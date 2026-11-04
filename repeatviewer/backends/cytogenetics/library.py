"""Library tab logic: one row per image, filters, facet counts and grouping (Qt-free)."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from ..cyto_images import CytoImage, load_images
from .model import Slide, image_details, load_slides
from .importer import natural_key

GROUP_BYS = ("Slide", "Probe", "Species", "Cell stage", "Nothing")
FACETS = ("probe", "species", "stage", "signals")
NO_SLIDE = "No slide"
SIGNALS_MARKED = "Marked"
SIGNALS_NOT_MARKED = "Not marked yet"


@dataclass
class LibraryItem:
    image: CytoImage
    slide: Slide | None
    probes: list[str]
    species: str
    stage: str
    marked: bool
    channels: int

    @property
    def title(self) -> str:
        return self.image.original_name or self.image.title()

    @property
    def signals_label(self) -> str:
        return SIGNALS_MARKED if self.marked else SIGNALS_NOT_MARKED

    @property
    def date(self) -> str:
        return self.image.metadata.get("date") or (self.slide.date if self.slide else "")


@dataclass
class LibraryFilter:
    text: str = ""
    probe: set[str] = field(default_factory=set)      # empty set = no restriction
    species: set[str] = field(default_factory=set)
    stage: set[str] = field(default_factory=set)
    signals: set[str] = field(default_factory=set)
    sequence: str = ""                                 # FASTA name an image must link to

    def is_empty(self) -> bool:
        return not (self.text or self.probe or self.species or self.stage or self.signals or self.sequence)


def build_items(folder: Path) -> list[LibraryItem]:
    slides = {s.id: s for s in load_slides(folder)}
    items = []
    for image in load_images(folder):
        details = image_details(image)
        probes: list[str] = []
        for p in [t.probe_label for t in image.targets] + [c.probe for c in details.channels] + [s.probe for s in details.signals]:
            if p and p not in probes:
                probes.append(p)
        slide = slides.get(details.slide_id)
        items.append(LibraryItem(
            image=image,
            slide=slide,
            probes=probes,
            species=image.metadata.get("species") or (slide.species if slide else ""),
            stage=image.metadata.get("stage", ""),
            marked=bool(details.signals),
            channels=len(details.channels),
        ))
    return items


def _haystack(item: LibraryItem) -> str:
    parts = [item.title, item.image.original_name, item.species, item.stage, item.slide.name if item.slide else "",
             item.slide.note if item.slide else "", " ".join(item.probes), " ".join(t.sequence for t in item.image.targets)]
    parts += [str(v) for v in item.image.metadata.values()]
    return " ".join(parts).lower()


def matches(item: LibraryItem, flt: LibraryFilter, *, skip: str = "") -> bool:
    """``skip`` leaves one facet out (used for the counts shown beside that facet)."""
    if flt.sequence and not any(t.sequence == flt.sequence for t in item.image.targets):
        return False
    if flt.text:
        hay = _haystack(item)
        if not all(word in hay for word in flt.text.lower().split()):
            return False
    if skip != "probe" and flt.probe and not (set(item.probes) & flt.probe):
        return False
    if skip != "species" and flt.species and item.species not in flt.species:
        return False
    if skip != "stage" and flt.stage and item.stage not in flt.stage:
        return False
    if skip != "signals" and flt.signals and item.signals_label not in flt.signals:
        return False
    return True


def apply_filter(items: list[LibraryItem], flt: LibraryFilter) -> list[LibraryItem]:
    return [i for i in items if matches(i, flt)]


def _facet_values(item: LibraryItem, facet: str) -> list[str]:
    if facet == "probe":
        return item.probes
    if facet == "species":
        return [item.species] if item.species else []
    if facet == "stage":
        return [item.stage] if item.stage else []
    if facet == "signals":
        return [item.signals_label]
    raise ValueError(facet)


def facet_counts(items: list[LibraryItem], flt: LibraryFilter, facet: str) -> list[tuple[str, int]]:
    """Value and image count for one facet, counted with the other facets applied.

    Values the folder has but the filter currently hides show 0, so they stay clickable.
    """
    counts: dict[str, int] = {}
    for item in items:
        values = _facet_values(item, facet)
        visible = matches(item, flt, skip=facet)
        for v in values:
            counts[v] = counts.get(v, 0) + (1 if visible else 0)
    order = sorted(counts, key=lambda v: (-counts[v] if facet == "signals" else 0, natural_key(v)))
    if facet == "signals":
        order = [v for v in (SIGNALS_MARKED, SIGNALS_NOT_MARKED) if v in counts]
    return [(v, counts[v]) for v in order]


@dataclass
class LibraryGroup:
    title: str
    subtitle: str
    items: list[LibraryItem]


def _slide_subtitle(slide: Slide | None, items: list[LibraryItem]) -> str:
    bits = []
    if slide is not None:
        if slide.date:
            bits.append(slide.date)
        if slide.species:
            bits.append(slide.species)
    probes: list[str] = []
    for i in items:
        for p in i.probes:
            if p not in probes:
                probes.append(p)
    if probes:
        bits.append(", ".join(probes))
    bits.append(f"{len(items)} image{'s' if len(items) != 1 else ''}")
    return " · ".join(bits)


def _sorted_items(items: list[LibraryItem]) -> list[LibraryItem]:
    return sorted(items, key=lambda i: natural_key(i.title))


def group_items(items: list[LibraryItem], by: str) -> list[LibraryGroup]:
    if by == "Nothing":
        return [LibraryGroup("", f"{len(items)} image{'s' if len(items) != 1 else ''}", _sorted_items(items))] if items else []
    buckets: dict[str, list[LibraryItem]] = {}
    slides: dict[str, Slide | None] = {}
    for item in items:
        if by == "Slide":
            keys = [item.slide.name if item.slide else NO_SLIDE]
            slides.setdefault(keys[0], item.slide)
        elif by == "Probe":
            keys = item.probes or ["No probe"]
        elif by == "Species":
            keys = [item.species or "No species"]
        elif by == "Cell stage":
            keys = [item.stage or "No stage"]
        else:
            raise ValueError(by)
        for key in keys:
            buckets.setdefault(key, []).append(item)
    last = {NO_SLIDE, "No probe", "No species", "No stage"}
    groups = []
    for key in sorted(buckets, key=lambda k: (k in last, natural_key(k))):
        members = _sorted_items(buckets[key])
        groups.append(LibraryGroup(key, _slide_subtitle(slides.get(key), members) if by == "Slide"
                                   else f"{len(members)} image{'s' if len(members) != 1 else ''}", members))
    return groups


def details_rows(item: LibraryItem) -> list[tuple[str, str]]:
    """(label, value) rows for the details panel; empty values are left out."""
    img, slide = item.image, item.slide
    details = image_details(img)
    meta = img.metadata
    channels = []
    for c in details.channels:
        channels.append(f"{c.dye or 'channel'}" + (f" ({c.probe})" if c.probe else ""))
    pixel = f"{details.pixel_size_um:.4g} µm per pixel" if details.pixel_size_um else ""
    rows = [
        ("Slide", slide.name if slide else ""),
        ("Species", item.species),
        ("Accession", meta.get("accession", "") or (slide.accession if slide else "")),
        ("Technique", meta.get("technique", "") or (slide.technique if slide else "")),
        ("Counterstain", meta.get("counterstain", "") or (slide.counterstain if slide else "")),
        ("Cell stage", item.stage),
        ("Date", item.date),
        ("Channels", " · ".join(channels)),
        ("Pixel size", pixel),
        ("Probes", ", ".join(item.probes)),
        ("Sequences", ", ".join(dict.fromkeys(t.sequence for t in img.targets))),
        ("Signals", f"{len(details.signals)} on {len(details.chromosomes)} chromosomes" if details.signals else ""),
        ("Notes", meta.get("notes", "")),
    ]
    return [(k, v) for k, v in rows if str(v).strip()]
