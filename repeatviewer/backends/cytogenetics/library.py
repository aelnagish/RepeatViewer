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
    channel_colours: list[str] = field(default_factory=list)
    probe_info: list[tuple[str, str, str, str]] = field(default_factory=list)  # probe, dye, colour, sequence
    n_marks: int = 0
    n_signals: int = 0
    n_expected: int = 0  # chromosomes in the karyotype (2n), 0 when unknown

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
    from .lists import dye_color, entry, load_lists

    from .figure import ordered_channels
    from .karyotype import list_karyotypes

    slides = {s.id: s for s in load_slides(folder)}
    lists = load_lists(folder)
    karyotypes = list_karyotypes(folder)

    def probe_dye_colour(probe: str, dye: str, colour: str) -> tuple[str, str]:
        e = entry(lists, "probes", probe) or {}
        dye = dye or str(e.get("dye") or "")
        if not colour or colour.lower() == "#ffffff":
            colour = dye_color(lists, dye, "") if dye else ""
        return dye, colour
    items = []
    for image in load_images(folder):
        details = image_details(image)
        probes: list[str] = []
        for p in [t.probe_label for t in image.targets] + [c.probe for c in details.channels] + [s.probe for s in details.signals]:
            if p and p not in probes:
                probes.append(p)
        slide = slides.get(details.slide_id)
        sequences = {t.probe_label: t.sequence for t in image.targets if t.probe_label and t.sequence}
        info: list[tuple[str, str, str, str]] = []
        channels = ordered_channels(details)
        for c in channels:
            if c.probe and all(c.probe != row[0] for row in info):
                dye, colour = probe_dye_colour(c.probe, c.dye, c.display.color)
                info.append((c.probe, dye, colour, sequences.get(c.probe, "")))
        for probe in probes:
            if all(probe != row[0] for row in info):
                dye, colour = probe_dye_colour(probe, "", "")
                info.append((probe, dye, colour, sequences.get(probe, "")))
        items.append(LibraryItem(
            image=image,
            slide=slide,
            probes=probes,
            species=image.metadata.get("species") or (slide.species if slide else ""),
            stage=image.metadata.get("stage", ""),
            marked=bool(details.signals),
            channels=len(details.channels),
            channel_colours=[c.display.color for c in channels],
            probe_info=info,
            n_marks=sum(1 for m in details.chromosomes if m.named),
            n_signals=len(details.signals),
            n_expected=_expected(karyotypes, details.karyotype,
                                 image.metadata.get("species") or (slide.species if slide else "")),
        ))
    return items


def _expected(karyotypes, name: str, species: str) -> int:
    k = next((k for k in karyotypes if name and k.name == name), None)
    if k is None and species:
        k = next((k for k in karyotypes if k.species.lower() == species.lower()), None)
    return k.total_count if k is not None and k.chromosomes else 0


_MONTHS = ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec")


def pretty_date(text: str) -> str:
    """"2026-10-14" -> "14 Oct 2026"; anything else unchanged."""
    import re

    m = re.match(r"^(\d{4})-(\d{2})-(\d{2})", str(text or ""))
    if not m:
        return str(text or "")
    y, mo, d = int(m.group(1)), int(m.group(2)), int(m.group(3))
    return f"{d} {_MONTHS[mo - 1]} {y}" if 1 <= mo <= 12 else str(text)


def short_species(name: str) -> str:
    """"Crocus sativus" -> "C. sativus"."""
    parts = str(name or "").split()
    if len(parts) >= 2 and parts[0][:1].isupper() and not parts[0].endswith("."):
        return f"{parts[0][0]}. {' '.join(parts[1:])}"
    return str(name or "")


def _probe_sort_key(value: str):
    return (value[:1].isdigit(), natural_key(value))


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


def facet_values(items: list[LibraryItem], facet: str) -> set[str]:
    return {v for item in items for v in _facet_values(item, facet)}


def include_all_but(items: list[LibraryItem], facet: str, excluded: set[str]) -> set[str]:
    """Filter set for "everything ticked except ``excluded``" (empty set = no restriction)."""
    if not excluded:
        return set()
    keep = facet_values(items, facet) - excluded
    return keep or {"\x00none"}


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
    order = sorted(counts, key=_probe_sort_key)
    if facet == "signals":
        order = [v for v in (SIGNALS_MARKED, SIGNALS_NOT_MARKED) if v in counts]
    return [(v, counts[v]) for v in order]


@dataclass
class LibraryGroup:
    title: str
    subtitle: str
    items: list[LibraryItem]


def _slide_subtitle(slide: Slide | None, items: list[LibraryItem]) -> str:
    """"14 Oct 2026 · C. sativus · CroSat2 (FITC) + 45S rDNA (Cy3) · 9 cells" (as in the mockup)."""
    bits = []
    if slide is not None:
        if slide.date:
            bits.append(pretty_date(slide.date))
        species = slide.species or next((i.species for i in items if i.species), "")
        if species:
            bits.append(short_species(species))
    probes: list[str] = []
    for i in items:
        for probe, dye, _colour, _seq in i.probe_info:
            text = f"{probe} ({dye})" if dye else probe
            if text not in probes:
                probes.append(text)
    if probes:
        bits.append(" + ".join(probes))
    word = "cell" if slide is not None else "image"
    bits.append(f"{len(items)} {word}{'s' if len(items) != 1 else ''}")
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
    def order(k: str):
        if by == "Slide":
            slide = slides.get(k)
            date = slide.date if slide is not None else ""
            return (k in last, "".join(chr(0x10FFFF - ord(ch)) for ch in date) if date else "\U0010FFFF", natural_key(k))
        return (k in last, _probe_sort_key(k))

    for key in sorted(buckets, key=order):
        members = _sorted_items(buckets[key])
        groups.append(LibraryGroup(key, _slide_subtitle(slides.get(key), members) if by == "Slide"
                                   else f"{len(members)} image{'s' if len(members) != 1 else ''}", members))
    return groups


def details_rows(item: LibraryItem) -> list[tuple[str, str]]:
    """(label, value) rows for the details panel; empty values are left out."""
    from .figure import ordered_channels

    img, slide = item.image, item.slide
    details = image_details(img)
    meta = img.metadata
    channels = ordered_channels(details)
    copied = sum(1 for c in channels if not Path(c.file).is_absolute())
    channel_text = " · ".join(c.dye or "channel" for c in channels)
    if channels:
        channel_text += f" ({len(channels)} TIFF{'s' if len(channels) != 1 else ''}{'' if copied == len(channels) else ', not copied'})"
    pixel = f"{details.pixel_size_um:.4g} µm" if details.pixel_size_um else ""
    rows = [
        ("Slide", slide.name if slide else ""),
        ("Species", item.species),
        ("Accession", meta.get("accession", "") or (slide.accession if slide else "")),
        ("Technique", meta.get("technique", "") or (slide.technique if slide else "")),
        ("Counterstain", meta.get("counterstain", "") or (slide.counterstain if slide else "")),
        ("Cell stage", item.stage),
        ("Channels", channel_text),
        ("Pixel size", pixel),
        ("Date", pretty_date(item.date)),
        ("Probes", ", ".join(item.probes)),
        ("Sequences", ", ".join(dict.fromkeys(t.sequence for t in img.targets))),
        ("Signals", f"{len(details.signals)} on {len(details.chromosomes)} chromosomes" if details.signals else ""),
        ("Notes", meta.get("notes", "")),
    ]
    return [(k, v) for k, v in rows if str(v).strip()]


def probe_colours(items: list[LibraryItem]) -> dict[str, str]:
    """Probe -> colour of the channel it was seen in (for the dots in the filters)."""
    out: dict[str, str] = {}
    for item in items:
        for probe, _dye, colour, _seq in item.probe_info:
            if colour and probe not in out:
                out[probe] = colour
    return out


def signals_text(item: LibraryItem) -> str:
    if not item.n_marks and not item.n_signals:
        return "Not marked yet"
    numbered = f"{item.n_marks} of {item.n_expected}" if item.n_expected else f"{item.n_marks}"
    return f"{numbered} chromosomes numbered, {item.n_signals} signals marked"


def cell_label(item: LibraryItem) -> str:
    """"cell 03" for cells of a slide (from the file name), else the image title."""
    import re

    m = re.search(r"(\d+)\D*$", item.image.original_name or "")
    if item.slide is not None and m:
        return f"cell {m.group(1).zfill(2)}"
    return item.title


def full_title(item: LibraryItem) -> str:
    return f"{item.slide.name} · {cell_label(item)}" if item.slide else item.title
