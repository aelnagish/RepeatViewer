"""Data model v2: slides, channels, chromosome marks and signals."""

from __future__ import annotations

import re

import json
import os
import uuid
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path

from ..cyto_images import CytoImage, store_dir

SLIDES_FILE = "slides.json"

ARMS = ("p", "cen", "q")
POSITIONS = ("pericentromeric", "interstitial", "subterminal", "terminal")
STRENGTHS = ("strong", "medium", "weak")
UNNAMED = "?"  # chromosome name of a detected chromosome that is not paired yet


def new_id() -> str:
    return uuid.uuid4().hex


def _build(cls, raw):
    """Dataclass from a dict, ignoring unknown keys (newer files) and bad input."""
    if not isinstance(raw, dict):
        return None
    names = {f.name for f in fields(cls)}
    try:
        return cls(**{k: v for k, v in raw.items() if k in names})
    except (TypeError, ValueError):
        return None


def _point(value) -> tuple[float, float] | None:
    if isinstance(value, (list, tuple)) and len(value) == 2:
        try:
            return (float(value[0]), float(value[1]))
        except (TypeError, ValueError):
            return None
    return None


@dataclass
class ChannelDisplay:
    color: str = "#ffffff"
    visible: bool = True
    low: float = 0.0
    high: float = 0.0  # high <= low means "not set yet": use auto levels
    gamma: float = 1.0

    @property
    def is_set(self) -> bool:
        return self.high > self.low


@dataclass
class Channel:
    id: str
    file: str  # file name inside the store directory
    dye: str = ""
    probe: str = ""
    bits: int = 8
    original_name: str = ""
    display: ChannelDisplay = field(default_factory=ChannelDisplay)
    index: int = -1  # plane in a multichannel or RGB file; -1 = the whole first page

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, raw) -> "Channel | None":
        if not isinstance(raw, dict) or not raw.get("id") or not raw.get("file"):
            return None
        data = dict(raw)
        data["display"] = _build(ChannelDisplay, raw.get("display")) or ChannelDisplay()
        return _build(cls, data)


@dataclass
class ChromosomeMark:
    """One numbered chromosome in one image (e.g. "2b")."""

    id: str
    chromosome: str  # karyotype chromosome name, e.g. "2"
    homologue: str = ""  # "a", "b", ...
    x: float = 0.0  # label anchor (image pixels)
    y: float = 0.0
    centromere: tuple[float, float] | None = None
    telomere_p: tuple[float, float] | None = None
    telomere_q: tuple[float, float] | None = None
    # found by detection (Signals step 1): outline and midline from p to q telomere
    midline: list | None = None  # [[x, y], ...]
    mask: dict | None = None  # segment.encode_mask output
    auto: bool = False
    check: str = ""  # why it should be looked at ("" = fine or accepted)

    @property
    def label(self) -> str:
        return f"{self.chromosome}{self.homologue}"

    @property
    def named(self) -> bool:
        """False for chromosomes found by detection and not yet paired with the karyotype ("?")."""
        return self.chromosome != UNNAMED

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, raw) -> "ChromosomeMark | None":
        if not isinstance(raw, dict) or not raw.get("id") or not str(raw.get("chromosome", "")).strip():
            return None
        data = dict(raw)
        for key in ("centromere", "telomere_p", "telomere_q"):
            data[key] = _point(raw.get(key))
        line = raw.get("midline")
        data["midline"] = [p for p in (_point(q) for q in line) if p] if isinstance(line, list) else None
        data["mask"] = raw.get("mask") if isinstance(raw.get("mask"), dict) else None
        return _build(cls, data)


@dataclass
class Signal:
    id: str
    probe: str
    x: float
    y: float
    chromosome_id: str = ""  # ChromosomeMark.id, "" when not assigned
    arm: str = ""  # p / cen / q
    position: str = ""  # see POSITIONS
    rel_distance: float | None = None  # 0 at the centromere, 1 at the telomere
    strength: str = ""
    manual: bool = False  # arm/position set by hand, never recomputed
    auto: bool = False  # found by detection (replaced when the probe is detected again)
    check: str = ""  # why it should be looked at
    peak: float | None = None  # detection strength, in thresholds

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, raw) -> "Signal | None":
        if not isinstance(raw, dict) or not raw.get("id"):
            return None
        return _build(cls, raw)


@dataclass
class ImageDetails:
    """The v2 part of one image record."""

    slide_id: str = ""
    karyotype: str = ""  # karyotype name; "" = the species' default
    pixel_size_um: float | None = None
    channels: list[Channel] = field(default_factory=list)
    chromosomes: list[ChromosomeMark] = field(default_factory=list)
    signals: list[Signal] = field(default_factory=list)
    adjustment_log: list[str] = field(default_factory=list)
    detection: dict = field(default_factory=dict)  # Signals detection settings and probe table of this cell

    @property
    def is_multichannel(self) -> bool:
        return bool(self.channels)

    def chromosome(self, mark_id: str) -> ChromosomeMark | None:
        return next((c for c in self.chromosomes if c.id == mark_id), None)

    def to_dict(self) -> dict:
        return {
            "slide_id": self.slide_id,
            "karyotype": self.karyotype,
            "pixel_size_um": self.pixel_size_um,
            "channels": [c.to_dict() for c in self.channels],
            "chromosomes": [c.to_dict() for c in self.chromosomes],
            "signals": [s.to_dict() for s in self.signals],
            "adjustment_log": list(self.adjustment_log),
            "detection": dict(self.detection),
        }

    @classmethod
    def from_extra(cls, extra: dict | None) -> "ImageDetails":
        extra = extra if isinstance(extra, dict) else {}
        try:
            pixel = float(extra["pixel_size_um"]) if extra.get("pixel_size_um") else None
        except (TypeError, ValueError):
            pixel = None
        return cls(
            slide_id=str(extra.get("slide_id") or ""),
            karyotype=str(extra.get("karyotype") or ""),
            pixel_size_um=pixel,
            channels=[c for c in (Channel.from_dict(x) for x in extra.get("channels") or []) if c],
            chromosomes=[c for c in (ChromosomeMark.from_dict(x) for x in extra.get("chromosomes") or []) if c],
            signals=[s for s in (Signal.from_dict(x) for x in extra.get("signals") or []) if s],
            adjustment_log=[str(x) for x in extra.get("adjustment_log") or []],
            detection=dict(extra.get("detection") or {}) if isinstance(extra.get("detection"), dict) else {},
        )


def image_details(image: CytoImage) -> ImageDetails:
    return ImageDetails.from_extra(image.extra)


def set_image_details(image: CytoImage, details: ImageDetails) -> None:
    image.extra.update(details.to_dict())


_DYE_RANK = {"dapi": 0, "hoechst": 0, "pi": 1, "fitc": 2, "gfp": 2, "alexa488": 2, "cy3": 3, "tritc": 3, "alexa555": 3,
             "texasred": 4, "alexa594": 4, "cy5": 5, "alexa647": 5, "cy7": 6}


def ordered_channels(details: ImageDetails) -> list:
    """Counterstain first, then by wavelength (DAPI, FITC, Cy3, Cy5), others last in their order."""
    def rank(item):
        i, c = item
        key = re.sub(r"[^a-z0-9]", "", c.dye.lower())
        return (_DYE_RANK.get(key, 50), i)
    return [c for _i, c in sorted(enumerate(details.channels), key=rank)]



def channel_path(folder: Path, channel: Channel) -> Path:
    """Copied channels live in the store; channels imported without copying keep their absolute path."""
    if Path(channel.file).is_absolute():
        return Path(channel.file)
    return store_dir(Path(folder).expanduser().resolve()) / Path(channel.file).name


# ---------------------------------------------------------------------------
# Slides
# ---------------------------------------------------------------------------


@dataclass
class Slide:
    id: str
    name: str
    date: str = ""
    species: str = ""
    accession: str = ""
    technique: str = ""
    counterstain: str = ""
    note: str = ""
    box: str = ""  # slide or box number
    probes: str = ""  # probes / hybridisation, free text
    hyb_temp: str = ""  # hybridisation temperature
    fixation_date: str = ""
    microscope: str = ""
    experimenter: str = ""
    protocol: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


# Slide values an image takes unless it was given its own value
SLIDE_INHERITED = ("species", "accession", "technique", "counterstain", "date", "microscope", "experimenter")


def _slides_path(folder: Path) -> Path:
    return store_dir(Path(folder).expanduser().resolve()) / SLIDES_FILE


def load_slides(folder: Path) -> list[Slide]:
    try:
        raw = json.loads(_slides_path(folder).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    items = raw.get("slides") if isinstance(raw, dict) else None
    out = []
    for item in items or []:
        slide = _build(Slide, item)
        if slide is not None and slide.id:
            out.append(slide)
    return out


def save_slides(folder: Path, slides: list[Slide]) -> None:
    path = _slides_path(folder)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    payload = {"version": 1, "slides": [s.to_dict() for s in slides]}
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def add_slide(folder: Path, name: str, **values) -> Slide:
    slides = load_slides(folder)
    slide = Slide(id=new_id(), name=str(name).strip() or "Slide", **values)
    slides.append(slide)
    save_slides(folder, slides)
    return slide


def _inherit(image: CytoImage, old: Slide | None, new: Slide | None) -> bool:
    """Give the image the new slide's values where it still has the old slide's value (or none)."""
    changed = False
    for key in SLIDE_INHERITED:
        before = getattr(old, key, "") if old is not None else ""
        value = getattr(new, key, "") if new is not None else ""
        current = image.metadata.get(key, "")
        if current in ("", before) and current != value:
            if value:
                image.metadata[key] = value
            else:
                image.metadata.pop(key, None)
            changed = True
    return changed


def update_slide(folder: Path, slide: Slide) -> int:
    """Save the slide's details; its images follow unless they were given their own value. Returns images changed."""
    from ..cyto_images import load_images, save_images

    slides = load_slides(folder)
    old = next((s for s in slides if s.id == slide.id), None)
    slides = [slide if s.id == slide.id else s for s in slides] if old else [*slides, slide]
    save_slides(folder, slides)
    images = load_images(folder)
    changed = [im for im in images if image_details(im).slide_id == slide.id and _inherit(im, old, slide)]
    if changed:
        save_images(folder, images)
    return len(changed)


def assign_to_slide(folder: Path, image_ids: list[str], slide_id: str) -> int:
    """Put images on a slide ("" = no slide). Empty image values are filled from the new slide."""
    from ..cyto_images import load_images, save_images

    slides = {s.id: s for s in load_slides(folder)}
    new = slides.get(slide_id)
    images = load_images(folder)
    n = 0
    for im in images:
        if im.id not in image_ids:
            continue
        d = image_details(im)
        if d.slide_id == (slide_id or ""):
            continue
        old = slides.get(d.slide_id)
        _inherit(im, old, new)
        d.slide_id = slide_id or ""
        set_image_details(im, d)
        n += 1
    if n:
        save_images(folder, images)
    return n


def slide_for_image(image: CytoImage, slides: list[Slide]) -> Slide | None:
    sid = image_details(image).slide_id
    return next((s for s in slides if s.id == sid), None) if sid else None
