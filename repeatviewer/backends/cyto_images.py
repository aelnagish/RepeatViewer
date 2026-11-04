"""Cytogenetic images (FISH, GISH, …) attached to sequences and probe annotations.

Qt-free storage layer. Each data folder keeps its images in a hidden
``.repeatviewer_cytogenetics/`` directory next to the FASTA files:

    my_folder/
        CroSat1.fasta
        .repeatviewer_cytogenetics/
            index.json
            3f2a…c1.png        (copy of the dropped image)

``index.json`` holds one record per image. A record can be linked to several
targets (multicolour FISH: one picture, several probes), each target being a
FASTA file in the same folder and, optionally, one annotation on it.
"""

from __future__ import annotations

import json
import os
import shutil
import uuid
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Iterable

STORE_DIR_NAME = ".repeatviewer_cytogenetics"
INDEX_FILE_NAME = "index.json"
# Version 2 adds per-image fields (channels, signals, ...) kept in CytoImage.extra;
# version 1 files load unchanged.
INDEX_VERSION = 2

IMAGE_SUFFIXES = (".png", ".jpg", ".jpeg", ".tif", ".tiff", ".bmp", ".gif", ".webp")

# Ordered metadata fields: (key, label). All optional.
METADATA_FIELDS: tuple[tuple[str, str], ...] = (
    ("species", "Species"),
    ("accession", "Accession / line"),
    ("technique", "Technique"),
    ("counterstain", "Counterstain"),
    ("stage", "Cell stage"),
    ("chromosomes", "Chromosomes with signal"),
    ("date", "Date"),
    ("experimenter", "Experimenter"),
    ("microscope", "Microscope / objective"),
    ("notes", "Notes"),
)
METADATA_KEYS = tuple(key for key, _label in METADATA_FIELDS)

TECHNIQUE_CHOICES = ("FISH", "GISH", "BAC-FISH", "Oligo-FISH / painting", "Immunostaining")
COUNTERSTAIN_CHOICES = ("DAPI", "PI", "Hoechst")
STAGE_CHOICES = ("Mitotic metaphase", "Prometaphase", "Pachytene", "Interphase", "Extended fibre")
DYE_CHOICES = (
    ("Cy3", "#ff3b30"),
    ("Texas Red", "#e0282e"),
    ("Alexa 594", "#ff453a"),
    ("FITC", "#34c759"),
    ("Alexa 488", "#30d158"),
    ("Cy5", "#bf5af2"),
    ("Atto 647N", "#af52de"),
    ("DEAC", "#0a84ff"),
    ("Aqua", "#5ac8fa"),
)


def is_image_path(path: Path) -> bool:
    return path.suffix.lower() in IMAGE_SUFFIXES


def default_dye_color(dye: str) -> str:
    text = str(dye or "").strip().lower()
    for name, color in DYE_CHOICES:
        if name.lower() == text:
            return color
    return ""


# ---------------------------------------------------------------------------
# Data model
# ---------------------------------------------------------------------------


@dataclass
class FeatureRef:
    """Annotation on a sequence, identified by name + 1-based coordinates."""

    name: str
    start: int
    end: int

    def __post_init__(self) -> None:
        self.name = str(self.name or "").strip()
        self.start = int(self.start)
        self.end = int(self.end)
        if self.start > self.end:
            self.start, self.end = self.end, self.start

    def to_dict(self) -> dict:
        return {"name": self.name, "start": self.start, "end": self.end}

    @classmethod
    def from_dict(cls, raw: dict | None) -> "FeatureRef | None":
        if not isinstance(raw, dict):
            return None
        try:
            return cls(raw.get("name", ""), int(raw["start"]), int(raw["end"]))
        except (KeyError, TypeError, ValueError):
            return None


@dataclass
class ImageTarget:
    """One probe shown in an image: a FASTA (and optionally one annotation)."""

    sequence: str  # FASTA file name, relative to the folder
    feature: FeatureRef | None = None
    probe: str = ""
    dye: str = ""
    color: str = ""

    def __post_init__(self) -> None:
        self.sequence = str(self.sequence or "").strip()
        self.probe = str(self.probe or "").strip()
        self.dye = str(self.dye or "").strip()
        self.color = str(self.color or "").strip()

    @property
    def probe_label(self) -> str:
        if self.probe:
            return self.probe
        if self.feature is not None and self.feature.name:
            return self.feature.name
        return Path(self.sequence).stem

    def to_dict(self) -> dict:
        out: dict = {"sequence": self.sequence}
        if self.feature is not None:
            out["feature"] = self.feature.to_dict()
        for key in ("probe", "dye", "color"):
            value = getattr(self, key)
            if value:
                out[key] = value
        return out

    @classmethod
    def from_dict(cls, raw: dict) -> "ImageTarget | None":
        if not isinstance(raw, dict) or not raw.get("sequence"):
            return None
        return cls(
            sequence=raw.get("sequence", ""),
            feature=FeatureRef.from_dict(raw.get("feature")),
            probe=raw.get("probe", ""),
            dye=raw.get("dye", ""),
            color=raw.get("color", ""),
        )


_KNOWN_IMAGE_KEYS = frozenset({"id", "file", "original_name", "added", "targets", "metadata"})


@dataclass
class CytoImage:
    """One stored image with its links and optional metadata."""

    id: str
    file: str  # file name inside the store directory
    original_name: str = ""
    added: str = ""
    targets: list[ImageTarget] = field(default_factory=list)
    metadata: dict[str, str] = field(default_factory=dict)
    # Fields this module does not know (v2 cytogenetics data: channels, signals, ...).
    # Kept as-is so saving here never drops them.
    extra: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return {
            **self.extra,
            "id": self.id,
            "file": self.file,
            "original_name": self.original_name,
            "added": self.added,
            "targets": [t.to_dict() for t in self.targets],
            "metadata": {k: v for k, v in self.metadata.items() if str(v).strip()},
        }

    @classmethod
    def from_dict(cls, raw: dict) -> "CytoImage | None":
        if not isinstance(raw, dict) or not raw.get("id") or not raw.get("file"):
            return None
        targets = [t for t in (ImageTarget.from_dict(x) for x in raw.get("targets") or []) if t]
        meta_raw = raw.get("metadata") or {}
        metadata = {str(k): str(v) for k, v in meta_raw.items()} if isinstance(meta_raw, dict) else {}
        return cls(
            id=str(raw["id"]),
            file=str(raw["file"]),
            original_name=str(raw.get("original_name", "")),
            added=str(raw.get("added", "")),
            targets=targets,
            metadata=metadata,
            extra={k: v for k, v in raw.items() if k not in _KNOWN_IMAGE_KEYS},
        )

    def targets_for_sequence(self, sequence_name: str) -> list[ImageTarget]:
        return [t for t in self.targets if t.sequence == sequence_name]

    def title(self) -> str:
        probes = []
        for target in self.targets:
            label = target.probe_label
            if label and label not in probes:
                probes.append(label)
        return " + ".join(probes) if probes else (self.original_name or self.file)


# ---------------------------------------------------------------------------
# Store
# ---------------------------------------------------------------------------


def store_dir(folder: Path) -> Path:
    return Path(folder) / STORE_DIR_NAME


def index_path(folder: Path) -> Path:
    return store_dir(folder) / INDEX_FILE_NAME


_CACHE: dict[str, tuple[float, list[CytoImage]]] = {}
# Sentinel mtime when index.json is missing — avoids re-stat'ing per FASTA row
# while opening large folders with no cytogenetic store.
_MISSING_INDEX_MTIME = -1.0


def _index_mtime(folder: Path) -> float | None:
    try:
        return index_path(folder).stat().st_mtime
    except OSError:
        return None


def load_images(folder: Path) -> list[CytoImage]:
    """All image records for ``folder`` (cached; re-read when index.json changes)."""
    folder = Path(folder).expanduser().resolve()
    key = str(folder)
    mtime = _index_mtime(folder)
    if mtime is None:
        cached = _CACHE.get(key)
        if cached is not None and cached[0] == _MISSING_INDEX_MTIME:
            return cached[1]
        _CACHE[key] = (_MISSING_INDEX_MTIME, [])
        return []
    cached = _CACHE.get(key)
    if cached is not None and cached[0] == mtime:
        return cached[1]
    try:
        raw = json.loads(index_path(folder).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        _CACHE[key] = (_MISSING_INDEX_MTIME, [])
        return []
    items = raw.get("images") if isinstance(raw, dict) else None
    images = [img for img in (CytoImage.from_dict(x) for x in items or []) if img]
    _CACHE[key] = (mtime, images)
    return images


def image_counts_by_sequence(folder: Path) -> dict[str, int]:
    """How many images link to each FASTA basename in ``folder`` (one index load)."""
    counts: dict[str, int] = {}
    for img in load_images(folder):
        seen: set[str] = set()
        for target in img.targets:
            name = str(getattr(target, "sequence", "") or "").strip()
            if not name or name in seen:
                continue
            seen.add(name)
            counts[name] = counts.get(name, 0) + 1
    return counts


def save_images(folder: Path, images: Iterable[CytoImage]) -> None:
    folder = Path(folder).expanduser().resolve()
    directory = store_dir(folder)
    directory.mkdir(parents=True, exist_ok=True)
    payload = {"version": INDEX_VERSION, "images": [img.to_dict() for img in images]}
    target = index_path(folder)
    tmp = target.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, target)
    _CACHE.pop(str(folder), None)


def image_file_path(folder: Path, image: CytoImage) -> Path:
    return store_dir(Path(folder).expanduser().resolve()) / image.file


def image_file_date(path: Path) -> str:
    """Best-effort acquisition date (YYYY-MM-DD) from the file timestamps."""
    try:
        st = Path(path).stat()
    except OSError:
        return ""
    stamp = getattr(st, "st_birthtime", None) or st.st_mtime
    return datetime.fromtimestamp(min(stamp, st.st_mtime)).strftime("%Y-%m-%d")


def add_image(
    folder: Path,
    source: Path,
    *,
    targets: list[ImageTarget],
    metadata: dict[str, str] | None = None,
) -> CytoImage:
    """Copy ``source`` into the folder's store and record it."""
    folder = Path(folder).expanduser().resolve()
    source = Path(source).expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(source)
    image_id = uuid.uuid4().hex
    suffix = source.suffix.lower() or ".png"
    directory = store_dir(folder)
    directory.mkdir(parents=True, exist_ok=True)
    shutil.copy2(source, directory / f"{image_id}{suffix}")
    image = CytoImage(
        id=image_id,
        file=f"{image_id}{suffix}",
        original_name=source.name,
        added=datetime.now().isoformat(timespec="seconds"),
        targets=list(targets),
        metadata={k: str(v) for k, v in (metadata or {}).items() if str(v).strip()},
    )
    images = list(load_images(folder))
    images.append(image)
    save_images(folder, images)
    return image


def update_image(folder: Path, image: CytoImage) -> None:
    images = list(load_images(folder))
    for i, existing in enumerate(images):
        if existing.id == image.id:
            images[i] = image
            break
    else:
        raise KeyError(image.id)
    save_images(folder, images)


def remove_image(folder: Path, image_id: str) -> bool:
    """Delete one record and its image file. True if something was removed."""
    folder = Path(folder).expanduser().resolve()
    images = list(load_images(folder))
    keep = [img for img in images if img.id != image_id]
    if len(keep) == len(images):
        return False
    directory = store_dir(folder)
    for img in images:
        if img.id != image_id:
            continue
        names = [img.file]
        for channel in img.extra.get("channels") or []:
            if isinstance(channel, dict) and channel.get("file"):
                names.append(str(channel["file"]))
        for name in names:
            if Path(name).is_absolute():
                continue  # channel left in place at import (not copied): never delete the user's file
            try:
                (directory / Path(name).name).unlink()
            except OSError:
                pass
    save_images(folder, keep)
    return True


def unlink_target(folder: Path, image_id: str, sequence_name: str) -> None:
    """Remove the links from one image to one sequence; drop the image if none remain."""
    images = list(load_images(folder))
    for img in images:
        if img.id != image_id:
            continue
        img.targets = [t for t in img.targets if t.sequence != sequence_name]
        if not img.targets:
            remove_image(folder, image_id)
            return
        save_images(folder, images)
        return


# ---------------------------------------------------------------------------
# Queries
# ---------------------------------------------------------------------------


def images_for_sequence(fasta_path: Path) -> list[CytoImage]:
    """Images linked to this FASTA (whole sequence or any of its annotations)."""
    fasta_path = Path(fasta_path).expanduser().resolve()
    name = fasta_path.name
    return [img for img in load_images(fasta_path.parent) if img.targets_for_sequence(name)]


def image_count_for_sequence(fasta_path: Path) -> int:
    return len(images_for_sequence(fasta_path))


def feature_matches(ref: FeatureRef, name: str, start: int, end: int) -> bool:
    return ref.name == str(name or "").strip() and ref.start == int(start) and ref.end == int(end)


def images_for_feature(fasta_path: Path, name: str, start: int, end: int) -> list[CytoImage]:
    """Images linked to one annotation.

    Exact name + coordinates first; if the annotation was moved or resized
    since the image was attached, fall back to a unique name match.
    """
    fasta_path = Path(fasta_path).expanduser().resolve()
    seq_images = images_for_sequence(fasta_path)
    exact: list[CytoImage] = []
    by_name: list[CytoImage] = []
    clean = str(name or "").strip()
    for img in seq_images:
        for target in img.targets_for_sequence(fasta_path.name):
            if target.feature is None:
                continue
            if feature_matches(target.feature, clean, start, end):
                exact.append(img)
                break
            if clean and target.feature.name == clean:
                by_name.append(img)
                break
    return exact or by_name


def feature_keys_with_images(fasta_path: Path) -> tuple[set[tuple[str, int, int]], set[str]]:
    """(exact keys, names) of annotations on this FASTA that carry images."""
    fasta_path = Path(fasta_path).expanduser().resolve()
    exact: set[tuple[str, int, int]] = set()
    names: set[str] = set()
    for img in images_for_sequence(fasta_path):
        for target in img.targets_for_sequence(fasta_path.name):
            if target.feature is not None:
                exact.add((target.feature.name, target.feature.start, target.feature.end))
                if target.feature.name:
                    names.add(target.feature.name)
    return exact, names


def rename_sequence(folder: Path, old_name: str, new_name: str) -> int:
    """Follow a FASTA rename. Returns the number of links updated."""
    images = list(load_images(folder))
    changed = 0
    for img in images:
        for target in img.targets:
            if target.sequence == old_name:
                target.sequence = new_name
                changed += 1
    if changed:
        save_images(folder, images)
    return changed


def rename_sequences(folder: Path, mapping: dict[str, str]) -> int:
    """Follow several FASTA renames in one folder at once (old name -> new name).

    All links change together, so swaps (a <-> b) and chains (a -> b, b -> c)
    end up right. Returns the number of links updated.
    """
    if not mapping:
        return 0
    images = list(load_images(folder))
    changed = 0
    for img in images:
        for target in img.targets:
            new_name = mapping.get(target.sequence)
            if new_name is not None and new_name != target.sequence:
                target.sequence = new_name
                changed += 1
    if changed:
        save_images(folder, images)
    return changed


# ---------------------------------------------------------------------------
# Caption
# ---------------------------------------------------------------------------


def build_caption(image: CytoImage) -> str:
    """Figure caption from the metadata, e.g.

    "CroSat1 (Cy3, red) and CroSat3 (FITC, green) on Crocus sativus mitotic
    metaphase chromosomes, DAPI counterstain."
    """
    parts: list[str] = []
    for target in image.targets:
        label = target.probe_label
        detail = ", ".join(x for x in (target.dye, _color_word(target.color)) if x)
        parts.append(f"{label} ({detail})" if detail else label)
    probes = _join_and(parts)
    meta = image.metadata
    technique = meta.get("technique", "").strip()
    species = meta.get("species", "").strip()
    accession = meta.get("accession", "").strip()
    stage = meta.get("stage", "").strip()
    subject_bits = []
    if species:
        subject_bits.append(f"{species} {accession}".strip() if accession else species)
    if stage:
        subject_bits.append(stage[0].lower() + stage[1:])
    subject = " ".join(subject_bits)
    if subject and "chromosome" not in subject.lower() and stage.lower() not in {"interphase", "extended fibre"}:
        subject += " chromosomes"
    sentence = probes or "Cytogenetic image"
    if technique:
        sentence = f"{technique} of {sentence}" if probes else technique
    if subject:
        sentence += f" on {subject}"
    extras = []
    if meta.get("counterstain", "").strip():
        extras.append(f"{meta['counterstain'].strip()} counterstain")
    if meta.get("chromosomes", "").strip():
        extras.append(f"signals on {meta['chromosomes'].strip()}")
    if extras:
        sentence += ", " + ", ".join(extras)
    return sentence.rstrip(".") + "."


def _join_and(items: list[str]) -> str:
    if not items:
        return ""
    if len(items) == 1:
        return items[0]
    return ", ".join(items[:-1]) + " and " + items[-1]


_COLOR_WORDS = (
    ("red", (255, 0, 0)),
    ("green", (0, 200, 0)),
    ("blue", (0, 90, 255)),
    ("yellow", (255, 220, 0)),
    ("magenta", (220, 0, 220)),
    ("cyan", (0, 210, 230)),
    ("orange", (255, 140, 0)),
    ("white", (245, 245, 245)),
)


def _color_word(color: str) -> str:
    text = str(color or "").strip().lstrip("#")
    if len(text) != 6:
        return ""
    try:
        r, g, b = (int(text[i : i + 2], 16) for i in (0, 2, 4))
    except ValueError:
        return ""
    best = min(_COLOR_WORDS, key=lambda item: sum((a - c) ** 2 for a, c in zip((r, g, b), item[1])))
    return best[0]
