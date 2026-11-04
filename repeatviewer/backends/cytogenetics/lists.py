"""Editable lists: probes, species, cell stages, techniques, counterstains, dyes.

Saved as ``lists.json`` in the folder's store. Defaults come from the constants in
``cyto_images`` so nothing is lost when the file does not exist yet. Renaming or
merging an entry updates every image (and slide) that uses it.
"""

from __future__ import annotations

import json
import os
import re
from pathlib import Path

from ..cyto_images import (
    COUNTERSTAIN_CHOICES,
    DYE_CHOICES,
    STAGE_CHOICES,
    TECHNIQUE_CHOICES,
    load_images,
    save_images,
    store_dir,
)
from .model import image_details, load_slides, save_slides, set_image_details

LISTS_FILE = "lists.json"
LIST_KEYS = ("probes", "species", "stages", "techniques", "counterstains", "dyes")
LIST_TITLES = {
    "probes": "Probes",
    "species": "Species",
    "stages": "Cell stages",
    "techniques": "Techniques",
    "counterstains": "Counterstains",
    "dyes": "Dyes and channels",
}
STAGE_MARKS = ("chromosomes", "bivalents", "nuclei", "fibres")

_DYE_WORDS = {
    "DAPI": "dapi, 405, blue",
    "Hoechst": "hoechst",
    "FITC": "fitc, gfp, 488, green",
    "Alexa 488": "a488, alexa488",
    "Cy3": "cy3, tritc, 555, red",
    "Texas Red": "texasred, txred",
    "Alexa 594": "a594, alexa594",
    "Cy5": "cy5, 647, farred",
    "Atto 647N": "atto647",
    "DEAC": "deac",
    "Aqua": "aqua",
}

Lists = dict[str, list[dict]]

# order of the lists in the Edit lists window: what you need first comes first
EDIT_ORDER = ("dyes", "probes", "species", "stages", "techniques", "counterstains")

# excitation and emission maxima (nm) of common dyes
DYE_WAVELENGTHS = {
    "dapi": (358, 461), "hoechst": (350, 461), "pi": (535, 617), "dy-415": (418, 467), "deac": (432, 472),
    "aqua": (433, 480), "fitc": (495, 519), "alexa 488": (490, 525), "gfp": (488, 507), "cy3": (550, 570),
    "tritc": (557, 576), "rhodamine": (552, 575), "texas red": (596, 615), "alexa 594": (590, 617),
    "cy5": (649, 666), "atto 647n": (646, 664), "alexa 647": (650, 668), "cy7": (743, 767),
}

# every cell stage a FISH lab meets, with what the Signals tab marks on it
STANDARD_STAGES = (
    ("Interphase", "nuclei"), ("Prophase", "chromosomes"), ("Prometaphase", "chromosomes"),
    ("Mitotic metaphase", "chromosomes"), ("Anaphase", "chromosomes"), ("Telophase", "chromosomes"),
    ("Leptotene", "chromosomes"), ("Zygotene", "bivalents"), ("Pachytene", "bivalents"), ("Diplotene", "bivalents"),
    ("Diakinesis", "bivalents"), ("Metaphase I", "bivalents"), ("Anaphase I", "chromosomes"),
    ("Telophase I", "chromosomes"), ("Metaphase II", "chromosomes"), ("Anaphase II", "chromosomes"),
    ("Telophase II", "chromosomes"), ("Tetrad", "nuclei"), ("Pollen mitosis", "chromosomes"),
    ("Extended fibre", "fibres"),
)


def complete_lists(lists: Lists) -> Lists:
    """Add the standard stages that are missing and the known wavelengths of dyes. Never removes anything."""
    have = {str(e["name"]).lower() for e in lists.get("stages", [])}
    for name, marks in STANDARD_STAGES:
        if name.lower() not in have:
            lists.setdefault("stages", []).append({"name": name, "marks": marks})
    for e in lists.get("dyes", []):
        known = DYE_WAVELENGTHS.get(str(e["name"]).lower())
        if known:
            e.setdefault("excitation", str(known[0]))
            e.setdefault("emission", str(known[1]))
    return lists


def default_lists() -> Lists:
    return complete_lists(_base_lists())


def _base_lists() -> Lists:
    stage_marks = {
        "Pachytene": "bivalents",
        "Interphase": "nuclei",
        "Extended fibre": "fibres",
    }
    dyes = [{"name": "DAPI", "color": "#5b74ff", "words": _DYE_WORDS["DAPI"]}]
    dyes += [{"name": name, "color": color, "words": _DYE_WORDS.get(name, name.lower())} for name, color in DYE_CHOICES]
    return {
        "probes": [],
        "species": [],
        "stages": [{"name": s, "marks": stage_marks.get(s, "chromosomes")} for s in STAGE_CHOICES],
        "techniques": [{"name": t} for t in TECHNIQUE_CHOICES],
        "counterstains": [{"name": c} for c in COUNTERSTAIN_CHOICES],
        "dyes": dyes,
    }


def _path(folder: Path) -> Path:
    return store_dir(Path(folder).expanduser().resolve()) / LISTS_FILE


def load_lists(folder: Path) -> Lists:
    lists = default_lists()
    try:
        raw = json.loads(_path(folder).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return lists
    saved = raw.get("lists") if isinstance(raw, dict) else None
    if not isinstance(saved, dict):
        return lists
    for key in LIST_KEYS:
        items = saved.get(key)
        if isinstance(items, list):
            lists[key] = [dict(x) for x in items if isinstance(x, dict) and str(x.get("name", "")).strip()]
    return complete_lists(lists)


def save_lists(folder: Path, lists: Lists) -> None:
    path = _path(folder)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".json.tmp")
    payload = {"version": 1, "lists": {k: lists.get(k, []) for k in LIST_KEYS}}
    tmp.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)


def names(lists: Lists, key: str) -> list[str]:
    return [str(e["name"]) for e in lists.get(key, [])]


def entry(lists: Lists, key: str, name: str) -> dict | None:
    return next((e for e in lists.get(key, []) if e.get("name") == name), None)


def add_entry(lists: Lists, key: str, name: str, **values) -> dict:
    """Add ``name`` (no-op if present, case-insensitive). Returns the entry."""
    name = str(name).strip()
    if not name:
        raise ValueError("Empty name")
    for e in lists.setdefault(key, []):
        if str(e["name"]).lower() == name.lower():
            return e
    new = {"name": name, **{k: v for k, v in values.items() if v not in (None, "")}}
    lists[key].append(new)
    return new


def ensure_used_values(folder: Path, lists: Lists) -> bool:
    """Add values used by images but missing from the lists (typed in a field). True if changed."""
    before = {k: len(v) for k, v in lists.items()}
    for key, counts in usage_counts(folder).items():
        for value in counts:
            add_entry(lists, key, value)
    return any(len(lists[k]) != before.get(k, 0) for k in lists)


# ---------------------------------------------------------------------------
# Where entries are used, rename and merge
# ---------------------------------------------------------------------------

_META_KEY = {"species": "species", "stages": "stage", "techniques": "technique", "counterstains": "counterstain"}
_SLIDE_KEY = {"species": "species", "techniques": "technique", "counterstains": "counterstain"}


def _image_values(image, key: str) -> list[str]:
    if key in _META_KEY:
        value = image.metadata.get(_META_KEY[key], "")
        return [value] if value else []
    details = image_details(image)
    if key == "probes":
        vals = [t.probe for t in image.targets] + [c.probe for c in details.channels] + [s.probe for s in details.signals]
    elif key == "dyes":
        vals = [t.dye for t in image.targets] + [c.dye for c in details.channels]
    else:
        vals = []
    return [v for v in vals if v]


def usage_counts(folder: Path) -> dict[str, dict[str, int]]:
    """Per list: value -> number of images using it."""
    out: dict[str, dict[str, int]] = {k: {} for k in LIST_KEYS}
    for image in load_images(folder):
        for key in LIST_KEYS:
            for value in set(_image_values(image, key)):
                out[key][value] = out[key].get(value, 0) + 1
    for slide in load_slides(folder):
        for key, attr in _SLIDE_KEY.items():
            value = getattr(slide, attr)
            if value:
                out[key].setdefault(value, 0)
    return out


def _replace_in_images(folder: Path, key: str, mapping: dict[str, str]) -> int:
    images = list(load_images(folder))
    changed = 0
    for image in images:
        touched = False
        if key in _META_KEY:
            mkey = _META_KEY[key]
            old = image.metadata.get(mkey, "")
            if old in mapping:
                image.metadata[mkey] = mapping[old]
                touched = True
        elif key in ("probes", "dyes"):
            attr = "probe" if key == "probes" else "dye"
            for target in image.targets:
                if getattr(target, attr) in mapping:
                    setattr(target, attr, mapping[getattr(target, attr)])
                    touched = True
            details = image_details(image)
            items = list(details.channels) + (list(details.signals) if key == "probes" else [])
            for item in items:
                if getattr(item, attr) in mapping:
                    setattr(item, attr, mapping[getattr(item, attr)])
                    touched = True
            if touched:
                set_image_details(image, details)
        changed += int(touched)
    if changed:
        save_images(folder, images)
    if key in _SLIDE_KEY:
        slides = load_slides(folder)
        attr = _SLIDE_KEY[key]
        hit = False
        for slide in slides:
            if getattr(slide, attr) in mapping:
                setattr(slide, attr, mapping[getattr(slide, attr)])
                hit = True
        if hit:
            save_slides(folder, slides)
    return changed


def rename_entry(folder: Path, lists: Lists, key: str, old: str, new: str) -> int:
    """Rename in the list and in every image. Returns how many images changed."""
    new = str(new).strip()
    if not new:
        raise ValueError("Empty name")
    if new != old and entry(lists, key, new) is not None:
        raise ValueError(f"{new} is already in the list; use Merge instead.")
    e = entry(lists, key, old)
    if e is not None:
        e["name"] = new
    return _replace_in_images(folder, key, {old: new})


def merge_entries(folder: Path, lists: Lists, key: str, merged: list[str], into: str) -> int:
    """Replace every name in ``merged`` by ``into`` and drop the merged entries."""
    mapping = {m: into for m in merged if m != into}
    lists[key] = [e for e in lists.get(key, []) if e.get("name") not in mapping]
    if entry(lists, key, into) is None:
        add_entry(lists, key, into)
    return _replace_in_images(folder, key, mapping)


def remove_entry(folder: Path, lists: Lists, key: str, name: str) -> None:
    """Remove from the list only when no image uses it."""
    if usage_counts(folder).get(key, {}).get(name):
        raise ValueError(f"{name} is still used by images; merge it into another entry instead.")
    lists[key] = [e for e in lists.get(key, []) if e.get("name") != name]


# ---------------------------------------------------------------------------
# Dyes from file names
# ---------------------------------------------------------------------------


def filename_tokens(name: str) -> list[str]:
    stem = Path(name).name
    for suffix in (".ome.tiff", ".ome.tif", ".tiff", ".tif", ".jpeg", ".jpg", ".png", ".bmp", ".gif", ".webp"):
        if stem.lower().endswith(suffix):
            stem = stem[: -len(suffix)]
            break
    return [t for t in re.split(r"[\s_.\-]+", stem) if t]


def dye_words(lists: Lists) -> dict[str, str]:
    """lower-case word -> dye name."""
    out: dict[str, str] = {}
    for e in lists.get("dyes", []):
        words = [e["name"]] + [w for w in str(e.get("words", "")).split(",")]
        for w in words:
            w = re.sub(r"[\s_.\-]+", "", w).lower()
            if w and w not in out:
                out[w] = str(e["name"])
    return out


def dye_for_filename(name: str, lists: Lists) -> tuple[str, int]:
    """(dye, token index) for the last token naming a dye; ("", -1) if none."""
    words = dye_words(lists)
    tokens = filename_tokens(name)
    for i in range(len(tokens) - 1, -1, -1):
        dye = words.get(tokens[i].lower())
        if dye:
            return dye, i
    return "", -1


def dye_color(lists: Lists, dye: str, fallback: str = "#ffffff") -> str:
    e = entry(lists, "dyes", dye)
    return str(e.get("color") or fallback) if e else fallback


# ---------------------------------------------------------------------------
# Editing session (Edit lists window): changes are kept until Save
# ---------------------------------------------------------------------------

import copy


class ListsEditor:
    """Working copy of the lists. Renames and merges are applied to the images on ``apply``."""

    def __init__(self, folder: Path, lists: Lists | None = None) -> None:
        self.folder = Path(folder)
        self.lists: Lists = copy.deepcopy(lists if lists is not None else load_lists(folder))
        self.usage = usage_counts(folder)
        self._ops: list[tuple] = []  # ("rename", key, old, new) / ("merge", key, [names], into)

    def names(self, key: str) -> list[str]:
        return names(self.lists, key)

    def used_by(self, key: str, name: str) -> int:
        return self.usage.get(key, {}).get(name, 0)

    def _exists(self, key: str, name: str, ignore: str = "") -> bool:
        return any(str(e["name"]).lower() == name.lower() and e["name"] != ignore for e in self.lists.get(key, []))

    def add(self, key: str, name: str, **values) -> dict:
        name = str(name).strip()
        if not name:
            raise ValueError("The name is empty.")
        if self._exists(key, name):
            raise ValueError(f"{name} is already in the list.")
        return add_entry(self.lists, key, name, **values)

    def move(self, key: str, src: int, dst: int) -> None:
        """Change the order of the entries (drag in the table)."""
        items = self.lists.get(key, [])
        if 0 <= src < len(items) and 0 <= dst < len(items) and src != dst:
            items.insert(dst, items.pop(src))

    def rename(self, key: str, old: str, new: str) -> None:
        new = str(new).strip()
        if not new:
            raise ValueError("The name is empty.")
        if new == old:
            return
        if self._exists(key, new, ignore=old):
            raise ValueError(f"{new} is already in the list. Use Merge to combine entries.")
        entry(self.lists, key, old)["name"] = new
        self.usage.setdefault(key, {})[new] = self.usage.get(key, {}).pop(old, 0)
        self._ops.append(("rename", key, old, new))

    def set_value(self, key: str, name: str, field_name: str, value: str) -> None:
        e = entry(self.lists, key, name)
        if e is None:
            raise KeyError(name)
        if value in (None, ""):
            e.pop(field_name, None)
        else:
            e[field_name] = value

    def merge(self, key: str, merged: list[str], into: str) -> None:
        merged = [m for m in merged if m != into]
        if not merged:
            return
        total = self.used_by(key, into) + sum(self.used_by(key, m) for m in merged)
        self.lists[key] = [e for e in self.lists.get(key, []) if e.get("name") not in merged]
        if entry(self.lists, key, into) is None:
            add_entry(self.lists, key, into)
        usage = self.usage.setdefault(key, {})
        for m in merged:
            usage.pop(m, None)
        usage[into] = total
        self._ops.append(("merge", key, merged, into))

    def can_remove(self, key: str, name: str) -> bool:
        return self.used_by(key, name) == 0

    def remove(self, key: str, name: str) -> None:
        if not self.can_remove(key, name):
            raise ValueError(f"{name} is used by {self.used_by(key, name)} image(s). Merge it into another entry instead.")
        self.lists[key] = [e for e in self.lists.get(key, []) if e.get("name") != name]

    def apply(self) -> int:
        """Update the images for every rename and merge, then save the lists. Returns images changed."""
        changed = 0
        for op in self._ops:
            if op[0] == "rename":
                changed += _replace_in_images(self.folder, op[1], {op[2]: op[3]})
            else:
                changed += _replace_in_images(self.folder, op[1], {m: op[3] for m in op[2]})
        save_lists(self.folder, self.lists)
        self._ops.clear()
        return changed
