"""The Cytogenetics workspace of one data folder.

The workspace page opens from a small marker file, ``workspace.rvcyto``, kept inside
the hidden store (``.repeatviewer_cytogenetics/``). It is not shown in the Files panel;
Tools > Cytogenetics creates it when needed and opens it in the Preview area, like a
Mapping result (``.rvmap``). Nothing but the format version is stored in it.
"""

from __future__ import annotations

import json
import os
from dataclasses import dataclass
from pathlib import Path

from ..cyto_images import STORE_DIR_NAME, load_images, store_dir
from .karyotype import list_karyotypes
from .model import image_details, load_slides

WORKSPACE_SUFFIX = ".rvcyto"
WORKSPACE_FILE = "workspace" + WORKSPACE_SUFFIX
TAB_NAMES = ("Library", "Image", "Signals", "Idiogram", "Figure")


def workspace_path(folder: Path) -> Path:
    return store_dir(Path(folder).expanduser().resolve()) / WORKSPACE_FILE


def is_workspace_path(path: Path) -> bool:
    return Path(path).suffix.lower() == WORKSPACE_SUFFIX


def workspace_folder(path: Path) -> Path:
    """The data folder a workspace marker belongs to."""
    path = Path(path).expanduser().resolve()
    return path.parent.parent if path.parent.name == STORE_DIR_NAME else path.parent


def ensure_workspace(folder: Path) -> Path:
    """Create the marker (and the store folder) if missing. Returns its path."""
    path = workspace_path(folder)
    if not path.is_file():
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(".tmp")
        tmp.write_text(json.dumps({"version": 1}) + "\n", encoding="utf-8")
        os.replace(tmp, path)
    return path


@dataclass
class WorkspaceSummary:
    images: int = 0
    slides: int = 0
    multichannel: int = 0
    with_signals: int = 0
    karyotypes: int = 0


def folder_summary(folder: Path) -> WorkspaceSummary:
    """Counts shown in the workspace's Data section."""
    images = load_images(folder)
    details = [image_details(i) for i in images]
    return WorkspaceSummary(
        images=len(images),
        slides=len(load_slides(folder)),
        multichannel=sum(1 for d in details if d.channels),
        with_signals=sum(1 for d in details if d.signals),
        karyotypes=len(list_karyotypes(folder)),
    )


def summary_lines(s: WorkspaceSummary) -> list[tuple[str, str]]:
    """(label, value) rows for the Data section."""
    return [
        ("Images", str(s.images)),
        ("Slides", str(s.slides)),
        ("With channels", str(s.multichannel)),
        ("With signals marked", str(s.with_signals)),
        ("Karyotypes", str(s.karyotypes)),
    ]


def normalize_tab(name: str) -> str:
    """A valid tab name (case-insensitive); the first tab when unknown."""
    wanted = str(name or "").strip().lower()
    return next((t for t in TAB_NAMES if t.lower() == wanted), TAB_NAMES[0])
