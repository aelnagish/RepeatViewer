"""Top-level folder roots hidden from the tree (not deleted on disk)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings

FOLDERS_SETTINGS_GROUP = "folders"
HIDDEN_ROOTS_KEY = "hidden_roots"


def _norm(path: Path) -> str:
    return str(path.expanduser().resolve())


def load_hidden_root_folders() -> set[str]:
    s = QSettings()
    s.beginGroup(FOLDERS_SETTINGS_GROUP)
    raw = s.value(HIDDEN_ROOTS_KEY, [])
    s.endGroup()
    if raw is None:
        return set()
    if isinstance(raw, str):
        raw = [raw]
    return {str(p) for p in raw if p}


def save_hidden_root_folders(paths: set[str]) -> None:
    s = QSettings()
    s.beginGroup(FOLDERS_SETTINGS_GROUP)
    s.setValue(HIDDEN_ROOTS_KEY, sorted(paths))
    s.endGroup()


def hide_root_folder(hidden: set[str], folder: Path) -> set[str]:
    hidden = set(hidden)
    hidden.add(_norm(folder))
    save_hidden_root_folders(hidden)
    return hidden


def unhide_root_folder(hidden: set[str], folder: Path) -> set[str]:
    hidden = set(hidden)
    hidden.discard(_norm(folder))
    save_hidden_root_folders(hidden)
    return hidden


def remove_hidden_root(hidden: set[str], folder: Path) -> set[str]:
    hidden = set(hidden)
    hidden.discard(_norm(folder))
    save_hidden_root_folders(hidden)
    return hidden


def replace_hidden_root_path(hidden: set[str], old: Path, new: Path) -> set[str]:
    old_key = _norm(old)
    if old_key not in hidden:
        return hidden
    hidden = set(hidden)
    hidden.discard(old_key)
    hidden.add(_norm(new))
    save_hidden_root_folders(hidden)
    return hidden
