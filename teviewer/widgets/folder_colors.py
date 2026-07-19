"""User-assigned colors for folders in the Folders tree."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QSettings, Qt
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import QTreeWidgetItem

FOLDER_COLORS_GROUP = "folders"
FOLDER_COLORS_KEY = "colors_json"

FOLDER_COLOR_PRESETS: tuple[tuple[str, str], ...] = (
    ("Blue", "#2563eb"),
    ("Teal", "#0d9488"),
    ("Green", "#16a34a"),
    ("Orange", "#ea580c"),
    ("Purple", "#7c3aed"),
    ("Brown", "#92400e"),
)


def _norm(path: Path) -> str:
    return str(path.expanduser().resolve())


def load_folder_colors() -> dict[str, str]:
    s = QSettings()
    s.beginGroup(FOLDER_COLORS_GROUP)
    raw = s.value(FOLDER_COLORS_KEY, "{}")
    s.endGroup()
    if not raw:
        return {}
    try:
        decoded = json.loads(str(raw))
    except (json.JSONDecodeError, TypeError):
        return {}
    if not isinstance(decoded, dict):
        return {}
    colors: dict[str, str] = {}
    for key, value in decoded.items():
        color = str(value).strip()
        if key and color.startswith("#") and len(color) >= 4:
            colors[str(key)] = color
    return colors


def save_folder_colors(colors: dict[str, str]) -> None:
    s = QSettings()
    s.beginGroup(FOLDER_COLORS_GROUP)
    s.setValue(FOLDER_COLORS_KEY, json.dumps(colors, indent=2, sort_keys=True))
    s.endGroup()


def get_folder_color(folder: Path) -> str | None:
    return load_folder_colors().get(_norm(folder))


def set_folder_color(folder: Path, color: str | None) -> None:
    folder_key = _norm(folder)
    colors = load_folder_colors()
    if color:
        colors[folder_key] = color
    else:
        colors.pop(folder_key, None)
    save_folder_colors(colors)


def replace_folder_color_paths(old: Path, new: Path) -> None:
    """Move stored colors when a folder is renamed on disk."""
    old_root = Path(_norm(old))
    new_root = Path(_norm(new))
    colors = load_folder_colors()
    if not colors:
        return
    updated: dict[str, str] = {}
    for path_key, color in colors.items():
        path = Path(path_key)
        try:
            rel = path.relative_to(old_root)
            updated[_norm(new_root / rel)] = color
        except ValueError:
            if path_key == _norm(old):
                updated[_norm(new)] = color
            else:
                updated[path_key] = color
    save_folder_colors(updated)


def remove_folder_colors_under(folder: Path) -> None:
    """Drop colors for a deleted folder and everything beneath it."""
    root = Path(_norm(folder))
    colors = load_folder_colors()
    if not colors:
        return
    kept: dict[str, str] = {}
    for path_key, color in colors.items():
        path = Path(path_key)
        if path == root:
            continue
        try:
            path.relative_to(root)
        except ValueError:
            kept[path_key] = color
    save_folder_colors(kept)


def apply_folder_colors_to_tree(tree) -> None:
    """Re-apply stored colors to every folder item in a tree widget."""
    from PySide6.QtWidgets import QTreeWidget, QTreeWidgetItem

    if not isinstance(tree, QTreeWidget):
        return

    def walk(item: QTreeWidgetItem) -> None:
        raw = item.data(0, Qt.ItemDataRole.UserRole)
        if raw:
            apply_folder_item_color(item, Path(str(raw)))
        for index in range(item.childCount()):
            walk(item.child(index))

    for index in range(tree.topLevelItemCount()):
        walk(tree.topLevelItem(index))


def apply_folder_item_color(item: QTreeWidgetItem, folder: Path) -> None:
    color = get_folder_color(folder)
    if color:
        item.setForeground(0, QBrush(QColor(color)))
    else:
        item.setForeground(0, QBrush())
