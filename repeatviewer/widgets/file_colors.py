"""User-assigned colors for files in the Files list."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QSettings
from PySide6.QtGui import QBrush, QColor
from PySide6.QtWidgets import QTableWidgetItem

from repeatviewer.widgets.folder_colors import FOLDER_COLOR_PRESETS

FILE_COLORS_GROUP = "files"
FILE_COLORS_KEY = "colors_json"

# Same palette as folders so the UI feels consistent.
FILE_COLOR_PRESETS = FOLDER_COLOR_PRESETS


def _norm(path: Path) -> str:
    return str(path.expanduser().resolve())


def load_file_colors() -> dict[str, str]:
    s = QSettings()
    s.beginGroup(FILE_COLORS_GROUP)
    raw = s.value(FILE_COLORS_KEY, "{}")
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


def save_file_colors(colors: dict[str, str]) -> None:
    s = QSettings()
    s.beginGroup(FILE_COLORS_GROUP)
    s.setValue(FILE_COLORS_KEY, json.dumps(colors, indent=2, sort_keys=True))
    s.endGroup()


def get_file_color(path: Path) -> str | None:
    return load_file_colors().get(_norm(path))


def set_file_color(path: Path, color: str | None) -> None:
    file_key = _norm(path)
    colors = load_file_colors()
    if color:
        colors[file_key] = color
    else:
        colors.pop(file_key, None)
    save_file_colors(colors)


def replace_file_color_path(old: Path, new: Path) -> None:
    """Move a stored color when a file is renamed or moved on disk."""
    old_key = _norm(old)
    new_key = _norm(new)
    colors = load_file_colors()
    color = colors.pop(old_key, None)
    if color is None:
        return
    colors[new_key] = color
    save_file_colors(colors)


def replace_file_color_paths(moves: list[tuple[Path, Path]]) -> None:
    """Move stored colors for several renames at once (swaps and chains included)."""
    colors = load_file_colors()
    picked: list[tuple[str, str]] = []
    for old, new in moves:
        color = colors.pop(_norm(old), None)
        if color is not None:
            picked.append((_norm(new), color))
    if not picked:
        return
    for key, color in picked:
        colors[key] = color
    save_file_colors(colors)


def remove_file_color(path: Path) -> None:
    set_file_color(path, None)


def color_label(hex_color: str) -> str:
    """Human-readable name for a preset, otherwise the hex value."""
    needle = hex_color.strip().lower()
    for label, preset in FILE_COLOR_PRESETS:
        if preset.lower() == needle:
            return label
    return hex_color.strip().upper()


def apply_name_foreground(name_item: QTableWidgetItem | None, path: Path) -> None:
    """Set or clear the Name cell text color from stored file color."""
    if name_item is None:
        return
    color = get_file_color(path)
    if color:
        name_item.setForeground(QBrush(QColor(color)))
    else:
        name_item.setForeground(QBrush())
