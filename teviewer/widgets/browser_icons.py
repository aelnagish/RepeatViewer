"""Browser icons loaded from bundled SVG assets."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QIcon, QPixmap
from PySide6.QtWidgets import QWidget

from teviewer.widgets.svg_icons import themed_svg_icon_modes

_BROWSER_ICON_SIZE = 16
_BROWSER_GRAY = QColor("#5f6368")
_BROWSER_WHITE = QColor("#ffffff")
_BROWSER_MODES = [
    (QIcon.Mode.Normal, _BROWSER_GRAY),
    (QIcon.Mode.Selected, _BROWSER_WHITE),
    (QIcon.Mode.Active, _BROWSER_WHITE),
]
_BROWSER_SUBDIR = "browser"

_IMAGE_SUFFIXES = {
    ".png",
    ".jpg",
    ".jpeg",
    ".gif",
    ".webp",
    ".bmp",
    ".tif",
    ".tiff",
    ".svg",
}
_FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".fna", ".ffn", ".faa"}

# Reuse rendered icons — building SVG pixmaps per row is expensive in large folders.
_ICON_CACHE: dict[str, QIcon] = {}


def _browser_icon(name: str, widget: QWidget | None) -> QIcon:
    cached = _ICON_CACHE.get(name)
    if cached is not None and not cached.isNull():
        return cached
    icon = themed_svg_icon_modes(
        name,
        _BROWSER_ICON_SIZE,
        widget=widget,
        subdir=_BROWSER_SUBDIR,
        mode_colors=_BROWSER_MODES,
    )
    if not icon.isNull():
        _ICON_CACHE[name] = icon
    return icon


def browser_folder_icon(widget: QWidget | None = None, *, open_: bool = False) -> QIcon:
    return _browser_icon("chevron-down" if open_ else "chevron-right", widget)


def browser_folder_spacer_icon(widget: QWidget | None = None) -> QIcon:
    """Transparent icon matching chevron size so leaf folder text stays aligned."""
    del widget  # same logical size regardless of widget
    cached = _ICON_CACHE.get("__spacer__")
    if cached is not None and not cached.isNull():
        return cached
    pixmap = QPixmap(_BROWSER_ICON_SIZE, _BROWSER_ICON_SIZE)
    pixmap.fill(Qt.GlobalColor.transparent)
    icon = QIcon(pixmap)
    _ICON_CACHE["__spacer__"] = icon
    return icon


def browser_file_icon(widget: QWidget | None, path: Path) -> QIcon:
    suffix = path.suffix.lower()
    if suffix in _FASTA_SUFFIXES or path.name.lower().endswith(".fa.gz"):
        name = "sequence"
    elif suffix in _IMAGE_SUFFIXES:
        name = "image"
    else:
        name = "file"
    return _browser_icon(name, widget)
