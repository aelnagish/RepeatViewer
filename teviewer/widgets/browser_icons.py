"""Browser icons loaded from bundled SVG assets."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QIcon, QPalette, QPixmap
from PySide6.QtWidgets import QApplication, QWidget

from teviewer.widgets.svg_icons import themed_svg_icon_modes

_BROWSER_ICON_SIZE = 16
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


def clear_browser_icon_cache() -> None:
    _ICON_CACHE.clear()


def _is_dark_palette() -> bool:
    app = QApplication.instance()
    if app is None:
        return False
    return app.palette().color(QPalette.ColorRole.Window).lightness() < 128


def _browser_modes() -> list[tuple[QIcon.Mode, QColor]]:
    if _is_dark_palette():
        normal = QColor("#c4c7c5")
        selected = QColor("#ffffff")
    else:
        normal = QColor("#5f6368")
        selected = QColor("#ffffff")
    return [
        (QIcon.Mode.Normal, normal),
        (QIcon.Mode.Selected, selected),
        (QIcon.Mode.Active, selected),
    ]


def _browser_icon(name: str, widget: QWidget | None) -> QIcon:
    theme_key = "dark" if _is_dark_palette() else "light"
    cache_key = f"{theme_key}:{name}"
    cached = _ICON_CACHE.get(cache_key)
    if cached is not None and not cached.isNull():
        return cached
    icon = themed_svg_icon_modes(
        name,
        _BROWSER_ICON_SIZE,
        widget=widget,
        subdir=_BROWSER_SUBDIR,
        mode_colors=_browser_modes(),
    )
    if not icon.isNull():
        _ICON_CACHE[cache_key] = icon
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
