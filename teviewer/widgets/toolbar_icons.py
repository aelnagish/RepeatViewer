"""Toolbar icons: bundled assets with Qt fallback."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QIcon, QPainter, QPixmap
from PySide6.QtWidgets import QStyle, QWidget

_DEFAULT_ICON_SIZE = 22
_ALIVIEW_ICON_SIZE = 40

# teviewer/assets/icons (preferred) and legacy backends/inco overrides
_ASSETS_DIR = Path(__file__).resolve().parent.parent / "assets" / "icons"
_LEGACY_DIR = Path(__file__).resolve().parent.parent / "backends" / "inco"

_NAME_FILES: dict[str, tuple[str, ...]] = {
    "aliview": ("aliview.png", "aliview.svg"),
    "align": ("align.svg",),
    "blast": ("blast.svg",),
    "te_genome_map": ("blast.svg",),
    "extend": ("extend.svg",),
    "teaid": ("teaid.svg", "TEAid.png"),
    "msa_clean": ("msa_clean.svg",),
    "consensus": ("consensus.svg",),
    "flexidot": ("flexidot.svg", "flexidot.png"),
}

_FALLBACK_PIXMAP: dict[str, QStyle.StandardPixmap] = {
    "aliview": QStyle.StandardPixmap.SP_FileDialogDetailedView,
    "align": QStyle.StandardPixmap.SP_BrowserReload,
    "blast": QStyle.StandardPixmap.SP_ArrowForward,
    "te_genome_map": QStyle.StandardPixmap.SP_FileDialogListView,
    "extend": QStyle.StandardPixmap.SP_ArrowDown,
    "teaid": QStyle.StandardPixmap.SP_FileDialogInfoView,
    "msa_clean": QStyle.StandardPixmap.SP_DialogResetButton,
    "consensus": QStyle.StandardPixmap.SP_DialogApplyButton,
    "flexidot": QStyle.StandardPixmap.SP_DesktopIcon,
}


def _resolve_icon_path(name: str) -> Path | None:
    for filename in _NAME_FILES.get(name, ()):
        for folder in (_ASSETS_DIR, _LEGACY_DIR):
            path = folder / filename
            if path.is_file():
                return path
    return None


def _pixmap_on_canvas(source: QPixmap, size: int) -> QPixmap:
    """Scale into a square canvas so wide icons (e.g. old TEAid PNG) fit the toolbar."""
    canvas = QPixmap(size, size)
    canvas.fill(Qt.GlobalColor.transparent)
    scaled = source.scaled(
        size,
        size,
        Qt.AspectRatioMode.KeepAspectRatio,
        Qt.TransformationMode.SmoothTransformation,
    )
    x = (size - scaled.width()) // 2
    y = (size - scaled.height()) // 2
    painter = QPainter(canvas)
    painter.drawPixmap(x, y, scaled)
    painter.end()
    return canvas


def function_icon(
    widget: QWidget,
    name: str,
    *,
    size: int | None = None,
) -> QIcon:
    """Return toolbar icon for a function button."""
    icon_size = size if size is not None else _DEFAULT_ICON_SIZE
    path = _resolve_icon_path(name)
    if path is not None:
        icon = QIcon(str(path))
        if not icon.isNull():
            pix = icon.pixmap(icon_size, icon_size)
            if not pix.isNull():
                return QIcon(_pixmap_on_canvas(pix, icon_size))
            return icon

    style = widget.style()
    pixmap = _FALLBACK_PIXMAP.get(name, QStyle.StandardPixmap.SP_FileIcon)
    return style.standardIcon(pixmap)


def aliview_icon(widget: QWidget) -> QIcon:
    """Larger AliView toolbar icon."""
    return function_icon(widget, "aliview", size=_ALIVIEW_ICON_SIZE)
