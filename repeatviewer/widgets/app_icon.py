"""Application icon loaded from bundled SVG (macOS-style rounded square)."""

from __future__ import annotations

from functools import lru_cache
from pathlib import Path

from PySide6.QtCore import QByteArray, Qt, QRectF
from PySide6.QtGui import QIcon, QPainter, QPixmap
from PySide6.QtSvg import QSvgRenderer

_ASSETS_ROOT = Path(__file__).resolve().parent.parent / "assets" / "icons"

# macOS dock icons leave ~18% margin around artwork (1024 canvas → ~832pt glyph).
_MACOS_ICON_ARTWORK_SCALE = 0.8125


def _resolve_app_icon_path() -> Path | None:
    for path in (
        _ASSETS_ROOT / "app_icon.svg",
        _ASSETS_ROOT / "browser" / "app_icon.svg",
    ):
        if path.is_file():
            return path
    return None


_MACOS_ICON_SIZES = (16, 32, 64, 128, 256, 512, 1024)


def _render_svg_pixmap(svg_path: Path, logical_size: int, *, device_ratio: float = 1.0) -> QPixmap | None:
    svg_bytes = QByteArray(svg_path.read_bytes())
    renderer = QSvgRenderer(svg_bytes)
    if not renderer.isValid():
        return None

    ratio = device_ratio if device_ratio > 0 else 1.0
    physical = max(1, int(round(logical_size * ratio)))
    pixmap = QPixmap(physical, physical)
    pixmap.setDevicePixelRatio(ratio)
    pixmap.fill(Qt.GlobalColor.transparent)

    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    artwork = logical_size * _MACOS_ICON_ARTWORK_SCALE
    margin = (logical_size - artwork) / 2.0
    renderer.render(painter, QRectF(margin, margin, artwork, artwork))
    painter.end()
    return pixmap


@lru_cache(maxsize=1)
def app_icon() -> QIcon:
    """Return the RepeatViewer application icon at macOS-friendly sizes."""
    svg_path = _resolve_app_icon_path()
    if svg_path is None:
        return QIcon()

    icon = QIcon()
    for size in _MACOS_ICON_SIZES:
        for ratio in (1.0, 2.0):
            pixmap = _render_svg_pixmap(svg_path, size, device_ratio=ratio)
            if pixmap is not None:
                icon.addPixmap(pixmap)
    return icon
