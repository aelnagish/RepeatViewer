"""Application icon loaded from bundled SVG (macOS-style rounded square)."""

from __future__ import annotations

import math
from functools import lru_cache
from pathlib import Path

from PySide6.QtCore import QByteArray, QRectF, Qt
from PySide6.QtGui import QIcon, QPainter, QPainterPath, QPixmap
from PySide6.QtSvg import QSvgRenderer

_ASSETS_ROOT = Path(__file__).resolve().parent.parent / "assets" / "icons"

# Plate size as a fraction of the icon canvas. Apple's grid is 824/1024 = 0.8125,
# but measuring a real Dock screenshot: neighbouring icons cover 90-93% of their
# 38px slot while this one covered 95.4%, because a rounded rect keeps more
# corner area than a squircle. Squircle + 0.80 would match Terminal exactly at
# 90.2%; 0.78 lands at ~86%, the extra trim offsetting the fact that a solid
# white plate reads larger than a dark one at equal size. Raise to 0.80 if this
# ever looks too small.
_MACOS_ICON_ARTWORK_SCALE = 0.78
# Superellipse exponent for the plate outline: |x|^n + |y|^n = 1. n=4 matches the
# macOS squircle; a plain rounded rectangle (which is what QPainterPath's
# addRoundedRect draws) keeps too much corner area and looks blockier.
_ICON_SQUIRCLE_EXPONENT = 4.0
# Points sampled along the outline; plenty at 1024px with antialiasing on.
_ICON_SQUIRCLE_SEGMENTS = 720


def _resolve_app_icon_path() -> Path | None:
    for path in (
        _ASSETS_ROOT / "app_icon.svg",
        _ASSETS_ROOT / "browser" / "app_icon.svg",
    ):
        if path.is_file():
            return path
    return None


_MACOS_ICON_SIZES = (16, 32, 64, 128, 256, 512, 1024)


def _squircle_path(rect: QRectF) -> QPainterPath:
    """Superellipse outline inscribed in *rect* (the macOS icon shape)."""
    cx, cy = rect.center().x(), rect.center().y()
    rx, ry = rect.width() / 2.0, rect.height() / 2.0
    power = 2.0 / _ICON_SQUIRCLE_EXPONENT

    path = QPainterPath()
    for index in range(_ICON_SQUIRCLE_SEGMENTS):
        angle = 2.0 * math.pi * index / _ICON_SQUIRCLE_SEGMENTS
        cos_a, sin_a = math.cos(angle), math.sin(angle)
        x = cx + rx * math.copysign(abs(cos_a) ** power, cos_a)
        y = cy + ry * math.copysign(abs(sin_a) ** power, sin_a)
        if index == 0:
            path.moveTo(x, y)
        else:
            path.lineTo(x, y)
    path.closeSubpath()
    return path


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
    plate = QRectF(margin, margin, artwork, artwork)

    # Clip to a squircle so the logo matches other Dock icons (not a sharp tile).
    painter.setClipPath(_squircle_path(plate))
    renderer.render(painter, plate)
    painter.end()
    return pixmap


def render_app_icon_pixmap(size: int) -> QPixmap | None:
    """Squircle app icon at *size* px, shared with the packaging scripts.

    The .icns/.ico/.png built for a release must use this same geometry, or the
    icon macOS shows before launch (from the bundle) differs from the one Qt
    puts in the Dock after launch, and Tahoe wraps the square one in a grey tile.
    """
    svg_path = _resolve_app_icon_path()
    if svg_path is None:
        return None
    return _render_svg_pixmap(svg_path, size)


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
