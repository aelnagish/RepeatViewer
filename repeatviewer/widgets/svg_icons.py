"""Flat monochrome SVG icons (Tabler-style, MIT) with palette-aware tinting."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QByteArray, Qt, QRectF
from PySide6.QtGui import QColor, QIcon, QPainter, QPixmap
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import QWidget

# Tabler Icons (https://tabler.io/icons) — MIT license
_ASSETS_ROOT = Path(__file__).resolve().parent.parent / "assets" / "icons"


def _resolve_svg_path(name: str, *, subdir: str = "") -> Path | None:
    candidates: list[Path] = []
    if subdir:
        candidates.append(_ASSETS_ROOT / subdir / f"{name}.svg")
    candidates.append(_ASSETS_ROOT / f"{name}.svg")
    for path in candidates:
        if path.is_file():
            return path
    return None


def _device_pixel_ratio(widget: QWidget | None) -> float:
    if widget is None:
        return 1.0
    ratio = widget.devicePixelRatioF()
    return ratio if ratio > 0 else 1.0


def _tint_svg(svg_text: str, color: QColor) -> str:
    tinted = svg_text.replace("currentColor", color.name())
    for token in ("#000", "#000000", "black"):
        tinted = tinted.replace(f'stroke="{token}"', f'stroke="{color.name()}"')
        tinted = tinted.replace(f"stroke='{token}'", f"stroke='{color.name()}'")
        tinted = tinted.replace(f'fill="{token}"', f'fill="{color.name()}"')
    return tinted


def _render_themed_pixmap(
    name: str,
    size: int,
    *,
    color: QColor,
    widget: QWidget | None,
    subdir: str,
) -> QPixmap | None:
    path = _resolve_svg_path(name, subdir=subdir)
    if path is None:
        return None
    svg_text = _tint_svg(path.read_text(encoding="utf-8"), color)
    renderer = QSvgRenderer(QByteArray(svg_text.encode("utf-8")))
    if not renderer.isValid():
        return None

    ratio = _device_pixel_ratio(widget)
    physical = max(1, int(round(size * ratio)))
    pixmap = QPixmap(physical, physical)
    pixmap.setDevicePixelRatio(ratio)
    pixmap.fill(Qt.GlobalColor.transparent)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    # Logical bounds — without this, Retina pixmaps only show a corner of the SVG.
    renderer.render(painter, QRectF(0, 0, size, size))
    painter.end()
    return pixmap


def themed_svg_icon_modes(
    name: str,
    size: int,
    *,
    widget: QWidget | None = None,
    subdir: str = "",
    mode_colors: list[tuple[QIcon.Mode, QColor]],
) -> QIcon:
    """SVG icon with explicit colors per QIcon mode (e.g. gray normal, white selected)."""
    icon = QIcon()
    for mode, color in mode_colors:
        pixmap = _render_themed_pixmap(
            name,
            size,
            color=color,
            widget=widget,
            subdir=subdir,
        )
        if pixmap is not None:
            icon.addPixmap(pixmap, mode, QIcon.State.Off)
    return icon
