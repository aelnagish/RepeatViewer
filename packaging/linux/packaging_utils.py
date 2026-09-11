#!/usr/bin/env python3
"""Build RepeatViewer packaging helpers (Linux)."""

from __future__ import annotations

import os
import sys
from pathlib import Path

PACKAGING = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGING.parent))

from bundle_datas import app_version, collect_datas_for  # noqa: E402

PKG = Path(__file__).resolve().parents[2] / "repeatviewer"
ASSETS = PKG / "assets" / "icons"


def collect_datas() -> list[tuple[str, str]]:
    """Datas for PyInstaller: Linux tools only (no mac/windows binaries)."""
    return collect_datas_for("linux")


def make_png_icon(out_png: Path, size: int = 256) -> Path | None:
    """Render app_icon.svg to a PNG via Qt (for .desktop / PyInstaller icon)."""
    svg = ASSETS / "browser" / "app_icon.svg"
    if not svg.is_file():
        svg = ASSETS / "app_icon.svg"
    if not svg.is_file():
        return None

    from PySide6.QtCore import QByteArray, QRectF, Qt
    from PySide6.QtGui import QGuiApplication, QPainter, QPixmap
    from PySide6.QtSvg import QSvgRenderer

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    app = QGuiApplication.instance() or QGuiApplication([])

    renderer = QSvgRenderer(QByteArray(svg.read_bytes()))
    if not renderer.isValid():
        return None

    pixmap = QPixmap(size, size)
    pixmap.fill(Qt.GlobalColor.transparent)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    art = size * 0.8125
    margin = (size - art) / 2.0
    renderer.render(painter, QRectF(margin, margin, art, art))
    painter.end()

    out_png.parent.mkdir(parents=True, exist_ok=True)
    if out_png.exists():
        out_png.unlink()
    if not pixmap.save(str(out_png), "PNG"):
        return None
    del app
    return out_png if out_png.is_file() else None


if __name__ == "__main__":
    target = PACKAGING / "RepeatViewer.png"
    path = make_png_icon(target)
    print(path or "failed")
    sys.exit(0 if path else 1)
