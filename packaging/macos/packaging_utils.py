#!/usr/bin/env python3
"""Build RepeatViewer.app / .dmg packaging helpers."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

PACKAGING = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGING.parent))

from bundle_datas import app_version, collect_datas_for  # noqa: E402

PKG = Path(__file__).resolve().parents[2] / "repeatviewer"
ASSETS = PKG / "assets" / "icons"


def collect_datas() -> list[tuple[str, str]]:
    """Datas for PyInstaller: macOS tools only (no linux/windows binaries)."""
    return collect_datas_for("macos")


def make_icns(out_icns: Path) -> Path | None:
    """Render app_icon.svg to a .icns via Qt (no extra deps)."""
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

    iconset = out_icns.with_suffix(".iconset")
    if iconset.exists():
        shutil.rmtree(iconset)
    iconset.mkdir(parents=True)

    svg_bytes = QByteArray(svg.read_bytes())
    renderer = QSvgRenderer(svg_bytes)
    if not renderer.isValid():
        return None

    sizes = (16, 32, 64, 128, 256, 512)
    for size in sizes:
        for scale, suffix in ((1, ""), (2, "@2x")):
            pixel = size * scale
            if pixel > 1024:
                continue
            pixmap = QPixmap(pixel, pixel)
            pixmap.fill(Qt.GlobalColor.transparent)
            painter = QPainter(pixmap)
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
            art = pixel * 0.8125
            margin = (pixel - art) / 2.0
            renderer.render(painter, QRectF(margin, margin, art, art))
            painter.end()
            name = f"icon_{size}x{size}{suffix}.png"
            pixmap.save(str(iconset / name), "PNG")

    pixmap = QPixmap(1024, 1024)
    pixmap.fill(Qt.GlobalColor.transparent)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    art = 1024 * 0.8125
    margin = (1024 - art) / 2.0
    renderer.render(painter, QRectF(margin, margin, art, art))
    painter.end()
    pixmap.save(str(iconset / "icon_512x512@2x.png"), "PNG")

    out_icns.parent.mkdir(parents=True, exist_ok=True)
    if out_icns.exists():
        out_icns.unlink()
    subprocess.run(
        ["iconutil", "-c", "icns", str(iconset), "-o", str(out_icns)],
        check=True,
    )
    shutil.rmtree(iconset, ignore_errors=True)
    del app
    return out_icns if out_icns.is_file() else None


if __name__ == "__main__":
    target = PACKAGING / "RepeatViewer.icns"
    path = make_icns(target)
    print(path or "failed")
