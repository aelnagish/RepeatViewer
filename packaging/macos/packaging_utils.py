#!/usr/bin/env python3
"""Build RepeatViewer.app / .dmg packaging helpers."""

from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
PKG = ROOT / "repeatviewer"
BACKENDS = PKG / "backends"
ASSETS = PKG / "assets" / "icons"
PACKAGING = Path(__file__).resolve().parent


def app_version() -> str:
    init = (PKG / "__init__.py").read_text(encoding="utf-8")
    for line in init.splitlines():
        if line.startswith("__version__"):
            return line.split("=", 1)[1].strip().strip("\"'")
    return "0.0.0"


def collect_datas() -> list[tuple[str, str]]:
    """Datas for PyInstaller: macOS tools + icons (skip linux/windows binaries)."""
    datas: list[tuple[str, str]] = []

    blast_dst = "repeatviewer/backends/blast"
    for name in (
        "NOTICE.md",
        "blastn_mac",
        "blastp_mac",
        "blastx_mac",
        "makeblastdb_mac",
        "makeprofiledb_mac",
        "rpstblastn_mac",
    ):
        path = BACKENDS / "blast" / name
        if path.is_file():
            datas.append((str(path), blast_dst))

    mafft_macos = BACKENDS / "mafft" / "macos"
    if mafft_macos.is_dir():
        datas.append((str(mafft_macos), "repeatviewer/backends/mafft/macos"))
    mafft_notice = BACKENDS / "mafft" / "NOTICE.md"
    if mafft_notice.is_file():
        datas.append((str(mafft_notice), "repeatviewer/backends/mafft"))

    muscle_macos = BACKENDS / "muscle" / "macos"
    if muscle_macos.is_dir():
        datas.append((str(muscle_macos), "repeatviewer/backends/muscle/macos"))
    muscle_notice = BACKENDS / "muscle" / "NOTICE.md"
    if muscle_notice.is_file():
        datas.append((str(muscle_notice), "repeatviewer/backends/muscle"))

    aliview = BACKENDS / "aliview"
    if aliview.is_dir():
        datas.append((str(aliview), "repeatviewer/backends/aliview"))

    icons = ASSETS
    if icons.is_dir():
        datas.append((str(icons), "repeatviewer/assets/icons"))

    return datas


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

    # QGuiApplication needed for offscreen pixmaps on some Qt builds.
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

    # Apple iconset naming: icon_16x16.png, icon_16x16@2x.png, …
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
            # Leave ~18% margin like macOS dock icons.
            art = pixel * 0.8125
            margin = (pixel - art) / 2.0
            renderer.render(painter, QRectF(margin, margin, art, art))
            painter.end()
            name = f"icon_{size}x{size}{suffix}.png"
            pixmap.save(str(iconset / name), "PNG")

    # Also required 512@2x = 1024
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
