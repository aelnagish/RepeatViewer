#!/usr/bin/env python3
"""Build RepeatViewer.exe packaging helpers (Windows)."""

from __future__ import annotations

import os
import struct
import sys
from pathlib import Path

PACKAGING = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGING.parent))

from bundle_datas import app_version, collect_datas_for  # noqa: E402

PKG = Path(__file__).resolve().parents[2] / "repeatviewer"
ASSETS = PKG / "assets" / "icons"


def collect_datas() -> list[tuple[str, str]]:
    """Datas for PyInstaller: Windows tools only (no mac/linux binaries)."""
    return collect_datas_for("windows")


def _png_to_ico(png_bytes: bytes, out_ico: Path, sizes: tuple[int, ...] = (16, 32, 48, 256)) -> None:
    """Write a multi-size .ico from a square RGBA PNG (Qt-rendered)."""
    from PySide6.QtCore import QByteArray, QBuffer, QIODevice
    from PySide6.QtGui import QImage

    images: list[tuple[int, bytes]] = []
    base = QImage.fromData(QByteArray(png_bytes), "PNG")
    if base.isNull():
        raise RuntimeError("failed to decode PNG for .ico")

    for size in sizes:
        scaled = base.scaled(size, size)
        buf = QBuffer()
        buf.open(QIODevice.OpenModeFlag.WriteOnly)
        scaled.save(buf, "PNG")
        images.append((size, bytes(buf.data())))

    count = len(images)
    offset = 6 + 16 * count
    entries = bytearray()
    payload = bytearray()
    for size, data in images:
        w = 0 if size >= 256 else size
        h = 0 if size >= 256 else size
        entries += struct.pack("<BBBBHHII", w, h, 0, 0, 1, 32, len(data), offset)
        payload += data
        offset += len(data)

    header = struct.pack("<HHH", 0, 1, count)
    out_ico.write_bytes(header + entries + payload)


def make_ico(out_ico: Path) -> Path | None:
    """Render app_icon.svg to a .ico via Qt (no extra deps)."""
    svg = ASSETS / "browser" / "app_icon.svg"
    if not svg.is_file():
        svg = ASSETS / "app_icon.svg"
    if not svg.is_file():
        return None

    from PySide6.QtCore import QByteArray, QBuffer, QIODevice, QRectF, Qt
    from PySide6.QtGui import QGuiApplication, QPainter, QPixmap
    from PySide6.QtSvg import QSvgRenderer

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    app = QGuiApplication.instance() or QGuiApplication([])

    renderer = QSvgRenderer(QByteArray(svg.read_bytes()))
    if not renderer.isValid():
        return None

    pixel = 256
    pixmap = QPixmap(pixel, pixel)
    pixmap.fill(Qt.GlobalColor.transparent)
    painter = QPainter(pixmap)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
    art = pixel * 0.8125
    margin = (pixel - art) / 2.0
    renderer.render(painter, QRectF(margin, margin, art, art))
    painter.end()

    buf_bytes = QByteArray()
    buffer = QBuffer(buf_bytes)
    buffer.open(QIODevice.OpenModeFlag.WriteOnly)
    pixmap.save(buffer, "PNG")
    png = bytes(buf_bytes)

    out_ico.parent.mkdir(parents=True, exist_ok=True)
    if out_ico.exists():
        out_ico.unlink()
    _png_to_ico(png, out_ico)
    del app
    return out_ico if out_ico.is_file() else None


if __name__ == "__main__":
    target = PACKAGING / "RepeatViewer.ico"
    path = make_ico(target)
    print(path or "failed")
    sys.exit(0 if path else 1)
