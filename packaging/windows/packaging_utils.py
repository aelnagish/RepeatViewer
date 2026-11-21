#!/usr/bin/env python3
"""Build RepeatViewer.exe packaging helpers (Windows)."""

from __future__ import annotations

import os
import struct
import sys
from pathlib import Path

PACKAGING = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGING.parent))

from bundle_datas import app_version, collect_datas_for, load_app_icon_renderer  # noqa: E402

PKG = Path(__file__).resolve().parents[2] / "repeatviewer"
ASSETS = PKG / "assets" / "icons"


def collect_datas() -> list[tuple[str, str]]:
    """Datas for PyInstaller: Windows tools only (no mac/linux binaries)."""
    return collect_datas_for("windows")


def collect_conda_runtime_dlls() -> list[tuple[str, str]]:
    """Collect conda ``Library\\bin`` DLLs that PyInstaller often fails to resolve.

    Conda Python loads ``_ctypes`` from ``ffi.dll`` (and similar stdlib
    extensions from liblzma / libbz2 / libexpat / sqlite3). Those live next
    to conda, not on the default DLL search path, so a frozen app crashes with:

        ImportError: DLL load failed while importing _ctypes
    """
    names = (
        "ffi.dll",
        "libffi.dll",
        "libffi-8.dll",
        "liblzma.dll",
        "libbz2.dll",
        "LIBBZ2.dll",
        "libexpat.dll",
        "sqlite3.dll",
    )
    folders = (
        Path(sys.prefix) / "Library" / "bin",
        Path(sys.base_prefix) / "Library" / "bin",
    )
    seen: set[str] = set()
    binaries: list[tuple[str, str]] = []
    for folder in folders:
        if not folder.is_dir():
            continue
        for name in names:
            path = folder / name
            key = path.name.lower()
            if path.is_file() and key not in seen:
                seen.add(key)
                binaries.append((str(path), "."))
                print(f"INFO: bundling conda DLL {path}")
    return binaries


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

    from PySide6.QtCore import QByteArray, QBuffer, QIODevice
    from PySide6.QtGui import QGuiApplication

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    app = QGuiApplication.instance() or QGuiApplication([])

    pixmap = load_app_icon_renderer()(256)
    if pixmap is None:
        return None

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
