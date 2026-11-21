#!/usr/bin/env python3
"""Build RepeatViewer packaging helpers (Linux)."""

from __future__ import annotations

import os
import sys
from pathlib import Path

PACKAGING = Path(__file__).resolve().parent
sys.path.insert(0, str(PACKAGING.parent))

from bundle_datas import app_version, collect_datas_for, load_app_icon_renderer  # noqa: E402

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

    from PySide6.QtGui import QGuiApplication

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    app = QGuiApplication.instance() or QGuiApplication([])

    pixmap = load_app_icon_renderer()(size)
    if pixmap is None:
        return None

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
