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

from bundle_datas import app_version, collect_datas_for, load_app_icon_renderer  # noqa: E402

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

    from PySide6.QtGui import QGuiApplication

    os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
    app = QGuiApplication.instance() or QGuiApplication([])
    render_icon = load_app_icon_renderer()

    iconset = out_icns.with_suffix(".iconset")
    if iconset.exists():
        shutil.rmtree(iconset)
    iconset.mkdir(parents=True)

    for size in (16, 32, 128, 256, 512):
        for scale, suffix in ((1, ""), (2, "@2x")):
            pixmap = render_icon(size * scale)
            if pixmap is None:
                return None
            pixmap.save(str(iconset / f"icon_{size}x{size}{suffix}.png"), "PNG")

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
