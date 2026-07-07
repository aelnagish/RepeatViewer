"""Resolve minimap2: bundled on Linux/Windows, Homebrew on macOS."""

from __future__ import annotations

import os
import platform
import shutil
from pathlib import Path


def _plugin_root() -> Path:
    return Path(__file__).resolve().parent


def _bundled_minimap2() -> Path | None:
    system = platform.system()
    if system == "Linux":
        candidate = _plugin_root() / "binaries" / "linux" / "minimap2"
    elif system == "Windows":
        candidate = _plugin_root() / "binaries" / "win" / "minimap2.exe"
    else:
        return None
    if candidate.is_file():
        return candidate
    return None


def minimap2_executable() -> str:
    system = platform.system()

    if system == "Darwin":
        found = shutil.which("minimap2")
        if found:
            return found
        raise FileNotFoundError(
            "minimap2 not found on macOS. Install with:\n  brew install minimap2"
        )

    bundled = _bundled_minimap2()
    if bundled is not None:
        if system == "Windows" or os.access(bundled, os.X_OK):
            return str(bundled)
        raise PermissionError(f"Bundled minimap2 is not executable: {bundled}")

    found = shutil.which("minimap2")
    if found:
        return found

    raise FileNotFoundError(
        "minimap2 not found. On Linux/Windows, use the bundled plugin binaries "
        "or add minimap2 to PATH."
    )
