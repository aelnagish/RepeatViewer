"""Bundled minimap2 2.31 shipped with RepeatViewer (like MAFFT / IQ-TREE).

Layout under ``backends/minimap2/``:

* ``linux/bin/minimap2`` — official lh3 x86_64 release
* ``windows/bin/minimap2.exe`` (+ DLLs) — win-ngs community UCRT64 build
* ``macos/arm64/minimap2``, ``macos/x86_64/minimap2`` — Homebrew bottles

All three are **2.31 (r1302)**. Tick *Use programs from the system PATH* to ignore
the bundle and use the user's own install.
"""

from __future__ import annotations

import os
import platform
import sys
from pathlib import Path

MINIMAP2_VERSION = "2.31"


def minimap2_bundle_root() -> Path:
    if hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS) / "repeatviewer" / "backends" / "minimap2"
    return Path(__file__).resolve().parent / "minimap2"


def _is_win() -> bool:
    return sys.platform == "win32"


def _bundled_executable_path(path: Path, target_system: str) -> bool:
    """True if the vendored binary is present (and executable on this host OS)."""
    if not path.is_file():
        return False
    if target_system != platform.system():
        return True
    return _is_win() or os.access(path, os.X_OK)


def _macos_arch(machine: str | None = None) -> str:
    m = (machine or platform.machine()).lower()
    return "arm64" if m in {"arm64", "aarch64"} else "x86_64"


def bundled_minimap2_executable(
    os_type: str | None = None,
    machine: str | None = None,
) -> str | None:
    """Path to the vendored minimap2 for this (or a given) platform, or None."""
    system = os_type or platform.system()
    root = minimap2_bundle_root()
    if system == "Darwin":
        path = root / "macos" / _macos_arch(machine) / "minimap2"
    elif system == "Windows":
        arch = (machine or platform.machine()).lower()
        if arch not in {"amd64", "x86_64"}:
            return None
        path = root / "windows" / "bin" / "minimap2.exe"
    elif system == "Linux":
        arch = (machine or platform.machine()).lower()
        if arch not in {"amd64", "x86_64"}:
            return None  # official release is x86_64 only
        path = root / "linux" / "bin" / "minimap2"
    else:
        return None
    return str(path) if _bundled_executable_path(path, system) else None
