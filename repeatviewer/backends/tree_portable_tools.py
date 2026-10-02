"""Bundled tree programs shipped with RepeatViewer (like MAFFT).

**IQ-TREE 3.1.4** binaries live under ``backends/iqtree/{macos,linux,windows}/bin/``:

* macOS universal (``iqtree-3.1.4-macOS.zip``)
* Linux universal (``iqtree-3.1.4-Linux.tar.gz`` — launcher + intel/arm builds)
* Windows 64-bit (``iqtree-3.1.4-Windows.zip``)

``engine_executable`` prefers these unless the user ticks *Use … from the system PATH*.

FastME publishes no ready multi-platform binaries (source / conda only).
FastTree publishes Windows and Linux x86-64 builds on the authors' site, but not
macOS, and that host is not reliable for vendoring — so FastTree is not bundled.
"""

from __future__ import annotations

import os
import platform
import sys
from pathlib import Path

ENGINE_IQTREE = "iqtree"
ENGINE_FASTME = "fastme"
ENGINE_FASTTREE = "fasttree"
BUNDLED_ENGINES = (ENGINE_IQTREE,)  # only IQ-TREE ships ready binaries we vendor
IQTREE_VERSION = "3.1.4"

ENGINE_TITLES = {
    ENGINE_IQTREE: "IQ-TREE",
    ENGINE_FASTME: "FastME",
    ENGINE_FASTTREE: "FastTree",
}

# Kept for older tests / call sites that still import these names.
PORTABLE_ENGINES = BUNDLED_ENGINES


def _is_win() -> bool:
    return sys.platform == "win32"


def _runnable(path: Path) -> bool:
    return path.is_file() and (_is_win() or os.access(path, os.X_OK))


def iqtree_bundle_root() -> Path:
    """Root of the vendored IQ-TREE trees (PyInstaller-aware)."""
    if hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS) / "repeatviewer" / "backends" / "iqtree"
    return Path(__file__).resolve().parent / "iqtree"


def _platform_iqtree_dir(
    os_type: str | None = None,
    machine: str | None = None,
) -> Path | None:
    system = os_type or platform.system()
    root = iqtree_bundle_root()
    if system == "Darwin":
        return root / "macos"
    if system == "Windows":
        arch = (machine or platform.machine()).lower()
        if arch not in {"amd64", "x86_64"}:
            return None
        return root / "windows"
    if system == "Linux":
        return root / "linux"
    return None


def bundled_iqtree_executable(
    os_type: str | None = None,
    machine: str | None = None,
) -> str | None:
    """Path to the bundled IQ-TREE for this (or a given) platform, or None."""
    base = _platform_iqtree_dir(os_type, machine)
    if base is None:
        return None
    system = os_type or platform.system()
    if system == "Windows":
        names = ("iqtree3.exe", "iqtree2.exe", "iqtree.exe")
    else:
        names = ("iqtree3", "iqtree2", "iqtree")
    for folder in (base / "bin", base):
        for name in names:
            path = folder / name
            if _runnable(path):
                return str(path)
    return None


def bundled_engine_executable(
    engine: str,
    *,
    os_type: str | None = None,
    machine: str | None = None,
) -> str | None:
    """Bundled program for ``engine``, or None when that engine is not vendored."""
    if engine == ENGINE_IQTREE:
        return bundled_iqtree_executable(os_type, machine)
    return None


# --------------------------------------------------------------------------
# Back-compat aliases used by older dialog / test code
# --------------------------------------------------------------------------


class PortableToolError(RuntimeError):
    """Legacy name; download/build helpers were removed in favour of bundling."""


def portable_tree_tools_dir() -> Path:
    """Deprecated: previously the download target; now the IQ-TREE bundle root."""
    return iqtree_bundle_root()


def portable_engine_executable(engine: str) -> str | None:
    """Deprecated alias for :func:`bundled_engine_executable`."""
    return bundled_engine_executable(engine)


def missing_portable_engines() -> list[str]:
    """Engines that should be bundled but are missing on this machine."""
    return [e for e in BUNDLED_ENGINES if bundled_engine_executable(e) is None]


def download_portable_engine(engine: str, progress_cb=None) -> Path:  # noqa: ANN001
    """Removed: binaries are bundled. Raises so old UI cannot silently no-op."""
    raise PortableToolError(
        f"{ENGINE_TITLES.get(engine, engine)} is no longer downloaded. "
        "IQ-TREE is bundled with RepeatViewer; FastME and FastTree need conda or PATH."
    )
