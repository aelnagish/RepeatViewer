"""Shared PyInstaller data collection: only ship binaries for the target OS/arch."""

from __future__ import annotations

import platform
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PKG = ROOT / "repeatviewer"
BACKENDS = PKG / "backends"
ASSETS = PKG / "assets" / "icons"

# Runtime AliView files only (skip installer / makeself helpers).
_ALIVIEW_RUNTIME = (
    "aliview.jar",
    "aliview",
    "aliicon_128x128.png",
    "README.txt",
)


def app_version() -> str:
    init = (PKG / "__init__.py").read_text(encoding="utf-8")
    for line in init.splitlines():
        if line.startswith("__version__"):
            return line.split("=", 1)[1].strip().strip("\"'")
    return "0.0.0"


def host_cpu() -> str:
    """Normalize machine to ``x86_64`` or ``aarch64`` when possible."""
    machine = platform.machine().lower()
    if machine in {"amd64", "x86_64"}:
        return "x86_64"
    if machine in {"arm64", "aarch64"}:
        return "aarch64"
    return machine


def _add_file(datas: list[tuple[str, str]], path: Path, dest_dir: str) -> None:
    if path.is_file():
        datas.append((str(path), dest_dir))


def _add_tree(datas: list[tuple[str, str]], path: Path, dest_dir: str) -> None:
    if path.is_dir():
        datas.append((str(path), dest_dir))


def _add_aliview(datas: list[tuple[str, str]]) -> None:
    base = BACKENDS / "aliview"
    dest = "repeatviewer/backends/aliview"
    if not base.is_dir():
        return
    for name in _ALIVIEW_RUNTIME:
        _add_file(datas, base / name, dest)


def _add_icons(datas: list[tuple[str, str]]) -> None:
    if ASSETS.is_dir():
        datas.append((str(ASSETS), "repeatviewer/assets/icons"))


def collect_datas_for(target: str) -> list[tuple[str, str]]:
    """Collect bundled tool datas for ``macos``, ``linux``, or ``windows``.

    Never includes other OS tool trees. On Linux, also skips wrong-CPU
    muscle/blast/mafft when those binaries are arch-specific.
    """
    target = target.lower().strip()
    if target not in {"macos", "linux", "windows"}:
        raise ValueError(f"Unknown packaging target: {target!r}")

    datas: list[tuple[str, str]] = []
    blast_dst = "repeatviewer/backends/blast"
    cpu = host_cpu()

    if target == "macos":
        for name in (
            "NOTICE.md",
            "blastn_mac",
            "blastp_mac",
            "blastx_mac",
            "makeblastdb_mac",
            "makeprofiledb_mac",
            "rpstblastn_mac",
        ):
            _add_file(datas, BACKENDS / "blast" / name, blast_dst)
        _add_tree(datas, BACKENDS / "mafft" / "macos", "repeatviewer/backends/mafft/macos")
        _add_file(datas, BACKENDS / "mafft" / "NOTICE.md", "repeatviewer/backends/mafft")
        # Bundled MUSCLE for macOS is arm64-only in this repo.
        muscle_arch = "arm64" if cpu == "aarch64" else None
        if muscle_arch:
            _add_tree(
                datas,
                BACKENDS / "muscle" / "macos" / muscle_arch,
                f"repeatviewer/backends/muscle/macos/{muscle_arch}",
            )
        _add_file(datas, BACKENDS / "muscle" / "NOTICE.md", "repeatviewer/backends/muscle")

    elif target == "linux":
        # Bundled BLAST + MAFFT are x86_64; omit on aarch64 hosts.
        if cpu == "x86_64":
            for name in (
                "NOTICE.md",
                "blastn_linux",
                "blastp_linux",
                "blastx_linux",
                "makeblastdb_linux",
                "makeprofiledb_linux",
                "rpstblastn_linux",
            ):
                _add_file(datas, BACKENDS / "blast" / name, blast_dst)
            _add_tree(
                datas,
                BACKENDS / "mafft" / "linux",
                "repeatviewer/backends/mafft/linux",
            )
            _add_file(datas, BACKENDS / "mafft" / "NOTICE.md", "repeatviewer/backends/mafft")
        else:
            print(
                f"WARNING: bundled Linux BLAST/MAFFT are x86_64; "
                f"skipping them for host CPU {cpu!r}. Use system tools or conda."
            )
        # MUSCLE: only the matching arch subdirectory (not both).
        if cpu in {"x86_64", "aarch64"}:
            _add_tree(
                datas,
                BACKENDS / "muscle" / "linux" / cpu,
                f"repeatviewer/backends/muscle/linux/{cpu}",
            )
        _add_file(datas, BACKENDS / "muscle" / "NOTICE.md", "repeatviewer/backends/muscle")

    else:  # windows
        for name in (
            "NOTICE.md",
            "blastn.exe",
            "blastp.exe",
            "blastx.exe",
            "makeblastdb.exe",
            "makeprofiledb.exe",
            "rpstblastn.exe",
        ):
            _add_file(datas, BACKENDS / "blast" / name, blast_dst)
        _add_tree(
            datas,
            BACKENDS / "mafft" / "windows",
            "repeatviewer/backends/mafft/windows",
        )
        _add_file(datas, BACKENDS / "mafft" / "NOTICE.md", "repeatviewer/backends/mafft")
        _add_tree(
            datas,
            BACKENDS / "muscle" / "windows",
            "repeatviewer/backends/muscle/windows",
        )
        _add_file(datas, BACKENDS / "muscle" / "NOTICE.md", "repeatviewer/backends/muscle")

    _add_aliview(datas)
    _add_icons(datas)
    return datas
