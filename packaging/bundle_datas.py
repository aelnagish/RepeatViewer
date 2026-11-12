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
    manuals = PKG / "assets" / "manuals"
    if manuals.is_dir():
        datas.append((str(manuals), "repeatviewer/assets/manuals"))


def _add_database(datas: list[tuple[str, str]]) -> None:
    """Ship bundled motif / reference FASTA under ``repeatviewer/database/``.

    Required for PBS feature detection (Neumann2019 motif FASTA) in packaged
    apps; without this the .app looks for the file and every sequence is skipped.
    """
    _add_tree(datas, PKG / "database", "repeatviewer/database")


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
        _add_tree(datas, BACKENDS / "iqtree" / "macos", "repeatviewer/backends/iqtree/macos")
        _add_file(datas, BACKENDS / "iqtree" / "NOTICE.md", "repeatviewer/backends/iqtree")
        # minimap2: ship only the host CPU bottle
        mm_arch = "arm64" if cpu == "aarch64" else "x86_64"
        _add_tree(
            datas,
            BACKENDS / "minimap2" / "macos" / mm_arch,
            f"repeatviewer/backends/minimap2/macos/{mm_arch}",
        )
        _add_file(datas, BACKENDS / "minimap2" / "NOTICE.md", "repeatviewer/backends/minimap2")
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
        # IQ-TREE Linux universal covers intel + arm via launcher.
        _add_tree(datas, BACKENDS / "iqtree" / "linux", "repeatviewer/backends/iqtree/linux")
        _add_file(datas, BACKENDS / "iqtree" / "NOTICE.md", "repeatviewer/backends/iqtree")
        # Official minimap2 Linux binary is x86_64 only.
        if cpu == "x86_64":
            _add_tree(
                datas,
                BACKENDS / "minimap2" / "linux",
                "repeatviewer/backends/minimap2/linux",
            )
            _add_file(
                datas, BACKENDS / "minimap2" / "NOTICE.md", "repeatviewer/backends/minimap2"
            )
        else:
            print(
                f"WARNING: bundled Linux minimap2 is x86_64; "
                f"skipping it for host CPU {cpu!r}."
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
            BACKENDS / "iqtree" / "windows",
            "repeatviewer/backends/iqtree/windows",
        )
        _add_file(datas, BACKENDS / "iqtree" / "NOTICE.md", "repeatviewer/backends/iqtree")
        _add_tree(
            datas,
            BACKENDS / "minimap2" / "windows",
            "repeatviewer/backends/minimap2/windows",
        )
        _add_file(datas, BACKENDS / "minimap2" / "NOTICE.md", "repeatviewer/backends/minimap2")
        _add_tree(
            datas,
            BACKENDS / "muscle" / "windows",
            "repeatviewer/backends/muscle/windows",
        )
        _add_file(datas, BACKENDS / "muscle" / "NOTICE.md", "repeatviewer/backends/muscle")

    _add_aliview(datas)
    _add_icons(datas)
    _add_database(datas)
    return datas


def collect_annotation_python_deps() -> tuple[
    list[tuple[str, str]],
    list,
    list[str],
]:
    """Collect orfipy + pyhmmer + primer3 (+ pysam on macOS/Linux) for PyInstaller.

    Binary extensions alone are not enough: a ``pyhmmer/`` folder of ``.so`` files
    shadows the pure-Python package in PYZ and ``import pyhmmer`` fails at runtime.
    ``collect_all`` also misses some Cython extensions (e.g. ``pyhmmer.platform.*``)
    and orfipy's ``ahocorasick`` dependency — those are added explicitly.
    ``pysam`` enables large indexed BAM browsing; skip on Windows (no wheels).
    """
    import sys

    from PyInstaller.utils.hooks import collect_all

    datas: list[tuple[str, str]] = []
    binaries: list = []
    hidden: list[str] = ["orfipy_core", "ahocorasick", "pyhmmer.platform",
                         # Tools > Primers (primer3-py C extensions)
                         "primer3.thermoanalysis", "primer3.p3helpers"]

    def _without_tests(entries: list[tuple[str, str]]) -> list[tuple[str, str]]:
        kept: list[tuple[str, str]] = []
        for src, dest in entries:
            if "tests" in Path(src).parts or "tests" in Path(dest).parts:
                continue
            kept.append((src, dest))
        return kept

    pkgs = ["pyhmmer", "orfipy", "primer3"]
    if sys.platform != "win32":
        pkgs.append("pysam")
        hidden.append("pysam")
    for pkg in pkgs:
        try:
            pkg_datas, pkg_binaries, pkg_hidden = collect_all(pkg)
        except Exception as exc:  # noqa: BLE001
            print(f"WARNING: collect_all({pkg}) failed: {exc}")
            continue
        datas += _without_tests(list(pkg_datas))
        binaries += list(pkg_binaries)
        hidden += list(pkg_hidden)

    # Explicit Cython / extension modules that Analysis often skips.
    try:
        import pyhmmer  # type: ignore[import-untyped]

        pyhmmer_root = Path(pyhmmer.__file__).resolve().parent
        site_packages = pyhmmer_root.parent
        for so in pyhmmer_root.rglob("*.so"):
            dest = str(so.relative_to(site_packages).parent).replace("\\", "/")
            binaries.append((str(so), dest))
            # e.g. pyhmmer/platform/linux*.so → pyhmmer.platform.linux
            rel = so.relative_to(pyhmmer_root)
            mod_parts = list(rel.parts[:-1]) + [so.name.split(".", 1)[0]]
            if mod_parts:
                hidden.append("pyhmmer." + ".".join(mod_parts))
        for pyd in pyhmmer_root.rglob("*.pyd"):
            dest = str(pyd.relative_to(site_packages).parent).replace("\\", "/")
            binaries.append((str(pyd), dest))
            rel = pyd.relative_to(pyhmmer_root)
            mod_parts = list(rel.parts[:-1]) + [pyd.name.split(".", 1)[0]]
            if mod_parts:
                hidden.append("pyhmmer." + ".".join(mod_parts))
    except Exception as exc:  # noqa: BLE001
        print(f"WARNING: pyhmmer extension collect failed: {exc}")
        # Fallbacks if import failed before discovery (build host still needs them).
        hidden += [
            "pyhmmer.platform.bsd",
            "pyhmmer.platform.linux",
            "pyhmmer.platform.win32",
        ]

    try:
        import ahocorasick  # type: ignore[import-untyped]

        path = getattr(ahocorasick, "__file__", None)
        if path:
            binaries.append((str(path), "."))
    except Exception as exc:  # noqa: BLE001
        print(f"WARNING: ahocorasick collect failed: {exc}")

    try:
        import orfipy_core  # type: ignore[import-untyped]

        path = getattr(orfipy_core, "__file__", None)
        if path:
            binaries.append((str(path), "."))
    except Exception as exc:  # noqa: BLE001
        print(f"WARNING: orfipy_core collect failed: {exc}")

    # De-dupe while preserving order.
    seen: set[str] = set()
    unique_hidden: list[str] = []
    for name in hidden:
        if name in seen:
            continue
        seen.add(name)
        unique_hidden.append(name)

    return datas, binaries, unique_hidden
