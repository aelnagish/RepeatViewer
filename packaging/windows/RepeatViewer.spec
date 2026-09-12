# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for RepeatViewer (Windows).

Build on Windows (from repo root, with venv active):

    python packaging/windows/packaging_utils.py   # optional .ico
    pyinstaller packaging/windows/RepeatViewer.spec

Or run:

    powershell -ExecutionPolicy Bypass -File packaging/windows/build_windows.ps1
"""

from __future__ import annotations

import sys
from pathlib import Path

from PyInstaller.utils.hooks import collect_all, collect_data_files, collect_submodules

SPECDIR = Path(SPEC).resolve().parent  # noqa: F821 — injected by PyInstaller
ROOT = SPECDIR.parents[1]
sys.path.insert(0, str(SPECDIR))
sys.path.insert(0, str(SPECDIR.parent))

from packaging_utils import collect_datas  # noqa: E402
from bundle_datas import app_version  # noqa: E402

block_cipher = None
version = app_version()

datas = collect_datas()
binaries: list = []
hiddenimports: list[str] = [
    "flexidot",
    "flexidot.app",
    "flexidot.plotting",
    "flexidot.utils",
    "Bio",
    "Bio.SeqIO",
    "sklearn",
    "sklearn.cluster",
    "pytrf",
    "plotly",
    "matplotlib",
    "matplotlib.backends.backend_agg",
    "numpy",
    "numpy._core",
    "numpy._core._multiarray_umath",
    "six",
    "six.moves",
]

for pkg in ("flexidot",):
    try:
        pkg_datas, pkg_binaries, pkg_hidden = collect_all(pkg)
        datas += pkg_datas
        binaries += pkg_binaries
        hiddenimports += pkg_hidden
    except Exception as exc:  # noqa: BLE001
        print(f"WARNING: collect_all({pkg}) failed: {exc}")

# Plotly: data + submodules only. collect_all() can pull a second copy of
# numpy/scipy binaries and trigger:
#   ImportError: cannot load module more than once per process
try:
    datas += collect_data_files("plotly")
    hiddenimports += collect_submodules("plotly")
except Exception as exc:  # noqa: BLE001
    print(f"WARNING: plotly collect failed: {exc}")

try:
    datas += collect_data_files(
        "PySide6",
        includes=["**/Qt/resources/**", "**/Qt/translations/**", "**/Qt/bin/**"],
    )
except Exception as exc:  # noqa: BLE001
    print(f"WARNING: PySide6 data collect failed: {exc}")

hiddenimports += [
    "PySide6.QtWebEngineWidgets",
    "PySide6.QtWebEngineCore",
    "PySide6.QtWebChannel",
    "PySide6.QtPrintSupport",
    "PySide6.QtSvg",
]

hiddenimports += collect_submodules("repeatviewer")

ico = SPECDIR / "RepeatViewer.ico"
icon_path = str(ico) if ico.is_file() else None

a = Analysis(  # noqa: F821
    [str(ROOT / "RepeatViewer.py")],
    pathex=[str(ROOT)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[str(SPECDIR / "hooks")],
    hooksconfig={},
    runtime_hooks=[
        str(SPECDIR / "pyi_rth_six_importer.py"),
        str(SPECDIR / "pyi_rth_numpy_guard.py"),
        str(SPECDIR / "pyi_rth_scipy_distn.py"),
    ],
    excludes=[
        "tkinter",
        "PyQt5",
        "PyQt6",
        "IPython",
        "jupyter",
        "notebook",
    ],
    win_no_prefer_redirects=False,
    win_private_assemblies=False,
    cipher=block_cipher,
    noarchive=False,
)


def _dedupe_numpy_binaries(entries: list) -> list:
    """Keep one copy of each numpy C-extension (numpy 2.4+ forbids re-init)."""
    markers = (
        "_multiarray_umath",
        "_multiarray_tests",
        "_umath_linalg",
        "_pocketfft_umath",
    )
    groups: dict[str, list] = {}
    others: list = []
    for entry in entries:
        name = str(entry[0]).replace("\\", "/")
        base = name.rsplit("/", 1)[-1].lower()
        key = next((marker for marker in markers if marker in base), None)
        if key is None:
            others.append(entry)
            continue
        groups.setdefault(key, []).append(entry)

    out = list(others)
    for key, group in groups.items():
        if len(group) == 1:
            out.extend(group)
            continue
        nested = [
            item
            for item in group
            if "/numpy/" in str(item[0]).replace("\\", "/").lower()
        ]
        chosen = nested[0] if nested else group[0]
        dropped = [str(item[0]) for item in group if item is not chosen]
        print(
            f"WARNING: duplicate numpy binary {key!r}; keeping {chosen[0]!r}, "
            f"dropping {dropped}"
        )
        out.append(chosen)
    return out


a.binaries = _dedupe_numpy_binaries(list(a.binaries))

pyz = PYZ(a.pure, a.zipped_data, cipher=block_cipher)  # noqa: F821

exe = EXE(  # noqa: F821
    pyz,
    a.scripts,
    [],
    exclude_binaries=True,
    name="RepeatViewer",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=icon_path,
    version=None,
)

coll = COLLECT(  # noqa: F821
    exe,
    a.binaries,
    a.zipfiles,
    a.datas,
    strip=False,
    upx=False,
    upx_exclude=[],
    name="RepeatViewer",
)
