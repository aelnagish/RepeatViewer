# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller spec for RepeatViewer.app (macOS).

Build (from repo root, with venv active):

    python packaging/macos/packaging_utils.py   # optional .icns
    pyinstaller packaging/macos/RepeatViewer.spec
    bash packaging/macos/build_macos_dmg.sh
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
from bundle_datas import app_version, collect_annotation_python_deps  # noqa: E402

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
    "orfipy_core",
    "pyhmmer",
]

# Qt WebEngine + Flexidot need explicit collection; PySide6 hooks cover the rest.
for pkg in ("flexidot", "plotly"):
    try:
        pkg_datas, pkg_binaries, pkg_hidden = collect_all(pkg)
        datas += pkg_datas
        binaries += pkg_binaries
        hiddenimports += pkg_hidden
    except Exception as exc:  # noqa: BLE001
        print(f"WARNING: collect_all({pkg}) failed: {exc}")

ann_datas, ann_bins, ann_hidden = collect_annotation_python_deps()
datas += ann_datas
binaries += ann_bins
hiddenimports += ann_hidden

# WebEngine resources (locales, pak files, helpers).
try:
    datas += collect_data_files("PySide6", includes=["**/Qt/lib/**", "**/Qt/resources/**", "**/Qt/translations/**"])
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

icns = SPECDIR / "RepeatViewer.icns"
icon_path = str(icns) if icns.is_file() else None

a = Analysis(  # noqa: F821
    [str(ROOT / "RepeatViewer.py")],
    pathex=[str(ROOT)],
    binaries=binaries,
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
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
    argv_emulation=True,  # macOS: pass argv to frozen GUI apps
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=icon_path,
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

app = BUNDLE(  # noqa: F821
    coll,
    name="RepeatViewer.app",
    icon=icon_path,
    bundle_identifier="com.tetrimmer.repeatviewer",
    version=version,
    info_plist={
        "CFBundleName": "RepeatViewer",
        "CFBundleDisplayName": "RepeatViewer",
        "CFBundleShortVersionString": version,
        "CFBundleVersion": version,
        "NSHighResolutionCapable": True,
        "NSRequiresAquaSystemAppearance": False,
        "LSMinimumSystemVersion": "12.0",
        "NSPrincipalClass": "NSApplication",
        "CFBundleDocumentTypes": [],
    },
)
