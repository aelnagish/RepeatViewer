"""Toolbar icons: flat monochrome Tabler SVGs with Qt fallback."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtGui import QIcon
from PySide6.QtWidgets import QStyle, QWidget

from teviewer.widgets.svg_icons import themed_svg_icon

_DEFAULT_ICON_SIZE = 20
_ALIVIEW_ICON_SIZE = 22

# teviewer/assets/icons (preferred) and legacy backends/inco overrides
_ASSETS_DIR = Path(__file__).resolve().parent.parent / "assets" / "icons"
_LEGACY_DIR = Path(__file__).resolve().parent.parent / "backends" / "inco"

_NAME_FILES: dict[str, tuple[str, ...]] = {
    "aliview": ("aliview.svg", "aliview.png"),
    "align": ("align.svg",),
    "blast": ("blast.svg",),
    "te_genome_map": ("te_genome_map.svg", "blast.svg"),
    "extend": ("extend.svg",),
    "teaid": ("teaid.svg", "TEAid.png"),
    "msa_clean": ("msa_clean.svg",),
    "consensus": ("consensus.svg",),
    "flexidot": ("flexidot.svg", "flexidot.png"),
    "te_benchmark": ("te_benchmark.svg",),
    "te_divergence_plot": ("te_divergence.svg",),
    "minimap2": ("minimap2.svg",),
    "te_annotation": ("browser/sequence.svg",),
    "tsd_search": ("blast.svg",),
}

_FALLBACK_PIXMAP: dict[str, QStyle.StandardPixmap] = {
    "aliview": QStyle.StandardPixmap.SP_FileDialogDetailedView,
    "align": QStyle.StandardPixmap.SP_BrowserReload,
    "blast": QStyle.StandardPixmap.SP_ArrowForward,
    "te_genome_map": QStyle.StandardPixmap.SP_FileDialogListView,
    "extend": QStyle.StandardPixmap.SP_ArrowDown,
    "teaid": QStyle.StandardPixmap.SP_FileDialogInfoView,
    "msa_clean": QStyle.StandardPixmap.SP_DialogResetButton,
    "consensus": QStyle.StandardPixmap.SP_DialogApplyButton,
    "flexidot": QStyle.StandardPixmap.SP_DesktopIcon,
    "te_benchmark": QStyle.StandardPixmap.SP_FileDialogListView,
    "te_divergence_plot": QStyle.StandardPixmap.SP_FileDialogListView,
    "minimap2": QStyle.StandardPixmap.SP_BrowserReload,
    "te_annotation": QStyle.StandardPixmap.SP_FileIcon,
    "tsd_search": QStyle.StandardPixmap.SP_ArrowForward,
}


def _resolve_icon_asset(name: str) -> tuple[str, str] | None:
    """Return (svg_stem, subdir) for a toolbar function id."""
    for filename in _NAME_FILES.get(name, (f"{name}.svg",)):
        if "/" in filename:
            path = _ASSETS_DIR / filename
            if path.is_file():
                parts = Path(filename).parts
                subdir = parts[0] if len(parts) > 1 else ""
                return parts[-1].removesuffix(".svg"), subdir
        for folder in (_ASSETS_DIR, _LEGACY_DIR):
            path = folder / filename
            if path.is_file():
                subdir = "browser" if path.parent.name == "browser" else ""
                return path.stem, subdir
        browser_path = _ASSETS_DIR / "browser" / filename
        if browser_path.is_file():
            return Path(filename).stem, "browser"
    root_path = _ASSETS_DIR / f"{name}.svg"
    if root_path.is_file():
        return name, ""
    return None


def function_icon(
    widget: QWidget,
    name: str,
    *,
    size: int | None = None,
) -> QIcon:
    """Return a monochrome toolbar icon."""
    icon_size = size if size is not None else _DEFAULT_ICON_SIZE
    asset = _resolve_icon_asset(name)
    if asset is not None:
        stem, subdir = asset
        icon = themed_svg_icon(stem, icon_size, widget=widget, subdir=subdir)
        if not icon.isNull():
            return icon
    style = widget.style()
    pixmap = _FALLBACK_PIXMAP.get(name, QStyle.StandardPixmap.SP_FileIcon)
    return style.standardIcon(pixmap)


def aliview_icon(widget: QWidget) -> QIcon:
    """AliView toolbar icon."""
    return function_icon(widget, "aliview", size=_ALIVIEW_ICON_SIZE)
