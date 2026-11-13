"""Lightweight HTML output kind checks (no matplotlib / plotly / pandas).

Preview and policy code only needs filename heuristics. Keep those helpers here
so importing ``html_preview`` does not pull the full plotting stack at startup.
"""

from __future__ import annotations

from pathlib import Path

GENOME_MAP_HTML_SUFFIX = "_den_plot.html"
GENOME_MAP_REGION_HTML_SUFFIX = "_den_plot_region.html"
_LEGACY_GENOME_MAP_SUFFIXES = (
    "_density_plot.html",
    "_density_plot_region.html",
    "_genome_map.html",
    "_genome_map_region.html",
)


def is_flexidot_html_file(path: Path) -> bool:
    """True for a Flexidot interactive dotplot HTML output (``*_flexidot*.html``)."""
    name = path.name.lower()
    if not name.endswith((".html", ".htm")):
        return False
    return "_flexidot" in name


def is_genome_map_html_file(path: Path) -> bool:
    name = path.name.lower()
    if not name.endswith((".html", ".htm")):
        return False
    if "_dens_plot" in name or "_den_plot" in name or "_density_plot" in name:
        return True
    return (
        name.endswith(GENOME_MAP_HTML_SUFFIX)
        or name.endswith(GENOME_MAP_REGION_HTML_SUFFIX)
        or any(name.endswith(suffix) for suffix in _LEGACY_GENOME_MAP_SUFFIXES)
        or "_blast_genome_map" in name
        or "_genome_map" in name
    )
