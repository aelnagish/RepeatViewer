"""Plotly genome map of BLAST hits (per-chromosome bigwig-style tracks)."""

from __future__ import annotations

import json
import logging
import math
import os
import re
import shutil
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from Bio import SeqIO

from repeatviewer.job_options import (
    CHROMOSOME_SORT_NATURAL,
    CHROMOSOME_SORT_SIZE_ASC,
    CHROMOSOME_SORT_SIZE_DESC,
    DEFAULT_CHROMOSOME_SORT,
    TE_GENOME_MAP_COLOR_FAMILY,
    TE_GENOME_MAP_COLOR_IDENTITY,
    TE_GENOME_MAP_COLOR_ORDER,
    TE_GENOME_MAP_COLOR_SUPERFAMILY,
    TE_GENOME_MAP_DENSITY_HIT_COUNT,
    TE_GENOME_MAP_DENSITY_TE_LENGTH,
)
from repeatviewer.backends.bismark_cov import (
    METH_SUBTRACK_CG_COVERAGE_PCT,
    METH_SUBTRACK_CG_COVERED,
    METH_SUBTRACK_CG_TOTAL,
    METH_SUBTRACK_METH_LEVEL,
    MethylationChromBins,
    MethylationChromCov,
    bin_methylation_chromosome_worker,
    filter_methylation_cov_by_region,
    is_bismark_cov_file,
    load_bismark_cov,
    load_methylation_cov_by_chrom,
    map_bismark_chromosomes,
    methylation_cov_arrays,
    partition_cov_by_chrom_arrays,
    methylation_cov_site_count,
    methylation_diff_bins_worker,
    precompute_methylation_cg_totals,
    shared_cg_count_axis_max,
    _load_methylation_bins_chrom_cache,
    _methylation_bins_cache_path,
    _methylation_cov_cache_path,
    _save_methylation_bins_cache,
)
from repeatviewer.backends.blast_output_names import (
    blast_table_candidate_names,
    blast_table_source_tag,
    local_blast_txt_path,
    strip_blast_table_suffix as _strip_blast_table_suffix,
)
from repeatviewer.backends.genome_map_binning import bin_chromosome_tracks_worker
from repeatviewer.backends.job_cancel import JobCancellation, JobCancelledError, check_cancelled
from repeatviewer.plot_preferences import (
    DEFAULT_DIFF_COLOR_A,
    DEFAULT_DIFF_COLOR_B,
    DEFAULT_WINDOW_REF_COLOR,
    ensure_default_prefs_file,
    load_view_limits,
    load_view_state,
    merge_view_state_into_scale_meta,
    sidecar_prefs_path,
)

logger = logging.getLogger(__name__)


class GenomeMapNothingToPlotError(ValueError):
    """No plottable features and no selected chromosomes to render as empty tracks."""

ensure_default_prefs_file()
_VIEW_LIMITS = load_view_limits()

BLAST_COLUMNS = [
    "qseqid",
    "sseqid",
    "pident",
    "length",
    "mismatch",
    "gapopen",
    "qstart",
    "qend",
    "sstart",
    "send",
    "sstrand",
    "evalue",
    "bitscore",
]

DEFAULT_BIN_SIZE = 100_000  # bp per density window (100 kb)
MAX_CHROMOSOMES = 60
# Parallel binning kicks in above these thresholds (avoids process spawn overhead).
_PAR_BIN_MIN_HITS = 15_000
_PAR_BIN_MIN_CHROMS = 2
_PAR_METH_MIN_HITS = 0
_PAR_METH_MIN_CHROMS = 1
_PLOTLY_JS_FILENAME = "plotly.min.js"
# Vertical layout: grey track per chromosome; density bars stay below the next row
ROW_Y_STEP = 2.45
TRACK_HALF = 0.10
BAR_GAP_BELOW_TRACK = 0.04
BAR_GAP_ABOVE_NEXT = 0.22
# Diff track: tighter step from last TE input, extra gap before next chromosome block.
DIFF_TRACK_BEFORE_STEP = 1.72
DIFF_TRACK_AFTER_GAP = 1.20
DIFF_TRACK_BAR_HALF_FRAC = 0.36
DIFF_AXIS_TICK_COUNT = 5
DIFF_AXIS_LABEL_PAD = 0.22
BP_PER_MBP = 1_000_000.0
TE_INPUT_KIND_MASKING = "masking"
TE_INPUT_KIND_CONSENSUS = "consensus"


def _track_y_positions(
    tracks: list["_GenomeMapTrack"],
    *,
    show_te_diff: bool = False,
    show_meth_diff: bool = False,
) -> dict[str, float]:
    """First listed track at top; optional compact diff row spacing."""
    if not tracks:
        return {}
    use_diff_spacing = show_te_diff or show_meth_diff
    if not use_diff_spacing:
        n = len(tracks)
        return {t.label: (n - 1 - idx) * ROW_Y_STEP for idx, t in enumerate(tracks)}

    y_bottom_up: dict[str, float] = {}
    y = 0.0
    for idx, track in enumerate(tracks):
        if idx > 0:
            prev = tracks[idx - 1]
            if track.kind in {"te_diff", "meth_diff"}:
                y += DIFF_TRACK_BEFORE_STEP
            elif prev.kind in {"te_diff", "meth_diff"}:
                y += ROW_Y_STEP + DIFF_TRACK_AFTER_GAP
            else:
                y += ROW_Y_STEP
        y_bottom_up[track.label] = y
    max_y = max(y_bottom_up.values())
    return {label: max_y - y_val for label, y_val in y_bottom_up.items()}


HOVER_TRACK_POINTS_PER_CHROM = 36

# Cached genome FASTA metadata: path -> (mtime, names, sizes)
_GENOME_META_CACHE: dict[str, tuple[float, list[str], dict[str, int]]] = {}

PLOTLY_HTML_CONFIG: dict[str, Any] = {
    "scrollZoom": True,
    "displayModeBar": True,
    "displaylogo": False,
    # responsive reflow can undo restyle/redraw bar-height updates in Qt WebEngine.
    "responsive": False,
}


def _ensure_plotly_js_asset(output_dir: Path) -> None:
    """Copy plotly.min.js beside the HTML once (avoids embedding ~3 MB per export)."""
    dest = output_dir / _PLOTLY_JS_FILENAME
    if dest.is_file():
        return
    src = Path(pio.__file__).resolve().parent / "package_data" / _PLOTLY_JS_FILENAME
    if not src.is_file():
        return
    try:
        shutil.copy2(src, dest)
    except OSError as exc:
        logger.debug("Could not copy Plotly JS to %s: %s", dest, exc)


def _figure_to_preview_html(fig: go.Figure, output_html: Path) -> str:
    """Serialize a Plotly figure to self-contained HTML (works in external browsers via file://)."""
    output_html = output_html.expanduser().resolve()
    _ensure_plotly_js_asset(output_html.parent)
    return fig.to_html(
        include_plotlyjs=True,
        full_html=True,
        config=PLOTLY_HTML_CONFIG,
    )

# Interactive HTML panel: density / genome-track bar height multiplier (from prefs file).
PLOT_SCALE_MIN = _VIEW_LIMITS.plot_scale_min
PLOT_SCALE_MAX = _VIEW_LIMITS.plot_scale_max
PLOT_SCALE_MAX_DIVERGENCE = _VIEW_LIMITS.plot_scale_max_divergence
PLOT_SCALE_MAX_GENOME = _VIEW_LIMITS.plot_scale_max_genome
# Divergence histogram: bin width on x-axis (1 - identity), 0–1 scale.
DEFAULT_DIVERGENCE_BIN_WIDTH = 0.005
PLOT_SCALE_STEP = _VIEW_LIMITS.plot_scale_step
PLOT_XZOOM_MIN = 1.0
PLOT_XZOOM_MAX_GENOME = 50.0
PLOT_XZOOM_MAX_DIVERGENCE = 20.0
PLOT_XZOOM_STEP = 0.5
ROW_SPACING_MIN = _VIEW_LIMITS.row_spacing_min
ROW_SPACING_MAX = _VIEW_LIMITS.row_spacing_max
ROW_SPACING_STEP = _VIEW_LIMITS.row_spacing_step
DIFF_TRACK_SCALE_MIN = _VIEW_LIMITS.diff_track_scale_min
DIFF_TRACK_SCALE_MAX = _VIEW_LIMITS.diff_track_scale_max
DIFF_TRACK_SCALE_STEP = _VIEW_LIMITS.diff_track_scale_step
# Left margin for chromosome-name column; density scale sits in the x-axis domain gutter.
GENOME_MAP_MARGIN_LEFT = 148
CHROM_LABEL_PAPER_X = 0.0
CHROM_LABEL_XSHIFT = -8
# Narrow left strip reserved for density y-axes (empty when zooming x).
XAXIS_DOMAIN_START = 0.042
DENSITY_TICK_W_PAPER = 0.005

# TE landscape plot defaults (TE_landscape_plot.R / TEtrimmer taxonomy).
TE_TYPE_CATEGORIES: dict[str, tuple[str, ...]] = {
    "LTR": (
        "LTR",
        "RLC",
        "RLB",
        "RLR",
        "RLE",
        "RLG",
        "RLX",
        "Gypsy",
        "Copia",
        "BEL",
    ),
    "LINE": ("LINE", "RIL", "RIT", "RIX", "Penelope"),
    "SINE": ("SINE", "RST", "RSX"),
    "DNA": (
        "DNA",
        "TIR",
        "hAT",
        "hAT-Ac",
        "MULE",
        "MLE",
        "MuDR",
        "Tourist",
        "CACTA",
        "PILE",
        "POLE",
        "Stowaway",
        "TcMar-Stowaway",
        "PIF",
        "Harbinger",
        "Tc1",
        "En-Spm",
        "EnSpm",
        "CMC-EnSpm",
        "PiggyBac",
        "Mirage",
        "P-element",
        "Transib",
        "DTA",
        "DTH",
        "DTT",
        "DTM",
        "DTC",
        "DTX",
        "DTR",
        "DTE",
        "Merlin",
        "DTP",
        "DTB",
        "polinton",
        "Academ",
        "Crypton",
    ),
    "MITE": ("MITE",),
    "Helitron": ("Helitron", "DHH", "DHX", "helitron"),
    "Unknown": ("Unknown", "unknown", "Unspecified", "Unclassified"),
    "simple_repeat": (
        "Low_complexity",
        "repeat_region",
        "Simple_repeat",
        "Low",
        "Satellite",
    ),
    "Other_non-LTR": ("Non-LTR", "Retroposon", "nonLTR", "YR"),
}
TE_CATEGORY_COLORS: dict[str, str] = {
    "LTR": "#97d672",
    "LINE": "#d62728",
    "SINE": "#ffbb78",
    "Other_non-LTR": "#17becf",
    "DNA": "#1f77b4",
    "MITE": "#9467bd",
    "Helitron": "#e377c2",
    "Unknown": "#919697",
    "simple_repeat": "#8c564b",
    "Non-TE": "#919697",
}
TE_CATEGORY_ORDER: tuple[str, ...] = (
    "LTR",
    "LINE",
    "SINE",
    "Other_non-LTR",
    "DNA",
    "MITE",
    "Helitron",
    "Unknown",
    "simple_repeat",
    "Non-TE",
)
_TE_TYPE_TO_CATEGORY: dict[str, str] = {}
for _parent, _types in TE_TYPE_CATEGORIES.items():
    for _te_type in _types:
        _TE_TYPE_TO_CATEGORY[_te_type] = _parent
        _TE_TYPE_TO_CATEGORY[_te_type.lower()] = _parent
_TE_CATEGORY_RANK: dict[str, int] = {
    name: idx for idx, name in enumerate(TE_CATEGORY_ORDER)
}
# Distinct non-grey colours for TE types outside the default taxonomy (e.g. RC).
TE_EXTRA_TYPE_COLORS: tuple[str, ...] = (
    "#bcbd22",
    "#aec7e8",
    "#ff9896",
    "#c5b0d5",
    "#dbdb8d",
    "#9edae5",
    "#f7b6d2",
    "#637939",
    "#393b79",
    "#ad494a",
    "#5254a3",
    "#e7969c",
    "#843c39",
    "#6b6ecf",
    "#cedb9c",
    "#e7ba52",
)


def _density_paper_x_in_gutter() -> float:
    """Paper x for density axis line, centred in the domain gutter."""
    return XAXIS_DOMAIN_START * 0.58


def _density_label_paper_x(axis_paper_x: float) -> float:
    return axis_paper_x - DENSITY_TICK_W_PAPER * 1.15


def _track_pad_mbp(x_span_mbp: float) -> float:
    """Small genomic x padding so grey bars start just right of the density gutter."""
    span = max(0.5, float(x_span_mbp))
    return max(0.32, min(0.9, span * 0.010))


def _genomic_xaxis_ticks(
    chrom_start: float,
    view_hi: float,
    pad: float,
) -> tuple[list[float], list[str]]:
    span = max(0.5, float(view_hi) - float(chrom_start))
    if span > 40:
        step = 10.0
    elif span > 15:
        step = 5.0
    else:
        step = max(0.5, round(span / 6.0, 1))
    start = math.floor(float(chrom_start) / step) * step
    genomic: list[float] = []
    g = start
    while g <= float(view_hi) + step * 0.001:
        genomic.append(float(g))
        g += step
    tickvals = [t + pad for t in genomic]
    ticktext = [
        str(int(t)) if abs(t - round(t)) < 0.05 else f"{t:.1f}"
        for t in genomic
    ]
    return tickvals, ticktext


def _plot_scale_slider_attrs(*, max_scale: float | None = None) -> str:
    cap = PLOT_SCALE_MAX if max_scale is None else max_scale
    return f'min="{PLOT_SCALE_MIN}" max="{cap}" step="{PLOT_SCALE_STEP}"'


def _diff_axis_tick_fractions(count: int = DIFF_AXIS_TICK_COUNT) -> list[float]:
    n = max(3, int(count))
    return [i / (n - 1) for i in range(n)]


def _diff_axis_tick_label(frac: float, track_max: float, *, te_length_mode: bool) -> str:
    if abs(frac - 0.5) < 1e-6:
        return "0"
    sign = 1.0 if frac > 0.5 else -1.0
    magnitude = track_max * abs(frac - 0.5) * 2.0
    text = _format_density_tick_value(magnitude, te_length_mode=te_length_mode)
    return f"-{text}" if sign < 0 else text


def _meth_diff_axis_tick_label(frac: float, track_max: float) -> str:
    if abs(frac - 0.5) < 1e-6:
        return "0"
    sign = 1.0 if frac > 0.5 else -1.0
    magnitude = track_max * abs(frac - 0.5) * 2.0
    text = _format_unit_interval_axis_tick_value(magnitude)
    return f"-{text}" if sign < 0 else text


def _format_fraction_axis_tick_value(value: float) -> str:
    """Format a 0–1 axis value as a percentage tick label."""
    if value <= 0:
        return "0%"
    pct = value * 100.0
    if pct >= 10:
        return f"{pct:.0f}%"
    if pct >= 1:
        return f"{pct:.0f}%" if abs(pct - round(pct)) < 0.05 else f"{pct:.1f}%"
    return f"{pct:.1f}%"


def _format_unit_interval_axis_tick_value(value: float) -> str:
    """Format a 0–1 axis value as a decimal tick label."""
    if value <= 0:
        return "0"
    if value < 1.0:
        return f"{value:.2f}"
    return f"{value:.1f}"


def _format_density_tick_value(value: float, *, te_length_mode: bool) -> str:
    if value <= 0:
        return "0"
    if te_length_mode:
        if value < 1000:
            return f"{int(round(value))} bp"
        kb = value / 1000.0
        if kb >= 100:
            return f"{kb:.0f} kb"
        if kb >= 10:
            return (
                f"{kb:.0f} kb"
                if abs(kb - round(kb)) < 0.05
                else f"{kb:.1f} kb"
            )
        if kb >= 1:
            return (
                f"{kb:.0f} kb"
                if abs(kb - round(kb)) < 0.05
                else f"{kb:.1f} kb"
            )
        return f"{kb:.2f} kb"
    return str(int(round(value)))


def _density_axis_tick_fractions(tick_count: int = 5) -> list[float]:
    n = max(2, int(tick_count))
    return [i / (n - 1) for i in range(n)]


def _density_axis_tick_count(axis_height: float, *, default: int = 5) -> int:
    """Fewer tick labels when the mini axis is short, to avoid overlap."""
    if axis_height < 0.55:
        return 3
    if axis_height < 1.0:
        return 4
    return default


def _density_axis_height_for_track(
    _track_label: str,
    _track_labels: list[str],
    _y_pos: dict[str, float],
    *,
    max_bar_height: float,
    top_first: bool = True,
) -> float:
    """Nominal density-axis height; bars may extend into the track above when scaled."""
    return max_bar_height


def _format_track_label_for_display(label: str, *, is_diff: bool = False) -> str:
    """Break long track labels before a parenthetical suffix for the left margin."""
    meth_diff_marker = " (methylation diff:"
    if meth_diff_marker in label:
        chrom = label.split(meth_diff_marker, 1)[0]
        rest = label.split(meth_diff_marker, 1)[1].rstrip(")")
        parts = [p.strip() for p in rest.split(" vs ", 1)]
        if len(parts) == 2:
            short_a = _compact_bismark_cov_label(parts[0], max_len=18)
            short_b = _compact_bismark_cov_label(parts[1], max_len=18)
            return f"{chrom}<br>meth diff<br>{short_a} vs {short_b}"
        return f"{chrom}<br>meth diff"
    if is_diff:
        diff_idx = label.find(" (diff:")
        if diff_idx > 0:
            chrom = label[:diff_idx]
            rest = label[diff_idx + len(" (diff:") :].rstrip(")")
            parts = [p.strip() for p in rest.split(" vs ", 1)]
            if len(parts) == 2:
                short_a = _compact_bismark_cov_label(parts[0], max_len=18)
                short_b = _compact_bismark_cov_label(parts[1], max_len=18)
                return f"{chrom}<br>TE diff<br>{short_a} vs {short_b}"
            return chrom
        paren = label.find(" (")
        if paren > 0:
            return label[:paren]
    idx = label.find(" (")
    if idx > 0:
        chrom = label[:idx]
        rest = label[idx + 2 :].rstrip(")")
        if ", " in rest:
            suffix, source = rest.rsplit(", ", 1)
            short_src = _compact_bismark_cov_label(source, max_len=20)
            return f"{chrom}<br>{suffix}<br>{short_src}"
        return f"{chrom}<br>{rest}"
    return label


def _genome_map_margin_left(track_labels: list[str]) -> int:
    """Left margin wide enough for density axes and the longest track label line."""
    max_line_chars = 18
    max_label_line_chars = 34
    for label in track_labels:
        display = _format_track_label_for_display(label)
        for line in display.split("<br>"):
            max_line_chars = max(
                max_line_chars,
                min(len(line), max_label_line_chars),
            )
    label_px = int(max_line_chars * 5.4) + 16
    gutter_px = 58
    return max(GENOME_MAP_MARGIN_LEFT, label_px + gutter_px)


def _genome_map_top_margin(title: str, caption: str, *, width: int = 1120) -> int:
    """Top margin tall enough for multi-line title + caption without overlapping tracks."""
    chars_per_line = max(48, int(width / 10.5))
    title_lines = max(1, (len(title) + chars_per_line - 1) // chars_per_line)
    caption_lines = max(1, (len(caption) + chars_per_line - 1) // chars_per_line)
    return max(136, 34 + title_lines * 24 + caption_lines * 15)


def _append_chromosome_labels(
    fig: go.Figure,
    *,
    chrom_order: list[str],
    y_pos: dict[str, float],
    track_kinds: dict[str, str] | None = None,
    track_chroms: dict[str, str] | None = None,
    chrom_sizes: dict[str, int] | None = None,
    diff_bar_half: float | None = None,
) -> list[int]:
    """Chromosome IDs in the left margin (paper x), separate from density tick labels."""
    ann_indices: list[int] = []
    seen_chrom_sizes: set[str] = set()
    for label in chrom_order:
        y_center = float(y_pos[label])
        is_diff = (track_kinds or {}).get(label) in {"te_diff", "meth_diff"}
        half_h = float(diff_bar_half) if is_diff and diff_bar_half is not None else 0.0
        label_y = y_center - half_h - 0.10 if is_diff else y_center
        yanchor = "top" if is_diff else "middle"
        chrom = (track_chroms or {}).get(label) or _chromosome_from_track_label(label)
        display = _format_track_label_for_display(label, is_diff=is_diff)
        if (
            chrom_sizes
            and chrom in chrom_sizes
            and chrom not in seen_chrom_sizes
            and not is_diff
        ):
            seen_chrom_sizes.add(chrom)
            size_line = _format_genome_length_bp(chrom_sizes[chrom])
            head, tail = (display.split("<br>", 1) + [""])[:2]
            display = f"{head}<br>{size_line}" + (f"<br>{tail}" if tail else "")
        ann_idx = len(fig.layout.annotations) if fig.layout.annotations else 0
        fig.add_annotation(
            x=CHROM_LABEL_PAPER_X,
            y=label_y,
            xref="paper",
            yref="y",
            text=display,
            showarrow=False,
            xanchor="right",
            yanchor=yanchor,
            xshift=CHROM_LABEL_XSHIFT - (6 if is_diff else 0),
            font={"size": 8 if is_diff else 9, "color": "#222222"},
            align="right",
        )
        ann_indices.append(ann_idx)
    return ann_indices


def _append_track_density_axes(
    fig: go.Figure,
    *,
    chrom_order: list[str],
    y_pos: dict[str, float],
    max_bar_height: float,
    density_max: float,
    window_kb: int,
    te_length_mode: bool,
    tick_count: int = 5,
    density_paper_x: float | None = None,
    gene_density_max: float | None = None,
    annotation_length_mode: bool = True,
    track_kinds: dict[str, str] | None = None,
    track_density_norm: dict[str, float] | None = None,
    track_input_kinds: dict[str, str] | None = None,
    track_methylation_subtracks: dict[str, str] | None = None,
    diff_bar_half: float | None = None,
) -> dict[str, Any]:
    """Mini y-axes in the left plot gutter (paper coords; empty strip when zooming x)."""
    axis_paper_x = float(
        density_paper_x if density_paper_x is not None else _density_paper_x_in_gutter()
    )
    tick_w_paper = DENSITY_TICK_W_PAPER
    label_paper_x = _density_label_paper_x(axis_paper_x)
    axis_line = {"color": "#555555", "width": 1.5}
    tick_line = {"color": "#555555", "width": 1.0}
    tick_fracs = _density_axis_tick_fractions(tick_count)
    tracks: list[dict[str, Any]] = []

    for chrom in chrom_order:
        is_gene = (track_kinds or {}).get(chrom) == "gene"
        is_methylation = (track_kinds or {}).get(chrom) == "methylation"
        is_te_diff = (track_kinds or {}).get(chrom) == "te_diff"
        is_meth_diff = (track_kinds or {}).get(chrom) == "meth_diff"
        is_diff = is_te_diff or is_meth_diff
        if is_diff:
            track_max = float((track_density_norm or {}).get(chrom, density_max))
            y_center = float(y_pos[chrom])
            half_h = float(
                diff_bar_half if diff_bar_half is not None else max_bar_height * 0.45
            )
            y0 = y_center - half_h
            y1 = y_center + half_h
            shape_idx = len(fig.layout.shapes) if fig.layout.shapes else 0
            fig.add_shape(
                type="line",
                xref="paper",
                yref="y",
                x0=axis_paper_x,
                x1=axis_paper_x,
                y0=y0,
                y1=y1,
                line=axis_line,
            )
            tick_meta: list[dict[str, Any]] = []
            for frac in _diff_axis_tick_fractions():
                y_tick = y0 + (y1 - y0) * frac
                tick_shape_idx = len(fig.layout.shapes) if fig.layout.shapes else 0
                ann_idx = len(fig.layout.annotations) if fig.layout.annotations else 0
                fig.add_shape(
                    type="line",
                    xref="paper",
                    yref="y",
                    x0=axis_paper_x - tick_w_paper,
                    x1=axis_paper_x,
                    y0=y_tick,
                    y1=y_tick,
                    line=tick_line,
                )
                label = (
                    _meth_diff_axis_tick_label(frac, track_max)
                    if is_meth_diff
                    else _diff_axis_tick_label(frac, track_max, te_length_mode=True)
                )
                fig.add_annotation(
                    x=label_paper_x,
                    y=y_tick,
                    xref="paper",
                    yref="y",
                    text=label,
                    showarrow=False,
                    xanchor="right",
                    font={"size": 6, "color": "#333333"},
                )
                tick_meta.append(
                    {
                        "frac": float(frac),
                        "shapeIdx": tick_shape_idx,
                        "annIdx": ann_idx,
                    }
                )
            tracks.append(
                {
                    "chrom": chrom,
                    "shapeIdx": shape_idx,
                    "y0": y0,
                    "yCenter": y_center,
                    "ticks": tick_meta,
                    "densityNormMax": track_max,
                    "diffNormMax": track_max,
                    "isDiffTrack": True,
                    "isMethDiffTrack": is_meth_diff,
                    "isGeneTrack": False,
                }
            )
            continue
        meth_subtrack = (track_methylation_subtracks or {}).get(chrom, "")
        if is_methylation and meth_subtrack == METH_SUBTRACK_CG_COVERAGE_PCT:
            per_track = (track_density_norm or {}).get(chrom)
            track_max = float(
                min(
                    METH_COVERAGE_PCT_AXIS_MAX,
                    per_track if per_track is not None else METH_COVERAGE_PCT_AXIS_MAX,
                )
            )
        elif is_methylation and meth_subtrack == METH_SUBTRACK_METH_LEVEL:
            track_max = float(METH_LEVEL_AXIS_MAX)
        elif is_gene and gene_density_max is not None and gene_density_max > 0:
            track_max = float(gene_density_max)
        else:
            per_track = (track_density_norm or {}).get(chrom)
            track_max = float(per_track if per_track is not None else density_max)
        y0 = float(y_pos[chrom]) + TRACK_HALF + BAR_GAP_BELOW_TRACK
        axis_height = _density_axis_height_for_track(
            chrom,
            chrom_order,
            y_pos,
            max_bar_height=max_bar_height,
            top_first=True,
        )
        tick_fracs_row = _density_axis_tick_fractions(
            _density_axis_tick_count(axis_height, default=tick_count),
        )
        y1 = y0 + axis_height
        shape_idx = len(fig.layout.shapes) if fig.layout.shapes else 0

        fig.add_shape(
            type="line",
            xref="paper",
            yref="y",
            x0=axis_paper_x,
            x1=axis_paper_x,
            y0=y0,
            y1=y1,
            line=axis_line,
        )

        tick_meta: list[dict[str, Any]] = []
        for frac in tick_fracs_row:
            y_tick = y0 + axis_height * frac
            tick_shape_idx = len(fig.layout.shapes) if fig.layout.shapes else 0
            ann_idx = len(fig.layout.annotations) if fig.layout.annotations else 0
            fig.add_shape(
                type="line",
                xref="paper",
                yref="y",
                x0=axis_paper_x - tick_w_paper,
                x1=axis_paper_x,
                y0=y_tick,
                y1=y_tick,
                line=tick_line,
            )
            if is_methylation and meth_subtrack == METH_SUBTRACK_CG_COVERAGE_PCT:
                label = _format_unit_interval_axis_tick_value(track_max * frac)
            elif is_methylation and meth_subtrack == METH_SUBTRACK_METH_LEVEL:
                label = _format_unit_interval_axis_tick_value(track_max * frac)
            elif is_methylation:
                label = str(int(round(track_max * frac)))
            elif is_gene:
                label = _format_density_tick_value(
                    track_max * frac,
                    te_length_mode=annotation_length_mode,
                )
            else:
                label = _format_density_tick_value(
                    track_max * frac,
                    te_length_mode=te_length_mode,
                )
            fig.add_annotation(
                x=label_paper_x,
                y=y_tick,
                xref="paper",
                yref="y",
                text=label,
                showarrow=False,
                xanchor="right",
                font={"size": 7, "color": "#333333"},
            )
            tick_meta.append(
                {
                    "frac": float(frac),
                    "shapeIdx": tick_shape_idx,
                    "annIdx": ann_idx,
                }
            )

        tracks.append(
            {
                "chrom": chrom,
                "shapeIdx": shape_idx,
                "y0": y0,
                "ticks": tick_meta,
                "axisHeight": float(axis_height),
                "densityNormMax": track_max,
                "isGeneTrack": is_gene,
                "isMethylationTrack": is_methylation,
                "methylationSubtrack": meth_subtrack or None,
                "isDiffTrack": is_diff,
                "teInputKind": (
                    "gene"
                    if is_gene
                    else (
                        "methylation"
                        if is_methylation
                        else (track_input_kinds or {}).get(chrom, TE_INPUT_KIND_MASKING)
                    )
                ),
            }
        )

    return {
        "tracks": tracks,
        "maxBarHeight": float(max_bar_height),
        "densityMax": float(density_max),
        "teLengthMode": te_length_mode,
        "annotationLengthMode": annotation_length_mode,
        "windowKb": int(window_kb),
        "tickCount": len(tick_fracs),
        "densityPaperAxis": True,
    }


def _window_ref_density_value(*, window_kb: int, te_length_mode: bool) -> float:
    """Density value for a full window on the mini y-axis (e.g. 100 kb TE per 100 kb window)."""
    if te_length_mode:
        return float(window_kb * 1000)
    return 1.0


def _append_window_ref_lines(
    fig: go.Figure,
    *,
    track_labels: list[str],
    track_chroms: dict[str, str],
    y_pos: dict[str, float],
    chrom_sizes: dict[str, int],
    track_x_off: float,
    max_bar_height: float,
    density_norm_max: float,
    window_kb: int,
    te_length_mode: bool,
    density_paper_x: float,
    line_color: str = DEFAULT_WINDOW_REF_COLOR,
    track_kinds: dict[str, str] | None = None,
    track_density_norm: dict[str, float] | None = None,
    show_lines: bool = True,
) -> dict[str, dict[str, int]]:
    """Dashed horizontal reference line per track at the density window size."""
    ref_value = _window_ref_density_value(window_kb=window_kb, te_length_mode=te_length_mode)
    if ref_value <= 0:
        return {}
    tick_w_paper = DENSITY_TICK_W_PAPER
    ref_color = str(line_color or DEFAULT_WINDOW_REF_COLOR)
    ref_line = {"color": ref_color, "width": 1.2, "dash": "dash"}
    ref_tick = {"color": ref_color, "width": 1.2}
    meta: dict[str, dict[str, int]] = {}

    for label in track_labels:
        chrom = track_chroms.get(label, label)
        if chrom not in chrom_sizes:
            continue
        if (track_kinds or {}).get(label) in ("te_diff", "meth_diff", "gene", "methylation"):
            continue
        per_track = (track_density_norm or {}).get(label)
        track_norm = float(per_track if per_track is not None else density_norm_max)
        if track_norm <= 0:
            continue
        density_scale = max_bar_height / track_norm
        ref_height = ref_value * density_scale
        y_base = float(y_pos[label])
        ground = y_base + TRACK_HALF + BAR_GAP_BELOW_TRACK
        y_line = ground + ref_height
        within_scale = ref_value <= track_norm * 1.001
        line_visible = bool(show_lines and within_scale)
        length_mbp = chrom_sizes[chrom] / BP_PER_MBP
        line_idx = len(fig.layout.shapes) if fig.layout.shapes else 0
        fig.add_shape(
            type="line",
            x0=track_x_off,
            x1=track_x_off + length_mbp,
            y0=y_line,
            y1=y_line,
            line=ref_line,
            layer="above",
            visible=line_visible,
        )
        tick_idx = len(fig.layout.shapes) if fig.layout.shapes else 0
        fig.add_shape(
            type="line",
            xref="paper",
            yref="y",
            x0=density_paper_x - tick_w_paper,
            x1=density_paper_x,
            y0=y_line,
            y1=y_line,
            line=ref_tick,
            layer="above",
            visible=line_visible,
        )
        meta[label] = {"lineShapeIdx": line_idx, "tickShapeIdx": tick_idx}
    return meta


def _plot_xzoom_slider_attrs(*, max_scale: float) -> str:
    return (
        f'min="{PLOT_XZOOM_MIN}" max="{max_scale}" step="{PLOT_XZOOM_STEP}"'
    )


def _plotly_interactive_controls_html(
    *,
    scale_label: str,
    panel_mode: str,
    prefs_path: str | None = None,
) -> str:
    """Inject draggable filter panel + bar-height scale slider into Plotly HTML."""
    genome_mode = panel_mode == "genome"
    prefs_path_js = json.dumps(prefs_path)
    scale_attrs = _plot_scale_slider_attrs(
        max_scale=PLOT_SCALE_MAX_GENOME if genome_mode else PLOT_SCALE_MAX_DIVERGENCE,
    )
    scale_max_js = PLOT_SCALE_MAX_GENOME if genome_mode else PLOT_SCALE_MAX_DIVERGENCE
    xzoom_max_js = PLOT_XZOOM_MAX_GENOME if genome_mode else PLOT_XZOOM_MAX_DIVERGENCE
    xzoom_attrs = _plot_xzoom_slider_attrs(max_scale=xzoom_max_js)
    filter_title = (
        "TE types — legend (show / colour / order)"
        if genome_mode
        else "TE classes — legend"
    )
    # Genome map: stacked category bars + plain identity-colour bars (no legendgroup).
    apply_body = (
        _GENOME_MAP_SCALE_APPLY_JS.replace("__ROW_SPACING_MIN__", str(ROW_SPACING_MIN))
        .replace("__ROW_SPACING_MAX__", str(ROW_SPACING_MAX))
        .replace("__DIFF_TRACK_SCALE_MIN__", str(DIFF_TRACK_SCALE_MIN))
        .replace("__DIFF_TRACK_SCALE_MAX__", str(DIFF_TRACK_SCALE_MAX))
        .replace("{{", "{")
        .replace("}}", "}")
        if genome_mode
        else _DIVERGENCE_SCALE_APPLY_JS
    )
    xzoom_panel_html = (
        ""
        if genome_mode
        else (
            '<div class="te-sf-scale">X axis zoom: <span id="te-sf-xzoom-val">1.0x</span>'
            f'<input id="te-sf-xzoom" type="range" {xzoom_attrs} value="1"></div>'
        )
    )
    genome_extra_html = (
        (
            '<div class="te-sf-scale">Track spacing: <span id="te-sf-rowsp-val">1.0x</span>'
            f'<input id="te-sf-rowsp" type="range" min="{ROW_SPACING_MIN}" max="{ROW_SPACING_MAX}"'
            f' step="{ROW_SPACING_STEP}" value="1"></div>'
            '<div id="te-sf-diffscale-wrap" class="te-sf-scale te-sf-hidden">'
            'TE difference track scale: <span id="te-sf-diffscale-val">1.0x</span>'
            f'<input id="te-sf-diffscale" type="range" min="{DIFF_TRACK_SCALE_MIN}"'
            f' max="{DIFF_TRACK_SCALE_MAX}" step="{DIFF_TRACK_SCALE_STEP}" value="1"></div>'
            '<div id="te-sf-diffcolors-wrap" class="te-sf-hidden">'
            '<div class="te-sf-item te-sf-diff-color">'
            f'<input type="color" id="te-sf-diffcolor-a" class="te-sf-color" value="{TE_DIFF_COLOR_A}">'
            '<span id="te-sf-diffcolor-a-lbl">TE difference track positive color</span></div>'
            '<div class="te-sf-item te-sf-diff-color">'
            f'<input type="color" id="te-sf-diffcolor-b" class="te-sf-color" value="{TE_DIFF_COLOR_B}">'
            '<span id="te-sf-diffcolor-b-lbl">TE difference track negative color</span></div></div>'
            '<div class="te-sf-item te-sf-window-ref">'
            '<input type="checkbox" id="te-sf-window-ref" checked>'
            f'<input type="color" id="te-sf-window-ref-color" class="te-sf-color" value="{DEFAULT_WINDOW_REF_COLOR}">'
            "<span>Show window-size reference line</span></div>"
            '<div id="te-sf-methdiffscale-wrap" class="te-sf-scale te-sf-hidden">'
            'Methylation difference track scale: <span id="te-sf-methdiffscale-val">1.0x</span>'
            f'<input id="te-sf-methdiffscale" type="range" min="{DIFF_TRACK_SCALE_MIN}"'
            f' max="{DIFF_TRACK_SCALE_MAX}" step="{DIFF_TRACK_SCALE_STEP}" value="1"></div>'
            '<div id="te-sf-methdiffcolors-wrap" class="te-sf-hidden">'
            '<div class="te-sf-item te-sf-diff-color">'
            f'<input type="color" id="te-sf-methdiffcolor-a" class="te-sf-color" value="{METH_DIFF_COLOR_A}">'
            '<span id="te-sf-methdiffcolor-a-lbl">Methylation difference track positive color</span></div>'
            '<div class="te-sf-item te-sf-diff-color">'
            f'<input type="color" id="te-sf-methdiffcolor-b" class="te-sf-color" value="{METH_DIFF_COLOR_B}">'
            '<span id="te-sf-methdiffcolor-b-lbl">Methylation difference track negative color</span></div></div>'
            '<div id="te-sf-methcolors-wrap" class="te-sf-hidden">'
            '<div class="te-sf-item te-sf-meth-color">'
            '<input type="checkbox" id="te-sf-methshow-cgtotal">'
            f'<input type="color" id="te-sf-methcolor-cgtotal" class="te-sf-color" value="{METH_COLOR_CG_TOTAL}">'
            "<span>Total C+G number</span></div>"
            '<div class="te-sf-item te-sf-meth-color">'
            '<input type="checkbox" id="te-sf-methshow-cgcovered">'
            f'<input type="color" id="te-sf-methcolor-cgcovered" class="te-sf-color" value="{METH_COLOR_CG_COVERED}">'
            "<span>Covered C+G number</span></div>"
            '<div class="te-sf-item te-sf-meth-color">'
            '<input type="checkbox" id="te-sf-methshow-cgpct" checked>'
            f'<input type="color" id="te-sf-methcolor-cgpct" class="te-sf-color" value="{METH_COLOR_CG_COVERAGE_PCT}">'
            "<span>Covered fraction</span></div>"
            '<div class="te-sf-item te-sf-meth-color">'
            '<input type="checkbox" id="te-sf-methshow-level" checked>'
            f'<input type="color" id="te-sf-methcolor-level" class="te-sf-color" value="{METH_COLOR_METH_LEVEL}">'
            "<span>Methylation level</span></div>"
            '<div class="te-sf-item te-sf-meth-color te-sf-hidden" id="te-sf-methshow-diff-row">'
            '<input type="checkbox" id="te-sf-methshow-diff" checked>'
            "<span>Methylation difference track</span></div></div>"
        )
        if genome_mode
        else ""
    )
    stack_hint_html = ""
    return f"""
<style>
.te-sf-panel{{position:fixed;right:14px;top:14px;z-index:100000;background:#fff;border:1px solid #d0d7de;border-radius:8px;padding:8px;max-height:70vh;overflow:auto;min-width:230px;font-family:Arial,sans-serif;font-size:12px;box-shadow:0 4px 16px rgba(0,0,0,.15)}}
.te-sf-title{{font-weight:600;margin-bottom:6px;cursor:move;user-select:none}}
.te-sf-actions{{display:flex;gap:6px;margin-bottom:6px}}
.te-sf-actions button{{font-size:11px;padding:3px 8px;cursor:pointer}}
.te-sf-item{{display:flex;align-items:center;gap:6px;margin:2px 0;padding:2px 4px;border-radius:4px}}
.te-sf-item.te-sf-track-group{{flex-direction:column;align-items:stretch;gap:2px}}
.te-sf-track-head{{display:flex;align-items:center;gap:6px;min-width:0}}
.te-sf-track-controls{{display:flex;align-items:center;gap:6px;padding-left:22px;flex-wrap:wrap}}
.te-sf-item.te-sf-dragging{{opacity:0.45}}
.te-sf-item.te-sf-drop-target{{box-shadow:inset 0 -2px 0 #0969da}}
.te-sf-drag-handle{{cursor:grab;color:#8b949e;font-size:11px;line-height:1;user-select:none;flex:0 0 auto;padding:0 2px}}
.te-sf-drag-handle:active{{cursor:grabbing}}
.te-sf-list{{margin-top:4px}}
.te-sf-color{{width:28px;height:20px;padding:0;border:1px solid #ccc;border-radius:3px;cursor:pointer;flex:0 0 auto;background:none}}
.te-sf-scale{{margin:6px 0}}
.te-sf-scale input[type="range"]{{width:100%;display:block}}
.te-sf-window-ref{{margin:6px 0}}
.te-sf-window-ref input[type="checkbox"]{{width:14px;height:14px;margin:0;flex:0 0 auto;cursor:pointer}}
.te-sf-window-ref span{{line-height:1.3;cursor:default}}
.te-sf-hidden{{display:none}}
.te-sf-hint{{font-size:10px;color:#57606a;line-height:1.3;margin:2px 0 4px}}
.te-sf-track-label{{flex:1 1 auto;min-width:0;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;font-size:11px}}
.te-sf-ymax{{display:flex;align-items:center;gap:4px}}
.te-sf-ymax input{{width:74px;height:20px;padding:0 4px;border:1px solid #ccc;border-radius:3px}}
.te-sf-ymax label{{font-size:10px;color:#57606a}}
.te-sf-ymax-actions{{margin:4px 0 6px;padding-left:2px}}
.te-sf-ymax-actions button{{font-size:11px;padding:3px 10px;cursor:pointer}}
.te-sf-ymax-status{{font-size:10px;color:#1a7f37;margin-left:6px}}
@media print{{
  .te-sf-panel,.modebar{{display:none !important}}
  html,body{{margin:0 !important;padding:0 !important}}
}}
</style>
<script>
(function(){{
  const GENOME_MODE = {str(genome_mode).lower()};
  const SCALE_MAX = {scale_max_js};
  const XZOOM_MAX = {xzoom_max_js};
  const PREFS_PATH = {prefs_path_js};
  let prefsSaveTimer = null;

  // Plotly.js v3 may store large numeric arrays as binary objects
  // ({{dtype, bdata, _inputArray}}) instead of plain Arrays. The panel must
  // accept both forms or it never initializes on big TE density plots.
  function traceNumericArray(v){{
    if(v == null){{ return []; }}
    if(Array.isArray(v)){{
      return v.map((x) => Number(x || 0));
    }}
    if(typeof ArrayBuffer !== 'undefined' && ArrayBuffer.isView(v)){{
      return Array.from(v, (x) => Number(x || 0));
    }}
    if(typeof v === 'object'){{
      if(v._inputArray != null){{
        return traceNumericArray(v._inputArray);
      }}
      if(typeof v.length === 'number' && v.length > 0){{
        const out = new Array(v.length);
        for(let i = 0; i < v.length; i += 1){{
          out[i] = Number(v[i] || 0);
        }}
        return out;
      }}
    }}
    return [];
  }}

  function hasTraceValues(v){{
    if(v == null){{ return false; }}
    if(Array.isArray(v) || (typeof ArrayBuffer !== 'undefined' && ArrayBuffer.isView(v))){{
      return v.length > 0;
    }}
    if(typeof v === 'object'){{
      if(v._inputArray != null){{
        return hasTraceValues(v._inputArray);
      }}
      if(typeof v.length === 'number'){{
        return v.length > 0;
      }}
      if(v.dtype && v.bdata){{
        return true;
      }}
    }}
    return false;
  }}

  function traceArrayLength(v){{
    if(v == null){{ return 0; }}
    if(Array.isArray(v) || (typeof ArrayBuffer !== 'undefined' && ArrayBuffer.isView(v))){{
      return v.length;
    }}
    if(typeof v === 'object'){{
      if(v._inputArray != null){{
        return traceArrayLength(v._inputArray);
      }}
      if(typeof v.length === 'number'){{
        return v.length;
      }}
      if(v.dtype && v.bdata){{
        return traceNumericArray(v).length;
      }}
    }}
    return 0;
  }}

  function collectPlotPrefs(){{
    const out = {{
      bar_scale: 1,
      row_spacing: 1,
      diff_track_scale: 1,
      diff_color_a: '{TE_DIFF_COLOR_A}',
      diff_color_b: '{TE_DIFF_COLOR_B}',
      show_window_ref_line: true,
      window_ref_color: '{DEFAULT_WINDOW_REF_COLOR}',
      category_order: [],
      category_colors: {{}},
      category_hidden: [],
      meth_diff_track_scale: 1,
      meth_diff_color_a: '{METH_DIFF_COLOR_A}',
      meth_diff_color_b: '{METH_DIFF_COLOR_B}',
      meth_show_cg_total: false,
      meth_show_cg_covered: false,
      meth_show_cg_pct: true,
      meth_show_level: true,
      meth_show_diff: true,
      meth_color_cg_total: '{METH_COLOR_CG_TOTAL}',
      meth_color_cg_covered: '{METH_COLOR_CG_COVERED}',
      meth_color_cg_pct: '{METH_COLOR_CG_COVERAGE_PCT}',
      meth_color_level: '{METH_COLOR_METH_LEVEL}',
      gene_source_hidden: [],
      te_source_hidden: [],
      track_group_order: [],
      track_group_hidden: [],
      track_ymax_by_group: {{}},
      chromosome_sort: 'natural',
    }};
    const scaleEl = document.querySelector('#te-sf-scale');
    if(scaleEl){{ out.bar_scale = Number(scaleEl.value || 1); }}
    const rowEl = document.querySelector('#te-sf-rowsp');
    if(rowEl){{ out.row_spacing = Number(rowEl.value || 1); }}
    const diffEl = document.querySelector('#te-sf-diffscale');
    if(diffEl){{ out.diff_track_scale = Number(diffEl.value || 1); }}
    const diffColorA = document.querySelector('#te-sf-diffcolor-a');
    if(diffColorA){{ out.diff_color_a = diffColorA.value; }}
    const diffColorB = document.querySelector('#te-sf-diffcolor-b');
    if(diffColorB){{ out.diff_color_b = diffColorB.value; }}
    const refCb = document.querySelector('#te-sf-window-ref');
    if(refCb){{ out.show_window_ref_line = !!refCb.checked; }}
    const refColorEl = document.querySelector('#te-sf-window-ref-color');
    if(refColorEl){{ out.window_ref_color = refColorEl.value || '{DEFAULT_WINDOW_REF_COLOR}'; }}
    const methDiffEl = document.querySelector('#te-sf-methdiffscale');
    if(methDiffEl){{ out.meth_diff_track_scale = Number(methDiffEl.value || 1); }}
    const methDiffColorA = document.querySelector('#te-sf-methdiffcolor-a');
    if(methDiffColorA){{ out.meth_diff_color_a = methDiffColorA.value; }}
    const methDiffColorB = document.querySelector('#te-sf-methdiffcolor-b');
    if(methDiffColorB){{ out.meth_diff_color_b = methDiffColorB.value; }}
    const methShowCgTotal = document.querySelector('#te-sf-methshow-cgtotal');
    if(methShowCgTotal){{ out.meth_show_cg_total = !!methShowCgTotal.checked; }}
    const methShowCgCovered = document.querySelector('#te-sf-methshow-cgcovered');
    if(methShowCgCovered){{ out.meth_show_cg_covered = !!methShowCgCovered.checked; }}
    const methShowCgPct = document.querySelector('#te-sf-methshow-cgpct');
    if(methShowCgPct){{ out.meth_show_cg_pct = !!methShowCgPct.checked; }}
    const methShowLevel = document.querySelector('#te-sf-methshow-level');
    if(methShowLevel){{ out.meth_show_level = !!methShowLevel.checked; }}
    const methShowDiff = document.querySelector('#te-sf-methshow-diff');
    if(methShowDiff){{ out.meth_show_diff = !!methShowDiff.checked; }}
    const methColorCgTotal = document.querySelector('#te-sf-methcolor-cgtotal');
    if(methColorCgTotal){{ out.meth_color_cg_total = methColorCgTotal.value; }}
    const methColorCgCovered = document.querySelector('#te-sf-methcolor-cgcovered');
    if(methColorCgCovered){{ out.meth_color_cg_covered = methColorCgCovered.value; }}
    const methColorCgPct = document.querySelector('#te-sf-methcolor-cgpct');
    if(methColorCgPct){{ out.meth_color_cg_pct = methColorCgPct.value; }}
    const methColorLevel = document.querySelector('#te-sf-methcolor-level');
    if(methColorLevel){{ out.meth_color_level = methColorLevel.value; }}
    const trackOrder = [];
    document.querySelectorAll('#te-sf-tracks-list .te-sf-track-group').forEach((row) => {{
      const key = row.dataset.trackGroup || '';
      if(key){{ trackOrder.push(key); }}
    }});
    out.track_group_order = trackOrder;
    const trackHidden = [];
    document.querySelectorAll('#te-sf-tracks-list input[type="checkbox"]').forEach((cb) => {{
      if(!cb.checked && cb.id && String(cb.id).startsWith('te-sf-track-')){{
        trackHidden.push(String(cb.id).replace('te-sf-track-', ''));
      }}
    }});
    out.track_group_hidden = trackHidden;
    out.te_source_hidden = trackHidden.filter((k) => String(k).startsWith('te_src_'));
    out.gene_source_hidden = trackHidden.filter((k) => String(k).startsWith('gene_src_'));
    const trackYMax = {{}};
    document.querySelectorAll('#te-sf-tracks-list .te-sf-track-group').forEach((row) => {{
      const key = String(row.dataset.trackGroup || '');
      if(!key){{ return; }}
      const mainEl = row.querySelector('input[data-ymax-main="1"]');
      if(mainEl){{
        const v = Number(mainEl.value);
        if(Number.isFinite(v) && v > 0){{ trackYMax[key] = v; }}
      }}
      const fracEl = row.querySelector('input[data-ymax-sub="cg_coverage_pct"]');
      const levelEl = row.querySelector('input[data-ymax-sub="meth_level"]');
      const hasFrac = !!(fracEl && fracEl.value !== '');
      const hasLevel = !!(levelEl && levelEl.value !== '');
      if(hasFrac || hasLevel){{
        const bySub = {{}};
        if(hasFrac){{
          const fv = Number(fracEl.value);
          if(Number.isFinite(fv) && fv > 0 && fv <= 1){{ bySub.cg_coverage_pct = fv; }}
        }}
        if(hasLevel){{
          const lv = Number(levelEl.value);
          if(Number.isFinite(lv) && lv > 0 && lv <= 1){{ bySub.meth_level = lv; }}
        }}
        if(Object.keys(bySub).length){{
          trackYMax[key] = bySub;
        }}
      }}
    }});
    out.track_ymax_by_group = trackYMax;
    const gdPrefs = document.querySelector('.plotly-graph-div');
    if(gdPrefs && gdPrefs._teScaleMeta && gdPrefs._teScaleMeta.chromLayout
      && gdPrefs._teScaleMeta.chromLayout.chromosomeSort){{
      out.chromosome_sort = String(gdPrefs._teScaleMeta.chromLayout.chromosomeSort);
    }}
    const gd = document.querySelector('.plotly-graph-div');
    if(gd && gd._teXRangeZoom && gd._teXRangeZoom.length >= 2){{
      out.x_range = gd._teXRangeZoom.slice();
    }} else if(gd && gd.layout && gd.layout.xaxis && Array.isArray(gd.layout.xaxis.range)){{
      const xr = gd.layout.xaxis.range;
      if(
        gd._teXRange0
        && (
          Math.abs(Number(xr[0]) - Number(gd._teXRange0[0])) > 0.05
          || Math.abs(Number(xr[1]) - Number(gd._teXRange0[1])) > 0.05
        )
      ){{
        out.x_range = xr.slice();
      }}
    }}
    const order = [];
    const hidden = [];
    const colors = {{}};
    document.querySelectorAll('.te-sf-item').forEach((row) => {{
      const g = row.dataset.group;
      if(!g){{ return; }}
      order.push(g);
      const cb = row.querySelector('input[type="checkbox"]');
      if(cb && !cb.checked){{ hidden.push(g); }}
      const colorInput = row.querySelector('input[type="color"]');
      if(colorInput && colorInput.value){{ colors[g] = colorInput.value; }}
    }});
    out.category_order = order;
    out.category_hidden = hidden;
    out.category_colors = colors;
    return out;
  }}

  function schedulePrefsSave(){{
    if(!PREFS_PATH){{ return; }}
    if(prefsSaveTimer){{ clearTimeout(prefsSaveTimer); }}
    prefsSaveTimer = setTimeout(() => {{
      prefsSaveTimer = null;
      const payload = collectPlotPrefs();
      if(typeof window.repeatviewerSavePlotPrefs === 'function'){{
        window.repeatviewerSavePlotPrefs(payload);
      }}
    }}, 400);
  }}

  window.repeatviewerCollectPlotPrefs = collectPlotPrefs;
  function readLayoutMeta(gd){{
    const fl = gd._fullLayout && gd._fullLayout.meta;
    const ly = gd.layout && gd.layout.meta;
    return {{
      teGenomeScale: (fl && fl.teGenomeScale) || (ly && ly.teGenomeScale) || null,
      teDivergenceScale: (fl && fl.teDivergenceScale) || (ly && ly.teDivergenceScale) || null,
    }};
  }}

  function loadScaleMeta(gd){{
    const meta = readLayoutMeta(gd);
    if(GENOME_MODE && meta.teGenomeScale){{
      gd._teScaleMeta = meta.teGenomeScale;
    }} else if(!GENOME_MODE && meta.teDivergenceScale){{
      gd._teDivergenceMeta = meta.teDivergenceScale;
      if(Array.isArray(gd._teDivergenceMeta.xRange) && gd._teDivergenceMeta.xRange.length >= 2){{
        gd._teXRange0 = gd._teDivergenceMeta.xRange.slice();
      }}
    }}
  }}

  function init(){{
    const gd = document.querySelector('.plotly-graph-div');
    if(!gd || !Array.isArray(gd.data) || gd.data.length === 0){{ return false; }}
    if(gd._teSfInitialized){{
      loadScaleMeta(gd);
      if(
        GENOME_MODE
        && gd._teScaleMeta
        && !gd._teSfMetaApplied
        && typeof gd._teSfScheduleApply === 'function'
      ){{
        gd._teSfMetaApplied = true;
        const needsApply = typeof gd._teGenomePanelNeedsInitialApply === 'function'
          ? gd._teGenomePanelNeedsInitialApply(gd._teScaleMeta)
          : true;
        if(needsApply){{
          gd._teSfScheduleApply();
        }}
      }}
      return true;
    }}
    loadScaleMeta(gd);

    function isGeneTrackLabel(label){{
      const s = String(label || '').trim();
      if(!s){{ return false; }}
      const meta = gd._teScaleMeta;
      const kinds = meta && meta.chromLayout && meta.chromLayout.trackKindsByLabel;
      if(kinds && kinds[s] === 'gene'){{ return true; }}
      if(s.indexOf('(genes') >= 0){{ return true; }}
      if(meta && Array.isArray(meta.tracks)){{
        const tm = meta.tracks.find((t) => t.chrom === s || String(t.chrom) === s);
        if(tm && tm.isGeneTrack){{ return true; }}
      }}
      return false;
    }}

    gd._teIsGeneTrackLabel = isGeneTrackLabel;

    function barTrackKindFromTrace(tr){{
      if(!GENOME_MODE){{ return ''; }}
      const meta = gd._teScaleMeta;
      const cd = Array.isArray(tr.customdata) ? tr.customdata[0] : null;
      const rawLabel = cd && cd[0];
      if(rawLabel === undefined || rawLabel === null){{ return ''; }}
      const label = String(rawLabel).trim();
      if(isGeneTrackLabel(label)){{ return 'gene'; }}
      const kinds = meta && meta.chromLayout && meta.chromLayout.trackKindsByLabel;
      if(kinds && kinds[label]){{ return kinds[label]; }}
      if(meta && Array.isArray(meta.tracks)){{
        const tm = meta.tracks.find((t) => t.chrom === label || String(t.chrom) === label);
        if(tm && tm.isGeneTrack){{ return 'gene'; }}
        if(tm && tm.isMethylationTrack){{ return 'methylation'; }}
      }}
      return '';
    }}

    function registerPlainBarTrace(idx, tr){{
      if(!plainBarIndices.includes(idx)){{ plainBarIndices.push(idx); }}
      const keys = [];
      const grounds = [];
      const cd = Array.isArray(tr.customdata) ? tr.customdata : [];
      const bases = originalBase.get(idx) || [];
      for(let i = 0; i < Math.max(cd.length, (originalY.get(idx) || []).length); i += 1){{
        const row = cd[i] || [];
        const chrom = String(row[0] ?? '');
        const binIdx = row.length > 1 ? String(row[1] ?? i) : String(i);
        keys.push(chrom + '|' + binIdx);
        const ground = row.length > 4
          ? Number(row[4] ?? 0)
          : Number(bases[i] ?? 0);
        grounds.push(ground);
      }}
      if(keys.length){{
        pointKeys.set(idx, keys);
        pointGround.set(idx, grounds);
      }}
    }}

    const barTraceIndices = [];
    const originalY = new Map();
    const originalBase = new Map();
    const groups = [];
    const seenGroups = new Set();
    const barsByGroup = new Map();
    const plainBarIndices = [];
    const pointKeys = new Map();
    const pointGround = new Map();

    function registerBarTraceFromGd(idx, tr){{
      if(tr.type !== 'bar' || !hasTraceValues(tr.y)){{ return; }}
      if(!barTraceIndices.includes(idx)){{ barTraceIndices.push(idx); }}
      if(!originalY.has(idx)){{
        originalY.set(idx, traceNumericArray(tr.y));
        originalBase.set(idx, hasTraceValues(tr.base) ? traceNumericArray(tr.base) : []);
      }}
      const g = tr.legendgroup || tr.name || '';
      const barTrackKind = barTrackKindFromTrace(tr);
      if(g === '__te_diff__' || g === '__meth_diff__'){{
        if(!plainBarIndices.includes(idx)){{ plainBarIndices.push(idx); }}
        if(pointKeys.has(idx)){{ return; }}
        const keys = [];
        const cd = Array.isArray(tr.customdata) ? tr.customdata : [];
        for(let i = 0; i < cd.length; i += 1){{
          const row = cd[i] || [];
          const chrom = String(row[0] ?? '');
          const binIdx = String(row[1] ?? i);
          keys.push(chrom + '|' + binIdx);
        }}
        if(keys.length === 0){{
          for(let i = 0; i < originalY.get(idx).length; i += 1){{
            keys.push(String(idx) + '|' + String(i));
          }}
        }}
        pointKeys.set(idx, keys);
        return;
      }}
      if(barTrackKind === 'gene' || barTrackKind === 'methylation'){{
        registerPlainBarTrace(idx, tr);
        return;
      }}
      if(g){{
        if(!seenGroups.has(g)){{ seenGroups.add(g); groups.push(g); }}
        if(!barsByGroup.has(g)){{ barsByGroup.set(g, []); }}
        if(!barsByGroup.get(g).includes(idx)){{ barsByGroup.get(g).push(idx); }}
        if(pointKeys.has(idx)){{ return; }}
        const keys = [];
        const grounds = [];
        const cd = Array.isArray(tr.customdata) ? tr.customdata : [];
        for(let i = 0; i < cd.length; i += 1){{
          const row = cd[i] || [];
          const chrom = String(row[0] ?? '');
          const binIdx = String(row[1] ?? i);
          keys.push(chrom + '|' + binIdx);
          const ground = row.length > 4 ? Number(row[4] ?? 0) : Number(originalBase.get(idx)[i] ?? 0);
          grounds.push(ground);
        }}
        if(keys.length === 0){{
          for(let i = 0; i < originalY.get(idx).length; i += 1){{
            keys.push(String(idx) + '|' + String(i));
            grounds.push(Number(originalBase.get(idx)[i] ?? 0));
          }}
        }}
        pointKeys.set(idx, keys);
        pointGround.set(idx, grounds);
        return;
      }}
      registerPlainBarTrace(idx, tr);
    }}

    function rescanBarTracesIfNeeded(){{
      gd.data.forEach((tr, idx) => registerBarTraceFromGd(idx, tr));
    }}

    gd.data.forEach((tr, idx) => registerBarTraceFromGd(idx, tr));
    if(barTraceIndices.length === 0){{ return false; }}

    const hasCategoryFilters = groups.length > 0;
    const panel = document.createElement('div');
    panel.className = 'te-sf-panel';
    const filterBlock = '<div class="te-sf-title" title="Drag to move">'
      + (hasCategoryFilters ? '{filter_title}' : 'Track display')
      + '</div>'
      + '<div class="te-sf-actions">'
      + '<button type="button" id="te-sf-save-svg" title="Download the current chart view as SVG">Save SVG</button>'
      + '<button type="button" id="te-sf-save-pdf" title="Open the print dialog to save the current chart view as PDF">Save PDF</button>'
      + '</div>';
    panel.innerHTML = filterBlock
      + '<div id="te-sf-tracks-wrap" class="te-sf-hidden">'
      + '<div class="te-sf-title" style="font-size:11px;margin-top:6px">Tracks</div>'
      + '<div id="te-sf-tracks-list" class="te-sf-list"></div>'
      + '<div id="te-sf-ymax-actions" class="te-sf-ymax-actions te-sf-hidden">'
      + '<button type="button" id="te-sf-ymax-apply">Apply Y max</button>'
      + '<span id="te-sf-ymax-status" class="te-sf-ymax-status"></span>'
      + '</div></div>'
      + '<div class="te-sf-scale">{scale_label}: <span id="te-sf-scale-val">1.0x</span>'
      + '<input id="te-sf-scale" type="range" {scale_attrs} value="1"></div>'
      + {json.dumps(xzoom_panel_html)}
      + {json.dumps(genome_extra_html)}
      + (hasCategoryFilters ? {json.dumps(stack_hint_html)} : '');
    if(hasCategoryFilters){{
      const teMasterRow = document.createElement('div');
      teMasterRow.className = 'te-sf-item te-sf-te-master';
      teMasterRow.innerHTML = '<input type="checkbox" id="te-sf-te-master" checked>'
        + '<span>Show all TE types</span>';
      panel.appendChild(teMasterRow);
    }}
    const list = document.createElement('div');
    list.className = hasCategoryFilters ? 'te-sf-list' : 'te-sf-hidden te-sf-list';
    panel.appendChild(list);

    const checkByGroup = new Map();
    const colorByGroup = new Map();
    const selectedOrder = [];
    const checkByTrackGroup = new Map();
    const trackGroupOrder = [];
    let tracksListEl = panel.querySelector('#te-sf-tracks-list');

    if(!gd._teYRange0 && gd.layout && gd.layout.yaxis && Array.isArray(gd.layout.yaxis.range)){{
      gd._teYRange0 = gd.layout.yaxis.range.slice();
    }}
    if(!gd._teXRange0 && gd.layout && gd.layout.xaxis && Array.isArray(gd.layout.xaxis.range)){{
      gd._teXRange0 = gd.layout.xaxis.range.slice();
    }}
    if(gd._teXRange0 && gd._teXCenter === undefined){{
      gd._teXCenter = (Number(gd._teXRange0[0]) + Number(gd._teXRange0[1])) / 2;
    }}

    const metaForOrder = GENOME_MODE ? gd._teScaleMeta : gd._teDivergenceMeta;
    if(metaForOrder && Array.isArray(metaForOrder.categoryOrder) && metaForOrder.categoryOrder.length){{
      const orderIdx = new Map(metaForOrder.categoryOrder.map((g, i) => [g, i]));
      groups.sort((a, b) => {{
        const ia = orderIdx.has(a) ? orderIdx.get(a) : 1e9;
        const ib = orderIdx.has(b) ? orderIdx.get(b) : 1e9;
        if(ia !== ib){{ return ia - ib; }}
        return a.localeCompare(b, undefined, {{numeric: true, sensitivity: 'base'}});
      }});
    }}
    if(hasCategoryFilters){{
      groups.reverse();
    }}

    function normalizeHexColor(raw){{
      const s = String(raw || '').trim();
      if(/^#[0-9a-fA-F]{{6}}$/.test(s)){{ return s.toLowerCase(); }}
      const m = s.match(/rgb\\(\\s*(\\d+)\\s*,\\s*(\\d+)\\s*,\\s*(\\d+)\\s*\\)/);
      if(m){{
        const hx = (n) => Math.max(0, Math.min(255, Number(n))).toString(16).padStart(2, '0');
        return ('#' + hx(m[1]) + hx(m[2]) + hx(m[3])).toLowerCase();
      }}
      return '#7f7f7f';
    }}

    function defaultColorForGroup(group){{
      const meta = GENOME_MODE
        ? (gd._teScaleMeta && gd._teScaleMeta.categoryColors)
        : (gd._teDivergenceMeta && gd._teDivergenceMeta.categoryColors);
      if(meta && meta[group]){{ return normalizeHexColor(meta[group]); }}
      for(let i = 0; i < gd.data.length; i += 1){{
        const tr = gd.data[i];
        const g = tr.legendgroup || tr.name || '';
        if(g !== group){{ continue; }}
        const c = tr.marker && tr.marker.color;
        if(typeof c === 'string' && c){{ return normalizeHexColor(c); }}
      }}
      return '#7f7f7f';
    }}

    function diffTrackColorInputs(){{
      return {{
        a: panel.querySelector('#te-sf-diffcolor-a'),
        b: panel.querySelector('#te-sf-diffcolor-b'),
      }};
    }}

    function currentDiffTrackColors(){{
      const m = gd._teScaleMeta || {{}};
      const inputs = diffTrackColorInputs();
      return {{
        a: normalizeHexColor(
          (inputs.a && inputs.a.value) || m.diffColorA || m.savedDiffColorA || '{TE_DIFF_COLOR_A}',
        ),
        b: normalizeHexColor(
          (inputs.b && inputs.b.value) || m.diffColorB || m.savedDiffColorB || '{TE_DIFF_COLOR_B}',
        ),
      }};
    }}

    function currentMethDiffTrackColors(){{
      const m = gd._teScaleMeta || {{}};
      const pick = (id, fallback) => {{
        const el = panel.querySelector(id);
        return normalizeHexColor((el && el.value) || fallback);
      }};
      return {{
        a: pick('#te-sf-methdiffcolor-a', m.methDiffColorA || m.savedMethDiffColorA || '{METH_DIFF_COLOR_A}'),
        b: pick('#te-sf-methdiffcolor-b', m.methDiffColorB || m.savedMethDiffColorB || '{METH_DIFF_COLOR_B}'),
      }};
    }}

    function currentMethylationTrackColors(){{
      const m = gd._teScaleMeta || {{}};
      const pick = (id, fallback) => {{
        const el = panel.querySelector(id);
        return normalizeHexColor((el && el.value) || fallback);
      }};
      return {{
        cg_total: pick('#te-sf-methcolor-cgtotal', m.methColorCgTotal || '{METH_COLOR_CG_TOTAL}'),
        cg_covered: pick('#te-sf-methcolor-cgcovered', m.methColorCgCovered || '{METH_COLOR_CG_COVERED}'),
        cg_coverage_pct: pick('#te-sf-methcolor-cgpct', m.methColorCgPct || '{METH_COLOR_CG_COVERAGE_PCT}'),
        meth_level: pick('#te-sf-methcolor-level', m.methColorLevel || '{METH_COLOR_METH_LEVEL}'),
      }};
    }}

    function applyMethylationTrackColors(){{
      if(!GENOME_MODE || !gd._teScaleMeta || !gd._teScaleMeta.hasMethylationTracks){{ return Promise.resolve(); }}
      const m = gd._teScaleMeta;
      const colors = currentMethylationTrackColors();
      const labelToSub = new Map();
      (m.tracks || []).forEach((t) => {{
        if(t.isMethylationTrack){{ labelToSub.set(t.chrom, t.methylationSubtrack || 'cg_total'); }}
      }});
      const indices = [];
      const barColors = [];
      gd.data.forEach((tr, idx) => {{
        if(tr.type !== 'bar' || !Array.isArray(tr.customdata) || !tr.customdata.length){{ return; }}
        const trackLabel = tr.customdata[0] && tr.customdata[0][0];
        const sub = labelToSub.get(String(trackLabel || ''));
        if(!sub){{ return; }}
        const col = colors[sub] || colors.cg_total;
        const n = originalY.has(idx) ? originalY.get(idx).length : traceArrayLength(tr.y);
        indices.push(idx);
        barColors.push(new Array(n).fill(col));
      }});
      if(!indices.length){{ return Promise.resolve(); }}
      return Plotly.restyle(gd, {{'marker.color': barColors}}, indices);
    }}

    function applyDiffTrackColors(){{
      if(!GENOME_MODE || !gd._teScaleMeta || !gd._teScaleMeta.hasDiffTrack){{ return Promise.resolve(); }}
      const {{ a: colorA, b: colorB }} = currentDiffTrackColors();
      const indices = [];
      const colors = [];
      gd.data.forEach((tr, idx) => {{
        if(tr.type !== 'bar' || tr.legendgroup !== '__te_diff__'){{ return; }}
        const ys = originalY.get(idx) || traceNumericArray(tr.y);
        let isUp = true;
        for(let i = 0; i < ys.length; i += 1){{
          const y = Number(ys[i] || 0);
          if(y < 0){{ isUp = false; break; }}
          if(y > 0){{ isUp = true; break; }}
        }}
        indices.push(idx);
        colors.push(isUp ? colorA : colorB);
      }});
      if(!indices.length){{ return Promise.resolve(); }}
      return Plotly.restyle(gd, {{'marker.color': colors}}, indices);
    }}

    function applyMethDiffTrackColors(){{
      if(!GENOME_MODE || !gd._teScaleMeta || !gd._teScaleMeta.hasMethDiffTrack){{ return Promise.resolve(); }}
      const {{ a: colorA, b: colorB }} = currentMethDiffTrackColors();
      const indices = [];
      const colors = [];
      gd.data.forEach((tr, idx) => {{
        if(tr.type !== 'bar' || tr.legendgroup !== '__meth_diff__'){{ return; }}
        const ys = originalY.get(idx) || traceNumericArray(tr.y);
        let isUp = true;
        for(let i = 0; i < ys.length; i += 1){{
          const y = Number(ys[i] || 0);
          if(y < 0){{ isUp = false; break; }}
          if(y > 0){{ isUp = true; break; }}
        }}
        indices.push(idx);
        colors.push(isUp ? colorA : colorB);
      }});
      if(!indices.length){{ return Promise.resolve(); }}
      return Plotly.restyle(gd, {{'marker.color': colors}}, indices);
    }}

    function applyCategoryColors(){{
      if(!hasCategoryFilters || colorByGroup.size === 0){{ return Promise.resolve(); }}
      const indices = [];
      const colors = [];
      gd.data.forEach((tr, idx) => {{
        if(plainBarIndices.includes(idx)){{ return; }}
        const g = tr.legendgroup || tr.name || '';
        if(!g){{ return; }}
        const picker = colorByGroup.get(g);
        if(!picker){{ return; }}
        const col = picker.value;
        if(tr.type === 'bar'){{
          indices.push(idx);
          colors.push(col);
        }} else if(tr.type === 'scatter' && tr.mode === 'markers' && tr.showlegend){{
          indices.push(idx);
          colors.push(col);
        }}
      }});
      if(!indices.length){{ return Promise.resolve(); }}
      return Plotly.restyle(gd, {{'marker.color': colors}}, indices);
    }}

    function canonicalCategoryOrder(){{
      const meta = GENOME_MODE ? gd._teScaleMeta : gd._teDivergenceMeta;
      if(meta && Array.isArray(meta.categoryOrder) && meta.categoryOrder.length){{
        return meta.categoryOrder.slice();
      }}
      return groups.slice();
    }}

    function insertGroupInCanonicalOrder(g){{
      if(!g){{ return; }}
      const pos = selectedOrder.indexOf(g);
      if(pos >= 0){{ selectedOrder.splice(pos, 1); }}
      const built = canonicalCategoryOrder();
      const builtIdx = built.indexOf(g);
      let insertAt = selectedOrder.length;
      if(builtIdx >= 0){{
        insertAt = 0;
        for(let i = 0; i < builtIdx; i += 1){{
          if(selectedOrder.includes(built[i])){{ insertAt += 1; }}
        }}
      }}
      selectedOrder.splice(insertAt, 0, g);
      reorderListRowsToGroups(selectedOrder.slice());
    }}

    function syncSelectedOrderFromList(){{
      selectedOrder.length = 0;
      list.querySelectorAll('.te-sf-item').forEach((row) => {{
        const g = row.dataset.group || '';
        if(g){{ selectedOrder.push(g); }}
      }});
    }}

    function syncTeMasterCheckbox(){{
      const teMaster = panel.querySelector('#te-sf-te-master');
      if(!teMaster || !hasCategoryFilters){{ return; }}
      let nOn = 0;
      checkByGroup.forEach((cb) => {{
        if(cb.checked){{ nOn += 1; }}
      }});
      const n = checkByGroup.size;
      teMaster.indeterminate = nOn > 0 && nOn < n;
      teMaster.checked = nOn === n;
    }}

    function reorderListRowsToGroups(order){{
      order.forEach((g) => {{
        list.querySelectorAll('.te-sf-item').forEach((row) => {{
          if(row.dataset.group === g){{ list.appendChild(row); }}
        }});
      }});
      syncSelectedOrderFromList();
    }}

    function initRowDragDrop(row){{
      const handle = row.querySelector('.te-sf-drag-handle');
      if(!handle){{ return; }}
      handle.addEventListener('dragstart', (e) => {{
        row.classList.add('te-sf-dragging');
        e.dataTransfer.effectAllowed = 'move';
        e.dataTransfer.setData('text/plain', row.dataset.group || '');
        if(e.dataTransfer.setDragImage){{
          e.dataTransfer.setDragImage(row, 12, 12);
        }}
      }});
      handle.addEventListener('dragend', () => {{
        row.classList.remove('te-sf-dragging');
        list.querySelectorAll('.te-sf-drop-target').forEach((el) => {{
          el.classList.remove('te-sf-drop-target');
        }});
      }});
      row.addEventListener('dragover', (e) => {{
        e.preventDefault();
        e.dataTransfer.dropEffect = 'move';
        const dragging = list.querySelector('.te-sf-item.te-sf-dragging');
        if(!dragging || dragging === row){{ return; }}
        row.classList.add('te-sf-drop-target');
      }});
      row.addEventListener('dragleave', () => {{
        row.classList.remove('te-sf-drop-target');
      }});
      row.addEventListener('drop', (e) => {{
        e.preventDefault();
        row.classList.remove('te-sf-drop-target');
        const dragging = list.querySelector('.te-sf-item.te-sf-dragging');
        if(!dragging || dragging === row){{ return; }}
        const rect = row.getBoundingClientRect();
        const before = e.clientY < rect.top + rect.height / 2;
        if(before){{
          list.insertBefore(dragging, row);
        }} else {{
          list.insertBefore(dragging, row.nextSibling);
        }}
        syncSelectedOrderFromList();
        scheduleApply();
        schedulePrefsSave();
      }});
    }}

    groups.forEach((g) => {{
      const row = document.createElement('div');
      row.className = 'te-sf-item';
      row.dataset.group = g;
      const handle = document.createElement('span');
      handle.className = 'te-sf-drag-handle';
      handle.textContent = '⋮⋮';
      handle.title = 'Drag to reorder';
      handle.draggable = true;
      const cb = document.createElement('input');
      cb.type = 'checkbox';
      cb.checked = true;
      cb.dataset.group = g;
      const colorInput = document.createElement('input');
      colorInput.type = 'color';
      colorInput.className = 'te-sf-color';
      colorInput.value = defaultColorForGroup(g);
      colorInput.title = 'Colour for ' + g;
      colorInput.addEventListener('input', () => {{
        applyCategoryColors();
      }});
      const txt = document.createElement('span');
      txt.textContent = g;
      row.appendChild(handle);
      row.appendChild(cb);
      row.appendChild(colorInput);
      row.appendChild(txt);
      list.appendChild(row);
      initRowDragDrop(row);
      checkByGroup.set(g, cb);
      colorByGroup.set(g, colorInput);
      selectedOrder.push(g);
    }});

    if(GENOME_MODE){{
      const m = gd._teScaleMeta;
      const xr = gd.layout && gd.layout.xaxis && gd.layout.xaxis.range;
      if(m && Array.isArray(m.xRange) && m.xRange.length >= 2){{
        gd._teXRange0 = m.xRange.slice();
      }} else if(Array.isArray(xr) && xr.length >= 2){{
        gd._teXRange0 = xr.slice();
      }}
      gd._teXRangeZoom = null;
      gd._teUserXZoom = false;
      gd._teRelayoutSync = false;
      gd._teApplyDepth = 0;
      if(m && Array.isArray(m.savedXRange) && m.savedXRange.length >= 2){{
        gd._teXRangeZoom = m.savedXRange.slice();
        gd._teUserXZoom = true;
      }}
    }}

    function xRangesClose(a, b){{
      if(!Array.isArray(a) || !Array.isArray(b) || a.length < 2 || b.length < 2){{ return false; }}
      return Math.abs(Number(a[0]) - Number(b[0])) < 0.05
        && Math.abs(Number(a[1]) - Number(b[1])) < 0.05;
    }}

    function readXZoomScale(){{
      const el = panel.querySelector('#te-sf-xzoom');
      let zoom = Number(el ? el.value : 1);
      if(!Number.isFinite(zoom) || zoom <= 0){{ zoom = 1; }}
      return Math.min(XZOOM_MAX, Math.max({PLOT_XZOOM_MIN}, zoom));
    }}

    function snapDivergenceRange(r0, r1, meta){{
      const bw = Number(meta && meta.binWidth);
      if(!Number.isFinite(bw) || bw <= 0){{ return [r0, r1]; }}
      const lo = Math.floor(Number(r0) / bw) * bw;
      const hi = Math.ceil(Number(r1) / bw) * bw;
      return [lo, Math.max(lo + bw, hi)];
    }}

    function computeXRangeFromZoom(zoomScale){{
      if(!gd._teXRange0 || gd._teXRange0.length < 2){{ return null; }}
      const full0 = Number(gd._teXRange0[0]);
      const full1 = Number(gd._teXRange0[1]);
      const fullSpan = full1 - full0;
      if(!Number.isFinite(fullSpan) || fullSpan <= 0){{ return null; }}
      if(zoomScale <= 1.001){{ return null; }}
      const visibleSpan = fullSpan / zoomScale;
      const center = Number.isFinite(gd._teXCenter)
        ? Number(gd._teXCenter)
        : (full0 + full1) / 2;
      let r0 = center - visibleSpan / 2;
      let r1 = center + visibleSpan / 2;
      if(r0 < full0){{
        r0 = full0;
        r1 = Math.min(full1, r0 + visibleSpan);
      }}
      if(r1 > full1){{
        r1 = full1;
        r0 = Math.max(full0, r1 - visibleSpan);
      }}
      if(GENOME_MODE && gd._teScaleMeta){{
        return snapDisplayRangeToBins(r0, r1, gd._teScaleMeta);
      }}
      if(!GENOME_MODE && gd._teDivergenceMeta){{
        return snapDivergenceRange(r0, r1, gd._teDivergenceMeta);
      }}
      return [r0, r1];
    }}

    function updateXRangeFromSlider(){{
      const el = panel.querySelector('#te-sf-xzoom');
      if(!el){{ return; }}
      const zoomScale = readXZoomScale();
      gd._teXRangeZoom = computeXRangeFromZoom(zoomScale);
    }}

    function syncXZoomSliderFromRange(r0, r1){{
      if(!gd._teXRange0 || gd._teXRange0.length < 2){{ return; }}
      const fullSpan = Number(gd._teXRange0[1]) - Number(gd._teXRange0[0]);
      const visSpan = Number(r1) - Number(r0);
      if(!Number.isFinite(fullSpan) || fullSpan <= 0 || !Number.isFinite(visSpan) || visSpan <= 0){{ return; }}
      gd._teXCenter = (Number(r0) + Number(r1)) / 2;
      const zoom = Math.min(XZOOM_MAX, Math.max({PLOT_XZOOM_MIN}, fullSpan / visSpan));
      const el = panel.querySelector('#te-sf-xzoom');
      const valEl = panel.querySelector('#te-sf-xzoom-val');
      if(el){{ el.value = String(zoom); }}
      if(valEl){{ valEl.textContent = zoom.toFixed(1) + 'x'; }}
    }}

    function applyXAxisZoomLayout(layoutUpdate){{
      if(panel.querySelector('#te-sf-xzoom')){{
        updateXRangeFromSlider();
      }}
      if(!gd._teXRange0){{ return; }}
      const range = gd._teXRangeZoom ? gd._teXRangeZoom.slice() : gd._teXRange0.slice();
      layoutUpdate['xaxis.range'] = range;
      layoutUpdate['xaxis.autorange'] = false;
      if(GENOME_MODE){{
        layoutUpdate['xaxis.visible'] = true;
        layoutUpdate['xaxis.showticklabels'] = true;
        layoutUpdate['xaxis.showline'] = true;
        layoutUpdate['xaxis.title.text'] = 'Mbp';
        layoutUpdate['margin.b'] = 56;
        const mLeft = gd._teScaleMeta && gd._teScaleMeta.marginLeft;
        if(mLeft){{ layoutUpdate['margin.l'] = mLeft; }}
        const m = gd._teScaleMeta;
        const domStart = m && Number.isFinite(m.xDomainStart) ? Number(m.xDomainStart) : 0.042;
        layoutUpdate['xaxis.domain'] = [domStart, 1];
        if(!gd._teXRangeZoom && m && Array.isArray(m.xTickvals) && Array.isArray(m.xTicktext)){{
          layoutUpdate['xaxis.tickmode'] = 'array';
          layoutUpdate['xaxis.tickvals'] = m.xTickvals;
          layoutUpdate['xaxis.ticktext'] = m.xTicktext;
        }} else if(gd._teXRangeZoom && m){{
          const ticks = zoomBinTickvals(range[0], range[1], m);
          layoutUpdate['xaxis.tickmode'] = 'array';
          layoutUpdate['xaxis.tickvals'] = ticks.tickvals;
          layoutUpdate['xaxis.ticktext'] = ticks.ticktext;
        }}
      }}
    }}

    function snapDisplayRangeToBins(r0, r1, meta){{
      const pad = Number(meta && meta.trackPadMbp) || 0;
      const binSpan = Number(meta && meta.binSpanMbp);
      if(!Number.isFinite(binSpan) || binSpan <= 0){{ return [r0, r1]; }}
      const chromStart = Number(meta.chromStartMbp) || 0;
      const viewHi = Number(meta.viewHiMbp);
      let gen0 = Number(r0) - pad;
      let gen1 = Number(r1) - pad;
      const snappedGen0 = Math.floor(gen0 / binSpan) * binSpan;
      const snappedGen1 = Math.ceil(gen1 / binSpan) * binSpan;
      const clampLo = chromStart;
      const clampHi = Number.isFinite(viewHi) ? viewHi : snappedGen1;
      const g0 = Math.max(clampLo, snappedGen0);
      const g1 = Math.min(clampHi, Math.max(g0 + binSpan, snappedGen1));
      return [g0 + pad, g1 + pad];
    }}

    function zoomBinTickvals(r0, r1, meta){{
      const pad = Number(meta && meta.trackPadMbp) || 0;
      const gen0 = Math.max(Number(meta && meta.chromStartMbp) || 0, Number(r0) - pad);
      const gen1 = Math.max(gen0, Number(r1) - pad);
      const span = Math.max(0.5, gen1 - gen0);
      let step;
      if(span > 40){{ step = 10; }}
      else if(span > 15){{ step = 5; }}
      else {{ step = Math.max(0.5, Math.round((span / 6) * 10) / 10); }}
      const start = Math.floor(gen0 / step) * step;
      const genomic = [];
      for(let g = start; g <= gen1 + step * 0.001; g += step){{
        genomic.push(g);
      }}
      return {{
        tickvals: genomic.map((t) => t + pad),
        ticktext: genomic.map((t) => (
          Math.abs(t - Math.round(t)) < 0.05 ? String(Math.round(t)) : t.toFixed(1)
        )),
      }};
    }}

    function readXRangeFromRelayoutEvent(ev){{
      if(!ev){{ return null; }}
      let r0;
      let r1;
      if(Array.isArray(ev['xaxis.range']) && ev['xaxis.range'].length >= 2){{
        r0 = Number(ev['xaxis.range'][0]);
        r1 = Number(ev['xaxis.range'][1]);
      }} else if(ev['xaxis.range[0]'] !== undefined || ev['xaxis.range[1]'] !== undefined){{
        r0 = Number(ev['xaxis.range[0]'] ?? gd.layout?.xaxis?.range?.[0]);
        r1 = Number(ev['xaxis.range[1]'] ?? gd.layout?.xaxis?.range?.[1]);
      }} else {{
        return null;
      }}
      if(!Number.isFinite(r0) || !Number.isFinite(r1)){{ return null; }}
      return [r0, r1];
    }}

    function captureXRangeFromRelayoutEvent(ev){{
      const newRange = readXRangeFromRelayoutEvent(ev);
      if(!newRange){{ return; }}
      if(gd._teXRange0 && xRangesClose(newRange, gd._teXRange0)){{
        gd._teXRangeZoom = null;
        gd._teUserXZoom = false;
        return;
      }}
      const r0 = newRange[0];
      const r1 = newRange[1];
      gd._teUserXZoom = true;
      gd._teXCenter = (r0 + r1) / 2;
      if(GENOME_MODE && gd._teScaleMeta){{
        gd._teXRangeZoom = snapDisplayRangeToBins(r0, r1, gd._teScaleMeta).slice();
      }} else if(gd._teDivergenceMeta){{
        gd._teXRangeZoom = snapDivergenceRange(r0, r1, gd._teDivergenceMeta).slice();
      }} else {{
        gd._teXRangeZoom = [r0, r1];
      }}
      syncXZoomSliderFromRange(gd._teXRangeZoom[0], gd._teXRangeZoom[1]);
    }}

    function preserveXRange(layoutUpdate){{
      if(gd._teUserXZoom && gd._teXRangeZoom && gd._teXRangeZoom.length >= 2){{
        const range = gd._teXRangeZoom.slice();
        layoutUpdate['xaxis.range'] = range;
        layoutUpdate['xaxis.autorange'] = false;
        const m = gd._teScaleMeta;
        if(GENOME_MODE && m){{
          const ticks = zoomBinTickvals(range[0], range[1], m);
          layoutUpdate['xaxis.tickmode'] = 'array';
          layoutUpdate['xaxis.tickvals'] = ticks.tickvals;
          layoutUpdate['xaxis.ticktext'] = ticks.ticktext;
        }}
        return;
      }}
      applyXAxisZoomLayout(layoutUpdate);
    }}

    function refreshPlotLayout(){{
      if(!GENOME_MODE || typeof Plotly === 'undefined'){{ return Promise.resolve(); }}
      const layoutUpdate = {{}};
      preserveXRange(layoutUpdate);
      if(!Object.keys(layoutUpdate).length){{ return Promise.resolve(); }}
      return safeRelayout(layoutUpdate).then(() => {{
        if(Plotly.Plots && typeof Plotly.Plots.resize === 'function'){{
          Plotly.Plots.resize(gd);
        }}
        window.dispatchEvent(new Event('resize'));
      }});
    }}

    function safeRelayout(layoutUpdate){{
      gd._teRelayoutSync = true;
      return Plotly.relayout(gd, layoutUpdate).finally(() => {{
        gd._teRelayoutSync = false;
      }});
    }}

    function formatDensityTick(value, teLengthMode){{
      if(!Number.isFinite(value) || value <= 0){{ return '0'; }}
      if(teLengthMode){{
        if(value < 1000){{ return String(Math.round(value)) + ' bp'; }}
        const kb = value / 1000;
        if(kb >= 100){{ return kb.toFixed(0) + ' kb'; }}
        if(kb >= 10){{
          return (Math.abs(kb - Math.round(kb)) < 0.05 ? kb.toFixed(0) : kb.toFixed(1)) + ' kb';
        }}
        if(kb >= 1){{
          return (Math.abs(kb - Math.round(kb)) < 0.05 ? kb.toFixed(0) : kb.toFixed(1)) + ' kb';
        }}
        return kb.toFixed(2) + ' kb';
      }}
      return String(Math.round(value));
    }}

    function formatFractionAxisTick(value){{
      if(!Number.isFinite(value) || value <= 0){{ return '0%'; }}
      const pct = value * 100;
      if(pct >= 10){{ return pct.toFixed(0) + '%'; }}
      if(pct >= 1){{
        return (Math.abs(pct - Math.round(pct)) < 0.05 ? pct.toFixed(0) : pct.toFixed(1)) + '%';
      }}
      return pct.toFixed(1) + '%';
    }}

    function formatUnitIntervalAxisTick(value){{
      if(!Number.isFinite(value) || value <= 0){{ return '0'; }}
      if(value < 1){{ return value.toFixed(2); }}
      return value.toFixed(1);
    }}

    function formatDiffAxisTick(frac, diffMax, teLengthMode){{
      if(Math.abs(frac - 0.5) < 1e-6){{ return '0'; }}
      const sign = frac > 0.5 ? 1 : -1;
      const mag = diffMax * Math.abs(frac - 0.5) * 2.0;
      const text = formatDensityTick(mag, !!teLengthMode);
      return sign < 0 ? '-' + text : text;
    }}

    function formatMethDiffAxisTick(frac, diffMax){{
      if(Math.abs(frac - 0.5) < 1e-6){{ return '0'; }}
      const sign = frac > 0.5 ? 1 : -1;
      const mag = diffMax * Math.abs(frac - 0.5) * 2.0;
      const text = formatUnitIntervalAxisTick(mag);
      return sign < 0 ? '-' + text : text;
    }}

    function diffTrackHalf(cl, rowMult, diffScale){{
      const m = gd._teScaleMeta;
      const base = m && Number(m.diffBarHalf) > 0
        ? Number(m.diffBarHalf)
        : Number(m && m.maxBarHeight || 0.36) * 0.36;
      const mult = Number.isFinite(diffScale) && diffScale > 0 ? diffScale : 1;
      return Math.max(0.12, base * (rowMult || 1) * mult);
    }}

    function greyTrackHalf(label, idx, cl, rowMult, diffScale, methDiffScale){{
      const kinds = cl && cl.trackKindsByLabel ? cl.trackKindsByLabel : {{}};
      if(kinds[label] === 'te_diff'){{
        return diffTrackHalf(cl, rowMult, diffScale);
      }}
      if(kinds[label] === 'meth_diff'){{
        return methDiffTrackHalf(cl, rowMult, methDiffScale);
      }}
      return cl ? cl.trackHalf : 0.1;
    }}

    function methylationSubtrackShowId(sub){{
      const map = {{
        cg_total: 'te-sf-methshow-cgtotal',
        cg_covered: 'te-sf-methshow-cgcovered',
        cg_coverage_pct: 'te-sf-methshow-cgpct',
        meth_level: 'te-sf-methshow-level',
      }};
      return map[String(sub || '')] || '';
    }}

    function methylationSubtrackTypeVisible(sub){{
      const id = methylationSubtrackShowId(sub);
      if(!id){{ return true; }}
      const cb = panel.querySelector('#' + id);
      if(!cb){{ return true; }}
      return !!cb.checked;
    }}

    function chromFromTrackLabel(label){{
      const s = String(label || '');
      const idx = s.indexOf(' (');
      return idx >= 0 ? s.slice(0, idx) : s;
    }}

    function naturalSortKey(name){{
      return String(name).split(/(\\d+)/).map((part, i) => (
        i % 2 ? Number(part) : part.toLowerCase()
      ));
    }}

    function compareNaturalChromNames(a, b){{
      const ka = naturalSortKey(a);
      const kb = naturalSortKey(b);
      const n = Math.min(ka.length, kb.length);
      for(let i = 0; i < n; i += 1){{
        if(ka[i] < kb[i]){{ return -1; }}
        if(ka[i] > kb[i]){{ return 1; }}
      }}
      return ka.length - kb.length;
    }}

    function sortChromosomeList(chroms, sizes, mode){{
      const list = chroms.slice();
      if(mode === 'size_desc'){{
        return list.sort((a, b) => {{
          const da = Number(sizes[a] || 0);
          const db = Number(sizes[b] || 0);
          if(db !== da){{ return db - da; }}
          return compareNaturalChromNames(a, b);
        }});
      }}
      if(mode === 'size_asc'){{
        return list.sort((a, b) => {{
          const da = Number(sizes[a] || 0);
          const db = Number(sizes[b] || 0);
          if(da !== db){{ return da - db; }}
          return compareNaturalChromNames(a, b);
        }});
      }}
      return list.sort(compareNaturalChromNames);
    }}

    function applyChromosomeSort(mode){{
      const m = gd._teScaleMeta;
      if(!m || !m.chromLayout){{ return; }}
      const cl = m.chromLayout;
      const sizes = cl.chromSizes || {{}};
      const baseChroms = Array.isArray(cl.naturalChromOrder) && cl.naturalChromOrder.length
        ? cl.naturalChromOrder.slice()
        : (Array.isArray(cl.chromOrder) ? [...new Set(cl.chromOrder)] : []);
      cl.chromOrder = sortChromosomeList(baseChroms, sizes, mode);
      cl.chromosomeSort = mode;
      syncChromLayoutOrderFromTrackGroups();
      scheduleApply();
      schedulePrefsSave();
    }}

    function trackGroupVisible(groupKey){{
      if(!groupKey){{ return true; }}
      const cb = checkByTrackGroup.get(String(groupKey));
      if(cb){{ return !!cb.checked; }}
      const el = panel.querySelector('#te-sf-track-' + String(groupKey));
      if(el){{ return !!el.checked; }}
      return true;
    }}

    function syncTrackGroupOrderFromList(){{
      trackGroupOrder.length = 0;
      if(!tracksListEl){{ return; }}
      tracksListEl.querySelectorAll('.te-sf-track-group').forEach((row) => {{
        const key = row.dataset.trackGroup || '';
        if(key){{ trackGroupOrder.push(key); }}
      }});
    }}

    function reorderTrackLabelsByGroupOrder(baseOrder, cl, order){{
      if(!cl || !Array.isArray(baseOrder) || !order.length){{ return baseOrder.slice(); }}
      const byGroup = cl.trackGroupKeyByLabel || {{}};
      const chromOrder = Array.isArray(cl.chromOrder) && cl.chromOrder.length
        ? cl.chromOrder.slice()
        : [...new Set(baseOrder.map((lbl) => chromFromTrackLabel(lbl)))];
      const newOrder = [];
      chromOrder.forEach((chrom) => {{
        order.forEach((groupKey) => {{
          baseOrder.forEach((lbl) => {{
            if(chromFromTrackLabel(lbl) === chrom && byGroup[lbl] === groupKey){{
              newOrder.push(lbl);
            }}
          }});
        }});
      }});
      baseOrder.forEach((lbl) => {{
        if(!newOrder.includes(lbl)){{ newOrder.push(lbl); }}
      }});
      return newOrder;
    }}

    function syncChromLayoutOrderFromTrackGroups(){{
      const m = gd._teScaleMeta;
      if(!m || !m.chromLayout || !Array.isArray(m._baseTrackOrder)){{ return; }}
      m.chromLayout.order = reorderTrackLabelsByGroupOrder(
        m._baseTrackOrder,
        m.chromLayout,
        trackGroupOrder,
      );
    }}

    function initTrackGroupDragDrop(row){{
      if(!tracksListEl){{ return; }}
      const handle = row.querySelector('.te-sf-drag-handle');
      if(!handle){{ return; }}
      handle.addEventListener('dragstart', (e) => {{
        row.classList.add('te-sf-dragging');
        if(e.dataTransfer){{
          e.dataTransfer.effectAllowed = 'move';
          e.dataTransfer.setData('text/plain', row.dataset.trackGroup || '');
        }}
      }});
      handle.addEventListener('dragend', () => {{
        row.classList.remove('te-sf-dragging');
        tracksListEl.querySelectorAll('.te-sf-drop-target').forEach((el) => {{
          el.classList.remove('te-sf-drop-target');
        }});
      }});
      row.addEventListener('dragover', (e) => {{
        e.preventDefault();
        const dragging = tracksListEl.querySelector('.te-sf-item.te-sf-dragging');
        if(!dragging || dragging === row){{ return; }}
        row.classList.add('te-sf-drop-target');
      }});
      row.addEventListener('dragleave', () => {{
        row.classList.remove('te-sf-drop-target');
      }});
      row.addEventListener('drop', (e) => {{
        e.preventDefault();
        row.classList.remove('te-sf-drop-target');
        const dragging = tracksListEl.querySelector('.te-sf-item.te-sf-dragging');
        if(!dragging || dragging === row){{ return; }}
        const rect = row.getBoundingClientRect();
        const before = e.clientY < rect.top + rect.height / 2;
        if(before){{
          tracksListEl.insertBefore(dragging, row);
        }} else {{
          tracksListEl.insertBefore(dragging, row.nextSibling);
        }}
        syncTrackGroupOrderFromList();
        syncChromLayoutOrderFromTrackGroups();
        scheduleApply();
        schedulePrefsSave();
      }});
    }}

    function shortTrackLabel(label){{
      const s = String(label || '');
      const parts = s.split(/[/\\\\]/);
      const base = parts[parts.length - 1] || s;
      if(base.length > 52){{ return base.slice(0, 49) + '…'; }}
      return base;
    }}

    function collectTrackYMaxFromPanel(){{
      const byGroup = {{}};
      if(!tracksListEl){{ return byGroup; }}
      tracksListEl.querySelectorAll('.te-sf-track-group').forEach((row) => {{
        const key = String(row.dataset.trackGroup || '');
        if(!key){{ return; }}
        const mainEl = row.querySelector('input[data-ymax-main="1"]');
        if(mainEl && mainEl.value !== ''){{
          const v = Number(mainEl.value);
          if(Number.isFinite(v) && v > 0){{
            if(key === 'meth_diff'){{
              byGroup[key] = Math.min(1, Math.max(0.0001, v));
            }} else {{
              byGroup[key] = Math.max(1, Math.round(v));
            }}
            return;
          }}
        }}
        const fracEl = row.querySelector('input[data-ymax-sub="cg_coverage_pct"]');
        const levelEl = row.querySelector('input[data-ymax-sub="meth_level"]');
        const sub = {{}};
        if(fracEl && fracEl.value !== ''){{
          const v = Number(fracEl.value);
          if(Number.isFinite(v) && v > 0 && v <= 1){{ sub.cg_coverage_pct = v; }}
        }}
        if(levelEl && levelEl.value !== ''){{
          const v = Number(levelEl.value);
          if(Number.isFinite(v) && v > 0 && v <= 1){{ sub.meth_level = v; }}
        }}
        if(Object.keys(sub).length){{
          byGroup[key] = sub;
        }}
      }});
      return byGroup;
    }}

    function commitTrackYMaxFromPanel(){{
      const m = gd._teScaleMeta;
      if(!m){{ return false; }}
      const byGroup = collectTrackYMaxFromPanel();
      m._appliedTrackYMaxByGroup = byGroup;
      m._panelTrackYMaxByGroup = byGroup;
      m.trackYMaxByGroup = byGroup;
      return true;
    }}

    function applyTrackYMax(){{
      if(!commitTrackYMaxFromPanel()){{ return; }}
      const m = gd._teScaleMeta;
      if(m){{ m._yMaxRescaleActive = true; }}
      gd._teApplyFastPath = '';
      invalidateLayoutCache();
      scheduleApply();
      schedulePrefsSave();
      const status = panel.querySelector('#te-sf-ymax-status');
      if(status){{
        status.textContent = 'Applied';
        window.setTimeout(() => {{ status.textContent = ''; }}, 2000);
      }}
    }}

    function buildTracksPanel(){{
      const m = gd._teScaleMeta;
      const wrap = panel.querySelector('#te-sf-tracks-wrap');
      tracksListEl = panel.querySelector('#te-sf-tracks-list');
      if(!wrap || !tracksListEl || !m || !Array.isArray(m.panelTrackGroups) || !m.panelTrackGroups.length){{
        if(wrap){{ wrap.classList.add('te-sf-hidden'); }}
        return;
      }}
      wrap.classList.remove('te-sf-hidden');
      const ymaxActions = panel.querySelector('#te-sf-ymax-actions');
      if(ymaxActions){{ ymaxActions.classList.remove('te-sf-hidden'); }}
      tracksListEl.innerHTML = '';
      checkByTrackGroup.clear();
      trackGroupOrder.length = 0;

      const hidden = new Set();
      if(Array.isArray(m.savedTrackGroupHidden) && m.savedTrackGroupHidden.length){{
        m.savedTrackGroupHidden.forEach((k) => {{ hidden.add(String(k)); }});
      }} else {{
        (m.savedTeSourceHidden || []).forEach((k) => {{ hidden.add(String(k)); }});
        (m.savedGeneSourceHidden || []).forEach((k) => {{ hidden.add(String(k)); }});
      }}

      const validKeys = new Set(m.panelTrackGroups.map((g) => String(g.key)));
      let order = Array.isArray(m.savedTrackGroupOrder) && m.savedTrackGroupOrder.length
        ? m.savedTrackGroupOrder.map(String).filter((k) => validKeys.has(k))
        : m.panelTrackGroups.map((g) => String(g.key));
      m.panelTrackGroups.forEach((g) => {{
        const k = String(g.key);
        if(!order.includes(k)){{ order.push(k); }}
      }});

      const groupByKey = new Map(m.panelTrackGroups.map((g) => [String(g.key), g]));
      order.forEach((key) => {{
        const src = groupByKey.get(key);
        if(!src){{ return; }}
        const row = document.createElement('div');
        row.className = 'te-sf-item te-sf-track-group';
        row.dataset.trackGroup = key;
        const head = document.createElement('div');
        head.className = 'te-sf-track-head';
        const controls = document.createElement('div');
        controls.className = 'te-sf-track-controls';
        const handle = document.createElement('span');
        handle.className = 'te-sf-drag-handle';
        handle.textContent = '⋮⋮';
        handle.title = 'Drag to reorder tracks';
        handle.draggable = true;
        const cb = document.createElement('input');
        cb.type = 'checkbox';
        cb.id = 'te-sf-track-' + key;
        cb.checked = !hidden.has(key);
        const fullLabel = String(src.label || key);
        const txt = document.createElement('span');
        txt.className = 'te-sf-track-label';
        txt.textContent = shortTrackLabel(fullLabel);
        txt.title = fullLabel;
        const savedYMaxByGroup = m.savedTrackYMaxByGroup || m.trackYMaxByGroup || {{}};
        const savedYMax = savedYMaxByGroup[key];
        const defaultYMax = src.defaultYMax;
        function resolveYMaxValue(saved, fallback){{
          if(saved !== undefined && saved !== null && saved !== ''){{ return saved; }}
          if(fallback !== undefined && fallback !== null && fallback !== ''){{ return fallback; }}
          return '';
        }}
        function makeYMaxInput(opts){{
          const wrap = document.createElement('span');
          wrap.className = 'te-sf-ymax';
          const lbl = document.createElement('label');
          lbl.textContent = opts.label;
          const input = document.createElement('input');
          input.type = 'number';
          input.placeholder = 'max Y';
          input.step = opts.step;
          if(opts.min !== undefined){{ input.min = String(opts.min); }}
          if(opts.max !== undefined){{ input.max = String(opts.max); }}
          input.dataset.ymaxMain = opts.main ? '1' : '0';
          if(opts.subKey){{ input.dataset.ymaxSub = opts.subKey; }}
          const resolved = resolveYMaxValue(opts.saved, opts.fallback);
          if(resolved !== ''){{
            input.value = String(resolved);
          }}
          wrap.appendChild(lbl);
          wrap.appendChild(input);
          return wrap;
        }}
        head.appendChild(handle);
        head.appendChild(cb);
        head.appendChild(txt);
        row.appendChild(head);
        if(src.kind === 'te' || src.kind === 'gene' || src.kind === 'te_diff'){{
          const saved = (typeof savedYMax === 'number') ? Math.max(1, Math.round(savedYMax)) : '';
          const fallback = (typeof defaultYMax === 'number') ? Math.max(1, Math.round(defaultYMax)) : '';
          const ymaxLabel = src.kind === 'te_diff' ? 'Y max ±' : 'Y max';
          controls.appendChild(makeYMaxInput({{ label: ymaxLabel, step: '1', min: 1, main: true, saved: saved, fallback: fallback }}));
        }} else if(src.kind === 'meth_diff'){{
          const saved = (typeof savedYMax === 'number') ? savedYMax : '';
          const fallback = (typeof defaultYMax === 'number') ? defaultYMax : '';
          controls.appendChild(makeYMaxInput({{
            label: 'Y max ±',
            step: '0.01',
            min: 0.01,
            max: 1,
            main: true,
            saved: saved,
            fallback: fallback,
          }}));
        }} else if(src.kind === 'methylation'){{
          const defaults = (defaultYMax && typeof defaultYMax === 'object') ? defaultYMax : {{}};
          const savedFrac = (savedYMax && typeof savedYMax === 'object') ? savedYMax.cg_coverage_pct : '';
          const savedLevel = (savedYMax && typeof savedYMax === 'object') ? savedYMax.meth_level : '';
          controls.appendChild(
            makeYMaxInput({{
              label: 'Y max fraction',
              step: '0.01',
              min: 0,
              max: 1,
              subKey: 'cg_coverage_pct',
              saved: savedFrac,
              fallback: defaults.cg_coverage_pct,
            }}),
          );
          controls.appendChild(
            makeYMaxInput({{
              label: 'Y max level',
              step: '0.01',
              min: 0,
              max: 1,
              subKey: 'meth_level',
              saved: savedLevel,
              fallback: defaults.meth_level,
            }}),
          );
        }}
        if(controls.childNodes.length){{
          row.appendChild(controls);
        }}
        tracksListEl.appendChild(row);
        initTrackGroupDragDrop(row);
        checkByTrackGroup.set(key, cb);
        trackGroupOrder.push(key);
      }});

      if(!m._baseTrackOrder){{
        if(Array.isArray(m.baseTrackOrder) && m.baseTrackOrder.length){{
          m._baseTrackOrder = m.baseTrackOrder.slice();
        }} else if(m.chromLayout && Array.isArray(m.chromLayout.order)){{
          m._baseTrackOrder = m.chromLayout.order.slice();
        }}
      }}
      syncChromLayoutOrderFromTrackGroups();

      const hasMethDiffTrack = m.panelTrackGroups.some((g) => String(g.key) === 'meth_diff');
      const methShowDiffRow = panel.querySelector('#te-sf-methshow-diff-row');
      if(methShowDiffRow){{
        methShowDiffRow.classList.toggle('te-sf-hidden', !hasMethDiffTrack);
      }}
    }}

    function isTrackLabelVisible(label, cl){{
      const kinds = cl && cl.trackKindsByLabel ? cl.trackKindsByLabel : {{}};
      const kind = kinds[label];
      const byGroup = cl.trackGroupKeyByLabel || {{}};
      const groupKey = byGroup[label];
      if(groupKey && !trackGroupVisible(groupKey)){{ return false; }}
      if(kind === 'methylation'){{
        const bySub = cl.methylationSubtrackByLabel || {{}};
        return methylationSubtrackTypeVisible(bySub[label] || 'cg_total');
      }}
      if(kind === 'meth_diff'){{
        return trackGroupVisible('meth_diff');
      }}
      if(kind === 'te_diff'){{
        return trackGroupVisible('te_diff');
      }}
      return true;
    }}

    function anyTrackLabelHidden(cl){{
      if(!cl || !Array.isArray(cl.order)){{ return false; }}
      return cl.order.some((label) => !isTrackLabelVisible(label, cl));
    }}

    function visibleTrackLabels(cl){{
      if(!cl || !Array.isArray(cl.order)){{ return []; }}
      return cl.order.filter((label) => isTrackLabelVisible(label, cl));
    }}

    function computeVisibleYPosByLabel(cl, rowMult){{
      const order = visibleTrackLabels(cl);
      const yByLabel = {{}};
      if(!order.length || !cl){{ return yByLabel; }}
      const mult = rowMult || 1;
      const hasDiff = order.some((lbl) => {{
        const k = (cl.trackKindsByLabel || {{}})[lbl];
        return k === 'te_diff' || k === 'meth_diff';
      }});
      if(!hasDiff){{
        const n = order.length;
        order.forEach((label, idx) => {{
          yByLabel[label] = (n - 1 - idx) * cl.rowYStep * mult;
        }});
        return yByLabel;
      }}
      let y = 0;
      const bottomUp = {{}};
      const diffBefore = Number(cl.diffBeforeStep || {DIFF_TRACK_BEFORE_STEP});
      const diffAfter = Number(cl.diffAfterGap || {DIFF_TRACK_AFTER_GAP});
      order.forEach((label, idx) => {{
        if(idx > 0){{
          const kind = (cl.trackKindsByLabel || {{}})[label];
          const prevKind = (cl.trackKindsByLabel || {{}})[order[idx - 1]];
          if(kind === 'te_diff' || kind === 'meth_diff'){{
            y += diffBefore;
          }} else if(prevKind === 'te_diff' || prevKind === 'meth_diff'){{
            y += cl.rowYStep + diffAfter;
          }} else {{
            y += cl.rowYStep;
          }}
        }}
        bottomUp[label] = y;
      }});
      const maxY = Math.max(...Object.values(bottomUp));
      order.forEach((label) => {{
        yByLabel[label] = (maxY - bottomUp[label]) * mult;
      }});
      return yByLabel;
    }}

    function trackYCenter(label, idx, cl, rowMult){{
      rowMult = rowMult || 1;
      if(!cl){{ return (idx >= 0 ? idx : 0) * 2.45 * rowMult; }}
      const yMap = computeVisibleYPosByLabel(cl, rowMult);
      if(yMap[label] !== undefined){{ return yMap[label]; }}
      return -1e6;
    }}

    function trackBarBase(label, idx, cl, rowMult){{
      return trackYCenter(label, idx, cl, rowMult) + cl.trackHalf + cl.barGapBelow;
    }}

    function methDiffTrackHalf(cl, rowMult, methDiffScale){{
      const m = gd._teScaleMeta;
      const base = m && Number(m.methDiffBarHalf) > 0
        ? Number(m.methDiffBarHalf)
        : Number(m && m.maxBarHeight || 0.36) * 0.36;
      const mult = Number.isFinite(methDiffScale) && methDiffScale > 0 ? methDiffScale : 1;
      return Math.max(0.12, base * (rowMult || 1) * mult);
    }}

    function densityAxisHeightForTrack(idx, y0, cl, rowMult, scale, maxBarHeight){{
      return maxBarHeight * scale;
    }}

    function categoriesFiltered(){{
      if(!hasCategoryFilters){{ return false; }}
      return groups.some((g) => {{
        const cb = checkByGroup.get(g);
        return !!(cb && !cb.checked);
      }});
    }}

    function panelGroupDefaultYMax(key, subKey){{
      const m = gd._teScaleMeta;
      if(!m || !Array.isArray(m.panelTrackGroups)){{ return undefined; }}
      const group = m.panelTrackGroups.find((g) => String(g.key) === String(key));
      if(!group){{ return undefined; }}
      const def = group.defaultYMax;
      if(subKey){{
        return def && typeof def === 'object' ? def[subKey] : undefined;
      }}
      return def;
    }}

    function appliedTeDiffYMax(){{
      const m = gd._teScaleMeta;
      if(!m){{ return null; }}
      const byGroup = m._appliedTrackYMaxByGroup || m._panelTrackYMaxByGroup || m.trackYMaxByGroup || {{}};
      const v = Number(byGroup.te_diff);
      if(!Number.isFinite(v) || v <= 0){{ return null; }}
      return v;
    }}

    function appliedMethDiffYMax(){{
      const m = gd._teScaleMeta;
      if(!m){{ return null; }}
      const byGroup = m._appliedTrackYMaxByGroup || m._panelTrackYMaxByGroup || m.trackYMaxByGroup || {{}};
      const v = Number(byGroup.meth_diff);
      if(!Number.isFinite(v) || v <= 0){{ return null; }}
      return v;
    }}

    function hasTeDiffYMaxOverride(){{
      const v = appliedTeDiffYMax();
      if(v === null){{ return false; }}
      const def = panelGroupDefaultYMax('te_diff');
      if(typeof def === 'number' && Math.abs(v - def) <= 0.5){{ return false; }}
      return true;
    }}

    function hasMethDiffYMaxOverride(){{
      const v = appliedMethDiffYMax();
      if(v === null){{ return false; }}
      const def = panelGroupDefaultYMax('meth_diff');
      if(typeof def === 'number' && Math.abs(v - def) <= 1e-4){{ return false; }}
      return true;
    }}

    function teDiffNormMax(){{
      const m = gd._teScaleMeta;
      if(hasTeDiffYMaxOverride()){{
        return Math.max(appliedTeDiffYMax(), 1e-9);
      }}
      return Math.max(Number(m && m.diffDensityNormMax) || 1, 1e-9);
    }}

    function methDiffNormMax(){{
      const m = gd._teScaleMeta;
      if(hasMethDiffYMaxOverride()){{
        return Math.max(appliedMethDiffYMax(), 1e-9);
      }}
      return Math.max(Number(m && m.methDiffDensityNormMax || 0.1), 1e-9);
    }}

    function isTeDiffTrackLabel(label){{
      const m = gd._teScaleMeta;
      const cl = m && m.chromLayout;
      const kinds = cl && cl.trackKindsByLabel ? cl.trackKindsByLabel : {{}};
      const resolved = resolveTrackLabel(String(label || ''));
      return kinds[resolved] === 'te_diff' || kinds[String(label || '')] === 'te_diff';
    }}

    function isMethDiffTrackLabel(label){{
      const m = gd._teScaleMeta;
      const cl = m && m.chromLayout;
      const kinds = cl && cl.trackKindsByLabel ? cl.trackKindsByLabel : {{}};
      const resolved = resolveTrackLabel(String(label || ''));
      return kinds[resolved] === 'meth_diff' || kinds[String(label || '')] === 'meth_diff';
    }}

    function trackYMaxOverrideForMeta(t, m){{
      if(!t || !m || !m.chromLayout || !m.chromLayout.trackGroupKeyByLabel){{ return null; }}
      const groupKey = String(m.chromLayout.trackGroupKeyByLabel[t.chrom] || '');
      if(!groupKey){{ return null; }}
      if(groupKey === 'te_diff'){{
        const v = appliedTeDiffYMax();
        if(v !== null && hasTeDiffYMaxOverride()){{ return v; }}
        return null;
      }}
      if(groupKey === 'meth_diff'){{
        const v = appliedMethDiffYMax();
        if(v !== null && hasMethDiffYMaxOverride()){{ return v; }}
        return null;
      }}
      const byGroup = m._appliedTrackYMaxByGroup || m._panelTrackYMaxByGroup || m.trackYMaxByGroup || m.savedTrackYMaxByGroup || {{}};
      const raw = byGroup[groupKey];
      if(raw == null){{ return null; }}
      if(t.isMethylationTrack && raw && typeof raw === 'object'){{
        const sub = t.methylationSubtrack || '';
        if(sub === 'cg_coverage_pct' || sub === 'meth_level'){{
          const v = Number(raw[sub]);
          if(Number.isFinite(v) && v > 0){{ return v; }}
        }}
        return null;
      }}
      const v = Number(raw);
      if(Number.isFinite(v) && v > 0){{ return v; }}
      return null;
    }}

    function hasTrackYMaxOverride(t){{
      const m = gd._teScaleMeta;
      if(!t || !m || !m.chromLayout || !m.chromLayout.trackGroupKeyByLabel){{ return false; }}
      const groupKey = String(m.chromLayout.trackGroupKeyByLabel[t.chrom] || '');
      if(!groupKey){{ return false; }}
      if(groupKey === 'te_diff'){{
        return hasTeDiffYMaxOverride();
      }}
      if(groupKey === 'meth_diff'){{
        return hasMethDiffYMaxOverride();
      }}
      const override = trackYMaxOverrideForMeta(t, m);
      if(override === null){{ return false; }}
      if(t.isMethylationTrack){{
        const sub = t.methylationSubtrack || '';
        const def = panelGroupDefaultYMax(groupKey, sub);
        if(Number.isFinite(def) && Math.abs(Number(override) - Number(def)) <= 1e-4){{ return false; }}
        return true;
      }}
      const def = panelGroupDefaultYMax(groupKey);
      if(typeof def === 'number' && Math.abs(Number(override) - def) <= 0.5){{ return false; }}
      return true;
    }}

    function anyTrackGroupYMaxOverride(){{
      const m = gd._teScaleMeta;
      const byGroup = (m && (m._appliedTrackYMaxByGroup || m._panelTrackYMaxByGroup || m.trackYMaxByGroup)) || {{}};
      return Object.keys(byGroup).some((key) => {{
        const raw = byGroup[key];
        if(typeof raw === 'number'){{
          const def = panelGroupDefaultYMax(key);
          if(key === 'meth_diff'){{
            if(Number.isFinite(def) && Math.abs(raw - def) <= 1e-4){{ return false; }}
          }} else if(typeof def === 'number' && Math.abs(raw - def) <= 0.5){{ return false; }}
          return Number.isFinite(raw) && raw > 0;
        }}
        if(raw && typeof raw === 'object'){{
          return Object.keys(raw).some((subKey) => {{
            const v = Number(raw[subKey]);
            const def = panelGroupDefaultYMax(key, subKey);
            if(Number.isFinite(def) && Math.abs(v - def) <= 1e-4){{ return false; }}
            return Number.isFinite(v) && v > 0;
          }});
        }}
        return false;
      }});
    }}

    function autoDensityNormForTrackMeta(t, m){{
      if(!t || !m){{ return 1; }}
      if(Number(t.densityNormMax) > 0){{ return Number(t.densityNormMax); }}
      if(t.isMethylationTrack){{
        const sub = t.methylationSubtrack || 'cg_total';
        const bySub = m.methylationDensityNormMaxBySubtrack || {{}};
        if(sub === 'cg_total' || sub === 'cg_covered'){{
          return Number(bySub.cg_count_shared || bySub.cg_total || bySub.cg_covered || 1);
        }}
        if(sub === 'cg_coverage_pct'){{
          return Number(bySub.cg_coverage_pct || 1);
        }}
        if(sub === 'meth_level'){{
          return Number(bySub.meth_level || 1);
        }}
        return Number(bySub[sub] || bySub.cg_total || 1);
      }}
      if(t.isGeneTrack){{
        return Number(m.geneDensityNormMax || 1);
      }}
      const kind = t.teInputKind || 'masking';
      if(kind === 'consensus'){{
        return Number(m.consensusDensityNormMax || m.densityNormMax || m.densityMax || 1);
      }}
      return Number(m.maskingDensityNormMax || m.densityNormMax || m.densityMax || 1);
    }}

    function densityNormForTrackMeta(t, m){{
      const override = trackYMaxOverrideForMeta(t, m);
      if(override !== null && hasTrackYMaxOverride(t)){{
        return Number(override);
      }}
      if(t.isMethylationTrack){{
        const sub = t.methylationSubtrack || 'cg_total';
        const bySub = m.methylationDensityNormMaxBySubtrack || {{}};
        if(sub === 'cg_total' || sub === 'cg_covered'){{
          return Number(bySub.cg_count_shared || bySub.cg_total || bySub.cg_covered || 1);
        }}
        if(sub === 'cg_coverage_pct'){{
          return Number(bySub.cg_coverage_pct || 1);
        }}
        if(sub === 'meth_level'){{
          return Number(bySub.meth_level || 1);
        }}
        return Number(bySub[sub] || bySub.cg_total || 1);
      }}
      if(t.isGeneTrack){{
        return Number(m.geneDensityNormMax || 1);
      }}
      const kind = t.teInputKind || 'masking';
      if(kind === 'consensus'){{
        return Number(m.consensusDensityNormMax || m.densityNormMax || m.densityMax || 1);
      }}
      return Number(m.maskingDensityNormMax || m.densityNormMax || m.densityMax || 1);
    }}

    const originalRawY = new Map();

    function resolveTrackLabel(label){{
      const s = String(label ?? '').trim();
      const cl = gd._teScaleMeta && gd._teScaleMeta.chromLayout;
      if(!cl || !Array.isArray(cl.order)){{ return s; }}
      const idx = cl.order.indexOf(s);
      if(idx >= 0){{ return cl.order[idx]; }}
      return s;
    }}

    function trackMetaForLabel(label){{
      const m = gd._teScaleMeta;
      if(!m || !Array.isArray(m.tracks)){{ return null; }}
      const s = resolveTrackLabel(label);
      return m.tracks.find((t) => t.chrom === s || String(t.chrom) === s) || null;
    }}

    function rawDensityFromRendered(yVal, t){{
      const m = gd._teScaleMeta;
      if(!m || !t){{ return Number(yVal || 0); }}
      const normMax = autoDensityNormForTrackMeta(t, m);
      const axisH = Number(t.axisHeight) > 0 ? Number(t.axisHeight) : Number(m.maxBarHeight || 1);
      if(axisH <= 0 || normMax <= 0){{ return Number(yVal || 0); }}
      return Number(yVal || 0) * normMax / axisH;
    }}

    function populateOriginalRawY(){{
      originalRawY.clear();
      const m = gd._teScaleMeta;
      function rawTeFromRendered(yVal, chrom){{
        if(!m || !m.maxBarHeight){{ return Number(yVal || 0); }}
        const t = trackMetaForLabel(chrom);
        const normMax = t ? autoDensityNormForTrackMeta(t, m) : Number(m.densityNormMax || 1);
        if(!normMax){{ return Number(yVal || 0); }}
        return Number(yVal || 0) * normMax / m.maxBarHeight;
      }}
      function rawDiffFromRendered(yVal, tr, i){{
        const cd = tr && Array.isArray(tr.customdata) ? tr.customdata : null;
        const row = cd && cd[i];
        if(row && row.length > 5){{
          const v = Number(row[5]);
          if(Number.isFinite(v)){{ return Math.abs(v); }}
        }}
        const initialDiffMax = Math.max(Number(m.diffDensityNormMax || 1), 1e-9);
        const diffHalf0 = Math.max(0.12, Number(m.diffBarHalf || m.maxBarHeight * 0.36));
        return Math.abs(Number(yVal || 0)) * initialDiffMax / diffHalf0;
      }}
      barTraceIndices.forEach((idx) => {{
        const tr = gd.data[idx];
        if(!tr){{ return; }}
        const g = tr.legendgroup || '';
        if(g === '__te_diff__' || g === '__meth_diff__'){{
          const ys0 = originalY.get(idx) || [];
          const raw = new Array(ys0.length);
          for(let i = 0; i < ys0.length; i += 1){{
            raw[i] = rawDiffFromRendered(ys0[i], tr, i);
          }}
          if(raw.length){{ originalRawY.set(idx, raw); }}
          return;
        }}
        const keys = pointKeys.get(idx) || [];
        const ys0 = originalY.get(idx) || [];
        const cd = Array.isArray(tr.customdata) ? tr.customdata : [];
        const raw = new Array(ys0.length);
        for(let i = 0; i < ys0.length; i += 1){{
          const row = cd[i] || [];
          if(row.length > 5){{
            const v = Number(row[5]);
            if(Number.isFinite(v)){{
              raw[i] = v;
              continue;
            }}
          }}
          const trackLabel = String(keys[i] || '').split('|')[0] || '';
          const t = trackMetaForLabel(trackLabel);
          if(g && barsByGroup.has(g)){{
            raw[i] = rawTeFromRendered(ys0[i], trackLabel);
          }} else if(t && (t.isGeneTrack || t.isMethylationTrack)){{
            raw[i] = rawDensityFromRendered(ys0[i], t);
          }} else {{
            raw[i] = Number(ys0[i] || 0);
          }}
        }}
        if(raw.length){{ originalRawY.set(idx, raw); }}
      }});
    }}

    function categoryStackTraceIndices(){{
      const out = [];
      barsByGroup.forEach((arr) => {{
        arr.forEach((idx) => {{ out.push(idx); }});
      }});
      return out;
    }}

    function activeTickMaxForTrack(t, m, filtered){{
      const kind = t.isMethylationTrack
        ? 'methylation'
        : (t.isGeneTrack ? 'gene' : (t.teInputKind || 'masking'));
      const groupNorm = densityNormForTrackMeta(t, m);
      if(!filtered){{
        return groupNorm;
      }}
      const byKind = m._activeTickMaxByKind;
      if(byKind && Number(byKind[kind]) > 0){{
        return Number(byKind[kind]);
      }}
      return groupNorm;
    }}

    function updateWindowRefLineShapes(layoutUpdate, t, y0, axisHeight, tickMax, showRefLine, refColor){{
      if(t.isGeneTrack || t.isDiffTrack || t.isMethylationTrack){{ return; }}
      if(t.windowRefLineShapeIdx === undefined){{ return; }}
      const refVal = Number(gd._teScaleMeta.windowRefDensity);
      const showRef = !!showRefLine && refVal > 0 && refVal <= tickMax * 1.001;
      if(!showRef){{
        layoutUpdate['shapes[' + t.windowRefLineShapeIdx + '].visible'] = false;
        if(t.windowRefTickShapeIdx !== undefined){{
          layoutUpdate['shapes[' + t.windowRefTickShapeIdx + '].visible'] = false;
        }}
        return;
      }}
      const refFrac = Math.min(1, refVal / tickMax);
      const yRef = y0 + axisHeight * refFrac;
      layoutUpdate['shapes[' + t.windowRefLineShapeIdx + '].y0'] = yRef;
      layoutUpdate['shapes[' + t.windowRefLineShapeIdx + '].y1'] = yRef;
      layoutUpdate['shapes[' + t.windowRefLineShapeIdx + '].line.color'] = refColor;
      layoutUpdate['shapes[' + t.windowRefLineShapeIdx + '].visible'] = true;
      if(t.windowRefTickShapeIdx !== undefined){{
        layoutUpdate['shapes[' + t.windowRefTickShapeIdx + '].y0'] = yRef;
        layoutUpdate['shapes[' + t.windowRefTickShapeIdx + '].y1'] = yRef;
        layoutUpdate['shapes[' + t.windowRefTickShapeIdx + '].line.color'] = refColor;
        layoutUpdate['shapes[' + t.windowRefTickShapeIdx + '].visible'] = true;
      }}
    }}

    function updateTrackDensityAxes(scale, layoutUpdate, visibleMax, rowMult, diffScale, methDiffScale){{
      const m = gd._teScaleMeta;
      if(!m || !Array.isArray(m.tracks)){{ return; }}
      const diffMult = Number.isFinite(diffScale) && diffScale > 0 ? diffScale : 1;
      const methDiffMult = Number.isFinite(methDiffScale) && methDiffScale > 0 ? methDiffScale : 1;
      const cl = m.chromLayout;
      m.tracks.forEach((t) => {{
        const trackVisible = !cl || isTrackLabelVisible(t.chrom, cl);
        if(t.isDiffTrack){{
          const kind = cl && cl.trackKindsByLabel ? cl.trackKindsByLabel[t.chrom] : '';
          const isMethDiff = kind === 'meth_diff' || !!t.isMethDiffTrack;
          let yCenter = Number(t.yCenter || 0);
          let idx = -1;
          if(cl && Array.isArray(cl.order)){{
            idx = cl.order.indexOf(t.chrom);
            if(idx >= 0){{
              yCenter = trackYCenter(t.chrom, idx, cl, rowMult);
            }}
          }}
          const halfH = isMethDiff
            ? methDiffTrackHalf(cl, rowMult, methDiffMult)
            : Math.max(
              0.12,
              Number(m.diffBarHalf || m.maxBarHeight * 0.36) * rowMult * diffMult,
            );
          const y0 = yCenter - halfH;
          const y1 = yCenter + halfH;
          const diffMax = isMethDiff
            ? methDiffNormMax()
            : teDiffNormMax();
          layoutUpdate['shapes[' + t.shapeIdx + '].visible'] = trackVisible;
          layoutUpdate['shapes[' + t.shapeIdx + '].y0'] = y0;
          layoutUpdate['shapes[' + t.shapeIdx + '].y1'] = y1;
          (t.ticks || []).forEach((tick) => {{
            layoutUpdate['shapes[' + tick.shapeIdx + '].visible'] = trackVisible;
            layoutUpdate['annotations[' + tick.annIdx + '].visible'] = trackVisible;
            if(!trackVisible){{ return; }}
            const frac = Number(tick.frac || 0);
            const yTick = y0 + (y1 - y0) * frac;
            layoutUpdate['shapes[' + tick.shapeIdx + '].y0'] = yTick;
            layoutUpdate['shapes[' + tick.shapeIdx + '].y1'] = yTick;
            layoutUpdate['annotations[' + tick.annIdx + '].y'] = yTick;
            layoutUpdate['annotations[' + tick.annIdx + '].text'] = isMethDiff
              ? formatMethDiffAxisTick(frac, diffMax)
              : formatDiffAxisTick(frac, diffMax, !!m.teLengthMode);
          }});
          return;
        }}
        let y0 = Number(t.y0 || 0);
        let idx = -1;
        if(cl && Array.isArray(cl.order)){{
          idx = cl.order.indexOf(t.chrom);
          if(idx >= 0){{
            y0 = trackBarBase(t.chrom, idx, cl, rowMult);
          }}
        }}
        const axisHeight = Number(t.axisHeight) > 0
          ? Math.min(
              Number(t.axisHeight) * rowMult * scale,
              densityAxisHeightForTrack(idx, y0, cl, rowMult, scale, m.maxBarHeight),
            )
          : densityAxisHeightForTrack(idx, y0, cl, rowMult, scale, m.maxBarHeight);
        const nominalAxisH = Math.max(1e-12, m.maxBarHeight * scale);
        const filtered = categoriesFiltered();
        const useSharedTe = (
          !t.isGeneTrack && !t.isMethylationTrack && !t.isDiffTrack
          && m.sharedTeDensityAxis !== false
          && !hasTrackYMaxOverride(t)
        );
        const tickMax = useSharedTe
          ? activeTickMaxForTrack(t, m, filtered)
          : (
            hasTrackYMaxOverride(t)
              ? densityNormForTrackMeta(t, m)
              : Number(t.densityNormMax || densityNormForTrackMeta(t, m))
                * (axisHeight / nominalAxisH)
          );
        layoutUpdate['shapes[' + t.shapeIdx + '].visible'] = trackVisible;
        layoutUpdate['shapes[' + t.shapeIdx + '].y0'] = y0;
        layoutUpdate['shapes[' + t.shapeIdx + '].y1'] = y0 + axisHeight;
        (t.ticks || []).forEach((tick) => {{
          layoutUpdate['shapes[' + tick.shapeIdx + '].visible'] = trackVisible;
          layoutUpdate['annotations[' + tick.annIdx + '].visible'] = trackVisible;
          if(!trackVisible){{ return; }}
          const frac = Number(tick.frac || 0);
          const yTick = y0 + axisHeight * frac;
          layoutUpdate['shapes[' + tick.shapeIdx + '].y0'] = yTick;
          layoutUpdate['shapes[' + tick.shapeIdx + '].y1'] = yTick;
          layoutUpdate['annotations[' + tick.annIdx + '].y'] = yTick;
          let tickText = formatDensityTick(tickMax * frac, !!m.teLengthMode);
          if(t.isGeneTrack){{
            tickText = formatDensityTick(tickMax * frac, !!m.annotationLengthMode);
          }} else if(t.isMethylationTrack && t.methylationSubtrack === 'cg_coverage_pct'){{
            tickText = formatUnitIntervalAxisTick(tickMax * frac);
          }} else if(t.isMethylationTrack && t.methylationSubtrack === 'meth_level'){{
            tickText = formatUnitIntervalAxisTick(tickMax * frac);
          }} else if(t.isMethylationTrack){{
            tickText = formatDensityTick(tickMax * frac, false);
          }}
          layoutUpdate['annotations[' + tick.annIdx + '].text'] = tickText;
        }});
        const showRefLine = trackVisible && !!m.showWindowRefLine && Number(m.windowRefDensity) > 0;
        const refColor = String(m.windowRefColor || '{DEFAULT_WINDOW_REF_COLOR}');
        updateWindowRefLineShapes(layoutUpdate, t, y0, axisHeight, tickMax, showRefLine, refColor);
      }});
    }}

    function invalidateLayoutCache(){{
      gd._teCachedYMin = undefined;
      gd._teCachedAxisTop = undefined;
      gd._teCachedMaxStackTop = undefined;
      gd._teYRangeLocked = null;
      gd._teCategoryFastRelayout = false;
    }}

    function updateScatterTrackY(rowMult){{
      const m = gd._teScaleMeta;
      const cl = m && m.chromLayout;
      if(!cl || !Array.isArray(cl.scatterIdx)){{ return; }}
      cl.order.forEach((chrom, i) => {{
        if(!isTrackLabelVisible(chrom, cl)){{ return; }}
        const ti = cl.scatterIdx[i];
        if(ti === undefined || !gd.data[ti]){{ return; }}
        const yc = trackYCenter(chrom, i, cl, rowMult);
        const n = traceArrayLength(gd.data[ti].y);
        gd.data[ti].y = new Array(n).fill(yc);
      }});
    }}

    let lastVisibleMask = groups.map(() => true).join('|');
    function checkedGroupsInOrder(){{
      return selectedOrder.filter((g) => {{
        const cb = checkByGroup.get(g);
        return !!(cb && cb.checked);
      }});
    }}

    function stackGroupsInOrder(){{
      const ordered = checkedGroupsInOrder();
      if(hasCategoryFilters){{
        return ordered.slice().reverse();
      }}
      return ordered;
    }}

    function writeBarsToDom(yByIdx, baseByIdx){{
      barTraceIndices.forEach((idx) => {{
        if(yByIdx.has(idx)){{ gd.data[idx].y = yByIdx.get(idx); }}
        if(baseByIdx.has(idx)){{ gd.data[idx].base = baseByIdx.get(idx); }}
      }});
    }}

    function applyBarDataViaRestyle(yByIdx, baseByIdx, opts){{
      opts = opts || {{}};
      const traceSubset = opts.traceSubset;
      const visibleByIdx = opts.visibleByIdx;
      const indices = traceSubset || barTraceIndices;
      if(!indices.length){{ return Promise.resolve(); }}
      const ys = [];
      const bases = [];
      const vis = visibleByIdx ? [] : null;
      indices.forEach((idx) => {{
        ys.push(yByIdx.has(idx) ? yByIdx.get(idx) : (originalY.get(idx) || []));
        bases.push(baseByIdx.has(idx) ? baseByIdx.get(idx) : (originalBase.get(idx) || []));
        if(vis){{
          if(plainBarIndices.includes(idx)){{
            vis.push(true);
          }} else {{
            vis.push(visibleByIdx.has(idx) ? visibleByIdx.get(idx) : !!gd.data[idx].visible);
          }}
        }}
      }});
      const patch = {{y: ys, base: bases}};
      if(vis){{ patch.visible = vis; }}
      return Plotly.restyle(gd, patch, indices);
    }}

    function buildTraceVisibilityMap(){{
      const byIdx = new Map();
      gd.data.forEach((tr, idx) => {{
        if(plainBarIndices.includes(idx)){{ return; }}
        const g = tr.legendgroup || tr.name || '';
        if(!g){{ return; }}
        const cb = checkByGroup.get(g);
        if(!cb){{ return; }}
        byIdx.set(idx, !!cb.checked);
      }});
      return byIdx;
    }}

    function applyScatterViaRestyle(rowMult){{
      const m = gd._teScaleMeta;
      const cl = m && m.chromLayout;
      if(!cl || !Array.isArray(cl.scatterIdx)){{ return Promise.resolve(); }}
      const indices = [];
      const ys = [];
      cl.order.forEach((chrom, i) => {{
        if(!isTrackLabelVisible(chrom, cl)){{ return; }}
        const ti = cl.scatterIdx[i];
        if(ti === undefined || !gd.data[ti]){{ return; }}
        const yc = trackYCenter(chrom, i, cl, rowMult);
        const n = traceArrayLength(gd.data[ti].y);
        indices.push(ti);
        ys.push(new Array(n).fill(yc));
      }});
      if(!indices.length){{ return Promise.resolve(); }}
      return Plotly.restyle(gd, {{y: ys}}, indices);
    }}

    function applyScale(){{
      loadScaleMeta(gd);
      if(GENOME_MODE && !gd._teScaleMeta){{ return Promise.resolve(); }}
      if(GENOME_MODE && gd._teScaleMeta && trackGroupOrder.length){{
        syncChromLayoutOrderFromTrackGroups();
      }}
      rescanBarTracesIfNeeded();
      populateOriginalRawY();
      if(GENOME_MODE && gd._teScaleMeta){{
        const refCb = panel.querySelector('#te-sf-window-ref');
        if(refCb){{
          gd._teScaleMeta.showWindowRefLine = !!refCb.checked;
        }}
        const refColorEl = panel.querySelector('#te-sf-window-ref-color');
        if(refColorEl){{
          gd._teScaleMeta.windowRefColor = refColorEl.value || '{DEFAULT_WINDOW_REF_COLOR}';
        }}
      }}
      const scaleInput = panel.querySelector('#te-sf-scale');
      let scale = Number(scaleInput ? scaleInput.value : 1);
      if(!Number.isFinite(scale) || scale <= 0){{ scale = 1; }}
      scale = Math.min(SCALE_MAX, Math.max({PLOT_SCALE_MIN}, scale));

      let categoryOnly = false;
      let methylationOnly = false;
      {apply_body}

      const visibleByIdx = hasCategoryFilters ? buildTraceVisibilityMap() : null;
      let visibleMask = '';
      if(hasCategoryFilters){{
        const parts = [];
        gd.data.forEach((tr) => {{
          const g = tr.legendgroup || tr.name || '';
          if(!g){{ return; }}
          const cb = checkByGroup.get(g);
          if(!cb){{ return; }}
          parts.push(cb.checked ? '1' : '0');
        }});
        visibleMask = parts.join('|');
      }}

      function applyDivergenceBars(yByIdx, baseByIdx){{
        writeBarsToDom(yByIdx, baseByIdx || new Map());
        const indices = barTraceIndices.slice();
        const ys = indices.map((idx) => yByIdx.get(idx) || originalY.get(idx) || []);
        const bases = indices.map((idx) => {{
          if(baseByIdx && baseByIdx.has(idx)){{ return baseByIdx.get(idx); }}
          return originalBase.get(idx) || [];
        }});
        if(!indices.length){{ return Promise.resolve(); }}
        return Plotly.restyle(gd, {{y: ys, base: bases}}, indices);
      }}

      function applyDivergenceLayout(layoutUpdate){{
        gd._teRelayoutSync = true;
        return Plotly.relayout(gd, layoutUpdate).finally(() => {{
          gd._teRelayoutSync = false;
        }});
      }}

      function redrawAndRelayout(){{
        if(GENOME_MODE && gd._teScaleMeta){{
          preserveXRange(layoutUpdate);
          const catStackIdx = categoryStackTraceIndices();
          const restyleOpts = {{
            traceSubset: categoryOnly ? catStackIdx : null,
            visibleByIdx: visibleByIdx,
          }};
          const skipColors = categoryOnly || methylationOnly;
          const skipFullAxis = categoryOnly && !!gd._teCategoryFastRelayout;
          if(!skipFullAxis){{
            updateTrackDensityAxes(
              scale,
              layoutUpdate,
              gd._teVisibleDensityMax,
              gd._teRowMult || 1,
              gd._teDiffScale || 1,
              gd._teMethDiffScale || 1,
            );
          }}
          if(categoryOnly){{ gd._teCategoryFastRelayout = false; }}
          writeBarsToDom(yByIdx, baseByIdx);
          let barPromise;
          if(categoryOnly){{
            if(visibleByIdx){{
              catStackIdx.forEach((idx) => {{
                if(visibleByIdx.has(idx)){{
                  gd.data[idx].visible = visibleByIdx.get(idx);
                }}
              }});
            }}
            const plainSubset = plainBarIndices.filter((idx) => catStackIdx.indexOf(idx) < 0);
            barPromise = applyBarDataViaRestyle(yByIdx, baseByIdx, {{
              traceSubset: catStackIdx.concat(plainSubset),
              visibleByIdx: visibleByIdx,
            }});
          }} else {{
            barPromise = applyBarDataViaRestyle(yByIdx, baseByIdx, restyleOpts);
          }}
          const scatterPromise = categoryOnly
            ? Promise.resolve()
            : applyScatterViaRestyle(gd._teRowMult || 1);
          return barPromise
            .then(() => scatterPromise)
            .then(() => {{
              if(hasCategoryFilters){{
                lastVisibleMask = visibleMask;
              }}
              if(!Object.keys(layoutUpdate).length){{ return Promise.resolve(); }}
              return safeRelayout(layoutUpdate);
            }})
            .then(() => (skipColors ? Promise.resolve() : applyCategoryColors()))
            .then(() => (skipColors ? Promise.resolve() : applyDiffTrackColors()))
            .then(() => (skipColors ? Promise.resolve() : applyMethDiffTrackColors()))
            .then(() => (skipColors ? Promise.resolve() : applyMethylationTrackColors()));
        }}
        return applyDivergenceBars(yByIdx, baseByIdx)
          .then(() => {{
            applyXAxisZoomLayout(layoutUpdate);
            return applyDivergenceLayout(layoutUpdate);
          }})
          .then(() => {{
            if(typeof Plotly.redraw === 'function'){{
              return Plotly.redraw(gd);
            }}
          }})
          .then(() => applyCategoryColors())
          .then(() => applyDiffTrackColors())
          .then(() => applyMethDiffTrackColors())
          .then(() => applyMethylationTrackColors());
      }}

      gd._teApplyDepth = (gd._teApplyDepth || 0) + 1;
      const applyPromise = redrawAndRelayout();
      return applyPromise.finally(() => {{
        gd._teApplyDepth = Math.max(0, (gd._teApplyDepth || 1) - 1);
        if(GENOME_MODE && gd._teScaleMeta){{
          gd._teScaleMeta._yMaxRescaleActive = false;
        }}
        if(gd._teApplyDepth === 0 && !(categoryOnly || methylationOnly)){{
          return refreshPlotLayout();
        }}
      }});
    }}

    let applyRaf = 0;
    let delayedApply = 0;
    function scheduleApply(){{
      if(applyRaf){{ return; }}
      applyRaf = window.requestAnimationFrame(() => {{
        applyRaf = 0;
        applyScale();
      }});
    }}
    function scheduleApplyDebounced(ms){{
      if(delayedApply){{ window.clearTimeout(delayedApply); }}
      delayedApply = window.setTimeout(() => {{
        delayedApply = 0;
        scheduleApply();
      }}, ms);
    }}

    panel.addEventListener('change', (e) => {{
      const t = e.target;
      if(t && t.matches('input[type="checkbox"]')){{
        if(t.id === 'te-sf-te-master'){{
          const on = !!t.checked;
          checkByGroup.forEach((cb) => {{ cb.checked = on; }});
          syncSelectedOrderFromList();
          gd._teApplyFastPath = '';
          scheduleApplyDebounced(80);
          schedulePrefsSave();
          return;
        }}
        const g = t.dataset.group || '';
        if(g){{
          syncSelectedOrderFromList();
        }}
        if(t.id && String(t.id).startsWith('te-sf-methshow')){{
          gd._teApplyFastPath = 'methylation';
        }} else if(t.id && String(t.id).startsWith('te-sf-track-')){{
          gd._teApplyFastPath = '';
          invalidateLayoutCache();
          syncChromLayoutOrderFromTrackGroups();
        }} else if(t.dataset && t.dataset.group){{
          gd._teApplyFastPath = '';
          syncTeMasterCheckbox();
        }} else {{
          gd._teApplyFastPath = '';
        }}
        scheduleApplyDebounced(80);
        schedulePrefsSave();
      }}
      if(t && t.matches('input[type="color"]')){{
        scheduleApply();
        schedulePrefsSave();
      }}
    }});

    function applySavedPanelState(){{
      const m = GENOME_MODE ? gd._teScaleMeta : gd._teDivergenceMeta;
      if(!m){{ return; }}
      const scaleInput = panel.querySelector('#te-sf-scale');
      if(scaleInput && Number.isFinite(m.savedBarScale) && m.savedBarScale > 0){{
        scaleInput.value = String(m.savedBarScale);
        panel.querySelector('#te-sf-scale-val').textContent = Number(m.savedBarScale).toFixed(1) + 'x';
      }}
      const rowInput = panel.querySelector('#te-sf-rowsp');
      if(rowInput && Number.isFinite(m.savedRowSpacing) && m.savedRowSpacing > 0){{
        rowInput.value = String(m.savedRowSpacing);
        panel.querySelector('#te-sf-rowsp-val').textContent = Number(m.savedRowSpacing).toFixed(1) + 'x';
      }}
      if(m.chromLayout){{
        const savedSort = m.savedChromosomeSort || m.chromLayout.chromosomeSort || 'natural';
        if(savedSort !== (m.chromLayout.chromosomeSort || 'natural')){{
          applyChromosomeSort(savedSort);
        }}
      }}
      const diffWrap = panel.querySelector('#te-sf-diffscale-wrap');
      const diffColorsWrap = panel.querySelector('#te-sf-diffcolors-wrap');
      const diffInput = panel.querySelector('#te-sf-diffscale');
      if(m.hasDiffTrack){{
        if(diffWrap){{ diffWrap.classList.remove('te-sf-hidden'); }}
        if(diffColorsWrap){{ diffColorsWrap.classList.remove('te-sf-hidden'); }}
        const lblA = panel.querySelector('#te-sf-diffcolor-a-lbl');
        const lblB = panel.querySelector('#te-sf-diffcolor-b-lbl');
        if(lblA){{ lblA.textContent = 'TE difference track positive color'; }}
        if(lblB){{ lblB.textContent = 'TE difference track negative color'; }}
        const diffColorA = panel.querySelector('#te-sf-diffcolor-a');
        const diffColorB = panel.querySelector('#te-sf-diffcolor-b');
        if(diffColorA){{
          diffColorA.value = normalizeHexColor(m.savedDiffColorA || m.diffColorA || '{TE_DIFF_COLOR_A}');
        }}
        if(diffColorB){{
          diffColorB.value = normalizeHexColor(m.savedDiffColorB || m.diffColorB || '{TE_DIFF_COLOR_B}');
        }}
      }}
      const methColorsWrap = panel.querySelector('#te-sf-methcolors-wrap');
      if(m.hasMethylationTracks){{
        if(methColorsWrap){{ methColorsWrap.classList.remove('te-sf-hidden'); }}
        const methCgTotal = panel.querySelector('#te-sf-methcolor-cgtotal');
        const methCgCovered = panel.querySelector('#te-sf-methcolor-cgcovered');
        const methCgPct = panel.querySelector('#te-sf-methcolor-cgpct');
        const methLevel = panel.querySelector('#te-sf-methcolor-level');
        if(methCgTotal){{
          methCgTotal.value = normalizeHexColor(m.savedMethColorCgTotal || m.methColorCgTotal || '{METH_COLOR_CG_TOTAL}');
        }}
        if(methCgCovered){{
          methCgCovered.value = normalizeHexColor(m.savedMethColorCgCovered || m.methColorCgCovered || '{METH_COLOR_CG_COVERED}');
        }}
        if(methCgPct){{
          methCgPct.value = normalizeHexColor(m.savedMethColorCgPct || m.methColorCgPct || '{METH_COLOR_CG_COVERAGE_PCT}');
        }}
        if(methLevel){{
          methLevel.value = normalizeHexColor(m.savedMethColorLevel || m.methColorLevel || '{METH_COLOR_METH_LEVEL}');
        }}
        const methShowDiffRow = panel.querySelector('#te-sf-methshow-diff-row');
        if(methShowDiffRow){{
          methShowDiffRow.classList.toggle('te-sf-hidden', !m.hasMethDiffTrack);
        }}
        if(m.hasMethDiffTrack){{
          const methDiffWrap = panel.querySelector('#te-sf-methdiffscale-wrap');
          const methDiffColorsWrap = panel.querySelector('#te-sf-methdiffcolors-wrap');
          if(methDiffWrap){{ methDiffWrap.classList.remove('te-sf-hidden'); }}
          if(methDiffColorsWrap){{ methDiffColorsWrap.classList.remove('te-sf-hidden'); }}
          const lblA = panel.querySelector('#te-sf-methdiffcolor-a-lbl');
          const lblB = panel.querySelector('#te-sf-methdiffcolor-b-lbl');
          if(lblA){{ lblA.textContent = 'Methylation difference track positive color'; }}
          if(lblB){{ lblB.textContent = 'Methylation difference track negative color'; }}
        }}
      }}
      if(diffInput && Number.isFinite(m.savedDiffTrackScale) && m.savedDiffTrackScale > 0){{
        diffInput.value = String(m.savedDiffTrackScale);
        panel.querySelector('#te-sf-diffscale-val').textContent =
          Number(m.savedDiffTrackScale).toFixed(1) + 'x';
      }}
      const methDiffInput = panel.querySelector('#te-sf-methdiffscale');
      if(methDiffInput && Number.isFinite(m.savedMethDiffTrackScale) && m.savedMethDiffTrackScale > 0){{
        methDiffInput.value = String(m.savedMethDiffTrackScale);
        const methDiffVal = panel.querySelector('#te-sf-methdiffscale-val');
        if(methDiffVal){{
          methDiffVal.textContent = Number(m.savedMethDiffTrackScale).toFixed(1) + 'x';
        }}
      }}
      if(m.hasMethDiffTrack){{
        const methDiffColorA = panel.querySelector('#te-sf-methdiffcolor-a');
        const methDiffColorB = panel.querySelector('#te-sf-methdiffcolor-b');
        if(methDiffColorA){{
          methDiffColorA.value = normalizeHexColor(
            m.savedMethDiffColorA || m.methDiffColorA || '{METH_DIFF_COLOR_A}',
          );
        }}
        if(methDiffColorB){{
          methDiffColorB.value = normalizeHexColor(
            m.savedMethDiffColorB || m.methDiffColorB || '{METH_DIFF_COLOR_B}',
          );
        }}
      }}
      if(m.hasMethylationTracks){{
        const showMap = [
          ['#te-sf-methshow-cgtotal', m.savedMethShowCgTotal],
          ['#te-sf-methshow-cgcovered', m.savedMethShowCgCovered],
          ['#te-sf-methshow-cgpct', m.savedMethShowCgPct],
          ['#te-sf-methshow-level', m.savedMethShowLevel],
          ['#te-sf-methshow-diff', m.savedMethShowDiff],
        ];
        showMap.forEach(([sel, saved]) => {{
          if(saved === undefined){{ return; }}
          const cb = panel.querySelector(sel);
          if(cb){{ cb.checked = !!saved; }}
        }});
      }}
      const refInput = panel.querySelector('#te-sf-window-ref');
      if(refInput && m.savedShowWindowRefLine !== undefined){{
        refInput.checked = !!m.savedShowWindowRefLine;
      }}
      const refColorInput = panel.querySelector('#te-sf-window-ref-color');
      if(refColorInput){{
        const savedRefColor = m.savedWindowRefColor || m.windowRefColor;
        if(savedRefColor){{ refColorInput.value = normalizeHexColor(savedRefColor); }}
      }}
      if(Array.isArray(m.savedCategoryHidden) && m.savedCategoryHidden.length){{
        const hidden = new Set(m.savedCategoryHidden);
        checkByGroup.forEach((cb, g) => {{
          if(hidden.has(g)){{ cb.checked = false; }}
        }});
      }}
      if(m.savedCategoryColors && typeof m.savedCategoryColors === 'object'){{
        Object.entries(m.savedCategoryColors).forEach(([g, c]) => {{
          const colorInput = colorByGroup.get(g);
          if(colorInput && c){{ colorInput.value = normalizeHexColor(c); }}
        }});
      }}
      if(Array.isArray(m.savedCategoryOrder) && m.savedCategoryOrder.length){{
        reorderListRowsToGroups(m.savedCategoryOrder.filter((g) => checkByGroup.has(g)));
      }}
      buildTracksPanel();
      if(m.savedTrackYMaxByGroup && typeof m.savedTrackYMaxByGroup === 'object'
        && Object.keys(m.savedTrackYMaxByGroup).length){{
        m._appliedTrackYMaxByGroup = JSON.parse(JSON.stringify(m.savedTrackYMaxByGroup));
        m._panelTrackYMaxByGroup = m._appliedTrackYMaxByGroup;
        m.trackYMaxByGroup = m._appliedTrackYMaxByGroup;
      }}
      syncTeMasterCheckbox();
    }}

    const scaleEl = panel.querySelector('#te-sf-scale');
    scaleEl.max = String(SCALE_MAX);
    scaleEl.addEventListener('input', (e) => {{
      const v = Number(e.target.value || 1);
      panel.querySelector('#te-sf-scale-val').textContent = v.toFixed(1) + 'x';
      scheduleApplyDebounced(60);
      schedulePrefsSave();
    }});
    scaleEl.addEventListener('change', () => {{
      scheduleApply();
      schedulePrefsSave();
    }});
    const xZoomEl = panel.querySelector('#te-sf-xzoom');
    if(xZoomEl){{
      xZoomEl.max = String(XZOOM_MAX);
      function applyXZoomRelayout(){{
        const layoutUpdate = {{}};
        applyXAxisZoomLayout(layoutUpdate);
        if(!Object.keys(layoutUpdate).length){{ return Promise.resolve(); }}
        return safeRelayout(layoutUpdate);
      }}
      xZoomEl.addEventListener('input', (e) => {{
        const v = Number(e.target.value || 1);
        const valEl = panel.querySelector('#te-sf-xzoom-val');
        if(valEl){{ valEl.textContent = v.toFixed(1) + 'x'; }}
        applyXZoomRelayout();
        schedulePrefsSave();
      }});
      xZoomEl.addEventListener('change', () => {{
        applyXZoomRelayout();
        schedulePrefsSave();
      }});
    }}
    if(GENOME_MODE){{
      const rowEl = panel.querySelector('#te-sf-rowsp');
      if(rowEl){{
        rowEl.addEventListener('input', (e) => {{
          const v = Number(e.target.value || 1).toFixed(1);
          panel.querySelector('#te-sf-rowsp-val').textContent = v + 'x';
          scheduleApplyDebounced(60);
          schedulePrefsSave();
        }});
        rowEl.addEventListener('change', () => {{
          scheduleApply();
          schedulePrefsSave();
        }});
      }}
      const diffEl = panel.querySelector('#te-sf-diffscale');
      if(diffEl){{
        diffEl.addEventListener('input', (e) => {{
          const v = Number(e.target.value || 1).toFixed(1);
          panel.querySelector('#te-sf-diffscale-val').textContent = v + 'x';
          scheduleApplyDebounced(60);
          schedulePrefsSave();
        }});
        diffEl.addEventListener('change', () => {{
          scheduleApply();
          schedulePrefsSave();
        }});
      }}
      const methDiffEl = panel.querySelector('#te-sf-methdiffscale');
      if(methDiffEl){{
        methDiffEl.addEventListener('input', (e) => {{
          const v = Number(e.target.value || 1).toFixed(1);
          const valEl = panel.querySelector('#te-sf-methdiffscale-val');
          if(valEl){{ valEl.textContent = v + 'x'; }}
          scheduleApplyDebounced(60);
          schedulePrefsSave();
        }});
        methDiffEl.addEventListener('change', () => {{
          scheduleApply();
          schedulePrefsSave();
        }});
      }}
      const refCb = panel.querySelector('#te-sf-window-ref');
      if(refCb){{
        if(gd._teScaleMeta && gd._teScaleMeta.savedShowWindowRefLine === undefined){{
          refCb.checked = !!gd._teScaleMeta.showWindowRefLine;
        }}
        refCb.addEventListener('change', () => {{
          scheduleApply();
          schedulePrefsSave();
        }});
      }}
      const refColorCb = panel.querySelector('#te-sf-window-ref-color');
      if(refColorCb){{
        refColorCb.addEventListener('input', () => {{
          scheduleApply();
          schedulePrefsSave();
        }});
        refColorCb.addEventListener('change', () => {{
          scheduleApply();
          schedulePrefsSave();
        }});
      }}
      ['#te-sf-methshow-cgtotal', '#te-sf-methshow-cgcovered', '#te-sf-methshow-cgpct', '#te-sf-methshow-level', '#te-sf-methshow-diff'].forEach((sel) => {{
        const el = panel.querySelector(sel);
        if(el){{
          el.addEventListener('change', () => {{
            if(sel === '#te-sf-methshow-diff'){{
              const trackCb = panel.querySelector('#te-sf-track-meth_diff');
              if(trackCb){{ trackCb.checked = !!el.checked; }}
            }}
            scheduleApply();
            schedulePrefsSave();
          }});
        }}
      }});
      panel.addEventListener('change', (e) => {{
        const t = e.target;
        if(t && t.id && String(t.id).startsWith('te-sf-track-')){{
          if(t.id === 'te-sf-track-meth_diff'){{
            const mdCb = panel.querySelector('#te-sf-methshow-diff');
            if(mdCb){{ mdCb.checked = !!t.checked; }}
          }}
          invalidateLayoutCache();
          syncChromLayoutOrderFromTrackGroups();
          scheduleApplyDebounced(80);
          schedulePrefsSave();
        }}
      }});
      panel.addEventListener('keydown', (e) => {{
        const t = e.target;
        if(!(t && t.matches && t.matches('input[data-ymax-main="1"], input[data-ymax-sub]'))){{ return; }}
        if(e.key === 'Enter'){{
          e.preventDefault();
          applyTrackYMax();
        }}
      }});
      const ymaxApplyBtn = panel.querySelector('#te-sf-ymax-apply');
      if(ymaxApplyBtn && !ymaxApplyBtn._teBound){{
        ymaxApplyBtn._teBound = true;
        ymaxApplyBtn.addEventListener('click', () => {{
          applyTrackYMax();
        }});
      }}
      ['#te-sf-methcolor-cgtotal', '#te-sf-methcolor-cgcovered', '#te-sf-methcolor-cgpct', '#te-sf-methcolor-level'].forEach((sel) => {{
        const el = panel.querySelector(sel);
        if(!el){{ return; }}
        el.addEventListener('input', () => {{
          applyMethylationTrackColors();
          schedulePrefsSave();
        }});
        el.addEventListener('change', () => {{
          applyMethylationTrackColors();
          schedulePrefsSave();
        }});
      }});
      ['#te-sf-methdiffcolor-a', '#te-sf-methdiffcolor-b'].forEach((sel) => {{
        const el = panel.querySelector(sel);
        if(!el){{ return; }}
        el.addEventListener('input', () => {{
          applyMethDiffTrackColors();
          schedulePrefsSave();
        }});
        el.addEventListener('change', () => {{
          applyMethDiffTrackColors();
          schedulePrefsSave();
        }});
      }});
    }}
    const title = panel.querySelector('.te-sf-title');
    let dragging = false;
    let dragDx = 0;
    let dragDy = 0;
    title.addEventListener('mousedown', (e) => {{
      dragging = true;
      const rect = panel.getBoundingClientRect();
      panel.style.left = rect.left + 'px';
      panel.style.top = rect.top + 'px';
      panel.style.right = 'auto';
      dragDx = e.clientX - rect.left;
      dragDy = e.clientY - rect.top;
      e.preventDefault();
    }});
    window.addEventListener('mousemove', (e) => {{
      if(!dragging){{ return; }}
      panel.style.left = Math.max(4, e.clientX - dragDx) + 'px';
      panel.style.top = Math.max(4, e.clientY - dragDy) + 'px';
    }});
    window.addEventListener('mouseup', () => {{ dragging = false; }});
    document.body.appendChild(panel);
    panel.style.zIndex = '100000';
    const exportName = (document.title || (GENOME_MODE ? 'te-genome-map' : 'te-divergence-plot'))
      .replace(/[^a-z0-9._-]+/gi, '_')
      .replace(/^_+|_+$/g, '') || 'te-plot';
    panel.querySelector('#te-sf-save-svg').addEventListener('click', () => {{
      Plotly.downloadImage(gd, {{format: 'svg', filename: exportName}});
    }});
    panel.querySelector('#te-sf-save-pdf').addEventListener('click', () => {{
      window.print();
    }});
    gd.on('plotly_relayout', (ev) => {{
      if(gd._teRelayoutSync || (gd._teApplyDepth || 0) > 0){{ return; }}
      if(GENOME_MODE){{
        if(!gd._teScaleMeta || !Array.isArray(gd._teScaleMeta.tracks)){{ return; }}
      }} else if(!gd._teDivergenceMeta){{
        return;
      }}

      const xAutorange = ev && ev['xaxis.autorange'] === true;
      const yAutorange = ev && ev['yaxis.autorange'] === true;
      let yNeedsRestore = GENOME_MODE && yAutorange;
      if(
        GENOME_MODE
        && !yNeedsRestore
        && ev
        && (ev['yaxis.range[0]'] !== undefined || ev['yaxis.range[1]'] !== undefined)
      ){{
        const yr0 = Number(ev['yaxis.range[0]'] ?? gd.layout?.yaxis?.range?.[0]);
        const yr1 = Number(ev['yaxis.range[1]'] ?? gd.layout?.yaxis?.range?.[1]);
        if(
          gd._teYRangeLocked
          && Number.isFinite(yr0)
          && Number.isFinite(yr1)
          && !xRangesClose([yr0, yr1], gd._teYRangeLocked)
        ){{
          yNeedsRestore = true;
        }}
      }}
      if(xAutorange || yNeedsRestore){{
        if(xAutorange){{
          gd._teXRangeZoom = null;
          gd._teUserXZoom = false;
        }} else if(yNeedsRestore){{
          captureXRangeFromRelayoutEvent(ev);
        }}
        if(GENOME_MODE && yNeedsRestore){{
          scheduleApply();
        }}
        schedulePrefsSave();
        if(xAutorange && !GENOME_MODE){{
          const layoutUpdate = {{}};
          applyXAxisZoomLayout(layoutUpdate);
          if(Object.keys(layoutUpdate).length){{ safeRelayout(layoutUpdate); }}
        }}
        return;
      }}

      const newRange = readXRangeFromRelayoutEvent(ev);
      if(!newRange){{ return; }}
      if(gd._teXRange0 && xRangesClose(newRange, gd._teXRange0)){{
        gd._teXRangeZoom = null;
        gd._teUserXZoom = false;
        syncXZoomSliderFromRange(gd._teXRange0[0], gd._teXRange0[1]);
        schedulePrefsSave();
        return;
      }}
      if(gd._teXRangeZoom && xRangesClose(newRange, gd._teXRangeZoom)){{ return; }}
      captureXRangeFromRelayoutEvent(ev);
      const snapped = gd._teXRangeZoom ? gd._teXRangeZoom.slice() : newRange;
      const layoutUpdate = {{
        'xaxis.range': snapped,
        'xaxis.autorange': false,
      }};
      if(GENOME_MODE && gd._teScaleMeta){{
        const ticks = zoomBinTickvals(snapped[0], snapped[1], gd._teScaleMeta);
        layoutUpdate['xaxis.tickmode'] = 'array';
        layoutUpdate['xaxis.tickvals'] = ticks.tickvals;
        layoutUpdate['xaxis.ticktext'] = ticks.ticktext;
      }}
      safeRelayout(layoutUpdate);
      schedulePrefsSave();
    }});
    gd._teSfScheduleApply = scheduleApply;
    gd._teSfInitialized = true;

    function genomePanelNeedsInitialApply(m){{
      if(!m){{ return false; }}
      const scaleEl = panel.querySelector('#te-sf-scale');
      const scale = Number(scaleEl ? scaleEl.value : (m.savedBarScale || 1));
      if(Number.isFinite(scale) && Math.abs(scale - 1) > 1e-6){{ return true; }}
      const rowEl = panel.querySelector('#te-sf-rowsp');
      if(rowEl){{
        const rowMult = Number(rowEl.value || 1);
        if(Number.isFinite(rowMult) && Math.abs(rowMult - 1) > 1e-6){{ return true; }}
      }}
      if(m.hasDiffTrack){{
        const diffEl = panel.querySelector('#te-sf-diffscale');
        const ds = Number(diffEl ? diffEl.value : (m.savedDiffTrackScale || 1));
        if(Number.isFinite(ds) && Math.abs(ds - 1) > 1e-6){{ return true; }}
      }}
      if(m.hasMethDiffTrack){{
        const mdsEl = panel.querySelector('#te-sf-methdiffscale');
        const mds = Number(mdsEl ? mdsEl.value : (m.savedMethDiffTrackScale || 1));
        if(Number.isFinite(mds) && Math.abs(mds - 1) > 1e-6){{ return true; }}
      }}
      const cl = m.chromLayout;
      if(cl && anyTrackLabelHidden(cl)){{ return true; }}
      if(hasCategoryFilters && Array.isArray(m.savedCategoryHidden) && m.savedCategoryHidden.length){{
        return true;
      }}
      if(Array.isArray(m.savedCategoryOrder) && m.savedCategoryOrder.length){{ return true; }}
      if(m.savedCategoryColors && typeof m.savedCategoryColors === 'object'
        && Object.keys(m.savedCategoryColors).length){{ return true; }}
      if(Array.isArray(m.savedTeSourceHidden) && m.savedTeSourceHidden.length){{ return true; }}
      if(Array.isArray(m.savedTrackGroupHidden) && m.savedTrackGroupHidden.length){{ return true; }}
      if(Array.isArray(m.savedTrackGroupOrder) && m.savedTrackGroupOrder.length){{ return true; }}
      if(m.savedTrackYMaxByGroup && typeof m.savedTrackYMaxByGroup === 'object'
        && Object.keys(m.savedTrackYMaxByGroup).length){{ return true; }}
      if(Number.isFinite(m.savedBarScale) && Math.abs(m.savedBarScale - 1) > 1e-6){{ return true; }}
      if(Number.isFinite(m.savedRowSpacing) && Math.abs(m.savedRowSpacing - 1) > 1e-6){{ return true; }}
      return false;
    }}

    gd._teGenomePanelNeedsInitialApply = genomePanelNeedsInitialApply;

    populateOriginalRawY();
    applySavedPanelState();
    if(GENOME_MODE){{
      if(gd._teScaleMeta){{
        gd._teSfMetaApplied = true;
        if(genomePanelNeedsInitialApply(gd._teScaleMeta)){{
          scheduleApply();
        }}
      }}
    }} else if(hasCategoryFilters){{
      scheduleApply();
    }}
    return true;
  }}
  function boot(retries){{
    const gd = document.querySelector('.plotly-graph-div');
    if(gd && typeof gd.on === 'function' && !gd._teSfAfterPlotHook){{
      gd._teSfAfterPlotHook = true;
      gd.on('plotly_afterplot', () => {{
        init();
      }});
    }}
    if(init()){{ return; }}
    if(retries <= 0){{ return; }}
    setTimeout(() => boot(retries - 1), 250);
  }}
  if(document.readyState === 'loading'){{
    document.addEventListener('DOMContentLoaded', () => boot(80));
  }} else {{
    boot(80);
  }}
}})();
</script>
"""


_DIVERGENCE_SCALE_APPLY_JS = """
      function divergenceVisibleStackMax(){
        const nBins = barTraceIndices.length > 0 ? (originalY.get(barTraceIndices[0]) || []).length : 0;
        if(nBins <= 0){ return 1; }
        const stackSums = new Array(nBins).fill(0);
        const groupsInUse = hasCategoryFilters ? stackGroupsInOrder() : groups;
        groupsInUse.forEach((g) => {
          (barsByGroup.get(g) || []).forEach((idx) => {
            const ys0 = originalY.get(idx) || [];
            for(let i = 0; i < ys0.length; i += 1){
              stackSums[i] += ys0[i];
            }
          });
        });
        let maxStack = 0;
        for(let i = 0; i < stackSums.length; i += 1){
          if(stackSums[i] > maxStack){ maxStack = stackSums[i]; }
        }
        return Math.max(maxStack, 1);
      }

      const yByIdx = new Map();
      const baseByIdx = new Map();
      const rawMax = divergenceVisibleStackMax();
      const baselineTop = Math.max(1, rawMax * 1.08);
      const top = Math.max(1, baselineTop / scale);
      const layoutUpdate = {
        'yaxis.range': [0, top],
        'yaxis.autorange': false,
      };

      if(hasCategoryFilters){
        const groupsInUse = stackGroupsInOrder();
        const nBins = barTraceIndices.length > 0 ? (originalY.get(barTraceIndices[0]) || []).length : 0;
        const cumulative = new Array(nBins).fill(0);
        groupsInUse.forEach((g) => {
          (barsByGroup.get(g) || []).forEach((idx) => {
            const ys0 = originalY.get(idx) || [];
            const ys = new Array(ys0.length);
            const bases = new Array(ys0.length);
            for(let i = 0; i < ys0.length; i += 1){
              bases[i] = cumulative[i];
              const v = Number(ys0[i] || 0);
              ys[i] = v;
              cumulative[i] += v;
            }
            yByIdx.set(idx, ys);
            baseByIdx.set(idx, bases);
          });
        });
        barTraceIndices.forEach((idx) => {
          const g = gd.data[idx].legendgroup || gd.data[idx].name || '';
          if(!g){ return; }
          const cb = checkByGroup.get(g);
          const on = !!(cb && cb.checked);
          if(on){ return; }
          const n = (originalY.get(idx) || []).length;
          yByIdx.set(idx, new Array(n).fill(0));
          baseByIdx.set(idx, new Array(n).fill(0));
        });
      } else {
        barTraceIndices.forEach((idx) => {
          const ys0 = originalY.get(idx) || [];
          yByIdx.set(idx, ys0.slice());
          baseByIdx.set(idx, new Array(ys0.length).fill(0));
        });
      }
"""

_GENOME_MAP_SCALE_APPLY_JS = """
      const layoutUpdate = {};
      const rowInput = panel.querySelector('#te-sf-rowsp');
      let rowMult = Number(rowInput ? rowInput.value : 1);
      if(!Number.isFinite(rowMult) || rowMult <= 0){ rowMult = 1; }
      rowMult = Math.max(__ROW_SPACING_MIN__, Math.min(__ROW_SPACING_MAX__, rowMult));
      gd._teRowMult = rowMult;
      const m = gd._teScaleMeta;
      if(m){{
        const byGroup = m._appliedTrackYMaxByGroup || m._panelTrackYMaxByGroup || m.trackYMaxByGroup || {{}};
        m._panelTrackYMaxByGroup = byGroup;
        m.trackYMaxByGroup = byGroup;
      }}
      let diffScale = 1;
      if(m && m.hasDiffTrack){{
        const diffInput = panel.querySelector('#te-sf-diffscale');
        diffScale = Number(diffInput ? diffInput.value : 1);
        if(!Number.isFinite(diffScale) || diffScale <= 0){{ diffScale = 1; }}
        diffScale = Math.max(
          __DIFF_TRACK_SCALE_MIN__,
          Math.min(__DIFF_TRACK_SCALE_MAX__, diffScale),
        );
      }}
      gd._teDiffScale = diffScale;
      let methDiffScale = 1;
      if(m && m.hasMethDiffTrack){{
        const methDiffInput = panel.querySelector('#te-sf-methdiffscale');
        methDiffScale = Number(methDiffInput ? methDiffInput.value : 1);
        if(!Number.isFinite(methDiffScale) || methDiffScale <= 0){{ methDiffScale = 1; }}
        methDiffScale = Math.max(
          __DIFF_TRACK_SCALE_MIN__,
          Math.min(__DIFF_TRACK_SCALE_MAX__, methDiffScale),
        );
      }}
      gd._teMethDiffScale = methDiffScale;
      const cl = gd._teScaleMeta && gd._teScaleMeta.chromLayout;
      const fastPath = String(gd._teApplyFastPath || '');
      gd._teApplyFastPath = '';
      categoryOnly = fastPath === 'category' && hasCategoryFilters;
      methylationOnly = fastPath === 'methylation';
      const heightScaleCache = new Map();
      function cachedHeightScale(chrom){{
        if(!heightScaleCache.has(chrom)){{
          heightScaleCache.set(chrom, heightScaleForChrom(chrom));
        }}
        return heightScaleCache.get(chrom);
      }}

      function teInputKindForChrom(chrom){{
        const t = trackMetaForChrom(chrom);
        if(!t){{ return 'masking'; }}
        if(t.isGeneTrack){{ return 'gene'; }}
        if(t.isMethylationTrack){{ return 'methylation'; }}
        return t.teInputKind || 'masking';
      }}

      function trackMetaForChrom(chrom){{
        const label = resolveTrackLabel(chrom);
        if(!m || !Array.isArray(m.tracks)){{ return null; }}
        return m.tracks.find((t) => t.chrom === label || String(t.chrom) === label) || null;
      }}

      function chromHasYMaxOverride(chrom){{
        const t = trackMetaForChrom(chrom);
        return !!(t && hasTrackYMaxOverride(t));
      }}

      function isGenePlainBarIdx(idx){{
        if(!plainBarIndices.includes(idx)){{ return false; }}
        const tr = gd.data[idx];
        if(!tr || tr.legendgroup === '__te_diff__' || tr.legendgroup === '__meth_diff__'){{ return false; }}
        const cd = tr.customdata && tr.customdata[0];
        const label = cd && cd[0];
        if(label == null){{ return false; }}
        if(typeof gd._teIsGeneTrackLabel === 'function' && gd._teIsGeneTrackLabel(String(label))){{
          return true;
        }}
        const resolved = resolveTrackLabel(String(label));
        const kinds = cl && cl.trackKindsByLabel;
        if(kinds && kinds[resolved] === 'gene'){{ return true; }}
        const tm = trackMetaForLabel(resolved);
        return !!(tm && tm.isGeneTrack);
      }}

      function normMaxForInputKind(kind){{
        if(kind === 'methylation'){{
          return Number(
            (m.methylationDensityNormMaxBySubtrack && m.methylationDensityNormMaxBySubtrack.cg_count_shared)
            || (m.methylationDensityNormMaxBySubtrack && m.methylationDensityNormMaxBySubtrack.cg_total)
            || 1,
          );
        }}
        if(kind === 'gene'){{
          return Number(m.geneDensityNormMax || 1);
        }}
        if(kind === 'consensus'){{
          return Number(
            m.consensusDensityNormMax || m.densityNormMax || m.densityMax || 1,
          );
        }}
        return Number(
          m.maskingDensityNormMax || m.densityNormMax || m.densityMax || 1,
        );
      }}

      function normMaxForChrom(chrom){{
        const t = trackMetaForChrom(chrom);
        if(t && hasTrackYMaxOverride(t)){{
          return densityNormForTrackMeta(t, m);
        }}
        return normMaxForInputKind(teInputKindForChrom(chrom));
      }}

      function rawFromY(yVal, chrom){{
        const t = trackMetaForChrom(chrom);
        const normMax = t
          ? densityNormForTrackMeta(t, m)
          : normMaxForChrom(chrom);
        if(!m || !normMax || !m.maxBarHeight){{ return Number(yVal || 0); }}
        return Number(yVal || 0) * normMax / m.maxBarHeight;
      }}

      function rawBinForPlainPoint(tr, idx, i, ys0i, trackLabel){{
        const cd = tr && Array.isArray(tr.customdata) ? tr.customdata : null;
        const row = cd && cd[i];
        if(row && row.length > 5){{
          const v = Number(row[5]);
          if(Number.isFinite(v)){{ return v; }}
        }}
        const store = originalRawY.get(idx);
        if(store && store[i] !== undefined){{ return store[i]; }}
        const t = trackMetaForLabel(trackLabel);
        if(t && (t.isGeneTrack || t.isMethylationTrack)){{
          return rawDensityFromRendered(ys0i, t);
        }}
        return rawFromY(ys0i, trackLabel);
      }}

      function parseChromFromKey(key){{
        return resolveTrackLabel(String(key || '').split('|')[0] || '');
      }}

      function chromGround(chrom){{
        if(!cl || !Array.isArray(cl.order)){{ return 0; }}
        const label = resolveTrackLabel(chrom);
        const idx = cl.order.indexOf(label);
        if(idx < 0){{ return 0; }}
        return trackBarBase(label, idx, cl, rowMult);
      }}

      function diffNormMaxForLabel(trackLabel){{
        if(isTeDiffTrackLabel(trackLabel)){{
          return teDiffNormMax();
        }}
        if(isMethDiffTrackLabel(trackLabel)){{
          return methDiffNormMax();
        }}
        if(m && m.sharedDiffAxis !== false){{
          return Math.max(Number(m.diffDensityNormMax || 1), 1e-9);
        }}
        if(!m || !Array.isArray(m.tracks)){{
          return Math.max(Number(m && m.diffDensityNormMax) || 1, 1e-9);
        }}
        const entry = m.tracks.find((tr) => tr.chrom === trackLabel && tr.isDiffTrack);
        return Math.max(
          Number(entry && (entry.diffNormMax || entry.densityNormMax)) || 0,
          Number(m.diffDensityNormMax) || 1,
          1e-9,
        );
      }}

      function computeVisibleMaxByChrom(){{
        const byChrom = new Map();
        function bump(chrom, raw){{
          if(!chrom){{ return; }}
          byChrom.set(chrom, Math.max(byChrom.get(chrom) || 0, raw));
        }}
        const groupsInUse = hasCategoryFilters ? stackGroupsInOrder() : groups;
        if(hasCategoryFilters && !methylationOnly){{
          const binTotals = new Map();
          groupsInUse.forEach((g) => {{
            (barsByGroup.get(g) || []).forEach((idx) => {{
              const ys0 = originalY.get(idx) || [];
              const keys = pointKeys.get(idx) || [];
              for(let i = 0; i < ys0.length; i += 1){{
                const key = keys[i] || (String(idx) + '|' + String(i));
                const raw = rawFromY(ys0[i], parseChromFromKey(key));
                binTotals.set(key, (binTotals.get(key) || 0) + raw);
              }}
            }});
          }});
          binTotals.forEach((v, key) => {{
            const chrom = parseChromFromKey(key);
            if(cl && !isTrackLabelVisible(chrom, cl)){{ return; }}
            bump(chrom, v);
          }});
        }} else {{
          groupsInUse.forEach((g) => {{
            (barsByGroup.get(g) || []).forEach((idx) => {{
              const ys0 = originalY.get(idx) || [];
              const keys = pointKeys.get(idx) || [];
              for(let i = 0; i < ys0.length; i += 1){{
                const chrom = parseChromFromKey(keys[i] || '');
                if(cl && !isTrackLabelVisible(chrom, cl)){{ continue; }}
                bump(chrom, rawFromY(ys0[i], chrom));
              }}
            }});
          }});
        }}
        plainBarIndices.forEach((idx) => {{
          const tr = gd.data[idx];
          const ys0 = originalY.get(idx) || [];
          const keys = pointKeys.get(idx) || [];
          for(let i = 0; i < ys0.length; i += 1){{
            const chrom = parseChromFromKey(keys[i] || '');
            if(cl && !isTrackLabelVisible(chrom, cl)){{ continue; }}
            bump(chrom, rawBinForPlainPoint(tr, idx, i, ys0[i], chrom));
          }}
        }});
        byChrom.forEach((v, chrom) => {{
          const floor = normMaxForChrom(chrom) * 1e-6;
          byChrom.set(chrom, Math.max(v, floor));
        }});
        return byChrom;
      }}

      function computeVisibleStackRawMax(){{
        const byChrom = computeVisibleMaxByChrom();
        let maxV = 0;
        byChrom.forEach((v) => {{ if(v > maxV){{ maxV = v; }} }});
        return Math.max(maxV, 1e-9);
      }}

      const visibleMaxByChrom = computeVisibleMaxByChrom();

      function allCategoryGroupsVisible(){{
        if(!hasCategoryFilters){{ return true; }}
        return groups.every((g) => {{
          const cb = checkByGroup.get(g);
          return !cb || cb.checked;
        }});
      }}

      function globalTeNormMaxForKind(kind){{
        return normMaxForInputKind(kind);
      }}

      function globalTeVisibleMaxForKind(kind){{
        let maxV = 0;
        if(m && Array.isArray(m.tracks)){{
          m.tracks.forEach((tr) => {{
            if(tr.isGeneTrack || tr.isDiffTrack){{ return; }}
            const trKind = tr.teInputKind || 'masking';
            if(trKind !== kind){{ return; }}
            const v = visibleMaxByChrom.get(tr.chrom) || 0;
            if(v > maxV){{ maxV = v; }}
          }});
        }}
        return Math.max(maxV, globalTeNormMaxForKind(kind) * 1e-6);
      }}

      function globalGeneVisibleMax(){{
        let maxV = 0;
        if(m && Array.isArray(m.tracks)){{
          m.tracks.forEach((tr) => {{
            if(!tr.isGeneTrack){{ return; }}
            const v = visibleMaxByChrom.get(tr.chrom) || 0;
            if(v > maxV){{ maxV = v; }}
          }});
        }}
        return Math.max(maxV, Number(m.geneDensityNormMax || 1) * 1e-6);
      }}

      function activeDensityMaxForKind(kind){{
        if(kind === 'gene'){{
          if(allCategoryGroupsVisible()){{
            return Number(m.geneDensityNormMax || 1);
          }}
          return globalGeneVisibleMax();
        }}
        if(allCategoryGroupsVisible()){{
          return globalTeNormMaxForKind(kind);
        }}
        return globalTeVisibleMaxForKind(kind);
      }}

      function globalTeNormMax(){{
        return globalTeNormMaxForKind('masking');
      }}

      function globalTeVisibleMax(){{
        return globalTeVisibleMaxForKind('masking');
      }}

      function activeTeDensityMax(){{
        return activeDensityMaxForKind('masking');
      }}

      gd._teVisibleDensityMax = activeTeDensityMax();
      gd._teActiveDensityMax = gd._teVisibleDensityMax;
      if(m){{
        m._activeTickMaxByKind = {{
          masking: activeDensityMaxForKind('masking'),
          consensus: activeDensityMaxForKind('consensus'),
          gene: activeDensityMaxForKind('gene'),
        }};
      }}

      function teAxisHeightForTrackLabel(label){{
        if(!cl || !Array.isArray(cl.order) || !m.maxBarHeight){{ return m.maxBarHeight * scale; }}
        const idx = cl.order.indexOf(label);
        if(idx < 0){{ return m.maxBarHeight * scale; }}
        const y0 = trackBarBase(label, idx, cl, rowMult);
        return densityAxisHeightForTrack(idx, y0, cl, rowMult, scale, m.maxBarHeight);
      }}

      function cappedTrackHeight(raw, chrom){{
        const label = resolveTrackLabel(chrom);
        const axisH = teAxisHeightForTrackLabel(label);
        const h = Number(raw || 0) * cachedHeightScale(label);
        if(!Number.isFinite(h)){{ return 0; }}
        return Math.max(0, Math.min(h, axisH));
      }}

      function heightScaleForChrom(chrom){{
        if(!m || !m.maxBarHeight){{ return scale; }}
        const label = resolveTrackLabel(chrom);
        const kind = teInputKindForChrom(label);
        const t = trackMetaForChrom(label);
        if(t && hasTrackYMaxOverride(t)){{
          const axisH = teAxisHeightForTrackLabel(label);
          const normMax = densityNormForTrackMeta(t, m);
          return axisH / Math.max(normMax, 1e-9);
        }}
        if(kind === 'gene' && Number(m.geneDensityNormMax) > 0){{
          const axisH = teAxisHeightForTrackLabel(label);
          return axisH / Number(m.geneDensityNormMax);
        }}
        if(kind === 'methylation' && t){{
          const sub = t.methylationSubtrack || 'cg_total';
          const bySub = m.methylationDensityNormMaxBySubtrack || {{}};
          let normMax = 1;
          if(sub === 'cg_total' || sub === 'cg_covered'){{
            normMax = Number(bySub.cg_count_shared || bySub.cg_total || bySub.cg_covered || 1);
          }} else if(sub === 'cg_coverage_pct'){{
            normMax = Number(bySub.cg_coverage_pct || 1);
          }} else if(sub === 'meth_level'){{
            normMax = Number(bySub.meth_level || 1);
          }} else {{
            normMax = Number(bySub[sub] || bySub.cg_total || 1);
          }}
          const axisH = teAxisHeightForTrackLabel(label);
          return axisH / Math.max(normMax, 1e-9);
        }}
        const useSharedTe = !m || m.sharedTeDensityAxis !== false;
        if(
          useSharedTe && t && !t.isGeneTrack && !t.isMethylationTrack && !t.isDiffTrack
          && !hasTrackYMaxOverride(t)
        ){{
          const axisH = teAxisHeightForTrackLabel(label);
          return axisH / activeDensityMaxForKind(kind);
        }}
        const normMax = normMaxForChrom(label);
        let visMax = normMax;
        if(!hasTrackYMaxOverride(t)){{
          visMax = visibleMaxByChrom.get(label) || normMax;
        }}
        const axisH = teAxisHeightForTrackLabel(label);
        return axisH / Math.max(visMax, 1e-9);
      }}

      const yByIdx = new Map();
      const baseByIdx = new Map();
      let maxStackTop = 0;
      const useSharedTeAxis = !m || m.sharedTeDensityAxis !== false;
      function builtStackOrder(){{
        const built = (m && Array.isArray(m.categoryOrder) && m.categoryOrder.length)
          ? m.categoryOrder.slice()
          : groups.slice();
        return built.slice().reverse().filter((g) => {{
          const cb = checkByGroup.get(g);
          return !!(cb && cb.checked);
        }});
      }}

      function stackOrderMatchesBuilt(){{
        if(!hasCategoryFilters){{ return true; }}
        const current = stackGroupsInOrder();
        const built = builtStackOrder();
        if(current.length !== built.length){{ return false; }}
        for(let i = 0; i < current.length; i += 1){{
          if(current[i] !== built[i]){{ return false; }}
        }}
        return true;
      }}

      const preserveOriginalHeights =
        useSharedTeAxis &&
        scale === 1 &&
        rowMult === 1 &&
        diffScale === 1 &&
        methDiffScale === 1 &&
        allCategoryGroupsVisible() &&
        stackOrderMatchesBuilt() &&
        !anyTrackLabelHidden(cl) &&
        !anyTrackGroupYMaxOverride() &&
        !(m && m._yMaxRescaleActive);

      let minStackExtent = Number.POSITIVE_INFINITY;
      function noteStackTop(yVal, baseVal){{
        const base = Number(baseVal || 0);
        const y = Number(yVal || 0);
        const top = base + y;
        if(top > maxStackTop){{ maxStackTop = top; }}
        const bottom = base + Math.min(0, y);
        if(bottom < minStackExtent){{ minStackExtent = bottom; }}
      }}

      if(preserveOriginalHeights){{
        barTraceIndices.forEach((idx) => {{
          const ys0 = originalY.get(idx) || [];
          const bs0 = originalBase.get(idx) || [];
          yByIdx.set(idx, ys0.slice());
          baseByIdx.set(idx, bs0.slice());
          for(let i = 0; i < ys0.length; i += 1){{
            noteStackTop(ys0[i], bs0[i]);
          }}
        }});
      }} else if(hasCategoryFilters){{
      const cumulative = new Map();

        const groupsInUse = stackGroupsInOrder();
        groupsInUse.forEach((g) => {{
          (barsByGroup.get(g) || []).forEach((idx) => {{
            const ys0 = originalY.get(idx) || [];
            const keys = pointKeys.get(idx) || [];
            const ys = new Array(ys0.length);
            const bases = new Array(ys0.length);
            for(let i = 0; i < ys0.length; i += 1){{
              const key = keys[i] || (String(idx) + '|' + String(i));
              const chrom = parseChromFromKey(key);
              if(cl && !isTrackLabelVisible(chrom, cl)){{
                ys[i] = 0;
                bases[i] = 0;
                continue;
              }}
              const ground = cumulative.has(key) ? cumulative.get(key) : chromGround(chrom);
              const rawArr = originalRawY.get(idx);
              const raw = rawArr && rawArr[i] !== undefined
                ? rawArr[i]
                : rawFromY(ys0[i], chrom);
              const h = cappedTrackHeight(raw, chrom);
              ys[i] = h;
              bases[i] = ground;
              cumulative.set(key, ground + h);
              if(ground + h > maxStackTop){{ maxStackTop = ground + h; }}
            }}
            yByIdx.set(idx, ys);
            baseByIdx.set(idx, bases);
          }});
        }});
        barTraceIndices.forEach((idx) => {{
          if(plainBarIndices.includes(idx)){{ return; }}
          const g = gd.data[idx].legendgroup || gd.data[idx].name || '';
          if(!g){{ return; }}
          const cb = checkByGroup.get(g);
          const on = !!(cb && cb.checked);
          if(!on){{
            const n = (originalY.get(idx) || []).length;
            yByIdx.set(idx, new Array(n).fill(0));
            baseByIdx.set(idx, (originalBase.get(idx) || []).slice());
          }}
        }});
        cumulative.forEach((top) => {{
          if(top > maxStackTop){{ maxStackTop = top; }}
        }});
      }} else if(!preserveOriginalHeights){{
        barTraceIndices.forEach((idx) => {{
          if(plainBarIndices.includes(idx)){{ return; }}
          const tr = gd.data[idx];
          if(!tr){{ return; }}
          const g = tr.legendgroup || '';
          if(!g || !barsByGroup.has(g)){{ return; }}
          const ys0 = originalY.get(idx) || [];
          const keys = pointKeys.get(idx) || [];
          const ys = new Array(ys0.length);
          const bases = new Array(ys0.length);
          for(let i = 0; i < ys0.length; i += 1){{
            const chrom = parseChromFromKey(keys[i] || '');
            if(cl && !isTrackLabelVisible(chrom, cl)){{
              ys[i] = 0;
              bases[i] = 0;
              continue;
            }}
            const ground = chromGround(chrom);
            const rawArr = originalRawY.get(idx);
            const raw = rawArr && rawArr[i] !== undefined
              ? rawArr[i]
              : rawBinForPlainPoint(tr, idx, i, ys0[i], chrom);
            const h = cappedTrackHeight(raw, chrom);
            ys[i] = h;
            bases[i] = ground;
            if(ground + h > maxStackTop){{ maxStackTop = ground + h; }}
          }}
          yByIdx.set(idx, ys);
          baseByIdx.set(idx, bases);
        }});
      }}

      plainBarIndices.forEach((idx) => {{
        const tr = gd.data[idx];
        if(tr && tr.legendgroup === '__te_diff__'){{
          const ys0 = originalY.get(idx) || [];
          const keys = pointKeys.get(idx) || [];
          const diffHalf0 = Math.max(0.12, Number(m.diffBarHalf || m.maxBarHeight * 0.36));
          const ys = new Array(ys0.length);
          const bases = new Array(ys0.length);
          for(let i = 0; i < ys0.length; i += 1){{
            const trackLabel = parseChromFromKey(keys[i] || '');
            if(cl && !isTrackLabelVisible(trackLabel, cl)){{
              ys[i] = 0;
              bases[i] = 0;
              continue;
            }}
            const diffMax = diffNormMaxForLabel(trackLabel);
            const resolved = resolveTrackLabel(trackLabel);
            const chromIdx = cl && cl.order ? cl.order.indexOf(resolved) : -1;
            const center = trackYCenter(resolved, chromIdx, cl, rowMult);
            const rawArr = originalRawY.get(idx);
            let rawBp = rawArr && rawArr[i] !== undefined
              ? Number(rawArr[i] || 0)
              : Math.abs(Number(ys0[i] || 0)) * Math.max(Number(m.diffDensityNormMax || 1), 1e-9) / diffHalf0;
            rawBp = Math.min(rawBp, diffMax);
            const sign = Number(ys0[i] || 0) < 0 ? -1 : 1;
            ys[i] = sign * (rawBp / diffMax) * diffHalf0 * rowMult * diffScale;
            bases[i] = center;
            noteStackTop(Math.abs(ys[i]), center - Math.abs(ys[i]));
            noteStackTop(ys[i], center);
          }}
          yByIdx.set(idx, ys);
          baseByIdx.set(idx, bases);
          return;
        }}
        if(tr && tr.legendgroup === '__meth_diff__'){{
          const ys0 = originalY.get(idx) || [];
          const keys = pointKeys.get(idx) || [];
          const diffHalf0 = Math.max(0.12, Number(m.methDiffBarHalf || m.maxBarHeight * 0.36));
          const ys = new Array(ys0.length);
          const bases = new Array(ys0.length);
          for(let i = 0; i < ys0.length; i += 1){{
            const trackLabel = parseChromFromKey(keys[i] || '');
            if(cl && !isTrackLabelVisible(trackLabel, cl)){{
              ys[i] = 0;
              bases[i] = 0;
              continue;
            }}
            const diffMax = methDiffNormMax();
            const resolved = resolveTrackLabel(trackLabel);
            const chromIdx = cl && cl.order ? cl.order.indexOf(resolved) : -1;
            const center = trackYCenter(resolved, chromIdx, cl, rowMult);
            const rawVal = Math.abs(Number(ys0[i] || 0)) * diffMax / diffHalf0;
            const sign = Number(ys0[i] || 0) < 0 ? -1 : 1;
            ys[i] = sign * (rawVal / diffMax) * diffHalf0 * rowMult * methDiffScale;
            bases[i] = center;
            noteStackTop(Math.abs(ys[i]), center - Math.abs(ys[i]));
            noteStackTop(ys[i], center);
          }}
          yByIdx.set(idx, ys);
          baseByIdx.set(idx, bases);
          return;
        }}
        if(
          isGenePlainBarIdx(idx)
          && scale === 1
          && rowMult === 1
          && diffScale === 1
          && methDiffScale === 1
          && !anyTrackGroupYMaxOverride()
          && !(m && m._yMaxRescaleActive)
        ){{
          const ys0 = originalY.get(idx) || [];
          const keys = pointKeys.get(idx) || [];
          const ys = new Array(ys0.length);
          const bases = new Array(ys0.length);
          for(let i = 0; i < ys0.length; i += 1){{
            const chrom = parseChromFromKey(keys[i] || '');
            if(cl && !isTrackLabelVisible(chrom, cl)){{
              ys[i] = 0;
              bases[i] = 0;
              continue;
            }}
            const ground = chromGround(chrom);
            const origBase = Number((originalBase.get(idx) || [])[i] ?? 0);
            ys[i] = Number(ys0[i] || 0);
            bases[i] = ground !== 0 ? ground : origBase;
            if(bases[i] + ys[i] > maxStackTop){{ maxStackTop = bases[i] + ys[i]; }}
          }}
          yByIdx.set(idx, ys);
          baseByIdx.set(idx, bases);
          return;
        }}
        if(isGenePlainBarIdx(idx)){{
          const ys0 = originalY.get(idx) || [];
          const keys = pointKeys.get(idx) || [];
          const ys = new Array(ys0.length);
          const bases = new Array(ys0.length);
          for(let i = 0; i < ys0.length; i += 1){{
            const chrom = parseChromFromKey(keys[i] || '');
            if(cl && !isTrackLabelVisible(chrom, cl)){{
              ys[i] = 0;
              bases[i] = 0;
              continue;
            }}
            const ground = chromGround(chrom);
            const raw = rawBinForPlainPoint(tr, idx, i, ys0[i], chrom);
            const tmeta = trackMetaForLabel(chrom);
            const normMax = tmeta
              ? densityNormForTrackMeta(tmeta, m)
              : Number(m.geneDensityNormMax || 1);
            const h = raw * (Number(m.maxBarHeight || 1) * scale / Math.max(normMax, 1e-9));
            ys[i] = h;
            bases[i] = ground !== 0 ? ground : Number((originalBase.get(idx) || [])[i] ?? 0);
            if(bases[i] + ys[i] > maxStackTop){{ maxStackTop = bases[i] + ys[i]; }}
          }}
          yByIdx.set(idx, ys);
          baseByIdx.set(idx, bases);
          return;
        }}
        const ys0 = originalY.get(idx) || [];
        const keys = pointKeys.get(idx) || [];
        const ys = new Array(ys0.length);
        const bases = new Array(ys0.length);
        for(let i = 0; i < ys0.length; i += 1){{
          const chrom = parseChromFromKey(keys[i] || '');
          const ground = chromGround(chrom);
          if(cl && !isTrackLabelVisible(chrom, cl)){{
            ys[i] = 0;
            bases[i] = 0;
            continue;
          }}
          const raw = rawBinForPlainPoint(tr, idx, i, ys0[i], chrom);
          const h = cappedTrackHeight(raw, chrom);
          ys[i] = h;
          bases[i] = ground;
          if(ground + h > maxStackTop){{ maxStackTop = ground + h; }}
        }}
        yByIdx.set(idx, ys);
        baseByIdx.set(idx, bases);
      }});

      const yRange0 = gd._teYRange0 || (gd.layout.yaxis && gd.layout.yaxis.range) || [-0.35, 10];
      let yMin;
      let axisTop;
      if(categoryOnly){{
        const prevRange = gd._teYRangeLocked
          || (gd.layout.yaxis && gd.layout.yaxis.range)
          || [Number(yRange0[0] ?? -0.35), maxStackTop + 0.35];
        yMin = Number(gd._teCachedYMin ?? prevRange[0] ?? -0.35);
        axisTop = Math.max(
          Number(gd._teCachedAxisTop ?? prevRange[1] ?? maxStackTop + 0.35),
          maxStackTop + 0.35,
        );
        const prevTickMax = Number(gd._teCachedVisibleDensityMax ?? gd._teVisibleDensityMax ?? 0);
        const tickDelta = Math.abs(prevTickMax - Number(gd._teVisibleDensityMax || 0));
        const tickMaxChanged = tickDelta > Math.max(prevTickMax, 1) * 1e-5;
        gd._teCachedVisibleDensityMax = gd._teVisibleDensityMax;
        gd._teCategoryFastRelayout = !tickMaxChanged;
        if(Math.abs(axisTop - Number(prevRange[1] || 0)) > 0.05){{
          layoutUpdate['yaxis.range'] = [yMin, axisTop];
          layoutUpdate['yaxis.autorange'] = false;
          layoutUpdate['yaxis.fixedrange'] = true;
          gd._teYRangeLocked = [yMin, axisTop];
          gd._teCategoryFastRelayout = false;
        }}
        gd._teCachedAxisTop = axisTop;
        gd._teCachedYMin = yMin;
        gd._teCachedMaxStackTop = maxStackTop;
        preserveXRange(layoutUpdate);
      }} else {{
      yMin = Number((cl && cl.yMin !== undefined) ? cl.yMin : (yRange0[0] ?? -0.35));
      if(cl && Array.isArray(cl.order)){{
        let lowestDiff = Number.POSITIVE_INFINITY;
        cl.order.forEach((label, i) => {{
          const kind = cl.trackKindsByLabel ? cl.trackKindsByLabel[label] : '';
          if(kind !== 'te_diff' && kind !== 'meth_diff'){{ return; }}
          if(!isTrackLabelVisible(label, cl)){{ return; }}
          const halfDiff = greyTrackHalf(label, i, cl, rowMult, diffScale, methDiffScale);
          const yc = trackYCenter(label, i, cl, rowMult);
          lowestDiff = Math.min(lowestDiff, yc - halfDiff);
        }});
        if(Number.isFinite(lowestDiff) && lowestDiff > -1e5){{
          yMin = Math.min(yMin, lowestDiff - (cl.diffAxisLabelPad || 0.18));
        }}
      }}
      if(minStackExtent < Number.POSITIVE_INFINITY){{
        yMin = Math.min(yMin, minStackExtent - 0.18);
      }}
      let axisTop = maxStackTop + 0.35;
      if(cl && Array.isArray(cl.order) && cl.order.length > 0){{
        const visibleOrder = visibleTrackLabels(cl);
        const n = visibleOrder.length || cl.order.length;
        const firstLabel = visibleOrder[0] || cl.order[0];
        const firstBase = trackBarBase(firstLabel, 0, cl, rowMult);
        const firstAxisH = densityAxisHeightForTrack(
          0, firstBase, cl, rowMult, scale, m.maxBarHeight,
        );
        const lastLabel = visibleOrder[visibleOrder.length - 1] || cl.order[cl.order.length - 1];
        const lastY = trackYCenter(lastLabel, cl.order.length - 1, cl, rowMult);
        axisTop = Math.max(
          maxStackTop + 0.35,
          firstBase + firstAxisH + 0.28 * rowMult,
          lastY + cl.trackHalf + cl.barGapBelow + m.maxBarHeight * scale + cl.yTopPad,
        );
        layoutUpdate['yaxis.showticklabels'] = false;
        layoutUpdate['yaxis.showline'] = false;
        layoutUpdate['yaxis.ticks'] = '';
        cl.order.forEach((chrom, i) => {{
          const visible = isTrackLabelVisible(chrom, cl);
          const yc = trackYCenter(chrom, i, cl, rowMult);
          const si = cl.greyShapeIdx ? cl.greyShapeIdx[i] : undefined;
          if(si !== undefined){{
            layoutUpdate['shapes[' + si + '].visible'] = visible;
            if(visible){{
              const half = greyTrackHalf(chrom, i, cl, rowMult, diffScale, methDiffScale);
              layoutUpdate['shapes[' + si + '].y0'] = yc - half;
              layoutUpdate['shapes[' + si + '].y1'] = yc + half;
            }}
          }}
          const liMap = cl.chromLabelAnnIdxByLabel || {{}};
          const li = liMap[chrom] !== undefined
            ? liMap[chrom]
            : (cl.chromLabelAnnIdx ? cl.chromLabelAnnIdx[i] : undefined);
          if(li !== undefined){{
            layoutUpdate['annotations[' + li + '].visible'] = visible;
            if(!visible){{ return; }}
            const kind = cl.trackKindsByLabel && cl.trackKindsByLabel[chrom];
            const isDiff = kind === 'te_diff' || kind === 'meth_diff';
            if(isDiff){{
              const half = greyTrackHalf(chrom, i, cl, rowMult, diffScale, methDiffScale);
              layoutUpdate['annotations[' + li + '].y'] = yc - half - 0.10 * rowMult;
              layoutUpdate['annotations[' + li + '].yanchor'] = 'top';
            }} else {{
              layoutUpdate['annotations[' + li + '].y'] = yc;
              layoutUpdate['annotations[' + li + '].yanchor'] = 'middle';
            }}
          }}
        }});
        layoutUpdate.height = Math.max(380, Math.round(68 * n * rowMult + 130));
      }}
      layoutUpdate['yaxis.range'] = [yMin, axisTop];
      layoutUpdate['yaxis.autorange'] = false;
      layoutUpdate['yaxis.fixedrange'] = true;
      gd._teYRangeLocked = [yMin, axisTop];
      gd._teCachedYMin = yMin;
      gd._teCachedAxisTop = axisTop;
      gd._teCachedMaxStackTop = maxStackTop;
      gd._teCachedVisibleDensityMax = gd._teVisibleDensityMax;
      preserveXRange(layoutUpdate);
      }}
"""


def _figure_genome_scale_meta(fig: go.Figure) -> dict[str, Any]:
    raw_meta = getattr(fig.layout, "meta", None) if fig.layout else None
    if not raw_meta:
        return {}
    if isinstance(raw_meta, dict):
        te_scale = raw_meta.get("teGenomeScale")
    elif hasattr(raw_meta, "get"):
        te_scale = raw_meta.get("teGenomeScale")  # type: ignore[union-attr]
    else:
        te_scale = getattr(raw_meta, "teGenomeScale", None)
    return dict(te_scale) if isinstance(te_scale, dict) else {}


def _inject_plotly_controls(
    html: str,
    *,
    scale_label: str,
    panel_mode: str,
    prefs_path: Path | None = None,
) -> str:
    prefs_str = str(prefs_path.expanduser().resolve()) if prefs_path else None
    controls = _plotly_interactive_controls_html(
        scale_label=scale_label,
        panel_mode=panel_mode,
        prefs_path=prefs_str,
    )
    if "</body>" in html:
        return html.replace("</body>", controls + "\n</body>")
    return html + controls


def blast_txt_path_for_query(query_fasta: Path, program: str = "blastn") -> Path:
    """Path to ``{fasta_basename}_{bn|bp|bx}.txt`` next to the query FASTA."""
    return local_blast_txt_path(query_fasta, program)


def _legacy_repeatviewer_dir(query_fasta: Path) -> Path:
    """Older runs wrote outputs under ``{stem}.repeatviewer/`` next to the FASTA."""
    query_fasta = query_fasta.expanduser().resolve()
    return query_fasta.parent / f"{query_fasta.stem}.repeatviewer"


def blast_table_name_suffixes() -> tuple[str, ...]:
    """Filename endings recognized as BLAST hit tables (canonical first)."""
    from repeatviewer.backends.blast_output_names import all_blast_table_txt_suffixes

    return all_blast_table_txt_suffixes()


def strip_blast_table_suffix(name: str) -> str | None:
    """If ``name`` is a BLAST table, return the query basename; else ``None``."""
    return _strip_blast_table_suffix(name)


def blast_table_candidates_for_query(query_fasta: Path) -> list[Path]:
    """Possible BLAST table paths for a query FASTA (canonical tags first)."""
    query_fasta = query_fasta.expanduser().resolve()
    parent = query_fasta.parent
    legacy = _legacy_repeatviewer_dir(query_fasta)
    names = blast_table_candidate_names(query_fasta.name, stem_only=True)
    seen: set[str] = set()
    candidates: list[Path] = []
    for base in (parent, legacy):
        for name in names:
            path = base / name
            key = str(path)
            if key in seen:
                continue
            seen.add(key)
            candidates.append(path)
    return candidates


def find_blast_table_for_query(query_fasta: Path) -> Path | None:
    """Return an existing BLAST table for the query FASTA, if present."""
    for path in blast_table_candidates_for_query(query_fasta):
        if path.is_file():
            return path
    return None


GENOME_MAP_HTML_SUFFIX = "_den_plot.html"
GENOME_MAP_REGION_HTML_SUFFIX = "_den_plot_region.html"
# Legacy suffixes kept for preview detection of older outputs.
_LEGACY_GENOME_MAP_SUFFIXES = (
    "_density_plot.html",
    "_density_plot_region.html",
    "_genome_map.html",
    "_genome_map_region.html",
)


def te_hit_source_tag(path: Path) -> str:
    """Return ``RM`` or a BLAST tag (``bn``, ``b``, …) for TE hit-table naming."""
    name = path.name.lower()
    if name.endswith(".out") and strip_blast_table_suffix(path.name) is None:
        return "RM"
    return blast_table_source_tag(path.name)


def te_plot_output_suffix(
    te_input_files: list[Path] | Path,
    *,
    kind: str,
    region: bool = False,
) -> str:
    """Build a plot HTML suffix: ``_{b|RM|multi}_{kind}_plot[.region].html``."""
    te_paths = (
        [te_input_files.expanduser().resolve()]
        if isinstance(te_input_files, Path)
        else [p.expanduser().resolve() for p in te_input_files if p]
    )
    if len(te_paths) > 1:
        tag = "multi"
    elif te_paths:
        tag = te_hit_source_tag(te_paths[0])
    else:
        tag = "b"
    kind = str(kind).strip("_") or "plot"
    # Multi-source density plots use ``dens`` (``_multi_dens_plot.html``).
    if tag == "multi" and kind in {"den", "density"}:
        kind = "dens"
    if region:
        return f"_{tag}_{kind}_plot_region.html"
    return f"_{tag}_{kind}_plot.html"


def is_blast_table_file(path: Path) -> bool:
    """True if ``path`` looks like a BLAST hit table from RepeatViewer/Blast."""
    return strip_blast_table_suffix(path.name) is not None


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


def is_divergence_plot_html_file(path: Path) -> bool:
    name = path.name.lower()
    if not name.endswith((".html", ".htm")):
        return False
    return (
        "_div_plot" in name
        or "_divergence_plot" in name
        or "_blast_divergence" in name
    )


def blast_genome_map_html_path_for_blast_table(blast_txt: Path) -> Path:
    """Expected density-plot HTML path for a BLAST table."""
    blast_txt = blast_txt.expanduser().resolve()
    stem = te_hit_table_stem(blast_txt)
    suffix = te_plot_output_suffix(blast_txt, kind="den")
    return blast_txt.parent / f"{stem}{suffix}"


def genome_map_html_path_for_te_inputs(
    te_input_files: list[Path] | Path,
    *,
    region: bool = False,
    gene_gff: Path | None = None,
) -> Path:
    """Expected density-plot HTML path for one or more TE hit-table inputs."""
    te_paths = (
        [te_input_files.expanduser().resolve()]
        if isinstance(te_input_files, Path)
        else [p.expanduser().resolve() for p in te_input_files if p]
    )
    if te_paths:
        primary = te_paths[0]
        stem = _genome_map_output_stem(te_paths)
        if gene_gff is not None:
            stem = f"{stem}_genes"
        suffix = te_plot_output_suffix(te_paths, kind="den", region=region)
    elif gene_gff is not None:
        primary = gene_gff.expanduser().resolve()
        stem = f"{primary.stem}_genes"
        suffix = (
            GENOME_MAP_REGION_HTML_SUFFIX if region else GENOME_MAP_HTML_SUFFIX
        )
    else:
        raise ValueError("At least one TE or gene input path is required")
    return primary.parent / f"{stem}{suffix}"


def blast_genome_map_html_path_for_query(query_fasta: Path) -> Path:
    query_fasta = query_fasta.expanduser().resolve()
    return query_fasta.parent / f"{query_fasta.name}{GENOME_MAP_HTML_SUFFIX}"


def divergence_plot_html_path_for_te_inputs(
    te_input_files: list[Path] | Path,
) -> Path:
    """Expected divergence-plot HTML path for one or more TE hit tables."""
    te_paths = (
        [te_input_files.expanduser().resolve()]
        if isinstance(te_input_files, Path)
        else [p.expanduser().resolve() for p in te_input_files if p]
    )
    if not te_paths:
        raise ValueError("At least one TE input path is required")
    primary = te_paths[0]
    stem = _genome_map_output_stem(te_paths)
    suffix = te_plot_output_suffix(te_paths, kind="div")
    return primary.parent / f"{stem}{suffix}"


def read_blast_table(path: Path) -> pd.DataFrame:
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)

    with path.open(encoding="utf-8", errors="replace") as handle:
        first = handle.readline().strip()
    if not first:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    if first.lower().startswith("qseqid"):
        df = pd.read_csv(path, sep="\t")
        missing = [c for c in BLAST_COLUMNS if c not in df.columns]
        if missing:
            raise ValueError(f"BLAST table missing columns: {missing}")
        return df[BLAST_COLUMNS].copy()

    return pd.read_csv(path, sep="\t", header=None, names=BLAST_COLUMNS)


REPEATMASKER_OUT_COLUMNS = [
    "score",
    "perc_div",
    "perc_del",
    "perc_ins",
    "query_name",
    "query_start",
    "query_end",
    "query_left",
    "strand",
    "repeat_name",
    "repeat_class",
    "repeat_start",
    "repeat_end",
    "repeat_left",
    "rm_id",
]


def _looks_like_repeatmasker_out_text(text: str) -> bool:
    """Detect RepeatMasker ``.out`` even when the header splits ``perc`` and ``div.`` across lines."""
    low = text.lower()
    has_div = "div." in low or "perc div" in low
    has_query = "query" in low or "sequence" in low
    has_repeat = "repeat" in low or "matching" in low
    if has_div and has_query and has_repeat:
        return True
    return bool(
        re.search(
            r"^\s*\d+\s+\d+(?:\.\d+)?\s+\d+(?:\.\d+)?\s+\d+(?:\.\d+)?\s+\S+",
            text,
            re.MULTILINE,
        )
    )


def is_repeatmasker_out_file(path: Path) -> bool:
    """True when ``path`` looks like a RepeatMasker ``*.out`` annotation table."""
    if not path.name.lower().endswith(".out"):
        return False
    try:
        with path.expanduser().resolve().open(encoding="utf-8", errors="replace") as handle:
            head = "".join(handle.readline() for _ in range(60))
    except OSError:
        return False
    if not head.strip():
        return False
    if head.lstrip().lower().startswith("qseqid"):
        return False
    return _looks_like_repeatmasker_out_text(head)


def is_te_hit_table_file(path: Path) -> bool:
    """BLAST hit table or RepeatMasker ``.out`` usable for TE plots."""
    return is_blast_table_file(path) or is_repeatmasker_out_file(path)


def te_hit_table_stem(path: Path) -> str:
    """Base name for plot HTML outputs derived from a hit table path."""
    name = path.name
    stripped = strip_blast_table_suffix(name)
    if stripped is not None:
        return stripped
    if name.lower().endswith(".out"):
        return name[: -len(".out")]
    return path.stem


def _compact_bismark_cov_label(name: str, *, max_len: int = 24) -> str:
    """Shorten Bismark coverage filenames for track labels and margins."""
    s = str(name).strip()
    if not s:
        return s
    lower = s.lower()
    for suffix in (
        ".deduplicated.bismark.cov.gz",
        ".deduplicated.bismark.cov",
        ".bismark.cov.gz",
        ".bismark.cov",
        ".cov.gz",
        ".cov",
    ):
        if lower.endswith(suffix):
            s = s[: -len(suffix)]
            lower = s.lower()
            break
    s = s.rstrip("._-")
    if len(s) > max_len:
        s = s[: max_len - 1].rstrip("._-") + "…"
    return s or str(name).strip()


def _disambiguate_path_labels(
    paths: list[Path],
    *,
    label_for_path: Callable[[Path], str],
) -> list[str]:
    """Return display labels that are unique even when basenames repeat."""
    raw = [str(label_for_path(p)) for p in paths]
    if len(set(raw)) == len(raw):
        return raw
    out: list[str] = []
    for path, label in zip(paths, raw, strict=True):
        if raw.count(label) == 1:
            out.append(label)
            continue
        parent = path.parent.name
        candidate = f"{label} ({parent})" if parent and parent not in {".", ""} else f"{label} [2]"
        suffix = 2
        while candidate in out:
            suffix += 1
            candidate = (
                f"{label} ({parent}) #{suffix}"
                if parent and parent not in {".", ""}
                else f"{label} #{suffix}"
            )
        out.append(candidate)
    return out


def _repeatmasker_data_start_line(lines: list[str]) -> int:
    for index, line in enumerate(lines):
        if re.match(
            r"^\s*\d+\s+\d+(?:\.\d+)?\s+\d+(?:\.\d+)?\s+\d+(?:\.\d+)?\s+\S+",
            line,
        ):
            return index
    for index, line in enumerate(lines):
        lower = line.lower()
        if ("div." in lower or "perc div" in lower) and (
            "query" in lower or "sequence" in lower
        ):
            return index + 1
    return 3


def _repeatmasker_data_start_line_from_file(path: Path, max_header_lines: int = 400) -> int:
    """Locate the first data row without loading the whole RepeatMasker ``.out`` file."""
    try:
        with path.open(encoding="utf-8", errors="replace") as handle:
            header_lines: list[str] = []
            for index, line in enumerate(handle):
                if re.match(
                    r"^\s*\d+\s+\d+(?:\.\d+)?\s+\d+(?:\.\d+)?\s+\d+(?:\.\d+)?\s+\S+",
                    line,
                ):
                    return index
                if index < max_header_lines:
                    header_lines.append(line)
                elif header_lines:
                    break
    except OSError:
        return 3
    return _repeatmasker_data_start_line(header_lines)


def read_repeatmasker_out(path: Path) -> pd.DataFrame:
    """Read RepeatMasker ``.out`` and convert to BLAST-shaped columns for ``prepare_hits``."""
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    if not is_repeatmasker_out_file(path):
        raise ValueError(f"Not a RepeatMasker .out file: {path}")

    data_start = _repeatmasker_data_start_line_from_file(path)
    try:
        rm = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            skiprows=data_start,
            usecols=range(15),
            on_bad_lines="skip",
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    if rm.empty:
        return pd.DataFrame(columns=BLAST_COLUMNS)

    rm.columns = REPEATMASKER_OUT_COLUMNS
    query_start = pd.to_numeric(rm["query_start"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    query_end = pd.to_numeric(rm["query_end"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    starts = np.minimum(query_start, query_end)
    ends = np.maximum(query_start, query_end)
    te_bp = (np.abs(query_end - query_start) + 1).astype(np.int64)
    te_bp = np.clip(te_bp, 1, None)
    hit_bp = np.maximum(ends - starts + 1, 1)
    query_left = pd.to_numeric(
        rm["query_left"].astype(str).str.strip("()"), errors="coerce"
    ).fillna(0).to_numpy(dtype=np.int64)
    query_size = np.maximum(query_start, query_end) + np.maximum(query_left, 0)
    # Zero-based field index 9: RepeatMasker's ``matching repeat`` name.
    repeat_name = rm.iloc[:, 9].astype(str).str.strip()
    repeat_class = rm["repeat_class"].astype(str).str.strip()
    repeat_class = repeat_class.replace({"": "Unclassified", "nan": "Unclassified"})
    class_parts = repeat_class.str.split("/", n=1, expand=True)
    order = class_parts[0].fillna("Unclassified").replace("", "Unclassified")
    if 1 in class_parts.columns:
        raw_superfamily = class_parts[1]
        superfamily = raw_superfamily.where(
            raw_superfamily.notna() & raw_superfamily.ne(""),
            order,
        )
    else:
        superfamily = order.copy()
    family = repeat_name.fillna("Unclassified").replace(
        {"": "Unclassified", "nan": "Unclassified"}
    )
    qseqid = repeat_name + "#" + repeat_class
    pident = 100.0 - pd.to_numeric(rm["perc_div"], errors="coerce").fillna(0.0)

    return pd.DataFrame(
        {
            "qseqid": qseqid,
            "sseqid": rm["query_name"].astype(str).str.strip(),
            "pident": pident,
            "length": te_bp,
            "mismatch": 0,
            "gapopen": 0,
            "qstart": 1,
            "qend": te_bp,
            "sstart": starts,
            "send": ends,
            "sstrand": rm["strand"].astype(str).str.strip().replace({"C": "-"}),
            "evalue": 0.0,
            "bitscore": pd.to_numeric(rm["score"], errors="coerce").fillna(0.0),
            # Prepared columns (skip slow qseqid re-parse in prepare_hits).
            "chrom": rm["query_name"].astype(str).str.strip(),
            "start": starts,
            "end": ends,
            # Original RepeatMasker column 10 (matching repeat name).
            "repeat_name": repeat_name,
            "order": order,
            "superfamily": superfamily,
            "family": family,
            "hit_bp": hit_bp,
            "te_bp": te_bp,
            "query_size": query_size,
        }
    )


def read_te_hit_table(path: Path) -> pd.DataFrame:
    """Read BLAST table or RepeatMasker ``.out`` in BLAST column layout."""
    path = path.expanduser().resolve()
    if is_repeatmasker_out_file(path):
        return read_repeatmasker_out(path)
    return read_blast_table(path)


def te_hit_plot_source_label(path: Path) -> str:
    return "RepeatMasker" if is_repeatmasker_out_file(path) else "BLAST"


def _natural_sort_key(name: str) -> list:
    return [int(part) if part.isdigit() else part.lower() for part in re.split(r"(\d+)", name)]


def _format_genome_length_bp(bp: int) -> str:
    """Human-readable chromosome length for map labels."""
    bp = max(0, int(bp))
    if bp >= 1_000_000:
        val = bp / BP_PER_MBP
        return f"{val:.2f} Mb" if val < 100 else f"{val:.1f} Mb"
    if bp >= 10_000:
        return f"{bp / 1000:.0f} kb"
    if bp >= 1000:
        return f"{bp / 1000:.1f} kb"
    return f"{bp:,} bp"


format_genome_length_bp = _format_genome_length_bp


def _chromosome_from_track_label(label: str) -> str:
    """Base chromosome/scaffold name from a composite track label."""
    s = str(label).strip()
    for marker in (" (methylation diff:", " (diff:", " ("):
        idx = s.find(marker)
        if idx > 0:
            return s[:idx]
    return s


def _normalize_coords(start: int, end: int) -> tuple[int, int]:
    s, e = int(start), int(end)
    return (min(s, e), max(s, e))


def invalidate_genome_chromosome_cache(genome_fasta: Path | None = None) -> None:
    """Drop cached chromosome metadata (call when the project genome FASTA changes)."""
    if genome_fasta is None:
        _GENOME_META_CACHE.clear()
        return
    _GENOME_META_CACHE.pop(str(genome_fasta.expanduser().resolve()), None)


def _load_genome_meta_from_fasta(genome_fasta: Path) -> tuple[list[str], dict[str, int]]:
    """Read chromosome names and lengths once; result is cached by path + mtime."""
    path = genome_fasta.expanduser().resolve()
    key = str(path)
    try:
        mtime = path.stat().st_mtime
    except OSError:
        return [], {}
    from repeatviewer.backends.genome_prepare import resolve_genome_chrom_sizes_file

    sizes_path = resolve_genome_chrom_sizes_file(path)
    if sizes_path is not None and sizes_path.is_file():
        try:
            mtime = max(mtime, sizes_path.stat().st_mtime)
        except OSError:
            pass
    else:
        sizes_path = None

    cached = _GENOME_META_CACHE.get(key)
    if cached is not None and cached[0] == mtime:
        return cached[1], cached[2]

    sizes: dict[str, int] = {}
    if sizes_path is not None and sizes_path.is_file():
        try:
            for line in sizes_path.read_text(encoding="utf-8", errors="replace").splitlines():
                if not line.strip():
                    continue
                parts = line.split()
                if len(parts) < 2:
                    continue
                sizes[str(parts[0])] = int(parts[1])
        except Exception as exc:  # noqa: BLE001
            logger.warning("Could not read genome sizes from %s: %s", sizes_path, exc)
            sizes = {}
        if sizes:
            names = sorted(sizes.keys(), key=_natural_sort_key)
            _GENOME_META_CACHE[key] = (mtime, names, sizes)
            return names, sizes

    try:
        for record in SeqIO.parse(str(path), "fasta"):
            sizes[record.id] = len(record.seq)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Could not read genome sizes from %s: %s", path, exc)
        return [], {}

    names = sorted(sizes.keys(), key=_natural_sort_key)
    _GENOME_META_CACHE[key] = (mtime, names, sizes)
    return names, sizes


def _ensure_te_bp(hits: pd.DataFrame) -> pd.DataFrame:
    """TE element length (bp) for divergence weighting; genome span stays in ``hit_bp``."""
    if hits.empty or "te_bp" in hits.columns:
        return hits
    out = hits.copy()
    if "hit_bp" not in out.columns and {"start", "end"}.issubset(out.columns):
        out["hit_bp"] = (out["end"] - out["start"] + 1).clip(lower=1)
    if {"qstart", "qend"}.issubset(out.columns):
        qs = pd.to_numeric(out["qstart"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
        qe = pd.to_numeric(out["qend"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
        out["te_bp"] = np.clip(np.maximum(qs, qe) - np.minimum(qs, qe) + 1, 1, None)
    elif "length" in out.columns:
        out["te_bp"] = pd.to_numeric(out["length"], errors="coerce").fillna(1).clip(lower=1)
    elif "hit_bp" in out.columns:
        out["te_bp"] = out["hit_bp"].clip(lower=1)
    else:
        out["te_bp"] = 1
    return out


REPEATVIEWER_HITS_CACHE_SUFFIX = ".repeatviewer_hits.pkl"
REPEATVIEWER_HITS_CLASSIFICATION_SCHEMA = 2


def is_repeatviewer_hits_cache_file(path: Path) -> bool:
    """True for internal TE hit-table pickle caches (hidden from the file list)."""
    return path.name.endswith(REPEATVIEWER_HITS_CACHE_SUFFIX)


def _hits_cache_path(path: Path) -> Path:
    return path.parent / f"{path.name}{REPEATVIEWER_HITS_CACHE_SUFFIX}"


def _concat_dataframes(frames: list[pd.DataFrame]) -> pd.DataFrame:
    """Concatenate frames, returning empty DataFrame when there is nothing to join."""
    non_empty = [df for df in frames if df is not None and not df.empty]
    if not non_empty:
        return pd.DataFrame()
    if len(non_empty) == 1:
        return non_empty[0]
    return pd.concat(non_empty, ignore_index=True)


def load_combined_te_hits(paths: list[Path]) -> pd.DataFrame:
    """Load and concatenate multiple TE hit tables."""
    frames: list[pd.DataFrame] = []
    for path in paths:
        path = path.expanduser().resolve()
        if not path.is_file():
            continue
        frames.append(load_prepared_te_hits(path))
    if not frames:
        return pd.DataFrame()
    if len(frames) == 1:
        return frames[0]
    return _concat_dataframes(frames)


def _genome_map_output_stem(te_input_files: list[Path]) -> str:
    """Base stem for plot HTML; multi-source is marked in the file suffix instead."""
    if not te_input_files:
        return "genes"
    return te_hit_table_stem(te_input_files[0])


def _annotation_feature_types_for_path(
    path: Path,
    *,
    annotation_feature_types: list[str] | None,
    annotation_feature_types_by_file: dict[str, list[str]] | None,
) -> list[str] | None:
    """Resolve GFF/GTF feature types for one annotation file."""
    if annotation_feature_types_by_file:
        resolved = path.expanduser().resolve()
        key = str(resolved)
        if key in annotation_feature_types_by_file:
            picked = annotation_feature_types_by_file[key]
            return picked if picked else None
        name = resolved.name
        for raw_key, picked in annotation_feature_types_by_file.items():
            if Path(raw_key).expanduser().resolve().name == name:
                return picked if picked else None
    return annotation_feature_types


def gene_hits_from_gff(
    gff_path: Path,
    chrom_sizes: dict[str, int],
    *,
    feature_types: list[str] | None = None,
) -> pd.DataFrame:
    """Convert GFF/GTF feature rows or BED intervals into a hit table for genome-map binning."""
    from repeatviewer.backends.gff_features import _is_vcf_path, parse_annotation

    gff_path = gff_path.expanduser().resolve()
    if not gff_path.is_file():
        return pd.DataFrame()
    is_bed = gff_path.name.lower().endswith(".bed")
    is_vcf = _is_vcf_path(gff_path)
    features = parse_annotation(gff_path)
    if not features:
        return pd.DataFrame()

    if is_bed:
        picked = features
    elif is_vcf:
        want = {str(t).strip().lower() for t in (feature_types or []) if str(t).strip()}
        picked = (
            [f for f in features if f.feature_type.lower() in want]
            if want
            else features
        )
    else:
        want = {str(t).strip().lower() for t in (feature_types or ["gene"]) if str(t).strip()}
        picked = [f for f in features if f.feature_type.lower() in want]
        if not picked:
            return pd.DataFrame()

    genome_names = list(chrom_sizes.keys()) if chrom_sizes else []
    rows: list[dict[str, object]] = []
    for feat in picked:
        chrom: str | None = None
        if genome_names:
            chrom = match_genome_chromosome(feat.seqid, genome_names)
        if chrom is None and not genome_names:
            chrom = str(feat.seqid).strip()
        if chrom is None:
            continue
        ft = feat.feature_type.lower()
        strand_norm = _normalize_hit_strand(feat.strand)
        rows.append(
            {
                "chrom": chrom,
                "start": int(feat.start),
                "end": int(feat.end),
                "pident": 100.0,
                "order": ft,
                "superfamily": ft,
                "family": ft,
                "te_bp": int(feat.length),
                "sstrand": feat.strand,
                "_strand_norm": strand_norm,
            }
        )
    if not rows:
        return pd.DataFrame()
    return _ensure_te_bp(pd.DataFrame(rows))


def load_prepared_te_hits(path: Path) -> pd.DataFrame:
    """Read and normalize TE hits, using a sidecar pickle cache when the source is unchanged."""
    path = path.expanduser().resolve()
    cache = _hits_cache_path(path)
    try:
        src_mtime = path.stat().st_mtime
        if cache.is_file() and cache.stat().st_mtime >= src_mtime:
            cached = pd.read_pickle(cache)
            # Rebuild caches created under older classification rules.
            if (
                {"order", "superfamily", "family"}.issubset(cached.columns)
                and cached.attrs.get("repeatviewer_classification_schema")
                == REPEATVIEWER_HITS_CLASSIFICATION_SCHEMA
            ):
                return cached
    except OSError:
        pass
    hits = prepare_hits(read_te_hit_table(path))
    hits.attrs["repeatviewer_classification_schema"] = REPEATVIEWER_HITS_CLASSIFICATION_SCHEMA
    try:
        hits.to_pickle(cache)
    except OSError as exc:
        logger.debug("TE hits cache write skipped for %s: %s", path, exc)
    return hits


def prepare_hits(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df
    if {"start", "end", "chrom", "pident", "order", "superfamily", "family"}.issubset(df.columns):
        return _ensure_te_bp(df)
    hits = df.copy()
    sstart = pd.to_numeric(hits["sstart"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    send = pd.to_numeric(hits["send"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    hits["start"] = np.minimum(sstart, send)
    hits["end"] = np.maximum(sstart, send)
    hits["pident"] = pd.to_numeric(hits["pident"], errors="coerce").fillna(0.0)
    hits["chrom"] = hits["sseqid"].astype(str)
    hits["hit_bp"] = (hits["end"] - hits["start"] + 1).clip(lower=1)
    qnames = hits["qseqid"].astype(str)
    name_parts = qnames.str.split("#", n=1, expand=True)
    family = name_parts[0].fillna("Unclassified").replace("", "Unclassified")
    after_hash = name_parts[1] if 1 in name_parts.columns else pd.Series("", index=hits.index)
    has_hash = qnames.str.contains("#", regex=False)
    class_parts = after_hash.where(has_hash, "").str.split("/", n=1, expand=True)
    order = class_parts[0].fillna("Unclassified").replace("", "Unclassified")
    if 1 in class_parts.columns:
        superfamily = class_parts[1].fillna("Unclassified").replace("", "Unclassified")
    else:
        superfamily = pd.Series("Unclassified", index=hits.index, dtype=str)
    hits["order"] = order
    hits["superfamily"] = superfamily.fillna("Unclassified").replace("", "Unclassified")
    hits["family"] = family.fillna("Unclassified").replace("", "Unclassified")
    return _ensure_te_bp(hits)


@dataclass(frozen=True)
class _GenomeMapTrack:
    chrom: str
    strand: str | None
    label: str
    kind: str = "te"
    source_key: str | None = None
    diff_pair: tuple[str, str] | None = None
    methylation_subtrack: str | None = None


TE_DIFF_LEGENDGROUP = "__te_diff__"
METH_DIFF_LEGENDGROUP = "__meth_diff__"
TE_DIFF_COLOR_A = DEFAULT_DIFF_COLOR_A
TE_DIFF_COLOR_B = DEFAULT_DIFF_COLOR_B
METH_DIFF_COLOR_A = "#c97b4a"
METH_DIFF_COLOR_B = "#4a8fc9"
METH_COLOR_CG_TOTAL = "#7c9a6e"
METH_COLOR_CG_COVERED = "#3d7ea6"
METH_COLOR_CG_COVERAGE_PCT = "#8b6fad"
METH_COLOR_METH_LEVEL = "#b84c4c"
METH_COVERAGE_PCT_AXIS_MAX = 1.0
METH_LEVEL_AXIS_MAX = 1.0


def _te_diff_tracks_enabled(
    te_sources: list[_TeSourceSpec] | None,
    *,
    show_te_diff_track: bool,
) -> bool:
    """True when a pairwise exclusive-coverage diff track should be drawn."""
    if not show_te_diff_track or not te_sources or len(te_sources) != 2:
        return False
    return not any(spec.split_strand for spec in te_sources)


def _hit_coord_arrays(hits: pd.DataFrame) -> tuple[np.ndarray, np.ndarray]:
    if hits.empty or "start" not in hits.columns or "end" not in hits.columns:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    starts = pd.to_numeric(hits["start"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    ends = pd.to_numeric(hits["end"], errors="coerce").fillna(0).to_numpy(dtype=np.int64)
    lo = np.minimum(starts, ends)
    hi = np.maximum(starts, ends)
    valid = hi >= lo
    return lo[valid], hi[valid]


def _merged_interval_arrays(
    starts: np.ndarray,
    ends: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Merge overlapping/adjacent intervals (sorted single pass)."""
    if starts.size == 0:
        return starts, ends
    order = np.argsort(starts, kind="mergesort")
    starts = starts[order]
    ends = ends[order]
    merged_s: list[int] = []
    merged_e: list[int] = []
    cur_s = int(starts[0])
    cur_e = int(ends[0])
    for start, end in zip(starts[1:], ends[1:], strict=False):
        s, e = int(start), int(end)
        if s <= cur_e + 1:
            cur_e = max(cur_e, e)
        else:
            merged_s.append(cur_s)
            merged_e.append(cur_e)
            cur_s, cur_e = s, e
    merged_s.append(cur_s)
    merged_e.append(cur_e)
    return np.asarray(merged_s, dtype=np.int64), np.asarray(merged_e, dtype=np.int64)


def _subtract_interval_arrays(
    base_st: np.ndarray,
    base_en: np.ndarray,
    sub_st: np.ndarray,
    sub_en: np.ndarray,
) -> tuple[np.ndarray, np.ndarray]:
    """Return merged ``base`` intervals with merged ``sub`` removed."""
    if base_st.size == 0:
        return base_st, base_en
    if sub_st.size == 0:
        return base_st, base_en
    out_s: list[int] = []
    out_e: list[int] = []
    j = 0
    sub_len = int(sub_st.size)
    for i in range(int(base_st.size)):
        seg_s = int(base_st[i])
        seg_e = int(base_en[i])
        cursor = seg_s
        while cursor <= seg_e and j < sub_len:
            bs = int(sub_st[j])
            be = int(sub_en[j])
            if be < cursor:
                j += 1
                continue
            if bs > seg_e:
                break
            if bs > cursor:
                out_s.append(cursor)
                out_e.append(min(bs - 1, seg_e))
            cursor = max(cursor, be + 1)
            if cursor > seg_e:
                break
        if cursor <= seg_e:
            out_s.append(cursor)
            out_e.append(seg_e)
    if not out_s:
        return np.array([], dtype=np.int64), np.array([], dtype=np.int64)
    return np.asarray(out_s, dtype=np.int64), np.asarray(out_e, dtype=np.int64)


def _bin_disjoint_interval_arrays(
    starts: np.ndarray,
    ends: np.ndarray,
    chrom_length: int,
    bin_size: int,
) -> np.ndarray:
    """Sum bp per bin for non-overlapping intervals (vectorized edges + slice fills)."""
    bin_size_i = max(5_000, int(bin_size))
    n_bins = max(1, int(np.ceil(chrom_length / bin_size_i)))
    counts = np.zeros(n_bins, dtype=np.float64)
    if starts.size == 0:
        return counts
    max_pos = max(0, int(chrom_length))
    starts = starts.astype(np.int64, copy=False)
    ends = np.minimum(ends.astype(np.int64, copy=False), max(max_pos - 1, 0))
    valid = ends >= starts
    starts, ends = starts[valid], ends[valid]
    if starts.size == 0:
        return counts

    b0 = starts // bin_size_i
    b1 = np.minimum(ends // bin_size_i, n_bins - 1)
    same = b0 == b1
    if np.any(same):
        np.add.at(counts, b0[same], (ends[same] - starts[same] + 1).astype(np.float64))

    diff = ~same
    if not np.any(diff):
        return counts

    s0 = starts[diff]
    e0 = ends[diff]
    bb0 = b0[diff]
    bb1 = b1[diff]
    first_hi = (bb0 + 1) * bin_size_i - 1
    np.add.at(counts, bb0, (first_hi - s0 + 1).astype(np.float64))
    last_lo = bb1 * bin_size_i
    np.add.at(counts, bb1, (e0 - last_lo + 1).astype(np.float64))
    for bin_lo, bin_hi in zip(bb0, bb1, strict=False):
        bl = int(bin_lo)
        bh = int(bin_hi)
        if bh > bl + 1:
            counts[bl + 1 : bh] += bin_size_i
    return counts


def _compute_te_diff_bins_for_chrom(
    hits_a: pd.DataFrame,
    hits_b: pd.DataFrame,
    chrom_length: int,
    bin_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Exclusive TE coverage (bp) per bin for two merged occupancy maps."""
    a_st, a_en = _merged_interval_arrays(*_hit_coord_arrays(hits_a))
    b_st, b_en = _merged_interval_arrays(*_hit_coord_arrays(hits_b))
    ex_a_st, ex_a_en = _subtract_interval_arrays(a_st, a_en, b_st, b_en)
    ex_b_st, ex_b_en = _subtract_interval_arrays(b_st, b_en, a_st, a_en)
    return (
        _bin_disjoint_interval_arrays(ex_a_st, ex_a_en, chrom_length, bin_size),
        _bin_disjoint_interval_arrays(ex_b_st, ex_b_en, chrom_length, bin_size),
    )


def _te_diff_bins_worker(payload: dict[str, Any]) -> tuple[str, np.ndarray, np.ndarray]:
    """ProcessPool worker: one chromosome of diff bins."""
    starts_a = np.asarray(payload["starts_a"], dtype=np.int64)
    ends_a = np.asarray(payload["ends_a"], dtype=np.int64)
    starts_b = np.asarray(payload["starts_b"], dtype=np.int64)
    ends_b = np.asarray(payload["ends_b"], dtype=np.int64)
    chrom_length = int(payload["chrom_length"])
    bin_size = int(payload["bin_size"])
    label = str(payload["label"])
    bins_a, bins_b = _compute_te_diff_bins_for_chrom_arrays(
        starts_a,
        ends_a,
        starts_b,
        ends_b,
        chrom_length,
        bin_size,
    )
    return label, bins_a, bins_b


def _compute_te_diff_bins_for_chrom_arrays(
    starts_a: np.ndarray,
    ends_a: np.ndarray,
    starts_b: np.ndarray,
    ends_b: np.ndarray,
    chrom_length: int,
    bin_size: int,
) -> tuple[np.ndarray, np.ndarray]:
    a_st, a_en = _merged_interval_arrays(starts_a, ends_a)
    b_st, b_en = _merged_interval_arrays(starts_b, ends_b)
    ex_a_st, ex_a_en = _subtract_interval_arrays(a_st, a_en, b_st, b_en)
    ex_b_st, ex_b_en = _subtract_interval_arrays(b_st, b_en, a_st, a_en)
    return (
        _bin_disjoint_interval_arrays(ex_a_st, ex_a_en, chrom_length, bin_size),
        _bin_disjoint_interval_arrays(ex_b_st, ex_b_en, chrom_length, bin_size),
    )


def _compute_all_track_diff_bins(
    *,
    chrom_order: list[str],
    plot_chrom_sizes: dict[str, int],
    te_hits_by_chrom_source: dict[str, dict[str, pd.DataFrame]],
    te_sources: list[_TeSourceSpec],
    bin_size: int,
    num_workers: int | None,
    cancel: JobCancellation | None = None,
) -> dict[str, tuple[np.ndarray, np.ndarray, str, str]]:
    check_cancelled(cancel)
    key_a, key_b = te_sources[0].key, te_sources[1].key
    label_a, label_b = te_sources[0].label, te_sources[1].label
    empty = pd.DataFrame()
    payloads: list[dict[str, Any]] = []
    for chrom in chrom_order:
        hits_a = te_hits_by_chrom_source.get(key_a, {}).get(chrom, empty)
        hits_b = te_hits_by_chrom_source.get(key_b, {}).get(chrom, empty)
        sa, ea = _hit_coord_arrays(hits_a)
        sb, eb = _hit_coord_arrays(hits_b)
        payloads.append(
            {
                "label": _diff_track_label(chrom, label_a, label_b),
                "chrom_length": plot_chrom_sizes[chrom],
                "bin_size": bin_size,
                "starts_a": sa,
                "ends_a": ea,
                "starts_b": sb,
                "ends_b": eb,
                "label_a": label_a,
                "label_b": label_b,
            }
        )

    out: dict[str, tuple[np.ndarray, np.ndarray, str, str]] = {}
    workers = _effective_genome_map_workers(
        num_workers,
        num_chromosomes=len(payloads),
        num_hits=sum(int(p["starts_a"].size + p["starts_b"].size) for p in payloads),
    )
    t0 = time.perf_counter()
    if workers <= 1 or len(payloads) < 2:
        for payload in payloads:
            check_cancelled(cancel)
            label, bins_a, bins_b = _te_diff_bins_worker(payload)
            out[label] = (bins_a, bins_b, label_a, label_b)
    else:
        logger.info("TE diff track: comparing %d chromosomes with %d workers", len(payloads), workers)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_te_diff_bins_worker, payload) for payload in payloads]
            try:
                for future in as_completed(futures):
                    check_cancelled(cancel)
                    label, bins_a, bins_b = future.result()
                    out[label] = (bins_a, bins_b, label_a, label_b)
            except JobCancelledError:
                for future in futures:
                    future.cancel()
                pool.shutdown(wait=False, cancel_futures=True)
                raise
    logger.info("TE diff track: finished in %.1fs", time.perf_counter() - t0)
    return out


def _global_te_diff_max(
    track_diff_bins: dict[str, tuple[np.ndarray, np.ndarray, str, str]],
) -> float:
    """Peak exclusive TE bp in any diff window across all chromosomes."""
    peak = 1.0
    for bins_a, bins_b, _, _ in track_diff_bins.values():
        peak = max(
            peak,
            float(bins_a.max(initial=0)),
            float(bins_b.max(initial=0)),
        )
    return peak


def _shared_te_diff_scale_max(
    track_diff_bins: dict[str, tuple[np.ndarray, np.ndarray, str, str]],
) -> tuple[float, float]:
    """Shared diff axis/bar scale (95th pct.) and absolute peak for metadata."""
    values: list[float] = []
    for bins_a, bins_b, _, _ in track_diff_bins.values():
        for arr in (bins_a, bins_b):
            if arr.size:
                positive = arr[arr > 0]
                if positive.size:
                    values.extend(positive.astype(float).tolist())
    if not values:
        return 1.0, 1.0
    peak = float(max(values))
    p95 = float(np.percentile(values, 95))
    scale_max = max(p95, peak * 0.2, 1000.0)
    return scale_max, peak


def _diff_track_label(chrom: str, label_a: str, label_b: str) -> str:
    return f"{chrom} (diff: {label_a} vs {label_b})"


def _meth_diff_track_label(chrom: str, label_a: str, label_b: str) -> str:
    short_a = _compact_bismark_cov_label(label_a, max_len=20)
    short_b = _compact_bismark_cov_label(label_b, max_len=20)
    return f"{chrom} (methylation diff: {short_a} vs {short_b})"


def _meth_diff_tracks_enabled(
    methylation_sources: list[_MethylationSourceSpec] | None,
    *,
    show_meth_diff_track: bool,
) -> bool:
    if not show_meth_diff_track or not methylation_sources or len(methylation_sources) != 2:
        return False
    return True


def _compute_methylation_diff_bins(
    *,
    chrom_order: list[str],
    plot_chrom_sizes: dict[str, int],
    methylation_hits_by_chrom_source: dict[str, dict[str, MethylationChromCov]],
    methylation_sources: list[_MethylationSourceSpec],
    bin_size: int,
    genome_fasta: Path | None,
    region: tuple[float, float] | None,
    cg_by_chrom: dict[str, np.ndarray] | None = None,
    num_workers: int | None = None,
    cancel: JobCancellation | None = None,
) -> dict[str, tuple[np.ndarray, np.ndarray, str, str]]:
    check_cancelled(cancel)
    key_a, key_b = methylation_sources[0].key, methylation_sources[1].key
    label_a, label_b = methylation_sources[0].label, methylation_sources[1].label
    path_a = methylation_sources[0].path
    path_b = methylation_sources[1].path
    cov_cache_a = (
        str(_methylation_cov_cache_path(path_a))
        if path_a and _methylation_cov_cache_path(path_a).is_file()
        else ""
    )
    cov_cache_b = (
        str(_methylation_cov_cache_path(path_b))
        if path_b and _methylation_cov_cache_path(path_b).is_file()
        else ""
    )
    use_process_pool = bool(cov_cache_a and cov_cache_b)
    payloads: list[dict[str, Any]] = []
    num_hits = 0
    cpus = os.cpu_count() or 4
    array_workers = max(1, min(8, int(num_workers or cpus), cpus))
    for chrom in chrom_order:
        cov_a = (
            {}
            if use_process_pool
            else methylation_cov_arrays(
                methylation_hits_by_chrom_source.get(key_a, {}).get(chrom),
            )
        )
        cov_b = (
            {}
            if use_process_pool
            else methylation_cov_arrays(
                methylation_hits_by_chrom_source.get(key_b, {}).get(chrom),
            )
        )
        num_hits += methylation_cov_site_count(cov_a) + methylation_cov_site_count(cov_b)
        payload: dict[str, Any] = {
            "label": _meth_diff_track_label(chrom, label_a, label_b),
            "chrom": chrom,
            "chrom_length": plot_chrom_sizes[chrom],
            "bin_size": bin_size,
            "cg_total": (cg_by_chrom or {}).get(
                chrom,
                np.zeros(
                    max(1, int(np.ceil(plot_chrom_sizes[chrom] / max(5_000, int(bin_size))))),
                    dtype=np.float64,
                ),
            ),
            "region": region,
            "label_a": label_a,
            "label_b": label_b,
            "cov_a": cov_a,
            "cov_b": cov_b,
            "array_workers": array_workers,
        }
        if use_process_pool:
            payload["source_cov_cache_paths"] = {"a": cov_cache_a, "b": cov_cache_b}
        payloads.append(payload)

    out: dict[str, tuple[np.ndarray, np.ndarray, str, str]] = {}
    workers = _effective_methylation_workers(
        num_workers,
        num_chromosomes=len(payloads),
        num_hits=num_hits,
    )
    t0 = time.perf_counter()
    if workers <= 1 or len(payloads) < 2:
        for payload in payloads:
            check_cancelled(cancel)
            label, diff_a, diff_b, src_a, src_b = methylation_diff_bins_worker(payload)
            out[label] = (diff_a, diff_b, src_a, src_b)
    else:
        logger.info(
            "Methylation diff track: comparing %d chromosomes with %d workers",
            len(payloads),
            workers,
        )
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(methylation_diff_bins_worker, payload) for payload in payloads]
            try:
                for future in as_completed(futures):
                    check_cancelled(cancel)
                    label, diff_a, diff_b, src_a, src_b = future.result()
                    out[label] = (diff_a, diff_b, src_a, src_b)
            except JobCancelledError:
                for future in futures:
                    future.cancel()
                pool.shutdown(wait=False, cancel_futures=True)
                raise
    logger.info("Methylation diff track: finished in %.1fs", time.perf_counter() - t0)
    return out


def _shared_meth_diff_scale_max(
    track_meth_diff_bins: dict[str, tuple[np.ndarray, np.ndarray, str, str]],
) -> tuple[float, float]:
    values: list[float] = []
    for bins_a, bins_b, _, _ in track_meth_diff_bins.values():
        for arr in (bins_a, bins_b):
            if arr.size:
                positive = arr[arr > 0]
                if positive.size:
                    values.extend(positive.astype(float).tolist())
    if not values:
        return 0.1, 0.1
    peak = float(max(values))
    p95 = float(np.percentile(values, 95))
    scale_max = min(1.0, max(p95, peak * 0.2, 0.01))
    return scale_max, peak


@dataclass(frozen=True)
class _TeSourceSpec:
    key: str
    label: str
    split_strand: bool = False
    input_kind: str = TE_INPUT_KIND_MASKING


@dataclass(frozen=True)
class _GeneSourceSpec:
    key: str
    label: str
    split_strand: bool = False


@dataclass(frozen=True)
class _MethylationSourceSpec:
    key: str
    label: str
    path: Path | None = None


def te_source_spec_for_path(
    path: Path,
    index: int,
    *,
    multi_file: bool,
    split_strand: bool = False,
    label: str | None = None,
    input_kind: str = TE_INPUT_KIND_MASKING,
) -> _TeSourceSpec:
    """Unique source key/label for one TE hit table in a multi-file genome map."""
    display = label or (path.name if multi_file else te_hit_table_stem(path))
    return _TeSourceSpec(
        key=f"te_src_{index}",
        label=display,
        split_strand=split_strand,
        input_kind=input_kind,
    )


def _te_input_kind_for_source(
    te_sources: list[_TeSourceSpec] | None,
    source_key: str | None,
) -> str:
    if not te_sources or not source_key:
        return TE_INPUT_KIND_MASKING
    for spec in te_sources:
        if spec.key == source_key:
            return spec.input_kind
    return TE_INPUT_KIND_MASKING


def _density_max_for_te_tracks(
    tracks: list["_GenomeMapTrack"],
    track_density_max: dict[str, float],
    te_sources: list[_TeSourceSpec] | None,
    *,
    use_category_density: bool,
    track_cat_bins: dict[tuple[str, str], tuple],
    categories_for_legend: list[str],
) -> tuple[float, float]:
    """Peak density for masking vs consensus TE track groups."""
    masking_max = 1.0
    consensus_max = 1.0
    for track in tracks:
        if track.kind != "te":
            continue
        input_kind = _te_input_kind_for_source(te_sources, track.source_key)
        vmax = float(track_density_max.get(track.label, 1.0))
        if use_category_density and categories_for_legend:
            sample = track_cat_bins.get((track.label, categories_for_legend[0]))
            if sample is not None:
                totals = np.zeros(len(sample[1]), dtype=float)
                for cat in categories_for_legend:
                    cat_bins = track_cat_bins.get((track.label, cat))
                    if cat_bins is None:
                        continue
                    totals += cat_bins[1]
                if totals.size:
                    vmax = max(vmax, float(totals.max()))
        if input_kind == TE_INPUT_KIND_CONSENSUS:
            consensus_max = max(consensus_max, vmax)
        else:
            masking_max = max(masking_max, vmax)
    return masking_max, consensus_max


def _track_input_kind_map(
    tracks: list["_GenomeMapTrack"],
    te_sources: list[_TeSourceSpec] | None,
) -> dict[str, str]:
    out: dict[str, str] = {}
    for track in tracks:
        if track.kind == "gene":
            out[track.label] = "gene"
        elif track.kind == "methylation":
            out[track.label] = "methylation"
        elif track.kind == "te":
            out[track.label] = _te_input_kind_for_source(te_sources, track.source_key)
    return out


def _bin_methylation_tracks_for_genome_map(
    tracks: list[_GenomeMapTrack],
    methylation_hits_by_chrom_source: dict[str, dict[str, MethylationChromCov]],
    chrom_sizes: dict[str, int],
    *,
    bin_size: int,
    coverage_min: int,
    region: tuple[float, float] | None,
    genome_fasta: Path | None = None,
    cg_by_chrom: dict[str, np.ndarray] | None = None,
    methylation_sources: list[_MethylationSourceSpec] | None = None,
    num_workers: int | None = None,
    cancel: JobCancellation | None = None,
) -> dict[str, tuple]:
    check_cancelled(cancel)
    chrom_track_specs: dict[str, list[dict[str, str]]] = {}
    for track in tracks:
        if track.kind != "methylation" or not track.source_key or not track.methylation_subtrack:
            continue
        chrom_track_specs.setdefault(track.chrom, []).append(
            {
                "track_label": track.label,
                "source_key": track.source_key,
                "subtrack": track.methylation_subtrack,
            }
        )
    if not chrom_track_specs:
        return {}

    num_hits = 0
    for chrom, specs in chrom_track_specs.items():
        source_keys = {spec["source_key"] for spec in specs}
        for source_key in source_keys:
            cov = methylation_hits_by_chrom_source.get(source_key, {}).get(chrom)
            num_hits += methylation_cov_site_count(cov)

    workers = _effective_methylation_workers(
        num_workers,
        num_chromosomes=len(chrom_track_specs),
        num_hits=num_hits,
    )
    cpus = os.cpu_count() or 4
    array_workers = max(1, min(8, int(num_workers or cpus), cpus))
    source_paths: dict[str, Path] = {}
    if methylation_sources:
        for spec in methylation_sources:
            if spec.path and spec.path.is_file():
                source_paths[spec.key] = spec.path
    source_cov_cache_paths = {
        key: str(_methylation_cov_cache_path(path))
        for key, path in source_paths.items()
        if _methylation_cov_cache_path(path).is_file()
    }
    source_bins_cache_paths = {
        key: str(_methylation_bins_cache_path(path, bin_size))
        for key, path in source_paths.items()
        if _methylation_bins_cache_path(path, bin_size).is_file()
    }
    use_process_pool = workers > 1 and bool(source_cov_cache_paths)
    payloads: list[dict[str, Any]] = []
    for chrom, track_specs in chrom_track_specs.items():
        source_keys = {spec["source_key"] for spec in track_specs}
        source_cov = (
            {}
            if use_process_pool
            else {
                source_key: methylation_cov_arrays(
                    methylation_hits_by_chrom_source.get(source_key, {}).get(chrom),
                )
                for source_key in source_keys
            }
        )
        payloads.append(
            {
                "chrom": chrom,
                "chrom_length": chrom_sizes[chrom],
                "bin_size": bin_size,
                "coverage_min": max(1, int(coverage_min)),
                "region": region,
                "array_workers": array_workers,
                "cg_total": (cg_by_chrom or {}).get(
                    chrom,
                    np.zeros(
                        max(1, int(np.ceil(chrom_sizes[chrom] / max(5_000, int(bin_size))))),
                        dtype=np.float64,
                    ),
                ),
                "source_cov": source_cov,
                "source_cov_cache_paths": {
                    key: source_cov_cache_paths[key]
                    for key in source_keys
                    if key in source_cov_cache_paths
                },
                "source_bins_cache_paths": {
                    key: source_bins_cache_paths[key]
                    for key in source_keys
                    if key in source_bins_cache_paths
                },
                "track_specs": track_specs,
            }
        )

    track_bins: dict[str, tuple] = {}
    binned_by_source_chrom: dict[str, dict[str, Any]] = {}
    t0 = time.perf_counter()

    def _consume_worker_result(
        chrom_name: str,
        chrom_bins_map: dict[str, tuple],
        binned: dict[str, MethylationChromBins],
    ) -> None:
        track_bins.update(chrom_bins_map)
        for source_key, chrom_bins in binned.items():
            binned_by_source_chrom.setdefault(source_key, {})[chrom_name] = chrom_bins

    if workers <= 1:
        for payload in payloads:
            check_cancelled(cancel)
            chrom_name, chrom_bins_map, binned = bin_methylation_chromosome_worker(payload)
            _consume_worker_result(chrom_name, chrom_bins_map, binned)
    elif use_process_pool:
        logger.info(
            "Methylation genome map: binning %d chromosomes with %d processes",
            len(payloads),
            workers,
        )
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(bin_methylation_chromosome_worker, payload) for payload in payloads]
            try:
                for future in as_completed(futures):
                    check_cancelled(cancel)
                    chrom_name, chrom_bins_map, binned = future.result()
                    _consume_worker_result(chrom_name, chrom_bins_map, binned)
            except JobCancelledError:
                for future in futures:
                    future.cancel()
                pool.shutdown(wait=False, cancel_futures=True)
                raise
    else:
        logger.info(
            "Methylation genome map: binning %d chromosomes with %d threads",
            len(payloads),
            workers,
        )
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(bin_methylation_chromosome_worker, payload) for payload in payloads]
            try:
                for future in as_completed(futures):
                    check_cancelled(cancel)
                    chrom_name, chrom_bins_map, binned = future.result()
                    _consume_worker_result(chrom_name, chrom_bins_map, binned)
            except JobCancelledError:
                for future in futures:
                    future.cancel()
                raise
    for source_key, chrom_bins_map in binned_by_source_chrom.items():
        cov_path = source_paths.get(source_key)
        if not cov_path or not cov_path.is_file() or not chrom_bins_map:
            continue
        try:
            source_mtime = cov_path.stat().st_mtime
        except OSError:
            source_mtime = 0.0
        cache_path = _methylation_bins_cache_path(cov_path, bin_size)
        merged: dict[str, MethylationChromBins] = {}
        if cache_path.is_file():
            for chrom_name in chrom_track_specs:
                cached = _load_methylation_bins_chrom_cache(
                    cache_path,
                    chrom_name,
                    bin_size=bin_size,
                    cg_total=(cg_by_chrom or {}).get(
                        chrom_name,
                        np.zeros(1, dtype=np.float64),
                    ),
                )
                if cached is not None:
                    merged[chrom_name] = cached
        merged.update(chrom_bins_map)
        _save_methylation_bins_cache(
            cache_path,
            merged,
            bin_size=bin_size,
            source_mtime=source_mtime,
        )
    logger.info("Methylation genome map: binning finished in %.1fs", time.perf_counter() - t0)
    return track_bins


def gene_source_spec_for_path(
    path: Path,
    index: int,
    *,
    multi_file: bool,
    split_strand: bool = False,
    label: str | None = None,
) -> _GeneSourceSpec:
    display = label or (path.name if multi_file else path.stem)
    return _GeneSourceSpec(key=f"gene_src_{index}", label=display, split_strand=split_strand)


def methylation_source_spec_for_path(
    path: Path,
    index: int,
    *,
    multi_file: bool,
    label: str | None = None,
) -> _MethylationSourceSpec:
    raw = label or (path.name if multi_file else path.stem)
    display = _compact_bismark_cov_label(raw, max_len=28 if multi_file else 24)
    return _MethylationSourceSpec(
        key=f"meth_src_{index}",
        label=display,
        path=path.expanduser().resolve(),
    )


def _track_display_label(
    chrom: str,
    src: _TeSourceSpec,
    strand: str | None,
    *,
    multi_te: bool,
) -> str:
    if strand == "+":
        if multi_te:
            return f"{chrom} (forward, {src.label})"
        return f"{chrom} (forward)"
    if strand == "-":
        if multi_te:
            return f"{chrom} (reverse, {src.label})"
        return f"{chrom} (reverse)"
    if multi_te:
        return f"{chrom} ({src.label})"
    return chrom


def _normalize_hit_strand(raw: object) -> str:
    s = str(raw).strip().lower()
    if s in {"-", "minus", "c", "reverse", "rev", "complement"}:
        return "-"
    return "+"


def _methylation_track_display_label(
    chrom: str,
    src: _MethylationSourceSpec,
    subtrack: str,
    *,
    multi_meth: bool,
) -> str:
    suffix = {
        METH_SUBTRACK_CG_TOTAL: "Total C+G number",
        METH_SUBTRACK_CG_COVERED: "Covered C+G number",
        METH_SUBTRACK_CG_COVERAGE_PCT: "Covered fraction",
        METH_SUBTRACK_METH_LEVEL: "Methylation level",
    }.get(subtrack, subtrack)
    if multi_meth:
        return f"{chrom} ({suffix}, {src.label})"
    return f"{chrom} ({suffix})"


def _track_group_key_for_track(track: _GenomeMapTrack) -> str | None:
    if track.kind == "te_diff":
        return "te_diff"
    if track.kind == "meth_diff":
        return "meth_diff"
    if track.source_key:
        return track.source_key
    return None


def _panel_track_groups_for_map(
    *,
    te_sources: list[_TeSourceSpec],
    show_te_diff: bool,
    resolved_gene_sources: list[_GeneSourceSpec],
    include_gene_tracks: bool,
    resolved_methylation_sources: list[_MethylationSourceSpec],
    include_methylation_tracks: bool,
    show_meth_diff: bool,
) -> list[dict[str, str]]:
    """Logical track groups shown in the genome-map panel (one row per input / diff track)."""
    groups: list[dict[str, str]] = []
    for src in te_sources:
        groups.append({"key": src.key, "kind": "te", "label": src.label})
    if show_te_diff and len(te_sources) == 2:
        groups.append({"key": "te_diff", "kind": "te_diff", "label": "TE difference"})
    if include_gene_tracks:
        for src in resolved_gene_sources:
            groups.append({"key": src.key, "kind": "gene", "label": src.label})
    if include_methylation_tracks:
        for src in resolved_methylation_sources:
            groups.append({"key": src.key, "kind": "methylation", "label": src.label})
    if show_meth_diff and len(resolved_methylation_sources) == 2:
        groups.append(
            {"key": "meth_diff", "kind": "meth_diff", "label": "Methylation difference"}
        )
    return groups


def _enrich_panel_track_groups_with_default_ymax(
    panel_groups: list[dict[str, object]],
    *,
    tracks: list[_GenomeMapTrack],
    track_density_norm: dict[str, float],
    te_sources: list[_TeSourceSpec],
    masking_density_max: float,
    consensus_density_max: float,
    gene_density_max: float,
    methylation_norm_max: dict[str, float],
    te_diff_density_max: float | None = None,
    meth_diff_density_max: float | None = None,
) -> list[dict[str, object]]:
    """Attach per-group default Y-max values used to pre-fill the tracks panel."""
    enriched: list[dict[str, object]] = []
    for group in panel_groups:
        entry = dict(group)
        key = str(group.get("key") or "")
        kind = str(group.get("kind") or "")
        if kind == "te":
            vmax = 0.0
            for track in tracks:
                if track.kind == "te" and track.source_key == key:
                    vmax = max(vmax, float(track_density_norm.get(track.label, 0.0)))
            if vmax <= 0:
                input_kind = _te_input_kind_for_source(te_sources, key)
                vmax = (
                    consensus_density_max
                    if input_kind == TE_INPUT_KIND_CONSENSUS
                    else masking_density_max
                )
            entry["defaultYMax"] = int(max(1, round(vmax)))
        elif kind == "gene":
            vmax = max(
                (
                    float(track_density_norm.get(track.label, 0.0))
                    for track in tracks
                    if track.kind == "gene" and track.source_key == key
                ),
                default=0.0,
            )
            if vmax <= 0:
                vmax = gene_density_max
            entry["defaultYMax"] = int(max(1, round(vmax)))
        elif kind == "methylation":
            entry["defaultYMax"] = {
                "cg_coverage_pct": round(
                    float(methylation_norm_max.get(METH_SUBTRACK_CG_COVERAGE_PCT, 1.0)),
                    4,
                ),
                "meth_level": round(
                    float(methylation_norm_max.get(METH_SUBTRACK_METH_LEVEL, 1.0)),
                    4,
                ),
            }
        elif kind == "te_diff" and te_diff_density_max is not None:
            entry["defaultYMax"] = int(max(1, round(float(te_diff_density_max))))
        elif kind == "meth_diff" and meth_diff_density_max is not None:
            entry["defaultYMax"] = round(float(meth_diff_density_max), 4)
        enriched.append(entry)
    return enriched


def _gene_track_display_label(
    chrom: str,
    src: _GeneSourceSpec,
    strand: str | None,
    *,
    multi_gene: bool,
) -> str:
    if strand == "+":
        if multi_gene:
            return f"{chrom} (genes, forward, {src.label})"
        return f"{chrom} (genes, forward)"
    if strand == "-":
        if multi_gene:
            return f"{chrom} (genes, reverse, {src.label})"
        return f"{chrom} (genes, reverse)"
    if multi_gene:
        return f"{chrom} (genes, {src.label})"
    return f"{chrom} (genes)"


def _ensure_unique_track_labels(tracks: list[_GenomeMapTrack]) -> list[_GenomeMapTrack]:
    """Guarantee unique track labels (required for y-pos and binning dict keys)."""
    seen: dict[str, int] = {}
    unique: list[_GenomeMapTrack] = []
    for track in tracks:
        label = track.label
        if label not in seen:
            seen[label] = 1
            unique.append(track)
            continue
        seen[label] += 1
        unique.append(replace(track, label=f"{label} [{seen[label]}]"))
    return unique


def _genome_map_tracks(
    chrom_order: list[str],
    *,
    te_sources: list[_TeSourceSpec],
    gene_sources: list[_GeneSourceSpec] | None = None,
    include_gene_tracks: bool = False,
    methylation_sources: list[_MethylationSourceSpec] | None = None,
    include_methylation_tracks: bool = False,
    split_strand_gene_tracks: bool = False,
    show_te_diff: bool = False,
    show_meth_diff: bool = False,
) -> list[_GenomeMapTrack]:
    """Build tracks: per TE input file on each chromosome (top = first file), optional diff, then genes."""
    tracks: list[_GenomeMapTrack] = []
    multi_te = len(te_sources) > 1
    gene_sources = list(gene_sources or [])
    methylation_sources = list(methylation_sources or [])
    multi_gene = len(gene_sources) > 1
    multi_meth = len(methylation_sources) > 1
    diff_pair: tuple[str, str] | None = None
    meth_diff_pair: tuple[str, str] | None = None
    if show_te_diff and len(te_sources) == 2:
        diff_pair = (te_sources[0].key, te_sources[1].key)
    if show_meth_diff and len(methylation_sources) == 2:
        meth_diff_pair = (methylation_sources[0].key, methylation_sources[1].key)
    for chrom in chrom_order:
        for src in te_sources:
            if src.split_strand:
                tracks.append(
                    _GenomeMapTrack(
                        chrom,
                        "+",
                        _track_display_label(chrom, src, "+", multi_te=multi_te),
                        "te",
                        src.key,
                    )
                )
                tracks.append(
                    _GenomeMapTrack(
                        chrom,
                        "-",
                        _track_display_label(chrom, src, "-", multi_te=multi_te),
                        "te",
                        src.key,
                    )
                )
            else:
                tracks.append(
                    _GenomeMapTrack(
                        chrom,
                        None,
                        _track_display_label(chrom, src, None, multi_te=multi_te),
                        "te",
                        src.key,
                    )
                )
        if diff_pair is not None:
            tracks.append(
                _GenomeMapTrack(
                    chrom,
                    None,
                    _diff_track_label(chrom, te_sources[0].label, te_sources[1].label),
                    "te_diff",
                    diff_pair=diff_pair,
                )
            )
        if include_gene_tracks:
            sources = gene_sources or [
                _GeneSourceSpec(key="gene_src_0", label="genes", split_strand=split_strand_gene_tracks)
            ]
            for src in sources:
                split_gene = split_strand_gene_tracks or src.split_strand
                if split_gene:
                    tracks.append(
                        _GenomeMapTrack(
                            chrom,
                            "+",
                            _gene_track_display_label(chrom, src, "+", multi_gene=multi_gene),
                            "gene",
                            src.key,
                        )
                    )
                    tracks.append(
                        _GenomeMapTrack(
                            chrom,
                            "-",
                            _gene_track_display_label(chrom, src, "-", multi_gene=multi_gene),
                            "gene",
                            src.key,
                        )
                    )
                else:
                    tracks.append(
                        _GenomeMapTrack(
                            chrom,
                            None,
                            _gene_track_display_label(chrom, src, None, multi_gene=multi_gene),
                            "gene",
                            src.key,
                        )
                    )
        if include_methylation_tracks:
            sources = methylation_sources or [
                _MethylationSourceSpec(key="meth_src_0", label="methylation")
            ]
            for src in sources:
                for subtrack in (
                    METH_SUBTRACK_CG_TOTAL,
                    METH_SUBTRACK_CG_COVERED,
                    METH_SUBTRACK_CG_COVERAGE_PCT,
                    METH_SUBTRACK_METH_LEVEL,
                ):
                    tracks.append(
                        _GenomeMapTrack(
                            chrom,
                            None,
                            _methylation_track_display_label(
                                chrom,
                                src,
                                subtrack,
                                multi_meth=multi_meth,
                            ),
                            "methylation",
                            src.key,
                            methylation_subtrack=subtrack,
                        )
                    )
        if meth_diff_pair is not None:
            tracks.append(
                _GenomeMapTrack(
                    chrom,
                    None,
                    _meth_diff_track_label(
                        chrom,
                        methylation_sources[0].label,
                        methylation_sources[1].label,
                    ),
                    "meth_diff",
                    diff_pair=meth_diff_pair,
                )
            )
    return tracks


def _filter_hits_for_track(hits: pd.DataFrame, track: _GenomeMapTrack) -> pd.DataFrame:
    if hits.empty:
        return hits
    if track.strand is not None:
        if "_strand_norm" in hits.columns:
            return hits[hits["_strand_norm"] == track.strand]
        if "sstrand" in hits.columns:
            strands = hits["sstrand"].map(_normalize_hit_strand)
            return hits[strands == track.strand]
        return hits.iloc[0:0]
    if track.kind in {"gene", "methylation"}:
        return hits
    return hits


def _track_strand_caption(track: _GenomeMapTrack) -> str:
    if track.kind == "methylation":
        sub = track.methylation_subtrack or ""
        if sub == METH_SUBTRACK_CG_TOTAL:
            return "total in window (genome)"
        if sub == METH_SUBTRACK_CG_COVERED:
            return "covered lines in .cov per window"
        if sub == METH_SUBTRACK_CG_COVERAGE_PCT:
            return "covered fraction (track 2 / track 1)"
        if sub == METH_SUBTRACK_METH_LEVEL:
            return "mean methylation level (0–1)"
        return "methylation"
    if track.kind == "gene":
        if track.strand == "+":
            return "genes, forward strand"
        if track.strand == "-":
            return "genes, reverse strand"
        return "genes"
    if track.strand == "+":
        return "forward strand"
    if track.strand == "-":
        return "reverse strand"
    return ""


def chromosome_names_from_genome(genome_fasta: Path | None) -> list[str]:
    """Chromosome IDs from genome FASTA (natural sort: chr1, chr2, …, chr10)."""
    if genome_fasta is None:
        return []
    genome_fasta = genome_fasta.expanduser().resolve()
    if not genome_fasta.is_file():
        return []
    names, _ = _load_genome_meta_from_fasta(genome_fasta)
    return names


def genome_chrom_sizes_from_fasta(genome_fasta: Path | None) -> dict[str, int]:
    """Chromosome lengths (bp) from genome FASTA or sidecar ``*.chrom.sizes``."""
    if genome_fasta is None:
        return {}
    genome_fasta = genome_fasta.expanduser().resolve()
    if not genome_fasta.is_file():
        return {}
    _, sizes = _load_genome_meta_from_fasta(genome_fasta)
    return dict(sizes)


def superfamily_names_from_blast(blast_txt: Path) -> list[str]:
    """Unique TE superfamilies parsed from qseqid format: after the class '/'."""
    blast_txt = blast_txt.expanduser().resolve()
    if not blast_txt.is_file():
        return []
    try:
        hits = load_prepared_te_hits(blast_txt)
    except Exception:
        return []
    if hits.empty or "superfamily" not in hits.columns:
        return []
    return _ordered_te_types({str(v) for v in hits["superfamily"].astype(str) if str(v)})


def chromosome_sizes_from_genome(
    genome_fasta: Path | None,
    hits: pd.DataFrame,
) -> dict[str, int]:
    sizes: dict[str, int] = {}
    if genome_fasta is not None:
        genome_fasta = genome_fasta.expanduser().resolve()
        if genome_fasta.is_file():
            _, sizes = _load_genome_meta_from_fasta(genome_fasta)
            sizes = dict(sizes)

    if not hits.empty and "chrom" in hits.columns and "end" in hits.columns:
        genome_keys = set(sizes.keys())
        ends = hits.groupby("chrom", sort=False)["end"].max()
        for chrom, max_end in ends.items():
            chrom_s = str(chrom)
            if genome_keys and chrom_s not in genome_keys:
                continue
            sizes[chrom_s] = max(sizes.get(chrom_s, 0), int(max_end))

    return sizes


def _ordered_chromosomes(
    chrom_sizes: dict[str, int],
    _hits: pd.DataFrame,
    *,
    chromosome_sort: str = DEFAULT_CHROMOSOME_SORT,
) -> list[str]:
    """Chromosomes in display order (natural or by size), capped at MAX_CHROMOSOMES."""
    keys = list(chrom_sizes.keys())
    if chromosome_sort == CHROMOSOME_SORT_SIZE_DESC:
        ordered = sorted(keys, key=lambda c: (-chrom_sizes[c], _natural_sort_key(c)))
    elif chromosome_sort == CHROMOSOME_SORT_SIZE_ASC:
        ordered = sorted(keys, key=lambda c: (chrom_sizes[c], _natural_sort_key(c)))
    else:
        ordered = sorted(keys, key=_natural_sort_key)
    return ordered[:MAX_CHROMOSOMES]


_CHROM_BP_SUFFIX_RE = re.compile(r"_bp\d+$", re.IGNORECASE)


def _chrom_base_name(name: str) -> str:
    return _CHROM_BP_SUFFIX_RE.sub("", str(name).strip())


def match_genome_chromosome(hit_chrom: str, genome_names: list[str]) -> str | None:
    """Map a hit / RepeatMasker chromosome ID to a genome FASTA name when uniquely possible."""
    hc = str(hit_chrom).strip()
    if not hc or not genome_names:
        return None
    if hc in genome_names:
        return hc
    prefix_hits = [
        g for g in genome_names if g.startswith(hc + "_") or g.startswith(hc + ".")
    ]
    if len(prefix_hits) == 1:
        return prefix_hits[0]
    prefix_genome = [
        g for g in genome_names if hc.startswith(g + "_") or hc.startswith(g + ".")
    ]
    if len(prefix_genome) == 1:
        return prefix_genome[0]
    hc_lower = hc.lower()
    suffix_matches = [
        g
        for g in genome_names
        if g.lower() == hc_lower
        or g.lower().endswith("." + hc_lower)
        or g.lower().endswith("_" + hc_lower)
    ]
    if len(suffix_matches) == 1:
        return suffix_matches[0]
    hc_base = _chrom_base_name(hc)
    if hc_base:
        base_matches = [g for g in genome_names if _chrom_base_name(g) == hc_base]
        if len(base_matches) == 1:
            return base_matches[0]
    return None


def harmonize_hit_chromosomes(hits: pd.DataFrame, genome_sizes: dict[str, int]) -> pd.DataFrame:
    """Rewrite hit ``chrom`` values to canonical genome FASTA IDs when names differ slightly."""
    if hits.empty or not genome_sizes or "chrom" not in hits.columns:
        return hits
    genome_names = list(genome_sizes.keys())
    mapping: dict[str, str] = {}
    for chrom in hits["chrom"].astype(str).unique():
        canon = match_genome_chromosome(chrom, genome_names)
        if canon and canon != chrom:
            mapping[chrom] = canon
    if not mapping:
        return hits
    out = hits.copy()
    out["chrom"] = out["chrom"].astype(str).replace(mapping)
    return out


def _map_chromosome_filter_to_genome(
    chromosomes_only: list[str],
    genome_names: list[str],
) -> set[str]:
    mapped: set[str] = set()
    for name in chromosomes_only:
        n = str(name).strip()
        if not n:
            continue
        if n in genome_names:
            mapped.add(n)
            continue
        canon = match_genome_chromosome(n, genome_names)
        if canon:
            mapped.add(canon)
    return mapped


def _resolve_chrom_order(
    chrom_sizes: dict[str, int],
    hits: pd.DataFrame,
    chromosomes_only: list[str] | None,
    *,
    chromosome_sort: str = DEFAULT_CHROMOSOME_SORT,
) -> tuple[list[str], str | None]:
    chrom_order = _ordered_chromosomes(
        chrom_sizes,
        hits,
        chromosome_sort=chromosome_sort,
    )
    if not chromosomes_only:
        return chrom_order, None

    genome_names = list(chrom_sizes.keys())
    allowed = _map_chromosome_filter_to_genome(chromosomes_only, genome_names)
    if not allowed:
        hit_chroms = set(hits["chrom"].astype(str)) if not hits.empty else set()
        with_hits = [c for c in chrom_order if c in hit_chroms]
        if with_hits:
            return with_hits, (
                "Chromosome filter does not match this genome — showing chromosomes with TE hits."
            )
        return chrom_order, (
            "Chromosome filter does not match this genome — showing all genome chromosomes."
        )

    filtered = [c for c in chrom_order if c in allowed]
    if filtered:
        return filtered, None

    hit_chroms = set(hits["chrom"].astype(str)) if not hits.empty else set()
    with_hits = [c for c in chrom_order if c in hit_chroms]
    if with_hits:
        return with_hits, (
            "Chromosome filter matched no tracks — showing chromosomes with TE hits."
        )
    return chrom_order, (
        "Chromosome filter matched no tracks — showing all genome chromosomes."
    )


def _chromosomes_with_plottable_data(
    hits: pd.DataFrame,
    gene_df: pd.DataFrame,
    methylation_df: pd.DataFrame,
    methylation_cov_by_chrom_source: dict[str, dict[str, Any]] | None,
) -> set[str]:
    chroms: set[str] = set()
    for df in (hits, gene_df, methylation_df):
        if not df.empty and "chrom" in df.columns:
            chroms |= set(df["chrom"].astype(str))
    if methylation_cov_by_chrom_source:
        for per_source in methylation_cov_by_chrom_source.values():
            chroms |= {str(k) for k in per_source.keys()}
    return chroms


def _prune_chromosomes_without_data(
    chrom_order: list[str],
    chroms_with_data: set[str],
) -> tuple[list[str], list[str]]:
    if not chrom_order or not chroms_with_data:
        return chrom_order, []
    kept = [c for c in chrom_order if c in chroms_with_data]
    dropped = [c for c in chrom_order if c not in chroms_with_data]
    return kept, dropped


def _empty_chrom_order_hint(
    chrom_sizes: dict[str, int],
    hits: pd.DataFrame,
    chromosomes_only: list[str] | None,
) -> str:
    hints: list[str] = []
    if hits.empty:
        hints.append("input has no TE hits")
    elif chrom_sizes:
        genome_keys = set(chrom_sizes.keys())
        hit_chroms = set(hits["chrom"].astype(str).unique())
        if genome_keys and hit_chroms and not (genome_keys & hit_chroms):
            sample_h = ", ".join(sorted(hit_chroms, key=_natural_sort_key)[:3])
            sample_g = ", ".join(sorted(genome_keys, key=_natural_sort_key)[:3])
            hints.append(
                f"TE file uses {sample_h}… but genome FASTA has {sample_g}… "
                "(RepeatMasker output must be from the same genome)"
            )
    if chromosomes_only and chrom_sizes:
        allowed = _map_chromosome_filter_to_genome(chromosomes_only, list(chrom_sizes.keys()))
        if not allowed:
            hints.append("saved chromosome filter does not match this genome")
    if not hints:
        return ""
    return "<br>".join(hints)


def _nothing_to_plot_message(
    *,
    te_paths: list[Path],
    gene_paths: list[Path],
    methylation_paths: list[Path],
    chromosomes_only: list[str] | None,
    chrom_sizes: dict[str, int],
) -> str:
    lines = ["There is nothing to plot for the current inputs."]
    if te_paths and not gene_paths and not methylation_paths:
        lines.append("The TE input file(s) contain no genome hits.")
    elif not te_paths and gene_paths and not methylation_paths:
        lines.append("The annotation file(s) contain no features mapped to this genome.")
    elif not te_paths and not gene_paths and methylation_paths:
        lines.append("The methylation file(s) contain no sites mapped to this genome.")
    else:
        lines.append("No TE hits, gene features, or methylation sites were loaded.")
    if chromosomes_only and chrom_sizes:
        allowed = _map_chromosome_filter_to_genome(chromosomes_only, list(chrom_sizes.keys()))
        if not allowed:
            lines.append(
                "The selected chromosome name(s) do not match the genome FASTA. "
                "Clear the chromosome filter or pick names from the list."
            )
        else:
            names = ", ".join(sorted(allowed, key=_natural_sort_key)[:4])
            if len(allowed) > 4:
                names += ", …"
            lines.append(
                f"No data on the selected chromosome(s) ({names}). "
                "Try another chromosome, run Blast on the TE consensus first, "
                "or confirm the genome FASTA matches the hit table."
            )
    elif te_paths:
        lines.append(
            "If you use a TE consensus FASTA, run Blast against the genome first. "
            "RepeatMasker or masking output must be from the same genome as the project FASTA."
        )
    return "\n\n".join(lines)


def _pident_colorscale(z: np.ndarray, zmin: float, zmax: float) -> list[str]:
    if z.size == 0:
        return []
    span = max(zmax - zmin, 1.0)
    colors = []
    for value in z:
        t = (float(value) - zmin) / span
        t = max(0.0, min(1.0, t))
        r = int(230 - 90 * t)
        g = int(90 + 140 * t)
        b = int(60 + 120 * t)
        colors.append(f"rgb({r},{g},{b})")
    return colors


def _bin_chromosome(
    chrom_hits: pd.DataFrame,
    chrom_length: int,
    bin_size: int,
    *,
    density_mode: str = TE_GENOME_MAP_DENSITY_HIT_COUNT,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    bin_size = max(5_000, int(bin_size))
    n_bins = max(1, int(np.ceil(chrom_length / bin_size)))
    counts = np.zeros(n_bins, dtype=float)
    pident_sum = np.zeros(n_bins, dtype=float)
    bp_cover = np.zeros(n_bins, dtype=float)

    if chrom_hits.empty:
        pass
    else:
        starts = chrom_hits["start"].to_numpy(dtype=np.int64)
        ends = chrom_hits["end"].to_numpy(dtype=np.int64)
        pidents = chrom_hits["pident"].to_numpy(dtype=np.float64)
        bin_size_i = int(bin_size)
        for start, end, pident in zip(starts, ends, pidents):
            if end < start:
                start, end = end, start
            b0 = int(start // bin_size_i)
            b1 = int(min(end // bin_size_i, n_bins - 1))
            for bin_idx in range(b0, b1 + 1):
                bin_lo = bin_idx * bin_size_i
                bin_hi = min((bin_idx + 1) * bin_size_i - 1, int(end))
                seg_lo = max(int(start), bin_lo)
                seg_hi = min(int(end), bin_hi)
                overlap = seg_hi - seg_lo + 1
                if overlap <= 0:
                    continue
                counts[bin_idx] += 1.0
                pident_sum[bin_idx] += pident * overlap
                bp_cover[bin_idx] += overlap

    centers_mbp = (np.arange(n_bins) + 0.5) * bin_size / BP_PER_MBP
    mean_pident = np.divide(
        pident_sum,
        bp_cover,
        out=np.full(n_bins, np.nan, dtype=np.float64),
        where=bp_cover > 0,
    )
    if density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH:
        density = bp_cover
    else:
        density = counts
    return centers_mbp, density, mean_pident, counts, bp_cover


def _bin_chromosome_by_category(
    chrom_hits: pd.DataFrame,
    chrom_length: int,
    bin_size: int,
    categories: list[str],
    cat_col: str,
    *,
    density_mode: str = TE_GENOME_MAP_DENSITY_HIT_COUNT,
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    """Bin all TE categories on one track in a single pass over hits."""
    bin_size_i = max(5_000, int(bin_size))
    n_bins = max(1, int(np.ceil(chrom_length / bin_size_i)))
    cat_to_idx = {str(cat): idx for idx, cat in enumerate(categories)}
    n_cats = len(categories)
    counts = np.zeros((n_cats, n_bins), dtype=np.float64)
    pident_sum = np.zeros((n_cats, n_bins), dtype=np.float64)
    bp_cover = np.zeros((n_cats, n_bins), dtype=np.float64)

    if not chrom_hits.empty and cat_col in chrom_hits.columns:
        starts = chrom_hits["start"].to_numpy(dtype=np.int64)
        ends = chrom_hits["end"].to_numpy(dtype=np.int64)
        pidents = chrom_hits["pident"].to_numpy(dtype=np.float64)
        cat_vals = chrom_hits[cat_col].astype(str).to_numpy()
        for start, end, pident, cat in zip(starts, ends, pidents, cat_vals):
            ci = cat_to_idx.get(cat)
            if ci is None:
                continue
            if end < start:
                start, end = end, start
            b0 = int(start // bin_size_i)
            b1 = int(min(end // bin_size_i, n_bins - 1))
            for bin_idx in range(b0, b1 + 1):
                bin_lo = bin_idx * bin_size_i
                bin_hi = min((bin_idx + 1) * bin_size_i - 1, int(end))
                seg_lo = max(int(start), bin_lo)
                seg_hi = min(int(end), bin_hi)
                overlap = seg_hi - seg_lo + 1
                if overlap <= 0:
                    continue
                counts[ci, bin_idx] += 1.0
                pident_sum[ci, bin_idx] += pident * overlap
                bp_cover[ci, bin_idx] += overlap

    centers_mbp = (np.arange(n_bins) + 0.5) * bin_size_i / BP_PER_MBP
    out: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    for idx, cat in enumerate(categories):
        mean_pident = np.divide(
            pident_sum[idx],
            bp_cover[idx],
            out=np.full(n_bins, np.nan, dtype=np.float64),
            where=bp_cover[idx] > 0,
        )
        density = bp_cover[idx] if density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH else counts[idx]
        out[cat] = (centers_mbp, density, mean_pident, counts[idx], bp_cover[idx])
    return out


def _effective_genome_map_workers(
    requested: int | None,
    *,
    num_chromosomes: int,
    num_hits: int,
) -> int:
    cpus = os.cpu_count() or 4
    cap = max(1, min(int(requested or cpus), cpus))
    if num_hits < _PAR_BIN_MIN_HITS or num_chromosomes < _PAR_BIN_MIN_CHROMS:
        return 1
    return max(1, min(cap, num_chromosomes))


def _effective_methylation_workers(
    requested: int | None,
    *,
    num_chromosomes: int,
    num_hits: int,
) -> int:
    """Worker count for methylation binning (lower thresholds than TE binning)."""
    cpus = os.cpu_count() or 4
    cap = max(1, min(int(requested or cpus), cpus))
    if num_hits < _PAR_METH_MIN_HITS or num_chromosomes < _PAR_METH_MIN_CHROMS:
        return 1
    return max(1, min(cap, num_chromosomes))


def _bin_arrays_single_track(
    starts: np.ndarray,
    ends: np.ndarray,
    pidents: np.ndarray,
    chrom_length: int,
    bin_size: int,
    *,
    density_mode: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    bin_size_i = max(5_000, int(bin_size))
    n_bins = max(1, int(np.ceil(chrom_length / bin_size_i)))
    counts = np.zeros(n_bins, dtype=np.float64)
    pident_sum = np.zeros(n_bins, dtype=np.float64)
    bp_cover = np.zeros(n_bins, dtype=np.float64)
    if starts.size:
        _accumulate_hit_bins(
            starts,
            ends,
            pidents,
            counts,
            pident_sum,
            bp_cover,
            bin_size_i=bin_size_i,
            n_bins=n_bins,
        )
    centers_mbp = (np.arange(n_bins) + 0.5) * bin_size_i / BP_PER_MBP
    mean_pident = np.divide(
        pident_sum,
        bp_cover,
        out=np.full(n_bins, np.nan, dtype=np.float64),
        where=bp_cover > 0,
    )
    density = bp_cover if density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH else counts
    return centers_mbp, density, mean_pident, counts, bp_cover


def _bin_arrays_by_category(
    starts: np.ndarray,
    ends: np.ndarray,
    pidents: np.ndarray,
    cat_values: np.ndarray,
    chrom_length: int,
    bin_size: int,
    categories: list[str],
    *,
    density_mode: str,
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    bin_size_i = max(5_000, int(bin_size))
    n_bins = max(1, int(np.ceil(chrom_length / bin_size_i)))
    n_cats = len(categories)
    cat_to_idx = {str(cat): idx for idx, cat in enumerate(categories)}
    counts = np.zeros((n_cats, n_bins), dtype=np.float64)
    pident_sum = np.zeros((n_cats, n_bins), dtype=np.float64)
    bp_cover = np.zeros((n_cats, n_bins), dtype=np.float64)
    if starts.size:
        cat_indices = np.fromiter(
            (cat_to_idx.get(str(c), -1) for c in cat_values),
            dtype=np.int16,
            count=cat_values.size,
        )
        _accumulate_hit_bins(
            starts,
            ends,
            pidents,
            counts,
            pident_sum,
            bp_cover,
            bin_size_i=bin_size_i,
            n_bins=n_bins,
            cat_indices=cat_indices,
        )
    centers_mbp = (np.arange(n_bins) + 0.5) * bin_size_i / BP_PER_MBP
    out: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    for idx, cat in enumerate(categories):
        mean_pident = np.divide(
            pident_sum[idx],
            bp_cover[idx],
            out=np.full(n_bins, np.nan, dtype=np.float64),
            where=bp_cover[idx] > 0,
        )
        density = bp_cover[idx] if density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH else counts[idx]
        out[cat] = (centers_mbp, density, mean_pident, counts[idx], bp_cover[idx])
    return out


def _accumulate_hit_bins(
    starts: np.ndarray,
    ends: np.ndarray,
    pidents: np.ndarray,
    counts: np.ndarray,
    pident_sum: np.ndarray,
    bp_cover: np.ndarray,
    *,
    bin_size_i: int,
    n_bins: int,
    cat_indices: np.ndarray | None = None,
) -> None:
    """Accumulate per-bin hit statistics (optionally per TE category)."""
    multi_cat = cat_indices is not None and counts.ndim == 2
    for i in range(starts.size):
        start = int(starts[i])
        end = int(ends[i])
        if end < start:
            start, end = end, start
        if multi_cat:
            ci = int(cat_indices[i])
            if ci < 0:
                continue
            row_counts = counts[ci]
            row_pident = pident_sum[ci]
            row_bp = bp_cover[ci]
        else:
            row_counts = counts
            row_pident = pident_sum
            row_bp = bp_cover
        pident = float(pidents[i])
        b0 = start // bin_size_i
        b1 = min(end // bin_size_i, n_bins - 1)
        for bin_idx in range(b0, b1 + 1):
            bin_lo = bin_idx * bin_size_i
            seg_lo = max(start, bin_lo)
            seg_hi = min(end, min((bin_idx + 1) * bin_size_i - 1, end))
            overlap = seg_hi - seg_lo + 1
            if overlap <= 0:
                continue
            row_counts[bin_idx] += 1.0
            row_pident[bin_idx] += pident * overlap
            row_bp[bin_idx] += overlap


def _maybe_filter_bins_region(
    binned: tuple[np.ndarray, ...],
    bin_size: int,
    region: tuple[float, float] | None,
) -> tuple[np.ndarray, ...]:
    if region is None:
        return binned
    return _filter_bins_to_mbp_window(*binned, bin_size, region[0], region[1])


def _track_hit_arrays(
    track_hits: pd.DataFrame,
    cat_col: str | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray | None]:
    if track_hits.empty:
        empty_i = np.array([], dtype=np.int64)
        empty_f = np.array([], dtype=np.float64)
        return empty_i, empty_i, empty_f, None
    starts = track_hits["start"].to_numpy(dtype=np.int64, copy=False)
    ends = track_hits["end"].to_numpy(dtype=np.int64, copy=False)
    pidents = track_hits["pident"].to_numpy(dtype=np.float64, copy=False)
    cats = None
    if cat_col and cat_col in track_hits.columns:
        cats = track_hits[cat_col].astype(str).to_numpy()
    return starts, ends, pidents, cats


def _bin_all_tracks_for_genome_map(
    tracks: list[_GenomeMapTrack],
    hits_by_chrom: dict[str, pd.DataFrame],
    chrom_sizes: dict[str, int],
    *,
    bin_size: int,
    density_mode: str,
    use_category_density: bool,
    categories_for_legend: list[str],
    cat_col: str | None,
    region: tuple[float, float] | None,
    num_workers: int | None,
    num_hits: int,
    collapse_overlaps: bool = False,
    annotation_density_mode: str = TE_GENOME_MAP_DENSITY_TE_LENGTH,
    gene_hits_by_chrom: dict[str, pd.DataFrame] | None = None,
    te_hits_by_chrom_source: dict[str, dict[str, pd.DataFrame]] | None = None,
    gene_hits_by_chrom_source: dict[str, dict[str, pd.DataFrame]] | None = None,
    cancel: JobCancellation | None = None,
) -> tuple[dict[tuple[str, str], tuple], dict[str, tuple], float]:
    """Bin TE density per track, using multiple processes when worthwhile."""
    check_cancelled(cancel)
    track_cat_bins: dict[tuple[str, str], tuple] = {}
    track_bins: dict[str, tuple] = {}
    global_max_count = 1.0
    gene_hits_by_chrom = gene_hits_by_chrom or {}
    te_hits_by_chrom_source = te_hits_by_chrom_source or {}
    gene_hits_by_chrom_source = gene_hits_by_chrom_source or {}
    chrom_to_tracks: dict[str, list[_GenomeMapTrack]] = {}
    for track in tracks:
        chrom_to_tracks.setdefault(track.chrom, []).append(track)

    workers = _effective_genome_map_workers(
        num_workers,
        num_chromosomes=len(chrom_to_tracks),
        num_hits=num_hits,
    )
    payloads: list[dict[str, Any]] = []
    empty_hits = pd.DataFrame()
    for chrom, chrom_tracks in chrom_to_tracks.items():
        chrom_te_hits = hits_by_chrom.get(chrom, empty_hits)
        track_specs: list[dict[str, Any]] = []
        for track in chrom_tracks:
            if track.kind == "te_diff":
                continue
            if track.kind == "methylation":
                continue
            if track.kind == "gene":
                if track.source_key and track.source_key in gene_hits_by_chrom_source:
                    source_chrom = gene_hits_by_chrom_source[track.source_key]
                    track_hits = _filter_hits_for_track(
                        source_chrom.get(chrom, empty_hits),
                        track,
                    )
                else:
                    track_hits = _filter_hits_for_track(
                        gene_hits_by_chrom.get(chrom, empty_hits),
                        track,
                    )
                track_use_category = False
            elif track.source_key:
                source_chrom = te_hits_by_chrom_source.get(track.source_key, {})
                track_hits = _filter_hits_for_track(
                    source_chrom.get(chrom, empty_hits),
                    track,
                )
                track_use_category = use_category_density
            else:
                track_hits = _filter_hits_for_track(chrom_te_hits, track)
                track_use_category = use_category_density
            starts, ends, pidents, cats = _track_hit_arrays(
                track_hits,
                cat_col if track_use_category else None,
            )
            spec: dict[str, Any] = {
                "track_label": track.label,
                "starts": starts,
                "ends": ends,
                "pidents": pidents,
                "use_category": track_use_category,
            }
            if track_use_category and cats is not None:
                spec["cats"] = cats
            if track.kind == "gene":
                spec["density_mode"] = annotation_density_mode
            track_specs.append(spec)
        payloads.append(
            {
                "chrom": chrom,
                "chrom_length": chrom_sizes[chrom],
                "bin_size": bin_size,
                "density_mode": density_mode,
                "collapse_overlaps": collapse_overlaps,
                "use_category": use_category_density,
                "categories": tuple(categories_for_legend),
                "region": region,
                "tracks": track_specs,
            }
        )

    if workers <= 1:
        for payload in payloads:
            check_cancelled(cancel)
            result = bin_chromosome_tracks_worker(payload)
            track_cat_bins.update(result["track_cat_bins"])
            track_bins.update(result["track_bins"])
            global_max_count = max(global_max_count, float(result["global_max"]))
    else:
        logger.info("TE genome map: binning %d chromosomes with %d workers", len(payloads), workers)
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(bin_chromosome_tracks_worker, payload) for payload in payloads]
            try:
                for future in as_completed(futures):
                    check_cancelled(cancel)
                    result = future.result()
                    track_cat_bins.update(result["track_cat_bins"])
                    track_bins.update(result["track_bins"])
                    global_max_count = max(global_max_count, float(result["global_max"]))
            except JobCancelledError:
                for future in futures:
                    future.cancel()
                pool.shutdown(wait=False, cancel_futures=True)
                raise

    return track_cat_bins, track_bins, global_max_count


def _filter_bins_to_mbp_window(
    centers_mbp: np.ndarray,
    density: np.ndarray,
    mean_pident: np.ndarray,
    hit_counts: np.ndarray,
    te_bp: np.ndarray,
    bin_size: int,
    lo_mbp: float,
    hi_mbp: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    half_w = (bin_size / BP_PER_MBP) / 2.0
    keep = (centers_mbp + half_w >= lo_mbp) & (centers_mbp - half_w <= hi_mbp)
    return (
        centers_mbp[keep],
        density[keep],
        mean_pident[keep],
        hit_counts[keep],
        te_bp[keep],
    )


def _category_series(hits: pd.DataFrame, color_by: str) -> pd.Series:
    if color_by == TE_GENOME_MAP_COLOR_ORDER and "order" in hits.columns:
        return hits["order"].astype(str)
    if color_by == TE_GENOME_MAP_COLOR_FAMILY and "family" in hits.columns:
        return hits["family"].astype(str)
    if color_by == TE_GENOME_MAP_COLOR_SUPERFAMILY and "superfamily" in hits.columns:
        return hits["superfamily"].astype(str)
    return pd.Series(["% identity"] * len(hits), index=hits.index, dtype=str)


def _is_default_te_type(te_type: str) -> bool:
    name = str(te_type).strip()
    if not name:
        return False
    return name in _TE_TYPE_TO_CATEGORY or name.lower() in _TE_TYPE_TO_CATEGORY


def _te_type_parent_category(te_type: str) -> str:
    name = str(te_type).strip()
    if not name:
        return "Unknown"
    if name in _TE_TYPE_TO_CATEGORY:
        return _TE_TYPE_TO_CATEGORY[name]
    lowered = name.lower()
    if lowered in _TE_TYPE_TO_CATEGORY:
        return _TE_TYPE_TO_CATEGORY[lowered]
    return "Unknown"


def _te_type_display_color(te_type: str) -> str:
    parent = _te_type_parent_category(te_type)
    return TE_CATEGORY_COLORS.get(parent, TE_CATEGORY_COLORS["Unknown"])


def _extra_te_type_color(te_type: str, extra_index: int) -> str:
    if TE_EXTRA_TYPE_COLORS:
        return TE_EXTRA_TYPE_COLORS[extra_index % len(TE_EXTRA_TYPE_COLORS)]
    # Stable hash fallback if palette is empty.
    digest = sum(ord(ch) for ch in str(te_type)) + extra_index * 17
    hue = digest % 360
    return f"hsl({hue}, 55%, 52%)"


def _ordered_te_types(types: list[str] | set[str]) -> list[str]:
    unique = sorted({str(t) for t in types if str(t)}, key=_natural_sort_key)

    def _sort_key(te_type: str) -> tuple[int, tuple[Any, ...]]:
        if _is_default_te_type(te_type):
            parent = _te_type_parent_category(te_type)
            return (_TE_CATEGORY_RANK.get(parent, len(TE_CATEGORY_ORDER)), _natural_sort_key(te_type))
        return (len(TE_CATEGORY_ORDER) + 1, _natural_sort_key(te_type))

    return sorted(unique, key=_sort_key)


def _family_level_palette() -> list[str]:
    seen: set[str] = set()
    palette: list[str] = []
    for color in (*TE_CATEGORY_COLORS.values(), *TE_EXTRA_TYPE_COLORS):
        if color not in seen:
            seen.add(color)
            palette.append(color)
    return palette


def _family_level_color_map(categories: list[str]) -> dict[str, str]:
    """Assign a distinct colour to each TE family/subfamily element."""
    palette = _family_level_palette()
    colors: dict[str, str] = {}
    for index, cat in enumerate(categories):
        if index < len(palette):
            colors[cat] = palette[index]
        else:
            hue = (index * 37) % 360
            colors[cat] = f"hsl({hue}, 58%, 50%)"
    return colors


def _ordered_te_families(hits: pd.DataFrame) -> list[str]:
    """Stable family order grouped by superfamily taxonomy, then family name."""
    if hits.empty or "family" not in hits.columns:
        return []
    family_rank: dict[str, tuple[int, tuple[Any, ...], tuple[Any, ...]]] = {}
    has_superfamily = "superfamily" in hits.columns
    for fam_raw, sf_raw in zip(
        hits["family"].astype(str),
        hits["superfamily"].astype(str) if has_superfamily else [""] * len(hits),
    ):
        fam = fam_raw.strip()
        if not fam:
            continue
        if has_superfamily and sf_raw.strip():
            sf = sf_raw.strip()
            rank = _TE_CATEGORY_RANK.get(
                _te_type_parent_category(sf), len(TE_CATEGORY_ORDER)
            )
            sf_key = _natural_sort_key(sf)
        else:
            rank = len(TE_CATEGORY_ORDER)
            sf_key = ()
        key = (rank, sf_key, _natural_sort_key(fam))
        if fam not in family_rank or key < family_rank[fam]:
            family_rank[fam] = key
    return sorted(family_rank.keys(), key=lambda name: family_rank[name])


def _legend_categories(hits: pd.DataFrame, color_by: str) -> list[str]:
    if hits.empty:
        return []
    if color_by == TE_GENOME_MAP_COLOR_FAMILY and "family" in hits.columns:
        return _ordered_te_families(hits)
    cat_col = (
        "order"
        if color_by == TE_GENOME_MAP_COLOR_ORDER
        else "superfamily"
        if color_by == TE_GENOME_MAP_COLOR_SUPERFAMILY
        else "family"
    )
    if cat_col not in hits.columns:
        return []
    return _ordered_te_types({str(v) for v in hits[cat_col].astype(str) if str(v)})


def _category_color_map(
    categories: list[str],
    *,
    color_by: str = TE_GENOME_MAP_COLOR_SUPERFAMILY,
) -> dict[str, str]:
    if color_by == TE_GENOME_MAP_COLOR_FAMILY:
        return _family_level_color_map(categories)
    colors: dict[str, str] = {}
    extra_index = 0
    for cat in categories:
        if _is_default_te_type(cat):
            colors[cat] = _te_type_display_color(cat)
        else:
            colors[cat] = _extra_te_type_color(cat, extra_index)
            extra_index += 1
    return colors


def _prepare_te_hits_for_plot(
    hits: pd.DataFrame,
    chrom_sizes: dict[str, int],
    *,
    split_strand: bool = False,
    superfamilies_only: list[str] | None = None,
) -> pd.DataFrame:
    if hits.empty:
        return hits
    if "start" not in hits.columns:
        hits = prepare_hits(hits)
    hits = harmonize_hit_chromosomes(hits, chrom_sizes)
    if (
        split_strand
        and "sstrand" in hits.columns
        and "_strand_norm" not in hits.columns
    ):
        hits = hits.assign(_strand_norm=hits["sstrand"].map(_normalize_hit_strand))
    if superfamilies_only is not None and "superfamily" in hits.columns:
        allowed_sf = {str(v) for v in superfamilies_only}
        hits = hits[hits["superfamily"].astype(str).isin(allowed_sf)]
    return hits


def _per_track_density_max(
    tracks: list[_GenomeMapTrack],
    track_bins: dict[str, tuple],
    track_cat_bins: dict[tuple[str, str], tuple],
    *,
    use_category_density: bool,
    categories_for_legend: list[str],
) -> dict[str, float]:
    """Maximum binned density per track label (TE and gene scales computed separately)."""
    out: dict[str, float] = {}
    for track in tracks:
        vmax = 1.0
        if use_category_density and categories_for_legend and track.kind not in {
            "gene",
            "te_diff",
            "methylation",
        }:
            sample = track_cat_bins.get((track.label, categories_for_legend[0]))
            if sample is not None:
                totals = np.zeros(len(sample[1]), dtype=float)
                for cat in categories_for_legend:
                    cat_bins = track_cat_bins.get((track.label, cat))
                    if cat_bins is None:
                        continue
                    totals += cat_bins[1]
                if totals.size:
                    vmax = max(vmax, float(totals.max()))
        else:
            binned = track_bins.get(track.label)
            if binned is not None and binned[1].size:
                vmax = max(vmax, float(np.max(binned[1])))
        out[track.label] = vmax
    return out


def build_blast_genome_map_figure(
    hits: pd.DataFrame,
    chrom_sizes: dict[str, int],
    *,
    bin_size: int = DEFAULT_BIN_SIZE,
    density_mode: str = TE_GENOME_MAP_DENSITY_HIT_COUNT,
    color_by: str = TE_GENOME_MAP_COLOR_SUPERFAMILY,
    query_label: str = "",
    region: tuple[str, float, float] | None = None,
    chromosomes_only: list[str] | None = None,
    superfamilies_only: list[str] | None = None,
    chromosome_sort: str = DEFAULT_CHROMOSOME_SORT,
    split_strand_tracks: bool = False,
    collapse_overlaps: bool = False,
    show_te_diff_track: bool = False,
    show_meth_diff_track: bool = False,
    show_window_ref_line: bool = True,
    annotation_density_mode: str = TE_GENOME_MAP_DENSITY_TE_LENGTH,
    gene_hits: pd.DataFrame | None = None,
    gene_sources: list[_GeneSourceSpec] | None = None,
    gene_hits_by_source: dict[str, pd.DataFrame] | None = None,
    methylation_hits: pd.DataFrame | None = None,
    methylation_sources: list[_MethylationSourceSpec] | None = None,
    methylation_hits_by_source: dict[str, pd.DataFrame] | None = None,
    methylation_cov_by_chrom_source: dict[str, dict[str, MethylationChromCov]] | None = None,
    methylation_coverage_min: int = 5,
    genome_fasta: Path | None = None,
    te_sources: list[_TeSourceSpec] | None = None,
    te_hits_by_source: dict[str, pd.DataFrame] | None = None,
    plot_source: str = "BLAST",
    num_workers: int | None = None,
    cancel: JobCancellation | None = None,
    progress: Callable[[int, str], None] | None = None,
) -> go.Figure:
    """
    region: optional (chrom_name, start_mbp, end_mbp) to zoom one chromosome and finer context.
    """
    def _emit(message: str) -> None:
        if progress is not None:
            progress(3, message)

    check_cancelled(cancel)
    prepared_by_source: dict[str, pd.DataFrame] = {}
    if te_hits_by_source and te_sources:
        for spec in te_sources:
            prepared_by_source[spec.key] = _prepare_te_hits_for_plot(
                te_hits_by_source.get(spec.key, pd.DataFrame()),
                chrom_sizes,
                split_strand=spec.split_strand,
                superfamilies_only=superfamilies_only,
            )
        if prepared_by_source:
            hits = _concat_dataframes(
                [df for df in prepared_by_source.values() if not df.empty],
            )
        else:
            hits = pd.DataFrame()
    else:
        hits = _prepare_te_hits_for_plot(
            hits,
            chrom_sizes,
            split_strand=split_strand_tracks,
            superfamilies_only=superfamilies_only,
        )
        if te_sources is None:
            te_sources = [
                _TeSourceSpec(
                    key="__single__",
                    label="",
                    split_strand=split_strand_tracks,
                )
            ]
            prepared_by_source = {"__single__": hits}
        elif te_sources:
            for spec in te_sources:
                prepared_by_source[spec.key] = hits
        else:
            hits = pd.DataFrame()
            prepared_by_source = {}

    gene_df = gene_hits.copy() if gene_hits is not None else pd.DataFrame()
    methylation_df = (
        methylation_hits.copy() if methylation_hits is not None else pd.DataFrame()
    )
    prepared_gene_by_source: dict[str, pd.DataFrame] = {}
    prepared_methylation_by_source: dict[str, pd.DataFrame] = {}
    prepared_methylation_cov_by_source: dict[str, dict[str, MethylationChromCov]] = {}
    if methylation_cov_by_chrom_source:
        prepared_methylation_cov_by_source = {
            key: dict(src) for key, src in methylation_cov_by_chrom_source.items()
        }
    if gene_hits_by_source:
        for key, src_df in gene_hits_by_source.items():
            prepared_gene_by_source[key] = src_df.copy() if not src_df.empty else pd.DataFrame()
        if not gene_df.empty:
            gene_df = _concat_dataframes(
                [df for df in prepared_gene_by_source.values() if not df.empty],
            )
        else:
            gene_df = _concat_dataframes(
                [df for df in prepared_gene_by_source.values() if not df.empty],
            )
    if not gene_df.empty:
        gene_df = harmonize_hit_chromosomes(gene_df, chrom_sizes)
        if prepared_gene_by_source and chrom_sizes:
            allowed = set(chrom_sizes.keys())
            for key in list(prepared_gene_by_source.keys()):
                df = harmonize_hit_chromosomes(prepared_gene_by_source[key], chrom_sizes)
                if not df.empty and "chrom" in df.columns:
                    df = df[df["chrom"].astype(str).isin(allowed)]
                prepared_gene_by_source[key] = df
    if methylation_hits_by_source:
        for key, src_df in methylation_hits_by_source.items():
            prepared_methylation_by_source[key] = (
                src_df.copy() if not src_df.empty else pd.DataFrame()
            )
        if not methylation_df.empty:
            methylation_df = _concat_dataframes(
                [df for df in prepared_methylation_by_source.values() if not df.empty],
            )
        else:
            methylation_df = _concat_dataframes(
                [df for df in prepared_methylation_by_source.values() if not df.empty],
            )

    if not chrom_sizes:
        fig = go.Figure()
        fig.update_layout(
            title="Genome BLAST map — no genome/chromosome data",
            template="plotly_white",
        )
        return fig

    x_range: tuple[float, float] | None = None
    plot_chrom_sizes = chrom_sizes
    if region is not None:
        chrom_r, lo_m, hi_m = region[0], float(region[1]), float(region[2])
        if chrom_r not in chrom_sizes:
            fig = go.Figure()
            fig.update_layout(
                title=f"Region plot — chromosome {chrom_r!r} not in genome FASTA",
                template="plotly_white",
            )
            return fig
        if not hits.empty:
            hits = hits[hits["chrom"] == chrom_r]
            lo_bp = int(lo_m * BP_PER_MBP)
            hi_bp = int(hi_m * BP_PER_MBP)
            hits = hits[(hits["end"] >= lo_bp) & (hits["start"] <= hi_bp)]
            for key in list(prepared_by_source.keys()):
                src_df = prepared_by_source[key]
                if src_df.empty:
                    continue
                src_df = src_df[src_df["chrom"] == chrom_r]
                src_df = src_df[(src_df["end"] >= lo_bp) & (src_df["start"] <= hi_bp)]
                prepared_by_source[key] = src_df
        if not gene_df.empty:
            gene_df = gene_df[gene_df["chrom"] == chrom_r]
            lo_bp = int(lo_m * BP_PER_MBP)
            hi_bp = int(hi_m * BP_PER_MBP)
            gene_df = gene_df[(gene_df["end"] >= lo_bp) & (gene_df["start"] <= hi_bp)]
            for key in list(prepared_gene_by_source.keys()):
                src_df = prepared_gene_by_source[key]
                if src_df.empty:
                    continue
                src_df = src_df[src_df["chrom"] == chrom_r]
                src_df = src_df[(src_df["end"] >= lo_bp) & (src_df["start"] <= hi_bp)]
                prepared_gene_by_source[key] = src_df
        if not methylation_df.empty:
            methylation_df = methylation_df[methylation_df["chrom"] == chrom_r]
            methylation_df = methylation_df[
                (methylation_df["pos"] >= lo_bp) & (methylation_df["pos"] <= hi_bp)
            ]
            for key in list(prepared_methylation_by_source.keys()):
                src_df = prepared_methylation_by_source[key]
                if src_df.empty:
                    continue
                src_df = src_df[src_df["chrom"] == chrom_r]
                src_df = src_df[(src_df["pos"] >= lo_bp) & (src_df["pos"] <= hi_bp)]
                prepared_methylation_by_source[key] = src_df
        if prepared_methylation_cov_by_source:
            for key in list(prepared_methylation_cov_by_source.keys()):
                src_cov = prepared_methylation_cov_by_source[key]
                filtered = filter_methylation_cov_by_region(
                    src_cov,
                    chrom=chrom_r,
                    lo_bp=lo_bp,
                    hi_bp=hi_bp,
                )
                prepared_methylation_cov_by_source[key] = filtered
        chrom_order = [chrom_r]
        plot_chrom_sizes = {chrom_r: chrom_sizes[chrom_r]}
        x_range = (max(0.0, lo_m - (hi_m - lo_m) * 0.02), hi_m + (hi_m - lo_m) * 0.02)
    else:
        if prepared_methylation_cov_by_source:
            pseudo_chroms = sorted(
                {
                    chrom
                    for src in prepared_methylation_cov_by_source.values()
                    for chrom in src.keys()
                }
            )
            order_hits = (
                pd.DataFrame({"chrom": pseudo_chroms, "pos": [1] * len(pseudo_chroms)})
                if pseudo_chroms
                else methylation_df
            )
        else:
            order_hits = hits if not hits.empty else (gene_df if not gene_df.empty else methylation_df)
        chrom_order, filter_note = _resolve_chrom_order(
            chrom_sizes,
            order_hits,
            chromosomes_only,
            chromosome_sort=chromosome_sort,
        )
        if filter_note:
            logger.warning("TE genome map: %s", filter_note)
        data_chroms = _chromosomes_with_plottable_data(
            hits,
            gene_df,
            methylation_df,
            prepared_methylation_cov_by_source,
        )
        # Keep explicitly selected chromosomes even when they have no TE/annotation hits.
        if data_chroms and len(chrom_order) > 1 and not chromosomes_only:
            chrom_order, dropped_empty = _prune_chromosomes_without_data(
                chrom_order,
                data_chroms,
            )
            if dropped_empty:
                sample = ", ".join(dropped_empty[:4])
                if len(dropped_empty) > 4:
                    sample += f", … (+{len(dropped_empty) - 4} more)"
                logger.info(
                    "TE genome map: skipped %d chromosome(s) with no TE/annotation/methylation data: %s",
                    len(dropped_empty),
                    sample,
                )
        plot_chrom_sizes = {c: chrom_sizes[c] for c in chrom_order}

    if not chrom_order:
        hint = _empty_chrom_order_hint(chrom_sizes, hits if not hits.empty else gene_df, chromosomes_only)
        title = "Genome BLAST map — no chromosomes to display"
        if hint:
            title += f"<br><sup>{hint}</sup>"
        fig = go.Figure()
        fig.update_layout(
            title=title,
            template="plotly_white",
        )
        return fig

    include_gene_tracks = not gene_df.empty
    gene_hits_by_chrom: dict[str, pd.DataFrame] = {}
    gene_hits_by_chrom_source: dict[str, dict[str, pd.DataFrame]] = {}
    if include_gene_tracks:
        gene_hits_by_chrom = {str(k): v for k, v in gene_df.groupby("chrom", sort=False)}
        for key, src_df in prepared_gene_by_source.items():
            if src_df.empty:
                gene_hits_by_chrom_source[key] = {}
            else:
                gene_hits_by_chrom_source[key] = {
                    str(k): v for k, v in src_df.groupby("chrom", sort=False)
                }

    split_strand_gene_tracks = bool(split_strand_tracks and include_gene_tracks)
    resolved_gene_sources = list(gene_sources or [])
    if include_gene_tracks and not resolved_gene_sources:
        resolved_gene_sources = [
            _GeneSourceSpec(
                key="gene_src_0",
                label="genes",
                split_strand=split_strand_gene_tracks,
            )
        ]
    elif split_strand_gene_tracks and resolved_gene_sources:
        resolved_gene_sources = [
            replace(src, split_strand=True) for src in resolved_gene_sources
        ]

    include_methylation_tracks = bool(prepared_methylation_cov_by_source) or not methylation_df.empty
    methylation_hits_by_chrom_source: dict[str, dict[str, MethylationChromCov]] = {}
    if include_methylation_tracks:
        if prepared_methylation_cov_by_source:
            methylation_hits_by_chrom_source = prepared_methylation_cov_by_source
        else:
            for key, src_df in prepared_methylation_by_source.items():
                if src_df.empty:
                    methylation_hits_by_chrom_source[key] = {}
                else:
                    methylation_hits_by_chrom_source[key] = partition_cov_by_chrom_arrays(
                        src_df,
                        coverage_min=max(1, int(methylation_coverage_min)),
                    )
    resolved_methylation_sources = list(methylation_sources or [])
    if include_methylation_tracks and not resolved_methylation_sources:
        resolved_methylation_sources = [
            _MethylationSourceSpec(key="meth_src_0", label="methylation")
        ]

    show_te_diff = _te_diff_tracks_enabled(
        te_sources,
        show_te_diff_track=show_te_diff_track,
    )
    show_meth_diff = _meth_diff_tracks_enabled(
        resolved_methylation_sources if include_methylation_tracks else None,
        show_meth_diff_track=show_meth_diff_track,
    )
    tracks = _ensure_unique_track_labels(
        _genome_map_tracks(
            chrom_order,
            te_sources=te_sources or [],
            gene_sources=resolved_gene_sources if include_gene_tracks else None,
            include_gene_tracks=include_gene_tracks,
            methylation_sources=resolved_methylation_sources if include_methylation_tracks else None,
            include_methylation_tracks=include_methylation_tracks,
            split_strand_gene_tracks=split_strand_gene_tracks,
            show_te_diff=show_te_diff,
            show_meth_diff=show_meth_diff,
        )
    )
    track_labels = [t.label for t in tracks]
    track_methylation_subtracks = {
        t.label: t.methylation_subtrack for t in tracks if t.methylation_subtrack
    }
    n_tracks = len(tracks)
    y_pos = _track_y_positions(tracks, show_te_diff=show_te_diff, show_meth_diff=show_meth_diff)
    any_split_tracks = any(spec.split_strand for spec in (te_sources or [])) or split_strand_gene_tracks
    multi_te_sources = len(te_sources or []) > 1
    gene_only = include_gene_tracks and hits.empty

    max_mbp = max(plot_chrom_sizes[c] / BP_PER_MBP for c in chrom_order)
    chrom_start_x = float(x_range[0]) if x_range is not None else 0.0
    view_hi_mbp = float(x_range[1]) if x_range is not None else max_mbp * 1.02
    track_pad_mbp = _track_pad_mbp(view_hi_mbp - chrom_start_x)
    track_x_off = track_pad_mbp

    if hits.empty:
        zmin, zmax = 70.0, 100.0
    else:
        zmin = float(hits["pident"].min())
        zmax = max(100.0, float(hits["pident"].max()))

    max_bar_height = (
        ROW_Y_STEP - 2 * TRACK_HALF - BAR_GAP_BELOW_TRACK - BAR_GAP_ABOVE_NEXT
    )
    max_bar_height = max(0.32, max_bar_height)

    hits_by_chrom: dict[str, pd.DataFrame] = {}
    if not hits.empty:
        hits_by_chrom = {str(k): v for k, v in hits.groupby("chrom", sort=False)}

    te_hits_by_chrom_source: dict[str, dict[str, pd.DataFrame]] = {}
    for key, src_df in prepared_by_source.items():
        if src_df.empty:
            te_hits_by_chrom_source[key] = {}
        else:
            te_hits_by_chrom_source[key] = {
                str(k): v for k, v in src_df.groupby("chrom", sort=False)
            }

    track_diff_bins: dict[str, tuple[np.ndarray, np.ndarray, str, str]] = {}
    if show_te_diff and te_sources and len(te_sources) == 2:
        track_diff_bins = _compute_all_track_diff_bins(
            chrom_order=chrom_order,
            plot_chrom_sizes=plot_chrom_sizes,
            te_hits_by_chrom_source=te_hits_by_chrom_source,
            te_sources=te_sources,
            bin_size=bin_size,
            num_workers=num_workers,
            cancel=cancel,
        )

    te_diff_peak = 1.0
    if track_diff_bins:
        te_diff_max, te_diff_peak = _shared_te_diff_scale_max(track_diff_bins)
    else:
        te_diff_max = 1.0
    shared_diff_norm = float(te_diff_peak) if show_te_diff else 1.0
    diff_bar_half = max(0.28, max_bar_height * DIFF_TRACK_BAR_HALF_FRAC) if show_te_diff else 0.06

    use_category_density = color_by in {
        TE_GENOME_MAP_COLOR_ORDER,
        TE_GENOME_MAP_COLOR_SUPERFAMILY,
        TE_GENOME_MAP_COLOR_FAMILY,
    }
    categories_for_legend: list[str] = []
    if use_category_density and not hits.empty:
        categories_for_legend = _legend_categories(hits, color_by)
    cat_col = (
        "order"
        if color_by == TE_GENOME_MAP_COLOR_ORDER
        else "superfamily"
        if color_by == TE_GENOME_MAP_COLOR_SUPERFAMILY
        else "family"
        if use_category_density
        else None
    )
    category_label = {
        TE_GENOME_MAP_COLOR_ORDER: "Order",
        TE_GENOME_MAP_COLOR_SUPERFAMILY: "Superfamily",
        TE_GENOME_MAP_COLOR_FAMILY: "Family",
    }.get(color_by, "TE category")
    region_window = (x_range[0], x_range[1]) if region is not None and x_range is not None else None
    methylation_cg_by_chrom: dict[str, np.ndarray] = {}
    if include_methylation_tracks or show_meth_diff:
        _emit("Precomputing methylation reference C+G bins…")
        methylation_cg_by_chrom = precompute_methylation_cg_totals(
            genome_fasta,
            plot_chrom_sizes,
            bin_size,
            chromosomes=chrom_order,
            num_workers=num_workers,
        )
    track_meth_diff_bins: dict[str, tuple[np.ndarray, np.ndarray, str, str]] = {}
    if show_meth_diff and resolved_methylation_sources and len(resolved_methylation_sources) == 2:
        track_meth_diff_bins = _compute_methylation_diff_bins(
            chrom_order=chrom_order,
            plot_chrom_sizes=plot_chrom_sizes,
            methylation_hits_by_chrom_source=methylation_hits_by_chrom_source,
            methylation_sources=resolved_methylation_sources,
            bin_size=bin_size,
            genome_fasta=genome_fasta,
            region=region_window,
            cg_by_chrom=methylation_cg_by_chrom,
            num_workers=num_workers,
            cancel=cancel,
        )
    meth_diff_peak = 0.1
    if track_meth_diff_bins:
        meth_diff_max, meth_diff_peak = _shared_meth_diff_scale_max(track_meth_diff_bins)
    else:
        meth_diff_max = 0.1
    shared_meth_diff_norm = float(meth_diff_peak) if show_meth_diff else 0.1
    meth_diff_bar_half = (
        max(0.28, max_bar_height * DIFF_TRACK_BAR_HALF_FRAC) if show_meth_diff else 0.06
    )
    track_cat_bins, track_bins, global_max_count = _bin_all_tracks_for_genome_map(
        tracks,
        hits_by_chrom,
        plot_chrom_sizes,
        bin_size=bin_size,
        density_mode=density_mode,
        use_category_density=use_category_density,
        categories_for_legend=categories_for_legend,
        cat_col=cat_col,
        region=region_window,
        num_workers=num_workers,
        num_hits=len(hits)
        + len(gene_df)
        + (
            sum(
                methylation_cov_site_count(cov)
                for src in methylation_hits_by_chrom_source.values()
                for cov in src.values()
            )
            if methylation_hits_by_chrom_source
            else len(methylation_df)
        ),
        collapse_overlaps=collapse_overlaps,
        annotation_density_mode=annotation_density_mode,
        gene_hits_by_chrom=gene_hits_by_chrom,
        te_hits_by_chrom_source=te_hits_by_chrom_source,
        gene_hits_by_chrom_source=gene_hits_by_chrom_source,
        cancel=cancel,
    )
    if include_methylation_tracks:
        _emit("Binning methylation tracks…")
        track_bins.update(
            _bin_methylation_tracks_for_genome_map(
                tracks,
                methylation_hits_by_chrom_source,
                plot_chrom_sizes,
                bin_size=bin_size,
                coverage_min=max(1, int(methylation_coverage_min)),
                region=region_window,
                genome_fasta=genome_fasta,
                cg_by_chrom=methylation_cg_by_chrom,
                methylation_sources=resolved_methylation_sources,
                num_workers=num_workers,
                cancel=cancel,
            )
        )

    global_max_stack = float(global_max_count)
    if use_category_density and categories_for_legend:
        global_max_stack = 1.0
        for track in tracks:
            if track.kind in {"gene", "te_diff", "meth_diff", "methylation"}:
                continue
            sample = track_cat_bins.get((track.label, categories_for_legend[0]))
            if sample is None:
                continue
            totals = np.zeros(len(sample[1]), dtype=float)
            for cat in categories_for_legend:
                cat_bins = track_cat_bins.get((track.label, cat))
                if cat_bins is None:
                    continue
                totals += cat_bins[1]
            if totals.size:
                global_max_stack = max(global_max_stack, float(totals.max()))

    track_density_max = _per_track_density_max(
        tracks,
        track_bins,
        track_cat_bins,
        use_category_density=use_category_density,
        categories_for_legend=categories_for_legend,
    )
    masking_density_max, consensus_density_max = _density_max_for_te_tracks(
        tracks,
        track_density_max,
        te_sources,
        use_category_density=use_category_density,
        track_cat_bins=track_cat_bins,
        categories_for_legend=categories_for_legend,
    )
    te_density_max = max(masking_density_max, consensus_density_max)
    gene_density_max = 1.0
    for track in tracks:
        if track.kind == "gene":
            vmax = track_density_max.get(track.label, 1.0)
            gene_density_max = max(gene_density_max, vmax)
    methylation_norm_max: dict[str, float] = {
        METH_SUBTRACK_CG_TOTAL: 1.0,
        METH_SUBTRACK_CG_COVERED: 1.0,
        METH_SUBTRACK_CG_COVERAGE_PCT: METH_COVERAGE_PCT_AXIS_MAX,
        METH_SUBTRACK_METH_LEVEL: METH_LEVEL_AXIS_MAX,
    }
    for track in tracks:
        if track.kind != "methylation" or not track.methylation_subtrack:
            continue
        vmax = track_density_max.get(track.label, 1.0)
        sub = track.methylation_subtrack
        if sub == METH_SUBTRACK_CG_COVERAGE_PCT:
            vmax = min(float(METH_COVERAGE_PCT_AXIS_MAX), max(vmax, 1e-6))
        elif sub == METH_SUBTRACK_METH_LEVEL:
            vmax = min(float(METH_LEVEL_AXIS_MAX), max(vmax, 1e-6))
        methylation_norm_max[sub] = max(methylation_norm_max.get(sub, 1.0), vmax)
    cg_count_shared = shared_cg_count_axis_max(methylation_norm_max)
    methylation_norm_max[METH_SUBTRACK_CG_TOTAL] = cg_count_shared
    methylation_norm_max[METH_SUBTRACK_CG_COVERED] = cg_count_shared

    def _density_norm_for(track: _GenomeMapTrack) -> float:
        if track.kind == "te_diff":
            return shared_diff_norm
        if track.kind == "meth_diff":
            return shared_meth_diff_norm
        if track.kind == "methylation":
            sub = track.methylation_subtrack or METH_SUBTRACK_CG_TOTAL
            return methylation_norm_max.get(sub, 1.0)
        if track.kind == "gene":
            return gene_density_max
        input_kind = _te_input_kind_for_source(te_sources, track.source_key)
        if input_kind == TE_INPUT_KIND_CONSENSUS:
            return consensus_density_max
        return masking_density_max

    def _density_scale_for(track: _GenomeMapTrack) -> float:
        bar_h = max_bar_height
        if track.kind == "te":
            bar_h = _density_axis_height_for_track(
                track.label,
                track_labels,
                y_pos,
                max_bar_height=max_bar_height,
                top_first=True,
            )
        return bar_h / max(_density_norm_for(track), 1e-12)

    def _track_axis_norm_max(track: _GenomeMapTrack) -> float:
        if track.kind == "te_diff":
            return shared_diff_norm
        if track.kind == "meth_diff":
            return shared_meth_diff_norm
        if track.kind == "methylation":
            sub = track.methylation_subtrack or METH_SUBTRACK_CG_TOTAL
            return methylation_norm_max.get(sub, 1.0)
        if track.kind == "gene":
            return gene_density_max
        if track.kind == "te":
            return _density_norm_for(track)
        return masking_density_max

    methylation_colors = {
        METH_SUBTRACK_CG_TOTAL: METH_COLOR_CG_TOTAL,
        METH_SUBTRACK_CG_COVERED: METH_COLOR_CG_COVERED,
        METH_SUBTRACK_CG_COVERAGE_PCT: METH_COLOR_CG_COVERAGE_PCT,
        METH_SUBTRACK_METH_LEVEL: METH_COLOR_METH_LEVEL,
    }
    track_density_norm = {t.label: _track_axis_norm_max(t) for t in tracks}
    track_input_kinds = _track_input_kind_map(tracks, te_sources)

    bin_span_mbp = bin_size / BP_PER_MBP
    cat_colors = _category_color_map(categories_for_legend, color_by=color_by)

    fig = go.Figure()
    chrom_grey_shape_idx: list[int] = []
    chrom_scatter_idx: list[int] = []

    for track in tracks:
        chrom = track.chrom
        chrom_len = plot_chrom_sizes[chrom]
        length_mbp = chrom_len / BP_PER_MBP
        y_base = float(y_pos[track.label])
        strand_note = _track_strand_caption(track)

        track_vis_half = (
            meth_diff_bar_half
            if track.kind == "meth_diff"
            else diff_bar_half
            if track.kind == "te_diff"
            else TRACK_HALF
        )
        chrom_grey_shape_idx.append(len(fig.layout.shapes) if fig.layout.shapes else 0)
        fig.add_shape(
            type="rect",
            x0=track_x_off,
            x1=track_x_off + length_mbp,
            y0=y_base - track_vis_half,
            y1=y_base + track_vis_half,
            fillcolor="#e0e0e0",
            line={"color": "#a8a8a8", "width": 1},
            layer="below",
        )

        n_hover = min(HOVER_TRACK_POINTS_PER_CHROM, max(4, int(length_mbp * 4)))
        x_track = np.linspace(track_x_off, track_x_off + length_mbp, n_hover)
        track_hover = [
            (
                f"<b>{track.label}</b><br>"
                f"Position: {(x - track_x_off):.3f} Mbp ({int((x - track_x_off) * BP_PER_MBP):,} bp)<br>"
                f"Chromosome length: {length_mbp:.3f} Mbp"
                + (f"<br>Strand: {strand_note}" if strand_note else "")
            )
            for x in x_track
        ]
        chrom_scatter_idx.append(len(fig.data))
        fig.add_trace(
            go.Scatter(
                x=x_track,
                y=[y_base] * n_hover,
                mode="markers",
                marker={"size": 14, "opacity": 0, "color": "#e0e0e0"},
                hovertext=track_hover,
                hoverinfo="text",
                showlegend=False,
            )
        )

        window_kb = max(5, bin_size // 1000)
        track_density_scale = _density_scale_for(track)
        if track.kind == "te_diff":
            diff_entry = track_diff_bins.get(track.label)
            if diff_entry is None:
                continue
            bins_a, bins_b, label_a, label_b = diff_entry
            half_track = diff_bar_half
            diff_scale = half_track / max(shared_diff_norm, 1e-12)
            n_bins = len(bins_a)
            lo_mbp_all = np.arange(n_bins, dtype=float) * bin_span_mbp
            hi_mbp_all = np.minimum((np.arange(n_bins, dtype=float) + 1.0) * bin_span_mbp, length_mbp)

            def _add_diff_bars(
                values: np.ndarray,
                *,
                upward: bool,
                sample_label: str,
                color: str,
            ) -> None:
                active = values > 0
                if not np.any(active):
                    return
                bin_indices = np.where(active)[0]
                heights = values[active] * diff_scale
                if not upward:
                    heights = -heights
                lo_mbp = lo_mbp_all[active]
                hi_mbp = hi_mbp_all[active]
                x_vals = (bin_indices.astype(float) + 0.5) * bin_span_mbp + track_x_off
                bar_widths = np.full(len(bin_indices), bin_span_mbp * 0.98)
                hover = [
                    (
                        f"<b>{track.label}</b><br>"
                        f"Only in {sample_label}<br>"
                        f"{lo:.2f}–{hi:.2f} Mbp ({window_kb} kb window)<br>"
                        f"Exclusive TE bp: {int(val):,}"
                    )
                    for lo, hi, val in zip(lo_mbp, hi_mbp, values[active], strict=False)
                ]
                custom = np.column_stack(
                    [
                        np.full(len(bin_indices), track.label),
                        bin_indices,
                        lo_mbp,
                        hi_mbp,
                        np.full(len(bin_indices), y_base),
                        values[active],
                    ]
                )
                fig.add_trace(
                    go.Bar(
                        x=x_vals,
                        y=heights,
                        base=[y_base] * len(x_vals),
                        width=bar_widths,
                        marker={"color": color, "line": {"width": 0}},
                        customdata=custom,
                        hovertext=hover,
                        hoverinfo="text",
                        name=sample_label,
                        legendgroup=TE_DIFF_LEGENDGROUP,
                        alignmentgroup=track.label,
                        offsetgroup=track.label,
                        showlegend=False,
                    )
                )

            _add_diff_bars(
                bins_a,
                upward=True,
                sample_label=label_a,
                color=TE_DIFF_COLOR_A,
            )
            _add_diff_bars(
                bins_b,
                upward=False,
                sample_label=label_b,
                color=TE_DIFF_COLOR_B,
            )
            continue
        if track.kind == "meth_diff":
            diff_entry = track_meth_diff_bins.get(track.label)
            if diff_entry is None:
                continue
            bins_a, bins_b, label_a, label_b = diff_entry
            half_track = meth_diff_bar_half
            diff_scale = half_track / max(shared_meth_diff_norm, 1e-12)
            n_bins = len(bins_a)
            lo_mbp_all = np.arange(n_bins, dtype=float) * bin_span_mbp
            hi_mbp_all = np.minimum((np.arange(n_bins, dtype=float) + 1.0) * bin_span_mbp, length_mbp)

            def _add_meth_diff_bars(
                values: np.ndarray,
                *,
                upward: bool,
                sample_label: str,
                color: str,
            ) -> None:
                active = values > 0
                if not np.any(active):
                    return
                bin_indices = np.where(active)[0]
                heights = values[active] * diff_scale
                if not upward:
                    heights = -heights
                lo_mbp = lo_mbp_all[active]
                hi_mbp = hi_mbp_all[active]
                x_vals = (bin_indices.astype(float) + 0.5) * bin_span_mbp + track_x_off
                bar_widths = np.full(len(bin_indices), bin_span_mbp * 0.98)
                hover = [
                    (
                        f"<b>{track.label}</b><br>"
                        f"Higher in {sample_label}<br>"
                        f"{lo:.2f}–{hi:.2f} Mbp ({window_kb} kb window)<br>"
                        f"Methylation level diff: {float(val):.3f}"
                    )
                    for lo, hi, val in zip(lo_mbp, hi_mbp, values[active], strict=False)
                ]
                custom = np.column_stack(
                    [
                        np.full(len(bin_indices), track.label),
                        bin_indices,
                        lo_mbp,
                        hi_mbp,
                        np.full(len(bin_indices), y_base),
                    ]
                )
                fig.add_trace(
                    go.Bar(
                        x=x_vals,
                        y=heights,
                        base=[y_base] * len(x_vals),
                        width=bar_widths,
                        marker={"color": color, "line": {"width": 0}},
                        customdata=custom,
                        hovertext=hover,
                        hoverinfo="text",
                        name=sample_label,
                        legendgroup=METH_DIFF_LEGENDGROUP,
                        alignmentgroup=track.label,
                        offsetgroup=track.label,
                        showlegend=False,
                    )
                )

            _add_meth_diff_bars(
                bins_a,
                upward=True,
                sample_label=label_a,
                color=METH_DIFF_COLOR_A,
            )
            _add_meth_diff_bars(
                bins_b,
                upward=False,
                sample_label=label_b,
                color=METH_DIFF_COLOR_B,
            )
            continue
        if use_category_density and categories_for_legend and track.kind not in {
            "gene",
            "methylation",
        }:
            stacked_base = None
            for cat in categories_for_legend:
                binned = track_cat_bins.get((track.label, cat))
                if binned is None:
                    continue
                centers_mbp, bin_values, mean_pident, hit_counts, te_bp = binned
                active = bin_values > 0
                if not np.any(active):
                    continue
                density = bin_values[active] * track_density_scale
                bin_indices = np.where(active)[0]
                lo_mbp = bin_indices * bin_span_mbp
                hi_mbp = np.minimum((bin_indices + 1) * bin_span_mbp, length_mbp)
                x_vals = (bin_indices.astype(float) + 0.5) * bin_span_mbp + track_x_off
                bar_widths = np.full(len(bin_indices), bin_span_mbp * 0.98)
                if stacked_base is None:
                    stacked_base = {idx: (y_base + TRACK_HALF + BAR_GAP_BELOW_TRACK) for idx in bin_indices}
                hover = []
                row_ground = y_base + TRACK_HALF + BAR_GAP_BELOW_TRACK
                for idx, lo_m, hi_m, hits_n, te_len in zip(
                    bin_indices,
                    lo_mbp,
                    hi_mbp,
                    hit_counts[active],
                    te_bp[active],
                ):
                    hover.append(
                        f"<b>{track.label}</b><br>"
                        f"{lo_m:.2f}–{hi_m:.2f} Mbp ({window_kb} kb window)"
                        + (f"<br>Strand: {strand_note}" if strand_note else "")
                        + f"<br>{category_label}: {cat}<br>"
                        f"Hits: {int(hits_n)}<br>"
                        f"Total TE length: {int(te_len):,} bp ({te_len / 1000:.1f} kb)"
                    )
                bases = [stacked_base.get(idx, y_base + TRACK_HALF + BAR_GAP_BELOW_TRACK) for idx in bin_indices]
                for idx, d, b in zip(bin_indices, density, bases):
                    stacked_base[idx] = b + d
                custom = np.column_stack(
                    [
                        np.full(len(bin_indices), track.label),
                        bin_indices,
                        lo_mbp,
                        hi_mbp,
                        np.full(len(bin_indices), row_ground),
                        bin_values[active],
                    ]
                )
                fig.add_trace(
                    go.Bar(
                        x=x_vals,
                        y=density,
                        base=bases,
                        width=bar_widths,
                        marker={"color": cat_colors.get(cat, TE_EXTRA_TYPE_COLORS[0]), "line": {"width": 0}},
                        customdata=custom,
                        hovertext=hover,
                        hoverinfo="text",
                        name=cat,
                        legendgroup=cat,
                        alignmentgroup=track.label,
                        offsetgroup=track.label,
                        showlegend=False,
                    )
                )
        else:
            binned = track_bins.get(track.label)
            if binned is None:
                continue
            centers_mbp, bin_values, mean_pident, hit_counts, te_bp = binned
            active = bin_values > 0
            if not np.any(active):
                continue

            density = bin_values[active] * track_density_scale
            if track.kind == "methylation":
                sub = track.methylation_subtrack or METH_SUBTRACK_CG_TOTAL
                bar_color = methylation_colors.get(sub, METH_COLOR_CG_TOTAL)
                colors = [bar_color] * int(np.count_nonzero(active))
            elif track.kind == "gene":
                colors = ["#4a90d9"] * int(np.count_nonzero(active))
            else:
                colors = _pident_colorscale(mean_pident[active], zmin, zmax)
            bar_base = y_base + TRACK_HALF + BAR_GAP_BELOW_TRACK
            bin_indices = np.where(active)[0]
            lo_mbp = bin_indices * bin_span_mbp
            hi_mbp = np.minimum((bin_indices + 1) * bin_span_mbp, length_mbp)
            x_vals = (bin_indices.astype(float) + 0.5) * bin_span_mbp + track_x_off
            bar_widths = np.full(len(bin_indices), bin_span_mbp * 0.98)
            hover = []
            if track.kind == "methylation":
                sub = track.methylation_subtrack or ""
                for lo_m, hi_m, hits_n, counted_n, meth_val in zip(
                    lo_mbp,
                    hi_mbp,
                    hit_counts[active],
                    te_bp[active],
                    bin_values[active],
                ):
                    if sub == METH_SUBTRACK_CG_TOTAL:
                        detail = f"Total C+G number: {int(hits_n)}"
                    elif sub == METH_SUBTRACK_CG_COVERED:
                        detail = f"Covered C+G number: {int(hits_n)}"
                    elif sub == METH_SUBTRACK_CG_COVERAGE_PCT:
                        detail = (
                            f"Covered fraction: {float(meth_val):.3f} "
                            f"({int(counted_n)} / {int(hits_n)})"
                        )
                    else:
                        detail = (
                            f"Methylation level: {float(meth_val):.3f} "
                            f"(from {int(hits_n)} covered lines)"
                        )
                    hover.append(
                        f"<b>{track.label}</b><br>"
                        f"{lo_m:.2f}–{hi_m:.2f} Mbp ({window_kb} kb window)"
                        + (f"<br>{strand_note}" if strand_note else "")
                        + f"<br>{detail}"
                    )
            else:
                length_label = "gene" if track.kind == "gene" else "TE"
                feature_label = "Genes" if track.kind == "gene" else "Features"
                for lo_m, hi_m, hits_n, te_len, pid in zip(
                    lo_mbp,
                    hi_mbp,
                    hit_counts[active],
                    te_bp[active],
                    np.nan_to_num(mean_pident[active], nan=0.0),
                ):
                    hover.append(
                        f"<b>{track.label}</b><br>"
                        f"{lo_m:.2f}–{hi_m:.2f} Mbp ({window_kb} kb window)"
                        + (f"<br>Strand: {strand_note}" if strand_note else "")
                        + f"<br>{feature_label}: {int(hits_n)}<br>"
                        f"Total {length_label} length: {int(te_len):,} bp ({te_len / 1000:.1f} kb)"
                        + (f"<br>Mean identity: {pid:.1f}%" if track.kind != "gene" else "")
                    )
            custom = np.column_stack(
                [
                    np.full(len(bin_indices), track.label),
                    bin_indices,
                    lo_mbp,
                    hi_mbp,
                    np.full(len(bin_indices), bar_base),
                    bin_values[active],
                ]
            )

            fig.add_trace(
                go.Bar(
                    x=x_vals,
                    y=density,
                    base=[bar_base] * len(x_vals),
                    width=bar_widths,
                    marker={"color": colors, "line": {"width": 0}},
                    customdata=custom,
                    hovertext=hover,
                    hoverinfo="text",
                    showlegend=False,
                )
            )

    if use_category_density and categories_for_legend:
        for cat in categories_for_legend:
            fig.add_trace(
                go.Scatter(
                    x=[None],
                    y=[None],
                    mode="markers",
                    marker={
                        "size": 10,
                        "symbol": "square",
                        "color": cat_colors.get(cat, TE_EXTRA_TYPE_COLORS[0]),
                    },
                    name=cat,
                    legendgroup=cat,
                    showlegend=True,
                    hoverinfo="skip",
                )
            )
    elif not hits.empty:
        fig.add_trace(
            go.Scatter(
                x=[None],
                y=[None],
                mode="markers",
                marker={
                    "colorscale": [
                        [0.0, "rgb(230,90,60)"],
                        [0.5, "rgb(200,180,80)"],
                        [1.0, "rgb(60,180,140)"],
                    ],
                    "cmin": zmin,
                    "cmax": zmax,
                    "colorbar": {
                        "title": "% ident",
                        "thickness": 14,
                        "len": 0.5,
                        "x": 1.02,
                    },
                    "showscale": True,
                },
                hoverinfo="skip",
                showlegend=False,
            )
        )

    hit_n = len(hits)
    gene_n = len(gene_df)
    if gene_only:
        title = f"Gene density in genome ({plot_source})"
    else:
        title = f"TE positions in genome ({plot_source})"
    if query_label:
        title += f" — {query_label}"
    if region is not None:
        title += f" — {region[0]} {region[1]:.2f}–{region[2]:.2f} Mbp"
    summary_parts: list[str] = []
    if multi_te_sources and te_sources:
        summary_parts.append(f"{len(te_sources)} TE inputs")
    if hit_n:
        summary_parts.append(f"{hit_n} TE hits")
    if gene_n:
        summary_parts.append(f"{gene_n} gene features")
    if summary_parts:
        title += f" ({', '.join(summary_parts)}, {len(chrom_order)} chr."
        if n_tracks > len(chrom_order):
            title += f", {n_tracks} tracks"
        title += ")"
    window_kb = max(5, bin_size // 1000)
    y_top = max(y_pos.values()) + TRACK_HALF + max_bar_height + 0.24 if y_pos else (
        (n_tracks - 1) * ROW_Y_STEP + TRACK_HALF + max_bar_height + 0.24
    )
    diff_half_for_margin = max(
        diff_bar_half if show_te_diff else 0.0,
        meth_diff_bar_half if show_meth_diff else 0.0,
    )
    if y_pos and (show_te_diff or show_meth_diff):
        y_bottom = min(y_pos.values()) - diff_half_for_margin - DIFF_AXIS_LABEL_PAD
    else:
        y_bottom = -0.35
    if gene_only:
        ann_mode = annotation_density_mode
        if ann_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH:
            density_caption = f"total element length per {window_kb} kb window"
        else:
            density_caption = f"element hit number per {window_kb} kb window"
        colour_caption = "colour = annotation density"
    elif density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH:
        density_caption = f"total TE length per {window_kb} kb window"
    else:
        density_caption = f"TE hits number per {window_kb} kb window"
    if not gene_only:
        if use_category_density:
            colour_caption = (
                "colour = TE order"
                if color_by == TE_GENOME_MAP_COLOR_ORDER
                else "colour = TE superfamily"
                if color_by == TE_GENOME_MAP_COLOR_SUPERFAMILY
                else "colour = TE family"
            )
        else:
            colour_caption = "colour = % ident"

    x_lo = 0.0
    x_hi = view_hi_mbp + track_pad_mbp
    x_tickvals, x_ticktext = _genomic_xaxis_ticks(chrom_start_x, view_hi_mbp, track_pad_mbp)
    density_paper_x = _density_paper_x_in_gutter()
    if gene_only:
        if annotation_density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH:
            density_unit = f"element bp / {window_kb} kb window"
        else:
            density_unit = f"elements / {window_kb} kb window"
    elif density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH:
        density_unit = f"TE bp / {window_kb} kb window"
    else:
        density_unit = f"hits / {window_kb} kb window"

    margin_left = _genome_map_margin_left(track_labels)
    track_kinds = {t.label: t.kind for t in tracks}
    track_chroms = {t.label: t.chrom for t in tracks}
    map_caption = (
        f"Grey bar = chromosome (hover for position); each coloured bar = one "
        f"{window_kb} kb window ({density_caption}); {colour_caption}."
        + (
            " Optional diff track (two TE inputs): exclusive TE bp per window "
            "(above = first input only, below = second only; shared diff-track scale)."
            if show_te_diff and te_sources and len(te_sources) == 2
            else ""
        )
        + (
            " RepeatMasker .out: separate forward and reverse tracks per TE input."
            if any_split_tracks and not split_strand_gene_tracks
            else ""
        )
        + (
            " GFF/GTF/BED/VCF annotations: separate forward and reverse tracks per chromosome."
            if split_strand_gene_tracks and not any(spec.split_strand for spec in (te_sources or []))
            else ""
        )
        + (
            " RepeatMasker .out and GFF/GTF/BED/VCF annotations: separate forward and reverse tracks."
            if split_strand_gene_tracks
            and any(spec.split_strand for spec in (te_sources or []))
            else ""
        )
        + (
            " Annotation tracks use features from GFF/GTF (gene rows), BED intervals, "
            "or VCF variants (SNP, INDEL, MNP, SV)."
            "with an independent y-axis scale."
            if include_gene_tracks and not gene_only
            else ""
        )
        + (
            " Bismark methylation: four tracks per file per chromosome "
            "(Total, Covered, Covered fraction, Methylation level 0–1)."
            if include_methylation_tracks
            else ""
        )
        + f" Mini y-axis in left padding ({density_unit}). "
        "Track labels in left margin. TE types legend on the right; "
        "floating panel (top right) for show/hide, colour, and order."
        + (
            f" Dashed orange line = one window ({window_kb} kb"
            + (" TE" if density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH else "")
            + ")."
            if show_window_ref_line
            else ""
        )
    )

    plot_width = 1120
    top_margin = _genome_map_top_margin(title, map_caption, width=plot_width)
    fig.update_layout(
        title={
            "text": (
                f"{title}<br><span style='font-size:10px;color:#57606a'>{map_caption}</span>"
            ),
            "x": 0.01,
            "xanchor": "left",
            "pad": {"t": 8, "b": 12},
        },
        template="plotly_white",
        height=max(380, int(68 * n_tracks) + 130),
        width=plot_width,
        autosize=False,
        bargap=0,
        # Manual base stacking per bin; overlay keeps category segments at the same x.
        barmode="overlay",
        xaxis={
            "title": "Mbp",
            "range": [x_lo, x_hi],
            "domain": [XAXIS_DOMAIN_START, 1.0],
            "autorange": False,
            "fixedrange": False,
            "visible": True,
            "showticklabels": True,
            "showline": True,
            "tickmode": "array",
            "tickvals": x_tickvals,
            "ticktext": x_ticktext,
            "showgrid": True,
            "gridcolor": "#f0f0f0",
            "zeroline": False,
            "rangeslider": {"visible": False},
        },
        dragmode="zoom",
        uirevision="te-genome-map",
        yaxis={
            "showticklabels": False,
            "showline": False,
            "ticks": "",
            "title": "",
            "range": [y_bottom, y_top],
            "autorange": False,
            "fixedrange": True,
            "showgrid": False,
        },
        margin={"l": margin_left, "r": 200 if use_category_density else 100, "t": top_margin, "b": 56},
        plot_bgcolor="white",
        paper_bgcolor="white",
        showlegend=bool(use_category_density and categories_for_legend),
        legend=(
            {
                "title": {"text": "TE types"},
                "orientation": "v",
                "yanchor": "top",
                "y": 1.0,
                "xanchor": "left",
                "x": 1.02,
                "font": {"size": 10},
            }
            if use_category_density and categories_for_legend
            else {}
        ),
    )
    chrom_label_ann_idx = _append_chromosome_labels(
        fig,
        chrom_order=track_labels,
        y_pos=y_pos,
        track_kinds=track_kinds,
        track_chroms=track_chroms,
        chrom_sizes=plot_chrom_sizes,
        diff_bar_half=diff_half_for_margin if (show_te_diff or show_meth_diff) else None,
    )
    te_length_mode = density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH
    annotation_length_mode = annotation_density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH
    scale_meta = _append_track_density_axes(
        fig,
        chrom_order=track_labels,
        y_pos=y_pos,
        max_bar_height=max_bar_height,
        density_max=te_density_max,
        window_kb=window_kb,
        te_length_mode=te_length_mode,
        tick_count=5,
        density_paper_x=density_paper_x,
        gene_density_max=gene_density_max if include_gene_tracks else None,
        annotation_length_mode=annotation_length_mode,
        track_kinds=track_kinds,
        track_density_norm=track_density_norm,
        track_input_kinds=track_input_kinds,
        track_methylation_subtracks=track_methylation_subtracks,
        diff_bar_half=diff_half_for_margin if (show_te_diff or show_meth_diff) else None,
    )
    scale_meta["densityNormMax"] = float(te_density_max)
    scale_meta["maskingDensityNormMax"] = float(masking_density_max)
    scale_meta["consensusDensityNormMax"] = float(consensus_density_max)
    if include_gene_tracks:
        scale_meta["geneDensityNormMax"] = float(gene_density_max)
        scale_meta["annotationLengthMode"] = bool(annotation_length_mode)
    if resolved_gene_sources:
        scale_meta["geneSources"] = [
            {"key": src.key, "label": src.label} for src in resolved_gene_sources
        ]
        scale_meta["multiGeneSources"] = len(resolved_gene_sources) > 1
    if te_sources:
        scale_meta["teSources"] = [
            {"key": src.key, "label": src.label} for src in te_sources
        ]
    panel_groups = _panel_track_groups_for_map(
        te_sources=list(te_sources or []),
        show_te_diff=bool(show_te_diff),
        resolved_gene_sources=resolved_gene_sources,
        include_gene_tracks=include_gene_tracks,
        resolved_methylation_sources=resolved_methylation_sources,
        include_methylation_tracks=include_methylation_tracks,
        show_meth_diff=bool(show_meth_diff),
    )
    if panel_groups:
        scale_meta["panelTrackGroups"] = _enrich_panel_track_groups_with_default_ymax(
            panel_groups,
            tracks=tracks,
            track_density_norm=track_density_norm,
            te_sources=list(te_sources or []),
            masking_density_max=float(masking_density_max),
            consensus_density_max=float(consensus_density_max),
            gene_density_max=float(gene_density_max),
            methylation_norm_max=methylation_norm_max,
            te_diff_density_max=float(shared_diff_norm) if show_te_diff else None,
            meth_diff_density_max=float(shared_meth_diff_norm) if show_meth_diff else None,
        )
    scale_meta["baseTrackOrder"] = list(track_labels)
    if include_methylation_tracks:
        scale_meta["hasMethylationTracks"] = True
        scale_meta["methylationCoverageMin"] = max(1, int(methylation_coverage_min))
        scale_meta["methylationDensityNormMaxBySubtrack"] = {
            "cg_count_shared": float(cg_count_shared),
            METH_SUBTRACK_CG_TOTAL: float(methylation_norm_max[METH_SUBTRACK_CG_TOTAL]),
            METH_SUBTRACK_CG_COVERED: float(methylation_norm_max[METH_SUBTRACK_CG_COVERED]),
            METH_SUBTRACK_CG_COVERAGE_PCT: float(
                methylation_norm_max[METH_SUBTRACK_CG_COVERAGE_PCT]
            ),
            METH_SUBTRACK_METH_LEVEL: float(methylation_norm_max[METH_SUBTRACK_METH_LEVEL]),
        }
        scale_meta["methColorCgTotal"] = METH_COLOR_CG_TOTAL
        scale_meta["methColorCgCovered"] = METH_COLOR_CG_COVERED
        scale_meta["methColorCgPct"] = METH_COLOR_CG_COVERAGE_PCT
        scale_meta["methColorLevel"] = METH_COLOR_METH_LEVEL
    scale_meta["sharedTeDensityAxis"] = True
    if show_te_diff and te_sources and len(te_sources) == 2:
        scale_meta["hasDiffTrack"] = True
        scale_meta["diffColorA"] = str(scale_meta.get("diffColorA") or TE_DIFF_COLOR_A)
        scale_meta["diffColorB"] = str(scale_meta.get("diffColorB") or TE_DIFF_COLOR_B)
        scale_meta["diffLabelA"] = te_sources[0].label
        scale_meta["diffLabelB"] = te_sources[1].label
        scale_meta["diffDensityNormMax"] = float(shared_diff_norm)
        scale_meta["diffDensityRefMax"] = float(te_diff_max)
        scale_meta["diffAxisTickCount"] = DIFF_AXIS_TICK_COUNT
        scale_meta["diffDensityPeak"] = float(te_diff_peak)
        scale_meta["diffBarHalf"] = float(diff_bar_half)
        scale_meta["sharedDiffAxis"] = True
    scale_meta["showWindowRefLine"] = bool(
        show_window_ref_line
    )
    scale_meta["windowRefDensity"] = _window_ref_density_value(
        window_kb=window_kb,
        te_length_mode=te_length_mode,
    )
    scale_meta["windowRefColor"] = str(
        scale_meta.get("windowRefColor") or DEFAULT_WINDOW_REF_COLOR
    )
    ref_meta = _append_window_ref_lines(
        fig,
        track_labels=track_labels,
        track_chroms=track_chroms,
        y_pos=y_pos,
        chrom_sizes=plot_chrom_sizes,
        track_x_off=track_x_off,
        max_bar_height=max_bar_height,
        density_norm_max=float(te_density_max),
        window_kb=window_kb,
        te_length_mode=te_length_mode,
        density_paper_x=density_paper_x,
        line_color=str(scale_meta.get("windowRefColor") or DEFAULT_WINDOW_REF_COLOR),
        track_kinds=track_kinds,
        track_density_norm=track_density_norm,
        show_lines=bool(scale_meta["showWindowRefLine"]),
    )
    for track_entry in scale_meta["tracks"]:
        ref_shapes = ref_meta.get(track_entry["chrom"], {})
        if ref_shapes:
            track_entry["windowRefLineShapeIdx"] = ref_shapes.get("lineShapeIdx")
            track_entry["windowRefTickShapeIdx"] = ref_shapes.get("tickShapeIdx")
    scale_meta["xRange"] = [float(x_lo), float(x_hi)]
    scale_meta["trackPadMbp"] = float(track_pad_mbp)
    scale_meta["binSpanMbp"] = float(bin_span_mbp)
    scale_meta["chromStartMbp"] = float(chrom_start_x)
    scale_meta["viewHiMbp"] = float(view_hi_mbp)
    scale_meta["xDomainStart"] = float(XAXIS_DOMAIN_START)
    scale_meta["densityPaperX"] = float(density_paper_x)
    scale_meta["xTickvals"] = x_tickvals
    scale_meta["xTicktext"] = x_ticktext
    scale_meta["marginLeft"] = margin_left
    scale_meta["chromLayout"] = {
        "order": list(track_labels),
        "chromOrder": list(chrom_order),
        "chromSizes": {str(k): int(v) for k, v in plot_chrom_sizes.items()},
        "naturalChromOrder": sorted(chrom_order, key=_natural_sort_key),
        "chromosomeSort": chromosome_sort,
        "splitStrandTracks": any_split_tracks,
        "multiTeSources": multi_te_sources,
        "rowYStep": ROW_Y_STEP,
        "yPosByLabel": {label: float(y_pos[label]) for label in track_labels},
        "topFirstOrder": True,
        "trackHalf": TRACK_HALF,
        "barGapBelow": BAR_GAP_BELOW_TRACK,
        "rowGapAbove": BAR_GAP_ABOVE_NEXT,
        "yMin": float(y_bottom),
        "yTopPad": 0.24,
        "diffAxisLabelPad": DIFF_AXIS_LABEL_PAD,
        "greyShapeIdx": chrom_grey_shape_idx,
        "scatterIdx": chrom_scatter_idx,
        "chromLabelAnnIdx": chrom_label_ann_idx,
        "chromLabelAnnIdxByLabel": {
            label: chrom_label_ann_idx[i] for i, label in enumerate(track_labels)
        },
    }
    scale_meta["chromLayout"]["trackKindsByLabel"] = {t.label: t.kind for t in tracks}
    scale_meta["chromLayout"]["trackGroupKeyByLabel"] = {
        t.label: key
        for t in tracks
        if (key := _track_group_key_for_track(t))
    }
    scale_meta["chromLayout"]["teSourceKeyByLabel"] = {
        t.label: t.source_key
        for t in tracks
        if t.kind == "te" and t.source_key
    }
    scale_meta["chromLayout"]["geneSourceKeyByLabel"] = {
        t.label: t.source_key
        for t in tracks
        if t.kind == "gene" and t.source_key
    }
    scale_meta["chromLayout"]["methylationSubtrackByLabel"] = {
        t.label: t.methylation_subtrack
        for t in tracks
        if t.methylation_subtrack
    }
    scale_meta["chromLayout"]["diffBeforeStep"] = float(DIFF_TRACK_BEFORE_STEP)
    scale_meta["chromLayout"]["diffAfterGap"] = float(DIFF_TRACK_AFTER_GAP)
    if show_te_diff:
        scale_meta["chromLayout"]["diffTrackHalf"] = float(diff_bar_half)
    if show_meth_diff and resolved_methylation_sources and len(resolved_methylation_sources) == 2:
        scale_meta["hasMethDiffTrack"] = True
        scale_meta["methDiffLabelA"] = resolved_methylation_sources[0].label
        scale_meta["methDiffLabelB"] = resolved_methylation_sources[1].label
        scale_meta["methDiffDensityNormMax"] = float(shared_meth_diff_norm)
        scale_meta["methDiffDensityRefMax"] = float(meth_diff_max)
        scale_meta["methDiffDensityPeak"] = float(meth_diff_peak)
        scale_meta["methDiffBarHalf"] = float(meth_diff_bar_half)
        scale_meta["methDiffColorA"] = METH_DIFF_COLOR_A
        scale_meta["methDiffColorB"] = METH_DIFF_COLOR_B
        scale_meta["chromLayout"]["methDiffTrackHalf"] = float(meth_diff_bar_half)
    if use_category_density and categories_for_legend:
        scale_meta["categoryColors"] = {
            cat: cat_colors.get(cat, TE_EXTRA_TYPE_COLORS[0]) for cat in categories_for_legend
        }
        scale_meta["categoryOrder"] = list(categories_for_legend)
    fig.update_layout(meta={"teGenomeScale": scale_meta})
    return fig


def write_blast_divergence_plot_html(
    te_input_files: list[Path] | Path,
    output_html: Path | None = None,
    *,
    query_label: str | None = None,
    bin_width: float = DEFAULT_DIVERGENCE_BIN_WIDTH,
) -> Path | None:
    """Write divergence histogram HTML from one or more TE hit tables.

    X axis: divergence fraction (1 - identity), starts at 0 on the left.
    Y axis: total aligned TE length (bp) summed in each divergence bin.
    """
    te_paths = (
        [te_input_files.expanduser().resolve()]
        if isinstance(te_input_files, Path)
        else [p.expanduser().resolve() for p in te_input_files]
    )
    te_paths = [p for p in te_paths if p.is_file()]
    if not te_paths:
        return None
    primary = te_paths[0]
    if output_html is None:
        output_html = divergence_plot_html_path_for_te_inputs(te_paths)
    else:
        output_html = output_html.expanduser().resolve()

    try:
        hits = load_combined_te_hits(te_paths)
    except Exception as exc:  # noqa: BLE001
        logger.warning("TE divergence plot: cannot read TE inputs: %s", exc)
        return None

    source = te_hit_plot_source_label(primary)
    if hits.empty:
        fig = go.Figure()
        fig.update_layout(
            title=f"TE divergence from {source} — no hits",
            xaxis={"title": "Divergence", "range": [0, 1]},
            yaxis={"title": "Total TE length (bp)"},
            template="plotly_white",
        )
    else:
        bw = max(0.001, float(bin_width))
        n_bins = int(np.ceil(1.0 / bw))
        edges = np.linspace(0.0, 1.0, n_bins + 1)
        centers = (edges[:-1] + edges[1:]) / 2.0
        # Stack / colour by TE class (order), matching landscape-plot legends
        # (LTR, LINE, DNA, …) rather than fine-grained superfamily names.
        if "order" in hits.columns:
            cat_series = hits["order"].astype(str)
        elif "superfamily" in hits.columns:
            cat_series = hits["superfamily"].astype(str)
        else:
            cat_series = pd.Series(["Unclassified"] * len(hits), index=hits.index, dtype=str)
        categories = _ordered_te_types({str(v) for v in cat_series if str(v)})
        if not categories:
            categories = ["Unclassified"]
            hits = hits.copy()
            cat_series = pd.Series("Unclassified", index=hits.index, dtype=str)
        cat_colors = _category_color_map(
            categories,
            color_by=TE_GENOME_MAP_COLOR_ORDER,
        )
        stack_totals = np.zeros(n_bins, dtype=float)
        per_cat_bins: dict[str, np.ndarray] = {}
        for cat in categories:
            part = hits[cat_series == cat]
            part_div = ((100.0 - part["pident"].astype(float)) / 100.0).clip(lower=0.0, upper=1.0)
            part_bp = part["te_bp"].astype(float).clip(lower=0.0)
            bp_per_bin, _ = np.histogram(
                part_div.to_numpy(dtype=float),
                bins=edges,
                weights=part_bp.to_numpy(dtype=float),
            )
            per_cat_bins[cat] = bp_per_bin
            stack_totals += bp_per_bin
        active = stack_totals > 0
        if not np.any(active):
            fig = go.Figure()
            fig.update_layout(
                title=f"TE divergence from {source} — no hits in divergence bins",
                xaxis={"title": "Divergence", "range": [0, 1]},
                yaxis={"title": "Total TE length (bp)"},
                template="plotly_white",
            )
        else:
            active_idx = np.where(active)[0]
            plot_centers = centers[active]
            plot_edges_lo = edges[active_idx]
            plot_edges_hi = edges[active_idx + 1]
            bar_width = bw
            global_stack_max = float(stack_totals[active].max()) if stack_totals[active].size else 1.0
            y_top = max(1.0, global_stack_max * 1.08)
            fig = go.Figure()
            stack_base = np.zeros(len(plot_centers), dtype=float)
            for cat in categories:
                bp_active = per_cat_bins[cat][active]
                bar_base = stack_base.copy()
                hover = [
                    (
                        f"Divergence: {lo:.3f} - {hi:.3f}<br>"
                        f"TE class: {cat}<br>"
                        f"Total TE length: {int(bp):,} bp"
                    )
                    for lo, hi, bp in zip(plot_edges_lo, plot_edges_hi, bp_active)
                ]
                fig.add_trace(
                    go.Bar(
                        x=plot_centers,
                        y=bp_active,
                        base=bar_base,
                        width=[bar_width] * len(plot_centers),
                        marker={
                            "color": cat_colors.get(cat, TE_EXTRA_TYPE_COLORS[0]),
                            "line": {"width": 0},
                        },
                        hovertext=hover,
                        hoverinfo="text",
                        name=cat,
                        legendgroup=cat,
                        showlegend=True,
                    )
                )
                stack_base += bp_active
            label = query_label or (
                te_hit_table_stem(primary)
                if len(te_paths) == 1
                else f"{te_hit_table_stem(primary)} + {len(te_paths) - 1} more"
            )
            fig.update_layout(
                title=f"TE divergence from {source} — {label}",
                showlegend=True,
                legend={
                    "title": {"text": ""},
                    "orientation": "v",
                    "yanchor": "top",
                    "y": 1.0,
                    "xanchor": "left",
                    "x": 1.02,
                    "bgcolor": "rgba(255,255,255,0.92)",
                    "bordercolor": "#d0d7de",
                    "borderwidth": 1,
                    "font": {"size": 11, "family": "Arial, sans-serif"},
                    "itemsizing": "constant",
                    "traceorder": "normal",
                },
                xaxis={
                    "title": "Divergence",
                    "range": [0.0, 1.0],
                    "autorange": False,
                    "showgrid": True,
                    "gridcolor": "#f0f0f0",
                    "zeroline": False,
                },
                yaxis={
                    "title": "Total TE length (bp)",
                    "range": [0, y_top],
                    "autorange": False,
                    "showgrid": True,
                    "gridcolor": "#f4f4f4",
                    "zeroline": False,
                },
                barmode="overlay",
                dragmode="zoom",
                uirevision="te-divergence-plot",
                template="plotly_white",
                bargap=0,
                margin={"l": 72, "r": 160, "t": 96, "b": 60},
                width=1120,
                height=540,
                meta={
                    "teDivergenceScale": {
                        "stackMax": global_stack_max,
                        "yMax": y_top,
                        "xRange": [0.0, 1.0],
                        "binWidth": bw,
                        "categoryColors": {
                            cat: cat_colors.get(cat, TE_EXTRA_TYPE_COLORS[0])
                            for cat in categories
                        },
                        "categoryOrder": list(categories),
                    }
                },
            )

    html = _figure_to_preview_html(fig, output_html)
    html = _inject_plotly_controls(
        html,
        scale_label="Y axis zoom",
        panel_mode="divergence",
    )
    output_html.write_text(html, encoding="utf-8")
    return output_html


def write_blast_genome_map_html(
    te_input_files: list[Path],
    output_html: Path | None = None,
    *,
    gene_gff: Path | list[Path] | None = None,
    methylation_files: Path | list[Path] | None = None,
    methylation_coverage_min: int = 5,
    genome_fasta: Path | None = None,
    bin_size: int = DEFAULT_BIN_SIZE,
    density_mode: str = TE_GENOME_MAP_DENSITY_HIT_COUNT,
    color_by: str = TE_GENOME_MAP_COLOR_SUPERFAMILY,
    query_label: str | None = None,
    region: tuple[str, float, float] | None = None,
    chromosomes_only: list[str] | None = None,
    superfamilies_only: list[str] | None = None,
    chromosome_sort: str = DEFAULT_CHROMOSOME_SORT,
    split_strand_tracks: bool = False,
    collapse_overlaps: bool = False,
    show_te_diff_track: bool = False,
    show_meth_diff_track: bool = False,
    show_window_ref_line: bool = True,
    annotation_density_mode: str = TE_GENOME_MAP_DENSITY_TE_LENGTH,
    annotation_feature_types: list[str] | None = None,
    annotation_feature_types_by_file: dict[str, list[str]] | None = None,
    te_input_kinds: list[str] | None = None,
    num_workers: int | None = None,
    cancel: JobCancellation | None = None,
    progress: Callable[[int, str], None] | None = None,
) -> Path | None:
    def _emit_progress(step: int, message: str) -> None:
        if progress is not None:
            progress(step, message)

    check_cancelled(cancel)
    _emit_progress(0, "Preparing TE genome map inputs…")
    te_paths = [p.expanduser().resolve() for p in te_input_files if p]
    gene_paths: list[Path] = []
    if gene_gff is not None:
        if isinstance(gene_gff, Path):
            gene_paths = [gene_gff.expanduser().resolve()]
        else:
            gene_paths = [p.expanduser().resolve() for p in gene_gff if p]
    methylation_paths: list[Path] = []
    if methylation_files is not None:
        if isinstance(methylation_files, Path):
            methylation_paths = [methylation_files.expanduser().resolve()]
        else:
            methylation_paths = [p.expanduser().resolve() for p in methylation_files if p]
    if not te_paths and not gene_paths and not methylation_paths:
        logger.warning("TE genome map: no TE, gene, or methylation input files")
        return None
    primary = (
        te_paths[0]
        if te_paths
        else (gene_paths[0] if gene_paths else methylation_paths[0])
    )
    if output_html is None:
        if te_paths:
            output_html = genome_map_html_path_for_te_inputs(
                te_paths,
                region=region is not None,
                gene_gff=gene_paths[0] if gene_paths else None,
            )
        elif gene_paths:
            output_html = genome_map_html_path_for_te_inputs(
                [],
                region=region is not None,
                gene_gff=gene_paths[0],
            )
        elif methylation_paths:
            stem = methylation_paths[0].stem
            suffix = GENOME_MAP_REGION_HTML_SUFFIX if region is not None else GENOME_MAP_HTML_SUFFIX
            output_html = methylation_paths[0].parent / f"{stem}{suffix}"
    else:
        output_html = output_html.expanduser().resolve()

    t0 = time.perf_counter()
    _emit_progress(1, "Reading TE hit table(s)…")
    te_hits_by_source: dict[str, pd.DataFrame] = {}
    te_source_specs: list[_TeSourceSpec] = []
    multi_te_files = len(te_paths) > 1
    te_labels = _disambiguate_path_labels(
        te_paths,
        label_for_path=lambda p: p.name if multi_te_files else te_hit_table_stem(p),
    )
    for index, path in enumerate(te_paths):
        check_cancelled(cancel)
        input_kind = TE_INPUT_KIND_MASKING
        if te_input_kinds and index < len(te_input_kinds):
            raw_kind = str(te_input_kinds[index]).strip().lower()
            if raw_kind == TE_INPUT_KIND_CONSENSUS:
                input_kind = TE_INPUT_KIND_CONSENSUS
        try:
            spec = te_source_spec_for_path(
                path,
                index,
                multi_file=multi_te_files,
                split_strand=bool(split_strand_tracks) and is_repeatmasker_out_file(path),
                label=te_labels[index],
                input_kind=input_kind,
            )
            te_hits_by_source[spec.key] = load_prepared_te_hits(path)
        except JobCancelledError:
            raise
        except Exception as exc:  # noqa: BLE001
            logger.warning("TE genome map: cannot read %s: %s", path, exc)
            return None
        te_source_specs.append(spec)
    hits = pd.DataFrame()
    t_read = time.perf_counter()

    genome_sizes: dict[str, int] = {}
    _emit_progress(2, "Reading annotation/methylation tracks…")
    if genome_fasta is not None:
        genome_path = genome_fasta.expanduser().resolve()
        if genome_path.is_file():
            _, genome_sizes = _load_genome_meta_from_fasta(genome_path)
    if genome_sizes and te_hits_by_source:
        allowed = set(genome_sizes.keys())
        for key in list(te_hits_by_source.keys()):
            df = harmonize_hit_chromosomes(te_hits_by_source[key], genome_sizes)
            if not df.empty and "chrom" in df.columns:
                df = df[df["chrom"].astype(str).isin(allowed)]
            te_hits_by_source[key] = df
        hits = _concat_dataframes(
            [df for df in te_hits_by_source.values() if not df.empty],
        )
    elif te_hits_by_source:
        hits = _concat_dataframes(
            [df for df in te_hits_by_source.values() if not df.empty],
        )
    chrom_sizes = chromosome_sizes_from_genome(genome_fasta, hits)
    if not chrom_sizes and genome_sizes:
        chrom_sizes = dict(genome_sizes)
    gene_hits = pd.DataFrame()
    gene_hits_by_source: dict[str, pd.DataFrame] = {}
    gene_source_specs: list[_GeneSourceSpec] = []
    multi_gene_files = len(gene_paths) > 1
    if gene_paths:
        gene_labels = _disambiguate_path_labels(
            gene_paths,
            label_for_path=lambda p: p.name if multi_gene_files else p.stem,
        )
        for index, gff_path in enumerate(gene_paths):
            check_cancelled(cancel)
            spec = gene_source_spec_for_path(
                gff_path,
                index,
                multi_file=multi_gene_files,
                label=gene_labels[index],
            )
            src_hits = gene_hits_from_gff(
                gff_path,
                chrom_sizes,
                feature_types=_annotation_feature_types_for_path(
                    gff_path,
                    annotation_feature_types=annotation_feature_types,
                    annotation_feature_types_by_file=annotation_feature_types_by_file,
                ),
            )
            if src_hits.empty:
                logger.warning("TE genome map: no gene features mapped from %s", gff_path)
            gene_hits_by_source[spec.key] = src_hits
            gene_source_specs.append(spec)
        gene_hits = _concat_dataframes(
            [df for df in gene_hits_by_source.values() if not df.empty],
        )
        if hits.empty and not gene_hits.empty:
            chrom_sizes = chromosome_sizes_from_genome(genome_fasta, gene_hits) or chrom_sizes
    methylation_hits = pd.DataFrame()
    methylation_hits_by_source: dict[str, pd.DataFrame] = {}
    methylation_cov_by_chrom_source: dict[str, dict[str, MethylationChromCov]] = {}
    methylation_source_specs: list[_MethylationSourceSpec] = []
    multi_meth_files = len(methylation_paths) > 1
    cov_min = max(1, int(methylation_coverage_min))
    if methylation_paths:
        meth_labels = _disambiguate_path_labels(
            methylation_paths,
            label_for_path=lambda p: _compact_bismark_cov_label(
                p.name if multi_meth_files else p.stem,
                max_len=28,
            ),
        )

        def _load_methylation_source(index: int, meth_path: Path) -> tuple[_MethylationSourceSpec, dict[str, MethylationChromCov]]:
            check_cancelled(cancel)
            if not is_bismark_cov_file(meth_path):
                logger.warning("TE genome map: not a Bismark coverage file: %s", meth_path)
                return (
                    methylation_source_spec_for_path(
                        meth_path,
                        index,
                        multi_file=multi_meth_files,
                        label=meth_labels[index],
                    ),
                    {},
                )
            spec = methylation_source_spec_for_path(
                meth_path,
                index,
                multi_file=multi_meth_files,
                label=meth_labels[index],
            )
            by_chrom = load_methylation_cov_by_chrom(
                meth_path,
                chrom_sizes,
                coverage_min=cov_min,
            )
            if not by_chrom:
                logger.warning("TE genome map: no methylation sites mapped from %s", meth_path)
            return spec, by_chrom

        meth_tasks = list(enumerate(methylation_paths))
        if len(meth_tasks) == 1:
            index, meth_path = meth_tasks[0]
            _emit_progress(2, f"Loading methylation: {meth_path.name}…")
            loaded = [_load_methylation_source(index, meth_path)]
        else:
            _emit_progress(2, f"Loading {len(meth_tasks)} methylation files…")
            max_io_workers = min(max(1, int(num_workers or (os.cpu_count() or 4))), len(meth_tasks))
            with ThreadPoolExecutor(max_workers=max_io_workers) as pool:
                loaded = list(pool.map(lambda item: _load_methylation_source(*item), meth_tasks))
        for spec, by_chrom in loaded:
            methylation_cov_by_chrom_source[spec.key] = by_chrom
            methylation_source_specs.append(spec)
        methylation_site_total = sum(
            methylation_cov_site_count(cov)
            for src in methylation_cov_by_chrom_source.values()
            for cov in src.values()
        )
        if methylation_site_total > 0:
            methylation_hits = pd.DataFrame({"chrom": ["__meth__"], "pos": [1]})
        if hits.empty and gene_hits.empty and methylation_cov_by_chrom_source:
            pseudo_chroms = sorted(
                {
                    chrom
                    for src in methylation_cov_by_chrom_source.values()
                    for chrom in src.keys()
                }
            )
            if pseudo_chroms:
                chrom_sizes = (
                    chromosome_sizes_from_genome(
                        genome_fasta,
                        pd.DataFrame({"chrom": pseudo_chroms, "pos": [1] * len(pseudo_chroms)}),
                    )
                    or chrom_sizes
                )
    has_methylation_data = any(
        methylation_cov_site_count(cov) > 0
        for src in methylation_cov_by_chrom_source.values()
        for cov in src.values()
    )
    if hits.empty and gene_hits.empty and not has_methylation_data:
        selected_chroms: set[str] = set()
        if chromosomes_only and chrom_sizes:
            selected_chroms = _map_chromosome_filter_to_genome(
                chromosomes_only,
                list(chrom_sizes.keys()),
            )
        has_track_sources = bool(
            te_source_specs or gene_source_specs or methylation_source_specs
        )
        if selected_chroms and has_track_sources and chrom_sizes:
            logger.info(
                "TE genome map: no features in inputs; showing %d selected "
                "chromosome(s) as empty tracks",
                len(selected_chroms),
            )
        else:
            message = _nothing_to_plot_message(
                te_paths=te_paths,
                gene_paths=gene_paths,
                methylation_paths=methylation_paths,
                chromosomes_only=chromosomes_only,
                chrom_sizes=chrom_sizes,
            )
            logger.info("TE genome map: %s", message.replace("\n\n", " "))
            raise GenomeMapNothingToPlotError(message)
    if query_label:
        label = query_label
    elif te_paths:
        if len(te_paths) == 1:
            label = te_hit_table_stem(te_paths[0])
        else:
            label = f"{te_hit_table_stem(te_paths[0])} + {len(te_paths) - 1} more"
    elif gene_paths:
        label = gene_paths[0].stem if len(gene_paths) == 1 else f"{gene_paths[0].stem} + {len(gene_paths) - 1} more"
    elif methylation_paths:
        label = (
            methylation_paths[0].stem
            if len(methylation_paths) == 1
            else f"{methylation_paths[0].stem} + {len(methylation_paths) - 1} more"
        )
    else:
        label = "genome map"
    source = te_hit_plot_source_label(te_paths[0]) if te_paths else "GFF/GTF"
    logger.info(
        "TE genome map: %d TE input file(s) -> %d track source(s)",
        len(te_paths),
        len(te_source_specs),
    )
    t_prep = time.perf_counter()
    check_cancelled(cancel)
    _emit_progress(3, "Building chromosome tracks…")
    fig = build_blast_genome_map_figure(
        hits,
        chrom_sizes,
        bin_size=bin_size,
        density_mode=density_mode,
        color_by=color_by,
        query_label=label,
        region=region,
        chromosomes_only=chromosomes_only,
        superfamilies_only=superfamilies_only,
        chromosome_sort=chromosome_sort,
        collapse_overlaps=collapse_overlaps,
        show_te_diff_track=show_te_diff_track,
        show_meth_diff_track=show_meth_diff_track,
        show_window_ref_line=show_window_ref_line,
        annotation_density_mode=annotation_density_mode,
        gene_hits=gene_hits if not gene_hits.empty else None,
        gene_sources=gene_source_specs or None,
        gene_hits_by_source=gene_hits_by_source or None,
        methylation_hits=methylation_hits if not methylation_hits.empty else None,
        methylation_sources=methylation_source_specs or None,
        methylation_hits_by_source=methylation_hits_by_source or None,
        methylation_cov_by_chrom_source=methylation_cov_by_chrom_source or None,
        methylation_coverage_min=cov_min,
        genome_fasta=genome_fasta,
        te_sources=te_source_specs,
        te_hits_by_source=te_hits_by_source,
        plot_source=source,
        split_strand_tracks=split_strand_tracks,
        num_workers=num_workers,
        cancel=cancel,
        progress=_emit_progress,
    )
    t_fig = time.perf_counter()
    _emit_progress(4, "Applying saved view settings…")
    view_state = load_view_state(output_html)
    scale_meta = _figure_genome_scale_meta(fig)
    merge_view_state_into_scale_meta(scale_meta, view_state)
    fig.update_layout(meta={"teGenomeScale": scale_meta})
    if view_state.x_range and len(view_state.x_range) >= 2:
        fig.update_layout(xaxis={"range": view_state.x_range, "autorange": False})
    html = _figure_to_preview_html(fig, output_html)
    html = _inject_plotly_controls(
        html,
        scale_label="Genome track bar scale",
        panel_mode="genome",
        prefs_path=sidecar_prefs_path(output_html),
    )
    _emit_progress(5, "Writing HTML output…")
    output_html.write_text(html, encoding="utf-8")
    _emit_progress(6, "TE genome map ready")
    t_done = time.perf_counter()
    logger.info(
        "TE genome map %s: read %.1fs, prep %.1fs, figure %.1fs, html %.1fs (total %.1fs, %d hits)",
        output_html.name,
        t_read - t0,
        t_prep - t_read,
        t_fig - t_prep,
        t_done - t_fig,
        t_done - t0,
        len(hits),
    )
    return output_html
