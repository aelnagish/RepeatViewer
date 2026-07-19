"""File-based TE genome map view preferences (global defaults + per-HTML sidecar)."""

from __future__ import annotations

import json
import logging
import os
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

PREFS_VERSION = 1

# Defaults (used when no JSON file exists yet).
DEFAULT_BAR_SCALE = 1.0
DEFAULT_ROW_SPACING = 1.0
DEFAULT_SHOW_WINDOW_REF_LINE = True
DEFAULT_WINDOW_REF_COLOR = "#b45309"
DEFAULT_PLOT_SCALE_MIN = 0.5
DEFAULT_PLOT_SCALE_MAX = 100.0
DEFAULT_PLOT_SCALE_MAX_GENOME = 15.0
DEFAULT_PLOT_SCALE_MAX_DIVERGENCE = 10.0
DEFAULT_PLOT_SCALE_STEP = 0.5
DEFAULT_ROW_SPACING_MIN = 0.6
DEFAULT_ROW_SPACING_MAX = 3.0
DEFAULT_ROW_SPACING_STEP = 0.1
DEFAULT_DIFF_TRACK_SCALE = 1.0
DEFAULT_DIFF_TRACK_SCALE_MIN = 0.2
DEFAULT_DIFF_TRACK_SCALE_MAX = 5.0
DEFAULT_DIFF_TRACK_SCALE_STEP = 0.1
DEFAULT_DIFF_COLOR_A = "#c0392b"
DEFAULT_DIFF_COLOR_B = "#2980b9"
DEFAULT_METH_DIFF_TRACK_SCALE = 1.0
DEFAULT_METH_DIFF_COLOR_A = "#8e44ad"
DEFAULT_METH_DIFF_COLOR_B = "#16a085"
DEFAULT_METH_COLOR_CG_TOTAL = "#7f7f7f"
DEFAULT_METH_COLOR_CG_COVERED = "#1f77b4"
DEFAULT_METH_COLOR_CG_COVERAGE_PCT = "#2ca02c"
DEFAULT_METH_COLOR_METH_LEVEL = "#d62728"
DEFAULT_METH_SHOW_CG_TOTAL = False
DEFAULT_METH_SHOW_CG_COVERED = False


@dataclass
class TeGenomeMapViewLimits:
    """Slider limits for the interactive genome map panel."""

    plot_scale_min: float = DEFAULT_PLOT_SCALE_MIN
    plot_scale_max: float = DEFAULT_PLOT_SCALE_MAX
    plot_scale_max_genome: float = DEFAULT_PLOT_SCALE_MAX_GENOME
    plot_scale_max_divergence: float = DEFAULT_PLOT_SCALE_MAX_DIVERGENCE
    plot_scale_step: float = DEFAULT_PLOT_SCALE_STEP
    row_spacing_min: float = DEFAULT_ROW_SPACING_MIN
    row_spacing_max: float = DEFAULT_ROW_SPACING_MAX
    row_spacing_step: float = DEFAULT_ROW_SPACING_STEP
    diff_track_scale_min: float = DEFAULT_DIFF_TRACK_SCALE_MIN
    diff_track_scale_max: float = DEFAULT_DIFF_TRACK_SCALE_MAX
    diff_track_scale_step: float = DEFAULT_DIFF_TRACK_SCALE_STEP


@dataclass
class TeGenomeMapViewState:
    """Per-map interactive view state saved beside the HTML file."""

    bar_scale: float = DEFAULT_BAR_SCALE
    row_spacing: float = DEFAULT_ROW_SPACING
    diff_track_scale: float = DEFAULT_DIFF_TRACK_SCALE
    diff_color_a: str = DEFAULT_DIFF_COLOR_A
    diff_color_b: str = DEFAULT_DIFF_COLOR_B
    show_window_ref_line: bool = DEFAULT_SHOW_WINDOW_REF_LINE
    window_ref_color: str = DEFAULT_WINDOW_REF_COLOR
    x_range: list[float] | None = None
    category_order: list[str] = field(default_factory=list)
    category_colors: dict[str, str] = field(default_factory=dict)
    category_hidden: list[str] = field(default_factory=list)
    meth_diff_track_scale: float = DEFAULT_METH_DIFF_TRACK_SCALE
    meth_diff_color_a: str = DEFAULT_METH_DIFF_COLOR_A
    meth_diff_color_b: str = DEFAULT_METH_DIFF_COLOR_B
    meth_show_cg_total: bool = DEFAULT_METH_SHOW_CG_TOTAL
    meth_show_cg_covered: bool = DEFAULT_METH_SHOW_CG_COVERED
    meth_show_cg_pct: bool = True
    meth_show_level: bool = True
    meth_show_diff: bool = True
    meth_color_cg_total: str = DEFAULT_METH_COLOR_CG_TOTAL
    meth_color_cg_covered: str = DEFAULT_METH_COLOR_CG_COVERED
    meth_color_cg_pct: str = DEFAULT_METH_COLOR_CG_COVERAGE_PCT
    meth_color_level: str = DEFAULT_METH_COLOR_METH_LEVEL
    gene_source_hidden: list[str] = field(default_factory=list)
    te_source_hidden: list[str] = field(default_factory=list)
    track_group_order: list[str] = field(default_factory=list)
    track_group_hidden: list[str] = field(default_factory=list)
    track_ymax_by_group: dict[str, Any] = field(default_factory=dict)
    chromosome_sort: str = "natural"


def _config_dir() -> Path:
    base = os.environ.get("XDG_CONFIG_HOME", str(Path.home() / ".config"))
    return Path(base) / "TEviewer"


def global_view_limits_path() -> Path:
    return _config_dir() / "te_genome_map_view.json"


def sidecar_prefs_path(html_path: Path) -> Path:
    html_path = html_path.expanduser().resolve()
    return html_path.with_name(f"{html_path.stem}.prefs.json")


def _coerce_float(value: Any, default: float) -> float:
    try:
        return float(value)
    except (TypeError, ValueError):
        return default


def _coerce_bool(value: Any, default: bool) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return default
    return str(value).lower() in {"1", "true", "yes"}


def load_view_limits() -> TeGenomeMapViewLimits:
    path = global_view_limits_path()
    if not path.is_file():
        return TeGenomeMapViewLimits()
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("Could not read %s: %s", path, exc)
        return TeGenomeMapViewLimits()
    if not isinstance(raw, dict):
        return TeGenomeMapViewLimits()
    limits = raw.get("limits", raw)
    if not isinstance(limits, dict):
        limits = {}
    return TeGenomeMapViewLimits(
        plot_scale_min=_coerce_float(limits.get("plot_scale_min"), DEFAULT_PLOT_SCALE_MIN),
        plot_scale_max=_coerce_float(limits.get("plot_scale_max"), DEFAULT_PLOT_SCALE_MAX),
        plot_scale_max_genome=_coerce_float(
            limits.get("plot_scale_max_genome"), DEFAULT_PLOT_SCALE_MAX_GENOME
        ),
        plot_scale_max_divergence=_coerce_float(
            limits.get("plot_scale_max_divergence"), DEFAULT_PLOT_SCALE_MAX_DIVERGENCE
        ),
        plot_scale_step=_coerce_float(limits.get("plot_scale_step"), DEFAULT_PLOT_SCALE_STEP),
        row_spacing_min=_coerce_float(limits.get("row_spacing_min"), DEFAULT_ROW_SPACING_MIN),
        row_spacing_max=_coerce_float(limits.get("row_spacing_max"), DEFAULT_ROW_SPACING_MAX),
        row_spacing_step=_coerce_float(limits.get("row_spacing_step"), DEFAULT_ROW_SPACING_STEP),
        diff_track_scale_min=_coerce_float(
            limits.get("diff_track_scale_min"), DEFAULT_DIFF_TRACK_SCALE_MIN
        ),
        diff_track_scale_max=_coerce_float(
            limits.get("diff_track_scale_max"), DEFAULT_DIFF_TRACK_SCALE_MAX
        ),
        diff_track_scale_step=_coerce_float(
            limits.get("diff_track_scale_step"), DEFAULT_DIFF_TRACK_SCALE_STEP
        ),
    )


def save_view_limits(limits: TeGenomeMapViewLimits) -> None:
    path = global_view_limits_path()
    path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "version": PREFS_VERSION,
        "limits": asdict(limits),
    }
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _parse_x_range(raw: Any) -> list[float] | None:
    if not isinstance(raw, (list, tuple)) or len(raw) < 2:
        return None
    try:
        r0 = float(raw[0])
        r1 = float(raw[1])
    except (TypeError, ValueError):
        return None
    if r1 <= r0:
        return None
    return [r0, r1]


def _parse_str_list(raw: Any) -> list[str]:
    if not isinstance(raw, list):
        return []
    return [str(x).strip() for x in raw if str(x).strip()]


def _parse_color_map(raw: Any) -> dict[str, str]:
    if not isinstance(raw, dict):
        return {}
    out: dict[str, str] = {}
    for key, value in raw.items():
        name = str(key).strip()
        color = str(value).strip()
        if name and color:
            out[name] = color
    return out


def _parse_track_ymax_map(raw: Any) -> dict[str, Any]:
    if not isinstance(raw, dict):
        return {}
    out: dict[str, Any] = {}
    for key, value in raw.items():
        group = str(key).strip()
        if not group:
            continue
        if isinstance(value, dict):
            sub_out: dict[str, float] = {}
            for sub_key, sub_value in value.items():
                sub_name = str(sub_key).strip()
                if not sub_name:
                    continue
                try:
                    num = float(sub_value)
                except (TypeError, ValueError):
                    continue
                if num > 0:
                    sub_out[sub_name] = num
            if sub_out:
                out[group] = sub_out
            continue
        try:
            num = float(value)
        except (TypeError, ValueError):
            continue
        if num > 0:
            out[group] = num
    return out


def load_view_state(html_path: Path | None = None) -> TeGenomeMapViewState:
    if html_path is None:
        return TeGenomeMapViewState()
    path = sidecar_prefs_path(html_path)
    if not path.is_file():
        return TeGenomeMapViewState()
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        logger.warning("Could not read %s: %s", path, exc)
        return TeGenomeMapViewState()
    if not isinstance(raw, dict):
        return TeGenomeMapViewState()
    state = raw.get("state", raw)
    if not isinstance(state, dict):
        state = {}
    return TeGenomeMapViewState(
        bar_scale=_coerce_float(state.get("bar_scale"), DEFAULT_BAR_SCALE),
        row_spacing=_coerce_float(state.get("row_spacing"), DEFAULT_ROW_SPACING),
        diff_track_scale=_coerce_float(state.get("diff_track_scale"), DEFAULT_DIFF_TRACK_SCALE),
        diff_color_a=str(state.get("diff_color_a") or DEFAULT_DIFF_COLOR_A),
        diff_color_b=str(state.get("diff_color_b") or DEFAULT_DIFF_COLOR_B),
        show_window_ref_line=_coerce_bool(
            state.get("show_window_ref_line"), DEFAULT_SHOW_WINDOW_REF_LINE
        ),
        window_ref_color=str(state.get("window_ref_color") or DEFAULT_WINDOW_REF_COLOR),
        x_range=_parse_x_range(state.get("x_range")),
        category_order=_parse_str_list(state.get("category_order")),
        category_colors=_parse_color_map(state.get("category_colors")),
        category_hidden=_parse_str_list(state.get("category_hidden")),
        meth_diff_track_scale=_coerce_float(
            state.get("meth_diff_track_scale"), DEFAULT_METH_DIFF_TRACK_SCALE
        ),
        meth_diff_color_a=str(state.get("meth_diff_color_a") or DEFAULT_METH_DIFF_COLOR_A),
        meth_diff_color_b=str(state.get("meth_diff_color_b") or DEFAULT_METH_DIFF_COLOR_B),
        meth_show_cg_total=_coerce_bool(
            state.get("meth_show_cg_total"), DEFAULT_METH_SHOW_CG_TOTAL
        ),
        meth_show_cg_covered=_coerce_bool(
            state.get("meth_show_cg_covered"), DEFAULT_METH_SHOW_CG_COVERED
        ),
        meth_show_cg_pct=_coerce_bool(state.get("meth_show_cg_pct"), True),
        meth_show_level=_coerce_bool(state.get("meth_show_level"), True),
        meth_show_diff=_coerce_bool(state.get("meth_show_diff"), True),
        meth_color_cg_total=str(state.get("meth_color_cg_total") or DEFAULT_METH_COLOR_CG_TOTAL),
        meth_color_cg_covered=str(
            state.get("meth_color_cg_covered") or DEFAULT_METH_COLOR_CG_COVERED
        ),
        meth_color_cg_pct=str(state.get("meth_color_cg_pct") or DEFAULT_METH_COLOR_CG_COVERAGE_PCT),
        meth_color_level=str(state.get("meth_color_level") or DEFAULT_METH_COLOR_METH_LEVEL),
        gene_source_hidden=_parse_str_list(state.get("gene_source_hidden")),
        te_source_hidden=_parse_str_list(state.get("te_source_hidden")),
        track_group_order=_parse_str_list(state.get("track_group_order")),
        track_group_hidden=_parse_str_list(state.get("track_group_hidden")),
        track_ymax_by_group=_parse_track_ymax_map(state.get("track_ymax_by_group")),
        chromosome_sort=str(state.get("chromosome_sort") or "natural"),
    )


def save_view_state(html_path: Path, state: TeGenomeMapViewState) -> None:
    path = sidecar_prefs_path(html_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    payload: dict[str, Any] = {
        "version": PREFS_VERSION,
        "state": {
            "bar_scale": state.bar_scale,
            "row_spacing": state.row_spacing,
            "diff_track_scale": state.diff_track_scale,
            "diff_color_a": state.diff_color_a,
            "diff_color_b": state.diff_color_b,
            "show_window_ref_line": state.show_window_ref_line,
            "window_ref_color": state.window_ref_color,
            "category_order": list(state.category_order),
            "category_colors": dict(state.category_colors),
            "category_hidden": list(state.category_hidden),
            "meth_diff_track_scale": state.meth_diff_track_scale,
            "meth_diff_color_a": state.meth_diff_color_a,
            "meth_diff_color_b": state.meth_diff_color_b,
            "meth_show_cg_total": state.meth_show_cg_total,
            "meth_show_cg_covered": state.meth_show_cg_covered,
            "meth_show_cg_pct": state.meth_show_cg_pct,
            "meth_show_level": state.meth_show_level,
            "meth_show_diff": state.meth_show_diff,
            "meth_color_cg_total": state.meth_color_cg_total,
            "meth_color_cg_covered": state.meth_color_cg_covered,
            "meth_color_cg_pct": state.meth_color_cg_pct,
            "meth_color_level": state.meth_color_level,
            "gene_source_hidden": list(state.gene_source_hidden),
            "te_source_hidden": list(state.te_source_hidden),
            "track_group_order": list(state.track_group_order),
            "track_group_hidden": list(state.track_group_hidden),
            "track_ymax_by_group": dict(state.track_ymax_by_group),
            "chromosome_sort": state.chromosome_sort,
        },
    }
    if state.x_range and len(state.x_range) >= 2:
        payload["state"]["x_range"] = [float(state.x_range[0]), float(state.x_range[1])]
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def view_state_from_dict(data: dict[str, Any]) -> TeGenomeMapViewState:
    return TeGenomeMapViewState(
        bar_scale=_coerce_float(data.get("bar_scale"), DEFAULT_BAR_SCALE),
        row_spacing=_coerce_float(data.get("row_spacing"), DEFAULT_ROW_SPACING),
        diff_track_scale=_coerce_float(data.get("diff_track_scale"), DEFAULT_DIFF_TRACK_SCALE),
        diff_color_a=str(data.get("diff_color_a") or DEFAULT_DIFF_COLOR_A),
        diff_color_b=str(data.get("diff_color_b") or DEFAULT_DIFF_COLOR_B),
        show_window_ref_line=_coerce_bool(
            data.get("show_window_ref_line"), DEFAULT_SHOW_WINDOW_REF_LINE
        ),
        window_ref_color=str(data.get("window_ref_color") or DEFAULT_WINDOW_REF_COLOR),
        x_range=_parse_x_range(data.get("x_range")),
        category_order=_parse_str_list(data.get("category_order")),
        category_colors=_parse_color_map(data.get("category_colors")),
        category_hidden=_parse_str_list(data.get("category_hidden")),
        meth_diff_track_scale=_coerce_float(
            data.get("meth_diff_track_scale"), DEFAULT_METH_DIFF_TRACK_SCALE
        ),
        meth_diff_color_a=str(data.get("meth_diff_color_a") or DEFAULT_METH_DIFF_COLOR_A),
        meth_diff_color_b=str(data.get("meth_diff_color_b") or DEFAULT_METH_DIFF_COLOR_B),
        meth_show_cg_total=_coerce_bool(
            data.get("meth_show_cg_total"), DEFAULT_METH_SHOW_CG_TOTAL
        ),
        meth_show_cg_covered=_coerce_bool(
            data.get("meth_show_cg_covered"), DEFAULT_METH_SHOW_CG_COVERED
        ),
        meth_show_cg_pct=_coerce_bool(data.get("meth_show_cg_pct"), True),
        meth_show_level=_coerce_bool(data.get("meth_show_level"), True),
        meth_show_diff=_coerce_bool(data.get("meth_show_diff"), True),
        meth_color_cg_total=str(data.get("meth_color_cg_total") or DEFAULT_METH_COLOR_CG_TOTAL),
        meth_color_cg_covered=str(
            data.get("meth_color_cg_covered") or DEFAULT_METH_COLOR_CG_COVERED
        ),
        meth_color_cg_pct=str(data.get("meth_color_cg_pct") or DEFAULT_METH_COLOR_CG_COVERAGE_PCT),
        meth_color_level=str(data.get("meth_color_level") or DEFAULT_METH_COLOR_METH_LEVEL),
        gene_source_hidden=_parse_str_list(data.get("gene_source_hidden")),
        te_source_hidden=_parse_str_list(data.get("te_source_hidden")),
        track_group_order=_parse_str_list(data.get("track_group_order")),
        track_group_hidden=_parse_str_list(data.get("track_group_hidden")),
        track_ymax_by_group=_parse_track_ymax_map(data.get("track_ymax_by_group")),
        chromosome_sort=str(data.get("chromosome_sort") or "natural"),
    )


def merge_view_state_into_scale_meta(scale_meta: dict[str, Any], state: TeGenomeMapViewState) -> None:
    """Apply saved view state into Plotly layout meta consumed by embedded JS."""
    if state.bar_scale > 0:
        scale_meta["savedBarScale"] = float(state.bar_scale)
    if state.row_spacing > 0:
        scale_meta["savedRowSpacing"] = float(state.row_spacing)
    if state.diff_track_scale > 0:
        scale_meta["savedDiffTrackScale"] = float(state.diff_track_scale)
    if state.diff_color_a:
        scale_meta["savedDiffColorA"] = str(state.diff_color_a)
        scale_meta["diffColorA"] = str(state.diff_color_a)
    if state.diff_color_b:
        scale_meta["savedDiffColorB"] = str(state.diff_color_b)
        scale_meta["diffColorB"] = str(state.diff_color_b)
    scale_meta["savedShowWindowRefLine"] = bool(state.show_window_ref_line)
    if state.window_ref_color:
        scale_meta["savedWindowRefColor"] = str(state.window_ref_color)
        scale_meta["windowRefColor"] = str(state.window_ref_color)
    if state.x_range and len(state.x_range) >= 2:
        scale_meta["savedXRange"] = [float(state.x_range[0]), float(state.x_range[1])]
    if state.category_order:
        scale_meta["savedCategoryOrder"] = list(state.category_order)
    if state.category_colors:
        scale_meta["savedCategoryColors"] = dict(state.category_colors)
    if state.category_hidden:
        scale_meta["savedCategoryHidden"] = list(state.category_hidden)
    if state.meth_diff_track_scale > 0:
        scale_meta["savedMethDiffTrackScale"] = float(state.meth_diff_track_scale)
    if state.meth_diff_color_a:
        scale_meta["savedMethDiffColorA"] = str(state.meth_diff_color_a)
        scale_meta["methDiffColorA"] = str(state.meth_diff_color_a)
    if state.meth_diff_color_b:
        scale_meta["savedMethDiffColorB"] = str(state.meth_diff_color_b)
        scale_meta["methDiffColorB"] = str(state.meth_diff_color_b)
    scale_meta["savedMethShowCgTotal"] = bool(state.meth_show_cg_total)
    scale_meta["savedMethShowCgCovered"] = bool(state.meth_show_cg_covered)
    scale_meta["savedMethShowCgPct"] = bool(state.meth_show_cg_pct)
    scale_meta["savedMethShowLevel"] = bool(state.meth_show_level)
    scale_meta["savedMethShowDiff"] = bool(state.meth_show_diff)
    if state.meth_color_cg_total:
        scale_meta["savedMethColorCgTotal"] = str(state.meth_color_cg_total)
        scale_meta["methColorCgTotal"] = str(state.meth_color_cg_total)
    if state.meth_color_cg_covered:
        scale_meta["savedMethColorCgCovered"] = str(state.meth_color_cg_covered)
        scale_meta["methColorCgCovered"] = str(state.meth_color_cg_covered)
    if state.meth_color_cg_pct:
        scale_meta["savedMethColorCgPct"] = str(state.meth_color_cg_pct)
        scale_meta["methColorCgPct"] = str(state.meth_color_cg_pct)
    if state.meth_color_level:
        scale_meta["savedMethColorLevel"] = str(state.meth_color_level)
        scale_meta["methColorLevel"] = str(state.meth_color_level)
    if state.gene_source_hidden:
        scale_meta["savedGeneSourceHidden"] = list(state.gene_source_hidden)
    if state.te_source_hidden:
        scale_meta["savedTeSourceHidden"] = list(state.te_source_hidden)
    if state.track_group_order:
        scale_meta["savedTrackGroupOrder"] = list(state.track_group_order)
    if state.track_group_hidden:
        scale_meta["savedTrackGroupHidden"] = list(state.track_group_hidden)
    if state.track_ymax_by_group:
        scale_meta["savedTrackYMaxByGroup"] = dict(state.track_ymax_by_group)
        scale_meta["trackYMaxByGroup"] = dict(state.track_ymax_by_group)
    if state.chromosome_sort:
        scale_meta["savedChromosomeSort"] = str(state.chromosome_sort)


def ensure_default_prefs_file() -> None:
    """Write global limits JSON once if missing (migrates defaults out of Python)."""
    path = global_view_limits_path()
    if path.is_file():
        return
    save_view_limits(TeGenomeMapViewLimits())
