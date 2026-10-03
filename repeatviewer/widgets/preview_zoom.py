"""Persisted zoom / font / fit settings for Preview (text, HTML, PDF, images)."""

from __future__ import annotations

from enum import Enum

from PySide6.QtCore import QSettings
from PySide6.QtGui import QWheelEvent

PREVIEW_SETTINGS_GROUP = "preview"
DEFAULT_PREVIEW_ZOOM_PERCENT = 100
MIN_PREVIEW_ZOOM_PERCENT = 25
MAX_PREVIEW_ZOOM_PERCENT = 500
PREVIEW_WHEEL_ZOOM_STEP_PERCENT = 2
PREVIEW_WHEEL_FONT_STEP_PT = 1
# Accumulated scroll units per zoom step. Classic mouse notch = 120.
PREVIEW_WHEEL_DELTA_PER_STEP = 80
DEFAULT_PREVIEW_FONT_SIZE = 11
MIN_PREVIEW_FONT_SIZE = 8
MAX_PREVIEW_FONT_SIZE = 72


class PreviewFitMode(str, Enum):
    """How images are scaled before applying zoom percent."""

    FIT_WINDOW = "fit_window"
    FIT_WIDTH = "fit_width"
    ORIGINAL = "original"


class DefaultFastaPreviewTab(str, Enum):
    """Which sub-tab to show when opening a single-sequence FASTA in Preview."""

    TEXT = "text"
    ANNOTATION = "annotation"
    DOTPLOT = "dotplot"


class PreviewZoomKind(str, Enum):
    """Preview content kind for chrome / persistence.

    HTML, PDF, and image each keep an independent zoom %. Text uses font size
    (pt), not a percent scale.
    """

    TEXT = "text"
    HTML = "html"
    PDF = "pdf"
    IMAGE = "image"


DEFAULT_FASTA_PREVIEW_TAB = DefaultFastaPreviewTab.ANNOTATION

_PREVIEW_ZOOM_KIND_KEYS: dict[PreviewZoomKind, str] = {
    PreviewZoomKind.HTML: "zoom_percent_html",
    PreviewZoomKind.PDF: "zoom_percent_pdf",
    PreviewZoomKind.IMAGE: "zoom_percent_image",
}


def _clamp_preview_zoom_percent(percent: int) -> int:
    return max(MIN_PREVIEW_ZOOM_PERCENT, min(MAX_PREVIEW_ZOOM_PERCENT, int(percent)))


def _coerce_preview_zoom_kind(
    kind: PreviewZoomKind | str | None,
) -> PreviewZoomKind | None:
    if kind is None:
        return None
    if isinstance(kind, PreviewZoomKind):
        return kind
    try:
        return PreviewZoomKind(str(kind))
    except ValueError:
        return PreviewZoomKind.TEXT


def load_preview_zoom_percent(kind: PreviewZoomKind | str | None = None) -> int:
    """Load zoom for a media preview kind (HTML / PDF / image).

    ``kind=None`` returns the legacy shared value. Text preview uses font size
    (see ``load_preview_font_size``), not a percent scale.
    """
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    legacy = _clamp_preview_zoom_percent(
        int(s.value("zoom_percent", DEFAULT_PREVIEW_ZOOM_PERCENT))
    )
    zoom_kind = _coerce_preview_zoom_kind(kind)
    if zoom_kind is None:
        s.endGroup()
        return legacy
    key = _PREVIEW_ZOOM_KIND_KEYS.get(zoom_kind)
    if key is None:
        s.endGroup()
        return legacy
    if s.contains(key):
        value = int(s.value(key, legacy))
    else:
        value = legacy
    s.endGroup()
    return _clamp_preview_zoom_percent(value)


def save_preview_zoom_percent(
    percent: int,
    kind: PreviewZoomKind | str | None = None,
) -> None:
    """Save zoom for one media preview kind (or the legacy shared key)."""
    percent = _clamp_preview_zoom_percent(percent)
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    zoom_kind = _coerce_preview_zoom_kind(kind)
    if zoom_kind is None:
        s.setValue("zoom_percent", percent)
    else:
        key = _PREVIEW_ZOOM_KIND_KEYS.get(zoom_kind)
        if key is not None:
            s.setValue(key, percent)
    s.endGroup()


def load_all_preview_zoom_percents() -> dict[PreviewZoomKind, int]:
    return {
        kind: load_preview_zoom_percent(kind)
        for kind in (PreviewZoomKind.HTML, PreviewZoomKind.PDF, PreviewZoomKind.IMAGE)
    }


def save_all_preview_zoom_percents(values: dict[PreviewZoomKind, int]) -> None:
    for kind, percent in values.items():
        if kind in _PREVIEW_ZOOM_KIND_KEYS:
            save_preview_zoom_percent(percent, kind)


def _clamp_preview_font_size(size: int) -> int:
    return max(MIN_PREVIEW_FONT_SIZE, min(MAX_PREVIEW_FONT_SIZE, int(size)))


def load_preview_fit_mode() -> PreviewFitMode:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    raw = str(s.value("fit_mode", PreviewFitMode.FIT_WINDOW.value))
    s.endGroup()
    try:
        return PreviewFitMode(raw)
    except ValueError:
        return PreviewFitMode.FIT_WINDOW


def wheel_zoom_delta_from_event(event: QWheelEvent) -> int:
    """Signed raw scroll delta for zoom (positive = zoom in).

    On macOS, Shift+two-finger scroll is remapped to a horizontal delta, so
    use the dominant axis (prefer Y, fall back to X). Trackpad events often
    only set ``pixelDelta``; scale those so they accumulate like angle notches.
    """
    angle = event.angleDelta()
    pixel = event.pixelDelta()
    ax, ay = int(angle.x()), int(angle.y())
    if ay or ax:
        return ay if abs(ay) >= abs(ax) else ax
    px, py = int(pixel.x()), int(pixel.y())
    if py or px:
        # ~30 px ≈ one classic 120-unit notch after scaling.
        scaled_y = py * 4
        scaled_x = px * 4
        return scaled_y if abs(scaled_y) >= abs(scaled_x) else scaled_x
    return 0


def wheel_zoom_steps_from_event(event: QWheelEvent) -> int:
    """Return signed zoom steps from a single wheel event (no accumulation).

    Prefer :func:`wheel_zoom_delta_from_event` with a running accumulator for
    trackpad gestures so the same drag distance does not jump too quickly.
    """
    delta = wheel_zoom_delta_from_event(event)
    if delta == 0:
        return 0
    steps = int(round(delta / float(PREVIEW_WHEEL_DELTA_PER_STEP)))
    if steps == 0:
        steps = 1 if delta > 0 else -1
    return steps


def save_preview_fit_mode(mode: PreviewFitMode) -> None:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    s.setValue("fit_mode", mode.value)
    s.endGroup()


def load_preview_font_size() -> int:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    value = int(s.value("text_font_size", DEFAULT_PREVIEW_FONT_SIZE))
    s.endGroup()
    return _clamp_preview_font_size(value)


def save_preview_font_size(size: int) -> None:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    s.setValue("text_font_size", _clamp_preview_font_size(size))
    s.endGroup()


def load_default_fasta_preview_tab() -> DefaultFastaPreviewTab:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    raw = str(s.value("default_fasta_tab", DEFAULT_FASTA_PREVIEW_TAB.value))
    s.endGroup()
    try:
        return DefaultFastaPreviewTab(raw)
    except ValueError:
        return DEFAULT_FASTA_PREVIEW_TAB


def save_default_fasta_preview_tab(tab: DefaultFastaPreviewTab) -> None:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    s.setValue("default_fasta_tab", tab.value)
    s.endGroup()


DEFAULT_SHOW_SEQUENCE_SUMMARY = True


def load_show_sequence_summary() -> bool:
    """Whether to show sequence length / GC% for single-sequence FASTA files.

    Computing GC means scanning every base, so on very large sequences (a
    500 Mbp scaffold) it costs a few seconds even after optimisation. Turning
    this off skips that work entirely.
    """
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    raw = s.value("show_sequence_summary", DEFAULT_SHOW_SEQUENCE_SUMMARY)
    s.endGroup()
    if isinstance(raw, bool):
        return raw
    return str(raw).strip().lower() in {"1", "true", "yes", "on"}


def save_show_sequence_summary(enabled: bool) -> None:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    s.setValue("show_sequence_summary", bool(enabled))
    s.endGroup()


# Wong / Paul Tol–style nucleotide colors (color-blind friendly).
_COLORBLIND_BASE_COLORS: dict[str, tuple[str, str]] = {
    "A": ("#E69F00", "#000000"),
    "C": ("#0072B2", "#FFFFFF"),
    "G": ("#000000", "#FFFFFF"),
    "T": ("#009E73", "#FFFFFF"),
    "U": ("#009E73", "#FFFFFF"),
    "N": ("#C8C8C8", "#000000"),
}

# Same MSA "Default" palette used by Element colors / MSA plot (figure 2).
_STANDARD_BASE_COLORS: dict[str, tuple[str, str]] = {
    "A": ("#00C800", "#000000"),
    "C": ("#7040E8", "#FFFFFF"),
    "G": ("#000000", "#FFFFFF"),
    "T": ("#F06090", "#FFFFFF"),
    "U": ("#F06090", "#FFFFFF"),
    "N": ("#C8C8C8", "#000000"),
}

# Pre–MSA-default AliView ACGT colours; ignore if still stored as "custom".
_LEGACY_ALIVIEW_ACGT: dict[str, str] = {
    "A": "#00FF00",
    "C": "#0000FF",
    "G": "#000000",
    "T": "#FF0000",
    "U": "#FF0000",
}


# Hot path (Annotation color bar) must not reopen QSettings per base.
_nucleotide_color_map_cache: dict[str, tuple[str, str]] | None = None
_amino_acid_color_map_cache: dict[str, tuple[str, str]] | None = None


def invalidate_residue_color_maps() -> None:
    """Drop cached nucleotide / amino-acid palettes (call after color settings change)."""
    global _nucleotide_color_map_cache, _amino_acid_color_map_cache
    _nucleotide_color_map_cache = None
    _amino_acid_color_map_cache = None


def load_color_blind_nucleotide_scheme() -> bool:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    enabled = bool(s.value("color_blind_nucleotide_scheme", False))
    s.endGroup()
    return enabled


def save_color_blind_nucleotide_scheme(enabled: bool) -> None:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    s.setValue("color_blind_nucleotide_scheme", bool(enabled))
    s.endGroup()
    invalidate_residue_color_maps()


def _scheme_base_colors() -> dict[str, tuple[str, str]]:
    if load_color_blind_nucleotide_scheme():
        return dict(_COLORBLIND_BASE_COLORS)
    return dict(_STANDARD_BASE_COLORS)


def _strip_legacy_aliview_overrides(colors: dict[str, str]) -> dict[str, str]:
    """Drop old AliView ACGT custom values so MSA Default is the true default."""
    cleaned = dict(colors)
    for base, legacy in _LEGACY_ALIVIEW_ACGT.items():
        if cleaned.get(base, "").upper() == legacy.upper():
            cleaned.pop(base, None)
    return cleaned


def load_custom_nucleotide_backgrounds() -> dict[str, str]:
    """Load custom nucleotide colours (full alphabet JSON).

    Read-only: do not write/remove settings here — the Annotation color bar
    calls into this path thousands of times when the palette is uncached.
    Legacy per-base keys are stripped in ``save_custom_nucleotide_backgrounds``.
    """
    import json

    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    raw = s.value("nucleotide_colors_json", "")
    s.endGroup()

    custom: dict[str, str] = {}
    text = str(raw or "").strip()
    if text:
        try:
            decoded = json.loads(text)
        except (json.JSONDecodeError, TypeError):
            decoded = None
        if isinstance(decoded, dict):
            for key, value in decoded.items():
                color = str(value).strip()
                if key and color.startswith("#") and len(color) >= 4:
                    custom[str(key).upper()] = color
    custom = _strip_legacy_aliview_overrides(custom)
    if "U" in custom:
        custom["T"] = custom["U"]
    elif "T" in custom:
        custom["U"] = custom["T"]
    return custom


def save_custom_nucleotide_backgrounds(colors: dict[str, str]) -> None:
    import json

    cleaned: dict[str, str] = {}
    for key, value in colors.items():
        color = str(value).strip()
        if key and color.startswith("#") and len(color) >= 4:
            cleaned[str(key).upper()] = color
    if "U" in cleaned:
        cleaned["T"] = cleaned["U"]
    elif "T" in cleaned:
        cleaned["U"] = cleaned["T"]

    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    for base in ("A", "C", "G", "T"):
        s.remove(f"nucleotide_{base}_bg")
    if cleaned:
        s.setValue("nucleotide_colors_json", json.dumps(cleaned, sort_keys=True))
    else:
        s.remove("nucleotide_colors_json")
    s.endGroup()
    invalidate_residue_color_maps()


def clear_custom_nucleotide_backgrounds() -> None:
    save_custom_nucleotide_backgrounds({})


def default_nucleotide_backgrounds() -> dict[str, str]:
    """MSA-style Default palette (figure 2 / Element colors Default preset)."""
    from repeatviewer.backends.msa_plot import DEFAULT_NUCLEOTIDE_COLORS

    return dict(DEFAULT_NUCLEOTIDE_COLORS)


def nucleotide_base_color_map() -> dict[str, tuple[str, str]]:
    """Return bg/fg hex colors for nucleotide letters in the Annotation view."""
    global _nucleotide_color_map_cache
    if _nucleotide_color_map_cache is not None:
        return _nucleotide_color_map_cache

    from repeatviewer.backends.msa_plot import (
        DEFAULT_NUCLEOTIDE_COLORS,
        NUCLEOTIDE_COLORS_COLORBLIND,
    )

    if load_color_blind_nucleotide_scheme():
        merged = dict(NUCLEOTIDE_COLORS_COLORBLIND)
    else:
        merged = dict(DEFAULT_NUCLEOTIDE_COLORS)
    merged.update(load_custom_nucleotide_backgrounds())
    if "U" in merged:
        merged["T"] = merged["U"]
    elif "T" in merged:
        merged["U"] = merged["T"]
    _nucleotide_color_map_cache = {
        base: (bg, _contrast_fg_for_bg(bg)) for base, bg in merged.items()
    }
    return _nucleotide_color_map_cache


def _contrast_fg_for_bg(bg_hex: str) -> str:
    from PySide6.QtGui import QColor

    from repeatviewer.widgets.residue_color_grid import ResidueColorGrid

    return ResidueColorGrid.contrast_text_color(QColor(bg_hex))


def default_amino_acid_backgrounds() -> dict[str, str]:
    from repeatviewer.backends.msa_plot import DEFAULT_AMINO_ACID_COLORS

    return dict(DEFAULT_AMINO_ACID_COLORS)


def load_custom_amino_acid_backgrounds() -> dict[str, str]:
    import json

    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    raw = s.value("amino_acid_colors_json", "")
    s.endGroup()
    text = str(raw or "").strip()
    if not text:
        return {}
    try:
        decoded = json.loads(text)
    except (json.JSONDecodeError, TypeError):
        return {}
    if not isinstance(decoded, dict):
        return {}
    colors: dict[str, str] = {}
    for key, value in decoded.items():
        color = str(value).strip()
        if key and color.startswith("#") and len(color) >= 4:
            colors[str(key).upper()] = color
    return colors


def save_custom_amino_acid_backgrounds(colors: dict[str, str]) -> None:
    import json

    cleaned: dict[str, str] = {}
    for key, value in colors.items():
        color = str(value).strip()
        if key and color.startswith("#") and len(color) >= 4:
            cleaned[str(key).upper()] = color
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    if cleaned:
        s.setValue("amino_acid_colors_json", json.dumps(cleaned, sort_keys=True))
    else:
        s.remove("amino_acid_colors_json")
    s.endGroup()
    invalidate_residue_color_maps()


def clear_custom_amino_acid_backgrounds() -> None:
    save_custom_amino_acid_backgrounds({})


def amino_acid_color_map() -> dict[str, tuple[str, str]]:
    """Return bg/fg hex colors for amino-acid residues (Annotation protein view)."""
    global _amino_acid_color_map_cache
    if _amino_acid_color_map_cache is not None:
        return _amino_acid_color_map_cache
    merged = default_amino_acid_backgrounds()
    merged.update(load_custom_amino_acid_backgrounds())
    _amino_acid_color_map_cache = {
        aa: (bg, _contrast_fg_for_bg(bg)) for aa, bg in merged.items()
    }
    return _amino_acid_color_map_cache
