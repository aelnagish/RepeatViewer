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


# Remembered inline Split layout per FASTA shape (multi-seq MSA vs single-seq).
PreviewSplitKind = str  # "multi" | "single"


def load_preview_split_preference(
    kind: PreviewSplitKind,
) -> tuple[bool, int, int] | None:
    """Return ``(enabled, left_tab, right_tab)`` for ``kind``, or None if unset.

    ``enabled`` means the next open of that FASTA shape should auto-split.
    Tab indices: 0=Text, 1=Overview, 2=Dotplot, 3=Cytogenetics.
    """
    key = str(kind or "").strip().lower()
    if key not in {"multi", "single"}:
        return None
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    s.beginGroup("split_layout")
    if not s.contains(f"{key}/enabled"):
        s.endGroup()
        s.endGroup()
        return None
    enabled = bool(s.value(f"{key}/enabled", False))
    try:
        left = int(s.value(f"{key}/left", 1))
        right = int(s.value(f"{key}/right", 2))
    except (TypeError, ValueError):
        left, right = 1, 2
    s.endGroup()
    s.endGroup()
    left = max(0, min(3, left))
    right = max(0, min(3, right))
    return enabled, left, right


def save_preview_split_preference(
    kind: PreviewSplitKind,
    *,
    enabled: bool,
    left: int,
    right: int,
) -> None:
    """Persist whether / how to auto-split the next file of this FASTA shape."""
    key = str(kind or "").strip().lower()
    if key not in {"multi", "single"}:
        return
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    s.beginGroup("split_layout")
    s.setValue(f"{key}/enabled", bool(enabled))
    s.setValue(f"{key}/left", max(0, min(3, int(left))))
    s.setValue(f"{key}/right", max(0, min(3, int(right))))
    s.endGroup()
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
# G matches Default charcoal so black letters stay readable on fill.
_COLORBLIND_BASE_COLORS: dict[str, tuple[str, str]] = {
    "A": ("#E69F00", "#000000"),
    "C": ("#0072B2", "#FFFFFF"),
    "G": ("#707070", "#000000"),
    "T": ("#009E73", "#FFFFFF"),
    "U": ("#009E73", "#FFFFFF"),
    "N": ("#F0C000", "#000000"),
}

# Same MSA "Default" palette used by Element colors / MSA plot.
# G is charcoal grey (not pure black) so black letters remain readable.
_STANDARD_BASE_COLORS: dict[str, tuple[str, str]] = {
    "A": ("#00C800", "#000000"),
    "C": ("#7040E8", "#000000"),
    "G": ("#707070", "#000000"),
    "T": ("#F06090", "#000000"),
    "U": ("#F06090", "#000000"),
    "N": ("#F0C000", "#000000"),
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
# Per-context caches: single-sequence map, multi-sequence overview, MSA overview.
_nucleotide_color_map_cache: dict[str, dict[str, tuple[str, str]]] = {}
_amino_acid_color_map_cache: dict[str, tuple[str, str]] | None = None

# Independent nucleotide palettes for each sequence display mode.
NUCLEOTIDE_COLOR_CONTEXTS: tuple[str, ...] = ("single", "multi", "msa")
_NUCLEOTIDE_SCHEME_NAMES: frozenset[str] = frozenset(
    {"Default", "Color-blind", "All black", "Custom"}
)
_NUCLEOTIDE_NAMED_PRESETS: frozenset[str] = frozenset(
    {"Default", "Color-blind", "All black"}
)
_NUCLEOTIDE_SCHEME_CORE: tuple[str, ...] = ("A", "C", "G", "T", "U")
_NUCLEOTIDE_CONTEXTS_MIGRATED_KEY = "nucleotide_contexts_migrated_v1"


def normalize_nucleotide_color_context(context: str | None) -> str:
    key = str(context or "single").strip().lower()
    if key in NUCLEOTIDE_COLOR_CONTEXTS:
        return key
    return "single"


def invalidate_residue_color_maps(context: str | None = None) -> None:
    """Drop cached nucleotide / amino-acid palettes (call after color settings change)."""
    global _nucleotide_color_map_cache, _amino_acid_color_map_cache
    if context is None:
        _nucleotide_color_map_cache = {}
    else:
        _nucleotide_color_map_cache.pop(normalize_nucleotide_color_context(context), None)
    _amino_acid_color_map_cache = None


def _scheme_base_colors() -> dict[str, tuple[str, str]]:
    if load_color_blind_nucleotide_scheme():
        return dict(_COLORBLIND_BASE_COLORS)
    return dict(_STANDARD_BASE_COLORS)


def _strip_legacy_aliview_overrides(colors: dict[str, str]) -> dict[str, str]:
    """Drop a full legacy AliView ACGT set so MSA Default is the true default.

    Only strip when A/C/G/T all match the old AliView palette together.
    Per-base stripping wrongly removed G from the All-black preset
    (legacy AliView G was also ``#000000``).
    """
    cleaned = dict(colors)
    cores = ("A", "C", "G", "T")
    if all(
        cleaned.get(base, "").upper() == _LEGACY_ALIVIEW_ACGT[base].upper()
        for base in cores
    ):
        for base in (*cores, "U"):
            cleaned.pop(base, None)
    return cleaned


def _normalize_tu_colors(colors: dict[str, str]) -> dict[str, str]:
    cleaned = dict(colors)
    if "U" in cleaned:
        cleaned["T"] = cleaned["U"]
    elif "T" in cleaned:
        cleaned["U"] = cleaned["T"]
    return cleaned


def _parse_nucleotide_colors_json(raw: object) -> dict[str, str]:
    import json

    custom: dict[str, str] = {}
    text = str(raw or "").strip()
    if not text:
        return custom
    try:
        decoded = json.loads(text)
    except (json.JSONDecodeError, TypeError):
        return custom
    if not isinstance(decoded, dict):
        return custom
    for key, value in decoded.items():
        color = str(value).strip()
        if key and color.startswith("#") and len(color) >= 4:
            custom[str(key).upper()] = color
    return _normalize_tu_colors(_strip_legacy_aliview_overrides(custom))


def _clean_nucleotide_color_dict(colors: dict[str, str]) -> dict[str, str]:
    cleaned: dict[str, str] = {}
    for key, value in colors.items():
        color = str(value).strip()
        if key and color.startswith("#") and len(color) >= 4:
            cleaned[str(key).upper()] = color
    return _normalize_tu_colors(cleaned)


def match_nucleotide_scheme_name(colors: dict[str, str]) -> str | None:
    """Return a preset name when ``colors`` match that preset on A/C/G/T/U."""
    from repeatviewer.backends.msa_plot import NUCLEOTIDE_COLOR_PRESETS

    cleaned = _clean_nucleotide_color_dict(colors)
    if not cleaned:
        return "Default"
    for name, palette in NUCLEOTIDE_COLOR_PRESETS:
        if all(
            cleaned.get(base, "").upper() == str(palette.get(base, "")).upper()
            for base in _NUCLEOTIDE_SCHEME_CORE
        ):
            return name
    return None


def _nucleotide_preset_palette(name: str) -> dict[str, str]:
    from repeatviewer.backends.msa_plot import (
        DEFAULT_NUCLEOTIDE_COLORS,
        NUCLEOTIDE_COLOR_PRESETS,
    )

    for preset_name, palette in NUCLEOTIDE_COLOR_PRESETS:
        if preset_name == name:
            return dict(palette)
    return dict(DEFAULT_NUCLEOTIDE_COLORS)


def _nucleotide_scheme_key(context: str) -> str:
    return f"nucleotide_scheme_{context}"


def _nucleotide_colors_key(context: str) -> str:
    return f"nucleotide_colors_json_{context}"


def _nucleotide_preset_key_token(preset: str) -> str:
    return str(preset or "Default").strip().lower().replace(" ", "_")


def _nucleotide_preset_overrides_key(context: str, preset: str) -> str:
    ctx = normalize_nucleotide_color_context(context)
    return (
        f"nucleotide_preset_overrides_json_{ctx}_"
        f"{_nucleotide_preset_key_token(preset)}"
    )


def load_nucleotide_preset_overrides(
    context: str, preset: str
) -> dict[str, str]:
    """Custom colours stored for one named preset in a display context."""
    name = str(preset or "").strip()
    if name not in _NUCLEOTIDE_NAMED_PRESETS:
        return {}
    _ensure_nucleotide_context_migration()
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    raw = s.value(_nucleotide_preset_overrides_key(context, name), "")
    s.endGroup()
    return _parse_nucleotide_colors_json(raw)


def _nucleotide_palette_matches_stock(
    cleaned: dict[str, str], stock: dict[str, str]
) -> bool:
    """True when ``cleaned`` agrees with ``stock`` on every shared key."""
    if not cleaned:
        return True
    for key, value in cleaned.items():
        stock_val = stock.get(key)
        if stock_val is None:
            return False
        if str(value).upper() != str(stock_val).upper():
            return False
    for key, value in stock.items():
        if key in cleaned and str(cleaned[key]).upper() != str(value).upper():
            return False
    return True


def save_nucleotide_preset_overrides(
    context: str,
    preset: str,
    colors: dict[str, str] | None,
) -> None:
    """Persist (or clear) per-preset custom colours for a display context."""
    import json

    name = str(preset or "").strip()
    if name not in _NUCLEOTIDE_NAMED_PRESETS:
        return
    _ensure_nucleotide_context_migration()
    ctx = normalize_nucleotide_color_context(context)
    cleaned = _clean_nucleotide_color_dict(colors or {})
    stock = _nucleotide_preset_palette(name)
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    key = _nucleotide_preset_overrides_key(ctx, name)
    if cleaned and not _nucleotide_palette_matches_stock(cleaned, stock):
        s.setValue(key, json.dumps(cleaned, sort_keys=True))
    else:
        s.remove(key)
    s.endGroup()
    invalidate_residue_color_maps(ctx)


def _write_nucleotide_context_settings(
    context: str,
    scheme: str,
    colors: dict[str, str] | None,
) -> None:
    import json

    ctx = normalize_nucleotide_color_context(context)
    name = str(scheme or "Default").strip()
    if name not in _NUCLEOTIDE_SCHEME_NAMES:
        name = "Default"
    cleaned = _clean_nucleotide_color_dict(colors or {})
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    s.setValue(_nucleotide_scheme_key(ctx), name)
    for base in ("A", "C", "G", "T"):
        s.remove(f"nucleotide_{base}_bg")
        s.remove(f"nucleotide_{base}_bg_{ctx}")
    if name == "Custom" and cleaned:
        s.setValue(_nucleotide_colors_key(ctx), json.dumps(cleaned, sort_keys=True))
    else:
        s.remove(_nucleotide_colors_key(ctx))
    s.endGroup()
    invalidate_residue_color_maps(ctx)


def _infer_legacy_nucleotide_scheme(
    *,
    color_blind: bool,
    custom: dict[str, str],
) -> tuple[str, dict[str, str]]:
    """Map pre-context settings to (scheme name, custom colors to store)."""
    from repeatviewer.backends.msa_plot import (
        DEFAULT_NUCLEOTIDE_COLORS,
        NUCLEOTIDE_COLORS_COLORBLIND,
    )

    if color_blind:
        base = dict(NUCLEOTIDE_COLORS_COLORBLIND)
    else:
        base = dict(DEFAULT_NUCLEOTIDE_COLORS)
    if not custom:
        return ("Color-blind" if color_blind else "Default"), {}
    merged = dict(base)
    merged.update(custom)
    matched = match_nucleotide_scheme_name(merged)
    if matched is not None:
        return matched, {}
    return "Custom", merged


def _ensure_nucleotide_context_migration() -> None:
    """Copy legacy global nucleotide settings into per-context keys once."""
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    migrated = bool(s.value(_NUCLEOTIDE_CONTEXTS_MIGRATED_KEY, False))
    if migrated:
        s.endGroup()
        return
    raw = s.value("nucleotide_colors_json", "")
    color_blind = bool(s.value("color_blind_nucleotide_scheme", False))
    s.endGroup()

    custom = _parse_nucleotide_colors_json(raw)
    scheme, store = _infer_legacy_nucleotide_scheme(
        color_blind=color_blind, custom=custom
    )
    for ctx in NUCLEOTIDE_COLOR_CONTEXTS:
        _write_nucleotide_context_settings(ctx, scheme, store)

    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    s.setValue(_NUCLEOTIDE_CONTEXTS_MIGRATED_KEY, True)
    s.endGroup()


def load_nucleotide_scheme_name(context: str = "single") -> str:
    """Saved preset name for a display context (Default / Color-blind / All black / Custom)."""
    _ensure_nucleotide_context_migration()
    ctx = normalize_nucleotide_color_context(context)
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    name = str(s.value(_nucleotide_scheme_key(ctx), "Default") or "Default").strip()
    s.endGroup()
    if name in _NUCLEOTIDE_SCHEME_NAMES:
        return name
    return "Default"


def save_nucleotide_scheme(
    context: str,
    name: str,
    colors: dict[str, str] | None = None,
) -> None:
    """Persist the active named preset or Custom palette for one context.

    Switching between Default / Color-blind / All black keeps each preset's
    own colour overrides; only the active scheme name changes.
    """
    _ensure_nucleotide_context_migration()
    label = str(name or "Default").strip()
    if label == "Custom":
        cleaned = _clean_nucleotide_color_dict(colors or {})
        matched = match_nucleotide_scheme_name(cleaned) if cleaned else "Default"
        if matched is not None:
            _write_nucleotide_context_settings(context, matched, None)
        else:
            _write_nucleotide_context_settings(context, "Custom", cleaned)
        return
    if label not in _NUCLEOTIDE_SCHEME_NAMES:
        label = "Default"
    # Named presets: optional ``colors`` update that preset's overrides.
    if label in _NUCLEOTIDE_NAMED_PRESETS and colors is not None:
        save_nucleotide_preset_overrides(context, label, colors)
    _write_nucleotide_context_settings(context, label, None)


def load_color_blind_nucleotide_scheme(context: str = "single") -> bool:
    """Back-compat: True when the context scheme is Color-blind."""
    return load_nucleotide_scheme_name(context) == "Color-blind"


def save_color_blind_nucleotide_scheme(
    enabled: bool, context: str = "single"
) -> None:
    """Back-compat: set Color-blind or Default for one context."""
    save_nucleotide_scheme(
        context, "Color-blind" if enabled else "Default"
    )


def load_custom_nucleotide_backgrounds(context: str = "single") -> dict[str, str]:
    """Load Custom-palette colours for a display context (empty for named presets).

    Read-only: do not write/remove settings here — the Annotation color bar
    calls into this path thousands of times when the palette is uncached.
    """
    _ensure_nucleotide_context_migration()
    ctx = normalize_nucleotide_color_context(context)
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    raw = s.value(_nucleotide_colors_key(ctx), "")
    s.endGroup()
    return _parse_nucleotide_colors_json(raw)


def save_custom_nucleotide_backgrounds(
    colors: dict[str, str], context: str = "single"
) -> None:
    """Save the current palette under the active scheme for one context.

    Named presets keep their own overrides so switching Default ↔ Color-blind
    restores each preset's custom colours. A palette whose A/C/G/T/U match a
    stock preset activates that preset (and stores any non-stock extras).
    """
    _ensure_nucleotide_context_migration()
    cleaned = _clean_nucleotide_color_dict(colors)
    if not cleaned:
        save_nucleotide_scheme(context, "Default")
        return
    matched = match_nucleotide_scheme_name(cleaned)
    scheme = load_nucleotide_scheme_name(context)
    if matched is not None:
        save_nucleotide_preset_overrides(context, matched, cleaned)
        _write_nucleotide_context_settings(context, matched, None)
        return
    if scheme in _NUCLEOTIDE_NAMED_PRESETS:
        save_nucleotide_preset_overrides(context, scheme, cleaned)
        _write_nucleotide_context_settings(context, scheme, None)
        return
    _write_nucleotide_context_settings(context, "Custom", cleaned)


def clear_custom_nucleotide_backgrounds(context: str = "single") -> None:
    """Reset the active context to stock Default (clears Default overrides)."""
    save_nucleotide_preset_overrides(context, "Default", None)
    save_nucleotide_scheme(context, "Default")


def default_nucleotide_backgrounds() -> dict[str, str]:
    """MSA-style Default palette (figure 2 / Element colors Default preset)."""
    from repeatviewer.backends.msa_plot import DEFAULT_NUCLEOTIDE_COLORS

    return dict(DEFAULT_NUCLEOTIDE_COLORS)


def nucleotide_palette_colors(context: str = "single") -> dict[str, str]:
    """Resolved background hex colours for one display context."""
    ctx = normalize_nucleotide_color_context(context)
    scheme = load_nucleotide_scheme_name(ctx)
    if scheme == "Custom":
        merged = default_nucleotide_backgrounds()
        merged.update(load_custom_nucleotide_backgrounds(ctx))
    else:
        merged = _nucleotide_preset_palette(scheme)
        merged.update(load_nucleotide_preset_overrides(ctx, scheme))
    return _normalize_tu_colors(merged)


def nucleotide_base_color_map(context: str = "single") -> dict[str, tuple[str, str]]:
    """Return bg/fg hex colors for nucleotide letters in one display context."""
    ctx = normalize_nucleotide_color_context(context)
    cached = _nucleotide_color_map_cache.get(ctx)
    if cached is not None:
        return cached

    merged = nucleotide_palette_colors(ctx)
    resolved = {
        base: (bg, _contrast_fg_for_bg(bg)) for base, bg in merged.items()
    }
    _nucleotide_color_map_cache[ctx] = resolved
    return resolved


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
