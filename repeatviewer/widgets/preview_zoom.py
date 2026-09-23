"""Persisted zoom / fit settings for TE-Aid HTML and Flexidot image previews."""

from __future__ import annotations

from enum import Enum

from PySide6.QtCore import QSettings
from PySide6.QtGui import QWheelEvent

PREVIEW_SETTINGS_GROUP = "preview"
DEFAULT_PREVIEW_ZOOM_PERCENT = 100
MIN_PREVIEW_ZOOM_PERCENT = 25
MAX_PREVIEW_ZOOM_PERCENT = 200
PREVIEW_WHEEL_ZOOM_STEP_PERCENT = 5
DEFAULT_PREVIEW_FONT_SIZE = 11
MIN_PREVIEW_FONT_SIZE = 8
MAX_PREVIEW_FONT_SIZE = 24


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


DEFAULT_FASTA_PREVIEW_TAB = DefaultFastaPreviewTab.ANNOTATION


def load_preview_zoom_percent() -> int:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    value = int(s.value("zoom_percent", DEFAULT_PREVIEW_ZOOM_PERCENT))
    s.endGroup()
    return max(MIN_PREVIEW_ZOOM_PERCENT, min(MAX_PREVIEW_ZOOM_PERCENT, value))


def save_preview_zoom_percent(percent: int) -> None:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    s.setValue(
        "zoom_percent",
        max(MIN_PREVIEW_ZOOM_PERCENT, min(MAX_PREVIEW_ZOOM_PERCENT, int(percent))),
    )
    s.endGroup()


def load_preview_fit_mode() -> PreviewFitMode:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    raw = str(s.value("fit_mode", PreviewFitMode.FIT_WINDOW.value))
    s.endGroup()
    try:
        return PreviewFitMode(raw)
    except ValueError:
        return PreviewFitMode.FIT_WINDOW


def wheel_zoom_steps_from_event(event: QWheelEvent) -> int:
    """Return signed zoom steps from a wheel event (positive = zoom in)."""
    delta = event.angleDelta().y()
    if delta == 0:
        delta = event.pixelDelta().y()
    if delta == 0:
        return 0
    steps = int(round(delta / 120.0))
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
    return max(MIN_PREVIEW_FONT_SIZE, min(MAX_PREVIEW_FONT_SIZE, value))


def save_preview_font_size(size: int) -> None:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    s.setValue(
        "text_font_size",
        max(MIN_PREVIEW_FONT_SIZE, min(MAX_PREVIEW_FONT_SIZE, int(size))),
    )
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
    """Load custom nucleotide colours (full alphabet JSON)."""
    import json

    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    raw = s.value("nucleotide_colors_json", "")
    # Drop legacy per-base keys so old AliView saves do not override MSA Default.
    for base in ("A", "C", "G", "T"):
        s.remove(f"nucleotide_{base}_bg")
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


def clear_custom_nucleotide_backgrounds() -> None:
    save_custom_nucleotide_backgrounds({})


def default_nucleotide_backgrounds() -> dict[str, str]:
    """MSA-style Default palette (figure 2 / Element colors Default preset)."""
    from repeatviewer.backends.msa_plot import DEFAULT_NUCLEOTIDE_COLORS

    return dict(DEFAULT_NUCLEOTIDE_COLORS)


def nucleotide_base_color_map() -> dict[str, tuple[str, str]]:
    """Return bg/fg hex colors for nucleotide letters in the Annotation view."""
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
    return {base: (bg, _contrast_fg_for_bg(bg)) for base, bg in merged.items()}


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


def clear_custom_amino_acid_backgrounds() -> None:
    save_custom_amino_acid_backgrounds({})


def amino_acid_color_map() -> dict[str, tuple[str, str]]:
    """Return bg/fg hex colors for amino-acid residues (Annotation protein view)."""
    merged = default_amino_acid_backgrounds()
    merged.update(load_custom_amino_acid_backgrounds())
    return {aa: (bg, _contrast_fg_for_bg(bg)) for aa, bg in merged.items()}
