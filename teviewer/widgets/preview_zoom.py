"""Persisted zoom / fit settings for TEAid HTML and Flexidot image previews."""

from __future__ import annotations

from enum import Enum

from PySide6.QtCore import QSettings

PREVIEW_SETTINGS_GROUP = "preview"
DEFAULT_PREVIEW_ZOOM_PERCENT = 100
MIN_PREVIEW_ZOOM_PERCENT = 25
MAX_PREVIEW_ZOOM_PERCENT = 200
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


DEFAULT_FASTA_PREVIEW_TAB = DefaultFastaPreviewTab.TEXT


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


# Wong / Paul Tol–style nucleotide colors (color-blind friendly).
_COLORBLIND_BASE_COLORS: dict[str, tuple[str, str]] = {
    "A": ("#E69F00", "#000000"),
    "C": ("#0072B2", "#FFFFFF"),
    "G": ("#000000", "#FFFFFF"),
    "T": ("#009E73", "#FFFFFF"),
    "U": ("#009E73", "#FFFFFF"),
    "N": ("#FFFFFF", "#000000"),
}

_STANDARD_BASE_COLORS: dict[str, tuple[str, str]] = {
    "A": ("#00FF00", "#000000"),
    "C": ("#0000FF", "#FFFFFF"),
    "G": ("#000000", "#FFFFFF"),
    "T": ("#FF0000", "#000000"),
    "U": ("#FF0000", "#000000"),
    "N": ("#FFFFFF", "#000000"),
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


def nucleotide_base_color_map() -> dict[str, tuple[str, str]]:
    """Return bg/fg hex colors for A/C/G/T/U/N."""
    if load_color_blind_nucleotide_scheme():
        return dict(_COLORBLIND_BASE_COLORS)
    return dict(_STANDARD_BASE_COLORS)
