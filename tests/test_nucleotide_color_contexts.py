"""Per-context nucleotide colour schemes (single / multi / MSA)."""

from __future__ import annotations

from PySide6.QtCore import QSettings

from repeatviewer.backends.msa_plot import (
    DEFAULT_NUCLEOTIDE_COLORS,
    NUCLEOTIDE_COLORS_ALL_BLACK,
    NUCLEOTIDE_COLORS_COLORBLIND,
)
from repeatviewer.widgets.preview_zoom import (
    NUCLEOTIDE_COLOR_CONTEXTS,
    PREVIEW_SETTINGS_GROUP,
    _NUCLEOTIDE_CONTEXTS_MIGRATED_KEY,
    invalidate_residue_color_maps,
    load_nucleotide_scheme_name,
    match_nucleotide_scheme_name,
    nucleotide_base_color_map,
    nucleotide_palette_colors,
    save_custom_nucleotide_backgrounds,
    save_nucleotide_scheme,
)


def _clear_nucleotide_settings() -> None:
    s = QSettings()
    s.beginGroup(PREVIEW_SETTINGS_GROUP)
    for key in list(s.allKeys()):
        if "nucleotide" in key:
            s.remove(key)
    s.endGroup()
    invalidate_residue_color_maps()


def test_nucleotide_contexts_are_independent() -> None:
    _clear_nucleotide_settings()
    try:
        save_nucleotide_scheme("msa", "All black")
        save_nucleotide_scheme("multi", "Color-blind")
        save_nucleotide_scheme("single", "Default")

        assert load_nucleotide_scheme_name("msa") == "All black"
        assert load_nucleotide_scheme_name("multi") == "Color-blind"
        assert load_nucleotide_scheme_name("single") == "Default"

        assert nucleotide_palette_colors("msa")["A"] == "#000000"
        assert (
            nucleotide_palette_colors("multi")["A"]
            == NUCLEOTIDE_COLORS_COLORBLIND["A"]
        )
        assert (
            nucleotide_palette_colors("single")["A"]
            == DEFAULT_NUCLEOTIDE_COLORS["A"]
        )
        assert nucleotide_base_color_map("msa")["G"][0] == "#000000"
    finally:
        _clear_nucleotide_settings()


def test_all_black_scheme_roundtrip_matches_combo_name() -> None:
    """Regression: All black must not reopen as Default in the combo."""
    _clear_nucleotide_settings()
    try:
        save_custom_nucleotide_backgrounds(
            dict(NUCLEOTIDE_COLORS_ALL_BLACK), "msa"
        )
        assert load_nucleotide_scheme_name("msa") == "All black"
        assert match_nucleotide_scheme_name(NUCLEOTIDE_COLORS_ALL_BLACK) == (
            "All black"
        )
    finally:
        _clear_nucleotide_settings()


def test_named_preset_custom_colors_survive_scheme_switch() -> None:
    """Editing Default then switching to Color-blind must not wipe Default."""
    _clear_nucleotide_settings()
    try:
        save_nucleotide_scheme("msa", "Default")
        custom_default = dict(DEFAULT_NUCLEOTIDE_COLORS)
        custom_default["A"] = "#112233"
        save_custom_nucleotide_backgrounds(custom_default, "msa")
        assert load_nucleotide_scheme_name("msa") == "Default"
        assert nucleotide_palette_colors("msa")["A"].upper() == "#112233"

        save_nucleotide_scheme("msa", "Color-blind")
        assert load_nucleotide_scheme_name("msa") == "Color-blind"
        assert (
            nucleotide_palette_colors("msa")["A"].upper()
            == NUCLEOTIDE_COLORS_COLORBLIND["A"].upper()
        )

        save_nucleotide_scheme("msa", "Default")
        assert load_nucleotide_scheme_name("msa") == "Default"
        assert nucleotide_palette_colors("msa")["A"].upper() == "#112233"
    finally:
        _clear_nucleotide_settings()


def test_legacy_all_black_migrates_to_named_scheme() -> None:
    import json

    _clear_nucleotide_settings()
    try:
        s = QSettings()
        s.beginGroup(PREVIEW_SETTINGS_GROUP)
        s.setValue(
            "nucleotide_colors_json",
            json.dumps(NUCLEOTIDE_COLORS_ALL_BLACK, sort_keys=True),
        )
        s.setValue("color_blind_nucleotide_scheme", False)
        s.remove(_NUCLEOTIDE_CONTEXTS_MIGRATED_KEY)
        s.endGroup()
        invalidate_residue_color_maps()

        for ctx in NUCLEOTIDE_COLOR_CONTEXTS:
            assert load_nucleotide_scheme_name(ctx) == "All black"
            assert nucleotide_palette_colors(ctx)["T"] == "#000000"
    finally:
        _clear_nucleotide_settings()
