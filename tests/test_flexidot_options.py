from repeatviewer.job_options import (
    FLEXIDOT_RIGHT_CLICK_BUILTIN,
    FLEXIDOT_RIGHT_CLICK_DIALOG,
    FLEXIDOT_RIGHT_CLICK_FLEXIDOT,
    FlexidotOptions,
    load_flexidot_options,
    normalize_flexidot_output_format,
    normalize_flexidot_right_click_action,
    normalize_flexidot_sequence_type,
    save_flexidot_options,
)


def test_normalize_flexidot_output_format() -> None:
    assert normalize_flexidot_output_format("png") == "png"
    assert normalize_flexidot_output_format("PNG") == "png"
    assert normalize_flexidot_output_format("html") == "html"
    # Unknown values fall back to the project default (png).
    assert normalize_flexidot_output_format("unknown") == "png"


def test_normalize_flexidot_sequence_type() -> None:
    assert normalize_flexidot_sequence_type("nuc") == "nuc"
    assert normalize_flexidot_sequence_type("DNA") == "nuc"
    assert normalize_flexidot_sequence_type("aa") == "aa"
    assert normalize_flexidot_sequence_type("protein") == "aa"
    assert normalize_flexidot_sequence_type("") == "nuc"


def test_normalize_flexidot_right_click_action() -> None:
    assert (
        normalize_flexidot_right_click_action("builtin")
        == FLEXIDOT_RIGHT_CLICK_BUILTIN
    )
    assert (
        normalize_flexidot_right_click_action("flexidot")
        == FLEXIDOT_RIGHT_CLICK_FLEXIDOT
    )
    assert (
        normalize_flexidot_right_click_action("dialog")
        == FLEXIDOT_RIGHT_CLICK_DIALOG
    )
    # Legacy “run with saved parameters” → Built-in.
    assert (
        normalize_flexidot_right_click_action("run") == FLEXIDOT_RIGHT_CLICK_BUILTIN
    )
    assert (
        normalize_flexidot_right_click_action("") == FLEXIDOT_RIGHT_CLICK_BUILTIN
    )


def test_save_and_load_flexidot_right_click(monkeypatch) -> None:
    saved: dict = {}

    def _save(group, values):
        saved[group] = values

    def _load(group):
        return saved.get(group, {})

    monkeypatch.setattr("repeatviewer.job_options._save_settings_group", _save)
    monkeypatch.setattr("repeatviewer.job_options._load_settings_group", _load)

    save_flexidot_options(
        FlexidotOptions(right_click_action=FLEXIDOT_RIGHT_CLICK_FLEXIDOT)
    )
    assert load_flexidot_options().right_click_action == FLEXIDOT_RIGHT_CLICK_FLEXIDOT

    saved.clear()
    save_flexidot_options(FlexidotOptions(right_click_action="run"))
    assert load_flexidot_options().right_click_action == FLEXIDOT_RIGHT_CLICK_BUILTIN


def test_save_and_load_flexidot_output_format(monkeypatch) -> None:
    saved: dict = {}

    def _save(group, values):
        saved[group] = values

    def _load(group):
        return saved.get(group, {})

    monkeypatch.setattr("repeatviewer.job_options._save_settings_group", _save)
    monkeypatch.setattr("repeatviewer.job_options._load_settings_group", _load)

    save_flexidot_options(
        FlexidotOptions(
            wordsize=5,
            sequence_type="aa",
            lcs_shading=False,
            use_gff=True,
            output_format="png",
            substitution_count=2,
            wobble_conversion=True,
            max_n_percentage=25.0,
            include_reverse=False,
            only_vs_first=True,
            line_width=2.5,
        )
    )
    loaded = load_flexidot_options()
    assert loaded.output_format == "png"
    assert loaded.wordsize == 5
    assert loaded.sequence_type == "aa"
    assert loaded.use_gff is True
    assert loaded.substitution_count == 2
    assert loaded.wobble_conversion is True
    assert loaded.max_n_percentage == 25.0
    assert loaded.include_reverse is False
    assert loaded.only_vs_first is True
    assert loaded.line_width == 2.5
