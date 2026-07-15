from teviewer.job_options import (
    FlexidotOptions,
    load_flexidot_options,
    normalize_flexidot_output_format,
    save_flexidot_options,
)


def test_normalize_flexidot_output_format() -> None:
    assert normalize_flexidot_output_format("png") == "png"
    assert normalize_flexidot_output_format("PNG") == "png"
    assert normalize_flexidot_output_format("html") == "html"
    assert normalize_flexidot_output_format("unknown") == "html"


def test_save_and_load_flexidot_output_format(monkeypatch) -> None:
    saved: dict = {}

    def _save(group, values):
        saved[group] = values

    def _load(group):
        return saved.get(group, {})

    monkeypatch.setattr("teviewer.job_options._save_settings_group", _save)
    monkeypatch.setattr("teviewer.job_options._load_settings_group", _load)

    save_flexidot_options(
        FlexidotOptions(wordsize=5, lcs_shading=False, use_gff=True, output_format="png")
    )
    loaded = load_flexidot_options()
    assert loaded.output_format == "png"
    assert loaded.wordsize == 5
    assert loaded.use_gff is True
