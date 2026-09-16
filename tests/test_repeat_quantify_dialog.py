"""Repeat quantify dialog: BLAST/Align-style layout, saved options, path fields.

Only behaviour is checked here; the visual layout was compared by rendering
the dialog offscreen next to the Align dialog.
"""

from pathlib import Path

import pytest
from PySide6.QtWidgets import QApplication

from repeatviewer import job_options
from repeatviewer.job_options import (
    RepeatQuantifyOptions,
    load_repeat_quantify_options,
    save_repeat_quantify_options,
)
from repeatviewer.views.tool_dialogs import RepeatQuantifyDialog

app = QApplication.instance() or QApplication([])


@pytest.fixture
def settings(monkeypatch):
    """In-memory stand-in for QSettings, so tests never touch real prefs."""
    store: dict[str, dict] = {}
    monkeypatch.setattr(
        job_options, "_load_settings_group", lambda g: dict(store.get(g, {}))
    )
    monkeypatch.setattr(
        job_options,
        "_save_settings_group",
        lambda g, values: store.setdefault(g, {}).update(values),
    )
    return store


def _dialog(**kwargs) -> RepeatQuantifyDialog:
    return RepeatQuantifyDialog(candidates=[], **kwargs)


def test_initial_options_fill_every_tab():
    opts = RepeatQuantifyOptions(
        engine="repeatmasker", blast_word_size=15, rm_div=25,
        rm_sensitivity="quick",
    )
    d = _dialog(initial=opts)
    assert d.engine() == "repeatmasker"
    assert d._blast_c["word"].value() == 15  # "Both" copy follows
    assert d._rm_d["div"].value() == 25
    assert d._rm_d["sens"].currentData() == "quick"
    assert d.options() == opts


def test_both_tab_and_single_tabs_stay_in_sync():
    d = _dialog()
    d._blast_c["mincov"].setValue(50)
    d._rm_b["sens"].setCurrentIndex(2)
    assert d._blast_a["mincov"].value() == 50
    assert d._rm_d["sens"].currentData() == "quick"


def test_restore_defaults_resets_parameters_not_files(tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">s\nACGT\n")
    d = _dialog(initial=RepeatQuantifyOptions(engine="blastn", rm_div=10))
    d._genome_input.add_row(fasta)
    d._restore_defaults()
    assert d.options() == RepeatQuantifyOptions()
    assert d.targets() == [fasta.resolve()]


def test_save_closes_without_files_and_does_not_run():
    d = _dialog()
    d._close_action = "save"
    d._on_footer_accept()
    assert d.result() == d.DialogCode.Accepted
    assert not d.should_run()


def test_run_without_files_keeps_dialog_open(monkeypatch):
    shown = []
    monkeypatch.setattr(
        "repeatviewer.views.tool_dialogs.QMessageBox.information",
        lambda *a, **k: shown.append(a),
    )
    d = _dialog()
    d._close_action = "run"
    d._on_footer_accept()
    assert shown
    assert d.result() != d.DialogCode.Accepted


def test_path_field_is_framed_and_shows_full_path_on_hover():
    long_path = Path("/Users/someone/Downloads/run2/RepeatQuantify/scaffold_1.fasta")
    d = RepeatQuantifyDialog(candidates=[long_path, long_path.with_name("lib.fa")])
    field = d._library_input._rows[0]._path
    assert field.objectName() == "PathField"
    assert field.toolTip() == str(long_path.with_name("lib.fa").resolve())
    field.setText("/elsewhere/other.fa")
    assert field.toolTip() == "/elsewhere/other.fa"


def test_options_round_trip_through_settings(settings):
    opts = RepeatQuantifyOptions(
        engine="blastn", blast_word_size=9, blast_evalue=0.001,
        blast_min_identity=80.0, blast_min_coverage=50.0, rm_div=30,
        rm_sensitivity="default",
    )
    save_repeat_quantify_options(opts)
    assert load_repeat_quantify_options() == opts


def test_bad_stored_values_fall_back_to_defaults(settings):
    settings["repeat_quantify"] = {
        "engine": "hmmer", "blast_word_size": "abc", "rm_div": "999",
        "rm_sensitivity": "turbo",
    }
    loaded = load_repeat_quantify_options()
    assert loaded.engine == "both"
    assert loaded.blast_word_size == 11
    assert loaded.rm_div == 100  # clamped to the spin box range
    assert loaded.rm_sensitivity == "sensitive"


def test_params_match_quantify_library_keywords():
    import inspect

    from repeatviewer.backends.repeat_quantify import quantify_library

    accepted = set(inspect.signature(quantify_library).parameters)
    assert set(_dialog().params()) <= accepted
