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
    return RepeatQuantifyDialog(**kwargs)


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
    d = _dialog()
    d._library_input._rows[0].set_path(long_path.with_name("lib.fa"))
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


# --- Output folder ------------------------------------------------------------


def test_output_folder_defaults_to_auto_and_is_not_saved(settings):
    d = _dialog()
    assert d.output_dir() is None
    assert "<genome folder>/<genome>_repeat_quant_" in d._output_preview.full_text()
    d.set_output_dir(Path("/tmp/results"))
    save_repeat_quantify_options(d.options())
    assert not any("output" in k for k in settings["repeat_quantify"])
    assert _dialog().output_dir() is None


def test_output_preview_follows_genome_and_chosen_folder(tmp_path):
    genome = tmp_path / "Csat_v1.fa"
    genome.write_text(">s\nACGT\n")
    d = _dialog()
    d._genome_input.add_row(genome)
    d._update_output_preview()
    text = d._output_preview.full_text()
    assert str(tmp_path.resolve() / "Csat_v1_repeat_quant_") in text
    results = tmp_path / "results"
    d.set_output_dir(results)
    assert str(results.resolve() / "Csat_v1_repeat_quant_") in d._output_preview.full_text()
    d.set_output_dir(None)
    assert d.output_dir() is None


def test_folder_click_fills_output_only_when_it_had_focus(tmp_path):
    fasta = tmp_path / "genome.fa"
    fasta.write_text(">s\nACGT\n")
    d = _dialog()
    assert not d.apply_folder_from_panel(tmp_path)  # output field never focused
    d._output_focused = True
    assert not d.apply_path_from_panel(fasta)  # a file is not a folder
    assert d.apply_path_from_panel(tmp_path)
    assert d.output_dir() == tmp_path.resolve()


def _run_with_output(tmp_path, monkeypatch, output: str):
    lib = tmp_path / "lib.fa"
    genome = tmp_path / "genome.fa"
    for f in (lib, genome):
        f.write_text(">s\nACGT\n")
    warned = []
    monkeypatch.setattr(
        "repeatviewer.views.tool_dialogs.QMessageBox.warning",
        lambda *a, **k: warned.append(a),
    )
    d = _dialog()
    d._library_input._rows[0].set_path(lib)
    d._genome_input._rows[0].set_path(genome)
    d._output_path.setText(output)
    d._close_action = "run"
    d._on_footer_accept()
    return d, warned


def test_run_rejects_relative_or_file_output_path(tmp_path, monkeypatch):
    d, warned = _run_with_output(tmp_path, monkeypatch, "results")
    assert warned and d.result() != d.DialogCode.Accepted
    d, warned = _run_with_output(tmp_path, monkeypatch, str(tmp_path / "lib.fa"))
    assert warned and d.result() != d.DialogCode.Accepted


def test_run_accepts_a_new_output_folder(tmp_path, monkeypatch):
    d, warned = _run_with_output(tmp_path, monkeypatch, str(tmp_path / "new"))
    assert not warned
    assert d.result() == d.DialogCode.Accepted
    assert d.should_run()


def test_output_preview_is_one_line_and_elided_in_the_middle(tmp_path):
    deep = tmp_path.joinpath(*["very_long_folder_name"] * 12)
    deep.mkdir(parents=True)
    genome = deep / "scaffold_5.fa"
    genome.write_text(">s\nACGT\n")
    d = _dialog()
    d._genome_input.add_row(genome)
    d._update_output_preview()
    label = d._output_preview
    label.resize(420, label.sizeHint().height())
    label._refit()  # hidden widgets get no resize event
    shown = label.text()
    assert "\n" not in shown and not label.wordWrap()
    assert "…" in shown
    assert shown.startswith("→ /")
    run_folder = Path(label.full_text()).name
    assert shown.endswith("/" + run_folder)  # run folder name stays whole
    assert label.toolTip() == label.full_text()
    assert str(deep.resolve()) in label.full_text()


def test_output_preview_keeps_the_more_genomes_note(tmp_path):
    a, b = tmp_path / "a.fa", tmp_path / "b.fa"
    for f in (a, b):
        f.write_text(">s\nACGT\n")
    d = _dialog()
    d._genome_input.add_row(a)
    d._genome_input.add_row(b)
    d._update_output_preview()
    d._output_preview.resize(200, 20)
    d._output_preview._refit()
    assert d._output_preview.text().endswith("(+1 more)")


def test_library_and_genome_start_empty():
    d = _dialog(start_dir="/some/folder")
    assert d.library() is None
    assert d.targets() == []


def test_save_index_box_is_on_by_default_and_synced_across_tabs():
    d = _dialog()
    assert d.save_blast_index() is True
    d._blast_c["save"].setChecked(False)
    assert d._blast_a["save"].isChecked() is False
    assert d.options().blast_save_index is False
    d._restore_defaults()
    assert d.save_blast_index() is True


def test_save_index_choice_round_trips_through_settings(settings):
    save_repeat_quantify_options(RepeatQuantifyOptions(blast_save_index=False))
    assert load_repeat_quantify_options().blast_save_index is False
    settings["repeat_quantify"]["blast_save_index"] = "true"  # QSettings string
    assert load_repeat_quantify_options().blast_save_index is True


def test_status_says_which_genomes_reuse_a_saved_database(tmp_path):
    from repeatviewer.project import BlastDatabase

    saved = tmp_path / "saved.fa"
    fresh = tmp_path / "fresh.fa"
    for f in (saved, fresh):
        f.write_text(">s\nACGT\n")
    prefix = tmp_path / "dbs" / "saved" / "saved"
    prefix.parent.mkdir(parents=True)
    Path(f"{prefix}.nin").write_text("")
    db = BlastDatabase(name="Saffron v3", source_fasta=saved, db_prefix_path=prefix)
    d = _dialog(blast_databases=[db])
    d._genome_input._rows[0].set_path(saved)
    d._genome_input.add_row()
    d._genome_input._rows[1].set_path(fresh)
    lines = d.blast_database_status()
    assert lines[0].startswith("saved.fa: uses saved BLAST database")
    assert "Saffron v3" in lines[0]
    assert "indexed once and added to Add/Manage databases" in lines[1]
    d._blast_a["save"].setChecked(False)
    assert "temporary index" in d.blast_database_status()[1]


def test_status_is_hidden_for_repeatmasker_only(tmp_path):
    genome = tmp_path / "g.fa"
    genome.write_text(">s\nACGT\n")
    d = _dialog()
    d._tabs.setCurrentIndex(0)
    d._genome_input._rows[0].set_path(genome)
    assert d._db_status.text()
    d._tabs.setCurrentIndex(1)  # RepeatMasker
    assert d._db_status.text() == ""
