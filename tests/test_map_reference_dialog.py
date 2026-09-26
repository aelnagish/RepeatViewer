"""Map to Reference dialog: rules (Cancel / Save / Run, Options tab, "?" help),
option round trip, Files-panel fill of the reference, the summary line."""

from __future__ import annotations

from dataclasses import replace

from PySide6.QtWidgets import QPushButton, QTabWidget

from repeatviewer.job_options import (
    MAP_ENGINE_MINIMAP2,
    MAP_MULTI_ALL,
    MAP_NO_COVERAGE_SKIP,
    MAP_RIGHT_CLICK_DIALOG,
    MAP_SENS_HIGH,
    MapOptions,
)


def _dialog(qapp, tmp_path, **kw):
    from repeatviewer.views.map_reference_dialog import MapReferenceDialog

    ref = tmp_path / "CroSat1.fa"
    ref.write_text(">sat\n" + "ACGT" * 200 + "\n")
    q = tmp_path / "variants.fa"
    q.write_text(">v1\nACGTACGTAC\n>v2\nACGTACGTAA\n")
    return MapReferenceDialog(files=[q], **kw), ref, q


def test_buttons_follow_the_rules(qapp, tmp_path):
    dlg, _ref, _q = _dialog(qapp, tmp_path)
    texts = [b.text() for b in dlg.findChildren(QPushButton)]
    assert "OK" not in texts
    order = [t for t in texts if t in ("Cancel", "Save", "Run")]
    assert order == ["Cancel", "Save", "Run"]
    tabs = dlg.findChild(QTabWidget)
    assert [tabs.tabText(i) for i in range(tabs.count())] == ["Mapping", "Advanced", "Output",
                                                              "Options"]


def test_options_round_trip(qapp, tmp_path):
    opts = MapOptions(engine=MAP_ENGINE_MINIMAP2, sensitivity=MAP_SENS_HIGH,
                      multiple_hits=MAP_MULTI_ALL, minimap2_preset="map-hifi",
                      min_identity=0.85, min_aligned_bp=40, min_mapq=5, threads=3,
                      consensus_threshold=0.75, no_coverage=MAP_NO_COVERAGE_SKIP,
                      min_variant_freq=0.1, save_unused=False, lift_annotations=False,
                      max_display_reads=500, right_click_action=MAP_RIGHT_CLICK_DIALOG,
                      show_process_window=True, use_system_path=True)
    dlg, ref, _q = _dialog(qapp, tmp_path, initial=replace(opts))
    dlg.apply_path_from_panel(ref)
    assert dlg.options() == replace(opts, reference_fasta=str(ref))


def test_panel_click_fills_reference_only_with_fasta(qapp, tmp_path):
    dlg, ref, _q = _dialog(qapp, tmp_path)
    other = tmp_path / "notes.txt"
    other.write_text("x")
    assert dlg.apply_path_from_panel(other) is False
    assert dlg.apply_path_from_panel(ref) is True
    assert dlg.options().reference_fasta == str(ref)
    assert "2 sequences (variants.fa) will be mapped to CroSat1.fa" in dlg._summary.text()


def test_restore_defaults(qapp, tmp_path):
    dlg, _ref, _q = _dialog(qapp, tmp_path, initial=MapOptions(min_identity=0.5, threads=7))
    dlg._restore_defaults()
    got = dlg.options()
    d = MapOptions()
    assert (got.min_identity, got.engine, got.multiple_hits) == (d.min_identity, d.engine,
                                                                 d.multiple_hits)


def test_summary_without_selection(qapp, tmp_path):
    from repeatviewer.views.map_reference_dialog import input_summary_text

    assert "Select the FASTA" in input_summary_text([], None)
