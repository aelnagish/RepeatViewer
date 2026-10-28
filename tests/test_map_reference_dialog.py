"""Mapping dialog: Align-style Minimap2 / HISAT2 / Options tabs."""

from __future__ import annotations

from dataclasses import replace

import pytest
from PySide6.QtWidgets import QPushButton, QTabBar

from repeatviewer.job_options import (
    MAP_ENGINE_HISAT2,
    MAP_ENGINE_MINIMAP2,
    MAP_HISAT2_DNA,
    MAP_HISAT2_RNA,
    MAP_HISAT2_READS_PAIRED,
    MAP_HISAT2_STRAND_RF,
    MAP_MULTI_ALL,
    MAP_NO_COVERAGE_SKIP,
    MAP_RIGHT_CLICK_DIALOG,
    MapOptions,
)


def _dialog(qapp, tmp_path, **kw):
    from repeatviewer.views.map_reference_dialog import MapReferenceDialog

    del qapp
    ref = tmp_path / "CroSat1.fa"
    ref.write_text(">sat\n" + "ACGT" * 200 + "\n")
    q = tmp_path / "variants.fa"
    q.write_text(">v1\nACGTACGTAC\n>v2\nACGTACGTAA\n")
    return MapReferenceDialog(files=[q], **kw), ref, q


@pytest.fixture(autouse=True)
def _silence_hisat2_reminder(monkeypatch):
    """HISAT2 tab open schedules an install reminder; keep tests non-interactive."""
    monkeypatch.setattr(
        "repeatviewer.views.map_reference_dialog.MapReferenceDialog._remind_hisat2_if_needed",
        lambda self: None,
    )


def test_buttons_and_align_style_tabs(qapp, tmp_path):
    dlg, _ref, _q = _dialog(qapp, tmp_path)
    texts = [b.text() for b in dlg.findChildren(QPushButton)]
    assert "OK" not in texts
    assert "Restore Defaults" in texts
    assert any(t.startswith("Advanced Settings") for t in texts)
    order = [t for t in texts if t in ("Cancel", "Save", "Run")]
    assert order == ["Cancel", "Save", "Run"]
    tabs = dlg.findChildren(QTabBar)[0]
    assert [tabs.tabText(i) for i in range(tabs.count())] == [
        "Minimap2",
        "HISAT2",
        "Options",
    ]
    assert tabs.expanding() is False  # BLAST-style compact tabs
    assert tabs.documentMode() is True
    assert tabs.currentIndex() == dlg._TAB_MINIMAP2
    assert dlg.width() == dlg._DIALOG_WIDTH


def test_reopens_last_map_tab(qapp, tmp_path):
    from repeatviewer.job_options import DIALOG_TAB_OPTIONS

    dlg, _ref, _q = _dialog(
        qapp,
        tmp_path,
        initial=MapOptions(
            engine=MAP_ENGINE_HISAT2,
            last_dialog_tab=MAP_ENGINE_HISAT2,
        ),
    )
    assert dlg._tab_bar.currentIndex() == dlg._TAB_HISAT2
    assert dlg.current_dialog_tab() == MAP_ENGINE_HISAT2
    assert dlg.options().engine == MAP_ENGINE_HISAT2

    dlg2, _r2, _q2 = _dialog(
        qapp,
        tmp_path,
        initial=MapOptions(last_dialog_tab=DIALOG_TAB_OPTIONS),
    )
    assert dlg2._tab_bar.currentIndex() == dlg2._TAB_MINIMAP2
    assert dlg2.current_dialog_tab() == MAP_ENGINE_MINIMAP2


def test_default_reference_uses_project_genome(qapp, tmp_path):
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\n" + "ACGT" * 100 + "\n")
    dlg, _ref, _q = _dialog(qapp, tmp_path, project_genome=genome)
    assert dlg.options().reference_fasta == str(genome.resolve())


def test_saved_reference_overrides_project_genome(qapp, tmp_path):
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nACGT\n")
    saved = tmp_path / "sat.fa"
    saved.write_text(">sat\nACGTACGT\n")
    dlg, _ref, _q = _dialog(
        qapp,
        tmp_path,
        initial=MapOptions(reference_fasta=str(saved)),
        project_genome=genome,
    )
    assert dlg.options().reference_fasta == str(saved)


def test_engine_follows_tab(qapp, tmp_path):
    dlg, ref, _q = _dialog(qapp, tmp_path)
    dlg.apply_path_from_panel(ref)
    dlg._tab_bar.setCurrentIndex(dlg._TAB_MINIMAP2)
    assert dlg.options().engine == MAP_ENGINE_MINIMAP2
    dlg._tab_bar.setCurrentIndex(dlg._TAB_HISAT2)
    assert dlg.options().engine == MAP_ENGINE_HISAT2


def test_options_round_trip(qapp, tmp_path):
    opts = MapOptions(
        engine=MAP_ENGINE_MINIMAP2,
        multiple_hits=MAP_MULTI_ALL,
        minimap2_preset="map-hifi",
        minimap2_output="sam",
        minimap2_include_cigar=False,
        minimap2_include_cs=False,
        minimap2_kmer=19,
        minimap2_window=19,
        hisat2_library=MAP_HISAT2_DNA,
        hisat2_reads=MAP_HISAT2_READS_PAIRED,
        hisat2_strandness=MAP_HISAT2_STRAND_RF,
        hisat2_dta=True,
        hisat2_min_intron=30,
        hisat2_max_intron=100_000,
        hisat2_no_softclip=True,
        min_identity=0.85,
        min_aligned_bp=40,
        min_mapq=5,
        threads=3,
        consensus_threshold=0.75,
        no_coverage=MAP_NO_COVERAGE_SKIP,
        min_variant_freq=0.1,
        save_unused=False,
        lift_annotations=False,
        max_display_reads=500,
        right_click_action=MAP_RIGHT_CLICK_DIALOG,
        show_process_window=True,
        use_system_path=True,
    )
    dlg, ref, _q = _dialog(qapp, tmp_path, initial=replace(opts))
    dlg.apply_path_from_panel(ref)
    dlg._tab_bar.setCurrentIndex(dlg._TAB_MINIMAP2)
    got = dlg.options()
    assert got.reference_fasta == str(ref)
    assert got.engine == MAP_ENGINE_MINIMAP2
    assert got.minimap2_preset == "map-hifi"
    assert got.minimap2_output == "sam"
    assert got.minimap2_kmer == 19
    assert got.minimap2_window == 19
    assert got.multiple_hits == MAP_MULTI_ALL
    assert got.min_identity == 0.85
    # CIGAR / CS live under Advanced (defaults still apply).
    assert dlg._mm_cigar.isChecked() is False
    assert dlg._mm_cs.isChecked() is False
    assert dlg._mm_cigar_field.parent() is dlg._advanced_panel
    # HISAT2 tab hides minimap2 Output / Advanced rows; shows library settings.
    dlg._tab_bar.setCurrentIndex(dlg._TAB_HISAT2)
    assert dlg._hisat2_library_field.isVisibleTo(dlg)
    assert dlg._hisat2_strand_field.isVisibleTo(dlg)
    assert not dlg._preset_field.isVisibleTo(dlg)
    assert not dlg._mm_format_field.isVisibleTo(dlg)
    hs = dlg.options()
    assert hs.engine == MAP_ENGINE_HISAT2
    assert hs.hisat2_library == MAP_HISAT2_DNA
    assert hs.hisat2_reads == MAP_HISAT2_READS_PAIRED
    assert hs.hisat2_strandness == MAP_HISAT2_STRAND_RF
    assert hs.hisat2_dta is True
    assert hs.hisat2_min_intron == 30
    assert hs.hisat2_max_intron == 100_000
    assert hs.hisat2_no_softclip is True
    # DNA mode disables RNA-only controls.
    assert dlg._hisat2_strand.isEnabled() is False
    dlg._set_combo(dlg._hisat2_library, MAP_HISAT2_RNA)
    dlg._sync_hisat2_library()
    assert dlg._hisat2_strand.isEnabled() is True
    # Changing preset fills concrete -k/-w (including RNA/cDNA splice).
    dlg._tab_bar.setCurrentIndex(dlg._TAB_MINIMAP2)
    dlg._set_combo(dlg._preset, "sr")
    dlg._on_preset_changed()
    assert (dlg._mm_kmer.value(), dlg._mm_window.value()) == (21, 11)
    dlg._set_combo(dlg._preset, "splice")
    dlg._on_preset_changed()
    assert (dlg._mm_kmer.value(), dlg._mm_window.value()) == (15, 5)
    assert "splice" in {dlg._preset.itemData(i) for i in range(dlg._preset.count())}
    assert "splice:hq" in {dlg._preset.itemData(i) for i in range(dlg._preset.count())}
    assert "splice:sr" in {dlg._preset.itemData(i) for i in range(dlg._preset.count())}


def test_reference_browse_aligns_with_define_inputs_browse(qapp, tmp_path):
    """Reference Browse must share the X of Define-inputs Browse (Minimap2 / HISAT2)."""
    dlg, _ref, _q = _dialog(qapp, tmp_path)
    dlg.show()
    qapp.processEvents()
    dlg._input_define.setChecked(True)
    dlg._sync_inputs()
    dlg._sync_field_column()
    qapp.processEvents()
    define_browse = dlg._active_define_inputs()._rows[0]._row._browse_btn
    ref_browse = dlg._reference._browse_btn
    define_x = define_browse.mapTo(dlg, define_browse.rect().topLeft()).x()
    ref_x = ref_browse.mapTo(dlg, ref_browse.rect().topLeft()).x()
    assert abs(define_x - ref_x) <= 2
    # HISAT2 has its own Define-inputs list; Browse still aligns with Reference.
    dlg._tab_bar.setCurrentIndex(dlg._TAB_HISAT2)
    qapp.processEvents()
    hs_browse = dlg._active_define_inputs()._rows[0]._row._browse_btn
    define_x = hs_browse.mapTo(dlg, hs_browse.rect().topLeft()).x()
    ref_x = ref_browse.mapTo(dlg, ref_browse.rect().topLeft()).x()
    assert abs(define_x - ref_x) <= 2


def test_minimap2_and_hisat2_define_inputs_are_independent(qapp, tmp_path):
    """Adding a Define-inputs row on Minimap2 must not add one on HISAT2."""
    dlg, _ref, _q = _dialog(qapp, tmp_path)
    dlg._input_define.setChecked(True)
    dlg._sync_inputs()
    assert len(dlg._mm_define_inputs._rows) == 1
    assert len(dlg._hs_define_inputs._rows) == 1
    dlg._mm_define_inputs.add_row(focus=True)
    assert len(dlg._mm_define_inputs._rows) == 2
    assert len(dlg._hs_define_inputs._rows) == 1
    dlg._tab_bar.setCurrentIndex(dlg._TAB_HISAT2)
    dlg._sync_inputs()
    assert dlg._active_define_inputs() is dlg._hs_define_inputs
    assert len(dlg._active_define_inputs()._rows) == 1
    dlg._hs_define_inputs.add_row(focus=True)
    assert len(dlg._hs_define_inputs._rows) == 2
    assert len(dlg._mm_define_inputs._rows) == 2
    assert dlg._mm_define_inputs._add_btn.text() == "Add another input"



@pytest.mark.parametrize("tab", ["minimap2", "hisat2"])
def test_define_inputs_panel_click_follows_clicked_field(qapp, tmp_path, tab):
    """Define inputs on: a Files click fills the field clicked last, Reference included."""
    from PySide6.QtCore import QEvent, Qt
    from PySide6.QtGui import QFocusEvent
    from PySide6.QtWidgets import QApplication

    def _click_into(line_edit) -> None:
        QApplication.sendEvent(
            line_edit, QFocusEvent(QEvent.Type.FocusIn, Qt.FocusReason.MouseFocusReason)
        )

    dlg, ref, q = _dialog(qapp, tmp_path)
    dlg._reference.clear()
    # Focus landing in Reference on open (not a user click) must not steal clicks.
    _click_into(dlg._reference.line_edit())
    dlg._clear_startup_path_focus()
    if tab == "hisat2":
        dlg._tab_bar.setCurrentIndex(dlg._TAB_HISAT2)
    dlg._input_define.setChecked(True)
    dlg._sync_inputs()
    define = dlg._active_define_inputs()

    assert dlg.apply_path_from_panel(q)
    assert define.paths() == [q]

    _click_into(dlg._reference.line_edit())
    assert dlg.apply_path_from_panel(ref)
    assert dlg._reference.path() == ref
    assert define.paths() == [q]

    other = tmp_path / "reads2.fa"
    other.write_text(">r\nACGT\n")
    _click_into(define._rows[0]._row.line_edit())
    assert dlg.apply_path_from_panel(other)
    assert define.paths() == [other]
    assert dlg._reference.path() == ref


def test_settings_use_align_field_width(qapp, tmp_path):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QLabel

    from repeatviewer.views.tool_dialogs import (
        TOOL_PARAMS_FIELD_WIDTH,
        TOOL_PARAMS_LABEL_COLUMN_REF,
    )

    dlg, _ref, _q = _dialog(qapp, tmp_path)
    dlg._sync_field_column()
    for w in (
        dlg._preset,
        dlg._multi_mm,
        dlg._mm_format,
        dlg._threads,
        dlg._mm_kmer,
        dlg._mm_window,
        dlg._hisat2_library,
        dlg._hisat2_strand,
    ):
        assert w.width() == TOOL_PARAMS_FIELD_WIDTH
    # Inputs radios indented to the shared Settings field column.
    assert dlg._inputs_radios_label.width() == dlg.fontMetrics().horizontalAdvance(
        TOOL_PARAMS_LABEL_COLUMN_REF
    )
    # Short labels sit on the right of that column (next to the combo).
    form = dlg._basic_form
    map_mode_lab = None
    for row in range(form.rowCount()):
        item = form.itemAt(row, form.ItemRole.LabelRole)
        w = item.widget() if item is not None else None
        if isinstance(w, QLabel) and w.text() == "Map mode:":
            map_mode_lab = w
            break
    assert map_mode_lab is not None
    assert map_mode_lab.width() == dlg._inputs_radios_label.width()
    assert int(map_mode_lab.alignment()) & int(Qt.AlignmentFlag.AlignRight)
    assert dlg._preset.itemText(0).startswith("map-ont")
    assert "Primary only" in {
        dlg._multi_mm.itemText(i) for i in range(dlg._multi_mm.count())
    }
    assert "Up to 50 places (-N 50)" in {
        dlg._multi_mm.itemText(i) for i in range(dlg._multi_mm.count())
    }


def test_restore_defaults(qapp, tmp_path):
    dlg, _ref, _q = _dialog(
        qapp, tmp_path, initial=MapOptions(min_identity=0.5, threads=7)
    )
    dlg._restore_defaults()
    got = dlg.options()
    d = MapOptions()
    assert (got.min_identity, got.multiple_hits) == (d.min_identity, d.multiple_hits)
    assert got.minimap2_preset == "map-ont"
    assert got.hisat2_library == MAP_HISAT2_RNA
    assert got.hisat2_dta is False
