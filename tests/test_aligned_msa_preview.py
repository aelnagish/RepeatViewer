"""Aligned FASTA detection and MSA residue-row overview load."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest
from PySide6.QtWidgets import QApplication

from repeatviewer.backends.sequence_features import (
    aligned_fasta_width,
    is_aligned_fasta,
    load_aligned_fasta_rows,
)
from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel
from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail


def test_is_aligned_fasta_equal_length_gapped(tmp_path: Path) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">s1\nACGT--GT\n>s2\nACGTACGT\n")
    assert is_aligned_fasta(path) is True
    assert aligned_fasta_width(path) == 8


def test_is_aligned_fasta_false_for_unequal_lengths(tmp_path: Path) -> None:
    path = tmp_path / "multi.fa"
    path.write_text(">s1\nACGT\n>s2\nACGTACGT\n")
    assert is_aligned_fasta(path) is False
    assert aligned_fasta_width(path) is None


def test_is_aligned_fasta_false_for_single_record(tmp_path: Path) -> None:
    path = tmp_path / "one.fa"
    path.write_text(">only\nACGTACGT\n")
    assert is_aligned_fasta(path) is False
    assert aligned_fasta_width(path) is None


def test_load_aligned_fasta_rows_uppercases(tmp_path: Path) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">a\nacgt--gt\n>b\nACGTACGT\n")
    ids, seqs = load_aligned_fasta_rows(path)
    assert ids == ["a", "b"]
    assert seqs == ["ACGT--GT", "ACGTACGT"]


def test_overview_load_path_enters_msa_mode(qapp: QApplication, tmp_path: Path) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(
        ">seqA\nACGT--GTAC\n"
        ">seqB\nACGTACGTAC\n"
        ">seqC\nACGTAAGTAC\n"
    )
    panel = MultiSequenceOverviewPanel(path)
    panel.resize(800, 400)
    panel.load_path(path, msa_mode=True)
    assert panel.is_msa_mode() is True
    assert panel._bars.record_ids() == ["seqA", "seqB", "seqC"]
    assert panel._bars.max_length() == 10
    # Force a paint pass; must not raise.
    panel._bars.repaint()
    qapp.processEvents()
    panel.deleteLater()


def test_name_multi_drag_keeps_horizontal_scroll_when_zoomed(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Dragging several name rows must not jump a zoomed MSA back to bp 1."""
    long_cols = "ACGT" * 800
    path = tmp_path / "aln.fa"
    rows = "\n".join(f">s{i}\n{long_cols}" for i in range(8))
    path.write_text(rows + "\n", encoding="utf-8")

    panel = MultiSequenceOverviewPanel(path)
    panel.resize(900, 500)
    panel.show()
    panel.load_path(path, msa_mode=True)
    qapp.processEvents()
    panel._fit_bars()
    panel._bars.zoom_by(10.0, anchor_x=100.0)
    qapp.processEvents()
    panel._bars_scroll.horizontalScrollBar().setValue(500)
    before = panel._bars_scroll.horizontalScrollBar().value()
    assert before > 0

    # Same call path as name-column multi-row drag preview.
    panel._on_name_rows_pixel_preview(40.0, 120.0)
    qapp.processEvents()

    assert panel._bars_scroll.horizontalScrollBar().value() == before
    panel.deleteLater()


def test_overview_switching_msa_resets_scroll_to_origin(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Loading a different MSA must show bp 1 / first row, not the prior scroll."""
    long_cols = "ACGT" * 500
    first = tmp_path / "first.fa"
    first.write_text(f">a\n{long_cols}\n>b\n{long_cols}\n", encoding="utf-8")
    second = tmp_path / "second.fa"
    second.write_text(f">x\n{long_cols}\n>y\n{long_cols}\n", encoding="utf-8")

    panel = MultiSequenceOverviewPanel(first)
    panel.resize(800, 400)
    panel.show()
    panel.load_path(first, msa_mode=True)
    qapp.processEvents()
    panel._fit_bars()
    # Simulate user zoom + scroll into the middle of the first MSA.
    panel._bars.zoom_by(8.0, anchor_x=200.0)
    qapp.processEvents()
    panel._bars_scroll.horizontalScrollBar().setValue(400)
    panel._bars_scroll.verticalScrollBar().setValue(30)
    assert panel._bars_scroll.horizontalScrollBar().value() > 0

    panel.load_path(second, msa_mode=True)
    qapp.processEvents()
    panel._try_initial_fit()
    panel._reset_scrollbars_to_origin()
    qapp.processEvents()

    assert panel._bars_scroll.horizontalScrollBar().value() == 0
    assert panel._bars_scroll.verticalScrollBar().value() == 0
    assert panel._names_scroll.verticalScrollBar().value() == 0
    # Fit zoom for the new file (not the previous zoomed-in ppb).
    fit_ppb, _ = panel._bars.zoom_bounds(panel._bars_viewport_width())
    assert panel._bars.pixels_per_base() == pytest.approx(fit_ppb, rel=0.05)
    panel.deleteLater()


def test_overview_fit_bars_keeps_vertical_scroll(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Fit / viewport refit must not jump back to the first sequence row."""
    long_cols = "ACGT" * 200
    rows = "\n".join(f">s{i}\n{long_cols}" for i in range(40))
    path = tmp_path / "tall.fa"
    path.write_text(rows + "\n", encoding="utf-8")

    panel = MultiSequenceOverviewPanel(path)
    panel.resize(900, 360)
    panel.show()
    panel.load_path(path, msa_mode=True)
    qapp.processEvents()
    panel._fit_bars()
    qapp.processEvents()

    vbar = panel._bars_scroll.verticalScrollBar()
    assert vbar.maximum() > 0
    vbar.setValue(min(180, vbar.maximum()))
    panel._names_scroll.verticalScrollBar().setValue(vbar.value())
    before = int(vbar.value())
    assert before > 0

    # Same path as window / split-pane resize while at Fit zoom.
    panel._refit_after_viewport_change()
    qapp.processEvents()
    assert int(vbar.value()) == before
    assert int(panel._names_scroll.verticalScrollBar().value()) == before

    panel._fit_bars()
    qapp.processEvents()
    assert int(vbar.value()) == before
    panel.deleteLater()


def test_overview_capture_restore_view_state_keeps_row(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Split panes clone the source Overview scroll via capture/restore."""
    long_cols = "ACGT" * 200
    rows = "\n".join(f">s{i}\n{long_cols}" for i in range(40))
    path = tmp_path / "split_src.fa"
    path.write_text(rows + "\n", encoding="utf-8")

    panel = MultiSequenceOverviewPanel(path)
    panel.resize(900, 360)
    panel.show()
    panel.load_path(path, msa_mode=True)
    qapp.processEvents()
    panel._fit_bars()
    qapp.processEvents()

    vbar = panel._bars_scroll.verticalScrollBar()
    assert vbar.maximum() > 0
    vbar.setValue(min(220, vbar.maximum()))
    panel._names_scroll.verticalScrollBar().setValue(vbar.value())
    state = panel.capture_view_state()
    before_anchor = state.get("row_anchor")
    assert before_anchor is not None

    # Simulate a fresh split pane load that would otherwise pin to row 0.
    panel._reset_scrollbars_to_origin()
    qapp.processEvents()
    assert int(vbar.value()) == 0

    panel.set_pending_view_restore(state)
    assert panel.restore_view_state(state) is True
    qapp.processEvents()
    restored = panel.linked_viewport_row_anchor()
    assert restored is not None
    assert int(restored[0]) == int(before_anchor[0])
    panel.deleteLater()


def test_overview_load_path_unaligned_stays_bars(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "multi.fa"
    path.write_text(">a\nACGT\n>b\nACGTACGTACGT\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=False)
    assert panel.is_msa_mode() is False
    assert panel._bars.max_length() == 12
    assert panel._bars.msa_sequences() == ["ACGT", "ACGTACGTACGT"]
    assert panel._bars.zoom_bounds()[1] >= 8.0
    panel._bars.set_pixels_per_base(8.0)
    assert panel._bars._bar_shows_letters() is True
    panel.deleteLater()


def test_unaligned_overview_hides_msa_only_controls(
    qapp: QApplication, tmp_path: Path
) -> None:
    from repeatviewer.widgets.preview_controls import AnnotationZoomBar

    path = tmp_path / "multi.fa"
    path.write_text(">a\nACGT\n>b\nACGTACGT\n")
    panel = MultiSequenceOverviewPanel(path)
    bar = AnnotationZoomBar()
    panel.bind_annotation_zoom_bar(bar)
    panel.load_path(path, msa_mode=False)
    assert bar._clean_btn.isHidden()
    assert bar._consensus_btn.isHidden()
    assert panel._rail.msa_coloring_header.isHidden()
    assert panel._rail.msa_coloring_mode_host.isHidden()
    assert panel._rail.msa_divergence_row.isHidden()
    assert panel._rail._column_distance_label.isHidden()
    assert panel._rail.column_distance.isHidden()
    assert not panel._rail.color_multi_sequence_background.isHidden()
    assert panel._rail.color_multi_sequence_background.isChecked()
    with patch.object(panel, "_save_overview_display_settings"):
        panel._rail.color_multi_sequence_background.setChecked(False)
    assert panel._bars.color_unaligned_background() is False
    panel.deleteLater()
    bar.deleteLater()


def test_aligned_overview_shows_msa_only_controls(
    qapp: QApplication, tmp_path: Path
) -> None:
    from repeatviewer.widgets.preview_controls import AnnotationZoomBar

    path = tmp_path / "aln.fa"
    path.write_text(">a\nACGTACGT\n>b\nACGT--GT\n")
    panel = MultiSequenceOverviewPanel(path)
    bar = AnnotationZoomBar()
    panel.bind_annotation_zoom_bar(bar)
    panel.load_path(path, msa_mode=True)
    assert not bar._clean_btn.isHidden()
    assert not bar._consensus_btn.isHidden()
    assert not panel._rail.msa_coloring_header.isHidden()
    assert not panel._rail.msa_coloring_mode_host.isHidden()
    assert not panel._rail.msa_divergence_row.isHidden()
    assert panel._rail._column_distance_label.isHidden()
    assert panel._rail.column_distance.isHidden()
    assert panel._rail.color_multi_sequence_background.isHidden()
    panel.deleteLater()
    bar.deleteLater()


def test_unaligned_overview_sequence_distance_defaults_to_two(
    qapp: QApplication, tmp_path: Path
) -> None:
    from PySide6.QtCore import QSettings

    settings = QSettings("RepeatViewer", "RepeatViewer")
    settings.beginGroup("annotation_overview_display")
    had_value = settings.contains("multi_sequence_distance")
    previous = settings.value("multi_sequence_distance")
    settings.remove("multi_sequence_distance")
    settings.endGroup()
    try:
        path = tmp_path / "multi.fa"
        path.write_text(">a\nACGT\n>b\nACGTACGT\n")
        panel = MultiSequenceOverviewPanel(path)
        panel.load_path(path, msa_mode=False)
        assert panel._rail.sequence_distance.value() == 2
        assert panel._bars.sequence_row_gap() == 2
        panel.deleteLater()
    finally:
        settings.beginGroup("annotation_overview_display")
        if had_value:
            settings.setValue("multi_sequence_distance", previous)
        else:
            settings.remove("multi_sequence_distance")
        settings.endGroup()


def test_msa_overview_sequence_distance_defaults_to_one(
    qapp: QApplication, tmp_path: Path
) -> None:
    from PySide6.QtCore import QSettings

    settings = QSettings("RepeatViewer", "RepeatViewer")
    settings.beginGroup("annotation_overview_display")
    had_value = settings.contains("sequence_distance")
    previous = settings.value("sequence_distance")
    settings.remove("sequence_distance")
    settings.endGroup()
    try:
        path = tmp_path / "aln.fa"
        path.write_text(">a\nACGT\n>b\nACGT\n")
        panel = MultiSequenceOverviewPanel(path)
        panel.load_path(path, msa_mode=True)
        assert panel._rail.sequence_distance.value() == 1
        assert panel._bars.sequence_row_gap() == 1
        panel.deleteLater()
    finally:
        settings.beginGroup("annotation_overview_display")
        if had_value:
            settings.setValue("sequence_distance", previous)
        else:
            settings.remove("sequence_distance")
        settings.endGroup()


def test_overview_auto_detects_aligned_msa(qapp: QApplication, tmp_path: Path) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">x\nAAAA----\n>y\nAAAATTTT\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path)  # msa_mode=None → auto
    assert panel.is_msa_mode() is True
    assert isinstance(panel._rail, SequenceAnnotationRail)
    panel.deleteLater()


def test_overview_annotation_tab_visible_and_wired(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "multi.fa"
    path.write_text(">a\nACGT\n>b\nTTTT\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=False)
    rail = panel._rail
    assert not rail._tab_buttons[SequenceAnnotationRail.TAB_ANNOTATIONS].isHidden()
    tip = rail.auto_ann_run_btn.toolTip().lower()
    assert "selected" in tip
    assert "msa" in tip or "open" in tip
    assert not rail._annotation_sections["manual"].isHidden()
    assert rail.show_annotation_track.isHidden()
    assert not rail.displayed_annotations_label.isHidden()
    assert not rail._annotation_kinds_host.isHidden()
    assert panel._auto_ann_checks
    assert callable(panel._on_auto_run_clicked)
    assert callable(panel._on_similarity_run_clicked)
    # No name selection → Feature detection / Similarity target the whole MSA.
    assert panel._selected_record_ids() == []
    assert panel._target_record_ids_for_msa_job() == ["a", "b"]
    tip_sim = rail.similarity_ann_run_btn.toolTip().lower()
    assert "selected" in tip_sim
    assert "msa" in tip_sim or "open" in tip_sim
    panel._bars.set_selected_rows([1])
    assert panel._target_record_ids_for_msa_job() == ["b"]
    # Regional MSA row box also scopes Feature detection / Similarity.
    panel._bars.set_regional_selection(1, 4, {0})
    assert panel._target_record_ids_for_msa_job() == ["a"]
    # Same control-panel width as single-sequence Annotation.
    assert rail._module_width == 260
    panel.deleteLater()


def test_overview_zoom_bar_matches_map_position_toolbar(
    qapp: QApplication, tmp_path: Path
) -> None:
    from PySide6.QtWidgets import QPushButton

    from repeatviewer.widgets.preview_controls import AnnotationZoomBar

    path = tmp_path / "aln.fa"
    path.write_text(">a\nACGT--GT\n>b\nACGTACGT\n")
    bar = AnnotationZoomBar()
    bar.set_overview_mode(True)
    assert not bar.pos_start.isHidden()
    assert bar.name_filter.isHidden()
    labels = {btn.text() for btn in bar._overview_host.findChildren(QPushButton)}
    assert "Select" in labels
    assert "Fit" in labels
    assert "Extract" not in labels

    panel = MultiSequenceOverviewPanel(path)
    panel.bind_annotation_zoom_bar(bar)
    panel.load_path(path, msa_mode=True)
    assert bar.pos_start.maximum() == 8
    assert bar.pos_end.value() == 8
    bar.pos_start.setValue(2)
    bar.pos_end.setValue(5)
    panel._select_by_position()
    assert panel._bars.pixels_per_base() > 0
    panel.deleteLater()
    bar.deleteLater()


def test_rail_tab_reclick_toggles_panel(qapp: QApplication) -> None:
    rail = SequenceAnnotationRail()
    rail.set_panel_collapsed(False, persist=False)
    rail.set_current_tab(SequenceAnnotationRail.TAB_DISPLAY)
    assert rail.is_panel_collapsed() is False
    assert rail._stack.isHidden() is False

    rail._on_tab_clicked(SequenceAnnotationRail.TAB_DISPLAY)
    assert rail.is_panel_collapsed() is True
    assert rail._stack.isHidden() is True
    assert rail.width() < 40

    rail._on_tab_clicked(SequenceAnnotationRail.TAB_DISPLAY)
    assert rail.is_panel_collapsed() is False
    assert rail._stack.isHidden() is False
    assert rail._stack.currentIndex() == SequenceAnnotationRail.TAB_DISPLAY
    rail.deleteLater()


def test_rail_switch_tab_while_expanded_keeps_panel(qapp: QApplication) -> None:
    rail = SequenceAnnotationRail()
    rail.set_panel_collapsed(False, persist=False)
    rail.set_current_tab(SequenceAnnotationRail.TAB_DISPLAY)
    rail._on_tab_clicked(SequenceAnnotationRail.TAB_FIND)
    assert rail.is_panel_collapsed() is False
    assert rail._stack.currentIndex() == SequenceAnnotationRail.TAB_FIND
    rail.deleteLater()


def test_msa_overview_find_tab_locates_motif(qapp: QApplication, tmp_path: Path) -> None:
    """Find next/prev on the overview rail must search multi/MSA FASTA."""
    path = tmp_path / "find_msa.fa"
    path.write_text(
        ">seqA\nACGTGTACGT\n"
        ">seqB\nTTTTTTTTTT\n"
        ">seqC\nCCCTGTCCCC\n"
    )
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=True)
    qapp.processEvents()
    assert panel._bars.is_msa_mode()
    rail = panel._rail
    rail.find_query.setText("TGT")
    rail.find_in_names.setChecked(False)
    rail.find_in_annotations.setChecked(False)
    rail.find_reverse_complement.setChecked(False)
    rail.findRequested.emit(1)
    qapp.processEvents()
    assert len(panel._find_items) >= 2
    assert panel._find_index == 0
    assert panel._bars._find_current is not None
    assert panel._bars._find_current[0] in {"seqA", "seqC"}
    first = panel._find_index
    rail.findRequested.emit(1)
    qapp.processEvents()
    assert panel._find_index != first or len(panel._find_items) == 1
    panel.close()


def test_msa_overview_find_annotation_selects_features_not_sequences(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Annotation-name Find selects matching features, not sequence rows."""
    from repeatviewer.backends.sequence_features import SequenceFeature

    path = tmp_path / "find_ann.fa"
    path.write_text(">seqA\nACGTACGTACGTACGTACGT\n>seqB\nACGTACGTACGTACGTACGT\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=True)
    qapp.processEvents()
    panel._features = {
        "seqA": [SequenceFeature(start=5, end=10, name="Gvch2", feature_type="domain")],
        "seqB": [SequenceFeature(start=2, end=4, name="other", feature_type="ORF")],
    }
    panel._bars.set_features_in_place(panel._features)
    panel._bars.set_selected_rows([1])  # prior sequence selection should be cleared
    qapp.processEvents()

    rail = panel._rail
    rail.find_query.setText("Gvch2")
    rail.find_in_sequence.setChecked(False)
    rail.find_in_names.setChecked(False)
    rail.find_in_annotations.setChecked(True)
    rail.findRequested.emit(1)
    qapp.processEvents()

    assert panel._find_items
    assert panel._find_items[0].kind == "annotation"
    assert panel._bars.selected_rows() == []
    assert panel._bars.selected_features() == [(0, 0)]
    # Annotation Find must not paint sequence underlays or a column band.
    assert panel._bars._find_hits == {}
    assert panel._bars.region_selection() is None
    panel.close()


def test_msa_overview_find_all_annotations_selects_every_match(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Find all on annotation names selects every match for bulk delete."""
    from repeatviewer.backends.sequence_features import SequenceFeature

    path = tmp_path / "find_all_ann.fa"
    path.write_text(
        ">seqA\nACGTACGTACGTACGTACGT\n"
        ">seqB\nACGTACGTACGTACGTACGT\n"
        ">seqC\nACGTACGTACGTACGTACGT\n"
    )
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=True)
    qapp.processEvents()
    panel._features = {
        "seqA": [SequenceFeature(start=2, end=5, name="LTR5", feature_type="LTR")],
        "seqB": [SequenceFeature(start=3, end=8, name="LTR3", feature_type="LTR")],
        "seqC": [SequenceFeature(start=1, end=4, name="ORF1", feature_type="ORF")],
    }
    panel._bars.set_features_in_place(panel._features)
    qapp.processEvents()

    rail = panel._rail
    rail.find_query.setText("LTR")
    rail.find_in_sequence.setChecked(False)
    rail.find_in_names.setChecked(False)
    rail.find_in_annotations.setChecked(True)
    rail.findRequested.emit(0)  # Find all
    qapp.processEvents()

    selected = set(panel._bars.selected_features())
    assert selected == {(0, 0), (1, 0)}
    assert panel._bars._find_hits == {}
    assert panel._bars.region_selection() is None
    panel.close()


def test_msa_overview_context_menu_has_rename_for_sequence_selection(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Right-click menu: Rename when a sequence row is selected."""
    path = tmp_path / "ctx.fa"
    path.write_text(">seqA\nACGTACGT\n>seqB\nGGGGGGGG\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=True)
    panel._bars.set_selected_rows([0])
    qapp.processEvents()

    menu = panel._build_overview_context_menu()
    texts = [a.text() for a in menu.actions() if not a.isSeparator()]
    assert "Add feature…" in texts
    assert "Rename sequence" in texts
    assert "Set as reference" in texts
    assert "Copy Sequence" in texts
    assert "Extract Sequence" in texts
    assert "Copy Sequence without gap" in texts
    assert "Extract Sequence without gap" in texts
    assert "Show self-dotplot" in texts
    # Copy grouped above its without-gap twin; same for Extract.
    assert texts.index("Copy Sequence") < texts.index("Copy Sequence without gap")
    assert texts.index("Copy Sequence without gap") < texts.index("Extract Sequence")
    assert texts.index("Extract Sequence") < texts.index(
        "Extract Sequence without gap"
    )
    rename = next(a for a in menu.actions() if a.text() == "Rename sequence")
    assert rename.isEnabled() is True
    set_std = next(a for a in menu.actions() if a.text() == "Set as reference")
    assert set_std.isEnabled() is True
    show_dot = next(a for a in menu.actions() if a.text() == "Show self-dotplot")
    assert show_dot.isEnabled() is True

    # Region-only selection: no Rename / Set as reference; self-dotplot disabled.
    panel._bars.set_selected_rows([])
    panel._bars.set_region_selection(2, 4)
    menu2 = panel._build_overview_context_menu()
    texts2 = [a.text() for a in menu2.actions() if not a.isSeparator()]
    assert "Rename sequence" not in texts2
    assert "Set as reference" not in texts2
    assert "Copy Sequence without gap" in texts2
    show_dot2 = next(a for a in menu2.actions() if a.text() == "Show self-dotplot")
    assert show_dot2.isEnabled() is False
    panel.close()


def test_msa_color_vs_standard_defaults_to_first_until_pinned(
    qapp: QApplication, tmp_path: Path
) -> None:
    from repeatviewer.backends.msa_coordinates import MSA_COLOR_DIFF_FROM_SELECTED

    path = tmp_path / "std.fa"
    path.write_text(">seqA\nACGT\n>seqB\nTGCA\n>seqC\nAAAA\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=True)
    panel._bars.set_msa_color_mode(MSA_COLOR_DIFF_FROM_SELECTED)
    panel._sync_standard_name_emphasis()
    assert panel._bars.msa_diff_standard_row() == 0
    assert panel._bars.msa_diff_standard_is_pinned() is False
    assert panel._names.standard_row() == 0

    panel._bars.set_selected_rows([2])
    assert panel._bars.msa_diff_standard_row() == 0

    panel._bars.set_selected_rows([1])
    panel._set_selected_as_msa_standard()
    assert panel._bars.msa_diff_standard_is_pinned() is True
    assert panel._bars.msa_diff_standard_id() == "seqB"
    assert panel._bars.msa_diff_standard_row() == 1
    assert panel._names.standard_row() == 1

    panel._bars.set_selected_rows([0])
    assert panel._bars.msa_diff_standard_row() == 1

    panel._bars.set_msa_color_mode("all_backgrounds")
    panel._rail.show_standard_boundaries.setChecked(True)
    panel._sync_standard_name_emphasis()
    # Show standard boundaries still emphasizes the name without Color vs standard.
    assert panel._names.standard_row() == 1
    panel._rail.show_standard_boundaries.setChecked(False)
    panel._sync_standard_name_emphasis()
    assert panel._names.standard_row() is None
    panel._bars.set_msa_color_mode(MSA_COLOR_DIFF_FROM_SELECTED)
    panel._sync_standard_name_emphasis()
    assert panel._bars.msa_diff_standard_row() == 1
    assert panel._names.standard_row() == 1
    panel.close()

    # Reopen the same MSA — pinned standard is restored from settings.
    panel2 = MultiSequenceOverviewPanel(path)
    panel2.load_path(path, msa_mode=True)
    panel2._bars.set_msa_color_mode(MSA_COLOR_DIFF_FROM_SELECTED)
    assert panel2._bars.msa_diff_standard_id() == "seqB"
    assert panel2._bars.msa_diff_standard_row() == 1
    panel2.close()


def test_msa_color_vs_standard_conserved_diff_uses_default_palette(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Conserved diffs → Default colours; rare diffs + standard → All black."""
    from repeatviewer.backends.msa_coordinates import MSA_COLOR_DIFF_FROM_SELECTED
    from repeatviewer.backends.msa_plot import (
        DEFAULT_NUCLEOTIDE_COLORS,
        NUCLEOTIDE_COLORS_ALL_BLACK,
    )

    path = tmp_path / "vs_std.fa"
    # Col1: majority T vs standard A → conserved T diff.
    # Col2: single T vs standard A → rare T diff.
    path.write_text(
        ">std\nAAA\n"
        ">s2\nATT\n"
        ">s3\nATA\n"
        ">s4\nATA\n"
        ">s5\nATA\n"
    )
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=True)
    bars = panel._bars
    bars.set_msa_color_mode(MSA_COLOR_DIFF_FROM_SELECTED)
    bars.set_msa_divergence_threshold(0.4)
    assert bars.msa_diff_standard_row() == 0

    std_a = bars._msa_residue_paint_pair(0, "A", row=0)
    assert std_a is not None
    assert std_a[0].name().lower() == NUCLEOTIDE_COLORS_ALL_BLACK["A"].lower()

    conserved_t = bars._msa_residue_paint_pair(1, "T", row=1)
    assert conserved_t is not None
    assert conserved_t[0].name().lower() == DEFAULT_NUCLEOTIDE_COLORS["T"].lower()

    rare_t = bars._msa_residue_paint_pair(2, "T", row=1)
    assert rare_t is not None
    assert rare_t[0].name().lower() == NUCLEOTIDE_COLORS_ALL_BLACK["T"].lower()

    # Match to standard: black letter colour (pair[0] used when unfilled).
    match_a = bars._msa_residue_paint_pair(0, "A", row=1)
    assert match_a is not None
    assert match_a[0].name().lower() == "#000000"

    # Gap that differs from standard stays the red accent fill.
    gap = bars._msa_overview_fill_color(0, "-", row=1)
    assert gap is not None
    assert gap.name().lower() in {"#e74c3c", "#c0392b"}

    # N/X always use the fixed yellow highlight under Color vs reference.
    n_fill = bars._msa_residue_paint_pair(0, "N", row=1)
    assert n_fill is not None
    assert n_fill[0].name().lower() == "#f0c000"
    assert bars._msa_should_fill_background(0, "N", row=1) is True
    panel.close()


def test_msa_standard_extent_for_boundary_guides(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Standard-boundary guides use leftmost/rightmost nucleotides in any scheme."""
    from repeatviewer.backends.msa_coordinates import (
        MSA_COLOR_BG_WHEN_CONSERVED,
        MSA_COLOR_DIFF_FROM_SELECTED,
    )

    path = tmp_path / "islands.fa"
    path.write_text(">std\nAC----GTAC\n>other\nACGTACGTAC\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=True)
    panel._rail.show_standard_boundaries.setChecked(True)
    panel._bars.set_msa_color_mode(MSA_COLOR_DIFF_FROM_SELECTED)
    assert panel._bars.msa_diff_standard_row() == 0
    # First non-gap col 0, last non-gap col 9 — ignore the internal gap.
    assert panel._bars.msa_standard_extent_columns() == (0, 9)
    assert panel._bars.show_msa_standard_boundaries() is True
    panel._sync_scale_standard_boundaries()
    assert panel._scale._standard_boundary_cols == (0, 9)
    # Guides stay available under Color conserved (and other schemes).
    panel._bars.set_msa_color_mode(MSA_COLOR_BG_WHEN_CONSERVED)
    assert panel._bars.msa_standard_extent_columns() == (0, 9)
    panel._sync_scale_standard_boundaries()
    assert panel._scale._standard_boundary_cols == (0, 9)
    assert panel._rail.show_standard_boundaries.text() == "Show reference boundaries"
    panel._rail.show_standard_boundaries.setChecked(False)
    assert panel._bars.show_msa_standard_boundaries() is False
    assert panel._scale._standard_boundary_cols is None
    panel.close()


def test_show_reference_boundaries_persists_across_new_msa(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Unchecking Show reference boundaries must stay off on the next MSA."""
    from PySide6.QtCore import QSettings

    a = tmp_path / "a.fa"
    b = tmp_path / "b.fa"
    a.write_text(">x\nACGT\n>y\nACGT\n")
    b.write_text(">p\nTTTT\n>q\nTTTT\n")
    first = MultiSequenceOverviewPanel(a)
    first.load_path(a, msa_mode=True)
    first._rail.show_standard_boundaries.setChecked(False)
    first._save_overview_display_settings()
    first.close()

    s = QSettings("RepeatViewer", "RepeatViewer")
    s.beginGroup("annotation_overview_display")
    assert s.value("show_standard_boundaries", True, type=bool) is False
    s.endGroup()

    second = MultiSequenceOverviewPanel(b)
    second.load_path(b, msa_mode=True)
    assert second._rail.show_standard_boundaries.isChecked() is False
    assert second._bars.show_msa_standard_boundaries() is False
    # Restore default for other tests in this process.
    second._rail.show_standard_boundaries.setChecked(True)
    second._save_overview_display_settings()
    second.close()


def test_msa_copy_sequence_without_gap_strips_dashes(
    qapp: QApplication, tmp_path: Path, monkeypatch
) -> None:
    path = tmp_path / "gaps.fa"
    path.write_text(">seqA\nAC--GT\n>seqB\nACGTGT\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=True)
    panel._bars.set_selected_rows([0])
    clipped: list[str] = []
    monkeypatch.setattr(
        QApplication,
        "clipboard",
        lambda: type("Clip", (), {"setText": lambda _s, t: clipped.append(t)})(),
    )
    panel._copy_selected_sequences(without_gaps=True)
    assert clipped == [">seqA\nACGT"]
    panel.close()


def test_msa_overview_rename_sequence_writes_fasta(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "rename.fa"
    path.write_text(">oldName\nACGTACGT\n>other\nTTTTTTTT\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=True)
    panel._bars.set_selected_rows([0])
    with patch(
        "repeatviewer.widgets.file_conflict.prompt_rename_text",
        return_value=("newName", True),
    ):
        panel._rename_selected_sequence()
    text = path.read_text()
    assert ">newName\n" in text
    assert ">oldName\n" not in text
    assert panel.record_ids()[0] == "newName"
    panel.close()


def test_overview_splitter_keeps_names_flush_with_tracks(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Widening the names|tracks splitter must not leave a blank gap."""
    path = tmp_path / "aln.fa"
    path.write_text(">seqA\nACGTACGTAC\n>seqB\nACGTACGTAC\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.resize(900, 400)
    panel.show()
    panel.load_path(path, msa_mode=True)
    qapp.processEvents()
    panel._fit_bars()
    qapp.processEvents()

    total = sum(panel._splitter.sizes()) or 800
    panel._splitter.setSizes([360, max(200, total - 360)])
    panel._on_splitter_moved()
    qapp.processEvents()

    assert panel._names.width() == panel._names_viewport_width()
    panel.close()
