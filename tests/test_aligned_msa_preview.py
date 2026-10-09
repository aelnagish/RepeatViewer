"""Aligned FASTA detection and MSA residue-row overview load."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication

from repeatviewer.backends.sequence_features import (
    is_aligned_fasta,
    load_aligned_fasta_rows,
)
from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel
from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail


def test_is_aligned_fasta_equal_length_gapped(tmp_path: Path) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">s1\nACGT--GT\n>s2\nACGTACGT\n")
    assert is_aligned_fasta(path) is True


def test_is_aligned_fasta_false_for_unequal_lengths(tmp_path: Path) -> None:
    path = tmp_path / "multi.fa"
    path.write_text(">s1\nACGT\n>s2\nACGTACGT\n")
    assert is_aligned_fasta(path) is False


def test_is_aligned_fasta_false_for_single_record(tmp_path: Path) -> None:
    path = tmp_path / "one.fa"
    path.write_text(">only\nACGTACGT\n")
    assert is_aligned_fasta(path) is False


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


def test_overview_load_path_unaligned_stays_bars(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "multi.fa"
    path.write_text(">a\nACGT\n>b\nACGTACGTACGT\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=False)
    assert panel.is_msa_mode() is False
    assert panel._bars.max_length() == 12
    panel.deleteLater()


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
    assert not rail.show_annotation_track.isHidden()
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
