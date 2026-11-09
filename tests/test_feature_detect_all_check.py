"""Feature detection "All" ticks every detector, TSD and PBS included."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication

from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel


def test_all_ticks_tsd_and_pbs(qapp: QApplication, tmp_path: Path) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">a\nACGTACGT\n>b\nACGT--GT\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=True)
    rail = panel._rail
    labels = {"LTR", "TIR", "Tandem", "ORF", "TSD", "PBS"}
    rail.auto_ann_all_check.setChecked(False)
    rail.auto_ann_all_check.setChecked(True)
    assert all(rail.auto_ann_checks[k].isChecked() for k in labels)
    rail.auto_ann_checks["PBS"].setChecked(False)
    assert not rail.auto_ann_all_check.isChecked()
    rail.auto_ann_checks["PBS"].setChecked(True)
    assert rail.auto_ann_all_check.isChecked()
    panel.deleteLater()
