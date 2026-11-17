"""Feature detection "All" ticks every detector, TSD and PBS included."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication

from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel
from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail


def test_all_ticks_tsd_and_pbs(qapp: QApplication, tmp_path: Path) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">a\nACGTACGT\n>b\nACGT--GT\n")
    panel = MultiSequenceOverviewPanel(path)
    panel.load_path(path, msa_mode=True)
    rail = panel._rail
    assert not rail.auto_ann_checks["TSD"].isHidden()
    labels = {"LTR", "TIR", "Tandem", "ORF", "TSD", "PBS"}
    rail.auto_ann_all_check.setChecked(False)
    rail.auto_ann_all_check.setChecked(True)
    assert all(rail.auto_ann_checks[k].isChecked() for k in labels)
    rail.auto_ann_checks["PBS"].setChecked(False)
    assert not rail.auto_ann_all_check.isChecked()
    rail.auto_ann_checks["PBS"].setChecked(True)
    assert rail.auto_ann_all_check.isChecked()
    panel.deleteLater()


def test_tsd_hidden_unless_msa(qapp: QApplication, tmp_path: Path) -> None:
    """TSD is MSA-only — hide on single-seq and multi unaligned overviews."""
    rail = SequenceAnnotationRail()
    assert rail.auto_ann_checks["TSD"].isHidden()
    rail.set_overview_display_mode(True)
    rail.set_msa_display_mode(False)
    assert rail.auto_ann_checks["TSD"].isHidden()
    rail.set_msa_display_mode(True)
    assert not rail.auto_ann_checks["TSD"].isHidden()
    rail.deleteLater()

    multi = tmp_path / "multi.fa"
    multi.write_text(">a\nACGTACGTAC\n>b\nTTTT\n")
    panel = MultiSequenceOverviewPanel(multi)
    panel.load_path(multi, msa_mode=False)
    assert panel._rail.auto_ann_checks["TSD"].isHidden()
    panel.deleteLater()


def test_pbs_aligns_with_ltr_when_tsd_hidden(qapp: QApplication) -> None:
    """Hidden TSD must not indent PBS on the wrapped chip row."""
    rail = SequenceAnnotationRail()
    rail.resize(280, 600)
    rail.show()
    qapp.processEvents()
    rail.set_msa_display_mode(False)
    qapp.processEvents()
    ltr = rail.auto_ann_checks["LTR"]
    pbs = rail.auto_ann_checks["PBS"]
    assert rail.auto_ann_checks["TSD"].isHidden()
    assert ltr.x() == pbs.x()
    rail.deleteLater()
