"""Overview toolbar: Cov hidden; "Ann" shows/hides all annotation tracks."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication

from repeatviewer.backends.sequence_features import SequenceFeature
from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel
from repeatviewer.widgets.preview_controls import AnnotationZoomBar


def test_ann_button_toggles_all_tracks_and_cov_hidden(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">a\nACGTACGTAC\n>b\nACGTACGTAC\n")
    panel = MultiSequenceOverviewPanel(path)
    bar = AnnotationZoomBar()
    bar.set_overview_mode(True)
    panel.bind_annotation_zoom_bar(bar)
    panel.load_path(path, msa_mode=True)
    panel._features = {
        "a": [SequenceFeature(start=2, end=5, name="ORF1", feature_type="ORF")],
        "b": [],
    }
    panel._refresh_features_after_edit()
    # Saved "Displayed annotations" picks may hide ORF; show every type.
    panel._rail.annotation_kinds_all_check.setChecked(False)
    panel._rail.annotation_kinds_all_check.setChecked(True)
    assert bar._coverage_btn.isHidden()
    assert bar._annotations_btn.isChecked()
    shown_h = panel._bars.row_heights()[0]

    bar._annotations_btn.setChecked(False)
    assert not panel._bars._show_annotation_track
    assert panel._bars.row_heights()[0] < shown_h
    # Per-type picks are untouched.
    assert panel._rail.show_annotation_track.isChecked()

    bar._annotations_btn.setChecked(True)
    assert panel._bars._show_annotation_track
    assert panel._bars.row_heights()[0] == shown_h
    panel.deleteLater()
    bar.deleteLater()
