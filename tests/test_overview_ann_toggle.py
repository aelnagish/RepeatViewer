"""Overview toolbar: Cov hidden; "Feat" shows/hides all feature tracks."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication

from repeatviewer.backends.sequence_features import SequenceFeature
from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel
from repeatviewer.widgets.preview_controls import AnnotationZoomBar


def test_feat_button_toggles_all_tracks_and_cov_hidden(
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
    assert bar._annotations_btn.text() == "Feat"
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


def test_reveal_features_overrides_feat_and_display_off(
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
        "a": [SequenceFeature(start=2, end=5, name="LTR1", feature_type="LTR")],
        "b": [],
    }
    panel._refresh_features_after_edit()
    panel._rail.ensure_annotation_kinds_visible(["LTR"])
    panel._rail.annotation_kind_checks["LTR"].setChecked(False)
    bar._annotations_btn.setChecked(False)
    assert not panel._bars._show_annotation_track

    panel._reveal_features_on_canvas(["LTR"])
    assert bar._annotations_btn.isChecked()
    assert not panel._annotations_hidden
    assert panel._rail.annotation_kind_checks["LTR"].isChecked()
    assert panel._bars._show_annotation_track
    panel.deleteLater()
    bar.deleteLater()
