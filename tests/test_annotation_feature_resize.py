"""Annotation track: edge drag updates feature start/end bounds."""

from __future__ import annotations

from PySide6.QtCore import QPointF
from PySide6.QtWidgets import QApplication

from repeatviewer.backends.sequence_features import SequenceFeature
from repeatviewer.widgets.sequence_map_viewer import (
    AnnotationTrackWidget,
    _feature_with_bounds,
)


def test_feature_with_bounds_keeps_metadata() -> None:
    feat = SequenceFeature(
        start=10,
        end=20,
        name="LTR5 (90%)",
        feature_type="LTR",
        color="#c0392b",
        strand="-",
        track=2,
        source="repeatviewer_internal",
        pident=90.0,
    )
    updated = _feature_with_bounds(feat, 5, 25)
    assert (updated.start, updated.end) == (5, 25)
    assert updated.name == feat.name
    assert updated.strand == "-"
    assert updated.track == 2
    assert updated.pident == 90.0


def test_overview_edge_hit_and_resize(qapp: QApplication, tmp_path) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    path = tmp_path / "aln.fa"
    path.write_text(">a\nACGTACGTACGTACGTACGT\n>b\nACGTACGTACGTACGTACGT\n")
    canvas = _BarsCanvas()
    feat = SequenceFeature(start=5, end=12, name="LTR1", feature_type="LTR", strand="+")
    canvas.set_msa_data(
        ["a", "b"],
        ["ACGTACGTACGTACGTACGT", "ACGTACGTACGTACGTACGT"],
        features={"a": [feat], "b": []},
        fit_bars_width=400,
    )
    canvas.set_pixels_per_base(10.0)
    canvas.set_show_annotation_track(True)
    canvas.resize(500, 200)
    canvas.show()

    rect = canvas._feature_rect(0, 0)
    assert rect is not None
    y = rect.center().y()
    assert canvas._edge_hit(0, rect.left(), y) == (0, "start")
    assert canvas._edge_hit(0, rect.right(), y) == (0, "end")
    assert canvas._edge_hit(0, rect.center().x(), y) is None

    resolved = canvas._store_index_for_flat(0, 0)
    assert resolved is not None
    rid, store_i, _ = resolved
    canvas._resize_active = True
    canvas._resize_rid = rid
    canvas._resize_store_index = store_i
    canvas._resize_flat = 0
    canvas._resize_row = 0
    canvas._resize_edge = "end"
    canvas._resize_orig_start = 5
    canvas._resize_orig_end = 12
    # bp 16 → x = pad + 15 * ppb
    canvas._apply_resize_at(canvas._bar_pad + 15 * canvas._ppb)
    assert canvas._features["a"][0].start == 5
    assert canvas._features["a"][0].end == 16
    canvas.deleteLater()


def test_edge_hit_and_resize_drag(qapp: QApplication) -> None:
    track = AnnotationTrackWidget()
    feat = SequenceFeature(start=10, end=20, name="ORF1", feature_type="ORF", strand="+")
    track.set_view(
        features=[feat],
        feature_rows=[0],
        scroll_px=0.0,
        ppb=10.0,
        viewport_w=400,
        selected_indices={0},
        seq_length=100,
    )
    # Feature spans x = left_margin + (10-1)*10 .. left_margin + 20*10
    from repeatviewer.widgets.sequence_map_viewer import _LEFT_MARGIN

    left = _LEFT_MARGIN + 9 * 10
    right = _LEFT_MARGIN + 20 * 10
    y = track._feature_rect(0).center().y()

    assert track._edge_hit(QPointF(left, y)) == (0, "start")
    assert track._edge_hit(QPointF(right, y)) == (0, "end")
    assert track._edge_hit(QPointF((left + right) / 2, y)) is None

    track._resize_active = True
    track._resize_idx = 0
    track._resize_edge = "end"
    track._resize_orig_start = 10
    track._resize_orig_end = 20
    # Drag end to bp 30 → x = margin + 29*10 (bp = (x-margin)/ppb + 1)
    track._apply_resize_at(_LEFT_MARGIN + 29 * 10)
    assert track._features[0].start == 10
    assert track._features[0].end == 30

    track._resize_edge = "start"
    track._resize_orig_start = 10
    track._resize_orig_end = 30
    track._apply_resize_at(_LEFT_MARGIN + 14 * 10)  # bp 15
    assert track._features[0].start == 15
    assert track._features[0].end == 30

    track.deleteLater()
