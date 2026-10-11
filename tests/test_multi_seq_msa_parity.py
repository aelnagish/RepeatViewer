"""Multi-seq / MSA Annotation parity: sizes, selection, tracks, row height."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest
from PySide6.QtCore import QPoint, QPointF, QRect, Qt
from PySide6.QtGui import QMouseEvent, QPalette
from PySide6.QtWidgets import QApplication, QMessageBox

from repeatviewer.backends.sequence_features import SequenceFeature
from repeatviewer.widgets.preview_controls import (
    ANNOTATION_FIT_BUTTON_RENDERED_WIDTH,
    ANNOTATION_SELECT_BUTTON_RENDERED_WIDTH,
    ANNOTATION_TOOLBAR_RENDERED_HEIGHT,
    ZOOM_SLIDER_HEIGHT,
    AnnotationZoomBar,
    PreviewHeaderBar,
    PreviewZoomBar,
    _OVERVIEW_POSITION_SPIN_WIDTH,
    _OVERVIEW_SELECT_FIT_HEIGHT,
    _ZOOM_BTN_SIDE,
)


def test_overview_zoom_bar_matches_single_seq_sizes(qapp: QApplication) -> None:
    bar = AnnotationZoomBar()
    bar.set_overview_mode(True)
    bar.show()
    qapp.processEvents()
    assert bar.height() == ANNOTATION_TOOLBAR_RENDERED_HEIGHT == 26
    assert bar._zoom_slider.height() == ZOOM_SLIDER_HEIGHT == 22
    assert bar._pos_start.width() == _OVERVIEW_POSITION_SPIN_WIDTH == 108
    assert _ZOOM_BTN_SIDE == ZOOM_SLIDER_HEIGHT
    # Actions lives on the Preview title bar in overview mode (not this bar).
    assert bar._actions_btn.isHidden()
    # Select / Fit / Clean / Cons live in overview host as QPushButtons with fixed height 20.
    from PySide6.QtWidgets import QPushButton

    labels = [
        w.text()
        for w in bar._overview_host.findChildren(QPushButton)
        if w.text() in {"Select", "Fit", "Clean", "Cons"}
    ]
    assert labels == ["Select", "Fit", "Clean", "Cons"]
    heights = [
        w.height()
        for w in bar._overview_host.findChildren(QPushButton)
        if w.text() in {"Select", "Fit", "Clean", "Cons"}
    ]
    assert heights == [_OVERVIEW_SELECT_FIT_HEIGHT] * 4
    bar.close()


def test_msa_toolbar_matches_rendered_single_sequence_control_sizes(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Compare final styled geometry, not only pre-stylesheet fixed values."""
    from PySide6.QtWidgets import QPushButton, QToolButton

    from repeatviewer.ui.theme import application_stylesheet
    from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel
    from repeatviewer.widgets.sequence_map_viewer import SequenceMapWidget

    previous_stylesheet = qapp.styleSheet()
    qapp.setStyleSheet(application_stylesheet(toolbar_font_pt=12))
    single = SequenceMapWidget(embedded=True)
    msa = AnnotationZoomBar()
    fasta = tmp_path / "msa.fa"
    fasta.write_text(">a\nACGT\n>b\nACGT\n", encoding="utf-8")
    overview = MultiSequenceOverviewPanel(fasta)
    try:
        overview.bind_annotation_zoom_bar(msa)
        single.show()
        overview.show()
        qapp.processEvents()

        single_minus = next(
            w for w in single._toolbar.findChildren(QToolButton) if w.text() == "−"
        )
        msa_minus = next(
            w for w in msa._overview_host.findChildren(QToolButton) if w.text() == "−"
        )
        single_select = next(
            w
            for w in single._toolbar.findChildren(QPushButton)
            if w.text() == "Select"
        )
        single_fit = next(
            w for w in single._toolbar.findChildren(QPushButton) if w.text() == "Fit"
        )

        # Shared Select / Fit widths recorded from multi-seq Overview.
        assert single._pos_start.size().width() == 108
        assert single._pos_start.size().height() == 20
        assert single_minus.size().width() == 22
        assert single_minus.size().height() == 22
        assert single_select.height() == 22
        assert (
            single_select.width()
            == msa._select_btn.width()
            == ANNOTATION_SELECT_BUTTON_RENDERED_WIDTH
            == 55
        )
        assert (
            single_fit.width()
            == msa._fit_btn.width()
            == ANNOTATION_FIT_BUTTON_RENDERED_WIDTH
            == 37
        )

        assert msa._pos_start.size() == single._pos_start.size()
        assert msa._pos_end.size() == single._pos_end.size()
        assert msa_minus.size() == single_minus.size()
        assert msa._select_btn.height() == single_select.height()
        assert msa._fit_btn.height() == single_select.height()
        assert msa._clean_btn.height() == single_select.height()
        assert msa._consensus_btn.height() == single_select.height()
        assert (
            single._toolbar.height()
            == overview._toolbar_slot.height()
            == msa.height()
            == ANNOTATION_TOOLBAR_RENDERED_HEIGHT
            == 26
        )
        assert msa_minus.mapTo(overview._toolbar_slot, QPoint(0, 0)) == (
            single_minus.mapTo(single._toolbar, QPoint(0, 0))
        ) == QPoint(6, 2)
        assert overview._toolbar_slot.palette().color(
            QPalette.ColorRole.Window
        ) == single._toolbar.palette().color(QPalette.ColorRole.Window)
    finally:
        overview._unbind_header_zoom_bar()
        overview.close()
        single.close()
        msa.close()
        qapp.setStyleSheet(previous_stylesheet)


def test_overview_mounts_zoom_bar_beside_rail(qapp: QApplication, tmp_path: Path) -> None:
    from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel

    fasta = tmp_path / "m.fa"
    fasta.write_text(">a\nACGT\n>b\nTTTT\n", encoding="utf-8")
    panel = MultiSequenceOverviewPanel(fasta)
    bar = AnnotationZoomBar()
    panel.show()
    bar.show()
    qapp.processEvents()
    panel.bind_annotation_zoom_bar(bar)
    qapp.processEvents()
    assert bar.parent() is panel._toolbar_slot
    assert not panel._toolbar_slot.isHidden()
    assert bar._actions_btn.isHidden()
    # Toolbar spans names + tracks and its first control aligns with name text.
    from PySide6.QtCore import QPoint
    from PySide6.QtWidgets import QPushButton, QToolButton

    slot_x = panel._toolbar_slot.mapTo(panel, QPoint(0, 0)).x()
    assert slot_x == 0
    zoom_out = next(
        w for w in bar._overview_host.findChildren(QToolButton) if w.text() == "−"
    )
    assert abs(zoom_out.mapTo(panel, QPoint(0, 0)).x() - 6) <= 1
    select = next(
        w
        for w in bar._overview_host.findChildren(QPushButton)
        if w.text() == "Select"
    )
    assert select.size().height() == _OVERVIEW_SELECT_FIT_HEIGHT
    assert select.size().width() <= 80
    panel._unbind_header_zoom_bar()
    assert panel._toolbar_slot.isHidden()
    panel.close()
    bar.close()


def test_header_refresh_keeps_overview_zoom_bar_visible(
    qapp: QApplication, tmp_path: Path
) -> None:
    """A late preview-header sync must not hide the MSA's mounted toolbar."""
    from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel

    fasta = tmp_path / "msa.fa"
    fasta.write_text(">a\nACGT\n>b\nACGT\n", encoding="utf-8")
    preview_bar = PreviewZoomBar()
    bar = AnnotationZoomBar()
    header = PreviewHeaderBar(preview_bar, annotation_zoom_bar=bar)
    panel = MultiSequenceOverviewPanel(fasta)
    panel.show()
    panel.bind_annotation_zoom_bar(bar)
    qapp.processEvents()

    assert bar.parent() is panel._toolbar_slot
    assert not bar.isHidden()

    # This is the refresh that used to call bar.setVisible(False) after the
    # Overview had already reparented and shown the shared control bar.
    header.set_header_zoom_mode("preview")
    qapp.processEvents()

    assert bar.parent() is panel._toolbar_slot
    assert not bar.isHidden()
    assert not panel._toolbar_slot.isHidden()

    panel._unbind_header_zoom_bar()
    panel.close()
    header.close()
    preview_bar.close()
    bar.close()


def test_overview_rail_keeps_add_transfer_and_track_toggle(qapp: QApplication) -> None:
    from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail

    rail = SequenceAnnotationRail()
    rail.show()
    qapp.processEvents()
    rail.set_overview_display_mode(True)
    assert not rail.show_annotation_track.isHidden()
    assert rail.show_annotation_track.isChecked()
    # Sort options: multi-seq has original/name; MSA also has selected columns.
    assert not rail.sequence_sort_row.isHidden()
    assert rail.sequence_sort.findData("selected_columns") < 0
    rail.set_msa_display_mode(True)
    assert rail.sequence_sort.findData("selected_columns") >= 0
    assert rail.show_complement.isHidden()
    sections = rail._annotation_sections
    assert not sections["manual"].isHidden()
    assert not sections["transfer"].isHidden()
    tip = rail.auto_ann_run_btn.toolTip().lower()
    assert "selected" in tip
    assert "ctrl+a" in tip or "cmd+a" in tip
    assert not rail.delete_btn.isHidden()
    assert not rail.restore_btn.isHidden()
    assert not rail.export_btn.isHidden()
    rail.close()


def test_msa_overview_always_paints_gap_backbone(qapp: QApplication) -> None:
    from PySide6.QtCore import QPoint
    from PySide6.QtGui import QColor, QImage, QPainter

    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    canvas.resize(400, 80)
    seq = "AAAA" + ("-" * 20) + "TTTT"
    canvas.set_msa_data(["a"], [seq])
    canvas.set_pixels_per_base(4.0)
    assert not canvas._msa_shows_letters()

    img = QImage(canvas.size(), QImage.Format.Format_ARGB32)
    img.fill(QColor("white"))
    painter = QPainter(img)
    canvas.render(painter, QPoint())
    painter.end()

    x = int(canvas.bar_pad() + 14 * canvas.pixels_per_base())
    y = int(canvas.row_tops()[0] + canvas._msa_letter_band_height() * 0.5)
    pix = img.pixelColor(x, y)
    assert pix.red() < 40 and pix.green() < 40 and pix.blue() < 40


def test_msa_sort_by_name_and_original(qapp: QApplication) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    canvas.set_msa_data(
        ["c_seq", "a_seq", "b_seq"],
        ["AAA", "CCC", "GGG"],
    )
    canvas.sort_rows("name")
    assert canvas.record_ids() == ["a_seq", "b_seq", "c_seq"]
    canvas.sort_rows("original")
    assert canvas.record_ids() == ["c_seq", "a_seq", "b_seq"]


def test_msa_sort_by_selected_columns(qapp: QApplication) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    canvas.set_msa_data(
        ["r1", "r2", "r3"],
        ["AAA", "CCC", "AAC"],
    )
    canvas.set_region_selection(1, 2)
    canvas.sort_rows("selected_columns")
    assert canvas.record_ids()[0] == "r1"
    assert canvas.record_ids()[1] == "r3"
    assert canvas.record_ids()[2] == "r2"


def test_sorting_label_and_selected_column_is_oneshot(
    qapp: QApplication, tmp_path: Path
) -> None:
    from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel
    from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail

    rail = SequenceAnnotationRail()
    rail.set_overview_display_mode(True)
    rail.set_msa_display_mode(True)
    assert rail.sequence_sort_label.text() == "Sorting:"
    assert rail.sequence_sort.findText("Selected column") >= 0
    assert rail.sequence_sort.findText("Selected columns") < 0
    rail.close()

    fasta = tmp_path / "sort_oneshot.fa"
    fasta.write_text(">r1\nAAA\n>r2\nCCC\n>r3\nAAC\n")
    panel = MultiSequenceOverviewPanel(fasta)
    panel.load_path(fasta, msa_mode=True)
    panel._rail.set_sequence_sort_mode("selected_columns")
    panel._bars.set_region_selection(1, 2)
    panel._bars.sort_rows("selected_columns")
    order_after_sort = list(panel._bars.record_ids())
    # Column click / region change must not re-sort while mode stays selected.
    panel._on_region_selection_changed((3, 3))
    assert panel._bars.record_ids() == order_after_sort
    panel.deleteLater()


def test_bars_select_all_and_region_without_zoom(qapp: QApplication) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    canvas.resize(600, 200)
    canvas.show()
    qapp.processEvents()
    canvas.set_data(
        ["a", "b", "c"],
        {"a": 500, "b": 400, "c": 300},
        {},
        fit_bars_width=500,
    )
    zoom_calls: list[tuple[float, float]] = []
    canvas.regionZoomRequested.connect(lambda s, e: zoom_calls.append((s, e)))

    canvas.select_all_rows()
    assert canvas.selected_rows() == [0, 1, 2]

    canvas.set_region_selection(10, 40)
    assert canvas.region_selection() == (10, 40)
    assert zoom_calls == []

    # Simulated rubber-band release path: set region, no zoom emit.
    canvas._drag_origin = QPointF(20, 10)
    canvas._rubber_active = True
    canvas._rubber = canvas.rect().adjusted(30, 5, -200, -100)
    release = QMouseEvent(
        QMouseEvent.Type.MouseButtonRelease,
        QPointF(120, 40),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    canvas.mouseReleaseEvent(release)
    qapp.processEvents()
    assert canvas.region_selection() is not None
    assert zoom_calls == []

    # Click inside (or outside) the selected band clears it.
    sel = canvas.region_selection()
    assert sel is not None
    inside_x = canvas.bp_to_x((sel[0] + sel[1]) / 2.0)
    canvas._drag_origin = QPointF(inside_x, 20)
    canvas._rubber_active = False
    canvas._rubber = QRect()
    click = QMouseEvent(
        QMouseEvent.Type.MouseButtonRelease,
        QPointF(inside_x, 20),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    canvas.mouseReleaseEvent(click)
    qapp.processEvents()
    assert canvas.region_selection() is None
    canvas.close()


def test_msa_residue_size_and_annotation_track(qapp: QApplication) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    canvas.show()
    qapp.processEvents()
    feats = {
        "s1": [
            SequenceFeature(start=2, end=8, name="f1", feature_type="LTR"),
        ],
        "s2": [],
    }
    canvas.set_msa_data(
        ["s1", "s2"],
        ["ACGTACGTAC", "ACGTACGTAC"],
        features=feats,
    )
    canvas.set_msa_residue_size(10)
    h_small = canvas.msa_residue_height()
    tops_small = list(canvas.row_heights())
    canvas.set_msa_residue_size(20)
    assert canvas.msa_residue_height() > h_small
    assert canvas.row_heights()[0] > tops_small[0]

    # Height keeps shrinking down to the minimum label size (2pt).
    canvas.set_msa_residue_size(8)
    h8 = canvas.msa_residue_height()
    canvas.set_msa_residue_size(2)
    h2 = canvas.msa_residue_height()
    assert h2 < h8
    assert h2 >= 1

    # Sequence distance opens a gap between rows (tops move apart).
    canvas.set_sequence_row_gap(0)
    top1_tight = canvas.row_tops()[1]
    canvas.set_sequence_row_gap(12)
    assert canvas.sequence_row_gap() == 12
    assert canvas.row_tops()[1] == top1_tight + 12

    with_track = canvas.row_heights()[0]
    canvas.set_show_annotation_track(False)
    assert canvas.row_heights()[0] < with_track
    assert canvas._hit_feature(0, canvas.bp_to_x(3), canvas._feature_band_y0(0) + 2) is None
    canvas.set_show_annotation_track(True)
    assert canvas._tracks["s1"]
    canvas.close()


def test_msa_name_label_tracks_dna_protein_size(qapp: QApplication, tmp_path: Path) -> None:
    from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel

    fasta = tmp_path / "aln.fa"
    fasta.write_text(">s1\nACGTACGT\n>s2\nACGTACGT\n", encoding="utf-8")
    panel = MultiSequenceOverviewPanel(fasta)
    panel.show()
    qapp.processEvents()
    panel.load_path(fasta, msa_mode=True)
    qapp.processEvents()
    panel._rail.bases_zoomed_in.setValue(12)
    qapp.processEvents()
    assert panel._names.label_point_size() == 13  # nucleotide size + 1
    panel._rail.bases_zoomed_in.setValue(4)
    qapp.processEvents()
    assert panel._names.label_point_size() == 5
    panel.close()


def test_add_features_keeps_msa_mode(qapp: QApplication, tmp_path: Path) -> None:
    from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel

    fasta = tmp_path / "aln.fa"
    fasta.write_text(
        ">s1\nACGTACGTACGTACGT\n>s2\nACGTACGTACGTACGT\n",
        encoding="utf-8",
    )
    panel = MultiSequenceOverviewPanel(fasta)
    panel.show()
    qapp.processEvents()
    panel.load_path(fasta, msa_mode=True)
    qapp.processEvents()
    assert panel.is_msa_mode()
    panel._bars.set_selected_rows([0, 1])
    panel._bars.set_region_selection(2, 5)
    ok = panel.add_features_and_save(
        {
            "s1": [SequenceFeature(start=2, end=5, name="x", feature_type="LTR")],
            "s2": [SequenceFeature(start=2, end=5, name="x", feature_type="LTR")],
        }
    )
    assert ok
    assert panel.is_msa_mode()
    assert len(panel._features["s1"]) >= 1
    panel.close()


def test_overview_ignores_nested_ltr_progress_counters(qapp: QApplication, tmp_path: Path) -> None:
    """Nested '0/1 LTR' must not reset the per-sequence 0/N progress bar."""
    from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel

    fasta = tmp_path / "m.fa"
    fasta.write_text(">a\nACGT\n>b\nACGT\n", encoding="utf-8")
    panel = MultiSequenceOverviewPanel(fasta)
    panel._ann_job_kind = "auto"
    panel._ann_progress_total = 41
    updates: list[tuple[int, str]] = []
    panel._update_window_job_progress = (  # type: ignore[method-assign]
        lambda completed, detail="": updates.append((int(completed), str(detail)))
    )
    statuses: list[str] = []
    panel._set_annotation_job_status = (  # type: ignore[method-assign]
        lambda text: statuses.append(str(text))
    )

    panel._on_batch_annotation_progress("0/41 lib.fa: sequence 1/41 (a)")
    panel._on_batch_annotation_progress("detail:0/1 LTR")
    panel._on_batch_annotation_progress("0/1 LTR")  # unprefixed nested — ignore
    panel._on_batch_annotation_progress("1/41 lib.fa: sequence 1/41 (a)")

    assert updates == [
        (0, "lib.fa: sequence 1/41 (a)"),
        (1, "lib.fa: sequence 1/41 (a)"),
    ]
    assert any("0/41" in s for s in statuses)
    assert any("1/41" in s for s in statuses)
    assert not any("0/1" in s for s in statuses)
    panel.close()


def test_run_te_annotation_filters_record_ids(tmp_path: Path) -> None:
    from repeatviewer.backends.auto_annotation import AutoAnnotationMode
    from repeatviewer.backends.te_annotation_batch import (
        TeAnnotationOptions,
        run_te_annotation_on_fasta,
    )

    fasta = tmp_path / "multi.fa"
    fasta.write_text(
        ">a\n" + "ACGT" * 40 + "\n>b\n" + "ACGT" * 40 + "\n>c\n" + "ACGT" * 40 + "\n",
        encoding="utf-8",
    )
    annotated: list[str] = []

    def _fake_annotate(record_fasta, options, **kwargs):
        del options, kwargs
        # Temp FASTA id is stripped of #suffix; track by reading the file header.
        text = Path(record_fasta).read_text(encoding="utf-8")
        header = text.splitlines()[0].lstrip(">")
        annotated.append(header)
        return []

    with patch(
        "repeatviewer.backends.te_annotation_batch._annotate_one_record_fasta",
        side_effect=_fake_annotate,
    ):
        run_te_annotation_on_fasta(
            fasta,
            TeAnnotationOptions(auto_modes=(AutoAnnotationMode.ORF,)),
            record_ids=["b"],
        )
    assert annotated == ["b"]


def test_name_click_selects_sequence_msa_click_clears(qapp: QApplication) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    canvas.set_msa_data(["a", "b"], ["ACGTACGT", "ACGTACGT"])
    canvas._select_row_click(1, Qt.KeyboardModifier.NoModifier)
    assert canvas.selected_rows() == [1]
    # Residue click (no drag) clears name-selected sequences.
    canvas._drag_origin = QPointF(10, 2)
    canvas._drag_kind = "region"
    canvas._press_feat = None
    canvas._rubber_active = False
    release = QMouseEvent(
        QMouseEvent.Type.MouseButtonRelease,
        QPointF(10, 2),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    canvas.mouseReleaseEvent(release)
    assert canvas.selected_rows() == []
    canvas.close()


def test_feature_click_sets_column_selection_and_shift_arrows(
    qapp: QApplication,
) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    feats = {
        "a": [SequenceFeature(start=4, end=6, name="f1", feature_type="LTR")],
        "b": [SequenceFeature(start=2, end=8, name="f2", feature_type="ORF")],
    }
    canvas.set_msa_data(["a", "b"], ["ACGTACGTAC", "ACGTACGTAC"], features=feats)
    # Select feature on row 0 (flat index 0) → columns 4–6.
    canvas._set_feature_selection({(0, 0)}, set_columns=True)
    assert canvas.region_selection() == (4, 6)
    assert canvas.selected_feature_column_span() == (4, 6)

    assert canvas.apply_shift_edge_expand("left")
    assert canvas.region_selection() == (1, 3)  # exclude left boundary
    assert canvas.apply_shift_edge_expand("left")
    assert canvas.region_selection() == (1, 6)  # include right boundary

    canvas._set_feature_selection({(0, 0)}, set_columns=True)
    assert canvas.apply_shift_edge_expand("right")
    assert canvas.region_selection() == (7, 10)  # exclude right boundary
    assert canvas.apply_shift_edge_expand("right")
    assert canvas.region_selection() == (4, 10)  # include left boundary

    # Multiple features → leftmost / rightmost.
    canvas._set_feature_selection({(0, 0), (1, 0)}, set_columns=True)
    assert canvas.selected_feature_column_span() == (2, 8)
    assert canvas.region_selection() == (2, 8)
    canvas.close()


def test_annotation_drag_selects_on_partial_overlap(qapp: QApplication) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    canvas.resize(800, 200)
    canvas.show()
    qapp.processEvents()
    feats = {
        "a": [SequenceFeature(start=2, end=5, name="f1", feature_type="LTR")],
        "b": [SequenceFeature(start=4, end=8, name="f2", feature_type="ORF")],
    }
    canvas.set_msa_data(["a", "b"], ["ACGTACGTAC", "ACGTACGTAC"], features=feats)
    canvas.set_msa_residue_size(12)
    canvas.set_pixels_per_base(20.0)
    qapp.processEvents()

    # Partial overlap of feature a still selects it.
    partial = QRect(
        int(canvas.bp_to_x(3.0)),
        int(canvas._feature_band_y0(0)),
        int(canvas.bp_to_x(4.5) - canvas.bp_to_x(3.0)),
        12,
    )
    assert (0, 0) in canvas._features_intersecting_rect(partial)

    full = QRect(
        int(canvas.bp_to_x(1.0)),
        int(canvas._feature_band_y0(0) - 1),
        int(canvas.bp_to_x(9.0) - canvas.bp_to_x(1.0)),
        int(canvas._feature_band_y0(1) + 14 - canvas._feature_band_y0(0)),
    )
    canvas._drag_origin = QPointF(full.left(), full.top())
    canvas._drag_kind = "features"
    canvas._rubber_active = True
    canvas._rubber = full
    release = QMouseEvent(
        QMouseEvent.Type.MouseButtonRelease,
        QPointF(full.right(), full.bottom()),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    canvas.mouseReleaseEvent(release)
    qapp.processEvents()
    assert (0, 0) in canvas.selected_features()
    assert (1, 0) in canvas.selected_features()
    assert canvas.region_selection() == (2, 8)
    assert canvas.regional_rows() is None
    assert canvas.selected_rows() == []  # column select clears name rows
    canvas.close()


def test_delete_selected_features_not_column_region(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Feature picks set a column span; Delete must remove only those features."""
    from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel

    fasta = tmp_path / "aln.fa"
    fasta.write_text(
        ">a\nACGTACGTAC\n>b\nACGTACGTAC\n>c\nACGTACGTAC\n",
        encoding="utf-8",
    )
    panel = MultiSequenceOverviewPanel(fasta)
    panel.show()
    qapp.processEvents()
    panel.load_path(fasta, msa_mode=True)
    qapp.processEvents()
    panel._features = {
        "a": [
            SequenceFeature(start=2, end=5, name="ORF1", feature_type="ORF"),
            SequenceFeature(start=6, end=8, name="ORF2", feature_type="ORF"),
        ],
        "b": [
            SequenceFeature(start=3, end=7, name="ORF3", feature_type="ORF"),
            SequenceFeature(start=4, end=6, name="Pfam", feature_type="domain"),
        ],
        "c": [
            SequenceFeature(start=2, end=8, name="keep", feature_type="ORF"),
        ],
    }
    panel._refresh_features_after_edit()
    # Selecting features also mirrors their column span onto region_selection.
    panel._bars._set_feature_selection({(0, 0), (1, 0), (1, 1)}, set_columns=True)
    assert panel._bars.region_selection() is not None
    assert len(panel._bars.selected_features()) == 3

    panel._on_delete_clicked()
    qapp.processEvents()

    assert [f.name for f in panel._features["a"]] == ["ORF2"]
    assert panel._features["b"] == []
    assert [f.name for f in panel._features["c"]] == ["keep"]
    assert panel._bars.msa_sequences() == ["ACGTACGTAC", "ACGTACGTAC", "ACGTACGTAC"]
    assert panel._bars.selected_features() == []
    panel.close()


def test_msa_pane_to_matrix_gap_does_not_jump_to_last_row(
    qapp: QApplication,
) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    canvas.resize(400, 400)
    canvas.set_msa_data(
        ["a", "b", "c", "d", "e"],
        ["ACGTACGTAC"] * 5,
    )
    canvas.set_sequence_row_gap(12)
    canvas.set_pixels_per_base(4.0)
    # y in the gap between row 1 and row 2 — must not clamp to last row.
    gap_y = canvas.row_tops()[1] + canvas.row_heights()[1] + 4
    assert canvas._row_at_y(gap_y) < 0
    _col, row = canvas._pane_to_matrix(canvas.bar_pad() + 8, gap_y)
    assert row in (1, 2)
    assert row != len(canvas.record_ids()) - 1


def test_msa_regional_select_snaps_box_and_keeps_matrix_rows(
    qapp: QApplication,
) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    canvas.resize(600, 200)
    canvas.show()
    qapp.processEvents()
    canvas.set_msa_data(
        ["a", "b", "c"],
        ["ACGTACGTAC", "ACGTACGTAC", "ACGTACGTAC"],
    )
    canvas.set_pixels_per_base(15.0)
    canvas._select_row_click(2, Qt.KeyboardModifier.NoModifier)
    assert canvas.selected_rows() == [2]
    qapp.processEvents()
    top0 = canvas.row_tops()[0]
    bot1 = canvas.row_tops()[1] + canvas.row_heights()[1] - 1
    x0 = canvas.bp_to_x(2.0)
    x1 = canvas.bp_to_x(6.0)
    # Partial-pixel rubber still snaps to full sequence bands of rows 0–1.
    canvas._set_live_region_rubber(
        QRect(int(x0), int(top0 + 2), int(x1 - x0), int((bot1 - top0) // 2))
    )
    assert canvas._temp_matrix_sel is not None
    c0, c1, r0, r1 = canvas._temp_matrix_sel
    snapped = canvas._matrix_selection_box(c0, c1, r0, r1)
    assert canvas._rubber == snapped
    canvas._drag_origin = QPointF(x0, top0)
    canvas._drag_kind = "region"
    canvas._rubber_active = True
    release = QMouseEvent(
        QMouseEvent.Type.MouseButtonRelease,
        QPointF(x1, bot1),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    canvas.mouseReleaseEvent(release)
    qapp.processEvents()
    assert canvas.is_regional_selection()
    assert canvas.region_selection() is not None
    assert canvas.regional_rows() == {0, 1}
    # Regional select clears prior whole-row selection (names use regional_rows).
    assert canvas.selected_rows() == []
    canvas.close()


def test_msa_regional_select_highlights_names_and_scale(
    qapp: QApplication, tmp_path: Path
) -> None:
    from repeatviewer.widgets.multi_sequence_overview import (
        MultiSequenceOverviewPanel,
    )

    fasta = tmp_path / "reg.fa"
    fasta.write_text(">a\nACGTACGTAC\n>b\nACGTACGTAC\n>c\nACGTACGTAC\n")
    panel = MultiSequenceOverviewPanel(fasta)
    panel.load_path(fasta, msa_mode=True)
    panel.show()
    qapp.processEvents()
    panel._bars.set_regional_selection(3, 7, {0, 1})
    qapp.processEvents()
    panel._sync_names_column()
    assert panel._names._selected_rows == {0, 1}
    assert panel._scale._region_sel == (3, 7)
    panel.close()


def test_scale_shows_column_selection_band(qapp: QApplication) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _ScaleRuler

    ruler = _ScaleRuler()
    ruler.set_scale(ppb=10.0, max_length=100, bar_pad=8)
    assert ruler._region_sel is None
    ruler.set_region_selection(10, 40)
    assert ruler._region_sel == (10, 40)
    ruler.set_region_selection(None)
    assert ruler._region_sel is None
    ruler.close()


def test_scale_drag_selects_columns_only(qapp: QApplication) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _ScaleRuler

    ruler = _ScaleRuler()
    ruler.resize(400, 28)
    ruler.set_scale(ppb=10.0, max_length=100, bar_pad=8)
    ruler.set_scroll_x(0)
    got: list[tuple[int, int]] = []
    ruler.columnRegionSelected.connect(lambda a, b: got.append((a, b)))
    ruler._drag_origin_x = 8.0 + 10.0 * 4  # near column 5
    ruler._rubber_active = True
    release = QMouseEvent(
        QMouseEvent.Type.MouseButtonRelease,
        QPointF(8.0 + 10.0 * 20, 10),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    ruler.mouseReleaseEvent(release)
    assert got
    assert got[0][0] < got[0][1]
    ruler.close()


def test_delete_removes_shift_delete_fills_gaps(
    qapp: QApplication, tmp_path: Path
) -> None:
    from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel

    fasta = tmp_path / "aln.fa"
    fasta.write_text(
        ">s1\nACGTACGTAC\n>s2\nACGTACGTAC\n",
        encoding="utf-8",
    )
    panel = MultiSequenceOverviewPanel(fasta)
    panel.show()
    qapp.processEvents()
    panel.load_path(fasta, msa_mode=True)
    qapp.processEvents()
    panel._bars.set_region_selection(3, 5)
    panel._features = {
        "s1": [
            SequenceFeature(start=2, end=8, name="overlap", feature_type="LTR"),
            SequenceFeature(start=9, end=10, name="keep", feature_type="ORF"),
        ],
        "s2": [],
    }
    # Plain Delete → remove columns; overlapping feature dropped + coords remapped.
    panel._delete_selected_msa_edit(hard_delete=True)
    qapp.processEvents()
    assert panel._bars.msa_sequences() == ["ACCGTAC", "ACCGTAC"]
    assert [f.name for f in panel._features["s1"]] == ["keep"]
    assert panel._features["s1"][0].start == 6  # 9 → 6 after deleting 3 cols
    assert panel._features["s1"][0].end == 7

    panel._restore_features()
    qapp.processEvents()
    assert panel._bars.msa_sequences() == ["ACGTACGTAC", "ACGTACGTAC"]

    # Shift+Delete → gaps; overlapping feature dropped.
    panel._bars.set_region_selection(3, 5)
    panel._features = {
        "s1": [
            SequenceFeature(start=2, end=8, name="overlap", feature_type="LTR"),
            SequenceFeature(start=9, end=10, name="keep", feature_type="ORF"),
        ],
        "s2": [],
    }
    panel._delete_selected_msa_edit(hard_delete=False)
    qapp.processEvents()
    assert panel._bars.msa_sequences() == ["AC---CGTAC", "AC---CGTAC"]
    assert [f.name for f in panel._features["s1"]] == ["keep"]
    assert panel._features["s1"][0].start == 9
    panel.close()


def test_scale_click_selects_column_double_click_zooms(
    qapp: QApplication,
) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _ScaleRuler

    ruler = _ScaleRuler()
    ruler.resize(400, 28)
    ruler.set_scale(ppb=10.0, max_length=100, bar_pad=8)
    ruler.set_scroll_x(20)
    zooms: list[tuple[float, bool]] = []
    selects: list[tuple[int, int]] = []
    ruler.zoomAtRequested.connect(lambda x, inn: zooms.append((x, inn)))
    ruler.columnRegionSelected.connect(lambda a, b: selects.append((a, b)))
    # Click without drag — select the column under the cursor.
    click_x = 8.0 + 9.5 * 10
    ruler._drag_origin_x = click_x
    ruler._rubber_active = False
    release = QMouseEvent(
        QMouseEvent.Type.MouseButtonRelease,
        QPointF(click_x, 10),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    ruler.mouseReleaseEvent(release)
    assert zooms == []
    expected = ruler._columns_from_x_range(click_x, click_x)
    assert selects == [expected]
    assert expected[0] == expected[1]

    # Double-click → zoom to max at content x (scale_x + scroll).
    zooms.clear()
    selects.clear()
    dbl = QMouseEvent(
        QMouseEvent.Type.MouseButtonDblClick,
        QPointF(click_x, 10),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    ruler.mouseDoubleClickEvent(dbl)
    assert selects == []
    assert len(zooms) == 1
    assert zooms[0][1] is True
    assert abs(zooms[0][0] - (click_x + 20)) < 1e-6

    # Shift+click → fit (zoom_in=False).
    zooms.clear()
    ruler._suppress_click_select = False
    ruler._drag_origin_x = click_x
    ruler._rubber_active = False
    release_out = QMouseEvent(
        QMouseEvent.Type.MouseButtonRelease,
        QPointF(click_x, 10),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.ShiftModifier,
    )
    ruler.mouseReleaseEvent(release_out)
    assert zooms == [(click_x + 20, False)]
    ruler.close()


def test_shift_arrow_expands_column_selection(qapp: QApplication) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    canvas.set_msa_data(["a", "b"], ["ACGTACGTAC", "ACGTACGTAC"])
    canvas.set_region_selection(4, 6)
    assert canvas.apply_shift_edge_expand("left")
    assert canvas.region_selection() == (1, 3)
    assert canvas.apply_shift_edge_expand("left")
    assert canvas.region_selection() == (1, 6)
    canvas.set_region_selection(4, 6)
    assert canvas.apply_shift_edge_expand("right")
    assert canvas.region_selection() == (7, 10)
    assert canvas.apply_shift_edge_expand("right")
    assert canvas.region_selection() == (4, 10)
    canvas.close()


def test_overview_consensus_button_runs_con_generater(
    qapp: QApplication, tmp_path: Path
) -> None:
    from repeatviewer.job_options import ConsensusOptions
    from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel

    fasta = tmp_path / "aln.fa"
    fasta.write_text(">a\nACGT\n>b\nACCT\n", encoding="utf-8")
    panel = MultiSequenceOverviewPanel(fasta)
    got: list[Path] = []
    panel.sequenceExtracted.connect(got.append)
    panel._run_consensus(
        ConsensusOptions(threshold=0.5, ambiguous="N", sequence_type="nucleotide")
    )
    assert len(got) == 1
    assert got[0].is_file()
    assert "_0.5_co.fa" in got[0].name
    text = got[0].read_text(encoding="utf-8")
    assert ">consensus" in text
    panel.close()


def test_overview_msa_clean_emits_additional_preview_tab(
    qapp: QApplication, tmp_path: Path
) -> None:
    from repeatviewer.job_options import MsaCleanOptions
    from repeatviewer.project import RepeatViewerProject
    from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel

    fasta = tmp_path / "aln.fa"
    # Column 2 is mostly gaps → gappy clean should drop it with default settings.
    fasta.write_text(
        ">a\nA-GT\n>b\nA-GT\n>c\nA-GT\n>d\nA-GT\n>e\nACGT\n",
        encoding="utf-8",
    )
    panel = MultiSequenceOverviewPanel(fasta)

    # Attach a parent so window() exposes _project for run_on_sequence.
    from PySide6.QtWidgets import QWidget

    host = QWidget()
    host._project = RepeatViewerProject()  # type: ignore[attr-defined]
    panel.setParent(host)

    extracted: list[Path] = []
    tabs: list[Path] = []
    panel.sequenceExtracted.connect(extracted.append)
    panel.openInAdditionalPreviewTab.connect(tabs.append)
    panel._run_msa_clean(
        "msa_clean_gappy",
        MsaCleanOptions(use_gappy_columns=True, use_divergence=False, use_gap_row=False),
    )
    assert len(extracted) == 1
    assert extracted == tabs
    assert extracted[0].is_file()
    assert extracted[0].name.endswith("_gap.fa") or "_gap" in extracted[0].name
    panel.close()
    host.close()


def test_scale_column_selection_shift_arrows_match_feature(
    qapp: QApplication,
) -> None:
    """Scale column select uses the same Shift+Left/Right two-step cycle."""
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    feats = {
        "a": [SequenceFeature(start=4, end=6, name="f1", feature_type="LTR")],
        "b": [],
    }
    canvas.set_msa_data(["a", "b"], ["ACGTACGTAC", "ACGTACGTAC"], features=feats)
    canvas._set_feature_selection({(0, 0)}, set_columns=True)
    assert canvas.selected_features() == [(0, 0)]

    # Pure column select (as from the scale) clears feature chrome.
    canvas.clear_feature_selection()
    canvas.set_region_selection(4, 6)
    assert canvas.selected_features() == []
    assert canvas.region_selection() == (4, 6)

    assert canvas.apply_shift_edge_expand("left")
    assert canvas.region_selection() == (1, 3)
    assert canvas.apply_shift_edge_expand("left")
    assert canvas.region_selection() == (1, 6)

    canvas.set_region_selection(4, 6)
    assert canvas.apply_shift_edge_expand("right")
    assert canvas.region_selection() == (7, 10)
    assert canvas.apply_shift_edge_expand("right")
    assert canvas.region_selection() == (4, 10)
    canvas.close()
