"""Multi-seq / MSA Annotation parity: sizes, selection, tracks, row height."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest
from PySide6.QtCore import QPointF, QRect, Qt
from PySide6.QtGui import QMouseEvent
from PySide6.QtWidgets import QApplication, QMessageBox

from repeatviewer.backends.sequence_features import SequenceFeature
from repeatviewer.widgets.preview_controls import (
    OVERVIEW_TOOLBAR_HEIGHT,
    ZOOM_SLIDER_HEIGHT,
    AnnotationZoomBar,
    _OVERVIEW_POSITION_SPIN_WIDTH,
    _OVERVIEW_SELECT_FIT_HEIGHT,
    _ZOOM_BTN_SIDE,
)


def test_overview_zoom_bar_matches_single_seq_sizes(qapp: QApplication) -> None:
    bar = AnnotationZoomBar()
    bar.set_overview_mode(True)
    bar.show()
    qapp.processEvents()
    assert bar.height() == OVERVIEW_TOOLBAR_HEIGHT == 28
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
    assert not panel._toolbar_shim.isHidden()
    assert bar._actions_btn.isHidden()
    # Toolbar sits above tracks only (left edge aligned with the sequence map).
    from PySide6.QtCore import QPoint
    from PySide6.QtWidgets import QPushButton

    slot_x = panel._toolbar_slot.mapTo(panel, QPoint(0, 0)).x()
    shim_w = panel._toolbar_shim.width()
    assert slot_x >= shim_w - 1
    select = next(
        w
        for w in bar._overview_host.findChildren(QPushButton)
        if w.text() == "Select"
    )
    assert select.size().height() == _OVERVIEW_SELECT_FIT_HEIGHT
    assert select.size().width() <= 80
    panel._unbind_header_zoom_bar()
    assert panel._toolbar_shim.isHidden()
    panel.close()
    bar.close()


def test_overview_rail_keeps_add_transfer_and_track_toggle(qapp: QApplication) -> None:
    from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail

    rail = SequenceAnnotationRail()
    rail.show()
    qapp.processEvents()
    rail.set_overview_display_mode(True)
    assert not rail.show_annotation_track.isHidden()
    assert rail.show_annotation_track.isChecked()
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


def test_msa_drag_keeps_regional_box_without_row_highlight(
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
    canvas._drag_origin = QPointF(x0, top0)
    canvas._drag_kind = "region"
    canvas._rubber_active = True
    canvas._rubber = QRect(int(x0), int(top0), int(x1 - x0), int(bot1 - top0))
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
    # Regional select clears prior name-column row selection.
    assert canvas.selected_rows() == []
    canvas.close()


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


def test_scale_click_selects_single_column(qapp: QApplication) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _ScaleRuler

    ruler = _ScaleRuler()
    ruler.resize(400, 28)
    ruler.set_scale(ppb=10.0, max_length=100, bar_pad=8)
    ruler.set_scroll_x(0)
    got: list[tuple[int, int]] = []
    ruler.columnRegionSelected.connect(lambda a, b: got.append((a, b)))
    # Click without drag near column 10 (x ≈ 8 + 9.5*10).
    ruler._drag_origin_x = 8.0 + 9.5 * 10
    ruler._rubber_active = False
    release = QMouseEvent(
        QMouseEvent.Type.MouseButtonRelease,
        QPointF(8.0 + 9.5 * 10, 10),
        Qt.MouseButton.LeftButton,
        Qt.MouseButton.LeftButton,
        Qt.KeyboardModifier.NoModifier,
    )
    ruler.mouseReleaseEvent(release)
    assert got == [(10, 10)]
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
