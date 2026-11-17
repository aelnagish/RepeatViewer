"""Dotplot axis ticks and equal X/Y scale helpers."""

from __future__ import annotations

import pytest
from PySide6.QtCore import QRectF

from repeatviewer.widgets.dotplot_view import (
    _axis_entry_at,
    _axis_tick_positions,
    _ensure_min_bp_rect,
    _format_bp,
    _nice_tick_step,
    _non_overlapping_tick_labels,
    _scene_to_bp,
    _scene_to_local_bp,
)


def test_dotplot_rail_range_boxes_align_hyphen(qapp) -> None:
    """X/Y range fields use ' - ' and share a start column so hyphens line up."""
    from repeatviewer.widgets.dotplot_rail import DotplotRail

    rail = DotplotRail()
    assert rail.sel_x_label.text() == "X"
    assert rail.sel_y_label.text() == "Y"
    assert rail.sel_x.start_text() == ""
    assert not rail.sel_x._sep.isVisibleTo(rail.sel_x)
    rail.set_selection_ranges(x=(256, 1698), y=(5517, 6553))
    assert rail.sel_x._sep.text() == " - "
    assert rail.sel_y._sep.text() == " - "
    assert rail.sel_x.start_text() == "256"
    assert rail.sel_y.start_text() == "5,517"
    # Longer start wins the shared column width.
    assert rail.sel_x._start.width() == rail.sel_y._start.width()
    assert rail.sel_x._start.width() >= rail.sel_y._start.fontMetrics().horizontalAdvance(
        "5,517"
    )
    assert rail.view_x_btn.isEnabled()
    assert rail.view_y_btn.isEnabled()
    # Selection block sits after the parameter host (bottom of the panel).
    page = rail.sel_row.parentWidget()
    page_layout = page.layout()
    assert page_layout.indexOf(rail.sel_row) > page_layout.indexOf(rail._params_host)
    # Show reverse sits below Line width in the form.
    form = rail.line_width.parentWidget().layout().itemAt(0).layout()
    widgets = []
    for i in range(form.rowCount()):
        field = form.itemAt(i, form.ItemRole.FieldRole)
        if field is not None and field.widget() is not None:
            widgets.append(field.widget())
    assert widgets.index(rail.line_width) < widgets.index(rail.reverse)
    rail.deleteLater()


def test_zoom_keeps_selection_visible(qapp) -> None:
    """Zooming with a selection anchors on it and keeps it in the viewport."""
    from PySide6.QtCore import QRectF

    from repeatviewer.backends.dotplot import (
        STRAND_FORWARD,
        DotplotResult,
        DotplotSegment,
    )
    from repeatviewer.widgets.dotplot_view import DotplotCanvas

    canvas = DotplotCanvas()
    canvas.resize(900, 700)
    canvas.show()
    qapp.processEvents()
    canvas.show_result(
        DotplotResult(
            sequence_id="t",
            sequence_length=10_000,
            word_size=6,
            segments=[DotplotSegment(0, 0, 100, STRAND_FORWARD)],
            include_reverse=True,
        )
    )
    qapp.processEvents()
    canvas.reset_zoom()
    qapp.processEvents()

    canvas._selection = QRectF(7000.0, 7000.0, 200.0, 200.0)
    # Zoom in several steps — selection must remain inside the visible scene.
    for _ in range(8):
        canvas._apply_zoom(1.5, canvas.viewport().rect().topLeft())
        qapp.processEvents()
        vis = canvas.visible_scene_rect()
        sel = canvas._selection
        assert sel is not None
        assert vis.intersects(sel), f"selection left the view: vis={vis} sel={sel}"
        # Centre of selection should stay on-screen after zoom-about-selection.
        assert vis.contains(sel.center())
    canvas.close()


def test_preserve_view_word_size_does_not_drift(qapp) -> None:
    """Zoomed word-size recomputes must keep scroll/hit mapping stable.

    Restoring via centerOn() rounded through integer scrollbars and crept ~1 px
    per change, so hits walked away from the crosshair after spinning Word size.
    """
    from PySide6.QtCore import QPointF, QRectF

    from repeatviewer.backends.dotplot import (
        STRAND_FORWARD,
        DotplotResult,
        DotplotSegment,
    )
    from repeatviewer.widgets.dotplot_view import DotplotCanvas

    canvas = DotplotCanvas()
    canvas.resize(900, 700)
    canvas.show()
    qapp.processEvents()
    seg = DotplotSegment(500, 800, 100, STRAND_FORWARD)
    canvas.show_result(
        DotplotResult(
            sequence_id="t",
            sequence_length=2_000,
            word_size=8,
            segments=[seg],
            include_reverse=False,
        )
    )
    qapp.processEvents()
    canvas.reset_zoom()
    qapp.processEvents()
    canvas._selection = QRectF(450.0, 750.0, 200.0, 200.0)
    assert canvas.zoom_to_selection()
    qapp.processEvents()

    sx, sy = 550.0, 850.0
    before_map = canvas.mapFromScene(QPointF(sx, sy))
    before_scroll = (
        canvas.horizontalScrollBar().value(),
        canvas.verticalScrollBar().value(),
    )
    for ws in (10, 15, 20, 25, 30, 35, 40, 45):
        canvas.show_result(
            DotplotResult(
                sequence_id="t",
                sequence_length=2_000,
                word_size=ws,
                segments=[seg],
                include_reverse=False,
            ),
            preserve_view=True,
        )
        qapp.processEvents()

    after_map = canvas.mapFromScene(QPointF(sx, sy))
    after_scroll = (
        canvas.horizontalScrollBar().value(),
        canvas.verticalScrollBar().value(),
    )
    assert after_scroll == before_scroll, (
        f"scroll drifted across word-size recomputes: "
        f"{before_scroll} -> {after_scroll}"
    )
    assert after_map == before_map, (
        f"scene→view mapping drifted: {before_map} -> {after_map}"
    )
    canvas.close()


def test_ensure_min_bp_rect_expands_sub_base_drag() -> None:
    """At high zoom a few-pixel drag is < 1 bp; still become a 1×1 selection."""
    bounds = QRectF(0.0, 0.0, 10_000.0, 10_000.0)
    tiny = QRectF(3768.2, 3343.1, 0.3, 0.4)
    grown = _ensure_min_bp_rect(tiny, bounds)
    assert grown.width() >= 1.0
    assert grown.height() >= 1.0
    assert bounds.contains(grown)


def test_ensure_min_bp_rect_clamps_near_edge() -> None:
    bounds = QRectF(0.0, 0.0, 100.0, 100.0)
    near = QRectF(99.7, 0.1, 0.2, 0.2)
    grown = _ensure_min_bp_rect(near, bounds)
    assert grown.width() >= 1.0
    assert grown.height() >= 1.0
    assert grown.right() <= bounds.right() + 1e-9
    assert grown.left() >= bounds.left() - 1e-9


def test_nice_tick_step_for_6kb_uses_1k_not_5k() -> None:
    """Fully visible ~6.6 kb plot should get ~1 kb steps, not a single 5k tick."""
    step = _nice_tick_step(6614, target_ticks=10)
    assert step == 1000
    ticks = list(range(0, 6615, step))
    assert 0 in ticks
    assert 5000 in ticks
    assert len(ticks) >= 6


def test_nice_tick_step_short_sequence() -> None:
    step = _nice_tick_step(182, target_ticks=10)
    assert step <= 50
    assert 182 / step >= 3


def test_format_bp_compact_labels() -> None:
    assert _format_bp(1) == "1"
    assert _format_bp(500) == "500"
    assert _format_bp(1000) == "1,000"
    assert _format_bp(5000) == "5,000"
    assert _format_bp(1_000_000) == "1,000,000"
    assert _format_bp(1500) == "1,500"


def test_scene_to_bp_is_1_based() -> None:
    assert _scene_to_bp(0, length=100, offset=0) == 1
    assert _scene_to_bp(50, length=100, offset=0) == 50
    assert _scene_to_bp(100, length=100, offset=0) == 100
    assert _scene_to_bp(0, length=50, offset=10) == 11
    assert _scene_to_bp(50, length=50, offset=10) == 60


def test_shared_tick_step_ignores_letterbox_span(qapp) -> None:
    """Empty space beside a fitted square must not inflate the tick interval."""
    from repeatviewer.widgets.dotplot_view import DotplotCanvas

    canvas = DotplotCanvas()
    canvas._length = 6614
    # Pretend the viewport mapped to a scene wider than the sequence (letterbox).
    canvas.visible_scene_rect = lambda: QRectF(-100, 0, 20000, 6614)  # type: ignore[method-assign]
    step = canvas.shared_tick_step()
    assert step == 1000
    canvas.deleteLater()


def test_axis_end_tick_replaces_crowding_major() -> None:
    """Length 810 with step 100 must keep 810 and drop 800 (no 800/810 overlap)."""
    ticks = _axis_tick_positions(lo=0.0, hi=810.0, length=810, step=100)
    assert 810 in ticks
    assert 800 not in ticks
    assert 0 in ticks
    assert 700 in ticks


def test_axis_exact_step_end_keeps_single_end() -> None:
    ticks = _axis_tick_positions(lo=0.0, hi=800.0, length=800, step=100)
    assert ticks.count(800) == 1
    assert 800 in ticks


def test_multi_axis_local_bp_restarts_at_one() -> None:
    """Each concatenated block reports 1-based positions from its own start."""
    from repeatviewer.backends.paf_dotplot import AxisEntry

    entries = [
        AxisEntry("a", 0, 100),
        AxisEntry("b", 100, 100),
    ]
    assert _axis_entry_at(entries, 0).name == "a"
    assert _axis_entry_at(entries, 99).name == "a"
    assert _axis_entry_at(entries, 100).name == "b"
    assert _axis_entry_at(entries, 200).name == "b"
    assert _scene_to_local_bp(0, entries) == 1
    assert _scene_to_local_bp(50, entries) == 50
    assert _scene_to_local_bp(100, entries) == 1
    assert _scene_to_local_bp(150, entries) == 50
    assert _scene_to_local_bp(200, entries) == 100
    # Scene 0 and scene 100 (starts of a and b) both label as local 1.
    assert _format_bp(_scene_to_local_bp(0, entries)) == "1"
    assert _format_bp(_scene_to_local_bp(100, entries)) == "1"


def test_multi_axes_widen_rulers_and_display_local_xy(qapp) -> None:
    from repeatviewer.backends.paf_dotplot import AxisEntry
    from repeatviewer.widgets.dotplot_view import (
        DotplotCanvas,
        _RULER_LEFT,
        _RULER_LEFT_MULTI,
        _RULER_TOP,
        _RULER_TOP_MULTI,
    )

    canvas = DotplotCanvas()
    assert canvas.ruler_top() == _RULER_TOP
    assert canvas.ruler_left() == _RULER_LEFT
    # Left margin must fit comma-formatted multi-Mb positions (was clipping
    # "38,110,309" down to "8,110,309" on the Y ruler).
    assert _RULER_LEFT >= 70
    assert _RULER_LEFT_MULTI >= _RULER_LEFT
    canvas.set_axes(
        [AxisEntry("a", 0, 100), AxisEntry("b", 100, 100)],
        [AxisEntry("a", 0, 100), AxisEntry("b", 100, 100)],
    )
    canvas._length = 200
    canvas._height = 200
    assert canvas.ruler_top() == _RULER_TOP_MULTI
    assert canvas.ruler_left() == _RULER_LEFT_MULTI
    assert canvas._display_xy(0, 100) == (1, 1)
    assert canvas._display_xy(50, 150) == (50, 50)
    canvas.set_axes([], [])
    assert canvas.ruler_top() == _RULER_TOP
    canvas.deleteLater()


def test_multi_dot_identity_spans_per_sequence_block(qapp) -> None:
    """All-sequences / multi_dot: thick identity only inside each self block."""
    from repeatviewer.backends.dotplot import DotplotResult
    from repeatviewer.backends.paf_dotplot import AxisEntry
    from repeatviewer.widgets.dotplot_view import DotplotCanvas

    canvas = DotplotCanvas()
    result = DotplotResult(
        sequence_id="2 x 2 sequences", sequence_length=200, word_size=0
    )
    canvas.set_axes(
        [AxisEntry("a", 0, 100), AxisEntry("b", 100, 100)],
        [AxisEntry("a", 0, 100), AxisEntry("b", 100, 100)],
    )
    spans, mirror = canvas._identity_spans_for(result)
    assert spans == [
        (0.0, 0.0, 100.0, 100.0),
        (100.0, 100.0, 200.0, 200.0),
    ]
    assert mirror is False
    canvas.set_axes([], [])
    cross, cross_mirror = canvas._identity_spans_for(
        DotplotResult(sequence_id="a vs b", sequence_length=100, word_size=0)
    )
    assert cross == [] and cross_mirror is False
    self_spans, self_mirror = canvas._identity_spans_for(
        DotplotResult(sequence_id="a", sequence_length=50, word_size=0)
    )
    assert self_spans == [(0.0, 0.0, 50.0, 50.0)] and self_mirror is True
    # Pairwise two-sequence labels: names on each axis; cross has no identity.
    canvas.set_axes(
        [AxisEntry("query", 0, 276)],
        [AxisEntry("target", 0, 755)],
    )
    cross_pair, cross_pair_mirror = canvas._identity_spans_for(
        DotplotResult(
            sequence_id="query vs target",
            sequence_length=276,
            word_size=0,
            target_id="target",
            target_length=755,
        )
    )
    assert cross_pair == [] and cross_pair_mirror is False
    canvas.set_axes(
        [AxisEntry("same", 0, 100)],
        [AxisEntry("same", 0, 100)],
    )
    self_pair, self_pair_mirror = canvas._identity_spans_for(
        DotplotResult(sequence_id="same", sequence_length=100, word_size=0)
    )
    assert self_pair == [(0.0, 0.0, 100.0, 100.0)] and self_pair_mirror is True
    canvas.set_axes(
        [AxisEntry("a", 0, 100), AxisEntry("b", 100, 100)],
        [AxisEntry("a", 0, 100), AxisEntry("b", 100, 100)],
    )
    canvas._length = 200
    canvas._height = 200
    seps = canvas._sequence_separator_segments(200.0, 200.0)
    assert len(seps) == 2
    assert seps[0][0].x() == 100.0 and seps[0][1].y() == 200.0
    assert seps[1][0].y() == 100.0 and seps[1][1].x() == 200.0
    canvas.deleteLater()


def test_non_overlapping_prefers_sequence_end() -> None:
    ticks = [0, 100, 200, 800, 810]
    positions = {t: float(t) for t in ticks}

    def extent(pos: float, label: str) -> tuple[float, float]:
        # Pretend each digit is 10 px so 800 and 810 collide hard.
        return pos, pos + 10.0 * len(label)

    shown = _non_overlapping_tick_labels(
        ticks,
        length=810,
        offset=0,
        positions=positions,
        label_extent=extent,
        pad=4.0,
    )
    assert 810 in shown
    assert 800 not in shown
    assert 0 in shown


def test_export_line_width_scales_with_dpi() -> None:
    """Display Line width must thicken on high-DPI PDF/PNG pages."""
    from repeatviewer.widgets.dotplot_view import (
        _export_line_width_px,
        _line_width_px,
    )

    assert _line_width_px(22) == pytest.approx(2.2)
    # SVG / 72 dpi: same as screen mapping, floored at 1 px.
    assert _export_line_width_px(8, 72) == pytest.approx(1.0)
    assert _export_line_width_px(22, 72) == pytest.approx(2.2)
    # Default PDF dpi: thicker strokes so step changes stay visible on A4.
    assert _export_line_width_px(8, 400) == pytest.approx(8 * 0.1 * 400 / 72.0)
    assert _export_line_width_px(22, 400) == pytest.approx(22 * 0.1 * 400 / 72.0)
    assert _export_line_width_px(22, 400) > _export_line_width_px(8, 400)


def test_pdf_export_applies_line_width_and_clears_cache(qapp, tmp_path) -> None:
    """PDF used to skip export-pen setup and keep DeviceCoordinateCache."""
    from PySide6.QtWidgets import QGraphicsItem

    from repeatviewer.backends.dotplot import (
        STRAND_FORWARD,
        DotplotResult,
        DotplotSegment,
    )
    from repeatviewer.widgets.dotplot_view import (
        DotplotCanvas,
        _export_line_width_px,
    )

    canvas = DotplotCanvas()
    canvas.resize(400, 400)
    canvas.show()
    qapp.processEvents()
    canvas.show_result(
        DotplotResult(
            sequence_id="seq",
            sequence_length=100,
            word_size=10,
            segments=[DotplotSegment(10, 10, 40, STRAND_FORWARD)],
            include_reverse=False,
        )
    )
    qapp.processEvents()
    canvas.set_line_width_step(22)
    assert canvas._item is not None
    canvas._item.setCacheMode(QGraphicsItem.CacheMode.DeviceCoordinateCache)

    widths: list[float | None] = []
    original = canvas._item.set_export_pen_px

    def _record(width: float | None) -> None:
        widths.append(width)
        original(width)

    canvas._item.set_export_pen_px = _record  # type: ignore[method-assign]
    out = tmp_path / "dotplot.pdf"
    canvas.export_to_path(out, dpi=400, title="t")
    assert out.is_file() and out.stat().st_size > 0
    expected = _export_line_width_px(22, 400)
    assert any(
        w is not None and abs(float(w) - expected) < 1e-6 for w in widths
    )
    assert widths[-1] is None
    assert canvas._item.cacheMode() == QGraphicsItem.CacheMode.DeviceCoordinateCache
    canvas.deleteLater()
