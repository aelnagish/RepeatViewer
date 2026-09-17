"""Dotplot axis ticks and equal X/Y scale helpers."""

from __future__ import annotations

from PySide6.QtCore import QRectF

from repeatviewer.widgets.dotplot_view import (
    _axis_tick_positions,
    _format_bp,
    _nice_tick_step,
    _non_overlapping_tick_labels,
    _scene_to_bp,
)


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
