"""TE-Aid consensus axis tick spacing stays readable for long sequences."""

from __future__ import annotations

from repeatviewer.backends.teaid_plotter import (
    _consensus_axis_tick_kwargs,
    _consensus_position_ticks,
    _format_consensus_tick_label,
    _nice_consensus_tick_step,
)


def test_long_consensus_uses_sparse_nice_steps() -> None:
    # ~27 kb TE previously used step=1000 → ~27 overlapping labels.
    ticks = _consensus_position_ticks(27_000)
    assert ticks[0] == 1
    assert ticks[-1] < 27_000
    assert len(ticks) <= 12
    assert _nice_consensus_tick_step(27_000) >= 2000


def test_short_consensus_still_has_useful_ticks() -> None:
    ticks = _consensus_position_ticks(500)
    assert ticks[0] == 1
    assert len(ticks) >= 3
    assert len(ticks) <= 12


def test_tick_kwargs_include_compact_labels() -> None:
    kwargs = _consensus_axis_tick_kwargs(20_000)
    assert kwargs["tickmode"] == "array"
    assert len(kwargs["tickvals"]) == len(kwargs["ticktext"])
    assert _format_consensus_tick_label(15_000) == "15k"
    assert _format_consensus_tick_label(500) == "500"
    assert "15k" in kwargs["ticktext"] or any(
        t.endswith("k") for t in kwargs["ticktext"]
    )
