"""Data-dependent zoom limits for the multi-sequence overview."""

from __future__ import annotations

import pytest

from repeatviewer.widgets.multi_sequence_overview import (
    _MAX_PPB,
    _QT_MAX_WIDGET_SIZE,
    _BarsCanvas,
)


def _canvas_with_lengths(qapp, lengths: dict[str, int], width: int = 800) -> _BarsCanvas:
    canvas = _BarsCanvas()
    canvas.set_viewport_min_size(width, 300)
    canvas.set_data(list(lengths), lengths, {})
    return canvas


def test_fit_uses_longest_genome_sequence_without_fixed_minimum(qapp) -> None:
    canvas = _canvas_with_lengths(
        qapp,
        {"chr1": 50_000_000, "chr2": 35_000_000},
    )

    canvas.fit_to_bars_width(800)

    expected = (800 - 16) / 50_000_000
    assert canvas.pixels_per_base() == pytest.approx(expected)
    assert canvas.width() == 800
    canvas.deleteLater()


def test_fit_can_exceed_legacy_maximum_for_short_sequences(qapp) -> None:
    canvas = _canvas_with_lengths(qapp, {"short-a": 118, "short-b": 114})

    canvas.fit_to_bars_width(800)

    assert canvas.pixels_per_base() == pytest.approx((800 - 16) / 118)
    assert canvas.pixels_per_base() > _MAX_PPB
    assert canvas.width() == 800
    canvas.deleteLater()


def test_loading_genome_after_short_sequences_keeps_widget_size_safe(qapp) -> None:
    canvas = _canvas_with_lengths(qapp, {"short": 118})
    canvas.fit_to_bars_width(800)
    assert canvas.pixels_per_base() > _MAX_PPB

    canvas.set_data(["chr1"], {"chr1": 50_000_000}, {})

    assert canvas.width() <= _QT_MAX_WIDGET_SIZE
    fit_ppb, max_ppb = canvas.zoom_bounds(800)
    assert fit_ppb <= canvas.pixels_per_base() <= max_ppb
    canvas.deleteLater()
