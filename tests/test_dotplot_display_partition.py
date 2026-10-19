"""Display partitioning: keep LTR-scale diagonals above short-word noise."""

from __future__ import annotations

from repeatviewer.backends.dotplot import (
    STRAND_FORWARD,
    DotplotSegment,
    partition_segments_for_display,
    structural_length_floor,
)


def test_structural_floor_grows_with_word_and_length() -> None:
    assert structural_length_floor(4, 11_000) >= 12
    assert structural_length_floor(10, 500) >= 30


def test_partition_keeps_long_runs_and_thins_noise() -> None:
    # Many short noise hits + a few longer LTR-like diagonals.
    segments = [
        DotplotSegment(x=i, y=i + 10, length=4, strand=STRAND_FORWARD)
        for i in range(5_000)
    ]
    segments.append(DotplotSegment(x=100, y=500, length=200, strand=STRAND_FORWARD))
    segments.append(DotplotSegment(x=300, y=800, length=120, strand=STRAND_FORWARD))

    noise, structural = partition_segments_for_display(
        segments,
        word_size=4,
        sequence_length=11_000,
        max_noise=500,
    )
    assert len(noise) <= 500
    assert len(structural) == 2
    assert [s.length for s in structural] == [120, 200]
    assert all(s.length < structural_length_floor(4, 11_000) for s in noise)
