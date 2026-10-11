"""orfipy ORF detection respects selected reading frames."""

from __future__ import annotations

import pytest

from repeatviewer.backends.auto_annotation_settings import OrfSettings, normalize_orf_frames
from repeatviewer.backends.orf_detect import detect_orf_features, orfipy_available

pytestmark = pytest.mark.skipif(not orfipy_available(), reason="orfipy not installed")


def test_normalize_orf_frames_defaults_and_order() -> None:
    assert normalize_orf_frames(()) == (1, 2, 3, -1, -2, -3)
    assert normalize_orf_frames((-2, 1, 9, 3)) == (1, 3, -2)


def test_detect_orf_features_filters_frames() -> None:
    # Forward ORF in frame +1 and a reverse ORF in frame -1 (from orfipy).
    seq = (
        "ATGGGGGTAGACATTCAGATGAATATATATTAGATGTTTTTTTAG"
        + "N" * 30
        + "CTAAAAAAAACATCTAATATATATTCATCTGAATGTCTACCCCCAT"
    )
    all_frames = detect_orf_features(
        "s",
        seq,
        settings=OrfSettings(min_length=9, max_length=10_000, frames=(1, 2, 3, -1, -2, -3)),
    )
    fwd_only = detect_orf_features(
        "s",
        seq,
        settings=OrfSettings(min_length=9, max_length=10_000, frames=(1, 2, 3)),
    )
    rev_only = detect_orf_features(
        "s",
        seq,
        settings=OrfSettings(min_length=9, max_length=10_000, frames=(-1, -2, -3)),
    )
    assert all_frames.features
    assert all(f.strand == "+" for f in fwd_only.features)
    assert all(f.strand == "-" for f in rev_only.features)
    assert len(fwd_only.features) + len(rev_only.features) == len(all_frames.features)
    assert "orfipy" in all_frames.message.lower()
