"""PBS motif path + LTR-pair selection for Feature detection."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.auto_annotation import AutoAnnotationMode
from repeatviewer.backends.pbs_detect import (
    _extract_ltrs,
    default_motif_fasta,
    pbs_motif_database_available,
)
from repeatviewer.backends.sequence_features import SequenceFeature


def test_bundled_pbs_motif_fasta_exists() -> None:
    path = default_motif_fasta()
    assert path.is_file()
    assert pbs_motif_database_available(path)
    assert path.name == "PBS_unique_motifs_Neumann2019.fasta"


def test_extract_ltrs_prefers_primary_pair() -> None:
    features = [
        SequenceFeature(1, 100, "LTR5", "LTR"),
        SequenceFeature(500, 600, "LTR3", "LTR"),
        SequenceFeature(50, 80, "LTR5_2", "LTR"),
        SequenceFeature(700, 800, "LTR3_2", "LTR"),
    ]
    left, right = _extract_ltrs(features)
    assert left.name == "LTR5"
    assert right.name == "LTR3"


def test_extract_ltrs_outermost_fallback() -> None:
    features = [
        SequenceFeature(10, 50, "left", "LTR"),
        SequenceFeature(900, 950, "right", "LTR"),
        SequenceFeature(100, 200, "mid", "ORF"),
    ]
    left, right = _extract_ltrs(features)
    assert (left.start, right.start) == (10, 900)


def test_extract_ltrs_requires_two() -> None:
    assert _extract_ltrs([SequenceFeature(1, 10, "LTR5", "LTR")]) is None


def test_pbs_mode_value() -> None:
    assert AutoAnnotationMode.PBS.value == "PBS"
