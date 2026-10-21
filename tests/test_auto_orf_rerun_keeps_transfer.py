"""ORF Feature detection re-run must keep Transfer ORFs."""

from __future__ import annotations

from repeatviewer.backends.auto_annotation import (
    AutoAnnotationMode,
    _is_auto_generated_orf,
    _remove_features_for_auto_rerun,
)
from repeatviewer.backends.sequence_features import (
    REPEATVIEWER_INTERNAL_SOURCE,
    SequenceFeature,
)


def test_transfer_orfs_are_not_auto_generated() -> None:
    transferred = SequenceFeature(
        start=10,
        end=100,
        name="ORF1",
        feature_type="ORF",
        source="Transfer (Tork0002192.fa)",
    )
    detected = SequenceFeature(
        start=20,
        end=200,
        name="ORF1",
        feature_type="ORF",
        source=REPEATVIEWER_INTERNAL_SOURCE,
    )
    legacy = SequenceFeature(
        start=30,
        end=90,
        name="ORF2",
        feature_type="ORF",
        source="",
    )
    imported = SequenceFeature(
        start=40,
        end=80,
        name="ORF3",
        feature_type="ORF",
        source="genes.gff3",
    )
    assert _is_auto_generated_orf(transferred) is False
    assert _is_auto_generated_orf(detected) is True
    assert _is_auto_generated_orf(legacy) is True
    assert _is_auto_generated_orf(imported) is False


def test_orf_rerun_removes_only_feature_detection_orfs() -> None:
    features = [
        SequenceFeature(
            start=1,
            end=50,
            name="LTR5",
            feature_type="LTR",
            source="Transfer (src.fa)",
        ),
        SequenceFeature(
            start=100,
            end=400,
            name="ORF1",
            feature_type="ORF",
            source="Transfer (src.fa)",
        ),
        SequenceFeature(
            start=120,
            end=500,
            name="ORF1",
            feature_type="ORF",
            source=REPEATVIEWER_INTERNAL_SOURCE,
        ),
        SequenceFeature(
            start=200,
            end=300,
            name="ORF2",
            feature_type="ORF",
            source="",
        ),
    ]
    kept = _remove_features_for_auto_rerun(features, AutoAnnotationMode.ORF)
    names_sources = {(f.name, f.source) for f in kept}
    assert ("LTR5", "Transfer (src.fa)") in names_sources
    assert ("ORF1", "Transfer (src.fa)") in names_sources
    assert ("ORF1", REPEATVIEWER_INTERNAL_SOURCE) not in names_sources
    assert ("ORF2", "") not in names_sources
