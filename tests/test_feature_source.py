"""Per-feature source field and split-by-source packing."""

from __future__ import annotations

import pytest

from repeatviewer.backends.sequence_features import (
    SequenceFeature,
    feature_display_pack_key,
    ordered_feature_pack_keys,
    pack_feature_display_rows,
)


def _feat(
    name: str,
    start: int,
    end: int,
    *,
    track: int = 0,
    source: str = "",
) -> SequenceFeature:
    return SequenceFeature(
        start=start,
        end=end,
        name=name,
        feature_type="orf",
        color="#27ae60",
        track=track,
        source=source,
    )


def test_source_round_trip_and_missing_defaults_empty() -> None:
    feat = _feat("ORF1", 10, 200, source="genes.gff3")
    payload = feat.to_dict()
    assert payload["source"] == "genes.gff3"
    restored = SequenceFeature.from_dict(payload)
    assert restored.source == "genes.gff3"

    legacy = SequenceFeature.from_dict(
        {
            "start": 1,
            "end": 5,
            "name": "x",
            "feature_type": "others",
            "color": "#888888",
        }
    )
    assert legacy.source == ""


def test_transfer_meta_round_trip_and_hover_text() -> None:
    from repeatviewer.backends.sequence_features import (
        feature_has_hover_meta,
        feature_hover_detail_text,
    )

    feat = SequenceFeature(
        start=10,
        end=100,
        name="LTR5 (93.8%)",
        feature_type="LTR",
        source="Transfer (Tork0002192.fa)",
        pident=37.3,
        coverage=0.806,
        subject_id="Tork0002192",
        source_start=120,
        source_end=210,
    )
    payload = feat.to_dict()
    assert payload["coverage"] == 0.806
    assert payload["pident"] == 37.3
    assert payload["source_start"] == 120
    assert payload["source_end"] == 210
    restored = SequenceFeature.from_dict(payload)
    assert restored.coverage == 0.806
    assert restored.pident == 37.3
    assert restored.source_start == 120
    assert restored.source_end == 210
    # Percent-style coverage in JSON is normalized to 0–1.
    assert SequenceFeature.from_dict({**payload, "coverage": 80.6}).coverage == pytest.approx(
        0.806
    )

    tip = feature_hover_detail_text(feat, primary=feat.name)
    assert "LTR5 (93.8%)" in tip
    assert "Source: Transfer (Tork0002192.fa)" in tip
    assert "Source position: 120–210" in tip
    assert "Identity: 37.3%" in tip
    assert "Coverage: 80.6%" in tip
    assert feature_has_hover_meta(feat) is True
    assert feature_has_hover_meta(SequenceFeature(start=1, end=2, name="x")) is False


def test_transfer_hover_omits_source_position_for_non_transfer() -> None:
    from repeatviewer.backends.sequence_features import feature_hover_detail_text

    tip = feature_hover_detail_text(
        SequenceFeature(
            start=1,
            end=10,
            name="ORF1",
            source="genes.gff3",
            source_start=5,
            source_end=15,
        ),
        primary="ORF1",
    )
    assert "Source position:" not in tip


def test_pack_rows_split_by_source_separates_same_track() -> None:
    features = [
        _feat("a", 1, 50, source=""),
        _feat("b", 1, 50, source="genes.gff3"),
        _feat("c", 60, 80, source="repeatviewer_internal"),
    ]
    rows_split, count_split = pack_feature_display_rows(features, split_by_source=True)
    assert count_split == 3
    assert len(set(rows_split)) == 3

    rows_combined, count_combined = pack_feature_display_rows(
        features, split_by_source=False
    )
    # Same track + overlapping a/b pack into 2 rows within one track band.
    assert count_combined == 2
    assert rows_combined[0] != rows_combined[1]
    assert rows_combined[2] in {0, 1}


def test_pack_key_order_empty_source_first() -> None:
    features = [
        _feat("z", 1, 2, track=1, source="b.gff"),
        _feat("y", 1, 2, track=0, source=""),
        _feat("x", 1, 2, track=0, source="a.gff"),
    ]
    keys = ordered_feature_pack_keys(features, split_by_source=True)
    assert keys[0] == ("", 0)
    assert keys[1] == ("a.gff", 0)
    assert keys[2] == ("b.gff", 1)
    assert feature_display_pack_key(features[1], split_by_source=False) == ("", 0)
