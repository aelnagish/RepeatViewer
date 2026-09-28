"""Per-feature source field and split-by-source packing."""

from __future__ import annotations

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
