"""Track-join / overlap-displace layout for annotation vertical moves."""

from __future__ import annotations

from repeatviewer.backends.sequence_features import SequenceFeature
from repeatviewer.widgets.sequence_map_viewer import reassign_tracks_after_row_join


def _feat(name: str, start: int, end: int, track: int = 0) -> SequenceFeature:
    return SequenceFeature(
        start=start,
        end=end,
        name=name,
        feature_type="domain",
        color="#888888",
        track=track,
    )


def test_moving_int_up_displaces_only_overlapping_hhcc() -> None:
    # Figure 1 layout: domains share track 1; INT is alone on track 2.
    features = [
        _feat("gag-pol", 1000, 6000, track=0),
        _feat("PROT", 2800, 3100, track=1),
        _feat("RT", 3200, 3800, track=1),
        _feat("RNH", 4000, 4360, track=1),
        _feat("HHCC", 4500, 4650, track=1),
        _feat("INT", 4500, 5600, track=2),
    ]
    int_idx = 5
    updated = reassign_tracks_after_row_join(
        features, moving={int_idx}, join_track=1, local_row=0
    )
    assert updated is not None
    by_name = {feat.name: feat for feat in updated}

    # INT shares a level with non-overlapping PROT / RT / RNH.
    shared = {by_name[name].track for name in ("PROT", "RT", "RNH", "INT")}
    assert len(shared) == 1

    # Only HHCC (overlap) is pushed to a lower track — not the whole row.
    assert by_name["HHCC"].track == by_name["INT"].track + 1
    assert by_name["PROT"].track == by_name["INT"].track
    assert by_name["RNH"].track == by_name["INT"].track
    assert by_name["gag-pol"].track < by_name["INT"].track


def test_join_without_overlap_keeps_single_row() -> None:
    features = [
        _feat("A", 100, 200, track=0),
        _feat("B", 300, 400, track=0),
        _feat("C", 500, 600, track=1),
    ]
    updated = reassign_tracks_after_row_join(
        features, moving={2}, join_track=0, local_row=0
    )
    assert updated is not None
    assert {feat.track for feat in updated} == {0}
