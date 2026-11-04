"""LTR / TIR / TSD partner selection helpers."""

from __future__ import annotations

from repeatviewer.backends.sequence_features import (
    SequenceFeature,
    expand_terminal_pair_indices,
    parse_terminal_pair_key,
    partner_terminal_pair_key,
)


def test_parse_ltr_tir_tsd_pair_keys() -> None:
    assert parse_terminal_pair_key("LTR5") == ("LTR", "5", "")
    assert parse_terminal_pair_key("LTR3 (99.5%)") == ("LTR", "3", "")
    assert parse_terminal_pair_key("LTR5_2 (s)") == ("LTR", "5", "_2")
    assert parse_terminal_pair_key("TIR3_3") == ("TIR", "3", "_3")
    assert parse_terminal_pair_key("TSD5 ATAT") == ("TSD", "5", "ATAT")
    assert parse_terminal_pair_key("TSD3 atat") == ("TSD", "3", "ATAT")
    assert parse_terminal_pair_key("TSD_ATGCA") == ("TSD", "", "ATGCA")
    assert parse_terminal_pair_key("ORF1") is None


def test_partner_swaps_five_and_three() -> None:
    assert partner_terminal_pair_key(("LTR", "5", "")) == ("LTR", "3", "")
    assert partner_terminal_pair_key(("TIR", "3", "_2")) == ("TIR", "5", "_2")
    assert partner_terminal_pair_key(("TSD", "5", "ATAT")) == ("TSD", "3", "ATAT")
    assert partner_terminal_pair_key(("TSD", "", "ATGCA")) == ("TSD", "", "ATGCA")


def test_expand_selects_ltr_partner() -> None:
    features = [
        SequenceFeature(1, 10, "LTR5 (99.5%)", feature_type="LTR"),
        SequenceFeature(20, 30, "ORF1", feature_type="ORF"),
        SequenceFeature(90, 100, "LTR3 (99.5%)", feature_type="LTR"),
    ]
    assert expand_terminal_pair_indices(features, {0}) == {0, 2}
    assert expand_terminal_pair_indices(features, {2}) == {0, 2}
    assert expand_terminal_pair_indices(features, {1}) == {1}


def test_expand_nested_ltr_pairs_stay_separate() -> None:
    features = [
        SequenceFeature(1, 10, "LTR5", feature_type="LTR"),
        SequenceFeature(90, 100, "LTR3", feature_type="LTR"),
        SequenceFeature(15, 25, "LTR5_2", feature_type="LTR"),
        SequenceFeature(70, 80, "LTR3_2", feature_type="LTR"),
    ]
    assert expand_terminal_pair_indices(features, {0}) == {0, 1}
    assert expand_terminal_pair_indices(features, {2}) == {2, 3}


def test_expand_tsd5_tsd3_and_shared_name() -> None:
    features = [
        SequenceFeature(1, 4, "TSD5 ATAT", feature_type="TSD"),
        SequenceFeature(50, 53, "TSD3 ATAT", feature_type="TSD"),
        SequenceFeature(5, 9, "TSD_ATGCA", feature_type="TSD"),
        SequenceFeature(40, 44, "TSD_ATGCA", feature_type="TSD"),
    ]
    assert expand_terminal_pair_indices(features, {0}) == {0, 1}
    assert expand_terminal_pair_indices(features, {2}) == {2, 3}
