"""Find in the Tree view: name / trait search and stepping (no Qt)."""

from __future__ import annotations

from repeatviewer.backends import tree_find as tf

TIPS = ["Knob_chr1_100", "Knob_chr2_200", "CroSat2_x", "LTR1_chr1"]
TRAITS = {
    "Knob_chr1_100": {"chromosome": "chr1", "array": "A1"},
    "Knob_chr2_200": {"chromosome": "chr2", "array": "A7"},
    "LTR1_chr1": {"chromosome": "chr1", "lineage": "Ale"},
}


def test_name_match_keeps_tree_order_and_ignores_case():
    assert tf.find_tips(TIPS, "knob") == ["Knob_chr1_100", "Knob_chr2_200"]
    assert tf.find_tips(TIPS, "knob", case_sensitive=True) == []


def test_several_words_match_any():
    assert tf.find_tips(TIPS, "crosat2, ltr1") == ["CroSat2_x", "LTR1_chr1"]


def test_trait_values_can_be_switched_off():
    assert tf.find_tips(TIPS, "Ale", traits=TRAITS) == ["LTR1_chr1"]
    assert tf.find_tips(TIPS, "Ale", traits=TRAITS, in_traits=False) == []
    assert tf.find_tips(TIPS, "A7", traits=TRAITS, in_names=False) == ["Knob_chr2_200"]


def test_empty_query_finds_nothing():
    assert tf.find_tips(TIPS, "  ,; ") == []


def test_step_index_wraps_both_ways():
    assert tf.step_index(-1, 3, 1) == 0
    assert tf.step_index(-1, 3, -1) == 2
    assert tf.step_index(2, 3, 1) == 0
    assert tf.step_index(0, 3, -1) == 2
    assert tf.step_index(-1, 0, 1) == -1
