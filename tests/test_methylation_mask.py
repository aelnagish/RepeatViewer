"""CG/CHG masking sensitivity run (pure)."""

from __future__ import annotations

import io

from Bio import Phylo

from repeatviewer.backends import methylation_mask as mm


def test_cg_marks_both_bases():
    assert mm.methylation_columns(["ACGT"], mm.MASK_CG) == [1, 2]


def test_chg_forward_and_reverse_strand():
    # CAG: the C is forward CHG; the G is the C of CTG on the other strand.
    assert mm.methylation_columns(["TCAGT"], mm.MASK_CG_CHG) == [1, 3]
    # CCG: C-C-G -> first C is CHG (H = C), second C and G are CG.
    assert mm.methylation_columns(["CCG"], mm.MASK_CG_CHG) == [0, 1, 2]
    # CHG only counts when asked for
    assert mm.methylation_columns(["TCAGT"], mm.MASK_CG) == []


def test_chh_is_not_masked():
    assert mm.methylation_columns(["TCATT", "ACTAA"], mm.MASK_CG_CHG) == []


def test_context_is_read_across_gaps():
    # C and G separated by a gap in the alignment are still a CG in the sequence.
    assert mm.methylation_columns(["AC--GT"], mm.MASK_CG) == [1, 4]


def test_any_sequence_flags_the_column():
    aligned = ["ACGT", "ATGT"]           # TG in the second = C→T at a CG
    cols = mm.methylation_columns(aligned, mm.MASK_CG)
    assert cols == [1, 2]
    assert mm.remove_columns(aligned, cols) == ["AT", "AT"]


def test_support_thresholds():
    assert mm.is_well_supported("92/97", "iqtree")
    assert not mm.is_well_supported("92/90", "iqtree")
    assert not mm.is_well_supported("75/99", "iqtree")
    assert mm.is_well_supported("96", "iqtree")
    assert mm.is_well_supported("70", "fastme")
    assert mm.is_well_supported("0.95", "fasttree")
    assert not mm.is_well_supported("", "iqtree")


def _tree(text):
    return Phylo.read(io.StringIO(text), "newick")


def test_identical_trees():
    t = "((A,B)95/99,(C,D)90/98,(E,F)85/96);"
    cmp = mm.compare_trees(_tree(t), _tree(t), "iqtree")
    assert cmp.rf == 0 and cmp.normalized_rf == 0.0
    assert cmp.supported_full == 3 and not cmp.lost_supported


def test_lost_supported_grouping_is_reported():
    full = _tree("(((A,B)95/99,C)50/60,(D,E)90/98,F);")
    masked = _tree("(((A,C),B),(D,E),F);")
    cmp = mm.compare_trees(full, masked, "iqtree")
    assert cmp.rf > 0
    assert cmp.lost_supported == [["A", "B"]]
    assert cmp.kept_supported == 1
    text = mm.comparison_sentence(cmp, 12, 100, mm.MASK_CG_CHG)
    assert "12 of 100 columns, 12%" in text and "1 of 2 well-supported" in text


def test_sentence_says_differences_are_only_weak_branches():
    full = _tree("(((A,B)95/99,C)50/60,((D,E)40/50,F)90/98,(G,H)99/100);")
    masked = _tree("(((A,B),C),((D,F),E),(G,H));")
    cmp = mm.compare_trees(full, masked, "iqtree")
    assert cmp.rf > 0 and not cmp.lost_supported
    text = mm.comparison_sentence(cmp, 5, 50, mm.MASK_CG)
    assert "differences concern only branches without strong support" in text
