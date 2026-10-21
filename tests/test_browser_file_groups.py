"""Tests for Files-panel source+output family grouping."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.browser_file_groups import file_family_key, family_band_flags


def test_file_family_key_groups_blast_outputs() -> None:
    assert file_family_key("reads.fa") == "reads.fa"
    assert file_family_key("reads.fa_bn.txt") == "reads.fa"
    assert file_family_key("reads.fa_bn.txt.fa") == "reads.fa"
    assert file_family_key("reads.fa_bp.txt") == "reads.fa"
    assert file_family_key("reads.fa_obn.html") == "reads.fa"
    assert file_family_key("reads.fa_obtn.txt") == "reads.fa"


def test_file_family_key_groups_align_and_ops() -> None:
    assert file_family_key("reads.fa_mafft.fa") == "reads.fa"
    assert file_family_key("reads.fa_mafft_1.fa") == "reads.fa"
    assert file_family_key("reads.fa_rc.fa") == "reads.fa"
    assert file_family_key("reads.fa_orfs.fa") == "reads.fa"
    assert file_family_key("reads.fa_translated.faa") == "reads.fa"
    assert file_family_key("reads.fa_TEAid.html") == "reads.fa"
    assert file_family_key("reads.fa_pTEAid.html") == "reads.fa"
    assert file_family_key("reads.fa_pTEAid_2.html") == "reads.fa"
    assert file_family_key("reads.fa_flexidot.png") == "reads.fa"
    assert file_family_key("reads.fa_multi_flexidot.html") == "reads.fa"
    assert file_family_key("reads.fa_multi_flexidot_2.png") == "reads.fa"
    assert file_family_key("reads.fa_self_dot.paf") == "reads.fa"
    assert file_family_key("reads.fa_multi_dot.paf") == "reads.fa"
    assert file_family_key("reads.fa_self_dot_2.paf") == "reads.fa"
    assert file_family_key("reads.fa_hmmer.txt") == "reads.fa"
    assert file_family_key("reads.fa_hmmer2.txt") == "reads.fa"
    assert file_family_key("reads.fa_L1000_R1000.fa") == "reads.fa"
    assert file_family_key("reads.fa_0.8_co.fa") == "reads.fa"


def test_family_band_flags_alternate_by_group() -> None:
    names = [
        "a.fa",
        "a.fa_bn.txt",
        "a.fa_bn.txt.fa",
        "b.fa",
        "b.fa_bn.txt",
        "c.fa",
    ]
    assert family_band_flags(names) == [False, False, False, True, True, False]


def test_row_zebra_colors_differ_from_family() -> None:
    from repeatviewer.browser_file_groups import (
        FAMILY_BAND_LIGHT,
        ROW_ZEBRA_LIGHT,
    )

    assert ROW_ZEBRA_LIGHT != FAMILY_BAND_LIGHT


def test_unrelated_files_keep_own_keys() -> None:
    assert file_family_key("notes.txt") == "notes.txt"
    assert file_family_key("genome.fasta") == "genome.fasta"
