"""Gappy-column filter helpers (Clean preview + remove_gaps)."""

from __future__ import annotations

from pathlib import Path

from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from repeatviewer.backends.remove_gap import (
    coalesce_column_indices_to_ranges,
    gappy_columns_to_keep,
    gappy_columns_to_remove,
    remove_gaps_with_similarity_check,
)


def _five_row_msa() -> list[str]:
    # Col0 A (keep); col1 mostly gaps (remove); col2 G (keep);
    # col3 all gaps (remove); col4 A (keep).
    return [
        "A-G-A",
        "A-G-A",
        "A-G-A",
        "A-G-A",
        "AAG-A",
    ]


def test_gappy_columns_to_remove_marks_high_gap_columns() -> None:
    seqs = _five_row_msa()
    # gap_threshold 0.5 drops columns with >50% gaps.
    removed = gappy_columns_to_remove(
        seqs,
        gap_threshold=0.5,
        simi_check_gap_thr=0.4,
        similarity_thr=0.7,
        min_nucleotide=2,
    )
    assert 1 in removed
    assert 3 in removed
    kept = gappy_columns_to_keep(
        seqs,
        gap_threshold=0.5,
        simi_check_gap_thr=0.4,
        similarity_thr=0.7,
        min_nucleotide=2,
    )
    assert 0 in kept
    assert 2 in kept
    assert 4 in kept
    assert set(removed).isdisjoint(kept)
    assert len(kept) + len(removed) == len(seqs[0])


def test_gappy_columns_to_remove_empty_when_too_few_sequences() -> None:
    assert gappy_columns_to_remove(["AAA", "AAA", "AAA"]) == []


def test_coalesce_column_indices_to_ranges() -> None:
    assert coalesce_column_indices_to_ranges([0, 1, 2, 5, 7, 8]) == [
        (1, 3),
        (6, 6),
        (8, 9),
    ]
    assert coalesce_column_indices_to_ranges([]) == []


def test_remove_gaps_file_matches_helper(tmp_path: Path) -> None:
    seqs = [
        "AAAA",
        "AAAA",
        "AAAA",
        "AAAA",
        "A-AA",
    ]
    # Col1 is 20% gaps — with gap_threshold 0.1 it is removed immediately.
    records = [
        SeqRecord(Seq(s), id=f"s{i}", description="") for i, s in enumerate(seqs)
    ]
    src = tmp_path / "in.fa"
    out = tmp_path / "out.fa"
    AlignIO.write(MultipleSeqAlignment(records), src, "fasta")
    remove_gaps_with_similarity_check(
        str(src),
        str(out),
        gap_threshold=0.1,
        simi_check_gap_thr=0.05,
        similarity_thr=0.7,
        min_nucleotide=2,
    )
    filtered = [str(rec.seq) for rec in AlignIO.read(out, "fasta")]
    removed = gappy_columns_to_remove(
        seqs,
        gap_threshold=0.1,
        simi_check_gap_thr=0.05,
        similarity_thr=0.7,
        min_nucleotide=2,
    )
    keep = [i for i in range(4) if i not in set(removed)]
    expected = ["".join(s[i] for i in keep) for s in seqs]
    assert filtered == expected
