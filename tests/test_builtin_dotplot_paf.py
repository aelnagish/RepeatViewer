"""Built-in self-dotplot → PAF export (reopenable in Preview)."""

from __future__ import annotations

from pathlib import Path

import pytest

from repeatviewer.backends.dotplot import (
    DotplotError,
    compute_cross_dotplot,
    compute_self_dotplot,
    segments_to_paf_lines,
    self_dotplot_paf_path,
    write_pairwise_dotplot_paf,
    write_self_dotplot_paf,
)
from repeatviewer.backends.paf_dotplot import parse_paf


def test_self_dotplot_paf_path_naming(tmp_path: Path) -> None:
    fa = tmp_path / "seq.fa"
    fa.write_text(">s\nACGT\n", encoding="utf-8")
    assert self_dotplot_paf_path(fa).name == "seq.fa_self_dot.paf"
    fasta = tmp_path / "seq.fasta"
    fasta.write_text(">s\nACGT\n", encoding="utf-8")
    assert self_dotplot_paf_path(fasta).name == "seq.fasta_self_dot.paf"
    gz = tmp_path / "seq.fa.gz"
    gz.write_bytes(b"")
    assert self_dotplot_paf_path(gz).name == "seq.fa.gz_self_dot.paf"


def test_segments_to_paf_round_trip() -> None:
    seq = ("ACGT" * 40) + ("TGCA" * 40) + ("ACGT" * 40)
    result = compute_self_dotplot(seq, sequence_id="demo", word_size=8)
    lines = segments_to_paf_lines(result)
    assert lines
    assert all(line.startswith("demo\t") for line in lines)


def test_write_self_dotplot_paf_reopens(tmp_path: Path) -> None:
    fasta = tmp_path / "te.fa"
    # Repeat motif so self-matches exist at modest word size.
    unit = "ACGTACGGTT"
    fasta.write_text(f">te1\n{unit * 30}\n", encoding="utf-8")
    out = write_self_dotplot_paf(fasta, word_size=8, include_reverse=True)
    assert out.name == "te.fa_self_dot.paf"
    assert out.is_file()
    text = out.read_text(encoding="utf-8")
    assert text.startswith("# RepeatViewer-dotplot")
    assert "word_size=8" in text.splitlines()[0]
    aln = parse_paf(out)
    assert aln.records
    assert aln.query_lengths.get("te1", 0) == 300
    assert aln.word_size == 8
    assert aln.include_reverse is True


def test_write_self_dotplot_paf_empty_raises(tmp_path: Path) -> None:
    fasta = tmp_path / "empty.fa"
    fasta.write_text(">x\n\n", encoding="utf-8")
    with pytest.raises(DotplotError):
        write_self_dotplot_paf(fasta, word_size=8)


def test_compute_cross_dotplot_paf_lines() -> None:
    unit = "ACGTACGGTT"
    a = unit * 20
    b = ("TTTT" + unit * 18 + "GGGG")
    result = compute_cross_dotplot(
        a, b, query_id="q1", target_id="t1", word_size=8
    )
    lines = segments_to_paf_lines(result)
    assert lines
    assert all(line.startswith("q1\t") for line in lines)
    assert all("\tt1\t" in line for line in lines)


def test_write_pairwise_dotplot_paf_two_files(tmp_path: Path) -> None:
    unit = "ACGTACGGTT"
    a = tmp_path / "a.fa"
    b = tmp_path / "b.fa"
    a.write_text(f">a1\n{unit * 25}\n", encoding="utf-8")
    b.write_text(f">b1\n{unit * 25}\n", encoding="utf-8")
    out = write_pairwise_dotplot_paf([a, b], word_size=8, include_reverse=False)
    assert out.name == "a.fa_multi_dot.paf"
    assert out.is_file()
    aln = parse_paf(out)
    assert aln.records
    assert aln.word_size == 8
    assert aln.include_reverse is False
    names = {r[0] for r in aln.records} | {r[5] for r in aln.records}
    assert "a1" in names and "b1" in names
    # Exactly two sequences: cross-only (no A×A / B×B).
    self_ids = {r[0] for r in aln.records if r[0] == r[5]}
    assert not self_ids


def test_write_pairwise_dotplot_paf_needs_two_sequences(tmp_path: Path) -> None:
    fasta = tmp_path / "solo.fa"
    fasta.write_text(">only\nACGTACGTACGTACGTACGT\n", encoding="utf-8")
    with pytest.raises(DotplotError, match="two or more"):
        write_pairwise_dotplot_paf([fasta], word_size=8)


def test_write_pairwise_dotplot_paf_three_files_includes_self(tmp_path: Path) -> None:
    unit = "ACGTACGGTT"
    paths = []
    for name in ("a", "b", "c"):
        path = tmp_path / f"{name}.fa"
        path.write_text(f">{name}1\n{unit * 25}\n", encoding="utf-8")
        paths.append(path)
    out = write_pairwise_dotplot_paf(paths, word_size=8, include_reverse=False)
    aln = parse_paf(out)
    self_ids = {r[0] for r in aln.records if r[0] == r[5]}
    assert self_ids == {"a1", "b1", "c1"}


def test_result_for_all_shared_axis_order_makes_straight_diagonal(
    tmp_path: Path,
) -> None:
    """Pairwise PAF: same sequences on X and Y at the same offsets.

    Combination pairs only list a name as query *or* target; without a unified
    axis the identity blocks stagger. Shared layout keeps one straight diagonal.
    """
    from repeatviewer.backends.paf_dotplot import parse_paf, result_for_all

    unit = "ACGTACGGTT"
    paths = []
    for name in ("s1", "s2", "s3"):
        path = tmp_path / f"{name}.fa"
        path.write_text(f">{name}\n{unit * 20}\n", encoding="utf-8")
        paths.append(path)
    out = write_pairwise_dotplot_paf(paths, word_size=8, include_reverse=False)
    aln = parse_paf(out)
    # Asymmetric membership like combinations: not every name is both q and t.
    assert set(aln.query_lengths) != set(aln.target_lengths) or len(aln.query_lengths) >= 2
    _result, x_entries, y_entries = result_for_all(aln)
    assert [e.name for e in x_entries] == [e.name for e in y_entries]
    assert [e.offset for e in x_entries] == [e.offset for e in y_entries]
    assert [e.length for e in x_entries] == [e.length for e in y_entries]
    # Identity spans (viewer) would be continuous: each block starts where the
    # previous ended on both axes.
    for i, entry in enumerate(x_entries):
        assert entry.offset == sum(e.length for e in x_entries[:i])
