"""Large / multi-contig PAF handling (no UI freeze on pair listing)."""

from __future__ import annotations

from pathlib import Path

import pytest

from repeatviewer.backends.dotplot import DotplotSegment, STRAND_FORWARD
from repeatviewer.backends.paf_dotplot import (
    MAX_ALL_SEQUENCES,
    MAX_PAIR_COMBO_ITEMS,
    MAX_PLOT_ALIGNMENTS,
    PafAlignments,
    PafError,
    keep_longest_segments,
    parse_paf,
    result_for_all,
    result_for_pair,
)


def _record(
    qname: str,
    tname: str,
    *,
    qlen: int = 1000,
    tlen: int = 1000,
) -> tuple:
    return (qname, qlen, 0, 100, "+", tname, tlen, 0, 100)


def test_pairs_scales_with_observed_pairs_not_name_product() -> None:
    """Hundreds of names must not trigger an N² name×name walk."""
    n = 800
    aln = PafAlignments(
        records=[_record(f"s{i}", f"s{i}", qlen=n - i, tlen=n - i) for i in range(n)],
        query_lengths={f"s{i}": n - i for i in range(n)},
        target_lengths={f"s{i}": n - i for i in range(n)},
    )
    pairs = aln.pairs()
    assert len(pairs) == n
    assert all(q == t for q, t in pairs)


def test_pairs_for_combo_prefers_longest_self_and_caps() -> None:
    records = []
    qlen: dict[str, int] = {}
    tlen: dict[str, int] = {}
    for i in range(MAX_PAIR_COMBO_ITEMS + 50):
        name = f"ctg{i}"
        length = 1000 + i
        records.append(_record(name, name, qlen=length, tlen=length))
        qlen[name] = length
        tlen[name] = length
    # One cross pair that should rank below long self-pairs.
    records.append(_record("ctg0", "ctg1", qlen=1000, tlen=1001))
    aln = PafAlignments(records=records, query_lengths=qlen, target_lengths=tlen)
    menu, total = aln.pairs_for_combo()
    assert total == MAX_PAIR_COMBO_ITEMS + 51
    assert len(menu) == MAX_PAIR_COMBO_ITEMS
    assert menu[0][0] == menu[0][1]
    assert menu[0][0] == f"ctg{MAX_PAIR_COMBO_ITEMS + 49}"


def test_supports_all_sequences_view_limit() -> None:
    names = [f"s{i}" for i in range(MAX_ALL_SEQUENCES + 1)]
    aln = PafAlignments(
        records=[_record(n, n) for n in names],
        query_lengths={n: 100 for n in names},
        target_lengths={n: 100 for n in names},
    )
    assert not aln.supports_all_sequences_view()
    with pytest.raises(PafError, match="All-sequences view"):
        result_for_all(aln)


def test_parse_paf_keeps_longest_not_first(tmp_path: Path) -> None:
    """Global cap must prefer late long hits over early short ones."""
    path = tmp_path / "long_tail.paf"
    lines = [
        f"early\t1000\t0\t10\t+\tearly\t1000\t0\t10\t10\t10\t60\n"
        for _ in range(20)
    ]
    lines.append("late\t5000\t0\t4000\t+\tlate\t5000\t0\t4000\t4000\t4000\t60\n")
    path.write_text("".join(lines), encoding="utf-8")
    aln = parse_paf(path, max_records=5)
    assert aln.truncated
    assert any(r[0] == "late" for r in aln.records)
    assert aln.pair_spans.get(("late", "late"), 0) >= 4000


def test_parse_paf_pair_rescans_one_pair(tmp_path: Path) -> None:
    from repeatviewer.backends.paf_dotplot import parse_paf_pair

    path = tmp_path / "multi.paf"
    lines = [
        "a\t1000\t0\t100\t+\tb\t2000\t0\t100\t100\t100\t60\n",
        "a\t1000\t10\t200\t+\tb\t2000\t50\t240\t190\t190\t60\n",
        "c\t800\t0\t50\t+\td\t900\t0\t50\t50\t50\t60\n",
    ]
    path.write_text("".join(lines), encoding="utf-8")
    pair = parse_paf_pair(path, "a", "b")
    assert len(pair.records) == 2
    assert all(r[0] == "a" and r[5] == "b" for r in pair.records)


def test_parse_paf_cancel(tmp_path: Path) -> None:
    path = tmp_path / "tiny.paf"
    lines = [
        f"q{i}\t100\t0\t50\t+\tt{i}\t100\t0\t50\t50\t50\t60\n" for i in range(20)
    ]
    path.write_text("".join(lines), encoding="utf-8")
    calls = {"n": 0}

    def cancel_after_few() -> bool:
        calls["n"] += 1
        return calls["n"] > 5

    with pytest.raises(PafError, match="Cancelled"):
        parse_paf(path, is_cancelled=cancel_after_few)


def test_is_single_pair_without_full_sort() -> None:
    aln = PafAlignments(
        records=[_record("a", "a"), _record("a", "a")],
        query_lengths={"a": 100},
        target_lengths={"a": 100},
    )
    assert aln.is_single_pair()
    aln.records.append(_record("a", "b"))
    aln.target_lengths["b"] = 100
    assert not aln.is_single_pair()


def test_result_for_pair_uses_rectangular_axes() -> None:
    """Unequal query/target lengths must not pad both axes to max(q, t)."""
    aln = PafAlignments(
        records=[_record("q", "t", qlen=2843, tlen=5372)],
        query_lengths={"q": 2843},
        target_lengths={"t": 5372},
    )
    result = result_for_pair(aln, "q", "t")
    assert result.sequence_length == 2843
    assert result.effective_target_length == 5372
    assert result.target_id == "t"


def test_result_for_all_mirrors_cross_pair_triangle() -> None:
    """Pairwise PAFs store A→B only; All-sequences should fill B→A too."""
    aln = PafAlignments(
        records=[
            ("a", 100, 10, 40, "+", "b", 100, 20, 50),
        ],
        query_lengths={"a": 100},
        target_lengths={"b": 100},
    )
    # lengths dict needs both names for a 2×2 layout
    aln.query_lengths["b"] = 100
    aln.target_lengths["a"] = 100
    result, x_entries, y_entries = result_for_all(aln)
    assert [e.name for e in x_entries] == [e.name for e in y_entries]
    # One stored hit + its transpose across y=x.
    assert result.segment_count == 2
    xs = sorted((s.x, s.y) for s in result.segments)
    assert xs[0] != xs[1]


def test_keep_longest_segments_matches_alnplot_budget() -> None:
    """FastGA ALNplot keeps the longest N hits so structure stays visible."""
    segs = [
        DotplotSegment(x=i, y=i, length=(i % 50) + 1, strand=STRAND_FORWARD)
        for i in range(500)
    ]
    kept, total, capped = keep_longest_segments(segs, limit=100)
    assert total == 500 and capped and len(kept) == 100
    assert kept[0].length <= kept[-1].length
    assert kept[-1].length == max(s.length for s in segs)
    assert MAX_PLOT_ALIGNMENTS == 100_000


def test_result_for_all_two_sequences_skips_self_hits() -> None:
    aln = PafAlignments(
        records=[
            ("a", 100, 0, 50, "+", "a", 100, 0, 50),
            ("a", 100, 10, 40, "+", "b", 100, 20, 50),
            ("b", 100, 0, 40, "+", "b", 100, 0, 40),
        ],
        query_lengths={"a": 100, "b": 100},
        target_lengths={"a": 100, "b": 100},
    )
    result, _x, _y = result_for_all(aln)
    # Self rows dropped; one cross hit + mirror.
    assert result.segment_count == 2
