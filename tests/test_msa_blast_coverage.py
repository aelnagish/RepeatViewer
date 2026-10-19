"""Unit tests for TE-Aid-style MSA BLAST coverage helpers."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.msa_blast_coverage import (
    blast_hit_coverage,
    gapped_coverage_from_blast_table,
    interpolate_coverage_across_gaps,
    parse_teaid_blast_qspans,
    project_coverage_to_gapped,
)
from repeatviewer.backends.teaid_plotter import _alignment_dumb_consensus


def test_blast_hit_coverage_overlaps() -> None:
    # Two hits overlapping on positions 2–3 (1-based) → depths [1,2,2,1].
    depth = blast_hit_coverage(4, [(1, 3), (2, 4)])
    assert depth == [1, 2, 2, 1]
    assert blast_hit_coverage(0, [(1, 2)]) == []
    assert blast_hit_coverage(3, []) == [0, 0, 0]


def test_interpolate_coverage_across_gaps() -> None:
    # A=10, gap, gap, T=40 → gaps become 20 and 30.
    seq = "A--T"
    cov = [10, 0, 0, 40]
    out = interpolate_coverage_across_gaps(seq, cov)
    assert abs(out[0] - 10.0) < 1e-9
    assert abs(out[1] - 20.0) < 1e-9
    assert abs(out[2] - 30.0) < 1e-9
    assert abs(out[3] - 40.0) < 1e-9
    # Leading / trailing gaps hold nearest nucleotide.
    out2 = interpolate_coverage_across_gaps("-A-", [0, 5, 0])
    assert abs(out2[0] - 5.0) < 1e-9
    assert abs(out2[1] - 5.0) < 1e-9
    assert abs(out2[2] - 5.0) < 1e-9


def test_project_coverage_to_gapped_gaps_are_zero() -> None:
    # Ungapped ACGT with depths 1,2,3,4 → gapped A-C-GT keeps gaps at 0.
    gapped = "A-C-GT"
    depth = [1, 2, 3, 4]
    out = project_coverage_to_gapped(gapped, depth)
    assert out == [1, 0, 2, 0, 3, 4]
    assert project_coverage_to_gapped("AC", [9, 8]) == [9, 8]
    assert project_coverage_to_gapped("---", []) == [0, 0, 0]


def test_parse_teaid_blast_and_gapped_from_table(tmp_path: Path) -> None:
    table = tmp_path / "hits.csv"
    table.write_text(
        "V1,V2,V3,V4,V5,V6,V7,V8,V9,V10,V11,V12,V13\n"
        "q,s,100,2,0,0,1,2,1,2,plus,0,10\n"
        "q,s,100,2,0,0,2,3,1,2,plus,0,10\n",
        encoding="utf-8",
    )
    spans = parse_teaid_blast_qspans(table)
    assert spans == [(1, 2), (2, 3)]
    # Ungapped AT: depths [1, 2] → gapped A-T- → [1, 0, 2, 0].
    cov = gapped_coverage_from_blast_table("A-T-", table)
    assert cov == [1, 0, 2, 0]


def test_consensus_threshold_zero_matches_msa_width() -> None:
    # Threshold 0: majority letter wins even with sparse support; width preserved.
    seqs = ["A-CGT", "ATCT-", "A-C-T"]
    cons = _alignment_dumb_consensus(seqs, threshold=0.0, ambiguous="N")
    assert len(cons) == 5
    assert "-" not in cons
    assert cons[0] == "A"


def test_coverage_axis_tick_helpers() -> None:
    from repeatviewer.widgets.multi_sequence_overview import (
        _coverage_axis_tick_values,
        _coverage_log_y,
        _format_coverage_tick,
    )

    assert _format_coverage_tick(0) == "0"
    assert _format_coverage_tick(42) == "42"
    assert _format_coverage_tick(1500) == "1.5k"
    assert _coverage_axis_tick_values(0) == [0.0]
    ticks = _coverage_axis_tick_values(100)
    assert ticks[0] == 0.0
    assert ticks[-1] == 100.0
    assert len(ticks) >= 2
    # Top of the band maps peak → near y=0; bottom maps 0 → bottom_y.
    assert abs(_coverage_log_y(0, peak=100, bottom_y=40, usable=38) - 40.0) < 1e-6
    assert abs(_coverage_log_y(100, peak=100, bottom_y=40, usable=38) - 2.0) < 1e-6


def test_blast_coverage_json_roundtrip_and_splice(tmp_path: Path) -> None:
    from repeatviewer.backends.msa_blast_coverage import (
        blast_coverage_sidecar_path,
        delete_blast_coverage_json,
        load_blast_coverage_json,
        save_blast_coverage_json,
        splice_coverage_columns,
    )

    fasta = tmp_path / "aln.fa"
    fasta.write_text(">a\nACGT\n", encoding="utf-8")
    cov = [1, 2, 3, 4, 5, 6]
    path = save_blast_coverage_json(fasta, cov)
    assert path == blast_coverage_sidecar_path(fasta)
    assert path.is_file()
    assert load_blast_coverage_json(fasta) == cov
    # Delete columns 2–4 (1-based) → keep 1 and 5–6.
    assert splice_coverage_columns(cov, 2, 4) == [1, 5, 6]
    delete_blast_coverage_json(fasta)
    assert load_blast_coverage_json(fasta) is None


def test_coverage_viewport_peak_safe_when_scrolled_past_end(qapp) -> None:
    """Horizontal scroll past the MSA end must not IndexError on peak calc."""
    from repeatviewer.widgets.multi_sequence_overview import _MsaTopCoverageProfile

    profile = _MsaTopCoverageProfile()
    profile.resize(400, 48)
    profile.set_coverage([1, 2, 3, 4, 5])
    profile.set_scale(ppb=10.0, max_length=5, bar_pad=8)
    # Scroll far past the 5-column alignment (content width ~ 8+50).
    profile.set_scroll_x(10_000)
    assert profile.viewport_peak() >= 0.0
    profile.set_scroll_x(-100)
    assert profile.viewport_peak() >= 0.0
    profile.close()
