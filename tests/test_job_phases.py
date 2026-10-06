"""Logical progress phases for BLAST / Align / Flexidot / TE-Aid / Mapping."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.job_phases import (
    ALIGN_PHASES,
    BLAST_PHASES,
    FLEXIDOT_PHASES,
    MAP_PHASES,
    ONLINE_BLAST_PHASES,
    TEAID_PHASES,
    blast_units_for_fasta,
    emit_phase,
    tool_units_for_fasta,
)


def test_blast_units_single_vs_multi(tmp_path: Path) -> None:
    single = tmp_path / "one.fa"
    single.write_text(">a\nACGT\n", encoding="utf-8")
    multi = tmp_path / "many.fa"
    multi.write_text(">a\nACGT\n>b\nTTTT\n>c\nGGGG\n", encoding="utf-8")
    assert blast_units_for_fasta(single) == BLAST_PHASES
    assert blast_units_for_fasta(multi) == 3
    assert blast_units_for_fasta(single, online=True) == ONLINE_BLAST_PHASES


def test_tool_units_for_fasta(tmp_path: Path) -> None:
    fasta = tmp_path / "q.fa"
    fasta.write_text(">a\nACGT\n", encoding="utf-8")
    assert tool_units_for_fasta("flexidot", fasta) == FLEXIDOT_PHASES
    assert tool_units_for_fasta("teaid", fasta) == TEAID_PHASES
    assert tool_units_for_fasta("blast", fasta) == BLAST_PHASES
    assert ALIGN_PHASES == 3
    assert MAP_PHASES == 5


def test_emit_phase_reports_steps() -> None:
    events: list[tuple[int, int, str]] = []

    def progress(done: int, total: int, detail: str) -> None:
        events.append((done, total, detail))

    emit_phase(progress, 1, BLAST_PHASES, "running")
    emit_phase(None, 1, BLAST_PHASES, "ignored")
    assert events == [(1, BLAST_PHASES, "running")]
