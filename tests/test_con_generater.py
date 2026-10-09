"""Consensus generation (con_generater): nucleotide and protein MSAs."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO

from repeatviewer.backends.teaid_plotter import (
    _alignment_dumb_consensus,
    con_generater,
)
from repeatviewer.job_options import MAFFT_SEQUENCE_NUCLEOTIDE, MAFFT_SEQUENCE_PROTEIN


def test_alignment_dumb_consensus_matches_biopython_rules() -> None:
    # Unanimous → letters; below threshold / ties / all-gap → ambiguous.
    assert _alignment_dumb_consensus(
        ["ACGT", "ACGT", "ACGT"], threshold=0.8, ambiguous="N"
    ) == "ACGT"
    assert _alignment_dumb_consensus(
        ["ACGT", "ACGA", "ACGC"], threshold=0.8, ambiguous="N"
    ) == "ACGN"
    assert _alignment_dumb_consensus(
        ["A-GT", "A-GT", "--GT"], threshold=0.5, ambiguous="N"
    ) == "ANGT"
    assert _alignment_dumb_consensus(
        ["AC", "TG"], threshold=0.5, ambiguous="N"
    ) == "NN"


def test_con_generater_nucleotide_majority(tmp_path: Path) -> None:
    fasta = tmp_path / "aln.fa"
    fasta.write_text(
        ">a\nACGT\n>b\nACGT\n>c\nACGT\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    out_path, length = con_generater(
        str(fasta),
        str(out_dir),
        threshold=0.8,
        ambiguous="N",
        sequence_type=MAFFT_SEQUENCE_NUCLEOTIDE,
    )
    assert length == 4
    rec = next(SeqIO.parse(out_path, "fasta"))
    assert str(rec.seq) == "ACGT"


def test_con_generater_protein_threshold(tmp_path: Path) -> None:
    fasta = tmp_path / "prot.fa"
    fasta.write_text(
        ">a\nMKTA\n>b\nMKTA\n>c\nMKTV\n",
        encoding="utf-8",
    )
    out_dir = tmp_path / "out"
    out_dir.mkdir()
    out_path, _length = con_generater(
        str(fasta),
        str(out_dir),
        threshold=0.8,
        ambiguous="X",
        sequence_type=MAFFT_SEQUENCE_PROTEIN,
    )
    rec = next(SeqIO.parse(out_path, "fasta"))
    assert str(rec.seq) == "MKTX"
