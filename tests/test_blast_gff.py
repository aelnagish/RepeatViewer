"""BLAST table to GFF conversion."""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from teviewer.backends.blast_gff import blast_hits_to_gff_rows, write_blast_gff


def _write_blast_table(path: Path) -> None:
    df = pd.DataFrame(
        [
            {
                "qseqid": "TE1",
                "sseqid": "chr1",
                "pident": 95.0,
                "length": 100,
                "mismatch": 5,
                "gapopen": 0,
                "qstart": 1,
                "qend": 100,
                "sstart": 1000,
                "send": 1100,
                "sstrand": "plus",
                "evalue": 1e-50,
                "bitscore": 200.0,
            },
            {
                "qseqid": "TE1",
                "sseqid": "chr1",
                "pident": 90.0,
                "length": 80,
                "mismatch": 8,
                "gapopen": 0,
                "qstart": 1,
                "qend": 80,
                "sstart": 2000,
                "send": 2080,
                "sstrand": "minus",
                "evalue": 1e-40,
                "bitscore": 50.0,
            },
        ]
    )
    df.to_csv(path, sep="\t", index=False, header=False)


def test_blast_hits_to_gff_rows_filters_bit_score(tmp_path: Path) -> None:
    blast_txt = tmp_path / "reads.fasta_blast.txt"
    _write_blast_table(blast_txt)
    rows = blast_hits_to_gff_rows(blast_txt, minimum_bit_score=100.0)
    assert len(rows) == 1
    assert rows[0][2] == "TE1"
    assert rows[0][6] == "+"


def test_write_blast_gff(tmp_path: Path) -> None:
    blast_txt = tmp_path / "reads.fasta_blast.txt"
    _write_blast_table(blast_txt)
    gff_path = write_blast_gff(blast_txt, minimum_bit_score=0.0)
    text = gff_path.read_text(encoding="utf-8")
    assert text.startswith("##gff-version 3")
    assert "TE1" in text
    assert "chr1" in text
