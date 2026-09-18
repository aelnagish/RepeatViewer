"""Tests for Sequence Annotation batch helpers."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.te_annotation_batch import existing_annotation_gff_paths


def test_existing_annotation_gff_paths_finds_rv_gff3(tmp_path: Path) -> None:
    fasta = tmp_path / "TE1.fa"
    fasta.write_text(">TE1\nACGT\n", encoding="utf-8")
    gff = tmp_path / "TE1.RV.gff3"
    gff.write_text("##gff-version 3\n", encoding="utf-8")
    other = tmp_path / "TE2.fa"
    other.write_text(">TE2\nACGT\n", encoding="utf-8")

    found = existing_annotation_gff_paths([fasta, other])
    assert found == [gff.resolve()]


def test_existing_annotation_gff_paths_empty_when_missing(tmp_path: Path) -> None:
    fasta = tmp_path / "TE1.fa"
    fasta.write_text(">TE1\nACGT\n", encoding="utf-8")
    assert existing_annotation_gff_paths([fasta]) == []
