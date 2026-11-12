"""Browser file table metadata."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.browser_file_metadata import (
    enrich_fasta_metadata,
    metadata_for_file,
    quick_metadata_for_file,
)


def test_fasta_metadata(tmp_path: Path) -> None:
    fasta = tmp_path / "GPTTR.fa"
    fasta.write_text(">GPTTR\n" + "ACGT" * 100 + "\n")
    meta = metadata_for_file(fasta)
    assert meta.molecule_type == "DNA"
    assert meta.sequence_count_text == "1"
    assert meta.length_text == "400"
    assert meta.size_text.endswith("B")
    assert meta.description_text == "GPTTR"
    assert meta.gc_text == "50.0"
    assert "20" in meta.modified_text or meta.modified_text != "—"


def test_quick_metadata_skips_fasta_read(tmp_path: Path) -> None:
    fasta = tmp_path / "x.fa"
    fasta.write_text(">seq\nACGT\n")
    quick = quick_metadata_for_file(fasta)
    assert quick.molecule_type == "DNA"
    assert quick.sequence_count_text == "…"
    assert quick.length_text == "…"
    assert quick.description_text == "…"
    rich = enrich_fasta_metadata(fasta)
    assert rich.sequence_count_text == "1"
    assert rich.length_text == "4"
    assert rich.description_text == "seq"


def test_protein_fasta_metadata(tmp_path: Path) -> None:
    fasta = tmp_path / "proteins.faa"
    fasta.write_text(">protein_1\nMKWVTFISLL\n")

    meta = metadata_for_file(fasta)

    assert meta.molecule_type == "Protein"
    assert meta.sequence_count_text == "1"
    assert meta.length_text == "10"
    assert meta.gc_text == "—"


def test_non_fasta_metadata(tmp_path: Path) -> None:
    path = tmp_path / "notes.txt"
    path.write_text("hello")
    meta = metadata_for_file(path)
    assert meta.molecule_type == "Text"
    assert meta.sequence_count_text == "—"
    assert meta.gc_text == "—"
    assert meta.length_text == "—"
    assert meta.size_text.endswith("B")
    assert meta.modified_text != "—"


def test_quick_metadata_fills_size_for_gff(tmp_path: Path) -> None:
    path = tmp_path / "features.gff"
    path.write_text("##gff-version 3\n")
    meta = quick_metadata_for_file(path)
    assert meta.molecule_type == "Annotation"
    assert meta.size_text.endswith("B")
    assert meta.modified_text != "—"
    assert meta.length_text == "—"
    assert meta.sequence_count_text == "—"


def test_msa_fasta_metadata_shows_column_count(tmp_path: Path) -> None:
    """Equal-length multi-FASTA (MSA): Length = alignment width."""
    fasta = tmp_path / "aln.fa"
    fasta.write_text(">a\nACGT--ACGT\n>b\nACGTA-ACGT\n>c\nACGTACACGT\n")
    meta = enrich_fasta_metadata(fasta)
    assert meta.sequence_count_text == "3"
    assert meta.length_text == "10"
    assert meta.sort_length == 10.0
    assert meta.gc_text == "—"


def test_unequal_multi_fasta_metadata_omits_length(tmp_path: Path) -> None:
    """Unaligned multi-FASTA keeps Length empty (no single column count)."""
    fasta = tmp_path / "multi.fa"
    fasta.write_text(">a\nACGT\n>b\nACGTACGTACGT\n")
    meta = enrich_fasta_metadata(fasta)
    assert meta.sequence_count_text == "2"
    assert meta.length_text == "—"
    assert meta.sort_length == -1.0