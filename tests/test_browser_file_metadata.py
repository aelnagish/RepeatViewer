"""Browser file table metadata."""

from __future__ import annotations

from pathlib import Path

from teviewer.browser_file_metadata import (
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
