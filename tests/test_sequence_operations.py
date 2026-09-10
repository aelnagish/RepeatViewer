from pathlib import Path

from Bio import SeqIO

from repeatviewer.backends.orf_detect import OrfDetectResult
from repeatviewer.backends.sequence_features import SequenceFeature
from repeatviewer.backends.sequence_operations import (
    write_predicted_orfs_fasta,
    write_reverse_complement_fasta,
    write_translated_orfs_fasta,
)


def _write_fasta(path: Path, sequence: str = "ATGAAATAA") -> None:
    path.write_text(f">seq1\n{sequence}\n", encoding="utf-8")


def test_reverse_complement_fasta(tmp_path: Path) -> None:
    source = tmp_path / "input.fa"
    _write_fasta(source, "ATGC")

    output = write_reverse_complement_fasta(source)

    record = next(SeqIO.parse(output, "fasta"))
    assert output.name == "input_rc.fa"
    assert str(record.seq) == "GCAT"


def test_orf_prediction_and_translation_use_detected_strand(
    tmp_path: Path, monkeypatch
) -> None:
    source = tmp_path / "input.fa"
    _write_fasta(source, "ATGAAATAA")
    feature = SequenceFeature(1, 9, "ORF1", "ORF", color="#27ae60", strand="+")
    monkeypatch.setattr(
        "repeatviewer.backends.sequence_operations.detect_orf_features",
        lambda *_args, **_kwargs: OrfDetectResult((feature,), "Detected 1 ORF."),
    )

    dna_output = write_predicted_orfs_fasta(source)
    protein_output = write_translated_orfs_fasta(source)

    dna = next(SeqIO.parse(dna_output, "fasta"))
    protein = next(SeqIO.parse(protein_output, "fasta"))
    assert str(dna.seq) == "ATGAAATAA"
    assert str(protein.seq) == "MK"
    assert protein_output.name == "input_translated.faa"
