"""Standalone nucleotide transformation, ORF prediction, and translation tools."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from teviewer.backends.auto_annotation_settings import load_auto_annotation_settings
from teviewer.backends.orf_detect import detect_orf_features


def _records(path: Path) -> list[SeqRecord]:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise RuntimeError(f"No FASTA sequences found in {path.name}.")
    return records


def write_reverse_complement_fasta(path: Path) -> Path:
    """Reverse-complement every FASTA record and preserve its identifier."""
    source = path.expanduser().resolve()
    output = source.with_name(f"{source.stem}_reverse_complement.fa")
    transformed = [
        SeqRecord(
            record.seq.reverse_complement(),
            id=record.id,
            description=f"{record.description} reverse_complement".strip(),
        )
        for record in _records(source)
    ]
    SeqIO.write(transformed, str(output), "fasta")
    return output


def predicted_orf_records(path: Path, *, translate: bool = False) -> list[SeqRecord]:
    """Predict ORFs with the Annotation panel settings and return FASTA records."""
    source = path.expanduser().resolve()
    settings = load_auto_annotation_settings().orf
    output_records: list[SeqRecord] = []
    for record in _records(source):
        sequence = str(record.seq).upper()
        result = detect_orf_features(record.id, sequence, settings=settings)
        for feature in result.features:
            dna = Seq(sequence[feature.start - 1 : feature.end])
            if feature.strand == "-":
                dna = dna.reverse_complement()
            seq = dna.translate(to_stop=False).rstrip("*") if translate else dna
            strand = feature.strand or "+"
            molecule = "protein" if translate else "DNA"
            output_records.append(
                SeqRecord(
                    Seq(str(seq)),
                    id=f"{record.id}_{feature.name}",
                    description=(
                        f"source={record.id} start={feature.start} end={feature.end} "
                        f"strand={strand} molecule={molecule}"
                    ),
                )
            )
    if not output_records:
        raise RuntimeError(f"No ORFs were predicted in {source.name}.")
    return output_records


def write_predicted_orfs_fasta(path: Path) -> Path:
    """Write predicted ORF nucleotide sequences to a FASTA file."""
    source = path.expanduser().resolve()
    output = source.with_name(f"{source.stem}_orfs.fa")
    SeqIO.write(predicted_orf_records(source), str(output), "fasta")
    return output


def write_translated_orfs_fasta(path: Path) -> Path:
    """Predict ORFs, then translate them with Biopython into protein FASTA."""
    source = path.expanduser().resolve()
    output = source.with_name(f"{source.stem}_translated.faa")
    SeqIO.write(predicted_orf_records(source, translate=True), str(output), "fasta")
    return output
