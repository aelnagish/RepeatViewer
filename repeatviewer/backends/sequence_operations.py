"""Standalone nucleotide transformation, ORF prediction, and translation tools."""

from __future__ import annotations

from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from repeatviewer.backends.auto_annotation_settings import OrfSettings
from repeatviewer.backends.orf_detect import detect_orf_features


def _records(path: Path) -> list[SeqRecord]:
    records = list(SeqIO.parse(str(path), "fasta"))
    if not records:
        raise RuntimeError(f"No FASTA sequences found in {path.name}.")
    return records


def write_reverse_complement_fasta(path: Path) -> Path:
    """Reverse-complement every FASTA record and preserve its identifier."""
    from repeatviewer.tool_output_paths import kept_input_name

    source = path.expanduser().resolve()
    output = source.with_name(f"{kept_input_name(source)}_rc.fa")
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


def predicted_orf_records(
    path: Path,
    *,
    translate: bool = False,
    settings: OrfSettings | None = None,
) -> list[SeqRecord]:
    """Predict ORFs with Sequence → ORF Prediction settings and return FASTA records."""
    source = path.expanduser().resolve()
    orf_settings = settings or OrfSettings()
    output_records: list[SeqRecord] = []
    for record in _records(source):
        sequence = str(record.seq).upper()
        result = detect_orf_features(record.id, sequence, settings=orf_settings)
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


def write_predicted_orfs_fasta(
    path: Path,
    *,
    settings: OrfSettings | None = None,
) -> Path:
    """Write predicted ORF nucleotide sequences to a FASTA file."""
    from repeatviewer.tool_output_paths import kept_input_name

    source = path.expanduser().resolve()
    output = source.with_name(f"{kept_input_name(source)}_orfs.fa")
    SeqIO.write(predicted_orf_records(source, settings=settings), str(output), "fasta")
    return output


def write_translated_orfs_fasta(
    path: Path,
    *,
    settings: OrfSettings | None = None,
) -> Path:
    """Predict ORFs, then translate them with Biopython into protein FASTA."""
    from repeatviewer.tool_output_paths import kept_input_name

    source = path.expanduser().resolve()
    output = source.with_name(f"{kept_input_name(source)}_translated.faa")
    SeqIO.write(
        predicted_orf_records(source, translate=True, settings=settings),
        str(output),
        "fasta",
    )
    return output
