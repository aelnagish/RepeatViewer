"""Sequence alignment via external CLI tools (MAFFT, Clustal Omega)."""

from __future__ import annotations

import logging
import platform
import shutil
import subprocess
from pathlib import Path

from Bio import SeqIO

from teviewer.backends.gui_functions import align_sequences

OS_TYPE = platform.system()
logger = logging.getLogger(__name__)


def alignment_backends() -> list[str]:
    """Installed aligners detected on PATH."""
    found: list[str] = []
    if OS_TYPE != "Windows" and shutil.which("mafft"):
        found.append("mafft")
    if shutil.which("clustalo"):
        found.append("clustalo")
    return found


def combined_alignment_output_stem(seq_files: list[Path]) -> str:
    """Basename for combined MSA: ``{first_stem}_msa_multi`` (same pattern as Flexidot multi)."""
    stem = seq_files[0].expanduser().resolve().stem
    return f"{stem}_msa_multi"


def _next_available_fasta(parent: Path, stem_base: str) -> Path:
    """``{stem}.fa``, then ``{stem}_2.fa``, ``{stem}_3.fa``, … (same index rule as Flexidot PNGs)."""
    candidate = parent / f"{stem_base}.fa"
    if not candidate.exists():
        return candidate
    index = 2
    while True:
        candidate = parent / f"{stem_base}_{index}.fa"
        if not candidate.exists():
            return candidate
        index += 1


def run_alignment(
    input_fasta: Path,
    *,
    backend: str = "auto",
    num_threads: int = 5,
) -> Path:
    """
    Align a FASTA/MSA file. Writes ``<stem>_msa.fa`` next to the input.

    ``backend``: ``auto`` (MAFFT then Clustal Omega), ``mafft``, or ``clustalo``.
    """
    input_fasta = input_fasta.expanduser().resolve()
    if input_fasta.suffix.lower() not in {".fa", ".fasta", ".fas"}:
        raise ValueError("Align only works on .fa / .fasta files.")

    output_file = input_fasta.parent / f"{input_fasta.stem}_msa.fa"
    order = ["mafft", "clustalo"] if backend == "auto" else [backend]
    errors: list[str] = []

    for tool in order:
        if tool == "mafft":
            if OS_TYPE == "Windows":
                errors.append("MAFFT is not supported on Windows in TEviewer.")
                continue
            if not shutil.which("mafft"):
                errors.append("MAFFT not found on PATH.")
                continue
            try:
                align_sequences(str(input_fasta), str(output_file), num_thread=num_threads)
                return output_file
            except Exception as exc:  # noqa: BLE001
                errors.append(f"MAFFT: {exc}")
        elif tool == "clustalo":
            if not shutil.which("clustalo"):
                errors.append("Clustal Omega (clustalo) not found on PATH.")
                continue
            try:
                _run_clustalo(input_fasta, output_file, num_threads=num_threads)
                return output_file
            except Exception as exc:  # noqa: BLE001
                errors.append(f"Clustal Omega: {exc}")
        else:
            raise ValueError(f"Unknown aligner: {tool}")

    hint = (
        "Install MAFFT (macOS/Linux: conda install -c bioconda mafft) "
        "or Clustal Omega (conda install -c bioconda clustalo)."
    )
    raise RuntimeError("No aligner available.\n" + "\n".join(errors) + f"\n\n{hint}")


def run_alignment_combined(
    seq_files: list[Path],
    *,
    backend: str = "auto",
    num_threads: int = 5,
) -> Path:
    """Merge selected FASTA files into one input, then align once."""
    if not seq_files:
        raise ValueError("No FASTA files to align.")
    resolved = [p.expanduser().resolve() for p in seq_files]
    parent = resolved[0].parent
    records = []
    used_ids: set[str] = set()
    for path in resolved:
        for record in SeqIO.parse(str(path), "fasta"):
            rid = record.id
            if rid in used_ids:
                base = rid
                n = 2
                while f"{base}_{n}" in used_ids:
                    n += 1
                rid = f"{base}_{n}"
                record.id = rid
                record.description = ""
            used_ids.add(rid)
            records.append(record)
    if not records:
        raise ValueError("Selected FASTA files contain no sequences.")
    if len(records) < 2:
        raise ValueError("Combined alignment needs at least two sequences.")
    output_stem = combined_alignment_output_stem(resolved)
    desired = _next_available_fasta(parent, output_stem)
    combined = parent / f".teviewer_{desired.stem}_input.fa"
    SeqIO.write(records, str(combined), "fasta")
    try:
        aligned = run_alignment(combined, backend=backend, num_threads=num_threads)
        if aligned.resolve() != desired.resolve() and aligned.is_file():
            if desired.is_file():
                desired.unlink()
            aligned.rename(desired)
            return desired
        return aligned
    finally:
        try:
            combined.unlink(missing_ok=True)
        except OSError:
            pass


def _run_clustalo(input_fasta: Path, output_file: Path, *, num_threads: int) -> None:
    cmd = [
        "clustalo",
        "-i",
        str(input_fasta),
        "-o",
        str(output_file),
        "--outfmt=fa",
        "--force",
        "--threads",
        str(num_threads),
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)
    if not output_file.is_file():
        raise RuntimeError("Clustal Omega did not produce an output file.")
    logger.info("Clustal Omega alignment finished: %s", output_file.name)
