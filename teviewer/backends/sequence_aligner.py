"""Sequence alignment via external CLI tools (MAFFT, Clustal Omega)."""

from __future__ import annotations

import logging
import platform
import shutil
import subprocess
import tempfile
from pathlib import Path

from Bio import SeqIO

from teviewer.backends.gui_functions import align_sequences

OS_TYPE = platform.system()
logger = logging.getLogger(__name__)

_clustalo_openmp: bool | None = None


def alignment_backends() -> list[str]:
    """Installed aligners detected on PATH."""
    found: list[str] = []
    if OS_TYPE != "Windows" and shutil.which("mafft"):
        found.append("mafft")
    if shutil.which("clustalo"):
        found.append("clustalo")
    return found


def _fasta_sequence_count(path: Path) -> int:
    count = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def clustalo_supports_openmp() -> bool:
    """Whether the installed ``clustalo`` binary was built with OpenMP threading."""
    global _clustalo_openmp
    if _clustalo_openmp is not None:
        return _clustalo_openmp
    exe = shutil.which("clustalo")
    if not exe:
        _clustalo_openmp = False
        return False
    with tempfile.TemporaryDirectory(prefix="teviewer_clustalo_probe_") as tmp:
        inp = Path(tmp) / "probe.fa"
        out = Path(tmp) / "probe_out.fa"
        inp.write_text(">a\nACGT\n>b\nACGT\n", encoding="utf-8")
        try:
            proc = subprocess.run(
                [exe, "-i", str(inp), "-o", str(out), "--outfmt=fa", "--force", "--threads", "2"],
                capture_output=True,
                text=True,
                timeout=30,
                check=False,
            )
        except (OSError, subprocess.TimeoutExpired):
            _clustalo_openmp = False
            return False
        if proc.returncode == 0 and out.is_file():
            _clustalo_openmp = True
            return True
        _clustalo_openmp = False
        return False


def clustalo_effective_threads(requested: int) -> int:
    """Clamp thread count when Clustal Omega lacks OpenMP (common on macOS conda builds)."""
    requested = max(1, int(requested))
    if clustalo_supports_openmp():
        return requested
    return 1


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
    if _fasta_sequence_count(input_fasta) < 2:
        raise ValueError(
            f"{input_fasta.name} has fewer than 2 sequences. "
            "Alignment needs at least two sequences."
        )

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
                errors.append(str(exc))
        else:
            raise ValueError(f"Unknown aligner: {tool}")

    hint = (
        "Install MAFFT (macOS/Linux: conda install -c bioconda mafft) "
        "or Clustal Omega (conda install -c bioconda clustalo)."
    )
    detail = "\n".join(errors)
    if not detail:
        detail = "No aligner could run."
    raise RuntimeError(f"{detail}\n\n{hint}")


def _format_subprocess_error(tool: str, exc: subprocess.CalledProcessError) -> str:
    detail = (exc.stderr or exc.stdout or "").strip()
    if detail:
        first_line = detail.splitlines()[0].strip()
        return f"{tool}: {first_line}"
    return f"{tool} failed (exit {exc.returncode})."


def _run_clustalo(input_fasta: Path, output_file: Path, *, num_threads: int) -> None:
    threads = clustalo_effective_threads(num_threads)
    cmd = [
        "clustalo",
        "-i",
        str(input_fasta),
        "-o",
        str(output_file),
        "--outfmt=fa",
        "--force",
        "--threads",
        str(threads),
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except subprocess.CalledProcessError as exc:
        if threads > 1 and "openmp" in f"{exc.stderr or ''}{exc.stdout or ''}".lower():
            cmd[-1] = "1"
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        else:
            raise RuntimeError(_format_subprocess_error("Clustal Omega", exc)) from exc
    if not output_file.is_file():
        raise RuntimeError("Clustal Omega did not produce an output file.")
    if threads < num_threads:
        logger.info(
            "Clustal Omega ran with 1 thread (installed build lacks OpenMP support)."
        )
    logger.info("Clustal Omega alignment finished: %s", output_file.name)


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
