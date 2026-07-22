"""Sequence alignment using TEviewer's bundled or system MAFFT/MUSCLE."""

from __future__ import annotations

import logging
import os
import platform
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path

from Bio import SeqIO

from teviewer.job_options import (
    ALIGNER_MAFFT,
    ALIGNER_MUSCLE,
    MAFFT_SEQUENCE_NUCLEOTIDE,
    MAFFT_SEQUENCE_PROTEIN,
    MAFFT_STRATEGY_AUTO,
    MAFFT_STRATEGY_CUSTOM,
    MAFFT_STRATEGY_EINSI,
    MAFFT_STRATEGY_FFT_NS_1,
    MAFFT_STRATEGY_FFT_NS_2,
    MAFFT_STRATEGY_GINSI,
    MAFFT_STRATEGY_LINSI,
)
from teviewer.external_process import run_logged_process

OS_TYPE = platform.system()
logger = logging.getLogger(__name__)

def _mafft_bundle_root() -> Path:
    if hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS) / "teviewer" / "backends" / "mafft"
    return Path(__file__).resolve().parent / "mafft"


def _muscle_bundle_root() -> Path:
    if hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS) / "teviewer" / "backends" / "muscle"
    return Path(__file__).resolve().parent / "muscle"


def bundled_mafft_command(
    os_type: str | None = None,
    machine: str | None = None,
) -> tuple[list[str], dict[str, str]] | None:
    """Command prefix and environment for TEviewer's bundled MAFFT."""
    system = os_type or OS_TYPE
    architecture = (machine or platform.machine()).lower()
    root = _mafft_bundle_root()
    env = os.environ.copy()
    if system == "Darwin":
        executable = root / "macos" / "bin" / "mafft"
        libexec = root / "macos" / "libexec"
        if executable.is_file() and libexec.is_dir():
            env["MAFFT_BINARIES"] = str(libexec)
            return [str(executable)], env
    elif system == "Windows":
        if architecture not in {"amd64", "x86_64"}:
            return None
        wrapper = root / "windows" / "mafft.bat"
        if wrapper.is_file():
            return [str(wrapper)], env
    elif system == "Linux":
        if architecture not in {"amd64", "x86_64"}:
            return None
        wrapper = root / "linux" / "mafft.bat"
        if wrapper.is_file():
            return [str(wrapper)], env
    return None


def mafft_command(*, use_system_mafft: bool = False) -> tuple[list[str], dict[str, str]] | None:
    if use_system_mafft:
        executable = shutil.which("mafft")
        return ([executable], os.environ.copy()) if executable else None
    return bundled_mafft_command()


def bundled_muscle_command(
    os_type: str | None = None,
    machine: str | None = None,
) -> tuple[list[str], dict[str, str]] | None:
    """Command prefix and environment for TEviewer's bundled MUSCLE 5."""
    system = os_type or OS_TYPE
    architecture = (machine or platform.machine()).lower()
    root = _muscle_bundle_root()
    env = os.environ.copy()
    if system == "Darwin":
        if architecture not in {"arm64", "aarch64"}:
            return None
        executable = root / "macos" / "arm64" / "muscle"
    elif system == "Linux":
        arch_dir = "aarch64" if architecture in {"arm64", "aarch64"} else "x86_64"
        if architecture not in {"arm64", "aarch64", "amd64", "x86_64"}:
            return None
        executable = root / "linux" / arch_dir / "muscle"
    elif system == "Windows":
        if architecture not in {"amd64", "x86_64"}:
            return None
        executable = root / "windows" / "muscle.exe"
    else:
        return None
    return ([str(executable)], env) if executable.is_file() else None


def muscle_command(*, use_system_muscle: bool = False) -> tuple[list[str], dict[str, str]] | None:
    if use_system_muscle:
        executable = shutil.which("muscle")
        return ([executable], os.environ.copy()) if executable else None
    return bundled_muscle_command()


def alignment_backends(
    *,
    use_system_mafft: bool = False,
    use_system_muscle: bool = False,
) -> list[str]:
    """Return the selected bundled/system alignment engines that are available."""
    engines = []
    if mafft_command(use_system_mafft=use_system_mafft):
        engines.append(ALIGNER_MAFFT)
    if muscle_command(use_system_muscle=use_system_muscle):
        engines.append(ALIGNER_MUSCLE)
    return engines


def _fasta_sequence_count(path: Path) -> int:
    count = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


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
    num_threads: int = 5,
    use_system_mafft: bool = False,
    use_system_muscle: bool = False,
    engine: str = ALIGNER_MAFFT,
    strategy: str = MAFFT_STRATEGY_FFT_NS_1,
    sequence_type: str = "auto",
    retree: int = 1,
    max_iterate: int = 0,
    gap_open_penalty: float = 1.53,
    offset: float = 0.0,
    reorder: bool = False,
    adjust_direction: bool = False,
    muscle_mode: str = "align",
    muscle_num_threads: int = 5,
) -> Path:
    """
    Align a FASTA/MSA file. Writes ``<stem>_msa.fa`` next to the input.

    ``engine`` selects MAFFT or MUSCLE 5.
    """
    input_fasta = input_fasta.expanduser().resolve()
    if input_fasta.suffix.lower() not in {".fa", ".fasta", ".fas", ".faa"}:
        raise ValueError("Align only works on .fa / .fasta files.")
    if _fasta_sequence_count(input_fasta) < 2:
        raise ValueError(
            f"{input_fasta.name} has fewer than 2 sequences. "
            "Alignment needs at least two sequences."
        )

    output_file = input_fasta.parent / f"{input_fasta.stem}_msa.fa"
    if engine == ALIGNER_MUSCLE:
        resolved = muscle_command(use_system_muscle=use_system_muscle)
        system_selected = use_system_muscle
        tool_name = "MUSCLE"
    elif engine == ALIGNER_MAFFT:
        resolved = mafft_command(use_system_mafft=use_system_mafft)
        system_selected = use_system_mafft
        tool_name = "MAFFT"
    else:
        raise ValueError(f"Unknown alignment engine: {engine}")
    if resolved is None:
        source = "system PATH" if system_selected else "TEviewer bundle"
        raise RuntimeError(
            f"{tool_name} was not found in the {source}. Disable 'Use system "
            f"{tool_name}' in Preferences to use TEviewer's bundled copy, or install "
            f"{tool_name} on PATH."
        )
    if engine == ALIGNER_MUSCLE:
        _run_muscle(
            input_fasta,
            output_file,
            command=resolved[0],
            env=resolved[1],
            mode=muscle_mode,
            num_threads=muscle_num_threads,
        )
    else:
        _run_mafft(
            input_fasta,
            output_file,
            command=resolved[0],
            env=resolved[1],
            num_threads=num_threads,
            strategy=strategy,
            sequence_type=sequence_type,
            retree=retree,
            max_iterate=max_iterate,
            gap_open_penalty=gap_open_penalty,
            offset=offset,
            reorder=reorder,
            adjust_direction=adjust_direction,
        )
    return output_file


def _format_subprocess_error(tool: str, exc: subprocess.CalledProcessError) -> str:
    detail = (exc.stderr or exc.stdout or "").strip()
    if detail:
        first_line = detail.splitlines()[0].strip()
        return f"{tool}: {first_line}"
    return f"{tool} failed (exit {exc.returncode})."


def _run_mafft(
    input_fasta: Path,
    output_file: Path,
    *,
    command: list[str],
    env: dict[str, str],
    num_threads: int,
    strategy: str,
    sequence_type: str,
    retree: int,
    max_iterate: int,
    gap_open_penalty: float,
    offset: float,
    reorder: bool,
    adjust_direction: bool,
) -> None:
    arguments = _mafft_arguments(
        num_threads=num_threads,
        strategy=strategy,
        sequence_type=sequence_type,
        retree=retree,
        max_iterate=max_iterate,
        gap_open_penalty=gap_open_penalty,
        offset=offset,
        reorder=reorder,
        adjust_direction=adjust_direction,
    )
    with tempfile.NamedTemporaryFile(
        prefix=".teviewer_mafft_progress_",
        suffix=".log",
        dir=input_fasta.parent,
        delete=False,
    ) as progress_handle:
        progress_path = Path(progress_handle.name)
    cmd = [
        *command,
        *arguments,
        "--progress",
        str(progress_path),
        str(input_fasta),
    ]
    try:
        try:
            result = run_logged_process(
                cmd,
                tool="MAFFT",
                inputs=[input_fasta],
                cwd=input_fasta.parent,
                check=True,
                env=env,
                stdout_is_result=True,
                output_path=output_file,
                progress_file=progress_path,
            )
        finally:
            progress_path.unlink(missing_ok=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(_format_subprocess_error("MAFFT", exc)) from exc
    if not result.stdout.strip():
        raise RuntimeError("MAFFT completed without producing an alignment.")
    output_file.write_text(result.stdout, encoding="utf-8")
    logger.info("MAFFT alignment finished: %s", output_file.name)


def _run_muscle(
    input_fasta: Path,
    output_file: Path,
    *,
    command: list[str],
    env: dict[str, str],
    mode: str,
    num_threads: int,
) -> None:
    if mode not in {"align", "super5"}:
        raise ValueError(f"Unknown MUSCLE mode: {mode}")
    cmd = [
        *command,
        f"-{mode}",
        str(input_fasta),
        "-output",
        str(output_file),
        "-threads",
        str(max(1, num_threads)),
    ]
    try:
        result = run_logged_process(
            cmd,
            tool="MUSCLE",
            inputs=[input_fasta],
            cwd=input_fasta.parent,
            check=True,
            env=env,
            output_path=output_file,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(_format_subprocess_error("MUSCLE", exc)) from exc
    if not output_file.is_file() or not output_file.stat().st_size:
        raise RuntimeError("MUSCLE completed without producing an alignment.")
    logger.info("MUSCLE alignment finished: %s", output_file.name)


def _mafft_arguments(
    *,
    num_threads: int,
    strategy: str,
    sequence_type: str,
    retree: int,
    max_iterate: int,
    gap_open_penalty: float,
    offset: float,
    reorder: bool,
    adjust_direction: bool,
) -> list[str]:
    """Translate persisted alignment options into MAFFT CLI arguments."""
    strategy_args = {
        MAFFT_STRATEGY_AUTO: ["--auto"],
        MAFFT_STRATEGY_FFT_NS_1: ["--retree", "1", "--maxiterate", "0"],
        MAFFT_STRATEGY_FFT_NS_2: ["--retree", "2", "--maxiterate", "0"],
        MAFFT_STRATEGY_LINSI: ["--localpair", "--maxiterate", "1000"],
        MAFFT_STRATEGY_GINSI: ["--globalpair", "--maxiterate", "1000"],
        MAFFT_STRATEGY_EINSI: ["--genafpair", "--maxiterate", "1000"],
        MAFFT_STRATEGY_CUSTOM: [
            "--retree",
            str(max(1, retree)),
            "--maxiterate",
            str(max(0, max_iterate)),
        ],
    }
    if strategy not in strategy_args:
        raise ValueError(f"Unknown MAFFT strategy: {strategy}")
    args = [*strategy_args[strategy]]
    if sequence_type == MAFFT_SEQUENCE_NUCLEOTIDE:
        args.append("--nuc")
    elif sequence_type == MAFFT_SEQUENCE_PROTEIN:
        args.append("--amino")
    elif sequence_type != "auto":
        raise ValueError(f"Unknown MAFFT sequence type: {sequence_type}")
    args.extend(
        [
            "--thread",
            str(max(1, num_threads)),
            "--op",
            f"{gap_open_penalty:g}",
            "--ep",
            f"{offset:g}",
        ]
    )
    if reorder:
        args.append("--reorder")
    if adjust_direction:
        if sequence_type == MAFFT_SEQUENCE_PROTEIN:
            raise ValueError("Reverse-complement adjustment is only valid for DNA/RNA.")
        args.append("--adjustdirectionaccurately")
    return args


def run_alignment_combined(
    seq_files: list[Path],
    *,
    num_threads: int = 5,
    use_system_mafft: bool = False,
    use_system_muscle: bool = False,
    engine: str = ALIGNER_MAFFT,
    strategy: str = MAFFT_STRATEGY_FFT_NS_1,
    sequence_type: str = "auto",
    retree: int = 1,
    max_iterate: int = 0,
    gap_open_penalty: float = 1.53,
    offset: float = 0.0,
    reorder: bool = False,
    adjust_direction: bool = False,
    muscle_mode: str = "align",
    muscle_num_threads: int = 5,
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
        aligned = run_alignment(
            combined,
            num_threads=num_threads,
            use_system_mafft=use_system_mafft,
            use_system_muscle=use_system_muscle,
            engine=engine,
            strategy=strategy,
            sequence_type=sequence_type,
            retree=retree,
            max_iterate=max_iterate,
            gap_open_penalty=gap_open_penalty,
            offset=offset,
            reorder=reorder,
            adjust_direction=adjust_direction,
            muscle_mode=muscle_mode,
            muscle_num_threads=muscle_num_threads,
        )
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
