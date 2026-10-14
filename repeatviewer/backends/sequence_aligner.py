"""Sequence alignment using RepeatViewer's bundled or system MAFFT/MUSCLE."""

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

from repeatviewer.job_options import (
    ALIGNER_MAFFT,
    ALIGNER_MUSCLE,
    MAFFT_ADD_METHOD_ADD,
    MAFFT_ADD_METHOD_ADDFRAGMENTS,
    MAFFT_ADD_METHOD_ADDFULL,
    MAFFT_ADD_METHOD_ADDLONG,
    MAFFT_ADD_METHODS,
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
from repeatviewer.external_process import run_logged_process

OS_TYPE = platform.system()
logger = logging.getLogger(__name__)

def _mafft_bundle_root() -> Path:
    if hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS) / "repeatviewer" / "backends" / "mafft"
    return Path(__file__).resolve().parent / "mafft"


def _muscle_bundle_root() -> Path:
    if hasattr(sys, "_MEIPASS"):
        return Path(sys._MEIPASS) / "repeatviewer" / "backends" / "muscle"
    return Path(__file__).resolve().parent / "muscle"


def bundled_mafft_command(
    os_type: str | None = None,
    machine: str | None = None,
) -> tuple[list[str], dict[str, str]] | None:
    """Command prefix and environment for RepeatViewer's bundled MAFFT."""
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
    """Command prefix and environment for RepeatViewer's bundled MUSCLE 5."""
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


def alignment_output_stem(
    input_fasta: Path | str,
    *,
    engine: str = ALIGNER_MAFFT,
) -> str:
    """Basename stem for Align outputs: ``{input.name}_mafft`` or ``_muscle``."""
    from repeatviewer.tool_output_paths import kept_input_name

    tool = "muscle" if engine == ALIGNER_MUSCLE else "mafft"
    if engine not in (ALIGNER_MAFFT, ALIGNER_MUSCLE):
        raise ValueError(f"Unknown alignment engine: {engine}")
    return f"{kept_input_name(input_fasta)}_{tool}"


def combined_alignment_output_stem(
    seq_files: list[Path],
    *,
    engine: str = ALIGNER_MAFFT,
) -> str:
    """Basename stem for combined MSA: ``{first.name}_cmafft`` / ``_cmuscle``."""
    if not seq_files:
        raise ValueError("No FASTA files to align.")
    from repeatviewer.tool_output_paths import kept_input_name

    if engine == ALIGNER_MUSCLE:
        tool = "cmuscle"
    elif engine == ALIGNER_MAFFT:
        tool = "cmafft"
    else:
        raise ValueError(f"Unknown alignment engine: {engine}")
    return f"{kept_input_name(seq_files[0])}_{tool}"


def _next_available_fasta(
    parent: Path,
    stem_base: str,
    *,
    first_conflict_index: int = 1,
) -> Path:
    """``{stem}.fa``, then ``{stem}_1.fa``, ``{stem}_2.fa``, … on conflict."""
    candidate = parent / f"{stem_base}.fa"
    if not candidate.exists():
        return candidate
    index = max(1, int(first_conflict_index))
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
    muscle_perm: str = "none",
    muscle_perturb: int = 0,
    cancel=None,
) -> Path:
    """
    Align a FASTA/MSA file. Writes ``<name>_mafft.fa`` or ``<name>_muscle.fa``
    next to the input (``_1``, ``_2``, … before ``.fa`` if the name exists).

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

    output_file = _next_available_fasta(
        input_fasta.parent, alignment_output_stem(input_fasta, engine=engine)
    )
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
        source = "system PATH" if system_selected else "RepeatViewer bundle"
        raise RuntimeError(
            f"{tool_name} was not found in the {source}. Disable 'Use system "
            f"{tool_name}' in Preferences to use RepeatViewer's bundled copy, or install "
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
            perm=muscle_perm,
            perturb=muscle_perturb,
            cancel=cancel,
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
            cancel=cancel,
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
    cancel=None,
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
        prefix=".repeatviewer_mafft_progress_",
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
                cancel=cancel,
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
    perm: str = "none",
    perturb: int = 0,
    cancel=None,
) -> None:
    if mode not in {"align", "super5"}:
        raise ValueError(f"Unknown MUSCLE mode: {mode}")
    cmd = [
        *command,
        *_muscle_arguments(
            mode=mode,
            num_threads=num_threads,
            perm=perm,
            perturb=perturb,
            input_fasta=input_fasta,
            output_file=output_file,
        ),
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
            cancel=cancel,
        )
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(_format_subprocess_error("MUSCLE", exc)) from exc
    if not output_file.is_file() or not output_file.stat().st_size:
        raise RuntimeError("MUSCLE completed without producing an alignment.")
    logger.info("MUSCLE alignment finished: %s", output_file.name)


def _muscle_arguments(
    *,
    mode: str,
    num_threads: int,
    perm: str,
    perturb: int,
    input_fasta: Path | str,
    output_file: Path | str,
) -> list[str]:
    """Translate persisted MUSCLE options into MUSCLE 5 CLI arguments."""
    if mode not in {"align", "super5"}:
        raise ValueError(f"Unknown MUSCLE mode: {mode}")
    args = [
        f"-{mode}",
        str(input_fasta),
        "-output",
        str(output_file),
        "-threads",
        str(max(1, num_threads)),
    ]
    perm_value = str(perm or "none").lower()
    if perm_value not in {"none", "abc", "acb", "bca"}:
        perm_value = "none"
    if perm_value != "none":
        args.extend(["-perm", perm_value])
    perturb_value = max(0, int(perturb))
    if perturb_value > 0:
        args.extend(["-perturb", str(perturb_value)])
    return args


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


def add_to_msa_output_stem(
    msa_fasta: Path | str,
    *,
    engine: str = ALIGNER_MAFFT,
) -> str:
    """Basename stem for add-to-MSA outputs: ``{msa.name}_add_mafft`` / ``_add_muscle``."""
    from repeatviewer.tool_output_paths import kept_input_name

    if engine == ALIGNER_MUSCLE:
        tool = "add_muscle"
    elif engine == ALIGNER_MAFFT:
        tool = "add_mafft"
    else:
        raise ValueError(f"Unknown alignment engine: {engine}")
    return f"{kept_input_name(msa_fasta)}_{tool}"


def run_alignment_add_to_msa(
    sequences_fasta: Path,
    msa_fasta: Path,
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
    muscle_perm: str = "none",
    muscle_perturb: int = 0,
    mafft_add_method: str = MAFFT_ADD_METHOD_ADD,
    mafft_add_keeplength: bool = False,
    cancel=None,
) -> Path:
    """Add unaligned sequences into an existing MSA.

    MAFFT uses ``--add`` / ``--addfragments`` / ``--addfull`` / ``--addlong``
    (optional ``--keeplength``). MUSCLE 5 has no profile-add; it degaps the MSA,
    concatenates the new sequences, and re-aligns the combined set.
    """
    sequences_fasta = sequences_fasta.expanduser().resolve()
    msa_fasta = msa_fasta.expanduser().resolve()
    if sequences_fasta.resolve() == msa_fasta.resolve():
        raise ValueError("Sequences to add and the existing MSA must be different files.")
    if _fasta_sequence_count(sequences_fasta) < 1:
        raise ValueError(f"{sequences_fasta.name} contains no sequences to add.")
    if _fasta_sequence_count(msa_fasta) < 1:
        raise ValueError(f"{msa_fasta.name} contains no sequences.")

    output_file = _next_available_fasta(
        msa_fasta.parent, add_to_msa_output_stem(msa_fasta, engine=engine)
    )
    if engine == ALIGNER_MUSCLE:
        return _run_muscle_add_to_msa(
            sequences_fasta,
            msa_fasta,
            output_file,
            use_system_muscle=use_system_muscle,
            mode=muscle_mode,
            num_threads=muscle_num_threads,
            perm=muscle_perm,
            perturb=muscle_perturb,
            cancel=cancel,
        )
    if engine != ALIGNER_MAFFT:
        raise ValueError(f"Unknown alignment engine: {engine}")
    resolved = mafft_command(use_system_mafft=use_system_mafft)
    if resolved is None:
        source = "system PATH" if use_system_mafft else "RepeatViewer bundle"
        raise RuntimeError(
            f"MAFFT was not found in the {source}. Disable 'Use system "
            "MAFFT' in Preferences to use RepeatViewer's bundled copy, or install "
            "MAFFT on PATH."
        )
    _run_mafft_add_to_msa(
        sequences_fasta,
        msa_fasta,
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
        add_method=mafft_add_method,
        keeplength=mafft_add_keeplength,
        cancel=cancel,
    )
    return output_file


def _run_mafft_add_to_msa(
    sequences_fasta: Path,
    msa_fasta: Path,
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
    add_method: str,
    keeplength: bool,
    cancel=None,
) -> None:
    method = str(add_method or MAFFT_ADD_METHOD_ADD).lower()
    if method not in MAFFT_ADD_METHODS:
        method = MAFFT_ADD_METHOD_ADD
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
    add_flag = f"--{method}"
    with tempfile.NamedTemporaryFile(
        prefix=".repeatviewer_mafft_progress_",
        suffix=".log",
        dir=msa_fasta.parent,
        delete=False,
    ) as progress_handle:
        progress_path = Path(progress_handle.name)
    cmd = [
        *command,
        *arguments,
        add_flag,
        str(sequences_fasta),
    ]
    if keeplength and method in {
        MAFFT_ADD_METHOD_ADD,
        MAFFT_ADD_METHOD_ADDFRAGMENTS,
        MAFFT_ADD_METHOD_ADDFULL,
        MAFFT_ADD_METHOD_ADDLONG,
    }:
        cmd.append("--keeplength")
    cmd.extend(
        [
            "--progress",
            str(progress_path),
            str(msa_fasta),
        ]
    )
    try:
        try:
            result = run_logged_process(
                cmd,
                tool="MAFFT",
                inputs=[sequences_fasta, msa_fasta],
                cwd=msa_fasta.parent,
                check=True,
                env=env,
                stdout_is_result=True,
                output_path=output_file,
                progress_file=progress_path,
                cancel=cancel,
            )
        finally:
            progress_path.unlink(missing_ok=True)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(_format_subprocess_error("MAFFT", exc)) from exc
    if not result.stdout.strip():
        raise RuntimeError("MAFFT completed without producing an alignment.")
    output_file.write_text(result.stdout, encoding="utf-8")
    logger.info("MAFFT add-to-MSA finished: %s", output_file.name)


def _run_muscle_add_to_msa(
    sequences_fasta: Path,
    msa_fasta: Path,
    output_file: Path,
    *,
    use_system_muscle: bool,
    mode: str,
    num_threads: int,
    perm: str,
    perturb: int,
    cancel=None,
) -> Path:
    """MUSCLE 5 has no profile-add: degap MSA + new sequences, then re-align."""
    resolved = muscle_command(use_system_muscle=use_system_muscle)
    if resolved is None:
        source = "system PATH" if use_system_muscle else "RepeatViewer bundle"
        raise RuntimeError(
            f"MUSCLE was not found in the {source}. Disable 'Use system "
            "MUSCLE' in Preferences to use RepeatViewer's bundled copy, or install "
            "MUSCLE on PATH."
        )
    records = []
    used_ids: set[str] = set()
    for path in (msa_fasta, sequences_fasta):
        for record in SeqIO.parse(str(path), "fasta"):
            # Degap so MUSCLE re-aligns from sequences (MSA columns are not kept).
            record.seq = record.seq.replace("-", "").replace(".", "")
            if not str(record.seq):
                continue
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
    if len(records) < 2:
        raise ValueError(
            "Add-to-MSA with MUSCLE needs at least two sequences after merging "
            "the MSA and the sequences to add."
        )
    combined = msa_fasta.parent / f".repeatviewer_{output_file.stem}_input.fa"
    SeqIO.write(records, str(combined), "fasta")
    try:
        _run_muscle(
            combined,
            output_file,
            command=resolved[0],
            env=resolved[1],
            mode=mode,
            num_threads=num_threads,
            perm=perm,
            perturb=perturb,
            cancel=cancel,
        )
    finally:
        try:
            combined.unlink(missing_ok=True)
        except OSError:
            pass
    return output_file


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
    muscle_perm: str = "none",
    muscle_perturb: int = 0,
    cancel=None,
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
    output_stem = combined_alignment_output_stem(resolved, engine=engine)
    desired = _next_available_fasta(parent, output_stem)
    combined = parent / f".repeatviewer_{desired.stem}_input.fa"
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
            muscle_perm=muscle_perm,
            muscle_perturb=muscle_perturb,
            cancel=cancel,
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
