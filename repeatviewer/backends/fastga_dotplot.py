"""Self-comparison dotplots for large sequences, via FastGA.

The built-in k-mer engine in ``repeatviewer.backends.dotplot`` is instant on TE-sized
sequences but its match count grows faster than the sequence, so it refuses
anything much beyond a few megabases. FastGA (Gene Myers) is a whole-genome
aligner built for exactly this scale -- it compares two genomes, or a genome
against itself, using adaptive seeds and the DALIGNER wave-based local aligner.

``FastGA <fasta>`` performs a self-comparison "carefully avoiding self matches",
which is what we want: only the off-diagonal repeat structure is reported, and
the viewer draws the main diagonal itself.

Output is parsed from PAF, and returned as the same ``DotplotResult`` the k-mer
engine produces, so the canvas needs no changes.

FastGA is optional and detected at runtime, like BLAST and MAFFT::

    conda install -c bioconda fastga

Only single-sequence FASTA files are supported here. A multi-scaffold assembly
would need both axes laid out as concatenated scaffolds with boundary markers,
which the current viewer does not do.
"""

from __future__ import annotations

import shutil
import logging
import subprocess
import tempfile
import threading
import time
from contextlib import contextmanager
from repeatviewer.external_process import log_command, log_process_output

process_logger = logging.getLogger("repeatviewer.commands")


def _stream_stderr(
    proc: subprocess.Popen,
    tool: str,
    *,
    on_line=None,
) -> tuple[list[str], threading.Thread]:
    lines: list[str] = []

    def _read() -> None:
        if proc.stderr is None:
            return
        for line in iter(proc.stderr.readline, ""):
            lines.append(line)
            text = line.rstrip("\r\n")
            process_logger.info("%s stderr | %s", tool, text)
            if on_line is not None and text:
                try:
                    on_line(text)
                except Exception:  # noqa: BLE001 - never break the reader
                    pass

    thread = threading.Thread(target=_read, daemon=True)
    thread.start()
    return lines, thread


def _format_command_line(cmd: list[str]) -> str:
    import shlex

    try:
        return shlex.join(str(part) for part in cmd)
    except Exception:  # noqa: BLE001
        return " ".join(str(part) for part in cmd)
from pathlib import Path

from repeatviewer.backends.dotplot import (
    STRAND_FORWARD,
    STRAND_REVERSE,
    DotplotResult,
    DotplotSegment,
)

_FASTGA_NAMES = ("FastGA", "fastga")

DEFAULT_THREADS = 8
DEFAULT_MIN_IDENTITY = 0.7
DEFAULT_MIN_LENGTH = 100
# FastGA chaining defaults from the tool itself (-c 85, -s 1000).
DEFAULT_CHAIN_COVERAGE = 85
DEFAULT_CHAIN_GAP = 1000

# FastGA's -f: a seed occurring more than this many times is deemed
# "repetitive" and discarded before alignment. FastGA's own default is 10,
# which is tuned for synteny between genomes -- and which silently removes the
# entire repeat structure of a satellite array, where every monomer k-mer
# occurs hundreds of times. For self-comparison of repeats we default much
# higher; the cost is many more seed hits, so it is slower.
FASTGA_DEFAULT_SEED_FREQ = 10
DEFAULT_SEED_FREQ = 200
MAX_SEED_FREQ = 100_000

# Above this the built-in k-mer engine is impractical and FastGA takes over.
FASTGA_PREFERRED_ABOVE_BYTES = 2_000_000


class FastgaError(RuntimeError):
    """FastGA is missing, or failed to run."""


@contextmanager
def _scratch_dir(near: Path, explicit: Path | None = None):
    """A writable scratch directory for FastGA's intermediate .las files.

    FastGA defaults to $TMPDIR. Launched from a desktop app on macOS that is a
    sandboxed per-process directory, which is often small and is not the one
    you get in a shell -- producing "Cannot write overlap gather file". Placing
    scratch beside the input instead is predictable and has the same free space
    as the data. Removed afterwards, even on failure.
    """
    if explicit is not None:
        explicit.mkdir(parents=True, exist_ok=True)
        yield explicit
        return
    path = Path(tempfile.mkdtemp(prefix=".repeatviewer_fastga_", dir=str(near)))
    try:
        yield path
    finally:
        shutil.rmtree(path, ignore_errors=True)


def fastga_executable() -> str | None:
    for name in _FASTGA_NAMES:
        found = shutil.which(name)
        if found:
            return found
    return None


def fastga_available() -> bool:
    return fastga_executable() is not None


def install_hint() -> str:
    return "Install FastGA with:  conda install -c bioconda fastga"


def parse_paf_segments(
    lines,
    *,
    min_length: int = 0,
) -> tuple[list[DotplotSegment], int, str]:
    """Turn PAF records into dotplot segments.

    PAF columns used here are 0-based, half-open::

        0 qname  1 qlen  2 qstart  3 qend  4 strand
        5 tname  6 tlen  7 tstart  8 tend

    For a '+' record the alignment runs down-right, so it starts at
    (qstart, tstart). For a '-' record the query is reverse-complemented
    relative to the target, so it runs up-right: it starts at (qstart, tend)
    and descends to (qend, tstart). Getting that inversion wrong silently
    mirrors every inverted repeat, so it is covered by tests.

    Returns (segments, sequence_length, sequence_name).
    """
    segments: list[DotplotSegment] = []
    seq_length = 0
    seq_name = ""

    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        try:
            qname = parts[0]
            qlen = int(parts[1])
            qstart = int(parts[2])
            qend = int(parts[3])
            strand = parts[4]
            tstart = int(parts[7])
            tend = int(parts[8])
        except (ValueError, IndexError):
            continue

        if not seq_name:
            seq_name = qname
        seq_length = max(seq_length, qlen)

        span = min(qend - qstart, abs(tend - tstart))
        if span <= 0 or span < min_length:
            continue

        if strand == "-":
            segments.append(
                DotplotSegment(
                    x=qstart, y=tend, length=span, strand=STRAND_REVERSE
                )
            )
        else:
            segments.append(
                DotplotSegment(
                    x=qstart, y=tstart, length=span, strand=STRAND_FORWARD
                )
            )

    return segments, seq_length, seq_name


def _append_core_flags(
    cmd: list[str],
    *,
    threads: int,
    min_identity: float,
    min_length: int,
    seed_freq: int,
    chain_coverage: int,
    chain_gap: int,
    symmetric: bool,
    keep_index: bool,
    no_soft_mask: bool,
    temp_dir: Path | None,
) -> None:
    cmd.append(f"-T{max(1, int(threads))}")
    if temp_dir is not None:
        cmd.append(f"-P{temp_dir}")
    if symmetric:
        cmd.append("-S")
    if keep_index:
        cmd.append("-k")
    if no_soft_mask:
        cmd.append("-N")
    cmd.append(f"-i{float(min_identity):.3f}")
    cmd.append(f"-l{max(1, int(min_length))}")
    cmd.append(f"-f{max(1, int(seed_freq))}")
    cmd.append(f"-c{max(1, int(chain_coverage))}")
    cmd.append(f"-s{max(1, int(chain_gap))}")


def build_fastga_command(
    fasta: Path,
    *,
    threads: int = DEFAULT_THREADS,
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_length: int = DEFAULT_MIN_LENGTH,
    seed_freq: int = DEFAULT_SEED_FREQ,
    chain_coverage: int = DEFAULT_CHAIN_COVERAGE,
    chain_gap: int = DEFAULT_CHAIN_GAP,
    symmetric: bool = True,
    keep_index: bool = True,
    no_soft_mask: bool = False,
    temp_dir: Path | None = None,
) -> list[str]:
    """Command line for a FastGA self-comparison streaming PAF to stdout.

    ``symmetric`` (-S) is on by default: the docs recommend it "when
    understanding the repeat structure ... is of interest", which is the point
    of a self dotplot. ``keep_index`` (-k) keeps the .1gdb/.gidx beside the
    FASTA so reopening the same file skips index building.
    """
    exe = fastga_executable()
    if not exe:
        raise FastgaError(f"FastGA is not installed.\n\n{install_hint()}")

    cmd = [exe]
    _append_core_flags(
        cmd,
        threads=threads,
        min_identity=min_identity,
        min_length=min_length,
        seed_freq=seed_freq,
        chain_coverage=chain_coverage,
        chain_gap=chain_gap,
        symmetric=symmetric,
        keep_index=keep_index,
        no_soft_mask=no_soft_mask,
        temp_dir=temp_dir,
    )
    cmd.append("-paf")
    cmd.append(str(fasta))
    return cmd


FORMAT_PAF = "paf"
FORMAT_PAFX = "pafx"
FORMAT_PAFM = "pafm"
FORMAT_PAFS = "pafs"
FORMAT_PAFS_LONG = "pafS"
FORMAT_PSL = "psl"
FORMAT_ALN = "1aln"

# Streaming output flags FastGA accepts (excluding -v / -L).
FORMAT_STREAM_FLAGS: dict[str, str] = {
    FORMAT_PAF: "-paf",
    FORMAT_PAFX: "-pafx",
    FORMAT_PAFM: "-pafm",
    FORMAT_PAFS: "-pafs",
    FORMAT_PAFS_LONG: "-pafS",
    FORMAT_PSL: "-psl",
}

FORMAT_LABELS = (
    ("-paf: PAF text file", FORMAT_PAF),
    ("-pafx: PAF text with CIGAR (X mismatches)", FORMAT_PAFX),
    ("-pafm: PAF text with CIGAR (= matches)", FORMAT_PAFM),
    ("-pafs: PAF text with CS string (short)", FORMAT_PAFS),
    ("-pafS: PAF text with CS string (long)", FORMAT_PAFS_LONG),
    ("-psl: PSL text file", FORMAT_PSL),
    ("-1: .1aln binary (ONEcode)", FORMAT_ALN),
)


def output_suffix_for_format(fmt: str) -> str:
    if fmt == FORMAT_ALN:
        return ".1aln"
    if fmt == FORMAT_PSL:
        return ".psl"
    return ".paf"


def append_output_format_flag(cmd: list[str], fmt: str, out_path: Path) -> None:
    if fmt == FORMAT_ALN:
        cmd.append(f"-1:{out_path}")
        return
    flag = FORMAT_STREAM_FLAGS.get(fmt, "-paf")
    cmd.append(flag)


def run_fastga_alignment(
    source1: Path,
    source2: Path | None,
    out_path: Path,
    *,
    fmt: str = FORMAT_PAF,
    threads: int = DEFAULT_THREADS,
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_length: int = DEFAULT_MIN_LENGTH,
    seed_freq: int = DEFAULT_SEED_FREQ,
    chain_coverage: int = DEFAULT_CHAIN_COVERAGE,
    chain_gap: int = DEFAULT_CHAIN_GAP,
    symmetric: bool = True,
    keep_index: bool = True,
    no_soft_mask: bool = False,
    temp_dir: Path | None = None,
    progress=None,
    on_log=None,
    is_cancelled=None,
) -> Path:
    """Align two genomes, or one against itself, writing ``out_path``.

    With a single source FastGA does a self-comparison "carefully avoiding self
    matches", which is how repetitive regions are found. PAF and PSL are
    streamed to stdout; .1aln is binary and written by FastGA itself via -1.
    """
    exe = fastga_executable()
    if not exe:
        raise FastgaError(f"FastGA is not installed.\n\n{install_hint()}")

    source1 = source1.expanduser().resolve()
    if not source1.is_file():
        raise FastgaError(f"No such file: {source1}")
    if source2 is not None:
        source2 = source2.expanduser().resolve()
        if not source2.is_file():
            raise FastgaError(f"No such file: {source2}")

    out_path.parent.mkdir(parents=True, exist_ok=True)

    def _emit(message: str) -> None:
        if progress is not None:
            progress(message)
        if on_log is not None:
            on_log(message)

    with _scratch_dir(source1.parent, temp_dir) as scratch:
        cmd = [exe]
        _append_core_flags(
            cmd,
            threads=threads,
            min_identity=min_identity,
            min_length=min_length,
            seed_freq=seed_freq,
            chain_coverage=chain_coverage,
            chain_gap=chain_gap,
            symmetric=symmetric,
            keep_index=keep_index,
            no_soft_mask=no_soft_mask,
            temp_dir=scratch,
        )
        append_output_format_flag(cmd, fmt, out_path)

        cmd.append(str(source1))
        if source2 is not None:
            cmd.append(str(source2))

        if on_log is not None:
            on_log("Command:")
            on_log(_format_command_line(cmd))
            on_log(f"Working directory: {source1.parent}")
            on_log(f"Output: {out_path}")
            on_log("")
            on_log("—— program output ——")

        label = (
            source1.name if source2 is None
            else f"{source1.name} vs {source2.name}"
        )
        _emit(f"FastGA: {label} (indexing on first run may take a while)…")

        stream_to_file = fmt != FORMAT_ALN
        handle = out_path.open("w", encoding="utf-8") if stream_to_file else None
        try:
            log_command(
                "FastGA",
                cmd,
                inputs=[p for p in (source1, source2) if p is not None],
                cwd=source1.parent,
            )
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE if stream_to_file else subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                cwd=str(source1.parent),
            )
            process_started = time.monotonic()
            process_logger.info("START tool=FastGA pid=%s", proc.pid)
            if on_log is not None:
                on_log(f"START pid={proc.pid}")
            stderr_lines, stderr_thread = _stream_stderr(
                proc, "FastGA", on_line=on_log
            )
        except OSError as exc:
            if handle is not None:
                handle.close()
            raise FastgaError(f"Could not start FastGA: {exc}") from exc

        try:
            if stream_to_file and proc.stdout is not None and handle is not None:
                for count, line in enumerate(proc.stdout, 1):
                    handle.write(line)
                    if count % 50_000 == 0:
                        _emit(f"FastGA: {count:,} alignments…")
                    if is_cancelled is not None and is_cancelled():
                        proc.terminate()
                        raise FastgaError("Cancelled.")
            else:
                while proc.poll() is None:
                    if is_cancelled is not None and is_cancelled():
                        proc.terminate()
                        raise FastgaError("Cancelled.")
                    try:
                        proc.wait(timeout=0.5)
                    except subprocess.TimeoutExpired:
                        continue
        finally:
            if handle is not None:
                handle.close()
            proc.wait()
            stderr_thread.join()
            stderr = "".join(stderr_lines)

    elapsed = time.monotonic() - process_started
    process_logger.log(
        logging.INFO if proc.returncode == 0 else logging.ERROR,
        "FINISH tool=FastGA pid=%s exit=%s elapsed=%.3fs",
        proc.pid,
        proc.returncode,
        elapsed,
    )
    if on_log is not None:
        on_log(f"FINISH exit={proc.returncode} elapsed={elapsed:.3f}s")

    log_process_output("FastGA", output_path=out_path)

    if proc.returncode not in (0, None):
        raise FastgaError(
            f"FastGA failed (exit {proc.returncode}).\n{stderr.strip()[:600]}"
        )
    if not out_path.exists():
        raise FastgaError("FastGA produced no output file.")
    return out_path


def run_fastga_self_dotplot(
    fasta: Path,
    *,
    threads: int = DEFAULT_THREADS,
    min_identity: float = DEFAULT_MIN_IDENTITY,
    min_length: int = DEFAULT_MIN_LENGTH,
    seed_freq: int = DEFAULT_SEED_FREQ,
    chain_coverage: int = DEFAULT_CHAIN_COVERAGE,
    chain_gap: int = DEFAULT_CHAIN_GAP,
    symmetric: bool = True,
    keep_index: bool = True,
    no_soft_mask: bool = False,
    temp_dir: Path | None = None,
    progress=None,
    is_cancelled=None,
) -> DotplotResult:
    """Run a FastGA self-comparison and return it as a DotplotResult.

    ``progress`` is called with status strings; ``is_cancelled`` is polled and
    the subprocess terminated when it returns True. This runs for minutes on a
    large genome, so it must be driven from a worker thread.
    """
    fasta = fasta.expanduser().resolve()
    if not fasta.is_file():
        raise FastgaError(f"No such file: {fasta}")

    paf_lines: list[str] = []
    with _scratch_dir(fasta.parent, temp_dir) as scratch:
        cmd = build_fastga_command(
            fasta,
            threads=threads,
            min_identity=min_identity,
            min_length=min_length,
            seed_freq=seed_freq,
            chain_coverage=chain_coverage,
            chain_gap=chain_gap,
            symmetric=symmetric,
            keep_index=keep_index,
            no_soft_mask=no_soft_mask,
            temp_dir=scratch,
        )

        if progress is not None:
            progress("Running FastGA (indexing on first run may take a while)…")

        try:
            log_command("FastGA", cmd, inputs=[fasta], cwd=fasta.parent)
            proc = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=str(fasta.parent),
            )
            process_started = time.monotonic()
            process_logger.info("START tool=FastGA pid=%s", proc.pid)
            stderr_lines, stderr_thread = _stream_stderr(proc, "FastGA")
        except OSError as exc:
            raise FastgaError(f"Could not start FastGA: {exc}") from exc

        try:
            assert proc.stdout is not None
            for count, line in enumerate(proc.stdout, 1):
                paf_lines.append(line)
                if progress is not None and count % 50_000 == 0:
                    progress(f"FastGA: {count:,} alignments…")
                if is_cancelled is not None and is_cancelled():
                    proc.terminate()
                    raise FastgaError("Cancelled.")
        finally:
            proc.wait()
            stderr_thread.join()
            stderr = "".join(stderr_lines)

    process_logger.log(
        logging.INFO if proc.returncode == 0 else logging.ERROR,
        "FINISH tool=FastGA pid=%s exit=%s elapsed=%.3fs",
        proc.pid,
        proc.returncode,
        time.monotonic() - process_started,
    )

    if proc.returncode not in (0, None):
        raise FastgaError(
            f"FastGA failed (exit {proc.returncode}).\n{stderr.strip()[:600]}"
        )

    if progress is not None:
        progress("Parsing alignments…")

    segments, seq_length, seq_name = parse_paf_segments(
        paf_lines, min_length=min_length
    )

    if seq_length <= 0:
        raise FastgaError(
            "FastGA returned no alignments, so the sequence length is unknown.\n"
            "This usually means no repeats were found above the current "
            "identity and length thresholds."
        )

    note = ""
    if not segments:
        note = "No repeats found above the current thresholds."

    return DotplotResult(
        sequence_id=seq_name,
        sequence_length=seq_length,
        word_size=0,  # not applicable: FastGA uses adaptive seeds
        segments=segments,
        include_reverse=True,
        raw_match_count=len(segments),
        note=note,
    )
