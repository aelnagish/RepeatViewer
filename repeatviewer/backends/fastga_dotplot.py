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
import sys
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
    """Pretty command for the process window (one option per line)."""
    from repeatviewer.external_process import command_text_multiline

    try:
        return command_text_multiline(cmd)
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
    from repeatviewer.conda_tools import which_tool

    for name in _FASTGA_NAMES:
        found = which_tool(name)
        if found:
            return found
    return None


def fastga_available() -> bool:
    return fastga_executable() is not None


_SUPPORTED_VERSION_HTML = "<b>Supported version: FastGA 1.5 or newer</b>"
_SUPPORTED_VERSION_PLAIN = "Supported version: FastGA 1.5 or newer"

_RECOMMENDED_ENV = "rv_env"
_RECOMMENDED_FASTGA = "fastga=1.5.1"
_CONDA_CREATE_CMD = (
    f"conda create -n {_RECOMMENDED_ENV} -c bioconda {_RECOMMENDED_FASTGA}"
)
_CONDA_INSTALL_CMD = (
    f"conda install -n {_RECOMMENDED_ENV} -c bioconda {_RECOMMENDED_FASTGA}"
)


def _fastga_conda_install_block(*, rich: bool) -> str:
    """Create-or-install into the shared recommended env (may already exist)."""
    if rich:
        return (
            "Open a Terminal window.<br><br>"
            "Install FastGA with conda into RepeatViewer:<br>"
            f"<b>If</b> Conda environment <code>{_RECOMMENDED_ENV}</code> exists:<br>"
            f"<code>{_CONDA_INSTALL_CMD}</code><br>"
            "<b>else:</b><br>"
            f"<code>{_CONDA_CREATE_CMD}</code>"
        )
    return (
        "Open a Terminal window.\n\n"
        "Install FastGA with conda into RepeatViewer:\n"
        f"If Conda environment {_RECOMMENDED_ENV} exists:\n"
        f"  {_CONDA_INSTALL_CMD}\n"
        "else:\n"
        f"  {_CONDA_CREATE_CMD}"
    )


def fastga_configuration_message(*, rich: bool = True) -> str:
    """How to install and point RepeatViewer at FastGA (platform-specific).

    Packed builds have no shell conda env; users must install FastGA separately
    and bind it under File → Setting → External tools.
    """
    version = _SUPPORTED_VERSION_HTML if rich else _SUPPORTED_VERSION_PLAIN
    install = _fastga_conda_install_block(rich=rich)
    env_name = _RECOMMENDED_ENV
    if sys.platform == "win32":
        if rich:
            return (
                "FastGA is not bundled with RepeatViewer and is not configured yet.<br><br>"
                f"{version}<br><br>"
                "On Windows, FastGA runs inside WSL2 (not as a native Windows "
                "binary).<br><br>"
                "1. Install WSL2 if needed (PowerShell as Administrator):"
                "<br><code>wsl --install</code><br><br>"
                f"2. Inside your WSL Linux distro:<br>{install}<br><br>"
                "3. <b>Then</b>, in RepeatViewer:<br>"
                "File → Setting → External tools (conda / WSL2) → enter "
                f"<code>{env_name}</code> (or the Linux path to that env), "
                "then click Test to confirm FastGA resolves."
            )
        return (
            "FastGA is not bundled with RepeatViewer and is not configured yet.\n\n"
            f"{version}\n\n"
            "On Windows, FastGA runs inside WSL2 (not as a native Windows binary).\n\n"
            "1. Install WSL2 if needed (PowerShell as Administrator):\n"
            "   wsl --install\n\n"
            f"2. Inside your WSL Linux distro:\n{install}\n\n"
            "3. Then, in RepeatViewer:\n"
            "   File → Setting → External tools (conda / WSL2) → enter "
            f"{env_name} (or the Linux path to that env), then Test."
        )
    if rich:
        return (
            "FastGA is not bundled with RepeatViewer and is not configured yet.<br><br>"
            f"{version}<br><br>"
            f"{install}<br><br>"
            "<b>Then</b>, in RepeatViewer:<br>"
            "File → Setting → External tools (conda) → enter "
            f"<code>{env_name}</code> (or the path to that env), "
            "then click Test to confirm FastGA resolves."
        )
    return (
        "FastGA is not bundled with RepeatViewer and is not configured yet.\n\n"
        f"{version}\n\n"
        f"{install}\n\n"
        "Then, in RepeatViewer:\n"
        "  File → Setting → External tools (conda) → enter "
        f"{env_name} (or the path to that env), then Test."
    )


def install_hint() -> str:
    """Plain-text install / bind instructions for the current platform."""
    return fastga_configuration_message(rich=False)


def show_fastga_not_configured_dialog(parent=None) -> None:
    """Remind the user how to install FastGA when it is missing."""
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QMessageBox

    box = QMessageBox(parent)
    box.setWindowTitle("FastGA not configured")
    box.setIcon(QMessageBox.Icon.Information)
    box.setTextFormat(Qt.TextFormat.RichText)
    # QMessageBox often renders setText() at bold weight; keep the body normal
    # so only explicit <b>…</b> (supported version) stands out — same as help_button.
    box.setText(
        f"<div style='font-weight:400'>{fastga_configuration_message(rich=True)}</div>"
    )
    box.setStandardButtons(QMessageBox.StandardButton.Ok)
    # Native macOS alerts ignore HTML weight / links.
    box.setOption(QMessageBox.Option.DontUseNativeDialog, True)
    box.exec()


# --------------------------------------------------------------- soft masking
#
# FastGA changed how soft masking is switched. Up to 1.4 it soft-masked
# lowercase bases by default and -N turned that off. From 1.5 (Dec 2025) it
# does not mask by default, -M turns masking on, and -N is rejected as an
# illegal option -- so a fixed flag breaks one version or the other. The usage
# text FastGA prints with no arguments tells us which it is.

MASK_STYLE_ON_FLAG = "M"    # 1.5+: -M turns soft masking on
MASK_STYLE_OFF_FLAG = "N"   # up to 1.4: -N turns soft masking off

_mask_style_cache: dict[str, str | None] = {}


def mask_style_from_usage(usage: str) -> str | None:
    """Which soft-mask switch a FastGA usage text documents (pure)."""
    head = usage.split("<source1", 1)[0]
    import re

    for match in re.finditer(r"\[-([A-Za-z]+)\]", head):
        letters = match.group(1)
        if "M" in letters:
            return MASK_STYLE_ON_FLAG
        if "N" in letters:
            return MASK_STYLE_OFF_FLAG
    if re.search(r"^\s*-M:", usage, re.MULTILINE):
        return MASK_STYLE_ON_FLAG
    if re.search(r"^\s*-N:", usage, re.MULTILINE):
        return MASK_STYLE_OFF_FLAG
    return None


def fastga_mask_style(exe: str | None = None) -> str | None:
    """Ask the installed FastGA once which soft-mask switch it takes."""
    exe = exe or fastga_executable()
    if not exe:
        return None
    if exe in _mask_style_cache:
        return _mask_style_cache[exe]
    style: str | None = None
    try:
        from repeatviewer.conda_tools import prepare_command
        from repeatviewer.external_process import windows_hidden_process_kwargs

        result = subprocess.run(
            prepare_command([exe]),
            capture_output=True,
            text=True,
            timeout=20,
            **windows_hidden_process_kwargs(),
        )
        style = mask_style_from_usage((result.stdout or "") + (result.stderr or ""))
    except (OSError, subprocess.SubprocessError):
        style = None
    _mask_style_cache[exe] = style
    return style


def soft_mask_flags(soft_mask: bool | None, style: str | None) -> list[str]:
    """Flags asking for soft masking on/off; None keeps FastGA's default."""
    if soft_mask is None or style is None:
        return []
    if style == MASK_STYLE_ON_FLAG:
        return ["-M"] if soft_mask else []
    if style == MASK_STYLE_OFF_FLAG:
        return [] if soft_mask else ["-N"]
    return []


# ------------------------------------------------------ against a reference
#
# "Selected files against a reference" runs FastGA once per query file with
# the same reference. FastGA would rebuild the reference index (.gix: about
# 14 GB per gigabase) on every run, so it is built once with GIXmake and the
# .gix is passed as the second source. With "Keep genome index" it is kept in
# a hidden .repeatviewer_fastga/ folder beside the reference (hidden folders
# are not listed in the Files panel) and reused while it is newer than the
# FASTA; otherwise it lives in scratch space and is removed afterwards.

FASTGA_MODE_SELF = "self"
FASTGA_MODE_PAIRS = "pairs"
FASTGA_MODE_REFERENCE = "reference"

REFERENCE_INDEX_DIRNAME = ".repeatviewer_fastga"


def gixmake_executable() -> str | None:
    from repeatviewer.conda_tools import which_tool

    return which_tool("GIXmake")


def _index_core_name(reference: Path) -> str:
    """GIX root name for a reference: its full file name, dots made safe."""
    return reference.name.replace(".", "_")


def cached_reference_index(reference: Path) -> Path | None:
    """The kept .gix for ``reference`` if it exists and is not older than it."""
    folder = reference.parent / REFERENCE_INDEX_DIRNAME
    core = _index_core_name(reference)
    gix = folder / f"{core}.gix"
    gdb = folder / f"{core}.1gdb"
    try:
        if gix.is_file() and gdb.is_file() and (
            gix.stat().st_mtime >= reference.stat().st_mtime
        ):
            return gix
    except OSError:
        return None
    return None


@contextmanager
def reference_index(
    reference: Path,
    *,
    threads: int = DEFAULT_THREADS,
    keep: bool = True,
    soft_mask: bool | None = None,
    temp_dir: Path | None = None,
    on_log=None,
    progress=None,
    is_cancelled=None,
):
    """Yield what to pass to FastGA as the reference: a .gix, or the FASTA.

    Falls back to the FASTA itself (FastGA then indexes it on every run) when
    GIXmake is missing or soft masking is requested -- a mask is baked into a
    GIX, and FastGA 1.5 only applies one it builds itself from the FASTA case.
    """
    reference = reference.expanduser().resolve()

    def _say(message: str) -> None:
        if progress is not None:
            progress(message)
        if on_log is not None:
            on_log(message)

    if soft_mask:
        _say("Soft masking is on: FastGA indexes the reference on every run.")
        yield reference
        return
    if keep:
        cached = cached_reference_index(reference)
        if cached is not None:
            _say(f"Reusing reference index {cached}")
            yield cached
            return
    gixmake = gixmake_executable()
    if not gixmake:
        _say("GIXmake not found: FastGA indexes the reference on every run.")
        yield reference
        return

    if keep:
        folder = reference.parent / REFERENCE_INDEX_DIRNAME
        try:
            folder.mkdir(parents=True, exist_ok=True)
            _remove_index(folder / _index_core_name(reference))  # stale copy
        except OSError:
            keep = False
    with _scratch_dir(reference.parent, temp_dir) as scratch:
        if not keep:
            folder = Path(tempfile.mkdtemp(prefix=".repeatviewer_fastga_ref_",
                                           dir=str(scratch)))
        core = folder / _index_core_name(reference)
        cmd = [gixmake, f"-T{max(1, int(threads))}", f"-P{scratch}",
               str(reference), str(core)]
        from repeatviewer.conda_tools import prepare_command
        from repeatviewer.external_process import windows_hidden_process_kwargs

        launch = prepare_command(cmd)
        _say(f"Indexing reference {reference.name} once for all runs…")
        if on_log is not None:
            from repeatviewer.external_process import PROCESS_LOG_COMMAND_HEADER

            on_log(PROCESS_LOG_COMMAND_HEADER)
            for line in _format_command_line(launch).splitlines():
                on_log(line)
        log_command("GIXmake", launch, inputs=[reference], cwd=reference.parent)
        try:
            proc = subprocess.Popen(
                launch,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                cwd=str(reference.parent),
                **windows_hidden_process_kwargs(),
            )
        except OSError as exc:
            raise FastgaError(f"Could not start GIXmake: {exc}") from exc
        stderr_lines, stderr_thread = _stream_stderr(proc, "GIXmake", on_line=on_log)
        while proc.poll() is None:
            if is_cancelled is not None and is_cancelled():
                proc.terminate()
                proc.wait()
                stderr_thread.join()
                _remove_index(core)
                raise FastgaError("Cancelled.")
            try:
                proc.wait(timeout=0.5)
            except subprocess.TimeoutExpired:
                continue
        stderr_thread.join()
        gix = core.with_name(core.name + ".gix")
        if proc.returncode != 0 or not gix.is_file():
            _remove_index(core)
            raise FastgaError(
                f"GIXmake failed (exit {proc.returncode}).\n"
                + "".join(stderr_lines).strip()[:600]
            )
        try:
            yield gix
        finally:
            if not keep:
                _remove_index(core)


def _remove_index(core: Path) -> None:
    """Delete a GDB/GIX pair and its hidden parts (.bps, .ktab.N, .post.N)."""
    folder = core.parent
    stem = core.name
    patterns = (f"{stem}.gix", f"{stem}.1gdb", f"{stem}.1ano", f".{stem}.*")
    for pattern in patterns:
        for path in folder.glob(pattern):
            try:
                path.unlink()
            except OSError:
                pass


def _same_file(a: Path, b: Path) -> bool:
    try:
        return a.expanduser().resolve() == b.expanduser().resolve()
    except OSError:
        return str(a) == str(b)


def reference_output_path(query: Path, reference: Path, fmt: str) -> Path:
    """``<query name>_FastGA_vs_<reference name><ext>`` next to the query.

    Rule 1 (keep the full input name, then a suffix). "_vs_<reference>" keeps
    it apart from the self-comparison result ``<query name>_FastGA<ext>``.
    Numbers ``_2``, ``_3``… when the path already exists.
    """
    from repeatviewer.backends.output_paths import next_numbered_path

    return next_numbered_path(
        query.parent
        / f"{query.name}_FastGA_vs_{reference.name}{output_suffix_for_format(fmt)}"
    )


def reference_combined_output_path(first_query: Path, reference: Path) -> Path:
    """One file for all queries (always PAF text: .1aln cannot be joined)."""
    from repeatviewer.backends.output_paths import next_numbered_path

    return next_numbered_path(
        first_query.parent
        / f"{first_query.name}_multi_FastGA_vs_{reference.name}.paf"
    )


def reference_queries(files, reference: Path | None) -> list[Path]:
    """Selected FASTA files minus the reference itself, order kept."""
    out: list[Path] = []
    for path in files:
        path = Path(path)
        if reference is not None and _same_file(path, reference):
            continue
        if any(_same_file(path, seen) for seen in out):
            continue
        out.append(path)
    return out


def pairwise_comparisons(files) -> list[tuple[Path, Path]]:
    """Unique unordered pairs among ``files`` (A vs B, A vs C, B vs C, …)."""
    import itertools

    paths = [Path(p) for p in files]
    return list(itertools.combinations(paths, 2))


def resolve_pairwise_comparisons(files, options) -> list[tuple[Path, Path]]:
    """All pairs, or only those listed in ``options.selected_pairs``."""
    pairs = pairwise_comparisons(files)
    selected = getattr(options, "selected_pairs", None)
    if not selected:
        return pairs
    wanted: set[tuple[str, str]] = set()
    for a, b in selected:
        wanted.add((str(Path(a)), str(Path(b))))
        wanted.add((str(Path(b)), str(Path(a))))
    return [(a, b) for a, b in pairs if (str(a), str(b)) in wanted]


_FASTGA_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".faa", ".fsa"}


def is_fastga_fasta_path(path: Path) -> bool:
    """True for uncompressed FASTA-like suffixes FastGA accepts as genomes."""
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return Path(name).suffix in _FASTGA_FASTA_SUFFIXES


def expand_fastga_input_paths(paths) -> list[Path]:
    """Resolve Define-inputs entries: FASTA files and folders of FASTA files."""
    out: list[Path] = []
    seen: set[Path] = set()
    for raw in paths or []:
        if raw is None or str(raw).strip() == "":
            continue
        try:
            path = Path(raw).expanduser().resolve()
        except OSError:
            continue
        candidates: list[Path] = []
        if path.is_file():
            if is_fastga_fasta_path(path):
                candidates.append(path)
        elif path.is_dir():
            try:
                children = sorted(path.iterdir(), key=lambda p: p.name.lower())
            except OSError:
                children = []
            for child in children:
                if child.is_file() and is_fastga_fasta_path(child):
                    candidates.append(child)
        for candidate in candidates:
            if candidate in seen:
                continue
            seen.add(candidate)
            out.append(candidate)
    return out


def fastga_mode_summary(
    files,
    mode: str,
    reference: Path | None = None,
    *,
    combined: bool = False,
    selected_pairs: list[tuple[str, str]] | None = None,
) -> str:
    """The grey line in the FastGA tab: what Run will do (pure)."""
    files = [Path(p) for p in files]
    n = len(files)
    if mode == FASTGA_MODE_REFERENCE:
        if reference is None:
            return "Choose the reference genome: click it in Files, or Browse."
        queries = reference_queries(files, reference)
        if not queries:
            if n:
                return (
                    f"Only the reference genome is selected. Select the query "
                    f"genomes to compare with {reference.name}."
                )
            return (
                f"No query genomes selected: select FASTA files to compare with "
                f"{reference.name}."
            )
        runs = len(queries)
        where = (
            "into one file"
            if combined and runs > 1
            else "one result file per query"
        )
        return (
            f"{runs} query genome{'s' if runs != 1 else ''} → {runs} run"
            f"{'s' if runs != 1 else ''} against {reference.name} ({where})."
        )
    if mode == FASTGA_MODE_PAIRS:
        if n < 2:
            return "Pairwise comparison needs two or more selected genomes."
        from types import SimpleNamespace

        pairs = resolve_pairwise_comparisons(
            files, SimpleNamespace(selected_pairs=selected_pairs)
        )
        total = n * (n - 1) // 2
        if not pairs:
            return f"{n} genomes selected → no comparisons checked."
        return (
            f"{n} genomes selected → {len(pairs)} comparison"
            f"{'s' if len(pairs) != 1 else ''}"
            + (f" (of {total} possible)" if len(pairs) != total else "")
            + ", written into one file."
        )
    if n == 0:
        return "Select one or more genomes: each is compared against itself."
    return (
        f"{n} genome{'s' if n != 1 else ''} selected → "
        f"{n} self-comparison{'s' if n != 1 else ''}."
    )


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
    soft_mask: bool | None,
    temp_dir: Path | None,
    mask_style: str | None = None,
) -> None:
    cmd.append(f"-T{max(1, int(threads))}")
    if temp_dir is not None:
        cmd.append(f"-P{temp_dir}")
    if symmetric:
        cmd.append("-S")
    if keep_index:
        cmd.append("-k")
    cmd.extend(soft_mask_flags(soft_mask, mask_style))
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
        # no_soft_mask asks for "off"; otherwise keep FastGA's own default.
        soft_mask=False if no_soft_mask else None,
        temp_dir=temp_dir,
        mask_style=fastga_mask_style(exe),
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
    soft_mask: bool | None = None,
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
            soft_mask=soft_mask,
            temp_dir=scratch,
            mask_style=fastga_mask_style(exe),
        )
        append_output_format_flag(cmd, fmt, out_path)

        cmd.append(str(source1))
        if source2 is not None:
            cmd.append(str(source2))

        from repeatviewer.conda_tools import prepare_command
        from repeatviewer.external_process import windows_hidden_process_kwargs

        launch_cmd = prepare_command(cmd)

        if on_log is not None:
            from repeatviewer.external_process import (
                PROCESS_LOG_COMMAND_HEADER,
                PROCESS_LOG_OUTPUT_HEADER,
            )

            on_log(PROCESS_LOG_COMMAND_HEADER)
            for line in _format_command_line(launch_cmd).splitlines():
                on_log(line)
            on_log(f"Working directory: {source1.parent}")
            on_log(f"Output: {out_path}")
            on_log("")
            on_log(PROCESS_LOG_OUTPUT_HEADER)

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
                launch_cmd,
                inputs=[p for p in (source1, source2) if p is not None],
                cwd=source1.parent,
            )
            proc = subprocess.Popen(
                launch_cmd,
                stdout=subprocess.PIPE if stream_to_file else subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                text=True,
                cwd=str(source1.parent),
                **windows_hidden_process_kwargs(),
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
            from repeatviewer.conda_tools import prepare_command
            from repeatviewer.external_process import windows_hidden_process_kwargs

            launch_cmd = prepare_command(cmd)
            log_command("FastGA", launch_cmd, inputs=[fasta], cwd=fasta.parent)
            proc = subprocess.Popen(
                launch_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                cwd=str(fasta.parent),
                **windows_hidden_process_kwargs(),
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


# ------------------------------------------------------------------ batches


def fastga_step_count(files, options) -> int:
    """Progress steps a batch reports (index build counts as one)."""
    mode = getattr(options, "mode", FASTGA_MODE_SELF)
    files = list(files)
    if mode == FASTGA_MODE_PAIRS and len(files) >= 2:
        return max(1, len(resolve_pairwise_comparisons(files, options)))
    if mode == FASTGA_MODE_REFERENCE:
        ref = Path(options.reference_fasta) if options.reference_fasta else None
        return 1 + len(reference_queries(files, ref))
    return max(1, len(files))


def run_fastga_batch(
    files,
    options,
    *,
    progress=None,
    on_log=None,
    is_cancelled=None,
    on_step=None,
) -> list[Path]:
    """Run FastGA over the selected files as ``options.mode`` says.

    Returns the files written. ``on_step(n)`` is called before each FastGA
    run (and the reference index build) so a progress bar can advance.
    """
    import itertools

    files = [Path(p) for p in files]
    mode = getattr(options, "mode", FASTGA_MODE_SELF)
    fmt = options.output_format
    temp_raw = str(getattr(options, "temp_dir", "") or "").strip()
    temp_dir = Path(temp_raw).expanduser() if temp_raw else None
    if temp_dir is not None:
        temp_dir.mkdir(parents=True, exist_ok=True)
    soft_mask = bool(getattr(options, "soft_mask", False))

    def _log(message: str) -> None:
        if on_log is not None:
            on_log(message)

    def _step(n: int) -> None:
        if on_step is not None:
            on_step(n)

    def _run(source1: Path, source2: Path | None, out: Path, run_fmt: str,
             *, keep_index: bool) -> Path:
        return run_fastga_alignment(
            source1,
            source2,
            out,
            fmt=run_fmt,
            threads=options.threads,
            min_identity=options.min_identity,
            min_length=options.min_length,
            seed_freq=options.seed_freq,
            chain_coverage=options.chain_coverage,
            chain_gap=options.chain_gap,
            symmetric=options.symmetric,
            keep_index=keep_index,
            soft_mask=soft_mask,
            temp_dir=temp_dir,
            progress=progress,
            on_log=on_log,
            is_cancelled=is_cancelled,
        )

    def _run_into_combined(jobs, out_path: Path, *, keep_index: bool,
                           first_step: int, label) -> Path:
        """Several runs streamed into one PAF/PSL text file.

        Callers must not pass ``FORMAT_ALN`` — binary ``.1aln`` cannot be joined.
        """
        if fmt == FORMAT_ALN:
            raise FastgaError(
                "FastGA .1aln output cannot be combined into one file. "
                "Use PAF/PSL for a single multi file, or leave .1aln and "
                "write one file per comparison."
            )
        chunks: list[str] = []
        for index, (source1, source2) in enumerate(jobs):
            _step(first_step + index)
            _log(label(index, source1, source2))
            with tempfile.NamedTemporaryFile(
                prefix=".repeatviewer_fastga_part_",
                suffix=output_suffix_for_format(fmt),
                dir=str(out_path.parent),
                delete=False,
            ) as tmp:
                tmp_path = Path(tmp.name)
            try:
                written = _run(source1, source2, tmp_path, fmt,
                               keep_index=keep_index)
                text = written.read_text(encoding="utf-8")
                if text and not text.endswith("\n"):
                    text += "\n"
                chunks.append(text)
            finally:
                try:
                    tmp_path.unlink(missing_ok=True)
                except OSError:
                    pass
        out_path.write_text("".join(chunks), encoding="utf-8")
        return out_path

    if mode == FASTGA_MODE_PAIRS:
        if len(files) < 2:
            raise FastgaError("Pairwise comparison needs two or more genomes.")
        pairs = resolve_pairwise_comparisons(files, options)
        if not pairs:
            raise FastgaError("No pairwise comparisons are checked.")
        from repeatviewer.backends.output_paths import next_numbered_path

        # Binary .1aln cannot be concatenated — one file per pair (like
        # reference mode without "Write all into one file").
        if fmt == FORMAT_ALN:
            created: list[Path] = []
            for index, (source1, source2) in enumerate(pairs):
                _step(1 + index)
                out_path = next_numbered_path(
                    source1.parent
                    / f"{source1.name}_FastGA_vs_{source2.name}.1aln"
                )
                _log(f"Pair {index + 1}/{len(pairs)}: {source1.name} vs {source2.name}")
                created.append(
                    _run(
                        source1,
                        source2,
                        out_path,
                        fmt,
                        keep_index=options.keep_index,
                    )
                )
            return created

        first = files[0]
        out_path = next_numbered_path(
            first.parent
            / f"{first.name}_multi_FastGA{output_suffix_for_format(fmt)}"
        )
        _run_into_combined(
            pairs, out_path, keep_index=options.keep_index, first_step=1,
            label=lambda i, a, b: f"Pair {i + 1}/{len(pairs)}: {a.name} vs {b.name}",
        )
        return [out_path]

    if mode == FASTGA_MODE_REFERENCE:
        if not options.reference_fasta:
            raise FastgaError("Choose a reference FASTA first (Align → FastGA tab).")
        reference = Path(options.reference_fasta).expanduser()
        if not reference.is_file():
            raise FastgaError(f"Reference not found: {reference}")
        queries = reference_queries(files, reference)
        if not queries:
            raise FastgaError(
                "Select the FASTA files to compare with the reference "
                "(the reference itself is skipped)."
            )
        _step(1)
        with reference_index(
            reference,
            threads=options.threads,
            keep=options.keep_index,
            soft_mask=soft_mask,
            temp_dir=temp_dir,
            on_log=on_log,
            progress=progress,
            is_cancelled=is_cancelled,
        ) as ref_source:
            # .1aln is binary — never join into one file; write one per query.
            combine = (
                bool(getattr(options, "reference_combined", False))
                and len(queries) > 1
                and fmt != FORMAT_ALN
            )
            if combine:
                out_path = reference_combined_output_path(queries[0], reference)
                _run_into_combined(
                    [(q, ref_source) for q in queries],
                    out_path,
                    keep_index=False,
                    first_step=2,
                    label=lambda i, a, _b: (
                        f"Run {i + 1}/{len(queries)}: {a.name} vs {reference.name}"
                    ),
                )
                return [out_path]
            created = []
            for index, query in enumerate(queries):
                _step(2 + index)
                if len(queries) > 1:
                    _log(f"Run {index + 1}/{len(queries)}: "
                         f"{query.name} vs {reference.name}")
                # Query indexes are temporary: "Keep genome index" is about the
                # reference, which is the one reused.
                created.append(_run(query, ref_source,
                                    reference_output_path(query, reference, fmt),
                                    fmt, keep_index=False))
            return created

    # Self: each file against itself.
    from repeatviewer.backends.output_paths import next_numbered_path

    created = []
    suffix = output_suffix_for_format(fmt)
    for index, source in enumerate(files, start=1):
        _step(index)
        out_path = next_numbered_path(
            source.parent / f"{source.name}_FastGA{suffix}"
        )
        if len(files) > 1:
            _log(f"Job {index}/{len(files)}: {out_path.name}")
        created.append(_run(source, None, out_path, fmt,
                            keep_index=options.keep_index))
    return created
