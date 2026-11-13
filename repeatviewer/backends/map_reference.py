"""Map to Reference (Tools → Map to Reference).

Places sequences or reads on an (annotated) reference and reports what differs,
like Geneious Prime's "Map to Reference", with additions for repeat work:

* **Engines.** minimap2 for long reads or assemblies; HISAT2 for short DNA/RNA-seq.
  ``auto`` picks between them by input. minimap2 is bundled where platform
  builds exist; HISAT2 comes from conda or PATH. A missing engine the user
  chose is an error with an install hint; ``auto`` records any fallback in the
  report.
* **No pysam / samtools.** minimap2 writes PAF with the long ``cs`` tag and
  HISAT2 writes SAM; both are parsed here in pure Python.
* **Annotation lift-over.** Every reference feature (``.RV.gff3``) is checked in
  every mapped sequence: intact, changed, partial, lost or not covered, with its
  position on the mapped sequence. Mapped sequences and the consensus are
  written with the lifted annotations, so they open annotated in the viewer.

Results follow the RepeatViewer naming rule (full input name, then a suffix):
the selected aligner output beside the input
(``…_minimap.paf`` / ``…_minimap.sam`` / ``…_hisat.sam``), a contig-view
document ``…_minimap.rvmap`` / ``…_hisat.rvmap``, and a matching ``_details/``
folder with the report, consensus, mapped and unused sequences, and variant /
feature / read tables. A second run becomes ``…_minimap_2…`` / ``…_hisat_2…``;
nothing is overwritten. Older ``<alignment>.rv.json`` sidecars are still
recognised if present, but new runs do not write them.

Pure functions (parsers, projection, consensus, variants, lift-over) have no Qt
dependency and are tested in ``tests/test_map_reference.py``.
"""

from __future__ import annotations

import gzip
import hashlib
import html
import json
import re
import shutil
import statistics
import subprocess
import sys
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable, Iterable, Iterator

import numpy as np

from repeatviewer.job_options import (
    MAP_REPEAT_SEEDS_ARGS,
    MAP_REPEAT_SEEDS_DEFAULT,
    MAP_ENGINE_AUTO,
    MAP_ENGINE_HISAT2,
    MAP_ENGINE_MAFFT,
    MAP_ENGINE_MINIMAP2,
    MAP_HISAT2_DNA,
    MAP_HISAT2_RNA,
    MAP_HISAT2_READS_PAIRED,
    MAP_HISAT2_READS_SINGLE,
    MAP_MINIMAP2_OUT_PAF,
    MAP_MINIMAP2_OUT_SAM,
    MAP_MINIMAP2_PRESET_AUTO,
    MAP_MULTI_ALL,
    MAP_MULTI_BEST,
    MAP_MULTI_NONE,
    MAP_NO_COVERAGE_REFERENCE,
    MAP_NO_COVERAGE_SKIP,
    MAP_SENS_FAST,
    MAP_SENS_HIGH,
    MapOptions,
)
from repeatviewer.tool_output_paths import kept_input_name

RVMAP_SUFFIX = ".rvmap"
RVMAP_FORMAT = "repeatviewer-map"
RVMAP_VERSION = 1
RUN_JSON_SUFFIX = ".rv.json"          # ``reads.fa_minimap.sam`` → ``….sam.rv.json``
RUN_STAMP_FORMAT = "%Y-%m-%d_%H%M"

# Auto engine: HISAT2 for many short reads, otherwise minimap2.
SHORT_READ_LENGTH = 400
SHORT_READ_MIN_QUERIES = 500
# Per-base detail (consensus, variants, lift-over) up to this reference length;
# longer contigs get binned coverage and the reads only.
DETAIL_LIMIT = 5_000_000
COVERAGE_BINS = 4_000
MAX_VARIANT_ROWS = 20_000
MAX_FEATURE_COLUMNS = 200
MAX_READ_TABLE = 100_000
JSON_BASE_BUDGET = 40_000_000    # bases of drawn placements stored in one .rvmap
JSON_MAX_FEATURES = 20_000       # reference features per sequence stored in the .rvmap

ENGINE_LABELS = {
    MAP_ENGINE_AUTO: "Auto (choose by input)",
    MAP_ENGINE_MINIMAP2: "minimap2 (long reads, assemblies)",
    MAP_ENGINE_HISAT2: "HISAT2 (short DNA / RNA-seq)",
}
ENGINE_NAMES = {
    MAP_ENGINE_MINIMAP2: "minimap2",
    MAP_ENGINE_HISAT2: "HISAT2",
}
_EXECUTABLES = {
    MAP_ENGINE_MINIMAP2: ("minimap2",),
    MAP_ENGINE_HISAT2: ("hisat2",),
}
_CONDA_PACKAGE = {MAP_ENGINE_MINIMAP2: "minimap2", MAP_ENGINE_HISAT2: "hisat2"}
# Sibling-file tag before ``.rvmap`` (minimap2 → minimap, HISAT2 → hisat).
OUTPUT_ENGINE_TAG = {
    MAP_ENGINE_MINIMAP2: "minimap",
    MAP_ENGINE_HISAT2: "hisat",
}

STATUS_INTACT = "intact"
STATUS_CHANGED = "changed"
STATUS_PARTIAL = "partial"
STATUS_LOST = "lost"
STATUS_NOT_COVERED = "not covered"

_FASTQ_SUFFIXES = {".fq", ".fastq"}
_IUPAC = {
    frozenset("A"): "A", frozenset("C"): "C", frozenset("G"): "G", frozenset("T"): "T",
    frozenset("AG"): "R", frozenset("CT"): "Y", frozenset("CG"): "S", frozenset("AT"): "W",
    frozenset("GT"): "K", frozenset("AC"): "M", frozenset("CGT"): "B", frozenset("AGT"): "D",
    frozenset("ACT"): "H", frozenset("ACG"): "V", frozenset("ACGT"): "N",
}
BASES = "ACGTN-"          # count columns: A C G T other gap
_BASE_INDEX = np.full(256, 4, dtype=np.int8)
for _i, _b in enumerate("ACGT"):
    _BASE_INDEX[ord(_b)] = _i
    _BASE_INDEX[ord(_b.lower())] = _i
_BASE_INDEX[ord("U")] = 3
_BASE_INDEX[ord("u")] = 3
_BASE_INDEX[ord("-")] = 5


class MapReferenceError(RuntimeError):
    """Mapping cannot run; the message is shown to the user as is."""


# --------------------------------------------------------------------------
# Sequence input
# --------------------------------------------------------------------------


def _open_text(path: Path):
    path = Path(path)
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open(encoding="utf-8", errors="replace")


def is_fastq(path: Path) -> bool:
    suffixes = [s.lower() for s in Path(path).suffixes]
    if suffixes and suffixes[-1] == ".gz":
        suffixes = suffixes[:-1]
    return bool(suffixes) and suffixes[-1] in _FASTQ_SUFFIXES


def iter_sequences(path: Path) -> Iterator[tuple[str, str]]:
    """(name, sequence) from FASTA or FASTQ (optionally gzipped). Name = first word."""
    with _open_text(path) as handle:
        if is_fastq(path):
            while True:
                header = handle.readline()
                if not header:
                    return
                if not header.strip():
                    continue
                seq = handle.readline().strip()
                handle.readline()
                handle.readline()
                yield header[1:].split()[0] if header[1:].split() else "", seq
            return
        name, chunks = None, []
        for line in handle:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                words = line[1:].split()
                name, chunks = (words[0] if words else ""), []
            elif name is not None:
                chunks.append(line.strip())
        if name is not None:
            yield name, "".join(chunks)


@dataclass
class InputInfo:
    n: int = 0
    total: int = 0
    lengths_sample: list[int] = field(default_factory=list)
    longest: int = 0
    names: list[str] = field(default_factory=list)
    complete: bool = True            # False when only the first records were read

    @property
    def median(self) -> float:
        return float(statistics.median(self.lengths_sample)) if self.lengths_sample else 0.0


def scan_inputs(files: Iterable[Path], *, keep_names: int = 2_000_000,
                sample: int = 20_000, max_records: int | None = None) -> InputInfo:
    """Count and sample the input. With ``max_records`` only that many records
    are read (large read sets); ``complete`` then says whether all were seen."""
    info = InputInfo()
    for path in files:
        for name, seq in iter_sequences(Path(path)):
            if max_records is not None and info.n >= max_records:
                info.complete = False
                return info
            n = len(seq)
            info.n += 1
            info.total += n
            info.longest = max(info.longest, n)
            if len(info.lengths_sample) < sample:
                info.lengths_sample.append(n)
            if len(info.names) < keep_names:
                info.names.append(name)
    return info


def reference_lengths(path: Path, duplicates: list[str] | None = None) -> dict[str, int]:
    """Record lengths, from a current ``.fai`` index when there is one, otherwise
    by reading the FASTA line by line (no whole-chromosome strings).

    Names that occur more than once are appended to ``duplicates`` when given:
    the mappers report hits by name only, so such a reference cannot be used.
    """
    path = Path(path)
    seen: set[str] = set()

    def _note(name: str) -> None:
        if name in seen and duplicates is not None and name not in duplicates:
            duplicates.append(name)
        seen.add(name)

    fai = path.with_name(path.name + ".fai")
    try:
        if fai.is_file() and fai.stat().st_mtime >= path.stat().st_mtime:
            out = {}
            for line in fai.read_text(encoding="utf-8").splitlines():
                parts = line.split("\t")
                if len(parts) >= 2:
                    _note(parts[0])
                    out[parts[0]] = int(parts[1])
            if out:
                return out
    except (OSError, ValueError):
        pass
    out: dict[str, int] = {}
    name = None
    with _open_text(path) as handle:
        for line in handle:
            if line.startswith(">"):
                words = line[1:].split()
                name = words[0] if words else ""
                _note(name)
                out[name] = 0
            elif name is not None:
                out[name] += len(line.strip())
    return out


def load_reference_records(path: Path, wanted: set[str] | None = None,
                           max_length: int | None = None) -> dict[str, str]:
    out: dict[str, str] = {}
    for name, seq in iter_sequences(path):
        if wanted is not None and name not in wanted:
            continue
        if max_length is not None and len(seq) > max_length:
            continue
        out[name] = seq.upper()
    return out


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTUNacgtunRYKMSWBDHVrykmswbdhv",
                                       "TGCAANtgcaanYRMKSWVHDByrmkswvhdb"))[::-1]


# --------------------------------------------------------------------------
# Engines
# --------------------------------------------------------------------------


def engine_executable(engine: str, *, system_path_only: bool = False) -> str | None:
    """minimap2 / HISAT2 executable.

    ``system_path_only``: PATH only. Otherwise prefer the bundled binary
    (minimap2), then the conda env from Settings, then PATH.
    """
    if not system_path_only and engine == MAP_ENGINE_MINIMAP2:
        from repeatviewer.backends.minimap2_tools import bundled_minimap2_executable

        bundled = bundled_minimap2_executable()
        if bundled:
            return bundled
    for name in _EXECUTABLES.get(engine, ()):
        if system_path_only:
            found = shutil.which(name)
        else:
            from repeatviewer.conda_tools import which_tool

            found = which_tool(name) or shutil.which(name)
        if found:
            return found
    return None


def hisat2_build_executable(hisat2: str) -> str | None:
    """``hisat2-build`` next to ``hisat2`` (same install), else PATH."""
    sibling = Path(hisat2).with_name("hisat2-build")
    if sibling.is_file():
        return str(sibling)
    for name in ("hisat2-build", "hisat2-build.exe", "hisat2-build.bat"):
        sib = Path(hisat2).parent / name
        if sib.is_file():
            return str(sib)
    from repeatviewer.conda_tools import which_tool

    return which_tool("hisat2-build") or shutil.which("hisat2-build")


# Bioconda current stable; Mapping calls the HISAT2 2.2 CLI (-p, -k, --dta, …).
HISAT2_SUPPORTED_VERSION = "2.2"
HISAT2_RECOMMENDED_VERSION = "2.2.3"
_HISAT2_SUPPORTED_VERSION_HTML = (
    f"<b>Supported version: HISAT2 {HISAT2_SUPPORTED_VERSION} or newer</b>"
)
_HISAT2_SUPPORTED_VERSION_PLAIN = (
    f"Supported version: HISAT2 {HISAT2_SUPPORTED_VERSION} or newer"
)

_HISAT2_RECOMMENDED_ENV = "rv_env"
_HISAT2_CONDA_PKG = f"hisat2={HISAT2_RECOMMENDED_VERSION}"
_HISAT2_CONDA_INSTALL_CMD = (
    f"conda install -n {_HISAT2_RECOMMENDED_ENV} -c conda-forge -c bioconda "
    f"{_HISAT2_CONDA_PKG}"
)
_HISAT2_CONDA_CREATE_CMD = (
    f"conda create -n {_HISAT2_RECOMMENDED_ENV} -c conda-forge -c bioconda "
    f"{_HISAT2_CONDA_PKG}"
)


def _hisat2_conda_install_block(*, rich: bool) -> str:
    """Create-or-install into the shared recommended env (same pattern as FastGA)."""
    if rich:
        return (
            "Open a Terminal window.<br><br>"
            "Install HISAT2 with conda into RepeatViewer:<br>"
            f"<b>If</b> Conda environment <code>{_HISAT2_RECOMMENDED_ENV}</code> "
            "exists:<br>"
            f"<code>{_HISAT2_CONDA_INSTALL_CMD}</code><br>"
            "<b>else:</b><br>"
            f"<code>{_HISAT2_CONDA_CREATE_CMD}</code>"
        )
    return (
        "Open a Terminal window.\n\n"
        "Install HISAT2 with conda into RepeatViewer:\n"
        f"If Conda environment {_HISAT2_RECOMMENDED_ENV} exists:\n"
        f"  {_HISAT2_CONDA_INSTALL_CMD}\n"
        "else:\n"
        f"  {_HISAT2_CONDA_CREATE_CMD}"
    )


def hisat2_configuration_message(*, rich: bool = True) -> str:
    """How to install and point RepeatViewer at HISAT2 (FastGA reminder style)."""
    version = (
        _HISAT2_SUPPORTED_VERSION_HTML if rich else _HISAT2_SUPPORTED_VERSION_PLAIN
    )
    install = _hisat2_conda_install_block(rich=rich)
    env_name = _HISAT2_RECOMMENDED_ENV
    if sys.platform == "win32":
        if rich:
            return (
                "HISAT2 is not bundled with RepeatViewer and is not configured "
                "yet.<br><br>"
                f"{version}<br><br>"
                "On Windows, HISAT2 typically runs inside WSL2 (not as a native "
                "Windows binary).<br><br>"
                "1. Install WSL2 if needed (PowerShell as Administrator):"
                "<br><code>wsl --install</code><br><br>"
                f"2. Inside your WSL Linux distro:<br>{install}<br><br>"
                "3. <b>Then</b>, in RepeatViewer:<br>"
                "File → Setting → External tools (conda / WSL2) → enter "
                f"<code>{env_name}</code> (or the Linux path to that env), "
                "then click Test to confirm HISAT2 resolves."
            )
        return (
            "HISAT2 is not bundled with RepeatViewer and is not configured yet.\n\n"
            f"{version}\n\n"
            "On Windows, HISAT2 typically runs inside WSL2 (not as a native "
            "Windows binary).\n\n"
            "1. Install WSL2 if needed (PowerShell as Administrator):\n"
            "   wsl --install\n\n"
            f"2. Inside your WSL Linux distro:\n{install}\n\n"
            "3. Then, in RepeatViewer:\n"
            "   File → Setting → External tools (conda / WSL2) → enter "
            f"{env_name} (or the Linux path to that env), then Test."
        )
    if rich:
        return (
            "HISAT2 is not bundled with RepeatViewer and is not configured "
            "yet.<br><br>"
            f"{version}<br><br>"
            f"{install}<br><br>"
            "<b>Then</b>, in RepeatViewer:<br>"
            "File → Setting → External tools (conda) → enter "
            f"<code>{env_name}</code> (or the path to that env), "
            "then click Test to confirm HISAT2 resolves."
        )
    return (
        "HISAT2 is not bundled with RepeatViewer and is not configured yet.\n\n"
        f"{version}\n\n"
        f"{install}\n\n"
        "Then, in RepeatViewer:\n"
        "  File → Setting → External tools (conda) → enter "
        f"{env_name} (or the path to that env), then Test."
    )


def show_hisat2_not_configured_dialog(parent=None) -> None:
    """Remind the user how to install HISAT2 when it is missing (FastGA style)."""
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QMessageBox

    box = QMessageBox(parent)
    box.setWindowTitle("HISAT2 not configured")
    box.setIcon(QMessageBox.Icon.Information)
    box.setTextFormat(Qt.TextFormat.RichText)
    # Keep body normal weight so only explicit <b>…</b> stands out — same as FastGA.
    box.setText(
        f"<div style='font-weight:400'>{hisat2_configuration_message(rich=True)}</div>"
    )
    box.setStandardButtons(QMessageBox.StandardButton.Ok)
    box.setOption(QMessageBox.Option.DontUseNativeDialog, True)
    box.exec()


def engine_install_hint(engine: str) -> str:
    if engine == MAP_ENGINE_MINIMAP2:
        from repeatviewer.backends.minimap2_tools import (
            MINIMAP2_VERSION,
            bundled_minimap2_executable,
        )

        if bundled_minimap2_executable():
            return (
                f"minimap2 {MINIMAP2_VERSION} is bundled with RepeatViewer. "
                "Leave \"Use Minimap2 from the system PATH or defined Conda environment\" "
                "off to use it, or tick that box to use your own install."
            )
        return (
            f"minimap2 {MINIMAP2_VERSION} should be bundled with RepeatViewer but was not found. "
            "Reinstall the app, or install it yourself:\n"
            "  conda install -c conda-forge -c bioconda minimap2"
        )
    if engine == MAP_ENGINE_HISAT2:
        return hisat2_configuration_message(rich=False)
    pkg = _CONDA_PACKAGE.get(engine, engine)
    return ("Install it into the conda env set in File → Settings → External tools:\n"
            f"  conda install -c conda-forge -c bioconda {pkg}")


def engine_available(engine: str, *, system_path_only: bool = False) -> bool:
    if engine == MAP_ENGINE_MAFFT:
        # MAFFT is no longer a Mapping engine.
        return False
    return engine_executable(engine, system_path_only=system_path_only) is not None


def choose_engine(ref_lengths: dict[str, int], queries: InputInfo) -> str:
    """The engine ``auto`` uses for this input (minimap2 or HISAT2)."""
    del ref_lengths  # kept for call-site compatibility
    if (
        queries.median
        and queries.median < SHORT_READ_LENGTH
        and queries.n > SHORT_READ_MIN_QUERIES
    ):
        return MAP_ENGINE_HISAT2
    return MAP_ENGINE_MINIMAP2


def auto_minimap2_preset(queries: InputInfo) -> str:
    """``-x`` value for Map Options → Preset = Auto.

    * median length &lt; 400 bp → ``sr`` (short reads)
    * ≤200 sequences and median ≥100 kb → ``asm20`` (assemblies / BACs)
    * otherwise → ``map-ont`` (tolerant long-read default)
    """
    if queries.median and queries.median < SHORT_READ_LENGTH:
        return "sr"
    if queries.n <= 200 and queries.median >= 100_000:
        return "asm20"
    return "map-ont"


def preset_kmer_window(preset: str) -> tuple[int, int]:
    """Concrete ``-k`` / ``-w`` for a resolved minimap2 ``-x`` preset."""
    from repeatviewer.job_options import MAP_MINIMAP2_PRESET_KW

    return MAP_MINIMAP2_PRESET_KW.get(preset, (15, 10))


def minimap2_arguments(opts: MapOptions, preset: str) -> list[str]:
    """CLI args for minimap2 (after the executable; before reference and query)."""
    out = getattr(opts, "minimap2_output", MAP_MINIMAP2_OUT_PAF) or MAP_MINIMAP2_OUT_PAF
    args: list[str] = ["-x", preset, "-t", str(max(1, opts.threads))]
    if out == MAP_MINIMAP2_OUT_SAM:
        args = ["-a", "--MD", *args]
        if getattr(opts, "minimap2_include_cs", False):
            args.append("--cs=long")
    else:
        # PAF (default). Contig view needs long cs; CIGAR is optional (-c).
        if getattr(opts, "minimap2_include_cigar", True):
            args.append("-c")
        if getattr(opts, "minimap2_include_cs", True):
            args.append("--cs=long")
        args.append("--paf-no-hit")
    k = int(getattr(opts, "minimap2_kmer", 0) or 0)
    w = int(getattr(opts, "minimap2_window", 0) or 0)
    # Auto always uses the resolved preset's -k/-w; explicit presets use the spinbox values.
    if getattr(opts, "minimap2_preset", "") == MAP_MINIMAP2_PRESET_AUTO or k <= 0 or w <= 0:
        dk, dw = preset_kmer_window(preset)
        if getattr(opts, "minimap2_preset", "") == MAP_MINIMAP2_PRESET_AUTO or k <= 0:
            k = dk
        if getattr(opts, "minimap2_preset", "") == MAP_MINIMAP2_PRESET_AUTO or w <= 0:
            w = dw
    args += ["-k", str(min(28, k)), "-w", str(w)]
    args += MAP_REPEAT_SEEDS_ARGS.get(
        getattr(opts, "minimap2_repeat_seeds", MAP_REPEAT_SEEDS_DEFAULT) or
        MAP_REPEAT_SEEDS_DEFAULT, [])
    if opts.multiple_hits == MAP_MULTI_ALL:
        args += ["-N", "50", "-p", "0.95", "--secondary=yes"]
    else:
        args += ["--secondary=no"]
    return args


def hisat2_arguments(opts: MapOptions) -> list[str]:
    """CLI args for HISAT2 (after the executable; before -x / -U / -S)."""
    args = ["-p", str(max(1, opts.threads))]
    library = getattr(opts, "hisat2_library", MAP_HISAT2_RNA) or MAP_HISAT2_RNA
    if library == MAP_HISAT2_DNA:
        args.append("--no-spliced-alignment")
    else:
        strand = getattr(opts, "hisat2_strandness", "") or ""
        if strand:
            args += ["--rna-strandness", strand]
        if getattr(opts, "hisat2_dta", False):
            args.append("--dta")
        min_i = int(getattr(opts, "hisat2_min_intron", 20) or 20)
        max_i = int(getattr(opts, "hisat2_max_intron", 500_000) or 500_000)
        if min_i > 0:
            args += ["--min-intronlen", str(min_i)]
        if max_i > 0:
            args += ["--max-intronlen", str(max(min_i, max_i))]
    if getattr(opts, "hisat2_no_softclip", False):
        args.append("--no-softclip")
    if opts.multiple_hits == MAP_MULTI_ALL:
        args += ["-k", "20", "--secondary"]
    elif opts.multiple_hits == MAP_MULTI_BEST:
        args += ["-k", "1"]
    return args


def mafft_arguments(opts: MapOptions) -> list[str]:
    args = ["--nuc", "--preservecase", "--thread", str(max(1, opts.threads))]
    args += ["--6merpair"] if opts.sensitivity == MAP_SENS_FAST else ["--multipair"]
    args += ["--adjustdirectionaccurately"]
    return args


def tool_version(executable: str) -> str:
    try:
        from repeatviewer.conda_tools import prepare_command

        result = subprocess.run(prepare_command([executable, "--version"]),
                                capture_output=True, text=True, timeout=20)
    except (OSError, subprocess.SubprocessError):
        return ""
    text = (result.stdout or "") + "\n" + (result.stderr or "")
    m = re.search(r"version\s*v?(\d+\.\d+(?:\.\d+)*)", text) or re.search(
        r"\bv?(\d+\.\d+(?:[.-]r?\d+)*)\b", text)
    return m.group(1) if m else ""


# --------------------------------------------------------------------------
# Mapped sequences and the parsers that make them
# --------------------------------------------------------------------------


@dataclass
class MappedRead:
    """One placement of a sequence on a reference record.

    ``row`` holds the sequence's bases on reference positions ``start`` to
    ``end`` (0-based, end exclusive), already in reference orientation;
    ``-`` marks a reference base deleted in the sequence. ``insertions`` lists
    (p, bases) inserted before reference position p. ``ref_row`` is the
    reference over the same span when the aligner reports it ('' otherwise).
    """

    name: str
    ref_id: str
    start: int
    end: int
    strand: str
    row: str
    insertions: list[tuple[int, str]] = field(default_factory=list)
    ref_row: str = ""
    mapq: int = -1
    kind: str = "primary"            # primary / secondary / supplementary
    query_length: int = 0
    clipped: int = 0
    clip_left: int = 0               # unaligned bases before the placed part (reference orientation)
    matches: int = 0
    mismatches: int = 0
    deleted: int = 0
    inserted: int = 0

    @property
    def aligned_bp(self) -> int:
        """Sequence bases placed on the reference (matches, mismatches, insertions)."""
        return len(self.row) - self.deleted + self.inserted

    @property
    def identity(self) -> float:
        cols = len(self.row) + self.inserted
        return self.matches / cols if cols else 0.0

    @property
    def display_name(self) -> str:
        return f"{self.name} (reversed)" if self.strand == "-" else self.name


def compute_stats(read: MappedRead, ref_row: str | None = None) -> MappedRead:
    """Fill matches / mismatches / deleted / inserted from ``ref_row``."""
    ref = (ref_row if ref_row is not None else read.ref_row).upper()
    row = read.row.upper()
    matches = mismatches = deleted = 0
    if ref and len(ref) == len(row):
        a = np.frombuffer(row.encode("ascii", "replace"), dtype=np.uint8)
        b = np.frombuffer(ref.encode("ascii", "replace"), dtype=np.uint8)
        gap = a == ord("-")
        deleted = int(gap.sum())
        same = (a == b) & ~gap
        matches = int(same.sum())
        mismatches = int((~gap & ~same).sum())
    else:
        deleted = row.count("-")
        mismatches = 0
        matches = len(row) - deleted
    read.matches, read.mismatches, read.deleted = matches, mismatches, deleted
    read.inserted = sum(len(s) for _p, s in read.insertions)
    return read


def project_gapped_pair(ref_gapped: str, query_gapped: str) -> tuple | None:
    """Project a pairwise gapped alignment onto reference coordinates.

    Returns (start, end, row, ref_row, insertions, clipped, clip_left) or None when the
    query shares no column with the reference. Query bases before the first or
    after the last shared column count as clipped, not as insertions.
    """
    if len(ref_gapped) != len(query_gapped):
        raise ValueError("alignment rows differ in length")
    shared = [i for i, (r, q) in enumerate(zip(ref_gapped, query_gapped))
              if r != "-" and q != "-"]
    if not shared:
        return None
    first, last = shared[0], shared[-1]
    ref_pos = sum(1 for c in ref_gapped[:first] if c != "-")
    clip_left = sum(1 for c in query_gapped[:first] if c != "-")
    clipped = clip_left + sum(1 for c in query_gapped[last + 1:] if c != "-")
    start = ref_pos
    row: list[str] = []
    ref_row: list[str] = []
    insertions: list[tuple[int, str]] = []
    pending: list[str] = []
    for i in range(first, last + 1):
        r, q = ref_gapped[i], query_gapped[i]
        if r == "-":
            if q != "-":
                pending.append(q)
            continue
        if pending:
            insertions.append((ref_pos, "".join(pending).upper()))
            pending = []
        row.append(q.upper() if q != "-" else "-")
        ref_row.append(r.upper())
        ref_pos += 1
    return start, ref_pos, "".join(row), "".join(ref_row), insertions, clipped, clip_left


_CS_TOKEN = re.compile(r"(=[A-Za-z]+|:\d+|\*[A-Za-z][A-Za-z]|\+[A-Za-z]+|-[A-Za-z]+|"
                       r"~[A-Za-z]{2}\d+[A-Za-z]{2})")


def parse_cs(cs: str, tstart: int) -> tuple[str, str, list[tuple[int, str]]]:
    """(row, ref_row, insertions) from a minimap2 ``cs`` string (long form).

    The short form ``:n`` gives no bases; those columns get ``?`` in both rows
    and are filled from the reference later."""
    row: list[str] = []
    ref: list[str] = []
    ins: list[tuple[int, str]] = []
    pos = tstart
    for token in _CS_TOKEN.findall(cs):
        op, body = token[0], token[1:]
        if op == "=":
            seq = body.upper()
            row.append(seq)
            ref.append(seq)
            pos += len(seq)
        elif op == ":":
            n = int(body)
            row.append("?" * n)
            ref.append("?" * n)
            pos += n
        elif op == "*":
            ref.append(body[0].upper())
            row.append(body[1].upper())
            pos += 1
        elif op == "+":
            ins.append((pos, body.upper()))
        elif op == "-":
            row.append("-" * len(body))
            ref.append(body.upper())
            pos += len(body)
        elif op == "~":
            n = int(body[2:-2])
            row.append("-" * n)
            ref.append("N" * n)
            pos += n
    return "".join(row), "".join(ref), ins


def parse_paf_line(line: str) -> tuple[str, int, MappedRead | None]:
    """(query name, query length, placement or None for an unmapped line)."""
    f = line.rstrip("\n").split("\t")
    qname, qlen = f[0], int(f[1])
    if len(f) < 12 or f[4] == "*":
        return qname, qlen, None
    qstart, qend, strand = int(f[2]), int(f[3]), f[4]
    tname, tstart = f[5], int(f[7])
    mapq = int(f[11])
    tags = dict(t.split(":", 2)[0::2] for t in f[12:] if t.count(":") >= 2)
    cs = tags.get("cs", "")
    if not cs:
        raise MapReferenceError("minimap2 output has no cs tag (run with -c --cs=long).")
    row, ref_row, ins = parse_cs(cs, tstart)
    kind = {"P": "primary", "S": "secondary"}.get(tags.get("tp", "P"), "primary")
    read = MappedRead(name=qname, ref_id=tname, start=tstart, end=tstart + len(row),
                      strand=strand, row=row, insertions=ins, ref_row=ref_row, mapq=mapq,
                      kind=kind, query_length=qlen, clipped=qlen - (qend - qstart),
                      clip_left=qstart if strand == "+" else qlen - qend)
    return qname, qlen, compute_stats(read)


_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_MD = re.compile(r"(\d+)|(\^[A-Za-z]+)|([A-Za-z])")


def _ref_row_from_md(row: str, md: str) -> str:
    ref = list(row)
    j = 0
    for num, dele, mis in _MD.findall(md):
        if num:
            j += int(num)
        elif dele:
            for c in dele[1:]:
                if j < len(ref):
                    ref[j] = c.upper()
                j += 1
        elif mis:
            if j < len(ref):
                ref[j] = mis.upper()
            j += 1
    return "".join(ref)


def parse_sam_line(line: str) -> tuple[str, int, MappedRead | None]:
    """(query name, query length, placement or None when unmapped). Header lines
    must be skipped by the caller."""
    f = line.rstrip("\n").split("\t")
    qname, flag = f[0], int(f[1])
    seq = f[9] if f[9] != "*" else ""
    cigar = f[5]
    hard = sum(int(n) for n, op in _CIGAR.findall(cigar) if op == "H")
    qlen = len(seq) + hard
    if flag & 4 or f[2] == "*" or cigar == "*":
        return qname, qlen, None
    if not seq:
        return qname, qlen, None          # secondary without bases: nothing to draw
    pos = int(f[3]) - 1
    start = pos
    row: list[str] = []
    ins: list[tuple[int, str]] = []
    q = 0
    clipped = 0
    clip_left = 0
    seen_aligned = False
    for n_text, op in _CIGAR.findall(cigar):
        if op in "M=XID":
            seen_aligned = True
        elif op in "SH" and not seen_aligned:
            clip_left += int(n_text)
        n = int(n_text)
        if op in "M=X":
            row.append(seq[q:q + n].upper())
            q += n
            pos += n
        elif op == "I":
            ins.append((pos, seq[q:q + n].upper()))
            q += n
        elif op in "DN":
            row.append("-" * n)
            pos += n
        elif op == "S":
            clipped += n
            q += n
        elif op == "H":
            clipped += n
    row_text = "".join(row)
    tags = {t[:2]: t[5:] for t in f[11:] if len(t) > 5}
    ref_row = _ref_row_from_md(row_text, tags["MD"]) if "MD" in tags else ""
    kind = "secondary" if flag & 256 else ("supplementary" if flag & 2048 else "primary")
    read = MappedRead(name=qname, ref_id=f[2], start=start, end=pos,
                      strand="-" if flag & 16 else "+", row=row_text, insertions=ins,
                      ref_row=ref_row, mapq=int(f[4]), kind=kind, query_length=qlen,
                      clipped=clipped, clip_left=clip_left)
    return qname, qlen, compute_stats(read)


def parse_mafft_alignment(pairs: list[tuple[str, str]], ref_name: str, ref_id: str,
                          safe_to_name: dict[str, str]) -> tuple[list[MappedRead], list[str]]:
    """Placements from MAFFT ``--addfragments`` output (reference row + fragments).

    ``_R_`` in front of a name means MAFFT reversed that sequence."""
    rows = {name: seq for name, seq in pairs}
    ref_gapped = rows.get(ref_name)
    if ref_gapped is None:
        raise MapReferenceError("MAFFT output has no reference row.")
    reads: list[MappedRead] = []
    unmapped: list[str] = []
    for name, gapped in pairs:
        if name == ref_name:
            continue
        strand = "+"
        safe = name
        if safe.startswith("_R_"):
            safe, strand = safe[3:], "-"
        original = safe_to_name.get(safe, safe)
        qlen = sum(1 for c in gapped if c != "-")
        projected = project_gapped_pair(ref_gapped, gapped)
        if projected is None:
            unmapped.append(original)
            continue
        start, end, row, ref_row, ins, clipped, clip_left = projected
        reads.append(compute_stats(MappedRead(
            name=original, ref_id=ref_id, start=start, end=end, strand=strand, row=row,
            insertions=ins, ref_row=ref_row, query_length=qlen, clipped=clipped,
            clip_left=clip_left)))
    return reads, unmapped


def fill_unknown_bases(read: MappedRead, reference: str) -> None:
    """Replace ``?`` (cs short form) with reference bases and recount."""
    if "?" not in read.row and read.ref_row and "?" not in read.ref_row:
        return
    ref = reference[read.start:read.end]
    read.row = "".join(r if c == "?" else c for c, r in zip(read.row, ref))
    read.ref_row = ref
    compute_stats(read)


# --------------------------------------------------------------------------
# Filtering and multiple placements
# --------------------------------------------------------------------------


def placement_verdict(read: MappedRead, opts: MapOptions, engine: str) -> str | None:
    """None when the placement is kept; otherwise why it is not (the Advanced-tab
    limits and the multiple-placement choice). Secondary placements are dropped
    silently unless 'every equally good place' is chosen."""
    if read.kind == "secondary" and opts.multiple_hits != MAP_MULTI_ALL:
        return ""
    ambiguous_mapq = 1 if engine == MAP_ENGINE_HISAT2 else 0
    if (opts.multiple_hits == MAP_MULTI_NONE and read.kind == "primary"
            and 0 <= read.mapq <= ambiguous_mapq):
        return "ambiguous: fits equally well in several places"
    if read.aligned_bp < opts.min_aligned_bp:
        return f"aligned part shorter than {opts.min_aligned_bp} bp"
    if read.identity < opts.min_identity:
        return f"identity {read.identity:.0%} below {opts.min_identity:.0%}"
    if read.mapq >= 0 and read.mapq < opts.min_mapq:
        return f"mapping quality {read.mapq} below {opts.min_mapq}"
    return None


@dataclass
class FilterResult:
    kept: list[MappedRead]
    unused: dict[str, str]           # name → reason


def filter_placements(reads: list[MappedRead], all_names: Iterable[str], opts: MapOptions,
                      engine: str, unmapped: Iterable[str] = ()) -> FilterResult:
    """Apply :func:`placement_verdict` to a list (the run streams instead)."""
    unused: dict[str, str] = {}
    for name in unmapped:
        unused.setdefault(name, "no alignment found")
    kept: list[MappedRead] = []
    for read in reads:
        reason = placement_verdict(read, opts, engine)
        if reason is None:
            kept.append(read)
        elif reason:
            unused.setdefault(read.name, reason)
    placed = {r.name for r in kept}
    for name in all_names:
        if name not in placed:
            unused.setdefault(name, "no alignment found")
    for name in placed:
        unused.pop(name, None)
    return FilterResult(kept, unused)


# --------------------------------------------------------------------------
# Per-base counts, consensus, variants
# --------------------------------------------------------------------------


def row_indices(row: str) -> np.ndarray:
    """Column index in BASES (A C G T other gap) for every character of ``row``."""
    return _BASE_INDEX[np.frombuffer(row.encode("ascii", "replace"), dtype=np.uint8)]


def add_to_counts(counts: np.ndarray, read: MappedRead) -> None:
    idx = row_indices(read.row)
    positions = np.arange(read.start, read.start + len(idx))
    ok = positions < counts.shape[0]
    np.add.at(counts, (positions[ok], idx[ok].astype(np.intp)), 1)


_CODE_BITS = {"A": 1, "C": 2, "G": 4, "T": 8}


def _code_table() -> list[str]:
    table = []
    for bits in range(64):
        if bits & 16:                                   # N among the chosen: unknown
            table.append("N")
            continue
        bases = frozenset(b for b, v in _CODE_BITS.items() if bits & v)
        if not bases:
            table.append("" if bits & 32 else "N")      # only the gap: position deleted
        else:
            table.append(_IUPAC[bases])
    return table


_CODES = _code_table()


def consensus_calls(counts: np.ndarray, threshold: float = 0.0) -> list[str]:
    """Consensus character per reference position ('' = deleted, None-like '?' =
    no coverage).

    ``threshold`` 0: the most common state (ties give an IUPAC code). Otherwise
    the fewest states that together reach ``threshold`` of the depth, as an
    IUPAC code (Geneious style). A gap wins only on its own."""
    depth = counts.sum(axis=1)
    if threshold <= 0:
        top = counts.max(axis=1, keepdims=True)
        chosen = (counts == top) & (counts > 0)
    else:
        order = np.argsort(-counts, axis=1, kind="stable")
        ranked = np.take_along_axis(counts, order, axis=1)
        cum = np.cumsum(ranked, axis=1)
        need = (cum < (threshold * depth)[:, None] - 1e-9).sum(axis=1) + 1
        rank = np.argsort(order, axis=1)
        chosen = (rank < need[:, None]) & (counts > 0)
    weights = np.array([1, 2, 4, 8, 16, 32])
    bits = (chosen * weights).sum(axis=1)
    calls = [_CODES[b] for b in bits.tolist()]
    for i in np.flatnonzero(depth == 0).tolist():
        calls[i] = "?"
    return calls


def build_consensus(reference: str, counts: np.ndarray, insertions: dict[int, Counter],
                    *, threshold: float = 0.0, no_coverage: str = "N"
                    ) -> tuple[str, np.ndarray]:
    """(consensus, ref_to_cons) where ref_to_cons[i] is the 0-based consensus
    position of reference position i, or -1 when it was deleted / skipped."""
    calls = consensus_calls(counts, threshold)
    depth = counts.sum(axis=1)
    out: list[str] = []
    ref_to_cons = np.full(len(calls), -1, dtype=np.int64)
    n = 0
    need_ins = max(0.5, threshold)
    for i, call in enumerate(calls):
        if i in insertions:
            counter = insertions[i]
            spanning = min(int(depth[i - 1]) if i > 0 else 0, int(depth[i]))
            n_ins = sum(counter.values())
            if spanning and n_ins / spanning >= need_ins:
                seq = counter.most_common(1)[0][0]
                out.append(seq)
                n += len(seq)
        if call == "?":
            if no_coverage == MAP_NO_COVERAGE_SKIP:
                continue
            call = reference[i] if no_coverage == MAP_NO_COVERAGE_REFERENCE and i < len(reference) \
                else "N"
        if call:
            ref_to_cons[i] = n
            out.append(call)
            n += 1
    return "".join(out), ref_to_cons


def methylation_context(reference: str, pos: int) -> str:
    """'CG', 'CHG' or '' for the cytosine context at ``pos`` (either strand)."""
    ref = reference
    b = ref[pos] if 0 <= pos < len(ref) else ""
    if b == "C":
        if ref[pos + 1:pos + 2] == "G":
            return "CG"
        if ref[pos + 2:pos + 3] == "G":
            return "CHG"
    elif b == "G":
        if pos >= 1 and ref[pos - 1] == "C":
            return "CG"
        if pos >= 2 and ref[pos - 2] == "C":
            return "CHG"
    return ""


@dataclass
class Variant:
    pos: int            # 1-based reference position
    ref: str
    alt: str            # base, '-' (deletion) or '+SEQ' (insertion before pos)
    count: int
    depth: int
    features: str = ""
    context: str = ""

    @property
    def freq(self) -> float:
        return self.count / self.depth if self.depth else 0.0


def call_variants(reference: str, counts: np.ndarray, insertions: dict[int, Counter],
                  features: list, *, min_freq: float = 0.0, min_count: int = 1,
                  limit: int = MAX_VARIANT_ROWS) -> list[Variant]:
    """Every state that differs from the reference, with its count and depth."""
    L = min(len(reference), counts.shape[0])
    if L == 0:
        return []
    ref_idx = row_indices(reference[:L]).astype(np.intp)
    depth = counts[:L].sum(axis=1)
    out: list[Variant] = []
    for j, base in enumerate(BASES):
        col = counts[:L, j]
        mask = (col >= max(1, min_count)) & (ref_idx != j) & (depth > 0)
        if min_freq > 0:
            mask &= col >= min_freq * depth
        for i in np.flatnonzero(mask).tolist():
            out.append(Variant(i + 1, reference[i], base, int(col[i]), int(depth[i])))
    for p, counter in insertions.items():
        if p <= 0 or p >= L:
            continue
        spanning = int(min(depth[p - 1], depth[p]))
        for seq, c in counter.items():
            if c >= max(1, min_count) and spanning and c / spanning >= min_freq:
                out.append(Variant(p + 1, reference[p], "+" + seq, c, spanning))
    out.sort(key=lambda v: (v.pos, v.alt))
    if len(out) > limit:
        out = sorted(out, key=lambda v: -v.count)[:limit]
        out.sort(key=lambda v: (v.pos, v.alt))
    starts = [f.start for f in features]
    for v in out:
        v.context = methylation_context(reference, v.pos - 1)
        names = [f.name for f in features[:_bisect_right(starts, v.pos)] if f.end >= v.pos]
        v.features = ", ".join(names[:6])
    return out


def _bisect_right(values: list[int], x: int) -> int:
    import bisect

    return bisect.bisect_right(values, x)


# --------------------------------------------------------------------------
# Annotation lift-over
# --------------------------------------------------------------------------


@dataclass
class LiftedFeature:
    index: int          # position in the contig's feature list
    status: str
    differences: int    # mismatched + deleted + inserted bases inside the feature
    identity: float
    start: int          # 1-based on the mapped sequence (0 when not placed)
    end: int


def read_offsets(read: MappedRead) -> tuple[np.ndarray, np.ndarray]:
    """(bases_before, present): for each reference column of the read, how many
    sequence bases (including insertions) come before it, and whether the
    column holds a base."""
    n = len(read.row)
    present = (np.frombuffer(read.row.encode("ascii", "replace"), dtype=np.uint8)
               != ord("-")).astype(np.int64)
    ins_len = np.zeros(n + 1, dtype=np.int64)
    for p, seq in read.insertions:
        k = p - read.start
        if 0 <= k <= n:
            ins_len[k] += len(seq)
    before = np.concatenate(([0], np.cumsum(present)))[:n] + np.cumsum(ins_len)[:n]
    return before, present


@dataclass
class FeatureIndex:
    """Features sorted by start, for finding those overlapping a placement fast."""

    features: list
    starts: list[int] = field(default_factory=list)
    longest: int = 0

    def __post_init__(self) -> None:
        self.starts = [f.start for f in self.features]
        self.longest = max((f.end - f.start + 1 for f in self.features), default=0)

    def overlapping(self, start0: int, end0: int) -> range:
        """Indices of features that may overlap 0-based [start0, end0)."""
        import bisect

        lo = bisect.bisect_left(self.starts, start0 + 1 - self.longest)
        hi = bisect.bisect_right(self.starts, end0)
        return range(max(0, lo), hi)


def lift_features(read: MappedRead, features, reference: str = "") -> list[LiftedFeature]:
    """Status of each reference feature in this placement (features not
    overlapping the placement are left out). Coordinates are on the whole
    mapped sequence in reference orientation (clipped ends included).
    ``features`` is a list sorted by start, or a :class:`FeatureIndex`."""
    index = features if isinstance(features, FeatureIndex) else FeatureIndex(list(features or []))
    if not index.features:
        return []
    ref = read.ref_row or reference[read.start:read.end]
    before, present = read_offsets(read)
    out: list[LiftedFeature] = []
    ins_positions = [(p, len(s)) for p, s in read.insertions]
    for idx in index.overlapping(read.start, read.end):
        feat = index.features[idx]
        fs, fe = feat.start - 1, feat.end - 1
        lo, hi = max(fs, read.start), min(fe, read.end - 1)
        if lo > hi:
            continue
        flen = fe - fs + 1
        a, b = lo - read.start, hi - read.start
        seg = read.row[a:b + 1]
        rseg = ref[a:b + 1] if len(ref) == len(read.row) else ""
        dels = seg.count("-")
        if rseg:
            mism = sum(1 for x, y in zip(seg, rseg) if x != "-" and x != y)
        else:
            mism = 0
        ins_bp = sum(n for p, n in ins_positions if lo < p <= hi)
        covered = hi - lo + 1
        diffs = mism + dels + ins_bp
        matches = covered - dels - mism
        identity = matches / (covered + ins_bp) if covered + ins_bp else 0.0
        if covered < flen:
            status = STATUS_PARTIAL
        elif dels >= 0.5 * flen:
            status = STATUS_LOST
        elif diffs == 0:
            status = STATUS_INTACT
        else:
            status = STATUS_CHANGED
        s_on = int(before[a]) + 1
        e_on = int(before[b]) + int(present[b])
        if e_on < s_on:
            s_on = e_on = 0
        else:
            s_on += read.clip_left
            e_on += read.clip_left
        out.append(LiftedFeature(idx, status, diffs, round(identity, 4), s_on, e_on))
    return out


# --------------------------------------------------------------------------
# Reference annotations
# --------------------------------------------------------------------------


def load_reference_features(reference: Path, lengths: dict[str, int] | None = None
                            ) -> dict[str, list]:
    """Features of every reference record (``.RV.gff3`` plus the colour sidecar),
    sorted by start. Empty lists when the reference has no annotations."""
    from repeatviewer.backends.sequence_features import (
        _merge_feature_colors_from_json,
        load_all_features_by_id,
        load_features_by_id_from_json,
        resolve_features_gff_path,
    )

    ref = Path(reference).expanduser().resolve()
    lengths = lengths if lengths is not None else reference_lengths(ref)
    ids = list(lengths)
    try:
        gff = resolve_features_gff_path(ref)
        by_id = (load_all_features_by_id(gff, ids, lengths_by_id=lengths) if gff is not None
                 else {rid: [] for rid in ids})
        for rid, js in load_features_by_id_from_json(ref, ids).items():
            if js:
                by_id[rid] = (_merge_feature_colors_from_json(by_id[rid], js)
                              if by_id.get(rid) else js)
    except Exception:  # noqa: BLE001 - annotations are optional
        by_id = {rid: [] for rid in ids}
    return {rid: sorted(feats, key=lambda f: (f.start, f.end)) for rid, feats in by_id.items()}


# --------------------------------------------------------------------------
# Streaming accumulation
# --------------------------------------------------------------------------


@dataclass
class ReadRecord:
    """A kept placement without its bases (for the tables and the FASTA)."""

    name: str
    ref_id: str
    start: int
    end: int
    strand: str
    identity: float
    mismatches: int
    deleted: int
    inserted: int
    aligned_bp: int
    mapq: int
    kind: str
    query_length: int
    clipped: int
    clip_left: int
    lifted: list[LiftedFeature] = field(default_factory=list)
    placement: int = 1              # 1 for the first placement of this name, 2, 3, …

    @property
    def display_name(self) -> str:
        base = f"{self.name} (reversed)" if self.strand == "-" else self.name
        return base if self.placement == 1 else f"{base} #{self.placement}"


class ContigAccumulator:
    def __init__(self, ref_id: str, length: int, detailed: bool, features: list,
                 max_display: int, seed: int = 12345) -> None:
        import random

        self.ref_id = ref_id
        self.length = length
        self.detailed = detailed
        self.features = features
        self.max_display = max_display
        self.rng = random.Random(seed)
        self.counts = np.zeros((length, 6), dtype=np.int32) if detailed else None
        self.insertions: dict[int, Counter] = {}
        self.bin_size = 1 if detailed else max(1, -(-length // COVERAGE_BINS))
        self.bins = None if detailed else np.zeros(-(-length // self.bin_size), dtype=np.int64)
        self.n = 0
        self.display: list[tuple[MappedRead, list[LiftedFeature], int]] = []
        self.identities: list[float] = []

    def add(self, read: MappedRead, lifted: list[LiftedFeature], placement: int) -> None:
        self.n += 1
        if len(self.identities) < 200_000:
            self.identities.append(read.identity)
        if self.detailed:
            add_to_counts(self.counts, read)
            for p, seq in read.insertions:
                self.insertions.setdefault(p, Counter())[seq] += 1
        else:
            idx = row_indices(read.row)
            pos = np.arange(read.start, read.start + len(idx))[idx != 5] // self.bin_size
            pos = pos[pos < len(self.bins)]
            self.bins += np.bincount(pos, minlength=len(self.bins))[:len(self.bins)]
        item = (read, lifted, placement)
        if len(self.display) < self.max_display:
            self.display.append(item)
        else:
            j = self.rng.randrange(self.n)
            if j < self.max_display:
                self.display[j] = item

    def coverage(self) -> np.ndarray:
        """Bases per reference position (detailed) or mean per bin."""
        if self.detailed:
            return self.counts[:, :5].sum(axis=1)
        return self.bins / float(self.bin_size)


def rle(values: np.ndarray) -> list[list[int]]:
    """[[value, run], …] for integer arrays (coverage in the .rvmap file)."""
    v = np.asarray(values)
    if v.size == 0:
        return []
    change = np.flatnonzero(np.diff(v)) + 1
    starts = np.concatenate(([0], change))
    runs = np.diff(np.concatenate((starts, [v.size])))
    return [[int(v[s]), int(r)] for s, r in zip(starts.tolist(), runs.tolist())]


def unrle(pairs: list) -> np.ndarray:
    if not pairs:
        return np.zeros(0, dtype=np.int64)
    vals = np.array([p[0] for p in pairs], dtype=np.float64)
    runs = np.array([p[1] for p in pairs], dtype=np.int64)
    return np.repeat(vals, runs)


# --------------------------------------------------------------------------
# Output layout
# --------------------------------------------------------------------------


def run_stamp(now: datetime | None = None) -> str:
    return (now or datetime.now()).strftime(RUN_STAMP_FORMAT)


def output_base_name(inputs: list[Path]) -> str:
    first = kept_input_name(inputs[0])
    return first if len(inputs) == 1 else f"{first}_multi"


@dataclass
class OutputLayout:
    parent: Path
    base: str
    details: Path

    @property
    def rvmap(self) -> Path:
        return self.parent / f"{self.base}{RVMAP_SUFFIX}"

    def path(self, suffix: str) -> Path:
        return self.details / f"{self.base}{suffix}"

    def alignment(self, suffix: str) -> Path:
        """Selected aligner format (``.paf`` / ``.sam``) next to the run JSON."""
        ext = suffix if suffix.startswith(".") else f".{suffix}"
        return self.parent / f"{self.base}{ext}"


def run_json_path(alignment: Path) -> Path:
    """Run record beside an alignment: ``path.sam`` → ``path.sam.rv.json``."""
    return Path(str(Path(alignment)) + RUN_JSON_SUFFIX)


_RUN_JSON_TAILS = (
    ".sam.gz.rv.json", ".sam.rv.json", ".bam.rv.json", ".cram.rv.json", ".paf.rv.json",
)
_ALIGNMENT_TAILS = (".sam.gz", ".cram", ".bam", ".sam", ".paf")


def alignment_run_stem(path: Path) -> str:
    """Filename without ``.sam`` / ``.bam`` / ``.cram`` / ``.paf`` / ``.sam.gz``."""
    name = Path(path).name
    lower = name.lower()
    for tail in _ALIGNMENT_TAILS:
        if lower.endswith(tail):
            return name[: -len(tail)]
    if lower.endswith(RUN_JSON_SUFFIX):
        return alignment_run_stem(Path(name[: -len(RUN_JSON_SUFFIX)]))
    return Path(path).stem


def find_rvmap_for_alignment(path: Path) -> Path | None:
    """Contig document ``{stem}.rvmap`` beside an alignment (SAM / BAM / PAF)."""
    path = Path(path)
    cand = path.parent / f"{alignment_run_stem(path)}{RVMAP_SUFFIX}"
    return cand if cand.is_file() else None


def find_run_json(path: Path) -> Path | None:
    """Legacy ``*.rv.json`` beside an alignment (older Mapping runs).

    Looks in the same directory for ``{path}.rv.json``, then
    ``{stem}.sam.rv.json``, ``.bam.rv.json``, ``.cram.rv.json``, ``.paf.rv.json``,
    and ``.sam.gz.rv.json``. Prefer :func:`find_mapping_document` for new code.
    """
    path = Path(path)
    parent = path.parent
    stem = alignment_run_stem(path)
    ordered = [Path(str(path) + RUN_JSON_SUFFIX)]
    ordered.extend(parent / f"{stem}{tail}" for tail in _RUN_JSON_TAILS)
    seen: set[str] = set()
    for cand in ordered:
        key = str(cand)
        if key in seen:
            continue
        seen.add(key)
        if cand.is_file():
            return cand
    return None


def find_mapping_document(path: Path) -> Path | None:
    """Mapping snapshot for an alignment: ``.rvmap`` first, else legacy ``.rv.json``."""
    return find_rvmap_for_alignment(path) or find_run_json(path)


def is_alignment_path(path: Path) -> bool:
    """True for SAM / BAM / CRAM / PAF, including ``.sam.gz``."""
    name = Path(path).name.lower()
    return any(name.endswith(tail) for tail in _ALIGNMENT_TAILS)


def is_run_json_file(path: Path) -> bool:
    """True for a mapping run record (``*.rv.json``), hidden beside its alignment."""
    return Path(path).name.lower().endswith(RUN_JSON_SUFFIX)


def alignment_run_json_companions(path: Path) -> list[Path]:
    """The run record that belongs with this alignment, if it is still on disk."""
    path = Path(path)
    if not is_alignment_path(path):
        return []
    found = find_run_json(path)
    return [found] if found is not None else []


def alignment_run_json_renames(old: Path, new: Path) -> list[tuple[Path, Path]]:
    """Keep the run record named ``<new alignment>.rv.json``.

    Renaming, moving, or saving SAM as BAM / CRAM moves the one record with the
    alignment. A ``reads.sam.rv.json`` becomes ``reads.bam.rv.json``.
    """
    old, new = Path(old), Path(new)
    if not is_alignment_path(old):
        return []
    found = find_run_json(old)
    if found is None:
        return []
    dest = run_json_path(new)
    try:
        same = found.resolve() == dest.resolve()
    except OSError:
        same = found == dest
    return [] if same else [(found, dest)]


def alignment_for_run_json(path: Path) -> Path | None:
    """Alignment file that a ``*.rv.json`` record belongs to, if it is still on disk."""
    path = Path(path)
    name = path.name
    if not name.lower().endswith(RUN_JSON_SUFFIX):
        return path if path.is_file() else None
    stem = name[: -len(RUN_JSON_SUFFIX)]
    direct = path.parent / stem
    if direct.is_file():
        return direct
    base = alignment_run_stem(Path(stem))
    for tail in _ALIGNMENT_TAILS:
        cand = path.parent / f"{base}{tail}"
        if cand.is_file():
            return cand
    return None


def _base_outputs_exist(parent: Path, base: str) -> bool:
    """True when a previous run already claimed this base name."""
    if (parent / f"{base}{RVMAP_SUFFIX}").exists():
        return True
    for ext in (".paf", ".sam"):
        if (parent / f"{base}{ext}").exists():
            return True
        if (parent / f"{base}{ext}{RUN_JSON_SUFFIX}").exists():
            return True
    return False


def plan_outputs(inputs: list[Path], *, engine: str | None = None) -> OutputLayout:
    """``<input>_minimap`` / ``_hisat`` (or ``_2``, ``_3`` …) and details folder."""
    parent = Path(inputs[0]).expanduser().resolve().parent
    engine_key = engine or MAP_ENGINE_MINIMAP2
    engine_tag = OUTPUT_ENGINE_TAG.get(engine_key, "minimap")
    stem = f"{output_base_name(inputs)}_{engine_tag}"
    n = 1
    while True:
        base = stem if n == 1 else f"{stem}_{n}"
        details = parent / f"{base}_details"
        if not _base_outputs_exist(parent, base):
            try:
                details.mkdir(exist_ok=False)
                return OutputLayout(parent, base, details)
            except FileExistsError:
                pass
        n += 1


# --------------------------------------------------------------------------
# Running the engines
# --------------------------------------------------------------------------


def _check_cancel(cancel) -> None:
    if cancel is not None and getattr(cancel, "requested", False):
        from repeatviewer.backends.job_cancel import JobCancelledError

        raise JobCancelledError("Job stopped by user")


def _write_fasta(path: Path, items: Iterable[tuple[str, str]], width: int = 80) -> None:
    with path.open("w", encoding="utf-8", newline="\n") as out:
        for name, seq in items:
            out.write(f">{name}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")
            if not seq:
                out.write("\n")


def _parse_fasta_text(text: str) -> list[tuple[str, str]]:
    pairs: list[tuple[str, str]] = []
    name, chunks = None, []
    for line in text.splitlines():
        if line.startswith(">"):
            if name is not None:
                pairs.append((name, "".join(chunks)))
            name, chunks = line[1:].strip(), []
        elif name is not None:
            chunks.append(line.strip())
    if name is not None:
        pairs.append((name, "".join(chunks)))
    return pairs


def _single_query_file(queries: list[Path], work: Path, *, fasta_only: bool = False
                       ) -> Path:
    """The query file for a read mapper: the input itself when there is one file
    (FASTA / FASTQ, gzipped or not), otherwise all inputs combined as FASTA."""
    if len(queries) == 1 and not (fasta_only and is_fastq(queries[0])):
        return queries[0]
    combined = work / "queries.fa"
    _write_fasta(combined, (item for q in queries for item in iter_sequences(q)))
    return combined


# Illumina-style mate tags in the filename (before .fastq / .fq / .gz).
_HISAT2_MATE_RE = re.compile(
    r"(?P<pre>.*?)(?P<tag>"
    r"_R1|_R2"           # sample_R1.fastq.gz / sample_R1_001.fastq
    r"|_1|_2"            # sample_1.fastq / sample_2.fastq
    r"|\.1|\.2"          # sample.1.fq / sample.2.fq (older)
    r")(?P<suf>(?:_\d+)?(?:\.(?:fastq|fq|fasta|fa|fna))+)$",
    re.IGNORECASE,
)


def _hisat2_mate_parts(path: Path) -> tuple[str, int, str] | None:
    """``(pair_key, mate 1|2, suffix)`` from an Illumina-style filename, else None."""
    name = path.name
    lower = name.lower()
    if lower.endswith(".gz"):
        name = name[:-3]
    m = _HISAT2_MATE_RE.match(name)
    if not m:
        return None
    tag = m.group("tag").upper().replace(".", "_")
    if tag in ("_R1", "_1"):
        mate = 1
    elif tag in ("_R2", "_2"):
        mate = 2
    else:
        return None
    key = f"{m.group('pre').lower()}|{m.group('suf').lower()}"
    return key, mate, m.group("suf")


def pair_hisat2_queries(queries: list[Path]) -> tuple[list[Path], list[Path]]:
    """Split Illumina R1/R2 FASTQ paths into mate-1 and mate-2 lists for ``-1``/``-2``.

    Pairs by filename: ``*_R1*``/``*_R2*``, ``*_1.*``/``*_2.*``, or ``*.1.*``/``*.2.*``.
    Raises ``MapReferenceError`` when files cannot be paired.
    """
    if len(queries) < 2:
        raise MapReferenceError(
            "Paired-end HISAT2 needs both mates (R1 and R2).\n\n"
            "Select two FASTQ files whose names differ by R1/R2 (or _1/_2), "
            "for example sample_R1.fastq.gz and sample_R2.fastq.gz."
        )
    if len(queries) % 2 != 0:
        raise MapReferenceError(
            "Paired-end HISAT2 needs an even number of FASTQ files "
            f"(got {len(queries)}).\n\n"
            "Each R1 file must have a matching R2 file."
        )
    buckets: dict[str, dict[int, Path]] = {}
    unpaired: list[Path] = []
    for path in queries:
        parts = _hisat2_mate_parts(path)
        if parts is None:
            unpaired.append(path)
            continue
        key, mate, _suf = parts
        slot = buckets.setdefault(key, {})
        if mate in slot:
            raise MapReferenceError(
                f"Two files look like mate {mate} for the same sample:\n"
                f"  {slot[mate].name}\n  {path.name}"
            )
        slot[mate] = path
    if unpaired:
        names = "\n".join(f"  {p.name}" for p in unpaired)
        raise MapReferenceError(
            "These files do not look like Illumina R1/R2 pairs:\n"
            f"{names}\n\n"
            "Rename them so mates share a name and differ by R1/R2 "
            "(e.g. sample_R1.fastq.gz / sample_R2.fastq.gz), "
            "or switch Reads to Single-end."
        )
    mates1: list[Path] = []
    mates2: list[Path] = []
    for key, slot in sorted(buckets.items()):
        if set(slot) != {1, 2}:
            have = ", ".join(p.name for p in slot.values())
            missing = "R2" if 1 in slot else "R1"
            raise MapReferenceError(
                f"Incomplete pair (missing {missing}):\n  {have}"
            )
        mates1.append(slot[1])
        mates2.append(slot[2])
    return mates1, mates2


def _mafft_version(command: list[str], env: dict) -> str:
    try:
        from repeatviewer.conda_tools import prepare_command

        r = subprocess.run(prepare_command([*command, "--version"]), capture_output=True,
                           text=True, timeout=30, env=env)
    except (OSError, subprocess.SubprocessError):
        return ""
    m = re.search(r"v?(\d+\.\d+)", (r.stdout or "") + (r.stderr or ""))
    return m.group(1) if m else ""


def _summarize_error(tool: str, returncode: int, output: str) -> str:
    lines = [ln for ln in (output or "").splitlines() if ln.strip()]
    tail = "\n".join(lines[-12:])
    return f"{tool} failed (exit {returncode})." + (f"\n\n{tail}" if tail else "")


def hisat2_index_dir(reference: Path) -> Path:
    ref = Path(reference).expanduser().resolve()
    st = ref.stat()
    sig = hashlib.md5(f"{ref.name}|{st.st_size}|{st.st_mtime_ns}".encode()).hexdigest()[:10]
    return ref.parent / ".repeatviewer_hisat2" / f"{ref.name}_{sig}"


def ensure_hisat2_index(reference: Path, hisat2: str, threads: int, commands: list[str],
                        say, cancel) -> Path:
    """Build (once) and return the index prefix; kept next to the reference."""
    from repeatviewer.external_process import command_text, run_logged_process

    folder = hisat2_index_dir(reference)
    try:
        folder.mkdir(parents=True, exist_ok=True)
    except OSError:
        # reference folder is read-only: keep the index in the user's home
        folder = Path.home() / ".repeatviewer" / "hisat2" / folder.name
    prefix = folder / "index"
    done = folder / "BUILT"
    # Small indexes use .1.ht2; large genomes use .1.ht2l.
    if done.is_file() and (
        Path(f"{prefix}.1.ht2").is_file() or Path(f"{prefix}.1.ht2l").is_file()
    ):
        return prefix
    build = hisat2_build_executable(hisat2)
    if not build:
        raise MapReferenceError("hisat2-build was not found next to hisat2.\n\n"
                                + engine_install_hint(MAP_ENGINE_HISAT2))
    folder.mkdir(parents=True, exist_ok=True)
    say("Indexing the reference for HISAT2 (once)…")
    cmd = [build, "-p", str(max(1, threads)), str(reference), str(prefix)]
    commands.append(command_text(cmd))
    result = run_logged_process(cmd, tool="hisat2-build", cwd=folder, cancel=cancel)
    _check_cancel(cancel)
    if result.returncode != 0:
        raise MapReferenceError(_summarize_error("hisat2-build", result.returncode,
                                                 result.stderr or result.stdout))
    done.write_text(f"{reference}\n", encoding="utf-8")
    return prefix


@dataclass
class EngineOutput:
    records: Iterator[tuple[str, int, MappedRead | None]]
    raw_path: Path | None
    raw_name: str
    version: str


def run_mafft_engine(reference: Path, ref_id: str, ref_seq: str, queries: list[Path],
                     opts: MapOptions, work: Path, commands: list[str], cancel
                     ) -> tuple[list[MappedRead], list[str], list[tuple[str, str]], str]:
    """MAFFT --addfragments: (placements, unmapped names, alignment rows with
    original names, MAFFT version)."""
    from repeatviewer.backends.sequence_aligner import mafft_command
    from repeatviewer.external_process import command_text, run_logged_process

    cmd = mafft_command(use_system_mafft=opts.use_system_path)
    if cmd is None:
        raise MapReferenceError("MAFFT was not found (bundled or on PATH).")
    command, env = cmd
    ref_file = work / "reference.fa"
    _write_fasta(ref_file, [("REF", ref_seq.replace("-", ""))])
    safe_to_name: dict[str, str] = {}
    frags: list[tuple[str, str]] = []
    empty: list[str] = []
    for q in queries:
        for name, seq in iter_sequences(q):
            seq = seq.replace("-", "").replace(".", "")
            if not seq:
                empty.append(name)
                continue
            safe = f"Q{len(frags) + 1:06d}"
            safe_to_name[safe] = name
            frags.append((safe, seq))
    if not frags:
        raise MapReferenceError("The selected files contain no sequences.")
    frag_file = work / "fragments.fa"
    _write_fasta(frag_file, frags)
    full = [*command, *mafft_arguments(opts), "--addfragments", str(frag_file), str(ref_file)]
    commands.append(command_text(full))
    result = run_logged_process(full, tool="MAFFT", cwd=work, env=env, stdout_is_result=True,
                                cancel=cancel)
    _check_cancel(cancel)
    if result.returncode != 0 or not (result.stdout or "").strip():
        raise MapReferenceError(_summarize_error("MAFFT", result.returncode, result.stderr))
    pairs = _parse_fasta_text(result.stdout)
    reads, unmapped = parse_mafft_alignment(pairs, "REF", ref_id, safe_to_name)
    named = []
    for name, row in pairs:
        if name == "REF":
            named.append((ref_id, row.upper()))
        else:
            rev = name.startswith("_R_")
            base = safe_to_name.get(name[3:] if rev else name, name)
            named.append((f"{base} (reversed)" if rev else base, row.upper()))
    return reads, unmapped + empty, named, _mafft_version(command, env)


def run_minimap2_engine(reference: Path, queries: list[Path], opts: MapOptions,
                        info: InputInfo, work: Path, commands: list[str], cancel
                        ) -> tuple[Path, str, str, str]:
    """Run minimap2; returns (alignment path, preset, version, format ``paf``|``sam``)."""
    from repeatviewer.external_process import command_text, run_logged_process

    exe = engine_executable(MAP_ENGINE_MINIMAP2, system_path_only=opts.use_system_path)
    if not exe:
        raise MapReferenceError("minimap2 was not found.\n\n"
                                + engine_install_hint(MAP_ENGINE_MINIMAP2))
    preset = (auto_minimap2_preset(info) if opts.minimap2_preset == MAP_MINIMAP2_PRESET_AUTO
              else opts.minimap2_preset)
    out_fmt = getattr(opts, "minimap2_output", MAP_MINIMAP2_OUT_PAF) or MAP_MINIMAP2_OUT_PAF
    # PAF without cs cannot build the contig view.
    if out_fmt == MAP_MINIMAP2_OUT_PAF and not getattr(opts, "minimap2_include_cs", True):
        raise MapReferenceError(
            "PAF output needs Include CS tag (--cs=long) for the contig view.\n"
            "Tick it under Minimap2 → Output, or choose SAM."
        )
    query = _single_query_file(queries, work)
    if out_fmt == MAP_MINIMAP2_OUT_SAM:
        out_path = work / "minimap2.sam"
    else:
        out_fmt = MAP_MINIMAP2_OUT_PAF
        out_path = work / "minimap2.paf"
    cmd = [exe, *minimap2_arguments(opts, preset), "-o", str(out_path),
           str(reference), str(query)]
    commands.append(command_text(cmd))
    result = run_logged_process(cmd, tool="minimap2", cwd=work, cancel=cancel)
    _check_cancel(cancel)
    if result.returncode != 0 or not out_path.is_file():
        raise MapReferenceError(_summarize_error("minimap2", result.returncode,
                                                 result.stderr or result.stdout))
    return out_path, preset, tool_version(exe), out_fmt


def run_hisat2_engine(reference: Path, queries: list[Path], opts: MapOptions, work: Path,
                      commands: list[str], say, cancel) -> tuple[Path, str]:
    """Run HISAT2 (building the index once); returns (SAM path, version)."""
    from repeatviewer.external_process import command_text, run_logged_process

    exe = engine_executable(MAP_ENGINE_HISAT2, system_path_only=opts.use_system_path)
    if not exe:
        raise MapReferenceError(
            "HISAT2 is not configured.\n\n"
            + engine_install_hint(MAP_ENGINE_HISAT2)
        )
    prefix = ensure_hisat2_index(reference, exe, opts.threads, commands, say, cancel)
    paired = getattr(opts, "hisat2_reads", MAP_HISAT2_READS_SINGLE) == MAP_HISAT2_READS_PAIRED
    sam = work / "hisat2.sam"
    cmd = [exe, *hisat2_arguments(opts), "-x", str(prefix)]
    if paired:
        if not all(is_fastq(q) for q in queries):
            raise MapReferenceError(
                "Paired-end HISAT2 expects FASTQ mates (R1/R2).\n\n"
                "FASTA inputs are treated as single-end — switch Reads to Single-end."
            )
        mates1, mates2 = pair_hisat2_queries(queries)
        cmd += [
            "-1", ",".join(str(p) for p in mates1),
            "-2", ",".join(str(p) for p in mates2),
        ]
        say(f"Mapping {len(mates1)} Illumina pair(s) with HISAT2…")
    else:
        fastq = all(is_fastq(q) for q in queries)
        query = queries[0] if len(queries) == 1 else _single_query_file(
            queries, work, fasta_only=True
        )
        if len(queries) > 1:
            fastq = False
        cmd += ["-U", str(query), *([] if fastq else ["-f"])]
        say("Mapping with HISAT2…")
    cmd += ["-S", str(sam)]
    commands.append(command_text(cmd))
    result = run_logged_process(cmd, tool="HISAT2", cwd=work, cancel=cancel)
    _check_cancel(cancel)
    if result.returncode != 0 or not sam.is_file():
        raise MapReferenceError(_summarize_error("HISAT2", result.returncode,
                                                 result.stderr or result.stdout))
    return sam, tool_version(exe)


def iter_paf(path: Path) -> Iterator[tuple[str, int, MappedRead | None]]:
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.strip():
                yield parse_paf_line(line)


def iter_sam(path: Path) -> Iterator[tuple[str, int, MappedRead | None]]:
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("@") or not line.strip():
                continue
            yield parse_sam_line(line)


# --------------------------------------------------------------------------
# The whole run
# --------------------------------------------------------------------------


@dataclass
class ContigResult:
    ref_id: str
    length: int
    detailed: bool
    n_reads: int
    reference: str = ""
    consensus: str = ""
    ref_to_cons: np.ndarray | None = None
    coverage: np.ndarray | None = None
    bin_size: int = 1
    variants: list[Variant] = field(default_factory=list)
    features: list = field(default_factory=list)
    feature_summary: list[dict] = field(default_factory=list)
    display: list[tuple[MappedRead, list[LiftedFeature], int]] = field(default_factory=list)
    mean_identity: float = 0.0


@dataclass
class MapRun:
    reference: Path
    queries: list[Path]
    engine: str
    engine_version: str
    preset: str
    opts: MapOptions
    n_queries: int
    contigs: list[ContigResult]
    records: list[ReadRecord]
    n_placements: int
    placed_names: set[str]
    reversed_names: set[str]
    unused: dict[str, str]
    commands: list[str]
    warnings: list[str]
    stamp: str


def fasta_id(record: ReadRecord) -> str:
    rid = record.name + ("_reversed" if record.strand == "-" else "")
    return rid if record.placement == 1 else f"{rid}_p{record.placement}"


def feature_summary(features: list, lifted_lists: Iterable[list[LiftedFeature]]) -> list[dict]:
    """Per reference feature: how many placements have it intact / changed / … ."""
    rows = [{"intact": 0, "changed": 0, "partial": 0, "lost": 0, "identity_sum": 0.0,
             "n": 0} for _ in features]
    key = {STATUS_INTACT: "intact", STATUS_CHANGED: "changed", STATUS_PARTIAL: "partial",
           STATUS_LOST: "lost"}
    for lifted in lifted_lists:
        for lf in lifted:
            if 0 <= lf.index < len(rows):
                r = rows[lf.index]
                r[key[lf.status]] += 1
                if lf.status in (STATUS_INTACT, STATUS_CHANGED):
                    r["identity_sum"] += lf.identity
                    r["n"] += 1
    out = []
    for feat, r in zip(features, rows):
        full = r["n"]
        out.append({"name": feat.name, "type": feat.feature_type, "start": feat.start,
                    "end": feat.end, "intact": r["intact"], "changed": r["changed"],
                    "partial": r["partial"], "lost": r["lost"],
                    "mean_identity": round(r["identity_sum"] / full, 4) if full else None})
    return out


def run_mapping(
    queries: list[Path],
    opts: MapOptions,
    *,
    progress: Callable[[str], None] | None = None,
    cancel=None,
    stamp: str | None = None,
) -> tuple[Path, list[Path], Path]:
    """Map ``queries`` to ``opts.reference_fasta``.

    Returns (``.rvmap`` when ``opts.generate_rvmap``, else the alignment path,
    created files, report HTML)."""
    say = progress or (lambda _m: None)
    queries = [Path(q).expanduser() for q in queries]
    ref_text = str(opts.reference_fasta or "").strip()
    if not ref_text:
        raise MapReferenceError("Choose a reference FASTA first (Mapping tab → Reference).")
    reference = Path(ref_text).expanduser()
    if not reference.is_file():
        raise MapReferenceError(f"Reference not found:\n{reference}")
    reference = reference.resolve()
    warnings: list[str] = []
    kept_queries = [q.resolve() for q in queries if q.is_file() and q.resolve() != reference]
    if len(kept_queries) < len([q for q in queries if q.is_file()]):
        warnings.append("The reference was also selected as input; it was not mapped to itself.")
    queries = kept_queries
    if not queries:
        raise MapReferenceError("Select the FASTA / FASTQ files to map in the file panel.")

    say("Reading the input…")
    duplicate_refs: list[str] = []
    ref_lengths = reference_lengths(reference, duplicate_refs)
    if not ref_lengths or not any(ref_lengths.values()):
        raise MapReferenceError(f"{reference.name} contains no sequences.")
    if duplicate_refs:
        # minimap2 / HISAT2 report hits by name only: with two records of the
        # same name the result would be drawn against the wrong sequence.
        shown = ", ".join(duplicate_refs[:5]) + (" …" if len(duplicate_refs) > 5 else "")
        raise MapReferenceError(
            f"{reference.name} has more than one sequence named {shown}.\n\n"
            "Give every sequence in the reference its own name (edit the headers, for "
            "example with Open in Text Editor) and run Mapping again.")
    info = scan_inputs(queries, max_records=50_000)
    if info.n == 0:
        raise MapReferenceError("The selected files contain no sequences.")
    if info.complete and len(set(info.names)) < len(info.names):
        warnings.append("Some sequence names occur more than once; their results are merged. "
                        "Rename them for separate results.")
    engine = opts.engine
    if engine == MAP_ENGINE_MAFFT:
        warnings.append(
            "MAFFT is no longer used for Mapping; Auto chose the mapper instead."
        )
        engine = MAP_ENGINE_AUTO
    if engine == MAP_ENGINE_AUTO:
        engine = choose_engine(ref_lengths, info)
        if not engine_available(engine, system_path_only=opts.use_system_path):
            fallback = next((e for e in (MAP_ENGINE_MINIMAP2, MAP_ENGINE_HISAT2)
                             if e != engine and engine_available(
                                 e, system_path_only=opts.use_system_path)), None)
            if fallback is None:
                raise MapReferenceError(
                    f"Auto chose {ENGINE_NAMES[engine]} for this input, but it was not found.\n\n"
                    + engine_install_hint(engine))
            warnings.append(f"Auto would use {ENGINE_NAMES[engine]} for this input, which was "
                            f"not found; {ENGINE_NAMES[fallback]} was used instead.")
            engine = fallback
    elif not engine_available(engine, system_path_only=opts.use_system_path):
        if engine == MAP_ENGINE_HISAT2:
            raise MapReferenceError(
                "HISAT2 is not configured.\n\n" + engine_install_hint(engine)
            )
        raise MapReferenceError(
            f"{ENGINE_NAMES.get(engine, engine)} was not found.\n\n"
            + engine_install_hint(engine)
        )

    features_by_id = (load_reference_features(reference, ref_lengths)
                      if opts.lift_annotations else {})
    stamp = stamp or run_stamp()
    layout = plan_outputs(queries, engine=engine)
    created: list[Path] = []
    from repeatviewer.backends.job_cancel import JobCancelledError

    def _cleanup() -> None:
        layout.rvmap.unlink(missing_ok=True)
        for ext in (".paf", ".sam"):
            layout.alignment(ext).unlink(missing_ok=True)
            run_json_path(layout.alignment(ext)).unlink(missing_ok=True)
        shutil.rmtree(layout.details, ignore_errors=True)

    try:
        run = _map(layout, reference, ref_lengths, features_by_id, queries, info, engine,
                   opts, warnings, stamp, say, cancel)
        say("Writing results…")
        return _write_outputs(run, layout, info, created, say, cancel)
    except JobCancelledError:
        _cleanup()
        raise
    except Exception as exc:  # noqa: BLE001 - reported with a saved log
        message = str(exc) if isinstance(exc, MapReferenceError) else f"{type(exc).__name__}: {exc}"
        log = _write_failure_log(queries, message)
        _cleanup()
        raise MapReferenceError(message + (f"\n\nFull log: {log}" if log else "")) from exc


def _write_failure_log(queries: list[Path], error: str) -> Path | None:
    import traceback

    path = Path(queries[0]).expanduser().resolve().parent / \
        f"{output_base_name(queries)}_map_failed.log"
    try:
        path.write_text(f"Mapping failed {datetime.now():%Y-%m-%d %H:%M}\n\n{error}\n\n"
                        + traceback.format_exc(), encoding="utf-8")
        return path
    except OSError:
        return None


def _map(layout: OutputLayout, reference: Path, ref_lengths: dict[str, int],
         features_by_id: dict[str, list], queries: list[Path], info: InputInfo, engine: str,
         opts: MapOptions, warnings: list[str], stamp: str, say, cancel) -> MapRun:
    work = layout.details / "work"
    work.mkdir(parents=True, exist_ok=True)
    commands: list[str] = []
    preset = ""
    ref_seqs: dict[str, str] = {}
    unmapped_extra: list[str] = []
    if engine == MAP_ENGINE_MINIMAP2:
        say("Mapping with minimap2…")
        aln, preset, version, out_fmt = run_minimap2_engine(
            reference, queries, opts, info, work, commands, cancel
        )
        if out_fmt == MAP_MINIMAP2_OUT_SAM:
            raw = layout.alignment(".sam")
            shutil.move(str(aln), str(raw))
            records = iter_sam(raw)
        else:
            raw = layout.alignment(".paf")
            shutil.move(str(aln), str(raw))
            records = iter_paf(raw)
    elif engine == MAP_ENGINE_HISAT2:
        sam, version = run_hisat2_engine(reference, queries, opts, work, commands, say, cancel)
        raw = layout.alignment(".sam")
        shutil.move(str(sam), str(raw))
        records = iter_sam(raw)
        library = getattr(opts, "hisat2_library", MAP_HISAT2_RNA) or MAP_HISAT2_RNA
        preset = "dna" if library == MAP_HISAT2_DNA else "rna"
    else:
        raise MapReferenceError(
            f"Unknown Mapping engine: {engine!r}. Use minimap2 or HISAT2."
        )

    say("Collecting the placements…")
    detail_budget = 50_000_000
    accs: dict[str, ContigAccumulator] = {}
    indexes: dict[str, FeatureIndex] = {}
    placements: Counter = Counter()
    unused: dict[str, str] = {n: "no alignment found" for n in unmapped_extra}
    placed: set[str] = set()
    reversed_names: set[str] = set()
    table: list[ReadRecord] = []
    n_placements = 0
    seen: set[str] = set()
    for count, (qname, _qlen, read) in enumerate(records):
        if count % 2000 == 0:
            _check_cancel(cancel)
        seen.add(qname)
        if read is None:
            unused.setdefault(qname, "no alignment found")
            continue
        reason = placement_verdict(read, opts, engine)
        if reason is not None:
            if reason:
                unused.setdefault(qname, reason)
            continue
        acc = accs.get(read.ref_id)
        if acc is None:
            length = ref_lengths.get(read.ref_id, read.end)
            detailed = length <= DETAIL_LIMIT and detail_budget - length >= 0
            if detailed:
                detail_budget -= length
            feats = features_by_id.get(read.ref_id, [])
            acc = ContigAccumulator(read.ref_id, length, detailed, feats,
                                    max(1, opts.max_display_reads))
            accs[read.ref_id] = acc
            indexes[read.ref_id] = FeatureIndex(feats)
        lifted = (lift_features(read, indexes[read.ref_id])
                  if acc.features and n_placements < MAX_READ_TABLE else [])
        placements[qname] += 1
        acc.add(read, lifted, placements[qname])
        n_placements += 1
        placed.add(qname)
        if read.strand == "-":
            reversed_names.add(qname)
        if len(table) < MAX_READ_TABLE:
            table.append(ReadRecord(
                read.name, read.ref_id, read.start, read.end, read.strand,
                round(read.identity, 4), read.mismatches, read.deleted, read.inserted,
                read.aligned_bp, read.mapq, read.kind, read.query_length, read.clipped,
                read.clip_left, lifted, placements[qname]))
    # minimap2 (--paf-no-hit) and HISAT2 list every read in their output.
    for name in seen:
        if name not in placed:
            unused.setdefault(name, "no alignment found")
    for name in placed:
        unused.pop(name, None)
    n_queries = len(seen | set(unmapped_extra))
    if n_placements > MAX_READ_TABLE:
        warnings.append(f"{n_placements:,} placements: the read and feature tables list the "
                        f"first {MAX_READ_TABLE:,}; coverage, consensus and variants use all.")
    if any(not a.detailed for a in accs.values()):
        warnings.append("Some reference sequences are too long for base-level detail "
                        f"(over {DETAIL_LIMIT:,} bp, or {50:,} Mb in total): they show binned "
                        "coverage and the mapped reads, without consensus or variants.")

    say("Consensus and variants…")
    wanted = {rid for rid, a in accs.items() if a.detailed and rid not in ref_seqs}
    if wanted:
        ref_seqs.update(load_reference_records(reference, wanted=wanted))
    contigs: list[ContigResult] = []
    for rid, acc in sorted(accs.items(), key=lambda kv: -kv[1].n):
        _check_cancel(cancel)
        c = ContigResult(rid, acc.length, acc.detailed, acc.n, features=acc.features,
                         display=acc.display, bin_size=acc.bin_size,
                         mean_identity=float(np.mean(acc.identities)) if acc.identities else 0.0)
        c.coverage = acc.coverage()
        if acc.detailed:
            ref = ref_seqs.get(rid, "")
            c.reference = ref
            for read, _lifted, _pl in acc.display:
                if "?" in read.row:
                    fill_unknown_bases(read, ref)
            c.consensus, c.ref_to_cons = build_consensus(
                ref, acc.counts, acc.insertions, threshold=opts.consensus_threshold,
                no_coverage=opts.no_coverage)
            c.variants = call_variants(ref, acc.counts, acc.insertions, acc.features,
                                       min_freq=opts.min_variant_freq)
        c.feature_summary = feature_summary(
            acc.features, (r.lifted for r in table if r.ref_id == rid))
        contigs.append(c)
    shutil.rmtree(work, ignore_errors=True)
    return MapRun(reference, queries, engine, version, preset, opts, n_queries, contigs, table,
                  n_placements, placed, reversed_names, unused, commands, warnings, stamp)


# --------------------------------------------------------------------------
# Outputs
# --------------------------------------------------------------------------

CITATIONS = {
    MAP_ENGINE_MINIMAP2: [
        "Li H (2018) Minimap2: pairwise alignment for nucleotide sequences. "
        "Bioinformatics 34:3094–3100.",
    ],
    MAP_ENGINE_HISAT2: [
        "Kim D, Paggi JM, Park C, Bennett C, Salzberg SL (2019) Graph-based genome "
        "alignment and genotyping with HISAT2 and HISAT-genotype. Nat Biotechnol "
        "37:907–915.",
    ],
}


def threshold_label(threshold: float) -> str:
    return "majority" if threshold <= 0 else f"{threshold:.0%} threshold (IUPAC codes)"


def methods_paragraph(run: MapRun) -> str:
    o = run.opts
    name = ENGINE_NAMES[run.engine] + (f" v{run.engine_version}" if run.engine_version else "")
    if run.engine == MAP_ENGINE_MINIMAP2:
        out = getattr(o, "minimap2_output", MAP_MINIMAP2_OUT_PAF) or MAP_MINIMAP2_OUT_PAF
        fmt = "SAM (-a --MD)" if out == MAP_MINIMAP2_OUT_SAM else "PAF (-c --cs=long)"
        seeds = getattr(o, "minimap2_repeat_seeds", MAP_REPEAT_SEEDS_DEFAULT)
        seed_note = {"more": ", keeping more repetitive seeds with -f 0.00002",
                     "all": ", keeping all repetitive seeds with -f 1000000"}.get(seeds, "")
        how = f"with {name} (preset {run.preset}, {fmt}{seed_note})"
    elif run.engine == MAP_ENGINE_HISAT2:
        how = f"with {name} ({run.preset or 'rna'})"
    else:
        how = f"with {name}"
    multi = {"best": "Sequences fitting equally well in several places were placed once, "
                     "in the place chosen by the mapper.",
             "all": "Sequences fitting equally well in several places were placed in every "
                    "such place.",
             "none": "Sequences fitting equally well in several places were not placed."}
    text = (f"{run.n_queries} sequence(s) were mapped to {run.reference.name} {how}. "
            f"Placements with less than {o.min_identity:.0%} identity or fewer than "
            f"{o.min_aligned_bp} aligned bases were not used"
            + (f", nor those with mapping quality below {o.min_mapq}"
               if o.min_mapq else "") + ". "
            + multi[o.multiple_hits] + " "
            + f"{len(run.placed_names)} sequence(s) were placed ({run.n_placements} placements"
            f"{', ' + str(len(run.reversed_names)) + ' reverse-complemented' if run.reversed_names else ''}). "
            f"The consensus was called by {threshold_label(o.consensus_threshold)}. ")
    if any(c.features for c in run.contigs):
        text += ("The reference annotations were transferred to every placed sequence and "
                 "classified as intact, changed, partial or lost. ")
    text += "Mapping was run and summarised in RepeatViewer."
    return text


REPORT_STATUS_COLORS = {STATUS_INTACT: "#5fae66", STATUS_CHANGED: "#e0b43a",
                        STATUS_PARTIAL: "#b9b9b5", STATUS_LOST: "#d0574a"}
REPORT_VARIABLE_SHARE = 0.25     # same threshold as the contig view's identity bar
REPORT_MAX_DIFFERENCES = 20


def _svg_coverage(c: ContigResult, width: int = 940, height: int = 84) -> str:
    """Coverage profile with the reference features as a band below; scales to the
    page width so it lines up with the tables."""
    cov = np.asarray(c.coverage if c.coverage is not None else [], dtype=float)
    if cov.size == 0:
        return ""
    n = min(width, cov.size)
    edges = np.linspace(0, cov.size, n + 1).astype(int)
    vals = np.array([cov[a:b].mean() if b > a else 0 for a, b in zip(edges[:-1], edges[1:])])
    top = max(1.0, float(vals.max()))
    base = height - 22
    pts = " ".join(f"{i * width / n:.1f},{base - v / top * (base - 6):.1f}"
                   for i, v in enumerate(vals))
    feats = ""
    L = max(1, c.length)
    for f in c.features[:400]:
        x = (f.start - 1) / L * width
        w = max(1.0, (f.end - f.start + 1) / L * width)
        feats += (f'<rect x="{x:.1f}" y="{base + 4}" width="{w:.1f}" height="6" rx="1" '
                  f'fill="{html.escape(f.color or "#888")}"><title>{html.escape(f.name)} '
                  f'{f.start:,}–{f.end:,}</title></rect>')
    return (f'<svg class="cov" viewBox="0 0 {width} {height}" preserveAspectRatio="none">'
            f'<polyline points="0,{base} {pts} {width},{base}" fill="#dbe7f6" '
            f'stroke="#3b6fb6" stroke-width="1" vector-effect="non-scaling-stroke"/>{feats}'
            f'</svg><div class="axis"><span>1</span><span>max {top:.0f}×</span>'
            f'<span>{c.length:,} bp</span></div>')


def _stamp_text(stamp: str) -> str:
    """'2026-10-27_0812' → 'Run on 2026-10-27 at 08:12'."""
    try:
        day, hm = stamp.split("_", 1)
        return f"Run on {day} at {hm[:2]}:{hm[2:4]}"
    except (ValueError, AttributeError):
        return ""


def _status_bar(fs: dict) -> str:
    total = sum(int(fs.get(k, 0) or 0) for k in REPORT_STATUS_COLORS)
    if not total:
        return ""
    parts = "".join(
        f'<i style="width:{100 * int(fs.get(k, 0) or 0) / total:.2f}%;background:{col}" '
        f'title="{int(fs.get(k, 0) or 0)} {k}"></i>'
        for k, col in REPORT_STATUS_COLORS.items() if int(fs.get(k, 0) or 0))
    return f'<span class="bar">{parts}</span>'


def _type_legend(c: ContigResult) -> str:
    seen: dict[str, str] = {}
    for f in c.features:
        seen.setdefault(str(f.feature_type or "other"), f.color or "#888")
    return "".join(f'<span><b style="background:{html.escape(col)}"></b>{html.escape(t)}</span>'
                   for t, col in seen.items())


def report_html(run: MapRun, files: list[Path], title: str) -> str:
    e = html.escape
    o = run.opts
    ident = [r.identity for r in run.records]
    n_var_pos = sum(len({v.pos for v in c.variants}) for c in run.contigs)
    tiles = [
        ("Placed", f"{len(run.placed_names):,} / {run.n_queries:,}",
         f"{run.n_placements:,} placements"),
        ("Reversed to fit", f"{len(run.reversed_names):,}",
         f"{len(run.unused):,} not used"),
        ("Mean identity", f"{statistics.mean(ident):.1%}" if ident else "–",
         f"min {min(ident):.1%}" if ident else ""),
        ("Variant positions", f"{n_var_pos:,}",
         f"{len(run.contigs):,} reference sequence(s)"),
    ]
    tile_html = "".join(f'<div class="tile"><div class="k">{e(k)}</div><div class="v">{e(v)}'
                        f'</div><div class="s">{e(sub)}</div></div>' for k, v, sub in tiles)
    seeds = getattr(o, "minimap2_repeat_seeds", MAP_REPEAT_SEEDS_DEFAULT)
    rows = [
        ("Reference", run.reference.name),
        ("Input", ", ".join(q.name for q in run.queries)),
        ("Mapper", ENGINE_NAMES[run.engine]
         + (f" {run.engine_version}" if run.engine_version else "")
         + (f" · {run.preset}" if run.preset else "")
         + ({"more": " · more repetitive seeds", "all": " · all repetitive seeds"}
            .get(seeds, "") if run.engine == MAP_ENGINE_MINIMAP2 else "")),
        ("Limits", f"identity ≥ {o.min_identity:.0%} · ≥ {o.min_aligned_bp} aligned bp"
         + (f" · mapping quality ≥ {o.min_mapq}" if o.min_mapq else "")),
        ("Consensus", threshold_label(o.consensus_threshold)),
    ]
    details = "".join(f"<tr><th>{e(k)}</th><td>{e(v)}</td></tr>" for k, v in rows)
    warn = "".join(f"<li>{e(w)}</li>" for w in run.warnings)
    one = len(run.contigs) == 1
    contig_rows = ""
    for c in run.contigs[:50]:
        mean_cov = float(np.mean(c.coverage)) if c.coverage is not None and len(c.coverage) else 0
        contig_rows += (f"<tr><td>{e(c.ref_id)}</td><td class=n>{c.length:,}</td>"
                        f"<td class=n>{c.n_reads:,}</td><td class=n>{mean_cov:.1f}×</td>"
                        f"<td class=n>{c.mean_identity:.1%}</td>"
                        f"<td class=n>{len(c.variants):,}</td></tr>")
    plots = "".join(
        (f"<h3>{e(c.ref_id)}</h3>" if not one else "") + _svg_coverage(c)
        + (f'<div class="legend">{_type_legend(c)}</div>' if c.features else "")
        for c in run.contigs[:5])
    feat_rows = ""
    for c in run.contigs[:5]:
        for fs in c.feature_summary[:MAX_FEATURE_COLUMNS]:
            ident_txt = f"{fs['mean_identity']:.1%}" if fs["mean_identity"] is not None else "–"
            feat_rows += "<tr>" + ("" if one else f"<td>{e(c.ref_id)}</td>")
            feat_rows += (f"<td>{e(fs['name'])}</td><td class=m>{e(fs['type'])}</td>"
                          f"<td class=n>{fs['start']:,}–{fs['end']:,}</td>"
                          f"<td class=n>{fs['intact']}</td><td class=n>{fs['changed']}</td>"
                          f"<td class=n>{fs['partial']}</td><td class=n>{fs['lost']}</td>"
                          f"<td class=n>{ident_txt}</td><td>{_status_bar(fs)}</td></tr>")
    all_vars = [(c, v) for c in run.contigs for v in c.variants]
    common = sorted((cv for cv in all_vars if cv[1].freq >= REPORT_VARIABLE_SHARE),
                    key=lambda cv: (-cv[1].freq, cv[1].pos))
    var_note = ""
    if not common and all_vars:
        common = sorted(all_vars, key=lambda cv: (-cv[1].freq, cv[1].pos))[:10]
        var_note = ("No position differs in a quarter or more of the sequences; "
                    "the most frequent differences are listed.")
    shown_vars = common[:REPORT_MAX_DIFFERENCES]
    var_rows = "".join(
        "<tr>" + (f"<td>{e(c.ref_id)}</td>" if not one else "")
        + f"<td class=n>{v.pos:,}</td><td class=c>{e(v.ref)} → {e(v.alt)}</td>"
          f"<td class=n>{v.count}/{v.depth}</td><td class=n>{v.freq:.0%}</td>"
          f"<td class=m>{e(v.context) or '–'}</td><td class=m>{e(v.features) or '–'}</td></tr>"
        for c, v in shown_vars)
    var_tsv = next((p.name for p in files if p.name.endswith("_variants.tsv")), "")
    rest = len(all_vars) - len(shown_vars)
    var_foot = (f'<p class="note">{e(var_note)} {rest:,} more in {e(var_tsv)}.</p>'
                if rest > 0 or var_note else "")
    reasons = Counter(run.unused.values())
    unused_rows = "".join(f"<li>{n:,} × {e(r)}</li>" for r, n in reasons.most_common())
    main = [p for p in files if not p.parent.name.endswith("_details")]
    extra = [p for p in files if p not in main]
    flist = "".join(f"<li>{e(p.name)}</li>" for p in main)
    if extra:
        flist += (f"<li>{e(extra[0].parent.name)}/<ul>"
                  + "".join(f"<li>{e(p.name)}</li>" for p in extra) + "</ul></li>")
    cites = "".join(f"<li>{e(c)}</li>" for c in CITATIONS[run.engine])
    cmds = e("\n".join(run.commands))
    ref_col = "" if one else "<th>Sequence</th>"
    legend_status = "".join(f'<span><b style="background:{col}"></b>{k}</span>'
                            for k, col in REPORT_STATUS_COLORS.items())
    return f"""<!doctype html><html><head><meta charset="utf-8"><title>{e(title)}</title>
<style>
:root{{--fg:#1d1d1b;--muted:#6b6b67;--line:#e4e4e0;--soft:#f6f6f4}}
body{{font-family:-apple-system,"Segoe UI",Helvetica,Arial,sans-serif;margin:0;color:var(--fg);background:#fff;font-size:13px}}
main{{max-width:940px;margin:0 auto;padding:28px 24px 40px}}
h1{{font-size:19px;font-weight:600;margin:0 0 4px}} .sub{{color:var(--muted);margin:0 0 20px}}
h2{{font-size:14px;font-weight:600;margin:28px 0 8px;padding-bottom:6px;border-bottom:1px solid var(--line)}}
h3{{font-size:12px;font-weight:600;color:var(--muted);margin:14px 0 4px}}
.tiles{{display:grid;grid-template-columns:repeat(4,minmax(0,1fr));gap:12px}}
.tile{{background:var(--soft);border-radius:8px;padding:12px 14px}}
.tile .k{{color:var(--muted);font-size:12px}} .tile .v{{font-size:21px;font-weight:600;margin:2px 0}}
.tile .s{{color:var(--muted);font-size:12px}}
table{{border-collapse:collapse;width:100%}} th,td{{padding:5px 8px;text-align:left;vertical-align:top}}
.grid th{{color:var(--muted);font-weight:500;border-bottom:1px solid var(--line)}}
.grid td{{border-bottom:1px solid #f0f0ec}} td.n,th.n{{text-align:right;font-variant-numeric:tabular-nums}}
td.c{{font-family:Menlo,Consolas,monospace;font-size:12px}} td.m{{color:var(--muted)}}
.details th{{width:120px;color:var(--muted);font-weight:500;padding-left:0}} .details td{{padding-left:0}}
svg.cov{{display:block;width:100%;height:84px}}
.axis{{display:flex;justify-content:space-between;color:var(--muted);font-size:11px;margin-top:2px}}
.legend{{display:flex;flex-wrap:wrap;gap:14px;color:var(--muted);font-size:12px;margin-top:8px}}
.legend b{{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:5px;vertical-align:-1px}}
.bar{{display:flex;width:120px;height:8px;border-radius:4px;overflow:hidden;background:var(--soft);margin-top:4px}}
.bar i{{display:block;height:100%}}
.warn{{background:#fff6e5;border:1px solid #f0c36d;border-radius:8px;padding:10px 14px;margin-bottom:16px}}
.methods{{line-height:1.55;background:var(--soft);border-radius:8px;padding:12px 14px}}
.note{{color:var(--muted);font-size:12px;margin:6px 0 0}}
pre{{font-size:11px;background:var(--soft);border-radius:8px;padding:10px 12px;white-space:pre-wrap;word-break:break-all;margin:0}}
.two{{display:grid;grid-template-columns:repeat(2,minmax(0,1fr));gap:24px}}
ul{{margin:4px 0 0 18px;padding:0}} li{{margin:2px 0}}
@media (max-width:700px){{.tiles{{grid-template-columns:repeat(2,minmax(0,1fr))}}.two{{grid-template-columns:1fr}}}}
</style></head><body><main>
<h1>{e(title)}</h1>
<p class="sub">{e(_stamp_text(getattr(run, 'stamp', '')))}</p>
{f'<div class="warn"><b>Check before use</b><ul>{warn}</ul></div>' if warn else ''}
<div class="tiles">{tile_html}</div>
<h2>Run</h2><table class="details">{details}</table>
<h2>Coverage</h2>{plots}
{'' if one else f'<h2>Reference sequences</h2><table class="grid"><tr><th>Sequence</th><th class=n>Length</th><th class=n>Placements</th><th class=n>Mean coverage</th><th class=n>Mean identity</th><th class=n>Variants</th></tr>{contig_rows}</table>'}
{f'<h2>Reference features in the placed sequences</h2><table class="grid"><tr>{ref_col}<th>Feature</th><th>Type</th><th class=n>Reference</th><th class=n>Intact</th><th class=n>Changed</th><th class=n>Partial</th><th class=n>Lost</th><th class=n>Identity</th><th>Status</th></tr>{feat_rows}</table><div class="legend">{legend_status}</div>' if feat_rows else ''}
{f'<h2>Most frequent differences</h2><table class="grid"><tr>{ref_col}<th class=n>Position</th><th>Change</th><th class=n>Sequences</th><th class=n>Share</th><th>Context</th><th>In feature</th></tr>{var_rows}</table>{var_foot}' if var_rows else ''}
{f'<h2>Not used</h2><ul>{unused_rows}</ul>' if unused_rows else ''}
<h2>Methods text</h2><div class="methods">{e(methods_paragraph(run))}</div>
<h2>Commands</h2><pre>{cmds}</pre>
<div class="two"><div><h2>Cite</h2><ul>{cites}</ul></div><div><h2>Files written</h2><ul>{flist}</ul></div></div>
</main></body></html>"""


def _lifted_feature_objects(record: ReadRecord, features: list) -> list:
    from repeatviewer.backends.sequence_features import SequenceFeature

    out = []
    for lf in record.lifted:
        if lf.start <= 0 or not (0 <= lf.index < len(features)):
            continue
        src = features[lf.index]
        out.append(SequenceFeature(start=lf.start, end=lf.end, name=src.name,
                                   feature_type=src.feature_type, color=src.color,
                                   strand=src.strand, track=src.track, parent=src.parent))
    return out


def _consensus_features(c: ContigResult) -> list:
    from repeatviewer.backends.sequence_features import SequenceFeature

    out = []
    if c.ref_to_cons is None:
        return out
    m = c.ref_to_cons
    for f in c.features:
        seg = m[f.start - 1:f.end]
        placed = seg[seg >= 0]
        if placed.size == 0:
            continue
        out.append(SequenceFeature(start=int(placed.min()) + 1, end=int(placed.max()) + 1,
                                   name=f.name, feature_type=f.feature_type, color=f.color,
                                   strand=f.strand, track=f.track, parent=f.parent))
    return out


def _read_json(read: MappedRead, lifted: list[LiftedFeature], placement: int,
               with_marks: bool = False) -> dict:
    """One drawn placement in the .rvmap file. ``with_marks`` adds the mismatch
    positions for reference sequences stored without their bases."""
    extra = {}
    if with_marks and read.ref_row and len(read.ref_row) == len(read.row):
        extra["mism"] = [[read.start + k, c] for k, (c, r) in enumerate(zip(read.row, read.ref_row))
                         if c != "-" and c != r and r != "?"]
    return {**extra,
        "name": read.name, "start": read.start, "end": read.end, "strand": read.strand,
        "row": read.row, "ins": [[p, s] for p, s in read.insertions],
        "identity": round(read.identity, 4), "mismatches": read.mismatches,
        "deleted": read.deleted, "inserted": read.inserted, "mapq": read.mapq,
        "kind": read.kind, "qlen": read.query_length, "clip_left": read.clip_left,
        "placement": placement,
        "lift": [[lf.index, lf.status, lf.differences, lf.identity, lf.start, lf.end]
                 for lf in lifted],
    }


def _write_outputs(run: MapRun, layout: OutputLayout, info: InputInfo, created: list[Path],
                   say, cancel) -> tuple[Path, list[Path], Path]:
    from repeatviewer.backends.sequence_features import save_annotations_multi

    o = run.opts
    # Prefer the selected format next to the .rvmap; fall back to older runs that
    # kept the aligner file only inside ``_details/``.
    raw = next(
        (
            path
            for path in (
                *(layout.alignment(s) for s in (".paf", ".sam")),
                *(layout.path(s) for s in (".aln.fa", ".paf", ".sam")),
            )
            if path.is_file()
        ),
        None,
    )
    details: dict[str, str] = {}

    # Consensus, with the reference annotations moved to consensus coordinates.
    detailed = [c for c in run.contigs if c.detailed and c.consensus]
    if o.save_consensus and detailed:
        path = layout.path("_consensus.fa")
        _write_fasta(path, [(f"{c.ref_id}_consensus", c.consensus) for c in detailed])
        feats = {f"{c.ref_id}_consensus": _consensus_features(c) for c in detailed}
        if any(feats.values()):
            save_annotations_multi(path, features_by_id=feats,
                                   lengths_by_id={f"{c.ref_id}_consensus": len(c.consensus)
                                                  for c in detailed})
        created.append(path)
        details["consensus"] = path.name

    # Placed and unused sequences: one more pass over the input.
    by_name: dict[str, list[ReadRecord]] = {}
    for rec in run.records:
        by_name.setdefault(rec.name, []).append(rec)
    feats_of = {c.ref_id: c.features for c in run.contigs}
    mapped_path = layout.path("_mapped.fa") if o.save_mapped and run.records else None
    unused_path = layout.path("_unused.fa") if o.save_unused and run.unused else None
    if mapped_path or unused_path:
        say("Writing the placed and unused sequences…")
        mapped_out = mapped_path.open("w", encoding="utf-8", newline="\n") if mapped_path else None
        unused_out = unused_path.open("w", encoding="utf-8", newline="\n") if unused_path else None
        lifted_feats: dict[str, list] = {}
        lengths: dict[str, int] = {}
        annotate = len(run.records) <= 5000
        try:
            for q in run.queries:
                for name, seq in iter_sequences(q):
                    if mapped_out is not None and name in by_name:
                        for rec in by_name[name]:
                            oriented = revcomp(seq) if rec.strand == "-" else seq
                            rid = fasta_id(rec)
                            mapped_out.write(f">{rid} mapped to {rec.ref_id}:{rec.start + 1}-"
                                             f"{rec.end} ({rec.strand})\n")
                            for i in range(0, len(oriented), 80):
                                mapped_out.write(oriented[i:i + 80] + "\n")
                            if annotate:
                                lf = _lifted_feature_objects(rec, feats_of.get(rec.ref_id, []))
                                if lf:
                                    lifted_feats[rid] = lf
                                lengths[rid] = len(oriented)
                    if unused_out is not None and name in run.unused:
                        unused_out.write(f">{name} {run.unused[name]}\n")
                        for i in range(0, len(seq), 80):
                            unused_out.write(seq[i:i + 80] + "\n")
                _check_cancel(cancel)
        finally:
            if mapped_out:
                mapped_out.close()
            if unused_out:
                unused_out.close()
        if mapped_path:
            if lifted_feats:
                save_annotations_multi(mapped_path, features_by_id=lifted_feats,
                                       lengths_by_id=lengths)
            created.append(mapped_path)
            details["mapped"] = mapped_path.name
        if unused_path:
            created.append(unused_path)
            details["unused"] = unused_path.name

    if o.save_variants and any(c.variants for c in run.contigs):
        path = layout.path("_variants.tsv")
        with path.open("w", encoding="utf-8", newline="\n") as out:
            out.write("sequence\tposition\tref\talt\tcount\tdepth\tfrequency\tcontext\tfeatures\n")
            for c in run.contigs:
                for v in c.variants:
                    out.write(f"{c.ref_id}\t{v.pos}\t{v.ref}\t{v.alt}\t{v.count}\t{v.depth}\t"
                              f"{v.freq:.4f}\t{v.context}\t{v.features}\n")
        created.append(path)
        details["variants"] = path.name

    if any(r.lifted for r in run.records):
        path = layout.path("_features.tsv")
        with path.open("w", encoding="utf-8", newline="\n") as out:
            out.write("sequence\treference\tfeature\ttype\tref_start\tref_end\tstatus\t"
                      "differences\tidentity\tstart_on_sequence\tend_on_sequence\n")
            for r in run.records:
                feats = feats_of.get(r.ref_id, [])
                for lf in r.lifted:
                    f = feats[lf.index]
                    out.write(f"{fasta_id(r)}\t{r.ref_id}\t{f.name}\t{f.feature_type}\t{f.start}\t"
                              f"{f.end}\t{lf.status}\t{lf.differences}\t{lf.identity:.4f}\t"
                              f"{lf.start or ''}\t{lf.end or ''}\n")
        created.append(path)
        details["features"] = path.name

    path = layout.path("_reads.tsv")
    with path.open("w", encoding="utf-8", newline="\n") as out:
        out.write("sequence\treference\tstart\tend\tstrand\tidentity\tmismatches\tdeleted\t"
                  "inserted\taligned_bp\tmapq\tkind\tlength\tunaligned_ends\n")
        for r in run.records:
            out.write(f"{fasta_id(r)}\t{r.ref_id}\t{r.start + 1}\t{r.end}\t{r.strand}\t"
                      f"{r.identity:.4f}\t{r.mismatches}\t{r.deleted}\t{r.inserted}\t"
                      f"{r.aligned_bp}\t{r.mapq if r.mapq >= 0 else ''}\t{r.kind}\t"
                      f"{r.query_length}\t{r.clipped}\n")
    created.append(path)
    details["reads"] = path.name

    if raw is not None:
        created.append(raw)
        details["raw"] = raw.name
    methods = methods_paragraph(run)
    path = layout.path("_methods.txt")
    path.write_text(methods + "\n\nReferences\n" + "\n".join(CITATIONS[run.engine]) + "\n",
                    encoding="utf-8")
    created.append(path)
    path = layout.path("_commands.txt")
    path.write_text(f"RepeatViewer Mapping, run {run.stamp}\n"
                    f"{ENGINE_NAMES[run.engine]} {run.engine_version}\n\n"
                    + "\n".join(run.commands) + "\n", encoding="utf-8")
    created.append(path)

    report = layout.path("_report.html")
    details["report"] = report.name
    write_rvmap = bool(getattr(o, "generate_rvmap", True))
    rvmap = layout.rvmap
    title = f"{output_base_name(run.queries)} mapped to {run.reference.name}"
    listed = ([rvmap] if write_rvmap else []) + [*created, report]
    report.write_text(report_html(run, listed, title), encoding="utf-8")
    created.append(report)

    if not write_rvmap:
        primary = raw if raw is not None else report
        return primary, list(created), report

    payload = {
        "format": RVMAP_FORMAT, "version": RVMAP_VERSION, "created": run.stamp,
        "title": title,
        "reference": str(run.reference), "queries": [str(q) for q in run.queries],
        "engine": run.engine, "engine_version": run.engine_version, "preset": run.preset,
        "options": {"sensitivity": o.sensitivity, "multiple_hits": o.multiple_hits,
                    "minimap2_preset": o.minimap2_preset,
                    "minimap2_output": getattr(o, "minimap2_output", "paf"),
                    "min_identity": o.min_identity, "min_aligned_bp": o.min_aligned_bp,
                    "min_mapq": o.min_mapq, "consensus_threshold": o.consensus_threshold,
                    "no_coverage": o.no_coverage},
        "summary": {"n_sequences": run.n_queries, "n_placed": len(run.placed_names),
                    "n_placements": run.n_placements, "n_reversed": len(run.reversed_names),
                    "n_unused": len(run.unused),
                    "mean_identity": round(statistics.mean([r.identity for r in run.records]), 4)
                    if run.records else None},
        "warnings": run.warnings,
        "alignment": raw.name if raw is not None else "",
        "details_dir": layout.details.name,
        "details": details,
        "unused": [[n, r] for n, r in list(run.unused.items())[:10_000]],
        "contigs": [],
    }
    budget = JSON_BASE_BUDGET
    for c in run.contigs[:500]:
        shown = []
        for item in sorted(c.display, key=lambda t: (t[0].start, t[0].end)):
            if budget <= 0:
                break
            budget -= len(item[0].row)
            shown.append(item)
        entry = {
            "id": c.ref_id, "length": c.length, "detailed": c.detailed, "n_reads": c.n_reads,
            "mean_identity": round(c.mean_identity, 4),
            "reference": c.reference, "consensus": c.consensus,
            "ref_to_cons": rle(c.ref_to_cons) if c.ref_to_cons is not None else [],
            "coverage": rle(np.rint(np.asarray(c.coverage)).astype(np.int64))
            if c.coverage is not None else [],
            "bin_size": c.bin_size,
            "features": [f.to_dict() for f in c.features[:JSON_MAX_FEATURES]],
            "feature_summary": c.feature_summary[:JSON_MAX_FEATURES],
            "variants": [[v.pos, v.ref, v.alt, v.count, v.depth, v.context, v.features]
                         for v in c.variants[:MAX_VARIANT_ROWS]],
            "reads": [_read_json(r, lf, pl, with_marks=not c.detailed) for r, lf, pl in shown],
        }
        payload["contigs"].append(entry)
    rvmap.write_text(json.dumps(payload, separators=(",", ":")), encoding="utf-8")
    return rvmap, [rvmap, *created], report


# --------------------------------------------------------------------------
# Reading a result back (the contig view)
# --------------------------------------------------------------------------


def is_rvmap_file(path: Path) -> bool:
    return Path(path).suffix.lower() == RVMAP_SUFFIX


def load_rvmap(path: Path) -> dict:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("format") != RVMAP_FORMAT:
        raise MapReferenceError(f"{Path(path).name} is not a RepeatViewer mapping result.")
    return data
