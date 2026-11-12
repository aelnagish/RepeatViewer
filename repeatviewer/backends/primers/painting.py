"""Chromosome painting probes with Chorus2 (Tools > Primers and Probes > Design).

Chorus2 (Zhang et al. 2021, Plant Biotechnol J 19:1967; MIT licence) finds
single-copy oligos along a genome. RepeatViewer either runs it here, in its
own conda env (Chorus2 needs ``pandas<2``), or imports a BED it made on a
server, then selects a pool (spacing, strand, density or pool size, regions)
and writes an order file and a report.

Coordinates: Chorus2 BED rows are ``chrom start end sequence [score strand]``
with a 1-based ``start`` and an inclusive ``end``. Inside this module probes
keep that convention; the BED we write is standard (0-based start).

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import html
import json
import math
import os
import re
import shlex
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable

from repeatviewer.backends.primers.thermo import PrimerError

SOURCE_RUN = "run"
SOURCE_IMPORT = "import"
SOURCES = (SOURCE_RUN, SOURCE_IMPORT)
STRAND_BOTH = "both"     # alternate + / - along the chromosome (Chorus2 default)
STRAND_PLUS = "plus"
STRANDS = (STRAND_BOTH, STRAND_PLUS)
TARGET_ALL = "all"
TARGET_DENSITY = "density"   # probes per kb
TARGET_POOL = "pool"         # total number of probes
TARGETS = (TARGET_ALL, TARGET_DENSITY, TARGET_POOL)
PAINTING_TAG = "painting"
DEFAULT_ENV = "chorus2"

INSTALL_COMMANDS = (
    "conda create -n chorus2 -c conda-forge -c bioconda python=3.9 \"pandas<2\" numpy "
    "cython pyfasta \"primer3-py<2\" pybigwig pybedtools kmer-jellyfish bwa",
    "conda run -n chorus2 pip install --no-deps git+https://github.com/zhangtaolab/Chorus2.git",
)
CITATION = ("Zhang T, Liu G, Zhao H, Braz GT, Jiang J (2021) Chorus2: design of genome-scale "
            "oligonucleotide-based probes for fluorescence in situ hybridization. "
            "Plant Biotechnol J 19:1967-1978.")

_COMP = str.maketrans("ACGTNacgtnRYKMSWrykmsw", "TGCANtgcanYRMKSWyrmksw")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# --------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------


@dataclass
class PaintingSettings:
    source: str = SOURCE_RUN
    genome_fasta: str = ""
    input_fasta: str = ""          # chromosome / region FASTA; "" = whole genome
    import_bed: str = ""           # Chorus2 BED made elsewhere
    conda_env: str = DEFAULT_ENV   # env name or path holding Chorus2, bwa, jellyfish
    threads: int = 4
    length: int = 45
    homology: float = 75.0
    dtm: float = 10.0
    step: int = 5
    ploidy: int = 2
    use_ngs: bool = False
    ngs_fastq: str = ""            # comma-separated FASTQ files (gz)
    ngs_kmer: int = 17
    ngs_low_quantile: float = 0.1
    ngs_high_quantile: float = 0.9
    spacing: int = 25              # minimum gap (bp) between neighbouring probes
    strand: str = STRAND_BOTH
    target: str = TARGET_ALL
    density_per_kb: float = 1.0
    pool_size: int = 20000
    regions: str = ""              # "chr1, chr2:1-5000000"; "" = everything
    primer5: str = ""              # 5' adapter added to every ordered oligo
    primer3: str = ""              # 3' adapter added to every ordered oligo
    pool_name: str = "Pool1"

    def __post_init__(self) -> None:
        if self.source not in SOURCES:
            self.source = SOURCE_RUN
        if self.strand not in STRANDS:
            self.strand = STRAND_BOTH
        if self.target not in TARGETS:
            self.target = TARGET_ALL
        self.threads = max(1, int(self.threads))
        self.length = max(20, min(100, int(self.length)))
        self.step = max(1, int(self.step))
        self.spacing = max(0, int(self.spacing))
        self.pool_size = max(1, int(self.pool_size))
        self.density_per_kb = max(0.001, float(self.density_per_kb))
        self.pool_name = (self.pool_name or "Pool1").strip() or "Pool1"

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict | None) -> PaintingSettings:
        data = dict(data or {})
        return cls(**{k: data[k] for k in cls.__dataclass_fields__ if k in data})


# --------------------------------------------------------------------------
# Probes
# --------------------------------------------------------------------------


@dataclass
class PaintingProbe:
    chrom: str
    start: int          # 1-based
    end: int            # inclusive
    sequence: str
    score: int | None = None
    strand: str = "+"

    @property
    def length(self) -> int:
        return self.end - self.start + 1


def parse_probe_bed(path: Path) -> list[PaintingProbe]:
    """Read a Chorus2 BED (4 columns, or 6 after ChorusNGSfilter / select)."""
    path = Path(path).expanduser()
    if not path.is_file():
        raise PrimerError(f"Probe BED not found: {path}")
    probes: list[PaintingProbe] = []
    with open(path, encoding="utf-8", errors="replace") as fh:
        for n, line in enumerate(fh, 1):
            line = line.strip()
            if not line or line.startswith(("#", "track", "browser")):
                continue
            cols = line.split()
            if len(cols) < 4:
                raise PrimerError(f"{path.name} line {n}: expected chrom, start, end, sequence "
                                  "(a Chorus2 BED).")
            try:
                start, end = int(cols[1]), int(cols[2])
            except ValueError:
                raise PrimerError(f"{path.name} line {n}: start / end are not numbers.") from None
            seq = cols[3].upper()
            if not re.fullmatch(r"[ACGTN]+", seq):
                raise PrimerError(f"{path.name} line {n}: column 4 is not a DNA sequence.")
            score = None
            if len(cols) >= 5:
                try:
                    score = int(float(cols[4]))
                except ValueError:
                    score = None
            strand = cols[5] if len(cols) >= 6 and cols[5] in "+-" else "+"
            probes.append(PaintingProbe(cols[0], start, end, seq, score, strand))
    if not probes:
        raise PrimerError(f"{path.name} contains no probes.")
    return probes


_REGION_RE = re.compile(r"^(?P<chrom>[^:\s]+)(?::(?P<a>[\d,_]+)\s*-\s*(?P<b>[\d,_]+))?$")


def parse_regions(text: str) -> list[tuple[str, int, int]]:
    """``"chr1, chr2:1-5,000,000"`` -> ``[("chr1", 1, 0), ("chr2", 1, 5000000)]``.

    End 0 means "to the end of the chromosome".
    """
    out = []
    # Items are split on ";", new lines, ", " and "," before a name; "5,000,000" stays whole.
    for raw in re.split(r"[;\n]|,\s+|,(?=[^\d\s])", text or ""):
        item = raw.strip().strip(",")
        if not item:
            continue
        m = _REGION_RE.match(item)
        if not m:
            raise PrimerError(f"Region not understood: {item!r}. Use chr1 or chr1:1000-2000.")
        a = int(m["a"].replace(",", "").replace("_", "")) if m["a"] else 1
        b = int(m["b"].replace(",", "").replace("_", "")) if m["b"] else 0
        if b and b < a:
            raise PrimerError(f"Region {item!r}: the end is before the start.")
        out.append((m["chrom"], max(1, a), b))
    return out


def fasta_lengths(path: Path) -> dict[str, int]:
    """Sequence lengths, from ``.fai`` when present, else by reading the FASTA."""
    path = Path(path).expanduser()
    fai = Path(str(path) + ".fai")
    lengths: dict[str, int] = {}
    if fai.is_file():
        with open(fai, encoding="utf-8") as fh:
            for line in fh:
                cols = line.split("\t")
                if len(cols) >= 2 and cols[1].strip().isdigit():
                    lengths[cols[0]] = int(cols[1])
        if lengths:
            return lengths
    import gzip

    opener = gzip.open if path.name.endswith(".gz") else open
    name, n = None, 0
    with opener(path, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    lengths[name] = n
                name, n = line[1:].split(None, 1)[0] if line[1:].strip() else "", 0
            else:
                n += len(line.strip())
    if name is not None:
        lengths[name] = n
    return lengths


def _read_len_file(path: Path) -> dict[str, int]:
    out = {}
    try:
        for line in Path(path).read_text(encoding="utf-8").splitlines():
            cols = line.split()
            if len(cols) >= 2 and cols[1].isdigit():
                out[cols[0]] = int(cols[1])
    except OSError:
        pass
    return out


def select_probes(probes: list[PaintingProbe], s: PaintingSettings,
                  lengths: dict[str, int] | None = None) -> tuple[list[PaintingProbe], dict]:
    """Choose the pool: regions, k-mer score window, spacing, density / size, strand.

    Returns the chosen probes (sorted by chromosome order of ``lengths`` or of
    the input, then position) and a dict of counts for the report.
    """
    stats: dict = {"input": len(probes)}
    regions = parse_regions(s.regions)
    if regions:
        def inside(p: PaintingProbe) -> bool:
            return any(p.chrom == c and p.start >= a and (not b or p.end <= b)
                       for c, a, b in regions)
        probes = [p for p in probes if inside(p)]
    stats["in_regions"] = len(probes)
    if s.use_ngs:
        scored = sorted(p.score for p in probes if p.score is not None)
        if scored:
            lo = _quantile(scored, s.ngs_low_quantile)
            hi = _quantile(scored, s.ngs_high_quantile)
            probes = [p for p in probes if p.score is not None and lo < p.score < hi]
            stats["score_window"] = [lo, hi]
    stats["after_score"] = len(probes)

    order = {c: i for i, c in enumerate(lengths or {})}
    for p in probes:
        order.setdefault(p.chrom, len(order))
    by_chrom: dict[str, list[PaintingProbe]] = {}
    for p in sorted(probes, key=lambda x: (order[x.chrom], x.start)):
        by_chrom.setdefault(p.chrom, []).append(p)

    # Spacing first: no two probes closer than ``spacing`` bp (gap between them).
    spaced: dict[str, list[PaintingProbe]] = {}
    for chrom, rows in by_chrom.items():
        keep, last_end = [], -10**12
        for p in rows:
            if p.start > last_end + s.spacing:
                keep.append(p)
                last_end = p.end
        spaced[chrom] = keep
    stats["after_spacing"] = sum(len(v) for v in spaced.values())

    # Then thin evenly along each chromosome to the density / pool size.
    available = stats["after_spacing"]
    if s.target == TARGET_DENSITY and available:
        step_bp = 1000.0 / s.density_per_kb
        for chrom, rows in spaced.items():
            spaced[chrom] = _thin(rows, step_bp)
    elif s.target == TARGET_POOL and available > s.pool_size:
        # Probes cluster (gaps at centromeres and repeats), so the step is searched:
        # the widest even step that still leaves at least the pool size.
        def count(step: float) -> int:
            return sum(len(_thin(rows, step)) for rows in spaced.values())

        lo, hi = 0.0, 1.0
        while count(hi) >= s.pool_size and hi < 1e9:
            lo, hi = hi, hi * 2
        for _ in range(24):
            mid = (lo + hi) / 2
            if count(mid) >= s.pool_size:
                lo = mid
            else:
                hi = mid
        for chrom, rows in spaced.items():
            spaced[chrom] = _thin(rows, lo)
        # Drop evenly from the whole pool down to the size asked for.
        flat = [p for rows in spaced.values() for p in rows]
        if len(flat) > s.pool_size:
            keep_ids = {id(flat[int(i * len(flat) / s.pool_size)]) for i in range(s.pool_size)}
            for chrom in spaced:
                spaced[chrom] = [p for p in spaced[chrom] if id(p) in keep_ids]
    chosen = [p for rows in spaced.values() for p in rows]

    if s.strand == STRAND_BOTH:
        # Chorus2 style: every second probe on a chromosome is the other strand.
        out = []
        for rows in spaced.values():
            for i, p in enumerate(rows):
                if i % 2 == 1:
                    out.append(PaintingProbe(p.chrom, p.start, p.end, revcomp(p.sequence),
                                             p.score, "-" if p.strand == "+" else "+"))
                else:
                    out.append(p)
        chosen = out
    stats["selected"] = len(chosen)
    if s.target == TARGET_POOL and len(chosen) < s.pool_size:
        stats["warning"] = (f"Only {len(chosen):,} probes are available for a pool of "
                            f"{s.pool_size:,}: lower the spacing, widen the regions or the "
                            "NGS score window, or run Chorus2 with a smaller k-mer step.")
    return chosen, stats


def _quantile(sorted_values: list[int], q: float) -> float:
    if not sorted_values:
        return 0.0
    q = min(1.0, max(0.0, q))
    pos = q * (len(sorted_values) - 1)
    lo = math.floor(pos)
    hi = min(len(sorted_values) - 1, lo + 1)
    return sorted_values[lo] + (sorted_values[hi] - sorted_values[lo]) * (pos - lo)


def _thin(rows: list[PaintingProbe], step_bp: float) -> list[PaintingProbe]:
    """Keep about one probe per ``step_bp``, as evenly spaced as the probes allow."""
    if step_bp <= 0 or not rows:
        return rows
    keep = [rows[0]]
    next_at = rows[0].start + step_bp
    for p in rows[1:]:
        if p.start >= next_at:
            keep.append(p)
            next_at = p.start + step_bp
    return keep


def density_table(probes: list[PaintingProbe], lengths: dict[str, int],
                  window: int = 1_000_000) -> dict[str, list[int]]:
    """Probes per ``window`` bp along each chromosome."""
    out: dict[str, list[int]] = {}
    for chrom, n in lengths.items():
        out[chrom] = [0] * max(1, math.ceil(n / window))
    for p in probes:
        bins = out.setdefault(p.chrom, [0])
        i = (p.start - 1) // window
        if i >= len(bins):
            bins.extend([0] * (i + 1 - len(bins)))
        bins[i] += 1
    return out


# --------------------------------------------------------------------------
# Running Chorus2
# --------------------------------------------------------------------------


@dataclass
class ChorusInstall:
    ok: bool
    message: str
    prefix: Path | None = None
    bin_dir: Path | None = None
    python: Path | None = None
    chorus2: Path | None = None
    ngs_filter: Path | None = None
    bwa: Path | None = None
    jellyfish: Path | None = None
    missing: list[str] = field(default_factory=list)


def find_chorus2(env: str = DEFAULT_ENV) -> ChorusInstall:
    """Locate the conda env with Chorus2, bwa and jellyfish."""
    if sys.platform == "win32":
        return ChorusInstall(False, "Chorus2 runs on Linux and macOS. On Windows, run it on a "
                                    "Linux server and use “Import results from the server”.")
    from repeatviewer.conda_tools import resolve_conda_env_prefix

    prefix = resolve_conda_env_prefix(env or DEFAULT_ENV)
    if prefix is None:
        return ChorusInstall(False, f"The conda env “{env or DEFAULT_ENV}” was not found.",
                             missing=["conda env"])
    bin_dir = prefix / "bin"

    def pick(*names: str) -> Path | None:
        for n in names:
            p = bin_dir / n
            if p.is_file() and os.access(p, os.X_OK):
                return p
        return None

    inst = ChorusInstall(True, "", prefix=prefix, bin_dir=bin_dir,
                         python=pick("python", "python3"),
                         chorus2=pick("Chorus2", "Chorus2.py"),
                         ngs_filter=pick("ChorusNGSfilter", "ChorusNGSfilter.py"),
                         bwa=pick("bwa"), jellyfish=pick("jellyfish"))
    for label, value in (("python", inst.python), ("Chorus2", inst.chorus2),
                         ("bwa", inst.bwa), ("jellyfish", inst.jellyfish)):
        if value is None:
            inst.missing.append(label)
    if inst.missing:
        inst.ok = False
        inst.message = (f"The env {prefix} lacks: {', '.join(inst.missing)}.")
    else:
        inst.message = f"Chorus2 found in {prefix}"
    return inst


def install_help(rich: bool = False) -> str:
    cmds = "\n".join(INSTALL_COMMANDS)
    if rich:
        return ("Install Chorus2 once, in its own conda env (it needs pandas&lt;2):<br>"
                + "<br>".join(f"<code>{html.escape(c)}</code>" for c in INSTALL_COMMANDS))
    return f"Install Chorus2 once, in its own conda env (it needs pandas<2):\n{cmds}"


def chorus_out_dir(s: PaintingSettings, details: Path) -> Path:
    return Path(details) / "chorus2"


def chorus_commands(s: PaintingSettings, out_dir: Path, inst: ChorusInstall | None = None,
                    *, server: bool = False) -> list[list[str]]:
    """Chorus2 (and ChorusNGSfilter) command lines.

    ``server=True`` gives plain ``Chorus2`` / ``ChorusNGSfilter`` names (as on
    a server where the env is active) and file names relative to the folder.
    """
    genome = Path(s.genome_fasta).expanduser()
    query = Path(s.input_fasta).expanduser() if s.input_fasta else genome
    if server:
        g, q, o = genome.name, query.name, "chorus2_out"
        pre = [["Chorus2"]]
        ngs = ["ChorusNGSfilter"]
    else:
        if inst is None or not inst.ok:
            raise PrimerError("Chorus2 is not installed.\n\n" + install_help())
        g, q, o = str(genome), str(query), str(out_dir)
        pre = [[str(inst.python), str(inst.chorus2)]]
        ngs = [str(inst.python), str(inst.ngs_filter)] if inst.ngs_filter else ["ChorusNGSfilter"]
    cmd = pre[0] + ["-g", g, "-i", q, "-s", o, "-t", str(s.threads), "-l", str(s.length),
                    "--homology", f"{s.homology:g}", "-d", f"{s.dtm:g}", "--step", str(s.step),
                    "--ploidy", str(s.ploidy), "--docker", "True"]
    if not server and inst is not None:
        cmd += ["-j", str(inst.jellyfish), "-b", str(inst.bwa)]
    cmds = [cmd]
    if s.use_ngs:
        fq = [x.strip() for x in s.ngs_fastq.split(",") if x.strip()]
        if not fq:
            raise PrimerError("NGS filter: add the FASTQ file(s) of genomic reads.")
        if server:
            fq = [Path(x).name for x in fq]
        gz = "gz" if all(x.endswith(".gz") for x in fq) else "text"
        bed = f"{o}/{q if server else Path(q).name}.bed"
        n = ngs + ["-g", g, "-i", ",".join(fq), "-z", gz, "-t", str(s.threads),
                   "-k", str(s.ngs_kmer), "-p", bed, "-o", f"{o}/probes_kmer_score.bed"]
        if not server and inst is not None:
            n += ["-j", str(inst.jellyfish)]
        cmds.append(n)
    return cmds


def _server_command_or_empty(s: PaintingSettings) -> str:
    if not s.genome_fasta:
        return ""
    try:
        return server_command_text(s)
    except PrimerError:
        return ""


def server_command_text(s: PaintingSettings) -> str:
    """Copy-paste text for running Chorus2 on a server, then importing here."""
    cmds = chorus_commands(s, Path("chorus2_out"), server=True)
    lines = ["conda activate chorus2"] + [shlex.join(c) for c in cmds]
    return "\n".join(lines)


def chorus_result_bed(s: PaintingSettings, out_dir: Path) -> Path:
    query = Path(s.input_fasta or s.genome_fasta).expanduser()
    if s.use_ngs:
        return Path(out_dir) / "probes_kmer_score.bed"
    return Path(out_dir) / f"{query.name}.bed"


def run_chorus2(s: PaintingSettings, out_dir: Path, *, cancel=None,
                log: Callable[[str], None] | None = None) -> Path:
    """Run Chorus2 (and the NGS filter); return the BED to select from."""
    from repeatviewer.external_process import run_logged_process

    inst = find_chorus2(s.conda_env)
    if not inst.ok:
        raise PrimerError(inst.message + "\n\n" + install_help())
    for label, value in (("Genome", s.genome_fasta), ("Input", s.input_fasta)):
        if value and not Path(value).expanduser().is_file():
            raise PrimerError(f"{label} FASTA not found: {value}")
    if not s.genome_fasta:
        raise PrimerError("Choose the genome FASTA (Chorus2 tab).")
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ)
    env["PATH"] = f"{inst.bin_dir}{os.pathsep}{env.get('PATH', '')}"
    env.pop("PYTHONPATH", None)
    env.pop("PYTHONHOME", None)
    for i, cmd in enumerate(chorus_commands(s, out_dir, inst)):
        tool = "Chorus2" if i == 0 else "ChorusNGSfilter"
        result = run_logged_process(cmd, tool=tool, cwd=out_dir, env=env, cancel=cancel)
        if cancel is not None and getattr(cancel, "requested", False):
            from repeatviewer.backends.primers.binding import PrimerCancelled

            raise PrimerCancelled()
        if result.returncode != 0:
            tail = "\n".join((result.stderr or result.stdout or "").strip().splitlines()[-12:])
            raise PrimerError(f"{tool} stopped (exit code {result.returncode}).\n\n{tail}")
    bed = chorus_result_bed(s, out_dir)
    if not bed.is_file():
        raise PrimerError(f"Chorus2 finished but {bed.name} was not written. See the process "
                          "window for its messages.")
    return bed


# --------------------------------------------------------------------------
# Job: run or import, select, write
# --------------------------------------------------------------------------


def run_painting(s: PaintingSettings, *, report_path: Path | None = None, cancel=None,
                 progress: Callable[[str], None] | None = None) -> tuple[Path, list[Path], dict]:
    from repeatviewer.backends.primers import jobs as pj

    def say(msg: str) -> None:
        if progress:
            progress(msg)

    if s.source == SOURCE_IMPORT:
        if not s.import_bed:
            raise PrimerError("Choose the Chorus2 BED to import (Chorus2 tab).")
        anchor = Path(s.import_bed).expanduser()
        if not anchor.is_file():
            raise PrimerError(f"Probe BED not found: {anchor}")
    else:
        if not s.genome_fasta:
            raise PrimerError("Choose the genome FASTA (Chorus2 tab).")
        anchor = Path(s.genome_fasta).expanduser()
        if not anchor.is_file():
            raise PrimerError(f"Genome FASTA not found: {anchor}")
        inst = find_chorus2(s.conda_env)
        if not inst.ok:
            raise PrimerError(inst.message + "\n\n" + install_help())
    layout = pj.plan_outputs(anchor, PAINTING_TAG, report_path=report_path)
    try:
        return _run_painting_into(s, layout, anchor, cancel, say)
    except BaseException:
        _remove_if_empty(layout.details)
        raise


def _remove_if_empty(folder: Path) -> None:
    try:
        if folder.is_dir() and not any(folder.iterdir()):
            folder.rmdir()
    except OSError:
        pass


def _run_painting_into(s: PaintingSettings, layout, anchor: Path, cancel, say):
    started = time.time()
    if s.source == SOURCE_IMPORT:
        bed = anchor
        say("Reading the Chorus2 probes…")
    else:
        say("Running Chorus2 (this takes a while on a large genome)…")
        bed = run_chorus2(s, chorus_out_dir(s, layout.details), cancel=cancel)
        say("Reading the Chorus2 probes…")
    probes = parse_probe_bed(bed)
    lengths: dict[str, int] = {}
    # Chorus2 writes <input>.len next to its BED: quicker than reading a large genome.
    for cand in (bed.with_suffix(".len"), bed.parent / (bed.stem + ".len")):
        if cand.is_file():
            lengths = _read_len_file(cand)
            if lengths:
                break
    if not lengths and s.genome_fasta and Path(s.genome_fasta).expanduser().is_file():
        say("Reading chromosome lengths…")
        lengths = fasta_lengths(Path(s.genome_fasta))
    if not lengths:
        for p in probes:
            lengths[p.chrom] = max(lengths.get(p.chrom, 0), p.end)
    say("Selecting the pool…")
    chosen, stats = select_probes(probes, s, lengths)
    if not chosen:
        raise PrimerError("No probe is left after selection. Check Regions, the NGS score window "
                          "and the spacing (Selection tab).")
    say("Writing results…")
    files = write_painting_outputs(layout, chosen, probes, stats, s, lengths, bed,
                                   time.time() - started)
    payload = {"kind": "painting", "selected": len(chosen), "input": len(probes),
               "stats": stats, "report": str(files[0]),
               "warnings": [stats["warning"]] if stats.get("warning") else []}
    return files[0], files, payload


def order_name(s: PaintingSettings, i: int, p: PaintingProbe) -> str:
    return f"{s.pool_name}_{p.chrom}_{p.start}{'r' if p.strand == '-' else ''}_{i:06d}"


def write_painting_outputs(layout, chosen: list[PaintingProbe], all_probes: list[PaintingProbe],
                           stats: dict, s: PaintingSettings, lengths: dict[str, int],
                           source_bed: Path, seconds: float) -> list[Path]:
    from repeatviewer.backends.primers import jobs as pj

    d = layout.details
    base = layout.base
    bed = d / f"{base}_probes.bed"
    fa = d / f"{base}_probes.fa"
    order = d / f"{base}_order.tsv"
    js = d / f"{base}.json"
    with open(bed, "w", encoding="utf-8") as fh:
        fh.write(f'track name="{s.pool_name}" description="Chorus2 painting probes"\n')
        for i, p in enumerate(chosen, 1):
            fh.write(f"{p.chrom}\t{p.start - 1}\t{p.end}\t{order_name(s, i, p)}\t"
                     f"{p.score if p.score is not None else 0}\t{p.strand}\n")
    with open(fa, "w", encoding="utf-8") as fh:
        for i, p in enumerate(chosen, 1):
            fh.write(f">{order_name(s, i, p)} {p.chrom}:{p.start}-{p.end}({p.strand})\n"
                     f"{p.sequence}\n")
    with open(order, "w", encoding="utf-8") as fh:
        fh.write("pool\tname\tsequence\tlength\tchrom\tstart\tend\tstrand\tscore\n")
        for i, p in enumerate(chosen, 1):
            full = f"{s.primer5.strip()}{p.sequence}{s.primer3.strip()}"
            fh.write(f"{s.pool_name}\t{order_name(s, i, p)}\t{full}\t{len(full)}\t{p.chrom}\t"
                     f"{p.start}\t{p.end}\t{p.strand}\t"
                     f"{p.score if p.score is not None else ''}\n")
    window = _density_window(lengths)
    dens = density_table(chosen, lengths, window)
    per_chrom = []
    for chrom, n in lengths.items():
        k = sum(1 for p in chosen if p.chrom == chrom)
        per_chrom.append({"chrom": chrom, "length": n, "probes": k,
                          "per_kb": round(k / (n / 1000), 3) if n else 0})
    data = {
        "format": pj.RESULT_FORMAT, "kind": "painting",
        "created": time.strftime("%Y-%m-%d %H:%M"),
        "settings": s.to_dict(), "source_bed": str(source_bed), "stats": stats,
        "per_chromosome": per_chrom, "density_window": window, "density": dens,
        "seconds": round(seconds, 1),
        "probe_length": (sorted(p.length for p in chosen)[len(chosen) // 2] if chosen else 0),
        "server_command": _server_command_or_empty(s),
    }
    js.write_text(json.dumps(data, indent=1), encoding="utf-8")
    report = layout.parent / f"{base}{pj.REPORT_SUFFIX}"
    report.write_text(painting_report_html(data, len(all_probes), [bed, fa, order]),
                      encoding="utf-8")
    return [report, bed, fa, order, js]


def _density_window(lengths: dict[str, int]) -> int:
    longest = max(lengths.values(), default=0)
    for w in (10_000, 100_000, 500_000, 1_000_000, 5_000_000):
        if longest / w <= 400:
            return w
    return 10_000_000


def _fmt_bp(n: int) -> str:
    if n >= 1_000_000:
        return f"{n / 1e6:.1f} Mb"
    if n >= 1000:
        return f"{n / 1e3:.1f} kb"
    return f"{n} bp"


def _density_svg(data: dict) -> str:
    dens = data["density"]
    window = data["density_window"]
    chroms = [c for c in dens if any(dens[c])] or list(dens)[:40]
    if not chroms:
        return ""
    chroms = chroms[:60]
    peak = max((max(v) for c, v in dens.items() if c in chroms and v), default=1) or 1
    longest = max(len(dens[c]) for c in chroms)
    label_w, row_h, plot_w = 110, 22, 600
    h = row_h * len(chroms) + 24
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {label_w + plot_w + 10} {h}" '
             f'width="100%" role="img" aria-label="Probe density per chromosome">']
    bw = plot_w / max(1, longest)
    for r, c in enumerate(chroms):
        y = r * row_h + 4
        parts.append(f'<text x="{label_w - 8}" y="{y + 13}" text-anchor="end" '
                     f'class="lab">{html.escape(c[:16])}</text>')
        n = len(dens[c])
        parts.append(f'<rect x="{label_w}" y="{y + 2}" width="{bw * n:.1f}" height="{row_h - 6}" '
                     f'rx="3" class="chr"/>')
        for i, v in enumerate(dens[c]):
            if v:
                hh = (row_h - 6) * v / peak
                parts.append(f'<rect x="{label_w + i * bw:.2f}" y="{y + 2 + (row_h - 6) - hh:.1f}" '
                             f'width="{max(bw, 0.6):.2f}" height="{hh:.1f}" class="bar"/>')
    parts.append(f'<text x="{label_w}" y="{h - 4}" class="lab">bar = probes per '
                 f'{_fmt_bp(window)} (max {peak})</text></svg>')
    return "".join(parts)


def painting_report_html(data: dict, n_input: int, files: list[Path]) -> str:
    e = html.escape
    s = data["settings"]
    st = data["stats"]
    rows = "".join(
        f"<tr><td>{e(r['chrom'])}</td><td>{_fmt_bp(r['length'])}</td><td>{r['probes']:,}</td>"
        f"<td>{r['per_kb']}</td></tr>" for r in data["per_chromosome"] if r["probes"] or len(
            data["per_chromosome"]) <= 60)
    target = {TARGET_ALL: "all that pass the spacing",
              TARGET_DENSITY: f"{s['density_per_kb']:g} probes per kb",
              TARGET_POOL: f"a pool of {s['pool_size']:,}"}[s["target"]]
    steps = [("Probes from Chorus2", n_input), ("In the regions", st.get("in_regions", n_input))]
    if s["use_ngs"]:
        steps.append(("Inside the NGS k-mer score window", st.get("after_score", 0)))
    steps += [(f"After spacing ≥ {s['spacing']} bp", st.get("after_spacing", 0)),
              ("Selected", st.get("selected", 0))]
    funnel = "".join(f"<tr><td>{e(a)}</td><td>{b:,}</td></tr>" for a, b in steps)
    imported = s["source"] == SOURCE_IMPORT
    length = data.get("probe_length") or s["length"]
    methods = (
        f"Oligo probes for chromosome painting were designed with Chorus2 v2.1 (Zhang et al. "
        f"2021){' on a server' if imported else ''}: {length} nt, at most {s['homology']:g} % homology to other genome "
        f"positions, dTm ≥ {s['dtm']:g} °C, {s['step']} bp k-mer step, ploidy {s['ploidy']}"
        + (f"; probes were scored with {s['ngs_kmer']}-mer counts from genomic reads "
           f"(ChorusNGSfilter) and kept between the {s['ngs_low_quantile']:.0%} and "
           f"{s['ngs_high_quantile']:.0%} score quantiles" if s["use_ngs"] else "")
        + f". From these, {st.get('selected', 0):,} probes were selected in RepeatViewer "
        f"({target}, at least {s['spacing']} bp apart"
        + (", alternating strands" if s["strand"] == STRAND_BOTH else ", plus strand") + ")."
    )
    files_html = "".join(f"<li><code>{e(Path(f).parent.name)}/{e(Path(f).name)}</code></li>"
                         for f in files)
    server = data.get("server_command") or ""
    server_html = (f"<h2>Run on a server</h2><p>Copy, run in the folder with the genome, then "
                   f"use <b>Import results from the server</b> with "
                   f"<code>chorus2_out/…bed</code>.</p><pre>{e(server)}</pre>") if server else ""
    return f"""<!doctype html><html><head><meta charset="utf-8">
<title>{e(s['pool_name'])}: chromosome painting probes</title>
<style>
:root{{--bg:#fff;--fg:#1f2328;--muted:#5f6670;--line:#d9dce0;--soft:#f3f4f6;--bar:#2f6fb2;--chr:#e6e9ee}}
@media (prefers-color-scheme: dark){{:root{{--bg:#1c1d1f;--fg:#e6e6e6;--muted:#a3a8ae;--line:#3a3d42;--soft:#26282b;--bar:#7fb2e8;--chr:#30343a}}}}
body{{margin:0;background:var(--bg);color:var(--fg);font:14px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Roboto,sans-serif}}
.wrap{{max-width:980px;margin:0 auto;padding:22px 18px 50px}}
h1{{font-size:22px;margin:0 0 4px}} h2{{font-size:17px;margin:28px 0 8px}}
.sub{{color:var(--muted)}}
.tiles{{display:flex;gap:10px;flex-wrap:wrap;margin:14px 0}}
.tile{{background:var(--soft);border-radius:8px;padding:10px 14px;min-width:130px}}
.tile b{{display:block;font-size:20px}}
table{{border-collapse:collapse;font-size:13px}} td,th{{padding:4px 10px;border-bottom:1px solid var(--line);text-align:left;white-space:nowrap}}
th{{color:var(--muted)}}
pre{{background:var(--soft);padding:10px;border-radius:6px;overflow-x:auto;font-size:12px}}
.warn{{border-left:3px solid #b45309;background:var(--soft);padding:8px 12px;border-radius:0 6px 6px 0}}
.methods{{background:var(--soft);border-radius:8px;padding:12px 14px;line-height:1.55}}
svg .lab{{font-size:11px;fill:var(--muted)}} svg .chr{{fill:var(--chr)}} svg .bar{{fill:var(--bar)}}
</style></head><body><div class="wrap">
<h1>Chromosome painting probes</h1>
<div class="sub">Pool {e(s['pool_name'])} · {e(data['created'])} · Chorus2 BED: {e(Path(data['source_bed']).name)}</div>
<div class="tiles">
<div class="tile"><b>{st.get('selected', 0):,}</b>probes selected</div>
<div class="tile"><b>{n_input:,}</b>from Chorus2</div>
<div class="tile"><b>{data.get('probe_length') or s['length']} nt</b>probe length</div>
<div class="tile"><b>{len([r for r in data['per_chromosome'] if r['probes']])}</b>chromosomes covered</div>
</div>
{f'<p class="warn">{e(st["warning"])}</p>' if st.get("warning") else ""}
<h2>Density along the chromosomes</h2>
{_density_svg(data)}
<h2>Selection</h2><table>{funnel}</table>
<h2>Per chromosome</h2>
<table><tr><th>Chromosome</th><th>Length</th><th>Probes</th><th>Probes per kb</th></tr>{rows}</table>
<h2>Files</h2><ul>{files_html}</ul>
<p class="sub">The order file holds the full oligo (5′ adapter + probe + 3′ adapter). The BED
is standard (0-based start), for IGV or any genome browser.</p>
<h2>Methods text</h2><div class="methods">{e(methods)}</div>
<p class="sub">{e(CITATION)}</p>
{server_html}
</div></body></html>"""
