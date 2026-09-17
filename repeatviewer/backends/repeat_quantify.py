"""Quantify a tandem-repeat / satellite library across target sequences.

Given a library of known monomers/families and a target (scaffolds or a whole
assembly), report where each family occurs and how much of the target it
covers. Two engines:

- **BLASTn (light):** library-vs-target search, aggregated per family into copy
  count, merged-array count, covered bp, % of target, and mean identity, plus a
  GFF track of the arrays. Reuses the bundled BLAST binaries.
- **RepeatMasker (standard):** ``-lib library`` run for the genome-wide landscape
  (added in a later step; see ``quantify_library_repeatmasker``).

This is *quantification against a known library* — not de-novo discovery.
"""

from __future__ import annotations

import hashlib
import heapq
import json
import logging
import math
import os
import platform
import shutil
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from Bio import SeqIO

from repeatviewer.app_preferences import cache_dir
from repeatviewer.backends.gui_functions import (
    blast_database_basename,
    check_database,
    get_original_file_path,
)
from repeatviewer.backends.job_cancel import run_subprocess_with_cancel

logger = logging.getLogger(__name__)

# blastn tabular fields we request (outfmt 6 + qlen/slen).
_BLAST_FIELDS = (
    "qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen"
)


def _blastn_executable(use_system_blast: bool = False) -> str:
    """Resolve the blastn binary (bundled per-OS, or system PATH)."""
    if use_system_blast:
        import shutil

        found = shutil.which("blastn")
        if found:
            return found
    blast_dir = os.path.join(get_original_file_path(), "blast")
    system = platform.system()
    if system == "Linux":
        return os.path.join(blast_dir, "blastn_linux")
    if system == "Darwin":
        return os.path.join(blast_dir, "blastn_mac")
    return os.path.join(blast_dir, "blastn.exe")


@dataclass
class FamilyQuant:
    """Per-family quantification against the target."""

    family: str
    hits: int = 0
    arrays: int = 0            # merged overlapping/adjacent hit intervals
    covered_bp: int = 0        # bp of target covered (merged, non-double-counted)
    pct_of_target: float = 0.0
    mean_identity: float = 0.0
    # RepeatMasker only (0 for blastn): CpG-adjusted Kimura when the run has
    # a .align file, else the raw .out value -- see QuantifyResult.divergence_metric.
    mean_divergence: float = 0.0
    scaffolds: int = 0         # distinct target sequences hit
    monomer_len: int = 0       # library monomer length (period size)
    copy_number: float = 0.0   # covered_bp / monomer_len (fractional below 1)
    mean_array_bp: float = 0.0  # covered_bp / arrays
    max_array_bp: int = 0      # largest single contiguous array


@dataclass
class QuantifyResult:
    """Outcome of a library-quantification run."""

    engine: str
    target_length: int
    families: list[FamilyQuant] = field(default_factory=list)
    summary_tsv: Path | None = None
    gff: Path | None = None
    # RepeatMasker landscape: (family, divergence_bin_start, bp) rows.
    landscape: list[tuple[str, float, int]] = field(default_factory=list)
    landscape_tsv: Path | None = None
    landscape_bin_width: float = 1.0
    # "kimura": CpG-adjusted Kimura from the RepeatMasker .align file;
    # "raw": the .out "perc div." column (no .align available).
    landscape_metric: str = "raw"
    # Hits with no matching .align record (Jukes-Cantor-corrected instead).
    landscape_unmatched_hits: int = 0
    rm_out: Path | None = None

    @property
    def divergence_metric(self) -> str:
        """Metric behind ``mean_divergence`` and the landscape ("" for blastn)."""
        return self.landscape_metric if self.engine == "repeatmasker" else ""

    @property
    def total_covered_bp(self) -> int:
        return sum(f.covered_bp for f in self.families)


def _format_copy_number(value: float, *, thousands: bool = False) -> str:
    """Render copy number: whole copies as an integer, sub-copy as a fraction.

    A family covering less than one monomer -- e.g. a partial hit against the
    18 kb 35S_rDNA locus consensus -- used to round to a bare ``0``, which
    reads like "not found" rather than "a fifth of one copy". Keep two
    decimals below 1 and round to whole copies at or above it.
    """
    if 0 < value < 1:
        return f"{value:.2f}".rstrip("0").rstrip(".")
    return f"{value:,.0f}" if thousands else f"{value:.0f}"


def _fasta_total_length(fasta: Path) -> int:
    return sum(len(rec.seq) for rec in SeqIO.parse(str(fasta), "fasta"))


def _fasta_lengths(fasta: Path) -> dict[str, int]:
    """Map each record id to its length (used for satellite copy number).

    Library headers carry the class as a suffix (``TRC_7#Simple_repeat/TRC_7``)
    and blastn reports that whole string back as ``sseqid``. RepeatMasker does
    not: its ``.out`` splits the header on ``#`` and reports the bare name in
    the "matching repeat" column, so a full-id lookup misses and monomer_len /
    copy_number come out 0. Register the pre-``#`` stem as an extra key so both
    engines resolve.

    A stem shared by two different library entries is ambiguous and is left
    out rather than resolved arbitrarily; an exact id always wins over a stem.
    """
    lengths: dict[str, int] = {}
    stem_ids: dict[str, set[str]] = {}
    stem_lengths: dict[str, int] = {}
    for rec in SeqIO.parse(str(fasta), "fasta"):
        # Last duplicate wins, matching what library_duplicate_warnings() reports.
        lengths[rec.id] = len(rec.seq)
        stem = rec.id.split("#", 1)[0]
        if stem == rec.id:
            continue
        stem_ids.setdefault(stem, set()).add(rec.id)
        stem_lengths[stem] = len(rec.seq)
    for stem, ids in stem_ids.items():
        if len(ids) == 1 and stem not in lengths:
            lengths[stem] = stem_lengths[stem]
    return lengths


def library_duplicate_warnings(fasta: Path) -> list[str]:
    """Warn when the library has repeated record ids (they silently merge)."""
    lengths: dict[str, list[int]] = {}
    for rec in SeqIO.parse(str(fasta), "fasta"):
        lengths.setdefault(rec.id, []).append(len(rec.seq))
    warnings: list[str] = []
    for name, lens in lengths.items():
        if len(lens) > 1:
            shown = ", ".join(f"{length:,}" for length in lens)
            warnings.append(
                f"'{name}' appears {len(lens)}× in the library ({shown} bp) — "
                f"hits are merged into one row and copy_number uses {lens[-1]:,} bp. "
                "Give each entry a unique name to keep them separate."
            )
    return warnings


def _array_stats(by_scaffold: dict) -> tuple[int, int, int]:
    """Return (arrays, covered_bp, max_array_bp) from per-scaffold intervals."""
    arrays = 0
    covered = 0
    max_array = 0
    for _scaf, ivals in by_scaffold.items():
        for start, end in _merge_intervals(ivals):
            arrays += 1
            length = end - start + 1
            covered += length
            if length > max_array:
                max_array = length
    return arrays, covered, max_array


def _merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    """Merge overlapping/adjacent [start, end] intervals (1-based inclusive)."""
    if not intervals:
        return []
    ordered = sorted((min(a, b), max(a, b)) for a, b in intervals)
    merged = [ordered[0]]
    for start, end in ordered[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end + 1:  # overlap or directly adjacent
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def _aggregate_blast_tab(
    tab_path: Path,
    target_length: int,
    *,
    monomer_lengths: dict[str, int] | None = None,
) -> list[FamilyQuant]:
    """Aggregate a blastn outfmt-6 file into per-family quantification."""
    monomer_lengths = monomer_lengths or {}
    # family -> {"pident_len": [(pident, length), ...],
    #            "by_scaffold": {sseqid: [(sstart, send), ...]}}
    acc: dict[str, dict] = {}
    with open(tab_path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 8:
                continue
            qseqid, sseqid = parts[0], parts[1]
            try:
                pident = float(parts[2])
                length = int(parts[3])
                sstart = int(parts[6])
                send = int(parts[7])
            except ValueError:
                continue
            fam = acc.setdefault(
                qseqid, {"pident_len": [], "by_scaffold": {}}
            )
            fam["pident_len"].append((pident, length))
            fam["by_scaffold"].setdefault(sseqid, []).append((sstart, send))

    families: list[FamilyQuant] = []
    for fam_name, data in acc.items():
        hits = len(data["pident_len"])
        arrays, covered, max_array = _array_stats(data["by_scaffold"])
        total_len = sum(length for _pi, length in data["pident_len"]) or 1
        mean_id = sum(pi * length for pi, length in data["pident_len"]) / total_len
        pct = (covered / target_length * 100.0) if target_length else 0.0
        mono = int(monomer_lengths.get(fam_name, 0))
        families.append(
            FamilyQuant(
                family=fam_name,
                hits=hits,
                arrays=arrays,
                covered_bp=covered,
                pct_of_target=round(pct, 4),
                mean_identity=round(mean_id, 2),
                scaffolds=len(data["by_scaffold"]),
                monomer_len=mono,
                copy_number=round(covered / mono, 2) if mono else 0.0,
                mean_array_bp=round(covered / arrays, 1) if arrays else 0.0,
                max_array_bp=max_array,
            )
        )
    families.sort(key=lambda f: f.covered_bp, reverse=True)
    return families


def _write_summary_tsv(result: QuantifyResult, path: Path) -> Path:
    lines = [
        "family\tmonomer_len\tcopy_number\thits\tarrays\tcovered_bp\t"
        "pct_of_target\tmean_array_bp\tmax_array_bp\tmean_identity\t"
        "mean_divergence\tscaffolds\tdivergence_metric"
    ]
    metric = _LANDSCAPE_METRIC_LABELS.get(
        result.divergence_metric, result.divergence_metric
    )
    for f in result.families:
        lines.append(
            f"{f.family}\t{f.monomer_len}\t{_format_copy_number(f.copy_number)}\t{f.hits}\t"
            f"{f.arrays}\t{f.covered_bp}\t{f.pct_of_target}\t"
            f"{round(f.mean_array_bp)}\t{f.max_array_bp}\t{f.mean_identity}\t"
            f"{f.mean_divergence}\t{f.scaffolds}\t{metric}"
        )
    path.write_text("\n".join(lines) + "\n")
    return path


def _write_arrays_gff(tab_path: Path, path: Path) -> Path:
    """Write a GFF3 of merged per-(family, scaffold) arrays for track viewing."""
    # (family, scaffold) -> list of (start, end, best_pident)
    groups: dict[tuple[str, str], list[tuple[int, int, float]]] = {}
    with open(tab_path) as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            try:
                pident = float(parts[2])
                sstart, send = int(parts[6]), int(parts[7])
            except ValueError:
                continue
            groups.setdefault((parts[0], parts[1]), []).append(
                (min(sstart, send), max(sstart, send), pident)
            )

    with open(path, "w") as out:
        out.write("##gff-version 3\n")
        for (family, scaffold), ivals in groups.items():
            merged = _merge_intervals([(s, e) for s, e, _p in ivals])
            for start, end in merged:
                score = max(
                    (p for s, e, p in ivals if s <= end and e >= start),
                    default=0.0,
                )
                out.write(
                    f"{scaffold}\tRepeatViewer\ttandem_repeat\t{start}\t{end}\t"
                    f"{score:.1f}\t.\t.\tName={family}\n"
                )
    return path


def _blast_db_cache_dir(target_fasta: Path) -> Path:
    """Per-genome folder for the BLAST index, keyed by path, size and mtime.

    Indexing a large genome is the slowest part of a BLASTn quantify run and the
    result is identical every time, so it is kept instead of rebuilt. The key
    includes size and mtime, so an edited genome lands in a new folder rather
    than silently reusing a stale index. Lives under the app cache directory,
    which "Clear cache" in Settings already removes.
    """
    target_fasta = Path(target_fasta).expanduser().resolve()
    try:
        stat = target_fasta.stat()
        fingerprint = f"{target_fasta}|{stat.st_size}|{stat.st_mtime_ns}"
    except OSError:
        fingerprint = str(target_fasta)
    digest = hashlib.sha1(fingerprint.encode("utf-8", "replace")).hexdigest()[:16]
    safe_stem = blast_database_basename(target_fasta.name)
    return cache_dir() / "blast_db" / f"{safe_stem}-{digest}"


def _blast_index_present(db_path: str) -> bool:
    """True when makeblastdb output for ``db_path`` is already on disk.

    ``.nin`` is a single-volume index; ``.nal`` is the alias file written for
    genomes large enough to be split across volumes.
    """
    return os.path.isfile(f"{db_path}.nin") or os.path.isfile(f"{db_path}.nal")


def quantify_library_blastn(
    library_fasta: Path,
    target_fasta: Path,
    out_dir: Path,
    *,
    e_value: float = 1e-5,
    min_identity: float = 0.0,
    min_coverage: float = 0.0,
    word_size: int = 11,
    threads: int = 1,
    use_system_blast: bool = False,
    cancel=None,
) -> QuantifyResult:
    """Quantify ``library_fasta`` occurrences across ``target_fasta`` via blastn.

    Builds a nucleotide DB from the target (cached between runs, see
    ``_blast_db_cache_dir``), searches the library against it, then aggregates
    per family. Returns a QuantifyResult and writes a summary TSV and an arrays
    GFF into ``out_dir``.
    """
    library_fasta = Path(library_fasta).expanduser().resolve()
    target_fasta = Path(target_fasta).expanduser().resolve()
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    for path in (library_fasta, target_fasta):
        if not path.is_file():
            raise FileNotFoundError(f"FASTA not found: {path}")

    os_type = platform.system()
    # The index lives in the app cache, not the results folder: it keeps the ~8
    # index files out of the way (the reason this used to be a scratch dir) and
    # lets repeat runs against the same genome skip makeblastdb entirely.
    db_dir = _blast_db_cache_dir(target_fasta)
    db_dir.mkdir(parents=True, exist_ok=True)
    db_path = os.path.join(str(db_dir), blast_database_basename(str(target_fasta)))
    reused = _blast_index_present(db_path)
    logger.info(
        "repeat quantify blastn: %s BLAST index for %s (%s)",
        "reusing cached" if reused else "building",
        target_fasta.name,
        db_dir,
    )
    tab_path = out_dir / f"{target_fasta.stem}_vs_{library_fasta.stem}.blastn.tsv"
    db = check_database(
        str(target_fasta),
        str(db_dir),
        os_type=os_type,
        use_system_blast=use_system_blast,
    )
    if db in {"makeblastdb_not_found", "makeblastdb_got_error"}:
        raise RuntimeError(
            "Could not build a BLAST database for the target sequences "
            f"({db}). Check that makeblastdb is available."
        )

    cmd = [
        _blastn_executable(use_system_blast),
        "-task", "blastn",
        "-query", str(library_fasta),
        "-db", db_path,
        "-evalue", str(e_value),
        "-word_size", str(int(word_size)),
        "-num_threads", str(max(1, int(threads))),
        "-dust", "no",
        "-outfmt", f"6 {_BLAST_FIELDS}",
        "-out", str(tab_path),
    ]
    if min_identity and min_identity > 0:
        cmd.extend(["-perc_identity", str(float(min_identity))])
    if min_coverage and min_coverage > 0:
        cmd.extend(["-qcov_hsp_perc", str(float(min_coverage))])
    # Cancellable runner so Stop actually terminates BLASTn on big inputs.
    run_subprocess_with_cancel(cmd, cancel=cancel, check=True)

    target_length = _fasta_total_length(target_fasta)
    result = QuantifyResult(
        engine="blastn",
        target_length=target_length,
        families=_aggregate_blast_tab(
            tab_path, target_length, monomer_lengths=_fasta_lengths(library_fasta)
        ),
    )
    result.summary_tsv = _write_summary_tsv(
        result, out_dir / f"{target_fasta.stem}_repeat_quant.tsv"
    )
    result.gff = _write_arrays_gff(
        tab_path, out_dir / f"{target_fasta.stem}_repeat_quant.gff"
    )
    return result


# ---------------------------------------------------------------------------
# RepeatMasker engine (standard: -lib library, genome-wide + divergence landscape)
# ---------------------------------------------------------------------------

# RepeatMasker .out whitespace columns (0-based): see REPEATMASKER_OUT_COLUMNS
# in blast_genome_plot.py. We only need divergence, query coords, and the
# matching repeat (library family) name.
# "perc div." is raw percent substitution (mismatched / aligned bases, no
# indels, no multiple-hit correction; rmhubley/RepeatMasker#216). Kimura
# divergence exists only in the -a .align file -- see
# parse_repeatmasker_align_kimura().
_RM_PERC_DIV = 1
_RM_QUERY_NAME = 4
_RM_QUERY_START = 5
_RM_QUERY_END = 6
_RM_REPEAT_NAME = 9


def _parse_repeatmasker_out(out_path: Path) -> list[dict]:
    """Parse a RepeatMasker ``.out`` into per-hit dicts on the target."""
    rows: list[dict] = []
    with open(out_path) as handle:
        for line in handle:
            parts = line.split()
            if len(parts) < 10:
                continue
            try:
                div = float(parts[_RM_PERC_DIV])
                qstart = int(parts[_RM_QUERY_START])
                qend = int(parts[_RM_QUERY_END])
            except ValueError:
                continue  # header / non-data lines
            rows.append(
                {
                    "family": parts[_RM_REPEAT_NAME],
                    "scaffold": parts[_RM_QUERY_NAME],
                    "start": min(qstart, qend),
                    "end": max(qstart, qend),
                    "divergence": div,
                }
            )
    return rows


def _aggregate_rm_rows(
    rows: list[dict],
    target_length: int,
    *,
    monomer_lengths: dict[str, int] | None = None,
) -> list[FamilyQuant]:
    monomer_lengths = monomer_lengths or {}
    acc: dict[str, dict] = {}
    for r in rows:
        fam = acc.setdefault(
            r["family"], {"div_len": [], "by_scaffold": {}}
        )
        length = r["end"] - r["start"] + 1
        # (raw perc div, reported divergence, length). Identity stays
        # 100 - raw substitution %; divergence uses Kimura when present.
        fam["div_len"].append(
            (r["divergence"], r.get("landscape_div", r["divergence"]), length)
        )
        fam["by_scaffold"].setdefault(r["scaffold"], []).append((r["start"], r["end"]))

    families: list[FamilyQuant] = []
    for fam_name, data in acc.items():
        arrays, covered, max_array = _array_stats(data["by_scaffold"])
        total_len = sum(length for _r, _d, length in data["div_len"]) or 1
        mean_raw = sum(r * length for r, _d, length in data["div_len"]) / total_len
        mean_div = sum(d * length for _r, d, length in data["div_len"]) / total_len
        pct = (covered / target_length * 100.0) if target_length else 0.0
        mono = int(monomer_lengths.get(fam_name, 0))
        families.append(
            FamilyQuant(
                family=fam_name,
                hits=len(data["div_len"]),
                arrays=arrays,
                covered_bp=covered,
                pct_of_target=round(pct, 4),
                mean_identity=round(100.0 - mean_raw, 2),
                mean_divergence=round(mean_div, 2),
                scaffolds=len(data["by_scaffold"]),
                monomer_len=mono,
                copy_number=round(covered / mono, 2) if mono else 0.0,
                mean_array_bp=round(covered / arrays, 1) if arrays else 0.0,
                max_array_bp=max_array,
            )
        )
    families.sort(key=lambda f: f.covered_bp, reverse=True)
    return families


# Trailer line RepeatMasker writes under each alignment in the ``-a`` file.
_ALIGN_KIMURA_PREFIX = "Kimura (with divCpGMod)"


def align_path_for_out(out_path: Path) -> Path | None:
    """Return the ``.align`` written next to a RepeatMasker ``.out`` (if any)."""
    out_path = Path(out_path)
    if out_path.name.endswith(".out"):
        candidate = out_path.with_name(out_path.name[: -len(".out")] + ".align")
        if candidate.is_file():
            return candidate
    return None


def _align_header_key(parts: list[str]) -> tuple[str, int, int, str] | None:
    """Key an ``.align`` header line the same way as a parsed ``.out`` row.

    Header layout (whitespace separated)::

        score div del ins query qstart qend (qleft) [C] repeat[#class] ...

    Returns ``(query, start, end, bare_repeat_name)`` or ``None`` for any
    other line.
    """
    if len(parts) < 9:
        return None
    try:
        int(parts[0])
        float(parts[1])
        qstart = int(parts[5])
        qend = int(parts[6])
    except ValueError:
        return None
    if not parts[7].startswith("("):
        return None
    idx = 9 if parts[8] == "C" else 8
    if idx >= len(parts):
        return None
    repeat = parts[idx].split("#", 1)[0]
    return (parts[4], min(qstart, qend), max(qstart, qend), repeat)


def parse_repeatmasker_align_kimura(
    align_path: Path,
) -> dict[tuple[str, int, int, str], float]:
    """Map each alignment to its CpG-adjusted Kimura divergence (percent).

    RepeatMasker (``-a``) writes ``Kimura (with divCpGMod) = X`` under every
    alignment: Kimura two-parameter distance with transitions at CpG sites
    down-weighted, the metric used for published repeat landscapes. Keys match
    ``_parse_repeatmasker_out`` rows: ``(query, start, end, repeat)``.
    """
    kimura: dict[tuple[str, int, int, str], float] = {}
    key: tuple[str, int, int, str] | None = None
    with open(align_path) as handle:
        for line in handle:
            if line[:1].isdigit():
                header = _align_header_key(line.split())
                if header is not None:
                    key = header
                    continue
            if key is not None and line.startswith(_ALIGN_KIMURA_PREFIX):
                try:
                    value = float(line.split("=", 1)[1].split()[0])
                except (IndexError, ValueError):
                    key = None
                    continue
                kimura.setdefault(key, value)
                key = None
    return kimura


def jukes_cantor_percent(p_distance_percent: float, *, cap: float = 50.0) -> float:
    """Jukes-Cantor-correct a raw mismatch percentage (fallback only)."""
    p = max(p_distance_percent, 0.0) / 100.0
    inner = 1.0 - 4.0 * p / 3.0
    if inner <= 0:
        return cap
    return min(-0.75 * math.log(inner) * 100.0, cap)


def apply_kimura_to_rows(
    rows: list[dict], kimura: dict[tuple[str, int, int, str], float]
) -> int:
    """Set ``row["landscape_div"]`` from the .align Kimura values.

    Must run while ``row["family"]`` still holds the name RepeatMasker saw.
    Rows without an alignment record get the Jukes-Cantor-corrected ``.out``
    value, so every hit stays in the landscape on a comparable scale.
    Returns the number of unmatched rows.
    """
    unmatched = 0
    for row in rows:
        key = (row["scaffold"], row["start"], row["end"], row["family"])
        value = kimura.get(key)
        if value is None:
            unmatched += 1
            value = jukes_cantor_percent(row["divergence"])
        row["landscape_div"] = value
    return unmatched


def _bin_start(div: float, bin_width: float, div_max: float) -> float:
    div = min(max(div, 0.0), div_max)
    b = math.floor(div / bin_width + 1e-9) * bin_width
    # Values at the ceiling land in the last bin, not a bin beyond the range.
    b = min(b, div_max - bin_width)
    return round(b, 6)


def _least_diverged_segments(hits: list[tuple[int, int, float]]):
    """Yield ``(divergence, bp)`` counting each base once.

    ``hits`` are ``(start, end, divergence)`` on one sequence (1-based,
    inclusive). Where hits overlap, the base goes to the least-diverged hit --
    the copy's best match -- so the total equals the merged covered bp.
    """
    if not hits:
        return
    points = sorted({h[0] for h in hits} | {h[1] + 1 for h in hits})
    order = sorted(hits)
    heap: list[tuple[float, int]] = []
    k = 0
    for left, right in zip(points, points[1:]):
        while k < len(order) and order[k][0] <= left:
            start, end, div = order[k]
            heapq.heappush(heap, (div, end))
            k += 1
        while heap and heap[0][1] < left:
            heapq.heappop(heap)
        if heap:
            yield heap[0][0], right - left


def _landscape_from_rows(
    rows: list[dict],
    *,
    bin_width: float = 1.0,
    div_max: float = 50.0,
    div_key: str = "divergence",
) -> list[tuple[str, float, int]]:
    """Bin covered bp by divergence, per family (repeat landscape).

    Overlapping hits of the same family are not double-counted: each base is
    assigned to the least-diverged hit covering it, so a family's landscape
    total equals its ``covered_bp`` in the summary table.
    """
    groups: dict[tuple[str, str], list[tuple[int, int, float]]] = {}
    for r in rows:
        div = r.get(div_key, r["divergence"])
        groups.setdefault((r["family"], r["scaffold"]), []).append(
            (r["start"], r["end"], div)
        )
    bins: dict[tuple[str, float], int] = {}
    for (family, _scaffold), hits in groups.items():
        for div, bp in _least_diverged_segments(hits):
            key = (family, _bin_start(div, bin_width, div_max))
            bins[key] = bins.get(key, 0) + bp
    return sorted(
        ((fam, b, bp) for (fam, b), bp in bins.items()),
        key=lambda t: (t[0], t[1]),
    )


_LANDSCAPE_METRIC_LABELS = {
    "kimura": "kimura_cpg_adjusted",
    "raw": "raw_perc_div",
}


def _write_landscape_tsv(
    landscape: list[tuple[str, float, int]],
    path: Path,
    *,
    bin_width: float = 1.0,
    metric: str = "raw",
) -> Path:
    label = _LANDSCAPE_METRIC_LABELS.get(metric, metric)
    lines = ["family\tdivergence_bin\tbp\tdivergence_bin_end\tmetric"]
    for fam, b, bp in landscape:
        lines.append(f"{fam}\t{b}\t{bp}\t{round(b + bin_width, 6)}\t{label}")
    path.write_text("\n".join(lines) + "\n")
    return path


# RepeatMasker's FastaDB rejects sequence ids longer than this.
_RM_MAX_ID_LEN = 50


def _sanitize_fasta_ids(src: Path, dest: Path, prefix: str) -> dict[str, str]:
    """Rewrite ``src`` to ``dest`` with short unique ids (RepeatMasker <=50).

    Returns a ``{short_id: original_id}`` map so output names can be restored.
    """
    from Bio.SeqRecord import SeqRecord

    mapping: dict[str, str] = {}
    records = []
    for index, rec in enumerate(SeqIO.parse(str(src), "fasta")):
        short = f"{prefix}{index}"
        mapping[short] = rec.id
        records.append(SeqRecord(rec.seq, id=short, description=""))
    SeqIO.write(records, str(dest), "fasta")
    return mapping


def quantify_library_repeatmasker(
    library_fasta: Path,
    target_fasta: Path,
    out_dir: Path,
    *,
    threads: int = 4,
    div: int = 40,
    sensitivity: str = "sensitive",
    rm_out: Path | None = None,
    cancel=None,
) -> QuantifyResult:
    """Quantify ``library_fasta`` across ``target_fasta`` with RepeatMasker ``-lib``.

    Runs RepeatMasker (target masked with the library) unless an existing
    ``.out`` is supplied, then aggregates per family and builds a per-family
    divergence landscape. Returns the same QuantifyResult shape as the blastn
    engine so the UI can treat both uniformly.
    """
    from repeatviewer.backends.te_library_benchmark import run_repeatmasker_lib_benchmark

    library_fasta = Path(library_fasta).expanduser().resolve()
    target_fasta = Path(target_fasta).expanduser().resolve()
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    for path in (library_fasta, target_fasta):
        if not path.is_file():
            raise FileNotFoundError(f"FASTA not found: {path}")

    library_map: dict[str, str] = {}
    if rm_out is None:
        # RepeatMasker's FastaDB caps ids at 50 chars, so feed it short ids and
        # map the family (library) names back afterward.
        short_lib = out_dir / f"{target_fasta.stem}_rm_lib.fa"
        short_query = out_dir / f"{target_fasta.stem}_rm.fa"
        library_map = _sanitize_fasta_ids(library_fasta, short_lib, "L")
        _sanitize_fasta_ids(target_fasta, short_query, "q")
        rm_out = run_repeatmasker_lib_benchmark(
            query_fasta=short_query,
            standard_fasta=short_lib,
            output_dir=out_dir,
            threads=threads,
            div=div,
            sensitivity=sensitivity,
            write_alignments=True,
        )
    rm_out = Path(rm_out).expanduser().resolve()

    rows = _parse_repeatmasker_out(rm_out)
    # Divergence (table + landscape): CpG-adjusted Kimura from the .align file when
    # present (matched before family names are restored, since the .align
    # carries the names RepeatMasker saw); otherwise the raw .out column.
    metric = "raw"
    unmatched = 0
    align_path = align_path_for_out(rm_out)
    if align_path is not None:
        kimura = parse_repeatmasker_align_kimura(align_path)
        if kimura:
            unmatched = apply_kimura_to_rows(rows, kimura)
            metric = "kimura"
            if unmatched:
                logger.warning(
                    "Repeat landscape: %d of %d hits had no .align record; "
                    "used Jukes-Cantor-corrected .out divergence for those.",
                    unmatched, len(rows),
                )
    if library_map:  # restore original family names
        for row in rows:
            row["family"] = library_map.get(row["family"], row["family"])
    target_length = _fasta_total_length(target_fasta)
    result = QuantifyResult(
        engine="repeatmasker",
        target_length=target_length,
        families=_aggregate_rm_rows(
            rows, target_length, monomer_lengths=_fasta_lengths(library_fasta)
        ),
        rm_out=rm_out,
        landscape_metric=metric,
        landscape_unmatched_hits=unmatched,
    )
    result.summary_tsv = _write_summary_tsv(
        result, out_dir / f"{target_fasta.stem}_repeat_quant_rm.tsv"
    )
    bin_width = 1.0
    result.landscape = _landscape_from_rows(
        rows,
        bin_width=bin_width,
        div_key="landscape_div" if metric == "kimura" else "divergence",
    )
    result.landscape_bin_width = bin_width
    result.landscape_tsv = _write_landscape_tsv(
        result.landscape,
        out_dir / f"{target_fasta.stem}_repeat_landscape.tsv",
        bin_width=bin_width,
        metric=metric,
    )
    return result


# ---------------------------------------------------------------------------
# HTML result rendering (for the preview pane)
# ---------------------------------------------------------------------------

def _html_escape(text: str) -> str:
    return (
        str(text)
        .replace("&", "&amp;")
        .replace("<", "&lt;")
        .replace(">", "&gt;")
    )


# Column display headers and their short definitions (shown as a side panel).
_COLUMN_LABELS = {
    "family": "family",
    "monomer_len": "monomer_len",
    "copy_number": "copy_number",
    "hits": "hits",
    "arrays": "arrays",
    "covered_bp": "covered_bp",
    "pct_of_target": "%target (%)",
    "mean_array_bp": "mean_array_bp",
    "max_array_bp": "max_array_bp",
    "mean_identity": "mean_identity",
    "mean_divergence": "mean_divergence",
    "scaffolds": "scaffolds",
}
_COLUMN_DEFS = [
    ("family", "Library family that was matched."),
    ("monomer_len", "Length of the library monomer (period size), bp."),
    ("copy_number", "Estimated copies = covered_bp ÷ monomer length. Shown as a "
                    "fraction when a family covers less than one full monomer."),
    ("hits", "Number of alignment fragments (raw hit count)."),
    ("arrays", "Hits merged into contiguous blocks (≈ distinct array loci)."),
    ("covered_bp", "Target bp covered by the family (overlaps counted once)."),
    ("%target (%)", "covered_bp as a percentage of the target length."),
    ("mean_array_bp", "Average tandem-array length = covered_bp ÷ arrays."),
    ("max_array_bp", "Largest single contiguous array, bp."),
    ("mean_identity", "Length-weighted mean % identity to the consensus "
                      "(RepeatMasker: 100 − raw substitution %; "
                      "high = young / conserved copies)."),
    ("mean_divergence", "Length-weighted mean divergence from the consensus. "
                        "Kimura: CpG-adjusted Kimura 2-parameter % from the "
                        "RepeatMasker .align file. raw: .out substitution % "
                        "(runs without a .align). High = old / diverged copies."),
    ("scaffolds", "Number of distinct target sequences hit."),
]


def _global_family_order(results: list[QuantifyResult]) -> list[str]:
    """One family order used by every table: by max covered_bp across engines."""
    max_cov: dict[str, int] = {}
    for result in results:
        for fam in result.families:
            max_cov[fam.family] = max(max_cov.get(fam.family, 0), fam.covered_bp)
    return sorted(max_cov, key=lambda name: max_cov[name], reverse=True)


def write_summary_html(
    results: list[QuantifyResult],
    path: Path,
    *,
    target_label: str = "",
    warnings: list[str] | None = None,
) -> Path:
    """Render one or more QuantifyResults as a styled HTML summary page."""
    path = Path(path)
    engine_display = {"blastn": "BLASTn", "repeatmasker": "RepeatMasker"}
    family_order = _global_family_order(results)
    blocks: list[str] = []
    for result in results:
        pct = (
            result.total_covered_bp / result.target_length * 100.0
            if result.target_length else 0.0
        )
        header = (
            "<div class='qhead'>"
            f"<span class='qname'>Quantifier: "
            f"{engine_display.get(result.engine, result.engine)}</span>"
            f"<span class='qcov'>Covered: {result.total_covered_bp:,} bp "
            f"({pct:.2f}%)</span></div>"
        )
        show_div = result.engine == "repeatmasker"
        cols = ["family", "monomer_len", "copy_number", "hits", "arrays",
                "covered_bp", "pct_of_target", "mean_array_bp", "max_array_bp",
                "mean_identity"] + (["mean_divergence"] if show_div else []) + ["scaffolds"]
        div_label = {"kimura": "mean_divergence (Kimura)",
                     "raw": "mean_divergence (raw)"}.get(
            result.divergence_metric, _COLUMN_LABELS["mean_divergence"])
        thead = "".join(
            f"<th>{div_label if c == 'mean_divergence' else _COLUMN_LABELS[c]}</th>"
            for c in cols
        )
        by_name = {f.family: f for f in result.families}
        body_rows = []
        for fam_name in family_order:  # same order in every table
            f = by_name.get(fam_name)
            if f is None:
                continue
            cells = [
                _html_escape(f.family), f"{f.monomer_len:,}",
                _format_copy_number(f.copy_number, thousands=True),
                f"{f.hits}", f"{f.arrays}", f"{f.covered_bp:,}",
                f"{f.pct_of_target:.3f}", f"{f.mean_array_bp:,.0f}",
                f"{f.max_array_bp:,}", f"{f.mean_identity:.1f}",
            ]
            if show_div:
                cells.append(f"{f.mean_divergence:.1f}")
            cells.append(f"{f.scaffolds}")
            body_rows.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
        blocks.append(
            f"{header}<table><thead><tr>{thead}</tr></thead>"
            f"<tbody>{''.join(body_rows)}</tbody></table>"
        )
    definitions = "".join(
        f"<dt>{_html_escape(k)}</dt><dd>{_html_escape(v)}</dd>"
        for k, v in _COLUMN_DEFS
    )
    aside = f"<aside><h4>Column definitions</h4><dl>{definitions}</dl></aside>"
    style = (
        "body{font-family:Arial,sans-serif;margin:16px;color:#111;}"
        ".wrap{display:flex;align-items:flex-start;gap:28px;}"
        ".tables{flex:1;min-width:0;}"
        ".genome{font-size:14px;margin:0 0 16px;color:#333;}"
        ".qhead{margin:16px 0 4px;}"
        ".qname{font-weight:bold;font-size:14px;}"
        ".qcov{color:#555;margin-left:14px;}"
        "table{border-collapse:collapse;margin:4px 0 20px;font-size:13px;}"
        "th,td{border:1px solid #ccc;padding:4px 8px;text-align:right;}"
        "th{background:#f2f2f2;} td:first-child,th:first-child{text-align:left;}"
        "tr:nth-child(even){background:#fafafa;}h2{margin:0 0 6px;}"
        "aside{width:300px;flex:none;font-size:12px;color:#333;"
        "border-left:2px solid #eee;padding-left:16px;}"
        "aside h4{margin:0 0 8px;} dt{font-weight:bold;margin-top:8px;}"
        "dd{margin:2px 0 0;color:#555;}"
        ".warn{background:#fff8e1;border:1px solid #f0c36d;border-radius:6px;"
        "padding:8px 12px;margin:0 0 16px;font-size:13px;color:#5c4400;}"
        ".warn ul{margin:6px 0 0;padding-left:20px;} .warn b{color:#7a5b00;}"
    )
    glen = results[0].target_length if results else 0
    genome_line = (
        f"Genome: {_html_escape(target_label)}, {glen:,} bp"
        if target_label else f"Genome: {glen:,} bp"
    )
    banner = ""
    if warnings:
        items = "".join(f"<li>{_html_escape(w)}</li>" for w in warnings)
        banner = f"<div class='warn'><b>⚠ Library warning</b><ul>{items}</ul></div>"
    path.write_text(
        f"<!doctype html><html><head><meta charset='utf-8'><style>{style}</style>"
        f"</head><body><h2>Repeat quantification</h2>"
        f"<p class='genome'>{genome_line}</p>{banner}"
        f"<div class='wrap'><div class='tables'>{''.join(blocks)}</div>{aside}</div>"
        f"</body></html>"
    )
    return path


def write_landscape_html(
    result: QuantifyResult, path: Path, *, target_label: str = ""
) -> Path:
    """Stacked-bar divergence landscape (bp vs divergence, per family).

    Bars are centred on their bins (a 1 % bar for 8-9 % sits at 8.5).
    """
    import plotly.graph_objects as go
    import plotly.io as pio

    path = Path(path)
    families = sorted({fam for fam, _b, _bp in result.landscape})
    bins = sorted({b for _f, b, _bp in result.landscape})
    lookup = {(fam, b): bp for fam, b, bp in result.landscape}
    width = result.landscape_bin_width or 1.0
    centres = [round(b + width / 2, 6) for b in bins]
    ranges = [f"{b:g}–{b + width:g}%" for b in bins]
    fig = go.Figure()
    for fam in families:
        fig.add_trace(
            go.Bar(
                x=centres,
                y=[lookup.get((fam, b), 0) for b in bins],
                width=width,
                name=fam,
                customdata=ranges,
                hovertemplate=(
                    f"{_html_escape(fam)}<br>%{{customdata}}: %{{y:,}} bp"
                    "<extra></extra>"
                ),
            )
        )
    title = (
        f"Repeat divergence landscape — {target_label}"
        if target_label else "Repeat divergence landscape"
    )
    fig.update_layout(
        barmode="stack",
        bargap=0.1,
        title=title,
        xaxis_title=(
            "Kimura divergence, CpG-adjusted (%)"
            if result.landscape_metric == "kimura"
            else "Divergence from consensus, raw (%)"
        ),
        yaxis_title="Covered bp",
        plot_bgcolor="white",
        paper_bgcolor="white",
        font={"family": "Arial, sans-serif", "size": 12, "color": "black"},
    )
    # Draw solid x- and y-axis lines.
    axis_line = {
        "showline": True,
        "linecolor": "black",
        "linewidth": 1,
        "ticks": "outside",
    }
    fig.update_xaxes(**axis_line)
    fig.update_yaxes(**axis_line)
    pio.write_html(fig, file=str(path), auto_open=False, config={"responsive": True})
    return path


def repeatmasker_executable_available() -> bool:
    """True if RepeatMasker is on PATH (needed for the RepeatMasker engine)."""
    from repeatviewer.conda_tools import which_tool

    return which_tool("RepeatMasker") is not None


# --- Output folders -----------------------------------------------------------
# Every run writes into its own new folder, so a repeat run never overwrites or
# mixes with an earlier one: <base>/<genome stem>_repeat_quant_<YYYY-MM-DD_HHMM>,
# where <base> defaults to the genome's own folder.

RUN_STAMP_FORMAT = "%Y-%m-%d_%H%M"
RUN_PARAMETERS_FILENAME = "run_parameters.json"
_RM_PARAM_KEYS = ("rm_div", "rm_sensitivity")


def run_stamp(now: datetime | None = None) -> str:
    """Timestamp used in run folder names, e.g. ``2026-09-17_1432``."""
    return (now or datetime.now()).strftime(RUN_STAMP_FORMAT)


def repeat_quant_output_dir(
    target_fasta: Path, stamp: str, base_dir: Path | None = None
) -> Path:
    """Planned run folder for one genome (not created, no clash handling)."""
    target_fasta = Path(target_fasta).expanduser()
    parent = Path(base_dir).expanduser() if base_dir else target_fasta.parent
    return parent / f"{target_fasta.stem}_repeat_quant_{stamp}"


def create_repeat_quant_output_dir(
    target_fasta: Path, stamp: str, base_dir: Path | None = None
) -> Path:
    """Create and return a fresh run folder; adds ``_2``, ``_3``... on a clash.

    ``mkdir(exist_ok=False)`` claims the name atomically, so two genomes with
    the same file name (or two runs in the same minute) never share a folder.
    """
    planned = repeat_quant_output_dir(target_fasta, stamp, base_dir).resolve()
    planned.parent.mkdir(parents=True, exist_ok=True)
    candidate, n = planned, 1
    while True:
        try:
            candidate.mkdir(exist_ok=False)
            return candidate
        except FileExistsError:
            n += 1
            candidate = planned.with_name(f"{planned.name}_{n}")


def engine_params(engine: str, params: dict) -> dict:
    """Only the parameters the chosen engine actually uses."""
    use_blast = engine in ("blastn", "both")
    use_rm = engine in ("repeatmasker", "both")
    out = {}
    for key, value in params.items():
        if key.startswith("blast_") and not use_blast:
            continue
        if key in _RM_PARAM_KEYS and not use_rm:
            continue
        out[key] = value
    return out


def write_run_parameters(
    out_dir: Path,
    *,
    library_fasta: Path,
    target_fasta: Path,
    engine: str,
    threads: int,
    params: dict,
    started: datetime,
    status: str = "running",
    finished: datetime | None = None,
    outputs: list[Path] | None = None,
    error: str = "",
) -> Path:
    """Write (or rewrite) ``run_parameters.json`` describing one run.

    Written once with status "running" before the tools start and again at the
    end ("completed", "failed" or "cancelled"), so an interrupted run still
    leaves a record of how it was started.
    """
    from repeatviewer import __version__

    out_dir = Path(out_dir)
    record = {
        "tool": "RepeatViewer - Repeat quantify",
        "repeatviewer_version": __version__,
        "status": status,
        "started": started.isoformat(timespec="seconds"),
        "finished": finished.isoformat(timespec="seconds") if finished else None,
        "library": str(Path(library_fasta).expanduser().resolve()),
        "genome": str(Path(target_fasta).expanduser().resolve()),
        "output_dir": str(out_dir.resolve()),
        "engine": engine,
        "threads": int(threads),
        "parameters": engine_params(engine, params),
        "outputs": sorted(
            str(Path(p).resolve().relative_to(out_dir.resolve()))
            if Path(p).resolve().is_relative_to(out_dir.resolve())
            else str(Path(p).resolve())
            for p in (outputs or [])
        ),
        "platform": platform.platform(),
    }
    if error:
        record["error"] = error
    path = out_dir / RUN_PARAMETERS_FILENAME
    tmp = path.with_suffix(".json.tmp")
    tmp.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    os.replace(tmp, path)
    return path


def quantify_library(
    library_fasta: Path,
    target_fasta: Path,
    out_dir: Path,
    *,
    engine: str = "both",
    threads: int = 4,
    blast_word_size: int = 11,
    blast_evalue: float = 1e-5,
    blast_min_identity: float = 0.0,
    blast_min_coverage: float = 0.0,
    rm_div: int = 40,
    rm_sensitivity: str = "sensitive",
    cancel=None,
) -> list[QuantifyResult]:
    """Run the requested engine(s): ``blastn``, ``repeatmasker``, or ``both``."""
    results: list[QuantifyResult] = []
    if engine in ("blastn", "both"):
        results.append(
            quantify_library_blastn(
                library_fasta, target_fasta, out_dir,
                e_value=blast_evalue,
                min_identity=blast_min_identity,
                min_coverage=blast_min_coverage,
                word_size=blast_word_size,
                threads=threads,
                cancel=cancel,
            )
        )
    if engine in ("repeatmasker", "both"):
        results.append(
            quantify_library_repeatmasker(
                library_fasta, target_fasta, out_dir,
                threads=threads,
                div=rm_div,
                sensitivity=rm_sensitivity,
                cancel=cancel,
            )
        )
    return results
