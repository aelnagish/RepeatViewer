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

import os
import platform
import shutil
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO

from teviewer.backends.gui_functions import (
    blast_database_basename,
    check_database,
    get_original_file_path,
)
from teviewer.backends.job_cancel import run_subprocess_with_cancel

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
    mean_divergence: float = 0.0   # RepeatMasker only (0 for blastn)
    scaffolds: int = 0         # distinct target sequences hit
    monomer_len: int = 0       # library monomer length (period size)
    copy_number: int = 0       # covered_bp / monomer_len
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
    rm_out: Path | None = None

    @property
    def total_covered_bp(self) -> int:
        return sum(f.covered_bp for f in self.families)


def _fasta_total_length(fasta: Path) -> int:
    return sum(len(rec.seq) for rec in SeqIO.parse(str(fasta), "fasta"))


def _fasta_lengths(fasta: Path) -> dict[str, int]:
    """Map each record id to its length (used for satellite copy number)."""
    return {rec.id: len(rec.seq) for rec in SeqIO.parse(str(fasta), "fasta")}


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
                copy_number=round(covered / mono) if mono else 0,
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
        "mean_divergence\tscaffolds"
    ]
    for f in result.families:
        lines.append(
            f"{f.family}\t{f.monomer_len}\t{f.copy_number}\t{f.hits}\t"
            f"{f.arrays}\t{f.covered_bp}\t{f.pct_of_target}\t"
            f"{round(f.mean_array_bp)}\t{f.max_array_bp}\t{f.mean_identity}\t"
            f"{f.mean_divergence}\t{f.scaffolds}"
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


def quantify_library_blastn(
    library_fasta: Path,
    target_fasta: Path,
    out_dir: Path,
    *,
    e_value: float = 1e-5,
    min_identity: float = 0.0,
    word_size: int = 11,
    threads: int = 1,
    use_system_blast: bool = False,
    cancel=None,
) -> QuantifyResult:
    """Quantify ``library_fasta`` occurrences across ``target_fasta`` via blastn.

    Builds a nucleotide DB from the target, searches the library against it,
    then aggregates per family. Returns a QuantifyResult and writes a summary
    TSV and an arrays GFF into ``out_dir``.
    """
    library_fasta = Path(library_fasta).expanduser().resolve()
    target_fasta = Path(target_fasta).expanduser().resolve()
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    for path in (library_fasta, target_fasta):
        if not path.is_file():
            raise FileNotFoundError(f"FASTA not found: {path}")

    os_type = platform.system()
    # Build the BLAST database in a throwaway scratch dir so the ~8 index files
    # don't clutter the results folder; removed after the search.
    scratch = Path(tempfile.mkdtemp(prefix="repeat_quant_db_", dir=str(out_dir)))
    tab_path = out_dir / f"{target_fasta.stem}_vs_{library_fasta.stem}.blastn.tsv"
    try:
        db = check_database(
            str(target_fasta),
            str(scratch),
            os_type=os_type,
            use_system_blast=use_system_blast,
        )
        if db in {"makeblastdb_not_found", "makeblastdb_got_error"}:
            raise RuntimeError(
                "Could not build a BLAST database for the target sequences "
                f"({db}). Check that makeblastdb is available."
            )
        db_path = os.path.join(str(scratch), blast_database_basename(str(target_fasta)))

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
        # Cancellable runner so Stop actually terminates BLASTn on big inputs.
        run_subprocess_with_cancel(cmd, cancel=cancel, check=True)
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

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
        fam["div_len"].append((r["divergence"], length))
        fam["by_scaffold"].setdefault(r["scaffold"], []).append((r["start"], r["end"]))

    families: list[FamilyQuant] = []
    for fam_name, data in acc.items():
        arrays, covered, max_array = _array_stats(data["by_scaffold"])
        total_len = sum(length for _d, length in data["div_len"]) or 1
        mean_div = sum(d * length for d, length in data["div_len"]) / total_len
        pct = (covered / target_length * 100.0) if target_length else 0.0
        mono = int(monomer_lengths.get(fam_name, 0))
        families.append(
            FamilyQuant(
                family=fam_name,
                hits=len(data["div_len"]),
                arrays=arrays,
                covered_bp=covered,
                pct_of_target=round(pct, 4),
                mean_identity=round(100.0 - mean_div, 2),
                mean_divergence=round(mean_div, 2),
                scaffolds=len(data["by_scaffold"]),
                monomer_len=mono,
                copy_number=round(covered / mono) if mono else 0,
                mean_array_bp=round(covered / arrays, 1) if arrays else 0.0,
                max_array_bp=max_array,
            )
        )
    families.sort(key=lambda f: f.covered_bp, reverse=True)
    return families


def _landscape_from_rows(
    rows: list[dict], *, bin_width: float = 1.0, div_max: float = 50.0
) -> list[tuple[str, float, int]]:
    """Bin covered bp by divergence, per family (for a Kimura landscape)."""
    bins: dict[tuple[str, float], int] = {}
    for r in rows:
        div = min(max(r["divergence"], 0.0), div_max)
        bin_start = float(int(div / bin_width) * bin_width)
        length = r["end"] - r["start"] + 1
        key = (r["family"], bin_start)
        bins[key] = bins.get(key, 0) + length
    return sorted(
        ((fam, b, bp) for (fam, b), bp in bins.items()),
        key=lambda t: (t[0], t[1]),
    )


def _write_landscape_tsv(
    landscape: list[tuple[str, float, int]], path: Path
) -> Path:
    lines = ["family\tdivergence_bin\tbp"]
    for fam, b, bp in landscape:
        lines.append(f"{fam}\t{b}\t{bp}")
    path.write_text("\n".join(lines) + "\n")
    return path


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
    from teviewer.backends.te_library_benchmark import run_repeatmasker_lib_benchmark

    library_fasta = Path(library_fasta).expanduser().resolve()
    target_fasta = Path(target_fasta).expanduser().resolve()
    out_dir = Path(out_dir).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    for path in (library_fasta, target_fasta):
        if not path.is_file():
            raise FileNotFoundError(f"FASTA not found: {path}")

    if rm_out is None:
        rm_out = run_repeatmasker_lib_benchmark(
            query_fasta=target_fasta,
            standard_fasta=library_fasta,
            output_dir=out_dir,
            threads=threads,
            div=div,
            sensitivity=sensitivity,
        )
    rm_out = Path(rm_out).expanduser().resolve()

    rows = _parse_repeatmasker_out(rm_out)
    target_length = _fasta_total_length(target_fasta)
    result = QuantifyResult(
        engine="repeatmasker",
        target_length=target_length,
        families=_aggregate_rm_rows(
            rows, target_length, monomer_lengths=_fasta_lengths(library_fasta)
        ),
        rm_out=rm_out,
    )
    result.summary_tsv = _write_summary_tsv(
        result, out_dir / f"{target_fasta.stem}_repeat_quant_rm.tsv"
    )
    result.landscape = _landscape_from_rows(rows)
    result.landscape_tsv = _write_landscape_tsv(
        result.landscape, out_dir / f"{target_fasta.stem}_repeat_landscape.tsv"
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
    ("copy_number", "Estimated copies = covered_bp ÷ monomer length."),
    ("hits", "Number of alignment fragments (raw hit count)."),
    ("arrays", "Hits merged into contiguous blocks (≈ distinct array loci)."),
    ("covered_bp", "Target bp covered by the family (overlaps counted once)."),
    ("%target (%)", "covered_bp as a percentage of the target length."),
    ("mean_array_bp", "Average tandem-array length = covered_bp ÷ arrays."),
    ("max_array_bp", "Largest single contiguous array, bp."),
    ("mean_identity", "Length-weighted mean % identity to the consensus "
                      "(high = young / conserved copies)."),
    ("mean_divergence", "RepeatMasker Kimura substitution % "
                        "(high = old / diverged copies)."),
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
        thead = "".join(f"<th>{_COLUMN_LABELS[c]}</th>" for c in cols)
        by_name = {f.family: f for f in result.families}
        body_rows = []
        for fam_name in family_order:  # same order in every table
            f = by_name.get(fam_name)
            if f is None:
                continue
            cells = [
                _html_escape(f.family), f"{f.monomer_len:,}", f"{f.copy_number:,}",
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
    """Stacked-bar Kimura divergence landscape (bp vs divergence, per family)."""
    import plotly.graph_objects as go
    import plotly.io as pio

    path = Path(path)
    families = sorted({fam for fam, _b, _bp in result.landscape})
    bins = sorted({b for _f, b, _bp in result.landscape})
    lookup = {(fam, b): bp for fam, b, bp in result.landscape}
    fig = go.Figure()
    for fam in families:
        fig.add_trace(
            go.Bar(
                x=bins,
                y=[lookup.get((fam, b), 0) for b in bins],
                name=fam,
            )
        )
    title = (
        f"Repeat divergence landscape — {target_label}"
        if target_label else "Repeat divergence landscape"
    )
    fig.update_layout(
        barmode="stack",
        title=title,
        xaxis_title="Divergence (%)",
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
    import shutil

    return shutil.which("RepeatMasker") is not None


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
