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
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO

from teviewer.backends.gui_functions import (
    blast_database_basename,
    check_database,
    get_original_file_path,
)
from teviewer.external_process import run_logged_process

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


def _aggregate_blast_tab(tab_path: Path, target_length: int) -> list[FamilyQuant]:
    """Aggregate a blastn outfmt-6 file into per-family quantification."""
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
        arrays = 0
        covered = 0
        for _sseqid, ivals in data["by_scaffold"].items():
            merged = _merge_intervals(ivals)
            arrays += len(merged)
            covered += sum(end - start + 1 for start, end in merged)
        total_len = sum(length for _pi, length in data["pident_len"]) or 1
        mean_id = sum(pi * length for pi, length in data["pident_len"]) / total_len
        pct = (covered / target_length * 100.0) if target_length else 0.0
        families.append(
            FamilyQuant(
                family=fam_name,
                hits=hits,
                arrays=arrays,
                covered_bp=covered,
                pct_of_target=round(pct, 4),
                mean_identity=round(mean_id, 2),
                scaffolds=len(data["by_scaffold"]),
            )
        )
    families.sort(key=lambda f: f.covered_bp, reverse=True)
    return families


def _write_summary_tsv(result: QuantifyResult, path: Path) -> Path:
    lines = [
        "family\thits\tarrays\tcovered_bp\tpct_of_target\tmean_identity\tscaffolds"
    ]
    for f in result.families:
        lines.append(
            f"{f.family}\t{f.hits}\t{f.arrays}\t{f.covered_bp}\t"
            f"{f.pct_of_target}\t{f.mean_identity}\t{f.scaffolds}"
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
    word_size: int = 7,
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
    db = check_database(
        str(target_fasta),
        str(out_dir),
        os_type=os_type,
        use_system_blast=use_system_blast,
    )
    if db in {"makeblastdb_not_found", "makeblastdb_got_error"}:
        raise RuntimeError(
            "Could not build a BLAST database for the target sequences "
            f"({db}). Check that makeblastdb is available."
        )
    db_path = os.path.join(str(out_dir), blast_database_basename(str(target_fasta)))

    tab_path = out_dir / f"{target_fasta.stem}_vs_{library_fasta.stem}.blastn.tsv"
    cmd = [
        _blastn_executable(use_system_blast),
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
    run_logged_process(cmd, tool="Repeat quantify (blastn)", check=True, cwd=out_dir)

    target_length = _fasta_total_length(target_fasta)
    result = QuantifyResult(
        engine="blastn",
        target_length=target_length,
        families=_aggregate_blast_tab(tab_path, target_length),
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


def _aggregate_rm_rows(rows: list[dict], target_length: int) -> list[FamilyQuant]:
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
        arrays = 0
        covered = 0
        for _scaf, ivals in data["by_scaffold"].items():
            merged = _merge_intervals(ivals)
            arrays += len(merged)
            covered += sum(end - start + 1 for start, end in merged)
        total_len = sum(length for _d, length in data["div_len"]) or 1
        mean_div = sum(d * length for d, length in data["div_len"]) / total_len
        pct = (covered / target_length * 100.0) if target_length else 0.0
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
        )
    rm_out = Path(rm_out).expanduser().resolve()

    rows = _parse_repeatmasker_out(rm_out)
    target_length = _fasta_total_length(target_fasta)
    result = QuantifyResult(
        engine="repeatmasker",
        target_length=target_length,
        families=_aggregate_rm_rows(rows, target_length),
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


def write_summary_html(
    results: list[QuantifyResult],
    path: Path,
    *,
    target_label: str = "",
) -> Path:
    """Render one or more QuantifyResults as a styled HTML summary page."""
    path = Path(path)
    blocks: list[str] = []
    for result in results:
        header = (
            f"Engine: {result.engine} &nbsp;|&nbsp; "
            f"target {result.target_length:,} bp &nbsp;|&nbsp; "
            f"covered {result.total_covered_bp:,} bp "
            f"({(result.total_covered_bp / result.target_length * 100.0) if result.target_length else 0:.2f}%)"
        )
        show_div = result.engine == "repeatmasker"
        cols = ["family", "hits", "arrays", "covered_bp", "pct_of_target",
                "mean_identity"] + (["mean_divergence"] if show_div else []) + ["scaffolds"]
        thead = "".join(f"<th>{c}</th>" for c in cols)
        body_rows = []
        for f in result.families:
            cells = [
                _html_escape(f.family), f"{f.hits}", f"{f.arrays}",
                f"{f.covered_bp:,}", f"{f.pct_of_target:.3f}",
                f"{f.mean_identity:.1f}",
            ]
            if show_div:
                cells.append(f"{f.mean_divergence:.1f}")
            cells.append(f"{f.scaffolds}")
            body_rows.append("<tr>" + "".join(f"<td>{c}</td>" for c in cells) + "</tr>")
        blocks.append(
            f"<h3>{header}</h3><table><thead><tr>{thead}</tr></thead>"
            f"<tbody>{''.join(body_rows)}</tbody></table>"
        )
    style = (
        "body{font-family:Arial,sans-serif;margin:16px;color:#111;}"
        "table{border-collapse:collapse;margin:8px 0 24px;font-size:13px;}"
        "th,td{border:1px solid #ccc;padding:4px 8px;text-align:right;}"
        "th{background:#f2f2f2;} td:first-child,th:first-child{text-align:left;}"
        "tr:nth-child(even){background:#fafafa;}h2{margin:0 0 4px;}"
    )
    title = f"Repeat quantification — {_html_escape(target_label)}" if target_label else "Repeat quantification"
    path.write_text(
        f"<!doctype html><html><head><meta charset='utf-8'><style>{style}</style>"
        f"</head><body><h2>{title}</h2>{''.join(blocks)}</body></html>"
    )
    return path


def write_landscape_html(result: QuantifyResult, path: Path) -> Path:
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
    fig.update_layout(
        barmode="stack",
        title="Repeat divergence landscape",
        xaxis_title="Divergence (%)",
        yaxis_title="Covered bp",
        plot_bgcolor="white",
        paper_bgcolor="white",
        font={"family": "Arial, sans-serif", "size": 12, "color": "black"},
    )
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
    cancel=None,
) -> list[QuantifyResult]:
    """Run the requested engine(s): ``blastn``, ``repeatmasker``, or ``both``."""
    results: list[QuantifyResult] = []
    if engine in ("blastn", "both"):
        results.append(
            quantify_library_blastn(
                library_fasta, target_fasta, out_dir, threads=threads, cancel=cancel
            )
        )
    if engine in ("repeatmasker", "both"):
        results.append(
            quantify_library_repeatmasker(
                library_fasta, target_fasta, out_dir, threads=threads, cancel=cancel
            )
        )
    return results
