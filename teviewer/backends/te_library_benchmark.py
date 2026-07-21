"""TE consensus-library benchmark (Perfect / Good / Present / Not found).

Python port of the RepeatModeler2-style ``get_family_summary_paper.sh`` workflow.
"""

from __future__ import annotations

import re
import shutil
import subprocess
from teviewer.external_process import run_logged_process
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

from teviewer.backends.blast_genome_plot import read_repeatmasker_out

_LTR_RE = re.compile(r"-LTR|_LTR", re.IGNORECASE)


@dataclass(frozen=True)
class LibraryBenchmarkCounts:
    perfect: int
    good: int
    all_present: int
    only_present: int
    not_found: int
    total_te: int
    total_ltr: int
    perfect_no_ltr: int
    good_no_ltr: int
    all_present_no_ltr: int
    all_present_only_ltr: int
    perfect_only_ltr: int
    good_only_ltr: int


def _fasta_sequence_count(path: Path) -> int:
    count = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">"):
                count += 1
    return count


def _fasta_ltr_count(path: Path) -> int:
    count = 0
    with path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith(">") and _LTR_RE.search(line):
                count += 1
    return count


def _normalize_strand(raw: object) -> str:
    s = str(raw).strip()
    if s in {"C", "c", "-"}:
        return "C"
    return "+"


def _alignment_metrics(row: pd.Series) -> tuple[float, float, float, float]:
    """Return family_length, target_cov, query_cov, divergence."""
    rs = float(row["repeat_start"])
    re = float(row["repeat_end"])
    rl = float(row["repeat_left"])
    qs = float(row["query_start"])
    qe = float(row["query_end"])
    ql = float(row["query_left"])
    div = float(row["perc_div"])
    strand = _normalize_strand(row["strand"])

    if strand == "+":
        family_length = rl + re
        align_family = re - rs
        query_length = qe + ql
        align_query = qe - qs
    else:
        family_length = rs + re
        align_family = re - rl
        query_length = qe + ql
        align_query = qe - qs

    family_length = max(family_length, 1.0)
    query_length = max(query_length, 1.0)
    target_cov = max(0.0, align_family) / family_length
    query_cov = max(0.0, align_query) / query_length
    return family_length, target_cov, query_cov, div


def _merged_interval_length(intervals: list[tuple[int, int]], *, max_gap: int = 10) -> int:
    if not intervals:
        return 0
    cleaned = [(min(a, b), max(a, b)) for a, b in intervals if max(a, b) > min(a, b)]
    if not cleaned:
        return 0
    cleaned.sort(key=lambda x: x[0])
    merged: list[tuple[int, int]] = [cleaned[0]]
    for start, end in cleaned[1:]:
        prev_start, prev_end = merged[-1]
        if start <= prev_end + max_gap:
            merged[-1] = (prev_start, max(prev_end, end))
        else:
            merged.append((start, end))
    return sum(end - start for start, end in merged)


def _family_interval(row: pd.Series) -> tuple[int, int]:
    rs = int(row["repeat_start"])
    re = int(row["repeat_end"])
    if _normalize_strand(row["strand"]) == "+":
        return rs, re
    rl = int(row["repeat_left"])
    return rs, re if re >= rs else rl


def run_library_benchmark_from_out(
    rm_out: Path,
    *,
    standard_fasta: Path | None = None,
) -> tuple[LibraryBenchmarkCounts, pd.DataFrame]:
    """Benchmark a lib-vs-lib RepeatMasker ``.out`` (query masked with ``-lib standard``)."""
    rm_out = rm_out.expanduser().resolve()
    hits = read_repeatmasker_out(rm_out)
    if hits.empty:
        raise ValueError(f"No RepeatMasker hits found in {rm_out}")

    rows: list[dict[str, object]] = []
    for _, row in hits.iterrows():
        family_length, target_cov, query_cov, div = _alignment_metrics(row)
        if query_cov <= 0.1:
            continue
        rows.append(
            {
                "family": str(row["repeat_name"]),
                "repeat_class": str(row["repeat_class"]),
                "family_length": family_length,
                "target_cov": target_cov,
                "query_cov": query_cov,
                "divergence": div,
                "strand": _normalize_strand(row["strand"]),
                "repeat_start": int(row["repeat_start"]),
                "repeat_end": int(row["repeat_end"]),
                "repeat_left": int(row["repeat_left"]),
            }
        )
    if not rows:
        raise ValueError("No alignments passed the minimum query-coverage filter (>10%).")

    df = pd.DataFrame(rows)
    families = sorted(df["family"].astype(str).unique())

    perfect: set[str] = set()
    for _, row in df.iterrows():
        if row["target_cov"] > 0.95 and row["query_cov"] > 0.95 and row["divergence"] < 5.0:
            perfect.add(str(row["family"]))

    present80: set[str] = set()
    good95: set[str] = set()
    for family in families:
        part = df[df["family"] == family]
        fam_len = float(part["family_length"].iloc[0])

        p_rows = part[part["divergence"] < 20.0]
        if not p_rows.empty:
            intervals = [_family_interval(r) for _, r in p_rows.iterrows()]
            cov = _merged_interval_length(intervals, max_gap=10)
            if fam_len > 0 and cov / fam_len > 0.8:
                present80.add(family)

        g_rows = part[part["divergence"] < 5.0]
        if not g_rows.empty:
            intervals = [_family_interval(r) for _, r in g_rows.iterrows()]
            cov = _merged_interval_length(intervals, max_gap=10)
            if fam_len > 0 and cov / fam_len > 0.95:
                good95.add(family)

    good = good95 - perfect
    present = present80
    only_present = present - perfect - good

    total_te = _fasta_sequence_count(standard_fasta) if standard_fasta else len(families)
    total_ltr = _fasta_ltr_count(standard_fasta) if standard_fasta else sum(
        1 for f in families if _LTR_RE.search(f)
    )

    def _count_no_ltr(names: set[str]) -> int:
        return sum(1 for n in names if not _LTR_RE.search(n))

    def _count_only_ltr(names: set[str]) -> int:
        return sum(1 for n in names if _LTR_RE.search(n))

    counts = LibraryBenchmarkCounts(
        perfect=len(perfect),
        good=len(good),
        all_present=len(present),
        only_present=len(only_present),
        not_found=max(0, total_te - len(present)),
        total_te=total_te,
        total_ltr=total_ltr,
        perfect_no_ltr=_count_no_ltr(perfect),
        good_no_ltr=_count_no_ltr(good),
        all_present_no_ltr=_count_no_ltr(present),
        all_present_only_ltr=_count_only_ltr(present),
        perfect_only_ltr=_count_only_ltr(perfect),
        good_only_ltr=_count_only_ltr(good),
    )
    return counts, df


def write_library_benchmark_summary(
    *,
    query_fasta: Path,
    standard_fasta: Path,
    counts: LibraryBenchmarkCounts,
    output_dir: Path,
) -> Path:
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    summary_path = output_dir / "results_summary.txt"
    lines = [
        query_fasta.name,
        standard_fasta.name,
        f"Perfect_number\t{counts.perfect}",
        f"Perfect_without_single_LTR\t{counts.perfect_no_ltr}",
        f"Perfect_only_single_LTR\t{counts.perfect_only_ltr}",
        f"Good_number\t{counts.good}",
        f"Good_without_single_LTR\t{counts.good_no_ltr}",
        f"Good_only_single_LTR\t{counts.good_only_ltr}",
        f"All_present\t{counts.all_present}",
        f"All_present_without_single_LTR\t{counts.all_present_no_ltr}",
        f"All_present_only_single_LTR\t{counts.all_present_only_ltr}",
        f"Only_present\t{counts.only_present}",
        f"Total_TE_number\t{counts.total_te}",
        f"Total_LTR_number\t{counts.total_ltr}",
        f"Not_found_number\t{counts.not_found}",
        "",
    ]
    summary_path.write_text("\n".join(lines), encoding="utf-8")
    return summary_path


def write_library_benchmark_plot_html(
    counts: LibraryBenchmarkCounts,
    output_path: Path,
    *,
    query_label: str = "",
    standard_label: str = "",
) -> Path:
    """Simple bar chart of Perfect / Good / Present / Not found."""
    import plotly.graph_objects as go

    categories = ["Perfect", "Good", "Only present", "Not found"]
    values = [
        counts.perfect,
        counts.good,
        counts.only_present,
        counts.not_found,
    ]
    colors = ["#f066b2", "#33a02c", "#ff7f00", "#919697"]

    fig = go.Figure(
        data=[
            go.Bar(
                x=categories,
                y=values,
                marker={"color": colors},
                text=values,
                textposition="outside",
            )
        ]
    )
    title = "TE library benchmark"
    if query_label and standard_label:
        title += f"<br><sup>{query_label} vs {standard_label}</sup>"
    fig.update_layout(
        title=title,
        yaxis={"title": "Families"},
        template="plotly_white",
        height=480,
        width=720,
    )
    output_path = output_path.expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.write_html(str(output_path), include_plotlyjs="cdn")
    return output_path


def repeatmasker_executable() -> str | None:
    return shutil.which("RepeatMasker")


def run_repeatmasker_lib_benchmark(
    *,
    query_fasta: Path,
    standard_fasta: Path,
    output_dir: Path,
    threads: int = 4,
) -> Path:
    """Run RepeatMasker ``-lib standard query`` and return the ``.out`` path."""
    exe = repeatmasker_executable()
    if not exe:
        raise RuntimeError(
            "RepeatMasker is not on PATH. Install RepeatMasker or provide an existing .out file."
        )
    query_fasta = query_fasta.expanduser().resolve()
    standard_fasta = standard_fasta.expanduser().resolve()
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    cmd = [
        exe,
        "-pa",
        str(max(1, int(threads))),
        "-s",
        "-div",
        "40",
        "-gff",
        "-no_is",
        "-norna",
        "-nolow",
        "-lib",
        str(standard_fasta),
        str(query_fasta),
        "-dir",
        str(output_dir),
    ]
    run_logged_process(
        cmd,
        tool="RepeatMasker",
        inputs=[standard_fasta, query_fasta],
        cwd=output_dir,
        check=True,
        output_path=output_dir,
    )
    out_path = output_dir / f"{query_fasta.name}.out"
    if not out_path.is_file():
        candidates = sorted(output_dir.glob("*.out"))
        if not candidates:
            raise RuntimeError("RepeatMasker finished but no .out file was found.")
        out_path = candidates[0]
    return out_path


def run_te_library_benchmark(
    *,
    query_fasta: Path,
    standard_fasta: Path,
    output_dir: Path | None = None,
    rm_out: Path | None = None,
    run_repeatmasker: bool = True,
    threads: int = 4,
) -> list[Path]:
    """Full library benchmark: optional RepeatMasker, then summary + plot."""
    query_fasta = query_fasta.expanduser().resolve()
    standard_fasta = standard_fasta.expanduser().resolve()
    work_dir = (output_dir or query_fasta.parent / f"{query_fasta.stem}_lib_benchmark").resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    out_file = rm_out.expanduser().resolve() if rm_out else None
    if out_file is None:
        if not run_repeatmasker:
            raise ValueError("Provide a RepeatMasker .out file or enable Run RepeatMasker.")
        out_file = run_repeatmasker_lib_benchmark(
            query_fasta=query_fasta,
            standard_fasta=standard_fasta,
            output_dir=work_dir,
            threads=threads,
        )

    counts, _detail = run_library_benchmark_from_out(
        out_file,
        standard_fasta=standard_fasta,
    )
    summary = write_library_benchmark_summary(
        query_fasta=query_fasta,
        standard_fasta=standard_fasta,
        counts=counts,
        output_dir=work_dir,
    )
    plot_html = write_library_benchmark_plot_html(
        counts,
        work_dir / f"{query_fasta.stem}_library_benchmark.html",
        query_label=query_fasta.name,
        standard_label=standard_fasta.name,
    )
    combined = work_dir / "final_combined_output.txt"
    combined.write_text(
        summary.read_text(encoding="utf-8"),
        encoding="utf-8",
    )
    return [summary, plot_html, combined, out_file]
