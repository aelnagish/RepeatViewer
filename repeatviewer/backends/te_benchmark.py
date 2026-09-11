"""TE annotation confusion-matrix benchmarking (TEtrimmer TE_anno_benchmarker)."""

from __future__ import annotations

import importlib.util
import os
from pathlib import Path
from typing import Callable


def default_tetrimmer_home() -> Path | None:
    env = os.environ.get("TETRIMMER_HOME", "").strip()
    if env:
        path = Path(env).expanduser()
        if (path / "tetrimmer" / "TE_anno_benchmarker.py").is_file():
            return path
    candidates = [
        Path.home() / "Documents" / "TE_Trimmer" / "TEtrimmer",
        Path.home() / "Documents" / "GitHub" / "TEtrimmer",
    ]
    for path in candidates:
        if (path / "tetrimmer" / "TE_anno_benchmarker.py").is_file():
            return path
    return None


def load_calculate_overlap(tetrimmer_home: Path) -> Callable:
    """Import ``calculate_overlap`` from a TEtrimmer clone."""
    home = tetrimmer_home.expanduser().resolve()
    script = home / "tetrimmer" / "TE_anno_benchmarker.py"
    if not script.is_file():
        raise FileNotFoundError(
            f"TE_anno_benchmarker.py not found under {home}/tetrimmer/"
        )
    spec = importlib.util.spec_from_file_location("te_anno_benchmarker", script)
    if spec is None or spec.loader is None:
        raise ImportError(f"Cannot load module from {script}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    fn = getattr(module, "calculate_overlap", None)
    if fn is None:
        raise AttributeError("calculate_overlap not found in TE_anno_benchmarker.py")
    return fn


def require_bedtools() -> None:
    from repeatviewer.conda_tools import which_tool

    if not which_tool("bedtools"):
        raise RuntimeError(
            "bedtools is not on PATH. Install BEDtools (required by TE_anno_benchmarker)."
        )


def confusion_matrix_output_paths(query_file: Path, output_dir: Path) -> list[Path]:
    """Return paths TE_anno_benchmarker creates for one query."""
    base = query_file.name
    names = [
        f"{base}_sankey.html",
        f"{base}_heatmap_proportion_plot.pdf",
        f"{base}_TE_length.txt",
        f"{base}_confusion_matrix.txt",
        f"{base}_confusion_matrix_TP_TN_FP_FN.txt",
        f"{base}_spider_plot.pdf",
    ]
    return [output_dir / n for n in names if (output_dir / n).is_file()]


def run_confusion_matrix_benchmark(
    *,
    query_out: Path,
    reference_out: Path,
    genome_fasta: Path,
    output_dir: Path | None = None,
    tetrimmer_home: Path,
    less_classification: bool = False,
    debug: bool = False,
) -> list[Path]:
    """Run TEtrimmer ``calculate_overlap`` for one query RepeatMasker ``.out`` file."""
    require_bedtools()
    query_out = query_out.expanduser().resolve()
    reference_out = reference_out.expanduser().resolve()
    genome_fasta = genome_fasta.expanduser().resolve()
    out_dir = (output_dir or query_out.parent).expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    calculate_overlap = load_calculate_overlap(tetrimmer_home)
    calculate_overlap(
        query=str(query_out),
        reference=str(reference_out),
        output_dir=str(out_dir),
        genome_dir=str(genome_fasta),
        genome_length=None,
        less_classification=bool(less_classification),
        debug=bool(debug),
        query2=None,
    )
    created = confusion_matrix_output_paths(query_out, out_dir)
    if not created:
        raise RuntimeError("Benchmark finished but no expected output files were found.")
    return created
