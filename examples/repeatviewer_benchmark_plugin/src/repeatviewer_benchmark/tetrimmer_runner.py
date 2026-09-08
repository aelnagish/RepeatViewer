"""Load and run TEtrimmer ``TE_anno_benchmarker.calculate_overlap``."""

from __future__ import annotations

import importlib.util
import os
import shutil
from pathlib import Path
from typing import Callable


def load_calculate_overlap(tetrimmer_home: Path) -> Callable:
    """Import ``calculate_overlap`` from a TEtrimmer clone (no pip install required)."""
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


def default_tetrimmer_home() -> Path | None:
    env = os.environ.get("TETRIMMER_HOME", "").strip()
    if env:
        path = Path(env).expanduser()
        if (path / "tetrimmer" / "TE_anno_benchmarker.py").is_file():
            return path
    return None


def require_bedtools() -> None:
    if not shutil.which("bedtools"):
        raise RuntimeError(
            "bedtools is not on PATH. Install BEDtools (required by TE_anno_benchmarker)."
        )


def benchmark_outputs(query_file: Path, output_dir: Path) -> list[Path]:
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
