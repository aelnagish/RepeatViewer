from __future__ import annotations

import subprocess
from pathlib import Path

from repeatviewer_minimap2.executable import minimap2_executable
from repeatviewer_minimap2.options import Minimap2Options


def build_command(query: Path, reference: Path, opts: Minimap2Options) -> list[str]:
    exe = minimap2_executable()
    return [
        exe,
        "-t",
        str(max(1, int(opts.threads))),
        "-x",
        opts.preset,
        str(reference),
        str(query),
    ]


def run_minimap2(
    query_fasta: Path,
    reference_fasta: Path,
    opts: Minimap2Options,
) -> list[Path]:
    query = query_fasta.expanduser().resolve()
    reference = reference_fasta.expanduser().resolve()
    if not query.is_file():
        raise FileNotFoundError(query)
    if not reference.is_file():
        raise FileNotFoundError(reference)

    out_path = query.parent / f"{query.stem}.minimap2.paf"
    cmd = build_command(query, reference, opts)

    with out_path.open("w", encoding="utf-8") as handle:
        completed = subprocess.run(
            cmd,
            check=False,
            stdout=handle,
            stderr=subprocess.PIPE,
            text=True,
        )
    if completed.returncode != 0:
        detail = (completed.stderr or "").strip()
        raise RuntimeError(detail or f"minimap2 exited with code {completed.returncode}")
    if not out_path.is_file() or out_path.stat().st_size == 0:
        raise RuntimeError(f"minimap2 produced no output at {out_path}")
    return [out_path]
