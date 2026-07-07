from __future__ import annotations

import subprocess
from pathlib import Path

from teviewer_minimap2.executable import minimap2_executable
from teviewer_minimap2.options import Minimap2Options


def build_command(query: Path, reference: Path, opts: Minimap2Options) -> list[str]:
    exe = minimap2_executable()
    cmd = [
        exe,
        "-t",
        str(max(1, int(opts.threads))),
        "-x",
        opts.preset,
    ]
    if opts.output_format == "sam":
        cmd.append("-a")
    cmd.extend([str(reference), str(query)])
    return cmd


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

    suffix = "sam" if opts.output_format == "sam" else "paf"
    out_path = query.parent / f"{query.stem}.minimap2.{suffix}"
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
        raise RuntimeError(
            f"minimap2 failed (exit {completed.returncode})"
            + (f":\n{detail}" if detail else "")
        )
    if not out_path.is_file() or out_path.stat().st_size == 0:
        raise RuntimeError("minimap2 finished but produced no output.")
    return [out_path]
