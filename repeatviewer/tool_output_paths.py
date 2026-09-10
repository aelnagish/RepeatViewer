"""Helpers for tool output paths beside input FASTA files."""

from __future__ import annotations

import os
from pathlib import Path


def write_genome_reference_link(
    out_dir: Path,
    genome_fasta: Path,
    *,
    genome_label: str = "",
) -> Path:
    """Symlink to the BLAST database FASTA for Flexidot multi-select in the output folder.

    Returns the link/copy path beside the outputs (not the resolved genome target),
    so downstream tools do not treat the genome as a BLAST-hit FASTA.
    """
    from repeatviewer.pan_genome import genome_output_suffix

    genome_fasta = genome_fasta.expanduser().resolve()
    out_dir = out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    suffix = genome_output_suffix(genome_label)
    if suffix:
        link_name = f"{suffix.lstrip('_')}.lnk.fa"
    else:
        link_name = f"{genome_fasta.stem}.lnk.fa"

    link_path = out_dir / link_name
    if link_path.exists() or link_path.is_symlink():
        link_path.unlink()

    try:
        link_path.symlink_to(genome_fasta)
    except OSError:
        # Relative symlink helps when source/target share a parent directory.
        try:
            link_path.symlink_to(os.path.relpath(genome_fasta, out_dir))
        except OSError:
            import shutil

            shutil.copy2(genome_fasta, link_path)
    # Do not Path.resolve() here — that follows the symlink to the genome FASTA.
    return link_path
