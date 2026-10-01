"""Helpers for tool output paths beside input FASTA files."""

from __future__ import annotations

import os
from pathlib import Path


def kept_input_name(path: Path | str) -> str:
    """Basename kept in derived tool outputs (same convention as BLAST).

    Example: ``reads.fa`` → ``reads.fa_mafft.fa``, ``reads.fa_flexidot.png``.
    Always use the full filename (including extension), never ``Path.stem``.
    """
    return Path(path).name


def write_genome_reference_link(
    out_dir: Path,
    genome_fasta: Path,
    *,
    genome_label: str = "",
) -> Path:
    """Symlink to the BLAST database FASTA for Flexidot multi-select in the output folder.

    Returns the link/copy path beside the outputs (not the resolved genome target),
    so downstream tools do not treat the genome as a BLAST-hit FASTA.

    Safe for parallel BLAST workers sharing one ``out_dir``: if the link already
    points at the same genome, it is left alone (avoids ``SameFileError`` races).
    """
    from repeatviewer.pan_genome import genome_output_suffix

    genome_fasta = genome_fasta.expanduser().resolve()
    out_dir = out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    suffix = genome_output_suffix(genome_label)
    if suffix:
        link_name = f"{suffix.lstrip('_')}.lnk.fa"
    else:
        link_name = f"{kept_input_name(genome_fasta)}.lnk.fa"

    link_path = out_dir / link_name

    def _points_to_genome() -> bool:
        if not (link_path.exists() or link_path.is_symlink()):
            return False
        try:
            return os.path.samefile(link_path, genome_fasta)
        except OSError:
            return False

    # Parallel jobs share this path — reuse if another worker already linked it.
    if _points_to_genome():
        return link_path

    if link_path.exists() or link_path.is_symlink():
        try:
            link_path.unlink()
        except FileNotFoundError:
            pass

    try:
        link_path.symlink_to(genome_fasta)
    except OSError:
        if _points_to_genome():
            return link_path
        # Relative symlink helps when source/target share a parent directory.
        try:
            link_path.symlink_to(os.path.relpath(genome_fasta, out_dir))
        except OSError:
            if _points_to_genome():
                return link_path
            import shutil

            try:
                shutil.copy2(genome_fasta, link_path)
            except shutil.SameFileError:
                # Another worker created a symlink/hardlink to the genome.
                pass
            except OSError:
                if not _points_to_genome():
                    raise
    # Do not Path.resolve() here — that follows the symlink to the genome FASTA.
    return link_path
