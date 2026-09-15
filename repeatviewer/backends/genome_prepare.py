"""Prepare genome FASTA sidecar files (chromosome sizes + BLAST index)."""

from __future__ import annotations

import logging
import os
import platform
from pathlib import Path

from repeatviewer.backends.curation import genome_chrom_sizes
from repeatviewer.backends.gui_functions import blast_database_basename, check_database

logger = logging.getLogger(__name__)

_BLAST_NUCL_EXTENSIONS = (
    ".nhr",
    ".nin",
    ".nsq",
    ".ndb",
    ".not",
    ".nto",
    ".ntf",
    ".nog",
    ".nos",
)


def genome_chrom_sizes_path(genome_fasta: Path) -> Path:
    """Tab-separated chromosome lengths next to the FASTA (``*.chrom.sizes``).

    The sidecar basename uses underscores instead of spaces so tools that
    mishandle spaced paths can still open it.
    """
    genome_fasta = genome_fasta.expanduser().resolve()
    safe_name = genome_fasta.name.replace(" ", "_")
    return genome_fasta.parent / f"{safe_name}.chrom.sizes"


def _legacy_genome_chrom_sizes_path(genome_fasta: Path) -> Path:
    """Pre-sanitization sidecar name (kept for reading older files)."""
    genome_fasta = genome_fasta.expanduser().resolve()
    return genome_fasta.parent / f"{genome_fasta.name}.chrom.sizes"


def resolve_genome_chrom_sizes_file(genome_fasta: Path) -> Path | None:
    """Return an existing chrom.sizes sidecar, preferring the space-free name."""
    preferred = genome_chrom_sizes_path(genome_fasta)
    if preferred.is_file():
        return preferred
    legacy = _legacy_genome_chrom_sizes_path(genome_fasta)
    if legacy != preferred and legacy.is_file():
        return legacy
    return None


def _sidecar_is_stale(sidecar: Path, genome_fasta: Path) -> bool:
    if not sidecar.is_file():
        return True
    try:
        return sidecar.stat().st_mtime < genome_fasta.stat().st_mtime
    except OSError:
        return True


def write_genome_chrom_sizes_file(genome_fasta: Path) -> Path:
    """Write or refresh ``{genome}.chrom.sizes`` from FASTA sequence lengths."""
    genome_fasta = genome_fasta.expanduser().resolve()
    if not genome_fasta.is_file():
        raise FileNotFoundError(genome_fasta)
    sizes_path = genome_chrom_sizes_path(genome_fasta)
    sizes = genome_chrom_sizes(genome_fasta)
    lines = [f"{chrom}\t{length}" for chrom, length in sorted(sizes.items(), key=lambda kv: kv[0])]
    sizes_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    # Drop an older spaced sidecar so it is not listed as a problem path.
    legacy = _legacy_genome_chrom_sizes_path(genome_fasta)
    if legacy != sizes_path and legacy.is_file():
        try:
            legacy.unlink()
        except OSError:
            pass
    return sizes_path


def ensure_genome_chrom_sizes_file(genome_fasta: Path) -> Path:
    genome_fasta = genome_fasta.expanduser().resolve()
    preferred = genome_chrom_sizes_path(genome_fasta)
    existing = resolve_genome_chrom_sizes_file(genome_fasta)
    if existing is not None and not _sidecar_is_stale(existing, genome_fasta):
        if existing == preferred:
            return preferred
        # Migrate legacy spaced name → underscore name.
        try:
            existing.replace(preferred)
            return preferred
        except OSError:
            pass
    return write_genome_chrom_sizes_file(genome_fasta)


def _remove_blast_nucleotide_index(db_prefix: Path) -> None:
    base = str(db_prefix.expanduser().resolve())
    for ext in _BLAST_NUCL_EXTENSIONS:
        path = Path(base + ext)
        if path.is_file():
            path.unlink()


def ensure_genome_blast_index(
    genome_fasta: Path,
    *,
    use_system_blast: bool = False,
) -> str:
    """Return BLAST database prefix path, building the index when missing or stale."""
    genome_fasta = genome_fasta.expanduser().resolve()
    if not genome_fasta.is_file():
        raise FileNotFoundError(genome_fasta)
    # Index files use a sanitized basename (spaces → ``_``); do not look for
    # ``{original name with spaces}.nin``.
    db_prefix = genome_fasta.parent / blast_database_basename(genome_fasta)
    nin = Path(f"{db_prefix}.nin")
    if nin.is_file() and not _sidecar_is_stale(nin, genome_fasta):
        return str(db_prefix)
    _remove_blast_nucleotide_index(db_prefix)
    os_type = platform.system()
    result = check_database(
        str(genome_fasta),
        output_dir=None,
        os_type=os_type,
        use_system_blast=use_system_blast,
    )
    if result == "makeblastdb_not_found":
        raise RuntimeError(
            "makeblastdb was not found. Install BLAST+ or enable bundled BLAST in Settings."
        )
    if result == "makeblastdb_got_error":
        raise RuntimeError(
            "makeblastdb failed while indexing the genome FASTA. "
            "See the process log for the command output."
        )
    if not result:
        raise RuntimeError("Could not prepare BLAST index for the genome FASTA.")
    return str(result)


def ensure_genome_artifacts(
    genome_fasta: Path,
    *,
    use_system_blast: bool = False,
) -> tuple[bool, str]:
    """Ensure chromosome sizes file and BLAST nucleotide index exist for a genome."""
    genome_fasta = genome_fasta.expanduser().resolve()
    if not genome_fasta.is_file():
        return False, f"Genome file not found: {genome_fasta}"
    try:
        sizes_path = ensure_genome_chrom_sizes_file(genome_fasta)
        db_path = ensure_genome_blast_index(
            genome_fasta,
            use_system_blast=use_system_blast,
        )
    except Exception as exc:  # noqa: BLE001
        logger.warning("Genome prepare failed for %s: %s", genome_fasta, exc)
        return False, str(exc)
    return (
        True,
        f"Genome ready: {sizes_path.name} and BLAST index ({Path(db_path).name}.*)",
    )
