"""Resolve companion GFF files and Flexidot color configs for dotplot shading."""

from __future__ import annotations

import logging
from pathlib import Path

from Bio import SeqIO

from teviewer.backends.annotation_colors import feature_type_color
from teviewer.backends.gff_features import (
    _seqid_matches,
    find_companion_gff,
    list_feature_types_in_annotation,
    parse_gff,
)
from teviewer.backends.gui_functions import sanitize_name
from teviewer.backends.sequence_features import resolve_features_gff_path

logger = logging.getLogger(__name__)

_FLEXIDOT_GFF_ALPHA = 0.75
_FLEXIDOT_GFF_ALPHA_OTHER = 0.5
_FLEXIDOT_GFF_ZOOM = 0


def resolve_gff_for_fasta(fasta_path: Path) -> Path | None:
    """Return a companion GFF3/GFF for a FASTA, preferring sequence-map exports."""
    resolved = fasta_path.expanduser().resolve()
    gff = resolve_features_gff_path(resolved)
    if gff is not None:
        return gff
    return find_companion_gff(resolved)


def resolve_gff_paths_for_fastas(fasta_paths: list[Path]) -> list[Path]:
    """Unique existing GFF paths for the given FASTA files (stable order)."""
    found: list[Path] = []
    seen: set[str] = set()
    for fasta_path in fasta_paths:
        gff = resolve_gff_for_fasta(fasta_path)
        if gff is None:
            continue
        key = str(gff.resolve())
        if key in seen:
            continue
        seen.add(key)
        found.append(gff.resolve())
    return found


def build_flexidot_seqid_remap(fasta_paths: list[Path]) -> dict[str, str]:
    """
    Map original FASTA record IDs to Flexidot-safe IDs.

    Uses the same sanitization rules as combined Flexidot runs and
    ``prepare_flexidot_input_fasta``.
    """
    remap: dict[str, str] = {}
    seen: dict[str, int] = {}
    for fasta_path in fasta_paths:
        resolved = fasta_path.expanduser().resolve()
        if not resolved.is_file():
            continue
        for record in SeqIO.parse(str(resolved), "fasta"):
            original = str(record.id)
            safe = sanitize_name(original)
            count = seen.get(safe, 0)
            seen[safe] = count + 1
            final = f"{safe}_{count + 1}" if count else safe
            remap[original] = final
            base = original.split()[0]
            if base and base not in remap:
                remap[base] = final
    return remap


def _remap_gff_seqid(seqid: str, id_remap: dict[str, str]) -> str:
    if seqid in id_remap:
        return id_remap[seqid]
    for original, mapped in id_remap.items():
        if _seqid_matches(original, seqid):
            return mapped
    return seqid


def _gff_needs_remap(gff_path: Path, id_remap: dict[str, str]) -> bool:
    if not id_remap:
        return False
    try:
        with gff_path.open(encoding="utf-8", errors="replace") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 3:
                    continue
                seqid = cols[0]
                mapped = _remap_gff_seqid(seqid, id_remap)
                if mapped != seqid:
                    return True
    except OSError:
        return False
    return False


def write_remapped_gff(gff_path: Path, *, id_remap: dict[str, str], out_path: Path) -> Path:
    """Write a copy of ``gff_path`` with sequence IDs replaced for Flexidot."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    lines_out: list[str] = []
    with gff_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            stripped = line.rstrip("\n")
            if not stripped or stripped.startswith("#"):
                lines_out.append(stripped)
                continue
            cols = stripped.split("\t")
            if len(cols) < 9:
                lines_out.append(stripped)
                continue
            cols[0] = _remap_gff_seqid(cols[0], id_remap)
            lines_out.append("\t".join(cols))
    out_path.write_text("\n".join(lines_out) + "\n", encoding="utf-8")
    return out_path


def _feature_color_from_gff_row(feat_type: str, attributes: dict[str, str]) -> str:
    raw = str(attributes.get("color", "")).strip()
    if raw.startswith("#") and len(raw) >= 4:
        return raw
    return feature_type_color(feat_type)


def _flexidot_lookup_keys(feat_type: str, strand: str) -> list[str]:
    """Keys Flexidot uses in ``read_gffs`` (column 3 lowercased, ``_rev`` on minus)."""
    base = feat_type.lower()
    if strand == "-":
        return [f"{base}_rev", base]
    return [base, f"{base}_rev"]


def _collect_flexidot_gff_colors(gff_paths: list[Path]) -> dict[str, str]:
    """
    Build Flexidot ``-G`` color keys from GFF rows.

    Uses each feature's ``color=`` attribute when present (as written by TEviewer),
    otherwise the annotation color scheme. Reverse-strand rows use ``*_rev`` keys
  so Flexidot's built-in ``orf_rev`` defaults do not override TEviewer colors.
    """
    color_by_key: dict[str, str] = {}
    for gff_path in gff_paths:
        for feat in parse_gff(gff_path):
            feat_type = feat.attributes.get("type", feat.feature_type)
            color = _feature_color_from_gff_row(feat_type, feat.attributes)
            for key in _flexidot_lookup_keys(feat_type, feat.strand):
                color_by_key[key] = color

    if not color_by_key:
        for gff_path in gff_paths:
            for feat_type in list_feature_types_in_annotation(gff_path):
                if not feat_type:
                    continue
                color = feature_type_color(feat_type)
                for key in _flexidot_lookup_keys(feat_type, "+"):
                    color_by_key[key] = color

    return color_by_key


def write_flexidot_color_config(
    gff_paths: list[Path],
    out_path: Path,
    *,
    alpha: float = _FLEXIDOT_GFF_ALPHA,
    alpha_other: float = _FLEXIDOT_GFF_ALPHA_OTHER,
) -> Path:
    """
    Write a Flexidot ``-G`` tab-delimited color file using TEviewer feature colors.

    Colors are taken from GFF ``color=`` attributes when available, with matching
    ``*_rev`` entries for reverse-strand features.
    """
    color_by_key = _collect_flexidot_gff_colors(gff_paths)

    lines = [
        "# TEviewer Flexidot GFF color config",
        "# annotation_type\tcolor\ttransparency\tzoom",
    ]
    for key in sorted(color_by_key.keys(), key=str.lower):
        lines.append(f"{key}\t{color_by_key[key]}\t{alpha}\t{_FLEXIDOT_GFF_ZOOM}")
    others_color = feature_type_color("others")
    lines.append(f"others\t{others_color}\t{alpha_other}\t{_FLEXIDOT_GFF_ZOOM}")
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out_path


def prepare_flexidot_gff_inputs(
    fasta_paths: list[Path],
    out_dir: Path,
    *,
    output_stem: str,
) -> tuple[list[str], str, list[Path]] | None:
    """
    Resolve GFF files and write Flexidot ``-g`` / ``-G`` inputs.

    Returns ``(gff_cli_paths, color_config_path, temp_paths)`` or ``None`` when
    no companion GFF exists for any input FASTA.
  """
    gff_sources = resolve_gff_paths_for_fastas(fasta_paths)
    if not gff_sources:
        return None

    out_dir = out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    id_remap = build_flexidot_seqid_remap(fasta_paths)

    gff_cli_paths: list[str] = []
    temp_paths: list[Path] = []
    prepared_gff_paths: list[Path] = []

    for index, source in enumerate(gff_sources, start=1):
        if _gff_needs_remap(source, id_remap):
            dest = out_dir / f".teviewer_flexidot_{output_stem}_{index}.gff3"
            write_remapped_gff(source, id_remap=id_remap, out_path=dest)
            temp_paths.append(dest)
            prepared_gff_paths.append(dest)
            gff_cli_paths.append(str(dest))
        else:
            prepared_gff_paths.append(source)
            gff_cli_paths.append(str(source))

    color_config = out_dir / f".teviewer_flexidot_{output_stem}.gff.colors.config"
    write_flexidot_color_config(prepared_gff_paths, color_config)
    temp_paths.append(color_config)
    logger.info(
        "Flexidot GFF shading: %d file(s), temporary color config %s",
        len(gff_cli_paths),
        color_config.name,
    )
    return gff_cli_paths, str(color_config), temp_paths
