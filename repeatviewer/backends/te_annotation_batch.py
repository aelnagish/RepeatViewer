"""Batch TE annotation: run auto/similarity detectors and write companion JSON."""

from __future__ import annotations

import tempfile
from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from repeatviewer.backends.auto_annotation import AutoAnnotationMode, run_auto_annotations
from repeatviewer.backends.auto_annotation_settings import AutoAnnotationSettings, load_auto_annotation_settings
from repeatviewer.backends.job_cancel import JobCancellation, check_cancelled
from repeatviewer.backends.sequence_features import (
    load_features_by_id_from_json,
    load_sequence_map_data,
    load_single_sequence_fasta,
    resolve_features_gff_path,
    resolve_sidecar_path,
    save_annotations_multi,
)
from repeatviewer.backends.similarity_annotation import (
    SimilarityAnnotationSettings,
    SimilarityDatabaseRun,
    load_use_system_blast_from_settings,
    run_similarity_annotations,
)
from repeatviewer.backends.similarity_annotation_settings import load_similarity_annotation_settings


@dataclass(frozen=True)
class TeAnnotationOptions:
    """Selected annotation modules for batch TE annotation."""

    auto_modes: tuple[AutoAnnotationMode, ...] = ()
    similarity_databases: tuple[SimilarityDatabaseRun, ...] = ()


def existing_annotation_gff_paths(fasta_paths: list[Path]) -> list[Path]:
    """Return existing companion annotation GFF files (``.RV.gff3`` / legacy)."""
    found: list[Path] = []
    seen: set[Path] = set()
    for fasta in fasta_paths:
        gff = resolve_features_gff_path(fasta)
        if gff is None:
            continue
        resolved = gff.expanduser().resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        found.append(resolved)
    return found


def existing_annotation_json_paths(fasta_paths: list[Path]) -> list[Path]:
    """Return existing JSON annotation sidecars next to the FASTAs."""
    found: list[Path] = []
    seen: set[Path] = set()
    for fasta in fasta_paths:
        sidecar = resolve_sidecar_path(fasta)
        if sidecar is None:
            continue
        resolved = sidecar.expanduser().resolve()
        if resolved in seen:
            continue
        seen.add(resolved)
        found.append(resolved)
    return found


def is_single_sequence_fasta(path: Path) -> bool:
    """True when the FASTA contains exactly one sequence."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        return False
    try:
        load_single_sequence_fasta(resolved)
    except (OSError, ValueError):
        return False
    return True


def is_annotatable_fasta(path: Path) -> bool:
    """True when the path is a readable FASTA with at least one sequence."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        return False
    try:
        data = load_sequence_map_data(resolved)
    except (OSError, ValueError):
        return False
    return bool(data.record_ids)


def _write_single_record_temp_fasta(
    record: SeqRecord,
    *,
    directory: Path,
    index: int,
) -> Path:
    """Write one FASTA record for detectors that require a single-sequence file."""
    safe = "".join(c if c.isalnum() or c in "-_" else "_" for c in str(record.id))[:80]
    out = directory / f"seq_{index}_{safe or 'record'}.fa"
    SeqIO.write(record, str(out), "fasta")
    return out


def _annotate_one_record_fasta(
    record_fasta: Path,
    options: TeAnnotationOptions,
    *,
    auto_cfg: AutoAnnotationSettings,
    sim_cfg: SimilarityAnnotationSettings,
    use_blast: bool,
    base_features: tuple,
    cancel: JobCancellation | None,
    progress: Callable[[str], None] | None,
    label: str,
) -> list:
    """Run selected modules on one single-sequence FASTA; return feature list."""
    features: tuple = tuple(base_features)
    if options.auto_modes:
        if progress is not None:
            progress(f"Auto annotation: {label}")
        auto_result = run_auto_annotations(
            record_fasta,
            options.auto_modes,
            cancel=cancel,
            settings=auto_cfg,
            base_features=features,
            progress=progress,
        )
        features = auto_result.features

    if options.similarity_databases:
        if progress is not None:
            progress(f"Similarity annotation: {label}")
        sim_result = run_similarity_annotations(
            record_fasta,
            options.similarity_databases,
            settings=sim_cfg,
            use_system_blast=use_blast,
            orf_settings=auto_cfg.orf,
            cancel=cancel,
            base_features=features,
            progress=progress,
        )
        features = sim_result.features
    return list(features)


def run_te_annotation_on_fasta(
    fasta_path: Path,
    options: TeAnnotationOptions,
    *,
    auto_settings: AutoAnnotationSettings | None = None,
    similarity_settings: SimilarityAnnotationSettings | None = None,
    use_system_blast: bool | None = None,
    cancel: JobCancellation | None = None,
    progress: Callable[[str], None] | None = None,
) -> list[Path]:
    """
    Run selected annotation modules on a FASTA (one or many sequences).

    Existing JSON sidecar features are kept and new detections are added on top.
    Returns paths to the JSON sidecar (and companion GFF3 when written).
    """
    if not options.auto_modes and not options.similarity_databases:
        raise ValueError("Select at least one annotation module.")

    resolved = fasta_path.expanduser().resolve()
    if not resolved.is_file():
        raise FileNotFoundError(resolved)

    auto_cfg = auto_settings or load_auto_annotation_settings()
    sim_cfg = similarity_settings or load_similarity_annotation_settings()
    use_blast = (
        load_use_system_blast_from_settings()
        if use_system_blast is None
        else use_system_blast
    )

    data = load_sequence_map_data(resolved)
    record_ids = list(data.record_ids)
    lengths = dict(data.lengths_by_id)
    if not record_ids:
        raise ValueError(f"No sequences in {resolved.name}")

    # Prefer existing JSON; fall back to whatever load_sequence_map_data already merged.
    json_feats = load_features_by_id_from_json(resolved, record_ids)
    existing_by_id = {
        rid: list(json_feats.get(rid) or data.features_by_id.get(rid) or [])
        for rid in record_ids
    }

    # Single-sequence path: annotate the file directly (no temp split).
    if len(record_ids) == 1:
        rid = record_ids[0]
        features = _annotate_one_record_fasta(
            resolved,
            options,
            auto_cfg=auto_cfg,
            sim_cfg=sim_cfg,
            use_blast=use_blast,
            base_features=tuple(existing_by_id.get(rid) or []),
            cancel=cancel,
            progress=progress,
            label=resolved.name,
        )
        features_by_id = {rid: features}
    else:
        features_by_id = {rid: list(existing_by_id.get(rid) or []) for rid in record_ids}
        records = list(SeqIO.parse(str(resolved), "fasta"))
        if len(records) != len(record_ids):
            # Uniquified ids vs raw headers: zip by index.
            records = records[: len(record_ids)]
        total = len(record_ids)
        with tempfile.TemporaryDirectory(prefix="repeatviewer_teann_") as tmp:
            tmp_dir = Path(tmp)
            for index, (rid, record) in enumerate(
                zip(record_ids, records, strict=False)
            ):
                check_cancelled(cancel)
                if progress is not None:
                    progress(
                        f"{resolved.name}: sequence {index + 1}/{total} ({rid})"
                    )
                original_id = rid.split("#", 1)[0]
                temp_record = SeqRecord(
                    record.seq,
                    id=original_id,
                    name=original_id,
                    description=original_id,
                )
                temp_fa = _write_single_record_temp_fasta(
                    temp_record, directory=tmp_dir, index=index
                )
                new_feats = _annotate_one_record_fasta(
                    temp_fa,
                    options,
                    auto_cfg=auto_cfg,
                    sim_cfg=sim_cfg,
                    use_blast=use_blast,
                    base_features=tuple(features_by_id.get(rid) or []),
                    cancel=cancel,
                    progress=progress,
                    label=f"{resolved.name} [{rid}]",
                )
                features_by_id[rid] = new_feats

    json_path, gff_path = save_annotations_multi(
        resolved,
        features_by_id=features_by_id,
        lengths_by_id=lengths,
        record_ids=record_ids,
    )
    outputs = [json_path.resolve()]
    if gff_path is not None:
        outputs.append(gff_path.resolve())
    return outputs
