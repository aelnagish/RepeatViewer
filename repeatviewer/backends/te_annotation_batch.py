"""Batch TE annotation: run auto/similarity detectors and write companion GFF3 files."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass
from pathlib import Path

from repeatviewer.backends.auto_annotation import AutoAnnotationMode, run_auto_annotations
from repeatviewer.backends.auto_annotation_settings import AutoAnnotationSettings, load_auto_annotation_settings
from repeatviewer.backends.job_cancel import JobCancellation
from repeatviewer.backends.sequence_features import (
    load_sequence_map_data,
    load_single_sequence_fasta,
    save_annotations,
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
    Run the selected annotation modules on one single-sequence FASTA and save GFF3.

    Returns paths to the JSON sidecar and companion GFF3 (when features were found).
    """
    if not options.auto_modes and not options.similarity_databases:
        raise ValueError("Select at least one annotation module.")

    resolved = fasta_path.expanduser().resolve()
    if not is_single_sequence_fasta(resolved):
        raise ValueError(f"{resolved.name} is not a single-sequence FASTA file.")

    auto_cfg = auto_settings or load_auto_annotation_settings()
    sim_cfg = similarity_settings or load_similarity_annotation_settings()
    use_blast = (
        load_use_system_blast_from_settings()
        if use_system_blast is None
        else use_system_blast
    )

    features: tuple = ()
    if options.auto_modes:
        if progress is not None:
            progress(f"Auto annotation: {resolved.name}")
        auto_result = run_auto_annotations(
            resolved,
            options.auto_modes,
            cancel=cancel,
            settings=auto_cfg,
            base_features=(),
            progress=progress,
        )
        features = auto_result.features

    if options.similarity_databases:
        if progress is not None:
            progress(f"Similarity annotation: {resolved.name}")
        sim_result = run_similarity_annotations(
            resolved,
            options.similarity_databases,
            settings=sim_cfg,
            use_system_blast=use_blast,
            orf_settings=auto_cfg.orf,
            cancel=cancel,
            base_features=features,
            progress=progress,
        )
        features = sim_result.features

    data = load_sequence_map_data(resolved)
    json_path, gff_path = save_annotations(
        resolved,
        sequence_id=data.sequence_id,
        sequence_length=data.length,
        features=list(features),
    )
    outputs = [json_path.resolve()]
    if gff_path is not None:
        outputs.append(gff_path.resolve())
    return outputs
