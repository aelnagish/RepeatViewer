"""Dispatch auto-annotation jobs (LTR, TIR, tandem repeat, ORF, TSD)."""

from __future__ import annotations

import re
from collections.abc import Callable
from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from Bio import SeqIO

from repeatviewer.backends.auto_annotation_settings import AutoAnnotationSettings
from repeatviewer.backends.job_cancel import JobCancellation
from repeatviewer.backends.orf_detect import detect_orf_features
from repeatviewer.backends.sequence_features import (
    REPEATVIEWER_INTERNAL_SOURCE,
    SequenceFeature,
)
from repeatviewer.backends.tandem_repeat_detect import detect_tandem_repeat_features
from repeatviewer.backends.terminal_repeat_detect import auto_annotate_terminal_repeats


class AutoAnnotationMode(str, Enum):
    LTR = "LTR"
    TIR = "TIR"
    TANDEM = "tandem repeat"
    ORF = "ORF"
    TSD = "TSD"


_AUTO_ORF_NAME = re.compile(r"^ORF\d+$")
_AUTO_TSD_NAME = re.compile(r"^TSD_")


@dataclass(frozen=True)
class AutoAnnotationResult:
    features: tuple[SequenceFeature, ...]
    mode: AutoAnnotationMode
    message: str


@dataclass(frozen=True)
class MultiAutoAnnotationResult:
    features: tuple[SequenceFeature, ...]
    modes: tuple[AutoAnnotationMode, ...]
    message: str


def _next_track(features: list[SequenceFeature]) -> int:
    if not features:
        return 0
    return max(feat.track for feat in features) + 1


def _auto_mode_feature_type(mode: AutoAnnotationMode) -> str:
    if mode == AutoAnnotationMode.TANDEM:
        return "tandem repeat"
    return mode.value


def _feature_signature(feat: SequenceFeature) -> tuple[int, int, str, str]:
    return (feat.start, feat.end, feat.name, feat.strand)


def _features_for_mode(
    features: list[SequenceFeature],
    mode: AutoAnnotationMode,
) -> list[SequenceFeature]:
    feature_type = _auto_mode_feature_type(mode)
    return [feat for feat in features if feat.feature_type == feature_type]


def _remove_features_for_mode(
    features: list[SequenceFeature],
    mode: AutoAnnotationMode,
) -> list[SequenceFeature]:
    feature_type = _auto_mode_feature_type(mode)
    return [feat for feat in features if feat.feature_type != feature_type]


def _is_auto_generated_orf(feat: SequenceFeature) -> bool:
    return feat.feature_type == "ORF" and bool(_AUTO_ORF_NAME.match(feat.name or ""))


def _is_auto_generated_tsd(feat: SequenceFeature) -> bool:
    return feat.feature_type == "TSD" and bool(_AUTO_TSD_NAME.match(feat.name or ""))


def _remove_features_for_auto_rerun(
    features: list[SequenceFeature],
    mode: AutoAnnotationMode,
) -> list[SequenceFeature]:
    if mode == AutoAnnotationMode.ORF:
        return [feat for feat in features if not _is_auto_generated_orf(feat)]
    if mode == AutoAnnotationMode.TSD:
        return [feat for feat in features if not _is_auto_generated_tsd(feat)]
    feature_type = _auto_mode_feature_type(mode)
    has_manual = any(feat.track == 0 for feat in _features_for_mode(features, mode))
    if has_manual:
        return [
            feat
            for feat in features
            if not (feat.feature_type == feature_type and feat.track != 0)
        ]
    return _remove_features_for_mode(features, mode)


def _target_track_for_auto_rerun(
    features: list[SequenceFeature],
    mode: AutoAnnotationMode,
    *,
    next_track: int,
) -> int:
    existing = _features_for_mode(features, mode)
    if mode == AutoAnnotationMode.ORF:
        auto_existing = [feat for feat in existing if _is_auto_generated_orf(feat)]
        manual_existing = [feat for feat in existing if not _is_auto_generated_orf(feat)]
        if auto_existing:
            return min(feat.track for feat in auto_existing)
        if manual_existing:
            return _next_track(features)
        return next_track
    if mode == AutoAnnotationMode.TSD:
        auto_existing = [feat for feat in existing if _is_auto_generated_tsd(feat)]
        manual_existing = [feat for feat in existing if not _is_auto_generated_tsd(feat)]
        if auto_existing:
            return min(feat.track for feat in auto_existing)
        if manual_existing:
            return _next_track(features)
        return next_track
    has_manual = any(feat.track == 0 for feat in existing)
    if has_manual:
        auto_existing = [feat for feat in existing if feat.track != 0]
        if auto_existing:
            return min(feat.track for feat in auto_existing)
        return _next_track(features)
    if existing:
        return min(feat.track for feat in existing)
    return next_track


def _stamp_track(
    features: list[SequenceFeature] | tuple[SequenceFeature, ...],
    track: int,
) -> list[SequenceFeature]:
    stamped: list[SequenceFeature] = []
    for feat in features:
        stamped.append(
            SequenceFeature(
                start=feat.start,
                end=feat.end,
                name=feat.name,
                feature_type=feat.feature_type,
                color=feat.color,
                strand=feat.strand,
                track=track,
                parent=feat.parent,
                source=REPEATVIEWER_INTERNAL_SOURCE,
            )
        )
    return stamped


def run_auto_annotations(
    fasta_path: Path,
    modes: tuple[AutoAnnotationMode, ...],
    *,
    cancel: JobCancellation | None = None,
    settings: AutoAnnotationSettings | None = None,
    base_features: tuple[SequenceFeature, ...] | None = None,
    progress: Callable[[str], None] | None = None,
    num_threads: int | None = None,
) -> MultiAutoAnnotationResult:
    """Run one or more auto-annotation modules and append results on new tracks.

    MSA rows with gaps are ungapped for detection; coordinates are projected
    back to alignment columns in the returned features.
    """
    ordered = tuple(modes)
    if not ordered:
        return MultiAutoAnnotationResult((), (), "No auto-annotation modules selected.")

    from repeatviewer.backends.msa_coordinates import prepare_msa_annotation_workspace
    from repeatviewer.backends.sequence_features import load_sequence_map_data

    cfg = settings or AutoAnnotationSettings()
    threads = max(
        1,
        int(num_threads if num_threads is not None else cfg.num_threads),
    )
    data = load_sequence_map_data(fasta_path)
    incoming = tuple(base_features if base_features is not None else data.features)
    workspace = prepare_msa_annotation_workspace(fasta_path, incoming)
    try:
        work_path = workspace.work_fasta
        features = list(workspace.base_features)
        messages: list[str] = []
        if workspace.keep is not None and progress is not None:
            progress("MSA: removed gaps for feature detection")
        next_track = _next_track(features)
        total = len(ordered)
        if progress is not None:
            progress(f"0/{total} starting")
        for index, mode in enumerate(ordered, start=1):
            if progress is not None:
                # Completed count before this mode: index-1.
                progress(f"{index - 1}/{total} {mode.value}")
            result = run_auto_annotation(
                work_path,
                mode,
                cancel=cancel,
                settings=cfg,
                num_threads=threads,
            )
            if progress is not None:
                progress(f"{index}/{total} {mode.value}")
            if result.message.strip():
                messages.append(result.message.strip())
            existing = _features_for_mode(features, mode)
            new_sigs = {_feature_signature(feat) for feat in result.features}
            if existing and new_sigs == {_feature_signature(feat) for feat in existing}:
                messages.append(
                    f"{mode.value}: unchanged ({len(result.features)} feature(s))."
                )
                continue
            target_track = _target_track_for_auto_rerun(
                features,
                mode,
                next_track=next_track,
            )
            features = _remove_features_for_auto_rerun(features, mode)
            if result.features:
                features = merge_auto_annotation_features(
                    features,
                    result.features,
                    result.mode,
                    track=target_track,
                )
                if not existing:
                    next_track += 1
        summary = "\n".join(messages) if messages else "Auto annotation complete."
        return MultiAutoAnnotationResult(
            features=tuple(workspace.project_features(features)),
            modes=ordered,
            message=summary,
        )
    finally:
        workspace.close()


def merge_auto_annotation_features(
    existing: list[SequenceFeature],
    detected: list[SequenceFeature] | tuple[SequenceFeature, ...],
    mode: AutoAnnotationMode,
    *,
    track: int | None = None,
) -> list[SequenceFeature]:
    """Append detected features on a new track without removing existing annotations."""
    del mode  # kept for call-site clarity
    if not detected:
        return list(existing)
    target_track = _next_track(existing) if track is None else max(0, int(track))
    return list(existing) + _stamp_track(detected, target_track)


def run_auto_annotation(
    fasta_path: Path,
    mode: AutoAnnotationMode,
    *,
    cancel: JobCancellation | None = None,
    settings: AutoAnnotationSettings | None = None,
    num_threads: int | None = None,
) -> AutoAnnotationResult:
    """Run the selected auto-annotation module on a single-sequence FASTA."""
    fasta_path = fasta_path.expanduser().resolve()
    cfg = settings or AutoAnnotationSettings()
    threads = max(
        1,
        int(num_threads if num_threads is not None else cfg.num_threads),
    )

    if mode == AutoAnnotationMode.TANDEM:
        record = SeqIO.read(str(fasta_path), "fasta")
        tandem = detect_tandem_repeat_features(
            record.id,
            str(record.seq),
            settings=cfg.tandem,
        )
        return AutoAnnotationResult(
            features=tandem.features,
            mode=mode,
            message=tandem.message,
        )

    if mode == AutoAnnotationMode.ORF:
        record = SeqIO.read(str(fasta_path), "fasta")
        orfs = detect_orf_features(
            record.id,
            str(record.seq),
            settings=cfg.orf,
        )
        return AutoAnnotationResult(
            features=orfs.features,
            mode=mode,
            message=orfs.message,
        )

    if mode == AutoAnnotationMode.TSD:
        # TSD needs the full MSA; batch annotation runs it on the multi-FASTA.
        return AutoAnnotationResult(
            features=(),
            mode=mode,
            message="TSD detection runs on the full MSA (see batch annotation).",
        )

    terminal = auto_annotate_terminal_repeats(
        fasta_path,
        cancel=cancel,
        repeat_mode=mode.value,
        ltr_settings=cfg.ltr,
        tir_settings=cfg.tir,
        num_threads=threads,
    )
    return AutoAnnotationResult(
        features=terminal.features,
        mode=mode,
        message=terminal.message,
    )
