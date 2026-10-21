"""Batch TE annotation: run auto/similarity detectors (optionally write JSON)."""

from __future__ import annotations

import logging
import tempfile
import threading
from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from repeatviewer.backends.auto_annotation import (
    AutoAnnotationMode,
    _is_auto_generated_tsd,
    _next_track,
    _target_track_for_auto_rerun,
    merge_auto_annotation_features,
    run_auto_annotations,
)
from repeatviewer.backends.auto_annotation_settings import (
    AutoAnnotationSettings,
    default_annotation_threads,
    load_auto_annotation_settings,
)
from repeatviewer.backends.job_cancel import (
    JobCancelledError,
    JobCancellation,
    check_cancelled,
)
from repeatviewer.backends.sequence_features import (
    SequenceFeature,
    load_features_by_id_from_json,
    load_sequence_map_data,
    load_single_sequence_fasta,
    resolve_features_gff_path,
    resolve_sidecar_path,
    save_annotations_multi,
)
from repeatviewer.backends.thread_budget import load_batch_workers
from repeatviewer.backends.tsd_search import detect_tsd_features_by_id
from repeatviewer.job_options import TsdSearchOptions, load_tsd_search_options
from repeatviewer.backends.similarity_annotation import (
    SimilarityAnnotationSettings,
    SimilarityDatabaseRun,
    load_use_system_blast_from_settings,
    run_similarity_annotations,
)
from repeatviewer.backends.similarity_annotation_settings import load_similarity_annotation_settings


def _annotation_parallel_sequences(num_sequences: int) -> int:
    """Like BLAST: Max. parallel jobs capped by selection count and CPUs."""
    import os

    requested = load_batch_workers()
    sequences = max(1, int(num_sequences))
    if sequences <= 1 or requested <= 1:
        return 1
    cpus = os.cpu_count() or 4
    return max(1, min(requested, sequences, cpus))


@dataclass(frozen=True)
class TeAnnotationOptions:
    """Selected annotation modules for batch TE annotation."""

    auto_modes: tuple[AutoAnnotationMode, ...] = ()
    similarity_databases: tuple[SimilarityDatabaseRun, ...] = ()


@dataclass
class TeAnnotationBatchResult:
    """In-memory annotation result; sidecar paths only when ``persist=True``."""

    features_by_id: dict[str, list[SequenceFeature]]
    record_ids: list[str]
    lengths_by_id: dict[str, int]
    json_path: Path | None = None
    gff_path: Path | None = None
    # Sequences that failed detection but did not abort the batch: (id, reason).
    skipped_sequences: tuple[tuple[str, str], ...] = ()
    # Features newly present after this run (not in the pre-job baseline).
    identified_count: int = 0

    @property
    def output_paths(self) -> list[Path]:
        paths: list[Path] = []
        if self.json_path is not None:
            paths.append(self.json_path.resolve())
        if self.gff_path is not None:
            paths.append(self.gff_path.resolve())
        return paths

    def skipped_reminder_parts(
        self, *, kind: str = "Feature detection"
    ) -> tuple[str, str]:
        """Return ``(headline, body)`` for a skip reminder dialog.

        Use ``headline`` as ``QMessageBox`` text (native bold key line) and
        ``body`` as informative text (regular weight on macOS).
        """
        skipped = list(self.skipped_sequences)
        if not skipped:
            return "", ""
        headline = (
            f"{kind} finished, but {len(skipped)} sequence(s) were skipped"
        )
        lines: list[str] = []
        for rid, reason in skipped[:20]:
            detail = str(reason or "unknown error").strip() or "unknown error"
            if len(detail) > 180:
                detail = detail[:177] + "…"
            lines.append(f"• {rid}: {detail}")
        if len(skipped) > 20:
            lines.append(f"• …and {len(skipped) - 20} more")
        lines.append("")
        lines.append("Other sequences were annotated normally.")
        return headline, "\n".join(lines)

    def skipped_reminder_text(self, *, kind: str = "Feature detection") -> str:
        """Plain-text summary (tests / logs); prefer :meth:`skipped_reminder_parts` in UI."""
        headline, body = self.skipped_reminder_parts(kind=kind)
        if not headline:
            return ""
        return f"{headline}:\n\n{body}" if body else headline

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
    auto_threads: int | None = None,
    run_auto: bool = True,
    run_similarity: bool = True,
) -> list:
    """Run selected modules on one single-sequence FASTA; return feature list.

    Gap stripping and MSA coordinate projection are handled inside
    ``run_auto_annotations`` / ``run_similarity_annotations``.
    TSD is skipped here — it needs the full MSA (see ``_apply_msa_tsd``).
    """
    features: tuple = tuple(base_features)
    per_seq_modes = tuple(
        mode for mode in options.auto_modes if mode != AutoAnnotationMode.TSD
    )
    if run_auto and per_seq_modes:
        if progress is not None:
            progress(f"Auto annotation: {label}")
        auto_result = run_auto_annotations(
            record_fasta,
            per_seq_modes,
            cancel=cancel,
            settings=auto_cfg,
            base_features=features,
            progress=progress,
            num_threads=auto_threads,
        )
        features = auto_result.features
        if progress is not None and auto_result.message.strip():
            for line in auto_result.message.splitlines():
                text = line.strip()
                if text:
                    progress(text)

    if run_similarity and options.similarity_databases:
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
        if progress is not None and sim_result.message.strip():
            for line in sim_result.message.splitlines():
                text = line.strip()
                if text:
                    progress(text)
    return list(features)


def run_te_annotation_on_fasta(
    fasta_path: Path,
    options: TeAnnotationOptions,
    *,
    auto_settings: AutoAnnotationSettings | None = None,
    similarity_settings: SimilarityAnnotationSettings | None = None,
    use_system_blast: bool | None = None,
    cancel: JobCancellation | None = None,
    progress: Callable[[int, int, str], None] | None = None,
    record_ids: Sequence[str] | None = None,
    persist: bool = True,
    base_features_by_id: Mapping[str, Sequence[SequenceFeature]] | None = None,
    tsd_search_options: TsdSearchOptions | None = None,
) -> TeAnnotationBatchResult:
    """
    Run selected annotation modules on a FASTA (one or many sequences).

    Existing features are kept and new detections are added on top. By default
    results are written to the JSON sidecar (and companion GFF3 when written).
    Pass ``persist=False`` to keep results in memory only (Preview UI until Save).

    When ``base_features_by_id`` is set, that in-memory map is preferred over the
    on-disk JSON (so unsaved Preview edits are not discarded).

    When ``record_ids`` is set, only those sequences are annotated; other
    records keep their existing features unchanged.

    ``progress(completed, total, label)`` is called after each sequence finishes
    (and once with completed=0 when work starts).
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
    all_record_ids = list(data.record_ids)
    lengths = dict(data.lengths_by_id)
    if not all_record_ids:
        raise ValueError(f"No sequences in {resolved.name}")

    wanted: set[str] | None = None
    if record_ids is not None:
        wanted = {str(rid) for rid in record_ids if str(rid)}
        annotate_ids = [rid for rid in all_record_ids if rid in wanted]
        if not annotate_ids:
            raise ValueError("None of the selected sequences were found in the FASTA.")
    else:
        annotate_ids = list(all_record_ids)

    # Prefer caller-supplied in-memory features, then JSON, then load merge.
    json_feats = load_features_by_id_from_json(resolved, all_record_ids)
    existing_by_id: dict[str, list[SequenceFeature]] = {}
    for rid in all_record_ids:
        if base_features_by_id is not None and rid in base_features_by_id:
            existing_by_id[rid] = list(base_features_by_id.get(rid) or [])
        else:
            existing_by_id[rid] = list(
                json_feats.get(rid) or data.features_by_id.get(rid) or []
            )

    def _detail_progress(message: str) -> None:
        if progress is None:
            return
        # Keep the last completed count; detail-only updates use completed=-1 sentinel
        # via a nested helper — callers that only care about sequence steps ignore this.
        progress(-1, len(annotate_ids), message)

    features_by_id = {
        rid: list(existing_by_id.get(rid) or []) for rid in all_record_ids
    }
    skipped: list[tuple[str, str]] = []

    def _apply_msa_tsd() -> None:
        """Scan TSD on each MSA row and merge features (no cluster step)."""
        if AutoAnnotationMode.TSD not in options.auto_modes:
            return
        if len(all_record_ids) < 2:
            return
        check_cancelled(cancel)
        if progress is not None:
            progress(-1, len(annotate_ids), f"{resolved.name}: TSD scan")
        tsd_opts = tsd_search_options or load_tsd_search_options()
        tsd_by_id = detect_tsd_features_by_id(
            resolved,
            record_ids=all_record_ids,
            left_boundary=tsd_opts.left_boundary,
            right_boundary=tsd_opts.right_boundary,
            nt_radius=tsd_opts.nt_radius,
            k_min=tsd_opts.min_tsd_length,
            k_max=tsd_opts.max_tsd_length,
            min_distance=tsd_opts.min_distance,
            max_distance=tsd_opts.max_distance,
        )
        for rid in annotate_ids:
            existing = list(features_by_id.get(rid) or [])
            track = _target_track_for_auto_rerun(
                existing,
                AutoAnnotationMode.TSD,
                next_track=_next_track(existing),
            )
            existing = [feat for feat in existing if not _is_auto_generated_tsd(feat)]
            detected = list(tsd_by_id.get(rid) or [])
            if not detected:
                features_by_id[rid] = existing
                continue
            features_by_id[rid] = merge_auto_annotation_features(
                existing,
                detected,
                AutoAnnotationMode.TSD,
                track=track,
            )

    needs_per_sequence = bool(options.similarity_databases) or any(
        mode != AutoAnnotationMode.TSD for mode in options.auto_modes
    )

    # TSD-only on a multi-sequence MSA: skip per-row temp FASTA split.
    if not needs_per_sequence and len(all_record_ids) >= 2:
        if progress is not None:
            progress(0, 1, resolved.name)
        _apply_msa_tsd()
        if progress is not None:
            progress(1, 1, resolved.name)
    # Single-sequence path: annotate the file directly (no temp split).
    elif len(all_record_ids) == 1 and len(annotate_ids) == 1:
        rid = annotate_ids[0]
        if progress is not None:
            progress(0, 1, resolved.name)
        try:
            features = _annotate_one_record_fasta(
                resolved,
                options,
                auto_cfg=auto_cfg,
                sim_cfg=sim_cfg,
                use_blast=use_blast,
                base_features=tuple(existing_by_id.get(rid) or []),
                cancel=cancel,
                progress=_detail_progress,
                label=resolved.name,
                auto_threads=max(
                    1,
                    int(getattr(auto_cfg, "num_threads", default_annotation_threads())),
                ),
            )
            features_by_id[rid] = features
        except JobCancelledError:
            raise
        except Exception as exc:
            logging.exception("Annotation failed for %s", rid)
            skipped.append((rid, str(exc) or exc.__class__.__name__))
        if progress is not None:
            progress(1, 1, resolved.name)
        _apply_msa_tsd()
    else:
        records = list(SeqIO.parse(str(resolved), "fasta"))
        if len(records) != len(all_record_ids):
            # Uniquified ids vs raw headers: zip by index.
            records = records[: len(all_record_ids)]
        id_to_record = {
            rid: rec
            for rid, rec in zip(all_record_ids, records, strict=False)
        }
        total = len(annotate_ids)
        auto_threads = max(
            1, int(getattr(auto_cfg, "num_threads", default_annotation_threads()))
        )
        parallel = _annotation_parallel_sequences(len(annotate_ids))
        with tempfile.TemporaryDirectory(prefix="repeatviewer_teann_") as tmp:
            tmp_dir = Path(tmp)
            temp_by_id: dict[str, Path] = {}
            for index, rid in enumerate(annotate_ids):
                check_cancelled(cancel)
                record = id_to_record.get(rid)
                if record is None:
                    continue
                original_id = rid.split("#", 1)[0]
                temp_record = SeqRecord(
                    record.seq,
                    id=original_id,
                    name=original_id,
                    description=original_id,
                )
                temp_by_id[rid] = _write_single_record_temp_fasta(
                    temp_record, directory=tmp_dir, index=index
                )

            # BLAST-style: up to Max. parallel jobs sequences at once;
            # each sequence job uses its Threads setting for tools.
            completed = 0
            progress_lock = threading.Lock()

            def _run_one(
                rid: str, temp_fa: Path
            ) -> tuple[str, list[SequenceFeature] | None, str | None]:
                check_cancelled(cancel)
                try:
                    feats = _annotate_one_record_fasta(
                        temp_fa,
                        options,
                        auto_cfg=auto_cfg,
                        sim_cfg=sim_cfg,
                        use_blast=use_blast,
                        base_features=tuple(features_by_id.get(rid) or []),
                        cancel=cancel,
                        progress=None,
                        label=f"{resolved.name} [{rid}]",
                        auto_threads=auto_threads,
                    )
                    return rid, feats, None
                except JobCancelledError:
                    raise
                except Exception as exc:
                    logging.exception("Annotation failed for %s", rid)
                    return rid, None, str(exc) or exc.__class__.__name__

            if progress is not None:
                progress(
                    0,
                    total,
                    f"{resolved.name}: up to {parallel} sequence(s) in parallel "
                    f"({auto_threads} feature-detection thread(s) each)",
                )
            with ThreadPoolExecutor(max_workers=parallel) as pool:
                futures = [
                    pool.submit(_run_one, rid, temp_fa)
                    for rid, temp_fa in temp_by_id.items()
                ]
                for future in as_completed(futures):
                    check_cancelled(cancel)
                    try:
                        rid, feats, error = future.result()
                    except JobCancelledError:
                        for pending in futures:
                            pending.cancel()
                        raise
                    if error:
                        skipped.append((rid, error))
                    elif feats is not None:
                        features_by_id[rid] = feats
                    with progress_lock:
                        completed += 1
                        done = completed
                    if progress is not None:
                        if error:
                            progress(
                                done,
                                total,
                                f"{resolved.name}: skipped {done}/{total} ({rid})",
                            )
                        else:
                            progress(
                                done,
                                total,
                                f"{resolved.name}: sequence {done}/{total} ({rid})",
                            )
        _apply_msa_tsd()

    def _feature_signature(feat: SequenceFeature) -> tuple:
        return (
            str(getattr(feat, "feature_type", "") or ""),
            int(getattr(feat, "start", 0) or 0),
            int(getattr(feat, "end", 0) or 0),
            str(getattr(feat, "strand", "") or ""),
            str(getattr(feat, "name", "") or ""),
            int(getattr(feat, "track", 0) or 0),
        )

    baseline_sigs = {
        rid: {_feature_signature(feat) for feat in (existing_by_id.get(rid) or [])}
        for rid in annotate_ids
    }
    identified_count = 0
    for rid in annotate_ids:
        before = baseline_sigs.get(rid) or set()
        for feat in features_by_id.get(rid) or []:
            if _feature_signature(feat) not in before:
                identified_count += 1

    # Stream into the process monitor (Worker.log_line / FunctionWorker.log_line).
    from repeatviewer.external_process import emit_process_log

    emit_process_log(f"Identified features: {identified_count}")

    json_path: Path | None = None
    gff_path: Path | None = None
    if persist:
        json_path, gff_path = save_annotations_multi(
            resolved,
            features_by_id=features_by_id,
            lengths_by_id=lengths,
            record_ids=all_record_ids,
        )
        if json_path is not None:
            json_path = json_path.resolve()
        if gff_path is not None:
            gff_path = gff_path.resolve()
    return TeAnnotationBatchResult(
        features_by_id=features_by_id,
        record_ids=all_record_ids,
        lengths_by_id=lengths,
        json_path=json_path,
        gff_path=gff_path,
        skipped_sequences=tuple(skipped),
        identified_count=identified_count,
    )
