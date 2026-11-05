"""Batch TE annotation: run auto/similarity detectors (optionally write JSON)."""

from __future__ import annotations

import logging
import tempfile
import threading
from collections import Counter
from collections.abc import Callable, Mapping, Sequence
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from repeatviewer.backends.auto_annotation import (
    AutoAnnotationMode,
    _feature_coord_key,
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
from repeatviewer.backends.tsd_search import (
    detect_tsd_features_by_id,
    tsd_cluster_report_path,
)
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
    tsd_report_path: Path | None = None
    # Sequences that failed detection but did not abort the batch: (id, reason).
    skipped_sequences: tuple[tuple[str, str], ...] = ()
    # Features newly present after this run (not in the pre-job baseline).
    identified_count: int = 0
    # New-feature counts by ``feature_type`` (e.g. LTR, PBS, domain).
    identified_by_type: dict[str, int] = field(default_factory=dict)
    # Detector totals this run (includes hits that already matched existing).
    found_by_type: dict[str, int] = field(default_factory=dict)
    overlap_by_type: dict[str, int] = field(default_factory=dict)
    new_by_type: dict[str, int] = field(default_factory=dict)
    # How many sequences were actually annotated in this job.
    annotated_sequence_count: int = 0

    @property
    def output_paths(self) -> list[Path]:
        paths: list[Path] = []
        if self.json_path is not None:
            paths.append(self.json_path.resolve())
        if self.gff_path is not None:
            paths.append(self.gff_path.resolve())
        if self.tsd_report_path is not None:
            paths.append(self.tsd_report_path.resolve())
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
) -> tuple[list, dict[str, int], dict[str, int], dict[str, int]]:
    """Run selected modules on one single-sequence FASTA.

    Returns ``(features, found_by_type, overlap_by_type, new_by_type)``.
    Gap stripping and MSA coordinate projection are handled inside
    ``run_auto_annotations`` / ``run_similarity_annotations``.
    TSD is skipped here — it needs the full MSA (see ``_apply_msa_tsd``).
    """
    features: tuple = tuple(base_features)
    found_by_type: Counter[str] = Counter()
    overlap_by_type: Counter[str] = Counter()
    new_by_type: Counter[str] = Counter()
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
        found_by_type.update(auto_result.found_by_type)
        overlap_by_type.update(auto_result.overlap_by_type)
        new_by_type.update(auto_result.new_by_type)
        if progress is not None and auto_result.message.strip():
            for line in auto_result.message.splitlines():
                text = line.strip()
                if text:
                    progress(text)

    if run_similarity and options.similarity_databases:
        if progress is not None:
            progress(f"Similarity annotation: {label}")
        before_coords = {_feature_coord_key(feat) for feat in features}
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
        for feat in features:
            key = _feature_coord_key(feat)
            if key in before_coords:
                continue
            ftype = key[0]
            found_by_type[ftype] += 1
            new_by_type[ftype] += 1
        if progress is not None and sim_result.message.strip():
            for line in sim_result.message.splitlines():
                text = line.strip()
                if text:
                    progress(text)
    return (
        list(features),
        dict(found_by_type),
        dict(overlap_by_type),
        dict(new_by_type),
    )


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
    found_by_type: Counter[str] = Counter()
    overlap_by_type: Counter[str] = Counter()
    new_by_type: Counter[str] = Counter()
    stats_lock = threading.Lock()

    def _accumulate_stats(
        found: Mapping[str, int],
        overlap: Mapping[str, int],
        new: Mapping[str, int],
    ) -> None:
        with stats_lock:
            found_by_type.update(
                {k: v for k, v in dict(found or {}).items() if str(k).strip()}
            )
            overlap_by_type.update(
                {k: v for k, v in dict(overlap or {}).items() if str(k).strip()}
            )
            new_by_type.update(
                {k: v for k, v in dict(new or {}).items() if str(k).strip()}
            )

    tsd_report_holder: list[Path] = []

    def _apply_msa_tsd(*, drive_progress: bool = False) -> None:
        """Scan TSD on the MSA, annotate top cluster, write clustered text report.

        When ``drive_progress`` is True (TSD-only Feature detection), report
        per-row scan progress against ``len(annotate_ids)`` so the Overview bar
        advances. Otherwise emit detail-only chatter (TSD runs first; other
        detectors then own the sequence-level progress channel).
        """
        if AutoAnnotationMode.TSD not in options.auto_modes:
            return
        if len(all_record_ids) < 2:
            return
        check_cancelled(cancel)
        expected_total = max(1, len(annotate_ids))

        def _tsd_progress(done: int, total: int, label: str) -> None:
            if progress is None:
                return
            if drive_progress:
                # Map scan rows onto the Overview's expected N (annotate_ids).
                mapped = (
                    int(round(done * expected_total / max(1, total)))
                    if total > 0
                    else 0
                )
                mapped = max(0, min(expected_total, mapped))
                progress(mapped, expected_total, f"TSD {label}")
            else:
                progress(
                    -1,
                    expected_total,
                    f"{resolved.name}: TSD {done}/{total} {label}",
                )

        if progress is not None and not drive_progress:
            # Mixed Feature detection: TSD runs before per-sequence detectors.
            progress(
                -1,
                expected_total,
                f"{resolved.name}: TSD scan (before other detectors)",
            )
        tsd_opts = tsd_search_options or load_tsd_search_options()
        # Optional clustered text report (off by default; TSD settings).
        report_path = (
            tsd_cluster_report_path(resolved)
            if bool(getattr(tsd_opts, "write_text_report", False))
            else None
        )
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
            max_sim_threshold=tsd_opts.max_sim_threshold,
            nt_threshold=tsd_opts.nt_threshold,
            cluster_eps=tsd_opts.cluster_eps,
            cluster_report_path=report_path,
            progress=_tsd_progress,
        )
        if report_path is not None and report_path.is_file():
            tsd_report_holder.append(report_path.resolve())
        for rid in annotate_ids:
            existing = list(features_by_id.get(rid) or [])
            track = _target_track_for_auto_rerun(
                existing,
                AutoAnnotationMode.TSD,
                next_track=_next_track(existing),
            )
            # Overlap vs what was on the map before this TSD pass.
            prior_coords = {
                _feature_coord_key(feat)
                for feat in existing
                if str(getattr(feat, "feature_type", "") or "") == "TSD"
            }
            existing = [feat for feat in existing if not _is_auto_generated_tsd(feat)]
            detected = list(tsd_by_id.get(rid) or [])
            tsd_found: Counter[str] = Counter()
            tsd_overlap: Counter[str] = Counter()
            tsd_new: Counter[str] = Counter()
            for feat in detected:
                key = _feature_coord_key(feat)
                tsd_found["TSD"] += 1
                if key in prior_coords:
                    tsd_overlap["TSD"] += 1
                else:
                    tsd_new["TSD"] += 1
            _accumulate_stats(tsd_found, tsd_overlap, tsd_new)
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

    # MSA-wide TSD first when selected (needs full alignment; other detectors
    # then run per sequence and keep these TSD features as base).
    tsd_selected = AutoAnnotationMode.TSD in options.auto_modes and len(all_record_ids) >= 2
    if tsd_selected and not needs_per_sequence:
        n = max(1, len(annotate_ids))
        if progress is not None:
            progress(0, n, f"{resolved.name}: TSD scan")
        _apply_msa_tsd(drive_progress=True)
        if progress is not None:
            progress(n, n, resolved.name)
    elif tsd_selected:
        _apply_msa_tsd(drive_progress=False)

    # Single-sequence path: annotate the file directly (no temp split).
    if needs_per_sequence and len(all_record_ids) == 1 and len(annotate_ids) == 1:
        rid = annotate_ids[0]
        if progress is not None:
            progress(0, 1, resolved.name)
        try:
            features, found, overlap, new = _annotate_one_record_fasta(
                resolved,
                options,
                auto_cfg=auto_cfg,
                sim_cfg=sim_cfg,
                use_blast=use_blast,
                base_features=tuple(features_by_id.get(rid) or []),
                cancel=cancel,
                progress=_detail_progress,
                label=resolved.name,
                auto_threads=max(
                    1,
                    int(getattr(auto_cfg, "num_threads", default_annotation_threads())),
                ),
            )
            features_by_id[rid] = features
            _accumulate_stats(found, overlap, new)
        except JobCancelledError:
            raise
        except Exception as exc:
            logging.exception("Annotation failed for %s", rid)
            skipped.append((rid, str(exc) or exc.__class__.__name__))
        if progress is not None:
            progress(1, 1, resolved.name)
    elif needs_per_sequence:
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
            ) -> tuple[
                str,
                list[SequenceFeature] | None,
                dict[str, int],
                dict[str, int],
                dict[str, int],
                str | None,
            ]:
                check_cancelled(cancel)
                try:
                    feats, found, overlap, new = _annotate_one_record_fasta(
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
                    return rid, feats, found, overlap, new, None
                except JobCancelledError:
                    raise
                except Exception as exc:
                    logging.exception("Annotation failed for %s", rid)
                    return rid, None, {}, {}, {}, str(exc) or exc.__class__.__name__

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
                        rid, feats, found, overlap, new, error = future.result()
                    except JobCancelledError:
                        for pending in futures:
                            pending.cancel()
                        raise
                    if error:
                        skipped.append((rid, error))
                    elif feats is not None:
                        features_by_id[rid] = feats
                        _accumulate_stats(found, overlap, new)
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
    identified_by_type: Counter[str] = Counter()
    for rid in annotate_ids:
        before = baseline_sigs.get(rid) or set()
        for feat in features_by_id.get(rid) or []:
            if _feature_signature(feat) not in before:
                identified_count += 1
                ftype = str(getattr(feat, "feature_type", "") or "").strip() or "others"
                identified_by_type[ftype] += 1

    # If detectors reported no found map (e.g. unexpected path), fall back to
    # baseline-diff counts so the summary still has numbers.
    if not found_by_type and identified_by_type:
        found_by_type.update(identified_by_type)
        new_by_type.update(identified_by_type)
    elif not new_by_type and identified_by_type:
        # Prefer detector "new" when present; else baseline-diff.
        new_by_type.update(identified_by_type)

    preferred_types = [
        str(mode.value) for mode in options.auto_modes if str(mode.value or "").strip()
    ]
    # Similarity hits are usually ``domain`` (or similar); include when that ran.
    if options.similarity_databases and "domain" not in preferred_types:
        if (
            identified_by_type.get("domain", 0)
            or found_by_type.get("domain", 0)
            or options.similarity_databases
        ):
            preferred_types.append("domain")

    # Stream into the process monitor (Worker.log_line / FunctionWorker.log_line).
    from repeatviewer.external_process import emit_process_log
    from repeatviewer.widgets.process_monitor_dialog import format_annotation_run_summary

    for line in format_annotation_run_summary(
        n_sequences=len(annotate_ids),
        counts_by_type=dict(new_by_type),
        preferred_types=preferred_types,
        total_new=sum(new_by_type.values()),
        found_by_type=dict(found_by_type),
        overlap_by_type=dict(overlap_by_type),
    ):
        emit_process_log(line)

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
        tsd_report_path=tsd_report_holder[0] if tsd_report_holder else None,
        skipped_sequences=tuple(skipped),
        identified_count=identified_count,
        identified_by_type=dict(identified_by_type),
        found_by_type=dict(found_by_type),
        overlap_by_type=dict(overlap_by_type),
        new_by_type=dict(new_by_type),
        annotated_sequence_count=len(annotate_ids),
    )
