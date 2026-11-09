"""TSD (target site duplication) search in aligned FASTA MSAs (from TEtrimmer TSD_search.py)."""

from __future__ import annotations

import collections
import os
from collections import defaultdict
from collections.abc import Callable
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path

import numpy as np
from Bio import AlignIO
from sklearn.cluster import DBSCAN

from repeatviewer.backends.sequence_features import (
    REPEATVIEWER_INTERNAL_SOURCE,
    SequenceFeature,
)

# ``progress(completed, total, label)`` — completed rows / phase label.
TsdProgressCallback = Callable[[int, int, str], None]

# Allow TE edges (next nt after left TSD / prev nt before right TSD) to sit
# within this many MSA columns of the defined left/right boundaries.
TSD_BOUNDARY_MATCH_TOLERANCE_BP = 5


def _te_edges_match_boundaries(
    te_start_0: int,
    te_end_0: int,
    left_boundary_1: int,
    right_boundary_1: int,
    *,
    tolerance: int = TSD_BOUNDARY_MATCH_TOLERANCE_BP,
) -> bool:
    """True when TE body edges abut the defined boundaries (± ``tolerance`` bp).

    ``te_start_0`` / ``te_end_0`` are 0-based MSA columns (first TE nt after the
    left TSD, last TE nt before the right TSD). Boundaries are 1-based.
    """
    left_1 = int(left_boundary_1)
    right_1 = int(right_boundary_1)
    if left_1 <= 0 or right_1 <= 0:
        return False
    tol = max(0, int(tolerance))
    left_ok = abs((int(te_start_0) + 1) - left_1) <= tol
    right_ok = abs((int(te_end_0) + 1) - right_1) <= tol
    return left_ok and right_ok


def _prefer_tsd_candidate(
    prev: dict[str, object] | None,
    candidate: dict[str, object],
) -> dict[str, object]:
    """Keep the longer motif; tie-break toward the leftmost TE start."""
    if prev is None:
        return candidate
    prev_len = len(str(prev["tsd_seq"]))
    new_len = len(str(candidate["tsd_seq"]))
    if new_len > prev_len or (
        new_len == prev_len and int(candidate["te_start"]) < int(prev["te_start"])
    ):
        return candidate
    return prev


def _tsd_features_from_hit(
    sequence: str,
    *,
    te_start: int,
    te_end: int,
    tsd_seq: str,
    k_min: int,
) -> list[SequenceFeature] | None:
    """Build left/right TSD annotation features for one hit, or ``None``."""
    motif = str(tsd_seq or "").upper()
    if len(motif) < k_min:
        return None
    left_span = _tsd_nucleotide_span(sequence, int(te_start), len(motif), "left")
    right_span = _tsd_nucleotide_span(sequence, int(te_end), len(motif), "right")
    if left_span is None or right_span is None:
        return None
    if left_span[2] != motif or right_span[2] != motif:
        return None
    name = f"TSD_{motif}"
    return [
        SequenceFeature(
            start=left_span[0],
            end=left_span[1],
            name=name,
            feature_type="TSD",
            strand="+",
            source=REPEATVIEWER_INTERNAL_SOURCE,
        ),
        SequenceFeature(
            start=right_span[0],
            end=right_span[1],
            name=name,
            feature_type="TSD",
            strand="+",
            source=REPEATVIEWER_INTERNAL_SOURCE,
        ),
    ]


def _tsd_motif_passes_basic_quality(
    tsd_seq: str,
    *,
    k_min: int,
    complexity_threshold: float = 0.8,
) -> bool:
    """Per-copy quality used when keeping boundary-flush hits outside diversity."""
    motif = str(tsd_seq or "").upper()
    if len(motif) < max(1, int(k_min)):
        return False
    counts = collections.Counter(motif)
    _, most_freq = counts.most_common(1)[0]
    if (most_freq / len(motif)) > float(complexity_threshold):
        return False
    if len(motif) >= 5 and len(set(motif.lower())) == 1:
        return False
    return True


def _collect_boundary_flush_by_id(
    hit_maps: list[dict],
    *,
    left_boundary: int,
    right_boundary: int,
    k_min: int,
) -> dict[str, dict[str, object]]:
    """Best per-sequence hit whose TE edges abut defined boundaries (± tol).

    Used by the TEtrimmer annotation path in addition to top-cluster picks.
    ``hit_maps`` are ``{(te_start, te_end): [(seq_id, tsd_seq), ...]}`` from
    raw identify and/or post-filter maps.
    """
    best: dict[str, dict[str, object]] = {}
    left_1 = int(left_boundary)
    right_1 = int(right_boundary)
    if left_1 <= 0 or right_1 <= 0:
        return best
    for hit_map in hit_maps:
        if not hit_map:
            continue
        for (te_start, te_end), matches in hit_map.items():
            if not _te_edges_match_boundaries(
                int(te_start), int(te_end), left_1, right_1
            ):
                continue
            for raw_id, raw_seq in matches or ():
                seq_id = str(raw_id).split()[0]
                tsd_seq = str(raw_seq or "").upper()
                if not seq_id:
                    continue
                if not _tsd_motif_passes_basic_quality(tsd_seq, k_min=k_min):
                    continue
                candidate = {
                    "te_start": int(te_start),
                    "te_end": int(te_end),
                    "tsd_seq": tsd_seq,
                }
                best[seq_id] = _prefer_tsd_candidate(best.get(seq_id), candidate)
    return best


def _tsd_scan_workers(num_sequences: int) -> int:
    """Parallel processes for per-row TSD scans (Preferences → Max. parallel jobs)."""
    sequences = max(1, int(num_sequences))
    if sequences <= 1:
        return 1
    try:
        from repeatviewer.backends.thread_budget import (
            analysis_thread_pool,
            load_batch_workers,
        )

        requested = load_batch_workers()
        cpus = analysis_thread_pool()
    except Exception:
        requested = 5
        cpus = os.cpu_count() or 4
    return max(1, min(int(requested), sequences, max(1, int(cpus))))


def _scan_record_job(
    payload: tuple[str, str, int, int, int, int, int, int, int, bool],
) -> list[TsdFlankHit]:
    """Picklable worker for ``ProcessPoolExecutor`` (one MSA row)."""
    (
        seq_str,
        seq_id,
        left_boundary,
        right_boundary,
        nt_radius,
        k_min,
        k_max,
        min_distance,
        max_distance,
        whole_msa,
    ) = payload
    return _scan_record_tsd_hits(
        seq_str,
        seq_id=seq_id,
        left_boundary=left_boundary,
        right_boundary=right_boundary,
        nt_radius=nt_radius,
        k_min=k_min,
        k_max=k_max,
        min_distance=min_distance,
        max_distance=max_distance,
        whole_msa=whole_msa,
    )


def get_nt_radius_indices(sequence, center_idx, nt_radius, direction="both"):
    """MSA column window of ``nt_radius`` positions around ``center_idx``.

    Gaps count as positions. By default (``direction='both'``) the window is
    ``[center - radius, center + radius]`` inclusive (half-open end exclusive
    as ``center + radius + 1``). ``left`` / ``right`` keep one-sided windows
    for callers that need them.
    """
    seq_len = len(sequence)
    radius = max(0, int(nt_radius))
    if seq_len <= 0:
        return 0, 0
    center = max(0, min(int(center_idx), seq_len - 1))
    start_col = center
    end_col = center + 1

    if direction in ("left", "both"):
        start_col = max(0, center - radius)
    if direction in ("right", "both"):
        end_col = min(seq_len, center + radius + 1)
    elif direction == "left":
        end_col = center
    if direction == "right":
        start_col = center + 1
        end_col = min(seq_len, center + 1 + radius)
    return start_col, end_col


def get_sub_kmers(sequence, k=2):
    """Break a TSD into small overlapping k-mers (default 2-mers)."""
    return {sequence[i : i + k] for i in range(len(sequence) - k + 1)}


def calculate_jaccard_similarity(seq1, seq2):
    """Calculate similarity based on shared sub-k-mers."""
    set1 = get_sub_kmers(seq1.upper())
    set2 = get_sub_kmers(seq2.upper())
    if not set1 or not set2:
        return 0
    return len(set1.intersection(set2)) / len(set1.union(set2))


@dataclass(frozen=True)
class TsdFlankHit:
    """One TSD pair on a single MSA row (0-based inclusive MSA columns)."""

    seq_id: str
    left_start: int
    left_end: int
    right_start: int
    right_end: int
    tsd_seq: str
    te_start: int
    te_end: int


def resolve_tsd_boundaries(
    aln_len: int,
    left_boundary: int,
    right_boundary: int,
    nt_radius: int,
) -> tuple[int, int]:
    """Return MSA TE-body edge columns (0-based inclusive).

    ``left_boundary`` / ``right_boundary`` are **1-based** inclusive UI columns
    (``0`` = unset). Unset both → midpoint (whole-MSA half split).
    """
    del nt_radius
    if aln_len <= 1:
        return 0, max(0, aln_len - 1)
    left_1 = int(left_boundary)
    right_1 = int(right_boundary)
    if left_1 <= 0 and right_1 <= 0:
        mid = max(1, (aln_len - 1) // 2)
        return mid, mid
    l = max(0, min(left_1 - 1, aln_len - 1)) if left_1 > 0 else 0
    r = max(0, min(right_1 - 1, aln_len - 1)) if right_1 > 0 else aln_len - 1
    if r < l:
        l, r = r, l
    return l, r


TSD_BOUNDARY_DEFINITION_HELP = (
    "Sets the TE left and right edges; TSDs are searched within the search "
    "range (± columns) outside each edge.<br><br>"
    "<b>Auto detect</b><br>Uses Boundary prediction (Annotation tab). The "
    "first available boundary is used, in this order:<br>"
    "1. LTR or TIR<br>"
    "2. Genome Cov<br>"
    "3. MSA identity<br>"
    "Results already shown in Boundary prediction (or saved by TEtrimmer as "
    "*.boundaries_RV.json) are reused. In an open MSA without them, Boundary "
    "prediction runs first, then the TSD search starts. If no boundary is "
    "found, the TSD search splits the MSA at its middle column.<br><br>"
    "<b>User define</b><br>A window opens at run time; type the left and right "
    "columns or drag across the MSA to fill them."
)


def auto_tsd_span(prediction) -> tuple[int, int] | None:
    """Auto detect edge span (0-based, exclusive end): LTR/TIR, Genome Cov, MSA identity.

    A TEtrimmer sidecar that only stores ``final`` falls back to it.
    """
    if prediction is None:
        return None
    for span in (
        prediction.terminal_repeat,
        prediction.blast_coverage,
        prediction.msa_identity,
        prediction.final,
    ):
        if span is not None and int(span[1]) > int(span[0]):
            return int(span[0]), int(span[1])
    return None


def boundaries_from_prediction_1based(prediction) -> tuple[int, int] | None:
    """Auto detect edges as 1-based inclusive columns."""
    span = auto_tsd_span(prediction)
    if span is None:
        return None
    return span[0] + 1, span[1]  # exclusive end → 1-based last inclusive


def boundaries_from_sidecar_1based(
    file_path: str | Path,
    aln_len: int,
) -> tuple[int, int] | None:
    """Auto detect edges from ``*.boundaries_RV.json`` (1-based inclusive)."""
    from repeatviewer.backends.boundary_predict import load_boundary_prediction_json

    prediction = load_boundary_prediction_json(file_path, msa_width=aln_len)
    if prediction is None:
        return None
    span = auto_tsd_span(prediction)
    if span is None:
        return None
    start_0, end_excl = int(span[0]), int(span[1])
    if end_excl <= start_0:
        return None
    left_1 = start_0 + 1
    right_1 = end_excl  # exclusive end → 1-based last inclusive
    if right_1 < left_1:
        return None
    return left_1, right_1


def boundaries_from_sequence_extent_1based(seq: str) -> tuple[int, int] | None:
    """First/last non-gap columns of ``seq`` as 1-based inclusive edges."""
    first: int | None = None
    last: int | None = None
    for index, char in enumerate(seq):
        if char in "-.":
            continue
        if first is None:
            first = index
        last = index
    if first is None or last is None:
        return None
    return first + 1, last + 1


def preferred_tsd_boundaries_1based(
    file_path: str | Path,
    aln_len: int,
    *,
    boundary_source: str,
    saved_left: int = 0,
    saved_right: int = 0,
    reference_extent_1based: tuple[int, int] | None = None,
    reference_sequence: str | None = None,
    prediction=None,
) -> tuple[int, int]:
    """Prefill left/right (1-based) for the run-time TSD boundary dialog.

    Auto detect uses ``prediction`` (Boundary prediction in the open MSA)
    first, then ``*.boundaries_RV.json``.

    ``boundary_source`` is one of ``repeatviewer`` / ``reference`` / ``user``.
    Falls back to saved values, then ``0, 0`` (whole-MSA midpoint search).
    """
    from repeatviewer.job_options import (
        TSD_BOUNDARY_REFERENCE,
        TSD_BOUNDARY_REPEATVIEWER,
        TSD_BOUNDARY_USER,
    )

    source = str(boundary_source or TSD_BOUNDARY_REPEATVIEWER).strip().lower()
    saved = (max(0, int(saved_left)), max(0, int(saved_right)))

    if source == TSD_BOUNDARY_USER:
        return saved

    if source == TSD_BOUNDARY_REFERENCE:
        if reference_extent_1based is not None:
            left, right = reference_extent_1based
            return max(0, int(left)), max(0, int(right))
        if reference_sequence:
            extent = boundaries_from_sequence_extent_1based(reference_sequence)
            if extent is not None:
                return extent
        return saved

    # Default / Auto detect: Boundary prediction, then *.boundaries_RV.json
    if source != TSD_BOUNDARY_REPEATVIEWER:
        source = TSD_BOUNDARY_REPEATVIEWER
    del source
    auto = boundaries_from_prediction_1based(prediction)
    if auto is not None:
        return auto
    rv = boundaries_from_sidecar_1based(file_path, aln_len)
    if rv is not None:
        return rv
    return saved


def resolve_tsd_search_scope(
    file_path: str | Path,
    aln_len: int,
    left_boundary: int,
    right_boundary: int,
    nt_radius: int,
) -> tuple[int, int, int, str]:
    """Resolve TE edges (0-based) + search radius.

    Returns ``(left, right, radius, mode)`` with ``mode`` ``explicit`` when
    1-based boundaries are set, else ``whole_msa`` (midpoint half split).
    """
    del file_path
    radius = max(1, int(nt_radius))
    left_1 = int(left_boundary)
    right_1 = int(right_boundary)
    if left_1 > 0 or right_1 > 0:
        edges = resolve_tsd_boundaries(aln_len, left_1, right_1, radius)
        return edges[0], edges[1], radius, "explicit"
    edges = resolve_tsd_boundaries(aln_len, 0, 0, radius)
    return edges[0], edges[1], max(radius, aln_len), "whole_msa"


def _msa_midpoint_halves(seq_str: str) -> tuple[list[int], list[int]]:
    """Split one MSA row at the middle column into left / right nt columns."""
    aln_len = len(seq_str)
    if aln_len <= 1:
        return [], []
    mid = max(1, (aln_len - 1) // 2)
    left = [c for c in range(0, mid) if seq_str[c] != "-"]
    right = [c for c in range(mid, aln_len) if seq_str[c] != "-"]
    return left, right


def _tsd_nucleotide_span(
    sequence: str,
    anchor: int,
    length: int,
    side: str,
) -> tuple[int, int, str] | None:
    """1-based inclusive MSA columns of ``length`` nucleotides beside ``anchor``.

    ``side='left'`` walks bases immediately before ``anchor``.
    ``side='right'`` walks bases immediately after it.
    """
    if length <= 0:
        return None
    columns: list[int] = []
    if side == "left":
        index = anchor - 1
        while index >= 0 and len(columns) < length:
            if sequence[index] != "-":
                columns.append(index)
            index -= 1
        columns.reverse()
    else:
        index = anchor + 1
        while index < len(sequence) and len(columns) < length:
            if sequence[index] != "-":
                columns.append(index)
            index += 1
    if len(columns) != length:
        return None
    gathered = "".join(sequence[col] for col in columns).upper()
    return columns[0] + 1, columns[-1] + 1, gathered


def _scan_record_tsd_hits(
    seq_str: str,
    *,
    seq_id: str,
    left_boundary: int = 0,
    right_boundary: int = 0,
    nt_radius: int = 100,
    k_min: int,
    k_max: int = 20,
    min_distance: int = 0,
    max_distance: int = 0,
    whole_msa: bool = True,
) -> list[TsdFlankHit]:
    """Per-sequence TSD scan (no cross-sequence filtering or clustering).

    When ``whole_msa`` is False, ``left_boundary`` / ``right_boundary`` are
    0-based TE edges and ``nt_radius`` MSA columns on both sides of each
    edge are searched (gaps count). When True, split at the MSA midpoint.
    """
    k_min = max(2, int(k_min))
    k_max = max(k_min, int(k_max))
    min_distance = max(0, int(min_distance))
    max_distance = max(0, int(max_distance))
    all_nt_pos = [i for i, char in enumerate(seq_str) if char != "-"]
    if len(all_nt_pos) < 2 * k_min:
        return []
    nt_to_msa = {i: col for i, col in enumerate(all_nt_pos)}
    msa_to_nt = {col: i for i, col in enumerate(all_nt_pos)}

    if whole_msa:
        l_nt_pos_map, r_nt_pos_map = _msa_midpoint_halves(seq_str)
    else:
        # Search range around each edge: both sides of the left boundary and
        # both sides of the right boundary (MSA columns; gaps count).
        l_col_start, l_col_end = get_nt_radius_indices(
            seq_str, left_boundary, nt_radius, direction="both"
        )
        r_col_start, r_col_end = get_nt_radius_indices(
            seq_str, right_boundary, nt_radius, direction="both"
        )
        l_nt_pos_map = [c for c in range(l_col_start, l_col_end) if seq_str[c] != "-"]
        r_nt_pos_map = [c for c in range(r_col_start, r_col_end) if seq_str[c] != "-"]
    if len(l_nt_pos_map) < k_min or len(r_nt_pos_map) < k_min:
        return []

    k = k_min
    left_kmers_seeds: dict[str, list[int]] = defaultdict(list)
    for i in range(len(l_nt_pos_map) - k + 1):
        kmer_cols = l_nt_pos_map[i : i + k]
        kmer = "".join(seq_str[c] for c in kmer_cols)
        left_kmers_seeds[kmer].append(i)

    record_unique_hits: set[tuple[int, int, int, int]] = set()
    # te_coords -> (flank cols, tsd_seq)
    longest: dict[tuple[int, int], tuple[tuple[int, int, int, int], str]] = {}

    for j in range(len(r_nt_pos_map) - k + 1):
        kmer_cols_r = r_nt_pos_map[j : j + k]
        kmer = "".join(seq_str[c] for c in kmer_cols_r)
        if kmer not in left_kmers_seeds:
            continue

        for seed_idx_l in left_kmers_seeds[kmer]:
            # Extend maximally (TEtrimmer): no k_max cap during extension.
            l_shift = 0
            while seed_idx_l - l_shift - 1 >= 0 and j - l_shift - 1 >= 0:
                char_l = seq_str[l_nt_pos_map[seed_idx_l - l_shift - 1]]
                char_r = seq_str[r_nt_pos_map[j - l_shift - 1]]
                if char_l == char_r:
                    l_shift += 1
                else:
                    break

            r_shift = 0
            while (
                seed_idx_l + k + r_shift < len(l_nt_pos_map)
                and j + k + r_shift < len(r_nt_pos_map)
            ):
                char_l = seq_str[l_nt_pos_map[seed_idx_l + k + r_shift]]
                char_r = seq_str[r_nt_pos_map[j + k + r_shift]]
                if char_l == char_r:
                    r_shift += 1
                else:
                    break

            l_indices = l_nt_pos_map[seed_idx_l - l_shift : seed_idx_l + k + r_shift]
            r_indices = r_nt_pos_map[j - l_shift : j + k + r_shift]
            if not l_indices or not r_indices:
                continue
            hit_coords = (l_indices[0], l_indices[-1], r_indices[0], r_indices[-1])
            if hit_coords in record_unique_hits:
                continue
            record_unique_hits.add(hit_coords)

            last_nt_l_col = l_indices[-1]
            nt_idx_l = msa_to_nt[last_nt_l_col]
            if nt_idx_l + 1 not in nt_to_msa:
                continue
            te_start_col = nt_to_msa[nt_idx_l + 1]
            first_nt_r_col = r_indices[0]
            nt_idx_r = msa_to_nt[first_nt_r_col]
            if nt_idx_r - 1 not in nt_to_msa:
                continue
            te_end_col = nt_to_msa[nt_idx_r - 1]
            # Optional RV distance bounds (not in TEtrimmer).
            te_nt_span = max(0, nt_idx_r - nt_idx_l - 1)
            if te_nt_span < min_distance:
                continue
            if max_distance > 0 and te_nt_span > max_distance:
                continue
            final_kmer_seq = "".join(seq_str[c] for c in l_indices)
            # UI Max TSD length (TEtrimmer signature has k_max but does not apply it).
            if len(final_kmer_seq) > k_max:
                continue
            # Homopolymer motif filter (TEtrimmer).
            if len(final_kmer_seq) >= 5 and len(set(final_kmer_seq.lower())) == 1:
                continue

            # Longest motif per (te_start, te_end) — TEtrimmer GUI logic.
            te_key = (te_start_col, te_end_col)
            prev = longest.get(te_key)
            if prev is None or len(final_kmer_seq) > len(prev[1]):
                longest[te_key] = (hit_coords, final_kmer_seq)

    hits: list[TsdFlankHit] = []
    for (te_start, te_end), (flanks, tsd_seq) in longest.items():
        l0, l1, r0, r1 = flanks
        hits.append(
            TsdFlankHit(
                seq_id=seq_id,
                left_start=l0,
                left_end=l1,
                right_start=r0,
                right_end=r1,
                tsd_seq=tsd_seq,
                te_start=te_start,
                te_end=te_end,
            )
        )
    return hits


def scan_tsds_per_sequence(
    file_path: str | Path,
    left_boundary: int,
    right_boundary: int,
    *,
    nt_radius: int = 100,
    k_min: int = 4,
    k_max: int = 20,
    min_distance: int = 0,
    max_distance: int = 0,
    record_ids: list[str] | None = None,
    num_workers: int | None = None,
    progress: TsdProgressCallback | None = None,
) -> list[TsdFlankHit]:
    """Scan each MSA row for TSDs; no cross-sequence filter or clustering.

    With 1-based ``left_boundary`` / ``right_boundary`` set, search
    ``nt_radius`` MSA columns outside each TE edge. Otherwise midpoint split.

    Rows are scanned in parallel processes (``num_workers`` or Preferences →
    Max. parallel jobs). Filter/cluster steps stay single-process afterward.
    ``progress(done, total, label)`` is called after each row finishes.
    """
    alignment = AlignIO.read(str(file_path), "fasta")
    if not alignment:
        return []
    aln_len = len(alignment[0])
    left_0, right_0, radius, mode = resolve_tsd_search_scope(
        file_path, aln_len, left_boundary, right_boundary, nt_radius
    )
    whole_msa = mode == "whole_msa"
    jobs: list[tuple[str, str, int, int, int, int, int, int, int, bool]] = []
    for index, record in enumerate(alignment):
        if record_ids is not None and index < len(record_ids):
            seq_id = record_ids[index]
        else:
            seq_id = str(record.id)
        jobs.append(
            (
                str(record.seq),
                seq_id,
                int(left_0),
                int(right_0),
                int(radius),
                int(k_min),
                int(k_max),
                int(min_distance),
                int(max_distance),
                bool(whole_msa),
            )
        )
    total = len(jobs)
    workers = (
        max(1, int(num_workers))
        if num_workers is not None
        else _tsd_scan_workers(total)
    )

    def _emit(done: int, label: str) -> None:
        if progress is None or total <= 0:
            return
        progress(int(done), total, label)

    def _run_sequential() -> list[TsdFlankHit]:
        out: list[TsdFlankHit] = []
        for index, job in enumerate(jobs, start=1):
            out.extend(_scan_record_job(job))
            _emit(index, str(job[1]))
        return out

    # Process spawn overhead dominates tiny MSAs / test fixtures.
    if workers <= 1 or total < 4:
        return _run_sequential()

    hits: list[TsdFlankHit] = []
    try:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = {
                pool.submit(_scan_record_job, job): str(job[1]) for job in jobs
            }
            done = 0
            for fut in as_completed(futures):
                hits.extend(fut.result())
                done += 1
                _emit(done, futures[fut])
    except (PermissionError, OSError):
        # Restricted environments (some sandboxes) cannot spawn processes.
        return _run_sequential()
    return hits


def detect_tsd_features_by_id(
    file_path: str | Path,
    *,
    record_ids: list[str] | None = None,
    left_boundary: int = 0,
    right_boundary: int = 0,
    nt_radius: int = 100,
    k_min: int = 4,
    k_max: int = 20,
    min_distance: int = 0,
    max_distance: int = 0,
    max_sim_threshold: float = 0.3,
    nt_threshold: float = 0.95,
    cluster_eps: float = 1.0,
    cluster_report_path: str | Path | None = None,
    progress: TsdProgressCallback | None = None,
) -> dict[str, list[SequenceFeature]]:
    """TSD flanks as annotation features (1-based MSA columns).

    TEtrimmer ``TSD_identifier`` / ``detect_tsd_features`` pipeline:
    flank search -> process_tsd_hits (diversity / complexity) ->
    cluster_tsds_ml, then annotate:

    1. Top (largest) cluster — one best pair per sequence (TEtrimmer default).
    2. Boundary-flush pairs (when left/right boundaries are set): among
       TEtrimmer-filtered hits (diversity / complexity already applied), keep
       locations where the next nucleotide after the left TSD and the previous
       nucleotide before the right TSD sit within
       ``TSD_BOUNDARY_MATCH_TOLERANCE_BP`` of the defined start/end — even when
       that location is not in cluster 1.

    Boundaries are 1-based inclusive TE edges (``0`` = unset -> midpoint
    whole-MSA search). With edges set, ``nt_radius`` MSA columns on both sides
    of each edge are scanned. If ``cluster_report_path`` is set, write the
    full clustered text report (Tools -> TSD search).
    """
    alignment = AlignIO.read(str(file_path), "fasta")
    if not alignment:
        return {}
    sequences = {
        (
            record_ids[index]
            if record_ids is not None and index < len(record_ids)
            else str(record.id)
        ): str(record.seq)
        for index, record in enumerate(alignment)
    }
    n_rows = max(1, len(sequences))

    location_hits = identify_tsds_with_boundary_logic(
        file_path,
        left_boundary,
        right_boundary,
        nt_radius=nt_radius,
        k_min=k_min,
        k_max=k_max,
        min_distance=min_distance,
        max_distance=max_distance,
        record_ids=record_ids,
        progress=progress,
    )

    if progress is not None:
        progress(n_rows, n_rows, "TSD filter")
    filtered = process_tsd_hits(
        location_hits,
        max_sim_threshold=max_sim_threshold,
        nt_threshold=nt_threshold,
        min_length=k_min,
    )
    if progress is not None:
        progress(n_rows, n_rows, "TSD cluster")
    clusters = cluster_tsds_ml(filtered, eps=cluster_eps)
    if cluster_report_path is not None:
        write_tsd_cluster_report(clusters, cluster_report_path)

    # TEtrimmer selection: top cluster + boundary-flush (defined start/end).
    selected: dict[tuple[str, int, int, str], dict[str, object]] = {}

    def _remember(seq_id: str, te_start: int, te_end: int, tsd_seq: str) -> None:
        motif = str(tsd_seq or "").upper()
        if not seq_id or len(motif) < k_min:
            return
        key = (seq_id, int(te_start), int(te_end), motif)
        candidate = {
            "te_start": int(te_start),
            "te_end": int(te_end),
            "tsd_seq": motif,
        }
        prev = selected.get(key)
        selected[key] = _prefer_tsd_candidate(prev, candidate)

    if clusters:
        # Best hit per sequence from the largest cluster (TEtrimmer default).
        best_top: dict[str, dict[str, object]] = {}
        for item in clusters[0].get("hits") or []:
            seq_id = str(item.get("id", "")).split()[0]
            tsd_seq = str(item.get("seq") or "").upper()
            coords = item.get("coords")
            if not seq_id or not tsd_seq or not coords or len(coords) != 2:
                continue
            if len(tsd_seq) < k_min:
                continue
            candidate = {
                "te_start": int(coords[0]),
                "te_end": int(coords[1]),
                "tsd_seq": tsd_seq,
            }
            best_top[seq_id] = _prefer_tsd_candidate(best_top.get(seq_id), candidate)
        for seq_id, hit in best_top.items():
            _remember(seq_id, int(hit["te_start"]), int(hit["te_end"]), str(hit["tsd_seq"]))

    # Boundary-flush on TEtrimmer-filtered findings (not raw identify): next nt
    # after left TSD ≈ defined start, previous nt before right TSD ≈ defined
    # end (± tolerance), even when that location is outside cluster 1.
    for seq_id, hit in _collect_boundary_flush_by_id(
        [filtered],
        left_boundary=left_boundary,
        right_boundary=right_boundary,
        k_min=k_min,
    ).items():
        _remember(seq_id, int(hit["te_start"]), int(hit["te_end"]), str(hit["tsd_seq"]))

    if not selected:
        return {}

    by_id: dict[str, list[SequenceFeature]] = {}
    seen_coords: dict[str, set[tuple[int, int, int, int]]] = {}
    for (seq_id, _te_s, _te_e, _motif), hit in selected.items():
        sequence = sequences.get(seq_id)
        if sequence is None:
            continue
        feats = _tsd_features_from_hit(
            sequence,
            te_start=int(hit["te_start"]),
            te_end=int(hit["te_end"]),
            tsd_seq=str(hit["tsd_seq"]),
            k_min=k_min,
        )
        if not feats:
            continue
        coord_key = (
            int(feats[0].start),
            int(feats[0].end),
            int(feats[1].start),
            int(feats[1].end),
        )
        used = seen_coords.setdefault(seq_id, set())
        if coord_key in used:
            continue
        used.add(coord_key)
        by_id.setdefault(seq_id, []).extend(feats)
    return by_id


def identify_tsds_with_boundary_logic(
    file_path,
    left_boundary,
    right_boundary,
    nt_radius=100,
    k_min=4,
    k_max=20,
    min_distance=0,
    max_distance=0,
    record_ids=None,
    progress: TsdProgressCallback | None = None,
):
    """Legacy aggregator used by Tools → TSD search (keeps cluster pipeline)."""
    location_hits: dict[tuple[int, int], list[tuple[str, str]]] = defaultdict(list)
    for hit in scan_tsds_per_sequence(
        file_path,
        left_boundary,
        right_boundary,
        nt_radius=nt_radius,
        k_min=k_min,
        k_max=k_max,
        min_distance=min_distance,
        max_distance=max_distance,
        record_ids=record_ids,
        progress=progress,
    ):
        location_hits[(hit.te_start, hit.te_end)].append((hit.seq_id, hit.tsd_seq))
    return location_hits


def process_tsd_hits(
    hits,
    max_sim_threshold=0.35,
    nt_threshold=0.9,
    complexity_threshold=0.8,
    min_length: int = 0,
):
    """Filter TSD location hits (TEtrimmer diversity / complexity / poly-N trim).

    ``min_length`` drops motifs shortened below the caller's Min TSD length
    (poly-N edge trimming must not produce ``TSD_AA`` when min is 4).
    """
    final_processed_hits = defaultdict(list)
    min_length = max(0, int(min_length))

    for (te_start, te_end), matches in hits.items():
        if not matches:
            continue

        kmers = [m[1].upper() for m in matches]
        hit_n = len(kmers)

        if hit_n >= 2:
            # Diversity filter (TEtrimmer): true insertion-site TSDs differ
            # across copies. Identical / near-identical motifs at one MSA
            # location (avg Jaccard > threshold) are usually conserved TE
            # terminal sequence, not genomic TSDs — drop the location.
            sim_scores = [
                calculate_jaccard_similarity(s1, s2) for s1, s2 in combinations(kmers, 2)
            ]
            avg_similarity = sum(sim_scores) / len(sim_scores)
            if avg_similarity > max_sim_threshold:
                continue

        cleaned_matches = []
        for rec_id, seq in matches:
            seq_upper = seq.upper()
            seq_len = len(seq_upper)
            if seq_len == 0:
                continue
            if min_length and seq_len < min_length:
                continue

            counts = collections.Counter(seq_upper)
            _, most_freq_count = counts.most_common(1)[0]
            if (most_freq_count / seq_len) > complexity_threshold:
                continue

            cleaned_matches.append((rec_id, seq))

        hit_n = len(cleaned_matches)
        if hit_n == 0:
            continue

        first_letters = [m[1][0].upper() for m in cleaned_matches if len(m[1]) > 0]
        last_letters = [m[1][-1].upper() for m in cleaned_matches if len(m[1]) > 0]

        shift_start = 0
        shift_end = 0

        # Poly-N edge trim (TEtrimmer GUI): when ≥4 copies remain.
        if hit_n >= 4:
            if first_letters:
                counts_f = collections.Counter(first_letters)
                _, freq_f = counts_f.most_common(1)[0]
                if (freq_f / hit_n) >= nt_threshold:
                    shift_end = 1

            if last_letters:
                counts_l = collections.Counter(last_letters)
                _, freq_l = counts_l.most_common(1)[0]
                if (freq_l / hit_n) >= nt_threshold:
                    shift_start = -1

            new_loc = (te_start + shift_start, te_end + shift_end)
            final_matches = []
            for rec_id, seq in cleaned_matches:
                t_seq = seq
                if shift_start != 0 and len(t_seq) > 1:
                    t_seq = t_seq[1:]
                if shift_end != 0 and len(t_seq) > 1:
                    t_seq = t_seq[:-1]
                # RV Min TSD length: drop motifs shortened below the UI minimum.
                if min_length and len(t_seq) < min_length:
                    continue
                final_matches.append((rec_id, t_seq))

            if final_matches:
                final_processed_hits[new_loc].extend(final_matches)
        else:
            kept = [
                (rec_id, seq)
                for rec_id, seq in cleaned_matches
                if not min_length or len(seq) >= min_length
            ]
            if kept:
                final_processed_hits[te_start, te_end].extend(kept)

    return final_processed_hits


def cluster_tsds_ml(filtered_hits, eps=1):
    if not filtered_hits:
        return []

    locs = list(filtered_hits.keys())
    data = np.array(locs)

    db = DBSCAN(eps=eps, min_samples=1, metric="manhattan").fit(data)

    clusters = defaultdict(list)
    for idx, label in enumerate(db.labels_):
        clusters[label].append(locs[idx])

    results = []
    for label, member_locs in clusters.items():
        all_matches_with_coords = []
        for loc in member_locs:
            for match in filtered_hits[loc]:
                all_matches_with_coords.append(
                    {
                        "id": match[0],
                        "seq": match[1],
                        "coords": loc,
                    }
                )

        results.append(
            {
                "cluster_id": label,
                "all_locations": member_locs,
                "hits": all_matches_with_coords,
                "unique_count": len({m["id"] for m in all_matches_with_coords}),
            }
        )

    results.sort(key=lambda x: x["unique_count"], reverse=True)
    return results


def write_tsd_cluster_report(
    clusters: list[dict[str, object]],
    output_file: str | Path,
) -> Path:
    """Write DBSCAN TSD clusters to the Tools-style ``.tsd_search.txt`` report."""
    output_path = Path(output_file).expanduser().resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for index, group in enumerate(clusters, start=1):
            handle.write(f"### CLUSTER {index} ###\n")
            handle.write(
                "Unique Sequences: "
                f"{group['unique_count']} | Locations in Cluster: "
                f"{len(group['all_locations'])}\n"
            )
            handle.write("All Elements:\n")
            for item in group.get("hits") or []:
                start, end = item["coords"]
                handle.write(
                    f" TE_Start: {start:<6} | TE_End: {end:<6} | "
                    f"ID: {item['id']:<25} |  TSD: {item['seq']}\n"
                )
            handle.write("-" * 90 + "\n\n")
    return output_path


def tsd_cluster_report_path(msa_path: str | Path) -> Path:
    """``{fasta_name}.tsd_search.txt`` beside the MSA (Tools / Feature detection)."""
    path = Path(msa_path).expanduser().resolve()
    return path.parent / f"{path.name}.tsd_search.txt"


def TSD_identifier(
    left_n: int,
    right_n: int,
    extend_n: int,
    input_file: str | Path,
    output_file: str | Path,
    *,
    k_min: int = 4,
    k_max: int = 20,
    min_distance: int = 0,
    max_distance: int = 0,
    max_sim_threshold: float = 0.3,
    nt_threshold: float = 0.95,
    cluster_eps: float = 1.0,
) -> list[list[int]]:
    """Search TSDs near MSA boundaries and write clustered results to a text file."""
    input_path = Path(input_file).expanduser().resolve()

    location_hits = identify_tsds_with_boundary_logic(
        str(input_path),
        left_n,
        right_n,
        nt_radius=extend_n,
        k_min=k_min,
        k_max=k_max,
        min_distance=min_distance,
        max_distance=max_distance,
    )
    filtered_location_hits = process_tsd_hits(
        location_hits,
        max_sim_threshold=max_sim_threshold,
        nt_threshold=nt_threshold,
        min_length=k_min,
    )
    final_clusters = cluster_tsds_ml(filtered_location_hits, eps=cluster_eps)
    write_tsd_cluster_report(final_clusters, output_file)

    left_pos: list[int] = []
    right_pos: list[int] = []
    for group in final_clusters:
        for l_coord, r_coord in group["all_locations"]:
            left_pos.append(l_coord)
            right_pos.append(r_coord)

    return [left_pos, right_pos]
