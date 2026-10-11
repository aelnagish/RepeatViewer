"""TSD (target site duplication) search in aligned FASTA MSAs (from TEtrimmer TSD_search.py)."""

from __future__ import annotations

import collections
from collections import defaultdict
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


def get_nt_radius_indices(sequence, center_idx, nt_radius, direction="both"):
    """Find MSA column indices containing nt_radius nucleotides, skipping gaps."""
    seq_len = len(sequence)
    start_col = center_idx
    end_col = center_idx

    if direction in ("left", "both"):
        count = 0
        while start_col > 0 and count < nt_radius:
            start_col -= 1
            if sequence[start_col] != "-":
                count += 1
    if direction in ("right", "both"):
        count = 0
        while end_col < seq_len and count < nt_radius:
            if sequence[end_col] != "-":
                count += 1
            end_col += 1
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
    """Return MSA TE-body columns. ``0, 0`` means auto (leave flank windows)."""
    if aln_len <= 1:
        return 0, max(0, aln_len - 1)
    left = int(left_boundary)
    right = int(right_boundary)
    if left > 0 or right > 0:
        l = max(0, min(left, aln_len - 1))
        r = right if right > 0 else aln_len - 1
        r = max(0, min(r, aln_len - 1))
        if r < l:
            l, r = r, l
        return l, r
    pad = max(1, int(nt_radius))
    if aln_len <= 2 * pad + 8:
        pad = max(1, aln_len // 10)
    return pad, max(pad, aln_len - pad - 1)


def _scan_record_tsd_hits(
    seq_str: str,
    *,
    seq_id: str,
    left_boundary: int,
    right_boundary: int,
    nt_radius: int,
    k_min: int,
) -> list[TsdFlankHit]:
    """Per-sequence TSD scan (no cross-sequence filtering or clustering)."""
    all_nt_pos = [i for i, char in enumerate(seq_str) if char != "-"]
    if len(all_nt_pos) < 2 * k_min:
        return []
    nt_to_msa = {i: col for i, col in enumerate(all_nt_pos)}
    msa_to_nt = {col: i for i, col in enumerate(all_nt_pos)}

    l_col_start, l_col_end = get_nt_radius_indices(
        seq_str, left_boundary, nt_radius, direction="left"
    )
    r_col_start, r_col_end = get_nt_radius_indices(
        seq_str, right_boundary, nt_radius, direction="right"
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
            final_kmer_seq = "".join(seq_str[c] for c in l_indices)
            if len(final_kmer_seq) >= 5 and len(set(final_kmer_seq.lower())) == 1:
                continue

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
    record_ids: list[str] | None = None,
) -> list[TsdFlankHit]:
    """Scan each MSA row for TSDs; no cross-sequence filter or clustering."""
    alignment = AlignIO.read(str(file_path), "fasta")
    if not alignment:
        return []
    aln_len = len(alignment[0])
    left, right = resolve_tsd_boundaries(
        aln_len, left_boundary, right_boundary, nt_radius
    )
    hits: list[TsdFlankHit] = []
    for index, record in enumerate(alignment):
        if record_ids is not None and index < len(record_ids):
            seq_id = record_ids[index]
        else:
            seq_id = str(record.id)
        hits.extend(
            _scan_record_tsd_hits(
                str(record.seq),
                seq_id=seq_id,
                left_boundary=left,
                right_boundary=right,
                nt_radius=nt_radius,
                k_min=k_min,
            )
        )
    return hits


def detect_tsd_features_by_id(
    file_path: str | Path,
    *,
    record_ids: list[str] | None = None,
    left_boundary: int = 0,
    right_boundary: int = 0,
    nt_radius: int = 100,
    k_min: int = 4,
) -> dict[str, list[SequenceFeature]]:
    """Per-sequence TSD flanks as annotation features (1-based MSA columns).

    No TSD cluster calculation — each sequence's hits are annotated independently.
    """
    hits = scan_tsds_per_sequence(
        file_path,
        left_boundary,
        right_boundary,
        nt_radius=nt_radius,
        k_min=k_min,
        record_ids=record_ids,
    )
    by_id: dict[str, list[SequenceFeature]] = {}
    for hit in hits:
        name = f"TSD_{hit.tsd_seq.upper()}"
        left = SequenceFeature(
            start=hit.left_start + 1,
            end=hit.left_end + 1,
            name=name,
            feature_type="TSD",
            strand="+",
            source=REPEATVIEWER_INTERNAL_SOURCE,
        )
        right = SequenceFeature(
            start=hit.right_start + 1,
            end=hit.right_end + 1,
            name=name,
            feature_type="TSD",
            strand="+",
            source=REPEATVIEWER_INTERNAL_SOURCE,
        )
        by_id.setdefault(hit.seq_id, []).extend([left, right])
    return by_id


def identify_tsds_with_boundary_logic(
    file_path,
    left_boundary,
    right_boundary,
    nt_radius=100,
    k_min=4,
    k_max=7,
):
    """Legacy aggregator used by Tools → TSD search (keeps cluster pipeline)."""
    del k_max  # retained for TEtrimmer API compatibility
    location_hits: dict[tuple[int, int], list[tuple[str, str]]] = defaultdict(list)
    for hit in scan_tsds_per_sequence(
        file_path,
        left_boundary,
        right_boundary,
        nt_radius=nt_radius,
        k_min=k_min,
    ):
        location_hits[(hit.te_start, hit.te_end)].append((hit.seq_id, hit.tsd_seq))
    return location_hits


def process_tsd_hits(
    hits,
    max_sim_threshold=0.35,
    nt_threshold=0.9,
    complexity_threshold=0.8,
):
    final_processed_hits = defaultdict(list)

    for (te_start, te_end), matches in hits.items():
        if not matches:
            continue

        kmers = [m[1].upper() for m in matches]
        hit_n = len(kmers)

        if hit_n >= 2:
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
                final_matches.append((rec_id, t_seq))

            final_processed_hits[new_loc].extend(final_matches)
        else:
            final_processed_hits[te_start, te_end].extend(cleaned_matches)

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


def TSD_identifier(
    left_n: int,
    right_n: int,
    extend_n: int,
    input_file: str | Path,
    output_file: str | Path,
    *,
    max_sim_threshold: float = 0.3,
    nt_threshold: float = 0.95,
    cluster_eps: float = 1.0,
) -> list[list[int]]:
    """Search TSDs near MSA boundaries and write clustered results to a text file."""
    input_path = Path(input_file).expanduser().resolve()
    output_path = Path(output_file).expanduser().resolve()

    location_hits = identify_tsds_with_boundary_logic(
        str(input_path),
        left_n,
        right_n,
        nt_radius=extend_n,
    )
    filtered_location_hits = process_tsd_hits(
        location_hits,
        max_sim_threshold=max_sim_threshold,
        nt_threshold=nt_threshold,
    )
    final_clusters = cluster_tsds_ml(filtered_location_hits, eps=cluster_eps)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8") as handle:
        for index, group in enumerate(final_clusters, start=1):
            handle.write(f"### CLUSTER {index} ###\n")
            handle.write(
                "Unique Sequences: "
                f"{group['unique_count']} | Locations in Cluster: "
                f"{len(group['all_locations'])}\n"
            )
            handle.write("All Elements:\n")
            for item in group["hits"]:
                start, end = item["coords"]
                handle.write(
                    f" TE_Start: {start:<6} | TE_End: {end:<6} | "
                    f"ID: {item['id']:<25} |  TSD: {item['seq']}\n"
                )
            handle.write("-" * 90 + "\n\n")

    left_pos: list[int] = []
    right_pos: list[int] = []
    for group in final_clusters:
        for l_coord, r_coord in group["all_locations"]:
            left_pos.append(l_coord)
            right_pos.append(r_coord)

    return [left_pos, right_pos]
