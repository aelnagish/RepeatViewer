"""Numpy-only genome-map binning (ProcessPool-safe: no Plotly/pandas imports)."""

from __future__ import annotations

from typing import Any

import numpy as np

from repeatviewer.job_options import TE_GENOME_MAP_DENSITY_TE_LENGTH

BP_PER_MBP = 1_000_000.0


def _filter_bins_to_mbp_window(
    centers_mbp: np.ndarray,
    density: np.ndarray,
    mean_pident: np.ndarray,
    hit_counts: np.ndarray,
    te_bp: np.ndarray,
    bin_size: int,
    lo_mbp: float,
    hi_mbp: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    half_w = (bin_size / BP_PER_MBP) / 2.0
    keep = (centers_mbp + half_w >= lo_mbp) & (centers_mbp - half_w <= hi_mbp)
    return (
        centers_mbp[keep],
        density[keep],
        mean_pident[keep],
        hit_counts[keep],
        te_bp[keep],
    )


def _union_coverage_and_mean_pident(
    segments: list[tuple[int, int, float]],
) -> tuple[float, float]:
    """Union bp length and length-weighted mean identity for overlapping segments."""
    if not segments:
        return 0.0, float("nan")
    segs = sorted((int(lo), int(hi), float(p)) for lo, hi, p in segments if hi >= lo)
    if not segs:
        return 0.0, float("nan")

    boundaries = sorted({x for lo, hi, _ in segs for x in (lo, hi + 1)})
    total_len = 0
    pident_len_sum = 0.0
    for i in range(len(boundaries) - 1):
        seg_lo = boundaries[i]
        seg_hi = boundaries[i + 1] - 1
        if seg_hi < seg_lo:
            continue
        covering = [p for lo, hi, p in segs if lo <= seg_lo and hi >= seg_hi]
        if not covering:
            continue
        length = seg_hi - seg_lo + 1
        mean_p = sum(covering) / len(covering)
        total_len += length
        pident_len_sum += mean_p * length
    if total_len <= 0:
        return 0.0, float("nan")
    return float(total_len), float(pident_len_sum / total_len)


def _accumulate_hit_bins(
    starts: np.ndarray,
    ends: np.ndarray,
    pidents: np.ndarray,
    counts: np.ndarray,
    pident_sum: np.ndarray,
    bp_cover: np.ndarray,
    *,
    bin_size_i: int,
    n_bins: int,
    cat_indices: np.ndarray | None = None,
    collapse_overlaps: bool = False,
) -> None:
    """Accumulate per-bin hit statistics (optionally per TE category)."""
    multi_cat = cat_indices is not None and counts.ndim == 2
    collapsed_segments: dict[tuple[int, int], list[tuple[int, int, float]]] | None = None
    if collapse_overlaps:
        collapsed_segments = {}

    for i in range(starts.size):
        start = int(starts[i])
        end = int(ends[i])
        if end < start:
            start, end = end, start
        if multi_cat:
            ci = int(cat_indices[i])
            if ci < 0:
                continue
            row_counts = counts[ci]
            row_pident = pident_sum[ci]
            row_bp = bp_cover[ci]
        else:
            ci = 0
            row_counts = counts
            row_pident = pident_sum
            row_bp = bp_cover
        pident = float(pidents[i])
        b0 = start // bin_size_i
        b1 = min(end // bin_size_i, n_bins - 1)
        for bin_idx in range(b0, b1 + 1):
            bin_lo = bin_idx * bin_size_i
            seg_lo = max(start, bin_lo)
            seg_hi = min(end, min((bin_idx + 1) * bin_size_i - 1, end))
            overlap = seg_hi - seg_lo + 1
            if overlap <= 0:
                continue
            row_counts[bin_idx] += 1.0
            if collapse_overlaps:
                assert collapsed_segments is not None
                key = (ci, bin_idx)
                collapsed_segments.setdefault(key, []).append((seg_lo, seg_hi, pident))
            else:
                row_pident[bin_idx] += pident * overlap
                row_bp[bin_idx] += overlap

    if collapse_overlaps and collapsed_segments:
        for (ci, bin_idx), segments in collapsed_segments.items():
            union_len, mean_pident = _union_coverage_and_mean_pident(segments)
            if multi_cat:
                bp_cover[ci, bin_idx] = union_len
                pident_sum[ci, bin_idx] = (
                    mean_pident * union_len if union_len > 0 else 0.0
                )
            else:
                bp_cover[bin_idx] = union_len
                pident_sum[bin_idx] = mean_pident * union_len if union_len > 0 else 0.0


def _bin_arrays_single_track(
    starts: np.ndarray,
    ends: np.ndarray,
    pidents: np.ndarray,
    chrom_length: int,
    bin_size: int,
    *,
    density_mode: str,
    collapse_overlaps: bool = False,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    bin_size_i = max(5_000, int(bin_size))
    n_bins = max(1, int(np.ceil(chrom_length / bin_size_i)))
    counts = np.zeros(n_bins, dtype=np.float64)
    pident_sum = np.zeros(n_bins, dtype=np.float64)
    bp_cover = np.zeros(n_bins, dtype=np.float64)
    use_collapse = bool(collapse_overlaps) and density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH
    if starts.size:
        _accumulate_hit_bins(
            starts,
            ends,
            pidents,
            counts,
            pident_sum,
            bp_cover,
            bin_size_i=bin_size_i,
            n_bins=n_bins,
            collapse_overlaps=use_collapse,
        )
    centers_mbp = (np.arange(n_bins) + 0.5) * bin_size_i / BP_PER_MBP
    mean_pident = np.divide(
        pident_sum,
        bp_cover,
        out=np.full(n_bins, np.nan, dtype=np.float64),
        where=bp_cover > 0,
    )
    density = bp_cover if density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH else counts
    return centers_mbp, density, mean_pident, counts, bp_cover


def _bin_arrays_by_category(
    starts: np.ndarray,
    ends: np.ndarray,
    pidents: np.ndarray,
    cat_values: np.ndarray,
    chrom_length: int,
    bin_size: int,
    categories: list[str],
    *,
    density_mode: str,
    collapse_overlaps: bool = False,
) -> dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]]:
    bin_size_i = max(5_000, int(bin_size))
    n_bins = max(1, int(np.ceil(chrom_length / bin_size_i)))
    n_cats = len(categories)
    cat_to_idx = {str(cat): idx for idx, cat in enumerate(categories)}
    counts = np.zeros((n_cats, n_bins), dtype=np.float64)
    pident_sum = np.zeros((n_cats, n_bins), dtype=np.float64)
    bp_cover = np.zeros((n_cats, n_bins), dtype=np.float64)
    use_collapse = bool(collapse_overlaps) and density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH
    if starts.size:
        cat_indices = np.fromiter(
            (cat_to_idx.get(str(c), -1) for c in cat_values),
            dtype=np.int16,
            count=cat_values.size,
        )
        _accumulate_hit_bins(
            starts,
            ends,
            pidents,
            counts,
            pident_sum,
            bp_cover,
            bin_size_i=bin_size_i,
            n_bins=n_bins,
            cat_indices=cat_indices,
            collapse_overlaps=use_collapse,
        )
    centers_mbp = (np.arange(n_bins) + 0.5) * bin_size_i / BP_PER_MBP
    out: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]] = {}
    for idx, cat in enumerate(categories):
        mean_pident = np.divide(
            pident_sum[idx],
            bp_cover[idx],
            out=np.full(n_bins, np.nan, dtype=np.float64),
            where=bp_cover[idx] > 0,
        )
        density = bp_cover[idx] if density_mode == TE_GENOME_MAP_DENSITY_TE_LENGTH else counts[idx]
        out[cat] = (centers_mbp, density, mean_pident, counts[idx], bp_cover[idx])
    return out


def _maybe_filter_bins_region(
    binned: tuple[np.ndarray, ...],
    bin_size: int,
    region: tuple[float, float] | None,
) -> tuple[np.ndarray, ...]:
    if region is None:
        return binned
    return _filter_bins_to_mbp_window(*binned, bin_size, region[0], region[1])


def bin_chromosome_tracks_worker(payload: dict[str, Any]) -> dict[str, Any]:
    """ProcessPool worker: bin all tracks on one chromosome."""
    track_cat_bins: dict[tuple[str, str], tuple] = {}
    track_bins: dict[str, tuple] = {}
    global_max = 1.0
    chrom_length = int(payload["chrom_length"])
    bin_size = int(payload["bin_size"])
    density_mode = str(payload["density_mode"])
    collapse_overlaps = bool(payload.get("collapse_overlaps", False))
    use_category = bool(payload["use_category"])
    categories = list(payload.get("categories") or [])
    region = payload.get("region")

    for spec in payload["tracks"]:
        track_label = str(spec["track_label"])
        track_density_mode = str(spec.get("density_mode", density_mode))
        starts = np.asarray(spec["starts"], dtype=np.int64)
        ends = np.asarray(spec["ends"], dtype=np.int64)
        pidents = np.asarray(spec["pidents"], dtype=np.float64)
        use_category = bool(spec.get("use_category", payload.get("use_category")))
        if use_category:
            cats = np.asarray(spec.get("cats", []), dtype=str)
            per_cat = _bin_arrays_by_category(
                starts,
                ends,
                pidents,
                cats,
                chrom_length,
                bin_size,
                categories,
                density_mode=track_density_mode,
                collapse_overlaps=collapse_overlaps,
            )
            for cat in categories:
                binned = per_cat.get(cat)
                if binned is None:
                    continue
                binned = _maybe_filter_bins_region(binned, bin_size, region)
                track_cat_bins[(track_label, cat)] = binned
            if categories:
                n_bins = 0
                stack_totals: np.ndarray | None = None
                for cat in categories:
                    binned = track_cat_bins.get((track_label, cat))
                    if binned is None or not binned[1].size:
                        continue
                    if stack_totals is None:
                        n_bins = len(binned[1])
                        stack_totals = np.zeros(n_bins, dtype=float)
                    stack_totals += binned[1]
                if stack_totals is not None and stack_totals.size:
                    global_max = max(global_max, float(stack_totals.max()))
        else:
            binned = _bin_arrays_single_track(
                starts,
                ends,
                pidents,
                chrom_length,
                bin_size,
                density_mode=track_density_mode,
                collapse_overlaps=collapse_overlaps,
            )
            binned = _maybe_filter_bins_region(binned, bin_size, region)
            if binned[1].size and float(binned[1].max()) > global_max:
                global_max = float(binned[1].max())
            track_bins[track_label] = binned

    return {
        "track_cat_bins": track_cat_bins,
        "track_bins": track_bins,
        "global_max": global_max,
    }
