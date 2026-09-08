"""Load Bismark coverage (.cov / .cov.gz) files for methylation genome-map tracks."""

from __future__ import annotations

import gzip
import logging
import os
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from Bio import SeqIO

from repeatviewer.backends.gff_features import _seqid_matches

logger = logging.getLogger(__name__)

# Per-chromosome .cov arrays keyed by "pos" and "meth_pct".
MethylationChromCov = dict[str, np.ndarray]

_ARRAY_BIN_MIN_SITES = 500_000
_ARRAY_BIN_MAX_THREADS = 8
_COV_CACHE_VERSION = 1
_CG_BINS_CACHE_VERSION = 1
_METH_BINS_CACHE_VERSION = 1

METH_SUBTRACK_CG_TOTAL = "cg_total"
METH_SUBTRACK_CG_COVERED = "cg_covered"
METH_SUBTRACK_CG_COVERAGE_PCT = "cg_coverage_pct"
METH_SUBTRACK_METH_LEVEL = "meth_level"

# Legacy aliases (older HTML / pickles).
METH_SUBTRACK_C_TOTAL = METH_SUBTRACK_CG_TOTAL
METH_SUBTRACK_C_COUNTED = METH_SUBTRACK_CG_COVERED
METH_SUBTRACK_LEVEL = METH_SUBTRACK_METH_LEVEL

BP_PER_MBP = 1_000_000.0

_CG_BIN_CACHE: dict[tuple[str, int, str, int], np.ndarray] = {}
_GENOME_FASTA_INDEX: dict[tuple[str, int], object] = {}
_CHROM_RECORD_ID_CACHE: dict[tuple[str, int], dict[str, str]] = {}


def is_bismark_cov_file(path: Path) -> bool:
    """True for Bismark coverage files (plain or gzipped)."""
    name = path.name.lower()
    if name.endswith(".cov") or name.endswith(".cov.gz"):
        return True
    if name.endswith(".gz") and "bismark" in name and ".cov" in name:
        return True
    if name.endswith(".gz") and name.endswith("bismark.cov.gz"):
        return True
    return "bismark" in name and (name.endswith(".cov") or ".cov." in name)


def _open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open(encoding="utf-8", errors="replace")


def load_bismark_cov(path: Path) -> pd.DataFrame:
    """Read Bismark coverage: chrom, pos, meth_pct, methylated, unmethylated."""
    path = path.expanduser().resolve()
    compression = "gzip" if path.name.lower().endswith(".gz") else None
    try:
        df = pd.read_csv(
            path,
            sep=r"\s+",
            header=None,
            comment="#",
            compression=compression,
            usecols=[0, 1, 3, 4, 5],
            names=["chrom", "pos", "meth_pct", "methylated", "unmethylated"],
            dtype={
                "chrom": str,
                "pos": np.int64,
                "methylated": np.int32,
                "unmethylated": np.int32,
            },
            engine="c",
            on_bad_lines="skip",
        )
    except (OSError, ValueError, pd.errors.ParserError) as exc:
        logger.warning("Bismark cov: fast read failed for %s (%s); falling back", path, exc)
        return _load_bismark_cov_slow(path)
    if df.empty:
        return df
    df = df[~df["chrom"].astype(str).str.lower().str.startswith("track")]
    df["meth_pct"] = pd.to_numeric(df["meth_pct"], errors="coerce").fillna(0.0)
    df["coverage"] = df["methylated"] + df["unmethylated"]
    return df.reset_index(drop=True)


def _load_bismark_cov_slow(path: Path) -> pd.DataFrame:
    rows: list[tuple[str, int, float, int, int]] = []
    try:
        with _open_text(path) as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#") or line.startswith("track"):
                    continue
                parts = line.split("\t")
                if len(parts) < 6:
                    parts = line.split()
                if len(parts) < 6:
                    continue
                chrom = str(parts[0]).strip()
                try:
                    start = int(float(parts[1]))
                    meth_pct = float(parts[3])
                    methylated = int(float(parts[4]))
                    unmethylated = int(float(parts[5]))
                except (TypeError, ValueError):
                    continue
                rows.append((chrom, start, meth_pct, methylated, unmethylated))
    except OSError as exc:
        logger.warning("Bismark cov: cannot read %s: %s", path, exc)
        return pd.DataFrame()
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(
        rows,
        columns=["chrom", "pos", "meth_pct", "methylated", "unmethylated"],
    )
    df["coverage"] = df["methylated"] + df["unmethylated"]
    return df


def _build_bismark_chrom_mapping(
    unique_raw_chroms: np.ndarray | list[str],
    chrom_sizes: dict[str, int],
) -> dict[str, str]:
    genome_names = list(chrom_sizes.keys())
    mapping: dict[str, str] = {}
    for raw in unique_raw_chroms:
        raw_s = str(raw)
        if raw_s in mapping:
            continue
        for gname in genome_names:
            if _seqid_matches(gname, raw_s):
                mapping[raw_s] = gname
                break
    return mapping


def map_bismark_chromosomes(df: pd.DataFrame, chrom_sizes: dict[str, int]) -> pd.DataFrame:
    """Map Bismark seqids to genome chromosome names."""
    if df.empty or not chrom_sizes:
        return pd.DataFrame()
    mapping = _build_bismark_chrom_mapping(df["chrom"].astype(str).unique(), chrom_sizes)
    if not mapping:
        return pd.DataFrame()
    raw = df["chrom"].astype(str).to_numpy()
    mapped = np.fromiter((mapping.get(str(r)) for r in raw), dtype=object, count=len(raw))
    valid = mapped != None  # noqa: E711
    if not np.any(valid):
        return pd.DataFrame()
    out = df.loc[valid].copy()
    out["chrom"] = mapped[valid]
    max_pos = out["chrom"].map(chrom_sizes).to_numpy(dtype=np.int64, copy=False)
    pos = out["pos"].to_numpy(dtype=np.int64, copy=False)
    in_range = (pos >= 1) & (pos <= max_pos)
    out = out.loc[in_range]
    return out.reset_index(drop=True)


def partition_cov_by_chrom_arrays(
    df: pd.DataFrame,
    *,
    coverage_min: int = 1,
) -> dict[str, MethylationChromCov]:
    """Split mapped .cov rows into per-chrom numpy arrays (no per-chrom DataFrames)."""
    if df.empty:
        return {}
    chrom = df["chrom"].astype(str).to_numpy()
    pos = df["pos"].to_numpy(dtype=np.int64, copy=False)
    meth = df["meth_pct"].to_numpy(dtype=np.float64, copy=False)
    if coverage_min > 1 and "coverage" in df.columns:
        cov = df["coverage"].to_numpy(dtype=np.int32, copy=False)
        keep = cov >= int(coverage_min)
        chrom = chrom[keep]
        pos = pos[keep]
        meth = meth[keep]
    if chrom.size == 0:
        return {}
    order = np.argsort(chrom, kind="mergesort")
    chrom_s = chrom[order]
    pos_s = pos[order]
    meth_s = meth[order]
    unique, starts = np.unique(chrom_s, return_index=True)
    ends = np.append(starts[1:], chrom_s.size)
    out: dict[str, MethylationChromCov] = {}
    for i, name in enumerate(unique):
        sl = slice(int(starts[i]), int(ends[i]))
        out[str(name)] = {
            "pos": np.ascontiguousarray(pos_s[sl]),
            "meth_pct": np.ascontiguousarray(meth_s[sl]),
        }
    return out


def methylation_cov_site_count(cov: MethylationChromCov | pd.DataFrame | None) -> int:
    if cov is None:
        return 0
    if isinstance(cov, dict):
        pos = cov.get("pos")
        return int(pos.size) if isinstance(pos, np.ndarray) else 0
    return len(cov)


def methylation_cov_arrays(cov: MethylationChromCov | pd.DataFrame | None) -> MethylationChromCov:
    if cov is None:
        return {}
    if isinstance(cov, dict) and "pos" in cov:
        return cov
    if isinstance(cov, pd.DataFrame):
        return methylation_cov_payload(cov)
    return {}


def _methylation_cov_cache_path(path: Path) -> Path:
    path = path.expanduser().resolve()
    return path.parent / f"{path.name}.repeatviewer.npz"


def _methylation_bins_cache_path(path: Path, bin_size: int) -> Path:
    path = path.expanduser().resolve()
    bin_size_i = max(5_000, int(bin_size))
    return path.parent / f"{path.name}.bins_{bin_size_i}.repeatviewer.npz"


def _save_methylation_cov_cache(
    cache_path: Path,
    by_chrom: dict[str, MethylationChromCov],
    *,
    coverage_min: int,
) -> None:
    payload: dict[str, np.ndarray] = {
        "_version": np.array([_COV_CACHE_VERSION], dtype=np.int32),
        "_coverage_min": np.array([int(coverage_min)], dtype=np.int32),
        "_chrom_names": np.asarray(list(by_chrom.keys()), dtype=str),
    }
    for chrom, cov in by_chrom.items():
        safe = chrom.replace("/", "_")
        payload[f"{safe}__pos"] = cov["pos"]
        payload[f"{safe}__meth_pct"] = cov["meth_pct"]
    try:
        np.savez_compressed(cache_path, **payload)
    except OSError as exc:
        logger.warning("Bismark cov: could not write cache %s: %s", cache_path, exc)


def _load_methylation_cov_cache(
    cache_path: Path,
    *,
    coverage_min: int,
) -> dict[str, MethylationChromCov] | None:
    try:
        with np.load(cache_path, allow_pickle=False) as data:
            if int(data["_version"][0]) != _COV_CACHE_VERSION:
                return None
            cached_min = int(data["_coverage_min"][0])
            if cached_min < int(coverage_min):
                return None
            chrom_names = [str(x) for x in data["_chrom_names"].tolist()]
            out: dict[str, MethylationChromCov] = {}
            for chrom in chrom_names:
                safe = chrom.replace("/", "_")
                pos_key = f"{safe}__pos"
                meth_key = f"{safe}__meth_pct"
                if pos_key not in data or meth_key not in data:
                    return None
                pos = np.asarray(data[pos_key], dtype=np.int64)
                meth = np.asarray(data[meth_key], dtype=np.float64)
                out[chrom] = {
                    "pos": np.ascontiguousarray(pos),
                    "meth_pct": np.ascontiguousarray(meth),
                }
            return out
    except (OSError, KeyError, ValueError, TypeError) as exc:
        logger.warning("Bismark cov: could not read cache %s: %s", cache_path, exc)
        return None


def load_methylation_cov_chrom_from_cache(
    cache_path: Path,
    chrom: str,
) -> MethylationChromCov:
    """Load one chromosome slice from a RepeatViewer .cov cache (ProcessPool-safe)."""
    cache_path = cache_path.expanduser().resolve()
    try:
        with np.load(cache_path, allow_pickle=False) as data:
            if int(data["_version"][0]) != _COV_CACHE_VERSION:
                return {}
            chrom_names = [str(x) for x in data["_chrom_names"].tolist()]
            if chrom not in chrom_names:
                return {}
            safe = chrom.replace("/", "_")
            pos_key = f"{safe}__pos"
            meth_key = f"{safe}__meth_pct"
            if pos_key not in data or meth_key not in data:
                return {}
            return {
                "pos": np.asarray(data[pos_key], dtype=np.int64),
                "meth_pct": np.asarray(data[meth_key], dtype=np.float64),
            }
    except (OSError, KeyError, ValueError, TypeError) as exc:
        logger.debug("Bismark cov: chrom cache read failed %s (%s): %s", cache_path, chrom, exc)
        return {}


def _load_methylation_bins_chrom_cache(
    cache_path: Path,
    chrom: str,
    *,
    bin_size: int,
    cg_total: np.ndarray,
) -> MethylationChromBins | None:
    try:
        with np.load(cache_path, allow_pickle=False) as data:
            if int(data["_version"][0]) != _METH_BINS_CACHE_VERSION:
                return None
            if int(data["_bin_size"][0]) != max(5_000, int(bin_size)):
                return None
            names = [str(x) for x in data["_chrom_names"].tolist()]
            if chrom not in names:
                return None
            safe = chrom.replace("/", "_")
            counts_key = f"{safe}__counts"
            meth_key = f"{safe}__meth_sum"
            if counts_key not in data or meth_key not in data:
                return None
            line_counts = np.asarray(data[counts_key], dtype=np.float64)
            meth_sum = np.asarray(data[meth_key], dtype=np.float64)
            n_bins = line_counts.size
            bin_size_i = max(5_000, int(bin_size))
            centers_mbp = (np.arange(n_bins) + 0.5) * bin_size_i / BP_PER_MBP
            return MethylationChromBins(
                centers_mbp=centers_mbp,
                cg_total=_align_cg_total(cg_total, n_bins),
                line_counts=line_counts,
                meth_sum=meth_sum,
                bin_size=bin_size_i,
            )
    except (OSError, KeyError, ValueError, TypeError) as exc:
        logger.debug("Bismark cov: bins cache read failed %s (%s): %s", cache_path, chrom, exc)
        return None


def _save_methylation_bins_cache(
    cache_path: Path,
    by_chrom: dict[str, MethylationChromBins],
    *,
    bin_size: int,
    source_mtime: float,
) -> None:
    if not by_chrom:
        return
    bin_size_i = max(5_000, int(bin_size))
    payload: dict[str, np.ndarray] = {
        "_version": np.array([_METH_BINS_CACHE_VERSION], dtype=np.int32),
        "_bin_size": np.array([bin_size_i], dtype=np.int32),
        "_source_mtime": np.array([float(source_mtime)], dtype=np.float64),
        "_chrom_names": np.asarray(list(by_chrom.keys()), dtype=str),
    }
    for chrom, bins in by_chrom.items():
        safe = chrom.replace("/", "_")
        payload[f"{safe}__counts"] = np.ascontiguousarray(bins.line_counts, dtype=np.float64)
        payload[f"{safe}__meth_sum"] = np.ascontiguousarray(bins.meth_sum, dtype=np.float64)
    try:
        np.savez_compressed(cache_path, **payload)
    except OSError as exc:
        logger.warning("Bismark cov: could not write bins cache %s: %s", cache_path, exc)


def try_load_methylation_bins_cache(
    cov_path: Path,
    chromosomes: list[str],
    *,
    bin_size: int,
    cg_by_chrom: dict[str, np.ndarray],
) -> dict[str, MethylationChromBins] | None:
    """Return cached per-chromosome binned methylation when cov file and bin size match."""
    cov_path = cov_path.expanduser().resolve()
    cache_path = _methylation_bins_cache_path(cov_path, bin_size)
    if not cache_path.is_file():
        return None
    try:
        cov_mtime = cov_path.stat().st_mtime
        if cache_path.stat().st_mtime < cov_mtime:
            return None
    except OSError:
        return None
    out: dict[str, MethylationChromBins] = {}
    for chrom in chromosomes:
        bins = _load_methylation_bins_chrom_cache(
            cache_path,
            chrom,
            bin_size=bin_size,
            cg_total=cg_by_chrom.get(
                chrom,
                np.zeros(1, dtype=np.float64),
            ),
        )
        if bins is None:
            return None
        out[chrom] = bins
    logger.info(
        "Bismark cov: loaded binned cache %s (%d chromosomes)",
        cache_path.name,
        len(out),
    )
    return out


def load_methylation_cov_by_chrom(
    path: Path,
    chrom_sizes: dict[str, int],
    *,
    coverage_min: int = 1,
    use_cache: bool = True,
) -> dict[str, MethylationChromCov]:
    """Load a Bismark .cov file into per-chromosome numpy arrays."""
    path = path.expanduser().resolve()
    coverage_min = max(1, int(coverage_min))
    cache_path = _methylation_cov_cache_path(path)
    if use_cache and cache_path.is_file():
        try:
            if cache_path.stat().st_mtime >= path.stat().st_mtime:
                cached = _load_methylation_cov_cache(cache_path, coverage_min=coverage_min)
                if cached is not None:
                    logger.info("Bismark cov: loaded cache %s (%d chromosomes)", cache_path.name, len(cached))
                    return cached
        except OSError:
            pass
    t0 = time.perf_counter()
    raw = load_bismark_cov(path)
    mapped = map_bismark_chromosomes(raw, chrom_sizes)
    by_chrom = partition_cov_by_chrom_arrays(mapped, coverage_min=coverage_min)
    logger.info(
        "Bismark cov: loaded %s -> %d chromosomes, %d sites in %.1fs",
        path.name,
        len(by_chrom),
        sum(methylation_cov_site_count(c) for c in by_chrom.values()),
        time.perf_counter() - t0,
    )
    if use_cache and by_chrom:
        _save_methylation_cov_cache(cache_path, by_chrom, coverage_min=coverage_min)
    return by_chrom


def filter_methylation_cov_by_region(
    by_chrom: dict[str, MethylationChromCov],
    *,
    chrom: str,
    lo_bp: int,
    hi_bp: int,
) -> dict[str, MethylationChromCov]:
    cov = by_chrom.get(chrom)
    if not cov:
        return {}
    pos = cov.get("pos")
    meth = cov.get("meth_pct")
    if pos is None or meth is None or pos.size == 0:
        return {}
    keep = (pos >= lo_bp) & (pos <= hi_bp)
    if not np.any(keep):
        return {}
    return {
        chrom: {
            "pos": np.ascontiguousarray(pos[keep]),
            "meth_pct": np.ascontiguousarray(meth[keep]),
        }
    }


def _genome_fasta_index(genome_fasta: Path):
    path = genome_fasta.expanduser().resolve()
    try:
        mtime = int(path.stat().st_mtime)
    except OSError:
        mtime = 0
    key = (str(path), mtime)
    cached = _GENOME_FASTA_INDEX.get(key)
    if cached is not None:
        return cached
    _GENOME_FASTA_INDEX.clear()
    try:
        index = SeqIO.index(str(path), "fasta")
    except (OSError, ValueError) as exc:
        logger.warning("Bismark cov: cannot index genome %s: %s", path, exc)
        return None
    _GENOME_FASTA_INDEX[key] = index
    return index


def _index_record_id_for_chrom(index, chrom: str) -> str | None:
    for rec_id in index:  # type: ignore[union-attr]
        if _seqid_matches(chrom, rec_id):
            return rec_id
    return None


def _chrom_record_ids_for_genome(
    genome_fasta: Path,
    chrom_sizes: dict[str, int],
) -> dict[str, str]:
    path = genome_fasta.expanduser().resolve()
    try:
        mtime = int(path.stat().st_mtime)
    except OSError:
        mtime = 0
    cache_key = (str(path), mtime)
    cached = _CHROM_RECORD_ID_CACHE.get(cache_key)
    if cached is not None:
        return cached
    index = _genome_fasta_index(path)
    lookup: dict[str, str] = {}
    if index is not None:
        for chrom in chrom_sizes:
            rec_id = _index_record_id_for_chrom(index, chrom)
            if rec_id is not None:
                lookup[chrom] = rec_id
    _CHROM_RECORD_ID_CACHE[cache_key] = lookup
    return lookup


def _chrom_sequence_from_index(genome_fasta: Path, chrom: str) -> str:
    index = _genome_fasta_index(genome_fasta)
    if index is None:
        return ""
    rec_id = _index_record_id_for_chrom(index, chrom)
    if rec_id is None:
        return ""
    return str(index[rec_id].seq)  # type: ignore[index]


def count_cg_in_bins(seq: str, bin_size: int, n_bins: int) -> np.ndarray:
    """Count C and G bases per bin (0-based sequence indexing)."""
    if not seq:
        return np.zeros(n_bins, dtype=np.float64)
    arr = np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)
    cg_mask = (arr == ord("C")) | (arr == ord("G"))
    positions = np.flatnonzero(cg_mask)
    if positions.size == 0:
        return np.zeros(n_bins, dtype=np.float64)
    bin_indices = np.minimum(positions // bin_size, n_bins - 1)
    return np.bincount(bin_indices, minlength=n_bins).astype(np.float64)


def count_cg_bins_from_fasta(
    genome_fasta: Path,
    chrom: str,
    chrom_length: int,
    bin_size: int,
) -> np.ndarray:
    """Count reference C+G bases per window, reading the FASTA in chunks."""
    bin_size_i = max(5_000, int(bin_size))
    n_bins = max(1, int(np.ceil(chrom_length / bin_size_i)))
    path = genome_fasta.expanduser().resolve()
    if not path.is_file():
        return np.zeros(n_bins, dtype=np.float64)
    try:
        mtime = int(path.stat().st_mtime)
    except OSError:
        mtime = 0
    cache_key = (str(path), mtime, chrom, bin_size_i)
    cached = _CG_BIN_CACHE.get(cache_key)
    if cached is not None and cached.size == n_bins:
        return cached

    index = _genome_fasta_index(path)
    if index is None:
        return np.zeros(n_bins, dtype=np.float64)
    rec_lookup = _chrom_record_ids_for_genome(path, {chrom: chrom_length})
    rec_id = rec_lookup.get(chrom)
    if rec_id is None:
        rec_id = _index_record_id_for_chrom(index, chrom)
    if rec_id is None:
        return np.zeros(n_bins, dtype=np.float64)

    record = index[rec_id]
    rec_len = len(record)
    counts = np.zeros(n_bins, dtype=np.float64)
    if rec_len <= 0:
        _CG_BIN_CACHE[cache_key] = counts
        return counts

    chunk_bp = 4_000_000
    c_ord = ord("C")
    g_ord = ord("G")
    for start in range(0, rec_len, chunk_bp):
        end = min(start + chunk_bp, rec_len)
        chunk = str(record.seq[start:end]).upper()
        if not chunk:
            continue
        arr = np.frombuffer(chunk.encode("ascii"), dtype=np.uint8)
        cg_mask = (arr == c_ord) | (arr == g_ord)
        if not np.any(cg_mask):
            continue
        positions = np.flatnonzero(cg_mask) + start
        bin_indices = np.minimum(positions // bin_size_i, n_bins - 1)
        np.add.at(counts, bin_indices, 1.0)

    _CG_BIN_CACHE[cache_key] = counts
    return counts


def bin_cg_total_from_genome(
    genome_fasta: Path | None,
    chrom: str,
    chrom_length: int,
    bin_size: int,
) -> np.ndarray:
    """Track 1: total C+G count in each window from the reference genome."""
    bin_size_i = max(5_000, int(bin_size))
    n_bins = max(1, int(np.ceil(chrom_length / bin_size_i)))
    if genome_fasta is None:
        return np.zeros(n_bins, dtype=np.float64)
    path = genome_fasta.expanduser().resolve()
    if not path.is_file():
        return np.zeros(n_bins, dtype=np.float64)
    try:
        mtime = int(path.stat().st_mtime)
    except OSError:
        mtime = 0
    cache_key = (str(path), mtime, chrom, bin_size_i)
    cached = _CG_BIN_CACHE.get(cache_key)
    if cached is not None and cached.size == n_bins:
        return cached
    counts = count_cg_bins_from_fasta(path, chrom, chrom_length, bin_size_i)
    _CG_BIN_CACHE[cache_key] = counts
    return counts


def _count_cg_bins_worker(payload: dict[str, Any]) -> tuple[str, np.ndarray]:
    chrom = str(payload["chrom"])
    length = int(payload["chrom_length"])
    bin_size = int(payload["bin_size"])
    fasta = Path(str(payload["fasta"]))
    return chrom, count_cg_bins_from_fasta(fasta, chrom, length, bin_size)


def _cg_bins_cache_path(genome_fasta: Path, bin_size: int) -> Path:
    path = genome_fasta.expanduser().resolve()
    bin_size_i = max(5_000, int(bin_size))
    return path.parent / f"{path.name}.cg_bins_{bin_size_i}.repeatviewer.npz"


def _save_cg_bins_cache(cache_path: Path, by_chrom: dict[str, np.ndarray]) -> None:
    payload: dict[str, np.ndarray] = {
        "_version": np.array([_CG_BINS_CACHE_VERSION], dtype=np.int32),
        "_chrom_names": np.asarray(list(by_chrom.keys()), dtype=str),
    }
    for chrom, counts in by_chrom.items():
        safe = chrom.replace("/", "_")
        payload[f"{safe}__cg"] = counts
    try:
        np.savez_compressed(cache_path, **payload)
    except OSError as exc:
        logger.warning("Bismark cov: could not write C+G bin cache %s: %s", cache_path, exc)


def _load_cg_bins_cache(
    cache_path: Path,
    *,
    chromosomes: list[str],
) -> dict[str, np.ndarray] | None:
    try:
        with np.load(cache_path, allow_pickle=False) as data:
            if int(data["_version"][0]) != _CG_BINS_CACHE_VERSION:
                return None
            names = [str(x) for x in data["_chrom_names"].tolist()]
            out: dict[str, np.ndarray] = {}
            for chrom in chromosomes:
                if chrom not in names:
                    return None
                safe = chrom.replace("/", "_")
                key = f"{safe}__cg"
                if key not in data:
                    return None
                out[chrom] = np.asarray(data[key], dtype=np.float64)
            return out
    except (OSError, KeyError, ValueError, TypeError) as exc:
        logger.debug("Bismark cov: C+G bin cache unreadable %s: %s", cache_path, exc)
        return None


def precompute_methylation_cg_totals(
    genome_fasta: Path | None,
    chrom_sizes: dict[str, int],
    bin_size: int,
    *,
    chromosomes: list[str] | None = None,
    num_workers: int | None = None,
) -> dict[str, np.ndarray]:
    """Precompute reference C+G counts per chromosome once for all methylation sources."""
    chroms = chromosomes or list(chrom_sizes.keys())
    bin_size_i = max(5_000, int(bin_size))
    if genome_fasta is None:
        return {
            chrom: np.zeros(max(1, int(np.ceil(chrom_sizes[chrom] / bin_size_i))), dtype=np.float64)
            for chrom in chroms
            if chrom in chrom_sizes and chrom_sizes[chrom] > 0
        }
    fasta_path = genome_fasta.expanduser().resolve()
    cache_path = _cg_bins_cache_path(fasta_path, bin_size_i)
    if cache_path.is_file():
        try:
            if cache_path.stat().st_mtime >= fasta_path.stat().st_mtime:
                cached = _load_cg_bins_cache(cache_path, chromosomes=chroms)
                if cached is not None:
                    logger.info(
                        "Methylation: loaded reference C+G bin cache %s (%d chromosomes)",
                        cache_path.name,
                        len(cached),
                    )
                    return cached
        except OSError:
            pass
    t0 = time.perf_counter()
    tasks: list[dict[str, Any]] = []
    for chrom in chroms:
        length = int(chrom_sizes.get(chrom, 0))
        if length <= 0:
            continue
        tasks.append(
            {
                "chrom": chrom,
                "chrom_length": length,
                "bin_size": bin_size_i,
                "fasta": str(fasta_path),
            }
        )
    out: dict[str, np.ndarray] = {}
    cpus = os.cpu_count() or 4
    cap = max(1, min(int(num_workers or cpus), cpus))
    workers = max(1, min(cap, len(tasks)))
    if workers <= 1 or len(tasks) < 2:
        for task in tasks:
            chrom, counts = _count_cg_bins_worker(task)
            out[chrom] = counts
    else:
        logger.info(
            "Methylation: precomputing reference C+G bins for %d chromosomes with %d workers",
            len(tasks),
            workers,
        )
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_count_cg_bins_worker, task) for task in tasks]
            for future in as_completed(futures):
                chrom, counts = future.result()
                out[chrom] = counts
    if out:
        _save_cg_bins_cache(cache_path, out)
    logger.info(
        "Methylation: precomputed reference C+G bins for %d chromosomes in %.1fs",
        len(out),
        time.perf_counter() - t0,
    )
    return out


def _bin_cov_arrays(
    positions: np.ndarray,
    meth_pct: np.ndarray,
    chrom_length: int,
    bin_size: int,
    *,
    array_workers: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bin .cov positions into per-window line counts and methylation % sums."""
    bin_size_i = max(5_000, int(bin_size))
    n_bins = max(1, int(np.ceil(chrom_length / bin_size_i)))
    centers_mbp = (np.arange(n_bins) + 0.5) * bin_size_i / BP_PER_MBP
    line_counts = np.zeros(n_bins, dtype=np.float64)
    meth_sum = np.zeros(n_bins, dtype=np.float64)
    if not positions.size:
        return centers_mbp, line_counts, meth_sum
    valid = positions >= 1
    pos = positions[valid]
    pct = meth_pct[valid]
    if not pos.size:
        return centers_mbp, line_counts, meth_sum
    bin_indices = np.minimum((pos - 1) // bin_size_i, n_bins - 1)
    threads = max(1, min(int(array_workers), _ARRAY_BIN_MAX_THREADS))
    if pos.size >= _ARRAY_BIN_MIN_SITES and threads > 1:
        splits = np.array_split(np.arange(pos.size, dtype=np.int64), threads)

        def _accumulate_chunk(idxs: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
            bi = bin_indices[idxs]
            weights = pct[idxs]
            return (
                np.bincount(bi, minlength=n_bins).astype(np.float64),
                np.bincount(bi, weights=weights, minlength=n_bins).astype(np.float64),
            )

        with ThreadPoolExecutor(max_workers=threads) as pool:
            for local_counts, local_meth in pool.map(_accumulate_chunk, splits):
                line_counts += local_counts
                meth_sum += local_meth
    else:
        line_counts = np.bincount(bin_indices, minlength=n_bins).astype(np.float64)
        meth_sum = np.bincount(bin_indices, weights=pct, minlength=n_bins).astype(np.float64)
    return centers_mbp, line_counts, meth_sum


def _bin_cov_rows(
    df: pd.DataFrame,
    chrom_length: int,
    bin_size: int,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Bin .cov rows: line counts and sum of methylation % (column 4) per window."""
    if df.empty:
        bin_size_i = max(5_000, int(bin_size))
        n_bins = max(1, int(np.ceil(chrom_length / bin_size_i)))
        centers_mbp = (np.arange(n_bins) + 0.5) * bin_size_i / BP_PER_MBP
        empty = np.zeros(n_bins, dtype=np.float64)
        return centers_mbp, empty, empty
    return _bin_cov_arrays(
        df["pos"].to_numpy(dtype=np.int64, copy=False),
        df["meth_pct"].to_numpy(dtype=np.float64, copy=False),
        chrom_length,
        bin_size,
    )


def _align_cg_total(cg_total: np.ndarray, n_bins: int) -> np.ndarray:
    if cg_total.size == n_bins:
        return cg_total.astype(np.float64, copy=False)
    aligned = np.zeros(n_bins, dtype=np.float64)
    if cg_total.size:
        aligned[: min(n_bins, cg_total.size)] = cg_total[:n_bins]
    return aligned


@dataclass
class MethylationChromBins:
    centers_mbp: np.ndarray
    cg_total: np.ndarray
    line_counts: np.ndarray
    meth_sum: np.ndarray
    bin_size: int


def prepare_methylation_chrom_bins_from_arrays(
    positions: np.ndarray,
    meth_pct: np.ndarray,
    chrom_length: int,
    bin_size: int,
    *,
    cg_total: np.ndarray,
    array_workers: int = 1,
) -> MethylationChromBins:
    """Build prepared methylation bins from numpy arrays and precomputed C+G totals."""
    bin_size_i = max(5_000, int(bin_size))
    centers_mbp, line_counts, meth_sum = _bin_cov_arrays(
        np.asarray(positions, dtype=np.int64),
        np.asarray(meth_pct, dtype=np.float64),
        chrom_length,
        bin_size_i,
        array_workers=array_workers,
    )
    return MethylationChromBins(
        centers_mbp=centers_mbp,
        cg_total=_align_cg_total(np.asarray(cg_total, dtype=np.float64), centers_mbp.size),
        line_counts=line_counts,
        meth_sum=meth_sum,
        bin_size=bin_size_i,
    )


def prepare_methylation_chrom_bins(
    df: pd.DataFrame,
    chrom_length: int,
    bin_size: int,
    *,
    genome_fasta: Path | None = None,
    chrom: str = "",
    cg_total: np.ndarray | None = None,
) -> MethylationChromBins:
    """Bin .cov and genome C+G once per chromosome (reuse for all four subtracks)."""
    bin_size_i = max(5_000, int(bin_size))
    if cg_total is None:
        cg_total = bin_cg_total_from_genome(genome_fasta, chrom, chrom_length, bin_size_i)
    if df.empty:
        positions = np.array([], dtype=np.int64)
        meth_pct = np.array([], dtype=np.float64)
    else:
        positions = df["pos"].to_numpy(dtype=np.int64, copy=False)
        meth_pct = df["meth_pct"].to_numpy(dtype=np.float64, copy=False)
    return prepare_methylation_chrom_bins_from_arrays(
        positions,
        meth_pct,
        chrom_length,
        bin_size_i,
        cg_total=cg_total,
    )


def _density_from_chrom_bins(
    bins: MethylationChromBins,
    subtrack: str,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    mean_pident = np.full(bins.centers_mbp.size, np.nan, dtype=np.float64)
    density = np.zeros(bins.centers_mbp.size, dtype=np.float64)
    if subtrack == METH_SUBTRACK_CG_TOTAL:
        density = bins.cg_total.copy()
        hit_counts = bins.cg_total.copy()
        te_bp = bins.line_counts.copy()
    elif subtrack == METH_SUBTRACK_CG_COVERED:
        density = bins.line_counts.copy()
        hit_counts = bins.line_counts.copy()
        te_bp = bins.cg_total.copy()
    elif subtrack == METH_SUBTRACK_CG_COVERAGE_PCT:
        with np.errstate(divide="ignore", invalid="ignore"):
            density = np.where(
                bins.cg_total > 0,
                bins.line_counts / bins.cg_total,
                0.0,
            )
        mean_pident = density.copy()
        hit_counts = bins.cg_total.copy()
        te_bp = bins.line_counts.copy()
    elif subtrack == METH_SUBTRACK_METH_LEVEL:
        with np.errstate(divide="ignore", invalid="ignore"):
            density = np.where(
                bins.line_counts > 0,
                bins.meth_sum / (100.0 * bins.line_counts),
                0.0,
            )
        mean_pident = density.copy()
        hit_counts = bins.line_counts.copy()
        te_bp = bins.meth_sum.copy()
    else:
        hit_counts = bins.line_counts.copy()
        te_bp = bins.cg_total.copy()
    return density, mean_pident, hit_counts, te_bp


def bin_methylation_subtrack(
    df: pd.DataFrame,
    chrom_length: int,
    bin_size: int,
    *,
    subtrack: str,
    coverage_min: int = 1,
    region: tuple[float, float] | None = None,
    genome_fasta: Path | None = None,
    chrom: str = "",
    cg_total_bins: np.ndarray | None = None,
    chrom_bins: MethylationChromBins | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Bin one methylation subtrack for a chromosome.

    Returns (centers_mbp, density, mean_pident, hit_counts, te_bp) matching TE track bins.
    """
    if chrom_bins is None:
        chrom_bins = prepare_methylation_chrom_bins(
            df,
            chrom_length,
            bin_size,
            genome_fasta=genome_fasta,
            chrom=chrom,
        )
        if cg_total_bins is not None and cg_total_bins.size == chrom_bins.cg_total.size:
            chrom_bins = MethylationChromBins(
                centers_mbp=chrom_bins.centers_mbp,
                cg_total=cg_total_bins.astype(np.float64, copy=False),
                line_counts=chrom_bins.line_counts,
                meth_sum=chrom_bins.meth_sum,
                bin_size=chrom_bins.bin_size,
            )
    density, mean_pident, hit_counts, te_bp = _density_from_chrom_bins(chrom_bins, subtrack)
    if region is not None:
        from repeatviewer.backends.genome_map_binning import _filter_bins_to_mbp_window

        return _filter_bins_to_mbp_window(
            chrom_bins.centers_mbp,
            density,
            mean_pident,
            hit_counts,
            te_bp,
            chrom_bins.bin_size,
            region[0],
            region[1],
        )
    return chrom_bins.centers_mbp, density, mean_pident, hit_counts, te_bp


def methylation_level_density(bins: MethylationChromBins) -> np.ndarray:
    """Per-bin methylation level (0–1) from prepared chromosome bins."""
    density, _, _, _ = _density_from_chrom_bins(bins, METH_SUBTRACK_METH_LEVEL)
    return density


def methylation_cov_payload(df: pd.DataFrame) -> dict[str, np.ndarray]:
    """Serialize a chromosome .cov slice for ProcessPool workers."""
    if df.empty:
        return {}
    return {
        "pos": np.ascontiguousarray(df["pos"].to_numpy(dtype=np.int64, copy=False)),
        "meth_pct": np.ascontiguousarray(df["meth_pct"].to_numpy(dtype=np.float64, copy=False)),
    }


def _cg_total_from_payload(payload: dict[str, Any], key: str = "cg_total") -> np.ndarray:
    raw = payload.get(key)
    if raw is None:
        return np.array([], dtype=np.float64)
    return np.asarray(raw, dtype=np.float64)


def _methylation_cov_from_payload(data: dict[str, np.ndarray] | None) -> tuple[np.ndarray, np.ndarray]:
    if data is None or len(data) == 0:
        return np.array([], dtype=np.int64), np.array([], dtype=np.float64)
    pos = data.get("pos")
    meth_pct = data.get("meth_pct")
    if pos is None or meth_pct is None:
        return np.array([], dtype=np.int64), np.array([], dtype=np.float64)
    return np.asarray(pos, dtype=np.int64), np.asarray(meth_pct, dtype=np.float64)


# Backward-compatible alias.
methylation_df_payload = methylation_cov_payload


def _methylation_chrom_bins_from_cov(
    cov_payload: dict[str, np.ndarray] | None,
    *,
    chrom_length: int,
    bin_size: int,
    cg_total: np.ndarray,
    array_workers: int = 1,
) -> MethylationChromBins:
    positions, meth_pct = _methylation_cov_from_payload(cov_payload)
    return prepare_methylation_chrom_bins_from_arrays(
        positions,
        meth_pct,
        chrom_length,
        bin_size,
        cg_total=cg_total,
        array_workers=array_workers,
    )


def bin_methylation_chromosome_worker(
    payload: dict[str, Any],
) -> tuple[str, dict[str, tuple], dict[str, MethylationChromBins]]:
    """ProcessPool worker: bin all methylation subtracks on one chromosome."""
    chrom = str(payload.get("chrom") or "")
    chrom_length = int(payload["chrom_length"])
    bin_size = int(payload["bin_size"])
    coverage_min = int(payload.get("coverage_min", 1))
    region = payload.get("region")
    array_workers = int(payload.get("array_workers", 1))
    source_cov: dict[str, dict[str, np.ndarray]] = payload.get("source_cov") or {}
    source_cov_cache_paths: dict[str, str] = payload.get("source_cov_cache_paths") or {}
    source_bins_cache_paths: dict[str, str] = payload.get("source_bins_cache_paths") or {}
    track_specs: list[dict[str, str]] = list(payload.get("track_specs") or [])
    cg_total = _cg_total_from_payload(payload)

    by_source: dict[str, list[dict[str, str]]] = {}
    for spec in track_specs:
        source_key = str(spec["source_key"])
        by_source.setdefault(source_key, []).append(spec)

    track_bins: dict[str, tuple] = {}
    binned_by_source: dict[str, MethylationChromBins] = {}
    for source_key, specs in by_source.items():
        chrom_bins: MethylationChromBins | None = None
        bins_cache_path = source_bins_cache_paths.get(source_key)
        if bins_cache_path and chrom:
            chrom_bins = _load_methylation_bins_chrom_cache(
                Path(bins_cache_path),
                chrom,
                bin_size=bin_size,
                cg_total=cg_total,
            )
        if chrom_bins is None:
            cov_payload = source_cov.get(source_key)
            if (not cov_payload) and source_key in source_cov_cache_paths and chrom:
                cov_payload = load_methylation_cov_chrom_from_cache(
                    Path(source_cov_cache_paths[source_key]),
                    chrom,
                )
            chrom_bins = _methylation_chrom_bins_from_cov(
                cov_payload,
                chrom_length=chrom_length,
                bin_size=bin_size,
                cg_total=cg_total,
                array_workers=array_workers,
            )
            if chrom_bins is not None:
                binned_by_source[source_key] = chrom_bins
        for spec in specs:
            track_bins[str(spec["track_label"])] = bin_methylation_subtrack(
                pd.DataFrame(),
                chrom_length,
                bin_size,
                subtrack=str(spec["subtrack"]),
                coverage_min=coverage_min,
                region=region,
                chrom_bins=chrom_bins,
            )
    return chrom, track_bins, binned_by_source


def methylation_diff_bins_worker(
    payload: dict[str, Any],
) -> tuple[str, np.ndarray, np.ndarray, str, str]:
    """ProcessPool worker: methylation level diff for one chromosome."""
    chrom = str(payload["chrom"])
    chrom_length = int(payload["chrom_length"])
    bin_size = int(payload["bin_size"])
    region = payload.get("region")
    array_workers = int(payload.get("array_workers", 1))
    label_a = str(payload["label_a"])
    label_b = str(payload["label_b"])
    label = str(payload["label"])
    cg_total = _cg_total_from_payload(payload)
    source_cov_cache_paths: dict[str, str] = payload.get("source_cov_cache_paths") or {}
    cov_a = payload.get("cov_a")
    cov_b = payload.get("cov_b")
    if (not cov_a) and "a" in source_cov_cache_paths and chrom:
        cov_a = load_methylation_cov_chrom_from_cache(Path(source_cov_cache_paths["a"]), chrom)
    if (not cov_b) and "b" in source_cov_cache_paths and chrom:
        cov_b = load_methylation_cov_chrom_from_cache(Path(source_cov_cache_paths["b"]), chrom)
    bins_a = _methylation_chrom_bins_from_cov(
        cov_a,
        chrom_length=chrom_length,
        bin_size=bin_size,
        cg_total=cg_total,
        array_workers=array_workers,
    )
    bins_b = _methylation_chrom_bins_from_cov(
        cov_b,
        chrom_length=chrom_length,
        bin_size=bin_size,
        cg_total=cg_total,
        array_workers=array_workers,
    )
    level_a = methylation_level_density(bins_a)
    level_b = methylation_level_density(bins_b)
    diff_a = np.maximum(level_a - level_b, 0.0)
    diff_b = np.maximum(level_b - level_a, 0.0)
    if region is not None:
        from repeatviewer.backends.genome_map_binning import _filter_bins_to_mbp_window

        centers = bins_a.centers_mbp
        _, diff_a, _, _, _ = _filter_bins_to_mbp_window(
            centers,
            diff_a,
            diff_a,
            diff_a,
            diff_a,
            bins_a.bin_size,
            region[0],
            region[1],
        )
        _, diff_b, _, _, _ = _filter_bins_to_mbp_window(
            centers,
            diff_b,
            diff_b,
            diff_b,
            diff_b,
            bins_a.bin_size,
            region[0],
            region[1],
        )
    return label, diff_a, diff_b, label_a, label_b


def shared_cg_count_axis_max(
    methylation_norm_max: dict[str, float],
) -> float:
    """Tracks 1 and 2 share one y-axis scale."""
    return max(
        float(methylation_norm_max.get(METH_SUBTRACK_CG_TOTAL, 1.0)),
        float(methylation_norm_max.get(METH_SUBTRACK_CG_COVERED, 1.0)),
        1.0,
    )
