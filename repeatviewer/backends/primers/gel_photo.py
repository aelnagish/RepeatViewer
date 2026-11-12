"""Read a photo of a real agarose gel: lanes, bands, ladder fit, which band is which.

Qt-free; works on a 2-D numpy array (grey values, any range). Bands are made
bright (a light-background scan is inverted), lanes are found from the column
profile, bands from each lane's row profile. The ladder lane gives a fit of
log10(size) against distance (quadratic), which turns every other band into a
size; those sizes are then compared with the bands the in-silico gel expects.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np


def _smooth(x: np.ndarray, width: int) -> np.ndarray:
    width = max(1, int(width))
    if width <= 1:
        return x.astype(float)
    k = np.ones(width) / width
    return np.convolve(x.astype(float), k, mode="same")


def _peaks(x: np.ndarray, *, min_dist: int, threshold: float) -> list[int]:
    """Local maxima above ``threshold``, at least ``min_dist`` apart (strongest kept)."""
    if len(x) < 3:
        return []
    cand = [i for i in range(1, len(x) - 1) if x[i] >= x[i - 1] and x[i] > x[i + 1]
            and x[i] >= threshold]
    cand.sort(key=lambda i: -x[i])
    kept: list[int] = []
    for i in cand:
        if all(abs(i - j) >= min_dist for j in kept):
            kept.append(i)
    return sorted(kept)


def normalise(img: np.ndarray) -> np.ndarray:
    """Grey 0..1 with bands bright (inverts light-background scans)."""
    a = np.asarray(img, dtype=float)
    if a.ndim == 3:
        a = a[..., :3].mean(axis=2)
    lo, hi = np.percentile(a, 1), np.percentile(a, 99.5)
    a = np.clip((a - lo) / max(1e-9, hi - lo), 0, 1)
    if np.median(a) > 0.5:
        a = 1.0 - a
    return a


@dataclass
class Lane:
    x: int          # centre column
    half: int       # half width


def find_lanes(img: np.ndarray, *, max_lanes: int = 30) -> list[Lane]:
    """Lanes as runs of bright columns (separated by darker gaps)."""
    a = normalise(img)
    h, w = a.shape
    body = a[int(h * 0.05):int(h * 0.95)]
    col = _smooth(np.percentile(body, 99.5, axis=0), max(3, w // 100))
    lo, hi = float(np.percentile(col, 5)), float(col.max())
    if hi - lo <= 1e-6:
        return []
    thr = lo + 0.3 * (hi - lo)
    above = col >= thr
    runs: list[list[int]] = []
    start = None
    for i, on in enumerate(above):
        if on and start is None:
            start = i
        elif not on and start is not None:
            runs.append([start, i - 1])
            start = None
    if start is not None:
        runs.append([start, w - 1])
    # merge tiny gaps, drop slivers
    merged: list[list[int]] = []
    for r in runs:
        if merged and r[0] - merged[-1][1] <= max(2, w // 150):
            merged[-1][1] = r[1]
        else:
            merged.append(r)
    lanes = [Lane(x=(a0 + b0) // 2, half=max(2, int((b0 - a0) * 0.4)))
             for a0, b0 in merged if b0 - a0 >= max(3, w // 70)]
    return lanes[:max_lanes]


@dataclass
class Band:
    y: int
    strength: float


def lane_profile(img: np.ndarray, lane: Lane) -> np.ndarray:
    a = normalise(img)
    x0 = max(0, lane.x - lane.half)
    x1 = min(a.shape[1], lane.x + lane.half + 1)
    prof = a[:, x0:x1].mean(axis=1)
    # remove the slowly varying background (smear, uneven light)
    win = max(15, a.shape[0] // 12)
    base = np.array([np.percentile(prof[max(0, i - win):i + win], 20) for i in range(len(prof))])
    return _smooth(prof - base, max(2, a.shape[0] // 250))


def find_bands(img: np.ndarray, lane: Lane, *, min_strength: float = 0.04) -> list[Band]:
    prof = lane_profile(img, lane)
    if prof.size == 0:
        return []
    mad = float(np.median(np.abs(prof - np.median(prof)))) or 1e-3
    thr = max(min_strength, np.median(prof) + 6 * mad, 0.08 * float(prof.max()))
    peaks = _peaks(prof, min_dist=max(3, len(prof) // 100), threshold=thr)
    return [Band(y=int(p), strength=float(prof[p])) for p in peaks]


@dataclass
class LadderFit:
    """Size from distance: point-to-point on log size between the ladder bands."""

    used: list[tuple[int, int]] = field(default_factory=list)   # (y, size), y increasing
    max_error: float = 0.0               # relative error of a smooth fit (quality)

    def _arrays(self):
        ys = np.array([y for y, _s in self.used], dtype=float)
        ts = np.log10(np.array([s for _y, s in self.used], dtype=float))
        return ys, ts

    def size_at(self, y: float) -> float:
        ys, ts = self._arrays()
        if y <= ys[0]:
            k = (ts[1] - ts[0]) / (ys[1] - ys[0])
            return float(10 ** (ts[0] + k * (y - ys[0])))
        if y >= ys[-1]:
            k = (ts[-1] - ts[-2]) / (ys[-1] - ys[-2])
            return float(10 ** (ts[-1] + k * (y - ys[-1])))
        return float(10 ** np.interp(y, ys, ts))

    def y_at(self, size: float) -> float:
        ys, ts = self._arrays()
        t = math.log10(max(1.0, size))
        # ts falls with y: interpolate on the reversed arrays
        if t >= ts[0]:
            k = (ys[1] - ys[0]) / (ts[1] - ts[0])
            return float(ys[0] + k * (t - ts[0]))
        if t <= ts[-1]:
            k = (ys[-1] - ys[-2]) / (ts[-1] - ts[-2])
            return float(ys[-1] + k * (t - ts[-1]))
        return float(np.interp(t, ts[::-1], ys[::-1]))


def fit_ladder(band_ys: list[int], sizes: list[int]) -> LadderFit | None:
    """Match bands (top to bottom) to ladder sizes (large to small).

    Ladder bands may be missing at either end (run off, merged near the wells):
    every contiguous window of the size list is tried; the one a smooth cubic
    fits best (relative error) wins, with a small cost for each band left out.
    """
    ys = sorted(band_ys)
    sz = sorted(sizes, reverse=True)
    if len(ys) < 3 or len(sz) < 3:
        return None
    n = min(len(ys), len(sz))
    best, best_score = None, None
    for k in range(n, max(2, n - 4), -1):
        for ys_off in range(0, len(ys) - k + 1):
            for sz_off in range(0, len(sz) - k + 1):
                yy = ys[ys_off:ys_off + k]
                ss = sz[sz_off:sz_off + k]
                y = np.array(yy, dtype=float)
                t = np.log10(np.array(ss, dtype=float))
                deg = 3 if k >= 6 else (2 if k >= 4 else 1)
                c = np.polyfit(y, t, deg)
                pred = np.polyval(c, y)
                err = float(np.max(np.abs(10 ** (pred - t) - 1)))
                score = err + 0.015 * (n - k)
                if best_score is None or score < best_score:
                    best, best_score = LadderFit(used=list(zip(yy, ss)), max_error=err), score
    return best


@dataclass
class BandCall:
    lane: int
    y: int
    size: float
    status: str           # matched / brighter / unexplained / missing
    expected: int = 0     # matched simulated size
    note: str = ""


def compare(observed: list[float], expected: list[tuple[int, str, float]], *,
            tolerance: float = 0.07) -> list[tuple[float | None, int | None, str, str]]:
    """Pair observed sizes with expected bands (size, label, relative amount).

    Returns rows (observed size or None, expected size or None, status, label).
    """
    rows = []
    used: set[int] = set()
    for obs in sorted(observed, reverse=True):
        best, best_err = None, None
        for i, (size, _label, _amt) in enumerate(expected):
            if i in used:
                continue
            err = abs(obs / size - 1)
            if err <= tolerance and (best_err is None or err < best_err):
                best, best_err = i, err
        if best is None:
            rows.append((obs, None, "unexplained", ""))
        else:
            used.add(best)
            size, label, amount = expected[best]
            status = "matched"
            if label == "off-target":
                status = "off-target"
            rows.append((obs, size, status, label))
    for i, (size, label, amount) in enumerate(expected):
        if i not in used and amount >= 0.05:
            rows.append((None, size, "not seen", label))
    return rows
