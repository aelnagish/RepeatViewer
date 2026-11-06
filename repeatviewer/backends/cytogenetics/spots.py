"""Find FISH signals in one probe channel and place them on detected chromosomes. No Qt.

A signal is a local maximum of a difference-of-Gaussians image (spot-sized blur minus the
background around it) above a threshold set from the channel's own noise (median + k x MAD).
Sensitivity moves k; ``smallest_px`` drops tiny specks; spots closer than ``join_px``
count as one. Each spot is then given to the chromosome it sits on.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

try:
    from scipy import ndimage as ndi
except ImportError:  # pragma: no cover
    ndi = None


@dataclass
class SpotSettings:
    sensitivity: float = 0.5  # 0..1, higher finds weaker spots
    smallest_px: float = 2.0  # spot diameter in pixels
    join_px: float = 3.0  # spots closer than this are one

    def to_dict(self) -> dict:
        return {"sensitivity": self.sensitivity, "smallest_px": self.smallest_px, "join_px": self.join_px}

    @classmethod
    def from_dict(cls, raw) -> "SpotSettings":
        raw = raw if isinstance(raw, dict) else {}
        out = cls()
        for key, lo, hi in (("sensitivity", 0.0, 1.0), ("smallest_px", 0.5, 50.0), ("join_px", 0.0, 100.0)):
            try:
                setattr(out, key, max(lo, min(hi, float(raw.get(key, getattr(out, key))))))
            except (TypeError, ValueError):
                pass
        return out


@dataclass
class Spot:
    x: float
    y: float
    peak: float  # height above the threshold, in thresholds (1 = just found)
    intensity: float  # summed signal of the spot
    area: int


def _grey(array: np.ndarray) -> np.ndarray:
    a = np.asarray(array, dtype=np.float64)
    if a.ndim == 3:
        a = a[..., :3].max(axis=2)
    return a


def response(array: np.ndarray, spot_px: float) -> np.ndarray:
    """Difference of Gaussians tuned to the spot size."""
    g = _grey(array)
    s1 = max(0.7, spot_px / 2.5)
    return ndi.gaussian_filter(g, s1) - ndi.gaussian_filter(g, s1 * 4.0)


def threshold(dog: np.ndarray, sensitivity: float, where: np.ndarray | None = None) -> float:
    vals = dog[where] if where is not None and where.any() else dog
    vals = vals[:: max(1, vals.size // 400_000)] if vals.ndim == 1 else vals.ravel()[:: max(1, vals.size // 400_000)]
    med = float(np.median(vals))
    mad = float(np.median(np.abs(vals - med))) * 1.4826 or float(vals.std()) or 1.0
    k = 3.0 + (1.0 - max(0.0, min(1.0, sensitivity))) * 14.0  # 0.5 -> 10 MAD
    return med + k * mad


def find_spots(array: np.ndarray, settings: SpotSettings | None = None, chromosomes: np.ndarray | None = None) -> list[Spot]:
    """Spots of one channel. ``chromosomes`` (bool or label map) limits them to chromosomes and their edge."""
    if ndi is None:
        raise RuntimeError("Signal detection needs scipy (pip install scipy)")
    settings = settings or SpotSettings()
    dog = response(array, settings.smallest_px)
    near = None
    if chromosomes is not None and np.any(chromosomes):
        near = ndi.binary_dilation(np.asarray(chromosomes) > 0, iterations=4)
    thr = threshold(dog, settings.sensitivity)
    if thr <= 0:
        thr = float(dog.max()) * 0.5 or 1.0
    above = dog > thr
    if near is not None:
        above &= near
    if not above.any():
        return []
    size = max(3, int(round(2 * max(settings.join_px, settings.smallest_px / 2))) | 1)
    peaks = above & (dog == ndi.maximum_filter(dog, size=size))
    regions, n = ndi.label(dog > thr * 0.5)
    min_area = max(1, int(round(math.pi * (settings.smallest_px / 2) ** 2)))
    areas = np.bincount(regions.ravel(), minlength=n + 1)
    sums = ndi.sum(np.clip(dog, 0, None), regions, np.arange(n + 1)) if n else np.zeros(1)
    out: list[Spot] = []
    for y, x in np.argwhere(peaks):
        k = regions[y, x]
        if k == 0 or areas[k] < min_area:
            continue
        # centre of mass in a small window
        r = max(2, int(settings.smallest_px))
        y0, y1, x0, x1 = max(0, y - r), y + r + 1, max(0, x - r), x + r + 1
        w = np.clip(dog[y0:y1, x0:x1] - thr * 0.5, 0, None)
        if w.sum() > 0:
            yy, xx = np.mgrid[y0:y0 + w.shape[0], x0:x0 + w.shape[1]]
            cx, cy = float((w * xx).sum() / w.sum()), float((w * yy).sum() / w.sum())
        else:
            cx, cy = float(x), float(y)
        out.append(Spot(cx, cy, float(dog[y, x] / thr), float(sums[k]), int(areas[k])))
    # one spot per region and per join distance: keep the strongest
    out.sort(key=lambda s: -s.peak)
    kept: list[Spot] = []
    seen_regions: set[int] = set()
    for s in out:
        k = int(regions[int(round(s.y)), int(round(s.x))]) if 0 <= int(round(s.y)) < regions.shape[0] and 0 <= int(round(s.x)) < regions.shape[1] else 0
        if any(math.hypot(s.x - o.x, s.y - o.y) < max(settings.join_px, 1.0) for o in kept):
            continue
        if k and k in seen_regions and areas[k] < 12 * min_area:
            continue
        seen_regions.add(k)
        kept.append(s)
    return kept


def strength(spot: Spot) -> str:
    if spot.peak >= 3.0:
        return "strong"
    if spot.peak >= 1.6:
        return "medium"
    return "weak"


def assign(spot_xy: tuple[float, float], labels: np.ndarray, max_distance: float = 6.0) -> int:
    """Index (0-based) of the chromosome in a label map under or nearest to a point; -1 = none."""
    x, y = spot_xy
    h, w = labels.shape
    xi, yi = int(round(x)), int(round(y))
    if 0 <= yi < h and 0 <= xi < w and labels[yi, xi] > 0:
        return int(labels[yi, xi]) - 1
    r = int(math.ceil(max_distance))
    y0, y1, x0, x1 = max(0, yi - r), min(h, yi + r + 1), max(0, xi - r), min(w, xi + r + 1)
    win = labels[y0:y1, x0:x1]
    ys, xs = np.nonzero(win)
    if not len(xs):
        return -1
    d = np.hypot(xs + x0 - x, ys + y0 - y)
    k = int(np.argmin(d))
    return int(win[ys[k], xs[k]]) - 1 if d[k] <= max_distance else -1
