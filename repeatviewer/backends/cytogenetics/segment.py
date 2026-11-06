"""Find chromosomes in a counterstain (DAPI) image and measure them. No Qt.

``detect_chromosomes`` thresholds the counterstain, separates touching chromosomes
(distance-map watershed, then merges pieces that together form one unbranched
chromosome of normal size), and measures each object: its midline (the longest path
of its skeleton, extended to the ends), length, width along the midline and the
centromere (the narrowest place between 15 and 85 % of the length).

Objects are stored compactly as a bounding box plus run lengths (``encode_mask``).
The edit functions (split along a line, merge, brush) work on those masks and are
followed by ``measure`` again.
"""

from __future__ import annotations

import math
from collections import deque
from dataclasses import dataclass, field

import numpy as np

try:
    from scipy import ndimage as ndi
except ImportError:  # pragma: no cover - scipy comes with scikit-learn
    ndi = None

CHECK_BRANCHED = "branched, maybe two chromosomes"
CHECK_LARGE = "large, maybe two touching chromosomes"
CHECK_SMALL = "small, maybe a piece"
CHECK_CENTROMERE = "centromere not clear"
CHECK_EDGE = "cut by the image edge"


@dataclass
class DetectSettings:
    sensitivity: float = 0.5  # 0..1, higher finds fainter chromosomes
    separate: float = 0.5  # 0..1, higher splits touching chromosomes harder
    smallest_px: int = 0  # objects smaller than this area are dropped; 0 = automatic

    def to_dict(self) -> dict:
        return {"sensitivity": self.sensitivity, "separate": self.separate, "smallest_px": self.smallest_px}

    @classmethod
    def from_dict(cls, raw) -> "DetectSettings":
        raw = raw if isinstance(raw, dict) else {}
        out = cls()
        for key in ("sensitivity", "separate"):
            try:
                out.__setattr__(key, max(0.0, min(1.0, float(raw.get(key, getattr(out, key))))))
            except (TypeError, ValueError):
                pass
        try:
            out.smallest_px = max(0, int(raw.get("smallest_px", 0)))
        except (TypeError, ValueError):
            pass
        return out


@dataclass
class Measured:
    """One chromosome: its mask (full-image coordinates through the box) and measurements."""

    mask: np.ndarray  # bool crop
    x0: int
    y0: int
    midline: list[tuple[float, float]] = field(default_factory=list)  # (x, y) from p telomere to q telomere
    centromere: int = -1  # index into midline, -1 = none
    widths: list[float] = field(default_factory=list)
    checks: list[str] = field(default_factory=list)

    @property
    def area(self) -> int:
        return int(self.mask.sum())

    @property
    def length(self) -> float:
        return polyline_length(self.midline)

    def centre(self) -> tuple[float, float]:
        ys, xs = np.nonzero(self.mask)
        if not len(xs):
            return float(self.x0), float(self.y0)
        return float(xs.mean() + self.x0), float(ys.mean() + self.y0)


# ---------------------------------------------------------------------------
# small helpers
# ---------------------------------------------------------------------------


def polyline_length(points) -> float:
    if len(points) < 2:
        return 0.0
    a = np.asarray(points, dtype=np.float64)
    return float(np.hypot(*(a[1:] - a[:-1]).T).sum())


def normalise(grey: np.ndarray) -> np.ndarray:
    """Float image scaled so the background is near 0 and bright chromosomes near 1."""
    g = np.asarray(grey, dtype=np.float64)
    if g.ndim == 3:
        g = g.max(axis=2)
    sample = g[:: max(1, g.shape[0] // 512), :: max(1, g.shape[1] // 512)]
    lo, hi = np.percentile(sample, (1.0, 99.7))
    if hi <= lo:
        hi = lo + 1.0
    return np.clip((g - lo) / (hi - lo), 0.0, 1.5)


def otsu(values: np.ndarray) -> float:
    data = np.asarray(values, dtype=np.float64).ravel()
    lo, hi = float(data.min()), float(data.max())
    if hi <= lo:
        return hi
    hist, edges = np.histogram(data, bins=256, range=(lo, hi))
    centers = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * centers)
    mean0 = m0 / np.maximum(w0, 1)
    mean1 = (m0[-1] - m0) / np.maximum(w1, 1)
    between = w0 * w1 * (mean0 - mean1) ** 2
    return float(centers[int(np.argmax(between))])


def encode_mask(mask: np.ndarray, x0: int, y0: int) -> dict:
    """Bounding box + run lengths (starting with a run of False) of a bool crop."""
    m = np.asarray(mask, dtype=bool)
    flat = m.ravel()
    if not flat.size:
        return {"x": int(x0), "y": int(y0), "w": 0, "h": 0, "rle": []}
    change = np.flatnonzero(flat[1:] != flat[:-1]) + 1
    bounds = np.concatenate(([0], change, [flat.size]))
    runs = np.diff(bounds).tolist()
    if flat[0]:
        runs = [0] + runs
    return {"x": int(x0), "y": int(y0), "w": int(m.shape[1]), "h": int(m.shape[0]), "rle": [int(r) for r in runs]}


def decode_mask(raw) -> tuple[np.ndarray, int, int] | None:
    """(bool crop, x0, y0) from ``encode_mask`` output; None for bad input."""
    if not isinstance(raw, dict):
        return None
    try:
        w, h, x0, y0 = int(raw["w"]), int(raw["h"]), int(raw["x"]), int(raw["y"])
        runs = [int(r) for r in raw.get("rle") or []]
    except (KeyError, TypeError, ValueError):
        return None
    if w <= 0 or h <= 0 or sum(runs) != w * h:
        return None
    values = np.zeros(len(runs), dtype=bool)
    values[1::2] = True
    flat = np.repeat(values, runs)
    return flat.reshape(h, w), x0, y0


def _crop(mask: np.ndarray, pad: int = 2) -> tuple[np.ndarray, int, int]:
    ys, xs = np.nonzero(mask)
    if not len(xs):
        return np.zeros((0, 0), dtype=bool), 0, 0
    y0, y1 = max(0, ys.min() - pad), min(mask.shape[0], ys.max() + pad + 1)
    x0, x1 = max(0, xs.min() - pad), min(mask.shape[1], xs.max() + pad + 1)
    return mask[y0:y1, x0:x1].copy(), int(x0), int(y0)


def paste(shape: tuple[int, int], crop: np.ndarray, x0: int, y0: int) -> np.ndarray:
    """Full-image bool mask from a crop (parts outside the image are cut)."""
    out = np.zeros(shape, dtype=bool)
    h, w = crop.shape
    ya, xa = max(0, y0), max(0, x0)
    yb, xb = min(shape[0], y0 + h), min(shape[1], x0 + w)
    if yb > ya and xb > xa:
        out[ya:yb, xa:xb] = crop[ya - y0:yb - y0, xa - x0:xb - x0]
    return out


# ---------------------------------------------------------------------------
# skeleton and midline
# ---------------------------------------------------------------------------


def thin(mask: np.ndarray) -> np.ndarray:
    """Zhang-Suen thinning (vectorised)."""
    img = np.pad(np.asarray(mask, dtype=np.uint8), 1)
    while True:
        changed = False
        for step in (0, 1):
            c = img[1:-1, 1:-1]
            p2, p3, p4 = img[:-2, 1:-1], img[:-2, 2:], img[1:-1, 2:]
            p5, p6, p7 = img[2:, 2:], img[2:, 1:-1], img[2:, :-2]
            p8, p9 = img[1:-1, :-2], img[:-2, :-2]
            seq = (p2, p3, p4, p5, p6, p7, p8, p9, p2)
            b = p2.astype(np.int16) + p3 + p4 + p5 + p6 + p7 + p8 + p9
            a = sum(((seq[i] == 0) & (seq[i + 1] == 1)).astype(np.int16) for i in range(8))
            if step == 0:
                m = (p2 * p4 * p6 == 0) & (p4 * p6 * p8 == 0)
            else:
                m = (p2 * p4 * p8 == 0) & (p2 * p6 * p8 == 0)
            remove = (c == 1) & (b >= 2) & (b <= 6) & (a == 1) & m
            if remove.any():
                c[remove] = 0
                changed = True
        if not changed:
            break
    return img[1:-1, 1:-1].astype(bool)


_NEIGH = [(-1, -1), (-1, 0), (-1, 1), (0, -1), (0, 1), (1, -1), (1, 0), (1, 1)]


def _bfs(start, pixels: set) -> tuple[dict, tuple]:
    prev = {start: None}
    q = deque([start])
    last = start
    while q:
        cur = q.popleft()
        last = cur
        r, c = cur
        for dr, dc in _NEIGH:
            nb = (r + dr, c + dc)
            if nb in pixels and nb not in prev:
                prev[nb] = cur
                q.append(nb)
    return prev, last


def longest_path(skeleton: np.ndarray) -> list[tuple[int, int]]:
    """(row, col) pixels of the longest path through a skeleton (two breadth-first searches)."""
    pts = np.argwhere(skeleton)
    if not len(pts):
        return []
    pixels = {(int(r), int(c)) for r, c in pts}
    _prev, far = _bfs((int(pts[0][0]), int(pts[0][1])), pixels)
    prev, end = _bfs(far, pixels)
    path = []
    cur = end
    while cur is not None:
        path.append(cur)
        cur = prev[cur]
    return path[::-1]


def _smooth(points: np.ndarray, window: int) -> np.ndarray:
    if len(points) < 3 or window < 2:
        return points
    window = min(window, len(points))
    pad = window // 2
    padded = np.concatenate([np.repeat(points[:1], pad, axis=0), points, np.repeat(points[-1:], pad, axis=0)])
    kernel = np.ones(window) / window
    out = np.column_stack([np.convolve(padded[:, i], kernel, mode="valid") for i in range(2)])
    out = out[: len(points)]
    out[0], out[-1] = points[0], points[-1]
    return out


def _resample(points: np.ndarray, step: float = 1.0) -> np.ndarray:
    if len(points) < 2:
        return points
    seg = np.hypot(*(points[1:] - points[:-1]).T)
    s = np.concatenate(([0.0], np.cumsum(seg)))
    if s[-1] <= 0:
        return points[:1]
    n = max(2, int(round(s[-1] / step)) + 1)
    t = np.linspace(0.0, s[-1], n)
    return np.column_stack([np.interp(t, s, points[:, 0]), np.interp(t, s, points[:, 1])])


def _extend(points: np.ndarray, mask: np.ndarray, back: int = 6) -> np.ndarray:
    """Lengthen both ends along their direction until they leave the mask (skeletons stop short)."""
    h, w = mask.shape

    def grow(pts: np.ndarray) -> list:
        if len(pts) < 2:
            return []
        end = pts[-1]
        ref = pts[max(0, len(pts) - 1 - back)]
        d = end - ref
        n = math.hypot(*d)
        if n == 0:
            return []
        d = d / n
        extra = []
        p = end.copy()
        for _ in range(200):
            p = p + d * 0.5
            x, y = int(round(p[0])), int(round(p[1]))
            if not (0 <= x < w and 0 <= y < h) or not mask[y, x]:
                break
            extra.append(p.copy())
        return extra

    tail = grow(points)
    head = grow(points[::-1])
    parts = [np.array(head[::-1]).reshape(-1, 2), points, np.array(tail).reshape(-1, 2)]
    return np.concatenate(parts)


def measure(mask: np.ndarray, x0: int, y0: int, image_shape: tuple[int, int] | None = None,
            grey: np.ndarray | None = None) -> Measured:
    """Midline, widths and centromere of one object (crop + offset).

    The centromere is the place where the chromosome is narrowest and (with ``grey``, the full counterstain)
    dimmest, between 15 and 85 % of its length.
    """
    m = np.asarray(mask, dtype=bool)
    out = Measured(mask=m, x0=int(x0), y0=int(y0))
    if m.sum() < 6:
        return out
    if ndi is not None:
        m = ndi.binary_fill_holes(m)
        out.mask = m
        dist = ndi.distance_transform_edt(np.pad(m, 1))[1:-1, 1:-1]
    else:
        dist = m.astype(np.float64)
    skel = thin(m)
    path = longest_path(skel)
    if len(path) < 2:
        ys, xs = np.nonzero(m)
        pts = np.column_stack([xs, ys]).astype(np.float64)
        centre = pts.mean(axis=0)
        vals, vecs = np.linalg.eigh(np.cov((pts - centre).T)) if len(pts) > 2 else (np.ones(2), np.eye(2))
        v = vecs[:, int(np.argmax(vals))]
        t = (pts - centre) @ v
        line = np.array([centre + v * t.min(), centre + v * t.max()])
    else:
        line = np.array([(c, r) for r, c in path], dtype=np.float64)
        line = _smooth(line, max(3, int(round(float(np.median(dist[skel])) if skel.any() else 3)) * 2 + 1))
        line = _resample(line, 1.0)
        line = _extend(line, m)
        line = _resample(line, 1.0)
    rows = np.clip(np.round(line[:, 1]).astype(int), 0, m.shape[0] - 1)
    cols = np.clip(np.round(line[:, 0]).astype(int), 0, m.shape[1] - 1)
    widths = 2.0 * dist[rows, cols]
    n = len(line)
    shape_score = widths / (float(np.median(widths)) or 1.0)
    if grey is not None:
        gy = np.clip(rows + y0, 0, grey.shape[0] - 1)
        gx = np.clip(cols + x0, 0, grey.shape[1] - 1)
        bright = np.asarray(grey, dtype=np.float64)[gy, gx]
        shape_score = shape_score * np.sqrt(bright / (float(np.median(bright)) or 1.0))
    cen = -1
    unclear = True
    if n >= 10:
        win = max(3, n // 12)
        smooth_w = np.convolve(np.pad(shape_score, win // 2, mode="edge"), np.ones(win) / win, mode="valid")[:n]
        lo, hi = int(n * 0.15), int(n * 0.85)
        if hi > lo:
            k = lo + int(np.argmin(smooth_w[lo:hi]))
            cen = k
            unclear = smooth_w[k] >= 0.88 * float(np.median(smooth_w[lo:hi]))
    if cen >= 0 and cen > n - 1 - cen:  # p (the shorter arm) first
        line = line[::-1]
        widths = widths[::-1]
        cen = n - 1 - cen
    out.midline = [(float(x + x0), float(y + y0)) for x, y in line]
    out.widths = [round(float(v), 2) for v in widths]
    out.centromere = cen
    hw = float(np.median(widths)) / 2 if len(widths) else 2.0
    if len(path) > 8 and branch_count(skel, 1.5 * max(hw, 1.5)) > 2:
        out.checks.append(CHECK_BRANCHED)
    if unclear:
        out.checks.append(CHECK_CENTROMERE)
    if image_shape is not None:
        ys, xs = np.nonzero(m)
        if len(xs) and (ys.min() + y0 <= 0 or xs.min() + x0 <= 0 or ys.max() + y0 >= image_shape[0] - 1 or xs.max() + x0 >= image_shape[1] - 1):
            out.checks.append(CHECK_EDGE)
    return out


# ---------------------------------------------------------------------------
# detection
# ---------------------------------------------------------------------------


def prepared(grey: np.ndarray) -> np.ndarray:
    """Normalised counterstain with slow background changes removed, lightly smoothed."""
    g = normalise(grey)
    if ndi is None:
        return g
    background = ndi.gaussian_filter(ndi.grey_opening(g[::4, ::4], size=(31, 31)), sigma=8)
    background = np.kron(background, np.ones((4, 4)))[: g.shape[0], : g.shape[1]]
    g = np.clip(g - np.minimum(background, np.percentile(g, 60)), 0.0, None)
    return ndi.gaussian_filter(g, sigma=1.0)


def foreground(grey: np.ndarray, sensitivity: float = 0.5, *, ready: np.ndarray | None = None) -> np.ndarray:
    """Bool mask of chromosome pixels: the prepared counterstain above an Otsu threshold.

    Sensitivity 0.5 uses Otsu's threshold; higher values lower it (fainter chromosomes), lower values raise it
    (tighter outlines, touching chromosomes come apart more easily).
    """
    g = prepared(grey) if ready is None else ready
    sample = g[:: max(1, g.shape[0] // 500), :: max(1, g.shape[1] // 500)]
    s = max(0.0, min(1.0, sensitivity))
    thr = otsu(sample) * (1.45 - 0.9 * s)
    mask = g > thr
    if ndi is not None:
        # faint chromosomes far from the bright ones: brighter than their own surroundings
        around = ndi.gaussian_filter(g, sigma=15)
        mask |= (g > around * (1.6 - 0.5 * s)) & (g > thr * (0.45 - 0.2 * s))
        mask = ndi.binary_opening(mask, structure=np.ones((2, 2), dtype=bool))
        holes = ndi.binary_fill_holes(mask) & ~mask
        lab, n = ndi.label(holes)
        if n:
            sizes = np.bincount(lab.ravel())
            small = sizes < 40
            small[0] = False
            mask |= small[lab]  # fill only small holes; gaps between chromosomes stay open
    return mask


def half_width(mask: np.ndarray, dist: np.ndarray | None = None) -> float:
    """Typical half-width of the chromosomes: distance to the edge along the skeleton."""
    if dist is None:
        dist = ndi.distance_transform_edt(mask)
    ridge = thin(mask)
    return max(1.5, float(np.percentile(dist[ridge], 60)) if ridge.any() else 2.0)


def _segments(mask: np.ndarray, g: np.ndarray, separate: float, hw: float) -> tuple[np.ndarray, float]:
    """Over-segmentation and the chromosome half-width.

    Cores are pixels brighter than their surroundings (local contrast shows the dark seams between touching
    chromosomes), opened to break thin bridges; each core then grows back over the mask by a watershed on
    the intensity.
    """
    lab0, n0 = ndi.label(mask)
    if n0 == 0:
        return lab0, hw
    sep = max(0.0, min(1.0, separate))
    local = ndi.gaussian_filter(g, sigma=float(np.clip(hw, 4.0, 12.0)))
    core = mask & (g > local * (0.92 + 0.26 * sep))
    hwc = half_width(core) if core.any() else hw
    core = ndi.binary_opening(core, iterations=max(1, int(round((0.2 + 0.3 * sep) * hwc))))
    markers, nc = ndi.label(core)
    if nc:
        sizes = np.bincount(markers.ravel())
        small = sizes < max(6.0, hwc * hwc)
        small[0] = False
        markers[small[markers]] = 0
        markers, nc = ndi.label(markers > 0)
    if nc == 0:
        return lab0, hwc
    markers = markers.astype(np.int32)
    markers[~mask] = nc + 1
    cost = g.max() - g
    cost = (cost / (cost.max() or 1.0) * 65000).astype(np.uint16)
    res = ndi.watershed_ift(cost, markers)
    res[~mask] = 0
    res[res == nc + 1] = 0
    missing = mask & (res == 0)
    if missing.any():
        lm, _nm = ndi.label(missing)
        res[missing] = lm[missing] + nc + 1
    _labels, inverse = np.unique(res, return_inverse=True)
    return inverse.reshape(res.shape).astype(np.int32), hwc


def _neighbours(skel: np.ndarray) -> np.ndarray:
    pad = np.pad(skel.astype(np.uint8), 1)
    total = np.zeros(skel.shape, dtype=np.int16)
    for dr, dc in _NEIGH:
        total += pad[1 + dr:1 + dr + skel.shape[0], 1 + dc:1 + dc + skel.shape[1]]
    return total


def branch_count(skel: np.ndarray, spur: float) -> int:
    """Ends of the skeleton after removing spurs shorter than ``spur`` pixels (2 = one simple line)."""
    s = np.asarray(skel, dtype=bool).copy()
    for _ in range(max(1, int(spur))):
        ends = s & (_neighbours(s) <= 1)
        if not ends.any() or s.sum() - ends.sum() < 2:
            break
        s &= ~ends
    return int((s & (_neighbours(s) == 1)).sum())


def _unbranched(both: np.ndarray, hw: float) -> bool:
    crop, _x0, _y0 = _crop(both)
    return branch_count(thin(crop), 1.5 * hw) <= 2


def _merge_pieces(labels: np.ndarray, g: np.ndarray, separate: float, hw: float) -> np.ndarray:
    """Join neighbouring pieces when the seam between them is nearly as bright as the pieces (no dark gap)
    and together they still form one unbranched chromosome; the brightest seams first."""
    n = int(labels.max())
    if n < 2:
        return labels
    need = 0.88 + 0.1 * max(0.0, min(1.0, separate))
    flat_lab = labels.ravel()
    order = np.argsort(flat_lab, kind="stable")
    bounds = np.searchsorted(flat_lab[order], np.arange(n + 2))
    gv = g.ravel()
    level = np.ones(n + 1)
    for k in range(1, n + 1):
        vals = gv[order[bounds[k]:bounds[k + 1]]]
        level[k] = float(np.percentile(vals, 60)) if vals.size else 1.0
    grown = ndi.grey_dilation(labels, size=(3, 3))
    seam = (labels > 0) & (grown > 0) & (grown != labels)
    pairs: dict[tuple[int, int], list[float]] = {}
    for a, b, v in zip(labels[seam], grown[seam], g[seam]):
        key = (int(min(a, b)), int(max(a, b)))
        pairs.setdefault(key, []).append(float(v))
    candidates = []
    for (a, b), vals in pairs.items():
        if len(vals) >= 2:
            ratio = float(np.median(vals)) / max(1e-6, min(level[a], level[b]))
            if ratio >= need:
                candidates.append((ratio, a, b))
    areas = np.bincount(flat_lab, minlength=n + 1).astype(np.float64)
    cap = 1.8 * float(np.median(areas[1:][areas[1:] > 0])) if n > 2 else float("inf")
    parent = list(range(n + 1))

    def root(k: int) -> int:
        while parent[k] != k:
            parent[k] = parent[parent[k]]
            k = parent[k]
        return k

    slices = ndi.find_objects(labels)
    for _ratio, a, b in sorted(candidates, reverse=True):
        ra, rb = root(a), root(b)
        if ra == rb:
            continue
        members = [k for k in range(1, n + 1) if root(k) in (ra, rb)]
        if areas[members].sum() > cap:
            continue
        ys = slice(min(slices[k - 1][0].start for k in members), max(slices[k - 1][0].stop for k in members))
        xs = slice(min(slices[k - 1][1].start for k in members), max(slices[k - 1][1].stop for k in members))
        both = np.isin(labels[ys, xs], members)
        if _unbranched(both, hw):
            parent[rb] = ra
    lut = np.array([root(k) for k in range(n + 1)], dtype=np.int32)
    merged = lut[labels]
    _labels, inverse = np.unique(merged, return_inverse=True)
    return inverse.reshape(labels.shape).astype(np.int32)


DETECT_SIDE = 1400  # larger images are searched on a reduced copy (same result, seconds instead of minutes)


def _reduce(g: np.ndarray, step: int) -> np.ndarray:
    h, w = (g.shape[0] // step) * step, (g.shape[1] // step) * step
    return g[:h, :w].reshape(h // step, step, w // step, step).mean(axis=(1, 3))


def _enlarge(obj: Measured, step: int, shape: tuple[int, int]) -> Measured:
    mask = np.repeat(np.repeat(obj.mask, step, axis=0), step, axis=1)
    x0, y0 = obj.x0 * step, obj.y0 * step
    mask = mask[: max(0, shape[0] - y0), : max(0, shape[1] - x0)]
    off = (step - 1) / 2.0
    return Measured(mask=mask, x0=x0, y0=y0,
                    midline=[(x * step + off, y * step + off) for x, y in obj.midline],
                    centromere=obj.centromere, widths=[w * step for w in obj.widths], checks=list(obj.checks))


def detect_chromosomes(grey: np.ndarray, settings: DetectSettings | None = None) -> list[Measured]:
    """All chromosomes of a counterstain image, measured, sorted longest first."""
    if ndi is None:
        raise RuntimeError("Chromosome detection needs scipy (pip install scipy)")
    settings = settings or DetectSettings()
    g0 = np.asarray(grey, dtype=np.float64)
    if g0.ndim == 3:
        g0 = g0.max(axis=2)
    step = max(1, int(math.ceil(max(g0.shape) / DETECT_SIDE)))
    if step == 1:
        return _detect(g0, settings)
    small = DetectSettings(settings.sensitivity, settings.separate,
                           int(round(settings.smallest_px / (step * step))) if settings.smallest_px else 0)
    return [_enlarge(o, step, g0.shape) for o in _detect(_reduce(g0, step), small)]


def _detect(g0: np.ndarray, settings: DetectSettings) -> list[Measured]:
    g = prepared(g0)
    mask = foreground(g0, settings.sensitivity, ready=g)
    if not mask.any():
        return []
    hw = half_width(mask)
    labels, hw = _segments(mask, g, settings.separate, hw)
    labels = _merge_pieces(labels, g, settings.separate, hw)
    n = int(labels.max())
    areas = np.bincount(labels.ravel(), minlength=n + 1)
    smallest = settings.smallest_px
    if smallest <= 0:
        smallest = max(20, int(round(3.0 * hw * hw)))
    objects = []
    slices = ndi.find_objects(labels)
    for k, sl in enumerate(slices, start=1):
        if sl is None or areas[k] < smallest:
            continue
        ys, xs = sl
        y0, y1 = max(0, ys.start - 2), min(g.shape[0], ys.stop + 2)
        x0, x1 = max(0, xs.start - 2), min(g.shape[1], xs.stop + 2)
        crop = labels[y0:y1, x0:x1] == k
        objects.append(measure(crop, x0, y0, g.shape, g))
    if len(objects) >= 5 and settings.smallest_px <= 0:  # specks: far smaller than the big chromosomes
        big = float(np.percentile([o.area for o in objects], 90))
        objects = [o for o in objects if o.area >= 0.06 * big]
    flag_sizes(objects)
    objects.sort(key=lambda o: -o.length)
    return objects


def flag_sizes(objects: list[Measured]) -> None:
    """Add the large / small checks against the median area of the cell."""
    if len(objects) < 3:
        return
    med = float(np.median([o.area for o in objects]))
    for o in objects:
        o.checks = [c for c in o.checks if c not in (CHECK_LARGE, CHECK_SMALL)]
        if o.area > 1.9 * med:
            o.checks.append(CHECK_LARGE)
        elif o.area < 0.25 * med:
            o.checks.append(CHECK_SMALL)


# ---------------------------------------------------------------------------
# edits
# ---------------------------------------------------------------------------


def split_along(mask: np.ndarray, x0: int, y0: int, a: tuple[float, float], b: tuple[float, float], width: float = 1.5) -> list[tuple[np.ndarray, int, int]]:
    """Cut an object along the line a-b (image pixels); returns the pieces (crop, x0, y0), largest first.

    The line is lengthened a little at both ends so a cut drawn just across the chromosome goes all the way.
    """
    m = np.asarray(mask, dtype=bool).copy()
    h, w = m.shape
    ax, ay = a[0] - x0, a[1] - y0
    bx, by = b[0] - x0, b[1] - y0
    length = math.hypot(bx - ax, by - ay)
    if length < 1:
        return [(m, x0, y0)]
    dx, dy = (bx - ax) / length, (by - ay) / length
    ax, ay, bx, by = ax - dx * 4, ay - dy * 4, bx + dx * 4, by + dy * 4
    yy, xx = np.mgrid[0:h, 0:w]
    t = np.clip(((xx - ax) * dx + (yy - ay) * dy) / (length + 8), 0, 1) * (length + 8)
    px, py = ax + dx * t, ay + dy * t
    cut = np.hypot(xx - px, yy - py) <= width
    pieces_mask = m & ~cut
    lab, n = ndi.label(pieces_mask)
    if n < 2:
        return [(m, x0, y0)]
    # give the cut pixels back to the nearest piece
    _d, (ir, ic) = ndi.distance_transform_edt(lab == 0, return_indices=True)
    filled = lab[ir, ic]
    filled[~m] = 0
    sizes = np.bincount(filled.ravel())
    order = [k for k in np.argsort(-sizes) if k != 0 and sizes[k] > 0]
    out = []
    for k in order:
        piece, px0, py0 = _crop(filled == k, pad=1)
        out.append((piece, x0 + px0, y0 + py0))
    return out


def union(parts: list[tuple[np.ndarray, int, int]]) -> tuple[np.ndarray, int, int]:
    """One crop covering several (crop, x0, y0) masks."""
    x0 = min(p[1] for p in parts)
    y0 = min(p[2] for p in parts)
    x1 = max(p[1] + p[0].shape[1] for p in parts)
    y1 = max(p[2] + p[0].shape[0] for p in parts)
    out = np.zeros((y1 - y0, x1 - x0), dtype=bool)
    for crop, px, py in parts:
        out[py - y0:py - y0 + crop.shape[0], px - x0:px - x0 + crop.shape[1]] |= crop
    if ndi is not None:
        out = ndi.binary_closing(out, structure=np.ones((3, 3), dtype=bool))
    return out, x0, y0


def brush(mask: np.ndarray, x0: int, y0: int, cx: float, cy: float, radius: float, add: bool) -> tuple[np.ndarray, int, int]:
    """Paint (add) or erase a disc on an object; the crop grows when painting past its edge."""
    r = max(1.0, float(radius))
    m = np.asarray(mask, dtype=bool)
    if add:
        nx0 = min(x0, int(math.floor(cx - r)) - 1)
        ny0 = min(y0, int(math.floor(cy - r)) - 1)
        nx1 = max(x0 + m.shape[1], int(math.ceil(cx + r)) + 2)
        ny1 = max(y0 + m.shape[0], int(math.ceil(cy + r)) + 2)
        grown = np.zeros((ny1 - ny0, nx1 - nx0), dtype=bool)
        grown[y0 - ny0:y0 - ny0 + m.shape[0], x0 - nx0:x0 - nx0 + m.shape[1]] = m
        m, x0, y0 = grown, nx0, ny0
    else:
        m = m.copy()
    yy, xx = np.mgrid[0:m.shape[0], 0:m.shape[1]]
    disc = (xx + x0 - cx) ** 2 + (yy + y0 - cy) ** 2 <= r * r
    m[disc] = bool(add)
    return m, x0, y0


def object_at(grey: np.ndarray, x: float, y: float, sensitivity: float = 0.5, radius: int = 200) -> tuple[np.ndarray, int, int] | None:
    """The bright object under a click (for adding a chromosome the detection missed)."""
    g = np.asarray(grey)
    if g.ndim == 3:
        g = g.max(axis=2)
    h, w = g.shape
    xa, xb = max(0, int(x) - radius), min(w, int(x) + radius + 1)
    ya, yb = max(0, int(y) - radius), min(h, int(y) + radius + 1)
    if xb - xa < 5 or yb - ya < 5:
        return None
    win = foreground(g[ya:yb, xa:xb], sensitivity)
    lab, _n = ndi.label(win)
    cx, cy = int(x) - xa, int(y) - ya
    k = lab[cy, cx] if 0 <= cy < lab.shape[0] and 0 <= cx < lab.shape[1] else 0
    if not k:
        d, (ir, ic) = ndi.distance_transform_edt(lab == 0, return_indices=True)
        if d[cy, cx] > 6:
            return None
        k = lab[ir[cy, cx], ic[cy, cx]]
    crop, px, py = _crop(lab == k)
    return crop, xa + px, ya + py


def label_image(objects: list[tuple[np.ndarray, int, int]], shape: tuple[int, int]) -> np.ndarray:
    """int32 map: 0 background, i + 1 for object i (later objects win where they overlap)."""
    out = np.zeros(shape, dtype=np.int32)
    for i, (crop, x0, y0) in enumerate(objects):
        out[paste(shape, crop, x0, y0)] = i + 1
    return out


def outline(mask: np.ndarray) -> np.ndarray:
    """Bool edge pixels of a mask."""
    m = np.asarray(mask, dtype=bool)
    if ndi is None:
        return m
    return m & ~ndi.binary_erosion(m)
