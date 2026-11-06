"""Signals on marked chromosomes: arm, position along the arm, summaries."""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass

from .importer import natural_key
from .model import ChromosomeMark, ImageDetails, Signal

CENTROMERE_ZONE = 0.12  # signals closer to the centromere than this share of the arm count as "cen"


def _sub(a, b):
    return (a[0] - b[0], a[1] - b[1])


def _dot(a, b):
    return a[0] * b[0] + a[1] * b[1]


def _seg_distance(p, a, b) -> float:
    ab = _sub(b, a)
    denom = _dot(ab, ab)
    t = 0.0 if denom == 0 else max(0.0, min(1.0, _dot(_sub(p, a), ab) / denom))
    q = (a[0] + t * ab[0], a[1] + t * ab[1])
    return math.hypot(p[0] - q[0], p[1] - q[1])


def _project(p, line) -> tuple[float, float, float]:
    """(distance to the polyline, arc length of the nearest point, total length)."""
    best_d, best_s, run = float("inf"), 0.0, 0.0
    for a, b in zip(line, line[1:]):
        ab = _sub(b, a)
        seg = math.hypot(*ab)
        denom = seg * seg
        t = 0.0 if denom == 0 else max(0.0, min(1.0, _dot(_sub(p, a), ab) / denom))
        q = (a[0] + t * ab[0], a[1] + t * ab[1])
        d = math.hypot(p[0] - q[0], p[1] - q[1])
        if d < best_d:
            best_d, best_s = d, run + t * seg
        run += seg
    return best_d, best_s, run


def arc_position(line, point) -> float:
    """Arc length along a polyline of the point nearest to ``point``."""
    return _project(point, line)[1]


def distance_to_mark(x: float, y: float, mark: ChromosomeMark) -> float:
    """Distance from a point to the chromosome axis (midline, or telomere-centromere-telomere), or its label."""
    p = (x, y)
    if mark.midline and len(mark.midline) >= 2:
        return _project(p, mark.midline)[0]
    pts = [q for q in (mark.telomere_p, mark.centromere, mark.telomere_q) if q is not None]
    if len(pts) >= 2:
        return min(_seg_distance(p, pts[i], pts[i + 1]) for i in range(len(pts) - 1))
    anchor = pts[0] if pts else (mark.x, mark.y)
    return math.hypot(p[0] - anchor[0], p[1] - anchor[1])


def nearest_chromosome(x: float, y: float, marks: list[ChromosomeMark], max_distance: float | None = None) -> ChromosomeMark | None:
    best = min(marks, key=lambda m: distance_to_mark(x, y, m), default=None)
    if best is None or (max_distance is not None and distance_to_mark(x, y, best) > max_distance):
        return None
    return best


def position_class(arm: str, rel: float | None) -> str:
    if arm == "cen" or rel is None or rel < CENTROMERE_ZONE:
        return "pericentromeric"
    if rel < 0.8:
        return "interstitial"
    if rel < 0.95:
        return "subterminal"
    return "terminal"


def classify(signal: Signal, mark: ChromosomeMark | None) -> None:
    """Set arm, rel_distance and position from the chromosome's centromere and telomeres.

    Signals set by hand (``manual``) are left alone; without a centromere nothing is set.
    """
    if signal.manual or mark is None or mark.centromere is None:
        return
    if mark.midline and len(mark.midline) >= 2:
        _d, at, total = _project((signal.x, signal.y), mark.midline)
        cen = arc_position(mark.midline, mark.centromere)
        if total <= 0:
            return
        if at <= cen:
            arm, rel = "p", (cen - at) / cen if cen > 0 else 0.0
        else:
            arm, rel = "q", (at - cen) / (total - cen) if total > cen else 0.0
        rel = max(0.0, min(1.0, rel))
        signal.arm = "cen" if rel < CENTROMERE_ZONE else arm
        signal.rel_distance = round(rel, 3)
        signal.position = position_class(signal.arm, rel)
        return
    c = mark.centromere
    s = (signal.x, signal.y)
    best_arm, best_t = "", -1.0
    for arm, tel in (("p", mark.telomere_p), ("q", mark.telomere_q)):
        if tel is None:
            continue
        v = _sub(tel, c)
        denom = _dot(v, v)
        if denom == 0:
            continue
        t = _dot(_sub(s, c), v) / denom
        if t > best_t:
            best_arm, best_t = arm, t
    if not best_arm:
        return
    rel = max(0.0, min(1.0, best_t))
    signal.arm = "cen" if rel < CENTROMERE_ZONE else best_arm
    signal.rel_distance = round(rel, 3)
    signal.position = position_class(signal.arm, rel)


def assign_and_classify(details: ImageDetails, signal: Signal, max_distance: float | None = None) -> None:
    """Attach a new signal to the nearest numbered chromosome and classify it."""
    if not signal.chromosome_id:
        mark = nearest_chromosome(signal.x, signal.y, details.chromosomes, max_distance)
        signal.chromosome_id = mark.id if mark else ""
    classify(signal, details.chromosome(signal.chromosome_id))


@dataclass
class SummaryRow:
    chromosome: str
    probe: str
    arm: str
    position: str
    seen: int  # homologues with the signal
    total: int  # homologues marked
    cells: int  # cells contributing


def summarize(images: list[ImageDetails], order: list[str] | None = None) -> list[SummaryRow]:
    """Per chromosome and probe: how many marked homologues carry the signal, over all cells."""
    total: Counter = Counter()
    cells_per_chrom: dict[str, set[int]] = {}
    seen: dict[tuple[str, str], set[tuple[int, str]]] = {}
    places: dict[tuple[str, str], Counter] = {}
    for idx, details in enumerate(images):
        marks = {m.id: m for m in details.chromosomes if m.named}
        for m in marks.values():
            total[m.chromosome] += 1
            cells_per_chrom.setdefault(m.chromosome, set()).add(idx)
        for s in details.signals:
            mark = marks.get(s.chromosome_id)
            if mark is None or not s.probe:
                continue
            key = (mark.chromosome, s.probe)
            seen.setdefault(key, set()).add((idx, mark.id))
            places.setdefault(key, Counter())[(s.arm, s.position)] += 1
    rank = {name: i for i, name in enumerate(order or [])}
    rows = []
    for (chrom, probe), hits in seen.items():
        (arm, position), _n = places[(chrom, probe)].most_common(1)[0]
        rows.append(SummaryRow(chrom, probe, arm, position, len(hits), total[chrom], len(cells_per_chrom.get(chrom, ()))))
    rows.sort(key=lambda r: (rank.get(r.chromosome, len(rank)), natural_key(r.chromosome), r.probe.lower()))
    return rows


def unsignalled_chromosomes(images: list[ImageDetails], rows: list[SummaryRow]) -> list[str]:
    marked = {m.chromosome for d in images for m in d.chromosomes if m.named}
    with_signal = {r.chromosome for r in rows}
    return sorted(marked - with_signal, key=natural_key)


def summary_tsv(rows: list[SummaryRow]) -> str:
    lines = ["chromosome\tprobe\tarm\tposition\thomologues_with_signal\thomologues_marked\tcells"]
    for r in rows:
        lines.append(f"{r.chromosome}\t{r.probe}\t{r.arm}\t{r.position}\t{r.seen}\t{r.total}\t{r.cells}")
    return "\n".join(lines) + "\n"
