"""Signals step 3: straightened chromosomes, homologue pairing, profiles across cells, karyotype
from measured cells. No Qt.

Straightening samples the image along lines at right angles to the chromosome's midline
(as ImageJ's Straighten does). Pairing names the detected chromosomes: with a karyotype each
chromosome is matched to a karyotype chromosome by relative length, centromere index and
signal pattern (an assignment over the homologue slots); without one, chromosomes are paired
by the same measures and numbered from the longest pair down.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import segment as sm
from . import signals as sg
from .karyotype import SIGNAL_KIND_PREFIX, Band, Karyotype, KaryoChromosome, new_karyotype
from .model import UNNAMED, ChromosomeMark, ImageDetails

try:
    from scipy import ndimage as ndi
    from scipy.optimize import linear_sum_assignment
except ImportError:  # pragma: no cover
    ndi = None
    linear_sum_assignment = None

CHECK_PAIR = "homologues differ"
PROFILE_BINS = 100


# ---------------------------------------------------------------------------
# geometry
# ---------------------------------------------------------------------------


def _normals(line: np.ndarray) -> np.ndarray:
    d = np.gradient(line, axis=0)
    n = np.hypot(d[:, 0], d[:, 1])
    n[n == 0] = 1.0
    d = d / n[:, None]
    return np.column_stack([-d[:, 1], d[:, 0]])


def half_width(mark: ChromosomeMark) -> float:
    """Half the width of a detected chromosome (from its mask area and length)."""
    dec = sm.decode_mask(mark.mask) if mark.mask else None
    length = sm.polyline_length(mark.midline or [])
    if dec is None or length <= 0:
        return 6.0
    return max(2.0, float(dec[0].sum()) / length / 2.0)


def straighten(array: np.ndarray, midline, width: float) -> np.ndarray:
    """Image along the midline: rows from the p to the q telomere, ``2 * width + 1`` columns across."""
    a = np.asarray(array, dtype=np.float64)
    if a.ndim == 3:
        a = a[..., :3].max(axis=2)
    line = np.asarray(midline, dtype=np.float64)
    if len(line) < 2:
        return np.zeros((0, 0))
    normals = _normals(line)
    w = int(math.ceil(width))
    offsets = np.arange(-w, w + 1, dtype=np.float64)
    xs = line[:, 0][:, None] + normals[:, 0][:, None] * offsets[None, :]
    ys = line[:, 1][:, None] + normals[:, 1][:, None] * offsets[None, :]
    if ndi is None:
        xi = np.clip(np.round(xs).astype(int), 0, a.shape[1] - 1)
        yi = np.clip(np.round(ys).astype(int), 0, a.shape[0] - 1)
        return a[yi, xi]
    return ndi.map_coordinates(a, [ys, xs], order=1, mode="constant", cval=0.0)


def profile(array: np.ndarray, midline, width: float, bins: int = PROFILE_BINS) -> np.ndarray:
    """Mean intensity across the chromosome at ``bins`` evenly spaced places from p to q telomere."""
    strip = straighten(array, midline, width)
    if strip.size == 0:
        return np.zeros(bins)
    core = strip[:, max(0, strip.shape[1] // 4): strip.shape[1] - strip.shape[1] // 4 or None]
    values = core.mean(axis=1)
    t = np.linspace(0, 1, len(values))
    return np.interp(np.linspace(0, 1, bins), t, values)


def signal_fraction(mark: ChromosomeMark, x: float, y: float) -> float | None:
    """Position of a point along a chromosome, 0 = p telomere, 1 = q telomere."""
    if not mark.midline or len(mark.midline) < 2:
        return None
    total = sm.polyline_length(mark.midline)
    if total <= 0:
        return None
    return max(0.0, min(1.0, sg.arc_position(mark.midline, (x, y)) / total))


# ---------------------------------------------------------------------------
# measured chromosomes
# ---------------------------------------------------------------------------


@dataclass
class Measure:
    mark_id: str
    length: float  # pixels
    ci: float  # centromere index (percent, short arm / total)
    signals: dict[str, list[float]] = field(default_factory=dict)  # probe -> fractions from the p telomere


def measures(details: ImageDetails) -> list[Measure]:
    out = []
    for m in details.chromosomes:
        if not m.midline or len(m.midline) < 2:
            continue
        total = sm.polyline_length(m.midline)
        at = sg.arc_position(m.midline, m.centromere) if m.centromere else total / 2
        ci = 100.0 * min(at, total - at) / total if total > 0 else 50.0
        meas = Measure(m.id, total, ci)
        for s in details.signals:
            if s.chromosome_id == m.id and s.probe:
                f = signal_fraction(m, s.x, s.y)
                if f is not None:
                    meas.signals.setdefault(s.probe, []).append(f)
        out.append(meas)
    return out


def _pattern_distance(a: dict[str, list[float]], b: dict[str, list[float]]) -> float:
    probes = set(a) | set(b)
    if not probes:
        return 0.0
    cost = 0.0
    for p in probes:
        if (p in a) != (p in b):
            cost += 1.0
        else:
            cost += min(1.0, abs(float(np.median(a[p])) - float(np.median(b[p]))) * 3)
    return cost / len(probes)


def _pair_cost(a: Measure, b: Measure, total: float) -> float:
    la, lb = a.length / total, b.length / total
    return abs(la - lb) / max(la, lb, 1e-9) + abs(a.ci - b.ci) / 50.0 + 0.6 * _pattern_distance(a.signals, b.signals)


@dataclass
class Pairing:
    names: dict[str, tuple[str, str]]  # mark id -> (chromosome, homologue suffix)
    checks: dict[str, str]  # chromosome name -> why to look at it
    unplaced: list[str]  # mark ids left without a name


def pair_without_karyotype(details: ImageDetails, ploidy: int = 2) -> Pairing:
    """Group chromosomes into sets of ``ploidy`` similar ones and number the sets from the longest down."""
    ms = measures(details)
    total = sum(m.length for m in ms) or 1.0
    left = list(ms)
    groups: list[list[Measure]] = []
    k = max(1, int(ploidy))
    while len(left) >= k:
        left.sort(key=lambda m: -m.length)
        seed = left.pop(0)
        left.sort(key=lambda m: _pair_cost(seed, m, total))
        groups.append([seed] + left[: k - 1])
        del left[: k - 1]
    if left:
        groups.append(left)
    groups.sort(key=lambda g: -np.mean([m.length for m in g]))
    names, checks = {}, {}
    letters = "abcdefghijkl"
    for i, g in enumerate(groups, start=1):
        g.sort(key=lambda m: -m.length)
        for j, m in enumerate(g):
            names[m.mark_id] = (str(i), letters[j] if j < len(letters) else f".{j + 1}")
        lengths = [m.length for m in g]
        if len(g) < k or (max(lengths) - min(lengths)) / max(lengths) > 0.25:
            checks[str(i)] = CHECK_PAIR
    return Pairing(names, checks, [])


def _expected_fish(chrom: KaryoChromosome) -> dict[str, list[float]]:
    out: dict[str, list[float]] = {}
    for b in chrom.bands:
        if b.kind.startswith(SIGNAL_KIND_PREFIX):
            out.setdefault(b.kind[len(SIGNAL_KIND_PREFIX):], []).append((b.start + b.end) / 2)
    return out


def pair_with_karyotype(details: ImageDetails, karyotype: Karyotype) -> Pairing:
    """Match each detected chromosome to a homologue slot of the karyotype (Hungarian assignment)."""
    ms = measures(details)
    if not ms or not karyotype.chromosomes:
        return Pairing({}, {}, [m.mark_id for m in ms])
    total = sum(m.length for m in ms) or 1.0
    rel = karyotype.relative_lengths()
    haploid_total = sum(rel.values()) or 1.0
    slots: list[tuple[KaryoChromosome, str]] = []
    for c in karyotype.chromosomes:
        labels = karyotype.homologue_labels(c)
        for lab in labels:
            slots.append((c, lab[len(c.name):]))
    copies = max(1, len(slots) / max(1, len(karyotype.chromosomes)))
    has_lengths = any(c.length > 0 for c in karyotype.chromosomes)
    cost = np.zeros((len(ms), len(slots)))
    for i, m in enumerate(ms):
        frac = m.length / total * copies  # share of the haploid set
        for j, (c, _suffix) in enumerate(slots):
            exp = rel.get(c.name, 0.0) / haploid_total if has_lengths else 1.0 / len(karyotype.chromosomes)
            cl = abs(frac - exp) / max(exp, 1e-6)
            cc = abs(m.ci - c.centromere_index) / 50.0 if c.centromere_index > 0 else 0.0
            cp = 0.6 * _pattern_distance(m.signals, _expected_fish(c)) if _expected_fish(c) else 0.0
            cost[i, j] = cl + cc + cp
    if linear_sum_assignment is not None:
        rows, cols = linear_sum_assignment(cost)
    else:  # greedy fallback
        rows, cols, used = [], [], set()
        for i in np.argsort(cost.min(axis=1)):
            j = next(j for j in np.argsort(cost[i]) if j not in used)
            used.add(j)
            rows.append(i)
            cols.append(j)
    names, unplaced = {}, []
    by_chrom: dict[str, list[tuple[int, int]]] = {}
    placed = set()
    for i, j in zip(rows, cols):
        c, suffix = slots[j]
        by_chrom.setdefault(c.name, []).append((i, j))
        placed.add(i)
    # inside one chromosome the longest copy gets the first suffix
    for cname, items in by_chrom.items():
        suffixes = [slots[j][1] for _i, j in items]
        suffixes.sort(key=lambda s: [x for c, x in slots if c.name == cname].index(s))
        items.sort(key=lambda it: -ms[it[0]].length)
        for (i, _j), suffix in zip(items, suffixes):
            names[ms[i].mark_id] = (cname, suffix)
    checks = {}
    for cname, items in by_chrom.items():
        lengths = [ms[i].length for i, _j in items]
        worst = max(cost[i, j] for i, j in items)
        if (len(lengths) > 1 and (max(lengths) - min(lengths)) / max(lengths) > 0.25) or worst > 0.9:
            checks[cname] = CHECK_PAIR
    unplaced = [ms[i].mark_id for i in range(len(ms)) if i not in placed]
    return Pairing(names, checks, unplaced)


def apply_pairing(details: ImageDetails, pairing: Pairing) -> None:
    for m in details.chromosomes:
        if m.id in pairing.names:
            m.chromosome, m.homologue = pairing.names[m.id]
        elif m.id in pairing.unplaced:
            m.chromosome, m.homologue = UNNAMED, m.homologue if not m.named else "0"


def swap_names(details: ImageDetails, a_id: str, b_id: str) -> None:
    a, b = details.chromosome(a_id), details.chromosome(b_id)
    if a is not None and b is not None:
        (a.chromosome, a.homologue), (b.chromosome, b.homologue) = (b.chromosome, b.homologue), (a.chromosome, a.homologue)


# ---------------------------------------------------------------------------
# across cells
# ---------------------------------------------------------------------------


@dataclass
class ChromStats:
    name: str
    cells: int  # cells where this chromosome is named
    homologues: int  # named homologues over those cells
    length_px: float  # mean
    length_um: float | None
    rel_length: float  # percent of the haploid set (mean over cells)
    ci: float
    seen: dict[str, int] = field(default_factory=dict)  # probe -> cells with the signal on any homologue
    positions: dict[str, list[float]] = field(default_factory=dict)  # probe -> fractions from p telomere


def stats_over_cells(cells: list[ImageDetails], order: list[str] | None = None) -> list[ChromStats]:
    """Per chromosome name: mean length, CI, and in how many cells each probe was seen on it."""
    acc: dict[str, dict] = {}
    for d in cells:
        named = [m for m in d.chromosomes if m.named and m.midline and len(m.midline) >= 2]
        if not named:
            continue
        lengths = {m.id: sm.polyline_length(m.midline) for m in named}
        per_chrom: dict[str, list[ChromosomeMark]] = {}
        for m in named:
            per_chrom.setdefault(m.chromosome, []).append(m)
        haploid = sum(max(lengths[m.id] for m in ms) for ms in per_chrom.values()) or 1.0
        for cname, ms in per_chrom.items():
            a = acc.setdefault(cname, {"cells": 0, "hom": 0, "len": [], "um": [], "rel": [], "ci": [], "seen": {}, "pos": {}})
            a["cells"] += 1
            a["hom"] += len(ms)
            for m in ms:
                L = lengths[m.id]
                a["len"].append(L)
                if d.pixel_size_um:
                    a["um"].append(L * d.pixel_size_um)
                at = sg.arc_position(m.midline, m.centromere) if m.centromere else L / 2
                a["ci"].append(100.0 * min(at, L - at) / L if L > 0 else 50.0)
            a["rel"].append(100.0 * float(np.mean([lengths[m.id] for m in ms])) / haploid)
            ids = {m.id: m for m in ms}
            probes_here = set()
            for s in d.signals:
                if s.chromosome_id in ids and s.probe:
                    probes_here.add(s.probe)
                    f = signal_fraction(ids[s.chromosome_id], s.x, s.y)
                    if f is not None:
                        a["pos"].setdefault(s.probe, []).append(f)
            for p in probes_here:
                a["seen"][p] = a["seen"].get(p, 0) + 1
    rank = {n: i for i, n in enumerate(order or [])}
    out = []
    for cname, a in acc.items():
        out.append(ChromStats(cname, a["cells"], a["hom"], float(np.mean(a["len"])),
                              float(np.mean(a["um"])) if a["um"] else None, float(np.mean(a["rel"])),
                              float(np.mean(a["ci"])), dict(a["seen"]), dict(a["pos"])))
    from .importer import natural_key

    out.sort(key=lambda s: (rank.get(s.name, len(rank)), natural_key(s.name)))
    return out


def karyotype_from_cells(stats: list[ChromStats], name: str, *, species: str = "", ploidy: int = 2,
                         seen_at_least: float = 0.5, colours: dict[str, str] | None = None,
                         base: Karyotype | None = None) -> Karyotype:
    """A karyotype from measured cells: mean lengths and centromere indices, and a FISH band for every probe
    seen in at least ``seen_at_least`` of the cells, at its median position (± 0.015).

    With ``base`` its chromosomes keep their other bands and settings and only lengths, centromere and FISH
    bands are replaced.
    """
    colours = colours or {}
    if base is not None:
        k = Karyotype.from_dict(base.to_dict())
        k.name = name or base.name
    else:
        k = new_karyotype(name, max(1, len(stats)), ploidy, species)
        k.chromosomes = []
    use_um = all(s.length_um is not None for s in stats) and stats
    k.unit = "µm" if use_um else "%"
    for s in stats:
        chrom = k.chromosome(s.name)
        if chrom is None:
            chrom = KaryoChromosome(name=s.name)
            k.chromosomes.append(chrom)
        chrom.length = round(s.length_um, 3) if use_um else round(s.rel_length, 3)
        chrom.centromere_index = round(s.ci, 2)
        chrom.bands = [b for b in chrom.bands if not b.kind.startswith(SIGNAL_KIND_PREFIX)]
        for probe, n in s.seen.items():
            if s.cells and n / s.cells >= seen_at_least and s.positions.get(probe):
                f = float(np.median(s.positions[probe]))
                chrom.bands.append(Band(f"{SIGNAL_KIND_PREFIX}{probe}", max(0.0, f - 0.015), min(1.0, f + 0.015),
                                        colours.get(probe, "")))
    if base is None:
        k.basic_number = len([c for c in k.chromosomes if not c.b_chromosome]) or 1
    return k
