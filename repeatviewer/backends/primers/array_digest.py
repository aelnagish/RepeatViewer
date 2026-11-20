"""Digest of a tandem array: monomer ladders, site loss, methylation and partial cutting.

Every site of the enzyme on the sequence is cut with probability

    p = completeness × (1 − m)

where ``m`` is the share of molecules methylated at that site in a context
that blocks the enzyme (the highest blocking cytosine with enough reads in the
user's file; 0 without data or without reads). Sites are independent, so a
fragment between sites i and j (sequence ends always "cut") has the molar
amount p_i · p_j · Π(1 − p_k) over the sites in between. Copies that lost the
site by mutation give the 2×, 3× … monomer ladder; methylated or uncut sites
give long fragments, up to the uncut array.

Site loss "From the copies" uses the sites as they are in the sequence.
"A site in every copy" adds the missing sites at the array period (same phase
as most sites), to show the ladder expected without mutation.

The period comes from monomer annotations (a feature name with three or more
copies of similar length) or, without them, from the most common distance
between repeated 12-mers.

Pure Python, no Qt; unit-tested in ``tests/test_primers_array_digest.py``.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

from collections import Counter
from dataclasses import dataclass, field

from repeatviewer.backends.primers import enzymes as ez

LOSS_COPIES = "copies"
LOSS_EVERY = "every"
LOSS_LABELS = {LOSS_COPIES: "From the copies", LOSS_EVERY: "A site in every copy"}

MIN_AMOUNT = 1e-4      # molar amounts below this are left out
MAX_BANDS = 80


@dataclass
class ArrayInfo:
    period: int = 0
    copies: int = 0
    source: str = ""          # "annotation", "estimated" or ""

    def label(self) -> str:
        if not self.period:
            return "no tandem array found"
        how = "from the annotation" if self.source == "annotation" else "estimated"
        return f"{self.copies} copies of {self.period:,} bp ({how})"


def estimate_period(seq: str, *, k: int = 12, lo: int = 20, hi: int = 5000) -> int:
    """Most common distance between consecutive copies of the same k-mer (0 if none)."""
    seq = seq.upper()
    last: dict[str, int] = {}
    dist: Counter = Counter()
    step = max(1, len(seq) // 200_000)
    for i in range(0, len(seq) - k + 1, step):
        kmer = seq[i:i + k]
        if "N" in kmer:
            continue
        j = last.get(kmer)
        if j is not None and lo <= i - j <= hi:
            dist[i - j] += 1
        last[kmer] = i
    if not dist:
        return 0
    best, n = dist.most_common(1)[0]
    # merge near values (indels shift the distance by a few bases)
    near = sum(c for d, c in dist.items() if abs(d - best) <= max(2, best // 50))
    total = sum(dist.values())
    if near < max(10, total * 0.2):
        return 0
    # a multiple of a shorter period? prefer the shortest that is common too
    for d in sorted(dist):
        if d < best and best % d == 0 and dist[d] >= n * 0.5:
            return d
    return best


def array_info(seq: str, features=()) -> ArrayInfo:
    """Period and copy number from monomer annotations, else estimated."""
    by_name: dict[str, list[int]] = {}
    for f in features or ():
        name = str(getattr(f, "name", "") or getattr(f, "label", "") or "")
        try:
            start, end = int(getattr(f, "start")), int(getattr(f, "end"))
        except (TypeError, ValueError, AttributeError):
            continue
        if name and end > start:
            by_name.setdefault(name, []).append(end - start + 1)
    best: list[int] = []
    for name, lengths in by_name.items():
        if len(lengths) < 3:
            continue
        mid = sorted(lengths)[len(lengths) // 2]
        similar = [x for x in lengths if abs(x - mid) <= mid * 0.15]
        monomerish = "monomer" in name.lower() or "repeat" in name.lower()
        if len(similar) >= 3 and (monomerish or len(similar) > len(best)):
            best = similar
            if monomerish:
                break
    if best:
        period = sorted(best)[len(best) // 2]
        return ArrayInfo(period=period, copies=len(best), source="annotation")
    period = estimate_period(seq)
    if not period:
        return ArrayInfo()
    return ArrayInfo(period=period, copies=max(1, round(len(seq) / period)), source="estimated")


@dataclass
class Band:
    size: int
    molar: float
    mass_share: float = 0.0
    label: str = ""


@dataclass
class DigestLane:
    label: str
    enzyme: str
    use_data: bool
    bands: list[Band] = field(default_factory=list)
    sites: int = 0
    blocked: int = 0
    added: int = 0            # sites added for "a site in every copy"
    why: str = ""


@dataclass
class DigestSettings:
    site_loss: str = LOSS_COPIES
    completeness: float = 0.95
    min_percent: float = 50.0
    min_reads: int = 5
    agarose: float = 1.2
    ladder: str = ""


def cut_probabilities(sites: list[ez.Site], completeness: float, min_reads: int,
                      use_data: bool) -> list[float]:
    out = []
    for s in sites:
        m = ez.blocked_fraction(s, min_reads) if use_data else None
        out.append(max(0.0, min(1.0, completeness * (1.0 - (m or 0.0)))))
    return out


def partial_fragments(length: int, cuts: list[int], probs: list[float], *,
                      min_amount: float = MIN_AMOUNT) -> dict[int, float]:
    """Fragment size → molar amount for independent partial cutting (ends always cut)."""
    pts = sorted(zip(cuts, probs))
    pos = [1] + [c for c, _p in pts if 1 < c <= length] + [length + 1]
    pr = [1.0] + [p for c, p in pts if 1 < c <= length] + [1.0]
    out: dict[int, float] = {}
    for i in range(len(pos) - 1):
        if pr[i] <= 0:
            continue
        run = pr[i]
        for j in range(i + 1, len(pos)):
            amount = run * pr[j]
            if amount >= min_amount:
                size = pos[j] - pos[i]
                out[size] = out.get(size, 0.0) + amount
            run *= (1.0 - pr[j])
            if run < min_amount:
                break
    return out


def _label(size: int, length: int, period: int) -> str:
    if size >= length * 0.98:
        return "uncut"
    if period:
        n = round(size / period)
        if n >= 1 and abs(size - n * period) <= max(3, period * 0.06):
            return f"{n}×"
    return ""


def _fill_every_copy(sites: list[ez.Site], length: int, period: int, name: str,
                     size: int) -> list[ez.Site]:
    if not period or len(sites) < 2:
        return list(sites)
    phase = Counter(s.cut % period for s in sites).most_common(1)[0][0]
    have = sorted(s.cut for s in sites)
    tol = max(3, period // 10)
    out = list(sites)
    import bisect

    c = phase if phase > 1 else phase + period
    lo, hi = min(have) - period, max(have) + period
    while c <= length:
        if lo <= c <= hi:
            i = bisect.bisect_left(have, c)
            near = any(abs(have[k] - c) <= tol for k in (i - 1, i) if 0 <= k < len(have))
            if not near:
                out.append(ez.Site(enzyme=name, start=c, end=min(length, c + size - 1), cut=c,
                                   reason="added: a site in every copy"))
        c += period
    return sorted(out, key=lambda s: s.cut)


def digest_lane(sequence: str, row: ez.EnzymeRow, info: ArrayInfo, s: DigestSettings, *,
                use_data: bool, label: str = "") -> DigestLane:
    n = len(sequence)
    sites = list(row.sites)
    real = len(sites)
    if s.site_loss == LOSS_EVERY:
        try:
            size = int(ez.enzyme(row.name).size)
        except ValueError:
            size = 4
        sites = _fill_every_copy(sites, n, info.period, row.name, size)
    probs = cut_probabilities(sites, s.completeness, s.min_reads, use_data)
    frags = partial_fragments(n, [x.cut for x in sites], probs)
    total = sum(size * m for size, m in frags.items()) or 1.0
    bands = [Band(size=size, molar=m, mass_share=size * m / total,
                  label=_label(size, n, info.period)) for size, m in frags.items()]
    bands.sort(key=lambda b: -b.size * b.molar)
    bands = sorted(bands[:MAX_BANDS], key=lambda b: b.size)
    blocked = sum(1 for x in sites[:] if use_data and x.blocked)
    lane = DigestLane(label=label or row.name, enzyme=row.name, use_data=use_data, bands=bands,
                      sites=real, blocked=blocked, added=len(sites) - real)
    lane.why = _why(lane, info, s)
    return lane


def _why(lane: DigestLane, info: ArrayInfo, s: DigestSettings) -> str:
    parts = []
    if lane.use_data and lane.blocked:
        parts.append(f"{lane.blocked} of {lane.sites} {lane.enzyme} sites are methylated in your "
                     f"data (≥ {s.min_percent:.0f} %, ≥ {s.min_reads} reads).")
    elif lane.use_data:
        parts.append(f"None of the {lane.sites} {lane.enzyme} sites is blocked in your data.")
    else:
        parts.append(f"{lane.sites} {lane.enzyme} sites, no methylation.")
    if info.copies and lane.sites < info.copies and lane.added == 0:
        parts.append(f"Only {lane.sites} of about {info.copies} copies carry the site: the "
                     "others lost it, which makes the 2×, 3× … ladder.")
    if lane.added:
        parts.append(f"{lane.added} sites added at the period (a site in every copy).")
    if s.completeness < 1.0:
        parts.append(f"Each site cut in {s.completeness * 100:.0f} % of molecules.")
    return " ".join(parts)


def digest_lanes(sequence: str, rows: list[ez.EnzymeRow], info: ArrayInfo, s: DigestSettings,
                 *, has_data: bool) -> list[DigestLane]:
    """One lane per enzyme (with data when there is data), plus a "no data" lane for enzymes
    with sites blocked in the data."""
    lanes = []
    extra = []
    for row in rows:
        lane = digest_lane(sequence, row, info, s, use_data=has_data,
                           label=f"{row.name} (your data)" if has_data and row.blocked else row.name)
        lanes.append(lane)
        if has_data and row.blocked:
            extra.append(digest_lane(sequence, row, info, s, use_data=False,
                                     label=f"{row.name} (no data)"))
    return lanes + extra


def gel_doc(lanes: list[DigestLane], length: int, s: DigestSettings, *, title: str = ""):
    from repeatviewer.backends.primers import gel_doc as gd

    sizes = [b.size for ln in lanes for b in ln.bands] or [length]
    ladder = s.ladder if s.ladder in gd.LADDERS else gd.pick_ladder_name(sizes)
    doc = gd.GelDoc(ladder=ladder, agarose=float(s.agarose), title=title)
    doc.lanes = [gd.ladder_lane(ladder)]
    for ln in lanes:
        doc.lanes.append(gd.GelLane(
            label=ln.label, bands=[gd.GelBand(size=b.size, amount=b.molar,
                                              label=b.label or f"{b.size} bp") for b in ln.bands],
            note=ln.why))
    return doc
