"""Regions not suitable for hybridisation probes: low complexity, simple
repeats and inverted repeats (the dustmasker / einverted step of probe design).

Pure Python, no external tools. Coordinates are 0-based half-open.
"""

from __future__ import annotations

import re
from dataclasses import asdict, dataclass

_COMP = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def _revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclass
class Region:
    start: int          # 0-based
    end: int            # exclusive
    kind: str           # "low complexity" / "simple repeat" / "inverted repeat"
    note: str = ""

    def to_dict(self) -> dict:
        return asdict(self)

    def label(self) -> str:
        return f"{self.start + 1}-{self.end} {self.kind}" + (f" ({self.note})" if self.note else "")


def dust_regions(seq: str, *, window: int = 64, level: int = 30) -> list[Region]:
    """DUST low-complexity scan (triplet score of Morgulis et al. 2006).

    A window scores sum c(c-1)/2 / (l-1) over its l triplets; windows above
    ``level``/10 are low complexity, and in those the bases of triplets seen
    five or more times are masked. ``level`` 30 = ``dustmasker -level 30``.
    """
    s = seq.upper()
    n = len(s)
    if n < 3:
        return []
    w = min(window, n)
    trip = [s[i:i + 3] for i in range(n - 2)]
    lt = w - 2
    if lt < 2:
        return []
    thr = level / 10.0
    masked = bytearray(n)
    counts: dict[str, int] = {}
    pairs = 0
    for i in range(lt):
        c = counts.get(trip[i], 0)
        pairs += c
        counts[trip[i]] = c + 1
    for start in range(0, n - w + 1):
        if start:
            old = trip[start - 1]
            c = counts[old] - 1
            pairs -= c
            counts[old] = c
            new = trip[start + lt - 1]
            c = counts.get(new, 0)
            pairs += c
            counts[new] = c + 1
        if pairs / (lt - 1) > thr:
            hot = {t for t, c in counts.items() if c >= 5 and "N" not in t}
            for j in range(start, start + lt):
                if trip[j] in hot:
                    masked[j:j + 3] = b"\x01\x01\x01"
    # chance triplet repeats next to a real one leave short scraps: keep runs ≥ 10 nt
    return [Region(a, b, "low complexity") for a, b in _runs(masked) if b - a >= 10]


_SSR = re.compile(r"([ACGT]{1,6}?)\1{2,}")


def simple_repeat_regions(seq: str, *, min_length: int = 12) -> list[Region]:
    """Microsatellites: a 1 to 6 nt unit repeated at least three times, ≥ ``min_length`` nt."""
    out = []
    s = seq.upper()
    for m in _SSR.finditer(s):
        if m.end() - m.start() >= min_length:
            out.append(Region(m.start(), m.end(), "simple repeat", f"({m.group(1)})n"))
    return out


def inverted_repeat_regions(seq: str, *, min_stem: int = 12, max_span: int = 300,
                            min_loop: int = 0) -> list[Region]:
    """Inverted repeats (stems of a possible hairpin, like EMBOSS einverted).

    A ``min_stem``-mer whose reverse complement occurs again downstream, the
    two arms at most ``max_span`` nt apart (start to end), is a stem; both arms
    are reported. Perfect palindromes (the k-mer is its own reverse
    complement) count too.
    """
    s = seq.upper()
    n = len(s)
    k = min_stem
    if n < k:
        return []
    where: dict[str, list[int]] = {}
    for i in range(n - k + 1):
        kmer = s[i:i + k]
        if "N" not in kmer:
            where.setdefault(kmer, []).append(i)
    masked = bytearray(n)
    for i in range(n - k + 1):
        kmer = s[i:i + k]
        if "N" in kmer:
            continue
        rc = _revcomp(kmer)
        for j in where.get(rc, ()):
            if j == i or (j >= i + k + min_loop and j + k - i <= max_span):
                masked[i:i + k] = b"\x01" * k
                masked[j:j + k] = b"\x01" * k
    return [Region(a, b, "inverted repeat") for a, b in _runs(masked)]


def unsuitable_regions(seq: str, *, dust_level: int = 30, min_stem: int = 12,
                       max_span: int = 300) -> list[Region]:
    """All regions a probe should avoid, sorted, overlapping ones of one kind merged."""
    regs = (dust_regions(seq, level=dust_level) + simple_repeat_regions(seq)
            + inverted_repeat_regions(seq, min_stem=min_stem, max_span=max_span))
    regs.sort(key=lambda r: (r.start, r.end))
    return regs


def overlaps(start: int, end: int, regions: list[Region]) -> bool:
    """True when [start, end) touches any region."""
    return any(r.start < end and start < r.end for r in regions)


def _runs(mask: bytearray) -> list[tuple[int, int]]:
    out = []
    i, n = 0, len(mask)
    while i < n:
        if mask[i]:
            j = i
            while j < n and mask[j]:
                j += 1
            out.append((i, j))
            i = j
        else:
            i += 1
    return out
