"""Satellite monomer phase rotation (Tools → Tree, Alignment tab).

Monomers cut from a tandem array start at arbitrary positions: two copies of
the same monomer can be circular permutations of each other (``ACGTTG`` and
``TTGACG``). Aligned as they are, the shared sequence ends up split across
both ends of the alignment and the tree reflects where the cut fell, not how
the monomers are related. Rotating every monomer to a common start (and the
same strand) before alignment removes that artefact.

Method: one monomer of typical length (closest to the median) is the
reference. For every other monomer, and for its reverse complement, each
k-mer shared with the reference votes for the rotation offset that would put
the two in register (diagonal voting, as in dotplots). The strand and offset
with the most votes win. Monomers with too few shared k-mers are left
unchanged and reported.

Pure Python, no Qt; unit-tested in ``tests/test_monomer_phase.py``.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import math
import statistics
from collections import Counter
from dataclasses import dataclass, field

_COMPLEMENT = str.maketrans("ACGTUNRYSWKMBDHVacgtunryswkmbdhv-.",
                            "TGCAANYRSWMKVHDBtgcaanyrswmkvhdb-.")
DEFAULT_K = 8
SHORT_K = 6               # for monomers shorter than 60 bp
FALLBACK_KS = (6, 5)      # tried in turn for divergent monomers
MIN_VOTES = 3             # fewer shared k-mers than this: leave the monomer as it is
P_RANDOM = 0.001          # chance that an unrelated sequence passes (all offsets, both strands)


def reverse_complement(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


@dataclass
class RotationRow:
    name: str
    strand: str      # "+" as given, "-" reverse-complemented
    offset: int      # bases moved from the start to the end (0 = unchanged)
    votes: int       # shared k-mers supporting the chosen register
    matched: bool


@dataclass
class RotationResult:
    sequences: list[tuple[str, str]]
    reference: str
    k: int
    rows: list[RotationRow] = field(default_factory=list)

    @property
    def n_rotated(self) -> int:
        return sum(1 for r in self.rows if r.matched and r.offset)

    @property
    def n_reversed(self) -> int:
        return sum(1 for r in self.rows if r.matched and r.strand == "-")

    @property
    def unmatched(self) -> list[str]:
        return [r.name for r in self.rows if not r.matched]

    def table_tsv(self) -> str:
        lines = ["name\tstrand\toffset\tshared_kmers\tmatched"]
        for r in self.rows:
            lines.append(f"{r.name}\t{r.strand}\t{r.offset}\t{r.votes}\t{'yes' if r.matched else 'no'}")
        return "\n".join(lines) + "\n"


def choose_reference(seqs: list[str]) -> int:
    """Index of the monomer whose length is closest to the median (first on ties)."""
    lengths = [len(s) for s in seqs]
    median = statistics.median(lengths)
    return min(range(len(seqs)), key=lambda i: (abs(lengths[i] - median), i))


def pick_k(seqs: list[str]) -> int:
    return DEFAULT_K if min(len(s) for s in seqs) >= 60 else SHORT_K


def _poisson_sf(v: int, lam: float) -> float:
    """P(X >= v) for X ~ Poisson(lam)."""
    if v <= 0:
        return 1.0
    term = math.exp(-lam)
    cdf = term
    for i in range(1, v):
        term *= lam / i
        cdf += term
    return max(0.0, 1.0 - cdf)


def vote_threshold(seq_len: int, ref_len: int, k: int, min_votes: int = MIN_VOTES) -> int:
    """Pooled votes a register needs to count as real for this ``k``.

    For an unrelated sequence, chance k-mer hits (about ``seq_len * ref_len /
    4**k``) are spread over ``seq_len`` offsets, and the pooling window collects
    ``2w + 1`` of them: a Poisson count per window. The threshold is the
    smallest count that the best of all ``2 * seq_len`` windows (every offset,
    both strands) reaches by chance with probability below ``P_RANDOM``.
    """
    window = max(2, seq_len // 50)
    lam = ref_len * (2 * window + 1) / 4 ** k
    tests = 2 * max(1, seq_len)
    v = max(1, min_votes)
    while tests * _poisson_sf(v, lam) >= P_RANDOM:
        v += 1
    return v


def _kmer_positions(ref: str, k: int) -> dict[str, list[int]]:
    """k-mers of the reference read circularly (so the join also counts)."""
    circ = ref + ref[: k - 1]
    index: dict[str, list[int]] = {}
    for j in range(len(ref)):
        kmer = circ[j : j + k]
        if "N" not in kmer:
            index.setdefault(kmer, []).append(j)
    return index


def best_offset(seq: str, ref_index: dict[str, list[int]], k: int) -> tuple[int, int]:
    """(offset, votes): rotating ``seq`` left by ``offset`` puts it in register.

    Votes within a small window are pooled, so monomers that differ by a few
    indels (slightly shifted diagonals) still agree on one register. Inside
    that window the diagonal supported by k-mers near the reference's start
    wins, because the start is what has to line up.
    """
    n = len(seq)
    if n < k:
        return 0, 0
    circ = seq + seq[: k - 1]
    votes: Counter[int] = Counter()
    start_votes: Counter[int] = Counter()
    ref_len = max((max(p) for p in ref_index.values()), default=0) + 1
    start_zone = max(4 * k, ref_len // 4)
    # Count only non-overlapping k-mer hits per diagonal: one chance match of
    # k + 2 bases would otherwise cast three votes, and the random background
    # would no longer follow the Poisson model that sets the threshold.
    last_hit: dict[int, int] = {}
    for i in range(n):
        for j in ref_index.get(circ[i : i + k], ()):
            offset = (i - j) % n
            if i < last_hit.get(offset, -k) + k:
                continue
            last_hit[offset] = i
            votes[offset] += 1
            if j < start_zone:
                start_votes[offset] += 1
    if not votes:
        return 0, 0
    window = max(2, n // 50)

    def pooled(offset: int) -> int:
        return sum(votes.get((offset + d) % n, 0) for d in range(-window, window + 1))

    best = max(votes, key=lambda o: (pooled(o), votes[o]))
    score = pooled(best)
    near = [(best + d) % n for d in range(-window, window + 1)]
    if any(start_votes.get(o, 0) for o in near):
        best = max(near, key=lambda o: (start_votes.get(o, 0), votes.get(o, 0)))
    return best, score


def rotate(seq: str, offset: int) -> str:
    offset %= max(1, len(seq))
    return seq[offset:] + seq[:offset]


def rotate_to_common_start(
    records: list[tuple[str, str]], *, k: int | None = None, min_votes: int = MIN_VOTES
) -> RotationResult:
    """Rotate (and orient) every monomer to the register of a reference monomer.

    ``records`` are (name, ungapped sequence). Order is kept.
    """
    if not records:
        return RotationResult([], "", k or DEFAULT_K)
    seqs = [s.upper().replace("-", "").replace(".", "") for _n, s in records]
    k = k or pick_k(seqs)
    ref_i = choose_reference(seqs)
    ref = seqs[ref_i]
    ks = [k] + [kk for kk in FALLBACK_KS if kk < k]
    indexes = {kk: _kmer_positions(ref, kk) for kk in ks}
    out: list[tuple[str, str]] = []
    rows: list[RotationRow] = []
    for i, ((name, _orig), seq) in enumerate(zip(records, seqs)):
        if i == ref_i:
            out.append((name, seq))
            rows.append(RotationRow(name, "+", 0, len(seq), True))
            continue
        rc = reverse_complement(seq)
        chosen = None
        best_votes = 0
        for kk in ks:
            fwd_off, fwd_votes = best_offset(seq, indexes[kk], kk)
            rev_off, rev_votes = best_offset(rc, indexes[kk], kk)
            best_votes = max(best_votes, fwd_votes, rev_votes)
            need = vote_threshold(len(seq), len(ref), kk, min_votes)
            if max(fwd_votes, rev_votes) < need:
                continue
            if rev_votes > fwd_votes:
                chosen = ("-", rev_off, rev_votes)
            else:
                chosen = ("+", fwd_off, fwd_votes)
            break
        if chosen is None:
            out.append((name, seq))
            rows.append(RotationRow(name, "+", 0, best_votes, False))
            continue
        strand, offset, votes = chosen
        out.append((name, rotate(rc if strand == "-" else seq, offset)))
        rows.append(RotationRow(name, strand, offset, votes, True))
    return RotationResult(out, records[ref_i][0], k, rows)
