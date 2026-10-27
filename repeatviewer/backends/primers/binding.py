"""Where primers bind, and what they amplify (in-silico PCR).

Bases are stored as 4-bit masks (A=1, C=2, G=4, T=8; IUPAC codes are their
unions), so a primer base matches a template base when the masks share a bit.
Template N and gaps are 0 and never match, so runs of N cannot attract hits.

Two search strategies give the same answer where both apply:

* **Exhaustive** (inputs up to :data:`EXHAUSTIVE_MAX_BP`): every position is
  scored with numpy; mismatches are allowed anywhere, subject to the
  "last bases must match" rule.
* **Seeded** (genomes): the primer's 3' end (its last *k* bases, which must
  match anyway) is found with a regular expression, then each hit is scored.
  This is how a 3 Gb genome is searched in seconds per primer.

Coordinates in results are 1-based, inclusive, on the template's top strand.
"""

from __future__ import annotations

import gzip
import re
from collections.abc import Callable, Iterable, Iterator
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from repeatviewer.backends.primers.thermo import (
    IUPAC_BASES,
    PrimerError,
    binding_part,
    reverse_complement,
)

EXHAUSTIVE_MAX_BP = 3_000_000
MIN_SEED = 6
DEFAULT_SEED = 8
MAX_SITES_PER_PRIMER = 50_000
MAX_PRODUCTS = 20_000
CHUNK_BP = 4_000_000

_BIT = {"A": 1, "C": 2, "G": 4, "T": 8}
_MASK = np.zeros(256, dtype=np.uint8)
for _code, _bases in IUPAC_BASES.items():
    _m = 0
    for _b in _bases:
        _m |= _BIT[_b]
    _MASK[ord(_code)] = _m
    _MASK[ord(_code.lower())] = _m
_MASK[ord("N")] = 0     # template N / gaps never match
_MASK[ord("n")] = 0
_MASK[ord("U")] = _BIT["T"]
_MASK[ord("u")] = _BIT["T"]
# Primer N matches anything.
_PRIMER_MASK = _MASK.copy()
_PRIMER_MASK[ord("N")] = 15
_PRIMER_MASK[ord("n")] = 15


def encode_template(seq: str | bytes) -> np.ndarray:
    data = seq.encode("ascii", "replace") if isinstance(seq, str) else bytes(seq)
    return _MASK[np.frombuffer(data, dtype=np.uint8)]


def encode_primer(seq: str) -> np.ndarray:
    return _PRIMER_MASK[np.frombuffer(seq.upper().encode("ascii"), dtype=np.uint8)]


@dataclass
class SearchSettings:
    max_mismatches: int = 2
    # The last N bases at the primer's 3' end must match exactly (0 = no rule).
    three_prime_exact: int = 3
    both_strands: bool = True
    max_sites_per_primer: int = MAX_SITES_PER_PRIMER

    def __post_init__(self) -> None:
        self.max_mismatches = max(0, int(self.max_mismatches))
        self.three_prime_exact = max(0, int(self.three_prime_exact))

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict | None) -> SearchSettings:
        data = dict(data or {})
        return cls(**{k: data[k] for k in cls.__dataclass_fields__ if k in data})


@dataclass
class BindingSite:
    primer: str           # primer name
    record: str           # template record id
    strand: str           # "+" primer copies the top strand (points right), "-" points left
    start: int            # 1-based inclusive, top strand
    end: int
    mismatches: int
    mismatch_positions: list[int] = field(default_factory=list)  # 1-based in the primer, 5'→3'
    site: str = ""        # template stretch in the primer's orientation
    tm: float = 0.0       # Tm of the primer on this site (0 = not computed)

    @property
    def three_prime(self) -> int:
        """Template coordinate of the primer's 3' base."""
        return self.end if self.strand == "+" else self.start

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> BindingSite:
        return cls(**{k: d[k] for k in cls.__dataclass_fields__ if k in d})


@dataclass
class Product:
    record: str
    forward: str          # primer at the left end (binds +)
    reverse: str          # primer at the right end (binds -)
    start: int
    end: int
    mismatches: int
    forward_start: int = 0
    reverse_end: int = 0
    sequence: str = ""

    @property
    def size(self) -> int:
        return self.end - self.start + 1

    def to_dict(self) -> dict:
        d = asdict(self)
        d["size"] = self.size
        return d

    @classmethod
    def from_dict(cls, d: dict) -> Product:
        return cls(**{k: d[k] for k in cls.__dataclass_fields__ if k in d})


# --------------------------------------------------------------------------
# FASTA reading (streaming, fast enough for genomes)
# --------------------------------------------------------------------------


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    """(id, sequence) for each record; ids are the header up to the first space."""
    p = Path(path)
    opener = gzip.open if p.name.lower().endswith(".gz") else open
    name = None
    chunks: list[str] = []
    with opener(p, "rt", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name = line[1:].strip().split()[0] if line[1:].strip() else "unnamed"
                chunks = []
            elif name is not None:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def fasta_total_length(path: Path, limit: int | None = None) -> int:
    total = 0
    for _name, seq in iter_fasta(path):
        total += len(seq)
        if limit is not None and total > limit:
            return total
    return total


# --------------------------------------------------------------------------
# Matching one primer
# --------------------------------------------------------------------------


def _mismatch_rows(window: np.ndarray, pmask: np.ndarray) -> np.ndarray:
    """Boolean (n_windows, L): True where the template base does not fit."""
    return (window & pmask[None, :]) == 0


def _score_positions(tmask: np.ndarray, pmask: np.ndarray, settings: SearchSettings,
                     three_prime_at_end: bool) -> tuple[np.ndarray, np.ndarray]:
    """Exhaustive: start offsets whose window fits, and their mismatch counts."""
    L = len(pmask)
    n = len(tmask) - L + 1
    if n <= 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int32)
    mism = np.zeros(n, dtype=np.int16)
    exact_k = min(settings.three_prime_exact, L)
    if three_prime_at_end:
        exact_idx = range(L - exact_k, L)
    else:
        exact_idx = range(0, exact_k)
    ok = np.ones(n, dtype=bool)
    for i in exact_idx:
        ok &= (tmask[i:i + n] & pmask[i]) != 0
    if not ok.any():
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int32)
    for i in range(L):
        mism += ((tmask[i:i + n] & pmask[i]) == 0)
        if i % 4 == 3:
            ok &= mism <= settings.max_mismatches
    ok &= mism <= settings.max_mismatches
    idx = np.nonzero(ok)[0]
    return idx, mism[idx].astype(np.int32)


def _seed_regex(seed: str) -> re.Pattern:
    parts = []
    for ch in seed.upper():
        bases = IUPAC_BASES.get(ch, "ACGT")
        # Template letters (ambiguity codes too, never N) sharing a base with the primer.
        allowed = {code for code, members in IUPAC_BASES.items()
                   if code != "N" and set(members) & set(bases)}
        cls = "".join(sorted(allowed | {c.lower() for c in allowed}))
        parts.append(f"[{cls}]")
    return re.compile("(?=" + "".join(parts) + ")")


def _seeded_positions(seq: str, tmask: np.ndarray, oriented: str, settings: SearchSettings,
                      three_prime_at_end: bool) -> tuple[np.ndarray, np.ndarray]:
    L = len(oriented)
    k = max(MIN_SEED, settings.three_prime_exact, 1)
    k = min(k, L)
    if three_prime_at_end:
        seed, offset = oriented[L - k:], L - k
    else:
        seed, offset = oriented[:k], 0
    rx = _seed_regex(seed)
    starts = []
    cap = settings.max_sites_per_primer * 4
    for m in rx.finditer(seq):
        s = m.start() - offset
        if 0 <= s <= len(seq) - L:
            starts.append(s)
            if len(starts) >= cap:
                break
    if not starts:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int32)
    st = np.asarray(starts, dtype=np.int64)
    pmask = encode_primer(oriented)
    win = tmask[st[:, None] + np.arange(L)[None, :]]
    mm = _mismatch_rows(win, pmask)
    counts = mm.sum(axis=1)
    keep = counts <= settings.max_mismatches
    exact_k = min(settings.three_prime_exact, L)
    if exact_k:
        tail = mm[:, L - exact_k:] if three_prime_at_end else mm[:, :exact_k]
        keep &= ~tail.any(axis=1)
    return st[keep], counts[keep].astype(np.int32)


def find_sites_in_sequence(primer_name: str, primer_seq: str, record_id: str, seq: str,
                           settings: SearchSettings, *, exhaustive: bool | None = None,
                           offset: int = 0) -> list[BindingSite]:
    """Binding sites of one primer on one template sequence."""
    bind = binding_part(primer_seq)
    L = len(bind)
    if L == 0 or len(seq) < L:
        return []
    tmask = encode_template(seq)
    use_exhaustive = exhaustive if exhaustive is not None else len(seq) <= EXHAUSTIVE_MAX_BP
    sites: list[BindingSite] = []
    orientations = [("+", bind, True)]
    if settings.both_strands:
        # On the top strand a "-" primer reads as its reverse complement, whose
        # first base is the primer's 3' end.
        orientations.append(("-", reverse_complement(bind), False))
    for strand, oriented, at_end in orientations:
        if use_exhaustive:
            pmask = encode_primer(oriented)
            idx, mism = _score_positions(tmask, pmask, settings, at_end)
        else:
            idx, mism = _seeded_positions(seq, tmask, oriented, settings, at_end)
        if len(idx) > settings.max_sites_per_primer:
            idx, mism = idx[: settings.max_sites_per_primer], mism[: settings.max_sites_per_primer]
        if len(idx) == 0:
            continue
        pm = encode_primer(oriented)
        for s, m in zip(idx.tolist(), mism.tolist()):
            window = seq[s:s + L]
            positions: list[int] = []
            if m:
                bad = np.nonzero((tmask[s:s + L] & pm) == 0)[0].tolist()
                if strand == "+":
                    positions = [b + 1 for b in bad]
                else:
                    positions = sorted(L - b for b in bad)
            site = window.upper() if strand == "+" else reverse_complement(window.upper())
            sites.append(BindingSite(
                primer=primer_name, record=record_id, strand=strand,
                start=offset + s + 1, end=offset + s + L, mismatches=int(m),
                mismatch_positions=positions, site=site,
            ))
    sites.sort(key=lambda b: (b.start, b.strand))
    return sites


def find_sites(primers: list[tuple[str, str]], records: Iterable[tuple[str, str]],
               settings: SearchSettings, *, exhaustive: bool | None = None,
               progress: Callable[[str], None] | None = None,
               cancel: Callable[[], bool] | None = None) -> list[BindingSite]:
    """All binding sites of several primers on several records."""
    out: list[BindingSite] = []
    per_primer: dict[str, int] = {}
    for rid, seq in records:
        if cancel is not None and cancel():
            raise PrimerCancelled()
        if progress is not None:
            progress(f"Searching {rid} ({len(seq):,} bp)")
        for name, pseq in primers:
            if per_primer.get(name, 0) >= settings.max_sites_per_primer:
                continue
            hits = find_sites_in_sequence(name, pseq, rid, seq, settings, exhaustive=exhaustive)
            per_primer[name] = per_primer.get(name, 0) + len(hits)
            out.extend(hits)
    return out


class PrimerCancelled(Exception):
    pass


# --------------------------------------------------------------------------
# In-silico PCR
# --------------------------------------------------------------------------


def predict_products(sites: list[BindingSite], *, min_size: int = 1, max_size: int = 5000,
                     sequences: dict[str, str] | None = None,
                     allowed_pairs: set[tuple[str, str]] | None = None,
                     max_products: int = MAX_PRODUCTS,
                     keep_sequence_bp: int = 20_000) -> list[Product]:
    """Products between a "+" site and a downstream "-" site on the same record.

    ``allowed_pairs`` limits which (left primer, right primer) combinations
    count (e.g. only F+R of a designed pair, plus each primer on its own).
    Without it every combination counts, including one primer at both ends,
    which is what happens in a real tube.
    """
    by_record: dict[str, tuple[list[BindingSite], list[BindingSite]]] = {}
    for s in sites:
        fw, rv = by_record.setdefault(s.record, ([], []))
        (fw if s.strand == "+" else rv).append(s)
    products: list[Product] = []
    for rec, (fw, rv) in by_record.items():
        if not fw or not rv:
            continue
        fw.sort(key=lambda b: b.start)
        rv.sort(key=lambda b: b.end)
        rv_ends = np.asarray([b.end for b in rv], dtype=np.int64)
        for f in fw:
            lo = f.start + max(1, min_size) - 1
            hi = f.start + max_size - 1
            i0 = int(np.searchsorted(rv_ends, lo, side="left"))
            i1 = int(np.searchsorted(rv_ends, hi, side="right"))
            for r in rv[i0:i1]:
                if r.start < f.start:
                    continue
                if allowed_pairs is not None and (f.primer, r.primer) not in allowed_pairs:
                    continue
                size = r.end - f.start + 1
                if size < min_size or size > max_size:
                    continue
                seq = ""
                if sequences is not None and rec in sequences and size <= keep_sequence_bp:
                    seq = sequences[rec][f.start - 1:r.end].upper()
                products.append(Product(
                    record=rec, forward=f.primer, reverse=r.primer, start=f.start,
                    end=r.end, mismatches=f.mismatches + r.mismatches,
                    forward_start=f.start, reverse_end=r.end, sequence=seq,
                ))
                if len(products) >= max_products:
                    return products
    products.sort(key=lambda p: (p.record, p.start, p.end))
    return products


def products_for_pair(sites: list[BindingSite], forward: str, reverse: str, **kw) -> list[Product]:
    """Products a pair makes: F…R, plus F…F and R…R (single-primer products)."""
    allowed = {(forward, reverse), (forward, forward), (reverse, reverse), (reverse, forward)}
    own = [s for s in sites if s.primer in (forward, reverse)]
    return predict_products(own, allowed_pairs=allowed, **kw)


# --------------------------------------------------------------------------
# Alignment helpers (consensus, gap-free coordinates)
# --------------------------------------------------------------------------


def is_aligned(rows: list[str]) -> bool:
    if len(rows) < 2:
        return False
    lengths = {len(r) for r in rows}
    return len(lengths) == 1 and any("-" in r or "." in r for r in rows)


def alignment_masks(rows: list[str]) -> np.ndarray:
    """(n_seq, n_col) uint8 masks; gaps and N are 0."""
    width = max((len(r) for r in rows), default=0)
    mat = np.zeros((len(rows), width), dtype=np.uint8)
    for i, r in enumerate(rows):
        enc = encode_template(r.replace(".", "-"))
        mat[i, : len(enc)] = enc
    return mat


def base_counts(masks: np.ndarray) -> np.ndarray:
    """(n_col, 5) counts of A C G T gap (ambiguity spread evenly)."""
    n_col = masks.shape[1]
    counts = np.zeros((n_col, 5), dtype=np.float64)
    for bit_i, bit in enumerate((1, 2, 4, 8)):
        has = (masks & bit) != 0
        nbits = np.zeros_like(masks, dtype=np.float64)
        for b in (1, 2, 4, 8):
            nbits += (masks & b) != 0
        with np.errstate(divide="ignore", invalid="ignore"):
            share = np.where(has, 1.0 / np.maximum(nbits, 1), 0.0)
        counts[:, bit_i] = share.sum(axis=0)
    counts[:, 4] = (masks == 0).sum(axis=0)
    return counts


def consensus(rows: list[str], *, threshold: float = 0.0, ignore_gaps: bool = True) -> str:
    """Consensus of an alignment, gap columns dropped.

    ``threshold`` 0 = majority base. Above 0, the smallest set of bases that
    together reach the threshold is written as one IUPAC code (Geneious'
    "Threshold" consensus). With ``ignore_gaps`` gaps do not count toward the
    column total; a column that is mostly gaps is left out either way.
    """
    from repeatviewer.backends.primers.thermo import iupac_code

    masks = alignment_masks(rows)
    counts = base_counts(masks)
    out = []
    letters = "ACGT"
    for col in counts:
        bases = col[:4]
        gaps = col[4]
        total = bases.sum() + (0 if ignore_gaps else gaps)
        if bases.sum() == 0 or gaps > bases.sum():
            continue
        if threshold <= 0:
            out.append(letters[int(np.argmax(bases))])
            continue
        order = np.argsort(-bases)
        chosen, acc = [], 0.0
        for j in order:
            if bases[j] <= 0:
                break
            chosen.append(letters[j])
            acc += bases[j]
            if acc / max(total, 1e-9) >= threshold - 1e-9:
                break
        out.append(iupac_code(chosen))
    return "".join(out)


def check_primer_text(seq: str, name: str = "primer") -> str:
    b = binding_part(seq)
    if len(b) < 6:
        raise PrimerError(f"{name}: at least 6 binding bases are needed.")
    return seq


# --------------------------------------------------------------------------
# Genome scan: one pass over the sequence for all primers (k-mer seeds)
# --------------------------------------------------------------------------

_CODE2 = np.full(256, 255, dtype=np.uint8)
for _i, _b in enumerate("ACGT"):
    _CODE2[ord(_b)] = _i
    _CODE2[ord(_b.lower())] = _i


def _kmer_hashes(codes: np.ndarray, k: int) -> tuple[np.ndarray, np.ndarray]:
    """2-bit rolling hashes of every k-mer, and whether the k-mer is plain ACGT."""
    n = len(codes) - k + 1
    if n <= 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=bool)
    h = np.zeros(n, dtype=np.int64)
    bad = (codes == 255)
    c = np.where(bad, 0, codes).astype(np.int64)
    for j in range(k):
        h = (h << 2) | c[j:j + n]
    badc = np.concatenate([[0], np.cumsum(bad, dtype=np.int64)])
    valid = (badc[k:k + n] - badc[:n]) == 0
    return h, valid


def _hash_of(kmer: str) -> int:
    h = 0
    for ch in kmer:
        h = (h << 2) | "ACGT".index(ch)
    return h


@dataclass
class _SeedEntry:
    primer: str
    oriented: str
    strand: str
    offset: int          # seed start inside the oriented primer
    at_end: bool


class GenomeScanner:
    """Find all primers' binding sites in one pass per sequence chunk.

    Each primer is seeded by its 3' end (its last ``seed`` bases, or the first
    ones of its reverse complement), which must match exactly. Seeds are found
    for all primers at once through 2-bit k-mer hashes, then each candidate is
    scored like the exhaustive search.
    """

    def __init__(self, primers: list[tuple[str, str]], settings: SearchSettings,
                 seed: int = DEFAULT_SEED) -> None:
        self.settings = settings
        self.k = max(MIN_SEED, min(12, max(int(seed), settings.three_prime_exact)))
        self.entries: dict[int, list[_SeedEntry]] = {}
        self.max_len = 1
        for name, pseq in primers:
            bind = binding_part(pseq)
            if len(bind) < self.k:
                continue
            self.max_len = max(self.max_len, len(bind))
            orient = [("+", bind, True)]
            if settings.both_strands:
                orient.append(("-", reverse_complement(bind), False))
            for strand, oriented, at_end in orient:
                L = len(oriented)
                seed_seq, offset = ((oriented[L - self.k:], L - self.k) if at_end
                                    else (oriented[: self.k], 0))
                for variant in _expand_seed(seed_seq):
                    self.entries.setdefault(_hash_of(variant), []).append(
                        _SeedEntry(name, oriented, strand, offset, at_end))
        self._keys = np.asarray(sorted(self.entries), dtype=np.int64)
        self.counts: dict[str, int] = {name: 0 for name, _ in primers}

    def scan(self, record_id: str, seq: str) -> list[BindingSite]:
        if len(self._keys) == 0 or not seq:
            return []
        sites: list[BindingSite] = []
        n = len(seq)
        overlap = self.max_len
        pos = 0
        while pos < n:
            core_end = min(n, pos + CHUNK_BP)
            lo = max(0, pos - overlap)
            hi = min(n, core_end + overlap)
            chunk = seq[lo:hi]
            data = chunk.encode("ascii", "replace")
            raw = np.frombuffer(data, dtype=np.uint8)
            codes = _CODE2[raw]
            tmask = _MASK[raw]
            hashes, valid = _kmer_hashes(codes, self.k)
            idx = np.nonzero(valid & np.isin(hashes, self._keys, assume_unique=False))[0]
            if len(idx):
                hv = hashes[idx]
                order = np.argsort(hv, kind="stable")
                idx, hv = idx[order], hv[order]
                uniq, starts = np.unique(hv, return_index=True)
                bounds = list(starts) + [len(hv)]
                for u_i, key in enumerate(uniq.tolist()):
                    seed_pos = idx[bounds[u_i]:bounds[u_i + 1]]
                    for entry in self.entries.get(key, []):
                        sites.extend(self._score(entry, record_id, chunk, tmask, seed_pos,
                                                 lo, pos, core_end))
            pos = core_end
        sites.sort(key=lambda b: (b.start, b.strand))
        return sites

    def _score(self, entry: _SeedEntry, rid: str, chunk: str, tmask: np.ndarray,
               seed_pos: np.ndarray, lo: int, core_lo: int, core_hi: int) -> list[BindingSite]:
        name = entry.primer
        cap = self.settings.max_sites_per_primer
        if self.counts.get(name, 0) >= cap:
            return []
        L = len(entry.oriented)
        st = seed_pos.astype(np.int64) - entry.offset
        ok = (st >= 0) & (st + L <= len(tmask))
        st = st[ok]
        # only windows that start inside this chunk's core (no duplicates)
        abs_start = st + lo
        st = st[(abs_start >= core_lo) & (abs_start < core_hi)]
        if len(st) == 0:
            return []
        pm = encode_primer(entry.oriented)
        win = tmask[st[:, None] + np.arange(L)[None, :]]
        mm = (win & pm[None, :]) == 0
        counts = mm.sum(axis=1)
        keep = counts <= self.settings.max_mismatches
        ek = min(self.settings.three_prime_exact, L)
        if ek:
            tail = mm[:, L - ek:] if entry.at_end else mm[:, :ek]
            keep &= ~tail.any(axis=1)
        st, counts, mm = st[keep], counts[keep], mm[keep]
        room = cap - self.counts.get(name, 0)
        st, counts, mm = st[:room], counts[:room], mm[:room]
        self.counts[name] = self.counts.get(name, 0) + len(st)
        out = []
        for s, c, row in zip(st.tolist(), counts.tolist(), mm):
            window = chunk[s:s + L].upper()
            bad = np.nonzero(row)[0].tolist() if c else []
            if entry.strand == "+":
                positions = [b + 1 for b in bad]
                site = window
            else:
                positions = sorted(L - b for b in bad)
                site = reverse_complement(window)
            out.append(BindingSite(primer=name, record=rid, strand=entry.strand,
                                   start=lo + s + 1, end=lo + s + L, mismatches=int(c),
                                   mismatch_positions=positions, site=site))
        return out


def _expand_seed(seed: str, limit: int = 64) -> list[str]:
    from repeatviewer.backends.primers.thermo import expand_iupac

    return expand_iupac(seed.upper(), limit)


def search_records(primers: list[tuple[str, str]], records: Iterable[tuple[str, str]],
                   settings: SearchSettings, *, large: bool,
                   seed: int = DEFAULT_SEED,
                   progress: Callable[[str], None] | None = None,
                   cancel: Callable[[], bool] | None = None,
                   on_record: Callable[[str, str, list[BindingSite]], None] | None = None,
                   ) -> list[BindingSite]:
    """Binding sites on many records: exhaustive for small inputs, seeded scan for large.

    ``on_record(id, seq, sites)`` is called per record (e.g. to predict products
    while the sequence is in memory); results are also returned unless the
    callback is given.
    """
    scanner = GenomeScanner(primers, settings, seed) if large else None
    out: list[BindingSite] = []
    for rid, seq in records:
        if cancel is not None and cancel():
            raise PrimerCancelled()
        if progress is not None:
            progress(f"Searching {rid} ({len(seq):,} bp)")
        if scanner is not None:
            sites = scanner.scan(rid, seq)
        else:
            sites = []
            for name, pseq in primers:
                sites.extend(find_sites_in_sequence(name, pseq, rid, seq, settings,
                                                    exhaustive=True))
            sites.sort(key=lambda b: (b.start, b.strand))
        if on_record is not None:
            on_record(rid, seq, sites)
        else:
            out.extend(sites)
    return out
