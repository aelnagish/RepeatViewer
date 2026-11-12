"""Find a DNA or protein motif in sequences: the engine behind Edit > Find.

No Qt here, so everything is unit tested. Positions are 1-based and inclusive,
in the coordinates of the string that was searched: alignment rows keep their
gap columns, so a hit can be selected directly in the viewer.

Options follow Geneious' Find in Document:
  * both strands: also report hits of the reverse complement (strand "-")
  * IUPAC codes in the query: R matches A or G, N matches anything, ...
  * IUPAC codes in the sequence: an N (or R, Y, ...) in the sequence matches
  * allow gaps: ignore '-' and '.' in the sequence, so a motif is found even
    when an alignment gap splits it
  * mismatches: up to this many positions may differ
  * protein: the query is a peptide, searched in all six reading frames

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import re
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field

import numpy as np

IUPAC_DNA: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
_DNA_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVN", "TGCAAYRSWMKVHDBN"
)
# Amino acids, the stop '*', and the ambiguity codes B (D/N), Z (E/Q), J (I/L), X (any).
_AMINO = "ACDEFGHIKLMNPQRSTVWY"
IUPAC_PROTEIN: dict[str, str] = {aa: aa for aa in _AMINO}
IUPAC_PROTEIN.update({"*": "*", "B": "DN", "Z": "EQ", "J": "IL", "X": _AMINO, "U": "U", "O": "O"})
GAP_CHARS = "-."
_CHUNK = 4_000_000  # windows per numpy chunk (bounds memory on long scaffolds)


class FindCancelled(Exception):
    """The user pressed Stop."""


@dataclass(frozen=True)
class FindOptions:
    protein: bool = False
    both_strands: bool = True
    iupac_in_query: bool = True
    iupac_in_sequence: bool = False
    allow_gaps: bool = False
    max_mismatches: int = 0
    genetic_code: int = 1
    max_hits: int = 100_000


@dataclass(frozen=True)
class MotifHit:
    start: int  # 1-based, inclusive
    end: int  # 1-based, inclusive
    strand: str  # "+" or "-"
    mismatches: int = 0
    frame: int = 0  # protein search: +1..+3 / -1..-3; DNA: 0
    seq_id: str = ""


@dataclass
class FindResult:
    hits: list[MotifHit] = field(default_factory=list)
    truncated: bool = False  # stopped at FindOptions.max_hits


# ------------------------------------------------------------------ queries


def clean_query(query: str) -> str:
    """Upper case, without spaces, line breaks or position numbers (pasted text)."""
    return re.sub(r"[\s\d]+", "", query or "").upper()


def query_problem(query: str, options: FindOptions = FindOptions()) -> str | None:
    """A message for the user if ``query`` cannot be searched, else None."""
    q = clean_query(query)
    if not q:
        return "Type a sequence motif to search for."
    alphabet = IUPAC_PROTEIN if options.protein else IUPAC_DNA
    bad = sorted({c for c in q if c not in alphabet})
    if bad:
        kind = "amino acid" if options.protein else "DNA or IUPAC"
        return f"Not {kind} letters: {', '.join(bad)}"
    if options.max_mismatches < 0:
        return "Mismatches cannot be negative."
    if options.max_mismatches >= len(q):
        return "Allow fewer mismatches than the motif is long."
    return None


def reverse_complement(seq: str) -> str:
    return seq.upper().translate(_DNA_COMPLEMENT)[::-1]


def _allowed_sets(query: str, options: FindOptions) -> list[frozenset[str]]:
    """For each query position, the sequence letters that count as a match."""
    table = IUPAC_PROTEIN if options.protein else IUPAC_DNA
    sets = []
    for q in query:
        q_bases = set(table.get(q, q)) if options.iupac_in_query else {table.get(q, q) if len(table.get(q, q)) == 1 else q}
        allowed = set(q_bases)
        if not options.iupac_in_query:
            allowed.add(q)
        if options.iupac_in_sequence:
            for code, code_bases in table.items():
                if set(code_bases) & q_bases:
                    allowed.add(code)
        if not options.protein and "T" in allowed:
            allowed.add("U")
        sets.append(frozenset(allowed))
    return sets


# ------------------------------------------------------------------ matching


def _match_positions(
    text: str,
    allowed: list[frozenset[str]],
    max_mismatches: int,
    limit: int,
    cancelled: Callable[[], bool] | None,
) -> tuple[list[tuple[int, int]], bool]:
    """0-based starts (and mismatch counts) where ``allowed`` matches ``text``.

    Overlapping hits are all reported. Returns (hits, truncated).
    """
    length = len(allowed)
    n = len(text)
    hits: list[tuple[int, int]] = []
    if length == 0 or n < length:
        return hits, False
    if max_mismatches <= 0:
        pattern = "".join(
            "[" + "".join(re.escape(c) for c in sorted(chars)) + "]" for chars in allowed
        )
        regex = re.compile(f"(?=({pattern}))")
        # Scan in chunks (each overlapping the next by length - 1) so Stop is
        # noticed on long scaffolds even when there are no hits.
        step = _CHUNK * 2
        for chunk_start in range(0, n - length + 1, step):
            if cancelled is not None and cancelled():
                raise FindCancelled()
            chunk_stop = min(n - length + 1, chunk_start + step)
            for match in regex.finditer(text, chunk_start, chunk_stop + length - 1):
                if match.start() >= chunk_stop:
                    break
                if len(hits) >= limit:
                    return hits, True
                hits.append((match.start(), 0))
        return hits, False

    data = np.frombuffer(text.encode("latin-1", "replace"), dtype=np.uint8)
    lookups = []
    for chars in allowed:
        ok = np.zeros(256, dtype=bool)
        for c in chars:
            ok[ord(c)] = True
            ok[ord(c.lower())] = True
        lookups.append(ok)
    windows = n - length + 1
    for chunk_start in range(0, windows, _CHUNK):
        if cancelled is not None and cancelled():
            raise FindCancelled()
        chunk_end = min(windows, chunk_start + _CHUNK)
        size = chunk_end - chunk_start
        misses = np.zeros(size, dtype=np.uint16)
        for i, ok in enumerate(lookups):
            misses += ~ok[data[chunk_start + i : chunk_start + i + size]]
        found = np.flatnonzero(misses <= max_mismatches)
        for pos in found:
            if len(hits) >= limit:
                return hits, True
            hits.append((chunk_start + int(pos), int(misses[pos])))
    return hits, False


def _ungap(sequence: str, allow_gaps: bool) -> tuple[str, list[int] | None]:
    """Sequence without gap characters, and the original index of each kept letter."""
    if not allow_gaps or not any(g in sequence for g in GAP_CHARS):
        return sequence, None
    keep = [i for i, c in enumerate(sequence) if c not in GAP_CHARS]
    return "".join(sequence[i] for i in keep), keep


# ------------------------------------------------------------------ translation


_CODON_LETTERS = "ACGTRYSWKMBDHVN"  # index 0..14; anything else (gaps, ...) is 15
# Codons whose readings are exactly one of these pairs get the ambiguity letter.
_AMBIGUOUS_AMINO = {frozenset("DN"): "B", frozenset("EQ"): "Z", frozenset("IL"): "J"}


def _codon_table(genetic_code: int) -> np.ndarray:
    """Amino acid for every codon written with IUPAC letters (16 x 16 x 16 entries).

    Index c1*256 + c2*16 + c3 with each letter's position in _CODON_LETTERS.
    An ambiguous codon still gets an amino acid when all its readings agree
    (GGN is always G, TAR is always a stop), B/Z/J when they are D-or-N,
    E-or-Q or I-or-L, and X otherwise, exactly as Biopython translates.
    """
    from itertools import product

    from Bio.Data import CodonTable

    table = CodonTable.unambiguous_dna_by_id[genetic_code]

    def amino(codon: str) -> str:
        return "*" if codon in table.stop_codons else table.forward_table.get(codon, "X")

    out = np.full(16 * 16 * 16, ord("X"), dtype=np.uint8)
    for i1, l1 in enumerate(_CODON_LETTERS):
        for i2, l2 in enumerate(_CODON_LETTERS):
            for i3, l3 in enumerate(_CODON_LETTERS):
                readings = {
                    amino(a + b + c)
                    for a, b, c in product(IUPAC_DNA[l1], IUPAC_DNA[l2], IUPAC_DNA[l3])
                }
                if len(readings) == 1:
                    out[i1 * 256 + i2 * 16 + i3] = ord(readings.pop())
                elif frozenset(readings) in _AMBIGUOUS_AMINO:
                    out[i1 * 256 + i2 * 16 + i3] = ord(_AMBIGUOUS_AMINO[frozenset(readings)])
    return out


def translate_frame(dna: str, offset: int, codon_table: np.ndarray) -> str:
    """Translate ``dna`` from ``offset`` (0, 1 or 2). Unreadable codons give X."""
    code = np.full(256, 15, dtype=np.uint16)
    for i, letter in enumerate(_CODON_LETTERS):
        code[ord(letter)] = i
        code[ord(letter.lower())] = i
    code[ord("U")] = code[ord("u")] = 3
    data = np.frombuffer(dna.encode("latin-1", "replace"), dtype=np.uint8)[offset:]
    usable = (len(data) // 3) * 3
    if usable == 0:
        return ""
    codons = code[data[:usable]].reshape(-1, 3)
    index = codons[:, 0] * 256
    index += codons[:, 1] * 16
    index += codons[:, 2]
    return codon_table[index].tobytes().decode("ascii")


# ------------------------------------------------------------------ search


def find_in_sequence(
    sequence: str,
    query: str,
    options: FindOptions = FindOptions(),
    *,
    seq_id: str = "",
    cancelled: Callable[[], bool] | None = None,
    limit: int | None = None,
) -> FindResult:
    """All hits of ``query`` in one sequence, sorted by position."""
    q = clean_query(query)
    problem = query_problem(q, options)
    if problem:
        raise ValueError(problem)
    limit = options.max_hits if limit is None else limit
    text, keep = _ungap(sequence.upper(), options.allow_gaps or options.protein)
    n = len(text)

    def original(pos0: int) -> int:
        return (keep[pos0] if keep is not None else pos0) + 1

    hits: list[MotifHit] = []
    truncated = False
    if not options.protein:
        allowed = _allowed_sets(q, options)
        strands = [("+", allowed)]
        if options.both_strands:
            rc = reverse_complement(q)
            if rc != q:  # palindromes (e.g. GAATTC) would only repeat the + hits
                strands.append(("-", _allowed_sets(rc, options)))
        for strand, sets in strands:
            found, cut = _match_positions(text, sets, options.max_mismatches, limit - len(hits), cancelled)
            truncated = truncated or cut
            for start0, misses in found:
                end0 = start0 + len(q) - 1
                hits.append(MotifHit(original(start0), original(end0), strand, misses, 0, seq_id))
    else:
        allowed = _allowed_sets(q, options)
        codons = _codon_table(options.genetic_code)
        frames = [("+", text)]
        if options.both_strands:
            frames.append(("-", reverse_complement(text)))
        for strand, dna in frames:
            for offset in range(3):
                if cancelled is not None and cancelled():
                    raise FindCancelled()
                peptide = translate_frame(dna, offset, codons)
                found, cut = _match_positions(
                    peptide, allowed, options.max_mismatches, limit - len(hits), cancelled
                )
                truncated = truncated or cut
                for start0, misses in found:
                    nt_start = offset + 3 * start0
                    nt_end = offset + 3 * (start0 + len(q)) - 1
                    if strand == "-":
                        nt_start, nt_end = n - 1 - nt_end, n - 1 - nt_start
                    frame = (offset + 1) if strand == "+" else -(offset + 1)
                    hits.append(
                        MotifHit(original(nt_start), original(nt_end), strand, misses, frame, seq_id)
                    )
    hits.sort(key=lambda h: (h.start, h.end, h.strand != "+"))
    return FindResult(hits[:limit], truncated or len(hits) > limit)


def find_in_records(
    records: Iterable[tuple[str, str]],
    query: str,
    options: FindOptions = FindOptions(),
    *,
    cancelled: Callable[[], bool] | None = None,
) -> FindResult:
    """Search every (id, sequence) record; hits keep their record id."""
    result = FindResult()
    for seq_id, sequence in records:
        left = options.max_hits - len(result.hits)
        if left <= 0:
            result.truncated = True
            break
        part = find_in_sequence(sequence, query, options, seq_id=seq_id, cancelled=cancelled, limit=left)
        result.hits.extend(part.hits)
        result.truncated = result.truncated or part.truncated
    return result


# ------------------------------------------------------------------ positions


_POSITION = re.compile(
    r"^\s*(?:(?P<name>[^\s:]+)\s*:\s*)?(?P<start>[\d,_ ]*\d)\s*(?:(?:-|\.\.|–|to)\s*(?P<end>[\d,_ ]*\d))?\s*$",
    re.IGNORECASE,
)


def parse_position(text: str) -> tuple[str | None, int, int] | None:
    """Read a typed position for Go To Position: "114590001", "114,590,001",
    "1000-2000", "1,000..2,000", "scaffold_5:1000-2000". Returns (name, start,
    end) with start <= end, or None if the text is not a position."""
    match = _POSITION.match(text or "")
    if not match:
        return None

    def number(raw: str) -> int:
        return int(re.sub(r"[,_ ]", "", raw))

    start = number(match.group("start"))
    end = number(match.group("end")) if match.group("end") else start
    if start < 1 or end < 1:
        return None
    if end < start:
        start, end = end, start
    return match.group("name"), start, end
