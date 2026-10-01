"""Edit > Find engine: DNA/protein motifs, strands, IUPAC, gaps, mismatches, positions."""

from __future__ import annotations

import random
import time

import pytest

from repeatviewer.backends.motif_search import (
    IUPAC_DNA,
    FindCancelled,
    FindOptions,
    clean_query,
    find_in_records,
    find_in_sequence,
    parse_position,
    query_problem,
    reverse_complement,
    translate_frame,
    _codon_table,
)


def spans(result, strand=None):
    return [(h.start, h.end) for h in result.hits if strand is None or h.strand == strand]


# ------------------------------------------------------------------ DNA basics


def test_forward_hits_are_one_based_and_overlapping():
    r = find_in_sequence("GGAAAAGG", "AA", FindOptions(both_strands=False))
    assert spans(r) == [(3, 4), (4, 5), (5, 6)]


def test_reverse_complement_hit_is_on_minus_strand():
    # ACGTT reverse-complemented is AACGT
    r = find_in_sequence("GGGAACGTGGG", "ACGTT")
    assert [(h.start, h.end, h.strand) for h in r.hits] == [(4, 8, "-")]
    assert spans(find_in_sequence("GGGAACGTGGG", "ACGTT", FindOptions(both_strands=False))) == []


def test_palindrome_is_reported_once():
    r = find_in_sequence("TTGAATTCTT", "GAATTC")
    assert [(h.start, h.end, h.strand) for h in r.hits] == [(3, 8, "+")]


def test_soft_masked_lowercase_and_rna_letters_match():
    assert spans(find_in_sequence("ccacgtcc", "ACGT", FindOptions(both_strands=False))) == [(3, 6)]
    assert spans(find_in_sequence("CCACGUCC", "ACGT", FindOptions(both_strands=False))) == [(3, 6)]
    assert spans(find_in_sequence("CCACGTCC", "ACGU", FindOptions(both_strands=False))) == [(3, 6)]


def test_iupac_in_query():
    on = FindOptions(both_strands=False)
    assert spans(find_in_sequence("AAGTT", "RGT", on)) == [(2, 4)]          # R = A/G
    off = FindOptions(both_strands=False, iupac_in_query=False)
    assert spans(find_in_sequence("AAGTT", "RGT", off)) == []
    assert spans(find_in_sequence("ARGTT", "RGT", off)) == [(2, 4)]         # literal R


def test_iupac_in_sequence():
    seq = "TTANGTT"
    off = FindOptions(both_strands=False)
    assert spans(find_in_sequence(seq, "ACG", off)) == []
    on = FindOptions(both_strands=False, iupac_in_sequence=True)
    assert spans(find_in_sequence(seq, "ACG", on)) == [(3, 5)]              # N matches C


def test_gaps_in_alignment_rows():
    row = "TTAC--G-TTT"
    off = FindOptions(both_strands=False)
    assert spans(find_in_sequence(row, "ACGT", off)) == []
    on = FindOptions(both_strands=False, allow_gaps=True)
    assert spans(find_in_sequence(row, "ACGT", on)) == [(3, 9)]             # gapped columns


def test_mismatches_are_counted():
    opts = FindOptions(both_strands=False, max_mismatches=1)
    r = find_in_sequence("AAACCTTGGAAA", "CCATGG", opts)
    assert [(h.start, h.end, h.mismatches) for h in r.hits] == [(4, 9, 1)]
    assert find_in_sequence("AAACCTTGGAAA", "CCATGG", FindOptions(both_strands=False)).hits == []


# ------------------------------------------------------------------ vs brute force


def _naive(seq, query, opts):
    """Deliberately simple reference: IUPAC in query only, N etc. in the sequence mismatch."""
    upper = seq.upper().replace("U", "T")
    keep = [i for i, c in enumerate(upper) if not (opts.allow_gaps and c in "-.")]
    text = "".join(upper[i] for i in keep)
    found = set()
    strands = [("+", query)]
    if opts.both_strands and reverse_complement(query) != query:
        strands.append(("-", reverse_complement(query)))
    for strand, q in strands:
        for s in range(len(text) - len(q) + 1):
            misses = sum(text[s + i] not in IUPAC_DNA[q[i]] for i in range(len(q)))
            if misses <= opts.max_mismatches:
                found.add((keep[s] + 1, keep[s + len(q) - 1] + 1, strand, misses))
    return found


@pytest.mark.parametrize("seed", range(40))
def test_matches_brute_force_on_random_sequences(seed):
    rng = random.Random(seed)
    alphabet = "ACGT" * 6 + "N" + "acgt" + ("-" if seed % 2 else "")
    seq = "".join(rng.choice(alphabet) for _ in range(rng.randint(50, 400)))
    query = "".join(rng.choice("ACGT" * 4 + "RYN") for _ in range(rng.randint(3, 8)))
    opts = FindOptions(
        both_strands=bool(seed % 3),
        allow_gaps=bool(seed % 2),
        max_mismatches=rng.randint(0, min(2, len(query) - 1)),
    )
    got = {(h.start, h.end, h.strand, h.mismatches) for h in find_in_sequence(seq, query, opts).hits}
    assert got == _naive(seq, query, opts)


# ------------------------------------------------------------------ protein


def test_translation_frames_stops_and_unknown_codons():
    table = _codon_table(1)
    assert translate_frame("ATGAAATGGTAA", 0, table) == "MKW*"
    assert translate_frame("CATGAAA", 1, table) == "MK"
    assert translate_frame("ATGNAATGG", 0, table) == "MXW"   # NAA: K, Q, E or stop
    assert translate_frame("GGNTAR", 0, table) == "G*"       # all readings agree
    assert translate_frame("RAYSARMTH", 0, table) == "BZJ"   # D/N, E/Q, I/L like Biopython


def test_protein_hits_on_both_strands():
    opts = FindOptions(protein=True)
    plus = find_in_sequence("CC" + "ATGAAATGG" + "TT", "MKW", opts)
    assert [(h.start, h.end, h.strand, h.frame) for h in plus.hits] == [(3, 11, "+", 3)]
    minus_seq = "G" + reverse_complement("ATGAAATGG") + "AA"
    minus = find_in_sequence(minus_seq, "MKW", opts)
    assert [(h.start, h.end, h.strand) for h in minus.hits] == [(2, 10, "-")]


def test_protein_wildcard_and_mismatch():
    seq = "ATGAAATGG"
    assert spans(find_in_sequence(seq, "MXW", FindOptions(protein=True, both_strands=False))) == [(1, 9)]
    one = FindOptions(protein=True, both_strands=False, max_mismatches=1)
    assert [(h.start, h.mismatches) for h in find_in_sequence(seq, "MRW", one).hits] == [(1, 1)]


# ------------------------------------------------------------------ records, limits, stop


def test_records_keep_their_ids_and_the_hit_limit_holds():
    records = [("s1", "ACGTACGT"), ("s2", "TTACGTT")]
    r = find_in_records(records, "ACGT", FindOptions(both_strands=False))
    assert [(h.seq_id, h.start) for h in r.hits] == [("s1", 1), ("s1", 5), ("s2", 3)]
    capped = find_in_records(records, "ACGT", FindOptions(both_strands=False, max_hits=2))
    assert len(capped.hits) == 2 and capped.truncated


def test_stop_interrupts_a_long_search():
    seq = "ACGT" * 5_000_000
    with pytest.raises(FindCancelled):
        find_in_sequence(seq, "GGGGGGGG", FindOptions(), cancelled=lambda: True)
    with pytest.raises(FindCancelled):
        find_in_sequence(seq, "GGGGGGGG", FindOptions(max_mismatches=1), cancelled=lambda: True)


def test_twenty_megabase_scaffold_is_fast():
    rng = random.Random(1)
    seq = "".join(rng.choices("ACGT", k=20_000_000))
    t0 = time.perf_counter()
    exact = find_in_sequence(seq, "ACGTACGTACGT")
    t1 = time.perf_counter()
    fuzzy = find_in_sequence(seq, "ACGTACGTACGTAC", FindOptions(max_mismatches=2))
    t2 = time.perf_counter()
    assert exact.hits and fuzzy.hits
    assert t1 - t0 < 30 and t2 - t1 < 60, (t1 - t0, t2 - t1)


# ------------------------------------------------------------------ user input


def test_query_cleanup_and_problems():
    assert clean_query(" acg t\n12 NN ") == "ACGTNN"
    assert query_problem("") == "Type a sequence motif to search for."
    assert query_problem("ACGQ") == "Not DNA or IUPAC letters: Q"
    assert query_problem("ACG", FindOptions(max_mismatches=3)) == "Allow fewer mismatches than the motif is long."
    assert query_problem("MKW*", FindOptions(protein=True)) is None
    with pytest.raises(ValueError):
        find_in_sequence("ACGT", "ACGQ")


@pytest.mark.parametrize(
    "text, expected",
    [
        ("114590001", (None, 114590001, 114590001)),
        ("114,590,001", (None, 114590001, 114590001)),
        ("1000-2000", (None, 1000, 2000)),
        ("1,000 .. 2,000", (None, 1000, 2000)),
        ("2000-1000", (None, 1000, 2000)),
        ("scaffold_5:1000-2000", ("scaffold_5", 1000, 2000)),
        ("ACGT", None),
        ("0", None),
        ("", None),
    ],
)
def test_parse_position(text, expected):
    assert parse_position(text) == expected
