"""Satellite monomer phase rotation (pure)."""

from __future__ import annotations

import random

from repeatviewer.backends import monomer_phase as mp


def _monomer(length=156, seed=1):
    rng = random.Random(seed)
    return "".join(rng.choice("ACGT") for _ in range(length))


def _mutate(seq, n, seed):
    rng = random.Random(seed)
    s = list(seq)
    for _ in range(n):
        i = rng.randrange(len(s))
        s[i] = rng.choice("ACGT".replace(s[i], ""))
    return "".join(s)


def test_rotations_and_strand_are_undone():
    base = _monomer()
    records = [("m0", base)]
    for i, off in enumerate((17, 60, 101, 150), 1):
        seq = mp.rotate(_mutate(base, 8, i), off)
        if i % 2 == 0:
            seq = mp.reverse_complement(seq)
        records.append((f"m{i}", seq))
    result = mp.rotate_to_common_start(records)
    assert result.reference == "m0"
    assert not result.unmatched
    for (_name, seq), (_n, orig) in zip(result.sequences, [("m0", base)] + [
            (f"m{i}", _mutate(base, 8, i)) for i in range(1, 5)]):
        assert seq == orig
    assert result.n_reversed == 2
    assert result.n_rotated == 4


def test_reference_is_the_median_length_monomer():
    seqs = ["A" * 100, "A" * 150, "A" * 160, "A" * 170, "A" * 400]
    assert mp.choose_reference(seqs) == 2


def test_unrelated_sequence_is_left_unchanged_and_reported():
    base = _monomer(seed=2)
    other = _monomer(seed=99)
    result = mp.rotate_to_common_start([("a", base), ("b", mp.rotate(base, 30)), ("junk", other)])
    names = dict(result.sequences)
    assert names["junk"] == other
    assert result.unmatched == ["junk"]


def test_small_indels_still_find_the_register():
    base = _monomer(seed=5)
    shifted = mp.rotate(base[:40] + base[43:], 70)   # 3 bp deletion, then rotation
    result = mp.rotate_to_common_start([("ref", base), ("del", shifted), ("ref2", base)])
    rotated = dict(result.sequences)["del"]
    assert rotated.startswith(base[:30])


def test_short_monomers_use_a_smaller_k():
    assert mp.pick_k(["ACGT" * 10]) == mp.SHORT_K
    assert mp.pick_k(["ACGT" * 40]) == mp.DEFAULT_K


def test_table_lists_every_monomer():
    base = _monomer(seed=7)
    result = mp.rotate_to_common_start([("a", base), ("b", mp.rotate(base, 12))])
    lines = result.table_tsv().strip().splitlines()
    assert lines[0].startswith("name\tstrand\toffset")
    assert lines[2] == "b\t+\t144\t" + lines[2].split("\t")[3] + "\tyes"


def test_unrelated_random_sequences_are_not_rotated():
    rng = random.Random(11)
    base = _monomer(seed=12)
    junk = [(f"j{i}", "".join(rng.choice("ACGT") for _ in range(156))) for i in range(300)]
    result = mp.rotate_to_common_start([("b", base), ("b2", mp.rotate(base, 40))] + junk)
    wrongly = [n for n in (r.name for r in result.rows if r.matched) if n.startswith("j")]
    assert len(wrongly) <= 2          # expected about 0.3 at P_RANDOM = 0.001 per k
    assert "b2" not in result.unmatched


def test_threshold_grows_as_k_shrinks():
    assert mp.vote_threshold(156, 156, 8) < mp.vote_threshold(156, 156, 6) < mp.vote_threshold(156, 156, 5)


def test_poisson_tail():
    assert mp._poisson_sf(0, 2.0) == 1.0
    assert abs(mp._poisson_sf(1, 2.0) - (1 - 2.718281828 ** -2)) < 1e-9
