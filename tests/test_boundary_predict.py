"""TEtrimmer boundary prediction port (cleaning + MSA identity, no genome)."""

from __future__ import annotations

import random

import numpy as np

from repeatviewer.backends import boundary_predict as bp


def _random_msa(rows=20, core=300, flank=200, seed=1):
    rng = random.Random(seed)
    te = "".join(rng.choice("acgt") for _ in range(core))
    seqs = []
    for _ in range(rows):
        left = "".join(rng.choice("acgt") for _ in range(flank))
        right = "".join(rng.choice("acgt") for _ in range(flank))
        body = list(te)
        for _ in range(core // 50):  # ~2% divergence
            body[rng.randrange(core)] = rng.choice("acgt")
        seqs.append(left + "".join(body) + right)
    return seqs


def test_clean_column_turns_rare_nucleotides_into_gaps():
    # 101 rows: one "t" among 100 "a" is below 1% → gap; uppercase is untouched.
    col = ["a"] * 100 + ["t"]
    m = np.array([[c] for c in col])
    out = bp.clean_column(m, threshold=0.01)
    assert out[-1, 0] == "-"
    assert (out[:-1, 0] == "a").all()
    upper = np.array([[c.upper()] for c in col])
    assert bp.clean_column(upper)[-1, 0] == "T"


def test_remove_gaps_needs_fifty_columns():
    m = bp._to_matrix(["acgt" * 10] * 6)
    assert bp.remove_gaps_with_similarity_check(m) == (None, None)


def test_remove_gaps_mapping_has_exclusive_end():
    seqs = ["a" * 60 + "-" * 5] * 6
    gs, cmap = bp.remove_gaps_with_similarity_check(bp._to_matrix(seqs))
    assert gs.shape == (6, 60)
    assert cmap[0] == 0 and cmap[59] == 59 and cmap[60] == 60


def test_boundary_rows_skip_only_tetrimmer_prefixed_rows():
    ids = [
        "TEtrimmer_output_fam#LTR/Copia",
        "TEtrimmer_input_fam",
        "chr1:1-10(+)",
        "fam#LTR/Copia",  # no prefix: a real row
        "Input_fam",
    ]
    assert bp.boundary_rows(ids, [""] * 5) == [2, 3, 4]


def test_msa_identity_finds_conserved_core_without_genome():
    seqs = _random_msa()
    ids = [f"chr1:{i}-{i + 700}(+)" for i in range(len(seqs))]
    result = bp.predict_boundaries(ids, seqs, genome_path=None)
    start, end = result.msa_identity
    # Random flanks are not conserved; the 300-column core sits at 200–500.
    assert abs(start - 200) <= 5
    assert abs(end - 500) <= 5
    assert result.final is None
    assert any("No project genome" in n for n in result.notes)
