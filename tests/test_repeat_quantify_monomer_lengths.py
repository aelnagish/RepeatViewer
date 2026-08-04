"""Monomer-length lookup must resolve names from both quantifier engines.

blastn reports the full library header as ``sseqid``
(``run2_TRC_7#Simple_repeat/TRC_7``); RepeatMasker's ``.out`` splits that on
``#`` and reports only ``run2_TRC_7``. A full-id-only lookup silently yields
monomer_len 0, which zeroes copy_number too.
"""

from pathlib import Path

import pytest

from teviewer.backends.repeat_quantify import (
    _aggregate_rm_rows,
    _fasta_lengths,
)


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> Path:
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
    return path


def test_full_id_lookup_still_works(tmp_path):
    """The blastn path keys on the whole header and must be unaffected."""
    lib = _write_fasta(
        tmp_path / "lib.fa",
        [("run2_TRC_7#Simple_repeat/TRC_7", "A" * 302)],
    )
    lengths = _fasta_lengths(lib)
    assert lengths["run2_TRC_7#Simple_repeat/TRC_7"] == 302


def test_stem_is_registered_for_repeatmasker_names(tmp_path):
    """RepeatMasker's bare name resolves via the pre-'#' stem alias."""
    lib = _write_fasta(
        tmp_path / "lib.fa",
        [
            ("run2_TRC_7#Simple_repeat/TRC_7", "A" * 302),
            ("run1_TRC_4#Simple_repeat/TRC_4", "C" * 82),
        ],
    )
    lengths = _fasta_lengths(lib)
    assert lengths["run2_TRC_7"] == 302
    assert lengths["run1_TRC_4"] == 82


def test_headers_without_class_suffix_add_no_alias(tmp_path):
    lib = _write_fasta(tmp_path / "lib.fa", [("CroSat2_mono", "G" * 150)])
    assert _fasta_lengths(lib) == {"CroSat2_mono": 150}


def test_ambiguous_stem_is_skipped(tmp_path):
    """Two entries sharing a stem must not resolve to an arbitrary length."""
    lib = _write_fasta(
        tmp_path / "lib.fa",
        [("sat#Satellite/A", "A" * 100), ("sat#Satellite/B", "T" * 250)],
    )
    lengths = _fasta_lengths(lib)
    assert "sat" not in lengths
    assert lengths["sat#Satellite/A"] == 100
    assert lengths["sat#Satellite/B"] == 250


def test_exact_id_wins_over_stem(tmp_path):
    """A real record named 'sat' must not be shadowed by the 'sat#X' alias."""
    lib = _write_fasta(
        tmp_path / "lib.fa",
        [("sat", "A" * 11), ("sat#Satellite/A", "T" * 99)],
    )
    assert _fasta_lengths(lib)["sat"] == 11


def test_duplicate_ids_keep_last_length(tmp_path):
    """library_duplicate_warnings() documents last-wins; preserve that."""
    lib = _write_fasta(
        tmp_path / "lib.fa",
        [("dup#Simple_repeat/d", "A" * 10), ("dup#Simple_repeat/d", "T" * 40)],
    )
    lengths = _fasta_lengths(lib)
    assert lengths["dup#Simple_repeat/d"] == 40
    assert lengths["dup"] == 40


def test_rm_aggregate_populates_monomer_len_and_copy_number(tmp_path):
    """End-to-end: the RepeatMasker aggregate now fills both columns."""
    lib = _write_fasta(
        tmp_path / "lib.fa",
        [("run2_TRC_7#Simple_repeat/TRC_7", "A" * 100)],
    )
    rows = [
        {"family": "run2_TRC_7", "scaffold": "chr1",
         "start": 1, "end": 500, "divergence": 5.0},
        {"family": "run2_TRC_7", "scaffold": "chr1",
         "start": 1000, "end": 1499, "divergence": 5.0},
    ]
    families = _aggregate_rm_rows(
        rows, target_length=10_000, monomer_lengths=_fasta_lengths(lib)
    )
    assert len(families) == 1
    fam = families[0]
    assert fam.covered_bp == 1000
    assert fam.monomer_len == 100
    assert fam.copy_number == 10  # covered_bp / monomer_len


def test_rm_aggregate_without_a_match_stays_zero(tmp_path):
    """An unknown family still degrades gracefully rather than dividing by 0."""
    lib = _write_fasta(tmp_path / "lib.fa", [("other#Satellite/o", "A" * 100)])
    rows = [
        {"family": "not_in_library", "scaffold": "chr1",
         "start": 1, "end": 500, "divergence": 5.0},
    ]
    fam = _aggregate_rm_rows(
        rows, target_length=10_000, monomer_lengths=_fasta_lengths(lib)
    )[0]
    assert fam.monomer_len == 0
    assert fam.copy_number == 0


if __name__ == "__main__":
    raise SystemExit(pytest.main([__file__, "-v"]))
