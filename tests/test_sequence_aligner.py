"""Clustal Omega threading and alignment helpers."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from teviewer.backends.sequence_aligner import (
    _fasta_sequence_count,
    clustalo_effective_threads,
    run_alignment,
)


def test_fasta_sequence_count(tmp_path: Path) -> None:
    fasta = tmp_path / "t.fa"
    fasta.write_text(">a\nAC\n>b\nGT\n")
    assert _fasta_sequence_count(fasta) == 2


def test_run_alignment_requires_two_sequences(tmp_path: Path) -> None:
    fasta = tmp_path / "one.fa"
    fasta.write_text(">only\nACGT\n")
    with pytest.raises(ValueError, match="fewer than 2 sequences"):
        run_alignment(fasta)


def test_clustalo_effective_threads_without_openmp() -> None:
    with patch(
        "teviewer.backends.sequence_aligner.clustalo_supports_openmp",
        return_value=False,
    ):
        assert clustalo_effective_threads(5) == 1


def test_clustalo_effective_threads_with_openmp() -> None:
    with patch(
        "teviewer.backends.sequence_aligner.clustalo_supports_openmp",
        return_value=True,
    ):
        assert clustalo_effective_threads(5) == 5


def test_run_alignment_clustalo_on_real_input(tmp_path: Path) -> None:
    if not __import__("shutil").which("clustalo"):
        pytest.skip("clustalo not installed")
    fasta = tmp_path / "pair.fa"
    fasta.write_text(">a\nACGT\n>b\nACGT\n")
    out = run_alignment(fasta, backend="clustalo", num_threads=5)
    assert out.is_file()
    assert out.name == "pair_msa.fa"
