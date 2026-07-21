"""Clustal Omega threading and alignment helpers."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest

from teviewer.backends.sequence_aligner import (
    _fasta_sequence_count,
    bundled_mafft_command,
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


@pytest.mark.parametrize("os_type", ["Darwin", "Linux", "Windows"])
def test_bundled_mafft_exists_for_supported_platform(os_type: str) -> None:
    resolved = bundled_mafft_command(os_type, machine="x86_64")
    assert resolved is not None
    assert Path(resolved[0][0]).is_file()


def test_bundled_linux_mafft_rejects_unsupported_architecture() -> None:
    assert bundled_mafft_command("Linux", machine="aarch64") is None


def test_run_alignment_with_bundled_mafft(tmp_path: Path) -> None:
    if __import__("platform").system() not in {"Darwin", "Linux", "Windows"}:
        pytest.skip("No bundled MAFFT for this platform")
    fasta = tmp_path / "pair.fa"
    fasta.write_text(">a\nACGTACGT\n>b\nACGTTCGT\n")
    out = run_alignment(fasta, backend="mafft", num_threads=1)
    assert out.is_file()
    assert out.read_text().count(">") == 2


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
