"""Bundled MAFFT selection and alignment helpers."""

from __future__ import annotations

from pathlib import Path
import pytest

from teviewer.backends.sequence_aligner import (
    _fasta_sequence_count,
    _mafft_arguments,
    bundled_mafft_command,
    run_alignment,
)
from teviewer.job_options import (
    MAFFT_SEQUENCE_NUCLEOTIDE,
    MAFFT_STRATEGY_CUSTOM,
    MAFFT_STRATEGY_LINSI,
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
    out = run_alignment(
        fasta,
        num_threads=1,
        strategy=MAFFT_STRATEGY_LINSI,
        sequence_type=MAFFT_SEQUENCE_NUCLEOTIDE,
    )
    assert out.is_file()
    assert out.read_text().count(">") == 2


def test_custom_mafft_arguments_include_advanced_options() -> None:
    args = _mafft_arguments(
        num_threads=7,
        strategy=MAFFT_STRATEGY_CUSTOM,
        sequence_type=MAFFT_SEQUENCE_NUCLEOTIDE,
        retree=4,
        max_iterate=25,
        gap_open_penalty=2.25,
        offset=0.2,
        reorder=True,
        adjust_direction=True,
    )
    assert args == [
        "--quiet",
        "--retree", "4",
        "--maxiterate", "25",
        "--nuc",
        "--thread", "7",
        "--op", "2.25",
        "--ep", "0.2",
        "--reorder",
        "--adjustdirectionaccurately",
    ]
