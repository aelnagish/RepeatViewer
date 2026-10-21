"""Pipeline path, output naming, and cleanup logic."""

from __future__ import annotations

from pathlib import Path

import pytest

from repeatviewer.backends.blast_gff import blast_gff_output_path
from repeatviewer.backends.blast_genome_plot import (
    blast_table_candidates_for_query,
    blast_txt_path_for_query,
)
from repeatviewer.backends.curation import (
    cleanup_blast_scratch_files,
    cleanup_flexidot_intermediates,
    effective_batch_workers,
    flexidot_output_stem,
)
from repeatviewer.backends.flexidot_gff import (
    persistent_flexidot_color_config_path,
    persistent_flexidot_legend_path,
)
from repeatviewer.browser_file_filter import is_tool_output_file


@pytest.mark.parametrize(
    ("workers", "sequences", "expected"),
    [
        (4, 1, 1),
        (1, 10, 1),
        (8, 3, 3),
    ],
)
def test_effective_batch_workers(workers: int, sequences: int, expected: int) -> None:
    assert effective_batch_workers(workers, sequences) == expected


def test_flexidot_output_stem_single(tmp_path: Path) -> None:
    fasta = tmp_path / "reads.fasta"
    fasta.write_text(">a\nA\n")
    assert flexidot_output_stem([fasta]) == "reads.fasta_flexidot"


def test_flexidot_output_stem_multi(tmp_path: Path) -> None:
    a = tmp_path / "a.fa"
    b = tmp_path / "b.fa"
    a.write_text(">a\nA\n")
    b.write_text(">b\nA\n")
    assert flexidot_output_stem([a, b]) == "a.fa_multi_flexidot"


def test_flexidot_sidecar_paths_use_output_stem_once(tmp_path: Path) -> None:
    stem = "reads.fasta_flexidot"
    assert persistent_flexidot_color_config_path(tmp_path, stem).name == (
        "reads.fasta_flexidot.gff.colors.config"
    )
    assert persistent_flexidot_legend_path(tmp_path, stem).name == (
        "reads.fasta_flexidot_legend.png"
    )


def test_blast_paths_next_to_query(tmp_path: Path) -> None:
    fasta = tmp_path / "reads.fasta"
    fasta.write_text(">a\nA\n")
    assert blast_txt_path_for_query(fasta) == tmp_path / "reads.fasta_bn.txt"
    assert blast_txt_path_for_query(fasta, "blastp") == tmp_path / "reads.fasta_bp.txt"
    candidates = blast_table_candidates_for_query(fasta)
    assert candidates[0] == tmp_path / "reads.fasta_bn.txt"
    assert tmp_path / "reads.fasta_b.txt" in candidates
    assert tmp_path / "reads.fasta_blast.txt" in candidates


def test_blast_paths_include_legacy_repeatviewer_folder(tmp_path: Path) -> None:
    fasta = tmp_path / "reads.fasta"
    fasta.write_text(">a\nA\n")
    legacy = tmp_path / "reads.repeatviewer" / "reads.fasta_blast.txt"
    legacy.parent.mkdir()
    legacy.write_text("q\ts\n")
    found = blast_table_candidates_for_query(fasta)
    assert legacy in found


def test_blast_gff_output_path() -> None:
    blast_txt = Path("/tmp/reads.fasta_b.txt")
    assert blast_gff_output_path(blast_txt).name == "reads.fasta_b.gff"
    legacy = Path("/tmp/reads.fasta_blast.txt")
    assert blast_gff_output_path(legacy).name == "reads.fasta_blast.gff"


def test_kept_input_name_used_by_align_and_combine(tmp_path: Path) -> None:
    from repeatviewer.backends.sequence_aligner import (
        add_to_msa_output_stem,
        alignment_output_stem,
        combined_alignment_output_stem,
        _next_available_fasta,
    )
    from repeatviewer.backends.sequence_compose import default_group_name
    from repeatviewer.tool_output_paths import kept_input_name

    a = tmp_path / "reads.fa"
    b = tmp_path / "other.fasta"
    a.write_text(">a\nA\n")
    b.write_text(">b\nA\n")
    assert kept_input_name(a) == "reads.fa"
    assert alignment_output_stem(a, engine="mafft") == "reads.fa_mafft"
    assert alignment_output_stem(a, engine="muscle") == "reads.fa_muscle"
    assert combined_alignment_output_stem([a, b], engine="mafft") == "reads.fa_cmafft"
    assert combined_alignment_output_stem([a, b], engine="muscle") == "reads.fa_cmuscle"
    assert add_to_msa_output_stem(a, engine="mafft") == "reads.fa_amafft"
    assert add_to_msa_output_stem(a, engine="muscle") == "reads.fa_amuscle"
    first = _next_available_fasta(tmp_path, "reads.fa_cmafft")
    assert first.name == "reads.fa_cmafft.fa"
    first.write_text(">x\nA\n")
    second = _next_available_fasta(tmp_path, "reads.fa_cmafft")
    assert second.name == "reads.fa_cmafft_1.fa"
    assert default_group_name([a, b]) == "reads.fa_and_other.fasta_comb"


def test_cleanup_flexidot_intermediates(tmp_path: Path) -> None:
    keep = tmp_path / "reads.fasta_flexidot.png"
    keep.write_text("png")
    scratch_dir = tmp_path / "flexidot_output-1"
    scratch_dir.mkdir()
    (scratch_dir / "plot.png").write_text("x")
    temp_gff = tmp_path / ".repeatviewer_flexidot_reads.fasta_flexidot_1.gff3"
    temp_gff.write_text("##gff")
    cleanup_flexidot_intermediates(tmp_path, [keep])
    assert keep.is_file()
    assert not scratch_dir.exists()
    assert not temp_gff.exists()


def test_cleanup_blast_scratch_files(tmp_path: Path) -> None:
    fasta = tmp_path / "reads.fasta"
    fasta.write_text(">a\nA\n")
    scratch = tmp_path / "reads.fasta_bn_no_header.b"
    scratch.write_text("x")
    legacy = tmp_path / "reads.fasta_b_no_header.b"
    legacy.write_text("x")
    cleanup_blast_scratch_files(tmp_path, fasta)
    assert not scratch.exists()
    assert not legacy.exists()


@pytest.mark.parametrize(
    "name",
    [
        "reads.fasta_bn.txt",
        "reads.fasta_bn.gff",
        "reads.fasta_b.txt",
        "reads.fasta_b.gff",
        "reads.fasta_flexidot.png",
        "reads.fasta_flexidot.gff.colors.config",
        "reads.fasta_b.txt.fa",
        "reads.fa.tsd_search.txt",
        "reads.fasta_blast.txt",
        "reads.fasta_blast.gff",
    ],
)
def test_is_tool_output_file(name: str, tmp_path: Path) -> None:
    path = tmp_path / name
    path.write_text("x")
    assert is_tool_output_file(path)
