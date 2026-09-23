"""Flexidot companion GFF resolution and sidecar finalization."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.curation import finalize_flexidot_outputs
from repeatviewer.backends.flexidot_gff import (
    prepare_flexidot_gff_inputs,
    resolve_gff_for_fasta,
)
from repeatviewer.backends.gff_features import find_companion_gff


def test_find_companion_gff_prefers_rv_gff3(tmp_path: Path) -> None:
    fasta = tmp_path / "reads.fa"
    fasta.write_text(">a\nACGT\n")
    blast_gff = tmp_path / "reads.fa_b.gff"
    blast_gff.write_text("##gff-version 3\n")
    companion = tmp_path / "reads.RV.gff3"
    companion.write_text("##gff-version 3\n")
    assert find_companion_gff(fasta) == companion


def test_find_companion_gff_ignores_blast_hit_table(tmp_path: Path) -> None:
    fasta = tmp_path / "reads.fa"
    fasta.write_text(">a\nACGT\n")
    (tmp_path / "reads.fa_b.gff").write_text("##gff-version 3\n")
    (tmp_path / "other_thing.gff").write_text("##gff-version 3\n")
    assert find_companion_gff(fasta) is None
    assert resolve_gff_for_fasta(fasta) is None


def test_prepare_flexidot_gff_returns_none_without_companion(tmp_path: Path) -> None:
    fasta = tmp_path / "reads.fa"
    fasta.write_text(">a\nACGT\n")
    (tmp_path / "reads.fa_b.gff").write_text("##gff-version 3\n")
    assert prepare_flexidot_gff_inputs([fasta], tmp_path, output_stem="reads.fa_flexidot") is None


def test_finalize_flexidot_skips_sidecars_without_gff(tmp_path: Path) -> None:
    fasta = tmp_path / "reads.fa"
    fasta.write_text(">a\nACGT\n")
    plot = tmp_path / "raw.pdf"
    plot.write_bytes(b"%PDF")
    before: set[Path] = set()
    # Snapshot before would exclude raw.pdf if taken earlier; simulate post-run
    # by passing empty before so raw.pdf is the new file.
    outputs = finalize_flexidot_outputs(
        tmp_path,
        [fasta],
        before,
        ext="pdf",
        use_gff=True,
        gff_source_fastas=[fasta],
    )
    assert len(outputs) == 1
    assert Path(outputs[0]).name == "reads.fa_flexidot.pdf"
    assert not (tmp_path / "reads.fa_flexidot.gff.colors.config").exists()
    assert not (tmp_path / "reads.fa_flexidot_legend.png").exists()
