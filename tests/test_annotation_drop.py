"""Annotation files dropped on a FASTA row in the Files panel."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.sequence_features import (
    annotation_import_paths,
    import_annotation_as_rv_gff3,
    is_importable_annotation_path,
)


def test_importable_suffixes() -> None:
    for name in ("a.gff", "a.GFF3", "a.gtf", "a.bed", "a.RV.gff3"):
        assert is_importable_annotation_path(Path(name))
    for name in ("a.gff.gz", "a.fasta", "a.txt", "a.png", "a"):
        assert not is_importable_annotation_path(Path(name))


def test_annotation_import_paths_filters_and_dedupes(tmp_path: Path) -> None:
    gff = tmp_path / "x.gff3"
    bed = tmp_path / "y.bed"
    for path in (gff, bed, tmp_path / "z.gff.gz", tmp_path / "pic.png"):
        path.write_text("")
    picked = annotation_import_paths(
        [gff, tmp_path / "pic.png", tmp_path / "z.gff.gz", gff, bed, tmp_path / "gone.gff"]
    )
    assert [p.name for p in picked] == ["x.gff3", "y.bed"]


def test_dropped_gff_with_other_name_matches_by_id(tmp_path: Path) -> None:
    fasta = tmp_path / "genome.fasta"
    fasta.write_text(">chr1 desc\n" + "ACGT" * 8 + "\n>chr2\n" + "GC" * 8 + "\n")
    gff = tmp_path / "completely_different.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "chr1\tsrc\tgene\t2\t10\t.\t+\t.\tID=g1\n"
        "chr2\tsrc\trepeat\t3\t8\t.\t-\t.\tID=r1\n"
        "chrX\tsrc\tgene\t1\t5\t.\t+\t.\tID=skip\n"
    )
    first = import_annotation_as_rv_gff3(fasta, gff)
    assert (first.feature_count, first.sequence_count) == (2, 2)
    # Dropping the same file again must not double the features.
    again = import_annotation_as_rv_gff3(fasta, gff)
    assert again.feature_count == 2
    assert gff.is_file()
