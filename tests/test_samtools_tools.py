"""samtools Sort / Index / Convert backend helpers."""

from __future__ import annotations

from pathlib import Path

import pytest

from repeatviewer.backends import samtools_tools as st


def _write_mini_sam(path: Path) -> Path:
    path.write_text(
        "@HD\tVN:1.6\tSO:unsorted\n"
        "@SQ\tSN:chr1\tLN:1000\n"
        "r1\t0\tchr1\t1\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n"
        "r2\t0\tchr1\t50\t60\t10M\t*\t0\t0\tACGTACGTAC\t*\n",
        encoding="utf-8",
    )
    return path


def test_default_output_paths(tmp_path: Path) -> None:
    sam = tmp_path / "reads_minimap.sam"
    assert st.default_sorted_bam_path(sam).name == "reads_minimap.sorted.bam"
    assert st.default_converted_bam_path(sam).name == "reads_minimap.bam"
    assert st.default_converted_path(sam, "cram").name == "reads_minimap.cram"
    bam = tmp_path / "x.bam"
    assert st.default_converted_path(bam, "sam").name == "x.sam"
    assert st.default_converted_path(bam, "bam").name == "x.bam"


def test_convert_skips_same_format(tmp_path: Path) -> None:
    sam = _write_mini_sam(tmp_path / "reads.sam")
    assert st.convert_alignment(sam, output_format="sam") == []
    with pytest.raises(st.SamtoolsError, match="already SAM"):
        st.run_samtools_operation("convert", [sam], convert_format="sam")


def test_is_alignment_path() -> None:
    assert st.is_alignment_path(Path("a.sam"))
    assert st.is_alignment_path(Path("a.bam"))
    assert st.is_alignment_path(Path("a.cram"))
    assert not st.is_alignment_path(Path("a.fa"))


def test_alignment_index_companions_and_renames(tmp_path: Path) -> None:
    bam = tmp_path / "reads.sorted.bam"
    bam.write_bytes(b"BAM")
    bai = tmp_path / "reads.sorted.bam.bai"
    bai.write_bytes(b"BAI")
    assert st.is_alignment_index_file(bai)
    assert not st.is_alignment_index_file(bam)
    assert st.alignment_index_companions(bam) == [bai]
    new_bam = tmp_path / "other.bam"
    assert st.alignment_index_renames(bam, new_bam) == [
        (bai, tmp_path / "other.bam.bai")
    ]
    # SAM has no index companions.
    sam = tmp_path / "reads.sam"
    sam.write_text("@HD\tVN:1.6\n", encoding="utf-8")
    assert st.alignment_index_companions(sam) == []
    assert st.alignment_index_renames(sam, tmp_path / "x.bam") == []


def test_sort_index_convert_roundtrip(tmp_path: Path) -> None:
    pysam = pytest.importorskip("pysam")
    del pysam
    sam = _write_mini_sam(tmp_path / "reads.sam")
    converted = st.convert_to_bam(sam)
    assert converted and converted[0].is_file()
    to_sam = st.convert_alignment(converted[0], output_format="sam")
    assert to_sam and to_sam[0].is_file() and to_sam[0].suffix == ".sam"
    sorted_files = st.sort_alignment(sam, index_after=True)
    bam = sorted_files[0]
    assert bam.is_file()
    assert any(p.suffix == ".bai" or p.name.endswith(".bai") for p in sorted_files)
    indexed = st.index_alignment(bam)
    assert indexed and indexed[0].is_file()


def test_run_samtools_operation_sort(tmp_path: Path) -> None:
    pytest.importorskip("pysam")
    sam = _write_mini_sam(tmp_path / "q.sam")
    result = st.run_samtools_operation(
        "sort", [sam], threads=1, index_after_sort=True
    )
    assert result["operation"] == "sort"
    assert any(Path(p).suffix == ".bam" for p in result["files"])
