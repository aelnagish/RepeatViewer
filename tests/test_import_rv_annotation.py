"""Import external annotations into the JSON sidecar (GFF export is separate)."""

from __future__ import annotations

from pathlib import Path

import pytest

from repeatviewer.backends.sequence_features import (
    features_gff_path,
    import_annotation_as_rv_gff3,
    load_features_by_id_from_json,
    load_sequence_map_data,
    resolve_features_gff_path,
    sidecar_path,
)


def _write_fasta(path: Path, seq_id: str = "TE1", sequence: str = "ACGT" * 50) -> Path:
    path.write_text(f">{seq_id}\n{sequence}\n", encoding="utf-8")
    return path


def test_features_gff_path_uses_rv_suffix(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "demo.fa")
    assert features_gff_path(fasta).name == "demo.RV.gff3"


def test_import_gff3_writes_json_not_gff(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "demo.fa", "TE1", "A" * 200)
    src = tmp_path / "external.gff3"
    src.write_text(
        "##gff-version 3\n"
        "TE1\tother\tLTR\t1\t50\t.\t+\t.\tID=ltr5;Name=LTR5\n"
        "TE1\tother\tLTR\t151\t200\t.\t-\t.\tID=ltr3;Name=LTR3\n",
        encoding="utf-8",
    )

    out, count = import_annotation_as_rv_gff3(fasta, src, merge_existing=False)
    assert count == 2
    assert out == sidecar_path(fasta)
    assert out.is_file()
    assert not features_gff_path(fasta).exists()

    feats = load_features_by_id_from_json(fasta, ["TE1"]).get("TE1") or []
    assert len(feats) == 2
    assert {f.name for f in feats} == {"LTR5", "LTR3"}
    assert {f.source for f in feats} == {"external.gff3"}


def test_import_rejects_mismatched_seqid_for_single_fasta(tmp_path: Path) -> None:
    """Do not silently attach Elbe4 annotations onto an L1 (or other) FASTA."""
    fasta = _write_fasta(tmp_path / "L1.fa", "L1_Belline1", "ACGT" * 25)
    gff = tmp_path / "Elbe4.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "Elbe4_Beta_vulgaris\tRepeatViewer\tLTR\t1\t50\t.\t+\t.\tName=LTR5\n",
        encoding="utf-8",
    )
    with pytest.raises(ValueError, match="matched any sequence id"):
        import_annotation_as_rv_gff3(fasta, gff, merge_existing=False)


def test_import_matches_fasta_id_with_hash_and_slash(tmp_path: Path) -> None:
    seq_id = "TE/Class_I#LINE"
    fasta = _write_fasta(tmp_path / "te.fa", seq_id, "A" * 100)
    src = tmp_path / "te.gff3"
    src.write_text(
        f"##gff-version 3\n"
        f"{seq_id}\tx\tORF\t1\t30\t.\t+\t.\tName=orf1\n",
        encoding="utf-8",
    )
    out, count = import_annotation_as_rv_gff3(fasta, src, merge_existing=False)
    assert count == 1
    feats = load_features_by_id_from_json(fasta, [seq_id]).get(seq_id) or []
    assert len(feats) == 1
    assert feats[0].name == "orf1"
    assert out == sidecar_path(fasta)


def test_import_merges_with_existing_json(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "demo.fa", "TE1", "N" * 100)
    first = tmp_path / "a.gff3"
    first.write_text(
        "TE1\tx\tORF\t1\t30\t.\t+\t.\tName=orf1\n",
        encoding="utf-8",
    )
    import_annotation_as_rv_gff3(fasta, first, merge_existing=False)

    second = tmp_path / "b.gff3"
    second.write_text(
        "TE1\tx\tLTR\t70\t100\t.\t+\t.\tName=ltr\n",
        encoding="utf-8",
    )
    out, count = import_annotation_as_rv_gff3(fasta, second, merge_existing=True)
    assert count == 2
    assert out == sidecar_path(fasta)
    feats = load_features_by_id_from_json(fasta, ["TE1"]).get("TE1") or []
    assert {f.name for f in feats} == {"orf1", "ltr"}
    assert not features_gff_path(fasta).exists()


def test_resolve_still_finds_legacy_features_gff3(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "old.fa", "X", "ACGT")
    legacy = tmp_path / "old.features.gff3"
    legacy.write_text(
        "X\tRepeatViewer\tothers\t1\t4\t.\t+\t.\tName=a;type=others\n",
        encoding="utf-8",
    )
    found = resolve_features_gff_path(fasta)
    assert found == legacy


def test_import_multi_sequence_matches_by_seqid(tmp_path: Path) -> None:
    fasta = tmp_path / "lib.fa"
    fasta.write_text(
        ">TE1\n" + ("A" * 100) + "\n"
        ">TE2\n" + ("C" * 80) + "\n"
        ">TE3\n" + ("G" * 60) + "\n",
        encoding="utf-8",
    )
    src = tmp_path / "external.gff3"
    src.write_text(
        "##gff-version 3\n"
        "##sequence-region TE1 1 100\n"
        "##sequence-region TE2 1 80\n"
        "TE1\to\tLTR\t1\t20\t.\t+\t.\tName=a\n"
        "TE2\to\tORF\t5\t40\t.\t+\t.\tName=b\n"
        "OTHER\to\tLTR\t1\t10\t.\t+\t.\tName=orphan\n",
        encoding="utf-8",
    )

    result = import_annotation_as_rv_gff3(fasta, src, merge_existing=False)
    assert result.feature_count == 2
    assert result.sequence_count == 2
    assert result.json_path == sidecar_path(fasta)
    assert not features_gff_path(fasta).exists()

    data = load_sequence_map_data(fasta)
    assert data.is_multi_sequence
    assert len(data.record_ids) == 3
    assert data.sequence_id in {"TE1", "TE2"}  # first with features
    data.switch_to("TE2")
    assert data.sequence_id == "TE2"
    assert len(data.features) == 1
    assert data.features[0].name == "b"


def test_import_rejects_unsupported_suffix(tmp_path: Path) -> None:
    fasta = _write_fasta(tmp_path / "demo.fa")
    bad = tmp_path / "notes.txt"
    bad.write_text("nope\n", encoding="utf-8")
    with pytest.raises(ValueError, match="Unsupported"):
        import_annotation_as_rv_gff3(fasta, bad)
