"""Seqmap JSON naming and paired delete / rename (JSON only; not GFF export)."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.sequence_features import (
    SequenceFeature,
    companion_annotation_files_for_delete,
    companion_annotation_renames_for_fasta,
    features_gff_path,
    is_sequence_map_hidden_file,
    load_sequence_map_data,
    resolve_sidecar_path,
    save_annotations,
    sidecar_path,
)


def test_sidecar_path_uses_seqmap_rv_suffix(tmp_path: Path) -> None:
    fasta = tmp_path / "Beetle1.fasta"
    fasta.write_text(">TE1\nACGT\n", encoding="utf-8")
    assert sidecar_path(fasta).name == "Beetle1.seqmap_RV.json"


def test_legacy_seqmap_json_migrates_to_seqmap_rv(tmp_path: Path) -> None:
    fasta = tmp_path / "Beetle1.fasta"
    fasta.write_text(">TE1\nACGT\n", encoding="utf-8")
    legacy = tmp_path / "Beetle1.seqmap.json"
    legacy.write_text(
        '{"sequence_id":"TE1","features":[{"start":1,"end":4,"name":"x"}]}\n',
        encoding="utf-8",
    )
    resolved = resolve_sidecar_path(fasta)
    assert resolved is not None
    assert resolved.name == "Beetle1.seqmap_RV.json"
    assert resolved.is_file()
    assert not legacy.exists()


def test_save_annotations_writes_seqmap_rv_and_cleans_legacy(tmp_path: Path) -> None:
    fasta = tmp_path / "Beetle1.fasta"
    fasta.write_text(">TE1\n" + ("A" * 50) + "\n", encoding="utf-8")
    legacy = tmp_path / "Beetle1.seqmap.json"
    legacy.write_text('{"sequence_id":"TE1","features":[]}\n', encoding="utf-8")
    features = [
        SequenceFeature(
            start=1, end=10, name="orf1", feature_type="ORF", color="#27ae60"
        )
    ]
    json_path, gff_path = save_annotations(
        fasta,
        sequence_id="TE1",
        sequence_length=50,
        features=features,
    )
    assert json_path.name == "Beetle1.seqmap_RV.json"
    assert json_path.is_file()
    assert gff_path is None
    assert not (tmp_path / "Beetle1.RV.gff3").exists()
    assert not legacy.exists()


def test_companion_delete_fasta_takes_json_not_gff(tmp_path: Path) -> None:
    fasta = tmp_path / "Beetle1.fasta"
    fasta.write_text(">TE1\nACGT\n", encoding="utf-8")
    gff = tmp_path / "Beetle1.RV.gff3"
    gff.write_text("##gff-version 3\n", encoding="utf-8")
    json_new = tmp_path / "Beetle1.seqmap_RV.json"
    json_new.write_text("{}\n", encoding="utf-8")
    companions = companion_annotation_files_for_delete(fasta)
    names = {p.name for p in companions}
    assert names == {"Beetle1.seqmap_RV.json"}


def test_companion_delete_gff_does_not_take_json(tmp_path: Path) -> None:
    gff = tmp_path / "Beetle1.RV.gff3"
    gff.write_text("##gff-version 3\n", encoding="utf-8")
    json_new = tmp_path / "Beetle1.seqmap_RV.json"
    json_new.write_text("{}\n", encoding="utf-8")
    assert companion_annotation_files_for_delete(gff) == []


def test_seqmap_json_files_are_hidden_from_file_list() -> None:
    assert is_sequence_map_hidden_file(Path("Beetle1.seqmap_RV.json"))
    assert is_sequence_map_hidden_file(Path("Beetle1.seqmap.json"))


def test_empty_json_sidecar_does_not_restore_gff_features(tmp_path: Path) -> None:
    """Deleting all features must stick: empty JSON wins over leftover .RV.gff3."""
    fasta = tmp_path / "Elbe4.fa"
    fasta.write_text(">Elbe4_Beta_vulgaris\n" + ("ATGC" * 2500) + "\n", encoding="utf-8")
    gff = features_gff_path(fasta)
    gff.write_text(
        "##gff-version 3\n"
        "##sequence-region Elbe4_Beta_vulgaris 1 10000\n"
        "Elbe4_Beta_vulgaris\tRepeatViewer\tLTR\t1\t1404\t.\t+\t.\t"
        "ID=Elbe4_5-LTR;Name=Elbe4_5-LTR;type=LTR\n"
        "Elbe4_Beta_vulgaris\tRepeatViewer\tORF\t1962\t6844\t.\t+\t.\t"
        "ID=Elbe4_gag-pol-ORF;Name=Elbe4_gag-pol-ORF;type=ORF\n",
        encoding="utf-8",
    )
    data = load_sequence_map_data(fasta)
    assert len(data.features) == 2

    save_annotations(
        fasta,
        sequence_id=data.sequence_id,
        sequence_length=data.length,
        features=[],
    )
    assert sidecar_path(fasta).is_file()
    assert gff.is_file()

    reloaded = load_sequence_map_data(fasta)
    assert reloaded.features == []


def test_companion_renames_follow_fasta_json_only(tmp_path: Path) -> None:
    old_fasta = tmp_path / "Beetle1.fasta"
    new_fasta = tmp_path / "Beetle1_copy.fasta"
    old_fasta.write_text(">TE1\nACGT\n", encoding="utf-8")
    gff = tmp_path / "Beetle1.RV.gff3"
    gff.write_text("##gff-version 3\n", encoding="utf-8")
    json_path = tmp_path / "Beetle1.seqmap_RV.json"
    json_path.write_text("{}\n", encoding="utf-8")

    renames = companion_annotation_renames_for_fasta(old_fasta, new_fasta)
    by_src = {src.name: dst.name for src, dst in renames}
    assert "Beetle1.RV.gff3" not in by_src
    assert by_src["Beetle1.seqmap_RV.json"] == "Beetle1_copy.seqmap_RV.json"
    assert gff.is_file()
