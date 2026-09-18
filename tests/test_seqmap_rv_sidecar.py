"""Seqmap JSON naming and paired delete with RV.gff3."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.sequence_features import (
    SequenceFeature,
    companion_annotation_files_for_delete,
    is_sequence_map_hidden_file,
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
    features = [SequenceFeature(start=1, end=10, name="orf1", feature_type="ORF")]
    json_path, gff_path = save_annotations(
        fasta,
        sequence_id="TE1",
        sequence_length=50,
        features=features,
    )
    assert json_path.name == "Beetle1.seqmap_RV.json"
    assert json_path.is_file()
    assert gff_path is not None and gff_path.name == "Beetle1.RV.gff3"
    assert not legacy.exists()


def test_companion_delete_pairs_rv_gff3_with_seqmap_json(tmp_path: Path) -> None:
    gff = tmp_path / "Beetle1.RV.gff3"
    gff.write_text("##gff-version 3\n", encoding="utf-8")
    json_new = tmp_path / "Beetle1.seqmap_RV.json"
    json_new.write_text("{}\n", encoding="utf-8")
    json_old = tmp_path / "Beetle1.seqmap.json"
    json_old.write_text("{}\n", encoding="utf-8")
    companions = companion_annotation_files_for_delete(gff)
    names = {p.name for p in companions}
    assert names == {"Beetle1.seqmap_RV.json", "Beetle1.seqmap.json"}


def test_seqmap_json_files_are_hidden_from_file_list() -> None:
    assert is_sequence_map_hidden_file(Path("Beetle1.seqmap_RV.json"))
    assert is_sequence_map_hidden_file(Path("Beetle1.seqmap.json"))
