"""GFF3 Parent links round-trip and remount when a parent feature is renamed."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.sequence_features import (
    SequenceFeature,
    export_features_gff,
    load_features_from_gff,
    remount_parent_ids,
)


def test_export_gff_keeps_full_fasta_id_before_first_space(tmp_path: Path) -> None:
    """Export column 1 = full FASTA token (including # and /), not stripped at #."""
    from repeatviewer.backends.sequence_features import (
        SequenceFeature,
        export_features_gff,
        fasta_id_token,
    )

    seq_id = "Elbe4/Beta_vulgaris#Athila"
    assert fasta_id_token(f">{seq_id} description here") == seq_id
    fasta = tmp_path / "demo.fa"
    fasta.write_text(f">{seq_id} notes\n" + ("A" * 200) + "\n", encoding="utf-8")
    features = [
        SequenceFeature(start=1, end=50, name="LTR5", feature_type="LTR", color="#c0392b"),
    ]
    gff = export_features_gff(
        fasta,
        sequence_id=seq_id,
        sequence_length=200,
        features=features,
    )
    text = gff.read_text(encoding="utf-8")
    assert f"##sequence-region {seq_id} 1 200" in text
    assert text.splitlines()[2].startswith(f"{seq_id}\t")
    assert "Elbe4/Beta_vulgaris\t" not in text  # must not strip at #


def test_export_and_load_preserves_parent(tmp_path: Path) -> None:
    fasta = tmp_path / "demo.fa"
    fasta.write_text(">TE1\n" + ("A" * 200) + "\n", encoding="utf-8")
    features = [
        SequenceFeature(start=1, end=100, name="gene1", feature_type="gene"),
        SequenceFeature(
            start=10,
            end=40,
            name="exon1",
            feature_type="exon",
            parent="gene1",
        ),
    ]
    gff = export_features_gff(
        fasta,
        sequence_id="TE1",
        sequence_length=200,
        features=features,
    )
    text = gff.read_text(encoding="utf-8")
    assert "Parent=gene1" in text
    loaded = load_features_from_gff(gff, expected_id="TE1")
    by_name = {feat.name: feat for feat in loaded}
    assert by_name["exon1"].parent == "gene1"
    assert by_name["gene1"].parent == ""


def test_remount_parent_ids_updates_children() -> None:
    features = [
        SequenceFeature(start=1, end=100, name="geneA", feature_type="gene"),
        SequenceFeature(
            start=10,
            end=40,
            name="exon1",
            feature_type="exon",
            parent="geneA",
        ),
        SequenceFeature(
            start=50,
            end=80,
            name="exon2",
            feature_type="exon",
            parent="geneA",
        ),
        SequenceFeature(
            start=120,
            end=140,
            name="other",
            feature_type="ORF",
            parent="geneB",
        ),
    ]
    updated = remount_parent_ids(features, old_id="geneA", new_id="geneA_renamed")
    assert updated[0].name == "geneA"
    assert updated[1].parent == "geneA_renamed"
    assert updated[2].parent == "geneA_renamed"
    assert updated[3].parent == "geneB"
