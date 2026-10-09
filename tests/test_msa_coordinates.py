"""MSA gap strip / feature projection helpers."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from repeatviewer.backends.msa_coordinates import (
    gapped_to_ungapped_1based,
    project_feature_to_gapped,
    project_feature_to_ungapped,
    project_features_to_gapped,
    sequence_has_gaps,
    ungap_sequence,
    ungapped_to_gapped_1based,
)
from repeatviewer.backends.sequence_features import SequenceFeature


def test_ungap_sequence_builds_keep_map() -> None:
    ungapped, keep = ungap_sequence("AC--GT-A")
    assert ungapped == "ACGTA"
    assert keep == [0, 1, 4, 5, 7]
    assert sequence_has_gaps("AC--GT")
    assert not sequence_has_gaps("ACGT")
    plain, none_keep = ungap_sequence("ACGT")
    assert plain == "ACGT"
    assert none_keep is None


def test_round_trip_feature_coords() -> None:
    # Gapped: A C - - G T - A
    # Index:  1 2 3 4 5 6 7 8
    # Ungap:  A C G T A  → positions 1..5
    _, keep = ungap_sequence("AC--GT-A")
    assert keep is not None
    assert ungapped_to_gapped_1based(1, keep) == 1
    assert ungapped_to_gapped_1based(3, keep) == 5  # G
    assert ungapped_to_gapped_1based(5, keep) == 8  # A
    assert gapped_to_ungapped_1based(5, keep) == 3
    assert gapped_to_ungapped_1based(6, keep) == 4

    feat_u = SequenceFeature(start=2, end=4, name="x", feature_type="LTR")
    feat_g = project_feature_to_gapped(feat_u, keep)
    assert feat_g.start == 2
    assert feat_g.end == 6  # T at ungapped 4 → gapped 6
    back = project_feature_to_ungapped(feat_g, keep)
    assert back is not None
    assert back.start == 2
    assert back.end == 4


def test_annotate_one_record_projects_msa_features(tmp_path: Path) -> None:
    from repeatviewer.backends.auto_annotation import AutoAnnotationMode
    from repeatviewer.backends.te_annotation_batch import (
        TeAnnotationOptions,
        _annotate_one_record_fasta,
    )
    from repeatviewer.backends.auto_annotation_settings import (
        load_auto_annotation_settings,
    )
    from repeatviewer.backends.similarity_annotation_settings import (
        load_similarity_annotation_settings,
    )

    gapped = "ATG---CCC"
    fasta = tmp_path / "msa_row.fa"
    SeqIO.write(
        SeqRecord(Seq(gapped), id="s1", name="s1", description="s1"),
        str(fasta),
        "fasta",
    )

    def _fake_auto(path, modes, **kwargs):
        del modes, kwargs
        # Detector saw ungapped ATGCCC (6 bp) — report ORF covering 1..6.
        text = Path(path).read_text(encoding="utf-8")
        assert "ATG---CCC" not in text
        assert "ATGCCC" in text

        class _R:
            features = (
                SequenceFeature(start=1, end=6, name="orf1", feature_type="ORF"),
            )

        return _R()

    with patch(
        "repeatviewer.backends.te_annotation_batch.run_auto_annotations",
        side_effect=_fake_auto,
    ):
        out = _annotate_one_record_fasta(
            fasta,
            TeAnnotationOptions(auto_modes=(AutoAnnotationMode.ORF,)),
            auto_cfg=load_auto_annotation_settings(),
            sim_cfg=load_similarity_annotation_settings(),
            use_blast=False,
            base_features=(),
            cancel=None,
            progress=None,
            label="s1",
        )
    assert len(out) == 1
    # Ungapped 1..6 → gapped ATG---CCC columns 1..9
    assert out[0].start == 1
    assert out[0].end == 9


def test_project_features_list_noop_without_gaps() -> None:
    feats = [SequenceFeature(start=3, end=5, name="a", feature_type="LTR")]
    assert project_features_to_gapped(feats, None)[0].start == 3
