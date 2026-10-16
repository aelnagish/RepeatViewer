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


def test_run_auto_annotations_projects_msa_features(tmp_path: Path) -> None:
    from repeatviewer.backends.auto_annotation import (
        AutoAnnotationMode,
        AutoAnnotationResult,
        run_auto_annotations,
    )

    gapped = "ATG---CCC"
    fasta = tmp_path / "msa_row.fa"
    SeqIO.write(
        SeqRecord(Seq(gapped), id="s1", name="s1", description="s1"),
        str(fasta),
        "fasta",
    )

    def _fake_one(path, mode, **kwargs):
        del mode, kwargs
        # Detector must see ungapped ATGCCC (6 bp).
        text = Path(path).read_text(encoding="utf-8")
        assert "ATG---CCC" not in text
        assert "ATGCCC" in text
        return AutoAnnotationResult(
            features=(
                SequenceFeature(start=1, end=6, name="orf1", feature_type="ORF"),
            ),
            mode=AutoAnnotationMode.ORF,
            message="orf ok",
        )

    with patch(
        "repeatviewer.backends.auto_annotation.run_auto_annotation",
        side_effect=_fake_one,
    ):
        result = run_auto_annotations(
            fasta,
            (AutoAnnotationMode.ORF,),
            base_features=(),
        )
    assert len(result.features) == 1
    # Ungapped 1..6 → gapped ATG---CCC columns 1..9
    assert result.features[0].start == 1
    assert result.features[0].end == 9


def test_project_features_list_noop_without_gaps() -> None:
    feats = [SequenceFeature(start=3, end=5, name="a", feature_type="LTR")]
    assert project_features_to_gapped(feats, None)[0].start == 3


def test_column_residue_frequencies_and_threshold() -> None:
    from repeatviewer.backends.msa_coordinates import (
        MSA_COLOR_ALL_BACKGROUNDS,
        MSA_COLOR_BG_WHEN_CONSERVED,
        MSA_COLOR_DIFF_FROM_SELECTED,
        MSA_COLOR_LETTERS_ONLY,
        MSA_COLOR_MODE_CHOICES,
        MSA_COLOR_NO_BG_WHEN_CONSERVED,
        column_conservation_mask,
        column_residue_frequencies,
        letter_is_conserved,
        msa_residue_fill_background,
        normalize_msa_color_mode,
        normalize_msa_divergence_threshold,
        residue_differs_from_reference,
    )

    freqs = column_residue_frequencies(["AA-", "AA-", "AC-"])
    assert freqs[0] == {"A": 1.0}
    assert abs(freqs[1]["A"] - 2 / 3) < 1e-9
    assert abs(freqs[1]["C"] - 1 / 3) < 1e-9
    assert freqs[2] == {}
    assert letter_is_conserved(freqs[1], "A", threshold=0.4)
    assert not letter_is_conserved(freqs[1], "C", threshold=0.4)
    assert column_conservation_mask(
        ["AA-", "AA-", "AC-"], threshold=0.4
    ) == [True, True, False]
    # Unanimous columns at threshold 1.0.
    assert column_conservation_mask(
        ["AAT-", "ACT-", "A-T-"],
        threshold=1.0,
    ) == [True, False, True, False]

    assert normalize_msa_color_mode(None) == MSA_COLOR_ALL_BACKGROUNDS
    assert normalize_msa_divergence_threshold(0.4) == 0.4

    freqs_mix = {"A": 0.7, "T": 0.3}
    assert msa_residue_fill_background(
        MSA_COLOR_ALL_BACKGROUNDS,
        residue="A",
        column_freqs=freqs_mix,
        threshold=0.4,
    )
    assert not msa_residue_fill_background(
        MSA_COLOR_LETTERS_ONLY,
        residue="A",
        column_freqs=freqs_mix,
        threshold=0.4,
    )
    assert msa_residue_fill_background(
        MSA_COLOR_BG_WHEN_CONSERVED,
        residue="A",
        column_freqs=freqs_mix,
        threshold=0.4,
    )
    assert not msa_residue_fill_background(
        MSA_COLOR_BG_WHEN_CONSERVED,
        residue="T",
        column_freqs=freqs_mix,
        threshold=0.4,
    )
    assert msa_residue_fill_background(
        MSA_COLOR_NO_BG_WHEN_CONSERVED,
        residue="T",
        column_freqs=freqs_mix,
        threshold=0.4,
    )
    assert not msa_residue_fill_background(
        mode=MSA_COLOR_ALL_BACKGROUNDS,
        residue="-",
        column_freqs=freqs_mix,
        threshold=0.4,
    )
    # Retired conservation_grey settings migrate to Color all.
    assert normalize_msa_color_mode("conservation_grey") == (
        MSA_COLOR_ALL_BACKGROUNDS
    )
    assert normalize_msa_color_mode(MSA_COLOR_DIFF_FROM_SELECTED) == (
        MSA_COLOR_DIFF_FROM_SELECTED
    )
    assert dict(MSA_COLOR_MODE_CHOICES)[MSA_COLOR_DIFF_FROM_SELECTED] == (
        "Color vs reference"
    )

    assert residue_differs_from_reference("T", "A")
    assert not residue_differs_from_reference("A", "a")
    assert residue_differs_from_reference("-", "A")
    assert residue_differs_from_reference("A", "-")
    assert not residue_differs_from_reference("-", "-")
    assert not residue_differs_from_reference(".", "-")
    assert not residue_differs_from_reference("A", None)
    # Standard row: fully coloured.
    assert msa_residue_fill_background(
        MSA_COLOR_DIFF_FROM_SELECTED,
        residue="A",
        reference_residue="A",
        is_reference_row=True,
    )
    assert msa_residue_fill_background(
        MSA_COLOR_DIFF_FROM_SELECTED,
        residue="T",
        reference_residue="A",
        is_reference_row=True,
    )
    assert not msa_residue_fill_background(
        MSA_COLOR_DIFF_FROM_SELECTED,
        residue="-",
        reference_residue="A",
        is_reference_row=True,
    )
    # Other rows: only differences (including gap vs standard letter).
    assert msa_residue_fill_background(
        MSA_COLOR_DIFF_FROM_SELECTED,
        residue="T",
        reference_residue="A",
    )
    assert not msa_residue_fill_background(
        MSA_COLOR_DIFF_FROM_SELECTED,
        residue="A",
        reference_residue="A",
    )
    assert msa_residue_fill_background(
        MSA_COLOR_DIFF_FROM_SELECTED,
        residue="-",
        reference_residue="A",
    )
    assert not msa_residue_fill_background(
        MSA_COLOR_DIFF_FROM_SELECTED,
        residue="-",
        reference_residue="-",
    )
    assert not msa_residue_fill_background(
        MSA_COLOR_DIFF_FROM_SELECTED,
        residue="T",
        reference_residue=None,
    )


def test_fill_rows_region_and_normalize_edit_residue() -> None:
    from repeatviewer.backends.msa_coordinates import (
        fill_rows_region,
        mask_rows_region,
        normalize_msa_edit_residue,
    )

    assert normalize_msa_edit_residue("a", protein=False) == "A"
    assert normalize_msa_edit_residue("-", protein=False) == "-"
    assert normalize_msa_edit_residue("Q", protein=False) is None
    assert normalize_msa_edit_residue("q", protein=True) == "Q"
    assert normalize_msa_edit_residue("AC", protein=False) is None

    seqs = ["ACGTAC", "ACGTAC"]
    assert fill_rows_region(seqs, 2, 4, fill="N") == ["ANNNAC", "ANNNAC"]
    assert fill_rows_region(seqs, 2, 4, {1}, fill="T") == ["ACGTAC", "ATTTAC"]
    assert mask_rows_region(seqs, 3, 3) == ["AC-TAC", "AC-TAC"]
