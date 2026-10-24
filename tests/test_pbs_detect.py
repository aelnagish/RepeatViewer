"""PBS motif path + LTR-pair selection for Feature detection."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.auto_annotation import AutoAnnotationMode, _stamp_track
from repeatviewer.backends.pbs_detect import (
    PBSCandidate,
    _extract_ltrs,
    default_motif_fasta,
    motif_database_label,
    pbs_feature_from_candidate,
    pbs_motif_database_available,
    pbs_motif_source,
)
from repeatviewer.backends.sequence_features import (
    FEATURE_TYPE_COLORS,
    SequenceFeature,
    feature_has_hover_meta,
    feature_hover_detail_text,
)


def test_bundled_pbs_motif_fasta_exists() -> None:
    path = default_motif_fasta()
    assert path.is_file()
    assert pbs_motif_database_available(path)
    assert path.name == "PBS_unique_motifs_Neumann2019.fasta"
    assert motif_database_label(path) == "Neumann2019"
    assert pbs_motif_source(path) == "PBS motif: Neumann2019"


def test_pbs_feature_shows_identity_coverage_and_db() -> None:
    hit = PBSCandidate(
        motif_id="motif1_tRNA=Met",
        trna="Met",
        start=120,
        end=140,
        strand="+",
        identity_percent=96.5,
        alignment_length=18,
        motif_length=18,
        motif_coverage=1.0,
        distance_from_ltr5_nt=2,
        bitscore=30.0,
        evalue=1e-5,
    )
    feat = pbs_feature_from_candidate(hit, motif_database=default_motif_fasta())
    assert feat.name == "PBS_Met"
    assert feat.feature_type == "PBS"
    assert feat.pident == 96.5
    assert feat.coverage == 1.0
    assert feat.source == "PBS motif: Neumann2019"
    assert feat.subject_id == "motif1_tRNA=Met"
    assert feat.color == FEATURE_TYPE_COLORS["PBS"]
    assert FEATURE_TYPE_COLORS["PBS"] != FEATURE_TYPE_COLORS["domain"]

    tip = feature_hover_detail_text(feat, primary=feat.name)
    assert "Source: PBS motif: Neumann2019" in tip
    assert "Identity: 96.5%" in tip
    assert "Coverage: 100.0%" in tip
    assert feature_has_hover_meta(feat) is True

    stamped = _stamp_track([feat], track=3)[0]
    assert stamped.track == 3
    assert stamped.pident == 96.5
    assert stamped.coverage == 1.0
    assert stamped.source == "PBS motif: Neumann2019"
    assert stamped.subject_id == "motif1_tRNA=Met"


def test_extract_ltrs_prefers_primary_pair() -> None:
    features = [
        SequenceFeature(1, 100, "LTR5", "LTR"),
        SequenceFeature(500, 600, "LTR3", "LTR"),
        SequenceFeature(50, 80, "LTR5_2", "LTR"),
        SequenceFeature(700, 800, "LTR3_2", "LTR"),
    ]
    left, right = _extract_ltrs(features)
    assert left.name == "LTR5"
    assert right.name == "LTR3"


def test_extract_ltrs_outermost_fallback() -> None:
    features = [
        SequenceFeature(10, 50, "left", "LTR"),
        SequenceFeature(900, 950, "right", "LTR"),
        SequenceFeature(100, 200, "mid", "ORF"),
    ]
    left, right = _extract_ltrs(features)
    assert (left.start, right.start) == (10, 900)


def test_extract_ltrs_requires_two() -> None:
    assert _extract_ltrs([SequenceFeature(1, 10, "LTR5", "LTR")]) is None


def test_pbs_mode_value() -> None:
    assert AutoAnnotationMode.PBS.value == "PBS"


def test_pbs_defaults_match_refined_filters() -> None:
    from repeatviewer.backends.auto_annotation_settings import (
        DEFAULT_PBS_MAX_DISTANCE,
        DEFAULT_PBS_MIN_COVERAGE,
        DEFAULT_PBS_MIN_IDENTITY,
        DEFAULT_PBS_WINDOW,
        PbsSettings,
    )
    from repeatviewer.backends.pbs_detect import PBSSettings

    runtime = PBSSettings()
    ui = PbsSettings()
    assert runtime.window == DEFAULT_PBS_WINDOW == 40
    assert runtime.min_identity == DEFAULT_PBS_MIN_IDENTITY == 95.0
    assert runtime.min_coverage == DEFAULT_PBS_MIN_COVERAGE == 0.90
    assert runtime.max_distance_from_ltr5_nt == DEFAULT_PBS_MAX_DISTANCE == 20
    assert ui.max_distance_from_ltr5_nt == 20


def test_pbs_ranking_prefers_identity_and_length() -> None:
    from repeatviewer.backends.pbs_detect import _ranking

    short_full = PBSCandidate(
        motif_id="short",
        trna="Cys",
        start=101,
        end=112,
        strand="+",
        identity_percent=100.0,
        alignment_length=12,
        motif_length=12,
        motif_coverage=1.0,
        distance_from_ltr5_nt=0,
        bitscore=20.0,
        evalue=1e-3,
    )
    longer = PBSCandidate(
        motif_id="long",
        trna="Met",
        start=101,
        end=118,
        strand="+",
        identity_percent=100.0,
        alignment_length=18,
        motif_length=18,
        motif_coverage=1.0,
        distance_from_ltr5_nt=0,
        bitscore=30.0,
        evalue=1e-5,
    )
    assert _ranking(longer) < _ranking(short_full)


def test_candidate_rows_reject_gapped_and_distant(tmp_path: Path) -> None:
    from repeatviewer.backends.pbs_detect import PBSSettings, _candidate_rows

    # outfmt 6: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    table = tmp_path / "hits.tsv"
    table.write_text(
        "\n".join(
            [
                # ungapped, close — keep
                "q\tmotif_tRNA=Met\t100.0\t18\t0\t0\t1\t18\t1\t18\t1e-5\t36.2",
                # gapped — reject
                "q\tmotif_tRNA=Gap\t100.0\t18\t0\t1\t1\t18\t1\t18\t1e-5\t30.0",
                # too far from LTR5 (distance 25 with start at 126, anchor 100)
                "q\tmotif_tRNA=Far\t100.0\t18\t0\t0\t26\t43\t1\t18\t1e-5\t36.2",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    settings = PBSSettings(
        window=40,
        min_length=15,
        min_identity=95.0,
        min_coverage=0.90,
        max_distance_from_ltr5_nt=20,
    )
    hits = _candidate_rows(
        table,
        {"motif_tRNA=Met": 18, "motif_tRNA=Gap": 18, "motif_tRNA=Far": 18},
        "+",
        101,
        140,
        100,
        settings,
    )
    assert len(hits) == 1
    assert hits[0].trna == "Met"
    assert hits[0].distance_from_ltr5_nt == 0
