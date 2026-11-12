"""TSD feature detection (MSA annotation with TEtrimmer-style filter/cluster)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from repeatviewer.backends.auto_annotation import AutoAnnotationMode
from repeatviewer.backends.sequence_features import load_features_by_id_from_json
from repeatviewer.backends.te_annotation_batch import (
    TeAnnotationOptions,
    run_te_annotation_on_fasta,
)
from repeatviewer.backends.tsd_search import (
    detect_tsd_features_by_id,
    process_tsd_hits,
    resolve_tsd_boundaries,
    scan_tsds_per_sequence,
)
from repeatviewer.job_options import TsdSearchOptions

# Divergent insertion-site TSDs (TEtrimmer requires cross-copy diversity).
# Layout per row: TSD(5) + TE(20) + TSD(5)  → length 30
# 1-based TE edges: left=6, right=25
_TE_BODY = "CGTACGTACGTACGTACGTA"
_ROW_A = f"ATGCA{_TE_BODY}ATGCA"
_ROW_B = f"GCTAA{_TE_BODY}GCTAA"
_TE_LEFT_1BASED = 6
_TE_RIGHT_1BASED = 25


def _write_msa_with_tsd(path: Path) -> None:
    path.write_text(f">a\n{_ROW_A}\n>b\n{_ROW_B}\n", encoding="utf-8")


def test_te_edges_match_boundaries_within_tolerance() -> None:
    from repeatviewer.backends.tsd_search import _te_edges_match_boundaries

    # Exact: te_start 0-based 5 → 1-based 6; te_end 0-based 24 → 1-based 25.
    assert _te_edges_match_boundaries(5, 24, 6, 25, tolerance=5)
    assert _te_edges_match_boundaries(5, 24, 6, 25, tolerance=0)
    # Within 5 bp.
    assert _te_edges_match_boundaries(5, 24, 8, 22, tolerance=5)
    # Outside tolerance.
    assert not _te_edges_match_boundaries(5, 24, 15, 25, tolerance=5)
    # Unset boundaries.
    assert not _te_edges_match_boundaries(5, 24, 0, 0, tolerance=5)


def test_resolve_tsd_boundaries_zero_means_full_msa() -> None:
    left, right = resolve_tsd_boundaries(30, 0, 0, nt_radius=10)
    assert left == right == 14


def test_resolve_tsd_boundaries_converts_1based_ui() -> None:
    left, right = resolve_tsd_boundaries(
        30, _TE_LEFT_1BASED, _TE_RIGHT_1BASED, nt_radius=10
    )
    assert left == 5
    assert right == 24


def test_resolve_tsd_search_scope_explicit_when_boundaries_set(tmp_path: Path) -> None:
    from repeatviewer.backends.tsd_search import resolve_tsd_search_scope

    fasta = tmp_path / "demo.raw.fa"
    fasta.write_text(">a\n" + ("A" * 100) + "\n>b\n" + ("A" * 100) + "\n", encoding="utf-8")
    left, right, radius, mode = resolve_tsd_search_scope(fasta, 100, 20, 80, 50)
    assert mode == "explicit"
    assert left == 19
    assert right == 79
    assert radius == 50


def test_resolve_tsd_search_scope_whole_msa_when_unset(tmp_path: Path) -> None:
    from repeatviewer.backends.tsd_search import resolve_tsd_search_scope

    fasta = tmp_path / "demo.raw.fa"
    fasta.write_text(">a\n" + ("A" * 100) + "\n>b\n" + ("A" * 100) + "\n", encoding="utf-8")
    left, right, radius, mode = resolve_tsd_search_scope(fasta, 100, 0, 0, 50)
    assert mode == "whole_msa"
    assert left == right == 49
    assert radius == 50


def test_msa_midpoint_halves_split() -> None:
    from repeatviewer.backends.tsd_search import _msa_midpoint_halves

    # 10 columns → mid=4; left [0,4), right [4,10)
    left, right = _msa_midpoint_halves("ACGTACGTAA")
    assert left == [0, 1, 2, 3]
    assert right == [4, 5, 6, 7, 8, 9]
    # Gaps are skipped inside each half.
    left, right = _msa_midpoint_halves("AC-TACGTAA")
    assert left == [0, 1, 3]
    assert right == [4, 5, 6, 7, 8, 9]


def test_msa_whole_search_windows_use_terminal_radius() -> None:
    from repeatviewer.backends.tsd_search import _msa_whole_search_windows

    left, right = _msa_whole_search_windows("ACGTACGTAA", nt_radius=2)
    assert left == [0, 1]
    assert right == [8, 9]


def test_whole_msa_scan_stays_fast_on_long_alignment() -> None:
    import time

    from repeatviewer.backends.tsd_search import _scan_record_tsd_hits

    ltr = ("TGTTGG" + "AT" * 40 + "CGTA" * 20) * 12
    body = ("ACGT" * 800 + "AAAA" * 80 + "TATA" * 80) * 3
    seq = ltr + body + ltr
    assert len(seq) > 15_000
    t0 = time.perf_counter()
    hits = _scan_record_tsd_hits(
        seq, seq_id="x", nt_radius=100, k_min=4, k_max=20, whole_msa=True
    )
    assert time.perf_counter() - t0 < 2.0
    assert isinstance(hits, list)


def test_preferred_boundaries_uses_existing_ltr_prediction(tmp_path: Path) -> None:
    """Auto detect must reuse Boundary prediction LTR without re-running."""
    from repeatviewer.backends.boundary_predict import BoundaryPrediction
    from repeatviewer.backends.tsd_search import (
        auto_tsd_span_source,
        preferred_tsd_boundaries_1based,
    )
    from repeatviewer.job_options import TSD_BOUNDARY_REPEATVIEWER

    fasta = tmp_path / "msa.fa"
    fasta.write_text(">a\nACGT\n>b\nACGT\n", encoding="utf-8")
    # 0-based half-open [4740, 12150) → 1-based inclusive 4741..12150
    pred = BoundaryPrediction(
        terminal_repeat=(4740, 12150),
        terminal_repeat_type="LTR",
        blast_coverage=(4719, 12150),
        msa_identity=(4742, 12150),
    )
    assert auto_tsd_span_source(pred) == "LTR"
    left, right = preferred_tsd_boundaries_1based(
        fasta,
        18000,
        boundary_source=TSD_BOUNDARY_REPEATVIEWER,
        prediction=pred,
    )
    assert (left, right) == (4741, 12150)


def test_scan_tsds_per_sequence_finds_flanks(tmp_path: Path) -> None:
    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)
    hits = scan_tsds_per_sequence(
        fasta,
        left_boundary=_TE_LEFT_1BASED,
        right_boundary=_TE_RIGHT_1BASED,
        nt_radius=10,
        record_ids=["a", "b"],
    )
    assert hits
    assert {h.seq_id for h in hits} == {"a", "b"}
    by_id = {h.seq_id: h.tsd_seq.upper() for h in hits}
    assert by_id["a"] == "ATGCA"
    assert by_id["b"] == "GCTAA"


def test_detect_tsd_features_by_id_annotates_both_flanks(tmp_path: Path) -> None:
    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)
    by_id = detect_tsd_features_by_id(
        fasta,
        record_ids=["a", "b"],
        left_boundary=_TE_LEFT_1BASED,
        right_boundary=_TE_RIGHT_1BASED,
        nt_radius=10,
    )
    assert set(by_id) == {"a", "b"}
    assert {f.name for f in by_id["a"]} == {"TSD_ATGCA"}
    assert {f.name for f in by_id["b"]} == {"TSD_GCTAA"}
    for feats in by_id.values():
        assert len(feats) >= 2
        assert all(f.feature_type == "TSD" for f in feats)
        starts = sorted(f.start for f in feats)
        assert starts[0] <= 5
        assert starts[-1] >= 26


def test_detect_keeps_boundary_flush_besides_top_cluster(tmp_path: Path) -> None:
    """TEtrimmer path: top cluster plus hits abutting defined boundaries.

    Layout (length 40): noise TSD pair mid-MSA (large divergent cluster) and
    true boundary-flush TSD at columns 1-5 / 36-40 (TE edges 6..35).
    """
    te = "CGTACGTACGTACGTACGTACGTACGTACGTA"  # 32 nt
    # Mid noise flanks (divergent) sit inside the TE body window for clustering.
    # Boundary flanks: ATGCA / GCTAA / TTGGC / CAAGT at the true edges.
    rows = [
        ("a", f"ATGCA{te}ATGCA"),
        ("b", f"GCTAA{te}GCTAA"),
        ("c", f"TTGGC{te}TTGGC"),
        ("d", f"CAAGT{te}CAAGT"),
    ]
    # Pad so mid identical-ish conserved motifs can form a competing cluster
    # further from the boundaries (still divergent enough to pass diversity).
    assert all(len(seq) == 42 for _, seq in rows)
    fasta = tmp_path / "msa.fa"
    fasta.write_text(
        "".join(f">{rid}\n{seq}\n" for rid, seq in rows),
        encoding="utf-8",
    )
    by_id = detect_tsd_features_by_id(
        fasta,
        record_ids=[rid for rid, _ in rows],
        left_boundary=6,
        right_boundary=37,
        nt_radius=20,
        k_min=4,
        k_max=20,
        min_distance=0,
        max_sim_threshold=0.3,
    )
    assert set(by_id) == {"a", "b", "c", "d"}
    # Boundary-flush motifs at the defined edges must be annotated.
    assert {f.name for f in by_id["a"]} >= {"TSD_ATGCA"}
    assert {f.name for f in by_id["b"]} >= {"TSD_GCTAA"}
    for rid, feats in by_id.items():
        lefts = [f for f in feats if f.start <= 5]
        rights = [f for f in feats if f.start >= 36]
        assert lefts, rid
        assert rights, rid


def test_detect_tsd_without_boundaries_searches_full_msa(tmp_path: Path) -> None:
    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)
    by_id = detect_tsd_features_by_id(
        fasta,
        record_ids=["a", "b"],
        left_boundary=0,
        right_boundary=0,
        nt_radius=10,
    )
    assert set(by_id) == {"a", "b"}
    assert all(len(feats) >= 2 for feats in by_id.values())


def test_run_te_annotation_tsd_only_writes_sidecar(tmp_path: Path, monkeypatch) -> None:
    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)

    monkeypatch.setattr(
        "repeatviewer.backends.te_annotation_batch.load_tsd_search_options",
        lambda: TsdSearchOptions(
            left_boundary=_TE_LEFT_1BASED,
            right_boundary=_TE_RIGHT_1BASED,
            nt_radius=10,
            min_distance=0,
            write_text_report=True,
        ),
    )

    with patch(
        "repeatviewer.backends.te_annotation_batch._annotate_one_record_fasta"
    ) as per_seq:
        result = run_te_annotation_on_fasta(
            fasta,
            TeAnnotationOptions(auto_modes=(AutoAnnotationMode.TSD,)),
        )
        per_seq.assert_not_called()

    assert result.output_paths
    assert result.json_path is not None
    assert result.tsd_report_path is not None
    assert result.tsd_report_path.is_file()
    report = result.tsd_report_path.read_text(encoding="utf-8")
    assert "### CLUSTER 1 ###" in report
    assert "TSD:" in report
    features = load_features_by_id_from_json(fasta, ["a", "b"])
    assert features["a"]
    assert features["b"]
    assert all(f.feature_type == "TSD" for f in features["a"])
    assert all(f.feature_type == "TSD" for f in features["b"])


def test_run_te_annotation_persist_false_skips_sidecar(tmp_path: Path, monkeypatch) -> None:
    """Preview UI path: detections stay in memory until the user clicks Save."""
    from repeatviewer.backends.sequence_features import resolve_sidecar_path

    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)

    monkeypatch.setattr(
        "repeatviewer.backends.te_annotation_batch.load_tsd_search_options",
        lambda: TsdSearchOptions(
            left_boundary=_TE_LEFT_1BASED,
            right_boundary=_TE_RIGHT_1BASED,
            nt_radius=10,
            min_distance=0,
            write_text_report=True,
        ),
    )

    with patch(
        "repeatviewer.backends.te_annotation_batch._annotate_one_record_fasta"
    ) as per_seq:
        result = run_te_annotation_on_fasta(
            fasta,
            TeAnnotationOptions(auto_modes=(AutoAnnotationMode.TSD,)),
            persist=False,
        )
        per_seq.assert_not_called()

    assert result.json_path is None
    assert resolve_sidecar_path(fasta) is None
    # Clustered TSD text report is written even on Preview (persist=False).
    assert result.tsd_report_path is not None
    assert result.tsd_report_path.is_file()
    assert result.output_paths == [result.tsd_report_path.resolve()]
    assert result.features_by_id["a"]
    assert result.features_by_id["b"]
    assert all(f.feature_type == "TSD" for f in result.features_by_id["a"])


def test_scan_respects_min_tsd_length(tmp_path: Path) -> None:
    """Whole-MSA scan: no hit shorter than k_min (and none if only shorter exist)."""
    fasta = tmp_path / "msa.fa"
    # Halves share no ≥6-mer; 5-bp flanks alone stay below k_min=6.
    fasta.write_text(
        ">a\nAAACGTAAAAACGTTGCATTTTTGCAAA\n"
        ">b\nGGCGTACCCCCGTAATGCGGGGGTGCGG\n",
        encoding="utf-8",
    )
    hits = scan_tsds_per_sequence(
        fasta,
        left_boundary=0,
        right_boundary=0,
        k_min=6,
        k_max=20,
        record_ids=["a", "b"],
    )
    assert hits == []
    assert all(len(h.tsd_seq) >= 6 for h in hits)


def test_process_tsd_hits_rejects_identical_motifs() -> None:
    """Same TSD across all copies fails TEtrimmer diversity (avg sim = 1)."""
    hits = {
        (5, 25): [
            ("a", "TAAT"),
            ("b", "TAAT"),
            ("c", "TAAT"),
            ("d", "TAAT"),
        ]
    }
    assert process_tsd_hits(hits, max_sim_threshold=0.3, min_length=4) == {}


def test_process_tsd_hits_keeps_divergent_motifs() -> None:
    hits = {
        (5, 25): [
            ("a", "ATGCA"),
            ("b", "GCTAA"),
            ("c", "TTGGC"),
            ("d", "CAAGT"),
        ]
    }
    out = process_tsd_hits(hits, max_sim_threshold=0.3, min_length=4)
    motifs = {seq.upper() for matches in out.values() for _, seq in matches}
    assert motifs == {"ATGCA", "GCTAA", "TTGGC", "CAAGT"}


def test_process_tsd_hits_does_not_shrink_below_min_length() -> None:
    """After TEtrimmer poly-N trim (n≥4), drop motifs shorter than Min TSD length."""
    # Divergent 5-mers sharing only the last letter → trim one end → still ≥4.
    hits = {
        (5, 25): [
            ("a", "ATGCT"),
            ("b", "CGGGT"),
            ("c", "TTCCT"),
            ("d", "GACAT"),
        ]
    }
    out = process_tsd_hits(hits, min_length=4, nt_threshold=0.9)
    motifs = [seq for matches in out.values() for _, seq in matches]
    assert motifs
    assert all(len(seq) >= 4 for seq in motifs)


def test_detect_rejects_identical_conserved_flanks(tmp_path: Path) -> None:
    """All-identical TAAT across rows is not annotated (not divergent TSDs)."""
    fasta = tmp_path / "msa.fa"
    te = "CGTACGTACGTACGTACGTA"
    row = f"TAAT{te}TAAT"
    fasta.write_text(
        f">a\n{row}\n>b\n{row}\n>c\n{row}\n>d\n{row}\n",
        encoding="utf-8",
    )
    by_id = detect_tsd_features_by_id(
        fasta,
        record_ids=["a", "b", "c", "d"],
        left_boundary=5,
        right_boundary=24,
        nt_radius=10,
        k_min=4,
        k_max=20,
        max_sim_threshold=0.3,
    )
    assert by_id == {}


def test_detect_rejects_short_motifs(tmp_path: Path) -> None:
    """Annotated TSD motifs are never shorter than Min TSD length."""
    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)
    by_id = detect_tsd_features_by_id(
        fasta,
        record_ids=["a", "b"],
        left_boundary=0,
        right_boundary=0,
        k_min=4,
        k_max=20,
        min_distance=0,
    )
    assert by_id
    for feats in by_id.values():
        for feat in feats:
            motif = str(feat.name).removeprefix("TSD_")
            assert len(motif) >= 4


def test_scan_respects_min_distance(tmp_path: Path) -> None:
    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)
    hits = scan_tsds_per_sequence(
        fasta,
        left_boundary=_TE_LEFT_1BASED,
        right_boundary=_TE_RIGHT_1BASED,
        nt_radius=10,
        min_distance=50,
        record_ids=["a", "b"],
    )
    assert hits == []


def test_preferred_boundaries_from_sidecar(tmp_path: Path) -> None:
    import json

    from repeatviewer.backends.tsd_search import preferred_tsd_boundaries_1based
    from repeatviewer.job_options import TSD_BOUNDARY_REPEATVIEWER

    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)
    # 0-based half-open final span [5, 25) → 1-based inclusive 6..25
    sidecar = tmp_path / "msa.boundaries_RV.json"
    sidecar.write_text(
        json.dumps(
            {
                "msa_width": 30,
                "final": {"start": 5, "end": 25},
            }
        ),
        encoding="utf-8",
    )
    left, right = preferred_tsd_boundaries_1based(
        fasta,
        30,
        boundary_source=TSD_BOUNDARY_REPEATVIEWER,
        saved_left=1,
        saved_right=2,
    )
    assert (left, right) == (6, 25)


def test_tool_tsd_search_skips_text_report_by_default(tmp_path: Path) -> None:
    """Tools → TSD search: annotation JSON only; no .tsd_search.txt by default."""
    from repeatviewer.backends.curation import run_on_sequence
    from repeatviewer.backends.sequence_features import load_features_by_id_from_json
    from repeatviewer.project import RepeatViewerProject

    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)
    outs = run_on_sequence(
        "tsd_search",
        fasta,
        RepeatViewerProject(),
        tsd_search_options=TsdSearchOptions(
            left_boundary=_TE_LEFT_1BASED,
            right_boundary=_TE_RIGHT_1BASED,
            nt_radius=10,
            min_distance=0,
        ),
    )
    text_files = [p for p in outs if str(p).endswith(".tsd_search.txt")]
    assert not text_files
    assert not list(tmp_path.glob("*.tsd_search.txt"))
    features = load_features_by_id_from_json(fasta, ["a", "b"])
    assert features["a"]
    assert all(f.feature_type == "TSD" for f in features["a"])


def test_tool_tsd_search_writes_text_when_enabled(tmp_path: Path) -> None:
    """Tools → TSD search: clustered text report when write_text_report is on."""
    from repeatviewer.backends.curation import run_on_sequence
    from repeatviewer.backends.sequence_features import load_features_by_id_from_json
    from repeatviewer.project import RepeatViewerProject

    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)
    outs = run_on_sequence(
        "tsd_search",
        fasta,
        RepeatViewerProject(),
        tsd_search_options=TsdSearchOptions(
            left_boundary=_TE_LEFT_1BASED,
            right_boundary=_TE_RIGHT_1BASED,
            nt_radius=10,
            min_distance=0,
            write_text_report=True,
        ),
    )
    text_files = [p for p in outs if str(p).endswith(".tsd_search.txt")]
    assert text_files
    assert Path(text_files[0]).is_file()
    features = load_features_by_id_from_json(fasta, ["a", "b"])
    assert features["a"]
    assert all(f.feature_type == "TSD" for f in features["a"])
