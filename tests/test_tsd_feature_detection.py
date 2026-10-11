"""Per-sequence TSD feature detection (MSA annotation, no clustering)."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

from repeatviewer.backends.auto_annotation import AutoAnnotationMode
from repeatviewer.backends.sequence_features import load_features_by_id_from_json
from repeatviewer.backends.te_annotation_batch import (
    TeAnnotationOptions,
    run_te_annotation_on_fasta,
)
from repeatviewer.backends.tsd_search import detect_tsd_features_by_id, scan_tsds_per_sequence
from repeatviewer.job_options import TsdSearchOptions


def _write_msa_with_tsd(path: Path) -> None:
    # Distinct 5-mer TSD "ATGCA" flanking a short TE body on both rows.
    # Layout: ATGCA + TE(20) + ATGCA  (length 30)
    te = "CGTACGTACGTACGTACGTA"
    row = f"ATGCA{te}ATGCA"
    path.write_text(f">a\n{row}\n>b\n{row}\n", encoding="utf-8")


def test_scan_tsds_per_sequence_finds_flanks(tmp_path: Path) -> None:
    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)
    hits = scan_tsds_per_sequence(
        fasta,
        left_boundary=5,
        right_boundary=24,
        nt_radius=10,
        record_ids=["a", "b"],
    )
    assert hits
    assert {h.seq_id for h in hits} == {"a", "b"}
    assert all(h.tsd_seq.upper() == "ATGCA" for h in hits)


def test_detect_tsd_features_by_id_annotates_both_flanks(tmp_path: Path) -> None:
    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)
    by_id = detect_tsd_features_by_id(
        fasta,
        record_ids=["a", "b"],
        left_boundary=5,
        right_boundary=24,
        nt_radius=10,
    )
    assert set(by_id) == {"a", "b"}
    for feats in by_id.values():
        assert len(feats) >= 2
        assert all(f.feature_type == "TSD" for f in feats)
        assert all(str(f.name).startswith("TSD_") for f in feats)
        starts = sorted(f.start for f in feats)
        # Left flank near start, right flank near end (1-based).
        assert starts[0] <= 5
        assert starts[-1] >= 26


def test_run_te_annotation_tsd_only_writes_sidecar(tmp_path: Path, monkeypatch) -> None:
    fasta = tmp_path / "msa.fa"
    _write_msa_with_tsd(fasta)

    monkeypatch.setattr(
        "repeatviewer.backends.te_annotation_batch.load_tsd_search_options",
        lambda: TsdSearchOptions(left_boundary=5, right_boundary=24, nt_radius=10),
    )

    # Ensure we do not fall into per-sequence ORF/LTR path.
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
        lambda: TsdSearchOptions(left_boundary=5, right_boundary=24, nt_radius=10),
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
    assert result.output_paths == []
    assert resolve_sidecar_path(fasta) is None
    assert result.features_by_id["a"]
    assert result.features_by_id["b"]
    assert all(f.feature_type == "TSD" for f in result.features_by_id["a"])
