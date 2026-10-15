"""Visibility keys and filters for Display → Tracks annotation kinds."""

from __future__ import annotations

from PySide6.QtWidgets import QApplication

from repeatviewer.backends.sequence_features import (
    SequenceFeature,
    annotation_visibility_key,
    annotation_visibility_label,
    collect_annotation_visibility_keys,
    domain_kind_from_feature_name,
    domain_visibility_group,
    filter_features_by_annotation_visibility,
)


def test_domain_kind_strips_evalue_and_percent() -> None:
    assert domain_kind_from_feature_name("Pfam: RVT_1 (E=1e-20)") == "Pfam: RVT_1"
    assert (
        domain_kind_from_feature_name("Pfam: Integrase (E=0.001; 92%)")
        == "Pfam: Integrase"
    )
    assert domain_kind_from_feature_name("RVT_1") == "RVT_1"


def test_domain_visibility_groups_by_database() -> None:
    assert domain_visibility_group("Pfam: RVT_1 (E=1e-20)") == "Pfam"
    assert domain_visibility_group("Pfam: Integrase") == "Pfam"
    assert domain_visibility_group("CDD: rve") == "CDD"
    assert domain_visibility_group("RVT_1") == "domain"


def test_annotation_visibility_keys_group_pfam() -> None:
    feats = [
        SequenceFeature(1, 10, "left", feature_type="LTR"),
        SequenceFeature(20, 40, "orf1", feature_type="ORF"),
        SequenceFeature(50, 80, "Pfam: RVT_1 (E=1e-20)", feature_type="domain"),
        SequenceFeature(90, 110, "Pfam: Integrase (E=1e-5)", feature_type="domain"),
        SequenceFeature(120, 140, "Pfam: Chromo (E=1e-8)", feature_type="domain"),
    ]
    assert annotation_visibility_key(feats[0]) == "LTR"
    assert annotation_visibility_key(feats[2]) == "domain::Pfam"
    assert annotation_visibility_label("domain::Pfam") == "Pfam"
    keys = collect_annotation_visibility_keys(feats)
    assert keys[:4] == ["LTR", "TIR", "ORF", "TSD"]
    assert keys.count("domain::Pfam") == 1
    assert "domain::Pfam: RVT_1" not in keys


def test_filter_features_by_annotation_visibility() -> None:
    feats = [
        SequenceFeature(1, 10, "left", feature_type="LTR"),
        SequenceFeature(20, 40, "orf1", feature_type="ORF"),
        SequenceFeature(50, 80, "Pfam: RVT_1 (E=1e-20)", feature_type="domain"),
        SequenceFeature(90, 100, "Pfam: rve", feature_type="domain"),
    ]
    assert filter_features_by_annotation_visibility(feats, None) == feats
    only_ltr = filter_features_by_annotation_visibility(feats, {"LTR"})
    assert [f.feature_type for f in only_ltr] == ["LTR"]
    domains = filter_features_by_annotation_visibility(feats, {"domain::Pfam"})
    assert len(domains) == 2
    assert all(f.feature_type == "domain" for f in domains)
    assert filter_features_by_annotation_visibility(feats, set()) == []


def test_rail_nested_kind_checks(qapp: QApplication) -> None:
    from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail

    rail = SequenceAnnotationRail()
    rail.show()
    qapp.processEvents()
    rail.set_overview_display_mode(True)
    assert rail.show_annotation_track.isHidden()
    assert not rail.displayed_annotations_label.isHidden()
    assert not rail._annotation_kinds_host.isHidden()
    assert not rail.annotation_kinds_all_check.isHidden()
    assert "LTR" in rail.annotation_kind_checks
    assert "TIR" in rail.annotation_kind_checks
    assert "ORF" in rail.annotation_kind_checks
    assert "TSD" in rail.annotation_kind_checks

    feats = [
        SequenceFeature(1, 10, "left", feature_type="LTR"),
        SequenceFeature(50, 80, "Pfam: RVT_1 (E=1e-20)", feature_type="domain"),
        SequenceFeature(55, 70, "Pfam: Chromo (E=1e-3)", feature_type="domain"),
        SequenceFeature(90, 100, "tandem", feature_type="tandem repeat"),
    ]
    rail._hidden_annotation_kinds.clear()
    rail._rebuild_annotation_kind_checks(
        [
            "LTR",
            "TIR",
            "ORF",
            "TSD",
            "tandem repeat",
            "domain::Pfam",
        ]
    )
    assert "domain::Pfam" in rail.annotation_kind_checks
    assert "domain::Pfam: RVT_1" not in rail.annotation_kind_checks
    assert "tandem repeat" in rail.annotation_kind_checks
    assert rail.annotation_kind_checks["domain::Pfam"].text() == "Pfam"
    assert rail.annotation_kinds_all_check.isChecked()

    rail.annotation_kind_checks["LTR"].setChecked(False)
    visible = rail.visible_annotation_kinds()
    assert "LTR" not in visible
    assert "domain::Pfam" in visible
    assert not rail.annotation_kinds_all_check.isChecked()
    assert rail.show_annotation_track.isChecked()

    rail.annotation_kinds_all_check.setChecked(True)
    assert rail.annotation_kinds_all_check.isChecked()
    rail.annotation_kinds_all_check.setChecked(False)
    assert rail.visible_annotation_kinds() == frozenset()
    assert not rail.show_annotation_track.isChecked()
    assert not rail._annotation_kinds_host.isHidden()
    # Restore defaults so later tests are not polluted via QSettings.
    rail.annotation_kinds_all_check.setChecked(True)
    rail.close()


def test_bars_filter_packed_tracks(qapp: QApplication) -> None:
    from repeatviewer.widgets.multi_sequence_overview import _BarsCanvas

    canvas = _BarsCanvas()
    feats = {
        "seq1": [
            SequenceFeature(1, 10, "left", feature_type="LTR"),
            SequenceFeature(20, 40, "orf1", feature_type="ORF"),
            SequenceFeature(50, 80, "Pfam: RVT_1 (E=1e-20)", feature_type="domain"),
            SequenceFeature(90, 100, "Pfam: rve", feature_type="domain"),
        ]
    }
    canvas.set_data(["seq1"], {"seq1": 100}, feats)
    assert sum(len(t) for t in canvas._tracks["seq1"]) == 4

    canvas.set_visible_annotation_kinds(frozenset({"ORF", "domain::Pfam"}))
    packed = [f for track in canvas._tracks["seq1"] for f in track]
    assert {f.feature_type for f in packed} == {"ORF", "domain"}
    assert len(packed) == 3
    canvas.close()
