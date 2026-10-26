"""Visibility keys and filters for Display → Tracks annotation kinds."""

from __future__ import annotations

from PySide6.QtWidgets import QApplication, QToolButton

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
    # Only types present on the features (not empty TIR/TSD/PBS placeholders).
    assert keys == ["LTR", "ORF", "domain::Pfam"]
    assert keys.count("domain::Pfam") == 1
    assert "domain::Pfam: RVT_1" not in keys
    assert "TIR" not in keys
    assert "TSD" not in keys
    assert "PBS" not in keys


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
    # Empty until features are provided — no placeholder LTR/TIR/…
    assert rail.displayed_annotations_label.isHidden()
    assert rail._annotation_kinds_host.isHidden()
    assert rail.annotation_kind_checks == {}
    # All sits on its own row above the wrapping type flow.
    assert rail._annotation_kinds_flow.indexOf(rail.annotation_kinds_all_check) < 0
    assert rail.annotation_kinds_all_check.parentWidget() is rail._annotation_kinds_host
    assert rail.show_standard_boundaries.text() == "Show reference boundaries"
    assert not rail.show_standard_boundaries.font().bold()
    assert rail.show_blast_coverage.text() == "Show BLAST coverage"
    assert rail.show_blast_coverage_row is not None
    help_btns = [
        w
        for w in rail.show_blast_coverage_row.findChildren(QToolButton)
        if w.text() == "?"
    ]
    assert len(help_btns) == 1
    assert "Show BLAST coverage" in help_btns[0].accessibleName()
    tip = help_btns[0].toolTip()
    assert "reference genome" in tip
    assert "Manage Database" in tip
    assert "project genome" not in tip
    assert "Define Genome" not in tip
    assert "Reading the plot" not in tip
    assert rail.sync_split_views.text() == "Sync left and right view"
    assert rail.sync_split_views.isHidden()
    rail.set_overview_display_mode(True)
    rail.set_msa_display_mode(True)
    rail.set_split_sync_available(True)
    assert not rail.sync_split_views.isHidden()
    rail.set_split_sync_available(False)
    assert rail.sync_split_views.isHidden()

    feats = [
        SequenceFeature(1, 10, "left", feature_type="LTR"),
        SequenceFeature(50, 80, "Pfam: RVT_1 (E=1e-20)", feature_type="domain"),
        SequenceFeature(55, 70, "Pfam: Chromo (E=1e-3)", feature_type="domain"),
        SequenceFeature(90, 100, "tandem", feature_type="tandem repeat"),
    ]
    rail._hidden_annotation_kinds.clear()
    rail.set_available_annotation_kinds(feats)
    assert not rail.displayed_annotations_label.isHidden()
    assert not rail._annotation_kinds_host.isHidden()
    assert list(rail.annotation_kind_checks) == [
        "LTR",
        "tandem repeat",
        "domain::Pfam",
    ]
    assert "TIR" not in rail.annotation_kind_checks
    assert "ORF" not in rail.annotation_kind_checks
    assert "TSD" not in rail.annotation_kind_checks
    assert "PBS" not in rail.annotation_kind_checks
    assert "domain::Pfam" in rail.annotation_kind_checks
    assert "domain::Pfam: RVT_1" not in rail.annotation_kind_checks
    assert rail.annotation_kind_checks["domain::Pfam"].text() == "Pfam"
    assert rail.annotation_kinds_all_check.isChecked()

    # Displayed annotations sits below Split annotation by source in the Tracks
    # section widget order.
    tracks_layout = rail._tracks_section._body_layout
    split_idx = -1
    label_idx = -1
    for i in range(tracks_layout.count()):
        item = tracks_layout.itemAt(i)
        w = item.widget() if item is not None else None
        if w is rail.split_tracks_by_source:
            split_idx = i
        if w is rail.displayed_annotations_label:
            label_idx = i
    assert split_idx >= 0 and label_idx > split_idx

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


def test_ensure_annotation_kinds_visible_after_run(qapp: QApplication) -> None:
    """Feature/similarity runs re-tick Displayed annotations that were off."""
    from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail

    rail = SequenceAnnotationRail()
    rail.show()
    qapp.processEvents()
    rail.set_overview_display_mode(True)
    rail._hidden_annotation_kinds.clear()
    rail._rebuild_annotation_kind_checks(
        ["LTR", "TIR", "ORF", "TSD", "domain::Pfam"]
    )
    rail.annotation_kind_checks["LTR"].setChecked(False)
    rail.annotation_kind_checks["domain::Pfam"].setChecked(False)
    assert "LTR" not in rail.visible_annotation_kinds()
    assert "domain::Pfam" not in rail.visible_annotation_kinds()

    rail.ensure_annotation_kinds_visible(["LTR", "domain::Pfam"])
    assert "LTR" in rail.visible_annotation_kinds()
    assert "domain::Pfam" in rail.visible_annotation_kinds()
    assert rail.annotation_kind_checks["LTR"].isChecked()
    assert rail.annotation_kind_checks["domain::Pfam"].isChecked()
    assert "LTR" not in rail._hidden_annotation_kinds
    assert "domain::Pfam" not in rail._hidden_annotation_kinds

    # New domain group appears and is ticked when first predicted.
    rail.ensure_annotation_kinds_visible(["domain::CDD"])
    assert "domain::CDD" in rail.annotation_kind_checks
    assert rail.annotation_kind_checks["domain::CDD"].isChecked()

    rail.annotation_kinds_all_check.setChecked(True)
    rail.close()


def test_rail_apply_visible_annotation_kinds_sync(qapp: QApplication) -> None:
    from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail

    rail = SequenceAnnotationRail()
    rail.show()
    rail._stack.setCurrentIndex(0)
    rail._tracks_section.set_expanded(True)
    feats = [
        SequenceFeature(1, 10, "a", feature_type="LTR"),
        SequenceFeature(1, 10, "b", feature_type="ORF"),
        SequenceFeature(1, 10, "c", feature_type="PBS"),
    ]
    rail.set_available_annotation_kinds(feats)
    assert rail.visible_annotation_kinds() == frozenset({"LTR", "ORF", "PBS"})
    rail.apply_visible_annotation_kinds(frozenset({"LTR", "PBS"}))
    assert rail.visible_annotation_kinds() == frozenset({"LTR", "PBS"})
    assert rail.annotation_kind_checks["ORF"].isChecked() is False
    rail.close()


def test_annotation_kinds_more_overflow(qapp: QApplication) -> None:
    """Tall type lists stay capped; extras are reachable from More."""
    from repeatviewer.widgets.sequence_annotation_rail import (
        SequenceAnnotationRail,
        _ANNOTATION_KINDS_FLOW_MAX_HEIGHT,
    )

    rail = SequenceAnnotationRail()
    rail.show()
    rail._stack.setCurrentIndex(0)
    rail._tracks_section.set_expanded(True)
    qapp.processEvents()

    feats = [
        SequenceFeature(1, 10, "a", feature_type="LTR"),
        SequenceFeature(1, 10, "b", feature_type="TIR"),
        SequenceFeature(1, 10, "c", feature_type="ORF"),
        SequenceFeature(1, 10, "d", feature_type="TSD"),
        SequenceFeature(1, 10, "e", feature_type="PBS"),
        SequenceFeature(1, 10, "f", feature_type="extracted region"),
        SequenceFeature(1, 10, "g", feature_type="primer_binding_site"),
        SequenceFeature(1, 10, "h", feature_type="protein_domain"),
    ]
    rail.set_available_annotation_kinds(feats)
    rail._annotation_kinds_flow_host.setFixedWidth(200)
    rail._sync_annotation_kinds_overflow()
    qapp.processEvents()

    assert rail._annotation_kinds_flow_host.maximumHeight() == (
        _ANNOTATION_KINDS_FLOW_MAX_HEIGHT
    )
    assert rail._annotation_kinds_overflow_keys
    assert not rail._annotation_kinds_more_btn.isHidden()
    assert "More" in rail._annotation_kinds_more_btn.text()
    for key in rail._annotation_kinds_overflow_keys:
        assert rail.annotation_kind_checks[key].isHidden()
    # Overflow types still participate in visibility filtering.
    assert "extracted region" in rail.visible_annotation_kinds()
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
