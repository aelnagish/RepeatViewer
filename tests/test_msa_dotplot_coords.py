"""MSA Dotplot (Self): ungapped plot coords map back to gapped Overview columns."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication

from repeatviewer.backends.msa_coordinates import (
    gapped_to_ungapped_1based,
    ungap_sequence,
    ungapped_to_gapped_1based,
)
from repeatviewer.widgets.dotplot_view import DotplotPanel
from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel


def test_ungap_round_trip_maps_dotplot_selection_to_msa_columns() -> None:
    gapped = "AC--GT"
    ungapped, keep = ungap_sequence(gapped)
    assert ungapped == "ACGT"
    assert keep is not None
    # Plot residue 3 (G) is gapped column 5.
    assert ungapped_to_gapped_1based(3, keep) == 5
    assert gapped_to_ungapped_1based(5, keep) == 3


def test_dotplot_panel_ungaps_msa_records_and_maps_selection(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">a\nAC--GT\n>b\nACGTGT\n")
    panel = DotplotPanel(path)
    panel.set_record_choices(
        [("a", "AC--GT"), ("b", "ACGTGT")],
        ungap=True,
        preferred_id="a",
    )
    assert panel._active_sequence == "ACGT"
    assert panel._coord_keep is not None
    # Numbered labels match Overview order; userData stays the raw id.
    assert panel._sequence_combo.itemText(0) == "1. a"
    assert panel._sequence_combo.itemText(1) == "2. b"
    assert panel._sequence_combo.itemData(0) == "a"
    assert panel._sequence_combo.currentData() == "a"
    # Plot coords 1-4 → MSA columns 1-6 (gaps included at ends of span).
    assert panel._plot_to_source_range(1, 4) == (1, 6)
    assert panel._source_to_plot_range(1, 6) == (1, 4)
    panel.close()


def test_overview_select_by_position_sets_region(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">a\nACGTAC\n>b\nACGTAC\n")
    overview = MultiSequenceOverviewPanel(path)
    overview.load_path(path, msa_mode=True)
    overview.select_by_position(2, 4, sequence_id="b")
    assert overview._bars.region_selection() == (2, 4)
    # Only the target sequence — not a full-height column across all rows.
    assert overview._bars.regional_rows() == {1}
    assert overview.selected_record_id() == "b"
    overview.close()


def test_overview_sort_updates_dotplot_sequence_order(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">zeta\nACGTAC\n>alpha\nACGTAC\n>beta\nACGTAC\n")
    overview = MultiSequenceOverviewPanel(path)
    overview.load_path(path, msa_mode=True)
    assert overview.record_ids_for_dotplot() == ["zeta", "alpha", "beta"]
    ordered: list[list[str]] = []
    overview.sequenceOrderChanged.connect(
        lambda: ordered.append(overview.record_ids_for_dotplot())
    )
    overview._on_sequence_sort_changed("name")
    assert overview.record_ids() == ["alpha", "beta", "zeta"]
    assert overview.record_ids_for_dotplot() == ["alpha", "beta", "zeta"]
    assert ordered == [["alpha", "beta", "zeta"]]
    overview.close()


def test_dotplot_set_record_choices_loads_sequence_lazily(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">a\nAC--GT\n>b\nACGTGT\n>c\nAAAAAA\n")
    loads: list[str] = []

    def loader(sid: str) -> str:
        loads.append(sid)
        return {"a": "AC--GT", "b": "ACGTGT", "c": "AAAAAA"}[sid]

    panel = DotplotPanel(path)
    panel.set_record_choices(
        ["a", "b", "c"],
        ungap=True,
        preferred_id="b",
        sequence_loader=loader,
    )
    assert panel._sequence_combo.count() == 3
    assert panel._sequence_combo.itemText(2) == "3. c"
    assert panel.sequence_id() == "b"
    assert panel._active_sequence == "ACGTGT"
    # Only the preferred sequence is fetched when the list is built.
    assert loads == ["b"]
    panel._sequence_combo.setCurrentIndex(0)
    qapp.processEvents()
    assert loads == ["b", "a"]
    assert panel._active_sequence == "ACGT"
    panel.close()


def test_dotplot_numbered_sequence_label(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">a\nACGTAC\n>b\nACGTAC\n")
    panel = DotplotPanel(path)
    panel.set_record_choices(
        [("a", "ACGTAC"), ("b", "ACGTAC")],
        preferred_id="b",
    )
    assert panel._numbered_sequence_label("b") == "2. b"
    assert panel._numbered_sequence_label("a") == "1. a"
    assert panel._numbered_sequence_label("") == ""
    panel.close()


def test_dotplot_sequence_prev_next_buttons(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "aln.fa"
    path.write_text(">a\nACGTAC\n>b\nACGTAC\n>c\nACGTAC\n")
    panel = DotplotPanel(path)
    panel.set_record_choices(
        [("a", "ACGTAC"), ("b", "ACGTAC"), ("c", "ACGTAC")],
        preferred_id="a",
    )
    assert panel._sequence_combo.currentData() == "a"
    assert panel._sequence_prev_btn.isEnabled()
    assert panel._sequence_next_btn.isEnabled()

    panel._rail.step_sequence(1)
    assert panel._sequence_combo.currentData() == "b"
    panel._rail.step_sequence(1)
    assert panel._sequence_combo.currentData() == "c"
    # Wraps last → first.
    panel._rail.step_sequence(1)
    assert panel._sequence_combo.currentData() == "a"
    # Previous wraps first → last.
    panel._rail.step_sequence(-1)
    assert panel._sequence_combo.currentData() == "c"
    panel.close()


def test_overview_region_for_dotplot_requires_matching_single_sequence(
    qapp: QApplication, tmp_path: Path
) -> None:
    from types import SimpleNamespace

    from repeatviewer.widgets.preview_tabs import FilePreviewPage

    path = tmp_path / "aln.fa"
    path.write_text(">a\nACGTAC\n>b\nACGTAC\n>c\nACGTAC\n")
    page = FilePreviewPage(path)
    overview = MultiSequenceOverviewPanel(path)
    overview.load_path(path, msa_mode=True)
    page._annotation_overview = overview
    page._multi_fasta_mode = True
    page._msa_aligned_mode = True
    panel = SimpleNamespace(sequence_id=lambda: "b")

    # Single regional row matching Dotplot sequence → highlight.
    overview.select_by_position(2, 4, sequence_id="b")
    assert page._overview_region_for_dotplot(panel) == (2, 4)

    # Same region on a different sequence → no highlight.
    panel = SimpleNamespace(sequence_id=lambda: "a")
    assert page._overview_region_for_dotplot(panel) is None

    # Multi-row regional selection → no highlight.
    overview._bars.set_regional_selection(2, 4, [0, 1])
    panel = SimpleNamespace(sequence_id=lambda: "a")
    assert page._overview_region_for_dotplot(panel) is None

    # Full column with no single matching row → no highlight.
    overview._bars.set_selected_rows([])
    overview._bars.set_region_selection(2, 4)
    assert page._overview_region_for_dotplot(panel) is None

    # Full column + exactly one matching selected row → highlight.
    # set_selected_rows clears region, so restore the column band after.
    overview._bars.set_selected_rows([0])
    overview._bars.set_region_selection(2, 4)
    assert page._overview_region_for_dotplot(panel) == (2, 4)

    overview.close()
    page.deleteLater()
    qapp.processEvents()
