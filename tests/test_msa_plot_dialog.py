"""MSA plot dialog: Align-style Inputs / Settings / Advanced colors."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QApplication,
    QGroupBox,
    QLabel,
    QPushButton,
    QRadioButton,
)

from repeatviewer.job_options import (
    ALIGN_INPUT_DEFINE,
    ALIGN_INPUT_SELECTED,
    DEFAULT_MSA_PLOT_DPI,
    MsaPlotJobOptions,
)
from repeatviewer.views.tool_dialogs import MsaPlotDialog
from repeatviewer.widgets.residue_color_grid import ResidueColorGrid


def _msa_plot_dialog(qapp: QApplication, **kw) -> MsaPlotDialog:
    del qapp
    kw.setdefault("initial", MsaPlotJobOptions())
    dlg = MsaPlotDialog(**kw)
    dlg.show()
    dlg._resize_for_msa_advanced(False)
    return dlg


def test_msa_plot_dialog_align_style_sections(qapp: QApplication) -> None:
    dlg = _msa_plot_dialog(qapp)
    labels = [lab.text() for lab in dlg.findChildren(QLabel)]
    assert "<b>Inputs</b>" in labels
    assert "<b>Settings</b>" in labels
    assert any(lab.text() == "Sequence type:" for lab in dlg.findChildren(QLabel))
    assert any(lab.text() == "Sequence height:" for lab in dlg.findChildren(QLabel))
    radios = [r.text() for r in dlg.findChildren(QRadioButton)]
    assert "Use selected sequences" in radios
    assert "Define inputs" in radios
    assert dlg.input_source() == ALIGN_INPUT_SELECTED
    texts = [b.text() for b in dlg.findChildren(QPushButton)]
    assert any(t.startswith("Advanced Settings") for t in texts)
    assert "Restore Defaults" in texts
    assert "▶ Color" not in texts
    assert dlg.width() == MsaPlotDialog._DIALOG_WIDTH
    assert dlg._DIALOG_WIDTH == 608
    assert dlg._seq_type.width() == 260
    assert dlg._msa_advanced_panel.isHidden()
    assert any(g.title() == "Nucleotide colors" for g in dlg.findChildren(QGroupBox))
    assert any(g.title() == "Amino-acid colors" for g in dlg.findChildren(QGroupBox))

    dlg._advanced_settings_btn.setChecked(True)
    qapp.processEvents()
    assert dlg._msa_advanced_panel.isVisible()
    assert dlg._advanced_settings_btn.text() == "Advanced Settings ▾"
    assert dlg._colors_scroll.maximumHeight() == MsaPlotDialog._COLORS_SCROLL_MAX_HEIGHT
    # Ambiguity is split across two rows (label on first only).
    from repeatviewer.backends.msa_plot import NUCLEOTIDE_CATEGORIES

    ambiguity_rows = [c for c in NUCLEOTIDE_CATEGORIES if c[0] in {"Ambiguity", ""}]
    assert len(ambiguity_rows) == 2
    assert sum(len(group) for _, group in ambiguity_rows) == 11
    # Expanded dialog stays within a usable height (scroll, not full color dump).
    assert dlg.height() < 900


def test_msa_plot_define_inputs_and_restore(qapp: QApplication, tmp_path: Path) -> None:
    fasta = tmp_path / "aln.fa"
    fasta.write_text(">s\nACGT\n")
    dlg = _msa_plot_dialog(qapp)
    dlg._input_define.setChecked(True)
    assert dlg.input_source() == ALIGN_INPUT_DEFINE
    assert dlg.apply_path_from_panel(fasta)
    assert dlg.run_files() == [fasta.resolve()]
    dlg._dpi.setValue(200)
    dlg._restore_defaults()
    assert dlg.input_source() == ALIGN_INPUT_SELECTED
    assert dlg._dpi.value() == DEFAULT_MSA_PLOT_DPI


def test_residue_color_grid_category_labels_align(qapp: QApplication) -> None:
    del qapp
    grid = ResidueColorGrid(
        None,
        colors={"A": "#111111", "R": "#222222", "N": "#333333"},
        categories=[
            ("Bases", ("A",)),
            ("Ambiguity", ("R",)),
            ("Other / gap", ("N",)),
        ],
    )
    labels = [
        lab
        for lab in grid.findChildren(QLabel)
        if lab.text() in {"Bases", "Ambiguity", "Other / gap"}
    ]
    assert len(labels) == 3
    widths = {lab.width() or lab.sizeHint().width() for lab in labels}
    assert all(lab.minimumWidth() == ResidueColorGrid.LABEL_WIDTH for lab in labels)
    assert all(lab.maximumWidth() == ResidueColorGrid.LABEL_WIDTH for lab in labels)
    assert len(widths) == 1 or all(
        lab.minimumWidth() == ResidueColorGrid.LABEL_WIDTH for lab in labels
    )
