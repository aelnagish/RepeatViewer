"""Flexidot parameters dialog: Align-style Mode / Inputs / Outputs / Settings."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication, QCheckBox, QLabel, QPushButton, QRadioButton

from repeatviewer.job_options import (
    ALIGN_INPUT_DEFINE,
    ALIGN_INPUT_SELECTED,
    FlexidotOptions,
)
from repeatviewer.views.tool_dialogs import FlexidotParamsDialog


def _flexidot_dialog(qapp: QApplication, **kw) -> FlexidotParamsDialog:
    del qapp
    kw.setdefault("initial", FlexidotOptions())
    dlg = FlexidotParamsDialog(**kw)
    dlg.show()
    dlg._resize_for_flexidot_advanced(False)
    return dlg


def test_flexidot_dialog_align_style_sections(qapp: QApplication) -> None:
    dlg = _flexidot_dialog(qapp)
    labels = [lab.text() for lab in dlg.findChildren(QLabel)]
    assert "<b>Mode</b>" in labels
    assert "<b>Inputs</b>" in labels
    assert "<b>Outputs</b>" in labels
    assert "<b>Settings</b>" in labels
    assert any(lab.text() == "k-mer:" for lab in dlg.findChildren(QLabel))
    assert any(lab.text() == "Sequence type:" for lab in dlg.findChildren(QLabel))
    assert any(
        isinstance(w, QCheckBox) and w.text() == "Show GFF features on dotplot"
        for w in dlg.findChildren(QCheckBox)
    )
    texts = [b.text() for b in dlg.findChildren(QPushButton)]
    assert any(t.startswith("Advanced Settings") for t in texts)
    assert dlg.width() == FlexidotParamsDialog._DIALOG_WIDTH
    assert dlg._flexidot_advanced_panel.isHidden()
    assert dlg.input_source() == ALIGN_INPUT_SELECTED


def test_flexidot_define_inputs_session_only(
    qapp: QApplication, tmp_path: Path
) -> None:
    fasta = tmp_path / "a.fa"
    fasta.write_text(">s\nACGT\n")
    dlg = _flexidot_dialog(qapp)
    dlg._input_define.setChecked(True)
    assert dlg.input_source() == ALIGN_INPUT_DEFINE
    assert dlg._define_inputs.isVisible()
    assert dlg.apply_path_from_panel(fasta)
    assert dlg.run_files() == [fasta.resolve()]
    opts = dlg.options()
    assert opts.combined == dlg._combine.isChecked()
    assert "input_paths" not in opts.__dataclass_fields__


def test_flexidot_mode_combine_maps_to_options(qapp: QApplication) -> None:
    dlg = _flexidot_dialog(qapp, initial=FlexidotOptions(combined=True))
    assert dlg._combine.isChecked()
    dlg._separate.setChecked(True)
    assert dlg.options().combined is False
    radios = [r.text() for r in dlg.findChildren(QRadioButton)]
    assert "Use selected sequences" in radios
    assert "Define inputs" in radios


def test_flexidot_options_process_output_stays_compact(qapp: QApplication) -> None:
    dlg = _flexidot_dialog(qapp)
    dlg._resize_for_flexidot_advanced(False)
    qapp.processEvents()
    dlg._tabs.setCurrentIndex(0)
    dlg._resize_for_flexidot_advanced(False)
    qapp.processEvents()
    dlg._tabs.setCurrentIndex(1)
    qapp.processEvents()
    show = dlg._show_process_window
    close = dlg._close_process_window
    gap = close.mapTo(dlg, close.rect().topLeft()).y() - show.mapTo(
        dlg, show.rect().bottomLeft()
    ).y()
    assert gap <= 16
