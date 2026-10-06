"""TE-Aid parameters dialog: BLAST-style Inputs / Settings / Advanced / Options."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QApplication,
    QGroupBox,
    QLabel,
    QPushButton,
    QRadioButton,
    QToolButton,
)

from repeatviewer.job_options import (
    ALIGN_INPUT_DEFINE,
    ALIGN_INPUT_SELECTED,
    TeaAidOptions,
)
from repeatviewer.project import RepeatViewerProject
from repeatviewer.views.tool_dialogs import TeaAidParamsDialog


def _teaid_dialog(qapp: QApplication, **kw) -> TeaAidParamsDialog:
    del qapp
    kw.setdefault("initial", TeaAidOptions())
    kw.setdefault("project", RepeatViewerProject())
    dlg = TeaAidParamsDialog(**kw)
    dlg.show()
    dlg._resize_for_teaid_advanced(False)
    return dlg


def test_teaid_dialog_blast_style_sections(qapp: QApplication) -> None:
    dlg = _teaid_dialog(qapp)
    labels = [lab.text() for lab in dlg.findChildren(QLabel)]
    assert "<b>Mode</b>" in labels
    assert "<b>Inputs</b>" in labels
    assert "<b>Settings</b>" in labels
    assert not any(lab.text() == "Genome:" for lab in dlg.findChildren(QLabel))
    assert any(lab.text() == "Threads:" for lab in dlg.findChildren(QLabel))
    assert any(lab.text() == "Line transparency:" for lab in dlg.findChildren(QLabel))
    assert any(lab.text() == "Full-length hits (≥):" for lab in dlg.findChildren(QLabel))
    assert any(lab.text() == "Mid-length hits (≥):" for lab in dlg.findChildren(QLabel))
    assert any(lab.text() == "of sequence length" for lab in dlg.findChildren(QLabel))
    assert any(lab.text() == "but below full-length" for lab in dlg.findChildren(QLabel))
    assert not any(g.title() == "Plot labels" for g in dlg.findChildren(QGroupBox))
    assert not any(g.title() == "Genome BLAST hit colors" for g in dlg.findChildren(QGroupBox))
    texts = [b.text() for b in dlg.findChildren(QPushButton)]
    radios = [r.text() for r in dlg.findChildren(QRadioButton)]
    assert "Single genome" in radios
    assert "Pangenome" in radios
    assert "Use selected sequences" in radios
    assert "Manage Database" in texts
    assert "Add/Manage databases" not in texts
    assert any(t.startswith("Advanced Settings") for t in texts)
    assert dlg.width() == TeaAidParamsDialog._DIALOG_WIDTH
    assert dlg._evalue.width() == 260
    assert dlg._hit_high_pct.width() == 100
    assert dlg._hit_mid_pct.width() == 100
    assert dlg._teaid_advanced_panel.isHidden()
    assert dlg.input_source() == ALIGN_INPUT_SELECTED
    assert dlg._mode_single.isChecked()
    assert dlg._genome_row.isVisible()
    assert not dlg._pan_list.isVisible()
    # No "?" on the Current genome combo; mode help is on the radios.
    assert not any(
        child.text() == "?" for child in dlg._genome_row.findChildren(QToolButton)
    )
    pan_help = None
    for child in dlg.findChildren(QToolButton):
        if child.text() != "?":
            continue
        # Help next to the Pangenome radio.
        if child.parentWidget() is not None and dlg._mode_pan in child.parentWidget().findChildren(
            QRadioButton
        ):
            pan_help = child
            break
    assert pan_help is not None

    dlg._mode_pan.setChecked(True)
    qapp.processEvents()
    assert dlg._pan_list.isVisible()
    assert not dlg._genome_row.isVisible()
    assert "Add another genome" in texts or any(
        b.text() == "Add another genome" for b in dlg.findChildren(QPushButton)
    )
    # Pangenome uses genome dropdowns (not Name / Genome FASTA path rows).
    assert not any(lab.text() == "Name:" for lab in dlg.findChildren(QLabel))
    assert not any(lab.text() == "Genome FASTA:" for lab in dlg.findChildren(QLabel))
    assert len(dlg._pan_list._rows) >= 2
    opts = dlg.options()
    assert opts.use_pan_genome is True

    dlg._advanced_settings_btn.setChecked(True)
    qapp.processEvents()
    assert dlg._teaid_advanced_panel.isVisible()
    assert dlg._hit_high_pct.x() < dlg._hit_color_high_btn.x()
    assert dlg._hit_mid_pct.x() < dlg._hit_color_mid_btn.x()

def test_teaid_options_has_process_output(qapp: QApplication) -> None:
    dlg = _teaid_dialog(qapp)
    dlg._tabs.setCurrentIndex(1)
    qapp.processEvents()
    assert any(g.title() == "Process output" for g in dlg.findChildren(QGroupBox))
    assert not dlg._advanced_settings_btn.isVisible()
    dlg._show_process_window.setChecked(True)
    opts = dlg.options()
    assert opts.show_process_window is True
    assert opts.num_threads == dlg._threads.value()


def test_teaid_define_inputs_session_only(qapp: QApplication, tmp_path: Path) -> None:
    fasta = tmp_path / "a.fa"
    fasta.write_text(">s\nACGT\n")
    dlg = _teaid_dialog(qapp)
    dlg._input_define.setChecked(True)
    assert dlg.input_source() == ALIGN_INPUT_DEFINE
    assert dlg.apply_path_from_panel(fasta)
    assert dlg.run_files() == [fasta.resolve()]
    radios = [r.text() for r in dlg.findChildren(QRadioButton)]
    assert "Use selected sequences" in radios
    assert "Define inputs" in radios
