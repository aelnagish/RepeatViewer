"""Dotplot parameters dialog: Built-in | Flexidot | Options tabs."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication, QCheckBox, QLabel, QPushButton, QRadioButton

from repeatviewer.job_options import (
    ALIGN_INPUT_DEFINE,
    ALIGN_INPUT_SELECTED,
    DOTPLOT_DIALOG_TAB_BUILTIN,
    DOTPLOT_DIALOG_TAB_FLEXIDOT,
    FLEXIDOT_RIGHT_CLICK_BUILTIN,
    FLEXIDOT_RIGHT_CLICK_FLEXIDOT,
    FlexidotOptions,
)
from repeatviewer.views.tool_dialogs import (
    DotplotParamsDialog,
    FlexidotParamsDialog,
)


def _flexidot_dialog(qapp: QApplication, **kw) -> FlexidotParamsDialog:
    del qapp
    kw.setdefault("initial", FlexidotOptions())
    dlg = FlexidotParamsDialog(**kw)
    dlg.show()
    dlg._fit_active_dotplot_tab()
    return dlg


def test_dotplot_dialog_has_builtin_and_flexidot_tabs(qapp: QApplication) -> None:
    dlg = _flexidot_dialog(qapp)
    assert dlg.windowTitle() == "Dotplot parameters"
    titles = [dlg._tabs.tabText(i) for i in range(dlg._tabs.count())]
    assert titles == ["Built-in", "Flexidot", "Options"]
    assert dlg._tabs.currentWidget() is dlg._builtin_tab
    assert dlg.active_engine_tab() == DOTPLOT_DIALOG_TAB_BUILTIN
    assert DotplotParamsDialog is FlexidotParamsDialog


def test_flexidot_dialog_align_style_sections(qapp: QApplication) -> None:
    dlg = _flexidot_dialog(qapp)
    dlg._tabs.setCurrentWidget(dlg._params_tab)
    qapp.processEvents()
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


def test_builtin_tab_mode_inputs_settings(qapp: QApplication) -> None:
    dlg = _flexidot_dialog(qapp)
    assert dlg._tabs.currentWidget() is dlg._builtin_tab
    assert any(lab.text() == "Word size:" for lab in dlg.findChildren(QLabel))
    assert any(
        isinstance(w, QCheckBox) and "reverse-complement" in w.text().lower()
        for w in dlg.findChildren(QCheckBox)
    )
    assert dlg._builtin_self.isChecked()
    assert "Self comparison" in [r.text() for r in dlg.findChildren(QRadioButton)]
    assert "Pairwise comparison" in [r.text() for r in dlg.findChildren(QRadioButton)]
    dlg._builtin_pairs.setChecked(True)
    dlg._builtin_word_size.setValue(12)
    dlg._builtin_include_reverse.setChecked(False)
    opts = dlg.options()
    assert opts.builtin_combined is True  # Pairwise
    assert opts.builtin_word_size == 12
    assert opts.builtin_include_reverse is False
    assert opts.last_dialog_tab == DOTPLOT_DIALOG_TAB_BUILTIN
    dlg._builtin_self.setChecked(True)
    assert dlg.options().builtin_combined is False


def test_flexidot_define_inputs_session_only(
    qapp: QApplication, tmp_path: Path
) -> None:
    fasta = tmp_path / "a.fa"
    fasta.write_text(">s\nACGT\n")
    dlg = _flexidot_dialog(qapp)
    dlg._tabs.setCurrentWidget(dlg._params_tab)
    qapp.processEvents()
    dlg._input_define.setChecked(True)
    assert dlg.input_source() == ALIGN_INPUT_DEFINE
    assert dlg._define_inputs.isVisible()
    assert dlg.apply_path_from_panel(fasta)
    assert dlg.run_files() == [fasta.resolve()]
    opts = dlg.options()
    assert opts.combined == dlg._combine.isChecked()
    assert opts.last_dialog_tab == DOTPLOT_DIALOG_TAB_FLEXIDOT
    assert "input_paths" not in opts.__dataclass_fields__


def test_flexidot_mode_self_pairwise_maps_to_options(qapp: QApplication) -> None:
    dlg = _flexidot_dialog(
        qapp,
        initial=FlexidotOptions(
            combined=True, last_dialog_tab=DOTPLOT_DIALOG_TAB_FLEXIDOT
        ),
    )
    assert dlg._tabs.currentWidget() is dlg._params_tab
    assert dlg._flexidot_pairs.isChecked()
    dlg._flexidot_self.setChecked(True)
    assert dlg.options().combined is False
    radios = [r.text() for r in dlg.findChildren(QRadioButton)]
    assert "Self comparison" in radios
    assert "Pairwise comparison" in radios
    assert "Use selected sequences" in radios
    assert "Define inputs" in radios


def test_dotplot_right_click_builtin_or_flexidot(qapp: QApplication) -> None:
    dlg = _flexidot_dialog(
        qapp,
        initial=FlexidotOptions(right_click_action=FLEXIDOT_RIGHT_CLICK_FLEXIDOT),
    )
    labels = [r.text() for r in dlg._right_click_buttons.values()]
    assert "Run Built-in Dotplot" in labels
    assert "Run Flexidot" in labels
    assert "Open this dialog" in labels
    assert dlg._right_click_buttons[FLEXIDOT_RIGHT_CLICK_FLEXIDOT].isChecked()
    dlg._right_click_buttons[FLEXIDOT_RIGHT_CLICK_BUILTIN].setChecked(True)
    assert dlg.options().right_click_action == FLEXIDOT_RIGHT_CLICK_BUILTIN


def test_flexidot_options_process_output_stays_compact(qapp: QApplication) -> None:
    dlg = _flexidot_dialog(qapp)
    dlg._fit_active_dotplot_tab()
    qapp.processEvents()
    dlg._tabs.setCurrentWidget(dlg._params_tab)
    dlg._resize_for_flexidot_advanced(False)
    qapp.processEvents()
    dlg._tabs.setCurrentWidget(dlg._options_tab)
    qapp.processEvents()
    show = dlg._show_process_window
    close = dlg._close_process_window
    gap = close.mapTo(dlg, close.rect().topLeft()).y() - show.mapTo(
        dlg, show.rect().bottomLeft()
    ).y()
    assert gap <= 16
