"""samtools parameters dialog: BLAST-style tabs / Inputs / footer."""

from __future__ import annotations

from PySide6.QtWidgets import QApplication, QGroupBox, QLabel, QPushButton, QTabBar

from repeatviewer.job_options import ALIGN_INPUT_DEFINE, SamtoolsOptions
from repeatviewer.views.tool_dialogs import (
    TOOL_PARAMS_DIALOG_WIDTH,
    TOOL_PARAMS_DIALOG_WIDTH_WIDE,
)
from repeatviewer.views.samtools_dialog import SamtoolsParamsDialog


def _dialog(qapp: QApplication, **kw) -> SamtoolsParamsDialog:
    del qapp
    kw.setdefault("initial", SamtoolsOptions())
    dlg = SamtoolsParamsDialog(**kw)
    dlg.show()
    return dlg


def test_samtools_dialog_tabs_and_footer(qapp: QApplication) -> None:
    dlg = _dialog(qapp)
    dlg._fit_to_content()
    qapp.processEvents()
    tabs = dlg.findChild(QTabBar)
    assert [tabs.tabText(i) for i in range(tabs.count())] == [
        "Sort",
        "Index",
        "Convert",
        "Options",
    ]
    labels = [lab.text() for lab in dlg.findChildren(QLabel)]
    assert "<b>Inputs</b>" in labels
    assert "<b>Settings</b>" in labels
    texts = [b.text() for b in dlg.findChildren(QPushButton)]
    assert [t for t in texts if t in ("Cancel", "Save", "Run")] == [
        "Cancel",
        "Save",
        "Run",
    ]
    assert "Restore Defaults" in texts
    assert dlg.width() == TOOL_PARAMS_DIALOG_WIDTH
    # No manual executable PATH picker — backend is automatic.
    assert not any(g.title() == "Executable" for g in dlg.findChildren(QGroupBox))


def test_samtools_dialog_define_inputs_widens(qapp: QApplication) -> None:
    dlg = _dialog(qapp)
    dlg._fit_to_content()
    qapp.processEvents()
    assert dlg.width() == TOOL_PARAMS_DIALOG_WIDTH
    dlg._input_define.setChecked(True)
    dlg._fit_to_content()
    qapp.processEvents()
    assert dlg.width() == TOOL_PARAMS_DIALOG_WIDTH_WIDE
    assert dlg.input_source() == ALIGN_INPUT_DEFINE


def test_samtools_dialog_height_follows_tab(qapp: QApplication) -> None:
    """Options → Index must shrink; blank band under Index Settings is a bug."""
    dlg = _dialog(qapp)
    dlg._tabs.setCurrentIndex(3)  # Options
    dlg._fit_to_content()
    qapp.processEvents()
    options_h = dlg.height()
    dlg._tabs.setCurrentIndex(1)  # Index
    dlg._fit_to_content()
    qapp.processEvents()
    index_h = dlg.height()
    assert index_h < options_h
    # Index is Inputs + short Settings — should stay compact.
    assert index_h < 420


def test_samtools_dialog_options_roundtrip(qapp: QApplication) -> None:
    dlg = _dialog(
        qapp,
        initial=SamtoolsOptions(
            operation="index",
            num_threads=8,
            index_after_sort=False,
            convert_format="sam",
        ),
    )
    dlg._tabs.setCurrentIndex(0)  # Sort
    qapp.processEvents()
    opts = dlg.options()
    assert opts.operation == "sort"
    assert opts.num_threads == 8
    assert opts.index_after_sort is False
    assert opts.convert_format == "sam"
    dlg._tabs.setCurrentIndex(1)  # Index — threads shared
    assert dlg._threads_index.value() == 8
    dlg._tabs.setCurrentIndex(2)  # Convert
    assert dlg.options().operation == "convert"
    assert dlg._threads_convert.value() == 8
    assert dlg._convert_output.currentData() == "sam"
    # Output format choices match samtools view (BAM / SAM / CRAM).
    assert [
        dlg._convert_output.itemData(i) for i in range(dlg._convert_output.count())
    ] == ["bam", "sam", "cram"]
    dlg._convert_output.setCurrentIndex(dlg._convert_output.findData("bam"))
    assert dlg.options().convert_format == "bam"
