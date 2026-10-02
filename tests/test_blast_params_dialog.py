"""BLAST parameters dialog: Align-style Inputs / Settings / Advanced."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication, QGroupBox, QLabel, QPushButton, QTabBar

from repeatviewer.job_options import (
    ALIGN_INPUT_DEFINE,
    ALIGN_INPUT_SELECTED,
    BlastOptions,
    load_blast_options,
    save_blast_options,
)
from repeatviewer.project import RepeatViewerProject
from repeatviewer.views.tool_dialogs import BlastParamsDialog


def _blast_dialog(qapp: QApplication, **kw) -> BlastParamsDialog:
    del qapp
    kw.setdefault("project", RepeatViewerProject())
    kw.setdefault("initial", BlastOptions())
    dlg = BlastParamsDialog(**kw)
    dlg.show()
    return dlg


def test_blast_dialog_align_style_sections_and_footer(qapp: QApplication) -> None:
    dlg = _blast_dialog(qapp)
    dlg._resize_for_blast_advanced(False)
    qapp.processEvents()
    tabs = dlg.findChild(QTabBar)
    assert [tabs.tabText(i) for i in range(tabs.count())] == [
        "BLASTn",
        "BLASTp",
        "BLASTx",
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
    assert any(t.startswith("Advanced Settings") for t in texts)
    assert not any(g.title() == "Post-processing" for g in dlg.findChildren(QGroupBox))
    assert any("Max sequence number" in lab.text() for lab in dlg.findChildren(QLabel))
    assert any("Top longest sequences" in lab.text() for lab in dlg.findChildren(QLabel))
    assert not any("Max MSA lines" in lab.text() for lab in dlg.findChildren(QLabel))
    assert not any("Top longest lines" in lab.text() for lab in dlg.findChildren(QLabel))
    assert any(t == "Manage BLAST databases" for t in texts)
    assert dlg.width() == BlastParamsDialog._DIALOG_WIDTH
    assert dlg.width() == 608
    # Inputs radio dots share the Settings field column (Database combo).
    dlg._sync_blast_inputs_field_column()
    qapp.processEvents()
    radio_x = dlg._input_selected.mapTo(dlg, dlg._input_selected.rect().topLeft()).x()
    db_x = dlg._blastn_database.mapTo(dlg, dlg._blastn_database.rect().topLeft()).x()
    assert abs(radio_x - db_x) <= 8
    # Compact open: footer sits under Settings (no tall blank band).
    content_bottom = dlg._stack.geometry().bottom()
    buttons_top = dlg._button_box.geometry().top()
    assert buttons_top - content_bottom <= 20
    # Shared Inputs panel stays put; visible on program tabs, hidden on Options.
    dlg._tabs.setCurrentIndex(0)
    assert dlg._inputs_panel.isVisible()
    assert dlg._inputs_radios_box.parent() is dlg._inputs_panel
    assert dlg._define_inputs.parent() is dlg._inputs_panel
    dlg._tabs.setCurrentIndex(1)
    assert dlg._inputs_panel.isVisible()
    assert dlg._inputs_radios_box.parent() is dlg._inputs_panel
    dlg._tabs.setCurrentIndex(3)
    assert dlg._inputs_panel.isHidden()
    assert dlg._advanced_settings_btn.isHidden()


def test_blast_options_process_output_stays_compact_after_program_tab(
    qapp: QApplication,
) -> None:
    """BLASTn → Options must not stretch Process output checkboxes apart."""
    dlg = _blast_dialog(qapp)
    dlg._resize_for_blast_advanced(False)
    qapp.processEvents()
    dlg._tabs.setCurrentIndex(0)
    dlg._resize_for_blast_advanced(False)
    qapp.processEvents()
    dlg._tabs.setCurrentIndex(3)
    qapp.processEvents()
    show = dlg._show_process_window
    close = dlg._close_process_window
    gap = close.mapTo(dlg, close.rect().topLeft()).y() - show.mapTo(
        dlg, show.rect().bottomLeft()
    ).y()
    assert gap <= 16
    # Process box itself should stay near its sizeHint (not dialog-tall).
    box = show.parentWidget()
    assert box is not None
    assert box.height() <= box.sizeHint().height() + 24


def test_blast_define_inputs_are_session_only_and_fill_from_panel(
    qapp: QApplication, tmp_path: Path, monkeypatch
) -> None:
    monkeypatch.setattr(
        "repeatviewer.job_options._load_settings_group",
        lambda _g: {
            "max_msa_lines": 250,
            "top_msa_lines": 80,
        },
    )
    loaded = load_blast_options()
    assert loaded.input_source == ALIGN_INPUT_SELECTED
    assert loaded.input_paths == []
    assert loaded.max_msa_lines == 250

    q1 = tmp_path / "q1.fa"
    q2 = tmp_path / "q2.fa"
    q1.write_text(">a\nACGT\n")
    q2.write_text(">b\nTTTT\n")
    dlg = _blast_dialog(qapp, initial=loaded)
    assert dlg.input_source() == ALIGN_INPUT_SELECTED
    assert dlg.apply_path_from_panel(q1) is False

    dlg._input_define.setChecked(True)
    assert dlg.input_source() == ALIGN_INPUT_DEFINE
    assert dlg._define_inputs.isVisible()
    assert dlg.apply_path_from_panel(q1) is True
    assert dlg.apply_folder_from_panel(tmp_path) is True
    run_files = dlg.run_files()
    assert q1.resolve() in run_files
    assert q2.resolve() in run_files

    # Saved options never keep Define-input paths.
    opts = dlg.options()
    assert opts.input_source == ALIGN_INPUT_SELECTED
    assert opts.input_paths == []
    assert opts.max_msa_lines == 250

    saved: dict = {}

    def _capture(_group: str, values: dict) -> None:
        saved.update(values)

    monkeypatch.setattr("repeatviewer.job_options._save_settings_group", _capture)
    save_blast_options(opts)
    assert "input_source" not in saved
    assert "input_paths" not in saved
    assert saved["max_msa_lines"] == 250


def test_blast_advanced_settings_toggle(qapp: QApplication) -> None:
    dlg = _blast_dialog(qapp)
    assert dlg._blastn_advanced.isHidden()
    dlg._advanced_settings_btn.setChecked(True)
    assert not dlg._blastn_advanced.isHidden()
    assert "▾" in dlg._advanced_settings_btn.text()
    dlg._tabs.setCurrentIndex(1)
    assert dlg._advanced_settings_btn.isChecked() is False
    assert dlg._blastp_advanced.isHidden()
