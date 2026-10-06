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
        "Online BLAST",
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
    assert any(t == "Manage Database" for t in texts)
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
    # Online BLAST: Inputs visible, no Advanced Settings.
    dlg._tabs.setCurrentIndex(3)
    assert dlg._inputs_panel.isVisible()
    assert dlg._advanced_settings_btn.isHidden()
    dlg._tabs.setCurrentIndex(4)
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
    dlg._tabs.setCurrentIndex(4)
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


def test_online_blast_tab_options_round_trip(qapp: QApplication, monkeypatch) -> None:
    from repeatviewer.job_options import (
        BLAST_PROGRAM_BLASTP,
        BLAST_PROGRAM_ONLINE,
        BlastOptions,
        DIALOG_TAB_OPTIONS,
    )

    dlg = _blast_dialog(
        qapp,
        initial=BlastOptions(
            program=BLAST_PROGRAM_ONLINE,
            last_dialog_tab=BLAST_PROGRAM_ONLINE,
            online_program=BLAST_PROGRAM_BLASTP,
            online_database="swissprot",
            online_e_value=0.01,
            online_hitlist_size=50,
            online_write_text=False,
        ),
    )
    assert dlg._tabs.currentIndex() == 3
    assert dlg._advanced_settings_btn.isHidden()
    assert dlg._inputs_panel.isVisible()
    opts = dlg.options()
    assert opts.program == BLAST_PROGRAM_ONLINE
    assert opts.last_dialog_tab == BLAST_PROGRAM_ONLINE
    assert opts.online_program == BLAST_PROGRAM_BLASTP
    assert opts.online_database == "swissprot"
    assert opts.online_e_value == 0.01
    assert opts.online_hitlist_size == 50
    assert opts.online_write_text is False
    assert opts.active_database_name is None
    assert opts.is_online is True

    saved: dict = {}
    monkeypatch.setattr(
        "repeatviewer.job_options._save_settings_group",
        lambda _g, values: saved.update(values),
    )
    save_blast_options(opts)
    assert saved["program"] == BLAST_PROGRAM_ONLINE
    assert saved["online_program"] == BLAST_PROGRAM_BLASTP
    assert saved["online_database"] == "swissprot"
    assert saved["last_dialog_tab"] == BLAST_PROGRAM_ONLINE


def test_blast_dialog_reopens_last_tab(qapp: QApplication) -> None:
    from repeatviewer.job_options import (
        BLAST_PROGRAM_BLASTN,
        BLAST_PROGRAM_BLASTX,
        BlastOptions,
        DIALOG_TAB_OPTIONS,
    )

    dlg = _blast_dialog(
        qapp,
        initial=BlastOptions(
            program=BLAST_PROGRAM_BLASTN,
            last_dialog_tab=BLAST_PROGRAM_BLASTX,
        ),
    )
    assert dlg._tabs.currentIndex() == 2
    assert dlg.current_dialog_tab() == BLAST_PROGRAM_BLASTX

    dlg._tabs.setCurrentIndex(4)
    assert dlg.current_dialog_tab() == DIALOG_TAB_OPTIONS
    assert dlg.options().last_dialog_tab == BLAST_PROGRAM_BLASTN
    # Options tab does not change the run program.
    assert dlg.options().program == BLAST_PROGRAM_BLASTN


def test_blast_and_align_share_inputs_indent(qapp: QApplication) -> None:
    """BLAST, MAFFT, and FastGA Inputs radios share one label-column width."""
    from PySide6.QtCore import Qt

    from repeatviewer.views.tool_dialogs import (
        AlignParamsDialog,
        TOOL_PARAMS_LABEL_COLUMN_REF,
    )

    blast = _blast_dialog(qapp)
    blast._sync_blast_inputs_field_column()
    align = AlignParamsDialog(parent=None)
    align._sync_align_field_columns()
    qapp.processEvents()
    expected = blast.fontMetrics().horizontalAdvance(TOOL_PARAMS_LABEL_COLUMN_REF)
    assert blast._inputs_radios_label.width() == expected
    assert align._mafft_inputs_radios_label.width() == expected
    assert align._mafft_mode_label.width() == expected
    # FastGA must not widen past MAFFT / BLAST (previously used a longer ref_label).
    panel = align.fastga_panel()
    panel._sync_field_column()
    assert panel._inputs_radios_label.width() == expected
    assert panel._mode_label.width() == expected
    # Settings labels are right-aligned in that column so short text meets the field.
    strategy_lab = None
    form = align._mafft_basic.layout()
    for row in range(form.rowCount()):
        item = form.itemAt(row, form.ItemRole.LabelRole)
        w = item.widget() if item is not None else None
        if isinstance(w, QLabel) and w.text() == "Strategy:":
            strategy_lab = w
            break
    assert strategy_lab is not None
    assert strategy_lab.width() == expected
    assert int(strategy_lab.alignment()) & int(Qt.AlignmentFlag.AlignRight)

