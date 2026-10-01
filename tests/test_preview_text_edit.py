"""Preview > Edit text: Save text state, silent reloads and the unsaved-changes prompt."""

from __future__ import annotations

from pathlib import Path

import pytest


def _process(qapp, n: int = 20) -> None:
    import time

    for _ in range(n):
        qapp.processEvents()
        time.sleep(0.005)


@pytest.fixture
def preview(qapp, tmp_path):
    from repeatviewer.widgets.preview_controls import PreviewZoomBar
    from repeatviewer.widgets.preview_tabs import PreviewTabWidget

    a = tmp_path / "A.fasta"
    a.write_text(">A\nACGT\n")
    b = tmp_path / "B.fasta"
    b.write_text(">B\nTTTT\n")
    bar = PreviewZoomBar()
    tabs = PreviewTabWidget(zoom_bar=bar)
    tabs.show()
    tabs.open_file(a, replace=True)
    _process(qapp)
    return tabs, bar, a, b


def _start_editing(qapp, tabs, bar, text: str):
    from PySide6.QtTest import QTest

    from repeatviewer.widgets.preview_tabs import DefaultFastaPreviewTab

    page = tabs.currentWidget()
    page._activate_fasta_tab(DefaultFastaPreviewTab.TEXT)
    bar._edit_text_action.trigger()
    editor = page._active_text_editor()
    cursor = editor.textCursor()
    cursor.movePosition(cursor.MoveOperation.End)
    editor.setTextCursor(cursor)
    QTest.keyClicks(editor, text)
    _process(qapp, 5)
    return page, editor


def _answer_prompt(monkeypatch, choice: str) -> list[str]:
    from PySide6.QtWidgets import QMessageBox

    asked: list[str] = []

    def fake_exec(box):
        if box.windowTitle() == "Unsaved Changes":
            asked.append(box.text())
            box.button(getattr(QMessageBox.StandardButton, choice)).click()
        return 0

    monkeypatch.setattr(QMessageBox, "exec", fake_exec)
    return asked


def test_save_text_turns_on_after_typing_and_saves(qapp, preview):
    tabs, bar, a, _b = preview
    page, _editor = _start_editing(qapp, tabs, bar, "GG")
    assert bar._save_text_action.isEnabled()
    bar._save_text_action.trigger()
    _process(qapp, 5)
    assert a.read_text().endswith("GG")
    assert not page.has_unsaved_text()
    assert not bar._save_text_action.isEnabled()


def test_silent_reload_keeps_unsaved_edits(qapp, preview):
    tabs, bar, a, _b = preview
    page, editor = _start_editing(qapp, tabs, bar, "GG")
    page.reload_from_disk()                 # folder watcher / same file clicked again
    tabs.open_file(a, replace=True)
    assert page.has_unsaved_text() and editor.toPlainText().endswith("GG")


def test_undo_back_to_saved_text_is_not_unsaved(qapp, preview):
    tabs, bar, _a, _b = preview
    page, editor = _start_editing(qapp, tabs, bar, "G")
    editor.undo()
    assert not page.has_unsaved_text()


def test_cancel_keeps_the_edited_file_open(qapp, preview, monkeypatch):
    tabs, bar, a, b = preview
    page, _editor = _start_editing(qapp, tabs, bar, "GG")
    asked = _answer_prompt(monkeypatch, "Cancel")
    assert tabs.open_file(b, replace=True) is False
    assert asked == ['Save changes to "A.fasta"?']
    assert tabs.currentWidget().path == a.resolve() and page.has_unsaved_text()


def test_dont_save_switches_and_leaves_the_file_alone(qapp, preview, monkeypatch):
    tabs, bar, a, b = preview
    _start_editing(qapp, tabs, bar, "GG")
    _answer_prompt(monkeypatch, "Discard")
    assert tabs.open_file(b, replace=True) is True
    assert tabs.currentWidget().path == b.resolve()
    assert a.read_text() == ">A\nACGT\n"


def test_save_at_the_prompt_writes_the_file(qapp, preview, monkeypatch):
    tabs, bar, a, b = preview
    _start_editing(qapp, tabs, bar, "CC")
    _answer_prompt(monkeypatch, "Save")
    assert tabs.open_file(b, replace=True) is True
    assert a.read_text().endswith("CC")


def test_closing_the_tab_asks_first(qapp, preview, monkeypatch):
    tabs, bar, _a, _b = preview
    page, _editor = _start_editing(qapp, tabs, bar, "GG")
    asked = _answer_prompt(monkeypatch, "Cancel")
    tabs._close_tab(tabs.indexOf(page))
    assert asked and tabs.indexOf(page) >= 0
