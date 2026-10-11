"""Actions → Duplicate window opens a floating Preview that ignores file-list clicks."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication

from repeatviewer.widgets.preview_duplicate_window import PreviewDuplicateWindow
from repeatviewer.widgets.preview_tabs import FilePreviewPage, PreviewTabWidget


def test_duplicate_window_stays_when_main_preview_changes(
    qapp: QApplication, tmp_path: Path
) -> None:
    a = tmp_path / "a.fa"
    a.write_text(">a\nACGTACGT\n", encoding="utf-8")
    b = tmp_path / "b.fa"
    b.write_text(">b\nTTTTAAAA\n", encoding="utf-8")

    tabs = PreviewTabWidget()
    assert tabs.open_file(a, replace=False) is True
    assert tabs.count() == 1
    assert tabs.open_duplicate_window() is True
    assert len(tabs._duplicate_windows) == 1

    win = tabs._duplicate_windows[0]
    assert isinstance(win, PreviewDuplicateWindow)
    win_page = win._tabs.currentWidget()
    assert isinstance(win_page, FilePreviewPage)
    assert win_page.path == a.resolve()

    # Main Preview follows selection; the floating window stays on A.
    assert tabs.open_file(b, replace=True) is True
    qapp.processEvents()
    main_page = tabs.currentWidget()
    assert isinstance(main_page, FilePreviewPage)
    assert main_page.path == b.resolve()
    assert win_page.path == a.resolve()
    assert win.isVisible()

    win.close()
    qapp.processEvents()
    assert tabs._duplicate_windows == []
    tabs.close()
