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

    # Floating Overview/Annotation jobs need a local progress API (other screen).
    assert hasattr(win, "_show_job_progress")
    assert hasattr(win, "_on_job_progress")
    assert hasattr(win, "_hide_job_progress")
    assert win._task_progress is not None
    win._show_job_progress("Similarity annotation", 3)
    assert win._job_running is True
    assert win._status_bar.isVisible()
    win._on_job_progress(1, "1/3 seqA")
    win._hide_job_progress()
    assert win._job_running is False
    assert win._status_bar.isHidden()

    win.close()
    qapp.processEvents()
    assert tabs._duplicate_windows == []
    tabs.close()


def test_iter_job_progress_hosts_includes_main_and_duplicate(
    qapp: QApplication,
) -> None:
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QLabel, QMainWindow, QWidget

    from repeatviewer.widgets.task_progress import iter_job_progress_hosts

    class _Main(QMainWindow):
        def _show_job_progress(self, *_a, **_k) -> None:
            return None

        def _set_job_actions_enabled(self, *_a) -> None:
            return None

    main = _Main()
    tool = QWidget(
        main,
        Qt.WindowType.Tool | Qt.WindowType.WindowTitleHint,
    )

    def _show(*_a, **_k) -> None:
        return None

    tool._show_job_progress = _show  # type: ignore[attr-defined]
    inner = QLabel("x", tool)
    hosts = iter_job_progress_hosts(inner)
    assert any(host is tool for host in hosts)
    assert any(host is main for host in hosts)
    main.close()
    tool.close()
