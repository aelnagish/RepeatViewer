"""Actions → Duplicate window opens a floating Preview that ignores file-list clicks."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication

from repeatviewer.widgets.preview_duplicate_window import PreviewDuplicateWindow
from repeatviewer.widgets.preview_split_window import PreviewSplitWindow
from repeatviewer.widgets.preview_tabs import FilePreviewPage, PreviewTabWidget


def test_split_window_opens_left_and_right_panes(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "split.fa"
    path.write_text(">s\nACGTACGT\n", encoding="utf-8")

    tabs = PreviewTabWidget()
    assert tabs.open_file(path, replace=False) is True
    assert tabs.open_split_window() is True
    assert len(tabs._duplicate_windows) == 1

    win = tabs._duplicate_windows[0]
    assert isinstance(win, PreviewSplitWindow)
    assert win._left is not None and win._right is not None
    left_page = win._left._tabs.currentWidget()
    right_page = win._right._tabs.currentWidget()
    assert isinstance(left_page, FilePreviewPage)
    assert isinstance(right_page, FilePreviewPage)
    assert left_page.path == path.resolve()
    assert right_page.path == path.resolve()
    assert win.isVisible()

    win.close()
    qapp.processEvents()
    assert tabs._duplicate_windows == []
    tabs.close()


def test_split_window_allows_narrow_size_and_free_ratio(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "split_ratio.fa"
    path.write_text(">a\nACGT\n>b\nTTTT\n", encoding="utf-8")

    tabs = PreviewTabWidget()
    assert tabs.open_file(path, replace=False) is True
    assert tabs.open_split_window() is True
    win = tabs._duplicate_windows[0]
    assert isinstance(win, PreviewSplitWindow)

    assert win.minimumWidth() <= 320
    assert win._left.minimumWidth() <= 120
    assert win._right.minimumWidth() <= 120

    win.resize(700, 500)
    qapp.processEvents()
    win._splitter.setSizes([180, 520])
    qapp.processEvents()
    sizes = win._splitter.sizes()
    assert sizes[0] < sizes[1]
    assert sizes[0] >= win._left.minimumWidth() - 2

    win._splitter.setSizes([520, 180])
    qapp.processEvents()
    sizes = win._splitter.sizes()
    assert sizes[1] < sizes[0]

    win.close()
    qapp.processEvents()
    tabs.close()


def test_annotation_zoom_bar_overflow_hides_when_narrow(
    qapp: QApplication,
) -> None:
    from repeatviewer.widgets.preview_controls import AnnotationZoomBar

    bar = AnnotationZoomBar()
    bar.set_overview_mode(True)
    bar.set_msa_mode(True)
    bar.set_sequence_count(12)
    bar.resize(900, bar.height())
    bar.show()
    qapp.processEvents()
    bar._overview_overflow.sync()
    assert not bar._consensus_btn.isHidden()

    bar.resize(200, bar.height())
    qapp.processEvents()
    bar._overview_overflow.sync()
    # Right-side units drop first; −/+ remain.
    assert bar._consensus_btn.isHidden()
    assert bar._clean_btn.isHidden()
    zoom_out = bar._overview_host.layout().itemAt(0).widget()
    assert zoom_out is not None and not zoom_out.isHidden()

    bar.close()


def test_annotation_zoom_bar_overflow_uses_narrow_map_slot(
    qapp: QApplication,
) -> None:
    """Split overview: measure the map slot, not the preferred toolbar width."""
    from PySide6.QtWidgets import QHBoxLayout, QLabel, QVBoxLayout, QWidget

    from repeatviewer.widgets.preview_controls import AnnotationZoomBar

    pane = QWidget()
    pane.setFixedWidth(420)
    row = QHBoxLayout(pane)
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(0)
    map_host = QWidget()
    map_host.setMinimumWidth(0)
    map_col = QVBoxLayout(map_host)
    map_col.setContentsMargins(0, 0, 0, 0)
    map_col.setSpacing(0)
    slot = QWidget()
    slot.setFixedHeight(26)
    slot_layout = QHBoxLayout(slot)
    slot_layout.setContentsMargins(0, 0, 0, 0)
    bar = AnnotationZoomBar()
    bar.set_overview_mode(True)
    bar.set_msa_mode(True)
    slot_layout.addWidget(bar, 1)
    map_col.addWidget(slot)
    map_col.addWidget(QLabel("scale"))
    rail = QLabel("Display")
    rail.setFixedWidth(24)
    row.addWidget(map_host, 1)
    row.addWidget(rail, 0)
    pane.show()
    qapp.processEvents()

    hider = bar._overview_overflow
    hider.add_width_source(slot)
    hider.add_width_source(map_host)
    hider.sync()
    qapp.processEvents()

    assert map_host.width() < 420
    # Right-to-left: Cons/Clean gone; Fit usually gone too in a ~400px map column.
    assert bar._consensus_btn.isHidden()
    assert bar._clean_btn.isHidden()
    assert bar._fit_btn.isHidden()
    assert hider._available_width() <= map_host.width()
    pane.close()


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
