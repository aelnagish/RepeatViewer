"""Actions → Duplicate / Split Preview behaviors."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QApplication, QVBoxLayout, QWidget

from repeatviewer.widgets.preview_duplicate_window import PreviewDuplicateWindow
from repeatviewer.widgets.preview_split_window import PreviewSplitWindow
from repeatviewer.widgets.preview_tabs import FilePreviewPage, PreviewTabWidget


def test_split_window_opens_inline_in_preview_host(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "split.fa"
    path.write_text(">s\nACGTACGT\n", encoding="utf-8")

    host = QWidget()
    layout = QVBoxLayout(host)
    tabs = PreviewTabWidget()
    layout.addWidget(tabs, stretch=1)
    host.resize(900, 500)
    host.show()
    qapp.processEvents()

    assert tabs.open_file(path, replace=False) is True
    assert tabs.open_split_window() is True
    assert tabs.is_split_active()
    assert tabs.isHidden()
    split = tabs._inline_split
    assert isinstance(split, PreviewSplitWindow)
    assert split._floating is False
    left_page = split._left._tabs.currentWidget()
    right_page = split._right._tabs.currentWidget()
    assert isinstance(left_page, FilePreviewPage)
    assert isinstance(right_page, FilePreviewPage)
    assert left_page.path == path.resolve()
    assert right_page.path == path.resolve()

    other = tmp_path / "other.fa"
    other.write_text(">t\nTTTTAAAA\n", encoding="utf-8")
    assert tabs.open_file(other, replace=True) is True
    assert tabs.is_split_active()
    assert split.current_path() == other.resolve()

    assert tabs.open_split_window() is True  # toggle off
    qapp.processEvents()
    assert not tabs.is_split_active()
    assert not tabs.isHidden()
    from PySide6.QtTest import QTest

    QTest.qWait(300)
    host.close()
    tabs.close()


def test_split_window_allows_narrow_size_and_free_ratio(
    qapp: QApplication, tmp_path: Path
) -> None:
    path = tmp_path / "split_ratio.fa"
    path.write_text(">a\nACGT\n>b\nTTTT\n", encoding="utf-8")

    host = QWidget()
    layout = QVBoxLayout(host)
    tabs = PreviewTabWidget()
    layout.addWidget(tabs, stretch=1)
    host.resize(700, 500)
    host.show()
    from PySide6.QtTest import QTest

    QTest.qWait(50)

    assert tabs.open_file(path, replace=False) is True
    assert tabs.open_split_window() is True
    win = tabs._inline_split
    assert isinstance(win, PreviewSplitWindow)

    assert win._left.minimumWidth() <= 120
    assert win._right.minimumWidth() <= 120

    QTest.qWait(50)
    win._splitter.setSizes([180, 520])
    QTest.qWait(50)
    sizes = win._splitter.sizes()
    assert sizes[0] < sizes[1]
    assert sizes[0] >= win._left.minimumWidth() - 2

    win._splitter.setSizes([520, 180])
    QTest.qWait(50)
    sizes = win._splitter.sizes()
    assert sizes[1] < sizes[0]

    tabs._close_inline_split()
    QTest.qWait(300)
    host.close()
    tabs.close()


def test_split_stays_two_panes_when_opening_another_file(
    qapp: QApplication, tmp_path: Path
) -> None:
    """File-list click must not re-show the single Overview above the split."""
    from PySide6.QtTest import QTest

    a = tmp_path / "a.fa"
    a.write_text(">s1\nACGTACGTACGTACGT\n>s2\nTTTTAAAACCCCGGGG\n", encoding="utf-8")
    b = tmp_path / "b.fa"
    b.write_text(">t1\nGGGGCCCCAAAATTTT\n>t2\nACGTACGTACGTACGT\n", encoding="utf-8")

    host = QWidget()
    layout = QVBoxLayout(host)
    tabs = PreviewTabWidget()
    layout.addWidget(tabs, stretch=1)
    host.resize(1000, 600)
    host.show()
    QTest.qWait(50)

    assert tabs.open_file(a, replace=False) is True
    assert tabs.open_split_window() is True
    assert tabs.is_split_active()
    assert tabs.isHidden()

    assert tabs.open_file(b, replace=True) is True
    # Mimic MainWindow._ensure_preview_visible after a file-list click.
    if tabs.is_split_active():
        tabs.keep_split_chrome_hidden()
    else:
        tabs.show()
    QTest.qWait(200)

    assert tabs.is_split_active()
    assert tabs.isHidden()
    header = tabs._header_bar
    if header is not None:
        assert header.isHidden()
    split = tabs._inline_split
    assert isinstance(split, PreviewSplitWindow)
    assert split.isVisible()
    assert split.current_path() == b.resolve()

    tabs._close_inline_split()
    QTest.qWait(300)
    host.close()
    tabs.close()


def test_open_additional_tab_restores_split_on_new_tab(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Open in Additional Tab adds a tab and re-opens Overview|Dotplot on it."""
    from PySide6.QtTest import QTest

    a = tmp_path / "a.fa"
    a.write_text(
        ">s1\nACGTACGTACGTACGT\n>s2\nTTTTAAAACCCCGGGG\n", encoding="utf-8"
    )
    b = tmp_path / "b.fa"
    b.write_text(
        ">t1\nGGGGCCCCAAAATTTT\n>t2\nACGTACGTACGTACGT\n", encoding="utf-8"
    )

    host = QWidget()
    layout = QVBoxLayout(host)
    tabs = PreviewTabWidget()
    layout.addWidget(tabs, stretch=1)
    host.resize(1000, 600)
    host.show()
    QTest.qWait(50)

    assert tabs.open_file(a, replace=False) is True
    assert tabs.open_split_with_tabs(1, 2) is True
    QTest.qWait(500)
    assert tabs.is_split_active()
    assert tabs.isHidden()

    assert tabs.open_file(b, replace=False) is True
    QTest.qWait(500)

    assert tabs.count() == 2
    current = tabs.currentWidget()
    assert isinstance(current, FilePreviewPage)
    assert current.path == b.resolve()
    assert tabs.is_split_active()
    assert tabs.isHidden()
    split = tabs._inline_split
    assert isinstance(split, PreviewSplitWindow)
    assert split.current_path() == b.resolve()
    QTest.qWait(200)
    assert split._left.current_fasta_tab() == 1
    assert split._right.current_fasta_tab() == 2
    left_page = split._left.preview_page()
    right_page = split._right.preview_page()
    assert isinstance(left_page, FilePreviewPage)
    assert isinstance(right_page, FilePreviewPage)
    assert left_page.is_annotation_preview_active()
    assert right_page.is_dotplot_preview_active()

    # Switching Preview tabs while split loads that file into both panes.
    tabs.setCurrentIndex(0)
    QTest.qWait(500)
    assert tabs._inline_split.current_path() == a.resolve()

    tabs._close_inline_split()
    QTest.qWait(200)
    assert tabs.count() == 2
    assert not tabs.isHidden()

    host.close()
    tabs.close()


def test_msa_split_layout_remembered_for_next_msa(
    qapp: QApplication, tmp_path: Path, monkeypatch
) -> None:
    """Overview|Dotplot split on an MSA is restored when opening another MSA."""
    from PySide6.QtCore import QSettings
    from PySide6.QtTest import QTest

    from repeatviewer.widgets import preview_zoom as pz

    conf = tmp_path / "settings.ini"
    monkeypatch.setattr(
        pz,
        "QSettings",
        lambda *a, **k: QSettings(str(conf), QSettings.Format.IniFormat),
    )

    a = tmp_path / "a.fa"
    a.write_text(
        ">s1\nACGTACGTACGTACGT\n>s2\nTTTTAAAACCCCGGGG\n", encoding="utf-8"
    )
    b = tmp_path / "b.fa"
    b.write_text(
        ">t1\nGGGGCCCCAAAATTTT\n>t2\nACGTACGTACGTACGT\n", encoding="utf-8"
    )

    host = QWidget()
    layout = QVBoxLayout(host)
    tabs = PreviewTabWidget()
    layout.addWidget(tabs, stretch=1)
    host.resize(1000, 600)
    host.show()
    QTest.qWait(50)

    assert tabs.open_file(a, replace=False) is True
    assert tabs.open_split_with_tabs(1, 2) is True
    QTest.qWait(400)
    assert tabs.is_split_active()

    # Simulate leaving the dock without turning Split off (preference stays on).
    # Directly close chrome without clearing the remembered MSA layout.
    tabs._close_inline_split(clear_preference=False)
    QTest.qWait(100)
    assert not tabs.is_split_active()

    assert tabs.open_file(b, replace=True) is True
    QTest.qWait(500)
    assert tabs.is_split_active()
    split = tabs._inline_split
    assert isinstance(split, PreviewSplitWindow)
    assert split.current_path() == b.resolve()
    assert split._left.current_fasta_tab() == 1
    assert split._right.current_fasta_tab() == 2

    # Turning Split off clears auto-restore for the next MSA.
    tabs._close_inline_split(clear_preference=True)
    QTest.qWait(100)
    c = tmp_path / "c.fa"
    c.write_text(
        ">u1\nAAAACCCCGGGGTTTT\n>u2\nACGTACGTACGTACGT\n", encoding="utf-8"
    )
    assert tabs.open_file(c, replace=True) is True
    QTest.qWait(200)
    assert not tabs.is_split_active()

    host.close()
    tabs.close()


def test_open_split_with_tabs_overview_and_dotplot(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Right-click path: Overview left, Dotplot right."""
    from PySide6.QtTest import QTest

    path = tmp_path / "split_tabs.fa"
    path.write_text(
        ">seq_a\nACGTACGTACGTACGTACGT\n"
        ">seq_b\nTTTTAAAACCCCGGGGAAAA\n",
        encoding="utf-8",
    )

    host = QWidget()
    layout = QVBoxLayout(host)
    tabs = PreviewTabWidget()
    layout.addWidget(tabs, stretch=1)
    host.resize(1000, 600)
    host.show()
    QTest.qWait(50)

    assert tabs.open_file(path, replace=False) is True
    page = tabs.currentWidget()
    assert isinstance(page, FilePreviewPage)
    page._open_annotation_tab()
    QTest.qWait(200)

    assert tabs.open_split_with_tabs(1, 2) is True
    split = tabs._inline_split
    assert isinstance(split, PreviewSplitWindow)
    QTest.qWait(500)

    left = split._left.preview_page()
    right = split._right.preview_page()
    assert isinstance(left, FilePreviewPage)
    assert isinstance(right, FilePreviewPage)
    assert left.is_annotation_preview_active()
    assert right.is_dotplot_preview_active()

    # Open-in-panel while already split: switch right pane to Text.
    split.set_pane_fasta_tab("right", 0)
    QTest.qWait(200)
    assert right.is_text_preview_active()

    tabs._close_inline_split()
    QTest.qWait(300)
    host.close()
    tabs.close()


def test_split_overview_dotplot_sequence_sync(qapp: QApplication, tmp_path: Path) -> None:
    """Selecting a sequence on Overview highlights the peer Dotplot, and vice versa."""
    from PySide6.QtTest import QTest

    path = tmp_path / "split_sync.fa"
    path.write_text(
        ">seq_a\nACGTACGTACGTACGTACGT\n"
        ">seq_b\nTTTTAAAACCCCGGGGAAAA\n"
        ">seq_c\nGGGGCCCCAAAATTTTGGGG\n",
        encoding="utf-8",
    )

    host = QWidget()
    layout = QVBoxLayout(host)
    tabs = PreviewTabWidget()
    layout.addWidget(tabs, stretch=1)
    host.resize(1000, 600)
    host.show()
    QTest.qWait(50)

    assert tabs.open_file(path, replace=False) is True
    assert tabs.open_split_with_tabs(1, 2) is True
    split = tabs._inline_split
    assert isinstance(split, PreviewSplitWindow)
    QTest.qWait(500)
    split._setup_view_link()

    left = split._left.preview_page()
    right = split._right.preview_page()
    assert isinstance(left, FilePreviewPage)
    assert isinstance(right, FilePreviewPage)
    assert left.is_annotation_preview_active()
    assert right.is_dotplot_preview_active()

    assert left.select_records(["seq_b"]) is True
    QTest.qWait(150)
    assert right.current_split_sequence_id() == "seq_b"

    right.focus_sequence_for_split("seq_c")
    right.splitSequenceActivated.emit("seq_c")
    QTest.qWait(150)
    assert left.current_split_sequence_id() == "seq_c"

    tabs._close_inline_split()
    QTest.qWait(300)
    host.close()
    tabs.close()


def test_split_dotplot_view_highlights_peer_overview(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Dotplot View applies the region on the Overview pane; Dotplot stays open."""
    from PySide6.QtTest import QTest

    path = tmp_path / "split_view.fa"
    path.write_text(
        ">seq_a\nACGTACGTACGTACGTACGTACGT\n"
        ">seq_b\nTTTTAAAACCCCGGGGAAAATTTT\n",
        encoding="utf-8",
    )

    host = QWidget()
    layout = QVBoxLayout(host)
    tabs = PreviewTabWidget()
    layout.addWidget(tabs, stretch=1)
    host.resize(1000, 600)
    host.show()
    QTest.qWait(50)

    assert tabs.open_file(path, replace=False) is True
    assert tabs.open_split_with_tabs(1, 2) is True
    split = tabs._inline_split
    assert isinstance(split, PreviewSplitWindow)
    QTest.qWait(500)
    split._setup_view_link()

    left = split._left.preview_page()
    right = split._right.preview_page()
    assert isinstance(left, FilePreviewPage)
    assert isinstance(right, FilePreviewPage)
    assert left.is_annotation_preview_active()
    assert right.is_dotplot_preview_active()

    assert right.focus_sequence_for_split("seq_b") is True
    QTest.qWait(100)

    # Mimic Dotplot Display → View for a selected x/y range.
    right._on_dotplot_view_in_annotation(3, 10)
    QTest.qWait(200)

    assert right.is_dotplot_preview_active()
    assert left.is_annotation_preview_active()
    overview = left._annotation_overview
    assert overview is not None
    region = overview._bars.region_selection()
    assert region == (3, 10)
    assert overview._bars.is_regional_selection()

    tabs._close_inline_split()
    QTest.qWait(300)
    host.close()
    tabs.close()


def test_fasta_tab_strip_context_menu_signals(qapp: QApplication) -> None:
    from repeatviewer.widgets.preview_controls import FastaTabStrip

    strip = FastaTabStrip()
    seen: list[tuple[int, str]] = []
    closed: list[bool] = []
    strip.openInPanelRequested.connect(lambda i, s: seen.append((i, s)))
    strip.closeSplitRequested.connect(lambda: closed.append(True))
    strip.set_split_active(True)
    assert strip._split_active is True
    assert strip._buttons[2].contextMenuPolicy().name == "CustomContextMenu"

    strip.openInPanelRequested.emit(2, "right")
    strip.openInPanelRequested.emit(1, "left")
    strip.closeSplitRequested.emit()
    assert seen == [(2, "right"), (1, "left")]
    assert closed == [True]
    strip.close()


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
    from PySide6.QtWidgets import QHBoxLayout, QLabel, QVBoxLayout

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
    from PySide6.QtWidgets import QLabel, QMainWindow

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


def test_overflow_hider_survives_deleted_overview_width_source(
    qapp: QApplication,
) -> None:
    """Switching Overview files must not crash on a deleted toolbar slot."""
    from PySide6.QtWidgets import QHBoxLayout

    from repeatviewer.widgets.preview_controls import AnnotationZoomBar

    host = QWidget()
    layout = QHBoxLayout(host)
    bar = AnnotationZoomBar()
    layout.addWidget(bar)
    slot = QWidget()
    hider = bar._overview_overflow
    hider.add_width_source(slot)
    assert slot in hider._width_sources

    slot.deleteLater()
    qapp.processEvents()
    # Deferred sync + prune must not raise ``already deleted``.
    hider.sync()
    hider.clear_width_sources()
    assert hider._width_sources == []
    host.close()


def test_open_pdf_leaves_split_without_clearing_fasta_preference(
    qapp: QApplication, tmp_path: Path, monkeypatch
) -> None:
    """Split is remembered per FASTA shape; PDF opens single-pane."""
    from repeatviewer.widgets import preview_tabs as pt

    saved: list[dict] = []
    monkeypatch.setattr(
        pt,
        "save_preview_split_preference",
        lambda kind, *, enabled, left, right: saved.append(
            {"kind": kind, "enabled": enabled, "left": left, "right": right}
        ),
    )
    monkeypatch.setattr(
        pt, "load_preview_split_preference", lambda kind: (True, 1, 2)
    )

    fasta = tmp_path / "msa.fa"
    fasta.write_text(">a\nACGTACGTAC\n>b\nTGCATGCATG\n", encoding="utf-8")
    pdf = tmp_path / "report.pdf"
    pdf.write_bytes(b"%PDF-1.4\n%\xe2\xe3\xcf\xd3\n1 0 obj<<>>endobj\ntrailer<<>>\n")

    host = QWidget()
    layout = QVBoxLayout(host)
    tabs = PreviewTabWidget()
    layout.addWidget(tabs, stretch=1)
    host.resize(900, 500)
    host.show()
    qapp.processEvents()

    assert tabs.open_file(fasta, replace=False) is True
    assert tabs.open_split_with_tabs(1, 2) is True
    assert tabs.is_split_active()
    assert any(row.get("enabled") for row in saved)
    saved.clear()

    assert tabs.open_file(pdf, replace=True) is True
    assert not tabs.is_split_active()
    # Leaving for PDF must not write enabled=False for the FASTA shape.
    assert not any(row.get("enabled") is False for row in saved)

    assert tabs.open_file(fasta, replace=True) is True
    assert tabs.is_split_active()  # remembered multi split comes back

    tabs._close_inline_split()
    from PySide6.QtTest import QTest

    QTest.qWait(200)
    host.close()
    tabs.close()


def test_sync_does_not_force_matching_displayed_annotations(
    qapp: QApplication, tmp_path: Path
) -> None:
    """Sync left/right keeps Displayed annotations independent per pane."""
    from repeatviewer.backends.sequence_features import SequenceFeature

    path = tmp_path / "sync_ann.fa"
    path.write_text(">a\nACGTACGTAC\n>b\nTGCATGCATG\n", encoding="utf-8")

    host = QWidget()
    layout = QVBoxLayout(host)
    tabs = PreviewTabWidget()
    layout.addWidget(tabs, stretch=1)
    host.resize(1000, 560)
    host.show()
    qapp.processEvents()

    assert tabs.open_file(path, replace=False) is True
    assert tabs.open_split_with_tabs(1, 1) is True
    split = tabs._inline_split
    assert isinstance(split, PreviewSplitWindow)
    from PySide6.QtTest import QTest

    QTest.qWait(300)
    split._setup_view_link()
    left = split._left.overview_panel()
    right = split._right.overview_panel()
    assert left is not None and right is not None
    assert not hasattr(split, "_on_linked_annotation_kinds")

    feats = [SequenceFeature(1, 4, "x", feature_type="LTR")]
    for panel in (left, right):
        panel._rail.set_available_annotation_kinds(feats)
        panel._rail.apply_visible_annotation_kinds(frozenset({"LTR"}))

    left._rail.sync_split_views.setChecked(True)
    qapp.processEvents()
    assert split._link_enabled

    # Turn annotations off on the left only — right must stay on.
    left._rail.annotation_kinds_all_check.setChecked(False)
    qapp.processEvents()
    assert left._rail.visible_annotation_kinds() == frozenset()
    assert right._rail.visible_annotation_kinds() == frozenset({"LTR"})

    tabs._close_inline_split()
    QTest.qWait(200)
    host.close()
    tabs.close()
