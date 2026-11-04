"""Left/right split Preview panes (inline in the Preview dock, or floating)."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING

from PySide6.QtCore import QSize, Qt, QTimer, Signal, Slot
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import (
    QHBoxLayout,
    QPushButton,
    QSizePolicy,
    QSplitter,
    QStatusBar,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.ui.window_geometry import ensure_widget_on_screen
from repeatviewer.widgets.preview_controls import (
    AnnotationZoomBar,
    PreviewHeaderBar,
    PreviewZoomBar,
)
from repeatviewer.widgets.preview_tabs import FilePreviewPage, PreviewTabWidget
from repeatviewer.widgets.task_progress import TaskProgressPanel

if TYPE_CHECKING:
    from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel

# Narrow enough for free splitter ratios; toolbar overflow hides controls.
_PANE_MIN_WIDTH = 120
_WINDOW_MIN_WIDTH = 320
_WINDOW_MIN_HEIGHT = 240


class _SplitPreviewPane(QWidget):
    """One independent Preview copy inside a split view."""

    closeSplitRequested = Signal()

    def __init__(
        self,
        source: FilePreviewPage,
        *,
        source_tabs: PreviewTabWidget,
        parent: QWidget | None = None,
        initial_fasta_tab: int = -1,
    ) -> None:
        super().__init__(parent)
        self.setMinimumWidth(_PANE_MIN_WIDTH)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.setSizeConstraint(QVBoxLayout.SizeConstraint.SetNoConstraint)

        self._zoom_bar = PreviewZoomBar()
        self._annotation_zoom_bar = AnnotationZoomBar()
        # Actions stay on the Preview dock title bar (not duplicated per pane).
        self._zoom_bar.set_actions_visible(False)
        self._annotation_zoom_bar.set_actions_visible(False)
        try:
            self._zoom_bar.set_values(
                source_tabs._zoom_bar.current_zoom_percent(),
                source_tabs._zoom_bar.current_fit_mode(),
            )
        except Exception:
            pass

        self._zoom_host = QWidget(self)
        self._zoom_host.setObjectName("PreviewFloatingZoomHost")
        self._zoom_host.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        self._zoom_host.setMinimumWidth(0)
        zoom_layout = QHBoxLayout(self._zoom_host)
        zoom_layout.setContentsMargins(6, 0, 6, 0)
        zoom_layout.setSpacing(6)
        zoom_layout.addStretch(1)
        zoom_layout.addWidget(self._zoom_bar)
        self._collapse_zoom_host()

        self._header = PreviewHeaderBar(
            self._zoom_bar,
            self,
            annotation_zoom_bar=self._annotation_zoom_bar,
        )
        self._header.set_split_active(True)
        self._tabs = PreviewTabWidget(
            zoom_bar=self._zoom_bar,
            header_bar=self._header,
            annotation_zoom_bar=self._annotation_zoom_bar,
            allow_inline_split=False,
        )
        self._tabs.setObjectName("PreviewSplitPaneTabs")
        self._tabs.setMinimumWidth(0)
        self._tabs.tabBar().hide()
        self._tabs.currentChanged.connect(lambda _i: self._sync_title_zoom_host())
        self._header.fasta_tab_changed.connect(lambda _i: self._sync_title_zoom_host())
        self._header.closeSplitRequested.connect(self.closeSplitRequested.emit)

        layout.addWidget(self._zoom_host)
        layout.addWidget(self._header)
        layout.addWidget(self._tabs, stretch=1)

        self._job_running = False
        self._job_work_steps = 1
        self._status_bar = QStatusBar(self)
        self._status_bar.setSizeGripEnabled(False)
        self._status_bar.setFixedHeight(26)
        self._task_progress = TaskProgressPanel()
        self._status_bar.addPermanentWidget(self._task_progress)
        self._stop_btn = QPushButton("Stop")
        self._stop_btn.setToolTip("Stop the currently running job")
        self._stop_btn.setMaximumHeight(22)
        self._stop_btn.setEnabled(False)
        self._stop_btn.clicked.connect(self._stop_current_job)
        self._status_bar.addPermanentWidget(self._stop_btn)
        self._status_bar.hide()
        layout.addWidget(self._status_bar)

        self._source_overview_features = self._snapshot_source_features(source)
        fasta_tab = (
            int(initial_fasta_tab)
            if initial_fasta_tab >= 0
            else source.current_fasta_tab_index()
        )
        self._desired_fasta_tab = int(fasta_tab)
        self._tabs.open_file(source.path, replace=True)
        # File load is async; re-apply until the FASTA stack is ready (default
        # tab restore must not overwrite Overview↔Dotplot asymmetric splits).
        for delay in (0, 50, 120, 250, 500):
            QTimer.singleShot(delay, lambda: self._finish_open(self._desired_fasta_tab))

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return QSize(_PANE_MIN_WIDTH, 200)

    def preview_page(self) -> FilePreviewPage | None:
        page = self._tabs.currentWidget()
        return page if isinstance(page, FilePreviewPage) else None

    def open_path(self, path: Path, *, fasta_tab: int = -1) -> None:
        """Load ``path`` into this pane (keeps split alive across file clicks)."""
        resolved = Path(path).expanduser().resolve()
        if not resolved.is_file():
            return
        self._source_overview_features = None
        if fasta_tab >= 0:
            self._desired_fasta_tab = int(fasta_tab)
        self._tabs.open_file(resolved, replace=True)
        for delay in (0, 80, 200, 400):
            QTimer.singleShot(delay, lambda: self._finish_open(self._desired_fasta_tab))

    def set_fasta_tab(self, fasta_tab: int) -> None:
        """Switch this pane to a Text / Overview / Dotplot / Cytogenetics index."""
        self._desired_fasta_tab = int(fasta_tab)
        self._finish_open(self._desired_fasta_tab)

    def current_fasta_tab(self) -> int:
        page = self.preview_page()
        if page is None:
            return -1
        return page.current_fasta_tab_index()

    def focus_sequence(self, sequence_id: str) -> bool:
        page = self.preview_page()
        if page is None:
            return False
        return bool(page.focus_sequence_for_split(sequence_id))

    @staticmethod
    def _snapshot_source_features(source: FilePreviewPage) -> dict | None:
        overview = getattr(source, "_annotation_overview", None)
        if overview is None:
            return None
        features = getattr(overview, "_features", None)
        if not isinstance(features, dict) or not any(features.values()):
            return None
        try:
            from repeatviewer.widgets.multi_sequence_overview import (
                _clone_features_by_id,
            )

            return _clone_features_by_id(features)
        except Exception:
            return {rid: list(feats or []) for rid, feats in features.items()}

    def _collapse_zoom_host(self) -> None:
        self._zoom_host.hide()
        self._zoom_host.setFixedHeight(0)
        self._zoom_bar.hide()

    def _expand_zoom_host(self) -> None:
        self._zoom_host.setFixedHeight(22)
        self._zoom_host.show()
        self._zoom_bar.show()

    @staticmethod
    def _apply_fasta_tab(page: FilePreviewPage, fasta_tab: int) -> None:
        if fasta_tab < 0:
            return
        from repeatviewer.widgets.preview_zoom import DefaultFastaPreviewTab

        if fasta_tab == 1:
            page._open_annotation_tab()
        elif fasta_tab == 2:
            page._activate_fasta_tab(DefaultFastaPreviewTab.DOTPLOT)
        elif fasta_tab == 3:
            page.open_cyto_tab()
        else:
            page.set_fasta_tab_index(0)
        if page._fasta_stack.currentIndex() != fasta_tab:
            page.set_fasta_tab_index(fasta_tab)

    def _finish_open(self, fasta_tab: int) -> None:
        page = self._tabs.currentWidget()
        if not isinstance(page, FilePreviewPage):
            self._sync_title_zoom_host()
            return
        want = int(fasta_tab)
        # Always restore Text/Overview/Dotplot — including multi-seq MSA, where
        # file load defaults to Overview and would otherwise wipe a Dotplot pane.
        if want >= 0:
            self._apply_fasta_tab(page, want)
        if page.is_multi_fasta_preview() and page.is_annotation_preview_active():
            overview = getattr(page, "_annotation_overview", None)
            if overview is not None:
                overview.bind_annotation_zoom_bar(self._annotation_zoom_bar)
            page._ensure_annotation_data_loaded(
                progress_detail="Loading annotations…"
            )
            self._apply_source_features_if_needed(page)
        elif (
            page.is_single_record_fasta_preview()
            and page.is_annotation_preview_active()
        ):
            panel = getattr(page, "_annotation_panel", None)
            if panel is not None:
                panel.bind_annotation_zoom_bar(self._annotation_zoom_bar)
            page._ensure_annotation_data_loaded(
                progress_detail="Loading annotations…"
            )
        self._tabs._sync_zoom_bar_enabled()
        if hasattr(self._header, "annotation_controls_host"):
            self._header.annotation_controls_host().hide()
            if hasattr(self._header, "_update_header_height"):
                self._header._update_header_height()
        self._sync_title_zoom_host()

    def _apply_source_features_if_needed(self, page: FilePreviewPage) -> None:
        snapshot = self._source_overview_features
        if not snapshot:
            return
        overview = getattr(page, "_annotation_overview", None)
        if overview is None:
            return
        current = getattr(overview, "_features", None) or {}
        if any(current.values()):
            return
        overview._features = {
            rid: list(feats or []) for rid, feats in snapshot.items()
        }
        bars = getattr(overview, "_bars", None)
        if bars is not None and hasattr(bars, "set_features_in_place"):
            bars.set_features_in_place(snapshot)

    def _sync_title_zoom_host(self) -> None:
        page = self._tabs.currentWidget()
        if not isinstance(page, FilePreviewPage):
            self._collapse_zoom_host()
            return
        if page.is_annotation_preview_active() or page.is_dotplot_preview_active():
            self._collapse_zoom_host()
            return
        self._expand_zoom_host()
        if hasattr(self._zoom_bar, "set_chrome_mode"):
            self._zoom_bar.set_chrome_mode("text")

    @staticmethod
    def _running_job_title(job_title: str) -> str:
        title = str(job_title or "job").strip() or "job"
        if title.lower().startswith("running "):
            return title if title.endswith("…") else f"{title}…"
        return f"Running {title}…"

    @staticmethod
    def _parse_progress_filename(label: str) -> str:
        if ": " in label:
            part = label.split(": ", 1)[1]
        elif label.lower().startswith("aligning "):
            part = label[9:]
        else:
            return ""
        part = part.removesuffix("…").strip()
        if not part or part.lower() in {"finished", "alignment finished."}:
            return ""
        return part

    def _show_job_progress(self, job_title: str, total_steps: int) -> None:
        self._job_running = True
        self._job_work_steps = max(1, int(total_steps))
        self._status_bar.show()
        if self._task_progress is not None:
            self._task_progress.start(
                self._running_job_title(job_title),
                self._job_work_steps + 1,
            )
        if self._stop_btn is not None:
            self._stop_btn.setEnabled(True)

    @Slot(int, str)
    def _on_job_progress(self, value: int, label: str) -> None:
        if self._task_progress is None:
            return
        work = getattr(self, "_job_work_steps", self._task_progress._total)
        capped = min(max(0, int(value)), max(0, int(work)))
        detail = self._parse_progress_filename(label)
        if not detail and label:
            detail = label.split(": ", 1)[-1].strip() if ": " in label else label
            detail = detail.removesuffix("…").strip()
        self._task_progress.update(capped, detail)

    def _hide_job_progress(self) -> None:
        self._job_running = False
        if self._task_progress is not None:
            self._task_progress.reset()
        if self._stop_btn is not None:
            self._stop_btn.setEnabled(False)
        self._status_bar.hide()

    def _stop_current_job(self) -> None:
        from repeatviewer.widgets.multi_sequence_overview import (
            MultiSequenceOverviewPanel,
        )
        from repeatviewer.widgets.sequence_map_viewer import SequenceMapPanel

        for panel in self.findChildren(SequenceMapPanel):
            if (
                panel.cancel_auto_annotation()
                or panel.cancel_similarity_annotation()
                or panel.cancel_transfer_annotation()
            ):
                if self._stop_btn is not None:
                    self._stop_btn.setEnabled(False)
                return

        for panel in self.findChildren(MultiSequenceOverviewPanel):
            if panel.cancel_auto_annotation() or panel.cancel_similarity_annotation():
                if self._stop_btn is not None:
                    self._stop_btn.setEnabled(False)
                return

        if self._stop_btn is not None:
            self._stop_btn.setEnabled(False)

    def overview_panel(self) -> MultiSequenceOverviewPanel | None:
        """Multi-FASTA Overview panel for this pane, if loaded."""
        page = self.preview_page()
        if page is None:
            return None
        overview = getattr(page, "_annotation_overview", None)
        if overview is None:
            return None
        if not page.is_annotation_preview_active():
            return None
        return overview


class PreviewSplitWindow(QWidget):
    """Two side-by-side Preview panes of the same file.

    Default: embedded in the Preview dock (replaces the single-pane tab view).
    Optional ``floating=True`` keeps the legacy Tool window.
    ``Sync left and right view`` links vertical position, Displayed annotations,
    and selection. Overview↔Dotplot always mirrors the active sequence
    selection both ways while the split is open.
    """

    closeSplitRequested = Signal()

    _WINDOW_FLAGS = (
        Qt.WindowType.Tool
        | Qt.WindowType.CustomizeWindowHint
        | Qt.WindowType.WindowTitleHint
        | Qt.WindowType.WindowSystemMenuHint
        | Qt.WindowType.WindowMinMaxButtonsHint
        | Qt.WindowType.WindowCloseButtonHint
    )

    def __init__(
        self,
        source: FilePreviewPage,
        *,
        source_tabs: PreviewTabWidget,
        parent: QWidget | None = None,
        floating: bool = False,
        left_tab: int = -1,
        right_tab: int = -1,
    ) -> None:
        if floating:
            super().__init__(parent, self._WINDOW_FLAGS)
            self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
            self.setAttribute(Qt.WidgetAttribute.WA_MacAlwaysShowToolWindow, True)
            self.setWindowTitle(f"Split Preview — {source.path.name}")
            self.setMinimumSize(_WINDOW_MIN_WIDTH, _WINDOW_MIN_HEIGHT)
        else:
            super().__init__(parent)
            self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
            self.setMinimumSize(_PANE_MIN_WIDTH * 2, 200)
            self.setSizePolicy(
                QSizePolicy.Policy.Expanding,
                QSizePolicy.Policy.Expanding,
            )

        self._floating = bool(floating)
        self._source_tabs = source_tabs
        self._current_path = Path(source.path).expanduser().resolve()
        default_tab = source.current_fasta_tab_index()
        left_idx = int(left_tab) if left_tab >= 0 else default_tab
        right_idx = int(right_tab) if right_tab >= 0 else default_tab

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.setSizeConstraint(QVBoxLayout.SizeConstraint.SetNoConstraint)

        self._splitter = QSplitter(Qt.Orientation.Horizontal, self)
        self._splitter.setObjectName("PreviewSplitSplitter")
        self._splitter.setChildrenCollapsible(False)
        self._splitter.setHandleWidth(10)
        self._splitter.setOpaqueResize(True)
        # Clear vertical band so left/right panes read as separate windows.
        self._splitter.setStyleSheet(
            """
            QSplitter#PreviewSplitSplitter::handle:horizontal {
                background: #c8cdd4;
                margin: 0px 2px;
                border-left: 1px solid #aeb5bf;
                border-right: 1px solid #aeb5bf;
                width: 10px;
            }
            QSplitter#PreviewSplitSplitter::handle:horizontal:hover {
                background: #9aa3af;
                border-left: 1px solid #7f8894;
                border-right: 1px solid #7f8894;
            }
            """
        )

        self._left = _SplitPreviewPane(
            source,
            source_tabs=source_tabs,
            parent=self._splitter,
            initial_fasta_tab=left_idx,
        )
        self._right = _SplitPreviewPane(
            source,
            source_tabs=source_tabs,
            parent=self._splitter,
            initial_fasta_tab=right_idx,
        )

        self._splitter.addWidget(self._left)
        self._splitter.addWidget(self._right)
        self._splitter.setStretchFactor(0, 1)
        self._splitter.setStretchFactor(1, 1)
        self._splitter.setSizes([1000, 1000])
        layout.addWidget(self._splitter, stretch=1)

        self._link_enabled = False
        self._link_applying = False
        self._link_wired = False
        self._seq_link_wired = False
        self._overview_link_wired = False
        self._closed = False

        self._left.closeSplitRequested.connect(self.closeSplitRequested.emit)
        self._right.closeSplitRequested.connect(self.closeSplitRequested.emit)
        self._left._header.openInPanelRequested.connect(
            lambda tab, side: self._on_open_in_panel(tab, side)
        )
        self._right._header.openInPanelRequested.connect(
            lambda tab, side: self._on_open_in_panel(tab, side)
        )

        if floating:
            self._size_to_screen()
            QTimer.singleShot(0, self._apply_equal_split)
            ensure_widget_on_screen(self)
        else:
            QTimer.singleShot(0, self._apply_equal_split)

        QTimer.singleShot(0, self._setup_view_link)
        QTimer.singleShot(80, self._setup_view_link)
        QTimer.singleShot(250, self._setup_view_link)

    def mark_closed(self) -> None:
        """Stop deferred sync/open work before the widget is destroyed."""
        self._closed = True

    def open_path(self, path: Path, *, fasta_tab: int = -1) -> None:
        """Load the same file into both panes (file-list click while split)."""
        if self._closed:
            return
        resolved = Path(path).expanduser().resolve()
        if not resolved.is_file():
            return
        self._current_path = resolved
        if self._floating:
            self.setWindowTitle(f"Split Preview — {resolved.name}")
        # Prefer each pane's desired tab (survives async reload) over the live
        # stack index, which briefly falls back to Overview during file load.
        left_tab = int(getattr(self._left, "_desired_fasta_tab", -1))
        right_tab = int(getattr(self._right, "_desired_fasta_tab", -1))
        if left_tab < 0:
            left_tab = self._left.current_fasta_tab()
        if right_tab < 0:
            right_tab = self._right.current_fasta_tab()
        if fasta_tab >= 0:
            # Keep each pane's current tab when possible; seed only if unknown.
            if left_tab < 0:
                left_tab = fasta_tab
            if right_tab < 0:
                right_tab = fasta_tab
        self._left.open_path(resolved, fasta_tab=left_tab)
        self._right.open_path(resolved, fasta_tab=right_tab)
        self._link_wired = False
        self._seq_link_wired = False
        self._overview_link_wired = False
        QTimer.singleShot(0, self._setup_view_link)
        QTimer.singleShot(80, self._setup_view_link)
        QTimer.singleShot(250, self._setup_view_link)

    def current_path(self) -> Path:
        return self._current_path

    def set_pane_fasta_tab(self, side: str, fasta_tab: int) -> None:
        """Switch left or right pane to ``fasta_tab`` (0–3)."""
        if self._closed:
            return
        pane = self._left if str(side) == "left" else self._right
        pane.set_fasta_tab(int(fasta_tab))
        QTimer.singleShot(0, self._setup_view_link)
        QTimer.singleShot(80, self._setup_view_link)

    def _on_open_in_panel(self, tab_index: int, side: str) -> None:
        self.set_pane_fasta_tab(str(side or "right"), int(tab_index))

    def _setup_view_link(self) -> None:
        """Wire Overview↔Overview Sync when both exist; always wire sequence sync."""
        if self._closed:
            return
        left_page = self._left.preview_page()
        right_page = self._right.preview_page()
        if left_page is None or right_page is None:
            return

        if not self._seq_link_wired:
            self._seq_link_wired = True
            left_page.splitSequenceActivated.connect(
                lambda sid, src=self._left: self._on_pane_sequence(src, sid)
            )
            right_page.splitSequenceActivated.connect(
                lambda sid, src=self._right: self._on_pane_sequence(src, sid)
            )

        left = self._left.overview_panel()
        right = self._right.overview_panel()
        for overview in (left, right):
            if overview is None:
                continue
            rail = getattr(overview, "_rail", None)
            if rail is not None and hasattr(rail, "set_split_sync_available"):
                rail.set_split_sync_available(True)
            if self._link_enabled and rail is not None:
                check = getattr(rail, "sync_split_views", None)
                if check is not None and not check.isChecked():
                    check.blockSignals(True)
                    check.setChecked(True)
                    check.blockSignals(False)

        if left is None or right is None:
            return
        if self._overview_link_wired:
            return
        self._overview_link_wired = True
        self._link_wired = True

        for overview in (left, right):
            rail = overview._rail
            rail.syncSplitViewsChanged.connect(self._on_sync_toggled)
            rail.annotationTrackKindsChanged.connect(
                lambda keys, src=overview: self._on_linked_annotation_kinds(src, keys)
            )
            overview._bars.selectionChanged.connect(
                lambda src=overview: self._on_linked_selection(src)
            )
            overview._bars_scroll.verticalScrollBar().valueChanged.connect(
                lambda _v, src=overview: self._on_linked_vertical_scroll(src)
            )

    def _peer_pane(self, source: _SplitPreviewPane) -> _SplitPreviewPane | None:
        if source is self._left:
            return self._right
        if source is self._right:
            return self._left
        return None

    def _on_pane_sequence(self, source: _SplitPreviewPane, sequence_id: str) -> None:
        """Mirror sequence selection between Overview and Dotplot panes."""
        if self._closed or self._link_applying:
            return
        sid = str(sequence_id or "").strip()
        if not sid:
            return
        peer = self._peer_pane(source)
        if peer is None:
            return
        src_page = source.preview_page()
        peer_page = peer.preview_page()
        if src_page is None or peer_page is None:
            return
        src_dot = src_page.is_dotplot_preview_active()
        peer_dot = peer_page.is_dotplot_preview_active()
        # Overview↔Overview selection stays behind the Sync checkbox.
        if not src_dot and not peer_dot:
            return
        src_seq = src_page.is_annotation_preview_active() or src_dot
        peer_seq = peer_page.is_annotation_preview_active() or peer_dot
        if not src_seq or not peer_seq:
            return
        if src_dot == peer_dot:
            return
        current = peer_page.current_split_sequence_id()
        if current == sid and peer_dot:
            self._link_applying = True
            try:
                peer.focus_sequence(sid)
            finally:
                self._link_applying = False
            return
        if current == sid:
            return
        self._link_applying = True
        try:
            peer.focus_sequence(sid)
        finally:
            self._link_applying = False

    def try_apply_dotplot_view(
        self, source_page: FilePreviewPage, start: int, end: int
    ) -> bool:
        """Dotplot View → highlight peer Overview; keep source on Dotplot.

        Returns True when the peer Annotation/Overview pane accepted the region.
        """
        if self._closed:
            return False
        if not source_page.is_dotplot_preview_active():
            return False
        source_pane: _SplitPreviewPane | None = None
        if self._left.preview_page() is source_page:
            source_pane = self._left
        elif self._right.preview_page() is source_page:
            source_pane = self._right
        if source_pane is None:
            return False
        peer = self._peer_pane(source_pane)
        if peer is None:
            return False
        peer_page = peer.preview_page()
        if peer_page is None or not peer_page.is_annotation_preview_active():
            return False
        sid = source_page.current_split_sequence_id()
        self._link_applying = True
        try:
            return bool(
                peer_page.apply_region_highlight(
                    int(start), int(end), sequence_id=sid
                )
            )
        finally:
            self._link_applying = False

    def _peer_overview(
        self, source: MultiSequenceOverviewPanel
    ) -> MultiSequenceOverviewPanel | None:
        left = self._left.overview_panel()
        right = self._right.overview_panel()
        if source is left:
            return right
        if source is right:
            return left
        return None

    def _on_sync_toggled(self, enabled: bool) -> None:
        self._link_enabled = bool(enabled)
        left = self._left.overview_panel()
        right = self._right.overview_panel()
        for overview in (left, right):
            if overview is None:
                continue
            check = overview._rail.sync_split_views
            check.blockSignals(True)
            check.setChecked(self._link_enabled)
            check.blockSignals(False)
        if self._link_enabled and left is not None and right is not None:
            self._mirror_all(left, right)

    def _mirror_all(
        self,
        source: MultiSequenceOverviewPanel,
        target: MultiSequenceOverviewPanel,
    ) -> None:
        if self._link_applying:
            return
        self._link_applying = True
        try:
            anchor = source.linked_viewport_row_anchor()
            if anchor is not None:
                target.apply_linked_viewport_row_anchor(*anchor)
            target.apply_linked_annotation_kinds(
                source._rail.visible_annotation_kinds()
            )
            rows = source._bars.selected_rows()
            target.apply_linked_selected_rows(rows)
            self._align_peer_to_selected_row(source, target, rows)
        finally:
            self._link_applying = False

    @staticmethod
    def _align_peer_to_selected_row(
        source: MultiSequenceOverviewPanel,
        target: MultiSequenceOverviewPanel,
        rows: list[int],
    ) -> None:
        """Keep the selected sequence on the same vertical screen level."""
        if not rows:
            return
        row = int(rows[0])
        frac = source.linked_row_viewport_fraction(row)
        if frac is None:
            return
        target.apply_linked_row_at_viewport_fraction(row, frac)

    def _on_linked_vertical_scroll(
        self, source: MultiSequenceOverviewPanel
    ) -> None:
        if not self._link_enabled or self._link_applying:
            return
        target = self._peer_overview(source)
        if target is None:
            return
        anchor = source.linked_viewport_row_anchor()
        if anchor is None:
            return
        self._link_applying = True
        try:
            target.apply_linked_viewport_row_anchor(*anchor)
        finally:
            self._link_applying = False

    def _on_linked_selection(self, source: MultiSequenceOverviewPanel) -> None:
        if not self._link_enabled or self._link_applying:
            return
        target = self._peer_overview(source)
        if target is None:
            return
        rows = source._bars.selected_rows()
        self._link_applying = True
        try:
            target.apply_linked_selected_rows(rows)
            self._align_peer_to_selected_row(source, target, rows)
        finally:
            self._link_applying = False

    def _on_linked_annotation_kinds(
        self, source: MultiSequenceOverviewPanel, keys: object
    ) -> None:
        if not self._link_enabled or self._link_applying:
            return
        target = self._peer_overview(source)
        if target is None:
            return
        if keys is None:
            visible = source._rail.visible_annotation_kinds()
        else:
            visible = frozenset(str(k) for k in keys)
        self._link_applying = True
        try:
            target.apply_linked_annotation_kinds(visible)
        finally:
            self._link_applying = False

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        if self._floating:
            return QSize(_WINDOW_MIN_WIDTH, _WINDOW_MIN_HEIGHT)
        return QSize(_PANE_MIN_WIDTH * 2, 200)

    def _apply_equal_split(self) -> None:
        if self._closed:
            return
        total = max(2, self._splitter.width())
        half = total // 2
        self._splitter.setSizes([half, total - half])

    def _size_to_screen(self) -> None:
        """Open nearly full-screen so left/right panes are usable."""
        host = self.parentWidget()
        screen = None
        if host is not None:
            screen = host.screen()
        if screen is None:
            screen = QGuiApplication.primaryScreen()
        if screen is None:
            self.resize(1200, 720)
            return
        available = screen.availableGeometry()
        width = max(800, int(available.width() * 0.92))
        height = max(480, int(available.height() * 0.88))
        x = available.x() + max(0, (available.width() - width) // 2)
        y = available.y() + max(0, (available.height() - height) // 2)
        self.setGeometry(x, y, width, height)
