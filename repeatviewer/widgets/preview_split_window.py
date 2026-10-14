"""Floating Preview window with two side-by-side panes of the same file."""

from __future__ import annotations

from PySide6.QtCore import QSize, Qt, QTimer, Slot
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

# Narrow enough for free splitter ratios; toolbar overflow hides controls.
_PANE_MIN_WIDTH = 120
_WINDOW_MIN_WIDTH = 320
_WINDOW_MIN_HEIGHT = 240


class _SplitPreviewPane(QWidget):
    """One independent Preview copy inside a split window."""

    def __init__(
        self,
        source: FilePreviewPage,
        *,
        source_tabs: PreviewTabWidget,
        parent: QWidget | None = None,
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
        self._tabs = PreviewTabWidget(
            zoom_bar=self._zoom_bar,
            header_bar=self._header,
            annotation_zoom_bar=self._annotation_zoom_bar,
        )
        self._tabs.setObjectName("PreviewSplitPaneTabs")
        self._tabs.setMinimumWidth(0)
        self._tabs.tabBar().hide()
        self._tabs.currentChanged.connect(lambda _i: self._sync_title_zoom_host())
        self._header.fasta_tab_changed.connect(lambda _i: self._sync_title_zoom_host())

        layout.addWidget(self._zoom_host)
        layout.addWidget(self._header)
        layout.addWidget(self._tabs, stretch=1)

        self._job_running = False
        self._job_work_steps = 1
        self._status_bar = QStatusBar(self)
        self._status_bar.setSizeGripEnabled(False)
        self._status_bar.setFixedHeight(26)
        self._task_progress = TaskProgressPanel()
        self._status_bar.addPermanentWidget(self._task_progress, 1)
        self._stop_btn = QPushButton("Stop")
        self._stop_btn.setToolTip("Stop the currently running job")
        self._stop_btn.setMaximumHeight(22)
        self._stop_btn.setEnabled(False)
        self._stop_btn.clicked.connect(self._stop_current_job)
        self._status_bar.addPermanentWidget(self._stop_btn)
        self._status_bar.hide()
        layout.addWidget(self._status_bar)

        self._source_overview_features = self._snapshot_source_features(source)
        fasta_tab = source.current_fasta_tab_index()
        self._tabs.open_file(source.path, replace=True)
        QTimer.singleShot(0, lambda idx=fasta_tab: self._finish_open(idx))
        QTimer.singleShot(50, lambda idx=fasta_tab: self._finish_open(idx))

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return QSize(_PANE_MIN_WIDTH, 200)

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

    def _finish_open(self, fasta_tab: int) -> None:
        page = self._tabs.currentWidget()
        if not isinstance(page, FilePreviewPage):
            self._sync_title_zoom_host()
            return
        if fasta_tab >= 0:
            if fasta_tab == 1 or page.is_annotation_preview_active():
                page._open_annotation_tab()
            if page._fasta_stack.currentIndex() != fasta_tab:
                page.set_fasta_tab_index(fasta_tab)
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


class PreviewSplitWindow(QWidget):
    """Floating Preview with two independent left/right panes of the same file.

    Parented to the main window as a Tool window (not global always-on-top).
    """

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
    ) -> None:
        super().__init__(parent, self._WINDOW_FLAGS)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self.setAttribute(Qt.WidgetAttribute.WA_MacAlwaysShowToolWindow, True)
        self.setWindowTitle(f"Split Preview — {source.path.name}")
        self.setMinimumSize(_WINDOW_MIN_WIDTH, _WINDOW_MIN_HEIGHT)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.setSizeConstraint(QVBoxLayout.SizeConstraint.SetNoConstraint)

        self._splitter = QSplitter(Qt.Orientation.Horizontal, self)
        self._splitter.setObjectName("PreviewSplitSplitter")
        # Free left/right ratio down to each pane's minimum width.
        self._splitter.setChildrenCollapsible(False)
        self._splitter.setHandleWidth(6)
        self._splitter.setOpaqueResize(True)

        self._left = _SplitPreviewPane(
            source, source_tabs=source_tabs, parent=self._splitter
        )
        self._right = _SplitPreviewPane(
            source, source_tabs=source_tabs, parent=self._splitter
        )

        self._splitter.addWidget(self._left)
        self._splitter.addWidget(self._right)
        self._splitter.setStretchFactor(0, 1)
        self._splitter.setStretchFactor(1, 1)
        # Ignore child size hints so the handle can move freely within mins.
        self._splitter.setSizes([1000, 1000])
        layout.addWidget(self._splitter, stretch=1)

        self._size_to_screen()
        # Equal halves after the window has a real width.
        QTimer.singleShot(0, self._apply_equal_split)
        ensure_widget_on_screen(self)

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return QSize(_WINDOW_MIN_WIDTH, _WINDOW_MIN_HEIGHT)

    def _apply_equal_split(self) -> None:
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
