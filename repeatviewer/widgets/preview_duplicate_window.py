"""Standalone floating window that copies the current FASTA Preview."""

from __future__ import annotations

from PySide6.QtCore import Qt, QTimer
from PySide6.QtWidgets import QHBoxLayout, QSizePolicy, QVBoxLayout, QWidget

from repeatviewer.ui.window_geometry import ensure_widget_on_screen
from repeatviewer.widgets.preview_controls import (
    AnnotationZoomBar,
    PreviewHeaderBar,
    PreviewZoomBar,
)
from repeatviewer.widgets.preview_tabs import FilePreviewPage, PreviewTabWidget


class PreviewDuplicateWindow(QWidget):
    """Independent Preview window; file-list clicks in the main app do not change it.

    Uses a Tool window parented to the main window so it stays above RepeatViewer
    (without ``WindowStaysOnTopHint``, which would cover other apps).
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
        # Keep the tool window visible with the main app on macOS.
        self.setAttribute(Qt.WidgetAttribute.WA_MacAlwaysShowToolWindow, True)
        self.setWindowTitle(f"Preview — {source.path.name}")
        self.resize(max(640, source.width()), max(420, source.height()))

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self._zoom_bar = PreviewZoomBar()
        self._annotation_zoom_bar = AnnotationZoomBar()
        # Duplicate windows are sticky copies — no Actions / Duplicate again.
        self._zoom_bar.set_actions_visible(False)
        self._annotation_zoom_bar.set_actions_visible(False)
        try:
            self._zoom_bar.set_values(
                source_tabs._zoom_bar.current_zoom_percent(),
                source_tabs._zoom_bar.current_fit_mode(),
            )
        except Exception:
            pass

        # Title strip for Text font size only. Overview mounts AnnotationZoomBar
        # in its own toolbar (Select/Fit/…); keep this host collapsed then so it
        # does not leave a blank gray band above the tabs.
        self._zoom_host = QWidget(self)
        self._zoom_host.setObjectName("PreviewFloatingZoomHost")
        self._zoom_host.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
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
        self._tabs.setObjectName("PreviewDuplicateTabs")
        # Hide open-file tab strip — this window is a single sticky preview.
        self._tabs.tabBar().hide()
        self._tabs.currentChanged.connect(lambda _i: self._sync_title_zoom_host())
        self._header.fasta_tab_changed.connect(lambda _i: self._sync_title_zoom_host())

        layout.addWidget(self._zoom_host)
        layout.addWidget(self._header)
        layout.addWidget(self._tabs, stretch=1)

        # Keep a weak snapshot of source overview features so the duplicate can
        # mirror in-memory annotations if disk sidecars are empty.
        self._source_overview_features = self._snapshot_source_features(source)

        fasta_tab = source.current_fasta_tab_index()
        self._tabs.open_file(source.path, replace=True)
        # After open_file's own deferred tab open, bind chrome + load features.
        QTimer.singleShot(0, lambda idx=fasta_tab: self._finish_open(idx))
        # Second pass catches races where the first run bound before Overview
        # finished constructing its toolbar slot.
        QTimer.singleShot(50, lambda idx=fasta_tab: self._finish_open(idx))

        ensure_widget_on_screen(self)

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
        """Apply source tab, bind Overview chrome, and hide empty title strip."""
        page = self._tabs.currentWidget()
        if not isinstance(page, FilePreviewPage):
            self._sync_title_zoom_host()
            return
        if fasta_tab >= 0:
            # Force Overview/Annotation content even if the stack is already there
            # (set_fasta_tab_index no-ops when the index is unchanged).
            if fasta_tab == 1 or page.is_annotation_preview_active():
                page._open_annotation_tab()
            if page._fasta_stack.currentIndex() != fasta_tab:
                page.set_fasta_tab_index(fasta_tab)
        # Re-bind zoom bar after the page is visible so Select/Fit/Cons appear.
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
        # Keep header controls host collapsed — Overview owns the map toolbar.
        if hasattr(self._header, "annotation_controls_host"):
            self._header.annotation_controls_host().hide()
            if hasattr(self._header, "_update_header_height"):
                self._header._update_header_height()
        self._sync_title_zoom_host()

    def _apply_source_features_if_needed(self, page: FilePreviewPage) -> None:
        """If disk reload found no features, copy the source Overview's tracks."""
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
        """Show title zoom only for Text; Overview owns its toolbar below the tabs."""
        page = self._tabs.currentWidget()
        if not isinstance(page, FilePreviewPage):
            self._collapse_zoom_host()
            return
        if page.is_annotation_preview_active() or page.is_dotplot_preview_active():
            # Annotation/Overview chrome lives in the map toolbar, not here.
            self._collapse_zoom_host()
            return
        self._expand_zoom_host()
        if hasattr(self._zoom_bar, "set_chrome_mode"):
            self._zoom_bar.set_chrome_mode("text")
