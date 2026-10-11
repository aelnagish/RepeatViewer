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
        self._annotation_zoom_bar.hide()
        try:
            self._zoom_bar.set_values(
                source_tabs._zoom_bar.current_zoom_percent(),
                source_tabs._zoom_bar.current_fit_mode(),
            )
        except Exception:
            pass

        # Match main Preview dock: Actions / zoom on a slim title strip.
        zoom_host = QWidget(self)
        zoom_host.setObjectName("PreviewFloatingZoomHost")
        zoom_host.setFixedHeight(22)
        zoom_host.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        zoom_layout = QHBoxLayout(zoom_host)
        zoom_layout.setContentsMargins(6, 0, 6, 0)
        zoom_layout.setSpacing(6)
        zoom_layout.addStretch(1)
        zoom_layout.addWidget(self._zoom_bar)

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

        layout.addWidget(zoom_host)
        layout.addWidget(self._header)
        layout.addWidget(self._tabs, stretch=1)

        fasta_tab = source.current_fasta_tab_index()
        self._tabs.open_file(source.path, replace=True)
        page = self._tabs.currentWidget()
        if isinstance(page, FilePreviewPage) and fasta_tab >= 0:
            QTimer.singleShot(0, lambda idx=fasta_tab: page.set_fasta_tab_index(idx))

        ensure_widget_on_screen(self)
