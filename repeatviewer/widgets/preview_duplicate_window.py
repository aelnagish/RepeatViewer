"""Standalone duplicate of the FASTA preview panel."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QHBoxLayout, QSizePolicy, QVBoxLayout, QWidget

from repeatviewer.ui.window_geometry import ensure_widget_on_screen
from repeatviewer.widgets.preview_controls import (
    AnnotationZoomBar,
    PreviewHeaderBar,
    PreviewZoomBar,
)
from repeatviewer.widgets.preview_tabs import FilePreviewPage, PreviewTabWidget


class PreviewDuplicateWindow(QWidget):
    """Floating preview window mirroring one open FASTA preview tab."""

    def __init__(
        self,
        source: FilePreviewPage,
        *,
        source_tabs: PreviewTabWidget,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self.setWindowTitle(f"Preview — {source.path.name}")
        self.resize(max(480, source.width()), max(360, source.height()))

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self._zoom_bar = PreviewZoomBar()
        self._annotation_zoom_bar = AnnotationZoomBar()
        self._zoom_bar.set_values(
            source_tabs._zoom_bar.current_zoom_percent(),
            source_tabs._zoom_bar.current_fit_mode(),
        )

        zoom_host = QWidget(self)
        zoom_host.setObjectName("PreviewFloatingZoomHost")
        zoom_host.setFixedHeight(18)
        zoom_host.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        zoom_layout = QHBoxLayout(zoom_host)
        zoom_layout.setContentsMargins(6, 0, 6, 0)
        zoom_layout.setSpacing(6)
        zoom_layout.addStretch(1)
        zoom_layout.addWidget(self._zoom_bar)
        zoom_layout.addWidget(self._annotation_zoom_bar)

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
        self._tabs.open_file(source.path, replace=True)
        page = self._tabs.currentWidget()
        if isinstance(page, FilePreviewPage):
            page.set_fasta_tab_index(source.current_fasta_tab_index())
            if source.is_split_duplicate_active():
                page.set_split_duplicate(True)

        layout.addWidget(zoom_host)
        layout.addWidget(self._header)
        layout.addWidget(self._tabs, stretch=1)

        ensure_widget_on_screen(self)
