"""Zoom toolbar for plot/image previews."""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QComboBox,
    QDockWidget,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QSlider,
    QWidget,
)

from teviewer.widgets.preview_zoom import (
    DEFAULT_PREVIEW_ZOOM_PERCENT,
    MAX_PREVIEW_ZOOM_PERCENT,
    MIN_PREVIEW_ZOOM_PERCENT,
    PreviewFitMode,
    load_preview_fit_mode,
    load_preview_zoom_percent,
    save_preview_fit_mode,
    save_preview_zoom_percent,
)

ZOOM_BAR_HEIGHT = 24
ZOOM_CONTROL_HEIGHT = 20
ZOOM_SLIDER_HEIGHT = 16


def zoom_slider_stylesheet(object_name: str) -> str:
    """Shared horizontal zoom slider styling (preview + annotation)."""
    return (
        f"QSlider#{object_name}::groove:horizontal {{"
        "  height: 3px;"
        "  background: palette(mid);"
        "  border-radius: 1px;"
        "}"
        f"QSlider#{object_name}::sub-page:horizontal {{"
        "  background: #4a90d9;"
        "  border-radius: 1px;"
        "}"
        f"QSlider#{object_name}::handle:horizontal {{"
        "  width: 9px;"
        "  height: 9px;"
        "  margin: -4px 0;"
        "  background: palette(base);"
        "  border: 1px solid palette(mid);"
        "  border-radius: 4px;"
        "}"
        f"QSlider#{object_name}::handle:horizontal:hover {{"
        "  border-color: #4a90d9;"
        "}"
        f"QSlider#{object_name}:disabled::groove:horizontal {{"
        "  background: palette(midlight);"
        "}"
        f"QSlider#{object_name}:disabled::sub-page:horizontal {{"
        "  background: palette(mid);"
        "}"
        f"QSlider#{object_name}:disabled::handle:horizontal {{"
        "  background: palette(midlight);"
        "  border-color: palette(mid);"
        "}"
    )


def zoom_bar_container_stylesheet(object_name: str, *, embedded: bool = False) -> str:
    if embedded:
        return (
            f"QWidget#{object_name} {{"
            "  background: transparent;"
            "  border: none;"
            "}"
        )
    return (
        f"QWidget#{object_name} {{"
        "  background: palette(alternate-base);"
        "  border: 1px solid palette(mid);"
        "  border-radius: 4px;"
        "}"
    )


class PreviewDockTitleBar(QWidget):
    """Preview dock header: title on the left, zoom controls on the right."""

    def __init__(self, dock: QDockWidget, zoom_bar: "PreviewZoomBar") -> None:
        super().__init__(dock)
        self._dock = dock
        self.setObjectName("previewDockTitleBar")
        self.setStyleSheet(
            "QWidget#previewDockTitleBar {"
            "  background: palette(midlight);"
            "  border-bottom: 1px solid palette(mid);"
            "}"
            "QLabel#previewDockTitleLabel {"
            "  font-weight: 600;"
            "  padding: 5px 0px;"
            "}"
        )
        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 0, 2, 0)
        layout.setSpacing(6)

        self._title = QLabel(dock.windowTitle())
        self._title.setObjectName("previewDockTitleLabel")
        layout.addWidget(self._title, 0, Qt.AlignmentFlag.AlignVCenter)
        zoom_bar.setParent(self)
        zoom_bar.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        layout.addWidget(zoom_bar, 1, Qt.AlignmentFlag.AlignVCenter)
        self.setFixedHeight(max(ZOOM_BAR_HEIGHT, 22))

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802 (Qt API)
        if event.button() == Qt.MouseButton.LeftButton:
            self._dock.setFloating(not self._dock.isFloating())
        super().mouseDoubleClickEvent(event)


def attach_preview_zoom_to_dock_title(
    dock: QDockWidget,
    zoom_bar: "PreviewZoomBar",
) -> PreviewDockTitleBar:
    """Place preview zoom controls in the dock title row."""
    title_bar = PreviewDockTitleBar(dock, zoom_bar)
    dock.setTitleBarWidget(title_bar)
    return title_bar


class PreviewZoomBar(QWidget):
    """Adjust preview scale for HTML plots, images, and text files."""

    zoom_changed = Signal(int)
    fit_mode_changed = Signal(str)

    def __init__(self, parent: QWidget | None = None, *, show_heading: bool = True) -> None:
        super().__init__(parent)
        self._embedded_in_tabs = not show_heading
        if self._embedded_in_tabs:
            self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(4, 0, 4, 0)
        layout.setSpacing(4)

        if show_heading:
            preview_lbl = QLabel("Preview:")
            layout.addWidget(preview_lbl)

        self._fit_combo = QComboBox()
        self._fit_combo.setFixedHeight(20)
        for label, mode in (
            ("Fit window", PreviewFitMode.FIT_WINDOW),
            ("Fit width", PreviewFitMode.FIT_WIDTH),
            ("Original size", PreviewFitMode.ORIGINAL),
        ):
            index = self._fit_combo.count()
            self._fit_combo.addItem(label)
            self._fit_combo.setItemData(index, mode.value, Qt.ItemDataRole.UserRole)
        layout.addWidget(self._fit_combo)

        self._slider = QSlider(Qt.Orientation.Horizontal)
        self._slider.setObjectName("previewZoomSlider")
        self._slider.setFixedHeight(ZOOM_SLIDER_HEIGHT)
        self._slider.setRange(MIN_PREVIEW_ZOOM_PERCENT, MAX_PREVIEW_ZOOM_PERCENT)
        self._slider.valueChanged.connect(self._on_slider)
        self._slider.setStyleSheet(zoom_slider_stylesheet("previewZoomSlider"))
        self._slider.setMinimumWidth(160 if self._embedded_in_tabs else 72)
        layout.addWidget(self._slider, stretch=1)

        self._pct_label = QLabel()
        layout.addWidget(self._pct_label)

        self._reset_btn = QPushButton("Reset")
        self._reset_btn.setFixedHeight(20)
        self._reset_btn.setToolTip("Restore default preview zoom (100%)")
        self._reset_btn.clicked.connect(self._reset_zoom)
        layout.addWidget(self._reset_btn)

        self.setMinimumWidth(220)
        self.setFixedHeight(ZOOM_BAR_HEIGHT)
        self.setStyleSheet(
            zoom_bar_container_stylesheet("previewZoomBar", embedded=self._embedded_in_tabs)
        )
        self.setObjectName("previewZoomBar")
        self._block = False
        self._fit_combo.activated.connect(self._on_fit_activated)
        self.set_values(load_preview_zoom_percent(), load_preview_fit_mode())

    def set_values(self, percent: int, fit_mode: PreviewFitMode) -> None:
        self._block = True
        self._slider.setValue(percent)
        self._update_label(percent)
        index = self._fit_combo.findData(fit_mode.value, Qt.ItemDataRole.UserRole)
        if index >= 0:
            self._fit_combo.setCurrentIndex(index)
        self._block = False

    def current_zoom_percent(self) -> int:
        return self._slider.value()

    def current_fit_mode(self) -> PreviewFitMode:
        index = self._fit_combo.currentIndex()
        raw = self._fit_combo.itemData(index, Qt.ItemDataRole.UserRole)
        if raw is None:
            raw = self._fit_combo.itemData(index)
        try:
            return PreviewFitMode(str(raw))
        except ValueError:
            return PreviewFitMode.FIT_WINDOW

    def _update_label(self, percent: int) -> None:
        self._pct_label.setText(f"{percent}%")

    def _on_slider(self, value: int) -> None:
        self._update_label(value)
        if self._block:
            return
        save_preview_zoom_percent(value)
        self.zoom_changed.emit(value)

    def _on_fit_activated(self, _index: int) -> None:
        if self._block:
            return
        mode = self.current_fit_mode()
        save_preview_fit_mode(mode)
        self.fit_mode_changed.emit(mode.value)

    def _reset_zoom(self) -> None:
        self._slider.setValue(DEFAULT_PREVIEW_ZOOM_PERCENT)

    def set_controls_enabled(self, enabled: bool, *, reason: str = "") -> None:
        """Enable or grey out all zoom controls (bar stays visible)."""
        self.setEnabled(enabled)
        self.setToolTip("" if enabled else reason)
