"""Zoom toolbar and unified preview header for the browser panel."""

from __future__ import annotations

from PySide6.QtCore import QEvent, QObject, Qt, Signal
from PySide6.QtGui import QWheelEvent
from PySide6.QtWidgets import (
    QButtonGroup,
    QComboBox,
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSlider,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from teviewer.widgets.preview_zoom import (
    DEFAULT_PREVIEW_ZOOM_PERCENT,
    MAX_PREVIEW_ZOOM_PERCENT,
    MIN_PREVIEW_ZOOM_PERCENT,
    PREVIEW_WHEEL_ZOOM_STEP_PERCENT,
    PreviewFitMode,
    load_preview_fit_mode,
    load_preview_zoom_percent,
    save_preview_fit_mode,
    save_preview_zoom_percent,
)


ZOOM_BAR_HEIGHT = 24
ZOOM_CONTROL_HEIGHT = 20
ZOOM_SLIDER_HEIGHT = 18
ZOOM_SLIDER_MAX_WIDTH = 140


class PreviewZoomBar(QWidget):
    """Compact preview scale controls (embedded on the right of the header bar)."""

    zoom_changed = Signal(int)
    fit_mode_changed = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        preview_lbl = QLabel("Preview:")
        layout.addWidget(preview_lbl)

        self._fit_combo = QComboBox()
        self._fit_combo.setFixedHeight(ZOOM_CONTROL_HEIGHT)
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
        self._slider.setMaximumWidth(ZOOM_SLIDER_MAX_WIDTH)
        self._slider.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._slider.setRange(MIN_PREVIEW_ZOOM_PERCENT, MAX_PREVIEW_ZOOM_PERCENT)
        self._slider.valueChanged.connect(self._on_slider)
        layout.addWidget(self._slider)

        self._pct_label = QLabel()
        self._pct_label.setMinimumWidth(34)
        layout.addWidget(self._pct_label)

        self._reset_btn = QPushButton("Reset")
        self._reset_btn.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        self._reset_btn.setToolTip("Restore default preview zoom (100%)")
        self._reset_btn.clicked.connect(self._reset_zoom)
        layout.addWidget(self._reset_btn)

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

    def step_zoom(self, steps: int) -> None:
        if steps == 0:
            return
        value = self._slider.value() + steps * PREVIEW_WHEEL_ZOOM_STEP_PERCENT
        value = max(MIN_PREVIEW_ZOOM_PERCENT, min(MAX_PREVIEW_ZOOM_PERCENT, value))
        self._slider.setValue(value)


class AnnotationZoomBar(QWidget):
    """Sequence-map zoom controls for the preview header (Annotation tab)."""

    fit_mode_changed = Signal(str)
    reset_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        preview_lbl = QLabel("Preview:")
        layout.addWidget(preview_lbl)

        self._fit_combo = QComboBox()
        self._fit_combo.setFixedHeight(ZOOM_CONTROL_HEIGHT)
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
        self._slider.setMaximumWidth(ZOOM_SLIDER_MAX_WIDTH)
        self._slider.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._slider.setRange(0, 100)
        self._slider.setValue(0)
        self._slider.setToolTip(
            "Left: fit full sequence in window. Right: maximum zoom-in."
        )
        self._slider.valueChanged.connect(self._on_slider)
        layout.addWidget(self._slider)

        self._pct_label = QLabel()
        self._pct_label.setMinimumWidth(34)
        layout.addWidget(self._pct_label)

        self._reset_btn = QPushButton("Reset")
        self._reset_btn.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        self._reset_btn.setToolTip("Fit the full sequence in the window")
        self._reset_btn.clicked.connect(self.reset_requested.emit)
        layout.addWidget(self._reset_btn)

        self._block = False
        self._fit_combo.activated.connect(self._on_fit_activated)
        self._update_label(0)

    @property
    def slider(self) -> QSlider:
        return self._slider

    def sync_display(self, value: int) -> None:
        """Update slider and percent label without emitting map zoom."""
        self._block = True
        self._slider.blockSignals(True)
        self._slider.setValue(value)
        self._slider.blockSignals(False)
        self._update_label(value)
        self._block = False

    def _update_label(self, percent: int) -> None:
        self._pct_label.setText(f"{percent}%")

    def _on_slider(self, value: int) -> None:
        self._update_label(value)
        if self._block:
            return

    def _on_fit_activated(self, _index: int) -> None:
        if self._block:
            return
        index = self._fit_combo.currentIndex()
        raw = self._fit_combo.itemData(index, Qt.ItemDataRole.UserRole)
        if raw is None:
            raw = self._fit_combo.itemData(index)
        self.fit_mode_changed.emit(str(raw))


class PreviewWheelZoomFilter(QObject):
    """Scroll-wheel zoom for preview content (image/dotplot: wheel; text/html: Ctrl/Cmd+wheel)."""

    def __init__(self, tabs: QWidget, zoom_bar: PreviewZoomBar) -> None:
        super().__init__(tabs)
        self._tabs = tabs
        self._zoom_bar = zoom_bar
        self._watched: set[int] = set()

    def register_page(self, page: QWidget) -> None:
        self._watch_widget(page)

    def _watch_widget(self, widget: QWidget) -> None:
        key = id(widget)
        if key in self._watched:
            return
        self._watched.add(key)
        widget.installEventFilter(self)
        for child in widget.findChildren(QWidget):
            self._watch_widget(child)

    @staticmethod
    def _is_descendant_of(widget: QWidget, ancestor: QWidget | None) -> bool:
        if ancestor is None:
            return False
        current: QWidget | None = widget
        while current is not None:
            if current is ancestor:
                return True
            current = current.parentWidget()
        return False

    def _current_page(self):
        from teviewer.widgets.preview_tabs import FilePreviewPage, PreviewTabWidget

        if not isinstance(self._tabs, PreviewTabWidget):
            return None
        widget = self._tabs.currentWidget()
        return widget if isinstance(widget, FilePreviewPage) else None

    def _allows_plain_wheel(self, page, obj: QWidget) -> bool:
        if self._is_descendant_of(obj, page._image):
            return True
        dotplot = getattr(page, "_dotplot_panel", None)
        if dotplot is not None and page.is_dotplot_preview_active():
            image_preview = getattr(dotplot, "_image_preview", None)
            if image_preview is not None and self._is_descendant_of(obj, image_preview):
                return True
        return False

    def _can_zoom(self, page, obj: QWidget, wheel: QWheelEvent) -> bool:
        from PySide6.QtWidgets import QSpinBox

        if not self._zoom_bar.isVisible() or not self._zoom_bar._slider.isEnabled():
            return False
        if page.is_annotation_preview_active():
            return False
        if isinstance(obj, QSpinBox) or self._is_descendant_of(obj, page._sequence_info_bar):
            return False
        mod = wheel.modifiers() & (
            Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier
        )
        if mod:
            return True
        return self._allows_plain_wheel(page, obj)

    def eventFilter(self, watched, event) -> bool:  # noqa: N802 (Qt API)
        if event.type() != QEvent.Type.Wheel:
            return False
        wheel = event
        if not isinstance(wheel, QWheelEvent):
            return False
        page = self._current_page()
        if page is None:
            return False
        obj = watched if isinstance(watched, QWidget) else None
        if obj is None or not self._can_zoom(page, obj, wheel):
            return False
        steps = wheel_zoom_steps_from_event(wheel)
        if steps == 0:
            return False
        self._zoom_bar.step_zoom(steps)
        wheel.accept()
        return True


class FastaTabStrip(QWidget):
    """Text / Annotation / Dotplot switcher for the preview header."""

    tab_changed = Signal(int)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("FastaTabStrip")
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self._group = QButtonGroup(self)
        self._group.setExclusive(True)
        self._buttons: list[QPushButton] = []
        for index, label in enumerate(("Text", "Annotation", "Dotplot")):
            btn = QPushButton(label)
            btn.setObjectName("FastaPreviewTab")
            btn.setCheckable(True)
            btn.setFlat(True)
            btn.setCursor(Qt.CursorShape.PointingHandCursor)
            self._group.addButton(btn, index)
            layout.addWidget(btn)
            self._buttons.append(btn)
        self._buttons[0].setChecked(True)
        self._group.idClicked.connect(self.tab_changed.emit)

    def set_current_index(self, index: int) -> None:
        if 0 <= index < len(self._buttons):
            self._buttons[index].setChecked(True)

    def current_index(self) -> int:
        return self._group.checkedId()


class PreviewHeaderBar(QWidget):
    """Grey header: FASTA tabs left, zoom right, hairline divider below."""

    fasta_tab_changed = Signal(int)

    def __init__(
        self,
        preview_zoom_bar: PreviewZoomBar,
        parent: QWidget | None = None,
        *,
        annotation_zoom_bar: AnnotationZoomBar | None = None,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("PreviewHeaderBar")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        row = QHBoxLayout()
        row.setContentsMargins(8, 4, 8, 4)
        row.setSpacing(8)

        self._fasta_tabs = FastaTabStrip(self)
        self._fasta_tabs.hide()
        self._fasta_tabs.tab_changed.connect(self.fasta_tab_changed.emit)
        row.addWidget(self._fasta_tabs)

        row.addStretch(1)
        self._preview_zoom_bar = preview_zoom_bar
        self._annotation_zoom_bar = annotation_zoom_bar
        row.addWidget(preview_zoom_bar)
        if annotation_zoom_bar is not None:
            row.addWidget(annotation_zoom_bar)
            annotation_zoom_bar.hide()
        outer.addLayout(row)

        divider = QFrame(self)
        divider.setObjectName("PreviewHeaderDivider")
        divider.setFrameShape(QFrame.Shape.HLine)
        divider.setFixedHeight(1)
        outer.addWidget(divider)

    def set_header_zoom_mode(self, mode: str) -> None:
        """Show preview zoom, annotation zoom, or hide both."""
        show_preview = mode == "preview"
        show_annotation = mode == "annotation"
        self._preview_zoom_bar.setVisible(show_preview)
        if self._annotation_zoom_bar is not None:
            self._annotation_zoom_bar.setVisible(show_annotation)

    def set_fasta_tabs_visible(self, visible: bool) -> None:
        self._fasta_tabs.setVisible(visible)

    def sync_fasta_tab(self, index: int) -> None:
        self._fasta_tabs.set_current_index(index)

    def current_fasta_tab(self) -> int:
        return self._fasta_tabs.current_index()
