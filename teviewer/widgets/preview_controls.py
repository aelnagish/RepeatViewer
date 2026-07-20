"""Zoom toolbar and unified preview header for the browser panel."""

from __future__ import annotations

from PySide6.QtCore import QEvent, QObject, Qt, Signal
from PySide6.QtGui import QAction, QWheelEvent
from PySide6.QtWidgets import (
    QButtonGroup,
    QFrame,
    QHBoxLayout,
    QMenu,
    QPushButton,
    QSizePolicy,
    QSlider,
    QSpinBox,
    QToolButton,
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
    save_preview_zoom_percent,
    wheel_zoom_steps_from_event,
)
from teviewer.widgets.spin_boxes import NoWheelSpinBox

ZOOM_BAR_HEIGHT = 18
ZOOM_CONTROL_HEIGHT = 18
ZOOM_SPIN_WIDTH = 68
# Kept for sequence_map_viewer imports (embedded map toolbar height).
ZOOM_SLIDER_HEIGHT = 22
ANNOTATION_ZOOM_SLIDER_WIDTH = 240


def apply_annotation_zoom_slider_style(slider: QSlider) -> None:
    """Accent groove + small circular thumb for the annotation map zoom slider."""
    from teviewer.ui.theme_mode import load_accent_color

    accent = load_accent_color()
    slider.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
    slider.setStyleSheet(
        f"""
        QSlider#annotationZoomSlider::groove:horizontal {{
            border: 1px solid palette(mid);
            background: palette(base);
            height: 4px;
            border-radius: 2px;
        }}
        QSlider#annotationZoomSlider::sub-page:horizontal {{
            background: {accent};
            border: 1px solid {accent};
            height: 4px;
            border-radius: 2px;
        }}
        QSlider#annotationZoomSlider::add-page:horizontal {{
            background: palette(mid);
            border: 1px solid palette(mid);
            height: 4px;
            border-radius: 2px;
        }}
        QSlider#annotationZoomSlider::handle:horizontal {{
            background: {accent};
            border: 1px solid palette(text);
            width: 10px;
            height: 10px;
            margin: -4px 0;
            border-radius: 5px;
        }}
        QSlider#annotationZoomSlider::handle:horizontal:hover {{
            background: {accent};
            border: 1px solid palette(highlight);
        }}
        """
    )


def _make_zoom_percent_spin(
    *,
    minimum: int,
    maximum: int,
    tooltip: str,
) -> NoWheelSpinBox:
    """Rounded percent spinbox (figure-2 style); value shows with a ``%`` suffix."""
    spin = NoWheelSpinBox()
    spin.setObjectName("previewZoomSpin")
    spin.setRange(minimum, maximum)
    spin.setSuffix("%")
    spin.setFrame(False)
    spin.setAttribute(Qt.WidgetAttribute.WA_MacShowFocusRect, False)
    spin.setFixedHeight(ZOOM_CONTROL_HEIGHT)
    spin.setFixedWidth(ZOOM_SPIN_WIDTH)
    spin.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
    spin.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    spin.setButtonSymbols(QSpinBox.ButtonSymbols.UpDownArrows)
    spin.setToolTip(tooltip)
    spin.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
    spin.setSingleStep(PREVIEW_WHEEL_ZOOM_STEP_PERCENT)
    spin.setKeyboardTracking(False)
    return spin


def _make_preview_tools_button(*, tooltip: str) -> QToolButton:
    """Bordered Tools menu control matching the preview zoom spinbox."""
    btn = QToolButton()
    btn.setObjectName("PreviewToolsButton")
    btn.setText("Tools")
    btn.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
    btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
    btn.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
    btn.setAttribute(Qt.WidgetAttribute.WA_MacShowFocusRect, False)
    btn.setFixedHeight(ZOOM_CONTROL_HEIGHT)
    btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
    btn.setToolTip(tooltip)
    btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
    return btn


class PreviewZoomBar(QWidget):
    """Compact preview scale controls (% spin / Tools) for the Preview title bar."""

    zoom_changed = Signal(int)
    fit_mode_changed = Signal(str)
    edit_text_requested = Signal()
    save_text_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("PreviewZoomBar")
        self.setFixedHeight(ZOOM_BAR_HEIGHT)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        self._spin = _make_zoom_percent_spin(
            minimum=MIN_PREVIEW_ZOOM_PERCENT,
            maximum=MAX_PREVIEW_ZOOM_PERCENT,
            tooltip="Preview zoom",
        )
        self._spin.valueChanged.connect(self._on_spin)
        layout.addWidget(self._spin)

        self._tools_btn = _make_preview_tools_button(
            tooltip="Preview tools",
        )
        self._tools_menu = QMenu(self._tools_btn)
        self._reset_scale_action = QAction("Reset scale", self._tools_menu)
        self._reset_scale_action.setToolTip("Restore default preview zoom (100%)")
        self._reset_scale_action.triggered.connect(self._reset_zoom)
        self._edit_text_action = QAction("Edit text", self._tools_menu)
        self._edit_text_action.setToolTip("Allow editing the text shown in Preview")
        self._edit_text_action.triggered.connect(self.edit_text_requested.emit)
        self._save_text_action = QAction("Save text", self._tools_menu)
        self._save_text_action.setToolTip("Save Preview text changes to disk")
        self._save_text_action.triggered.connect(self.save_text_requested.emit)
        self._tools_menu.addAction(self._reset_scale_action)
        self._tools_menu.addSeparator()
        self._tools_menu.addAction(self._edit_text_action)
        self._tools_menu.addAction(self._save_text_action)
        self._tools_btn.setMenu(self._tools_menu)
        layout.addWidget(self._tools_btn)

        self._block = False
        self._fit_mode = load_preview_fit_mode()
        self.set_values(load_preview_zoom_percent(), self._fit_mode)
        self.set_text_actions_enabled(can_edit=False, can_save=False)

    def set_text_actions_enabled(self, *, can_edit: bool, can_save: bool) -> None:
        self._edit_text_action.setEnabled(can_edit)
        self._save_text_action.setEnabled(can_save)

    def set_values(self, percent: int, fit_mode: PreviewFitMode) -> None:
        self._block = True
        self._spin.setValue(percent)
        self._fit_mode = fit_mode
        self._block = False

    def current_zoom_percent(self) -> int:
        return self._spin.value()

    def current_fit_mode(self) -> PreviewFitMode:
        return self._fit_mode

    def _on_spin(self, value: int) -> None:
        if self._block:
            return
        save_preview_zoom_percent(value)
        self.zoom_changed.emit(value)

    def _reset_zoom(self) -> None:
        self._spin.setValue(DEFAULT_PREVIEW_ZOOM_PERCENT)

    def step_zoom(self, steps: int) -> None:
        if steps == 0:
            return
        value = self._spin.value() + steps * PREVIEW_WHEEL_ZOOM_STEP_PERCENT
        value = max(MIN_PREVIEW_ZOOM_PERCENT, min(MAX_PREVIEW_ZOOM_PERCENT, value))
        self._spin.setValue(value)

    def refresh_slider_style(self) -> None:
        """Kept for theme-reload callers; spinbox uses the global stylesheet."""


class AnnotationZoomBar(QWidget):
    """Sequence-map zoom controls for the Preview title bar (Annotation tab)."""

    fit_mode_changed = Signal(str)
    reset_requested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("PreviewZoomBar")
        self.setFixedHeight(ZOOM_BAR_HEIGHT)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        self._spin = _make_zoom_percent_spin(
            minimum=0,
            maximum=100,
            tooltip="0%: fit full sequence in window. 100%: maximum zoom-in.",
        )
        self._spin.setSingleStep(5)
        self._spin.setValue(0)
        self._spin.valueChanged.connect(self._on_spin)
        layout.addWidget(self._spin)

        self._tools_btn = _make_preview_tools_button(
            tooltip="Annotation tools",
        )
        self._tools_menu = QMenu(self._tools_btn)
        self._reset_scale_action = QAction("Reset scale", self._tools_menu)
        self._reset_scale_action.setToolTip("Fit the full sequence in the window")
        self._reset_scale_action.triggered.connect(self.reset_requested.emit)
        self._tools_menu.addAction(self._reset_scale_action)
        self._tools_btn.setMenu(self._tools_menu)
        layout.addWidget(self._tools_btn)

        self._block = False

    @property
    def slider(self) -> QSpinBox:
        """Zoom percent control (name kept for sequence-map binding)."""
        return self._spin

    def sync_display(self, value: int) -> None:
        """Update percent spin without emitting map zoom."""
        self._block = True
        self._spin.blockSignals(True)
        self._spin.setValue(value)
        self._spin.blockSignals(False)
        self._block = False

    def refresh_slider_style(self) -> None:
        """Kept for theme-reload callers; spinbox uses the global stylesheet."""

    def _on_spin(self, _value: int) -> None:
        if self._block:
            return


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
        return False

    def _can_zoom(self, page, obj: QWidget, wheel: QWheelEvent) -> bool:
        if not self._zoom_bar.isVisible() or not self._zoom_bar._spin.isEnabled():
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


FASTA_TAB_STRIP_HEIGHT = 22
PREVIEW_HEADER_HEIGHT = 24


class FastaTabStrip(QWidget):
    """Text / Annotation switcher for the preview header."""

    tab_changed = Signal(int)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("FastaTabStrip")
        self.setFixedHeight(FASTA_TAB_STRIP_HEIGHT)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self._group = QButtonGroup(self)
        self._group.setExclusive(True)
        self._buttons: list[QPushButton] = []
        for index, label in enumerate(("Text", "Annotation", "Dotplot (Self)")):
            btn = QPushButton(label)
            btn.setObjectName("FastaPreviewTab")
            btn.setCheckable(True)
            btn.setFlat(True)
            btn.setFixedHeight(FASTA_TAB_STRIP_HEIGHT)
            btn.setSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)
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
    """FASTA Text/Annotation switcher. Zoom controls live on the Preview dock title."""

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
        self.setFixedHeight(PREVIEW_HEADER_HEIGHT)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        row = QHBoxLayout()
        row.setContentsMargins(6, 0, 6, 0)
        row.setSpacing(4)

        self._fasta_tabs = FastaTabStrip(self)
        self._fasta_tabs.hide()
        self._fasta_tabs.tab_changed.connect(self.fasta_tab_changed.emit)
        row.addWidget(self._fasta_tabs)
        row.addStretch(1)

        # Zoom bars are mounted on the Preview dock title bar (same row as
        # "Preview"); keep refs so set_header_zoom_mode() can show/hide them.
        self._preview_zoom_bar = preview_zoom_bar
        self._annotation_zoom_bar = annotation_zoom_bar

        outer.addLayout(row)

        divider = QFrame(self)
        divider.setObjectName("PreviewHeaderDivider")
        divider.setFrameShape(QFrame.Shape.HLine)
        divider.setFixedHeight(1)
        outer.addWidget(divider)
        self.hide()

    def _update_bar_visibility(self) -> None:
        # Header only shows for Text/Annotation tabs; zoom lives on the title bar.
        self.setVisible(not self._fasta_tabs.isHidden())

    def set_header_zoom_mode(self, mode: str) -> None:
        """Show preview zoom, annotation zoom, or hide both (on the dock title)."""
        show_preview = mode == "preview"
        show_annotation = mode == "annotation"
        self._preview_zoom_bar.setVisible(show_preview)
        if self._annotation_zoom_bar is not None:
            self._annotation_zoom_bar.setVisible(show_annotation)
        self._update_bar_visibility()

    def set_fasta_tabs_visible(self, visible: bool) -> None:
        self._fasta_tabs.setVisible(visible)
        self._update_bar_visibility()

    def sync_fasta_tab(self, index: int) -> None:
        self._fasta_tabs.set_current_index(index)

    def current_fasta_tab(self) -> int:
        return self._fasta_tabs.current_index()
