"""Zoom toolbar and unified preview header for the browser panel."""

from __future__ import annotations

from PySide6.QtCore import QEvent, QObject, QPointF, Qt, Signal
from PySide6.QtGui import (
    QAction,
    QActionGroup,
    QFont,
    QFontMetrics,
    QGuiApplication,
    QMouseEvent,
    QWheelEvent,
)
from PySide6.QtWidgets import (
    QButtonGroup,
    QFrame,
    QHBoxLayout,
    QLineEdit,
    QMenu,
    QPushButton,
    QSizePolicy,
    QSlider,
    QSpinBox,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.preview_zoom import (
    DEFAULT_PREVIEW_FONT_SIZE,
    DEFAULT_PREVIEW_ZOOM_PERCENT,
    MAX_PREVIEW_FONT_SIZE,
    MAX_PREVIEW_ZOOM_PERCENT,
    MIN_PREVIEW_FONT_SIZE,
    MIN_PREVIEW_ZOOM_PERCENT,
    PREVIEW_WHEEL_DELTA_PER_STEP,
    PREVIEW_WHEEL_FONT_STEP_PT,
    PREVIEW_WHEEL_ZOOM_STEP_PERCENT,
    DefaultFastaPreviewTab,
    PreviewFitMode,
    PreviewZoomKind,
    load_default_fasta_preview_tab,
    load_preview_fit_mode,
    load_preview_font_size,
    load_preview_zoom_percent,
    save_default_fasta_preview_tab,
    save_preview_font_size,
    save_preview_zoom_percent,
    wheel_zoom_delta_from_event,
)
from repeatviewer.widgets.spin_boxes import NoWheelSpinBox

_FASTA_DEFAULT_TAB_CHOICES: tuple[tuple[str, DefaultFastaPreviewTab], ...] = (
    ("Text", DefaultFastaPreviewTab.TEXT),
    ("Annotation", DefaultFastaPreviewTab.ANNOTATION),
    ("Dotplot (Self)", DefaultFastaPreviewTab.DOTPLOT),
)


def _add_default_fasta_tab_menu(parent_menu: QMenu) -> QAction:
    """Add Settings-equivalent default FASTA sub-tab chooser to an Actions menu."""
    submenu = parent_menu.addMenu("Default preview tab")
    submenu.setToolTip(
        "Sub-tab shown when you single-click a single-sequence FASTA file."
    )
    group = QActionGroup(submenu)
    group.setExclusive(True)
    for label, tab in _FASTA_DEFAULT_TAB_CHOICES:
        action = QAction(label, submenu)
        action.setCheckable(True)
        action.setData(tab.value)
        group.addAction(action)
        submenu.addAction(action)

    def _sync_checks() -> None:
        current = load_default_fasta_preview_tab().value
        for action in group.actions():
            action.setChecked(action.data() == current)

    def _on_chosen(action: QAction) -> None:
        try:
            tab = DefaultFastaPreviewTab(str(action.data()))
        except ValueError:
            return
        save_default_fasta_preview_tab(tab)

    submenu.aboutToShow.connect(_sync_checks)
    group.triggered.connect(_on_chosen)
    _sync_checks()
    return submenu.menuAction()

ZOOM_BAR_HEIGHT = 18
ZOOM_CONTROL_HEIGHT = 18
ZOOM_SPIN_WIDTH = 68
# Kept for sequence_map_viewer imports (embedded map toolbar height).
ZOOM_SLIDER_HEIGHT = 22
ANNOTATION_ZOOM_SLIDER_WIDTH = 240


def apply_annotation_zoom_slider_style(slider: QSlider) -> None:
    """Accent groove + small circular thumb for the annotation map zoom slider."""
    from repeatviewer.ui.theme_mode import load_accent_color

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


def _make_font_size_spin(
    *,
    minimum: int,
    maximum: int,
    tooltip: str,
) -> NoWheelSpinBox:
    """Font-size spinbox for text preview (``pt`` suffix)."""
    spin = NoWheelSpinBox()
    spin.setObjectName("previewZoomSpin")
    spin.setRange(minimum, maximum)
    spin.setSuffix(" pt")
    spin.setFrame(False)
    spin.setAttribute(Qt.WidgetAttribute.WA_MacShowFocusRect, False)
    spin.setFixedHeight(ZOOM_CONTROL_HEIGHT)
    spin.setFixedWidth(ZOOM_SPIN_WIDTH)
    spin.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
    spin.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    spin.setButtonSymbols(QSpinBox.ButtonSymbols.UpDownArrows)
    spin.setToolTip(tooltip)
    spin.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
    spin.setSingleStep(PREVIEW_WHEEL_FONT_STEP_PT)
    spin.setKeyboardTracking(False)
    return spin


def _make_preview_actions_button(*, tooltip: str) -> QToolButton:
    """Bordered Actions menu control matching the preview zoom spinbox."""
    btn = QToolButton()
    btn.setObjectName("PreviewActionsButton")
    btn.setText("Actions")
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
    """Compact preview controls for the Preview title bar.

    Text preview uses a font-size (pt) spin; HTML / PDF / image use zoom %.
    """

    zoom_changed = Signal(int)
    font_size_changed = Signal(int)
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
            tooltip="Preview zoom (Shift+scroll, or trackpad pinch)",
        )
        self._spin.valueChanged.connect(self._on_spin)
        layout.addWidget(self._spin)

        self._font_spin = _make_font_size_spin(
            minimum=MIN_PREVIEW_FONT_SIZE,
            maximum=MAX_PREVIEW_FONT_SIZE,
            tooltip="Preview font size (Shift+scroll, or trackpad pinch)",
        )
        self._font_spin.valueChanged.connect(self._on_font_spin)
        layout.addWidget(self._font_spin)

        self._actions_btn = _make_preview_actions_button(
            tooltip="Preview actions",
        )
        self._actions_menu = QMenu(self._actions_btn)
        self._reset_scale_action = QAction("Reset scale", self._actions_menu)
        self._reset_scale_action.setToolTip("Restore default preview zoom (100%)")
        self._reset_scale_action.triggered.connect(self._reset_zoom)
        self._edit_text_action = QAction("Edit text", self._actions_menu)
        self._edit_text_action.setToolTip("Allow editing the text shown in Preview")
        self._edit_text_action.triggered.connect(self.edit_text_requested.emit)
        self._save_text_action = QAction("Save text", self._actions_menu)
        self._save_text_action.setToolTip("Save Preview text changes to disk")
        self._save_text_action.triggered.connect(self.save_text_requested.emit)
        self._actions_menu.addAction(self._reset_scale_action)
        self._text_actions_sep = self._actions_menu.addSeparator()
        self._actions_menu.addAction(self._edit_text_action)
        self._actions_menu.addAction(self._save_text_action)
        self._default_tab_sep = self._actions_menu.addSeparator()
        self._default_tab_action = _add_default_fasta_tab_menu(self._actions_menu)
        self._actions_btn.setMenu(self._actions_menu)
        layout.addWidget(self._actions_btn)

        self._block = False
        self._fit_mode = load_preview_fit_mode()
        self._zoom_kind = PreviewZoomKind.TEXT
        self._chrome_mode = "text"
        self.set_values(DEFAULT_PREVIEW_ZOOM_PERCENT, self._fit_mode)
        self._block = True
        self._font_spin.setValue(load_preview_font_size())
        self._block = False
        self.set_text_actions_enabled(can_edit=False, can_save=False)
        self.set_chrome_mode("text")
        self.set_default_preview_tab_visible(False)

    def set_chrome_mode(self, mode: str) -> None:
        """Configure title-bar chrome for the active preview kind.

        ``text`` — font-size (pt) spin + Reset / Edit text / Save text.
        ``media`` — % spin + Reset scale (HTML / PDF / image).
        ``dotplot`` — Actions only; no size controls.
        """
        self._chrome_mode = mode
        show_percent = mode == "media"
        show_font = mode == "text"
        show_reset = mode in {"text", "media"}
        show_text_tools = mode == "text"
        self._spin.setVisible(show_percent)
        self._font_spin.setVisible(show_font)
        self._reset_scale_action.setVisible(show_reset)
        self._text_actions_sep.setVisible(show_text_tools)
        self._edit_text_action.setVisible(show_text_tools)
        self._save_text_action.setVisible(show_text_tools)
        if mode == "text":
            self._reset_scale_action.setText("Reset font size")
            self._reset_scale_action.setToolTip(
                f"Restore default preview font size ({DEFAULT_PREVIEW_FONT_SIZE} pt)"
            )
        else:
            self._reset_scale_action.setText("Reset scale")
            self._reset_scale_action.setToolTip("Restore default preview zoom (100%)")

    def set_default_preview_tab_visible(self, visible: bool) -> None:
        """Show Default preview tab only for FASTA files (not plain text/PDF/PNG/HTML)."""
        show = bool(visible)
        self._default_tab_action.setVisible(show)
        self._default_tab_sep.setVisible(show)

    def set_zoom_kind(self, kind: PreviewZoomKind | str) -> None:
        """Load the saved font size (text) or zoom % (media) for this kind."""
        try:
            zoom_kind = kind if isinstance(kind, PreviewZoomKind) else PreviewZoomKind(str(kind))
        except ValueError:
            zoom_kind = PreviewZoomKind.TEXT
        if zoom_kind == PreviewZoomKind.TEXT:
            size = load_preview_font_size()
            if zoom_kind == self._zoom_kind and not self._block:
                if size != self._font_spin.value():
                    self._block = True
                    self._font_spin.setValue(size)
                    self._block = False
                self._zoom_kind = zoom_kind
                return
            self._zoom_kind = zoom_kind
            self._block = True
            self._font_spin.setValue(size)
            self._block = False
            return
        if zoom_kind == self._zoom_kind and not self._block:
            # Still refresh from disk in case Settings changed the value.
            percent = load_preview_zoom_percent(zoom_kind)
            if percent != self._spin.value():
                self._block = True
                self._spin.setValue(percent)
                self._block = False
            return
        self._zoom_kind = zoom_kind
        percent = load_preview_zoom_percent(zoom_kind)
        self._block = True
        self._spin.setValue(percent)
        self._block = False

    def current_zoom_kind(self) -> PreviewZoomKind:
        return self._zoom_kind

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

    def current_font_size(self) -> int:
        return self._font_spin.value()

    def current_fit_mode(self) -> PreviewFitMode:
        return self._fit_mode

    def _size_control_active(self) -> bool:
        """True when text font or media % control should handle Shift+zoom."""
        if self._chrome_mode == "dotplot":
            return False
        if self._chrome_mode == "text":
            return self._font_spin.isEnabled()
        if self._chrome_mode == "media":
            return self._spin.isEnabled()
        return (
            (self._spin.isVisible() and self._spin.isEnabled())
            or (self._font_spin.isVisible() and self._font_spin.isEnabled())
        )

    def _on_spin(self, value: int) -> None:
        if self._block:
            return
        if self._zoom_kind == PreviewZoomKind.TEXT:
            return
        save_preview_zoom_percent(value, self._zoom_kind)
        self.zoom_changed.emit(value)
        # Do not clearFocus here — that stops keyboard ↑/↓ auto-repeat.
        # Stepper mouse clicks clear focus in NoWheelSpinBox.

    def _on_font_spin(self, value: int) -> None:
        if self._block:
            return
        save_preview_font_size(value)
        self.font_size_changed.emit(value)

    def _reset_zoom(self) -> None:
        if self._chrome_mode == "text" or self._zoom_kind == PreviewZoomKind.TEXT:
            self._font_spin.setValue(DEFAULT_PREVIEW_FONT_SIZE)
        else:
            self._spin.setValue(DEFAULT_PREVIEW_ZOOM_PERCENT)

    def step_zoom(self, steps: int) -> None:
        if steps == 0:
            return
        if self._font_spin.isVisible():
            value = self._font_spin.value() + steps * PREVIEW_WHEEL_FONT_STEP_PT
            value = max(MIN_PREVIEW_FONT_SIZE, min(MAX_PREVIEW_FONT_SIZE, value))
            self._font_spin.setValue(value)
            return
        value = self._spin.value() + steps * PREVIEW_WHEEL_ZOOM_STEP_PERCENT
        value = max(MIN_PREVIEW_ZOOM_PERCENT, min(MAX_PREVIEW_ZOOM_PERCENT, value))
        self._spin.setValue(value)

    def refresh_slider_style(self) -> None:
        """Kept for theme-reload callers; spinbox uses the global stylesheet."""


class AnnotationZoomBar(QWidget):
    """Controls beside Text/Annotation tabs.

    Single-sequence Annotation: Actions menu (zoom lives on the map toolbar).
    Multi-sequence overview: −/+/slider, Fit, Extract, and name filter.
    """

    fit_mode_changed = Signal(str)
    reset_requested = Signal()
    fit_requested = Signal()
    extract_requested = Signal()
    filter_text_changed = Signal(str)
    filter_return_pressed = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("PreviewZoomBar")
        self.setFixedHeight(ZOOM_BAR_HEIGHT)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)

        # Hidden spin kept for sequence-map bind API (embedded maps ignore it).
        self._spin = _make_zoom_percent_spin(
            minimum=0,
            maximum=100,
            tooltip="0%: fit full sequence in window. 100%: maximum zoom-in.",
        )
        self._spin.setSingleStep(5)
        self._spin.setValue(0)
        self._spin.hide()
        self._spin.valueChanged.connect(self._on_spin)

        self._overview_host = QWidget(self)
        self._overview_host.setObjectName("AnnotationOverviewControls")
        ov = QHBoxLayout(self._overview_host)
        ov.setContentsMargins(0, 0, 0, 0)
        ov.setSpacing(6)

        zoom_out = QToolButton(self._overview_host)
        zoom_out.setText("−")
        zoom_out.setToolTip("Zoom out (Shift+scroll)")
        zoom_out.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        zoom_out.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        zoom_out.clicked.connect(lambda: self._nudge_zoom_slider(-8))
        zoom_in = QToolButton(self._overview_host)
        zoom_in.setText("+")
        zoom_in.setToolTip("Zoom in (Shift+scroll)")
        zoom_in.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        zoom_in.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        zoom_in.clicked.connect(lambda: self._nudge_zoom_slider(8))

        self._zoom_slider = QSlider(Qt.Orientation.Horizontal, self._overview_host)
        self._zoom_slider.setObjectName("annotationZoomSlider")
        self._zoom_slider.setRange(0, 100)
        self._zoom_slider.setValue(0)
        self._zoom_slider.setFixedWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
        self._zoom_slider.setMinimumWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
        self._zoom_slider.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        self._zoom_slider.setToolTip(
            "Left: fit longest sequence in window. Right: maximum zoom-in."
        )
        apply_annotation_zoom_slider_style(self._zoom_slider)

        fit_btn = QPushButton("Fit", self._overview_host)
        fit_btn.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        fit_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        fit_btn.setToolTip(
            "Fit the longest sequence to the visible width · "
            "Drag on tracks to zoom a region"
        )
        fit_btn.clicked.connect(self.fit_requested.emit)

        extract_btn = QPushButton("Extract", self._overview_host)
        extract_btn.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        extract_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        extract_btn.setToolTip(
            "Extract selected sequences (click / Shift / ⌘, or drag on names)"
        )
        extract_btn.clicked.connect(self.extract_requested.emit)

        self._name_filter = QLineEdit(self._overview_host)
        self._name_filter.setObjectName("FileFilter")
        self._name_filter.setPlaceholderText("Filter by name…")
        self._name_filter.setClearButtonEnabled(True)
        self._name_filter.setFrame(False)
        self._name_filter.setAttribute(Qt.WidgetAttribute.WA_MacShowFocusRect, False)
        self._name_filter.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        self._name_filter.setMaximumHeight(ZOOM_CONTROL_HEIGHT)
        self._name_filter.setMinimumWidth(140)
        self._name_filter.setMaximumWidth(280)
        self._name_filter.setSizePolicy(
            QSizePolicy.Policy.Preferred,
            QSizePolicy.Policy.Fixed,
        )
        self._name_filter.setToolTip(
            "Type part of a sequence name to jump to it. "
            "Press Enter for the next match."
        )
        self._name_filter.textChanged.connect(self.filter_text_changed.emit)
        self._name_filter.returnPressed.connect(self.filter_return_pressed.emit)

        for w in (
            zoom_out,
            zoom_in,
            self._zoom_slider,
            fit_btn,
            extract_btn,
            self._name_filter,
        ):
            ov.addWidget(w, 0, Qt.AlignmentFlag.AlignVCenter)
        layout.addWidget(self._overview_host)

        self._actions_btn = _make_preview_actions_button(
            tooltip="Annotation actions",
        )
        self._actions_menu = QMenu(self._actions_btn)
        self._default_tab_action = _add_default_fasta_tab_menu(self._actions_menu)
        self._actions_btn.setMenu(self._actions_menu)
        layout.addWidget(self._actions_btn)

        self._block = False
        self._overview_mode = False
        self.set_overview_mode(False)

    def set_overview_mode(self, enabled: bool) -> None:
        """Show multi-seq overview controls, or the single-seq Actions menu."""
        self._overview_mode = bool(enabled)
        self._overview_host.setVisible(self._overview_mode)
        self._actions_btn.setVisible(not self._overview_mode)

    def is_overview_mode(self) -> bool:
        return self._overview_mode

    @property
    def slider(self) -> QSpinBox:
        """Zoom percent control (name kept for sequence-map binding)."""
        return self._spin

    @property
    def zoom_slider(self) -> QSlider:
        """Horizontal zoom slider used by the multi-sequence overview."""
        return self._zoom_slider

    @property
    def name_filter(self) -> QLineEdit:
        return self._name_filter

    def sync_display(self, value: int) -> None:
        """Update percent spin without emitting map zoom."""
        self._block = True
        self._spin.blockSignals(True)
        self._spin.setValue(value)
        self._spin.blockSignals(False)
        self._block = False

    def refresh_slider_style(self) -> None:
        """Re-apply accent colors after a theme change."""
        apply_annotation_zoom_slider_style(self._zoom_slider)

    def _nudge_zoom_slider(self, delta: int) -> None:
        self._zoom_slider.setValue(
            max(0, min(100, self._zoom_slider.value() + int(delta)))
        )

    def _on_spin(self, _value: int) -> None:
        if self._block:
            return


class PreviewWheelZoomFilter(QObject):
    """Shift+scroll / pinch / Shift+drag zooms; plain two-finger scroll navigates.

    Watches only an explicit whitelist of text/PDF/image widgets. Never installs
    on QApplication or WebEngine — app-level filters segfault on macOS when an
    HTML preview starts Chromium.
    """

    # Pixels of vertical Shift+drag per one zoom step (percent or pt).
    _DRAG_PIXELS_PER_STEP = 28.0
    _INTERESTING = (
        QEvent.Type.Wheel,
        QEvent.Type.MouseButtonPress,
        QEvent.Type.MouseMove,
        QEvent.Type.MouseButtonRelease,
    )

    def __init__(self, tabs: QWidget, zoom_bar: PreviewZoomBar) -> None:
        # Parent the filter to the zoom bar so its lifetime matches the title-bar
        # controls (not the tab widget, which can be destroyed first).
        super().__init__(zoom_bar)
        self._tabs: QWidget | None = tabs
        self._zoom_bar = zoom_bar
        self._enabled = True
        self._watched: set[int] = set()
        self._zoom_accum = 0
        self._drag_origin: QPointF | None = None
        self._drag_accum = 0.0
        tabs.destroyed.connect(self._disable)
        zoom_bar.destroyed.connect(self._disable)
        self._watch(zoom_bar)

    def register_page(self, page: QWidget) -> None:
        """Attach to safe preview surfaces for this page (never HTML/WebEngine)."""
        if not self._alive() or not isinstance(page, QWidget):
            return
        # Explicit safe targets only — do not walk findChildren / ChildAdded.
        targets: list[QWidget] = [page]
        for attr in (
            "_text",
            "_fasta_text",
            "_text_scroll",
            "_fasta_text_scroll",
            "_placeholder",
        ):
            widget = getattr(page, attr, None)
            if isinstance(widget, QWidget):
                targets.append(widget)
                self._append_viewport(targets, widget)
        # Image/PDF panels wrap a QScrollArea — watch that viewport/label so
        # Shift+scroll/drag reaches us instead of panning the scroll area alone.
        for panel_attr in ("_pdf", "_image"):
            panel = getattr(page, panel_attr, None)
            if not isinstance(panel, QWidget):
                continue
            targets.append(panel)
            scroll = getattr(panel, "_scroll", None)
            if isinstance(scroll, QWidget):
                targets.append(scroll)
                self._append_viewport(targets, scroll)
            for nested_attr in ("_label", "_content"):
                nested = getattr(panel, nested_attr, None)
                if isinstance(nested, QWidget):
                    targets.append(nested)
        for widget in targets:
            self._watch(widget)

    @staticmethod
    def _append_viewport(targets: list[QWidget], widget: QWidget) -> None:
        viewport = getattr(widget, "viewport", None)
        if not callable(viewport):
            return
        try:
            vp = viewport()
        except RuntimeError:
            return
        if isinstance(vp, QWidget):
            targets.append(vp)

    def _disable(self, *_args) -> None:
        if not self._enabled:
            return
        self._enabled = False
        self._tabs = None

    def _alive(self) -> bool:
        if not self._enabled or self._tabs is None:
            return False
        try:
            from shiboken6 import isValid

            return bool(isValid(self) and isValid(self._tabs) and isValid(self._zoom_bar))
        except Exception:
            return False

    def _watch(self, widget: QWidget) -> None:
        if not self._alive() or not isinstance(widget, QWidget):
            return
        name = type(widget).__name__
        if "WebEngine" in name or "RenderWidget" in name or "Html" in name:
            return
        key = id(widget)
        if key in self._watched:
            return
        self._watched.add(key)
        try:
            widget.installEventFilter(self)
        except RuntimeError:
            self._watched.discard(key)

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
        if not self._alive():
            return None
        from repeatviewer.widgets.preview_tabs import FilePreviewPage, PreviewTabWidget

        tabs = self._tabs
        try:
            if not isinstance(tabs, PreviewTabWidget):
                return None
            widget = tabs.currentWidget()
        except RuntimeError:
            self._disable()
            return None
        return widget if isinstance(widget, FilePreviewPage) else None

    def _zoom_modifiers_held(self, event) -> bool:
        mods = event.modifiers() | QGuiApplication.keyboardModifiers()
        if mods & Qt.KeyboardModifier.ShiftModifier:
            return True
        if mods & (
            Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier
        ):
            return True
        return False

    def _page_uses_own_zoom(self, page) -> bool:
        try:
            if page.is_annotation_preview_active() or page.is_dotplot_preview_active():
                return True
            if page.is_html_preview_active():
                return True
        except RuntimeError:
            return True
        return False

    def _can_zoom(self, page, obj: QWidget, event) -> bool:
        if not self._alive():
            return False
        try:
            if not self._zoom_bar.isVisible() or not self._zoom_bar._size_control_active():
                return False
            if self._page_uses_own_zoom(page):
                return False
        except RuntimeError:
            return False
        if isinstance(obj, QSpinBox) and not self._is_descendant_of(
            obj, self._zoom_bar
        ):
            return False
        try:
            if self._is_descendant_of(obj, page._sequence_info_bar):
                return False
            # Never hijack events under the HTML/WebEngine stack entry.
            html = getattr(page, "_html", None)
            if isinstance(html, QWidget) and self._is_descendant_of(obj, html):
                return False
        except RuntimeError:
            return False
        return self._zoom_modifiers_held(event)

    def _apply_zoom_delta(self, delta: int) -> None:
        if delta == 0 or not self._alive():
            return
        self._zoom_accum += delta
        threshold = PREVIEW_WHEEL_DELTA_PER_STEP
        steps = int(self._zoom_accum / threshold)
        if steps != 0:
            self._zoom_accum -= steps * threshold
            try:
                self._zoom_bar.step_zoom(steps)
            except RuntimeError:
                pass

    def eventFilter(self, watched, event) -> bool:  # noqa: N802 (Qt API)
        try:
            if event.type() not in self._INTERESTING:
                return False
        except RuntimeError:
            return False
        if not self._alive():
            return False
        try:
            return self._event_filter_impl(watched, event)
        except RuntimeError:
            return False

    def _event_filter_impl(self, watched, event) -> bool:
        page = self._current_page()
        obj = watched if isinstance(watched, QWidget) else None

        if isinstance(event, QMouseEvent):
            if page is None or obj is None:
                self._drag_origin = None
                self._drag_accum = 0.0
                return False
            if self._page_uses_own_zoom(page):
                return False
            if event.type() == QEvent.Type.MouseButtonPress:
                if (
                    event.button() == Qt.MouseButton.LeftButton
                    and self._can_zoom(page, obj, event)
                ):
                    self._drag_origin = QPointF(event.globalPosition())
                    self._drag_accum = 0.0
                    # Consume press so QScrollArea does not start a pan/drag.
                    event.accept()
                    return True
                self._drag_origin = None
                self._drag_accum = 0.0
                return False
            if event.type() == QEvent.Type.MouseMove:
                if (
                    self._drag_origin is not None
                    and event.buttons() & Qt.MouseButton.LeftButton
                ):
                    # Keep zooming for the active drag even if Shift flickers;
                    # never let the scroll area pan during a zoom drag.
                    pos = QPointF(event.globalPosition())
                    dy = self._drag_origin.y() - pos.y()
                    self._drag_accum += dy
                    self._drag_origin = pos
                    steps = int(self._drag_accum / self._DRAG_PIXELS_PER_STEP)
                    if steps != 0:
                        self._drag_accum -= steps * self._DRAG_PIXELS_PER_STEP
                        self._zoom_bar.step_zoom(steps)
                    event.accept()
                    return True
                return False
            if event.type() == QEvent.Type.MouseButtonRelease:
                if self._drag_origin is not None:
                    self._drag_origin = None
                    self._drag_accum = 0.0
                    event.accept()
                    return True
                self._drag_origin = None
                self._drag_accum = 0.0
                return False

        if event.type() != QEvent.Type.Wheel or not isinstance(event, QWheelEvent):
            return False
        if page is None or obj is None:
            self._zoom_accum = 0
            return False
        if self._page_uses_own_zoom(page):
            self._zoom_accum = 0
            return False
        if not self._can_zoom(page, obj, event):
            self._zoom_accum = 0
            return False
        delta = wheel_zoom_delta_from_event(event)
        self._apply_zoom_delta(delta)
        # Consume fully — on macOS Shift+wheel is also horizontal scroll.
        event.accept()
        return True


FASTA_TAB_STRIP_HEIGHT = 22
PREVIEW_HEADER_HEIGHT = 24
# Tabs row stays FASTA_TAB_STRIP_HEIGHT; controls sit under the divider.
PREVIEW_HEADER_WITH_CONTROLS_HEIGHT = (
    FASTA_TAB_STRIP_HEIGHT + 1 + ZOOM_BAR_HEIGHT + 4
)  # tabs + divider + padded controls


class FastaTabStrip(QWidget):
    """Text / Annotation / Dotplot / Cytogenetics switcher for the preview header."""

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
        for index, label in enumerate(("Text", "Annotation", "Dotplot (Self)", "Cytogenetics")):
            btn = QPushButton(label)
            btn.setObjectName("FastaPreviewTab")
            btn.setCheckable(True)
            btn.setFlat(True)
            btn.setFixedHeight(FASTA_TAB_STRIP_HEIGHT)
            btn.setSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)
            # The selected tab is demi-bold in the stylesheet.  Measure that
            # wider state up front so labels such as "Annotation" never get
            # clipped when the font weight changes after layout.
            selected_font = btn.font()
            selected_font.setWeight(QFont.Weight.DemiBold)
            text_width = QFontMetrics(selected_font).horizontalAdvance(label)
            btn.setMinimumWidth(text_width + 24)
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

    def set_tab_visible(self, index: int, visible: bool) -> None:
        if 0 <= index < len(self._buttons):
            self._buttons[index].setVisible(bool(visible))

    def set_tab_label(self, index: int, label: str) -> None:
        if not (0 <= index < len(self._buttons)):
            return
        btn = self._buttons[index]
        btn.setText(label)
        selected_font = btn.font()
        selected_font.setWeight(QFont.Weight.DemiBold)
        text_width = QFontMetrics(selected_font).horizontalAdvance(label)
        btn.setMinimumWidth(text_width + 24)


class PreviewHeaderBar(QWidget):
    """FASTA Text / Annotation (or Overview) tabs; controls sit below the divider."""

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

        # Tab row keeps the original compact height (never grows with controls).
        tabs_host = QWidget(self)
        tabs_host.setFixedHeight(FASTA_TAB_STRIP_HEIGHT)
        tabs_row = QHBoxLayout(tabs_host)
        tabs_row.setContentsMargins(6, 0, 6, 0)
        tabs_row.setSpacing(4)

        self._fasta_tabs = FastaTabStrip(tabs_host)
        self._fasta_tabs.hide()
        self._fasta_tabs.tab_changed.connect(self.fasta_tab_changed.emit)
        tabs_row.addWidget(self._fasta_tabs)
        tabs_row.addStretch(1)
        outer.addWidget(tabs_host)

        divider = QFrame(self)
        divider.setObjectName("PreviewHeaderDivider")
        divider.setFrameShape(QFrame.Shape.HLine)
        divider.setFixedHeight(1)
        outer.addWidget(divider)

        # −/+/Fit/Extract/filter sit under the line; same control styles as before.
        self._preview_zoom_bar = preview_zoom_bar
        self._annotation_zoom_bar = annotation_zoom_bar
        self._controls_host = QWidget(self)
        self._controls_host.setObjectName("PreviewHeaderControls")
        self._controls_host.setFixedHeight(ZOOM_BAR_HEIGHT + 4)
        controls_row = QHBoxLayout(self._controls_host)
        controls_row.setContentsMargins(6, 2, 6, 2)
        controls_row.setSpacing(4)
        controls_row.addStretch(1)
        if annotation_zoom_bar is not None:
            controls_row.addWidget(annotation_zoom_bar, 0, Qt.AlignmentFlag.AlignVCenter)
        self._controls_host.hide()
        outer.addWidget(self._controls_host)

        self.hide()

    def _update_header_height(self) -> None:
        if self._controls_host.isVisible():
            self.setFixedHeight(PREVIEW_HEADER_WITH_CONTROLS_HEIGHT)
        else:
            self.setFixedHeight(PREVIEW_HEADER_HEIGHT)

    def _update_bar_visibility(self) -> None:
        # Header only shows for Text/Annotation tabs; text zoom lives on the title bar.
        self.setVisible(not self._fasta_tabs.isHidden())
        self._update_header_height()

    def set_header_zoom_mode(self, mode: str) -> None:
        """Show preview zoom, annotation/overview controls, or hide both."""
        show_preview = mode == "preview"
        show_annotation = mode == "annotation"
        self._preview_zoom_bar.setVisible(show_preview)
        if self._annotation_zoom_bar is not None:
            self._annotation_zoom_bar.setVisible(show_annotation)
        self._controls_host.setVisible(
            show_annotation and self._annotation_zoom_bar is not None
        )
        self._update_bar_visibility()

    def set_fasta_tabs_visible(self, visible: bool) -> None:
        self._fasta_tabs.setVisible(visible)
        self._update_bar_visibility()

    def sync_fasta_tab(self, index: int) -> None:
        self._fasta_tabs.set_current_index(index)

    def set_fasta_tab_visible(self, index: int, visible: bool) -> None:
        self._fasta_tabs.set_tab_visible(index, visible)

    def set_fasta_tab_label(self, index: int, label: str) -> None:
        self._fasta_tabs.set_tab_label(index, label)

    def current_fasta_tab(self) -> int:
        return self._fasta_tabs.current_index()
