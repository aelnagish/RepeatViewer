"""Zoom toolbar and unified preview header for the browser panel."""

from __future__ import annotations

from PySide6.QtCore import QEvent, QObject, QPoint, QPointF, QSize, Qt, QTimer, Signal
from PySide6.QtGui import (
    QAction,
    QActionGroup,
    QFont,
    QFontMetrics,
    QGuiApplication,
    QMouseEvent,
    QNativeGestureEvent,
    QWheelEvent,
)
from PySide6.QtWidgets import (
    QButtonGroup,
    QFrame,
    QHBoxLayout,
    QLabel,
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
from repeatviewer.widgets.spin_boxes import NoWheelSlider, NoWheelSpinBox
from repeatviewer.widgets.toolbar_overflow import ToolbarOverflowHider

_FASTA_DEFAULT_TAB_CHOICES: tuple[tuple[str, DefaultFastaPreviewTab], ...] = (
    ("Text", DefaultFastaPreviewTab.TEXT),
    ("Overview", DefaultFastaPreviewTab.ANNOTATION),
    ("Dotplot (Self)", DefaultFastaPreviewTab.DOTPLOT),
)
_MULTI_FASTA_DEFAULT_TAB_CHOICES: tuple[tuple[str, DefaultFastaPreviewTab], ...] = (
    ("Text", DefaultFastaPreviewTab.TEXT),
    ("Overview", DefaultFastaPreviewTab.ANNOTATION),
    ("Dotplot (Self)", DefaultFastaPreviewTab.DOTPLOT),
)


def _add_default_fasta_tab_menu(parent_menu: QMenu) -> QAction:
    """Add Settings-equivalent default FASTA sub-tab chooser to an Actions menu."""
    submenu = parent_menu.addMenu("Default preview tab")
    group = QActionGroup(submenu)
    group.setExclusive(True)
    state = {"multi": False}

    def _rebuild_actions() -> None:
        for action in list(group.actions()):
            group.removeAction(action)
            submenu.removeAction(action)
            action.deleteLater()
        choices = (
            _MULTI_FASTA_DEFAULT_TAB_CHOICES
            if state["multi"]
            else _FASTA_DEFAULT_TAB_CHOICES
        )
        for label, tab in choices:
            action = QAction(label, submenu)
            action.setCheckable(True)
            action.setData(tab.value)
            group.addAction(action)
            submenu.addAction(action)
        if state["multi"]:
            submenu.setToolTip(
                "Sub-tab shown when you single-click a multi-sequence FASTA file."
            )
        else:
            submenu.setToolTip(
                "Sub-tab shown when you single-click a single-sequence FASTA file."
            )
        _sync_checks()

    def _sync_checks() -> None:
        current = load_default_fasta_preview_tab().value
        matched = False
        for action in group.actions():
            checked = action.data() == current
            action.setChecked(checked)
            matched = matched or checked
        if not matched and group.actions():
            # Fall back to Overview/Annotation when saved tab is unavailable.
            for action in group.actions():
                if action.data() == DefaultFastaPreviewTab.ANNOTATION.value:
                    action.setChecked(True)
                    break

    def _on_chosen(action: QAction) -> None:
        try:
            tab = DefaultFastaPreviewTab(str(action.data()))
        except ValueError:
            return
        save_default_fasta_preview_tab(tab)

    def set_multi_fasta_choices(enabled: bool) -> None:
        enabled = bool(enabled)
        if state["multi"] == enabled and group.actions():
            _sync_checks()
            return
        state["multi"] = enabled
        _rebuild_actions()

    submenu.aboutToShow.connect(_sync_checks)
    group.triggered.connect(_on_chosen)
    _rebuild_actions()
    menu_action = submenu.menuAction()
    menu_action.set_multi_fasta_choices = set_multi_fasta_choices  # type: ignore[attr-defined]
    return menu_action

# Title-bar PreviewZoomBar stays compact; overview map toolbar is taller.
ZOOM_BAR_HEIGHT = 18
OVERVIEW_TOOLBAR_HEIGHT = 28
# With the application stylesheet, the single-sequence AnnotationMapToolbar
# renders at 26 px (22 px controls plus 2 px top/bottom margins).
ANNOTATION_TOOLBAR_RENDERED_HEIGHT = 26
# Shared Select / Fit widths (recorded from multi-seq Overview toolbar).
# 55 clipped "Select" on macOS once the buttons became 18 px tall (the small
# bezel has wider side insets): 64 leaves room for the text
ANNOTATION_SELECT_BUTTON_RENDERED_WIDTH = 64
ANNOTATION_FIT_BUTTON_RENDERED_WIDTH = 37
ZOOM_CONTROL_HEIGHT = 18
ZOOM_SPIN_WIDTH = 56
# Kept for sequence_map_viewer imports (embedded map toolbar height).
ZOOM_SLIDER_HEIGHT = 22
ANNOTATION_ZOOM_SLIDER_WIDTH = 120
_OVERVIEW_POSITION_SPIN_WIDTH = 108
_OVERVIEW_SELECT_FIT_HEIGHT = 18
# Public alias for other toolbars (dotplot rail, etc.).
OVERVIEW_SELECT_FIT_HEIGHT = _OVERVIEW_SELECT_FIT_HEIGHT
# −/+ as tall as the Select / Fit buttons and the Position spins (was 22)
_ZOOM_BTN_SIDE = ZOOM_CONTROL_HEIGHT


# macOS draws a native push button with its own bezel at least ~21 px tall, whatever
# fixed height it is given, so toolbar buttons looked taller than the 18 px −/+.
# A style sheet takes the native bezel away: the button is then exactly its height.
COMPACT_BUTTON_QSS = (
    "QPushButton { border: 1px solid palette(mid); border-radius: 4px;"
    " background: palette(base); padding: 0px 8px; margin: 0px;"
    " min-height: 16px; max-height: 16px; }"
    "QPushButton:hover { background: palette(midlight); }"
    "QPushButton:pressed { background: palette(mid); }"
    "QPushButton:disabled { color: palette(mid); }"
)


def style_compact_button(btn: QPushButton) -> None:
    """Toolbar push button (Select, Fit, Clean, …) exactly as tall as the −/+."""
    btn.setStyleSheet(COMPACT_BUTTON_QSS)


def style_compact_spin(spin: QWidget) -> None:
    """Position spin box exactly as tall as the −/+ (the app style sheet's
    min-height made it 2 px taller)."""
    spin.setStyleSheet("QSpinBox { min-height: 16px; max-height: 16px; }")


def make_overview_zoom_button(text: str, parent: QWidget) -> QToolButton:
    """Square −/+ control aligned with the overview zoom slider."""
    return _make_overview_zoom_button(text, parent)


def _make_overview_zoom_button(text: str, parent: QWidget) -> QToolButton:
    """Square −/+ control aligned with the overview zoom slider."""
    btn = QToolButton(parent)
    btn.setText(text)
    btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
    btn.setAutoRaise(False)
    btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
    side = int(_ZOOM_BTN_SIDE)
    btn.setFixedSize(side, side)
    btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
    btn.setStyleSheet("QToolButton { padding: 0px; margin: 0px; }")
    return btn


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


def apply_preview_scale_slider_style(slider: QSlider) -> None:
    """Blue filled track + circular thumb (same look as the annotation zoom bar)."""
    from repeatviewer.ui.theme_mode import load_accent_color

    accent = load_accent_color()
    slider.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
    slider.setStyleSheet(
        f"""
        QSlider#previewScaleSlider::groove:vertical {{
            border: 1px solid palette(mid);
            background: palette(mid);
            width: 4px;
            border-radius: 2px;
        }}
        QSlider#previewScaleSlider::sub-page:vertical {{
            background: {accent};
            border: 1px solid {accent};
            width: 4px;
            border-radius: 2px;
        }}
        QSlider#previewScaleSlider::add-page:vertical {{
            background: palette(mid);
            border: 1px solid palette(mid);
            width: 4px;
            border-radius: 2px;
        }}
        QSlider#previewScaleSlider::handle:vertical {{
            background: {accent};
            border: 1px solid palette(text);
            height: 10px;
            width: 10px;
            margin: 0 -4px;
            border-radius: 5px;
        }}
        QSlider#previewScaleSlider::handle:vertical:hover {{
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


def _make_preview_actions_button(
    *, tooltip: str, parent: QWidget | None = None
) -> QToolButton:
    """Bordered Actions menu control matching the preview zoom spinbox.

    Always pass ``parent`` when possible — a parentless shown QToolButton
    becomes a tiny floating window on macOS.
    """
    btn = QToolButton(parent)
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


class ScaleSliderPopup(QFrame):
    """Vertical slider popup; range matches the linked spin (Shift+drag bounds)."""

    closed = Signal()

    def __init__(self, spin: QSpinBox, anchor: QWidget) -> None:
        super().__init__(
            None,
            Qt.WindowType.Popup | Qt.WindowType.FramelessWindowHint,
        )
        self.setObjectName("previewScaleSliderPopup")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        # Lifetime is owned by ScaleSpinCluster (deleteLater on close) — do not
        # use WA_DeleteOnClose; that raced with Popup dismiss and crashed.
        self._spin = spin
        self._block = False
        self._closed_emitted = False

        self._slider = QSlider(Qt.Orientation.Vertical, self)
        self._slider.setObjectName("previewScaleSlider")
        self._slider.setRange(spin.minimum(), spin.maximum())
        self._slider.setSingleStep(max(1, spin.singleStep()))
        self._slider.setPageStep(max(1, spin.singleStep() * 5))
        # Top = smallest scale / font; bottom = largest.
        self._slider.setInvertedAppearance(True)
        self._slider.setInvertedControls(True)
        self._slider.setMinimumHeight(140)
        self._slider.setFixedWidth(28)
        self._slider.setValue(spin.value())
        apply_preview_scale_slider_style(self._slider)
        self._slider.valueChanged.connect(self._on_slider)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(6, 8, 6, 8)
        layout.setSpacing(0)
        layout.addWidget(self._slider, alignment=Qt.AlignmentFlag.AlignHCenter)

        spin.valueChanged.connect(self._on_spin)
        self._reposition(anchor)

    def paintEvent(self, event) -> None:  # noqa: N802
        # Rounded like every other popup (theme_mode makes the window transparent).
        from repeatviewer.ui.theme_mode import paint_rounded_popup

        paint_rounded_popup(self)

    def hideEvent(self, event) -> None:  # noqa: N802
        super().hideEvent(event)
        # Defer: Popup dismiss runs inside Qt's click handling; syncing the arrow
        # synchronously (or deleteLater here) can exit the process on macOS.
        if not self._closed_emitted:
            self._closed_emitted = True
            QTimer.singleShot(0, self.closed.emit)

    def _reposition(self, anchor: QWidget) -> None:
        self.adjustSize()
        # Drop just below the scale cluster, left-aligned with the arrow.
        origin = anchor.mapToGlobal(QPoint(0, anchor.height() + 2))
        self.move(origin)

    def _on_slider(self, value: int) -> None:
        if self._block:
            return
        self._block = True
        try:
            self._spin.setValue(int(value))
        finally:
            self._block = False

    def _on_spin(self, value: int) -> None:
        if self._block or not self.isVisible():
            return
        self._block = True
        try:
            self._slider.setRange(self._spin.minimum(), self._spin.maximum())
            self._slider.setValue(int(value))
        finally:
            self._block = False


class ScaleSpinCluster(QWidget):
    """Scale spin with a Color-style ▶/▼ arrow that opens a vertical slider."""

    def __init__(self, spin: NoWheelSpinBox, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("previewScaleCluster")
        self.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        self._spin = spin
        self._popup: ScaleSliderPopup | None = None
        # Qt dismisses Popup on press; the same click must not reopen it.
        self._suppress_reopen = False

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        self._arrow = QToolButton(self)
        self._arrow.setObjectName("previewScaleMenuButton")
        # Same filled triangles as the Color disclosure button (▶ / ▼).
        self._arrow.setText("▶")
        self._arrow.setAutoRaise(True)
        self._arrow.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._arrow.setFixedSize(14, ZOOM_CONTROL_HEIGHT)
        self._arrow.setToolTip("Show scale slider")
        self._arrow.installEventFilter(self)
        self._arrow.clicked.connect(self._toggle_slider)
        layout.addWidget(self._arrow)

        spin.setParent(self)
        layout.addWidget(spin)

    @property
    def spin(self) -> NoWheelSpinBox:
        return self._spin

    def eventFilter(self, watched, event) -> bool:  # noqa: N802
        if watched is self._arrow and event.type() == QEvent.Type.MouseButtonPress:
            if (
                isinstance(event, QMouseEvent)
                and event.button() == Qt.MouseButton.LeftButton
                and self._suppress_reopen
            ):
                self._suppress_reopen = False
                return True
        return super().eventFilter(watched, event)

    def _set_arrow_expanded(self, expanded: bool) -> None:
        arrow = getattr(self, "_arrow", None)
        if arrow is None:
            return
        try:
            from shiboken6 import isValid

            if not isValid(self) or not isValid(arrow):
                return
        except ImportError:
            pass
        try:
            arrow.setText("▼" if expanded else "▶")
            arrow.setToolTip(
                "Hide scale slider" if expanded else "Show scale slider"
            )
        except RuntimeError:
            return

    def _clear_suppress_reopen(self) -> None:
        self._suppress_reopen = False

    def _popup_is_open(self) -> bool:
        popup = self._popup
        if popup is None:
            return False
        try:
            from shiboken6 import isValid

            if not isValid(popup):
                self._popup = None
                return False
        except ImportError:
            pass
        try:
            return bool(popup.isVisible())
        except RuntimeError:
            self._popup = None
            return False

    def _close_popup(self) -> None:
        popup = self._popup
        self._popup = None
        self._set_arrow_expanded(False)
        if popup is None:
            return
        try:
            from shiboken6 import isValid

            if not isValid(popup):
                return
        except ImportError:
            pass
        try:
            popup.closed.disconnect(self._on_popup_closed)
        except (RuntimeError, TypeError):
            pass
        try:
            popup.hide()
            popup.deleteLater()
        except RuntimeError:
            pass

    def _open_popup(self) -> None:
        self._close_popup()
        popup = ScaleSliderPopup(self._spin, self)
        self._popup = popup
        popup.closed.connect(self._on_popup_closed)
        self._set_arrow_expanded(True)
        popup.show()
        popup._slider.setFocus(Qt.FocusReason.PopupFocusReason)

    def _on_popup_closed(self, *_args) -> None:
        # Popup already hidden by Qt (outside click / Escape) or about to be.
        popup = self._popup
        self._popup = None
        self._set_arrow_expanded(False)
        self._suppress_reopen = True
        QTimer.singleShot(0, self._clear_suppress_reopen)
        if popup is None:
            return
        try:
            from shiboken6 import isValid

            if not isValid(popup):
                return
        except ImportError:
            pass
        try:
            popup.deleteLater()
        except RuntimeError:
            pass

    def _toggle_slider(self) -> None:
        if self._suppress_reopen:
            self._suppress_reopen = False
            return
        if self._popup_is_open():
            self._close_popup()
            return
        self._open_popup()

    def hideEvent(self, event) -> None:  # noqa: N802
        self._close_popup()
        super().hideEvent(event)


class PreviewZoomBar(QWidget):
    """Compact preview controls for the Preview title bar.

    Text preview uses a font-size (pt) spin; HTML / PDF / image use zoom %.
    """

    zoom_changed = Signal(int)
    font_size_changed = Signal(int)
    fit_mode_changed = Signal(str)
    edit_text_requested = Signal()
    save_text_requested = Signal()
    duplicate_window_requested = Signal()
    split_window_requested = Signal()

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
            tooltip="Preview zoom (Shift+scroll, trackpad pinch, or scale slider)",
        )
        self._spin.valueChanged.connect(self._on_spin)
        self._spin_cluster = ScaleSpinCluster(self._spin, self)
        layout.addWidget(self._spin_cluster)

        self._font_spin = _make_font_size_spin(
            minimum=MIN_PREVIEW_FONT_SIZE,
            maximum=MAX_PREVIEW_FONT_SIZE,
            tooltip="Preview font size (Shift+scroll, trackpad pinch, or scale slider)",
        )
        self._font_spin.valueChanged.connect(self._on_font_spin)
        self._font_cluster = ScaleSpinCluster(self._font_spin, self)
        layout.addWidget(self._font_cluster)

        self._actions_btn = _make_preview_actions_button(
            tooltip="Preview actions",
            parent=self,
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
        self._duplicate_sep = self._actions_menu.addSeparator()
        self._duplicate_action = QAction("Duplicate window", self._actions_menu)
        self._duplicate_action.setToolTip(
            "Open a separate Preview window for this file. "
            "It stays open until you close it and is not replaced when you "
            "click other sequences."
        )
        self._duplicate_action.triggered.connect(self.duplicate_window_requested.emit)
        self._actions_menu.addAction(self._duplicate_action)
        self._split_action = QAction("Split window", self._actions_menu)
        self._split_action.setCheckable(True)
        self._split_action.setToolTip(
            "Show this file in two side-by-side panes inside Preview. "
            "Right-click Text / Overview / Dotplot to open that tab in the "
            "other pane. Choose again to close the split."
        )
        self._split_action.triggered.connect(self.split_window_requested.emit)
        self._actions_menu.addAction(self._split_action)
        self._actions_btn.setMenu(self._actions_menu)
        layout.addWidget(self._actions_btn)
        self._actions_visible = True

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

    def set_actions_visible(self, visible: bool) -> None:
        """Show or hide the Actions button (off in duplicated Preview windows)."""
        self._actions_visible = bool(visible)
        self._actions_btn.setVisible(self._actions_visible)

    def set_split_checked(self, checked: bool) -> None:
        action = getattr(self, "_split_action", None)
        if action is None:
            return
        action.blockSignals(True)
        action.setChecked(bool(checked))
        action.blockSignals(False)

    def set_chrome_mode(self, mode: str) -> None:
        """Configure title-bar chrome for the active preview kind.

        ``text`` — font-size (pt) spin + Reset / Edit text / Save text.
        ``media`` — % spin + Reset scale (HTML / PDF / image).
        ``dotplot`` — Actions only; no size controls (Default preview tab when shown).
        """
        self._chrome_mode = mode
        show_percent = mode == "media"
        show_font = mode == "text"
        show_reset = mode in {"text", "media"}
        show_text_tools = mode == "text"
        self._spin_cluster.setVisible(show_percent)
        self._font_cluster.setVisible(show_font)
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
        self._sync_default_tab_separator()

    def set_default_preview_tab_visible(self, visible: bool) -> None:
        """Show Default preview tab only for FASTA files (not plain text/PDF/PNG/HTML)."""
        show = bool(visible)
        self._default_tab_action.setVisible(show)
        self._sync_default_tab_separator()

    def _sync_default_tab_separator(self) -> None:
        """Hide separators when neighboring Actions items are hidden."""
        show_tab = self._default_tab_action.isVisible()
        something_above = (
            self._reset_scale_action.isVisible()
            or self._edit_text_action.isVisible()
            or self._save_text_action.isVisible()
        )
        self._default_tab_sep.setVisible(show_tab and something_above)
        self._duplicate_sep.setVisible(show_tab or something_above)

    def set_multi_fasta_default_tabs(self, enabled: bool) -> None:
        """Multi-seq: Default preview tab is Text / Overview / Dotplot (Self)."""
        setter = getattr(self._default_tab_action, "set_multi_fasta_choices", None)
        if callable(setter):
            setter(bool(enabled))

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
    """Controls beside Text/Overview tabs.

    Single-sequence Annotation: Actions menu (zoom lives on the map toolbar).
    Multi-sequence overview: same chrome as the map toolbar — −/+/slider,
    Position, Select, Fit, Clean, Cons.
    """

    fit_mode_changed = Signal(str)
    reset_requested = Signal()
    fit_requested = Signal()
    select_requested = Signal()
    clean_requested = Signal()
    consensus_requested = Signal()
    coverage_requested = Signal()
    extract_requested = Signal()  # kept for callers; Extract is not on this bar
    filter_text_changed = Signal(str)
    filter_return_pressed = Signal()
    duplicate_window_requested = Signal()
    split_window_requested = Signal()

    _POSITION_SPIN_WIDTH = _OVERVIEW_POSITION_SPIN_WIDTH

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("PreviewZoomBar")
        # Height matches single-seq AnnotationMapToolbar when hosting overview controls.
        self.setFixedHeight(OVERVIEW_TOOLBAR_HEIGHT)
        # Allow split/duplicate panes to shrink; overflow hides controls instead.
        self.setMinimumWidth(0)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(6, 2, 8, 2)
        layout.setSpacing(6)
        layout.setAlignment(Qt.AlignmentFlag.AlignVCenter)
        layout.setSizeConstraint(QHBoxLayout.SizeConstraint.SetNoConstraint)

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
        # Same padding as SequenceMapWidget AnnotationMapToolbar.
        ov.setContentsMargins(6, 2, 8, 2)
        ov.setSpacing(6)
        ov.setAlignment(Qt.AlignmentFlag.AlignVCenter)

        # Finer than the old ±8 so MSA / multi-seq −/+ are easier to fine-tune
        # (slider range is 0–100).
        zoom_out = _make_overview_zoom_button("−", self._overview_host)
        zoom_out.setToolTip("Zoom out (Shift+scroll)")
        zoom_out.clicked.connect(lambda: self._nudge_zoom_slider(-3))
        zoom_in = _make_overview_zoom_button("+", self._overview_host)
        zoom_in.setToolTip("Zoom in (Shift+scroll)")
        zoom_in.clicked.connect(lambda: self._nudge_zoom_slider(3))

        self._zoom_slider = NoWheelSlider(Qt.Orientation.Horizontal, self._overview_host)
        self._zoom_slider.setObjectName("annotationZoomSlider")
        self._zoom_slider.setRange(0, 100)
        self._zoom_slider.setValue(0)
        self._zoom_slider.setFixedWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
        self._zoom_slider.setMinimumWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
        self._zoom_slider.setFixedHeight(ZOOM_SLIDER_HEIGHT)
        self._zoom_slider.setToolTip(
            "Click the bar: left = fit, right = maximum zoom. "
            "Drag to adjust (keeps a selected region in place)."
        )
        apply_annotation_zoom_slider_style(self._zoom_slider)

        self._pos_label = QLabel("Position:", self._overview_host)
        self._pos_start = NoWheelSpinBox(self._overview_host)
        self._pos_start.setMinimum(1)
        self._pos_start.setMaximum(1)
        self._pos_start.setValue(1)
        self._pos_start.setKeyboardTracking(False)
        self._pos_start.setFixedWidth(self._POSITION_SPIN_WIDTH)
        self._pos_start.setFixedHeight(18)  # same as AnnotationMapToolbar spins
        style_compact_spin(self._pos_start)
        self._pos_start.setAlignment(Qt.AlignmentFlag.AlignRight)
        self._pos_start.setGroupSeparatorShown(True)
        self._pos_start.setToolTip("1-based start position (column)")
        self._pos_start.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        self._to_label = QLabel("to", self._overview_host)
        self._pos_end = NoWheelSpinBox(self._overview_host)
        self._pos_end.setMinimum(1)
        self._pos_end.setMaximum(1)
        self._pos_end.setValue(1)
        self._pos_end.setKeyboardTracking(False)
        self._pos_end.setFixedWidth(self._POSITION_SPIN_WIDTH)
        self._pos_end.setFixedHeight(18)
        style_compact_spin(self._pos_end)
        self._pos_end.setAlignment(Qt.AlignmentFlag.AlignRight)
        self._pos_end.setGroupSeparatorShown(True)
        self._pos_end.setToolTip("1-based end position (column)")
        self._pos_end.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        # Shared Select / Fit widths with single-seq AnnotationMapToolbar.
        select_btn = QPushButton("Select", self._overview_host)
        select_btn.setFixedWidth(ANNOTATION_SELECT_BUTTON_RENDERED_WIDTH)
        select_btn.setFixedHeight(_OVERVIEW_SELECT_FIT_HEIGHT)
        style_compact_button(select_btn)
        select_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        select_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        select_btn.setToolTip(
            "Zoom to the Position range and center it in the window"
        )
        select_btn.clicked.connect(self.select_requested.emit)
        self._select_btn = select_btn

        fit_btn = QPushButton("Fit", self._overview_host)
        fit_btn.setFixedWidth(ANNOTATION_FIT_BUTTON_RENDERED_WIDTH)
        fit_btn.setFixedHeight(_OVERVIEW_SELECT_FIT_HEIGHT)
        style_compact_button(fit_btn)
        fit_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        fit_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        fit_btn.setToolTip(
            "Fit the longest sequence to the visible width · "
            "Drag on tracks to select a region"
        )
        fit_btn.clicked.connect(self.fit_requested.emit)
        self._fit_btn = fit_btn

        clean_btn = QPushButton("Clean", self._overview_host)
        clean_btn.setFixedHeight(_OVERVIEW_SELECT_FIT_HEIGHT)
        style_compact_button(clean_btn)
        clean_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        clean_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        clean_btn.setToolTip(
            "Clean this alignment (opens MSA clean parameters)"
        )
        clean_btn.clicked.connect(self.clean_requested.emit)
        self._clean_btn = clean_btn

        consensus_btn = QPushButton("Cons", self._overview_host)
        consensus_btn.setFixedHeight(_OVERVIEW_SELECT_FIT_HEIGHT)
        style_compact_button(consensus_btn)
        consensus_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        consensus_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        consensus_btn.setToolTip(
            "Build a consensus sequence for this alignment "
            "(opens Consensus parameters)"
        )
        consensus_btn.clicked.connect(self.consensus_requested.emit)
        self._consensus_btn = consensus_btn

        coverage_btn = QPushButton("Cov", self._overview_host)
        coverage_btn.setFixedHeight(_OVERVIEW_SELECT_FIT_HEIGHT)
        style_compact_button(coverage_btn)
        coverage_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        coverage_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        coverage_btn.setToolTip(
            "Build a consensus (threshold 0), BLAST it against the project "
            "genome, and plot TE-Aid-style coverage on top of the MSA."
        )
        coverage_btn.clicked.connect(self.coverage_requested.emit)
        self._coverage_btn = coverage_btn

        # Hidden: overview name jump still available via bind API / Find.
        self._name_filter = QLineEdit(self._overview_host)
        self._name_filter.setObjectName("FileFilter")
        self._name_filter.setPlaceholderText("Filter by name…")
        self._name_filter.setClearButtonEnabled(True)
        self._name_filter.hide()
        self._name_filter.textChanged.connect(self.filter_text_changed.emit)
        self._name_filter.returnPressed.connect(self.filter_return_pressed.emit)

        self._sequence_count_label = QLabel("", self._overview_host)
        self._sequence_count_label.setObjectName("OverviewSequenceCount")
        self._sequence_count_label.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        self._sequence_count_label.setToolTip("Number of sequences in this file")
        self._sequence_count_label.hide()

        for w in (
            zoom_out,
            zoom_in,
            self._zoom_slider,
            self._pos_label,
            self._pos_start,
            self._to_label,
            self._pos_end,
            select_btn,
            fit_btn,
            clean_btn,
            consensus_btn,
            coverage_btn,
        ):
            ov.addWidget(w, 0, Qt.AlignmentFlag.AlignVCenter)
        # Left-pack controls; sequence count sits on the right after stretch.
        ov.addStretch(1)
        ov.addWidget(
            self._sequence_count_label, 0, Qt.AlignmentFlag.AlignVCenter
        )
        self._overview_host.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._overview_host.setMinimumWidth(0)
        self._overview_host.setFixedHeight(OVERVIEW_TOOLBAR_HEIGHT)
        # Hide whole controls from the right when the bar is too narrow
        # (−/+ stay longest; never show a half-cut button).
        self._overview_overflow = ToolbarOverflowHider(
            self._overview_host,
            [
                [zoom_out],
                [zoom_in],
                [self._zoom_slider],
                [
                    self._pos_label,
                    self._pos_start,
                    self._to_label,
                    self._pos_end,
                ],
                [select_btn],
                [fit_btn],
                [clean_btn],
                [consensus_btn],
                [coverage_btn],
                [self._sequence_count_label],
            ],
            spacing=6,
            layout=ov,
        )
        # Count / Clean / Cons / Cov start disallowed; set_overview_mode / set_msa_mode enable.
        self._overview_overflow.set_unit_allowed_by_widget(
            self._sequence_count_label, False
        )
        self._overview_overflow.set_unit_allowed_by_widget(clean_btn, False)
        self._overview_overflow.set_unit_allowed_by_widget(consensus_btn, False)
        self._overview_overflow.set_unit_allowed_by_widget(coverage_btn, False)
        layout.addWidget(self._overview_host, 1)

        self._actions_btn = _make_preview_actions_button(
            tooltip="Annotation actions",
            parent=self,
        )
        self._actions_menu = QMenu(self._actions_btn)
        self._default_tab_action = _add_default_fasta_tab_menu(self._actions_menu)
        self._actions_menu.addSeparator()
        self._duplicate_action = QAction("Duplicate window", self._actions_menu)
        self._duplicate_action.setToolTip(
            "Open a separate Preview window for this file. "
            "It stays open until you close it and is not replaced when you "
            "click other sequences."
        )
        self._duplicate_action.triggered.connect(self.duplicate_window_requested.emit)
        self._actions_menu.addAction(self._duplicate_action)
        self._split_action = QAction("Split window", self._actions_menu)
        self._split_action.setCheckable(True)
        self._split_action.setToolTip(
            "Show this file in two side-by-side panes inside Preview. "
            "Right-click Text / Overview / Dotplot to open that tab in the "
            "other pane. Choose again to close the split."
        )
        self._split_action.triggered.connect(self.split_window_requested.emit)
        self._actions_menu.addAction(self._split_action)
        self._actions_btn.setMenu(self._actions_menu)
        layout.addWidget(self._actions_btn, 0, Qt.AlignmentFlag.AlignVCenter)
        self._actions_visible = True
        # Split panes: keep Actions on the Overview strip (title-bar Actions
        # still work; this mirrors Split window next to Fit/Cons).
        self._force_actions_in_overview = False

        self._block = False
        self._overview_mode = False
        self._msa_mode = False
        self.set_overview_mode(False)

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        # Do not let Position/Cons/… fixed widths set the window/splitter floor;
        # ToolbarOverflowHider drops whole controls when this bar is narrow.
        return QSize(0, self.height() or OVERVIEW_TOOLBAR_HEIGHT)

    def sizeHint(self) -> QSize:  # noqa: N802
        height = self.height() or OVERVIEW_TOOLBAR_HEIGHT
        width = max(0, self.width()) if self.width() > 0 else 240
        return QSize(width, height)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        hider = getattr(self, "_overview_overflow", None)
        if hider is not None:
            hider.sync()

    def set_actions_visible(self, visible: bool) -> None:
        """Show or hide the Actions button (off in duplicated Preview windows)."""
        self._actions_visible = bool(visible)
        self._sync_actions_btn_visibility()

    def set_force_actions_in_overview(self, enabled: bool) -> None:
        """Keep Actions on the Overview strip (Split window panes)."""
        self._force_actions_in_overview = bool(enabled)
        self._sync_actions_btn_visibility()

    def _sync_actions_btn_visibility(self) -> None:
        show = bool(self._actions_visible) and (
            not self._overview_mode or bool(self._force_actions_in_overview)
        )
        self._actions_btn.setVisible(show)

    def set_split_checked(self, checked: bool) -> None:
        action = getattr(self, "_split_action", None)
        if action is None:
            return
        action.blockSignals(True)
        action.setChecked(bool(checked))
        action.blockSignals(False)

    def set_msa_mode(self, enabled: bool) -> None:
        """Show alignment-only actions only for a genuine aligned MSA."""
        self._msa_mode = bool(enabled)
        visible = self._overview_mode and self._msa_mode
        # Overflow hider owns visibility; mark Clean/Cons/Cov as allowed only for MSA.
        hider = getattr(self, "_overview_overflow", None)
        if hider is not None:
            hider.set_unit_allowed_by_widget(self._clean_btn, visible)
            hider.set_unit_allowed_by_widget(self._consensus_btn, visible)
            hider.set_unit_allowed_by_widget(self._coverage_btn, visible)
        else:
            self._clean_btn.setVisible(visible)
            self._consensus_btn.setVisible(visible)
            self._coverage_btn.setVisible(visible)

    def set_sequence_count(self, count: int) -> None:
        """Show total sequences on the right (same row as Cons)."""
        n = max(0, int(count))
        hider = getattr(self, "_overview_overflow", None)
        if n <= 0:
            self._sequence_count_label.clear()
            if hider is not None:
                hider.set_unit_allowed_by_widget(self._sequence_count_label, False)
            else:
                self._sequence_count_label.hide()
            return
        label = "1 sequence" if n == 1 else f"{n:,} sequences"
        self._sequence_count_label.setText(label)
        if hider is not None:
            hider.set_unit_allowed_by_widget(
                self._sequence_count_label, self._overview_mode
            )
        else:
            self._sequence_count_label.setVisible(self._overview_mode)

    def set_overview_mode(self, enabled: bool) -> None:
        """Show multi-seq −/+/Position/Select/Fit/Clean/Cons; Actions stay on Preview title."""
        self._overview_mode = bool(enabled)
        self._overview_host.setVisible(self._overview_mode)
        self.set_msa_mode(self._msa_mode)
        if self._sequence_count_label.text().strip():
            self._overview_overflow.set_unit_allowed_by_widget(
                self._sequence_count_label, self._overview_mode
            )
        else:
            self._overview_overflow.set_unit_allowed_by_widget(
                self._sequence_count_label, False
            )
        # Title-bar PreviewZoomBar owns Actions (aligned with "Preview"), matching
        # single-seq Annotation. Split panes force Actions onto this strip too.
        self._sync_actions_btn_visibility()
        if self._overview_mode:
            # Match AnnotationMapToolbar chrome identity for theme / sizing.
            self.setObjectName("AnnotationMapToolbar")
            self.setFixedHeight(ANNOTATION_TOOLBAR_RENDERED_HEIGHT)
            self._overview_host.setFixedHeight(ANNOTATION_TOOLBAR_RENDERED_HEIGHT)
            # Zero outer margins — overview_host already matches map toolbar padding.
            lay = self.layout()
            if lay is not None:
                lay.setContentsMargins(0, 0, 0, 0)
                lay.setSpacing(0)
        else:
            self.setObjectName("PreviewZoomBar")
            self.setFixedHeight(OVERVIEW_TOOLBAR_HEIGHT)
            self._overview_host.setFixedHeight(OVERVIEW_TOOLBAR_HEIGHT)
            lay = self.layout()
            if lay is not None:
                lay.setContentsMargins(6, 2, 8, 2)
                lay.setSpacing(6)
        setter = getattr(self._default_tab_action, "set_multi_fasta_choices", None)
        if callable(setter):
            setter(self._overview_mode)
        # The MSA bar becomes ``AnnotationMapToolbar`` only after construction
        # and reparenting. Qt does not automatically re-evaluate descendant
        # stylesheet selectors after that object-name change, so its controls
        # otherwise retain different native sizes from the single-sequence
        # toolbar. Re-polish the whole subtree to apply the same toolbar rules.
        for widget in (self, *self.findChildren(QWidget)):
            style = widget.style()
            style.unpolish(widget)
            style.polish(widget)
            widget.updateGeometry()
        self._overview_overflow.sync()
        self.update()

    def set_multi_fasta_default_tabs(self, enabled: bool) -> None:
        """Multi-seq: Default preview tab is Text / Overview / Dotplot (Self)."""
        setter = getattr(self._default_tab_action, "set_multi_fasta_choices", None)
        if callable(setter):
            setter(bool(enabled))

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

    @property
    def pos_start(self) -> QSpinBox:
        return self._pos_start

    @property
    def pos_end(self) -> QSpinBox:
        return self._pos_end

    def set_position_range(self, maximum: int, *, start: int = 1, end: int | None = None) -> None:
        """Update Position spin limits for the loaded alignment / longest sequence."""
        maximum = max(1, int(maximum))
        start = max(1, min(maximum, int(start)))
        end_val = maximum if end is None else max(1, min(maximum, int(end)))
        if end_val < start:
            end_val = start
        for spin in (self._pos_start, self._pos_end):
            spin.blockSignals(True)
            spin.setMaximum(maximum)
            spin.blockSignals(False)
        self._pos_start.blockSignals(True)
        self._pos_start.setValue(start)
        self._pos_start.blockSignals(False)
        self._pos_end.blockSignals(True)
        self._pos_end.setValue(end_val)
        self._pos_end.blockSignals(False)

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
    # Trackpad pinch magnification accumulated per zoom step.
    _PINCH_DELTA_PER_STEP = 0.04
    _INTERESTING = (
        QEvent.Type.Wheel,
        QEvent.Type.MouseButtonPress,
        QEvent.Type.MouseMove,
        QEvent.Type.MouseButtonRelease,
        QEvent.Type.NativeGesture,
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
        self._pinch_accum = 0.0
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

    def _can_pinch_zoom(self, page, obj: QWidget | None) -> bool:
        """Pinch does not need Shift; still skip Annotation / HTML / foreign widgets."""
        if not self._alive() or obj is None:
            return False
        try:
            if not self._zoom_bar.isVisible() or not self._zoom_bar._size_control_active():
                return False
            if self._page_uses_own_zoom(page):
                return False
        except RuntimeError:
            return False
        try:
            if self._is_descendant_of(obj, page._sequence_info_bar):
                return False
            html = getattr(page, "_html", None)
            if isinstance(html, QWidget) and self._is_descendant_of(obj, html):
                return False
        except RuntimeError:
            return False
        return True

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

    def _apply_pinch_delta(self, magnification: float) -> None:
        if magnification == 0.0 or not self._alive():
            return
        self._pinch_accum += float(magnification)
        steps = int(self._pinch_accum / self._PINCH_DELTA_PER_STEP)
        if steps != 0:
            self._pinch_accum -= steps * self._PINCH_DELTA_PER_STEP
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

        if event.type() == QEvent.Type.NativeGesture and isinstance(
            event, QNativeGestureEvent
        ):
            if page is None or obj is None:
                return False
            if event.gestureType() != Qt.NativeGestureType.ZoomNativeGesture:
                return False
            if not self._can_pinch_zoom(page, obj):
                return False
            # value() is the magnification delta (positive = zoom in).
            self._apply_pinch_delta(float(event.value()))
            event.accept()
            return True

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
    """Text / Overview / Dotplot / Cytogenetics switcher for the preview header."""

    tab_changed = Signal(int)
    # Open this tab index in the left or right split pane (``"left"`` / ``"right"``).
    openInPanelRequested = Signal(int, str)
    closeSplitRequested = Signal()

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
        self._split_active = False
        for index, label in enumerate(("Text", "Overview", "Dotplot (Self)", "Cytogenetics")):
            btn = QPushButton(label)
            btn.setObjectName("FastaPreviewTab")
            btn.setCheckable(True)
            btn.setFlat(True)
            btn.setFixedHeight(FASTA_TAB_STRIP_HEIGHT)
            btn.setSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)
            # The selected tab is demi-bold in the stylesheet.  Measure that
            # wider state up front so labels such as "Overview" never get
            # clipped when the font weight changes after layout.
            selected_font = btn.font()
            selected_font.setWeight(QFont.Weight.DemiBold)
            text_width = QFontMetrics(selected_font).horizontalAdvance(label)
            btn.setMinimumWidth(text_width + 24)
            btn.setCursor(Qt.CursorShape.PointingHandCursor)
            btn.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
            btn.customContextMenuRequested.connect(
                lambda pos, idx=index, b=btn: self._show_tab_context_menu(b, pos, idx)
            )
            self._group.addButton(btn, index)
            layout.addWidget(btn)
            self._buttons.append(btn)
        self._buttons[0].setChecked(True)
        self._group.idClicked.connect(self.tab_changed.emit)

    def set_split_active(self, active: bool) -> None:
        """Show Close split in the tab context menu when the Preview is split."""
        self._split_active = bool(active)

    def _show_tab_context_menu(
        self, button: QPushButton, pos: QPoint, tab_index: int
    ) -> None:
        if not button.isVisible():
            return
        menu = QMenu(button)
        open_right = menu.addAction("Open in right panel")
        open_left = menu.addAction("Open in left panel")
        close_split = None
        if self._split_active:
            menu.addSeparator()
            close_split = menu.addAction("Close split")
        chosen = menu.exec(button.mapToGlobal(pos))
        if chosen is open_right:
            self.openInPanelRequested.emit(int(tab_index), "right")
        elif chosen is open_left:
            self.openInPanelRequested.emit(int(tab_index), "left")
        elif close_split is not None and chosen is close_split:
            self.closeSplitRequested.emit()

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
    """FASTA Text / Overview tabs; controls sit below the divider."""

    fasta_tab_changed = Signal(int)
    openInPanelRequested = Signal(int, str)
    closeSplitRequested = Signal()

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
        self._fasta_tabs.openInPanelRequested.connect(self.openInPanelRequested.emit)
        self._fasta_tabs.closeSplitRequested.connect(self.closeSplitRequested.emit)
        tabs_row.addWidget(self._fasta_tabs)
        tabs_row.addStretch(1)
        outer.addWidget(tabs_host)

        divider = QFrame(self)
        divider.setObjectName("PreviewHeaderDivider")
        divider.setFrameShape(QFrame.Shape.HLine)
        divider.setFixedHeight(1)
        outer.addWidget(divider)

        # −/+/Position/Select/Fit sit under the line; same control styles as before.
        self._preview_zoom_bar = preview_zoom_bar
        self._annotation_zoom_bar = annotation_zoom_bar
        self._controls_host = QWidget(self)
        self._controls_host.setObjectName("PreviewHeaderControls")
        self._controls_host.setFixedHeight(ZOOM_BAR_HEIGHT + 4)
        controls_row = QHBoxLayout(self._controls_host)
        controls_row.setContentsMargins(6, 2, 6, 2)
        controls_row.setSpacing(4)
        # Left-align −/+/Position/Select/Fit (same side as Text/Overview tabs).
        if annotation_zoom_bar is not None:
            controls_row.addWidget(annotation_zoom_bar, 0, Qt.AlignmentFlag.AlignVCenter)
        controls_row.addStretch(1)
        self._controls_host.hide()
        outer.addWidget(self._controls_host)

        self.hide()

    def _update_header_height(self) -> None:
        if self._controls_host.isVisible():
            self.setFixedHeight(PREVIEW_HEADER_WITH_CONTROLS_HEIGHT)
        else:
            self.setFixedHeight(PREVIEW_HEADER_HEIGHT)

    def _update_bar_visibility(self) -> None:
        # Header only shows for Text/Overview tabs; text zoom lives on the title bar.
        self.setVisible(not self._fasta_tabs.isHidden())
        self._update_header_height()

    def set_header_zoom_mode(self, mode: str) -> None:
        """Show preview zoom, annotation/overview controls, or hide both.

        ``tabs`` — FASTA tabs only; AnnotationZoomBar is hosted inside the
        multi-seq overview panel (next to the rail), matching single-seq layout.
        ``hidden`` — hide title zoom chrome but keep Text/Overview tabs (used by
        duplicated Overview windows that have no Actions button).
        """
        show_preview = mode == "preview"
        show_annotation = mode == "annotation"
        tabs_only = mode in {"tabs", "hidden"}
        self._preview_zoom_bar.setVisible(show_preview)
        bar = self._annotation_zoom_bar
        bar_alive = False
        if bar is not None:
            try:
                from shiboken6 import isValid

                bar_alive = bool(isValid(bar))
            except Exception:
                bar_alive = False
            if not bar_alive:
                self._annotation_zoom_bar = None
                bar = None
        if bar is not None and not tabs_only:
            # Multi-sequence Overview reparents the shared bar into its own
            # toolbar slot.  Header refreshes still run after that mount (for
            # example while an MSA is opening); hiding the bar here would then
            # hide the controls inside the active Overview.  Only the header
            # owns the bar's visibility while no Overview is bound to it.
            if getattr(bar, "_bound_overview", None) is None:
                bar.setVisible(show_annotation)
        self._controls_host.setVisible(
            show_annotation and bar is not None and not tabs_only
        )
        self._update_bar_visibility()

    def annotation_controls_host(self) -> QWidget:
        """Parent that normally holds AnnotationZoomBar under the tab divider."""
        return self._controls_host

    def set_fasta_tabs_visible(self, visible: bool) -> None:
        self._fasta_tabs.setVisible(visible)
        self._update_bar_visibility()

    def set_split_active(self, active: bool) -> None:
        self._fasta_tabs.set_split_active(active)

    def sync_fasta_tab(self, index: int) -> None:
        self._fasta_tabs.set_current_index(index)

    def set_fasta_tab_visible(self, index: int, visible: bool) -> None:
        self._fasta_tabs.set_tab_visible(index, visible)

    def set_fasta_tab_label(self, index: int, label: str) -> None:
        self._fasta_tabs.set_tab_label(index, label)

    def current_fasta_tab(self) -> int:
        return self._fasta_tabs.current_index()
