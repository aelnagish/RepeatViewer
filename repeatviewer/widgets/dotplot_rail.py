"""Geneious-style right rail for the Dotplot preview (Display parameters)."""

from __future__ import annotations

from PySide6.QtCore import QPoint, QRectF, QSettings, QSize, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPalette, QPen
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QStackedWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.dotplot import (
    DEFAULT_WORD_SIZE,
    MAX_WORD_SIZE,
    MIN_WORD_SIZE,
)
from repeatviewer.backends.fastga_dotplot import (
    DEFAULT_MIN_IDENTITY,
    DEFAULT_MIN_LENGTH,
    DEFAULT_SEED_FREQ,
    MAX_SEED_FREQ,
)
from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox
from repeatviewer.widgets.preview_controls import (
    OVERVIEW_SELECT_FIT_HEIGHT,
    ZOOM_SLIDER_HEIGHT,
    apply_annotation_zoom_slider_style,
    make_overview_zoom_button,
)
from repeatviewer.widgets.spin_boxes import (
    NoWheelDoubleSpinBox as QDoubleSpinBox,
    NoWheelSlider,
    NoWheelSpinBox as QSpinBox,
)

_RAIL_TAB_WIDTH = 24
_MODULE_WIDTH = 260
_TAB_BUTTON_HEIGHT_MIN = 36
_TAB_PILL_RADIUS = 7
_TAB_TEXT_PAD = 4
_COLLAPSED_RAIL_WIDTH = _RAIL_TAB_WIDTH

_SETTINGS_GROUP = "dotplot_rail"
_SETTINGS_COLLAPSED_KEY = "panel_collapsed"

_DEFAULT_LINE_WIDTH_STEP = 3
_MIN_LINE_WIDTH_STEP = 1
_MAX_LINE_WIDTH_STEP = 15

def _section_label(text: str, parent: QWidget | None = None) -> QLabel:
    label = QLabel(text, parent)
    font = label.font()
    font.setBold(True)
    label.setFont(font)
    return label


def _h_divider(parent: QWidget | None = None) -> QFrame:
    line = QFrame(parent)
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Sunken)
    line.setFixedHeight(2)
    return line


def _page_title_with_rule(title: str, parent: QWidget) -> QWidget:
    host = QWidget(parent)
    host.setMinimumWidth(0)
    column = QVBoxLayout(host)
    column.setContentsMargins(0, 0, 0, 0)
    column.setSpacing(0)
    title_wrap = QWidget(host)
    title_layout = QVBoxLayout(title_wrap)
    title_layout.setContentsMargins(0, 0, 0, 6)
    title_layout.setSpacing(0)
    title_layout.addWidget(_section_label(title, title_wrap))
    column.addWidget(title_wrap)
    column.addWidget(_h_divider(host))
    return host


class _RotatedTabButton(QToolButton):
    """Vertical pill tab (label rotated −90°)."""

    def __init__(self, label: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._label = label
        self.setObjectName("DotplotRailTabButton")
        self.setText(label)
        self.setCheckable(True)
        self.setAutoRaise(False)
        self.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setMouseTracking(True)
        self.toggled.connect(lambda _checked: self.update())
        fm = QFontMetrics(self.font())
        height = max(
            _TAB_BUTTON_HEIGHT_MIN,
            fm.horizontalAdvance(label) + 2 * _TAB_TEXT_PAD,
        )
        self.setFixedSize(_RAIL_TAB_WIDTH - 2, height)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

    def sizeHint(self) -> QSize:  # noqa: N802
        return self.minimumSize()

    def enterEvent(self, event) -> None:  # noqa: N802
        super().enterEvent(event)
        self.update()

    def leaveEvent(self, event) -> None:  # noqa: N802
        super().leaveEvent(event)
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        del event
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
        pal = self.palette()
        rect = QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5)
        if self.isChecked():
            fill = pal.color(pal.ColorRole.Highlight)
            text = pal.color(pal.ColorRole.HighlightedText)
            border = fill.darker(120)
            pen_w = 1.5
        elif self.underMouse():
            fill = pal.color(pal.ColorRole.AlternateBase)
            text = pal.color(pal.ColorRole.ButtonText)
            border = pal.color(pal.ColorRole.Mid)
            pen_w = 1.0
        else:
            fill = pal.color(pal.ColorRole.Base)
            text = pal.color(pal.ColorRole.ButtonText)
            border = pal.color(pal.ColorRole.Mid)
            pen_w = 1.0
        if not border.isValid() or border.alpha() == 0:
            border = QColor(140, 140, 140)
        painter.setPen(QPen(border, pen_w))
        painter.setBrush(fill)
        painter.drawRoundedRect(rect, _TAB_PILL_RADIUS, _TAB_PILL_RADIUS)
        font = QFont(self.font())
        font.setBold(self.isChecked())
        painter.setPen(text)
        painter.setFont(font)
        fm = QFontMetrics(font)
        painter.translate(self.width() / 2.0, self.height() / 2.0)
        painter.rotate(-90)
        text_w = fm.horizontalAdvance(self._label)
        text_h = fm.ascent()
        painter.drawText(
            QPoint(int(-text_w / 2), int(text_h / 2) - 1),
            self._label,
        )
        painter.end()


class _AlignedRangeField(QFrame):
    """Read-only range box styled like a spin/line edit, with a shared start column.

    Layout inside the box is ``{start} - {end}`` (single spaces around ``-``).
    The start label is right-aligned at a width shared with its sibling field so
    the hyphens on the x and y rows line up vertically.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("dotplotRangeField")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.setFixedHeight(OVERVIEW_SELECT_FIT_HEIGHT)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setStyleSheet(
            "QFrame#dotplotRangeField {"
            "  background-color: palette(base);"
            "  border: 1px solid palette(mid);"
            "  border-radius: 4px;"
            "}"
            "QFrame#dotplotRangeField QLabel {"
            "  background: transparent;"
            "  border: none;"
            "  color: palette(text);"
            "}"
        )
        row = QHBoxLayout(self)
        row.setContentsMargins(6, 0, 6, 0)
        row.setSpacing(0)
        self._start = QLabel(self)
        self._start.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        self._sep = QLabel(" - ", self)
        self._sep.setAlignment(
            Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignVCenter
        )
        self._end = QLabel(self)
        self._end.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        )
        row.addWidget(self._start, 0)
        row.addWidget(self._sep, 0)
        row.addWidget(self._end, 1)
        self._has_range = False
        self.clear()

    def clear(self) -> None:
        """Leave the box blank when there is no selection (no placeholder dash)."""
        self._has_range = False
        self._start.clear()
        self._end.clear()
        self._start.hide()
        self._sep.hide()
        self._end.hide()

    def set_range(self, start: int, end: int) -> None:
        self._has_range = True
        self._start.setText(f"{int(start):,}")
        self._end.setText(f"{int(end):,}")
        self._start.show()
        self._sep.show()
        self._end.show()

    def start_text(self) -> str:
        return self._start.text() if self._has_range else ""

    def set_start_width(self, width: int) -> None:
        """Fix the start column so sibling fields share the same hyphen x-position."""
        width = max(0, int(width))
        if width <= 0:
            self._start.setMinimumWidth(0)
            self._start.setMaximumWidth(16777215)
            return
        self._start.setFixedWidth(width)


class DotplotRail(QWidget):
    """Vertical Display tab + flat parameter panel for dotplot controls."""

    panelCollapsedChanged = Signal(bool)
    computeClicked = Signal()
    cancelClicked = Signal()
    zoomSliderChanged = Signal(int)
    viewModeChanged = Signal(str)  # "single" | "grid"
    galleryPageStep = Signal(int)  # -1 / +1
    galleryPageSelected = Signal(int)  # 0-based page index from dropdown
    galleryPageSizeChanged = Signal(int)
    galleryThumbSizeChanged = Signal(int)
    returnToGridClicked = Signal()

    TAB_DISPLAY = 0

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("dotplotRail")
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Expanding)
        self._module_width = _MODULE_WIDTH
        self._panel_collapsed = False
        self._active_tab = self.TAB_DISPLAY
        self._updating_tabs = False
        self._updating_zoom_slider = False
        self.setFixedWidth(self._expanded_rail_width())

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        self._stack = QStackedWidget(self)
        self._stack.setFixedWidth(self._module_width)
        self._stack.setSizePolicy(
            QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Expanding
        )
        self._stack.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self._stack.setStyleSheet("QStackedWidget { background: transparent; }")

        self._display_page = self._build_display_page()
        self._stack.addWidget(self._display_page)

        self._divider = QFrame(self)
        self._divider.setFrameShape(QFrame.Shape.VLine)
        self._divider.setFrameShadow(QFrame.Shadow.Plain)
        self._divider.setFixedWidth(1)
        self._divider.setStyleSheet("color: palette(mid);")

        self._tab_strip = QWidget(self)
        self._tab_strip.setFixedWidth(_RAIL_TAB_WIDTH)
        self._tab_strip.setObjectName("dotplotRailTabs")
        tab_layout = QVBoxLayout(self._tab_strip)
        tab_layout.setContentsMargins(1, 2, 1, 2)
        tab_layout.setSpacing(2)

        self._tab_group = QButtonGroup(self)
        self._tab_group.setExclusive(False)
        self._tab_buttons: list[_RotatedTabButton] = []
        btn = _RotatedTabButton("Display", self._tab_strip)
        btn.setToolTip("Display — plot parameters (click again to hide)")
        self._tab_group.addButton(btn, self.TAB_DISPLAY)
        tab_layout.addWidget(btn, 0, Qt.AlignmentFlag.AlignHCenter)
        self._tab_buttons.append(btn)
        btn.clicked.connect(
            lambda _checked=False, i=self.TAB_DISPLAY: self._on_tab_clicked(i)
        )
        tab_layout.addStretch(1)

        root.addWidget(self._stack, 0)
        root.addWidget(self._divider, 0)
        root.addWidget(self._tab_strip, 0)
        self._restore_rail_ui_state()

    def _expanded_rail_width(self) -> int:
        return int(self._module_width) + _RAIL_TAB_WIDTH + 1

    def _fill_rail_surface(self, widget: QWidget) -> None:
        widget.setAutoFillBackground(True)
        color = self.palette().color(QPalette.ColorRole.Window)
        pal = widget.palette()
        pal.setColor(QPalette.ColorRole.Window, color)
        pal.setColor(QPalette.ColorRole.Base, color)
        widget.setPalette(pal)

    def _wrap_scroll(self, body: QWidget) -> QScrollArea:
        width = self._module_width
        body.setMinimumWidth(0)
        body.setMaximumWidth(width)
        body.setSizePolicy(
            QSizePolicy.Policy.Ignored,
            QSizePolicy.Policy.Minimum,
        )
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setMinimumWidth(0)
        scroll.setFixedWidth(width)
        scroll.setSizePolicy(
            QSizePolicy.Policy.Fixed,
            QSizePolicy.Policy.Expanding,
        )
        scroll.setWidget(body)
        self._fill_rail_surface(scroll)
        self._fill_rail_surface(scroll.viewport())
        self._fill_rail_surface(body)
        scroll.viewport().setMinimumWidth(0)
        return scroll

    def _form_layout(self) -> QFormLayout:
        form = QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        form.setHorizontalSpacing(8)
        form.setVerticalSpacing(6)
        form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        return form

    def _expanding_spin(self, spin: QWidget) -> QWidget:
        spin.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        spin.setMinimumWidth(0)
        return spin

    def _action_button(self, text: str, parent: QWidget) -> QPushButton:
        """Select / Fit / Export — same height as MSA overview toolbar buttons."""
        btn = QPushButton(text, parent)
        btn.setFixedHeight(OVERVIEW_SELECT_FIT_HEIGHT)
        btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        btn.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        btn.setMinimumWidth(0)
        return btn

    def _build_display_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)
        layout.addWidget(_page_title_with_rule("Display", page))

        spin_expand = QSizePolicy.Policy.Expanding

        # --- Zoom: − / + / slider (no section header) ---------------------
        self._view_host = QWidget(page)
        view_col = QVBoxLayout(self._view_host)
        view_col.setContentsMargins(0, 0, 0, 0)
        view_col.setSpacing(6)

        zoom_row = QHBoxLayout()
        zoom_row.setContentsMargins(0, 0, 0, 0)
        zoom_row.setSpacing(6)
        self.zoom_out_btn = make_overview_zoom_button("−", self._view_host)
        self.zoom_out_btn.setToolTip(
            "Zoom out (Shift+scroll, Shift+drag, or pinch)"
        )
        self.zoom_in_btn = make_overview_zoom_button("+", self._view_host)
        self.zoom_in_btn.setToolTip(
            "Zoom in (Shift+scroll, Shift+drag, or pinch)"
        )
        self.zoom_slider = NoWheelSlider(
            Qt.Orientation.Horizontal, self._view_host
        )
        self.zoom_slider.setObjectName("annotationZoomSlider")
        self.zoom_slider.setRange(0, 100)
        self.zoom_slider.setValue(0)
        self.zoom_slider.setMinimumWidth(72)
        self.zoom_slider.setFixedHeight(ZOOM_SLIDER_HEIGHT)
        self.zoom_slider.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self.zoom_slider.setToolTip(
            "Click the bar: left = fit, right = maximum zoom. "
            "Drag to adjust."
        )
        apply_annotation_zoom_slider_style(self.zoom_slider)
        self.zoom_slider.valueChanged.connect(self._on_zoom_slider)
        zoom_row.addWidget(self.zoom_out_btn, 0)
        zoom_row.addWidget(self.zoom_in_btn, 0)
        zoom_row.addWidget(self.zoom_slider, 1)
        view_col.addLayout(zoom_row)

        action_row = QHBoxLayout()
        action_row.setContentsMargins(0, 0, 0, 0)
        action_row.setSpacing(6)
        self.select_btn = self._action_button("Select", self._view_host)
        self.select_btn.setToolTip(
            "Zoom to the selected region (drag a rectangle on the plot first)"
        )
        self.fit_btn = self._action_button("Fit", self._view_host)
        self.fit_btn.setToolTip("Fit the whole plot in the view")
        self.export_btn = self._action_button("Export", self._view_host)
        self.export_btn.setToolTip(
            "Export the current view as PDF, PNG, or SVG"
        )
        action_row.addWidget(self.select_btn, 1)
        action_row.addWidget(self.fit_btn, 1)
        action_row.addWidget(self.export_btn, 1)
        view_col.addLayout(action_row)
        layout.addWidget(self._view_host)

        # --- Sequence / parameter / appearance controls (flat) ------------
        self._params_divider = _h_divider(page)
        layout.addWidget(self._params_divider)

        self._params_host = QWidget(page)
        params_col = QVBoxLayout(self._params_host)
        params_col.setContentsMargins(0, 0, 0, 0)
        params_col.setSpacing(6)
        form = self._form_layout()

        self.pair_label = QLabel("Compare:", self._params_host)
        self.pair_combo = QComboBox(self._params_host)
        self.pair_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.pair_combo.setMinimumContentsLength(8)
        self.pair_combo.setSizePolicy(spin_expand, QSizePolicy.Policy.Fixed)
        self.pair_combo.setToolTip(
            "Which sequence pair to plot. 'All sequences' lays every query and "
            "target end to end, with separators between them."
        )
        form.addRow(self.pair_label, self.pair_combo)

        self.sequence_label = QLabel("Sequence:", self._params_host)
        self.sequence_combo = QComboBox(self._params_host)
        self.sequence_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.sequence_combo.setMinimumContentsLength(8)
        self.sequence_combo.setSizePolicy(spin_expand, QSizePolicy.Policy.Fixed)
        self.sequence_combo.setToolTip(
            "Which sequence to show as a self-dotplot."
        )
        self.sequence_prev_btn = QToolButton(self._params_host)
        self.sequence_prev_btn.setText("◀")
        self.sequence_prev_btn.setAutoRaise(True)
        self.sequence_prev_btn.setFixedSize(22, 22)
        self.sequence_prev_btn.setToolTip("Previous sequence")
        self.sequence_prev_btn.clicked.connect(
            lambda: self.step_sequence(-1)
        )
        self.sequence_next_btn = QToolButton(self._params_host)
        self.sequence_next_btn.setText("▶")
        self.sequence_next_btn.setAutoRaise(True)
        self.sequence_next_btn.setFixedSize(22, 22)
        self.sequence_next_btn.setToolTip("Next sequence")
        self.sequence_next_btn.clicked.connect(
            lambda: self.step_sequence(1)
        )
        sequence_row = QWidget(self._params_host)
        sequence_row.setMinimumWidth(0)
        sequence_layout = QHBoxLayout(sequence_row)
        sequence_layout.setContentsMargins(0, 0, 0, 0)
        sequence_layout.setSpacing(2)
        sequence_layout.addWidget(self.sequence_combo, 1)
        sequence_layout.addWidget(self.sequence_prev_btn, 0)
        sequence_layout.addWidget(self.sequence_next_btn, 0)
        self.sequence_combo.currentIndexChanged.connect(
            self.sync_sequence_step_buttons
        )
        form.addRow(self.sequence_label, sequence_row)
        self.sync_sequence_step_buttons()

        # Always directly under Sequence (Grid and Single).
        self.word_label = QLabel("Word size:", self._params_host)
        self.word = QSpinBox(self._params_host)
        self.word.setRange(MIN_WORD_SIZE, MAX_WORD_SIZE)
        self.word.setValue(DEFAULT_WORD_SIZE)
        self._expanding_spin(self.word)
        self.word.setToolTip(
            "Length of the exact match used to place a dot. Larger = stricter "
            "and faster; smaller = more sensitive and noisier."
        )
        form.addRow(self.word_label, self.word)

        self.line_width_label = QLabel("Line width:", self._params_host)
        self.line_width = QSpinBox(self._params_host)
        self.line_width.setRange(_MIN_LINE_WIDTH_STEP, _MAX_LINE_WIDTH_STEP)
        self.line_width.setValue(_DEFAULT_LINE_WIDTH_STEP)
        self._expanding_spin(self.line_width)
        self.line_width.setToolTip(
            "Match and self-diagonal thickness for Grid thumbnails and "
            "the Single plot.\n"
            "Steps 1–15 use a fine gradient (1–4.5 screen pixels)."
        )
        form.addRow(self.line_width_label, self.line_width)

        # Multi-sequence: Grid (default) vs Single (one plot).
        self.view_mode_label = QLabel("View:", self._params_host)
        self.view_mode_combo = QComboBox(self._params_host)
        self.view_mode_combo.addItem("Grid", "grid")
        self.view_mode_combo.addItem("Single", "single")
        self.view_mode_combo.setToolTip(
            "Grid: many sequences at a glance (paged; default for MSA / multi-FASTA).\n"
            "Single: one self-dotplot at a time (◀/▶ to step)."
        )
        self.view_mode_combo.currentIndexChanged.connect(
            self._on_view_mode_changed
        )
        form.addRow(self.view_mode_label, self.view_mode_combo)

        self.page_size_label = QLabel("Per page:", self._params_host)
        self.page_size = QSpinBox(self._params_host)
        self.page_size.setRange(4, 100)
        self.page_size.setValue(20)
        self.page_size.setToolTip(
            "How many self-dotplots to show per Grid page (default 20, max 100)."
        )
        self._expanding_spin(self.page_size)
        form.addRow(self.page_size_label, self.page_size)

        self.thumb_size_label = QLabel("Plot size:", self._params_host)
        self.thumb_size = QSpinBox(self._params_host)
        self.thumb_size.setRange(80, 280)
        self.thumb_size.setSingleStep(10)
        self.thumb_size.setValue(140)
        self.thumb_size.setSuffix(" px")
        self.thumb_size.setToolTip(
            "Side length of each Grid plot (width and height). "
            "The grid fills the window width with square tiles."
        )
        self._expanding_spin(self.thumb_size)
        form.addRow(self.thumb_size_label, self.thumb_size)

        # Page picker mirrors Sequence: dropdown + ◀/▶.
        self.gallery_page_label = QLabel("Page:", self._params_host)
        self.gallery_page_combo = QComboBox(self._params_host)
        self.gallery_page_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self.gallery_page_combo.setMinimumContentsLength(8)
        self.gallery_page_combo.setSizePolicy(spin_expand, QSizePolicy.Policy.Fixed)
        self.gallery_page_combo.setToolTip(
            "Which Grid page (sequence range) to display."
        )
        self.gallery_page_combo.currentIndexChanged.connect(
            self._on_gallery_page_combo_changed
        )
        self.gallery_prev_btn = QToolButton(self._params_host)
        self.gallery_prev_btn.setText("◀")
        self.gallery_prev_btn.setAutoRaise(True)
        self.gallery_prev_btn.setFixedSize(22, 22)
        self.gallery_prev_btn.setToolTip("Previous Grid page")
        self.gallery_prev_btn.clicked.connect(
            lambda: self.galleryPageStep.emit(-1)
        )
        self.gallery_next_btn = QToolButton(self._params_host)
        self.gallery_next_btn.setText("▶")
        self.gallery_next_btn.setAutoRaise(True)
        self.gallery_next_btn.setFixedSize(22, 22)
        self.gallery_next_btn.setToolTip("Next Grid page")
        self.gallery_next_btn.clicked.connect(
            lambda: self.galleryPageStep.emit(1)
        )
        gallery_page_row = QWidget(self._params_host)
        gallery_page_row.setMinimumWidth(0)
        gallery_page_layout = QHBoxLayout(gallery_page_row)
        gallery_page_layout.setContentsMargins(0, 0, 0, 0)
        gallery_page_layout.setSpacing(2)
        gallery_page_layout.addWidget(self.gallery_page_combo, 1)
        gallery_page_layout.addWidget(self.gallery_prev_btn, 0)
        gallery_page_layout.addWidget(self.gallery_next_btn, 0)
        form.addRow(self.gallery_page_label, gallery_page_row)
        self.page_size.valueChanged.connect(self.galleryPageSizeChanged.emit)
        self.thumb_size.valueChanged.connect(self.galleryThumbSizeChanged.emit)
        self._updating_gallery_page_combo = False

        # Created here so visibility helpers can run; placed in the form after
        # "Show reverse" below.
        self.return_grid_btn = QPushButton("◀ Back to Grid", self._params_host)
        self.return_grid_btn.setToolTip(
            "Return to the Grid view at the same page / sequence."
        )
        self.return_grid_btn.clicked.connect(self.returnToGridClicked.emit)
        self.return_grid_btn.hide()

        self.set_view_mode_visible(False)
        self.set_gallery_nav_visible(False)
        self.set_return_to_grid_visible(False)

        self.filter_len_label = QLabel("Min length:", self._params_host)
        self.filter_len = QSpinBox(self._params_host)
        self.filter_len.setRange(0, 10_000_000)
        self.filter_len.setValue(0)
        self.filter_len.setSingleStep(100)
        self.filter_len.setGroupSeparatorShown(True)
        self._expanding_spin(self.filter_len)
        self.filter_len.setToolTip(
            "Hide alignments shorter than this, in bp.\n"
            "0 shows everything. Raise it to strip short, noisy matches the "
            "way a larger word size would when computing a dotplot."
        )
        form.addRow(self.filter_len_label, self.filter_len)

        self.region_start_label = QLabel("From:", self._params_host)
        self.region_start = QSpinBox(self._params_host)
        self.region_start.setRange(1, 2_000_000_000)
        self.region_start.setValue(1)
        self.region_start.setGroupSeparatorShown(True)
        self._expanding_spin(self.region_start)
        self.region_start.setToolTip("First base of the window to plot.")
        form.addRow(self.region_start_label, self.region_start)

        self.region_size_label = QLabel("Window:", self._params_host)
        self.region_size = QSpinBox(self._params_host)
        self.region_size.setRange(1_000, 5_000_000)
        self.region_size.setValue(200_000)
        self.region_size.setSingleStep(50_000)
        self.region_size.setGroupSeparatorShown(True)
        self._expanding_spin(self.region_size)
        self.region_size.setToolTip(
            "How many bases to plot from the start position.\n"
            "A whole scaffold cannot be plotted at small word sizes -- the "
            "match count grows with the square of the repeat copy number -- "
            "but a window of a few hundred kb is instant."
        )
        form.addRow(self.region_size_label, self.region_size)

        self.reverse = QCheckBox("Show reverse", self._params_host)
        self.reverse.setChecked(True)
        self.reverse.setToolTip(
            "Show reverse-strand hits (inverted repeats, TIRs, hairpins) "
            "as anti-diagonal lines."
        )

        self.identity_label = QLabel("Min identity:", self._params_host)
        self.identity = QDoubleSpinBox(self._params_host)
        self.identity.setRange(0.5, 1.0)
        self.identity.setSingleStep(0.05)
        self.identity.setDecimals(2)
        self.identity.setValue(DEFAULT_MIN_IDENTITY)
        self._expanding_spin(self.identity)
        self.identity.setToolTip("Minimum alignment identity reported by FastGA.")
        form.addRow(self.identity_label, self.identity)

        self.min_len_label = QLabel("Min length:", self._params_host)
        self.min_len = QSpinBox(self._params_host)
        self.min_len.setRange(100, 100_000)
        self.min_len.setSingleStep(50)
        self.min_len.setValue(DEFAULT_MIN_LENGTH)
        self._expanding_spin(self.min_len)
        self.min_len.setToolTip(
            "Shortest alignment FastGA will report, in bp.\n"
            "Cannot go below 100: alignments must contain a seed, and FastGA "
            "rejects shorter limits.\n"
            "Raise it on satellite arrays, where short matches swamp the plot."
        )
        form.addRow(self.min_len_label, self.min_len)

        self.seed_freq_label = QLabel("Repeat cutoff:", self._params_host)
        self.seed_freq = QSpinBox(self._params_host)
        self.seed_freq.setRange(1, MAX_SEED_FREQ)
        self.seed_freq.setValue(DEFAULT_SEED_FREQ)
        self.seed_freq.setSingleStep(50)
        self._expanding_spin(self.seed_freq)
        self.seed_freq.setToolTip(
            "FastGA's -f: seeds occurring more often than this are discarded as "
            "'repetitive' before alignment.\n"
            "FastGA's own default is 10, which removes satellite arrays entirely "
            "-- every monomer k-mer occurs far more than 10 times.\n"
            "Raise it to see repeat structure; higher values are slower."
        )
        form.addRow(self.seed_freq_label, self.seed_freq)

        # Appear last among the common self-dotplot controls.
        form.addRow(self.reverse)
        # Directly under "Show reverse" (Single mode after opening a thumb).
        form.addRow(self.return_grid_btn)

        params_col.addLayout(form)

        self.compute_btn = QPushButton("Run FastGA", self._params_host)
        self.compute_btn.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self.compute_btn.setToolTip(
            "Run FastGA to build a self-dotplot of this long sequence. "
            "This can take several minutes on large chromosomes."
        )
        self.compute_btn.hide()
        self.compute_btn.clicked.connect(self.computeClicked.emit)
        params_col.addWidget(self.compute_btn)

        self.cancel_btn = QPushButton("Cancel", self._params_host)
        self.cancel_btn.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self.cancel_btn.hide()
        self.cancel_btn.clicked.connect(self.cancelClicked.emit)
        params_col.addWidget(self.cancel_btn)

        layout.addWidget(self._params_host)

        # Selection ranges at the bottom: axis label + aligned range box + View.
        self.sel_row = QWidget(page)
        sel_col = QVBoxLayout(self.sel_row)
        sel_col.setContentsMargins(0, 4, 0, 0)
        sel_col.setSpacing(4)

        x_row = QHBoxLayout()
        x_row.setContentsMargins(0, 0, 0, 0)
        x_row.setSpacing(6)
        axis_fm = QFontMetrics(self.sel_row.font())
        axis_w = max(axis_fm.horizontalAdvance("X"), axis_fm.horizontalAdvance("Y"))
        self.sel_x_label = QLabel("X", self.sel_row)
        self.sel_x_label.setFixedWidth(axis_w)
        self.sel_x_label.setFixedHeight(OVERVIEW_SELECT_FIT_HEIGHT)
        self.sel_x_label.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        )
        self.sel_x = _AlignedRangeField(self.sel_row)
        self.view_x_btn = self._action_button("View", self.sel_row)
        self.view_x_btn.setToolTip(
            "Open Annotation and highlight the selected X (horizontal) range"
        )
        x_row.addWidget(self.sel_x_label, 0)
        x_row.addWidget(self.sel_x, 1)
        x_row.addWidget(self.view_x_btn, 0)
        sel_col.addLayout(x_row)

        y_row = QHBoxLayout()
        y_row.setContentsMargins(0, 0, 0, 0)
        y_row.setSpacing(6)
        self.sel_y_label = QLabel("Y", self.sel_row)
        self.sel_y_label.setFixedWidth(axis_w)
        self.sel_y_label.setFixedHeight(OVERVIEW_SELECT_FIT_HEIGHT)
        self.sel_y_label.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        )
        self.sel_y = _AlignedRangeField(self.sel_row)
        self.view_y_btn = self._action_button("View", self.sel_row)
        self.view_y_btn.setToolTip(
            "Open Annotation and highlight the selected Y (vertical) range"
        )
        y_row.addWidget(self.sel_y_label, 0)
        y_row.addWidget(self.sel_y, 1)
        y_row.addWidget(self.view_y_btn, 0)
        sel_col.addLayout(y_row)

        self.view_x_btn.setEnabled(False)
        self.view_y_btn.setEnabled(False)
        layout.addWidget(self.sel_row)
        layout.addStretch(1)
        return self._wrap_scroll(page)

    def _on_zoom_slider(self, value: int) -> None:
        if self._updating_zoom_slider:
            return
        self.zoomSliderChanged.emit(int(value))

    def set_zoom_slider_value(self, value: int) -> None:
        value = max(0, min(100, int(value)))
        if self.zoom_slider.value() == value:
            return
        self._updating_zoom_slider = True
        try:
            self.zoom_slider.setValue(value)
        finally:
            self._updating_zoom_slider = False

    def set_selection_ranges(
        self,
        *,
        x: tuple[int, int] | None = None,
        y: tuple[int, int] | None = None,
    ) -> None:
        """Update x/y range boxes; enable View when a range is present.

        Start numbers are right-aligned in a shared column so the `` - ``
        separators on both rows line up vertically.
        """
        if x is None:
            self.sel_x.clear()
        else:
            self.sel_x.set_range(int(x[0]), int(x[1]))
        if y is None:
            self.sel_y.clear()
        else:
            self.sel_y.set_range(int(y[0]), int(y[1]))

        fm = QFontMetrics(self.sel_x.font())
        start_w = 0
        for field in (self.sel_x, self.sel_y):
            text = field.start_text()
            if text:
                start_w = max(start_w, fm.horizontalAdvance(text))
        self.sel_x.set_start_width(start_w)
        self.sel_y.set_start_width(start_w)

        self.view_x_btn.setEnabled(x is not None)
        self.view_y_btn.setEnabled(y is not None)

    def step_sequence(self, delta: int) -> None:
        """Move the Sequence combo by ``delta`` entries (wraps at the ends)."""
        count = int(self.sequence_combo.count())
        if count <= 0:
            return
        current = int(self.sequence_combo.currentIndex())
        if current < 0:
            current = 0
        next_index = (current + int(delta)) % count
        if next_index == current:
            return
        self.sequence_combo.setCurrentIndex(next_index)

    def sync_sequence_step_buttons(self, *_args) -> None:
        count = int(self.sequence_combo.count())
        enabled = count > 1
        self.sequence_prev_btn.setEnabled(enabled)
        self.sequence_next_btn.setEnabled(enabled)

    def view_mode(self) -> str:
        data = self.view_mode_combo.currentData()
        return str(data or "single")

    def set_view_mode(self, mode: str) -> None:
        want = "grid" if str(mode).strip().lower() == "grid" else "single"
        idx = self.view_mode_combo.findData(want)
        if idx < 0 or self.view_mode_combo.currentIndex() == idx:
            return
        self.view_mode_combo.blockSignals(True)
        self.view_mode_combo.setCurrentIndex(idx)
        self.view_mode_combo.blockSignals(False)

    def set_view_mode_visible(self, visible: bool) -> None:
        visible = bool(visible)
        self.view_mode_label.setVisible(visible)
        self.view_mode_combo.setVisible(visible)
        if not visible:
            self.set_gallery_nav_visible(False)
            self.set_return_to_grid_visible(False)

    def set_gallery_nav_visible(self, visible: bool) -> None:
        visible = bool(visible)
        self.page_size_label.setVisible(visible)
        self.page_size.setVisible(visible)
        self.thumb_size_label.setVisible(visible)
        self.thumb_size.setVisible(visible)
        self.gallery_page_label.setVisible(visible)
        self.gallery_page_combo.setVisible(visible)
        self.gallery_prev_btn.setVisible(visible)
        self.gallery_next_btn.setVisible(visible)

    def set_return_to_grid_visible(self, visible: bool) -> None:
        self.return_grid_btn.setVisible(bool(visible))

    def set_gallery_page_choices(
        self,
        *,
        page: int,
        page_size: int,
        total: int,
    ) -> None:
        """Fill the Page dropdown like Sequence (ranges + ◀/▶)."""
        page = max(0, int(page))
        page_size = max(1, int(page_size))
        total = max(0, int(total))
        page_count = max(1, (total + page_size - 1) // page_size) if total else 1
        page = min(page, page_count - 1)
        self._updating_gallery_page_combo = True
        try:
            self.gallery_page_combo.clear()
            if total <= 0:
                self.gallery_page_combo.addItem("—", -1)
                self.gallery_prev_btn.setEnabled(False)
                self.gallery_next_btn.setEnabled(False)
                return
            for i in range(page_count):
                start = i * page_size + 1
                end = min(total, (i + 1) * page_size)
                self.gallery_page_combo.addItem(
                    f"{start}–{end} of {total}", i
                )
            self.gallery_page_combo.setCurrentIndex(page)
            self.gallery_prev_btn.setEnabled(page > 0)
            self.gallery_next_btn.setEnabled(page + 1 < page_count)
        finally:
            self._updating_gallery_page_combo = False

    def set_gallery_page_text(self, text: str) -> None:
        # Kept for callers; Page dropdown now carries the range label.
        del text

    def _on_gallery_page_combo_changed(self, index: int) -> None:
        if self._updating_gallery_page_combo:
            return
        if index < 0:
            return
        page = self.gallery_page_combo.itemData(index)
        if page is None:
            return
        self.galleryPageSelected.emit(int(page))

    def _on_view_mode_changed(self, _index: int = -1) -> None:
        mode = self.view_mode()
        self.set_gallery_nav_visible(mode == "grid")
        self.viewModeChanged.emit(mode)

    def sync_control_visibility(
        self,
        *,
        view: bool = True,
        params: bool = True,
    ) -> None:
        """Show or hide the zoom block vs. parameter fields as a group."""
        self._view_host.setVisible(view)
        self._params_host.setVisible(params)
        self._params_divider.setVisible(view and params)

    def is_panel_collapsed(self) -> bool:
        return bool(self._panel_collapsed)

    def set_panel_collapsed(self, collapsed: bool, *, persist: bool = True) -> None:
        collapsed = bool(collapsed)
        if collapsed == self._panel_collapsed:
            if collapsed:
                self._sync_tab_checked_state(None)
            return
        self._panel_collapsed = collapsed
        self._stack.setVisible(not collapsed)
        self._divider.setVisible(not collapsed)
        if collapsed:
            self.setFixedWidth(_COLLAPSED_RAIL_WIDTH)
            self._sync_tab_checked_state(None)
        else:
            self.setFixedWidth(self._expanded_rail_width())
            self._stack.setCurrentIndex(self._active_tab)
            self._sync_tab_checked_state(self._active_tab)
        self.panelCollapsedChanged.emit(collapsed)
        if persist:
            self._save_rail_ui_state()

    def _sync_tab_checked_state(self, active: int | None) -> None:
        self._updating_tabs = True
        try:
            for i, btn in enumerate(self._tab_buttons):
                btn.setChecked(active is not None and i == active)
                btn.update()
        finally:
            self._updating_tabs = False

    def _on_tab_clicked(self, index: int) -> None:
        if self._updating_tabs:
            return
        index = max(0, min(index, self._stack.count() - 1))
        if not self._panel_collapsed and index == self._active_tab:
            self.set_panel_collapsed(True)
            return
        self._active_tab = index
        was_collapsed = self._panel_collapsed
        if was_collapsed:
            self._panel_collapsed = False
            self._stack.setVisible(True)
            self._divider.setVisible(True)
            self.setFixedWidth(self._expanded_rail_width())
        self._stack.setCurrentIndex(index)
        self._sync_tab_checked_state(index)
        if was_collapsed:
            self.panelCollapsedChanged.emit(False)
        self._save_rail_ui_state()

    def _settings(self) -> QSettings:
        settings = QSettings()
        settings.beginGroup(_SETTINGS_GROUP)
        return settings

    def _restore_rail_ui_state(self) -> None:
        settings = self._settings()
        try:
            collapsed_raw = settings.value(_SETTINGS_COLLAPSED_KEY, False)
            if isinstance(collapsed_raw, str):
                collapsed = collapsed_raw.strip().lower() in {"1", "true", "yes"}
            else:
                collapsed = bool(collapsed_raw)
            self._panel_collapsed = False
            self.set_panel_collapsed(collapsed, persist=False)
            if not collapsed:
                self._sync_tab_checked_state(self.TAB_DISPLAY)
        finally:
            settings.endGroup()

    def _save_rail_ui_state(self, *_args) -> None:
        settings = self._settings()
        try:
            settings.setValue(_SETTINGS_COLLAPSED_KEY, bool(self._panel_collapsed))
        finally:
            settings.endGroup()
