"""Geneious-style right rail for the Annotation preview (Display / Annotation / Find / Others)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QPoint, QRect, QRectF, QSettings, QSize, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPen
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLayout,
    QLayoutItem,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QStackedWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox

from repeatviewer.backends.annotation_colors import load_add_feature_buttons
from repeatviewer.backends.genetic_codes import (
    genetic_code_choices,
    normalize_genetic_code_id,
)
from repeatviewer.backends.msa_plot import (
    AMINO_ACID_COLOR_PRESETS,
    NUCLEOTIDE_COLOR_PRESETS,
)
from repeatviewer.widgets.help_button import field_with_help, make_help_button
from repeatviewer.widgets.preview_zoom import (
    load_color_blind_nucleotide_scheme,
    save_color_blind_nucleotide_scheme,
    save_custom_amino_acid_backgrounds,
    save_custom_nucleotide_backgrounds,
)
from repeatviewer.widgets.path_browse import PathBrowseRow
from repeatviewer.widgets.preview_controls import apply_annotation_zoom_slider_style
from repeatviewer.widgets.spin_boxes import (
    NoWheelDoubleSpinBox,
    NoWheelSlider,
    NoWheelSpinBox,
)

_NCBI_GENETIC_CODES_URL = (
    "https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi"
)
_READING_FRAMES_HELP = (
    "Tick which frames to show under the sequence when translation is on.<br><br>"
    "<b>Frame 1 / 2 / 3</b> — forward reading frames on the displayed strand "
    "(offset by 0, 1, or 2 bases).<br><br>"
    "<b>Reverse 1 / 2 / 3</b> — reverse reading frames on the complementary strand "
    "(reverse-complement translation)."
)
_FIND_IUPAC_HELP = (
    "<b>On</b>: IUPAC ambiguity codes are treated as wildcards.<br><br>"
    "<b>Query</b> — e.g. R matches A or G, Y matches C or T, N matches any base.<br>"
    "<b>Sequence</b> — ambiguous letters in the sequence (N, R, Y, …) also count "
    "as a match for the corresponding motif base.<br><br>"
    "<b>Off</b>: each letter matches only itself (an R in the motif only matches "
    "an R in the sequence)."
)
_GENETIC_CODE_HELP = (
    "The genetic code table controls how DNA codons are translated into amino acids "
    "(same IDs as GenBank <code>/transl_table</code>). Use <b>1. Standard</b> for most "
    "nuclear genomes; choose a mitochondrial, plastid, or other table when your sequence "
    "comes from those systems.<br><br>"
    "For more detailed explanation please visit:<br>"
    f'<a href="{_NCBI_GENETIC_CODES_URL}">{_NCBI_GENETIC_CODES_URL}</a>'
)

_RAIL_TAB_WIDTH = 24
_MODULE_WIDTH = 260
_TAB_BUTTON_HEIGHT_MIN = 36
_TAB_PILL_RADIUS = 7
_TAB_TEXT_PAD = 4  # vertical padding around rotated label

# Persist Annotation rail tab + collapsible section expand state.
_SETTINGS_GROUP = "sequence_annotation_rail"
_SETTINGS_TAB_KEY = "last_tab"
_SETTINGS_SECTIONS_KEY = "expanded_sections"
_SETTINGS_COLLAPSED_KEY = "panel_collapsed"
# Default first open: Annotation tab with all sections collapsed.
_DEFAULT_TAB_INDEX = 1
_COLLAPSED_RAIL_WIDTH = _RAIL_TAB_WIDTH
_ANNOTATION_SECTION_IDS = (
    "manual",
    "feature_detection",
    "similarity",
    "transfer",
)

# Translation frames: +1..+3 forward, -1..-3 reverse.
TRANSLATION_FRAMES: tuple[int, ...] = (1, 2, 3, -1, -2, -3)


class TranslationSettingsDialog(QDialog):
    """Frames, three-letter amino acids, and NCBI genetic code for translation tracks."""

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        frames: frozenset[int] | set[int] | None = None,
        three_letter: bool = False,
        genetic_code: int = 1,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Translation settings")
        self.setModal(True)
        self.setMinimumWidth(400)

        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(10)

        frames_section = QVBoxLayout()
        frames_section.setContentsMargins(0, 0, 0, 0)
        frames_section.setSpacing(4)
        frames_header = QHBoxLayout()
        frames_header.setContentsMargins(0, 0, 0, 0)
        frames_header.setSpacing(6)
        frames_header.addWidget(_section_label("Reading frames", self))
        frames_header.addWidget(
            make_help_button(self, "Reading frames", _READING_FRAMES_HELP)
        )
        frames_header.addStretch(1)
        frames_section.addLayout(frames_header)

        frames_host = QWidget(self)
        frames_grid = QGridLayout(frames_host)
        frames_grid.setContentsMargins(0, 0, 0, 0)
        frames_grid.setHorizontalSpacing(12)
        frames_grid.setVerticalSpacing(4)
        selected = set(frames or {1})
        self._frame_checks: dict[int, QCheckBox] = {}
        for col, frame in enumerate((1, 2, 3)):
            fwd = QCheckBox(f"Frame {frame}", frames_host)
            fwd.setToolTip(f"Forward reading frame {frame}")
            fwd.setChecked(frame in selected)
            self._frame_checks[frame] = fwd
            frames_grid.addWidget(fwd, 0, col)

            rev_frame = -frame
            rev = QCheckBox(f"Reverse {frame}", frames_host)
            rev.setToolTip(f"Reverse reading frame {frame}")
            rev.setChecked(rev_frame in selected)
            self._frame_checks[rev_frame] = rev
            frames_grid.addWidget(rev, 1, col)
        frames_section.addWidget(frames_host)
        root.addLayout(frames_section)

        root.addWidget(_h_divider(self))

        aa_section = QVBoxLayout()
        aa_section.setContentsMargins(0, 0, 0, 0)
        aa_section.setSpacing(4)
        aa_section.addWidget(_section_label("Amino acids", self))
        self.three_letter = QCheckBox("Show three-letter codes (e.g. Met)", self)
        self.three_letter.setChecked(bool(three_letter))
        self.three_letter.setToolTip(
            "When checked, translation tracks show Ala/Met/Ter instead of A/M/*."
        )
        aa_section.addWidget(self.three_letter)
        root.addLayout(aa_section)

        root.addWidget(_h_divider(self))

        code_section = QVBoxLayout()
        code_section.setContentsMargins(0, 0, 0, 0)
        code_section.setSpacing(4)
        code_section.addWidget(_section_label("Genetic code", self))
        code_row = QHBoxLayout()
        code_row.setContentsMargins(0, 0, 0, 0)
        code_row.setSpacing(6)
        self.genetic_code = QComboBox(self)
        self.genetic_code.setSizePolicy(
            QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed
        )
        self.genetic_code.setFixedWidth(220)
        for table_id, label in genetic_code_choices():
            self.genetic_code.addItem(label, table_id)
        wanted = normalize_genetic_code_id(genetic_code)
        idx = self.genetic_code.findData(wanted)
        self.genetic_code.setCurrentIndex(max(0, idx))
        code_row.addWidget(self.genetic_code)
        code_row.addWidget(
            make_help_button(self, "Genetic code", _GENETIC_CODE_HELP)
        )
        code_row.addStretch(1)
        code_section.addLayout(code_row)
        root.addLayout(code_section)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel,
            self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        root.addWidget(buttons)

    def selected_frames(self) -> frozenset[int]:
        frames = frozenset(
            frame for frame, check in self._frame_checks.items() if check.isChecked()
        )
        return frames or frozenset({1})

    def show_three_letter(self) -> bool:
        return bool(self.three_letter.isChecked())

    def genetic_code_id(self) -> int:
        return normalize_genetic_code_id(self.genetic_code.currentData() or 1)


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


class _FlowLayout(QLayout):
    """Left-to-right wrapping row layout for compact checkbox chips."""

    def __init__(self, parent: QWidget | None = None, *, spacing: int = 6) -> None:
        super().__init__(parent)
        self._items: list[QLayoutItem] = []
        self.setContentsMargins(0, 0, 0, 0)
        self.setSpacing(spacing)

    def addItem(self, item: QLayoutItem) -> None:  # noqa: N802
        self._items.append(item)

    def count(self) -> int:
        return len(self._items)

    def itemAt(self, index: int) -> QLayoutItem | None:  # noqa: N802
        if 0 <= index < len(self._items):
            return self._items[index]
        return None

    def takeAt(self, index: int) -> QLayoutItem | None:  # noqa: N802
        if 0 <= index < len(self._items):
            return self._items.pop(index)
        return None

    def expandingDirections(self) -> Qt.Orientation:  # noqa: N802
        return Qt.Orientation(0)

    def hasHeightForWidth(self) -> bool:  # noqa: N802
        return True

    def heightForWidth(self, width: int) -> int:  # noqa: N802
        return self._do_layout(QRect(0, 0, width, 0), test_only=True)

    def setGeometry(self, rect: QRect) -> None:  # noqa: N802
        super().setGeometry(rect)
        self._do_layout(rect, test_only=False)

    def sizeHint(self) -> QSize:  # noqa: N802
        return self.minimumSize()

    def minimumSize(self) -> QSize:  # noqa: N802
        size = QSize()
        for item in self._items:
            size = size.expandedTo(item.minimumSize())
        left, top, right, bottom = self.getContentsMargins()
        size += QSize(left + right, top + bottom)
        return size

    def _do_layout(self, rect: QRect, *, test_only: bool) -> int:
        left, top, right, bottom = self.getContentsMargins()
        effective = rect.adjusted(left, top, -right, -bottom)
        x = effective.x()
        y = effective.y()
        line_height = 0
        space = self.spacing()
        for item in self._items:
            hint = item.sizeHint()
            next_x = x + hint.width() + space
            if line_height and next_x - space > effective.right() + 1:
                x = effective.x()
                y += line_height + space
                next_x = x + hint.width() + space
                line_height = 0
            if not test_only:
                item.setGeometry(QRect(QPoint(x, y), hint))
            x = next_x
            line_height = max(line_height, hint.height())
        return y + line_height - rect.y() + bottom


class _CollapsibleSection(QWidget):
    """Section title matching Annotations headers; click (or right arrow) to expand."""

    expandedChanged = Signal(bool)

    def __init__(
        self,
        title: str,
        parent: QWidget | None = None,
        *,
        expanded: bool = True,
    ) -> None:
        super().__init__(parent)
        self.setMinimumWidth(0)
        self._expanded = bool(expanded)
        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(6)

        self._header = QWidget(self)
        self._header.setObjectName("annotationSectionHeader")
        self._header.setCursor(Qt.CursorShape.PointingHandCursor)
        self._header.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        header_row = QHBoxLayout(self._header)
        header_row.setContentsMargins(0, 2, 0, 2)
        header_row.setSpacing(6)

        self._title = QLabel(title, self._header)
        title_font = self._title.font()
        title_font.setBold(True)
        self._title.setFont(title_font)
        self._title.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred
        )
        header_row.addWidget(self._title, 1)

        self._arrow = QLabel(self._header)
        self._arrow.setFixedWidth(12)
        self._arrow.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        arrow_font = self._arrow.font()
        # 3/4 of the title size (half was too small; 3/2 of that).
        arrow_font.setPointSizeF(max(6.0, title_font.pointSizeF() * 0.75))
        self._arrow.setFont(arrow_font)
        header_row.addWidget(self._arrow, 0)
        self._sync_arrow()
        self._header.mousePressEvent = self._on_header_pressed  # type: ignore[method-assign]
        root.addWidget(self._header)

        self._body = QWidget(self)
        self._body_layout = QVBoxLayout(self._body)
        self._body_layout.setContentsMargins(0, 0, 0, 0)
        self._body_layout.setSpacing(6)
        self._body.setVisible(self._expanded)
        root.addWidget(self._body)

    def body_layout(self) -> QVBoxLayout:
        return self._body_layout

    def add_widget(self, widget: QWidget) -> None:
        self._body_layout.addWidget(widget)

    def add_layout(self, layout: QHBoxLayout | QVBoxLayout | QFormLayout | QGridLayout) -> None:
        self._body_layout.addLayout(layout)

    def is_expanded(self) -> bool:
        return self._expanded

    def set_expanded(self, expanded: bool, *, notify: bool = True) -> None:
        expanded = bool(expanded)
        if expanded == self._expanded:
            return
        self._expanded = expanded
        self._body.setVisible(expanded)
        self._sync_arrow()
        if notify:
            self.expandedChanged.emit(expanded)

    def _sync_arrow(self) -> None:
        # Right side disclosure: ▼ expanded, ▶ collapsed.
        self._arrow.setText("▼" if self._expanded else "▶")

    def _on_header_pressed(self, event) -> None:  # noqa: ANN001
        if event.button() == Qt.MouseButton.LeftButton:
            self.set_expanded(not self._expanded)
        event.accept()


class _RotatedTabButton(QToolButton):
    """Vertical pill tab (Actions-style chrome, label rotated −90°)."""

    def __init__(self, label: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._label = label
        self.setObjectName("AnnotationRailTabButton")
        self.setText(label)
        self.setCheckable(True)
        self.setAutoRaise(False)
        self.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        self.setMouseTracking(True)
        self.toggled.connect(lambda _checked: self.update())
        fm = QFontMetrics(self.font())
        # Rotated text height ≈ label width + minimal inset (no tall empty caps).
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
        text_label = self._label
        painter.translate(self.width() / 2.0, self.height() / 2.0)
        painter.rotate(-90)
        text_w = fm.horizontalAdvance(text_label)
        text_h = fm.ascent()
        painter.drawText(
            QPoint(int(-text_w / 2), int(text_h / 2) - 1),
            text_label,
        )
        painter.end()


class SequenceAnnotationRail(QWidget):
    """Vertical pill-tab strip + stacked Display / Annotation / Find / Others panels."""

    addFeatureClicked = Signal()  # Add with type from manual_feature_combo
    addFeatureSettingsClicked = Signal()  # Manage feature types in the dropdown
    deleteClicked = Signal()
    restoreClicked = Signal()
    saveClicked = Signal()
    exportClicked = Signal()
    colorsClicked = Signal()
    residueColorsChanged = Signal()  # quick nucleotide / protein scheme applied
    featureDetectClicked = Signal()
    featureDetectSettingsClicked = Signal()
    similarityAnnotateClicked = Signal()
    similaritySettingsClicked = Signal()
    transferRequested = Signal()
    transferThresholdChanged = Signal()
    transferSourceChanged = Signal()
    manageDatabaseClicked = Signal()
    tabChanged = Signal(int)
    panelCollapsedChanged = Signal(bool)
    findRequested = Signal(int)  # +1 next, -1 previous, 0 find all

    # Display / size signals (live apply).
    showAnnotationTrackChanged = Signal(bool)
    showComplementChanged = Signal(bool)
    showTranslationChanged = Signal(bool)
    translationFramesChanged = Signal(object)  # frozenset[int]
    translationOptionsChanged = Signal()  # three-letter / genetic code
    basesZoomedInChanged = Signal(int)
    basesZoomedOutChanged = Signal(int)
    annotationLabelSizeChanged = Signal(int)
    annotationThicknessChanged = Signal(int)
    annotationDistanceChanged = Signal(int)
    sequenceNameSizeChanged = Signal(int)
    baseNumberSizeChanged = Signal(int)
    splitTracksBySourceChanged = Signal(bool)

    TAB_DISPLAY = 0
    TAB_SHOW = 0  # backward-compatible alias
    TAB_ANNOTATIONS = 1
    TAB_FIND = 2
    TAB_OTHERS = 3

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("sequenceAnnotationRail")
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Expanding)
        self._module_width = _MODULE_WIDTH
        self._overview_mode = False
        self._panel_collapsed = False
        self._active_tab = _DEFAULT_TAB_INDEX
        self._updating_tabs = False
        self.setFixedWidth(self._expanded_rail_width())

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        self._stack = QStackedWidget(self)
        self._stack.setFixedWidth(self._module_width)
        self._stack.setSizePolicy(
            QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Expanding
        )
        # Clip overflowing page content so it cannot paint over the tab strip.
        self._stack.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self._stack.setStyleSheet(
            "QStackedWidget { background: transparent; }"
        )

        self._show_page = self._build_show_page()
        self._annotations_page = self._build_annotations_page()
        self._find_page = self._build_find_page()
        self._others_page = self._build_others_page()
        self._stack.addWidget(self._show_page)
        self._stack.addWidget(self._annotations_page)
        self._stack.addWidget(self._find_page)
        self._stack.addWidget(self._others_page)

        self._divider = QFrame(self)
        self._divider.setFrameShape(QFrame.Shape.VLine)
        self._divider.setFrameShadow(QFrame.Shadow.Plain)
        self._divider.setFixedWidth(1)
        self._divider.setStyleSheet("color: palette(mid);")

        self._tab_strip = QWidget(self)
        self._tab_strip.setFixedWidth(_RAIL_TAB_WIDTH)
        self._tab_strip.setObjectName("annotationRailTabs")
        tab_layout = QVBoxLayout(self._tab_strip)
        tab_layout.setContentsMargins(1, 2, 1, 2)
        tab_layout.setSpacing(2)

        # Non-exclusive: re-clicking the active tab collapses the panel.
        self._tab_group = QButtonGroup(self)
        self._tab_group.setExclusive(False)
        self._tab_buttons: list[_RotatedTabButton] = []
        self._tab_tips = (
            "Display — sequence display and colors (click again to hide)",
            "Annotations — detect, edit, and transfer features (click again to hide)",
            "Find — motif search in the active sequence (click again to hide)",
            "Others — sizes and extra display options (click again to hide)",
        )
        for index, (label, tip) in enumerate(
            (
                ("Display", self._tab_tips[0]),
                ("Annotation", self._tab_tips[1]),
                ("Find", self._tab_tips[2]),
                ("Others", self._tab_tips[3]),
            )
        ):
            btn = _RotatedTabButton(label, self._tab_strip)
            btn.setToolTip(tip)
            self._tab_group.addButton(btn, index)
            tab_layout.addWidget(btn, 0, Qt.AlignmentFlag.AlignHCenter)
            self._tab_buttons.append(btn)
            btn.clicked.connect(lambda _checked=False, i=index: self._on_tab_clicked(i))
        tab_layout.addStretch(1)

        root.addWidget(self._stack, 0)
        root.addWidget(self._divider, 0)
        root.addWidget(self._tab_strip, 0)
        self._restore_rail_ui_state()

    def _expanded_rail_width(self) -> int:
        return int(self._module_width) + _RAIL_TAB_WIDTH + 1

    def _apply_module_geometry(self) -> None:
        """Resize the module column and every scroll page to ``_module_width``."""
        width = max(200, int(self._module_width))
        self._module_width = width
        self._stack.setFixedWidth(width)
        for index in range(self._stack.count()):
            page = self._stack.widget(index)
            if page is None:
                continue
            page.setFixedWidth(width)
            page.setMinimumWidth(0)
            page.setMaximumWidth(width)
            if isinstance(page, QScrollArea):
                body = page.widget()
                if body is not None:
                    body.setMinimumWidth(0)
                    body.setMaximumWidth(width)
        if not self._panel_collapsed:
            self.setFixedWidth(self._expanded_rail_width())

    def _wrap_scroll(self, body: QWidget) -> QScrollArea:
        # Never let page contents force the rail wider than the module column.
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
        # Keep page content as wide as the viewport (minus scrollbar when shown).
        scroll.viewport().setAutoFillBackground(False)
        scroll.viewport().setMinimumWidth(0)
        return scroll

    def _full_width_button(self, text: str, page: QWidget) -> QPushButton:
        btn = QPushButton(text, page)
        btn.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        # macOS QPushButton sizeHints are often wider than half the rail.
        btn.setMinimumWidth(0)
        btn.setMaximumWidth(self._module_width)
        return btn

    def _build_show_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)
        layout.addWidget(_section_label("Display", page))

        # Keep sequence switcher widgets for multi-record wiring, but do not show
        # them on the Display panel (overview / Find cover that job).
        self.sequence_combo = QComboBox(page)
        self.sequence_combo.setEditable(False)
        self.sequence_combo.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self.sequence_combo.hide()
        self.sequence_search = QLineEdit(page)
        self.sequence_search.setPlaceholderText("Filter by name…")
        self.sequence_search.setClearButtonEnabled(True)
        self.sequence_search.hide()
        self.seq_count_label = QLabel("", page)
        self.seq_count_label.setStyleSheet("color: palette(mid);")
        self.seq_count_label.hide()

        self._tracks_divider = _h_divider(page)
        self._tracks_label = _section_label("Tracks", page)
        layout.addWidget(self._tracks_divider)
        layout.addWidget(self._tracks_label)

        self.show_annotation_track = QCheckBox("Show annotation track", page)
        self.show_annotation_track.setChecked(True)
        self.show_annotation_track.setToolTip(
            "Show feature annotation tracks under each sequence "
            "(including MSA residue rows)."
        )
        self.show_annotation_track.toggled.connect(
            self.showAnnotationTrackChanged.emit
        )
        layout.addWidget(self.show_annotation_track)

        self.show_complement = QCheckBox("Show complement", page)
        self.show_complement.setToolTip(
            "Show the complementary strand under the forward sequence "
            "(only when zoomed in to individual nucleotides)."
        )
        self.show_complement.toggled.connect(self.showComplementChanged.emit)
        layout.addWidget(self.show_complement)

        translation_row = QHBoxLayout()
        translation_row.setContentsMargins(0, 0, 0, 0)
        translation_row.setSpacing(8)
        self.show_translation = QCheckBox("Show translation", page)
        self.show_translation.setToolTip(
            "Show amino-acid translation rows under the sequence "
            "(only when zoomed in to individual nucleotides)."
        )
        self.show_translation.toggled.connect(self._on_show_translation_toggled)
        translation_row.addWidget(self.show_translation)
        self.translation_more_btn = QPushButton("Setting", page)
        self.translation_more_btn.setEnabled(False)
        self.translation_more_btn.setSizePolicy(
            QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed
        )
        self.translation_more_btn.setToolTip(
            "Choose reading frames, three-letter amino acids, and genetic code."
        )
        self.translation_more_btn.clicked.connect(self._open_translation_settings)
        translation_row.addWidget(self.translation_more_btn)
        translation_row.addStretch(1)
        layout.addLayout(translation_row)

        # Hidden frame checkboxes keep settings load/save and map wiring stable.
        self.frame_checks: dict[int, QCheckBox] = {}
        for frame in TRANSLATION_FRAMES:
            check = QCheckBox(page)
            check.hide()
            check.setChecked(frame == 1)
            self.frame_checks[frame] = check
        self._translation_three_letter = False
        self._translation_genetic_code = 1

        self.split_tracks_by_source = QCheckBox("Split tracks by source", page)
        self.split_tracks_by_source.setToolTip(
            "When checked, features from different sources use separate visual "
            "tracks even if they share the same track number."
        )
        self.split_tracks_by_source.toggled.connect(
            self.splitTracksBySourceChanged.emit
        )
        layout.addWidget(self.split_tracks_by_source)

        layout.addWidget(_h_divider(page))
        layout.addWidget(_section_label("Color", page))
        color_form = QFormLayout()
        color_form.setContentsMargins(0, 0, 0, 0)
        color_form.setHorizontalSpacing(8)
        color_form.setVerticalSpacing(6)
        color_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        color_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )

        scheme_width = 120
        self.nucleotide_scheme = QComboBox(page)
        self.nucleotide_scheme.setSizePolicy(
            QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed
        )
        self.nucleotide_scheme.setFixedWidth(scheme_width)
        self.nucleotide_scheme.setToolTip(
            "Quick nucleotide letter colours for the sequence view."
        )
        for name, _palette in NUCLEOTIDE_COLOR_PRESETS:
            self.nucleotide_scheme.addItem(name)
        self.nucleotide_scheme.currentIndexChanged.connect(
            self._on_nucleotide_scheme_changed
        )

        self.protein_scheme = QComboBox(page)
        self.protein_scheme.setSizePolicy(
            QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed
        )
        self.protein_scheme.setFixedWidth(scheme_width)
        self.protein_scheme.setToolTip(
            "Quick amino-acid colour scheme for protein sequences / translations."
        )
        for name, _palette in AMINO_ACID_COLOR_PRESETS:
            self.protein_scheme.addItem(name)
        self.protein_scheme.currentIndexChanged.connect(
            self._on_protein_scheme_changed
        )

        color_form.addRow("Nucleotide:", self.nucleotide_scheme)
        color_form.addRow("Protein:", self.protein_scheme)
        layout.addLayout(color_form)

        self.colors_btn = QPushButton("Setting", page)
        self.colors_btn.setSizePolicy(
            QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed
        )
        self.colors_btn.setToolTip(
            "Open the full colour editor for nucleotides, amino acids, "
            "feature types, and similarity databases."
        )
        self.colors_btn.clicked.connect(self.colorsClicked.emit)
        layout.addWidget(self.colors_btn, 0, Qt.AlignmentFlag.AlignLeft)

        self.status_label = QLabel("", page)
        self.status_label.setWordWrap(True)
        self.status_label.setStyleSheet("color: palette(mid);")
        layout.addWidget(self.status_label)
        layout.addStretch(1)
        self.sync_color_scheme_combos()
        return self._wrap_scroll(page)

    def _build_others_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)
        layout.addWidget(_section_label("Others", page))
        layout.addWidget(_section_label("Size", page))

        sizes = QFormLayout()
        sizes.setContentsMargins(0, 0, 0, 0)
        sizes.setHorizontalSpacing(8)
        sizes.setVerticalSpacing(6)
        sizes.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        sizes.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )

        spin_expand = QSizePolicy.Policy.Expanding
        self.bases_font_label = QLabel("DNA/protein labels:", page)
        self.bases_zoomed_in = NoWheelSpinBox(page)
        self.bases_zoomed_in.setRange(2, 36)
        self.bases_zoomed_in.setValue(6)  # single-seq default; MSA overview uses 4
        self.bases_zoomed_in.setSizePolicy(spin_expand, QSizePolicy.Policy.Fixed)
        self.bases_zoomed_in.setToolTip(
            "Point size of sequence letters when zoomed in enough to show "
            "individual bases (DNA) or residues (protein)."
        )
        self.bases_zoomed_in.valueChanged.connect(self.basesZoomedInChanged.emit)

        # Kept for settings / map wiring; not shown on Others.
        self.bases_zoomed_out = NoWheelSpinBox(page)
        self.bases_zoomed_out.setRange(1, 24)
        self.bases_zoomed_out.setValue(4)
        self.bases_zoomed_out.hide()
        self.bases_zoomed_out.valueChanged.connect(self.basesZoomedOutChanged.emit)

        self.annotation_label_size = NoWheelSpinBox(page)
        self.annotation_label_size.setRange(6, 24)
        self.annotation_label_size.setValue(7)
        self.annotation_label_size.setSizePolicy(spin_expand, QSizePolicy.Policy.Fixed)
        self.annotation_label_size.setToolTip(
            "Font size for labels on annotation bars (e.g. ORF / LTR / Pfam names)."
        )
        self.annotation_label_size.valueChanged.connect(
            self.annotationLabelSizeChanged.emit
        )

        self.annotation_thickness = NoWheelSpinBox(page)
        self.annotation_thickness.setRange(8, 40)
        self.annotation_thickness.setValue(18)
        self.annotation_thickness.setSizePolicy(spin_expand, QSizePolicy.Policy.Fixed)
        self.annotation_thickness.setToolTip(
            "Vertical height of annotation feature bars (pixels)."
        )
        self.annotation_thickness.valueChanged.connect(
            self.annotationThicknessChanged.emit
        )

        self.annotation_distance = NoWheelSpinBox(page)
        self.annotation_distance.setRange(0, 40)
        self.annotation_distance.setValue(0)
        self.annotation_distance.setSizePolicy(spin_expand, QSizePolicy.Policy.Fixed)
        self.annotation_distance.setToolTip(
            "Vertical gap between annotation rows (pixels)."
        )
        self.annotation_distance.valueChanged.connect(
            self.annotationDistanceChanged.emit
        )

        self.sequence_name_size = NoWheelSpinBox(page)
        self.sequence_name_size.setRange(8, 24)
        self.sequence_name_size.setValue(12)
        self.sequence_name_size.hide()
        self.sequence_name_size.valueChanged.connect(self.sequenceNameSizeChanged.emit)

        self.base_number_size = NoWheelSpinBox(page)
        self.base_number_size.setRange(6, 20)
        self.base_number_size.setValue(7)
        self.base_number_size.hide()
        self.base_number_size.valueChanged.connect(self.baseNumberSizeChanged.emit)

        sizes.addRow(self.bases_font_label, self.bases_zoomed_in)
        sizes.addRow("Annotation labels:", self.annotation_label_size)
        sizes.addRow("Annotation thickness:", self.annotation_thickness)
        self._annotation_distance_label = QLabel("Annotation distance:", page)
        sizes.addRow(self._annotation_distance_label, self.annotation_distance)
        layout.addLayout(sizes)
        layout.addStretch(1)
        return self._wrap_scroll(page)

    def _two_column_buttons(
        self, buttons: list[QPushButton], page: QWidget
    ) -> QWidget:
        """Lay out action buttons two per row, filling the rail width."""
        host = QWidget(page)
        host.setMinimumWidth(0)
        host.setSizePolicy(
            QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed
        )
        column = QVBoxLayout(host)
        column.setContentsMargins(0, 0, 0, 0)
        column.setSpacing(6)
        for index in range(0, len(buttons), 2):
            row = QHBoxLayout()
            row.setContentsMargins(0, 0, 0, 0)
            row.setSpacing(6)
            left = buttons[index]
            left.setMinimumWidth(0)
            left.setMaximumWidth(self._module_width)
            left.setSizePolicy(
                QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed
            )
            if index + 1 < len(buttons):
                right = buttons[index + 1]
                right.setMinimumWidth(0)
                right.setMaximumWidth(self._module_width)
                right.setSizePolicy(
                    QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed
                )
                row.addWidget(left, 1)
                row.addWidget(right, 1)
            else:
                # Odd last button spans the full row.
                row.addWidget(left, 1)
            column.addLayout(row)
        return host

    def _compact_run_more_row(
        self,
        page: QWidget,
        *,
        run_tooltip: str,
        more_tooltip: str,
        run_clicked,
        more_clicked,
    ) -> tuple[QWidget, QPushButton, QPushButton]:
        row = QWidget(page)
        row.setMinimumWidth(0)
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        run_btn = QPushButton("Run", row)
        run_btn.setToolTip(run_tooltip)
        run_btn.setFixedWidth(52)
        run_btn.clicked.connect(run_clicked)
        more_btn = QPushButton("Setting", row)
        more_btn.setToolTip(more_tooltip)
        more_btn.setFixedWidth(64)
        more_btn.clicked.connect(more_clicked)
        layout.addWidget(run_btn, 0)
        layout.addWidget(more_btn, 0)
        layout.addStretch(1)
        return row, run_btn, more_btn

    def _build_annotations_page(self) -> QWidget:
        page = QWidget()
        page.setMinimumWidth(0)
        layout = QVBoxLayout(page)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)
        layout.addWidget(_section_label("Annotations", page))
        layout.addWidget(_h_divider(page))

        # Hidden compatibility alias — Add moved to "Manually add feature" below.
        self.add_feature_btn = QPushButton("Add", page)
        self.add_feature_btn.hide()

        self.delete_btn = self._full_width_button("Delete", page)
        self.delete_btn.clicked.connect(self.deleteClicked.emit)
        self.restore_btn = self._full_width_button("Restore", page)
        self.restore_btn.clicked.connect(self.restoreClicked.emit)
        self.export_btn = self._full_width_button("Export", page)
        self.export_btn.clicked.connect(self.exportClicked.emit)
        self.save_btn = self._full_width_button("Save", page)
        self.save_btn.clicked.connect(self.saveClicked.emit)
        layout.addWidget(
            self._two_column_buttons(
                [
                    self.delete_btn,
                    self.restore_btn,
                    self.export_btn,
                    self.save_btn,
                ],
                page,
            )
        )

        self.annotation_job_status = QLabel("", page)
        self.annotation_job_status.setWordWrap(True)
        self.annotation_job_status.setStyleSheet("color: palette(mid);")
        self.annotation_job_status.hide()
        layout.addWidget(self.annotation_job_status)

        layout.addWidget(_h_divider(page))
        manual_section = _CollapsibleSection(
            "Manually add feature", page, expanded=False
        )
        self.manual_feature_host = QWidget(manual_section)
        self.manual_feature_host.setMinimumWidth(0)
        self.manual_feature_flow = _FlowLayout(self.manual_feature_host, spacing=6)
        self.manual_feature_host.setLayout(self.manual_feature_flow)
        self.manual_feature_checks: dict[str, QCheckBox] = {}
        self._manual_feature_exclusive = False
        manual_section.add_widget(self.manual_feature_host)
        self.manual_add_btn = QPushButton("Add", manual_section)
        self.manual_add_btn.setFixedWidth(52)
        self.manual_add_btn.setToolTip(
            "Add the checked feature type on the current sequence selection."
        )
        self.manual_add_btn.clicked.connect(self.addFeatureClicked.emit)
        self.manual_feature_settings_btn = QPushButton("Setting", manual_section)
        self.manual_feature_settings_btn.setFixedWidth(64)
        self.manual_feature_settings_btn.setToolTip(
            "Manage which feature types appear in the list."
        )
        self.manual_feature_settings_btn.clicked.connect(
            self.addFeatureSettingsClicked.emit
        )
        manual_btn_row = QHBoxLayout()
        manual_btn_row.setContentsMargins(0, 0, 0, 0)
        manual_btn_row.setSpacing(6)
        manual_btn_row.addWidget(self.manual_add_btn, 0)
        manual_btn_row.addWidget(self.manual_feature_settings_btn, 0)
        manual_btn_row.addStretch(1)
        manual_section.add_layout(manual_btn_row)
        layout.addWidget(manual_section)
        # Hidden combo kept for tests / legacy callers.
        self.manual_feature_combo = QComboBox(page)
        self.manual_feature_combo.hide()
        self.populate_manual_feature_types()

        layout.addWidget(_h_divider(page))
        feature_section = _CollapsibleSection(
            "Feature detection", page, expanded=False
        )
        self.auto_ann_all_check = QCheckBox("All", feature_section)
        feature_section.add_widget(self.auto_ann_all_check)
        self.auto_ann_checks: dict[object, QCheckBox] = {}
        auto_host = QWidget(feature_section)
        auto_host.setMinimumWidth(0)
        auto_flow = _FlowLayout(auto_host, spacing=6)
        auto_host.setLayout(auto_flow)
        for label in ("LTR", "TIR", "Tandem", "ORF"):
            check = QCheckBox(label, auto_host)
            check.setChecked(True)
            check.setMinimumWidth(0)
            auto_flow.addWidget(check)
            # Keys filled by SequenceMapPanel with AutoAnnotationMode enums.
            self.auto_ann_checks[label] = check
        feature_section.add_widget(auto_host)
        auto_btn_row, self.auto_ann_run_btn, self.auto_ann_settings_btn = (
            self._compact_run_more_row(
                feature_section,
                run_tooltip="Run feature detection",
                more_tooltip="Feature detection settings",
                run_clicked=self.featureDetectClicked.emit,
                more_clicked=self.featureDetectSettingsClicked.emit,
            )
        )
        feature_section.add_widget(auto_btn_row)
        layout.addWidget(feature_section)

        layout.addWidget(_h_divider(page))
        similarity_section = _CollapsibleSection(
            "Similarity annotation", page, expanded=False
        )
        self.similarity_ann_all_check = QCheckBox(
            "All databases", similarity_section
        )
        similarity_section.add_widget(self.similarity_ann_all_check)
        self.similarity_db_host = QWidget(similarity_section)
        self.similarity_db_host.setMinimumWidth(0)
        self.similarity_db_layout = QVBoxLayout(self.similarity_db_host)
        self.similarity_db_layout.setContentsMargins(0, 0, 0, 0)
        self.similarity_db_layout.setSpacing(2)
        similarity_section.add_widget(self.similarity_db_host)
        sim_btn_row, self.similarity_ann_run_btn, self.similarity_ann_settings_btn = (
            self._compact_run_more_row(
                similarity_section,
                run_tooltip="Run similarity annotation",
                more_tooltip=(
                    "Similarity annotation settings "
                    "(also used by Sequence Annotation)"
                ),
                run_clicked=self.similarityAnnotateClicked.emit,
                more_clicked=self.similaritySettingsClicked.emit,
            )
        )
        similarity_section.add_widget(sim_btn_row)
        layout.addWidget(similarity_section)

        layout.addWidget(_h_divider(page))
        transfer_section = _CollapsibleSection(
            "Transfer annotation", page, expanded=False
        )

        self.transfer_source_row = PathBrowseRow(
            transfer_section,
            placeholder="FASTA file or folder with annotations…",
            browse_title="Choose transfer source FASTA",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.ffn *.faa);;All files (*)"
            ),
            # Browse opens a menu: FASTA file... or Folder...
            mode="file_or_directory",
            minimum_width=120,
            browse_below=True,
            elide_middle=True,
        )
        self.transfer_source_row.set_default_tooltip(
            "Click a FASTA in Files or a folder in Folders to fill this path "
            "(click the field first to replace an existing path). Or use Browse. "
            "Only annotated sequences (*.seqmap_RV.json) are used."
        )
        self.transfer_source_row.path_changed.connect(self.transferSourceChanged.emit)
        self.transfer_source_row.text_changed.connect(
            lambda _text: self.transferSourceChanged.emit()
        )
        # Compatibility aliases (older code expected a combo + manage button).
        self.transfer_source = QComboBox(transfer_section)
        self.transfer_source.hide()
        self.transfer_manage_db_btn = QPushButton("Browse", transfer_section)
        self.transfer_manage_db_btn.hide()
        self.transfer_manage_db_btn.clicked.connect(
            self.transfer_source_row._on_browse_clicked
        )
        # Hidden target widgets: transfer always uses the currently open sequence.
        self.transfer_target = QComboBox(transfer_section)
        self.transfer_target.hide()
        self.transfer_target_browse = QPushButton("Browse…", transfer_section)
        self.transfer_target_browse.hide()

        self.transfer_selected_only = QCheckBox(
            "Selected features only", transfer_section
        )
        self.transfer_selected_only.hide()
        self.transfer_selected_only.setToolTip(
            "When checked, only features selected on the source track are "
            "transferred. Otherwise all features on the source are used."
        )
        self.transfer_allow_reverse = QCheckBox(
            "Allow reverse strand", transfer_section
        )
        self.transfer_allow_reverse.setChecked(True)

        self.transfer_min_identity_slider, self.transfer_min_identity_value_label = (
            self._transfer_threshold_slider_row(transfer_section)
        )
        self.transfer_min_coverage_slider, self.transfer_min_coverage_value_label = (
            self._transfer_threshold_slider_row(transfer_section)
        )
        # Legacy alias for code that read spin-box values.
        self.transfer_min_identity = NoWheelDoubleSpinBox(transfer_section)
        self.transfer_min_identity.hide()
        self.transfer_min_coverage = NoWheelDoubleSpinBox(transfer_section)
        self.transfer_min_coverage.hide()

        source_label = QLabel("Source", transfer_section)
        transfer_section.add_widget(source_label)
        transfer_section.add_widget(self.transfer_source_row)

        xfer_sliders = QVBoxLayout()
        xfer_sliders.setContentsMargins(0, 0, 0, 0)
        xfer_sliders.setSpacing(6)
        xfer_sliders.addWidget(
            self._threshold_slider_field(
                self.transfer_min_identity_slider,
                self.transfer_min_identity_value_label,
                label="Identity",
            )
        )
        xfer_sliders.addWidget(
            self._threshold_slider_field(
                self.transfer_min_coverage_slider,
                self.transfer_min_coverage_value_label,
                label="Coverage",
            )
        )
        transfer_section.add_layout(xfer_sliders)
        transfer_section.add_widget(self.transfer_selected_only)
        transfer_section.add_widget(self.transfer_allow_reverse)

        self.transfer_run_btn = QPushButton("Run", transfer_section)
        self.transfer_run_btn.setFixedWidth(52)
        self.transfer_run_btn.setToolTip(
            "Hybrid: whole-seq alignment + feature lift (MAFFT fast if <10 kb, "
            "else minimap2), then blastn for leftovers. Matches become source "
            "“Transfer (…)”."
        )
        self.transfer_run_btn.clicked.connect(
            lambda: self.transferRequested.emit()
        )
        transfer_run_row = QHBoxLayout()
        transfer_run_row.setContentsMargins(0, 0, 0, 0)
        transfer_run_row.addWidget(self.transfer_run_btn, 0)
        transfer_run_row.addStretch(1)
        transfer_section.add_layout(transfer_run_row)
        layout.addWidget(transfer_section)

        self._annotation_sections: dict[str, _CollapsibleSection] = {
            "manual": manual_section,
            "feature_detection": feature_section,
            "similarity": similarity_section,
            "transfer": transfer_section,
        }
        for section in self._annotation_sections.values():
            section.expandedChanged.connect(self._save_rail_ui_state)

        layout.addStretch(1)
        return self._wrap_scroll(page)

    def _build_find_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)
        layout.addWidget(_section_label("Find", page))
        layout.addWidget(_h_divider(page))

        self.find_query = QLineEdit(page)
        self.find_query.setPlaceholderText("Sequence or name…")
        self.find_query.setClearButtonEnabled(True)
        self.find_query.returnPressed.connect(
            lambda: self.findRequested.emit(1)
        )
        find_form = QFormLayout()
        find_form.setContentsMargins(0, 0, 0, 0)
        find_form.setHorizontalSpacing(8)
        find_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        find_form.addRow("Find:", self.find_query)
        layout.addLayout(find_form)

        layout.addWidget(_section_label("Search in", page))
        self.find_in_sequence = QCheckBox("Sequence", page)
        self.find_in_sequence.setChecked(True)
        self.find_in_sequence.toggled.connect(
            lambda *_a: self.sync_find_option_enabled()
        )
        layout.addWidget(self.find_in_sequence)

        self.find_nucleotides = QRadioButton("Nucleotides", page)
        self.find_nucleotides.setChecked(True)
        self.find_protein = QRadioButton("Protein (6 frames)", page)
        mode_group = QButtonGroup(page)
        mode_group.addButton(self.find_nucleotides)
        mode_group.addButton(self.find_protein)
        mode_host = QWidget(page)
        mode_layout = QVBoxLayout(mode_host)
        mode_layout.setContentsMargins(16, 0, 0, 0)
        mode_layout.setSpacing(2)
        mode_layout.addWidget(self.find_nucleotides)
        mode_layout.addWidget(self.find_protein)
        layout.addWidget(mode_host)

        opts = QWidget(page)
        opts_layout = QVBoxLayout(opts)
        opts_layout.setContentsMargins(16, 0, 0, 0)
        opts_layout.setSpacing(2)
        self.find_reverse_complement = QCheckBox("Include reverse complement", opts)
        self.find_reverse_complement.setChecked(True)
        self.find_iupac = QCheckBox("Interpret IUPAC codes", opts)
        self.find_iupac.setChecked(True)
        self.find_allow_gaps = QCheckBox("Allow gaps in sequence", opts)
        opts_layout.addWidget(self.find_reverse_complement)
        opts_layout.addWidget(
            field_with_help(
                self.find_iupac,
                "Interpret IUPAC codes",
                _FIND_IUPAC_HELP,
                parent=opts,
            )
        )
        opts_layout.addWidget(self.find_allow_gaps)
        self.find_max_mismatches = NoWheelSpinBox(opts)
        self.find_max_mismatches.setRange(0, 20)
        self.find_max_mismatches.setValue(0)
        mm_row = QHBoxLayout()
        mm_row.setContentsMargins(0, 0, 0, 0)
        mm_row.addWidget(QLabel("Max mismatches:", opts))
        mm_row.addWidget(self.find_max_mismatches, 1)
        opts_layout.addLayout(mm_row)
        layout.addWidget(opts)

        self.find_in_names = QCheckBox("Sequence names", page)
        self.find_in_names.setChecked(True)
        self.find_in_annotations = QCheckBox("Annotation names", page)
        self.find_in_annotations.setChecked(True)
        self.find_case_sensitive = QCheckBox("Case sensitive", page)
        layout.addWidget(self.find_in_names)
        layout.addWidget(self.find_in_annotations)
        layout.addWidget(self.find_case_sensitive)

        self.find_next_btn = self._full_width_button("Find next", page)
        self.find_prev_btn = self._full_width_button("Find previous", page)
        self.find_all_btn = self._full_width_button("Find all", page)
        self.find_next_btn.clicked.connect(lambda: self.findRequested.emit(1))
        self.find_prev_btn.clicked.connect(lambda: self.findRequested.emit(-1))
        self.find_all_btn.clicked.connect(lambda: self.findRequested.emit(0))
        find_nav = QHBoxLayout()
        find_nav.setContentsMargins(0, 0, 0, 0)
        find_nav.setSpacing(6)
        find_nav.addWidget(self.find_prev_btn, 1)
        find_nav.addWidget(self.find_next_btn, 1)
        layout.addLayout(find_nav)
        layout.addWidget(self.find_all_btn)

        # Status used to show grey match text under the buttons; keep a hidden
        # label so set_find_status callers stay valid without showing it.
        self.find_status = QLabel("", page)
        self.find_status.hide()
        layout.addStretch(1)
        self.sync_find_option_enabled()
        return self._wrap_scroll(page)

    def find_request(self):
        """Build a FindRequest from the Find tab controls."""
        from repeatviewer.backends.find_in_file import FindRequest
        from repeatviewer.backends.motif_search import FindOptions

        iupac = self.find_iupac.isChecked()
        options = FindOptions(
            protein=self.find_protein.isChecked(),
            both_strands=self.find_reverse_complement.isChecked(),
            iupac_in_query=iupac,
            iupac_in_sequence=iupac,
            allow_gaps=self.find_allow_gaps.isChecked(),
            max_mismatches=int(self.find_max_mismatches.value()),
        )
        return FindRequest(
            query=self.find_query.text().strip(),
            options=options,
            in_sequence=self.find_in_sequence.isChecked(),
            in_names=self.find_in_names.isChecked(),
            in_annotations=self.find_in_annotations.isChecked(),
            case_sensitive=self.find_case_sensitive.isChecked(),
        )

    def set_find_status(self, text: str) -> None:
        # Find panel no longer shows status text; keep API for map-viewer calls.
        self.find_status.setText(text or "")

    def sync_find_option_enabled(self) -> None:
        on = self.find_in_sequence.isChecked()
        for widget in (
            self.find_nucleotides,
            self.find_protein,
            self.find_reverse_complement,
            self.find_iupac,
            self.find_allow_gaps,
            self.find_max_mismatches,
        ):
            widget.setEnabled(on)

    def _on_show_translation_toggled(self, checked: bool) -> None:
        self.translation_more_btn.setEnabled(checked)
        self.showTranslationChanged.emit(checked)
        self._emit_translation_frames()

    def _emit_translation_frames(self, *_args) -> None:
        frames = frozenset(
            frame for frame, check in self.frame_checks.items() if check.isChecked()
        )
        self.translationFramesChanged.emit(frames)

    def selected_translation_frames(self) -> frozenset[int]:
        if not self.show_translation.isChecked():
            return frozenset()
        frames = frozenset(
            frame for frame, check in self.frame_checks.items() if check.isChecked()
        )
        return frames or frozenset({1})

    def translation_three_letter(self) -> bool:
        return bool(self._translation_three_letter)

    def translation_genetic_code(self) -> int:
        return normalize_genetic_code_id(self._translation_genetic_code)

    def set_translation_options(
        self,
        *,
        frames: frozenset[int] | set[int] | None = None,
        three_letter: bool | None = None,
        genetic_code: int | None = None,
    ) -> None:
        """Update stored translation options without opening the dialog."""
        if frames is not None:
            selected = set(frames) or {1}
            for frame, check in self.frame_checks.items():
                check.blockSignals(True)
                check.setChecked(frame in selected)
                check.blockSignals(False)
        if three_letter is not None:
            self._translation_three_letter = bool(three_letter)
        if genetic_code is not None:
            self._translation_genetic_code = normalize_genetic_code_id(genetic_code)

    def _open_translation_settings(self) -> None:
        current = frozenset(
            frame for frame, check in self.frame_checks.items() if check.isChecked()
        ) or frozenset({1})
        dialog = TranslationSettingsDialog(
            self,
            frames=current,
            three_letter=self._translation_three_letter,
            genetic_code=self._translation_genetic_code,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self.set_translation_options(
            frames=dialog.selected_frames(),
            three_letter=dialog.show_three_letter(),
            genetic_code=dialog.genetic_code_id(),
        )
        self._emit_translation_frames()
        self.translationOptionsChanged.emit()

    def sync_color_scheme_combos(self) -> None:
        """Match scheme dropdowns to the currently saved residue palettes."""
        nuc_name = (
            "Color-blind"
            if load_color_blind_nucleotide_scheme()
            else "Default"
        )
        idx = self.nucleotide_scheme.findText(nuc_name)
        self.nucleotide_scheme.blockSignals(True)
        self.nucleotide_scheme.setCurrentIndex(max(0, idx))
        self.nucleotide_scheme.blockSignals(False)

    def _on_nucleotide_scheme_changed(self, index: int) -> None:
        if index < 0 or index >= len(NUCLEOTIDE_COLOR_PRESETS):
            return
        name, palette = NUCLEOTIDE_COLOR_PRESETS[index]
        save_custom_nucleotide_backgrounds(dict(palette))
        save_color_blind_nucleotide_scheme(name == "Color-blind")
        self.residueColorsChanged.emit()

    def _on_protein_scheme_changed(self, index: int) -> None:
        if index < 0 or index >= len(AMINO_ACID_COLOR_PRESETS):
            return
        _name, palette = AMINO_ACID_COLOR_PRESETS[index]
        save_custom_amino_acid_backgrounds(dict(palette))
        self.residueColorsChanged.emit()

    def set_sequence_kind(self, is_protein: bool) -> None:
        """Compatibility hook; size label is fixed as DNA/protein labels."""
        del is_protein

    def set_overview_display_mode(self, enabled: bool = True) -> None:
        """Adapt Display/Annotation for multi-sequence overview.

        Keeps Annotation (manual add, transfer, feature detection / similarity).
        Hides single-map-only complement/translation toggles.
        Module width matches the single-sequence Annotation rail.
        """
        overview = bool(enabled)
        self._overview_mode = overview
        for widget in (
            self.show_complement,
            self.show_translation,
            self.translation_more_btn,
            self.split_tracks_by_source,
        ):
            widget.setVisible(not overview)
        # Keep Tracks header + annotation-track toggle in overview.
        self._tracks_divider.setVisible(True)
        self._tracks_label.setVisible(True)
        self.show_annotation_track.setVisible(True)
        self._tab_buttons[self.TAB_ANNOTATIONS].setVisible(True)
        sections = getattr(self, "_annotation_sections", None) or {}
        for section_id in ("manual", "transfer"):
            section = sections.get(section_id)
            if section is not None:
                section.setVisible(True)
        if overview:
            self.auto_ann_run_btn.setToolTip(
                "Run feature detection on selected sequences (Ctrl+A for all)"
            )
            self.similarity_ann_run_btn.setToolTip(
                "Run similarity annotation on selected sequences (Ctrl+A for all)"
            )
            # Same spin as single-seq annotation-row gap; in overview it spaces sequences.
            self._annotation_distance_label.setText("Sequence distance:")
            self.annotation_distance.setToolTip(
                "Vertical gap between sequences (pixels)."
            )
        else:
            self.auto_ann_run_btn.setToolTip("Run feature detection")
            self.similarity_ann_run_btn.setToolTip("Run similarity annotation")
            self._annotation_distance_label.setText("Annotation distance:")
            self.annotation_distance.setToolTip(
                "Vertical gap between annotation rows (pixels)."
            )
        for btn in (self.delete_btn, self.restore_btn, self.export_btn):
            btn.setVisible(True)
        self._module_width = _MODULE_WIDTH
        self._apply_module_geometry()

    def is_panel_collapsed(self) -> bool:
        return bool(self._panel_collapsed)

    def set_panel_collapsed(self, collapsed: bool, *, persist: bool = True) -> None:
        """Hide or show the control pages; tab strip always stays visible."""
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
        # Re-click active tab → collapse; click while collapsed or another tab → show.
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
        self.tabChanged.emit(index)
        self._save_rail_ui_state()

    def set_current_tab(self, index: int) -> None:
        index = max(0, min(index, self._stack.count() - 1))
        self._active_tab = index
        if self._panel_collapsed:
            self.set_panel_collapsed(False)
        else:
            self._stack.setCurrentIndex(index)
            self._sync_tab_checked_state(index)
            self._save_rail_ui_state()

    def _settings(self) -> QSettings:
        settings = QSettings()
        settings.beginGroup(_SETTINGS_GROUP)
        return settings

    def _restore_rail_ui_state(self) -> None:
        """Default: Annotation tab, all sections collapsed; restore last state."""
        settings = self._settings()
        try:
            raw_tab = settings.value(_SETTINGS_TAB_KEY, _DEFAULT_TAB_INDEX)
            try:
                tab = int(raw_tab)  # type: ignore[arg-type]
            except (TypeError, ValueError):
                tab = _DEFAULT_TAB_INDEX
            tab = max(0, min(tab, self._stack.count() - 1))
            self._active_tab = tab
            self._stack.setCurrentIndex(tab)

            collapsed_raw = settings.value(_SETTINGS_COLLAPSED_KEY, False)
            if isinstance(collapsed_raw, str):
                collapsed = collapsed_raw.strip().lower() in {"1", "true", "yes"}
            else:
                collapsed = bool(collapsed_raw)
            self._panel_collapsed = False  # force apply below
            self.set_panel_collapsed(collapsed, persist=False)
            if not collapsed:
                self._sync_tab_checked_state(tab)

            sections = getattr(self, "_annotation_sections", None) or {}
            saved = settings.value(_SETTINGS_SECTIONS_KEY, [])
            if isinstance(saved, str):
                expanded_ids = {saved} if saved else set()
            elif isinstance(saved, (list, tuple)):
                expanded_ids = {str(item) for item in saved}
            else:
                expanded_ids = set()
            for section_id, section in sections.items():
                section.set_expanded(
                    section_id in expanded_ids, notify=False
                )
        finally:
            settings.endGroup()

    def _save_rail_ui_state(self, *_args) -> None:
        settings = self._settings()
        try:
            settings.setValue(_SETTINGS_TAB_KEY, int(self._active_tab))
            settings.setValue(_SETTINGS_COLLAPSED_KEY, bool(self._panel_collapsed))
            sections = getattr(self, "_annotation_sections", None) or {}
            expanded = [
                section_id
                for section_id, section in sections.items()
                if section.is_expanded()
            ]
            settings.setValue(_SETTINGS_SECTIONS_KEY, expanded)
        finally:
            settings.endGroup()

    def _transfer_threshold_slider_row(
        self, parent: QWidget
    ) -> tuple[NoWheelSlider, QLabel]:
        slider = NoWheelSlider(Qt.Orientation.Horizontal, parent)
        slider.setObjectName("annotationZoomSlider")
        slider.setRange(0, 100)
        slider.setValue(50)
        slider.setMinimumWidth(80)
        slider.setFixedHeight(18)
        slider.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        apply_annotation_zoom_slider_style(slider)
        value_label = QLabel("0.50", parent)
        value_label.setMinimumWidth(36)
        value_label.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )

        def _sync(value: int) -> None:
            frac = max(0, min(100, int(value))) / 100.0
            value_label.setText(f"{frac:.2f}")
            self.transferThresholdChanged.emit()

        slider.valueChanged.connect(_sync)
        return slider, value_label

    @staticmethod
    def _threshold_slider_field(
        slider: NoWheelSlider,
        value_label: QLabel,
        *,
        label: str = "",
    ) -> QWidget:
        row = QWidget(slider.parentWidget())
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        if label:
            name = QLabel(label, row)
            # Fixed column + right align so "Identity"/"Coverage" share an end
            # edge and both slider tracks stay the same length.
            fm = QFontMetrics(name.font())
            name.setFixedWidth(fm.horizontalAdvance("Coverage") + 2)
            name.setAlignment(
                Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
            )
            name.setToolTip(f"Minimum {label.lower()} threshold")
            layout.addWidget(name, 0)
        layout.addWidget(slider, 1)
        layout.addWidget(value_label, 0)
        return row

    def transfer_min_identity_value(self) -> float:
        return float(self.transfer_min_identity_slider.value()) / 100.0

    def transfer_min_coverage_value(self) -> float:
        return float(self.transfer_min_coverage_slider.value()) / 100.0

    def populate_manual_feature_types(
        self, buttons: list[tuple[str, str]] | None = None
    ) -> None:
        """Fill the Manually add feature type list (label, feature_type)."""
        entries = list(buttons) if buttons is not None else load_add_feature_buttons()
        previous = str(self.manual_feature_combo.currentData() or "").strip()
        if not previous:
            _label, previous = self.selected_manual_feature()
            if previous == "others" and not self.manual_feature_checks:
                previous = ""

        self.manual_feature_combo.blockSignals(True)
        try:
            self.manual_feature_combo.clear()
            for label, feature_type in entries:
                self.manual_feature_combo.addItem(str(label), str(feature_type))
            if previous:
                idx = self.manual_feature_combo.findData(previous)
                if idx >= 0:
                    self.manual_feature_combo.setCurrentIndex(idx)
        finally:
            self.manual_feature_combo.blockSignals(False)

        while self.manual_feature_flow.count():
            item = self.manual_feature_flow.takeAt(0)
            if item.widget():
                item.widget().deleteLater()
        self.manual_feature_checks.clear()
        self._manual_feature_exclusive = True
        restore_type = previous or (
            str(entries[0][1]) if entries else "others"
        )
        for label, feature_type in entries:
            check = QCheckBox(str(label), self.manual_feature_host)
            check.setProperty("feature_type", str(feature_type))
            check.setChecked(str(feature_type) == restore_type)
            check.setMinimumWidth(0)
            check.toggled.connect(self._on_manual_feature_check_toggled)
            self.manual_feature_flow.addWidget(check)
            self.manual_feature_checks[str(feature_type)] = check
        self.manual_feature_host.updateGeometry()

    def _on_manual_feature_check_toggled(self, checked: bool) -> None:
        if not self._manual_feature_exclusive:
            return
        sender = self.sender()
        if not isinstance(sender, QCheckBox):
            return
        if checked:
            for check in self.manual_feature_checks.values():
                if check is not sender and check.isChecked():
                    check.blockSignals(True)
                    check.setChecked(False)
                    check.blockSignals(False)
            return
        if any(c.isChecked() for c in self.manual_feature_checks.values()):
            return
        sender.blockSignals(True)
        sender.setChecked(True)
        sender.blockSignals(False)

    def selected_manual_feature(self) -> tuple[str, str]:
        """Return ``(label, feature_type)`` for the checked manual feature type."""
        for check in self.manual_feature_checks.values():
            if check.isChecked():
                feature_type = str(check.property("feature_type") or "").strip()
                label = str(check.text() or "").strip() or feature_type
                if feature_type:
                    return label, feature_type
        combo = self.manual_feature_combo
        if combo.count() > 0:
            label = str(combo.currentText() or "").strip() or "Others"
            feature_type = str(combo.currentData() or "").strip() or "others"
            return label, feature_type
        return ("Others", "others")

    def set_status(self, text: str) -> None:
        self.status_label.setText(text or "")

    def set_annotation_job_status(self, text: str) -> None:
        """Live progress for Feature detection / Similarity / Transfer on Annotation tab."""
        msg = str(text or "").strip()
        if self.annotation_job_status.text() == msg and (
            self.annotation_job_status.isVisible() == bool(msg)
        ):
            return
        self.annotation_job_status.setText(msg)
        self.annotation_job_status.setVisible(bool(msg))
        # Keep Display-tab status in sync for consistency.
        self.status_label.setText(msg)

    def set_transfer_summary(self, text: str) -> None:
        """Transfer progress is shown on the map preview, not in the rail."""
        del text

    def set_sequence_name_font_point_size(self, point_size: int) -> None:
        font = self.sequence_combo.font()
        font.setPointSize(max(8, min(24, int(point_size))))
        self.sequence_combo.setFont(font)
        self.sequence_search.setFont(font)

    def browse_target_fasta(self, parent: QWidget | None = None) -> Path | None:
        path, _ = QFileDialog.getOpenFileName(
            parent or self,
            "Choose target FASTA",
            "",
            "FASTA (*.fa *.fasta *.fna *.fas);;All files (*)",
        )
        if not path:
            return None
        return Path(path)

    def clear_similarity_db_checks(self) -> None:
        while self.similarity_db_layout.count():
            item = self.similarity_db_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()

    def add_similarity_db_widget(self, widget: QWidget) -> None:
        widget.setParent(self.similarity_db_host)
        self.similarity_db_layout.addWidget(widget)

    def add_similarity_db_check(self, check: QCheckBox) -> None:
        self.add_similarity_db_widget(check)
