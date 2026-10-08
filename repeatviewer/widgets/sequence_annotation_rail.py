"""Geneious-style right rail for the Annotation preview (Display / Annotation / Find / Others)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QPoint, QRectF, QSize, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPen
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
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

from repeatviewer.backends.genetic_codes import (
    genetic_code_choices,
    normalize_genetic_code_id,
)
from repeatviewer.backends.msa_plot import (
    AMINO_ACID_COLOR_PRESETS,
    NUCLEOTIDE_COLOR_PRESETS,
)
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.preview_zoom import (
    load_color_blind_nucleotide_scheme,
    save_color_blind_nucleotide_scheme,
    save_custom_amino_acid_backgrounds,
    save_custom_nucleotide_backgrounds,
)
from repeatviewer.widgets.spin_boxes import NoWheelDoubleSpinBox, NoWheelSpinBox

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
_GENETIC_CODE_HELP = (
    "The genetic code table controls how DNA codons are translated into amino acids "
    "(same IDs as GenBank <code>/transl_table</code>). Use <b>1. Standard</b> for most "
    "nuclear genomes; choose a mitochondrial, plastid, or other table when your sequence "
    "comes from those systems.<br><br>"
    "For more detailed explanation please visit:<br>"
    f'<a href="{_NCBI_GENETIC_CODES_URL}">{_NCBI_GENETIC_CODES_URL}</a>'
)

_RAIL_TAB_WIDTH = 24
_MODULE_WIDTH = 300
_TAB_BUTTON_HEIGHT_MIN = 36
_TAB_PILL_RADIUS = 7
_TAB_TEXT_PAD = 4  # vertical padding around rotated label

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

    addFeatureClicked = Signal()
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
    tabChanged = Signal(int)
    findRequested = Signal(int)  # +1 next, -1 previous, 0 find all

    # Display / size signals (live apply).
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
        self.setFixedWidth(_MODULE_WIDTH + _RAIL_TAB_WIDTH + 1)

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        self._stack = QStackedWidget(self)
        self._stack.setFixedWidth(_MODULE_WIDTH)
        self._stack.setSizePolicy(
            QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Expanding
        )

        self._show_page = self._build_show_page()
        self._annotations_page = self._build_annotations_page()
        self._find_page = self._build_find_page()
        self._others_page = self._build_others_page()
        self._stack.addWidget(self._show_page)
        self._stack.addWidget(self._annotations_page)
        self._stack.addWidget(self._find_page)
        self._stack.addWidget(self._others_page)

        divider = QFrame(self)
        divider.setFrameShape(QFrame.Shape.VLine)
        divider.setFrameShadow(QFrame.Shadow.Plain)
        divider.setFixedWidth(1)
        divider.setStyleSheet("color: palette(mid);")

        self._tab_strip = QWidget(self)
        self._tab_strip.setFixedWidth(_RAIL_TAB_WIDTH)
        self._tab_strip.setObjectName("annotationRailTabs")
        tab_layout = QVBoxLayout(self._tab_strip)
        tab_layout.setContentsMargins(1, 2, 1, 2)
        tab_layout.setSpacing(2)

        self._tab_group = QButtonGroup(self)
        self._tab_group.setExclusive(True)
        self._tab_buttons: list[_RotatedTabButton] = []
        for index, (label, tip) in enumerate(
            (
                ("Display", "Display — sequence display and colors"),
                ("Annotation", "Annotations — detect, edit, and transfer features"),
                ("Find", "Find — motif search in the active sequence"),
                ("Others", "Others — sizes and extra display options"),
            )
        ):
            btn = _RotatedTabButton(label, self._tab_strip)
            btn.setToolTip(tip)
            self._tab_group.addButton(btn, index)
            tab_layout.addWidget(btn, 0, Qt.AlignmentFlag.AlignHCenter)
            self._tab_buttons.append(btn)
        tab_layout.addStretch(1)
        self._tab_group.idClicked.connect(self._on_tab_clicked)
        self._tab_buttons[0].setChecked(True)
        # Ensure the initially selected tab paints with the active highlight.
        self._tab_buttons[0].update()

        root.addWidget(self._stack, 0)
        root.addWidget(divider, 0)
        root.addWidget(self._tab_strip, 0)

    def _wrap_scroll(self, body: QWidget) -> QScrollArea:
        body.setMinimumWidth(0)
        body.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Minimum,
        )
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        scroll.setWidget(body)
        # Keep page content as wide as the viewport (minus scrollbar when shown).
        scroll.viewport().setAutoFillBackground(False)
        return scroll

    def _full_width_button(self, text: str, page: QWidget) -> QPushButton:
        btn = QPushButton(text, page)
        btn.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        btn.setMinimumWidth(0)
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

        layout.addWidget(_h_divider(page))
        layout.addWidget(_section_label("Tracks", page))

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
        self.translation_more_btn = QPushButton("More…", page)
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

        self.colors_btn = QPushButton("Advanced color setting…", page)
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
        self.bases_zoomed_in.setRange(6, 36)
        self.bases_zoomed_in.setValue(13)
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
        sizes.addRow("Annotation distance:", self.annotation_distance)
        layout.addLayout(sizes)
        layout.addStretch(1)
        return self._wrap_scroll(page)

    def _build_annotations_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)
        layout.addWidget(_section_label("Annotations", page))

        hint = QLabel(
            "Select a region on the sequence, then Add feature. "
            "Or run feature / similarity detection below.",
            page,
        )
        hint.setWordWrap(True)
        layout.addWidget(hint)

        self.add_feature_btn = self._full_width_button("Add feature", page)
        self.add_feature_btn.clicked.connect(self.addFeatureClicked.emit)
        layout.addWidget(self.add_feature_btn)

        self.delete_btn = self._full_width_button("Delete", page)
        self.delete_btn.clicked.connect(self.deleteClicked.emit)
        self.restore_btn = self._full_width_button("Restore", page)
        self.restore_btn.clicked.connect(self.restoreClicked.emit)
        self.export_btn = self._full_width_button("Export", page)
        self.export_btn.clicked.connect(self.exportClicked.emit)
        self.save_btn = self._full_width_button("Save", page)
        self.save_btn.clicked.connect(self.saveClicked.emit)
        for btn in (
            self.delete_btn,
            self.restore_btn,
            self.export_btn,
            self.save_btn,
        ):
            layout.addWidget(btn)

        layout.addWidget(_section_label("Feature detection", page))
        self.auto_ann_all_check = QCheckBox("All", page)
        layout.addWidget(self.auto_ann_all_check)
        self.auto_ann_checks: dict[object, QCheckBox] = {}
        auto_row = QWidget(page)
        auto_layout = QHBoxLayout(auto_row)
        auto_layout.setContentsMargins(0, 0, 0, 0)
        auto_layout.setSpacing(6)
        for label in ("LTR", "TIR", "Tandem", "ORF"):
            check = QCheckBox(label, auto_row)
            check.setChecked(True)
            auto_layout.addWidget(check)
            # Keys filled by SequenceMapPanel with AutoAnnotationMode enums.
            self.auto_ann_checks[label] = check
        auto_layout.addStretch(1)
        layout.addWidget(auto_row)
        auto_btn_row = QWidget(page)
        auto_btn_layout = QHBoxLayout(auto_btn_row)
        auto_btn_layout.setContentsMargins(0, 0, 0, 0)
        auto_btn_layout.setSpacing(6)
        self.auto_ann_run_btn = self._full_width_button(
            "Run feature detection", auto_btn_row
        )
        self.auto_ann_run_btn.clicked.connect(self.featureDetectClicked.emit)
        self.auto_ann_settings_btn = QToolButton(auto_btn_row)
        self.auto_ann_settings_btn.setText("⚙")
        self.auto_ann_settings_btn.setToolTip("Feature detection settings")
        self.auto_ann_settings_btn.clicked.connect(self.featureDetectSettingsClicked.emit)
        auto_btn_layout.addWidget(self.auto_ann_run_btn, 1)
        auto_btn_layout.addWidget(self.auto_ann_settings_btn, 0)
        layout.addWidget(auto_btn_row)

        layout.addWidget(_section_label("Similarity annotation", page))
        self.similarity_ann_all_check = QCheckBox("All databases", page)
        layout.addWidget(self.similarity_ann_all_check)
        self.similarity_db_host = QWidget(page)
        self.similarity_db_layout = QVBoxLayout(self.similarity_db_host)
        self.similarity_db_layout.setContentsMargins(0, 0, 0, 0)
        self.similarity_db_layout.setSpacing(2)
        layout.addWidget(self.similarity_db_host)
        sim_btn_row = QWidget(page)
        sim_btn_layout = QHBoxLayout(sim_btn_row)
        sim_btn_layout.setContentsMargins(0, 0, 0, 0)
        sim_btn_layout.setSpacing(6)
        self.similarity_ann_run_btn = self._full_width_button(
            "Run similarity", sim_btn_row
        )
        self.similarity_ann_run_btn.clicked.connect(self.similarityAnnotateClicked.emit)
        self.similarity_ann_settings_btn = QToolButton(sim_btn_row)
        self.similarity_ann_settings_btn.setText("⚙")
        self.similarity_ann_settings_btn.setToolTip("Similarity annotation settings")
        self.similarity_ann_settings_btn.clicked.connect(
            self.similaritySettingsClicked.emit
        )
        sim_btn_layout.addWidget(self.similarity_ann_run_btn, 1)
        sim_btn_layout.addWidget(self.similarity_ann_settings_btn, 0)
        layout.addWidget(sim_btn_row)

        layout.addWidget(_section_label("Transfer from another sequence", page))
        intro = QLabel(
            "Lift annotations from a source sequence onto a target via "
            "minimap2 alignment.",
            page,
        )
        intro.setWordWrap(True)
        layout.addWidget(intro)

        self.transfer_source = QComboBox(page)
        self.transfer_target = QComboBox(page)
        self.transfer_target_browse = QPushButton("Browse…", page)
        self.transfer_target_browse.setToolTip(
            "Choose another FASTA as the transfer target"
        )
        target_row = QWidget(page)
        target_layout = QHBoxLayout(target_row)
        target_layout.setContentsMargins(0, 0, 0, 0)
        target_layout.setSpacing(4)
        target_layout.addWidget(self.transfer_target, 1)
        target_layout.addWidget(self.transfer_target_browse, 0)

        self.transfer_selected_only = QCheckBox("Selected features only", page)
        self.transfer_selected_only.setToolTip(
            "When checked, only features selected on the track are transferred. "
            "Otherwise all features on the source are used."
        )
        self.transfer_allow_reverse = QCheckBox("Allow reverse strand", page)
        self.transfer_allow_reverse.setChecked(True)

        self.transfer_min_identity = NoWheelDoubleSpinBox(page)
        self.transfer_min_identity.setRange(0.0, 1.0)
        self.transfer_min_identity.setDecimals(2)
        self.transfer_min_identity.setSingleStep(0.05)
        self.transfer_min_identity.setValue(0.50)
        self.transfer_min_identity.setMinimumWidth(0)

        self.transfer_min_coverage = NoWheelDoubleSpinBox(page)
        self.transfer_min_coverage.setRange(0.0, 1.0)
        self.transfer_min_coverage.setDecimals(2)
        self.transfer_min_coverage.setSingleStep(0.05)
        self.transfer_min_coverage.setValue(0.50)
        self.transfer_min_coverage.setMinimumWidth(0)

        xfer_form = QFormLayout()
        xfer_form.setContentsMargins(0, 0, 0, 0)
        xfer_form.setHorizontalSpacing(8)
        xfer_form.setVerticalSpacing(6)
        xfer_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        xfer_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        xfer_form.addRow("Source:", self.transfer_source)
        xfer_form.addRow("Target:", target_row)
        xfer_form.addRow("Min identity:", self.transfer_min_identity)
        xfer_form.addRow("Min coverage:", self.transfer_min_coverage)
        layout.addLayout(xfer_form)
        layout.addWidget(self.transfer_selected_only)
        layout.addWidget(self.transfer_allow_reverse)

        self.transfer_run_btn = self._full_width_button("Transfer", page)
        self.transfer_run_btn.clicked.connect(self.transferRequested.emit)
        layout.addWidget(self.transfer_run_btn)

        self.transfer_summary = QLabel("", page)
        self.transfer_summary.setWordWrap(True)
        layout.addWidget(self.transfer_summary)
        layout.addStretch(1)
        return self._wrap_scroll(page)

    def _build_find_page(self) -> QWidget:
        page = QWidget()
        layout = QVBoxLayout(page)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)
        layout.addWidget(_section_label("Find", page))

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
        self.find_reverse_complement = QCheckBox("Also find reverse complement", opts)
        self.find_reverse_complement.setChecked(True)
        self.find_iupac_query = QCheckBox("Interpret IUPAC codes in query", opts)
        self.find_iupac_query.setChecked(True)
        self.find_iupac_sequence = QCheckBox("Interpret IUPAC codes in sequence", opts)
        self.find_allow_gaps = QCheckBox("Allow gaps in sequence", opts)
        for check in (
            self.find_reverse_complement,
            self.find_iupac_query,
            self.find_iupac_sequence,
            self.find_allow_gaps,
        ):
            opts_layout.addWidget(check)
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
        layout.addWidget(self.find_next_btn)
        layout.addWidget(self.find_prev_btn)
        layout.addWidget(self.find_all_btn)

        self.find_status = QLabel("", page)
        self.find_status.setWordWrap(True)
        self.find_status.setStyleSheet("color: palette(mid);")
        layout.addWidget(self.find_status)
        layout.addStretch(1)
        self.sync_find_option_enabled()
        return self._wrap_scroll(page)

    def find_request(self):
        """Build a FindRequest from the Find tab controls."""
        from repeatviewer.backends.find_in_file import FindRequest
        from repeatviewer.backends.motif_search import FindOptions

        options = FindOptions(
            protein=self.find_protein.isChecked(),
            both_strands=self.find_reverse_complement.isChecked(),
            iupac_in_query=self.find_iupac_query.isChecked(),
            iupac_in_sequence=self.find_iupac_sequence.isChecked(),
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
        self.find_status.setText(text or "")

    def sync_find_option_enabled(self) -> None:
        on = self.find_in_sequence.isChecked()
        for widget in (
            self.find_nucleotides,
            self.find_protein,
            self.find_reverse_complement,
            self.find_iupac_query,
            self.find_iupac_sequence,
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

    def _on_tab_clicked(self, index: int) -> None:
        self._stack.setCurrentIndex(index)
        self.tabChanged.emit(index)

    def set_current_tab(self, index: int) -> None:
        index = max(0, min(index, self._stack.count() - 1))
        self._tab_buttons[index].setChecked(True)
        self._stack.setCurrentIndex(index)
        for btn in self._tab_buttons:
            btn.update()

    def set_status(self, text: str) -> None:
        self.status_label.setText(text or "")

    def set_transfer_summary(self, text: str) -> None:
        self.transfer_summary.setText(text or "")

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
