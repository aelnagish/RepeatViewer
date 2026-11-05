"""Interactive single-sequence map: scroll, zoom to nucleotide, feature tracks."""

from __future__ import annotations

import logging
import math
import platform
import re
from bisect import bisect_left, bisect_right
from pathlib import Path

from PySide6.QtCore import (
    QEvent,
    QPoint,
    QPointF,
    QRect,
    QRectF,
    QSettings,
    QSize,
    Qt,
    QThread,
    QTimer,
    Signal,
    Slot,
)
from PySide6.QtGui import (
    QColor,
    QContextMenuEvent,
    QCursor,
    QFont,
    QFontMetrics,
    QGuiApplication,
    QKeyEvent,
    QKeySequence,
    QMouseEvent,
    QNativeGestureEvent,
    QPainter,
    QPainterPath,
    QImage,
    QPalette,
    QPen,
    QPixmap,
    QPolygonF,
    QRegion,
    QShortcut,
    QPageSize,
    QPdfWriter,
    QWheelEvent,
)
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QColorDialog,
    QComboBox,
    QCompleter,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMessageBox,
    QMenu,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QScrollBar,
    QSizePolicy,
    QSlider,
    QStackedWidget,
    QStyle,
    QTextEdit,
    QToolButton,
    QToolTip,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.spin_boxes import (
    NoWheelDoubleSpinBox as QDoubleSpinBox,
    NoWheelSlider,
    NoWheelSpinBox as QSpinBox,
)
from repeatviewer.widgets.help_button import field_with_help, make_help_button

# Match TOOL_PARAMS_FIELD_WIDTH (BLAST / Align Settings).
_SETTINGS_FIELD_WIDTH = 260
# Shared label column so LTR / TIR / Tandem / ORF forms line up vertically.
_SETTINGS_LABEL_WIDTH = 210
# Label + field + ``?`` help + margins (long Tandem/LTR labels need room).
_SETTINGS_DIALOG_MIN_WIDTH = 620

# Shift+scroll / trackpad pinch zoom (same factor as MSA / multi-sequence).
_MAP_ZOOM_STEP = 1.08
_MAP_PINCH_DELTA_PER_STEP = 0.04

from repeatviewer.backends.auto_annotation import (
    AutoAnnotationMode,
    AutoAnnotationResult,
    MultiAutoAnnotationResult,
    merge_auto_annotation_features,
    run_auto_annotation,
    run_auto_annotations,
)
from repeatviewer.backends.auto_annotation_settings import (
    AutoAnnotationSettings,
    DEFAULT_MAX_TERMINAL_PAIRS,
    DEFAULT_PBS_MAX_DISTANCE,
    DEFAULT_PBS_MIN_COVERAGE,
    DEFAULT_PBS_MIN_IDENTITY,
    DEFAULT_PBS_MIN_LENGTH,
    DEFAULT_PBS_ORIENTATION,
    DEFAULT_PBS_WINDOW,
    DEFAULT_TIR_MAX_LENGTH,
    DEFAULT_TIR_MIN_LENGTH,
    DEFAULT_TIR_MIN_PAIR_DISTANCE,
    DEFAULT_TIR_SIMILARITY_THRESHOLD,
    MAX_TERMINAL_PAIRS_LIMIT,
    OrfSettings,
    PbsSettings,
    TerminalRepeatSettings,
    TandemRepeatSettings,
    load_auto_annotation_settings,
    save_auto_annotation_settings,
)
from repeatviewer.backends.pbs_detect import default_motif_fasta
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.path_browse import PathBrowseRow
from repeatviewer.branding import APP_DISPLAY_NAME
from repeatviewer.backends.job_cancel import JobCancellation
from repeatviewer.backends.annotation_colors import (
    DEFAULT_ADD_FEATURE_BUTTONS,
    feature_type_color,
    load_add_feature_buttons,
    similarity_database_color,
)
from repeatviewer.backends.database_prepare import database_search_readiness
from repeatviewer.backends.dotplot import reverse_complement
from repeatviewer.backends.sequence_features import (
    SequenceFeature,
    SequenceMapData,
    SingleSequenceSummary,
    REPEATVIEWER_INTERNAL_SOURCE,
    annotation_visibility_key,
    expand_terminal_pair_indices,
    features_gff_path,
    format_single_sequence_summary_parts,
    format_single_sequence_summary_text,
    load_sequence_map_data,
    pack_feature_display_rows,
    remount_parent_ids,
    summarize_loaded_sequence,
    resolve_features_gff_path,
    save_annotations,
    save_annotations_multi,
    summarize_single_sequence_fasta,
)
from repeatviewer.database_settings import load_similarity_annotation_databases
from repeatviewer.backends.similarity_annotation_settings import (
    SimilarityAnnotationSettings,
    load_similarity_annotation_settings,
    save_similarity_annotation_settings,
    similarity_feature_track_label,
)
from repeatviewer.job_options import (
    configure_report_evalue_spinbox,
    load_tsd_search_options,
    save_tsd_search_options,
    TSD_BOUNDARY_REPEATVIEWER,
    TSD_BOUNDARY_SOURCES,
    TSD_BOUNDARY_USER,
    TsdSearchOptions,
)
from repeatviewer.backends.orf_detect import orfipy_available
from repeatviewer.backends.similarity_annotation import (
    SimilarityAnnotationResult,
    SimilarityDatabaseRun,
    load_use_system_blast_from_settings,
    run_similarity_annotations,
    similarity_search_needs_orfs,
)
from repeatviewer.project import UserDatabase
from repeatviewer.workers.base import Worker
from repeatviewer.widgets.annotation_colors_dialog import (
    AddFeatureButtonsDialog,
    AnnotationElementColorsDialog,
)
from repeatviewer.widgets.sequence_annotation_rail import (
    SequenceAnnotationRail,
    SimilarityDatabaseCheck,
)
from repeatviewer.widgets.preview_controls import (
    ANNOTATION_FIT_BUTTON_RENDERED_WIDTH,
    ANNOTATION_SELECT_BUTTON_RENDERED_WIDTH,
    ANNOTATION_ZOOM_SLIDER_WIDTH,
    ZOOM_SLIDER_HEIGHT,
    apply_annotation_zoom_slider_style,
)
from repeatviewer.widgets.preview_zoom import (
    load_show_sequence_summary,
    nucleotide_base_color_map,
)

def _annotation_panel_separator(parent: QWidget | None = None) -> QFrame:
    line = QFrame(parent)
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Plain)
    line.setLineWidth(1)
    line.setFixedHeight(1)
    return line


_SPLIT_TRACKS_BY_SOURCE_KEY = "annotation/split_tracks_by_source"
_SHOW_ONLINE_BLAST_KEY = "annotation/show_online_blast"


def _load_split_tracks_by_source() -> bool:
    value = QSettings().value(_SPLIT_TRACKS_BY_SOURCE_KEY, True)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in ("1", "true", "yes", "on")
    return bool(value)


def _save_split_tracks_by_source(enabled: bool) -> None:
    QSettings().setValue(_SPLIT_TRACKS_BY_SOURCE_KEY, bool(enabled))


def _load_show_online_blast() -> bool:
    """Online BLAST overlay is retired; always hide legacy track features."""
    return False


def _save_show_online_blast(enabled: bool) -> None:
    del enabled
    QSettings().setValue(_SHOW_ONLINE_BLAST_KEY, False)


def filter_features_for_online_blast_display(
    features: list[SequenceFeature],
    *,
    show_online_blast: bool,
) -> tuple[list[SequenceFeature], list[int]]:
    """Return (display_features, full_index_for_each_display_row)."""
    from repeatviewer.backends.ncbi_blast_online import ONLINE_BLAST_SOURCE

    if show_online_blast:
        return list(features), list(range(len(features)))
    display: list[SequenceFeature] = []
    mapping: list[int] = []
    for index, feat in enumerate(features):
        if str(feat.source or "") == ONLINE_BLAST_SOURCE:
            continue
        display.append(feat)
        mapping.append(index)
    return display, mapping


_TRANSFER_PREVIEW_SOURCE = "transfer_preview"


def _transfer_preview_alpha(
    identity: float,
    coverage: float,
    min_identity: float,
    min_coverage: float,
) -> int | None:
    """Alpha for a preview feature, or None when below thresholds."""
    if identity < min_identity or coverage < min_coverage:
        return None
    quality = min(float(identity), float(coverage))
    floor = max(float(min_identity), float(min_coverage))
    if quality <= floor:
        return 72
    span = max(1e-6, 1.0 - floor)
    frac = min(1.0, (quality - floor) / span)
    return int(72 + frac * (255 - 72))


def _configure_annotation_controls_panel(panel: QWidget) -> None:
    """Tag annotation controls panel for layout hooks."""
    panel.setObjectName("annotationControlsPanel")
    panel.setAutoFillBackground(True)


_ADD_FEATURE_BUTTONS: tuple[tuple[str, str], ...] = DEFAULT_ADD_FEATURE_BUTTONS
_AUTO_ANNOTATION_MODES: tuple[tuple[str, AutoAnnotationMode], ...] = (
    ("LTR", AutoAnnotationMode.LTR),
    ("TIR", AutoAnnotationMode.TIR),
    ("Tandem", AutoAnnotationMode.TANDEM),
    ("ORF", AutoAnnotationMode.ORF),
    ("TSD", AutoAnnotationMode.TSD),
    ("PBS", AutoAnnotationMode.PBS),
)
_AUTO_MODE_LABEL = {
    AutoAnnotationMode.LTR: "LTR",
    AutoAnnotationMode.TIR: "TIR",
    AutoAnnotationMode.TANDEM: "Tandem",
    AutoAnnotationMode.ORF: "ORF",
    AutoAnnotationMode.TSD: "TSD",
    AutoAnnotationMode.PBS: "PBS",
}
# Shared column indices for Add feature / Auto annotation alignment.
_ANNOTATION_FEATURE_COLUMN = {
    "LTR": 0,
    "TIR": 1,
    "ORF": 2,
    "Tandem": 3,
    "PBS": 4,
    "PPT": 5,
    "Others": 6,
}
_ANNOTATION_FEATURE_COLUMN_COUNT = 7
_AUTO_ANNOTATION_COLUMN_COUNT = 4
_AUTO_CHECK_COLUMN = {
    AutoAnnotationMode.LTR: _ANNOTATION_FEATURE_COLUMN["LTR"],
    AutoAnnotationMode.TIR: _ANNOTATION_FEATURE_COLUMN["TIR"],
    AutoAnnotationMode.TANDEM: _ANNOTATION_FEATURE_COLUMN["Tandem"],
    AutoAnnotationMode.ORF: _ANNOTATION_FEATURE_COLUMN["ORF"],
}
_ANNOTATION_ROW_LABEL_WIDTH = 118
_ANNOTATION_RUN_SETTINGS_WIDTH = 128
_CONTROL_BAND_HEIGHT = 30
_AUTO_ANN_CHECK_EXTRA_SPACING = 4
_MIN_FEATURE_COLUMN_WIDTH = 54
_SIMILARITY_DB_CHIP_SPACING = 4
_MARQUEE_DRAG_THRESHOLD = 3.0
_FEATURE_EDGE_HIT_PX = 6.0
_ANNOT_TRACK_VPAD = 8
_ANNOT_FEATURE_HEIGHT = 18
logger = logging.getLogger(__name__)
_IMAGE_MARKER_SIZE = 12  # cytogenetic-image dot on annotation features
_ANNOT_MIN_VIEWPORT_HEIGHT = 56
_OVERVIEW_SYNC_MS = 16
_RULER_HEIGHT = 20
# Room for number above a short tick centered on each nucleotide.
_VIEWPORT_RULER_HEIGHT = 18
# Ignore tiny moves on the scale (click = select; drag = select range).
_RULER_RUBBER_MIN_PX = 4
_OVERVIEW_HEIGHT = 16
# Bar-mode (zoomed-out) strip height. Letter mode sizes from the Display font.
_SEQUENCE_ROW = 14
_TRACK_HEIGHT = 12
_TRACK_GAP = 0
_SELECTION_HIGHLIGHT_MARGIN = 2
_FEATURE_UNDO_LIMIT = 32
# Flush sequence to the left edge so bp 1 lines up with the ruler (no leading gap).
_LEFT_MARGIN = 0
_RIGHT_MARGIN = 4
# Annotation map "paper" — always light so track colors stay readable in dark mode.
# Only the zoom/Position toolbar follows the app chrome color.
_MAP_PAPER = QColor("#fafafa")
_ANNOT_PAPER = QColor("#f5f5f5")
_MAP_RAIL_DIVIDER = QColor("#9a9a9a")


def _draw_selection_side_guides(painter: QPainter, rect: QRectF) -> None:
    """Dashed left/right edges only — no top/bottom rail on range selections."""
    if rect.width() < 0.5 or rect.height() < 0.5:
        return
    pen = QPen(_SELECTION_BORDER, 2, Qt.PenStyle.DashLine)
    pen.setCosmetic(True)
    painter.setPen(pen)
    painter.setBrush(Qt.BrushStyle.NoBrush)
    x1 = rect.left() + 0.5
    x2 = rect.right() - 0.5
    y1 = rect.top()
    y2 = rect.bottom()
    painter.drawLine(QPointF(x1, y1), QPointF(x1, y2))
    painter.drawLine(QPointF(x2, y1), QPointF(x2, y2))
    painter.setPen(Qt.PenStyle.NoPen)


def _app_chrome_background() -> QColor:
    """Live application window color (for the Annotation toolbar strip only)."""
    app = QApplication.instance()
    if app is None:
        return QColor(_MAP_PAPER)
    return QColor(app.palette().color(QPalette.ColorRole.Window))


def _apply_map_paper_background(widget: QWidget, color: QColor | None = None) -> None:
    """Force a fixed fill so dark theme does not bleed under the map or toolbar."""
    fill = color if color is not None else _MAP_PAPER
    widget.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
    widget.setAutoFillBackground(True)
    palette = widget.palette()
    palette.setColor(QPalette.ColorRole.Window, fill)
    palette.setColor(QPalette.ColorRole.Base, fill)
    widget.setPalette(palette)


operation_logger = logging.getLogger("repeatviewer.operations")
_MIN_PPB = 0.05
_FIT_MIN_PPB = 1e-6
_MAX_PPB = 48.0
_INTERACTIVE_SEQUENCE_LEN = 35_000
_LARGE_SEQUENCE_LEN = 150_000
# Hard Annotation display limit (sequence length). Viewport-only letters keep
# long sequences usable; this is a safety gate, not a render window size.
_HUGE_SEQUENCE_LEN = 100_000_000
_MAX_ANNOTATION_DISPLAY_LEN = _HUGE_SEQUENCE_LEN
_LAYOUT_SYNC_MS = 40


def _annotation_too_long_message(length: int = 0) -> str:
    """Status text when Annotation refuses to render a long sequence."""
    del length  # kept for call-site compatibility
    return (
        "This sequence is too long for Annotation tab "
        "(limit 100 Mbp of sequence length)."
    )
_BUCKET_DRAW_PPB = 1.0
_WRAP_MIN_PPB = 4.0
_WRAP_LINE_RULER = 12
_WRAP_LINE_GAP = 6
_NO_COLOR_FG = QColor("#222222")
_NO_COLOR_BG = QColor("#fafafa")
_SELECTION_FILL = QColor(10, 40, 120, 100)
_SELECTION_BORDER = QColor(5, 20, 70, 255)
# Unselected annotation-region veil (not applied on the nucleotide strip).
_SELECTION_DIM = QColor(210, 210, 210, 90)
# Edit > Find hits (the current hit is also selected, so it shows in blue on top).
_FIND_HIT_FILL = QColor(245, 158, 11, 110)
_FIND_HIT_BORDER = QColor(180, 83, 9, 230)
_FEATURE_SELECTION_BORDER = QColor("#0b5cab")
_FEATURE_SELECTION_OUTLINE = QColor("#ffffff")
_FEATURE_SELECTION_OVERLAY = QColor(11, 92, 171, 120)
_MIN_NUC_FONT_PT = 2
_MAX_NUC_FONT_PT = 36
_DEFAULT_NUC_FONT_PT = 6
# Square −/+ chrome matching zoom-slider height; glyph keeps default font size.
_ZOOM_BTN_SIDE = ZOOM_SLIDER_HEIGHT
_ANNOTATION_SIDE_PANEL_WIDTH = 108
# Wide enough for "2,000,000" plus spin arrows without clipping.
_POSITION_SPIN_WIDTH = 108


def _make_zoom_button(text: str, parent: QWidget) -> QToolButton:
    """Square −/+ control aligned with the zoom slider; default font size kept."""
    btn = QToolButton(parent)
    btn.setText(text)
    btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
    btn.setAutoRaise(False)
    btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
    side = int(_ZOOM_BTN_SIDE)
    btn.setFixedSize(side, side)
    btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
    # padding/margin only — min/max width/height inflate size by ~2px on macOS.
    btn.setStyleSheet("QToolButton { padding: 0px; margin: 0px; }")
    return btn


def _configure_position_spin(spin: QSpinBox) -> None:
    """Size Position spins so the Annotation max coordinate is fully visible."""
    spin.setFixedWidth(_POSITION_SPIN_WIDTH)
    spin.setAlignment(Qt.AlignmentFlag.AlignRight)
    spin.setGroupSeparatorShown(True)
    spin.setToolTip(
        f"1-based position (up to {_MAX_ANNOTATION_DISPLAY_LEN:,} bp on Annotation)"
    )


def _zoom_ppb_max() -> float:
    """Maximum zoom (pixels per base) for the annotation slider."""
    return max(float(_WRAP_MIN_PPB), float(_MAX_PPB))


def _zoom_ppb_max_for(seq_len: int) -> float:
    """Cap zoom-in for very long sequences to keep scrolling responsive.

    Keep the cap high enough that letter mode (≈ font width px/base) is still
    reachable: letters are painted only for the visible viewport, so deep zoom
    on a multi-Mb sequence stays cheap.
    """
    length = max(0, int(seq_len))
    if length <= _INTERACTIVE_SEQUENCE_LEN:
        return _MAX_PPB
    if length <= _LARGE_SEQUENCE_LEN:
        return 16.0
    if length <= _HUGE_SEQUENCE_LEN:
        return 16.0
    return 8.0


def _annotation_row_label(text: str, parent: QWidget) -> QLabel:
    label = QLabel(text, parent)
    label.setFixedWidth(_ANNOTATION_ROW_LABEL_WIDTH)
    label.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
    return label


def _toolbar_stat_separator(parent: QWidget) -> QLabel:
    label = QLabel("|", parent)
    label.setAlignment(Qt.AlignmentFlag.AlignCenter)
    return label


def _annotation_run_settings_box(
    run_btn: QPushButton | None,
    settings_btn: QPushButton,
    parent: QWidget,
) -> QWidget:
    """Fixed-width Run/Setting column shared by annotation control rows.

    Pass ``run_btn=None`` to keep Setting aligned with rows that have Run.
    """
    box = QWidget(parent)
    box.setObjectName("annotationRunSettingsBox")
    box.setFixedWidth(_ANNOTATION_RUN_SETTINGS_WIDTH)
    layout = QHBoxLayout(box)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(14)
    run_spacer: QWidget | None = None
    if run_btn is not None:
        layout.addWidget(run_btn, 0, Qt.AlignmentFlag.AlignVCenter)
    else:
        run_spacer = QWidget(box)
        run_spacer.setObjectName("annotationRunSettingsSpacer")
        run_spacer.setFixedWidth(1)
        layout.addWidget(run_spacer, 0, Qt.AlignmentFlag.AlignVCenter)
    layout.addWidget(settings_btn, 0, Qt.AlignmentFlag.AlignVCenter)
    layout.addStretch()
    box._run_spacer = run_spacer  # type: ignore[attr-defined]
    return box


def _control_band(parent: QWidget) -> tuple[QWidget, QHBoxLayout]:
    band = QWidget(parent)
    band.setObjectName("annotationControlBand")
    band.setFixedHeight(_CONTROL_BAND_HEIGHT)
    layout = QHBoxLayout(band)
    layout.setContentsMargins(8, 0, 8, 0)
    layout.setSpacing(6)
    layout.setAlignment(Qt.AlignmentFlag.AlignVCenter)
    return band, layout


class _ColumnCell(QWidget):
    """Fixed-width slot with a right-aligned child so label text lines up across rows."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("annotationFeatureColumnCell")
        self.setAutoFillBackground(False)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._layout = QHBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._layout.setSpacing(0)
        self._child: QWidget | None = None

    def set_widget(self, widget: QWidget) -> None:
        if self._child is not None:
            self._layout.removeWidget(self._child)
            self._child.setParent(None)
        self._child = widget
        self._layout.addWidget(
            widget,
            0,
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
        )
        self._layout.addStretch(1)

    def set_column_width(self, width: int) -> None:
        width = max(_MIN_FEATURE_COLUMN_WIDTH, int(width))
        self.setFixedWidth(width)


class _FeatureColumnRow(QWidget):
    """One row of feature controls in fixed-width, right-aligned columns."""

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        column_count: int = _ANNOTATION_FEATURE_COLUMN_COUNT,
    ) -> None:
        super().__init__(parent)
        self.setObjectName("annotationFeatureColumnRow")
        self._column_count = column_count
        self.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        self._grid = QGridLayout(self)
        self._grid.setContentsMargins(0, 0, 0, 0)
        self._grid.setHorizontalSpacing(_AUTO_ANN_CHECK_EXTRA_SPACING)
        self._grid.setVerticalSpacing(0)
        self._cells: dict[int, _ColumnCell] = {}

    def place(self, column: int, widget: QWidget) -> None:
        cell = self._cells.get(column)
        if cell is None:
            cell = _ColumnCell(self)
            self._cells[column] = cell
            self._grid.addWidget(
                cell,
                0,
                column,
                Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
            )
        cell.set_widget(widget)
        self._column_count = max(self._column_count, column + 1)

    def clear(self) -> None:
        for cell in list(self._cells.values()):
            self._grid.removeWidget(cell)
            cell.setParent(None)
            cell.deleteLater()
        self._cells.clear()
        self._column_count = 0
        self.setFixedWidth(0)

    def set_column_count(self, count: int) -> None:
        self._column_count = max(0, int(count))

    def sync_column_widths(self, widths: list[int]) -> None:
        spacing = _AUTO_ANN_CHECK_EXTRA_SPACING
        active_widths = widths[: self._column_count]
        for index, width in enumerate(active_widths):
            self._grid.setColumnMinimumWidth(index, max(0, int(width)))
            self._grid.setColumnStretch(index, 0)
            cell = self._cells.get(index)
            if cell is not None:
                cell.set_column_width(width)
        total = sum(active_widths) + spacing * max(0, len(active_widths) - 1)
        self.setFixedWidth(max(0, total))


def _column_content_width(widget: QWidget) -> int:
    return max(_MIN_FEATURE_COLUMN_WIDTH, widget.sizeHint().width())


class _AddFeatureButtonsRow(QWidget):
    """Horizontally spaced Add feature buttons with uniform gaps."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("addFeatureButtonsRow")
        self.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        self._layout = QHBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._layout.setSpacing(6)
        self._layout.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        )
        self._controls: list[_AddFeatureControl] = []

    def set_controls(self, controls: list[_AddFeatureControl]) -> None:
        while self._layout.count():
            item = self._layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
        self._controls = list(controls)
        for control in self._controls:
            control.setParent(self)
            control.set_label_offset(0)
            self._layout.addWidget(
                control, 0, Qt.AlignmentFlag.AlignVCenter
            )
        self._update_size()

    def clear(self) -> None:
        self.set_controls([])

    def _update_size(self) -> None:
        if not self._controls:
            self.setFixedWidth(0)
            self.setFixedHeight(_CONTROL_BAND_HEIGHT)
            return
        width = sum(control.sizeHint().width() for control in self._controls)
        width += self._layout.spacing() * max(0, len(self._controls) - 1)
        height = max(
            _CONTROL_BAND_HEIGHT,
            max(control.sizeHint().height() for control in self._controls),
        )
        self.setFixedWidth(max(0, width))
        self.setFixedHeight(height)

    def content_width(self) -> int:
        return self.width() if self._controls else 0


class _AddFeatureControl(QWidget):
    """Add-feature button hosted in a shared annotation column."""

    clicked = Signal()

    def __init__(self, label: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("addFeatureControl")
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        self.setAutoFillBackground(False)
        self.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self._offset = QWidget(self)
        self._offset.setAutoFillBackground(False)
        self._offset.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        self._offset.setFixedWidth(0)
        layout.addWidget(self._offset)
        self._btn = QPushButton(label, self)
        self._btn.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        self._btn.clicked.connect(self.clicked.emit)
        layout.addWidget(self._btn)

    def set_label_offset(self, width: int) -> None:
        self._offset.setFixedWidth(max(0, int(width)))

    def label_offset_width(self) -> int:
        return self._offset.width()

    def button(self) -> QPushButton:
        return self._btn


class _SimilarityDbFlowWidget(QWidget):
    """Single-line horizontal row for similarity database checkboxes."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        self.setMinimumWidth(0)
        self._row = QHBoxLayout(self)
        self._row.setContentsMargins(0, 0, 0, 0)
        self._row.setSpacing(_SIMILARITY_DB_CHIP_SPACING)
        self._row.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self._chips: list[QWidget] = []

    def set_chips(self, chips: list[QWidget]) -> None:
        while self._row.count():
            item = self._row.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
        self._chips = list(chips)
        for chip in self._chips:
            chip.setParent(self)
            if isinstance(chip, QLabel):
                h_policy = QSizePolicy.Policy.Expanding
            else:
                h_policy = QSizePolicy.Policy.Maximum
            chip.setSizePolicy(
                h_policy,
                QSizePolicy.Policy.Fixed,
            )
            self._row.addWidget(chip, 0, Qt.AlignmentFlag.AlignVCenter)
        self.setFixedHeight(max(_CONTROL_BAND_HEIGHT, self.sizeHint().height()))
        self.updateGeometry()

    def sizeHint(self) -> QSize:  # noqa: N802
        if not self._chips:
            return QSize(0, _CONTROL_BAND_HEIGHT)
        height = max(chip.sizeHint().height() for chip in self._chips)
        width = sum(chip.sizeHint().width() for chip in self._chips)
        width += _SIMILARITY_DB_CHIP_SPACING * max(0, len(self._chips) - 1)
        return QSize(width, max(_CONTROL_BAND_HEIGHT, height))

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return self.sizeHint()

    def _reflow(self) -> None:
        self.updateGeometry()

    def _available_width(self) -> int:
        return self.width()


def _clone_features(features: list[SequenceFeature]) -> list[SequenceFeature]:
    return [SequenceFeature.from_dict(feat.to_dict()) for feat in features]


def _settings_separator() -> QFrame:
    line = QFrame()
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Sunken)
    return line


def _monospace_font(*, point_size: int | None = None) -> QFont:
    font = QFont()
    font.setFamilies(["Menlo", "Monaco", "Courier New", "monospace"])
    font.setStyleHint(QFont.StyleHint.Monospace)
    font.setFixedPitch(True)
    if point_size is not None:
        font.setPointSize(point_size)
    return font


def _letter_cell_width_for_pt(point_size: int) -> int:
    font = _monospace_font(point_size=point_size)
    fm = QFontMetrics(font)
    return max(1, max(fm.horizontalAdvance(base) for base in "ACGTUNWYF"))


# QColor objects for the Annotation bar — rebuilt when residue palettes change.
_base_qcolor_cache: dict[bool, dict[str, tuple[QColor, QColor]]] = {}


def invalidate_base_qcolor_cache() -> None:
    _base_qcolor_cache.clear()


def _base_colors(base: str, *, is_protein: bool = False) -> tuple[QColor, QColor]:
    """Return (bg, fg) QColors; palette is cached (no QSettings on the paint path)."""
    cache = _base_qcolor_cache.get(is_protein)
    if cache is None:
        if is_protein:
            from repeatviewer.widgets.preview_zoom import amino_acid_color_map

            palette = amino_acid_color_map()
        else:
            palette = nucleotide_base_color_map()
        cache = {
            key: (QColor(bg), QColor(fg)) for key, (bg, fg) in palette.items()
        }
        _base_qcolor_cache[is_protein] = cache
    key = (base or ("X" if is_protein else "N")).upper()
    if key in cache:
        return cache[key]
    fallback = "X" if is_protein else "N"
    if fallback in cache:
        return cache[fallback]
    return _NO_COLOR_BG, _NO_COLOR_FG


def _draw_directed_feature(
    painter: QPainter,
    rect: QRectF,
    color: QColor,
    *,
    forward: bool,
    label: str = "",
    label_point_size: int = 7,
) -> None:
    """Draw a feature bar with a directional arrowhead on the leading end."""
    if rect.width() < 2 or rect.height() < 2:
        return
    arrow_w = min(rect.height() * 0.9, max(4.0, rect.width() * 0.18))
    body = QRectF(rect)
    if forward:
        tip_x = rect.right()
        base_x = rect.right() - arrow_w
        body.setRight(base_x)
        arrow = QPolygonF(
            [
                QPointF(tip_x, rect.center().y()),
                QPointF(base_x, rect.top()),
                QPointF(base_x, rect.bottom()),
            ]
        )
    else:
        tip_x = rect.left()
        base_x = rect.left() + arrow_w
        body.setLeft(base_x)
        arrow = QPolygonF(
            [
                QPointF(tip_x, rect.center().y()),
                QPointF(base_x, rect.top()),
                QPointF(base_x, rect.bottom()),
            ]
        )
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(color)
    if body.width() > 0.5:
        painter.drawRect(body)
    painter.drawPolygon(arrow)
    pt = max(0, min(24, int(label_point_size)))
    if label and pt > 0 and rect.width() > 28:
        painter.setPen(QColor("#ffffff"))
        label_font = QFont(painter.font())
        label_font.setPointSize(pt)
        painter.setFont(label_font)
        text_rect = rect.adjusted(3, 0, -3, 0)
        if not forward:
            text_rect = text_rect.adjusted(arrow_w * 0.4, 0, 0, 0)
        else:
            text_rect = text_rect.adjusted(0, 0, -arrow_w * 0.4, 0)
        painter.drawText(text_rect, Qt.AlignmentFlag.AlignVCenter, label)


# Single-base complement (aligned under the forward strand — not a full RC string).
_COMPLEMENT_BASE = str.maketrans(
    "ACGTURYSWKMBDHVacgturyswkmbdhvnN",
    "TGCAAYRSWMKVHDBtgcaayrswmkvhdbnN",
)


def _complement_base(base: str) -> str:
    if not base:
        return "N"
    return base.translate(_COMPLEMENT_BASE)


def _translate_dna_fragment(dna: str, *, table: int = 1) -> str:
    """Translate a codon-aligned DNA fragment (empty-safe)."""
    from repeatviewer.backends.genetic_codes import translate_dna_codons

    return translate_dna_codons(dna, table=table)


def _translate_dna_frame(sequence: str, frame: int, *, table: int = 1) -> str:
    """Translate one reading frame; ``frame`` in {±1, ±2, ±3}.

    Prefer :func:`_iter_visible_frame_residues` for painting — this full-sequence
    path is only for small strings / tests.
    """
    text = (sequence or "").upper().replace("U", "T")
    if not text or frame == 0:
        return ""
    if frame > 0:
        dna = text[frame - 1 :]
    else:
        dna = reverse_complement(text)[abs(frame) - 1 :]
    return _translate_dna_fragment(dna, table=table)


def _iter_visible_frame_residues(
    sequence: str,
    frame: int,
    visible_start: int,
    visible_end: int,
    *,
    table: int = 1,
):
    """Yield ``(fwd_start0, fwd_end0, aa)`` for residues overlapping the viewport.

    Only the visible window (+ codon padding) is reverse-complemented / translated,
    so enabling translation on multi-megabase sequences stays interactive.
    """
    seq_len = len(sequence or "")
    if seq_len <= 0 or frame == 0 or visible_end <= visible_start:
        return
    vs = max(0, int(visible_start))
    ve = min(seq_len, int(visible_end))
    if ve <= vs:
        return

    if frame > 0:
        offset = frame - 1
        if offset >= seq_len:
            return
        # First codon that can overlap the viewport.
        ai0 = max(0, (vs - offset) // 3)
        dna_start = offset + ai0 * 3
        while dna_start + 3 <= vs and dna_start + 3 <= seq_len:
            dna_start += 3
        if dna_start >= ve or dna_start >= seq_len:
            return
        dna_end = min(seq_len, ve + 2)
        length = dna_end - dna_start
        length -= length % 3
        if length < 3:
            return
        aa = _translate_dna_fragment(
            sequence[dna_start : dna_start + length], table=table
        )
        for ai, residue in enumerate(aa):
            start0 = dna_start + ai * 3
            end0 = start0 + 3
            if end0 <= vs or start0 >= ve:
                continue
            yield start0, end0, residue
        return

    # Reverse frames: translate a short RC window covering the viewport.
    offset = abs(frame) - 1
    rc_lo = max(0, seq_len - ve - 2)
    rc_hi = min(seq_len, seq_len - vs + 2)
    if rc_lo < offset:
        rc_lo = offset
    else:
        rc_lo = offset + ((rc_lo - offset) // 3) * 3
    if rc_hi <= rc_lo:
        return
    fwd_lo = seq_len - rc_hi
    fwd_hi = seq_len - rc_lo
    if fwd_hi <= fwd_lo:
        return
    fragment_rc = reverse_complement(sequence[fwd_lo:fwd_hi])
    length = len(fragment_rc)
    length -= length % 3
    if length < 3:
        return
    aa = _translate_dna_fragment(fragment_rc[:length], table=table)
    for ai, residue in enumerate(aa):
        c0 = rc_lo + ai * 3
        c1 = c0 + 3
        f_start = seq_len - c1
        f_end = seq_len - c0
        if f_end <= vs or f_start >= ve:
            continue
        yield f_start, f_end, residue


def _draw_cell_separator(painter: QPainter, cell: QRectF) -> None:
    """Thin right-edge hairline so adjacent coloured cells read as a grid."""
    painter.setPen(QPen(QColor(0, 0, 0, 60), 1))
    painter.drawLine(cell.topRight(), cell.bottomRight())


def _draw_feature_selection_chrome(painter: QPainter, rect: QRectF) -> None:
    """Bold outline + tint so selected features are easy to see."""
    if rect.width() < 2 or rect.height() < 2:
        return
    frame = rect.adjusted(-2.5, -2.5, 2.5, 2.5)
    painter.setPen(QPen(_FEATURE_SELECTION_OUTLINE, 1))
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawRect(frame)
    painter.setPen(QPen(_FEATURE_SELECTION_BORDER, 2))
    painter.drawRect(rect.adjusted(-1.5, -1.5, 1.5, 1.5))
    painter.setPen(Qt.PenStyle.NoPen)
    painter.setBrush(_FEATURE_SELECTION_OVERLAY)
    painter.drawRect(rect)


def _max_cell_w() -> int:
    return max(1, int(_zoom_ppb_max()))


def _pack_interval_rows(
    segments: list[tuple[int, int, int]],
) -> tuple[list[int], int]:
    """Assign stacked rows to (index, start, end) segments on one sequence line."""
    if not segments:
        return [], 0
    row_for_index = [0] * len(segments)
    row_ends: list[int] = []
    for seg_i, (_idx, start, end) in sorted(
        enumerate(segments),
        key=lambda item: item[1][1],
    ):
        for row, last_end in enumerate(row_ends):
            if start > last_end:
                row_ends[row] = end
                row_for_index[seg_i] = row
                break
        else:
            row_for_index[seg_i] = len(row_ends)
            row_ends.append(end)
    return row_for_index, max(0, len(row_ends))


def _pack_feature_row_group(features: list[SequenceFeature]) -> tuple[list[int], int]:
    """Assign stacked rows for one track, packing overlaps vertically."""
    if not features:
        return [], 0
    row_for_index = [0] * len(features)
    row_ends: list[int] = []
    for idx, feat in sorted(enumerate(features), key=lambda item: item[1].start):
        for row, last_end in enumerate(row_ends):
            if feat.start > last_end:
                row_ends[row] = feat.end
                row_for_index[idx] = row
                break
        else:
            row_for_index[idx] = len(row_ends)
            row_ends.append(feat.end)
    return row_for_index, max(0, len(row_ends))


def _pack_feature_rows(
    features: list[SequenceFeature],
    *,
    split_by_source: bool = True,
) -> tuple[list[int], int]:
    """Assign rows per (source, track) or track, then stack bands vertically."""
    return pack_feature_display_rows(features, split_by_source=split_by_source)


class OverviewWidget(QWidget):
    """Mini-map showing features and the visible viewport."""

    viewportChanged = Signal(float)  # fraction 0..1 of sequence start

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedHeight(_OVERVIEW_HEIGHT)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setMinimumWidth(0)
        self._length = 1
        self._features: list[SequenceFeature] = []
        self._view_start = 0.0
        self._view_fraction = 1.0
        self._selected_indices: set[int] = set()
        self._dragging = False
        self._mouse_grabbed = False

    def _release_drag_mouse(self) -> None:
        if self._mouse_grabbed:
            self.releaseMouse()
            self._mouse_grabbed = False

    def hideEvent(self, event) -> None:  # noqa: N802
        self._dragging = False
        self._release_drag_mouse()
        super().hideEvent(event)

    def set_data(
        self,
        length: int,
        features: list[SequenceFeature],
        *,
        view_start_bp: float,
        view_end_bp: float,
        selected_indices: set[int] | frozenset[int] | None = None,
    ) -> None:
        self._length = max(1, int(length))
        self._features = list(features)
        self._selected_indices = set(selected_indices or [])
        self._view_start = max(0.0, view_start_bp) / self._length
        span = max(1.0, view_end_bp - view_start_bp)
        self._view_fraction = min(1.0, span / self._length)
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor("#ececec"))
        w = self.width()
        h = self.height()
        # Rectangular track (not pill-shaped) — keep annotation overview crisp.
        inner = QRectF(0, 2, max(1.0, float(w)), h - 4)
        painter.setPen(QPen(QColor("#999999")))
        painter.drawRect(inner)

        for idx, feat in enumerate(self._features):
            x1 = inner.left() + inner.width() * (feat.start - 1) / self._length
            x2 = inner.left() + inner.width() * feat.end / self._length
            bar = QRectF(x1, inner.top(), max(2.0, x2 - x1), inner.height())
            color = QColor(feat.color or feature_type_color(feat.feature_type))
            painter.fillRect(bar, color)
            if idx in self._selected_indices:
                _draw_feature_selection_chrome(painter, bar)

        vx = inner.left() + inner.width() * self._view_start
        vw = max(4.0, inner.width() * self._view_fraction)
        painter.setPen(QPen(QColor("#2c3e50"), 2))
        painter.setBrush(QColor(44, 62, 80, 40))
        painter.drawRect(QRectF(vx, inner.top(), vw, inner.height()))
        painter.end()

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._dragging = True
            self.grabMouse()
            self._mouse_grabbed = True
            self._emit_fraction(event.position().x())
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if self._dragging:
            self._emit_fraction(event.position().x())
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._dragging = False
            self._release_drag_mouse()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def _emit_fraction(self, x: float) -> None:
        inner_w = max(1.0, float(self.width()))
        frac = x / inner_w - self._view_fraction / 2
        frac = max(0.0, min(1.0 - self._view_fraction, frac))
        self.viewportChanged.emit(frac)


def _bp_view_x(bp: int, scroll_px: float, ppb: float) -> float:
    """Map 1-based bp to x in the viewport-aligned ruler/annotation strip."""
    return _LEFT_MARGIN + (bp - 1) * ppb - scroll_px


def _horizontal_wheel_delta(event: QWheelEvent, *, scroll_step: int) -> int:
    """Pixel delta for horizontal scrolling (trackpad swipe or Shift+wheel)."""
    pixel = event.pixelDelta()
    if pixel.x() != 0:
        return pixel.x()

    angle = event.angleDelta()
    if angle.x() != 0:
        return int((angle.x() / 120.0) * max(1, scroll_step))
    if event.modifiers() & Qt.KeyboardModifier.ShiftModifier and angle.y() != 0:
        return int((angle.y() / 120.0) * max(1, scroll_step))
    return 0


class SequenceScrollArea(QScrollArea):
    """Scroll area that pans horizontally from trackpad swipe (not Shift+wheel)."""

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802
        # Shift+wheel is zoom (handled by SequenceMapWidget); do not remap to pan.
        mods = event.modifiers() | QGuiApplication.keyboardModifiers()
        if mods & Qt.KeyboardModifier.ShiftModifier:
            super().wheelEvent(event)
            return
        hbar = self.horizontalScrollBar()
        delta_x = _horizontal_wheel_delta(event, scroll_step=hbar.singleStep() * 3)
        if delta_x != 0 and hbar.maximum() > 0:
            hbar.setValue(hbar.value() - delta_x)
            event.accept()
            return
        super().wheelEvent(event)


class ViewportRulerWidget(QWidget):
    """Position scale aligned with the visible sequence (scroll + zoom synced).

    Click → select position/column; double-click → zoom to max; drag → range;
    Shift/right-click → fit.
    """

    # Viewport-local x under the cursor + zoom-in flag (True=in, False=out).
    zoomAtRequested = Signal(float, bool)
    # 1-based inclusive range from a completed scale drag.
    regionSelected = Signal(int, int)
    # Live range while dragging (for map/Position preview).
    regionPreview = Signal(int, int)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedHeight(_VIEWPORT_RULER_HEIGHT)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setCursor(Qt.CursorShape.SizeHorCursor)
        self.setToolTip(
            "Click to select a position · "
            "Double-click to zoom to maximum at this position · "
            "Drag to select a region · "
            "Right-click or Shift+click to fit"
        )
        self._sequence = ""
        self._scroll_px = 0.0
        self._ppb = 1.0
        self._viewport_w = 200
        self._font_point_size = 8
        self._region_sel: tuple[int, int] | None = None
        self._drag_origin_x: float | None = None
        self._drag_current_x: float | None = None
        self._rubber_active = False
        self._press_modifiers = Qt.KeyboardModifier.NoModifier
        self._suppress_click_select = False

    def set_font_point_size(self, point_size: int) -> None:
        self._font_point_size = max(4, min(24, int(point_size)))
        self.update()

    def set_view(
        self,
        *,
        sequence: str,
        scroll_px: float,
        ppb: float,
        viewport_w: int,
    ) -> None:
        self._sequence = sequence.upper()
        self._scroll_px = max(0.0, float(scroll_px))
        self._ppb = max(_FIT_MIN_PPB, float(ppb))
        self._viewport_w = max(1, int(viewport_w))
        self.update()

    def set_region_selection(self, start: int | None, end: int | None = None) -> None:
        """Mirror the map selection band onto the scale."""
        if start is None:
            if self._region_sel is None:
                return
            self._region_sel = None
            self.update()
            return
        seq_len = max(1, len(self._sequence))
        a = max(1, min(seq_len, int(start)))
        b = max(1, min(seq_len, int(end if end is not None else start)))
        if b < a:
            a, b = b, a
        new = (a, b)
        if new == self._region_sel:
            return
        self._region_sel = new
        self.update()

    def _x_to_bp(self, x: float) -> int:
        seq_len = max(1, len(self._sequence))
        ppb = max(_FIT_MIN_PPB, float(self._ppb))
        bp = int((float(x) + float(self._scroll_px) - float(_LEFT_MARGIN)) / ppb) + 1
        return max(1, min(seq_len, bp))

    def _bp_to_x(self, bp: float) -> float:
        ppb = max(_FIT_MIN_PPB, float(self._ppb))
        return float(_LEFT_MARGIN) + (float(bp) - 1.0) * ppb - float(self._scroll_px)

    def _range_from_x(self, x0: float, x1: float) -> tuple[int, int]:
        a = self._x_to_bp(x0)
        b = self._x_to_bp(x1)
        if b < a:
            a, b = b, a
        return a, b

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._suppress_click_select = False
            self._drag_origin_x = float(event.position().x())
            self._drag_current_x = self._drag_origin_x
            self._rubber_active = False
            self._press_modifiers = event.modifiers()
            self.grabMouse()
            event.accept()
            return
        if event.button() == Qt.MouseButton.RightButton:
            self.zoomAtRequested.emit(float(event.position().x()), False)
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if (
            self._drag_origin_x is not None
            and event.buttons() & Qt.MouseButton.LeftButton
        ):
            # Prefer live cursor (coalesce backlog); fall back to the event x.
            mapped = float(self.mapFromGlobal(QCursor.pos()).x())
            self._drag_current_x = (
                mapped if self.rect().adjusted(-80, -80, 80, 80).contains(
                    int(mapped), int(self.height() // 2)
                )
                else float(event.position().x())
            )
            if abs(self._drag_current_x - self._drag_origin_x) >= _RULER_RUBBER_MIN_PX:
                if not self._rubber_active:
                    self._region_sel = None
                self._rubber_active = True
                start_i, end_i = self._range_from_x(
                    self._drag_origin_x, self._drag_current_x
                )
                self.regionPreview.emit(start_i, end_i)
                self.update()
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseDoubleClickEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._suppress_click_select = True
            self._drag_origin_x = None
            self._drag_current_x = None
            self._rubber_active = False
            if QWidget.mouseGrabber() is self:
                self.releaseMouse()
            self.zoomAtRequested.emit(float(event.position().x()), True)
            event.accept()
            return
        super().mouseDoubleClickEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton and self._drag_origin_x is not None:
            if QWidget.mouseGrabber() is self:
                self.releaseMouse()
            x0 = self._drag_origin_x
            x1 = float(event.position().x())
            active = self._rubber_active
            mods = self._press_modifiers
            self._drag_origin_x = None
            self._drag_current_x = None
            self._rubber_active = False
            if self._suppress_click_select:
                self._suppress_click_select = False
                self.update()
                event.accept()
                return
            if active and abs(x1 - x0) >= _RULER_RUBBER_MIN_PX:
                start_i, end_i = self._range_from_x(x0, x1)
                self._region_sel = (start_i, end_i)
                self.regionSelected.emit(start_i, end_i)
            elif mods & Qt.KeyboardModifier.ShiftModifier:
                # Shift+click → fit (zoom out) at this position.
                self.zoomAtRequested.emit(x1, False)
            else:
                # Click (no drag) → select the position/column under the cursor.
                start_i, end_i = self._range_from_x(x1, x1)
                self._region_sel = (start_i, end_i)
                self.regionSelected.emit(start_i, end_i)
            self.update()
            event.accept()
            return
        if (
            event.button() == Qt.MouseButton.LeftButton
            and self._suppress_click_select
        ):
            self._suppress_click_select = False
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def paintEvent(self, event) -> None:  # noqa: N802
        """Short ticks centered on nucleotides; numbers centered above the ticks."""
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.fillRect(self.rect(), _MAP_PAPER)
        seq_len = len(self._sequence)
        if seq_len <= 0:
            painter.end()
            return
        ppb = max(_FIT_MIN_PPB, float(self._ppb))
        scroll = float(self._scroll_px)
        vw = max(1, int(self._viewport_w))
        # Selection / live rubber band under the ticks.
        band: tuple[int, int] | None = self._region_sel
        if (
            self._rubber_active
            and self._drag_origin_x is not None
            and self._drag_current_x is not None
        ):
            band = self._range_from_x(self._drag_origin_x, self._drag_current_x)
        if band is not None:
            x0 = self._bp_to_x(band[0])
            x1 = self._bp_to_x(band[1] + 1)
            if x1 < x0:
                x0, x1 = x1, x0
            painter.fillRect(
                QRectF(x0, 0.0, max(1.0, x1 - x0), float(self.height())),
                _SELECTION_FILL,
            )
        ruler_font = _monospace_font(
            point_size=max(4, min(24, int(self._font_point_size)))
        )
        painter.setFont(ruler_font)
        fm = QFontMetrics(ruler_font)
        fg = QColor("#666666")

        origin = float(_LEFT_MARGIN)
        start_bp = max(1.0, (scroll - origin) / ppb + 1.0)
        end_bp = min(float(seq_len), (scroll - origin + float(vw)) / ppb + 1.0)
        span = max(1.0, end_bp - start_bp)
        sample = f"{int(end_bp)}"
        min_label_px = fm.horizontalAdvance(sample) + 12
        step = max(1, SequenceMapCanvas._nice_tick_step(max(1, int(span)), ppb, fm))
        est = span * float(min_label_px) / max(1.0, float(vw))
        if step < est:
            step = max(
                step,
                SequenceMapCanvas._nice_tick_step(max(1, int(est * 4)), ppb, fm),
            )

        # Number (top) → short tick → nucleotide.
        tick_len = 3
        tick_bottom = self.height() - 2
        tick_top = tick_bottom - tick_len
        label_baseline = max(fm.ascent() + 1, tick_top - 2)

        def _center_x(pos: float) -> float:
            return origin + (pos - 0.5) * ppb - scroll

        label_positions: list[int] = []
        if start_bp <= 1.0 <= end_bp:
            label_positions.append(1)
        first = int(math.ceil(start_bp / step) * step)
        if first < step:
            first = step
        for pos in range(first, int(end_bp) + 1, step):
            if 1 <= pos <= seq_len and pos != 1:
                label_positions.append(pos)

        tick_start = max(1, int(math.floor(start_bp)))
        tick_end = min(seq_len, int(math.ceil(end_bp)))
        show_every_base = ppb >= 3.5
        tick_positions = (
            range(tick_start, tick_end + 1)
            if show_every_base
            else label_positions
        )
        painter.setPen(QPen(fg, 1))
        for pos in tick_positions:
            cx = _center_x(float(pos))
            if cx < -2 or cx > vw + 2:
                continue
            painter.drawLine(int(cx), tick_top, int(cx), tick_bottom)

        last_label_right = -1.0e9
        for pos in label_positions:
            cx = _center_x(float(pos))
            if cx < -2 or cx > vw + 2:
                continue
            text = str(pos)
            text_w = fm.horizontalAdvance(text)
            label_left = cx - text_w / 2
            if label_left < 0:
                label_left = 0
            label_right = label_left + text_w
            if label_right > vw:
                label_left = vw - text_w
                label_right = label_left + text_w
            if label_left >= last_label_right + 4 and label_right <= vw + 0.5:
                painter.drawText(int(label_left), int(label_baseline), text)
                last_label_right = label_right
        painter.end()


def _is_delete_key(key: int) -> bool:
    return key in (int(Qt.Key.Key_Delete), int(Qt.Key.Key_Backspace))


def _safe_filename_part(text: str) -> str:
    cleaned = "".join(
        ch if ch.isalnum() or ch in "-_." else "_" for ch in str(text or "").strip()
    )
    return cleaned.strip("._") or ""


def _feature_export_tag(feat: SequenceFeature) -> str:
    """Type/database tag used in extracted FASTA names (``seq_LTR.fa``)."""
    if feat.feature_type == "domain":
        database, sep, _rest = str(feat.name or "").partition(": ")
        if sep and database.strip():
            return _safe_filename_part(database.strip()) or "domain"
        return "domain"
    for label, feature_type in DEFAULT_ADD_FEATURE_BUTTONS:
        if feature_type == feat.feature_type:
            return _safe_filename_part(label) or "feature"
    return _safe_filename_part(feat.feature_type.replace(" ", "_")) or "feature"


def _feature_export_stem(sequence_id: str, feat: SequenceFeature) -> str:
    base = _safe_filename_part(sequence_id) or "sequence"
    return f"{base}_{_feature_export_tag(feat)}"


def _unique_feature_fa_path(directory: Path, sequence_id: str, feat: SequenceFeature) -> Path:
    """``name_LTR.fa``, then ``name_LTR2.fa``, ``name_LTR3.fa``, …"""
    stem = _feature_export_stem(sequence_id, feat)
    candidate = directory / f"{stem}.fa"
    if not candidate.exists():
        return candidate
    index = 2
    while True:
        candidate = directory / f"{stem}{index}.fa"
        if not candidate.exists():
            return candidate
        index += 1


def _feature_sequence_fragment(sequence: str, feat: SequenceFeature) -> str:
    lo = max(1, min(int(feat.start), int(feat.end)))
    hi = min(len(sequence), max(int(feat.start), int(feat.end)))
    if hi < lo:
        return ""
    fragment = sequence[lo - 1 : hi]
    if not feat.is_forward():
        fragment = reverse_complement(fragment)
    return fragment


def _feature_with_track(feat: SequenceFeature, track: int) -> SequenceFeature:
    return SequenceFeature(
        start=feat.start,
        end=feat.end,
        name=feat.name,
        feature_type=feat.feature_type,
        color=feat.color,
        strand=feat.strand,
        track=max(0, int(track)),
        parent=feat.parent,
        source=feat.source,
        pident=feat.pident,
        subject_id=feat.subject_id,
        coverage=feat.coverage,
        source_start=feat.source_start,
        source_end=feat.source_end,
    )


def _feature_with_bounds(
    feat: SequenceFeature, start: int, end: int
) -> SequenceFeature:
    a, b = int(start), int(end)
    if b < a:
        a, b = b, a
    return SequenceFeature(
        start=a,
        end=b,
        name=feat.name,
        feature_type=feat.feature_type,
        color=feat.color,
        strand=feat.strand,
        track=feat.track,
        parent=feat.parent,
        source=feat.source,
        pident=feat.pident,
        subject_id=feat.subject_id,
        coverage=feat.coverage,
        source_start=feat.source_start,
        source_end=feat.source_end,
    )


def _feature_with_color(feat: SequenceFeature, color: str) -> SequenceFeature:
    return SequenceFeature(
        start=feat.start,
        end=feat.end,
        name=feat.name,
        feature_type=feat.feature_type,
        color=color,
        strand=feat.strand,
        track=feat.track,
        parent=feat.parent,
        source=feat.source,
        pident=feat.pident,
        subject_id=feat.subject_id,
        coverage=feat.coverage,
        source_start=feat.source_start,
        source_end=feat.source_end,
    )


def _features_overlap(a: SequenceFeature, b: SequenceFeature) -> bool:
    return a.start <= b.end and b.start <= a.end


def _normalize_feature_track_ids(features: list[SequenceFeature]) -> list[SequenceFeature]:
    """Keep relative track order; remap to contiguous non-negative ids."""
    ordered = sorted({feat.track for feat in features})
    if not ordered:
        return list(features)
    remap = {old: idx for idx, old in enumerate(ordered)}
    return [_feature_with_track(feat, remap[feat.track]) for feat in features]


def reassign_tracks_after_row_join(
    features: list[SequenceFeature],
    moving: set[int],
    join_track: int,
    local_row: int,
) -> list[SequenceFeature] | None:
    """Join ``moving`` onto a packed row of ``join_track``.

    Non-overlapping neighbours stay on the joined level. Residents on that row
    that overlap the moved feature(s) are displaced to the next track only
    (e.g. HHCC moves down when INT joins PROT/RT/RNH — not the whole row).
    Other packed rows of ``join_track`` keep their relative order as separate
    tracks. Returns None when nothing would change.
    """
    if not moving or not features:
        return None
    moving = {idx for idx in moving if 0 <= idx < len(features)}
    if not moving:
        return None

    residents = [
        idx
        for idx, feat in enumerate(features)
        if feat.track == join_track and idx not in moving
    ]
    group = [features[idx] for idx in residents]
    group_rows, group_count = _pack_feature_row_group(group)
    group_count = max(1, int(group_count) or 1)
    target_row = max(0, min(group_count - 1, int(local_row)))

    row_members: list[list[int]] = [[] for _ in range(group_count)]
    for local_i, idx in enumerate(residents):
        row = group_rows[local_i] if local_i < len(group_rows) else 0
        row = max(0, min(group_count - 1, int(row)))
        row_members[row].append(idx)

    moving_feats = [features[idx] for idx in sorted(moving)]

    def overlaps_moving(idx: int) -> bool:
        feat = features[idx]
        return any(_features_overlap(feat, other) for other in moving_feats)

    keep_on_target = [idx for idx in row_members[target_row] if not overlaps_moving(idx)]
    displace = [idx for idx in row_members[target_row] if overlaps_moving(idx)]

    ordered_tracks = sorted({feat.track for feat in features})
    new_track_for_index = [feat.track for feat in features]
    next_track = 0
    target_new_track = 0
    for track in ordered_tracks:
        if track != join_track:
            assigned = False
            for idx, feat in enumerate(features):
                if feat.track == track and idx not in moving:
                    new_track_for_index[idx] = next_track
                    assigned = True
            if assigned:
                next_track += 1
            continue
        for row in range(group_count):
            if row == target_row:
                target_new_track = next_track
                for idx in keep_on_target:
                    new_track_for_index[idx] = next_track
                next_track += 1
                if displace:
                    for idx in displace:
                        new_track_for_index[idx] = next_track
                    next_track += 1
            elif row_members[row]:
                for idx in row_members[row]:
                    new_track_for_index[idx] = next_track
                next_track += 1
    for idx in moving:
        new_track_for_index[idx] = target_new_track

    updated = [
        _feature_with_track(feat, new_track_for_index[idx])
        for idx, feat in enumerate(features)
    ]
    normalized = _normalize_feature_track_ids(updated)
    current = _normalize_feature_track_ids(list(features))
    if [feat.track for feat in normalized] == [feat.track for feat in current]:
        return None
    return normalized


class AnnotationTrackWidget(QWidget):
    """Feature annotations below the sequence; stacked rows for overlaps."""

    selectionChanged = Signal(object)  # frozenset[int]
    clearRangeSelectionRequested = Signal()
    deleteRequested = Signal()
    copySequenceRequested = Signal()
    extractSequenceRequested = Signal()
    modifyRequested = Signal()
    changeColorRequested = Signal()
    tracksMoved = Signal(object)  # list[SequenceFeature]
    featuresResized = Signal(object)  # list[SequenceFeature]
    attachImageRequested = Signal()  # selected feature
    openImagesRequested = Signal(int)  # feature index with images
    imagesDropped = Signal(int, list)  # feature index, image paths

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMouseTracking(True)
        self.setAcceptDrops(True)
        # Returns indices of features that carry cytogenetic images (see
        # SequenceMapPanel._marked_feature_indices); evaluated once per paint.
        self._image_marker_lookup = None
        self._marked_indices: set[int] = set()
        self._drop_hover_idx: int | None = None
        self._features: list[SequenceFeature] = []
        self._feature_rows: list[int] = []
        self._scroll_px = 0.0
        self._ppb = 1.0
        self._viewport_w = 200
        self._seq_length = 1
        self._selected_indices: set[int] = set()
        self._row_count = 0
        self._marquee_active = False
        self._marquee_origin: QPointF | None = None
        self._marquee_current: QPointF | None = None
        self._press_pos: QPointF | None = None
        self._press_hit: int | None = None
        self._press_edge: str | None = None  # "start" | "end"
        self._marquee_additive_base: set[int] = set()
        self._mouse_grabbed = False
        self._show_similarity_database_name = True
        self._move_active = False
        self._move_indices: set[int] = set()
        self._move_preview_track: int | None = None
        self._move_ghost_y: float | None = None
        self._resize_active = False
        self._resize_idx: int | None = None
        self._resize_edge: str | None = None
        self._resize_orig_start = 0
        self._resize_orig_end = 0
        self._feature_bar_height = float(_ANNOT_FEATURE_HEIGHT)
        self._track_gap = float(_TRACK_GAP)
        self._label_point_size = 7
        self._hover_idx: int | None = None
        self._feature_alphas: list[int] = []
        self._range_selection: tuple[int, int] | None = None
        self.setFixedHeight(0)

    def set_range_selection(self, start: int | None, end: int | None = None) -> None:
        """Highlight a bp range vertically across the annotation track."""
        if start is None or end is None:
            new: tuple[int, int] | None = None
        else:
            a, b = int(start), int(end)
            if a > b:
                a, b = b, a
            new = (a, b)
        if new == self._range_selection:
            return
        self._range_selection = new
        self.update()

    def set_feature_bar_height(self, height: int) -> None:
        self._feature_bar_height = float(max(2, min(40, int(height))))
        self._refresh_geometry_from_layout()

    def set_track_gap(self, gap: int) -> None:
        self._track_gap = float(max(0, min(40, int(gap))))
        self._refresh_geometry_from_layout()

    def _refresh_geometry_from_layout(self) -> None:
        if self._features:
            self.set_view(
                features=self._features,
                feature_rows=self._feature_rows,
                scroll_px=self._scroll_px,
                ppb=self._ppb,
                viewport_w=self._viewport_w,
                selected_indices=self._selected_indices,
                seq_length=self._seq_length,
            )
        else:
            self.update()

    def set_label_point_size(self, point_size: int) -> None:
        self._label_point_size = max(0, min(24, int(point_size)))
        self.update()

    def set_show_similarity_database_name(self, show: bool) -> None:
        show = bool(show)
        if show == self._show_similarity_database_name:
            return
        self._show_similarity_database_name = show
        self.update()

    def _feature_track_label(self, feat: SequenceFeature, rect_width: float, min_width: float) -> str:
        if self._label_point_size <= 0 or rect_width <= min_width:
            return ""
        base = feat.name or feat.feature_type
        return similarity_feature_track_label(
            base,
            feat.feature_type,
            show_database_name=self._show_similarity_database_name,
        )

    def _release_drag_mouse(self) -> None:
        if self._mouse_grabbed:
            self.releaseMouse()
            self._mouse_grabbed = False
        self.unsetCursor()

    def hideEvent(self, event) -> None:  # noqa: N802
        self._press_pos = None
        self._press_edge = None
        self._marquee_active = False
        self._move_active = False
        self._move_indices = set()
        self._move_preview_track = None
        self._move_ghost_y = None
        self._resize_active = False
        self._resize_idx = None
        self._resize_edge = None
        self._hover_idx = None
        self._release_drag_mouse()
        super().hideEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802
        if (
            self._hover_idx is not None
            and not self._move_active
            and not self._marquee_active
            and not self._resize_active
        ):
            self._hover_idx = None
            self.unsetCursor()
            self.update()
        super().leaveEvent(event)

    def _row_height(self) -> float:
        return float(self._feature_bar_height + self._track_gap)

    def _feature_full_label(self, feat: SequenceFeature) -> str:
        base = feat.name or feat.feature_type
        return similarity_feature_track_label(
            base,
            feat.feature_type,
            show_database_name=self._show_similarity_database_name,
        )

    def _label_overflows(self, feat: SequenceFeature, rect: QRectF) -> bool:
        label = self._feature_full_label(feat)
        if not label or self._label_point_size <= 0:
            return False
        if rect.width() <= 28:
            return True
        label_font = QFont(self.font())
        label_font.setPointSize(self._label_point_size)
        fm = QFontMetrics(label_font)
        available = max(0.0, rect.width() - 10.0)
        return fm.horizontalAdvance(label) > available

    def _track_band_ranges(self) -> list[tuple[int, float, float]]:
        """Return (track_id, y_top, y_bottom) for each track in display order."""
        if not self._features:
            return []
        tracks = sorted({feat.track for feat in self._features})
        bands: list[tuple[int, float, float]] = []
        row_offset = 0
        row_h = self._row_height()
        for track in tracks:
            group = [feat for feat in self._features if feat.track == track]
            _, group_count = _pack_feature_row_group(group)
            group_count = max(1, group_count)
            y0 = _ANNOT_TRACK_VPAD + row_offset * row_h
            y1 = y0 + group_count * row_h
            bands.append((track, y0, y1))
            row_offset += group_count
        return bands

    def _target_track_at_y(self, y: float) -> int:
        """Map a vertical cursor position to a durable feature track id."""
        bands = self._track_band_ranges()
        if not bands:
            return 0
        if y < bands[0][1]:
            return bands[0][0]
        for track, y0, y1 in bands:
            if y0 <= y < y1:
                return track
        return bands[-1][0]

    def _features_after_track_move(self, y: float) -> list[SequenceFeature] | None:
        """Return updated features after a vertical track move, or None if unchanged.

        Anywhere inside a track band joins that band's visual row under the
        cursor. Only residents that overlap the moved feature(s) are displaced
        downward (HHCC), so non-overlapping neighbours (PROT / RT / RNH) stay
        put. A new empty track is created only when dropping above the first
        band or below the last — not on in-band edges (that used to park INT on
        its own level and shove the whole grey band down).
        """
        moving = {idx for idx in self._move_indices if 0 <= idx < len(self._features)}
        if not moving:
            return None
        bands = self._track_band_ranges()
        if not bands:
            updated = [
                _feature_with_track(feat, 0) if idx in moving else feat
                for idx, feat in enumerate(self._features)
            ]
            return _normalize_feature_track_ids(updated)

        if y < bands[0][1]:
            insert_before = 0
        elif y >= bands[-1][2]:
            insert_before = len(bands)
        else:
            for track, y0, y1 in bands:
                if y0 <= y < y1:
                    row_h = max(1.0, float(self._row_height()))
                    local_row = int((y - y0) // row_h)
                    return reassign_tracks_after_row_join(
                        list(self._features), moving, track, local_row
                    )
            return None

        old_to_new = {
            track: (index if index < insert_before else index + 1)
            for index, (track, _, _) in enumerate(bands)
        }
        new_track = insert_before
        updated = []
        for idx, feat in enumerate(self._features):
            if idx in moving:
                updated.append(_feature_with_track(feat, new_track))
            else:
                updated.append(
                    _feature_with_track(feat, old_to_new.get(feat.track, feat.track))
                )
        return _normalize_feature_track_ids(updated)

    def _join_features_to_track_row(
        self,
        moving: set[int],
        join_track: int,
        local_row: int,
    ) -> list[SequenceFeature] | None:
        return reassign_tracks_after_row_join(
            list(self._features), moving, join_track, local_row
        )

    def _normalize_track_ids(self, features: list[SequenceFeature]) -> list[SequenceFeature]:
        return _normalize_feature_track_ids(features)

    def _can_start_move(self, hit: int | None) -> bool:
        return (
            hit is not None
            and hit in self._selected_indices
            and bool(self._selected_indices)
        )

    def _feature_rect(self, idx: int) -> QRectF:
        feat = self._features[idx]
        row = self._feature_rows[idx] if idx < len(self._feature_rows) else 0
        row_h = self._row_height()
        x1 = _bp_view_x(feat.start, self._scroll_px, self._ppb)
        x2 = _bp_view_x(feat.end + 1, self._scroll_px, self._ppb)
        y = _ANNOT_TRACK_VPAD + row * row_h + 0.5
        return QRectF(
            x1, y, max(2.0, x2 - x1), max(2.0, float(self._feature_bar_height) - 1.0)
        )

    def _feature_row_rect(self, idx: int) -> QRectF:
        feat = self._features[idx]
        row = self._feature_rows[idx] if idx < len(self._feature_rows) else 0
        row_h = self._row_height()
        x1 = _bp_view_x(feat.start, self._scroll_px, self._ppb)
        x2 = _bp_view_x(feat.end + 1, self._scroll_px, self._ppb)
        y = _ANNOT_TRACK_VPAD + row * row_h
        return QRectF(x1, y, max(2.0, x2 - x1), row_h)

    def _additive_modifier(self, modifiers: Qt.KeyboardModifier) -> bool:
        return bool(
            modifiers
            & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.ControlModifier)
        )

    def _hit_test_feature(self, pos: QPointF) -> int | None:
        for idx in range(len(self._features) - 1, -1, -1):
            if self._feature_row_rect(idx).contains(pos):
                return idx
        return None

    def _x_to_bp(self, x: float) -> int:
        if self._ppb <= 0:
            return 1
        bp = int((float(x) + self._scroll_px - _LEFT_MARGIN) / self._ppb) + 1
        return max(1, min(max(1, self._seq_length), bp))

    def _feature_resizable(self, idx: int) -> bool:
        if not (0 <= idx < len(self._features)):
            return False
        if idx < len(self._feature_alphas) and int(self._feature_alphas[idx]) < 255:
            return False  # transfer / ghost preview overlay
        return True

    def _edge_hit(self, pos: QPointF) -> tuple[int, str] | None:
        """Return ``(feature_index, 'start'|'end')`` when near a feature edge."""
        idx = self._hit_test_feature(pos)
        if idx is None or not self._feature_resizable(idx):
            return None
        rect = self._feature_rect(idx)
        threshold = min(
            max(4.0, _FEATURE_EDGE_HIT_PX),
            max(2.0, rect.width() * 0.35),
        )
        dist_left = abs(float(pos.x()) - float(rect.left()))
        dist_right = abs(float(pos.x()) - float(rect.right()))
        if dist_left <= threshold and dist_left <= dist_right:
            return idx, "start"
        if dist_right <= threshold:
            return idx, "end"
        return None

    def _apply_resize_at(self, x: float) -> None:
        idx = self._resize_idx
        edge = self._resize_edge
        if idx is None or edge is None or not (0 <= idx < len(self._features)):
            return
        bp = self._x_to_bp(x)
        start, end = self._resize_orig_start, self._resize_orig_end
        if edge == "start":
            start = min(bp, end)
        else:
            end = max(bp, start)
        feat = self._features[idx]
        if int(feat.start) == start and int(feat.end) == end:
            return
        self._features[idx] = _feature_with_bounds(feat, start, end)
        self.update()

    def _row_band_rect(self, row: int) -> QRectF:
        row_h = self._row_height()
        y = _ANNOT_TRACK_VPAD + row * row_h
        return QRectF(0.0, y, float(max(1, self.width())), row_h)

    def _row_index_at(self, y: float) -> int | None:
        row_h = self._row_height()
        if row_h <= 0 or self._row_count <= 0:
            return None
        rel = y - _ANNOT_TRACK_VPAD
        if rel < 0:
            return None
        row = int(rel // row_h)
        if row < 0 or row >= self._row_count:
            return None
        return row

    def _marquee_rect(self) -> QRectF:
        if self._marquee_origin is None or self._marquee_current is None:
            return QRectF()
        rect = QRectF(self._marquee_origin, self._marquee_current).normalized()
        if rect.width() < 1.0:
            rect.setWidth(1.0)
        if rect.height() < 1.0:
            rect.setHeight(1.0)
        # Mostly horizontal drag: limit selection to the row where the drag started.
        if rect.height() <= _MARQUEE_DRAG_THRESHOLD + 2 and self._press_pos is not None:
            row = self._row_index_at(self._press_pos.y())
            if row is not None:
                band = self._row_band_rect(row)
                return QRectF(rect.left(), band.top(), rect.width(), band.height())
        return rect

    def _indices_in_marquee(self) -> set[int]:
        rect = self._marquee_rect()
        if rect.isNull():
            return set()
        return {
            idx
            for idx in range(len(self._features))
            if self._feature_row_rect(idx).intersects(rect)
        }

    def _apply_marquee_preview(self, *, additive: bool) -> None:
        picked = self._indices_in_marquee()
        if additive:
            self._selected_indices = set(self._marquee_additive_base) | picked
        else:
            self._selected_indices = picked

    def _emit_selection(self) -> None:
        self.selectionChanged.emit(frozenset(self._selected_indices))

    def selected_indices(self) -> set[int]:
        return set(self._selected_indices)

    def set_view(
        self,
        *,
        features: list[SequenceFeature],
        feature_rows: list[int],
        scroll_px: float,
        ppb: float,
        viewport_w: int,
        selected_indices: set[int] | frozenset[int] | None = None,
        feature_alphas: list[int] | None = None,
        seq_length: int | None = None,
    ) -> None:
        if not self._resize_active:
            self._features = list(features)
            self._feature_rows = list(feature_rows)
            if feature_alphas is not None and len(feature_alphas) == len(features):
                self._feature_alphas = list(feature_alphas)
            else:
                self._feature_alphas = [255] * len(features)
        self._scroll_px = max(0.0, float(scroll_px))
        self._ppb = max(_FIT_MIN_PPB, float(ppb))
        self._viewport_w = max(1, int(viewport_w))
        if seq_length is not None:
            self._seq_length = max(1, int(seq_length))
        if (
            not self._marquee_active
            and not self._move_active
            and not self._resize_active
            and self._press_pos is None
        ):
            self._selected_indices = set(selected_indices or [])
        valid = set(range(len(self._features)))
        self._selected_indices &= valid
        if self._hover_idx is not None and self._hover_idx not in valid:
            self._hover_idx = None
        self._row_count = (
            max((r for r in self._feature_rows), default=-1) + 1
            if self._features
            else 0
        )
        track_w = max(1, int(viewport_w))
        bar_h = int(self._feature_bar_height)
        if self._row_count > 0:
            h = int(self._row_count * self._row_height()) + 2 * _ANNOT_TRACK_VPAD
            h = max(h, bar_h + 2 * _ANNOT_TRACK_VPAD)
            self.setFixedSize(track_w, h)
        else:
            # Keep a visible empty track so the annotation panel is not blank.
            self.setFixedSize(
                track_w, max(_ANNOT_MIN_VIEWPORT_HEIGHT, bar_h + 2 * _ANNOT_TRACK_VPAD)
            )
        self.show()
        self.update()

    def has_content(self) -> bool:
        return self._row_count > 0

    def content_height(self) -> int:
        """Height of feature rows only (not viewport-stretched empty fill)."""
        bar_h = int(self._feature_bar_height)
        if self._row_count > 0:
            h = int(self._row_count * self._row_height()) + 2 * _ANNOT_TRACK_VPAD
            return max(h, bar_h + 2 * _ANNOT_TRACK_VPAD)
        return max(_ANNOT_MIN_VIEWPORT_HEIGHT, bar_h + 2 * _ANNOT_TRACK_VPAD)

    # -- cytogenetic image markers ---------------------------------------

    def set_image_marker_lookup(self, lookup) -> None:
        """``lookup(features) -> set[int]`` of feature indices with images."""
        self._image_marker_lookup = lookup
        self.update()

    def _refresh_marked_indices(self) -> None:
        lookup = self._image_marker_lookup
        if lookup is None or not self._features:
            self._marked_indices = set()
            return
        try:
            self._marked_indices = set(lookup(self._features))
        except Exception:  # noqa: BLE001 - never break painting over images
            logger.exception("Could not look up cytogenetic image markers")
            self._marked_indices = set()

    def has_images(self, idx: int) -> bool:
        return idx in self._marked_indices

    def _image_marker_rect(self, idx: int) -> QRectF:
        rect = self._feature_rect(idx)
        d = float(_IMAGE_MARKER_SIZE)
        x = rect.right() - d - 2.0 if rect.width() >= d + 6.0 else rect.right() + 2.0
        return QRectF(x, rect.center().y() - d / 2.0, d, d)

    def _image_marker_hit(self, pos: QPointF) -> int | None:
        for idx in sorted(self._marked_indices, reverse=True):
            if 0 <= idx < len(self._features) and self._image_marker_rect(idx).adjusted(-2, -2, 2, 2).contains(pos):
                return idx
        return None

    def _draw_image_marker(self, painter: QPainter, idx: int) -> None:
        rect = self._image_marker_rect(idx)
        painter.save()
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        painter.setPen(QPen(QColor("#1f1f1f"), 1.2))
        painter.setBrush(QColor("#ffffff"))
        painter.drawEllipse(rect)
        inner = rect.adjusted(3.0, 3.0, -3.0, -3.0)
        painter.setPen(Qt.PenStyle.NoPen)
        painter.setBrush(QColor("#e0282e"))
        painter.drawEllipse(inner)
        painter.restore()

    def dragEnterEvent(self, event) -> None:  # noqa: N802
        from repeatviewer.widgets.cyto_images_ui import image_paths_from_mime

        if self._features and image_paths_from_mime(event.mimeData()):
            event.acceptProposedAction()
        else:
            event.ignore()

    def dragMoveEvent(self, event) -> None:  # noqa: N802
        from repeatviewer.widgets.cyto_images_ui import image_paths_from_mime

        if not image_paths_from_mime(event.mimeData()):
            event.ignore()
            return
        hit = self._hit_test_feature(event.position())
        if hit != self._drop_hover_idx:
            self._drop_hover_idx = hit
            self.update()
        if hit is None:
            event.ignore()
        else:
            event.acceptProposedAction()

    def dragLeaveEvent(self, event) -> None:  # noqa: N802
        if self._drop_hover_idx is not None:
            self._drop_hover_idx = None
            self.update()
        super().dragLeaveEvent(event)

    def dropEvent(self, event) -> None:  # noqa: N802
        from repeatviewer.widgets.cyto_images_ui import image_paths_from_mime

        paths = image_paths_from_mime(event.mimeData())
        hit = self._hit_test_feature(event.position())
        self._drop_hover_idx = None
        self.update()
        if not paths or hit is None:
            event.ignore()
            return
        event.acceptProposedAction()
        self._selected_indices = {hit}
        self._emit_selection()
        QTimer.singleShot(0, lambda i=hit, p=list(paths): self.imagesDropped.emit(i, p))

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.fillRect(self.rect(), _ANNOT_PAPER)
        ppb = self._ppb
        scroll = self._scroll_px
        vw = self._viewport_w
        # Annotation selection: light-grey unselected + dashed L/R guides.
        # Selected band keeps the normal paper (no blue fill here).
        if self._range_selection is not None and ppb > 0:
            a, b = self._range_selection
            if a != b:
                x1 = _bp_view_x(a, scroll, ppb)
                x2 = _bp_view_x(b + 1, scroll, ppb)
                h = float(self.height())
                clip_r = QRectF(event.rect())
                left = QRectF(0.0, 0.0, max(0.0, x1), h)
                right = QRectF(x2, 0.0, max(0.0, float(vw) - x2), h)
                for dim in (left, right):
                    if dim.width() > 0.5 and dim.intersects(clip_r):
                        painter.fillRect(dim, _SELECTION_DIM)
                band = QRectF(x1, 0.0, max(2.0, x2 - x1), h)
                if band.intersects(clip_r):
                    _draw_selection_side_guides(painter, band)
        if not self._features:
            self._marked_indices = set()
            painter.end()
            return
        self._refresh_marked_indices()
        row_h = self._row_height()
        painter.setPen(Qt.PenStyle.NoPen)

        from repeatviewer.backends.annotation_colors import (
            pident_color_hex,
            pident_scale_bounds,
        )
        from repeatviewer.backends.ncbi_blast_online import ONLINE_BLAST_SOURCE

        online_pidents = [
            float(feat.pident)
            for feat in self._features
            if str(feat.source or "") == ONLINE_BLAST_SOURCE and float(feat.pident or 0) > 0
        ]
        zmin, zmax = pident_scale_bounds(online_pidents)

        # Precompute feature rects for arc connectors.
        rects: list[QRectF | None] = [None] * len(self._features)
        for idx, feat in enumerate(self._features):
            x1 = _bp_view_x(feat.start, scroll, ppb)
            x2 = _bp_view_x(feat.end + 1, scroll, ppb)
            if x2 < _LEFT_MARGIN or x1 > vw:
                continue
            row = self._feature_rows[idx] if idx < len(self._feature_rows) else 0
            y = _ANNOT_TRACK_VPAD + row * row_h + 0.5
            rects[idx] = QRectF(
                x1, y, max(2.0, x2 - x1), max(2.0, float(self._feature_bar_height) - 1.0)
            )

        online_bar_h = 10.0
        for idx, feat in enumerate(self._features):
            if str(feat.source or "") != ONLINE_BLAST_SOURCE:
                continue
            row = self._feature_rows[idx] if idx < len(self._feature_rows) else 0
            x1 = _bp_view_x(feat.start, scroll, ppb)
            x2 = _bp_view_x(feat.end + 1, scroll, ppb)
            if x2 < _LEFT_MARGIN or x1 > vw:
                continue
            y_row = _ANNOT_TRACK_VPAD + row * row_h
            y_off = (float(self._feature_bar_height) - online_bar_h) / 2.0
            rects[idx] = QRectF(
                x1, y_row + y_off, max(2.0, x2 - x1), max(2.0, online_bar_h)
            )

        # NCBI graphic summary: thin horizontal connectors between HSPs on one row.
        groups: dict[str, list[int]] = {}
        for idx, feat in enumerate(self._features):
            if str(feat.source or "") != ONLINE_BLAST_SOURCE:
                continue
            key = str(feat.subject_id or feat.parent or "").strip()
            if not key:
                continue
            groups.setdefault(key, []).append(idx)
        painter.save()
        conn_pen = QPen(QColor("#444444"))
        conn_pen.setWidthF(1.0)
        painter.setPen(conn_pen)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        for _key, members in groups.items():
            if len(members) < 2:
                continue
            ordered = sorted(members, key=lambda i: self._features[i].start)
            row0 = self._feature_rows[ordered[0]] if ordered else 0
            y_line = _ANNOT_TRACK_VPAD + row0 * row_h + float(self._feature_bar_height) / 2.0
            for left_i, right_i in zip(ordered, ordered[1:]):
                if self._feature_rows[left_i] != self._feature_rows[right_i]:
                    continue
                left = rects[left_i]
                right = rects[right_i]
                if left is None or right is None:
                    continue
                x1 = left.right()
                x2 = right.left()
                if x2 > x1 + 0.5:
                    painter.drawLine(QPointF(x1, y_line), QPointF(x2, y_line))
        painter.restore()
        painter.setPen(Qt.PenStyle.NoPen)

        for idx, feat in enumerate(self._features):
            rect = rects[idx]
            if rect is None:
                continue
            is_online = str(feat.source or "") == ONLINE_BLAST_SOURCE
            if is_online and float(feat.pident or 0) > 0:
                color = QColor(pident_color_hex(float(feat.pident), zmin, zmax))
            else:
                color = QColor(feat.color or feature_type_color(feat.feature_type))
            if self._feature_alphas and idx < len(self._feature_alphas):
                color.setAlpha(int(self._feature_alphas[idx]))
            if is_online:
                label = ""
            else:
                label = self._feature_track_label(feat, rect.width(), 28.0)
            selected = idx in self._selected_indices
            if is_online:
                # NCBI graphic summary style: plain bars (arrows make overlaps look messy).
                painter.setPen(Qt.PenStyle.NoPen)
                painter.setBrush(color)
                painter.drawRect(rect)
            else:
                _draw_directed_feature(
                    painter,
                    rect,
                    color,
                    forward=feat.is_forward(),
                    label=label,
                    label_point_size=self._label_point_size,
                )
            if selected:
                _draw_feature_selection_chrome(painter, rect)
            if idx == self._drop_hover_idx:
                painter.save()
                painter.setPen(QPen(QColor("#2c6fad"), 2, Qt.PenStyle.DashLine))
                painter.setBrush(Qt.BrushStyle.NoBrush)
                painter.drawRect(rect.adjusted(-1, -1, 1, 1))
                painter.restore()
                painter.setPen(Qt.PenStyle.NoPen)

        for idx in sorted(self._marked_indices):
            if 0 <= idx < len(self._features):
                marker = self._image_marker_rect(idx)
                if marker.right() >= _LEFT_MARGIN and marker.left() <= vw:
                    self._draw_image_marker(painter, idx)

        if self._move_active and self._move_ghost_y is not None and self._move_indices:
            for idx in sorted(self._move_indices):
                if not (0 <= idx < len(self._features)):
                    continue
                feat = self._features[idx]
                x1 = _bp_view_x(feat.start, scroll, ppb)
                x2 = _bp_view_x(feat.end + 1, scroll, ppb)
                rect = QRectF(
                    x1,
                    self._move_ghost_y,
                    max(2.0, x2 - x1),
                    max(2.0, float(self._feature_bar_height) - 1.0),
                )
                if (
                    str(feat.source or "") == ONLINE_BLAST_SOURCE
                    and float(feat.pident or 0) > 0
                ):
                    color = QColor(pident_color_hex(float(feat.pident), zmin, zmax))
                else:
                    color = QColor(feat.color or feature_type_color(feat.feature_type))
                color.setAlpha(110)
                if str(feat.source or "") == ONLINE_BLAST_SOURCE:
                    painter.setPen(Qt.PenStyle.NoPen)
                    painter.setBrush(color)
                    painter.drawRect(rect)
                else:
                    _draw_directed_feature(
                        painter,
                        rect,
                        color,
                        forward=feat.is_forward(),
                        label="",
                    )
                painter.setPen(QPen(_FEATURE_SELECTION_BORDER, 1, Qt.PenStyle.DashLine))
                painter.setBrush(Qt.BrushStyle.NoBrush)
                painter.drawRect(rect)

        hover_idx = self._hover_idx
        if (
            hover_idx is not None
            and 0 <= hover_idx < len(self._features)
            and not self._marquee_active
            and not self._move_active
            and not self._resize_active
        ):
            from repeatviewer.backends.sequence_features import (
                feature_has_hover_meta,
                feature_hover_detail_text,
            )

            feat = self._features[hover_idx]
            rect = self._feature_rect(hover_idx)
            tip_text = feature_hover_detail_text(
                feat, primary=self._feature_full_label(feat)
            )
            show_tip = bool(tip_text) and (
                feature_has_hover_meta(feat) or self._label_overflows(feat, rect)
            )
            if show_tip:
                tip_font = QFont(painter.font())
                tip_font.setPointSize(8)
                painter.setFont(tip_font)
                fm = QFontMetrics(tip_font)
                lines = tip_text.split("\n")
                tip_w = min(
                    max(max(fm.horizontalAdvance(line) for line in lines) + 10, 40),
                    max(80, self.width() - 8),
                )
                tip_h = fm.height() * len(lines) + 6
                tip_x = min(
                    max(4.0, rect.left()),
                    max(4.0, float(self.width()) - tip_w - 4.0),
                )
                tip_y = rect.bottom() + 3.0
                if tip_y + tip_h > self.height() - 2:
                    tip_y = max(2.0, rect.top() - tip_h - 3.0)
                tip_rect = QRectF(tip_x, tip_y, tip_w, tip_h)
                painter.setPen(QPen(QColor("#555555"), 1))
                painter.setBrush(QColor(255, 255, 255, 235))
                painter.drawRoundedRect(tip_rect, 3, 3)
                painter.setPen(QColor("#222222"))
                for index, line in enumerate(lines):
                    line_rect = QRectF(
                        tip_rect.left() + 5.0,
                        tip_rect.top() + 3.0 + index * fm.height(),
                        tip_rect.width() - 10.0,
                        float(fm.height()),
                    )
                    painter.drawText(
                        line_rect,
                        Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                        line,
                    )

        if self._marquee_active:
            marquee = self._marquee_rect()
            painter.setPen(QPen(QColor("#2c6fad"), 1, Qt.PenStyle.DashLine))
            painter.setBrush(QColor(44, 111, 173, 50))
            painter.drawRect(marquee)
        painter.end()

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() != Qt.MouseButton.LeftButton or not self._features:
            super().mousePressEvent(event)
            return
        pos = event.position()
        marker_hit = self._image_marker_hit(pos)
        if marker_hit is not None:
            self._selected_indices = {marker_hit}
            self._emit_selection()
            self.update()
            self.openImagesRequested.emit(marker_hit)
            event.accept()
            return
        self._press_pos = pos
        edge = self._edge_hit(pos)
        if edge is not None:
            self._press_hit, self._press_edge = edge
        else:
            self._press_hit = self._hit_test_feature(pos)
            self._press_edge = None
        self._marquee_active = False
        self._marquee_origin = None
        self._marquee_current = None
        self._move_active = False
        self._move_indices = set()
        self._move_preview_track = None
        self._move_ghost_y = None
        self._resize_active = False
        self._resize_idx = None
        self._resize_edge = None
        self.setFocus(Qt.FocusReason.MouseFocusReason)
        event.accept()

    def mouseMoveEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        pos = event.position()
        if self._press_pos is None:
            edge = self._edge_hit(pos)
            hit = edge[0] if edge is not None else self._hit_test_feature(pos)
            if hit != self._hover_idx:
                self._hover_idx = hit
                self.update()
            if edge is not None:
                self.setCursor(Qt.CursorShape.SizeHorCursor)
            else:
                self.unsetCursor()
            super().mouseMoveEvent(event)
            return

        if self._resize_active:
            self._apply_resize_at(pos.x())
            event.accept()
            return

        if self._move_active:
            self._move_ghost_y = pos.y() - float(self._feature_bar_height) * 0.5
            self._move_preview_track = self._target_track_at_y(pos.y())
            self.update()
            event.accept()
            return

        if self._marquee_active:
            self._marquee_current = pos
            self._apply_marquee_preview(additive=self._additive_modifier(event.modifiers()))
            self.update()
            return

        delta = pos - self._press_pos
        if math.hypot(delta.x(), delta.y()) < _MARQUEE_DRAG_THRESHOLD:
            super().mouseMoveEvent(event)
            return

        # Edge grab → resize start/end (takes priority over track move / marquee).
        if self._press_edge is not None and self._press_hit is not None:
            idx = self._press_hit
            feat = self._features[idx]
            self._resize_active = True
            self._resize_idx = idx
            self._resize_edge = self._press_edge
            self._resize_orig_start = int(feat.start)
            self._resize_orig_end = int(feat.end)
            if idx not in self._selected_indices:
                self._selected_indices = {idx}
                self._emit_selection()
            self._hover_idx = None
            self.setCursor(Qt.CursorShape.SizeHorCursor)
            self.grabMouse()
            self._mouse_grabbed = True
            self._apply_resize_at(pos.x())
            event.accept()
            return

        # Drag from an already-selected feature (mostly vertical) → move track.
        # Otherwise → marquee multi-select.
        if (
            self._can_start_move(self._press_hit)
            and not self._additive_modifier(event.modifiers())
            and abs(delta.y()) >= abs(delta.x())
        ):
            self._move_active = True
            self._move_indices = set(self._selected_indices)
            self._move_ghost_y = pos.y() - float(self._feature_bar_height) * 0.5
            self._move_preview_track = self._target_track_at_y(pos.y())
            self._hover_idx = None
            self.setCursor(Qt.CursorShape.SizeVerCursor)
            self.grabMouse()
            self._mouse_grabbed = True
            self.update()
            event.accept()
            return

        self._marquee_active = True
        self._marquee_origin = self._press_pos
        self._marquee_current = pos
        if self._additive_modifier(event.modifiers()):
            self._marquee_additive_base = set(self._selected_indices)
        else:
            self._marquee_additive_base = set()
            self._selected_indices.clear()
        self._apply_marquee_preview(additive=self._additive_modifier(event.modifiers()))
        self.setCursor(Qt.CursorShape.CrossCursor)
        self.grabMouse()
        self._mouse_grabbed = True
        self.update()

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() != Qt.MouseButton.LeftButton:
            super().mouseReleaseEvent(event)
            return
        if (
            self._press_pos is None
            and not self._marquee_active
            and not self._move_active
            and not self._resize_active
        ):
            self._release_drag_mouse()
            super().mouseReleaseEvent(event)
            return
        pos = event.position()
        additive = self._additive_modifier(event.modifiers())
        try:
            if self._resize_active and self._resize_idx is not None:
                idx = self._resize_idx
                feat = self._features[idx] if 0 <= idx < len(self._features) else None
                changed = (
                    feat is not None
                    and (
                        int(feat.start) != self._resize_orig_start
                        or int(feat.end) != self._resize_orig_end
                    )
                )
                if changed:
                    self.featuresResized.emit(list(self._features))
                self._emit_selection()
            elif self._move_active and self._move_indices:
                updated = self._features_after_track_move(pos.y())
                if updated is not None:
                    self.tracksMoved.emit(updated)
                self._emit_selection()
            elif self._marquee_active:
                picked = self._indices_in_marquee()
                if additive:
                    self._selected_indices = set(self._marquee_additive_base) | picked
                else:
                    self._selected_indices = picked
                self._emit_selection()
            elif self._press_pos is not None:
                hit = self._press_hit
                if hit is None:
                    hit = self._hit_test_feature(pos)
                if hit is not None:
                    if additive:
                        if hit in self._selected_indices:
                            self._selected_indices.remove(hit)
                        else:
                            self._selected_indices.add(hit)
                    else:
                        self._selected_indices = {hit}
                elif not additive:
                    self._selected_indices.clear()
                    # Click on grey (outside range selection) clears bp selection.
                    if self._range_selection is not None and self._ppb > 0:
                        a, b = self._range_selection
                        if a != b:
                            bp = int((pos.x() + self._scroll_px - _LEFT_MARGIN) / self._ppb) + 1
                            if bp < a or bp > b:
                                self.clearRangeSelectionRequested.emit()
                self._emit_selection()
        finally:
            self._press_pos = None
            self._press_hit = None
            self._press_edge = None
            self._marquee_active = False
            self._marquee_origin = None
            self._marquee_current = None
            self._marquee_additive_base = set()
            self._move_active = False
            self._move_indices = set()
            self._move_preview_track = None
            self._move_ghost_y = None
            self._resize_active = False
            self._resize_idx = None
            self._resize_edge = None
            self._release_drag_mouse()
            self.setFocus(Qt.FocusReason.MouseFocusReason)
            self.update()
        event.accept()

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802
        if _is_delete_key(event.key()) and self._selected_indices:
            self.deleteRequested.emit()
            event.accept()
            return
        super().keyPressEvent(event)

    def contextMenuEvent(self, event: QContextMenuEvent) -> None:  # noqa: N802
        if not self._features:
            return
        hit = self._hit_test_feature(QPointF(event.pos()))
        if hit is None:
            return
        if hit not in self._selected_indices:
            self._selected_indices = {hit}
            self._emit_selection()
            self.update()
        menu = QMenu(self)
        menu.setObjectName("AnnotationFeatureMenu")
        delete_action = menu.addAction("Delete")
        delete_action.triggered.connect(self.deleteRequested.emit)
        copy_action = menu.addAction("Copy sequence")
        copy_action.triggered.connect(self.copySequenceRequested.emit)
        extract_action = menu.addAction("Extract Sequence")
        extract_action.triggered.connect(self.extractSequenceRequested.emit)
        if len(self._selected_indices) > 1:
            from repeatviewer.backends.ncbi_blast_online import ONLINE_BLAST_SOURCE

            has_online = any(
                0 <= i < len(self._features)
                and str(self._features[i].source or "") == ONLINE_BLAST_SOURCE
                for i in self._selected_indices
            )
            if not has_online:
                color_action = menu.addAction("Change color")
                color_action.triggered.connect(self.changeColorRequested.emit)
        else:
            modify_action = menu.addAction("Modify")
            modify_action.triggered.connect(
                lambda _checked=False, idx=hit: self._request_modify(idx)
            )
        if len(self._selected_indices) == 1:
            menu.addSeparator()
            attach_action = menu.addAction("Attach cytogenetic image…")
            attach_action.triggered.connect(self.attachImageRequested.emit)
            if hit in self._marked_indices:
                show_action = menu.addAction("Show cytogenetic images")
                show_action.triggered.connect(
                    lambda _checked=False, idx=hit: self.openImagesRequested.emit(idx)
                )
        menu.exec(event.globalPos())
        event.accept()

    def _request_modify(self, index: int) -> None:
        self._selected_indices = {int(index)}
        self._emit_selection()
        self.modifyRequested.emit()
        self.update()


class SequenceMapCanvas(QWidget):
    """Paints sequence + feature tracks; supports wrap and stacked feature rows."""

    selectionChanged = Signal(int, int)  # 1-based inclusive
    selectionCleared = Signal()
    featureClicked = Signal(int)  # index in features list
    addFeatureRequested = Signal(str, str)  # feature_type, label
    copySequenceRequested = Signal(int, int)  # 1-based inclusive range
    extractSequenceRequested = Signal(int, int)  # 1-based inclusive range
    featureActionRequested = Signal(str)  # delete | copy | extract | modify
    # Parent must re-size canvas_row from content height (not a clipped layout height).
    contentHeightChanged = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._sequence = ""
        self._features: list[SequenceFeature] = []
        self._feature_rows: list[int] = []
        self._feature_row_count = 1
        self._pixels_per_bp = 1.0
        self._viewport_width = 0
        self._scroll_px = 0.0
        self._color_bases = True
        self._is_protein = False
        self._font_point_size = _DEFAULT_NUC_FONT_PT
        self._selected_indices: set[int] = set()
        self._select_start_bp: int | None = None
        self._select_end_bp: int | None = None
        self._find_hits: list[tuple[int, int]] = []
        self._find_starts: list[int] = []
        self._find_max_len = 0
        self._drag_select = False
        self._mouse_grabbed = False
        self._pending_selection_emit: tuple[int, int] | None = None
        self._selection_emit_timer = QTimer(self)
        self._selection_emit_timer.setSingleShot(True)
        self._selection_emit_timer.setInterval(33)  # ~30 Hz spin/status updates
        self._selection_emit_timer.timeout.connect(self._flush_pending_selection_emit)
        self._seq_layer_cache: QPixmap | None = None
        # Set by export_annotation_view_dialog.vector_export(): paint live so
        # PDF/SVG get vector letters instead of the screen pixmap cache.
        self._vector_export = False
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._show_similarity_database_name = True
        self._split_tracks_by_source = True
        self._show_complement = False
        self._translation_frames: frozenset[int] = frozenset()
        self._translation_table = 1
        self._translation_three_letter = False
        self._letter_min_px = 4
        self._hover_track_tip = ""

    def sizeHint(self) -> QSize:  # noqa: N802
        return QSize(max(120, self.width()), max(1, int(self._content_height())))

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return QSize(120, max(1, int(self._content_height())))

    def set_show_complement(self, enabled: bool) -> None:
        enabled = bool(enabled)
        changed = enabled != self._show_complement
        self._show_complement = enabled
        # Always refresh geometry/cache: parent row may still be clipping a prior grow.
        self._invalidate_seq_layer_cache()
        self._update_geometry()
        self.update()
        if not changed:
            self.contentHeightChanged.emit()

    def set_translation_frames(self, frames: object) -> None:
        cleaned = frozenset(int(f) for f in (frames or []) if int(f) in {1, 2, 3, -1, -2, -3})
        changed = cleaned != self._translation_frames
        self._translation_frames = cleaned
        self._invalidate_seq_layer_cache()
        self._update_geometry()
        self.update()
        if not changed:
            self.contentHeightChanged.emit()

    def set_translation_table(self, table_id: int) -> None:
        from repeatviewer.backends.genetic_codes import normalize_genetic_code_id

        table_id = normalize_genetic_code_id(table_id)
        if table_id == self._translation_table:
            return
        self._translation_table = table_id
        self._invalidate_seq_layer_cache()
        self.update()

    def set_translation_three_letter(self, enabled: bool) -> None:
        enabled = bool(enabled)
        if enabled == self._translation_three_letter:
            return
        self._translation_three_letter = enabled
        self._invalidate_seq_layer_cache()
        self.update()

    def set_letter_min_px(self, pixels: int) -> None:
        pixels = max(1, min(24, int(pixels)))
        if pixels == self._letter_min_px:
            return
        self._letter_min_px = pixels
        self._invalidate_seq_layer_cache()
        self._update_geometry()
        self.update()

    def set_split_tracks_by_source(self, enabled: bool) -> None:
        enabled = bool(enabled)
        if enabled == self._split_tracks_by_source:
            return
        self._split_tracks_by_source = enabled
        self._feature_rows, self._feature_row_count = _pack_feature_rows(
            self._features,
            split_by_source=self._split_tracks_by_source,
        )
        self._invalidate_seq_layer_cache()
        self._update_geometry()
        self.update()

    def _invalidate_seq_layer_cache(self) -> None:
        self._seq_layer_cache = None

    def _seq_layer_cache_dpr(self) -> float:
        return max(1.0, float(self.devicePixelRatioF()))

    def _ensure_seq_layer_cache(self) -> None:
        """Lazy sequence (+ wrap features) snapshot without selection chrome.

        Rebuilt on zoom/scroll/sequence change — not on every mouse press — so
        selection drag only blits the cache and redraws the highlight.

        The pixmap is allocated at device-pixel resolution (Retina / HiDPI) so
        letters stay sharp when zoomed in; a 1× cache looks soft on 2× displays.
        """
        dpr = self._seq_layer_cache_dpr()
        logical_w = self.width()
        # Always paint the intended strip height; layout must match via _update_geometry.
        logical_h = max(1, int(self._content_height()))
        if self.height() != logical_h:
            self.setMinimumHeight(logical_h)
            self.setMaximumHeight(logical_h)
            self.setFixedHeight(logical_h)
        if (
            self._seq_layer_cache is not None
            and not self._seq_layer_cache.isNull()
            and abs(float(self._seq_layer_cache.devicePixelRatio()) - dpr) < 1e-6
            and self._seq_layer_cache.width() == int(math.ceil(logical_w * dpr))
            and self._seq_layer_cache.height() == int(math.ceil(logical_h * dpr))
        ):
            return
        if logical_w < 2 or logical_h < 2 or not self._sequence:
            self._seq_layer_cache = None
            return
        cache = QPixmap(
            int(math.ceil(logical_w * dpr)),
            int(math.ceil(logical_h * dpr)),
        )
        cache.setDevicePixelRatio(dpr)
        cache.fill(_MAP_PAPER)
        painter = QPainter(cache)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        # Smooth glyphs; fills stay aliased so base cells stay crisp.
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
        self._draw_sequence(painter, QRectF(0, 0, logical_w, logical_h))
        if self.uses_wrap():
            self._draw_features(painter)
            track_bottom = self._feature_track_top() + self._feature_block_height()
            painter.setPen(QColor("#bbbbbb"))
            painter.drawLine(
                _LEFT_MARGIN,
                track_bottom + 2,
                self._content_width() - _LEFT_MARGIN,
                track_bottom + 2,
            )
        painter.end()
        self._seq_layer_cache = cache

    def set_show_similarity_database_name(self, show: bool) -> None:
        show = bool(show)
        if show == self._show_similarity_database_name:
            return
        self._show_similarity_database_name = show
        self._invalidate_seq_layer_cache()
        self.update()

    def _feature_track_label(self, feat: SequenceFeature, rect_width: float, min_width: float) -> str:
        label_pt = getattr(self, "_label_point_size", None)
        if (label_pt is not None and int(label_pt) <= 0) or rect_width <= min_width:
            return ""
        base = feat.name or feat.feature_type
        return similarity_feature_track_label(
            base,
            feat.feature_type,
            show_database_name=self._show_similarity_database_name,
        )

    def _release_drag_mouse(self) -> None:
        if self._mouse_grabbed:
            self.releaseMouse()
            self._mouse_grabbed = False

    def hideEvent(self, event) -> None:  # noqa: N802
        self._drag_select = False
        self._release_drag_mouse()
        super().hideEvent(event)

    def resizeEvent(self, event) -> None:  # noqa: N802
        self._invalidate_seq_layer_cache()
        super().resizeEvent(event)

    def set_font_point_size(self, point_size: int) -> None:
        self._font_point_size = max(_MIN_NUC_FONT_PT, min(_MAX_NUC_FONT_PT, int(point_size)))
        self._invalidate_seq_layer_cache()
        self._update_geometry()
        self.update()

    def font_point_size(self) -> int:
        return self._font_point_size

    def _sequence_top(self) -> int:
        # Flush to the top so the selection box lines up with the sequence track.
        return 0

    def set_sequence(
        self,
        sequence: str,
        features: list[SequenceFeature],
        *,
        is_protein: bool = False,
    ) -> None:
        self._sequence = sequence.upper()
        self._features = list(features)
        self._is_protein = bool(is_protein)
        self._feature_rows, self._feature_row_count = _pack_feature_rows(
            self._features,
            split_by_source=self._split_tracks_by_source,
        )
        self._invalidate_seq_layer_cache()
        self._update_geometry()
        self.update()

    def set_pixels_per_bp(self, ppb: float) -> None:
        cap = _zoom_ppb_max_for(len(self._sequence))
        self._pixels_per_bp = max(_FIT_MIN_PPB, min(cap, float(ppb)))
        self._invalidate_seq_layer_cache()
        self._update_geometry()
        self.update()

    def zoom_ppb_max(self) -> float:
        return _zoom_ppb_max_for(len(self._sequence))

    def set_viewport_width(self, width: int) -> None:
        """Update paint viewport width (must not drive parent layout width)."""
        width = max(120, int(width))
        if width == self._viewport_width:
            return
        self._viewport_width = width
        self._invalidate_seq_layer_cache()
        self.update()

    def sync_viewport_from_widget(self) -> None:
        """Read the layout-assigned width; never expand the parent."""
        width = max(120, self.width())
        if width != self._viewport_width:
            self._viewport_width = width
            self._invalidate_seq_layer_cache()
            self.update()

    def set_scroll_px(self, scroll_px: float) -> None:
        self._scroll_px = max(0.0, float(scroll_px))
        self._invalidate_seq_layer_cache()
        self.update()
    def scroll_px(self) -> float:
        return self._scroll_px

    def content_pixel_width(self) -> int:
        """Full horizontal extent of the sequence at the current zoom (for scroll range)."""
        if not self._sequence:
            return max(120, self._viewport_width)
        if self.uses_wrap():
            return max(120, self._viewport_width)
        return _LEFT_MARGIN + self._sequence_pixel_width() + _RIGHT_MARGIN

    def bp_range_center_content_x(self, start: int, end: int) -> float:
        """Horizontal center of a bp range in content coordinates (scroll=0 origin)."""
        _, _, cell_w, _ = self._nucleotide_cell_geometry()
        a, b = sorted((int(start), int(end)))
        return _LEFT_MARGIN + (a + b - 1) * float(cell_w) / 2.0

    def pixels_per_bp(self) -> float:
        return self.effective_ppb()

    def uses_wrap(self) -> bool:
        return False

    def set_selection_range(self, start: int, end: int) -> None:
        """Select a 1-based inclusive bp interval."""
        if not self._sequence:
            return
        seq_len = len(self._sequence)
        a = max(1, min(seq_len, int(start)))
        b = max(1, min(seq_len, int(end)))
        if a > b:
            a, b = b, a
        self._select_start_bp = a
        self._select_end_bp = b
        self._selected_indices.clear()
        self._drag_select = False
        self.selectionChanged.emit(a, b)
        self.update()

    def set_color_bases(self, enabled: bool) -> None:
        self._color_bases = bool(enabled)
        self._invalidate_seq_layer_cache()
        self.update()

    def color_bases(self) -> bool:
        return self._color_bases

    def set_selected_features(self, indices: set[int] | frozenset[int] | None) -> None:
        self._selected_indices = set(indices or [])
        if self.uses_wrap():
            self._invalidate_seq_layer_cache()
        self.update()

    def set_selected_feature(self, index: int | None) -> None:
        self._selected_indices = {index} if index is not None else set()
        self.update()

    def selected_feature_indices(self) -> set[int]:
        return set(self._selected_indices)

    def selected_feature_index(self) -> int | None:
        if len(self._selected_indices) == 1:
            return next(iter(self._selected_indices))
        return None

    def clear_feature_selection(self) -> None:
        self._selected_indices.clear()
        self.update()

    def clear_selection(self) -> None:
        had_range = self.selection_range() is not None
        self._select_start_bp = None
        self._select_end_bp = None
        self._selected_indices.clear()
        self._drag_select = False
        if had_range:
            self.selectionCleared.emit()
        self.update()

    def set_find_hits(self, hits: list[tuple[int, int]]) -> None:
        """Edit > Find hits to highlight (1-based inclusive ranges); [] clears."""
        self._find_hits = sorted((min(int(a), int(b)), max(int(a), int(b))) for a, b in hits)
        self._find_starts = [a for a, _b in self._find_hits]
        self._find_max_len = max((b - a + 1 for a, b in self._find_hits), default=0)
        self.update()

    def find_hits(self) -> list[tuple[int, int]]:
        return list(self._find_hits)

    def selection_range(self) -> tuple[int, int] | None:
        if self._select_start_bp is None or self._select_end_bp is None:
            return None
        a, b = sorted((self._select_start_bp, self._select_end_bp))
        return a, b

    def _sequence_font(self) -> QFont:
        font = _monospace_font(point_size=self._font_point_size)
        return font

    def _sequence_font_for_track(self, cell_h: float) -> QFont:
        """Letter font that fits the track height (shrink only if the cell is short)."""
        font = self._sequence_font()
        max_h = max(1.0, float(cell_h) - 2.0)
        # Prefer tight bounding box — height() alone under-reports some macOS fonts
        # and left the bottom half of glyphs clipped.
        min_pt = int(_MIN_NUC_FONT_PT)
        for pt in range(max(min_pt, int(self._font_point_size)), min_pt - 1, -1):
            font.setPointSize(pt)
            fm = QFontMetrics(font)
            glyph_h = max(fm.height(), fm.ascent() + fm.descent())
            if glyph_h <= max_h:
                return font
        font.setPointSize(min_pt)
        return font

    def _font_letter_cell_size(self) -> tuple[int, int]:
        """Monospace cell size for nucleotide / amino-acid letters.

        Height always tracks the Display DNA/protein label size (no fixed floor),
        so shrinking the spin below ~8pt also shrinks the colored backgrounds.
        """
        fm = QFontMetrics(self._sequence_font())
        sample = "ACDEFGHIKLMNPQRSTVWYBXZ*" if self._is_protein else "ACGTUN"
        cell_w = max(1, max(fm.horizontalAdvance(base) for base in sample))
        glyph_h = max(fm.height(), fm.ascent() + fm.descent())
        cell_h = max(1, glyph_h + 4)
        return cell_w, cell_h

    def shows_letters(self) -> bool:
        """True when zoomed in enough to draw A/C/G/T (or AA) glyphs.

        Letters are painted only for the visible viewport, so this is safe for
        long sequences as well as short ones — no full-sequence letter pass.
        """
        if not self._sequence:
            return False
        font_w, _ = self._font_letter_cell_size()
        threshold = max(float(self._letter_min_px), float(font_w))
        return self._pixels_per_bp >= threshold

    def _extra_display_rows(self) -> list[str]:
        """Row kinds under the forward strand (complement / translation).

        Only at nucleotide (letter) zoom — zoomed-out bar mode stays a single strip.
        """
        if self._is_protein or not self.shows_letters():
            return []
        rows: list[str] = []
        if self._show_complement:
            rows.append("complement")
        for frame in sorted(self._translation_frames, key=lambda f: (f < 0, abs(f))):
            rows.append(f"frame:{frame}")
        return rows

    @staticmethod
    def _track_kind_label(kind: str, *, is_protein: bool = False) -> str:
        if kind == "forward":
            return "Protein sequence" if is_protein else "Forward strand"
        if kind == "complement":
            return "Reverse strand"
        if kind.startswith("frame:"):
            try:
                frame = int(kind.split(":", 1)[1])
            except ValueError:
                return "Translation"
            if frame > 0:
                return f"Forward frame {frame}"
            return f"Reverse frame {abs(frame)}"
        return ""

    def _display_row_kind_at_y(self, y: float) -> str | None:
        """Return ``forward`` / ``complement`` / ``frame:N`` under vertical ``y``, if any."""
        if not self._sequence:
            return None
        _origin_x, origin_y, _cell_w, cell_h = self._nucleotide_cell_geometry()
        if cell_h <= 0:
            return None
        y0 = float(origin_y)
        h = float(cell_h)
        if y0 <= y < y0 + h:
            return "forward"
        gap = 2.0
        cursor = y0 + h + gap
        for kind in self._extra_display_rows():
            if cursor <= y < cursor + h:
                return kind
            cursor += h + gap
        return None

    def _update_sequence_track_tooltip(self, pos: QPointF, global_pos) -> None:
        kind = self._display_row_kind_at_y(float(pos.y()))
        tip = (
            self._track_kind_label(kind, is_protein=self._is_protein)
            if kind
            else ""
        )
        if tip == self._hover_track_tip:
            return
        self._hover_track_tip = tip
        if tip:
            QToolTip.showText(global_pos, tip, self)
        else:
            QToolTip.hideText()

    def _extra_display_height(self) -> float:
        rows = self._extra_display_rows()
        if not rows:
            return 0.0
        _, _, _, cell_h = self._nucleotide_cell_geometry()
        # Matches _draw_extra_sequence_rows: 2px gap (includes hairline) + cell per row.
        gap = 2.0
        return len(rows) * (float(cell_h) + gap)

    def _letter_cell_size(self) -> tuple[int, int]:
        """Current nucleotide cell size when drawing letters."""
        cell_w = max(1, int(math.ceil(self._pixels_per_bp)))
        _, _, _, cell_h = self._nucleotide_cell_geometry()
        return cell_w, cell_h

    def _nucleotide_cell_geometry(self) -> tuple[int, int, float, int]:
        """Shared top-left origin and cell size for every nucleotide background.

        Bar mode keeps a compact strip. Letter mode sizes height from the Display
        font so glyphs are never clipped mid-letter.
        """
        origin_y = int(self._sequence_top())
        ppb = max(_FIT_MIN_PPB, float(self._pixels_per_bp))
        font_w, font_h = self._font_letter_cell_size()
        if self.shows_letters():
            cell_w = max(float(font_w), ppb)
            cell_h = int(font_h)
        else:
            cell_w = ppb
            cell_h = int(_SEQUENCE_ROW)
        return _LEFT_MARGIN, origin_y, cell_w, cell_h

    def effective_ppb(self) -> float:
        """Pixels per base on screen."""
        return max(_FIT_MIN_PPB, float(self._pixels_per_bp))

    def effective_cell_w(self) -> float:
        """Pixel width of one base on screen (matches drawing and hit-testing)."""
        _, _, cell_w, _ = self._nucleotide_cell_geometry()
        return float(cell_w)

    def bp_range_center_canvas_x(self, start: int, end: int) -> float:
        """Horizontal center of a bp interval in viewport coordinates."""
        return self.bp_range_center_content_x(start, end) - self._scroll_px

    def _cell_size(self) -> tuple[float, float, QFontMetrics]:
        font = self._sequence_font()
        fm = QFontMetrics(font)
        if self.shows_letters():
            cell_w, cell_h = self._letter_cell_size()
            return float(cell_w), float(cell_h), fm
        cell_w = max(1.0, float(fm.horizontalAdvance("A")))
        cell_h = max(1.0, float(fm.height()))
        return cell_w, cell_h, fm

    def _char_width(self) -> float:
        return self.effective_ppb()

    def _bases_per_line(self) -> int:
        if not self.uses_wrap():
            return max(1, len(self._sequence))
        inner = max(40, self._viewport_width - _LEFT_MARGIN * 2)
        return max(1, int(inner / self._char_width()))

    def _sequence_line_count(self) -> int:
        if not self._sequence:
            return 1
        if not self.uses_wrap():
            return 1
        return max(1, math.ceil(len(self._sequence) / self._bases_per_line()))

    def _features_on_line(self, line: int) -> list[tuple[int, int, int]]:
        bpl = self._bases_per_line()
        line_lo = line * bpl + 1
        line_hi = min(len(self._sequence), (line + 1) * bpl)
        segments: list[tuple[int, int, int]] = []
        for idx, feat in enumerate(self._features):
            if feat.end < line_lo or feat.start > line_hi:
                continue
            segments.append((idx, max(feat.start, line_lo), min(feat.end, line_hi)))
        return segments

    def _feature_rows_on_line(self, line: int) -> int:
        segments = self._features_on_line(line)
        if not segments:
            return 0
        _, rows = _pack_interval_rows(segments)
        return rows

    def _seq_row_height(self) -> float:
        # Same height in bar and letter mode (from Display font size only).
        _, _, _, cell_h = self._nucleotide_cell_geometry()
        return float(cell_h)

    def _line_block_height(self, line: int) -> float:
        if not self.uses_wrap():
            return 0.0
        height = _WRAP_LINE_RULER + self._seq_row_height()
        annot_rows = self._feature_rows_on_line(line)
        if annot_rows > 0:
            height += annot_rows * (_TRACK_HEIGHT + _TRACK_GAP)
        return height + _WRAP_LINE_GAP

    def _line_top(self, line: int) -> float:
        y = 4.0
        for ln in range(line):
            y += self._line_block_height(ln)
        return y

    def _sequence_block_height(self) -> int:
        if self.uses_wrap():
            total = 4.0
            for line in range(self._sequence_line_count()):
                total += self._line_block_height(line)
            return int(total)
        return int(self._seq_row_height() + self._extra_display_height())

    def _feature_block_height(self) -> int:
        if not self.uses_wrap():
            return 0
        if not self._features:
            return 0
        return self._feature_row_count * (_TRACK_HEIGHT + _TRACK_GAP) - _TRACK_GAP

    def _sequence_pixel_width(self) -> int:
        if not self._sequence:
            return 0
        ppb = max(_FIT_MIN_PPB, float(self._pixels_per_bp))
        return int(math.ceil(len(self._sequence) * ppb))

    def _content_width(self) -> int:
        return max(120, self.width())

    def _content_height(self) -> int:
        if self.uses_wrap():
            return int(self._sequence_block_height()) + 8
        return max(
            _SEQUENCE_ROW + 2,
            int(self._seq_row_height() + self._extra_display_height()) + 2,
        )

    def _update_geometry(self) -> None:
        height = max(1, int(self._content_height()))
        # Always pin min/max/fixed so a short parent layout cannot squeeze us
        # permanently and hide complement / translation rows.
        if (
            self.height() != height
            or self.minimumHeight() != height
            or self.maximumHeight() != height
        ):
            self.setMinimumHeight(height)
            self.setMaximumHeight(height)
            self.setFixedHeight(height)
            self.contentHeightChanged.emit()

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self.sync_viewport_from_widget()

    def changeEvent(self, event) -> None:  # noqa: N802
        super().changeEvent(event)
        # Moving between 1× / Retina screens must rebuild the letter cache.
        if event.type() == QEvent.Type.DevicePixelRatioChange:
            self._invalidate_seq_layer_cache()
            self.update()

    def _index_to_xy(self, index: int) -> tuple[float, float]:
        """0-based index -> top-left (x, y) of base cell."""
        if self.uses_wrap():
            bpl = self._bases_per_line()
            line = index // bpl
            col = index % bpl
            _, origin_y, cell_w, _ = self._nucleotide_cell_geometry()
            x = _LEFT_MARGIN + col * cell_w
            y = int(self._line_top(line) + _WRAP_LINE_RULER)
            return float(x), float(y)
        _, origin_y, cell_w, _ = self._nucleotide_cell_geometry()
        x = _LEFT_MARGIN + index * cell_w - self._scroll_px
        return float(x), float(origin_y)

    def _bp_at_pos(self, x: float, y: float) -> int:
        if not self._sequence:
            return 1
        seq_len = len(self._sequence)
        origin_x, _, cell_w, _ = self._nucleotide_cell_geometry()
        if self.uses_wrap():
            for line in range(self._sequence_line_count()):
                top = self._line_top(line)
                bottom = top + self._line_block_height(line) - _WRAP_LINE_GAP
                if y < top or y >= bottom:
                    continue
                rel_y = y - top
                bpl = self._bases_per_line()
                if rel_y < _WRAP_LINE_RULER + self._seq_row_height():
                    col = int((x - _LEFT_MARGIN) / cell_w)
                    if col < 0:
                        col = 0
                    bp = line * bpl + col + 1
                    return max(1, min(seq_len, bp))
                return max(1, min(seq_len, line * bpl + 1))
            return 1
        bp = int((x - origin_x + self._scroll_px) / cell_w) + 1
        return max(1, min(seq_len, bp))

    def _feature_track_top(self) -> int:
        if self.uses_wrap():
            return 0
        return self._sequence_top() + int(self._seq_row_height()) + _TRACK_GAP

    def _feature_at(self, pos: QPointF) -> int | None:
        if self.uses_wrap():
            cell_w, _, _ = self._cell_size()
            bpl = self._bases_per_line()
            for line in range(self._sequence_line_count()):
                segments = self._features_on_line(line)
                if not segments:
                    continue
                row_for_seg, _ = _pack_interval_rows(segments)
                track_base = self._line_top(line) + _WRAP_LINE_RULER + self._seq_row_height()
                row_h = _TRACK_HEIGHT + _TRACK_GAP
                for seg_i, (feat_idx, start_bp, end_bp) in enumerate(segments):
                    row = row_for_seg[seg_i]
                    y0 = track_base + row * row_h
                    col_start = start_bp - 1 - line * bpl
                    col_end = end_bp - 1 - line * bpl
                    x1 = _LEFT_MARGIN + col_start * cell_w
                    x2 = _LEFT_MARGIN + (col_end + 1) * cell_w
                    rect = QRectF(x1, y0, x2 - x1, _TRACK_HEIGHT)
                    if rect.contains(pos):
                        return feat_idx
            return None
        return None

    def _base_style(self, base: str) -> tuple[QColor, QColor]:
        if self._color_bases:
            return _base_colors(base, is_protein=self._is_protein)
        return _NO_COLOR_BG, _NO_COLOR_FG

    def _base_rect_at(self, index: int) -> QRectF:
        """Cell rect for one base; width can be fractional in bar mode."""
        origin_x, origin_y, cell_w, cell_h = self._nucleotide_cell_geometry()
        x = origin_x + index * cell_w - self._scroll_px
        if self.shows_letters():
            w = max(1.0, cell_w)
        else:
            w = max(0.05, cell_w)
        return QRectF(x, float(origin_y), w, float(cell_h))

    def _draw_base_cell(
        self,
        painter: QPainter,
        base: str,
        cell: QRectF,
        *,
        show_letters: bool,
    ) -> None:
        bg, fg = self._base_style(base)
        if self._color_bases:
            painter.fillRect(cell, bg)
            if show_letters:
                # thin separator so adjacent cells read as a clean uniform grid
                # (bright colors like green can otherwise look "off" next to
                # duller neighbors even though the cells are the same width)
                _draw_cell_separator(painter, cell)
        if show_letters:
            painter.save()
            painter.setClipRect(cell)
            painter.setFont(self._sequence_font_for_track(cell.height()))
            if self._color_bases:
                from repeatviewer.widgets.residue_color_grid import ResidueColorGrid

                painter.setPen(
                    QColor(ResidueColorGrid.contrast_text_color(bg))
                )
            else:
                painter.setPen(fg)
            painter.drawText(cell, Qt.AlignmentFlag.AlignCenter, base)
            painter.restore()

    def _draw_ruler(self, painter: QPainter, seq_len: int, ppb: float) -> None:
        painter.setPen(QColor("#666666"))
        ruler_font = _monospace_font(point_size=7)
        painter.setFont(ruler_font)
        fm = QFontMetrics(ruler_font)
        if self.uses_wrap():
            cell_w, _, _ = self._cell_size()
            bpl = self._bases_per_line()
            tick_step = self._nice_tick_step(bpl, cell_w, fm)
            for line in range(self._sequence_line_count()):
                line_start_bp = line * bpl + 1
                line_end_bp = min(seq_len, (line + 1) * bpl)
                ruler_y = self._line_top(line)
                seq_y = ruler_y + _WRAP_LINE_RULER
                last_label_x = -1e9
                for bp in range(line_start_bp, line_end_bp + 1):
                    idx = bp - 1
                    x, _ = self._index_to_xy(idx)
                    is_major = bp == line_start_bp or bp == line_end_bp or (bp - line_start_bp) % tick_step == 0
                    if not is_major:
                        continue
                    painter.drawLine(int(x), int(ruler_y + 2), int(x), int(seq_y - 2))
                    if x - last_label_x >= fm.horizontalAdvance("00000") + 4:
                        painter.drawText(int(x) + 2, int(ruler_y + fm.ascent() + 1), str(bp))
                        last_label_x = x
            return

        major_step = self._nice_tick_step(seq_len, ppb, fm)
        minor_step = max(1, major_step // 5) if major_step >= 10 else max(1, major_step // 2)
        min_label_w = fm.horizontalAdvance("00000") + 6
        last_label_x = -1e9

        for bp in range(1, seq_len + 1, minor_step):
            x = _LEFT_MARGIN + (bp - 1) * ppb
            is_major = bp == 1 or bp == seq_len or (bp - 1) % major_step == 0
            if bp == seq_len:
                continue
            if is_major:
                painter.drawLine(int(x), 0, int(x), _RULER_HEIGHT - 8)
            else:
                painter.drawLine(int(x), _RULER_HEIGHT - 10, int(x), _RULER_HEIGHT - 8)

        label_positions = list(range(1, seq_len + 1, major_step))
        if label_positions[-1] != seq_len:
            label_positions.append(seq_len)
        for bp in label_positions:
            x = _LEFT_MARGIN + (bp - 1) * ppb
            if x - last_label_x < min_label_w and bp != 1:
                continue
            painter.drawText(int(x) + 2, _RULER_HEIGHT - 6, str(bp))
            last_label_x = x

    def _sequence_viewport_x_extent(self) -> tuple[float, float]:
        """Viewport x range covered by bases 1..seq_len at the current scroll/zoom."""
        seq_len = len(self._sequence)
        if seq_len <= 0:
            return 0.0, 0.0
        origin_x, _, cell_w, _ = self._nucleotide_cell_geometry()
        scroll = self._scroll_px
        left = float(origin_x) - scroll
        right = float(origin_x) + seq_len * cell_w - scroll
        return left, right

    def _draw_sequence_bar_viewport(self, painter: QPainter, clip: QRectF) -> None:
        """Draw colored bars for the visible viewport (batched for long sequences)."""
        seq_len = len(self._sequence)
        origin_x, origin_y, cell_w, cell_h = self._nucleotide_cell_geometry()
        scroll = self._scroll_px
        if cell_w <= 0:
            return
        if not self._color_bases:
            painter.fillRect(
                QRectF(clip.left(), float(origin_y), clip.width(), float(cell_h)),
                _NO_COLOR_BG,
            )
            return

        seq_left, seq_right = self._sequence_viewport_x_extent()

        # Zoomed in (≈1 px/base or more): RLE by identical bases.
        if cell_w >= _BUCKET_DRAW_PPB:
            visible_start = max(0, int(math.floor((scroll + clip.left() - origin_x) / cell_w)))
            visible_end = min(
                seq_len,
                int(math.ceil((scroll + clip.right() - origin_x) / cell_w)),
            )
            i = visible_start
            while i < visible_end:
                base = self._sequence[i]
                bg, _ = self._base_style(base)
                run_start = i
                i += 1
                while i < visible_end and self._sequence[i] == base:
                    i += 1
                x1 = origin_x + run_start * cell_w - scroll
                x2 = origin_x + i * cell_w - scroll
                x1 = max(x1, seq_left)
                x2 = min(x2, seq_right)
                if x2 <= x1:
                    continue
                painter.fillRect(
                    QRectF(x1, float(origin_y), max(0.05, x2 - x1), float(cell_h)),
                    bg,
                )
            return

        # Fit / overview (sub-pixel bases): one sample per pixel, coalesce runs.
        # Always-RLE-by-base was the ~7 kb drag regression — DNA has ~1 base
        # per run, so it issued thousands of fillRects vs ~viewport-width here.
        px_left = max(0, int(math.floor(clip.left())), int(math.floor(seq_left)))
        px_right = min(
            self.width(),
            int(math.ceil(clip.right())) + 1,
            int(math.ceil(seq_right)),
        )
        if px_left >= px_right:
            return
        painter.setPen(Qt.PenStyle.NoPen)
        content_span = seq_len * cell_w
        px = px_left
        while px < px_right:
            content_x0 = scroll + px - origin_x
            if content_x0 < 0.0 or content_x0 >= content_span:
                px += 1
                continue
            bp0 = min(seq_len - 1, int(content_x0 / cell_w))
            bg, _ = self._base_style(self._sequence[bp0])
            run_start = px
            px += 1
            while px < px_right:
                content_x0 = scroll + px - origin_x
                if content_x0 < 0.0 or content_x0 >= content_span:
                    break
                bp1 = min(seq_len - 1, int(content_x0 / cell_w))
                if self._sequence[bp1] != self._sequence[bp0]:
                    break
                px += 1
            painter.fillRect(
                QRectF(float(run_start), float(origin_y), float(px - run_start), float(cell_h)),
                bg,
            )

    def _draw_sequence(self, painter: QPainter, clip: QRectF) -> None:
        if not self._sequence:
            return
        seq_len = len(self._sequence)
        show_letters = self.shows_letters()
        font = self._sequence_font()
        painter.setFont(font)

        if self.uses_wrap():
            letter_w, letter_h = self._letter_cell_size() if show_letters else (0, 0)
            for i in range(seq_len):
                x, y = self._index_to_xy(i)
                if show_letters:
                    cell = QRect(int(x), int(y), letter_w, letter_h)
                else:
                    _, _, cell_w, cell_h = self._nucleotide_cell_geometry()
                    cell = QRect(int(x), int(y), cell_w, cell_h)
                if not clip.intersects(QRectF(cell)):
                    continue
                self._draw_base_cell(
                    painter,
                    self._sequence[i],
                    cell,
                    show_letters=show_letters,
                )
            # Per-line complement / translation under each wrapped forward row.
            _, _, cell_w, _ = self._nucleotide_cell_geometry()
            if cell_w > 0 and self._extra_display_rows():
                for line in range(self._sequence_line_count()):
                    bpl = self._bases_per_line()
                    line_start = line * bpl
                    line_end = min(seq_len, line_start + bpl)
                    self._draw_extra_sequence_rows(
                        painter,
                        clip,
                        line_start,
                        line_end,
                        show_letters=show_letters,
                        origin_y_override=int(
                            self._line_top(line) + _WRAP_LINE_RULER
                        ),
                    )
            return

        origin_x, _, cell_w, _ = self._nucleotide_cell_geometry()
        if cell_w <= 0:
            return
        scroll = self._scroll_px
        visible_start = max(0, int(math.floor((scroll + clip.left() - origin_x) / cell_w)))
        visible_end = min(
            seq_len,
            int(math.ceil((scroll + clip.right() - origin_x) / cell_w)),
        )
        if not show_letters or visible_end - visible_start > 20_000:
            self._draw_sequence_bar_viewport(painter, clip)
            # Bars for complement/translation when zoomed out (height already reserved).
            if cell_w >= 1.0 and visible_end - visible_start <= 20_000:
                self._draw_extra_sequence_rows(
                    painter, clip, visible_start, visible_end, show_letters=False
                )
            return

        for i in range(visible_start, visible_end):
            cell = self._base_rect_at(i)
            self._draw_base_cell(
                painter,
                self._sequence[i],
                cell,
                show_letters=True,
            )
        self._draw_extra_sequence_rows(
            painter, clip, visible_start, visible_end, show_letters=True
        )

    def _draw_extra_sequence_rows(
        self,
        painter: QPainter,
        clip: QRectF,
        visible_start: int,
        visible_end: int,
        *,
        show_letters: bool,
        origin_y_override: int | None = None,
    ) -> None:
        """Complement and translation rows under the forward strand.

        Only the on-screen bases/codons are complemented or translated — never
        the full sequence (Geneious-style viewport work).
        """
        rows = self._extra_display_rows()
        if not rows:
            return
        origin_x, origin_y, cell_w, cell_h = self._nucleotide_cell_geometry()
        if origin_y_override is not None:
            origin_y = int(origin_y_override)
        if cell_w <= 0:
            return
        scroll = self._scroll_px
        seq = self._sequence
        gap = 2.0
        # Hairline under the forward strand so the complement row reads as its own track.
        sep_y = float(origin_y) + float(cell_h)
        sep = QRectF(clip.left(), sep_y, clip.width(), 1.0)
        if clip.intersects(sep):
            painter.fillRect(sep, QColor("#c8c8c8"))
        y = sep_y + gap
        letter_font = (
            self._sequence_font_for_track(float(cell_h)) if show_letters else None
        )

        from repeatviewer.backends.genetic_codes import amino_acid_display

        for kind in rows:
            if kind == "complement":
                for i in range(visible_start, visible_end):
                    # Same-position complement of the forward base (not a full RC).
                    letter = _complement_base(seq[i])
                    x = origin_x + i * cell_w - scroll
                    cell = QRectF(x, y, max(1.0, cell_w), float(cell_h))
                    if not clip.intersects(cell):
                        continue
                    bg, fg = _base_colors(letter, is_protein=False)
                    # Always fill — even if color-bases is off — so the second
                    # track cannot disappear into the paper background.
                    painter.fillRect(cell, bg if self._color_bases else QColor("#e8e8e8"))
                    if show_letters:
                        _draw_cell_separator(painter, cell)
                        if letter_font is not None:
                            painter.setFont(letter_font)
                            if self._color_bases:
                                from repeatviewer.widgets.residue_color_grid import (
                                    ResidueColorGrid,
                                )

                                painter.setPen(
                                    QColor(
                                        ResidueColorGrid.contrast_text_color(bg)
                                    )
                                )
                            else:
                                painter.setPen(QColor("#222222"))
                            painter.drawText(
                                cell, Qt.AlignmentFlag.AlignCenter, letter
                            )
            elif kind.startswith("frame:"):
                frame = int(kind.split(":", 1)[1])
                for start0, end0, residue in _iter_visible_frame_residues(
                    seq,
                    frame,
                    visible_start,
                    visible_end,
                    table=self._translation_table,
                ):
                    x = origin_x + start0 * cell_w - scroll
                    w = float(end0 - start0) * cell_w
                    cell = QRectF(x, y, max(1.0, w), float(cell_h))
                    if not clip.intersects(cell):
                        continue
                    bg, fg = _base_colors(residue, is_protein=True)
                    painter.fillRect(cell, bg if self._color_bases else QColor("#e8e8e8"))
                    if show_letters:
                        _draw_cell_separator(painter, cell)
                        if letter_font is not None:
                            label = amino_acid_display(
                                residue,
                                three_letter=self._translation_three_letter,
                            )
                            font = letter_font
                            if self._translation_three_letter and len(label) > 1:
                                # Fit Met / Ter into a 3-base-wide codon cell.
                                font = QFont(letter_font)
                                pt = max(5, int(letter_font.pointSizeF() * 0.72))
                                font.setPointSize(pt)
                            painter.setFont(font)
                            if self._color_bases:
                                from repeatviewer.widgets.residue_color_grid import (
                                    ResidueColorGrid,
                                )

                                painter.setPen(
                                    QColor(
                                        ResidueColorGrid.contrast_text_color(bg)
                                    )
                                )
                            else:
                                painter.setPen(QColor("#222222"))
                            painter.drawText(
                                cell, Qt.AlignmentFlag.AlignCenter, label
                            )
            y += float(cell_h) + gap

    def _draw_features(self, painter: QPainter) -> None:
        if not self._features:
            return
        ppb = self._pixels_per_bp
        cell_w, _, _ = self._cell_size()
        painter.setPen(Qt.PenStyle.NoPen)

        if self.uses_wrap():
            bpl = self._bases_per_line()
            for line in range(self._sequence_line_count()):
                segments = self._features_on_line(line)
                if not segments:
                    continue
                row_for_seg, _ = _pack_interval_rows(segments)
                track_base = self._line_top(line) + _WRAP_LINE_RULER + _SEQUENCE_ROW
                row_h = _TRACK_HEIGHT + _TRACK_GAP
                for seg_i, (feat_idx, start_bp, end_bp) in enumerate(segments):
                    row = row_for_seg[seg_i]
                    y = track_base + row * row_h
                    feat = self._features[feat_idx]
                    color = QColor(feat.color or feature_type_color(feat.feature_type))
                    col_start = start_bp - 1 - line * bpl
                    col_end = end_bp - 1 - line * bpl
                    x1 = _LEFT_MARGIN + col_start * cell_w
                    x2 = _LEFT_MARGIN + (col_end + 1) * cell_w
                    rect = QRectF(x1, y, max(2.0, x2 - x1), _TRACK_HEIGHT)
                    label = self._feature_track_label(feat, rect.width(), 40.0)
                    _draw_directed_feature(
                        painter,
                        rect,
                        color,
                        forward=feat.is_forward(),
                        label=label,
                    )
                    if feat_idx in self._selected_indices:
                        _draw_feature_selection_chrome(painter, rect)
            return

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        if not self._sequence:
            painter.fillRect(event.rect(), _MAP_PAPER)
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No sequence")
            painter.end()
            return

        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)

        # Blit cached sequence layer; rebuild only after zoom/scroll/etc.
        vector = self._vector_export
        if not vector:
            self._ensure_seq_layer_cache()
        if not vector and self._seq_layer_cache is not None:
            # IMPORTANT: do not pass event.rect() as the source on a HiDPI
            # QPixmap. Qt/PySide treats that source in *device* pixels, so only
            # the top half of the letter cache was stretched into the strip
            # ("half nucleotide" clipping). Draw the full logical pixmap instead.
            painter.drawPixmap(0, 0, self._seq_layer_cache)
        else:
            painter.fillRect(event.rect(), _MAP_PAPER)
            self._draw_sequence(painter, QRectF(event.rect()))
            if self.uses_wrap():
                self._draw_features(painter)
                track_bottom = self._feature_track_top() + self._feature_block_height()
                painter.setPen(QColor("#bbbbbb"))
                painter.drawLine(
                    _LEFT_MARGIN,
                    track_bottom + 2,
                    self._content_width() - _LEFT_MARGIN,
                    track_bottom + 2,
                )

        clip = QRectF(event.rect())
        if self._find_hits and not self.uses_wrap():
            self._paint_find_hits(painter, clip)
        self._paint_selection_overlay(painter, clip)
        painter.end()

    def _paint_selection_overlay(self, painter: QPainter, clip: QRectF) -> None:
        """Blue tint on selected nucleotides only (no grey veil on this strip)."""
        sel = self.selection_range()
        if not sel:
            return
        a, b = sel
        origin_x, origin_y, cell_w, cell_h = self._nucleotide_cell_geometry()
        if self.uses_wrap():
            cell_w_f, cell_h_f, _ = self._cell_size()
            span_h = float(cell_h_f + self._extra_display_height())
            seq_len = len(self._sequence)
            first = max(1, self._bp_at_pos(clip.left(), clip.top()) - 2)
            last = min(seq_len, self._bp_at_pos(clip.right(), clip.bottom()) + 2)
            for bp in range(max(a, first), min(b, last) + 1):
                x, y = self._index_to_xy(bp - 1)
                rect = QRectF(x, y, cell_w_f, span_h)
                if not clip.intersects(rect):
                    continue
                painter.fillRect(rect, _SELECTION_FILL)
            return

        # Nucleotide strip only — selected band blue; leave other bases untouched.
        top = float(origin_y)
        span_h = float(cell_h + self._extra_display_height())
        x1 = origin_x + (a - 1) * cell_w - self._scroll_px
        x2 = origin_x + b * cell_w - self._scroll_px
        rect = QRectF(x1, top, max(1.0, x2 - x1), span_h)
        if clip.intersects(rect):
            painter.fillRect(rect, _SELECTION_FILL)

    def _paint_find_hits(self, painter: QPainter, clip: QRectF) -> None:
        """Amber boxes for the Find hits that are on screen (bisect, not all hits)."""
        origin_x, origin_y, cell_w, cell_h = self._nucleotide_cell_geometry()
        if cell_w <= 0:
            return
        first_bp = int((clip.left() + self._scroll_px - origin_x) / cell_w) + 1 - self._find_max_len
        last_bp = int((clip.right() + self._scroll_px - origin_x) / cell_w) + 2
        lo = bisect_left(self._find_starts, max(1, first_bp))
        hi = bisect_right(self._find_starts, last_bp)
        top = float(origin_y)
        height = float(cell_h + self._extra_display_height())
        painter.setPen(QPen(_FIND_HIT_BORDER, 1))
        for a, b in self._find_hits[lo:hi]:
            x1 = origin_x + (a - 1) * cell_w - self._scroll_px
            x2 = origin_x + b * cell_w - self._scroll_px
            rect = QRectF(x1, top, max(2.0, x2 - x1), height)
            if clip.intersects(rect):
                painter.fillRect(rect, _FIND_HIT_FILL)
                painter.drawRect(rect)
        painter.setPen(Qt.PenStyle.NoPen)

    def _selection_paint_rect(self, sel: tuple[int, int]) -> QRect:
        """Pixel rect for a selection range (including highlight margin)."""
        a, b = sel
        origin_x, origin_y, cell_w, cell_h = self._nucleotide_cell_geometry()
        pad = 3
        if self.uses_wrap():
            cell_w_f, cell_h_f, _ = self._cell_size()
            span_h = float(cell_h_f + self._extra_display_height())
            x0, y0 = self._index_to_xy(a - 1)
            x1, y1 = self._index_to_xy(b - 1)
            left = min(x0, x1) - pad
            top = min(y0, y1) - pad
            right = max(x0 + cell_w_f, x1 + cell_w_f) + pad
            bottom = max(y0 + span_h, y1 + span_h) + pad
            return QRect(int(left), int(top), int(right - left), int(bottom - top)).intersected(
                self.rect()
            )
        x1 = origin_x + (a - 1) * cell_w - self._scroll_px
        x2 = origin_x + b * cell_w - self._scroll_px
        top = int(origin_y) - pad
        height = int(cell_h + self._extra_display_height()) + 2 * pad
        return QRect(
            int(min(x1, x2)) - pad,
            top,
            int(abs(x2 - x1)) + 2 * pad,
            height,
        ).intersected(self.rect())

    def _selection_change_damage(
        self, old: tuple[int, int] | None, new: tuple[int, int] | None
    ) -> QRect:
        """Damage for selection chrome.

        Range selections grey the outside, so the whole strip must refresh.
        Single-base carets still use the cheap edge-only path.
        """
        if old is None and new is None:
            return QRect()

        def _is_range(sel: tuple[int, int] | None) -> bool:
            return sel is not None and sel[0] != sel[1]

        if _is_range(old) or _is_range(new) or (old is not None and new is None):
            return self.rect()
        if old is None:
            return self._selection_paint_rect(new)  # type: ignore[arg-type]
        if new is None:
            return self._selection_paint_rect(old)
        oa, ob = old
        na, nb = new
        rect = QRect()
        if oa != na:
            rect = rect.united(self._selection_paint_rect((min(oa, na), max(oa, na))))
        if ob != nb:
            rect = rect.united(self._selection_paint_rect((min(ob, nb), max(ob, nb))))
        if not rect.isNull():
            return rect
        return self._selection_paint_rect(new)

    def _queue_selection_emit(self, start: int, end: int) -> None:
        self._pending_selection_emit = (int(start), int(end))
        if not self._selection_emit_timer.isActive():
            self._selection_emit_timer.start()

    def _flush_pending_selection_emit(self) -> None:
        pending = self._pending_selection_emit
        self._pending_selection_emit = None
        if pending is None:
            return
        self.selectionChanged.emit(pending[0], pending[1])

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() != Qt.MouseButton.LeftButton:
            super().mousePressEvent(event)
            return
        pos = event.position()
        if self.uses_wrap():
            feat_idx = self._feature_at(pos)
            if feat_idx is not None:
                self._selected_indices = {feat_idx}
                self._select_start_bp = None
                self._select_end_bp = None
                self.featureClicked.emit(feat_idx)
                self.setFocus(Qt.FocusReason.MouseFocusReason)
                self.update()
                return
        self._selected_indices.clear()
        bp = self._bp_at_pos(pos.x(), pos.y())
        old_sel = self.selection_range()
        # Always begin a new range here — even when clicking outside an existing
        # selection — so the user can drag another region without a clear click first.
        self._select_start_bp = bp
        self._select_end_bp = bp
        self._drag_select = True
        self.grabMouse()
        self._mouse_grabbed = True
        self.selectionChanged.emit(bp, bp)
        self.update(self._selection_change_damage(old_sel, (bp, bp)))
        event.accept()

    def mouseMoveEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if self._drag_select and self._select_start_bp is not None:
            if self._hover_track_tip:
                self._hover_track_tip = ""
                QToolTip.hideText()
            bp = self._bp_at_pos(event.position().x(), event.position().y())
            if bp == self._select_end_bp:
                event.accept()
                return
            old_sel = self.selection_range()
            self._select_end_bp = bp
            a, b = sorted((self._select_start_bp, bp))
            # Paint immediately; defer spin/status updates so drag stays smooth.
            self._queue_selection_emit(a, b)
            self.update(self._selection_change_damage(old_sel, (a, b)))
            event.accept()
            return
        self._update_sequence_track_tooltip(
            event.position(), event.globalPosition().toPoint()
        )
        super().mouseMoveEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802
        if self._hover_track_tip:
            self._hover_track_tip = ""
            QToolTip.hideText()
        super().leaveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._drag_select = False
            self._release_drag_mouse()
            self._selection_emit_timer.stop()
            self._flush_pending_selection_emit()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    @staticmethod
    def _nice_tick_step(seq_len: int, ppb: float, fm: QFontMetrics) -> int:
        sample_w = max(
            fm.horizontalAdvance("1"),
            fm.horizontalAdvance(str(max(1, seq_len))),
        )
        target_px = max(sample_w + 8, 56.0)
        # Do not floor ppb here — fit views for multi-Mb sequences are << 0.01.
        raw = max(1, int(target_px / max(float(ppb), _FIT_MIN_PPB)))
        for step in (
            1,
            2,
            5,
            10,
            20,
            50,
            100,
            200,
            500,
            1000,
            2000,
            5000,
            10_000,
            20_000,
            50_000,
            100_000,
            200_000,
            500_000,
            1_000_000,
            2_000_000,
        ):
            if step >= raw:
                return min(step, max(1, seq_len))
        return max(1, seq_len // 10)

    def contextMenuEvent(self, event: QContextMenuEvent) -> None:  # noqa: N802
        pos = QPointF(event.pos())
        feat_idx = self._feature_at(pos) if self.uses_wrap() else None
        if feat_idx is not None:
            self._selected_indices = {feat_idx}
            self._select_start_bp = None
            self._select_end_bp = None
            self.featureClicked.emit(feat_idx)
            self.update()
        elif self.selection_range() is None:
            bp = self._bp_at_pos(pos.x(), pos.y())
            if bp is None:
                return
            self._selected_indices.clear()
            self._select_start_bp = bp
            self._select_end_bp = bp
            self.selectionChanged.emit(bp, bp)
            self.update()

        menu = QMenu(self)
        menu.setObjectName("SequenceMapMenu")
        menu.setToolTipsVisible(False)

        if feat_idx is not None:
            menu.addAction("Delete").triggered.connect(
                lambda _checked=False: self.featureActionRequested.emit("delete")
            )
            menu.addAction("Copy sequence").triggered.connect(
                lambda _checked=False: self.featureActionRequested.emit("copy")
            )
            menu.addAction("Extract Sequence").triggered.connect(
                lambda _checked=False: self.featureActionRequested.emit("extract")
            )
            menu.addAction("Modify").triggered.connect(
                lambda _checked=False: self.featureActionRequested.emit("modify")
            )
            menu.addAction("Attach cytogenetic image…").triggered.connect(
                lambda _checked=False: self.featureActionRequested.emit("attach_image")
            )
            menu.addSeparator()

        add_action = menu.addAction("Add feature…")
        add_action.triggered.connect(
            lambda _checked=False: self.addFeatureRequested.emit("", "")
        )

        menu.addSeparator()
        selection = self.selection_range()
        copy_seq = menu.addAction("Copy Sequence")
        extract = menu.addAction("Extract Sequence")
        if selection is None:
            copy_seq.setEnabled(False)
            extract.setEnabled(False)
        else:
            start, end = selection
            copy_seq.triggered.connect(
                lambda _checked=False, a=start, b=end: self.copySequenceRequested.emit(a, b)
            )
            extract.triggered.connect(
                lambda _checked=False, a=start, b=end: self.extractSequenceRequested.emit(a, b)
            )

        menu.exec(event.globalPos())
        event.accept()

class _ExtractSequenceDialog(QDialog):
    """Name and adjust the range before extracting it as a new FASTA.

    The range starts at whatever was selected on the map but stays editable,
    so a selection made by dragging can be tidied up to exact coordinates
    without going back and re-selecting.

    When the Name still matches the auto ``{prefix}_{start}_{end}`` pattern,
    editing Start/End updates the Name so the output file tracks the range.
    """

    def __init__(
        self,
        *,
        sequence_id: str,
        sequence_length: int,
        start: int,
        end: int,
        default_name: str,
        reverse: bool = False,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Extract Sequence")
        self._sequence_length = max(1, int(sequence_length))

        self._start = QSpinBox(self)
        self._start.setRange(1, self._sequence_length)
        start_i = max(1, min(int(start), self._sequence_length))
        self._start.setValue(start_i)
        self._start.setGroupSeparatorShown(True)

        self._end = QSpinBox(self)
        self._end.setRange(1, self._sequence_length)
        end_i = max(1, min(int(end), self._sequence_length))
        self._end.setValue(end_i)
        self._end.setGroupSeparatorShown(True)

        for box in (self._start, self._end):
            box.setFixedWidth(160)

        self._strand = QComboBox(self)
        self._strand.addItem("Forward", False)
        self._strand.addItem("Reverse", True)
        self._strand.setCurrentIndex(1 if reverse else 0)
        self._strand.setFixedWidth(160)
        self._strand.setToolTip(
            "Reverse reverse-complements the extracted bases before writing the file"
        )

        self._name = QLineEdit(default_name, self)
        self._name.setFixedWidth(260)
        self._name.selectAll()
        # Auto-sync Name with Start/End when it still matches ``prefix_start_end``.
        self._name_prefix = self._prefix_from_ranged_name(default_name, start_i, end_i)
        self._last_auto_name = default_name if self._name_prefix is not None else ""
        self._name_user_edited = self._name_prefix is None

        self._summary = QLabel(self)
        self._summary.setStyleSheet("color: palette(mid);")

        source = QLabel(
            f"From {sequence_id or 'sequence'} ({self._sequence_length:,} bp)", self
        )
        source.setStyleSheet("color: palette(mid);")

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        form.addRow("Start:", self._start)
        form.addRow("End:", self._end)
        form.addRow("Strand:", self._strand)
        form.addRow("Name:", self._name)

        self._buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self._buttons.accepted.connect(self.accept)
        self._buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(12)
        layout.addWidget(source)
        layout.addLayout(form)
        layout.addWidget(self._summary)
        layout.addWidget(self._buttons)

        self._start.valueChanged.connect(self._on_range_changed)
        self._end.valueChanged.connect(self._on_range_changed)
        self._strand.currentIndexChanged.connect(self._sync)
        self._name.textEdited.connect(self._on_name_edited)
        self._name.textChanged.connect(self._sync)
        self._sync()

    @staticmethod
    def _ranged_name(prefix: str, start: int, end: int) -> str:
        return f"{prefix}_{int(start)}_{int(end)}"

    @classmethod
    def _prefix_from_ranged_name(
        cls, name: str, start: int, end: int
    ) -> str | None:
        text = str(name or "")
        suffix = f"_{int(start)}_{int(end)}"
        if text.endswith(suffix) and len(text) > len(suffix):
            return text[: -len(suffix)]
        return None

    def _on_name_edited(self, _text: str = "") -> None:
        self._name_user_edited = True

    def _on_range_changed(self, _value: int = 0) -> None:
        if (
            not self._name_user_edited
            and self._name_prefix is not None
            and self._name.text() == self._last_auto_name
        ):
            auto = self._ranged_name(
                self._name_prefix, self._start.value(), self._end.value()
            )
            self._name.blockSignals(True)
            self._name.setText(auto)
            self._name.blockSignals(False)
            self._last_auto_name = auto
        self._sync()

    def _sync(self) -> None:
        start = self._start.value()
        end = self._end.value()
        valid_range = end >= start
        has_name = bool(self._name.text().strip())
        if valid_range:
            strand = "reverse" if bool(self._strand.currentData()) else "forward"
            self._summary.setText(
                f"{end - start + 1:,} bp will be extracted ({strand})."
            )
        else:
            self._summary.setText("End must not be before start.")
        self._buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(
            valid_range and has_name
        )

    def values(self) -> tuple[int, int, str, bool]:
        return (
            self._start.value(),
            self._end.value(),
            self._name.text().strip(),
            bool(self._strand.currentData()),
        )


class SequenceMapWidget(QWidget):
    """Scrollable sequence map with overview and zoom controls."""

    featuresChanged = Signal()
    statusMessage = Signal(str)
    deleteFeaturesRequested = Signal()
    copyFeatureSequencesRequested = Signal()
    extractFeatureSequencesRequested = Signal()
    modifyFeatureRequested = Signal()
    changeFeatureColorRequested = Signal()
    annotationTracksMoveRequested = Signal(object)  # list[SequenceFeature]
    annotationFeaturesResizeRequested = Signal(object)  # list[SequenceFeature]
    attachFeatureImageRequested = Signal()
    openFeatureImagesRequested = Signal(int)  # feature index
    featureImagesDropped = Signal(int, list)  # feature index, image paths

    def __init__(self, parent: QWidget | None = None, *, embedded: bool = False) -> None:
        super().__init__(parent)
        self._embedded = embedded
        self._data: SequenceMapData | None = None
        self._dirty = False
        self._selected_feature_indices: set[int] = set()
        self._fit_retries = 0
        self._layout_intent_generation = 0
        self._show_similarity_database_name = True
        self._split_tracks_by_source = _load_split_tracks_by_source()
        self._show_online_blast = _load_show_online_blast()
        # None = no type filter; frozenset = Displayed annotations checkboxes.
        self._visible_annotation_kinds: frozenset[str] | None = None
        # Maps AnnotationTrackWidget indices → full ``_data.features`` indices.
        self._annot_display_to_full: list[int] = []
        self._transfer_preview_overlay: list[tuple[SequenceFeature, int]] = []
        if embedded:
            self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        else:
            self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        self._zoom_slider: QSlider | QSpinBox | None = None
        self._annotation_header_zoom_bar = None
        self._annotation_side_panel: QWidget | None = None
        self._toolbar: QWidget | None = None
        self._pinch_accum = 0.0
        zoom_row: QHBoxLayout | None = None

        # Shared position / font / stats controls (parented into toolbar or side panel).
        controls_parent = self
        self._sequence_stats_label = QLabel("", controls_parent)
        self._sequence_stats_label.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse
        )

        self._font_size = QSpinBox(controls_parent)
        self._font_size.setRange(_MIN_NUC_FONT_PT, _MAX_NUC_FONT_PT)
        self._font_size.setValue(_DEFAULT_NUC_FONT_PT)
        self._font_size.setSuffix(" pt")
        self._font_size.setFixedWidth(72)
        self._font_size.setToolTip("Nucleotide letter size (independent of zoom scale)")
        self._font_size.valueChanged.connect(self._on_font_size_changed)
        self._font_size.hide()

        self._pos_start = QSpinBox(controls_parent)
        self._pos_start.setMinimum(1)
        self._pos_start.setMaximum(1)
        self._pos_start.setKeyboardTracking(False)
        _configure_position_spin(self._pos_start)
        self._pos_end = QSpinBox(controls_parent)
        self._pos_end.setMinimum(1)
        self._pos_end.setMaximum(1)
        self._pos_end.setKeyboardTracking(False)
        _configure_position_spin(self._pos_end)
        self._pos_start.valueChanged.connect(self._on_position_spin_changed)
        self._pos_end.valueChanged.connect(self._on_position_spin_changed)
        self._pos_start.editingFinished.connect(self._on_position_editing_finished)
        self._pos_end.editingFinished.connect(self._on_position_editing_finished)
        select_btn = QPushButton("Select", controls_parent)
        select_btn.setToolTip(
            "Zoom to the selected range and center its midpoint in the window"
        )
        select_btn.clicked.connect(self._select_by_position)
        fit_btn = QPushButton("Fit", controls_parent)
        fit_btn.setToolTip("Fit the full sequence in the window")
        fit_btn.clicked.connect(self.fit_window)

        if embedded:
            # Horizontal toolbar matching the standalone annotation window layout:
            # Zoom | Position | Select | Fit | … | Length / GC
            toolbar = QWidget(self)
            toolbar.setObjectName("AnnotationMapToolbar")
            toolbar.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            toolbar.setMinimumWidth(0)
            toolbar.setFixedHeight(28)
            # Toolbar follows app chrome; sequence/annotation paper stays light.
            _apply_map_paper_background(toolbar, _app_chrome_background())
            tb = QHBoxLayout(toolbar)
            tb.setContentsMargins(6, 2, 8, 2)
            tb.setSpacing(6)
            tb.setAlignment(Qt.AlignmentFlag.AlignVCenter)
            tb.setSizeConstraint(QHBoxLayout.SizeConstraint.SetNoConstraint)

            zoom_out_btn = _make_zoom_button("−", toolbar)
            zoom_out_btn.setToolTip("Zoom out")
            zoom_out_btn.clicked.connect(self.zoom_out)
            zoom_in_btn = _make_zoom_button("+", toolbar)
            zoom_in_btn.setToolTip("Zoom in")
            zoom_in_btn.clicked.connect(self.zoom_in)
            self._zoom_slider = NoWheelSlider(Qt.Orientation.Horizontal, toolbar)
            self._zoom_slider.setObjectName("annotationZoomSlider")
            self._zoom_slider.setRange(0, 100)
            self._zoom_slider.setValue(0)
            self._zoom_slider.setFixedWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
            self._zoom_slider.setMinimumWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
            self._zoom_slider.setFixedHeight(ZOOM_SLIDER_HEIGHT)
            self._zoom_slider.setToolTip(
                "Click the bar: left = fit, right = maximum zoom. Drag to adjust."
            )
            apply_annotation_zoom_slider_style(self._zoom_slider)
            self._zoom_slider.valueChanged.connect(self._on_zoom_changed)

            for w in (self._pos_start, self._pos_end, select_btn, fit_btn):
                w.setParent(toolbar)
            self._sequence_stats_label.setParent(toolbar)
            _configure_position_spin(self._pos_start)
            _configure_position_spin(self._pos_end)
            for spin in (self._pos_start, self._pos_end):
                spin.setFixedHeight(18)
            # Match multi-seq Overview Select / Fit widths.
            select_btn.setFixedWidth(ANNOTATION_SELECT_BUTTON_RENDERED_WIDTH)
            select_btn.setFixedHeight(20)
            select_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
            fit_btn.setFixedWidth(ANNOTATION_FIT_BUTTON_RENDERED_WIDTH)
            fit_btn.setFixedHeight(20)
            fit_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

            pos_label = QLabel("Position:", toolbar)
            to_label = QLabel("to", toolbar)
            tb.addWidget(zoom_out_btn, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(zoom_in_btn, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(self._zoom_slider, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(pos_label, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(self._pos_start, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(to_label, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(self._pos_end, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(select_btn, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(fit_btn, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addStretch(1)
            self._sequence_stats_label.setAlignment(
                Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
            )
            self._sequence_stats_label.setWordWrap(False)
            tb.addWidget(self._sequence_stats_label, 0, Qt.AlignmentFlag.AlignVCenter)
            self._toolbar = toolbar
            from repeatviewer.widgets.toolbar_overflow import ToolbarOverflowHider

            self._toolbar_overflow = ToolbarOverflowHider(
                toolbar,
                [
                    [zoom_out_btn],
                    [zoom_in_btn],
                    [self._zoom_slider],
                    [pos_label, self._pos_start, to_label, self._pos_end],
                    [select_btn],
                    [fit_btn],
                    [self._sequence_stats_label],
                ],
                spacing=6,
                layout=tb,
            )
        else:
            zoom_row = QHBoxLayout()
            zoom_row.setContentsMargins(6, 4, 6, 4)
            zoom_row.setSpacing(6)
            zoom_row.setAlignment(Qt.AlignmentFlag.AlignVCenter)

            zoom_strip = QWidget(self)
            zoom_strip.setFixedHeight(_ZOOM_BTN_SIDE)
            zoom_strip.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            zoom_strip_layout = QHBoxLayout(zoom_strip)
            zoom_strip_layout.setContentsMargins(0, 0, 0, 0)
            zoom_strip_layout.setSpacing(4)
            zoom_strip_layout.setAlignment(Qt.AlignmentFlag.AlignVCenter)

            zoom_out_btn = _make_zoom_button("−", zoom_strip)
            zoom_out_btn.setToolTip("Zoom out")
            zoom_out_btn.clicked.connect(self.zoom_out)
            zoom_in_btn = _make_zoom_button("+", zoom_strip)
            zoom_in_btn.setToolTip("Zoom in")
            zoom_in_btn.clicked.connect(self.zoom_in)
            self._zoom_slider = NoWheelSlider(Qt.Orientation.Horizontal, zoom_strip)
            self._zoom_slider.setObjectName("annotationZoomSlider")
            self._zoom_slider.setRange(0, 100)
            self._zoom_slider.setValue(0)
            self._zoom_slider.setMinimumWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
            self._zoom_slider.setFixedHeight(ZOOM_SLIDER_HEIGHT)
            self._zoom_slider.setToolTip(
                "Click the bar: left = fit, right = maximum zoom. Drag to adjust."
            )
            self._zoom_slider.setSizePolicy(
                QSizePolicy.Policy.Expanding,
                QSizePolicy.Policy.Fixed,
            )
            apply_annotation_zoom_slider_style(self._zoom_slider)
            self._zoom_slider.valueChanged.connect(self._on_zoom_changed)

            zoom_strip_layout.addWidget(zoom_out_btn)
            zoom_strip_layout.addWidget(zoom_in_btn)
            zoom_strip_layout.addWidget(self._zoom_slider, 1)
            zoom_row.addWidget(zoom_strip, 1, Qt.AlignmentFlag.AlignVCenter)

            annotation_controls = QWidget(self)
            annotation_controls.setObjectName("AnnotationControls")
            annotation_controls.setFixedWidth(_ANNOTATION_SIDE_PANEL_WIDTH)
            controls_layout = QVBoxLayout(annotation_controls)
            controls_layout.setContentsMargins(8, 8, 8, 8)
            controls_layout.setSpacing(6)

            for w in (
                self._sequence_stats_label,
                self._pos_start,
                self._pos_end,
                select_btn,
                fit_btn,
            ):
                w.setParent(annotation_controls)

            position_label = QLabel("Position", annotation_controls)
            position_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
            controls_layout.addWidget(position_label)
            controls_layout.addWidget(self._pos_start)
            to_label = QLabel("to", annotation_controls)
            to_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
            controls_layout.addWidget(to_label)
            controls_layout.addWidget(self._pos_end)
            select_btn.setSizePolicy(
                QSizePolicy.Policy.Expanding,
                QSizePolicy.Policy.Fixed,
            )
            fit_btn.setSizePolicy(
                QSizePolicy.Policy.Expanding,
                QSizePolicy.Policy.Fixed,
            )
            controls_layout.addWidget(select_btn)
            controls_layout.addWidget(fit_btn)

            self._sequence_stats_label.setAlignment(
                Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop
            )
            self._sequence_stats_label.setWordWrap(True)
            controls_layout.addWidget(self._sequence_stats_label)
            controls_layout.addStretch()
            self._annotation_side_panel = annotation_controls
        self._canvas = SequenceMapCanvas(self)
        self._canvas.set_color_bases(True)
        self._canvas.set_split_tracks_by_source(self._split_tracks_by_source)
        self._canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._canvas.selectionChanged.connect(self._on_selection_changed)
        self._canvas.selectionCleared.connect(self._on_selection_cleared)
        self._canvas.featureClicked.connect(self._on_canvas_feature_clicked)
        self._canvas.featureActionRequested.connect(self._on_canvas_feature_action)

        self._ruler = ViewportRulerWidget(self)
        self._ruler.zoomAtRequested.connect(self._on_ruler_zoom_at)
        self._ruler.regionPreview.connect(self._on_ruler_region_preview)
        self._ruler.regionSelected.connect(self._on_ruler_region_selected)
        self._annot_track = AnnotationTrackWidget(self)
        self._annot_track.selectionChanged.connect(self._on_annot_track_selection_changed)
        self._annot_track.clearRangeSelectionRequested.connect(
            self._canvas.clear_selection
        )
        self._annot_track.deleteRequested.connect(self._request_delete_features)
        self._annot_track.copySequenceRequested.connect(self._request_copy_feature_sequences)
        self._annot_track.extractSequenceRequested.connect(
            self._request_extract_feature_sequences
        )
        self._annot_track.modifyRequested.connect(self._request_modify_feature)
        self._annot_track.changeColorRequested.connect(self._request_change_feature_color)
        self._annot_track.tracksMoved.connect(self._on_annot_tracks_moved)
        self._annot_track.featuresResized.connect(self._on_annot_features_resized)
        self._annot_track.attachImageRequested.connect(self._request_attach_feature_image)
        self._annot_track.openImagesRequested.connect(self._on_annot_open_images)
        self._annot_track.imagesDropped.connect(self._on_annot_images_dropped)
        self._overview = OverviewWidget(self)
        self._overview.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._overview.viewportChanged.connect(self._scroll_to_fraction)

        self._header = QWidget(self)
        self._header.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        header_layout = QVBoxLayout(self._header)
        header_layout.setContentsMargins(0, 0, 0, 0)
        header_layout.setSpacing(0)
        # Keep OverviewWidget for scroll/sync wiring, but do not show the minimap bar.
        self._overview.hide()
        header_layout.addWidget(self._ruler)

        self._canvas_row = QWidget(self)
        self._canvas_row.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        canvas_row_layout = QVBoxLayout(self._canvas_row)
        canvas_row_layout.setContentsMargins(0, 0, 0, 0)
        canvas_row_layout.setSpacing(0)
        # Canvas alone in this row — hbar is a sibling below so it can never
        # overlay and clip the bottom half of the nucleotide letters.
        canvas_row_layout.addWidget(self._canvas, 0)
        self._hbar = QScrollBar(Qt.Orientation.Horizontal, self)
        self._hbar.valueChanged.connect(self._on_horizontal_scroll)
        self._canvas.contentHeightChanged.connect(self._sync_header_views)

        self._overview_sync_timer = QTimer(self)
        self._overview_sync_timer.setSingleShot(True)
        self._overview_sync_timer.setInterval(_OVERVIEW_SYNC_MS)
        self._overview_sync_timer.timeout.connect(self._sync_overview_now)

        self._layout_sync_timer = QTimer(self)
        self._layout_sync_timer.setSingleShot(True)
        self._layout_sync_timer.setInterval(_LAYOUT_SYNC_MS)
        self._layout_sync_timer.timeout.connect(self._sync_header_views)

        self._sequence_block = QWidget(self)
        self._sequence_block.setObjectName("AnnotationSequenceBlock")
        self._sequence_block.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        _apply_map_paper_background(self._sequence_block, _ANNOT_PAPER)
        _apply_map_paper_background(self._canvas_row, _MAP_PAPER)
        _apply_map_paper_background(self._header, _MAP_PAPER)

        self._annot_scroll = QScrollArea(self._sequence_block)
        self._annot_scroll.setObjectName("AnnotationTrackScroll")
        self._annot_scroll.setWidget(self._annot_track)
        self._annot_scroll.setWidgetResizable(False)
        self._annot_scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._annot_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._annot_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._annot_scroll.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._annot_scroll.setFocusProxy(self._annot_track)
        # Keep a visible feature-track region even before annotations exist
        # (otherwise the expanding map looks like a blank/broken panel).
        self._annot_scroll.setMinimumHeight(_ANNOT_MIN_VIEWPORT_HEIGHT)
        _apply_map_paper_background(
            self._annot_scroll.viewport(), _ANNOT_PAPER
        )

        block_layout = QVBoxLayout(self._sequence_block)
        block_layout.setContentsMargins(0, 0, 0, 0)
        block_layout.setSpacing(0)
        block_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        block_layout.addWidget(self._canvas_row, 0, Qt.AlignmentFlag.AlignTop)
        block_layout.addWidget(self._hbar, 0)
        block_layout.addWidget(self._annot_scroll, 1)

        main_column = QWidget(self)
        main_column.setObjectName("AnnotationMapColumn")
        main_column.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._main_column = main_column
        _apply_map_paper_background(main_column, _MAP_PAPER)
        main_layout = QVBoxLayout(main_column)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        if self._toolbar is not None:
            main_layout.addWidget(self._toolbar)
        if zoom_row is not None:
            main_layout.addLayout(zoom_row)
        main_layout.addWidget(self._header)
        main_layout.addWidget(self._sequence_block, stretch=1)

        self._too_long_panel = QWidget(self)
        self._too_long_panel.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        too_long_layout = QVBoxLayout(self._too_long_panel)
        too_long_layout.setContentsMargins(8, 4, 8, 4)
        too_long_layout.setSpacing(0)
        self._too_long_label = QLabel("", self._too_long_panel)
        self._too_long_label.setWordWrap(True)
        self._too_long_label.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop
        )
        self._too_long_label.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse
        )
        self._too_long_label.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Minimum,
        )
        too_long_layout.addWidget(self._too_long_label, 0, Qt.AlignmentFlag.AlignTop)
        too_long_layout.addStretch(1)

        self._map_stack = QStackedWidget(self)
        self._map_stack.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._map_stack.addWidget(main_column)
        self._map_stack.addWidget(self._too_long_panel)
        self._too_long_active = False

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._map_stack, stretch=1)
        if self._annotation_side_panel is not None:
            layout.addWidget(self._annotation_side_panel)

        self._canvas.installEventFilter(self)
        self._canvas_row.installEventFilter(self)
        self._sequence_block.installEventFilter(self)
        self._annot_track.installEventFilter(self)
        self._annot_scroll.installEventFilter(self)
        self._annot_scroll.viewport().installEventFilter(self)
        self._ruler.installEventFilter(self)
        self._overview.installEventFilter(self)
        self._header.installEventFilter(self)
        # Needed for trackpad pinch (ZoomNativeGesture) to reach the filter.
        for _w in (
            self._canvas,
            self._canvas_row,
            self._sequence_block,
            self._annot_track,
            self._annot_scroll.viewport(),
            self._ruler,
            self._overview,
            self._header,
        ):
            _w.setAttribute(Qt.WidgetAttribute.WA_AcceptTouchEvents, True)

    def changeEvent(self, event) -> None:  # noqa: N802
        super().changeEvent(event)
        if event.type() in (
            QEvent.Type.PaletteChange,
            QEvent.Type.ApplicationPaletteChange,
            QEvent.Type.StyleChange,
        ):
            self._refresh_map_paper_backgrounds()

    def _refresh_map_paper_backgrounds(self) -> None:
        # Sequence + annotation tracks always stay on light paper.
        _apply_map_paper_background(self._sequence_block, _ANNOT_PAPER)
        _apply_map_paper_background(self._canvas_row, _MAP_PAPER)
        _apply_map_paper_background(self._header, _MAP_PAPER)
        _apply_map_paper_background(self._annot_scroll.viewport(), _ANNOT_PAPER)
        if getattr(self, "_main_column", None) is not None:
            _apply_map_paper_background(self._main_column, _MAP_PAPER)
        # Only the zoom/Position strip follows dark/light chrome.
        if getattr(self, "_toolbar", None) is not None:
            _apply_map_paper_background(self._toolbar, _app_chrome_background())
        self._canvas._invalidate_seq_layer_cache()
        self._canvas.update()
        self._annot_track.update()
        self._ruler.update()
        self._overview.update()

    def eventFilter(self, watched, event) -> bool:  # noqa: N802
        if event.type() == QEvent.Type.KeyPress and isinstance(event, QKeyEvent):
            if watched in (self._annot_track, self._annot_scroll, self._annot_scroll.viewport()):
                if _is_delete_key(event.key()) and self._annot_track.selected_indices():
                    self._sync_feature_selection_from_annot_track()
                    self._request_delete_features()
                    return True
        if event.type() == QEvent.Type.Resize and watched in (
            self._canvas,
            self._canvas_row,
            self._sequence_block,
            self._annot_scroll.viewport(),
        ):
            QTimer.singleShot(0, self._layout_sync_timer.start)

        zoom_surfaces = (
            self._canvas,
            self._canvas_row,
            self._annot_track,
            self._annot_scroll,
            self._annot_scroll.viewport(),
            self._ruler,
            self._overview,
            self._header,
            self._sequence_block,
        )

        # Trackpad pinch → zoom (same path as Preview / Dotplot).
        if (
            event.type() == QEvent.Type.NativeGesture
            and watched in zoom_surfaces
            and isinstance(event, QNativeGestureEvent)
            and event.gestureType() == Qt.NativeGestureType.ZoomNativeGesture
        ):
            if self._apply_pinch_zoom(float(event.value()), watched, event):
                event.accept()
                return True

        if event.type() == QEvent.Type.Wheel and watched in zoom_surfaces:
            wheel = event
            if isinstance(wheel, QWheelEvent):
                mods = wheel.modifiers() | QGuiApplication.keyboardModifiers()
                if mods & Qt.KeyboardModifier.ShiftModifier:
                    # Always consume so Shift+scroll never remaps to horizontal pan.
                    self._apply_shift_wheel_zoom(wheel, watched)
                    wheel.accept()
                    return True
                # Plain two-finger swipe: horizontal pan (not zoom).
                delta_x = _horizontal_wheel_delta(
                    wheel,
                    scroll_step=max(1, self._hbar.singleStep() * 3),
                )
                if delta_x != 0 and self._nudge_horizontal_scroll(delta_x):
                    wheel.accept()
                    return True
        return super().eventFilter(watched, event)

    def _event_to_canvas_x(self, watched, event) -> float:
        """Map an event position into SequenceMapCanvas local x."""
        if hasattr(event, "globalPosition"):
            try:
                global_pt = event.globalPosition().toPoint()
                return float(self._canvas.mapFromGlobal(global_pt).x())
            except Exception:
                pass
        if watched is self._canvas and hasattr(event, "position"):
            return float(event.position().x())
        return float(self._canvas.width()) * 0.5

    def _apply_shift_wheel_zoom(self, wheel: QWheelEvent, watched) -> bool:
        """Shift+scroll zoom: one ×1.08 step per event (MSA / multi-seq parity)."""
        if self._data is None or self.is_too_long_for_display():
            return False
        # Same delta extraction as MultiSequenceOverview / _BarsCanvas.
        delta = wheel.angleDelta().y()
        if delta == 0:
            delta = wheel.angleDelta().x()
        pixel = wheel.pixelDelta()
        if delta == 0 and not pixel.isNull():
            delta = int(pixel.y()) or int(pixel.x())
        if delta == 0:
            return False
        factor = _MAP_ZOOM_STEP if delta > 0 else 1.0 / _MAP_ZOOM_STEP
        self._zoom_by_at_view_x(factor, self._event_to_canvas_x(watched, wheel))
        return True

    def _apply_pinch_zoom(self, magnification: float, watched, event) -> bool:
        """Accumulate trackpad pinch magnification into discrete zoom steps."""
        if self._data is None or self.is_too_long_for_display():
            return False
        self._pinch_accum += float(magnification)
        steps = int(self._pinch_accum / _MAP_PINCH_DELTA_PER_STEP)
        if steps == 0:
            return True
        self._pinch_accum -= steps * _MAP_PINCH_DELTA_PER_STEP
        factor = _MAP_ZOOM_STEP ** steps
        self._zoom_by_at_view_x(factor, self._event_to_canvas_x(watched, event))
        return True

    def _zoom_by_at_view_x(self, factor: float, view_x: float) -> None:
        """Multiply current ppb by ``factor``, keeping the bp under ``view_x`` fixed."""
        if self._data is None or self.is_too_long_for_display():
            return
        if factor <= 0:
            return
        old_ppb = max(self._canvas.effective_ppb(), 1e-12)
        scroll = float(self._hbar.value())
        bp0 = (float(view_x) + scroll - float(_LEFT_MARGIN)) / old_ppb
        fit_ppb, max_ppb = self._slider_ppb_bounds()
        new_ppb = max(fit_ppb, min(max_ppb, old_ppb * float(factor)))
        if abs(new_ppb - old_ppb) < 1e-12:
            return
        self._apply_zoom_ppb(new_ppb, sync_slider=True, previous_ppb=old_ppb)
        new_ppb = max(self._canvas.effective_ppb(), 1e-12)
        new_content_x = float(_LEFT_MARGIN) + bp0 * new_ppb
        target = int(round(new_content_x - float(view_x)))
        self._hbar.setValue(max(0, min(self._hbar.maximum(), target)))
        self._sync_overview_now()

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._layout_sync_timer.start()
        if self._fit_retries > 0:
            QTimer.singleShot(0, self.fit_window)
        else:
            # Rail collapse/expand and window resize: re-fit when still at Fit zoom.
            QTimer.singleShot(0, self._ensure_fitted_to_viewport)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        QTimer.singleShot(0, self._sync_header_views)
        # First open often fits before the panel has its real width; re-check.
        QTimer.singleShot(0, self._ensure_fitted_to_viewport)

    def refresh_layout(self) -> None:
        self._sync_header_views()
        self._ensure_fitted_to_viewport()

    def _is_near_fit_zoom(self) -> bool:
        """True when the view is at (or past) whole-sequence Fit zoom."""
        if self._data is None or not self._data.sequence:
            return False
        if self._zoom_slider is not None and int(self._zoom_slider.value()) <= 1:
            return True
        fit_ppb, _ = self._slider_ppb_bounds()
        cur = float(self._canvas.effective_ppb())
        # After the viewport grows, current ppb is below the new fit target.
        return cur <= fit_ppb * 1.02

    def _sync_header_views(self) -> None:
        self._canvas.sync_viewport_from_widget()
        # Update scroll range/visibility before sizing the canvas row so the
        # horizontal bar gets a real height when zoomed in (not clipped to 0).
        self._update_horizontal_scroll()
        self._sync_sequence_block_layout()
        self._sync_overview_now()

    def _scroll_viewport_width(self) -> int:
        """Width available for sequence painting — never the full sequence pixel width."""
        for widget in (self._sequence_block, self._canvas, self):
            width = widget.width()
            if width > 20:
                return max(120, width - 2)
        parent = self.parentWidget()
        while parent is not None:
            if parent.width() > 100:
                return max(120, parent.width() - 48)
            parent = parent.parentWidget()
        return 120

    def _ensure_fitted_to_viewport(self) -> None:
        """Re-fit when Fit-zoom content no longer matches the viewport width."""
        if self._data is None or not self._data.sequence:
            return
        if not self.isVisible() or self.is_too_long_for_display():
            return
        viewport_w = self._scroll_viewport_width()
        panel_w = max(self.width(), self._sequence_block.width())
        if viewport_w <= 120 or (panel_w > 200 and viewport_w < panel_w * 0.45):
            self._fit_retries += 1
            if self._fit_retries < 40:
                QTimer.singleShot(50, self._ensure_fitted_to_viewport)
            return
        content_w = self._canvas.content_pixel_width()
        # At Fit zoom, keep the sequence edge-to-edge when the rail toggles or
        # the window resizes (content_w ≈ old viewport; new viewport differs).
        if (
            content_w > 0
            and self._is_near_fit_zoom()
            and (content_w < viewport_w * 0.98 or content_w > viewport_w * 1.02)
        ):
            self.fit_window()
            return
        # Tiny first-open fit: content fills only a sliver of the panel.
        if content_w > 0 and content_w < viewport_w * 0.55:
            self.fit_window()
            return
        self._fit_retries = 0

    def _horizontal_scroll_px(self) -> float:
        return float(self._hbar.value())

    def _update_horizontal_scroll(self) -> None:
        viewport_w = self._scroll_viewport_width()
        content_w = self._canvas.content_pixel_width()
        max_scroll = max(0, content_w - viewport_w)
        self._hbar.blockSignals(True)
        old = self._hbar.value()
        self._hbar.setRange(0, max_scroll)
        self._hbar.setPageStep(max(1, viewport_w // 2))
        cell_w = self._canvas.effective_cell_w()
        self._hbar.setSingleStep(max(1, int(cell_w * 5)))
        self._hbar.setValue(min(old, max_scroll))
        # Keep the bar in the layout whenever scrolling is possible so zoom-in
        # does not leave a 0-height groove under a fixed canvas row.
        need_bar = max_scroll > 0
        if self._hbar.isVisible() != need_bar:
            self._hbar.setVisible(need_bar)
        self._hbar.blockSignals(False)
        self._canvas.set_scroll_px(float(self._hbar.value()))

    def _nudge_horizontal_scroll(self, delta_x: int) -> bool:
        """Apply a horizontal pan; True when the scrollbar moved."""
        if delta_x == 0 or self._hbar.maximum() <= 0:
            return False
        new_val = max(0, min(self._hbar.maximum(), self._hbar.value() - delta_x))
        if new_val == self._hbar.value():
            return False
        self._hbar.setValue(new_val)
        return True

    def _sync_overview(self) -> None:
        self._overview_sync_timer.start()

    def _sync_overview_now(self) -> None:
        if self._data is None:
            return
        scroll_val = float(self._hbar.value())
        cell_w = self._canvas.effective_cell_w()
        vw = self._scroll_viewport_width()
        seq_len = self._data.length

        if cell_w > 0:
            start_bp = max(1, int((scroll_val - _LEFT_MARGIN) / cell_w) + 1)
            end_bp = min(seq_len, int((scroll_val + vw - _LEFT_MARGIN) / cell_w) + 1)
        else:
            start_bp = 1
            end_bp = seq_len
        display_features, _mapping, display_selected, feature_alphas = (
            self._annotation_display_bundle()
        )
        self._overview.set_data(
            seq_len,
            display_features,
            view_start_bp=start_bp,
            view_end_bp=end_bp,
            selected_indices=display_selected,
        )

        self._ruler.set_view(
            sequence=self._data.sequence,
            scroll_px=scroll_val,
            ppb=cell_w,
            viewport_w=vw,
        )

        feature_rows, _ = _pack_feature_rows(
            display_features,
            split_by_source=self._split_tracks_by_source,
        )
        self._annot_track.set_view(
            features=display_features,
            feature_rows=feature_rows,
            scroll_px=scroll_val,
            ppb=cell_w,
            viewport_w=vw,
            selected_indices=display_selected,
            feature_alphas=feature_alphas,
            seq_length=seq_len,
        )
        self._sync_sequence_block_layout()

    def _filter_display_by_annotation_kinds(
        self,
        display: list[SequenceFeature],
        mapping: list[int],
    ) -> tuple[list[SequenceFeature], list[int]]:
        """Apply Displayed annotations type filter (None = show all kinds)."""
        visible = self._visible_annotation_kinds
        if visible is None:
            return display, mapping
        kept_features: list[SequenceFeature] = []
        kept_mapping: list[int] = []
        for feat, full_idx in zip(display, mapping):
            if annotation_visibility_key(feat) in visible:
                kept_features.append(feat)
                kept_mapping.append(full_idx)
        return kept_features, kept_mapping

    def set_visible_annotation_kinds(
        self, keys: frozenset[str] | set[str] | None
    ) -> None:
        """Limit painted tracks to these visibility keys (None = all kinds)."""
        if keys is None:
            new_keys: frozenset[str] | None = None
        else:
            new_keys = frozenset(str(k) for k in keys)
        if new_keys == self._visible_annotation_kinds:
            return
        self._visible_annotation_kinds = new_keys
        if self._data is None:
            return
        display, mapping, display_selected, _alphas = self._annotation_display_bundle()
        # Nucleotide strip: real features only (no transfer-preview rows).
        canvas_features = [
            feat for feat, full_idx in zip(display, mapping) if full_idx >= 0
        ]
        self._canvas.set_sequence(
            self._data.sequence,
            canvas_features,
            is_protein=self._data.is_protein,
        )
        self._canvas.set_selected_features(display_selected)
        self._sync_overview()

    def visible_annotation_kinds(self) -> frozenset[str] | None:
        return self._visible_annotation_kinds

    def _annotation_display_bundle(
        self,
    ) -> tuple[list[SequenceFeature], list[int], set[int], list[int]]:
        """Visible features, display→full index map, selection, per-row alpha."""
        full = list(self._data.features) if self._data is not None else []
        display, mapping = filter_features_for_online_blast_display(
            full,
            show_online_blast=self._show_online_blast,
        )
        display, mapping = self._filter_display_by_annotation_kinds(display, mapping)
        alphas = [255] * len(display)
        preview = list(self._transfer_preview_overlay)
        if preview:
            for feat, alpha in preview:
                display.append(feat)
                mapping.append(-1)
                alphas.append(int(alpha))
            self._annot_display_to_full = mapping
        else:
            self._annot_display_to_full = mapping
        full_to_display = {
            full_idx: display_idx
            for display_idx, full_idx in enumerate(mapping)
            if full_idx >= 0
        }
        display_selected = {
            full_to_display[i]
            for i in self._selected_feature_indices
            if i in full_to_display
        }
        return display, mapping, display_selected, alphas

    def _map_annot_indices_to_full(self, indices: set[int]) -> set[int]:
        mapping = self._annot_display_to_full
        full: set[int] = set()
        for index in indices:
            if 0 <= int(index) < len(mapping):
                full.add(int(mapping[int(index)]))
        return full

    def _merge_display_features_into_full(
        self,
        display_features: list[SequenceFeature],
    ) -> list[SequenceFeature]:
        """Recombine a track-edited display list with features hidden from view."""
        from repeatviewer.backends.ncbi_blast_online import ONLINE_BLAST_SOURCE

        if self._data is None:
            return list(display_features)
        visible = self._visible_annotation_kinds
        hidden: list[SequenceFeature] = []
        for feat in self._data.features:
            source = str(feat.source or "")
            if not self._show_online_blast and source == ONLINE_BLAST_SOURCE:
                hidden.append(feat)
                continue
            if (
                visible is not None
                and annotation_visibility_key(feat) not in visible
            ):
                hidden.append(feat)
        return list(display_features) + hidden

    def _sync_sequence_block_layout(self) -> None:
        # Canvas row is exactly the nucleotide strip — scrollbar sits below it.
        canvas_h = max(1, int(self._canvas._content_height()))
        self._canvas.setMinimumHeight(canvas_h)
        self._canvas.setMaximumHeight(canvas_h)
        self._canvas.setFixedHeight(canvas_h)
        self._canvas_row.setMinimumHeight(canvas_h)
        self._canvas_row.setMaximumHeight(canvas_h)
        self._canvas_row.setFixedHeight(canvas_h)
        self._sync_annot_scroll_area()

    def _sync_annot_scroll_area(self) -> None:
        """Let annotation tracks fill space below the sequence; scroll when rows overflow."""
        self._annot_scroll.show()
        self._annot_scroll.setMinimumHeight(_ANNOT_MIN_VIEWPORT_HEIGHT)
        self._annot_scroll.setMaximumHeight(16777215)
        vw = max(1, self._scroll_viewport_width())
        vp_h = max(_ANNOT_MIN_VIEWPORT_HEIGHT, self._annot_scroll.viewport().height())
        # Grow the track to the viewport so empty space below features keeps the
        # same paper background as the annotation canvas (not the rail theme).
        content_h = max(1, self._annot_track.content_height())
        track_h = max(content_h, vp_h)
        self._annot_track.setFixedSize(vw, track_h)

    def _set_position_spin_limits(self, seq_len: int) -> None:
        seq_len = max(1, int(seq_len))
        for spin in (self._pos_start, self._pos_end):
            spin.setMaximum(seq_len)
            spin.setValue(min(spin.value(), seq_len) if spin.value() > 0 else 1)

    def _sync_position_spins(self, start: int, end: int) -> None:
        # Do not overwrite an in-progress typed value in either spin box.
        if self._pos_start.hasFocus() or self._pos_end.hasFocus():
            return
        a, b = sorted((int(start), int(end)))
        self._pos_start.blockSignals(True)
        self._pos_end.blockSignals(True)
        self._pos_start.setValue(a)
        self._pos_end.setValue(b)
        self._pos_start.blockSignals(False)
        self._pos_end.blockSignals(False)

    def _ensure_position_end_not_before_start(self) -> None:
        """If start > end, bump end to start+1 (clamped to sequence length)."""
        start = self._pos_start.value()
        end = self._pos_end.value()
        if start <= end:
            return
        seq_len = max(1, int(self._pos_end.maximum()))
        new_end = min(seq_len, start + 1)
        if new_end < start:
            new_end = start
        self._pos_end.blockSignals(True)
        self._pos_end.setValue(new_end)
        self._pos_end.blockSignals(False)

    def _apply_position_spins_to_selection(self) -> None:
        if self._data is None:
            return
        self._ensure_position_end_not_before_start()
        self._canvas.set_selection_range(self._pos_start.value(), self._pos_end.value())

    def _on_position_spin_changed(self, _value: int) -> None:
        # With keyboardTracking=False this runs after Enter/arrows/focus-out,
        # so typed multi-digit start values are not interrupted mid-entry.
        self._apply_position_spins_to_selection()

    def _on_position_editing_finished(self) -> None:
        self._apply_position_spins_to_selection()

    def _on_horizontal_scroll(self, value: int) -> None:
        self._canvas.set_scroll_px(float(value))
        self._overview_sync_timer.start()

    def _range_center_viewport_x(self, start: int, end: int) -> float:
        """Viewport x of the center of a bp range at the current scroll/zoom."""
        return self._canvas.bp_range_center_canvas_x(start, end)

    def _scroll_preserve_viewport_x(self, start: int, end: int, viewport_x: float) -> None:
        """Scroll so the range center stays at the given viewport x coordinate."""
        target = int(round(self._canvas.bp_range_center_content_x(start, end) - viewport_x))
        max_scroll = self._hbar.maximum()
        self._hbar.setValue(max(0, min(max_scroll, target)))
        self._sync_overview_now()

    def show_region(self, start: int, end: int, *, _attempt: int = 0) -> None:
        """Select [start, end] and bring it into view (Edit > Find, Go To Position).

        Unlike select_by_position, this keeps the current zoom when the region
        fits and the letters are readable; otherwise it zooms so the region is
        shown with some sequence on either side.
        """
        if self._data is None or self.is_too_long_for_display():
            return
        self._layout_intent_generation += 1
        viewport_w = self._scroll_viewport_width()
        if (not self.isVisible() or viewport_w <= 20) and _attempt < 40:
            QTimer.singleShot(
                50,
                lambda s=int(start), e=int(end), a=_attempt + 1: self.show_region(
                    s, e, _attempt=a
                ),
            )
            return
        a, b = sorted((int(start), int(end)))
        self._canvas.set_selection_range(a, b)
        self._sync_position_spins(a, b)
        span = max(1, b - a + 1)
        usable_w = max(1.0, float(viewport_w - _LEFT_MARGIN - _RIGHT_MARGIN))
        visible_bp = usable_w / max(self._canvas.pixels_per_bp(), 1e-9)
        if not self._canvas.shows_letters() or span > 0.9 * visible_bp:
            target_ppb = usable_w / max(3 * span, 60)  # the region with context
            letters_ppb = float(self._canvas._font_letter_cell_size()[0])
            if span * letters_ppb <= 0.9 * usable_w:
                target_ppb = max(target_ppb, letters_ppb)  # readable letters when it fits
            target_ppb = min(target_ppb, usable_w / span)  # the whole region stays on screen
            target_ppb = max(_FIT_MIN_PPB, min(self._canvas.zoom_ppb_max(), target_ppb))
            self._apply_zoom_ppb(target_ppb, sync_slider=True)
        self._scroll_to_bp_range(a, b)
        QTimer.singleShot(0, lambda lo=a, hi=b: self._scroll_to_bp_range(lo, hi))
        self.refresh_layout()

    def set_find_hits(self, hits: list[tuple[int, int]]) -> None:
        self._canvas.set_find_hits(hits)

    def select_by_position(self, start: int, end: int, *, _attempt: int = 0) -> None:
        if self._data is None or self.is_too_long_for_display():
            return
        # Marks this as the most recent explicit positioning request so a still-pending
        # fit_window() from the initial Annotation-tab auto-load (see load()/
        # _schedule_fit_window()) knows to back off instead of resetting a Go-to jump
        # back to whole-sequence view.
        self._layout_intent_generation += 1
        viewport_w = self._scroll_viewport_width()
        if (not self.isVisible() or viewport_w <= 20) and _attempt < 40:
            QTimer.singleShot(
                50,
                lambda s=int(start), e=int(end), a=_attempt + 1: self.select_by_position(
                    s, e, _attempt=a
                ),
            )
            return
        a, b = sorted((int(start), int(end)))
        self._canvas.set_selection_range(a, b)
        self._sync_position_spins(a, b)

        span = max(1, b - a + 1)
        usable_w = max(1.0, float(viewport_w - _LEFT_MARGIN - _RIGHT_MARGIN))
        if usable_w > 0:
            target_ppb = max(_FIT_MIN_PPB, min(self._canvas.zoom_ppb_max(), usable_w / span))
            self._apply_zoom_ppb(target_ppb, sync_slider=True)

        self._scroll_to_bp_range(a, b)
        QTimer.singleShot(0, lambda lo=a, hi=b: self._scroll_to_bp_range(lo, hi))
        self.refresh_layout()

    def _select_by_position(self) -> None:
        if self._data is None:
            return
        self._ensure_position_end_not_before_start()
        self.select_by_position(self._pos_start.value(), self._pos_end.value())

    def _scroll_to_bp_range(self, start: int, end: int) -> None:
        viewport_w = self._scroll_viewport_width()
        self._scroll_preserve_viewport_x(start, end, viewport_w / 2.0)

    def is_dirty(self) -> bool:
        return self._dirty

    def data(self) -> SequenceMapData | None:
        return self._data

    def is_too_long_for_display(self) -> bool:
        """True when Annotation is showing the too-long placeholder instead of the map."""
        return bool(getattr(self, "_too_long_active", False))

    def _set_too_long_display(self, active: bool, *, length: int = 0) -> None:
        """Swap between the map and the too-long message panel."""
        self._too_long_active = bool(active)
        if active:
            message = _annotation_too_long_message(length)
            self._too_long_label.setText(message)
            self._map_stack.setCurrentWidget(self._too_long_panel)
            if self._annotation_side_panel is not None:
                self._annotation_side_panel.hide()
            self.statusMessage.emit(message)
        else:
            self._map_stack.setCurrentIndex(0)
            if self._annotation_side_panel is not None:
                self._annotation_side_panel.show()

    def load(self, data: SequenceMapData, *, fit_after_load: bool = True) -> None:
        self._data = data
        self._dirty = False
        seq_len = int(data.length)
        if seq_len > _MAX_ANNOTATION_DISPLAY_LEN:
            # Keep metadata for multi-seq switching / summary, but do not paint.
            self._set_position_spin_limits(1)
            self._canvas.set_sequence("", [], is_protein=data.is_protein)
            self._canvas.clear_selection()
            self._annot_track.set_range_selection(None)
            self._sync_position_spins(1, 1)
            self._apply_sequence_stats_text("", None)
            self._set_too_long_display(True, length=seq_len)
            return

        self._set_too_long_display(False)
        self._set_position_spin_limits(data.length)
        if self._font_size.value() != _DEFAULT_NUC_FONT_PT:
            self._font_size.blockSignals(True)
            self._font_size.setValue(_DEFAULT_NUC_FONT_PT)
            self._font_size.blockSignals(False)
            self._canvas.set_font_point_size(_DEFAULT_NUC_FONT_PT)
        self._canvas.set_sequence(
            data.sequence,
            data.features,
            is_protein=data.is_protein,
        )
        self._canvas.clear_selection()
        self._annot_track.set_range_selection(None)
        self._sync_position_spins(1, 1)
        if seq_len > _INTERACTIVE_SEQUENCE_LEN:
            self.statusMessage.emit(
                f"Long sequence ({seq_len:,} bp): zoom in to see nucleotide letters "
                "(letters are drawn only for the visible region)."
            )
        if fit_after_load:
            self._schedule_fit_window()

    def _schedule_fit_window(self, delay_ms: int = 0) -> None:
        """Defer a fit_window() call, skipping it if select_by_position() sets a more
        specific view in the meantime (e.g. a Flexidot Go-to request racing the
        Annotation tab's own initial auto-load-and-fit)."""
        generation = self._layout_intent_generation

        def _maybe_fit() -> None:
            if self._layout_intent_generation == generation:
                self.fit_window()

        QTimer.singleShot(delay_ms, _maybe_fit)

    def features(self) -> list[SequenceFeature]:
        if self._data is None:
            return []
        return list(self._data.features)

    def set_transfer_preview(
        self, items: list[tuple[SequenceFeature, int]] | None
    ) -> None:
        self._transfer_preview_overlay = list(items or [])
        self._sync_overview()

    def clear_transfer_preview(self) -> None:
        self.set_transfer_preview(None)

    def set_features(self, features: list[SequenceFeature]) -> None:
        if self._data is None:
            return
        self._data.features = list(features)
        display, mapping = filter_features_for_online_blast_display(
            self._data.features,
            show_online_blast=self._show_online_blast,
        )
        display, mapping = self._filter_display_by_annotation_kinds(display, mapping)
        self._annot_display_to_full = mapping
        self._canvas.set_sequence(
            self._data.sequence,
            display,
            is_protein=self._data.is_protein,
        )
        valid = set(range(len(self._data.features)))
        self._selected_feature_indices &= valid
        self._dirty = True
        self.featuresChanged.emit()
        self._sync_overview()

    def set_show_online_blast(self, enabled: bool) -> None:
        enabled = bool(enabled)
        if enabled == self._show_online_blast:
            return
        self._show_online_blast = enabled
        _save_show_online_blast(enabled)
        if self._data is not None:
            display, mapping = filter_features_for_online_blast_display(
                self._data.features,
                show_online_blast=self._show_online_blast,
            )
            display, mapping = self._filter_display_by_annotation_kinds(
                display, mapping
            )
            self._annot_display_to_full = mapping
            self._canvas.set_sequence(
                self._data.sequence,
                display,
                is_protein=self._data.is_protein,
            )
        self._sync_overview()

    def show_online_blast(self) -> bool:
        return bool(self._show_online_blast)

    def reload_show_online_blast_pref(self) -> None:
        self.set_show_online_blast(_load_show_online_blast())

    def bind_annotation_zoom_bar(self, bar) -> None:
        """Attach header zoom controls when embedded in the preview panel.

        Embedded maps own a horizontal toolbar with zoom; ignore the header bar
        so it cannot steal the slider connection and leave the map unzoomable.
        """
        if self._embedded:
            return
        from repeatviewer.widgets.preview_controls import AnnotationZoomBar

        if not isinstance(bar, AnnotationZoomBar):
            return
        # Already bound — reconnecting would stack duplicate slots and make
        # every later zoom/tab switch progressively slower.
        if self._annotation_header_zoom_bar is bar and self._zoom_slider is bar.slider:
            return
        old = self._annotation_header_zoom_bar
        if old is not None:
            try:
                old.slider.valueChanged.disconnect(self._on_zoom_changed)
            except (RuntimeError, TypeError):
                pass
            try:
                old.fit_mode_changed.disconnect(self._on_annotation_header_fit_mode)
            except (RuntimeError, TypeError):
                pass
            try:
                old.reset_requested.disconnect(self.fit_window)
            except (RuntimeError, TypeError):
                pass
        self._annotation_header_zoom_bar = bar
        self._zoom_slider = bar.slider
        self._zoom_slider.valueChanged.connect(self._on_zoom_changed)
        bar.fit_mode_changed.connect(self._on_annotation_header_fit_mode)
        bar.reset_requested.connect(self.fit_window)

    def _on_annotation_header_fit_mode(self, mode: str) -> None:
        from repeatviewer.widgets.preview_zoom import PreviewFitMode

        if mode in (PreviewFitMode.FIT_WINDOW.value, PreviewFitMode.FIT_WIDTH.value):
            self.fit_window()
            return
        if mode == PreviewFitMode.ORIGINAL.value:
            self._apply_zoom_ppb(1.0, sync_slider=True)

    def mark_clean(self) -> None:
        self._dirty = False

    def zoom_in(self) -> None:
        if self._zoom_slider is None:
            return
        self._zoom_slider.setValue(min(100, self._zoom_slider.value() + 8))

    def zoom_out(self) -> None:
        if self._zoom_slider is None:
            return
        self._zoom_slider.setValue(max(0, self._zoom_slider.value() - 8))

    def _on_ruler_zoom_at(self, view_x: float, zoom_in: bool) -> None:
        """Scale-bar click: jump to max zoom (or fit) centered on the click."""
        if self._zoom_slider is None or self._data is None:
            return
        if self.is_too_long_for_display():
            return
        # Click → most zoom; right-click / Shift+click → fit (slider left).
        new_val = 100 if zoom_in else 0
        old_ppb = max(self._canvas.effective_ppb(), 1e-12)
        scroll = float(self._hbar.value())
        # 0-based bp under the click stays fixed in the viewport.
        bp0 = (float(view_x) + scroll - float(_LEFT_MARGIN)) / old_ppb
        new_ppb = self._slider_value_to_ppb(new_val)
        fit_ppb, _ = self._slider_ppb_bounds()
        reset_scroll = new_val <= 0 or new_ppb <= fit_ppb * 1.001
        self._apply_zoom_ppb(
            new_ppb,
            sync_slider=True,
            reset_scroll=reset_scroll,
            previous_ppb=old_ppb,
        )
        if reset_scroll:
            return
        new_content_x = float(_LEFT_MARGIN) + bp0 * max(
            self._canvas.effective_ppb(), 1e-12
        )
        target = int(round(new_content_x - float(view_x)))
        self._hbar.setValue(max(0, min(self._hbar.maximum(), target)))
        self._sync_overview_now()

    def _on_ruler_region_preview(self, start: int, end: int) -> None:
        """Live map + Position update while dragging on the scale."""
        if self._data is None or self.is_too_long_for_display():
            return
        self._canvas.set_selection_range(int(start), int(end))

    def _on_ruler_region_selected(self, start: int, end: int) -> None:
        """Commit a range selection from a scale drag."""
        if self._data is None or self.is_too_long_for_display():
            return
        self._canvas.set_selection_range(int(start), int(end))
        self._canvas.setFocus(Qt.FocusReason.MouseFocusReason)

    def _slider_ppb_bounds(self) -> tuple[float, float]:
        """Slider endpoints: left = fit entire sequence, right = max zoom-in."""
        if self._data is None or not self._data.sequence:
            return _FIT_MIN_PPB, float(_zoom_ppb_max())
        seq_len = self._data.length
        fit_ppb = self._ppb_to_fit_viewport(self._scroll_viewport_width(), seq_len)
        max_ppb = self._canvas.zoom_ppb_max()
        if max_ppb < fit_ppb:
            max_ppb = fit_ppb
        return fit_ppb, max_ppb

    def _slider_value_to_ppb(self, value: int) -> float:
        fit_ppb, max_ppb = self._slider_ppb_bounds()
        if max_ppb <= fit_ppb + 1e-15:
            return fit_ppb
        t = max(0, min(100, int(value))) / 100.0
        if max_ppb / fit_ppb > 50:
            log_fit = math.log(fit_ppb)
            log_max = math.log(max_ppb)
            return math.exp(log_fit + t * (log_max - log_fit))
        return fit_ppb + (max_ppb - fit_ppb) * t

    def _ppb_to_slider_value(self, ppb: float) -> int:
        fit_ppb, max_ppb = self._slider_ppb_bounds()
        ppb = max(fit_ppb, min(max_ppb, float(ppb)))
        if max_ppb <= fit_ppb + 1e-15:
            return 0
        if max_ppb / fit_ppb > 50:
            log_fit = math.log(fit_ppb)
            log_max = math.log(max_ppb)
            t = (math.log(ppb) - log_fit) / (log_max - log_fit)
        else:
            t = (ppb - fit_ppb) / (max_ppb - fit_ppb)
        return int(round(max(0.0, min(1.0, t)) * 100))

    def _ppb_to_fit_viewport(self, viewport_w: int, seq_len: int) -> float:
        """Scale so the full sequence fits exactly in the visible width."""
        if seq_len <= 0:
            return 1.0
        usable = max(1.0, float(viewport_w - _LEFT_MARGIN - _RIGHT_MARGIN))
        cap = self._canvas.zoom_ppb_max()
        return max(_FIT_MIN_PPB, min(cap, usable / seq_len))

    def _apply_zoom_ppb(
        self,
        ppb: float,
        *,
        sync_slider: bool = False,
        reset_scroll: bool = False,
        anchor_if_zoom_in: bool = False,
        previous_ppb: float | None = None,
    ) -> None:
        ppb = max(_FIT_MIN_PPB, min(self._canvas.zoom_ppb_max(), float(ppb)))
        prev_ppb = previous_ppb if previous_ppb is not None else self._canvas.effective_ppb()
        sel = self._canvas.selection_range()
        anchor_viewport_x = None
        if (
            anchor_if_zoom_in
            and sel is not None
            and sel[0] != sel[1]
            and ppb > prev_ppb
        ):
            anchor_viewport_x = self._range_center_viewport_x(sel[0], sel[1])

        if sync_slider:
            if self._zoom_slider is not None:
                slider_value = self._ppb_to_slider_value(ppb)
                self._zoom_slider.blockSignals(True)
                self._zoom_slider.setValue(slider_value)
                self._zoom_slider.blockSignals(False)
                header_bar = getattr(self, "_annotation_header_zoom_bar", None)
                if header_bar is not None:
                    header_bar.sync_display(slider_value)

        self._canvas.set_pixels_per_bp(ppb)
        if reset_scroll:
            self._hbar.setValue(0)
            self._canvas.set_scroll_px(0.0)
        self._sync_header_views()

        if anchor_viewport_x is not None and sel is not None:
            self._scroll_preserve_viewport_x(sel[0], sel[1], anchor_viewport_x)

    def fit_window(self) -> None:
        """Zoom out so the entire sequence fits in the visible sequence width."""
        if self._data is None or not self._data.sequence or self.is_too_long_for_display():
            return
        viewport_w = self._scroll_viewport_width()
        panel_w = max(self.width(), self._sequence_block.width())
        # Treat premature layout widths as not ready (not only <= 20).
        if viewport_w <= 120 or (panel_w > 240 and viewport_w < panel_w * 0.45):
            self._fit_retries += 1
            if self._fit_retries < 40:
                QTimer.singleShot(50, self.fit_window)
            return
        self._fit_retries = 0
        seq_len = len(self._data.sequence)
        ppb = self._ppb_to_fit_viewport(viewport_w, seq_len)
        # Keep any range highlight (including one opened from Dotplot → View).
        # Fit only changes zoom/scroll — same as Overview Fit.
        sel = self._canvas.selection_range()
        self._apply_zoom_ppb(ppb, sync_slider=True, reset_scroll=True)
        if sel is not None:
            self._annot_track.set_range_selection(sel[0], sel[1])
            self._ruler.set_region_selection(sel[0], sel[1])
            self._sync_position_spins(sel[0], sel[1])
        else:
            self._annot_track.set_range_selection(None)
            self._ruler.set_region_selection(None)
            self._sync_position_spins(1, 1)

    def fit_width(self) -> None:
        """Backward-compatible alias for fit_window."""
        self.fit_window()

    def refresh_nucleotide_colors(self) -> None:
        from repeatviewer.widgets.preview_zoom import invalidate_residue_color_maps

        invalidate_residue_color_maps()
        invalidate_base_qcolor_cache()
        self._canvas._invalidate_seq_layer_cache()
        self._canvas.update()

    def selection_range(self) -> tuple[int, int] | None:
        return self._canvas.selection_range()

    def view_state_for_session(self) -> dict | None:
        """Snapshot scroll/zoom/selection for session restore."""
        if self._data is None or self.is_too_long_for_display():
            return None
        sel = self.selection_range()
        has_selection = sel is not None
        if has_selection:
            start, end = int(sel[0]), int(sel[1])
        else:
            start, end = int(self._pos_start.value()), int(self._pos_end.value())
        return {
            "start": start,
            "end": end,
            "has_selection": bool(has_selection),
            "ppb": float(self._canvas.pixels_per_bp()),
            "scroll_px": float(self._canvas.scroll_px()),
            "sequence_id": str(self._data.sequence_id or ""),
        }

    def restore_view_state_for_session(
        self, state: dict, *, _attempt: int = 0
    ) -> None:
        """Reapply zoom/scroll from ``view_state_for_session`` (defers until visible).

        Restores the viewport the user left on, not the selection region. A saved
        selection is re-highlighted without changing zoom/scroll.
        """
        if self._data is None or self.is_too_long_for_display() or not state:
            return
        try:
            start = int(state.get("start", 1))
            end = int(state.get("end", 1))
        except (TypeError, ValueError):
            start, end = 1, 1
        seq_len = max(1, len(self._data.sequence))
        start = max(1, min(seq_len, start))
        end = max(1, min(seq_len, end))
        if end < start:
            start, end = end, start
        try:
            ppb = float(state.get("ppb") or 0.0)
        except (TypeError, ValueError):
            ppb = 0.0
        try:
            scroll_px = float(state.get("scroll_px") or 0.0)
        except (TypeError, ValueError):
            scroll_px = 0.0

        # Legacy snapshots without usable zoom: fall back to fitting start/end.
        if ppb <= 0:
            span = end - start + 1
            if span > 1 or (start == end and start > 1):
                self.select_by_position(start, end)
            return

        viewport_w = self._scroll_viewport_width()
        if (not self.isVisible() or viewport_w <= 20) and _attempt < 40:
            QTimer.singleShot(
                50,
                lambda s=dict(state), a=_attempt + 1: self.restore_view_state_for_session(
                    s, _attempt=a
                ),
            )
            return

        # Prevent a pending fit_window() from wiping the restored viewport.
        self._layout_intent_generation += 1
        self._apply_zoom_ppb(
            max(_FIT_MIN_PPB, min(self._canvas.zoom_ppb_max(), ppb)),
            sync_slider=True,
            reset_scroll=False,
        )
        self._update_horizontal_scroll()
        max_scroll = max(0, int(self._hbar.maximum()))
        scroll_px = max(0.0, min(float(max_scroll), scroll_px))
        self._hbar.setValue(int(scroll_px))
        self._canvas.set_scroll_px(scroll_px)
        self._sync_header_views()

        has_selection = state.get("has_selection")
        if has_selection is None:
            # Older saves: treat a non-trivial range as a selection to highlight.
            has_selection = (end - start + 1) > 1 or (start == end and start > 1)
        if has_selection:
            self._canvas.set_selection_range(start, end)
            self._sync_position_spins(start, end)
        else:
            self._sync_position_spins(start, end)
        self.refresh_layout()

    def selected_feature_index(self) -> int | None:
        return self._canvas.selected_feature_index()

    def selected_feature_indices(self) -> set[int]:
        self._sync_feature_selection_from_annot_track()
        return set(self._selected_feature_indices)

    def add_feature_from_selection(self, name: str, feature_type: str) -> bool:
        sel = self._canvas.selection_range()
        if sel is None or self._data is None:
            return False
        start, end = sel
        feat = SequenceFeature(
            start=start,
            end=end,
            name=name,
            feature_type=feature_type,
            source=REPEATVIEWER_INTERNAL_SOURCE,
        )
        features = list(self._data.features)
        features.append(feat)
        self.set_features(features)
        self.set_selected_features({len(features) - 1})
        self._sync_overview()
        return True

    def delete_selected_features(self) -> int:
        self._sync_feature_selection_from_annot_track()
        indices = self._selected_feature_indices
        if not indices or self._data is None:
            return 0
        features = [feat for i, feat in enumerate(self._data.features) if i not in indices]
        removed = len(self._data.features) - len(features)
        if removed <= 0:
            return 0
        self.set_features(features)
        self.clear_feature_selection()
        self._sync_overview()
        return removed

    def delete_selected_feature(self) -> bool:
        return self.delete_selected_features() > 0

    def clear_feature_selection(self) -> None:
        self._selected_feature_indices.clear()
        self._canvas.clear_feature_selection()
        self._annot_track._selected_indices = set()
        self._annot_track.update()
        self._sync_overview()
        self._sync_header_views()

    def set_show_similarity_database_name(self, show: bool) -> None:
        show = bool(show)
        if show == self._show_similarity_database_name:
            return
        self._show_similarity_database_name = show
        self._canvas.set_show_similarity_database_name(show)
        self._annot_track.set_show_similarity_database_name(show)
        self.update()

    def set_show_complement(self, enabled: bool) -> None:
        self._canvas.set_show_complement(enabled)
        self._sync_header_views()
        self._canvas._invalidate_seq_layer_cache()
        self._canvas.update()

    def set_translation_frames(self, frames: object) -> None:
        self._canvas.set_translation_frames(frames)
        self._sync_header_views()
        self._canvas._invalidate_seq_layer_cache()
        self._canvas.update()

    def set_translation_table(self, table_id: int) -> None:
        self._canvas.set_translation_table(table_id)
        self._canvas._invalidate_seq_layer_cache()
        self._canvas.update()

    def set_translation_three_letter(self, enabled: bool) -> None:
        self._canvas.set_translation_three_letter(enabled)
        self._canvas._invalidate_seq_layer_cache()
        self._canvas.update()

    def set_letter_min_px(self, pixels: int) -> None:
        self._canvas.set_letter_min_px(pixels)
        self._sync_header_views()

    def set_annotation_label_size(self, point_size: int) -> None:
        self._annot_track.set_label_point_size(point_size)
        self._sync_header_views()

    def set_annotation_thickness(self, height: int) -> None:
        self._annot_track.set_feature_bar_height(height)
        self._sync_header_views()

    def set_annotation_distance(self, gap: int) -> None:
        self._annot_track.set_track_gap(gap)
        self._sync_header_views()

    def set_base_number_size(self, point_size: int) -> None:
        self._ruler.set_font_point_size(point_size)
        self._sync_header_views()

    def set_split_tracks_by_source(self, enabled: bool) -> None:
        enabled = bool(enabled)
        if enabled == self._split_tracks_by_source:
            return
        self._split_tracks_by_source = enabled
        self._canvas.set_split_tracks_by_source(enabled)
        self._sync_overview()
        self.update()

    def split_tracks_by_source(self) -> bool:
        return bool(self._split_tracks_by_source)

    def reload_split_tracks_by_source_pref(self) -> None:
        """Reload the QSettings preference and re-pack display rows."""
        self.set_split_tracks_by_source(_load_split_tracks_by_source())

    def font_point_size(self) -> int:
        return int(self._font_size.value())

    def set_font_point_size(self, point_size: int) -> None:
        size = max(_MIN_NUC_FONT_PT, min(_MAX_NUC_FONT_PT, int(point_size)))
        if self._font_size.value() != size:
            self._font_size.blockSignals(True)
            self._font_size.setValue(size)
            self._font_size.blockSignals(False)
        self._canvas.set_font_point_size(size)

    def _apply_sequence_stats_text(self, length_part: str, gc_part: str | None) -> None:
        length_part = length_part.strip()
        if not length_part:
            self._sequence_stats_label.clear()
            self._sequence_stats_label.hide()
            return
        text = length_part
        if gc_part:
            gc_part = gc_part.strip()
            if gc_part:
                # Embedded toolbar is a single horizontal row (match reference UI).
                sep = "  " if self._embedded else "\n"
                text = f"{length_part}{sep}{gc_part}"
        self._sequence_stats_label.setText(text)
        self._sequence_stats_label.show()

    def set_sequence_summary_text(self, text: str) -> None:
        text = text.strip()
        if not text:
            self._apply_sequence_stats_text("", None)
            return
        if " | " in text:
            length_part, gc_part = text.split(" | ", 1)
            self._apply_sequence_stats_text(length_part, gc_part)
            return
        self._apply_sequence_stats_text(text, None)

    def set_sequence_summary(self, summary: SingleSequenceSummary) -> None:
        length_part, gc_part = format_single_sequence_summary_parts(summary)
        self._apply_sequence_stats_text(length_part, gc_part)

    def _on_font_size_changed(self, point_size: int) -> None:
        self._canvas.set_font_point_size(point_size)
        self._sync_header_views()

    def _on_zoom_changed(self, value: int) -> None:
        ppb = self._slider_value_to_ppb(value)
        fit_ppb, _ = self._slider_ppb_bounds()
        previous_ppb = self._canvas.effective_ppb()
        reset_scroll = value <= 0 or ppb <= fit_ppb * 1.001
        self._apply_zoom_ppb(
            ppb,
            anchor_if_zoom_in=not reset_scroll,
            previous_ppb=previous_ppb,
            reset_scroll=reset_scroll,
        )

    def _on_selection_changed(self, start: int, end: int) -> None:
        self._sync_position_spins(start, end)
        self._annot_track.set_range_selection(start, end)
        self._ruler.set_region_selection(start, end)
        if start == end:
            self.statusMessage.emit(f"Position: {start}")
            return
        a, b = sorted((start, end))
        unit = "aa" if (self._data is not None and self._data.is_protein) else "bp"
        self.statusMessage.emit(f"Selection: {a} – {b} ({b - a + 1} {unit})")

    def _on_selection_cleared(self) -> None:
        self._annot_track.set_range_selection(None)
        self._ruler.set_region_selection(None)
        self._sync_position_spins(1, 1)
        self.statusMessage.emit("Selection cleared")

    def set_selected_features(
        self,
        indices: set[int] | frozenset[int] | None,
        *,
        expand_pairs: bool = True,
    ) -> None:
        selected = set(indices or [])
        if expand_pairs and self._data is not None:
            selected = expand_terminal_pair_indices(self._data.features, selected)
        self._selected_feature_indices = selected
        self._canvas.set_selected_features(self._selected_feature_indices)
        # Annot track may be showing a filtered list — map full → display indices.
        _display, _mapping, display_selected, _alphas = self._annotation_display_bundle()
        self._annot_track._selected_indices = set(display_selected)
        self._annot_track.update()
        self._sync_overview()

    def _on_canvas_feature_clicked(self, index: int) -> None:
        self.set_selected_features({index})
        self._on_feature_clicked(index)

    def _on_annot_track_selection_changed(self, indices: object) -> None:
        if isinstance(indices, frozenset):
            selected = set(indices)
        elif isinstance(indices, set):
            selected = set(indices)
        else:
            selected = set()
        selected = self._map_annot_indices_to_full(selected)
        if self._data is not None:
            selected = expand_terminal_pair_indices(self._data.features, selected)
        self._selected_feature_indices = selected
        self._canvas.set_selected_features(selected)
        _display, _mapping, display_selected, _alphas = self._annotation_display_bundle()
        self._annot_track._selected_indices = set(display_selected)
        self._annot_track.update()
        self._sync_overview()
        if len(selected) == 1:
            self._on_feature_clicked(next(iter(selected)))
        elif selected and self._data:
            starts: list[int] = []
            ends: list[int] = []
            for index in selected:
                if 0 <= index < len(self._data.features):
                    feat = self._data.features[index]
                    starts.append(feat.start)
                    ends.append(feat.end)
            if starts:
                lo, hi = min(starts), max(ends)
                self._sync_position_spins(lo, hi)
                self._canvas.set_selection_range(lo, hi)
            self.statusMessage.emit(f"{len(selected)} features selected")

    def _on_feature_clicked(self, index: int) -> None:
        if self._data and 0 <= index < len(self._data.features):
            selected = expand_terminal_pair_indices(
                self._data.features, {index}
            )
            if len(selected) > 1:
                starts = [
                    self._data.features[i].start
                    for i in selected
                    if 0 <= i < len(self._data.features)
                ]
                ends = [
                    self._data.features[i].end
                    for i in selected
                    if 0 <= i < len(self._data.features)
                ]
                if starts and ends:
                    lo, hi = min(starts), max(ends)
                    self._sync_position_spins(lo, hi)
                    self._canvas.set_selection_range(lo, hi)
                names = ", ".join(
                    self._data.features[i].name
                    for i in sorted(selected)
                    if 0 <= i < len(self._data.features)
                )
                self.statusMessage.emit(f"Features: {names}")
                return
            feat = self._data.features[index]
            self._sync_position_spins(feat.start, feat.end)
            self._canvas.set_selection_range(feat.start, feat.end)
            self.statusMessage.emit(f"Feature: {feat.name} ({feat.start}–{feat.end})")

    def _sync_feature_selection_from_annot_track(self) -> None:
        """Prefer annot-track selection; keep existing map selection if track is empty."""
        annot_sel = self._annot_track.selected_indices()
        if annot_sel:
            mapped = self._map_annot_indices_to_full(annot_sel)
            if mapped != self._selected_feature_indices:
                self._selected_feature_indices = mapped
                self._canvas.set_selected_features(mapped)
                self._sync_overview()
            return
        # Annot track empty — keep canvas/map selection and mirror it onto the track.
        if self._selected_feature_indices:
            _display, _mapping, display_selected, _alphas = self._annotation_display_bundle()
            self._annot_track._selected_indices = set(display_selected)

    def _request_delete_features(self) -> None:
        self._sync_feature_selection_from_annot_track()
        if self._selected_feature_indices:
            self.deleteFeaturesRequested.emit()

    def _request_copy_feature_sequences(self) -> None:
        self._sync_feature_selection_from_annot_track()
        if self._selected_feature_indices:
            self.copyFeatureSequencesRequested.emit()

    def _request_extract_feature_sequences(self) -> None:
        self._sync_feature_selection_from_annot_track()
        if self._selected_feature_indices:
            self.extractFeatureSequencesRequested.emit()

    def _request_modify_feature(self) -> None:
        self._sync_feature_selection_from_annot_track()
        if self._selected_feature_indices:
            self.modifyFeatureRequested.emit()

    def _request_change_feature_color(self) -> None:
        self._sync_feature_selection_from_annot_track()
        if len(self._selected_feature_indices) > 1:
            self.changeFeatureColorRequested.emit()

    def _on_annot_open_images(self, display_index: int) -> None:
        mapped = self._map_annot_indices_to_full({int(display_index)})
        if mapped:
            self.openFeatureImagesRequested.emit(next(iter(mapped)))

    def _on_annot_images_dropped(self, display_index: int, paths: list) -> None:
        mapped = self._map_annot_indices_to_full({int(display_index)})
        if mapped:
            self.featureImagesDropped.emit(next(iter(mapped)), paths)

    def _on_annot_tracks_moved(self, features: object) -> None:
        updated = list(features or [])
        if not updated:
            return
        self.annotationTracksMoveRequested.emit(
            self._merge_display_features_into_full(updated)
        )

    def _on_annot_features_resized(self, features: object) -> None:
        updated = list(features or [])
        if not updated:
            return
        self.annotationFeaturesResizeRequested.emit(
            self._merge_display_features_into_full(updated)
        )

    def _on_canvas_feature_action(self, action: str) -> None:
        kind = str(action or "").strip().lower()
        if kind == "delete":
            self._request_delete_features()
        elif kind == "copy":
            self.copyFeatureSequencesRequested.emit()
        elif kind == "extract":
            self.extractFeatureSequencesRequested.emit()
        elif kind == "modify":
            self.modifyFeatureRequested.emit()
        elif kind == "attach_image":
            self.attachFeatureImageRequested.emit()

    def _request_attach_feature_image(self) -> None:
        self._sync_feature_selection_from_annot_track()
        self.attachFeatureImageRequested.emit()

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802
        if _is_delete_key(event.key()):
            self._sync_feature_selection_from_annot_track()
            if self._selected_feature_indices:
                self.deleteFeaturesRequested.emit()
                event.accept()
                return
        super().keyPressEvent(event)

    def _scroll_to_fraction(self, fraction: float) -> None:
        """Scroll so the overview viewport starts at ``fraction`` of the sequence.

        ``OverviewWidget`` emits the start of the visible window as a fraction of
        sequence length in ``[0, 1 - view_fraction]``. That must be stretched onto
        the full horizontal scrollbar range ``[0, max_scroll]``; multiplying by
        ``max_scroll`` alone left the rightmost bases unreachable.
        """
        if self._data is None:
            return
        max_scroll = self._hbar.maximum()
        if max_scroll <= 0:
            return
        view_frac = max(0.0, min(1.0, float(self._overview._view_fraction)))
        usable = max(1e-9, 1.0 - view_frac)
        t = max(0.0, min(1.0, float(fraction) / usable))
        self._hbar.setValue(int(round(t * max_scroll)))
        self._sync_overview_now()


class _AddFeatureColorField(QWidget):
    """Color swatch matching native combo row size."""

    def __init__(self, initial_color: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._color = QColor(initial_color)
        self._button = QPushButton(self)
        self._button.setCursor(Qt.CursorShape.PointingHandCursor)
        self._button.setToolTip("Choose feature color")
        self._button.clicked.connect(self._pick_color)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._button)
        self._refresh_button()

    def _refresh_button(self) -> None:
        hex_color = self._color.name(QColor.NameFormat.HexRgb)
        self._button.setFlat(False)
        self._button.setStyleSheet(f"background-color: {hex_color};")

    def _pick_color(self) -> None:
        chosen = QColorDialog.getColor(self._color, self, "Feature color")
        if not chosen.isValid():
            return
        self._color = chosen
        self._refresh_button()

    def set_field_size(self, width: int, height: int) -> None:
        self.setFixedSize(width, height)
        self._button.setFixedSize(width, height)

    def color_hex(self) -> str:
        return self._color.name(QColor.NameFormat.HexRgb)

    def set_color_hex(self, color: str) -> None:
        self._color = QColor(color)
        self._refresh_button()


class AddFeatureDialog(QDialog):
    """Add a sequence feature with uniform native macOS form fields."""

    _POSITION_GAP = 8

    def __init__(
        self,
        *,
        start: int,
        end: int,
        seq_length: int,
        parent: QWidget | None = None,
        feature_type: str | None = None,
        feature_label: str | None = None,
        default_name: str = "",
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Add feature")

        start_pos, end_pos = sorted((int(start), int(end)))

        self._type = QComboBox(self)
        for label, ftype in load_add_feature_buttons():
            self._type.addItem(label, ftype)
        if feature_type is not None:
            index = self._type.findData(feature_type)
            if index < 0 and feature_label:
                index = self._type.findText(feature_label)
            if index >= 0:
                self._type.setCurrentIndex(index)

        self._strand = QComboBox(self)
        self._strand.addItem("Forward (+)", "+")
        self._strand.addItem("Reverse (−)", "-")

        self._start = QSpinBox(self)
        self._start.setRange(1, seq_length)
        self._start.setValue(start_pos)

        self._end = QSpinBox(self)
        self._end.setRange(1, seq_length)
        self._end.setValue(end_pos)

        initial_type = str(self._type.currentData() or self._type.currentText())
        self._color = _AddFeatureColorField(feature_type_color(initial_type), self)

        self._name = QLineEdit(self)
        self._name.setPlaceholderText("e.g. LTR5")
        if default_name:
            self._name.setText(default_name)

        self._position_row = QWidget(self)
        position_layout = QHBoxLayout(self._position_row)
        position_layout.setContentsMargins(0, 0, 0, 0)
        position_layout.setSpacing(6)
        position_layout.addWidget(self._start)
        self._end_label = QLabel("End:", self._position_row)
        self._end_label.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        position_layout.addWidget(self._end_label)
        position_layout.addWidget(self._end)

        form = QFormLayout()
        form.addRow("Type:", self._type)
        form.addRow("Direction:", self._strand)
        form.addRow("Start:", self._position_row)
        form.addRow("Color:", self._color)
        form.addRow("Name:", self._name)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

        self._type.currentIndexChanged.connect(self._on_type_changed)

    def _on_type_changed(self) -> None:
        ftype = str(self._type.currentData() or self._type.currentText())
        self._color.set_color_hex(feature_type_color(ftype))

    def feature(self) -> SequenceFeature | None:
        name = self._name.text().strip()
        if not name:
            return None
        return SequenceFeature(
            start=self._start.value(),
            end=self._end.value(),
            name=name,
            feature_type=str(self._type.currentData() or self._type.currentText()),
            color=self._color.color_hex(),
            strand=str(self._strand.currentData()),
            source=REPEATVIEWER_INTERNAL_SOURCE,
        )


class ModifyFeatureDialog(QDialog):
    """Edit an existing feature; Type is free-form and Parent preserves GFF3 links."""

    _POSITION_GAP = 8

    def __init__(
        self,
        *,
        feature: SequenceFeature,
        seq_length: int,
        parent_candidates: list[str] | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Modify feature")
        self._track = max(0, int(getattr(feature, "track", 0) or 0))
        raw_source = str(getattr(feature, "source", "") or "").strip()
        # Empty / legacy Sequence Annotation → shared internal source.
        if not raw_source or raw_source == "Sequence Annotation":
            self._source = REPEATVIEWER_INTERNAL_SOURCE
        else:
            self._source = raw_source
        self._pident = float(getattr(feature, "pident", 0.0) or 0.0)
        self._subject_id = str(getattr(feature, "subject_id", "") or "").strip()
        self._coverage = float(getattr(feature, "coverage", 0.0) or 0.0)
        self._source_start = int(getattr(feature, "source_start", 0) or 0)
        self._source_end = int(getattr(feature, "source_end", 0) or 0)

        start_pos, end_pos = sorted((int(feature.start), int(feature.end)))

        self._type = QComboBox(self)
        self._type.setEditable(True)
        self._type.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        known_types = {ftype for _, ftype in load_add_feature_buttons()}
        for label, ftype in load_add_feature_buttons():
            self._type.addItem(label, ftype)
        current_type = str(feature.feature_type or "").strip() or "others"
        index = self._type.findData(current_type)
        if index < 0:
            index = self._type.findText(current_type)
        if index >= 0:
            self._type.setCurrentIndex(index)
        else:
            self._type.setEditText(current_type)
            if current_type not in known_types:
                self._type.addItem(current_type, current_type)
                self._type.setCurrentIndex(self._type.count() - 1)

        self._strand = QComboBox(self)
        self._strand.addItem("Forward (+)", "+")
        self._strand.addItem("Reverse (−)", "-")
        self._strand.setCurrentIndex(0 if feature.is_forward() else 1)

        self._start = QSpinBox(self)
        self._start.setRange(1, seq_length)
        self._start.setValue(start_pos)

        self._end = QSpinBox(self)
        self._end.setRange(1, seq_length)
        self._end.setValue(end_pos)

        self._color = _AddFeatureColorField(feature.color or feature_type_color(current_type), self)

        self._name = QLineEdit(self)
        self._name.setPlaceholderText("e.g. LTR5")
        self._name.setText(feature.name)

        self._parent = QLineEdit(self)
        self._parent.setPlaceholderText("GFF3 Parent ID (e.g. gene name)")
        self._parent.setText(feature.parent)
        if parent_candidates:
            completer = QCompleter(sorted({c for c in parent_candidates if c}), self)
            completer.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
            completer.setFilterMode(Qt.MatchFlag.MatchContains)
            self._parent.setCompleter(completer)

        self._position_row = QWidget(self)
        position_layout = QHBoxLayout(self._position_row)
        position_layout.setContentsMargins(0, 0, 0, 0)
        position_layout.setSpacing(6)
        position_layout.addWidget(self._start)
        self._end_label = QLabel("End:", self._position_row)
        self._end_label.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        position_layout.addWidget(self._end_label)
        position_layout.addWidget(self._end)

        form = QFormLayout()
        form.addRow("Type:", self._type)
        form.addRow("Direction:", self._strand)
        form.addRow("Start:", self._position_row)
        form.addRow("Color:", self._color)
        form.addRow("Name:", self._name)
        form.addRow("Parent:", self._parent)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def _selected_type(self) -> str:
        typed = self._type.currentText().strip()
        if not typed:
            return ""
        for i in range(self._type.count()):
            if self._type.itemText(i).strip() == typed:
                data = self._type.itemData(i)
                if data is not None and str(data).strip():
                    return str(data).strip()
                return typed
        return typed

    def feature(self) -> SequenceFeature | None:
        name = self._name.text().strip()
        if not name:
            return None
        ftype = self._selected_type()
        if not ftype:
            return None
        return SequenceFeature(
            start=self._start.value(),
            end=self._end.value(),
            name=name,
            feature_type=ftype,
            color=self._color.color_hex(),
            strand=str(self._strand.currentData()),
            parent=self._parent.text().strip(),
            track=self._track,
            source=self._source,
            pident=self._pident,
            subject_id=self._subject_id,
            coverage=self._coverage,
            source_start=self._source_start,
            source_end=self._source_end,
        )


class AutoAnnotationSettingsPanel(QWidget):
    """Embeddable LTR / TIR / tandem / ORF auto-annotation parameters."""

    runOptionsChanged = Signal()

    def __init__(
        self,
        settings: AutoAnnotationSettings,
        parent: QWidget | None = None,
        *,
        include_run_controls: bool = True,
    ) -> None:
        super().__init__(parent)
        self._settings = settings
        self._include_run_controls = bool(include_run_controls)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        layout.addWidget(self._section_title("LTR (long terminal repeat)"))
        layout.addLayout(self._terminal_repeat_form(settings.ltr or TerminalRepeatSettings(), "LTR"))
        layout.addWidget(_settings_separator())

        layout.addWidget(self._section_title("TIR (terminal inverted repeat)"))
        layout.addLayout(
            self._terminal_repeat_form(
                settings.tir
                or TerminalRepeatSettings(
                    min_length=DEFAULT_TIR_MIN_LENGTH,
                    max_length=DEFAULT_TIR_MAX_LENGTH,
                    min_pair_distance=DEFAULT_TIR_MIN_PAIR_DISTANCE,
                    similarity_threshold=DEFAULT_TIR_SIMILARITY_THRESHOLD,
                ),
                "TIR",
            )
        )
        layout.addWidget(_settings_separator())

        layout.addWidget(self._section_title("Tandem repeat"))
        layout.addLayout(self._tandem_repeat_form(settings.tandem or TandemRepeatSettings()))
        layout.addWidget(_settings_separator())

        layout.addWidget(self._section_title("ORF (orfipy)"))
        layout.addLayout(self._orf_form(settings.orf or OrfSettings()))
        # PBS / Others are built here but inserted after TSD by the host dialog
        # (see ``pbs_section_widget`` / ``others_section_widget``).
        self._pbs_section = QWidget(self)
        pbs_layout = QVBoxLayout(self._pbs_section)
        pbs_layout.setContentsMargins(0, 0, 0, 0)
        pbs_layout.setSpacing(8)
        pbs_layout.addWidget(_settings_separator())
        pbs_layout.addWidget(self._section_title("PBS (primer binding site)"))
        pbs_layout.addLayout(self._pbs_form(settings.pbs or PbsSettings()))
        self._others_section: QWidget | None
        if self._include_run_controls:
            self._others_section = QWidget(self)
            others_layout = QVBoxLayout(self._others_section)
            others_layout.setContentsMargins(0, 0, 0, 0)
            others_layout.setSpacing(8)
            others_layout.addWidget(_settings_separator())
            others_layout.addWidget(self._section_title("Others"))
            others_layout.addLayout(self._others_form(settings))
        else:
            self._others_section = None
            self._num_threads = None
            self._show_process_window = None
            self._close_process_window = None

        self._ltr_fields = self._terminal_fields[0]
        self._tir_fields = self._terminal_fields[1]
        self._tandem_fields = self._tandem_field_list
        self._orf_fields = self._orf_field_list
        self._pbs_fields = self._pbs_field_list

    def pbs_section_widget(self) -> QWidget:
        """PBS parameters; place after the TSD panel in host layouts."""
        return self._pbs_section

    def others_section_widget(self) -> QWidget | None:
        """Others (threads / process window); place after PBS when present."""
        return self._others_section

    def _section_title(self, text: str) -> QLabel:
        title = QLabel(text)
        font = title.font()
        font.setBold(True)
        title.setFont(font)
        return title

    def _sized_spin(self, spin: QSpinBox | QDoubleSpinBox) -> QSpinBox | QDoubleSpinBox:
        spin.setFixedWidth(_SETTINGS_FIELD_WIDTH)
        spin.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        return spin

    def _sized_line_edit(self, edit: QLineEdit) -> QLineEdit:
        edit.setFixedWidth(_SETTINGS_FIELD_WIDTH)
        edit.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        return edit

    def _field_help(self, field: QWidget, title: str, help_text: str) -> QWidget:
        """``?`` after the value control (BLAST Settings style)."""
        return field_with_help(field, title, help_text, parent=self)

    def _settings_form(self) -> QFormLayout:
        form = QFormLayout()
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        return form

    def _form_label(self, text: str) -> QLabel:
        """Fixed-width label so every section shares the same field-column edge."""
        lab = QLabel(text, self)
        lab.setFixedWidth(_SETTINGS_LABEL_WIDTH)
        lab.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        return lab

    def _terminal_repeat_form(self, cfg: TerminalRepeatSettings, repeat_kind: str) -> QFormLayout:
        form = self._settings_form()
        min_len = QSpinBox(self)
        min_len.setRange(1, 1_000_000)
        min_len.setValue(cfg.min_length)
        max_len = QSpinBox(self)
        max_len.setRange(1, 1_000_000)
        max_len.setValue(cfg.max_length)
        min_distance = QSpinBox(self)
        min_distance.setRange(0, 1_000_000)
        min_distance.setValue(cfg.min_pair_distance)
        max_distance_to_end = QSpinBox(self)
        max_distance_to_end.setRange(0, 1_000_000)
        max_distance_to_end.setValue(cfg.max_distance_to_end)
        max_pairs = QSpinBox(self)
        max_pairs.setRange(1, MAX_TERMINAL_PAIRS_LIMIT)
        max_pairs.setValue(cfg.max_pairs)
        similarity = QDoubleSpinBox(self)
        similarity.setRange(0.0, 1.0)
        similarity.setDecimals(2)
        similarity.setSingleStep(0.05)
        similarity.setValue(cfg.similarity_threshold)
        form.addRow(
            self._form_label("Min length (bp):"),
            self._field_help(
                self._sized_spin(min_len),
                f"{repeat_kind} min length",
                "Minimum repeat length (bp). Values below 11 also set the "
                "self-BLAST word size for this repeat type only (minimum word "
                "size 4; default 11).",
            ),
        )
        form.addRow(
            self._form_label("Max length (bp):"),
            self._field_help(
                self._sized_spin(max_len),
                f"{repeat_kind} max length",
                f"Maximum accepted {repeat_kind} length in bp. Longer candidate "
                "pairs are discarded.",
            ),
        )
        form.addRow(
            self._form_label(f"Min distance between {repeat_kind} (bp):"),
            self._field_help(
                self._sized_spin(min_distance),
                f"Min distance between {repeat_kind}",
                f"Minimum spacer between the two {repeat_kind} copies (bp), "
                f"not the distance from the sequence ends. "
                f"Pairs closer than this are ignored.",
            ),
        )
        form.addRow(
            self._form_label("Max distance to end (bp):"),
            self._field_help(
                self._sized_spin(max_distance_to_end),
                f"{repeat_kind} max distance to end",
                f"The maximum distance between the sequence ends to each of "
                f"the {repeat_kind}. The left (5′) copy must start within this "
                f"distance of the beginning, and the right (3′) copy must end "
                f"within this distance of the sequence end. Larger values "
                f"allow more internal candidates; default 10000.",
            ),
        )
        form.addRow(
            self._form_label("Max pairs:"),
            self._field_help(
                self._sized_spin(max_pairs),
                f"{repeat_kind} max pairs",
                f"How many {repeat_kind} pairs to keep "
                f"(1–{MAX_TERMINAL_PAIRS_LIMIT}). "
                f"1 keeps only the outermost pair; higher values also keep "
                f"inner nested pairs (named {repeat_kind}5_2 / {repeat_kind}3_2, "
                f"etc.). Default {DEFAULT_MAX_TERMINAL_PAIRS}.",
            ),
        )
        form.addRow(
            self._form_label("Similarity threshold:"),
            self._field_help(
                self._sized_spin(similarity),
                f"{repeat_kind} similarity threshold",
                "Minimum fraction identity (0–1) between the two terminal "
                f"repeats for a {repeat_kind} pair to be kept "
                "(e.g. 0.60 = 60%).",
            ),
        )
        if not hasattr(self, "_terminal_fields"):
            self._terminal_fields = []
        self._terminal_fields.append(
            (
                min_len,
                max_len,
                min_distance,
                max_distance_to_end,
                max_pairs,
                similarity,
            )
        )
        return form

    def _tandem_repeat_form(self, cfg: TandemRepeatSettings) -> QFormLayout:
        form = self._settings_form()

        fields: list[QSpinBox] = []
        for motif_size in (2, 3, 4, 5):
            spin = QSpinBox(self)
            spin.setRange(1, 1000)
            spin.setValue(getattr(cfg, f"min_repeat_{motif_size}"))
            form.addRow(
                self._form_label(f"Min repeat number ({motif_size}-motif):"),
                self._field_help(
                    self._sized_spin(spin),
                    f"Min repeat number ({motif_size}-motif)",
                    f"Minimum copy count for exact tandem repeats whose motif "
                    f"is {motif_size} bp long.",
                ),
            )
            fields.append(spin)
        spin_6 = QSpinBox(self)
        spin_6.setRange(1, 1000)
        spin_6.setValue(cfg.min_repeat_6)
        form.addRow(
            self._form_label("Min repeat number (≥ 6-motif):"),
            self._field_help(
                self._sized_spin(spin_6),
                "Min repeat number (≥ 6-motif)",
                "Minimum copy count for exact tandem repeats with motif length "
                "of 6 bp or more.",
            ),
        )
        fields.append(spin_6)
        self._tandem_field_list = fields

        atr_identity = QDoubleSpinBox(self)
        atr_identity.setRange(0.0, 1.0)
        atr_identity.setDecimals(2)
        atr_identity.setSingleStep(0.05)
        atr_identity.setValue(cfg.atr_min_identity)
        atr_min_motif = QSpinBox(self)
        atr_min_motif.setRange(1, 100)
        atr_min_motif.setValue(cfg.atr_min_motif)
        atr_max_errors = QSpinBox(self)
        atr_max_errors.setRange(0, 50)
        atr_max_errors.setValue(cfg.atr_max_errors)

        form.addRow(
            self._form_label("Approx. min identity:"),
            self._field_help(
                self._sized_spin(atr_identity),
                "Approx. min identity",
                "Minimum alignment identity (0–1) while extending an imperfect "
                "repeat. 0.75 = 75% similarity required.",
            ),
        )
        form.addRow(
            self._form_label("Approx. min motif (bp):"),
            self._field_help(
                self._sized_spin(atr_min_motif),
                "Approx. min motif",
                "Minimum motif length (bp) for approximate repeat seeds.",
            ),
        )
        form.addRow(
            self._form_label("Approx. max errors:"),
            self._field_help(
                self._sized_spin(atr_max_errors),
                "Approx. max errors",
                "Maximum consecutive mismatches allowed while extending an "
                "imperfect repeat.",
            ),
        )
        self._tandem_approximate_fields = (atr_identity, atr_min_motif, atr_max_errors)

        approximate_cb = QCheckBox("Enable approximate tandem repeat", self)
        approximate_cb.setChecked(cfg.enable_approximate)
        # Keep the enable tick after its Approx.* parameters.
        form.addRow(
            self._form_label(""),
            field_with_help(
                approximate_cb,
                "Approximate tandem repeat",
                "Run pytrf ATRFinder to detect imperfect repeats with "
                "mismatches. ATR runs before exact STR/GTR so long approximate "
                "blocks are kept.",
                parent=self,
                stretch=True,
            ),
        )
        self._tandem_approximate_cb = approximate_cb
        approximate_cb.toggled.connect(self._sync_tandem_approximate_fields)
        self._sync_tandem_approximate_fields(cfg.enable_approximate)

        return form

    def _sync_tandem_approximate_fields(self, enabled: bool) -> None:
        for widget in getattr(self, "_tandem_approximate_fields", ()):
            widget.setEnabled(enabled)

    def _orf_form(self, cfg: OrfSettings) -> QFormLayout:
        from repeatviewer.backends.auto_annotation_settings import (
            DEFAULT_ORF_FRAMES,
            normalize_orf_frames,
        )

        form = self._settings_form()
        min_len = QSpinBox(self)
        min_len.setRange(3, 1_000_000)
        min_len.setValue(cfg.min_length)
        max_len = QSpinBox(self)
        max_len.setRange(3, 1_000_000)
        max_len.setValue(cfg.max_length)
        alt_starts = QLineEdit(self)
        alt_starts.setPlaceholderText("e.g. TTG, CTG")
        alt_starts.setText(cfg.alternative_start_codons)
        form.addRow(
            self._form_label("Min ORF length (bp):"),
            self._field_help(
                self._sized_spin(min_len),
                "Min ORF length",
                "Minimum open reading frame length in bp (must be a multiple "
                "of 3 after translation rules).",
            ),
        )
        form.addRow(
            self._form_label("Max ORF length (bp):"),
            self._field_help(
                self._sized_spin(max_len),
                "Max ORF length",
                "Maximum open reading frame length in bp.",
            ),
        )
        form.addRow(
            self._form_label("Add other start codons:"),
            self._field_help(
                self._sized_line_edit(alt_starts),
                "Other start codons",
                "Optional extra start codons (comma- or space-separated), "
                "e.g. <b>TTG, CTG</b>, in addition to ATG.",
            ),
        )

        selected = set(normalize_orf_frames(cfg.frames))
        # Compact checkbox grid packed left inside a field-width host so the
        # ``?`` lines up with start-codon / spin rows (without stretching gaps).
        frames_inner = QWidget(self)
        frames_inner.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        frames_grid = QGridLayout(frames_inner)
        frames_grid.setContentsMargins(0, 0, 0, 0)
        frames_grid.setHorizontalSpacing(4)
        frames_grid.setVerticalSpacing(2)
        frame_checks: dict[int, QCheckBox] = {}
        for col, frame in enumerate((1, 2, 3)):
            frames_grid.setColumnStretch(col, 0)
            fwd = QCheckBox(f"Frame {frame}", frames_inner)
            fwd.setToolTip(f"Forward reading frame {frame}")
            fwd.setChecked(frame in selected)
            frame_checks[frame] = fwd
            frames_grid.addWidget(fwd, 0, col, Qt.AlignmentFlag.AlignLeft)
            rev = QCheckBox(f"Reverse {frame}", frames_inner)
            rev.setToolTip(f"Reverse reading frame {frame}")
            rev.setChecked((-frame) in selected)
            frame_checks[-frame] = rev
            frames_grid.addWidget(rev, 1, col, Qt.AlignmentFlag.AlignLeft)
        frames_host = QWidget(self)
        frames_host.setFixedWidth(_SETTINGS_FIELD_WIDTH)
        frames_host.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        frames_host_layout = QHBoxLayout(frames_host)
        frames_host_layout.setContentsMargins(0, 0, 0, 0)
        frames_host_layout.setSpacing(0)
        frames_host_layout.addWidget(
            frames_inner,
            0,
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop,
        )
        frames_host_layout.addStretch(1)
        form.addRow(
            self._form_label("Reading frames:"),
            self._field_help(
                frames_host,
                "Reading frames",
                "Keep ORFs from these orfipy reading frames "
                "(Frame 1–3 forward, Reverse 1–3). "
                "Default is all six frames.",
            ),
        )
        self._orf_frame_checks = frame_checks
        self._orf_default_frames = DEFAULT_ORF_FRAMES
        self._orf_field_list = (min_len, max_len, alt_starts)
        return form

    def _pbs_form(self, cfg: PbsSettings) -> QFormLayout:
        form = self._settings_form()
        window = QSpinBox(self)
        window.setRange(1, 10_000)
        window.setValue(cfg.window)
        min_len = QSpinBox(self)
        min_len.setRange(1, 500)
        min_len.setValue(cfg.min_length)
        min_ident = QDoubleSpinBox(self)
        min_ident.setRange(0.0, 100.0)
        min_ident.setDecimals(1)
        min_ident.setSingleStep(1.0)
        min_ident.setValue(cfg.min_identity)
        min_cov = QDoubleSpinBox(self)
        min_cov.setRange(0.0, 1.0)
        min_cov.setDecimals(2)
        min_cov.setSingleStep(0.05)
        min_cov.setValue(cfg.min_coverage)
        max_dist = QSpinBox(self)
        max_dist.setRange(0, 10_000)
        max_dist.setValue(cfg.max_distance_from_ltr5_nt)
        orient = NoWheelComboBox(self)
        orient.addItem("Auto (both strands)", "auto")
        orient.addItem("Plus (+)", "+")
        orient.addItem("Minus (−)", "-")
        orient_idx = orient.findData(cfg.orientation or DEFAULT_PBS_ORIENTATION)
        orient.setCurrentIndex(orient_idx if orient_idx >= 0 else 0)
        motif = PathBrowseRow(
            self,
            mode="file",
            file_filter="FASTA (*.fa *.fasta *.fna);;All files (*)",
            browse_title="Choose PBS motif FASTA",
            placeholder="Bundled Neumann 2019 motifs (default)",
            minimum_width=_SETTINGS_FIELD_WIDTH,
        )
        motif.set_path(cfg.motif_database)
        if not cfg.motif_database:
            bundled = default_motif_fasta()
            if bundled.is_file():
                motif.setToolTip(f"Default: {bundled}")

        form.addRow(
            self._form_label("Search window (bp):"),
            self._field_help(
                self._sized_spin(window),
                "PBS search window",
                "Bases immediately inside the 5′ LTR to search for a PBS "
                f"(default {DEFAULT_PBS_WINDOW}). "
                "Max distance from LTR5 further limits accepted hits.",
            ),
        )
        form.addRow(
            self._form_label("Min alignment length:"),
            self._field_help(
                self._sized_spin(min_len),
                "Min PBS alignment length",
                "Minimum ungapped blastn alignment length "
                f"(default {DEFAULT_PBS_MIN_LENGTH}).",
            ),
        )
        form.addRow(
            self._form_label("Min identity (%):"),
            self._field_help(
                self._sized_spin(min_ident),
                "Min PBS identity",
                f"Minimum percent identity (default {DEFAULT_PBS_MIN_IDENTITY:g}).",
            ),
        )
        form.addRow(
            self._form_label("Min motif coverage:"),
            self._field_help(
                self._sized_spin(min_cov),
                "Min PBS motif coverage",
                "Fraction of the motif that must align "
                f"(default {DEFAULT_PBS_MIN_COVERAGE:g}).",
            ),
        )
        form.addRow(
            self._form_label("Max distance from LTR5:"),
            self._field_help(
                self._sized_spin(max_dist),
                "Max PBS distance from LTR5",
                "Maximum gap (nt) between the biological 5′ LTR and the PBS "
                f"(default {DEFAULT_PBS_MAX_DISTANCE}).",
            ),
        )
        form.addRow(
            self._form_label("Orientation:"),
            self._field_help(
                orient,
                "PBS orientation",
                "Auto tries both sides of the LTR pair; + searches after the "
                "left LTR; − searches before the right LTR.",
            ),
        )
        form.addRow(
            self._form_label("Motif FASTA:"),
            self._field_help(
                motif,
                "PBS motif database",
                "Optional custom FASTA of PBS motifs. Leave empty to use the "
                "bundled Neumann et al. 2019 unique motifs.",
            ),
        )
        self._pbs_orient = orient
        self._pbs_motif = motif
        self._pbs_field_list = (window, min_len, min_ident, min_cov, max_dist)
        return form

    def _others_form(self, settings: AutoAnnotationSettings) -> QFormLayout:
        form = self._settings_form()
        threads = QSpinBox(self)
        threads.setRange(1, 256)
        threads.setValue(max(1, int(settings.num_threads)))
        form.addRow(
            self._form_label("Threads:"),
            self._field_help(
                self._sized_spin(threads),
                "Threads",
                "CPU threads inside each sequence job (like BLAST Threads). "
                "How many sequences run at once is Preferences → "
                "Max. parallel jobs. Peak use ≈ parallel jobs × Threads. "
                "Default is 4. Shared with Similarity annotation.",
            ),
        )
        self._num_threads = threads
        threads.valueChanged.connect(lambda _v: self.runOptionsChanged.emit())
        show_box = QCheckBox("Show process window while run", self)
        show_box.setChecked(bool(settings.show_process_window))
        show_box.setToolTip(
            "Open a window that streams Feature detection progress and tool "
            "messages. Useful for long runs."
        )
        close_box = QCheckBox("Close process window when finished", self)
        close_box.setChecked(bool(settings.close_process_window_when_finished))
        close_box.setEnabled(show_box.isChecked())
        close_box.setToolTip(
            "Automatically close the process window after a successful run. "
            "On failure the window stays open so you can read the log."
        )
        show_box.toggled.connect(close_box.setEnabled)
        show_box.toggled.connect(lambda _v: self.runOptionsChanged.emit())
        close_box.toggled.connect(lambda _v: self.runOptionsChanged.emit())
        form.addRow(
            self._form_label(""),
            field_with_help(
                show_box,
                "Show process window",
                "Open a live process window while Feature detection runs. "
                "Shared with Similarity annotation.",
                parent=self,
                stretch=True,
            ),
        )
        form.addRow(
            self._form_label(""),
            field_with_help(
                close_box,
                "Close process window when finished",
                "Close the process window automatically after a successful run. "
                "Shared with Similarity annotation.",
                parent=self,
                stretch=True,
            ),
        )
        self._show_process_window = show_box
        self._close_process_window = close_box
        return form

    def _selected_orf_frames(self) -> tuple[int, ...]:
        from repeatviewer.backends.auto_annotation_settings import normalize_orf_frames

        checks = getattr(self, "_orf_frame_checks", {})
        chosen = tuple(
            frame for frame, check in checks.items() if check.isChecked()
        )
        return normalize_orf_frames(chosen)

    def settings(self) -> AutoAnnotationSettings:
        ltr = self._read_terminal_settings(self._ltr_fields, default_min=100)
        tir = self._read_terminal_settings(self._tir_fields, default_min=20)
        tandem = TandemRepeatSettings(
            min_repeat_2=self._tandem_fields[0].value(),
            min_repeat_3=self._tandem_fields[1].value(),
            min_repeat_4=self._tandem_fields[2].value(),
            min_repeat_5=self._tandem_fields[3].value(),
            min_repeat_6=self._tandem_fields[4].value(),
            enable_approximate=self._tandem_approximate_cb.isChecked(),
            atr_min_identity=self._tandem_approximate_fields[0].value(),
            atr_min_motif=self._tandem_approximate_fields[1].value(),
            atr_max_errors=self._tandem_approximate_fields[2].value(),
            atr_max_motif=TandemRepeatSettings().atr_max_motif,
            atr_max_extend=TandemRepeatSettings().atr_max_extend,
            atr_min_seed_repeat=TandemRepeatSettings().atr_min_seed_repeat,
            atr_min_seed_length=TandemRepeatSettings().atr_min_seed_length,
        )
        min_len, max_len, alt_starts = self._orf_fields
        orf = OrfSettings(
            min_length=min_len.value(),
            max_length=max(max_len.value(), min_len.value()),
            alternative_start_codons=alt_starts.text().strip(),
            frames=self._selected_orf_frames(),
        )
        window, pbs_min_len, pbs_ident, pbs_cov, pbs_max_dist = self._pbs_fields
        orient_data = self._pbs_orient.currentData()
        pbs = PbsSettings(
            window=int(window.value()),
            min_length=int(pbs_min_len.value()),
            min_identity=float(pbs_ident.value()),
            min_coverage=float(pbs_cov.value()),
            max_distance_from_ltr5_nt=int(pbs_max_dist.value()),
            orientation=str(orient_data or DEFAULT_PBS_ORIENTATION),
            motif_database=str(self._pbs_motif.path() or "").strip(),
        )
        if self._include_run_controls and self._num_threads is not None:
            threads = max(1, int(self._num_threads.value()))
            show = bool(self._show_process_window.isChecked())
            close = bool(self._close_process_window.isChecked())
        else:
            from repeatviewer.backends.auto_annotation_settings import (
                load_annotation_run_options,
            )

            threads, show, close = load_annotation_run_options()
        return AutoAnnotationSettings(
            ltr=ltr,
            tir=tir,
            tandem=tandem,
            orf=orf,
            pbs=pbs,
            num_threads=threads,
            show_process_window=show,
            close_process_window_when_finished=close,
        )

    def set_run_options(
        self,
        *,
        num_threads: int,
        show_process_window: bool,
        close_process_window_when_finished: bool,
    ) -> None:
        """Update shared Threads / process-window controls (no signal echo)."""
        if not self._include_run_controls or self._num_threads is None:
            return
        self._num_threads.blockSignals(True)
        self._show_process_window.blockSignals(True)
        self._close_process_window.blockSignals(True)
        self._num_threads.setValue(max(1, int(num_threads)))
        self._show_process_window.setChecked(bool(show_process_window))
        self._close_process_window.setChecked(bool(close_process_window_when_finished))
        self._close_process_window.setEnabled(self._show_process_window.isChecked())
        self._num_threads.blockSignals(False)
        self._show_process_window.blockSignals(False)
        self._close_process_window.blockSignals(False)
    @staticmethod
    def _read_terminal_settings(
        fields: tuple[
            QSpinBox, QSpinBox, QSpinBox, QSpinBox, QSpinBox, QDoubleSpinBox
        ],
        *,
        default_min: int,
    ) -> TerminalRepeatSettings:
        (
            min_len,
            max_len,
            min_distance,
            max_distance_to_end,
            max_pairs,
            similarity,
        ) = fields
        return TerminalRepeatSettings(
            min_length=min_len.value() or default_min,
            max_length=max(max_len.value(), min_len.value()),
            min_pair_distance=min_distance.value(),
            max_distance_to_end=max_distance_to_end.value(),
            max_pairs=max_pairs.value(),
            similarity_threshold=similarity.value(),
        )


class TsdSearchSettingsPanel(QWidget):
    """TSD boundaries, search radius, length / distance."""

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        initial: TsdSearchOptions | None = None,
        note: str | None = None,
        allow_user_boundary: bool = True,
    ) -> None:
        super().__init__(parent)
        del note  # gray help blurb removed; kept for call-site compatibility
        opts = initial or load_tsd_search_options()
        self._allow_user_boundary = bool(allow_user_boundary)

        title = QLabel("TSD (target site duplication)", self)
        font = title.font()
        font.setBold(True)
        title.setFont(font)

        form = QFormLayout()
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)

        def _label(text: str) -> QLabel:
            lab = QLabel(text, self)
            lab.setFixedWidth(_SETTINGS_LABEL_WIDTH)
            lab.setAlignment(
                Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
            )
            return lab

        def _sized(spin: QSpinBox) -> QSpinBox:
            spin.setFixedWidth(_SETTINGS_FIELD_WIDTH)
            spin.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
            return spin

        self._boundary_source = NoWheelComboBox(self)
        self._boundary_source.addItem("RepeatViewer", TSD_BOUNDARY_REPEATVIEWER)
        if self._allow_user_boundary:
            self._boundary_source.addItem("User define", TSD_BOUNDARY_USER)
        source = (
            opts.boundary_source
            if opts.boundary_source in TSD_BOUNDARY_SOURCES
            else TSD_BOUNDARY_REPEATVIEWER
        )
        if not self._allow_user_boundary:
            source = TSD_BOUNDARY_REPEATVIEWER
        idx = self._boundary_source.findData(source)
        self._boundary_source.setCurrentIndex(max(0, idx))
        self._boundary_source.setFixedWidth(max(_SETTINGS_FIELD_WIDTH, 160))
        help_boundary = (
            "RepeatViewer: *.boundaries_RV.json."
            if not self._allow_user_boundary
            else (
                "RepeatViewer: *.boundaries_RV.json. "
                "User define: popup at run time to drag the MSA."
            )
        )
        form.addRow(
            _label("Boundary definition:"),
            field_with_help(
                self._boundary_source,
                "Boundary definition",
                help_boundary,
                parent=self,
            ),
        )

        self._tsd_radius = QSpinBox(self)
        self._tsd_radius.setRange(1, 10_000_000)
        self._tsd_radius.setValue(max(1, int(opts.nt_radius)))
        form.addRow(
            _label("Search range around boundary (bp):"),
            field_with_help(
                _sized(self._tsd_radius),
                "Search range around boundary",
                "MSA columns to search on both sides of each TE edge "
                "(left and right of the left boundary, and left and right "
                "of the right boundary; gaps count as positions).",
                parent=self,
            ),
        )

        self._tsd_min_len = QSpinBox(self)
        self._tsd_min_len.setRange(2, 200)
        self._tsd_min_len.setValue(opts.min_tsd_length)
        form.addRow(
            _label("Min TSD length (bp):"),
            field_with_help(
                _sized(self._tsd_min_len),
                "Min TSD length",
                "Minimum duplicated motif length.",
                parent=self,
            ),
        )

        self._tsd_max_len = QSpinBox(self)
        self._tsd_max_len.setRange(2, 200)
        self._tsd_max_len.setValue(opts.max_tsd_length)
        form.addRow(
            _label("Max TSD length (bp):"),
            field_with_help(
                _sized(self._tsd_max_len),
                "Max TSD length",
                "Maximum duplicated motif length.",
                parent=self,
            ),
        )

        self._tsd_min_dist = QSpinBox(self)
        self._tsd_min_dist.setRange(0, 10_000_000)
        self._tsd_min_dist.setValue(opts.min_distance)
        form.addRow(
            _label("Min distance between TSDs (bp):"),
            field_with_help(
                _sized(self._tsd_min_dist),
                "Min distance between TSDs",
                "Minimum TE span between the left and right TSD copies. "
                "Default 200 bp. 0 = no minimum.",
                parent=self,
            ),
        )

        self._write_text_report = QCheckBox(
            "Write .tsd_search.txt report", self
        )
        self._write_text_report.setChecked(bool(opts.write_text_report))
        form.addRow(
            _label(""),
            field_with_help(
                self._write_text_report,
                "Write .tsd_search.txt report",
                "When checked, write a clustered TSD text report beside the "
                "MSA (shows in Files). Off by default.",
                parent=self,
            ),
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.addWidget(_settings_separator())
        layout.addWidget(title)
        layout.addLayout(form)

    def options(self) -> TsdSearchOptions:
        min_len = int(self._tsd_min_len.value())
        max_len = max(min_len, int(self._tsd_max_len.value()))
        prev = load_tsd_search_options()
        source = self._boundary_source.currentData()
        if not self._allow_user_boundary:
            source = TSD_BOUNDARY_REPEATVIEWER
        elif source not in TSD_BOUNDARY_SOURCES:
            source = TSD_BOUNDARY_REPEATVIEWER
        return TsdSearchOptions(
            left_boundary=prev.left_boundary,
            right_boundary=prev.right_boundary,
            nt_radius=int(self._tsd_radius.value()),
            min_tsd_length=min_len,
            max_tsd_length=max_len,
            min_distance=self._tsd_min_dist.value(),
            max_distance=prev.max_distance,
            max_sim_threshold=prev.max_sim_threshold,
            nt_threshold=prev.nt_threshold,
            cluster_eps=prev.cluster_eps,
            boundary_source=str(source),
            write_text_report=self._write_text_report.isChecked(),
        )


class AutoAnnotationSettingsDialog(QDialog):
    """Configure LTR, TIR, tandem-repeat, ORF, PBS, and TSD parameters."""

    def __init__(
        self,
        settings: AutoAnnotationSettings,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Feature detection settings")
        self.setMinimumWidth(_SETTINGS_DIALOG_MIN_WIDTH)
        self._panel = AutoAnnotationSettingsPanel(settings, self)
        self._tsd_panel = TsdSearchSettingsPanel(self)
        body = QWidget(self)
        body_layout = QVBoxLayout(body)
        body_layout.setContentsMargins(0, 0, 0, 0)
        body_layout.setSpacing(12)
        body_layout.addWidget(self._panel)
        body_layout.addWidget(self._tsd_panel)
        body_layout.addWidget(self._panel.pbs_section_widget())
        others = self._panel.others_section_widget()
        if others is not None:
            body_layout.addWidget(others)
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(body)
        # Cap height so the dense LTR/TIR/Tandem/ORF/TSD form scrolls on small screens.
        screen = QApplication.primaryScreen()
        max_h = 640
        if screen is not None:
            max_h = max(420, min(720, int(screen.availableGeometry().height() * 0.75)))
        scroll.setMaximumHeight(max_h)

        layout = QVBoxLayout(self)
        layout.setSpacing(8)
        layout.addWidget(scroll, stretch=1)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def settings(self) -> AutoAnnotationSettings:
        return self._panel.settings()

    def tsd_options(self) -> TsdSearchOptions:
        return self._tsd_panel.options()


class SimilarityAnnotationSettingsPanel(QWidget):
    """Embeddable similarity-annotation search thresholds."""

    runOptionsChanged = Signal()

    def __init__(
        self,
        settings: SimilarityAnnotationSettings,
        parent: QWidget | None = None,
        *,
        include_run_controls: bool = True,
    ) -> None:
        super().__init__(parent)
        self._include_run_controls = bool(include_run_controls)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        form = QFormLayout()
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)

        self._blastn_evalue = QDoubleSpinBox(self)
        configure_report_evalue_spinbox(self._blastn_evalue)
        self._blastn_evalue.setValue(settings.blastn_evalue)
        self._blastn_evalue.setFixedWidth(_SETTINGS_FIELD_WIDTH)
        form.addRow(
            "BLASTN E-value:",
            field_with_help(
                self._blastn_evalue,
                "BLASTN E-value",
                "Reporting E-value for BLASTN (DNA) similarity searches.",
                parent=self,
            ),
        )

        self._blastp_evalue = QDoubleSpinBox(self)
        configure_report_evalue_spinbox(self._blastp_evalue)
        self._blastp_evalue.setValue(settings.blastp_evalue)
        self._blastp_evalue.setFixedWidth(_SETTINGS_FIELD_WIDTH)
        form.addRow(
            "BLASTP E-value:",
            field_with_help(
                self._blastp_evalue,
                "BLASTP E-value",
                "Reporting E-value for BLASTP (protein) similarity searches.",
                parent=self,
            ),
        )

        self._cdd_evalue = QDoubleSpinBox(self)
        configure_report_evalue_spinbox(self._cdd_evalue)
        self._cdd_evalue.setValue(settings.cdd_evalue)
        self._cdd_evalue.setFixedWidth(_SETTINGS_FIELD_WIDTH)
        form.addRow(
            "CDD E-value (rpstblastn):",
            field_with_help(
                self._cdd_evalue,
                "CDD E-value",
                "E-value cutoff for CDD rpstblastn similarity searches.",
                parent=self,
            ),
        )

        self._hmm_evalue = QDoubleSpinBox(self)
        configure_report_evalue_spinbox(self._hmm_evalue)
        self._hmm_evalue.setValue(settings.hmm_evalue)
        self._hmm_evalue.setFixedWidth(_SETTINGS_FIELD_WIDTH)
        form.addRow(
            "HMM E-value (pyhmmer):",
            field_with_help(
                self._hmm_evalue,
                "HMM E-value",
                "Reporting E-value for pyhmmer HMM similarity searches.",
                parent=self,
            ),
        )

        self._max_hits = QSpinBox(self)
        self._max_hits.setRange(1, 100)
        self._max_hits.setValue(settings.max_hits_per_database)
        self._max_hits.setFixedWidth(_SETTINGS_FIELD_WIDTH)
        form.addRow(
            "Maximum hits per database:",
            field_with_help(
                self._max_hits,
                "Maximum hits per database",
                (
                    f"Caps how many hits from <b>one database</b> are kept on the "
                    f"track after search (default 10). Protein and DNA searches "
                    f"use the same steps.<br><br>"
                    f"<b>How filtering works</b> (in order):<br>"
                    f"1. <b>Same span</b> — If two hits cover exactly the same "
                    f"coordinates, keep the one with the higher bit score.<br>"
                    f"2. <b>Fully nested</b> — If a shorter hit lies completely "
                    f"inside a longer one and its E-value is not better "
                    f"(≥ covering hit, including both E=0), drop the shorter "
                    f"hit.<br>"
                    f"3. <b>Partial overlap</b> (only if still over the limit) — "
                    f"Hits that overlap only partly are grouped; keep the "
                    f"highest bit-score hit in each group.<br>"
                    f"4. <b>Top N by score</b> (only if still over the limit) — "
                    f"Keep the best-scoring hits until this maximum is reached."
                    f"<br><br>"
                    f"<b>Examples</b> (max = 2):<br>"
                    f"• Nested: long hit 100–500 (E=1e−40) and short hit 200–250 "
                    f"(E=1e−10) → short is dropped (worse E and fully inside).<br>"
                    f"• Nested exception: long 100–500 (E=1e−5) and short 200–250 "
                    f"(E=1e−40) → both kept (short has a better E-value).<br>"
                    f"• Partial overlap: A 100–200 (score 80), B 150–250 "
                    f"(score 120), C 400–500 (score 50) and still above the "
                    f"limit → A and B form one group (keep B); C is separate → "
                    f"kept hits are B and C, then capped to the top 2 by score."
                ),
                parent=self,
            ),
        )

        if self._include_run_controls:
            self._num_threads = QSpinBox(self)
            self._num_threads.setRange(1, 256)
            self._num_threads.setValue(max(1, int(settings.num_threads)))
            self._num_threads.setFixedWidth(_SETTINGS_FIELD_WIDTH)
            form.addRow(
                "Threads:",
                field_with_help(
                    self._num_threads,
                    "Threads",
                    "CPU threads for BLAST / HMM / CDD inside each sequence job "
                    "(like BLAST Threads). How many sequences run at once is "
                    "Preferences → Max. parallel jobs. When several databases "
                    "run together on one sequence, this budget is split across "
                    "them. Peak use ≈ parallel jobs × Threads. Default is 4. "
                    "Shared with Feature detection.",
                    parent=self,
                ),
            )
            self._num_threads.valueChanged.connect(
                lambda _v: self.runOptionsChanged.emit()
            )
        else:
            self._num_threads = None

        # Tick aligns with field column (left edge of number boxes); text to the right.
        self._show_db_name = QCheckBox("Show database name on feature track", self)
        self._show_db_name.setChecked(settings.show_database_name_on_track)
        form.addRow(
            "",
            field_with_help(
                self._show_db_name,
                "Show database name",
                "When checked, similarity hits are labeled as "
                "“Database: target”. When unchecked, only the target or "
                "domain name is shown.",
                parent=self,
                stretch=True,
            ),
        )

        if self._include_run_controls:
            self._show_process_window = QCheckBox(
                "Show process window while run", self
            )
            self._show_process_window.setChecked(bool(settings.show_process_window))
            self._show_process_window.setToolTip(
                "Open a window that streams Similarity annotation progress and "
                "tool messages. Useful for long runs."
            )
            form.addRow(
                "",
                field_with_help(
                    self._show_process_window,
                    "Show process window",
                    "Open a live process window while Similarity annotation runs. "
                    "Shared with Feature detection.",
                    parent=self,
                    stretch=True,
                ),
            )

            self._close_process_window = QCheckBox(
                "Close process window when finished", self
            )
            self._close_process_window.setChecked(
                bool(settings.close_process_window_when_finished)
            )
            self._close_process_window.setEnabled(
                self._show_process_window.isChecked()
            )
            self._close_process_window.setToolTip(
                "Automatically close the process window after a successful run. "
                "On failure the window stays open so you can read the log."
            )
            self._show_process_window.toggled.connect(
                self._close_process_window.setEnabled
            )
            self._show_process_window.toggled.connect(
                lambda _v: self.runOptionsChanged.emit()
            )
            self._close_process_window.toggled.connect(
                lambda _v: self.runOptionsChanged.emit()
            )
            form.addRow(
                "",
                field_with_help(
                    self._close_process_window,
                    "Close process window when finished",
                    "Close the process window automatically after a successful run. "
                    "Shared with Feature detection.",
                    parent=self,
                    stretch=True,
                ),
            )
        else:
            self._show_process_window = None
            self._close_process_window = None
        layout.addLayout(form)

    def settings(self) -> SimilarityAnnotationSettings:
        if self._include_run_controls and self._num_threads is not None:
            threads = max(1, int(self._num_threads.value()))
            show = bool(self._show_process_window.isChecked())
            close = bool(self._close_process_window.isChecked())
        else:
            from repeatviewer.backends.auto_annotation_settings import (
                load_annotation_run_options,
            )

            threads, show, close = load_annotation_run_options()
        return SimilarityAnnotationSettings(
            blastn_evalue=self._blastn_evalue.value(),
            blastp_evalue=self._blastp_evalue.value(),
            hmm_evalue=self._hmm_evalue.value(),
            cdd_evalue=self._cdd_evalue.value(),
            max_hits_per_database=self._max_hits.value(),
            num_threads=threads,
            show_database_name_on_track=self._show_db_name.isChecked(),
            show_process_window=show,
            close_process_window_when_finished=close,
        )

    def set_run_options(
        self,
        *,
        num_threads: int,
        show_process_window: bool,
        close_process_window_when_finished: bool,
    ) -> None:
        """Update shared Threads / process-window controls (no signal echo)."""
        if not self._include_run_controls or self._num_threads is None:
            return
        self._num_threads.blockSignals(True)
        self._show_process_window.blockSignals(True)
        self._close_process_window.blockSignals(True)
        self._num_threads.setValue(max(1, int(num_threads)))
        self._show_process_window.setChecked(bool(show_process_window))
        self._close_process_window.setChecked(bool(close_process_window_when_finished))
        self._close_process_window.setEnabled(self._show_process_window.isChecked())
        self._num_threads.blockSignals(False)
        self._show_process_window.blockSignals(False)
        self._close_process_window.blockSignals(False)




class TransferAnnotationSettingsDialog(QDialog):
    """Similarity mode, threads, and process-window options for Transfer."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        from repeatviewer.backends.annotation_transfer import (
            TRANSFER_STRATEGY_CHOICES,
            load_transfer_strategy,
            mafft_liftover_max_kb,
        )
        from repeatviewer.backends.auto_annotation_settings import (
            load_annotation_run_options,
        )
        from repeatviewer.widgets.combo_boxes import NoWheelComboBox

        self.setWindowTitle("Transfer annotation settings")
        self.setMinimumWidth(420)
        threads, show, close = load_annotation_run_options()

        layout = QVBoxLayout(self)
        layout.setSpacing(8)
        form = QFormLayout()
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)

        mafft_kb = mafft_liftover_max_kb()
        self._strategy = NoWheelComboBox(self)
        for value, label in TRANSFER_STRATEGY_CHOICES:
            self._strategy.addItem(label, value)
        strategy_idx = self._strategy.findData(load_transfer_strategy())
        if strategy_idx >= 0:
            self._strategy.setCurrentIndex(strategy_idx)
        self._strategy.setFixedWidth(_SETTINGS_FIELD_WIDTH)
        form.addRow(
            "Target similar to source:",
            field_with_help(
                self._strategy,
                "Target similar to source",
                "<b>Target similar to source</b><br>"
                "Is the open target a related copy of the annotated source?<br>"
                "• <b>Yes</b> — whole-sequence align + lift features "
                f"(MAFFT if &lt;{mafft_kb} kb, else minimap2), then blastn "
                "for any leftovers.<br>"
                "• <b>No</b> — skip whole-seq lift; write all feature DNAs "
                "into one query FASTA and run a single blastn against the "
                "target (divergent, truncated, or rearranged copies).<br><br>"
                "Identity / Coverage sliders on the Transfer rail still "
                "filter the preview after Run.",
                parent=self,
            ),
        )

        self._num_threads = QSpinBox(self)
        self._num_threads.setRange(1, 256)
        self._num_threads.setValue(max(1, int(threads)))
        self._num_threads.setFixedWidth(_SETTINGS_FIELD_WIDTH)
        form.addRow(
            "Threads:",
            field_with_help(
                self._num_threads,
                "Threads",
                "<b>What this controls</b><br>"
                "CPU threads passed into each Transfer alignment / search "
                "tool for the current source→target pair:<br>"
                "• <b>MAFFT</b> (<code>--thread</code>) — whole-sequence "
                f"pairwise lift when both sequences are under {mafft_kb} kb "
                "(fast FFT-NS-1 mode)<br>"
                "• <b>minimap2</b> (<code>-t</code>) — whole-sequence lift "
                f"when either sequence is ≥{mafft_kb} kb, or when MAFFT is "
                "missing / too weak<br>"
                "• <b>blastn</b> (<code>-num_threads</code>) — one search of "
                "all feature DNAs (Yes leftovers, or No / blastn-only "
                "mode)<br><br>"
                "<b>How Transfer uses them</b><br>"
                "Each source→target pair gets this Threads budget. When the "
                "source is a multi-record FASTA or a folder, several pairs "
                "can run at once — concurrency comes from Preferences → "
                "<b>Max. parallel jobs</b> (capped by how many annotated "
                "sources you have and by CPUs).<br><br>"
                "Peak use ≈ parallel jobs × Threads.<br><br>"
                "<b>Shared setting</b><br>"
                "The same Threads number is stored for Feature detection, "
                "Similarity annotation, Sequence Annotation, and Transfer. "
                "Changing it here updates those panels too. Default is 4.",
                parent=self,
            ),
        )

        self._show_process_window = QCheckBox(
            "Show process window while run", self
        )
        self._show_process_window.setChecked(bool(show))
        form.addRow(
            "",
            field_with_help(
                self._show_process_window,
                "Show process window",
                "Open a live process window with Transfer steps and command "
                "lines (MAFFT / minimap2 / blastn). Shared with Feature "
                "detection and Similarity annotation.",
                parent=self,
                stretch=True,
            ),
        )

        self._close_process_window = QCheckBox(
            "Close process window when finished", self
        )
        self._close_process_window.setChecked(bool(close))
        self._close_process_window.setEnabled(
            self._show_process_window.isChecked()
        )
        self._show_process_window.toggled.connect(
            self._close_process_window.setEnabled
        )
        form.addRow(
            "",
            field_with_help(
                self._close_process_window,
                "Close process window when finished",
                "Close the process window automatically after a successful run. "
                "On failure the window stays open so you can read the log.",
                parent=self,
                stretch=True,
            ),
        )
        layout.addLayout(form)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok
            | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def run_options(self) -> tuple[int, bool, bool]:
        """Return ``(num_threads, show_process_window, close_when_finished)``."""
        return (
            max(1, int(self._num_threads.value())),
            bool(self._show_process_window.isChecked()),
            bool(self._close_process_window.isChecked()),
        )

    def transfer_strategy(self) -> str:
        """Selected Transfer strategy (`hybrid` / `liftover` / `blastn`)."""
        from repeatviewer.backends.annotation_transfer import (
            normalize_transfer_strategy,
        )

        return normalize_transfer_strategy(self._strategy.currentData())


class SimilarityAnnotationSettingsDialog(QDialog):
    """Configure similarity-annotation search thresholds."""

    manageDatabaseClicked = Signal()

    def __init__(
        self,
        settings: SimilarityAnnotationSettings,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Similarity annotation settings")
        self.setMinimumWidth(_SETTINGS_DIALOG_MIN_WIDTH)
        self._panel = SimilarityAnnotationSettingsPanel(settings, self)
        layout = QVBoxLayout(self)
        layout.setSpacing(8)
        layout.addWidget(self._panel)
        footer = QHBoxLayout()
        footer.setContentsMargins(0, 0, 0, 0)
        footer.setSpacing(8)
        manage_btn = QPushButton("Manage database", self)
        manage_btn.setToolTip(
            "Open Manage Database to add or edit similarity databases"
        )
        manage_btn.clicked.connect(self.manageDatabaseClicked.emit)
        footer.addWidget(manage_btn, 0, Qt.AlignmentFlag.AlignLeft)
        footer.addStretch(1)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        footer.addWidget(buttons, 0, Qt.AlignmentFlag.AlignRight)
        layout.addLayout(footer)

    def settings(self) -> SimilarityAnnotationSettings:
        return self._panel.settings()


class SequenceMapPanel(QWidget):
    """Embedded sequence map with annotation toolbar (preview Annotation tab)."""

    statusMessage = Signal(str)
    annotationJobProgress = Signal(str)
    annotationsSaved = Signal(object)  # Path to the written GFF3
    sequenceExtracted = Signal(object)  # Path to a newly written FASTA
    activeSequenceChanged = Signal(str)  # sequence id after multi-seq switch

    _MIN_PANEL_WIDTH = 160
    _MIN_PANEL_HEIGHT = 120

    def __init__(
        self,
        fasta_path: Path,
        parent: QWidget | None = None,
        *,
        autoload: bool = False,
    ) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setMinimumWidth(0)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._path = fasta_path.expanduser().resolve()
        # Sequence from the last load_path(), so the summary can be computed
        # without re-reading and re-parsing the FASTA from disk.
        self._loaded_data = None
        self._map = SequenceMapWidget(self, embedded=True)
        self._map.statusMessage.connect(self.statusMessage.emit)
        self._map._canvas.addFeatureRequested.connect(
            lambda _ft, _bl: self._on_add_feature_clicked()
        )
        self._map._canvas.copySequenceRequested.connect(self._on_copy_sequence_range)
        self._map._canvas.extractSequenceRequested.connect(self._on_extract_sequence)
        self._map.deleteFeaturesRequested.connect(
            lambda: self._delete_feature(confirm_empty=False)
        )
        self._map.copyFeatureSequencesRequested.connect(self._copy_selected_feature_sequences)
        self._map.extractFeatureSequencesRequested.connect(
            self._extract_selected_feature_sequences
        )
        self._map.modifyFeatureRequested.connect(self._modify_feature)
        self._map.changeFeatureColorRequested.connect(self._change_selected_feature_colors)
        self._map.annotationTracksMoveRequested.connect(self._apply_annotation_track_move)
        self._map.annotationFeaturesResizeRequested.connect(
            self._apply_annotation_feature_resize
        )
        self._map.attachFeatureImageRequested.connect(self._attach_image_to_selected_feature)
        self._map.openFeatureImagesRequested.connect(self._open_feature_images)
        self._map.featureImagesDropped.connect(self._on_feature_images_dropped)
        self._map._annot_track.set_image_marker_lookup(self._marked_feature_indices)
        from repeatviewer.widgets.cyto_images_ui import cyto_notifier

        cyto_notifier().changed.connect(self._on_cyto_images_changed)
        self.annotationJobProgress.connect(self._set_annotation_job_status)
        self._transfer_thread: QThread | None = None
        self._transfer_worker: Worker | None = None
        self._transfer_cancel = JobCancellation()
        self._transfer_preview_lifted = None
        self._transfer_preview_source_path: Path | None = None
        self._external_transfer_target: Path | None = None
        self._similarity_ann_settings = load_similarity_annotation_settings()
        self._auto_ann_settings = load_auto_annotation_settings()
        self._auto_ann_thread: QThread | None = None
        self._auto_ann_worker: Worker | None = None
        self._auto_ann_cancel = JobCancellation()
        self._auto_ann_reveal_kinds: tuple[str, ...] = ()
        self._sim_ann_thread: QThread | None = None
        self._sim_ann_worker: Worker | None = None
        self._sim_ann_cancel = JobCancellation()
        self._sim_ann_reveal_kinds: tuple[str, ...] = ()
        self._auto_all_updating = False
        self._similarity_all_updating = False
        self._similarity_db_rows: list[tuple] = []
        self._feature_undo_stack: list[list[SequenceFeature]] = []
        self._sequence_switch_updating = False
        self._sequence_search_completer: QCompleter | None = None

        # Geneious-style right rail (Display / Annotation / Find).
        self._rail = SequenceAnnotationRail(self)
        self._sequence_combo = self._rail.sequence_combo
        self._sequence_search = self._rail.sequence_search
        self._seq_count_label = self._rail.seq_count_label
        self._seq_row_container = self._rail
        self._controls_panel = self._rail
        self._controls_separator = QWidget(self)
        self._controls_separator.hide()
        self._add_row_container = self._rail
        self._auto_row_container = self._rail
        self._sim_row_container = self._rail
        self._action_row_container = self._rail
        self._auto_ann_all_check = self._rail.auto_ann_all_check
        self._auto_ann_run_btn = self._rail.auto_ann_run_btn
        self._auto_ann_settings_btn = self._rail.auto_ann_settings_btn
        self._similarity_ann_all_check = self._rail.similarity_ann_all_check
        self._similarity_ann_run_btn = self._rail.similarity_ann_run_btn
        self._similarity_ann_settings_btn = self._rail.similarity_ann_settings_btn
        self._auto_ann_checks: dict[AutoAnnotationMode, QCheckBox] = {}
        for label, mode in _AUTO_ANNOTATION_MODES:
            check = self._rail.auto_ann_checks[label]
            self._auto_ann_checks[mode] = check
            check.toggled.connect(self._on_auto_annotation_check_toggled)
        self._auto_ann_all_check.toggled.connect(self._on_auto_annotation_all_toggled)
        self._similarity_ann_all_check.toggled.connect(
            self._on_similarity_annotation_all_toggled
        )

        self._sequence_combo.currentIndexChanged.connect(self._on_sequence_combo_changed)
        self._sequence_search.returnPressed.connect(self._on_sequence_search_submitted)
        self._rail.addFeatureClicked.connect(self._on_add_feature_clicked)
        self._rail.addFeatureSettingsClicked.connect(
            self._open_add_feature_buttons_settings
        )
        self._rail.deleteClicked.connect(lambda: self._delete_feature(confirm_empty=True))
        self._rail.restoreClicked.connect(self._restore_features)
        self._rail.saveClicked.connect(self._save_annotations_from_button)
        self._rail.exportClicked.connect(self._export_annotation_view)
        self._rail.colorsClicked.connect(self._open_annotation_settings_dialog)
        self._rail.residueColorsChanged.connect(self._on_residue_colors_changed)
        self._rail.featureDetectClicked.connect(self._on_auto_run_clicked)
        self._rail.featureDetectSettingsClicked.connect(self._open_auto_annotation_settings)
        self._rail.similarityAnnotateClicked.connect(self._on_similarity_run_clicked)
        self._rail.similaritySettingsClicked.connect(self._open_similarity_annotation_settings)
        self._rail.transferSettingsClicked.connect(self._open_transfer_annotation_settings)
        self._rail.transferRequested.connect(self._on_transfer_annotations_clicked)
        self._rail.transferSourceChanged.connect(self._on_transfer_source_changed)
        self._rail.transferThresholdChanged.connect(self._apply_transfer_preview_filter)
        self._rail.transfer_allow_reverse.toggled.connect(
            self._on_transfer_reverse_changed
        )
        self._rail.transfer_target_browse.clicked.connect(self._on_transfer_browse_target)
        self._transfer_path_armed = False
        if hasattr(self._rail, "transfer_source_row"):
            self._rail.transfer_source_row.activated.connect(
                self._arm_transfer_source_path_fill
            )
            self._rail.transfer_source_row.set_default_tooltip(
                "Click a FASTA in Files or a folder in Folders to fill this path "
                "(click the field first to replace an existing path). Or use Browse."
            )
        self._rail.findRequested.connect(self._on_rail_find_requested)
        self._rail.panelCollapsedChanged.connect(self._on_rail_panel_collapsed_changed)
        self._map.featuresChanged.connect(self._sync_annotation_kind_filters)
        # "Save •" while annotation edits are not saved yet.
        self._map.featuresChanged.connect(self._sync_save_marker)
        self._find_items: list = []
        self._find_index: int = -1
        self._find_key: object | None = None
        self._wire_show_panel_controls()

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._map, stretch=1)
        map_rail_divider = QFrame(self)
        map_rail_divider.setObjectName("annotationMapRailDivider")
        map_rail_divider.setFrameShape(QFrame.Shape.VLine)
        map_rail_divider.setFrameShadow(QFrame.Shadow.Plain)
        map_rail_divider.setFixedWidth(1)
        map_rail_divider.setAutoFillBackground(True)
        divider_palette = map_rail_divider.palette()
        divider_palette.setColor(QPalette.ColorRole.Window, _MAP_RAIL_DIVIDER)
        divider_palette.setColor(QPalette.ColorRole.Dark, _MAP_RAIL_DIVIDER)
        map_rail_divider.setPalette(divider_palette)
        map_rail_divider.setStyleSheet(
            f"QFrame#annotationMapRailDivider {{ background-color: {_MAP_RAIL_DIVIDER.name()}; "
            f"border: none; max-width: 1px; }}"
        )
        layout.addWidget(map_rail_divider, stretch=0)
        layout.addWidget(self._rail, stretch=0)

        if autoload:
            self.load_path(self._path)

        self._map.set_show_similarity_database_name(
            self._similarity_ann_settings.show_database_name_on_track
        )
        self._install_annotation_shortcuts()
        self._refresh_transfer_blast_sources()
        self._sync_auto_annotation_all_check()
        self._refresh_similarity_database_rows()
        self._load_show_panel_settings()

    def _wire_show_panel_controls(self) -> None:
        rail = self._rail
        rail.showComplementChanged.connect(self._map.set_show_complement)
        rail.showTranslationChanged.connect(self._on_show_translation_toggled)
        rail.translationFramesChanged.connect(self._on_translation_frames_changed)
        rail.translationOptionsChanged.connect(self._on_translation_options_changed)
        rail.basesZoomedInChanged.connect(self._on_bases_zoomed_in_changed)
        rail.basesZoomedOutChanged.connect(self._map.set_letter_min_px)
        rail.annotationLabelSizeChanged.connect(self._map.set_annotation_label_size)
        rail.annotationThicknessChanged.connect(self._map.set_annotation_thickness)
        rail.annotationDistanceChanged.connect(self._map.set_annotation_distance)
        rail.sequenceNameSizeChanged.connect(self._on_sequence_name_size_changed)
        rail.baseNumberSizeChanged.connect(self._map.set_base_number_size)
        rail.splitTracksBySourceChanged.connect(self._on_split_tracks_by_source_changed)
        rail.annotationTrackKindsChanged.connect(
            self._on_annotation_track_kinds_changed
        )
        for signal in (
            rail.showComplementChanged,
            rail.showTranslationChanged,
            rail.translationFramesChanged,
            rail.translationOptionsChanged,
            rail.basesZoomedInChanged,
            rail.basesZoomedOutChanged,
            rail.annotationLabelSizeChanged,
            rail.annotationThicknessChanged,
            rail.annotationDistanceChanged,
            rail.sequenceNameSizeChanged,
            rail.baseNumberSizeChanged,
            rail.splitTracksBySourceChanged,
        ):
            signal.connect(lambda *_a: self._save_show_panel_settings())

    def _on_show_translation_toggled(self, checked: bool) -> None:
        if checked:
            self._map.set_translation_frames(self._rail.selected_translation_frames())
        else:
            self._map.set_translation_frames(frozenset())

    def _on_translation_frames_changed(self, frames: object) -> None:
        if self._rail.show_translation.isChecked():
            self._map.set_translation_frames(frames)
        else:
            self._map.set_translation_frames(frozenset())

    def _on_translation_options_changed(self) -> None:
        rail = self._rail
        self._map.set_translation_three_letter(rail.translation_three_letter())
        self._map.set_translation_table(rail.translation_genetic_code())
        if rail.show_translation.isChecked():
            self._map.set_translation_frames(rail.selected_translation_frames())

    def _on_bases_zoomed_in_changed(self, point_size: int) -> None:
        self._map.set_font_point_size(point_size)

    def _on_split_tracks_by_source_changed(self, enabled: bool) -> None:
        _save_split_tracks_by_source(bool(enabled))
        self._map.set_split_tracks_by_source(bool(enabled))

    def _on_annotation_track_kinds_changed(self, keys: object) -> None:
        if keys is None:
            visible: frozenset[str] | None = None
        else:
            visible = frozenset(str(k) for k in keys)
        self._map.set_visible_annotation_kinds(visible)

    def _sync_annotation_kind_filters(self) -> None:
        """Refresh Displayed annotations checkboxes and apply them to the map."""
        features = self._map.features()
        self._rail.set_available_annotation_kinds(features)
        self._map.set_visible_annotation_kinds(self._rail.visible_annotation_kinds())

    def _on_sequence_name_size_changed(self, point_size: int) -> None:
        self._rail.set_sequence_name_font_point_size(point_size)

    def _show_settings(self) -> QSettings:
        return QSettings("RepeatViewer", "RepeatViewer")

    def _load_show_panel_settings(self) -> None:
        s = self._show_settings()
        s.beginGroup("annotation_map_display")
        rail = self._rail
        blockers = [
            rail.show_complement,
            rail.show_translation,
            rail.split_tracks_by_source,
            rail.bases_zoomed_in,
            rail.bases_zoomed_out,
            rail.annotation_label_size,
            rail.annotation_thickness,
            rail.annotation_distance,
            rail.sequence_name_size,
            rail.base_number_size,
            *rail.frame_checks.values(),
        ]
        for widget in blockers:
            widget.blockSignals(True)
        try:
            rail.show_complement.setChecked(bool(s.value("show_complement", False, type=bool)))
            show_tr = bool(s.value("show_translation", False, type=bool))
            rail.show_translation.setChecked(show_tr)
            rail.translation_more_btn.setEnabled(show_tr)
            frames_raw = s.value("translation_frames", "1")
            try:
                frames = {
                    int(part)
                    for part in str(frames_raw).split(",")
                    if str(part).strip()
                }
            except ValueError:
                frames = {1}
            three_letter = bool(
                s.value("translation_three_letter", False, type=bool)
            )
            genetic_code = int(s.value("translation_genetic_code", 1))
            rail.set_translation_options(
                frames=frames or {1},
                three_letter=three_letter,
                genetic_code=genetic_code,
            )
            rail.split_tracks_by_source.setChecked(_load_split_tracks_by_source())
            rail.bases_zoomed_in.setValue(
                int(s.value("bases_zoomed_in", _DEFAULT_NUC_FONT_PT))
            )
            rail.bases_zoomed_out.setValue(int(s.value("bases_zoomed_out", 4)))
            rail.annotation_label_size.setValue(int(s.value("annotation_label_size", 7)))
            rail.annotation_thickness.setValue(int(s.value("annotation_thickness", 18)))
            rail.annotation_distance.setValue(int(s.value("annotation_distance", 0)))
            rail.sequence_name_size.setValue(int(s.value("sequence_name_size", 12)))
            rail.base_number_size.setValue(int(s.value("base_number_size", 8)))
        finally:
            for widget in blockers:
                widget.blockSignals(False)
        s.endGroup()
        self._map.set_show_complement(rail.show_complement.isChecked())
        self._on_show_translation_toggled(rail.show_translation.isChecked())
        self._on_translation_options_changed()
        self._map.set_split_tracks_by_source(rail.split_tracks_by_source.isChecked())
        self._map.set_font_point_size(rail.bases_zoomed_in.value())
        self._map.set_letter_min_px(rail.bases_zoomed_out.value())
        self._map.set_annotation_label_size(rail.annotation_label_size.value())
        self._map.set_annotation_thickness(rail.annotation_thickness.value())
        self._map.set_annotation_distance(rail.annotation_distance.value())
        self._rail.set_sequence_name_font_point_size(rail.sequence_name_size.value())
        self._map.set_base_number_size(rail.base_number_size.value())

    def _save_show_panel_settings(self) -> None:
        rail = self._rail
        s = self._show_settings()
        s.beginGroup("annotation_map_display")
        s.setValue("show_complement", rail.show_complement.isChecked())
        s.setValue("show_translation", rail.show_translation.isChecked())
        frames = sorted(
            frame
            for frame, check in rail.frame_checks.items()
            if check.isChecked()
        ) or [1]
        s.setValue("translation_frames", ",".join(str(f) for f in frames))
        s.setValue("translation_three_letter", rail.translation_three_letter())
        s.setValue("translation_genetic_code", rail.translation_genetic_code())
        s.setValue("bases_zoomed_in", rail.bases_zoomed_in.value())
        s.setValue("bases_zoomed_out", rail.bases_zoomed_out.value())
        s.setValue("annotation_label_size", rail.annotation_label_size.value())
        s.setValue("annotation_thickness", rail.annotation_thickness.value())
        s.setValue("annotation_distance", rail.annotation_distance.value())
        s.setValue("sequence_name_size", rail.sequence_name_size.value())
        s.setValue("base_number_size", rail.base_number_size.value())
        s.endGroup()

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return QSize(self._MIN_PANEL_WIDTH, self._MIN_PANEL_HEIGHT)

    def bind_annotation_zoom_bar(self, bar) -> None:
        self._map.bind_annotation_zoom_bar(bar)

    def showEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().showEvent(event)
        self._refresh_transfer_blast_sources()

    def _refresh_annotation_button_styles(self) -> None:
        return

    def _install_annotation_shortcuts(self) -> None:
        # Widget-scoped shortcuts only — never install a Python filter on
        # QApplication (PySide6 can segfault on macOS during setProperty).
        for key in (QKeySequence.StandardKey.Delete, QKeySequence(Qt.Key.Key_Backspace)):
            shortcut = QShortcut(key, self)
            shortcut.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            shortcut.activated.connect(lambda: self._delete_feature(confirm_empty=False))

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802
        if _is_delete_key(event.key()):
            self._map._sync_feature_selection_from_annot_track()
            if self._map.selected_feature_indices():
                self._delete_feature(confirm_empty=False)
                event.accept()
                return
        super().keyPressEvent(event)

    def hideEvent(self, event) -> None:  # noqa: N802
        super().hideEvent(event)

    def eventFilter(self, watched, event) -> bool:  # noqa: N802 (Qt API)
        return super().eventFilter(watched, event)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)

    def _sync_annotation_control_rows(self) -> None:
        return

    def _sync_run_settings_boxes(self) -> None:
        return

    def _sync_similarity_row_layout(self) -> None:
        return

    def _sync_add_feature_row(self) -> None:
        return

    def _sync_auto_feature_columns(self) -> None:
        return

    def _rebuild_add_feature_buttons(
        self, buttons: list[tuple[str, str]] | None = None
    ) -> None:
        if hasattr(self, "_rail"):
            self._rail.populate_manual_feature_types(buttons)

    def _open_add_feature_buttons_settings(self) -> None:
        dialog = AddFeatureButtonsDialog(parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._rebuild_add_feature_buttons(dialog.buttons())
        self._broadcast_shared_annotation_settings_changed()
        self.statusMessage.emit("Feature types updated.")

    def _sync_auto_annotation_row_alignment(self) -> None:
        return

    def _on_auto_annotation_check_toggled(self, _checked: bool) -> None:
        if self._auto_all_updating:
            return
        self._sync_auto_annotation_all_check()

    def _sync_auto_annotation_all_check(self) -> None:
        # "All" covers LTR/TIR/Tandem/ORF; TSD (MSA-only) and PBS (needs LTR +
        # BLAST) stay independent and off by default.
        checks = [
            check
            for mode, check in self._auto_ann_checks.items()
            if mode not in {AutoAnnotationMode.TSD, AutoAnnotationMode.PBS}
        ]
        if not checks:
            return
        all_checked = all(check.isChecked() for check in checks)
        self._auto_all_updating = True
        self._auto_ann_all_check.setChecked(all_checked)
        self._auto_all_updating = False

    def _on_auto_annotation_all_toggled(self, checked: bool) -> None:
        if self._auto_all_updating:
            return
        self._auto_all_updating = True
        for mode, check in self._auto_ann_checks.items():
            if mode in {AutoAnnotationMode.TSD, AutoAnnotationMode.PBS}:
                continue
            check.setChecked(checked)
        self._auto_all_updating = False

    def _on_similarity_db_check_toggled(self, _checked: bool) -> None:
        if self._similarity_all_updating:
            return
        self._sync_similarity_annotation_all_check()

    def _sync_similarity_annotation_all_check(self) -> None:
        checks = [check for _entry, check in self._similarity_db_rows]
        if not checks:
            self._similarity_all_updating = True
            self._similarity_ann_all_check.setChecked(False)
            self._similarity_all_updating = False
            return
        all_checked = all(check.isChecked() for check in checks)
        self._similarity_all_updating = True
        self._similarity_ann_all_check.setChecked(all_checked)
        self._similarity_all_updating = False

    def _on_similarity_annotation_all_toggled(self, checked: bool) -> None:
        if self._similarity_all_updating:
            return
        self._similarity_all_updating = True
        for _entry, check in self._similarity_db_rows:
            check.setChecked(checked)
        self._similarity_all_updating = False

    def _on_add_feature_overflow_activated(self, index: int) -> None:
        feature_type = self._add_feature_overflow.itemData(index)
        label = self._add_feature_overflow.itemText(index)
        if feature_type:
            self._add_feature_of_type(str(feature_type), feature_label=label)

    def refresh_similarity_databases(self) -> None:
        self._refresh_similarity_database_rows()
        self._refresh_transfer_blast_sources()

    def _on_show_online_blast_toggled(self, checked: bool) -> None:
        del checked

    def _on_residue_colors_changed(self) -> None:
        self._map.refresh_nucleotide_colors()
        self._map.update()
        self._broadcast_shared_annotation_settings_changed()

    def _open_annotation_settings_dialog(self) -> None:
        dialog = AnnotationElementColorsDialog(
            parent=self,
            font_point_size=self._map.font_point_size(),
            split_tracks_by_source=self._map.split_tracks_by_source(),
            min_font_pt=_MIN_NUC_FONT_PT,
            max_font_pt=_MAX_NUC_FONT_PT,
            nucleotide_color_context="single",
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        if dialog.apply_to_current_sequence():
            self._apply_element_colors_to_map(
                dialog.old_feature_type_colors(),
                dialog.new_feature_type_colors(),
                dialog.old_database_colors(),
                dialog.new_database_colors(),
            )
        self._rail.sync_color_scheme_combos()
        self._rail.refresh_feature_detect_check_colors()
        self._rail.refresh_manual_feature_check_colors()
        self._map.refresh_nucleotide_colors()
        self._map.update()
        self._broadcast_shared_annotation_settings_changed()
        self.statusMessage.emit("Annotation settings updated.")

    def _open_element_colors_dialog(self) -> None:
        # Back-compat alias for older call sites.
        self._open_annotation_settings_dialog()

    def _apply_element_colors_to_map(
        self,
        old_type_colors: dict[str, str],
        new_type_colors: dict[str, str],
        old_db_colors: dict[str, str],
        new_db_colors: dict[str, str],
    ) -> None:
        updated: list[SequenceFeature] = []
        changed = False
        for feat in self._map.features():
            color = feat.color
            if feat.feature_type == "domain":
                for db_name, old_color in old_db_colors.items():
                    if color == old_color:
                        color = new_db_colors.get(db_name, color)
                        break
            else:
                old_color = old_type_colors.get(feat.feature_type, "")
                if not color or color == old_color:
                    color = new_type_colors.get(
                        feat.feature_type,
                        feature_type_color(feat.feature_type),
                    )
            if color != feat.color:
                changed = True
            updated.append(
                SequenceFeature(
                    start=feat.start,
                    end=feat.end,
                    name=feat.name,
                    feature_type=feat.feature_type,
                    color=color,
                    strand=feat.strand,
                    track=feat.track,
                    parent=feat.parent,
                    source=feat.source,
                    pident=feat.pident,
                    subject_id=feat.subject_id,
                    coverage=feat.coverage,
                    source_start=feat.source_start,
                    source_end=feat.source_end,
                )
            )
        if changed:
            self._map.set_features(updated)

    def refresh_nucleotide_colors(self) -> None:
        self._map.refresh_nucleotide_colors()

    def open_annotation_tab(self, *, reload: bool = False) -> None:
        """Load annotations and fit the sequence to the visible window width."""
        if reload or self._map.data() is None:
            # load_path(fit_after_load=True) already schedules the fit_window() call
            # (see SequenceMapWidget.load()); do not schedule a second, redundant one
            # here - it used to race a Go-to jump landing on this same tab.
            self.load_path(self._path, fit_after_load=True)
            QTimer.singleShot(50, self._map.refresh_layout)
        else:
            # Already loaded — only re-sync scroll/viewport geometry (cheap).
            QTimer.singleShot(0, self._map.refresh_layout)

    # -- cytogenetic images ------------------------------------------------

    @staticmethod
    def _feature_ref(feat: SequenceFeature):
        from repeatviewer.backends.cyto_images import FeatureRef

        return FeatureRef(feat.name or feat.feature_type, feat.start, feat.end)

    def _marked_feature_indices(self, features: list[SequenceFeature]) -> set[int]:
        """Indices of ``features`` that carry cytogenetic images.

        Exact name + coordinates first; a link whose annotation was moved or
        resized still marks it when exactly one annotation has that name.
        """
        from repeatviewer.backends.cyto_images import images_for_sequence

        images = images_for_sequence(self._path)
        refs = [
            target.feature
            for image in images
            for target in image.targets_for_sequence(self._path.name)
            if target.feature is not None
        ]
        if not refs:
            return set()
        exact: dict[tuple[str, int, int], int] = {}
        by_name: dict[str, list[int]] = {}
        for idx, feat in enumerate(features):
            name = feat.name or feat.feature_type
            exact.setdefault((name, feat.start, feat.end), idx)
            by_name.setdefault(name, []).append(idx)
        marked: set[int] = set()
        for ref in refs:
            idx = exact.get((ref.name, ref.start, ref.end))
            if idx is not None:
                marked.add(idx)
                continue
            candidates = by_name.get(ref.name, [])
            if len(candidates) == 1:
                marked.add(candidates[0])
        return marked

    def _on_cyto_images_changed(self, folder: str) -> None:
        if folder == str(self._path.parent):
            self._map._annot_track.update()

    def _single_selected_feature(self) -> SequenceFeature | None:
        data = self._map.data()
        if data is None:
            return None
        self._map._sync_feature_selection_from_annot_track()
        indices = sorted(self._map.selected_feature_indices())
        if len(indices) != 1 or not (0 <= indices[0] < len(data.features)):
            return None
        return data.features[indices[0]]

    def _attach_image_to_selected_feature(self) -> None:
        from repeatviewer.widgets.cyto_images_ui import attach_images

        feat = self._single_selected_feature()
        if feat is None:
            QMessageBox.information(
                self, "Attach cytogenetic image", "Select a single annotation first."
            )
            return
        attach_images(self, self._path, [], feature=self._feature_ref(feat))

    def _on_feature_images_dropped(self, index: int, paths: list) -> None:
        from repeatviewer.widgets.cyto_images_ui import attach_images

        data = self._map.data()
        if data is None or not (0 <= index < len(data.features)):
            return
        feat = data.features[index]
        attach_images(self, self._path, [Path(p) for p in paths], feature=self._feature_ref(feat))

    def _open_feature_images(self, index: int) -> None:
        from repeatviewer.widgets.cyto_images_ui import open_images_for_feature

        data = self._map.data()
        if data is None or not (0 <= index < len(data.features)):
            return
        ref = self._feature_ref(data.features[index])
        open_images_for_feature(self, self._path, ref.name, ref.start, ref.end)

    def show_feature(self, ref) -> bool:
        """Select and jump to the annotation ``ref`` (a FeatureRef). True if found."""
        if self._map.data() is None:
            self.load_path(self._path, fit_after_load=False)
        data = self._map.data()
        if data is None or ref is None:
            return False
        marked_exact = None
        same_name = []
        for idx, feat in enumerate(data.features):
            name = feat.name or feat.feature_type
            if (name, feat.start, feat.end) == (ref.name, ref.start, ref.end):
                marked_exact = idx
                break
            if name == ref.name:
                same_name.append(idx)
        idx = marked_exact if marked_exact is not None else (same_name[0] if len(same_name) == 1 else None)
        if idx is None:
            self.jump_to_region(ref.start, ref.end)
            return False
        feat = data.features[idx]
        self.jump_to_region(feat.start, feat.end)
        self._map.set_selected_features({idx})
        self._map._annot_track._selected_indices = {idx}
        self._map._annot_track.update()
        return True

    def jump_to_region(self, start: int, end: int) -> None:
        """Zoom/scroll to a bp range without reloading or fit-to-window reset."""
        if self._map.data() is None:
            self.load_path(self._path, fit_after_load=False)
        self._map.select_by_position(start, end)

    def select_region(self, start: int, end: int) -> None:
        """Select and zoom/scroll the map to a (1-based, inclusive) bp range."""
        self.jump_to_region(start, end)

    # ------------------------------------------------ Edit > Find / Go To / Copy

    def active_sequence_id(self) -> str:
        data = self._loaded_data or self._map.data()
        return str(getattr(data, "sequence_id", "") or "")

    def active_sequence_length(self) -> int:
        data = self._map.data()
        return len(data.sequence) if data is not None else 0

    def is_too_long_for_display(self) -> bool:
        return self._map.is_too_long_for_display()

    def show_find_hit(
        self,
        seq_id: str,
        start: int,
        end: int,
        hits: list[tuple[int, int]],
        *,
        whole_sequence: bool = False,
        select_annotations: list | None = None,
    ) -> bool:
        """Highlight ``hits`` and select/show [start, end] of ``seq_id``. False if not shown."""
        if self._map.data() is None:
            self.load_path(self._path, fit_after_load=False)
        if seq_id and not self.ensure_active_sequence(seq_id):
            return False
        if self._map.data() is None or self._map.is_too_long_for_display():
            return False
        if select_annotations is not None:
            indices = self._feature_indices_for_annotation_hits(select_annotations)
            # Annotation search: select features only (no find underlays / base band).
            self._map.set_find_hits([])
            self._map._canvas.clear_selection()
            # Scroll/zoom without treating the hit as a base selection.
            a, b = sorted((int(start), int(end)))
            self._map._scroll_to_bp_range(a, b)
            span = max(1, b - a + 1)
            viewport_w = self._map._scroll_viewport_width()
            usable_w = max(1.0, float(viewport_w - _LEFT_MARGIN - _RIGHT_MARGIN))
            visible_bp = usable_w / max(self._map._canvas.pixels_per_bp(), 1e-9)
            if span > 0.9 * visible_bp:
                target_ppb = usable_w / max(3 * span, 60)
                target_ppb = max(_FIT_MIN_PPB, min(self._map._canvas.zoom_ppb_max(), target_ppb))
                self._map._apply_zoom_ppb(target_ppb, sync_slider=True)
                self._map._scroll_to_bp_range(a, b)
            self._map.set_selected_features(indices)
            return True
        self._map.set_find_hits(hits)
        if whole_sequence:
            self._map.fit_window()
            self._map._canvas.set_selection_range(start, end)
        else:
            self._map.show_region(start, end)
        return True

    def _feature_indices_for_annotation_hits(self, hits: list) -> set[int]:
        """Map Find annotation hits to indices in the active sequence feature list."""
        features = self._map.features()
        if not features or not hits:
            return set()
        indices: set[int] = set()
        for hit in hits:
            a, b = sorted(
                (int(getattr(hit, "start", 0)), int(getattr(hit, "end", 0)))
            )
            needle = str(getattr(hit, "label", "") or "").strip()
            for index, feat in enumerate(features):
                fa, fb = sorted((int(feat.start), int(feat.end)))
                if fa != a or fb != b:
                    continue
                name = str(feat.name or "")
                ftype = str(feat.feature_type or "")
                if (
                    not needle
                    or needle == name
                    or needle == ftype
                    or needle.lower() in name.lower()
                    or needle.lower() in ftype.lower()
                ):
                    indices.add(index)
        return indices

    def clear_find_hits(self) -> None:
        self._map.set_find_hits([])
        self._find_items = []
        self._find_index = -1
        self._find_key = None

    def _on_rail_panel_collapsed_changed(self, _collapsed: bool) -> None:
        """After the control rail hides/shows, re-fit if still at Fit zoom."""
        # Wait for the layout pass so the map has its new width.
        QTimer.singleShot(0, self._map.refresh_layout)

    def _on_rail_find_requested(self, direction: int) -> None:
        """Find tab: Geneious-style next / previous / find-all on this FASTA."""
        from repeatviewer.backends.find_in_file import find_in_fasta
        from repeatviewer.backends.motif_search import parse_position, query_problem

        query = self._rail.find_query.text().strip()
        if not query:
            self._rail.set_find_status("Type something to find.")
            self._rail.find_query.setFocus()
            return
        position = parse_position(query)
        if position is not None and (
            position[0] is not None or not any(c.isalpha() for c in query)
        ):
            name, start, end = position
            if name and not self.ensure_active_sequence(name):
                self._rail.set_find_status(f"No sequence named {name}.")
                return
            self._map.show_region(start, end)
            where = f"{start:,}" if start == end else f"{start:,}-{end:,}"
            self._rail.set_find_status(
                f"Went to position {where}"
                + (f" of {name}" if name else "")
                + "."
            )
            return

        request = self._rail.find_request()
        if request.in_sequence:
            problem = query_problem(query, request.options)
            if problem and not request.in_names and not request.in_annotations:
                self._rail.set_find_status(problem)
                return
        try:
            stat = self._path.stat()
        except OSError:
            self._rail.set_find_status("The file is no longer there.")
            return
        key = (str(self._path), stat.st_mtime_ns, request)
        if self._find_key != key or not self._find_items:
            features = None
            data = self._loaded_data or self._map.data()
            if request.in_annotations and data is not None:
                features = getattr(data, "features_by_id", None) or {}
                # Include the active sequence's in-memory features.
                active = str(getattr(data, "sequence_id", "") or "")
                if active:
                    features = dict(features)
                    features[active] = list(self._map.features())
            try:
                QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
                result = find_in_fasta(
                    self._path, request, features_by_id=features
                )
            except (OSError, ValueError) as exc:
                self._rail.set_find_status(f"Could not search: {exc}")
                return
            finally:
                QApplication.restoreOverrideCursor()
            self._find_items = list(result.items)
            self._find_key = key
            self._find_index = -1
            if result.sequence_skipped and not self._find_items:
                self._rail.set_find_status(result.sequence_skipped)
                self.clear_find_hits()
                return
            if not self._find_items:
                self._rail.set_find_status("No matches.")
                self.clear_find_hits()
                self._find_key = key
                return

        if direction == 0:
            self._find_index = 0
        else:
            step = 1 if direction >= 0 else -1
            if self._find_index < 0:
                self._find_index = 0 if step > 0 else len(self._find_items) - 1
            else:
                self._find_index = (self._find_index + step) % len(self._find_items)

        item = self._find_items[self._find_index]
        seq_hits = [
            (h.start, h.end)
            for h in self._find_items
            if h.seq_id == item.seq_id and h.kind in {"sequence", "annotation"}
        ]
        if item.kind == "annotation":
            ann_hits = [h for h in self._find_items if h.kind == "annotation"]
            ok = self.show_find_hit(
                item.seq_id,
                item.start,
                item.end,
                [],
                select_annotations=ann_hits,
            )
        else:
            ok = self.show_find_hit(
                item.seq_id,
                item.start,
                item.end,
                seq_hits or [(item.start, item.end)],
                whole_sequence=(item.kind == "name"),
            )
        if not ok:
            self._rail.set_find_status("Could not show that hit.")
            return
        kind = item.kind
        label = f" ({item.label})" if item.label else ""
        self._rail.set_find_status(
            f"{self._find_index + 1} of {len(self._find_items)} — "
            f"{item.seq_id}:{item.start:,}-{item.end:,} [{kind}]{label}"
        )

    def copy_selection_to_clipboard(self) -> str | None:
        """Cmd+C: selected annotations as before, else the selected bases as text."""
        data = self._loaded_data or self._map.data()
        if data is None or not getattr(data, "sequence", ""):
            return None
        if self._selected_features():
            self._copy_selected_feature_sequences()
            return "Copied the selected annotation sequence(s)"
        selection = self._map.selection_range()
        if selection is None:
            return None
        a, b = selection
        fragment = data.sequence[a - 1 : b]
        if not fragment:
            return None
        QApplication.clipboard().setText(fragment)
        return f"Copied {len(fragment):,} bp ({data.sequence_id}:{a:,}-{b:,})"

    def add_features_and_save(self, features: list[SequenceFeature]) -> bool:
        """Append features to the shown sequence and save the annotations."""
        if self._map.data() is None or not features:
            return False
        current = list(self._map.features())
        current.extend(features)
        self._map.set_features(current)
        return self.save_annotations(silent=True)

    def fasta_path(self) -> Path:
        return self._path

    def is_dirty(self) -> bool:
        return self._map.is_dirty()

    def load_path(self, path: Path, *, fit_after_load: bool = True) -> None:
        try:
            resolved = path.expanduser().resolve()
        except OSError:
            resolved = path.expanduser()
        if resolved != self._path:
            self._map.set_find_hits([])  # hits belong to the previous file
            self._stop_transfer_session(
                status="Transfer stopped (file changed)."
            )
        self._path = resolved
        QTimer.singleShot(0, self._sync_save_marker)
        try:
            data = load_sequence_map_data(self._path)
        except (OSError, ValueError) as exc:
            self._loaded_data = None
            self._populate_sequence_combo(None)
            self._sync_too_long_controls()
            self.statusMessage.emit(str(exc))
            return
        self._loaded_data = data
        self._populate_sequence_combo(data)
        # Patch legacy Online BLAST features missing pident from *_ob.txt.
        try:
            from repeatviewer.backends.ncbi_blast_online import (
                enrich_online_blast_features_from_ob,
                online_blast_output_path,
            )

            ob_path = online_blast_output_path(self._path)
            if ob_path.is_file() and data.features:
                enriched = enrich_online_blast_features_from_ob(list(data.features), ob_path)
                data.features = enriched
                rid = str(data.sequence_id or "")
                if rid and data.features_by_id is not None:
                    data.features_by_id[rid] = list(enriched)
        except Exception:  # noqa: BLE001
            pass
        from repeatviewer.session_restore import lookup_annotation_view

        saved_view = lookup_annotation_view(
            self._path,
            sequence_id=str(data.sequence_id or ""),
        )
        fit = fit_after_load and saved_view is None
        self._map.load(data, fit_after_load=fit)
        self._rail.set_sequence_kind(bool(data.is_protein))
        self._sync_annotation_kind_filters()
        if saved_view is not None:
            self._map.restore_view_state_for_session(saved_view)
        self._feature_undo_stack.clear()
        self._sync_too_long_controls()
        if not self._map.is_too_long_for_display():
            self.refresh_sequence_summary()

    def _sync_too_long_controls(self) -> None:
        """Hide the right rail when the map is blocked for a too-long sequence."""
        too_long = self._map.is_too_long_for_display()
        multi = bool(
            self._loaded_data is not None and len(self._loaded_data.record_ids) > 1
        )
        # Keep the rail for multi-FASTA so a shorter record can be selected.
        self._rail.setVisible(not too_long or multi)

    def _populate_sequence_combo(self, data: SequenceMapData | None) -> None:
        self._sequence_switch_updating = True
        try:
            self._sequence_combo.clear()
            self._sequence_search.clear()
            if data is None or not data.record_ids:
                self._seq_count_label.clear()
                self._sequence_search_completer = None
                self._refresh_transfer_sequence_lists()
                return
            for rid in data.record_ids:
                n_feat = len(data.features_by_id.get(rid) or [])
                label = rid if n_feat == 0 else f"{rid}  ({n_feat})"
                self._sequence_combo.addItem(label, rid)
            idx = self._sequence_combo.findData(data.sequence_id)
            if idx < 0:
                idx = 0
            self._sequence_combo.setCurrentIndex(idx)
            self._sequence_combo.setToolTip(data.sequence_id)
            total = len(data.record_ids)
            if total > 1:
                self._seq_count_label.setText(f"{total:,} sequences")
            else:
                self._seq_count_label.clear()
            self._refresh_transfer_sequence_lists()
            completer = QCompleter(list(data.record_ids), self._sequence_search)
            completer.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
            completer.setFilterMode(Qt.MatchFlag.MatchContains)
            completer.setCompletionMode(QCompleter.CompletionMode.PopupCompletion)
            completer.activated.connect(self._on_sequence_search_activated)
            self._sequence_search.setCompleter(completer)
            self._sequence_search_completer = completer
        finally:
            self._sequence_switch_updating = False

    @Slot(int)
    def _on_sequence_combo_changed(self, index: int) -> None:
        if self._sequence_switch_updating or index < 0:
            return
        seq_id = self._sequence_combo.itemData(index)
        if not seq_id:
            return
        self.ensure_active_sequence(str(seq_id))

    @Slot(str)
    def _on_sequence_search_activated(self, text: str) -> None:
        self._apply_sequence_search(str(text))

    @Slot()
    def _on_sequence_search_submitted(self) -> None:
        self._apply_sequence_search(self._sequence_search.text())

    def _apply_sequence_search(self, query: str) -> None:
        data = self._loaded_data
        if data is None:
            return
        text = str(query or "").strip()
        if not text:
            return
        # Exact id, or first contains match (case-insensitive).
        lowered = text.lower()
        target = ""
        if text in data.record_ids:
            target = text
        else:
            for rid in data.record_ids:
                if rid.lower() == lowered:
                    target = rid
                    break
            if not target:
                for rid in data.record_ids:
                    if lowered in rid.lower():
                        target = rid
                        break
        if not target:
            self.statusMessage.emit(f"No sequence matching “{text}”.")
            return
        if self.ensure_active_sequence(target):
            self._sequence_search.clear()

    def ensure_active_sequence(self, seq_id: str) -> bool:
        """Switch the map to ``seq_id`` (saving in-memory edits for the previous id)."""
        data = self._loaded_data
        if data is None:
            return False
        wanted = str(seq_id or "").strip()
        if not wanted:
            return False
        if data.sequence_id == wanted:
            return True
        # Soft-match against known ids.
        target = wanted
        if wanted not in data.record_ids:
            from repeatviewer.backends.gff_features import _seqid_matches

            target = next(
                (rid for rid in data.record_ids if _seqid_matches(rid, wanted)),
                "",
            )
            if not target:
                self.statusMessage.emit(f"Sequence “{wanted}” is not in this FASTA.")
                return False
        if data.sequence_id == target:
            return True

        # Keep unsaved canvas edits in the multi-seq store (disk write still via Save).
        data.sync_active_features(self._map.features())
        try:
            data.switch_to(target)
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, "Sequence", str(exc))
            return False
        self._loaded_data = data
        self._sequence_switch_updating = True
        try:
            idx = self._sequence_combo.findData(data.sequence_id)
            if idx >= 0:
                self._sequence_combo.setCurrentIndex(idx)
        finally:
            self._sequence_switch_updating = False
        self._map.load(data, fit_after_load=True)
        # Leaving this sequence ends an active Transfer preview session.
        self._stop_transfer_session(status="Transfer stopped (sequence changed).")
        self._rail.set_sequence_kind(bool(data.is_protein))
        self._feature_undo_stack.clear()
        self._sync_too_long_controls()
        if not self._map.is_too_long_for_display():
            self.refresh_sequence_summary()
        self._sequence_combo.setToolTip(data.sequence_id)
        self.activeSequenceChanged.emit(data.sequence_id)
        return True

    def _next_extraction_name(self, seq_id: str) -> str:
        """Default name: <seq>_extraction_N, skipping names already on disk."""
        base = seq_id or self._path.stem
        parent = self._path.parent
        index = 1
        while (parent / f"{base}_extraction_{index}.fasta").exists():
            index += 1
        return f"{base}_extraction_{index}"

    def _next_feature_extraction_name(self, seq_id: str, feat: SequenceFeature) -> str:
        """Default name: <seq>_<type>, then <seq>_<type>2, …"""
        stem = _feature_export_stem(seq_id, feat)
        parent = self._path.parent
        if not (parent / f"{stem}.fasta").exists() and not (parent / f"{stem}.fa").exists():
            return stem
        index = 2
        while True:
            candidate = f"{stem}{index}"
            if (
                not (parent / f"{candidate}.fasta").exists()
                and not (parent / f"{candidate}.fa").exists()
            ):
                return candidate
            index += 1

    def _write_extracted_fasta(
        self,
        *,
        sequence: str,
        sequence_id: str,
        start: int,
        end: int,
        name: str,
        reverse: bool = False,
    ) -> Path | None:
        lo = max(1, min(int(start), int(end)))
        hi = min(len(sequence), max(int(start), int(end)))
        fragment = sequence[lo - 1 : hi]
        if not fragment:
            QMessageBox.warning(self, "Extract Sequence", "The selection is empty.")
            return None
        if reverse:
            fragment = reverse_complement(fragment)
        safe = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in name)
        out_path = self._path.parent / f"{safe}.fasta"
        if out_path.exists():
            answer = QMessageBox.question(
                self,
                "Extract Sequence",
                f"{out_path.name} already exists. Overwrite it?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if answer != QMessageBox.StandardButton.Yes:
                return None
        strand = "-" if reverse else "+"
        header = f">{name} {sequence_id}:{lo}-{hi}({strand})"
        try:
            with out_path.open("w", encoding="utf-8") as handle:
                handle.write(header + "\n")
                for offset in range(0, len(fragment), 60):
                    handle.write(fragment[offset : offset + 60] + "\n")
        except OSError as exc:
            QMessageBox.warning(
                self, "Extract Sequence", f"Could not write {out_path.name}:\n{exc}"
            )
            return None
        self.statusMessage.emit(
            f"Extracted {lo:,}-{hi:,} ({len(fragment):,} bp) to {out_path.name}"
        )
        self.sequenceExtracted.emit(out_path)
        return out_path

    @Slot(int, int)
    def _on_copy_sequence_range(self, start: int, end: int) -> None:
        """Copy the selected nucleotide/protein range to the clipboard."""
        data = self._loaded_data or self._map.data()
        if data is None or not getattr(data, "sequence", ""):
            QMessageBox.warning(
                self, "Copy Sequence", "The sequence is not loaded yet."
            )
            return
        lo = min(int(start), int(end))
        hi = max(int(start), int(end))
        fragment = data.sequence[lo - 1 : hi]
        if not fragment:
            QMessageBox.warning(self, "Copy Sequence", "The selection is empty.")
            return
        QApplication.clipboard().setText(fragment)
        unit = "aa" if bool(getattr(data, "is_protein", False)) else "bp"
        self.statusMessage.emit(
            f"Copied {len(fragment):,} {unit} ({data.sequence_id}:{lo:,}-{hi:,})"
        )

    @Slot(int, int)
    def _on_extract_sequence(self, start: int, end: int) -> None:
        """Write the selected range out as a new single-record FASTA."""
        data = self._loaded_data
        if data is None or not getattr(data, "sequence", ""):
            QMessageBox.warning(
                self, "Extract Sequence", "The sequence is not loaded yet."
            )
            return

        sequence = data.sequence
        lo = min(int(start), int(end))
        hi = max(int(start), int(end))
        base = _safe_filename_part(data.sequence_id) or self._path.stem or "sequence"
        dialog = _ExtractSequenceDialog(
            sequence_id=data.sequence_id,
            sequence_length=len(sequence),
            start=lo,
            end=hi,
            default_name=f"{base}_{lo}_{hi}",
            reverse=False,
            parent=self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        chosen_start, chosen_end, name, reverse = dialog.values()
        self._write_extracted_fasta(
            sequence=sequence,
            sequence_id=data.sequence_id,
            start=chosen_start,
            end=chosen_end,
            name=name,
            reverse=reverse,
        )

    def select_by_position(self, start: int, end: int) -> None:
        """Select a 1-based inclusive range (delegates to the map widget)."""
        self._map.select_by_position(int(start), int(end))

    def selection_range(self) -> tuple[int, int] | None:
        """Current selection as a 1-based inclusive range, or None."""
        return self._map.selection_range()

    def refresh_sequence_summary(self) -> None:
        if self._map.is_too_long_for_display():
            self._map.set_sequence_summary_text("")
            self._sequence_summary_cache_key = None
            return
        if not load_show_sequence_summary():
            self._map.set_sequence_summary_text("")
            self._sequence_summary_cache_key = None
            return

        # Prefer the sequence we already hold: re-reading the file here meant a
        # 500 Mbp scaffold was parsed twice per open, doubling both the wait and
        # peak memory.
        data = self._loaded_data
        if data is not None and data.fasta_path == self._path:
            cache_key = (data.sequence_id, len(data.sequence), id(data.sequence))
            if cache_key == getattr(self, "_sequence_summary_cache_key", None):
                return
            self._map.set_sequence_summary(
                summarize_loaded_sequence(data.sequence_id, data.sequence)
            )
            self._sequence_summary_cache_key = cache_key
            return

        try:
            mtime = self._path.stat().st_mtime_ns
        except OSError:
            mtime = None
        cache_key = (str(self._path), mtime)
        if cache_key == getattr(self, "_sequence_summary_cache_key", None):
            return
        try:
            summary = summarize_single_sequence_fasta(self._path)
        except (OSError, ValueError):
            self._map.set_sequence_summary_text("")
            self._sequence_summary_cache_key = None
            return
        self._map.set_sequence_summary(summary)
        self._sequence_summary_cache_key = cache_key

    def reload_annotations(self) -> None:
        """Reload features from the companion GFF3 on disk (Annotation tab open)."""
        self.open_annotation_tab(reload=True)

    def reload_shared_annotation_settings(self) -> None:
        """Reload auto/similarity settings and colors shared with TE annotation."""
        self._auto_ann_settings = load_auto_annotation_settings()
        self._similarity_ann_settings = load_similarity_annotation_settings()
        self._map.set_show_similarity_database_name(
            self._similarity_ann_settings.show_database_name_on_track
        )
        split_on = _load_split_tracks_by_source()
        self._map.set_split_tracks_by_source(split_on)
        rail = getattr(self, "_rail", None)
        if rail is not None and hasattr(rail, "split_tracks_by_source"):
            rail.split_tracks_by_source.blockSignals(True)
            rail.split_tracks_by_source.setChecked(split_on)
            rail.split_tracks_by_source.blockSignals(False)
        self._map.set_show_online_blast(False)
        self._rebuild_add_feature_buttons()
        self.refresh_similarity_databases()
        self.refresh_nucleotide_colors()
        if rail is not None:
            rail.refresh_feature_detect_check_colors()
            rail.refresh_manual_feature_check_colors()
        self._map.update()

    def _broadcast_shared_annotation_settings_changed(self) -> None:
        """Push disk settings into other open annotation panels."""
        window = self.window()
        if window is None:
            return
        for panel in window.findChildren(SequenceMapPanel):
            if panel is self:
                continue
            panel.reload_shared_annotation_settings()

    def save_annotations(self, *, silent: bool = False) -> bool:
        return self._save_annotations(silent=silent)

    def has_unsaved_annotation_changes(self) -> bool:
        """Edits made in this view that Save has not written yet.

        Unlike ``_annotation_changes_pending`` this ignores "features loaded
        from a GFF but no JSON yet", which is not an edit and must not trigger
        the leave-file prompt.
        """
        try:
            if self._map.data() is None:
                return False
            return bool(
                self._map.is_dirty() or self._visible_transfer_prediction_count() > 0
            )
        except RuntimeError:
            return False

    def discard_annotation_changes(self) -> None:
        """Don't Save at the leave-file prompt: forget the unsaved edits.

        The next time this file is shown it is reloaded from the JSON on disk.
        """
        try:
            self._clear_transfer_preview_state()
        except Exception:  # noqa: BLE001 - discarding must not fail
            pass
        self._map.mark_clean()
        self._sync_save_marker()

    def _sync_save_marker(self) -> None:
        rail = getattr(self, "_rail", None)
        button = getattr(rail, "save_btn", None) if rail is not None else None
        if button is None:
            return
        try:
            dirty = self.has_unsaved_annotation_changes()
            button.setText("Save \u2022" if dirty else "Save")
            button.setToolTip(
                "Unsaved annotation changes. Save writes them to the "
                ".seqmap_RV.json next to the FASTA."
                if dirty
                else "Save annotations to the .seqmap_RV.json next to the FASTA."
            )
        except RuntimeError:
            return

    def _annotation_changes_pending(self) -> bool:
        if self._map.is_dirty():
            return True
        if self._visible_transfer_prediction_count() > 0:
            return True
        if not self._map.features():
            return False
        # Pending if there are features but no JSON sidecar yet.
        from repeatviewer.backends.sequence_features import resolve_sidecar_path

        return resolve_sidecar_path(self._path) is None

    def _save_annotations(self, *, silent: bool = False) -> bool:
        if self._map.data() is None:
            return False
        committed = self._commit_visible_transfer_predictions()
        if not self._annotation_changes_pending() and committed <= 0:
            if not silent:
                self.statusMessage.emit("No annotation changes to save.")
            return True
        if not self._persist_annotations(show_message=not silent):
            return False
        if not silent:
            from repeatviewer.backends.sequence_features import sidecar_path

            QMessageBox.information(
                self,
                "Saved",
                f"Annotations saved to:\n{sidecar_path(self._path).name}\n\n"
                "To get a GFF3 / GTF / BED file, use Export in the Annotations "
                "panel (or right-click the FASTA → Export Annotation…).",
            )
        return True

    def _persist_annotations(self, *, show_message: bool = True) -> bool:
        data = self._map.data()
        if data is None:
            return False
        # Sync canvas edits into the multi-seq store before writing.
        if self._loaded_data is not None:
            self._loaded_data.sync_active_features(self._map.features())
            data = self._loaded_data
        else:
            data.sync_active_features(self._map.features())
        try:
            if data.record_ids and (
                data.is_multi_sequence or len(data.features_by_id) > 1
            ):
                lengths = dict(data.lengths_by_id)
                lengths[data.sequence_id] = data.length
                json_path, _gff_path = save_annotations_multi(
                    self._path,
                    features_by_id=data.features_by_id,
                    lengths_by_id=lengths,
                    record_ids=list(data.record_ids),
                )
            else:
                json_path, _gff_path = save_annotations(
                    self._path,
                    sequence_id=data.sequence_id,
                    sequence_length=data.length,
                    features=self._map.features(),
                )
        except OSError as exc:
            QMessageBox.critical(self, "Save failed", str(exc))
            return False
        self._map.mark_clean()
        if self._loaded_data is not None:
            self._populate_sequence_combo(self._loaded_data)
            self._sync_too_long_controls()
        self._last_saved_json = json_path
        self._last_saved_gff = _gff_path
        self._sync_save_marker()
        if show_message:
            self.statusMessage.emit(f"Saved annotations: {json_path.name}")
        self.annotationsSaved.emit(json_path)
        if _gff_path is not None:
            self.annotationsSaved.emit(_gff_path)
        return True

    def _save_annotations_from_button(self) -> None:
        """Panel Save: write the JSON (and GFF3 if enabled) and say where."""
        self._last_saved_json = None
        self._last_saved_gff = None
        if not self._save_annotations(silent=True):
            return
        json_path = getattr(self, "_last_saved_json", None)
        if json_path is None:
            self.statusMessage.emit("No annotation changes to save.")
            return
        data = self._loaded_data if self._loaded_data is not None else self._map.data()
        count = 0
        if data is not None:
            count = sum(len(v) for v in data.features_by_id.values()) or len(
                self._map.features()
            )
        noun = "annotation" if count == 1 else "annotations"
        gff = getattr(self, "_last_saved_gff", None)
        if gff is not None:
            text = f"Saved {count:,} {noun} to {json_path.name} and {gff.name}"
        else:
            text = (
                f"Saved {count:,} {noun} to {json_path.name}. "
                "To get a GFF3: Export → GFF3"
            )
        self.statusMessage.emit(text)


    def _on_add_feature_clicked(self) -> None:
        """Add the feature type chosen in Manually add feature."""
        if hasattr(self, "_rail"):
            label, feature_type = self._rail.selected_manual_feature()
            self._add_feature_of_type(feature_type, feature_label=label)
            return
        self._add_feature_of_type("", feature_label=None)

    def _refresh_transfer_sequence_lists(self) -> None:
        """Keep hidden target combo pointed at the open sequence."""
        self._refresh_transfer_blast_sources()

    def _refresh_transfer_blast_sources(self) -> None:
        """Compatibility no-op — Transfer Source is a free path (file or folder)."""
        if not hasattr(self, "_rail"):
            return
        data = self._loaded_data
        tgt = self._rail.transfer_target
        tgt.blockSignals(True)
        try:
            tgt.clear()
            if data is not None and str(data.sequence_id or ""):
                rid = str(data.sequence_id)
                tgt.addItem(rid, rid)
                tgt.setCurrentIndex(0)
        finally:
            tgt.blockSignals(False)

    def _on_manage_database_clicked(self, host: QWidget | None = None) -> None:
        # Transfer Source uses Browse; Manage Database remains for Similarity DBs.
        # Pass the modal settings dialog as host so Manage is not frozen under it.
        window = self.window()
        opener = getattr(window, "_open_manage_references", None)
        if callable(opener):
            opener(role_filter="blast", host=host)
        self.refresh_similarity_databases()

    def _on_transfer_source_picked(self, payload) -> None:
        """Fill Transfer Source path from an external pick (file / DB / genome)."""
        from repeatviewer.project import BlastDatabase, ManagedGenome, UserDatabase

        path: Path | None = None
        if isinstance(payload, BlastDatabase):
            path = Path(payload.source_fasta).expanduser().resolve()
        elif isinstance(payload, ManagedGenome):
            path = Path(payload.path).expanduser().resolve()
        elif isinstance(payload, UserDatabase):
            path = Path(payload.path).expanduser().resolve()
        else:
            try:
                path = Path(payload).expanduser().resolve()
            except (TypeError, ValueError):
                path = None
        if path is None:
            return
        self.apply_path_from_panel(path, force=True)

    def _arm_transfer_source_path_fill(self, _widget=None) -> None:
        """Mark Transfer Source as the click-to-fill target (Files/Folders panel)."""
        self._transfer_path_armed = True

    def _transfer_source_accepts_panel_fill(self) -> bool:
        if not hasattr(self, "_rail") or not hasattr(self._rail, "transfer_source_row"):
            return False
        if bool(getattr(self, "_transfer_path_armed", False)):
            return True
        line = self._rail.transfer_source_row.line_edit()
        if line.hasFocus():
            return True
        # Empty Source with Transfer open: first Files/Folders click fills it
        # without requiring a prior focus (same idea as Define inputs).
        sections = getattr(self._rail, "_annotation_sections", None) or {}
        transfer = sections.get("transfer")
        if transfer is not None and transfer.is_expanded():
            return self._rail.transfer_source_row.path() is None
        return False

    def apply_path_from_panel(self, path: Path, *, force: bool = False) -> bool:
        """Fill Transfer Source from a Files-panel click (when the field is armed)."""
        if not hasattr(self, "_rail") or not hasattr(self._rail, "transfer_source_row"):
            return False
        if not force and not self._transfer_source_accepts_panel_fill():
            return False
        try:
            resolved = Path(path).expanduser().resolve()
        except OSError:
            resolved = Path(path).expanduser()
        if not resolved.is_file() and not resolved.is_dir():
            return False
        if resolved.is_file():
            from repeatviewer.backends.annotation_transfer import _is_fasta_path

            if not _is_fasta_path(resolved):
                return False
        self._rail.transfer_source_row.set_path(resolved)
        self._transfer_path_armed = False
        self.statusMessage.emit(f"Transfer source: {resolved.name}")
        return True

    def apply_folder_from_panel(self, folder: Path, *, force: bool = False) -> bool:
        """Fill Transfer Source from a Folders-panel click (when the field is armed)."""
        if not hasattr(self, "_rail") or not hasattr(self._rail, "transfer_source_row"):
            return False
        if not force and not self._transfer_source_accepts_panel_fill():
            return False
        try:
            resolved = Path(folder).expanduser().resolve()
        except OSError:
            resolved = Path(folder).expanduser()
        if not resolved.is_dir():
            return False
        self._rail.transfer_source_row.set_path(resolved)
        self._transfer_path_armed = False
        self.statusMessage.emit(f"Transfer source: {resolved.name}/")
        return True

    def _on_transfer_browse_target(self) -> None:
        # Target is always the open sequence; browse is hidden but kept for API compat.
        return

    def _sequence_text_for_id(self, seq_id: str) -> tuple[str, str, list]:
        """Return (id, sequence, features) for a record in the loaded FASTA."""
        from Bio import SeqIO

        data = self._loaded_data
        if data is None:
            raise ValueError("No sequence loaded.")
        wanted = str(seq_id or "").strip()
        if not wanted:
            raise ValueError("Empty sequence id.")
        if wanted == str(data.sequence_id or "") and getattr(data, "sequence", ""):
            return wanted, data.sequence, list(data.features)
        # Load from disk for another record in the same file.
        for rec in SeqIO.parse(str(self._path), "fasta"):
            rid = str(rec.id)
            if rid == wanted:
                feats = list((data.features_by_id or {}).get(rid) or [])
                return rid, str(rec.seq), feats
        raise ValueError(f"Sequence “{wanted}” was not found in {self._path.name}.")

    def _on_transfer_annotations_clicked(self) -> None:
        try:
            if self._transfer_session_active():
                self._stop_transfer_session()
                return
            # Run starts the hybrid search; while Stop is shown, sliders re-filter
            # stored hits. Save commits visible predictions to the annotation file.
            self._start_transfer_worker()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Transfer", str(exc))
            self._stop_transfer_session(status="")

    def _transfer_running(self) -> bool:
        thread = self._transfer_thread
        return thread is not None and thread.isRunning()

    def _transfer_session_active(self) -> bool:
        """True while searching or showing live Transfer predictions."""
        return self._transfer_running() or self._transfer_preview_lifted is not None

    def _set_transfer_run_button(self, *, stop_mode: bool) -> None:
        btn = self._rail.transfer_run_btn
        settings_btn = getattr(self._rail, "transfer_ann_settings_btn", None)
        if stop_mode:
            btn.setText("Stop")
            btn.setFixedWidth(52)
            btn.setEnabled(True)
        else:
            btn.setText("Run")
            btn.setFixedWidth(52)
            btn.setEnabled(True)
        if settings_btn is not None:
            # Allow Setting during live preview; disable only while searching.
            settings_btn.setEnabled(not self._transfer_running())
        self._sync_main_window_stop_for_transfer(stop_mode)

    def _sync_main_window_stop_for_transfer(self, stop_mode: bool) -> None:
        """Enable the status-bar Stop while a Transfer session is active."""
        win = self.window()
        stop_btn = getattr(win, "_stop_btn", None) if win is not None else None
        if stop_btn is None:
            return
        if stop_mode:
            stop_btn.setEnabled(True)
            return
        if getattr(win, "_job_running", False):
            return
        if self._auto_annotate_running() or self._similarity_annotation_running():
            return
        stop_btn.setEnabled(False)

    def cancel_transfer_annotation(self) -> bool:
        """Stop an active Transfer search/preview (rail or status-bar Stop)."""
        if not self._transfer_session_active():
            return False
        self._stop_transfer_session()
        return True

    def cancel_running_transfer_annotation(self) -> bool:
        """Stop only while the Transfer search worker is running (e.g. file click)."""
        if not self._transfer_running():
            return False
        self._stop_transfer_session(status="Transfer stopped.")
        return True

    def _stop_transfer_session(self, *, status: str | None = None) -> None:
        was_running = self._transfer_running()
        if was_running:
            self._transfer_cancel.request()
            self._rail.transfer_run_btn.setText("Stopping…")
            self._rail.transfer_run_btn.setFixedWidth(80)
            self._rail.transfer_run_btn.setEnabled(False)
            win = self._main_window_for_job_progress()
            if win is not None and hasattr(win, "_task_progress") and win._task_progress:
                win._task_progress.update(
                    getattr(win, "_job_work_steps", 0), "Stopping…"
                )
        self._clear_transfer_preview_state()
        if not was_running:
            # Preview-only stop: hide the held "Running Transfer annotation…" bar.
            self._set_transfer_run_button(stop_mode=False)
            self._end_window_job_progress()
        if status is None:
            status = "Transfer stopped."
        if status:
            self.statusMessage.emit(status)

    def _on_transfer_source_changed(self) -> None:
        # New source ends the current session; click Run again to search.
        if self._transfer_session_active():
            self._stop_transfer_session(
                status="Transfer cleared (source changed). Click Run to search again."
            )

    def _on_transfer_reverse_changed(self) -> None:
        if self._transfer_session_active():
            self._stop_transfer_session(
                status=(
                    "Transfer cleared (reverse-strand option changed). "
                    "Click Run to search again."
                )
            )

    def _clear_transfer_preview_state(self) -> None:
        self._transfer_preview_lifted = None
        self._transfer_preview_source_path = None
        self._map.clear_transfer_preview()

    def _transfer_source_path(self) -> Path | None:
        """User-chosen FASTA file or folder for Transfer Source."""
        row = getattr(self._rail, "transfer_source_row", None)
        if row is not None:
            path = row.path()
            if path is None:
                return None
            try:
                resolved = path.expanduser().resolve()
            except OSError:
                resolved = path.expanduser()
            if resolved.is_file() or resolved.is_dir():
                return resolved
            return None
        # Legacy combo fallback.
        source_path_raw = str(self._rail.transfer_source.currentData() or "").strip()
        if not source_path_raw:
            return None
        path = Path(source_path_raw).expanduser().resolve()
        return path if path.is_file() or path.is_dir() else None

    def _visible_transfer_predictions(self) -> list[SequenceFeature]:
        """Features that pass the current Min identity / Min coverage sliders."""
        from repeatviewer.backends.annotation_transfer import LiftedAnnotation

        lifted = self._transfer_preview_lifted
        if not lifted:
            return []
        min_id = self._rail.transfer_min_identity_value()
        min_cov = self._rail.transfer_min_coverage_value()
        out: list[SequenceFeature] = []
        for item in lifted:
            if not isinstance(item, LiftedAnnotation) or item.feature is None:
                continue
            if item.identity < min_id or item.coverage < min_cov:
                continue
            out.append(_clone_features([item.feature])[0])
        return out

    def _visible_transfer_prediction_count(self) -> int:
        return len(self._visible_transfer_predictions())

    def _commit_visible_transfer_predictions(self) -> int:
        """Merge live Transfer predictions into the map (disk write via Save)."""
        new_feats = self._visible_transfer_predictions()
        if not new_feats:
            if self._transfer_preview_lifted is not None:
                self._clear_transfer_preview_state()
                self._set_transfer_run_button(stop_mode=False)
                self._end_window_job_progress()
            return 0
        data = self._map.data()
        if data is None:
            return 0
        self._push_feature_undo_snapshot()
        merged = list(self._map.features()) + new_feats
        self._map.set_features(merged)
        if self._loaded_data is not None:
            self._loaded_data.sync_active_features(merged)
        count = len(new_feats)
        self._clear_transfer_preview_state()
        self._set_transfer_run_button(stop_mode=False)
        self._end_window_job_progress()
        self.statusMessage.emit(
            f"Kept {count} transferred feature(s). Saving…"
        )
        return count

    def _apply_transfer_preview_filter(self) -> None:
        """Re-filter stored Run results when Min identity / Min coverage change."""
        from repeatviewer.backends.annotation_transfer import LiftedAnnotation

        lifted = self._transfer_preview_lifted
        if not lifted:
            self._map.clear_transfer_preview()
            return
        min_id = self._rail.transfer_min_identity_value()
        min_cov = self._rail.transfer_min_coverage_value()
        overlay: list[tuple[SequenceFeature, int]] = []
        for item in lifted:
            if not isinstance(item, LiftedAnnotation) or item.feature is None:
                continue
            alpha = _transfer_preview_alpha(
                item.identity, item.coverage, min_id, min_cov
            )
            if alpha is None:
                continue
            feat = _clone_features([item.feature])[0]
            feat.source = _TRANSFER_PREVIEW_SOURCE
            overlay.append((feat, alpha))
        self._map.set_transfer_preview(overlay)
        self._set_transfer_run_button(stop_mode=True)
        # Slider changes re-run the status-bar progress (0 → full) while the
        # session stays active with "Running Transfer annotation…".
        if not self._transfer_running():
            self._replay_transfer_session_progress()
        kept = len(overlay)
        total = sum(1 for i in lifted if i.feature is not None)
        self.statusMessage.emit(
            f"Transfer preview: {kept} of {total} feature(s) — "
            "adjust sliders, then Save to keep."
        )

    def _start_transfer_worker(self) -> None:
        from repeatviewer.backends.annotation_transfer import (
            AnnotationTransferError,
            AnnotationTransferOptions,
            AnnotationTransferResult,
            discover_annotated_transfer_sources,
            load_transfer_strategy,
            load_use_system_mafft_from_settings,
            transfer_annotations,
        )
        from repeatviewer.backends.auto_annotation_settings import (
            load_annotation_run_options,
        )
        from repeatviewer.backends.sequence_features import load_sequence_map_data
        from repeatviewer.backends.similarity_annotation import (
            load_use_system_blast_from_settings,
        )
        from repeatviewer.widgets.process_monitor_dialog import (
            attach_annotation_worker_process_monitor,
        )

        data = self._loaded_data
        if data is None or not getattr(data, "sequence", ""):
            QMessageBox.information(self, "Transfer", "Load a sequence first.")
            return
        source_path = self._transfer_source_path()
        if source_path is None:
            QMessageBox.information(
                self,
                "Transfer",
                "Choose a source FASTA file or a folder of annotated FASTAs.",
            )
            return

        sources = discover_annotated_transfer_sources(source_path)
        if not sources:
            kind = "folder" if source_path.is_dir() else "file"
            QMessageBox.warning(
                self,
                "Transfer",
                f"No annotated sequences found in this {kind}.\n\n"
                "Only FASTA files with a RepeatViewer annotation JSON "
                "(*.seqmap_RV.json) and at least one feature are used.\n\n"
                f"Source:\n{source_path}",
            )
            return

        tgt_id = str(data.sequence_id or "")
        tgt_seq = str(data.sequence or "")
        if not tgt_id or not tgt_seq:
            QMessageBox.information(self, "Transfer", "Load a sequence first.")
            return

        source_label = source_path.name
        if source_path.is_dir():
            source_label = f"{source_path.name}/ ({len(sources)} annotated)"
        elif len(sources) > 1:
            source_label = f"{source_path.stem} ({len(sources)} records)"
        ann_source = f"Transfer ({source_label})"

        threads, show_process, close_process = load_annotation_run_options()
        # Thresholds applied live via sliders after Run; search keeps all hits.
        # BLAST/lift intermediates stay in tempfile dirs (never next to source).
        opts = AnnotationTransferOptions(
            min_identity=0.0,
            min_coverage=0.0,
            threads=max(1, int(threads)),
            allow_reverse=bool(self._rail.transfer_allow_reverse.isChecked()),
            strategy=load_transfer_strategy(),
            annotation_source=ann_source,
            use_system_blast=load_use_system_blast_from_settings(),
            use_system_mafft=load_use_system_mafft_from_settings(),
        )

        self._clear_transfer_preview_state()
        self._transfer_cancel = JobCancellation()
        self._set_transfer_run_button(stop_mode=True)
        self._log_annotation_job(
            "JOB START",
            "Transfer annotation",
            detail=f"sources={len(sources)} from={source_label}",
        )
        self.statusMessage.emit(
            f"Running Transfer annotation from {len(sources)} source(s)…"
        )
        self._begin_window_job_progress("Transfer annotation", len(sources))

        job_sources = list(sources)
        job_tgt_id = tgt_id
        job_tgt_seq = tgt_seq
        job_target_path = self._path
        job_opts = opts

        def _job(*, progress=None):
            from concurrent.futures import ThreadPoolExecutor, as_completed
            import threading

            from repeatviewer.backends.annotation_transfer import (
                transfer_parallel_workers,
            )
            from repeatviewer.backends.job_cancel import JobCancelledError

            # Prepare transferable pairs on this worker thread first.
            work: list[tuple] = []
            skipped: list[tuple[str, str]] = []
            for entry in job_sources:
                if self._transfer_cancel.requested:
                    raise JobCancelledError("Annotation transfer cancelled.")
                label = str(entry.label or entry.sequence_id or entry.fasta_path.stem)
                try:
                    same_file = entry.fasta_path.resolve() == job_target_path.resolve()
                except OSError:
                    same_file = False
                try:
                    source_data = load_sequence_map_data(
                        entry.fasta_path, active_id=entry.sequence_id
                    )
                except Exception as exc:  # noqa: BLE001 — keep batch going
                    skipped.append((label, str(exc) or exc.__class__.__name__))
                    continue
                src_id = str(source_data.sequence_id or entry.sequence_id)
                src_seq = str(source_data.sequence or "")
                features = list(source_data.features)
                if not src_seq or not features:
                    continue
                if same_file and src_id == job_tgt_id:
                    continue
                work.append((entry, label, src_id, src_seq, features))

            total = max(1, len(work) + len(skipped))
            parallel = transfer_parallel_workers(len(work))
            if progress is not None:
                progress(f"0/{total} starting")
                progress(
                    f"Threads={job_opts.threads} per pair; "
                    f"up to {parallel} pair(s) in parallel "
                    f"(Preferences → Max. parallel jobs); "
                    f"allow_reverse={job_opts.allow_reverse}"
                )

            all_lifted = []
            used = 0
            completed = len(skipped)
            progress_lock = threading.Lock()

            def _run_one(payload):
                entry, label, src_id, src_seq, features = payload
                per_opts = AnnotationTransferOptions(
                    min_identity=job_opts.min_identity,
                    min_coverage=job_opts.min_coverage,
                    threads=job_opts.threads,
                    allow_reverse=job_opts.allow_reverse,
                    evalue=job_opts.evalue,
                    strategy=job_opts.strategy,
                    annotation_source=f"Transfer ({entry.label or src_id})",
                    use_system_blast=job_opts.use_system_blast,
                    use_system_mafft=job_opts.use_system_mafft,
                )
                try:
                    result = transfer_annotations(
                        source_id=src_id,
                        source_sequence=src_seq,
                        target_id=job_tgt_id,
                        target_sequence=job_tgt_seq,
                        features=features,
                        opts=per_opts,
                        cancel=self._transfer_cancel,
                        source_fasta=entry.fasta_path,
                        progress=None if parallel > 1 else progress,
                    )
                    return label, result, None
                except JobCancelledError:
                    raise
                except Exception as exc:  # noqa: BLE001 — skip this source only
                    return label, None, str(exc) or exc.__class__.__name__

            with ThreadPoolExecutor(max_workers=max(1, parallel)) as pool:
                futures = [pool.submit(_run_one, item) for item in work]
                for future in as_completed(futures):
                    if self._transfer_cancel.requested:
                        for pending in futures:
                            pending.cancel()
                        raise JobCancelledError("Annotation transfer cancelled.")
                    try:
                        label, result, error = future.result()
                    except JobCancelledError:
                        for pending in futures:
                            pending.cancel()
                        raise
                    with progress_lock:
                        completed += 1
                        done = completed
                    if error:
                        skipped.append((label, error))
                        if progress is not None:
                            short = error if len(error) <= 80 else error[:77] + "…"
                            progress(f"{done}/{total} skipped {label}: {short}")
                        continue
                    assert result is not None
                    all_lifted.extend(result.lifted)
                    used += 1
                    if progress is not None:
                        progress(f"{done}/{total} {label}")

            if used <= 0 and not skipped:
                raise AnnotationTransferError(
                    "No transferable annotated sequences were found.\n"
                    "Source and target may be the same sequence, or sources "
                    "had no usable features."
                )
            if used <= 0 and skipped:
                raise AnnotationTransferError(
                    "Transfer could not place features from any source.\n\n"
                    + "\n".join(
                        f"• {lab}: {reason}" for lab, reason in skipped[:12]
                    )
                    + (
                        f"\n• …and {len(skipped) - 12} more"
                        if len(skipped) > 12
                        else ""
                    )
                )
            msg = (
                f"Transferred from {used} source sequence(s); "
                f"{sum(1 for item in all_lifted if item.feature is not None)} "
                f"placement(s) available for filtering."
            )
            if skipped:
                msg += f" Skipped {len(skipped)} source(s)."
            return (
                AnnotationTransferResult(
                    lifted=all_lifted,
                    message=msg,
                    skipped_sources=tuple(skipped),
                ),
                job_tgt_id,
                job_tgt_seq,
                None,
            )

        worker = Worker(_job, progress=True)
        attach_annotation_worker_process_monitor(
            self,
            worker,
            title="Transfer annotation",
            enabled=bool(show_process),
            close_when_finished=bool(close_process),
            tool_key="transfer_annotation",
        )
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.signals.progress.connect(self._on_transfer_annotation_progress)
        worker.signals.finished.connect(self._on_transfer_finished)
        worker.signals.error.connect(self._on_transfer_error)
        worker.signals.finished.connect(thread.quit)
        worker.signals.error.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._clear_transfer_annotation_thread)
        self._transfer_worker = worker
        self._transfer_thread = thread
        thread.start()
        settings_btn = getattr(self._rail, "transfer_ann_settings_btn", None)
        if settings_btn is not None:
            settings_btn.setEnabled(False)

    @Slot(str)
    def _on_transfer_annotation_progress(self, message: str) -> None:
        parsed = self._parse_annotation_step_progress(message)
        if parsed is None:
            return
        completed, total, detail = parsed
        self._update_window_job_progress(
            completed,
            "starting" if detail.lower() == "starting" else detail,
        )

    @Slot()
    def _clear_transfer_annotation_thread(self) -> None:
        self._transfer_thread = None
        self._transfer_worker = None
        # Keep status-bar Stop + full "Running Transfer annotation…" while live
        # predictions remain for slider filter / Save.
        if self._transfer_preview_lifted is not None:
            self._hold_transfer_session_progress()
            self._set_transfer_run_button(stop_mode=True)
        else:
            self._end_window_job_progress()
            self._set_transfer_run_button(stop_mode=False)

    @Slot(object)
    def _on_transfer_finished(self, payload: object) -> None:
        cancelled = bool(getattr(self._transfer_cancel, "requested", False))
        if cancelled:
            self._clear_transfer_preview_state()
            self._log_annotation_job(
                "JOB FINISH",
                "Transfer annotation",
                detail="stopped by user",
            )
            self.statusMessage.emit("Annotation transfer stopped.")
            return
        try:
            result, _tid, _tseq, _ext = payload  # type: ignore[misc]
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Transfer", f"Unexpected result: {exc}")
            self._clear_transfer_preview_state()
            return
        self._transfer_preview_lifted = list(result.lifted)
        self._transfer_preview_source_path = self._transfer_source_path()
        self._apply_transfer_preview_filter()
        self._log_annotation_job(
            "JOB FINISH",
            "Transfer annotation",
            detail=getattr(result, "message", ""),
        )
        skipped = list(getattr(result, "skipped_sources", ()) or ())
        if skipped:
            from repeatviewer.backends.annotation_transfer import (
                transfer_skipped_reminder_parts,
            )

            headline, body = transfer_skipped_reminder_parts(skipped)
            box = QMessageBox(self)
            box.setIcon(QMessageBox.Icon.Warning)
            box.setWindowTitle("Transfer — some sources skipped")
            box.setText(headline)
            box.setInformativeText(body)
            box.setStandardButtons(QMessageBox.StandardButton.Ok)
            box.exec()
            self.statusMessage.emit(
                f"Transfer finished with {len(skipped)} source(s) skipped."
            )

    @Slot(str)
    def _on_transfer_error(self, message: str) -> None:
        self._clear_transfer_preview_state()
        lower = message.lower()
        if "cancel" in lower or "stopped by user" in lower:
            self._log_annotation_job(
                "JOB FINISH",
                "Transfer annotation",
                detail="stopped by user",
            )
            self.statusMessage.emit("Annotation transfer stopped.")
            return
        self._log_annotation_job(
            "JOB ERROR",
            "Transfer annotation",
            detail=message,
            level=logging.ERROR,
        )
        QMessageBox.warning(self, "Transfer", message)
        self.statusMessage.emit(f"Transfer failed: {message}")

    def _persist_transferred_features(
        self,
        new_feats: list[SequenceFeature],
        target_id: str,
        target_seq: str,
        *,
        external_path: Path | None = None,
    ) -> None:
        from repeatviewer.backends.sequence_features import (
            load_features_by_id_from_json,
            save_annotations,
            save_annotations_multi,
        )

        if external_path is not None:
            path = Path(external_path)
            existing = load_features_by_id_from_json(path, [target_id]).get(target_id) or []
            merged = list(existing) + list(new_feats)
            save_annotations(
                path,
                sequence_id=target_id,
                sequence_length=len(target_seq),
                features=merged,
            )
            return

        data = self._loaded_data
        if data is None:
            return
        by_id = dict(data.features_by_id or {})
        current = list(by_id.get(target_id) or [])
        current.extend(new_feats)
        by_id[target_id] = current
        lengths = {rid: 0 for rid in data.record_ids}
        if target_id == data.sequence_id:
            lengths[target_id] = data.length
        else:
            lengths[target_id] = len(target_seq)
        save_annotations_multi(
            self._path,
            features_by_id=by_id,
            lengths_by_id=lengths,
            record_ids=list(data.record_ids),
        )
        data.features_by_id = by_id
        if target_id == data.sequence_id:
            self._push_feature_undo_snapshot()
            self._map.set_features(list(by_id[target_id]))
        self._populate_sequence_combo(data)

    def _apply_transfer_result(self, result, target_id: str, target_seq: str, external_path) -> None:
        from repeatviewer.backends.annotation_transfer import AnnotationTransferResult

        if not isinstance(result, AnnotationTransferResult):
            QMessageBox.warning(self, "Transfer", "Unexpected transfer result.")
            return
        new_feats = result.features
        if not new_feats:
            self.statusMessage.emit(
                result.message
                or "No features were transferred at the current thresholds."
            )
            return
        self._persist_transferred_features(
            new_feats, target_id, target_seq, external_path=external_path
        )
        self._clear_transfer_preview_state()
        self.statusMessage.emit(
            f"Transferred {len(new_feats)} feature(s). Click Save to write to disk."
        )

    def _add_feature_of_type(self, feature_type: str, *, feature_label: str | None = None) -> None:
        data = self._map.data()
        if data is None:
            return
        sel = self._map.selection_range()
        if sel is None:
            QMessageBox.information(
                self,
                "Add feature",
                "Drag on the sequence row to select a region, then choose a feature type.",
            )
            return
        start, end = sel
        label = feature_label or feature_type
        existing = sum(1 for feat in data.features if feat.feature_type == feature_type)
        if feature_type == "ORF":
            default_name = f"ORF_{existing + 1}"
        else:
            default_name = f"{feature_type.replace(' ', '_')}{existing + 1}"
        dialog = AddFeatureDialog(
            start=start,
            end=end,
            seq_length=data.length,
            parent=self,
            feature_type=feature_type,
            feature_label=label,
            default_name=default_name,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        feat = dialog.feature()
        if feat is None:
            QMessageBox.warning(self, "Add feature", "Enter a feature name.")
            return
        self._push_feature_undo_snapshot()
        features = list(data.features)
        features.append(feat)
        self._map.set_features(features)
        # Show the new type under Display → Displayed annotations (e.g. LTR).
        kind = str(feat.feature_type or "").strip()
        if kind:
            self._reveal_annotation_kinds_after_run((kind,))
        self.statusMessage.emit("Feature added. Click Save to write changes to disk.")

    def _push_feature_undo_snapshot(self) -> None:
        self._feature_undo_stack.append(_clone_features(self._map.features()))
        if len(self._feature_undo_stack) > _FEATURE_UNDO_LIMIT:
            self._feature_undo_stack.pop(0)

    def _delete_feature(self, *, confirm_empty: bool = True) -> None:
        self._map._sync_feature_selection_from_annot_track()
        if not self._map._selected_feature_indices:
            if confirm_empty:
                QMessageBox.information(
                    self,
                    "Delete",
                    "Select one or more features in the track below the sequence, "
                    "then press Delete or Backspace, or click Delete.\n\n"
                    "Drag an area on the feature track to select multiple features.",
                )
            return
        self._push_feature_undo_snapshot()
        removed = self._map.delete_selected_features()
        if removed <= 0:
            self._feature_undo_stack.pop()
            return
        self._persist_annotations()
        self.statusMessage.emit(f"Deleted {removed} feature(s).")

    def _modify_feature(self) -> None:
        data = self._map.data()
        if data is None:
            return
        self._map._sync_feature_selection_from_annot_track()
        indices = sorted(self._map.selected_feature_indices())
        if len(indices) != 1:
            QMessageBox.information(
                self,
                "Modify feature",
                "Select a single feature to modify.",
            )
            return
        index = indices[0]
        features = list(data.features)
        if not (0 <= index < len(features)):
            return
        current = features[index]
        parent_candidates = [
            feat.gff_id()
            for i, feat in enumerate(features)
            if i != index and feat.gff_id()
        ]
        dialog = ModifyFeatureDialog(
            feature=current,
            seq_length=data.length,
            parent_candidates=parent_candidates,
            parent=self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        edited = dialog.feature()
        if edited is None:
            QMessageBox.warning(self, "Modify feature", "Enter a feature name and type.")
            return
        if edited.parent and edited.parent == edited.gff_id():
            QMessageBox.warning(
                self,
                "Modify feature",
                "A feature cannot be its own GFF3 Parent.",
            )
            return
        self._push_feature_undo_snapshot()
        old_id = current.gff_id()
        features[index] = edited
        features = remount_parent_ids(features, old_id=old_id, new_id=edited.gff_id())
        self._map.set_features(features)
        self._map.set_selected_features({index})
        self._persist_annotations()
        self.statusMessage.emit(f"Modified feature “{edited.name}”.")

    def _change_selected_feature_colors(self) -> None:
        data = self._map.data()
        if data is None:
            return
        self._map._sync_feature_selection_from_annot_track()
        indices = sorted(self._map.selected_feature_indices())
        if len(indices) < 2:
            QMessageBox.information(
                self,
                "Change color",
                "Select two or more features to change their color together.",
            )
            return
        features = list(data.features)
        first = features[indices[0]]
        initial = QColor(first.color or feature_type_color(first.feature_type))
        chosen = QColorDialog.getColor(initial, self, "Change color")
        if not chosen.isValid():
            return
        color_hex = chosen.name(QColor.NameFormat.HexRgb)
        selected = set(indices)
        self._push_feature_undo_snapshot()
        updated = [
            _feature_with_color(feat, color_hex) if idx in selected else feat
            for idx, feat in enumerate(features)
        ]
        self._map.set_features(updated)
        self._map.set_selected_features(selected)
        self._persist_annotations()
        self.statusMessage.emit(f"Changed color for {len(selected)} feature(s).")

    def _apply_annotation_track_move(self, features: object) -> None:
        updated = list(features or [])
        if not updated:
            return
        selected = set(self._map.selected_feature_indices())
        self._push_feature_undo_snapshot()
        self._map.set_features(updated)
        self._map.set_selected_features(selected)
        self._persist_annotations()
        self.statusMessage.emit("Moved annotation track(s).")

    def _apply_annotation_feature_resize(self, features: object) -> None:
        updated = list(features or [])
        if not updated:
            return
        selected = set(self._map.selected_feature_indices())
        self._push_feature_undo_snapshot()
        self._map.set_features(updated)
        self._map.set_selected_features(selected)
        self._persist_annotations()
        picked = [
            updated[i] for i in sorted(selected) if 0 <= i < len(updated)
        ]
        if len(picked) == 1:
            feat = picked[0]
            self.statusMessage.emit(
                f"Resized “{feat.name or feat.feature_type}” "
                f"to {feat.start}–{feat.end}."
            )
        else:
            self.statusMessage.emit("Resized annotation feature(s).")

    def _selected_features(self) -> list[SequenceFeature]:
        self._map._sync_feature_selection_from_annot_track()
        features = self._map.features()
        indices = sorted(self._map.selected_feature_indices())
        return [features[index] for index in indices if 0 <= index < len(features)]

    def _copy_selected_feature_sequences(self) -> None:
        data = self._loaded_data
        if data is None or not getattr(data, "sequence", ""):
            QMessageBox.warning(
                self, "Copy sequence", "The sequence is not loaded yet."
            )
            return
        selected = self._selected_features()
        if not selected:
            QMessageBox.information(
                self, "Copy sequence", "Select a feature first."
            )
            return
        blocks: list[str] = []
        for feat in selected:
            fragment = _feature_sequence_fragment(data.sequence, feat)
            if not fragment:
                continue
            strand = "-" if not feat.is_forward() else "+"
            header = (
                f">{feat.name or feat.feature_type} "
                f"{data.sequence_id}:{feat.start}-{feat.end}({strand})"
            )
            if len(selected) == 1:
                blocks.append(fragment)
            else:
                lines = [header]
                for offset in range(0, len(fragment), 60):
                    lines.append(fragment[offset : offset + 60])
                blocks.append("\n".join(lines))
        if not blocks:
            QMessageBox.warning(self, "Copy sequence", "The selection is empty.")
            return
        QApplication.clipboard().setText("\n\n".join(blocks))
        total_bp = sum(
            len(_feature_sequence_fragment(data.sequence, feat)) for feat in selected
        )
        self.statusMessage.emit(
            f"Copied {len(selected)} feature sequence(s) ({total_bp:,} bp) to the clipboard"
        )

    def _extract_selected_feature_sequences(self) -> None:
        data = self._loaded_data
        if data is None or not getattr(data, "sequence", ""):
            QMessageBox.warning(
                self, "Extract Sequence", "The sequence is not loaded yet."
            )
            return
        selected = self._selected_features()
        if not selected:
            QMessageBox.information(
                self, "Extract Sequence", "Select a feature first."
            )
            return
        sequence = data.sequence
        for feat in selected:
            dialog = _ExtractSequenceDialog(
                sequence_id=data.sequence_id,
                sequence_length=len(sequence),
                start=int(feat.start),
                end=int(feat.end),
                default_name=self._next_feature_extraction_name(
                    data.sequence_id, feat
                ),
                reverse=not feat.is_forward(),
                parent=self,
            )
            if dialog.exec() != QDialog.DialogCode.Accepted:
                continue
            chosen_start, chosen_end, name, reverse = dialog.values()
            self._write_extracted_fasta(
                sequence=sequence,
                sequence_id=data.sequence_id,
                start=chosen_start,
                end=chosen_end,
                name=name,
                reverse=reverse,
            )

    def _restore_features(self) -> None:
        if not self._feature_undo_stack:
            QMessageBox.information(
                self,
                "Restore",
                "No feature changes to restore for this sequence.",
            )
            return
        features = self._feature_undo_stack.pop()
        self._map.set_features(features)
        self._map.clear_feature_selection()
        self._persist_annotations()
        self.statusMessage.emit("Restored previous feature state.")

    @Slot(str)
    def _set_annotation_job_status(self, text: str) -> None:
        msg = text.strip()
        if hasattr(self, "_rail") and hasattr(self._rail, "set_annotation_job_status"):
            self._rail.set_annotation_job_status(msg)

    def _clear_annotation_job_status(self) -> None:
        if hasattr(self, "_rail") and hasattr(self._rail, "set_annotation_job_status"):
            self._rail.set_annotation_job_status("")

    def _show_wide_information(self, title: str, message: str) -> None:
        """Multi-line summary dialog with text that fills the content width.

        Native macOS QMessageBox ignores stylesheet width and leaves a large
        empty gap beside short lines — use a plain QDialog instead.
        """
        dialog_w = 440
        margin = 16
        icon_w = 40
        gap = 12
        text_w = dialog_w - 2 * margin - icon_w - gap

        dlg = QDialog(self)
        dlg.setWindowTitle(title)
        dlg.setModal(True)
        dlg.setFixedWidth(dialog_w)

        root = QVBoxLayout(dlg)
        root.setContentsMargins(margin, margin, margin, 12)
        root.setSpacing(14)

        row = QHBoxLayout()
        row.setSpacing(gap)
        icon_label = QLabel(dlg)
        icon = dlg.style().standardIcon(QStyle.StandardPixmap.SP_MessageBoxInformation)
        icon_label.setPixmap(icon.pixmap(icon_w, icon_w))
        icon_label.setAlignment(Qt.AlignmentFlag.AlignTop)
        icon_label.setFixedWidth(icon_w)
        row.addWidget(icon_label, 0)

        text = QLabel(message, dlg)
        text.setWordWrap(True)
        text.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        text.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        text.setFixedWidth(text_w)
        text.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Preferred)
        row.addWidget(text, 0)
        root.addLayout(row)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok, dlg)
        buttons.accepted.connect(dlg.accept)
        root.addWidget(buttons, 0, Qt.AlignmentFlag.AlignRight)

        dlg.adjustSize()
        dlg.exec()

    def _set_auto_run_button_running(self, running: bool) -> None:
        btn = self._auto_ann_run_btn
        btn.setFixedWidth(52)
        btn.setText("Stop" if running else "Run")

    def _set_similarity_run_button_running(self, running: bool) -> None:
        btn = self._similarity_ann_run_btn
        btn.setFixedWidth(52)
        btn.setText("Stop" if running else "Run")

    def _on_auto_run_clicked(self) -> None:
        if self._auto_annotate_running():
            self._stop_auto_annotation()
            return
        self._auto_annotate()

    def _on_similarity_run_clicked(self) -> None:
        if self._similarity_annotation_running():
            self._stop_similarity_annotation()
            return
        self._similarity_based_annotate()

    def _log_annotation_job(
        self,
        event: str,
        name: str,
        *,
        detail: str = "",
        level: int = logging.INFO,
    ) -> None:
        operation_logger.log(
            level,
            "%s name=%s fasta=%s detail=%s",
            event,
            name,
            str(self._path) if self._path is not None else "",
            detail,
        )

    def _job_progress_hosts(self) -> list:
        from repeatviewer.widgets.task_progress import iter_job_progress_hosts

        return iter_job_progress_hosts(self)

    def _main_window_for_job_progress(self):
        """Prefer MainWindow (tool lock + status bar); else the local window."""
        hosts = self._job_progress_hosts()
        for win in hosts:
            if hasattr(win, "_set_job_actions_enabled"):
                return win
        return hosts[0] if hosts else None

    def _begin_window_job_progress(self, title: str, total_steps: int) -> None:
        steps = max(1, int(total_steps))
        for win in self._job_progress_hosts():
            win._show_job_progress(title, steps)
            if hasattr(win, "_set_job_actions_enabled"):
                win._set_job_actions_enabled(False)

    def _update_window_job_progress(self, completed: int, detail: str = "") -> None:
        """Update the status-bar progress bar (``completed`` = finished steps)."""
        value = max(0, int(completed))
        text = str(detail or "")
        for win in self._job_progress_hosts():
            if hasattr(win, "_on_job_progress"):
                win._on_job_progress(value, text)

    def _hold_transfer_session_progress(self) -> None:
        """Keep a full bar labeled Running Transfer annotation… after search."""
        win = self._main_window_for_job_progress()
        if win is None:
            return
        title = "Transfer annotation"
        if not getattr(win, "_job_running", False):
            self._begin_window_job_progress(title, 1)
        bar = getattr(win, "_task_progress", None)
        if bar is None:
            return
        # Prefer the canonical Running title; empty detail keeps it uncluttered.
        running = (
            win._running_job_title(title)
            if hasattr(win, "_running_job_title")
            else f"Running {title}…"
        )
        total = max(1, int(getattr(bar, "_total", 1)))
        bar._job_title = running
        bar.update(total, "")
        stop_btn = getattr(win, "_stop_btn", None)
        if stop_btn is not None:
            stop_btn.setEnabled(True)

    def _replay_transfer_session_progress(self) -> None:
        """Reset then refill the status bar when identity/coverage change."""
        win = self._main_window_for_job_progress()
        if win is None:
            return
        title = "Transfer annotation"
        self._begin_window_job_progress(title, 1)
        bar = getattr(win, "_task_progress", None)
        if bar is None:
            return
        # Start empty so the redo is visible, then fill on the next tick.
        bar.update(0, "")
        QTimer.singleShot(40, self._hold_transfer_session_progress)

    def _end_window_job_progress(self) -> None:
        # Fill the reserved final slot, then hide (same as BLAST after UI ready).
        for win in self._job_progress_hosts():
            if getattr(win, "_job_running", False) and hasattr(win, "_task_progress"):
                bar = win._task_progress
                if bar is not None:
                    total = max(1, int(getattr(bar, "_total", 1)))
                    bar.update(total, "Done")
            if hasattr(win, "_hide_job_progress"):
                win._hide_job_progress()
            # Unlock BLAST/TE-Aid unless a MainWindow QThread job is still active.
            if hasattr(win, "_set_job_actions_enabled"):
                worker_busy = (
                    getattr(win, "_worker_thread", None) is not None
                    or getattr(win, "_genome_map_thread", None) is not None
                    or getattr(win, "_benchmark_thread", None) is not None
                )
                if not worker_busy:
                    win._set_job_actions_enabled(True)

    @staticmethod
    def _parse_annotation_step_progress(message: str) -> tuple[int, int, str] | None:
        """Parse ``N/M label`` progress from feature detection / similarity."""
        text = str(message or "").strip()
        match = re.match(r"^(\d+)\s*/\s*(\d+)\s+(.+)$", text)
        if not match:
            return None
        return int(match.group(1)), int(match.group(2)), match.group(3).strip()

    @Slot(str)
    def _on_auto_annotation_progress(self, message: str) -> None:
        parsed = self._parse_annotation_step_progress(message)
        if parsed is None:
            self._set_annotation_job_status(message)
            return
        completed, total, detail = parsed
        if detail.lower() == "starting":
            status = f"Feature detection 0/{total}…"
        else:
            status = f"Feature detection {completed}/{total}: {detail}…"
        self._set_annotation_job_status(status)
        self._update_window_job_progress(completed, detail)

    @Slot(str)
    def _on_similarity_annotation_progress(self, message: str) -> None:
        parsed = self._parse_annotation_step_progress(message)
        if parsed is None:
            self._set_annotation_job_status(message)
            return
        completed, total, detail = parsed
        if detail.lower() == "starting":
            status = f"Similarity 0/{total}…"
        else:
            status = f"Similarity {completed}/{total}: {detail}…"
        self._set_annotation_job_status(status)
        self._update_window_job_progress(completed, detail)

    def _stop_auto_annotation(self) -> None:
        if not self._auto_annotate_running():
            return
        self._log_annotation_job("JOB STOP", "Auto annotation")
        self._auto_ann_cancel.request()
        self._set_annotation_job_status("Stopping auto annotation…")
        win = self._main_window_for_job_progress()
        if win is not None and hasattr(win, "_task_progress") and win._task_progress:
            win._task_progress.update(
                getattr(win, "_job_work_steps", 0), "Stopping…"
            )

    def _stop_similarity_annotation(self) -> None:
        if not self._similarity_annotation_running():
            return
        self._log_annotation_job("JOB STOP", "Similarity annotation")
        self._sim_ann_cancel.request()
        self._set_annotation_job_status("Stopping similarity annotation…")
        win = self._main_window_for_job_progress()
        if win is not None and hasattr(win, "_task_progress") and win._task_progress:
            win._task_progress.update(
                getattr(win, "_job_work_steps", 0), "Stopping…"
            )

    def _auto_annotate_running(self) -> bool:
        thread = self._auto_ann_thread
        return thread is not None and thread.isRunning()

    def _set_auto_annotation_enabled(self, enabled: bool) -> None:
        self._auto_ann_settings_btn.setEnabled(enabled)
        self._auto_ann_all_check.setEnabled(enabled)
        for check in self._auto_ann_checks.values():
            check.setEnabled(enabled)
        self._auto_ann_run_btn.setEnabled(True)
        self._set_auto_run_button_running(not enabled)

    def _open_auto_annotation_settings(self) -> None:
        self._auto_ann_settings = load_auto_annotation_settings()
        dialog = AutoAnnotationSettingsDialog(self._auto_ann_settings, parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._auto_ann_settings = dialog.settings()
        save_auto_annotation_settings(self._auto_ann_settings)
        save_tsd_search_options(dialog.tsd_options())
        self._broadcast_shared_annotation_settings_changed()
        self.statusMessage.emit("Auto annotation settings saved.")

    def _selected_auto_annotation_modes(self) -> list[AutoAnnotationMode]:
        return [mode for mode, check in self._auto_ann_checks.items() if check.isChecked()]

    def _open_similarity_annotation_settings(self) -> None:
        self._similarity_ann_settings = load_similarity_annotation_settings()
        dialog = SimilarityAnnotationSettingsDialog(self._similarity_ann_settings, parent=self)
        dialog.manageDatabaseClicked.connect(
            lambda: self._on_manage_database_clicked(host=dialog)
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            self.refresh_similarity_databases()
            return
        self._similarity_ann_settings = dialog.settings()
        save_similarity_annotation_settings(self._similarity_ann_settings)
        self._map.set_show_similarity_database_name(
            self._similarity_ann_settings.show_database_name_on_track
        )
        self._broadcast_shared_annotation_settings_changed()
        self.refresh_similarity_databases()
        self.statusMessage.emit("Similarity annotation settings saved.")

    def _open_transfer_annotation_settings(self) -> None:
        from repeatviewer.backends.annotation_transfer import save_transfer_strategy
        from repeatviewer.backends.auto_annotation_settings import (
            save_annotation_run_options,
        )

        dialog = TransferAnnotationSettingsDialog(parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        threads, show, close = dialog.run_options()
        save_annotation_run_options(
            num_threads=threads,
            show_process_window=show,
            close_process_window_when_finished=close,
        )
        save_transfer_strategy(dialog.transfer_strategy())
        self._rail.refresh_transfer_run_tooltip()
        self._broadcast_shared_annotation_settings_changed()
        self.statusMessage.emit("Transfer annotation settings saved.")

    def _current_similarity_annotation_settings(self) -> SimilarityAnnotationSettings:
        return self._similarity_ann_settings

    def _reflow_similarity_database_chips(self) -> None:
        return

    def _make_similarity_db_check(self, entry: UserDatabase, index: int) -> QCheckBox:
        color = similarity_database_color(entry.name, index)
        check = SimilarityDatabaseCheck(
            entry.name, index, color, self._rail.similarity_db_host
        )
        readiness = database_search_readiness(entry)
        tooltip = f"{index + 1}. {entry.name}\n{entry.display_type_label}\n{entry.path}"
        if readiness is not None:
            tooltip += f"\n\nNot ready: {readiness}"
            check.setChecked(False)
        else:
            check.setChecked(True)
        check.setToolTip(tooltip)
        check.toggled.connect(self._on_similarity_db_check_toggled)
        return check

    def _refresh_similarity_database_rows(self) -> None:
        self._rail.clear_similarity_db_checks()
        self._similarity_db_rows = []
        databases = load_similarity_annotation_databases()
        if not databases:
            empty = QLabel("No similarity databases configured.", self._rail.similarity_db_host)
            empty.setStyleSheet("color: palette(mid);")
            empty.setWordWrap(False)
            self._rail.add_similarity_db_widget(empty)
            self._sync_similarity_annotation_all_check()
            return
        for index, entry in enumerate(databases):
            check = self._make_similarity_db_check(entry, index)
            self._rail.add_similarity_db_check(check)
            self._similarity_db_rows.append((entry, check))
        self._sync_similarity_annotation_all_check()

    def _selected_similarity_databases(self) -> list[SimilarityDatabaseRun]:
        selected: list[SimilarityDatabaseRun] = []
        for index, (entry, check) in enumerate(self._similarity_db_rows):
            if check.isChecked():
                selected.append(
                    SimilarityDatabaseRun(
                        database=entry,
                        color=similarity_database_color(entry.name, index),
                    )
                )
        return selected

    def _similarity_annotation_running(self) -> bool:
        return self._sim_ann_thread is not None and self._sim_ann_thread.isRunning()

    def _set_similarity_annotation_enabled(self, enabled: bool) -> None:
        self._similarity_ann_settings_btn.setEnabled(enabled)
        self._similarity_ann_all_check.setEnabled(enabled)
        for _entry, check in self._similarity_db_rows:
            check.setEnabled(enabled)
        self._similarity_ann_run_btn.setEnabled(True)
        self._set_similarity_run_button_running(not enabled)

    def _similarity_based_annotate(self) -> None:
        if self._similarity_annotation_running():
            return
        if self._map.data() is None:
            return
        databases = self._selected_similarity_databases()
        if not databases:
            QMessageBox.information(
                self,
                "Similarity annotation",
                "Select at least one database to search.",
            )
            return
        readiness_issues = []
        for entry in databases:
            readiness = database_search_readiness(entry.database)
            if readiness is not None:
                readiness_issues.append(f"{entry.database.name}: {readiness}")
        if readiness_issues:
            QMessageBox.warning(
                self,
                "Similarity annotation",
                "These databases are not ready to search:\n\n"
                + "\n\n".join(readiness_issues),
            )
            return
        if similarity_search_needs_orfs(tuple(databases)):
            map_orfs = [
                feat
                for feat in self._map.features()
                if str(feat.feature_type).upper() == "ORF"
                or str(feat.name or "").upper() == "ORF"
            ]
            if not map_orfs and not orfipy_available():
                QMessageBox.warning(
                    self,
                    "Similarity annotation",
                    "Protein similarity search needs ORF intervals on the sequence.\n\n"
                    "Run auto annotation with ORF checked first, or install orfipy:\n"
                    "  pip install orfipy",
                )
                return

        self._push_feature_undo_snapshot()
        self._sim_ann_cancel.reset()
        self._sim_ann_reveal_kinds = tuple(
            f"domain::{run.database.name}"
            for run in databases
            if str(run.database.name or "").strip()
        )
        self._set_similarity_annotation_enabled(False)
        labels = ", ".join(run.database.name for run in databases)
        self._log_annotation_job(
            "JOB START",
            "Similarity annotation",
            detail=f"databases=[{labels}]",
        )
        self._set_annotation_job_status(
            f"Similarity annotation (BLAST/HMM/CDD): {labels}…"
        )
        self.statusMessage.emit(f"Running similarity annotation: {labels}…")
        self._begin_window_job_progress("Similarity annotation", len(databases))

        worker = Worker(
            run_similarity_annotations,
            self._path,
            tuple(databases),
            settings=self._current_similarity_annotation_settings(),
            use_system_blast=load_use_system_blast_from_settings(),
            orf_settings=self._auto_ann_settings.orf,
            cancel=self._sim_ann_cancel,
            base_features=tuple(self._map.features()),
            progress=True,
        )
        from repeatviewer.widgets.process_monitor_dialog import (
            attach_annotation_worker_process_monitor,
        )

        cfg = self._similarity_ann_settings
        attach_annotation_worker_process_monitor(
            self,
            worker,
            title="Similarity annotation",
            enabled=bool(cfg.show_process_window),
            close_when_finished=bool(cfg.close_process_window_when_finished),
            tool_key="similarity_annotation",
        )
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.signals.progress.connect(self._on_similarity_annotation_progress)
        worker.signals.finished.connect(self._on_similarity_annotation_finished)
        worker.signals.error.connect(self._on_similarity_annotation_error)
        worker.signals.finished.connect(thread.quit)
        worker.signals.error.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._clear_similarity_annotation_thread)
        self._sim_ann_worker = worker
        self._sim_ann_thread = thread
        thread.start()

    @Slot()
    def _clear_similarity_annotation_thread(self) -> None:
        self._sim_ann_thread = None
        self._sim_ann_worker = None
        self._set_similarity_annotation_enabled(True)
        self._clear_annotation_job_status()
        self._end_window_job_progress()

    @Slot(object)
    def _reveal_annotation_kinds_after_run(self, keys: tuple[str, ...]) -> None:
        """Enable Displayed annotations for kinds just predicted."""
        if not keys:
            return
        ensure = getattr(self._rail, "ensure_annotation_kinds_visible", None)
        if callable(ensure):
            ensure(keys)
        self._map.set_visible_annotation_kinds(self._rail.visible_annotation_kinds())

    def _on_similarity_annotation_finished(self, result: object) -> None:
        if not isinstance(result, SimilarityAnnotationResult):
            self._log_annotation_job(
                "JOB ERROR",
                "Similarity annotation",
                detail="unexpected result",
                level=logging.ERROR,
            )
            self.statusMessage.emit("Similarity annotation finished with an unexpected result.")
            self._clear_annotation_job_status()
            return
        data = self._map.data()
        if data is None:
            self._clear_annotation_job_status()
            return
        self._map.set_features(list(result.features))
        self._reveal_annotation_kinds_after_run(self._sim_ann_reveal_kinds)
        self._clear_annotation_job_status()
        self._log_annotation_job(
            "JOB FINISH",
            "Similarity annotation",
            detail=result.message,
        )
        self.statusMessage.emit(f"{result.message} Click Save to write changes to disk.")

    @Slot(str)
    def _on_similarity_annotation_error(self, message: str) -> None:
        self._clear_annotation_job_status()
        if "stopped by user" in message.lower() or "job stopped" in message.lower():
            self._log_annotation_job(
                "JOB FINISH",
                "Similarity annotation",
                detail="stopped by user",
            )
            self.statusMessage.emit("Similarity annotation stopped.")
            return
        self._log_annotation_job(
            "JOB ERROR",
            "Similarity annotation",
            detail=message,
            level=logging.ERROR,
        )
        QMessageBox.warning(self, "Similarity annotation failed", message)
        self.statusMessage.emit(f"Similarity annotation failed: {message}")

    def _auto_annotate(self) -> None:
        if self._auto_annotate_running():
            return
        data = self._map.data()
        if data is None:
            return

        modes = self._selected_auto_annotation_modes()
        if not modes:
            QMessageBox.information(
                self,
                "Auto annotation",
                "Select at least one detector to run.",
            )
            return

        self._push_feature_undo_snapshot()
        self._auto_ann_cancel.reset()
        self._auto_ann_reveal_kinds = tuple(str(mode.value) for mode in modes)
        self._set_auto_annotation_enabled(False)
        labels = ", ".join(
            _AUTO_MODE_LABEL.get(mode, str(mode.value)) for mode in modes
        )
        self._log_annotation_job(
            "JOB START",
            "Auto annotation",
            detail=f"modes=[{labels}]",
        )
        self._set_annotation_job_status(f"Feature detection: {labels}…")
        self.statusMessage.emit(f"Running feature detection: {labels}…")
        self._begin_window_job_progress("Feature detection", len(modes))

        worker = Worker(
            run_auto_annotations,
            self._path,
            tuple(modes),
            cancel=self._auto_ann_cancel,
            settings=self._auto_ann_settings,
            base_features=tuple(self._map.features()),
            progress=True,
        )
        from repeatviewer.widgets.process_monitor_dialog import (
            attach_annotation_worker_process_monitor,
        )

        cfg = self._auto_ann_settings
        attach_annotation_worker_process_monitor(
            self,
            worker,
            title="Feature detection",
            enabled=bool(cfg.show_process_window),
            close_when_finished=bool(cfg.close_process_window_when_finished),
            tool_key="feature_detection",
        )
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.signals.progress.connect(self._on_auto_annotation_progress)
        worker.signals.finished.connect(self._on_auto_annotation_finished)
        worker.signals.error.connect(self._on_auto_annotation_error)
        worker.signals.finished.connect(thread.quit)
        worker.signals.error.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._clear_auto_annotation_thread)
        self._auto_ann_worker = worker
        self._auto_ann_thread = thread
        thread.start()

    @Slot()
    def _clear_auto_annotation_thread(self) -> None:
        self._auto_ann_thread = None
        self._auto_ann_worker = None
        self._set_auto_annotation_enabled(True)
        self._clear_annotation_job_status()
        self._end_window_job_progress()

    @Slot(object)
    def _on_auto_annotation_finished(self, result: object) -> None:
        reveal_kinds = tuple(self._auto_ann_reveal_kinds)
        if isinstance(result, MultiAutoAnnotationResult):
            data = self._map.data()
            if data is None:
                self._clear_annotation_job_status()
                return
            features = list(result.features)
            self._map.set_features(features)
            self._reveal_annotation_kinds_after_run(reveal_kinds)
            self._clear_annotation_job_status()
            self._log_annotation_job(
                "JOB FINISH",
                "Auto annotation",
                detail=result.message,
            )
            if features:
                self.statusMessage.emit(f"{result.message} Click Save to write changes to disk.")
            else:
                self.statusMessage.emit(result.message)
            return

        if not isinstance(result, AutoAnnotationResult):
            self._log_annotation_job(
                "JOB ERROR",
                "Auto annotation",
                detail="unexpected result",
                level=logging.ERROR,
            )
            self.statusMessage.emit("Auto annotation finished with an unexpected result.")
            self._clear_annotation_job_status()
            return

        data = self._map.data()
        if data is None:
            self._clear_annotation_job_status()
            return

        if not result.features:
            self._show_wide_information("Auto annotation", result.message)
            self._reveal_annotation_kinds_after_run(reveal_kinds)
            self._log_annotation_job(
                "JOB FINISH",
                "Auto annotation",
                detail=result.message,
            )
            self.statusMessage.emit(result.message)
            self._clear_annotation_job_status()
            return

        merged = merge_auto_annotation_features(data.features, result.features, result.mode)
        self._map.set_features(merged)
        self._reveal_annotation_kinds_after_run(reveal_kinds)
        self._clear_annotation_job_status()
        self._log_annotation_job(
            "JOB FINISH",
            "Auto annotation",
            detail=result.message,
        )
        self.statusMessage.emit(f"{result.message} Click Save to write changes to disk.")

    @Slot(str)
    def _on_auto_annotation_error(self, message: str) -> None:
        self._clear_annotation_job_status()
        if "stopped by user" in message.lower() or "job stopped" in message.lower():
            self._log_annotation_job(
                "JOB FINISH",
                "Auto annotation",
                detail="stopped by user",
            )
            self.statusMessage.emit("Auto annotation stopped.")
            return
        self._log_annotation_job(
            "JOB ERROR",
            "Auto annotation",
            detail=message,
            level=logging.ERROR,
        )
        QMessageBox.warning(self, "Auto annotation failed", message)
        self.statusMessage.emit(f"Auto annotation failed: {message}")

    def cancel_auto_annotation(self) -> bool:
        """Request cancellation of a running auto-annotation job."""
        if not self._auto_annotate_running():
            return False
        self._stop_auto_annotation()
        return True

    def cancel_similarity_annotation(self) -> bool:
        """Request cancellation of a running similarity-annotation job."""
        if not self._similarity_annotation_running():
            return False
        self._stop_similarity_annotation()
        return True

    def _export_gff(self) -> None:
        self._save_annotations(silent=True)

    def _open_exported_file(self, path: Path) -> None:
        """Open a successful export in Preview (or the system viewer as fallback)."""
        resolved = path.expanduser().resolve()
        if not resolved.is_file():
            return
        window = self.window()
        opener = getattr(window, "_auto_open_outputs", None)
        if callable(opener):
            try:
                opener([resolved])
                return
            except Exception:
                pass
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices

        QDesktopServices.openUrl(QUrl.fromLocalFile(str(resolved)))

    def _annotation_view_export_layout(
        self,
    ) -> tuple[QWidget | None, QWidget | None, QWidget | None, int, int, int, int]:
        """Header, canvas row, annot track, and tight export geometry.

        Uses annotation *content* height (not the viewport-stretched empty fill
        below feature rows) so PDF/PNG/SVG exports have no large blank footer.
        """
        widget = self._map
        header = getattr(widget, "_header", None)
        canvas_row = getattr(widget, "_canvas_row", None)
        annot_track = getattr(widget, "_annot_track", None)
        if header is None or canvas_row is None or annot_track is None:
            return None, None, None, 0, 0, 0, 0

        def _h(w: QWidget) -> int:
            h = int(w.height())
            if h <= 0:
                h = int(w.sizeHint().height())
            return max(1, h)

        def _w(w: QWidget) -> int:
            width = int(w.width())
            if width <= 0:
                width = int(w.sizeHint().width())
            return max(1, width)

        header_h = _h(header)
        canvas_h = _h(canvas_row)
        annot_h = max(1, int(annot_track.content_height()))
        width = max(_w(header), _w(canvas_row), _w(annot_track))
        height = header_h + canvas_h + annot_h
        return header, canvas_row, annot_track, width, height, header_h, canvas_h

    def _export_annotation_view_to_file(
        self,
        out_path: Path,
        *,
        dpi: int = 400,
    ) -> bool:
        """Export the current SequenceMapWidget view to disk. Return True on success."""
        from repeatviewer.views.export_annotation_view_dialog import (
            vector_export,
            write_view_vector_pdf,
        )

        widget = self._map
        (
            header,
            canvas_row,
            annot_track,
            width,
            height,
            header_h,
            canvas_h,
        ) = self._annotation_view_export_layout()
        if header is None or width <= 0 or height <= 0:
            width = max(1, int(widget.rect().width()))
            height = max(1, int(widget.rect().height()))
            header = canvas_row = annot_track = None
            header_h = canvas_h = 0
        if width <= 0 or height <= 0:
            QMessageBox.information(self, "Export", "Nothing to export (empty view).")
            return False
        size = QSize(width, height)
        dpi = max(50, min(600, int(dpi)))
        scale = max(1.0, float(dpi) / 96.0)

        def _paint_view(painter: QPainter) -> None:
            if header is not None and canvas_row is not None and annot_track is not None:
                y = 0
                header.render(painter, QPoint(0, y))
                y += int(header_h)
                canvas_row.render(painter, QPoint(0, y))
                y += int(canvas_h)
                # Track widget may be taller than content (fills the Preview pane).
                source = QRect(0, 0, int(width), int(annot_track.content_height()))
                annot_track.render(painter, QPoint(0, y), QRegion(source))
            else:
                widget.render(painter, QPoint(0, 0))

        suffix = out_path.suffix.lower()
        if suffix == ".svg":
            try:
                from PySide6.QtSvg import QSvgGenerator  # type: ignore
            except Exception:  # noqa: BLE001
                QMessageBox.warning(
                    self,
                    "Export",
                    "SVG export needs PySide6 QtSvg module. Install it and retry.",
                )
                return False

            gen = QSvgGenerator()
            gen.setFileName(str(out_path))
            gen.setSize(size)
            gen.setViewBox(QRect(0, 0, size.width(), size.height()))
            gen.setResolution(dpi)
            painter = QPainter(gen)
            try:
                with vector_export(self):
                    _paint_view(painter)
            finally:
                painter.end()
            return True

        if suffix == ".pdf":
            # Vector PDF: paint the widgets straight onto the page (no bitmap).
            with vector_export(self):
                write_view_vector_pdf(
                    out_path,
                    _paint_view,
                    width=width,
                    height=height,
                    logical_dpi=float(widget.logicalDpiY()),
                )
            return True

        # PNG: render at DPI scale.
        pixel = QSize(
            max(1, int(round(size.width() * scale))),
            max(1, int(round(size.height() * scale))),
        )
        img = QImage(pixel, QImage.Format.Format_ARGB32)
        img.fill(Qt.GlobalColor.transparent)
        painter = QPainter(img)
        try:
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            painter.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
            painter.scale(scale, scale)
            _paint_view(painter)
        finally:
            painter.end()

        if suffix == ".png":
            if not img.save(str(out_path), "PNG"):
                QMessageBox.warning(self, "Export", "Could not save PNG.")
                return False
            return True

        QMessageBox.warning(
            self,
            "Export",
            f"Unsupported export format: {suffix!r}. Use .svg, .png, or .pdf.",
        )
        return False

    def _annotations_for_export(self) -> tuple[list[str], dict, dict]:
        """Features as shown now (incl. unsaved edits) for GFF3/GTF/BED export."""
        data = self._loaded_data if self._loaded_data is not None else self._map.data()
        if data is None:
            return [], {}, {}
        data.sync_active_features(self._map.features())
        lengths = dict(data.lengths_by_id)
        if data.sequence_id:
            lengths.setdefault(data.sequence_id, data.length)
        features = {
            rid: list(feats) for rid, feats in data.features_by_id.items()
        }
        if data.sequence_id and data.sequence_id not in features:
            features[data.sequence_id] = list(self._map.features())
        record_ids = list(data.record_ids) or list(features.keys())
        return record_ids, features, lengths

    def _export_annotation_view(self) -> None:
        from repeatviewer.views.export_annotation_view_dialog import (
            ExportAnnotationViewDialog,
        )

        dialog = ExportAnnotationViewDialog(
            source_path=self._path,
            title="Export",
            name_tag="annotation",
            default_format="pdf",
            annotation_provider=self._annotations_for_export,
            parent=self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        out = dialog.chosen_path()
        if out is None:
            return
        if dialog.wrote_annotation_file():
            self.statusMessage.emit(f"Exported annotations → {out.name}")
            # Shows the new GFF3/GTF/BED in Files without leaving this preview.
            self.annotationsSaved.emit(out)
            QTimer.singleShot(0, self._map.setFocus)
            return
        if self._export_annotation_view_to_file(out, dpi=dialog.dpi()):
            self.statusMessage.emit(f"Exported annotation view → {out.name}")
            self._open_exported_file(out)
        # Avoid letting focus fall back to the left Files panel after export.
        QTimer.singleShot(0, self._map.setFocus)
