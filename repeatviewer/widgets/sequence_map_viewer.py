"""Interactive single-sequence map: scroll, zoom to nucleotide, feature tracks."""

from __future__ import annotations

import logging
import math
import platform
from pathlib import Path

from PySide6.QtCore import (
    QEvent,
    QPoint,
    QPointF,
    QRect,
    QRectF,
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
    QFont,
    QFontMetrics,
    QKeyEvent,
    QKeySequence,
    QMouseEvent,
    QNativeGestureEvent,
    QPainter,
    QImage,
    QPalette,
    QPen,
    QPixmap,
    QPolygonF,
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
    QStyle,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.spin_boxes import (
    NoWheelDoubleSpinBox as QDoubleSpinBox,
    NoWheelSpinBox as QSpinBox,
)

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
    DEFAULT_TIR_SIMILARITY_THRESHOLD,
    OrfSettings,
    TerminalRepeatSettings,
    TandemRepeatSettings,
    load_auto_annotation_settings,
    save_auto_annotation_settings,
)
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
    features_gff_path,
    format_single_sequence_summary_parts,
    format_single_sequence_summary_text,
    load_sequence_map_data,
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
from repeatviewer.job_options import configure_report_evalue_spinbox
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
from repeatviewer.widgets.preview_controls import (
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
)
_AUTO_MODE_LABEL = {
    AutoAnnotationMode.LTR: "LTR",
    AutoAnnotationMode.TIR: "TIR",
    AutoAnnotationMode.TANDEM: "Tandem",
    AutoAnnotationMode.ORF: "ORF",
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
_ANNOT_TRACK_VPAD = 8
_ANNOT_FEATURE_HEIGHT = 18
logger = logging.getLogger(__name__)
_IMAGE_MARKER_SIZE = 12  # cytogenetic-image dot on annotation features
_ANNOT_MIN_VIEWPORT_HEIGHT = 56
_OVERVIEW_SYNC_MS = 16
_RULER_HEIGHT = 20
_VIEWPORT_RULER_HEIGHT = 14
_OVERVIEW_HEIGHT = 16
_SEQUENCE_ROW = 14
_TRACK_HEIGHT = 12
_TRACK_GAP = 0
_SELECTION_HIGHLIGHT_MARGIN = 2
_FEATURE_UNDO_LIMIT = 32
_LEFT_MARGIN = 8
_RIGHT_MARGIN = 8

operation_logger = logging.getLogger("repeatviewer.operations")
_MIN_PPB = 0.05
_FIT_MIN_PPB = 1e-6
_MAX_PPB = 48.0
_INTERACTIVE_SEQUENCE_LEN = 35_000
_LARGE_SEQUENCE_LEN = 150_000
_HUGE_SEQUENCE_LEN = 2_000_000
_LAYOUT_SYNC_MS = 40
_BUCKET_DRAW_PPB = 1.0
_WRAP_MIN_PPB = 4.0
_WRAP_LINE_RULER = 12
_WRAP_LINE_GAP = 6
_NO_COLOR_FG = QColor("#222222")
_NO_COLOR_BG = QColor("#fafafa")
_SELECTION_FILL = QColor(10, 40, 120, 100)
_SELECTION_BORDER = QColor(5, 20, 70, 255)
_FEATURE_SELECTION_BORDER = QColor("#0b5cab")
_FEATURE_SELECTION_OUTLINE = QColor("#ffffff")
_FEATURE_SELECTION_OVERLAY = QColor(11, 92, 171, 120)
_MIN_NUC_FONT_PT = 6
_MAX_NUC_FONT_PT = 36
_DEFAULT_NUC_FONT_PT = 13
_ZOOM_BTN_SIZE = 28
_ZOOM_BTN_HEIGHT = 24
_ANNOTATION_SIDE_PANEL_WIDTH = 108


def _zoom_ppb_max() -> float:
    """Maximum zoom (pixels per base) for the annotation slider."""
    return max(float(_WRAP_MIN_PPB), float(_MAX_PPB))


def _zoom_ppb_max_for(seq_len: int) -> float:
    """Cap zoom-in for very long sequences to keep scrolling responsive."""
    length = max(0, int(seq_len))
    if length <= _INTERACTIVE_SEQUENCE_LEN:
        return _MAX_PPB
    if length <= _LARGE_SEQUENCE_LEN:
        return 16.0
    if length <= _HUGE_SEQUENCE_LEN:
        return 8.0
    return 4.0


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


def _base_colors(base: str, *, is_protein: bool = False) -> tuple[QColor, QColor]:
    if is_protein:
        from repeatviewer.widgets.preview_zoom import amino_acid_color_map

        palette = amino_acid_color_map()
        key = base.upper()
        bg, fg = palette.get(key, palette.get("X", ("#C8C8C8", "#000000")))
        return QColor(bg), QColor(fg)
    palette = nucleotide_base_color_map()
    bg, fg = palette.get(base.upper(), palette["N"])
    return QColor(bg), QColor(fg)


def _draw_directed_feature(
    painter: QPainter,
    rect: QRectF,
    color: QColor,
    *,
    forward: bool,
    label: str = "",
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
    if label and rect.width() > 28:
        painter.setPen(QColor("#ffffff"))
        label_font = QFont(painter.font())
        label_font.setPointSize(7)
        painter.setFont(label_font)
        text_rect = rect.adjusted(3, 0, -3, 0)
        if not forward:
            text_rect = text_rect.adjusted(arrow_w * 0.4, 0, 0, 0)
        else:
            text_rect = text_rect.adjusted(0, 0, -arrow_w * 0.4, 0)
        painter.drawText(text_rect, Qt.AlignmentFlag.AlignVCenter, label)


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


def _pack_feature_rows(features: list[SequenceFeature]) -> tuple[list[int], int]:
    """Assign rows per track, then stack tracks vertically."""
    if not features:
        return [], 0
    tracks = sorted({feat.track for feat in features})
    row_for_index = [0] * len(features)
    row_offset = 0
    for track in tracks:
        indices = [idx for idx, feat in enumerate(features) if feat.track == track]
        group = [features[idx] for idx in indices]
        group_rows, group_count = _pack_feature_row_group(group)
        for local_idx, global_idx in enumerate(indices):
            row_for_index[global_idx] = group_rows[local_idx] + row_offset
        row_offset += group_count
    return row_for_index, row_offset


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
    """Scroll area that pans horizontally from trackpad / Shift+wheel input."""

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802
        hbar = self.horizontalScrollBar()
        delta_x = _horizontal_wheel_delta(event, scroll_step=hbar.singleStep() * 3)
        if delta_x != 0 and hbar.maximum() > 0:
            hbar.setValue(hbar.value() - delta_x)
            event.accept()
            return
        super().wheelEvent(event)


class ViewportRulerWidget(QWidget):
    """Position scale aligned with the visible sequence (scroll + zoom synced)."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedHeight(_VIEWPORT_RULER_HEIGHT)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._sequence = ""
        self._scroll_px = 0.0
        self._ppb = 1.0
        self._viewport_w = 200

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
        self._ppb = max(0.01, float(ppb))
        self._viewport_w = max(1, int(viewport_w))
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor("#fafafa"))
        seq_len = len(self._sequence)
        if seq_len <= 0:
            painter.end()
            return
        ppb = self._ppb
        scroll = self._scroll_px
        vw = self._viewport_w
        ruler_font = _monospace_font(point_size=7)
        painter.setFont(ruler_font)
        fm = QFontMetrics(ruler_font)
        painter.setPen(QColor("#666666"))

        start_bp = max(1, int((scroll - _LEFT_MARGIN) / ppb) + 1)
        end_bp = min(seq_len, int((scroll + vw - _LEFT_MARGIN) / ppb) + 1)
        major_step = SequenceMapCanvas._nice_tick_step(seq_len, ppb, fm)
        tick_step = major_step if major_step >= 5 else 5
        min_label_w = fm.horizontalAdvance("00000") + 4
        last_label_x = -1e9

        for bp in range(start_bp, end_bp + 1):
            if bp != 1 and bp != seq_len and bp % tick_step != 0:
                continue
            x = _bp_view_x(bp, scroll, ppb)
            if x < _LEFT_MARGIN or x > vw - 4:
                continue
            if bp != seq_len:
                painter.drawLine(int(x), 1, int(x), _VIEWPORT_RULER_HEIGHT - 3)
            text = str(bp)
            if x - last_label_x < min_label_w and bp != start_bp:
                continue
            painter.drawText(int(x + 2), _VIEWPORT_RULER_HEIGHT - 2, text)
            last_label_x = x
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
    deleteRequested = Signal()
    copySequenceRequested = Signal()
    extractSequenceRequested = Signal()
    modifyRequested = Signal()
    changeColorRequested = Signal()
    tracksMoved = Signal(object)  # list[SequenceFeature]
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
        self._selected_indices: set[int] = set()
        self._row_count = 0
        self._marquee_active = False
        self._marquee_origin: QPointF | None = None
        self._marquee_current: QPointF | None = None
        self._press_pos: QPointF | None = None
        self._press_hit: int | None = None
        self._marquee_additive_base: set[int] = set()
        self._mouse_grabbed = False
        self._show_similarity_database_name = True
        self._move_active = False
        self._move_indices: set[int] = set()
        self._move_preview_track: int | None = None
        self._move_ghost_y: float | None = None
        self._hover_idx: int | None = None
        self.setFixedHeight(0)

    def set_show_similarity_database_name(self, show: bool) -> None:
        show = bool(show)
        if show == self._show_similarity_database_name:
            return
        self._show_similarity_database_name = show
        self.update()

    def _feature_track_label(self, feat: SequenceFeature, rect_width: float, min_width: float) -> str:
        if rect_width <= min_width:
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
        self._marquee_active = False
        self._move_active = False
        self._move_indices = set()
        self._move_preview_track = None
        self._move_ghost_y = None
        self._hover_idx = None
        self._release_drag_mouse()
        super().hideEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802
        if self._hover_idx is not None and not self._move_active and not self._marquee_active:
            self._hover_idx = None
            self.update()
        super().leaveEvent(event)

    def _row_height(self) -> float:
        return float(_ANNOT_FEATURE_HEIGHT + _TRACK_GAP)

    def _feature_full_label(self, feat: SequenceFeature) -> str:
        base = feat.name or feat.feature_type
        return similarity_feature_track_label(
            base,
            feat.feature_type,
            show_database_name=self._show_similarity_database_name,
        )

    def _label_overflows(self, feat: SequenceFeature, rect: QRectF) -> bool:
        label = self._feature_full_label(feat)
        if not label or rect.width() <= 28:
            return bool(label) and rect.width() <= 28
        fm = QFontMetrics(self.font())
        # Match the 7pt label font used while painting features.
        label_font = QFont(self.font())
        label_font.setPointSize(7)
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
        return QRectF(x1, y, max(2.0, x2 - x1), max(2.0, float(_ANNOT_FEATURE_HEIGHT) - 1.0))

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
    ) -> None:
        self._features = list(features)
        self._feature_rows = list(feature_rows)
        self._scroll_px = max(0.0, float(scroll_px))
        self._ppb = max(0.01, float(ppb))
        self._viewport_w = max(1, int(viewport_w))
        if (
            not self._marquee_active
            and not self._move_active
            and self._press_pos is None
        ):
            self._selected_indices = set(selected_indices or [])
        valid = set(range(len(self._features)))
        self._selected_indices &= valid
        if self._hover_idx is not None and self._hover_idx not in valid:
            self._hover_idx = None
        self._row_count = max((r for r in feature_rows), default=-1) + 1 if features else 0
        track_w = max(1, int(viewport_w))
        if self._row_count > 0:
            h = int(self._row_count * self._row_height()) + 2 * _ANNOT_TRACK_VPAD
            h = max(h, _ANNOT_FEATURE_HEIGHT + 2 * _ANNOT_TRACK_VPAD)
            self.setFixedSize(track_w, h)
        else:
            # Keep a visible empty track so the annotation panel is not blank.
            self.setFixedSize(track_w, max(_ANNOT_MIN_VIEWPORT_HEIGHT, _ANNOT_FEATURE_HEIGHT + 2 * _ANNOT_TRACK_VPAD))
        self.show()
        self.update()

    def has_content(self) -> bool:
        return self._row_count > 0

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
        painter.fillRect(self.rect(), QColor("#f5f5f5"))
        if not self._features:
            self._marked_indices = set()
            painter.end()
            return
        self._refresh_marked_indices()
        ppb = self._ppb
        scroll = self._scroll_px
        vw = self._viewport_w
        row_h = self._row_height()
        painter.setPen(Qt.PenStyle.NoPen)

        for idx, feat in enumerate(self._features):
            x1 = _bp_view_x(feat.start, scroll, ppb)
            x2 = _bp_view_x(feat.end + 1, scroll, ppb)
            if x2 < _LEFT_MARGIN or x1 > vw:
                continue
            row = self._feature_rows[idx] if idx < len(self._feature_rows) else 0
            y = _ANNOT_TRACK_VPAD + row * row_h + 0.5
            color = QColor(feat.color or feature_type_color(feat.feature_type))
            rect = QRectF(x1, y, max(2.0, x2 - x1), max(2.0, float(_ANNOT_FEATURE_HEIGHT) - 1.0))
            label = self._feature_track_label(feat, rect.width(), 28.0)
            selected = idx in self._selected_indices
            _draw_directed_feature(
                painter,
                rect,
                color,
                forward=feat.is_forward(),
                label=label,
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
                    max(2.0, float(_ANNOT_FEATURE_HEIGHT) - 1.0),
                )
                color = QColor(feat.color or feature_type_color(feat.feature_type))
                color.setAlpha(110)
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
        ):
            feat = self._features[hover_idx]
            rect = self._feature_rect(hover_idx)
            full_label = self._feature_full_label(feat)
            if full_label and self._label_overflows(feat, rect):
                tip_font = QFont(painter.font())
                tip_font.setPointSize(8)
                painter.setFont(tip_font)
                fm = QFontMetrics(tip_font)
                tip_w = min(
                    max(fm.horizontalAdvance(full_label) + 10, 40),
                    max(80, self.width() - 8),
                )
                tip_h = fm.height() + 6
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
                painter.drawText(
                    tip_rect.adjusted(5, 0, -5, 0),
                    Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                    full_label,
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
        self._press_hit = self._hit_test_feature(pos)
        self._marquee_active = False
        self._marquee_origin = None
        self._marquee_current = None
        self._move_active = False
        self._move_indices = set()
        self._move_preview_track = None
        self._move_ghost_y = None
        self.setFocus(Qt.FocusReason.MouseFocusReason)
        event.accept()

    def mouseMoveEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        pos = event.position()
        if self._press_pos is None:
            hit = self._hit_test_feature(pos)
            if hit != self._hover_idx:
                self._hover_idx = hit
                self.update()
            super().mouseMoveEvent(event)
            return

        if self._move_active:
            self._move_ghost_y = pos.y() - float(_ANNOT_FEATURE_HEIGHT) * 0.5
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

        # Drag from an already-selected feature (mostly vertical) → move track.
        # Otherwise → marquee multi-select.
        if (
            self._can_start_move(self._press_hit)
            and not self._additive_modifier(event.modifiers())
            and abs(delta.y()) >= abs(delta.x())
        ):
            self._move_active = True
            self._move_indices = set(self._selected_indices)
            self._move_ghost_y = pos.y() - float(_ANNOT_FEATURE_HEIGHT) * 0.5
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
        if self._press_pos is None and not self._marquee_active and not self._move_active:
            self._release_drag_mouse()
            super().mouseReleaseEvent(event)
            return
        pos = event.position()
        additive = self._additive_modifier(event.modifiers())
        try:
            if self._move_active and self._move_indices:
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
                self._emit_selection()
        finally:
            self._press_pos = None
            self._press_hit = None
            self._marquee_active = False
            self._marquee_origin = None
            self._marquee_current = None
            self._marquee_additive_base = set()
            self._move_active = False
            self._move_indices = set()
            self._move_preview_track = None
            self._move_ghost_y = None
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
    featureClicked = Signal(int)  # index in features list
    addFeatureRequested = Signal(str, str)  # feature_type, label
    extractSequenceRequested = Signal(int, int)  # 1-based inclusive range
    featureActionRequested = Signal(str)  # delete | copy | extract | modify

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
        self._drag_select = False
        self._mouse_grabbed = False
        self._pending_selection_emit: tuple[int, int] | None = None
        self._selection_emit_timer = QTimer(self)
        self._selection_emit_timer.setSingleShot(True)
        self._selection_emit_timer.setInterval(33)  # ~30 Hz spin/status updates
        self._selection_emit_timer.timeout.connect(self._flush_pending_selection_emit)
        self._seq_layer_cache: QPixmap | None = None
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._show_similarity_database_name = True

    def _invalidate_seq_layer_cache(self) -> None:
        self._seq_layer_cache = None

    def _ensure_seq_layer_cache(self) -> None:
        """Lazy sequence (+ wrap features) snapshot without selection chrome.

        Rebuilt on zoom/scroll/sequence change — not on every mouse press — so
        selection drag only blits the cache and redraws the highlight.
        """
        if (
            self._seq_layer_cache is not None
            and self._seq_layer_cache.size() == self.size()
            and not self._seq_layer_cache.isNull()
        ):
            return
        if self.width() < 2 or self.height() < 2 or not self._sequence:
            self._seq_layer_cache = None
            return
        cache = QPixmap(self.size())
        cache.fill(QColor("#fafafa"))
        painter = QPainter(cache)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, False)
        self._draw_sequence(painter, QRectF(cache.rect()))
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
        if rect_width <= min_width:
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
        # In the single-line view, leave room for the range-selection border
        # above the bases so the DNA row can sit in the highlight's middle.
        return 0 if self.uses_wrap() else _SELECTION_HIGHLIGHT_MARGIN

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
        self._feature_rows, self._feature_row_count = _pack_feature_rows(self._features)
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
        self._select_start_bp = None
        self._select_end_bp = None
        self._selected_indices.clear()
        self.update()

    def selection_range(self) -> tuple[int, int] | None:
        if self._select_start_bp is None or self._select_end_bp is None:
            return None
        a, b = sorted((self._select_start_bp, self._select_end_bp))
        return a, b

    def _sequence_font(self) -> QFont:
        font = _monospace_font(point_size=self._font_point_size)
        return font

    def _font_letter_cell_size(self) -> tuple[int, int]:
        """Minimum monospace cell size for legible nucleotide / amino-acid letters."""
        fm = QFontMetrics(self._sequence_font())
        sample = "ACDEFGHIKLMNPQRSTVWYBXZ*" if self._is_protein else "ACGTUN"
        cell_w = max(1, max(fm.horizontalAdvance(base) for base in sample))
        cell_h = max(_SEQUENCE_ROW, fm.height())
        return cell_w, max(1, cell_h)

    def shows_letters(self) -> bool:
        if len(self._sequence) > _INTERACTIVE_SEQUENCE_LEN:
            return False
        font_w, _ = self._font_letter_cell_size()
        return self._pixels_per_bp >= font_w

    def _letter_cell_size(self) -> tuple[int, int]:
        """Current nucleotide cell size when drawing letters."""
        cell_w = max(1, int(math.ceil(self._pixels_per_bp)))
        _, cell_h = self._font_letter_cell_size()
        return cell_w, max(_SEQUENCE_ROW, cell_h)

    def _nucleotide_cell_geometry(self) -> tuple[int, int, float, int]:
        """Shared top-left origin and cell size for every nucleotide background."""
        origin_y = int(self._sequence_top())
        ppb = max(_FIT_MIN_PPB, float(self._pixels_per_bp))
        font_w, font_h = self._font_letter_cell_size()
        if self.shows_letters():
            cell_w = max(float(font_w), ppb)
            cell_h = max(_SEQUENCE_ROW, font_h)
        else:
            cell_w = ppb
            _, cell_h, _ = self._cell_size()
            cell_h = max(_SEQUENCE_ROW, int(math.ceil(cell_h)))
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
        if self.shows_letters():
            _, cell_h = self._letter_cell_size()
            return float(cell_h)
        _, cell_h, _ = self._cell_size()
        return max(float(_SEQUENCE_ROW), cell_h)

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
        return int(self._seq_row_height())

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
            int(self._seq_row_height()) + 2 * _SELECTION_HIGHLIGHT_MARGIN + 2,
        )

    def _update_geometry(self) -> None:
        height = self._content_height()
        if self.height() != height:
            self.setFixedHeight(height)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self.sync_viewport_from_widget()

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
                painter.setPen(QPen(QColor(0, 0, 0, 60), 1))
                painter.drawLine(cell.topRight(), cell.bottomRight())
        if show_letters:
            painter.setPen(fg)
            painter.drawText(cell, Qt.AlignmentFlag.AlignCenter, base)

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
            return

        if not show_letters:
            self._draw_sequence_bar_viewport(painter, clip)
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
        if show_letters:
            for i in range(visible_start, visible_end):
                cell = self._base_rect_at(i)
                self._draw_base_cell(
                    painter,
                    self._sequence[i],
                    cell,
                    show_letters=True,
                )

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
            painter.fillRect(event.rect(), QColor("#fafafa"))
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No sequence")
            painter.end()
            return

        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, False)

        # Blit cached sequence layer; rebuild only after zoom/scroll/etc.
        self._ensure_seq_layer_cache()
        if self._seq_layer_cache is not None:
            painter.drawPixmap(event.rect(), self._seq_layer_cache, event.rect())
        else:
            painter.fillRect(event.rect(), QColor("#fafafa"))
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
        sel = self.selection_range()
        if sel:
            a, b = sel
            origin_x, origin_y, cell_w, cell_h = self._nucleotide_cell_geometry()
            if self.uses_wrap():
                cell_w_f, cell_h_f, _ = self._cell_size()
                for bp in range(a, b + 1):
                    x, y = self._index_to_xy(bp - 1)
                    rect = QRectF(x, y, cell_w_f, cell_h_f)
                    if not clip.intersects(rect):
                        continue
                    painter.fillRect(rect, _SELECTION_FILL)
                    painter.setPen(QPen(_SELECTION_BORDER, 1))
                    painter.drawRect(rect)
            else:
                x1 = origin_x + (a - 1) * cell_w - self._scroll_px
                x2 = origin_x + b * cell_w - self._scroll_px
                rect = QRectF(
                    x1,
                    float(origin_y - _SELECTION_HIGHLIGHT_MARGIN),
                    x2 - x1,
                    float(cell_h + 2 * _SELECTION_HIGHLIGHT_MARGIN),
                )
                if clip.intersects(rect):
                    painter.fillRect(rect, _SELECTION_FILL)
                    painter.setPen(QPen(_SELECTION_BORDER, 1))
                    painter.drawRect(rect)
            painter.setPen(Qt.PenStyle.NoPen)
        painter.end()

    def _selection_paint_rect(self, sel: tuple[int, int]) -> QRect:
        """Pixel rect for a selection range (including highlight margin)."""
        a, b = sel
        origin_x, origin_y, cell_w, cell_h = self._nucleotide_cell_geometry()
        pad = 3
        if self.uses_wrap():
            cell_w_f, cell_h_f, _ = self._cell_size()
            x0, y0 = self._index_to_xy(a - 1)
            x1, y1 = self._index_to_xy(b - 1)
            left = min(x0, x1) - pad
            top = min(y0, y1) - pad
            right = max(x0 + cell_w_f, x1 + cell_w_f) + pad
            bottom = max(y0 + cell_h_f, y1 + cell_h_f) + pad
            return QRect(int(left), int(top), int(right - left), int(bottom - top)).intersected(
                self.rect()
            )
        x1 = origin_x + (a - 1) * cell_w - self._scroll_px
        x2 = origin_x + b * cell_w - self._scroll_px
        top = int(origin_y - _SELECTION_HIGHLIGHT_MARGIN) - pad
        height = int(cell_h + 2 * _SELECTION_HIGHLIGHT_MARGIN) + 2 * pad
        return QRect(
            int(min(x1, x2)) - pad,
            top,
            int(abs(x2 - x1)) + 2 * pad,
            height,
        ).intersected(self.rect())

    def _selection_change_damage(
        self, old: tuple[int, int] | None, new: tuple[int, int] | None
    ) -> QRect:
        """Repaint only edges that changed — not the full old∪new union.

        Combined with the persistent sequence-layer cache, drag paints avoid
        re-walking every base on each mouse move.
        """
        if old is None and new is None:
            return QRect()
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
        # Identical range (e.g. press): refresh the thin caret/selection chrome.
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
        super().mouseMoveEvent(event)

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
        raw = max(1, int(target_px / max(ppb, 0.01)))
        for step in (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 5000, 10000, 50000):
            if step >= raw:
                return min(step, seq_len)
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

        add_menu = menu.addMenu("Add annotation")
        for label, feature_type in load_add_feature_buttons():
            action = add_menu.addAction(label)
            action.triggered.connect(
                lambda _checked=False, ft=feature_type, bl=label: self.addFeatureRequested.emit(
                    ft, bl
                )
            )

        menu.addSeparator()
        selection = self.selection_range()
        extract = menu.addAction("Extract Sequence")
        if selection is None:
            extract.setEnabled(False)
        else:
            start, end = selection
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
        self._start.setValue(max(1, min(int(start), self._sequence_length)))
        self._start.setGroupSeparatorShown(True)

        self._end = QSpinBox(self)
        self._end.setRange(1, self._sequence_length)
        self._end.setValue(max(1, min(int(end), self._sequence_length)))
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

        self._start.valueChanged.connect(self._sync)
        self._end.valueChanged.connect(self._sync)
        self._strand.currentIndexChanged.connect(self._sync)
        self._name.textChanged.connect(self._sync)
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
        if embedded:
            self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        else:
            self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        self._zoom_slider: QSlider | QSpinBox | None = None
        self._annotation_header_zoom_bar = None
        self._annotation_side_panel: QWidget | None = None
        self._toolbar: QWidget | None = None
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

        self._pos_start = QSpinBox(controls_parent)
        self._pos_start.setMinimum(1)
        self._pos_start.setMaximum(1)
        self._pos_start.setFixedWidth(72)
        self._pos_start.setKeyboardTracking(False)
        self._pos_end = QSpinBox(controls_parent)
        self._pos_end.setMinimum(1)
        self._pos_end.setMaximum(1)
        self._pos_end.setFixedWidth(72)
        self._pos_end.setKeyboardTracking(False)
        self._pos_start.valueChanged.connect(self._on_position_spin_changed)
        self._pos_end.valueChanged.connect(self._on_position_spin_changed)
        self._pos_start.editingFinished.connect(self._on_position_editing_finished)
        self._pos_end.editingFinished.connect(self._on_position_editing_finished)
        select_btn = QPushButton("Select", controls_parent)
        select_btn.setToolTip(
            "Zoom to the selected range and center its midpoint in the window"
        )
        select_btn.clicked.connect(self._select_by_position)

        if embedded:
            # Horizontal toolbar matching the standalone annotation window layout:
            # Zoom | Position | Select | Font | … | Length / GC
            toolbar = QWidget(self)
            toolbar.setObjectName("AnnotationMapToolbar")
            toolbar.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            toolbar.setFixedHeight(28)
            tb = QHBoxLayout(toolbar)
            tb.setContentsMargins(6, 2, 8, 2)
            tb.setSpacing(6)
            tb.setAlignment(Qt.AlignmentFlag.AlignVCenter)

            zoom_out_btn = QPushButton("−", toolbar)
            zoom_out_btn.setFixedSize(_ZOOM_BTN_SIZE, _ZOOM_BTN_HEIGHT)
            zoom_out_btn.setToolTip("Zoom out")
            zoom_out_btn.clicked.connect(self.zoom_out)
            zoom_in_btn = QPushButton("+", toolbar)
            zoom_in_btn.setFixedSize(_ZOOM_BTN_SIZE, _ZOOM_BTN_HEIGHT)
            zoom_in_btn.setToolTip("Zoom in")
            zoom_in_btn.clicked.connect(self.zoom_in)
            self._zoom_slider = QSlider(Qt.Orientation.Horizontal, toolbar)
            self._zoom_slider.setObjectName("annotationZoomSlider")
            self._zoom_slider.setRange(0, 100)
            self._zoom_slider.setValue(0)
            self._zoom_slider.setFixedWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
            self._zoom_slider.setMinimumWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
            self._zoom_slider.setFixedHeight(ZOOM_SLIDER_HEIGHT)
            self._zoom_slider.setToolTip(
                "Left: fit full sequence in window. Right: maximum zoom-in."
            )
            apply_annotation_zoom_slider_style(self._zoom_slider)
            self._zoom_slider.valueChanged.connect(self._on_zoom_changed)

            for w in (self._pos_start, self._pos_end, self._font_size, select_btn):
                w.setParent(toolbar)
            self._sequence_stats_label.setParent(toolbar)
            self._pos_start.setFixedWidth(64)
            self._pos_end.setFixedWidth(64)
            self._font_size.setFixedWidth(70)
            for spin in (self._pos_start, self._pos_end, self._font_size):
                spin.setFixedHeight(18)
            select_btn.setFixedHeight(20)
            select_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

            tb.addWidget(zoom_out_btn, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(zoom_in_btn, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(self._zoom_slider, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(QLabel("Position:", toolbar), 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(self._pos_start, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(QLabel("to", toolbar), 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(self._pos_end, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(select_btn, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(QLabel("Font:", toolbar), 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addWidget(self._font_size, 0, Qt.AlignmentFlag.AlignVCenter)
            tb.addStretch(1)
            self._sequence_stats_label.setAlignment(
                Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
            )
            self._sequence_stats_label.setWordWrap(False)
            tb.addWidget(self._sequence_stats_label, 0, Qt.AlignmentFlag.AlignVCenter)
            self._toolbar = toolbar
        else:
            zoom_row = QHBoxLayout()
            zoom_row.setContentsMargins(6, 4, 6, 4)
            zoom_row.setSpacing(6)
            zoom_row.setAlignment(Qt.AlignmentFlag.AlignVCenter)

            zoom_strip = QWidget(self)
            zoom_strip.setFixedHeight(_ZOOM_BTN_HEIGHT)
            zoom_strip.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            zoom_strip_layout = QHBoxLayout(zoom_strip)
            zoom_strip_layout.setContentsMargins(0, 0, 0, 0)
            zoom_strip_layout.setSpacing(4)
            zoom_strip_layout.setAlignment(Qt.AlignmentFlag.AlignVCenter)

            zoom_out_btn = QPushButton("−", zoom_strip)
            zoom_out_btn.setFixedSize(_ZOOM_BTN_SIZE, _ZOOM_BTN_HEIGHT)
            zoom_out_btn.setToolTip("Zoom out")
            zoom_out_btn.clicked.connect(self.zoom_out)
            zoom_in_btn = QPushButton("+", zoom_strip)
            zoom_in_btn.setFixedSize(_ZOOM_BTN_SIZE, _ZOOM_BTN_HEIGHT)
            zoom_in_btn.setToolTip("Zoom in")
            zoom_in_btn.clicked.connect(self.zoom_in)
            self._zoom_slider = QSlider(Qt.Orientation.Horizontal, zoom_strip)
            self._zoom_slider.setObjectName("annotationZoomSlider")
            self._zoom_slider.setRange(0, 100)
            self._zoom_slider.setValue(0)
            self._zoom_slider.setMinimumWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
            self._zoom_slider.setFixedHeight(ZOOM_SLIDER_HEIGHT)
            self._zoom_slider.setToolTip(
                "Left: fit full sequence in window. Right: maximum zoom-in."
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
                self._font_size,
                self._pos_start,
                self._pos_end,
                select_btn,
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
            controls_layout.addWidget(select_btn)

            font_label = QLabel("Font", annotation_controls)
            font_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
            controls_layout.addWidget(font_label)
            controls_layout.addWidget(self._font_size)

            self._sequence_stats_label.setAlignment(
                Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop
            )
            self._sequence_stats_label.setWordWrap(True)
            controls_layout.addWidget(self._sequence_stats_label)
            controls_layout.addStretch()
            self._annotation_side_panel = annotation_controls
        self._canvas = SequenceMapCanvas(self)
        self._canvas.set_color_bases(True)
        self._canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._canvas.selectionChanged.connect(self._on_selection_changed)
        self._canvas.featureClicked.connect(self._on_canvas_feature_clicked)
        self._canvas.featureActionRequested.connect(self._on_canvas_feature_action)

        self._ruler = ViewportRulerWidget(self)
        self._annot_track = AnnotationTrackWidget(self)
        self._annot_track.selectionChanged.connect(self._on_annot_track_selection_changed)
        self._annot_track.deleteRequested.connect(self._request_delete_features)
        self._annot_track.copySequenceRequested.connect(self._request_copy_feature_sequences)
        self._annot_track.extractSequenceRequested.connect(
            self._request_extract_feature_sequences
        )
        self._annot_track.modifyRequested.connect(self._request_modify_feature)
        self._annot_track.changeColorRequested.connect(self._request_change_feature_color)
        self._annot_track.tracksMoved.connect(self._on_annot_tracks_moved)
        self._annot_track.attachImageRequested.connect(self._request_attach_feature_image)
        self._annot_track.openImagesRequested.connect(self.openFeatureImagesRequested.emit)
        self._annot_track.imagesDropped.connect(self.featureImagesDropped.emit)
        self._overview = OverviewWidget(self)
        self._overview.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._overview.viewportChanged.connect(self._scroll_to_fraction)

        self._header = QWidget(self)
        self._header.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        header_layout = QVBoxLayout(self._header)
        header_layout.setContentsMargins(0, 0, 0, 0)
        header_layout.setSpacing(0)
        header_layout.addWidget(self._overview)
        header_layout.addWidget(self._ruler)

        self._canvas_row = QWidget(self)
        self._canvas_row.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        canvas_row_layout = QVBoxLayout(self._canvas_row)
        canvas_row_layout.setContentsMargins(0, 0, 0, 0)
        canvas_row_layout.setSpacing(0)
        canvas_row_layout.addWidget(self._canvas, 1)
        self._hbar = QScrollBar(Qt.Orientation.Horizontal, self._canvas_row)
        self._hbar.valueChanged.connect(self._on_horizontal_scroll)
        canvas_row_layout.addWidget(self._hbar)

        self._overview_sync_timer = QTimer(self)
        self._overview_sync_timer.setSingleShot(True)
        self._overview_sync_timer.setInterval(_OVERVIEW_SYNC_MS)
        self._overview_sync_timer.timeout.connect(self._sync_overview_now)

        self._layout_sync_timer = QTimer(self)
        self._layout_sync_timer.setSingleShot(True)
        self._layout_sync_timer.setInterval(_LAYOUT_SYNC_MS)
        self._layout_sync_timer.timeout.connect(self._sync_header_views)

        self._sequence_block = QWidget(self)
        self._sequence_block.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
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
        self._annot_scroll.viewport().setAutoFillBackground(True)

        block_layout = QVBoxLayout(self._sequence_block)
        block_layout.setContentsMargins(0, 0, 0, 0)
        block_layout.setSpacing(0)
        block_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        block_layout.addWidget(self._canvas_row, 0, Qt.AlignmentFlag.AlignTop)
        block_layout.addWidget(self._annot_scroll, 1)

        main_column = QWidget(self)
        main_column.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        main_layout = QVBoxLayout(main_column)
        main_layout.setContentsMargins(0, 0, 0, 0)
        main_layout.setSpacing(0)
        if self._toolbar is not None:
            main_layout.addWidget(self._toolbar)
        if zoom_row is not None:
            main_layout.addLayout(zoom_row)
        main_layout.addWidget(self._header)
        main_layout.addWidget(self._sequence_block, stretch=1)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(main_column, stretch=1)
        if self._annotation_side_panel is not None:
            layout.addWidget(self._annotation_side_panel)

        self._canvas.installEventFilter(self)
        self._canvas_row.installEventFilter(self)
        self._sequence_block.installEventFilter(self)
        self._annot_track.installEventFilter(self)
        self._annot_scroll.installEventFilter(self)
        self._annot_scroll.viewport().installEventFilter(self)

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
        if event.type() == QEvent.Type.Wheel and watched in (self._canvas, self._canvas_row):
            wheel = event
            if isinstance(wheel, QWheelEvent):
                delta_x = _horizontal_wheel_delta(
                    wheel,
                    scroll_step=max(1, self._hbar.singleStep() * 3),
                )
                if delta_x != 0 and self._hbar.maximum() > 0:
                    self._hbar.setValue(
                        max(0, min(self._hbar.maximum(), self._hbar.value() - delta_x))
                    )
                    wheel.accept()
                    return True
        return super().eventFilter(watched, event)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._layout_sync_timer.start()
        if self._fit_retries > 0:
            QTimer.singleShot(0, self.fit_window)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        QTimer.singleShot(0, self._sync_header_views)
        # First open often fits before the panel has its real width; re-check.
        QTimer.singleShot(0, self._ensure_fitted_to_viewport)

    def refresh_layout(self) -> None:
        self._sync_header_views()
        self._ensure_fitted_to_viewport()

    def _sync_header_views(self) -> None:
        self._sync_sequence_block_layout()
        self._canvas.sync_viewport_from_widget()
        self._update_horizontal_scroll()
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
        """Re-fit when the sequence strip is much narrower than the panel."""
        if self._data is None or not self._data.sequence:
            return
        if not self.isVisible():
            return
        viewport_w = self._scroll_viewport_width()
        panel_w = max(self.width(), self._sequence_block.width())
        if viewport_w <= 120 or (panel_w > 200 and viewport_w < panel_w * 0.45):
            self._fit_retries += 1
            if self._fit_retries < 40:
                QTimer.singleShot(50, self._ensure_fitted_to_viewport)
            return
        content_w = self._canvas.content_pixel_width()
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
        self._hbar.setVisible(max_scroll > 0)
        self._hbar.blockSignals(False)
        self._canvas.set_scroll_px(float(self._hbar.value()))

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
        self._overview.set_data(
            seq_len,
            self._data.features,
            view_start_bp=start_bp,
            view_end_bp=end_bp,
            selected_indices=self._selected_feature_indices,
        )

        self._ruler.set_view(
            sequence=self._data.sequence,
            scroll_px=scroll_val,
            ppb=cell_w,
            viewport_w=vw,
        )

        feature_rows, _ = _pack_feature_rows(self._data.features)
        self._annot_track.set_view(
            features=self._data.features,
            feature_rows=feature_rows,
            scroll_px=scroll_val,
            ppb=cell_w,
            viewport_w=vw,
            selected_indices=self._selected_feature_indices,
        )
        self._sync_sequence_block_layout()

    def _sync_sequence_block_layout(self) -> None:
        canvas_h = max(1, self._canvas.height())
        hbar_h = self._hbar.sizeHint().height() if self._hbar.isVisible() else 0
        row_h = canvas_h + hbar_h
        self._canvas_row.setMinimumHeight(row_h)
        self._canvas_row.setMaximumHeight(row_h)
        self._canvas_row.setFixedHeight(row_h)
        self._sync_annot_scroll_area()

    def _sync_annot_scroll_area(self) -> None:
        """Let annotation tracks fill space below the sequence; scroll when rows overflow."""
        self._annot_scroll.show()
        self._annot_scroll.setMinimumHeight(_ANNOT_MIN_VIEWPORT_HEIGHT)
        self._annot_scroll.setMaximumHeight(16777215)
        vw = max(1, self._scroll_viewport_width())
        track_h = max(1, self._annot_track.height())
        if self._annot_track.has_content():
            self._annot_track.setFixedSize(vw, track_h)
        else:
            # Empty placeholder sized to the viewport so the panel is not a void.
            vp_h = max(_ANNOT_MIN_VIEWPORT_HEIGHT, self._annot_scroll.viewport().height())
            self._annot_track.setFixedSize(vw, max(1, vp_h))

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

    def select_by_position(self, start: int, end: int, *, _attempt: int = 0) -> None:
        if self._data is None:
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

    def load(self, data: SequenceMapData, *, fit_after_load: bool = True) -> None:
        self._data = data
        self._dirty = False
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
        self._sync_position_spins(1, 1)
        seq_len = data.length
        if seq_len > _HUGE_SEQUENCE_LEN:
            self.statusMessage.emit(
                f"Very long sequence ({seq_len:,} bp): zoom and pan are limited; "
                "annotation tracks remain available."
            )
        elif seq_len > _INTERACTIVE_SEQUENCE_LEN:
            self.statusMessage.emit(
                f"Long sequence ({seq_len:,} bp): letter mode disabled for performance."
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

    def set_features(self, features: list[SequenceFeature]) -> None:
        if self._data is None:
            return
        self._data.features = list(features)
        self._canvas.set_sequence(
            self._data.sequence,
            self._data.features,
            is_protein=self._data.is_protein,
        )
        valid = set(range(len(self._data.features)))
        self._selected_feature_indices &= valid
        self._dirty = True
        self.featuresChanged.emit()
        self._sync_overview()

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
        if self._data is None or not self._data.sequence:
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
        self._canvas.clear_selection()
        self._sync_position_spins(1, 1)
        self._apply_zoom_ppb(ppb, sync_slider=True, reset_scroll=True)

    def fit_width(self) -> None:
        """Backward-compatible alias for fit_window."""
        self.fit_window()

    def refresh_nucleotide_colors(self) -> None:
        self._canvas.update()

    def selection_range(self) -> tuple[int, int] | None:
        return self._canvas.selection_range()

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
        feat = SequenceFeature(start=start, end=end, name=name, feature_type=feature_type)
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
        if start == end:
            self.statusMessage.emit(f"Position: {start}")
            return
        a, b = sorted((start, end))
        unit = "aa" if (self._data is not None and self._data.is_protein) else "bp"
        self.statusMessage.emit(f"Selection: {a} – {b} ({b - a + 1} {unit})")

    def set_selected_features(self, indices: set[int] | frozenset[int] | None) -> None:
        self._selected_feature_indices = set(indices or [])
        self._canvas.set_selected_features(self._selected_feature_indices)
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
        self._selected_feature_indices = selected
        self._canvas.set_selected_features(selected)
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
            feat = self._data.features[index]
            self._sync_position_spins(feat.start, feat.end)
            self._canvas.set_selection_range(feat.start, feat.end)
            self.statusMessage.emit(f"Feature: {feat.name} ({feat.start}–{feat.end})")

    def _sync_feature_selection_from_annot_track(self) -> None:
        annot_sel = self._annot_track.selected_indices()
        if not annot_sel:
            return
        if annot_sel == self._selected_feature_indices:
            return
        self._selected_feature_indices = set(annot_sel)
        self._canvas.set_selected_features(annot_sel)
        self._sync_overview()

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

    def _on_annot_tracks_moved(self, features: object) -> None:
        updated = list(features or [])
        if not updated:
            return
        self.annotationTracksMoveRequested.emit(updated)

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
        )


class AutoAnnotationSettingsPanel(QWidget):
    """Embeddable LTR / TIR / tandem / ORF auto-annotation parameters."""

    def __init__(
        self,
        settings: AutoAnnotationSettings,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._settings = settings

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        layout.addWidget(self._section_title("LTR"))
        layout.addLayout(self._terminal_repeat_form(settings.ltr or TerminalRepeatSettings(), "LTR"))
        layout.addWidget(_settings_separator())

        layout.addWidget(self._section_title("TIR"))
        layout.addLayout(
            self._terminal_repeat_form(
                settings.tir
                or TerminalRepeatSettings(
                    min_length=20,
                    max_length=500,
                    similarity_threshold=DEFAULT_TIR_SIMILARITY_THRESHOLD,
                ),
                "TIR",
            )
        )
        layout.addWidget(_settings_separator())

        layout.addWidget(self._section_title("Tandem repeat"))
        layout.addLayout(self._tandem_repeat_form(settings.tandem or TandemRepeatSettings()))
        layout.addWidget(_settings_separator())

        layout.addWidget(self._section_title("ORF"))
        layout.addLayout(self._orf_form(settings.orf or OrfSettings()))

        self._ltr_fields = self._terminal_fields[0]
        self._tir_fields = self._terminal_fields[1]
        self._tandem_fields = self._tandem_field_list
        self._orf_fields = self._orf_field_list

    def _section_title(self, text: str) -> QLabel:
        title = QLabel(text)
        font = title.font()
        font.setBold(True)
        title.setFont(font)
        return title

    def _bp_spin_row(self, spin: QSpinBox) -> QWidget:
        row = QWidget(self)
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        layout.addWidget(spin)
        layout.addWidget(QLabel("bp", row))
        layout.addStretch()
        return row

    def _terminal_repeat_form(self, cfg: TerminalRepeatSettings, repeat_kind: str) -> QFormLayout:
        form = QFormLayout()
        min_len = QSpinBox(self)
        min_len.setRange(1, 1_000_000)
        min_len.setValue(cfg.min_length)
        min_len.setToolTip(
            "Minimum repeat length (bp). Values below 11 also set the self-BLAST "
            "word size for this repeat type only (minimum word size 4; default 11)."
        )
        max_len = QSpinBox(self)
        max_len.setRange(1, 1_000_000)
        max_len.setValue(cfg.max_length)
        min_distance = QSpinBox(self)
        min_distance.setRange(0, 1_000_000)
        min_distance.setValue(cfg.min_pair_distance)
        similarity = QDoubleSpinBox(self)
        similarity.setRange(0.0, 1.0)
        similarity.setDecimals(2)
        similarity.setSingleStep(0.05)
        similarity.setValue(cfg.similarity_threshold)
        form.addRow("Min length:", self._bp_spin_row(min_len))
        form.addRow("Max length:", self._bp_spin_row(max_len))
        form.addRow(f"Min distance between {repeat_kind}:", self._bp_spin_row(min_distance))
        form.addRow("Similarity threshold:", similarity)
        if not hasattr(self, "_terminal_fields"):
            self._terminal_fields = []
        self._terminal_fields.append((min_len, max_len, min_distance, similarity))
        return form

    def _tandem_repeat_form(self, cfg: TandemRepeatSettings) -> QFormLayout:
        form = QFormLayout()

        fields: list[QSpinBox] = []
        for motif_size in (2, 3, 4, 5):
            spin = QSpinBox(self)
            spin.setRange(1, 1000)
            spin.setValue(getattr(cfg, f"min_repeat_{motif_size}"))
            form.addRow(f"Min repeat number ({motif_size}-motif):", spin)
            fields.append(spin)
        spin_6 = QSpinBox(self)
        spin_6.setRange(1, 1000)
        spin_6.setValue(cfg.min_repeat_6)
        form.addRow("Min repeat number (≥ 6-motif):", spin_6)
        fields.append(spin_6)
        self._tandem_field_list = fields

        approximate_cb = QCheckBox("Enable approximate tandem repeat", self)
        approximate_cb.setChecked(cfg.enable_approximate)
        approximate_cb.setToolTip(
            "Run pytrf ATRFinder to detect imperfect repeats with mismatches. "
            "ATR runs before exact STR/GTR so long approximate blocks are kept."
        )
        form.addRow(approximate_cb)
        self._tandem_approximate_cb = approximate_cb

        atr_identity = QDoubleSpinBox(self)
        atr_identity.setRange(0.0, 1.0)
        atr_identity.setDecimals(2)
        atr_identity.setSingleStep(0.05)
        atr_identity.setValue(cfg.atr_min_identity)
        atr_identity.setToolTip(
            "Minimum alignment identity (0–1). 0.75 = 75% similarity required while extending."
        )
        atr_min_motif = QSpinBox(self)
        atr_min_motif.setRange(1, 100)
        atr_min_motif.setValue(cfg.atr_min_motif)
        atr_min_motif.setToolTip("Minimum motif length for approximate repeat seeds.")
        atr_max_errors = QSpinBox(self)
        atr_max_errors.setRange(0, 50)
        atr_max_errors.setValue(cfg.atr_max_errors)
        atr_max_errors.setToolTip(
            "Maximum consecutive mismatches allowed while extending an imperfect repeat."
        )

        form.addRow("Approx. min identity:", atr_identity)
        form.addRow("Approx. min motif:", self._bp_spin_row(atr_min_motif))
        form.addRow("Approx. max errors:", atr_max_errors)
        self._tandem_approximate_fields = (atr_identity, atr_min_motif, atr_max_errors)
        approximate_cb.toggled.connect(self._sync_tandem_approximate_fields)
        self._sync_tandem_approximate_fields(cfg.enable_approximate)

        return form

    def _sync_tandem_approximate_fields(self, enabled: bool) -> None:
        for widget in getattr(self, "_tandem_approximate_fields", ()):
            widget.setEnabled(enabled)

    def _orf_form(self, cfg: OrfSettings) -> QFormLayout:
        form = QFormLayout()
        min_len = QSpinBox(self)
        min_len.setRange(3, 1_000_000)
        min_len.setValue(cfg.min_length)
        max_len = QSpinBox(self)
        max_len.setRange(3, 1_000_000)
        max_len.setValue(cfg.max_length)
        alt_starts = QLineEdit(self)
        alt_starts.setPlaceholderText("e.g. TTG, CTG")
        alt_starts.setText(cfg.alternative_start_codons)
        alt_starts.setToolTip("Optional extra start codons (comma- or space-separated).")
        form.addRow("Min ORF length:", self._bp_spin_row(min_len))
        form.addRow("Max ORF length:", self._bp_spin_row(max_len))
        form.addRow("Add other start codons:", alt_starts)
        self._orf_field_list = (min_len, max_len, alt_starts)
        return form

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
        )
        return AutoAnnotationSettings(ltr=ltr, tir=tir, tandem=tandem, orf=orf)

    @staticmethod
    def _read_terminal_settings(
        fields: tuple[QSpinBox, QSpinBox, QSpinBox, QDoubleSpinBox],
        *,
        default_min: int,
    ) -> TerminalRepeatSettings:
        min_len, max_len, min_distance, similarity = fields
        return TerminalRepeatSettings(
            min_length=min_len.value() or default_min,
            max_length=max(max_len.value(), min_len.value()),
            min_pair_distance=min_distance.value(),
            similarity_threshold=similarity.value(),
        )


class AutoAnnotationSettingsDialog(QDialog):
    """Configure LTR, TIR, tandem-repeat, and ORF auto-annotation parameters."""

    def __init__(
        self,
        settings: AutoAnnotationSettings,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Auto annotation settings")
        self._panel = AutoAnnotationSettingsPanel(settings, self)
        layout = QVBoxLayout(self)
        layout.setSpacing(8)
        layout.addWidget(self._panel)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def settings(self) -> AutoAnnotationSettings:
        return self._panel.settings()


class SimilarityAnnotationSettingsPanel(QWidget):
    """Embeddable similarity-annotation search thresholds."""

    def __init__(
        self,
        settings: SimilarityAnnotationSettings,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)

        intro = QLabel(
            "These E-value and output limits apply only to Similarity annotation. "
            "They do not change TE-Aid, genome BLAST, or other tools in RepeatViewer.",
            self,
        )
        intro.setWordWrap(True)
        layout.addWidget(intro)
        layout.addWidget(_settings_separator())

        form = QFormLayout()
        self._blast_evalue = QDoubleSpinBox(self)
        configure_report_evalue_spinbox(self._blast_evalue)
        self._blast_evalue.setValue(settings.blast_evalue)
        self._blast_evalue.setToolTip(
            "Reporting E-value for BLASTP and BLASTN similarity searches."
        )
        form.addRow("BLAST E-value (BLASTP & BLASTN):", self._blast_evalue)

        self._cdd_evalue = QDoubleSpinBox(self)
        configure_report_evalue_spinbox(self._cdd_evalue)
        self._cdd_evalue.setValue(settings.cdd_evalue)
        self._cdd_evalue.setToolTip("E-value cutoff for CDD rpstblastn similarity searches.")
        form.addRow("CDD E-value (rpstblastn):", self._cdd_evalue)

        self._hmm_evalue = QDoubleSpinBox(self)
        configure_report_evalue_spinbox(self._hmm_evalue)
        self._hmm_evalue.setValue(settings.hmm_evalue)
        self._hmm_evalue.setToolTip("Reporting E-value for pyhmmer HMM similarity searches.")
        form.addRow("HMM E-value (pyhmmer):", self._hmm_evalue)

        self._max_hits = QSpinBox(self)
        self._max_hits.setRange(1, 100)
        self._max_hits.setValue(settings.max_hits_per_database)
        self._max_hits.setToolTip(
            "Maximum total hits shown per database. When more hits pass filtering, "
            f"{APP_DISPLAY_NAME} first removes shorter fully nested overlaps when their E-value "
            "is not better than the covering hit (including equal E=0), then merges "
            "partially overlapping groups (keeping the highest score in each group), "
            "then keeps the top hits by score. Protein and DNA searches use the same rule."
        )
        form.addRow("Maximum hits per database:", self._max_hits)

        self._show_db_name = QCheckBox(self)
        self._show_db_name.setChecked(settings.show_database_name_on_track)
        self._show_db_name.setToolTip(
            "When checked, similarity hits are labeled as “Database: target”. "
            "When unchecked, only the target or domain name is shown."
        )
        form.addRow("Show database name on feature track:", self._show_db_name)

        self._show_info_dialog = QCheckBox(self)
        self._show_info_dialog.setChecked(settings.show_info_window_after_run)
        self._show_info_dialog.setToolTip(
            f"When checked, {APP_DISPLAY_NAME} shows an information dialog with the similarity search summary "
            "(e.g. hits/no hits or database-readiness problems)."
        )
        form.addRow("Show info window after run:", self._show_info_dialog)
        layout.addLayout(form)

    def settings(self) -> SimilarityAnnotationSettings:
        return SimilarityAnnotationSettings(
            blast_evalue=self._blast_evalue.value(),
            hmm_evalue=self._hmm_evalue.value(),
            cdd_evalue=self._cdd_evalue.value(),
            max_hits_per_database=self._max_hits.value(),
            show_database_name_on_track=self._show_db_name.isChecked(),
            show_info_window_after_run=self._show_info_dialog.isChecked(),
        )


class SimilarityAnnotationSettingsDialog(QDialog):
    """Configure similarity-annotation search thresholds."""

    def __init__(
        self,
        settings: SimilarityAnnotationSettings,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Similarity annotation settings")
        self._panel = SimilarityAnnotationSettingsPanel(settings, self)
        layout = QVBoxLayout(self)
        layout.setSpacing(8)
        layout.addWidget(self._panel)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def settings(self) -> SimilarityAnnotationSettings:
        return self._panel.settings()


class SequenceMapPanel(QWidget):
    """Embedded sequence map with annotation toolbar (preview Annotation tab)."""

    statusMessage = Signal(str)
    annotationJobProgress = Signal(str)
    annotationsSaved = Signal(object)  # Path to the written GFF3
    sequenceExtracted = Signal(object)  # Path to a newly written FASTA
    activeSequenceChanged = Signal(str)  # sequence id after multi-seq switch

    def __init__(
        self,
        fasta_path: Path,
        parent: QWidget | None = None,
        *,
        autoload: bool = False,
    ) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._path = fasta_path.expanduser().resolve()
        # Sequence from the last load_path(), so the summary can be computed
        # without re-reading and re-parsing the FASTA from disk.
        self._loaded_data = None
        self._map = SequenceMapWidget(self, embedded=True)
        self._map.statusMessage.connect(self.statusMessage.emit)
        self._map._canvas.addFeatureRequested.connect(
            lambda ft, bl: self._add_feature_of_type(ft, feature_label=bl)
        )
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
        self._map.attachFeatureImageRequested.connect(self._attach_image_to_selected_feature)
        self._map.openFeatureImagesRequested.connect(self._open_feature_images)
        self._map.featureImagesDropped.connect(self._on_feature_images_dropped)
        self._map._annot_track.set_image_marker_lookup(self._marked_feature_indices)
        from repeatviewer.widgets.cyto_images_ui import cyto_notifier

        cyto_notifier().changed.connect(self._on_cyto_images_changed)
        self.annotationJobProgress.connect(self._set_annotation_job_status)
        self._auto_ann_thread: QThread | None = None
        self._auto_ann_worker: Worker | None = None
        self._auto_ann_cancel = JobCancellation()
        self._sim_ann_thread: QThread | None = None
        self._sim_ann_worker: Worker | None = None
        self._sim_ann_cancel = JobCancellation()
        self._similarity_db_rows: list[tuple[UserDatabase, QCheckBox]] = []
        self._similarity_db_empty: QLabel | None = None
        self._auto_ann_settings = load_auto_annotation_settings()
        self._similarity_ann_settings = load_similarity_annotation_settings()
        self._feature_undo_stack: list[list[SequenceFeature]] = []

        controls = QWidget(self)
        self._controls_panel = controls
        _configure_annotation_controls_panel(controls)
        controls.setMinimumWidth(0)
        controls.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Minimum,
        )
        controls_layout = QVBoxLayout(controls)
        controls_layout.setContentsMargins(0, 0, 0, 0)
        controls_layout.setSpacing(0)

        self._auto_all_updating = False
        self._similarity_all_updating = False
        self._add_feature_buttons: list[_AddFeatureControl] = []
        self._add_feature_button_defs: list[tuple[str, str]] = []
        self._sequence_switch_updating = False

        self._seq_row_container, seq_row = _control_band(controls)
        self._seq_label = _annotation_row_label("Sequence:", self._seq_row_container)
        seq_row.addWidget(self._seq_label, 0, Qt.AlignmentFlag.AlignVCenter)
        self._sequence_combo = QComboBox(self._seq_row_container)
        self._sequence_combo.setEditable(False)
        self._sequence_combo.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        self._sequence_combo.setMinimumWidth(180)
        self._sequence_combo.setToolTip(
            "Current sequence shown on the map. Use Search to jump by name."
        )
        self._sequence_combo.currentIndexChanged.connect(self._on_sequence_combo_changed)
        seq_row.addWidget(self._sequence_combo, 1, Qt.AlignmentFlag.AlignVCenter)

        self._seq_search_label = QLabel("Search:", self._seq_row_container)
        self._seq_search_label.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        seq_row.addWidget(self._seq_search_label, 0, Qt.AlignmentFlag.AlignVCenter)
        self._sequence_search = QLineEdit(self._seq_row_container)
        self._sequence_search.setPlaceholderText("Filter by sequence name…")
        self._sequence_search.setClearButtonEnabled(True)
        self._sequence_search.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        self._sequence_search.setMinimumWidth(140)
        self._sequence_search.setToolTip(
            "Type part of a sequence id, then press Enter or pick a suggestion."
        )
        self._sequence_search.returnPressed.connect(self._on_sequence_search_submitted)
        seq_row.addWidget(self._sequence_search, 1, Qt.AlignmentFlag.AlignVCenter)

        self._seq_count_label = QLabel("", self._seq_row_container)
        self._seq_count_label.setStyleSheet("color: palette(mid);")
        seq_row.addWidget(self._seq_count_label, 0, Qt.AlignmentFlag.AlignVCenter)
        controls_layout.addWidget(self._seq_row_container)
        self._seq_row_container.hide()
        self._sequence_search_completer: QCompleter | None = None

        self._add_row_container, add_row = _control_band(controls)
        self._add_feature_label = _annotation_row_label("Add feature:", self._add_row_container)
        add_row.addWidget(self._add_feature_label, 0, Qt.AlignmentFlag.AlignVCenter)
        self._add_feature_buttons_row = _AddFeatureButtonsRow(self._add_row_container)
        self._add_feature_overflow = QComboBox(self._add_row_container)
        self._add_feature_overflow.setMinimumWidth(110)
        self._add_feature_overflow.activated.connect(self._on_add_feature_overflow_activated)
        self._add_feature_overflow.hide()
        add_row.addWidget(self._add_feature_buttons_row, 0, Qt.AlignmentFlag.AlignVCenter)
        add_row.addWidget(self._add_feature_overflow, 0, Qt.AlignmentFlag.AlignVCenter)
        add_row.addStretch(1)
        self._add_feature_settings_btn = QPushButton("Setting", self._add_row_container)
        self._add_feature_settings_btn.setToolTip(
            "Create, delete, or reorder Add feature buttons"
        )
        self._add_feature_settings_btn.clicked.connect(self._open_add_feature_buttons_settings)
        self._add_run_settings_box = _annotation_run_settings_box(
            None,
            self._add_feature_settings_btn,
            self._add_row_container,
        )
        add_row.addWidget(self._add_run_settings_box, 0, Qt.AlignmentFlag.AlignVCenter)
        controls_layout.addWidget(self._add_row_container)
        self._rebuild_add_feature_buttons(load_add_feature_buttons())

        self._auto_row_container, auto_row = _control_band(controls)
        self._auto_ann_label = _annotation_row_label("Auto annotation:", self._auto_row_container)
        auto_row.addWidget(self._auto_ann_label, 0, Qt.AlignmentFlag.AlignVCenter)
        self._auto_feature_columns = _FeatureColumnRow(
            self._auto_row_container,
            column_count=_AUTO_ANNOTATION_COLUMN_COUNT,
        )
        self._auto_ann_checks: dict[AutoAnnotationMode, QCheckBox] = {}
        for label, mode in _AUTO_ANNOTATION_MODES:
            check = QCheckBox(label, self._auto_feature_columns)
            check.setChecked(True)
            check.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
            check.toggled.connect(self._on_auto_annotation_check_toggled)
            self._auto_ann_checks[mode] = check
            self._auto_feature_columns.place(_AUTO_CHECK_COLUMN[mode], check)
        auto_row.addWidget(self._auto_feature_columns, 0, Qt.AlignmentFlag.AlignVCenter)
        self._auto_ann_all_check = QCheckBox("All", self._auto_row_container)
        self._auto_ann_all_check.setToolTip("Select or clear all auto-annotation detectors")
        self._auto_ann_all_check.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        self._auto_ann_all_check.toggled.connect(self._on_auto_annotation_all_toggled)
        auto_row.addWidget(self._auto_ann_all_check, 0, Qt.AlignmentFlag.AlignVCenter)
        auto_row.addStretch(1)
        self._auto_ann_run_btn = QPushButton("Run", self._auto_row_container)
        self._auto_ann_run_btn.setToolTip(
            "Run the selected auto-annotation detectors on this sequence"
        )
        self._auto_ann_run_btn.clicked.connect(self._on_auto_run_clicked)
        self._auto_ann_settings_btn = QPushButton("Setting", self._auto_row_container)
        self._auto_ann_settings_btn.setToolTip("Configure LTR, TIR, tandem, and ORF detection")
        self._auto_ann_settings_btn.clicked.connect(self._open_auto_annotation_settings)
        self._auto_run_settings_box = _annotation_run_settings_box(
            self._auto_ann_run_btn,
            self._auto_ann_settings_btn,
            self._auto_row_container,
        )
        auto_row.addWidget(self._auto_run_settings_box, 0, Qt.AlignmentFlag.AlignVCenter)
        controls_layout.addWidget(self._auto_row_container)
        self._sync_auto_annotation_all_check()

        self._sim_row_container = QWidget(controls)
        self._sim_row_container.setObjectName("annotationControlBand")
        self._sim_row_container.setMinimumWidth(0)
        self._sim_row_container.setFixedHeight(_CONTROL_BAND_HEIGHT)
        self._sim_row_container.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        sim_row = QHBoxLayout(self._sim_row_container)
        sim_row.setContentsMargins(8, 0, 8, 0)
        sim_row.setSpacing(6)
        sim_row.setAlignment(Qt.AlignmentFlag.AlignVCenter)

        homology_title = _annotation_row_label("Similarity anno:", self._sim_row_container)
        homology_title.setToolTip(
            "Search user-defined databases; each enabled database uses its own color"
        )
        sim_row.addWidget(homology_title, 0, Qt.AlignmentFlag.AlignVCenter)

        self._similarity_db_flow = _SimilarityDbFlowWidget(self._sim_row_container)
        sim_row.addWidget(self._similarity_db_flow, 1)

        sim_row.addStretch(1)

        self._similarity_ann_all_check = QCheckBox("All", self._sim_row_container)
        self._similarity_ann_all_check.setToolTip("Select or clear all similarity databases")
        self._similarity_ann_all_check.setSizePolicy(
            QSizePolicy.Policy.Maximum,
            QSizePolicy.Policy.Fixed,
        )
        self._similarity_ann_all_check.toggled.connect(self._on_similarity_annotation_all_toggled)

        self._similarity_ann_run_btn = QPushButton("Run", self._sim_row_container)
        self._similarity_ann_run_btn.setToolTip(
            "Run similarity annotation with the selected databases"
        )
        self._similarity_ann_run_btn.clicked.connect(self._on_similarity_run_clicked)
        self._similarity_ann_settings_btn = QPushButton("Setting", self._sim_row_container)
        self._similarity_ann_settings_btn.setToolTip(
            "Configure E-value cutoffs and similarity-annotation search options"
        )
        self._similarity_ann_settings_btn.clicked.connect(self._open_similarity_annotation_settings)
        self._similarity_run_settings_box = _annotation_run_settings_box(
            self._similarity_ann_run_btn,
            self._similarity_ann_settings_btn,
            self._sim_row_container,
        )
        sim_row.addWidget(self._similarity_run_settings_box, 0, Qt.AlignmentFlag.AlignVCenter)
        controls_layout.addWidget(self._sim_row_container)
        self._refresh_similarity_database_rows()

        self._action_row_container, action_row = _control_band(controls)
        self._action_row_container.setObjectName("annotationActionBand")
        action_row.setSpacing(6)
        self._fit_btn = QPushButton("Fit window", self._action_row_container)
        self._fit_btn.setToolTip(
            "Scale the sequence so the full length fits in the visible width"
        )
        self._fit_btn.clicked.connect(self._map.fit_window)
        action_row.addWidget(self._fit_btn, 0, Qt.AlignmentFlag.AlignVCenter)
        del_btn = QPushButton("Delete", self._action_row_container)
        del_btn.setToolTip(
            "Delete selected feature(s). Select on the track below, then press Delete or Backspace."
        )
        del_btn.clicked.connect(lambda: self._delete_feature(confirm_empty=True))
        restore_btn = QPushButton("Restore", self._action_row_container)
        restore_btn.setToolTip(
            "Undo the last feature change on this sequence (delete, modify, color, move, …). "
            "History is cleared when you open another sequence."
        )
        restore_btn.clicked.connect(self._restore_features)
        color_btn = QPushButton("Color", self._action_row_container)
        color_btn.setToolTip(
            "Set colors for feature types, similarity databases, and nucleotides (A/C/G/T)"
        )
        color_btn.clicked.connect(self._open_element_colors_dialog)
        export_btn = QPushButton("Export", self._action_row_container)
        export_btn.setToolTip(
            "Export the current annotation view as SVG / PNG / PDF"
        )
        export_btn.clicked.connect(self._export_annotation_view)
        save_btn = QPushButton("Save", self._action_row_container)
        save_btn.setToolTip("Write annotation changes to the JSON sidecar and GFF3 file")
        save_btn.clicked.connect(lambda: self._save_annotations(silent=True))
        for btn in (del_btn, restore_btn, color_btn, export_btn, save_btn):
            action_row.addWidget(btn, 0, Qt.AlignmentFlag.AlignVCenter)
        action_row.addStretch()
        self._annotation_job_status = QLabel("", self._action_row_container)
        self._annotation_job_status.setObjectName("annotationJobStatus")
        self._annotation_job_status.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        self._annotation_job_status.setMinimumWidth(180)
        action_row.addWidget(self._annotation_job_status, 0, Qt.AlignmentFlag.AlignVCenter)
        controls_layout.addWidget(self._action_row_container)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._map, stretch=1)
        layout.addWidget(_annotation_panel_separator(self))
        layout.addWidget(controls)
        controls.installEventFilter(self)

        if autoload:
            self.load_path(self._path)

        self._map.set_show_similarity_database_name(
            self._similarity_ann_settings.show_database_name_on_track
        )
        self._install_annotation_shortcuts()
        QTimer.singleShot(0, self._sync_annotation_control_rows)

    def bind_annotation_zoom_bar(self, bar) -> None:
        self._map.bind_annotation_zoom_bar(bar)

    def showEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().showEvent(event)
        QTimer.singleShot(0, self._sync_annotation_control_rows)
        QTimer.singleShot(0, self._refresh_annotation_button_styles)
        # Rebuild DB chips only when empty — every Show used to recreate them
        # and made Text ↔ Annotation ↔ Dotplot switches feel sluggish.
        if not self._similarity_db_rows:
            QTimer.singleShot(0, self._refresh_similarity_database_rows)

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
        if watched is getattr(self, "_controls_panel", None) and event.type() in (
            QEvent.Type.FocusIn,
            QEvent.Type.FocusOut,
            QEvent.Type.Show,
            QEvent.Type.Resize,
        ):
            QTimer.singleShot(0, self._sync_annotation_control_rows)
            QTimer.singleShot(0, self._refresh_annotation_button_styles)
        return super().eventFilter(watched, event)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        QTimer.singleShot(0, self._sync_annotation_control_rows)

    def _sync_annotation_control_rows(self) -> None:
        self._sync_add_feature_row()
        self._sync_auto_annotation_row_alignment()
        self._sync_run_settings_boxes()
        self._sync_similarity_row_layout()

    def _sync_run_settings_boxes(self) -> None:
        if not hasattr(self, "_auto_run_settings_box") or not hasattr(
            self, "_similarity_run_settings_box"
        ):
            return
        run_width = max(
            self._auto_ann_run_btn.sizeHint().width(),
            self._similarity_ann_run_btn.sizeHint().width(),
        )
        add_spacer = getattr(self._add_run_settings_box, "_run_spacer", None)
        if isinstance(add_spacer, QWidget):
            add_spacer.setFixedWidth(max(1, run_width))
        boxes = (
            self._add_run_settings_box,
            self._auto_run_settings_box,
            self._similarity_run_settings_box,
        )
        width = max(
            _ANNOTATION_RUN_SETTINGS_WIDTH,
            *(box.sizeHint().width() for box in boxes),
        )
        for box in boxes:
            box.setFixedWidth(width)

    def _sync_similarity_row_layout(self) -> None:
        if not hasattr(self, "_similarity_db_flow"):
            return
        self._similarity_db_flow._reflow()
        self._sim_row_container.setFixedHeight(_CONTROL_BAND_HEIGHT)

    def _sync_add_feature_row(self) -> None:
        if not self._add_feature_buttons:
            return
        self._add_feature_buttons_row._update_size()
        available = max(
            self._add_row_container.width(),
            self._auto_row_container.width(),
            max(0, self.width() - 24),
        )
        label_width = _ANNOTATION_ROW_LABEL_WIDTH + 12
        trailing = self._add_run_settings_box.width() + 8
        budget = max(120, available - label_width - trailing)
        buttons_width = self._add_feature_buttons_row.content_width()
        if buttons_width <= budget:
            self._add_feature_buttons_row.show()
            for control in self._add_feature_buttons:
                control.show()
            self._add_feature_overflow.hide()
        else:
            self._add_feature_buttons_row.hide()
            for control in self._add_feature_buttons:
                control.hide()
            self._add_feature_overflow.show()
        self._sync_auto_feature_columns()
        self._add_feature_buttons_row.updateGeometry()
        self._add_row_container.updateGeometry()
        self._auto_row_container.updateGeometry()

    def _sync_auto_feature_columns(self) -> None:
        if not hasattr(self, "_auto_feature_columns"):
            return
        column_count = max(_AUTO_CHECK_COLUMN.values(), default=0) + 1
        widths = [0] * column_count
        for mode, check in self._auto_ann_checks.items():
            column = _AUTO_CHECK_COLUMN[mode]
            widths[column] = max(widths[column], _column_content_width(check))
            check.show()
        self._auto_feature_columns.show()
        self._auto_feature_columns.sync_column_widths(widths)
        self._auto_feature_columns.updateGeometry()

    def _rebuild_add_feature_buttons(
        self, buttons: list[tuple[str, str]] | None = None
    ) -> None:
        defs = list(buttons if buttons is not None else load_add_feature_buttons())
        if not defs:
            defs = list(DEFAULT_ADD_FEATURE_BUTTONS)
        self._add_feature_button_defs = defs
        self._add_feature_buttons.clear()
        self._add_feature_overflow.blockSignals(True)
        self._add_feature_overflow.clear()
        rebuilt: list[_AddFeatureControl] = []
        for label, feature_type in defs:
            control = _AddFeatureControl(label, self._add_feature_buttons_row)
            control.setToolTip(f"Add a {label} feature for the selected region")
            control.clicked.connect(
                lambda _checked=False, ft=feature_type, bl=label: self._add_feature_of_type(
                    ft, feature_label=bl
                )
            )
            rebuilt.append(control)
            self._add_feature_overflow.addItem(label, feature_type)
        self._add_feature_buttons = rebuilt
        self._add_feature_buttons_row.set_controls(rebuilt)
        self._add_feature_overflow.blockSignals(False)
        if hasattr(self, "_auto_feature_columns"):
            self._sync_add_feature_row()

    def _open_add_feature_buttons_settings(self) -> None:
        dialog = AddFeatureButtonsDialog(parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._rebuild_add_feature_buttons(dialog.buttons())
        self._broadcast_shared_annotation_settings_changed()
        self.statusMessage.emit("Add feature buttons updated.")

    def _sync_auto_annotation_row_alignment(self) -> None:
        self._sync_add_feature_row()

    def _on_auto_annotation_check_toggled(self, _checked: bool) -> None:
        if self._auto_all_updating:
            return
        self._sync_auto_annotation_all_check()

    def _sync_auto_annotation_all_check(self) -> None:
        checks = list(self._auto_ann_checks.values())
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
        for check in self._auto_ann_checks.values():
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

    def _open_element_colors_dialog(self) -> None:
        dialog = AnnotationElementColorsDialog(parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        if dialog.apply_to_current_sequence():
            self._apply_element_colors_to_map(
                dialog.old_feature_type_colors(),
                dialog.new_feature_type_colors(),
                dialog.old_database_colors(),
                dialog.new_database_colors(),
            )
        self._refresh_similarity_database_rows()
        self._map.refresh_nucleotide_colors()
        self._map.update()
        self._broadcast_shared_annotation_settings_changed()
        self.statusMessage.emit("Element colors updated.")

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

    def fasta_path(self) -> Path:
        return self._path

    def is_dirty(self) -> bool:
        return self._map.is_dirty()

    def load_path(self, path: Path, *, fit_after_load: bool = True) -> None:
        self._path = path.expanduser().resolve()
        try:
            data = load_sequence_map_data(self._path)
        except (OSError, ValueError) as exc:
            self._loaded_data = None
            self._populate_sequence_combo(None)
            self.statusMessage.emit(str(exc))
            return
        self._loaded_data = data
        self._populate_sequence_combo(data)
        self._map.load(data, fit_after_load=fit_after_load)
        self._feature_undo_stack.clear()
        self.refresh_sequence_summary()

    def _populate_sequence_combo(self, data: SequenceMapData | None) -> None:
        self._sequence_switch_updating = True
        try:
            self._sequence_combo.clear()
            self._sequence_search.clear()
            if data is None or not data.record_ids:
                self._seq_row_container.hide()
                self._seq_count_label.clear()
                self._sequence_search_completer = None
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
                self._seq_row_container.show()
            else:
                self._seq_count_label.clear()
                self._seq_row_container.hide()
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
        self._feature_undo_stack.clear()
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
    def _on_extract_sequence(self, start: int, end: int) -> None:
        """Write the selected range out as a new single-record FASTA."""
        data = self._loaded_data
        if data is None or not getattr(data, "sequence", ""):
            QMessageBox.warning(
                self, "Extract Sequence", "The sequence is not loaded yet."
            )
            return

        sequence = data.sequence
        dialog = _ExtractSequenceDialog(
            sequence_id=data.sequence_id,
            sequence_length=len(sequence),
            start=min(int(start), int(end)),
            end=max(int(start), int(end)),
            default_name=self._next_extraction_name(data.sequence_id),
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
        self._rebuild_add_feature_buttons()
        self.refresh_similarity_databases()
        self.refresh_nucleotide_colors()
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

    def _annotation_changes_pending(self) -> bool:
        if self._map.is_dirty():
            return True
        if not self._map.features():
            return False
        return resolve_features_gff_path(self._path) is None

    def _save_annotations(self, *, silent: bool = False) -> bool:
        if self._map.data() is None:
            return False
        if not self._annotation_changes_pending():
            if not silent:
                self.statusMessage.emit("No annotation changes to save.")
            return True
        if not self._persist_annotations(show_message=not silent):
            return False
        if not silent:
            gff_path = features_gff_path(self._path)
            if gff_path.is_file():
                QMessageBox.information(
                    self,
                    "Saved",
                    f"Annotations saved to:\n{gff_path}",
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
                _json_path, gff_path = save_annotations_multi(
                    self._path,
                    features_by_id=data.features_by_id,
                    lengths_by_id=lengths,
                )
            else:
                _json_path, gff_path = save_annotations(
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
        if gff_path is not None:
            if show_message:
                self.statusMessage.emit(f"Saved annotations: {gff_path}")
            self.annotationsSaved.emit(gff_path)
        elif show_message:
            self.statusMessage.emit("Saved annotations (JSON sidecar).")
        return True

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

    def _set_annotation_job_status(self, text: str) -> None:
        self._annotation_job_status.setText(text.strip())

    def _clear_annotation_job_status(self) -> None:
        self._annotation_job_status.clear()

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
        if running:
            self._auto_ann_run_btn.setText("Stop")
            self._auto_ann_run_btn.setToolTip("Stop the running auto-annotation job")
        else:
            self._auto_ann_run_btn.setText("Run")
            self._auto_ann_run_btn.setToolTip(
                "Run the selected auto-annotation detectors on this sequence"
            )

    def _set_similarity_run_button_running(self, running: bool) -> None:
        if running:
            self._similarity_ann_run_btn.setText("Stop")
            self._similarity_ann_run_btn.setToolTip("Stop the running similarity-annotation job")
        else:
            self._similarity_ann_run_btn.setText("Run")
            self._similarity_ann_run_btn.setToolTip(
                "Run similarity annotation with the selected databases"
            )

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

    def _stop_auto_annotation(self) -> None:
        if not self._auto_annotate_running():
            return
        self._log_annotation_job("JOB STOP", "Auto annotation")
        self._auto_ann_cancel.request()
        self._set_annotation_job_status("Stopping auto annotation…")

    def _stop_similarity_annotation(self) -> None:
        if not self._similarity_annotation_running():
            return
        self._log_annotation_job("JOB STOP", "Similarity annotation")
        self._sim_ann_cancel.request()
        self._set_annotation_job_status("Stopping similarity annotation…")

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
        self._broadcast_shared_annotation_settings_changed()
        self.statusMessage.emit("Auto annotation settings saved.")

    def _selected_auto_annotation_modes(self) -> list[AutoAnnotationMode]:
        return [mode for mode, check in self._auto_ann_checks.items() if check.isChecked()]

    def _open_similarity_annotation_settings(self) -> None:
        self._similarity_ann_settings = load_similarity_annotation_settings()
        dialog = SimilarityAnnotationSettingsDialog(self._similarity_ann_settings, parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._similarity_ann_settings = dialog.settings()
        save_similarity_annotation_settings(self._similarity_ann_settings)
        self._map.set_show_similarity_database_name(
            self._similarity_ann_settings.show_database_name_on_track
        )
        self._broadcast_shared_annotation_settings_changed()
        self.statusMessage.emit("Similarity annotation settings saved.")

    def _current_similarity_annotation_settings(self) -> SimilarityAnnotationSettings:
        return self._similarity_ann_settings

    def _reflow_similarity_database_chips(self) -> None:
        if hasattr(self, "_similarity_db_flow"):
            self._similarity_db_flow._reflow()

    def _make_similarity_db_check(self, entry: UserDatabase, index: int) -> QCheckBox:
        color = similarity_database_color(entry.name, index)
        check = QCheckBox(entry.name, self._similarity_db_flow)
        readiness = database_search_readiness(entry)
        tooltip = f"{entry.display_type_label}\n{entry.path}"
        if readiness is not None:
            tooltip += f"\n\nNot ready: {readiness}"
            check.setChecked(False)
        else:
            check.setChecked(True)
        check.setToolTip(tooltip)
        check.setStyleSheet(f"QCheckBox {{ color: {color}; }}")
        check.setSizePolicy(
            QSizePolicy.Policy.Maximum,
            QSizePolicy.Policy.Fixed,
        )
        check.toggled.connect(self._on_similarity_db_check_toggled)
        return check

    def _refresh_similarity_database_rows(self) -> None:
        if not hasattr(self, "_similarity_db_flow"):
            return
        self._similarity_db_rows.clear()
        if self._similarity_db_empty is not None:
            self._similarity_db_empty = None

        databases = load_similarity_annotation_databases()
        if not databases:
            self._similarity_db_empty = QLabel(
                "No databases defined. Add them in File → Define database.",
                self._similarity_db_flow,
            )
            self._similarity_db_empty.setWordWrap(False)
            self._similarity_db_empty.setAlignment(
                Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
            )
            self._similarity_db_flow.set_chips([self._similarity_db_empty])
            self._similarity_ann_all_check.setParent(self._sim_row_container)
            self._similarity_ann_all_check.hide()
            self._sync_similarity_annotation_all_check()
            return

        chips: list[QWidget] = []
        for index, entry in enumerate(databases):
            check = self._make_similarity_db_check(entry, index)
            chips.append(check)
            self._similarity_db_rows.append((entry, check))
        chips.append(self._similarity_ann_all_check)
        self._similarity_db_flow.set_chips(chips)
        self._similarity_ann_all_check.show()
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
        self._set_similarity_annotation_enabled(False)
        labels = ", ".join(run.database.name for run in databases)
        self._log_annotation_job(
            "JOB START",
            "Similarity annotation",
            detail=f"databases=[{labels}]",
        )
        self._set_annotation_job_status(f"Running similarity annotation: {labels}…")
        self.statusMessage.emit(f"Running similarity annotation: {labels}…")

        worker = Worker(
            run_similarity_annotations,
            self._path,
            tuple(databases),
            settings=self._current_similarity_annotation_settings(),
            use_system_blast=load_use_system_blast_from_settings(),
            orf_settings=self._auto_ann_settings.orf,
            cancel=self._sim_ann_cancel,
            base_features=tuple(self._map.features()),
            progress=self.annotationJobProgress.emit,
        )
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.signals.progress.connect(self._set_annotation_job_status)
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

    def _clear_similarity_annotation_thread(self) -> None:
        self._sim_ann_thread = None
        self._sim_ann_worker = None
        self._set_similarity_annotation_enabled(True)
        self._clear_annotation_job_status()

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
        before_domains = sum(
            1 for feat in self._map.features() if feat.feature_type == "domain"
        )
        self._map.set_features(list(result.features))
        after_domains = sum(
            1 for feat in result.features if feat.feature_type == "domain"
        )
        self._clear_annotation_job_status()
        self._log_annotation_job(
            "JOB FINISH",
            "Similarity annotation",
            detail=result.message,
        )
        self.statusMessage.emit(f"{result.message} Click Save to write changes to disk.")
        added_domains = after_domains - before_domains
        lower_message = result.message.lower()
        needs_dialog = (
            added_domains <= 0
            or "not indexed" in lower_message
            or "not prepared" in lower_message
            or "no blast" in lower_message
            or "no hmm" in lower_message
            or "no cdd" in lower_message
            or "no orf" in lower_message
            or "failed" in lower_message
            or "not installed" in lower_message
            or "unchanged (0 hit" in lower_message
        )
        if (
            self._similarity_ann_settings.show_info_window_after_run
            and needs_dialog
            and result.message.strip()
        ):
            self._show_wide_information("Similarity annotation", result.message)

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
        self._set_auto_annotation_enabled(False)
        labels = ", ".join(_AUTO_MODE_LABEL[mode] for mode in modes)
        self._log_annotation_job(
            "JOB START",
            "Auto annotation",
            detail=f"modes=[{labels}]",
        )
        self._set_annotation_job_status(f"Running auto annotation: {labels}…")
        self.statusMessage.emit(f"Running auto annotation: {labels}…")

        worker = Worker(
            run_auto_annotations,
            self._path,
            tuple(modes),
            cancel=self._auto_ann_cancel,
            settings=self._auto_ann_settings,
            base_features=tuple(self._map.features()),
            progress=self.annotationJobProgress.emit,
        )
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.signals.progress.connect(self._set_annotation_job_status)
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

    def _clear_auto_annotation_thread(self) -> None:
        self._auto_ann_thread = None
        self._auto_ann_worker = None
        self._set_auto_annotation_enabled(True)
        self._clear_annotation_job_status()

    def _on_auto_annotation_finished(self, result: object) -> None:
        if isinstance(result, MultiAutoAnnotationResult):
            data = self._map.data()
            if data is None:
                self._clear_annotation_job_status()
                return
            features = list(result.features)
            self._map.set_features(features)
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
        self._clear_annotation_job_status()
        self._log_annotation_job(
            "JOB FINISH",
            "Auto annotation",
            detail=result.message,
        )
        self.statusMessage.emit(f"{result.message} Click Save to write changes to disk.")

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

    def _export_annotation_view_to_file(
        self,
        out_path: Path,
    ) -> None:
        """Export the current SequenceMapWidget view to disk."""
        widget = self._map
        header = getattr(widget, "_header", None)
        seq_block = getattr(widget, "_sequence_block", None)
        # Compose: [header (overview+ruler)] + [main sequence/annotation block]
        if header is not None and seq_block is not None:
            header_size = header.rect().size()
            if header_size.width() <= 0 or header_size.height() <= 0:
                header_size = header.sizeHint()
            seq_size = seq_block.rect().size()
            if seq_size.width() <= 0 or seq_size.height() <= 0:
                seq_size = seq_block.sizeHint()
            width = max(1, int(max(header_size.width(), seq_size.width())))
            height = max(1, int(header_size.height() + seq_size.height()))
            export_sizes = (width, height, header_size.height())
        else:
            export_sizes = (max(1, int(widget.rect().width())), max(1, int(widget.rect().height())), 0)
        width, height, header_h = export_sizes
        if width <= 0 or height <= 0:
            QMessageBox.information(self, "Export", "Nothing to export (empty view).")
            return
        size = QSize(width, height)

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
                return

            gen = QSvgGenerator()
            gen.setFileName(str(out_path))
            gen.setSize(size)
            gen.setViewBox(QRect(0, 0, size.width(), size.height()))
            gen.setResolution(96)
            painter = QPainter(gen)
            try:
                if header is not None and seq_block is not None:
                    header.render(painter, QPoint(0, 0))
                    seq_block.render(painter, QPoint(0, int(header_h)))
                else:
                    widget.render(painter, QPoint(0, 0))
            finally:
                painter.end()
            return

        # Render once (PNG-friendly) then reuse for PNG/PDF.
        img = QImage(size, QImage.Format.Format_ARGB32)
        img.fill(Qt.GlobalColor.transparent)
        painter = QPainter(img)
        try:
            if header is not None and seq_block is not None:
                header.render(painter, QPoint(0, 0))
                seq_block.render(painter, QPoint(0, int(header_h)))
            else:
                widget.render(painter, QPoint(0, 0))
        finally:
            painter.end()

        if suffix == ".png":
            if not img.save(str(out_path), "PNG"):
                QMessageBox.warning(self, "Export", "Could not save PNG.")
            return

        if suffix == ".pdf":
            writer = QPdfWriter(str(out_path))
            writer.setPageSize(QPageSize(QPageSize.PageSizeId.A4))
            writer.setResolution(96)
            pdf_painter = QPainter(writer)
            try:
                page_rect = writer.pageRect()
                scaled = img.scaled(
                    page_rect.size(),
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation,
                )
                x = page_rect.x() + (page_rect.width() - scaled.width()) // 2
                y = page_rect.y() + (page_rect.height() - scaled.height()) // 2
                pdf_painter.drawImage(QPoint(x, y), scaled)
            finally:
                pdf_painter.end()
            return

        QMessageBox.warning(
            self,
            "Export",
            f"Unsupported export format: {suffix!r}. Use .svg, .png, or .pdf.",
        )

    def _export_annotation_view(self) -> None:
        out, _filter = QFileDialog.getSaveFileName(
            self,
            "Export annotation view",
            str(self._path.with_suffix(".svg")),
            "SVG (*.svg);;PNG (*.png);;PDF (*.pdf)",
        )
        if not out:
            return
        self._export_annotation_view_to_file(Path(out))
        # Avoid letting focus fall back to the left Files panel after export.
        QTimer.singleShot(0, self._map.setFocus)
