"""Interactive single-sequence map: scroll, zoom to nucleotide, feature tracks."""

from __future__ import annotations

import math
from pathlib import Path

from PySide6.QtCore import QEvent, QPoint, QPointF, QRect, QRectF, QSize, Qt, QThread, Signal, QTimer
from PySide6.QtGui import (
    QColor,
    QFont,
    QFontMetrics,
    QKeyEvent,
    QKeySequence,
    QMouseEvent,
    QNativeGestureEvent,
    QPainter,
    QPen,
    QPolygonF,
    QShortcut,
    QWheelEvent,
)
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QScrollBar,
    QSizePolicy,
    QSlider,
    QSpinBox,
    QStyle,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from teviewer.backends.auto_annotation import (
    AutoAnnotationMode,
    AutoAnnotationResult,
    MultiAutoAnnotationResult,
    merge_auto_annotation_features,
    run_auto_annotation,
    run_auto_annotations,
)
from teviewer.backends.auto_annotation_settings import (
    AutoAnnotationSettings,
    DEFAULT_TIR_SIMILARITY_THRESHOLD,
    OrfSettings,
    TerminalRepeatSettings,
    TandemRepeatSettings,
    load_auto_annotation_settings,
    save_auto_annotation_settings,
)
from teviewer.backends.job_cancel import JobCancellation
from teviewer.backends.annotation_colors import feature_type_color, similarity_database_color
from teviewer.backends.sequence_features import (
    DEFAULT_FEATURE_TYPES,
    SequenceFeature,
    SequenceMapData,
    features_gff_path,
    format_single_sequence_summary_text,
    load_sequence_map_data,
    resolve_features_gff_path,
    save_annotations,
    summarize_single_sequence_fasta,
)
from teviewer.database_settings import load_similarity_annotation_databases
from teviewer.backends.similarity_annotation_settings import (
    SimilarityAnnotationSettings,
    load_similarity_annotation_settings,
    save_similarity_annotation_settings,
    similarity_feature_track_label,
)
from teviewer.job_options import configure_report_evalue_spinbox
from teviewer.backends.similarity_annotation import (
    SimilarityAnnotationResult,
    SimilarityDatabaseRun,
    load_use_system_blast_from_settings,
    run_similarity_annotations,
)
from teviewer.project import UserDatabase
from teviewer.workers.base import Worker
from teviewer.widgets.annotation_colors_dialog import AnnotationElementColorsDialog
from teviewer.widgets.preview_controls import ZOOM_SLIDER_HEIGHT, zoom_slider_stylesheet
from teviewer.widgets.preview_zoom import nucleotide_base_color_map

_ANNOTATION_CONTROLS_STYLE = (
    "QWidget#annotationControlsPanel {"
    " background: palette(alternate-base);"
    " border-top: 1px solid palette(mid);"
    "}"
    " QWidget#annotationControlBand {"
    " border-bottom: 1px solid palette(mid);"
    " background: transparent;"
    "}"
    " QWidget#annotationActionBand {"
    " border-top: 1px solid palette(mid);"
    " background: transparent;"
    "}"
    " QWidget#annotationControlBand QLabel {"
    " background: transparent;"
    " border: none;"
    "}"
    " QWidget#annotationControlsPanel QPushButton {"
    " padding: 2px 0px;"
    " border: none;"
    " border-radius: 4px;"
    " background: transparent;"
    " min-height: 0px;"
    "}"
    " QWidget#annotationActionBand QPushButton {"
    " padding: 2px 12px;"
    " }"
    " QWidget#annotationRunSettingsBox QPushButton {"
    " padding: 2px 4px;"
    " }"
    " QWidget#annotationControlsPanel QPushButton:hover {"
    " background: palette(window);"
    "}"
    " QWidget#annotationControlsPanel QPushButton:pressed {"
    " background: palette(midlight);"
    "}"
    " QWidget#annotationControlsPanel QCheckBox {"
    " spacing: 4px;"
    " background: transparent;"
    " border: none;"
    " min-height: 22px;"
    "}"
    " QWidget#annotationControlsPanel QCheckBox::indicator {"
    " width: 14px;"
    " height: 14px;"
    "}"
    " QWidget#annotationFeatureColumnCell {"
    " background: transparent;"
    " min-height: 22px;"
    "}"
    " QLabel#annotationJobStatus {"
    " color: palette(placeholder-text);"
    " padding: 0 4px;"
    " min-height: 22px;"
    "}"
)
_ADD_FEATURE_BUTTONS: tuple[tuple[str, str], ...] = (
    ("LTR", "LTR"),
    ("TIR", "TIR"),
    ("ORF", "ORF"),
    ("Tandem", "tandem repeat"),
    ("PBS", "PBS"),
    ("PPT", "PPT"),
    ("Others", "others"),
)
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
_AUTO_ANN_CHECK_EXTRA_SPACING = 6
_MIN_FEATURE_COLUMN_WIDTH = 54
_SIMILARITY_DB_CHIP_SPACING = 4
_MARQUEE_DRAG_THRESHOLD = 3.0
_ANNOT_TRACK_VPAD = 8
_ANNOT_FEATURE_HEIGHT = 18
_ANNOT_MIN_VIEWPORT_HEIGHT = 56
_OVERVIEW_SYNC_MS = 16
_RULER_HEIGHT = 16
_OVERVIEW_HEIGHT = 16
_SEQUENCE_ROW = 14
_TRACK_HEIGHT = 12
_TRACK_GAP = 0
_FEATURE_UNDO_LIMIT = 32
_LEFT_MARGIN = 8
_RIGHT_MARGIN = 8
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
_SELECTION_FILL = QColor(25, 95, 175, 175)
_SELECTION_BORDER = QColor(15, 60, 120, 230)
_FEATURE_SELECTION_BORDER = QColor("#0b5cab")
_FEATURE_SELECTION_OUTLINE = QColor("#ffffff")
_FEATURE_SELECTION_OVERLAY = QColor(11, 92, 171, 120)
_MIN_NUC_FONT_PT = 6
_MAX_NUC_FONT_PT = 36
_DEFAULT_NUC_FONT_PT = 13
_ZOOM_BTN_SIZE = 28
_ZOOM_BTN_HEIGHT = 24


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


def _annotation_run_settings_box(
    run_btn: QPushButton,
    settings_btn: QPushButton,
    parent: QWidget,
) -> QWidget:
    """Fixed-width Run/Setting column shared by auto and similarity rows."""
    box = QWidget(parent)
    box.setObjectName("annotationRunSettingsBox")
    box.setFixedWidth(_ANNOTATION_RUN_SETTINGS_WIDTH)
    layout = QHBoxLayout(box)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(14)
    layout.addWidget(run_btn, 0, Qt.AlignmentFlag.AlignVCenter)
    layout.addWidget(settings_btn, 0, Qt.AlignmentFlag.AlignVCenter)
    layout.addStretch()
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
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
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


def _checkbox_label_offset(widget: QWidget) -> int:
    style = widget.style()
    return (
        style.pixelMetric(QStyle.PixelMetric.PM_IndicatorWidth, None, widget)
        + style.pixelMetric(QStyle.PixelMetric.PM_CheckBoxLabelSpacing, None, widget)
    )


class _AddFeatureControl(QWidget):
    """Add-feature button with leading offset so label text lines up with checkboxes."""

    clicked = Signal()

    def __init__(self, label: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        self._offset = QWidget(self)
        self._offset.setFixedWidth(0)
        layout.addWidget(self._offset)
        self._btn = QPushButton(label, self)
        self._btn.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        self._btn.clicked.connect(self.clicked.emit)
        layout.addWidget(self._btn)

    def set_label_offset(self, width: int) -> None:
        self._offset.setFixedWidth(max(0, int(width)))

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


def _letter_cell_width_for_pt(point_size: int) -> int:
    font = QFont("Menlo, Monaco, Courier New, monospace")
    font.setPointSize(point_size)
    font.setStyleHint(QFont.StyleHint.Monospace)
    font.setFixedPitch(True)
    fm = QFontMetrics(font)
    return max(1, max(fm.horizontalAdvance(base) for base in "ACGTUN"))


def _base_colors(base: str) -> tuple[QColor, QColor]:
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
        self.setFixedHeight(_RULER_HEIGHT)
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
        ruler_font = QFont("Menlo, Monaco, Courier New, monospace")
        ruler_font.setPointSize(7)
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
                painter.drawLine(int(x), 0, int(x), _RULER_HEIGHT - 3)
            text = str(bp)
            if x - last_label_x < min_label_w and bp != start_bp:
                continue
            painter.drawText(int(x + 2), _RULER_HEIGHT - 3, text)
            last_label_x = x
        painter.end()


def _is_delete_key(key: int) -> bool:
    return key in (int(Qt.Key.Key_Delete), int(Qt.Key.Key_Backspace))


class AnnotationTrackWidget(QWidget):
    """Feature annotations below the sequence; stacked rows for overlaps."""

    selectionChanged = Signal(object)  # frozenset[int]
    deleteRequested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMouseTracking(True)
        self.setStyleSheet("background: #f5f5f5; border: none; margin: 0px; padding: 0px;")
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
        self._release_drag_mouse()
        super().hideEvent(event)

    def _row_height(self) -> float:
        return float(_ANNOT_FEATURE_HEIGHT + _TRACK_GAP)

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
        if not self._marquee_active and self._press_pos is None:
            self._selected_indices = set(selected_indices or [])
        valid = set(range(len(self._features)))
        self._selected_indices &= valid
        self._row_count = max((r for r in feature_rows), default=-1) + 1 if features else 0
        track_w = max(1, int(viewport_w))
        if self._row_count > 0:
            h = int(self._row_count * self._row_height()) + 2 * _ANNOT_TRACK_VPAD
            h = max(h, _ANNOT_FEATURE_HEIGHT + 2 * _ANNOT_TRACK_VPAD)
            self.setFixedSize(track_w, h)
            self.show()
        else:
            self.setFixedSize(0, 0)
            self.hide()
        self.update()

    def has_content(self) -> bool:
        return self._row_count > 0

    def paintEvent(self, event) -> None:  # noqa: N802
        if not self._features:
            return
        painter = QPainter(self)
        painter.fillRect(self.rect(), QColor("#f5f5f5"))
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
        self._press_pos = pos
        self._press_hit = self._hit_test_feature(pos)
        self._marquee_active = False
        self._marquee_origin = None
        self._marquee_current = None
        self.setFocus(Qt.FocusReason.MouseFocusReason)
        event.accept()

    def mouseMoveEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if self._press_pos is None:
            super().mouseMoveEvent(event)
            return
        pos = event.position()
        if self._marquee_active:
            self._marquee_current = pos
            self._apply_marquee_preview(additive=self._additive_modifier(event.modifiers()))
            self.update()
            return
        delta = pos - self._press_pos
        if math.hypot(delta.x(), delta.y()) < _MARQUEE_DRAG_THRESHOLD:
            super().mouseMoveEvent(event)
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
        if self._press_pos is None and not self._marquee_active:
            self._release_drag_mouse()
            super().mouseReleaseEvent(event)
            return
        pos = event.position()
        additive = self._additive_modifier(event.modifiers())
        try:
            if self._marquee_active:
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


class SequenceMapCanvas(QWidget):
    """Paints sequence + feature tracks; supports wrap and stacked feature rows."""

    selectionChanged = Signal(int, int)  # 1-based inclusive
    featureClicked = Signal(int)  # index in features list

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
        self._font_point_size = _DEFAULT_NUC_FONT_PT
        self._selected_indices: set[int] = set()
        self._select_start_bp: int | None = None
        self._select_end_bp: int | None = None
        self._drag_select = False
        self._mouse_grabbed = False
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._show_similarity_database_name = True

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

    def hideEvent(self, event) -> None:  # noqa: N802
        self._drag_select = False
        self._release_drag_mouse()
        super().hideEvent(event)

    def set_font_point_size(self, point_size: int) -> None:
        self._font_point_size = max(_MIN_NUC_FONT_PT, min(_MAX_NUC_FONT_PT, int(point_size)))
        self._update_geometry()
        self.update()

    def font_point_size(self) -> int:
        return self._font_point_size

    def _sequence_top(self) -> int:
        return 0

    def set_sequence(self, sequence: str, features: list[SequenceFeature]) -> None:
        self._sequence = sequence.upper()
        self._features = list(features)
        self._feature_rows, self._feature_row_count = _pack_feature_rows(self._features)
        self._update_geometry()
        self.update()

    def set_pixels_per_bp(self, ppb: float) -> None:
        cap = _zoom_ppb_max_for(len(self._sequence))
        self._pixels_per_bp = max(_FIT_MIN_PPB, min(cap, float(ppb)))
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
        self.update()

    def sync_viewport_from_widget(self) -> None:
        """Read the layout-assigned width; never expand the parent."""
        width = max(120, self.width())
        if width != self._viewport_width:
            self._viewport_width = width
            self.update()

    def set_scroll_px(self, scroll_px: float) -> None:
        self._scroll_px = max(0.0, float(scroll_px))
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
        self.update()

    def color_bases(self) -> bool:
        return self._color_bases

    def set_selected_features(self, indices: set[int] | frozenset[int] | None) -> None:
        self._selected_indices = set(indices or [])
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
        font = QFont("Menlo, Monaco, Courier New, monospace")
        font.setPointSize(self._font_point_size)
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setFixedPitch(True)
        return font

    def _font_letter_cell_size(self) -> tuple[int, int]:
        """Minimum monospace cell size for legible nucleotide letters."""
        fm = QFontMetrics(self._sequence_font())
        cell_w = max(1, max(fm.horizontalAdvance(base) for base in "ACGTUN"))
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
        return max(_SEQUENCE_ROW + 2, int(self._seq_row_height()) + 2)

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
            return _base_colors(base)
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
            painter.setPen(fg)
            painter.drawText(cell, Qt.AlignmentFlag.AlignCenter, base)

    def _draw_ruler(self, painter: QPainter, seq_len: int, ppb: float) -> None:
        painter.setPen(QColor("#666666"))
        ruler_font = QFont("Menlo, Monaco, Courier New, monospace")
        ruler_font.setPointSize(7)
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
                painter.drawLine(int(x), 0, int(x), _RULER_HEIGHT - 3)
            else:
                painter.drawLine(int(x), _RULER_HEIGHT - 6, int(x), _RULER_HEIGHT - 3)

        label_positions = list(range(1, seq_len + 1, major_step))
        if label_positions[-1] != seq_len:
            label_positions.append(seq_len)
        for bp in label_positions:
            x = _LEFT_MARGIN + (bp - 1) * ppb
            if x - last_label_x < min_label_w and bp != 1:
                continue
            painter.drawText(int(x) + 2, _RULER_HEIGHT - 3, str(bp))
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

        if cell_w >= _BUCKET_DRAW_PPB:
            visible_start = max(0, int(math.floor((scroll + clip.left() - origin_x) / cell_w)))
            visible_end = min(
                seq_len,
                int(math.ceil((scroll + clip.right() - origin_x) / cell_w)),
            )
            seq_left, seq_right = self._sequence_viewport_x_extent()
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

        seq_left, seq_right = self._sequence_viewport_x_extent()
        px_left = max(0, int(math.floor(clip.left())), int(math.floor(seq_left)))
        px_right = min(
            self.width(),
            int(math.ceil(clip.right())) + 1,
            int(math.ceil(seq_right)),
        )
        if px_left >= px_right:
            return
        painter.setPen(Qt.PenStyle.NoPen)
        for px in range(px_left, px_right):
            content_x0 = scroll + px - origin_x
            if content_x0 < 0.0 or content_x0 >= seq_len * cell_w:
                continue
            bp0 = min(seq_len - 1, int(content_x0 / cell_w))
            bg, _ = self._base_style(self._sequence[bp0])
            painter.fillRect(QRectF(float(px), float(origin_y), 1.0, float(cell_h)), bg)

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
        painter.fillRect(self.rect(), QColor("#fafafa"))
        if not self._sequence:
            painter.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No sequence")
            painter.end()
            return

        seq_len = len(self._sequence)
        ppb = self.effective_ppb()
        clip = QRectF(event.rect())

        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, False)

        self._draw_sequence(painter, clip)

        sel = self.selection_range()
        if sel:
            a, b = sel
            origin_x, origin_y, cell_w, cell_h = self._nucleotide_cell_geometry()
            if self.uses_wrap():
                cell_w_f, cell_h_f, _ = self._cell_size()
                for bp in range(a, b + 1):
                    x, y = self._index_to_xy(bp - 1)
                    rect = QRectF(x, y, cell_w_f, cell_h_f)
                    painter.fillRect(rect, _SELECTION_FILL)
                    painter.setPen(QPen(_SELECTION_BORDER, 1))
                    painter.drawRect(rect)
            else:
                x1 = origin_x + (a - 1) * cell_w - self._scroll_px
                x2 = origin_x + b * cell_w - self._scroll_px
                rect = QRectF(x1, float(origin_y), x2 - x1, float(cell_h))
                painter.fillRect(rect, _SELECTION_FILL)
                painter.setPen(QPen(_SELECTION_BORDER, 1))
                painter.drawRect(rect)
            painter.setPen(Qt.PenStyle.NoPen)

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
        self._select_start_bp = bp
        self._select_end_bp = bp
        self._drag_select = True
        self.grabMouse()
        self._mouse_grabbed = True
        self.selectionChanged.emit(bp, bp)
        self.update()
        event.accept()

    def mouseMoveEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if self._drag_select and self._select_start_bp is not None:
            bp = self._bp_at_pos(event.position().x(), event.position().y())
            self._select_end_bp = bp
            a, b = sorted((self._select_start_bp, bp))
            self.selectionChanged.emit(a, b)
            self.update()
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._drag_select = False
            self._release_drag_mouse()
            event.accept()
            return
        super().mouseReleaseEvent(event)

class SequenceMapWidget(QWidget):
    """Scrollable sequence map with overview and zoom controls."""

    featuresChanged = Signal()
    statusMessage = Signal(str)
    deleteFeaturesRequested = Signal()

    def __init__(self, parent: QWidget | None = None, *, embedded: bool = False) -> None:
        super().__init__(parent)
        self._embedded = embedded
        self._data: SequenceMapData | None = None
        self._dirty = False
        self._selected_feature_indices: set[int] = set()
        self._fit_retries = 0
        self._show_similarity_database_name = True
        if embedded:
            self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        else:
            self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

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
        self._zoom_slider.setMinimumWidth(72)
        self._zoom_slider.setFixedHeight(ZOOM_SLIDER_HEIGHT)
        self._zoom_slider.setToolTip(
            "Left: fit full sequence in window. Right: maximum zoom-in."
        )
        self._zoom_slider.setStyleSheet(zoom_slider_stylesheet("annotationZoomSlider"))
        self._zoom_slider.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._zoom_slider.valueChanged.connect(self._on_zoom_changed)

        zoom_strip_layout.addWidget(zoom_out_btn)
        zoom_strip_layout.addWidget(zoom_in_btn)
        zoom_strip_layout.addWidget(self._zoom_slider, 1)

        self._sequence_stats_label = QLabel("", self)
        self._sequence_stats_label.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        self._sequence_stats_label.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse
        )
        self._sequence_stats_label.setStyleSheet(
            "QLabel { color: palette(text); padding: 0 2px; }"
        )

        self._font_size = QSpinBox()
        self._font_size.setRange(_MIN_NUC_FONT_PT, _MAX_NUC_FONT_PT)
        self._font_size.setValue(_DEFAULT_NUC_FONT_PT)
        self._font_size.setSuffix(" pt")
        self._font_size.setFixedWidth(58)
        self._font_size.setToolTip("Nucleotide letter size (independent of zoom scale)")
        self._font_size.valueChanged.connect(self._on_font_size_changed)

        self._pos_start = QSpinBox()
        self._pos_start.setMinimum(1)
        self._pos_start.setMaximum(1)
        self._pos_start.setFixedWidth(64)
        self._pos_end = QSpinBox()
        self._pos_end.setMinimum(1)
        self._pos_end.setMaximum(1)
        self._pos_end.setFixedWidth(64)
        self._pos_start.valueChanged.connect(self._on_position_spin_changed)
        self._pos_end.valueChanged.connect(self._on_position_spin_changed)
        select_btn = QPushButton("Select")
        select_btn.setToolTip(
            "Zoom to the selected range and center its midpoint in the window"
        )
        select_btn.clicked.connect(self._select_by_position)

        zoom_row.addWidget(zoom_strip, 1, Qt.AlignmentFlag.AlignVCenter)
        zoom_row.addWidget(QLabel("Position:"))
        zoom_row.addWidget(self._pos_start)
        zoom_row.addWidget(QLabel("to"))
        zoom_row.addWidget(self._pos_end)
        zoom_row.addWidget(select_btn)
        zoom_row.addWidget(QLabel("Font:"))
        zoom_row.addWidget(self._font_size)
        zoom_row.addWidget(self._sequence_stats_label)

        self._canvas = SequenceMapCanvas(self)
        self._canvas.set_color_bases(True)
        self._canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._canvas.selectionChanged.connect(self._on_selection_changed)
        self._canvas.featureClicked.connect(self._on_canvas_feature_clicked)

        self._ruler = ViewportRulerWidget(self)
        self._annot_track = AnnotationTrackWidget(self)
        self._annot_track.selectionChanged.connect(self._on_annot_track_selection_changed)
        self._annot_track.deleteRequested.connect(self._request_delete_features)
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

        block_layout = QVBoxLayout(self._sequence_block)
        block_layout.setContentsMargins(0, 0, 0, 0)
        block_layout.setSpacing(0)
        block_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        block_layout.addWidget(self._canvas_row, 0, Qt.AlignmentFlag.AlignTop)
        block_layout.addWidget(self._annot_scroll, 1)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addLayout(zoom_row)
        layout.addWidget(self._header)
        layout.addWidget(self._sequence_block, stretch=1)

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

    def refresh_layout(self) -> None:
        self._sync_header_views()

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
        return 120

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
        if not self._annot_track.has_content() or self._annot_track.height() <= 0:
            self._annot_scroll.hide()
            self._annot_scroll.setMinimumHeight(0)
            self._annot_scroll.setMaximumHeight(0)
            return

        self._annot_scroll.show()
        vw = max(1, self._scroll_viewport_width())
        track_h = self._annot_track.height()
        self._annot_track.setFixedSize(vw, track_h)
        self._annot_scroll.setMinimumHeight(_ANNOT_MIN_VIEWPORT_HEIGHT)
        self._annot_scroll.setMaximumHeight(16777215)

    def _set_position_spin_limits(self, seq_len: int) -> None:
        seq_len = max(1, int(seq_len))
        for spin in (self._pos_start, self._pos_end):
            spin.setMaximum(seq_len)
            spin.setValue(min(spin.value(), seq_len) if spin.value() > 0 else 1)

    def _sync_position_spins(self, start: int, end: int) -> None:
        a, b = sorted((int(start), int(end)))
        self._pos_start.blockSignals(True)
        self._pos_end.blockSignals(True)
        self._pos_start.setValue(a)
        self._pos_end.setValue(b)
        self._pos_start.blockSignals(False)
        self._pos_end.blockSignals(False)

    def _on_position_spin_changed(self, _value: int) -> None:
        if self._data is None:
            return
        self._canvas.set_selection_range(self._pos_start.value(), self._pos_end.value())

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

    def select_by_position(self, start: int, end: int) -> None:
        if self._data is None:
            return
        a, b = sorted((int(start), int(end)))
        self._canvas.set_selection_range(a, b)
        self._sync_position_spins(a, b)

        span = max(1, b - a + 1)
        viewport_w = self._scroll_viewport_width()
        usable_w = max(1.0, float(viewport_w - _LEFT_MARGIN - _RIGHT_MARGIN))
        if usable_w > 0:
            target_ppb = max(_FIT_MIN_PPB, min(self._canvas.zoom_ppb_max(), usable_w / span))
            self._apply_zoom_ppb(target_ppb, sync_slider=True)

        self._scroll_to_bp_range(a, b)
        QTimer.singleShot(0, lambda lo=a, hi=b: self._scroll_to_bp_range(lo, hi))

    def _select_by_position(self) -> None:
        if self._data is None:
            return
        sel = self._canvas.selection_range()
        if sel is not None:
            start, end = sel
        else:
            start, end = self._pos_start.value(), self._pos_end.value()
        self.select_by_position(start, end)

    def _scroll_to_bp_range(self, start: int, end: int) -> None:
        viewport_w = self._scroll_viewport_width()
        self._scroll_preserve_viewport_x(start, end, viewport_w / 2.0)

    def is_dirty(self) -> bool:
        return self._dirty

    def data(self) -> SequenceMapData | None:
        return self._data

    def load(self, data: SequenceMapData) -> None:
        self._data = data
        self._dirty = False
        self._set_position_spin_limits(data.length)
        if self._font_size.value() != _DEFAULT_NUC_FONT_PT:
            self._font_size.blockSignals(True)
            self._font_size.setValue(_DEFAULT_NUC_FONT_PT)
            self._font_size.blockSignals(False)
            self._canvas.set_font_point_size(_DEFAULT_NUC_FONT_PT)
        self._canvas.set_sequence(data.sequence, data.features)
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
        QTimer.singleShot(0, self.fit_window)

    def features(self) -> list[SequenceFeature]:
        if self._data is None:
            return []
        return list(self._data.features)

    def set_features(self, features: list[SequenceFeature]) -> None:
        if self._data is None:
            return
        self._data.features = list(features)
        self._canvas.set_sequence(self._data.sequence, self._data.features)
        valid = set(range(len(self._data.features)))
        self._selected_feature_indices &= valid
        self._dirty = True
        self.featuresChanged.emit()
        self._sync_overview()

    def mark_clean(self) -> None:
        self._dirty = False

    def zoom_in(self) -> None:
        self._zoom_slider.setValue(min(100, self._zoom_slider.value() + 8))

    def zoom_out(self) -> None:
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
            self._zoom_slider.blockSignals(True)
            self._zoom_slider.setValue(self._ppb_to_slider_value(ppb))
            self._zoom_slider.blockSignals(False)

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
        if viewport_w <= 20:
            self._fit_retries += 1
            if self._fit_retries < 30:
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

    def set_sequence_summary_text(self, text: str) -> None:
        self._sequence_stats_label.setText(text.strip())
        self._sequence_stats_label.setVisible(bool(text.strip()))

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
        else:
            a, b = sorted((start, end))
            self.statusMessage.emit(f"Selection: {a} – {b} ({b - a + 1} bp)")

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

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802
        if _is_delete_key(event.key()):
            self._sync_feature_selection_from_annot_track()
            if self._selected_feature_indices:
                self.deleteFeaturesRequested.emit()
                event.accept()
                return
        super().keyPressEvent(event)

    def _scroll_to_fraction(self, fraction: float) -> None:
        if self._data is None:
            return
        max_scroll = self._hbar.maximum()
        target = int(round(fraction * max_scroll))
        self._hbar.setValue(max(0, min(max_scroll, target)))
        self._sync_overview_now()


class AddFeatureDialog(QDialog):
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
        a, b = sorted((start, end))
        self._name = QLineEdit(self)
        self._name.setPlaceholderText("e.g. LTR5")
        if default_name:
            self._name.setText(default_name)
        self._type = QComboBox(self)
        if feature_type is not None:
            self._type.addItem(feature_label or feature_type, feature_type)
            self._type.setEnabled(False)
        else:
            for t in DEFAULT_FEATURE_TYPES:
                self._type.addItem(t)
        self._start = QSpinBox(self)
        self._start.setRange(1, seq_length)
        self._start.setValue(a)
        self._end = QSpinBox(self)
        self._end.setRange(1, seq_length)
        self._end.setValue(b)
        self._strand = QComboBox(self)
        self._strand.addItem("Forward (+)", "+")
        self._strand.addItem("Reverse (−)", "-")

        form = QFormLayout()
        form.addRow("Name:", self._name)
        form.addRow("Type:", self._type)
        form.addRow("Direction:", self._strand)
        form.addRow("Start (bp):", self._start)
        form.addRow("End (bp):", self._end)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def feature(self) -> SequenceFeature | None:
        name = self._name.text().strip()
        if not name:
            return None
        start = self._start.value()
        end = self._end.value()
        ftype = str(self._type.currentData() or self._type.currentText())
        strand = str(self._strand.currentData())
        return SequenceFeature(start=start, end=end, name=name, feature_type=ftype, strand=strand)


class AutoAnnotationSettingsDialog(QDialog):
    """Configure LTR, TIR, tandem-repeat, and ORF auto-annotation parameters."""

    def __init__(
        self,
        settings: AutoAnnotationSettings,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Auto annotation settings")
        self._settings = settings

        layout = QVBoxLayout(self)
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

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

        self._ltr_fields = self._terminal_fields[0]
        self._tir_fields = self._terminal_fields[1]
        self._tandem_fields = self._tandem_field_list
        self._orf_fields = self._orf_field_list

    def _section_title(self, text: str) -> QLabel:
        title = QLabel(text)
        title.setStyleSheet("font-weight: 600;")
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


class SimilarityAnnotationSettingsDialog(QDialog):
    """Configure similarity-annotation search thresholds."""

    def __init__(
        self,
        settings: SimilarityAnnotationSettings,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Similarity annotation settings")
        layout = QVBoxLayout(self)
        layout.setSpacing(8)

        intro = QLabel(
            "These E-value and output limits apply only to Similarity annotation on this tab. "
            "They do not change TEAid, genome BLAST, or other tools in TEviewer.",
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
            "TEviewer first removes shorter fully nested overlaps when their E-value "
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
        layout.addLayout(form)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)

    def settings(self) -> SimilarityAnnotationSettings:
        return SimilarityAnnotationSettings(
            blast_evalue=self._blast_evalue.value(),
            hmm_evalue=self._hmm_evalue.value(),
            cdd_evalue=self._cdd_evalue.value(),
            max_hits_per_database=self._max_hits.value(),
            show_database_name_on_track=self._show_db_name.isChecked(),
        )


class SequenceMapPanel(QWidget):
    """Embedded sequence map with annotation toolbar (preview Annotation tab)."""

    statusMessage = Signal(str)
    annotationJobProgress = Signal(str)
    annotationsSaved = Signal(object)  # Path to the written GFF3

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
        self._map = SequenceMapWidget(self, embedded=True)
        self._map.statusMessage.connect(self.statusMessage.emit)
        self._map.deleteFeaturesRequested.connect(
            lambda: self._delete_feature(confirm_empty=False)
        )
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
        self._app_event_filter_installed = False

        controls = QWidget(self)
        controls.setObjectName("annotationControlsPanel")
        self._controls_panel = controls
        controls.setMinimumWidth(0)
        controls.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Minimum,
        )
        controls.setStyleSheet(_ANNOTATION_CONTROLS_STYLE)
        controls_layout = QVBoxLayout(controls)
        controls_layout.setContentsMargins(0, 0, 0, 0)
        controls_layout.setSpacing(0)

        self._auto_all_updating = False
        self._similarity_all_updating = False
        self._add_feature_buttons: list[_AddFeatureControl] = []

        self._add_row_container, add_row = _control_band(controls)
        self._add_feature_label = _annotation_row_label("Add feature:", self._add_row_container)
        add_row.addWidget(self._add_feature_label, 0, Qt.AlignmentFlag.AlignVCenter)
        self._add_feature_columns = _FeatureColumnRow(self._add_row_container)
        for label, feature_type in _ADD_FEATURE_BUTTONS:
            control = _AddFeatureControl(label, self._add_feature_columns)
            control.setToolTip(f"Add a {label} feature for the selected region")
            control.clicked.connect(
                lambda _checked=False, ft=feature_type, bl=label: self._add_feature_of_type(
                    ft, feature_label=bl
                )
            )
            self._add_feature_buttons.append(control)
            self._add_feature_columns.place(_ANNOTATION_FEATURE_COLUMN.get(label, 0), control)
        add_row.addWidget(self._add_feature_columns, 0, Qt.AlignmentFlag.AlignVCenter)
        self._add_feature_overflow = QComboBox(self._add_row_container)
        for label, feature_type in _ADD_FEATURE_BUTTONS:
            self._add_feature_overflow.addItem(label, feature_type)
        self._add_feature_overflow.setMinimumWidth(110)
        self._add_feature_overflow.activated.connect(self._on_add_feature_overflow_activated)
        self._add_feature_overflow.hide()
        add_row.addWidget(self._add_feature_overflow, 0, Qt.AlignmentFlag.AlignVCenter)
        add_row.addStretch(1)
        controls_layout.addWidget(self._add_row_container)

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

        self._similarity_feature_columns = _FeatureColumnRow(self._sim_row_container)
        self._similarity_db_flow = _SimilarityDbFlowWidget(self._similarity_feature_columns)
        self._similarity_feature_columns.place(0, self._similarity_db_flow)
        sim_row.addWidget(self._similarity_feature_columns, 0, Qt.AlignmentFlag.AlignVCenter)

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
        self._action_row_label = _annotation_row_label("", self._action_row_container)
        action_row.addWidget(self._action_row_label, 0, Qt.AlignmentFlag.AlignVCenter)
        self._fit_control = _AddFeatureControl("Fit window", self._action_row_container)
        self._fit_control.setToolTip(
            "Scale the sequence so the full length fits in the visible width"
        )
        self._fit_control.clicked.connect(self._map.fit_window)
        action_row.addWidget(self._fit_control, 0, Qt.AlignmentFlag.AlignVCenter)
        del_btn = QPushButton("Delete feature", self._action_row_container)
        del_btn.setToolTip("Delete selected feature(s). Select on the track below, then press Delete or Backspace.")
        del_btn.clicked.connect(lambda: self._delete_feature(confirm_empty=True))
        restore_btn = QPushButton("Restore feature", self._action_row_container)
        restore_btn.setToolTip("Restore features deleted in this session")
        restore_btn.clicked.connect(self._restore_features)
        color_btn = QPushButton("Color", self._action_row_container)
        color_btn.setToolTip(
            "Set colors for feature types, similarity databases, and nucleotides (A/C/G/T)"
        )
        color_btn.clicked.connect(self._open_element_colors_dialog)
        save_btn = QPushButton("Save", self._action_row_container)
        save_btn.setToolTip("Write annotation changes to the JSON sidecar and GFF3 file")
        save_btn.clicked.connect(lambda: self._save_annotations(silent=True))
        for btn in (del_btn, restore_btn, color_btn, save_btn):
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
        layout.addWidget(controls)
        controls.installEventFilter(self)

        if autoload:
            self.load_path(self._path)

        self._map.set_show_similarity_database_name(
            self._similarity_ann_settings.show_database_name_on_track
        )
        self._install_annotation_shortcuts()
        QTimer.singleShot(0, self._sync_annotation_control_rows)

    def showEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().showEvent(event)
        QTimer.singleShot(0, self._sync_annotation_control_rows)

    def eventFilter(self, watched, event) -> bool:  # noqa: N802 (Qt API)
        if watched is getattr(self, "_controls_panel", None) and event.type() in (
            QEvent.Type.FocusIn,
            QEvent.Type.FocusOut,
            QEvent.Type.Show,
            QEvent.Type.Resize,
        ):
            QTimer.singleShot(0, self._sync_annotation_control_rows)
        return super().eventFilter(watched, event)

    def _install_annotation_shortcuts(self) -> None:
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

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        QTimer.singleShot(0, self._sync_annotation_control_rows)
        QTimer.singleShot(0, self._refresh_similarity_database_rows)
        app = QApplication.instance()
        if app is not None and not self._app_event_filter_installed:
            app.installEventFilter(self)
            self._app_event_filter_installed = True

    def hideEvent(self, event) -> None:  # noqa: N802
        if self._app_event_filter_installed:
            app = QApplication.instance()
            if app is not None:
                app.removeEventFilter(self)
            self._app_event_filter_installed = False
        super().hideEvent(event)

    def eventFilter(self, watched, event) -> bool:  # noqa: N802
        if not self.isVisible():
            return False
        if event.type() != QEvent.Type.KeyPress or not isinstance(event, QKeyEvent):
            return False
        if not _is_delete_key(event.key()):
            return False
        focus = QApplication.focusWidget()
        if focus is None or not self.isAncestorOf(focus):
            return False
        if isinstance(focus, (QLineEdit, QTextEdit, QPlainTextEdit)):
            return False
        self._map._sync_feature_selection_from_annot_track()
        if not self._map.selected_feature_indices():
            return False
        self._delete_feature(confirm_empty=False)
        return True

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
        width = max(
            _ANNOTATION_RUN_SETTINGS_WIDTH,
            self._auto_run_settings_box.sizeHint().width(),
            self._similarity_run_settings_box.sizeHint().width(),
        )
        self._auto_run_settings_box.setFixedWidth(width)
        self._similarity_run_settings_box.setFixedWidth(width)

    def _sync_similarity_row_layout(self) -> None:
        if not hasattr(self, "_similarity_db_flow"):
            return
        self._similarity_db_flow._reflow()
        self._sim_row_container.setFixedHeight(_CONTROL_BAND_HEIGHT)

    def _sync_add_feature_row(self) -> None:
        if not self._add_feature_buttons:
            return
        available = max(
            self._add_row_container.width(),
            self._auto_row_container.width(),
            max(0, self.width() - 24),
        )
        label_width = _ANNOTATION_ROW_LABEL_WIDTH + 12
        auto_trailing = self._auto_run_settings_box.width() + 18
        all_width = self._auto_ann_all_check.sizeHint().width() + 6
        grid_budget = max(
            120,
            available - label_width - auto_trailing - all_width,
        )
        widths = [0] * _ANNOTATION_FEATURE_COLUMN_COUNT
        label_offset = 0
        if self._auto_ann_checks:
            label_offset = _checkbox_label_offset(next(iter(self._auto_ann_checks.values())))
        for control, (label, _feature_type) in zip(self._add_feature_buttons, _ADD_FEATURE_BUTTONS):
            control.set_label_offset(label_offset)
            column = _ANNOTATION_FEATURE_COLUMN.get(label, 0)
            widths[column] = max(widths[column], _column_content_width(control))
        self._fit_control.set_label_offset(label_offset)
        for mode, check in self._auto_ann_checks.items():
            column = _AUTO_CHECK_COLUMN[mode]
            widths[column] = max(widths[column], _column_content_width(check))
        if hasattr(self, "_similarity_db_flow"):
            widths[0] = max(
                widths[0],
                _column_content_width(self._similarity_db_flow),
            )
        spacing = _AUTO_ANN_CHECK_EXTRA_SPACING
        add_total = sum(widths) + spacing * max(0, _ANNOTATION_FEATURE_COLUMN_COUNT - 1)
        if add_total <= grid_budget:
            for control in self._add_feature_buttons:
                control.show()
            for check in self._auto_ann_checks.values():
                check.show()
            self._add_feature_columns.show()
            self._auto_feature_columns.show()
            self._add_feature_overflow.hide()
            self._add_feature_columns.sync_column_widths(widths)
            self._auto_feature_columns.sync_column_widths(widths)
            self._similarity_feature_columns.sync_column_widths(widths)
        else:
            for control in self._add_feature_buttons:
                control.hide()
            for check in self._auto_ann_checks.values():
                check.show()
            self._add_feature_columns.show()
            self._auto_feature_columns.show()
            self._add_feature_overflow.show()
            self._add_feature_columns.sync_column_widths(widths)
            self._auto_feature_columns.sync_column_widths(widths)
            self._similarity_feature_columns.sync_column_widths(widths)
        self._add_feature_columns.updateGeometry()
        self._auto_feature_columns.updateGeometry()
        self._similarity_feature_columns.updateGeometry()
        self._add_row_container.updateGeometry()
        self._auto_row_container.updateGeometry()

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
                )
            )
        if changed:
            self._map.set_features(updated)

    def refresh_nucleotide_colors(self) -> None:
        self._map.refresh_nucleotide_colors()

    def open_annotation_tab(self) -> None:
        """Load/reload annotations and fit the sequence to the visible window width."""
        self.load_path(self._path)
        QTimer.singleShot(0, self._map.fit_window)
        QTimer.singleShot(50, self._map.refresh_layout)

    def fasta_path(self) -> Path:
        return self._path

    def is_dirty(self) -> bool:
        return self._map.is_dirty()

    def load_path(self, path: Path) -> None:
        self._path = path.expanduser().resolve()
        try:
            data = load_sequence_map_data(self._path)
        except (OSError, ValueError) as exc:
            self.statusMessage.emit(str(exc))
            return
        self._map.load(data)
        self._feature_undo_stack.clear()
        self.refresh_sequence_summary()

    def refresh_sequence_summary(self) -> None:
        try:
            summary = summarize_single_sequence_fasta(self._path)
        except (OSError, ValueError):
            self._map.set_sequence_summary_text("")
            return
        self._map.set_sequence_summary_text(format_single_sequence_summary_text(summary))

    def reload_annotations(self) -> None:
        """Reload features from the companion GFF3 on disk (Annotation tab open)."""
        self.open_annotation_tab()

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
        try:
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
                    "Delete feature",
                    "Select one or more features in the track below the sequence, "
                    "then press Delete or Backspace, or click Delete feature.\n\n"
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

    def _restore_features(self) -> None:
        if not self._feature_undo_stack:
            QMessageBox.information(
                self,
                "Restore feature",
                "No deleted features to restore.",
            )
            return
        features = self._feature_undo_stack.pop()
        self._map.set_features(features)
        self._map.clear_feature_selection()
        self._persist_annotations()
        self.statusMessage.emit(f"Restored {len(features)} feature(s).")

    def _set_annotation_job_status(self, text: str) -> None:
        self._annotation_job_status.setText(text.strip())

    def _clear_annotation_job_status(self) -> None:
        self._annotation_job_status.clear()

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

    def _stop_auto_annotation(self) -> None:
        if not self._auto_annotate_running():
            return
        self._auto_ann_cancel.request()
        self._set_annotation_job_status("Stopping auto annotation…")

    def _stop_similarity_annotation(self) -> None:
        if not self._similarity_annotation_running():
            return
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
        dialog = AutoAnnotationSettingsDialog(self._auto_ann_settings, parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._auto_ann_settings = dialog.settings()
        save_auto_annotation_settings(self._auto_ann_settings)
        self.statusMessage.emit("Auto annotation settings saved.")

    def _selected_auto_annotation_modes(self) -> list[AutoAnnotationMode]:
        return [mode for mode, check in self._auto_ann_checks.items() if check.isChecked()]

    def _open_similarity_annotation_settings(self) -> None:
        dialog = SimilarityAnnotationSettingsDialog(self._similarity_ann_settings, parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._similarity_ann_settings = dialog.settings()
        save_similarity_annotation_settings(self._similarity_ann_settings)
        self._map.set_show_similarity_database_name(
            self._similarity_ann_settings.show_database_name_on_track
        )
        self.statusMessage.emit("Similarity annotation settings saved.")

    def _current_similarity_annotation_settings(self) -> SimilarityAnnotationSettings:
        return self._similarity_ann_settings

    def _reflow_similarity_database_chips(self) -> None:
        if hasattr(self, "_similarity_db_flow"):
            self._similarity_db_flow._reflow()

    def _make_similarity_db_check(self, entry: UserDatabase, index: int) -> QCheckBox:
        color = similarity_database_color(entry.name, index)
        check = QCheckBox(entry.name, self._similarity_db_flow)
        check.setToolTip(f"{entry.display_type_label}\n{entry.path}")
        check.setChecked(True)
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

        self._push_feature_undo_snapshot()
        self._sim_ann_cancel.reset()
        self._set_similarity_annotation_enabled(False)
        labels = ", ".join(run.database.name for run in databases)
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
            self.statusMessage.emit("Similarity annotation finished with an unexpected result.")
            self._clear_annotation_job_status()
            return
        data = self._map.data()
        if data is None:
            self._clear_annotation_job_status()
            return
        self._map.set_features(list(result.features))
        self._clear_annotation_job_status()
        self.statusMessage.emit(f"{result.message} Click Save to write changes to disk.")

    def _on_similarity_annotation_error(self, message: str) -> None:
        self._clear_annotation_job_status()
        if "stopped by user" in message.lower() or "job stopped" in message.lower():
            self.statusMessage.emit("Similarity annotation stopped.")
            return
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
            if features:
                self.statusMessage.emit(f"{result.message} Click Save to write changes to disk.")
            else:
                self.statusMessage.emit(result.message)
            return

        if not isinstance(result, AutoAnnotationResult):
            self.statusMessage.emit("Auto annotation finished with an unexpected result.")
            self._clear_annotation_job_status()
            return

        data = self._map.data()
        if data is None:
            self._clear_annotation_job_status()
            return

        if not result.features:
            QMessageBox.information(self, "Auto annotation", result.message)
            self.statusMessage.emit(result.message)
            self._clear_annotation_job_status()
            return

        merged = merge_auto_annotation_features(data.features, result.features, result.mode)
        self._map.set_features(merged)
        self._clear_annotation_job_status()
        self.statusMessage.emit(f"{result.message} Click Save to write changes to disk.")

    def _on_auto_annotation_error(self, message: str) -> None:
        self._clear_annotation_job_status()
        if "stopped by user" in message.lower() or "job stopped" in message.lower():
            self.statusMessage.emit("Auto annotation stopped.")
            return
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
