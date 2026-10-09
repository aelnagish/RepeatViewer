"""Multi-sequence FASTA annotation overview (names | tracks, synced rows)."""

from __future__ import annotations

import math
import re
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
    QFont,
    QFontMetrics,
    QGuiApplication,
    QImage,
    QKeyEvent,
    QPageSize,
    QPainter,
    QPalette,
    QPdfWriter,
    QPen,
    QPixmap,
    QPolygonF,
    QWheelEvent,
)
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QScrollArea,
    QSizePolicy,
    QSlider,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.annotation_colors import similarity_database_color
from repeatviewer.backends.auto_annotation import AutoAnnotationMode
from repeatviewer.backends.auto_annotation_settings import (
    load_auto_annotation_settings,
    save_auto_annotation_settings,
)
from repeatviewer.backends.database_prepare import database_search_readiness
from repeatviewer.backends.job_cancel import JobCancellation
from repeatviewer.backends.sequence_features import (
    SequenceFeature,
    fasta_record_lengths,
    feature_display_pack_key,
    is_aligned_fasta,
    load_aligned_fasta_rows,
    load_all_features_by_id,
    load_features_by_id_from_json,
    ordered_feature_pack_keys,
    remount_parent_ids,
    resolve_features_gff_path,
    resolve_sidecar_path,
    save_annotations_multi,
)
from repeatviewer.backends.similarity_annotation import (
    SimilarityDatabaseRun,
    load_use_system_blast_from_settings,
)
from repeatviewer.backends.similarity_annotation_settings import (
    load_similarity_annotation_settings,
    save_similarity_annotation_settings,
)
from repeatviewer.backends.te_annotation_batch import (
    TeAnnotationOptions,
    run_te_annotation_on_fasta,
)
from repeatviewer.database_settings import load_similarity_annotation_databases
from repeatviewer.widgets.preview_controls import (
    OVERVIEW_TOOLBAR_HEIGHT,
    AnnotationZoomBar,
)
from repeatviewer.widgets.preview_zoom import (
    amino_acid_color_map,
    invalidate_residue_color_maps,
    nucleotide_base_color_map,
)
from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail
from repeatviewer.workers.base import Worker

_AUTO_ANNOTATION_MODES: tuple[tuple[str, AutoAnnotationMode], ...] = (
    ("LTR", AutoAnnotationMode.LTR),
    ("TIR", AutoAnnotationMode.TIR),
    ("Tandem", AutoAnnotationMode.TANDEM),
    ("ORF", AutoAnnotationMode.ORF),
)

_BAR_HEIGHT = 7
_FIND_HIT_FILL = QColor(245, 158, 11, 150)
_FIND_CURRENT_FILL = QColor(37, 99, 235, 220)
_FIND_HIT_BORDER = QColor(180, 83, 9, 230)
_BAR_TOP_PAD = 3
_FEAT_HEIGHT = 10
_FEAT_GAP = 2
_ROW_BOTTOM_PAD = 4
_MIN_ROW_HEIGHT = 18
_NAME_COLUMN_DEFAULT = 200
_NAME_COLUMN_MAX = 220
_NAME_PAD = 8
_SELECTED_BG_LIGHT = QColor(11, 92, 171, 36)
_SELECTED_BG_DARK = QColor(80, 160, 255, 48)
_DEFAULT_FEAT = QColor("#c0392b")
# Only guards division/painting against a zero scale.  The actual zoom-out
# limit is calculated per data set so the longest sequence always fits.
_ABS_MIN_PPB = 1e-12
_MAX_PPB = 4.0
_MSA_MAX_PPB = 16.0
_MSA_LETTER_MIN_PPB = 7.5
_MSA_ROW_HEIGHT = 18
_DEFAULT_MSA_FONT_PT = 4
_DEFAULT_MSA_SEQUENCE_DISTANCE = 1
_OVERVIEW_DISPLAY_SETTINGS_GROUP = "annotation_overview_display"
_FEATURE_UNDO_LIMIT = 32


def _clone_features_by_id(
    features_by_id: dict[str, list[SequenceFeature]],
) -> dict[str, list[SequenceFeature]]:
    return {
        rid: [SequenceFeature.from_dict(feat.to_dict()) for feat in feats]
        for rid, feats in features_by_id.items()
    }
_MSA_GAP_FILL = QColor(220, 220, 220, 90)
_REGION_BAND_LIGHT = QColor(11, 92, 171, 55)
_REGION_BAND_DARK = QColor(80, 160, 255, 70)
_MAX_ZOOM_OVER_FIT = 32.0
_QT_MAX_WIDGET_SIZE = 16_777_215
_BAR_RIGHT_PAD = 8
_ZOOM_STEP = 1.08
_SCALE_HEIGHT = 28
# Ignore tiny moves (click / accidental jitter); macOS-like drag threshold.
_RUBBER_MIN_PX = 14


def _guess_msa_is_protein(seqs: list[str]) -> bool:
    """Heuristic: protein when most alphabetic residues are not ACGTUN."""
    sample_parts: list[str] = []
    for seq in seqs[:30]:
        sample_parts.append(seq[:400])
        if sum(len(p) for p in sample_parts) >= 8_000:
            break
    letters = [c for c in "".join(sample_parts) if c.isalpha()]
    if not letters:
        return False
    nuc = sum(1 for c in letters if c in "ACGTUN")
    return (nuc / len(letters)) < 0.85


def _is_dark(widget: QWidget) -> bool:
    return widget.palette().color(QPalette.ColorRole.Window).lightness() < 128


def _wheel_delta_y(event: QWheelEvent) -> int:
    """Vertical pixels for a wheel/trackpad event (do not steal horizontal swipes)."""
    pixel = event.pixelDelta()
    if not pixel.isNull() and pixel.y() != 0:
        return int(pixel.y())
    angle = event.angleDelta()
    dy = int(angle.y())
    if dy == 0:
        return 0
    # Mouse wheel notches are typically ±120; map to ~one row, not 120px jumps.
    return int(round(dy * (_MIN_ROW_HEIGHT / 120.0)))


def _wheel_delta_x(event: QWheelEvent, *, scroll_step: int = 40) -> int:
    """Horizontal pixels for trackpad swipe (or Shift+wheel when not zooming)."""
    pixel = event.pixelDelta()
    if not pixel.isNull() and pixel.x() != 0:
        return int(pixel.x())
    angle = event.angleDelta()
    if angle.x() != 0:
        return int(round((angle.x() / 120.0) * max(1, scroll_step)))
    return 0


def _bar_color(widget: QWidget, *, selected: bool) -> QColor:
    if selected:
        return QColor("#6eb0ef") if _is_dark(widget) else QColor("#5b8fc7")
    # Soft slate — keep bars readable without dominating the features.
    return QColor("#9aa3ad") if _is_dark(widget) else QColor("#b0b6be")


def _split_tracks_by_source_pref() -> bool:
    value = QSettings().value("annotation/split_tracks_by_source", True)
    if isinstance(value, bool):
        return value
    if isinstance(value, str):
        return value.lower() in ("1", "true", "yes", "on")
    return bool(value)


def _pack_feature_tracks(
    features: list[SequenceFeature],
    *,
    split_by_source: bool | None = None,
) -> list[list[SequenceFeature]]:
    """Pack features into visual bands by (source, track) or track only.

    Within each band, overlapping features stack into additional rows.
    """
    if not features:
        return []
    if split_by_source is None:
        split_by_source = _split_tracks_by_source_pref()
    tracks: list[list[SequenceFeature]] = []
    keys = ordered_feature_pack_keys(features, split_by_source=split_by_source)
    for key in keys:
        group = [
            feat
            for feat in sorted(features, key=lambda f: (f.start, f.end, f.name))
            if feature_display_pack_key(feat, split_by_source=split_by_source) == key
        ]
        band_rows: list[list[SequenceFeature]] = []
        for feat in group:
            placed = False
            for row in band_rows:
                if all(feat.start > other.end or feat.end < other.start for other in row):
                    row.append(feat)
                    placed = True
                    break
            if not placed:
                band_rows.append([feat])
        tracks.extend(band_rows)
    return tracks


def _row_height(n_tracks: int) -> int:
    feat_block = 0
    if n_tracks > 0:
        feat_block = n_tracks * _FEAT_HEIGHT + (n_tracks - 1) * _FEAT_GAP + 2
    return max(
        _MIN_ROW_HEIGHT,
        _BAR_TOP_PAD + _BAR_HEIGHT + feat_block + _ROW_BOTTOM_PAD,
    )


def _draw_arrow_feature(
    painter: QPainter,
    rect: QRectF,
    color: QColor,
    *,
    forward: bool,
) -> None:
    if rect.width() < 1.5 or rect.height() < 2:
        return
    arrow_w = min(rect.height() * 0.85, max(3.0, rect.width() * 0.2))
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


def _tick_step_for_span(span_bp: float, bars_width_px: float, min_label_px: float) -> int:
    if bars_width_px <= 1 or span_bp <= 0:
        return max(1, int(span_bp) or 1)
    raw = max(1.0, float(span_bp) * float(min_label_px) / float(bars_width_px))
    exp = math.floor(math.log10(raw)) if raw > 0 else 0
    base = 10**exp
    for mult in (1, 2, 5, 10):
        step = int(mult * base)
        if step >= raw:
            return max(1, step)
    return max(1, int(10 * base))


def _elide_middle(fm: QFontMetrics, text: str, width: int) -> str:
    if width <= 0 or fm.horizontalAdvance(text) <= width:
        return text
    if width < fm.horizontalAdvance("…"):
        return ""
    ell = "…"
    left = max(1, len(text) // 2)
    right = left
    while left + right > 0:
        candidate = (
            text[:left] + ell + text[len(text) - right :] if right else text[:left] + ell
        )
        if fm.horizontalAdvance(candidate) <= width:
            return candidate
        if left >= right and left > 1:
            left -= 1
        elif right > 0:
            right -= 1
        else:
            left -= 1
    return ell


class _NamesColumn(QWidget):
    """Sequence names; row geometry matches the tracks canvas exactly."""

    rowClicked = Signal(int, object)  # row, keyboard modifiers
    rowsDragSelected = Signal(object)  # list[int]
    selectAllRequested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._record_ids: list[str] = []
        self._row_tops: list[int] = []
        self._row_heights: list[int] = []
        self._selected_rows: set[int] = set()
        self._label_point_size = _DEFAULT_MSA_FONT_PT
        self._drag_origin: QPointF | None = None
        self._dragging = False
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent, True)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Minimum)

    def set_label_point_size(self, point_size: int) -> None:
        """Name labels use DNA/protein size + 1 (clamped)."""
        self._label_point_size = max(2, min(36, int(point_size) + 1))
        self.update()

    def label_point_size(self) -> int:
        return int(self._label_point_size)

    def set_rows(
        self,
        record_ids: list[str],
        tops: list[int],
        heights: list[int],
        selected: set[int],
        *,
        width: int,
        content_height: int,
    ) -> None:
        self._record_ids = list(record_ids)
        self._row_tops = list(tops)
        self._row_heights = list(heights)
        self._selected_rows = set(selected)
        w = max(80, int(width))
        h = max(_MIN_ROW_HEIGHT, int(content_height))
        self.setMinimumWidth(80)
        self.setMinimumHeight(h)
        self.resize(w, h)
        self.update()

    def set_selected_rows(self, rows: set[int]) -> None:
        self._selected_rows = set(rows)
        self.update()

    def _label_font(self) -> QFont:
        """Name font = DNA/protein label size + 1."""
        font = QFont(self.font())
        font.setPointSize(max(2, int(self._label_point_size)))
        return font

    def _row_at_y(self, y: float) -> int:
        for index, (top, height) in enumerate(zip(self._row_tops, self._row_heights)):
            if top <= y < top + height:
                return index
        return -1

    def _rows_in_y_range(self, y0: float, y1: float) -> list[int]:
        lo, hi = (y0, y1) if y0 <= y1 else (y1, y0)
        rows: list[int] = []
        for index, (top, height) in enumerate(zip(self._row_tops, self._row_heights)):
            if top + height > lo and top < hi:
                rows.append(index)
        return rows

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._drag_origin = QPointF(event.position())
            self._dragging = False
            row = self._row_at_y(event.position().y())
            if row >= 0:
                self.rowClicked.emit(row, event.modifiers())
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if self._drag_origin is not None and event.buttons() & Qt.MouseButton.LeftButton:
            delta = event.position() - self._drag_origin
            if (
                abs(delta.x()) >= _RUBBER_MIN_PX
                or abs(delta.y()) >= _RUBBER_MIN_PX
            ):
                self._dragging = True
                rubber = QRect(
                    self._drag_origin.toPoint(), event.position().toPoint()
                ).normalized()
                rows = self._rows_in_y_range(rubber.top(), rubber.bottom())
                if rows:
                    self.rowsDragSelected.emit(rows)
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._drag_origin = None
            self._dragging = False
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802
        from PySide6.QtGui import QKeySequence

        if event.matches(QKeySequence.StandardKey.SelectAll) or (
            event.key() == Qt.Key.Key_A
            and event.modifiers()
            & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier)
        ):
            self.selectAllRequested.emit()
            event.accept()
            return
        super().keyPressEvent(event)

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        clip = event.rect()
        base = self.palette().color(QPalette.ColorRole.Base)
        painter.fillRect(clip, base)
        sel_bg = _SELECTED_BG_DARK if _is_dark(self) else _SELECTED_BG_LIGHT
        text_color = self.palette().color(QPalette.ColorRole.Text)
        mid = self.palette().color(QPalette.ColorRole.Mid)

        font = self._label_font()
        painter.setFont(font)
        fm = QFontMetrics(font)

        for index, rid in enumerate(self._record_ids):
            top = self._row_tops[index]
            height = self._row_heights[index]
            if top + height < clip.top() or top > clip.bottom():
                continue
            selected = index in self._selected_rows
            if selected:
                painter.fillRect(0, top, self.width(), height, sel_bg)
            label = _elide_middle(fm, rid, self.width() - 2 * _NAME_PAD)
            text_h = fm.height()
            # Center on the row band (MSA residue / bar strip), not the gap below.
            text_y = int(top + max(0, (height - text_h) / 2))
            painter.setPen(text_color)
            painter.drawText(
                QRect(_NAME_PAD, text_y, self.width() - 2 * _NAME_PAD, text_h),
                int(Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft),
                label,
            )

        painter.setPen(QPen(mid, 1))
        painter.drawLine(self.width() - 1, clip.top(), self.width() - 1, clip.bottom())
        painter.end()


class _BarsCanvas(QWidget):
    """Length bars + feature arrows, or aligned MSA residue rows (no name column)."""

    selectionChanged = Signal()
    featureActivated = Signal(int, int)
    zoomChanged = Signal(float)
    regionZoomRequested = Signal(float, float)  # kept for Fit/Select zoom helpers
    regionSelectionChanged = Signal(object)  # tuple[int,int] | None (1-based inclusive)
    rowsChanged = Signal()  # layout / selection updated — refresh names column

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._record_ids: list[str] = []
        self._lengths: dict[str, int] = {}
        self._features: dict[str, list[SequenceFeature]] = {}
        self._tracks: dict[str, list[list[SequenceFeature]]] = {}
        self._msa_seqs: list[str] = []
        self._msa_mode = False
        self._msa_protein = False
        self._msa_colors: dict[str, tuple[str, str]] = {}
        self._msa_qcolors: dict[str, tuple[QColor, QColor]] = {}
        self._msa_residue_height = _MSA_ROW_HEIGHT
        self._msa_font_pt = _DEFAULT_MSA_FONT_PT
        self._sequence_row_gap = 0  # pixels between sequences (Others → distance)
        self._show_annotation_track = True
        self._region_sel: tuple[int, int] | None = None
        self._row_tops: list[int] = []
        self._row_heights: list[int] = []
        self._max_length = 1
        self._selected_rows: set[int] = set()
        self._selected_feat: tuple[int, int] | None = None
        self._ppb = 0.05
        self._bar_pad = 8
        self._min_viewport = QSize(120, 40)
        self._drag_origin: QPointF | None = None
        self._rubber = QRect()
        self._rubber_active = False  # True only after drag exceeds threshold
        self._rubber_prev = QRect()  # last overlay rect (for cheap partial updates)
        self._drag_paint_cache: QPixmap | None = None
        self._pending_scroll_bp: float | None = None
        self._pending_anchor: tuple[float, float] | None = None
        # Edit > Find: hit ranges per record id, and the hit being shown.
        self._find_hits: dict[str, list[tuple[int, int]]] = {}
        self._find_current: tuple[str, int, int] | None = None
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Minimum)
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent, True)

    def is_msa_mode(self) -> bool:
        return bool(self._msa_mode)

    def set_data(
        self,
        record_ids: list[str],
        lengths: dict[str, int],
        features: dict[str, list[SequenceFeature]],
        *,
        fit_bars_width: int | None = None,
    ) -> None:
        self._msa_mode = False
        self._msa_seqs = []
        self._msa_colors = {}
        self._msa_qcolors = {}
        self._clear_drag_paint_cache()
        self._record_ids = list(record_ids)
        self._lengths = dict(lengths)
        self._features = dict(features)
        self._tracks = {
            rid: _pack_feature_tracks(list(features.get(rid) or []))
            for rid in record_ids
        }
        self._max_length = max((lengths.get(rid, 1) for rid in record_ids), default=1)
        self._selected_rows = set()
        self._selected_feat = None
        self._region_sel = None
        self._rubber = QRect()
        self._rubber_active = False
        self._drag_origin = None
        if fit_bars_width and fit_bars_width > 40 and self._max_length > 0:
            self._ppb = self.fit_pixels_per_base(fit_bars_width)
        fit_ppb, max_ppb = self.zoom_bounds()
        self._ppb = max(fit_ppb, min(max_ppb, self._ppb))
        self._recompute_layout()
        self.update()
        self.zoomChanged.emit(self._ppb)
        self.regionSelectionChanged.emit(None)
        self.rowsChanged.emit()

    def set_msa_data(
        self,
        record_ids: list[str],
        sequences: list[str],
        *,
        features: dict[str, list[SequenceFeature]] | None = None,
        fit_bars_width: int | None = None,
    ) -> None:
        """Show aligned residue rows (same length), optionally with feature tracks."""
        if len(record_ids) != len(sequences):
            raise ValueError("record_ids and sequences length mismatch")
        self._msa_mode = True
        self._msa_seqs = list(sequences)
        self._msa_protein = _guess_msa_is_protein(self._msa_seqs)
        self._msa_colors = (
            amino_acid_color_map()
            if self._msa_protein
            else nucleotide_base_color_map()
        )
        self._rebuild_msa_qcolor_cache()
        self._clear_drag_paint_cache()
        self._record_ids = list(record_ids)
        width = len(sequences[0]) if sequences else 1
        self._lengths = {rid: width for rid in record_ids}
        feat_map = dict(features or {})
        self._features = {
            rid: list(feat_map.get(rid) or []) for rid in record_ids
        }
        self._tracks = {
            rid: _pack_feature_tracks(list(self._features.get(rid) or []))
            for rid in record_ids
        }
        self._max_length = max(1, width)
        self._selected_rows = set()
        self._selected_feat = None
        self._region_sel = None
        self._rubber = QRect()
        self._rubber_active = False
        self._drag_origin = None
        if fit_bars_width and fit_bars_width > 40 and self._max_length > 0:
            self._ppb = self.fit_pixels_per_base(fit_bars_width)
        fit_ppb, max_ppb = self.zoom_bounds()
        self._ppb = max(fit_ppb, min(max_ppb, self._ppb))
        self._recompute_layout()
        self.update()
        self.zoomChanged.emit(self._ppb)
        self.regionSelectionChanged.emit(None)
        self.rowsChanged.emit()

    def set_features_in_place(
        self, features: dict[str, list[SequenceFeature]]
    ) -> None:
        """Update features without leaving MSA/bar mode or clearing selection."""
        self._features = {
            rid: list(features.get(rid) or []) for rid in self._record_ids
        }
        self._tracks = {
            rid: _pack_feature_tracks(list(self._features.get(rid) or []))
            for rid in self._record_ids
        }
        self._clear_drag_paint_cache()
        self._recompute_layout()
        self.update()

    def set_show_annotation_track(self, enabled: bool) -> None:
        self._show_annotation_track = bool(enabled)
        self._clear_drag_paint_cache()
        self._recompute_layout()
        self.update()

    def show_annotation_track(self) -> bool:
        return bool(self._show_annotation_track)

    def set_msa_residue_size(self, point_size: int) -> None:
        """DNA/protein label size → MSA residue band height and letter font.

        Height always matches the monospace glyph box (same idea as single-seq
        letter cells) — no fixed pixel floor.
        """
        pt = max(2, min(36, int(point_size)))
        self._msa_font_pt = pt
        font = QFont(self.font())
        font.setFamily("Menlo")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(pt)
        fm = QFontMetrics(font)
        glyph_h = max(fm.height(), fm.ascent() + fm.descent())
        self._msa_residue_height = max(1, int(glyph_h + 4))
        if self._msa_mode:
            self._clear_drag_paint_cache()
            self._recompute_layout()
            self.update()

    def msa_residue_height(self) -> int:
        return int(self._msa_residue_height)

    def set_sequence_row_gap(self, gap: int) -> None:
        """Vertical gap between sequences (Others → Annotation/Sequence distance)."""
        self._sequence_row_gap = max(0, min(40, int(gap)))
        self._recompute_layout()
        self.update()

    def sequence_row_gap(self) -> int:
        return int(self._sequence_row_gap)

    def region_selection(self) -> tuple[int, int] | None:
        return self._region_sel

    def set_region_selection(self, start: int | None, end: int | None = None) -> None:
        if start is None:
            if self._region_sel is None:
                return
            self._region_sel = None
            self.update()
            self.regionSelectionChanged.emit(None)
            return
        maximum = max(1, int(self._max_length))
        a = max(1, min(maximum, int(start)))
        b = max(1, min(maximum, int(end if end is not None else start)))
        if b < a:
            a, b = b, a
        new = (a, b)
        if new == self._region_sel:
            return
        self._region_sel = new
        self.update()
        self.regionSelectionChanged.emit(new)

    def select_all_rows(self) -> None:
        n = len(self._record_ids)
        self._selected_rows = set(range(n)) if n else set()
        self._selected_feat = None
        self.update()
        self.selectionChanged.emit()
        self.rowsChanged.emit()

    def refresh_msa_colors(self) -> None:
        """Reload Annotation residue palettes after Display color settings change."""
        if not self._msa_mode:
            return
        invalidate_residue_color_maps()
        self._msa_colors = (
            amino_acid_color_map()
            if self._msa_protein
            else nucleotide_base_color_map()
        )
        self._rebuild_msa_qcolor_cache()
        self._clear_drag_paint_cache()
        self.update()

    def _rebuild_msa_qcolor_cache(self) -> None:
        cache: dict[str, tuple[QColor, QColor]] = {}
        for key, pair in self._msa_colors.items():
            bg = QColor(pair[0])
            fg = QColor(pair[1])
            if bg.isValid() and fg.isValid():
                cache[key] = (bg, fg)
                cache[key.upper()] = (bg, fg)
                cache[key.lower()] = (bg, fg)
        self._msa_qcolors = cache

    def _clear_drag_paint_cache(self) -> None:
        self._drag_paint_cache = None
        self._rubber_prev = QRect()

    def _capture_drag_paint_cache(self) -> None:
        """Snapshot the MSA/bars under the rubber overlay for smooth dragging."""
        w = max(1, self.width())
        h = max(1, self.height())
        pm = QPixmap(w, h)
        pm.fill(self.palette().color(QPalette.ColorRole.Base))
        painter = QPainter(pm)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        clip = QRect(0, 0, w, h)
        if self._msa_mode:
            self._paint_msa_rows(painter, clip)
        else:
            self._paint_bar_rows(painter, clip)
        painter.end()
        self._drag_paint_cache = pm

    def _update_rubber_overlay(self, rubber: QRect) -> None:
        """Repaint only the band that changed while dragging a region."""
        old = QRect(self._rubber_prev)
        self._rubber = QRect(rubber)
        self._rubber_active = True
        new = QRect(rubber).normalized()
        dirty = old.united(new)
        # Full-height strips (dim overlay spans the viewport vertically).
        if dirty.isNull():
            self.update()
        else:
            left = max(0, dirty.left() - 2)
            right = min(self.width(), dirty.right() + 2)
            self.update(QRect(left, 0, max(1, right - left), self.height()))
        self._rubber_prev = QRect(new)

    def pixels_per_base(self) -> float:
        return self._ppb

    def max_length(self) -> int:
        return self._max_length

    def bar_pad(self) -> int:
        return self._bar_pad

    def row_tops(self) -> list[int]:
        return list(self._row_tops)

    def row_heights(self) -> list[int]:
        return list(self._row_heights)

    def content_height(self) -> int:
        if not self._row_tops:
            return _MIN_ROW_HEIGHT
        return self._row_tops[-1] + self._row_heights[-1] + 4

    def record_ids(self) -> list[str]:
        return list(self._record_ids)

    def selected_rows(self) -> list[int]:
        return sorted(self._selected_rows)

    def selected_row_set(self) -> set[int]:
        return set(self._selected_rows)

    def selected_feature(self) -> tuple[int, int] | None:
        """``(row, flat_feature_index)`` or ``None``."""
        return self._selected_feat

    def clear_feature_selection(self) -> None:
        if self._selected_feat is None:
            return
        self._selected_feat = None
        self.update()

    def set_viewport_min_size(self, width: int, height: int) -> None:
        self._min_viewport = QSize(max(40, int(width)), max(40, int(height)))
        fit_ppb, max_ppb = self.zoom_bounds()
        self._ppb = max(fit_ppb, min(max_ppb, self._ppb))
        self._recompute_layout()

    def fit_pixels_per_base(self, bars_width: int | None = None) -> float:
        """Scale that fits the longest sequence inside the track viewport."""
        width = max(
            40,
            int(self._min_viewport.width() if bars_width is None else bars_width),
        )
        usable_width = max(1, width - self._bar_pad - _BAR_RIGHT_PAD)
        return max(_ABS_MIN_PPB, usable_width / max(1, self._max_length))

    def zoom_bounds(self, bars_width: int | None = None) -> tuple[float, float]:
        """Return data-dependent zoom limits that remain safe for QWidget."""
        fit_ppb = self.fit_pixels_per_base(bars_width)
        geometry_limit = max(
            _ABS_MIN_PPB,
            (
                _QT_MAX_WIDGET_SIZE
                - self._bar_pad
                - _BAR_RIGHT_PAD
            )
            / max(1, self._max_length),
        )
        ceil = _MSA_MAX_PPB if self._msa_mode else _MAX_PPB
        desired_max = max(ceil, fit_ppb * _MAX_ZOOM_OVER_FIT)
        return fit_ppb, max(fit_ppb, min(desired_max, geometry_limit))

    def set_pixels_per_base(self, ppb: float) -> None:
        fit_ppb, max_ppb = self.zoom_bounds()
        self._ppb = max(fit_ppb, min(max_ppb, float(ppb)))
        self._clear_drag_paint_cache()
        self._recompute_layout()
        self.update()
        self.zoomChanged.emit(self._ppb)

    def zoom_by(self, factor: float, *, anchor_x: float | None = None) -> None:
        old = self._ppb
        fit_ppb, max_ppb = self.zoom_bounds()
        new = max(fit_ppb, min(max_ppb, old * factor))
        if abs(new - old) < 1e-12:
            return
        if anchor_x is not None and old > 0:
            self._pending_anchor = (
                (float(anchor_x) - self._bar_pad) / old,
                float(anchor_x),
            )
        else:
            self._pending_anchor = None
        self._pending_scroll_bp = None
        self._ppb = new
        self._clear_drag_paint_cache()
        self._recompute_layout()
        self.update()
        self.zoomChanged.emit(self._ppb)

    def fit_to_bars_width(self, bars_width: int) -> None:
        if self._max_length <= 0:
            return
        self._pending_scroll_bp = None
        self._pending_anchor = None
        self.set_pixels_per_base(self.fit_pixels_per_base(bars_width))

    def zoom_to_bp_range(
        self, start_bp: float, end_bp: float, bars_width: int
    ) -> None:
        span = max(1.0, float(end_bp) - float(start_bp))
        width = max(40, int(bars_width))
        self._pending_scroll_bp = float(start_bp)
        self._pending_anchor = None
        usable_width = max(1, width - self._bar_pad - _BAR_RIGHT_PAD)
        self.set_pixels_per_base(usable_width / span)

    def take_pending_scroll_bp(self) -> float | None:
        value = self._pending_scroll_bp
        self._pending_scroll_bp = None
        return value

    def take_pending_anchor(self) -> tuple[float, float] | None:
        value = self._pending_anchor
        self._pending_anchor = None
        return value

    def set_find_hits(
        self,
        hits_by_id: dict[str, list[tuple[int, int]]],
        current: tuple[str, int, int] | None = None,
    ) -> None:
        """Edit > Find hits drawn on the bars ({} clears)."""
        self._find_hits = {rid: list(ranges) for rid, ranges in (hits_by_id or {}).items()}
        self._find_current = current
        self.update()

    def set_selected_rows(self, rows: list[int] | set[int]) -> None:
        self._selected_rows = {int(r) for r in rows if 0 <= int(r) < len(self._record_ids)}
        self._selected_feat = None
        self.update()
        self.selectionChanged.emit()
        self.rowsChanged.emit()

    def _msa_track_block_height(self, rid: str) -> int:
        if not self._show_annotation_track:
            return 0
        n_tracks = len(self._tracks.get(rid) or [])
        if n_tracks <= 0:
            return 0
        return n_tracks * _FEAT_HEIGHT + (n_tracks - 1) * _FEAT_GAP + 2

    def _recompute_layout(self) -> None:
        self._row_tops = []
        self._row_heights = []
        y = 0
        n_rows = len(self._record_ids)
        gap = max(0, int(self._sequence_row_gap))
        for i, rid in enumerate(self._record_ids):
            if self._msa_mode:
                # Residue band tracks letter size; no _MIN_ROW_HEIGHT floor.
                h = int(self._msa_residue_height) + self._msa_track_block_height(rid)
            else:
                n_tracks = (
                    len(self._tracks.get(rid) or [])
                    if self._show_annotation_track
                    else 0
                )
                h = _row_height(n_tracks)
            self._row_tops.append(y)
            self._row_heights.append(h)
            y += h
            # Empty space between sequences (not part of the row hit-target).
            if i < n_rows - 1:
                y += gap
        content_w = math.ceil(
            self._bar_pad + self._max_length * self._ppb + _BAR_RIGHT_PAD
        )
        # Extra bottom pad so the last row clears the horizontal scrollbar edge.
        y = y + 4 if y else _MIN_ROW_HEIGHT
        w = min(
            _QT_MAX_WIDGET_SIZE,
            max(120, content_w, self._min_viewport.width()),
        )
        h = min(
            _QT_MAX_WIDGET_SIZE,
            max(_MIN_ROW_HEIGHT, y, self._min_viewport.height()),
        )
        self.setMinimumSize(w, h)
        self.resize(w, h)
        self.updateGeometry()
        self.rowsChanged.emit()

    def x_to_bp(self, x: float) -> float:
        return max(0.0, (float(x) - self._bar_pad) / max(self._ppb, 1e-12))

    def bp_to_x(self, bp: float) -> float:
        return self._bar_pad + float(bp) * self._ppb

    def _row_at_y(self, y: float) -> int:
        for index, (top, height) in enumerate(zip(self._row_tops, self._row_heights)):
            if top <= y < top + height:
                return index
        return -1

    def _rows_in_y_range(self, y0: float, y1: float) -> list[int]:
        lo, hi = (y0, y1) if y0 <= y1 else (y1, y0)
        rows: list[int] = []
        for index, (top, height) in enumerate(zip(self._row_tops, self._row_heights)):
            if top + height > lo and top < hi:
                rows.append(index)
        return rows

    def _feature_band_y0(self, row: int) -> int:
        top = self._row_tops[row]
        if self._msa_mode:
            return top + int(self._msa_residue_height) + 1
        return top + _BAR_TOP_PAD + _BAR_HEIGHT + 2

    def _hit_feature(self, row: int, x: float, y: float) -> int | None:
        if not self._show_annotation_track:
            return None
        if not (0 <= row < len(self._record_ids)):
            return None
        rid = self._record_ids[row]
        length = max(1, int(self._lengths.get(rid, 1)))
        feat_y0 = self._feature_band_y0(row)
        flat_index = 0
        for track_i, track in enumerate(self._tracks.get(rid) or []):
            ty = feat_y0 + track_i * (_FEAT_HEIGHT + _FEAT_GAP)
            for feat in track:
                start = max(1, min(int(feat.start), length))
                end = max(1, min(int(feat.end), length))
                if start > end:
                    start, end = end, start
                x0 = self._bar_pad + (start - 1) * self._ppb
                x1 = self._bar_pad + end * self._ppb
                if x0 <= x <= max(x0 + 2, x1) and ty <= y <= ty + _FEAT_HEIGHT:
                    return flat_index
                flat_index += 1
        return None

    def _select_row_click(self, row: int, modifiers: Qt.KeyboardModifier) -> None:
        if row < 0:
            return
        if modifiers & (
            Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier
        ):
            if row in self._selected_rows:
                self._selected_rows.discard(row)
            else:
                self._selected_rows.add(row)
        elif modifiers & Qt.KeyboardModifier.ShiftModifier and self._selected_rows:
            anchor = min(self._selected_rows)
            lo, hi = (anchor, row) if anchor <= row else (row, anchor)
            self._selected_rows = set(range(lo, hi + 1))
        else:
            self._selected_rows = {row}
        self.update()
        self.selectionChanged.emit()
        self.rowsChanged.emit()

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.position()
            self._drag_origin = QPointF(pos)
            self._rubber = QRect()
            self._rubber_active = False
            row = self._row_at_y(pos.y())
            if row >= 0:
                feat_i = self._hit_feature(row, pos.x(), pos.y())
                self._selected_feat = (row, feat_i) if feat_i is not None else None
                self._select_row_click(row, event.modifiers())
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if self._drag_origin is not None and event.buttons() & Qt.MouseButton.LeftButton:
            origin = self._drag_origin.toPoint()
            current = event.position().toPoint()
            rubber = QRect(origin, current).normalized()
            # Only activate region UI after a real drag distance.
            if rubber.width() >= _RUBBER_MIN_PX or rubber.height() >= _RUBBER_MIN_PX:
                if not self._rubber_active:
                    # Snapshot once; drag overlays blit the cache (smooth).
                    self._capture_drag_paint_cache()
                    self._selected_feat = None
                # Defer row/name sync until release — only update the overlay.
                rows = self._rows_in_y_range(rubber.top(), rubber.bottom())
                if rows:
                    self._selected_rows = set(rows)
                self._update_rubber_overlay(rubber)
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton and self._drag_origin is not None:
            rubber = QRect(self._rubber)
            active = self._rubber_active
            self._drag_origin = None
            self._rubber = QRect()
            self._rubber_active = False
            self._clear_drag_paint_cache()
            if active and (
                rubber.width() >= _RUBBER_MIN_PX or rubber.height() >= _RUBBER_MIN_PX
            ):
                rows = self._rows_in_y_range(rubber.top(), rubber.bottom())
                if rows:
                    self._selected_rows = set(rows)
                    self._selected_feat = None
                    self.selectionChanged.emit()
                    self.rowsChanged.emit()
                # Drag selects a column/bp region (does not zoom).
                start_bp = self.x_to_bp(rubber.left())
                end_bp = self.x_to_bp(rubber.right())
                if end_bp < start_bp:
                    start_bp, end_bp = end_bp, start_bp
                start_i = max(1, int(math.floor(start_bp)) + 1)
                end_i = max(start_i, int(math.ceil(end_bp)))
                self.set_region_selection(start_i, end_i)
            elif self._region_sel is not None:
                # Plain click (inside or outside the blue band) clears the region.
                self.set_region_selection(None)
            self.update()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802
        from PySide6.QtGui import QKeySequence

        if event.key() == Qt.Key.Key_Escape and self._region_sel is not None:
            self.set_region_selection(None)
            event.accept()
            return
        if event.matches(QKeySequence.StandardKey.SelectAll) or (
            event.key() == Qt.Key.Key_A
            and event.modifiers()
            & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier)
        ):
            self.select_all_rows()
            event.accept()
            return
        super().keyPressEvent(event)

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.position()
            row = self._row_at_y(pos.y())
            if row >= 0:
                feat_i = self._hit_feature(row, pos.x(), pos.y())
                if feat_i is not None:
                    self.featureActivated.emit(row, feat_i)
                    event.accept()
                    return
        super().mouseDoubleClickEvent(event)

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802
        mods = event.modifiers() | QGuiApplication.keyboardModifiers()
        if mods & Qt.KeyboardModifier.ShiftModifier:
            delta = event.angleDelta().y()
            if delta == 0:
                delta = event.angleDelta().x()
            pixel = event.pixelDelta()
            if delta == 0 and not pixel.isNull():
                delta = int(pixel.y()) or int(pixel.x())
            if delta != 0:
                factor = _ZOOM_STEP if delta > 0 else 1 / _ZOOM_STEP
                self.zoom_by(factor, anchor_x=event.position().x())
            event.accept()
            return
        super().wheelEvent(event)

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        clip = event.rect()
        base = self.palette().color(QPalette.ColorRole.Base)
        # During rubber-band drag, blit the cached MSA and only redraw the overlay.
        cache = self._drag_paint_cache
        if self._rubber_active and cache is not None and not cache.isNull():
            painter.drawPixmap(clip, cache, clip)
        else:
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, not self._msa_mode)
            painter.fillRect(clip, base)
            if self._msa_mode:
                self._paint_msa_rows(painter, clip)
            else:
                self._paint_bar_rows(painter, clip)

        # Region overlay only after a real drag (not on click).
        if self._rubber_active and not self._rubber.isNull():
            x0 = self._rubber.left()
            x1 = max(x0 + 1, self._rubber.right())
            dim = QColor(0, 0, 0, 90) if _is_dark(self) else QColor(40, 40, 40, 55)
            if x0 > 0:
                painter.fillRect(QRect(0, clip.top(), x0, clip.height()), dim)
            if x1 < self.width():
                painter.fillRect(
                    QRect(x1, clip.top(), self.width() - x1, clip.height()), dim
                )
            edge = QColor("#4ea1ff") if _is_dark(self) else QColor("#0b5cab")
            painter.setPen(QPen(edge, 1.5, Qt.PenStyle.DashLine))
            painter.drawLine(x0, clip.top(), x0, clip.bottom())
            painter.drawLine(x1, clip.top(), x1, clip.bottom())
        painter.end()

    def _paint_region_band(self, painter: QPainter, clip: QRect) -> None:
        sel = self._region_sel
        if sel is None:
            return
        start, end = sel
        x0 = self.bp_to_x(start - 1)
        x1 = self.bp_to_x(end)
        if x1 < clip.left() or x0 > clip.right():
            return
        band = _REGION_BAND_DARK if _is_dark(self) else _REGION_BAND_LIGHT
        painter.fillRect(
            QRectF(x0, clip.top(), max(1.0, x1 - x0), float(clip.height())),
            band,
        )

    def _paint_feature_tracks(
        self,
        painter: QPainter,
        *,
        row: int,
        rid: str,
        length: int,
        feat_y0: int,
    ) -> None:
        if not self._show_annotation_track:
            return
        flat_index = 0
        for track_i, track in enumerate(self._tracks.get(rid) or []):
            ty = feat_y0 + track_i * (_FEAT_HEIGHT + _FEAT_GAP)
            for feat in track:
                start = max(1, min(int(feat.start), length))
                end = max(1, min(int(feat.end), length))
                if start > end:
                    start, end = end, start
                x0 = self._bar_pad + (start - 1) * self._ppb
                x1 = self._bar_pad + end * self._ppb
                rect = QRectF(x0, ty, max(2.0, x1 - x0), float(_FEAT_HEIGHT))
                try:
                    fc = QColor(feat.color) if feat.color else _DEFAULT_FEAT
                    if not fc.isValid():
                        fc = _DEFAULT_FEAT
                except Exception:
                    fc = _DEFAULT_FEAT
                _draw_arrow_feature(
                    painter, rect, fc, forward=feat.is_forward()
                )
                if self._selected_feat == (row, flat_index):
                    painter.setBrush(Qt.BrushStyle.NoBrush)
                    painter.setPen(
                        QPen(QColor("#4ea1ff" if _is_dark(self) else "#0b5cab"), 2)
                    )
                    painter.drawRect(rect.adjusted(-1, -1, 1, 1))
                flat_index += 1

    def _paint_msa_rows(self, painter: QPainter, clip: QRect) -> None:
        """Column-aligned residue cells (visible rows/columns only)."""
        sel_bg = _SELECTED_BG_DARK if _is_dark(self) else _SELECTED_BG_LIGHT
        ppb = max(self._ppb, 1e-12)
        show_letters = ppb >= _MSA_LETTER_MIN_PPB and ppb >= 6
        c0 = max(0, int((clip.left() - self._bar_pad) / ppb) - 1)
        c1 = min(self._max_length, int((clip.right() - self._bar_pad) / ppb) + 2)
        if show_letters:
            font = QFont(self.font())
            font.setFamily("Menlo")
            font.setStyleHint(QFont.StyleHint.Monospace)
            font.setPointSize(max(2, min(36, int(self._msa_font_pt))))
            painter.setFont(font)
        gap_fill = (
            QColor(70, 70, 70, 70) if _is_dark(self) else _MSA_GAP_FILL
        )
        qcolors = self._msa_qcolors
        # Colored cell height matches letter metrics (set_msa_residue_size).
        residue_h = max(1, int(self._msa_residue_height))
        painter.setPen(Qt.PenStyle.NoPen)
        for index, rid in enumerate(self._record_ids):
            top = self._row_tops[index]
            height = self._row_heights[index]
            if top + height < clip.top() or top > clip.bottom():
                continue
            if index in self._selected_rows:
                painter.fillRect(0, top, self.width(), height, sel_bg)
            seq = self._msa_seqs[index] if index < len(self._msa_seqs) else ""
            y = top
            # Coalesce same-residue runs when zoomed out (far fewer draw calls).
            j = c0
            while j < c1:
                x = self._bar_pad + j * ppb
                if x + ppb < clip.left():
                    j += 1
                    continue
                if x > clip.right():
                    break
                ch = seq[j] if j < len(seq) else "-"
                run_end = j + 1
                if not show_letters:
                    while run_end < c1:
                        nxt = seq[run_end] if run_end < len(seq) else "-"
                        if nxt != ch:
                            break
                        run_end += 1
                span = run_end - j
                cell = QRectF(x, y, max(0.5, ppb * span), float(residue_h))
                if ch in "-.":
                    painter.fillRect(cell, gap_fill)
                else:
                    pair = qcolors.get(ch)
                    if pair is not None:
                        painter.fillRect(cell, pair[0])
                        if show_letters:
                            painter.setPen(pair[1])
                            painter.drawText(
                                cell.toRect(),
                                int(Qt.AlignmentFlag.AlignCenter),
                                ch,
                            )
                            painter.setPen(Qt.PenStyle.NoPen)
                j = run_end
            hits = self._find_hits.get(rid)
            if hits:
                painter.setPen(QPen(_FIND_HIT_BORDER, 1))
                current = self._find_current
                for hit_start, hit_end in hits:
                    x0 = self._bar_pad + (hit_start - 1) * ppb
                    x1 = self._bar_pad + hit_end * ppb
                    mark = QRectF(x0, y, max(2.0, x1 - x0), float(residue_h))
                    is_current = (
                        current is not None and current == (rid, hit_start, hit_end)
                    )
                    painter.fillRect(
                        mark, _FIND_CURRENT_FILL if is_current else _FIND_HIT_FILL
                    )
                    painter.drawRect(mark)
                painter.setPen(Qt.PenStyle.NoPen)
            length = max(1, int(self._lengths.get(rid, 1)))
            self._paint_feature_tracks(
                painter,
                row=index,
                rid=rid,
                length=length,
                feat_y0=self._feature_band_y0(index),
            )
        self._paint_region_band(painter, clip)

    def _paint_bar_rows(self, painter: QPainter, clip: QRect) -> None:
        sel_bg = _SELECTED_BG_DARK if _is_dark(self) else _SELECTED_BG_LIGHT

        for index, rid in enumerate(self._record_ids):
            top = self._row_tops[index]
            height = self._row_heights[index]
            if top + height < clip.top() or top > clip.bottom():
                continue
            selected = index in self._selected_rows
            if selected:
                painter.fillRect(0, top, self.width(), height, sel_bg)

            length = max(1, int(self._lengths.get(rid, 1)))
            bar_w = max(2.0, length * self._ppb)
            bar_y = top + _BAR_TOP_PAD
            painter.fillRect(
                QRectF(self._bar_pad, bar_y, bar_w, _BAR_HEIGHT),
                _bar_color(self, selected=selected),
            )
            hits = self._find_hits.get(rid)
            if hits:
                current = self._find_current
                painter.setPen(QPen(_FIND_HIT_BORDER, 1))
                for hit_start, hit_end in hits:
                    x0 = self._bar_pad + (hit_start - 1) * self._ppb
                    x1 = self._bar_pad + hit_end * self._ppb
                    mark = QRectF(x0, bar_y - 2, max(3.0, x1 - x0), _BAR_HEIGHT + 4)
                    is_current = current is not None and current == (rid, hit_start, hit_end)
                    painter.fillRect(mark, _FIND_CURRENT_FILL if is_current else _FIND_HIT_FILL)
                    painter.drawRect(mark)
                painter.setPen(Qt.PenStyle.NoPen)

            self._paint_feature_tracks(
                painter,
                row=index,
                rid=rid,
                length=length,
                feat_y0=self._feature_band_y0(index),
            )
        self._paint_region_band(painter, clip)


class _ScaleRuler(QWidget):
    """Horizontal bp scale under the tracks only."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._ppb = 0.05
        self._max_length = 1
        self._bar_pad = 8
        self._scroll_x = 0
        self.setFixedHeight(_SCALE_HEIGHT)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_scale(self, *, ppb: float, max_length: int, bar_pad: int) -> None:
        self._ppb = max(_ABS_MIN_PPB, float(ppb))
        self._max_length = max(1, int(max_length))
        self._bar_pad = int(bar_pad)
        self.update()

    def set_scroll_x(self, value: int) -> None:
        self._scroll_x = int(value)
        self.update()

    def paintEvent(self, _event) -> None:  # noqa: N802
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
        base = self.palette().color(QPalette.ColorRole.Base)
        painter.fillRect(self.rect(), base)
        fg = self.palette().color(QPalette.ColorRole.Text)
        if _is_dark(self):
            fg = QColor("#c8c8c8")
        font = QFont(self.font())
        font.setPointSize(max(8, font.pointSize() - 1))
        painter.setFont(font)
        fm = QFontMetrics(font)

        width = self.width()
        bars_left = 0
        bars_right = width
        bars_width = max(1.0, float(width))
        origin = self._bar_pad
        ppb = max(self._ppb, 1e-12)

        start_bp = max(0.0, (self._scroll_x - origin) / ppb)
        end_bp = min(float(self._max_length), start_bp + bars_width / ppb)
        span = max(1.0, end_bp - start_bp)
        sample = f"{int(end_bp):,}"
        min_label_px = fm.horizontalAdvance(sample) + 12
        step = _tick_step_for_span(span, bars_width, float(min_label_px))

        baseline = self.height() - 4
        painter.setPen(QPen(fg, 1))
        painter.drawLine(0, baseline, width, baseline)

        first = int(start_bp // step) * step
        last_label_right = -1.0
        bp = first
        guard = 0
        while bp <= end_bp + step and guard < 500:
            guard += 1
            x = origin + bp * ppb - self._scroll_x
            if -2 <= x <= width + 2:
                painter.drawLine(int(x), baseline - 6, int(x), baseline)
                label = f"{int(bp):,}"
                label_w = fm.horizontalAdvance(label)
                # Prefer centering; left-align near the left edge so "0" is visible.
                label_left = x - label_w / 2
                if label_left < bars_left + 2:
                    label_left = bars_left + 2
                label_right = label_left + label_w
                if label_right > bars_right - 2:
                    label_left = bars_right - 2 - label_w
                    label_right = label_left + label_w
                if label_left >= last_label_right + 6 and label_right <= bars_right:
                    painter.drawText(int(label_left), baseline - 8, label)
                    last_label_right = label_right
            bp += step
        painter.end()


class MultiSequenceOverviewPanel(QWidget):
    """Annotation overview: names (left) | tracks + scale (right)."""

    statusMessage = Signal(str)
    annotationsSaved = Signal(object)
    sequenceExtracted = Signal(object)

    def __init__(self, fasta_path: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self._path = fasta_path.expanduser().resolve()
        self._record_ids: list[str] = []
        self._lengths: dict[str, int] = {}
        self._features: dict[str, list[SequenceFeature]] = {}
        self._did_initial_fit = False
        self._name_width = _NAME_COLUMN_DEFAULT
        self._syncing_scroll = False
        self._search_needle = ""
        self._search_hits: list[int] = []
        self._search_hit_index = -1
        self._updating_zoom_slider = False

        self._names = _NamesColumn(self)
        self._names_scroll = QScrollArea(self)
        self._names_scroll.setWidget(self._names)
        self._names_scroll.setWidgetResizable(False)
        self._names_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        self._names_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        # Tracks own the vertical scrollbar; names stay pixel-synced.
        self._names_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self._names_scroll.setMinimumWidth(120)

        self._bars = _BarsCanvas(self)
        self._bars_scroll = QScrollArea(self)
        self._bars_scroll.setWidget(self._bars)
        self._bars_scroll.setWidgetResizable(False)
        self._bars_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        # Always reserve the horizontal scrollbar so viewport height matches names.
        self._bars_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOn
        )
        self._bars_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )

        self._scale = _ScaleRuler(self)

        # Spacer under names matching the tracks h-scrollbar so last rows stay aligned.
        self._hbar_shim = QWidget(self)
        self._hbar_shim.setFixedHeight(
            max(12, self._bars_scroll.horizontalScrollBar().sizeHint().height())
        )
        self._scale_shim = QWidget(self)
        self._scale_shim.setFixedHeight(_SCALE_HEIGHT)

        base = self.palette().color(QPalette.ColorRole.Base)
        for area in (self._names_scroll, self._bars_scroll):
            area.setStyleSheet(
                f"QScrollArea {{ background: {base.name()}; border: none; }}"
            )
            area.viewport().setAutoFillBackground(True)
            vp_pal = area.viewport().palette()
            vp_pal.setColor(QPalette.ColorRole.Window, base)
            area.viewport().setPalette(vp_pal)

        self._names_scroll.verticalScrollBar().valueChanged.connect(
            self._on_names_vscroll
        )
        self._bars_scroll.verticalScrollBar().valueChanged.connect(
            self._on_bars_vscroll
        )
        self._bars_scroll.horizontalScrollBar().valueChanged.connect(self._on_hscroll)
        self._bars.rowsChanged.connect(self._sync_names_column)
        self._bars.featureActivated.connect(self._on_feature_activated)
        self._bars.zoomChanged.connect(self._on_zoom_changed)
        self._bars.regionZoomRequested.connect(self._on_region_zoom)
        self._bars.regionSelectionChanged.connect(self._on_region_selection_changed)
        self._names.rowClicked.connect(self._on_name_row_clicked)
        self._names.rowsDragSelected.connect(self._on_name_rows_dragged)
        self._names.selectAllRequested.connect(self._select_all_sequences)
        # Bars owns the vertical scrollbar; names stays synced. Wheel over either
        # pane uses pixelDelta (trackpad) or a gentle mouse-wheel step.
        self._names_scroll.viewport().installEventFilter(self)
        self._bars_scroll.viewport().installEventFilter(self)
        self._bars_scroll.verticalScrollBar().setSingleStep(18)
        self._bars_scroll.verticalScrollBar().setPageStep(120)

        # Zoom / Position / Select / Fit live on AnnotationZoomBar (same row as
        # Text / Annotation). Keep hidden locals until bind, and for standalone use.
        self._header_zoom_bar = None
        self._zoom_slider = QSlider(Qt.Orientation.Horizontal, self)
        self._zoom_slider.setObjectName("annotationZoomSlider")
        self._zoom_slider.setRange(0, 100)
        self._zoom_slider.setValue(0)
        self._zoom_slider.hide()
        self._zoom_slider.valueChanged.connect(self._on_zoom_slider_changed)
        self._name_filter = QLineEdit(self)
        self._name_filter.setObjectName("FileFilter")
        self._name_filter.setPlaceholderText("Filter by name…")
        self._name_filter.hide()
        self._name_filter.textChanged.connect(self._on_name_filter_edited)
        self._name_filter.returnPressed.connect(self._on_name_filter_next)

        # Host for AnnotationZoomBar — above the sequence tracks only (same
        # left edge as single-seq AnnotationMapToolbar over the map).
        self._toolbar_slot = QWidget(self)
        self._toolbar_slot.setObjectName("OverviewMapToolbarSlot")
        self._toolbar_slot.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._toolbar_slot.setFixedHeight(OVERVIEW_TOOLBAR_HEIGHT)
        self._toolbar_slot_layout = QHBoxLayout(self._toolbar_slot)
        self._toolbar_slot_layout.setContentsMargins(0, 0, 0, 0)
        self._toolbar_slot_layout.setSpacing(0)
        self._toolbar_slot.hide()

        # Keep name rows level with bars when the map toolbar is shown.
        self._toolbar_shim = QWidget(self)
        self._toolbar_shim.setFixedHeight(OVERVIEW_TOOLBAR_HEIGHT)
        self._toolbar_shim.hide()

        names_col = QWidget(self)
        names_layout = QVBoxLayout(names_col)
        names_layout.setContentsMargins(0, 0, 0, 0)
        names_layout.setSpacing(0)
        names_layout.addWidget(self._toolbar_shim, 0)
        names_layout.addWidget(self._names_scroll, 1)
        names_layout.addWidget(self._hbar_shim)
        names_layout.addWidget(self._scale_shim)

        tracks_col = QWidget(self)
        tracks_layout = QVBoxLayout(tracks_col)
        tracks_layout.setContentsMargins(0, 0, 0, 0)
        tracks_layout.setSpacing(0)
        tracks_layout.addWidget(self._toolbar_slot, 0)
        tracks_layout.addWidget(self._bars_scroll, 1)
        tracks_layout.addWidget(self._scale)

        self._splitter = QSplitter(Qt.Orientation.Horizontal, self)
        self._splitter.addWidget(names_col)
        self._splitter.addWidget(tracks_col)
        self._splitter.setStretchFactor(0, 0)
        self._splitter.setStretchFactor(1, 1)
        self._splitter.setChildrenCollapsible(False)
        self._splitter.splitterMoved.connect(self._on_splitter_moved)

        self._auto_ann_settings = load_auto_annotation_settings()
        self._similarity_ann_settings = load_similarity_annotation_settings()
        self._ann_cancel = JobCancellation()
        self._ann_thread: QThread | None = None
        self._ann_worker: Worker | None = None
        self._ann_job_kind = ""
        self._ann_progress_total = 0
        self._ann_last_sequence_status = ""
        self._auto_all_updating = False
        self._similarity_all_updating = False
        self._similarity_db_rows: list[tuple] = []
        self._feature_undo_stack: list[dict[str, list[SequenceFeature]]] = []

        self._rail = SequenceAnnotationRail(self)
        self._rail.set_overview_display_mode(True)
        self._rail.residueColorsChanged.connect(self._on_residue_colors_changed)
        self._rail.colorsClicked.connect(self._on_colors_setting_clicked)
        self._transfer_cancel = JobCancellation()
        self._transfer_thread: QThread | None = None
        self._transfer_worker: Worker | None = None
        self._wire_overview_annotation_rail()
        self._bars.set_show_annotation_track(
            bool(self._rail.show_annotation_track.isChecked())
        )

        # Same 1px map↔rail divider as single-seq AnnotationMapWindow.
        from repeatviewer.widgets.sequence_map_viewer import _MAP_RAIL_DIVIDER

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

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._splitter, 1)
        layout.addWidget(map_rail_divider, 0)
        layout.addWidget(self._rail, 0)

    def _compute_name_width(self, viewport_width: int | None = None) -> int:
        """Keep names compact so the MSA + right control panel get more width."""
        font = QFont(self.font())
        font.setPointSize(max(2, int(self._names.label_point_size())))
        fm = QFontMetrics(font)
        widest = 0
        for rid in self._record_ids[:500]:
            widest = max(widest, fm.horizontalAdvance(rid))
        natural = widest + 2 * _NAME_PAD + 4
        cap = _NAME_COLUMN_MAX
        if viewport_width and viewport_width > 100:
            # At most ~22% of the overview content (names | tracks), not the rail.
            cap = min(cap, max(120, int(viewport_width * 0.22)))
        return max(120, min(cap, natural))

    def _current_name_width(self) -> int:
        sizes = self._splitter.sizes()
        if sizes and sizes[0] > 40:
            return int(sizes[0])
        return self._name_width

    def _on_splitter_moved(self, _pos: int = 0, _index: int = 0) -> None:
        self._name_width = self._current_name_width()
        self._sync_names_column()
        self._sync_viewport_fill()
        self._update_scale()

    def _sync_hbar_shim(self) -> None:
        h = max(12, self._bars_scroll.horizontalScrollBar().sizeHint().height())
        if self._hbar_shim.height() != h:
            self._hbar_shim.setFixedHeight(h)

    def _sync_names_column(self) -> None:
        self._sync_hbar_shim()
        width = max(80, self._names_scroll.viewport().width() or self._name_width)
        # Match bars content; both viewports share the same height thanks to shims.
        content_h = max(
            self._bars.content_height(),
            self._bars_scroll.viewport().height(),
        )
        self._names.set_rows(
            self._bars.record_ids(),
            self._bars.row_tops(),
            self._bars.row_heights(),
            self._bars.selected_row_set(),
            width=width,
            content_height=content_h,
        )

    @staticmethod
    def _safe_disconnect(signal, slot) -> None:
        """Disconnect without raising when the slot/receiver is already gone."""
        try:
            signal.disconnect(slot)
        except (RuntimeError, TypeError, SystemError):
            pass

    def _restore_zoom_bar_to_header(self, bar: AnnotationZoomBar) -> None:
        """Return AnnotationZoomBar to PreviewHeaderBar after leaving Overview."""
        host = getattr(self, "_zoom_bar_header_host", None)
        bar.set_overview_mode(False)
        bar.hide()
        self._toolbar_slot.hide()
        self._toolbar_shim.hide()
        if host is None:
            return
        layout = host.layout()
        bar.setParent(host)
        if layout is not None and layout.indexOf(bar) < 0:
            layout.insertWidget(0, bar, 0, Qt.AlignmentFlag.AlignVCenter)

    def _mount_zoom_bar_in_panel(self, bar: AnnotationZoomBar) -> None:
        """Place −/+/Position/Select/Fit above the tracks (single-seq map alignment)."""
        from repeatviewer.widgets.sequence_map_viewer import (
            _app_chrome_background,
            _apply_map_paper_background,
        )

        # Remember header host so we can put the bar back on unbind.
        if bar.parent() is not None and bar.parent() is not self._toolbar_slot:
            self._zoom_bar_header_host = bar.parent()
        layout = self._toolbar_slot_layout
        while layout.count():
            item = layout.takeAt(0)
            w = item.widget()
            if w is not None and w is not bar:
                w.setParent(None)
        bar.setParent(self._toolbar_slot)
        if layout.indexOf(bar) < 0:
            layout.addWidget(bar, 1)
        bar.set_overview_mode(True)
        bar.show()
        self._toolbar_slot.show()
        self._toolbar_shim.show()
        # Match single-seq AnnotationMapToolbar chrome (Select/Fit button fill).
        # Shim uses the same chrome so the strip reads continuous over names|tracks,
        # while −/+/Select/Fit start at the left edge of the sequence track.
        chrome = _app_chrome_background()
        _apply_map_paper_background(self._toolbar_slot, chrome)
        _apply_map_paper_background(self._toolbar_shim, chrome)
        _apply_map_paper_background(bar, chrome)

    def _unbind_header_zoom_bar(self) -> None:
        """Disconnect this panel from the shared header zoom bar (no-op if unbound)."""
        bar = self._header_zoom_bar
        if bar is None:
            return
        self._safe_disconnect(bar.zoom_slider.valueChanged, self._on_zoom_slider_changed)
        self._safe_disconnect(bar.fit_requested, self._fit_bars)
        self._safe_disconnect(bar.select_requested, self._select_by_position)
        self._safe_disconnect(bar.extract_requested, self._extract_selected)
        self._safe_disconnect(bar.filter_text_changed, self._on_name_filter_edited)
        self._safe_disconnect(bar.filter_return_pressed, self._on_name_filter_next)
        if getattr(bar, "_bound_overview", None) is self:
            bar._bound_overview = None
        self._restore_zoom_bar_to_header(bar)
        self._header_zoom_bar = None

    def bind_annotation_zoom_bar(self, bar) -> None:
        """Attach −/+/slider/Position/Select/Fit beside the rail (single-seq top)."""
        if not isinstance(bar, AnnotationZoomBar):
            return

        prev_owner = getattr(bar, "_bound_overview", None)
        if prev_owner is self and self._header_zoom_bar is bar:
            self._mount_zoom_bar_in_panel(bar)
            self._sync_header_position_spins()
            return
        if prev_owner is not None and prev_owner is not self:
            try:
                from shiboken6 import isValid

                if isValid(prev_owner):
                    prev_owner._unbind_header_zoom_bar()
                else:
                    bar._bound_overview = None
            except (RuntimeError, TypeError, AttributeError, SystemError):
                bar._bound_overview = None
        if self._header_zoom_bar is not None:
            self._unbind_header_zoom_bar()

        prev_value = self._zoom_slider.value()
        prev_filter = self._name_filter.text()

        if self._zoom_slider is not bar.zoom_slider:
            self._safe_disconnect(
                self._zoom_slider.valueChanged, self._on_zoom_slider_changed
            )
        if self._name_filter is not bar.name_filter:
            self._safe_disconnect(
                self._name_filter.textChanged, self._on_name_filter_edited
            )
            self._safe_disconnect(
                self._name_filter.returnPressed, self._on_name_filter_next
            )

        self._header_zoom_bar = bar
        bar._bound_overview = self
        self._mount_zoom_bar_in_panel(bar)
        self._zoom_slider = bar.zoom_slider
        self._name_filter = bar.name_filter
        self._updating_zoom_slider = True
        try:
            self._zoom_slider.setValue(prev_value)
        finally:
            self._updating_zoom_slider = False
        self._name_filter.blockSignals(True)
        self._name_filter.setText(prev_filter)
        self._name_filter.blockSignals(False)
        self._zoom_slider.valueChanged.connect(self._on_zoom_slider_changed)
        bar.fit_requested.connect(self._fit_bars)
        bar.select_requested.connect(self._select_by_position)
        bar.extract_requested.connect(self._extract_selected)
        bar.filter_text_changed.connect(self._on_name_filter_edited)
        bar.filter_return_pressed.connect(self._on_name_filter_next)
        self._sync_header_position_spins()

    def is_msa_mode(self) -> bool:
        return self._bars.is_msa_mode()

    def _wire_overview_annotation_rail(self) -> None:
        """Connect Annotation / Display actions for multi-seq overview."""
        rail = self._rail
        self._auto_ann_checks: dict[AutoAnnotationMode, QCheckBox] = {}
        for label, mode in _AUTO_ANNOTATION_MODES:
            check = rail.auto_ann_checks.get(label)
            if check is None:
                continue
            self._auto_ann_checks[mode] = check
            check.toggled.connect(self._on_auto_annotation_check_toggled)
        rail.auto_ann_all_check.toggled.connect(self._on_auto_annotation_all_toggled)
        rail.similarity_ann_all_check.toggled.connect(
            self._on_similarity_annotation_all_toggled
        )
        rail.featureDetectClicked.connect(self._on_auto_run_clicked)
        rail.featureDetectSettingsClicked.connect(self._open_auto_annotation_settings)
        rail.similarityAnnotateClicked.connect(self._on_similarity_run_clicked)
        rail.similaritySettingsClicked.connect(
            self._open_similarity_annotation_settings
        )
        rail.saveClicked.connect(self._save_overview_annotations)
        rail.deleteClicked.connect(self._delete_selected_feature)
        rail.restoreClicked.connect(self._restore_features)
        rail.exportClicked.connect(self._export_overview_view)
        rail.addFeatureClicked.connect(self._on_add_feature_clicked)
        rail.transferRequested.connect(self._on_transfer_annotations_clicked)
        rail.basesZoomedInChanged.connect(self._on_bases_zoomed_in_changed)
        rail.annotationDistanceChanged.connect(self._on_sequence_row_gap_changed)
        rail.showAnnotationTrackChanged.connect(self._on_show_annotation_track_changed)
        self._load_overview_display_settings()
        bases = int(rail.bases_zoomed_in.value())
        self._bars.set_msa_residue_size(bases)
        self._names.set_label_point_size(bases)
        self._bars.set_sequence_row_gap(int(rail.annotation_distance.value()))
        self._sync_auto_annotation_all_check()
        self._refresh_similarity_database_rows()

    def _overview_display_settings(self) -> QSettings:
        s = QSettings("RepeatViewer", "RepeatViewer")
        s.beginGroup(_OVERVIEW_DISPLAY_SETTINGS_GROUP)
        return s

    def _load_overview_display_settings(self) -> None:
        """MSA/overview defaults: DNA/protein labels 4, sequence distance 1."""
        rail = self._rail
        s = self._overview_display_settings()
        try:
            bases = int(s.value("bases_zoomed_in", _DEFAULT_MSA_FONT_PT))
            distance = int(
                s.value("sequence_distance", _DEFAULT_MSA_SEQUENCE_DISTANCE)
            )
        finally:
            s.endGroup()
        bases = max(2, min(36, bases))
        distance = max(0, min(40, distance))
        for widget in (rail.bases_zoomed_in, rail.annotation_distance):
            widget.blockSignals(True)
        try:
            rail.bases_zoomed_in.setValue(bases)
            rail.annotation_distance.setValue(distance)
        finally:
            for widget in (rail.bases_zoomed_in, rail.annotation_distance):
                widget.blockSignals(False)

    def _save_overview_display_settings(self) -> None:
        rail = self._rail
        s = self._overview_display_settings()
        try:
            s.setValue("bases_zoomed_in", int(rail.bases_zoomed_in.value()))
            s.setValue("sequence_distance", int(rail.annotation_distance.value()))
        finally:
            s.endGroup()

    def _on_bases_zoomed_in_changed(self, point_size: int) -> None:
        pt = int(point_size)
        self._bars.set_msa_residue_size(pt)
        self._names.set_label_point_size(pt)
        self._save_overview_display_settings()
        self._sync_names_column()

    def _on_sequence_row_gap_changed(self, gap: int) -> None:
        self._bars.set_sequence_row_gap(int(gap))
        self._save_overview_display_settings()
        self._sync_names_column()

    def _on_show_annotation_track_changed(self, enabled: bool) -> None:
        self._bars.set_show_annotation_track(bool(enabled))
        self._sync_names_column()

    def _select_all_sequences(self) -> None:
        self._bars.select_all_rows()

    def _on_region_selection_changed(self, region: object) -> None:
        bar = self._header_zoom_bar
        if bar is None or not hasattr(bar, "set_position_range"):
            return
        maximum = max(1, int(self._bars.max_length()))
        if region is None:
            # Cleared by click outside / Escape — reset Position spins to full span.
            bar.set_position_range(maximum, start=1, end=maximum)
            return
        try:
            start, end = int(region[0]), int(region[1])  # type: ignore[index]
        except (TypeError, ValueError, IndexError):
            return
        bar.set_position_range(maximum, start=start, end=end)

    def select_all_sequences(self) -> None:
        """Public hook for Edit → Select All."""
        self._select_all_sequences()

    def _on_residue_colors_changed(self) -> None:
        self._bars.refresh_msa_colors()

    def _on_colors_setting_clicked(self) -> None:
        from repeatviewer.widgets.annotation_colors_dialog import (
            AnnotationElementColorsDialog,
        )

        dialog = AnnotationElementColorsDialog(parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._on_residue_colors_changed()

    def _annotation_running(self) -> bool:
        thread = self._ann_thread
        return thread is not None and thread.isRunning()

    def _set_annotation_job_status(self, text: str) -> None:
        self._rail.set_annotation_job_status(text)
        if text:
            self.statusMessage.emit(text)

    def _clear_annotation_job_status(self) -> None:
        self._rail.set_annotation_job_status("")

    def _main_window_for_job_progress(self):
        win = self.window()
        if win is not None and hasattr(win, "_show_job_progress"):
            return win
        return None

    def _begin_window_job_progress(self, title: str, total_steps: int) -> None:
        win = self._main_window_for_job_progress()
        if win is None:
            return
        win._show_job_progress(title, max(1, int(total_steps)))

    def _update_window_job_progress(self, completed: int, detail: str = "") -> None:
        win = self._main_window_for_job_progress()
        if win is None or not hasattr(win, "_on_job_progress"):
            return
        win._on_job_progress(max(0, int(completed)), str(detail or ""))

    def _end_window_job_progress(self) -> None:
        win = self._main_window_for_job_progress()
        if win is None:
            return
        if getattr(win, "_job_running", False) and hasattr(win, "_task_progress"):
            bar = win._task_progress
            if bar is not None:
                total = max(1, int(getattr(bar, "_total", 1)))
                bar.update(total, "Done")
        if hasattr(win, "_hide_job_progress"):
            win._hide_job_progress()

    def _sync_auto_annotation_all_check(self) -> None:
        checks = list(self._auto_ann_checks.values())
        if not checks:
            return
        all_checked = all(check.isChecked() for check in checks)
        self._auto_all_updating = True
        self._rail.auto_ann_all_check.setChecked(all_checked)
        self._auto_all_updating = False

    def _on_auto_annotation_check_toggled(self, _checked: bool) -> None:
        if self._auto_all_updating:
            return
        self._sync_auto_annotation_all_check()

    def _on_auto_annotation_all_toggled(self, checked: bool) -> None:
        if self._auto_all_updating:
            return
        self._auto_all_updating = True
        try:
            for check in self._auto_ann_checks.values():
                check.setChecked(checked)
        finally:
            self._auto_all_updating = False

    def _sync_similarity_annotation_all_check(self) -> None:
        checks = [check for _entry, check in self._similarity_db_rows]
        if not checks:
            self._similarity_all_updating = True
            self._rail.similarity_ann_all_check.setChecked(False)
            self._similarity_all_updating = False
            return
        all_checked = all(check.isChecked() for check in checks)
        self._similarity_all_updating = True
        self._rail.similarity_ann_all_check.setChecked(all_checked)
        self._similarity_all_updating = False

    def _on_similarity_annotation_all_toggled(self, checked: bool) -> None:
        if self._similarity_all_updating:
            return
        self._similarity_all_updating = True
        try:
            for _entry, check in self._similarity_db_rows:
                check.setChecked(checked)
        finally:
            self._similarity_all_updating = False

    def _on_similarity_db_check_toggled(self, _checked: bool) -> None:
        if self._similarity_all_updating:
            return
        self._sync_similarity_annotation_all_check()

    def _make_similarity_db_check(self, entry, index: int) -> QCheckBox:
        color = similarity_database_color(entry.name, index)
        check = QCheckBox(entry.name, self._rail.similarity_db_host)
        readiness = database_search_readiness(entry)
        tooltip = f"{entry.display_type_label}\n{entry.path}"
        if readiness is not None:
            tooltip += f"\n\nNot ready: {readiness}"
            check.setChecked(False)
        else:
            check.setChecked(True)
        check.setToolTip(tooltip)
        check.setStyleSheet(f"QCheckBox {{ color: {color}; }}")
        check.toggled.connect(self._on_similarity_db_check_toggled)
        return check

    def _refresh_similarity_database_rows(self) -> None:
        self._rail.clear_similarity_db_checks()
        self._similarity_db_rows = []
        databases = load_similarity_annotation_databases()
        if not databases:
            empty = QLabel(
                "No similarity databases configured.", self._rail.similarity_db_host
            )
            empty.setStyleSheet("color: palette(mid);")
            empty.setWordWrap(True)
            self._rail.add_similarity_db_widget(empty)
            self._sync_similarity_annotation_all_check()
            return
        for index, entry in enumerate(databases):
            check = self._make_similarity_db_check(entry, index)
            self._rail.add_similarity_db_check(check)
            self._similarity_db_rows.append((entry, check))
        self._sync_similarity_annotation_all_check()

    def _selected_auto_annotation_modes(self) -> list[AutoAnnotationMode]:
        return [
            mode for mode, check in self._auto_ann_checks.items() if check.isChecked()
        ]

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

    def _set_annotation_controls_enabled(self, enabled: bool) -> None:
        rail = self._rail
        rail.auto_ann_settings_btn.setEnabled(enabled)
        rail.auto_ann_all_check.setEnabled(enabled)
        for check in self._auto_ann_checks.values():
            check.setEnabled(enabled)
        rail.similarity_ann_settings_btn.setEnabled(enabled)
        rail.similarity_ann_all_check.setEnabled(enabled)
        for _entry, check in self._similarity_db_rows:
            check.setEnabled(enabled)
        running = not enabled
        rail.auto_ann_run_btn.setEnabled(True)
        rail.similarity_ann_run_btn.setEnabled(True)
        rail.auto_ann_run_btn.setText("Stop" if running and self._ann_job_kind == "auto" else "Run")
        rail.similarity_ann_run_btn.setText(
            "Stop" if running and self._ann_job_kind == "similarity" else "Run"
        )

    def _open_auto_annotation_settings(self) -> None:
        from repeatviewer.widgets.sequence_map_viewer import AutoAnnotationSettingsDialog

        self._auto_ann_settings = load_auto_annotation_settings()
        dialog = AutoAnnotationSettingsDialog(self._auto_ann_settings, parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._auto_ann_settings = dialog.settings()
        save_auto_annotation_settings(self._auto_ann_settings)
        self.statusMessage.emit("Auto annotation settings saved.")

    def _open_similarity_annotation_settings(self) -> None:
        from repeatviewer.widgets.sequence_map_viewer import (
            SimilarityAnnotationSettingsDialog,
        )

        self._similarity_ann_settings = load_similarity_annotation_settings()
        dialog = SimilarityAnnotationSettingsDialog(
            self._similarity_ann_settings, parent=self
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        self._similarity_ann_settings = dialog.settings()
        save_similarity_annotation_settings(self._similarity_ann_settings)
        self.statusMessage.emit("Similarity annotation settings saved.")

    def _require_selected_record_ids(self, title: str) -> list[str] | None:
        ids = self._selected_record_ids()
        if ids:
            return ids
        QMessageBox.information(
            self,
            title,
            "Select one or more sequences first.\n\n"
            "Tip: Ctrl+A (Cmd+A on Mac) selects all sequences.",
        )
        return None

    def _on_auto_run_clicked(self) -> None:
        if self._annotation_running():
            if self._ann_job_kind == "auto":
                self._stop_annotation_job()
            return
        modes = self._selected_auto_annotation_modes()
        if not modes:
            QMessageBox.information(
                self,
                "Feature detection",
                "Select at least one detector to run.",
            )
            return
        target_ids = self._require_selected_record_ids("Feature detection")
        if not target_ids:
            return
        self._start_batch_annotation(
            TeAnnotationOptions(auto_modes=tuple(modes)),
            kind="auto",
            title="Feature detection",
            total_steps=max(1, len(target_ids)),
            record_ids=target_ids,
        )

    def _on_similarity_run_clicked(self) -> None:
        if self._annotation_running():
            if self._ann_job_kind == "similarity":
                self._stop_annotation_job()
            return
        databases = self._selected_similarity_databases()
        if not databases:
            QMessageBox.information(
                self,
                "Similarity annotation",
                "Select at least one database to search.",
            )
            return
        issues = []
        for entry in databases:
            readiness = database_search_readiness(entry.database)
            if readiness is not None:
                issues.append(f"{entry.database.name}: {readiness}")
        if issues:
            QMessageBox.warning(
                self,
                "Similarity annotation",
                "These databases are not ready to search:\n\n" + "\n\n".join(issues),
            )
            return
        target_ids = self._require_selected_record_ids("Similarity annotation")
        if not target_ids:
            return
        self._start_batch_annotation(
            TeAnnotationOptions(similarity_databases=tuple(databases)),
            kind="similarity",
            title="Similarity annotation",
            total_steps=max(1, len(target_ids)),
            record_ids=target_ids,
        )

    def _stop_annotation_job(self) -> None:
        if not self._annotation_running():
            return
        self._ann_cancel.request()
        self._set_annotation_job_status("Stopping annotation…")

    def _start_batch_annotation(
        self,
        options: TeAnnotationOptions,
        *,
        kind: str,
        title: str,
        total_steps: int,
        record_ids: list[str] | None = None,
    ) -> None:
        if self._annotation_running():
            return
        if not self._path.is_file():
            QMessageBox.warning(self, title, "No FASTA file loaded.")
            return
        target_ids = list(record_ids or self._selected_record_ids())
        if not target_ids:
            return
        self._ann_cancel.reset()
        self._ann_job_kind = kind
        self._ann_progress_total = max(1, int(total_steps))
        self._ann_last_sequence_status = ""
        self._set_annotation_controls_enabled(False)
        n = len(target_ids)
        self._set_annotation_job_status(
            f"{title}: running on {n:,} selected sequence{'s' if n != 1 else ''}…"
        )
        self.statusMessage.emit(
            f"Running {title.lower()} on {n:,} selected sequence(s) in {self._path.name}…"
        )
        self._begin_window_job_progress(title, total_steps)

        path = self._path
        auto_cfg = self._auto_ann_settings
        sim_cfg = self._similarity_ann_settings
        cancel = self._ann_cancel
        job_ids = list(target_ids)

        def _job(*, progress=None):
            def _prog(completed: int, total: int, label: str) -> None:
                if progress is None:
                    return
                if completed < 0:
                    # Nested tool chatter (e.g. "0/1 LTR") — prefix so it cannot
                    # steal the sequence-level N/M progress channel.
                    progress(f"detail:{label}")
                else:
                    progress(f"{max(0, int(completed))}/{max(1, int(total))} {label}")

            return run_te_annotation_on_fasta(
                path,
                options,
                auto_settings=auto_cfg,
                similarity_settings=sim_cfg,
                use_system_blast=load_use_system_blast_from_settings(),
                cancel=cancel,
                progress=_prog,
                record_ids=job_ids,
            )

        worker = Worker(_job, progress=True)
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.signals.progress.connect(self._on_batch_annotation_progress)
        worker.signals.finished.connect(self._on_batch_annotation_finished)
        worker.signals.error.connect(self._on_batch_annotation_error)
        worker.signals.finished.connect(thread.quit)
        worker.signals.error.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._clear_annotation_thread)
        self._ann_worker = worker
        self._ann_thread = thread
        thread.start()

    @Slot(str)
    def _on_batch_annotation_progress(self, message: str) -> None:
        text = str(message or "").strip()
        if text.startswith("detail:"):
            # Nested detector updates — keep sequence progress stable (no rail flash).
            return
        match = re.match(r"^(\d+)\s*/\s*(\d+)\s+(.+)$", text)
        if match is None:
            self._set_annotation_job_status(text)
            return
        completed, total, detail = (
            int(match.group(1)),
            int(match.group(2)),
            match.group(3).strip(),
        )
        expected = max(1, int(getattr(self, "_ann_progress_total", total) or total))
        # Ignore nested mode counters (e.g. "0/1 LTR") that are not per-sequence.
        if total != expected:
            return
        kind = "Feature detection" if self._ann_job_kind == "auto" else "Similarity"
        status = f"{kind} {completed}/{total}: {detail}"
        if status != getattr(self, "_ann_last_sequence_status", ""):
            self._ann_last_sequence_status = status
            self._set_annotation_job_status(status)
        self._update_window_job_progress(completed, detail)

    @Slot(object)
    def _on_batch_annotation_finished(self, result: object) -> None:
        self._clear_annotation_job_status()
        kind = "Feature detection" if self._ann_job_kind == "auto" else "Similarity"
        if isinstance(result, list) and result:
            self.statusMessage.emit(
                f"{kind} finished for selected sequences · saved {Path(result[0]).name}"
            )
        else:
            self.statusMessage.emit(f"{kind} finished for selected sequences.")
        was_msa = self.is_msa_mode()
        self.load_path(self._path, msa_mode=was_msa)
        self.annotationsSaved.emit(self._path)

    @Slot(str)
    def _on_batch_annotation_error(self, message: str) -> None:
        self._clear_annotation_job_status()
        lower = str(message or "").lower()
        if "stopped by user" in lower or "job stopped" in lower:
            self.statusMessage.emit("Annotation stopped.")
            return
        QMessageBox.warning(self, "Annotation failed", str(message))
        self.statusMessage.emit(f"Annotation failed: {message}")

    @Slot()
    def _clear_annotation_thread(self) -> None:
        self._ann_thread = None
        self._ann_worker = None
        self._ann_job_kind = ""
        self._ann_progress_total = 0
        self._ann_last_sequence_status = ""
        self._set_annotation_controls_enabled(True)
        self._end_window_job_progress()

    def _push_feature_undo_snapshot(self) -> None:
        self._feature_undo_stack.append(_clone_features_by_id(self._features))
        if len(self._feature_undo_stack) > _FEATURE_UNDO_LIMIT:
            self._feature_undo_stack.pop(0)

    def _refresh_features_after_edit(self) -> None:
        """Reload canvas features while preserving MSA mode, zoom, and selection."""
        keep_ppb = self._bars.pixels_per_base()
        keep_sel = self._bars.selected_rows()
        keep_region = self._bars.region_selection()
        if self.is_msa_mode():
            self._bars.set_features_in_place(self._features)
        else:
            self._bars.set_data(
                self._record_ids,
                self._lengths,
                self._features,
                fit_bars_width=None,
            )
        self._bars.set_pixels_per_base(keep_ppb)
        self._bars.set_selected_rows(keep_sel)
        if keep_region is not None:
            self._bars.set_region_selection(keep_region[0], keep_region[1])
        self._sync_names_column()
        self._update_scale()

    def _persist_overview_features(self) -> Path | None:
        try:
            _json, gff_path = save_annotations_multi(
                self._path,
                features_by_id=self._features,
                lengths_by_id=self._lengths,
                record_ids=list(self._record_ids),
            )
        except (OSError, ValueError) as exc:
            QMessageBox.critical(self, "Save failed", str(exc))
            return None
        if gff_path is not None:
            self.annotationsSaved.emit(gff_path)
        return gff_path

    def _delete_selected_feature(self) -> None:
        sel = self._bars.selected_feature()
        if sel is None:
            QMessageBox.information(
                self,
                "Delete",
                "Select a feature on a sequence track, then click Delete.",
            )
            return
        row, feat_index = sel
        if not (0 <= row < len(self._record_ids)):
            return
        rid = self._record_ids[row]
        flat: list[SequenceFeature] = []
        for track in _pack_feature_tracks(list(self._features.get(rid) or [])):
            flat.extend(track)
        if not (0 <= feat_index < len(flat)):
            return
        current = flat[feat_index]
        store = list(self._features.get(rid) or [])
        store_index = next(
            (
                i
                for i, feat in enumerate(store)
                if (
                    feat.start == current.start
                    and feat.end == current.end
                    and feat.name == current.name
                    and feat.feature_type == current.feature_type
                    and feat.strand == current.strand
                )
            ),
            -1,
        )
        if store_index < 0:
            return
        self._push_feature_undo_snapshot()
        del store[store_index]
        self._features[rid] = store
        self._bars.clear_feature_selection()
        self._refresh_features_after_edit()
        self._persist_overview_features()
        self.statusMessage.emit(f"Deleted feature on {rid}.")

    def _restore_features(self) -> None:
        if not self._feature_undo_stack:
            QMessageBox.information(
                self,
                "Restore",
                "No feature changes to restore.",
            )
            return
        self._features = self._feature_undo_stack.pop()
        self._bars.clear_feature_selection()
        self._refresh_features_after_edit()
        self._persist_overview_features()
        self.statusMessage.emit("Restored previous feature state.")

    def _export_overview_view(self) -> None:
        """Export the names + tracks overview (same formats as single-seq Export)."""
        out, _filter = QFileDialog.getSaveFileName(
            self,
            "Export overview",
            str(self._path.with_suffix(".svg")),
            "SVG (*.svg);;PNG (*.png);;PDF (*.pdf)",
        )
        if not out:
            return
        out_path = Path(out)
        widget = self._splitter
        size = widget.size()
        if size.width() <= 0 or size.height() <= 0:
            size = widget.sizeHint()
        width = max(1, int(size.width()))
        height = max(1, int(size.height()))
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
            gen.setSize(QSize(width, height))
            gen.setViewBox(QRect(0, 0, width, height))
            gen.setResolution(96)
            painter = QPainter(gen)
            try:
                widget.render(painter, QPoint(0, 0))
            finally:
                painter.end()
            self.statusMessage.emit(f"Exported overview → {out_path.name}")
            return

        img = QImage(QSize(width, height), QImage.Format.Format_ARGB32)
        img.fill(Qt.GlobalColor.transparent)
        painter = QPainter(img)
        try:
            widget.render(painter, QPoint(0, 0))
        finally:
            painter.end()

        if suffix == ".png":
            if not img.save(str(out_path), "PNG"):
                QMessageBox.warning(self, "Export", "Could not save PNG.")
                return
            self.statusMessage.emit(f"Exported overview → {out_path.name}")
            return

        if suffix == ".pdf":
            writer = QPdfWriter(str(out_path))
            writer.setPageSize(QPageSize(QPageSize.PageSizeId.A4))
            writer.setResolution(96)
            pdf_painter = QPainter(writer)
            try:
                page_rect = writer.pageRect()
                scaled = img.scaled(
                    page_rect.size().toSize(),
                    Qt.AspectRatioMode.KeepAspectRatio,
                    Qt.TransformationMode.SmoothTransformation,
                )
                x = page_rect.x() + (page_rect.width() - scaled.width()) // 2
                y = page_rect.y() + (page_rect.height() - scaled.height()) // 2
                pdf_painter.drawImage(QPoint(x, y), scaled)
            finally:
                pdf_painter.end()
            self.statusMessage.emit(f"Exported overview → {out_path.name}")
            return

        QMessageBox.warning(
            self,
            "Export",
            f"Unsupported export format: {suffix!r}. Use .svg, .png, or .pdf.",
        )

    def _save_overview_annotations(self) -> None:
        if not self._record_ids:
            return
        try:
            json_path, _gff = save_annotations_multi(
                self._path,
                features_by_id=self._features,
                lengths_by_id=self._lengths,
                record_ids=self._record_ids,
            )
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, "Save annotations", str(exc))
            return
        self.statusMessage.emit(f"Saved annotations · {json_path.name}")
        self.annotationsSaved.emit(self._path)

    def load_path(self, path: Path, *, msa_mode: bool | None = None) -> None:
        """Load multi-FASTA overview.

        When ``msa_mode`` is True (or auto-detected via equal lengths), paint
        aligned residue rows instead of length bars.
        """
        self._path = path.expanduser().resolve()
        self._bars.set_find_hits({}, None)  # hits belong to the previous file
        use_msa = bool(msa_mode) if msa_mode is not None else is_aligned_fasta(self._path)
        if use_msa:
            self._load_msa_path()
            return
        try:
            lengths = fasta_record_lengths(self._path)
        except (OSError, ValueError) as exc:
            self._record_ids = []
            self._lengths = {}
            self._features = {}
            self._bars.set_data([], {}, {})
            self._scale.set_scale(ppb=_ABS_MIN_PPB, max_length=1, bar_pad=8)
            self.statusMessage.emit(str(exc))
            return

        record_ids = list(lengths.keys())
        # Prefer JSON sidecar (canonical store); fall back to companion GFF3.
        json_by_id = load_features_by_id_from_json(self._path, record_ids)
        gff = resolve_features_gff_path(self._path)
        sidecar = resolve_sidecar_path(self._path)
        if json_by_id and any(json_by_id.values()):
            features = {rid: list(json_by_id.get(rid) or []) for rid in record_ids}
            source_label = sidecar.name if sidecar is not None else "JSON"
        elif gff is not None:
            features = load_all_features_by_id(
                gff, record_ids, lengths_by_id=lengths
            )
            if json_by_id:
                for rid, json_feats in json_by_id.items():
                    if not json_feats:
                        continue
                    current = features.get(rid) or []
                    features[rid] = json_feats if not current else current
            source_label = gff.name
        else:
            features = {rid: [] for rid in record_ids}
            source_label = None

        self._record_ids = record_ids
        self._lengths = lengths
        self._features = features
        self._did_initial_fit = False
        self._search_needle = ""
        self._search_hits = []
        self._search_hit_index = -1
        total_w = max(400, self.width() or 800)
        self._name_width = self._compute_name_width(total_w)
        tracks_w = max(200, total_w - self._name_width)
        self._splitter.setSizes([self._name_width, tracks_w])
        self._sync_viewport_fill()
        self._bars.set_data(record_ids, lengths, features, fit_bars_width=None)
        self._sync_names_column()
        self._update_scale()
        self._sync_header_position_spins()
        QTimer.singleShot(0, self._try_initial_fit)
        QTimer.singleShot(50, self._try_initial_fit)
        QTimer.singleShot(200, self._try_initial_fit)

        annotated = sum(1 for rid in record_ids if features.get(rid))
        feat_total = sum(len(features.get(rid) or []) for rid in record_ids)
        if source_label is not None and feat_total:
            self.statusMessage.emit(
                f"{len(record_ids):,} sequences · {annotated:,} with annotations · "
                f"{feat_total:,} features · {source_label}"
            )
        elif source_label is not None:
            self.statusMessage.emit(
                f"{len(record_ids):,} sequences · {source_label} (no features yet)"
            )
        else:
            self.statusMessage.emit(
                f"{len(record_ids):,} sequences · no annotation JSON / .RV.gff3 yet"
            )

    def _load_msa_path(self) -> None:
        """Load equal-length FASTA as colored residue rows with annotations."""
        try:
            record_ids, sequences = load_aligned_fasta_rows(self._path)
        except (OSError, ValueError) as exc:
            self._record_ids = []
            self._lengths = {}
            self._features = {}
            self._bars.set_data([], {}, {})
            self._scale.set_scale(ppb=_ABS_MIN_PPB, max_length=1, bar_pad=8)
            self.statusMessage.emit(str(exc))
            return
        width = len(sequences[0]) if sequences else 1
        self._record_ids = record_ids
        self._lengths = {rid: width for rid in record_ids}
        json_by_id = load_features_by_id_from_json(self._path, record_ids)
        gff = resolve_features_gff_path(self._path)
        if resolve_sidecar_path(self._path) is not None:
            self._features = {
                rid: list(json_by_id.get(rid) or []) for rid in record_ids
            }
        elif gff is not None:
            self._features = load_all_features_by_id(
                gff, record_ids, lengths_by_id=self._lengths
            )
            if json_by_id:
                for rid, json_feats in json_by_id.items():
                    if not json_feats:
                        continue
                    current = self._features.get(rid) or []
                    self._features[rid] = json_feats if not current else current
        else:
            self._features = {rid: [] for rid in record_ids}
        self._did_initial_fit = False
        self._search_needle = ""
        self._search_hits = []
        self._search_hit_index = -1
        total_w = max(400, self.width() or 800)
        self._name_width = self._compute_name_width(total_w)
        tracks_w = max(200, total_w - self._name_width)
        self._splitter.setSizes([self._name_width, tracks_w])
        self._sync_viewport_fill()
        self._bars.set_msa_data(
            record_ids,
            sequences,
            features=self._features,
            fit_bars_width=None,
        )
        self._sync_names_column()
        self._update_scale()
        self._sync_header_position_spins()
        QTimer.singleShot(0, self._try_initial_fit)
        QTimer.singleShot(50, self._try_initial_fit)
        QTimer.singleShot(200, self._try_initial_fit)
        kind = "protein" if self._bars._msa_protein else "nucleotide"
        n_feat = sum(len(v) for v in self._features.values())
        self.statusMessage.emit(
            f"{len(record_ids):,} aligned sequences · {width:,} columns · "
            f"{kind} MSA · {n_feat:,} features"
        )

    def _sync_viewport_fill(self) -> None:
        vp = self._bars_scroll.viewport()
        self._bars.set_viewport_min_size(vp.width(), vp.height())

    def _try_initial_fit(self) -> None:
        if self._did_initial_fit or not self._record_ids:
            return
        if self._bars_scroll.viewport().width() < 120:
            return
        self._fit_bars()
        self._did_initial_fit = True

    def _bars_viewport_width(self) -> int:
        return max(80, self._bars_scroll.viewport().width())

    def _sync_header_position_spins(self) -> None:
        bar = self._header_zoom_bar
        if bar is None or not hasattr(bar, "set_position_range"):
            return
        maximum = max(1, int(self._bars.max_length()))
        bar.set_position_range(maximum, start=1, end=maximum)

    def _select_by_position(self) -> None:
        """Zoom/scroll to the current region or Position spins (single-seq Select)."""
        bar = self._header_zoom_bar
        if bar is None:
            return
        maximum = max(1, int(self._bars.max_length()))
        # Prefer the blue MSA band when present — spins can lag or mismatch.
        region = self._bars.region_selection()
        if region is not None:
            start, end = int(region[0]), int(region[1])
        else:
            start = int(bar.pos_start.value())
            end = int(bar.pos_end.value())
        if end < start:
            start, end = end, start
        start = max(1, min(maximum, start))
        end = max(start, min(maximum, end))
        bar.set_position_range(maximum, start=start, end=end)
        self._bars.set_region_selection(start, end)
        # Inclusive columns → half-open window for zoom_to_bp_range.
        self._on_region_zoom(float(start - 1), float(end))
        self._sync_zoom_slider_from_ppb()
        # After zoom, scroll so the band is centered (pending scroll may be start-aligned).
        ppb = max(self._bars.pixels_per_base(), 1e-12)
        mid_x = self._bars.bar_pad() + ((start - 1 + end) / 2.0) * ppb
        vp_w = max(1, self._bars_scroll.viewport().width())
        target = int(round(mid_x - vp_w / 2.0))
        hbar = self._bars_scroll.horizontalScrollBar()
        hbar.setValue(max(0, min(hbar.maximum(), target)))
        self._update_scale()
        self._sync_names_column()

    def _fit_bars(self) -> None:
        # Keep the user's splitter ratio; only zoom tracks to the track viewport.
        self._name_width = self._current_name_width()
        self._sync_viewport_fill()
        self._bars.fit_to_bars_width(self._bars_viewport_width())
        self._bars_scroll.horizontalScrollBar().setValue(0)
        self._sync_names_column()
        self._update_scale()
        self._sync_zoom_slider_from_ppb()
        self._sync_header_position_spins()

    def _slider_ppb_bounds(self) -> tuple[float, float]:
        return self._bars.zoom_bounds(self._bars_viewport_width())

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

    def _sync_zoom_slider_from_ppb(self) -> None:
        if self._updating_zoom_slider:
            return
        self._updating_zoom_slider = True
        try:
            self._zoom_slider.setValue(
                self._ppb_to_slider_value(self._bars.pixels_per_base())
            )
        finally:
            self._updating_zoom_slider = False

    def _nudge_zoom_slider(self, delta: int) -> None:
        self._zoom_slider.setValue(
            max(0, min(100, self._zoom_slider.value() + int(delta)))
        )

    def _on_zoom_slider_changed(self, value: int) -> None:
        if self._updating_zoom_slider:
            return
        ppb = self._slider_value_to_ppb(value)
        self._sync_viewport_fill()
        self._bars.set_pixels_per_base(ppb)
        self._sync_names_column()
        self._update_scale()

    def _scroll_to_row(self, row: int) -> None:
        tops = self._bars.row_tops()
        heights = self._bars.row_heights()
        if not (0 <= row < len(tops)):
            return
        y = tops[row]
        h = heights[row]
        vbar = self._bars_scroll.verticalScrollBar()
        vp_h = max(1, self._bars_scroll.viewport().height())
        if y < vbar.value():
            vbar.setValue(y)
        elif y + h > vbar.value() + vp_h:
            vbar.setValue(max(0, y + h - vp_h))

    def _on_name_filter_edited(self, text: str) -> None:
        needle = text.strip()
        if not needle:
            self._search_needle = ""
            self._search_hits = []
            self._search_hit_index = -1
            return
        self._jump_to_name_filter(needle, advance=False)

    def _on_name_filter_next(self) -> None:
        needle = self._name_filter.text().strip()
        if not needle:
            return
        self._jump_to_name_filter(needle, advance=True)

    def _jump_to_name_filter(self, needle: str, *, advance: bool) -> None:
        lowered = needle.lower()
        hits = [
            i
            for i, rid in enumerate(self._record_ids)
            if lowered in rid.lower()
        ]
        if not hits:
            self._search_needle = needle
            self._search_hits = []
            self._search_hit_index = -1
            self.statusMessage.emit(f"No sequence name contains “{needle}”")
            return

        if (
            advance
            and needle.lower() == self._search_needle.lower()
            and hits == self._search_hits
            and self._search_hit_index >= 0
        ):
            self._search_hit_index = (self._search_hit_index + 1) % len(hits)
        else:
            self._search_needle = needle
            self._search_hits = hits
            self._search_hit_index = 0

        row = hits[self._search_hit_index]
        self._bars.set_selected_rows([row])
        self._scroll_to_row(row)
        self.statusMessage.emit(
            f"Filter “{needle}”: {self._search_hit_index + 1}/{len(hits)} — "
            f"{self._record_ids[row]}"
        )

    def _zoom_step(self, factor: float) -> None:
        vp = self._bars_scroll.viewport()
        hbar = self._bars_scroll.horizontalScrollBar()
        anchor_x = hbar.value() + vp.width() / 2
        self._bars.zoom_by(factor, anchor_x=anchor_x)

    def _update_scale(self) -> None:
        self._scale.set_scale(
            ppb=self._bars.pixels_per_base(),
            max_length=self._bars.max_length(),
            bar_pad=self._bars.bar_pad(),
        )
        self._scale.set_scroll_x(self._bars_scroll.horizontalScrollBar().value())

    def _on_hscroll(self, value: int) -> None:
        self._scale.set_scroll_x(value)

    def _on_names_vscroll(self, value: int) -> None:
        # Names bar is AlwaysOff; ignore user-driven changes (should not happen).
        if self._syncing_scroll:
            return
        self._syncing_scroll = True
        self._bars_scroll.verticalScrollBar().setValue(value)
        self._syncing_scroll = False

    def _on_bars_vscroll(self, value: int) -> None:
        if self._syncing_scroll:
            return
        self._syncing_scroll = True
        self._names_scroll.verticalScrollBar().setValue(value)
        self._syncing_scroll = False

    def _nudge_vertical_scroll(self, delta_y: int) -> None:
        if not delta_y:
            return
        bar = self._bars_scroll.verticalScrollBar()
        bar.setValue(bar.value() - int(delta_y))

    def _nudge_horizontal_scroll(self, delta_x: int) -> None:
        if not delta_x:
            return
        bar = self._bars_scroll.horizontalScrollBar()
        if bar.maximum() <= 0:
            return
        bar.setValue(bar.value() - int(delta_x))

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if event.type() != QEvent.Type.Wheel or not isinstance(event, QWheelEvent):
            return super().eventFilter(obj, event)
        names_vp = self._names_scroll.viewport()
        bars_vp = self._bars_scroll.viewport()
        if obj not in (names_vp, bars_vp):
            return super().eventFilter(obj, event)
        mods = event.modifiers() | QGuiApplication.keyboardModifiers()
        if mods & Qt.KeyboardModifier.ShiftModifier:
            # Handle zoom here — wheel is delivered to the viewport, not the canvas.
            delta = event.angleDelta().y()
            if delta == 0:
                delta = event.angleDelta().x()
            pixel = event.pixelDelta()
            if delta == 0 and not pixel.isNull():
                delta = int(pixel.y()) or int(pixel.x())
            if delta != 0:
                factor = _ZOOM_STEP if delta > 0 else 1.0 / _ZOOM_STEP
                # Map viewport x → content x (include horizontal scroll).
                anchor_x = float(event.position().x()) + float(
                    self._bars_scroll.horizontalScrollBar().value()
                )
                self._bars.zoom_by(factor, anchor_x=anchor_x)
            event.accept()
            return True
        # Trackpad: apply both axes. Prefer the dominant axis when both move.
        dx = _wheel_delta_x(
            event, scroll_step=self._bars_scroll.horizontalScrollBar().singleStep() * 3
        )
        dy = _wheel_delta_y(event)
        if dx and dy and abs(dx) >= abs(dy):
            dy = 0
        elif dx and dy and abs(dy) > abs(dx):
            dx = 0
        handled = False
        if dx:
            self._nudge_horizontal_scroll(dx)
            handled = True
        if dy:
            self._nudge_vertical_scroll(dy)
            handled = True
        if handled:
            event.accept()
            return True
        return False

    def _on_zoom_changed(self, _ppb: float) -> None:
        self._sync_viewport_fill()
        self._sync_names_column()
        hbar = self._bars_scroll.horizontalScrollBar()
        scroll_to = self._bars.take_pending_scroll_bp()
        anchor = self._bars.take_pending_anchor()
        if scroll_to is not None:
            target = int(self._bars.bp_to_x(float(scroll_to)) - 4)
            hbar.setValue(max(0, target))
        elif anchor is not None:
            anchor_bp, anchor_x = anchor
            new_x = self._bars.bp_to_x(float(anchor_bp))
            view_x = float(anchor_x) - hbar.value()
            hbar.setValue(max(0, int(new_x - view_x)))
        self._update_scale()
        self._sync_zoom_slider_from_ppb()

    # ------------------------------------------------ Edit > Find / Go To

    def record_ids(self) -> list[str]:
        return list(self._record_ids)

    def record_length(self, seq_id: str) -> int:
        return int(self._lengths.get(seq_id, 0))

    def selected_record_id(self) -> str | None:
        ids = self._selected_record_ids()
        return ids[0] if ids else None

    def features_by_id(self) -> dict[str, list[SequenceFeature]]:
        return {rid: list(feats) for rid, feats in self._features.items()}

    def _row_for_id(self, seq_id: str) -> int:
        wanted = str(seq_id or "").strip()
        if wanted in self._record_ids:
            return self._record_ids.index(wanted)
        from repeatviewer.backends.gff_features import _seqid_matches

        for row, rid in enumerate(self._record_ids):
            if _seqid_matches(rid, wanted):
                return row
        return -1

    def show_find_hit(
        self,
        seq_id: str,
        start: int,
        end: int,
        hits_by_id: dict[str, list[tuple[int, int]]],
        *,
        whole_sequence: bool = False,
    ) -> bool:
        """Select the record, scroll to it and zoom to [start, end]. False if unknown."""
        row = self._row_for_id(seq_id)
        if row < 0:
            return False
        rid = self._record_ids[row]
        self._bars.set_find_hits(hits_by_id, (rid, int(start), int(end)))
        self._bars.set_selected_rows([row])
        self._sync_names_column()
        self._scroll_to_row(row)
        if whole_sequence:
            self._on_region_zoom(0.0, float(max(1, self._lengths.get(rid, end))))
        else:
            span = max(1, int(end) - int(start) + 1)
            context = max(3 * span, 60)
            center = (int(start) + int(end)) / 2.0
            self._on_region_zoom(max(0.0, center - context / 2.0), center + context / 2.0)
        return True

    def clear_find_hits(self) -> None:
        self._bars.set_find_hits({}, None)

    def add_features_and_save(self, new_by_id: dict[str, list[SequenceFeature]]) -> bool:
        """Append features per record and save the annotations (JSON sidecar)."""
        added = {
            rid: feats
            for rid, feats in new_by_id.items()
            if feats and rid in self._lengths
        }
        if not added:
            return False
        self._push_feature_undo_snapshot()
        for rid, feats in added.items():
            self._features[rid] = list(self._features.get(rid) or []) + list(feats)
        try:
            json_path, _gff = save_annotations_multi(
                self._path,
                features_by_id=self._features,
                lengths_by_id=self._lengths,
                record_ids=list(self._record_ids),
            )
        except OSError as exc:
            QMessageBox.critical(self, "Save failed", str(exc))
            return False
        keep_ppb = self._bars.pixels_per_base()
        keep_sel = self._bars.selected_rows()
        keep_region = self._bars.region_selection()
        was_msa = self._bars.is_msa_mode()
        if was_msa:
            self._bars.set_features_in_place(self._features)
        else:
            self._bars.set_data(
                self._record_ids, self._lengths, self._features, fit_bars_width=None
            )
            self._bars.set_pixels_per_base(keep_ppb)
        self._bars.set_selected_rows(keep_sel)
        if keep_region is not None:
            self._bars.set_region_selection(keep_region[0], keep_region[1])
        self._sync_names_column()
        self._update_scale()
        self.annotationsSaved.emit(json_path)
        return True

    def _current_region_for_add(self) -> tuple[int, int] | None:
        """Region from canvas selection or header Position spins."""
        region = self._bars.region_selection()
        if region is not None:
            return region
        bar = self._header_zoom_bar
        if bar is None:
            return None
        start = int(bar.pos_start.value())
        end = int(bar.pos_end.value())
        if end < start:
            start, end = end, start
        return start, end

    def _on_add_feature_clicked(self) -> None:
        from repeatviewer.widgets.sequence_map_viewer import AddFeatureDialog

        target_ids = self._require_selected_record_ids("Add feature")
        if not target_ids:
            return
        region = self._current_region_for_add()
        if region is None:
            QMessageBox.information(
                self,
                "Add feature",
                "Drag on the tracks to select a region (or set Position), "
                "select one or more sequences, then Add.",
            )
            return
        start, end = region
        label, feature_type = self._rail.selected_manual_feature()
        existing = 0
        for rid in target_ids:
            existing += sum(
                1
                for feat in (self._features.get(rid) or [])
                if feat.feature_type == feature_type
            )
        if feature_type == "ORF":
            default_name = f"ORF_{existing + 1}"
        else:
            default_name = f"{(feature_type or label or 'feature').replace(' ', '_')}{existing + 1}"
        seq_len = max(
            (int(self._lengths.get(rid, 1)) for rid in target_ids),
            default=max(1, int(self._bars.max_length())),
        )
        dialog = AddFeatureDialog(
            start=start,
            end=end,
            seq_length=seq_len,
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
        new_by_id: dict[str, list[SequenceFeature]] = {}
        for rid in target_ids:
            rid_len = max(1, int(self._lengths.get(rid, 1)))
            a = max(1, min(rid_len, int(feat.start)))
            b = max(1, min(rid_len, int(feat.end)))
            if b < a:
                a, b = b, a
            new_by_id[rid] = [
                SequenceFeature(
                    name=feat.name,
                    feature_type=feat.feature_type,
                    start=a,
                    end=b,
                    strand=feat.strand,
                    color=feat.color,
                    track=feat.track,
                    source=feat.source,
                    parent=feat.parent,
                )
            ]
        if self.add_features_and_save(new_by_id):
            n = len(target_ids)
            self.statusMessage.emit(
                f"Added “{feat.name}” to {n:,} sequence{'s' if n != 1 else ''}."
            )

    def _sequence_text_for_id(self, seq_id: str) -> str:
        """Return sequence string for a record (MSA cache or disk)."""
        if self._bars.is_msa_mode():
            try:
                row = self._record_ids.index(seq_id)
            except ValueError:
                row = -1
            if 0 <= row < len(self._bars._msa_seqs):
                return str(self._bars._msa_seqs[row])
        from Bio import SeqIO

        for rec in SeqIO.parse(str(self._path), "fasta"):
            if str(rec.id) == seq_id or str(rec.id).split("#", 1)[0] == seq_id.split(
                "#", 1
            )[0]:
                return str(rec.seq)
        raise ValueError(f"Sequence “{seq_id}” was not found in {self._path.name}.")

    def _on_transfer_annotations_clicked(self) -> None:
        if self._transfer_thread is not None and self._transfer_thread.isRunning():
            self._transfer_cancel.request()
            self._rail.transfer_run_btn.setText("Stopping…")
            return
        try:
            self._start_overview_transfer_worker()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Transfer", str(exc))

    def _start_overview_transfer_worker(self) -> None:
        from repeatviewer.backends.annotation_transfer import (
            AnnotationTransferOptions,
            discover_annotated_transfer_sources,
            load_use_system_mafft_from_settings,
            transfer_annotations,
        )
        from repeatviewer.backends.similarity_annotation import (
            load_use_system_blast_from_settings,
        )

        target_ids = self._require_selected_record_ids("Transfer")
        if not target_ids:
            return
        source_path = None
        if hasattr(self._rail, "transfer_source_row"):
            source_path = self._rail.transfer_source_row.path()
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

        source_label = source_path.name
        if source_path.is_dir():
            source_label = f"{source_path.name}/ ({len(sources)} annotated)"
        ann_source = f"Transfer ({source_label})"
        opts = AnnotationTransferOptions(
            min_identity=float(self._rail.transfer_min_identity_value()),
            min_coverage=float(self._rail.transfer_min_coverage_value()),
            allow_reverse=bool(self._rail.transfer_allow_reverse.isChecked()),
            annotation_source=ann_source,
            use_system_blast=load_use_system_blast_from_settings(),
            use_system_mafft=load_use_system_mafft_from_settings(),
        )

        self._transfer_cancel = JobCancellation()
        self._rail.transfer_run_btn.setText("Stop")
        self._rail.transfer_run_btn.setFixedWidth(52)
        n_targets = len(target_ids)
        self.statusMessage.emit(
            f"Transferring annotations onto {n_targets:,} selected sequence(s)…"
        )
        self._begin_window_job_progress(
            "Transfer annotation", max(1, n_targets * max(1, len(sources)))
        )

        job_sources = list(sources)
        job_targets = list(target_ids)
        job_path = self._path
        job_opts = opts
        cancel = self._transfer_cancel
        # Preload target sequences on the UI thread (worker must not touch widgets).
        target_seqs: dict[str, str] = {}
        for rid in job_targets:
            try:
                target_seqs[rid] = self._sequence_text_for_id(rid)
            except ValueError as exc:
                QMessageBox.warning(self, "Transfer", str(exc))
                return

        def _job(*, progress=None):
            from Bio import SeqIO

            source_seqs: dict[tuple[str, str], tuple[str, list]] = {}
            for entry in job_sources:
                key = (str(entry.fasta_path), str(entry.sequence_id))
                if key in source_seqs:
                    continue
                feats = []
                seq = ""
                try:
                    from repeatviewer.backends.sequence_features import (
                        load_features_by_id_from_json as _load_json,
                    )

                    by_id = _load_json(entry.fasta_path, [entry.sequence_id])
                    feats = list(by_id.get(entry.sequence_id) or [])
                    for rec in SeqIO.parse(str(entry.fasta_path), "fasta"):
                        if str(rec.id) == entry.sequence_id:
                            seq = str(rec.seq)
                            break
                except Exception:
                    continue
                source_seqs[key] = (seq, feats)

            new_by_id: dict[str, list] = {rid: [] for rid in job_targets}
            step = 0
            total = max(1, len(job_targets) * max(1, len(job_sources)))
            for tgt_id in job_targets:
                if cancel.requested:
                    from repeatviewer.backends.job_cancel import JobCancelledError

                    raise JobCancelledError("Annotation transfer cancelled.")
                tgt_seq = target_seqs.get(tgt_id) or ""
                if not tgt_seq:
                    continue
                for entry in job_sources:
                    step += 1
                    if progress is not None:
                        progress(f"{step}/{total} {tgt_id} ← {entry.sequence_id}")
                    key = (str(entry.fasta_path), str(entry.sequence_id))
                    src_seq, src_feats = source_seqs.get(key, ("", []))
                    if not src_seq or not src_feats:
                        continue
                    try:
                        same_file = (
                            entry.fasta_path.resolve() == job_path.resolve()
                            and entry.sequence_id == tgt_id
                        )
                    except OSError:
                        same_file = False
                    if same_file:
                        continue
                    # MSA rows: strip gaps for transfer, project hits back to columns.
                    from repeatviewer.backends.msa_coordinates import (
                        project_features_to_gapped,
                        project_features_to_ungapped,
                        ungap_sequence,
                    )

                    src_u, src_keep = ungap_sequence(src_seq)
                    tgt_u, tgt_keep = ungap_sequence(tgt_seq)
                    src_feats_u = project_features_to_ungapped(
                        list(src_feats), src_keep
                    )
                    if not src_u or not tgt_u or not src_feats_u:
                        continue
                    result = transfer_annotations(
                        source_id=entry.sequence_id,
                        source_sequence=src_u,
                        target_id=tgt_id,
                        target_sequence=tgt_u,
                        features=list(src_feats_u),
                        opts=job_opts,
                        cancel=cancel,
                        source_fasta=entry.fasta_path,
                    )
                    for item in getattr(result, "lifted", []) or []:
                        feat = getattr(item, "feature", None)
                        if feat is None:
                            continue
                        projected = project_features_to_gapped([feat], tgt_keep)
                        new_by_id[tgt_id].extend(projected)
            return new_by_id

        worker = Worker(_job, progress=True)
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.signals.progress.connect(self._on_transfer_progress)
        worker.signals.finished.connect(self._on_transfer_finished)
        worker.signals.error.connect(self._on_transfer_error)
        worker.signals.finished.connect(thread.quit)
        worker.signals.error.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._clear_transfer_thread)
        self._transfer_worker = worker
        self._transfer_thread = thread
        thread.start()

    @Slot(str)
    def _on_transfer_progress(self, message: str) -> None:
        text = str(message or "").strip()
        match = re.match(r"^(\d+)\s*/\s*(\d+)\s+(.+)$", text)
        if match is None:
            self.statusMessage.emit(text)
            return
        completed, _total, detail = (
            int(match.group(1)),
            int(match.group(2)),
            match.group(3).strip(),
        )
        self._update_window_job_progress(completed, detail)
        self.statusMessage.emit(f"Transfer {completed}: {detail}")

    @Slot(object)
    def _on_transfer_finished(self, result: object) -> None:
        self._rail.transfer_run_btn.setText("Run")
        self._rail.transfer_run_btn.setFixedWidth(52)
        self._end_window_job_progress()
        if not isinstance(result, dict):
            return
        new_by_id = {
            rid: list(feats)
            for rid, feats in result.items()
            if feats and rid in self._lengths
        }
        if not new_by_id:
            self.statusMessage.emit("Transfer finished — no features placed.")
            return
        if self.add_features_and_save(new_by_id):
            n = sum(len(v) for v in new_by_id.values())
            self.statusMessage.emit(
                f"Transfer saved {n:,} feature(s) on {len(new_by_id):,} sequence(s)."
            )

    @Slot(str)
    def _on_transfer_error(self, message: str) -> None:
        self._rail.transfer_run_btn.setText("Run")
        self._rail.transfer_run_btn.setFixedWidth(52)
        self._end_window_job_progress()
        lower = str(message or "").lower()
        if "cancelled" in lower or "stopped" in lower:
            self.statusMessage.emit("Transfer stopped.")
            return
        QMessageBox.warning(self, "Transfer failed", str(message))
        self.statusMessage.emit(f"Transfer failed: {message}")

    @Slot()
    def _clear_transfer_thread(self) -> None:
        self._transfer_thread = None
        self._transfer_worker = None
        self._rail.transfer_run_btn.setText("Run")
        self._rail.transfer_run_btn.setFixedWidth(52)

    def _on_region_zoom(self, start_bp: float, end_bp: float) -> None:
        self._sync_viewport_fill()
        self._bars.zoom_to_bp_range(start_bp, end_bp, self._bars_viewport_width())

    def _on_name_row_clicked(self, row: int, modifiers: object) -> None:
        self._bars._select_row_click(row, modifiers)  # type: ignore[arg-type]

    def _on_name_rows_dragged(self, rows: object) -> None:
        row_list = [int(r) for r in (rows or [])]
        self._bars.set_selected_rows(row_list)

    def _on_feature_activated(self, row: int, feat_index: int) -> None:
        self._modify_feature_at(row, feat_index)

    def _selected_record_ids(self) -> list[str]:
        ids: list[str] = []
        for row in self._bars.selected_rows():
            if 0 <= row < len(self._record_ids):
                ids.append(self._record_ids[row])
        return ids

    def _extract_selected(self) -> None:
        from repeatviewer.backends.sequence_compose import (
            SequenceComposeError,
            extract_records,
            proposed_extract_path,
        )
        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        ids = self._selected_record_ids()
        if not ids:
            QMessageBox.information(
                self,
                "Extract",
                "Select one or more sequences, then Extract.",
            )
            return

        proposed = proposed_extract_path(
            self._path, ids, self._path.parent, single_file=True
        )
        confirmed = confirm_existing_file_path(
            self,
            proposed,
            title="Extract",
            informative=(
                "Overwrite the existing file, or save under a different name?"
            ),
        )
        if confirmed is None:
            return

        try:
            written = extract_records(
                self._path,
                set(ids),
                self._path.parent,
                single_file=True,
                output_path=confirmed,
            )
        except SequenceComposeError as exc:
            QMessageBox.warning(self, "Extract", str(exc))
            return

        for path in written:
            self.sequenceExtracted.emit(path)
        names = ", ".join(p.name for p in written[:4])
        if len(written) > 4:
            names += ", …"
        self.statusMessage.emit(f"Extracted {len(written)} file(s): {names}")

    def _modify_feature_at(self, row: int, feat_index: int) -> None:
        from repeatviewer.widgets.sequence_map_viewer import ModifyFeatureDialog

        if not (0 <= row < len(self._record_ids)):
            return
        rid = self._record_ids[row]
        flat: list[SequenceFeature] = []
        for track in _pack_feature_tracks(list(self._features.get(rid) or [])):
            flat.extend(track)
        if not (0 <= feat_index < len(flat)):
            return
        current = flat[feat_index]
        store = list(self._features.get(rid) or [])
        store_index = next(
            (
                i
                for i, feat in enumerate(store)
                if (
                    feat.start == current.start
                    and feat.end == current.end
                    and feat.name == current.name
                    and feat.feature_type == current.feature_type
                    and feat.strand == current.strand
                )
            ),
            -1,
        )
        if store_index < 0:
            return
        seq_length = max(1, int(self._lengths.get(rid, 1)))
        parent_candidates = [
            feat.gff_id()
            for i, feat in enumerate(store)
            if i != store_index and feat.gff_id()
        ]
        dialog = ModifyFeatureDialog(
            feature=store[store_index],
            seq_length=seq_length,
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
        old_id = store[store_index].gff_id()
        store[store_index] = edited
        store = remount_parent_ids(store, old_id=old_id, new_id=edited.gff_id())
        self._features[rid] = store
        self._refresh_features_after_edit()
        gff_path = self._persist_overview_features()
        if gff_path is not None:
            self.statusMessage.emit(f"Updated feature on {rid} → {gff_path.name}")
        else:
            self.statusMessage.emit(f"Updated feature on {rid}.")

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        QTimer.singleShot(0, self._try_initial_fit)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._sync_hbar_shim()
        self._sync_viewport_fill()
        if self._record_ids and not self._did_initial_fit:
            self._try_initial_fit()
        else:
            self._sync_names_column()
            self._update_scale()
