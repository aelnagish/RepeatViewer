"""Multi-sequence FASTA annotation overview (names | tracks, synced rows)."""

from __future__ import annotations

import math
from pathlib import Path

from PySide6.QtCore import QEvent, QPointF, QRect, QRectF, QSize, Qt, QTimer, Signal
from PySide6.QtGui import (
    QColor,
    QFont,
    QFontMetrics,
    QGuiApplication,
    QPainter,
    QPalette,
    QPen,
    QPolygonF,
    QWheelEvent,
)
from PySide6.QtWidgets import (
    QDialog,
    QLineEdit,
    QMessageBox,
    QScrollArea,
    QSizePolicy,
    QSlider,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.sequence_features import (
    SequenceFeature,
    fasta_record_lengths,
    load_all_features_by_id,
    remount_parent_ids,
    resolve_features_gff_path,
    save_annotations_multi,
)
from repeatviewer.widgets.preview_controls import (
    AnnotationZoomBar,
)

_BAR_HEIGHT = 7
_BAR_TOP_PAD = 3
_FEAT_HEIGHT = 10
_FEAT_GAP = 2
_ROW_BOTTOM_PAD = 4
_MIN_ROW_HEIGHT = 18
_NAME_COLUMN_DEFAULT = 280
_NAME_PAD = 8
_SELECTED_BG_LIGHT = QColor(11, 92, 171, 36)
_SELECTED_BG_DARK = QColor(80, 160, 255, 48)
_DEFAULT_FEAT = QColor("#c0392b")
_MIN_PPB = 0.002
_MAX_PPB = 4.0
_ZOOM_STEP = 1.08
_SCALE_HEIGHT = 28
# Ignore tiny moves (click / accidental jitter); macOS-like drag threshold.
_RUBBER_MIN_PX = 14


def _is_dark(widget: QWidget) -> bool:
    return widget.palette().color(QPalette.ColorRole.Window).lightness() < 128


def _wheel_delta_y(event: QWheelEvent) -> int:
    """Pixels to scroll for a wheel/trackpad event (smooth on macOS trackpads)."""
    pixel = event.pixelDelta()
    if not pixel.isNull() and pixel.y() != 0:
        return int(pixel.y())
    angle = event.angleDelta()
    dy = int(angle.y()) if angle.y() != 0 else int(angle.x())
    if dy == 0:
        return 0
    # Mouse wheel notches are typically ±120; map to ~one row, not 120px jumps.
    return int(round(dy * (_MIN_ROW_HEIGHT / 120.0)))


def _bar_color(widget: QWidget, *, selected: bool) -> QColor:
    if selected:
        return QColor("#6eb0ef") if _is_dark(widget) else QColor("#5b8fc7")
    # Soft slate — keep bars readable without dominating the features.
    return QColor("#9aa3ad") if _is_dark(widget) else QColor("#b0b6be")


def _pack_feature_tracks(
    features: list[SequenceFeature],
) -> list[list[SequenceFeature]]:
    tracks: list[list[SequenceFeature]] = []
    for feat in sorted(features, key=lambda f: (f.start, f.end, f.name)):
        placed = False
        for track in tracks:
            if all(feat.start > other.end or feat.end < other.start for other in track):
                track.append(feat)
                placed = True
                break
        if not placed:
            tracks.append([feat])
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

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._record_ids: list[str] = []
        self._row_tops: list[int] = []
        self._row_heights: list[int] = []
        self._selected_rows: set[int] = set()
        self._drag_origin: QPointF | None = None
        self._dragging = False
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent, True)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Minimum)

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

    def paintEvent(self, event) -> None:  # noqa: N802
        painter = QPainter(self)
        clip = event.rect()
        base = self.palette().color(QPalette.ColorRole.Base)
        painter.fillRect(clip, base)
        sel_bg = _SELECTED_BG_DARK if _is_dark(self) else _SELECTED_BG_LIGHT
        text_color = self.palette().color(QPalette.ColorRole.Text)
        mid = self.palette().color(QPalette.ColorRole.Mid)
        font = QFont(self.font())
        font.setPointSize(max(10, font.pointSize() - 1))
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
            bar_y = top + _BAR_TOP_PAD
            text_y = int(bar_y + (_BAR_HEIGHT - text_h) / 2)
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
    """Length bars + feature arrows (no name column)."""

    selectionChanged = Signal()
    featureActivated = Signal(int, int)
    zoomChanged = Signal(float)
    regionZoomRequested = Signal(float, float)
    rowsChanged = Signal()  # layout / selection updated — refresh names column

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._record_ids: list[str] = []
        self._lengths: dict[str, int] = {}
        self._features: dict[str, list[SequenceFeature]] = {}
        self._tracks: dict[str, list[list[SequenceFeature]]] = {}
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
        self._pending_scroll_bp: float | None = None
        self._pending_anchor: tuple[float, float] | None = None
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setSizePolicy(QSizePolicy.Policy.MinimumExpanding, QSizePolicy.Policy.Minimum)
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent, True)

    def set_data(
        self,
        record_ids: list[str],
        lengths: dict[str, int],
        features: dict[str, list[SequenceFeature]],
        *,
        fit_bars_width: int | None = None,
    ) -> None:
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
        self._rubber = QRect()
        self._rubber_active = False
        self._drag_origin = None
        if fit_bars_width and fit_bars_width > 40 and self._max_length > 0:
            self._ppb = max(_MIN_PPB, (fit_bars_width - 16) / self._max_length)
        self._recompute_layout()
        self.update()
        self.zoomChanged.emit(self._ppb)
        self.rowsChanged.emit()

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

    def set_viewport_min_size(self, width: int, height: int) -> None:
        self._min_viewport = QSize(max(40, int(width)), max(40, int(height)))
        self._recompute_layout()

    def set_pixels_per_base(self, ppb: float) -> None:
        self._ppb = max(_MIN_PPB, min(_MAX_PPB, float(ppb)))
        self._recompute_layout()
        self.update()
        self.zoomChanged.emit(self._ppb)

    def zoom_by(self, factor: float, *, anchor_x: float | None = None) -> None:
        old = self._ppb
        new = max(_MIN_PPB, min(_MAX_PPB, old * factor))
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
        self._recompute_layout()
        self.update()
        self.zoomChanged.emit(self._ppb)

    def fit_to_bars_width(self, bars_width: int) -> None:
        if self._max_length <= 0:
            return
        self._pending_scroll_bp = None
        self._pending_anchor = None
        self.set_pixels_per_base(
            max(_MIN_PPB, (max(40, bars_width) - 16) / self._max_length)
        )

    def zoom_to_bp_range(
        self, start_bp: float, end_bp: float, bars_width: int
    ) -> None:
        span = max(1.0, float(end_bp) - float(start_bp))
        width = max(40, int(bars_width))
        self._pending_scroll_bp = float(start_bp)
        self._pending_anchor = None
        self.set_pixels_per_base(max(_MIN_PPB, (width - 16) / span))

    def take_pending_scroll_bp(self) -> float | None:
        value = self._pending_scroll_bp
        self._pending_scroll_bp = None
        return value

    def take_pending_anchor(self) -> tuple[float, float] | None:
        value = self._pending_anchor
        self._pending_anchor = None
        return value

    def set_selected_rows(self, rows: list[int] | set[int]) -> None:
        self._selected_rows = {int(r) for r in rows if 0 <= int(r) < len(self._record_ids)}
        self._selected_feat = None
        self.update()
        self.selectionChanged.emit()
        self.rowsChanged.emit()

    def _recompute_layout(self) -> None:
        self._row_tops = []
        self._row_heights = []
        y = 0
        for rid in self._record_ids:
            n_tracks = len(self._tracks.get(rid) or [])
            h = _row_height(n_tracks)
            self._row_tops.append(y)
            self._row_heights.append(h)
            y += h
        content_w = int(self._bar_pad + self._max_length * self._ppb + 12)
        # Extra bottom pad so the last row clears the horizontal scrollbar edge.
        y = y + 4 if y else _MIN_ROW_HEIGHT
        w = max(120, content_w, self._min_viewport.width())
        h = max(_MIN_ROW_HEIGHT, y, self._min_viewport.height())
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

    def _hit_feature(self, row: int, x: float, y: float) -> int | None:
        if not (0 <= row < len(self._record_ids)):
            return None
        rid = self._record_ids[row]
        length = max(1, int(self._lengths.get(rid, 1)))
        top = self._row_tops[row]
        bar_y = top + _BAR_TOP_PAD
        feat_y0 = bar_y + _BAR_HEIGHT + 2
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
            # Only activate region-zoom UI after a real drag distance.
            if rubber.width() >= _RUBBER_MIN_PX or rubber.height() >= _RUBBER_MIN_PX:
                self._rubber = rubber
                self._rubber_active = True
                rows = self._rows_in_y_range(rubber.top(), rubber.bottom())
                if rows:
                    self._selected_rows = set(rows)
                    self._selected_feat = None
                    self.selectionChanged.emit()
                    self.rowsChanged.emit()
                self.update()
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
            if active and (
                rubber.width() >= _RUBBER_MIN_PX or rubber.height() >= _RUBBER_MIN_PX
            ):
                rows = self._rows_in_y_range(rubber.top(), rubber.bottom())
                if rows:
                    self._selected_rows = set(rows)
                    self._selected_feat = None
                    self.selectionChanged.emit()
                    self.rowsChanged.emit()
                start_bp = self.x_to_bp(rubber.left())
                end_bp = self.x_to_bp(rubber.right())
                if end_bp < start_bp:
                    start_bp, end_bp = end_bp, start_bp
                if end_bp - start_bp < 1:
                    end_bp = start_bp + 1
                self.regionZoomRequested.emit(start_bp, end_bp)
            self.update()
            event.accept()
            return
        super().mouseReleaseEvent(event)

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
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        clip = event.rect()
        base = self.palette().color(QPalette.ColorRole.Base)
        painter.fillRect(clip, base)
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

            feat_y0 = bar_y + _BAR_HEIGHT + 2
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
                    if self._selected_feat == (index, flat_index):
                        painter.setBrush(Qt.BrushStyle.NoBrush)
                        painter.setPen(
                            QPen(QColor("#4ea1ff" if _is_dark(self) else "#0b5cab"), 2)
                        )
                        painter.drawRect(rect.adjusted(-1, -1, 1, 1))
                    flat_index += 1

        # Region-zoom overlay only after a real drag (not on click).
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
        self._ppb = max(_MIN_PPB, float(ppb))
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
        self._names.rowClicked.connect(self._on_name_row_clicked)
        self._names.rowsDragSelected.connect(self._on_name_rows_dragged)
        # Bars owns the vertical scrollbar; names stays synced. Wheel over either
        # pane uses pixelDelta (trackpad) or a gentle mouse-wheel step.
        self._names_scroll.viewport().installEventFilter(self)
        self._bars_scroll.viewport().installEventFilter(self)
        self._bars_scroll.verticalScrollBar().setSingleStep(18)
        self._bars_scroll.verticalScrollBar().setPageStep(120)

        # Zoom / Fit / Extract / filter live on AnnotationZoomBar (same row as
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

        names_col = QWidget(self)
        names_layout = QVBoxLayout(names_col)
        names_layout.setContentsMargins(0, 0, 0, 0)
        names_layout.setSpacing(0)
        names_layout.addWidget(self._names_scroll, 1)
        names_layout.addWidget(self._hbar_shim)
        names_layout.addWidget(self._scale_shim)

        tracks_col = QWidget(self)
        tracks_layout = QVBoxLayout(tracks_col)
        tracks_layout.setContentsMargins(0, 0, 0, 0)
        tracks_layout.setSpacing(0)
        tracks_layout.addWidget(self._bars_scroll, 1)
        tracks_layout.addWidget(self._scale)

        self._splitter = QSplitter(Qt.Orientation.Horizontal, self)
        self._splitter.addWidget(names_col)
        self._splitter.addWidget(tracks_col)
        self._splitter.setStretchFactor(0, 0)
        self._splitter.setStretchFactor(1, 1)
        self._splitter.setChildrenCollapsible(False)
        self._splitter.splitterMoved.connect(self._on_splitter_moved)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._splitter, 1)

    def _compute_name_width(self, viewport_width: int | None = None) -> int:
        font = QFont(self.font())
        font.setPointSize(max(10, font.pointSize() - 1))
        fm = QFontMetrics(font)
        widest = 0
        for rid in self._record_ids[:500]:
            widest = max(widest, fm.horizontalAdvance(rid))
        natural = widest + 2 * _NAME_PAD + 4
        cap = _NAME_COLUMN_DEFAULT + 40
        if viewport_width and viewport_width > 100:
            cap = min(cap, max(140, int(viewport_width * 0.38)))
        return max(140, min(cap, natural))

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

    def _unbind_header_zoom_bar(self) -> None:
        """Disconnect this panel from the shared header zoom bar (no-op if unbound)."""
        bar = self._header_zoom_bar
        if bar is None:
            return
        self._safe_disconnect(bar.zoom_slider.valueChanged, self._on_zoom_slider_changed)
        self._safe_disconnect(bar.fit_requested, self._fit_bars)
        self._safe_disconnect(bar.extract_requested, self._extract_selected)
        self._safe_disconnect(bar.filter_text_changed, self._on_name_filter_edited)
        self._safe_disconnect(bar.filter_return_pressed, self._on_name_filter_next)
        if getattr(bar, "_bound_overview", None) is self:
            bar._bound_overview = None
        self._header_zoom_bar = None

    def bind_annotation_zoom_bar(self, bar) -> None:
        """Attach header −/+/slider/Fit/Extract/filter (same row as Text/Annotation)."""
        if not isinstance(bar, AnnotationZoomBar):
            return

        prev_owner = getattr(bar, "_bound_overview", None)
        if prev_owner is self and self._header_zoom_bar is bar:
            bar.set_overview_mode(True)
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
        bar.set_overview_mode(True)
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
        bar.extract_requested.connect(self._extract_selected)
        bar.filter_text_changed.connect(self._on_name_filter_edited)
        bar.filter_return_pressed.connect(self._on_name_filter_next)

    def load_path(self, path: Path) -> None:
        self._path = path.expanduser().resolve()
        try:
            lengths = fasta_record_lengths(self._path)
        except (OSError, ValueError) as exc:
            self._record_ids = []
            self._lengths = {}
            self._features = {}
            self._bars.set_data([], {}, {})
            self._scale.set_scale(ppb=_MIN_PPB, max_length=1, bar_pad=8)
            self.statusMessage.emit(str(exc))
            return

        record_ids = list(lengths.keys())
        gff = resolve_features_gff_path(self._path)
        features: dict[str, list[SequenceFeature]] = {rid: [] for rid in record_ids}
        if gff is not None:
            features = load_all_features_by_id(
                gff, record_ids, lengths_by_id=lengths
            )

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
        QTimer.singleShot(0, self._try_initial_fit)
        QTimer.singleShot(50, self._try_initial_fit)
        QTimer.singleShot(200, self._try_initial_fit)

        annotated = sum(1 for rid in record_ids if features.get(rid))
        feat_total = sum(len(features.get(rid) or []) for rid in record_ids)
        if gff is not None:
            self.statusMessage.emit(
                f"{len(record_ids):,} sequences · {annotated:,} with annotations · "
                f"{feat_total:,} features · {gff.name}"
            )
        else:
            self.statusMessage.emit(
                f"{len(record_ids):,} sequences · no companion .RV.gff3 / annotation yet"
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

    def _fit_bars(self) -> None:
        # Keep the user's splitter ratio; only zoom tracks to the track viewport.
        self._name_width = self._current_name_width()
        self._sync_viewport_fill()
        self._bars.fit_to_bars_width(self._bars_viewport_width())
        self._bars_scroll.horizontalScrollBar().setValue(0)
        self._sync_names_column()
        self._update_scale()
        self._sync_zoom_slider_from_ppb()

    def _slider_ppb_bounds(self) -> tuple[float, float]:
        max_len = max(1, self._bars.max_length())
        fit_ppb = max(
            _MIN_PPB,
            (max(40, self._bars_viewport_width()) - 16) / max_len,
        )
        max_ppb = _MAX_PPB
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
        delta = _wheel_delta_y(event)
        if delta:
            self._nudge_vertical_scroll(delta)
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
        old_id = store[store_index].gff_id()
        store[store_index] = edited
        store = remount_parent_ids(store, old_id=old_id, new_id=edited.gff_id())
        self._features[rid] = store
        try:
            _json, gff_path = save_annotations_multi(
                self._path,
                features_by_id=self._features,
                lengths_by_id=self._lengths,
            )
        except OSError as exc:
            QMessageBox.critical(self, "Save failed", str(exc))
            return
        keep_ppb = self._bars.pixels_per_base()
        keep_sel = self._bars.selected_rows()
        self._bars.set_data(
            self._record_ids,
            self._lengths,
            self._features,
            fit_bars_width=None,
        )
        self._bars.set_pixels_per_base(keep_ppb)
        self._bars.set_selected_rows(keep_sel)
        self._sync_names_column()
        self._update_scale()
        if gff_path is not None:
            self.annotationsSaved.emit(gff_path)
            self.statusMessage.emit(f"Updated feature on {rid} → {gff_path.name}")

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
