"""Interactive self-comparison dotplot ("Dotplot (Self)" preview tab).

Segments are computed once by ``repeatviewer.backends.dotplot`` and then only
redrawn: zooming and panning never recompute. A new pass happens only when the
word size or the reverse-complement setting changes.

All segments are painted by a single QGraphicsItem rather than one item per
segment. Tens of thousands of QGraphicsLineItems would each carry their own
bounding-box bookkeeping and make panning stutter; one item painting a list of
lines stays smooth. Segments shorter than a pixel at the current zoom are
skipped (level of detail), which keeps a dense satellite plot usable zoomed out.

Rulers are drawn hard against the edges of the plot itself, not the viewport,
so they stay beside the data the way they do in Geneious.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QObject, QPoint, QPointF, QRect, QRectF, Qt, QThread, Signal, Slot
from PySide6.QtGui import QColor, QFont, QPainter, QPalette, QPen, QTransform
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QGraphicsItem,
    QGraphicsScene,
    QGraphicsView,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.spin_boxes import (
    NoWheelDoubleSpinBox as QDoubleSpinBox,
    NoWheelSpinBox as QSpinBox,
)
from repeatviewer.widgets.preview_controls import ZOOM_CONTROL_HEIGHT

from repeatviewer.backends.dotplot import (
    DEFAULT_WORD_SIZE,
    MAX_WORD_SIZE,
    MIN_WORD_SIZE,
    STRAND_REVERSE,
    DotplotError,
    DotplotResult,
    compute_self_dotplot_for_fasta,
    compute_self_dotplot_for_region,
    sequence_length_of_fasta,
    suggest_word_size,
)
from repeatviewer.backends.paf_dotplot import (
    PafError,
    load_alignment_file,
    result_for_all,
    result_for_pair,
)
from repeatviewer.backends.fastga_dotplot import (
    DEFAULT_MIN_IDENTITY,
    DEFAULT_MIN_LENGTH,
    DEFAULT_SEED_FREQ,
    MAX_SEED_FREQ,
    FASTGA_PREFERRED_ABOVE_BYTES,
    FastgaError,
    fastga_available,
    install_hint,
    run_fastga_self_dotplot,
)


class _FastgaWorker(QObject):
    """Runs a FastGA self-comparison off the UI thread.

    Indexing a large genome takes minutes, so this must never run inline. The
    worker polls a cancel flag between output lines and terminates the
    subprocess when asked.
    """

    progress = Signal(str)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(
        self,
        path: Path,
        *,
        threads: int,
        min_identity: float,
        min_length: int,
        seed_freq: int = DEFAULT_SEED_FREQ,
    ) -> None:
        super().__init__()
        self._path = path
        self._threads = threads
        self._min_identity = min_identity
        self._min_length = min_length
        self._seed_freq = seed_freq
        self._cancelled = False

    def cancel(self) -> None:
        self._cancelled = True

    @Slot()
    def run(self) -> None:
        try:
            result = run_fastga_self_dotplot(
                self._path,
                threads=self._threads,
                min_identity=self._min_identity,
                min_length=self._min_length,
                seed_freq=self._seed_freq,
                progress=self.progress.emit,
                is_cancelled=lambda: self._cancelled,
            )
        except FastgaError as exc:
            self.failed.emit(str(exc))
            return
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(f"FastGA failed: {exc}")
            return
        self.finished.emit(result)

_FORWARD_COLOR = QColor(40, 90, 180)   # direct repeats (blue)
_REVERSE_COLOR = QColor(200, 60, 40)   # inverted repeats (red)
# Charcoal — must stay distinct from warm reverse red and cool forward blue.
_DIAGONAL_COLOR = QColor(45, 45, 45)

_DEFAULT_LINE_WIDTH_STEP = 1
_MIN_LINE_WIDTH_STEP = 1
_MAX_LINE_WIDTH_STEP = 15
# UI steps 1…15: same fine gradient as before (≈0.25 px per step), up to ~4.5 px.
_LINE_WIDTH_PX_MIN = 1.0
_LINE_WIDTH_PX_MAX = 4.5


def _line_width_px(step: int) -> float:
    """Map toolbar step (1–5) to cosmetic pen width in screen pixels."""
    step = max(_MIN_LINE_WIDTH_STEP, min(_MAX_LINE_WIDTH_STEP, int(step)))
    if _MAX_LINE_WIDTH_STEP <= _MIN_LINE_WIDTH_STEP:
        return _LINE_WIDTH_PX_MIN
    t = (step - _MIN_LINE_WIDTH_STEP) / (
        _MAX_LINE_WIDTH_STEP - _MIN_LINE_WIDTH_STEP
    )
    return _LINE_WIDTH_PX_MIN + t * (_LINE_WIDTH_PX_MAX - _LINE_WIDTH_PX_MIN)


# Zoom limits are relative to the scale that fits the whole plot, not absolute.
# A 190 kb window fits at ~0.004 px/base, so an absolute lower bound of 0.02
# silently blocked every zoom attempt.
_ZOOM_OUT_LIMIT = 1.0    # never smaller than fit → no empty margins from zoom-out
_MAX_PIXELS_PER_BASE = 40.0
_ZOOM_STEP = 1.12        # +/- button step
# Trackpad: factor = TRACKPAD_ZOOM_BASE ** pixel_delta_y (smooth continuous zoom).
_TRACKPAD_ZOOM_BASE = 1.0022

_RULER_TOP = 18
# Left margin: name strip + tick labels flush against the axis (same px/base as X).
_RULER_LEFT = 48
_Y_NAME_X = 9
_Y_TICK_PAD = 3  # gap between label and tick mark
_RULER_BG = QColor(250, 250, 250)
_RULER_FG = QColor(110, 110, 110)
_RULER_LINE = QColor(200, 200, 200)

_SELECT_FILL = QColor(70, 130, 220, 55)
_SELECT_EDGE = QColor(60, 110, 200, 190)

_CLICK_SLOP = 4


def _nice_tick_step(span: float, target_ticks: int = 10) -> int:
    """A round tick interval (1/2/5 x 10^n) giving roughly ``target_ticks``."""
    if span <= 0:
        return 1
    rough = span / max(target_ticks, 1)
    magnitude = 10 ** int(max(0, len(str(int(rough))) - 1)) if rough >= 1 else 1
    for multiplier in (1, 2, 5, 10):
        step = multiplier * magnitude
        if step >= rough:
            return int(step)
    return int(10 * magnitude)


def _format_bp(value: int) -> str:
    """Always use comma grouping (1,000) — never mix with compact 1k / 1M."""
    return f"{value:,}"


def _scene_to_bp(scene_coord: float | int, *, length: int, offset: int) -> int:
    """Map a scene coordinate in ``[0, length]`` to a 1-based sequence position.

    The plot square is 0-based internally; sequence positions are 1-based, so
    scene 0 is base ``1 + offset`` and scene ``length`` is the last base
    (``length + offset``). Interior ticks keep ``tick + offset`` so round
    labels like 2k stay on the nice step values.
    """
    tick = int(scene_coord)
    if tick <= 0:
        return 1 + offset
    if length > 0 and tick >= length:
        return length + offset
    return tick + offset


def _axis_tick_positions(
    *,
    lo: float,
    hi: float,
    length: int,
    step: int,
) -> list[int]:
    """Major ticks in ``[0, length]`` intersecting the visible span, plus ends.

    Sequence start/end always win over a nearby major tick so labels do not
    collide at the axis tip (e.g. length 810 with step 100 keeps 810, drops 800).
    """
    if step <= 0:
        step = 1
    length = max(int(length), 0)
    ticks: list[int] = []
    tick = int(lo // step) * step
    while tick <= hi + step:
        if 0 <= tick <= length and lo - step <= tick <= hi + step:
            ticks.append(int(tick))
        tick += step

    # Drop majors that would sit under an end label (~half a step).
    min_sep = max(int(step * 0.45), 1)
    for edge in (0, length):
        if not (lo - step <= edge <= hi + step):
            continue
        ticks = [t for t in ticks if abs(t - edge) >= min_sep]
        if edge not in ticks:
            ticks.append(edge)
    ticks.sort()
    return ticks


def _non_overlapping_tick_labels(
    ticks: list[int],
    *,
    length: int,
    offset: int,
    positions: dict[int, float],
    label_extent,
    pad: float = 4.0,
) -> set[int]:
    """Return tick values whose labels fit without overlapping.

    Prefers true sequence ends over intermediate majors. ``label_extent`` is
    ``(pos, label) -> (axis_lo, axis_hi)`` in ruler pixels.
    """
    preferred = [t for t in ticks if t in (0, length)]
    preferred.extend(t for t in ticks if t not in (0, length))
    chosen_spans: list[tuple[float, float]] = []
    show: set[int] = set()
    for tick in preferred:
        pos = positions.get(tick)
        if pos is None:
            continue
        label = _format_bp(_scene_to_bp(tick, length=length, offset=offset))
        a0, a1 = label_extent(pos, label)
        if any(not (a1 + pad <= b0 or b1 + pad <= a0) for b0, b1 in chosen_spans):
            continue
        chosen_spans.append((a0, a1))
        show.add(tick)
    return show


class _SegmentsItem(QGraphicsItem):
    """Paints every dotplot segment in one pass."""

    def __init__(
        self,
        result: DotplotResult,
        *,
        line_width_step: int = _DEFAULT_LINE_WIDTH_STEP,
    ) -> None:
        super().__init__()
        self._result = result
        self._n = max(result.sequence_length, 1)
        self._line_width_step = max(
            _MIN_LINE_WIDTH_STEP, min(_MAX_LINE_WIDTH_STEP, int(line_width_step))
        )
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemUsesExtendedStyleOption, True)

    def set_line_width_step(self, step: int) -> None:
        step = max(_MIN_LINE_WIDTH_STEP, min(_MAX_LINE_WIDTH_STEP, int(step)))
        if step == self._line_width_step:
            return
        self._line_width_step = step
        self.update()

    def boundingRect(self) -> QRectF:  # noqa: N802
        return QRectF(0, 0, float(self._n), float(self._n))

    def _pen(self, color: QColor) -> QPen:
        pen = QPen(color)
        pen.setCosmetic(True)
        pen.setWidthF(_line_width_px(self._line_width_step))
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        return pen

    def paint(self, painter: QPainter, option, widget=None) -> None:  # noqa: N802
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        # White plot face; viewport outside the square uses the ruler grey so
        # there is no framed "box" around the data.
        painter.fillRect(self.boundingRect(), QColor(255, 255, 255))
        scale = painter.transform().m11() or 1.0
        min_length = 0.0 if scale >= 1.0 else (1.0 / scale)
        visible = option.exposedRect if option is not None else self.boundingRect()

        painter.setPen(self._pen(_DIAGONAL_COLOR))
        painter.drawLine(QPointF(0, 0), QPointF(self._n, self._n))

        forward: list[tuple[float, float, float, float]] = []
        reverse: list[tuple[float, float, float, float]] = []
        for seg in self._result.segments:
            if seg.length < min_length:
                continue
            x1, y1 = float(seg.x), float(seg.y)
            x2, y2 = float(seg.end_x), float(seg.end_y)
            if not visible.isNull():
                lo_x, hi_x = (x1, x2) if x1 <= x2 else (x2, x1)
                lo_y, hi_y = (y1, y2) if y1 <= y2 else (y2, y1)
                off = (
                    hi_x < visible.left() or lo_x > visible.right()
                    or hi_y < visible.top() or lo_y > visible.bottom()
                )
                if off:
                    if seg.strand == STRAND_REVERSE:
                        continue
                    mirrored_off = (
                        hi_y < visible.left() or lo_y > visible.right()
                        or hi_x < visible.top() or lo_x > visible.bottom()
                    )
                    if mirrored_off:
                        continue
            if seg.strand == STRAND_REVERSE:
                reverse.append((x1, y1, x2, y2))
            else:
                forward.append((x1, y1, x2, y2))

        painter.setPen(self._pen(_FORWARD_COLOR))
        for x1, y1, x2, y2 in forward:
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
            painter.drawLine(QPointF(y1, x1), QPointF(y2, x2))  # self-symmetry

        painter.setPen(self._pen(_REVERSE_COLOR))
        for x1, y1, x2, y2 in reverse:
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))


class _Ruler(QWidget):
    """Coordinate ruler living in the view's reserved margin.

    A child widget rather than painting inside QGraphicsView.paintEvent,
    because that paint event belongs to the viewport and cannot reach the
    margins. This is the same pattern as Qt's CodeEditor line-number area.
    """

    def __init__(self, view: "DotplotCanvas", horizontal: bool) -> None:
        super().__init__(view)
        self._view = view
        self._horizontal = horizontal

    def paintEvent(self, event) -> None:  # noqa: N802
        view = self._view
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        painter.fillRect(self.rect(), _RULER_BG)
        if view.plot_length() <= 0:
            return

        font = QFont(self.font())
        font.setPointSize(8)
        painter.setFont(font)

        entries = view.axis_entries(self._horizontal)
        if entries:
            painter.setPen(QPen(_RULER_FG))
            for entry in entries:
                centre = entry.offset + entry.length / 2.0
                pos = view.scene_to_ruler(centre, self._horizontal)
                text = entry.name
                if self._horizontal:
                    if 0 <= pos <= self.width():
                        w = painter.fontMetrics().horizontalAdvance(text)
                        painter.drawText(int(pos - w / 2), self.height() - 5, text)
                else:
                    if 0 <= pos <= self.height():
                        painter.save()
                        painter.translate(self.width() - 4, pos)
                        painter.rotate(-90)
                        w = painter.fontMetrics().horizontalAdvance(text)
                        painter.drawText(int(-w / 2), 0, text)
                        painter.restore()
            painter.setPen(QPen(_RULER_LINE))
            if self._horizontal:
                painter.drawLine(0, self.height() - 1, self.width(), self.height() - 1)
            else:
                painter.drawLine(self.width() - 1, 0, self.width() - 1, self.height())
            return

        visible = view.visible_scene_rect()
        length = view.plot_length()
        # Same tick interval on X and Y so equal scale is obvious (a wide
        # viewport would otherwise pick step 100 on X and 20 on Y).
        step = view.shared_tick_step()
        lo = visible.left() if self._horizontal else visible.top()
        hi = visible.right() if self._horizontal else visible.bottom()

        painter.setPen(QPen(_RULER_FG))
        metrics = painter.fontMetrics()
        offset = view.coordinate_offset()
        ticks = _axis_tick_positions(lo=lo, hi=hi, length=length, step=step)
        positions = {
            tick: view.scene_to_ruler(tick, self._horizontal) for tick in ticks
        }

        def _label_extent(pos: float, label: str) -> tuple[float, float]:
            if self._horizontal:
                # Horizontal labels start just right of the tick mark.
                width = float(metrics.horizontalAdvance(label))
                return pos + 2.0, pos + 2.0 + width
            height = float(metrics.height())
            return pos - height / 2.0, pos + height / 2.0

        labeled = _non_overlapping_tick_labels(
            ticks,
            length=length,
            offset=offset,
            positions=positions,
            label_extent=_label_extent,
        )

        for tick in ticks:
            pos = positions[tick]
            label = _format_bp(_scene_to_bp(tick, length=length, offset=offset))
            if self._horizontal:
                if 0 <= pos <= self.width():
                    painter.drawLine(int(pos), self.height() - 4, int(pos), self.height())
                    if tick in labeled:
                        painter.drawText(
                            int(pos) + 2,
                            self.height() - 6,
                            label,
                        )
            else:
                if 0 <= pos <= self.height():
                    painter.drawLine(self.width() - 4, int(pos), self.width(), int(pos))
                    if tick not in labeled:
                        continue
                    # Right-align and vertically centre on the tick so the
                    # number sits on the mark (drawText y is baseline).
                    label_h = metrics.height()
                    text_rect = QRect(
                        0,
                        int(pos) - label_h // 2,
                        self.width() - 4 - _Y_TICK_PAD,
                        label_h,
                    )
                    painter.drawText(
                        text_rect,
                        int(
                            Qt.AlignmentFlag.AlignRight
                            | Qt.AlignmentFlag.AlignVCenter
                        ),
                        label,
                    )

        cursor = view.cursor_scene_pos()
        if cursor is not None:
            painter.setPen(QPen(QColor(40, 90, 180)))
            value = cursor.x() if self._horizontal else cursor.y()
            pos = view.scene_to_ruler(value, self._horizontal)
            real = _scene_to_bp(value, length=length, offset=offset)
            label = f"{real:,}"
            if self._horizontal and 0 <= pos <= self.width():
                painter.drawText(int(pos) + 2, self.height() - 6, label)
            elif not self._horizontal and 0 <= pos <= self.height():
                label_h = metrics.height()
                text_rect = QRect(
                    0,
                    int(pos) - label_h // 2,
                    self.width() - 4 - _Y_TICK_PAD,
                    label_h,
                )
                painter.drawText(
                    text_rect,
                    int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter),
                    label,
                )

        # Sequence name down the y axis, packed against the tick labels.
        if not self._horizontal and view.plot_title():
            painter.setPen(QPen(_RULER_FG))
            painter.save()
            painter.translate(_Y_NAME_X, self.height() / 2.0)
            painter.rotate(-90)
            painter.drawText(
                QRectF(-self.height() / 2.0, -8, self.height(), 16),
                int(Qt.AlignmentFlag.AlignCenter),
                view.plot_title(),
            )
            painter.restore()

        painter.setPen(QPen(_RULER_LINE))
        if self._horizontal:
            painter.drawLine(0, self.height() - 1, self.width(), self.height() - 1)
        else:
            painter.drawLine(self.width() - 1, 0, self.width() - 1, self.height())


class DotplotCanvas(QGraphicsView):
    """Zoom/pan canvas with rulers, rectangular selection and a crosshair."""

    positionChanged = Signal(int, int)
    rectSelected = Signal(int, int, int, int)  # x_start, x_end, y_start, y_end

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setScene(QGraphicsScene(self))
        self.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        self.setDragMode(QGraphicsView.DragMode.NoDrag)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        # Anchor the plot to the top-left like Geneious, instead of floating it
        # in the middle of a wide panel with empty space on both sides.
        self.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        # Reserve real estate for the rulers *outside* the viewport. Previously
        # the rulers were painted over the viewport and room was faked by
        # padding the fitInView() rectangle, which fought with the alignment
        # and left the plot in the wrong place. With margins the plot area and
        # the ruler area are separate regions and the geometry is exact.
        self.setViewportMargins(_RULER_LEFT, _RULER_TOP, 0, 0)
        # Scrollbars stay available: panning moves them, and with them disabled
        # a zoomed-in plot cannot be scrolled at all.
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setMouseTracking(True)
        self.setBackgroundBrush(_RULER_BG)
        self._item: _SegmentsItem | None = None
        self._length = 0
        self._line_width_step = _DEFAULT_LINE_WIDTH_STEP
        # Plot coordinates are 0..length; add this to show real positions when
        # only a window of a large sequence is plotted.
        self._offset = 0
        self._fit_scale = 0.0
        # Sequence boundaries for the concatenated (all sequences) layout.
        self._axis_x: list = []
        self._axis_y: list = []
        self._press_pos: QPoint | None = None
        self._panning = False
        self._pan_origin: QPoint | None = None
        self._drag_start: QPointF | None = None
        self._drag_now: QPointF | None = None
        self._cursor_pos: QPointF | None = None
        # Exactly one selection exists at a time, always a rectangle in scene
        # coordinates, whether it came from a drag here or from the Annotation
        # tab. Two different-looking selections that both persisted was the
        # source of much confusion.
        self._selection: QRectF | None = None
        self._needs_fit = False
        self._title = ""
        self._top_ruler = _Ruler(self, horizontal=True)
        self._left_ruler = _Ruler(self, horizontal=False)
        # The viewport margins leave a hole at (0,0)–(_RULER_LEFT,_RULER_TOP)
        # that neither ruler covers; fill it with the same grey as the scales.
        self._ruler_corner = QWidget(self)
        self._ruler_corner.setAutoFillBackground(True)
        corner_pal = self._ruler_corner.palette()
        corner_pal.setColor(QPalette.ColorRole.Window, _RULER_BG)
        self._ruler_corner.setPalette(corner_pal)
        # Repaint the rulers whenever the view scrolls, so ticks track the plot.
        self.horizontalScrollBar().valueChanged.connect(self._refresh_rulers)
        self.verticalScrollBar().valueChanged.connect(self._refresh_rulers)

    def set_title(self, title: str) -> None:
        """Sequence name, drawn rotated beside the y axis."""
        self._title = title or ""
        self._refresh_rulers()

    # -- accessors used by the ruler widgets -------------------------------

    def plot_length(self) -> int:
        return self._length

    def coordinate_offset(self) -> int:
        return self._offset

    def set_axes(self, x_entries: list, y_entries: list) -> None:
        """Sequence layout for the concatenated view (empty for a single pair)."""
        self._axis_x = list(x_entries or [])
        self._axis_y = list(y_entries or [])
        self.viewport().update()
        self._refresh_rulers()

    def axis_entries(self, horizontal: bool) -> list:
        return self._axis_x if horizontal else self._axis_y

    def plot_title(self) -> str:
        return self._title

    def cursor_scene_pos(self) -> QPointF | None:
        return self._cursor_pos

    def visible_scene_rect(self) -> QRectF:
        return self.mapToScene(self.viewport().rect()).boundingRect()

    def scene_to_ruler(self, value: float, horizontal: bool) -> float:
        """Scene coordinate -> position along the matching ruler widget.

        The rulers are laid out flush with the viewport on their long axis, so
        the viewport coordinate is already the ruler coordinate.
        """
        point = self.mapFromScene(QPointF(value, value))
        return point.x() if horizontal else point.y()

    def set_line_width_step(self, step: int) -> None:
        """Toolbar step 1–15 maps evenly onto ~1–4.5 screen pixels."""
        step = max(_MIN_LINE_WIDTH_STEP, min(_MAX_LINE_WIDTH_STEP, int(step)))
        self._line_width_step = step
        if self._item is not None:
            self._item.set_line_width_step(step)

    def line_width_step(self) -> int:
        return self._line_width_step

    def shared_tick_step(self) -> int:
        """One tick interval for both axes (equal scale stays visually clear)."""
        visible = self.visible_scene_rect()
        length = max(float(self._length), 1.0)
        # Intersect with the sequence square so empty letterbox beside a fitted
        # plot does not inflate the span (that produced only 1 and 5,000 on a 6 kb
        # sequence in a wide window).
        shown = visible.intersected(QRectF(0.0, 0.0, length, length))
        span = max(float(shown.width()), float(shown.height()), 1.0)
        return _nice_tick_step(span, target_ticks=10)

    def _refresh_rulers(self) -> None:
        self._top_ruler.update()
        self._left_ruler.update()

    def _position_rulers(self) -> None:
        vp = self.viewport().geometry()
        # Only the top-left margin corner — plot stays white.
        self._ruler_corner.setGeometry(0, 0, max(vp.x(), 1), max(vp.y(), 1))
        self._ruler_corner.raise_()
        self._top_ruler.setGeometry(vp.x(), 0, vp.width(), _RULER_TOP)
        self._left_ruler.setGeometry(0, vp.y(), _RULER_LEFT, vp.height())
        self._top_ruler.raise_()
        self._left_ruler.raise_()

    # ------------------------------------------------------------------
    # content
    # ------------------------------------------------------------------

    def show_result(self, result: DotplotResult, *, preserve_view: bool = False) -> None:
        same_length = (
            preserve_view
            and self._length == result.sequence_length
            and self._length > 0
            and self._item is not None
        )
        saved_transform = self.transform() if same_length else None
        saved_fit = self._fit_scale if same_length else 0.0
        saved_center = (
            self.mapToScene(self.viewport().rect().center())
            if same_length
            else None
        )

        scene = self.scene()
        scene.clear()
        self._item = _SegmentsItem(result, line_width_step=self._line_width_step)
        self._length = result.sequence_length
        self._offset = int(getattr(result, "coordinate_offset", 0) or 0)
        self._selection = None
        scene.addItem(self._item)
        scene.setSceneRect(0, 0, result.sequence_length, result.sequence_length)

        if saved_transform is not None and saved_center is not None:
            # Keep pan/zoom when only word size / reverse-complement toggled.
            self.setTransform(saved_transform)
            self._fit_scale = saved_fit
            self.centerOn(saved_center)
            self._needs_fit = False
            self.viewport().update()
            self._refresh_rulers()
            return

        # The tab is often still hidden here, so the viewport has no usable size
        # yet and fitInView() would scale to almost nothing. Defer until we are
        # actually shown or resized.
        self._needs_fit = True
        self._fit_if_needed()

    def clear_plot(self) -> None:
        self.scene().clear()
        self._item = None
        self._length = 0
        self._selection = None

    def _fit_if_needed(self) -> None:
        if not self._needs_fit or self._item is None:
            return
        viewport = self.viewport().rect()
        if viewport.width() < 40 or viewport.height() < 40:
            return  # not laid out yet; try again on the next resize/show
        self._needs_fit = False
        self.reset_zoom()

    def reset_zoom(self) -> None:
        """Fit the whole sequence square into the plot area.

        X and Y use the same pixels-per-base so a self-dotplot is always a
        true square (45° diagonal). ``fitInView`` is avoided: with viewport
        margins it can leave a non-uniform transform and stretch the plot.
        """
        if not self._length:
            return
        vp = self.viewport().rect()
        if vp.width() < 1 or vp.height() < 1:
            return
        n = float(self._length)
        scale = min(vp.width() / n, vp.height() / n)
        self._set_uniform_scale(scale)
        self._fit_scale = scale
        self._clamp_scroll_to_plot()
        self.viewport().update()
        self._refresh_rulers()

    def _set_uniform_scale(self, scale: float) -> None:
        """Identical X/Y scale, plot pinned to the top-left of the viewport."""
        scale = max(float(scale), 1e-12)
        self.setTransform(QTransform.fromScale(scale, scale))
        self.horizontalScrollBar().setValue(0)
        self.verticalScrollBar().setValue(0)

    def _uniform_scale(self) -> float:
        """Current pixels-per-base (forces X/Y equal if they ever drift)."""
        t = self.transform()
        sx = abs(t.m11()) or 1.0
        sy = abs(t.m22()) or sx
        if abs(sx - sy) > 1e-9 * max(sx, sy):
            # Repair a non-uniform matrix without changing the intended zoom.
            scale = min(sx, sy)
            self.setTransform(QTransform.fromScale(scale, scale))
            return scale
        return sx

    def zoom_to_selection(self) -> bool:
        """Fit the current selection into the view. Returns False if none."""
        rect = self._selection
        if rect is None or not self._length:
            return False
        if rect.width() < 1.0 or rect.height() < 1.0:
            return False
        pad = max(rect.width(), rect.height()) * 0.05
        target = rect.adjusted(-pad, -pad, pad, pad).intersected(
            QRectF(0, 0, float(self._length), float(self._length))
        )
        if target.width() < 1.0 or target.height() < 1.0:
            return False
        vp = self.viewport().rect()
        if vp.width() < 1 or vp.height() < 1:
            return False
        # Same px/base on both axes; size to the larger side of the selection.
        side = max(target.width(), target.height())
        scale = min(vp.width() / side, vp.height() / side)
        scale = min(scale, _MAX_PIXELS_PER_BASE)
        self.setTransform(QTransform.fromScale(scale, scale))
        self.centerOn(target.center())
        self._clamp_scroll_to_plot()
        self.viewport().update()
        self._refresh_rulers()
        return True

    def has_selection(self) -> bool:
        rect = self._selection
        return rect is not None and rect.width() >= 1.0 and rect.height() >= 1.0

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._position_rulers()
        if self._needs_fit:
            self._fit_if_needed()
        elif self._item is not None and self._fit_scale > 0:
            # Stay locked to a square fit when the user hasn't zoomed in.
            current = self._uniform_scale()
            if abs(current - self._fit_scale) <= max(1e-6, self._fit_scale * 1e-4):
                self.reset_zoom()
        self._refresh_rulers()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._position_rulers()
        self._fit_if_needed()
        self._refresh_rulers()

    # ------------------------------------------------------------------
    # selection
    # ------------------------------------------------------------------

    def set_range(self, start: int, end: int) -> None:
        """Show a sequence range coming from the Annotation tab.

        A 1-D range maps onto the diagonal block [start,end] x [start,end],
        which is that region compared against itself -- the same rectangle you
        would draw by hand to inspect it, so both sources look identical.
        """
        if not self._length:
            return
        lo, hi = (start, end) if start <= end else (end, start)
        lo = max(0, lo - 1 - self._offset)
        hi = min(self._length, hi - self._offset)
        if hi <= lo:
            return
        self._selection = QRectF(lo, lo, hi - lo, hi - lo)
        self.viewport().update()

    def clear_selection(self) -> None:
        self._selection = None
        self._drag_start = self._drag_now = None
        self.viewport().update()

    def selection_sequence_coords(self) -> tuple[int, int, int, int] | None:
        """Return the current drag-box as 1-based inclusive (x1, x2, y1, y2), or None."""
        if self._selection is None or self._length <= 0:
            return None
        rect = self._selection
        if rect.width() < 1.0 or rect.height() < 1.0:
            return None
        return (
            int(rect.left()) + 1 + self._offset,
            int(rect.right()) + self._offset,
            int(rect.top()) + 1 + self._offset,
            int(rect.bottom()) + self._offset,
        )

    # ------------------------------------------------------------------
    # painting
    # ------------------------------------------------------------------

    def drawForeground(self, painter: QPainter, rect: QRectF) -> None:  # noqa: N802
        """Scene-space overlays only: selection and crosshair lines.

        Anything that belongs in the ruler margins is painted in paintEvent
        instead, because the margins are outside this viewport.
        """
        super().drawForeground(painter, rect)
        if self._item is None or not self._length:
            return

        active = self._selection
        if self._drag_start is not None and self._drag_now is not None:
            active = QRectF(self._drag_start, self._drag_now).normalized()
        if active is not None:
            painter.fillRect(active, _SELECT_FILL)
            edge = QPen(_SELECT_EDGE)
            edge.setCosmetic(True)
            painter.setPen(edge)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(active)

        if self._axis_x or self._axis_y:
            sep = QPen(QColor(150, 150, 150))
            sep.setCosmetic(True)
            painter.setPen(sep)
            for entry in self._axis_x[1:]:
                painter.drawLine(
                    QPointF(entry.offset, 0), QPointF(entry.offset, self._length)
                )
            for entry in self._axis_y[1:]:
                painter.drawLine(
                    QPointF(0, entry.offset), QPointF(self._length, entry.offset)
                )

        if self._cursor_pos is not None:
            cross = QPen(QColor(90, 140, 220, 150))
            cross.setCosmetic(True)
            painter.setPen(cross)
            cx = self._cursor_pos.x()
            cy = self._cursor_pos.y()
            painter.drawLine(QPointF(cx, 0), QPointF(cx, self._length))
            painter.drawLine(QPointF(0, cy), QPointF(self._length, cy))

    def _zoom_allowed(self, target: float) -> bool:
        floor = (self._fit_scale or target) * _ZOOM_OUT_LIMIT
        return floor <= target <= _MAX_PIXELS_PER_BASE

    def _clamp_scroll_to_plot(self) -> None:
        """Keep empty scene space out of the viewport (pin plot to top-left)."""
        if not self._length:
            return
        n = float(self._length)
        vis = self.mapToScene(self.viewport().rect()).boundingRect()
        dx = dy = 0.0
        if vis.width() <= n + 1e-6:
            if vis.left() < 0:
                dx = -vis.left()
            elif vis.right() > n:
                dx = n - vis.right()
        else:
            # Viewport wider than the plot: keep x=0 on the left edge.
            dx = -vis.left()
        if vis.height() <= n + 1e-6:
            if vis.top() < 0:
                dy = -vis.top()
            elif vis.bottom() > n:
                dy = n - vis.bottom()
        else:
            dy = -vis.top()
        if abs(dx) < 1e-6 and abs(dy) < 1e-6:
            return
        origin = self.mapFromScene(QPointF(0, 0))
        shifted = self.mapFromScene(QPointF(dx, dy))
        self.horizontalScrollBar().setValue(
            self.horizontalScrollBar().value() + int(round(shifted.x() - origin.x()))
        )
        self.verticalScrollBar().setValue(
            self.verticalScrollBar().value() + int(round(shifted.y() - origin.y()))
        )

    def _apply_zoom(self, factor: float, anchor: QPoint | None = None) -> None:
        """Zoom by ``factor`` about ``anchor`` (viewport coords), then clamp."""
        if self._item is None or not self._length:
            return
        current = self._uniform_scale()
        floor = (self._fit_scale or current) * _ZOOM_OUT_LIMIT
        target = min(_MAX_PIXELS_PER_BASE, max(floor, current * factor))
        actual = target / current
        if abs(actual - 1.0) < 1e-9:
            return

        if anchor is None:
            anchor = self.viewport().rect().center()

        # Manual under-cursor zoom. Do not combine with AnchorUnderMouse —
        # that double-applies the correction and drifts the plot, leaving a
        # blank band at the top/left (especially with a Mac trackpad).
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.NoAnchor)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.NoAnchor)
        scene_pt = self.mapToScene(anchor)
        self.scale(actual, actual)
        offset = self.mapFromScene(scene_pt) - anchor
        self.horizontalScrollBar().setValue(
            self.horizontalScrollBar().value() + int(round(offset.x()))
        )
        self.verticalScrollBar().setValue(
            self.verticalScrollBar().value() + int(round(offset.y()))
        )
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self._clamp_scroll_to_plot()
        self.viewport().update()
        self._refresh_rulers()

    def zoom_by(self, factor: float) -> None:
        """Zoom about the centre of the view (used by the +/- buttons)."""
        self._apply_zoom(factor, self.viewport().rect().center())

    def zoom_at(self, factor: float, view_pos) -> None:
        """Zoom keeping the scene point under ``view_pos`` fixed."""
        self._apply_zoom(factor, view_pos)

    def render_to_pdf(self, path: Path, *, title: str = "", subtitle: str = "") -> None:
        """Write the plot to a PDF with axes, labels and a legend.

        The on-screen rulers are child widgets in the view margins, so they are
        not part of the scene and cannot simply be rendered. The axes here are
        drawn directly, using the same tick spacing and coordinate offset as the
        interactive rulers so the two agree.
        """
        from PySide6.QtGui import QPageSize, QPdfWriter

        writer = QPdfWriter(str(path))
        writer.setPageSize(QPageSize(QPageSize.PageSizeId.A4))
        writer.setResolution(300)
        writer.setTitle(title or "Dotplot")

        painter = QPainter(writer)
        try:
            page = painter.viewport()
            unit = page.width() / 210.0          # device units per mm (A4 width)
            margin = int(12 * unit)
            gutter = int(16 * unit)              # room for tick labels
            head = int(26 * unit)                # title block + top tick labels
            foot = int(10 * unit)                # legend

            side = min(
                page.width() - margin * 2 - gutter,
                page.height() - margin * 2 - head - foot - gutter,
            )
            plot = QRectF(margin + gutter, margin + head, side, side)

            # The app may be in a dark theme; a figure should not be. Paint the
            # page and the plot area white and draw everything in dark ink.
            painter.fillRect(page, QColor(255, 255, 255))

            base_font = QFont(painter.font())

            if title:
                f = QFont(base_font); f.setPointSizeF(11); f.setBold(True)
                painter.setFont(f)
                painter.setPen(QPen(QColor(20, 20, 20)))
                painter.drawText(int(margin), int(margin + 6 * unit), title)
            if subtitle:
                f = QFont(base_font); f.setPointSizeF(8)
                painter.setFont(f)
                painter.setPen(QPen(QColor(90, 90, 90)))
                painter.drawText(int(margin), int(margin + 12 * unit), subtitle)

            # --- the plot itself ---
            painter.fillRect(plot, QColor(255, 255, 255))
            painter.save()
            painter.setClipRect(plot)
            self.scene().render(painter, plot, self.sceneRect())
            painter.restore()

            # --- axes, matching the on-screen rulers ---
            length = max(self._length, 1)
            offset = self._offset
            step = self.shared_tick_step() if self._length else _nice_tick_step(float(length))
            tick_len = 3 * unit
            f = QFont(base_font); f.setPointSizeF(7)
            painter.setFont(f)
            painter.setPen(QPen(QColor(70, 70, 70)))
            metrics = painter.fontMetrics()

            ticks = _axis_tick_positions(
                lo=0.0, hi=float(length), length=length, step=step
            )
            x_positions = {
                tick: plot.left() + (tick / length) * plot.width() for tick in ticks
            }
            y_positions = {
                tick: plot.top() + (tick / length) * plot.height() for tick in ticks
            }

            def _export_x_extent(pos: float, label: str) -> tuple[float, float]:
                half = metrics.horizontalAdvance(label) / 2.0
                return pos - half, pos + half

            def _export_y_extent(pos: float, label: str) -> tuple[float, float]:
                half = metrics.height() / 2.0
                return pos - half, pos + half

            x_labeled = _non_overlapping_tick_labels(
                ticks,
                length=length,
                offset=offset,
                positions=x_positions,
                label_extent=_export_x_extent,
            )
            y_labeled = _non_overlapping_tick_labels(
                ticks,
                length=length,
                offset=offset,
                positions=y_positions,
                label_extent=_export_y_extent,
            )

            for tick in ticks:
                label = _format_bp(_scene_to_bp(tick, length=length, offset=offset))
                x = x_positions[tick]
                painter.drawLine(int(x), int(plot.top() - tick_len), int(x), int(plot.top()))
                if tick in x_labeled:
                    painter.drawText(
                        int(x - metrics.horizontalAdvance(label) / 2),
                        int(plot.top() - tick_len - 2 * unit),
                        label,
                    )

                y = y_positions[tick]
                painter.drawLine(int(plot.left() - tick_len), int(y), int(plot.left()), int(y))
                if tick in y_labeled:
                    painter.drawText(
                        int(
                            plot.left()
                            - tick_len
                            - 2 * unit
                            - metrics.horizontalAdvance(label)
                        ),
                        int(y + metrics.ascent() / 2),
                        label,
                    )

            # sequence name down the y axis
            if title:
                painter.save()
                painter.translate(margin + 4 * unit, plot.center().y())
                painter.rotate(-90)
                painter.drawText(
                    int(-metrics.horizontalAdvance(title) / 2), 0, title
                )
                painter.restore()

            # --- legend ---
            legend_y = plot.bottom() + 8 * unit
            swatch = 6 * unit
            entries = (
                (_DIAGONAL_COLOR, "self"),
                (_FORWARD_COLOR, "direct repeat"),
                (_REVERSE_COLOR, "inverted repeat"),
            )
            x = plot.left()
            for colour, label in entries:
                pen = QPen(colour); pen.setWidth(max(1, int(0.4 * unit)))
                painter.setPen(pen)
                painter.drawLine(int(x), int(legend_y), int(x + swatch), int(legend_y))
                painter.setPen(QPen(QColor(70, 70, 70)))
                painter.drawText(
                    int(x + swatch + 2 * unit),
                    int(legend_y + metrics.ascent() / 2),
                    label,
                )
                x += swatch + 4 * unit + metrics.horizontalAdvance(label)
        finally:
            painter.end()

    def wheelEvent(self, event) -> None:  # noqa: N802
        if self._item is None:
            return
        # Mac trackpads send smooth pixelDelta; mouse wheels use angleDelta notches.
        pixel = event.pixelDelta()
        angle = event.angleDelta()
        if not pixel.isNull() and (pixel.y() != 0 or pixel.x() != 0):
            dy = float(pixel.y() if pixel.y() != 0 else pixel.x())
            factor = _TRACKPAD_ZOOM_BASE ** dy
        elif angle.y() != 0:
            factor = _ZOOM_STEP ** (angle.y() / 120.0)
        else:
            return
        self._apply_zoom(factor, event.position().toPoint())
        event.accept()

    def mousePressEvent(self, event) -> None:  # noqa: N802
        pos = event.position().toPoint()
        self._press_pos = pos
        pan_requested = (
            event.button() == Qt.MouseButton.MiddleButton
            or bool(event.modifiers() & Qt.KeyboardModifier.ShiftModifier)
        )
        if self._item is not None and pan_requested:
            self._panning = True
            self._pan_origin = pos
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
            return
        if self._item is not None and event.button() == Qt.MouseButton.LeftButton:
            # Starting a new drag drops the previous selection straight away,
            # so only one is ever on screen.
            self._selection = None
            self._drag_start = self.mapToScene(pos)
            self._drag_now = self._drag_start
            self.viewport().update()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        pos = event.position().toPoint()
        if self._panning and self._pan_origin is not None:
            delta = pos - self._pan_origin
            self._pan_origin = pos
            h = self.horizontalScrollBar()
            v = self.verticalScrollBar()
            h.setValue(h.value() - delta.x())
            v.setValue(v.value() - delta.y())
            self._clamp_scroll_to_plot()
            self.viewport().update()
            self._refresh_rulers()
            return
        super().mouseMoveEvent(event)
        if self._item is None:
            return
        point = self.mapToScene(pos)
        self._cursor_pos = point
        if self._drag_start is not None:
            self._drag_now = point
        self.viewport().update()
        self._refresh_rulers()
        x, y = int(point.x()), int(point.y())
        if 0 <= x <= self._length and 0 <= y <= self._length:
            self.positionChanged.emit(
                _scene_to_bp(x, length=self._length, offset=self._offset),
                _scene_to_bp(y, length=self._length, offset=self._offset),
            )

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if self._panning:
            self._panning = False
            self._pan_origin = None
            self.unsetCursor()
            self._press_pos = None
            return
        super().mouseReleaseEvent(event)
        if self._item is None or self._press_pos is None:
            self._drag_start = self._drag_now = None
            return
        released = event.position().toPoint()
        moved = (released - self._press_pos).manhattanLength()
        self._press_pos = None

        if moved <= _CLICK_SLOP:
            # A plain click only clears the selection; it never navigates.
            self._drag_start = self._drag_now = None
            self._selection = None
            self.viewport().update()
            return

        if self._drag_start is not None and self._drag_now is not None:
            rect = QRectF(self._drag_start, self._drag_now).normalized()
            self._drag_start = self._drag_now = None
            bounds = QRectF(0, 0, self._length, self._length)
            rect = rect.intersected(bounds)
            if rect.width() >= 1 and rect.height() >= 1:
                self._selection = rect
                self.viewport().update()
                self.rectSelected.emit(
                    int(rect.left()) + 1 + self._offset,
                    int(rect.right()) + self._offset,
                    int(rect.top()) + 1 + self._offset,
                    int(rect.bottom()) + self._offset,
                )

    def leaveEvent(self, event) -> None:  # noqa: N802
        super().leaveEvent(event)
        self._cursor_pos = None
        self.viewport().update()
        self._refresh_rulers()


class DotplotPanel(QWidget):
    """Controls plus canvas, in the spirit of the Geneious 'Dotplot (Self)' tab."""

    statusMessage = Signal(str)
    # Rectangular drag-selection, 1-based inclusive on each axis.
    rangeSelected = Signal(int, int, int, int)
    # "View select" on x or y: 1-based inclusive (start, end) for that axis.
    viewInAnnotationRequested = Signal(int, int)

    def __init__(self, path: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("DotplotPanel")
        self._path = path.expanduser().resolve()
        self._loaded_for: tuple | None = None
        self._sequence_id = ""
        self._sequence_length = 0
        self._full_length = 0
        # Optional override for multi-record FASTA (active Annotation sequence).
        self._active_sequence_id = ""
        self._active_sequence = ""
        # Large sequences (>~2 Mb) show a status note instead of the instant
        # k-mer window. FastGA is available from the main toolbar, not this tab.
        self._engine_choice = "fastga"
        self._too_long_prompt_active = False
        # "fasta" computes a dotplot; "alignment" displays a PAF/.1aln file.
        self._mode = "fasta"
        self._alignments = None

        self._word = QSpinBox(self)
        self._word.setRange(MIN_WORD_SIZE, MAX_WORD_SIZE)
        self._word.setValue(DEFAULT_WORD_SIZE)
        self._word.setToolTip(
            "Length of the exact match used to place a dot. Larger = stricter "
            "and faster; smaller = more sensitive and noisier."
        )
        self._word.setFixedWidth(70)

        self._reverse = QCheckBox("Show reverse", self)
        self._reverse.setChecked(True)
        self._reverse.setToolTip(
            "Show reverse-strand hits (inverted repeats, TIRs, hairpins) "
            "as anti-diagonal lines."
        )

        self._line_width = QSpinBox(self)
        self._line_width.setRange(_MIN_LINE_WIDTH_STEP, _MAX_LINE_WIDTH_STEP)
        self._line_width.setValue(_DEFAULT_LINE_WIDTH_STEP)
        self._line_width.setFixedWidth(55)
        self._line_width.setToolTip(
            "Match and self-diagonal thickness.\n"
            "Steps 1–15 use a fine gradient "
            f"({_LINE_WIDTH_PX_MIN:g}–{_LINE_WIDTH_PX_MAX:g} screen pixels)."
        )
        self._line_width_label = QLabel("Line width:", self)

        # FastGA controls: only meaningful for the large-sequence engine.
        self._identity = QDoubleSpinBox(self)
        self._identity.setRange(0.5, 1.0)
        self._identity.setSingleStep(0.05)
        self._identity.setDecimals(2)
        self._identity.setValue(DEFAULT_MIN_IDENTITY)
        self._identity.setFixedWidth(70)
        self._identity.setToolTip("Minimum alignment identity reported by FastGA.")

        self._min_len = QSpinBox(self)
        # FastGA rejects values below its seed/chain thresholds: alignments must
        # contain a 40 bp seed and chains must cover at least 85 bases, so a
        # very small -l simply fails. 100 is FastGA's own default.
        self._min_len.setRange(100, 100_000)
        self._min_len.setSingleStep(50)
        self._min_len.setValue(DEFAULT_MIN_LENGTH)
        self._min_len.setFixedWidth(90)
        self._min_len.setToolTip(
            "Shortest alignment FastGA will report, in bp.\n"
            "Cannot go below 100: alignments must contain a seed, and FastGA "
            "rejects shorter limits.\n"
            "Raise it on satellite arrays, where short matches swamp the plot."
        )

        self._seed_freq = QSpinBox(self)
        self._seed_freq.setRange(1, MAX_SEED_FREQ)
        self._seed_freq.setValue(DEFAULT_SEED_FREQ)
        self._seed_freq.setSingleStep(50)
        self._seed_freq.setFixedWidth(90)
        self._seed_freq.setToolTip(
            "FastGA's -f: seeds occurring more often than this are discarded as "
            "'repetitive' before alignment.\n"
            "FastGA's own default is 10, which removes satellite arrays entirely "
            "-- every monomer k-mer occurs far more than 10 times.\n"
            "Raise it to see repeat structure; higher values are slower."
        )

        self._region_start = QSpinBox(self)
        self._region_start.setRange(1, 2_000_000_000)
        self._region_start.setValue(1)
        self._region_start.setGroupSeparatorShown(True)
        self._region_start.setFixedWidth(120)
        self._region_start.setToolTip("First base of the window to plot.")

        self._region_size = QSpinBox(self)
        self._region_size.setRange(1_000, 5_000_000)
        self._region_size.setValue(200_000)
        self._region_size.setSingleStep(50_000)
        self._region_size.setGroupSeparatorShown(True)
        self._region_size.setFixedWidth(120)
        self._region_size.setToolTip(
            "How many bases to plot from the start position.\n"
            "A whole scaffold cannot be plotted at small word sizes -- the "
            "match count grows with the square of the repeat copy number -- "
            "but a window of a few hundred kb is instant."
        )

        self._region_start_label = QLabel("From:", self)
        self._region_size_label = QLabel("Window:", self)

        self._pair_combo = QComboBox(self)
        # Long pair labels live in the popup; keep the closed control compact.
        self._pair_combo.setFixedWidth(180)
        self._pair_combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self._pair_combo.setMinimumContentsLength(8)
        self._pair_combo.setToolTip(
            "Which sequence pair to plot. 'All sequences' lays every query and "
            "target end to end, with separators between them."
        )
        self._pair_label = QLabel("Compare:", self)

        # For an alignment file there is no k-mer step to tune, so the
        # equivalent of word size is which alignments are drawn at all.
        self._filter_len = QSpinBox(self)
        self._filter_len.setRange(0, 10_000_000)
        self._filter_len.setValue(0)
        self._filter_len.setSingleStep(100)
        self._filter_len.setGroupSeparatorShown(True)
        self._filter_len.setFixedWidth(72)
        self._filter_len.setToolTip(
            "Hide alignments shorter than this, in bp.\n"
            "0 shows everything. Raise it to strip short, noisy matches the "
            "way a larger word size would when computing a dotplot."
        )
        self._filter_len_label = QLabel("Min length:", self)

        self._identity_label = QLabel("Min identity:", self)
        self._min_len_label = QLabel("Min length:", self)
        self._seed_freq_label = QLabel("Repeat cutoff:", self)

        self._zoom_out_btn = QToolButton(self)
        self._zoom_out_btn.setText("−")
        self._zoom_out_btn.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        self._zoom_out_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._zoom_out_btn.setToolTip("Zoom out (or scroll the wheel)")

        self._zoom_in_btn = QToolButton(self)
        self._zoom_in_btn.setText("+")
        self._zoom_in_btn.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        self._zoom_in_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._zoom_in_btn.setToolTip("Zoom in (or scroll the wheel)")

        self._reset_btn = QPushButton("Fit", self)
        self._reset_btn.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        self._reset_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._reset_btn.setToolTip("Fit the whole plot in the view")

        self._zoom_select_btn = QPushButton("Select", self)
        self._zoom_select_btn.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        self._zoom_select_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._zoom_select_btn.setToolTip(
            "Zoom to the selected region (drag a rectangle on the plot first)"
        )

        self._pdf_btn = QPushButton("Save", self)
        self._pdf_btn.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        self._pdf_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._pdf_btn.setToolTip("Save the current plot as a PDF")

        self._cancel_btn = QToolButton(self)
        self._cancel_btn.setText("Cancel")
        self._cancel_btn.setAutoRaise(True)
        self._cancel_btn.hide()

        # FastGA on a large genome takes minutes, so it is never started just
        # because the tab was opened — the user must click this button.
        self._compute_btn = QToolButton(self)
        self._compute_btn.setText("Run FastGA")
        self._compute_btn.setAutoRaise(True)
        self._compute_btn.setToolTip(
            "Run FastGA to build a self-dotplot of this long sequence. "
            "This can take several minutes on large chromosomes."
        )
        self._compute_btn.hide()

        self._thread: QThread | None = None
        self._worker: _FastgaWorker | None = None
        self._pending_key = None
        # Last drag-box as 1-based inclusive (x1, x2, y1, y2), if any.
        self._selection_xy: tuple[int, int, int, int] | None = None
        # Status line without the appended x/y readout (alignment mode).
        self._status_base = ""
        self._status_plain = ""

        self._status = QLabel("", self)
        self._status.setTextFormat(Qt.TextFormat.RichText)
        self._status.setContentsMargins(0, 2, 6, 2)

        self._sel_row = QWidget(self)
        sel_layout = QHBoxLayout(self._sel_row)
        sel_layout.setContentsMargins(0, 0, 0, 0)
        sel_layout.setSpacing(4)
        self._sel_x = QLabel("", self._sel_row)
        self._sel_x.setStyleSheet("color: palette(mid);")
        self._view_x_btn = QPushButton("View", self._sel_row)
        self._view_x_btn.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        self._view_x_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._view_x_btn.setToolTip(
            "Open Annotation and highlight the selected x (horizontal) range"
        )
        self._sel_sep = QLabel("·", self._sel_row)
        self._sel_sep.setStyleSheet("color: palette(mid);")
        self._sel_y = QLabel("", self._sel_row)
        self._sel_y.setStyleSheet("color: palette(mid);")
        self._view_y_btn = QPushButton("View", self._sel_row)
        self._view_y_btn.setFixedHeight(ZOOM_CONTROL_HEIGHT)
        self._view_y_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._view_y_btn.setToolTip(
            "Open Annotation and highlight the selected y (vertical) range"
        )
        sel_layout.addWidget(self._sel_x)
        sel_layout.addWidget(self._view_x_btn)
        sel_layout.addWidget(self._sel_sep)
        sel_layout.addWidget(self._sel_y)
        sel_layout.addWidget(self._view_y_btn)
        self._sel_row.hide()

        bar = QHBoxLayout()
        bar.setContentsMargins(6, 4, 6, 4)
        bar.setSpacing(8)
        self._word_label = QLabel("Word size:", self)
        bar.addWidget(self._pair_label)
        bar.addWidget(self._pair_combo)
        bar.addWidget(self._filter_len_label)
        bar.addWidget(self._filter_len)
        bar.addWidget(self._region_start_label)
        bar.addWidget(self._region_start)
        bar.addWidget(self._region_size_label)
        bar.addWidget(self._region_size)
        bar.addWidget(self._word_label)
        bar.addWidget(self._word)
        bar.addWidget(self._reverse)
        bar.addWidget(self._line_width_label)
        bar.addWidget(self._line_width)
        bar.addWidget(self._identity_label)
        bar.addWidget(self._identity)
        bar.addWidget(self._min_len_label)
        bar.addWidget(self._min_len)
        bar.addWidget(self._seed_freq_label)
        bar.addWidget(self._seed_freq)
        bar.addWidget(self._zoom_out_btn)
        bar.addWidget(self._zoom_in_btn)
        bar.addWidget(self._reset_btn)
        bar.addWidget(self._zoom_select_btn)
        bar.addWidget(self._pdf_btn)
        # Selection readout stays after Save (not in the title row).
        bar.addWidget(self._sel_row)
        bar.addWidget(self._compute_btn)
        bar.addWidget(self._cancel_btn)
        bar.addStretch(1)

        self._controls_bar = QWidget(self)
        self._controls_bar.setLayout(bar)

        self._title = QLabel("", self)
        self._title.setContentsMargins(_RULER_LEFT, 2, 6, 2)
        title_font = QFont(self.font())
        title_font.setBold(True)
        self._title.setFont(title_font)
        self._title.hide()

        self._canvas = DotplotCanvas(self)
        self._canvas.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        title_row = QHBoxLayout()
        title_row.setContentsMargins(0, 0, 8, 0)
        title_row.setSpacing(8)
        title_row.addWidget(self._title, 0, Qt.AlignmentFlag.AlignVCenter)
        title_row.addWidget(self._status, 0, Qt.AlignmentFlag.AlignVCenter)
        title_row.addStretch(1)

        layout.addWidget(self._controls_bar)
        layout.addLayout(title_row)
        layout.addWidget(self._canvas, stretch=1)

        self._word.valueChanged.connect(self._recompute)
        self._reverse.toggled.connect(self._recompute)
        self._line_width.valueChanged.connect(self._apply_line_width)
        self._identity.valueChanged.connect(self._recompute)
        self._min_len.valueChanged.connect(self._recompute)
        self._seed_freq.valueChanged.connect(self._recompute)
        self._region_start.valueChanged.connect(self._recompute)
        self._region_size.valueChanged.connect(self._recompute)
        self._pair_combo.currentIndexChanged.connect(self._show_selected_pair)
        self._filter_len.valueChanged.connect(self._recompute)
        self._reset_btn.clicked.connect(self._canvas.reset_zoom)
        self._zoom_select_btn.clicked.connect(self._zoom_to_selection)
        self._view_x_btn.clicked.connect(lambda: self._view_axis_in_annotation("x"))
        self._view_y_btn.clicked.connect(lambda: self._view_axis_in_annotation("y"))
        self._zoom_in_btn.clicked.connect(lambda: self._canvas.zoom_by(_ZOOM_STEP))
        self._zoom_out_btn.clicked.connect(lambda: self._canvas.zoom_by(1 / _ZOOM_STEP))
        self._pdf_btn.clicked.connect(self._save_pdf)
        self._cancel_btn.clicked.connect(self._cancel_running)
        self._compute_btn.clicked.connect(self._start_requested_compute)
        self._sync_engine_controls()
        self._canvas.positionChanged.connect(self._on_position)
        self._canvas.rectSelected.connect(self._on_rect_selected)

    def set_path(self, path: Path) -> None:
        resolved = path.expanduser().resolve()
        if resolved == self._path:
            return
        # Abandon any FastGA run for the previous file.
        self._cancel_running()
        self._path = resolved
        self._loaded_for = None
        self._mode = "fasta"
        self._alignments = None
        self._active_sequence_id = ""
        self._active_sequence = ""
        self._too_long_prompt_active = False
        self._canvas.set_axes([], [])
        self._sequence_id = ""
        self._canvas.clear_plot()
        self._canvas.set_title("")
        self._title.clear()
        self._title.hide()
        self.clear_selection()
        self._set_status_text("")
        self._compute_btn.hide()

    def set_active_sequence(self, sequence_id: str, sequence: str) -> None:
        """Use this record for Dotplot (Self) on multi-sequence FASTA files."""
        sid = str(sequence_id or "").strip()
        seq = str(sequence or "")
        if sid == self._active_sequence_id and seq == self._active_sequence:
            return
        self._active_sequence_id = sid
        self._active_sequence = seq
        self._full_length = len(seq)
        self._loaded_for = None
        self._too_long_prompt_active = False

    def _is_large(self) -> bool:
        if self._active_sequence:
            return len(self._active_sequence) > FASTGA_PREFERRED_ABOVE_BYTES
        try:
            return self._path.stat().st_size > FASTGA_PREFERRED_ABOVE_BYTES
        except OSError:
            return False

    def load_alignment(self, path: Path) -> None:
        """Display a PAF (or .1aln, converted) instead of computing anything."""
        self._mode = "alignment"
        self._path = path.expanduser().resolve()
        self._loaded_for = None
        self._set_status_text("Reading alignments…")
        try:
            alignments, source = load_alignment_file(self._path)
        except PafError as exc:
            self._alignments = None
            self._sync_engine_controls()
            self._fail(str(exc))
            return
        except Exception as exc:  # noqa: BLE001
            self._alignments = None
            self._sync_engine_controls()
            self._fail(f"Could not read {self._path.name}: {exc}")
            return

        self._alignments = alignments
        pairs = alignments.pairs()
        self._pair_combo.blockSignals(True)
        self._pair_combo.clear()
        if len(pairs) > 1:
            self._pair_combo.addItem("All sequences", ("__all__", "__all__"))
        for query, target in pairs:
            label = query if query == target else f"{query}  vs  {target}"
            self._pair_combo.addItem(label, (query, target))
        self._pair_combo.setCurrentIndex(0)
        self._pair_combo.blockSignals(False)
        longest = max(
            list(alignments.query_lengths.values())
            + list(alignments.target_lengths.values())
            + [1]
        )
        for box in (self._region_start, self._region_size):
            box.blockSignals(True)
        self._region_start.setMaximum(max(1, longest))
        self._region_start.setValue(1)
        self._region_size.setMaximum(max(1_000, longest))
        self._region_size.setValue(max(1_000, longest))  # whole sequence first
        for box in (self._region_start, self._region_size):
            box.blockSignals(False)

        self._sync_engine_controls()
        if source != self._path:
            self.statusMessage.emit(f"Converted {self._path.name} to {source.name}")
        self._show_selected_pair()

    def _show_selected_pair(self) -> None:
        if self._alignments is None:
            return
        data = self._pair_combo.currentData()
        try:
            min_len = int(self._filter_len.value())
            if data is None or data[0] == "__all__":
                result, x_entries, y_entries = result_for_all(
                    self._alignments, min_length=min_len
                )
            else:
                full = max(
                    self._alignments.query_lengths.get(data[0], 0),
                    self._alignments.target_lengths.get(data[1], 0),
                    1,
                )
                start = int(self._region_start.value())
                size = int(self._region_size.value())
                window = None
                if start > 1 or size < full:
                    window = (start, min(start + size - 1, full))
                result = result_for_pair(
                    self._alignments,
                    data[0],
                    data[1],
                    window=window,
                    min_length=min_len,
                )
                x_entries, y_entries = [], []
        except Exception as exc:  # noqa: BLE001
            self._fail(f"Could not build the plot: {exc}")
            return
        self._canvas.set_axes(x_entries, y_entries)
        self._sync_engine_controls()
        self._apply_result(("alignment", self._path, str(data)), result)

    def _use_fastga(self) -> bool:
        """Large sequences use FastGA (opt-in via Run FastGA), not the k-mer window."""
        return self._is_large() and fastga_available()

    def _use_region(self) -> bool:
        # Instant From/Window plotting is disabled for long sequences; users
        # are directed to FastGA instead.
        return False

    def _sync_engine_controls(self) -> None:
        if self._mode == "alignment":
            # Nothing is computed from an alignment file, so the compute
            # controls have nothing to act on.
            for widget in (
                self._word_label, self._word, self._reverse,
                self._identity_label, self._identity,
                self._min_len_label, self._min_len,
                self._seed_freq_label, self._seed_freq,
            ):
                widget.setVisible(False)
            multi = self._alignments is not None and not self._alignments.is_single_pair()
            self._pair_label.setVisible(multi)
            self._pair_combo.setVisible(multi)
            self._filter_len_label.setVisible(True)
            self._filter_len.setVisible(True)
            # FastGA / PAF previews: no From/Window; View select stays hidden
            # (x/y is appended to the status line for self / pair compares).
            for widget in (
                self._region_start_label, self._region_start,
                self._region_size_label, self._region_size,
            ):
                widget.setVisible(False)
            self._sel_row.hide()
            return

        self._pair_label.setVisible(False)
        self._pair_combo.setVisible(False)
        self._filter_len_label.setVisible(False)
        self._filter_len.setVisible(False)
        # Status-only note for long sequences: hide Line width / zoom / Save.
        hide_bar = self._too_long_prompt_active
        self._controls_bar.setVisible(not hide_bar)
        if hide_bar:
            self._compute_btn.hide()
            self._title.hide()
            return
        for widget in (self._word_label, self._word, self._reverse):
            widget.setVisible(not self._is_large())
        for widget in (
            self._region_start_label,
            self._region_start,
            self._region_size_label,
            self._region_size,
        ):
            widget.setVisible(False)
        for widget in (
            self._identity_label,
            self._identity,
            self._min_len_label,
            self._min_len,
            self._seed_freq_label,
            self._seed_freq,
        ):
            widget.setVisible(False)

    def _current_key(self):
        active = self._active_sequence_id or ""
        if self._use_fastga():
            return (
                self._path,
                active,
                "fastga",
                round(float(self._identity.value()), 3),
                int(self._min_len.value()),
                int(self._seed_freq.value()),
            )
        return (
            self._path,
            active,
            "kmer",
            int(self._word.value()),
            self._reverse.isChecked(),
        )

    def ensure_loaded(self) -> None:
        """Compute on first view only, never when a file is merely selected.

        Short sequences use the instant k-mer engine. Sequences longer than
        ~2 Mb need FastGA: show a prompt (never auto-start a multi-minute run).
        """
        if self._mode == "alignment":
            return  # already loaded from file; nothing to compute
        self._sync_engine_controls()
        key = self._current_key()
        if self._loaded_for == key or self._thread is not None:
            return
        if self._is_large():
            if self._use_fastga():
                self._offer_fastga(key)
            else:
                self._offer_too_long_no_fastga()
            return
        self._compute_kmer(key)

    def _too_long_status_message(self, *, fastga_available_now: bool) -> str:
        del fastga_available_now  # same copy whether FastGA is installed
        return (
            "This sequence is too long for instant self-dotplot "
            "(limit 2.0 Mbp of sequence length)."
        )

    def _offer_too_long_no_fastga(self) -> None:
        """Large sequence but FastGA is not installed."""
        self._pending_key = None
        self._canvas.clear_plot()
        self._clear_selection_chrome()
        self._compute_btn.hide()
        # Remember so ensure_loaded does not re-fire the same prompt every paint.
        self._loaded_for = ("too_long", self._path, self._active_sequence_id or "")
        self._too_long_prompt_active = True
        self._sync_engine_controls()
        message = self._too_long_status_message(fastga_available_now=False)
        self._set_status_text(message, tooltip="")
        self.statusMessage.emit(message)

    def _offer_fastga(self, key) -> None:
        """Status-only note for long sequences (no in-tab Run FastGA UI)."""
        self._pending_key = None
        self._loaded_for = key  # prompt shown; do not re-fire on every ensure_loaded
        self._too_long_prompt_active = True
        self._canvas.clear_plot()
        self._clear_selection_chrome()
        self._compute_btn.hide()
        self._sync_engine_controls()
        message = self._too_long_status_message(fastga_available_now=True)
        self._set_status_text(message, tooltip="")

    def _start_requested_compute(self) -> None:
        self._compute_btn.hide()
        key = self._current_key()
        if self._thread is None:
            self._compute_fastga(key)

    def _apply_line_width(self, step: int) -> None:
        self._canvas.set_line_width_step(int(step))

    def _recompute(self) -> None:
        if self._mode == "alignment":
            if self._alignments is not None:
                self._show_selected_pair()
            return
        if self._loaded_for is None:
            return
        # Long-sequence note is informational only — do not relaunch FastGA from here.
        if self._is_large():
            return
        self.ensure_loaded()

    def _compute_region(self, key) -> None:
        """Plot one window of a large sequence with the instant engine."""
        self._set_status_text("Reading region…")
        start = int(self._region_start.value())
        size = int(self._region_size.value())
        try:
            if self._full_length <= 0:
                _name, self._full_length = sequence_length_of_fasta(self._path)
                self._region_start.setMaximum(max(1, self._full_length))
            end = min(start + size - 1, self._full_length or (start + size - 1))
            result = compute_self_dotplot_for_region(
                self._path,
                start,
                end,
                word_size=int(self._word.value()),
                include_reverse=self._reverse.isChecked(),
                full_length=self._full_length,
            )
        except DotplotError as exc:
            self._fail(str(exc))
            return
        except Exception as exc:  # noqa: BLE001
            self._fail(f"Could not plot region: {exc}")
            return
        self._apply_result(key, result)

    # -- built-in k-mer engine (fast, synchronous) -------------------------

    def _compute_kmer(self, key) -> None:
        self._set_status_text("Computing…")
        try:
            kwargs: dict = {
                "word_size": int(self._word.value()),
                "include_reverse": self._reverse.isChecked(),
            }
            if self._active_sequence:
                kwargs["sequence"] = self._active_sequence
                kwargs["sequence_id"] = self._active_sequence_id
            result = compute_self_dotplot_for_fasta(self._path, **kwargs)
        except DotplotError as exc:
            message = str(exc)
            # The built-in engine refuses very large sequences. If FastGA is not
            # installed, say so here -- that is exactly the case it solves.
            if not fastga_available() and "limited to" in message:
                message += f"\n\n{install_hint()}"
            self._fail(message)
            return
        except Exception as exc:  # noqa: BLE001
            self._fail(f"Could not build dotplot: {exc}")
            return
        self._apply_result(key, result)

    # -- FastGA engine (slow, threaded) ------------------------------------

    def _compute_fastga(self, key) -> None:
        self._pending_key = key
        self._set_status_text("Starting FastGA…")
        self._cancel_btn.show()

        thread = QThread(self)
        worker = _FastgaWorker(
            self._path,
            threads=8,
            min_identity=float(self._identity.value()),
            min_length=int(self._min_len.value()),
            seed_freq=int(self._seed_freq.value()),
        )
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.progress.connect(self._set_status_text)
        worker.finished.connect(self._on_fastga_finished)
        worker.failed.connect(self._on_fastga_failed)
        # Tear the thread down from either outcome, and only delete the worker
        # once its thread has actually stopped.
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(self._on_thread_finished)
        self._thread = thread
        self._worker = worker
        thread.start()

    def _cancel_running(self) -> None:
        if self._worker is not None:
            self._worker.cancel()
            self._set_status_text("Cancelling…")

    @Slot()
    def _on_thread_finished(self) -> None:
        thread, self._thread = self._thread, None
        self._worker = None
        self._cancel_btn.hide()
        if thread is not None:
            thread.deleteLater()

    @Slot(object)
    def _on_fastga_finished(self, result) -> None:
        self._apply_result(getattr(self, "_pending_key", None), result)

    @Slot(str)
    def _on_fastga_failed(self, message: str) -> None:
        self._fail(message)

    # -- shared -----------------------------------------------------------

    @staticmethod
    def _condense_error(message: str) -> str:
        """One informative line from a multi-line tool error.

        The first line is often just "FastGA failed (exit 1)."; the actual
        reason is on a later line. Prefer whichever line names the tool or
        explains the problem, so the status bar is not merely "it failed".
        """
        lines = [ln.strip() for ln in message.splitlines() if ln.strip()]
        if not lines:
            return "Failed."
        for line in lines[1:]:
            if ":" in line or line.lower().startswith(("cannot", "could", "no ")):
                return line
        return lines[0]

    def _fail(self, message: str) -> None:
        self._canvas.clear_plot()
        self._loaded_for = None
        self._clear_selection_chrome()
        self._set_status_text(self._condense_error(message), tooltip=message)
        self.statusMessage.emit(message)

    def _apply_result(self, key, result: DotplotResult) -> None:
        previous_length = self._sequence_length
        had_plot = self._loaded_for is not None
        self._loaded_for = key
        self._too_long_prompt_active = False
        self._sync_engine_controls()
        self._sequence_id = result.sequence_id
        self._sequence_length = result.sequence_length
        name = result.sequence_id or self._path.stem
        self._title.setText(name)
        self._title.setVisible(bool(name))
        self._canvas.set_title(name)
        # Keep pan/zoom when recomputing the same-length sequence (word size /
        # reverse complement). New files / region sizes still fit to view.
        preserve = had_plot and previous_length == result.sequence_length
        self._canvas.show_result(result, preserve_view=preserve)
        if result.word_size:
            engine = f"word {result.word_size}"
        elif self._mode == "alignment":
            # The file may come from FastGA or any other PAF producer -- do
            # not claim to know which.
            engine = self._path.suffix.lstrip(".").upper() or "alignments"
        else:
            engine = "FastGA"
        offset = int(getattr(result, "coordinate_offset", 0) or 0)
        full = int(getattr(result, "full_sequence_length", 0) or 0)
        # Length / word size stay available for PDF captions; the title row
        # only shows the segment count in parentheses next to the name.
        if offset or (full and full > result.sequence_length):
            if full:
                where_plain = (
                    f"{offset + 1:,}-{offset + result.sequence_length:,}"
                    f" of {full:,} bp"
                )
            else:
                where_plain = f"{offset + 1:,}-{offset + result.sequence_length:,}"
        else:
            where_plain = f"{result.sequence_length:,} bp"
        engine_plain = (
            f"word {result.word_size}" if result.word_size else engine
        )
        suffix = f" — {result.note}" if result.note else ""
        # New plot data clears any previous drag box; chrome follows.
        self._clear_selection_chrome()
        self._set_status_text(
            f"(<b>{result.segment_count:,}</b> hits){suffix}",
            plain=(
                f"{where_plain} · {engine_plain} · "
                f"{result.segment_count:,} hits{suffix}"
            ),
        )

    def suggest_default_word_size(self) -> None:
        try:
            size = self._path.stat().st_size
        except OSError:
            return
        self._word.blockSignals(True)
        self._word.setValue(suggest_word_size(size))
        self._word.blockSignals(False)

    def show_selection(self, start: int, end: int) -> None:
        """Mirror a 1-D Annotation selection onto the plot as an x/y box."""
        lo, hi = (int(start), int(end)) if start <= end else (int(end), int(start))
        self._canvas.set_range(lo, hi)
        # Same range on both axes — show x and y explicitly (not "Selected N bp").
        self._show_selection_status(lo, hi, lo, hi)

    def clear_selection(self) -> None:
        self._canvas.clear_selection()
        self._clear_selection_chrome()

    def _clear_selection_chrome(self) -> None:
        self._selection_xy = None
        self._sel_x.clear()
        self._sel_y.clear()
        self._sel_row.hide()
        if self._mode == "alignment" and self._status_base:
            self._status.setTextFormat(Qt.TextFormat.RichText)
            self._status.setText(self._status_base)

    def _set_status_text(
        self,
        text: str,
        *,
        tooltip: str = "",
        plain: str | None = None,
    ) -> None:
        """Show status beside the sequence name (numbers may be bold HTML)."""
        self._status_base = text
        self._status_plain = plain if plain is not None else text
        self._status.show()
        # Transient messages are plain; meta stats use RichText with <b> numbers.
        if "<" in text and ">" in text:
            self._status.setTextFormat(Qt.TextFormat.RichText)
        else:
            self._status.setTextFormat(Qt.TextFormat.PlainText)
        self._status.setText(text)
        self._status.setToolTip(tooltip or self._status_plain)

    def _status_caption(self) -> str:
        """Plain-text status for PDF captions / tooltips."""
        base = self._status_plain or self._status_base
        if self._selection_xy is not None:
            x1, x2, y1, y2 = self._selection_xy
            return f"{base} · x {x1:,}-{x2:,} · y {y1:,}-{y2:,}"
        return base

    def _alignment_shows_xy_readout(self) -> bool:
        """True for self / two-sequence compares (not 'All sequences')."""
        if self._mode != "alignment" or self._alignments is None:
            return False
        if self._alignments.is_single_pair():
            return True
        data = self._pair_combo.currentData()
        return bool(data) and data[0] != "__all__"

    def _show_selection_status(self, x1: int, x2: int, y1: int, y2: int) -> None:
        self._selection_xy = (int(x1), int(x2), int(y1), int(y2))
        if self._mode == "alignment":
            # No View select buttons; append coords after the segment count.
            self._sel_row.hide()
            if self._alignment_shows_xy_readout():
                base = self._status_base
                self._status.setTextFormat(Qt.TextFormat.RichText)
                self._status.setText(
                    f"{base} · x <b>{x1:,}</b>–<b>{x2:,}</b> · "
                    f"y <b>{y1:,}</b>–<b>{y2:,}</b>"
                )
            return
        self._sel_x.setText(f"x {x1:,}-{x2:,}")
        self._sel_y.setText(f"y {y1:,}-{y2:,}")
        self._sel_row.show()

    def _zoom_to_selection(self) -> None:
        if self._canvas.zoom_to_selection():
            return
        self._set_status_text("Select a region on the plot first (drag a rectangle).")

    def _view_axis_in_annotation(self, axis: str) -> None:
        if self._mode == "alignment":
            return
        coords = self._selection_xy or self._canvas.selection_sequence_coords()
        if coords is None:
            self._set_status_text("Select a region on the plot first (drag a rectangle).")
            return
        x1, x2, y1, y2 = coords
        if axis == "y":
            start, end = (y1, y2) if y1 <= y2 else (y2, y1)
        else:
            start, end = (x1, x2) if x1 <= x2 else (x2, x1)
        self.viewInAnnotationRequested.emit(int(start), int(end))

    def _save_pdf(self) -> None:
        from PySide6.QtWidgets import QFileDialog, QMessageBox

        if self._sequence_length <= 0:
            QMessageBox.information(self, "Save PDF", "There is no plot to save yet.")
            return
        default = self._path.parent / f"{self._path.stem}_dotplot.pdf"
        chosen, _ = QFileDialog.getSaveFileName(
            self, "Save dotplot as PDF", str(default), "PDF (*.pdf)"
        )
        if not chosen:
            return
        out = Path(chosen)
        if out.suffix.lower() != ".pdf":
            out = out.with_suffix(".pdf")
        try:
            self._canvas.render_to_pdf(
                out,
                title=self._title.text() or self._path.stem,
                # The status line already says range, word size and segment
                # count -- exactly what a figure caption needs.
                subtitle=self._status_caption(),
            )
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Save PDF", f"Could not write {out.name}:\n{exc}")
            return
        self.statusMessage.emit(f"Saved {out.name}")

    def sequence_id(self) -> str:
        return self._sequence_id

    def _on_rect_selected(self, x1: int, x2: int, y1: int, y2: int) -> None:
        self._show_selection_status(x1, x2, y1, y2)
        if self._mode != "alignment":
            self.rangeSelected.emit(x1, x2, y1, y2)

    def _on_position(self, x: int, y: int) -> None:
        self.statusMessage.emit(f"x: {x:,}  y: {y:,}")
