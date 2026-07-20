"""Interactive self-comparison dotplot ("Dotplot (Self)" preview tab).

Segments are computed once by ``teviewer.backends.dotplot`` and then only
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

from PySide6.QtCore import QObject, QPoint, QPointF, QRectF, Qt, QThread, Signal, Slot
from PySide6.QtGui import QColor, QFont, QPainter, QPen
from PySide6.QtWidgets import (
    QCheckBox,
    QDoubleSpinBox,
    QGraphicsItem,
    QGraphicsScene,
    QGraphicsView,
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QSpinBox,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from teviewer.backends.dotplot import (
    DEFAULT_WORD_SIZE,
    MAX_WORD_SIZE,
    MIN_WORD_SIZE,
    STRAND_REVERSE,
    DotplotError,
    DotplotResult,
    compute_self_dotplot_for_fasta,
    suggest_word_size,
)
from teviewer.backends.fastga_dotplot import (
    DEFAULT_MIN_IDENTITY,
    DEFAULT_MIN_LENGTH,
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
    ) -> None:
        super().__init__()
        self._path = path
        self._threads = threads
        self._min_identity = min_identity
        self._min_length = min_length
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

_FORWARD_COLOR = QColor(40, 90, 180)
_REVERSE_COLOR = QColor(200, 60, 40)
_DIAGONAL_COLOR = QColor(230, 120, 40)

_MIN_SCALE = 0.02
_MAX_SCALE = 400.0

_RULER_TOP = 18
# Wide enough for the tick labels plus the rotated sequence name beside them.
_RULER_LEFT = 66
_Y_LABEL_STRIP = 16
_RULER_BG = QColor(250, 250, 250)
_RULER_FG = QColor(110, 110, 110)
_RULER_LINE = QColor(200, 200, 200)

_SELECT_FILL = QColor(70, 130, 220, 55)
_SELECT_EDGE = QColor(60, 110, 200, 190)

_CLICK_SLOP = 4


def _nice_tick_step(span: float, target_ticks: int = 8) -> int:
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
    if value >= 1_000_000 and value % 1_000_000 == 0:
        return f"{value // 1_000_000}M"
    if value >= 1_000 and value % 1_000 == 0:
        return f"{value // 1_000}k"
    return f"{value:,}"


class _SegmentsItem(QGraphicsItem):
    """Paints every dotplot segment in one pass."""

    def __init__(self, result: DotplotResult) -> None:
        super().__init__()
        self._result = result
        self._n = max(result.sequence_length, 1)
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemUsesExtendedStyleOption, True)

    def boundingRect(self) -> QRectF:  # noqa: N802
        return QRectF(0, 0, float(self._n), float(self._n))

    def paint(self, painter: QPainter, option, widget=None) -> None:  # noqa: N802
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        scale = painter.transform().m11() or 1.0
        min_length = 0.0 if scale >= 1.0 else (1.0 / scale)
        visible = option.exposedRect if option is not None else self.boundingRect()

        pen = QPen(_DIAGONAL_COLOR)
        pen.setCosmetic(True)
        painter.setPen(pen)
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

        pen = QPen(_FORWARD_COLOR)
        pen.setCosmetic(True)
        painter.setPen(pen)
        for x1, y1, x2, y2 in forward:
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
            painter.drawLine(QPointF(y1, x1), QPointF(y2, x2))  # self-symmetry

        pen = QPen(_REVERSE_COLOR)
        pen.setCosmetic(True)
        painter.setPen(pen)
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

        visible = view.visible_scene_rect()
        length = view.plot_length()
        span = visible.width() if self._horizontal else visible.height()
        step = _nice_tick_step(max(span, 1.0))
        lo = visible.left() if self._horizontal else visible.top()
        hi = visible.right() if self._horizontal else visible.bottom()

        painter.setPen(QPen(_RULER_FG))
        tick = int(lo // step) * step
        while tick <= hi:
            if 0 <= tick <= length:
                pos = view.scene_to_ruler(tick, self._horizontal)
                if self._horizontal:
                    if 0 <= pos <= self.width():
                        painter.drawLine(int(pos), self.height() - 4, int(pos), self.height())
                        painter.drawText(int(pos) + 2, self.height() - 6, _format_bp(int(tick)))
                else:
                    if 0 <= pos <= self.height():
                        painter.drawLine(self.width() - 4, int(pos), self.width(), int(pos))
                        painter.drawText(_Y_LABEL_STRIP, int(pos) - 2, _format_bp(int(tick)))
            tick += step

        cursor = view.cursor_scene_pos()
        if cursor is not None:
            painter.setPen(QPen(QColor(40, 90, 180)))
            value = cursor.x() if self._horizontal else cursor.y()
            pos = view.scene_to_ruler(value, self._horizontal)
            if self._horizontal and 0 <= pos <= self.width():
                painter.drawText(int(pos) + 2, self.height() - 6, f"{int(value):,}")
            elif not self._horizontal and 0 <= pos <= self.height():
                painter.drawText(_Y_LABEL_STRIP, int(pos) - 2, f"{int(value):,}")

        # Sequence name down the y axis, as Geneious labels its dotplot.
        if not self._horizontal and view.plot_title():
            painter.setPen(QPen(_RULER_FG))
            painter.save()
            painter.translate(11, self.height() / 2.0)
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
        self.setBackgroundBrush(QColor(255, 255, 255))
        self._item: _SegmentsItem | None = None
        self._length = 0
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

    def _refresh_rulers(self) -> None:
        self._top_ruler.update()
        self._left_ruler.update()

    def _position_rulers(self) -> None:
        vp = self.viewport().geometry()
        self._top_ruler.setGeometry(vp.x(), 0, vp.width(), _RULER_TOP)
        self._left_ruler.setGeometry(0, vp.y(), _RULER_LEFT, vp.height())

    # ------------------------------------------------------------------
    # content
    # ------------------------------------------------------------------

    def show_result(self, result: DotplotResult) -> None:
        scene = self.scene()
        scene.clear()
        self._item = _SegmentsItem(result)
        self._length = result.sequence_length
        self._selection = None
        scene.addItem(self._item)
        scene.setSceneRect(0, 0, result.sequence_length, result.sequence_length)
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

        Ruler space is reserved by setViewportMargins(), so the viewport *is*
        the plot area and this is a straight fit with no padding arithmetic.
        """
        if not self._length:
            return
        self.resetTransform()
        # The viewport is now exactly the plot area, so fit the sequence square
        # into it directly -- no padding tricks.
        n = float(self._length)
        self.fitInView(QRectF(0, 0, n, n), Qt.AspectRatioMode.KeepAspectRatio)
        self.viewport().update()
        self._refresh_rulers()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._position_rulers()
        self._fit_if_needed()
        self._refresh_rulers()

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
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
        lo = max(0, lo - 1)
        hi = min(self._length, hi)
        if hi <= lo:
            return
        self._selection = QRectF(lo, lo, hi - lo, hi - lo)
        self.viewport().update()

    def clear_selection(self) -> None:
        self._selection = None
        self._drag_start = self._drag_now = None
        self.viewport().update()

    # ------------------------------------------------------------------
    # painting
    # ------------------------------------------------------------------

    def drawForeground(self, painter: QPainter, rect: QRectF) -> None:  # noqa: N802
        """Scene-space overlays only: selection, plot border, crosshair lines.

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

        border = QPen(QColor(120, 120, 120))
        border.setCosmetic(True)
        painter.setPen(border)
        painter.setBrush(Qt.BrushStyle.NoBrush)
        painter.drawRect(QRectF(0, 0, self._length, self._length))

        if self._cursor_pos is not None:
            cross = QPen(QColor(90, 140, 220, 150))
            cross.setCosmetic(True)
            painter.setPen(cross)
            cx = self._cursor_pos.x()
            cy = self._cursor_pos.y()
            painter.drawLine(QPointF(cx, 0), QPointF(cx, self._length))
            painter.drawLine(QPointF(0, cy), QPointF(self._length, cy))

    def wheelEvent(self, event) -> None:  # noqa: N802
        if self._item is None:
            return
        factor = 1.25 if event.angleDelta().y() > 0 else 1 / 1.25
        target = self.transform().m11() * factor
        if target < _MIN_SCALE or target > _MAX_SCALE:
            return
        self.scale(factor, factor)
        self.viewport().update()

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
            self.positionChanged.emit(x, y)

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
                    int(rect.left()) + 1,
                    int(rect.right()),
                    int(rect.top()) + 1,
                    int(rect.bottom()),
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

    def __init__(self, path: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("DotplotPanel")
        self._path = path.expanduser().resolve()
        self._loaded_for: tuple[Path, int, bool] | None = None
        self._sequence_id = ""
        self._sequence_length = 0

        self._word = QSpinBox(self)
        self._word.setRange(MIN_WORD_SIZE, MAX_WORD_SIZE)
        self._word.setValue(DEFAULT_WORD_SIZE)
        self._word.setToolTip(
            "Length of the exact match used to place a dot. Larger = stricter "
            "and faster; smaller = more sensitive and noisier."
        )
        self._word.setFixedWidth(70)

        self._reverse = QCheckBox("Reverse complement", self)
        self._reverse.setChecked(True)
        self._reverse.setToolTip(
            "Show inverted repeats (terminal inverted repeats, hairpins) as "
            "anti-diagonal lines."
        )

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

        self._identity_label = QLabel("Min identity:", self)
        self._min_len_label = QLabel("Min length:", self)

        self._reset_btn = QToolButton(self)
        self._reset_btn.setText("Reset zoom")
        self._reset_btn.setAutoRaise(True)

        self._cancel_btn = QToolButton(self)
        self._cancel_btn.setText("Cancel")
        self._cancel_btn.setAutoRaise(True)
        self._cancel_btn.hide()

        # FastGA on a large genome takes minutes (measured ~13 min on a 648 Mb
        # scaffold), so it is never started just because the tab was opened.
        self._compute_btn = QToolButton(self)
        self._compute_btn.setText("Compute dotplot")
        self._compute_btn.setAutoRaise(True)
        self._compute_btn.hide()

        self._thread: QThread | None = None
        self._worker: _FastgaWorker | None = None
        self._pending_key = None

        self._status = QLabel("", self)
        self._status.setStyleSheet("color: palette(mid);")

        bar = QHBoxLayout()
        bar.setContentsMargins(6, 4, 6, 4)
        bar.setSpacing(8)
        self._word_label = QLabel("Word size:", self)
        bar.addWidget(self._word_label)
        bar.addWidget(self._word)
        bar.addWidget(self._reverse)
        bar.addWidget(self._identity_label)
        bar.addWidget(self._identity)
        bar.addWidget(self._min_len_label)
        bar.addWidget(self._min_len)
        bar.addWidget(self._reset_btn)
        bar.addWidget(self._compute_btn)
        bar.addWidget(self._cancel_btn)
        bar.addStretch(1)
        bar.addWidget(self._status)

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
        layout.addLayout(bar)
        layout.addWidget(self._title)
        layout.addWidget(self._canvas, stretch=1)

        self._word.valueChanged.connect(self._recompute)
        self._reverse.toggled.connect(self._recompute)
        self._identity.valueChanged.connect(self._recompute)
        self._min_len.valueChanged.connect(self._recompute)
        self._reset_btn.clicked.connect(self._canvas.reset_zoom)
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
        self._sequence_id = ""
        self._canvas.clear_plot()
        self._canvas.set_title("")
        self._title.clear()
        self._title.hide()
        self._status.clear()
        self._compute_btn.hide()

    def _use_fastga(self) -> bool:
        """Large sequences go to FastGA; small ones stay on the instant engine."""
        try:
            size = self._path.stat().st_size
        except OSError:
            return False
        return size > FASTGA_PREFERRED_ABOVE_BYTES and fastga_available()

    def _sync_engine_controls(self) -> None:
        fastga = self._use_fastga()
        for widget in (self._word_label, self._word, self._reverse):
            widget.setVisible(not fastga)
        for widget in (
            self._identity_label,
            self._identity,
            self._min_len_label,
            self._min_len,
        ):
            widget.setVisible(fastga)

    def _current_key(self):
        if self._use_fastga():
            return (
                self._path,
                "fastga",
                round(float(self._identity.value()), 3),
                int(self._min_len.value()),
            )
        return (
            self._path,
            "kmer",
            int(self._word.value()),
            self._reverse.isChecked(),
        )

    def ensure_loaded(self) -> None:
        """Compute on first view only, never when a file is merely selected.

        The small-sequence engine runs immediately (milliseconds). FastGA is
        offered behind a button instead: on a large genome it runs for minutes,
        which must never start merely because a tab was clicked.
        """
        self._sync_engine_controls()
        key = self._current_key()
        if self._loaded_for == key or self._thread is not None:
            return
        if self._use_fastga():
            self._offer_fastga(key)
        else:
            self._compute_kmer(key)

    def _offer_fastga(self, key) -> None:
        """Show what it will cost and wait for the user to ask for it."""
        self._pending_key = key
        self._canvas.clear_plot()
        try:
            mb = self._path.stat().st_size / 1024 / 1024
        except OSError:
            mb = 0
        self._compute_btn.show()
        self._status.setText(
            f"{mb:,.0f} MB — too large for the instant engine. "
            "FastGA can plot it, but takes minutes on a genome this size."
        )

    def _start_requested_compute(self) -> None:
        self._compute_btn.hide()
        key = self._current_key()
        if self._thread is None:
            self._compute_fastga(key)

    def _recompute(self) -> None:
        if self._loaded_for is None:
            return
        # Changing a FastGA setting must not silently relaunch a long run.
        if self._use_fastga():
            self._loaded_for = None
            self._offer_fastga(self._current_key())
            return
        self.ensure_loaded()

    # -- built-in k-mer engine (fast, synchronous) -------------------------

    def _compute_kmer(self, key) -> None:
        self._status.setText("Computing…")
        try:
            result = compute_self_dotplot_for_fasta(
                self._path,
                word_size=int(self._word.value()),
                include_reverse=self._reverse.isChecked(),
            )
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
        self._status.setText("Starting FastGA…")
        self._cancel_btn.show()

        thread = QThread(self)
        worker = _FastgaWorker(
            self._path,
            threads=8,
            min_identity=float(self._identity.value()),
            min_length=int(self._min_len.value()),
        )
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.progress.connect(self._status.setText)
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
            self._status.setText("Cancelling…")

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
        self._status.setText(self._condense_error(message))
        # Full text on hover: the condensed line cannot carry everything.
        self._status.setToolTip(message)
        self.statusMessage.emit(message)

    def _apply_result(self, key, result: DotplotResult) -> None:
        self._loaded_for = key
        self._sequence_id = result.sequence_id
        self._sequence_length = result.sequence_length
        name = result.sequence_id or self._path.stem
        self._title.setText(name)
        self._title.setVisible(bool(name))
        self._canvas.set_title(name)
        self._canvas.show_result(result)
        engine = (
            f"word {result.word_size}" if result.word_size else "FastGA"
        )
        suffix = f" — {result.note}" if result.note else ""
        self._status.setText(
            f"{result.sequence_length:,} bp · {engine} · "
            f"{result.segment_count:,} segments{suffix}"
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
        """Mirror a selection made in the Annotation tab onto the plot."""
        self._canvas.set_range(start, end)
        self._status.setText(
            f"Selected {start:,}-{end:,} ({end - start + 1:,} bp)"
        )

    def clear_selection(self) -> None:
        self._canvas.clear_selection()

    def sequence_id(self) -> str:
        return self._sequence_id

    def _on_rect_selected(self, x1: int, x2: int, y1: int, y2: int) -> None:
        self._status.setText(
            f"x {x1:,}-{x2:,}  ·  y {y1:,}-{y2:,}"
        )
        self.rangeSelected.emit(x1, x2, y1, y2)

    def _on_position(self, x: int, y: int) -> None:
        self.statusMessage.emit(f"x: {x:,}  y: {y:,}")
