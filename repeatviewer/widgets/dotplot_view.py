"""Interactive self-comparison dotplot ("Dotplot (Self)" preview tab).

Segments are computed once by ``repeatviewer.backends.dotplot`` and then only
redrawn: zooming and panning never recompute. A new pass happens only when the
word size or the reverse-complement setting changes.

All segments are painted by a single QGraphicsItem rather than one item per
segment. Tens of thousands of QGraphicsLineItems would each carry their own
bounding-box bookkeeping and make panning stutter; one item painting a list of
lines stays smooth. Segments shorter than a pixel at the current zoom are
drawn as points so a dense plot still shows where matches sit when zoomed out.

Rulers are drawn hard against the edges of the plot itself, not the viewport,
so they stay beside the data the way they do in Geneious.
"""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import QEvent, QObject, QPoint, QPointF, QRect, QRectF, Qt, QThread, Signal, Slot
from PySide6.QtGui import (
    QColor,
    QFont,
    QNativeGestureEvent,
    QPainter,
    QPalette,
    QPen,
    QStandardItem,
    QStandardItemModel,
    QTransform,
)
from PySide6.QtWidgets import (
    QGraphicsItem,
    QGraphicsScene,
    QGraphicsView,
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.dotplot_gallery import DotplotGallery
from repeatviewer.widgets.dotplot_rail import DotplotRail
from repeatviewer.widgets.preview_zoom import (
    PREVIEW_WHEEL_DELTA_PER_STEP,
    wheel_zoom_delta_from_event,
)

from repeatviewer.backends.dotplot import (
    DEFAULT_WORD_SIZE,
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
    MAX_ALL_SEQUENCES,
    AxisEntry,
    PafError,
    load_alignment_file,
    parse_paf_pair,
    result_for_all,
    result_for_pair,
)
from repeatviewer.backends.fastga_dotplot import (
    DEFAULT_SEED_FREQ,
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


class _PafLoadWorker(QObject):
    """Parse a PAF / convert .1aln off the UI thread.

    Whole-genome multi-FastGA PAFs are often hundreds of megabytes; reading
    them on the GUI thread freezes Preview until the file is fully scanned.
    """

    progress = Signal(str)
    finished = Signal(object, object, object)  # alignments, source Path, request Path
    failed = Signal(str)

    def __init__(self, path: Path) -> None:
        super().__init__()
        self._path = path
        self._cancelled = False

    def cancel(self) -> None:
        self._cancelled = True

    @Slot()
    def run(self) -> None:
        try:
            alignments, source = load_alignment_file(
                self._path,
                progress=self.progress.emit,
                is_cancelled=lambda: self._cancelled,
            )
        except PafError as exc:
            self.failed.emit(str(exc))
            return
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(f"Could not read {self._path.name}: {exc}")
            return
        # Include the requested path so a superseded load is ignored if the
        # user clicks another file before this worker finishes.
        self.finished.emit(alignments, source, self._path)


class _PafPairLoadWorker(QObject):
    """Re-scan a huge PAF for one Compare pair (FastGA selection-style)."""

    progress = Signal(str)
    finished = Signal(object, object, object)  # alignments, query, target
    failed = Signal(str)

    def __init__(self, path: Path, query: str, target: str) -> None:
        super().__init__()
        self._path = path
        self._query = query
        self._target = target
        self._cancelled = False

    def cancel(self) -> None:
        self._cancelled = True

    @Slot()
    def run(self) -> None:
        try:
            alignments = parse_paf_pair(
                self._path,
                self._query,
                self._target,
                progress=self.progress.emit,
                is_cancelled=lambda: self._cancelled,
            )
        except PafError as exc:
            self.failed.emit(str(exc))
            return
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(f"Could not load pair: {exc}")
            return
        self.finished.emit(alignments, self._query, self._target)

_FORWARD_COLOR = QColor(40, 90, 180)   # direct repeats (blue)
_REVERSE_COLOR = QColor(200, 60, 40)   # inverted repeats (red)
# Charcoal — must stay distinct from warm reverse red and cool forward blue.
_DIAGONAL_COLOR = QColor(45, 45, 45)

_DEFAULT_LINE_WIDTH_STEP = 3
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
# Shift+drag zoom: pixels of vertical motion per one +/- step (matches PNG preview).
_DRAG_ZOOM_PIXELS_PER_STEP = 28.0
# Trackpad pinch magnification accumulated per zoom step (matches PNG preview).
_PINCH_DELTA_PER_STEP = 0.04

_RULER_TOP = 18
# Left margin: room for tick labels flush against the axis (same px/base as X).
_RULER_LEFT = 44
# Multi-sequence All / multi_dot: room for sequence names + per-block ticks.
_RULER_TOP_MULTI = 30
_RULER_LEFT_MULTI = 52
_Y_TICK_PAD = 3  # gap between label and tick mark
_RULER_BG = QColor(250, 250, 250)
_RULER_FG = QColor(110, 110, 110)
_RULER_LINE = QColor(200, 200, 200)
# Multi-sequence grid: light grey alone disappears under dense small-word hits.
# White halo + dark core keep block boundaries readable at Fit zoom.
# Integer pixel widths only — fractional cosmetic pens (e.g. 1.75) round
# unevenly by subpixel position, so some dividers looked thicker than others.
_SEPARATOR_HALO = QColor(255, 255, 255, 230)
_SEPARATOR_LINE = QColor(25, 25, 25)
_SEPARATOR_HALO_PX = 4
_SEPARATOR_LINE_PX = 2

_SELECT_FILL = QColor(70, 130, 220, 55)
_SELECT_EDGE = QColor(60, 110, 200, 190)

_CLICK_SLOP = 4
# Selection is stored in scene (base) units. At high zoom a normal drag is often
# < 1 bp wide; expand to this minimum so the rubber band still commits.
_MIN_SELECTION_BP = 1.0


def _ensure_min_bp_rect(
    rect: QRectF, bounds: QRectF, *, min_size: float = _MIN_SELECTION_BP
) -> QRectF:
    """Grow ``rect`` to at least ``min_size``×``min_size`` bases, clamped to ``bounds``."""
    r = QRectF(rect)
    if r.width() < min_size:
        cx = r.center().x()
        r.setLeft(cx - min_size / 2.0)
        r.setWidth(min_size)
    if r.height() < min_size:
        cy = r.center().y()
        r.setTop(cy - min_size / 2.0)
        r.setHeight(min_size)
    r = r.intersected(bounds)
    if r.width() < min_size and bounds.width() >= min_size:
        if r.left() <= bounds.left() + 1e-9:
            r.setRect(bounds.left(), r.top(), min_size, r.height())
        else:
            r.setRect(bounds.right() - min_size, r.top(), min_size, r.height())
        r = r.intersected(bounds)
    if r.height() < min_size and bounds.height() >= min_size:
        if r.top() <= bounds.top() + 1e-9:
            r.setRect(r.left(), bounds.top(), r.width(), min_size)
        else:
            r.setRect(r.left(), bounds.bottom() - min_size, r.width(), min_size)
        r = r.intersected(bounds)
    return r


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


def _axis_entry_at(entries: list, scene_coord: float):
    """Return the ``AxisEntry`` whose block contains ``scene_coord``, if any."""
    if not entries:
        return None
    x = float(scene_coord)
    for index, entry in enumerate(entries):
        lo = float(entry.offset)
        hi = float(entry.offset + entry.length)
        # Half-open [lo, hi) except the last block includes its right end.
        if lo <= x < hi or (index == len(entries) - 1 and lo <= x <= hi):
            return entry
    return None


def _scene_to_local_bp(scene_coord: float | int, entries: list) -> int | None:
    """1-based position within the sequence block under ``scene_coord``."""
    entry = _axis_entry_at(entries, float(scene_coord))
    if entry is None:
        return None
    local = float(scene_coord) - float(entry.offset)
    return _scene_to_bp(local, length=int(entry.length), offset=0)


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
        identity_spans: list[tuple[float, float, float, float]] | None = None,
        mirror_symmetric: bool = True,
    ) -> None:
        super().__init__()
        self._result = result
        self._w = max(int(result.sequence_length), 1)
        self._h = max(int(result.effective_target_length), 1)
        self._line_width_step = max(
            _MIN_LINE_WIDTH_STEP, min(_MAX_LINE_WIDTH_STEP, int(line_width_step))
        )
        # Self-match diagonals: full (0,0)→(n,n) for one sequence, or one
        # span per shared sequence block in a multi_dot / All-sequences layout.
        # Empty list = cross comparison (no identity line).
        if identity_spans is None:
            n = float(min(self._w, self._h))
            self._identity_spans: list[tuple[float, float, float, float]] = [
                (0.0, 0.0, n, n)
            ]
        else:
            self._identity_spans = list(identity_spans)
        # Mirror forward hits across y=x — only valid for true self-dotplots.
        self._mirror_symmetric = bool(mirror_symmetric)
        # Optional cosmetic width override for high-res PNG/SVG export.
        self._export_pen_px: float | None = None
        self.setFlag(QGraphicsItem.GraphicsItemFlag.ItemUsesExtendedStyleOption, True)
        # Cache the dense hit layer in device pixels so mouse-move overlays
        # (crosshair) do not re-rasterize hundreds of thousands of segments.
        self.setCacheMode(QGraphicsItem.CacheMode.DeviceCoordinateCache)

    def set_line_width_step(self, step: int) -> None:
        step = max(_MIN_LINE_WIDTH_STEP, min(_MAX_LINE_WIDTH_STEP, int(step)))
        if step == self._line_width_step:
            return
        self._line_width_step = step
        self.update()

    def set_export_pen_px(self, width: float | None) -> None:
        """Temporary cosmetic pen width for export (``None`` restores UI width)."""
        self._export_pen_px = None if width is None else max(1.0, float(width))
        self.update()

    def boundingRect(self) -> QRectF:  # noqa: N802
        # Pad half a scene unit so the cosmetic right/bottom axis strokes at
        # the far edges are not clipped by the item bounds.
        return QRectF(0.0, 0.0, float(self._w), float(self._h)).adjusted(
            0.0, 0.0, 0.5, 0.5
        )

    def _pen(self, color: QColor, *, square_cap: bool = False) -> QPen:
        pen = QPen(color)
        pen.setCosmetic(True)
        if self._export_pen_px is not None:
            pen.setWidthF(self._export_pen_px)
        else:
            pen.setWidthF(_line_width_px(self._line_width_step))
        # Square caps give the chunky staircase identity look when zoomed.
        pen.setCapStyle(
            Qt.PenCapStyle.SquareCap if square_cap else Qt.PenCapStyle.RoundCap
        )
        return pen

    def paint(self, painter: QPainter, option, widget=None) -> None:  # noqa: N802
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        # White plot face; viewport outside the plot uses the ruler grey.
        # Top/left edges are the X/Y ruler axes; close the rectangle on the
        # right/bottom with the same light axis stroke (no scale labels).
        exposed = option.exposedRect if option is not None else self.boundingRect()
        painter.fillRect(exposed, QColor(255, 255, 255))
        w, h = float(self._w), float(self._h)
        axis = QPen(_RULER_LINE)
        axis.setCosmetic(True)
        axis.setWidth(1)
        painter.setPen(axis)
        painter.drawLine(QPointF(w, 0.0), QPointF(w, h))  # right
        painter.drawLine(QPointF(0.0, h), QPointF(w, h))  # bottom
        scale = painter.transform().m11() or 1.0
        # Scene length that maps to ~1 screen pixel. Shorter hits still get a
        # point so Fit zoom is not an empty diagonal on long sequences.
        min_length = 0.0 if scale >= 1.0 else (1.0 / scale)
        visible = exposed

        painter.setPen(self._pen(_DIAGONAL_COLOR, square_cap=True))
        for x1, y1, x2, y2 in self._identity_spans:
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))

        forward: list[tuple[float, float, float, float]] = []
        reverse: list[tuple[float, float, float, float]] = []
        mirror = self._mirror_symmetric
        # On multi-Mb axes at Fit, short hits are << 1 px — stretch each
        # stroke to a few screen pixels so sparse FastGA pairs stay visible.
        min_draw = max(min_length, 4.0 / scale)
        for seg in self._result.segments:
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
                    if seg.strand == STRAND_REVERSE or not mirror:
                        continue
                    mirrored_off = (
                        hi_y < visible.left() or lo_y > visible.right()
                        or hi_x < visible.top() or lo_x > visible.bottom()
                    )
                    if mirrored_off:
                        continue
            if seg.length < min_draw:
                if seg.strand == STRAND_REVERSE:
                    x2, y2 = x1 + min_draw, y1 - min_draw
                else:
                    x2, y2 = x1 + min_draw, y1 + min_draw
            if seg.strand == STRAND_REVERSE:
                reverse.append((x1, y1, x2, y2))
            else:
                forward.append((x1, y1, x2, y2))

        painter.setPen(self._pen(_FORWARD_COLOR))
        for x1, y1, x2, y2 in forward:
            painter.drawLine(QPointF(x1, y1), QPointF(x2, y2))
            if mirror:
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
        self.setMouseTracking(True)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        """Show the full sequence name when the elided axis label is hovered."""
        view = self._view
        entries = view.axis_entries(self._horizontal)
        tip = ""
        if entries and view.plot_extent(self._horizontal) > 0:
            pos = float(event.position().x() if self._horizontal else event.position().y())
            for entry in entries:
                start = float(view.scene_to_ruler(float(entry.offset), self._horizontal))
                end = float(
                    view.scene_to_ruler(
                        float(entry.offset + entry.length), self._horizontal
                    )
                )
                lo, hi = (start, end) if start <= end else (end, start)
                if lo <= pos <= hi:
                    tip = str(entry.name)
                    break
        if tip != self.toolTip():
            self.setToolTip(tip)
        super().mouseMoveEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802
        self.setToolTip("")
        super().leaveEvent(event)

    def _paint_sequence_axis_labels(
        self, painter: QPainter, view: "DotplotCanvas", entries: list
    ) -> None:
        """Draw multi-sequence names, elided to each block so they do not overlap.

        Full names are often longer than a panel at Fit zoom; clipping each
        label to its sequence's pixel span keeps neighbours readable. Zoom in
        on a block to see more of the name.
        """
        painter.setPen(QPen(_RULER_FG))
        metrics = painter.fontMetrics()
        # Below this many pixels a label is just noise — skip it.
        min_px = max(14, metrics.horizontalAdvance("…") + 4)
        pad = 4.0
        for entry in entries:
            start = float(view.scene_to_ruler(float(entry.offset), self._horizontal))
            end = float(
                view.scene_to_ruler(
                    float(entry.offset + entry.length), self._horizontal
                )
            )
            lo, hi = (start, end) if start <= end else (end, start)
            avail = hi - lo - pad
            if avail < min_px:
                continue
            centre = (lo + hi) / 2.0
            if self._horizontal:
                if centre < 0 or centre > self.width():
                    continue
            else:
                if centre < 0 or centre > self.height():
                    continue
            text = metrics.elidedText(
                str(entry.name),
                Qt.TextElideMode.ElideMiddle,
                max(1, int(avail)),
            )
            if not text:
                continue
            w = metrics.horizontalAdvance(text)
            if self._horizontal:
                # Outer edge of the strip (away from the plot / tick marks).
                painter.drawText(
                    int(centre - w / 2),
                    metrics.ascent() + 1,
                    text,
                )
            else:
                painter.save()
                # Outer (left) edge so tick numbers can sit near the plot.
                painter.translate(metrics.height() + 2, centre)
                painter.rotate(-90)
                painter.drawText(int(-w / 2), 0, text)
                painter.restore()

    def _paint_tick_marks(
        self,
        painter: QPainter,
        view: "DotplotCanvas",
        *,
        ticks: list[int],
        positions: dict[int, float],
        labeled: set[int],
        length: int,
        offset: int,
    ) -> None:
        """Draw tick marks and labels for one axis span (scene coords)."""
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
                    label_h = painter.fontMetrics().height()
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

    def _paint_cursor_label(
        self,
        painter: QPainter,
        view: "DotplotCanvas",
        *,
        font: QFont,
        length: int,
        offset: int,
        entries: list | None = None,
    ) -> None:
        cursor = view.cursor_scene_pos()
        if cursor is None:
            return
        cursor_font = QFont(font)
        cursor_font.setBold(True)
        painter.setFont(cursor_font)
        value = cursor.x() if self._horizontal else cursor.y()
        pos = view.scene_to_ruler(value, self._horizontal)
        if entries:
            local = _scene_to_local_bp(value, entries)
            if local is None:
                return
            real = local
        else:
            real = _scene_to_bp(value, length=length, offset=offset)
        label = f"{real:,}"
        cursor_metrics = painter.fontMetrics()
        pad = 2
        if self._horizontal and 0 <= pos <= self.width():
            x0 = int(pos) + 2
            y_base = self.height() - 6
            tw = cursor_metrics.horizontalAdvance(label)
            th = cursor_metrics.height()
            bg = QRect(
                x0 - pad,
                y_base - cursor_metrics.ascent() - pad,
                tw + 2 * pad,
                th + 2 * pad,
            )
            painter.fillRect(bg, _RULER_BG)
            painter.setPen(QPen(QColor(40, 90, 180)))
            painter.drawText(x0, y_base, label)
        elif not self._horizontal and 0 <= pos <= self.height():
            label_h = cursor_metrics.height()
            text_rect = QRect(
                0,
                int(pos) - label_h // 2,
                self.width() - 4 - _Y_TICK_PAD,
                label_h,
            )
            tw = cursor_metrics.horizontalAdvance(label)
            bg = QRect(
                text_rect.right() - tw - pad,
                text_rect.top() - pad,
                tw + 2 * pad,
                label_h + 2 * pad,
            )
            painter.fillRect(bg, _RULER_BG)
            painter.setPen(QPen(QColor(40, 90, 180)))
            painter.drawText(
                text_rect,
                int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter),
                label,
            )

    def paintEvent(self, event) -> None:  # noqa: N802
        view = self._view
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        painter.fillRect(self.rect(), _RULER_BG)
        if view.plot_extent(self._horizontal) <= 0:
            return

        font = QFont(self.font())
        font.setPointSize(8)
        painter.setFont(font)

        visible = view.visible_scene_rect()
        length = view.plot_extent(self._horizontal)
        step = view.shared_tick_step()
        lo = visible.left() if self._horizontal else visible.top()
        hi = visible.right() if self._horizontal else visible.bottom()
        metrics = painter.fontMetrics()
        offset = view.coordinate_offset()
        entries = view.axis_entries(self._horizontal)

        def _label_extent(pos: float, label: str) -> tuple[float, float]:
            if self._horizontal:
                width = float(metrics.horizontalAdvance(label))
                return pos + 2.0, pos + 2.0 + width
            height = float(metrics.height())
            return pos - height / 2.0, pos + height / 2.0

        painter.setPen(QPen(_RULER_FG))
        if entries:
            self._paint_sequence_axis_labels(painter, view, entries)
            # Per-sequence local ticks (1-based within each block).
            for index, entry in enumerate(entries):
                block_lo = float(entry.offset)
                block_hi = float(entry.offset + entry.length)
                if block_hi < lo or block_lo > hi:
                    continue
                local_lo = max(0.0, lo - block_lo)
                local_hi = min(float(entry.length), hi - block_lo)
                local_ticks = _axis_tick_positions(
                    lo=local_lo,
                    hi=local_hi,
                    length=int(entry.length),
                    step=step,
                )
                # Shared block boundary: next sequence's start (label 1) owns
                # the tick — drop this block's end so labels do not stack.
                if index + 1 < len(entries):
                    local_ticks = [
                        t for t in local_ticks if t < int(entry.length)
                    ]
                # Positions dict keyed by local tick for overlap logic.
                local_positions = {
                    tick: view.scene_to_ruler(
                        block_lo + float(tick), self._horizontal
                    )
                    for tick in local_ticks
                }
                labeled = _non_overlapping_tick_labels(
                    local_ticks,
                    length=int(entry.length),
                    offset=0,
                    positions=local_positions,
                    label_extent=_label_extent,
                )
                self._paint_tick_marks(
                    painter,
                    view,
                    ticks=local_ticks,
                    positions=local_positions,
                    labeled=labeled,
                    length=int(entry.length),
                    offset=0,
                )
            self._paint_cursor_label(
                painter,
                view,
                font=font,
                length=length,
                offset=offset,
                entries=entries,
            )
        else:
            ticks = _axis_tick_positions(lo=lo, hi=hi, length=length, step=step)
            positions = {
                tick: view.scene_to_ruler(tick, self._horizontal) for tick in ticks
            }
            labeled = _non_overlapping_tick_labels(
                ticks,
                length=length,
                offset=offset,
                positions=positions,
                label_extent=_label_extent,
            )
            self._paint_tick_marks(
                painter,
                view,
                ticks=ticks,
                positions=positions,
                labeled=labeled,
                length=length,
                offset=offset,
            )
            self._paint_cursor_label(
                painter,
                view,
                font=font,
                length=length,
                offset=offset,
            )

        painter.setPen(QPen(_RULER_LINE))
        if self._horizontal:
            painter.drawLine(0, self.height() - 1, self.width(), self.height() - 1)
        else:
            painter.drawLine(self.width() - 1, 0, self.width() - 1, self.height())


class _PlotOverlay(QWidget):
    """Draws crosshair / selection in viewport space without re-painting hits."""

    def __init__(self, view: "DotplotCanvas") -> None:
        super().__init__(view.viewport())
        self._view = view
        self.setAttribute(Qt.WidgetAttribute.WA_TransparentForMouseEvents)
        self.setAttribute(Qt.WidgetAttribute.WA_NoSystemBackground)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)

    def paintEvent(self, event) -> None:  # noqa: N802
        view = self._view
        if view._item is None or not view._length:
            return
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)

        # Live drag is tracked in viewport pixels so the rubber band stays under
        # the cursor after zoom (scene↔view mapping can lag or disagree).
        rect: QRect | None = None
        if view._drag_start_vp is not None and view._drag_now_vp is not None:
            rect = QRect(view._drag_start_vp, view._drag_now_vp).normalized()
        elif view._selection is not None:
            tl = view.mapFromScene(view._selection.topLeft())
            br = view.mapFromScene(view._selection.bottomRight())
            rect = QRect(tl, br).normalized()
        if rect is not None and rect.width() >= 1 and rect.height() >= 1:
            painter.fillRect(rect, _SELECT_FILL)
            edge = QPen(_SELECT_EDGE)
            edge.setCosmetic(True)
            painter.setPen(edge)
            painter.setBrush(Qt.BrushStyle.NoBrush)
            painter.drawRect(rect)

        # Crosshair from viewport mouse position — same space as the rubber band.
        if view._cursor_vp is not None:
            cross = QPen(QColor(90, 140, 220, 150))
            cross.setCosmetic(True)
            painter.setPen(cross)
            pt = view._cursor_vp
            painter.drawLine(pt.x(), 0, pt.x(), self.height())
            painter.drawLine(0, pt.y(), self.width(), pt.y())
        painter.end()


class DotplotCanvas(QGraphicsView):
    """Zoom/pan canvas with rulers, rectangular selection and a crosshair."""

    positionChanged = Signal(int, int)
    rectSelected = Signal(int, int, int, int)  # x_start, x_end, y_start, y_end
    zoomChanged = Signal()  # scale changed (slider / fit / wheel / pinch)

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
        self._length = 0  # X / query extent (bp)
        self._height = 0  # Y / target extent (bp); 0 → same as ``_length``
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
        self._zoom_dragging = False
        self._zoom_drag_origin: QPointF | None = None
        self._zoom_drag_accum = 0.0
        self._wheel_zoom_accum = 0
        self._pinch_accum = 0.0
        # Rubber-band corners in viewport pixels (follow the mouse exactly).
        self._drag_start_vp: QPoint | None = None
        self._drag_now_vp: QPoint | None = None
        # Cursor in viewport pixels (crosshair) and scene coords (rulers / status).
        self._cursor_vp: QPoint | None = None
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
        # Crosshair / selection overlay — mouse move must not re-paint all hits.
        self._overlay = _PlotOverlay(self)
        # Repaint the rulers whenever the view scrolls, so ticks track the plot.
        self.horizontalScrollBar().valueChanged.connect(self._refresh_rulers)
        self.verticalScrollBar().valueChanged.connect(self._refresh_rulers)
        self.horizontalScrollBar().valueChanged.connect(self._overlay.update)
        self.verticalScrollBar().valueChanged.connect(self._overlay.update)
        self.viewport().setAttribute(Qt.WidgetAttribute.WA_AcceptTouchEvents, True)

    def set_title(self, title: str) -> None:
        """Sequence name for status / PDF (not drawn on the y-axis)."""
        self._title = title or ""
        self._refresh_rulers()

    # -- accessors used by the ruler widgets -------------------------------

    def plot_length(self) -> int:
        """X-axis (query) length in bp."""
        return self._length

    def plot_height(self) -> int:
        """Y-axis (target) length in bp; equals ``plot_length`` for self plots."""
        return self._height if self._height > 0 else self._length

    def plot_extent(self, horizontal: bool) -> int:
        """Axis length for the given ruler orientation."""
        return self._length if horizontal else self.plot_height()

    def coordinate_offset(self) -> int:
        return self._offset

    def _display_xy(self, x: int, y: int) -> tuple[int, int]:
        """1-based bp for status/cursor; local within each multi-seq block."""
        height = self.plot_height()
        if self._axis_x:
            lx = _scene_to_local_bp(x, self._axis_x)
            bx = lx if lx is not None else _scene_to_bp(
                x, length=self._length, offset=self._offset
            )
        else:
            bx = _scene_to_bp(x, length=self._length, offset=self._offset)
        if self._axis_y:
            ly = _scene_to_local_bp(y, self._axis_y)
            by = ly if ly is not None else _scene_to_bp(
                y, length=height, offset=self._offset
            )
        else:
            by = _scene_to_bp(y, length=height, offset=self._offset)
        return bx, by

    def set_axes(self, x_entries: list, y_entries: list) -> None:
        """Sequence layout for the concatenated view (empty for a single pair)."""
        self._axis_x = list(x_entries or [])
        self._axis_y = list(y_entries or [])
        self._apply_ruler_margins()
        self.viewport().update()
        self._refresh_rulers()

    def axis_entries(self, horizontal: bool) -> list:
        return self._axis_x if horizontal else self._axis_y

    def has_multi_axes(self) -> bool:
        return bool(self._axis_x or self._axis_y)

    def ruler_top(self) -> int:
        return _RULER_TOP_MULTI if self.has_multi_axes() else _RULER_TOP

    def ruler_left(self) -> int:
        return _RULER_LEFT_MULTI if self.has_multi_axes() else _RULER_LEFT

    def _apply_ruler_margins(self) -> None:
        self.setViewportMargins(self.ruler_left(), self.ruler_top(), 0, 0)
        self._position_rulers()

    def _identity_spans_for(
        self, result: DotplotResult
    ) -> tuple[list[tuple[float, float, float, float]], bool]:
        """Return (identity line spans, mirror_symmetric).

        Multi-sequence ``*_multi_dot.paf`` / All-sequences: one thick diagonal
        per shared sequence block (figure-2 style when zoomed), not one line
        across the whole grid. Cross-only pairs get no identity line.
        """
        if self._axis_x or self._axis_y:
            y_by_name = {e.name: e for e in self._axis_y}
            spans: list[tuple[float, float, float, float]] = []
            for xe in self._axis_x:
                ye = y_by_name.get(xe.name)
                if ye is None:
                    continue
                length = min(int(xe.length), int(ye.length))
                if length <= 0:
                    continue
                spans.append(
                    (
                        float(xe.offset),
                        float(ye.offset),
                        float(xe.offset + length),
                        float(ye.offset + length),
                    )
                )
            # One shared sequence on both axes (self pair with name labels):
            # keep triangle mirroring. Multi-sequence matrices stay unmirrored
            # (pairwise files already store / reflect the other triangle).
            mirror = (
                len(spans) == 1
                and len(self._axis_x) == 1
                and len(self._axis_y) == 1
                and self._axis_x[0].name == self._axis_y[0].name
            )
            return spans, mirror
        label = str(getattr(result, "sequence_id", "") or "")
        tid = str(getattr(result, "target_id", "") or "").strip()
        if " vs " in label or (tid and tid != label):
            return [], False
        n = float(
            max(
                min(int(result.sequence_length), int(result.effective_target_length)),
                1,
            )
        )
        return [(0.0, 0.0, n, n)], True

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
        width = max(float(self._length), 1.0)
        height = max(float(self.plot_height()), 1.0)
        # Intersect with the plot rectangle so empty letterbox beside a fitted
        # plot does not inflate the span (that produced only 1 and 5,000 on a 6 kb
        # sequence in a wide window).
        shown = visible.intersected(QRectF(0.0, 0.0, width, height))
        span = max(float(shown.width()), float(shown.height()), 1.0)
        return _nice_tick_step(span, target_ticks=10)

    def _refresh_rulers(self) -> None:
        self._top_ruler.update()
        self._left_ruler.update()

    def _refresh_overlay(self) -> None:
        self._overlay.update()

    def _position_rulers(self) -> None:
        vp = self.viewport().geometry()
        top = self.ruler_top()
        left = self.ruler_left()
        # Only the top-left margin corner — plot stays white.
        self._ruler_corner.setGeometry(0, 0, max(vp.x(), 1), max(vp.y(), 1))
        self._ruler_corner.raise_()
        self._top_ruler.setGeometry(vp.x(), 0, vp.width(), top)
        self._left_ruler.setGeometry(0, vp.y(), left, vp.height())
        self._top_ruler.raise_()
        self._left_ruler.raise_()
        # Overlay covers the plot viewport only (child of viewport()).
        self._overlay.setGeometry(self.viewport().rect())
        self._overlay.raise_()
        self._overlay.show()

    # ------------------------------------------------------------------
    # content
    # ------------------------------------------------------------------

    def show_result(self, result: DotplotResult, *, preserve_view: bool = False) -> None:
        plot_h = int(result.effective_target_length)
        same_length = (
            preserve_view
            and self._length == result.sequence_length
            and self.plot_height() == plot_h
            and self._length > 0
            and self._item is not None
        )
        # Restore exact scroll offsets — not centerOn(). centerOn maps through
        # floating-point scene coords and integer scrollbars, drifting ~1 px per
        # word-size recompute; after a zoom that crawl makes hits miss the
        # crosshair (same bug on single-sequence and multi-sequence Single).
        saved_transform = self.transform() if same_length else None
        saved_fit = self._fit_scale if same_length else 0.0
        saved_scroll = (
            (
                self.horizontalScrollBar().value(),
                self.verticalScrollBar().value(),
            )
            if same_length
            else None
        )

        scene = self.scene()
        scene.clear()
        identity_spans, mirror = self._identity_spans_for(result)
        self._item = _SegmentsItem(
            result,
            line_width_step=self._line_width_step,
            identity_spans=identity_spans,
            mirror_symmetric=mirror,
        )
        self._length = int(result.sequence_length)
        self._height = plot_h
        self._offset = int(getattr(result, "coordinate_offset", 0) or 0)
        self._selection = None
        scene.addItem(self._item)
        scene.setSceneRect(0, 0, self._length, self._height)

        if saved_transform is not None and saved_scroll is not None:
            # Keep pan/zoom when only word size / reverse-complement toggled.
            self.setTransform(saved_transform)
            self._fit_scale = saved_fit
            self.horizontalScrollBar().setValue(saved_scroll[0])
            self.verticalScrollBar().setValue(saved_scroll[1])
            self._needs_fit = False
            self._position_rulers()
            self.viewport().update()
            self._refresh_overlay()
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
        self._height = 0
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
        """Fit the whole plot into the viewport.

        X and Y use the same pixels-per-base so homology diagonals stay at
        45°. Self-dotplots are square; pairwise plots with unequal sequence
        lengths are rectangular. ``fitInView`` is avoided: with viewport
        margins it can leave a non-uniform transform and stretch the plot.
        """
        if not self._length:
            return
        vp = self.viewport().rect()
        if vp.width() < 1 or vp.height() < 1:
            return
        width = float(self._length)
        height = float(self.plot_height())
        scale = min(vp.width() / width, vp.height() / height)
        self._set_uniform_scale(scale)
        self._fit_scale = scale
        self._clamp_scroll_to_plot()
        self._position_rulers()
        self.viewport().update()
        self._refresh_overlay()
        self._refresh_rulers()
        self.zoomChanged.emit()

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
            QRectF(0, 0, float(self._length), float(self.plot_height()))
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
        self._position_rulers()
        self.viewport().update()
        self._refresh_overlay()
        self._refresh_rulers()
        self.zoomChanged.emit()
        return True

    def zoom_to_hits(self, result: DotplotResult, *, pad_frac: float = 0.08) -> bool:
        """Zoom to the bounding box of plotted hits (sparse pairs on huge axes)."""
        segs = getattr(result, "segments", None) or []
        if not segs or not self._length:
            return False
        x0 = y0 = 1e30
        x1 = y1 = -1e30
        for seg in segs:
            xa, xb = float(seg.x), float(seg.end_x)
            ya, yb = float(seg.y), float(seg.end_y)
            x0 = min(x0, xa, xb)
            x1 = max(x1, xa, xb)
            y0 = min(y0, ya, yb)
            y1 = max(y1, ya, yb)
        if x1 <= x0 or y1 <= y0:
            return False
        width = float(self._length)
        height = float(self.plot_height())
        # Pad in scene units; keep at least ~1% of each axis so tiny boxes expand.
        pad_x = max((x1 - x0) * pad_frac, width * 0.01, 50.0)
        pad_y = max((y1 - y0) * pad_frac, height * 0.01, 50.0)
        target = QRectF(x0 - pad_x, y0 - pad_y, (x1 - x0) + 2 * pad_x, (y1 - y0) + 2 * pad_y)
        target = target.intersected(QRectF(0.0, 0.0, width, height))
        if target.width() < 1.0 or target.height() < 1.0:
            return False
        # Whole-axis hits: Fit is enough.
        if target.width() >= width * 0.85 and target.height() >= height * 0.85:
            return False
        vp = self.viewport().rect()
        if vp.width() < 1 or vp.height() < 1:
            return False
        side = max(target.width(), target.height())
        scale = min(vp.width() / side, vp.height() / side)
        scale = min(scale, _MAX_PIXELS_PER_BASE)
        self._needs_fit = False
        self.setTransform(QTransform.fromScale(scale, scale))
        self.centerOn(target.center())
        self._clamp_scroll_to_plot()
        self._position_rulers()
        self.viewport().update()
        self._refresh_overlay()
        self._refresh_rulers()
        self.zoomChanged.emit()
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
        self._refresh_overlay()

    def clear_selection(self) -> None:
        self._selection = None
        self._drag_start_vp = self._drag_now_vp = None
        self._refresh_overlay()

    def _viewport_event_pos(self, event) -> QPoint:
        """Mouse/wheel position in viewport coordinates for ``mapToScene``.

        Events arrive through the viewport (typical) with viewport-local coords.
        If the originating widget is the view itself, map into the viewport so
        ruler margins are not treated as plot pixels.
        """
        pos = event.position().toPoint()
        widget = event.widget() if hasattr(event, "widget") else None
        if widget is self:
            return self.viewport().mapFrom(self, pos)
        return pos

    def _clear_drag(self) -> None:
        self._drag_start_vp = self._drag_now_vp = None

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

    def _sequence_separator_segments(
        self, width: float, height: float
    ) -> list[tuple[QPointF, QPointF]]:
        """Internal grid lines between concatenated sequences (scene coords)."""
        if not (self._axis_x or self._axis_y):
            return []
        lines: list[tuple[QPointF, QPointF]] = []
        for entry in self._axis_x[1:]:
            x = float(entry.offset)
            lines.append((QPointF(x, 0.0), QPointF(x, height)))
        for entry in self._axis_y[1:]:
            y = float(entry.offset)
            lines.append((QPointF(0.0, y), QPointF(width, y)))
        return lines

    @staticmethod
    def _paint_separator_lines(
        painter: QPainter, lines: list[tuple[QPointF, QPointF]]
    ) -> None:
        """Draw sequence dividers with a white halo so they cut through dense hits.

        ``lines`` must already be in the painter's current coordinate system
        (prefer integer device/viewport pixels so every divider is the same
        thickness).
        """
        if not lines:
            return
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        halo = QPen(_SEPARATOR_HALO)
        halo.setCosmetic(True)
        halo.setWidth(_SEPARATOR_HALO_PX)
        halo.setCapStyle(Qt.PenCapStyle.FlatCap)
        core = QPen(_SEPARATOR_LINE)
        core.setCosmetic(True)
        core.setWidth(_SEPARATOR_LINE_PX)
        core.setCapStyle(Qt.PenCapStyle.FlatCap)
        painter.setPen(halo)
        for a, b in lines:
            painter.drawLine(a, b)
        painter.setPen(core)
        for a, b in lines:
            painter.drawLine(a, b)

    def _separator_lines_viewport_px(
        self,
    ) -> list[tuple[QPointF, QPointF]]:
        """Map scene separators to integer viewport pixels (uniform stroke)."""
        scene_lines = self._sequence_separator_segments(
            float(self._length), float(self.plot_height())
        )
        out: list[tuple[QPointF, QPointF]] = []
        for a, b in scene_lines:
            pa = self.mapFromScene(a)
            pb = self.mapFromScene(b)
            # Snap to the pixel grid so Fit zoom does not make some lines 1 px
            # and others 2 px depending on fractional scene→device mapping.
            out.append(
                (
                    QPointF(round(pa.x()) + 0.5, round(pa.y()) + 0.5),
                    QPointF(round(pb.x()) + 0.5, round(pb.y()) + 0.5),
                )
            )
        return out

    def drawForeground(self, painter: QPainter, rect: QRectF) -> None:  # noqa: N802
        """Overlays for multi-sequence separators (after hits).

        Crosshair and selection live on ``_PlotOverlay`` so mouse moves do not
        re-rasterize the hit layer. Separators use a white halo and are painted
        in viewport pixels so every divider has the same thickness.
        """
        super().drawForeground(painter, rect)
        if self._item is None or not self._length:
            return
        lines = self._separator_lines_viewport_px()
        if not lines:
            return
        # drawForeground's painter is in scene space; paint in viewport pixels.
        painter.save()
        painter.setTransform(QTransform())
        self._paint_separator_lines(painter, lines)
        painter.restore()

    def _zoom_allowed(self, target: float) -> bool:
        floor = (self._fit_scale or target) * _ZOOM_OUT_LIMIT
        return floor <= target <= _MAX_PIXELS_PER_BASE

    def _clamp_scroll_to_plot(self) -> None:
        """Keep empty scene space out of the viewport (pin plot to top-left)."""
        if not self._length:
            return
        width = float(self._length)
        height = float(self.plot_height())
        vis = self.mapToScene(self.viewport().rect()).boundingRect()
        dx = dy = 0.0
        if vis.width() <= width + 1e-6:
            if vis.left() < 0:
                dx = -vis.left()
            elif vis.right() > width:
                dx = width - vis.right()
        else:
            # Viewport wider than the plot: keep x=0 on the left edge.
            dx = -vis.left()
        if vis.height() <= height + 1e-6:
            if vis.top() < 0:
                dy = -vis.top()
            elif vis.bottom() > height:
                dy = height - vis.bottom()
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

    def _pan_scene_by(self, dx: float, dy: float) -> None:
        """Shift the view by ``dx``/``dy`` scene (bp) units."""
        if abs(dx) < 1e-9 and abs(dy) < 1e-9:
            return
        origin = self.mapFromScene(QPointF(0, 0))
        shifted = self.mapFromScene(QPointF(dx, dy))
        self.horizontalScrollBar().setValue(
            self.horizontalScrollBar().value() + int(round(shifted.x() - origin.x()))
        )
        self.verticalScrollBar().setValue(
            self.verticalScrollBar().value() + int(round(shifted.y() - origin.y()))
        )

    def _ensure_selection_visible(self) -> None:
        """Pan so the current selection stays inside the viewport when possible."""
        sel = self._selection
        if sel is None or not self._length:
            return
        if sel.width() < 1.0 or sel.height() < 1.0:
            return
        vis = self.mapToScene(self.viewport().rect()).boundingRect()
        if vis.width() < 1e-6 or vis.height() < 1e-6:
            return
        # Selection larger than the view: keep its centre on screen.
        if sel.width() >= vis.width() or sel.height() >= vis.height():
            self.centerOn(sel.center())
            return
        dx = dy = 0.0
        if sel.left() < vis.left():
            dx = sel.left() - vis.left()
        elif sel.right() > vis.right():
            dx = sel.right() - vis.right()
        if sel.top() < vis.top():
            dy = sel.top() - vis.top()
        elif sel.bottom() > vis.bottom():
            dy = sel.bottom() - vis.bottom()
        self._pan_scene_by(dx, dy)

    def _selection_anchor_vp(self) -> QPoint | None:
        """Viewport point under the selection centre, if a selection exists."""
        if not self.has_selection() or self._selection is None:
            return None
        pt = self.mapFromScene(self._selection.center())
        return QPoint(pt)

    def _apply_zoom(self, factor: float, anchor: QPoint | None = None) -> None:
        """Zoom by ``factor`` about ``anchor`` (viewport coords), then clamp.

        When a rectangular selection exists, zoom about its centre and pan so
        the selection stays visible in the window.
        """
        if self._item is None or not self._length:
            return
        current = self._uniform_scale()
        floor = (self._fit_scale or current) * _ZOOM_OUT_LIMIT
        target = min(_MAX_PIXELS_PER_BASE, max(floor, current * factor))
        actual = target / current
        if abs(actual - 1.0) < 1e-9:
            return

        # Prefer the selection so zooming in/out never loses the highlight.
        sel_anchor = self._selection_anchor_vp()
        if sel_anchor is not None:
            anchor = sel_anchor
        elif anchor is None:
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
        self._ensure_selection_visible()
        self._clamp_scroll_to_plot()
        # Clamping to the plot can push the selection off-screen again on edges.
        self._ensure_selection_visible()
        self._clamp_scroll_to_plot()
        self._position_rulers()
        self.viewport().update()
        self._refresh_overlay()
        self._refresh_rulers()
        self.zoomChanged.emit()

    def zoom_by(self, factor: float) -> None:
        """Zoom about the selection (if any) or the centre of the view."""
        self._apply_zoom(factor, self.viewport().rect().center())

    def zoom_at(self, factor: float, view_pos) -> None:
        """Zoom keeping the scene point under ``view_pos`` fixed."""
        self._apply_zoom(factor, view_pos)

    def zoom_slider_value(self) -> int:
        """Map current scale to 0–100 (0 = fit, 100 = max zoom)."""
        import math

        fit = float(self._fit_scale or 0.0)
        current = float(self._uniform_scale())
        if fit <= 0 or current <= 0:
            return 0
        if current <= fit * 1.001:
            return 0
        max_s = float(_MAX_PIXELS_PER_BASE)
        if current >= max_s * 0.999:
            return 100
        lo = math.log(fit)
        hi = math.log(max_s)
        if hi <= lo:
            return 0
        t = (math.log(current) - lo) / (hi - lo)
        return int(round(max(0.0, min(100.0, t * 100.0))))

    def set_zoom_slider_value(self, value: int) -> None:
        """Set scale from a 0–100 slider (0 = fit, 100 = max zoom)."""
        import math

        value = max(0, min(100, int(value)))
        if value <= 0 or self._fit_scale <= 0:
            self.reset_zoom()
            return
        if not self._length:
            return
        fit = float(self._fit_scale)
        max_s = float(_MAX_PIXELS_PER_BASE)
        t = value / 100.0
        target = math.exp(math.log(fit) + t * (math.log(max_s) - math.log(fit)))
        target = min(max_s, max(fit, target))
        current = self._uniform_scale()
        if current <= 0:
            return
        self._apply_zoom(target / current, self.viewport().rect().center())

    @staticmethod
    def _a4_page_rect(dpi: int) -> QRect:
        """A4 page in device pixels at ``dpi`` — same geometry as QPdfWriter."""
        dpi = max(50, min(600, int(dpi)))
        width = max(1, int(round(210.0 / 25.4 * dpi)))
        height = max(1, int(round(297.0 / 25.4 * dpi)))
        return QRect(0, 0, width, height)

    @staticmethod
    def _export_font(
        base: QFont, point_size: float, unit: float, *, bold: bool = False
    ) -> QFont:
        """Font sized in pixels from A4 ``unit`` so PNG/SVG match PDF points."""
        # pt → px at the page DPI implied by unit (px per mm): dpi = unit * 25.4
        px = max(1, int(round(float(point_size) * unit * 25.4 / 72.0)))
        font = QFont(base)
        font.setPixelSize(px)
        font.setBold(bool(bold))
        return font

    def _paint_export_figure(
        self,
        painter: QPainter,
        page: QRect,
        *,
        title: str = "",
        subtitle: str = "",
    ) -> None:
        """Paint axes, plot, and legend into ``page`` (shared by PDF/PNG/SVG)."""
        unit = page.width() / 210.0  # device units per mm (A4 width)
        margin = int(12 * unit)
        gutter = int(16 * unit)  # room for tick labels
        head = int(26 * unit)  # title block + top tick labels
        foot = int(10 * unit)  # legend

        width = max(self._length, 1)
        height = max(self.plot_height(), 1)
        offset = self._offset
        plot_bounds = QRectF(0.0, 0.0, float(width), float(height))
        source = self.visible_scene_rect().intersected(plot_bounds)
        if source.width() < 1.0 or source.height() < 1.0:
            source = QRectF(plot_bounds)

        avail_w = page.width() - margin * 2 - gutter
        avail_h = page.height() - margin * 2 - head - foot - gutter
        aspect = source.width() / max(source.height(), 1e-9)
        if avail_w / max(avail_h, 1) > aspect:
            side_h = float(avail_h)
            side_w = side_h * aspect
        else:
            side_w = float(avail_w)
            side_h = side_w / aspect
        plot = QRectF(margin + gutter, margin + head, side_w, side_h)

        # The app may be in a dark theme; a figure should not be. Paint the
        # page and the plot area white and draw everything in dark ink.
        painter.fillRect(page, QColor(255, 255, 255))

        base_font = QFont(painter.font())

        if title:
            painter.setFont(self._export_font(base_font, 11, unit, bold=True))
            painter.setPen(QPen(QColor(20, 20, 20)))
            painter.drawText(int(margin), int(margin + 6 * unit), title)
        if subtitle:
            painter.setFont(self._export_font(base_font, 8, unit))
            painter.setPen(QPen(QColor(90, 90, 90)))
            painter.drawText(int(margin), int(margin + 12 * unit), subtitle)
        # --- the plot itself (current viewport only) ---
        painter.fillRect(plot, QColor(255, 255, 255))
        painter.save()
        painter.setClipRect(plot)
        self.scene().render(painter, plot, source)
        # View drawForeground is not part of scene.render — redraw separators.
        sx = plot.width() / max(source.width(), 1e-9)
        sy = plot.height() / max(source.height(), 1e-9)

        def _map_pt(pt: QPointF) -> QPointF:
            x = plot.left() + (pt.x() - source.left()) * sx
            y = plot.top() + (pt.y() - source.top()) * sy
            return QPointF(round(x) + 0.5, round(y) + 0.5)

        export_seps = [
            (_map_pt(a), _map_pt(b))
            for a, b in self._sequence_separator_segments(
                float(width), float(height)
            )
        ]
        self._paint_separator_lines(painter, export_seps)
        painter.restore()

        # --- axes for the visible range ---
        span = max(float(source.width()), float(source.height()), 1.0)
        step = _nice_tick_step(span, target_ticks=10)
        tick_len = 3 * unit
        painter.setFont(self._export_font(base_font, 7, unit))
        painter.setPen(QPen(QColor(70, 70, 70)))
        metrics = painter.fontMetrics()

        x_ticks = _axis_tick_positions(
            lo=float(source.left()),
            hi=float(source.right()),
            length=width,
            step=step,
        )
        y_ticks = _axis_tick_positions(
            lo=float(source.top()),
            hi=float(source.bottom()),
            length=height,
            step=step,
        )
        x_positions = {
            tick: plot.left()
            + ((tick - source.left()) / max(source.width(), 1e-9)) * plot.width()
            for tick in x_ticks
        }
        y_positions = {
            tick: plot.top()
            + ((tick - source.top()) / max(source.height(), 1e-9)) * plot.height()
            for tick in y_ticks
        }

        def _export_x_extent(pos: float, label: str) -> tuple[float, float]:
            half = metrics.horizontalAdvance(label) / 2.0
            return pos - half, pos + half

        def _export_y_extent(pos: float, label: str) -> tuple[float, float]:
            half = metrics.height() / 2.0
            return pos - half, pos + half

        x_labeled = _non_overlapping_tick_labels(
            x_ticks,
            length=width,
            offset=offset,
            positions=x_positions,
            label_extent=_export_x_extent,
        )
        y_labeled = _non_overlapping_tick_labels(
            y_ticks,
            length=height,
            offset=offset,
            positions=y_positions,
            label_extent=_export_y_extent,
        )

        for tick in x_ticks:
            label = _format_bp(_scene_to_bp(tick, length=width, offset=offset))
            x = x_positions[tick]
            painter.drawLine(int(x), int(plot.top() - tick_len), int(x), int(plot.top()))
            if tick in x_labeled:
                painter.drawText(
                    int(x - metrics.horizontalAdvance(label) / 2),
                    int(plot.top() - tick_len - 2 * unit),
                    label,
                )

        for tick in y_ticks:
            label = _format_bp(_scene_to_bp(tick, length=height, offset=offset))
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
            pen = QPen(colour)
            pen.setWidth(max(1, int(0.4 * unit)))
            painter.setPen(pen)
            painter.drawLine(int(x), int(legend_y), int(x + swatch), int(legend_y))
            painter.setPen(QPen(QColor(70, 70, 70)))
            painter.drawText(
                int(x + swatch + 2 * unit),
                int(legend_y + metrics.ascent() / 2),
                label,
            )
            x += swatch + 4 * unit + metrics.horizontalAdvance(label)

    def render_to_pdf(
        self,
        path: Path,
        *,
        title: str = "",
        subtitle: str = "",
        dpi: int = 400,
    ) -> None:
        """Write the current view to a PDF with axes, labels and a legend.

        Exports the zoomed/panned viewport (not the full sequence square). The
        on-screen rulers are child widgets in the view margins, so axes are
        redrawn here using the same tick logic as the interactive rulers.
        """
        from PySide6.QtGui import QPageSize, QPdfWriter

        writer = QPdfWriter(str(path))
        writer.setPageSize(QPageSize(QPageSize.PageSizeId.A4))
        writer.setResolution(max(50, min(600, int(dpi))))
        writer.setTitle(title or "Dotplot")

        painter = QPainter(writer)
        try:
            self._paint_export_figure(
                painter,
                painter.viewport(),
                title=title,
                subtitle=subtitle,
            )
        finally:
            painter.end()

    def export_to_path(
        self,
        path: Path,
        *,
        dpi: int = 400,
        title: str = "",
        subtitle: str = "",
    ) -> None:
        """Export the current view as PDF, PNG, or SVG (by ``path`` suffix).

        PNG/SVG use the same A4 layout and paint path as PDF so figures match.
        """
        from PySide6.QtGui import QImage

        suffix = path.suffix.lower()
        dpi = max(50, min(600, int(dpi)))
        if suffix == ".pdf":
            self.render_to_pdf(path, title=title, subtitle=subtitle, dpi=dpi)
            return

        item = self._item
        prev_cache = None
        if item is not None:
            prev_cache = item.cacheMode()
            # DeviceCoordinateCache is for the on-screen viewport; painting to
            # a QImage/SVG with it active yields a sparse/stale plot.
            item.setCacheMode(QGraphicsItem.CacheMode.NoCache)

        try:
            if suffix == ".svg":
                try:
                    from PySide6.QtSvg import QSvgGenerator  # type: ignore
                except Exception as exc:  # noqa: BLE001
                    raise RuntimeError(
                        "SVG export needs the PySide6 QtSvg module."
                    ) from exc
                # Paint in A4 *points* (72 dpi logical page), not a 400-DPI pixel
                # grid. That keeps the SVG resolution-independent so Preview and
                # browsers can scale it without looking like an upscaled bitmap.
                page = QRect(0, 0, 595, 842)  # A4 at 72 pt/in
                if item is not None:
                    item.set_export_pen_px(_line_width_px(self._line_width_step))
                gen = QSvgGenerator()
                gen.setFileName(str(path))
                gen.setSize(page.size())
                gen.setViewBox(page)
                gen.setResolution(72)
                gen.setTitle(title or "Dotplot")
                painter = QPainter(gen)
                try:
                    self._paint_export_figure(
                        painter, page, title=title, subtitle=subtitle
                    )
                finally:
                    painter.end()
                return

            if suffix == ".png":
                page = self._a4_page_rect(dpi)
                # Same cosmetic width the PDF path uses on the A4 device.
                if item is not None:
                    item.set_export_pen_px(_line_width_px(self._line_width_step))
                img = QImage(page.size(), QImage.Format.Format_ARGB32)
                # Match QPdfWriter font metrics / physical size metadata.
                dpm = int(round(dpi / 0.0254))
                img.setDotsPerMeterX(dpm)
                img.setDotsPerMeterY(dpm)
                img.fill(Qt.GlobalColor.white)
                painter = QPainter(img)
                try:
                    # Match PDF: crisp hit lines (no AA on the plot layer).
                    painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
                    painter.setRenderHint(
                        QPainter.RenderHint.SmoothPixmapTransform, False
                    )
                    self._paint_export_figure(
                        painter, page, title=title, subtitle=subtitle
                    )
                finally:
                    painter.end()
                if not img.save(str(path), "PNG"):
                    raise RuntimeError(f"Could not save PNG: {path}")
                return
        finally:
            if item is not None:
                item.set_export_pen_px(None)
                if prev_cache is not None:
                    item.setCacheMode(prev_cache)

        raise RuntimeError(
            f"Unsupported export format: {suffix!r}. Use .svg, .png, or .pdf."
        )

    def wheelEvent(self, event) -> None:  # noqa: N802
        if self._item is None:
            return
        # Shift+scroll zooms (same accumulation as PNG preview). Plain scroll pans.
        if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            delta = wheel_zoom_delta_from_event(event)
            if delta == 0:
                return
            self._wheel_zoom_accum += delta
            steps = int(self._wheel_zoom_accum / PREVIEW_WHEEL_DELTA_PER_STEP)
            if steps != 0:
                self._wheel_zoom_accum -= steps * PREVIEW_WHEEL_DELTA_PER_STEP
                self._apply_zoom(
                    _ZOOM_STEP ** steps, self._viewport_event_pos(event)
                )
            event.accept()
            return

        pixel = event.pixelDelta()
        angle = event.angleDelta()
        if not pixel.isNull() and (pixel.x() != 0 or pixel.y() != 0):
            dx, dy = int(pixel.x()), int(pixel.y())
        elif angle.x() != 0 or angle.y() != 0:
            # ~30 px per classic 120-unit notch.
            dx = int(angle.x() // 4)
            dy = int(angle.y() // 4)
        else:
            return
        h = self.horizontalScrollBar()
        v = self.verticalScrollBar()
        h.setValue(h.value() - dx)
        v.setValue(v.value() - dy)
        self._clamp_scroll_to_plot()
        self._refresh_overlay()
        self._refresh_rulers()
        event.accept()

    def viewportEvent(self, event) -> bool:  # noqa: N802
        # Two-finger pinch on the trackpad (same path as PNG preview).
        if event.type() == QEvent.Type.NativeGesture and isinstance(
            event, QNativeGestureEvent
        ):
            if (
                event.gestureType() == Qt.NativeGestureType.ZoomNativeGesture
                and self._item is not None
            ):
                self._pinch_accum += float(event.value())
                steps = int(self._pinch_accum / _PINCH_DELTA_PER_STEP)
                if steps != 0:
                    self._pinch_accum -= steps * _PINCH_DELTA_PER_STEP
                    self._apply_zoom(
                        _ZOOM_STEP ** steps, self._viewport_event_pos(event)
                    )
                event.accept()
                return True
        return super().viewportEvent(event)

    def mousePressEvent(self, event) -> None:  # noqa: N802
        pos = self._viewport_event_pos(event)
        self._press_pos = pos
        if self._item is None:
            super().mousePressEvent(event)
            return

        mods = event.modifiers()
        shift = bool(mods & Qt.KeyboardModifier.ShiftModifier)

        if event.button() == Qt.MouseButton.MiddleButton:
            self._panning = True
            self._pan_origin = pos
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
            return

        if event.button() == Qt.MouseButton.LeftButton:
            if shift:
                # Shift+drag zooms (PNG preview behaviour).
                self._zoom_dragging = True
                self._zoom_drag_origin = QPointF(event.globalPosition())
                self._zoom_drag_accum = 0.0
                self.setCursor(Qt.CursorShape.SizeVerCursor)
                return
            # Click-drag draws a selection rectangle (viewport pixels).
            self._selection = None
            self._drag_start_vp = pos
            self._drag_now_vp = pos
            self._cursor_vp = pos
            self._cursor_pos = self.mapToScene(pos)
            self._refresh_overlay()
            return

        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        pos = self._viewport_event_pos(event)
        if self._zoom_dragging and self._zoom_drag_origin is not None:
            # Drag up → zoom in; drag down → zoom out (matches PNG preview).
            global_pos = QPointF(event.globalPosition())
            dy = self._zoom_drag_origin.y() - global_pos.y()
            self._zoom_drag_accum += dy
            self._zoom_drag_origin = global_pos
            steps = int(self._zoom_drag_accum / _DRAG_ZOOM_PIXELS_PER_STEP)
            if steps != 0:
                self._zoom_drag_accum -= steps * _DRAG_ZOOM_PIXELS_PER_STEP
                self._apply_zoom(_ZOOM_STEP ** steps, pos)
            return
        if self._panning and self._pan_origin is not None:
            delta = pos - self._pan_origin
            self._pan_origin = pos
            h = self.horizontalScrollBar()
            v = self.verticalScrollBar()
            h.setValue(h.value() - delta.x())
            v.setValue(v.value() - delta.y())
            self._clamp_scroll_to_plot()
            self._refresh_overlay()
            self._refresh_rulers()
            return
        if self._item is None:
            super().mouseMoveEvent(event)
            return
        self._cursor_vp = pos
        point = self.mapToScene(pos)
        self._cursor_pos = point
        if self._drag_start_vp is not None:
            # Keep selection drag in viewport pixels so the highlight tracks the
            # cursor even when zoomed; convert to scene only on release.
            self._drag_now_vp = pos
            self._refresh_overlay()
            self._refresh_rulers()
            x, y = int(point.x()), int(point.y())
            height = self.plot_height()
            if 0 <= x <= self._length and 0 <= y <= height:
                self.positionChanged.emit(*self._display_xy(x, y))
            return
        super().mouseMoveEvent(event)
        # Overlay only — do not viewport().update() (that re-paints all hits).
        self._refresh_overlay()
        self._refresh_rulers()
        x, y = int(point.x()), int(point.y())
        height = self.plot_height()
        if 0 <= x <= self._length and 0 <= y <= height:
            self.positionChanged.emit(*self._display_xy(x, y))

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if self._zoom_dragging:
            self._zoom_dragging = False
            self._zoom_drag_origin = None
            self._zoom_drag_accum = 0.0
            self.unsetCursor()
            self._press_pos = None
            return
        if self._panning:
            press = self._press_pos
            self._panning = False
            self._pan_origin = None
            self.unsetCursor()
            self._press_pos = None
            # Tiny motion → click clears the selection (no accidental clear on pan).
            if press is not None:
                moved = (self._viewport_event_pos(event) - press).manhattanLength()
                if moved <= _CLICK_SLOP:
                    self._selection = None
                    self._clear_drag()
                    self._refresh_overlay()
            return
        super().mouseReleaseEvent(event)
        if self._item is None or self._press_pos is None:
            self._clear_drag()
            return
        released = self._viewport_event_pos(event)
        moved = (released - self._press_pos).manhattanLength()
        self._press_pos = None

        if moved <= _CLICK_SLOP:
            # A plain click only clears the selection; it never navigates.
            self._clear_drag()
            self._selection = None
            self._refresh_overlay()
            return

        if self._drag_start_vp is not None and self._drag_now_vp is not None:
            # Convert the pixel rubber band into scene (bp) coordinates once.
            a = self.mapToScene(self._drag_start_vp)
            b = self.mapToScene(self._drag_now_vp)
            self._clear_drag()
            rect = QRectF(a, b).normalized()
            bounds = QRectF(
                0.0, 0.0, float(self._length), float(self.plot_height())
            )
            rect = _ensure_min_bp_rect(rect.intersected(bounds), bounds)
            if rect.width() >= _MIN_SELECTION_BP and rect.height() >= _MIN_SELECTION_BP:
                self._selection = rect
                self._refresh_overlay()
                self.rectSelected.emit(
                    int(rect.left()) + 1 + self._offset,
                    int(rect.right()) + self._offset,
                    int(rect.top()) + 1 + self._offset,
                    int(rect.bottom()) + self._offset,
                )
            return
        self._clear_drag()

    def leaveEvent(self, event) -> None:  # noqa: N802
        super().leaveEvent(event)
        self._cursor_pos = None
        self._cursor_vp = None
        self._refresh_overlay()
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

        # Parameter controls live on the right Display rail.
        self._rail = DotplotRail(self)
        self._word = self._rail.word
        self._reverse = self._rail.reverse
        self._line_width = self._rail.line_width
        self._line_width_label = self._rail.line_width_label
        self._identity = self._rail.identity
        self._min_len = self._rail.min_len
        self._seed_freq = self._rail.seed_freq
        self._region_start = self._rail.region_start
        self._region_size = self._rail.region_size
        self._region_start_label = self._rail.region_start_label
        self._region_size_label = self._rail.region_size_label
        self._pair_combo = self._rail.pair_combo
        self._pair_label = self._rail.pair_label
        self._sequence_combo = self._rail.sequence_combo
        self._sequence_label = self._rail.sequence_label
        self._sequence_prev_btn = self._rail.sequence_prev_btn
        self._sequence_next_btn = self._rail.sequence_next_btn
        self._filter_len = self._rail.filter_len
        self._filter_len_label = self._rail.filter_len_label
        self._word_label = self._rail.word_label
        self._identity_label = self._rail.identity_label
        self._min_len_label = self._rail.min_len_label
        self._seed_freq_label = self._rail.seed_freq_label
        self._compute_btn = self._rail.compute_btn
        self._cancel_btn = self._rail.cancel_btn
        self._zoom_out_btn = self._rail.zoom_out_btn
        self._zoom_in_btn = self._rail.zoom_in_btn
        self._reset_btn = self._rail.fit_btn
        self._zoom_select_btn = self._rail.select_btn
        self._pdf_btn = self._rail.export_btn
        self._sel_row = self._rail.sel_row
        self._sel_x = self._rail.sel_x
        self._sel_y = self._rail.sel_y
        self._view_x_btn = self._rail.view_x_btn
        self._view_y_btn = self._rail.view_y_btn

        # Ids for the Sequence picker; sequence text is loaded lazily.
        self._sequence_ids: list[str] = []
        self._sequence_cache: dict[str, str] = {}
        self._sequence_loader: Callable[[str], str] | None = None
        # Legacy alias used by engine-visibility checks (id list length).
        self._sequence_records: list[tuple[str, str]] = []
        self._ungap_for_plot = False
        # When ungapping for MSA: 0-based gapped index of each ungapped residue.
        self._coord_keep: list[int] | None = None

        self._thread: QThread | None = None
        self._worker: _FastgaWorker | _PafLoadWorker | _PafPairLoadWorker | None = None
        self._pending_key = None
        self._paf_load_path: Path | None = None
        self._pair_scan_key: tuple[str, str] | None = None
        # Last drag-box as 1-based inclusive (x1, x2, y1, y2), if any.
        self._selection_xy: tuple[int, int, int, int] | None = None
        # Status line without the appended x/y readout (alignment mode).
        self._status_base = ""
        self._status_plain = ""

        self._status = QLabel("", self)
        self._status.setTextFormat(Qt.TextFormat.RichText)
        self._status.setContentsMargins(0, 2, 6, 2)

        self._title = QLabel("", self)
        self._title.setContentsMargins(_RULER_LEFT, 2, 6, 2)
        title_font = QFont(self.font())
        title_font.setBold(True)
        self._title.setFont(title_font)
        self._title.hide()

        self._canvas = DotplotCanvas(self)
        self._canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        self._gallery = DotplotGallery(self)
        self._gallery.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        # Multi-seq / MSA default to Grid; Single is used for one-at-a-time.
        self._view_mode = "grid"
        self._gallery_return_sid: str | None = None

        single_page = QWidget(self)
        single_layout = QVBoxLayout(single_page)
        single_layout.setContentsMargins(0, 0, 0, 0)
        single_layout.setSpacing(0)
        title_row = QHBoxLayout()
        title_row.setContentsMargins(0, 0, 8, 0)
        title_row.setSpacing(8)
        title_row.addWidget(self._title, 0, Qt.AlignmentFlag.AlignVCenter)
        title_row.addWidget(self._status, 0, Qt.AlignmentFlag.AlignVCenter)
        title_row.addStretch(1)
        single_layout.addLayout(title_row)
        single_layout.addWidget(self._canvas, stretch=1)

        self._plot_stack = QStackedWidget(self)
        self._plot_stack.addWidget(single_page)  # 0 = Single
        self._plot_stack.addWidget(self._gallery)  # 1 = Grid
        self._plot_stack.setCurrentIndex(0)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._plot_stack, stretch=1)
        layout.addWidget(self._rail, 0)

        self._word.valueChanged.connect(self._recompute)
        self._reverse.toggled.connect(self._recompute)
        self._line_width.valueChanged.connect(self._apply_line_width)
        self._identity.valueChanged.connect(self._recompute)
        self._min_len.valueChanged.connect(self._recompute)
        self._seed_freq.valueChanged.connect(self._recompute)
        self._region_start.valueChanged.connect(self._recompute)
        self._region_size.valueChanged.connect(self._recompute)
        self._pair_combo.currentIndexChanged.connect(self._show_selected_pair)
        self._sequence_combo.currentIndexChanged.connect(
            self._on_sequence_choice_changed
        )
        self._filter_len.valueChanged.connect(self._recompute)
        self._reset_btn.clicked.connect(self._canvas.reset_zoom)
        self._zoom_select_btn.clicked.connect(self._zoom_to_selection)
        self._view_x_btn.clicked.connect(
            lambda: self._view_axis_in_annotation("x")
        )
        self._view_y_btn.clicked.connect(
            lambda: self._view_axis_in_annotation("y")
        )
        self._zoom_in_btn.clicked.connect(lambda: self._canvas.zoom_by(_ZOOM_STEP))
        self._zoom_out_btn.clicked.connect(
            lambda: self._canvas.zoom_by(1 / _ZOOM_STEP)
        )
        self._pdf_btn.clicked.connect(self._export_dotplot)
        self._rail.zoomSliderChanged.connect(self._on_zoom_slider_changed)
        self._rail.cancelClicked.connect(self._cancel_running)
        self._rail.computeClicked.connect(self._start_requested_compute)
        self._rail.viewModeChanged.connect(self._on_view_mode_changed)
        self._rail.galleryPageStep.connect(self._on_gallery_page_step)
        self._rail.galleryPageSelected.connect(self._on_gallery_page_selected)
        self._rail.galleryPageSizeChanged.connect(self._on_gallery_page_size_changed)
        self._rail.galleryThumbSizeChanged.connect(self._on_gallery_thumb_size_changed)
        self._rail.returnToGridClicked.connect(self._return_to_gallery)
        self._gallery.thumbClicked.connect(self._on_gallery_thumb_clicked)
        self._gallery.pageChanged.connect(self._sync_gallery_page_label)
        self._sync_engine_controls()
        self._canvas.positionChanged.connect(self._on_position)
        self._canvas.rectSelected.connect(self._on_rect_selected)
        self._canvas.zoomChanged.connect(self._sync_zoom_slider)

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
        self._sequence_ids = []
        self._sequence_cache = {}
        self._sequence_loader = None
        self._sequence_records = []
        self._ungap_for_plot = False
        self._coord_keep = None
        self._sequence_combo.blockSignals(True)
        self._sequence_combo.setModel(QStandardItemModel(self._sequence_combo))
        self._sequence_combo.blockSignals(False)
        self._too_long_prompt_active = False
        self._view_mode = "grid"
        self._gallery_return_sid = None
        self._rail.set_view_mode("grid")
        self._rail.set_return_to_grid_visible(False)
        self._plot_stack.setCurrentIndex(0)
        self._gallery.clear()
        self._canvas.set_axes([], [])
        self._sequence_id = ""
        self._canvas.clear_plot()
        self._restore_fasta_word_size_range()
        self._canvas.set_title("")
        self._title.clear()
        self._title.hide()
        self.clear_selection()
        self._set_status_text("")
        self._compute_btn.hide()
        self._sync_engine_controls()
        self._rail.set_selection_ranges()

    def refresh_sequence_summary(self) -> None:
        """No-op: Summary was removed from the Display rail."""
        return

    def _on_zoom_slider_changed(self, value: int) -> None:
        self._canvas.set_zoom_slider_value(int(value))

    def _sync_zoom_slider(self) -> None:
        self._rail.set_zoom_slider_value(self._canvas.zoom_slider_value())

    def set_record_choices(
        self,
        records: list[str] | list[tuple[str, str]],
        *,
        ungap: bool = False,
        preferred_id: str | None = None,
        sequence_loader: Callable[[str], str] | None = None,
    ) -> None:
        """Populate the Sequence picker for multi-FASTA / MSA Dotplot (Self).

        ``records`` may be ids only (preferred for large MSAs) or
        ``(id, sequence)`` pairs. Sequence text is loaded lazily via
        ``sequence_loader`` when needed so opening the tab with thousands of
        rows does not copy every alignment string into the combo.

        When ``ungap`` is True (MSA), gaps are stripped for plotting and a
        coordinate map is kept so selections can highlight the MSA.
        """
        ids: list[str] = []
        cache: dict[str, str] = {}
        for item in records:
            if isinstance(item, tuple):
                sid = str(item[0] or "").strip()
                if not sid:
                    continue
                ids.append(sid)
                cache[sid] = str(item[1] or "")
            else:
                sid = str(item or "").strip()
                if sid:
                    ids.append(sid)
        same_ids = ids == self._sequence_ids and bool(ungap) == self._ungap_for_plot
        self._sequence_ids = ids
        self._sequence_loader = sequence_loader
        # Keep any preloaded sequences; drop cache entries that left the list.
        if cache:
            self._sequence_cache = cache
        else:
            keep = set(ids)
            self._sequence_cache = {
                sid: seq for sid, seq in self._sequence_cache.items() if sid in keep
            }
        self._sequence_records = [(sid, self._sequence_cache.get(sid, "")) for sid in ids]
        self._ungap_for_plot = bool(ungap)
        if not same_ids:
            # Batch model fill — much faster than addItem() × N for large MSAs.
            model = QStandardItemModel(len(ids), 1, self._sequence_combo)
            for index, sid in enumerate(ids):
                item = QStandardItem(f"{index + 1}. {sid}")
                item.setData(sid, Qt.ItemDataRole.UserRole)
                model.setItem(index, 0, item)
            self._sequence_combo.blockSignals(True)
            self._sequence_combo.setModel(model)
            self._sequence_combo.blockSignals(False)
        want = str(preferred_id or "").strip()
        idx = self._sequence_combo.findData(want) if want else -1
        if idx < 0 and ids:
            idx = 0
        self._sequence_combo.blockSignals(True)
        if idx >= 0:
            self._sequence_combo.setCurrentIndex(idx)
        self._sequence_combo.blockSignals(False)
        self._sync_engine_controls()
        if ids:
            sid = str(self._sequence_combo.currentData() or ids[0])
            self.set_active_sequence(
                sid, self._sequence_for_id(sid), ungap=self._ungap_for_plot
            )
            if len(ids) > 1:
                # Multi-seq display is always Grid (no View dropdown).
                # Single is only via thumb click / Overview "Show self-dotplot".
                if not same_ids:
                    self._view_mode = "grid"
                    self._rail.set_view_mode("grid")
                    if self.isVisible():
                        self._show_gallery(reset_page=True)
                elif self._view_mode == "single" and self._gallery_return_sid:
                    # Keep an open Single detail view across a list resync.
                    self._rail.set_view_mode("single")
                    self._plot_stack.setCurrentIndex(0)
                    self._rail.set_return_to_grid_visible(True)
                else:
                    self._view_mode = "grid"
                    self._rail.set_view_mode("grid")
                    if self.isVisible():
                        self._show_gallery(reset_page=False)
            else:
                self._view_mode = "single"
                self._rail.set_view_mode("single")
                self._plot_stack.setCurrentIndex(0)
                self._rail.set_return_to_grid_visible(False)
        else:
            self.set_active_sequence("", "")
            self._view_mode = "single"
            self._rail.set_view_mode("single")
            self._plot_stack.setCurrentIndex(0)
            self._gallery.clear()
            self._rail.set_return_to_grid_visible(False)

    def _sequence_for_id(self, sequence_id: str) -> str:
        """Return sequence text for ``sequence_id``, loading it once if needed."""
        sid = str(sequence_id or "").strip()
        if not sid:
            return ""
        cached = self._sequence_cache.get(sid)
        if cached is not None:
            return cached
        if self._sequence_loader is None:
            return ""
        try:
            seq = str(self._sequence_loader(sid) or "")
        except (OSError, ValueError):
            seq = ""
        self._sequence_cache[sid] = seq
        return seq

    def _on_sequence_choice_changed(self, _index: int = -1) -> None:
        if self._mode == "alignment" or not self._sequence_ids:
            return
        sid = str(self._sequence_combo.currentData() or "").strip()
        if not sid:
            return
        self.set_active_sequence(
            sid, self._sequence_for_id(sid), ungap=self._ungap_for_plot
        )
        if self.isVisible():
            self.ensure_loaded()

    def _numbered_sequence_label(self, sequence_id: str) -> str:
        """Match the Sequence combo: ``N. id`` when multiple sequences are loaded."""
        sid = str(sequence_id or "").strip()
        if not sid:
            return ""
        ids = self._sequence_ids
        if len(ids) > 1 and sid in ids:
            return f"{ids.index(sid) + 1}. {sid}"
        return sid

    def set_active_sequence(
        self, sequence_id: str, sequence: str, *, ungap: bool | None = None
    ) -> None:
        """Use this record for Dotplot (Self) on multi-sequence FASTA files.

        When ``ungap`` is True, strip ``-``/``.`` for the plot and remember
        gapped coordinates for Overview highlighting.
        """
        sid = str(sequence_id or "").strip()
        seq = str(sequence or "")
        do_ungap = self._ungap_for_plot if ungap is None else bool(ungap)
        plot_seq = seq
        keep: list[int] | None = None
        if do_ungap and seq:
            from repeatviewer.backends.msa_coordinates import ungap_sequence

            plot_seq, keep = ungap_sequence(seq)
        if (
            sid == self._active_sequence_id
            and plot_seq == self._active_sequence
            and keep == self._coord_keep
        ):
            return
        self._active_sequence_id = sid
        self._active_sequence = plot_seq
        self._coord_keep = keep
        self._full_length = len(plot_seq)
        self._loaded_for = None
        self._too_long_prompt_active = False
        if sid and self._sequence_combo.count():
            idx = self._sequence_combo.findData(sid)
            if idx >= 0 and self._sequence_combo.currentIndex() != idx:
                self._sequence_combo.blockSignals(True)
                self._sequence_combo.setCurrentIndex(idx)
                self._sequence_combo.blockSignals(False)

    def _is_large(self) -> bool:
        if self._active_sequence:
            return len(self._active_sequence) > FASTGA_PREFERRED_ABOVE_BYTES
        try:
            return self._path.stat().st_size > FASTGA_PREFERRED_ABOVE_BYTES
        except OSError:
            return False

    def load_alignment(self, path: Path) -> None:
        """Display a PAF (or .1aln, converted) instead of computing anything.

        Parsing runs on a worker thread so large multi-FastGA PAFs do not freeze
        Preview. Huge files keep the longest hits in memory (FastGA ALNplot
        style); picking a Compare pair re-scans that pair. Stop cancels the read.
        """
        self._mode = "alignment"
        self._path = path.expanduser().resolve()
        self._loaded_for = None
        self._alignments = None
        self._pair_scan_key = None
        self._pair_combo.blockSignals(True)
        self._pair_combo.clear()
        self._pair_combo.blockSignals(False)
        # Show recorded run word size once the PAF loads (read-only here).
        self._word.blockSignals(True)
        self._word.setMinimum(0)
        self._word.setValue(0)
        self._word.setReadOnly(True)
        self._word.setButtonSymbols(
            self._word.ButtonSymbols.NoButtons
        )
        self._word.setToolTip(
            "Word size used when this PAF was computed (built-in Dotplot). "
            "Not edited here — the file is already computed. Use Min length "
            "to hide short hits. Older / FastGA PAFs may show 0 if the value "
            "was not recorded."
        )
        self._word.blockSignals(False)
        self._sync_engine_controls()
        self._set_status_text("Reading alignments…")
        self._start_paf_load(self._path)

    def _restore_fasta_word_size_range(self) -> None:
        """Put Word size back to compute-mode bounds after leaving a PAF."""
        self._word.blockSignals(True)
        self._word.setReadOnly(False)
        self._word.setButtonSymbols(self._word.ButtonSymbols.UpDownArrows)
        self._word.setMinimum(MIN_WORD_SIZE)
        if int(self._word.value()) < MIN_WORD_SIZE:
            self._word.setValue(DEFAULT_WORD_SIZE)
        self._word.setToolTip(
            "Length of the exact match used to place a dot. Larger = stricter "
            "and faster; smaller = more sensitive and noisier."
        )
        self._word.blockSignals(False)

    def _start_paf_load(self, path: Path) -> None:
        """Kick off (or replace) a background PAF parse for ``path``."""
        if self._worker is not None:
            self._worker.cancel()

        self._paf_load_path = path
        self._cancel_btn.show()
        thread = QThread(self)
        worker = _PafLoadWorker(path)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.progress.connect(self._set_status_text)
        worker.finished.connect(self._on_paf_load_finished)
        worker.failed.connect(self._on_paf_load_failed)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(self._on_thread_finished)
        self._thread = thread
        self._worker = worker
        thread.start()

    @Slot(object, object, object)
    def _on_paf_load_finished(self, alignments, source, request_path) -> None:
        if Path(request_path) != self._path:
            return
        self._apply_loaded_alignments(alignments, source)

    @Slot(str)
    def _on_paf_load_failed(self, message: str) -> None:
        if message.strip().lower() == "cancelled.":
            # A newer load may already be running; only clear the status when
            # nothing else replaced this request.
            if self._thread is None or not self._thread.isRunning():
                self._set_status_text("Cancelled.")
            return
        self._alignments = None
        self._sync_engine_controls()
        self._fail(message)

    def _apply_loaded_alignments(self, alignments, source) -> None:
        """Fill the Compare menu and plot after a successful background parse."""
        self._alignments = alignments
        self._pair_scan_key = None
        recorded = getattr(alignments, "word_size", None)
        self._word.blockSignals(True)
        if recorded is not None and int(recorded) > 0:
            self._word.setValue(int(recorded))
        else:
            self._word.setValue(0)
        self._word.blockSignals(False)
        menu_pairs, total_pairs = alignments.pairs_for_combo()
        offer_all = (
            total_pairs > 1 and alignments.supports_all_sequences_view()
        )
        self._pair_combo.blockSignals(True)
        self._pair_combo.clear()
        if offer_all:
            self._pair_combo.addItem("All sequences", ("__all__", "__all__"))
        for query, target in menu_pairs:
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
        notes: list[str] = []
        if source != self._path:
            self.statusMessage.emit(
                f"Converted {self._path.name} to {source.name}"
            )
        if total_pairs > len(menu_pairs):
            notes.append(
                f"Compare lists {len(menu_pairs)} of {total_pairs} pairs "
                "(longest first)."
            )
        if total_pairs > 1 and not offer_all:
            notes.append(
                f"{alignments.sequence_count()} sequences — pick a pair "
                f"(All-sequences view is limited to {MAX_ALL_SEQUENCES})."
            )
        if alignments.truncated:
            notes.append("Truncated: too many records.")
        if notes:
            self._set_status_text(" · ".join(notes))
        self._show_selected_pair()

    def _show_selected_pair(self) -> None:
        if self._alignments is None:
            return
        data = self._pair_combo.currentData()
        try:
            min_len = max(int(self._filter_len.value()), 0)
            if data is None or data[0] == "__all__":
                result, x_entries, y_entries = result_for_all(
                    self._alignments, min_length=min_len
                )
                zoom_hits = False
                pending_pair_scan: tuple[str, str] | None = None
            else:
                query, target = data[0], data[1]
                window = None
                # Region window is only meaningful for self comparisons (square).
                if query == target:
                    full = max(
                        self._alignments.query_lengths.get(query, 0),
                        self._alignments.target_lengths.get(target, 0),
                        1,
                    )
                    start = int(self._region_start.value())
                    size = int(self._region_size.value())
                    if start > 1 or size < full:
                        window = (start, min(start + size - 1, full))
                result = result_for_pair(
                    self._alignments,
                    query,
                    target,
                    window=window,
                    min_length=min_len,
                )
                # Label both axes with the sequence names (same style as
                # All-sequences), including when only two sequences are compared.
                x_entries = [
                    AxisEntry(
                        name=query,
                        offset=0,
                        length=max(int(result.sequence_length), 1),
                    )
                ]
                y_entries = [
                    AxisEntry(
                        name=target,
                        offset=0,
                        length=max(int(result.effective_target_length), 1),
                    )
                ]
                # Sparse hits on multi-Mb scaffolds vanish at Fit — frame them.
                zoom_hits = True
                # Huge PAFs only keep global longest hits in RAM. Re-scan for
                # this Compare pair in the background (FastGA ALNplot style).
                src = getattr(self._alignments, "source_path", None) or self._path
                pending_pair_scan = None
                if (
                    self._alignments.truncated
                    and src is not None
                    and Path(src).is_file()
                    and self._pair_scan_key != (query, target)
                ):
                    pending_pair_scan = (query, target)
        except Exception as exc:  # noqa: BLE001
            self._fail(f"Could not build the plot: {exc}")
            return
        self._canvas.set_axes(x_entries, y_entries)
        self._sync_engine_controls()
        self._apply_result(
            ("alignment", self._path, str(data)),
            result,
            zoom_to_hits=zoom_hits,
        )
        if pending_pair_scan is not None:
            q, t = pending_pair_scan
            src = getattr(self._alignments, "source_path", None) or self._path
            self._start_paf_pair_load(Path(src), q, t)

    def _start_paf_pair_load(self, path: Path, query: str, target: str) -> None:
        """Background re-scan of ``path`` for one query×target pair."""
        if self._worker is not None:
            self._worker.cancel()
        self._pair_scan_key = (query, target)
        self._set_status_text(f"Loading {query} vs {target}…")
        self._cancel_btn.show()
        thread = QThread(self)
        worker = _PafPairLoadWorker(path, query, target)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.progress.connect(self._set_status_text)
        worker.finished.connect(self._on_paf_pair_finished)
        worker.failed.connect(self._on_paf_pair_failed)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(self._on_thread_finished)
        self._thread = thread
        self._worker = worker
        thread.start()

    @Slot(object, object, object)
    def _on_paf_pair_finished(self, pair_aln, query, target) -> None:
        if self._pair_scan_key != (query, target):
            return
        if self._alignments is None:
            return
        # Merge pair hits into the cached alignments (replace same-pair rows).
        kept = [
            r
            for r in self._alignments.records
            if not (r[0] == query and r[5] == target)
        ]
        kept.extend(pair_aln.records)
        self._alignments.records = kept
        for name, n in pair_aln.query_lengths.items():
            self._alignments.query_lengths[name] = max(
                self._alignments.query_lengths.get(name, 0), int(n)
            )
        for name, n in pair_aln.target_lengths.items():
            self._alignments.target_lengths[name] = max(
                self._alignments.target_lengths.get(name, 0), int(n)
            )
        if pair_aln.pair_spans:
            self._alignments.pair_spans.update(pair_aln.pair_spans)
        # Build from the freshly scanned pair without triggering another scan.
        min_len = max(int(self._filter_len.value()), 0)
        try:
            result = result_for_pair(
                pair_aln, query, target, min_length=min_len
            )
        except Exception as exc:  # noqa: BLE001
            self._fail(f"Could not build the plot: {exc}")
            return
        x_entries = [
            AxisEntry(
                name=query,
                offset=0,
                length=max(int(result.sequence_length), 1),
            )
        ]
        y_entries = [
            AxisEntry(
                name=target,
                offset=0,
                length=max(int(result.effective_target_length), 1),
            )
        ]
        self._canvas.set_axes(x_entries, y_entries)
        self._sync_engine_controls()
        self._apply_result(
            ("alignment", self._path, (query, target)),
            result,
            zoom_to_hits=True,
        )

    @Slot(str)
    def _on_paf_pair_failed(self, message: str) -> None:
        failed_key = self._pair_scan_key
        self._pair_scan_key = None
        if message.strip().lower() == "cancelled.":
            if self._thread is None or not self._thread.isRunning():
                self._set_status_text("Cancelled.")
            return
        # Keep whatever in-memory plot we already showed; report the scan miss.
        detail = message.strip() or "Pair re-scan failed."
        if failed_key is not None:
            self._set_status_text(detail)
        else:
            self._fail(detail)

    def _use_fastga(self) -> bool:
        """Large sequences use FastGA (opt-in via Run FastGA), not the k-mer window."""
        return self._is_large() and fastga_available()

    def _use_region(self) -> bool:
        # Instant From/Window plotting is disabled for long sequences; users
        # are directed to FastGA instead.
        return False

    def _sync_engine_controls(self) -> None:
        if self._mode == "alignment":
            # PAF is precomputed; Word size / Min length only filter the plot.
            for widget in (
                self._reverse,
                self._identity_label, self._identity,
                self._min_len_label, self._min_len,
                self._seed_freq_label, self._seed_freq,
            ):
                widget.setVisible(False)
            multi = (
                self._alignments is not None
                and not self._alignments.is_single_pair()
            )
            self._pair_label.setVisible(multi)
            self._pair_combo.setVisible(multi)
            self._sequence_label.setVisible(False)
            self._sequence_combo.setVisible(False)
            self._sequence_prev_btn.setVisible(False)
            self._sequence_next_btn.setVisible(False)
            self._rail.set_view_mode_visible(False)
            if self._view_mode == "grid":
                self._view_mode = "single"
                self._rail.set_view_mode("single")
                self._plot_stack.setCurrentIndex(0)
            # Word size is shown in the title hits line, not as a Display field.
            self._word_label.setVisible(False)
            self._word.setVisible(False)
            self._filter_len_label.setVisible(True)
            self._filter_len.setVisible(True)
            # FastGA / PAF previews: no From/Window; x/y goes on the status line.
            for widget in (
                self._region_start_label, self._region_start,
                self._region_size_label, self._region_size,
            ):
                widget.setVisible(False)
            self._sel_row.hide()
            self._compute_btn.hide()
            self._line_width_label.setVisible(True)
            self._line_width.setVisible(True)
            self._rail.sync_control_visibility(view=True, params=True)
            self._rail.setVisible(True)
            return

        self._pair_label.setVisible(False)
        self._pair_combo.setVisible(False)
        multi_seq = self._mode != "alignment" and len(self._sequence_ids) > 1
        self._sequence_label.setVisible(multi_seq)
        self._sequence_combo.setVisible(multi_seq)
        self._sequence_prev_btn.setVisible(multi_seq)
        self._sequence_next_btn.setVisible(multi_seq)
        self._rail.set_view_mode_visible(multi_seq)
        if multi_seq:
            self._rail.sync_sequence_step_buttons()
            self._rail.set_gallery_nav_visible(self._view_mode == "grid")
            self._rail.set_return_to_grid_visible(
                self._view_mode == "single" and self._gallery_return_sid is not None
            )
        else:
            self._rail.set_return_to_grid_visible(False)
        self._filter_len_label.setVisible(False)
        self._filter_len.setVisible(False)
        grid_mode = multi_seq and self._view_mode == "grid"
        self._sel_row.setVisible(not grid_mode)
        # Status-only note for long sequences: hide zoom / Export / plot params.
        # Word size + Line width stay under Sequence whenever Display is shown.
        hide_bar = self._too_long_prompt_active and not grid_mode
        if hide_bar:
            self._compute_btn.hide()
            self._title.hide()
            for widget in (
                self._reverse,
                self._region_start_label, self._region_start,
                self._region_size_label, self._region_size,
                self._identity_label, self._identity,
                self._min_len_label, self._min_len,
                self._seed_freq_label, self._seed_freq,
            ):
                widget.setVisible(False)
            self._word_label.setVisible(True)
            self._word.setVisible(True)
            self._line_width_label.setVisible(True)
            self._line_width.setVisible(True)
            # Multi-record: keep Sequence + Word size + Line width.
            show_params = multi_seq
            self._rail.sync_control_visibility(view=False, params=show_params)
            self._rail.setVisible(show_params)
            return
        show_kmer = not self._is_large() or grid_mode
        self._reverse.setVisible(show_kmer)
        self._word_label.setVisible(True)
        self._word.setVisible(True)
        self._line_width_label.setVisible(True)
        self._line_width.setVisible(True)
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
        self._rail.sync_control_visibility(view=not grid_mode, params=True)
        self._rail.setVisible(True)

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
        if len(self._sequence_ids) > 1 and self._view_mode == "grid":
            self._show_gallery(reset_page=False)
            return
        if len(self._sequence_ids) > 1 and self._view_mode != "single":
            # First open of a multi-record file → Grid.
            self._show_gallery(reset_page=True)
            return
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

    def _on_view_mode_changed(self, mode: str) -> None:
        want = "grid" if str(mode).strip().lower() == "grid" else "single"
        if want == self._view_mode and (
            (want == "grid" and self._plot_stack.currentIndex() == 1)
            or (want == "single" and self._plot_stack.currentIndex() == 0)
        ):
            return
        if want == "grid" and len(self._sequence_ids) > 1:
            self._return_to_gallery()
        else:
            self._view_mode = "single"
            self._rail.set_view_mode("single")
            self._plot_stack.setCurrentIndex(0)
            self._sync_engine_controls()
            if self.isVisible():
                self.ensure_loaded()

    def _show_gallery(
        self,
        *,
        reset_page: bool = True,
        focus_sequence_id: str | None = None,
    ) -> None:
        """Show the paged multi-sequence thumbnail grid."""
        if len(self._sequence_ids) <= 1:
            self._view_mode = "single"
            self._rail.set_view_mode("single")
            self._plot_stack.setCurrentIndex(0)
            self._rail.set_return_to_grid_visible(False)
            return
        self._view_mode = "grid"
        self._rail.set_view_mode("grid")
        self._rail.set_return_to_grid_visible(False)
        self._plot_stack.setCurrentIndex(1)
        self._title.hide()
        self._set_status_text("")
        # Keep Grid pens in sync with the shared Line width control.
        self._gallery.set_line_width_step(int(self._line_width.value()))
        focus = str(focus_sequence_id or self._gallery_return_sid or "").strip() or None
        self._gallery.set_sequences(
            self._sequence_ids,
            sequence_loader=self._sequence_for_id,
            ungap=self._ungap_for_plot,
            word_size=int(self._word.value()),
            include_reverse=self._reverse.isChecked(),
            page_size=int(self._rail.page_size.value()),
            thumb_size=int(self._rail.thumb_size.value()),
            reset_page=reset_page,
            focus_sequence_id=None if reset_page else focus,
        )
        self._gallery_return_sid = None
        self._sync_gallery_page_label()
        self._sync_engine_controls()

    def _return_to_gallery(self) -> None:
        """Return to Grid, keeping page and focusing the last opened sequence."""
        focus = str(
            self._gallery_return_sid or self._active_sequence_id or ""
        ).strip() or None
        self._show_gallery(reset_page=False, focus_sequence_id=focus)

    def _on_gallery_page_step(self, delta: int) -> None:
        if self._view_mode != "grid":
            return
        self._gallery.set_page(self._gallery.page() + int(delta))

    def _on_gallery_page_selected(self, page: int) -> None:
        if self._view_mode != "grid":
            return
        self._gallery.set_page(int(page))

    def _on_gallery_page_size_changed(self, value: int) -> None:
        if self._view_mode != "grid":
            return
        self._gallery.set_page_size(int(value))
        self._sync_gallery_page_label()

    def _on_gallery_thumb_size_changed(self, value: int) -> None:
        if self._view_mode != "grid":
            return
        self._gallery.set_thumb_size(int(value))

    def _sync_gallery_page_label(self, _page: int = 0) -> None:
        self._rail.set_gallery_page_choices(
            page=self._gallery.page(),
            page_size=self._gallery.page_size(),
            total=len(self._sequence_ids),
        )

    def open_single_for_sequence(self, sequence_id: str) -> None:
        """Show one sequence in Single view (from Grid thumb or Overview menu)."""
        sid = str(sequence_id or "").strip()
        if not sid:
            return
        if sid not in self._sequence_ids and self._sequence_combo.findData(sid) < 0:
            return
        self._gallery_return_sid = sid
        self._view_mode = "single"
        self._rail.set_view_mode("single")
        self._plot_stack.setCurrentIndex(0)
        self.set_active_sequence(
            sid, self._sequence_for_id(sid), ungap=self._ungap_for_plot
        )
        self._sync_engine_controls()
        self._loaded_for = None
        self.ensure_loaded()

    def _on_gallery_thumb_clicked(self, sequence_id: str) -> None:
        """Leave Grid and open the clicked sequence as a full Single plot."""
        self.open_single_for_sequence(sequence_id)

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
        width = int(step)
        self._canvas.set_line_width_step(width)
        self._gallery.set_line_width_step(width)

    def _recompute(self) -> None:
        if self._mode == "alignment":
            if self._alignments is not None:
                self._show_selected_pair()
            return
        if self._view_mode == "grid" and len(self._sequence_ids) > 1:
            self._gallery.set_compute_params(
                word_size=int(self._word.value()),
                include_reverse=self._reverse.isChecked(),
            )
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
            self._fail(str(exc), key=key)
            return
        except Exception as exc:  # noqa: BLE001
            self._fail(f"Could not plot region: {exc}", key=key)
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
            self._fail(message, key=key)
            return
        except Exception as exc:  # noqa: BLE001
            self._fail(f"Could not build dotplot: {exc}", key=key)
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
        finished = self.sender()
        if finished is not None and finished is not self._thread:
            # A superseded load (user opened another file) — just clean up.
            finished.deleteLater()
            return
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

    def _fail(self, message: str, *, key=None) -> None:
        self._canvas.clear_plot()
        # Keep the attempted settings (not None). Clearing to None made
        # _recompute ignore later word-size / reverse edits until the user
        # switched sequences — the error stayed on screen while the spin
        # already showed the new value.
        if key is not None:
            self._loaded_for = key
        elif self._pending_key is not None:
            self._loaded_for = self._pending_key
        else:
            self._loaded_for = (
                "failed",
                self._path,
                self._active_sequence_id or "",
            )
        self._clear_selection_chrome()
        self._set_status_text(self._condense_error(message), tooltip=message)
        self.statusMessage.emit(message)

    def _apply_result(
        self, key, result: DotplotResult, *, zoom_to_hits: bool = False
    ) -> None:
        previous_length = self._sequence_length
        had_plot = self._loaded_for is not None
        self._loaded_for = key
        self._too_long_prompt_active = False
        self._sync_engine_controls()
        self._sequence_id = result.sequence_id
        self._sequence_length = result.sequence_length
        name = self._numbered_sequence_label(
            result.sequence_id or self._path.stem
        )
        self._title.setText(name)
        self._title.setVisible(bool(name))
        self._canvas.set_title(name)
        # Keep pan/zoom when recomputing the same-length sequence (word size /
        # reverse complement). New files / region sizes still fit to view.
        preserve = (
            had_plot
            and previous_length == result.sequence_length
            and not zoom_to_hits
        )
        self._canvas.show_result(result, preserve_view=preserve)
        if zoom_to_hits and result.segments:
            # Sparse hits on multi-Mb scaffolds vanish at Fit — frame them.
            self._canvas.zoom_to_hits(result)
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
        # Title line reads like: ``4 x 4 sequences (56,910 hits · word 6)``.
        paren_bits = [f"<b>{result.segment_count:,}</b> hits"]
        plain_bits = [f"{result.segment_count:,} hits"]
        word = int(result.word_size or 0)
        if word <= 0 and self._mode == "alignment":
            recorded = getattr(self._alignments, "word_size", None)
            if recorded is not None:
                word = int(recorded)
            else:
                word = int(self._word.value())
        if word > 0:
            paren_bits.append(f"word size {word}")
            plain_bits.append(f"word size {word}")
        self._set_status_text(
            f"({' · '.join(paren_bits)}){suffix}",
            plain=(
                f"{where_plain} · {engine_plain} · "
                f"{' · '.join(plain_bits)}{suffix}"
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

    def _source_to_plot_range(self, start: int, end: int) -> tuple[int, int]:
        """Map Overview/MSA (source) coords onto plot coords (ungapped when MSA)."""
        lo, hi = (int(start), int(end)) if start <= end else (int(end), int(start))
        keep = self._coord_keep
        if not keep:
            return lo, hi
        from repeatviewer.backends.msa_coordinates import gapped_to_ungapped_1based

        u1 = gapped_to_ungapped_1based(lo, keep) or 1
        u2 = gapped_to_ungapped_1based(hi, keep) or u1
        if u2 < u1:
            u1, u2 = u2, u1
        return int(u1), int(u2)

    def _plot_to_source_range(self, start: int, end: int) -> tuple[int, int]:
        """Map plot coords back to Overview/MSA (gapped) coordinates."""
        lo, hi = (int(start), int(end)) if start <= end else (int(end), int(start))
        keep = self._coord_keep
        if not keep:
            return lo, hi
        from repeatviewer.backends.msa_coordinates import ungapped_to_gapped_1based

        g1 = ungapped_to_gapped_1based(lo, keep)
        g2 = ungapped_to_gapped_1based(hi, keep)
        if g2 < g1:
            g1, g2 = g2, g1
        return int(g1), int(g2)

    def show_selection(self, start: int, end: int) -> None:
        """Mirror a 1-D Annotation/Overview selection onto the plot as an x/y box.

        ``start``/``end`` are source coordinates (gapped MSA columns when ungapping).
        """
        lo, hi = self._source_to_plot_range(start, end)
        self._canvas.set_range(lo, hi)
        # Same range on both axes — show x and y explicitly (not "Selected N bp").
        self._show_selection_status(lo, hi, lo, hi)

    def clear_selection(self) -> None:
        self._canvas.clear_selection()
        self._clear_selection_chrome()

    def _clear_selection_chrome(self) -> None:
        self._selection_xy = None
        self._rail.set_selection_ranges()
        if self._mode == "alignment":
            self._sel_row.hide()
            if self._status_base:
                self._status.setTextFormat(Qt.TextFormat.RichText)
                self._status.setText(self._status_base)
        else:
            self._sel_row.show()

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
        self._rail.set_selection_ranges(
            x=(int(x1), int(x2)),
            y=(int(y1), int(y2)),
        )
        self._sel_row.show()

    def _zoom_to_selection(self) -> None:
        if self._canvas.zoom_to_selection():
            return
        self._set_status_text(
            "Select a region on the plot first (drag a rectangle)."
        )

    def _view_axis_in_annotation(self, axis: str) -> None:
        if self._mode == "alignment":
            return
        coords = self._selection_xy or self._canvas.selection_sequence_coords()
        if coords is None:
            self._set_status_text(
            "Select a region on the plot first (drag a rectangle)."
        )
            return
        x1, x2, y1, y2 = coords
        if axis == "y":
            start, end = (y1, y2) if y1 <= y2 else (y2, y1)
        else:
            start, end = (x1, x2) if x1 <= x2 else (x2, x1)
        start, end = self._plot_to_source_range(int(start), int(end))
        self.viewInAnnotationRequested.emit(int(start), int(end))

    def _export_dotplot(self) -> None:
        from PySide6.QtWidgets import QDialog, QMessageBox

        from repeatviewer.views.export_annotation_view_dialog import (
            ExportAnnotationViewDialog,
        )

        if self._sequence_length <= 0:
            QMessageBox.information(
                self, "Export dotplot", "There is no plot to export yet."
            )
            return

        dialog = ExportAnnotationViewDialog(
            source_path=self._path,
            title="Export dotplot",
            name_tag="dotplot",
            default_format="pdf",
            parent=self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        out = dialog.chosen_path()
        if out is None:
            return

        try:
            self._canvas.export_to_path(
                out,
                dpi=dialog.dpi(),
                title=self._title.text() or self._path.stem,
                # Status line already has range, word size, and segment count.
                subtitle=self._status_caption(),
            )
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(
                self, "Export dotplot", f"Could not write {out.name}:\n{exc}"
            )
            return

        self.statusMessage.emit(f"Exported {out.name}")
        self._open_saved_file(out)

    def _open_saved_file(self, path: Path) -> None:
        """Refresh the folder list and open the PDF in the preview pane."""
        window = self.window()
        refresh = getattr(window, "_refresh_folders_after_job", None)
        if callable(refresh):
            try:
                refresh([path], [self._path])
            except Exception:
                pass
        opener = getattr(window, "_auto_open_outputs", None)
        if callable(opener):
            try:
                opener([path])
                return
            except Exception:
                pass
        from PySide6.QtCore import QUrl
        from PySide6.QtGui import QDesktopServices

        QDesktopServices.openUrl(QUrl.fromLocalFile(str(path)))

    def sequence_id(self) -> str:
        return self._active_sequence_id or self._sequence_id

    def _on_rect_selected(self, x1: int, x2: int, y1: int, y2: int) -> None:
        self._show_selection_status(x1, x2, y1, y2)
        if self._mode != "alignment":
            sx1, sx2 = self._plot_to_source_range(x1, x2)
            sy1, sy2 = self._plot_to_source_range(y1, y2)
            self.rangeSelected.emit(sx1, sx2, sy1, sy2)

    def _on_position(self, x: int, y: int) -> None:
        # Live coords are on the axis rulers; skip status-bar spam (was a lag source).
        del x, y