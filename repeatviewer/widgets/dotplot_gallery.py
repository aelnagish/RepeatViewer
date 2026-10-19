"""Paged grid of self-dotplot thumbnails for multi-sequence Dotplot (Self)."""

from __future__ import annotations

from PySide6.QtCore import QPointF, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QPainter, QPen
from PySide6.QtWidgets import (
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QScrollArea,
    QSizePolicy,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.dotplot import (
    STRAND_REVERSE,
    DotplotError,
    DotplotResult,
    compute_self_dotplot,
)
from repeatviewer.backends.fastga_dotplot import FASTGA_PREFERRED_ABOVE_BYTES

_FORWARD = QColor(40, 90, 180)
_REVERSE = QColor(200, 60, 40)
_DIAGONAL = QColor(45, 45, 45)
_BORDER = QColor(200, 200, 200)
_BG = QColor(255, 255, 255)
_LABEL = QColor(40, 40, 40)

_DEFAULT_PAGE_SIZE = 20
_MIN_PAGE_SIZE = 4
_MAX_PAGE_SIZE = 100
_DEFAULT_THUMB_SIZE = 140
_MIN_THUMB_SIZE = 80
_MAX_THUMB_SIZE = 280
_THUMB_LABEL_H = 22
_GRID_GAP = 8
_GRID_MARGIN = 4
# Match DotplotPanel Single-view Line width mapping (step 1–15 → 1–4.5 px).
_MIN_LINE_WIDTH_STEP = 1
_MAX_LINE_WIDTH_STEP = 15
_DEFAULT_LINE_WIDTH_STEP = 1
_LINE_WIDTH_PX_MIN = 1.0
_LINE_WIDTH_PX_MAX = 4.5


def _line_width_px(step: int) -> float:
    step = max(_MIN_LINE_WIDTH_STEP, min(_MAX_LINE_WIDTH_STEP, int(step)))
    if _MAX_LINE_WIDTH_STEP <= _MIN_LINE_WIDTH_STEP:
        return _LINE_WIDTH_PX_MIN
    t = (step - _MIN_LINE_WIDTH_STEP) / (
        _MAX_LINE_WIDTH_STEP - _MIN_LINE_WIDTH_STEP
    )
    return _LINE_WIDTH_PX_MIN + t * (_LINE_WIDTH_PX_MAX - _LINE_WIDTH_PX_MIN)


class DotplotThumb(QFrame):
    """One compact self-dotplot; click opens the full single-sequence view."""

    clicked = Signal(str)  # sequence id

    def __init__(
        self,
        sequence_id: str,
        index: int,
        *,
        plot_size: int = _DEFAULT_THUMB_SIZE,
        line_width_step: int = _DEFAULT_LINE_WIDTH_STEP,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._sequence_id = str(sequence_id or "")
        self._index = int(index)
        self._plot_size = int(plot_size)
        self._line_width_step = max(
            _MIN_LINE_WIDTH_STEP,
            min(_MAX_LINE_WIDTH_STEP, int(line_width_step)),
        )
        self._result: DotplotResult | None = None
        self._error = ""
        self._busy = False
        self.setObjectName("DotplotThumb")
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        # Fixed square plot area — Plot size must change both X and Y.
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self.setFrameShape(QFrame.Shape.StyledPanel)
        self.setStyleSheet(
            "QFrame#DotplotThumb {"
            "  background: palette(base);"
            "  border: 1px solid palette(mid);"
            "  border-radius: 4px;"
            "}"
            "QFrame#DotplotThumb:hover {"
            "  border: 1px solid palette(highlight);"
            "}"
        )
        self.setToolTip(self._sequence_id)
        self._apply_size()

    def sequence_id(self) -> str:
        return self._sequence_id

    def index(self) -> int:
        return int(self._index)

    def set_plot_size(self, plot_size: int) -> None:
        size = max(_MIN_THUMB_SIZE, min(_MAX_THUMB_SIZE, int(plot_size)))
        if size == self._plot_size:
            return
        self._plot_size = size
        self._apply_size()
        self.update()

    def set_line_width_step(self, step: int) -> None:
        step = max(_MIN_LINE_WIDTH_STEP, min(_MAX_LINE_WIDTH_STEP, int(step)))
        if step == self._line_width_step:
            return
        self._line_width_step = step
        self.update()

    def line_width_step(self) -> int:
        return int(self._line_width_step)

    def _apply_size(self) -> None:
        # Square plot face (side = plot_size) + label strip under it.
        side = max(_MIN_THUMB_SIZE, int(self._plot_size))
        w = side + 8
        h = side + _THUMB_LABEL_H + 10
        self.setFixedSize(w, h)

    def set_busy(self, busy: bool) -> None:
        self._busy = bool(busy)
        if busy:
            self._result = None
            self._error = ""
        self.update()

    def set_result(self, result: DotplotResult | None, *, error: str = "") -> None:
        self._busy = False
        self._result = result
        self._error = str(error or "")
        self.update()

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton and self._sequence_id:
            self.clicked.emit(self._sequence_id)
        super().mousePressEvent(event)

    def paintEvent(self, event) -> None:  # noqa: N802
        super().paintEvent(event)
        painter = QPainter(self)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        plot = self.rect().adjusted(4, 4, -4, -(_THUMB_LABEL_H + 2))
        painter.fillRect(plot, _BG)
        painter.setPen(QPen(_BORDER))
        painter.drawRect(plot.adjusted(0, 0, -1, -1))

        label = f"{self._index + 1}. {self._sequence_id}"
        metrics = painter.fontMetrics()
        elided = metrics.elidedText(
            label, Qt.TextElideMode.ElideMiddle, max(20, plot.width())
        )
        painter.setPen(_LABEL)
        painter.drawText(
            4,
            plot.bottom() + 2,
            plot.width(),
            _THUMB_LABEL_H,
            int(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter),
            elided,
        )

        if self._busy:
            painter.setPen(_LABEL)
            painter.drawText(plot, int(Qt.AlignmentFlag.AlignCenter), "…")
            return
        if self._error:
            painter.setPen(QColor(140, 60, 60))
            font = QFont(painter.font())
            font.setPointSize(max(8, font.pointSize() - 2))
            painter.setFont(font)
            painter.drawText(
                plot.adjusted(4, 4, -4, -4),
                int(Qt.AlignmentFlag.AlignCenter | Qt.TextFlag.TextWordWrap),
                self._error,
            )
            return
        if self._result is None or self._result.sequence_length <= 0:
            painter.setPen(_LABEL)
            painter.drawText(plot, int(Qt.AlignmentFlag.AlignCenter), "—")
            return

        self._paint_segments(
            painter, plot, self._result, line_width_step=self._line_width_step
        )

    @staticmethod
    def _paint_segments(
        painter: QPainter,
        plot,
        result: DotplotResult,
        *,
        line_width_step: int = _DEFAULT_LINE_WIDTH_STEP,
    ) -> None:
        n = max(1, int(result.sequence_length))
        pad = 2
        side = max(1, min(plot.width(), plot.height()) - 2 * pad)
        origin_x = plot.x() + pad + max(0, (plot.width() - 2 * pad - side) // 2)
        origin_y = plot.y() + pad + max(0, (plot.height() - 2 * pad - side) // 2)
        scale = side / float(n)
        width = _line_width_px(line_width_step)

        def map_pt(x: float, y: float) -> QPointF:
            return QPointF(origin_x + x * scale, origin_y + y * scale)

        pen = QPen(_DIAGONAL)
        pen.setCosmetic(True)
        pen.setWidthF(width)
        painter.setPen(pen)
        painter.drawLine(map_pt(0, 0), map_pt(n, n))

        fwd = QPen(_FORWARD)
        fwd.setCosmetic(True)
        fwd.setWidthF(width)
        rev = QPen(_REVERSE)
        rev.setCosmetic(True)
        rev.setWidthF(width)
        segments = result.segments
        if len(segments) > 4000:
            step = max(1, len(segments) // 4000)
            segments = segments[::step]
        for seg in segments:
            x1, y1 = float(seg.x), float(seg.y)
            x2, y2 = float(seg.end_x), float(seg.end_y)
            if seg.strand == STRAND_REVERSE:
                painter.setPen(rev)
                painter.drawLine(map_pt(x1, y1), map_pt(x2, y2))
            else:
                painter.setPen(fwd)
                painter.drawLine(map_pt(x1, y1), map_pt(x2, y2))
                painter.drawLine(map_pt(y1, x1), map_pt(y2, x2))


class DotplotGallery(QWidget):
    """Scrollable grid of thumbs with page navigation."""

    thumbClicked = Signal(str)
    pageChanged = Signal(int)  # 0-based page index

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._ids: list[str] = []
        self._page = 0
        self._page_size = _DEFAULT_PAGE_SIZE
        self._thumb_size = _DEFAULT_THUMB_SIZE
        self._line_width_step = _DEFAULT_LINE_WIDTH_STEP
        self._thumbs: list[DotplotThumb] = []
        self._cache: dict[tuple, DotplotResult | str] = {}
        self._pending: list[str] = []
        self._word = 10
        self._reverse = True
        self._loader = None
        self._ungap = False
        self._timer = QTimer(self)
        self._timer.setInterval(0)
        self._timer.timeout.connect(self._compute_next)
        self._layout_cols = 0
        self._layout_side = 0
        self._resize_timer = QTimer(self)
        self._resize_timer.setSingleShot(True)
        self._resize_timer.setInterval(80)
        self._resize_timer.timeout.connect(self._on_resize_relayout)

        root = QVBoxLayout(self)
        root.setContentsMargins(6, 4, 6, 4)
        root.setSpacing(4)

        nav = QHBoxLayout()
        nav.setContentsMargins(0, 0, 0, 0)
        nav.setSpacing(6)
        self._prev_btn = QToolButton(self)
        self._prev_btn.setText("◀")
        self._prev_btn.setAutoRaise(True)
        self._prev_btn.setToolTip("Previous page")
        self._prev_btn.clicked.connect(lambda: self.set_page(self._page - 1))
        self._next_btn = QToolButton(self)
        self._next_btn.setText("▶")
        self._next_btn.setAutoRaise(True)
        self._next_btn.setToolTip("Next page")
        self._next_btn.clicked.connect(lambda: self.set_page(self._page + 1))
        self._page_label = QLabel("", self)
        nav.addWidget(self._prev_btn, 0)
        nav.addWidget(self._next_btn, 0)
        nav.addWidget(self._page_label, 0)
        nav.addStretch(1)
        hint = QLabel("Click a plot to open it full size", self)
        hint.setStyleSheet("color: palette(mid);")
        nav.addWidget(hint, 0)
        root.addLayout(nav)

        self._scroll = QScrollArea(self)
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self._host = QWidget(self._scroll)
        self._host.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred
        )
        self._grid = QGridLayout(self._host)
        self._grid.setContentsMargins(
            _GRID_MARGIN, _GRID_MARGIN, _GRID_MARGIN, _GRID_MARGIN
        )
        self._grid.setSpacing(_GRID_GAP)
        # Left-align square tiles; column count fills the available width.
        self._grid.setAlignment(
            Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft
        )
        self._scroll.setWidget(self._host)
        root.addWidget(self._scroll, 1)

    @staticmethod
    def default_page_size() -> int:
        return _DEFAULT_PAGE_SIZE

    @staticmethod
    def max_page_size() -> int:
        return _MAX_PAGE_SIZE

    @staticmethod
    def default_thumb_size() -> int:
        return _DEFAULT_THUMB_SIZE

    @staticmethod
    def clamp_page_size(value: int) -> int:
        return max(_MIN_PAGE_SIZE, min(_MAX_PAGE_SIZE, int(value)))

    @staticmethod
    def clamp_thumb_size(value: int) -> int:
        return max(_MIN_THUMB_SIZE, min(_MAX_THUMB_SIZE, int(value)))

    def page_size(self) -> int:
        return int(self._page_size)

    def thumb_size(self) -> int:
        return int(self._thumb_size)

    def page(self) -> int:
        return int(self._page)

    def page_count(self) -> int:
        if not self._ids:
            return 1
        return max(1, (len(self._ids) + self._page_size - 1) // self._page_size)

    def set_sequences(
        self,
        ids: list[str],
        *,
        sequence_loader,
        ungap: bool = False,
        word_size: int = 10,
        include_reverse: bool = True,
        page_size: int | None = None,
        thumb_size: int | None = None,
        reset_page: bool = True,
        focus_sequence_id: str | None = None,
    ) -> None:
        new_ids = [str(x) for x in ids if str(x).strip()]
        new_page_size = (
            self.clamp_page_size(page_size)
            if page_size is not None
            else self._page_size
        )
        new_thumb_size = (
            self.clamp_thumb_size(thumb_size)
            if thumb_size is not None
            else self._thumb_size
        )
        word = int(word_size)
        reverse = bool(include_reverse)
        ungap_flag = bool(ungap)
        focus = str(focus_sequence_id or "").strip()
        params_same = (
            new_ids == self._ids
            and ungap_flag == self._ungap
            and word == self._word
            and reverse == self._reverse
            and new_page_size == self._page_size
            and new_thumb_size == self._thumb_size
        )
        # Tab re-entry with the same list: keep page and scroll; avoid rebuild.
        if params_same and not reset_page and self._thumbs and not focus:
            self._loader = sequence_loader
            return
        scroll_y = (
            int(self._scroll.verticalScrollBar().value())
            if params_same and not reset_page
            else 0
        )
        self._ids = new_ids
        self._loader = sequence_loader
        self._ungap = ungap_flag
        self._word = word
        self._reverse = reverse
        self._page_size = new_page_size
        self._thumb_size = new_thumb_size
        if not params_same:
            self._cache.clear()
        if reset_page:
            self._page = 0
        else:
            self._page = max(0, min(self._page, self.page_count() - 1))
        if focus and focus in self._ids and not reset_page:
            idx = self._ids.index(focus)
            self._page = idx // self._page_size
        self._rebuild_page()
        if focus:
            QTimer.singleShot(0, lambda: self.ensure_sequence_visible(focus))
        elif scroll_y > 0:

            def _restore_scroll() -> None:
                self._scroll.verticalScrollBar().setValue(scroll_y)

            QTimer.singleShot(0, _restore_scroll)

    def set_compute_params(self, *, word_size: int, include_reverse: bool) -> None:
        word = int(word_size)
        reverse = bool(include_reverse)
        if word == self._word and reverse == self._reverse:
            return
        self._word = word
        self._reverse = reverse
        self._cache.clear()
        self._rebuild_page()

    def set_page_size(self, page_size: int) -> None:
        size = self.clamp_page_size(page_size)
        if size == self._page_size:
            return
        first = self._page * self._page_size
        self._page_size = size
        self._page = first // self._page_size if self._ids else 0
        self._rebuild_page()
        self.pageChanged.emit(self._page)

    def set_thumb_size(self, thumb_size: int) -> None:
        size = self.clamp_thumb_size(thumb_size)
        if size == self._thumb_size:
            return
        self._thumb_size = size
        self._rebuild_page()

    def set_line_width_step(self, step: int) -> None:
        """Update match/diagonal thickness on all visible Grid thumbs."""
        step = max(_MIN_LINE_WIDTH_STEP, min(_MAX_LINE_WIDTH_STEP, int(step)))
        if step == self._line_width_step:
            return
        self._line_width_step = step
        for thumb in self._thumbs:
            thumb.set_line_width_step(step)

    def line_width_step(self) -> int:
        return int(self._line_width_step)

    def set_page(self, page: int) -> None:
        page = max(0, min(int(page), self.page_count() - 1))
        if page == self._page and self._thumbs:
            return
        self._page = page
        self._rebuild_page()
        self.pageChanged.emit(self._page)

    def ensure_sequence_visible(self, sequence_id: str) -> None:
        """Scroll so ``sequence_id`` is on the current page viewport if present."""
        sid = str(sequence_id or "").strip()
        if not sid:
            return
        if sid in self._ids:
            want_page = self._ids.index(sid) // max(1, self._page_size)
            if want_page != self._page:
                self.set_page(want_page)
        thumb = next((t for t in self._thumbs if t.sequence_id() == sid), None)
        if thumb is None:
            return
        self._scroll.ensureWidgetVisible(thumb, 12, 12)

    def clear(self) -> None:
        self._timer.stop()
        self._pending = []
        self._ids = []
        self._cache.clear()
        self._clear_thumbs()
        self._page_label.setText("")
        self._prev_btn.setEnabled(False)
        self._next_btn.setEnabled(False)

    def _clear_thumbs(self) -> None:
        while self._grid.count():
            item = self._grid.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        self._thumbs = []

    def _page_slice(self) -> list[tuple[int, str]]:
        start = self._page * self._page_size
        end = min(len(self._ids), start + self._page_size)
        return [(i, self._ids[i]) for i in range(start, end)]

    def _result_matches_params(self, result: DotplotResult | None) -> bool:
        """True when a prior/cached result was computed with current params."""
        if result is None:
            return False
        return (
            int(result.word_size) == int(self._word)
            and bool(result.include_reverse) == bool(self._reverse)
        )

    def _rebuild_page(self) -> None:
        self._timer.stop()
        self._pending = []
        # Preserve results already on thumbs only when word/reverse still match
        # (resize / page layout). After Word size changes, stale priors must
        # not be reused — that left Grid plots unchanged while Single updated.
        prior_results = {
            t.sequence_id(): (t._result, t._error)
            for t in self._thumbs
            if t.sequence_id()
        }
        self._clear_thumbs()
        entries = self._page_slice()
        total = len(self._ids)
        if total <= 0:
            self._page_label.setText("No sequences")
            self._prev_btn.setEnabled(False)
            self._next_btn.setEnabled(False)
            return
        start = self._page * self._page_size + 1
        end = min(total, (self._page + 1) * self._page_size)
        self._page_label.setText(
            f"{start}–{end} of {total}  (page {self._page + 1}/{self.page_count()})"
        )
        self._prev_btn.setEnabled(self._page > 0)
        self._next_btn.setEnabled(self._page + 1 < self.page_count())

        cols, cell_side = self._layout_metrics()
        self._layout_cols = cols
        self._layout_side = cell_side
        for c in range(max(cols, 12)):
            self._grid.setColumnStretch(c, 0)
        for pos, (index, sid) in enumerate(entries):
            thumb = DotplotThumb(
                sid,
                index,
                plot_size=cell_side,
                line_width_step=self._line_width_step,
                parent=self._host,
            )
            thumb.clicked.connect(self.thumbClicked.emit)
            cached = self._cache.get(self._cache_key(sid))
            prior = prior_results.get(sid)
            if isinstance(cached, DotplotResult) and self._result_matches_params(
                cached
            ):
                thumb.set_result(cached)
            elif isinstance(cached, str):
                thumb.set_result(None, error=cached)
            elif (
                prior is not None
                and prior[0] is not None
                and self._result_matches_params(prior[0])
            ):
                thumb.set_result(prior[0])
            elif prior is not None and prior[1] and prior[0] is None:
                # Keep hard errors (e.g. sequence too long) across resize only.
                thumb.set_result(None, error=prior[1])
            else:
                thumb.set_busy(True)
                self._pending.append(sid)
            self._thumbs.append(thumb)
            self._grid.addWidget(thumb, pos // cols, pos % cols)
        if self._pending:
            self._timer.start()

    def _layout_metrics(self) -> tuple[int, int]:
        """Return ``(columns, square_side_px)``.

        ``Plot size`` sets how large each square aims to be (both X and Y).
        Columns are chosen from that size; each tile's side then grows equally
        so the row fills the viewport width — X and Y always stay equal.
        """
        avail = max(64, self._scroll.viewport().width() - 2 * _GRID_MARGIN)
        preferred = max(_MIN_THUMB_SIZE, int(self._thumb_size))
        n = len(self._page_slice()) or 1
        outer_pref = preferred + 8
        cols = max(1, (avail + _GRID_GAP) // (outer_pref + _GRID_GAP))
        cols = max(1, min(cols, n, 12))
        outer = (avail - _GRID_GAP * max(0, cols - 1)) // cols
        side = max(_MIN_THUMB_SIZE, outer - 8)
        # Allow filling past the spin max a little so the last column does not
        # leave a large empty gutter; still keep a hard cap for sanity.
        side = min(side, max(_MAX_THUMB_SIZE, preferred) + 80)
        return cols, side

    def _column_count(self) -> int:
        return self._layout_metrics()[0]

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        if not self._thumbs:
            return
        self._resize_timer.start()

    def _on_resize_relayout(self) -> None:
        if not self._thumbs:
            return
        cols, side = self._layout_metrics()
        if cols != self._layout_cols or abs(side - self._layout_side) > 2:
            self._rebuild_page()

    def _cache_key(self, sid: str) -> tuple:
        return (sid, int(self._word), bool(self._reverse), bool(self._ungap))

    def _load_sequence(self, sid: str) -> str:
        if self._loader is None:
            return ""
        try:
            seq = str(self._loader(sid) or "")
        except (OSError, ValueError, TypeError):
            return ""
        if self._ungap:
            seq = seq.replace("-", "").replace(".", "")
        return seq

    def _compute_next(self) -> None:
        if not self._pending:
            self._timer.stop()
            return
        sid = self._pending.pop(0)
        thumb = next((t for t in self._thumbs if t.sequence_id() == sid), None)
        if thumb is None:
            return
        key = self._cache_key(sid)
        cached = self._cache.get(key)
        if isinstance(cached, DotplotResult):
            thumb.set_result(cached)
            return
        if isinstance(cached, str):
            thumb.set_result(None, error=cached)
            return
        seq = self._load_sequence(sid)
        if not seq:
            msg = "Empty"
            self._cache[key] = msg
            thumb.set_result(None, error=msg)
            return
        if len(seq) > FASTGA_PREFERRED_ABOVE_BYTES:
            msg = "Too long\n(use Single)"
            self._cache[key] = msg
            thumb.set_result(None, error=msg)
            return
        try:
            result = compute_self_dotplot(
                seq,
                sequence_id=sid,
                word_size=int(self._word),
                include_reverse=self._reverse,
            )
            self._cache[key] = result
            thumb.set_result(result)
        except DotplotError as exc:
            msg = str(exc).split(".")[0][:40]
            self._cache[key] = msg
            thumb.set_result(None, error=msg)
        except Exception as exc:  # noqa: BLE001
            msg = f"Error: {exc}"[:40]
            self._cache[key] = msg
            thumb.set_result(None, error=msg)
