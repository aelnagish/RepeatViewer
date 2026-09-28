"""Contig view: the result of Tools → Map to Reference (``*.rvmap``).

Like Geneious' contig view: the reference with its annotations on top, then the
consensus, the coverage and the reference bases, then every placed sequence,
with the bases that differ from the reference highlighted. Zoom in far enough
and the bases are written out.

Below, four tables: Sequences (placement, strand, identity…), Features (every
reference feature in every placed sequence: intact, changed, partial, lost),
Variants (every difference, its frequency, cytosine context and the feature it
falls in) and Not used (with the reason). Clicking a row moves the view there.

Pure helpers (lane packing, difference marks, feature cells) have no Qt
drawing and are tested in ``tests/test_contig_view.py``.
"""

from __future__ import annotations

import bisect
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from PySide6.QtCore import QAbstractTableModel, QModelIndex, QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QAction, QColor, QFont, QFontMetrics, QPainter, QPen, QPolygonF
from PySide6.QtWidgets import (
    QAbstractItemView,
    QAbstractScrollArea,
    QApplication,
    QComboBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMenu,
    QSplitter,
    QTableView,
    QTabWidget,
    QToolButton,
    QToolTip,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import map_reference as mr

LABEL_W = 150
RULER_H = 24
FEATURE_ROW_H = 14
TRACK_H = 16
COVERAGE_H = 44
READ_H = 12
READ_H_LETTERS = 16
LETTER_PPB = 7.0
MAX_PPB = 18.0
ONE_PER_ROW_LIMIT = 200

STATUS_COLORS = {
    mr.STATUS_INTACT: "#cdebd3",
    mr.STATUS_CHANGED: "#fbe3b4",
    mr.STATUS_PARTIAL: "#e4e4e4",
    mr.STATUS_LOST: "#f6c4c0",
}
_READ_FILL = QColor("#cfcfcf")
_READ_FILL_SEL = QColor("#9cc0ef")
_INSERT_COLOR = QColor("#8e44ad")
_DELETE_COLOR = QColor("#5a5a5a")


# --------------------------------------------------------------------------
# Pure helpers
# --------------------------------------------------------------------------


def pack_lanes(spans: list[tuple[int, int]], gap: int = 2) -> tuple[list[int], int]:
    """Lane per (start, end) span, greedily, in the given order sorted by start."""
    order = sorted(range(len(spans)), key=lambda i: (spans[i][0], spans[i][1]))
    lane_ends: list[int] = []
    lanes = [0] * len(spans)
    for i in order:
        s, e = spans[i]
        for lane, last in enumerate(lane_ends):
            if s >= last + gap:
                lane_ends[lane] = e
                lanes[i] = lane
                break
        else:
            lanes[i] = len(lane_ends)
            lane_ends.append(e)
    return lanes, len(lane_ends)


def difference_marks(row: str, start: int, reference: str) -> tuple[list[int], list[str],
                                                                   list[tuple[int, int]]]:
    """(mismatch positions, their bases, deletion runs [s, e)) of a placed row."""
    if not row:
        return [], [], []
    a = np.frombuffer(row.encode("ascii", "replace"), dtype=np.uint8)
    gap = a == ord("-")
    positions: list[int] = []
    bases: list[str] = []
    if reference:
        ref = reference[start:start + len(row)]
        b = np.frombuffer(ref.encode("ascii", "replace"), dtype=np.uint8)
        n = min(len(a), len(b))
        diff = np.flatnonzero((a[:n] != b[:n]) & ~gap[:n])
        positions = (diff + start).tolist()
        bases = [row[k] for k in diff.tolist()]
    runs: list[tuple[int, int]] = []
    if gap.any():
        edges = np.diff(np.concatenate(([0], gap.astype(np.int8), [0])))
        run_starts = np.flatnonzero(edges == 1)
        run_ends = np.flatnonzero(edges == -1)
        runs = [(int(s) + start, int(e) + start) for s, e in zip(run_starts, run_ends)]
    return positions, bases, runs


def feature_cell_text(status: str, differences: int) -> str:
    if status == mr.STATUS_INTACT:
        return "intact"
    if status == mr.STATUS_CHANGED:
        return f"{differences} diff"
    return status


def consensus_on_reference(consensus: str, ref_to_cons: np.ndarray, length: int) -> str:
    """Consensus character at every reference position ('-' where it was deleted)."""
    if not consensus or ref_to_cons is None or len(ref_to_cons) == 0:
        return ""
    m = np.asarray(ref_to_cons[:length], dtype=np.int64)
    chars = np.array(list(consensus) + ["-"])
    return "".join(chars[np.where(m >= 0, m, len(consensus))].tolist())


def nice_step(bp_per_px: float, min_px: int = 70) -> int:
    raw = bp_per_px * min_px
    for mag in (1, 10, 100, 1000, 10_000, 100_000, 1_000_000, 10_000_000):
        for m in (1, 2, 5):
            if m * mag >= raw:
                return m * mag
    return 100_000_000


@dataclass
class ReadView:
    """One drawn placement."""

    data: dict
    name: str
    start: int
    end: int
    strand: str
    identity: float
    mism_pos: list[int] = field(default_factory=list)
    mism_base: list[str] = field(default_factory=list)
    del_runs: list[tuple[int, int]] = field(default_factory=list)
    ins_pos: list[int] = field(default_factory=list)
    lane: int = 0

    @property
    def label(self) -> str:
        base = f"{self.name} (reversed)" if self.strand == "-" else self.name
        pl = int(self.data.get("placement", 1) or 1)
        return base if pl == 1 else f"{base} #{pl}"


@dataclass
class ContigModel:
    """What the canvas draws for one reference sequence."""

    ref_id: str
    length: int
    detailed: bool
    reference: str
    consensus: str
    cons_on_ref: str
    coverage: np.ndarray
    bin_size: int
    features: list
    reads: list[ReadView]
    variants: list
    feature_summary: list
    n_reads: int

    @classmethod
    def from_json(cls, entry: dict) -> "ContigModel":
        from repeatviewer.backends.sequence_features import SequenceFeature

        reference = entry.get("reference") or ""
        reads = []
        for r in entry.get("reads", []):
            pos, bases, runs = difference_marks(r["row"], r["start"], reference)
            if not reference and r.get("mism"):
                pos = [int(x[0]) for x in r["mism"]]
                bases = [str(x[1]) for x in r["mism"]]
            reads.append(ReadView(r, r["name"], r["start"], r["end"], r["strand"],
                                  float(r.get("identity", 0)), pos, bases, runs,
                                  sorted(p for p, _s in r.get("ins", []))))
        ref_to_cons = mr.unrle(entry.get("ref_to_cons") or []).astype(np.int64)
        features = []
        for d in entry.get("features", []):
            try:
                features.append(SequenceFeature.from_dict(d))
            except (TypeError, ValueError):
                continue
        return cls(
            ref_id=entry["id"], length=int(entry["length"]), detailed=bool(entry["detailed"]),
            reference=reference, consensus=entry.get("consensus") or "",
            cons_on_ref=consensus_on_reference(entry.get("consensus") or "", ref_to_cons,
                                               int(entry["length"])),
            coverage=mr.unrle(entry.get("coverage") or []), bin_size=int(entry.get("bin_size", 1)),
            features=features, reads=reads, variants=entry.get("variants", []),
            feature_summary=entry.get("feature_summary", []), n_reads=int(entry.get("n_reads", 0)))

    def coverage_at(self, bp: int) -> float:
        i = bp // max(1, self.bin_size)
        return float(self.coverage[i]) if 0 <= i < len(self.coverage) else 0.0


# --------------------------------------------------------------------------
# Canvas
# --------------------------------------------------------------------------


def _mono(size: int) -> QFont:
    font = QFont()
    font.setFamilies(["Menlo", "Monaco", "Consolas", "Courier New", "monospace"])
    font.setStyleHint(QFont.StyleHint.Monospace)
    font.setPointSize(max(6, size))
    return font


def _base_colors() -> dict[str, tuple[QColor, QColor]]:
    """Base → (background, text) from the app's nucleotide colour scheme."""
    try:
        from repeatviewer.widgets.preview_zoom import nucleotide_base_color_map

        return {k: (QColor(v[0]), QColor(v[1])) for k, v in nucleotide_base_color_map().items()}
    except Exception:  # noqa: BLE001 - drawing falls back to fixed colours
        black = QColor("#000000")
        return {"A": (QColor("#6fbf73"), black), "C": (QColor("#5b9bd5"), black),
                "G": (QColor("#f2b134"), black), "T": (QColor("#e0605e"), black),
                "N": (QColor("#bdbdbd"), black)}


_FALLBACK = (QColor("#dddddd"), QColor("#000000"))


class ContigCanvas(QAbstractScrollArea):
    readClicked = Signal(int)

    HIGHLIGHT_DIFF = "diff"
    HIGHLIGHT_ALL = "all"
    HIGHLIGHT_NONE = "none"

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.model: ContigModel | None = None
        self.ppb = 1.0
        self.highlight = self.HIGHLIGHT_DIFF
        self.one_per_row = True
        self.selected = -1
        self._lane_count = 0
        self._feature_rows: list[int] = []
        self._feature_row_count = 0
        self._colors = _base_colors()
        self.viewport().setMouseTracking(True)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.horizontalScrollBar().valueChanged.connect(self.viewport().update)
        self.verticalScrollBar().valueChanged.connect(self.viewport().update)
        self.setMinimumHeight(220)

    # ---------------------------------------------------------------- model

    def set_model(self, model: ContigModel | None) -> None:
        self.model = model
        self.selected = -1
        if model is not None:
            from repeatviewer.widgets.sequence_map_viewer import _pack_feature_rows

            self._feature_rows, self._feature_row_count = _pack_feature_rows(model.features)
            self.one_per_row = len(model.reads) <= ONE_PER_ROW_LIMIT
            self._relayout()
            self.fit()
        self.viewport().update()

    def set_one_per_row(self, on: bool) -> None:
        self.one_per_row = on
        self._relayout()
        self.viewport().update()

    def _relayout(self) -> None:
        m = self.model
        if m is None:
            return
        if self.one_per_row:
            order = sorted(range(len(m.reads)), key=lambda i: (m.reads[i].start, m.reads[i].end))
            for lane, i in enumerate(order):
                m.reads[i].lane = lane
            self._lane_count = len(m.reads)
        else:
            lanes, self._lane_count = pack_lanes([(r.start, r.end) for r in m.reads])
            for r, lane in zip(m.reads, lanes):
                r.lane = lane
        self._update_scrollbars()

    # ---------------------------------------------------------------- geometry

    def letters(self) -> bool:
        return self.ppb >= LETTER_PPB and bool(self.model and self.model.detailed)

    def read_h(self) -> int:
        return READ_H_LETTERS if (self.letters() or self.one_per_row) else READ_H

    def header_h(self) -> int:
        n_tracks = 2 if (self.model and self.model.detailed) else 0   # consensus + reference
        return (RULER_H + max(1, self._feature_row_count) * FEATURE_ROW_H + 6
                + n_tracks * TRACK_H + COVERAGE_H + 8)

    def _track_y(self) -> dict[str, int]:
        y = RULER_H
        out = {"features": y}
        y += max(1, self._feature_row_count) * FEATURE_ROW_H + 6
        if self.model and self.model.detailed:
            out["consensus"] = y
            y += TRACK_H
        out["coverage"] = y
        y += COVERAGE_H + 4
        if self.model and self.model.detailed:
            out["reference"] = y
            y += TRACK_H
        out["reads"] = y + 4
        return out

    def _data_w(self) -> int:
        return max(1, self.viewport().width() - LABEL_W)

    def _update_scrollbars(self) -> None:
        m = self.model
        if m is None:
            return
        content = int(m.length * self.ppb) + 20
        h = self.horizontalScrollBar()
        h.setRange(0, max(0, content - self._data_w()))
        h.setPageStep(self._data_w())
        h.setSingleStep(max(1, int(10 * self.ppb)))
        reads_h = self._lane_count * self.read_h() + 10
        v = self.verticalScrollBar()
        v.setRange(0, max(0, reads_h - (self.viewport().height() - self.header_h())))
        v.setPageStep(max(1, self.viewport().height() - self.header_h()))
        v.setSingleStep(self.read_h())

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._update_scrollbars()

    def x_of(self, bp: float) -> float:
        return LABEL_W + bp * self.ppb - self.horizontalScrollBar().value()

    def bp_at(self, x: float) -> int:
        return int((x - LABEL_W + self.horizontalScrollBar().value()) / self.ppb)

    # ---------------------------------------------------------------- zoom and scroll

    def fit(self) -> None:
        if self.model is None:
            return
        self.ppb = max(1e-6, (self._data_w() - 12) / max(1, self.model.length))
        self._update_scrollbars()
        self.horizontalScrollBar().setValue(0)
        self.viewport().update()

    def zoom(self, factor: float, anchor_x: float | None = None) -> None:
        if self.model is None:
            return
        anchor_x = LABEL_W + self._data_w() / 2 if anchor_x is None else anchor_x
        bp = (anchor_x - LABEL_W + self.horizontalScrollBar().value()) / self.ppb
        fit = (self._data_w() - 12) / max(1, self.model.length)
        self.ppb = min(MAX_PPB, max(fit, self.ppb * factor))
        self._update_scrollbars()
        self.horizontalScrollBar().setValue(int(bp * self.ppb - (anchor_x - LABEL_W)))
        self.viewport().update()

    def show_base(self, bp: int, *, zoom_to_letters: bool = False) -> None:
        if zoom_to_letters and self.ppb < LETTER_PPB:
            self.ppb = 10.0
            self._update_scrollbars()
        self.horizontalScrollBar().setValue(int(bp * self.ppb - self._data_w() / 2))
        self.viewport().update()

    def show_read(self, index: int) -> None:
        m = self.model
        if m is None or not 0 <= index < len(m.reads):
            return
        r = m.reads[index]
        self.selected = index
        span = (r.end - r.start) * self.ppb
        if span > self._data_w():
            self.horizontalScrollBar().setValue(int(r.start * self.ppb) - 10)
        else:
            self.horizontalScrollBar().setValue(int((r.start + r.end) / 2 * self.ppb
                                                    - self._data_w() / 2))
        v = self.verticalScrollBar()
        y = r.lane * self.read_h()
        visible = self.viewport().height() - self.header_h()
        if y < v.value() or y > v.value() + visible - self.read_h():
            v.setValue(max(0, y - visible // 2))
        self.viewport().update()

    def wheelEvent(self, event) -> None:  # noqa: N802
        mods = event.modifiers()
        if mods & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier):
            delta = event.angleDelta().y()
            if delta:
                self.zoom(1.25 if delta > 0 else 0.8, event.position().x())
            event.accept()
            return
        if mods & Qt.KeyboardModifier.ShiftModifier:
            h = self.horizontalScrollBar()
            h.setValue(h.value() - event.angleDelta().y())
            event.accept()
            return
        super().wheelEvent(event)

    # ---------------------------------------------------------------- hit testing

    def read_at(self, pos: QPointF) -> int:
        m = self.model
        if m is None:
            return -1
        top = self._track_y()["reads"]
        y = pos.y() - top + self.verticalScrollBar().value()
        if pos.y() < self.header_h() or y < 0:
            return -1
        lane = int(y // self.read_h())
        bp = self.bp_at(pos.x())
        for i, r in enumerate(m.reads):
            if r.lane == lane and r.start <= bp < r.end:
                return i
            if r.lane == lane and self.one_per_row and pos.x() < LABEL_W:
                return i
        return -1

    def mousePressEvent(self, event) -> None:  # noqa: N802
        i = self.read_at(event.position())
        if i >= 0:
            self.selected = i
            self.readClicked.emit(i)
            self.viewport().update()
        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        if self.read_at(event.position()) < 0:
            self.fit()
        super().mouseDoubleClickEvent(event)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        m = self.model
        if m is None or event.position().x() < LABEL_W:
            QToolTip.hideText()
            return
        bp = self.bp_at(event.position().x())
        if not 0 <= bp < m.length:
            QToolTip.hideText()
            return
        lines = [f"Position {bp + 1:,}"]
        if m.reference and bp < len(m.reference):
            lines[0] += f" · reference {m.reference[bp]}"
        if m.cons_on_ref and bp < len(m.cons_on_ref):
            lines[0] += f" · consensus {m.cons_on_ref[bp]}"
        lines.append(f"Coverage {m.coverage_at(bp):.0f}")
        feats = [f.name for f in m.features if f.start <= bp + 1 <= f.end]
        if feats:
            lines.append("In " + ", ".join(feats[:6]))
        i = self.read_at(event.position())
        if i >= 0:
            r = m.reads[i]
            d = r.data
            lines.append(f"{r.label}: {r.start + 1:,}–{r.end:,}, identity {r.identity:.1%}, "
                         f"{d.get('mismatches', 0)} mismatches, {d.get('deleted', 0)} deleted, "
                         f"{d.get('inserted', 0)} inserted")
            k = bp - r.start
            row = d.get("row", "")
            if 0 <= k < len(row):
                lines.append(f"Base here: {row[k]}")
        QToolTip.showText(event.globalPosition().toPoint(), "\n".join(lines), self)

    # ---------------------------------------------------------------- painting

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self.viewport())
        pal = self.palette()
        bg = pal.color(pal.ColorRole.Base)
        fg = pal.color(pal.ColorRole.Text)
        mid = pal.color(pal.ColorRole.Mid)
        p.fillRect(self.viewport().rect(), bg)
        m = self.model
        if m is None:
            p.setPen(mid)
            p.drawText(self.viewport().rect(), Qt.AlignmentFlag.AlignCenter, "No mapping loaded")
            return
        w, h = self.viewport().width(), self.viewport().height()
        bp0 = max(0, self.bp_at(LABEL_W) - 1)
        bp1 = min(m.length, self.bp_at(w) + 2)
        ys = self._track_y()
        letters = self.letters()
        font = QFont(p.font())
        font.setPointSize(8)
        p.setFont(font)
        header = self.header_h()
        p.setClipRect(QRectF(LABEL_W, 0, w - LABEL_W, header))
        self._paint_ruler(p, bp0, bp1, fg, mid)
        self._paint_features(p, ys["features"], bp0, bp1)
        if m.detailed:
            self._paint_sequence_track(p, ys["consensus"], m.cons_on_ref, bp0, bp1, letters,
                                       compare=m.reference, fg=fg)
            self._paint_sequence_track(p, ys["reference"], m.reference, bp0, bp1, letters,
                                       compare="", fg=fg)
        self._paint_coverage(p, ys["coverage"], bp0, bp1, w)
        p.setClipRect(QRectF(LABEL_W, header, w - LABEL_W, h - header))
        self._paint_reads(p, ys["reads"], bp0, bp1, letters, fg, h)
        # label column
        p.setClipping(False)
        p.fillRect(QRectF(0, 0, LABEL_W, h), pal.color(pal.ColorRole.Window))
        p.setPen(QPen(mid, 1))
        p.drawLine(LABEL_W - 1, 0, LABEL_W - 1, h)
        p.drawLine(0, header - 2, w, header - 2)
        p.setPen(fg)
        label_font = QFont(p.font())
        label_font.setPointSize(9)
        p.setFont(label_font)

        def label(y: float, text: str, height: float = TRACK_H) -> None:
            p.drawText(QRectF(6, y, LABEL_W - 10, height),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, text)

        label(0, m.ref_id, RULER_H)
        label(ys["features"], "Annotations" if m.features else "No annotations", FEATURE_ROW_H)
        if m.detailed:
            label(ys["consensus"], "Consensus")
            label(ys["reference"], "Reference")
        top = max(1.0, float(m.coverage[max(0, bp0 // m.bin_size):
                                         max(1, bp1 // m.bin_size + 1)].max(initial=0)))
        label(ys["coverage"], f"Coverage (max {top:.0f})", COVERAGE_H)
        if self.one_per_row:
            p.setClipRect(QRectF(0, header, LABEL_W - 1, h - header))
            rh = self.read_h()
            off = ys["reads"] - self.verticalScrollBar().value()
            fm = QFontMetrics(label_font)
            for i, r in enumerate(m.reads):
                y = off + r.lane * rh
                if y < header - rh or y > h:
                    continue
                if i == self.selected:
                    p.fillRect(QRectF(0, y, LABEL_W - 1, rh), _READ_FILL_SEL)
                p.setPen(fg)
                p.drawText(QRectF(6, y, LABEL_W - 10, rh),
                           Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                           fm.elidedText(r.label, Qt.TextElideMode.ElideRight, LABEL_W - 12))
        else:
            label(ys["reads"], f"{len(m.reads):,} of {m.n_reads:,} shown"
                  if len(m.reads) < m.n_reads else f"{m.n_reads:,} sequences")

    def _paint_ruler(self, p: QPainter, bp0: int, bp1: int, fg: QColor, mid: QColor) -> None:
        step = nice_step(1.0 / self.ppb)
        p.setPen(QPen(mid, 1))
        p.drawLine(QPointF(self.x_of(0), RULER_H - 4), QPointF(self.x_of(self.model.length),
                                                              RULER_H - 4))
        first = (bp0 // step) * step
        p.setPen(fg)
        for bp in range(int(first), int(bp1) + step, step):
            if bp > self.model.length:
                break
            x = self.x_of(bp)
            p.drawLine(QPointF(x, RULER_H - 8), QPointF(x, RULER_H - 4))
            p.drawText(QPointF(x + 2, RULER_H - 10), f"{max(1, bp):,}")

    def _paint_features(self, p: QPainter, y0: int, bp0: int, bp1: int) -> None:
        from repeatviewer.widgets.sequence_map_viewer import _draw_directed_feature

        m = self.model
        for f, row in zip(m.features, self._feature_rows):
            if f.end < bp0 or f.start - 1 > bp1:
                continue
            x0, x1 = self.x_of(f.start - 1), self.x_of(f.end)
            rect = QRectF(x0, y0 + row * FEATURE_ROW_H + 1, max(2.0, x1 - x0), FEATURE_ROW_H - 3)
            _draw_directed_feature(p, rect, QColor(f.color or "#888888"),
                                   forward=f.is_forward(), label=f.name)

    def _paint_sequence_track(self, p: QPainter, y: int, seq: str, bp0: int, bp1: int,
                              letters: bool, *, compare: str, fg: QColor) -> None:
        if not seq:
            return
        rect_h = TRACK_H - 3
        if not letters:
            p.fillRect(QRectF(self.x_of(0), y + 4, self.model.length * self.ppb, rect_h - 6),
                       QColor("#b9b9b9"))
            if compare:
                p.setPen(Qt.PenStyle.NoPen)
                for bp in range(bp0, min(bp1, len(seq), len(compare))):
                    if seq[bp] != compare[bp]:
                        p.fillRect(QRectF(self.x_of(bp), y + 1, max(1.0, self.ppb), rect_h),
                                   QColor("#d9534f"))
            return
        font = _mono(int(min(11, self.ppb * 0.85)))
        p.setFont(font)
        for bp in range(bp0, min(bp1, len(seq))):
            c = seq[bp]
            rect = QRectF(self.x_of(bp), y, self.ppb, rect_h)
            differs = bool(compare) and bp < len(compare) and c != compare[bp]
            bg_c, fg_c = self._colors.get(c, _FALLBACK)
            if self.highlight == self.HIGHLIGHT_ALL or (differs and self.highlight !=
                                                        self.HIGHLIGHT_NONE) or not compare:
                p.fillRect(rect, bg_c)
                p.setPen(fg_c)
            else:
                p.setPen(fg)
            p.drawText(rect, Qt.AlignmentFlag.AlignCenter, c)

    def _paint_coverage(self, p: QPainter, y: int, bp0: int, bp1: int, w: int) -> None:
        m = self.model
        cov = m.coverage
        if cov.size == 0:
            return
        x_start = int(max(LABEL_W, self.x_of(bp0)))
        x_end = int(min(w, self.x_of(bp1)))
        xs, vals = [], []
        for x in range(x_start, x_end + 1, 2):
            a = max(0, self.bp_at(x))
            b = max(a + 1, self.bp_at(x + 2))
            ia, ib = a // m.bin_size, max(a // m.bin_size + 1, b // m.bin_size)
            seg = cov[ia:ib]
            xs.append(x)
            vals.append(float(seg.max()) if seg.size else 0.0)
        if not vals:
            return
        top = max(1.0, max(vals))
        base = y + COVERAGE_H - 2
        poly = QPolygonF([QPointF(xs[0], base)] + [
            QPointF(x, base - v / top * (COVERAGE_H - 6)) for x, v in zip(xs, vals)]
            + [QPointF(xs[-1], base)])
        p.setPen(QPen(QColor("#3b6fb6"), 1))
        p.setBrush(QColor("#cfe0f5"))
        p.drawPolygon(poly)
        p.setBrush(Qt.BrushStyle.NoBrush)

    def _paint_reads(self, p: QPainter, top: int, bp0: int, bp1: int, letters: bool,
                     fg: QColor, h: int) -> None:
        from repeatviewer.widgets.sequence_map_viewer import _draw_directed_feature

        m = self.model
        rh = self.read_h()
        off = top - self.verticalScrollBar().value()
        header = self.header_h()
        font = _mono(int(min(11, self.ppb * 0.85)))
        for i, r in enumerate(m.reads):
            y = off + r.lane * rh
            if y + rh < header or y > h or r.end < bp0 or r.start > bp1:
                continue
            sel = i == self.selected
            if not letters:
                rect = QRectF(self.x_of(r.start), y + 1, max(2.0, (r.end - r.start) * self.ppb),
                              rh - 3)
                _draw_directed_feature(p, rect, _READ_FILL_SEL if sel else _READ_FILL,
                                       forward=r.strand != "-")
                if self.highlight != self.HIGHLIGHT_NONE:
                    lo = bisect.bisect_left(r.mism_pos, bp0)
                    hi = bisect.bisect_right(r.mism_pos, bp1)
                    for k in range(lo, hi):
                        p.fillRect(QRectF(self.x_of(r.mism_pos[k]), y + 1,
                                          max(1.0, self.ppb), rh - 3),
                                   self._colors.get(r.mism_base[k], _FALLBACK)[0])
                for s, e in r.del_runs:
                    if e < bp0 or s > bp1:
                        continue
                    x0, x1 = self.x_of(s), self.x_of(e)
                    p.fillRect(QRectF(x0, y + 1, max(1.0, x1 - x0), rh - 3), QColor("#ffffff"))
                    p.setPen(QPen(_DELETE_COLOR, 1))
                    p.drawLine(QPointF(x0, y + rh / 2 - 0.5), QPointF(x1, y + rh / 2 - 0.5))
            else:
                p.setFont(font)
                row = r.data["row"]
                ref = m.reference
                for bp in range(max(bp0, r.start), min(bp1, r.end)):
                    c = row[bp - r.start]
                    rect = QRectF(self.x_of(bp), y, self.ppb, rh - 2)
                    if c == "-":
                        p.setPen(QPen(_DELETE_COLOR, 1))
                        p.drawLine(QPointF(rect.left(), rect.center().y()),
                                   QPointF(rect.right(), rect.center().y()))
                        continue
                    differs = bp < len(ref) and c != ref[bp]
                    bg_c, fg_c = self._colors.get(c, _FALLBACK)
                    if self.highlight == self.HIGHLIGHT_ALL or (differs and self.highlight ==
                                                                self.HIGHLIGHT_DIFF):
                        p.fillRect(rect, bg_c)
                        p.setPen(fg_c)
                    else:
                        if sel:
                            p.fillRect(rect, _READ_FILL_SEL)
                        p.setPen(QColor("#000000") if differs else QColor("#6b6b6b"))
                    p.drawText(rect, Qt.AlignmentFlag.AlignCenter, c)
            p.setPen(QPen(_INSERT_COLOR, 2))
            lo = bisect.bisect_left(r.ins_pos, bp0)
            hi = bisect.bisect_right(r.ins_pos, bp1)
            for k in range(lo, hi):
                x = self.x_of(r.ins_pos[k])
                p.drawLine(QPointF(x, y), QPointF(x, y + rh - 2))


# --------------------------------------------------------------------------
# Tables
# --------------------------------------------------------------------------


class _RowsModel(QAbstractTableModel):
    """Read-only table from Python rows, with optional colours and tooltips."""

    def __init__(self, headers: list[str], rows: list[list], *, colors=None, tips=None,
                 header_tips: list[str] | None = None, parent=None) -> None:
        super().__init__(parent)
        self._headers = headers
        self._rows = rows
        self._colors = colors
        self._tips = tips
        self._header_tips = header_tips or []

    def rowCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(self._rows)

    def columnCount(self, parent=QModelIndex()) -> int:  # noqa: N802
        return 0 if parent.isValid() else len(self._headers)

    def data(self, index, role=Qt.ItemDataRole.DisplayRole):
        if not index.isValid():
            return None
        value = self._rows[index.row()][index.column()]
        if role == Qt.ItemDataRole.DisplayRole:
            if isinstance(value, float):
                return f"{value:.1%}" if value <= 1.0 else f"{value:,.1f}"
            return "" if value is None else str(value)
        if role == Qt.ItemDataRole.UserRole:
            return value
        if role == Qt.ItemDataRole.BackgroundRole and self._colors is not None:
            c = self._colors(index.row(), index.column())
            return QColor(c) if c else None
        if role == Qt.ItemDataRole.ForegroundRole and self._colors is not None:
            if self._colors(index.row(), index.column()):
                return QColor("#1d1d1b")
        if role == Qt.ItemDataRole.ToolTipRole and self._tips is not None:
            return self._tips(index.row(), index.column())
        if role == Qt.ItemDataRole.TextAlignmentRole and isinstance(value, (int, float)):
            return int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        return None

    def headerData(self, section, orientation, role=Qt.ItemDataRole.DisplayRole):  # noqa: N802
        if orientation == Qt.Orientation.Horizontal:
            if role == Qt.ItemDataRole.DisplayRole:
                return self._headers[section]
            if role == Qt.ItemDataRole.ToolTipRole and section < len(self._header_tips):
                return self._header_tips[section]
        return None

    def sort(self, column: int, order=Qt.SortOrder.AscendingOrder) -> None:
        self.layoutAboutToBeChanged.emit()
        key = lambda row: (row[column] is None, row[column] if not isinstance(row[column], str)  # noqa: E731
                           else row[column].lower())
        self._rows.sort(key=key, reverse=order == Qt.SortOrder.DescendingOrder)
        self.layoutChanged.emit()


def _table(parent: QWidget) -> QTableView:
    view = QTableView(parent)
    view.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    view.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
    view.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    view.verticalHeader().setVisible(False)
    view.verticalHeader().setDefaultSectionSize(20)
    view.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
    view.setSortingEnabled(True)
    view.setAlternatingRowColors(True)
    return view


# --------------------------------------------------------------------------
# Panel
# --------------------------------------------------------------------------


class ContigViewPanel(QWidget):
    """Toolbar + contig canvas + tables for one ``.rvmap`` result."""

    openRequested = Signal(object)     # Path of a result file to open (report, consensus…)
    fileWritten = Signal(object)       # Path of a file written from the view

    def __init__(self, path: Path | None = None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._path: Path | None = None
        self._data: dict = {}
        self._contig_models: dict[str, ContigModel] = {}
        self.canvas = ContigCanvas(self)
        self._title = QLabel(self)
        self._title.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._contig = QComboBox(self)
        self._contig.setMinimumWidth(180)
        self._contig.currentIndexChanged.connect(self._on_contig_changed)
        self._highlight = QComboBox(self)
        for label, value in (("Highlight differences", ContigCanvas.HIGHLIGHT_DIFF),
                             ("Colour all bases", ContigCanvas.HIGHLIGHT_ALL),
                             ("No highlighting", ContigCanvas.HIGHLIGHT_NONE)):
            self._highlight.addItem(label, value)
        self._highlight.currentIndexChanged.connect(self._on_highlight)
        self._layout_mode = QComboBox(self)
        self._layout_mode.addItem("One sequence per row", True)
        self._layout_mode.addItem("Packed", False)
        self._layout_mode.currentIndexChanged.connect(
            lambda _i: self.canvas.set_one_per_row(bool(self._layout_mode.currentData())))

        def tool(text: str, tip: str, slot) -> QToolButton:
            b = QToolButton(self)
            b.setText(text)
            b.setToolTip(tip)
            b.clicked.connect(slot)
            return b

        zoom_out = tool("−", "Zoom out (Ctrl/Cmd + wheel)", lambda: self.canvas.zoom(0.66))
        zoom_in = tool("+", "Zoom in (Ctrl/Cmd + wheel)", lambda: self.canvas.zoom(1.5))
        fit = tool("Fit", "Show the whole reference (double-click the view)", self.canvas.fit)
        bases = tool("Bases", "Zoom in until the bases are written out",
                     lambda: self.canvas.show_base(self.canvas.bp_at(
                         LABEL_W + self.canvas._data_w() / 2), zoom_to_letters=True))
        self._files_btn = QToolButton(self)
        self._files_btn.setText("Results")
        self._files_btn.setToolTip("Open the report, the consensus, the placed or unused "
                                   "sequences, or the tables")
        self._files_btn.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self._files_menu = QMenu(self._files_btn)
        self._files_btn.setMenu(self._files_menu)

        bar = QWidget(self)
        bl = QHBoxLayout(bar)
        bl.setContentsMargins(6, 4, 6, 4)
        bl.setSpacing(6)
        for wdg in (self._contig, zoom_out, zoom_in, fit, bases, self._highlight,
                    self._layout_mode):
            bl.addWidget(wdg)
        bl.addStretch(1)
        bl.addWidget(self._files_btn)

        self._tabs = QTabWidget(self)
        self._reads_view = _table(self)
        self._features_view = _table(self)
        self._variants_view = _table(self)
        self._unused_view = _table(self)
        self._tabs.addTab(self._reads_view, "Sequences")
        self._tabs.addTab(self._features_view, "Features")
        self._tabs.addTab(self._variants_view, "Variants")
        self._tabs.addTab(self._unused_view, "Not used")
        self._reads_view.clicked.connect(self._on_read_row)
        self._features_view.clicked.connect(self._on_feature_cell)
        self._variants_view.clicked.connect(self._on_variant_row)
        self.canvas.readClicked.connect(self._select_read_row)

        split = QSplitter(Qt.Orientation.Vertical, self)
        split.addWidget(self.canvas)
        split.addWidget(self._tabs)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 2)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(self._title)
        lay.addWidget(bar)
        lay.addWidget(split, 1)
        self._title.setContentsMargins(8, 6, 8, 0)
        if path is not None:
            self.load(path)

    # ---------------------------------------------------------------- loading

    def path(self) -> Path | None:
        return self._path

    def load(self, path: Path) -> None:
        self._path = Path(path)
        try:
            self._data = mr.load_rvmap(self._path)
        except (OSError, ValueError, mr.MapReferenceError) as exc:
            self._data = {}
            self._title.setText(f"Cannot open {self._path.name}: {exc}")
            self.canvas.set_model(None)
            return
        d = self._data
        s = d.get("summary", {})
        engine = mr.ENGINE_NAMES.get(d.get("engine", ""), d.get("engine", ""))
        mean = s.get("mean_identity")
        self._title.setText(
            f"<b>{d.get('title', self._path.name)}</b> &nbsp; {s.get('n_placed', 0):,} of "
            f"{s.get('n_sequences', 0):,} placed ({s.get('n_reversed', 0):,} reversed, "
            f"{s.get('n_unused', 0):,} not used)"
            + (f" · mean identity {mean:.1%}" if mean is not None else "")
            + f" · {engine} {d.get('engine_version', '')}")
        self._title.setToolTip("\n".join(d.get("warnings", [])))
        self._contig_models.clear()
        self._contig.blockSignals(True)
        self._contig.clear()
        for c in d.get("contigs", []):
            self._contig.addItem(f"{c['id']} ({c['n_reads']:,})", c["id"])
        self._contig.blockSignals(False)
        self._contig.setVisible(self._contig.count() > 1)
        self._build_files_menu()
        self._fill_unused()
        if self._contig.count():
            self._on_contig_changed(0)
        else:
            self.canvas.set_model(None)

    def _entry(self, ref_id: str) -> dict | None:
        return next((c for c in self._data.get("contigs", []) if c["id"] == ref_id), None)

    def current_model(self) -> ContigModel | None:
        return self.canvas.model

    def _on_contig_changed(self, _index: int) -> None:
        ref_id = self._contig.currentData()
        entry = self._entry(ref_id) if ref_id else None
        if entry is None:
            return
        model = self._contig_models.get(ref_id)
        if model is None:
            model = ContigModel.from_json(entry)
            self._contig_models[ref_id] = model
        self.canvas.set_model(model)
        self._layout_mode.blockSignals(True)
        self._layout_mode.setCurrentIndex(0 if self.canvas.one_per_row else 1)
        self._layout_mode.blockSignals(False)
        self._fill_reads(model)
        self._fill_features(model)
        self._fill_variants(model)

    def _on_highlight(self, _index: int) -> None:
        self.canvas.highlight = str(self._highlight.currentData())
        self.canvas.viewport().update()

    # ---------------------------------------------------------------- tables

    def _fill_reads(self, m: ContigModel) -> None:
        rows = []
        for i, r in enumerate(m.reads):
            d = r.data
            rows.append([r.label, "−" if r.strand == "-" else "+", r.start + 1, r.end,
                         r.identity, d.get("mismatches", 0), d.get("deleted", 0),
                         d.get("inserted", 0), d.get("mapq") if d.get("mapq", -1) >= 0 else None,
                         d.get("kind", ""), i])
        headers = ["Sequence", "Strand", "Start", "End", "Identity", "Mismatches", "Deleted",
                   "Inserted", "MAPQ", "Kind", "#"]
        self._reads_view.setModel(_RowsModel(headers, rows, parent=self))
        self._reads_view.setColumnHidden(len(headers) - 1, True)
        self._reads_view.setColumnWidth(0, 220)

    def _fill_features(self, m: ContigModel) -> None:
        feats = m.features[:mr.MAX_FEATURE_COLUMNS]
        headers = ["Sequence"] + [f.name for f in feats] + ["#"]
        summary = m.feature_summary
        tips = [""] + [
            (f"{f.name} ({f.feature_type}) {f.start:,}–{f.end:,}\n"
             f"intact {s.get('intact', 0)}, changed {s.get('changed', 0)}, "
             f"partial {s.get('partial', 0)}, lost {s.get('lost', 0)}")
            if (s := (summary[k] if k < len(summary) else {})) else f.name
            for k, f in enumerate(feats)] + [""]
        rows, statuses, details = [], [], []
        for i, r in enumerate(m.reads):
            cells = [""] * len(feats)
            st = [""] * len(feats)
            det = [""] * len(feats)
            for fi, status, ndiff, ident, s_on, e_on in r.data.get("lift", []):
                if fi < len(feats):
                    cells[fi] = feature_cell_text(status, ndiff)
                    st[fi] = status
                    where = f", {s_on:,}–{e_on:,} on the sequence" if s_on else ""
                    det[fi] = f"{feats[fi].name}: {status}, identity {ident:.1%}{where}"
            rows.append([r.label, *cells, i])
            statuses.append(st)
            details.append(det)

        def colors(row: int, col: int):
            if 1 <= col <= len(feats):
                return STATUS_COLORS.get(rows_status(row, col - 1))
            return None

        def rows_status(row: int, k: int) -> str:
            idx = model._rows[row][-1]
            return statuses[idx][k]

        def tip(row: int, col: int):
            if 1 <= col <= len(feats):
                idx = model._rows[row][-1]
                return details[idx][col - 1] or "not covered"
            return None

        model = _RowsModel(headers, rows, colors=colors, tips=tip, header_tips=tips, parent=self)
        self._features_view.setModel(model)
        self._features_view.setColumnHidden(len(headers) - 1, True)
        self._features_view.setColumnWidth(0, 200)
        for c in range(1, len(headers) - 1):
            self._features_view.setColumnWidth(c, 70)
        self._tabs.setTabText(1, f"Features ({len(feats)})" if feats else "Features")

    def _fill_variants(self, m: ContigModel) -> None:
        rows = [[v[0], v[1], v[2], v[3], v[4], (v[3] / v[4]) if v[4] else 0.0, v[5], v[6]]
                for v in m.variants]
        headers = ["Position", "Ref", "Alt", "Count", "Depth", "Frequency", "Context", "In feature"]
        self._variants_view.setModel(_RowsModel(
            headers, rows, parent=self,
            header_tips=["Reference position", "Reference base",
                         "Base in the sequences; - deleted; +XY inserted before the position",
                         "Sequences with this state", "Sequences covering the position",
                         "", "Cytosine context of the reference position (CG / CHG), "
                         "where methylation-driven C→T changes are expected", ""]))
        self._tabs.setTabText(2, f"Variants ({len(rows):,})")

    def _fill_unused(self) -> None:
        rows = [[n, r] for n, r in self._data.get("unused", [])]
        self._unused_view.setModel(_RowsModel(["Sequence", "Reason"], rows, parent=self))
        self._unused_view.setColumnWidth(0, 240)
        self._tabs.setTabText(3, f"Not used ({self._data.get('summary', {}).get('n_unused', 0):,})")

    # ---------------------------------------------------------------- interaction

    def _read_index_from(self, view: QTableView, index) -> int:
        model = view.model()
        return int(model.index(index.row(), model.columnCount() - 1).data(
            Qt.ItemDataRole.UserRole))

    def _on_read_row(self, index) -> None:
        self.canvas.show_read(self._read_index_from(self._reads_view, index))

    def _on_feature_cell(self, index) -> None:
        m = self.canvas.model
        if m is None:
            return
        read_idx = self._read_index_from(self._features_view, index)
        col = index.column()
        if 1 <= col <= min(len(m.features), mr.MAX_FEATURE_COLUMNS):
            f = m.features[col - 1]
            self.canvas.selected = read_idx
            self.canvas.show_read(read_idx)
            self.canvas.show_base((f.start + f.end) // 2)
        else:
            self.canvas.show_read(read_idx)

    def _on_variant_row(self, index) -> None:
        pos = self._variants_view.model().index(index.row(), 0).data(Qt.ItemDataRole.UserRole)
        if pos:
            self.canvas.show_base(int(pos) - 1, zoom_to_letters=True)

    def _select_read_row(self, read_idx: int) -> None:
        model = self._reads_view.model()
        if model is None:
            return
        last = model.columnCount() - 1
        for row in range(model.rowCount()):
            if model.index(row, last).data(Qt.ItemDataRole.UserRole) == read_idx:
                self._reads_view.selectRow(row)
                self._reads_view.scrollTo(model.index(row, 0))
                break

    def _build_files_menu(self) -> None:
        self._files_menu.clear()
        if self._path is None:
            return
        details = self._path.parent / self._data.get("details_dir", "")
        labels = [("report", "Report"), ("consensus", "Consensus"),
                  ("mapped", "Placed sequences"), ("unused", "Unused sequences"),
                  ("variants", "Variants table"), ("features", "Features table"),
                  ("reads", "Sequences table"), ("raw", "Aligner output")]
        for key, text in labels:
            name = self._data.get("details", {}).get(key)
            if not name:
                continue
            target = details / name
            act = QAction(text, self._files_menu)
            act.setEnabled(target.is_file())
            act.triggered.connect(lambda _c=False, t=target: self.openRequested.emit(t))
            self._files_menu.addAction(act)
        self._files_menu.addSeparator()
        ref = Path(self._data.get("reference", ""))
        act = QAction("Reference", self._files_menu)
        act.setEnabled(ref.is_file())
        act.triggered.connect(lambda _c=False, t=ref: self.openRequested.emit(t))
        self._files_menu.addAction(act)
        copy = QAction("Copy consensus", self._files_menu)
        copy.triggered.connect(self._copy_consensus)
        self._files_menu.addAction(copy)

    def _copy_consensus(self) -> None:
        m = self.canvas.model
        if m is not None and m.consensus:
            QApplication.clipboard().setText(f">{m.ref_id}_consensus\n{m.consensus}\n")
