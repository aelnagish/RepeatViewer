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
COVERAGE_H = 30
IDENTITY_H = 6
READ_H = 12
READ_H_LETTERS = 16
LETTER_PPB = 6.0
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
_REVERSE_TINT = QColor(123, 63, 181, 22)    # reverse-strand rows, light lavender
_MARK_FILL = QColor(31, 95, 191, 28)        # current difference column
_MARK_EDGE = QColor("#1f5fbf")
_NO_COVERAGE = QColor("#d9d9d6")
_IDENTITY_OK = QColor("#5fae66")
_IDENTITY_SOME = QColor("#e0b43a")
_IDENTITY_MANY = QColor("#d0574a")
_CONS_DIFF_MARK = QColor("#e0a020")

# IUPAC codes: which plain bases each letter stands for. A read base "matches"
# the reference when it is one of the bases the reference letter allows, so an
# ambiguity code in the reference (Y, R, N…) is not drawn as a difference.
IUPAC = {
    "A": "A", "C": "C", "G": "G", "T": "T", "U": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}


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


def _compat_table() -> np.ndarray:
    """256 x 256 table: True when read byte b agrees with reference byte a."""
    table = np.zeros((256, 256), dtype=bool)
    for i in range(256):
        table[i, i] = True
    for ref_code, allowed in IUPAC.items():
        for ref_c in (ref_code, ref_code.lower()):
            for base in allowed:
                for read_c in (base, base.lower()):
                    table[ord(ref_c), ord(read_c)] = True
            for read_c in (ref_code, ref_code.lower()):
                table[ord(ref_c), ord(read_c)] = True
    # case never makes a difference
    for c in range(ord("A"), ord("Z") + 1):
        table[c, c + 32] = table[c + 32, c] = True
    return table


_COMPAT = _compat_table()


def bases_match(ref_base: str, base: str) -> bool:
    """True when ``base`` agrees with ``ref_base`` (IUPAC aware, case blind)."""
    if not ref_base or not base:
        return False
    a, b = ord(ref_base[0]), ord(base[0])
    return bool(a < 256 and b < 256 and _COMPAT[a, b])


def difference_profile(reads: list, length: int) -> np.ndarray:
    """Number of placed sequences that differ from the reference at each position
    (a mismatch or a deleted base)."""
    out = np.zeros(max(0, int(length)), dtype=np.int32)
    for r in reads:
        if r.mism_pos:
            pos = np.asarray(r.mism_pos, dtype=np.int64)
            pos = pos[(pos >= 0) & (pos < length)]
            np.add.at(out, pos, 1)
        for s, e in r.del_runs:
            s, e = max(0, s), min(length, e)
            if e > s:
                out[s:e] += 1
    return out


def difference_positions(variants: list) -> list[int]:
    """Sorted 0-based reference positions that have a variant (Variants table)."""
    out = set()
    for v in variants:
        try:
            out.add(int(v[0]) - 1)
        except (TypeError, ValueError, IndexError):
            continue
    return sorted(p for p in out if p >= 0)


def step_difference(positions: list[int], current: int | None, direction: int) -> int | None:
    """Next (direction > 0) or previous difference from ``current``; wraps around."""
    if not positions:
        return None
    if current is None:
        return positions[0] if direction > 0 else positions[-1]
    if direction > 0:
        k = bisect.bisect_right(positions, current)
        return positions[k] if k < len(positions) else positions[0]
    k = bisect.bisect_left(positions, current) - 1
    return positions[k] if k >= 0 else positions[-1]


SORT_START = "start"
SORT_STRAND = "strand"
SORT_IDENTITY = "identity"
SORT_NAME = "name"
COLOR_STRAND = "strand"
COLOR_IDENTITY = "identity"
COLOR_MAPQ = "mapq"


def read_passes(data: dict, identity: float, *, min_identity: float = 0.0,
                min_mapq: int = 0, hide_secondary: bool = False) -> bool:
    """Whether a placement is shown under the current filters. MAPQ < 0 (not
    reported by the mapper) never fails the MAPQ filter."""
    if identity < min_identity:
        return False
    mapq = data.get("mapq", -1)
    if min_mapq > 0 and mapq is not None and mapq >= 0 and mapq < min_mapq:
        return False
    if hide_secondary and str(data.get("kind", "primary")) != "primary":
        return False
    return True


def sort_order(reads: list, key: str) -> list[int]:
    """Indices of ``reads`` in display order for one-per-row layout."""
    idx = list(range(len(reads)))
    if key == SORT_STRAND:
        return sorted(idx, key=lambda i: (reads[i].strand == "-", reads[i].start, reads[i].end))
    if key == SORT_IDENTITY:
        return sorted(idx, key=lambda i: (-reads[i].identity, reads[i].start))
    if key == SORT_NAME:
        return sorted(idx, key=lambda i: (reads[i].name.lower(), reads[i].start))
    return sorted(idx, key=lambda i: (reads[i].start, reads[i].end))


def _mix(a: tuple[int, int, int], b: tuple[int, int, int], t: float) -> str:
    t = min(1.0, max(0.0, t))
    return "#{:02x}{:02x}{:02x}".format(*(round(x + (y - x) * t) for x, y in zip(a, b)))


def read_color_hex(strand: str, identity: float, mapq: int | None, mode: str) -> str:
    """Bar colour of a placement when zoomed out."""
    if mode == COLOR_IDENTITY:
        # 80 % or lower = red, 100 % = green
        return _mix((208, 87, 74), (95, 174, 102), (identity - 0.8) / 0.2)
    if mode == COLOR_MAPQ:
        if mapq is None or mapq < 0:
            return "#cfcfcf"
        return _mix((224, 180, 58), (63, 120, 201), mapq / 60.0)
    return "#a98bd6" if strand == "-" else "#8fb3e6"


def feature_status_text(summary: dict) -> str:
    """Short status for a reference feature over all placements:
    "intact", "changed", "3/7 intact"…"""
    counts = {k: int(summary.get(k, 0) or 0) for k in
              (mr.STATUS_INTACT, mr.STATUS_CHANGED, mr.STATUS_PARTIAL, mr.STATUS_LOST)}
    total = sum(counts.values())
    if total == 0:
        return ""
    for status, n in counts.items():
        if n == total:
            return status
    return f"{counts[mr.STATUS_INTACT]}/{total} intact"


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
        # IUPAC aware: a read C or T under a reference Y is not a difference.
        diff = np.flatnonzero(~_COMPAT[b[:n], a[:n]] & ~gap[:n])
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
    hidden: bool = False    # filtered out (identity / MAPQ / secondary)

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
    diff_count: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int32))
    diff_pos: list[int] = field(default_factory=list)

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
        length = int(entry["length"])
        variants = entry.get("variants", [])
        return cls(
            ref_id=entry["id"], length=length, detailed=bool(entry["detailed"]),
            reference=reference, consensus=entry.get("consensus") or "",
            cons_on_ref=consensus_on_reference(entry.get("consensus") or "", ref_to_cons,
                                               int(entry["length"])),
            coverage=mr.unrle(entry.get("coverage") or []), bin_size=int(entry.get("bin_size", 1)),
            features=features, reads=reads, variants=variants,
            feature_summary=entry.get("feature_summary", []), n_reads=int(entry.get("n_reads", 0)),
            diff_count=difference_profile(reads, length),
            diff_pos=difference_positions(variants))

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
    viewChanged = Signal()       # zoom or scroll (the overview strip follows)

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
        self.mark_bp: int | None = None     # current difference column (0-based)
        self.sort_key = SORT_START
        self.color_by = COLOR_STRAND
        self.min_identity = 0.0
        self.min_mapq = 0
        self.hide_secondary = False
        self._lane_count = 0
        self._feature_rows: list[int] = []
        self._feature_row_count = 0
        self._colors = _base_colors()
        self.viewport().setMouseTracking(True)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.horizontalScrollBar().valueChanged.connect(self.viewport().update)
        self.horizontalScrollBar().valueChanged.connect(lambda _v: self.viewChanged.emit())
        self.horizontalScrollBar().rangeChanged.connect(lambda _a, _b: self.viewChanged.emit())
        self.verticalScrollBar().valueChanged.connect(self.viewport().update)
        self.setMinimumHeight(220)

    # ---------------------------------------------------------------- model

    def set_model(self, model: ContigModel | None) -> None:
        self.model = model
        self.selected = -1
        self.mark_bp = None
        if model is not None:
            from repeatviewer.widgets.sequence_map_viewer import _pack_feature_rows

            self._feature_rows, self._feature_row_count = _pack_feature_rows(model.features)
            self.one_per_row = len(model.reads) <= ONE_PER_ROW_LIMIT
            self._apply_filters()
            self._relayout()
            self.fit()
        self.viewport().update()

    def set_one_per_row(self, on: bool) -> None:
        self.one_per_row = on
        self._relayout()
        self.viewport().update()

    def set_sort(self, key: str) -> None:
        self.sort_key = key
        self._relayout()
        self.viewport().update()

    def set_color_by(self, mode: str) -> None:
        self.color_by = mode
        self.viewport().update()

    def set_filters(self, *, min_identity: float, min_mapq: int, hide_secondary: bool) -> None:
        self.min_identity = float(min_identity)
        self.min_mapq = int(min_mapq)
        self.hide_secondary = bool(hide_secondary)
        self._apply_filters()
        self._relayout()
        self.viewport().update()

    def hidden_count(self) -> int:
        return sum(1 for r in self.model.reads if r.hidden) if self.model else 0

    def _apply_filters(self) -> None:
        m = self.model
        if m is None:
            return
        for i, r in enumerate(m.reads):
            r.hidden = not read_passes(r.data, r.identity, min_identity=self.min_identity,
                                       min_mapq=self.min_mapq,
                                       hide_secondary=self.hide_secondary)
            if r.hidden and i == self.selected:
                self.selected = -1

    def _relayout(self) -> None:
        m = self.model
        if m is None:
            return
        shown = [i for i, r in enumerate(m.reads) if not r.hidden]
        for r in m.reads:
            r.lane = -1
        if self.one_per_row:
            order = sort_order([m.reads[i] for i in shown], self.sort_key)
            for lane, k in enumerate(order):
                m.reads[shown[k]].lane = lane
            self._lane_count = len(shown)
        else:
            lanes, self._lane_count = pack_lanes([(m.reads[i].start, m.reads[i].end)
                                                  for i in shown])
            for i, lane in zip(shown, lanes):
                m.reads[i].lane = lane
        self._update_scrollbars()

    def _read_fill(self, r: ReadView) -> QColor:
        mapq = r.data.get("mapq", -1)
        return QColor(read_color_hex(r.strand, r.identity,
                                     None if mapq is None else int(mapq), self.color_by))

    # ---------------------------------------------------------------- geometry

    def letters(self) -> bool:
        return self.ppb >= LETTER_PPB and bool(self.model and self.model.detailed)

    def read_h(self) -> int:
        return READ_H_LETTERS if (self.letters() or self.one_per_row) else READ_H

    def header_h(self) -> int:
        detailed = bool(self.model and self.model.detailed)
        n_tracks = 2 if detailed else 0   # consensus + reference
        return (RULER_H + max(1, self._feature_row_count) * FEATURE_ROW_H + 6
                + n_tracks * TRACK_H + (IDENTITY_H + 3 if detailed else 0) + COVERAGE_H + 8)

    def _track_y(self) -> dict[str, int]:
        y = RULER_H
        out = {"features": y}
        y += max(1, self._feature_row_count) * FEATURE_ROW_H + 6
        if self.model and self.model.detailed:
            out["consensus"] = y
            y += TRACK_H
            out["identity"] = y + 1
            y += IDENTITY_H + 3
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
        self.viewChanged.emit()

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
        self.viewChanged.emit()

    def show_base(self, bp: int, *, zoom_to_letters: bool = False, mark: bool = False) -> None:
        if mark:
            self.mark_bp = int(bp)
        if zoom_to_letters and self.ppb < LETTER_PPB:
            self.ppb = 12.0
            self._update_scrollbars()
        self.horizontalScrollBar().setValue(int(bp * self.ppb - self._data_w() / 2))
        self.viewport().update()
        self.viewChanged.emit()

    def show_read(self, index: int) -> None:
        m = self.model
        if m is None or not 0 <= index < len(m.reads):
            return
        r = m.reads[index]
        if r.hidden:
            return
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
            if r.hidden:
                continue
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
        cov = m.coverage_at(bp)
        if cov <= 0:
            lines.append("No coverage")
        else:
            line = f"Depth {cov:.0f}"
            if m.detailed and bp < len(m.diff_count) and m.diff_count[bp]:
                k = int(m.diff_count[bp])
                line += f" · {k} differ ({k / cov:.0%})"
            counts = self._base_counts(bp)
            if counts:
                line += " · " + " ".join(f"{b} {n}" for b, n in counts)
            lines.append(line)
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

    def _base_counts(self, bp: int, limit: int = 5000) -> list[tuple[str, int]]:
        """Bases of the placed sequences at ``bp`` (detailed results only)."""
        m = self.model
        if m is None or not m.detailed:
            return []
        counts: dict[str, int] = {}
        for r in m.reads[:limit]:
            if r.start <= bp < r.end:
                row = r.data.get("row", "")
                k = bp - r.start
                if 0 <= k < len(row):
                    c = row[k].upper()
                    counts[c] = counts.get(c, 0) + 1
        order = "ACGT-"
        return sorted(counts.items(), key=lambda kv: (order.find(kv[0]) % 9, kv[0]))

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
            self._paint_identity(p, ys["identity"], bp0, bp1, w)
            self._paint_sequence_track(p, ys["reference"], m.reference, bp0, bp1, letters,
                                       compare="", fg=fg)
        self._paint_coverage(p, ys["coverage"], bp0, bp1, w)
        p.setClipRect(QRectF(LABEL_W, header, w - LABEL_W, h - header))
        self._paint_reads(p, ys["reads"], bp0, bp1, letters, fg, h)
        if self.mark_bp is not None and bp0 <= self.mark_bp <= bp1:
            p.setClipRect(QRectF(LABEL_W, RULER_H, w - LABEL_W, h - RULER_H))
            x = self.x_of(self.mark_bp)
            col = QRectF(x, RULER_H, max(2.0, self.ppb), h - RULER_H)
            p.fillRect(col, _MARK_FILL)
            p.setPen(QPen(_MARK_EDGE, 1))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(col)
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

        label(0, f"{m.length:,} bp", RULER_H)
        label(ys["features"], "Annotations" if m.features else "No annotations", FEATURE_ROW_H)
        if m.detailed:
            label(ys["consensus"], "Consensus")
            small = QFont(label_font)
            small.setPointSize(7)
            p.setFont(small)
            p.drawText(QRectF(6, ys["identity"] - 3, LABEL_W - 10, IDENTITY_H + 6),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, "Identity")
            p.setFont(label_font)
            fm_ref = QFontMetrics(label_font)
            label(ys["reference"], fm_ref.elidedText(f"Reference · {m.ref_id}",
                                                     Qt.TextElideMode.ElideRight, LABEL_W - 12))
        top = max(1.0, float(m.coverage[max(0, bp0 // m.bin_size):
                                         max(1, bp1 // m.bin_size + 1)].max(initial=0)))
        label(ys["coverage"], f"Coverage (max {top:.0f})", COVERAGE_H)
        if self.one_per_row:
            p.setClipRect(QRectF(0, header, LABEL_W - 1, h - header))
            rh = self.read_h()
            off = ys["reads"] - self.verticalScrollBar().value()
            fm = QFontMetrics(label_font)
            for i, r in enumerate(m.reads):
                if r.hidden:
                    continue
                y = off + r.lane * rh
                if y < header - rh or y > h:
                    continue
                if i == self.selected:
                    p.fillRect(QRectF(0, y, LABEL_W - 1, rh), _READ_FILL_SEL)
                elif r.strand == "-":
                    p.fillRect(QRectF(0, y, LABEL_W - 1, rh), _REVERSE_TINT)
                p.setPen(fg)
                p.drawText(QRectF(6, y, LABEL_W - 10, rh),
                           Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                           fm.elidedText(r.label, Qt.TextElideMode.ElideRight, LABEL_W - 12))
        else:
            shown = len(m.reads) - self.hidden_count()
            label(ys["reads"], f"{shown:,} of {m.n_reads:,} shown"
                  if shown < m.n_reads else f"{m.n_reads:,} sequences")

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
        summary = m.feature_summary
        for k, (f, row) in enumerate(zip(m.features, self._feature_rows)):
            if f.end < bp0 or f.start - 1 > bp1:
                continue
            x0, x1 = self.x_of(f.start - 1), self.x_of(f.end)
            rect = QRectF(x0, y0 + row * FEATURE_ROW_H + 1, max(2.0, x1 - x0), FEATURE_ROW_H - 3)
            status = feature_status_text(summary[k]) if k < len(summary) else ""
            _draw_directed_feature(p, rect, QColor(f.color or "#888888"),
                                   forward=f.is_forward(),
                                   label=f"{f.name} · {status}" if status else f.name)

    def _hatch(self, p: QPainter, rect: QRectF) -> None:
        """No coverage: light diagonal hatching, never a difference colour."""
        from PySide6.QtGui import QBrush

        p.fillRect(rect, QBrush(_NO_COVERAGE, Qt.BrushStyle.BDiagPattern))

    def _uncovered_runs(self, bp0: int, bp1: int) -> list[tuple[int, int]]:
        """[start, end) runs without coverage between bp0 and bp1 (vectorised)."""
        m = self.model
        step = max(1, m.bin_size)
        i0, i1 = bp0 // step, min(len(m.coverage), bp1 // step + 1)
        if i1 <= i0:
            return []
        empty = (np.asarray(m.coverage[i0:i1]) <= 0).astype(np.int8)
        if not empty.any():
            return []
        edges = np.diff(np.concatenate(([0], empty, [0])))
        starts = np.flatnonzero(edges == 1)
        ends = np.flatnonzero(edges == -1)
        return [(max(bp0, (i0 + int(a)) * step), min(bp1, (i0 + int(b)) * step))
                for a, b in zip(starts, ends)]

    def _paint_sequence_track(self, p: QPainter, y: int, seq: str, bp0: int, bp1: int,
                              letters: bool, *, compare: str, fg: QColor) -> None:
        """Reference (compare == "") or consensus (compare = reference)."""
        if not seq:
            return
        m = self.model
        rect_h = TRACK_H - 3
        is_consensus = bool(compare)
        uncovered = self._uncovered_runs(bp0, bp1) if is_consensus else []

        def covered(bp: int) -> bool:
            return not is_consensus or m.coverage_at(bp) > 0

        if not letters:
            p.fillRect(QRectF(self.x_of(0), y + 4, m.length * self.ppb, rect_h - 6),
                       QColor("#b9b9b9"))
            for s0, e0 in uncovered:
                r = QRectF(self.x_of(s0), y + 1, max(1.0, (e0 - s0) * self.ppb), rect_h)
                p.fillRect(r, pal_base(self))
                self._hatch(p, r)
            if is_consensus:
                hi = min(bp1, len(seq), len(compare))
                if hi > bp0:
                    a = np.frombuffer(seq[bp0:hi].encode("ascii", "replace"), dtype=np.uint8)
                    b = np.frombuffer(compare[bp0:hi].encode("ascii", "replace"), dtype=np.uint8)
                    step = max(1, m.bin_size)
                    cov = np.repeat(np.asarray(m.coverage), step)[bp0:hi] if step > 1 else \
                        np.asarray(m.coverage)[bp0:hi]
                    n = min(len(a), len(cov))
                    bad = np.flatnonzero(~_COMPAT[b[:n], a[:n]] & (cov[:n] > 0)) + bp0
                    red = QColor("#d9534f")
                    for bp in bad.tolist():
                        p.fillRect(QRectF(self.x_of(bp), y + 1, max(1.0, self.ppb), rect_h), red)
            return
        font = _mono(int(min(11, max(7, self.ppb * 0.85))))
        p.setFont(font)
        for bp in range(bp0, min(bp1, len(seq))):
            rect = QRectF(self.x_of(bp), y, self.ppb, rect_h)
            if not covered(bp):
                self._hatch(p, rect)
                continue
            c = seq[bp]
            differs = is_consensus and bp < len(compare) and not bases_match(compare[bp], c)
            bg_c, fg_c = self._colors.get(c.upper(), _FALLBACK)
            ambiguous = c.upper() not in "ACGT-"
            if ambiguous and not is_consensus:
                # ambiguity code in the reference: pale fill, dashed outline
                p.fillRect(rect, QColor("#fff4d6"))
                p.setPen(QPen(QColor("#a35f00"), 1, Qt.PenStyle.DashLine))
                p.drawRect(rect.adjusted(0.5, 0.5, -0.5, -0.5))
                p.setPen(QColor("#6b4100"))
            elif (self.highlight == self.HIGHLIGHT_ALL or not is_consensus
                  or (differs and self.highlight != self.HIGHLIGHT_NONE)):
                p.fillRect(rect, bg_c)
                p.setPen(fg_c)
            else:
                p.setPen(fg)
            p.drawText(rect, Qt.AlignmentFlag.AlignCenter, c)
            if (is_consensus and bp < len(m.diff_count) and m.diff_count[bp]
                    and self.highlight != self.HIGHLIGHT_NONE):
                # some sequences differ here: amber underline
                p.fillRect(QRectF(rect.left(), rect.bottom() - 2, rect.width(), 2),
                           _CONS_DIFF_MARK)

    def _paint_identity(self, p: QPainter, y: int, bp0: int, bp1: int, w: int) -> None:
        """Green where all sequences agree with the reference, yellow where a few
        differ, red where many (>= 30 %) differ; nothing where there is no coverage."""
        m = self.model
        if m.diff_count.size == 0:
            return
        x_start = int(max(LABEL_W, self.x_of(bp0)))
        x_end = int(min(w, self.x_of(bp1)))
        step = max(1, int(round(self.ppb))) if self.ppb >= 1 else 1
        for x in range(x_start, x_end + 1, step):
            a = max(0, self.bp_at(x))
            b = min(m.length, max(a + 1, self.bp_at(x + step)))
            if a >= m.length:
                break
            cov = max(m.coverage_at(k) for k in range(a, b)) if b - a <= 64 else \
                m.coverage_at((a + b) // 2)
            if cov <= 0:
                continue
            share = float(m.diff_count[a:b].max(initial=0)) / cov
            color = _IDENTITY_OK if share == 0 else (_IDENTITY_SOME if share < 0.3
                                                     else _IDENTITY_MANY)
            p.fillRect(QRectF(x, y, step, IDENTITY_H), color)

    def _paint_coverage(self, p: QPainter, y: int, bp0: int, bp1: int, w: int) -> None:
        """Coverage on a log scale, so one or two sequences do not fill the track."""
        import math

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
        scale = math.log1p(max(top, 4.0))      # one read stays a low band, not a wall
        span = COVERAGE_H - 6
        base = y + COVERAGE_H - 2
        poly = QPolygonF([QPointF(xs[0], base)] + [
            QPointF(x, base - math.log1p(v) / scale * span) for x, v in zip(xs, vals)]
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
        w = self.viewport().width()
        font = _mono(int(min(11, max(7, self.ppb * 0.85))))
        dot_color = QColor("#9a9a96")
        for i, r in enumerate(m.reads):
            if r.hidden:
                continue
            y = off + r.lane * rh
            if y + rh < header or y > h:
                continue
            sel = i == self.selected
            if self.one_per_row and r.strand == "-" and not sel:
                p.fillRect(QRectF(LABEL_W, y, w - LABEL_W, rh), _REVERSE_TINT)
            if r.end < bp0 or r.start > bp1:
                continue
            if not letters:
                rect = QRectF(self.x_of(r.start), y + 1, max(2.0, (r.end - r.start) * self.ppb),
                              rh - 3)
                _draw_directed_feature(p, rect, _READ_FILL_SEL if sel else self._read_fill(r),
                                       forward=r.strand != "-")
                if sel:
                    p.setPen(QPen(QColor("#1d1d1f"), 1.5))
                    p.setBrush(Qt.BrushStyle.NoBrush)
                    p.drawRect(rect.adjusted(-1, -1, 1, 1))
                if self.highlight != self.HIGHLIGHT_NONE:
                    lo = bisect.bisect_left(r.mism_pos, bp0)
                    hi = bisect.bisect_right(r.mism_pos, bp1)
                    for k in range(lo, hi):
                        p.fillRect(QRectF(self.x_of(r.mism_pos[k]), y + 1,
                                          max(1.0, self.ppb), rh - 3),
                                   self._colors.get(r.mism_base[k].upper(), _FALLBACK)[0])
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
                if sel:
                    p.fillRect(QRectF(self.x_of(max(bp0, r.start)), y,
                                      (min(bp1, r.end) - max(bp0, r.start)) * self.ppb, rh - 2),
                               _READ_FILL_SEL)
                for bp in range(max(bp0, r.start), min(bp1, r.end)):
                    c = row[bp - r.start]
                    rect = QRectF(self.x_of(bp), y, self.ppb, rh - 2)
                    if c == "-":
                        # deleted base: dark cell with a dash
                        p.fillRect(rect.adjusted(0, 1, 0, -1), _DELETE_COLOR)
                        p.setPen(QColor("#ffffff"))
                        p.drawText(rect, Qt.AlignmentFlag.AlignCenter, "-")
                        continue
                    differs = bp < len(ref) and not bases_match(ref[bp], c)
                    bg_c, fg_c = self._colors.get(c.upper(), _FALLBACK)
                    if self.highlight == self.HIGHLIGHT_ALL:
                        p.fillRect(rect, bg_c)
                        p.setPen(fg_c)
                        p.drawText(rect, Qt.AlignmentFlag.AlignCenter, c)
                    elif self.highlight == self.HIGHLIGHT_DIFF:
                        if differs:
                            p.fillRect(rect, bg_c)
                            p.setPen(fg_c)
                            bold = QFont(font)
                            bold.setBold(True)
                            p.setFont(bold)
                            p.drawText(rect, Qt.AlignmentFlag.AlignCenter, c)
                            p.setFont(font)
                        else:
                            # agrees with the reference: a dot, so differences stand out
                            p.setPen(dot_color)
                            p.drawText(rect, Qt.AlignmentFlag.AlignCenter, "·")
                    else:
                        p.setPen(fg)
                        p.drawText(rect, Qt.AlignmentFlag.AlignCenter, c)
            p.setPen(QPen(_INSERT_COLOR, 2))
            lo = bisect.bisect_left(r.ins_pos, bp0)
            hi = bisect.bisect_right(r.ins_pos, bp1)
            for k in range(lo, hi):
                x = self.x_of(r.ins_pos[k])
                p.drawLine(QPointF(x, y), QPointF(x, y + rh - 2))


class ContigOverview(QWidget):
    """Whole reference in one strip: coverage (amber where sequences differ) and a
    box for the part shown below. Click or drag to move the view."""

    HEIGHT = 30

    def __init__(self, canvas: ContigCanvas, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.canvas = canvas
        self.setFixedHeight(self.HEIGHT)
        self.setMouseTracking(True)
        self.setToolTip("Whole reference. Click or drag to move the view.")
        canvas.viewChanged.connect(self.update)

    def _x0(self) -> float:
        return float(LABEL_W)

    def _scale(self) -> float:
        m = self.canvas.model
        return max(1.0, self.width() - LABEL_W - 8) / max(1, m.length if m else 1)

    def paintEvent(self, event) -> None:  # noqa: N802
        import math

        p = QPainter(self)
        pal = self.palette()
        p.fillRect(self.rect(), pal.color(pal.ColorRole.Window))
        m = self.canvas.model
        if m is None:
            return
        p.setPen(pal.color(pal.ColorRole.Text))
        f = QFont(p.font())
        f.setPointSize(8)
        p.setFont(f)
        p.drawText(QRectF(6, 0, LABEL_W - 10, self.HEIGHT),
                   Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                   f"Whole reference · {m.length:,} bp")
        x0, sc = self._x0(), self._scale()
        track = QRectF(x0, 4, m.length * sc, self.HEIGHT - 8)
        p.fillRect(track, QColor("#ececea"))
        cov = np.asarray(m.coverage, dtype=float)
        if cov.size:
            top = math.log1p(max(4.0, float(cov.max())))
            n_px = max(1, int(track.width()))
            edges = np.linspace(0, m.length, n_px + 1).astype(np.int64)
            diff = m.diff_count if m.diff_count.size else None
            ok, amber = QColor("#9dbbe6"), QColor("#e0a020")
            for k in range(n_px):
                a, b = int(edges[k]), max(int(edges[k]) + 1, int(edges[k + 1]))
                ia, ib = a // max(1, m.bin_size), max(a // max(1, m.bin_size) + 1,
                                                       b // max(1, m.bin_size))
                v = float(cov[ia:ib].max(initial=0))
                if v <= 0:
                    continue
                hgt = max(2.0, math.log1p(v) / top * track.height())
                differs = diff is not None and bool(diff[a:b].max(initial=0))
                p.fillRect(QRectF(x0 + k, track.bottom() - hgt, 1.0, hgt),
                           amber if differs else ok)
        c = self.canvas
        if c.ppb > 0:
            vis0 = max(0, c.bp_at(LABEL_W))
            vis1 = min(m.length, c.bp_at(c.viewport().width()))
            box = QRectF(x0 + vis0 * sc, 2, max(4.0, (vis1 - vis0) * sc), self.HEIGHT - 4)
            p.setPen(QPen(_MARK_EDGE, 2))
            p.setBrush(QColor(31, 95, 191, 26))
            p.drawRoundedRect(box, 3, 3)
        p.end()

    def _jump(self, x: float) -> None:
        m = self.canvas.model
        if m is None:
            return
        bp = int((x - self._x0()) / self._scale())
        self.canvas.show_base(min(max(0, bp), m.length - 1))

    def mousePressEvent(self, event) -> None:  # noqa: N802
        self._jump(event.position().x())

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if event.buttons() & Qt.MouseButton.LeftButton:
            self._jump(event.position().x())


def pal_base(widget: QWidget) -> QColor:
    pal = widget.palette()
    return pal.color(pal.ColorRole.Base)


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
        self.overview = ContigOverview(self.canvas, self)
        self._title = QLabel(self)
        self._title.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._contig = QComboBox(self)
        self._contig.setMinimumWidth(180)
        self._contig.currentIndexChanged.connect(self._on_contig_changed)
        self._highlight = QComboBox(self)
        for label, value in (("Dots where they match", ContigCanvas.HIGHLIGHT_DIFF),
                             ("Colour all bases", ContigCanvas.HIGHLIGHT_ALL),
                             ("Plain letters", ContigCanvas.HIGHLIGHT_NONE)):
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
        self._prev_diff = tool("◀ Previous difference", "Go to the previous variant position",
                               lambda: self._step_difference(-1))
        self._next_diff = tool("Next difference ▶", "Go to the next variant position",
                               lambda: self._step_difference(1))
        self._diff_label = QLabel(self)
        self._diff_label.setMinimumWidth(60)
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
                    self._layout_mode, self._prev_diff, self._diff_label, self._next_diff):
            bl.addWidget(wdg)
        bl.addStretch(1)
        bl.addWidget(self._files_btn)

        # second row: how reads are ordered, coloured and filtered
        from repeatviewer.widgets.combo_boxes import NoWheelComboBox
        from repeatviewer.widgets.spin_boxes import NoWheelDoubleSpinBox, NoWheelSpinBox
        from PySide6.QtWidgets import QCheckBox

        self._sort = NoWheelComboBox(self)
        for text, key in (("Start", SORT_START), ("Strand", SORT_STRAND),
                          ("Identity", SORT_IDENTITY), ("Name", SORT_NAME)):
            self._sort.addItem(text, key)
        self._sort.setToolTip("Order of the sequences in One sequence per row")
        self._sort.currentIndexChanged.connect(
            lambda _i: self.canvas.set_sort(str(self._sort.currentData())))
        self._color_by = NoWheelComboBox(self)
        for text, key in (("Strand", COLOR_STRAND), ("Identity", COLOR_IDENTITY),
                          ("MAPQ", COLOR_MAPQ)):
            self._color_by.addItem(text, key)
        self._color_by.setToolTip("Colour of the sequence bars when zoomed out: forward blue / "
                                  "reverse purple, identity red (80 %) to green (100 %), or "
                                  "MAPQ amber (0) to blue (60)")
        self._color_by.currentIndexChanged.connect(
            lambda _i: self.canvas.set_color_by(str(self._color_by.currentData())))
        self._min_identity = NoWheelDoubleSpinBox(self)
        self._min_identity.setRange(0.0, 100.0)
        self._min_identity.setDecimals(0)
        self._min_identity.setSingleStep(5.0)
        self._min_identity.setSuffix(" %")
        self._min_identity.setToolTip("Hide sequences placed with a lower identity")
        self._min_mapq = NoWheelSpinBox(self)
        self._min_mapq.setRange(0, 60)
        self._min_mapq.setToolTip("Hide placements with a lower mapping quality "
                                  "(placements without a MAPQ are kept)")
        self._hide_secondary = QCheckBox("Hide secondary placements", self)
        self._hide_secondary.setToolTip("Show only each sequence's primary placement")
        self._hidden_label = QLabel(self)
        for w_ in (self._min_identity, self._min_mapq):
            w_.valueChanged.connect(lambda _v: self._apply_filters())
        self._hide_secondary.toggled.connect(lambda _on: self._apply_filters())

        bar2 = QWidget(self)
        b2 = QHBoxLayout(bar2)
        b2.setContentsMargins(6, 0, 6, 4)
        b2.setSpacing(6)
        for text, wdg in (("Sort", self._sort), ("Colour by", self._color_by),
                          ("Min identity", self._min_identity), ("Min MAPQ", self._min_mapq)):
            b2.addWidget(QLabel(text, self))
            b2.addWidget(wdg)
            b2.addSpacing(6)
        b2.addWidget(self._hide_secondary)
        b2.addSpacing(8)
        b2.addWidget(self._hidden_label)
        b2.addStretch(1)

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

        top = QWidget(self)
        tl = QVBoxLayout(top)
        tl.setContentsMargins(0, 0, 0, 0)
        tl.setSpacing(0)
        tl.addWidget(self.overview)
        tl.addWidget(self.canvas, 1)
        split = QSplitter(Qt.Orientation.Vertical, self)
        split.addWidget(top)
        split.addWidget(self._tabs)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 2)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(self._title)
        lay.addWidget(bar)
        lay.addWidget(bar2)
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
        self._apply_filters()
        self._layout_mode.blockSignals(True)
        self._layout_mode.setCurrentIndex(0 if self.canvas.one_per_row else 1)
        self._layout_mode.blockSignals(False)
        self._fill_reads(model)
        self._fill_features(model)
        self._fill_variants(model)
        self._sync_difference_label()

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
            self.canvas.show_base(int(pos) - 1, zoom_to_letters=True, mark=True)
            self._sync_difference_label()

    def _apply_filters(self) -> None:
        self.canvas.set_filters(min_identity=self._min_identity.value() / 100.0,
                                min_mapq=self._min_mapq.value(),
                                hide_secondary=self._hide_secondary.isChecked())
        n = self.canvas.hidden_count()
        self._hidden_label.setText(f"{n:,} hidden by filters" if n else "")
        self.overview.update()

    # ---------------------------------------------------------------- differences

    def _step_difference(self, direction: int) -> None:
        m = self.canvas.model
        if m is None or not m.diff_pos:
            return
        target = step_difference(m.diff_pos, self.canvas.mark_bp, direction)
        if target is None:
            return
        self.canvas.show_base(target, zoom_to_letters=True, mark=True)
        self._sync_difference_label()
        self._select_variant_row(target + 1)

    def _sync_difference_label(self) -> None:
        m = self.canvas.model
        n = len(m.diff_pos) if m is not None else 0
        self._prev_diff.setEnabled(n > 0)
        self._next_diff.setEnabled(n > 0)
        mark = self.canvas.mark_bp
        if not n:
            self._diff_label.setText("No differences")
        elif mark is not None and mark in m.diff_pos:
            self._diff_label.setText(f"{m.diff_pos.index(mark) + 1} of {n:,}")
        else:
            self._diff_label.setText(f"{n:,} differences")

    def _select_variant_row(self, pos1: int) -> None:
        model = self._variants_view.model()
        if model is None:
            return
        for row in range(model.rowCount()):
            if model.index(row, 0).data(Qt.ItemDataRole.UserRole) == pos1:
                self._variants_view.selectRow(row)
                self._variants_view.scrollTo(model.index(row, 0))
                break

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
