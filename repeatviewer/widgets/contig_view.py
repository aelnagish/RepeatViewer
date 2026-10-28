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
import sys
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
from PySide6.QtCore import QAbstractTableModel, QEvent, QModelIndex, QPoint, QPointF, QRectF, Qt, Signal
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

LABEL_W = 150             # default width of the name column (draggable)
LABEL_W_MIN = 80
RULER_H = 24
FEATURE_ROW_H = 14
TRACK_H = 16
COVERAGE_H = 30
IDENTITY_H = 6
LOGO_H = 34
READ_FEAT_H = 10          # one row of lifted features under a sequence
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
_ARROW_FWD = QColor("#1f5fbf")
_ARROW_REV = QColor("#7b3fb5")
_MARK_FILL = QColor(31, 95, 191, 28)        # current difference column
_MARK_EDGE = QColor("#1f5fbf")
_NO_COVERAGE = QColor("#d9d9d6")
_IDENTITY_OK = QColor("#5fae66")
_IDENTITY_SOME = QColor("#e0b43a")
_IDENTITY_MANY = QColor("#d0574a")
_CONS_DIFF_MARK = QColor("#e0a020")
# Identity bar / overview: a position counts as variable when at least this
# share of the covering sequences differ (single-copy noise stays green).
VARIABLE_SHARE = 0.25
STRONG_SHARE = 0.5

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
    "intact", "changed", … or "mixed" (the details are in the tooltip)."""
    counts = {k: int(summary.get(k, 0) or 0) for k in
              (mr.STATUS_INTACT, mr.STATUS_CHANGED, mr.STATUS_PARTIAL, mr.STATUS_LOST)}
    total = sum(counts.values())
    if total == 0:
        return ""
    for status, n in counts.items():
        if n == total:
            return status
    return "mixed"


NO_COVERAGE_CALLS = ("N", "?", "-")
CONSENSUS_THRESHOLDS = ((None, "As saved"), (0.0, "Majority"), (0.5, "50 %"),
                        (0.75, "75 %"), (0.95, "95 %"), (1.0, "100 % (strict)"))


def count_matrix(reads: list, length: int) -> np.ndarray:
    """(length, 6) counts of A C G T other gap at every reference position, from the
    placed rows (shown sequences only)."""
    counts = np.zeros((max(0, int(length)), 6), dtype=np.int32)
    for r in reads:
        if getattr(r, "hidden", False):
            continue
        row = r.data.get("row", "") if isinstance(r.data, dict) else ""
        if not row:
            continue
        idx = mr.row_indices(row).astype(np.intp)
        pos = np.arange(r.start, r.start + len(idx))
        ok = (pos >= 0) & (pos < length)
        np.add.at(counts, (pos[ok], idx[ok]), 1)
    return counts


def strand_coverage(reads: list, length: int) -> tuple[np.ndarray, np.ndarray]:
    """Per-position coverage by forward and by reverse placements (shown only)."""
    fwd = np.zeros(max(0, int(length)) + 1, dtype=np.int32)
    rev = np.zeros_like(fwd)
    for r in reads:
        if getattr(r, "hidden", False):
            continue
        s, e = max(0, r.start), min(length, r.end)
        if e <= s:
            continue
        arr = rev if r.strand == "-" else fwd
        arr[s] += 1
        arr[e] -= 1
    return np.cumsum(fwd)[:length], np.cumsum(rev)[:length]


def consensus_from_counts(counts: np.ndarray, threshold: float, no_coverage: str = "N") -> str:
    """Consensus on reference coordinates at ``threshold`` (0 = majority, else the
    fewest bases reaching that share, as an IUPAC code). Deleted positions are
    '-', positions without coverage get ``no_coverage``. Insertions are not shown."""
    if counts.size == 0:
        return ""
    calls = mr.consensus_calls(counts, threshold)
    out = []
    for c in calls:
        if c == "?":
            out.append(no_coverage or "N")
        elif c == "":
            out.append("-")
        else:
            out.append(c)
    return "".join(out)


def feature_type_counts(features: list) -> list[tuple[str, int, str]]:
    """(type, number, colour of the first feature of that type), in first-seen order."""
    seen: dict[str, list] = {}
    for f in features:
        key = str(getattr(f, "feature_type", "") or "other")
        if key not in seen:
            seen[key] = [0, str(getattr(f, "color", "") or "#888888")]
        seen[key][0] += 1
    return [(k, v[0], v[1]) for k, v in seen.items()]


def contig_statistics(m: "ContigModel") -> list[tuple[str, str]]:
    """Statistics rows for one reference (like Geneious' Statistics panel)."""
    shown = [r for r in m.reads if not r.hidden]
    cov = np.asarray(m.coverage, dtype=float)
    step = max(1, m.bin_size)
    per_bp = np.repeat(cov, step)[:m.length] if step > 1 else cov[:m.length]
    covered = int((per_bp > 0).sum()) if per_bp.size else 0
    rows: list[tuple[str, str]] = [
        ("Reference length", f"{m.length:,} bp"),
        ("Sequences placed", f"{len(shown):,}" + (f" of {len(m.reads):,} (filters)"
                                                   if len(shown) < len(m.reads) else "")),
        ("Forward / reverse", f"{sum(r.strand != '-' for r in shown):,} / "
                              f"{sum(r.strand == '-' for r in shown):,}"),
        ("Reference covered", f"{covered:,} bp ({covered / max(1, m.length):.1%})"),
    ]
    if per_bp.size:
        rows.append(("Coverage mean ± SD", f"{per_bp.mean():.1f} ± {per_bp.std():.1f}"))
        rows.append(("Coverage min / max", f"{per_bp.min():.0f} / {per_bp.max():.0f}"))
    if m.detailed and m.diff_count.size and covered:
        mask = per_bp[:len(m.diff_count)] > 0
        identical = int(((m.diff_count[:len(mask)] == 0) & mask).sum())
        rows.append(("Identical sites", f"{identical:,} of {covered:,} covered "
                                        f"({identical / covered:.1%})"))
    if shown:
        ident = np.array([r.identity for r in shown], dtype=float)
        rows.append(("Identity to reference", f"mean {ident.mean():.1%}, "
                                              f"min {ident.min():.1%}, max {ident.max():.1%}"))
    rows.append(("Variant positions", f"{len(m.diff_pos):,}"))
    cons = "".join(c for c in (m.cons_on_ref or "").upper() if c in "ACGT")
    if cons:
        n = len(cons)
        freq = {b: cons.count(b) for b in "ACGT"}
        rows.append(("Consensus bases", "  ".join(f"{b} {freq[b] / n:.1%}" for b in "ACGT")))
        rows.append(("Consensus GC", f"{(freq['G'] + freq['C']) / n:.1%}"))
    totals = {k: 0 for k in (mr.STATUS_INTACT, mr.STATUS_CHANGED, mr.STATUS_PARTIAL,
                             mr.STATUS_LOST)}
    for s_ in m.feature_summary:
        for k in totals:
            totals[k] += int(s_.get(k, 0) or 0)
    if any(totals.values()):
        rows.append(("Features in sequences", ", ".join(f"{v:,} {k}" for k, v in totals.items())))
    return rows


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
    hidden: bool = False    # filtered out (identity / MAPQ / secondary) or hidden by the user
    user_hidden: bool = False   # right-click > Hide (this session only, not the file)

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
    cache: dict = field(default_factory=dict, repr=False)

    def counts(self) -> np.ndarray:
        if "counts" not in self.cache:
            self.cache["counts"] = count_matrix(self.reads, self.length) if self.detailed \
                else np.zeros((0, 6), dtype=np.int32)
        return self.cache["counts"]

    def strands(self) -> tuple[np.ndarray, np.ndarray]:
        if "strands" not in self.cache:
            self.cache["strands"] = strand_coverage(self.reads, self.length)
        return self.cache["strands"]

    def consensus_with(self, threshold: float | None, no_coverage: str) -> str:
        """Consensus on the reference: as saved (threshold None) or recomputed."""
        if threshold is None or not self.detailed:
            return self.cons_on_ref
        key = ("cons", threshold, no_coverage)
        if key not in self.cache:
            self.cache[key] = consensus_from_counts(self.counts(), threshold, no_coverage)
        return self.cache[key]

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
    INFO_H = 48          # details of the selected sequence, foot of the name column
    readClicked = Signal(int)
    readDoubleClicked = Signal(int)
    viewChanged = Signal()       # zoom or scroll (the overview strip follows)
    labelWidthChanged = Signal(int)   # the name column was resized by the user
    readContextMenu = Signal(int, object)  # sequence index (-1: none), global QPoint

    HIGHLIGHT_DIFF = "diff"
    HIGHLIGHT_ALL = "all"
    HIGHLIGHT_NONE = "none"
    # Highlighting (Geneious): what to colour, compared to what
    SHOW_DISAGREE = "disagree"
    SHOW_AGREE = "agree"
    SHOW_ALL = "all"
    SHOW_NONE = "none"
    TO_REFERENCE = "reference"
    TO_CONSENSUS = "consensus"
    # annotation names
    LABELS_INSIDE = "inside"
    LABELS_OUTSIDE = "outside"
    LABELS_NONE = "none"
    OUTSIDE_LABEL_H = 10
    WRAP_GAP = 18           # space between wrapped blocks

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.model: ContigModel | None = None
        self.label_w = LABEL_W          # name column width; drag its edge to change
        self._edge_drag = False
        self.ppb = 1.0
        self.highlight = self.HIGHLIGHT_DIFF
        self.hl_show = self.SHOW_DISAGREE
        self.hl_target = self.TO_REFERENCE
        self.hl_dots = True
        self.label_mode = self.LABELS_INSIDE
        self.type_colors: dict[str, str] = {}     # feature type -> colour chosen by the user
        self.wrap = False
        self._wrap_origin: int | None = None      # first base of the block being painted
        self.one_per_row = True
        self.selected = -1
        self.mark_bp: int | None = None     # current difference column (0-based)
        self.sort_key = SORT_START
        self.color_by = COLOR_STRAND
        self.min_identity = 0.0
        self.min_mapq = 0
        self.hide_secondary = False
        # Annotations & Tracks
        self.hidden_types: set[str] = set()
        self.show_read_features = True
        self._feat_row: list[int] = []        # row of each reference feature (-1 hidden)
        # Graphs
        self.coverage_log = True
        self.cov_low: int | None = None       # highlight coverage below
        self.cov_high: int | None = None      # highlight coverage above
        self.highlight_single_strand = False
        self.show_logo = False
        # Consensus (None = as saved by the run)
        self.consensus_threshold: float | None = None
        self.no_coverage_call = "N"
        self._lane_count = 0
        self._feature_rows: list[int] = []
        self._feature_row_count = 0
        self._colors = _base_colors()
        self.viewport().setMouseTracking(True)
        self.viewport().grabGesture(Qt.GestureType.PinchGesture)
        self.viewport().setAttribute(Qt.WidgetAttribute.WA_AcceptTouchEvents, True)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOn)
        self.horizontalScrollBar().valueChanged.connect(self.viewport().update)
        self.horizontalScrollBar().valueChanged.connect(lambda _v: self.viewChanged.emit())
        self.horizontalScrollBar().rangeChanged.connect(lambda _a, _b: self.viewChanged.emit())
        self.verticalScrollBar().valueChanged.connect(self.viewport().update)
        self.verticalScrollBar().valueChanged.connect(
            lambda _v: self.viewChanged.emit() if self.wrap_active() else None)
        self.setMinimumHeight(220)

    # ---------------------------------------------------------------- model

    def set_model(self, model: ContigModel | None) -> None:
        self.model = model
        self.selected = -1
        self.mark_bp = None
        if model is not None:
            self._repack_features()
            self.one_per_row = len(model.reads) <= ONE_PER_ROW_LIMIT
            self._apply_filters()
            self._relayout()
            self._frame_pending = True
            self._frame_if_ready()
        self.viewport().update()

    def _frame_if_ready(self) -> None:
        """Opening view (whole reference, or the region the reads cover). Waits
        until the view has its real width: framing a 100 px wide widget left the
        result squeezed to one pixel column."""
        if not getattr(self, "_frame_pending", False) or self.model is None:
            return
        if not self.isVisible() or self.viewport().width() < self.label_w + 80:
            return
        self._frame_pending = False
        self.fit()
        self.fit_to_reads()

    def reads_span(self) -> tuple[int, int] | None:
        """[start, end) covered by the shown placements, or None."""
        m = self.model
        shown = [r for r in (m.reads if m else []) if not r.hidden]
        if not shown:
            return None
        return min(r.start for r in shown), max(r.end for r in shown)

    def fit_to_reads(self) -> None:
        """Open on the part of the reference the sequences cover (Geneious-like):
        a 455 bp read on a 3 kb reference fills the view instead of a thin strip.
        Nothing changes when the reads already cover most of the reference."""
        m = self.model
        span = self.reads_span()
        if m is None or span is None:
            return
        s0, e0 = span
        if e0 - s0 >= 0.6 * m.length:
            return
        pad = max(5, int((e0 - s0) * 0.03))
        s0, e0 = max(0, s0 - pad), min(m.length, e0 + pad)
        self.ppb = min(MAX_PPB, max(1e-6, (self._data_w() - 12) / max(1, e0 - s0)))
        self._update_scrollbars()
        self.horizontalScrollBar().setValue(int(s0 * self.ppb))
        self.viewport().update()
        self.viewChanged.emit()

    # ---------------------------------------------------------------- tracks

    def _repack_features(self) -> None:
        """Rows for the shown feature types only, so hidden types leave no gaps."""
        from repeatviewer.widgets.sequence_map_viewer import _pack_feature_rows

        m = self.model
        if m is None:
            self._feat_row, self._feature_row_count = [], 0
            return
        shown = [i for i, f in enumerate(m.features)
                 if str(f.feature_type or "other") not in self.hidden_types]
        rows, count = _pack_feature_rows([m.features[i] for i in shown])
        self._feat_row = [-1] * len(m.features)
        for i, r in zip(shown, rows):
            self._feat_row[i] = r
        self._feature_rows = [r for r in self._feat_row]
        self._feature_row_count = count

    def set_hidden_types(self, types: set[str]) -> None:
        self.hidden_types = set(types)
        self._repack_features()
        self._update_scrollbars()
        self.viewport().update()
        self.viewChanged.emit()

    def set_show_read_features(self, on: bool) -> None:
        self.show_read_features = bool(on)
        self._update_scrollbars()
        self.viewport().update()

    def set_graph_options(self, *, log_scale: bool, low: int | None, high: int | None,
                          single_strand: bool, logo: bool) -> None:
        self.coverage_log = bool(log_scale)
        self.cov_low = low
        self.cov_high = high
        self.highlight_single_strand = bool(single_strand)
        self.show_logo = bool(logo)
        self._update_scrollbars()
        self.viewport().update()

    def set_consensus(self, threshold: float | None, no_coverage: str) -> None:
        self.consensus_threshold = threshold
        self.no_coverage_call = no_coverage or "N"
        self.viewport().update()

    def consensus_seq(self) -> str:
        m = self.model
        return m.consensus_with(self.consensus_threshold, self.no_coverage_call) if m else ""

    def read_features_on(self) -> bool:
        m = self.model
        return bool(m and self.show_read_features and self.one_per_row and m.features
                    and self._feature_row_count)

    def logo_on(self) -> bool:
        return bool(self.show_logo and self.model and self.model.detailed)

    # ---------------------------------------------------------------- highlighting

    def set_highlighting(self, show: str, target: str, dots: bool) -> None:
        self.hl_show = show
        self.hl_target = target
        self.hl_dots = bool(dots)
        # older single setting, still read by the legend
        self.highlight = (self.HIGHLIGHT_DIFF if show == self.SHOW_DISAGREE and dots else
                          self.HIGHLIGHT_ALL if show == self.SHOW_ALL else
                          self.HIGHLIGHT_NONE if show == self.SHOW_NONE else show)
        self.viewport().update()

    def compare_seq(self) -> str:
        """What the sequences are compared with: the reference or the consensus."""
        m = self.model
        if m is None:
            return ""
        return self.consensus_seq() if self.hl_target == self.TO_CONSENSUS else m.reference

    def diff_positions(self) -> list[int]:
        """Positions (0-based) where a shown sequence differs from the comparison; the
        toolbar's ◀ ▶ steps through these."""
        m = self.model
        if m is None:
            return []
        if self.hl_target != self.TO_CONSENSUS or not m.detailed:
            return m.diff_pos
        key = ("cons_diff", self.consensus_threshold)
        if key not in m.cache:
            cons = self.consensus_seq()
            found: set[int] = set()
            for r in m.reads:
                if r.hidden:
                    continue
                pos, _b, runs = difference_marks(r.data.get("row", ""), r.start, cons)
                found.update(pos)
                found.update(s for s, _e in runs)
            m.cache[key] = sorted(found)
        return m.cache[key]

    def _marks(self, index: int, r: "ReadView") -> tuple[list[int], list[str]]:
        """Difference marks of one sequence against the current comparison."""
        m = self.model
        if self.hl_target != self.TO_CONSENSUS or not m.detailed or not r.data.get("row"):
            return r.mism_pos, r.mism_base
        cache = m.cache.setdefault(("cons_marks", self.consensus_threshold), {})
        if index not in cache:
            pos, bases, _runs = difference_marks(r.data["row"], r.start, self.consensus_seq())
            cache[index] = (pos, bases)
        return cache[index]

    # ---------------------------------------------------------------- labels and colours

    def set_label_mode(self, mode: str) -> None:
        self.label_mode = mode
        self._update_scrollbars()
        self.viewport().update()

    def feature_row_h(self) -> int:
        return FEATURE_ROW_H + (self.OUTSIDE_LABEL_H if self.label_mode == self.LABELS_OUTSIDE
                                else 0)

    def feature_color(self, f) -> QColor:
        return QColor(self.type_colors.get(str(f.feature_type or "other"))
                      or f.color or "#888888")

    def set_type_colors(self, colors: dict[str, str]) -> None:
        self.type_colors = dict(colors)
        self.viewport().update()
        self.viewChanged.emit()

    # ---------------------------------------------------------------- wrap

    def set_wrap(self, on: bool) -> None:
        keep = self.view_center_bp()
        self.wrap = bool(on)
        self._update_scrollbars()
        if keep is not None:
            self.show_base(keep)
        self.viewport().update()
        self.viewChanged.emit()

    def wrap_active(self) -> bool:
        """Wrapping only applies while the bases are written out."""
        return bool(self.wrap and self.model is not None and self.letters())

    def _wrap_block_bp(self) -> int:
        return max(10, int((self._data_w() - 8) / self.ppb))

    def _wrap_block_h(self) -> int:
        return self.header_h() + self._lane_count * self.read_h() + 10 + self.WRAP_GAP

    def _wrap_blocks(self) -> int:
        m = self.model
        return 0 if m is None else -(-m.length // self._wrap_block_bp())

    def view_center_bp(self) -> int | None:
        span = self.visible_span()
        return None if span is None else (span[0] + span[1]) // 2

    def visible_span(self) -> tuple[int, int] | None:
        """[first, last) reference positions on screen (wrapped: all visible blocks)."""
        m = self.model
        if m is None:
            return None
        if self.wrap_active():
            bh, bbp = self._wrap_block_h(), self._wrap_block_bp()
            v = self.verticalScrollBar().value()
            first = v // bh
            last = (v + self.viewport().height()) // bh
            return first * bbp, min(m.length, (last + 1) * bbp)
        return (max(0, self.bp_at(self.label_w)),
                min(m.length, self.bp_at(self.viewport().width())))

    def _wrap_scroll_to(self, bp: int, lane: int | None = None) -> None:
        """Scroll so the block holding ``bp`` (and the given sequence row) is in view."""
        bh = self._wrap_block_h()
        block = max(0, int(bp) // self._wrap_block_bp())
        y = block * bh
        if lane is not None and lane >= 0:
            y += self._track_y()["reads"] + lane * self.read_h() - self.viewport().height() // 2
        v = self.verticalScrollBar()
        top, bottom = v.value(), v.value() + self.viewport().height()
        block_top, block_bottom = block * bh, block * bh + self.header_h() + 40
        if lane is None and top <= block_top and block_bottom <= bottom:
            return
        v.setValue(max(0, min(v.maximum(), y)))

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
            r.hidden = r.user_hidden or not read_passes(
                r.data, r.identity, min_identity=self.min_identity, min_mapq=self.min_mapq,
                hide_secondary=self.hide_secondary)
            if r.hidden:
                self.selection.discard(i)
                if i == self._selected:
                    self._selected = -1
        m.cache.clear()    # counts, strands and recomputed consensus follow the filters

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

    def read_base_h(self) -> int:
        """Height of the sequence bar itself."""
        return READ_H_LETTERS if (self.letters() or self.one_per_row) else READ_H

    def read_h(self) -> int:
        """Lane height: the bar plus its lifted features when they are shown."""
        extra = self._feature_row_count * READ_FEAT_H + 3 if self.read_features_on() else 0
        return self.read_base_h() + extra

    def header_h(self) -> int:
        detailed = bool(self.model and self.model.detailed)
        n_tracks = 2 if detailed else 0   # consensus + reference
        return (RULER_H + max(1, self._feature_row_count) * self.feature_row_h() + 6
                + n_tracks * TRACK_H + (IDENTITY_H + 3 if detailed else 0) + COVERAGE_H + 8
                + (LOGO_H + 4 if self.logo_on() else 0))

    def _track_y(self) -> dict[str, int]:
        y = RULER_H
        out = {"features": y}
        y += max(1, self._feature_row_count) * self.feature_row_h() + 6
        if self.model and self.model.detailed:
            out["consensus"] = y
            y += TRACK_H
            out["identity"] = y + 1
            y += IDENTITY_H + 3
        out["coverage"] = y
        y += COVERAGE_H + 4
        if self.logo_on():
            out["logo"] = y
            y += LOGO_H + 4
        if self.model and self.model.detailed:
            out["reference"] = y
            y += TRACK_H
        out["reads"] = y + 4
        return out

    def _data_w(self) -> int:
        return max(1, self.viewport().width() - self.label_w)

    def _update_scrollbars(self) -> None:
        m = self.model
        if m is None:
            return
        if self.wrap_active():
            # blocks stacked down the page: no sideways scrolling
            h = self.horizontalScrollBar()
            h.setRange(0, 0)
            total = self._wrap_blocks() * self._wrap_block_h() + self.INFO_H
            v = self.verticalScrollBar()
            v.setRange(0, max(0, total - self.viewport().height()))
            v.setPageStep(max(1, self.viewport().height()))
            v.setSingleStep(max(1, self.read_h()))
            return
        content = int(m.length * self.ppb) + 20
        h = self.horizontalScrollBar()
        h.setRange(0, max(0, content - self._data_w()))
        h.setPageStep(self._data_w())
        h.setSingleStep(max(1, int(10 * self.ppb)))
        reads_h = self._lane_count * self.read_h() + 10 + self.INFO_H
        v = self.verticalScrollBar()
        v.setRange(0, max(0, reads_h - (self.viewport().height() - self.header_h())))
        v.setPageStep(max(1, self.viewport().height() - self.header_h()))
        v.setSingleStep(self.read_h())

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._update_scrollbars()
        self._frame_if_ready()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._frame_if_ready()

    def x_of(self, bp: float) -> float:
        if self._wrap_origin is not None:
            return self.label_w + (bp - self._wrap_origin) * self.ppb
        return self.label_w + bp * self.ppb - self.horizontalScrollBar().value()

    def bp_at(self, x: float) -> int:
        if self._wrap_origin is not None:
            return int(self._wrap_origin + (x - self.label_w) / self.ppb)
        return int((x - self.label_w + self.horizontalScrollBar().value()) / self.ppb)

    def _wrap_locate(self, pos: QPointF) -> tuple[int, float] | None:
        """Wrapped view: (block, y inside the block) under ``pos``."""
        bh = self._wrap_block_h()
        y = pos.y() + self.verticalScrollBar().value()
        block = int(y // bh)
        if block < 0 or block >= self._wrap_blocks():
            return None
        return block, y - block * bh

    def bp_at_pos(self, pos: QPointF) -> int:
        """Reference position under the mouse, wrapped or not."""
        if self.wrap_active():
            loc = self._wrap_locate(pos)
            if loc is None:
                return -1
            return int(loc[0] * self._wrap_block_bp() + (pos.x() - self.label_w) / self.ppb)
        return self.bp_at(pos.x())

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
        anchor_x = self.label_w + self._data_w() / 2 if anchor_x is None else anchor_x
        if self.wrap_active():
            keep = self.view_center_bp()
            fit = (self._data_w() - 12) / max(1, self.model.length)
            self.ppb = min(MAX_PPB, max(fit, self.ppb * factor))
            self._update_scrollbars()
            if self.wrap_active() and keep is not None:
                self._wrap_scroll_to(keep)
            elif keep is not None:
                self.horizontalScrollBar().setValue(int(keep * self.ppb - self._data_w() / 2))
            self.viewport().update()
            self.viewChanged.emit()
            return
        was_letters = self.letters()
        bp = (anchor_x - self.label_w + self.horizontalScrollBar().value()) / self.ppb
        fit = (self._data_w() - 12) / max(1, self.model.length)
        self.ppb = min(MAX_PPB, max(fit, self.ppb * factor))
        self._update_scrollbars()
        if not was_letters and self.wrap_active():
            self._wrap_scroll_to(int(bp))          # just reached the bases: blocks start
        else:
            self.horizontalScrollBar().setValue(int(bp * self.ppb - (anchor_x - self.label_w)))
        self.viewport().update()
        self.viewChanged.emit()

    def show_base(self, bp: int, *, zoom_to_letters: bool = False, mark: bool = False) -> None:
        if mark:
            self.mark_bp = int(bp)
        if zoom_to_letters and self.ppb < LETTER_PPB:
            self.ppb = 12.0
            self._update_scrollbars()
        if self.wrap_active():
            self._wrap_scroll_to(bp)
            self.viewport().update()
            self.viewChanged.emit()
            return
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
        if self.wrap_active():
            self._wrap_scroll_to(r.start, lane=r.lane)
            self.viewport().update()
            return
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

    def _zoom_at(self, factor: float, x: float) -> None:
        if factor > 0 and abs(factor - 1.0) > 1e-4:
            self.zoom(min(4.0, max(0.25, factor)), x)

    def event(self, event) -> bool:  # noqa: N802
        if self._gesture(event):
            return True
        return super().event(event)

    def viewportEvent(self, event) -> bool:  # noqa: N802
        if self._gesture(event):
            return True
        return super().viewportEvent(event)

    def _gesture(self, event) -> bool:
        """Trackpad pinch zooms around the fingers.

        macOS delivers pinch as native gestures (Begin / Zoom / End); the begin and end
        events must be accepted too, or the zoom steps never arrive. Platforms that only
        offer Qt's gesture framework send a QPinchGesture instead."""
        et = event.type()
        if et == QEvent.Type.NativeGesture:
            kind = event.gestureType()
            if kind == Qt.NativeGestureType.ZoomNativeGesture:
                x = self.viewport().mapFromGlobal(event.globalPosition().toPoint()).x()
                self._zoom_at(1.0 + float(event.value()), x)
            elif kind == Qt.NativeGestureType.SmartZoomNativeGesture:
                self.fit()
            elif kind not in (Qt.NativeGestureType.BeginNativeGesture,
                              Qt.NativeGestureType.EndNativeGesture):
                return False
            event.accept()
            return True
        if et == QEvent.Type.Gesture:
            pinch = event.gesture(Qt.GestureType.PinchGesture)
            if pinch is None:
                return False
            x = self.viewport().mapFromGlobal(pinch.centerPoint().toPoint()).x()
            self._zoom_at(float(pinch.scaleFactor()), x)
            event.accept(pinch)
            return True
        return False

    def wheelEvent(self, event) -> None:  # noqa: N802
        """Mouse wheel: zoom around the mouse (Shift: up/down through the sequences,
        Option/Alt: along the sequence). A trackpad scrolls as usual (two fingers up/down
        and sideways) and zooms with a pinch or Cmd/Ctrl."""
        mods = event.modifiers()
        d = event.angleDelta()
        step = d.y() or d.x()           # macOS turns Shift + wheel into a sideways delta
        trackpad = event.phase() != Qt.ScrollPhase.NoScrollPhase
        zoom_key = mods & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier)
        if zoom_key or (not trackpad and not mods and d.x() == 0):
            # one notch = 25 %; small trackpad steps zoom smoothly
            delta = step or (event.pixelDelta().y() or event.pixelDelta().x()) * 2
            if delta:
                self._zoom_at(1.25 ** (delta / 120.0), event.position().x())
            event.accept()
            return
        if not trackpad and mods & Qt.KeyboardModifier.ShiftModifier:
            v = self.verticalScrollBar()
            v.setValue(v.value() - step)
            event.accept()
            return
        if not trackpad and mods & Qt.KeyboardModifier.AltModifier:
            h = self.horizontalScrollBar()
            h.setValue(h.value() - step)
            event.accept()
            return
        super().wheelEvent(event)

    # ---------------------------------------------------------------- hit testing

    def read_at(self, pos: QPointF) -> int:
        m = self.model
        if m is None:
            return -1
        top = self._track_y()["reads"]
        if self.wrap_active():
            loc = self._wrap_locate(pos)
            if loc is None:
                return -1
            y = loc[1] - top
            if y < 0:
                return -1
        else:
            y = pos.y() - top + self.verticalScrollBar().value()
            if pos.y() < self.header_h() or y < 0:
                return -1
        lane = int(y // self.read_h())
        bp = self.bp_at_pos(pos)
        for i, r in enumerate(m.reads):
            if r.hidden:
                continue
            if r.lane == lane and r.start <= bp < r.end:
                return i
            if r.lane == lane and self.one_per_row and pos.x() < self.label_w:
                return i
        return -1

    EDGE_GRIP = 4      # px either side of the name column's edge that grab it

    def _on_edge(self, x: float) -> bool:
        return abs(x - self.label_w) <= self.EDGE_GRIP

    def set_label_width(self, width: int, *, notify: bool = False) -> None:
        """Name column width, kept between LABEL_W_MIN and half the view."""
        hi = max(LABEL_W_MIN, int(self.viewport().width() * 0.5)) if self.isVisible() else 600
        width = int(min(hi, max(LABEL_W_MIN, width)))
        if width == self.label_w:
            return
        self.label_w = width
        self._update_scrollbars()
        self.viewport().update()
        self.viewChanged.emit()
        if notify:
            self.labelWidthChanged.emit(width)

    def fit_label_width(self) -> None:
        """Wide enough for the longest name (and the track labels), within the limits."""
        m = self.model
        f = QFont(self.font())
        f.setPointSize(9)
        fm = QFontMetrics(f)
        names = [r.label for r in (m.reads if m else []) if not r.hidden]
        widest = max([fm.horizontalAdvance(n) + 34 for n in names]
                     + [fm.horizontalAdvance("Coverage (max 10,000)") + 16, LABEL_W_MIN])
        self.set_label_width(widest, notify=True)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if self._edge_drag:
            self._edge_drag = False
            self.labelWidthChanged.emit(self.label_w)
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802
        if not self._edge_drag:
            self.viewport().unsetCursor()
        super().leaveEvent(event)

    # ---------------------------------------------------------------- selection

    @property
    def selected(self) -> int:
        """The current (last clicked) sequence, or -1."""
        return getattr(self, "_selected", -1)

    @selected.setter
    def selected(self, index: int) -> None:
        self._selected = int(index)
        self.selection = {self._selected} if self._selected >= 0 else set()

    def display_order(self) -> list[int]:
        """Shown sequences top to bottom (lane, then start), for Shift-click ranges."""
        m = self.model
        if m is None:
            return []
        return sorted((i for i, r in enumerate(m.reads) if not r.hidden),
                      key=lambda i: (m.reads[i].lane, m.reads[i].start))

    def select_at(self, i: int, *, add: bool = False, extend: bool = False) -> None:
        """Click: select one. Cmd/Ctrl-click: add or remove. Shift-click: the range
        from the current one."""
        if extend and self._selected >= 0:
            order = self.display_order()
            if i in order and self._selected in order:
                a, b = sorted((order.index(self._selected), order.index(i)))
                self.selection = set(self.selection) | set(order[a:b + 1])
                self.viewport().update()
                return
        if add:
            if i in self.selection:
                self.selection.discard(i)
                self._selected = next(iter(self.selection), -1) if self._selected == i \
                    else self._selected
            else:
                self.selection.add(i)
                self._selected = i
        else:
            self.selected = i
        self.viewport().update()

    def hide_reads(self, indices) -> None:
        """Right-click > Hide: out of the view until Show hidden (not the file)."""
        m = self.model
        if m is None:
            return
        for i in indices:
            if 0 <= i < len(m.reads):
                m.reads[i].user_hidden = True
        self.selected = -1
        self._apply_filters()
        self._relayout()
        self.viewport().update()

    def show_hidden_reads(self) -> None:
        m = self.model
        if m is None:
            return
        for r in m.reads:
            r.user_hidden = False
        self._apply_filters()
        self._relayout()
        self.viewport().update()

    def user_hidden_count(self) -> int:
        return sum(1 for r in self.model.reads if r.user_hidden) if self.model else 0

    def contextMenuEvent(self, event) -> None:  # noqa: N802
        i = self.read_at(QPointF(event.pos()))
        if i >= 0 and i not in self.selection:
            self.selected = i
            self.viewport().update()
        self.readContextMenu.emit(i, event.globalPos())
        event.accept()

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton and self._on_edge(event.position().x()):
            self._edge_drag = True
            event.accept()
            return
        i = self.read_at(event.position())
        mods = event.modifiers()
        if i >= 0:
            if event.button() == Qt.MouseButton.RightButton:
                if i not in self.selection:     # right-click outside the selection
                    self.selected = i
                    self.viewport().update()
            else:
                self.select_at(
                    i, add=bool(mods & (Qt.KeyboardModifier.ControlModifier
                                        | Qt.KeyboardModifier.MetaModifier)),
                    extend=bool(mods & Qt.KeyboardModifier.ShiftModifier))
            self.readClicked.emit(i)
        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        if self._on_edge(event.position().x()):
            self.fit_label_width()
            event.accept()
            return
        i = self.read_at(event.position())
        if i < 0:
            self.fit()
        else:
            self.readDoubleClicked.emit(i)
        super().mouseDoubleClickEvent(event)

    def show_range(self, start0: int, end0: int) -> None:
        """Zoom so [start0, end0) fills the view (a single position: zoom to bases)."""
        m = self.model
        if m is None:
            return
        start0, end0 = max(0, start0), min(m.length, max(end0, start0 + 1))
        if end0 - start0 <= 1:
            self.show_base(start0, zoom_to_letters=True, mark=True)
            return
        self.ppb = min(MAX_PPB, max(1e-6, (self._data_w() - 12) / (end0 - start0)))
        self._update_scrollbars()
        # centred: a short range shown at the maximum zoom sits in the middle
        mid = (start0 + end0) / 2 * self.ppb
        self.horizontalScrollBar().setValue(int(mid - self._data_w() / 2))
        self.viewport().update()
        self.viewChanged.emit()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        x = event.position().x()
        if self._edge_drag:
            self.set_label_width(int(x))
            return
        vp = self.viewport()
        if self._on_edge(x):
            vp.setCursor(Qt.CursorShape.SplitHCursor)
            QToolTip.showText(event.globalPosition().toPoint(),
                              "Drag to resize the names; double-click to fit the longest",
                              self)
            return
        vp.unsetCursor()
        m = self.model
        if m is not None and x < self.label_w:
            # the name column: the full name of the sequence under the mouse
            i = self.read_at(event.position())
            if i >= 0:
                QToolTip.showText(event.globalPosition().toPoint(), m.reads[i].label, self)
            else:
                QToolTip.hideText()
            return
        if m is None:
            QToolTip.hideText()
            return
        bp = self.bp_at_pos(event.position())
        if not 0 <= bp < m.length:
            QToolTip.hideText()
            return
        lines = [f"Position {bp + 1:,}"]
        if m.reference and bp < len(m.reference):
            ctx = mr.methylation_context(m.reference, bp)
            if ctx:
                lines[0] += f" · {ctx} context"
        if m.reference and bp < len(m.reference):
            lines[0] += f" · reference {m.reference[bp]}"
        cons = self.consensus_seq()
        if cons and bp < len(cons):
            lines[0] += f" · consensus {cons[bp]}"
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
        here = [k for k, f in enumerate(m.features) if f.start <= bp + 1 <= f.end
                and (k >= len(self._feat_row) or self._feat_row[k] >= 0)]
        for k in here[:6]:
            f = m.features[k]
            s_ = m.feature_summary[k] if k < len(m.feature_summary) else {}
            counts = ", ".join(f"{int(s_.get(st, 0) or 0)} {st}" for st in
                               (mr.STATUS_INTACT, mr.STATUS_CHANGED, mr.STATUS_PARTIAL,
                                mr.STATUS_LOST) if int(s_.get(st, 0) or 0))
            lines.append(f"{f.name} ({f.start:,}–{f.end:,})" + (f": {counts}" if counts else ""))
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
            for fi, status, ndiff, ident, _s, _e in d.get("lift", []):
                if 0 <= fi < len(m.features):
                    f = m.features[fi]
                    if f.start <= bp + 1 <= f.end and status != mr.STATUS_NOT_COVERED:
                        extra = f", {ndiff} differences" if status == mr.STATUS_CHANGED else ""
                        lines.append(f"{f.name} here: {status}{extra}, identity {ident:.1%}")
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
        mid = pal.color(pal.ColorRole.Mid)
        p.fillRect(self.viewport().rect(), bg)
        m = self.model
        if m is None:
            p.setPen(mid)
            p.drawText(self.viewport().rect(), Qt.AlignmentFlag.AlignCenter, "No mapping loaded")
            return
        w, h = self.viewport().width(), self.viewport().height()
        if self.wrap_active():
            self._paint_wrapped(p, w, h)
        else:
            bp0 = max(0, self.bp_at(self.label_w) - 1)
            bp1 = min(m.length, self.bp_at(w) + 2)
            self._paint_page(p, w, h, bp0, bp1, self.verticalScrollBar().value())
        self._paint_info(p, h)

    def info_lines(self) -> tuple[str, str, str]:
        """Three short lines for the foot of the name column: the selected sequence, or
        how many are shown."""
        m = self.model
        if m is None:
            return "", "", ""
        n_sel = len(self.selection)
        if n_sel > 1:
            return f"{n_sel:,} sequences selected", "Right-click to copy,", "extract, hide or BLAST"
        if 0 <= self._selected < len(m.reads):
            r = m.reads[self._selected]
            mapq = r.data.get("mapq", -1)
            third = "reverse" if r.strand == "-" else "forward"
            if mapq is not None and int(mapq) >= 0:
                third += f" · mapping quality {int(mapq)}"
            return (r.label.replace(" (reversed)", ""),
                    f"{r.start + 1:,}–{r.end:,} · identity {r.identity:.1%}", third)
        shown = len(m.reads) - self.hidden_count()
        first = (f"{shown:,} of {m.n_reads:,} shown" if shown < m.n_reads
                 else f"{m.n_reads:,} sequences")
        hidden = self.user_hidden_count()
        if hidden:
            return first, f"{hidden:,} hidden", "Right-click to show them"
        return first, "Click a sequence", "for its details"

    def _paint_info(self, p: QPainter, h: int) -> None:
        p.setClipping(False)
        pal = self.palette()
        rect = QRectF(0, h - self.INFO_H, self.label_w - 1, self.INFO_H)
        p.fillRect(rect, pal.color(pal.ColorRole.Window))
        p.setPen(QPen(pal.color(pal.ColorRole.Mid), 1))
        p.drawLine(QPointF(0, rect.top()), QPointF(self.label_w - 1, rect.top()))
        first, second, third = self.info_lines()
        f = QFont(self.font())
        f.setPointSize(9)
        small = QFont(f)
        small.setPointSize(8)
        width = self.label_w - 12
        p.setPen(pal.color(pal.ColorRole.Text))
        p.setFont(f)
        p.drawText(QRectF(6, rect.top() + 2, width, 16),
                   Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                   QFontMetrics(f).elidedText(first, Qt.TextElideMode.ElideMiddle, width))
        p.setPen(pal.color(pal.ColorRole.PlaceholderText))
        p.setFont(small)
        for k, text in enumerate((second, third)):
            p.drawText(QRectF(6, rect.top() + 17 + 14 * k, width, 14),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                       QFontMetrics(small).elidedText(text, Qt.TextElideMode.ElideRight, width))

    def _paint_wrapped(self, p: QPainter, w: int, h: int) -> None:
        """Wrap contig: the alignment in blocks one view wide, stacked down the page,
        each with its own ruler, tracks and sequences."""
        m = self.model
        bh, bbp = self._wrap_block_h(), self._wrap_block_bp()
        v = self.verticalScrollBar().value()
        pal = self.palette()
        try:
            for block in range(max(0, v // bh), self._wrap_blocks()):
                top = block * bh - v
                if top > h:
                    break
                bp0 = block * bbp
                bp1 = min(m.length, bp0 + bbp)
                self._wrap_origin = bp0
                p.save()
                p.translate(0, top)
                self._paint_page(p, w, bh - self.WRAP_GAP, bp0, bp1, 0)
                p.restore()
                # gap between blocks
                p.setClipping(False)
                p.fillRect(QRectF(0, top + bh - self.WRAP_GAP, w, self.WRAP_GAP),
                           pal.color(pal.ColorRole.Base))
                p.setPen(QPen(pal.color(pal.ColorRole.Midlight), 1))
                y_sep = top + bh - self.WRAP_GAP / 2
                p.drawLine(QPointF(self.label_w, y_sep), QPointF(w, y_sep))
        finally:
            self._wrap_origin = None

    def _paint_page(self, p: QPainter, w: int, h: int, bp0: int, bp1: int,
                    scroll: int) -> None:
        """Ruler, annotations, consensus / identity / coverage / logo / reference, the
        sequences (scrolled by ``scroll``) and the label column, for bp0..bp1."""
        m = self.model
        pal = self.palette()
        fg = pal.color(pal.ColorRole.Text)
        mid = pal.color(pal.ColorRole.Mid)
        ys = self._track_y()
        letters = self.letters()
        font = QFont(p.font())
        font.setPointSize(8)
        p.setFont(font)
        header = self.header_h()
        p.setClipRect(QRectF(self.label_w, 0, w - self.label_w, header))
        self._paint_ruler(p, bp0, bp1, fg, mid)
        self._paint_features(p, ys["features"], bp0, bp1)
        if m.detailed:
            to_cons = self.hl_target == self.TO_CONSENSUS
            cons = self.consensus_seq()
            self._paint_sequence_track(p, ys["consensus"], cons, bp0, bp1, letters,
                                       compare="" if to_cons else m.reference, fg=fg,
                                       consensus_track=True)
            self._paint_identity(p, ys["identity"], bp0, bp1, w)
            self._paint_sequence_track(p, ys["reference"], m.reference, bp0, bp1, letters,
                                       compare=cons if to_cons else "", fg=fg,
                                       consensus_track=False)
        self._paint_coverage(p, ys["coverage"], bp0, bp1, w)
        if "logo" in ys:
            self._paint_logo(p, ys["logo"], bp0, bp1, fg)
        p.setClipRect(QRectF(self.label_w, header, w - self.label_w, h - header))
        self._paint_reads(p, ys["reads"], bp0, bp1, letters, fg, h, scroll)
        if self.mark_bp is not None and bp0 <= self.mark_bp <= bp1:
            p.setClipRect(QRectF(self.label_w, RULER_H, w - self.label_w, h - RULER_H))
            x = self.x_of(self.mark_bp)
            col = QRectF(x, RULER_H, max(2.0, self.ppb), h - RULER_H)
            p.fillRect(col, _MARK_FILL)
            p.setPen(QPen(_MARK_EDGE, 1))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(col)
        # label column
        p.setClipRect(QRectF(0, 0, w, h))
        p.fillRect(QRectF(0, 0, self.label_w, h), pal.color(pal.ColorRole.Window))
        p.setPen(QPen(mid, 1))
        p.drawLine(self.label_w - 1, 0, self.label_w - 1, h)
        p.drawLine(0, header - 2, w, header - 2)
        p.setPen(fg)
        label_font = QFont(p.font())
        label_font.setPointSize(9)
        p.setFont(label_font)

        def label(y: float, text: str, height: float = TRACK_H) -> None:
            p.drawText(QRectF(6, y, self.label_w - 10, height),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, text)

        if self._wrap_origin is not None:
            label(0, f"{bp0 + 1:,}–{bp1:,}", RULER_H)
        else:
            label(0, f"{m.length:,} bp", RULER_H)
        label(ys["features"], "Annotations" if m.features else "No annotations",
              self.feature_row_h())
        if m.detailed:
            label(ys["consensus"], "Consensus")
            small = QFont(label_font)
            small.setPointSize(7)
            p.setFont(small)
            p.drawText(QRectF(6, ys["identity"] - 3, self.label_w - 10, IDENTITY_H + 6),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, "Identity")
            p.setFont(label_font)
            fm_ref = QFontMetrics(label_font)
            label(ys["reference"], fm_ref.elidedText(f"Reference · {m.ref_id}",
                                                     Qt.TextElideMode.ElideMiddle, self.label_w - 12))
        top = max(1.0, float(m.coverage[max(0, bp0 // m.bin_size):
                                         max(1, bp1 // m.bin_size + 1)].max(initial=0)))
        label(ys["coverage"], f"Coverage (max {top:.0f})", COVERAGE_H)
        if "logo" in ys:
            label(ys["logo"], "Logo", LOGO_H)
        if self.one_per_row:
            p.setClipRect(QRectF(0, header, self.label_w - 1, h - header))
            rh = self.read_h()
            off = ys["reads"] - scroll
            fm = QFontMetrics(label_font)
            for i, r in enumerate(m.reads):
                if r.hidden:
                    continue
                y = off + r.lane * rh
                if y < header - rh or y > h:
                    continue
                if i in self.selection:
                    p.fillRect(QRectF(0, y, self.label_w - 1, rh), _READ_FILL_SEL)
                elif r.strand == "-":
                    p.fillRect(QRectF(0, y, self.label_w - 1, rh), _REVERSE_TINT)
                bh_ = self.read_base_h()
                rev = r.strand == "-"
                p.setPen(_ARROW_REV if rev else _ARROW_FWD)
                p.drawText(QRectF(6, y, 14, bh_), Qt.AlignmentFlag.AlignCenter,
                           "←" if rev else "→")
                p.setPen(fg)
                pl = int(r.data.get("placement", 1) or 1)
                name = r.name if pl == 1 else f"{r.name} #{pl}"
                p.drawText(QRectF(24, y, self.label_w - 28, bh_),
                           Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                           fm.elidedText(name, Qt.TextElideMode.ElideMiddle, self.label_w - 30))
                if self.read_features_on() and y + rh < h:
                    # light rule between sequences when their features are shown
                    p.setPen(QPen(QColor(0, 0, 0, 28), 1))
                    p.drawLine(QPointF(0, y + rh - 1), QPointF(w, y + rh - 1))
        else:
            # Packed: several sequences share a line, so only the selected ones are named
            p.setClipRect(QRectF(0, header, self.label_w - 1, h - header))
            rh = self.read_h()
            off = ys["reads"] - scroll
            fm = QFontMetrics(label_font)
            bh_ = self.read_base_h()
            for i in sorted(self.selection)[:200]:
                if not 0 <= i < len(m.reads) or m.reads[i].hidden:
                    continue
                r = m.reads[i]
                y = off + r.lane * rh
                if y < header - rh or y > h:
                    continue
                p.fillRect(QRectF(0, y, self.label_w - 1, rh), _READ_FILL_SEL)
                rev = r.strand == "-"
                p.setPen(_ARROW_REV if rev else _ARROW_FWD)
                p.drawText(QRectF(6, y, 14, bh_), Qt.AlignmentFlag.AlignCenter,
                           "←" if rev else "→")
                p.setPen(fg)
                p.drawText(QRectF(24, y, self.label_w - 28, bh_),
                           Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                           fm.elidedText(r.label.replace(" (reversed)", ""),
                                         Qt.TextElideMode.ElideMiddle, self.label_w - 30))
        p.setClipping(False)

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
        m = self.model
        summary = m.feature_summary
        row_h = self.feature_row_h()
        outside = self.label_mode == self.LABELS_OUTSIDE
        pal = self.palette()
        text_c = pal.color(pal.ColorRole.Text)
        next_start: dict[int, int] = {}
        if outside:
            # an outside label may run past its bar, up to the next annotation in its row
            by_row: dict[int, list[int]] = {}
            for k, row in enumerate(self._feat_row):
                if row >= 0:
                    by_row.setdefault(row, []).append(k)
            for ks in by_row.values():
                ks.sort(key=lambda k: m.features[k].start)
                for a_, b_ in zip(ks, ks[1:]):
                    next_start[a_] = m.features[b_].start - 1
        for k, (f, row) in enumerate(zip(m.features, self._feat_row)):
            if row < 0 or f.end < bp0 or f.start - 1 > bp1:
                continue
            x0, x1 = self.x_of(f.start - 1), self.x_of(f.end)
            top = y0 + row * row_h
            bar_top = top + (self.OUTSIDE_LABEL_H if outside else 0) + 1
            rect = QRectF(x0, bar_top, max(2.0, x1 - x0), FEATURE_ROW_H - 3)
            status = feature_status_text(summary[k]) if k < len(summary) else ""
            text = f"{f.name} · {status}" if status else f.name
            _feature_bar(p, rect, self.feature_color(f), forward=f.is_forward(),
                         label=text if self.label_mode == self.LABELS_INSIDE else "",
                         point_size=7, left_edge=self.label_w, right_edge=self.viewport().width())
            if outside:
                limit = self.x_of(next_start[k]) - 6 if k in next_start else \
                    self.viewport().width()
                room = QRectF(rect.left(), rect.top(), max(rect.width(), limit - rect.left()),
                              rect.height())
                _outside_label(p, room, text, top, self.OUTSIDE_LABEL_H + 1, text_c,
                               left_edge=self.label_w, right_edge=self.viewport().width())

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

    def _base_style(self, differs: bool) -> str:
        """How one compared base is drawn under the Highlighting settings:
        'colour' (filled), 'dot', or 'plain' (letter only)."""
        show = self.hl_show
        if show == self.SHOW_ALL:
            return "colour"
        if show == self.SHOW_NONE:
            return "plain"
        if show == self.SHOW_AGREE:
            return "plain" if differs else "colour"
        # disagreements
        if differs:
            return "colour"
        return "dot" if self.hl_dots else "plain"

    def _paint_sequence_track(self, p: QPainter, y: int, seq: str, bp0: int, bp1: int,
                              letters: bool, *, compare: str, fg: QColor,
                              consensus_track: bool) -> None:
        """Consensus or reference. The track that is the comparison target (compare
        == "") is drawn fully coloured; the other one is compared with it."""
        if not seq:
            return
        m = self.model
        rect_h = TRACK_H - 3
        uncovered = self._uncovered_runs(bp0, bp1) if consensus_track else []

        def covered(bp: int) -> bool:
            return not consensus_track or m.coverage_at(bp) > 0

        if not letters:
            p.fillRect(QRectF(self.x_of(0), y + 4, m.length * self.ppb, rect_h - 6),
                       QColor("#b9b9b9"))
            for s0, e0 in uncovered:
                r = QRectF(self.x_of(s0), y + 1, max(1.0, (e0 - s0) * self.ppb), rect_h)
                p.fillRect(r, pal_base(self))
                self._hatch(p, r)
            if compare and self.hl_show in (self.SHOW_DISAGREE, self.SHOW_ALL):
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
        dot_color = QColor("#9a9a96")
        for bp in range(bp0, min(bp1, len(seq))):
            rect = QRectF(self.x_of(bp), y, self.ppb, rect_h)
            if not covered(bp):
                self._hatch(p, rect)
                continue
            c = seq[bp]
            bg_c, fg_c = self._colors.get(c.upper(), _FALLBACK)
            ambiguous = c.upper() not in "ACGT-"
            if ambiguous and not consensus_track:
                # ambiguity code in the reference: pale fill, dashed outline
                p.fillRect(rect, QColor("#fff4d6"))
                p.setPen(QPen(QColor("#a35f00"), 1, Qt.PenStyle.DashLine))
                p.drawRect(rect.adjusted(0.5, 0.5, -0.5, -0.5))
                p.setPen(QColor("#6b4100"))
                p.drawText(rect, Qt.AlignmentFlag.AlignCenter, c)
                continue
            if not compare:
                style = "colour"                         # the comparison target itself
            else:
                differs = bp < len(compare) and not bases_match(compare[bp], c)
                style = self._base_style(differs)
                if style == "dot":
                    style = "plain"                      # tracks keep their letters
            if style == "colour":
                p.fillRect(rect, bg_c)
                p.setPen(fg_c)
            else:
                p.setPen(fg if compare else dot_color)
            p.drawText(rect, Qt.AlignmentFlag.AlignCenter, c)
            if (consensus_track and bp < len(m.diff_count) and m.diff_count[bp]
                    and self.hl_show != self.SHOW_NONE):
                # some sequences differ here: amber underline
                p.fillRect(QRectF(rect.left(), rect.bottom() - 2, rect.width(), 2),
                           _CONS_DIFF_MARK)

    def _paint_read_features(self, p: QPainter, r: ReadView, y0: float, bp0: int,
                             bp1: int) -> None:
        """The reference features as they came out in this sequence: intact and
        changed in the feature colour ("n diff" in the label), partial faded and
        dashed, lost in grey."""
        _paint_read_features_impl(self, p, r, y0, bp0, bp1)

    def _paint_logo(self, p: QPainter, y: int, bp0: int, bp1: int, fg: QColor) -> None:
        """Frequency logo: at each position the bases stacked by their share of the
        shown sequences, most common on top (letters once there is room)."""
        m = self.model
        counts = m.counts()
        if counts.size == 0:
            return
        hi = min(bp1, counts.shape[0])
        if hi <= bp0:
            return
        if self.ppb < LETTER_PPB:
            p.setPen(QColor("#9a9a96"))
            p.drawText(QRectF(self.x_of(bp0) if bp0 > 0 else self.label_w, y, 400, LOGO_H),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                       "  Zoom in to see the logo")
            return
        sub = counts[bp0:hi, :4].astype(float)
        depth = sub.sum(axis=1)
        letters = True
        font = _mono(9)
        p.setFont(font)
        step = 1 if self.ppb >= 1 else max(1, int(1 / self.ppb))
        for k in range(0, hi - bp0, step):
            d = depth[k]
            if d <= 0:
                continue
            x = self.x_of(bp0 + k)
            wpx = max(1.0, self.ppb * step)
            base_y = y + LOGO_H
            for b in np.argsort(sub[k]):            # smallest first, at the bottom
                share = sub[k, b] / d
                if share <= 0:
                    continue
                hgt = share * (LOGO_H - 2)
                base = "ACGT"[b]
                bg_c, fg_c = self._colors.get(base, _FALLBACK)
                rect = QRectF(x, base_y - hgt, wpx, hgt)
                p.fillRect(rect, bg_c)
                if letters and hgt >= 9:
                    p.setPen(fg_c)
                    p.drawText(rect, Qt.AlignmentFlag.AlignCenter, base)
                base_y -= hgt

    def _paint_identity(self, p: QPainter, y: int, bp0: int, bp1: int, w: int) -> None:
        """Green where (almost) all sequences agree with the reference, yellow where a
        quarter or more differ, red where half or more differ; nothing where there is
        no coverage. Single-copy differences stay green so real variants stand out."""
        m = self.model
        if m.diff_count.size == 0:
            return
        x_start = int(max(self.label_w, self.x_of(bp0)))
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
            color = _IDENTITY_OK if share < VARIABLE_SHARE else (
                _IDENTITY_SOME if share < STRONG_SHARE else _IDENTITY_MANY)
            p.fillRect(QRectF(x, y, step, IDENTITY_H), color)

    def _paint_coverage(self, p: QPainter, y: int, bp0: int, bp1: int, w: int) -> None:
        """Coverage (log scale by default, so one or two sequences do not fill the
        track), with optional highlight bands along its foot: red below the low
        threshold, orange above the high one, purple where only one strand covers."""
        import math

        m = self.model
        cov = m.coverage
        if cov.size == 0:
            return
        x_start = int(max(self.label_w, self.x_of(bp0)))
        x_end = int(min(w, self.x_of(bp1)))
        xs, vals, spans = [], [], []
        for x in range(x_start, x_end + 1, 2):
            a = max(0, self.bp_at(x))
            b = max(a + 1, self.bp_at(x + 2))
            ia, ib = a // m.bin_size, max(a // m.bin_size + 1, b // m.bin_size)
            seg = cov[ia:ib]
            xs.append(x)
            vals.append(float(seg.max()) if seg.size else 0.0)
            spans.append((a, b))
        if not vals:
            return
        top = max(1.0, max(vals))
        span = COVERAGE_H - 9
        base = y + COVERAGE_H - 5
        if self.coverage_log:
            scale = math.log1p(max(top, 4.0))
            height = [math.log1p(v) / scale * span for v in vals]
        else:
            height = [v / top * span for v in vals]
        poly = QPolygonF([QPointF(xs[0], base)] + [
            QPointF(x, base - hv) for x, hv in zip(xs, height)] + [QPointF(xs[-1], base)])
        p.setPen(QPen(QColor("#3b6fb6"), 1))
        p.setBrush(QColor("#cfe0f5"))
        p.drawPolygon(poly)
        p.setBrush(Qt.BrushStyle.NoBrush)
        band_y = base + 1
        if self.cov_low is not None or self.cov_high is not None:
            for x, v in zip(xs, vals):
                if self.cov_low is not None and v < self.cov_low:
                    p.fillRect(QRectF(x, band_y, 2, 3), QColor("#d0574a"))
                elif self.cov_high is not None and v > self.cov_high:
                    p.fillRect(QRectF(x, band_y, 2, 3), QColor("#e8913a"))
        if self.highlight_single_strand:
            fwd, rev = m.strands()
            for x, (a, b) in zip(xs, spans):
                b = min(b, len(fwd))
                if b <= a:
                    continue
                f_, r_ = int(fwd[a:b].max(initial=0)), int(rev[a:b].max(initial=0))
                if (f_ > 0) != (r_ > 0):
                    p.fillRect(QRectF(x, band_y + 3, 2, 2), QColor("#7b3fb5"))

    def _paint_reads(self, p: QPainter, top: int, bp0: int, bp1: int, letters: bool,
                     fg: QColor, h: int, scroll: int | None = None) -> None:
        from repeatviewer.widgets.sequence_map_viewer import _draw_directed_feature

        m = self.model
        rh = self.read_h()
        bh = self.read_base_h()
        feats_on = self.read_features_on()
        if scroll is None:
            scroll = self.verticalScrollBar().value()
        off = top - scroll
        header = self.header_h()
        w = self.viewport().width()
        font = _mono(int(min(11, max(7, self.ppb * 0.85))))
        bold = QFont(font)
        bold.setBold(True)
        dot_color = QColor("#9a9a96")
        ref = self.compare_seq()
        show_marks = self.hl_show in (self.SHOW_DISAGREE, self.SHOW_ALL)
        for i, r in enumerate(m.reads):
            if r.hidden:
                continue
            y = off + r.lane * rh
            if y + rh < header or y > h:
                continue
            sel = i in self.selection
            if self.one_per_row and r.strand == "-" and not sel:
                p.fillRect(QRectF(self.label_w, y, w - self.label_w, rh), _REVERSE_TINT)
            if r.end < bp0 or r.start > bp1:
                continue
            if not letters:
                rect = QRectF(self.x_of(r.start), y + 1, max(2.0, (r.end - r.start) * self.ppb),
                              bh - 3)
                _draw_directed_feature(p, rect, _READ_FILL_SEL if sel else self._read_fill(r),
                                       forward=r.strand != "-")
                if sel:
                    p.setPen(QPen(QColor("#1d1d1f"), 1.5))
                    p.setBrush(Qt.BrushStyle.NoBrush)
                    p.drawRect(rect.adjusted(-1, -1, 1, 1))
                if show_marks:
                    mism_pos, mism_base = self._marks(i, r)
                    lo = bisect.bisect_left(mism_pos, bp0)
                    hi = bisect.bisect_right(mism_pos, bp1)
                    for k in range(lo, hi):
                        p.fillRect(QRectF(self.x_of(mism_pos[k]), y + 1,
                                          max(1.0, self.ppb), bh - 3),
                                   self._colors.get(mism_base[k].upper(), _FALLBACK)[0])
                for s, e in r.del_runs:
                    if e < bp0 or s > bp1:
                        continue
                    x0, x1 = self.x_of(s), self.x_of(e)
                    p.fillRect(QRectF(x0, y + 1, max(1.0, x1 - x0), bh - 3), QColor("#ffffff"))
                    p.setPen(QPen(_DELETE_COLOR, 1))
                    p.drawLine(QPointF(x0, y + bh / 2 - 0.5), QPointF(x1, y + bh / 2 - 0.5))
            else:
                p.setFont(font)
                row = r.data["row"]
                if sel:
                    p.fillRect(QRectF(self.x_of(max(bp0, r.start)), y,
                                      (min(bp1, r.end) - max(bp0, r.start)) * self.ppb, bh - 2),
                               _READ_FILL_SEL)
                for bp in range(max(bp0, r.start), min(bp1, r.end)):
                    c = row[bp - r.start]
                    rect = QRectF(self.x_of(bp), y, self.ppb, bh - 2)
                    if c == "-":
                        # deleted base: dark cell with a dash
                        p.fillRect(rect.adjusted(0, 1, 0, -1), _DELETE_COLOR)
                        p.setPen(QColor("#ffffff"))
                        p.drawText(rect, Qt.AlignmentFlag.AlignCenter, "-")
                        continue
                    differs = bp < len(ref) and not bases_match(ref[bp], c)
                    style = self._base_style(differs)
                    if style == "colour":
                        bg_c, fg_c = self._colors.get(c.upper(), _FALLBACK)
                        p.fillRect(rect, bg_c)
                        p.setPen(fg_c)
                        if differs and self.hl_show == self.SHOW_DISAGREE:
                            p.setFont(bold)
                        p.drawText(rect, Qt.AlignmentFlag.AlignCenter, c)
                        p.setFont(font)
                    elif style == "dot":
                        # agrees: a dot, so the differences stand out
                        p.setPen(dot_color)
                        p.drawText(rect, Qt.AlignmentFlag.AlignCenter, "·")
                    else:
                        p.setPen(fg)
                        p.drawText(rect, Qt.AlignmentFlag.AlignCenter, c)
            if feats_on:
                self._paint_read_features(p, r, y + bh + 1, bp0, bp1)
            p.setPen(QPen(_INSERT_COLOR, 2))
            lo = bisect.bisect_left(r.ins_pos, bp0)
            hi = bisect.bisect_right(r.ins_pos, bp1)
            for k in range(lo, hi):
                x = self.x_of(r.ins_pos[k])
                p.drawLine(QPointF(x, y), QPointF(x, y + bh - 2))


def _feature_bar(p: QPainter, rect: QRectF, color: QColor, *, forward: bool, label: str,
                 point_size: int, left_edge: float, right_edge: float) -> None:
    """Directed bar with its label kept inside the visible part, so a long feature
    that starts off-screen still says what it is."""
    from repeatviewer.widgets.sequence_map_viewer import _draw_directed_feature

    _draw_directed_feature(p, rect, color, forward=forward, label="")
    vis_l, vis_r = max(rect.left(), left_edge), min(rect.right(), right_edge)
    if not label or vis_r - vis_l < 28:
        return
    font = QFont(p.font())
    font.setPointSize(point_size)
    p.setFont(font)
    fm = QFontMetrics(font)
    head = min(rect.height() * 0.9, 10.0)
    text_rect = QRectF(vis_l + 4, rect.top(), vis_r - vis_l - 8 - head, rect.height())
    if text_rect.width() < 12:
        return
    p.setPen(QColor("#ffffff"))
    p.drawText(text_rect, Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
               fm.elidedText(label, Qt.TextElideMode.ElideRight, int(text_rect.width())))


def _outside_label(p: QPainter, bar: QRectF, text: str, top: float, height: float,
                   color: QColor, *, left_edge: float, right_edge: float) -> None:
    """Annotation name above its bar, within the bar's visible width (left out when
    it does not fit)."""
    vis_l, vis_r = max(bar.left(), left_edge), min(bar.right(), right_edge)
    if not text or vis_r - vis_l < 24:
        return
    font = QFont(p.font())
    font.setPointSize(7)
    p.setFont(font)
    fm = QFontMetrics(font)
    shown = fm.elidedText(text, Qt.TextElideMode.ElideRight, int(vis_r - vis_l - 2))
    if shown != text and len(shown) < 6:
        return                       # too little room: leave the label out
    p.setPen(color)
    p.drawText(QRectF(vis_l + 1, top, vis_r - vis_l - 2, height),
               Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, shown)


def _paint_read_features_impl(canvas: "ContigCanvas", p: QPainter, r: "ReadView", y0: float,
                              bp0: int, bp1: int) -> None:
    m = canvas.model
    for fi, status, ndiff, _ident, _s, _e in r.data.get("lift", []):
        if not 0 <= fi < len(m.features) or status == mr.STATUS_NOT_COVERED:
            continue
        row = canvas._feat_row[fi] if fi < len(canvas._feat_row) else -1
        if row < 0:
            continue
        f = m.features[fi]
        s0, e0 = max(f.start - 1, r.start), min(f.end, r.end)
        if e0 <= s0 or e0 < bp0 or s0 > bp1:
            continue
        x0, x1 = canvas.x_of(s0), canvas.x_of(e0)
        rect = QRectF(x0, y0 + row * READ_FEAT_H, max(2.0, x1 - x0), READ_FEAT_H - 2)
        color = canvas.feature_color(f)
        if status == mr.STATUS_LOST:
            color = QColor("#c4c4c0")
        elif status == mr.STATUS_PARTIAL:
            color.setAlpha(120)
        label = f.name if status == mr.STATUS_INTACT else f"{f.name} · " + (
            f"{ndiff} diff" if status == mr.STATUS_CHANGED else status)
        if canvas.label_mode == canvas.LABELS_NONE:
            label = ""
        _feature_bar(p, rect, color, forward=f.is_forward(), label=label, point_size=6,
                     left_edge=canvas.label_w, right_edge=canvas.viewport().width())
        if status == mr.STATUS_PARTIAL:
            p.setPen(QPen(canvas.feature_color(f), 1, Qt.PenStyle.DashLine))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(rect)


LEGEND_ITEMS = (
    ("dot", "matches reference"), ("T", "differs"), ("-", "deleted"),
    ("Y", "ambiguity code"), ("hatch", "no coverage"),
    ("fwd", "forward"), ("rev", "reverse"),
    ("partial", "feature partial"), ("lost", "feature lost"),
)
_AMBIGUOUS = set("RYSWKMBDHVNryswkmbdhvn")


def legend_items(model: "ContigModel | None", highlight: str = "",
                 target: str = "reference", show: str = "disagree") -> list[tuple[str, str]]:
    """The legend entries that apply to this result: e.g. no 'feature lost' without
    annotations, no 'reverse' when nothing was reversed. Wording follows the
    Highlighting settings (compared to the reference or the consensus)."""
    if model is None:
        return []
    reads = [r for r in model.reads if not r.hidden] or list(model.reads)
    keep = set()
    if show != ContigCanvas.SHOW_NONE:
        keep.add("T")
    if highlight in ("", ContigCanvas.HIGHLIGHT_DIFF) and show == ContigCanvas.SHOW_DISAGREE:
        keep.add("dot")
    if any(r.del_runs or r.data.get("deleted", 0) for r in reads):
        keep.add("-")
    if any(ch in _AMBIGUOUS for ch in (model.consensus or "")[:2_000_000]) or any(
            not _AMBIGUOUS.isdisjoint(str(r.data.get("row", ""))) for r in reads[:5000]):
        keep.add("Y")
    cov = np.asarray(model.coverage, dtype=float)
    if cov.size and float(cov.min()) <= 0:
        keep.add("hatch")
    if any(r.strand != "-" for r in reads):
        keep.add("fwd")
    if any(r.strand == "-" for r in reads):
        keep.add("rev")
    if model.features:
        statuses = {lf[1] for r in reads for lf in r.data.get("lift", [])}
        if "partial" in statuses:
            keep.add("partial")
        if "lost" in statuses:
            keep.add("lost")
    words = {
        "dot": f"matches {target}",
        "T": (f"agrees with {target}" if show == ContigCanvas.SHOW_AGREE else
              "base colour" if show == ContigCanvas.SHOW_ALL else f"differs from {target}"),
    }
    return [(k, words.get(k, t)) for k, t in LEGEND_ITEMS if k in keep]


_BAR_FWD, _BAR_REV = "#8fb3e6", "#a98bd6"


def _gradient(mode: str):
    """Colour stops of a Colour by scale (identity or mapping quality)."""
    from PySide6.QtGui import QLinearGradient

    g = QLinearGradient(0, 0, 1, 0)
    g.setCoordinateMode(QLinearGradient.CoordinateMode.ObjectBoundingMode)
    if mode == COLOR_IDENTITY:
        g.setColorAt(0, QColor(read_color_hex("+", 0.8, None, COLOR_IDENTITY)))
        g.setColorAt(1, QColor(read_color_hex("+", 1.0, None, COLOR_IDENTITY)))
    else:
        g.setColorAt(0, QColor(read_color_hex("+", 1.0, 0, COLOR_MAPQ)))
        g.setColorAt(1, QColor(read_color_hex("+", 1.0, 60, COLOR_MAPQ)))
    return g


def bar_colour_items(mode: str) -> list[tuple[str, str]]:
    """Legend rows for the bar colour of the current Colour by."""
    if mode == COLOR_IDENTITY:
        rows = [("grad_identity", "identity ≤ 80 % → 100 %")]
    elif mode == COLOR_MAPQ:
        rows = [("grad_mapq", "mapping quality 0 → 60"), ("bar_none", "no mapping quality")]
    else:
        return []          # Strand: the forward / reverse rows already show the colours
    return [("heading", "Bar colour, zoomed out")] + rows


class ContigLegend(QWidget):
    """The key to the view's marks, one entry per row; shown from the "?" button."""

    ROW_H = 22
    PAD = 10

    def __init__(self, items: list[tuple[str, str]], parent: QWidget | None = None, *,
                 strand_bars: bool = True) -> None:
        super().__init__(parent)
        self._items = list(items)
        self._strand_bars = strand_bars     # bars are coloured by strand (Colour by)
        f = QFont(self.font())
        fm = QFontMetrics(f)
        text_w = max((fm.horizontalAdvance(t) for _k, t in self._items), default=80)
        self.setFixedSize(self.PAD * 2 + 14 + 8 + text_w + 4,
                          self.PAD * 2 + self.ROW_H * max(1, len(self._items)))

    def paintEvent(self, event) -> None:  # noqa: N802
        from PySide6.QtGui import QBrush

        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        pal = self.palette()
        f = QFont(self.font())
        p.setFont(f)
        text_c = pal.color(pal.ColorRole.WindowText)
        mono = _mono(9)
        colors = _base_colors()
        sw = 14
        for i, (kind, text) in enumerate(self._items):
            y = self.PAD + i * self.ROW_H
            box = QRectF(self.PAD, y + (self.ROW_H - 12) / 2, sw, 12)
            if kind == "dot":
                p.setPen(QColor("#9a9a96"))
                p.setFont(mono)
                p.drawText(box.adjusted(-3, 0, 3, 0), Qt.AlignmentFlag.AlignCenter, "···")
            elif kind in ("T", "-"):
                bg_c = colors.get("T", _FALLBACK)[0] if kind == "T" else _DELETE_COLOR
                p.fillRect(box, bg_c)
                p.setPen(QColor("#ffffff") if kind == "-" else colors.get("T", _FALLBACK)[1])
                p.setFont(mono)
                p.drawText(box, Qt.AlignmentFlag.AlignCenter, kind)
            elif kind == "Y":
                p.fillRect(box, QColor("#fff4d6"))
                p.setPen(QPen(QColor("#a35f00"), 1, Qt.PenStyle.DashLine))
                p.drawRect(box.adjusted(0.5, 0.5, -0.5, -0.5))
            elif kind == "hatch":
                p.fillRect(box, QBrush(QColor("#9a9a96"), Qt.BrushStyle.BDiagPattern))
                p.setPen(QPen(QColor("#9a9a96"), 1))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRect(box.adjusted(0.5, 0.5, -0.5, -0.5))
            elif kind in ("fwd", "rev"):
                # arrow (name column) over a bar in the strand colour (zoomed out)
                p.setPen(_ARROW_REV if kind == "rev" else _ARROW_FWD)
                p.drawText(box.adjusted(0, -3, 0, -3), Qt.AlignmentFlag.AlignCenter,
                           "←" if kind == "rev" else "→")
                if self._strand_bars:
                    p.fillRect(QRectF(box.left(), box.bottom() - 1, box.width(), 3),
                               QColor(_BAR_REV if kind == "rev" else _BAR_FWD))
            elif kind == "partial":
                c = QColor("#4caf50")
                c.setAlpha(120)
                p.fillRect(box.adjusted(0, 3, 0, -3), c)
                p.setPen(QPen(QColor("#4caf50"), 1, Qt.PenStyle.DashLine))
                p.drawRect(box.adjusted(0, 3, 0, -3))
            elif kind == "lost":
                p.fillRect(box.adjusted(0, 3, 0, -3), QColor("#c4c4c0"))
            elif kind == "heading":
                bold = QFont(f)
                bold.setBold(True)
                p.setFont(bold)
                p.setPen(text_c)
                p.drawText(QRectF(self.PAD, y, self.width(), self.ROW_H),
                           Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignLeft, text)
                continue
            elif kind in ("bar_fwd", "bar_rev", "bar_none"):
                color = {"bar_fwd": _BAR_FWD, "bar_rev": _BAR_REV, "bar_none": "#cfcfcf"}[kind]
                p.fillRect(box.adjusted(0, 3, 0, -3), QColor(color))
            elif kind in ("grad_identity", "grad_mapq"):
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(_gradient(COLOR_IDENTITY if kind == "grad_identity" else COLOR_MAPQ))
                p.drawRect(box.adjusted(0, 3, 0, -3))
                p.setBrush(Qt.BrushStyle.NoBrush)
            p.setFont(f)
            p.setPen(text_c)
            p.drawText(QRectF(self.PAD + sw + 8, y, self.width(), self.ROW_H),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, text)
        p.end()


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
        # where the view's data area starts (its label column ends)
        vp = self.canvas.viewport()
        return float(self.mapFromGlobal(vp.mapToGlobal(QPoint(self.canvas.label_w, 0))).x())

    def _track_w(self) -> float:
        # to the right edge of the view's frame below, so both boxes end on one line
        c = self.canvas
        right = self.mapFromGlobal(c.mapToGlobal(QPoint(c.width(), 0))).x()
        return max(1.0, min(float(self.width()), float(right)) - 1.0 - self._x0())

    def _scale(self) -> float:
        m = self.canvas.model
        return self._track_w() / max(1, m.length if m else 1)

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
        p.drawText(QRectF(6, 0, self.canvas.label_w - 10, self.HEIGHT),
                   Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                   f"Whole reference · {m.length:,} bp")
        x0, sc = self._x0(), self._scale()
        track = QRectF(x0, 4, m.length * sc, self.HEIGHT - 8)
        # same background as the bars around it; a thin outline marks the track
        p.setPen(QPen(pal.color(pal.ColorRole.Mid), 0.6))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(track)
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
                differs = diff is not None and \
                    float(diff[a:b].max(initial=0)) / v >= VARIABLE_SHARE
                p.fillRect(QRectF(x0 + k, track.bottom() - hgt, 1.0, hgt),
                           amber if differs else ok)
        # reference features as thin coloured bands along the foot of the strip
        for k, f in enumerate(m.features):
            if k < len(self.canvas._feat_row) and self.canvas._feat_row[k] < 0:
                continue
            fx0, fx1 = x0 + (f.start - 1) * sc, x0 + f.end * sc
            p.fillRect(QRectF(fx0, track.bottom() - 3, max(1.0, fx1 - fx0), 3),
                       self.canvas.feature_color(f))
        c = self.canvas
        if c.ppb > 0:
            vis0, vis1 = c.visible_span() or (0, m.length)
            # keep the 2 px outline inside the track at both ends
            left = max(x0 + 1.0, x0 + vis0 * sc)
            right = min(x0 + m.length * sc - 1.0, max(left + 4.0, x0 + vis1 * sc))
            box = QRectF(left, 2, right - left, self.HEIGHT - 4)
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


def _even_columns(view: QTableView, first: int | None = None) -> None:
    """Spread the columns evenly over the width (optionally a fixed first column),
    so tables line up with the view above instead of ending in a ragged gap."""
    header = view.horizontalHeader()
    header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
    if first is not None:
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
        view.setColumnWidth(0, first)


def _table(parent: QWidget) -> QTableView:
    view = QTableView(parent)
    view.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
    view.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
    view.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
    view.verticalHeader().setVisible(False)
    view.verticalHeader().setDefaultSectionSize(20)
    view.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
    # last column takes the rest of the width: no ragged gap on the right
    view.horizontalHeader().setStretchLastSection(True)
    view.setSortingEnabled(True)
    view.setAlternatingRowColors(True)
    return view


# --------------------------------------------------------------------------
# Panel
# --------------------------------------------------------------------------


_SETTINGS = "contig_view"


def _settings():
    from PySide6.QtCore import QSettings

    return QSettings()


def _read_bool(key: str, default: bool) -> bool:
    try:
        v = _settings().value(f"{_SETTINGS}/{key}", default)
    except Exception:  # noqa: BLE001 - settings are a convenience only
        return default
    if isinstance(v, str):
        return v.lower() in ("1", "true", "yes")
    return bool(v)


def _write(key: str, value) -> None:
    try:
        _settings().setValue(f"{_SETTINGS}/{key}", value)
    except Exception:  # noqa: BLE001
        pass


def _read_list(key: str) -> list[str]:
    try:
        v = _settings().value(f"{_SETTINGS}/{key}", [])
    except Exception:  # noqa: BLE001
        return []
    if v is None:
        return []
    if isinstance(v, str):
        return [v] if v else []
    return [str(x) for x in v]


def _read_text(key: str, default: str) -> str:
    try:
        v = _settings().value(f"{_SETTINGS}/{key}", default)
    except Exception:  # noqa: BLE001
        return default
    return str(v) if v not in (None, "") else default


def _read_json(key: str, default):
    import json

    raw = _read_text(key, "")
    if not raw:
        return default
    try:
        value = json.loads(raw)
    except ValueError:
        return default
    return value if isinstance(value, type(default)) else default


def _write_json(key: str, value) -> None:
    import json

    _write(key, json.dumps(value))


def _swatch_icon(color: str):
    """Rounded colour swatch for a tool button."""
    from PySide6.QtGui import QIcon, QPixmap

    icon = QIcon()
    for scale in (1, 2):
        pix = QPixmap(16 * scale, 12 * scale)
        pix.setDevicePixelRatio(scale)
        pix.fill(Qt.GlobalColor.transparent)
        pt = QPainter(pix)
        pt.setRenderHint(QPainter.RenderHint.Antialiasing)
        pt.setPen(QPen(QColor(0, 0, 0, 60), 1))
        pt.setBrush(QColor(color))
        pt.drawRoundedRect(QRectF(1.5, 1.5, 13, 9), 2, 2)
        pt.end()
        icon.addPixmap(pix)
    return icon


def _swatch(color: str):
    from PySide6.QtGui import QIcon, QPixmap

    pix = QPixmap(10, 10)
    pix.fill(QColor(color))
    return QIcon(pix)


def _sidebar_icon(widget: QWidget):
    """A window outline with its right column filled: the usual 'right panel' toggle."""
    from PySide6.QtGui import QIcon, QPixmap

    icon = QIcon()
    color = widget.palette().color(widget.palette().ColorRole.WindowText)
    for scale in (1, 2):
        pix = QPixmap(18 * scale, 14 * scale)
        pix.setDevicePixelRatio(scale)
        pix.fill(Qt.GlobalColor.transparent)
        p = QPainter(pix)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        frame = QRectF(1.5, 1.5, 15, 11)
        p.setPen(QPen(color, 1.2))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRoundedRect(frame, 2, 2)
        fill = QColor(color)
        fill.setAlphaF(0.55)
        p.fillRect(QRectF(11, 2.2, 5, 9.6), fill)
        p.drawLine(QPointF(11, 1.5), QPointF(11, 12.5))
        p.end()
        icon.addPixmap(pix)
    return icon


class ContigViewPanel(QWidget):
    """A Map to Reference result (``.rvmap``), laid out like Geneious: document tabs along
    the top (Contig view and its tables), a thin navigation bar, the view, and the display
    options in collapsible sections on the right."""

    openRequested = Signal(object)     # Path of a result file to open (report, consensus…)
    fileWritten = Signal(object)       # Path of a file written from the view
    blastRequested = Signal(object)    # Path of a FASTA to BLAST (right-click > BLAST)

    SIDE_W = 270

    def __init__(self, path: Path | None = None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        from PySide6.QtWidgets import QCheckBox, QLineEdit
        # one background colour for the title, the tab bar, the toolbar and the strips
        self.setAutoFillBackground(True)
        self.setBackgroundRole(self.palette().ColorRole.Window)

        from repeatviewer.widgets.combo_boxes import NoWheelComboBox
        from repeatviewer.widgets.spin_boxes import NoWheelDoubleSpinBox, NoWheelSpinBox

        self._path: Path | None = None
        self._data: dict = {}
        self._contig_models: dict[str, ContigModel] = {}
        self._report_loaded: Path | None = None
        self.canvas = ContigCanvas(self)
        self.overview = ContigOverview(self.canvas, self)
        # name column width: dragged by the user, remembered between sessions
        try:
            saved_w = int(_settings().value(f"{_SETTINGS}/label_width", LABEL_W))
        except (TypeError, ValueError):
            saved_w = LABEL_W
        self.canvas.label_w = max(LABEL_W_MIN, min(600, saved_w))
        self.canvas.labelWidthChanged.connect(lambda w: _write("label_width", int(w)))
        self._title = QLabel(self)
        self._title.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        # ------------------------------------------------ navigation bar (thin)
        self._contig = NoWheelComboBox(self)
        self._contig.setMinimumWidth(180)
        self._contig.currentIndexChanged.connect(self._on_contig_changed)

        def tool(text: str, tip: str, slot) -> QToolButton:
            b = QToolButton(self)
            b.setText(text)
            b.setToolTip(tip)
            b.setAutoRaise(True)
            b.clicked.connect(slot)
            return b

        alt = "Option" if sys.platform == "darwin" else "Alt"
        zoom_out = tool("−", "Zoom out\nMouse wheel or trackpad pinch zooms around the mouse; "
                        f"Shift + wheel scrolls the sequences, {alt} + wheel moves along",
                        lambda: self.canvas.zoom(0.66))
        zoom_in = tool("+", "Zoom in\nMouse wheel or trackpad pinch zooms around the mouse; "
                       f"Shift + wheel scrolls the sequences, {alt} + wheel moves along",
                       lambda: self.canvas.zoom(1.5))
        fit = tool("Fit", "Show the whole reference (or double-click empty space in the view)",
                   self.canvas.fit)
        bases = tool("Zoom to bases", "Zoom in until every base is written as a letter, "
                     "centred on the middle of the view",
                     lambda: self.canvas.show_base(self.canvas.bp_at(
                         self.canvas.label_w + self.canvas._data_w() / 2), zoom_to_letters=True))
        self._prev_diff = tool("◀", "Previous difference", lambda: self._step_difference(-1))
        self._next_diff = tool("▶", "Next difference", lambda: self._step_difference(1))
        self._diff_label = QLabel(self)
        self._diff_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._diff_label.setMinimumWidth(140)
        self._goto = QLineEdit(self)
        self._goto.setPlaceholderText("Go to, e.g. 90-130")
        self._goto.setMinimumWidth(120)
        self._goto.setClearButtonEnabled(True)
        self._goto.setToolTip("Go to a position or a range on the reference "
                              "(e.g. 103, 90-130, or CroSat2:90-130). Enter to go.")
        self._goto.returnPressed.connect(self._go_to)
        # Files ▾: the arrow is part of the label (Qt's own menu arrow sits low on macOS),
        # and the menu opens under the button.
        self._files_btn = QToolButton(self)
        self._files_btn.setText("Files ▾")
        self._files_btn.setAutoRaise(True)
        self._files_btn.setToolTip("Files written by this run: consensus, placed and "
                                   "unused sequences, aligner output, reference")
        self._files_menu = QMenu(self._files_btn)
        self._files_btn.clicked.connect(
            lambda: self._files_menu.exec(
                self._files_btn.mapToGlobal(self._files_btn.rect().bottomLeft())))
        # panel toggle: sidebar icon, highlighted while the options panel is open
        self._side_btn = QToolButton(self)
        self._side_btn.setCheckable(True)
        self._side_btn.setAutoRaise(True)
        self._side_btn.setIcon(_sidebar_icon(self))
        self._side_btn.setToolTip("Show or hide the options panel")
        self._side_btn.toggled.connect(self._set_side_visible)

        bar = QWidget(self)
        bar.setObjectName("contigNavBar")
        bar.setAutoFillBackground(True)
        bar.setBackgroundRole(bar.palette().ColorRole.Window)
        bl = QHBoxLayout(bar)
        bl.setContentsMargins(6, 3, 6, 3)
        bl.setSpacing(4)
        bl.addWidget(self._contig)
        for wdg in (zoom_out, zoom_in, fit, bases):
            bl.addWidget(wdg)
        bl.addStretch(1)
        for wdg in (self._prev_diff, self._diff_label, self._next_diff):
            bl.addWidget(wdg)
        bl.addStretch(1)
        # Go to · Files · panel toggle: over the options panel, Go to starting where
        # its section titles start (see _sync_nav_right)
        self._nav_right = QWidget(bar)
        nr = QHBoxLayout(self._nav_right)
        nr.setContentsMargins(0, 0, 0, 0)
        nr.setSpacing(4)
        nr.addWidget(self._goto, 1)
        nr.addSpacing(2)
        nr.addWidget(self._files_btn)
        nr.addWidget(self._side_btn)
        # legend: the app's shared round "?"; click shows the key as a popup
        from repeatviewer.widgets.help_button import make_help_popup_button

        self._legend_btn = make_help_popup_button(
            self, "Legend", "What the colours and marks mean", self._show_legend)
        nr.addSpacing(2)
        nr.addWidget(self._legend_btn)
        bl.addWidget(self._nav_right)
        self._nav_bar = bar

        # ------------------------------------------------ side panel controls
        self._layout_mode = NoWheelComboBox(self)
        self._layout_mode.addItem("One per row", True)
        self._layout_mode.addItem("Packed", False)
        self._layout_mode.currentIndexChanged.connect(self._on_layout_mode)
        self._sort = NoWheelComboBox(self)
        for text, k in (("Start", SORT_START), ("Strand", SORT_STRAND),
                        ("Identity", SORT_IDENTITY), ("Name", SORT_NAME)):
            self._sort.addItem(text, k)
        self._sort.setToolTip("Order of the sequences when shown one per row")
        self._sort.currentIndexChanged.connect(
            lambda _i: self.canvas.set_sort(str(self._sort.currentData())))
        self._color_by = NoWheelComboBox(self)
        for text, k in (("Strand", COLOR_STRAND), ("Identity", COLOR_IDENTITY),
                        ("Mapping quality", COLOR_MAPQ)):
            self._color_by.addItem(text, k)
        self._color_by.setToolTip("Colour of the sequence bars when zoomed out (the scale is in "
                                  "the ? legend): forward blue / "
                                  "reverse purple, identity red (80 %) to green (100 %), or "
                                  "mapping quality amber (0) to blue (60)")
        self._color_by.currentIndexChanged.connect(
            lambda _i: self.canvas.set_color_by(str(self._color_by.currentData())))

        self._threshold = NoWheelComboBox(self)
        for value, label in CONSENSUS_THRESHOLDS:
            self._threshold.addItem(label, value)
        self._threshold.setToolTip("Share of the sequences a base needs to be called; "
                                   "below it the consensus writes an ambiguity code")
        self._threshold.currentIndexChanged.connect(
            lambda _i: self._set_threshold(self._threshold.currentData()))
        self._threshold_value: float | None = None
        self._no_cov_value = "N"

        self._min_identity = NoWheelDoubleSpinBox(self)
        self._min_identity.setRange(0.0, 100.0)
        self._min_identity.setDecimals(0)
        self._min_identity.setSingleStep(5.0)
        self._min_identity.setSuffix(" %")
        self._min_identity.setToolTip("Hide sequences placed with a lower identity")
        self._min_mapq = NoWheelSpinBox(self)
        self._min_mapq.setRange(0, 60)
        self._min_mapq.setPrefix("≥ ")
        self._min_mapq.setToolTip("Hide placements with a lower mapping quality (MAPQ, "
                                  "0–60: how sure the mapper is that this is the only "
                                  "good place). Placements without one are kept.")
        self._hide_secondary = QCheckBox("Hide secondary placements", self)
        self._hide_secondary.setToolTip("Show only each sequence's primary placement")
        self._hidden_label = QLabel(self)
        self._hidden_label.setEnabled(False)
        for w_ in (self._min_identity, self._min_mapq):
            w_.valueChanged.connect(lambda _v: self._apply_filters())
        self._hide_secondary.toggled.connect(lambda _on: self._apply_filters())

        self._show_ann = QCheckBox("Show annotations", self)
        self._show_ann.setChecked(_read_bool("annotations_on", True))
        self._show_ann.toggled.connect(self._on_show_annotations)
        self._type_boxes: dict[str, QCheckBox] = {}
        self._types_host = QWidget(self)
        self._types_lay = QVBoxLayout(self._types_host)
        self._types_lay.setContentsMargins(18, 0, 0, 0)
        self._types_lay.setSpacing(2)
        self._per_seq = QCheckBox("In each sequence", self)
        self._per_seq.setToolTip("Draw the reference features under every sequence: intact, "
                                 "changed (\"n diff\"), partial (faded) or lost (grey)")
        self._per_seq.setChecked(_read_bool("read_features", False))
        self._per_seq.toggled.connect(self._on_read_features)
        self._no_ann = QLabel("The reference has no annotations", self)
        self._no_ann.setEnabled(False)
        self._hidden_types: set[str] = set(_read_list("hidden_types"))

        self._g_log = QCheckBox("Coverage on a log scale", self)
        self._g_log.setChecked(True)
        self._g_logo = QCheckBox("Sequence logo", self)
        self._g_logo.setToolTip("Drawn when the bases are written out")

        def threshold_row(text: str, value: int, tip: str):
            host = QWidget(self)
            row = QHBoxLayout(host)
            row.setContentsMargins(0, 0, 0, 0)
            row.setSpacing(6)
            box = QCheckBox(text, host)
            box.setToolTip(tip)
            spin = NoWheelSpinBox(host)
            spin.setRange(0, 1_000_000)
            spin.setValue(value)
            spin.setEnabled(False)
            spin.setFixedWidth(84)
            box.toggled.connect(spin.setEnabled)
            row.addWidget(box, 1)
            row.addWidget(spin)
            return box, spin, host

        self._g_low_box, self._g_low, low_row = threshold_row(
            "Coverage below", 5, "Mark positions covered by fewer sequences (red)")
        self._g_high_box, self._g_high, high_row = threshold_row(
            "Coverage above", 100, "Mark positions covered by more sequences (orange)")
        for box in (self._g_log, self._g_logo, self._g_low_box, self._g_high_box):
            box.toggled.connect(lambda _on: self._apply_graphs())
        for spin in (self._g_low, self._g_high):
            spin.valueChanged.connect(lambda _v: self._apply_graphs())

        # Highlighting (Geneious): what to colour, compared to what
        self._hl_show = NoWheelComboBox(self)
        for text, key in (("Disagreements", ContigCanvas.SHOW_DISAGREE),
                          ("Agreements", ContigCanvas.SHOW_AGREE),
                          ("All bases", ContigCanvas.SHOW_ALL),
                          ("Nothing", ContigCanvas.SHOW_NONE)):
            self._hl_show.addItem(text, key)
        self._hl_show.setToolTip("Bases to colour: those that differ, those that agree, "
                                 "every base, or none (plain letters)")
        self._hl_target = NoWheelComboBox(self)
        self._hl_target.addItem("Reference", ContigCanvas.TO_REFERENCE)
        self._hl_target.addItem("Consensus", ContigCanvas.TO_CONSENSUS)
        self._hl_target.setToolTip("Compare each sequence with the reference, or with the "
                                   "consensus of all sequences (shows which copies deviate "
                                   "from the family). ◀ ▶ in the toolbar steps through "
                                   "these differences.")
        self._hl_dots = QCheckBox("Dots where same", self)
        self._hl_dots.setChecked(True)
        self._hl_dots.setToolTip("Write a dot instead of the letter where a base agrees")
        self._hl_show.setCurrentIndex(max(0, self._hl_show.findData(
            _read_text("highlight_show", ContigCanvas.SHOW_DISAGREE))))
        self._hl_target.setCurrentIndex(max(0, self._hl_target.findData(
            _read_text("highlight_target", ContigCanvas.TO_REFERENCE))))
        self._hl_dots.setChecked(_read_bool("highlight_dots", True))
        for combo in (self._hl_show, self._hl_target):
            combo.currentIndexChanged.connect(lambda _i: self._on_highlight())
        self._hl_dots.toggled.connect(lambda _on: self._on_highlight())

        # Wrap contig (Display)
        self._wrap = QCheckBox("Wrap contig", self)
        self._wrap.setToolTip("When the bases are written out, continue the alignment on "
                              "the next block instead of scrolling sideways")
        self._wrap.setChecked(_read_bool("wrap", False))
        self._wrap.toggled.connect(self._on_wrap)

        # Annotation labels
        self._labels = NoWheelComboBox(self)
        for text, key in (("Inside", ContigCanvas.LABELS_INSIDE),
                          ("Outside", ContigCanvas.LABELS_OUTSIDE),
                          ("None", ContigCanvas.LABELS_NONE)):
            self._labels.addItem(text, key)
        self._labels.setToolTip("Names inside the bars, above them, or not at all. "
                                "Names that do not fit are left out.")
        self._labels.setCurrentIndex(max(0, self._labels.findData(
            _read_text("labels", ContigCanvas.LABELS_INSIDE))))
        self._labels.currentIndexChanged.connect(self._on_labels)
        self._type_colors: dict[str, str] = _read_json("type_colors", {})

        self._side = self._build_side(low_row, high_row)

        # ------------------------------------------------ view page
        view_area = QWidget(self)
        va = QVBoxLayout(view_area)
        va.setContentsMargins(0, 0, 0, 0)
        va.setSpacing(0)
        va.addWidget(self.overview)
        va.addWidget(self.canvas, 1)
        self._view_split = QSplitter(Qt.Orientation.Horizontal, self)
        self._view_split.addWidget(view_area)
        self._view_split.addWidget(self._side)
        self._view_split.setStretchFactor(0, 1)
        self._view_split.setStretchFactor(1, 0)
        self._view_split.setCollapsible(0, False)
        self._view_split.setChildrenCollapsible(False)
        self._view_split.splitterMoved.connect(self._queue_nav_sync)
        self._view_page = QWidget(self)
        self._view_page.setAutoFillBackground(True)
        self._view_page.setBackgroundRole(self._view_page.palette().ColorRole.Window)
        vp = QVBoxLayout(self._view_page)
        vp.setContentsMargins(0, 0, 0, 0)
        vp.setSpacing(0)
        vp.addWidget(bar)
        vp.addWidget(self._view_split, 1)

        # ------------------------------------------------ document tabs
        from repeatviewer.widgets.html_preview import create_html_preview_widget

        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        self._reads_view = _table(self)
        self._features_view = _table(self)
        self._features_page = self._build_features_page()
        self._variants_view = _table(self)
        self._unused_view = _table(self)
        self._stats_view = _table(self)
        self._stats_view.setSortingEnabled(False)   # fixed order, like a report
        self._report_view = create_html_preview_widget(self)
        for page, text in ((self._view_page, "Contig view"), (self._reads_view, "Sequences"),
                           (self._features_page, "Features"), (self._variants_view, "Variants"),
                           (self._unused_view, "Not used"), (self._stats_view, "Statistics"),
                           (self._report_view, "Report")):
            self._tabs.addTab(page, text)
        self._tabs.setTabToolTip(1, "Click a row to show the sequence in the Contig view")
        self._tabs.setTabToolTip(2, "Click a cell to show the feature in that sequence")
        self._tabs.setTabToolTip(3, "Click a row to show the position in the Contig view")
        self._tabs.currentChanged.connect(self._on_tab_changed)
        self._reads_view.clicked.connect(self._on_read_row)
        self._features_view.clicked.connect(self._on_feature_cell)
        self._variants_view.clicked.connect(self._on_variant_row)
        self.canvas.readClicked.connect(self._select_read_row)
        self.canvas.readDoubleClicked.connect(self._open_read_in_rows)
        self.canvas.readContextMenu.connect(self._read_menu)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(self._title)
        lay.addWidget(self._tabs, 1)
        self._title.setContentsMargins(8, 6, 8, 4)
        self._side_btn.setChecked(_read_bool("side_panel", True))
        self._set_side_visible(self._side_btn.isChecked())
        self.canvas.type_colors = dict(self._type_colors)
        self.canvas.set_label_mode(str(self._labels.currentData()))
        self.canvas.set_wrap(self._wrap.isChecked())
        self._on_highlight()
        self._apply_annotations()
        if path is not None:
            self.load(path)

    # ---------------------------------------------------------------- side panel

    def _build_side(self, low_row: QWidget, high_row: QWidget) -> QWidget:
        from PySide6.QtWidgets import QFormLayout, QFrame, QScrollArea

        from repeatviewer.widgets.sequence_annotation_rail import _CollapsibleSection

        body = QWidget(self)
        col = QVBoxLayout(body)
        col.setContentsMargins(10, 8, 10, 8)
        col.setSpacing(12)
        self._sections: dict[str, QWidget] = {}

        def section(key: str, title: str, default: bool = True):
            s = _CollapsibleSection(title, body, expanded=_read_bool(f"section/{key}", default))
            s.expandedChanged.connect(lambda on, k=key: _write(f"section/{k}", bool(on)))
            self._sections[key] = s
            col.addWidget(s)
            return s

        def form(*rows) -> QFormLayout:
            f = QFormLayout()
            f.setContentsMargins(0, 0, 0, 0)
            f.setHorizontalSpacing(8)
            f.setVerticalSpacing(6)
            f.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
            for label, w in rows:
                lab = QLabel(label, body)
                labels.append(lab)
                f.addRow(lab, w)
            return f

        labels: list[QLabel] = []

        s = section("display", "Display")
        s.add_layout(form(("Layout", self._layout_mode), ("Sort by", self._sort),
                          ("Colour by", self._color_by)))
        s.add_widget(self._wrap)
        s = section("highlighting", "Highlighting")
        s.add_layout(form(("Show", self._hl_show), ("Compared to", self._hl_target)))
        s.add_widget(self._hl_dots)
        s = section("consensus", "Consensus")
        s.add_layout(form(("Threshold", self._threshold)))
        s = section("filters", "Filters")
        self._filters_form = form(("Min identity", self._min_identity),
                                  ("Mapping quality", self._min_mapq))
        s.add_layout(self._filters_form)
        from repeatviewer.widgets.help_button import field_with_help

        self._hide_secondary_row = field_with_help(
            self._hide_secondary, "Secondary placements",
            "A sequence can be placed in several spots when the Mapping run keeps more than "
            "one place (Mapping window, Minimap2 tab, <b>Secondaries: Up to 50 places</b>). "
            "Its best spot is the <b>primary</b> placement; the others are "
            "<b>secondary</b> and are shown with #2, #3… after the name.<br><br>"
            "Tick to keep only each sequence's primary placement.", parent=body)
        s.add_widget(self._hide_secondary_row)
        s.add_widget(self._hidden_label)
        s = section("annotations", "Annotations")
        s.add_widget(self._show_ann)
        self._labels_form = form(("Labels", self._labels))
        s.add_layout(self._labels_form)
        s.add_widget(self._types_host)
        s.add_widget(self._per_seq)
        s.add_widget(self._no_ann)
        s = section("graphs", "Graphs", default=False)
        for w in (self._g_log, low_row, high_row, self._g_logo):
            s.add_widget(w)
        col.addStretch(1)
        # column width from the everyday labels; the rare long one wraps onto two lines
        long_labels = [lab for lab in labels if lab.text() == "Mapping quality"]
        label_w = max(lab.sizeHint().width() for lab in labels if lab not in long_labels)
        for lab in long_labels:
            lab.setWordWrap(True)
        for combo in (self._layout_mode, self._sort, self._color_by, self._hl_show,
                      self._hl_target, self._threshold, self._labels):
            combo.setSizeAdjustPolicy(
                QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            combo.setMinimumContentsLength(6)
        for spin in (self._min_identity, self._min_mapq):
            spin.setMinimumWidth(60)
        for lab in labels:                       # fields line up across sections
            lab.setFixedWidth(label_w)

        scroll = QScrollArea(self)
        scroll.setWidget(body)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setMinimumWidth(self.SIDE_W - 30)
        scroll.setMaximumWidth(self.SIDE_W + 120)
        return scroll

    SIDE_PAD = 10       # left margin of the options panel's sections

    def _show_legend(self) -> None:
        """Popup under the "?" with only the marks this result can show. A QMenu, so it
        has the app's menu look (rounded, theme.py) like every other popup."""
        from PySide6.QtWidgets import QWidgetAction

        c = self.canvas
        items = legend_items(c.model, str(c.highlight), target=c.hl_target, show=c.hl_show)
        if items:
            items = items + bar_colour_items(c.color_by)
        pop = QMenu(self)
        holder = QWidgetAction(pop)
        holder.setDefaultWidget(ContigLegend(items, pop, strand_bars=c.color_by == COLOR_STRAND)
                                if items else
                                QLabel("Open a result to see its legend", pop))
        pop.addAction(holder)
        btn = self._legend_btn
        size = pop.sizeHint()
        self._legend_popup = pop
        pop.popup(btn.mapToGlobal(QPoint(btn.width() - size.width(), btn.height() + 4)))

    def _sync_nav_right(self) -> None:
        """Line Go to up with the options panel: its left edge at the section titles,
        the group ending at the bar's right margin. Without the panel: natural size."""
        right = self._nav_right
        if not self._side.isVisible() or self._side.width() <= 0:
            right.setMinimumWidth(0)
            right.setMaximumWidth(16777215)
            self._goto.setFixedWidth(190)
            return
        self._goto.setMinimumWidth(120)
        self._goto.setMaximumWidth(16777215)
        bar = self._nav_bar
        start = bar.mapFromGlobal(self._side.mapToGlobal(QPoint(0, 0))).x() + self.SIDE_PAD
        end = bar.width() - bar.layout().contentsMargins().right()
        if end - start > 160:
            right.setFixedWidth(end - start)

    def _queue_nav_sync(self, *_args) -> None:
        from PySide6.QtCore import QTimer

        QTimer.singleShot(0, self._sync_nav_right)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._queue_nav_sync()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._queue_nav_sync()

    def _set_side_visible(self, on: bool) -> None:
        self._side.setVisible(on)
        _write("side_panel", bool(on))
        if on:
            sizes = self._view_split.sizes()
            total = sum(sizes) or self.width()
            if len(sizes) == 2 and sizes[1] < self.SIDE_W - 30:
                self._view_split.setSizes([max(0, total - self.SIDE_W), self.SIDE_W])
        self._queue_nav_sync()

    # ---------------------------------------------------------------- loading

    def path(self) -> Path | None:
        return self._path

    def load(self, path: Path) -> None:
        self._path = Path(path)
        self._report_loaded = None
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
        # Only the title here; the numbers are in the Report and Statistics tabs.
        self._title.setText(f"<b>{d.get('title', self._path.name)}</b>")
        summary = (f"{s.get('n_placed', 0):,} of {s.get('n_sequences', 0):,} placed "
                   f"({s.get('n_reversed', 0):,} reversed, {s.get('n_unused', 0):,} not used)"
                   + (f" · mean identity {mean:.1%}" if mean is not None else "")
                   + f" · {engine} {d.get('engine_version', '')}")
        self._title.setToolTip("\n".join([summary.strip(), *d.get("warnings", [])]))
        self._contig_models.clear()
        self._contig.blockSignals(True)
        self._contig.clear()
        for c in d.get("contigs", []):
            self._contig.addItem(f"{c['id']} ({c['n_reads']:,})", c["id"])
        self._contig.blockSignals(False)
        self._contig.setVisible(self._contig.count() > 1)
        self._build_files_menu()
        self._fill_unused()
        report = self._report_path()
        self._tabs.setTabEnabled(self._tabs.indexOf(self._report_view), report is not None)
        if self._tabs.currentWidget() is self._report_view:
            self._tabs.setCurrentWidget(self._view_page)
        if self._contig.count():
            self._on_contig_changed(0)
        else:
            self.canvas.set_model(None)

    def _report_path(self) -> Path | None:
        if self._path is None:
            return None
        name = self._data.get("details", {}).get("report")
        if not name:
            return None
        target = self._path.parent / self._data.get("details_dir", "") / name
        return target if target.is_file() else None

    def _on_tab_changed(self, _index: int) -> None:
        if self._tabs.currentWidget() is not self._report_view:
            return
        report = self._report_path()
        if report is not None and report != self._report_loaded:
            self._report_loaded = report
            try:
                self._report_view.load_html_file(report)
            except Exception:  # noqa: BLE001 - the report stays in the details folder
                self._report_loaded = None

    def _show_view(self) -> None:
        if self._tabs.currentWidget() is not self._view_page:
            self._tabs.setCurrentWidget(self._view_page)
            QApplication.processEvents()

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
        self._sort.setEnabled(self.canvas.one_per_row)
        self._sync_filter_rows(model)
        self._fill_reads(model)
        self._fill_features(model)
        self._fill_variants(model)
        self._build_annotation_types()
        self._fill_statistics()
        self._sync_difference_label()

    def _on_layout_mode(self, _index: int = 0) -> None:
        one = bool(self._layout_mode.currentData())
        self.canvas.set_one_per_row(one)
        self._sort.setEnabled(one)
        self._sort.setToolTip("Order of the sequences when shown one per row" if one else
                              "Sorting applies to One per row; Packed places sequences "
                              "by position")

    def _sync_filter_rows(self, m: "ContigModel") -> None:
        """Show Min mapping quality only when the run's placements differ in mapping
        quality, and Hide secondary only when there are secondary placements."""
        qualities = {int(r.data.get("mapq", -1)) for r in m.reads
                     if int(r.data.get("mapq", -1)) >= 0}
        varied = len(qualities) > 1
        self._filters_form.setRowVisible(self._min_mapq, varied)
        if not varied and self._min_mapq.value():
            self._min_mapq.setValue(0)
        secondary = any(str(r.data.get("kind", "")) != "primary"
                        or int(r.data.get("placement", 1) or 1) > 1 for r in m.reads)
        self._hide_secondary_row.setVisible(secondary)
        if not secondary and self._hide_secondary.isChecked():
            self._hide_secondary.setChecked(False)

    def _on_highlight(self, *_args) -> None:
        show = str(self._hl_show.currentData())
        target = str(self._hl_target.currentData())
        self._hl_dots.setEnabled(show == ContigCanvas.SHOW_DISAGREE)
        self.canvas.set_highlighting(show, target, self._hl_dots.isChecked())
        _write("highlight_show", show)
        _write("highlight_target", target)
        _write("highlight_dots", self._hl_dots.isChecked())
        self._sync_difference_label()

    def _on_wrap(self, on: bool) -> None:
        _write("wrap", bool(on))
        self.canvas.set_wrap(bool(on))

    def _on_labels(self, _index: int = 0) -> None:
        mode = str(self._labels.currentData())
        _write("labels", mode)
        self.canvas.set_label_mode(mode)

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
                   "Inserted", "Mapping quality", "Kind", "#"]
        self._reads_view.setModel(_RowsModel(headers, rows, parent=self))
        self._reads_view.setColumnHidden(len(headers) - 1, True)
        _even_columns(self._reads_view, first=220)

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
        head_fm = self._features_view.horizontalHeader().fontMetrics()
        for c in range(1, len(headers) - 1):
            # wide enough for the feature name in the header (e.g. Subrepeat_1)
            self._features_view.setColumnWidth(
                c, max(70, head_fm.horizontalAdvance(headers[c]) + 22))
        self._tabs.setTabEnabled(self._tabs.indexOf(self._features_page), bool(feats))
        self._fill_annotation_list(m)

    def _build_features_page(self) -> QWidget:
        """Features tab: every annotation of the reference (filterable, click to show it)
        above the per-sequence table of how each feature came out."""
        from PySide6.QtCore import QSortFilterProxyModel
        from PySide6.QtWidgets import QLineEdit

        page = QWidget(self)
        lay = QVBoxLayout(page)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        bar = QWidget(page)
        bl = QHBoxLayout(bar)
        bl.setContentsMargins(6, 4, 6, 4)
        self._ann_filter = QLineEdit(bar)
        self._ann_filter.setPlaceholderText("Filter annotations by name or type")
        self._ann_filter.setClearButtonEnabled(True)
        self._ann_count = QLabel(bar)
        self._ann_count.setEnabled(False)
        bl.addWidget(self._ann_filter, 1)
        bl.addWidget(self._ann_count)
        self._ann_view = _table(page)
        self._ann_proxy = QSortFilterProxyModel(self)
        self._ann_proxy.setFilterCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
        self._ann_proxy.setFilterKeyColumn(-1)
        self._ann_proxy.setSortRole(Qt.ItemDataRole.UserRole)
        self._ann_view.setModel(self._ann_proxy)
        self._ann_filter.textChanged.connect(self._on_ann_filter)
        self._ann_view.clicked.connect(self._on_annotation_row)
        split = QSplitter(Qt.Orientation.Vertical, page)
        top = QWidget(split)
        tl = QVBoxLayout(top)
        tl.setContentsMargins(0, 0, 0, 0)
        tl.setSpacing(0)
        tl.addWidget(bar)
        tl.addWidget(self._ann_view, 1)
        split.addWidget(top)
        split.addWidget(self._features_view)
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 2)
        lay.addWidget(split, 1)
        return page

    def _fill_annotation_list(self, m: ContigModel) -> None:
        rows = []
        for k, f in enumerate(m.features):
            summary = m.feature_summary[k] if k < len(m.feature_summary) else {}
            rows.append([f.name, str(f.feature_type or "other"), f.start, f.end,
                         f.end - f.start + 1, "+" if f.is_forward() else "−",
                         feature_status_text(summary) if summary else "", k])
        headers = ["Name", "Type", "Start", "End", "Length", "Strand", "In the sequences", "#"]
        self._ann_rows = _RowsModel(headers, rows, parent=self)
        self._ann_proxy.setSourceModel(self._ann_rows)
        self._ann_view.setColumnHidden(len(headers) - 1, True)
        _even_columns(self._ann_view, first=200)
        self._on_ann_filter(self._ann_filter.text())

    def _on_ann_filter(self, text: str) -> None:
        self._ann_proxy.setFilterFixedString(text.strip())
        total = self._ann_rows.rowCount() if hasattr(self, "_ann_rows") else 0
        shown = self._ann_proxy.rowCount()
        self._ann_count.setText(f"{shown:,} of {total:,}" if shown != total else f"{total:,}")

    def _on_annotation_row(self, index) -> None:
        m = self.canvas.model
        if m is None:
            return
        src = self._ann_proxy.mapToSource(index)
        k = self._ann_rows.index(src.row(), self._ann_rows.columnCount() - 1).data(
            Qt.ItemDataRole.UserRole)
        if k is None or not 0 <= int(k) < len(m.features):
            return
        f = m.features[int(k)]
        t = str(f.feature_type or "other")
        if t in self._hidden_types or not self._show_ann.isChecked():
            self._show_ann.setChecked(True)
            self._toggle_type(t, True)
        self._show_view()
        self.canvas.show_range(f.start - 1, f.end)

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
        _even_columns(self._variants_view)

    def _fill_unused(self) -> None:
        rows = [[n, r] for n, r in self._data.get("unused", [])]
        self._unused_view.setModel(_RowsModel(["Sequence", "Reason"], rows, parent=self))
        _even_columns(self._unused_view, first=240)

    # ---------------------------------------------------------------- interaction

    def _read_index_from(self, view: QTableView, index) -> int:
        model = view.model()
        return int(model.index(index.row(), model.columnCount() - 1).data(
            Qt.ItemDataRole.UserRole))

    def _on_read_row(self, index) -> None:
        read_idx = self._read_index_from(self._reads_view, index)
        self._show_view()
        self.canvas.selected = read_idx
        self.canvas.show_read(read_idx)

    def _on_feature_cell(self, index) -> None:
        m = self.canvas.model
        if m is None:
            return
        read_idx = self._read_index_from(self._features_view, index)
        self._show_view()
        col = index.column()
        self.canvas.selected = read_idx
        self.canvas.show_read(read_idx)
        if 1 <= col <= min(len(m.features), mr.MAX_FEATURE_COLUMNS):
            f = m.features[col - 1]
            self.canvas.show_base((f.start + f.end) // 2)

    def _on_variant_row(self, index) -> None:
        pos = self._variants_view.model().index(index.row(), 0).data(Qt.ItemDataRole.UserRole)
        if pos:
            self._show_view()
            self.canvas.show_base(int(pos) - 1, zoom_to_letters=True, mark=True)
            self._sync_difference_label()

    def _apply_filters(self) -> None:
        self.canvas.set_filters(min_identity=self._min_identity.value() / 100.0,
                                min_mapq=self._min_mapq.value(),
                                hide_secondary=self._hide_secondary.isChecked())
        n = self.canvas.hidden_count()
        self._hidden_label.setText(f"{n:,} hidden by the filters" if n else "")
        self._hidden_label.setVisible(bool(n))
        self.overview.update()
        if hasattr(self, "_stats_view"):
            self._fill_statistics()

    # ---------------------------------------------------------------- go to / reads

    def _go_to(self) -> None:
        from repeatviewer.backends.motif_search import parse_position

        text = self._goto.text().strip()
        parsed = parse_position(text)
        if parsed is None:
            self._goto.setToolTip(f"\"{text}\" is not a position. Type e.g. 103 or 90-130.")
            self._goto.selectAll()
            return
        name, start, end = parsed
        if name:
            for i in range(self._contig.count()):
                if str(self._contig.itemData(i)) == name:
                    self._contig.setCurrentIndex(i)
                    break
        m = self.canvas.model
        if m is None:
            return
        if start > m.length:
            self._goto.setToolTip(f"The reference is {m.length:,} bp long.")
            self._goto.selectAll()
            return
        self._show_view()
        self.canvas.show_range(start - 1, min(end, m.length))
        self._sync_difference_label()

    def _open_read_in_rows(self, read_idx: int) -> None:
        """Double-click a read (packed view): show it in One sequence per row."""
        if not self.canvas.one_per_row:
            self._layout_mode.setCurrentIndex(0)       # One sequence per row
        self.canvas.selected = read_idx
        self.canvas.show_read(read_idx)
        self._select_read_row(read_idx)

    # ---------------------------------------------------------------- annotations

    def _all_types(self) -> list[tuple[str, int, str]]:
        m = self.canvas.model
        if m is None:
            return []
        return [(t, n, self._type_colors.get(t, c)) for t, n, c in feature_type_counts(m.features)]

    def _build_annotation_types(self) -> None:
        """One row per feature type: tick (show / hide), name and count, a colour
        swatch (click to change, remembered), and ◀ ▶ to jump to the previous / next
        annotation of that type."""
        from PySide6.QtWidgets import QCheckBox

        while self._types_lay.count():
            item = self._types_lay.takeAt(0)
            if item.widget() is not None:
                item.widget().deleteLater()
        self._type_boxes.clear()
        self._type_swatches: dict[str, QToolButton] = {}
        types = self._all_types()
        for name, n, color in types:
            row = QWidget(self._types_host)
            rl = QHBoxLayout(row)
            rl.setContentsMargins(0, 0, 0, 0)
            rl.setSpacing(2)
            box = QCheckBox(f"{name}  ({n:,})", row)
            box.setChecked(name not in self._hidden_types)
            box.toggled.connect(lambda on, t=name: self._toggle_type(t, on))
            self._type_boxes[name] = box
            swatch = QToolButton(row)
            swatch.setAutoRaise(True)
            swatch.setIcon(_swatch_icon(self._type_colors.get(name, color)))
            swatch.setToolTip(f"Colour of {name} annotations (click to change)")
            swatch.clicked.connect(lambda _c=False, t=name: self._pick_type_color(t))
            self._type_swatches[name] = swatch
            rl.addWidget(box, 1)
            rl.addWidget(swatch)
            for text, step in (("◀", -1), ("▶", 1)):
                jump = QToolButton(row)
                jump.setText(text)
                jump.setAutoRaise(True)
                jump.setToolTip(("Previous" if step < 0 else "Next") + f" {name}")
                jump.clicked.connect(lambda _c=False, t=name, d=step: self._jump_type(t, d))
                rl.addWidget(jump)
            self._types_lay.addWidget(row)
        has = bool(types)
        for w in (self._show_ann, self._types_host, self._per_seq):
            w.setVisible(has)
        self._labels_form.setRowVisible(self._labels, has)
        self._no_ann.setVisible(not has)
        self._apply_annotations()

    def _pick_type_color(self, name: str) -> None:
        from PySide6.QtWidgets import QColorDialog

        current = next((c for t, _n, c in self._all_types() if t == name), "#888888")
        chosen = QColorDialog.getColor(QColor(current), self, f"Colour for {name}")
        if not chosen.isValid():
            return
        self._type_colors[name] = chosen.name()
        _write_json("type_colors", self._type_colors)
        swatch = self._type_swatches.get(name)
        if swatch is not None:
            swatch.setIcon(_swatch_icon(chosen.name()))
        self.canvas.set_type_colors(self._type_colors)
        self.overview.update()

    def _jump_type(self, name: str, direction: int) -> None:
        """Show the previous / next annotation of one type, counted from the middle of
        the view, and switch the type on if it was hidden."""
        m = self.canvas.model
        if m is None:
            return
        feats = sorted((f for f in m.features if str(f.feature_type or "other") == name),
                       key=lambda f: (f.start, f.end))
        if not feats:
            return
        if name in self._hidden_types or not self._show_ann.isChecked():
            self._show_ann.setChecked(True)
            self._toggle_type(name, True)
        here = self.canvas.view_center_bp() or 0
        if direction > 0:
            nxt = next((f for f in feats if f.start - 1 > here), feats[0])
        else:
            nxt = next((f for f in reversed(feats) if f.end - 1 < here), feats[-1])
        self._show_view()
        self.canvas.show_range(nxt.start - 1, nxt.end)

    def _apply_annotations(self) -> None:
        on = self._show_ann.isChecked()
        names = {t for t, _n, _c in self._all_types()}
        hidden = set(names) if not on else (self._hidden_types & names)
        self._types_host.setEnabled(on)
        self._per_seq.setEnabled(on)
        self.canvas.set_hidden_types(hidden)
        self.canvas.set_show_read_features(on and self._per_seq.isChecked())
        self.overview.update()

    def _on_show_annotations(self, on: bool) -> None:
        _write("annotations_on", bool(on))
        self._apply_annotations()

    def _toggle_type(self, name: str, on: bool) -> None:
        (self._hidden_types.discard if on else self._hidden_types.add)(name)
        box = self._type_boxes.get(name)
        if box is not None and box.isChecked() != on:
            box.blockSignals(True)
            box.setChecked(on)
            box.blockSignals(False)
        _write("hidden_types", sorted(self._hidden_types))
        self._apply_annotations()

    def _set_all_types(self, on: bool) -> None:
        for name, _n, _c in self._all_types():
            self._toggle_type(name, on)

    def _on_read_features(self, on: bool) -> None:
        _write("read_features", bool(on))
        self._apply_annotations()

    # ---------------------------------------------------------------- graphs / consensus

    def _apply_graphs(self) -> None:
        self.canvas.set_graph_options(
            log_scale=self._g_log.isChecked(),
            low=self._g_low.value() if self._g_low_box.isChecked() else None,
            high=self._g_high.value() if self._g_high_box.isChecked() else None,
            single_strand=False,
            logo=self._g_logo.isChecked())

    def _set_threshold(self, value: float | None) -> None:
        self._threshold_value = value
        i = self._threshold.findData(value)
        if value is None:
            i = 0
        if i >= 0 and i != self._threshold.currentIndex():
            self._threshold.blockSignals(True)
            self._threshold.setCurrentIndex(i)
            self._threshold.blockSignals(False)
        self.canvas.set_consensus(value, self._no_cov_value)
        self._fill_statistics()

    def _fill_statistics(self) -> None:
        m = self.canvas.model
        rows = [[k, v] for k, v in contig_statistics(m)] if m is not None else []
        self._stats_view.setModel(_RowsModel(["Statistic", "Value"], rows, parent=self))
        self._stats_view.setColumnWidth(0, 220)
        self._stats_view.horizontalHeader().setStretchLastSection(True)

    # ---------------------------------------------------------------- differences

    def _step_difference(self, direction: int) -> None:
        positions = self.canvas.diff_positions()
        if not positions:
            return
        target = step_difference(positions, self.canvas.mark_bp, direction)
        if target is None:
            return
        self._show_view()
        self.canvas.show_base(target, zoom_to_letters=True, mark=True)
        self._sync_difference_label()
        self._select_variant_row(target + 1)

    def _sync_difference_label(self) -> None:
        positions = self.canvas.diff_positions() if self.canvas.model is not None else []
        n = len(positions)
        self._prev_diff.setEnabled(n > 0)
        self._next_diff.setEnabled(n > 0)
        to_cons = self.canvas.hl_target == ContigCanvas.TO_CONSENSUS
        what = "consensus" if to_cons else "reference"
        for b in (self._prev_diff, self._next_diff):
            b.setToolTip(("Previous" if b is self._prev_diff else "Next")
                         + f" difference to the {what}")
        mark = self.canvas.mark_bp
        if not n:
            self._diff_label.setText("No differences")
        elif mark is not None:
            k = bisect.bisect_left(positions, mark)
            if k < n and positions[k] == mark:
                self._diff_label.setText(f"Difference {k + 1} of {n:,}")
            else:
                self._diff_label.setText(f"{n:,} differences")
        else:
            self._diff_label.setText(f"{n:,} differences")
        self._diff_label.setToolTip(f"Positions where a shown sequence differs from the {what}")

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

    FILES_MENU_W = 260

    def _build_files_menu(self) -> None:
        """Files the run wrote: what each one is on the left, the end of its name on the
        right, full path in the tooltip; the same width whatever the run is called. Tables and the report are not listed: the tabs show them live."""
        import platform

        from repeatviewer.ui.platform import reveal_in_file_manager

        menu = self._files_menu
        menu.clear()
        menu.setMinimumWidth(self.FILES_MENU_W)
        menu.setToolTipsVisible(True)
        if self._path is None:
            return
        fm = QFontMetrics(menu.font())
        name_w = fm.horizontalAdvance("x" * 16)           # same for every run

        def add(label: str, target: Path, tip: str) -> None:
            short = fm.elidedText(target.name, Qt.TextElideMode.ElideLeft, name_w)
            act = QAction(f"{label}\t{short}", menu)      # tab: name in the right column
            act.setToolTip(f"{tip}\n{target}")
            act.setEnabled(target.is_file())
            act.triggered.connect(lambda _c=False, t=target: self.openRequested.emit(t))
            menu.addAction(act)

        head = menu.addAction("Written by this run")
        head.setEnabled(False)
        details = self._path.parent / self._data.get("details_dir", "")
        for key, label, tip in (
                ("consensus", "Consensus", "Consensus sequence with the reference features"),
                ("mapped", "Placed sequences", "The placed sequences, on the reference strand"),
                ("unused", "Not used", "Sequences that could not be placed"),
                ("raw", "Aligner output", "The aligner's own output, for other tools")):
            name = self._data.get("details", {}).get(key)
            if name:
                add(label, details / name, tip)
        ref = Path(self._data.get("reference", ""))
        if ref.name:
            menu.addSeparator()
            add("Reference", ref, "The reference the sequences were mapped to")
        menu.addSeparator()
        menu.addAction("Copy consensus", self._copy_consensus)
        folder = details if details.is_dir() else self._path
        menu.addAction("Reveal in Finder" if platform.system() == "Darwin"
                       else "Show in File Manager",
                       lambda t=folder: reveal_in_file_manager(t))

    # ---------------------------------------------------------------- right-click on reads

    def _selected_reads(self) -> list[int]:
        c = self.canvas
        order = c.display_order()
        rank = {i: k for k, i in enumerate(order)}
        return sorted((i for i in c.selection if i in rank), key=lambda i: rank[i])

    def _read_menu(self, index: int, global_pos) -> None:
        """Geneious-style right-click on sequences: copy, extract, zoom, show, hide, BLAST."""
        c = self.canvas
        m = c.model
        if m is None:
            return
        menu = QMenu(self)
        picked = self._selected_reads() if index >= 0 else []
        n = len(picked)
        if n:
            one = n == 1
            what = "sequence" if one else f"{n:,} sequences"
            menu.addAction("Copy name" if one else "Copy names", lambda: self._copy_names(picked))
            menu.addAction(f"Copy {what} (FASTA)", lambda: self._copy_fasta(picked))
            menu.addAction(f"Extract {what} to a file…", lambda: self._extract_reads(picked))
            menu.addSeparator()
            if one:
                menu.addAction("Zoom to this sequence", lambda: self._zoom_to_read(picked[0]))
                if not c.one_per_row:
                    menu.addAction("Show in One per row",
                                   lambda: self._open_read_in_rows(picked[0]))
                menu.addAction("Show in the Sequences table",
                               lambda: self._show_in_table(picked[0]))
            menu.addAction(f"Hide {what}", lambda: self._hide_reads(picked))
        hidden = c.user_hidden_count()
        if hidden:
            menu.addAction(f"Show hidden sequences ({hidden:,})", self._show_hidden)
        if n:
            menu.addSeparator()
            menu.addAction(f"BLAST {what}…", lambda: self._blast_reads(picked))
        if not menu.actions():
            return
        menu.exec(global_pos)

    def _mapped_sequences(self) -> dict[str, str]:
        """Placed sequences as written by the run (_mapped.fa: whole sequences, turned to
        the reference strand), read once on first use."""
        cache = getattr(self, "_mapped_cache", None)
        if cache is not None and cache[0] == self._path:
            return cache[1]
        seqs: dict[str, str] = {}
        name = self._data.get("details", {}).get("mapped")
        if self._path is not None and name:
            target = self._path.parent / self._data.get("details_dir", "") / name
            try:
                with open(target, encoding="utf-8", errors="replace") as fh:
                    key, parts = None, []
                    for line in fh:
                        if line.startswith(">"):
                            if key is not None:
                                seqs.setdefault(key, "".join(parts))
                            key, parts = line[1:].split()[0] if line[1:].split() else "", []
                        else:
                            parts.append(line.strip())
                    if key is not None:
                        seqs.setdefault(key, "".join(parts))
            except OSError:
                seqs = {}
        self._mapped_cache = (self._path, seqs)
        return seqs

    def _read_fasta(self, indices: list[int]) -> str:
        """FASTA of the sequences: from _mapped.fa when the run wrote it, else the aligned
        part from the view (without insertions, which the title says)."""
        m = self.canvas.model
        whole = self._mapped_sequences()
        out = []
        for i in indices:
            r = m.reads[i]
            # same id as the run's _mapped.fa (map_reference.fasta_id)
            pl = int(r.data.get("placement", 1) or 1)
            rid = r.name + ("_reversed" if r.strand == "-" else "")
            rid = rid if pl == 1 else f"{rid}_p{pl}"
            seq = whole.get(rid)
            title = f"{rid} mapped to {m.ref_id}:{r.start + 1}-{r.end} ({r.strand})"
            if not seq:
                seq = str(r.data.get("row", "")).replace("-", "")
                title += " aligned part only, insertions not included"
            out.append(f">{title}\n" + "\n".join(seq[k:k + 70] for k in range(0, len(seq), 70)))
        return "\n".join(out) + "\n"

    def _copy_names(self, indices: list[int]) -> None:
        m = self.canvas.model
        QApplication.clipboard().setText("\n".join(m.reads[i].name for i in indices))

    def _copy_fasta(self, indices: list[int]) -> None:
        QApplication.clipboard().setText(self._read_fasta(indices))

    def _extract_reads(self, indices: list[int]) -> None:
        """Selected sequences into <result>_selected.fa next to the result, opened in
        Preview (file-exists question per the app's rule)."""
        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        if self._path is None or not indices:
            return
        proposed = self._path.with_name(f"{self._path.stem}_selected.fa")
        target = confirm_existing_file_path(self, proposed, title="Extract sequences",
                                            required_suffix=".fa")
        if target is None:
            return
        try:
            target.write_text(self._read_fasta(indices), encoding="utf-8")
        except OSError as exc:
            from PySide6.QtWidgets import QMessageBox

            QMessageBox.warning(self, "Extract sequences", f"Could not write {target.name}: {exc}")
            return
        self.fileWritten.emit(target)
        self.openRequested.emit(target)

    def _blast_reads(self, indices: list[int]) -> None:
        """Write the sequences as a query next to the run's other files and hand it to the
        BLAST window (Define inputs)."""
        if self._path is None or not indices:
            return
        folder = self._path.parent / self._data.get("details_dir", "")
        if not folder.is_dir():
            folder = self._path.parent
        target = folder / f"{self._path.stem}_blast_query.fa"
        try:
            target.write_text(self._read_fasta(indices), encoding="utf-8")
        except OSError:
            return
        self.fileWritten.emit(target)
        self.blastRequested.emit(target)

    def _zoom_to_read(self, i: int) -> None:
        r = self.canvas.model.reads[i]
        self.canvas.show_range(r.start, r.end)
        self.canvas.selected = i
        self.canvas.show_read(i)

    def _show_in_table(self, i: int) -> None:
        self._tabs.setCurrentWidget(self._reads_view)
        self._select_read_row(i)

    def _hide_reads(self, indices: list[int]) -> None:
        self.canvas.hide_reads(indices)
        self._apply_filters()

    def _show_hidden(self) -> None:
        self.canvas.show_hidden_reads()
        self._apply_filters()

    def _copy_consensus(self) -> None:
        """The consensus as shown: at the threshold picked in the Consensus section."""
        m = self.canvas.model
        if m is None:
            return
        seq = self.canvas.consensus_seq() or m.consensus
        if not seq:
            return
        t = self._threshold_value
        tag = "" if t is None else ("_majority" if t == 0 else f"_{round(t * 100)}pct")
        QApplication.clipboard().setText(f">{m.ref_id}_consensus{tag}\n{seq}\n")
