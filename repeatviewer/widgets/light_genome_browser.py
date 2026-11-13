"""Light genome browser: a read-only Preview for SAM / BAM / CRAM and GFF / GTF / BED.

Files are read by the track backend vendored from GenomeFlow
(``backends/light_browser/genome_tracks.py``): indexed BAM is fetched window by
window in background threads, SAM and annotations are indexed once. This module
only draws: a ruler, the reference bases (when a reference matches), and the
file as one track (coverage + packed reads, or gene models). Nothing can be
edited; the full editor is GenomeFlow.

The reference is the project's genome when its sequence names match the file;
otherwise the user picks one (``Reference`` button), and the choice is
remembered for that file.
"""

from __future__ import annotations

import json
import math
import re
from pathlib import Path

import numpy as np
from PySide6.QtCore import QEvent, QPointF, QRectF, Qt, QThread, QTimer, Signal
from PySide6.QtGui import (
    QAction,
    QColor,
    QFont,
    QFontMetrics,
    QGuiApplication,
    QKeySequence,
    QContextMenuEvent,
    QPainter,
    QPainterPath,
    QPen,
    QPolygonF,
)
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMenu,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QScrollBar,
    QSizePolicy,
    QSplitter,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QTextBrowser,
    QToolButton,
    QToolTip,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.light_browser import alignment_style as style
from repeatviewer.backends.light_browser import genome_tracks as gt
from repeatviewer.backends.light_browser import reference as refmod

# ---------------------------------------------------------------- file types

ALIGNMENT_SUFFIXES = (".sam", ".sam.gz", ".bam", ".cram")
ANNOTATION_SUFFIXES = (".gff", ".gff3", ".gtf", ".bed", ".gff.gz", ".gff3.gz", ".gtf.gz",
                       ".bed.gz")
LIGHT_BROWSER_SUFFIXES = ALIGNMENT_SUFFIXES + ANNOTATION_SUFFIXES

# Colour-by / Group-by menus for alignments (full labels — no combo elision).
# Tag… stays out of this light UI; insert-size colour is omitted.
_COLOR_MENU_ENTRIES: tuple[tuple[str, str] | None, ...] = (
    (style.COLOR_STRAND, "Strand"),
    None,
    (style.COLOR_FIRST_STRAND, "First-of-pair strand"),
    None,
    (style.COLOR_MAPQ, "Mapping quality"),
    (style.COLOR_IDENTITY, "Identity"),
    None,
    (style.COLOR_NONE, "None"),
)
_COLOR_MODE_LABELS = {
    mode: label for entry in _COLOR_MENU_ENTRIES if entry for mode, label in (entry,)
}
_GROUP_MENU_ENTRIES: tuple[tuple[str, str], ...] = (
    (style.GROUP_NONE, "None"),
    (style.COLOR_STRAND, "Strand"),
    (style.COLOR_FIRST_STRAND, "First-of-pair strand"),
    ("pair_orientation", "Pair orientation"),
    (style.COLOR_MAPQ, "Mapping quality"),
    (style.COLOR_READ_GROUP, "Read group"),
    (style.COLOR_SAMPLE, "Sample"),
)


def is_light_browser_file(path: Path) -> bool:
    return Path(path).name.lower().endswith(LIGHT_BROWSER_SUFFIXES)


def _is_alignment(path: Path) -> bool:
    return Path(path).name.lower().endswith(ALIGNMENT_SUFFIXES)


# Indexing a SAM or GFF is slower than reading the same number of bytes of text,
# so the status-bar bar appears sooner than the 10 MB text-preview threshold.
_LOAD_PROGRESS_MIN_BYTES = 2 * 1024 * 1024
_LOAD_PROGRESS_DELAY_MS = 400


def _is_binary(path: Path) -> bool:
    return Path(path).name.lower().endswith((".bam", ".cram"))


# ---------------------------------------------------------------- settings

_SETTINGS_REFS = "LightGenomeBrowser/references"
_SETTINGS_MODE = "LightGenomeBrowser/text_mode"   # kinds last shown as Text
_SETTINGS_COLOR = "LightGenomeBrowser/color_by"


def _settings():
    from PySide6.QtCore import QSettings

    return QSettings()


def remembered_reference(path: Path) -> Path | None:
    try:
        data = json.loads(str(_settings().value(_SETTINGS_REFS, "{}") or "{}"))
    except (ValueError, TypeError):
        return None
    got = data.get(str(Path(path).resolve()))
    return Path(got) if got else None


def remember_reference(path: Path, fasta: Path | None) -> None:
    s = _settings()
    try:
        data = json.loads(str(s.value(_SETTINGS_REFS, "{}") or "{}"))
    except (ValueError, TypeError):
        data = {}
    key = str(Path(path).resolve())
    if fasta is None:
        data.pop(key, None)
    else:
        data[key] = str(Path(fasta).resolve())
    if len(data) > 300:                       # keep the newest few hundred
        data = dict(list(data.items())[-300:])
    s.setValue(_SETTINGS_REFS, json.dumps(data))


def _text_mode_kinds() -> set[str]:
    raw = _settings().value(_SETTINGS_MODE, "") or ""
    return {k for k in str(raw).split(",") if k}


def _set_text_mode(kind: str, on: bool) -> None:
    kinds = _text_mode_kinds()
    (kinds.add if on else kinds.discard)(kind)
    _settings().setValue(_SETTINGS_MODE, ",".join(sorted(kinds)))


def _load_project():
    try:
        from repeatviewer.views.settings_dialog import SettingsDialog

        return SettingsDialog.load_project_from_settings()
    except Exception:  # noqa: BLE001 - no project: browse without a reference
        return None


# ---------------------------------------------------------------- loading


# Worker threads are not parented to the panel: a preview page may be deleted while
# a large file is still being indexed, and a QThread must outlive its run().
_LIVE_THREADS: set = set()


def _start_thread(thread: QThread) -> None:
    _LIVE_THREADS.add(thread)

    def _done(t=thread) -> None:
        _LIVE_THREADS.discard(t)
        t.deleteLater()

    thread.finished.connect(_done)
    thread.start()


class _LoadResult:
    def __init__(self) -> None:
        self.track = None
        self.error = ""
        self.reference: refmod.ReferenceMatch | None = None
        self.tried: list[refmod.ReferenceMatch] = []
        self.candidates: list[tuple[str, Path]] = []
        self.lengths: dict[str, int] = {}
        self.first: dict[str, tuple[int, int]] = {}     # sequence -> first item span
        self.order: list[str] = []                       # sequences, those with data first
        self.forced = False                              # the reference was given


class _Loader(QThread):
    """Open the track, match a reference and find where the data starts, off the GUI
    thread (a large SAM or GFF takes seconds to index; a new .fai takes a pass)."""

    loaded = Signal(object, object)       # token, _LoadResult

    def __init__(self, token, path: Path, project, reference: Path | None,
                 forced: bool, parent=None) -> None:
        super().__init__(parent)
        self.token = token
        self.path = path
        self.project = project
        self.reference = reference
        self.forced = forced

    def run(self) -> None:  # noqa: D401 - QThread API
        res = _LoadResult()
        res.forced = self.forced
        try:
            track = gt.open_track(self.path)
            res.track = track
            names = refmod.track_sequence_names(track)
            track_lengths, auth_lengths = refmod.track_sequence_lengths(track)
            if self.forced:
                if self.reference is not None:
                    m = refmod.match_reference(
                        self.reference,
                        names,
                        "Chosen for this file",
                        track_lengths=track_lengths,
                        authoritative_lengths=auth_lengths,
                    )
                    res.reference = m
                    res.tried = [m]
            else:
                res.candidates = refmod.reference_candidates(self.project, self.reference)
                res.reference, res.tried = refmod.choose_reference(
                    names,
                    res.candidates,
                    track_lengths=track_lengths,
                    authoritative_lengths=auth_lengths,
                )
            res.lengths = refmod.sequence_lengths(track, res.reference)
            res.first = _first_items(track, names)
            res.order = refmod.genome_order(track, res.reference)
        except gt.TrackError as exc:
            res.error = str(exc)
        except Exception as exc:  # noqa: BLE001 - shown in the preview instead of a crash
            res.error = f"{type(exc).__name__}: {exc}"
        self.loaded.emit(self.token, res)


def _first_items(track, names: list[str]) -> dict[str, tuple[int, int]]:
    """Span of the first item of each sequence that has any."""
    first: dict[str, tuple[int, int]] = {}
    bam = getattr(track, "_bam", None)
    if bam is not None:                      # indexed BAM: index statistics + one read
        mapped: dict[str, int] = {}
        try:
            for st in bam.get_index_statistics():
                mapped[st.contig] = int(st.mapped)
        except (AttributeError, ValueError):
            pass
        with track._lock:
            for n in names:
                if mapped and not mapped.get(n):
                    continue
                try:
                    a = next(iter(bam.fetch(n)), None)
                except (ValueError, OSError):
                    a = None
                if a is not None:
                    first[n] = (a.reference_start, a.reference_end or a.reference_start + 1)
                if len(first) >= 1 and not mapped:
                    break
    else:
        for n in names:
            ix = track._chroms.get(n) if hasattr(track, "_chroms") else None
            if ix is not None and len(ix):
                first[n] = (int(ix.starts[0]), int(ix.ends[0]))
    return first


class _RefJob(QThread):
    """Match a newly chosen reference against the open track (may build a .fai)."""

    done = Signal(object, object)      # ReferenceMatch | None, error text

    def __init__(self, fasta: Path | None, track, parent=None) -> None:
        super().__init__(parent)
        self.fasta = fasta
        self.track = track

    def run(self) -> None:  # noqa: D401 - QThread API
        out = None
        err = ""
        if self.fasta is not None:
            try:
                names = refmod.track_sequence_names(self.track)
                track_lengths, auth_lengths = refmod.track_sequence_lengths(self.track)
                out = refmod.match_reference(
                    self.fasta,
                    names,
                    "Chosen for this file",
                    track_lengths=track_lengths,
                    authoritative_lengths=auth_lengths,
                )
            except (OSError, gt.TrackError) as exc:
                err = str(exc)
        self.done.emit(out, err)


# ---------------------------------------------------------------- drawing helpers

RULER_H = 28
SEQ_H = 16
COV_H = 22
JUNC_H = 36
READ_H = 9
READ_H_SQUISHED = 4
READ_H_COLLAPSED = 2
READ_GAP = 2
ANN_ROW_H = 26
ANN_ROW_H_SMALL = 12

LAYOUT_EXPANDED = "expanded"
LAYOUT_SQUISHED = "squished"
LAYOUT_COLLAPSED = "collapsed"
LAYOUT_READ_H = {
    LAYOUT_EXPANDED: READ_H,
    LAYOUT_SQUISHED: READ_H_SQUISHED,
    LAYOUT_COLLAPSED: READ_H_COLLAPSED,
}
MARGIN_TOP = 2
MAX_READS_DRAWN = 30_000
MAX_FEATURES_DRAWN = 25_000
ZOOM_STEPS = 100                # zoom slider: Fit (0) to single bases
ONE_PER_ROW_LIMIT = 20_000      # names + one row each; above this, stay packed
NAME_W = 168                    # sequence-name column in one-per-row
NAME_ROW_H = 18                 # row height when the name is shown
FEAT_LANE_H = 22                # reference features from the run JSON
SEQ_WINDOW_MAX = 200_000        # bases fetched for the sequence lane / mismatches
MIN_BPP = 1 / 24                # at most 24 px per base
BASE_COLORS = {"A": "#3ca34a", "C": "#2f6fd6", "G": "#d29a12", "T": "#d33a2f",
               "N": "#9a9a9a"}
GENE_COLOR = "#2f5f9e"
_CHEVRON_GAP = 18


def _nice_step(span_bp: float, px: int) -> int:
    """Ruler tick step: about one label per 110 px, rounded to 1/2/5 x 10^n."""
    raw = span_bp * 110.0 / max(1, px)
    exp = 10 ** math.floor(math.log10(max(raw, 1)))
    for m in (1, 2, 5, 10):
        if m * exp >= raw:
            return int(max(1, m * exp))
    return int(10 * exp)


def _fmt_len(n: int) -> str:
    """A sequence length: 2.28 Mb, 12.4 kb, or 800 bp."""
    n = max(0, int(n))
    if n >= 1_000_000:
        text = f"{n / 1_000_000:.2f}".rstrip("0").rstrip(".")
        return f"{text} Mb"
    if n >= 1_000:
        text = f"{n / 1_000:.1f}".rstrip("0").rstrip(".")
        return f"{text} kb"
    return f"{n:,} bp"


_PLACEMENT_SUFFIX = re.compile(r"_p\d+$")


def _has_bases(seq: str) -> bool:
    return any(c in "ACGTNacgtn" for c in seq)


def _placement_name(rid: str, strand: str) -> str:
    """Query name encoded in a ``_mapped.fa`` record id (``name``, ``name_reversed``, ``_pN``)."""
    rid = _PLACEMENT_SUFFIX.sub("", rid.strip().split()[0])
    if strand == "-" and rid.endswith("_reversed"):
        rid = rid[: -len("_reversed")]
    return rid


def _bases_from_row(read: dict) -> str:
    """Aligned query bases from a run row, with insertions put back in."""
    row = str(read.get("row") or "")
    if not _has_bases(row):
        return ""
    start = int(read.get("start", 0))
    pending: dict[int, list[str]] = {}
    for item in read.get("ins") or []:
        if isinstance(item, (list, tuple)) and len(item) >= 2 and _has_bases(str(item[1])):
            pending.setdefault(int(item[0]), []).append(str(item[1]))
    out: list[str] = []
    for i, base in enumerate(row):
        out.extend(pending.pop(start + i, ()))
        if base not in "-?":
            out.append(base)
    out.extend(base for bases in pending.values() for base in bases)
    return "".join(out)


def _fmt_bp(v: int, step: int) -> str:
    if step >= 1_000_000 and v % 1_000_000 == 0:
        return f"{v // 1_000_000:,} Mb"
    if step >= 1_000 and v % 1_000 == 0:
        return f"{v // 1_000:,} kb"
    return f"{v:,}"


def _mismatches(it, ref: str, ref0: int, a: int, b: int) -> list[tuple[int, str]]:
    """Read bases differing from the reference window ``ref`` (starting at ``ref0``)."""
    try:
        seq = it.seq
    except Exception:  # noqa: BLE001
        return []
    if not seq or not ref:
        return []
    out: list[tuple[int, str]] = []
    ref1 = ref0 + len(ref)
    for s, e, q in (it.blocks or [(it.start, it.end, 0)]):
        lo, hi = max(s, a, ref0), min(e, b, ref1)
        if hi <= lo:
            continue
        k0 = q + lo - s
        qs = seq[k0:k0 + hi - lo].upper()
        rs = ref[lo - ref0:hi - ref0]
        for j, (x, y) in enumerate(zip(qs, rs)):
            if x != y and x != "N" and y != "N":
                out.append((lo + j, x))
    return out


# ---------------------------------------------------------------- canvas


class LightBrowserCanvas(QWidget):
    """The drawn area. Coordinates: ``start`` (bp, float) at x = 0 and ``bpp`` bp per
    pixel; rows scroll vertically with ``vscroll`` (px)."""

    locusChanged = Signal()
    rowsChanged = Signal()
    selectionChanged = Signal()
    readNameMenu = Signal(object, object)   # TrackItem, global position
    readPicked = Signal(object)             # left-click on a sequence name
    _tileReady = Signal()

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setAttribute(Qt.WidgetAttribute.WA_OpaquePaintEvent, True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setMinimumHeight(120)
        self.track = None
        self.box: dict = {}          # plain holder: the track is released after deletion
        self.chrom = ""
        self.length = 1
        self.start = 0.0
        self.bpp = 1.0
        self.vscroll = 0
        self.color_by = style.COLOR_STRAND
        self.run_features: list = []
        self.one_per_row = False
        self.sort_key = "start"
        self.layout_mode = LAYOUT_EXPANDED
        self.group_by = style.GROUP_NONE
        self.min_identity = 0.0
        self.min_mapq = 0
        self.hide_secondary = False
        self.show_features = True
        self.hidden_feature_types: set[str] = set()
        self.show_coverage = True
        self.show_junctions = False
        self.cov_log = False
        self.cov_low: int | None = None
        self.cov_high: int | None = None
        self.show_mismatches = True
        self.reference: refmod.ReferenceMatch | None = None
        self._ref_find = None
        self._seq_cache: tuple[str, int, str] | None = None   # (chrom, start, bases)
        self._drag: tuple[QPointF, float, int] | None = None
        self.sel: tuple[int, int] | None = None        # region dragged on the ruler
        self._sel_from: float | None = None
        self._sel_click_new = False
        self.bookmarks: list[dict] = []
        self._pinch_anchor: float | None = None
        self._content_h = 0
        self._hits: list[tuple[QRectF, object]] = []
        self._name_hits: list[tuple[QRectF, object]] = []
        self._names_on = False
        self.picked = None
        self._hover = None
        self._ruler_hover: int | None = None
        self._context_at: QPointF | None = None
        self._menu_pending = False
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.DefaultContextMenu)
        self._user_hidden: set[tuple[str, str, int]] = set()
        self._font = QFont(self.font())
        self._font.setPointSizeF(max(8.0, self._font.pointSizeF() - 2))
        self._mono = QFont("Menlo")
        self._mono.setStyleHint(QFont.StyleHint.Monospace)
        self._mono.setPointSizeF(self._font.pointSizeF())
        self._repaint = QTimer(self)
        self._repaint.setSingleShot(True)
        self._repaint.setInterval(16)
        self._repaint.timeout.connect(self.update)
        self._tileReady.connect(self._repaint.start)
        self.grabGesture(Qt.GestureType.PinchGesture)

    # ------------------------------------------------------------ model
    def set_track(self, track, reference: refmod.ReferenceMatch | None) -> None:
        old = self.track
        if old is not None and old is not track:
            release_track(old)
        self.track = track
        self._user_hidden.clear()
        self._named_cache_key = None
        self.box["track"] = track
        if track is not None and hasattr(track, "on_tile_ready"):
            track.on_tile_ready = self._tileReady.emit
        self.set_reference(reference)

    def set_reference(self, reference: refmod.ReferenceMatch | None) -> None:
        self.reference = reference
        self._ref_find = refmod.NameMatcher(reference.lengths) if reference else None
        self._seq_cache = None
        self.update()

    def set_chrom(self, chrom: str, length: int, start: int | None = None,
                  end: int | None = None) -> None:
        if chrom != self.chrom:
            self.set_selection(None)
        self.chrom = chrom
        self.length = max(1, int(length))
        self.vscroll = 0
        if start is None or end is None:
            start, end = 0, self.length
        self.show_range(start, end)

    def show_range(self, start: float, end: float) -> None:
        w = max(1, self._data_w())
        span = max(w * MIN_BPP, min(float(end - start), self.length * 1.0))
        self.bpp = span / w
        mid = (start + end) / 2.0
        self.start = mid - span / 2.0
        self._clamp()
        self.update()
        self.locusChanged.emit()

    def set_selection(self, sel: tuple[int, int] | None) -> None:
        if sel != self.sel:
            self.sel = sel
            self.update()
            self.selectionChanged.emit()

    def go_to_selection(self) -> None:
        if self.sel is None:
            return
        a, b = self.sel
        self.set_selection(None)
        self.show_range(a, max(b, a + 1))

    def _in_ruler(self, pos: QPointF) -> bool:
        return pos.y() < MARGIN_TOP + RULER_H and pos.x() >= self._name_w()

    def _name_w(self) -> int:
        return NAME_W if self._names_on else 0

    def _data_w(self) -> int:
        return max(1, self.width() - self._name_w())

    def view_range(self) -> tuple[int, int]:
        a = max(0, int(math.floor(self.start)))
        b = min(self.length, int(math.ceil(self.start + self.bpp * self._data_w())))
        return a, max(a + 1, b)

    def _clamp(self) -> None:
        w = max(1, self._data_w())
        self.bpp = min(max(self.bpp, MIN_BPP), max(MIN_BPP, self.length / w))
        span = self.bpp * w
        self.start = min(max(0.0, self.start), max(0.0, self.length - span))

    def x_of(self, bp: float) -> float:
        return self._name_w() + (bp - self.start) / self.bpp

    def bp_at(self, x: float) -> float:
        return self.start + (x - self._name_w()) * self.bpp

    def _view_coverage(self, track, chrom: str):
        """Coverage for the visible pixels, edges matching ``bp_at`` / ``x_of``."""
        w = max(1, self._data_w())
        left = float(self.start)
        right = left + float(self.bpp) * w
        if right <= left:
            return None, left, right, w
        # ChromIndex profile: resample on the exact continuous window (not
        # floor/ceil view_range), so density lines up with mapping ticks.
        index = getattr(track, "index", None)
        if callable(index):
            ix = index(chrom)
            if ix is not None and hasattr(ix, "coverage_profile"):
                span = right - left
                if span / w < gt.COV_BIN:
                    a = max(0, int(math.floor(left)))
                    b = min(self.length, max(a + 1, int(math.ceil(right))))
                    return ix.window_coverage(a, b, w), left, right, w
                return (
                    gt.window_profile(ix.coverage_profile(), gt.COV_BIN, left, right, w),
                    left,
                    right,
                    w,
                )
        a = max(0, int(math.floor(left)))
        b = min(self.length, max(a + 1, int(math.ceil(right))))
        return track.coverage(chrom, a, b, w), left, right, w

    def _paint_coverage_bars(
        self,
        p: QPainter,
        cov,
        left: float,
        right: float,
        top: float,
        height: float,
        scale: float,
        *,
        raw: np.ndarray | None = None,
    ) -> None:
        """Draw coverage on the same X scale as reads (``x_of`` / ``bp_at``)."""
        n = len(cov)
        if n <= 0 or right <= left:
            return
        span = float(right - left)
        raw_vals = raw if raw is not None else cov
        for i, v in enumerate(cov):
            if v <= 0:
                continue
            x0 = self.x_of(left + i * span / n)
            x1 = self.x_of(left + (i + 1) * span / n)
            width = max(1.0, x1 - x0)
            hh = max(1.0, float(v) * scale)
            p.drawRect(QRectF(x0, top + height - hh, width, hh))
            depth = float(raw_vals[i]) if i < len(raw_vals) else 0.0
            if self.cov_low is not None and 0 < depth < self.cov_low:
                p.fillRect(QRectF(x0, top + height - 3, width, 3), QColor("#c0392b"))
            elif self.cov_high is not None and depth > self.cov_high:
                p.fillRect(QRectF(x0, top + height - 3, width, 3), QColor("#e67e22"))

    def _paint_junctions(self, p: QPainter, top: float, a: int, b: int, w: int,
                         fg: QColor, mid: QColor) -> None:
        """Splice-junction arcs from CIGAR N / long gaps in the loaded reads."""
        t = self.track
        if t is None or not hasattr(t, "junctions"):
            return
        try:
            rows = t.junctions(self.chrom, a, b)
        except Exception:  # noqa: BLE001 - missing tiles / no junctions
            rows = None
        if rows is None or len(rows) == 0:
            p.setPen(QPen(mid, 1))
            p.drawLine(self._name_w(), top + JUNC_H - 2, w, top + JUNC_H - 2)
            return
        base = top + JUNC_H - 4
        peak = max(1, int(np.max(rows[:, 2:5]))) if len(rows) else 1
        for js, je, plus, minus, unk in rows:
            total = int(plus) + int(minus) + int(unk)
            if total <= 0:
                continue
            x0, x1 = self.x_of(int(js)), self.x_of(int(je))
            if x1 < self._name_w() or x0 > w:
                continue
            span = max(4.0, x1 - x0)
            lift = 6.0 + (JUNC_H - 12) * min(1.0, total / peak)
            if plus >= minus and plus >= unk:
                col = QColor("#2a62c9")
            elif minus > plus and minus >= unk:
                col = QColor("#c0392b")
            else:
                col = QColor("#7f8c8d")
            col.setAlpha(160)
            pen = QPen(col, 1.0 + min(2.5, total / max(1.0, peak) * 2.5))
            p.setPen(pen)
            p.setBrush(Qt.BrushStyle.NoBrush)
            path = QPainterPath(QPointF(x0, base))
            path.cubicTo(
                QPointF(x0 + span * 0.15, base - lift),
                QPointF(x1 - span * 0.15, base - lift),
                QPointF(x1, base),
            )
            p.drawPath(path)
        p.setPen(QPen(mid, 1))
        p.drawLine(self._name_w(), top + JUNC_H - 2, w, top + JUNC_H - 2)
        p.setPen(fg)
        p.drawText(
            QRectF(self._name_w() + 4, top, 120, 14),
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
            "Junctions",
        )

    def _zoom_focus_bp(self) -> float | None:
        """Where toolbar zoom should stay: the selected region, else the picked sequence."""
        if self.sel is not None:
            a, b = self.sel
            return (a + b) / 2.0
        it = self.picked
        if it is None:
            return None
        return (int(getattr(it, "start", 0)) + int(getattr(it, "end", 0))) / 2.0

    def zoom(self, factor: float, x: float | None = None) -> None:
        """``factor`` > 1 zooms in. A pointer ``x`` stays put; otherwise zoom on the selection."""
        focus = self._zoom_focus_bp() if x is None else None
        if focus is not None:
            self.bpp /= factor
            self._clamp()
            self.start = focus - self.bpp * self._data_w() / 2.0
            self._clamp()
            self.update()
            self.locusChanged.emit()
            return
        if x is None:
            x = self._name_w() + self._data_w() / 2.0
        anchor = self.bp_at(x)
        self.bpp /= factor
        self._clamp()
        self.start = anchor - (x - self._name_w()) * self.bpp
        self._clamp()
        self.update()
        self.locusChanged.emit()

    def pan_px(self, dx: float) -> None:
        self.start -= dx * self.bpp
        self._clamp()
        self.update()
        self.locusChanged.emit()

    def scroll_rows(self, dy: int) -> None:
        self.set_vscroll(self.vscroll + dy)

    def set_vscroll(self, v: int) -> None:
        v = max(0, min(int(v), max(0, self._content_h - self._rows_area_h())))
        if v != self.vscroll:
            self.vscroll = v
            self.update()
            self.rowsChanged.emit()

    def resizeEvent(self, event) -> None:  # noqa: N802
        old_w = max(1, event.oldSize().width() - self._name_w())
        if event.oldSize().width() > 0:
            # keep the same bases in view
            span = self.bpp * old_w
            self.bpp = span / max(1, self._data_w())
        self._clamp()
        super().resizeEvent(event)
        self.rowsChanged.emit()

    # ------------------------------------------------------------ geometry
    def _seq_lane_h(self) -> int:
        return SEQ_H if self.reference is not None and self._ref_name() else 0

    def _hide_key(self, it) -> tuple[str, str, int]:
        return (self.chrom, str(getattr(it, "name", "") or ""), int(getattr(it, "start", -1)))

    def hide_read(self, it) -> None:
        """Right-click > Hide. This session only; it is not written to the file."""
        self._user_hidden.add(self._hide_key(it))
        if self._same_read(self.picked, it):
            self.picked = None
        if self._same_read(self._hover, it):
            self._hover = None
        self._named_cache_key = None
        self.update()

    def show_hidden_reads(self) -> None:
        if not self._user_hidden:
            return
        self._user_hidden.clear()
        self._named_cache_key = None
        self.update()

    def user_hidden_count(self) -> int:
        return len(self._user_hidden)

    def read_passes(self, it) -> bool:
        if self._hide_key(it) in self._user_hidden:
            return False
        if self.hide_secondary and getattr(it, "secondary", False):
            return False
        mq = getattr(it, "mapq", None)
        if self.min_mapq and mq is not None and int(mq) >= 0 and int(mq) < self.min_mapq:
            return False
        if self.min_identity > 0:
            score = getattr(it, "score", None)
            if score is not None and not callable(score) and float(score) < self.min_identity:
                return False
        return True

    def _scaffold_items(self):
        """Every alignment already in memory on this scaffold, or None for a streamed file."""
        track = self.track
        if track is None or not self.chrom:
            return None
        index = getattr(track, "index", None)
        if not callable(index):
            return None
        ix = index(self.chrom)
        if ix is None or not hasattr(ix, "items"):
            return None
        return ix.items

    def _names_wanted(self) -> bool:
        """One per row lists the whole scaffold, until that list passes 20,000."""
        if not self.one_per_row or self.track is None or not self.chrom:
            return False
        if getattr(self.track, "kind", None) != gt.KIND_ALIGNMENT:
            return False
        items = self._scaffold_items()
        n = len(items) if items is not None else self.track.count_in(self.chrom, 0, self.length)
        return 0 < n <= ONE_PER_ROW_LIMIT

    def one_per_row_availability(self) -> tuple[bool, str]:
        """Whether One per row can list this scaffold, and the reason shown on Layout."""
        limit = f"{ONE_PER_ROW_LIMIT:,}"
        about = (
            "Packed shares rows. One per row gives each sequence its own row and its name. "
            f"It is available for up to {limit} sequences on a scaffold."
        )
        if self.track is None or not self.chrom:
            return False, about
        if getattr(self.track, "kind", None) != gt.KIND_ALIGNMENT:
            return False, "One per row is for alignments."
        window = getattr(self.track, "read_window", None)
        if window and self.length > int(window):
            return False, (
                f"{about} This file only loads {int(window) // 1000:,} kb at a time, "
                "so the whole scaffold cannot be listed."
            )
        items = self._scaffold_items()
        if items is not None:
            n = len(items)
        else:
            n = int(self.track.count(self.chrom)) if hasattr(self.track, "count") else -1
        if n < 0:
            return False, f"{about} This file does not list every sequence."
        if n > ONE_PER_ROW_LIMIT:
            return False, f"{about} This scaffold has {n:,}."
        if n <= 0:
            return False, "No sequences on this scaffold."
        return True, about

    def _named_in_view(self, _a: int, _b: int) -> list:
        """Every sequence on this scaffold that passes the filters.

        One per row keeps the full list while the window moves. A sequence whose
        alignment is off screen still has its row."""
        track = self.track
        if track is None:
            return []
        key = (
            id(track), self.chrom, self.sort_key, self.min_identity,
            self.min_mapq, self.hide_secondary, frozenset(self._user_hidden),
        )
        if key == getattr(self, "_named_cache_key", None):
            return self._named_cache
        raw = self._scaffold_items()
        if raw is None:
            raw = [it for it, _row in track.fetch(self.chrom, 0, max(1, self.length))]
        items = [it for it in raw if self.read_passes(it)]
        if len(items) > ONE_PER_ROW_LIMIT:
            items = []
        else:
            items = self._sort_reads(items)
        self._named_cache_key = key
        self._named_cache = items
        return items

    def _sort_reads(self, items: list) -> list:
        key = self.sort_key

        def order(it):
            start = int(getattr(it, "start", 0))
            if key == "strand":
                return (0 if getattr(it, "strand", "+") == "+" else 1, start)
            if key == "identity":
                score = getattr(it, "score", None)
                score = 0.0 if score is None or callable(score) else float(score)
                return (-score, start)
            if key == "name":
                return (str(getattr(it, "name", "")), start)
            return (start, int(getattr(it, "end", start)))

        return sorted(items, key=order)

    def _feature_lane_h(self) -> int:
        if not self.show_features:
            return 0
        if not any(self._feature_visible(f) for f in self.run_features):
            return 0
        return FEAT_LANE_H

    def _feature_visible(self, feat: dict) -> bool:
        kind = str(feat.get("feature_type") or feat.get("type") or "other")
        return kind not in self.hidden_feature_types

    def _read_body_h(self) -> int:
        return LAYOUT_READ_H.get(self.layout_mode, READ_H)

    def _rows_top(self) -> int:
        top = MARGIN_TOP + RULER_H + self._feature_lane_h() + self._seq_lane_h()
        if (
            self.show_coverage
            and self.track is not None
            and self.track.kind == gt.KIND_ALIGNMENT
        ):
            top += COV_H + 4
        if (
            self.show_junctions
            and self.track is not None
            and self.track.kind == gt.KIND_ALIGNMENT
        ):
            top += JUNC_H + 4
        return top

    def _rows_area_h(self) -> int:
        return max(0, self.height() - self._rows_top())

    def content_height(self) -> int:
        return self._content_h

    def _ref_name(self) -> str | None:
        if self._ref_find is None or not self.chrom:
            return None
        return self._ref_find(self.chrom)

    def _ref_window(self, a: int, b: int) -> tuple[int, str] | None:
        """Reference bases covering [a, b) (cached with some padding)."""
        name = self._ref_name()
        if name is None or b - a > SEQ_WINDOW_MAX:
            return None
        c = self._seq_cache
        if c is not None and c[0] == name and c[1] <= a and c[1] + len(c[2]) >= b:
            return c[1], c[2]
        entry = None
        try:
            entry = gt.fasta_index(self.reference.path).get(name)
        except (OSError, gt.TrackError):
            entry = None
        if entry is None:
            return None
        pad = (b - a)
        s0, s1 = max(0, a - pad), min(entry.length, b + pad)
        try:
            bases = gt.fetch_sequence(self.reference.path, entry, s0, s1)
        except OSError:
            return None
        self._seq_cache = (name, s0, bases)
        return s0, bases

    # ------------------------------------------------------------ painting
    def paintEvent(self, _event) -> None:  # noqa: N802
        p = QPainter(self)
        pal = self.palette()
        bg = pal.base().color()
        fg = pal.text().color()
        mid = pal.mid().color()
        p.fillRect(self.rect(), bg)
        p.setFont(self._font)
        self._hits = []
        self._name_hits = []
        self._names_on = self._names_wanted()
        if self._names_on:
            p.fillRect(QRectF(0, 0, NAME_W, self.height()), pal.window().color())
            p.setPen(QPen(mid, 1))
            p.drawLine(NAME_W - 1, 0, NAME_W - 1, self.height())
        if self.track is None or not self.chrom:
            p.end()
            return
        w = self.width()
        a, b = self.view_range()
        y = MARGIN_TOP
        self._paint_ruler(p, y, w, fg, mid)
        self._paint_marks(p, y)
        y += RULER_H
        if self._feature_lane_h():
            self._paint_run_features(p, y, a, b)
            y += FEAT_LANE_H
        ref = self._ref_window(a, b) if self.bpp <= 2 else None
        if self._seq_lane_h():
            self._paint_sequence(p, y, a, b, ref, fg, mid)
            y += SEQ_H
        drawing = getattr(self.track, "drawing", None)
        ctx = drawing() if callable(drawing) else _NoCtx()
        before = self._content_h
        with ctx:
            if self.track.kind == gt.KIND_ALIGNMENT:
                self._paint_alignments(p, y, a, b, w, ref, fg, mid, bg)
            else:
                self._paint_annotations(p, y, a, b, w, fg, mid, bg)
        if self.sel is not None:
            # the region dragged on the ruler, carried down through the view
            p.save()
            p.setClipRect(QRectF(self._name_w(), 0, self._data_w(), self.height()))
            x0, x1 = self.x_of(self.sel[0]), self.x_of(self.sel[1])
            p.fillRect(QRectF(x0, MARGIN_TOP, max(1.0, x1 - x0), self.height()),
                       QColor(61, 123, 224, 28))
            p.setPen(QPen(QColor("#3d7be0"), 1))
            p.drawLine(QPointF(x0, MARGIN_TOP), QPointF(x0, self.height()))
            p.drawLine(QPointF(x1, MARGIN_TOP), QPointF(x1, self.height()))
            p.restore()
        if self._names_on:
            self._paint_left_header(p, fg)
        self._paint_ruler_hover(p)
        p.end()
        if self._content_h != before:
            # the scroll bar follows the rows now known (after the paint, not in it)
            if self.vscroll > max(0, self._content_h - self._rows_area_h()):
                QTimer.singleShot(0, lambda: self.set_vscroll(self.vscroll))
            QTimer.singleShot(0, self.rowsChanged.emit)

    def _paint_left_header(self, p: QPainter, fg: QColor) -> None:
        """Name column, when One per row is on: only how much of the scaffold is in view."""
        a, b = self.view_range()
        text = f"{max(0, b - a):,} bp in view"
        width = NAME_W - 14
        fm = QFontMetrics(self._font)
        p.setFont(self._font)
        p.setPen(fg)
        p.drawText(QRectF(6, MARGIN_TOP, width, RULER_H),
                   Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                   fm.elidedText(text, Qt.TextElideMode.ElideMiddle, width))

    def _paint_ruler(self, p: QPainter, y: int, w: int, fg: QColor, mid: QColor) -> None:
        left = self._name_w()
        span = self.bpp * max(1, w - left)
        step = _nice_step(span, max(1, w - left))
        a = int(self.start // step) * step
        p.setPen(QPen(mid, 1))
        base = y + RULER_H - 6
        p.drawLine(left, base, w, base)
        fm = QFontMetrics(self._font)
        p.setPen(fg)
        v = a
        while v <= self.start + span:
            x = self.x_of(v)
            if x >= left - 1:
                p.drawLine(QPointF(x, base - 4), QPointF(x, base))
                label = _fmt_bp(v, step) if v else "0"
                tw = fm.horizontalAdvance(label)
                p.drawText(QPointF(max(float(left + 2), x - tw / 2), base - 7), label)
            v += step
        # minor ticks
        p.setPen(QPen(mid, 1))
        minor = step / 5
        if minor >= 1 and minor / self.bpp >= 6:
            v = a
            while v <= self.start + span:
                x = self.x_of(v)
                p.drawLine(QPointF(x, base - 2), QPointF(x, base))
                v += minor

    def _set_ruler_hover(self, pos: QPointF | None) -> None:
        bp = None
        if pos is not None and self._in_ruler(pos) and self.length > 0:
            bp = max(0, min(self.length - 1, int(self.bp_at(pos.x()))))
        if bp == self._ruler_hover:
            return
        self._ruler_hover = bp
        self.update()

    def _paint_ruler_hover(self, p: QPainter) -> None:
        """Coordinate under the pointer, on the same line as the scale labels.

        While a range is being dragged, both ends are labelled.
        """
        if self._sel_from is not None and self.sel is not None and self.sel[1] > self.sel[0]:
            a, b = self.sel
            self._paint_ruler_coords(p, ((float(a), f"{a + 1:,}"), (float(b), f"{b:,}")))
            return
        if self._ruler_hover is None:
            return
        text = f"{self._ruler_hover + 1:,}"
        self._paint_ruler_coords(p, ((self._ruler_hover + 0.5, text),))

    def _paint_ruler_coords(self, p: QPainter, anchors: tuple[tuple[float, str], ...]) -> None:
        """Blue coordinates on the scale baseline, each on a white patch."""
        fm = QFontMetrics(self._font)
        baseline = MARGIN_TOP + RULER_H - 13
        pad_x = 2
        gap = 4
        placed: list[list] = []
        for bp, text in anchors:
            tw = fm.horizontalAdvance(text)
            placed.append([self.x_of(bp) - tw / 2.0, tw, text])
        if len(placed) == 2:
            placed.sort(key=lambda item: item[0])
            left, right = placed
            overlap = (left[0] + left[1] + gap) - right[0]
            if overlap > 0:
                left[0] -= overlap / 2.0
                right[0] += overlap / 2.0
        limit_l = float(self._name_w() + 1)
        limit_r = float(self.width() - 1)
        bg = self.palette().base().color()
        for item in placed:
            x, tw, _text = item
            x = max(limit_l, min(limit_r - tw, x))
            item[0] = x
            p.fillRect(QRectF(x - pad_x, baseline - fm.ascent(), tw + pad_x * 2, fm.height()), bg)
        p.setFont(self._font)
        p.setPen(QColor("#2a62c9"))
        for x, _tw, text in placed:
            p.drawText(QPointF(x, baseline), text)

    def _paint_marks(self, p, y) -> None:
        """Bookmarks of this sequence: a coloured strip at the foot of the ruler."""
        for bm in self.bookmarks:
            if bm.get("chrom") != self.chrom:
                continue
            x0, x1 = self.x_of(bm["start"]), self.x_of(bm["end"])
            if x1 < 0 or x0 > self.width():
                continue
            p.fillRect(QRectF(x0, y + RULER_H - 4, max(3.0, x1 - x0), 3),
                       QColor(217, 83, 79, 200))

    def _paint_sequence(self, p, y, a, b, ref, fg, mid) -> None:
        h = SEQ_H - 3
        if ref is None:
            p.setPen(QPen(mid, 1))
            p.drawLine(QPointF(self.x_of(0), y + h / 2), QPointF(self.x_of(self.length),
                                                                 y + h / 2))
            return
        r0, bases = ref
        ppb = 1.0 / self.bpp
        if ppb < 1:
            # more than one base per pixel: colours would only be noise
            p.setPen(QPen(mid, 1))
            p.drawLine(QPointF(0, y + h / 2), QPointF(self.width(), y + h / 2))
            return
        letters = ppb >= 8
        if letters:
            p.setFont(self._mono)
        fm = QFontMetrics(p.font())
        for pos in range(a, min(b, r0 + len(bases))):
            if pos < r0:
                continue
            ch = bases[pos - r0]
            col = QColor(BASE_COLORS.get(ch, "#9a9a9a"))
            x0 = self.x_of(pos)
            if letters:
                p.setPen(col)
                p.drawText(QPointF(x0 + (ppb - fm.horizontalAdvance(ch)) / 2,
                                   y + h - 2), ch)
            else:
                p.fillRect(QRectF(x0, y + 2, max(1.0, ppb), h - 3), col)
        p.setFont(self._font)

    def _paint_run_features(self, p: QPainter, y: int, a: int, b: int) -> None:
        """Reference features from the mapping run, on the genome-browser ruler."""
        for feat in self.run_features:
            if not self._feature_visible(feat):
                continue
            start = int(feat.get("start", 0))
            end = int(feat.get("end", start))
            if end < a or start > b:
                continue
            x0 = self.x_of(max(start, a))
            x1 = self.x_of(min(end, b))
            if x1 - x0 < 5:
                mid = (x0 + x1) / 2.0
                x0, x1 = mid - 2.5, mid + 2.5
            color = QColor(str(feat.get("color") or "#2f5f9e"))
            color.setAlpha(235)
            rect = QRectF(x0, y + 3, max(5.0, x1 - x0), FEAT_LANE_H - 6)
            p.fillRect(rect, color)
            edge = QColor(color)
            edge.setAlpha(255)
            p.setPen(QPen(edge.darker(140), 1))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(rect)
            # a tick on the ruler, so a short feature still shows at Fit
            p.setPen(QPen(edge.darker(120), 2))
            p.drawLine(QPointF(rect.center().x(), MARGIN_TOP + RULER_H - 8),
                       QPointF(rect.center().x(), y + 3))
            name = str(feat.get("name") or "")
            if name and rect.width() > 28:
                p.setPen(QColor("#ffffff"))
                p.drawText(QPointF(rect.left() + 3, y + FEAT_LANE_H - 6), name)

    # alignments ----------------------------------------------------------------
    def _paint_alignments(self, p, top, a, b, w, ref, fg, mid, bg) -> None:
        t = self.track
        visible = t.reads_visible(a, b)
        n_items = t.count_in(self.chrom, a, b) if visible else 0
        # coverage — same start/bpp X mapping as reads
        band = top
        if self.show_coverage:
            cov, left, right, _cov_w = self._view_coverage(t, self.chrom)
            cov_top = band + 2
            if cov is not None and len(cov):
                raw = np.asarray(cov, dtype=float)
                raw_peak = float(np.max(raw)) if len(raw) else 0.0
                drawn = np.log1p(raw) if self.cov_log else raw
                peak = float(np.max(drawn)) if len(drawn) else 0.0
                if peak > 0:
                    p.setPen(Qt.PenStyle.NoPen)
                    p.setBrush(QColor("#a9a9a9"))
                    scale = (COV_H - 2) / peak
                    self._paint_coverage_bars(
                        p, drawn, left, right, cov_top, COV_H, scale, raw=raw
                    )
                p.setPen(fg)
                label = f"0–{int(round(raw_peak)):,}"
                if self._names_on:
                    p.drawText(
                        QRectF(6, cov_top, NAME_W - 14, COV_H),
                        Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                        label,
                    )
                else:
                    p.drawText(
                        QRectF(self._name_w() + 4, cov_top, 72, COV_H),
                        Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                        label,
                    )
            p.setPen(QPen(mid, 1))
            p.drawLine(self._name_w(), band + COV_H + 2, w, band + COV_H + 2)
            band += COV_H + 4
        if self.show_junctions:
            self._paint_junctions(p, band + 2, a, b, w, fg, mid)
            band += JUNC_H + 4
        rows_top = self._rows_top()
        area_h = self._rows_area_h()
        if not visible:
            self._content_h = 0
            win = getattr(t, "read_window", None) or 0
            msg = (f"Zoom in to {win // 1000:,} kb or less to see reads"
                   if win else "Zoom in to see reads")
            self._note(p, rows_top + 16, msg, fg)
            return
        if not self._names_on and n_items > MAX_READS_DRAWN:
            self._content_h = 0
            self._note(p, rows_top + 16, f"{n_items:,} reads here: zoom in to see them", fg)
            return
        named = self._named_in_view(a, b) if self._names_on else None
        if named is not None and not named:
            named = None
        ppb = 1.0 / self.bpp
        show_mism = (
            self.show_mismatches
            and ref is not None
            and ppb >= 0.5
            and self.layout_mode == LAYOUT_EXPANDED
        )
        show_letters = show_mism and ppb >= 8
        # Tall enough that a differing base sits inside the read, not under it.
        body_h = self._read_body_h()
        if show_letters:
            body_h = max(body_h, QFontMetrics(self._mono).height())
        row_px = (NAME_ROW_H if named is not None else body_h) + READ_GAP
        first_row = self.vscroll // row_px
        shown_rows = area_h // row_px + 2
        if named is not None:
            placed = [(it, i) for i, it in enumerate(named)]
            rows = len(placed)
        elif hasattr(t, "fetch_shown"):
            placed = t.fetch_shown(self.chrom, a, b, first_row + shown_rows)
            rows = t.row_count(self.chrom)
        else:
            placed = t.fetch(self.chrom, a, b)
            rows = t.row_count(self.chrom)
        self._content_h = rows * row_px + 6
        if not placed and hasattr(t, "loaded_near") and not t.loaded_near(self.chrom, a, b):
            self._note(p, rows_top + 16, "Loading reads…", mid)
        p.save()
        p.setClipRect(QRectF(self._name_w(), rows_top, max(1, w - self._name_w()), area_h))
        from repeatviewer.widgets.sequence_map_viewer import _draw_directed_feature

        alpha_secondary = 90
        name_rows: list[tuple[object, float]] = []
        for it, row in placed:
            if not self.read_passes(it):
                continue
            if row < first_row or row >= first_row + shown_rows:
                continue
            y_row = rows_top + row * row_px - self.vscroll
            y = y_row + max(0, (row_px - (13 if named is not None else body_h)) // 2)
            if named is not None:
                name_rows.append((it, y_row))
            strand_s = str(getattr(it, "strand", "") or "")
            col = QColor(style.read_color(it, self.color_by, chrom=self.chrom))
            if named is not None and self.color_by == style.COLOR_STRAND:
                col = QColor("#c0392b" if strand_s == "-" else "#1f6fbf")
            elif getattr(it, "secondary", False):
                col.setAlpha(alpha_secondary)
            if named is not None:
                col.setAlpha(255)
            blocks = it.ref_blocks
            x_lo = self.x_of(it.start)
            x_hi = self.x_of(it.end)
            if named is not None:
                bar_h = 13
                y = y_row + max(0, (row_px - bar_h) // 2)
                rect = QRectF(x_lo, y, max(2.0, x_hi - x_lo), bar_h)
                _draw_directed_feature(p, rect, col, forward=strand_s != "-")
                if self._is_picked(it):
                    self._paint_read_selection(p, rect)
            elif x_hi - x_lo < 2.5 or len(blocks) == 1:
                rect = QRectF(x_lo, y, max(1.0, x_hi - x_lo), body_h)
                self._read_body(p, rect, it.strand, col)
                if self._is_picked(it):
                    self._paint_read_selection(p, rect)
            else:
                last = len(blocks) - 1
                span = QRectF(x_lo, y, max(2.0, x_hi - x_lo), body_h)
                for k, (s, e, _q) in enumerate(blocks):
                    rect = QRectF(self.x_of(s), y, max(1.0, (e - s) / self.bpp), body_h)
                    strand = it.strand if (k == last and it.strand == "+") or (
                        k == 0 and it.strand == "-") else ""
                    self._read_body(p, rect, strand, col)
                if self._is_picked(it):
                    self._paint_read_selection(p, span)
            hit_h = 13 if named is not None else body_h
            self._hits.append((QRectF(x_lo, y, max(2.0, x_hi - x_lo), hit_h), it))
            glyph_h = 13 if named is not None else body_h
            if ppb >= 2:
                for pos in it.ins:
                    if a <= pos <= b:
                        x = self.x_of(pos)
                        p.fillRect(QRectF(x - 1, y - 1, 2, glyph_h + 2), QColor("#8e44ad"))
            if show_mism:
                r0, bases = ref
                for pos, base in _mismatches(it, bases, r0, a, b):
                    x = self.x_of(pos)
                    c = QColor(BASE_COLORS.get(base, "#9a9a9a"))
                    cell = QRectF(x, y, max(1.0, ppb), glyph_h)
                    if show_letters:
                        p.setPen(c)
                        p.setFont(self._mono)
                        p.drawText(cell, Qt.AlignmentFlag.AlignCenter, base)
                        p.setFont(self._font)
                    else:
                        p.fillRect(cell, c)
        p.restore()
        if named is not None:
            self._paint_name_column(p, name_rows, row_px, fg)
        if getattr(t, "truncated", False):
            self._note(p, rows_top + area_h - 6, "Very deep region: reads are down-sampled",
                       mid)

    def _read_body(self, p, rect: QRectF, strand: str, col: QColor) -> None:
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(col)
        tip = min(4.0, rect.width() / 2)
        if strand and rect.width() > 6:
            y0, y1, ym = rect.top(), rect.bottom(), rect.center().y()
            if strand == "+":
                poly = QPolygonF([QPointF(rect.left(), y0), QPointF(rect.right() - tip, y0),
                                  QPointF(rect.right(), ym), QPointF(rect.right() - tip, y1),
                                  QPointF(rect.left(), y1)])
            else:
                poly = QPolygonF([QPointF(rect.left() + tip, y0), QPointF(rect.right(), y0),
                                  QPointF(rect.right(), y1), QPointF(rect.left() + tip, y1),
                                  QPointF(rect.left(), ym)])
            p.drawPolygon(poly)
        else:
            p.drawRect(rect)

    def _paint_read_selection(self, p, rect: QRectF) -> None:
        """Square outline around the selected read glyph (not the name-column arrow)."""
        highlight = self.palette().color(self.palette().ColorRole.Highlight)
        # Pad tiny glyphs so the box stays readable around the arrow tip.
        pad = 2.0
        side = max(rect.width(), rect.height()) + pad * 2
        box = QRectF(
            rect.center().x() - side / 2,
            rect.center().y() - side / 2,
            side,
            side,
        )
        if rect.width() >= rect.height() + 4:
            box = rect.adjusted(-pad, -pad, pad, pad)
        p.save()
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.setPen(QPen(highlight, 1.5))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRoundedRect(box, 2, 2)
        p.restore()

    def _paint_name_column(self, p, rows: list[tuple[object, float]], row_px: int,
                           fg: QColor) -> None:
        """Sequence name at the left of each one-per-row read."""
        if not rows:
            return
        from repeatviewer.widgets.folder_tree_delegate import _folder_hover_background

        fm = QFontMetrics(self._font)
        pal = self.palette()
        window = pal.color(pal.ColorRole.Window)
        hover_bg = _folder_hover_background(pal)
        p.setFont(self._font)
        for it, y in rows:
            rect = QRectF(0, y, NAME_W - 1, row_px)
            self._name_hits.append((rect, it))
            strand = str(getattr(it, "strand", "") or "")
            hovering = self._same_read(self._hover, it)
            p.fillRect(rect, window)
            if hovering:
                p.save()
                p.setRenderHint(QPainter.RenderHint.Antialiasing)
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(hover_bg)
                p.drawRoundedRect(rect.adjusted(1, 1, -1, -1), 4, 4)
                p.restore()
            arrow_box = QRectF(2, y + max(0.0, (row_px - 14) / 2), 14, 14)
            p.setPen(QColor("#c0392b") if strand == "-" else QColor("#1f6fbf"))
            arrow = "←" if strand == "-" else ("→" if strand == "+" else "")
            p.drawText(arrow_box,
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignCenter, arrow)
            p.setPen(fg)
            name = str(getattr(it, "name", "") or "")
            label = fm.elidedText(name, Qt.TextElideMode.ElideMiddle, NAME_W - 26)
            p.drawText(QRectF(18, y, NAME_W - 26, row_px),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, label)
        p.setPen(QPen(self.palette().mid().color(), 1))
        p.drawLine(NAME_W - 1, 0, NAME_W - 1, self.height())

    def _same_read(self, a, b) -> bool:
        if a is b:
            return True
        if a is None or b is None:
            return False
        return (getattr(a, "name", None) == getattr(b, "name", None)
                and int(getattr(a, "start", -1)) == int(getattr(b, "start", -2)))

    def _is_picked(self, it) -> bool:
        return self._same_read(self.picked, it)

    def _set_name_hover(self, it) -> None:
        if self._same_read(self._hover, it):
            return
        self._hover = it
        self.update()

    def _note(self, p, y, text: str, color: QColor) -> None:
        p.setPen(color)
        p.drawText(QPointF(self._name_w() + 10, y), text)

    # annotations ---------------------------------------------------------------
    def _paint_annotations(self, p, top, a, b, w, fg, mid, bg) -> None:
        t = self.track
        rows_top = self._rows_top()
        area_h = self._rows_area_h()
        n = t.count_in(self.chrom, a, b)
        if n > MAX_FEATURES_DRAWN:
            cov, left, right, _cov_w = self._view_coverage(t, self.chrom)
            self._content_h = 0
            if cov is not None and len(cov) and float(np.max(cov)) > 0:
                peak = float(np.max(cov))
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(QColor(GENE_COLOR))
                self._paint_coverage_bars(
                    p, cov, left, right, rows_top + 10, 50, 50.0 / peak
                )
            self._note(p, rows_top + 80, f"{n:,} features here: zoom in to see them", fg)
            return
        placed = t.fetch(self.chrom, a, b)
        rows = t.row_count(self.chrom)
        small = rows > 40 and self.bpp > 50
        row_h = ANN_ROW_H_SMALL if small else ANN_ROW_H
        self._content_h = rows * row_h + 8
        p.save()
        p.setClipRect(QRectF(0, rows_top, w, area_h))
        fm = QFontMetrics(self._font)
        label_end: dict[int, float] = {}      # per row: where the last label ended
        for it, row in placed:
            y = rows_top + row * row_h - self.vscroll + 4
            if y > rows_top + area_h or y + row_h < rows_top:
                continue
            body_h = 10 if not small else 7
            col = QColor(GENE_COLOR if it.kind in ("gene", "mRNA", "transcript", "")
                         or len(it.blocks) > 1 else style.category_color(it.kind))
            x0, x1 = self.x_of(it.start), self.x_of(it.end)
            self._draw_feature(p, it, y, body_h, col, x0, x1)
            self._hits.append((QRectF(x0, y, max(3.0, x1 - x0), row_h - 2), it))
            label = feature_label(it) if not small else ""
            if label:
                tw = fm.horizontalAdvance(label)
                lx = min(max(x0, 2), max(2, x1 - tw)) if x1 - x0 > tw else (x0 + x1 - tw) / 2
                if lx >= label_end.get(row, -1e9) + 6:
                    p.setPen(fg)
                    p.drawText(QPointF(lx, y + body_h + 11), label)
                    label_end[row] = lx + tw
        p.restore()

    def _draw_feature(self, p, it, y, body_h, col, x0, x1) -> None:
        mid_y = y + body_h / 2
        blocks = it.blocks or [(it.start, it.end, 0)]
        thick = list(it.thick or [])
        if len(blocks) > 1 or thick:
            # intron line with strand chevrons
            p.setPen(QPen(col, 1))
            p.drawLine(QPointF(x0, mid_y), QPointF(x1, mid_y))
            if it.strand in "+-" and it.strand and x1 - x0 > _CHEVRON_GAP:
                d = 3 if it.strand == "+" else -3
                x = max(x0, 0) + _CHEVRON_GAP / 2
                while x < min(x1, self.width()):
                    p.drawLine(QPointF(x - d, mid_y - 3), QPointF(x, mid_y))
                    p.drawLine(QPointF(x - d, mid_y + 3), QPointF(x, mid_y))
                    x += _CHEVRON_GAP
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(col)
        for s, e, _q in blocks:
            xs, xe = self.x_of(s), self.x_of(e)
            if thick:
                # UTR thin, CDS full height
                p.drawRect(QRectF(xs, mid_y - body_h / 4, max(1.0, xe - xs), body_h / 2))
                for cs, ce in thick:
                    lo, hi = max(s, cs), min(e, ce)
                    if hi > lo:
                        p.drawRect(QRectF(self.x_of(lo), y, max(1.0, (hi - lo) / self.bpp),
                                          body_h))
            else:
                p.drawRect(QRectF(xs, y, max(1.0, xe - xs), body_h))
        if len(blocks) == 1 and not thick and it.strand in ("+", "-") and x1 - x0 > 14:
            p.setPen(QPen(QColor(255, 255, 255, 200), 1))
            d = 3 if it.strand == "+" else -3
            x = max(x0, 0) + 8
            while x < min(x1, self.width()) - 4:
                p.drawLine(QPointF(x - d, mid_y - 2.5), QPointF(x, mid_y))
                p.drawLine(QPointF(x - d, mid_y + 2.5), QPointF(x, mid_y))
                x += _CHEVRON_GAP

    # ------------------------------------------------------------ interaction
    def hit_at(self, pos: QPointF):
        for rect, it in reversed(self._hits):
            if rect.adjusted(-1, -1, 1, 1).contains(pos):
                return it
        return None

    def name_at(self, pos: QPointF):
        if not self._names_on or pos.x() >= self._name_w():
            return None
        for rect, it in reversed(self._name_hits):
            if rect.adjusted(0, -1, 0, 1).contains(pos):
                return it
        return None

    def _row_at(self, pos: QPointF):
        """The one-per-row sequence whose row contains this y, even beside its bar."""
        if not self._names_on:
            return None
        y = pos.y()
        for rect, it in reversed(self._name_hits):
            if rect.top() - 1 <= y <= rect.bottom() + 1:
                return it
        return None

    def reveal_read(self, it) -> None:
        """Scroll the one-per-row list so this sequence's name is in view."""
        if self.track is None:
            return
        a, b = self.view_range()
        rows = self._named_in_view(a, b)
        row_px = NAME_ROW_H + READ_GAP
        want_name = getattr(it, "name", None)
        want_start = int(getattr(it, "start", -1))
        for i, row_it in enumerate(rows):
            if row_it is it or (
                getattr(row_it, "name", None) == want_name
                and int(getattr(row_it, "start", -1)) == want_start
            ):
                self.set_vscroll(max(0, i * row_px - 4))
                return

    def _context_click(self, event) -> bool:
        """Right-click, or Control-click on macOS."""
        if event.button() == Qt.MouseButton.RightButton:
            return True
        return bool(
            event.button() == Qt.MouseButton.LeftButton
            and event.modifiers() & Qt.KeyboardModifier.ControlModifier
        )

    def _schedule_context_menu(self, pos: QPointF, global_pos) -> None:
        self._context_at = None
        if self._menu_pending:
            return
        self._menu_pending = True

        def show() -> None:
            self._menu_pending = False
            self._popup_context_menu(pos, global_pos)

        QTimer.singleShot(0, show)

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if self._context_click(event):
            self._context_at = event.position()
            event.accept()
            return
        if event.button() == Qt.MouseButton.LeftButton and self._in_ruler(event.position()):
            # A bare click selects one column only when nothing is selected.
            # With a selection already up, the click clears it and does not start another.
            self._sel_click_new = self.sel is None
            self._sel_from = self.bp_at(event.position().x())
            self.set_selection(None)
            return
        if event.button() == Qt.MouseButton.LeftButton:
            # A click off the scale drops the marked region.
            if self.sel is not None:
                self.set_selection(None)
            named = self.name_at(event.position())
            if named is not None:
                self.picked = named
                self.readPicked.emit(named)
                self.update()
                event.accept()
                return
            self._drag = (event.position(), self.start, self.vscroll)
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        self._set_ruler_hover(event.position())
        if self._sel_from is not None:
            here = self.bp_at(event.position().x())
            a, b = sorted((self._sel_from, here))
            a, b = max(0, int(a)), min(self.length, int(math.ceil(b)))
            if b > a and abs(self.x_of(b) - self.x_of(a)) >= 3:
                self.set_selection((a, b))
            return
        if self._drag is not None:
            p0, s0, v0 = self._drag
            d = event.position() - p0
            self.start = s0 - d.x() * self.bpp
            self._clamp()
            self.set_vscroll(v0 - int(d.y()))
            self.update()
            self.locusChanged.emit()
            return
        named = self.name_at(event.position())
        self._set_name_hover(named)
        if self._in_ruler(event.position()):
            QToolTip.hideText()
            return
        if named is not None or (self._names_on and event.position().x() < self._name_w()):
            QToolTip.hideText()
            return
        it = self.hit_at(event.position())
        if it is not None:
            rows = style.details(it, chrom=self.chrom)[:14]
            html = "<table>" + "".join(
                f"<tr><td><b>{k}</b>&nbsp;</td><td>{_esc(v)}</td></tr>" for k, v in rows
            ) + "</table>"
            QToolTip.showText(event.globalPosition().toPoint(), html, self)
        else:
            # No floating chrom:position tooltip while moving the mouse.
            QToolTip.hideText()

    def leaveEvent(self, event) -> None:  # noqa: N802
        self._set_name_hover(None)
        self._set_ruler_hover(None)
        super().leaveEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if self._context_at is not None and (
            event.button() in (Qt.MouseButton.RightButton, Qt.MouseButton.LeftButton)
        ):
            pos = self._context_at
            global_pos = event.globalPosition().toPoint()
            self._schedule_context_menu(pos, global_pos)
            event.accept()
            return
        if (
            self._sel_from is not None
            and self.sel is None
            and self._sel_click_new
            and event.button() == Qt.MouseButton.LeftButton
        ):
            col = int(self.bp_at(event.position().x()))
            col = max(0, min(max(0, self.length - 1), col))
            self.set_selection((col, col + 1))
        self._sel_from = None
        self._sel_click_new = False
        self._drag = None
        self.unsetCursor()
        super().mouseReleaseEvent(event)

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        x = event.position().x()
        if event.modifiers() & Qt.KeyboardModifier.ShiftModifier:
            self.zoom(0.5, x)
            return
        if self._in_ruler(event.position()):
            # ruler: zoom to the bases there
            bp = self.bp_at(x)
            half = self.width() * MIN_BPP * 2
            self.show_range(bp - half, bp + half)
            return
        it = self.name_at(event.position()) or self.hit_at(event.position())
        if it is None and self._names_on:
            it = self._row_at(event.position())
        if it is not None:
            pad = (it.end - it.start) * 0.15 + 10
            self.picked = it
            self.readPicked.emit(it)
            self.show_range(it.start - pad, it.end + pad)
            self.update()
            return
        self.zoom(2.0, x)

    def wheelEvent(self, event) -> None:  # noqa: N802
        mods = event.modifiers()
        pix = event.pixelDelta()
        ang = event.angleDelta()
        dx = pix.x() if not pix.isNull() else ang.x() / 120 * 40
        dy = pix.y() if not pix.isNull() else ang.y() / 120 * 40
        if mods & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.AltModifier
                   | Qt.KeyboardModifier.MetaModifier):
            step = ang.y() / 120.0 if ang.y() else dy / 40.0
            self.zoom(1.25 ** step, event.position().x())
        elif mods & Qt.KeyboardModifier.ShiftModifier and not dx:
            self.pan_px(dy)
        else:
            if dx:
                self.pan_px(dx)
            if dy and self._content_h > self._rows_area_h():
                self.scroll_rows(-int(dy))
            # Plain two-finger up/down must not zoom; use +/−, slider, pinch, or Cmd+scroll.
        event.accept()

    def event(self, event) -> bool:
        et = event.type()
        if et == QEvent.Type.NativeGesture:
            if event.gestureType() == Qt.NativeGestureType.ZoomNativeGesture:
                x = self.mapFromGlobal(event.globalPosition().toPoint()).x()
                self.zoom(1.0 + float(event.value()), float(x))
                return True
        if et == QEvent.Type.Gesture:
            pinch = event.gesture(Qt.GestureType.PinchGesture)
            if pinch is not None:
                f = float(pinch.scaleFactor())
                if f > 0:
                    self.zoom(f, self.mapFromGlobal(pinch.hotSpot().toPoint()).x())
                return True
        return super().event(event)

    def keyPressEvent(self, event) -> None:  # noqa: N802
        k = event.key()
        w = self.width()
        if k in (Qt.Key.Key_Plus, Qt.Key.Key_Equal):
            self.zoom(2.0)
        elif k in (Qt.Key.Key_Minus, Qt.Key.Key_Underscore):
            self.zoom(0.5)
        elif k == Qt.Key.Key_Left:
            self.pan_px(w * 0.2)
        elif k == Qt.Key.Key_Right:
            self.pan_px(-w * 0.2)
        elif k == Qt.Key.Key_Up:
            self.scroll_rows(-40)
        elif k == Qt.Key.Key_Down:
            self.scroll_rows(40)
        elif k == Qt.Key.Key_Home:
            self.show_range(0, self.length)
        elif k == Qt.Key.Key_S:
            self.go_to_selection()
        elif k == Qt.Key.Key_Escape:
            self.set_selection(None)
        else:
            super().keyPressEvent(event)

    def contextMenuEvent(self, event) -> None:  # noqa: N802
        # Mouse menus open on release, so the click that summoned them does not
        # close them at once. The menu key still comes through here.
        event.accept()
        if event.reason() != QContextMenuEvent.Reason.Keyboard and self._context_at is None:
            return
        self._schedule_context_menu(QPointF(event.pos()), event.globalPos())

    def _popup_context_menu(self, pos: QPointF, global_pos) -> None:
        named = self.name_at(pos)
        if named is None and self._names_on:
            named = self.hit_at(pos)
        if named is not None:
            self.picked = named
            self.readPicked.emit(named)
            self.update()
            self.readNameMenu.emit(named, global_pos)
            return
        menu = QMenu(self)
        a, b = self.view_range()
        locus = f"{self.chrom}:{a + 1:,}-{b:,}"
        menu.addAction(f"Copy Locus ({locus})",
                       lambda: QGuiApplication.clipboard().setText(locus))
        it = named or self.hit_at(pos)
        if it is not None:
            name = it.name or ""
            if name:
                menu.addAction(f"Copy Name ({name[:40]})",
                               lambda: QGuiApplication.clipboard().setText(name))
            span = f"{self.chrom}:{it.start + 1}-{it.end}"
            menu.addAction("Copy Position", lambda: QGuiApplication.clipboard().setText(span))
            menu.addAction("Zoom to This", lambda: self.show_range(
                it.start - (it.end - it.start) * 0.15 - 10,
                it.end + (it.end - it.start) * 0.15 + 10))
        if self.bookmarks:
            panel = self.parent()
            while panel is not None and not hasattr(panel, "go_to_bookmark"):
                panel = panel.parent()
            if panel is not None:
                sub = menu.addMenu("Bookmarks")
                for bm in self.bookmarks:
                    sub.addAction(bm.get("name", "?"),
                                  lambda b_=bm: panel.go_to_bookmark(b_))
                rm = sub.addMenu("Remove")
                for bm in self.bookmarks:
                    rm.addAction(bm.get("name", "?"),
                                 lambda b_=bm: panel.remove_bookmark(b_))
        menu.addSeparator()
        menu.addAction("Zoom In", lambda: self.zoom(2.0, float(pos.x())))
        menu.addAction("Zoom Out", lambda: self.zoom(0.5, float(pos.x())))
        menu.addAction("Whole Sequence", lambda: self.show_range(0, self.length))
        menu.exec(global_pos)


class _NoCtx:
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        return False


def _alive(obj) -> bool:
    try:
        import shiboken6

        return shiboken6.isValid(obj)
    except Exception:  # noqa: BLE001
        return True


def feature_label(it) -> str:
    """A feature's display name: its name, or for generic names (RepeatMasker GTF:
    "similarity") the repeat in ``Target``, then Name / gene_name / ID."""
    name = str(getattr(it, "name", "") or "")
    kind = str(getattr(it, "kind", "") or "")
    if name and name != kind:
        return name
    info = getattr(it, "info", None) or {}
    target = str(info.get("Target") or "").strip().strip('"')
    if target:
        first = target.split()[0].strip('"')
        return first.split(":", 1)[1] if first.lower().startswith("motif:") else first
    for key in ("Name", "gene_name", "gene_id", "ID", "transcript_id"):
        if info.get(key):
            return str(info[key])
    return name


def _esc(text) -> str:
    return (str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;"))


def release_track(track) -> None:
    """Stop a BAM track's background decoding and close its file."""
    if track is None:
        return
    try:
        track.on_tile_ready = None
    except AttributeError:
        pass
    pool = getattr(track, "_pool", None)
    if pool is not None:
        try:
            pool.shutdown(wait=False, cancel_futures=True)
        except TypeError:
            pool.shutdown(wait=False)
    bam = getattr(track, "_bam", None)
    if bam is not None:
        lock = getattr(track, "_lock", None)
        try:
            if lock is not None:
                with lock:
                    bam.close()
            else:
                bam.close()
        except Exception:  # noqa: BLE001, S110
            pass


# ---------------------------------------------------------------- run panel


def _sidebar_icon(widget: QWidget):
    """Window outline with the right column filled: show or hide the side panel."""
    from PySide6.QtGui import QIcon, QPixmap

    icon = QIcon()
    color = widget.palette().color(widget.palette().ColorRole.WindowText)
    for scale in (1, 2):
        pix = QPixmap(18 * scale, 14 * scale)
        pix.setDevicePixelRatio(scale)
        pix.fill(Qt.GlobalColor.transparent)
        p = QPainter(pix)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.setPen(QPen(color, 1.2))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRoundedRect(QRectF(1.5, 1.5, 15, 11), 2, 2)
        fill = QColor(color)
        fill.setAlphaF(0.55)
        p.fillRect(QRectF(11, 2.2, 5, 9.6), fill)
        p.drawLine(QPointF(11, 1.5), QPointF(11, 12.5))
        p.end()
        icon.addPixmap(pix)
    return icon


class _RunSidePanel(QWidget):
    """Sequences, variants, unused reads, statistics and the report for one mapping run.

    The canvas stays the genome browser. This panel only lists what the run JSON recorded.
    """

    locusRequested = Signal(str, int, int)
    optionsChanged = Signal()
    openRequested = Signal(object)
    pagesChanged = Signal()

    PAGE_KEYS = ("Sequences", "Features", "Variants", "Not used", "Report")

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        from repeatviewer.widgets.combo_boxes import NoWheelComboBox

        self._data: dict = {}
        self._chrom = ""
        self._type_boxes: dict[str, QCheckBox] = {}
        self.setAutoFillBackground(True)
        self.setBackgroundRole(self.palette().ColorRole.Window)
        col = QVBoxLayout(self)
        col.setContentsMargins(0, 0, 0, 0)
        col.setSpacing(0)

        header = QHBoxLayout()
        header.setContentsMargins(0, 0, 0, 0)
        header.setSpacing(6)
        self._summary = QLabel(self)
        self._summary.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._summary.setWordWrap(True)
        self._files_btn = QToolButton(self)
        self._files_btn.setText("Files ▾")
        self._files_btn.setAutoRaise(True)
        self._files_btn.setToolTip(
            "Consensus, placed sequences, unused sequences, tables, report, reference")
        self._files_menu = QMenu(self._files_btn)
        self._files_btn.clicked.connect(
            lambda: self._files_menu.exec(
                self._files_btn.mapToGlobal(self._files_btn.rect().bottomLeft())))
        header.addWidget(self._summary, 1)
        header.addWidget(self._files_btn, 0, Qt.AlignmentFlag.AlignTop)
        self._summary.hide()
        self._files_btn.hide()
        self._plain_alignment = False
        col.addLayout(header)

        from repeatviewer.widgets.sequence_annotation_rail import (
            _CollapsibleSection,
            _h_divider,
        )

        body = QWidget(self)
        body_col = QVBoxLayout(body)
        # Same inset, gap, and hairline as Annotation (Feature detection, …).
        body_col.setContentsMargins(10, 6, 10, 6)
        body_col.setSpacing(6)

        def section(title: str, expanded: bool = True) -> tuple[_CollapsibleSection, QFrame]:
            # One hairline above every section, including Display.
            rule = _h_divider(body)
            rule.setFrameShadow(QFrame.Shadow.Plain)
            rule.setFixedHeight(1)
            body_col.addWidget(rule)
            block = _CollapsibleSection(title, body, expanded=expanded)
            body_col.addWidget(block)
            return block, rule

        def form(*rows) -> QFormLayout:
            laid = QFormLayout()
            laid.setContentsMargins(0, 0, 0, 0)
            laid.setHorizontalSpacing(8)
            laid.setVerticalSpacing(6)
            laid.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
            for label, widget in rows:
                laid.addRow(label, widget)
            return laid

        self._layout_mode = NoWheelComboBox(body)
        for label, key in (
            ("Expanded", LAYOUT_EXPANDED),
            ("Squished", LAYOUT_SQUISHED),
            ("Collapsed", LAYOUT_COLLAPSED),
        ):
            self._layout_mode.addItem(label, key)
        self._layout_mode.setToolTip(
            "Expanded: normal read height. Squished: thinner rows. "
            "Collapsed: compact lines."
        )
        self._color_by = NoWheelComboBox(body)
        for entry in _COLOR_MENU_ENTRIES:
            if entry is None:
                continue
            mode, label = entry
            self._color_by.addItem(label, mode)
        self._group_by = NoWheelComboBox(body)
        for mode, label in _GROUP_MENU_ENTRIES:
            self._group_by.addItem(label, mode)
        self._cov_show = QCheckBox("Show coverage", body)
        self._cov_show.setChecked(True)
        self._cov_scale = NoWheelComboBox(body)
        self._cov_scale.addItem("Linear", False)
        self._cov_scale.addItem("Log", True)
        self._cov_scale.setToolTip("Coverage depth scale when Show coverage is on.")
        self._cov_show.toggled.connect(self._cov_scale.setEnabled)
        self._junc_show = QCheckBox("Show splice junctions", body)
        self._junc_show.setChecked(False)
        self._junc_show.setToolTip(
            "Draw splice junctions (CIGAR N / long gaps) as arcs above the reads."
        )
        display, _display_rule = section("Display")
        display.add_layout(form(
            ("Layout", self._layout_mode),
            ("Colour by", self._color_by),
            ("Group by", self._group_by),
            ("Coverage scale", self._cov_scale),
        ))
        display.add_widget(self._cov_show)
        display.add_widget(self._junc_show)

        self._show_ann = QCheckBox("Show on browser", body)
        self._show_ann.setChecked(True)
        self._types_host = QWidget(body)
        self._types_lay = QVBoxLayout(self._types_host)
        self._types_lay.setContentsMargins(16, 0, 0, 0)
        self._types_lay.setSpacing(2)
        annotations, self._annotations_rule = section("Annotations")
        annotations.add_widget(self._show_ann)
        annotations.add_widget(self._types_host)
        self._annotations = annotations
        body_col.addStretch(1)

        for widget in (
            self._layout_mode, self._color_by, self._group_by, self._cov_scale,
            self._show_ann, self._cov_show, self._junc_show,
        ):
            if hasattr(widget, "currentIndexChanged"):
                widget.currentIndexChanged.connect(self._emit_options)
            elif hasattr(widget, "valueChanged"):
                widget.valueChanged.connect(self._emit_options)
            elif hasattr(widget, "toggled"):
                widget.toggled.connect(self._emit_options)
        for combo in (self._layout_mode, self._color_by, self._group_by, self._cov_scale):
            combo.setSizeAdjustPolicy(
                NoWheelComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            combo.setMinimumContentsLength(8)
        scroll = QScrollArea(self)
        scroll.setWidget(body)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        col.addWidget(scroll, 1)

        self._page_labels = {key: key for key in self.PAGE_KEYS}
        self._stack = QStackedWidget(self)
        self._sequences = self._table(["Name", "Start", "End", "Strand", "Identity", "MAPQ"])
        self._features = self._table(["Name", "Type", "Start", "End", "Strand"])
        self._variants = self._table(["Position", "Ref", "Alt", "Count", "Depth"])
        self._unused = self._table(["Name", "Reason"])
        for named in (self._sequences, self._unused):
            header = named.horizontalHeader()
            header.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
            named.setColumnWidth(0, 380)
        self._report = QTextBrowser(self)
        self._report.setOpenExternalLinks(False)
        for page in (
            self._sequences, self._features, self._variants, self._unused, self._report,
        ):
            self._stack.addWidget(page)
        self._sequences.cellClicked.connect(self._on_sequence)
        self._features.cellClicked.connect(self._on_feature)
        self._variants.cellClicked.connect(self._on_variant)
        self.setMinimumWidth(280)

    def _table(self, headers: list[str]) -> QTableWidget:
        table = QTableWidget(self)
        table.setColumnCount(len(headers))
        table.setHorizontalHeaderLabels(headers)
        table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        table.setAlternatingRowColors(True)
        table.setShowGrid(False)
        table.setWordWrap(False)
        table.verticalHeader().setVisible(False)
        table.verticalHeader().setDefaultSectionSize(22)
        table.horizontalHeader().setStretchLastSection(True)
        table.horizontalHeader().setHighlightSections(False)
        table.setFrameShape(QFrame.Shape.NoFrame)
        return table

    def set_plain_alignment(self, on: bool) -> None:
        """Display, filters, and coverage for a SAM/BAM that has no mapping record."""
        self._plain_alignment = bool(on)
        show = not on
        self._annotations.setVisible(show)
        self._annotations_rule.setVisible(show)

    def set_run(self, data: dict | None) -> None:
        self._data = data or {}
        summary = (data or {}).get("summary") or {}
        placed = int(summary.get("n_placed") or 0)
        total = int(summary.get("n_sequences") or 0)
        text = f"{placed:,} of {total:,} placed" if total else f"{placed:,} placed"
        mean = summary.get("mean_identity")
        if mean is not None:
            text += f"  ·  {float(mean):.1%}"
        self._summary.setText(text)
        self._summary.setToolTip(str((data or {}).get("title") or ""))
        self._fill_files()
        self._fill_unused()
        self._fill_report()
        self.set_chrom(self._chrom)

    def options(self) -> dict:
        hidden = {name for name, box in self._type_boxes.items() if not box.isChecked()}
        return {
            "layout_mode": str(self._layout_mode.currentData() or LAYOUT_EXPANDED),
            "color_by": str(self._color_by.currentData() or style.COLOR_STRAND),
            "group_by": str(self._group_by.currentData() or style.GROUP_NONE),
            "min_identity": 0.0,
            "min_mapq": 0,
            "hide_secondary": False,
            "show_features": self._show_ann.isChecked(),
            "hidden_types": hidden,
            "show_coverage": self._cov_show.isChecked(),
            "show_junctions": self._junc_show.isChecked(),
            "cov_log": bool(self._cov_scale.currentData()),
            "cov_low": None,
            "cov_high": None,
        }

    def _emit_options(self, *_args) -> None:
        entry = self._entry()
        self._fill_sequences(entry)
        self._fill_features(entry)
        self._sync_page_titles(entry)
        self.optionsChanged.emit()

    def _read_visible(self, read: dict) -> bool:
        del read
        return True

    def set_chrom(self, chrom: str) -> None:
        self._chrom = chrom
        entry = self._entry()
        self._fill_types(entry)
        self._fill_sequences(entry)
        self._fill_features(entry)
        self._fill_variants(entry)
        self._annotations.set_expanded(bool((entry or {}).get("features")))
        self._sync_page_titles(entry)
        self._emit_options()

    def set_one_per_row_allowed(self, allowed: bool, reason: str) -> None:
        """Kept for callers; Layout is Expanded / Squished / Collapsed now."""
        del allowed, reason

    def _sync_page_titles(self, entry: dict | None) -> None:
        counts = {
            "Sequences": self._sequences.rowCount(),
            "Features": self._features.rowCount(),
            "Variants": self._variants.rowCount(),
            "Not used": self._unused.rowCount(),
        }
        for key, count in counts.items():
            self._page_labels[key] = f"{key} ({count:,})"
        self.pagesChanged.emit()

    def _entry(self) -> dict | None:
        for contig in self._data.get("contigs", []):
            if contig.get("id") == self._chrom:
                return contig
        contigs = self._data.get("contigs") or []
        return contigs[0] if len(contigs) == 1 else None

    def show_page(self, key: str) -> None:
        try:
            self._stack.setCurrentIndex(self.PAGE_KEYS.index(key))
        except ValueError:
            return

    def _fill_types(self, entry: dict | None) -> None:
        while self._types_lay.count():
            item = self._types_lay.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        self._type_boxes.clear()
        seen: dict[str, int] = {}
        for feat in (entry or {}).get("features") or []:
            kind = str(feat.get("feature_type") or feat.get("type") or "other")
            seen[kind] = seen.get(kind, 0) + 1
        self._types_host.setVisible(bool(seen))
        for kind, count in sorted(seen.items()):
            box = QCheckBox(f"{kind} ({count})", self._types_host)
            box.setChecked(True)
            box.toggled.connect(self._emit_options)
            self._type_boxes[kind] = box
            self._types_lay.addWidget(box)

    def _fill_files(self) -> None:
        self._files_menu.clear()
        folder = Path(str(self._data.get("_folder") or ""))
        details_name = str(self._data.get("details_dir") or "")
        details = folder / details_name if folder else None
        names = self._data.get("details") or {}

        def add(label: str, path: Path | None) -> None:
            if path is None:
                return
            act = self._files_menu.addAction(f"{label} — {path.name}")
            act.setEnabled(path.is_file())
            act.triggered.connect(lambda _c=False, p=path: self.openRequested.emit(p))

        for key, label in (
            ("consensus", "Consensus"),
            ("mapped", "Placed sequences"),
            ("unused", "Not used"),
            ("raw", "Aligner output"),
            ("variants", "Variants table"),
            ("features", "Features table"),
            ("report", "Report"),
        ):
            name = names.get(key)
            if not name or details is None:
                continue
            sibling = folder / str(name)
            target = sibling if sibling.is_file() else details / str(name)
            add(label, target)
        ref = Path(str(self._data.get("reference") or ""))
        if ref.name:
            self._files_menu.addSeparator()
            add("Reference", ref)

    def _fill_sequences(self, entry: dict | None) -> None:
        reads = [r for r in ((entry or {}).get("reads") or []) if self._read_visible(r)]
        reads.sort(key=lambda read: (int(read.get("start", 0)), int(read.get("end", 0))))
        self._sequences.setRowCount(len(reads))
        for row, read in enumerate(reads):
            ident = read.get("identity")
            values = [
                str(read.get("name", "")),
                f"{int(read.get('start', 0)) + 1:,}",
                f"{int(read.get('end', 0)):,}",
                str(read.get("strand", "")),
                f"{float(ident):.1%}" if ident is not None else "",
                "" if read.get("mapq", -1) in (-1, None) else str(read.get("mapq")),
            ]
            for col, text in enumerate(values):
                self._sequences.setItem(row, col, QTableWidgetItem(text))
            self._sequences.item(row, 0).setData(Qt.ItemDataRole.UserRole, read)

    def _fill_features(self, entry: dict | None) -> None:
        feats = []
        for feat in (entry or {}).get("features") or []:
            kind = str(feat.get("feature_type") or feat.get("type") or "other")
            box = self._type_boxes.get(kind)
            if box is not None and not box.isChecked():
                continue
            if not self._show_ann.isChecked():
                continue
            feats.append(feat)
        self._features.setRowCount(len(feats))
        for row, feat in enumerate(feats):
            values = [
                str(feat.get("name", "")),
                str(feat.get("feature_type") or feat.get("type") or ""),
                f"{int(feat.get('start', 0)):,}",
                f"{int(feat.get('end', 0)):,}",
                str(feat.get("strand", "")),
            ]
            for col, text in enumerate(values):
                self._features.setItem(row, col, QTableWidgetItem(text))
            self._features.item(row, 0).setData(Qt.ItemDataRole.UserRole, feat)

    def _on_feature(self, row: int, _col: int) -> None:
        item = self._features.item(row, 0)
        feat = item.data(Qt.ItemDataRole.UserRole) if item is not None else None
        if not feat or not self._chrom:
            return
        start = max(0, int(feat.get("start", 1)) - 1)
        end = max(start + 1, int(feat.get("end", start + 1)))
        self.locusRequested.emit(self._chrom, start, end)

    def _fill_variants(self, entry: dict | None) -> None:
        variants = list((entry or {}).get("variants") or [])
        self._variants.setRowCount(len(variants))
        for row, item in enumerate(variants):
            pos, ref, alt, count, depth = (list(item) + ["", "", "", "", ""])[:5]
            for col, text in enumerate((f"{int(pos):,}", str(ref), str(alt), str(count), str(depth))):
                self._variants.setItem(row, col, QTableWidgetItem(text))
            self._variants.item(row, 0).setData(Qt.ItemDataRole.UserRole, int(pos))

    def _fill_unused(self) -> None:
        unused = list(self._data.get("unused") or [])
        self._unused.setRowCount(len(unused))
        for row, item in enumerate(unused):
            name, reason = (list(item) + ["", ""])[:2]
            self._unused.setItem(row, 0, QTableWidgetItem(str(name)))
            self._unused.setItem(row, 1, QTableWidgetItem(str(reason)))

    def _fill_report(self) -> None:
        self._report.clear()
        details = self._data.get("details") or {}
        name = details.get("report")
        folder = self._data.get("_folder")
        if not name or not folder:
            return
        path = Path(folder) / str(self._data.get("details_dir") or "") / str(name)
        if path.is_file():
            self._report.setHtml(path.read_text(encoding="utf-8", errors="replace"))

    def _on_sequence(self, row: int, _col: int) -> None:
        item = self._sequences.item(row, 0)
        read = item.data(Qt.ItemDataRole.UserRole) if item is not None else None
        if not read or not self._chrom:
            return
        start = int(read.get("start", 0))
        end = max(start + 1, int(read.get("end", start + 1)))
        self.locusRequested.emit(self._chrom, start, end)

    def _on_variant(self, row: int, _col: int) -> None:
        item = self._variants.item(row, 0)
        pos = item.data(Qt.ItemDataRole.UserRole) if item is not None else None
        if pos is None or not self._chrom:
            return
        pos = int(pos)
        self.locusRequested.emit(self._chrom, max(0, pos - 20), pos + 20)


# ---------------------------------------------------------------- scaffold bar


class _ScaffoldStrip(QWidget):
    """Thin position bar: the whole scaffold, with a mark for the window."""

    HEIGHT = 16

    def __init__(self, canvas: LightBrowserCanvas, parent=None) -> None:
        super().__init__(parent)
        self.canvas = canvas
        self.setFixedHeight(self.HEIGHT)
        self.setMouseTracking(True)
        self.setToolTip("Whole scaffold. Click or drag to move the view.")
        self._drag = False
        canvas.locusChanged.connect(self.update)

    def _track(self) -> QRectF:
        return QRectF(6, 4, max(1.0, self.width() - 12), 8)

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.fillRect(self.rect(), self.palette().color(self.palette().ColorRole.Base))
        bar = self._track()
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor("#d4dee8"))
        p.drawRect(bar)
        canvas = self.canvas
        if canvas.chrom and canvas.length > 0:
            # The bar is the whole scaffold. The red span is the window, as in
            # the genome browser overview: zoom in and it shrinks.
            a, b = canvas.view_range()
            start_f = max(0.0, min(1.0, a / canvas.length))
            end_f = max(start_f, min(1.0, b / canvas.length))
            left = bar.left() + bar.width() * start_f
            right = bar.left() + bar.width() * end_f
            # Square frame, a little taller than the bar, so a tiny window
            # stays a box instead of a dot.
            side = bar.height() + 4
            span_w = max(0.0, right - left)
            width = max(side, span_w)
            x = left + span_w / 2.0 - width / 2.0
            x = min(bar.right() - width, max(bar.left(), x))
            mark = QRectF(x, bar.top() - 2, width, side).adjusted(0.5, 0.5, -0.5, -0.5)
            p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
            p.setPen(QPen(QColor("#e53935"), 1))
            p.setBrush(QColor(255, 255, 255, 90))
            p.drawRect(mark)
        p.end()

    def _jump(self, x: float) -> None:
        canvas = self.canvas
        if not canvas.chrom or canvas.length <= 1:
            return
        bar = self._track()
        bp = (x - bar.left()) / max(1.0, bar.width()) * canvas.length
        a, b = canvas.view_range()
        span = max(1, b - a)
        start = bp - span / 2.0
        canvas.show_range(start, start + span)

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self.canvas.set_selection(None)
            self._drag = True
            self._jump(event.position().x())
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if self._drag and event.buttons() & Qt.MouseButton.LeftButton:
            self._jump(event.position().x())
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        self._drag = False
        super().mouseReleaseEvent(event)


# ---------------------------------------------------------------- panel


class LightGenomeBrowserPanel(QWidget):
    """Preview page body for a SAM / BAM / GFF file: toolbar, banner, canvas, and a
    Text view of the file (Browser / Text toggle)."""

    openRequested = Signal(object)
    fileWritten = Signal(object)
    blastRequested = Signal(object)
    generateCoverageRequested = Signal(object)  # Path → build .rvcov.npz

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.path: Path | None = None
        self._token = 0
        self._result: _LoadResult | None = None
        self._load_progress_on = False
        self._load_progress_path: Path | None = None
        self._load_progress_timer = QTimer(self)
        self._load_progress_timer.setSingleShot(True)
        self._load_progress_timer.timeout.connect(self._show_load_progress)
        self._run: dict | None = None
        self._run_ref: Path | None = None
        self._text_loaded_for: Path | None = None

        self.canvas = LightBrowserCanvas(self)
        self.canvas.readNameMenu.connect(self._on_read_name_menu)
        self._vbar = QScrollBar(Qt.Orientation.Vertical, self)
        self._vbar.valueChanged.connect(self.canvas.set_vscroll)
        self.canvas.rowsChanged.connect(self._sync_vbar)
        self.canvas.locusChanged.connect(self._sync_locus)

        # navigation bar, as in the genome browser (GenomeFlow): back / forward, the
        # genome (reference) chip, the sequence chip, the position box, then − + and the
        # zoom slider, Select (the region dragged on the ruler), and Fit
        from repeatviewer.widgets.light_browser_nav import BreadcrumbBar
        from repeatviewer.widgets.preview_controls import (
            ANNOTATION_FIT_BUTTON_RENDERED_WIDTH,
            ANNOTATION_SELECT_BUTTON_RENDERED_WIDTH,
            ANNOTATION_ZOOM_SLIDER_WIDTH,
            ZOOM_SLIDER_HEIGHT,
            apply_annotation_zoom_slider_style,
            make_overview_zoom_button,
        )
        from repeatviewer.widgets.spin_boxes import NoWheelSlider

        # Nav buttons a bit taller than the shared 18 px overview controls, still
        # under the locus QLineEdit / chips (~24–26 px).
        _NAV_BTN_H = 22
        _NAV_BTN_INNER = 20

        def _style_nav_button(btn: QPushButton) -> None:
            btn.setFixedHeight(_NAV_BTN_H)
            btn.setStyleSheet(
                "QPushButton { border: 1px solid palette(mid); border-radius: 4px;"
                f" background: palette(base); padding: 0px 8px; margin: 0px;"
                f" min-height: {_NAV_BTN_INNER}px; max-height: {_NAV_BTN_INNER}px; }}"
                "QPushButton:hover { background: palette(midlight); }"
                "QPushButton:pressed { background: palette(mid); }"
                "QPushButton:checked { background: palette(midlight); }"
                "QPushButton:disabled { color: palette(mid); }"
            )

        self._history: list[tuple[str, int, int]] = []
        self._history_at = -1
        self._history_timer = QTimer(self)
        self._history_timer.setSingleShot(True)
        self._history_timer.setInterval(700)        # a view that rests is a Back stop
        self._history_timer.timeout.connect(self._record_history)
        self._crumbs = BreadcrumbBar(self)
        self._crumbs.backRequested.connect(lambda: self._step_history(-1))
        self._crumbs.forwardRequested.connect(lambda: self._step_history(1))
        self._crumbs.chromRequested.connect(self._on_chrom_chosen)
        self._crumbs.positionEntered.connect(self._on_locus_entered)
        self._ref_menu = QMenu(self)
        self._ref_menu.aboutToShow.connect(self._fill_ref_menu)
        self._crumbs.set_genome("No reference", self._ref_menu)

        nav = QWidget(self)
        nl = QHBoxLayout(nav)
        nl.setContentsMargins(6, 0, 0, 0)
        nl.setSpacing(6)
        zoom_out = make_overview_zoom_button("−", nav)
        zoom_out.setToolTip("Zoom out (−; pinch or Cmd + scroll)")
        zoom_out.clicked.connect(lambda: self.canvas.zoom(0.66))
        zoom_in = make_overview_zoom_button("+", nav)
        zoom_in.setToolTip("Zoom in (+; pinch or Cmd + scroll)")
        zoom_in.clicked.connect(lambda: self.canvas.zoom(1.5))
        for btn in (zoom_out, zoom_in):
            btn.setAutoRepeat(True)
            btn.setAutoRepeatDelay(300)
            btn.setAutoRepeatInterval(60)
            btn.setFixedSize(_NAV_BTN_H, _NAV_BTN_H)
        self._zoom_slider = NoWheelSlider(Qt.Orientation.Horizontal, nav)
        self._zoom_slider.setObjectName("annotationZoomSlider")
        self._zoom_slider.setRange(0, ZOOM_STEPS)
        self._zoom_slider.setFixedWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
        self._zoom_slider.setFixedHeight(ZOOM_SLIDER_HEIGHT)
        self._zoom_slider.setToolTip("Zoom: left shows the whole sequence, right single bases")
        apply_annotation_zoom_slider_style(self._zoom_slider)
        self._zoom_slider.valueChanged.connect(self._on_zoom_slid)
        select = QPushButton("Select", nav)
        select.setFixedWidth(ANNOTATION_SELECT_BUTTON_RENDERED_WIDTH)
        select.setToolTip("Zoom to the region dragged on the ruler (S)\n"
                          "Ruler: double-click zooms to the bases there")
        select.clicked.connect(self.canvas.go_to_selection)
        fit = QPushButton("Fit", nav)
        fit.setFixedWidth(ANNOTATION_FIT_BUTTON_RENDERED_WIDTH)
        fit.setToolTip("Show the whole sequence (Home)")
        fit.clicked.connect(lambda: self.canvas.show_range(0, self.canvas.length))
        saved = str(_settings().value(_SETTINGS_COLOR, style.COLOR_STRAND) or "")
        if saved not in _COLOR_MODE_LABELS:
            saved = style.COLOR_STRAND
        self.canvas.color_by = saved
        for b in (select, fit):
            _style_nav_button(b)
            b.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        for w in (zoom_out, zoom_in, self._zoom_slider, select, fit):
            w.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
            nl.addWidget(w)
        self.canvas.selectionChanged.connect(
            lambda: select.setEnabled(self.canvas.sel is not None))
        select.setEnabled(False)

        self._view_key = "Browser"
        self._view_btn = QPushButton("Browser ▾", self)
        _style_nav_button(self._view_btn)
        self._view_btn.setToolTip("Browser, text, or a mapping list")
        self._view_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._view_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._view_menu = QMenu(self._view_btn)
        self._view_btn.clicked.connect(self._popup_view_menu)
        self._fit_view_button("Browser")

        bar_w = QWidget(self)
        bar_w.setObjectName("lightBrowserNavBar")
        bar_w.setAutoFillBackground(True)
        bar_w.setBackgroundRole(bar_w.palette().ColorRole.Window)
        bar = QHBoxLayout(bar_w)
        bar.setContentsMargins(6, 4, 6, 4)
        bar.setSpacing(4)
        self._browser_tools = [self._crumbs, nav]
        bar.addWidget(self._crumbs)
        bar.addWidget(nav)
        bar.addStretch(1)
        bar.addWidget(self._view_btn)
        self._side_btn = QToolButton(bar_w)
        self._side_btn.setCheckable(True)
        self._side_btn.setAutoRaise(True)
        self._side_btn.setIcon(_sidebar_icon(self))
        self._side_btn.setToolTip("Show or hide the mapping panel")
        self._side_btn.setFixedHeight(_NAV_BTN_H)
        self._side_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._side_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._side_btn.toggled.connect(self._set_side_visible)
        self._side_btn.hide()
        bar.addWidget(self._side_btn)
        self._nav_bar = bar_w
        self._bookmarks: list[dict] = []

        # banner (reference mismatch / errors / info)
        self._banner = QFrame(self)
        self._banner.setObjectName("lightBrowserBanner")
        self._banner.setStyleSheet(
            "QFrame#lightBrowserBanner { background: rgba(255, 196, 0, 45);"
            " border-bottom: 1px solid palette(mid); }")
        bl = QHBoxLayout(self._banner)
        bl.setContentsMargins(8, 4, 8, 4)
        bl.setSpacing(6)
        self._banner_label = QLabel(self._banner)
        self._banner_label.setWordWrap(False)
        bl.addWidget(self._banner_label, 1)
        self._banner_btn = QPushButton("Choose Reference…", self._banner)
        self._banner_btn.clicked.connect(self._choose_reference_file)
        bl.addWidget(self._banner_btn)
        self._banner_close = QToolButton(self._banner)
        self._banner_close.setText("✕")
        self._banner_close.setAutoRaise(True)
        self._banner_close.setCursor(Qt.CursorShape.PointingHandCursor)
        self._banner_close.setToolTip("Dismiss this message")
        self._banner_close.setFixedSize(22, 22)
        self._banner_close.clicked.connect(self._dismiss_banner)
        bl.addWidget(self._banner_close, 0, Qt.AlignmentFlag.AlignVCenter)
        self._banner_dismissed_text = ""
        self._banner.hide()

        self._status = QLabel(self)
        self._status.setContentsMargins(8, 2, 8, 2)
        self._status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._status.hide()

        canvas_row = QWidget(self)
        cl = QVBoxLayout(canvas_row)
        cl.setContentsMargins(0, 0, 0, 0)
        cl.setSpacing(0)
        self._scaffold_strip = _ScaffoldStrip(self.canvas, canvas_row)
        cl.addWidget(self._scaffold_strip)
        canvas_line = QWidget(canvas_row)
        ll = QHBoxLayout(canvas_line)
        ll.setContentsMargins(0, 0, 0, 0)
        ll.setSpacing(0)
        ll.addWidget(self.canvas, 1)
        ll.addWidget(self._vbar)
        cl.addWidget(canvas_line, 1)
        self._main_stack = QStackedWidget(self)
        self._main_stack.addWidget(canvas_row)
        self._run_panel = _RunSidePanel(self)
        self._main_stack.addWidget(self._run_panel._stack)
        self._run_panel.pagesChanged.connect(self._sync_view_menu)
        self.canvas.readPicked.connect(self._on_read_picked)
        color_idx = self._run_panel._color_by.findData(self.canvas.color_by)
        if color_idx >= 0:
            self._run_panel._color_by.setCurrentIndex(color_idx)
        self._run_panel.locusRequested.connect(self._go)
        self._run_panel.optionsChanged.connect(self._apply_run_options)
        self._run_panel.openRequested.connect(self.openRequested)
        self._run_panel.hide()
        self._view_split = QSplitter(Qt.Orientation.Horizontal, self)
        self._view_split.addWidget(self._main_stack)
        self._view_split.addWidget(self._run_panel)
        self._view_split.setStretchFactor(0, 1)
        self._view_split.setStretchFactor(1, 0)
        self._view_split.setCollapsible(0, False)
        self._message = QLabel(self)
        self._message.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._message.setWordWrap(True)
        self._browser_stack = QStackedWidget(self)
        self._browser_stack.addWidget(self._message)
        self._browser_stack.addWidget(self._view_split)

        self._text = QPlainTextEdit(self)
        self._text.setReadOnly(True)
        self._text.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        mono = QFont("Menlo")
        mono.setStyleHint(QFont.StyleHint.Monospace)
        self._text.setFont(mono)

        self._stack = QStackedWidget(self)
        self._stack.addWidget(self._browser_stack)
        self._stack.addWidget(self._text)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(self._nav_bar)
        lay.addWidget(self._banner)
        lay.addWidget(self._stack, 1)
        lay.addWidget(self._status)

        box = self.canvas.box
        # the page may be deleted while a BAM is open: stop its decoding threads then
        self.destroyed.connect(lambda *_a, b=box: release_track(b.get("track")))

        copy = QAction(self)
        copy.setShortcut(QKeySequence.StandardKey.Copy)
        copy.setShortcutContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        copy.triggered.connect(self._copy_locus)
        self.addAction(copy)

    # ------------------------------------------------------------ loading
    def load(self, path: Path, *, reference: Path | None = None, forced: bool = False) -> None:
        path = Path(path).expanduser().resolve()
        new_file = path != self.path
        self.path = path
        self._token += 1
        if new_file:
            self._text_loaded_for = None
            self.canvas.set_track(None, None)
            self._result = None
        kind = "alignment" if _is_alignment(path) else "annotation"
        text_first = kind in _text_mode_kinds()
        self._set_mode(text=text_first, remember=False)
        self._read_run(path)
        if not forced:
            reference = self._run_ref or remembered_reference(path)
        self._message.setText(f"Opening {path.name}…")
        self._browser_stack.setCurrentIndex(0)
        self._banner_dismissed_text = ""
        self._banner.hide()
        self._set_browser_status("")
        self._arm_load_progress(path)
        project = None if forced else _load_project()
        loader = _Loader(self._token, path, project, reference, forced)
        loader.loaded.connect(self._on_loaded)
        _start_thread(loader)

    def _open_progress_host(self):
        widget = self.window()
        if hasattr(widget, "begin_preview_open_progress"):
            return widget
        return None

    def _arm_load_progress(self, path: Path) -> None:
        """Show the status-bar bar if indexing is still going after a short wait."""
        self._end_load_progress()
        try:
            big = path.stat().st_size >= _LOAD_PROGRESS_MIN_BYTES
        except OSError:
            big = False
        if not big or self._open_progress_host() is None:
            return
        self._load_progress_path = path
        self._load_progress_timer.start(_LOAD_PROGRESS_DELAY_MS)

    def _show_load_progress(self) -> None:
        path = self._load_progress_path
        host = self._open_progress_host()
        if path is None or host is None:
            return
        if not host.begin_preview_open_progress(path.name, total=0):
            return
        if self._load_progress_path is not path:
            # Indexing finished while the bar was being shown.
            end = getattr(host, "end_preview_open_progress", None)
            if callable(end):
                end()
            return
        self._load_progress_on = True
        detail = "Indexing…" if _is_alignment(path) else "Reading…"
        host.update_preview_open_progress(0, detail)

    def _end_load_progress(self) -> None:
        self._load_progress_timer.stop()
        self._load_progress_path = None
        if not self._load_progress_on:
            return
        self._load_progress_on = False
        host = self._open_progress_host()
        end = getattr(host, "end_preview_open_progress", None)
        if callable(end):
            end()

    def _on_loaded(self, token, res: _LoadResult) -> None:
        if token != self._token:
            release_track(res.track)
            return
        self._end_load_progress()
        if res.error:
            self._message.setText(f"Cannot show {self.path.name} as a genome browser:\n\n"
                                  f"{res.error}\n\nUse Text to see the file.")
            self._browser_stack.setCurrentIndex(0)
            return
        keep_view = None
        if self._result is not None and self._result.track is not None and self.canvas.chrom:
            keep_view = (self.canvas.chrom, *self.canvas.view_range())
        self._result = res
        self._history, self._history_at = [], -1
        self._load_bookmarks()
        self._apply_loaded_result(keep_view=keep_view)
        for w in res.track.warnings if res.track is not None else []:
            self._set_browser_status(w)

    def _go(self, chrom: str, a: int, b: int) -> None:
        length = (self._result.lengths.get(chrom, 1) if self._result else 1)
        self.canvas.set_chrom(chrom, length, a, b)
        self._apply_track_grouping()
        self._apply_run_chrom(chrom)
        self._sync_one_per_row_choice()
        self._sync_vbar()

    def _read_run(self, path: Path) -> None:
        from repeatviewer.backends.map_reference import find_mapping_document, load_rvmap

        self._run = None
        self._run_ref = None
        self._mapped_hit_cache = {}
        self._query_seq_cache = {}
        self.canvas.run_features = []
        self.canvas.one_per_row = False
        found = find_mapping_document(path)
        if found is not None:
            try:
                data = load_rvmap(found)
            except (OSError, ValueError, Exception):  # noqa: BLE001 - not a run record
                data = None
            if data and "contigs" in data:
                self._run = dict(data)
                self._run["_folder"] = str(found.parent)
                ref = str(data.get("reference") or "")
                if ref and Path(ref).is_file():
                    self._run_ref = Path(ref)
        self._sync_run_panel()

    def _sync_run_panel(self) -> None:
        alignment = self.path is not None and _is_alignment(self.path)
        show = self._run is not None or alignment
        self._side_btn.setVisible(show)
        self._side_btn.setToolTip(
            "Show or hide the mapping panel" if self._run is not None
            else "Show or hide display and filters")
        self._run_panel.set_plain_alignment(alignment and self._run is None)
        self._run_panel.set_run(self._run)
        if not show:
            self._run_panel.hide()
            self._side_btn.blockSignals(True)
            self._side_btn.setChecked(False)
            self._side_btn.blockSignals(False)
            if self._view_key not in ("Browser", "Text"):
                self._set_view("Browser")
            self._sync_view_menu()
            return
        if not self._side_btn.isChecked():
            self._side_btn.blockSignals(True)
            self._side_btn.setChecked(True)
            self._side_btn.blockSignals(False)
        self._set_side_visible(True)
        self._sync_view_menu()

    def _set_side_visible(self, on: bool) -> None:
        # The mapping panel belongs to the genome browser. Report, Text, and the
        # lists use the main area on their own.
        alignment = self.path is not None and _is_alignment(self.path)
        visible = bool(on) and (self._run is not None or alignment) and self._view_key == "Browser"
        self._run_panel.setVisible(visible)
        if not visible:
            return
        sizes = self._view_split.sizes()
        total = sum(sizes) or self.width()
        if len(sizes) == 2 and sizes[1] < 220:
            self._view_split.setSizes([max(0, total - 280), 280])

    def _apply_run_options(self) -> None:
        opts = self._run_panel.options()
        canvas = self.canvas
        canvas.layout_mode = str(opts.get("layout_mode") or LAYOUT_EXPANDED)
        canvas.group_by = str(opts.get("group_by") or style.GROUP_NONE)
        canvas.min_identity = float(opts["min_identity"])
        canvas.min_mapq = int(opts["min_mapq"])
        canvas.hide_secondary = bool(opts["hide_secondary"])
        canvas.show_features = bool(opts["show_features"])
        canvas.show_coverage = bool(opts["show_coverage"])
        canvas.show_junctions = bool(opts["show_junctions"])
        canvas.hidden_feature_types = set(opts["hidden_types"])
        canvas.show_mismatches = True
        canvas.cov_log = bool(opts["cov_log"])
        canvas.cov_low = opts["cov_low"]
        canvas.cov_high = opts["cov_high"]
        mode = str(opts["color_by"])
        if mode in _COLOR_MODE_LABELS and canvas.color_by != mode:
            canvas.color_by = mode
            _settings().setValue(_SETTINGS_COLOR, mode)
        self._apply_track_grouping()
        self._sync_name_column()
        canvas.update()

    def _apply_track_grouping(self) -> None:
        """Push Display → Group by onto the open alignment track."""
        track = self.canvas.track
        if track is None or not hasattr(track, "set_grouping"):
            return
        mode = str(self.canvas.group_by or style.GROUP_NONE)
        fn = style.group_function(
            mode,
            read_groups=getattr(track, "read_groups", None) or {},
            chrom=self.canvas.chrom or "",
        )
        key = None if fn is None else mode
        track.set_grouping(key, fn)

    def _sync_name_column(self) -> None:
        """Turn the name column on with One per row, and keep a fitted view fitted."""
        canvas = self.canvas
        on = canvas._names_wanted()
        if on != canvas._names_on:
            a, b = canvas.view_range()
            whole = a <= 1 and b >= canvas.length - 1
            canvas._names_on = on
            if whole and canvas.length > 1:
                canvas.show_range(0, canvas.length)
            else:
                canvas._clamp()
        canvas.update()

    def _on_read_name_menu(self, it, global_pos) -> None:
        menu = self._read_name_menu(it)
        if menu.actions():
            menu.exec(global_pos)

    def _read_name_menu(self, it) -> QMenu:
        """Right-click on a one-per-row sequence: copy, extract, zoom, show, hide, BLAST."""
        menu = QMenu(self)
        name = str(getattr(it, "name", "") or "")
        if name:
            menu.addAction("Copy name", lambda n=name: QGuiApplication.clipboard().setText(n))
        bases = self._read_bases(it)
        text = self._fasta_text(it, bases) if bases else ""
        copy = menu.addAction(
            "Copy sequence (FASTA)", lambda t=text: QGuiApplication.clipboard().setText(t))
        copy.setEnabled(bool(bases))
        extract = menu.addAction(
            "Extract sequence to a file…", lambda item=it: self._extract_read(item))
        extract.setEnabled(bool(bases) and self.path is not None)
        menu.addSeparator()
        menu.addAction("Zoom to this sequence", lambda item=it: self._zoom_to_read(item))
        if self._run is not None:
            menu.addAction(
                "Show in the Sequences table", lambda item=it: self._show_read_in_sequences(item))
        menu.addAction("Hide sequence", lambda item=it: self.canvas.hide_read(item))
        hidden = self.canvas.user_hidden_count()
        if hidden:
            menu.addAction(
                f"Show hidden sequences ({hidden:,})", self.canvas.show_hidden_reads)
        menu.addSeparator()
        blast = menu.addAction("BLAST sequence…", lambda item=it: self._blast_read(item))
        blast.setEnabled(bool(bases) and self.path is not None)
        return menu

    def _fasta_text(self, it, bases: str) -> str:
        title = (f"{getattr(it, 'name', '')} {self.canvas.chrom}:{int(it.start) + 1:,}-"
                 f"{int(it.end):,} {getattr(it, 'strand', '')}").strip()
        body = "\n".join(bases[i:i + 70] for i in range(0, len(bases), 70))
        return f">{title}\n{body}\n"

    def _read_bases(self, it) -> str:
        """Query bases for copy, extract, and BLAST.

        The alignment record is used when it still holds the sequence. A PAF
        file does not, and a large SAM drops sequences past a memory budget, so
        the run's placed-sequence file, the original query file, and the stored
        alignment row are used in that order.
        """
        seq = str(getattr(it, "seq", "") or "")
        if _has_bases(seq):
            return seq
        if not self._run:
            return ""
        placed = self._sequence_from_mapped(it)
        if placed:
            return placed
        name = str(getattr(it, "name", "") or "")
        query = self._query_sequence(name)
        if query:
            if str(getattr(it, "strand", "+") or "+") == "-":
                from repeatviewer.backends.map_reference import revcomp

                return revcomp(query)
            return query
        return self._sequence_from_run_row(it)

    def _mapped_fasta_path(self) -> Path | None:
        run = self._run or {}
        name = str((run.get("details") or {}).get("mapped") or "")
        if not name:
            return None
        folder = Path(str(run.get("_folder") or ""))
        details = str(run.get("details_dir") or "")
        path = (folder / details / name) if details else (folder / name)
        return path if path.is_file() else None

    def _sequence_from_mapped(self, it) -> str:
        """This placement from the run's ``_mapped.fa`` (whole sequence, reference strand)."""
        path = self._mapped_fasta_path()
        if path is None:
            return ""
        chrom = str(self.canvas.chrom or "")
        start0 = int(getattr(it, "start", -1))
        end = int(getattr(it, "end", -1))
        strand = str(getattr(it, "strand", "+") or "+")
        name = str(getattr(it, "name", "") or "")
        cache = getattr(self, "_mapped_hit_cache", None)
        if cache is None:
            cache = {}
            self._mapped_hit_cache = cache
        key = (str(path), chrom, start0, end, strand, name)
        if key in cache:
            return cache[key]
        exact = loose = ""
        taking = ""
        parts: list[str] = []

        def finish() -> None:
            nonlocal exact, loose
            if taking == "exact" and not exact:
                exact = "".join(parts)
            elif taking == "loose" and not loose:
                loose = "".join(parts)

        with path.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    finish()
                    if exact:
                        break
                    taking = ""
                    parts = []
                    header = line[1:].strip()
                    rid, _, loc = header.rpartition(" mapped to ")
                    if not rid or (name and _placement_name(rid, strand) != name):
                        continue
                    if loc == f"{chrom}:{start0 + 1}-{end} ({strand})":
                        taking = "exact"
                    elif loc.startswith(f"{chrom}:{start0 + 1}-") and loc.endswith(f"({strand})"):
                        taking = "loose"
                elif taking:
                    parts.append(line.strip())
            finish()
        seq = exact or loose
        cache[key] = seq
        return seq

    def _query_sequence(self, name: str) -> str:
        """The original query record, when the run still points at that file."""
        if not name or not self._run:
            return ""
        cache = getattr(self, "_query_seq_cache", None)
        if cache is None:
            cache = {}
            self._query_seq_cache = cache
        if name in cache:
            return cache[name]
        from repeatviewer.backends.map_reference import iter_sequences

        found = ""
        for raw in self._run.get("queries") or []:
            path = Path(str(raw))
            if not path.is_file():
                continue
            try:
                for qname, seq in iter_sequences(path):
                    if qname == name and _has_bases(seq):
                        found = seq
                        break
            except OSError:
                continue
            if found:
                break
        cache[name] = found
        return found

    def _sequence_from_run_row(self, it) -> str:
        """Aligned bases stored on the run, including bases inserted against the reference."""
        name = str(getattr(it, "name", "") or "")
        start = int(getattr(it, "start", -1))
        chrom = self.canvas.chrom
        matches: list[dict] = []
        for contig in (self._run or {}).get("contigs") or []:
            if contig.get("id") != chrom:
                continue
            for read in contig.get("reads") or []:
                if str(read.get("name", "")) != name:
                    continue
                matches.append(read)
        chosen = None
        for read in matches:
            if int(read.get("start", -2)) == start:
                chosen = read
                break
        if chosen is None and len(matches) == 1:
            chosen = matches[0]
        if chosen is None:
            return ""
        return _bases_from_row(chosen)

    def _extract_read(self, it) -> None:
        """Selected sequence into <alignment>_selected.fa next to the alignment."""
        from PySide6.QtWidgets import QMessageBox

        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        bases = self._read_bases(it)
        if self.path is None or not bases:
            return
        proposed = self.path.with_name(f"{self.path.stem}_selected.fa")
        target = confirm_existing_file_path(
            self, proposed, title="Extract sequence", required_suffix=".fa")
        if target is None:
            return
        try:
            target.write_text(self._fasta_text(it, bases), encoding="utf-8")
        except OSError as exc:
            QMessageBox.warning(self, "Extract sequence", f"Could not write {target.name}: {exc}")
            return
        self.fileWritten.emit(target)
        self.openRequested.emit(target)

    def _blast_read(self, it) -> None:
        """Write the sequence as a query next to the alignment and open BLAST."""
        bases = self._read_bases(it)
        if self.path is None or not bases:
            return
        target = self.path.with_name(f"{self.path.stem}_blast_query.fa")
        try:
            target.write_text(self._fasta_text(it, bases), encoding="utf-8")
        except OSError:
            return
        self.fileWritten.emit(target)
        self.blastRequested.emit(target)

    def _zoom_to_read(self, it) -> None:
        pad = max(10, int((int(it.end) - int(it.start)) * 0.15))
        self.canvas.show_range(int(it.start) - pad, int(it.end) + pad)
        self.canvas.reveal_read(it)

    def _show_read_in_sequences(self, it) -> None:
        if self._run is None:
            return
        self._set_view("Sequences")
        name = str(getattr(it, "name", "") or "")
        start = int(getattr(it, "start", -1))
        table = self._run_panel._sequences
        for row in range(table.rowCount()):
            item = table.item(row, 0)
            read = item.data(Qt.ItemDataRole.UserRole) if item is not None else None
            if not isinstance(read, dict):
                continue
            if str(read.get("name", "")) == name and int(read.get("start", -2)) == start:
                table.selectRow(row)
                table.scrollToItem(item)
                return

    def _sync_one_per_row_choice(self) -> None:
        allowed, reason = self.canvas.one_per_row_availability()
        self._run_panel.set_one_per_row_allowed(allowed, reason)

    def _apply_run_chrom(self, chrom: str) -> None:
        if self._run is None:
            self.canvas.run_features = []
            return
        self._run_panel.set_chrom(chrom)
        entry = next((c for c in self._run.get("contigs", []) if c.get("id") == chrom), None)
        self.canvas.run_features = list((entry or {}).get("features") or [])
        self.canvas.update()

    # ------------------------------------------------------------ reference
    _REF_MISMATCH_MSG = "Please check if the reference genome is correct."

    def _reference_mismatch(self, res: _LoadResult | None = None) -> bool:
        """True when chromosome names / count / lengths do not fit the reference."""
        res = res if res is not None else self._result
        if res is None or res.track is None:
            return False
        ref = res.reference
        if ref is None:
            # Tried candidates but none fit — warn and do not draw.
            return bool(res.tried) and not res.forced
        return not ref.ok

    def _apply_loaded_result(
        self, *, keep_view: tuple[str, int, int] | None = None
    ) -> None:
        """Show the track when the reference fits; otherwise warn and hide it."""
        res = self._result
        if res is None or res.track is None:
            return
        mismatch = self._reference_mismatch(res)
        self._crumbs.set_chroms([(n, res.lengths.get(n, 0)) for n in res.order])
        if mismatch:
            self.canvas.set_track(None, None)
            self._message.setText(self._REF_MISMATCH_MSG)
            self._browser_stack.setCurrentIndex(0)
            self._update_reference_ui()
            self._show_banner(self._REF_MISMATCH_MSG)
            return
        self.canvas.set_track(res.track, res.reference)
        self._browser_stack.setCurrentIndex(1)
        if keep_view is not None and keep_view[0] in res.lengths:
            chrom, a, b = keep_view
            self._go(chrom, a, b)
        elif res.order:
            # Open at Fit (whole sequence) for alignments and annotations alike.
            chrom = next((n for n in res.order if n in res.first), res.order[0])
            length = res.lengths.get(chrom, 1)
            self._go(chrom, 0, length)
        self._update_reference_ui()
        self._maybe_offer_coverage_cache()

    def _maybe_offer_coverage_cache(self) -> None:
        """Remind that large BAMs need ``.rvcov.npz`` for overview coverage."""
        track = self.canvas.track
        if track is None or not hasattr(track, "needs_coverage_cache"):
            return
        if not track.needs_coverage_cache():
            return
        if self.path is None:
            return
        # Keep a reference-mismatch banner if it is already showing.
        if self._banner.isVisible():
            return
        self._show_coverage_banner()

    def _show_coverage_banner(self) -> None:
        text = (
            "Large BAM: zoom to 100 kb or less for reads. "
            "Generate a coverage cache for overview bars at wide zoom."
        )
        if text == self._banner_dismissed_text:
            return
        self._banner_label.setText(text)
        self._banner_btn.setText("Generate Coverage…")
        try:
            self._banner_btn.clicked.disconnect()
        except RuntimeError:
            pass
        self._banner_btn.clicked.connect(self._request_coverage_cache)
        self._banner_btn.show()
        self._banner.show()

    def _request_coverage_cache(self) -> None:
        if self.path is not None:
            self.generateCoverageRequested.emit(self.path)

    def reload_coverage_cache(self) -> bool:
        """Reload ``.rvcov.npz`` into the open track after Generate coverage."""
        track = self.canvas.track
        reload = getattr(track, "reload_coverage_cache", None)
        if not callable(reload):
            return False
        ok = bool(reload())
        if ok:
            self.canvas.update()
            self._banner.hide()
            self._set_browser_status("Coverage cache loaded — overview bars available.")
        return ok

    def _update_reference_ui(self) -> None:
        res = self._result
        if res is None or res.track is None:
            return
        ref = res.reference
        if ref is not None:
            self._crumbs.set_genome(ref.path.name, self._ref_menu)
            self._crumbs.genome.setToolTip(
                f"Reference: {ref.label}\n{ref.path}\n{ref.matched:,} of {ref.total:,} "
                f"sequences of this file found in it\nClick to choose another reference")
            if self._reference_mismatch(res):
                self._show_banner(self._REF_MISMATCH_MSG)
            else:
                self._banner.hide()
            return
        self._crumbs.set_genome("No reference", self._ref_menu)
        self._crumbs.genome.setToolTip("No reference FASTA: reads and features are drawn "
                                       "without reference bases.\nClick to choose one")
        if res.forced:
            self._banner.hide()
            return
        if self._reference_mismatch(res):
            self._show_banner(self._REF_MISMATCH_MSG)
        else:
            self._banner.hide()

    def _show_banner(self, text: str) -> None:
        message = str(text or "").strip()
        if not message:
            self._banner.hide()
            return
        if message == self._banner_dismissed_text:
            self._banner.hide()
            return
        self._banner_label.setText(message)
        self._banner_btn.setText("Choose Reference…")
        try:
            self._banner_btn.clicked.disconnect()
        except RuntimeError:
            pass
        self._banner_btn.clicked.connect(self._choose_reference_file)
        self._banner_btn.show()
        self._banner.show()

    def _dismiss_banner(self) -> None:
        self._banner_dismissed_text = self._banner_label.text().strip()
        self._banner.hide()

    def _fill_ref_menu(self) -> None:
        m = self._ref_menu
        m.clear()
        res = self._result
        current = res.reference.path if res is not None and res.reference else None
        project = _load_project()
        cands = refmod.reference_candidates(project, remembered_reference(self.path)
                                            if self.path else None)
        for label, path in cands:
            act = m.addAction(f"{label}: {path.name}")
            act.setToolTip(str(path))
            act.setCheckable(True)
            act.setChecked(current is not None and Path(current).resolve() == path.resolve())
            act.triggered.connect(lambda _c=False, p=path: self._use_reference(p))
        if cands:
            m.addSeparator()
        m.addAction("Choose FASTA…", self._choose_reference_file)
        none = m.addAction("No Reference", lambda: self._use_reference(None))
        none.setCheckable(True)
        none.setChecked(current is None)

    def _choose_reference_file(self) -> None:
        start = str(self.path.parent) if self.path else ""
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose reference FASTA", start,
            "FASTA (*.fa *.fasta *.fna *.fas *.fsa *.faa *.ffn *.frn *.mfa);;All files (*)")
        if path:
            self._use_reference(Path(path))

    def _use_reference(self, fasta: Path | None) -> None:
        if self.path is None:
            return
        remember_reference(self.path, fasta)
        res = self._result
        if res is None or res.track is None:
            self.load(self.path, reference=fasta, forced=True)
            return
        # keep the opened track: only the reference changes (fast, still off the GUI
        # thread for a new .fai)
        self._token += 1
        token = self._token
        track = res.track
        self._set_browser_status("Reading reference…" if fasta else "")

        job = _RefJob(fasta, track)

        def finish(match, err) -> None:
            if not _alive(self) or token != self._token:
                return
            if err:
                self._set_browser_status(f"Cannot read the reference: {err}")
                return
            self._set_browser_status("")
            res.reference = match
            res.forced = True
            res.tried = [match] if match is not None else []
            res.lengths = refmod.sequence_lengths(track, match)
            res.order = refmod.genome_order(track, match)
            keep = (
                (self.canvas.chrom, *self.canvas.view_range())
                if self.canvas.chrom
                else None
            )
            self._banner_dismissed_text = ""
            self._apply_loaded_result(keep_view=keep)

        job.done.connect(finish)
        _start_thread(job)

    # ------------------------------------------------------------ toolbar
    def _on_chrom_chosen(self, chrom: str) -> None:
        if not chrom or self._result is None:
            return
        # Always Fit the whole scaffold; keep the previous zoom only hurts when
        # jumping between chromosomes of very different lengths.
        length = self._result.lengths.get(chrom, 1)
        self._go(chrom, 0, length)

    def _on_locus_entered(self, text: str) -> None:
        text = text.strip()
        if not text or self._result is None:
            return
        lengths = self._result.lengths
        bare = text.replace(",", "").replace(" ", "")
        if self.canvas.chrom and bare and all(c.isdigit() or c in "-–." for c in bare):
            # the position box shows start-end on the current sequence
            text = f"{self.canvas.chrom}:{text}"
        got = gt.parse_locus(text)
        if got is not None:
            name, s1, e1 = got
            chrom = name if name in lengths else _match_name(lengths, name)
            if chrom is not None:
                if s1 == e1:
                    s1, e1 = max(1, s1 - 50), e1 + 50
                self._go(chrom, s1 - 1, e1)
                self.canvas.setFocus()
                return
        chrom = text if text in lengths else _match_name(lengths, text)
        if chrom is not None:
            self._go(chrom, 0, lengths[chrom])
            self.canvas.setFocus()
            return
        track = self._result.track
        try:
            hits = track.find(text, limit=5)
        except Exception:  # noqa: BLE001
            hits = []
        if hits:
            chrom, it = hits[0]
            pad = (it.end - it.start) * 0.15 + 20
            self._go(chrom, int(it.start - pad), int(it.end + pad))
            self.canvas.setFocus()
        else:
            self._set_browser_status(f"Not found: {text}")

    def _set_browser_status(self, text: str) -> None:
        self._status.setText(text)
        self._status.setVisible(bool(text.strip()))

    def _sync_locus(self) -> None:
        if not self.canvas.chrom:
            return
        a, b = self.canvas.view_range()
        self._crumbs.set_where(self.canvas.chrom, a, b, "")
        # Packed view has no span footer. One per row shows it in the name column.
        self._set_browser_status("")
        self._sync_zoom_slider()
        self._history_timer.start()

    def _sync_vbar(self) -> None:
        total = self.canvas.content_height()
        page = self.canvas._rows_area_h()
        self._vbar.blockSignals(True)
        self._vbar.setRange(0, max(0, total - page))
        self._vbar.setPageStep(max(1, page))
        step = (
            NAME_ROW_H if self.canvas._names_on else self.canvas._read_body_h()
        ) + READ_GAP
        self._vbar.setSingleStep(step)
        self._vbar.setValue(self.canvas.vscroll)
        self._vbar.blockSignals(False)
        self._vbar.setVisible(total > page)

    # ------------------------------------------------------------ zoom slider
    def _zoom_bounds(self) -> tuple[float, float]:
        # bp per pixel at Fit (left end of the slider) and at single bases (right end)
        fit = max(MIN_BPP, self.canvas.length / max(1, self.canvas.width()))
        return fit, MIN_BPP

    def _sync_zoom_slider(self) -> None:
        fit, top = self._zoom_bounds()
        if fit <= top:
            v = 0.0
        else:
            v = math.log(fit / max(top, min(fit, self.canvas.bpp))) / math.log(fit / top)
        self._zoom_slider.blockSignals(True)
        self._zoom_slider.setValue(round(v * ZOOM_STEPS))
        self._zoom_slider.blockSignals(False)

    def _on_zoom_slid(self, value: int) -> None:
        fit, top = self._zoom_bounds()
        target = fit * (top / fit) ** (value / ZOOM_STEPS) if fit > top else fit
        self.canvas.zoom(self.canvas.bpp / max(1e-12, target))

    # ------------------------------------------------------------ back / forward
    def _record_history(self) -> None:
        if not self.canvas.chrom:
            return
        a, b = self.canvas.view_range()
        here = (self.canvas.chrom, a, b)
        if 0 <= self._history_at < len(self._history) and self._history[self._history_at] == here:
            return
        del self._history[self._history_at + 1:]
        self._history.append(here)
        if len(self._history) > 100:
            del self._history[0]
        self._history_at = len(self._history) - 1
        self._crumbs.set_history(self._history_at > 0, False)

    def _step_history(self, step: int) -> None:
        at = self._history_at + step
        if not 0 <= at < len(self._history):
            return
        self._history_at = at
        chrom, a, b = self._history[at]
        self._go(chrom, a, b)
        self._history_timer.stop()          # moving through history is not a new stop
        self._crumbs.set_history(at > 0, at < len(self._history) - 1)

    # ------------------------------------------------------------ bookmarks (legacy)
    def _bookmark_key(self) -> str:
        import hashlib

        digest = hashlib.md5(str(self.path).encode()).hexdigest()[:16]
        return f"LightGenomeBrowser/bookmarks/{digest}"

    def _load_bookmarks(self) -> None:
        try:
            data = json.loads(str(_settings().value(self._bookmark_key(), "[]") or "[]"))
        except (ValueError, TypeError):
            data = []
        self._bookmarks = [d for d in data if isinstance(d, dict)
                           and d.get("path") == str(self.path)]
        self.canvas.bookmarks = self._bookmarks
        self.canvas.update()

    def _save_bookmarks(self) -> None:
        _settings().setValue(self._bookmark_key(), json.dumps(self._bookmarks))
        self.canvas.bookmarks = self._bookmarks
        self.canvas.update()

    def go_to_bookmark(self, bm: dict) -> None:
        self._go(bm["chrom"], bm["start"], bm["end"])

    def remove_bookmark(self, bm: dict) -> None:
        if bm in self._bookmarks:
            self._bookmarks.remove(bm)
            self._save_bookmarks()

    def _copy_locus(self) -> None:
        if self.canvas.chrom:
            a, b = self.canvas.view_range()
            QGuiApplication.clipboard().setText(f"{self.canvas.chrom}:{a + 1}-{b}")

    # ------------------------------------------------------------ view menu
    def _set_mode(self, *, text: bool, remember: bool) -> None:
        self._set_view("Text" if text else "Browser", remember=remember)

    def _popup_view_menu(self) -> None:
        self._sync_view_menu()
        self._view_menu.exec(self._view_btn.mapToGlobal(self._view_btn.rect().bottomLeft()))

    def _sync_view_menu(self) -> None:
        menu = self._view_menu
        menu.clear()
        entries = [("Browser", "Browser"), ("Text", "Text")]
        if self._run is not None:
            entries.append(None)
            for key in _RunSidePanel.PAGE_KEYS:
                entries.append((self._run_panel._page_labels.get(key, key), key))
        for entry in entries:
            if entry is None:
                menu.addSeparator()
                continue
            label, key = entry
            act = menu.addAction(label)
            act.setCheckable(True)
            act.setChecked(key == self._view_key)
            act.triggered.connect(lambda _c=False, k=key: self._set_view(k, remember=True))
        if self._view_key in _RunSidePanel.PAGE_KEYS:
            self._fit_view_button(self._run_panel._page_labels.get(self._view_key, self._view_key))

    def _fit_view_button(self, label: str) -> None:
        text = f"{label} ▾"
        self._view_btn.setText(text)
        width = self._view_btn.fontMetrics().horizontalAdvance(text) + 28
        self._view_btn.setFixedWidth(max(96, width))

    def _set_view(self, key: str, *, remember: bool = False) -> None:
        pages = _RunSidePanel.PAGE_KEYS
        if key in pages and self._run is None:
            key = "Browser"
        self._view_key = key
        text = key == "Text"
        if self.path is not None and remember and key in ("Browser", "Text"):
            _set_text_mode("alignment" if _is_alignment(self.path) else "annotation", text)
        for w in self._browser_tools:
            w.setEnabled(not text)
        if text:
            self._fill_text()
            self._stack.setCurrentIndex(1)
            self._fit_view_button("Text")
            self._set_side_visible(self._side_btn.isChecked())
            return
        self._stack.setCurrentIndex(0)
        if key == "Browser":
            self._main_stack.setCurrentIndex(0)
            self._fit_view_button("Browser")
            self._set_side_visible(self._side_btn.isChecked())
            self.canvas.setFocus()
            return
        self._run_panel.show_page(key)
        self._main_stack.setCurrentIndex(1)
        self._fit_view_button(self._run_panel._page_labels.get(key, key))
        self._set_side_visible(False)

    def _on_read_picked(self, it) -> None:
        """Keep the Sequences row in step with the name clicked on the left."""
        if self._run is None:
            return
        name = str(getattr(it, "name", "") or "")
        start = int(getattr(it, "start", -1))
        table = self._run_panel._sequences
        for row in range(table.rowCount()):
            item = table.item(row, 0)
            read = item.data(Qt.ItemDataRole.UserRole) if item is not None else None
            if not isinstance(read, dict):
                continue
            if str(read.get("name", "")) == name and int(read.get("start", -2)) == start:
                table.selectRow(row)
                table.scrollToItem(item)
                return

    def _fill_text(self) -> None:
        if self.path is None or self._text_loaded_for == self.path:
            return
        self._text_loaded_for = self.path
        self._text.setPlainText(read_text_head(self.path))

    def shutdown(self) -> None:
        """Called when the page goes away: stop background reads."""
        self._token += 1
        release_track(self.canvas.track)
        self.canvas.track = None
        self.canvas.box["track"] = None


def _match_name(lengths: dict[str, int], name: str) -> str | None:
    return refmod.NameMatcher(lengths)(name)


TEXT_MAX_CHARS = 500_000
BAM_TEXT_RECORDS = 2_000


def read_text_head(path: Path) -> str:
    """The start of the file as text; a BAM / CRAM is shown as SAM (header + reads)."""
    path = Path(path)
    if _is_binary(path):
        if not gt.pysam_available():
            return "Install pysam to see BAM / CRAM records as text."
        import pysam

        try:
            with pysam.AlignmentFile(str(path), "rc" if path.suffix.lower() == ".cram"
                                     else "rb") as bam:
                lines = [str(bam.header).rstrip("\n")]
                for n, a in enumerate(bam.fetch(until_eof=True)):
                    if n >= BAM_TEXT_RECORDS:
                        lines.append(f"\n… (first {BAM_TEXT_RECORDS:,} records shown)")
                        break
                    lines.append(a.to_string())
            return "\n".join(lines)
        except (OSError, ValueError) as exc:
            return f"Cannot read {path.name}: {exc}"
    try:
        with gt._open_text(path) as fh:
            text = fh.read(TEXT_MAX_CHARS + 1)
    except OSError as exc:
        return f"Cannot read {path.name}: {exc}"
    if len(text) > TEXT_MAX_CHARS:
        text = text[:TEXT_MAX_CHARS] + (f"\n\n… (first {TEXT_MAX_CHARS:,} characters shown; "
                                        "the file is larger)")
    return text


__all__ = ["LightGenomeBrowserPanel", "LightBrowserCanvas", "is_light_browser_file",
           "LIGHT_BROWSER_SUFFIXES", "read_text_head", "release_track"]
