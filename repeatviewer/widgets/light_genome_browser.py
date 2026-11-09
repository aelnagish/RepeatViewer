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
    QPainter,
    QPen,
    QPolygonF,
)
from PySide6.QtWidgets import (
    QButtonGroup,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QMenu,
    QPlainTextEdit,
    QPushButton,
    QScrollBar,
    QSizePolicy,
    QStackedWidget,
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

# Colour-by menu for alignments (full labels — no combo elision).
# Separators match GenomeFlow's grouping; Tag / RG options stay out of this light UI.
_COLOR_MENU_ENTRIES: tuple[tuple[str, str] | None, ...] = (
    (style.COLOR_STRAND, "Strand"),
    None,
    (style.COLOR_FIRST_STRAND, "First-of-pair strand"),
    (style.COLOR_INSERT, "Insert size and pair orientation"),
    None,
    (style.COLOR_MAPQ, "Mapping quality"),
    (style.COLOR_IDENTITY, "Identity"),
    None,
    (style.COLOR_NONE, "None"),
)
_COLOR_MODE_LABELS = {
    mode: label for entry in _COLOR_MENU_ENTRIES if entry for mode, label in (entry,)
}


def is_light_browser_file(path: Path) -> bool:
    return Path(path).name.lower().endswith(LIGHT_BROWSER_SUFFIXES)


def _is_alignment(path: Path) -> bool:
    return Path(path).name.lower().endswith(ALIGNMENT_SUFFIXES)


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
COV_H = 44
READ_H = 9
READ_GAP = 2
ANN_ROW_H = 26
ANN_ROW_H_SMALL = 12
MARGIN_TOP = 2
MAX_READS_DRAWN = 30_000
MAX_FEATURES_DRAWN = 25_000
ZOOM_STEPS = 100                # zoom slider: Fit (0) to single bases
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
        self.reference: refmod.ReferenceMatch | None = None
        self._ref_find = None
        self._seq_cache: tuple[str, int, str] | None = None   # (chrom, start, bases)
        self._drag: tuple[QPointF, float, int] | None = None
        self.sel: tuple[int, int] | None = None        # region dragged on the ruler
        self._sel_from: float | None = None
        self.bookmarks: list[dict] = []
        self._pinch_anchor: float | None = None
        self._content_h = 0
        self._hits: list[tuple[QRectF, object]] = []
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
        w = max(1, self.width())
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
        return pos.y() < MARGIN_TOP + RULER_H

    def view_range(self) -> tuple[int, int]:
        a = max(0, int(math.floor(self.start)))
        b = min(self.length, int(math.ceil(self.start + self.bpp * self.width())))
        return a, max(a + 1, b)

    def _clamp(self) -> None:
        w = max(1, self.width())
        self.bpp = min(max(self.bpp, MIN_BPP), max(MIN_BPP, self.length / w))
        span = self.bpp * w
        self.start = min(max(0.0, self.start), max(0.0, self.length - span))

    def x_of(self, bp: float) -> float:
        return (bp - self.start) / self.bpp

    def bp_at(self, x: float) -> float:
        return self.start + x * self.bpp

    def zoom(self, factor: float, x: float | None = None) -> None:
        """``factor`` > 1 zooms in, around pixel ``x`` (default: the middle)."""
        if x is None:
            x = self.width() / 2.0
        anchor = self.bp_at(x)
        self.bpp /= factor
        self._clamp()
        self.start = anchor - x * self.bpp
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
        old_w = max(1, event.oldSize().width())
        if event.oldSize().width() > 0:
            # keep the same bases in view
            span = self.bpp * old_w
            self.bpp = span / max(1, self.width())
        self._clamp()
        super().resizeEvent(event)
        self.rowsChanged.emit()

    # ------------------------------------------------------------ geometry
    def _seq_lane_h(self) -> int:
        return SEQ_H if self.reference is not None and self._ref_name() else 0

    def _rows_top(self) -> int:
        top = MARGIN_TOP + RULER_H + self._seq_lane_h()
        if self.track is not None and self.track.kind == gt.KIND_ALIGNMENT:
            top += COV_H + 4
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
        if self.track is None or not self.chrom:
            p.end()
            return
        w = self.width()
        a, b = self.view_range()
        y = MARGIN_TOP
        self._paint_ruler(p, y, w, fg, mid)
        self._paint_marks(p, y)
        y += RULER_H
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
            x0, x1 = self.x_of(self.sel[0]), self.x_of(self.sel[1])
            p.fillRect(QRectF(x0, MARGIN_TOP, max(1.0, x1 - x0), self.height()),
                       QColor(61, 123, 224, 28))
            p.setPen(QPen(QColor("#3d7be0"), 1))
            p.drawLine(QPointF(x0, MARGIN_TOP), QPointF(x0, self.height()))
            p.drawLine(QPointF(x1, MARGIN_TOP), QPointF(x1, self.height()))
        p.end()
        if self._content_h != before:
            # the scroll bar follows the rows now known (after the paint, not in it)
            if self.vscroll > max(0, self._content_h - self._rows_area_h()):
                QTimer.singleShot(0, lambda: self.set_vscroll(self.vscroll))
            QTimer.singleShot(0, self.rowsChanged.emit)

    def _paint_ruler(self, p: QPainter, y: int, w: int, fg: QColor, mid: QColor) -> None:
        span = self.bpp * w
        step = _nice_step(span, w)
        a = int(self.start // step) * step
        p.setPen(QPen(mid, 1))
        base = y + RULER_H - 6
        p.drawLine(0, base, w, base)
        fm = QFontMetrics(self._font)
        p.setPen(fg)
        v = a
        while v <= self.start + span:
            x = self.x_of(v)
            if x >= -1:
                p.drawLine(QPointF(x, base - 4), QPointF(x, base))
                label = _fmt_bp(v, step) if v else "0"
                tw = fm.horizontalAdvance(label)
                p.drawText(QPointF(max(2.0, x - tw / 2), base - 7), label)
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

    # alignments ----------------------------------------------------------------
    def _paint_alignments(self, p, top, a, b, w, ref, fg, mid, bg) -> None:
        t = self.track
        visible = t.reads_visible(a, b)
        n_items = t.count_in(self.chrom, a, b) if visible else 0
        # coverage
        cov = t.coverage(self.chrom, a, b, w)
        cov_top = top + 2
        if cov is not None and len(cov):
            peak = float(np.max(cov)) if len(cov) else 0.0
            if peak > 0:
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(QColor("#a9a9a9"))
                xs = np.arange(len(cov), dtype=float)
                scale = (COV_H - 12) / peak
                for x, v in zip(xs, cov):
                    if v > 0:
                        hh = max(1.0, v * scale)
                        p.drawRect(QRectF(x, cov_top + COV_H - hh, 1.0, hh))
            p.setPen(fg)
            p.drawText(QPointF(4, cov_top + 10), f"0–{int(round(peak)):,}")
        p.setPen(QPen(mid, 1))
        p.drawLine(0, top + COV_H + 2, w, top + COV_H + 2)
        rows_top = self._rows_top()
        area_h = self._rows_area_h()
        row_px = READ_H + READ_GAP
        if not visible:
            self._content_h = 0
            win = getattr(t, "read_window", None) or 0
            msg = (f"Zoom in to {win // 1000:,} kb or less to see reads"
                   if win else "Zoom in to see reads")
            self._note(p, rows_top + 16, msg, fg)
            return
        if n_items > MAX_READS_DRAWN:
            self._content_h = 0
            self._note(p, rows_top + 16, f"{n_items:,} reads here: zoom in to see them", fg)
            return
        first_row = self.vscroll // row_px
        shown_rows = area_h // row_px + 2
        if hasattr(t, "fetch_shown"):
            placed = t.fetch_shown(self.chrom, a, b, first_row + shown_rows)
        else:
            placed = t.fetch(self.chrom, a, b)
        rows = t.row_count(self.chrom)
        self._content_h = rows * row_px + 6
        if not placed and hasattr(t, "loaded_near") and not t.loaded_near(self.chrom, a, b):
            self._note(p, rows_top + 16, "Loading reads…", mid)
        p.save()
        p.setClipRect(QRectF(0, rows_top, w, area_h))
        ppb = 1.0 / self.bpp
        show_mism = ref is not None and ppb >= 0.5
        alpha_secondary = 90
        line_pen = QPen(QColor(120, 125, 135, 170), 1)
        for it, row in placed:
            if row < first_row or row >= first_row + shown_rows:
                continue
            y = rows_top + row * row_px - self.vscroll + 2
            col = QColor(style.read_color(it, self.color_by, chrom=self.chrom))
            if getattr(it, "secondary", False):
                col.setAlpha(alpha_secondary)
            blocks = it.ref_blocks
            x_lo = self.x_of(it.start)
            x_hi = self.x_of(it.end)
            if x_hi - x_lo < 2.5 or len(blocks) == 1:
                rect = QRectF(x_lo, y, max(1.0, x_hi - x_lo), READ_H)
                self._read_body(p, rect, it.strand, col)
            else:
                p.setPen(line_pen)
                p.drawLine(QPointF(x_lo, y + READ_H / 2), QPointF(x_hi, y + READ_H / 2))
                last = len(blocks) - 1
                for k, (s, e, _q) in enumerate(blocks):
                    rect = QRectF(self.x_of(s), y, max(1.0, (e - s) / self.bpp), READ_H)
                    strand = it.strand if (k == last and it.strand == "+") or (
                        k == 0 and it.strand == "-") else ""
                    self._read_body(p, rect, strand, col)
            self._hits.append((QRectF(x_lo, y, max(2.0, x_hi - x_lo), READ_H), it))
            if ppb >= 2:
                for pos in it.ins:
                    if a <= pos <= b:
                        x = self.x_of(pos)
                        p.fillRect(QRectF(x - 1, y - 1, 2, READ_H + 2), QColor("#8e44ad"))
            if show_mism:
                r0, bases = ref
                for pos, base in _mismatches(it, bases, r0, a, b):
                    x = self.x_of(pos)
                    c = QColor(BASE_COLORS.get(base, "#9a9a9a"))
                    if ppb >= 8:
                        p.setPen(c)
                        p.setFont(self._mono)
                        p.drawText(QPointF(x + ppb / 2 - 3, y + READ_H - 1), base)
                        p.setFont(self._font)
                    else:
                        p.fillRect(QRectF(x, y, max(1.0, ppb), READ_H), c)
        p.restore()
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

    def _note(self, p, y, text: str, color: QColor) -> None:
        p.setPen(color)
        p.drawText(QPointF(10, y), text)

    # annotations ---------------------------------------------------------------
    def _paint_annotations(self, p, top, a, b, w, fg, mid, bg) -> None:
        t = self.track
        rows_top = self._rows_top()
        area_h = self._rows_area_h()
        n = t.count_in(self.chrom, a, b)
        if n > MAX_FEATURES_DRAWN:
            cov = t.coverage(self.chrom, a, b, w)
            self._content_h = 0
            if cov is not None and len(cov) and float(np.max(cov)) > 0:
                peak = float(np.max(cov))
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(QColor(GENE_COLOR))
                for x, v in enumerate(cov):
                    if v > 0:
                        hh = max(1.0, v / peak * 50)
                        p.drawRect(QRectF(x, rows_top + 60 - hh, 1.0, hh))
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

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton and self._in_ruler(event.position()):
            # drag on the ruler: mark a region (Select zooms there)
            self._sel_from = self.bp_at(event.position().x())
            self.set_selection(None)
            return
        if event.button() == Qt.MouseButton.LeftButton:
            self._drag = (event.position(), self.start, self.vscroll)
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
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
        it = self.hit_at(event.position())
        if it is not None:
            rows = style.details(it, chrom=self.chrom)[:14]
            html = "<table>" + "".join(
                f"<tr><td><b>{k}</b>&nbsp;</td><td>{_esc(v)}</td></tr>" for k, v in rows
            ) + "</table>"
            QToolTip.showText(event.globalPosition().toPoint(), html, self)
        else:
            bp = int(self.bp_at(event.position().x())) + 1
            if 0 < bp <= self.length:
                self.setToolTip("")
                QToolTip.showText(event.globalPosition().toPoint(), f"{self.chrom}:{bp:,}",
                                  self)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        self._sel_from = None
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
        it = self.hit_at(event.position())
        if it is not None and (it.end - it.start) / self.bpp < self.width() * 0.5:
            pad = (it.end - it.start) * 0.15 + 10
            self.show_range(it.start - pad, it.end + pad)
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
            if dy:
                if self._content_h > self._rows_area_h():
                    self.scroll_rows(-int(dy))
                elif not dx:
                    # nothing to scroll: a vertical wheel zooms, like most genome browsers
                    self.zoom(1.0015 ** dy, event.position().x())
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
        menu = QMenu(self)
        a, b = self.view_range()
        locus = f"{self.chrom}:{a + 1:,}-{b:,}"
        menu.addAction(f"Copy Locus ({locus})",
                       lambda: QGuiApplication.clipboard().setText(locus))
        it = self.hit_at(QPointF(event.pos()))
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
        menu.addAction("Zoom In", lambda: self.zoom(2.0, float(event.pos().x())))
        menu.addAction("Zoom Out", lambda: self.zoom(0.5, float(event.pos().x())))
        menu.addAction("Whole Sequence", lambda: self.show_range(0, self.length))
        menu.exec(event.globalPos())


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


# ---------------------------------------------------------------- panel


class LightGenomeBrowserPanel(QWidget):
    """Preview page body for a SAM / BAM / GFF file: toolbar, banner, canvas, and a
    Text view of the file (Browser / Text toggle)."""

    def __init__(self, parent=None) -> None:
        super().__init__(parent)
        self.path: Path | None = None
        self._token = 0
        self._result: _LoadResult | None = None
        self._text_loaded_for: Path | None = None

        self.canvas = LightBrowserCanvas(self)
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
            OVERVIEW_SELECT_FIT_HEIGHT,
            ZOOM_SLIDER_HEIGHT,
            apply_annotation_zoom_slider_style,
            make_overview_zoom_button,
            style_compact_button,
        )
        from repeatviewer.widgets.spin_boxes import NoWheelSlider

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
        color_btn = QPushButton("Color", nav)
        color_btn.setFixedWidth(
            max(
                ANNOTATION_FIT_BUTTON_RENDERED_WIDTH,
                color_btn.fontMetrics().horizontalAdvance("Color") + 20,
            )
        )
        color_btn.setToolTip("Colour reads by strand, quality, identity, …")
        color_btn.clicked.connect(self._show_color_menu)
        self._color_btn = color_btn
        saved = str(_settings().value(_SETTINGS_COLOR, style.COLOR_STRAND) or "")
        if saved not in _COLOR_MODE_LABELS:
            saved = style.COLOR_STRAND
        self.canvas.color_by = saved
        for b in (select, fit, color_btn):
            b.setFixedHeight(OVERVIEW_SELECT_FIT_HEIGHT)
            style_compact_button(b)
            b.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        for w in (zoom_out, zoom_in, self._zoom_slider, select, fit, color_btn):
            w.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
            nl.addWidget(w)
        self.canvas.selectionChanged.connect(
            lambda: select.setEnabled(self.canvas.sel is not None))
        select.setEnabled(False)

        self._browser_btn = QPushButton("Browser", self)
        self._text_btn = QPushButton("Text", self)
        for btn, tip in (
            (self._browser_btn, "Show the genome browser view"),
            (self._text_btn, "Show the file as text"),
        ):
            btn.setCheckable(True)
            btn.setFixedHeight(OVERVIEW_SELECT_FIT_HEIGHT)
            btn.setFixedWidth(
                max(
                    ANNOTATION_FIT_BUTTON_RENDERED_WIDTH,
                    btn.fontMetrics().horizontalAdvance(btn.text()) + 20,
                )
            )
            style_compact_button(btn)
            btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
            btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
            btn.setToolTip(tip)
        group = QButtonGroup(self)
        group.setExclusive(True)
        group.addButton(self._browser_btn)
        group.addButton(self._text_btn)
        self._browser_btn.setChecked(True)
        self._browser_btn.clicked.connect(lambda: self._set_mode(text=False, remember=True))
        self._text_btn.clicked.connect(lambda: self._set_mode(text=True, remember=True))

        bar = QHBoxLayout()
        bar.setContentsMargins(6, 4, 6, 4)
        bar.setSpacing(4)
        self._browser_tools = [self._crumbs, nav]
        bar.addWidget(self._crumbs)
        bar.addWidget(nav)
        bar.addStretch(1)
        bar.addWidget(self._browser_btn)
        bar.addWidget(self._text_btn)
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

        canvas_row = QWidget(self)
        cl = QHBoxLayout(canvas_row)
        cl.setContentsMargins(0, 0, 0, 0)
        cl.setSpacing(0)
        cl.addWidget(self.canvas, 1)
        cl.addWidget(self._vbar)
        self._message = QLabel(self)
        self._message.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._message.setWordWrap(True)
        self._browser_stack = QStackedWidget(self)
        self._browser_stack.addWidget(self._message)
        self._browser_stack.addWidget(canvas_row)

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
        lay.addLayout(bar)
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
        self._color_btn.setVisible(kind == "alignment")
        text_first = kind in _text_mode_kinds()
        self._set_mode(text=text_first, remember=False)
        if not forced:
            reference = remembered_reference(path)
        self._message.setText(f"Opening {path.name}…")
        self._browser_stack.setCurrentIndex(0)
        self._banner_dismissed_text = ""
        self._banner.hide()
        self._status.setText("")
        project = None if forced else _load_project()
        loader = _Loader(self._token, path, project, reference, forced)
        loader.loaded.connect(self._on_loaded)
        _start_thread(loader)

    def _on_loaded(self, token, res: _LoadResult) -> None:
        if token != self._token:
            release_track(res.track)
            return
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
            self._status.setText(w)

    def _go(self, chrom: str, a: int, b: int) -> None:
        length = (self._result.lengths.get(chrom, 1) if self._result else 1)
        self.canvas.set_chrom(chrom, length, a, b)
        self._sync_vbar()

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
            chrom = next((n for n in res.order if n in res.first), res.order[0])
            span = res.first.get(chrom)
            if span is not None:
                s, e = span
                width = max(e - s, 20_000)
                self._go(chrom, s, s + width)
            else:
                self._go(chrom, 0, res.lengths.get(chrom, 1))
        self._update_reference_ui()

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
        self._status.setText("Reading reference…" if fasta else "")

        job = _RefJob(fasta, track)

        def finish(match, err) -> None:
            if not _alive(self) or token != self._token:
                return
            if err:
                self._status.setText(f"Cannot read the reference: {err}")
                return
            self._status.setText("")
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
        span = self._result.first.get(chrom)
        length = self._result.lengths.get(chrom, 1)
        if span is not None and _is_alignment(self.path):
            self._go(chrom, span[0], span[0] + 20_000)
        else:
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
            self._status.setText(f"Not found: {text}")

    def _show_color_menu(self) -> None:
        menu = QMenu(self)
        menu.setToolTipsVisible(True)
        current = self.canvas.color_by
        for entry in _COLOR_MENU_ENTRIES:
            if entry is None:
                menu.addSeparator()
                continue
            mode, label = entry
            act = menu.addAction(label)
            act.setCheckable(True)
            act.setChecked(mode == current)
            act.triggered.connect(lambda _c=False, m=mode: self._set_color_mode(m))
        menu.exec(self._color_btn.mapToGlobal(self._color_btn.rect().bottomLeft()))

    def _set_color_mode(self, mode: str) -> None:
        if mode not in _COLOR_MODE_LABELS:
            mode = style.COLOR_STRAND
        self.canvas.color_by = mode
        _settings().setValue(_SETTINGS_COLOR, mode)
        self.canvas.update()

    def _sync_locus(self) -> None:
        if not self.canvas.chrom:
            return
        a, b = self.canvas.view_range()
        self._crumbs.set_where(self.canvas.chrom, a, b, "")
        self._status.setText(f"{b - a:,} bp shown")
        self._sync_zoom_slider()
        self._history_timer.start()

    def _sync_vbar(self) -> None:
        total = self.canvas.content_height()
        page = self.canvas._rows_area_h()
        self._vbar.blockSignals(True)
        self._vbar.setRange(0, max(0, total - page))
        self._vbar.setPageStep(max(1, page))
        self._vbar.setSingleStep(READ_H + READ_GAP)
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

    # ------------------------------------------------------------ text mode
    def _set_mode(self, *, text: bool, remember: bool) -> None:
        if self.path is not None and remember:
            _set_text_mode("alignment" if _is_alignment(self.path) else "annotation", text)
        self._browser_btn.setChecked(not text)
        self._text_btn.setChecked(text)
        for w in self._browser_tools:
            w.setEnabled(not text)
        if text:
            self._fill_text()
            self._stack.setCurrentIndex(1)
        else:
            self._stack.setCurrentIndex(0)
            self.canvas.setFocus()

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
