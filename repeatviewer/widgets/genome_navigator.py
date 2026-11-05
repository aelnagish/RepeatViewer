"""Genome Browser navigation: the right-hand rail (every chromosome, bookmarks with
notes, recent places), the breadcrumb bar (back / forward, chromosome › range › gene)
and the history behind back / forward.

The rail and the bar only show and ask; ContigViewPanel moves the view."""

from __future__ import annotations

import time
from dataclasses import asdict, dataclass, field

from PySide6.QtCore import QPointF, QRectF, QSize, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPen
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QPlainTextEdit,
    QScrollArea,
    QSizePolicy,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

BOOKMARK_COLORS = ("#d9534f", "#6f5bd3", "#2a9d8f", "#e08a2c", "#3d7be0")
MAX_RECENT = 12
MAX_HISTORY = 100
RAIL_W = 300


@dataclass
class Bookmark:
    name: str
    chrom: str
    start: int                  # 0-based
    end: int                    # end exclusive; start + 1 for a single position
    note: str = ""
    color: str = BOOKMARK_COLORS[0]

    def locus(self) -> str:
        if self.end - self.start <= 1:
            return f"{self.chrom}:{self.start + 1:,}"
        return f"{self.chrom}:{self.start + 1:,}-{self.end:,}"

    def to_json(self) -> dict:
        return asdict(self)

    @classmethod
    def from_json(cls, d: dict) -> Bookmark | None:
        try:
            return cls(name=str(d.get("name", "")), chrom=str(d["chrom"]),
                       start=int(d["start"]), end=int(d.get("end", int(d["start"]) + 1)),
                       note=str(d.get("note", "")),
                       color=str(d.get("color", BOOKMARK_COLORS[0])))
        except (KeyError, TypeError, ValueError):
            return None


# ---------------------------------------------------------------- BED import / export

def _is_num(v: str) -> bool:
    try:
        float(v)
    except ValueError:
        return False
    return True


def _is_int_list(v: str) -> bool:
    parts = [x for x in v.split(",") if x != ""]
    return bool(parts) and all(x.lstrip("-").isdigit() for x in parts)


def _rgb_hex(v: str) -> str | None:
    parts = v.split(",")
    if len(parts) != 3 or not all(x.strip().isdigit() for x in parts):
        return None
    r, g, b = (max(0, min(255, int(x))) for x in parts)
    if r == g == b == 0:
        return None             # 0,0,0 is BED's "no colour"
    return f"#{r:02x}{g:02x}{b:02x}"


# what a standard BED column 5..12 may hold; the first one that does not fit starts
# the free text (a description) that becomes the bookmark note
_BED_COLUMN_OK = (
    lambda v: v == "." or _is_num(v),                       # score
    lambda v: v in ("+", "-", "."),                          # strand
    lambda v: v.isdigit(),                                   # thickStart
    lambda v: v.isdigit(),                                   # thickEnd
    lambda v: v in ("0", ".") or _rgb_hex(v) is not None or v == "0,0,0",   # itemRgb
    lambda v: v.isdigit(),                                   # blockCount
    _is_int_list,                                            # blockSizes
    _is_int_list,                                            # blockStarts
)


def bookmarks_from_bed(text: str, chroms: dict[str, int] | None = None,
                       first_color: int = 0) -> tuple[list[Bookmark], list[str]]:
    """Bookmarks from BED text: chrom, start (0-based), end (exclusive), then optional
    name, score, strand, thickStart, thickEnd, itemRgb (the bookmark colour) and
    blocks. Anything after the standard columns, or a column that is not what BED
    expects there (a description in column 5, say), becomes the note. ``track`` /
    ``browser`` / ``#`` lines are skipped. With ``chroms`` (name -> length), lines on
    other sequences are skipped and ends are clipped to the sequence. Returns the
    bookmarks and a short reason for each skipped line."""
    out: list[Bookmark] = []
    skipped: list[str] = []
    for ln_no, raw in enumerate(text.splitlines(), 1):
        line = raw.rstrip("\r\n")
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        cols = line.split("\t") if "\t" in line else line.split()
        if len(cols) < 3 or not cols[1].strip().isdigit() or not cols[2].strip().isdigit():
            skipped.append(f"line {ln_no}: needs chrom, start and end")
            continue
        chrom = cols[0].strip()
        start, end = int(cols[1]), int(cols[2])
        if chroms is not None:
            if chrom not in chroms:
                skipped.append(f"line {ln_no}: {chrom} is not in this genome")
                continue
            end = min(end, chroms[chrom])
            if start >= chroms[chrom]:
                skipped.append(f"line {ln_no}: {start + 1:,} is past the end of {chrom}")
                continue
        end = max(end, start + 1)
        name = cols[3].strip() if len(cols) > 3 and cols[3].strip() not in ("", ".") else ""
        color = None
        note_from = len(cols)
        for k, ok in enumerate(_BED_COLUMN_OK):
            idx = 4 + k
            if idx >= len(cols):
                break
            v = cols[idx].strip()
            if not ok(v):
                note_from = idx
                break
            if k == 4:
                color = _rgb_hex(v)
        else:
            note_from = min(note_from, 4 + len(_BED_COLUMN_OK))
        note = " ".join(c.strip() for c in cols[note_from:] if c.strip())
        bm = Bookmark(name=name, chrom=chrom, start=start, end=end, note=note,
                      color=color or BOOKMARK_COLORS[(first_color + len(out))
                                                    % len(BOOKMARK_COLORS)])
        if not bm.name:
            bm.name = bm.locus()
        out.append(bm)
    return out, skipped


def bookmarks_to_bed(bookmarks: list[Bookmark]) -> str:
    """BED9 (+ the note as a 10th column when any bookmark has one): it reads back
    with ``bookmarks_from_bed`` and loads in IGV / UCSC as a coloured track."""
    def clean(v: str) -> str:
        return " ".join(str(v).split())

    any_note = any(b.note.strip() for b in bookmarks)
    rows = ['track name="Bookmarks" itemRgb="On"']
    for b in bookmarks:
        c = QColor(b.color)
        rgb = f"{c.red()},{c.green()},{c.blue()}" if c.isValid() else "0"
        cols = [b.chrom, str(b.start), str(b.end), clean(b.name) or b.locus(), "0", ".",
                str(b.start), str(b.end), rgb]
        if any_note:
            cols.append(clean(b.note) or ".")
        rows.append("\t".join(cols))
    return "\n".join(rows) + "\n"


@dataclass
class Place:
    """A recent place: where the view went on purpose (search, gene, bookmark…)."""
    label: str
    chrom: str
    start: int
    end: int
    when: float = field(default_factory=time.time)

    def to_json(self) -> dict:
        return asdict(self)

    @classmethod
    def from_json(cls, d: dict) -> Place | None:
        try:
            return cls(label=str(d["label"]), chrom=str(d["chrom"]), start=int(d["start"]),
                       end=int(d["end"]), when=float(d.get("when", time.time())))
        except (KeyError, TypeError, ValueError):
            return None


class History:
    """Back / forward: the places the view settled on, newest last."""

    def __init__(self) -> None:
        self.stops: list[tuple[str, int, int]] = []
        self.index = -1

    def record(self, chrom: str, a: int, b: int) -> bool:
        """Remember where the view is now; False when it is (about) the same place."""
        if 0 <= self.index < len(self.stops):
            c0, a0, b0 = self.stops[self.index]
            span = max(1, b0 - a0)
            if c0 == chrom and abs(a - a0) < span * 0.1 and abs((b - a) - span) < span * 0.1:
                return False
        del self.stops[self.index + 1:]
        self.stops.append((chrom, a, b))
        if len(self.stops) > MAX_HISTORY:
            del self.stops[0]
        self.index = len(self.stops) - 1
        return True

    def can_back(self) -> bool:
        return self.index > 0

    def can_forward(self) -> bool:
        return self.index < len(self.stops) - 1

    def back(self) -> tuple[str, int, int] | None:
        if not self.can_back():
            return None
        self.index -= 1
        return self.stops[self.index]

    def forward(self) -> tuple[str, int, int] | None:
        if not self.can_forward():
            return None
        self.index += 1
        return self.stops[self.index]


def span_text(a: int, b: int) -> str:
    """'2.10–2.18 Mb', '12.4–15.0 kb' or '1,204–1,890 bp' for a region."""
    n = b - a
    if n >= 2_000_000 or b >= 10_000_000 and n >= 200_000:
        return f"{a / 1e6:.2f}–{b / 1e6:.2f} Mb"
    if n >= 2_000:
        return f"{a / 1e3:,.1f}–{b / 1e3:,.1f} kb"
    return f"{a + 1:,}–{b:,} bp"


def ago(when: float, now: float | None = None) -> str:
    s = max(0, int((now or time.time()) - when))
    if s < 60:
        return "now"
    if s < 3600:
        return f"{s // 60} min ago"
    if s < 86400:
        return f"{s // 3600} h ago"
    return f"{s // 86400} d ago"


def _palette_colors(w: QWidget) -> tuple[QColor, QColor, QColor, QColor]:
    pal = w.palette()
    base = pal.color(pal.ColorRole.Base)
    dark = base.lightness() < 128
    text = pal.color(pal.ColorRole.Text)
    muted = QColor("#9aa3ad" if dark else "#6b7480")
    line = QColor("#3a3f47" if dark else "#e3e6ea")
    return base, text, muted, line


# ---------------------------------------------------------------- the rail


class ChromosomeMap(QWidget):
    """Every chromosome as a bar (length to scale, feature density as ticks); the one
    in view in bold with a red box where the view is, bookmarks as dots. Click a bar
    to go there, a name for the whole chromosome."""

    goRequested = Signal(str, int, int)        # chrom, start, end (end 0: whole)
    ROW_H = 24
    NAME_W = 70

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.chroms: list[tuple[str, int]] = []
        self.density: dict[str, list[float]] = {}
        self.current = ""
        self.view: tuple[int, int] | None = None
        self.view_w = 0                            # current width in bp (for clicks)
        self.marks: list[tuple[str, int, str]] = []    # chrom, position, colour
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def sizeHint(self) -> QSize:
        return QSize(RAIL_W - 30, max(1, len(self.chroms)) * self.ROW_H + 4)

    def set_chroms(self, chroms: list[tuple[str, int]]) -> None:
        self.chroms = list(chroms)
        self.setFixedHeight(max(1, len(self.chroms)) * self.ROW_H + 4)
        self.update()

    def row_of(self, chrom: str) -> int:
        for i, (c, _n) in enumerate(self.chroms):
            if c == chrom:
                return i
        return -1

    def _bar(self, i: int) -> QRectF:
        longest = max((n for _c, n in self.chroms), default=1) or 1
        x0 = self.NAME_W
        full = max(10.0, self.width() - x0 - 6)
        n = self.chroms[i][1]
        return QRectF(x0, i * self.ROW_H + 7, max(6.0, full * n / longest), 10)

    def paintEvent(self, event) -> None:
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        _base, text, muted, _line = _palette_colors(self)
        accent = QColor("#3d7be0")
        bar_bg = QColor(accent)
        bar_bg.setAlpha(40)
        tick = QColor(accent)
        tick.setAlpha(150)
        f = QFont(self.font())
        f.setPointSize(9)
        bold = QFont(f)
        bold.setBold(True)
        fm = QFontMetrics(f)
        top = event.rect().top() // self.ROW_H
        last = event.rect().bottom() // self.ROW_H + 1
        for i in range(max(0, top), min(len(self.chroms), last + 1)):
            chrom, n = self.chroms[i]
            y = i * self.ROW_H
            here = chrom == self.current
            p.setFont(bold if here else f)
            p.setPen(text if here else muted)
            p.drawText(QRectF(0, y, self.NAME_W - 6, self.ROW_H),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                       fm.elidedText(chrom, Qt.TextElideMode.ElideMiddle, self.NAME_W - 8))
            bar = self._bar(i)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(bar_bg)
            p.drawRoundedRect(bar, 2, 2)
            dens = self.density.get(chrom)
            if dens:
                step = bar.width() / len(dens)
                for k, v in enumerate(dens):
                    if v > 0:
                        c = QColor(tick)
                        c.setAlpha(int(60 + 170 * min(1.0, v)))
                        p.fillRect(QRectF(bar.left() + k * step, bar.top() + 1,
                                          max(1.0, step - 0.3), bar.height() - 2), c)
            for c_, pos, color in self.marks:
                if c_ == chrom and n:
                    x = bar.left() + bar.width() * pos / n
                    p.setPen(QPen(QColor("#ffffff"), 1))
                    p.setBrush(QColor(color))
                    p.drawEllipse(QPointF(x, bar.top() - 1), 3.5, 3.5)
            if here and self.view is not None and n:
                a, b = self.view
                x0 = bar.left() + bar.width() * a / n
                x1 = bar.left() + bar.width() * b / n
                if x1 - x0 < 4:
                    x0, x1 = (x0 + x1) / 2 - 2, (x0 + x1) / 2 + 2
                p.setPen(QPen(QColor("#d9534f"), 1.5))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRect(QRectF(x0, bar.top() - 3, x1 - x0, bar.height() + 6))
        p.end()

    def _hit(self, pos: QPointF) -> tuple[int, int | None]:
        i = int(pos.y() // self.ROW_H)
        if not 0 <= i < len(self.chroms):
            return -1, None
        bar = self._bar(i)
        if pos.x() < bar.left() - 2 or pos.x() > bar.right() + 2:
            return i, None
        n = self.chroms[i][1]
        frac = min(1.0, max(0.0, (pos.x() - bar.left()) / max(1.0, bar.width())))
        return i, int(frac * n)

    def mouseMoveEvent(self, event) -> None:
        i, bp = self._hit(event.position())
        if i < 0:
            self.setToolTip("")
            return
        chrom, n = self.chroms[i]
        tip = f"{chrom} · {n:,} bp"
        if bp is not None:
            tip += f"\nClick to go to {bp + 1:,}"
        else:
            tip += "\nClick to show the whole sequence"
        self.setToolTip(tip)

    def mousePressEvent(self, event) -> None:
        if event.button() != Qt.MouseButton.LeftButton:
            return
        i, bp = self._hit(event.position())
        if i < 0:
            return
        chrom, n = self.chroms[i]
        if bp is None:
            self.goRequested.emit(chrom, 0, n)
            return
        half = max(50, (self.view_w or n // 20) // 2)
        a = max(0, min(n - 2 * half, bp - half))
        self.goRequested.emit(chrom, a, min(n, a + 2 * half))


class _ListRow(QFrame):
    """One row of the Bookmarks / Recent lists: a dot, a title, a place and a note."""

    clicked = Signal()
    menuRequested = Signal(object)

    def __init__(self, title: str, where: str, note: str = "", color: str | None = None,
                 mono_title: bool = False, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("navRow")
        self.setCursor(Qt.CursorShape.ArrowCursor)
        self.setAttribute(Qt.WidgetAttribute.WA_Hover, True)
        _base, text, muted, _line = _palette_colors(self)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 5, 8, 5)
        lay.setSpacing(1)
        head = QHBoxLayout()
        head.setSpacing(8)
        if color:
            dot = QLabel(self)
            dot.setFixedSize(10, 10)
            dot.setStyleSheet(f"background: {color}; border-radius: 5px;")
            head.addWidget(dot, 0, Qt.AlignmentFlag.AlignVCenter)
        self.title = QLabel(title, self)
        tf = QFont(self.font())
        if mono_title:
            tf.setFamily("Menlo")
            tf.setStyleHint(QFont.StyleHint.Monospace)
        else:
            tf.setWeight(QFont.Weight.DemiBold)
        self.title.setFont(tf)
        head.addWidget(self.title, 1)
        self.where = QLabel(where, self)
        wf = QFont(self.font())
        wf.setPointSize(max(8, wf.pointSize() - 2))
        if not mono_title:
            wf.setFamily("Menlo")
            wf.setStyleHint(QFont.StyleHint.Monospace)
        self.where.setFont(wf)
        self.where.setStyleSheet(f"color: {muted.name()};")
        head.addWidget(self.where, 0, Qt.AlignmentFlag.AlignRight)
        lay.addLayout(head)
        self.note = None
        if note:
            self.note = QLabel(note, self)
            self.note.setWordWrap(True)
            self.note.setStyleSheet(f"color: {muted.name()};")
            if color:
                self.note.setContentsMargins(18, 0, 0, 0)
            lay.addWidget(self.note)
        self.setStyleSheet("QFrame#navRow { border-radius: 6px; }"
                           "QFrame#navRow:hover { background: palette(alternate-base); }")
        del text

    def mousePressEvent(self, event) -> None:
        if event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit()
        elif event.button() == Qt.MouseButton.RightButton:
            self.menuRequested.emit(event.globalPosition().toPoint())

    def contextMenuEvent(self, event) -> None:
        event.accept()


# the "Density of" picker: a flat rounded field with a small chevron, like the
# breadcrumb chips
_COMBO_QSS = """
QComboBox#densityCombo {
    border: 1px solid palette(midlight); border-radius: 6px;
    padding: 0px 22px 0px 8px; background: palette(base);
}
QComboBox#densityCombo:hover { border-color: #3d7be0; }
QComboBox#densityCombo:on { border-color: #3d7be0; }
QComboBox#densityCombo::drop-down {
    subcontrol-origin: padding; subcontrol-position: center right;
    width: 20px; border: none; background: transparent;
}
QComboBox#densityCombo::down-arrow { image: none; width: 0px; height: 0px; }
"""


class _ChevronCombo(QComboBox):
    """A combo box whose arrow is a thin painted chevron (the style sheet hides the
    platform arrow)."""

    def paintEvent(self, event) -> None:
        super().paintEvent(event)
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        pen = QPen(self.palette().color(self.palette().ColorRole.PlaceholderText), 1.5)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        p.setPen(pen)
        cx, cy = self.width() - 13.0, self.height() / 2
        p.drawPolyline([QPointF(cx - 4, cy - 2), QPointF(cx, cy + 2), QPointF(cx + 4, cy - 2)])
        p.end()


class _Section(QWidget):
    """A rail section in the style of the main window's side panels: a small ▼ / ▶
    and the title in regular weight (click either to fold), a note or button on the
    right, and its content."""

    toggled = Signal(bool)

    def __init__(self, title: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(6)
        row = QHBoxLayout()
        row.setSpacing(6)
        self.arrow = QLabel(self)
        self.arrow.setFixedWidth(12)
        af = QFont(self.font())
        af.setPointSizeF(max(6.0, af.pointSizeF() * 0.75))
        self.arrow.setFont(af)
        self.arrow.setCursor(Qt.CursorShape.PointingHandCursor)
        self.arrow.mousePressEvent = lambda e: self.head.toggle()  # type: ignore[method-assign]
        row.addWidget(self.arrow)
        self.head = QToolButton(self)
        self.head.setCheckable(True)
        self.head.setChecked(True)
        self.head.setAutoRaise(True)
        self.head.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        self.head.setToolTip(f"Show or hide {title}")
        self.head.setCursor(Qt.CursorShape.PointingHandCursor)
        self.head.setStyleSheet("QToolButton { border: none; padding: 0px; background: transparent; }")
        self.head.setText(title)
        tf = QFont(self.font())
        tf.setPointSizeF(tf.pointSizeF() + 1)
        tf.setBold(True)
        self.head.setFont(tf)
        self._title = title
        self.head.toggled.connect(self._on_toggle)
        row.addWidget(self.head)
        row.addStretch(1)
        self.right = QHBoxLayout()
        row.addLayout(self.right)
        lay.addLayout(row)
        self.body = QWidget(self)
        self.content = QVBoxLayout(self.body)
        self.content.setContentsMargins(0, 0, 0, 0)
        self.content.setSpacing(6)
        lay.addWidget(self.body)
        self._on_toggle(True)

    def _on_toggle(self, on: bool) -> None:
        self.arrow.setText("▼" if on else "▶")
        self.body.setVisible(on)
        self.toggled.emit(on)

    def set_open(self, on: bool) -> None:
        self.head.setChecked(on)


class NavigatorRail(QWidget):
    """The right-hand rail of the Genome Browser: Genome (every sequence, shaded by
    the feature density of one chosen track) and Bookmarks. Each section folds."""

    goRequested = Signal(str, int, int, str)       # chrom, start, end, label
    addBookmarkRequested = Signal()
    editBookmarkRequested = Signal(int)
    deleteBookmarkRequested = Signal(int)
    importBookmarksRequested = Signal()
    exportBookmarksRequested = Signal()
    densityTrackChanged = Signal(str)              # a track path, "" for none
    sectionToggled = Signal(str, bool)             # "genome" / "bookmarks", open

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)


        self.setObjectName("navigatorRail")
        self.setMinimumWidth(220)
        self.bookmarks: list[Bookmark] = []
        self.recent: list[Place] = []
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        outer.addWidget(scroll)
        body = QWidget()
        body.setObjectName("navigatorBody")
        scroll.setWidget(body)
        lay = QVBoxLayout(body)
        lay.setContentsMargins(14, 12, 14, 12)
        lay.setSpacing(8)
        self._body = body
        muted = _palette_colors(self)[2].name()

        # Genome
        self.genome = _Section("Genome", body)
        self.genome.toggled.connect(lambda on: self.sectionToggled.emit("genome", on))
        self.genome_count = QLabel("", body)
        self.genome_count.setStyleSheet(f"color: {muted};")
        self.genome.right.addWidget(self.genome_count)
        dens = QHBoxLayout()
        dens.setSpacing(6)
        dl = QLabel("Density of", body)
        dl.setStyleSheet(f"color: {muted};")
        dens.addWidget(dl)
        self.density_track = _ChevronCombo(body)
        self.density_track.setObjectName("densityCombo")
        self.density_track.setFixedHeight(18)
        sf = QFont(self.density_track.font())
        sf.setPointSizeF(max(8.0, sf.pointSizeF() - 1.5))
        self.density_track.setFont(sf)
        self.density_track.setToolTip(
            "Which track's features shade the sequence bars: darker where that track has "
            "more features (genes, repeats…) per stretch of sequence")
        self.density_track.currentIndexChanged.connect(
            lambda _i: self.densityTrackChanged.emit(str(self.density_track.currentData() or "")))
        dens.addWidget(self.density_track, 1)
        self.genome.content.addLayout(dens)
        self.map = ChromosomeMap(body)
        self.map.goRequested.connect(
            lambda c, a, b: self.goRequested.emit(c, a, b, ""))
        self.map_scroll = QScrollArea(body)
        self.map_scroll.setWidget(self.map)
        self.map_scroll.setWidgetResizable(True)
        self.map_scroll.setFrameShape(QFrame.Shape.NoFrame)
        self.map_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.genome.content.addWidget(self.map_scroll)
        self.filter = QLineEdit(body)
        self.filter.setPlaceholderText("Filter sequences…")
        self.filter.setFixedHeight(18)
        self.filter.setFont(sf)
        self.filter.setClearButtonEnabled(True)
        self.filter.textChanged.connect(self._apply_filter)
        self.genome.content.addWidget(self.filter)
        self._all_chroms: list[tuple[str, int]] = []
        lay.addWidget(self.genome)
        lay.addWidget(self._rule(body))

        # Bookmarks
        self.bm_section = _Section("Bookmarks", body)
        self.bm_section.toggled.connect(lambda on: self.sectionToggled.emit("bookmarks", on))
        self.add_btn = QToolButton(body)
        self.add_btn.setText("+ Add (B)")
        self.add_btn.setAutoRaise(True)
        self.add_btn.setStyleSheet("QToolButton { color: #3d7be0; border: none; }")
        self.add_btn.clicked.connect(self.addBookmarkRequested)
        self.add_btn.setToolTip("Bookmark the selected region or the middle of the view "
                                "(B)")
        self.bm_section.right.addWidget(self.add_btn)
        self.bm_more = QToolButton(body)
        self.bm_more.setText("⋯")
        self.bm_more.setAutoRaise(True)
        self.bm_more.setToolTip("Import bookmarks from a BED file, or export them as BED")
        self.bm_more.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.bm_more.setStyleSheet(
            "QToolButton { border: none; padding: 0px 4px; }"
            "QToolButton::menu-indicator { image: none; width: 0px; }")
        more = QMenu(self.bm_more)
        self.import_action = more.addAction("Import from BED file…")
        self.import_action.triggered.connect(self.importBookmarksRequested)
        self.export_action = more.addAction("Export as BED file…")
        self.export_action.triggered.connect(self.exportBookmarksRequested)
        self.bm_more.setMenu(more)
        self.bm_section.right.addWidget(self.bm_more)
        self.bm_box = QVBoxLayout()
        self.bm_box.setSpacing(2)
        self.bm_section.content.addLayout(self.bm_box)
        self.bm_empty = QLabel("Press B to bookmark the place in view, with a note, or import a BED file from ⋯.", body)
        self.bm_empty.setWordWrap(True)
        self.bm_empty.setStyleSheet(f"color: {muted};")
        self.bm_section.content.addWidget(self.bm_empty)
        lay.addWidget(self.bm_section)
        lay.addWidget(self._rule(body))
        lay.addStretch(1)
        # one flat panel background, like the main window's side panels: the
        # sequence list sits on it with no box of its own
        self.map_scroll.viewport().setAutoFillBackground(False)
        self.map.setAutoFillBackground(False)
        scroll.viewport().setAutoFillBackground(False)
        self.setStyleSheet(
            "QWidget#navigatorRail, QWidget#navigatorBody { background: palette(window); }"
            "QScrollArea { background: transparent; }"
            "QLineEdit { border: 1px solid palette(midlight); border-radius: 6px;"
            " padding: 0px 8px; background: palette(base); }"
            "QLineEdit:focus { border-color: #3d7be0; }"
            + _COMBO_QSS)

    @staticmethod
    def _rule(parent: QWidget) -> QFrame:
        r = QFrame(parent)
        r.setFrameShape(QFrame.Shape.HLine)
        r.setStyleSheet("color: palette(mid);")
        r.setFixedHeight(1)
        return r

    # ------------------------------------------------ genome

    def set_chroms(self, chroms: list[tuple[str, int]]) -> None:
        self._all_chroms = list(chroms)
        n = len(chroms)
        self.genome_count.setText(f"{n:,} sequences" if n > 1 else "")
        self.filter.setVisible(n > 12)
        self._apply_filter()

    def _apply_filter(self, *_a) -> None:
        q = self.filter.text().strip().lower()
        chroms = [c for c in self._all_chroms if q in c[0].lower()] if q else self._all_chroms
        self.map.set_chroms(chroms)
        self.map_scroll.setFixedHeight(min(self.map.minimumHeight(), 10 * ChromosomeMap.ROW_H + 6))
        self._scroll_to_current()

    def set_view(self, chrom: str, a: int, b: int) -> None:
        changed = chrom != self.map.current
        self.map.current = chrom
        self.map.view = (a, b)
        self.map.view_w = b - a
        self.map.update()
        if changed:
            self._scroll_to_current()

    def _scroll_to_current(self) -> None:
        i = self.map.row_of(self.map.current)
        if i >= 0:
            self.map_scroll.ensureVisible(10, i * ChromosomeMap.ROW_H + ChromosomeMap.ROW_H // 2,
                                          0, ChromosomeMap.ROW_H)

    def set_density(self, density: dict[str, list[float]]) -> None:
        self.map.density = density
        self.map.update()

    def set_density_tracks(self, tracks: list[tuple[str, str]], current: str) -> None:
        """The annotation tracks that can shade the bars: (path, name)."""
        combo = self.density_track
        combo.blockSignals(True)
        combo.clear()
        combo.addItem("None", "")
        for path, name in tracks:
            combo.addItem(name, path)
        k = combo.findData(current)
        combo.setCurrentIndex(max(k, 0))
        combo.blockSignals(False)

    # ------------------------------------------------ bookmarks

    @staticmethod
    def _clear(box: QVBoxLayout) -> None:
        while box.count():
            w = box.takeAt(0).widget()
            if w is not None:
                w.deleteLater()

    def set_bookmarks(self, bookmarks: list[Bookmark]) -> None:
        self.bookmarks = list(bookmarks)
        self.export_action.setEnabled(bool(self.bookmarks))
        self._clear(self.bm_box)
        for i, bm in enumerate(self.bookmarks):
            row = _ListRow(bm.name or bm.locus(), bm.locus(), bm.note, bm.color,
                           parent=self._body)
            row.setToolTip("Click to go there · right-click to edit or delete")
            row.clicked.connect(lambda b=bm: self.goRequested.emit(
                b.chrom, b.start, b.end, b.name or b.locus()))
            row.menuRequested.connect(lambda pos, k=i: self._bookmark_menu(k, pos))
            self.bm_box.addWidget(row)
        self.bm_empty.setVisible(not self.bookmarks)
        self.map.marks = [(b.chrom, (b.start + b.end) // 2, b.color) for b in self.bookmarks]
        self.map.update()

    def _bookmark_menu(self, k: int, pos) -> None:
        menu = QMenu(self)
        bm = self.bookmarks[k]
        menu.addAction("Go there", lambda: self.goRequested.emit(
            bm.chrom, bm.start, bm.end, bm.name or bm.locus()))
        menu.addAction("Edit…", lambda: self.editBookmarkRequested.emit(k))
        menu.addAction("Delete", lambda: self.deleteBookmarkRequested.emit(k))
        menu.exec(pos)

    # recent places are kept (for the session) but not shown
    def set_recent(self, recent: list[Place]) -> None:
        self.recent = list(recent)

    def refresh_recent(self) -> None:
        pass


class ChromPicker(QFrame):
    """The sequence list under the chromosome chip: a filter box and one scrolling
    column, never taller than about 14 rows."""

    picked = Signal(str)

    def __init__(self, chroms: list[tuple[str, int]], current: str,
                 parent: QWidget | None = None) -> None:
        from PySide6.QtWidgets import QListWidget, QListWidgetItem

        super().__init__(parent, Qt.WindowType.Popup | Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        # a rounded card like the menus: only the card is painted, the window
        # around its corners stays transparent
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.setFrameShape(QFrame.Shape.NoFrame)
        pal = self.palette()
        base = pal.color(pal.ColorRole.Base)
        dark = base.lightness() < 128
        self._card = (QColor(base), QColor("#3a3f47" if dark else "#d5d9e0"))
        self.setStyleSheet(
            "QListWidget { border: none; background: transparent; }"
            "QListWidget::item { padding: 2px 6px; border-radius: 5px; }"
            "QListWidget::item:selected { background: #3d7be0; color: white; }"
            "QListWidget::item:hover:!selected { background: rgba(61, 123, 224, 40); }")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(4)
        self.filter = QLineEdit(self)
        self.filter.setPlaceholderText(f"Filter {len(chroms):,} sequences…")
        self.filter.setClearButtonEnabled(True)
        lay.addWidget(self.filter)
        self.list = QListWidget(self)
        self.list.setUniformItemSizes(True)
        mono = QFont(self.font())
        for name, n in chroms:
            it = QListWidgetItem(f"{name}   {n:,} bp")
            it.setData(Qt.ItemDataRole.UserRole, name)
            it.setFont(mono)
            self.list.addItem(it)
            if name == current:
                self.list.setCurrentItem(it)
        lay.addWidget(self.list)
        rows = min(14, max(3, len(chroms)))
        self.list.setFixedHeight(rows * (self.list.sizeHintForRow(0) or 20) + 6)
        self.setFixedWidth(max(240, self.list.sizeHintForColumn(0) + 40))
        self.filter.textChanged.connect(self._filter)
        self.filter.returnPressed.connect(self._pick_first)
        self.list.itemClicked.connect(self._pick)
        self.list.itemActivated.connect(self._pick)
        if self.list.currentItem() is not None:
            self.list.scrollToItem(self.list.currentItem(),
                                   self.list.ScrollHint.PositionAtCenter)
        self.filter.setFocus()

    def paintEvent(self, event) -> None:
        fill, edge = self._card
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.setPen(QPen(edge, 1))
        p.setBrush(fill)
        p.drawRoundedRect(QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5), 10, 10)
        p.end()

    def _filter(self, text: str) -> None:
        q = text.strip().lower()
        for i in range(self.list.count()):
            it = self.list.item(i)
            it.setHidden(bool(q) and q not in str(it.data(Qt.ItemDataRole.UserRole)).lower())

    def _pick_first(self) -> None:
        for i in range(self.list.count()):
            it = self.list.item(i)
            if not it.isHidden():
                self._pick(it)
                return

    def _pick(self, item) -> None:
        self.picked.emit(str(item.data(Qt.ItemDataRole.UserRole)))
        self.close()


# ---------------------------------------------------------------- breadcrumb bar


_CHIP_CSS = """
QToolButton#navChip {{ border: 1px solid {line}; border-radius: 7px; padding: 3px 10px;
                       background: palette(base); }}
QToolButton#navChip:hover {{ border-color: #9aa3ad; }}
QToolButton#navChip[current="true"] {{ border: 1px solid #3d7be0; color: #2a62c9;
                                        background: {sel}; font-weight: 600; }}
QToolButton#navRound {{ border: 1px solid {line}; border-radius: 13px; min-width: 24px;
                        min-height: 24px; max-width: 24px; max-height: 24px;
                        background: palette(base); }}
QToolButton#navRound:disabled {{ color: #b5bcc4; }}
QLabel#navSep {{ color: #9aa3ad; }}
"""


class _ArrowButton(QToolButton):
    """A round Back / Forward button whose chevron is drawn at the exact centre
    (a font's ‹ › glyphs sit off-centre)."""

    def __init__(self, direction: int, parent=None) -> None:
        super().__init__(parent)
        self.direction = direction

    def paintEvent(self, event) -> None:
        super().paintEvent(event)       # the circle, from the style sheet
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        col = self.palette().color(self.foregroundRole())
        if not self.isEnabled():
            col = QColor("#b5bcc4")
        pen = QPen(col, 1.8)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        p.setPen(pen)
        cx, cy = self.width() / 2, self.height() / 2
        w, h = 2.6, 5.0                 # half width and half height of the chevron
        tip = cx + self.direction * w
        back = cx - self.direction * w
        p.drawPolyline([QPointF(back, cy - h), QPointF(tip, cy), QPointF(back, cy + h)])
        p.end()


class BreadcrumbBar(QWidget):
    """‹ › chr4 › 2.10–2.18 Mb › BLGH_03307: back / forward and where the view is.
    The chromosome chip lists the sequences, the range chip edits the position, the
    gene chip zooms to the gene in the middle of the view."""

    backRequested = Signal()
    forwardRequested = Signal()
    chromRequested = Signal(str)
    rangeClicked = Signal()
    geneClicked = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(6)
        self.back = self._round("‹", "Back (⌘[)", self.backRequested)
        self.forward = self._round("›", "Forward (⌘])", self.forwardRequested)
        lay.addWidget(self.back)
        lay.addWidget(self.forward)
        lay.addSpacing(4)
        self.chrom = self._chip("Sequence: pick another one")
        self.chrom.clicked.connect(self._open_picker)
        self.picker: ChromPicker | None = None
        self.range = self._chip("The region in view: click to type another position")
        self.range.clicked.connect(self.rangeClicked)
        self.gene = self._chip("The feature in the middle of the view: click to zoom to it")
        self.gene.clicked.connect(self.geneClicked)
        self.gene.setProperty("current", True)
        self.seps = [self._sep(), self._sep()]
        lay.addWidget(self.chrom)
        lay.addWidget(self.seps[0])
        lay.addWidget(self.range)
        lay.addWidget(self.seps[1])
        lay.addWidget(self.gene)
        self._chroms: list[tuple[str, int]] = []
        line = "#3a3f47" if self.palette().color(
            self.palette().ColorRole.Base).lightness() < 128 else "#d5d9e0"
        self.setStyleSheet(_CHIP_CSS.format(line=line, sel="rgba(61,123,224,0.10)"))
        self.set_where("", 0, 0, "")
        self.set_history(False, False)

    def _round(self, text: str, tip: str, signal) -> QToolButton:
        b = _ArrowButton(-1 if text == "‹" else 1, self)
        b.setObjectName("navRound")
        b.setToolTip(tip)
        b.clicked.connect(signal)
        return b

    def _chip(self, tip: str) -> QToolButton:
        b = QToolButton(self)
        b.setObjectName("navChip")
        b.setToolTip(tip)
        b.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        return b

    def _sep(self) -> QLabel:
        s = QLabel("›", self)
        s.setObjectName("navSep")
        return s

    def set_chroms(self, chroms: list[tuple[str, int]]) -> None:
        self._chroms = list(chroms)

    def _open_picker(self) -> None:
        self.picker = ChromPicker(self._chroms, self.chrom.text(), self)
        self.picker.picked.connect(self.chromRequested)
        self.picker.move(self.chrom.mapToGlobal(self.chrom.rect().bottomLeft()))
        self.picker.show()

    def set_where(self, chrom: str, a: int, b: int, gene: str) -> None:
        self.chrom.setText(chrom or "—")
        self.chrom.setVisible(bool(chrom))
        has_range = bool(chrom) and b > a
        self.range.setText(span_text(a, b) if has_range else "")
        self.range.setVisible(has_range)
        self.seps[0].setVisible(has_range)
        self.gene.setText(gene)
        self.gene.setVisible(bool(gene))
        self.seps[1].setVisible(bool(gene))

    def set_history(self, can_back: bool, can_forward: bool) -> None:
        self.back.setEnabled(can_back)
        self.forward.setEnabled(can_forward)


# ---------------------------------------------------------------- bookmark dialog


class BookmarkDialog(QDialog):
    def __init__(self, bm: Bookmark, parent: QWidget | None = None, title: str = "Bookmark"
                 ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.color = bm.color
        form = QFormLayout(self)
        self.name = QLineEdit(bm.name, self)
        self.name.selectAll()
        form.addRow("Name", self.name)
        where = QLabel(bm.locus(), self)
        where.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        form.addRow("Place", where)
        self.note = QPlainTextEdit(bm.note, self)
        self.note.setPlaceholderText("A note, e.g. exon 3 junction looks short")
        self.note.setFixedHeight(64)
        form.addRow("Note", self.note)
        swatches = QHBoxLayout()
        self._swatches = []
        for c in BOOKMARK_COLORS:
            b = QToolButton(self)
            b.setCheckable(True)
            b.setChecked(c == bm.color)
            b.setFixedSize(22, 22)
            b.setStyleSheet(
                f"QToolButton {{ background: {c}; border-radius: 11px; border: 2px solid "
                f"transparent; }} QToolButton:checked {{ border: 2px solid palette(text); }}")
            b.clicked.connect(lambda _c=False, c_=c: self._pick(c_))
            swatches.addWidget(b)
            self._swatches.append((c, b))
        swatches.addStretch(1)
        form.addRow("Colour", swatches)
        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok
                                   | QDialogButtonBox.StandardButton.Cancel, self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        form.addRow(buttons)
        self._bm = bm

    def _pick(self, color: str) -> None:
        self.color = color
        for c, b in self._swatches:
            b.setChecked(c == color)

    def result_bookmark(self) -> Bookmark:
        bm = self._bm
        return Bookmark(name=self.name.text().strip() or bm.locus(), chrom=bm.chrom,
                        start=bm.start, end=bm.end,
                        note=self.note.toPlainText().strip(), color=self.color)
