"""Binding Sites on Chromosomes: where primers bind across a genome.

Each sequence is a strip on a shared base-pair scale; forward sites are ticks
above it, reverse sites below, coloured by mismatches (darker = fewer).
Products of the chosen pair are an orange band. The view zooms (scroll wheel, pinch, − +) and
moves (drag, trackpad swipe); double-click shows everything again.
The side panel counts sites per sequence for the primer(s) shown and the
products of the chosen pair (the expected value for a qPCR copy-number assay).
"""

from __future__ import annotations

import html
import math
import threading
from pathlib import Path

from PySide6.QtCore import QEvent, QObject, QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QPainter, QPen
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QTableWidget,
    QTableWidgetItem,
    QToolTip,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import binding as bd
from repeatviewer.backends.primers import pcr
from repeatviewer.backends.primers.oligos import Primer
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.path_browse import PathBrowseRow

MM_COLORS = ("#0f3a7a", "#4f86d6", "#a9c6ef", "#d6e3f5", "#e9eef7")
PRODUCT_COLOR = "#f2a65a"
WARN_COLOR = "#b26b00"
MAX_STRIPS = 60
NAME_W = 120
ROW_H = 46
RULER_H = 30
ALL_PRIMERS = "__all__"

GENOME_HELP = (
    "<b>Genome</b>: the FASTA to search (the project genome is filled in; click a file in "
    "the Files panel to change it).<br>"
    "Each sequence is drawn on the same base-pair scale. <b>Forward</b> sites are ticks above "
    "the bar, <b>reverse</b> sites below; darker blue means fewer mismatches. The <b>orange</b> "
    "band under a strip shows where the qPCR pair makes products.<br>"
    "<b>Zoom</b>: scroll wheel or pinch over the strips, or the − and + buttons. "
    "<b>Move</b>: drag, or swipe sideways on a trackpad. <b>Double-click</b> or <b>All</b> "
    "shows the whole sequences again; <b>Sites</b> zooms to where the primers bind. "
    "Click a row in the table to zoom to that sequence.<br>"
    "<b>Show</b> limits ticks and counts to one primer. <b>qPCR pair</b>: products of these "
    "two primers per haploid genome, the value to expect in a copy-number assay.<br>"
    "Only the 60 longest sequences are drawn; the table lists every sequence with a site."
)


class _Relay(QObject):
    done = Signal(object)
    failed = Signal(str)
    progress = Signal(str)


def nice_step(span: float, target_ticks: int = 8) -> float:
    raw = max(1.0, span / max(1, target_ticks))
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 5, 10):
        if m * mag >= raw:
            return m * mag
    return 10 * mag


def fmt_bp(v: float, step: float) -> str:
    if step >= 1_000_000:
        return f"{v / 1e6:g} Mb"
    if step >= 1_000:
        return f"{v / 1e3:g} kb"
    return f"{int(v):,}"


class GenomeStrips(QWidget):
    """Strips on a shared, zoomable base-pair window [lo, hi]."""

    viewChanged = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._records: list[pcr.RecordSites] = []
        self._shown: set[str] | None = None
        self._lo = 0.0
        self._hi = 1.0
        self._drag_x: float | None = None
        self._drag_lo = 0.0
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setMinimumWidth(500)
        self.grabGesture(Qt.GestureType.PinchGesture)

    # ---------------------------------------------------------------- data

    def set_records(self, records: list[pcr.RecordSites]) -> None:
        recs = sorted(records, key=lambda r: -r.length)[:MAX_STRIPS]
        order = {r.record: i for i, r in enumerate(records)}
        self._records = sorted(recs, key=lambda r: order[r.record])
        self.setFixedHeight(RULER_H + 10 + ROW_H * max(1, len(self._records)))
        self.fit_all()

    def set_shown(self, primers: set[str] | None) -> None:
        self._shown = primers
        self.update()

    def set_pair(self, left: str, right: str) -> None:
        self._pair = (left, right)
        self.update()

    def _pair_spans(self, rs: pcr.RecordSites) -> list[tuple[int, int]]:
        """Merged intervals covered by products of the chosen pair."""
        left, right = getattr(self, "_pair", ("", ""))
        spans = sorted((a, b) for a, b, _s, f, r in rs.products
                       if {f, r} <= {left, right} and (f != r or left == right))
        out: list[list[int]] = []
        for a, b in spans:
            if out and a <= out[-1][1]:
                out[-1][1] = max(out[-1][1], b)
            else:
                out.append([a, b])
        return [(a, b) for a, b in out]

    def records(self) -> list[pcr.RecordSites]:
        return list(self._records)

    def row_of(self, record: str) -> int:
        for i, r in enumerate(self._records):
            if r.record == record:
                return i
        return -1

    # ---------------------------------------------------------------- view

    def longest(self) -> int:
        return max((r.length for r in self._records), default=1) or 1

    def view(self) -> tuple[float, float]:
        return self._lo, self._hi

    def set_view(self, lo: float, hi: float) -> None:
        longest = float(self.longest())
        span = max(50.0, min(hi - lo, longest * 1.02))
        lo = max(-longest * 0.01, min(lo, longest * 1.01 - span))
        self._lo, self._hi = lo, lo + span
        self.update()
        self.viewChanged.emit()

    def fit_all(self) -> None:
        self.set_view(0.0, self.longest() * 1.01)

    def fit_sites(self, record: str | None = None) -> None:
        pos = []
        for r in self._records:
            if record is not None and r.record != record:
                continue
            pos.extend(p for name, _s, p, _m in r.sites if self._visible(name))
        if not pos:
            if record is not None:
                rec = next((r for r in self._records if r.record == record), None)
                if rec is not None:
                    self.set_view(0, rec.length)
            return
        a, b = min(pos), max(pos)
        pad = max(200.0, (b - a) * 0.08)
        self.set_view(a - pad, b + pad)

    def zoom(self, factor: float, at_x: float | None = None) -> None:
        left, width = self._plot()
        if at_x is None:
            at_x = left + width / 2
        frac = min(1.0, max(0.0, (at_x - left) / max(1.0, width)))
        anchor = self._lo + frac * (self._hi - self._lo)
        span = (self._hi - self._lo) / factor
        self.set_view(anchor - frac * span, anchor - frac * span + span)

    def pan_pixels(self, dx: float) -> None:
        _left, width = self._plot()
        bp = dx / max(1.0, width) * (self._hi - self._lo)
        self.set_view(self._lo - bp, self._hi - bp)

    def _plot(self) -> tuple[float, float]:
        left = NAME_W + 12.0
        return left, max(100.0, self.width() - left - 16.0)

    def _x(self, bp: float) -> float:
        left, width = self._plot()
        return left + (bp - self._lo) / max(1e-9, self._hi - self._lo) * width

    def _visible(self, primer: str) -> bool:
        return self._shown is None or primer in self._shown

    # ---------------------------------------------------------------- paint

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        text = self.palette().color(self.foregroundRole())
        muted = QColor(text)
        muted.setAlpha(140)
        if not self._records:
            p.setPen(text)
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter,
                       "Choose a genome and press Find Sites")
            p.end()
            return
        left, width = self._plot()
        fm = p.fontMetrics()
        # ruler
        span = self._hi - self._lo
        step = nice_step(span)
        first = math.ceil(max(0.0, self._lo) / step) * step
        p.setPen(muted)
        p.drawLine(QPointF(left, RULER_H - 6), QPointF(left + width, RULER_H - 6))
        v = first
        while v <= self._hi:
            x = self._x(v)
            if left <= x <= left + width:
                p.drawLine(QPointF(x, RULER_H - 10), QPointF(x, RULER_H - 6))
                label = fmt_bp(v, step)
                p.drawText(QRectF(x - 40, 2, 80, 16), Qt.AlignmentFlag.AlignHCenter, label)
            v += step
        p.setClipRect(QRectF(0, RULER_H, self.width(), self.height()))
        y = RULER_H + 8.0
        bp_per_px = span / max(1.0, width)
        for rs in self._records:
            p.setClipping(False)
            p.setPen(text)
            p.drawText(QRectF(0, y + 4, NAME_W, 20),
                       Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                       fm.elidedText(rs.record, Qt.TextElideMode.ElideMiddle, NAME_W))
            p.setClipRect(QRectF(left, y - 6, width, ROW_H))
            x0 = self._x(0)
            x1 = self._x(rs.length)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor("#e1e1e6"))
            p.drawRoundedRect(QRectF(x0, y + 9, max(3.0, x1 - x0), 10), 5, 5)
            # products of the qPCR pair: where they are, as one band under the strip
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(PRODUCT_COLOR))
            for a, b in self._pair_spans(rs):
                if b < self._lo or a > self._hi:
                    continue
                xa, xb = self._x(a), self._x(b)
                p.drawRoundedRect(QRectF(xa, y + 33, max(3.0, xb - xa), 5), 2, 2)
            # sites: one tick per pixel column is enough when zoomed out
            p.setPen(Qt.PenStyle.NoPen)
            seen: set[tuple[int, str, int]] = set()
            tick_w = 2.0 if bp_per_px > 2 else max(2.0, 20 / bp_per_px)
            for primer, strand, pos, mm in rs.sites:
                if not self._visible(primer) or pos < self._lo or pos > self._hi:
                    continue
                x = self._x(pos)
                key = (int(x), strand, min(mm, 4))
                if key in seen:
                    continue
                seen.add(key)
                p.setBrush(QColor(MM_COLORS[min(mm, len(MM_COLORS) - 1)]))
                if strand == "+":
                    p.drawRect(QRectF(x - 1, y - 2, min(tick_w, 8), 11))
                else:
                    p.drawRect(QRectF(x - 1, y + 19, min(tick_w, 8), 11))
            y += ROW_H
        p.end()

    # ---------------------------------------------------------------- input

    def wheelEvent(self, event) -> None:  # noqa: N802
        mods = event.modifiers()
        pixel = event.pixelDelta()
        angle = event.angleDelta()
        trackpad = event.phase() != Qt.ScrollPhase.NoScrollPhase or not pixel.isNull()
        x = event.position().x()
        if mods & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier):
            dy = angle.y() or pixel.y()
            self.zoom(1.0015 ** dy, x)
            event.accept()
            return
        if trackpad:
            if abs(pixel.x()) > abs(pixel.y()):
                self.pan_pixels(pixel.x())
                event.accept()
            else:
                event.ignore()      # vertical swipe scrolls the list of strips
            return
        if mods & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.AltModifier):
            self.pan_pixels((angle.y() or angle.x()) / 2)
            event.accept()
            return
        self.zoom(1.25 if angle.y() > 0 else 0.8, x)
        event.accept()

    def event(self, ev) -> bool:  # noqa: D401
        if ev.type() == QEvent.Type.NativeGesture:
            if ev.gestureType() == Qt.NativeGestureType.ZoomNativeGesture:
                self.zoom(1.0 + ev.value(), ev.position().x())
                return True
            if ev.gestureType() == Qt.NativeGestureType.SmartZoomNativeGesture:
                self.fit_all()
                return True
        if ev.type() == QEvent.Type.Gesture:
            pinch = ev.gesture(Qt.GestureType.PinchGesture)
            if pinch is not None:
                self.zoom(pinch.scaleFactor())
                return True
        return super().event(ev)

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self._drag_x = event.position().x()
            self._drag_lo = self._lo
            self.setCursor(Qt.CursorShape.ClosedHandCursor)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if self._drag_x is not None:
            _left, width = self._plot()
            dx = event.position().x() - self._drag_x
            span = self._hi - self._lo
            lo = self._drag_lo - dx / max(1.0, width) * span
            self.set_view(lo, lo + span)
            return
        hit = self._site_at(event.position())
        if hit:
            QToolTip.showText(event.globalPosition().toPoint(), hit, self)
        else:
            QToolTip.hideText()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        self._drag_x = None
        self.unsetCursor()

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        self.fit_all()

    def _site_at(self, pos: QPointF) -> str:
        row = int((pos.y() - RULER_H - 8 + 6) // ROW_H)
        if not (0 <= row < len(self._records)):
            return ""
        rs = self._records[row]
        bp_per_px = (self._hi - self._lo) / max(1.0, self._plot()[1])
        target = self._lo + (pos.x() - self._plot()[0]) * bp_per_px
        near = [s for s in rs.sites if self._visible(s[0]) and abs(s[2] - target) <= 4 * bp_per_px]
        if not near:
            return ""
        lines = [f"<b>{html.escape(rs.record)}</b>"]
        for primer, strand, p, mm in near[:6]:
            lines.append(f"{html.escape(primer)} {'→' if strand == '+' else '←'} {p:,} · {mm} mm")
        if len(near) > 6:
            lines.append(f"… {len(near) - 6} more here")
        return "<br>".join(lines)


class PrimerGenomeDialog(QDialog):
    fileWritten = Signal(object)

    def __init__(self, primers: list[Primer], *, genome: Path | None = None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Binding Sites on Chromosomes")
        self.setMinimumSize(1200, 760)
        self._primers = list(primers)
        self._records: list[pcr.RecordSites] = []
        self._relay = _Relay(self)
        self._relay.done.connect(self._on_done)
        self._relay.failed.connect(self._on_failed)
        self._relay.progress.connect(lambda t: self.status.setText(t))

        root = QVBoxLayout(self)
        root.setContentsMargins(16, 12, 16, 12)
        root.setSpacing(10)
        top = QHBoxLayout()
        top.setSpacing(8)
        top.addWidget(QLabel("Genome", self))
        self.genome_row = PathBrowseRow(self, initial=genome, placeholder="Genome FASTA",
                                        browse_title="Choose the genome FASTA",
                                        file_filter="FASTA (*.fa *.fasta *.fna *.fa.gz);;All files (*)")
        self.genome_row.setFixedWidth(420)
        top.addWidget(self.genome_row)
        top.addSpacing(8)
        top.addWidget(QLabel("Mismatches", self))
        self.mm_combo = QComboBox(self)
        for n in range(0, 4):
            self.mm_combo.addItem("none" if n == 0 else f"up to {n}", n)
        self.mm_combo.setCurrentIndex(2)
        self.mm_combo.setFixedWidth(110)
        top.addWidget(self.mm_combo)
        self.find_button = QPushButton("Find Sites", self)
        self.find_button.setDefault(True)
        top.addWidget(self.find_button)
        top.addStretch(1)
        top.addWidget(make_help_button(self, "Binding Sites on Chromosomes", GENOME_HELP))
        root.addLayout(top)

        second = QHBoxLayout()
        second.setSpacing(8)
        second.addWidget(QLabel("Show", self))
        self.show_combo = QComboBox(self)
        self.show_combo.addItem(f"All primers ({len(self._primers)})", ALL_PRIMERS)
        for p in self._primers:
            self.show_combo.addItem(p.name, p.name)
        self.show_combo.setFixedWidth(240)
        second.addWidget(self.show_combo)
        second.addSpacing(8)
        second.addWidget(QLabel("qPCR pair", self))
        self.left_combo = QComboBox(self)
        self.right_combo = QComboBox(self)
        for c in (self.left_combo, self.right_combo):
            for p in self._primers:
                c.addItem(p.name, p.name)
            c.setFixedWidth(240)
        if len(self._primers) > 1:
            self.right_combo.setCurrentIndex(1)
        second.addWidget(self.left_combo)
        second.addWidget(QLabel("+", self))
        second.addWidget(self.right_combo)
        second.addStretch(1)
        legend = QLabel(
            " ".join(f'<span style="color:{c}">■</span> {i}' for i, c in enumerate(MM_COLORS[:3]))
            + " mismatches &nbsp; "
            f'<span style="color:{PRODUCT_COLOR}">▬</span> products of the pair', self)
        second.addWidget(legend)
        root.addLayout(second)

        body = QHBoxLayout()
        body.setSpacing(12)
        left = QVBoxLayout()
        left.setSpacing(6)
        nav = QHBoxLayout()
        nav.setSpacing(6)
        self.zoom_out = QPushButton("−", self)
        self.zoom_in = QPushButton("+", self)
        for b in (self.zoom_out, self.zoom_in):
            b.setFixedWidth(34)
        self.fit_all_button = QPushButton("All", self)
        self.fit_sites_button = QPushButton("Sites", self)
        self.fit_all_button.setToolTip("Show the whole sequences (double-click does the same)")
        self.fit_sites_button.setToolTip("Zoom to where the shown primers bind")
        for b in (self.zoom_out, self.zoom_in, self.fit_all_button, self.fit_sites_button):
            nav.addWidget(b)
        nav.addSpacing(10)
        self.view_label = QLabel(self)
        self.view_label.setStyleSheet("color: palette(placeholder-text);")
        nav.addWidget(self.view_label, 1)
        left.addLayout(nav)
        self.strips = GenomeStrips(self)
        self.scroll = QScrollArea(self)
        self.scroll.setWidgetResizable(True)
        self.scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.scroll.setWidget(self.strips)
        left.addWidget(self.scroll, 1)
        body.addLayout(left, 1)

        side_host = QWidget(self)
        side_host.setFixedWidth(360)
        side = QVBoxLayout(side_host)
        side.setContentsMargins(0, 0, 0, 0)
        side.setSpacing(8)
        side.addWidget(QLabel("<b>Sites per sequence</b>", side_host))
        self.table = QTableWidget(0, 5, side_host)
        self.table.setHorizontalHeaderLabels(["Sequence", "0 mm", "1 mm", "2+ mm", "Pair products"])
        self.table.verticalHeader().hide()
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setShowGrid(False)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        h = self.table.horizontalHeader()
        h.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        h.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        h.setStretchLastSection(True)
        for c, w in enumerate((96, 46, 46, 50)):
            self.table.setColumnWidth(c, w)
        side.addWidget(self.table, 1)
        self.summary = QLabel(side_host)
        self.summary.setWordWrap(True)
        self.summary.setTextFormat(Qt.TextFormat.RichText)
        self.summary.setStyleSheet("QLabel { border: 1px solid palette(mid); border-radius: 8px;"
                                   " padding: 8px; }")
        side.addWidget(self.summary)
        self.gff_button = QPushButton("Save Sites as GFF3…", side_host)
        self.png_button = QPushButton("Save Figure (PNG)…", side_host)
        side.addWidget(self.gff_button)
        side.addWidget(self.png_button)
        body.addWidget(side_host)
        root.addLayout(body, 1)
        self.status = QLabel(self)
        self.status.setWordWrap(True)
        self.status.setStyleSheet("color: palette(placeholder-text);")
        root.addWidget(self.status)

        self.find_button.clicked.connect(self._find)
        self.gff_button.clicked.connect(self._save_gff)
        self.png_button.clicked.connect(self._save_png)
        self.zoom_in.clicked.connect(lambda: self.strips.zoom(2.0))
        self.zoom_out.clicked.connect(lambda: self.strips.zoom(0.5))
        self.fit_all_button.clicked.connect(self.strips.fit_all)
        self.fit_sites_button.clicked.connect(lambda: self.strips.fit_sites())
        self.strips.viewChanged.connect(self._sync_view_label)
        self.show_combo.currentIndexChanged.connect(lambda _i: self._refresh())
        self.left_combo.currentIndexChanged.connect(lambda _i: self._refresh())
        self.right_combo.currentIndexChanged.connect(lambda _i: self._refresh())
        self.table.itemSelectionChanged.connect(self._on_row)
        for w in (self.gff_button, self.png_button, self.zoom_in, self.zoom_out,
                  self.fit_all_button, self.fit_sites_button):
            w.setEnabled(False)
        self.summary.hide()

    # ---------------------------------------------------------------- input

    def apply_path_from_panel(self, path: Path) -> bool:
        self.genome_row.set_path(path)
        return True

    def _shown(self) -> set[str] | None:
        key = self.show_combo.currentData()
        return None if key == ALL_PRIMERS else {str(key)}

    def _find(self) -> None:
        genome = self.genome_row.path()
        if genome is None or not Path(genome).is_file():
            QMessageBox.information(self, "Binding Sites", "Choose the genome FASTA first.")
            return
        if not self._primers:
            QMessageBox.information(self, "Binding Sites",
                                    "No primers. Open this window from the Primers tab of a sequence.")
            return
        self.find_button.setEnabled(False)
        self.status.setText(f"Searching {Path(genome).name}…")
        mm = int(self.mm_combo.currentData())
        relay = self._relay
        primers = list(self._primers)

        def work() -> None:
            try:
                recs = pcr.genome_sites(Path(genome), primers, max_mismatches=mm,
                                        progress=relay.progress.emit)
                relay.done.emit(recs)
            except Exception as exc:  # noqa: BLE001
                relay.failed.emit(str(exc))

        threading.Thread(target=work, daemon=True).start()

    def _on_failed(self, error: str) -> None:
        self.find_button.setEnabled(True)
        self.status.setText("")
        QMessageBox.warning(self, "Binding Sites", error)

    # ---------------------------------------------------------------- results

    def _on_done(self, records) -> None:
        self.find_button.setEnabled(True)
        self._records = list(records)
        self.strips.set_records(self._records)
        self.strips.set_shown(self._shown())
        for w in (self.png_button, self.zoom_in, self.zoom_out, self.fit_all_button,
                  self.fit_sites_button):
            w.setEnabled(True)
        # Few sequences: start on the sites, the whole sequence says little.
        if len(self._records) <= 3:
            self.strips.fit_sites()
        self._refresh()

    def _refresh(self) -> None:
        shown = self._shown()
        self.strips.set_shown(shown)
        left = str(self.left_combo.currentData() or "")
        right = str(self.right_combo.currentData() or "")
        self.strips.set_pair(left, right)
        with_sites = [r for r in self._records if r.total_sites(shown)]
        self.table.setRowCount(len(with_sites))
        for i, r in enumerate(with_sites):
            c = r.counts(shown)
            vals = [r.record, c.get(0, 0), c.get(1, 0), sum(n for mm, n in c.items() if mm >= 2),
                    r.pair_products(left, right)]
            for j, v in enumerate(vals):
                item = QTableWidgetItem(f"{v:,}" if isinstance(v, int) else str(v))
                if j == 0:
                    item.setToolTip(f"{r.record} · {r.length:,} bp")
                self.table.setItem(i, j, item)
        total = sum(r.total_sites(shown) for r in self._records)
        pair = sum(r.pair_products(left, right) for r in self._records)
        which = "all primers" if shown is None else html.escape(next(iter(shown)))
        lines = [f"<b>{total:,}</b> sites of {which} on <b>{len(with_sites)}</b> of "
                 f"{len(self._records):,} sequences."]
        if left and right:
            sizes: dict[int, int] = {}
            for r in self._records:
                for _a, _b, size, f, rv in r.products:
                    if {f, rv} <= {left, right} and (f != rv or left == right):
                        sizes[size] = sizes.get(size, 0) + 1
            lines.append(f"qPCR pair <b>{html.escape(left)}</b> + <b>{html.escape(right)}</b>: "
                         f"<b>{pair:,}</b> products up to 5 kb.")
            if sizes:
                short = min(sizes)
                lines.append(f"Shortest product {short:,} bp, <b>{sizes[short]:,}</b> copies: "
                             "the copy number a qPCR with a short extension measures"
                             + (" (tandem repeat: longer products are multiples)."
                                if len(sizes) > 1 else "."))
        warn = []
        if any(r.products_truncated for r in self._records):
            warn.append(f"more than {pcr.MAX_PRODUCTS_PER_RECORD:,} products on one sequence: "
                        "counting stopped there")
        per_primer: dict[str, int] = {}
        for r in self._records:
            for name, by_mm in r.site_counts.items():
                per_primer[name] = per_primer.get(name, 0) + sum(by_mm.values())
        capped = [n for n, v in per_primer.items() if v >= bd.MAX_SITES_PER_PRIMER]
        if capped:
            warn.append(f"the search stops at {bd.MAX_SITES_PER_PRIMER:,} sites per primer "
                        f"({', '.join(html.escape(x) for x in capped[:3])}): real numbers are higher")
        if any(r.sites_truncated for r in self._records):
            warn.append(f"only the first {pcr.MAX_SITES_PER_RECORD:,} sites per sequence are drawn")
        if warn:
            lines.append(f'<span style="color:{WARN_COLOR}">Note: ' + "; ".join(warn) + ".</span>")
        self.summary.setText("<br>".join(lines))
        self.summary.setVisible(bool(self._records))
        self.gff_button.setEnabled(bool(total))
        drawn = min(len(self._records), MAX_STRIPS)
        n = len(self._records)
        self.status.setText(f"{n:,} sequence{'s' if n != 1 else ''} searched, {drawn} drawn. "
                            "Scroll or pinch to zoom, drag to move, double-click to see all.")
        self._sync_view_label()

    def _sync_view_label(self) -> None:
        lo, hi = self.strips.view()
        if not self._records:
            self.view_label.setText("")
            return
        self.view_label.setText(f"Showing {max(0, int(lo)):,} to {int(hi):,} bp "
                                f"({int(hi - lo):,} bp wide)")

    def _on_row(self) -> None:
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if not rows:
            return
        item = self.table.item(rows[0], 0)
        if item is None:
            return
        record = item.text()
        self.strips.fit_sites(record)
        r = self.strips.row_of(record)
        if r >= 0:
            self.scroll.ensureVisible(0, RULER_H + r * ROW_H + ROW_H // 2, 0, ROW_H)

    # ---------------------------------------------------------------- save

    def _save_gff(self) -> None:
        genome = self.genome_row.path()
        start = str(Path(genome).with_name(Path(genome).name + "_primer_sites.gff3")) if genome else ""
        fn, _ = QFileDialog.getSaveFileName(self, "Save Sites as GFF3", start, "GFF3 (*.gff3)")
        if not fn:
            return
        Path(fn).write_text(pcr.sites_gff3(self._records), encoding="utf-8")
        self.fileWritten.emit(Path(fn))

    def _save_png(self) -> None:
        from repeatviewer.widgets.gel_view import save_view_image

        genome = self.genome_row.path()
        start = str(Path(genome).with_name(Path(genome).name + "_primer_sites.png")) if genome else ""
        fn, _ = QFileDialog.getSaveFileName(self, "Save Figure", start, "PNG (*.png)")
        if not fn:
            return
        save_view_image(self.strips, Path(fn), "png")
        self.fileWritten.emit(Path(fn))
