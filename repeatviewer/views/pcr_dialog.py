"""In-silico PCR window (Tools > Primers and Probes > PCR).

Left: what to run (templates, primers, kind of PCR, conditions) with the
Cancel / Save / Run buttons. Middle: the products on the template and the
products table. Right: a quick gel and the selected product.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import html
from pathlib import Path

from PySide6.QtCore import QRectF, Qt, Signal
from PySide6.QtGui import QColor, QPainter, QPen
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QRadioButton,
    QSizePolicy,
    QSplitter,
    QStyledItemDelegate,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import gel_doc as gd
from repeatviewer.backends.primers import pcr
from repeatviewer.backends.primers import thermo as th
from repeatviewer.backends.primers.oligos import Primer, parse_primer_text
from repeatviewer.widgets.gel_view import GelView
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.primers_panel import (
    CHECK_COLOR,
    FORWARD_COLOR,
    REVERSE_COLOR,
    mono_font,
)

PCR_HELP = (
    "<b>Templates</b>: click FASTA files in the Files panel to add them (a genome works too). "
    "With several templates each one is a lane.<br>"
    "<b>Primers</b>: from the Primers tab of a sequence, the Primer Database, or pasted. "
    "All combinations are tested, as in a real tube (also one primer at both ends).<br>"
    "<b>Standard</b>: every product. <b>Satellite ladder</b>: products that step by a common "
    "size are labelled 1×, 2×, 3× (tandem repeat monomer). With several templates each one "
    "is a lane, so primers can be compared across accessions.<br>"
    "<b>Mismatches per primer</b> and <b>3′ end exact</b> decide which sites count. "
    "<b>Expected yield</b> is relative: fewer mismatches, shorter products and repeated "
    "sites make more DNA.<br>"
    "<b>The strip</b> above the table shows the product you click in the table: forward primer "
    "(blue), reverse primer (orange) and the stretch between them. Small arrows on the line "
    "are the other places the primers bind; thin bars are other sites that give the same size. Scroll or pinch to zoom, drag to move, double-click "
    "to go back to the product.<br>"
    "<b>Sites in Genome…</b> shows every place each primer binds in a genome, also sites "
    "without a partner primer (probes, primer background). No PCR run needed.<br>"
    "Products are saved next to the first template (<b>…_PCR_….fa</b>, with the primers and "
    "the template's annotations); the table and the run's history go into <b>…_details</b>."
)

MODE_NOTES = {
    pcr.MODE_STANDARD: "Every product the primers can make, with the mismatches per primer.",
    pcr.MODE_LADDER: "When the primers sit in a tandem repeat, copies give 1×, 2×, 3× … products. "
                     "Brightness falls with size.",
}


class _YieldDelegate(QStyledItemDelegate):
    def paint(self, painter: QPainter, option, index) -> None:  # noqa: N802
        value = index.data(Qt.ItemDataRole.UserRole)
        if value is None:
            super().paint(painter, option, index)
            return
        painter.save()
        r = option.rect.adjusted(4, option.rect.height() // 2 - 3, -4, 0)
        r.setHeight(6)
        w = max(2, int(min(1.0, float(value)) * min(90, r.width())))
        painter.setPen(Qt.PenStyle.NoPen)
        color = index.data(Qt.ItemDataRole.UserRole + 1) or "#3d6ea8"
        painter.setBrush(QColor(color))
        painter.drawRoundedRect(QRectF(r.left(), r.top(), w, 6), 3, 3)
        painter.restore()


def _parse_location(text: str) -> tuple[str, int, int] | None:
    """"rec:1,234-5,678" → (rec, 1234, 5678)."""
    try:
        rec, rng = str(text).rsplit(":", 1)
        lo, hi = rng.replace(",", "").split("-", 1)
        return rec, int(lo), int(hi)
    except ValueError:
        return None


class ProductMap(QWidget):
    """One zoomable strip: the selected product and its two primers on the template.

    Other products of the same sequence are drawn as faint bars for context.
    Scroll or pinch to zoom, drag to move, double-click to go back to the product.
    """

    viewChanged = Signal()
    RULER_H = 24

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._others: list[tuple[int, int, pcr.PcrProduct]] = []
        self._spans: list[tuple[int, int]] = []          # selected product, every site
        self._selected: pcr.PcrProduct | None = None
        self._record = ""
        self._lens: dict[str, int] = {}
        self._lo, self._hi = 0.0, 1.0
        self._drag_x: float | None = None
        self._drag_lo = 0.0
        self.setMinimumHeight(118)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setFixedHeight(118)
        self.setMouseTracking(True)
        self.grabGesture(Qt.GestureType.PinchGesture)
        self.setToolTip("Scroll or pinch to zoom, drag to move, double-click to go back to the product.")

    # ---------------------------------------------------------------- data

    def set_primer_lengths(self, primers: list[Primer]) -> None:
        self._lens = {}
        for pr in primers:
            try:
                self._lens[pr.name] = len(th.binding_part(th.clean_oligo(pr.sequence)))
            except Exception:  # noqa: BLE001
                self._lens[pr.name] = len(pr.sequence)

    def set_products(self, products: list[pcr.PcrProduct], selected: pcr.PcrProduct | None) -> None:
        real = [p for p in products if p.kind != "dimer" and p.record]
        if selected is not None and (selected.kind == "dimer" or not selected.record):
            selected = None
        if selected is None and real:
            selected = real[0]
        self._selected = selected
        self._record = selected.record if selected is not None else ""
        self._spans = []
        if selected is not None:
            for loc in selected.locations:
                got = _parse_location(loc)
                if got and got[0] == self._record:
                    self._spans.append((got[1], got[2]))
            if not self._spans:
                self._spans = [(selected.start, selected.end)]
        self._others = []
        for p in real:
            if p is selected or p.record != self._record:
                continue
            for loc in p.locations or [f"{p.record}:{p.start}-{p.end}"]:
                got = _parse_location(loc)
                if got and got[0] == self._record:
                    self._others.append((got[1], got[2], p))
        self.fit_product()

    # ---------------------------------------------------------------- view

    def _extent(self) -> tuple[float, float]:
        pts = [x for s in self._spans for x in s] + [x for o in self._others for x in o[:2]]
        if not pts:
            return 0.0, 1.0
        lo, hi = min(pts), max(pts)
        pad = max(20.0, (hi - lo) * 0.04)
        return lo - pad, hi + pad

    def set_view(self, lo: float, hi: float) -> None:
        span = max(20.0, hi - lo)
        ext_lo, ext_hi = self._extent()
        room = (ext_hi - ext_lo) * 0.5 + span
        lo = max(ext_lo - room, min(lo, ext_hi + room - span))
        self._lo, self._hi = lo, lo + span
        self.update()
        self.viewChanged.emit()

    def fit_product(self) -> None:
        if not self._spans:
            self._lo, self._hi = self._extent()
            self.update()
            return
        a, b = self._spans[0]
        pad = max(30.0, (b - a) * 0.35)
        self.set_view(a - pad, b + pad)

    def fit_all(self) -> None:
        self.set_view(*self._extent())

    def zoom(self, factor: float, at_x: float | None = None) -> None:
        if factor <= 0:
            return
        left, width = self._plot()
        frac = 0.5 if at_x is None else min(1.0, max(0.0, (at_x - left) / max(1.0, width)))
        span = self._hi - self._lo
        anchor = self._lo + frac * span
        new = max(20.0, span / factor)
        self.set_view(anchor - frac * new, anchor - frac * new + new)

    def pan_pixels(self, dx: float) -> None:
        _left, width = self._plot()
        shift = -dx / max(1.0, width) * (self._hi - self._lo)
        self.set_view(self._lo + shift, self._hi + shift)

    def _plot(self) -> tuple[float, float]:
        return 14.0, max(10.0, self.width() - 28.0)

    def _x(self, bp: float) -> float:
        left, width = self._plot()
        return left + (bp - self._lo) / max(1e-9, self._hi - self._lo) * width

    # ---------------------------------------------------------------- paint

    def paintEvent(self, event) -> None:  # noqa: N802
        from PySide6.QtCore import QPointF

        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        text = self.palette().color(self.foregroundRole())
        muted = QColor(text)
        muted.setAlpha(150)
        fm = p.fontMetrics()
        if self._selected is None:
            p.setPen(muted)
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No products to show")
            p.end()
            return
        left, width = self._plot()
        right = left + width
        span = self._hi - self._lo
        # ruler
        step = _nice_step(span)
        p.setPen(muted)
        p.drawLine(QPointF(left, self.RULER_H - 6), QPointF(right, self.RULER_H - 6))
        v = (int(self._lo // step) + 1) * step
        while v <= self._hi:
            x = self._x(v)
            p.drawLine(QPointF(x, self.RULER_H - 10), QPointF(x, self.RULER_H - 6))
            p.drawText(QRectF(x - 45, 2, 90, 14), Qt.AlignmentFlag.AlignHCenter, _fmt_bp(v, step))
            v += step
        p.setClipRect(QRectF(left, self.RULER_H, width, self.height()))
        sel = self._selected
        line_y = 62.0
        # template
        p.setPen(QPen(muted, 2))
        p.drawLine(QPointF(left, line_y), QPointF(right, line_y))
        # every other place a primer binds (from all products): small ticks on the template
        p.setPen(Qt.PenStyle.NoPen)
        fwd_c = QColor(FORWARD_COLOR)
        rev_c = QColor(REVERSE_COLOR)
        fwd_c.setAlpha(150)
        rev_c.setAlpha(150)
        starts = {a for a, _b, _p in self._others} | {a for a, _b in self._spans[1:]}
        ends = {b for _a, b, _p in self._others} | {b for _a, b in self._spans[1:]}
        p.setBrush(fwd_c)
        for pos in starts:
            if self._lo <= pos <= self._hi:
                x = self._x(pos)
                p.drawPolygon([QPointF(x, line_y - 9), QPointF(x + 7, line_y - 5.5),
                               QPointF(x, line_y - 2)])
        p.setBrush(rev_c)
        for pos in ends:
            if self._lo <= pos <= self._hi:
                x = self._x(pos)
                p.drawPolygon([QPointF(x, line_y + 2), QPointF(x - 7, line_y + 5.5),
                               QPointF(x, line_y + 9)])
        # other sites that give the same size: thin bars, stacked so they don't overlap
        rows_end: list[float] = []
        thin = QColor("#3d6ea8")
        thin.setAlpha(120)
        p.setBrush(thin)
        for a, b in sorted(self._spans[1:]):
            row = next((i for i, e in enumerate(rows_end) if e < a), None)
            if row is None:
                if len(rows_end) >= 3:
                    continue
                rows_end.append(b)
                row = len(rows_end) - 1
            rows_end[row] = b
            if b < self._lo or a > self._hi:
                continue
            xa, xb = self._x(a), self._x(b)
            p.drawRect(QRectF(xa, line_y + 36 + row * 6, max(2.0, xb - xa), 3))
        # the selected product (first site)
        lf = self._lens.get(sel.forward, 20)
        lr = self._lens.get(sel.reverse, 20)
        a, b = self._spans[0]
        xa, xb = self._x(a), self._x(b)
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor("#3d6ea8"))
        p.drawRect(QRectF(xa, line_y + 8, max(2.0, xb - xa), 8))
        fx2 = max(xa + 10, self._x(a + lf))
        rx1 = min(xb - 10, self._x(b - lr))
        p.setBrush(QColor(FORWARD_COLOR))
        p.drawPolygon([QPointF(xa, line_y - 14), QPointF(fx2 - 6, line_y - 14),
                       QPointF(fx2, line_y - 9), QPointF(fx2 - 6, line_y - 4),
                       QPointF(xa, line_y - 4)])
        p.setBrush(QColor(REVERSE_COLOR))
        p.drawPolygon([QPointF(xb, line_y - 14), QPointF(rx1 + 6, line_y - 14),
                       QPointF(rx1, line_y - 9), QPointF(rx1 + 6, line_y - 4),
                       QPointF(xb, line_y - 4)])
        p.setPen(text)
        fname = sel.forward + (f" · {sel.mm_forward} mm" if sel.mm_forward else "")
        rname = sel.reverse + (f" · {sel.mm_reverse} mm" if sel.mm_reverse else "")
        fw = fm.horizontalAdvance(fname) + 4
        rw = fm.horizontalAdvance(rname) + 4
        fxl = max(left, min(xa, right - fw))
        rxl = min(right - rw, max(xb - rw, left))
        p.drawText(QRectF(fxl, line_y - 34, fw, 16), Qt.AlignmentFlag.AlignLeft, fname)
        if rxl < fxl + fw:            # too close: put the reverse name under the bar
            p.drawText(QRectF(rxl, line_y + 18, rw, 16), Qt.AlignmentFlag.AlignLeft, rname)
        else:
            p.drawText(QRectF(rxl, line_y - 34, rw, 16), Qt.AlignmentFlag.AlignLeft, rname)
        lab = f"{sel.size:,} bp" + (f" · {sel.label}" if sel.label else "")
        lw = fm.horizontalAdvance(lab) + 6
        cx = (max(xa, left) + min(xb, right)) / 2
        lx = max(left, min(cx - lw / 2, right - lw))
        if rxl < fxl + fw:
            lx = max(left, min(rxl + rw + 8, right - lw))
        p.drawText(QRectF(lx, line_y + 18, lw, 16), Qt.AlignmentFlag.AlignHCenter, lab)
        p.end()

    # ---------------------------------------------------------------- input

    def wheelEvent(self, event) -> None:  # noqa: N802
        mods = event.modifiers()
        pixel = event.pixelDelta()
        angle = event.angleDelta()
        trackpad = event.phase() != Qt.ScrollPhase.NoScrollPhase or not pixel.isNull()
        x = event.position().x()
        if mods & (Qt.KeyboardModifier.ControlModifier | Qt.KeyboardModifier.MetaModifier):
            self.zoom(1.0015 ** (angle.y() or pixel.y()), x)
        elif trackpad:
            if abs(pixel.x()) > abs(pixel.y()):
                self.pan_pixels(pixel.x())
            else:
                self.zoom(1.003 ** pixel.y(), x)
        elif mods & (Qt.KeyboardModifier.ShiftModifier | Qt.KeyboardModifier.AltModifier):
            self.pan_pixels((angle.y() or angle.x()) / 2)
        else:
            self.zoom(1.25 if angle.y() > 0 else 0.8, x)
        event.accept()

    def event(self, ev) -> bool:  # noqa: D401
        from PySide6.QtCore import QEvent

        if ev.type() == QEvent.Type.NativeGesture:
            if ev.gestureType() == Qt.NativeGestureType.ZoomNativeGesture:
                self.zoom(1.0 + ev.value(), ev.position().x())
                return True
            if ev.gestureType() == Qt.NativeGestureType.SmartZoomNativeGesture:
                self.fit_product()
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
        if self._drag_x is None:
            return
        _left, width = self._plot()
        span = self._hi - self._lo
        lo = self._drag_lo - (event.position().x() - self._drag_x) / max(1.0, width) * span
        self.set_view(lo, lo + span)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        self._drag_x = None
        self.unsetCursor()

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        self.fit_product()


def _nice_step(span: float, target: int = 7) -> float:
    import math

    raw = max(1.0, span / max(1, target))
    mag = 10 ** math.floor(math.log10(raw))
    for m in (1, 2, 5, 10):
        if m * mag >= raw:
            return m * mag
    return 10 * mag


def _fmt_bp(v: float, step: float) -> str:
    if step >= 1_000_000:
        return f"{v / 1e6:g} Mb"
    if step >= 1_000:
        return f"{v / 1e3:g} kb"
    return f"{int(v):,}"


class PcrWindow(QDialog):
    runRequested = Signal(object, object, object)     # templates, primers, PcrSettings
    saveRequested = Signal(object)                     # PcrSettings
    openGelRequested = Signal(object)                  # PcrResult
    openFileRequested = Signal(object)                 # Path
    benchCardRequested = Signal(object)                # PcrResult
    sitesRequested = Signal(object, object)            # primers, genome Path | None

    def __init__(self, settings: pcr.PcrSettings, *, templates: list[Path] | None = None,
                 primers: list[Primer] | None = None, database_path: Path | None = None,
                 right_click_runs: bool = False, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("PCR")
        self.setMinimumSize(1280, 820)
        self._database_path = database_path
        self._primers: list[Primer] = list(primers or [])
        self._templates: list[Path] = [Path(t) for t in templates or []]
        self._result: pcr.PcrResult | None = None
        self._outputs: dict = {}

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        # ---------------- left: inputs
        left = QWidget(self)
        left.setFixedWidth(340)
        ll = QVBoxLayout(left)
        ll.setContentsMargins(14, 12, 14, 12)
        ll.setSpacing(10)
        self.tabs = QTabWidget(left)
        page = QWidget()
        pl = QVBoxLayout(page)
        pl.setContentsMargins(8, 10, 8, 8)
        pl.setSpacing(8)
        th_row = QHBoxLayout()
        th_row.addWidget(QLabel("<b>Templates</b>", page))
        th_row.addStretch(1)
        th_row.addWidget(make_help_button(page, "PCR", PCR_HELP))
        pl.addLayout(th_row)
        self.templates_list = QListWidget(page)
        self.templates_list.setFixedHeight(76)
        pl.addWidget(self.templates_list)
        trow = QHBoxLayout()
        self.templates_list.setToolTip("Click FASTA files in the Files panel to add them.")
        trow.addStretch(1)
        self.remove_template = QPushButton("Remove", page)
        trow.addWidget(self.remove_template)
        pl.addLayout(trow)
        pl.addWidget(QLabel("<b>Primers</b>", page))
        self.primers_list = QListWidget(page)
        self.primers_list.setFixedHeight(76)
        pl.addWidget(self.primers_list)
        prow = QHBoxLayout()
        self.db_button = QPushButton("From Database…", page)
        self.paste_button = QPushButton("Paste…", page)
        self.remove_primer = QPushButton("Remove", page)
        for b in (self.db_button, self.paste_button, self.remove_primer):
            prow.addWidget(b)
        pl.addLayout(prow)
        pl.addWidget(QLabel("<b>Kind of PCR</b>", page))
        grid = QHBoxLayout()
        self.mode_group = QButtonGroup(page)
        self.mode_group.setExclusive(True)
        for i, mode in enumerate(pcr.MODES):
            b = QPushButton(pcr.MODE_LABELS[mode], page)
            b.setCheckable(True)
            b.setProperty("mode", mode)
            b.setToolTip(MODE_NOTES.get(mode, ""))
            self.mode_group.addButton(b, i)
            grid.addWidget(b, 1)
        pl.addLayout(grid)
        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignLeft)
        self.mm_combo = QComboBox(page)
        for n in range(0, 5):
            self.mm_combo.addItem("none" if n == 0 else f"up to {n}", n)
        self.three_combo = QComboBox(page)
        for n in (0, 2, 3, 4, 5):
            self.three_combo.addItem("no rule" if n == 0 else f"last {n} bases", n)
        self.max_combo = QComboBox(page)
        for n in (1000, 2000, 5000, 10000, 20000):
            self.max_combo.addItem(f"{n // 1000} kb", n)
        self.pol_combo = QComboBox(page)
        self.pol_combo.addItems(list(pcr.POLYMERASES))
        self.anneal = QDoubleSpinBox(page)
        self.anneal.setRange(40.0, 75.0)
        self.anneal.setDecimals(1)
        self.anneal.setSuffix(" °C")
        for w in (self.mm_combo, self.three_combo, self.max_combo, self.pol_combo, self.anneal):
            w.setFixedWidth(150)
        form.addRow("Mismatches per primer", self.mm_combo)
        form.addRow("3′ end exact", self.three_combo)
        form.addRow("Longest product", self.max_combo)
        form.addRow("Polymerase", self.pol_combo)
        form.addRow("Annealing", self.anneal)
        pl.addLayout(form)
        self.mode_note = QLabel(page)
        self.mode_note.setWordWrap(True)
        self.mode_note.setStyleSheet("QLabel { border: 1px solid palette(mid); border-radius: 8px;"
                                     " padding: 8px; }")
        self.mode_note.hide()
        pl.addStretch(1)
        self.tabs.addTab(page, "PCR")

        opt = QWidget()
        ol = QVBoxLayout(opt)
        ol.setContentsMargins(8, 10, 8, 8)
        ol.setSpacing(8)
        ol.addWidget(QLabel("<b>Right-click action</b>", opt))
        self.rc_run = QRadioButton("Run PCR with saved settings", opt)
        self.rc_open = QRadioButton("Open this window", opt)
        (self.rc_run if right_click_runs else self.rc_open).setChecked(True)
        ol.addWidget(self.rc_run)
        ol.addWidget(self.rc_open)
        ol.addWidget(QLabel("<b>Products</b>", opt))
        self.dimer_check = QCheckBox("Look for primer-dimer", opt)
        ol.addWidget(self.dimer_check)
        ol.addStretch(1)
        self.tabs.addTab(opt, "Options")
        ll.addWidget(self.tabs, 1)
        brow = QHBoxLayout()
        self.cancel_button = QPushButton("Cancel", left)
        self.save_button = QPushButton("Save", left)
        self.run_button = QPushButton("Run", left)
        self.run_button.setDefault(True)
        for b in (self.cancel_button, self.save_button, self.run_button):
            brow.addWidget(b)
        ll.addLayout(brow)
        root.addWidget(left)
        line = QFrame(self)
        line.setFrameShape(QFrame.Shape.VLine)
        line.setFrameShadow(QFrame.Shadow.Plain)
        line.setStyleSheet("color: palette(mid);")
        root.addWidget(line)

        # ---------------- middle: results
        mid = QWidget(self)
        ml = QVBoxLayout(mid)
        ml.setContentsMargins(18, 14, 18, 14)
        ml.setSpacing(10)
        self.summary = QLabel(mid)
        self.summary.setTextFormat(Qt.TextFormat.RichText)
        f = self.summary.font()
        f.setPointSize(f.pointSize() + 3)
        self.summary.setFont(f)
        ml.addWidget(self.summary)
        self.sub = QLabel(mid)
        self.sub.setStyleSheet("color: palette(placeholder-text);")
        self.sub.setWordWrap(True)
        ml.addWidget(self.sub)
        lane_row = QHBoxLayout()
        self.lane_label = QLabel("Lane", mid)
        lane_row.addWidget(self.lane_label)
        self.lane_combo = QComboBox(mid)
        self.lane_combo.setMinimumWidth(220)
        lane_row.addWidget(self.lane_combo)
        lane_row.addStretch(1)
        self.fit_product_button = QPushButton("Zoom to Product", mid)
        self.fit_all_button = QPushButton("Show All Products", mid)
        lane_row.addWidget(self.fit_product_button)
        lane_row.addWidget(self.fit_all_button)
        ml.addLayout(lane_row)
        self.map = ProductMap(mid)
        ml.addWidget(self.map)
        self.fit_product_button.clicked.connect(self.map.fit_product)
        self.fit_all_button.clicked.connect(self.map.fit_all)
        cols = ["Product", "Size", "Where", "Mism. F / R", "Sites", "Expected yield", "Kind"]
        self.table = QTableWidget(0, len(cols), mid)
        self.table.setHorizontalHeaderLabels(cols)
        self.table.verticalHeader().hide()
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        h = self.table.horizontalHeader()
        h.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        h.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        h.setStretchLastSection(True)
        for c, w in enumerate((190, 80, 220, 90, 60, 110)):
            self.table.setColumnWidth(c, w)
        self.table.setItemDelegateForColumn(5, _YieldDelegate(self.table))
        ml.addWidget(self.table, 1)
        actions = QHBoxLayout()
        self.gel_button = QPushButton("Gel", mid)
        self.gel_button.setToolTip("Open the products in a gel document you can edit and save")
        self.products_button = QPushButton("Products File", mid)
        self.products_button.setToolTip("Open the FASTA file with the product sequences")
        self.table_button = QPushButton("Table", mid)
        self.table_button.setToolTip("Open the products table (TSV)")
        self.bench_button = QPushButton("Bench Card", mid)
        self.bench_button.setToolTip("Mastermix, PCR program and expected gel for the bench")
        self.sites_button = QPushButton("Sites in Genome", mid)
        self.sites_button.setToolTip(
            "Every place each primer binds in a genome (by default the first template), "
            "per chromosome. Useful for probes and single primers: no partner primer needed.")
        for b in (self.gel_button, self.products_button, self.table_button, self.bench_button,
                  self.sites_button):
            b.setMinimumWidth(b.sizeHint().width())     # never cut the label (Mac fonts are wider)
            actions.addWidget(b)
        actions.addStretch(1)
        ml.addLayout(actions)
        root.addWidget(mid, 1)
        line2 = QFrame(self)
        line2.setFrameShape(QFrame.Shape.VLine)
        line2.setStyleSheet("color: palette(mid);")
        root.addWidget(line2)

        # ---------------- right: quick gel + product
        right = QWidget(self)
        right.setFixedWidth(310)
        rl = QVBoxLayout(right)
        rl.setContentsMargins(14, 14, 14, 14)
        rl.setSpacing(10)
        self.gel_title = QLabel("<b>Quick gel</b>", right)
        rl.addWidget(self.gel_title)
        self.gel = GelView(right, compact=True)
        self.gel.setFixedHeight(380)
        rl.addWidget(self.gel)
        self.product_box = QLabel(right)
        self.product_box.setWordWrap(True)
        self.product_box.setTextFormat(Qt.TextFormat.RichText)
        self.product_box.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.product_box.setStyleSheet("QLabel { border: 1px solid palette(mid); border-radius: 8px;"
                                       " padding: 8px; }")
        self.product_box.hide()
        self.copy_button = QPushButton("Copy Product Sequence", right)
        rl.addWidget(self.copy_button, 0, Qt.AlignmentFlag.AlignLeft)
        rl.addStretch(1)
        root.addWidget(right)

        # ---------------- wiring
        self.mode_group.idClicked.connect(lambda _i: self._sync_mode())
        self.remove_template.clicked.connect(self._remove_template)
        self.remove_primer.clicked.connect(self._remove_primer)
        self.db_button.clicked.connect(self._from_database)
        self.paste_button.clicked.connect(self._paste)
        self.cancel_button.clicked.connect(self.reject)
        self.save_button.clicked.connect(lambda: self.saveRequested.emit(self.settings()))
        self.run_button.clicked.connect(self._run)
        self.lane_combo.currentIndexChanged.connect(lambda _i: self._fill_table())
        self.table.itemSelectionChanged.connect(self._on_row)
        self._pending_action = ""
        self.gel_button.clicked.connect(lambda: self._after_run("gel"))
        self.products_button.clicked.connect(lambda: self._after_run("products"))
        self.table_button.clicked.connect(lambda: self._after_run("table"))
        self.bench_button.clicked.connect(lambda: self._after_run("bench"))
        self.sites_button.clicked.connect(self._open_sites)
        self.copy_button.clicked.connect(self._copy_product)
        self.gel.laneClicked.connect(self._on_gel_lane)

        self.set_settings(settings)
        self._fill_inputs()
        self.show_result(None)

    # ---------------------------------------------------------------- settings

    def set_settings(self, s: pcr.PcrSettings) -> None:
        idx = pcr.MODES.index(s.mode) if s.mode in pcr.MODES else 0
        self.mode_group.button(idx).setChecked(True)
        self.mm_combo.setCurrentIndex(max(0, self.mm_combo.findData(s.max_mismatches)))
        self.three_combo.setCurrentIndex(max(0, self.three_combo.findData(s.three_prime_exact)))
        i = self.max_combo.findData(s.max_product)
        self.max_combo.setCurrentIndex(i if i >= 0 else 2)
        self.pol_combo.setCurrentText(s.polymerase)
        self.anneal.setValue(s.annealing_c)
        self.dimer_check.setChecked(s.dimers)
        self._sync_mode()

    def settings(self) -> pcr.PcrSettings:
        b = self.mode_group.checkedButton()
        mode = b.property("mode") if b is not None else pcr.MODE_STANDARD
        return pcr.PcrSettings(mode=mode, max_mismatches=int(self.mm_combo.currentData()),
                               three_prime_exact=int(self.three_combo.currentData()),
                               max_product=int(self.max_combo.currentData()),
                               polymerase=self.pol_combo.currentText(),
                               annealing_c=float(self.anneal.value()),
                               dimers=self.dimer_check.isChecked())

    def right_click_runs(self) -> bool:
        return self.rc_run.isChecked()

    def _sync_mode(self) -> None:
        mode = self.settings().mode
        self.mode_note.setText(MODE_NOTES.get(mode, ""))

    # ---------------------------------------------------------------- inputs

    def templates(self) -> list[Path]:
        return list(self._templates)

    def primers(self) -> list[Primer]:
        return list(self._primers)

    def _fill_inputs(self) -> None:
        self.templates_list.clear()
        for t in self._templates:
            item = QListWidgetItem(t.name)
            item.setToolTip(str(t))
            self.templates_list.addItem(item)
        self.primers_list.clear()
        for p in self._primers:
            item = QListWidgetItem(f"{p.name}  ·  {len(p.sequence)}-mer")
            item.setToolTip(p.sequence)
            self.primers_list.addItem(item)
        self.run_button.setEnabled(bool(self._templates) and bool(self._primers))

    def apply_path_from_panel(self, path: Path) -> bool:
        from repeatviewer.primer_actions import is_primer_fasta

        path = Path(path)
        if not path.is_file() or not is_primer_fasta(path):
            return False
        if path not in self._templates:
            self._templates.append(path)
            self._fill_inputs()
        return True

    def add_primers(self, primers: list[Primer]) -> None:
        names = {p.name for p in self._primers}
        for p in primers:
            if p.name not in names:
                self._primers.append(p)
                names.add(p.name)
        self._fill_inputs()

    def _remove_template(self) -> None:
        rows = sorted({i.row() for i in self.templates_list.selectedIndexes()}, reverse=True)
        for r in rows:
            del self._templates[r]
        self._fill_inputs()

    def _remove_primer(self) -> None:
        rows = sorted({i.row() for i in self.primers_list.selectedIndexes()}, reverse=True)
        for r in rows:
            del self._primers[r]
        self._fill_inputs()

    def _from_database(self) -> None:
        from repeatviewer.views.primer_edit_dialog import PickPrimersDialog

        dlg = PickPrimersDialog(self._database_path, parent=self)
        if dlg.exec():
            self.add_primers(dlg.primers())

    def _paste(self) -> None:
        text, ok = QInputDialog.getMultiLineText(
            self, "Paste Primers", "FASTA, or one primer per line (name and sequence):")
        if not ok or not text.strip():
            return
        try:
            self.add_primers(parse_primer_text(text))
        except th.PrimerError as exc:
            from PySide6.QtWidgets import QMessageBox

            QMessageBox.warning(self, "PCR", str(exc))

    def _run(self) -> None:
        if not self._templates or not self._primers:
            return
        self.summary.setText('<span style="color:gray">Running PCR…</span>')
        self.runRequested.emit(self.templates(), self.primers(), self.settings())

    # ---------------------------------------------------------------- results

    def show_result(self, result: pcr.PcrResult | None, outputs: dict | None = None) -> None:
        self._result = result
        self._outputs = dict(outputs or {})
        has = result is not None
        pending, self._pending_action = self._pending_action, ""
        if has and pending:
            from PySide6.QtCore import QTimer

            QTimer.singleShot(0, lambda a=pending: self._do_action(a))
        self.copy_button.setEnabled(False)
        self.lane_combo.blockSignals(True)
        self.lane_combo.clear()
        if not has:
            self.lane_label.hide()
            self.lane_combo.hide()
            self.summary.setText("<b>No run yet</b>")
            self.sub.setText("Add templates and primers, choose the kind of PCR, then Run.")
            self.lane_combo.blockSignals(False)
            self.table.setRowCount(0)
            self.map.set_products([], None)
            self.gel.set_doc(None)
            self.product_box.setText("")
            return
        for lane in result.lanes:
            self.lane_combo.addItem(lane)
        self.lane_combo.blockSignals(False)
        self.lane_combo.setVisible(len(result.lanes) > 1)
        self.lane_label.setVisible(len(result.lanes) > 1)
        self.summary.setText(f"<b>{html.escape(pcr.describe(result))}</b>")
        names = [p.name for p in result.primers]
        prim = ", ".join(names) if len(names) <= 3 else f"{len(names)} primers"
        sub = [f"{len(result.templates)} template(s)", prim,
               pcr.MODE_LABELS.get(result.settings.mode, "")]
        if self._outputs.get("fasta"):
            sub.append(f"saved as {Path(self._outputs['fasta']).name}")
        if result.warnings:
            sub.append("; ".join(result.warnings[:3]))
        self.sub.setText(" · ".join(x for x in sub if x))
        self.map.set_primer_lengths(result.primers)
        doc = gd.gel_from_pcr(result)
        self.gel.set_doc(doc)
        self.gel_title.setText(f"<b>Quick gel</b> · {doc.agarose:.1f} % agarose")
        self._fill_table()

    def _lane_products(self) -> list[pcr.PcrProduct]:
        if self._result is None:
            return []
        lane = self.lane_combo.currentText() or (self._result.lanes[0] if self._result.lanes else "")
        return self._result.lane_products(lane)

    def _fill_table(self) -> None:
        prods = self._lane_products()
        self.table.setRowCount(len(prods))
        for r, p in enumerate(prods):
            name = f"{p.label} " if p.label and p.kind != "dimer" else ""
            name += f"{p.forward} + {p.reverse}" if p.kind != "dimer" else p.label
            vals = [name, f"{p.size:,} bp",
                    p.locations[0] if p.locations else "", f"{p.mm_forward} / {p.mm_reverse}"
                    if p.kind != "dimer" else "–", str(p.sites), "", p.kind]
            for c, v in enumerate(vals):
                item = QTableWidgetItem(v)
                if c == 5:
                    item.setData(Qt.ItemDataRole.UserRole, p.yield_)
                    item.setData(Qt.ItemDataRole.UserRole + 1,
                                 "#e0a12b" if p.kind == "off-target" else "#3d6ea8")
                    item.setToolTip(f"{p.yield_:.2f} (relative)")
                if c == 6 and p.kind == "off-target":
                    item.setForeground(QColor(CHECK_COLOR))
                if c == 2 and len(p.locations) > 1:
                    item.setToolTip("\n".join(p.locations))
                self.table.setItem(r, c, item)
        if prods:
            self.table.selectRow(0)
        else:
            self.map.set_products([], None)

    def _on_row(self) -> None:
        prods = self._lane_products()
        rows = sorted({i.row() for i in self.table.selectedIndexes()})
        if not rows or rows[0] >= len(prods):
            return
        p = prods[rows[0]]
        self.map.set_products(prods, p)
        lane_idx = self._result.lanes.index(p.lane) + 1 if p.lane in self._result.lanes else 1
        self.gel.select_band(lane_idx, p.size)
        gc = th.gc_percent(p.sequence) if p.sequence else 0.0
        lines = [f"<b>{html.escape(p.label or p.kind)} · {p.size:,} bp</b>"]
        if p.sequence:
            lines.append(f"{gc:.0f} % GC · {p.forward} → … ← {p.reverse}")
            snippet = html.escape(p.sequence[:44]) + ("…" if len(p.sequence) > 44 else "")
            lines.append(f'<span style="font-family:Menlo,Consolas,Courier;color:gray">{snippet}</span>')
        if p.locations:
            lines.append(html.escape(", ".join(p.locations[:3])))
        if p.sites > 1:
            lines.append(f"{p.sites} sites give this size")
        self.product_box.setText("<br>".join(lines))
        self.copy_button.setEnabled(bool(p.sequence))
        self._current_product = p

    def _on_gel_lane(self, lane: int) -> None:
        if self._result is None:
            return
        i = lane - 1
        if 0 <= i < self.lane_combo.count():
            self.lane_combo.setCurrentIndex(i)

    def _copy_product(self) -> None:
        p = getattr(self, "_current_product", None)
        if p is not None and p.sequence:
            QApplication.clipboard().setText(f">{p.forward}+{p.reverse}_{p.size}bp\n{p.sequence}\n")

    def _missing(self, need_templates: bool = True) -> bool:
        miss = []
        if need_templates and not self._templates:
            miss.append("a template (click a FASTA file in the Files panel)")
        if not self._primers:
            miss.append("primers")
        if miss:
            self.summary.setText(f"<b>Add {' and '.join(miss)} first</b>")
            return True
        return False

    def _after_run(self, action: str) -> None:
        """Gel / Products / Table / Bench Card: run the PCR first when there is no result."""
        if self._result is not None:
            self._do_action(action)
            return
        if self._missing():
            return
        self._pending_action = action
        self._run()

    def _do_action(self, action: str) -> None:
        if self._result is None:
            return
        if action == "gel":
            self.openGelRequested.emit(self._result)
        elif action == "bench":
            self.benchCardRequested.emit(self._result)
        elif action == "products":
            self._open_products()
        elif action == "table":
            self._open_table()

    def _open_sites(self) -> None:
        if self._missing(need_templates=False):
            return
        self.sitesRequested.emit(list(self._primers),
                                 self._templates[0] if self._templates else None)

    def _open_products(self) -> None:
        f = self._outputs.get("fasta")
        if f:
            self.openFileRequested.emit(Path(f))

    def _open_table(self) -> None:
        d = self._outputs.get("details")
        base = self._outputs.get("base")
        if d and base:
            self.openFileRequested.emit(Path(d) / f"{base}_products.tsv")
