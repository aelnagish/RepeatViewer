"""Bench Card: mastermix, PCR program, primers with their freezer place, expected gel.

The card is a page (white, A4 proportions) that can be printed or saved as PDF;
the options on the right change it live.
"""

from __future__ import annotations

import datetime as _dt
import html
from pathlib import Path

from PySide6.QtCore import QPointF, QRectF, Qt
from PySide6.QtGui import QColor, QPainter, QPen
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLayout,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import bench
from repeatviewer.backends.primers import gel_doc as gd
from repeatviewer.backends.primers.oligos import Primer
from repeatviewer.widgets.gel_view import GelView, save_view_image
from repeatviewer.widgets.help_button import make_help_button

BENCH_HELP = (
    "<b>Mastermix</b>: volumes per reaction and for all reactions plus the extra, without the "
    "template (added to each tube).<br>"
    "<b>Program</b>: annealing from the lower primer Tm (Taq: 3 °C below; Phusion / Q5: 3 °C "
    "above), extension from the longest product you want (Taq 1 min per kb, Phusion / Q5 "
    "20 s per kb).<br>"
    "<b>Freezer places</b> come from the Primer Database.<br>"
    "Change the numbers on the right; <b>Save PDF</b> or <b>Print</b> the card."
)


class ProgramChart(QWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._steps: list[bench.Step] = []
        self._cycles = 30
        self.setFixedHeight(170)

    def set_program(self, steps: list[bench.Step], cycles: int) -> None:
        self._steps = steps
        self._cycles = cycles
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        if not self._steps:
            return
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.fillRect(self.rect(), QColor("#ffffff"))
        n = len(self._steps)
        w = (self.width() - 20) / n
        top, bottom = 44.0, self.height() - 12.0

        def Y(t: float) -> float:
            return bottom - (t - 4) / (100 - 4) * (bottom - top)

        cyc = [i for i, s in enumerate(self._steps) if s.in_cycle]
        if cyc:
            x1 = 10 + cyc[0] * w
            x2 = 10 + (cyc[-1] + 1) * w
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor("#f3f7fd"))
            p.drawRoundedRect(QRectF(x1, 2, x2 - x1, bottom), 6, 6)
            p.setPen(QColor("#1f4e9a"))
            p.drawText(QRectF(x1, 4, x2 - x1, 16), Qt.AlignmentFlag.AlignCenter,
                       f"× {self._cycles} cycles")
        pen = QPen(QColor("#c2410c"), 2.5)
        prev = None
        for i, s in enumerate(self._steps):
            x1 = 10 + i * w
            y = Y(s.temp_c)
            p.setPen(pen)
            if prev is not None:
                p.drawLine(QPointF(x1, prev), QPointF(x1, y))
            p.drawLine(QPointF(x1, y), QPointF(x1 + w, y))
            prev = y
            p.setPen(QColor("#1d1d1f"))
            label = f"{s.temp_c:g} °C · {bench.fmt_seconds(s.seconds)}"
            p.drawText(QRectF(x1 + 4, y - 19, w - 6, 16), Qt.AlignmentFlag.AlignLeft, label)
        p.end()


class _Page(QWidget):
    """White page holding the card."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setAutoFillBackground(True)
        pal = self.palette()
        pal.setColor(self.backgroundRole(), QColor("#ffffff"))
        pal.setColor(self.foregroundRole(), QColor("#1d1d1f"))
        self.setPalette(pal)
        self.setFixedWidth(760)


class BenchCardDialog(QDialog):
    def __init__(self, primers: list[Primer], *, result=None, database_path: Path | None = None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Bench Card")
        self.setMinimumSize(1180, 820)
        self._primers = list(primers)
        self._result = result
        self._places = self._lookup_places(database_path)

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        holder = QWidget(self)
        holder.setObjectName("benchHolder")
        holder.setStyleSheet("QWidget#benchHolder, QScrollArea, QScrollArea > QWidget > QWidget"
                             " { background: #dadbdf; }")
        hl = QVBoxLayout(holder)
        hl.setContentsMargins(0, 0, 0, 0)
        self.page = _Page(holder)
        pl = QVBoxLayout(self.page)
        # The page grows with its content (never squeezed), the scroll area scrolls.
        pl.setSizeConstraint(QLayout.SizeConstraint.SetMinimumSize)
        pl.setContentsMargins(36, 30, 36, 30)
        pl.setSpacing(14)
        head = QHBoxLayout()
        self.title = QLabel(self.page)
        f = self.title.font()
        f.setPointSize(f.pointSize() + 5)
        f.setBold(True)
        self.title.setFont(f)
        head.addWidget(self.title, 1)
        self.date = QLabel(_dt.date.today().strftime("%d %b %Y"), self.page)
        head.addWidget(self.date)
        pl.addLayout(head)
        self.subtitle = QLabel(self.page)
        self.subtitle.setWordWrap(True)
        pl.addWidget(self.subtitle)
        cols = QHBoxLayout()
        cols.setSpacing(24)
        self.primers_label = QLabel(self.page)
        self.primers_label.setTextFormat(Qt.TextFormat.RichText)
        self.primers_label.setWordWrap(True)
        self.primers_label.setAlignment(Qt.AlignmentFlag.AlignTop)
        cols.addWidget(self.primers_label, 1)
        self.mix_label = QLabel(self.page)
        self.mix_label.setTextFormat(Qt.TextFormat.RichText)
        self.mix_label.setAlignment(Qt.AlignmentFlag.AlignTop)
        cols.addWidget(self.mix_label, 1)
        pl.addLayout(cols)
        pl.addWidget(QLabel("<b>Program</b>", self.page))
        self.chart = ProgramChart(self.page)
        pl.addWidget(self.chart)
        self.why = QLabel(self.page)
        self.why.setWordWrap(True)
        self.why.setStyleSheet("color: #5f6368;")
        pl.addWidget(self.why)
        gel_row = QHBoxLayout()
        self.gel = GelView(self.page, compact=True)
        self.gel.setFixedSize(230, 260)
        gel_row.addWidget(self.gel)
        self.gel_note = QLabel(self.page)
        self.gel_note.setWordWrap(True)
        self.gel_note.setTextFormat(Qt.TextFormat.RichText)
        self.gel_note.setAlignment(Qt.AlignmentFlag.AlignTop)
        gel_row.addWidget(self.gel_note, 1)
        pl.addLayout(gel_row)
        self.labels = QLabel(self.page)
        self.labels.setTextFormat(Qt.TextFormat.RichText)
        pl.addWidget(self.labels)
        pl.addStretch(1)
        scroll = QScrollArea(holder)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        canvas = QWidget(scroll)
        cl = QHBoxLayout(canvas)
        cl.setContentsMargins(0, 0, 0, 0)
        cl.addStretch(1)
        cl.addWidget(self.page, 0, Qt.AlignmentFlag.AlignTop)
        cl.addStretch(1)
        scroll.setWidget(canvas)
        hl.addWidget(scroll)
        root.addWidget(holder, 1)

        side = QWidget(self)
        side.setFixedWidth(320)
        sl = QVBoxLayout(side)
        sl.setContentsMargins(16, 16, 16, 16)
        sl.setSpacing(10)
        sh = QHBoxLayout()
        sh.addWidget(QLabel("<b>Reaction</b>", side))
        sh.addStretch(1)
        sh.addWidget(make_help_button(side, "Bench Card", BENCH_HELP))
        sl.addLayout(sh)
        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        self.n_spin = QSpinBox(side)
        self.n_spin.setRange(1, 384)
        self.n_spin.setValue(10)
        self.extra = QComboBox(side)
        for v in (0, 5, 10, 15, 20):
            self.extra.addItem(f"{v} %", v)
        self.extra.setCurrentIndex(2)
        self.volume = QComboBox(side)
        for v in (10, 20, 25, 50):
            self.volume.addItem(f"{v} µl", v)
        self.volume.setCurrentIndex(2)
        self.pol = QComboBox(side)
        self.pol.addItems(list(bench.RECIPES))
        self.stock = QComboBox(side)
        for v in (5, 10, 20, 100):
            self.stock.addItem(f"{v} µM", v)
        self.stock.setCurrentIndex(1)
        self.cycles = QSpinBox(side)
        self.cycles.setRange(10, 45)
        self.cycles.setValue(30)
        self.anneal = QDoubleSpinBox(side)
        self.anneal.setRange(40, 75)
        self.anneal.setDecimals(0)
        self.anneal.setSuffix(" °C")
        self.longest = QSpinBox(side)
        self.longest.setRange(50, 50000)
        self.longest.setSingleStep(100)
        self.longest.setSuffix(" bp")
        for w in (self.n_spin, self.extra, self.volume, self.pol, self.stock, self.cycles,
                  self.anneal, self.longest):
            w.setFixedWidth(150)
        form.addRow("Reactions", self.n_spin)
        form.addRow("Extra", self.extra)
        form.addRow("Volume", self.volume)
        form.addRow("Polymerase", self.pol)
        form.addRow("Primer stock", self.stock)
        form.addRow("Cycles", self.cycles)
        form.addRow("Annealing", self.anneal)
        form.addRow("Longest product", self.longest)
        sl.addLayout(form)
        sl.addWidget(QLabel("<b>On the card</b>", side))
        self.show_gel = QCheckBox("Expected gel", side)
        self.show_places = QCheckBox("Freezer places", side)
        self.show_labels = QCheckBox("Tube labels", side)
        self.show_seqs = QCheckBox("Primer sequences", side)
        for c in (self.show_gel, self.show_places, self.show_labels, self.show_seqs):
            c.setChecked(True)
            sl.addWidget(c)
        sl.addStretch(1)
        brow = QHBoxLayout()
        self.pdf_button = QPushButton("Save PDF…", side)
        self.print_button = QPushButton("Print…", side)
        self.print_button.setDefault(True)
        brow.addWidget(self.pdf_button)
        brow.addWidget(self.print_button)
        sl.addLayout(brow)
        root.addWidget(side)

        longest = 1000
        if result is not None:
            wanted = [p.size for p in result.products if p.kind in ("expected", "ladder")]
            longest = max(wanted) if wanted else longest
            if result.settings.polymerase in bench.RECIPES:
                self.pol.setCurrentText(result.settings.polymerase)
        self.longest.setValue(int(longest))
        ta, _tms = bench.suggested_annealing(self._primers, self.pol.currentText())
        self.anneal.setValue(ta)
        for w in (self.n_spin, self.cycles, self.longest):
            w.valueChanged.connect(lambda _v: self._render())
        self.anneal.valueChanged.connect(lambda _v: self._render())
        for w in (self.extra, self.volume, self.stock):
            w.currentIndexChanged.connect(lambda _i: self._render())
        self.pol.currentIndexChanged.connect(self._on_pol)
        for c in (self.show_gel, self.show_places, self.show_labels, self.show_seqs):
            c.toggled.connect(lambda _c: self._render())
        self.pdf_button.clicked.connect(self._save_pdf)
        self.print_button.clicked.connect(self._print)
        self._render()

    # ---------------------------------------------------------------- data

    def _lookup_places(self, database_path: Path | None) -> dict[str, str]:
        try:
            from repeatviewer.backends.primers.inventory import open_database
            from repeatviewer.backends.primers.jobs import primer_location

            db = open_database(database_path)
        except Exception:  # noqa: BLE001 - no database
            return {}
        out = {}
        for p in self._primers:
            hit = db.get(p.name) or db.find_sequence(p.sequence)
            if hit is not None:
                place = primer_location(hit)
                status = (hit.extra or {}).get("status", "")
                out[p.name] = " · ".join(x for x in (place, status) if x)
        return out

    def _on_pol(self, _i: int) -> None:
        ta, _ = bench.suggested_annealing(self._primers, self.pol.currentText())
        self.anneal.blockSignals(True)
        self.anneal.setValue(ta)
        self.anneal.blockSignals(False)
        self._render()

    # ---------------------------------------------------------------- render

    def _render(self) -> None:
        names = " + ".join(p.name for p in self._primers[:2]) or "PCR"
        self.title.setText(f"PCR: {html.escape(names)}")
        sub = []
        if self._result is not None:
            from repeatviewer.backends.primers import pcr

            sub.append(", ".join(Path(t).name for t in self._result.templates))
            sub.append(pcr.describe(self._result))
        sub.append(f"{self.n_spin.value()} reactions")
        self.subtitle.setText(" · ".join(sub))
        rows = []
        for p in self._primers:
            line = f'<b>{html.escape(p.name)}</b>'
            if self.show_places.isChecked() and self._places.get(p.name):
                line += f' · {html.escape(self._places[p.name])}'
            if self.show_seqs.isChecked():
                line += (f'<br><span style="font-family:monospace;font-size:small">'
                         f"{html.escape(p.sequence)}</span>")
            rows.append(line)
        self.primers_label.setText("<b>Primers</b><br>" + "<br>".join(rows))
        mix = bench.mastermix(self._primers, reactions=self.n_spin.value(),
                              extra_percent=float(self.extra.currentData()),
                              volume_ul=float(self.volume.currentData()),
                              polymerase=self.pol.currentText(),
                              primer_stock_uM=float(self.stock.currentData()))
        n_total = self.n_spin.value() * (1 + float(self.extra.currentData()) / 100)
        trs = "".join(
            f"<tr><td>{html.escape(r.reagent)}</td><td align='right'>{r.per_reaction:g}</td>"
            f"<td align='right'>{r.total:g}</td><td>☐</td></tr>" for r in mix)
        per = sum(r.per_reaction for r in mix)
        self.mix_label.setText(
            f"<b>Mastermix</b> ({n_total:g} × {self.volume.currentData()} µl)"
            "<table cellspacing='0' cellpadding='3'>"
            "<tr><td style='color:#5f6368'>Reagent</td><td style='color:#5f6368'>1 × µl</td>"
            "<td style='color:#5f6368'>all µl</td><td></td></tr>"
            f"{trs}<tr><td><b>{per:g} µl mix + 1 µl DNA</b></td><td></td>"
            f"<td align='right'><b>{sum(r.total for r in mix):g}</b></td><td></td></tr></table>")
        steps, cycles, why = bench.program(self._primers, longest_bp=self.longest.value(),
                                           polymerase=self.pol.currentText(),
                                           cycles=self.cycles.value(),
                                           annealing_c=float(self.anneal.value()))
        self.chart.set_program(steps, cycles)
        self.why.setText(why)
        show_gel = self.show_gel.isChecked() and self._result is not None
        self.gel.setVisible(show_gel)
        self.gel_note.setVisible(show_gel)
        if show_gel:
            doc = gd.gel_from_pcr(self._result)
            self.gel.set_doc(doc)
            exp = [p for p in self._result.products if p.kind == "expected"]
            off = [p for p in self._result.products if p.kind == "off-target"]
            note = [f"<b>Expected gel</b> · {doc.agarose:.1f} % agarose, {html.escape(doc.ladder)}"]
            if exp:
                note.append(f"Strongest band {exp[0].size:,} bp.")
            exp_sizes = {p.size for p in exp}
            off_sizes = sorted({p.size for p in off} - exp_sizes)
            if off_sizes:
                shown = ", ".join(f"{x:,}" for x in off_sizes[:6])
                more = f" and {len(off_sizes) - 6} more" if len(off_sizes) > 6 else ""
                note.append(f"Possible off-target bands: {shown}{more} bp.")
            note.append("Notes:<br><br>_______________________________")
            self.gel_note.setText("<br>".join(note))
        if self.show_labels.isChecked():
            tag = (self._primers[0].name.split("_")[0] if self._primers else "PCR")[:10]
            today = _dt.date.today().strftime("%d.%m")
            cells = "".join(
                f"<td style='border:1px dashed #9b9ba3;padding:4px 6px;font-size:small'>"
                f"{html.escape(tag)} #{i}<br>{today}</td>"
                for i in range(1, min(self.n_spin.value(), 10) + 1))
            cells += ("<td style='border:1px dashed #9b9ba3;padding:4px 6px;font-size:small'>NTC"
                      f"<br>{today}</td>")
            self.labels.setText(f"<table cellspacing='6'><tr>{cells}</tr></table>")
        self.labels.setVisible(self.show_labels.isChecked())
        self.page.updateGeometry()

    # ---------------------------------------------------------------- output

    def _save_pdf(self) -> None:
        fn, _ = QFileDialog.getSaveFileName(self, "Save Bench Card", "bench_card.pdf", "PDF (*.pdf)")
        if not fn:
            return
        try:
            save_view_image(self.page, Path(fn), "pdf")
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Bench Card", str(exc))

    def _print(self) -> None:
        try:
            from PySide6.QtPrintSupport import QPrintDialog, QPrinter
        except ImportError:
            QMessageBox.information(self, "Bench Card", "Printing is not available; use Save PDF.")
            return
        printer = QPrinter(QPrinter.PrinterMode.HighResolution)
        dlg = QPrintDialog(printer, self)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        p = QPainter(printer)
        rect = printer.pageRect(QPrinter.Unit.DevicePixel)
        scale = min(rect.width() / self.page.width(), rect.height() / self.page.height())
        p.scale(scale, scale)
        self.page.render(p)
        p.end()
