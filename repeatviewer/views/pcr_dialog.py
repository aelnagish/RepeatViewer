"""In-silico PCR window (Tools > Primers and Probes > PCR).

Left: what to run (templates, primers, kind of PCR, conditions) with the
Cancel / Save / Run buttons. Middle: the products on the template and the
products table. Right: a quick gel and the selected product.
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
    "size are labelled 1×, 2×, 3× (tandem repeat monomer). <b>Bisulfite / MSP</b>: templates "
    "converted as methylated and unmethylated DNA (plant: CG and CHG stay C when methylated). "
    "<b>TE present / absent</b>: one lane per accession.<br>"
    "<b>Mismatches per primer</b> and <b>3′ end exact</b> decide which sites count. "
    "<b>Expected yield</b> is relative: fewer mismatches, shorter products and repeated "
    "sites make more DNA.<br>"
    "Products are saved next to the first template (<b>…_PCR_….fa</b>, with the primers and "
    "the template's annotations); the table and the run's history go into <b>…_details</b>."
)

MODE_NOTES = {
    pcr.MODE_STANDARD: "Every product the primers can make, with the mismatches per primer.",
    pcr.MODE_LADDER: "When the primers sit in a tandem repeat, copies give 1×, 2×, 3× … products. "
                     "Brightness falls with size.",
    pcr.MODE_BISULFITE: "Templates converted in silico, once methylated and once unmethylated, "
                        "on both strands. Use primers designed for the converted strand.",
    pcr.MODE_PRESENCE: "One lane per template (accession): a band at the expected size tells "
                       "whether the insertion is there.",
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


class ProductMap(QWidget):
    """Products of the selected lane drawn on their template stretch."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._products: list[pcr.PcrProduct] = []
        self._selected: pcr.PcrProduct | None = None
        self._title = ""
        self.setMinimumHeight(150)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

    def set_products(self, products: list[pcr.PcrProduct], selected: pcr.PcrProduct | None) -> None:
        real = [p for p in products if p.kind != "dimer" and p.record]
        rec = selected.record if selected is not None and selected.record else (
            real[0].record if real else "")
        self._products = [p for p in real if p.record == rec][:12]
        self._selected = selected
        self._title = rec
        rows = len(self._products)
        self.setFixedHeight(70 + max(1, rows) * 16)
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        text = self.palette().color(self.foregroundRole())
        muted = QColor(text)
        muted.setAlpha(150)
        if not self._products:
            p.setPen(muted)
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No products to show")
            p.end()
            return
        lo = min(x.start for x in self._products)
        hi = max(x.end for x in self._products)
        span = max(1, hi - lo)
        lo -= int(span * 0.04)
        hi += int(span * 0.04)
        span = max(1, hi - lo)
        left, right = 12.0, self.width() - 12.0
        w = right - left

        def X(pos: int) -> float:
            return left + (pos - lo) / span * w

        p.setPen(muted)
        p.drawText(QRectF(left, 2, w, 16), Qt.AlignmentFlag.AlignLeft,
                   f"{self._title} · {lo:,}..{hi:,}")
        p.setPen(QPen(text, 2))
        p.drawLine(int(left), 34, int(right), 34)
        # primer sites: starts (+) and ends (-) of the products
        starts = sorted({x.start for x in self._products})
        ends = sorted({x.end for x in self._products})
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor(FORWARD_COLOR))
        for s in starts:
            x = X(s)
            p.drawPolygon([self._pt(x, 22), self._pt(x + 9, 27), self._pt(x, 32)])
        p.setBrush(QColor(REVERSE_COLOR))
        for e in ends:
            x = X(e)
            p.drawPolygon([self._pt(x, 36), self._pt(x - 9, 41), self._pt(x, 46)])
        y = 56
        for prod in self._products:
            x1, x2 = X(prod.start), X(prod.end)
            c = QColor("#3d6ea8")
            c.setAlphaF(max(0.25, min(1.0, prod.yield_ + 0.2)))
            if prod.kind == "off-target":
                c = QColor("#e0a12b")
            p.setBrush(c)
            p.setPen(Qt.PenStyle.NoPen)
            p.drawRoundedRect(QRectF(x1, y, max(3.0, x2 - x1), 9), 2, 2)
            if prod is self._selected:
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.setPen(QPen(QColor("#7fb0ff"), 2))
                p.drawRoundedRect(QRectF(x1 - 2, y - 2, max(3.0, x2 - x1) + 4, 13), 3, 3)
            p.setPen(text)
            label = f"{prod.size:,} bp" + (f" ({prod.label})" if prod.label else "")
            if prod.sites > 1:
                label += f" ×{prod.sites} sites"
            p.drawText(QRectF(x2 + 6, y - 4, 220, 16), Qt.AlignmentFlag.AlignLeft, label)
            y += 16
        p.end()

    @staticmethod
    def _pt(x: float, y: float):
        from PySide6.QtCore import QPointF

        return QPointF(x, y)


class PcrWindow(QDialog):
    runRequested = Signal(object, object, object)     # templates, primers, PcrSettings
    saveRequested = Signal(object)                     # PcrSettings
    openGelRequested = Signal(object)                  # PcrResult
    openFileRequested = Signal(object)                 # Path
    benchCardRequested = Signal(object)                # PcrResult

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
        hint = QLabel("Click FASTA files in the Files panel to add.", page)
        hint.setStyleSheet("color: palette(placeholder-text);")
        hint.setWordWrap(True)
        trow.addWidget(hint, 1)
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
        col1, col2 = QVBoxLayout(), QVBoxLayout()
        for i, mode in enumerate(pcr.MODES):
            b = QPushButton(pcr.MODE_LABELS[mode], page)
            b.setCheckable(True)
            b.setProperty("mode", mode)
            self.mode_group.addButton(b, i)
            (col1 if i % 2 == 0 else col2).addWidget(b)
        grid.addLayout(col1)
        grid.addLayout(col2)
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
        self.meth_combo = QComboBox(page)
        self.meth_combo.addItem("Plant (CG, CHG)", "plant")
        self.meth_combo.addItem("Animal (CG)", "animal")
        for w in (self.mm_combo, self.three_combo, self.max_combo, self.pol_combo, self.anneal,
                  self.meth_combo):
            w.setFixedWidth(150)
        form.addRow("Mismatches per primer", self.mm_combo)
        form.addRow("3′ end exact", self.three_combo)
        form.addRow("Longest product", self.max_combo)
        form.addRow("Polymerase", self.pol_combo)
        form.addRow("Annealing", self.anneal)
        self.meth_label = QLabel("Methylation", page)
        form.addRow(self.meth_label, self.meth_combo)
        pl.addLayout(form)
        self.mode_note = QLabel(page)
        self.mode_note.setWordWrap(True)
        self.mode_note.setStyleSheet("QLabel { border: 1px solid palette(mid); border-radius: 8px;"
                                     " padding: 8px; }")
        pl.addWidget(self.mode_note)
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
        ml.addLayout(lane_row)
        self.map = ProductMap(mid)
        ml.addWidget(self.map)
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
        self.gel_button = QPushButton("Open in a Gel", mid)
        self.products_button = QPushButton("Open Products File", mid)
        self.table_button = QPushButton("Show Table (TSV)", mid)
        self.bench_button = QPushButton("Bench Card…", mid)
        for b in (self.gel_button, self.products_button, self.table_button, self.bench_button):
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
        rl.addWidget(self.product_box)
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
        self.gel_button.clicked.connect(lambda: self._result and self.openGelRequested.emit(self._result))
        self.products_button.clicked.connect(self._open_products)
        self.table_button.clicked.connect(self._open_table)
        self.bench_button.clicked.connect(lambda: self._result and self.benchCardRequested.emit(self._result))
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
        self.meth_combo.setCurrentIndex(max(0, self.meth_combo.findData(s.methylation)))
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
                               methylation=str(self.meth_combo.currentData()),
                               dimers=self.dimer_check.isChecked())

    def right_click_runs(self) -> bool:
        return self.rc_run.isChecked()

    def _sync_mode(self) -> None:
        mode = self.settings().mode
        self.mode_note.setText(MODE_NOTES.get(mode, ""))
        bis = mode == pcr.MODE_BISULFITE
        self.meth_combo.setVisible(bis)
        self.meth_label.setVisible(bis)

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
        for b in (self.gel_button, self.bench_button):
            b.setEnabled(has)
        self.products_button.setEnabled(bool(self._outputs.get("fasta")))
        self.table_button.setEnabled(bool(self._outputs.get("details")))
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
        sub = [f"{len(result.templates)} template(s)", ", ".join(p.name for p in result.primers),
               pcr.MODE_LABELS.get(result.settings.mode, "")]
        if self._outputs.get("fasta"):
            sub.append(f"saved as {Path(self._outputs['fasta']).name}")
        if result.warnings:
            sub.append("; ".join(result.warnings[:3]))
        self.sub.setText(" · ".join(x for x in sub if x))
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
            lines.append(f'<span style="font-family:monospace;color:gray">{snippet}</span>')
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

    def _open_products(self) -> None:
        f = self._outputs.get("fasta")
        if f:
            self.openFileRequested.emit(Path(f))

    def _open_table(self) -> None:
        d = self._outputs.get("details")
        base = self._outputs.get("base")
        if d and base:
            self.openFileRequested.emit(Path(d) / f"{base}_products.tsv")
