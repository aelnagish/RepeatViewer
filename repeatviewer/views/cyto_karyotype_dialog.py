"""Edit karyotype: chromosomes of a species (number, lengths, centromere index, copies, bands)."""

from __future__ import annotations

import copy
from pathlib import Path

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QFont, QPainter, QPen
from PySide6.QtWidgets import (
    QApplication,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.cytogenetics import karyotype as kt
from repeatviewer.backends.cytogenetics.lists import load_lists
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.help_button import make_help_button

HELP = (
    "<b>Karyotype</b><br>"
    "One set of chromosomes for a species: basic number x, ploidy and, per chromosome, length, centromere index, "
    "copies and bands.<br>"
    "<b>Centromere index</b> is the short arm as percent of the whole chromosome; the type (m, sm, st, t) follows Levan et al. (1964).<br>"
    "<b>Copies</b> can differ from the ploidy (aneuploidy, B chromosomes). Leave it at the ploidy for a normal set.<br>"
    "<b>Order</b>: drag the ⋮⋮ handle to move a chromosome, or use the Order buttons.<br>"
    "<b>Import measurements</b> reads a table with chromosome and length (and optionally centromere index or p and q), "
    "for example from KaryoScope.<br>"
    "<b>Take lengths from assembly</b> links Chr1, Chr01 … to chromosome 1 and uses the sequence lengths in Mb.<br>"
    "<b>Preview</b>: drag a centromere or a band edge to change it.<br>"
    "<b>Bands</b> are positions like \"q 0.92\": 92 % of the way from the centromere to the q telomere."
)

UNITS = (("µm", "µm"), ("% of genome (relative)", "%"), ("Mb (assembly)", "Mb"))
PLOIDIES = (1, 2, 3, 4, 5, 6, 8)
BAND_COLOURS = {
    "DAPI+ heterochromatin": "#5B74FF",
    "CMA+ band": "#E8A33D",
}
CUSTOM_KIND = "Custom…"
LINK_CSS = "QPushButton { color:#1F4E9A; border:none; } QPushButton:hover { text-decoration: underline; }"


class _Edit(QLineEdit):
    focused = Signal()

    def focusInEvent(self, event) -> None:  # noqa: N802
        self.focused.emit()
        super().focusInEvent(event)


class _ReorderTable(QTableWidget):
    rowMoved = Signal(int, int)

    def __init__(self, *args) -> None:
        super().__init__(*args)
        self.setDragEnabled(True)
        self.setAcceptDrops(True)
        self.viewport().setAcceptDrops(True)
        self.setDragDropMode(QTableWidget.DragDropMode.InternalMove)
        self.setDropIndicatorShown(True)
        self.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)

    def dropEvent(self, event) -> None:  # noqa: N802
        src = self.currentRow()
        dst = self.rowAt(int(event.position().y()))
        if dst < 0:
            dst = self.rowCount() - 1
        event.setDropAction(Qt.DropAction.IgnoreAction)
        event.ignore()
        if src >= 0 and dst != src:
            self.rowMoved.emit(src, dst)


class KaryoPreview(QWidget):
    """Chromosome bars side by side. Drag a centromere or a band edge."""

    selected = Signal(int)
    edited = Signal()

    BAR_W = 22
    MAX_H = 300

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumSize(300, 340)
        self.karyotype: kt.Karyotype | None = None
        self.index = -1
        self.band_index = -1
        self._drag: tuple[str, int] | None = None

    def _geometry(self, i: int, chrom: kt.KaryoChromosome) -> QRectF:
        k = self.karyotype
        n = max(1, len(k.chromosomes))
        slot = max(30.0, (self.width() - 16) / n)
        longest = max((c.length for c in k.chromosomes), default=0.0)
        h = self.MAX_H * (chrom.length / longest) if longest > 0 and chrom.length > 0 else self.MAX_H * 0.5
        h = max(30.0, min(h, self.height() - 40.0))
        x = 8 + slot * (i + 0.5) - self.BAR_W / 2
        return QRectF(x, self.height() - 26 - h, self.BAR_W, h)

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.fillRect(self.rect(), QColor("#FFFFFF"))
        k = self.karyotype
        if k is None or not k.chromosomes:
            p.setPen(QColor("#5A5D63"))
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No chromosomes")
            return
        font = QFont(self.font())
        font.setBold(True)
        p.setFont(font)
        for i, c in enumerate(k.chromosomes):
            r = self._geometry(i, c)
            if i == self.index:
                p.setPen(QPen(QColor("#1F4E9A"), 1.5, Qt.PenStyle.DashLine))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRoundedRect(r.adjusted(-8, -8, 8, 22), 8, 8)
            for b_i, b in enumerate(c.bands):
                top = r.top() + b.start * r.height()
                rect = QRectF(r.left(), top, r.width(), max(2.0, (b.end - b.start) * r.height()))
                if b.kind.startswith("Secondary"):
                    brush = QBrush(QColor("#8A8D93"), Qt.BrushStyle.BDiagPattern)
                else:
                    brush = QBrush(QColor(BAND_COLOURS.get(b.kind, "#8A8D93")))
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(brush)
                p.drawRect(rect)
            p.setPen(QPen(QColor("#6C6F75"), 1.5))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRoundedRect(r, self.BAR_W / 2, self.BAR_W / 2)
            cy = r.top() + c.centromere_fraction * r.height()
            p.setBrush(QColor("#FFFFFF"))
            p.drawRoundedRect(QRectF(r.left() - 3, cy - 3, r.width() + 6, 6), 3, 3)
            p.setPen(QColor("#1C1C1E"))
            p.drawText(QRectF(r.left() - 15, self.height() - 22, r.width() + 30, 18), Qt.AlignmentFlag.AlignCenter, c.name)
            if i == self.index and 0 <= self.band_index < len(c.bands):
                b = c.bands[self.band_index]
                p.setPen(QPen(QColor("#1F4E9A"), 2))
                for f in (b.start, b.end):
                    y = r.top() + f * r.height()
                    p.drawLine(QPointF(r.left() - 6, y), QPointF(r.right() + 6, y))

    def _hit(self, pos: QPointF) -> tuple[str, int] | None:
        k = self.karyotype
        if k is None:
            return None
        if 0 <= self.index < len(k.chromosomes):
            c = k.chromosomes[self.index]
            r = self._geometry(self.index, c)
            if 0 <= self.band_index < len(c.bands) and r.left() - 10 <= pos.x() <= r.right() + 10:
                b = c.bands[self.band_index]
                for edge, f in (("start", b.start), ("end", b.end)):
                    if abs(pos.y() - (r.top() + f * r.height())) <= 4:
                        return edge, self.index
        for i, c in enumerate(k.chromosomes):
            r = self._geometry(i, c)
            cy = r.top() + c.centromere_fraction * r.height()
            if r.left() - 8 <= pos.x() <= r.right() + 8 and abs(pos.y() - cy) <= 6:
                return "centromere", i
        for i, c in enumerate(k.chromosomes):
            if self._geometry(i, c).adjusted(-8, -8, 8, 8).contains(pos):
                return "select", i
        return None

    def mousePressEvent(self, event) -> None:  # noqa: N802
        hit = self._hit(event.position())
        if hit is None:
            return
        kind, i = hit
        if i != self.index:
            self.index = i
            self.band_index = 0 if self.karyotype.chromosomes[i].bands else -1
            self.selected.emit(i)
        self._drag = None if kind == "select" else (kind, i)
        self.update()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        k = self.karyotype
        if self._drag is None or k is None:
            hit = self._hit(event.position())
            self.setCursor(Qt.CursorShape.SizeVerCursor if hit and hit[0] != "select" else Qt.CursorShape.ArrowCursor)
            return
        kind, i = self._drag
        c = k.chromosomes[i]
        r = self._geometry(i, c)
        f = max(0.0, min(1.0, (event.position().y() - r.top()) / max(r.height(), 1.0)))
        if kind == "centromere":
            c.centromere_index = round(max(5.0, min(95.0, f * 100.0)), 1)
        elif 0 <= self.band_index < len(c.bands):
            b = c.bands[self.band_index]
            if kind == "start":
                b.start = min(f, b.end)
            else:
                b.end = max(f, b.start)
        self.edited.emit()
        self.update()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        self._drag = None


class KaryotypeDialog(QDialog):
    """Create or edit one karyotype. Changes are written on Save."""

    def __init__(self, folder: Path, name: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Edit karyotype")
        self.setMinimumSize(1100, 700)
        self._folder = Path(folder).expanduser().resolve()
        existing = kt.load_karyotype(self._folder, name) if name else None
        self._old_name = existing.name if existing else ""
        self._k: kt.Karyotype = copy.deepcopy(existing) if existing else kt.new_karyotype("New karyotype", 8, 2)
        self._listed = [c.name for c in self._k.chromosomes]
        self._saved = ""
        self._busy = False
        self._seq_names: list[str] = []

        # ---- header
        self._name = QLineEdit(self._k.name, self)
        self._name.setFixedWidth(150)
        self._species = NoWheelComboBox(self)
        self._species.setEditable(True)
        self._species.addItems([""] + [e.get("name", "") for e in load_lists(self._folder).get("species", [])])
        self._species.setCurrentText(self._k.species)
        self._species.setFixedWidth(170)
        self._species.setToolTip("Species this karyotype is for; images of this species pick it automatically")
        self._basic = QSpinBox(self)
        self._basic.setRange(1, 60)
        self._basic.setValue(self._k.basic_number)
        self._ploidy = NoWheelComboBox(self)
        for n in PLOIDIES:
            self._ploidy.addItem(f"{n}x", n)
        self._ploidy.setCurrentIndex(max(0, self._ploidy.findData(self._k.ploidy)))
        self._formula = QLabel(self)
        f = self._formula.font()
        f.setBold(True)
        self._formula.setFont(f)
        self._hom_names = QLabel(self)
        self._unit = NoWheelComboBox(self)
        for text, value in UNITS:
            self._unit.addItem(text, value)
        self._unit.setCurrentIndex(max(0, self._unit.findData(self._k.unit)))
        head = QHBoxLayout()
        head.setSpacing(10)
        for text, w in (("Name", self._name), ("Species", self._species), ("Basic number x", self._basic), ("Ploidy", self._ploidy)):
            head.addWidget(QLabel(text, self))
            head.addWidget(w)
        head.addWidget(self._formula)
        head.addWidget(self._hom_names)
        head.addStretch(1)
        head.addWidget(QLabel("Length unit", self))
        head.addWidget(self._unit)
        head.addWidget(make_help_button(self, "Karyotype", HELP))

        # ---- table side
        add_btn = QPushButton("+ Add chromosome", self)
        imp_btn = QPushButton("Import measurements (KaryoScope TSV)", self)
        asm_btn = QPushButton("Take lengths from assembly", self)
        add_btn.clicked.connect(self._add_chromosome)
        imp_btn.clicked.connect(self._import_measurements)
        asm_btn.clicked.connect(self._take_from_assembly)
        tools = QHBoxLayout()
        tools.addWidget(add_btn)
        tools.addWidget(imp_btn)
        tools.addWidget(asm_btn)
        tools.addStretch(1)
        tools.addWidget(QLabel("Order", self))
        self._order_btns: dict[str, QPushButton] = {}
        for key, text in (("listed", "As listed"), ("length", "Longest first"), ("ci", "By centromere index")):
            b = QPushButton(text, self)
            b.setCheckable(True)
            b.clicked.connect(lambda _c=False, k=key: self._order(k))
            self._order_btns[key] = b
            tools.addWidget(b)

        self._table = _ReorderTable(0, 9, self)
        self._table.setHorizontalHeaderLabels(["", "Name", f"Length ({self._k.unit})", "Centromere index (%)", "Type", "Copies",
                                               "Assembly sequence", "Bands", ""])
        self._table.verticalHeader().setVisible(False)
        self._table.setShowGrid(False)
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(6, QHeaderView.ResizeMode.Stretch)
        self._table.rowMoved.connect(self._move_row)
        self._table.currentCellChanged.connect(lambda r, _c, _pr, _pc: self._select_row(r))
        note = QLabel("Type follows Levan et al. (1964) from the centromere index. Drag ⋮⋮ to change the order. "
                      "Copies can differ from the ploidy (aneuploidy, B chromosomes).", self)
        note.setWordWrap(True)
        note.setStyleSheet("color: palette(dark); font-size: 12px;")
        left = QVBoxLayout()
        left.addLayout(tools)
        left.addWidget(self._table, 1)
        left.addWidget(note)

        # ---- preview side
        self._preview = KaryoPreview(self)
        self._preview.karyotype = self._k
        self._preview.selected.connect(self._on_preview_selected)
        self._preview.edited.connect(self._on_preview_edited)
        ph = QHBoxLayout()
        ph.addWidget(QLabel("<b>Preview</b>", self), 1)
        hint = QLabel("drag a centromere or a band edge", self)
        hint.setStyleSheet("color: palette(dark); font-size: 12px;")
        ph.addWidget(hint)
        self._band_title = QLabel(self)
        self._band_list = NoWheelComboBox(self)
        self._band_list.activated.connect(self._on_band_chosen)
        self._kind = NoWheelComboBox(self)
        for kind in kt.BAND_KINDS:
            self._kind.addItem(kind)
        self._kind.addItem(CUSTOM_KIND)
        self._kind.activated.connect(self._on_kind)
        self._from = QLineEdit(self)
        self._to = QLineEdit(self)
        for e in (self._from, self._to):
            e.setFixedWidth(80)
            e.editingFinished.connect(self._on_band_range)
        add_band = QPushButton("+ Add band", self)
        del_band = QPushButton("Remove band", self)
        add_band.clicked.connect(self._add_band)
        del_band.clicked.connect(self._remove_band)
        form = QVBoxLayout()
        form.addWidget(self._band_title)
        for text, w in (("Band", self._band_list), ("Kind", self._kind)):
            r = QHBoxLayout()
            lab = QLabel(text, self)
            lab.setFixedWidth(70)
            r.addWidget(lab)
            r.addWidget(w, 1)
            form.addLayout(r)
        r = QHBoxLayout()
        lab = QLabel("From / to", self)
        lab.setFixedWidth(70)
        r.addWidget(lab)
        r.addWidget(self._from)
        r.addWidget(self._to)
        r.addStretch(1)
        form.addLayout(r)
        br = QHBoxLayout()
        br.addWidget(add_band)
        br.addWidget(del_band)
        br.addStretch(1)
        form.addLayout(br)
        right = QVBoxLayout()
        right.addLayout(ph)
        right.addWidget(self._preview, 1)
        right.addLayout(form)
        right_host = QWidget(self)
        right_host.setLayout(right)
        right_host.setFixedWidth(430)

        body = QHBoxLayout()
        body.addLayout(left, 1)
        body.addWidget(right_host)

        # ---- footer
        self._footer = QLabel(self)
        self._footer.setStyleSheet("color: palette(dark); font-size: 12px;")
        cancel = QPushButton("Cancel", self)
        save = QPushButton("Save", self)
        save.setDefault(True)
        save.setStyleSheet("QPushButton { background:#1F4E9A; color:#FFFFFF; border:1px solid #1F4E9A; border-radius:6px; padding:6px 16px; }")
        cancel.clicked.connect(self.reject)
        save.clicked.connect(self._save)
        foot = QHBoxLayout()
        foot.addWidget(self._footer, 1)
        foot.addWidget(cancel)
        foot.addWidget(save)

        root = QVBoxLayout(self)
        root.addLayout(head)
        root.addLayout(body, 1)
        root.addLayout(foot)

        self._name.textChanged.connect(self._on_name)
        self._species.editTextChanged.connect(lambda t: setattr(self._k, "species", t.strip()))
        self._basic.editingFinished.connect(self._on_basic)
        self._ploidy.activated.connect(self._on_ploidy)
        self._unit.activated.connect(self._on_unit)
        self._row = 0 if self._k.chromosomes else -1
        self._preview.index = self._row
        self._fill_table()
        self._sync_head()
        self._sync_bands()

    # ------------------------------------------------------------------ results

    def saved_name(self) -> str:
        return self._saved

    # ------------------------------------------------------------------ header

    def _sync_head(self) -> None:
        k = self._k
        self._formula.setText(k.formula)
        letters = ", ".join(k.homologue_letters[: max(1, k.ploidy)])
        self._hom_names.setText(f"Homologues named {letters}")
        self._footer.setText(f"Saved as {kt.karyotype_path(self._folder, k.name).name} in the folder; used by Signals, Idiogram and Figure.")
        header = self._table.horizontalHeaderItem(2)
        if header is not None:
            header.setText(f"Length ({k.unit})")

    def _on_name(self, text: str) -> None:
        self._k.name = text.strip()
        self._sync_head()

    def _on_basic(self) -> None:
        n = self._basic.value()
        k = self._k
        if n == k.basic_number and n == len(k.chromosomes):
            return
        dropped = k.chromosomes[n:]
        if any(c.length > 0 or c.bands for c in dropped):
            answer = QMessageBox.question(self, "Basic number", f"Remove {len(dropped)} chromosome(s) with measurements from the end of the list?")
            if answer != QMessageBox.StandardButton.Yes:
                self._basic.setValue(len(k.chromosomes))
                return
        k.set_basic_number(n)
        self._relist()
        self._row = min(self._row, len(k.chromosomes) - 1)
        self._fill_table()
        self._sync_head()
        self._sync_bands()

    def _on_ploidy(self, index: int) -> None:
        self._k.ploidy = int(self._ploidy.itemData(index))
        self._fill_table()
        self._sync_head()

    def _on_unit(self, index: int) -> None:
        self._k.unit = str(self._unit.itemData(index))
        self._sync_head()

    # ------------------------------------------------------------------ table

    def _relist(self) -> None:
        self._listed = [c.name for c in self._k.chromosomes]
        self._order_btns["listed"].setChecked(True)
        self._order_btns["length"].setChecked(False)
        self._order_btns["ci"].setChecked(False)

    def _fill_table(self) -> None:
        self._busy = True
        k = self._k
        self._table.setRowCount(len(k.chromosomes))
        for r, c in enumerate(k.chromosomes):
            handle = QTableWidgetItem("⋮⋮")
            handle.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsDragEnabled)
            handle.setForeground(QColor("#8A8D93"))
            self._table.setItem(r, 0, handle)
            self._table.setCellWidget(r, 1, self._edit(r, c.name, 52, self._on_name_edit))
            self._table.setCellWidget(r, 2, self._edit(r, f"{c.length:g}", 70, self._on_length_edit))
            self._table.setCellWidget(r, 3, self._edit(r, f"{c.centromere_index:g}", 60, self._on_ci_edit))
            type_item = QTableWidgetItem(c.type)
            type_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
            self._table.setItem(r, 4, type_item)
            self._table.setCellWidget(r, 5, self._edit(r, str(k.copies_of(c)), 40, self._on_copies_edit))
            combo = NoWheelComboBox(self)
            combo.setEditable(True)
            combo.addItems([""] + self._seq_names)
            combo.setCurrentText(c.assembly_seq)
            combo.lineEdit().editingFinished.connect(lambda r=r, cb=combo: self._on_seq_edit(r, cb.currentText()))
            self._table.setCellWidget(r, 6, combo)
            bands_item = QTableWidgetItem(kt.band_summary(c))
            bands_item.setFlags(Qt.ItemFlag.ItemIsEnabled)
            self._table.setItem(r, 7, bands_item)
            remove = QPushButton("×", self)
            remove.setFixedSize(22, 22)
            remove.setToolTip("Remove chromosome")
            remove.clicked.connect(lambda _c=False, r=r: self._remove_chromosome(r))
            self._table.setCellWidget(r, 8, remove)
            self._table.setRowHeight(r, 30)
        if 0 <= self._row < len(k.chromosomes):
            self._table.selectRow(self._row)
        self._preview.index = self._row
        self._preview.update()
        self._busy = False

    def _edit(self, row: int, text: str, width: int, slot) -> _Edit:
        e = _Edit(text, self)
        e.setFixedWidth(width)
        e.focused.connect(lambda r=row: self._select_row(r))
        e.editingFinished.connect(lambda r=row, w=e: slot(r, w.text()))
        return e

    def _refresh_row(self, r: int) -> None:
        c = self._k.chromosomes[r]
        self._table.item(r, 4).setText(c.type)
        self._table.item(r, 7).setText(kt.band_summary(c))
        self._preview.update()

    def _on_name_edit(self, r: int, text: str) -> None:
        if self._busy or r >= len(self._k.chromosomes):
            return
        text = text.strip()
        other = [c.name for i, c in enumerate(self._k.chromosomes) if i != r]
        if not text or text in other:
            QMessageBox.information(self, "Chromosome name", "Every chromosome needs its own name.")
            self._table.cellWidget(r, 1).setText(self._k.chromosomes[r].name)
            return
        old = self._k.chromosomes[r].name
        self._k.chromosomes[r].name = text
        self._listed = [text if n == old else n for n in self._listed]
        self._preview.update()

    def _num(self, text: str) -> float | None:
        try:
            return float(text.replace(",", "."))
        except ValueError:
            return None

    def _on_length_edit(self, r: int, text: str) -> None:
        if self._busy or r >= len(self._k.chromosomes):
            return
        v = self._num(text)
        if v is None or v < 0:
            self._table.cellWidget(r, 2).setText(f"{self._k.chromosomes[r].length:g}")
            return
        self._k.chromosomes[r].length = v
        self._preview.update()

    def _on_ci_edit(self, r: int, text: str) -> None:
        if self._busy or r >= len(self._k.chromosomes):
            return
        v = self._num(text)
        if v is None or not 0 < v <= 100:
            self._table.cellWidget(r, 3).setText(f"{self._k.chromosomes[r].centromere_index:g}")
            return
        self._k.chromosomes[r].centromere_index = min(v, 100 - v) if v > 50 else v
        self._table.cellWidget(r, 3).setText(f"{self._k.chromosomes[r].centromere_index:g}")
        self._refresh_row(r)

    def _on_copies_edit(self, r: int, text: str) -> None:
        if self._busy or r >= len(self._k.chromosomes):
            return
        try:
            n = int(text)
            if n < 0:
                raise ValueError
        except ValueError:
            self._table.cellWidget(r, 5).setText(str(self._k.copies_of(self._k.chromosomes[r])))
            return
        self._k.chromosomes[r].copies = None if n == self._k.ploidy else n
        self._sync_head()

    def _on_seq_edit(self, r: int, text: str) -> None:
        if not self._busy and r < len(self._k.chromosomes):
            self._k.chromosomes[r].assembly_seq = text.strip()

    def _select_row(self, r: int) -> None:
        if self._busy or r < 0 or r >= len(self._k.chromosomes) or r == self._row and self._preview.index == r:
            return
        self._row = r
        self._preview.index = r
        self._preview.band_index = 0 if self._k.chromosomes[r].bands else -1
        self._table.selectRow(r)
        self._preview.update()
        self._sync_bands()

    def _on_preview_selected(self, i: int) -> None:
        self._row = i
        self._table.selectRow(i)
        self._sync_bands()

    def _on_preview_edited(self) -> None:
        if 0 <= self._preview.index < len(self._k.chromosomes):
            r = self._preview.index
            c = self._k.chromosomes[r]
            self._table.cellWidget(r, 3).setText(f"{c.centromere_index:g}")
            self._refresh_row(r)
            self._sync_bands(keep_index=True)

    def _add_chromosome(self) -> None:
        k = self._k
        used = {c.name for c in k.chromosomes}
        n = len(k.chromosomes) + 1
        while str(n) in used:
            n += 1
        k.chromosomes.append(kt.KaryoChromosome(name=str(n)))
        k.basic_number = len(k.chromosomes)
        self._basic.blockSignals(True)
        self._basic.setValue(k.basic_number)
        self._basic.blockSignals(False)
        self._row = len(k.chromosomes) - 1
        self._relist()
        self._fill_table()
        self._sync_head()
        self._sync_bands()

    def _remove_chromosome(self, r: int) -> None:
        k = self._k
        if len(k.chromosomes) <= 1:
            return
        c = k.chromosomes[r]
        if c.length > 0 or c.bands:
            answer = QMessageBox.question(self, "Remove chromosome", f"Remove chromosome {c.name} with its measurements?")
            if answer != QMessageBox.StandardButton.Yes:
                return
        del k.chromosomes[r]
        k.basic_number = len(k.chromosomes)
        self._basic.blockSignals(True)
        self._basic.setValue(k.basic_number)
        self._basic.blockSignals(False)
        self._row = min(r, len(k.chromosomes) - 1)
        self._relist()
        self._fill_table()
        self._sync_head()
        self._sync_bands()

    def _move_row(self, src: int, dst: int) -> None:
        chroms = self._k.chromosomes
        chroms.insert(dst, chroms.pop(src))
        self._row = dst
        self._relist()
        self._fill_table()

    def _order(self, key: str) -> None:
        k = self._k
        if key == "listed":
            rank = {name: i for i, name in enumerate(self._listed)}
            k.chromosomes.sort(key=lambda c: rank.get(c.name, len(rank)))
        else:
            k.chromosomes[:] = k.ordered(key)
        for name, b in self._order_btns.items():
            b.setChecked(name == key)
        self._row = 0 if k.chromosomes else -1
        self._fill_table()
        self._sync_bands()

    # ------------------------------------------------------------------ tools

    def _import_measurements(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Import measurements", str(self._folder), "Tables (*.tsv *.csv *.txt);;All files (*)")
        if not path:
            return
        try:
            rows = kt.parse_measurement_table(Path(path).read_text(encoding="utf-8", errors="replace"))
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Import measurements", str(exc))
            return
        if not rows:
            QMessageBox.information(self, "Import measurements",
                                    "No chromosome rows found. The table needs a chromosome column and a length "
                                    "(or p and q) column.")
            return
        unit = self._k.unit if self._k.unit in ("µm", "%") else "µm"
        missing = kt.apply_measurements(self._k, rows, unit)
        self._unit.setCurrentIndex(max(0, self._unit.findData(unit)))
        self._fill_table()
        self._sync_head()
        text = f"{len(rows) - len(missing)} chromosome(s) updated."
        if missing:
            text += "\nNot in this karyotype: " + ", ".join(missing)
        QMessageBox.information(self, "Import measurements", text)

    def _take_from_assembly(self) -> None:
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            lengths = kt.folder_sequence_lengths(self._folder)
        finally:
            QApplication.restoreOverrideCursor()
        if not lengths:
            QMessageBox.information(self, "Assembly", "No FASTA file with sequences found in this folder.")
            return
        self._seq_names = sorted(lengths)
        linked = kt.guess_assembly_links(self._k, self._seq_names)
        missing = kt.apply_assembly_lengths(self._k, lengths)
        self._unit.setCurrentIndex(max(0, self._unit.findData("Mb")))
        self._k.unit = "Mb"
        self._fill_table()
        self._sync_head()
        text = f"{linked} new link(s) set; lengths are now in Mb."
        if missing:
            text += "\nNo linked sequence for chromosome: " + ", ".join(missing) + "\nChoose the sequence in the Assembly column."
        QMessageBox.information(self, "Assembly", text)

    # ------------------------------------------------------------------ bands

    def _current(self) -> kt.KaryoChromosome | None:
        return self._k.chromosomes[self._row] if 0 <= self._row < len(self._k.chromosomes) else None

    def _sync_bands(self, keep_index: bool = False) -> None:
        c = self._current()
        self._band_title.setText(f"<b>Band on selected chromosome ({c.name})</b>" if c else "<b>No chromosome selected</b>")
        if c is None:
            self._band_list.clear()
            return
        if not keep_index:
            self._preview.band_index = min(max(self._preview.band_index, 0), len(c.bands) - 1) if c.bands else -1
        self._band_list.blockSignals(True)
        self._band_list.clear()
        for i, b in enumerate(c.bands):
            self._band_list.addItem(f"{i + 1}: {b.kind}", i)
        if not c.bands:
            self._band_list.addItem("no bands", -1)
        self._band_list.setCurrentIndex(max(self._preview.band_index, 0))
        self._band_list.blockSignals(False)
        has = 0 <= self._preview.band_index < len(c.bands)
        for w in (self._kind, self._from, self._to):
            w.setEnabled(has)
        if has:
            b = c.bands[self._preview.band_index]
            idx = self._kind.findText(b.kind)
            if idx < 0:
                self._kind.insertItem(self._kind.count() - 1, b.kind)
                idx = self._kind.findText(b.kind)
            self._kind.setCurrentIndex(idx)
            f = c.centromere_fraction
            self._from.setText(kt.fraction_to_arm_text(b.start, f))
            self._to.setText(kt.fraction_to_arm_text(b.end, f))
        else:
            self._from.clear()
            self._to.clear()
        self._preview.update()

    def _on_band_chosen(self, index: int) -> None:
        self._preview.band_index = self._band_list.itemData(index) if self._band_list.itemData(index) is not None else -1
        self._sync_bands(keep_index=True)

    def _band(self) -> kt.Band | None:
        c = self._current()
        if c is None or not 0 <= self._preview.band_index < len(c.bands):
            return None
        return c.bands[self._preview.band_index]

    def _on_kind(self, index: int) -> None:
        b = self._band()
        if b is None:
            return
        text = self._kind.itemText(index)
        if text == CUSTOM_KIND:
            text, ok = QInputDialog.getText(self, "Band kind", "Name of the band kind:")
            if not ok or not text.strip():
                self._sync_bands(keep_index=True)
                return
            text = text.strip()
        b.kind = text
        self._refresh_row(self._row)
        self._sync_bands(keep_index=True)

    def _on_band_range(self) -> None:
        c, b = self._current(), self._band()
        if c is None or b is None:
            return
        f = c.centromere_fraction
        start = kt.arm_text_to_fraction(self._from.text(), f)
        end = kt.arm_text_to_fraction(self._to.text(), f)
        if start is None or end is None:
            QMessageBox.information(self, "Band", "Write the position as p or q and a number from 0 to 1, for example q 0.92.")
        else:
            b.start, b.end = min(start, end), max(start, end)
        self._refresh_row(self._row)
        self._sync_bands(keep_index=True)

    def _add_band(self) -> None:
        c = self._current()
        if c is None:
            return
        f = c.centromere_fraction
        c.bands.append(kt.Band(kt.BAND_KINDS[0], min(1.0, f + (1 - f) * 0.9), 1.0))
        self._preview.band_index = len(c.bands) - 1
        self._refresh_row(self._row)
        self._sync_bands(keep_index=True)

    def _remove_band(self) -> None:
        c, b = self._current(), self._band()
        if c is None or b is None:
            return
        c.bands.remove(b)
        self._preview.band_index = min(self._preview.band_index, len(c.bands) - 1)
        self._refresh_row(self._row)
        self._sync_bands(keep_index=True)

    # ------------------------------------------------------------------ save

    def _save(self) -> None:
        k = self._k
        k.name = self._name.text().strip()
        if not k.name:
            QMessageBox.information(self, "Karyotype", "Give the karyotype a name.")
            return
        names = [c.name for c in k.chromosomes]
        if len(set(names)) != len(names) or not all(names):
            QMessageBox.information(self, "Karyotype", "Every chromosome needs its own name.")
            return
        if k.name != self._old_name and kt.karyotype_path(self._folder, k.name).exists():
            answer = QMessageBox.question(self, "Karyotype", f"A karyotype named \"{k.name}\" exists. Replace it?")
            if answer != QMessageBox.StandardButton.Yes:
                return
        k.basic_number = len(k.chromosomes)
        try:
            kt.save_karyotype(self._folder, k, old_name=self._old_name)
        except OSError as exc:
            QMessageBox.warning(self, "Karyotype", f"Could not save:\n{exc}")
            return
        from repeatviewer.widgets.cyto_images_ui import notify_changed

        self._saved = k.name
        notify_changed(self._folder)
        self.accept()
