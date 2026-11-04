"""Edit lists: dyes, probes, species, cell stages, techniques, counterstains (macOS Settings style)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSize, Qt, QTimer
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QColorDialog,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.cytogenetics import lists as ls
from repeatviewer.views.cyto_karyotype_dialog import _ReorderTable
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cyto_images_ui import _fasta_files_in, notify_changed
from repeatviewer.widgets.help_button import make_help_button

HELP = (
    "<b>Lists</b> hold the values offered when you import, in the Library filters and in the image details.<br>"
    "<b>+</b> adds an entry, <b>−</b> removes the selected ones (only while no image uses them).<br>"
    "<b>Rename</b>: double-click a name. Every image that uses it is updated when you click Save.<br>"
    "<b>Merge…</b>: select two or more rows; the images move to the entry you keep.<br>"
    "<b>Order</b>: drag a row by its first column. The order is the order of the menus everywhere.<br>"
    "<b>Dyes</b>: excitation and emission maxima in nm, the words that identify the dye in file names "
    "(comma separated) and the colour it is shown in.<br>"
    "<b>Probes</b>: the sequence a probe comes from links its images to that FASTA file; the dye is the one it is "
    "usually detected with. Choose <i>Add new dye…</i> in the Dye menu to add a dye on the spot."
)

ADD_DYE = "Add new dye…"
CHOOSE_FILE = "Choose a file…"

# (header, field, kind) per list; Name comes first and "Used in" last
_COLUMNS = {
    "dyes": [("Excitation (nm)", "excitation", "number"), ("Emission (nm)", "emission", "number"),
             ("Words in file names", "words", "text"), ("Colour", "color", "color")],
    "probes": [("Sequence", "sequence", "sequence"), ("Dye", "dye", "dye")],
    "species": [],
    "stages": [("Signals are marked on", "marks", "marks")],
    "techniques": [],
    "counterstains": [],
}
_DESCRIPTIONS = {
    "dyes": "Fluorochromes and their channels. Import uses the words to sort TIFFs into channels.",
    "probes": "Probes and the sequence they come from; images with a probe also show under its sequence.",
    "species": "Species of your slides. Each species can have its own karyotype (Signals > Edit karyotype).",
    "stages": "Cell stages. The stage decides what the Signals tab marks: chromosomes, bivalents, nuclei or fibres.",
    "techniques": "Cytogenetic techniques.",
    "counterstains": "Counterstains of the chromosomes.",
}
ROW_H = 30


class _NavRow(QWidget):
    def __init__(self, title: str, count: int, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(8, 0, 8, 0)
        name = QLabel(title, self)
        num = QLabel(str(count), self)
        self._num = num
        lay.addWidget(name, 1)
        lay.addWidget(num)
        self.set_selected(False)

    def set_selected(self, on: bool) -> None:
        self._num.setStyleSheet("color: palette(highlighted-text);" if on else "color: gray;")


class CytoListsDialog(QDialog):
    """Edit the lists of a data folder. Changes are applied on Save."""

    def __init__(self, folder: Path, parent: QWidget | None = None, *, start: str = "dyes") -> None:
        super().__init__(parent)
        self.setWindowTitle("Edit lists")
        self.setMinimumSize(980, 620)
        self._folder = Path(folder).expanduser().resolve()
        self._editor = ls.ListsEditor(self._folder)
        ls.ensure_used_values(self._folder, self._editor.lists)
        self._key = start if start in ls.EDIT_ORDER else "dyes"
        self._loading = False

        self._nav = QListWidget(self)
        self._nav.setFixedWidth(220)
        self._nav.setSpacing(1)
        self._nav.currentRowChanged.connect(self._on_nav)

        self._title = QLabel(self)
        f = self._title.font()
        f.setBold(True)
        f.setPointSizeF(f.pointSizeF() + 2)
        self._title.setFont(f)
        self._desc = QLabel(self)
        self._desc.setWordWrap(True)
        self._desc.setStyleSheet("color: gray;")

        self._table = _ReorderTable(0, 0, self)
        self._table.setSelectionMode(_ReorderTable.SelectionMode.ExtendedSelection)
        self._table.verticalHeader().setVisible(False)
        self._table.verticalHeader().setDefaultSectionSize(ROW_H)
        self._table.setShowGrid(False)
        self._table.setAlternatingRowColors(True)
        self._table.setWordWrap(False)
        self._table.itemChanged.connect(self._on_item_changed)
        self._table.rowMoved.connect(self._on_row_moved)

        add = QPushButton("+", self)
        remove = QPushButton("−", self)
        for b, tip in ((add, "Add an entry"), (remove, "Remove the selected entries")):
            b.setFixedWidth(32)
            b.setToolTip(tip)
        merge = QPushButton("Merge…", self)
        merge.setToolTip("Merge the selected entries into one")
        add.clicked.connect(self._add)
        remove.clicked.connect(self._remove)
        merge.clicked.connect(self._merge)
        under = QHBoxLayout()
        under.setSpacing(4)
        under.addWidget(add)
        under.addWidget(remove)
        under.addSpacing(8)
        under.addWidget(merge)
        under.addStretch(1)

        head = QHBoxLayout()
        head.addWidget(self._title, 1)
        head.addWidget(make_help_button(self, "Edit lists", HELP))
        right = QVBoxLayout()
        right.setSpacing(8)
        right.addLayout(head)
        right.addWidget(self._desc)
        right.addWidget(self._table, 1)
        right.addLayout(under)

        body = QHBoxLayout()
        body.setSpacing(16)
        body.addWidget(self._nav)
        body.addLayout(right, 1)

        note = QLabel("Lists are kept per folder. Renaming here updates every image that uses the entry.", self)
        note.setStyleSheet("color: gray;")
        cancel = QPushButton("Cancel", self)
        save = QPushButton("Save", self)
        save.setDefault(True)
        cancel.clicked.connect(self.reject)
        save.clicked.connect(self._save)
        buttons = QHBoxLayout()
        buttons.addWidget(note, 1)
        buttons.addWidget(cancel)
        buttons.addWidget(save)

        root = QVBoxLayout(self)
        root.setContentsMargins(20, 20, 20, 16)
        root.setSpacing(16)
        root.addLayout(body, 1)
        root.addLayout(buttons)

        self._fill_nav()
        self._nav.setCurrentRow(ls.EDIT_ORDER.index(self._key))
        self._fill_table()  # the row may already be current, so currentRowChanged does not fire

    # ------------------------------------------------------------------ navigation

    def _fill_nav(self) -> None:
        row = self._nav.currentRow()
        self._nav.blockSignals(True)
        self._nav.clear()
        for key in ls.EDIT_ORDER:
            item = QListWidgetItem()
            item.setSizeHint(QSize(200, 28))
            self._nav.addItem(item)
            self._nav.setItemWidget(item, _NavRow(ls.LIST_TITLES[key], len(self._editor.names(key)), self._nav))
        self._nav.setCurrentRow(row if row >= 0 else ls.EDIT_ORDER.index(self._key))
        self._nav.blockSignals(False)
        self._mark_nav()

    def _mark_nav(self) -> None:
        for i in range(self._nav.count()):
            w = self._nav.itemWidget(self._nav.item(i))
            if isinstance(w, _NavRow):
                w.set_selected(i == self._nav.currentRow())

    def _on_nav(self, row: int) -> None:
        self._mark_nav()
        if 0 <= row < len(ls.EDIT_ORDER):
            self._key = ls.EDIT_ORDER[row]
            self._fill_table()

    # ------------------------------------------------------------------ table

    def _fill_table(self) -> None:
        self._loading = True
        key = self._key
        cols = _COLUMNS[key]
        self._title.setText(ls.LIST_TITLES[key])
        self._desc.setText(_DESCRIPTIONS.get(key, ""))
        self._table.clear()
        self._table.setColumnCount(len(cols) + 2)
        self._table.setHorizontalHeaderLabels(["Name", *[c[0] for c in cols], "Used in"])
        entries = self._editor.lists.get(key, [])
        self._table.setRowCount(len(entries))
        for r, e in enumerate(entries):
            name = str(e["name"])
            item = QTableWidgetItem(name)
            item.setData(Qt.ItemDataRole.UserRole, name)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsDragEnabled)
            self._table.setItem(r, 0, item)
            for c, (_h, field, kind) in enumerate(cols, start=1):
                self._set_cell(r, c, name, e, field, kind)
            n = self._editor.used_by(key, name)
            used = QTableWidgetItem(f"{n} image{'s' if n != 1 else ''}" if n else "")
            used.setFlags(used.flags() & ~Qt.ItemFlag.ItemIsEditable)
            used.setForeground(QColor("gray"))
            self._table.setItem(r, len(cols) + 1, used)
        header = self._table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        widths = {"sequence": QHeaderView.ResizeMode.Stretch, "text": QHeaderView.ResizeMode.Stretch}
        for c, (_h, _f, kind) in enumerate(cols, start=1):
            if kind in widths:
                header.setSectionResizeMode(c, widths[kind])
        if not any(kind in widths for _h, _f, kind in cols):
            header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self._table.resizeColumnsToContents()
        if not any(kind in widths for _h, _f, kind in cols):
            pass
        elif self._table.columnWidth(0) < 170:
            header.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
            self._table.setColumnWidth(0, 170)
        self._loading = False

    def _put(self, row: int, col: int, widget: QWidget) -> None:
        holder = QWidget(self)
        lay = QHBoxLayout(holder)
        lay.setContentsMargins(4, 2, 4, 2)
        lay.addWidget(widget, 1)
        self._table.setCellWidget(row, col, holder)

    def _set_cell(self, row: int, col: int, name: str, entry: dict, field: str, kind: str) -> None:
        value = str(entry.get(field, ""))
        if kind == "sequence":
            box = NoWheelComboBox(self)
            box.setEditable(True)
            box.addItem("")
            box.addItems(_fasta_files_in(self._folder))
            if value and box.findText(value) < 0:
                box.addItem(value)
            box.insertSeparator(box.count())
            box.addItem(CHOOSE_FILE)
            box.setCurrentText(value)
            box.setToolTip(value or "FASTA file this probe comes from")
            box.setMinimumContentsLength(24)
            box.activated.connect(lambda _i, n=name, b=box: self._on_sequence(n, b))
            box.lineEdit().editingFinished.connect(lambda n=name, b=box: self._set(n, "sequence", b.currentText()))
            self._put(row, col, box)
        elif kind == "dye":
            box = NoWheelComboBox(self)
            box.addItem("")
            box.addItems(self._editor.names("dyes"))
            box.insertSeparator(box.count())
            box.addItem(ADD_DYE)
            box.setCurrentText(value)
            box.activated.connect(lambda _i, n=name, b=box: self._on_dye(n, b))
            self._put(row, col, box)
        elif kind == "marks":
            box = NoWheelComboBox(self)
            box.addItems(list(ls.STAGE_MARKS))
            box.setCurrentText(value or "chromosomes")
            box.activated.connect(lambda _i, n=name, b=box: self._set(n, "marks", b.currentText()))
            self._put(row, col, box)
        elif kind == "color":
            holder = QWidget(self)
            lay = QHBoxLayout(holder)
            lay.setContentsMargins(6, 3, 6, 3)
            btn = QPushButton(holder)
            btn.setFixedSize(44, 20)
            btn.setToolTip(value or "Choose a colour")
            btn.setStyleSheet(f"QPushButton {{ background: {value or '#ffffff'}; border: 1px solid #9a9a9a; border-radius: 4px; }}")
            btn.clicked.connect(lambda _c=False, n=name, b=btn: self._pick_color(n, b))
            lay.addWidget(btn)
            lay.addStretch(1)
            self._table.setCellWidget(row, col, holder)
        else:
            item = QTableWidgetItem(value)
            item.setData(Qt.ItemDataRole.UserRole + 1, field)
            if kind == "number":
                item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            self._table.setItem(row, col, item)

    def _set(self, name: str, field: str, value: str) -> None:
        if not self._loading:
            self._editor.set_value(self._key, name, field, str(value).strip())

    def _on_sequence(self, name: str, box: NoWheelComboBox) -> None:
        if box.currentText() != CHOOSE_FILE:
            self._set(name, "sequence", box.currentText())
            return
        path, _ = QFileDialog.getOpenFileName(self, "Sequence of the probe", str(self._folder),
                                              "FASTA (*.fa *.fasta *.fna *.fas);;All files (*)")
        if not path:
            box.setCurrentText(str((ls.entry(self._editor.lists, "probes", name) or {}).get("sequence", "")))
            return
        p = Path(path)
        value = p.name if p.parent.resolve() == self._folder else str(p)
        self._set(name, "sequence", value)
        QTimer.singleShot(0, self._fill_table)

    def _on_dye(self, name: str, box: NoWheelComboBox) -> None:
        if box.currentText() != ADD_DYE:
            self._set(name, "dye", box.currentText())
            return
        dye, ok = QInputDialog.getText(self, "Add new dye", "Name of the dye:")
        if not ok or not dye.strip():
            box.setCurrentText(str((ls.entry(self._editor.lists, "probes", name) or {}).get("dye", "")))
            return
        colour = QColorDialog.getColor(QColor("#33ff33"), self, f"Colour of {dye.strip()}")
        try:
            self._editor.add("dyes", dye.strip(), color=colour.name() if colour.isValid() else "#ffffff",
                             words=dye.strip().lower())
        except ValueError as exc:
            QMessageBox.information(self, "Add new dye", str(exc))
        self._set(name, "dye", dye.strip())
        self._fill_nav()
        QTimer.singleShot(0, self._fill_table)

    def _pick_color(self, name: str, button: QPushButton) -> None:
        e = ls.entry(self._editor.lists, self._key, name) or {}
        chosen = QColorDialog.getColor(QColor(e.get("color", "#ffffff")), self, "Colour")
        if chosen.isValid():
            self._set(name, "color", chosen.name())
            button.setToolTip(chosen.name())
            button.setStyleSheet(f"QPushButton {{ background: {chosen.name()}; border: 1px solid #9a9a9a; border-radius: 4px; }}")

    def _on_item_changed(self, item: QTableWidgetItem) -> None:
        if self._loading:
            return
        row, col = item.row(), item.column()
        old = self._table.item(row, 0).data(Qt.ItemDataRole.UserRole)
        if col == 0:
            try:
                self._editor.rename(self._key, old, item.text())
            except ValueError as exc:
                QMessageBox.warning(self, "Edit lists", str(exc))
                self._loading = True
                item.setText(old)
                self._loading = False
                return
            QTimer.singleShot(0, self._fill_table)  # other cells refer to the entry by its name
        else:
            field = item.data(Qt.ItemDataRole.UserRole + 1)
            if field:
                self._editor.set_value(self._key, old, field, item.text().strip())

    def _on_row_moved(self, src: int, dst: int) -> None:
        self._editor.move(self._key, src, dst)
        QTimer.singleShot(0, self._fill_table)

    def _selected_names(self) -> list[str]:
        rows = sorted({i.row() for i in self._table.selectedIndexes()})
        return [self._table.item(r, 0).data(Qt.ItemDataRole.UserRole) for r in rows]

    # ------------------------------------------------------------------ actions

    def _add(self) -> None:
        name, ok = QInputDialog.getText(self, "Add", f"New entry in {ls.LIST_TITLES[self._key]}:")
        if not ok or not name.strip():
            return
        try:
            self._editor.add(self._key, name)
        except ValueError as exc:
            QMessageBox.warning(self, "Edit lists", str(exc))
            return
        self._fill_nav()
        self._fill_table()
        self._table.selectRow(self._table.rowCount() - 1)

    def _merge(self) -> None:
        names = self._selected_names()
        if len(names) < 2:
            QMessageBox.information(self, "Merge", "Select two or more rows to merge.")
            return
        keep, ok = QInputDialog.getItem(self, "Merge", "Keep this entry (the others are merged into it):", names, 0, False)
        if not ok:
            return
        self._editor.merge(self._key, names, keep)
        self._fill_nav()
        self._fill_table()

    def _remove(self) -> None:
        names = self._selected_names()
        if not names:
            QMessageBox.information(self, "Remove", "Select the rows to remove.")
            return
        blocked = [n for n in names if not self._editor.can_remove(self._key, n)]
        for n in names:
            if n not in blocked:
                self._editor.remove(self._key, n)
        if blocked:
            QMessageBox.information(
                self, "Remove",
                "Not removed because images still use them: " + ", ".join(blocked) + ".\nMerge them into another entry instead.",
            )
        self._fill_nav()
        self._fill_table()

    def _save(self) -> None:
        try:
            self._editor.apply()
        except OSError as exc:
            QMessageBox.critical(self, "Edit lists", f"Could not save the lists:\n{exc}")
            return
        notify_changed(self._folder)
        self.accept()
