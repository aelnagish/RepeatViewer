"""Edit lists: probes, species, cell stages, techniques, counterstains, dyes."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QColorDialog,
    QDialog,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.cytogenetics import lists as ls
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cyto_images_ui import _fasta_files_in, notify_changed
from repeatviewer.widgets.help_button import make_help_button

HELP = (
    "<b>Lists</b> hold the values offered in the import, the Library filters and the image details.<br>"
    "<b>Rename</b> a value by editing its name: every image that uses it is updated when you click Save.<br>"
    "<b>Merge selected</b> combines entries (select two or more rows); images are moved to the entry you keep.<br>"
    "<b>Remove</b> only works while no image uses the value; merge it into another entry otherwise.<br>"
    "<b>Probes</b> can name a sequence: images with that probe are then linked to the FASTA file, "
    "so they also show under the sequence.<br>"
    "<b>Dyes</b> list the words in file names that identify a channel (comma separated) and its default colour."
)

# columns per list: (header, field, kind); the name column is always first, "Used in" last
_COLUMNS = {
    "probes": [("Sequence", "sequence", "sequence"), ("Default dye", "dye", "dye")],
    "species": [],
    "stages": [("Signals are marked on", "marks", "marks")],
    "techniques": [],
    "counterstains": [],
    "dyes": [("Words in file names", "words", "text"), ("Colour", "color", "color")],
}


_FOOTS = {
    "probes": "A probe can name the sequence it comes from; images with that probe then also show under the sequence.",
    "species": "Each species points to its karyotype (Signals tab > Edit karyotype), so Signals and Idiogram pick the right one.",
    "stages": "The stage decides what the Signals tab marks: chromosomes, nuclei or fibres.",
    "dyes": "Import uses these words to sort TIFFs into channels. Any number of channels per image.",
}


class CytoListsDialog(QDialog):
    """Edit the lists of a data folder. Changes are applied on Save."""

    def __init__(self, folder: Path, parent: QWidget | None = None, *, start: str = "probes") -> None:
        super().__init__(parent)
        self.setWindowTitle("Edit lists")
        self.setMinimumSize(820, 520)
        self._folder = Path(folder).expanduser().resolve()
        self._editor = ls.ListsEditor(self._folder)
        ls.ensure_used_values(self._folder, self._editor.lists)
        self._key = start if start in ls.LIST_KEYS else "probes"
        self._loading = False

        self._nav = QListWidget(self)
        self._nav.setFixedWidth(190)
        self._nav.currentRowChanged.connect(self._on_nav)
        self._table = QTableWidget(0, 0, self)
        self._table.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QTableWidget.SelectionMode.ExtendedSelection)
        self._table.verticalHeader().setVisible(False)
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.itemChanged.connect(self._on_item_changed)

        self._title = QLabel(self)
        font = self._title.font()
        font.setBold(True)
        self._title.setFont(font)
        add = QPushButton("Add", self)
        merge = QPushButton("Merge selected", self)
        remove = QPushButton("Remove", self)
        add.clicked.connect(self._add)
        merge.clicked.connect(self._merge)
        remove.clicked.connect(self._remove)
        bar = QHBoxLayout()
        bar.addWidget(self._title, 1)
        for b in (add, merge, remove):
            bar.addWidget(b)
        bar.addWidget(make_help_button(self, "Edit lists", HELP))
        self._foot = QLabel(self)
        self._foot.setObjectName("cytoMuted")
        self._foot.setWordWrap(True)
        right = QVBoxLayout()
        right.addLayout(bar)
        right.addWidget(self._table, 1)
        right.addWidget(self._foot)
        body = QHBoxLayout()
        body.addWidget(self._nav)
        body.addLayout(right, 1)

        cancel = QPushButton("Cancel", self)
        save = QPushButton("Save", self)
        save.setDefault(True)
        cancel.clicked.connect(self.reject)
        save.clicked.connect(self._save)
        buttons = QHBoxLayout()
        note = QLabel("Lists are kept per folder. Renaming here updates every image that uses the entry.", self)
        note.setObjectName("cytoMuted")
        note.setWordWrap(True)
        buttons.addWidget(note, 1)
        buttons.addWidget(cancel)
        buttons.addWidget(save)
        root = QVBoxLayout(self)
        root.addLayout(body, 1)
        root.addLayout(buttons)

        self._fill_nav()
        self._nav.setCurrentRow(ls.LIST_KEYS.index(self._key))
        self._fill_table()  # the row may already be current, so currentRowChanged does not fire

    # ------------------------------------------------------------------ navigation

    def _fill_nav(self) -> None:
        row = self._nav.currentRow()
        self._nav.blockSignals(True)
        self._nav.clear()
        for key in ls.LIST_KEYS:
            self._nav.addItem(QListWidgetItem(f"{ls.LIST_TITLES[key]}   ({len(self._editor.names(key))})"))
        self._nav.setCurrentRow(row if row >= 0 else ls.LIST_KEYS.index(self._key))
        self._nav.blockSignals(False)

    def _on_nav(self, row: int) -> None:
        if 0 <= row < len(ls.LIST_KEYS):
            self._key = ls.LIST_KEYS[row]
            self._fill_table()

    # ------------------------------------------------------------------ table

    def _fill_table(self) -> None:
        self._loading = True
        key = self._key
        cols = _COLUMNS[key]
        self._title.setText(ls.LIST_TITLES[key])
        self._foot.setText(_FOOTS.get(key, ""))
        self._table.clear()
        self._table.setColumnCount(len(cols) + 2)
        self._table.setHorizontalHeaderLabels(["Name", *[c[0] for c in cols], "Used in"])
        entries = self._editor.lists.get(key, [])
        self._table.setRowCount(len(entries))
        for r, e in enumerate(entries):
            name = str(e["name"])
            item = QTableWidgetItem(name)
            item.setData(Qt.ItemDataRole.UserRole, name)
            self._table.setItem(r, 0, item)
            for c, (_h, field, kind) in enumerate(cols, start=1):
                self._set_cell(r, c, name, e, field, kind)
            n = self._editor.used_by(key, name)
            used = QTableWidgetItem(f"{n} image{'s' if n != 1 else ''}" if n else "")
            used.setFlags(used.flags() & ~Qt.ItemFlag.ItemIsEditable)
            self._table.setItem(r, len(cols) + 1, used)
        self._table.setColumnWidth(0, 200)
        self._loading = False

    def _set_cell(self, row: int, col: int, name: str, entry: dict, field: str, kind: str) -> None:
        value = str(entry.get(field, ""))
        if kind in ("sequence", "dye", "marks"):
            choices = {
                "sequence": _fasta_files_in(self._folder),
                "dye": ls.names(self._editor.lists, "dyes"),
                "marks": list(ls.STAGE_MARKS),
            }[kind]
            box = NoWheelComboBox(self)
            box.setEditable(kind == "sequence")
            if kind != "marks":
                box.addItem("")
            box.addItems(choices)
            if value and box.findText(value) < 0:
                box.addItem(value)
            box.setCurrentText(value or ("chromosomes" if kind == "marks" else ""))
            box.currentTextChanged.connect(lambda t, n=name, f=field: self._editor.set_value(self._key, n, f, t.strip()))
            self._table.setCellWidget(row, col, box)
        elif kind == "color":
            btn = QPushButton(value or "choose", self)
            btn.setStyleSheet(f"background:{value}; color:{'#000' if QColor(value).lightness() > 128 else '#fff'};" if value else "")
            btn.clicked.connect(lambda _c=False, n=name, b=btn: self._pick_color(n, b))
            self._table.setCellWidget(row, col, btn)
        else:
            item = QTableWidgetItem(value)
            item.setData(Qt.ItemDataRole.UserRole + 1, field)
            self._table.setItem(row, col, item)

    def _pick_color(self, name: str, button: QPushButton) -> None:
        e = ls.entry(self._editor.lists, self._key, name) or {}
        chosen = QColorDialog.getColor(QColor(e.get("color", "#ffffff")), self, "Default colour")
        if chosen.isValid():
            hexcolor = chosen.name()
            self._editor.set_value(self._key, name, "color", hexcolor)
            button.setText(hexcolor)
            button.setStyleSheet(f"background:{hexcolor}; color:{'#000' if chosen.lightness() > 128 else '#fff'};")

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
            item.setData(Qt.ItemDataRole.UserRole, item.text().strip())
            # the other cells of this row refer to the entry by its name: rebuild after this signal ends
            QTimer.singleShot(0, self._fill_table)
        else:
            field = item.data(Qt.ItemDataRole.UserRole + 1)
            if field:
                self._editor.set_value(self._key, old, field, item.text().strip())

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
