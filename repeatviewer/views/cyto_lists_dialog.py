"""Edit lists: dyes, probes, species, cell stages, techniques, counterstains (macOS Settings style).

The table gives the overview; the form under it edits the selected entry. Text cells can also be edited in place.
Not modal: a click on a FASTA file in the Files panel fills the Sequence field of the selected probe.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSize, Qt, QTimer
from PySide6.QtGui import QColor, QIcon, QPixmap
from PySide6.QtWidgets import (
    QColorDialog,
    QDialog,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.cytogenetics import lists as ls
from repeatviewer.views.cyto_karyotype_dialog import _ReorderTable
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cyto_images_ui import notify_changed
from repeatviewer.widgets.cyto_panel_targets import PanelPathFields
from repeatviewer.widgets.cytogenetics_style import form_layout, rounded_box
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.path_browse import PathBrowseRow

HELP = (
    "<b>Lists</b> hold the values offered when you import, in the Library filters and in the image details.<br>"
    "<b>Edit</b>: select a row and change it in the form under the table. Names, numbers and words can also be "
    "edited in the table: double-click a cell, or select it and press Enter.<br>"
    "<b>Rename</b>: every image that uses the entry is updated when you click Save.<br>"
    "<b>+</b> adds an entry, <b>−</b> removes the selected ones (only while no image uses them).<br>"
    "<b>Merge…</b>: select two or more rows; the images move to the entry you keep.<br>"
    "<b>Order</b>: drag a row by its name. The order is the order of the menus everywhere.<br>"
    "<b>Columns</b>: drag a column edge to make it wider or narrower (never narrower than its title).<br>"
    "<b>Dyes</b>: excitation and emission maxima in nm, the words that identify the dye in file names "
    "(comma separated) and the colour it is shown in.<br>"
    "<b>Probes</b>: the sequence a probe comes from links its images to that FASTA file. Click the file in the "
    "Files panel, or Browse. The dye is the one it is usually detected with; <i>Add a dye…</i> adds one on the spot."
)

ADD_DYE = "Add a dye…"
NO_DYE = "None"

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
    "dyes": "Fluorochromes and their channels. Import uses the words to sort pictures into channels.",
    "probes": "Probes and the sequence they come from; images with a probe also show under its sequence.",
    "species": "Species of your slides. Each species can have its own karyotype.",
    "stages": "Cell stages. The stage decides what the Signals tab marks: chromosomes, bivalents, nuclei or fibres.",
    "techniques": "Cytogenetic techniques.",
    "counterstains": "Counterstains of the chromosomes.",
}
_INLINE = ("text", "number")  # kinds edited in the table cell itself
ROW_H = 26


def _swatch(colour: str) -> QIcon:
    pix = QPixmap(28, 14)
    pix.fill(QColor(colour or "#ffffff"))
    return QIcon(pix)


class _NavRow(QWidget):
    def __init__(self, title: str, count: int, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(8, 0, 8, 0)
        self._num = QLabel(str(count), self)
        lay.addWidget(QLabel(title, self), 1)
        lay.addWidget(self._num)
        self.set_selected(False)

    def set_selected(self, on: bool) -> None:
        self._num.setStyleSheet("color: palette(highlighted-text);" if on else "color: gray;")


class CytoListsDialog(PanelPathFields, QDialog):
    """Edit the lists of a data folder. Changes are applied on Save."""

    def __init__(self, folder: Path, parent: QWidget | None = None, *, start: str = "dyes") -> None:
        super().__init__(parent)
        self.setWindowTitle("Edit lists")
        self.setWindowModality(Qt.WindowModality.NonModal)
        self.setMinimumSize(760, 500)
        self._folder = Path(folder).expanduser().resolve()
        self._editor = ls.ListsEditor(self._folder)
        ls.ensure_used_values(self._folder, self._editor.lists)
        self._key = start if start in ls.EDIT_ORDER else "dyes"
        self._loading = False

        self._nav = QListWidget(self)
        self._nav.setSpacing(1)
        self._nav.currentRowChanged.connect(self._on_nav)
        nav_box = rounded_box(self._nav, self, margin=4)
        nav_box.setFixedWidth(200)

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
        self._table.setEditTriggers(_ReorderTable.EditTrigger.DoubleClicked | _ReorderTable.EditTrigger.EditKeyPressed
                                    | _ReorderTable.EditTrigger.SelectedClicked)
        self._table.verticalHeader().setVisible(False)
        self._table.verticalHeader().setDefaultSectionSize(ROW_H)
        self._table.setShowGrid(False)
        self._table.setAlternatingRowColors(True)
        self._table.setWordWrap(False)
        self._table.itemChanged.connect(self._on_item_changed)
        self._table.rowMoved.connect(self._on_row_moved)
        self._table.itemSelectionChanged.connect(self._fill_form)
        header = self._table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        table_box = rounded_box(self._table, self)

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
        hint = QLabel("Double-click a cell to edit it, or change the selected entry below.", self)
        hint.setStyleSheet("color: gray;")
        under.addWidget(hint)

        # ---- form for the selected entry
        self._form_host = QWidget(self)
        self._form = form_layout(self._form_host)
        self._f_name = QLineEdit(self)
        self._f_name.editingFinished.connect(self._on_form_name)
        self._f_ex = QLineEdit(self)
        self._f_em = QLineEdit(self)
        self._f_words = QLineEdit(self)
        self._f_words.setPlaceholderText("dapi, 405, blue")
        for w, field in ((self._f_ex, "excitation"), (self._f_em, "emission"), (self._f_words, "words")):
            w.editingFinished.connect(lambda w=w, fld=field: self._set_selected(fld, w.text()))
        for w in (self._f_ex, self._f_em):
            w.setPlaceholderText("nm")
        self._f_colour = QPushButton(self)
        self._f_colour.setText("Choose…")
        self._f_colour.setIconSize(QSize(28, 14))
        self._f_colour.clicked.connect(self._pick_colour)
        self._f_seq = PathBrowseRow(self, placeholder="Click a FASTA file in the Files panel, or Browse",
                                    browse_title="Sequence of the probe", start_dir=str(self._folder),
                                    file_filter="FASTA (*.fa *.fasta *.fna *.fas);;All files (*)", elide_middle=True)
        self._f_seq.panel_suffixes = (".fa", ".fasta", ".fna", ".fas")
        self._f_seq.path_changed.connect(self._on_form_sequence)
        self._f_dye = NoWheelComboBox(self)
        self._f_dye.setSizeAdjustPolicy(NoWheelComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self._f_dye.activated.connect(self._on_form_dye)
        self._f_marks = NoWheelComboBox(self)
        self._f_marks.addItems(list(ls.STAGE_MARKS))
        self._f_marks.activated.connect(lambda _i: self._set_selected("marks", self._f_marks.currentText()))
        self._rows = {}
        for label, w in (("Name", self._f_name), ("Excitation", self._f_ex), ("Emission", self._f_em),
                         ("Words in file names", self._f_words), ("Colour", self._f_colour),
                         ("Sequence", self._f_seq), ("Dye", self._f_dye), ("Signals are marked on", self._f_marks)):
            lab = QLabel(label, self._form_host)
            self._form.addRow(lab, w)
            self._rows[w] = lab
        for w in (self._f_name, self._f_ex, self._f_em, self._f_words, self._f_colour, self._f_seq, self._f_dye,
                  self._f_marks):
            w.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

        head = QHBoxLayout()
        head.addWidget(self._title, 1)
        head.addWidget(make_help_button(self, "Edit lists", HELP))
        right = QVBoxLayout()
        right.setSpacing(6)
        right.addLayout(head)
        right.addWidget(self._desc)
        right.addWidget(table_box, 1)
        right.addLayout(under)
        right.addSpacing(4)
        right.addWidget(self._form_host)

        body = QHBoxLayout()
        body.setSpacing(14)
        body.addWidget(nav_box)
        body.addLayout(right, 1)

        cancel = QPushButton("Cancel", self)
        save = QPushButton("Save", self)
        save.setDefault(True)
        cancel.clicked.connect(self.reject)
        save.clicked.connect(self._save)
        buttons = QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(cancel)
        buttons.addWidget(save)

        root = QVBoxLayout(self)
        root.setContentsMargins(16, 16, 16, 14)
        root.setSpacing(12)
        root.addLayout(body, 1)
        root.addLayout(buttons)

        self.watch_path_fields(self._f_seq)
        self._fill_nav()
        self._nav.setCurrentRow(ls.EDIT_ORDER.index(self._key))
        self._fill_table()
        self.resize(900, 600)

    # ------------------------------------------------------------------ navigation

    def _fill_nav(self) -> None:
        row = self._nav.currentRow()
        self._nav.blockSignals(True)
        self._nav.clear()
        for key in ls.EDIT_ORDER:
            item = QListWidgetItem()
            item.setSizeHint(QSize(180, 28))
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

    def _entry(self, name: str) -> dict:
        return ls.entry(self._editor.lists, self._key, name) or {}

    def _fill_table(self, select: str = "") -> None:
        self._loading = True
        key = self._key
        cols = _COLUMNS[key]
        self._title.setText(ls.LIST_TITLES[key])
        self._desc.setText(_DESCRIPTIONS.get(key, ""))
        self._table.clear()
        self._table.setColumnCount(len(cols) + 2)
        heads = ["Name", *[c[0] for c in cols], "Used in"]
        self._table.setHorizontalHeaderLabels(heads)
        entries = self._editor.lists.get(key, [])
        self._table.setRowCount(len(entries))
        for r, e in enumerate(entries):
            name = str(e["name"])
            item = QTableWidgetItem(name)
            item.setData(Qt.ItemDataRole.UserRole, name)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsDragEnabled | Qt.ItemFlag.ItemIsEditable)
            self._table.setItem(r, 0, item)
            for c, (_h, field, kind) in enumerate(cols, start=1):
                value = str(e.get(field, ""))
                cell = QTableWidgetItem()
                cell.setData(Qt.ItemDataRole.UserRole + 1, field)
                if kind == "color":
                    cell.setIcon(_swatch(value))
                    cell.setToolTip(value)
                elif kind == "sequence":
                    cell.setText(value)
                    cell.setToolTip(value)
                elif kind == "marks":
                    cell.setText(value or "chromosomes")
                else:
                    cell.setText(value)
                if kind == "number":
                    cell.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                flags = cell.flags()
                cell.setFlags(flags | Qt.ItemFlag.ItemIsEditable if kind in _INLINE else flags & ~Qt.ItemFlag.ItemIsEditable)
                self._table.setItem(r, c, cell)
            n = self._editor.used_by(key, name)
            used = QTableWidgetItem(f"{n} image{'s' if n != 1 else ''}" if n else "")
            used.setFlags(used.flags() & ~Qt.ItemFlag.ItemIsEditable)
            used.setForeground(QColor("gray"))
            self._table.setItem(r, len(cols) + 1, used)
        self._table.setIconSize(QSize(28, 14))
        header = self._table.horizontalHeader()
        for c in range(len(heads)):
            kind = cols[c - 1][2] if 0 < c <= len(cols) else "text"
            fit = kind in ("number", "color")
            header.setSectionResizeMode(c, QHeaderView.ResizeMode.ResizeToContents if fit
                                        else QHeaderView.ResizeMode.Stretch)
        self._loading = False
        row = next((r for r, e in enumerate(entries) if str(e["name"]) == select), 0 if entries else -1)
        if row >= 0:
            self._table.selectRow(row)
        self._fill_form()

    def _selected_names(self) -> list[str]:
        rows = sorted({i.row() for i in self._table.selectedIndexes()})
        return [self._table.item(r, 0).data(Qt.ItemDataRole.UserRole) for r in rows if self._table.item(r, 0)]

    def _current_name(self) -> str:
        names = self._selected_names()
        return names[0] if len(names) == 1 else ""

    def _on_item_changed(self, item: QTableWidgetItem) -> None:
        if self._loading:
            return
        row, col = item.row(), item.column()
        old = self._table.item(row, 0).data(Qt.ItemDataRole.UserRole)
        if col == 0:
            self._rename(old, item.text())
        else:
            field = item.data(Qt.ItemDataRole.UserRole + 1)
            if field:
                self._editor.set_value(self._key, old, field, item.text().strip())
                self._fill_form()

    def _rename(self, old: str, new: str) -> None:
        new = new.strip()
        if new == old:
            return
        try:
            self._editor.rename(self._key, old, new)
        except ValueError as exc:
            QMessageBox.warning(self, "Edit lists", str(exc))
            new = old
        QTimer.singleShot(0, lambda n=new: self._fill_table(select=n))

    def _on_row_moved(self, src: int, dst: int) -> None:
        name = self._table.item(src, 0).data(Qt.ItemDataRole.UserRole)
        self._editor.move(self._key, src, dst)
        QTimer.singleShot(0, lambda: self._fill_table(select=name))

    # ------------------------------------------------------------------ form

    def _fill_form(self) -> None:
        name = self._current_name()
        e = self._entry(name) if name else {}
        kinds = {kind for _h, _f, kind in _COLUMNS[self._key]}
        shown = {
            self._f_name: True,
            self._f_ex: "number" in kinds, self._f_em: "number" in kinds,
            self._f_words: "text" in kinds, self._f_colour: "color" in kinds,
            self._f_seq: "sequence" in kinds, self._f_dye: "dye" in kinds, self._f_marks: "marks" in kinds,
        }
        for w, on in shown.items():
            w.setVisible(on)
            self._rows[w].setVisible(on)
        self._form_host.setEnabled(bool(name))
        self._f_name.setText(name)
        self._f_ex.setText(str(e.get("excitation", "")))
        self._f_em.setText(str(e.get("emission", "")))
        self._f_words.setText(str(e.get("words", "")))
        colour = str(e.get("color", "")) or "#ffffff"
        self._f_colour.setIcon(_swatch(colour))
        self._f_colour.setText(colour)
        self._f_seq.set_text(str(e.get("sequence", "")), emit=False)
        self._f_dye.blockSignals(True)
        self._f_dye.clear()
        self._f_dye.addItem(NO_DYE)
        self._f_dye.addItems(self._editor.names("dyes"))
        self._f_dye.addItem(ADD_DYE)
        self._f_dye.setCurrentText(str(e.get("dye", "")) or NO_DYE)
        self._f_dye.blockSignals(False)
        self._f_marks.setCurrentText(str(e.get("marks", "")) or "chromosomes")

    def _set_selected(self, field: str, value: str) -> None:
        name = self._current_name()
        if not name or self._loading:
            return
        if str(self._entry(name).get(field, "")) == str(value).strip():
            return
        self._editor.set_value(self._key, name, field, str(value).strip())
        self._fill_table(select=name)

    def _on_form_name(self) -> None:
        name = self._current_name()
        if name and self._f_name.text().strip() and self._f_name.text().strip() != name:
            self._rename(name, self._f_name.text())

    def _on_form_sequence(self) -> None:
        path = self._f_seq.path()
        value = ""
        if path is not None:
            value = path.name if path.parent.resolve() == self._folder else str(path)
        self._set_selected("sequence", value)

    def _on_form_dye(self, _index: int) -> None:
        if self._f_dye.currentText() != ADD_DYE:
            text = self._f_dye.currentText()
            self._set_selected("dye", "" if text == NO_DYE else text)
            return
        dye, ok = QInputDialog.getText(self, "Add a dye", "Name of the dye:")
        if not ok or not dye.strip():
            self._fill_form()
            return
        colour = QColorDialog.getColor(QColor("#33ff33"), self, f"Colour of {dye.strip()}")
        try:
            self._editor.add("dyes", dye.strip(), color=colour.name() if colour.isValid() else "#ffffff",
                             words=dye.strip().lower())
        except ValueError as exc:
            QMessageBox.information(self, "Add a dye", str(exc))
        self._fill_nav()
        self._set_selected("dye", dye.strip())

    def _pick_colour(self) -> None:
        name = self._current_name()
        if not name:
            return
        chosen = QColorDialog.getColor(QColor(self._entry(name).get("color", "#ffffff")), self, "Colour")
        if chosen.isValid():
            self._set_selected("color", chosen.name())

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
        self._fill_table(select=name.strip())

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
        self._fill_table(select=keep)

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
