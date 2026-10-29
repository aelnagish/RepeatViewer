"""Tools > Primers and Probes > Database: the lab's primer inventory.

One window for everything the lab knows about a primer: sequence and checks
(Tm, hairpin, dimers), what it amplifies and in which organism, where the tube
is (box, freezer, room), whether it is in stock and whether it worked, who
ordered it and when, and the history of every change.

The database can be a folder on this computer or a shared folder (Sciebo):
see :mod:`repeatviewer.backends.primers.inventory`.
"""

from __future__ import annotations

import time
from pathlib import Path

from PySide6.QtCore import QSettings, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QDialog,
    QFileDialog,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSplitter,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import inventory as inv
from repeatviewer.backends.primers.oligos import (
    KIND_LABELS,
    PRIMER_KINDS,
    Primer,
    PrimerDatabase,
    load_primers,
    parse_primer_text,
    write_primers,
)
from repeatviewer.backends.primers.thermo import (
    PrimerError,
    ThermoSettings,
    analyze_oligo,
    heterodimer,
)
from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox

def _mono_font():
    from PySide6.QtGui import QFontDatabase

    return QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)


LAB_USER_KEY = "primers/lab_user"
_ALL = "__all__"
_COLUMNS: tuple[tuple[str, str], ...] = (
    ("name", "Name"), ("pair", "Pair"), ("sequence", "Sequence 5′→3′"), ("length", "Length"),
    ("tm", "Tm °C"), ("organism", "Organism"), ("target", "Target"), ("box", "Box"),
    ("freezer", "Freezer"), ("room", "Room"), ("status", "Status"), ("tested", "Tested"),
    ("ordered_by", "Ordered by"), ("ordered_on", "Ordered on"),
)
_COL = {k: i for i, (k, _l) in enumerate(_COLUMNS)}
# Form: (field, label, widget kind). Kinds: line, seq, combo (editable, filled
# from the values in use), kind, status, tested, check.
_FORM: tuple[tuple[str, str, str], ...] = (
    ("name", "Name:", "line"), ("sequence", "Sequence:", "seq"),
    ("kind", "Direction:", "kind"), ("pair", "Pairs with:", "line"),
    ("organism", "Organism:", "combo"), ("target", "Target:", "combo"),
    ("box", "Box:", "combo"), ("freezer", "Freezer:", "combo"),
    ("room", "Room:", "combo"), ("status", "Status:", "status"),
    ("tested", "Tested:", "tested"), ("anneal_temp", "Annealing °C:", "line"),
    ("tm_supplier", "Tm (supplier):", "line"), ("concentration", "Concentration:", "line"),
    ("ordered_by", "Ordered by:", "combo"), ("ordered_on", "Ordered on:", "line"),
    ("aliases", "Other names:", "line"), ("dissolved", "", "check"),
    ("note", "Note:", "line"),
)
_MULTI_SKIP = ("name", "sequence", "aliases")


def lab_user_setting() -> str:
    return str(QSettings().value(LAB_USER_KEY, "") or "").strip()


def apply_lab_user_setting() -> None:
    inv.set_lab_user(lab_user_setting())


class PrimerDatabaseDialog(QDialog):
    """The primer inventory window (non-modal)."""

    databaseChanged = Signal()
    databasePathChanged = Signal(str)      # "" = the database on this computer
    checkSequenceRequested = Signal(str)   # database path

    def __init__(self, path: Path | None = None, *, thermo: ThermoSettings | None = None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Primer Database")
        self.setModal(False)
        self.setWindowModality(Qt.WindowModality.NonModal)
        self.resize(1180, 720)
        apply_lab_user_setting()
        self._thermo = thermo or ThermoSettings()
        self._stats_cache: dict[tuple, object] = {}
        self._form_ids: list[str] = []          # records shown in the form
        self._form_original: dict = {}
        self._filling = False
        self._db: PrimerDatabase = self._open(path)

        # ---- top: filters and where the database is
        self._filter = QLineEdit(self)
        self._filter.setPlaceholderText("Filter: name, sequence, target, box, note, who …")
        self._filter.setClearButtonEnabled(True)
        self._filter.textChanged.connect(self._fill)
        self._organism = QComboBox(self)
        self._organism.setMinimumWidth(180)
        self._organism.currentIndexChanged.connect(self._fill)
        self._status = QComboBox(self)
        self._status.addItem("Any status", _ALL)
        for v in inv.STATUS_VALUES[1:]:
            self._status.addItem(inv.STATUS_LABELS[v], v)
        self._status.currentIndexChanged.connect(self._fill)
        top = QHBoxLayout()
        top.setSpacing(8)
        top.addWidget(self._filter, 1)
        top.addWidget(self._organism)
        top.addWidget(self._status)

        self._where = QLabel(self)
        self._where.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._where.setWordWrap(True)
        self._change_btn = QPushButton("Change…", self)
        menu = QMenu(self._change_btn)
        menu.addAction("Use a Shared Folder (Sciebo)…", self._use_shared_folder)
        menu.addAction("Open Another Database…", self._open_other)
        menu.addAction("Use the Database on This Computer", self._use_local)
        menu.addSeparator()
        menu.addAction("Your Name in the History…", lambda: self._ask_user(force=True))
        self._change_btn.setMenu(menu)
        where_row = QHBoxLayout()
        where_row.addWidget(self._where, 1)
        where_row.addWidget(self._change_btn)

        # ---- table
        self._table = QTableWidget(0, len(_COLUMNS), self)
        self._table.setHorizontalHeaderLabels([lab for _k, lab in _COLUMNS])
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.verticalHeader().setVisible(False)
        self._table.verticalHeader().setDefaultSectionSize(22)
        self._table.setAlternatingRowColors(True)
        self._table.setSortingEnabled(True)
        self._table.setWordWrap(False)
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        hh.setStretchLastSection(True)
        hh.setSortIndicator(0, Qt.SortOrder.AscendingOrder)
        self._table.itemSelectionChanged.connect(self._on_selection)

        # ---- form (left) and checks / history (right)
        self._form_host = QWidget(self)
        grid = QGridLayout(self._form_host)
        grid.setContentsMargins(0, 0, 8, 0)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(4)
        self._fields: dict[str, QWidget] = {}
        row = col = 0
        for key, label, kind in _FORM:
            w = self._make_field(key, kind)
            self._fields[key] = w
            wide = kind == "seq" or key in ("note", "aliases")
            if wide and col == 1:
                row, col = row + 1, 0
            lab = QLabel(label, self._form_host)
            lab.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            grid.addWidget(lab, row, col * 2)
            grid.addWidget(w, row, col * 2 + 1, 1, 3 if wide else 1)
            if wide:
                row, col = row + 1, 0
            else:
                col += 1
                if col == 2:
                    row, col = row + 1, 0
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(3, 1)
        self._form_title = QLabel(self)
        self._save_btn = QPushButton("Save Changes", self)
        self._save_btn.setDefault(False)
        self._save_btn.clicked.connect(self._save_form)
        self._revert_btn = QPushButton("Revert", self)
        self._revert_btn.clicked.connect(self._on_selection)
        fb = QHBoxLayout()
        fb.addWidget(self._form_title, 1)
        fb.addWidget(self._revert_btn)
        fb.addWidget(self._save_btn)
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setWidget(self._form_host)
        left = QWidget(self)
        ll = QVBoxLayout(left)
        ll.setContentsMargins(0, 0, 0, 0)
        ll.addLayout(fb)
        ll.addWidget(scroll, 1)

        self._details = QPlainTextEdit(self)
        self._details.setReadOnly(True)
        self._details.setFont(_mono_font())
        self._details.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
        self._details.setPlaceholderText(
            "Select a primer to see its Tm, GC, hairpin, self-dimer and history.\n"
            "Select two to see the dimer between them.\n"
            "Select several to change their box, status … at once.")
        bottom = QSplitter(Qt.Orientation.Horizontal, self)
        bottom.addWidget(left)
        bottom.addWidget(self._details)
        bottom.setStretchFactor(0, 3)
        bottom.setStretchFactor(1, 2)
        bottom.setSizes([660, 480])
        split = QSplitter(Qt.Orientation.Vertical, self)
        split.addWidget(self._table)
        split.addWidget(bottom)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 2)

        # ---- buttons
        add = QPushButton("Add…", self)
        add.setToolTip("Paste primers: name<TAB>sequence lines, FASTA, or bare sequences")
        add.clicked.connect(self._add)
        imp = QPushButton("Import…", self)
        imp.setToolTip("The lab sheet (CSV, TSV or Excel .xlsx, e.g. downloaded from Google "
                       "Sheets) or a primers FASTA")
        imp.clicked.connect(self._import)
        exp = QPushButton("Export…", self)
        exp.setToolTip("Table with every field (.tsv), FASTA, or an order sheet")
        exp.clicked.connect(self._export)
        self._delete = QPushButton("Delete", self)
        self._delete.clicked.connect(self._delete_selected)
        dups = QPushButton("Find Duplicates", self)
        dups.setToolTip("Same sequence under two names, reverse complements, one primer "
                        "inside another")
        dups.clicked.connect(self._find_duplicates)
        check = QPushButton("Check a Sequence…", self)
        check.setToolTip("Test every saved primer against a new sequence: which pairs "
                         "already give a product, where the tubes are")
        check.clicked.connect(lambda: self.checkSequenceRequested.emit(str(self._db.path)))
        close = QPushButton("Close", self)
        close.clicked.connect(self.close)
        bar = QHBoxLayout()
        bar.setSpacing(8)
        for b in (add, imp, exp, self._delete, dups, check):
            bar.addWidget(b)
        bar.addStretch(1)
        bar.addWidget(close)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.setSpacing(8)
        lay.addLayout(top)
        lay.addLayout(where_row)
        lay.addWidget(split, 1)
        lay.addLayout(bar)

        self._poll = QTimer(self)
        self._poll.setInterval(4000)
        self._poll.timeout.connect(self._poll_disk)
        self._poll.start()
        self._refresh_where()
        self._refill_choices()
        self._fill()
        self._on_selection()

    # ------------------------------------------------------------ database

    def database(self) -> PrimerDatabase:
        return self._db

    def _inventory(self) -> inv.PrimerInventory | None:
        return self._db if isinstance(self._db, inv.PrimerInventory) else None

    @staticmethod
    def _open(path) -> PrimerDatabase:
        try:
            return inv.open_database(path)
        except PrimerError:
            return inv.open_database(None)

    def _switch_to(self, path: str) -> None:
        self._db = self._open(path or None)
        self._stats_cache.clear()
        self._refresh_where()
        self._refill_choices()
        self._fill()
        self._on_selection()
        self.databasePathChanged.emit(path)

    def _refresh_where(self) -> None:
        db = self._inventory()
        if db is None:
            self._where.setText(
                f"<span style='color:palette(mid)'>Database: {self._db.path} (an older "
                "single FASTA file: only names and sequences). Change… › Use the Database "
                "on This Computer, then Import this FASTA, to keep locations and stock."
                "</span>")
            return
        shared = db.path.resolve() != inv.default_inventory_folder().resolve()
        people = db.people()
        who = (f" · {len(people)} people have made changes" if len(people) > 1 else "")
        kind = "shared folder" if shared else "on this computer only"
        self._where.setText(
            f"<span style='color:palette(mid)'>Database: {db.path} ({kind}{who}) · "
            f"you are <b>{db.user}</b></span>")

    def _poll_disk(self) -> None:
        db = self._inventory()
        if db is None or not self.isVisible():
            return
        if not db.changed_on_disk():
            return
        if self._form_dirty():
            self._form_title.setText("<span style='color:#b45309'>Someone changed the "
                                     "database; save or revert to see it.</span>")
            return
        keep = self._form_ids
        db.reload()
        self._stats_cache.clear()
        self._refresh_where()
        self._refill_choices()
        self._fill(select_ids=keep)

    def _ask_user(self, *, force: bool = False) -> bool:
        """The name written into the history; asked once."""
        name = lab_user_setting()
        if name and not force:
            return True
        text, ok = QInputDialog.getText(
            self, "Your Name",
            "Your name, as it should appear in the history and in “Ordered by”:",
            text=name or inv.default_user())
        if not ok or not text.strip():
            return bool(name)
        QSettings().setValue(LAB_USER_KEY, text.strip())
        apply_lab_user_setting()
        db = self._inventory()
        if db is not None:
            db.set_user(text.strip())
        self._refresh_where()
        return True

    def _before_write(self) -> bool:
        db = self._inventory()
        if db is None:
            return True
        if not self._ask_user():
            return False
        if db.user != inv.lab_user():
            db.set_user(inv.lab_user())
        return True

    # ------------------------------------------------------------ table

    def _stats(self, p: Primer):
        key = (p.sequence, repr(self._thermo.to_dict()))
        if key not in self._stats_cache:
            try:
                self._stats_cache[key] = analyze_oligo(p.sequence, self._thermo)
            except Exception:  # noqa: BLE001
                self._stats_cache[key] = None
        return self._stats_cache[key]

    def _rows(self) -> list[tuple[str, Primer, dict]]:
        """(id, primer, fields) for every primer in the database."""
        db = self._inventory()
        out = []
        if db is not None:
            for p in db.primers:
                rec = db.records.get(p.extra.get("id", ""))
                if rec is not None:
                    out.append((rec.rid, p, dict(rec.data)))
        else:
            for p in self._db.primers:
                out.append((p.name, p, {"name": p.name, "sequence": p.sequence,
                                        "kind": p.kind, "note": p.note}))
        return out

    def _refill_choices(self) -> None:
        rows = self._rows()
        current = self._organism.currentData()
        self._organism.blockSignals(True)
        self._organism.clear()
        self._organism.addItem("All organisms", _ALL)
        for org in sorted({str(d.get("organism", "")).strip() for _i, _p, d in rows} - {""},
                          key=str.lower):
            self._organism.addItem(org, org)
        idx = self._organism.findData(current)
        self._organism.setCurrentIndex(max(0, idx))
        self._organism.blockSignals(False)
        for key, _lab, kind in _FORM:
            if kind != "combo":
                continue
            box = self._fields[key]
            values = sorted({str(d.get(key, "")).strip() for _i, _p, d in rows} - {""},
                            key=str.lower)
            text = box.currentText()
            box.blockSignals(True)
            box.clear()
            box.addItems(values)
            box.setEditText(text)
            box.blockSignals(False)

    def _fill(self, *_a, select_ids: list[str] | None = None) -> None:
        if select_ids is None:
            select_ids = list(self._form_ids)
        text = self._filter.text().strip().lower()
        org = self._organism.currentData()
        status = self._status.currentData()
        rows = self._rows()
        pairs = inv.partners([p for _i, p, _d in rows])
        shown = []
        for rid, p, d in rows:
            if org not in (None, _ALL) and str(d.get("organism", "")).strip() != org:
                continue
            if status not in (None, _ALL) and d.get("status", "") != status:
                continue
            if text:
                hay = " ".join(str(v) for v in d.values()).lower() + " " + p.sequence.lower()
                if text not in hay:
                    continue
            shown.append((rid, p, d))
        self._filling = True
        self._table.blockSignals(True)
        self._table.setSortingEnabled(False)
        self._table.setRowCount(len(shown))
        mono = _mono_font()
        for i, (rid, p, d) in enumerate(shown):
            st = self._stats(p)
            vals = {
                "name": p.name, "pair": ", ".join(pairs.get(p.name, [])),
                "sequence": p.sequence, "length": len(p.binding),
                "tm": round(st.tm, 1) if st else "",
                "status": inv.STATUS_LABELS.get(d.get("status", ""), d.get("status", "")),
                "tested": inv.TESTED_LABELS.get(d.get("tested", ""), d.get("tested", "")),
            }
            for key, _lab in _COLUMNS:
                v = vals[key] if key in vals else str(d.get(key, "") or "")
                item = QTableWidgetItem()
                item.setData(Qt.ItemDataRole.DisplayRole, v)
                item.setData(Qt.ItemDataRole.UserRole, rid)
                if key == "sequence":
                    item.setFont(mono)
                    if st and st.warnings:
                        item.setToolTip("; ".join(st.warnings))
                if key in ("length", "tm"):
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                          | Qt.AlignmentFlag.AlignVCenter)
                if key == "status" and d.get("status") in ("used up", "running low"):
                    item.setForeground(Qt.GlobalColor.darkRed)
                if key == "tested" and d.get("tested") == "failed":
                    item.setForeground(Qt.GlobalColor.darkRed)
                self._table.setItem(i, _COL[key], item)
        self._table.setSortingEnabled(True)
        hh = self._table.horizontalHeader()
        self._table.sortByColumn(hh.sortIndicatorSection(), hh.sortIndicatorOrder())
        self._table.resizeColumnsToContents()
        self._table.setColumnWidth(_COL["sequence"],
                                   max(200, min(330, self._table.columnWidth(_COL["sequence"]))))
        for key in ("organism", "target", "freezer"):
            c = _COL[key]
            self._table.setColumnWidth(c, min(200, self._table.columnWidth(c)))
        wanted = set(select_ids or [])
        if wanted:
            self._table.clearSelection()
            mode = self._table.selectionMode()
            self._table.setSelectionMode(QAbstractItemView.SelectionMode.MultiSelection)
            for r in range(self._table.rowCount()):
                item = self._table.item(r, 0)
                if item is not None and item.data(Qt.ItemDataRole.UserRole) in wanted:
                    self._table.selectRow(r)
            self._table.setSelectionMode(mode)
        self._table.blockSignals(False)
        self._filling = False
        total = len(rows)
        self._delete.setEnabled(bool(shown))
        self.setWindowTitle(f"Primer Database ({len(shown)} of {total})" if len(shown) != total
                            else f"Primer Database ({total})")
        if wanted:
            self._on_selection()

    def _selected_ids(self) -> list[str]:
        ids = []
        for idx in self._table.selectionModel().selectedRows():
            item = self._table.item(idx.row(), 0)
            if item is not None:
                ids.append(item.data(Qt.ItemDataRole.UserRole))
        return ids

    def _primer_for(self, rid: str) -> Primer | None:
        db = self._inventory()
        if db is not None:
            rec = db.records.get(rid)
            return rec.to_primer() if rec is not None else None
        return self._db.get(rid)

    def _selected(self) -> list[Primer]:
        return [p for rid in self._selected_ids() if (p := self._primer_for(rid)) is not None]

    # ------------------------------------------------------------ form

    def _make_field(self, key: str, kind: str) -> QWidget:
        if kind == "check":
            w = QCheckBox("Dissolved (water added to the tube)", self)
            w.setTristate(False)
            w.toggled.connect(self._form_touched)
            return w
        if kind in ("kind", "status", "tested"):
            w = QComboBox(self)
            if kind == "kind":
                w.addItem("From the name (…_F / …_rev)", "")
                for k in PRIMER_KINDS:
                    w.addItem(KIND_LABELS[k], k)
            else:
                values = inv.STATUS_VALUES if kind == "status" else inv.TESTED_VALUES
                labels = inv.STATUS_LABELS if kind == "status" else inv.TESTED_LABELS
                for v in values:
                    w.addItem(labels[v] or "–", v)
            w.currentIndexChanged.connect(self._form_touched)
            return w
        if kind == "combo":
            w = QComboBox(self)
            w.setEditable(True)
            w.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
            w.lineEdit().setClearButtonEnabled(True)
            w.editTextChanged.connect(self._form_touched)
            return w
        w = QLineEdit(self)
        w.setClearButtonEnabled(True)
        if kind == "seq":
            w.setFont(_mono_font())
            w.setPlaceholderText("5′→3′; lower-case 5′ letters are a tail")
        if key == "pair":
            w.setPlaceholderText("Automatic from the names (…_for / …_rev)")
        if key == "ordered_on":
            w.setPlaceholderText(time.strftime("%Y-%m-%d"))
        w.textChanged.connect(self._form_touched)
        return w

    def _get_field(self, key: str):
        w = self._fields[key]
        if isinstance(w, QCheckBox):
            return w.isChecked()
        if isinstance(w, QComboBox):
            return str(w.currentData() or "") if not w.isEditable() else w.currentText().strip()
        return w.text().strip()

    def _set_field(self, key: str, value) -> None:
        w = self._fields[key]
        w.blockSignals(True)
        if isinstance(w, QCheckBox):
            w.setChecked(bool(value))
        elif isinstance(w, QComboBox) and not w.isEditable():
            idx = w.findData(value or "")
            w.setCurrentIndex(max(0, idx))
        elif isinstance(w, QComboBox):
            w.setEditText(str(value or ""))
        else:
            w.setText(str(value or ""))
            w.setCursorPosition(0)
        w.blockSignals(False)

    def _form_values(self) -> dict:
        return {key: self._get_field(key) for key, _l, _k in _FORM}

    def _form_dirty(self) -> bool:
        return bool(self._form_ids) and self._form_values() != self._form_original

    def _form_touched(self, *_a) -> None:
        dirty = self._form_dirty()
        self._save_btn.setEnabled(dirty)
        self._revert_btn.setEnabled(dirty)

    def _on_selection(self, *_a) -> None:
        if self._filling:
            return
        ids = self._selected_ids()
        self._form_ids = ids
        db = self._inventory()
        rows = {rid: d for rid, _p, d in self._rows()}
        datas = [rows[i] for i in ids if i in rows]
        for key, _l, _k in _FORM:
            vals = {repr(d.get(key, False if key == "dissolved" else "")) for d in datas}
            value = datas[0].get(key, "") if len(vals) == 1 and datas else ""
            self._set_field(key, value)
        full = db is not None
        for key, _l, _k in _FORM:
            editable = bool(datas) and (full or key in ("name", "sequence", "kind", "note"))
            if len(datas) > 1 and key in _MULTI_SKIP:
                editable = False
            self._fields[key].setEnabled(editable)
        if not datas:
            self._form_title.setText("<span style='color:palette(mid)'>Select a primer to "
                                     "edit it</span>")
        elif len(datas) == 1:
            self._form_title.setText(f"<b>{datas[0].get('name', '')}</b>")
        else:
            self._form_title.setText(f"<b>{len(datas)} primers</b> <span style='color:"
                                     "palette(mid)'>· fields you change are set on all of "
                                     "them</span>")
        self._form_original = self._form_values()
        self._form_touched()
        self._show_details()

    def _save_form(self) -> None:
        if not self._form_dirty():
            return
        if not self._before_write():
            return
        now = self._form_values()
        diff = {k: v for k, v in now.items() if v != self._form_original.get(k)}
        db = self._inventory()
        try:
            if db is not None:
                if len(self._form_ids) > 1:
                    diff = {k: v for k, v in diff.items() if k not in _MULTI_SKIP}
                for rid in self._form_ids:
                    db.update(rid, diff)
            else:
                p = self._db.get(self._form_ids[0])
                if p is None:
                    return
                if "name" in diff:
                    if self._db.get(diff["name"]) is not None:
                        raise PrimerError(f"Another primer is already called {diff['name']}.")
                    p.name = diff["name"]
                if "sequence" in diff:
                    p.sequence = Primer(p.name, diff["sequence"]).sequence
                if "kind" in diff:
                    p.kind = diff["kind"] or ""
                    p.__post_init__()
                if "note" in diff:
                    p.note = diff["note"]
                self._db.save()
                self._form_ids = [p.name]
        except PrimerError as exc:
            QMessageBox.warning(self, "Primer Database", str(exc))
            return
        self.databaseChanged.emit()
        self._stats_cache.clear()
        self._refill_choices()
        self._fill(select_ids=list(self._form_ids))

    # ------------------------------------------------------------ details

    def _show_details(self) -> None:
        sel = self._selected()
        if not sel:
            self._details.clear()
            return
        lines: list[str] = []
        db = self._inventory()
        pairs = inv.partners(self._db.primers)
        for p in sel[:2]:
            st = self._stats(p)
            if st is None:
                lines.append(f"{p.name}: cannot analyse {p.sequence}")
                continue
            lines.append(f"{p.name}  ({KIND_LABELS.get(p.kind, p.kind)})")
            lines.append(f"  5′-{p.sequence}-3′")
            lines.append(f"  {st.binding_length} nt binding"
                         + (f" + {len(p.tail)} nt tail" if p.tail else "")
                         + f" · Tm {st.tm:.1f} °C"
                         + (f" ({st.tm_min:.1f}–{st.tm_max:.1f})" if st.degeneracy > 1 else "")
                         + f" · GC {st.gc:.0f} % · 3′ ΔG {st.end_dg:.1f}")
            sup = p.extra.get("tm_supplier")
            if sup:
                lines.append(f"  Tm on the synthesis report: {sup} °C (other method)")
            if pairs.get(p.name):
                lines.append(f"  Pairs with: {', '.join(pairs[p.name])}")
            if st.warnings:
                lines.append("  ⚠ " + "; ".join(st.warnings))
            for s in (st.hairpin, st.self_dimer):
                if s.found:
                    tm = f", Tm {s.tm:.1f} °C" if s.tm > 0 else ""
                    lines.append(f"  {s.kind}: ΔG {s.dg:.1f} kcal/mol{tm}")
                    lines.extend("    " + x for x in s.lines)
                else:
                    lines.append(f"  {s.kind}: none")
            lines.append("")
        if len(sel) >= 2:
            d = heterodimer(sel[0].sequence, sel[1].sequence, self._thermo)
            lines.append(f"{sel[0].name} × {sel[1].name}: "
                         + (f"cross-dimer ΔG {d.dg:.1f} kcal/mol" if d.found else "no dimer"))
            lines.extend("    " + x for x in d.lines)
            ta, tb = self._stats(sel[0]), self._stats(sel[1])
            if ta and tb:
                lines.append(f"    Tm difference {abs(ta.tm - tb.tm):.1f} °C")
            lines.append("")
        if db is not None and len(sel) == 1:
            rec = db.records.get(sel[0].extra.get("id", ""))
            if rec is not None and rec.history:
                lines.append("History")
                lines.extend("  " + inv.describe_change(ch) for ch in reversed(rec.history[-30:]))
        self._details.setPlainText("\n".join(lines))

    # ------------------------------------------------------------ actions

    def add_primers(self, primers: list[Primer]) -> int:
        if not self._before_write():
            return 0
        added = self._db.add(primers)
        if added and not isinstance(self._db, inv.PrimerInventory):
            self._db.save()
        if added:
            self.databaseChanged.emit()
            self._refill_choices()
            ids = [p.extra.get("id", p.name) for p in added]
            self._fill(select_ids=ids)
        return len(added)

    def _add(self) -> None:
        dlg = QDialog(self)
        dlg.setWindowTitle("Add Primers")
        dlg.resize(560, 320)
        text = QPlainTextEdit(dlg)
        text.setFont(_mono_font())
        text.setPlaceholderText(
            "CroSat2_F\tGACTTGCAAGTTCGAGAAC\nCroSat2_R\tGTTGCATCCACTTAGCTTAG\n\n"
            "or FASTA, or one sequence per line.\n"
            "Lower-case 5′ letters are a tail (e.g. caggaaacagctatgacGTTGCATCC…).\n"
            "After adding, fill in box, organism … below the table.")
        kind = QComboBox(dlg)
        kind.addItem("Direction from the name (…_F / …_rev)", "")
        for k in PRIMER_KINDS:
            kind.addItem(KIND_LABELS[k], k)
        cancel = QPushButton("Cancel", dlg)
        save = QPushButton("Save", dlg)
        save.setDefault(True)
        cancel.clicked.connect(dlg.reject)
        save.clicked.connect(dlg.accept)
        bar = QHBoxLayout()
        bar.addWidget(kind)
        bar.addStretch(1)
        bar.addWidget(cancel)
        bar.addWidget(save)
        lay = QVBoxLayout(dlg)
        lay.addWidget(text, 1)
        lay.addLayout(bar)
        if dlg.exec() != int(QDialog.DialogCode.Accepted):
            return
        try:
            primers = parse_primer_text(text.toPlainText())
        except PrimerError as exc:
            QMessageBox.warning(self, "Add Primers", str(exc))
            return
        forced = kind.currentData()
        if forced:
            for p in primers:
                p.kind = forced
        n = self.add_primers(primers)
        if n < len(primers):
            QMessageBox.information(self, "Add Primers",
                                    f"Added {n} of {len(primers)}; the others were already "
                                    "in the database (same sequence).")

    def import_file(self, path: Path, *, confirm: bool = True) -> int:
        """Lab sheet (CSV / TSV / XLSX) or primers FASTA. Returns primers added."""
        path = Path(path)
        if path.suffix.lower() in (".csv", ".tsv", ".xlsx", ".xlsm"):
            plan = inv.plan_import(inv.read_table(path), self._db)
            if confirm:
                labels = ", ".join(f"{inv.FIELD_LABELS.get(k, k)}" for k in
                                   plan.columns.values())
                text = (f"{path.name}: {plan.summary()}.\n\nColumns read: {labels}.")
                details = []
                if plan.renamed:
                    details.append("Already in under another name (the sheet name is kept "
                                   "as Other names):")
                    details += [f"  {a} = {b}" for a, b in plan.renamed]
                if plan.problems:
                    details.append("Skipped:")
                    details += ["  " + x for x in plan.problems]
                if plan.notes:
                    details.append("To check:")
                    details += ["  " + x for x in plan.notes]
                box = QMessageBox(QMessageBox.Icon.Question, "Import", text,
                                  QMessageBox.StandardButton.Cancel, self)
                go = box.addButton("Import", QMessageBox.ButtonRole.AcceptRole)
                if details:
                    box.setDetailedText("\n".join(details))
                box.exec()
                if box.clickedButton() is not go:
                    return 0
            if not self._before_write():
                return 0
            n = inv.apply_import(plan, self._db)
            self.databaseChanged.emit()
            self._refill_choices()
            self._fill()
            return n
        primers = load_primers(path)
        return self.add_primers(primers)

    def _import(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Import primers", str(Path.home()),
            "Lab sheet or primers (*.csv *.tsv *.xlsx *.fa *.fasta *.txt);;All files (*)")
        if not path:
            return
        try:
            n = self.import_file(Path(path))
        except PrimerError as exc:
            QMessageBox.warning(self, "Import", str(exc))
            return
        if n:
            QMessageBox.information(self, "Import", f"Added {n} primer(s) from "
                                                    f"{Path(path).name}.")

    def _export(self) -> None:
        ids = self._selected_ids()
        path, chosen = QFileDialog.getSaveFileName(
            self, "Export primers", str(Path.home() / "primers.tsv"),
            "Table, every field (*.tsv);;FASTA (*.fa);;Order sheet (*.tsv)")
        if not path:
            return
        p = Path(path)
        sel = self._selected() if ids else list(self._db.primers)
        if chosen.startswith("FASTA") or p.suffix.lower() in (".fa", ".fasta"):
            write_primers(p, sel)
        elif chosen.startswith("Order"):
            lines = ["name\tsequence\tlength\tscale\tpurification"]
            lines += [f"{pr.name}\t{pr.sequence}\t{len(pr.sequence)}\t25 nmol\tdesalted"
                      for pr in sel]
            p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        else:
            rows = inv.export_rows(self._db)
            if ids:
                names = {pr.name for pr in sel}
                rows = [rows[0]] + [r for r in rows[1:] if r[0] in names]
            p.write_text("\n".join("\t".join(c.replace("\t", " ") for c in r) for r in rows)
                         + "\n", encoding="utf-8")

    def _delete_selected(self) -> None:
        sel = self._selected()
        if not sel:
            return
        names = ", ".join(p.name for p in sel[:5]) + (" …" if len(sel) > 5 else "")
        ok = QMessageBox.question(self, "Delete primers",
                                  f"Delete {len(sel)} primer(s) from the database?\n{names}")
        if ok != QMessageBox.StandardButton.Yes:
            return
        if not self._before_write():
            return
        self._db.remove([p.name for p in sel])
        if not isinstance(self._db, inv.PrimerInventory):
            self._db.save()
        self.databaseChanged.emit()
        self._form_ids = []
        self._fill(select_ids=[])
        self._on_selection()

    def _find_duplicates(self) -> None:
        dups = inv.find_duplicates(self._db.primers)
        if not dups:
            QMessageBox.information(self, "Find Duplicates", "No duplicates: every sequence "
                                    "is there once, and none lies inside another.")
            return
        text = "\n".join(f"{d.kind}: {d.a} and {d.b}" for d in dups)
        box = QMessageBox(QMessageBox.Icon.Information, "Find Duplicates",
                          f"{len(dups)} found. Select them in the table to compare; keep one "
                          "and note the other name under Other names.", parent=self)
        box.setDetailedText(text)
        box.exec()
        names = {d.a for d in dups} | {d.b for d in dups}
        self._filter.clear()
        self._fill(select_ids=[p.extra.get("id", p.name) for p in self._db.primers
                               if p.name in names])

    # ------------------------------------------------------------ location

    def _use_shared_folder(self) -> None:
        start = str(Path.home() / "Sciebo") if (Path.home() / "Sciebo").is_dir() \
            else str(Path.home())
        folder = QFileDialog.getExistingDirectory(
            self, "Choose the shared folder (for example in your Sciebo folder)", start)
        if not folder:
            return
        try:
            target = inv.create_inventory(Path(folder))
        except OSError as exc:
            QMessageBox.warning(self, "Shared Folder", f"Cannot use {folder}:\n{exc}")
            return
        self._move_to(target)

    def _open_other(self) -> None:
        folder = QFileDialog.getExistingDirectory(
            self, f"Choose a “{inv.FOLDER_NAME}” folder", str(Path.home()))
        if not folder:
            return
        if not inv.is_inventory_folder(folder):
            QMessageBox.information(self, "Open Database",
                                    f"{folder} is not a primer database folder. Use "
                                    "“Use a Shared Folder” to make one.")
            return
        self._switch_to(folder)

    def _use_local(self) -> None:
        self._switch_to("")

    def _move_to(self, target: Path) -> None:
        if not self._ask_user():
            return
        dst = inv.PrimerInventory(target)
        same = Path(self._db.path).resolve() == target.resolve()
        if not same and self._db.primers:
            ans = QMessageBox.question(
                self, "Shared Folder",
                f"Copy the {len(self._db.primers)} primers of the current database into "
                f"{target}?\n\nPrimers already there (same sequence) are not copied twice.")
            if ans == QMessageBox.StandardButton.Yes:
                inv.copy_into(self._db, dst)
        self._switch_to(str(target))
        QMessageBox.information(
            self, "Shared Folder",
            f"The Primer Database is now\n{target}\n\nEveryone in the lab: Database › "
            "Change… › Open Another Database…, and choose this folder in their own Sciebo "
            "folder. Changes appear on the other computers once Sciebo has synced.")

    def closeEvent(self, event) -> None:  # noqa: N802
        if self._form_dirty():
            ans = QMessageBox.question(
                self, "Primer Database", "Save your changes to the selected primer(s)?",
                QMessageBox.StandardButton.Save | QMessageBox.StandardButton.Discard
                | QMessageBox.StandardButton.Cancel)
            if ans == QMessageBox.StandardButton.Cancel:
                event.ignore()
                return
            if ans == QMessageBox.StandardButton.Save:
                self._save_form()
        self._poll.stop()
        super().closeEvent(event)


__all__ = ["PrimerDatabaseDialog", "apply_lab_user_setting", "lab_user_setting"]
