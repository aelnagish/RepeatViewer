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
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
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
    ("tm", "Tm Primer3 °C"), ("tm_supplier", "Tm supplier °C"), ("organism", "Organism"),
    ("target", "Target"), ("box", "Box"),
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
    ("ordered_by", "Ordered by:", "people"), ("ordered_on", "Ordered on:", "line"),
    ("aliases", "Other names:", "line"), ("dissolved", "", "check"),
    ("note", "Note:", "line"),
)
_MULTI_SKIP = ("name", "sequence", "aliases")


def _one_decimal(value) -> str:
    try:
        return f"{float(str(value).replace(',', '.')):.1f}"
    except (TypeError, ValueError):
        return str(value or "")


def lab_user_setting() -> str:
    return str(QSettings().value(LAB_USER_KEY, "") or "").strip()


def apply_lab_user_setting() -> None:
    inv.set_lab_user(lab_user_setting())


def _keeps_front(method):
    """After a file dialog or message box closes, macOS activates the main window
    and this (non-modal) window falls behind it: bring it back to the front."""
    import functools

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        try:
            return method(self, *args, **kwargs)
        finally:
            self.bring_to_front()

    return wrapper


def _keeps_front_of_parent(method):
    """Same as :func:`_keeps_front`, for small dialogs opened from this window."""
    import functools

    @functools.wraps(method)
    def wrapper(self, *args, **kwargs):
        try:
            return method(self, *args, **kwargs)
        finally:
            for delay in (0, 120, 350):
                QTimer.singleShot(delay, lambda w=self: (w.raise_(), w.activateWindow())
                                  if w.isVisible() else None)

    return wrapper


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
        self._person = QComboBox(self)
        self._person.setMinimumWidth(140)
        self._person.setToolTip("Show the primers ordered by one person")
        self._person.currentIndexChanged.connect(self._fill)
        top = QHBoxLayout()
        top.setSpacing(8)
        top.addWidget(self._filter, 1)
        top.addWidget(self._organism)
        top.addWidget(self._person)
        top.addWidget(self._status)

        self._where = QLabel(self)
        self._where.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._where.setWordWrap(False)
        self._where.setMinimumWidth(80)
        self._where.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self._change_btn = QPushButton("Change…", self)
        menu = QMenu(self._change_btn)
        # Run after the menu has closed: a file chooser opened while the menu is
        # still closing hands the focus back to the main window on macOS, and
        # this window drops behind it.
        later = lambda slot: (lambda: QTimer.singleShot(0, slot))  # noqa: E731
        menu.addAction("Use a Shared Folder (Sciebo)…", later(self._use_shared_folder))
        menu.addAction("Open Another Database…", later(self._open_other))
        menu.addAction("Use the Database on This Computer", later(self._use_local))
        menu.addSeparator()
        menu.addAction("Your Name in the History…",
                       later(lambda: self._ask_user(force=True)))
        self._change_btn.setMenu(menu)


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
        people = QPushButton("People…", self)
        people.setToolTip("The lab members: add people, set their role, mark who has left, "
                          "merge two entries for one person")
        people.clicked.connect(self._open_people)
        self._people_btn = people
        check = QPushButton("Check a Sequence…", self)
        check.setToolTip("Test every saved primer against a new sequence: which pairs "
                         "already give a product, where the tubes are")
        check.clicked.connect(lambda: self.checkSequenceRequested.emit(str(self._db.path)))
        close = QPushButton("Close", self)
        close.clicked.connect(self.close)
        bar = QHBoxLayout()
        bar.setSpacing(8)
        for b in (add, imp, exp, self._delete, dups, people, check):
            bar.addWidget(b)
        bar.addSpacing(12)
        bar.addWidget(self._where, 1)
        bar.addWidget(self._change_btn)
        bar.addSpacing(12)
        bar.addWidget(close)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.setSpacing(8)
        lay.addLayout(top)
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
            db = inv.open_database(path)
        except PrimerError:
            db = inv.open_database(None)
        if isinstance(db, inv.PrimerInventory):
            try:
                db.migrate_people()   # older "Anna & Nusrat" text -> two people
            except (OSError, PrimerError):
                pass
        return db

    def _switch_to(self, path: str) -> None:
        self._db = self._open(path or None)
        self._stats_cache.clear()
        self._refresh_where()
        self._refill_choices()
        self._fill()
        self._on_selection()
        self.databasePathChanged.emit(path)

    def _refresh_where(self) -> None:
        """Short line next to Change…; the full path is in the tooltip."""
        db = self._inventory()
        if db is None:
            self._where.setText("<span style='color:palette(mid)'>Database is an <b>older FASTA "
                                "file</b></span>")
            self._where.setToolTip(
                f"Database: {self._db.path}\nAn older single FASTA file: only names and "
                "sequences. Change… › Use the Database on This Computer, then Import this "
                "FASTA, to keep locations and stock.")
            return
        shared = db.path.resolve() != inv.default_inventory_folder().resolve()
        writers = db.writers()
        place = "a shared folder" if shared else "this computer"
        self._where.setText(f"<span style='color:palette(mid)'>Database is on <b>{place}</b>. "
                            f"User: <b>{db.user}</b></span>")
        self._where.setToolTip(
            f"Database: {db.path}\n"
            + ("A shared folder (Sciebo): changes from others appear after syncing."
               if shared else "On this computer only. Change… › Use a Shared Folder to "
               "share it with the lab.")
            + (f"\n{len(writers)} people have made changes." if len(writers) > 1 else ""))

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

    @_keeps_front
    def _ask_user(self, *, force: bool = False) -> bool:
        """Who is making changes (written into the history); asked once.

        Chosen from the lab's people; a new name is added as a person."""
        name = lab_user_setting()
        if name and not force:
            return True
        db = self._inventory()
        names = [p.name for p in db.people(include_former=False)] if db is not None else []
        start = name or inv.default_user()
        if start not in names:
            names = [start] + names
        text, ok = QInputDialog.getItem(
            self, "Who Are You?",
            "Your name, for the history and for \u201cOrdered by\u201d.\n"
            "Pick yourself, or type your name if you are new:",
            names, max(0, names.index(start)), True)
        text = (text or "").strip()
        if not ok or not text:
            return bool(name)
        QSettings().setValue(LAB_USER_KEY, text)
        apply_lab_user_setting()
        if db is not None:
            db.set_user(text)
            if db.person_by_name(text) is None:
                db.ensure_people([text])
            self._refill_choices()
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
                    d = dict(rec.data)
                    d["ordered_by"] = db.ordered_by_ids(rec)
                    d["_ordered_names"] = db.ordered_by_text(rec)
                    out.append((rec.rid, p, d))
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
        current = self._person.currentData()
        self._person.blockSignals(True)
        self._person.clear()
        self._person.addItem("Anyone", _ALL)
        db = self._inventory()
        for person in (db.people() if db is not None else []):
            n = len(db.primers_of(person.pid))
            self._person.addItem(f"{person.label} ({n})", person.pid)
        idx = self._person.findData(current)
        self._person.setCurrentIndex(max(0, idx))
        self._person.setEnabled(db is not None)
        self._person.blockSignals(False)
        self._people_btn.setEnabled(db is not None) if hasattr(self, "_people_btn") else None
        if "ordered_by" in self._fields:
            self._fields["ordered_by"].set_database(db)
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
        who = self._person.currentData()
        rows = self._rows()
        pairs = inv.partners([p for _i, p, _d in rows])
        shown = []
        for rid, p, d in rows:
            if org not in (None, _ALL) and str(d.get("organism", "")).strip() != org:
                continue
            if status not in (None, _ALL) and d.get("status", "") != status:
                continue
            if who not in (None, _ALL) and who not in (d.get("ordered_by") or []):
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
                # Text with one decimal (65.0, not 65); same width, so it sorts right.
                "tm": f"{st.tm:.1f}" if st else "",
                "tm_supplier": _one_decimal(d.get("tm_supplier", "")),
                "status": inv.STATUS_LABELS.get(d.get("status", ""), d.get("status", "")),
                "tested": inv.TESTED_LABELS.get(d.get("tested", ""), d.get("tested", "")),
                "ordered_by": d.get("_ordered_names", ""),
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
                if key in ("length", "tm", "tm_supplier"):
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                          | Qt.AlignmentFlag.AlignVCenter)
                if key == "tm":
                    item.setToolTip(f"Calculated here (Primer3, {self._thermo.describe()}).\n"
                                    "The supplier uses other conditions, so its Tm differs.")
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
        unit = "primer" if total == 1 else "primers"
        self.setWindowTitle(f"Primer Database: {len(shown)} of {total} {unit} shown"
                            if len(shown) != total else f"Primer Database: {total} {unit}")
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
        if kind == "people":
            w = PeopleField(self, before_add=self._before_write)
            w.changed.connect(self._form_touched)
            w.peopleAdded.connect(self._refill_choices)
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
        if isinstance(w, PeopleField):
            return w.value()
        if isinstance(w, QCheckBox):
            return w.isChecked()
        if isinstance(w, QComboBox):
            return str(w.currentData() or "") if not w.isEditable() else w.currentText().strip()
        return w.text().strip()

    def _set_field(self, key: str, value) -> None:
        w = self._fields[key]
        if isinstance(w, PeopleField):
            w.set_value(list(value) if isinstance(value, list) else [])
            return
        w.blockSignals(True)
        if isinstance(w, QCheckBox):
            w.setChecked(bool(value))
        elif isinstance(w, QComboBox) and not w.isEditable():
            idx = w.findData(value or "")
            w.setCurrentIndex(max(0, idx))
        elif isinstance(w, QComboBox):
            w.setEditText(str(value or ""))
            w.lineEdit().setCursorPosition(0)   # show the start of long values
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

    @_keeps_front
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
                         + f" · Tm {st.tm:.1f} °C (Primer3)"
                         + (f" ({st.tm_min:.1f}–{st.tm_max:.1f})" if st.degeneracy > 1 else "")
                         + f" · GC {st.gc:.0f} % · 3′ ΔG {st.end_dg:.1f}")
            sup = p.extra.get("tm_supplier")
            if sup:
                lines.append(f"  Tm on the synthesis report: {sup} °C (the supplier's "
                             "method and conditions)")
            lines.append(f"  Tm here: Primer3, {self._thermo.describe()}")
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
                lines.extend("  " + inv.describe_change(ch, db) for ch in reversed(rec.history[-30:]))
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

    @_keeps_front
    def _add(self) -> None:
        if not self._before_write():
            return
        dlg = AddPrimersDialog(self._db, self._thermo, self)
        if dlg.exec() != int(QDialog.DialogCode.Accepted):
            return
        primers = dlg.primers()
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

    @_keeps_front
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

    @_keeps_front
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

    @_keeps_front
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

    @_keeps_front
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

    @_keeps_front
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

    @_keeps_front
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

    @_keeps_front
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

    @_keeps_front
    def _open_people(self) -> None:
        db = self._inventory()
        if db is None:
            return
        if not self._before_write():
            return
        dlg = PeopleDialog(db, self)
        dlg.exec()
        self._refill_choices()
        self._fill()
        self._on_selection()

    def bring_to_front(self) -> None:
        def _go() -> None:
            try:
                if self.isVisible():
                    self.raise_()
                    self.activateWindow()
            except RuntimeError:
                pass

        _go()
        # macOS re-activates the main window a moment after a native panel closes
        # (its closing animation): raise again after it.
        for delay in (0, 120, 350, 700):
            QTimer.singleShot(delay, _go)

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


# --------------------------------------------------------------------------
# Add primers
# --------------------------------------------------------------------------


def _rows_from_text(text: str) -> list[tuple[str, str]]:
    """(name, sequence) pairs from FASTA, ``name<TAB>sequence`` lines or bare sequences.
    Sequences are not checked here."""
    body = str(text or "").strip()
    if not body:
        return []
    out: list[tuple[str, str]] = []
    if body.startswith(">"):
        name, seq = None, []
        for line in body.splitlines():
            line = line.strip()
            if line.startswith(">"):
                if name is not None:
                    out.append((name, "".join(seq)))
                name, seq = line[1:].split()[0] if line[1:].split() else "", []
            elif line:
                seq.append(line)
        if name is not None:
            out.append((name, "".join(seq)))
        return out
    import re as _re

    for raw in body.splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p for p in _re.split(r"\t+|;|\s{2,}", line) if p.strip()]
        if len(parts) == 1 and " " in line:
            head, _sp, tail = line.rpartition(" ")
            parts = [head, tail]
        if len(parts) >= 2:
            out.append((parts[0].strip(), parts[1].strip()))
        else:
            out.append(("", parts[0].strip()))
    return out


class AddPrimersDialog(QDialog):
    """Add one or more primers: a row per primer with live checks, and the
    details they share (organism, target, who ordered them, where they go)."""

    _COLS = ("Name", "Sequence 5′→3′", "Length", "Tm Primer3 °C", "GC %", "Check")

    def __init__(self, db: PrimerDatabase, thermo: ThermoSettings, parent: QWidget) -> None:
        super().__init__(parent)
        self.setWindowTitle("Add Primers")
        self.resize(860, 600)
        self._db = db
        self._thermo = thermo
        self._inv = db if isinstance(db, inv.PrimerInventory) else None
        self._checking = False

        intro = QLabel("One row per primer. Paste several at once (from a sheet, an order "
                       "confirmation or FASTA) with <b>Paste</b>. Lower-case 5′ letters are "
                       "a tail that does not bind (e.g. an M13 or T7 adapter).", self)
        intro.setWordWrap(True)
        self._table = QTableWidget(0, len(self._COLS), self)
        self._table.setHorizontalHeaderLabels(list(self._COLS))
        self._table.verticalHeader().setVisible(False)
        self._table.verticalHeader().setDefaultSectionSize(24)
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
        hh.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        for c in (2, 3, 4):
            hh.setSectionResizeMode(c, QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(5, QHeaderView.ResizeMode.Interactive)
        self._table.setColumnWidth(0, 170)
        self._table.setColumnWidth(5, 220)
        self._table.itemChanged.connect(self._check_row)
        add_row = QPushButton("Add Row", self)
        add_row.clicked.connect(lambda: self._add_row())
        add_pair = QPushButton("Add Pair", self)
        add_pair.setToolTip("Two rows: …_F and …_R")
        add_pair.clicked.connect(self._add_pair)
        paste = QPushButton("Paste", self)
        paste.setToolTip("Rows from the clipboard: name<TAB>sequence lines (copied from a "
                         "sheet), FASTA, or one sequence per line")
        paste.clicked.connect(self._paste)
        remove = QPushButton("Remove Row", self)
        remove.clicked.connect(self._remove_rows)
        rows_bar = QHBoxLayout()
        for b in (add_row, add_pair, paste, remove):
            rows_bar.addWidget(b)
        rows_bar.addStretch(1)

        # Shared details
        box = QGroupBox("Details for these primers", self)
        grid = QGridLayout(box)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(4)

        def combo(key: str) -> QComboBox:
            w = QComboBox(self)
            w.setEditable(True)
            w.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
            values = sorted({str((p.extra or {}).get(key, "")).strip() for p in db.primers}
                            - {""}, key=str.lower)
            w.addItems(values)
            w.setEditText("")
            return w

        self._organism = combo("organism")
        self._target = combo("target")
        self._box = combo("box")
        self._freezer = combo("freezer")
        self._room = combo("room")
        self._ordered_by = PeopleField(self)
        self._ordered_by.set_database(self._inv)
        if self._inv is not None:
            me = self._inv.person_by_name(inv.lab_user())
            if me is not None and me.active:
                self._ordered_by.set_value([me.pid])
        self._ordered_on = QLineEdit(time.strftime("%Y-%m-%d"), self)
        self._status = QComboBox(self)
        for v in inv.STATUS_VALUES:
            self._status.addItem(inv.STATUS_LABELS[v] or "–", v)
        self._status.setCurrentIndex(self._status.findData("ordered"))
        self._dissolved = QCheckBox("Dissolved (water added to the tube)", self)
        self._note = QLineEdit(self)
        self._note.setPlaceholderText("e.g. project, supplier, order number")
        fields = (("Organism:", self._organism), ("Target:", self._target),
                  ("Ordered by:", self._ordered_by), ("Ordered on:", self._ordered_on),
                  ("Status:", self._status), ("Box:", self._box),
                  ("Freezer:", self._freezer), ("Room:", self._room))
        for i, (label, w) in enumerate(fields):
            r, c = divmod(i, 2)
            lab = QLabel(label, self)
            lab.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            grid.addWidget(lab, r, c * 2)
            grid.addWidget(w, r, c * 2 + 1)
        r = len(fields) // 2
        lab = QLabel("Note:", self)
        lab.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        grid.addWidget(lab, r, 0)
        grid.addWidget(self._note, r, 1, 1, 3)
        grid.addWidget(self._dissolved, r + 1, 1, 1, 3)
        grid.setColumnStretch(1, 1)
        grid.setColumnStretch(3, 1)
        if self._inv is None:
            box.setEnabled(False)
            box.setTitle("Details (only in a database folder, not in an older FASTA)")

        self._summary = QLabel(self)
        cancel = QPushButton("Cancel", self)
        cancel.clicked.connect(self.reject)
        self._save = QPushButton("Save", self)
        self._save.setDefault(True)
        self._save.clicked.connect(self.accept)
        bar = QHBoxLayout()
        bar.addWidget(self._summary, 1)
        bar.addWidget(cancel)
        bar.addWidget(self._save)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.setSpacing(8)
        lay.addWidget(intro)
        lay.addWidget(self._table, 1)
        lay.addLayout(rows_bar)
        lay.addWidget(box)
        lay.addLayout(bar)
        self._add_pair()

    # ---------------------------------------------------------------- rows

    def _add_row(self, name: str = "", seq: str = "") -> int:
        r = self._table.rowCount()
        self._checking = True
        self._table.insertRow(r)
        for c in range(len(self._COLS)):
            item = QTableWidgetItem("")
            if c >= 2:
                item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
            if c == 1:
                item.setFont(_mono_font())
            if c in (2, 3, 4):
                item.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                      | Qt.AlignmentFlag.AlignVCenter)
            self._table.setItem(r, c, item)
        self._table.item(r, 0).setText(name)
        self._table.item(r, 1).setText(seq)
        self._checking = False
        self._check_row(self._table.item(r, 1))
        return r

    def _add_pair(self) -> None:
        base = ""
        names = [self._table.item(r, 0).text().strip() for r in range(self._table.rowCount())]
        if names and names[-1]:
            base = inv.pair_stem(names[-1])[0]
        r = self._add_row(f"{base}_F" if base else "", "")
        self._add_row(f"{base}_R" if base else "", "")
        self._table.setCurrentCell(r, 0 if not base else 1)
        self._table.editItem(self._table.item(r, 0 if not base else 1))

    def _remove_rows(self) -> None:
        rows = sorted({i.row() for i in self._table.selectedIndexes()}, reverse=True)
        if not rows and self._table.rowCount():
            rows = [self._table.rowCount() - 1]
        for r in rows:
            self._table.removeRow(r)
        self._update_summary()

    def _paste(self) -> None:
        from PySide6.QtWidgets import QApplication

        self.paste_text(QApplication.clipboard().text())

    def paste_text(self, text: str) -> int:
        """Rows from pasted text; a line with a bad sequence still becomes a row,
        marked in Check, so nothing pasted is silently lost."""
        rows = _rows_from_text(text)
        if not rows:
            QMessageBox.information(self, "Paste", "The clipboard holds no primers.")
            return 0
        empty = [r for r in range(self._table.rowCount())
                 if not self._table.item(r, 0).text().strip()
                 and not self._table.item(r, 1).text().strip()]
        for name, seq in rows:
            if empty:
                r = empty.pop(0)
                self._table.item(r, 0).setText(name)
                self._table.item(r, 1).setText(seq)
            else:
                self._add_row(name, seq)
        for r in reversed(empty):
            self._table.removeRow(r)
        return len(rows)

    def _check_row(self, item: QTableWidgetItem | None) -> None:
        if self._checking or item is None or item.column() > 1:
            return
        r = item.row()
        self._checking = True
        try:
            name = self._table.item(r, 0).text().strip()
            raw = self._table.item(r, 1).text().strip()
            problems, stats = [], None
            if raw:
                try:
                    p = Primer(name or "x", raw)
                    stats = analyze_oligo(p.sequence, self._thermo)
                    same = self._db.find_sequence(p.sequence)
                    if same is not None:
                        problems.append(f"already saved as {same.name}")
                    if stats.warnings:
                        problems.append(stats.warnings[0])
                except PrimerError as exc:
                    problems.append(str(exc).split(": ", 1)[-1])
            if raw and not name:
                problems.insert(0, "needs a name")
            if name and self._db.get(name) is not None:
                problems.insert(0, f"name taken: saved as {name}_2")
            self._table.item(r, 2).setText(str(stats.binding_length) if stats else "")
            self._table.item(r, 3).setText(f"{stats.tm:.1f}" if stats else "")
            self._table.item(r, 4).setText(f"{stats.gc:.0f}" if stats else "")
            check = self._table.item(r, 5)
            bad = any(x.startswith(("needs", "already")) for x in problems) or bool(
                raw and stats is None)
            check.setText("; ".join(problems) if problems else ("OK" if raw else ""))
            check.setToolTip("; ".join(problems))
            check.setForeground(Qt.GlobalColor.darkRed if bad else
                                (Qt.GlobalColor.darkYellow if problems else
                                 Qt.GlobalColor.darkGreen))
        finally:
            self._checking = False
        self._update_summary()

    def _valid_rows(self) -> list[tuple[str, str]]:
        out = []
        for r in range(self._table.rowCount()):
            name = self._table.item(r, 0).text().strip()
            seq = self._table.item(r, 1).text().strip()
            if not name or not seq:
                continue
            try:
                p = Primer(name, seq)
            except PrimerError:
                continue
            if self._db.find_sequence(p.sequence) is not None:
                continue
            out.append((name, seq))
        return out

    def _update_summary(self) -> None:
        n = len(self._valid_rows())
        filled = sum(1 for r in range(self._table.rowCount())
                     if self._table.item(r, 1).text().strip())
        skip = filled - n
        self._summary.setText(f"<b>{n}</b> primer{'s' if n != 1 else ''} to add"
                              + (f" <span style='color:#b45309'>· {skip} row(s) not added: "
                                 "see Check</span>" if skip > 0 else ""))
        self._save.setEnabled(n > 0)
        self._save.setText(f"Add {n}" if n else "Save")

    # ---------------------------------------------------------------- result

    def details(self) -> dict:
        d = {"organism": self._organism.currentText().strip(),
             "target": self._target.currentText().strip(),
             "box": self._box.currentText().strip(),
             "freezer": self._freezer.currentText().strip(),
             "room": self._room.currentText().strip(),
             "ordered_by": self._ordered_by.value(),
             "ordered_on": self._ordered_on.text().strip(),
             "status": str(self._status.currentData() or ""),
             "dissolved": self._dissolved.isChecked(),
             "note": self._note.text().strip()}
        return {k: v for k, v in d.items() if v not in ("", [], None, False)}

    def primers(self) -> list[Primer]:
        shared = self.details() if self._inv is not None else {}
        note = shared.pop("note", "") if self._inv is not None else self._note.text().strip()
        out = []
        for name, seq in self._valid_rows():
            out.append(Primer(name, seq, note=note, source="added", extra=dict(shared)))
        return out


# --------------------------------------------------------------------------
# People
# --------------------------------------------------------------------------


class PeopleField(QWidget):
    """``Ordered by``: one or more lab members, chosen from the list of people."""

    changed = Signal()
    peopleAdded = Signal()

    def __init__(self, parent: QWidget, *, before_add=None) -> None:
        super().__init__(parent)
        self._db: inv.PrimerInventory | None = None
        self._ids: list[str] = []
        self._before_add = before_add
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)
        self._text = QLineEdit(self)
        self._text.setReadOnly(True)
        self._text.setPlaceholderText("Nobody chosen")
        self._btn = QPushButton("Choose…", self)
        self._btn.clicked.connect(self._choose)
        lay.addWidget(self._text, 1)
        lay.addWidget(self._btn)

    def set_database(self, db: inv.PrimerInventory | None) -> None:
        self._db = db
        self._show()

    def value(self) -> list[str]:
        return list(self._ids)

    def set_value(self, ids: list[str]) -> None:
        self._ids = list(ids)
        self._show()

    def _show(self) -> None:
        names = self._db.names_of(self._ids) if self._db is not None else list(self._ids)
        self._text.setText(", ".join(names))
        self._text.setCursorPosition(0)

    def _choose(self) -> None:
        if self._db is None:
            return
        dlg = PeoplePicker(self._db, self._ids, self, before_add=self._before_add)
        if dlg.exec() != int(QDialog.DialogCode.Accepted):
            if dlg.added:
                self.peopleAdded.emit()
            return
        if dlg.added:
            self.peopleAdded.emit()
        ids = dlg.chosen()
        if ids != self._ids:
            self._ids = ids
            self._show()
            self.changed.emit()


class PeoplePicker(QDialog):
    """Tick the people who ordered the primer; add someone new."""

    def __init__(self, db: inv.PrimerInventory, chosen: list[str], parent: QWidget,
                 *, before_add=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Ordered by")
        self.resize(320, 420)
        self._db = db
        self._before_add = before_add
        self.added = False
        self._chosen = list(chosen)
        self._filter = QLineEdit(self)
        self._filter.setPlaceholderText("Filter")
        self._filter.setClearButtonEnabled(True)
        self._filter.textChanged.connect(self._fill)
        self._list = QListWidget(self)
        self._list.itemChanged.connect(self._ticked)
        add = QPushButton("Add Person…", self)
        add.clicked.connect(self._add)
        cancel = QPushButton("Cancel", self)
        save = QPushButton("Save", self)
        save.setDefault(True)
        cancel.clicked.connect(self.reject)
        save.clicked.connect(self.accept)
        bar = QHBoxLayout()
        bar.addWidget(add)
        bar.addStretch(1)
        bar.addWidget(cancel)
        bar.addWidget(save)
        lay = QVBoxLayout(self)
        lay.addWidget(self._filter)
        lay.addWidget(self._list, 1)
        lay.addLayout(bar)
        self._fill()

    def _fill(self, *_a) -> None:
        text = self._filter.text().strip().lower()
        self._list.blockSignals(True)
        self._list.clear()
        for person in self._db.people():
            if not person.active and person.pid not in self._chosen:
                continue
            if text and text not in person.name.lower():
                continue
            item = QListWidgetItem(person.label + (f"  ·  {person.role}" if person.role else ""))
            item.setData(Qt.ItemDataRole.UserRole, person.pid)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked if person.pid in self._chosen
                               else Qt.CheckState.Unchecked)
            self._list.addItem(item)
        self._list.blockSignals(False)

    def _ticked(self, item: QListWidgetItem) -> None:
        pid = item.data(Qt.ItemDataRole.UserRole)
        if item.checkState() == Qt.CheckState.Checked and pid not in self._chosen:
            self._chosen.append(pid)
        elif item.checkState() != Qt.CheckState.Checked and pid in self._chosen:
            self._chosen.remove(pid)

    def _add(self) -> None:
        if self._before_add is not None and not self._before_add():
            return
        name, ok = QInputDialog.getText(self, "Add Person", "Name of the new lab member:")
        if not ok or not name.strip():
            return
        try:
            pid = self._db.add_person(name)
        except PrimerError as exc:
            existing = self._db.person_by_name(name)
            if existing is None:
                QMessageBox.warning(self, "Add Person", str(exc))
                return
            pid = existing.pid
        self.added = True
        if pid not in self._chosen:
            self._chosen.append(pid)
        self._filter.clear()
        self._fill()

    def chosen(self) -> list[str]:
        return list(self._chosen)


class PeopleDialog(QDialog):
    """Manage the lab's people: add, edit, mark as former member, merge, remove."""

    _COLS = ("Name", "Initials", "Role", "Email", "Primers", "Member")

    def __init__(self, db: inv.PrimerInventory, parent: QWidget) -> None:
        super().__init__(parent)
        self.setWindowTitle("People")
        self.resize(720, 520)
        self._db = db
        self._table = QTableWidget(0, len(self._COLS), self)
        self._table.setHorizontalHeaderLabels(list(self._COLS))
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.verticalHeader().setVisible(False)
        self._table.setAlternatingRowColors(True)
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.itemSelectionChanged.connect(self._on_select)

        form = QGridLayout()
        form.setHorizontalSpacing(8)
        self._name = QLineEdit(self)
        self._initials = QLineEdit(self)
        self._initials.setMaxLength(6)
        self._role = QComboBox(self)
        self._role.setEditable(True)
        self._role.addItems(list(inv.PERSON_ROLES))
        self._email = QLineEdit(self)
        self._active = QCheckBox("Current lab member", self)
        for r, (label, w) in enumerate((("Name:", self._name), ("Initials:", self._initials),
                                        ("Role:", self._role), ("Email:", self._email))):
            lab = QLabel(label, self)
            lab.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            form.addWidget(lab, r, 0)
            form.addWidget(w, r, 1)
        form.addWidget(self._active, 4, 1)
        self._primers = QLabel(self)
        self._primers.setWordWrap(True)
        self._primers.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        form.addWidget(self._primers, 5, 1)
        self._save = QPushButton("Save Changes", self)
        self._save.clicked.connect(self._save_person)
        form.addWidget(self._save, 6, 1, alignment=Qt.AlignmentFlag.AlignRight)

        add = QPushButton("Add Person…", self)
        add.clicked.connect(self._add)
        self._merge = QPushButton("Merge…", self)
        self._merge.setToolTip("Two entries for one person (e.g. \u201cAnna\u201d and "
                               "\u201cAnna K.\u201d): select both, keep one")
        self._merge.clicked.connect(self._merge_selected)
        self._remove = QPushButton("Remove", self)
        self._remove.setToolTip("Someone with primers is kept as a former member, so the "
                                "record of who ordered what stays")
        self._remove.clicked.connect(self._remove_selected)
        close = QPushButton("Close", self)
        close.clicked.connect(self.accept)
        bar = QHBoxLayout()
        for b in (add, self._merge, self._remove):
            bar.addWidget(b)
        bar.addStretch(1)
        bar.addWidget(close)
        lay = QVBoxLayout(self)
        lay.addWidget(self._table, 1)
        lay.addLayout(form)
        lay.addLayout(bar)
        self._fill()

    def _fill(self, select: list[str] | None = None) -> None:
        people = self._db.people()
        self._table.blockSignals(True)
        self._table.setRowCount(len(people))
        for i, p in enumerate(people):
            n = len(self._db.primers_of(p.pid))
            for j, v in enumerate((p.name, p.initials, p.role, p.email, n,
                                   "Yes" if p.active else "Former")):
                item = QTableWidgetItem()
                item.setData(Qt.ItemDataRole.DisplayRole, v)
                item.setData(Qt.ItemDataRole.UserRole, p.pid)
                if not p.active:
                    item.setForeground(Qt.GlobalColor.gray)
                self._table.setItem(i, j, item)
        self._table.resizeColumnsToContents()
        if select:
            for r in range(self._table.rowCount()):
                if self._table.item(r, 0).data(Qt.ItemDataRole.UserRole) in select:
                    self._table.selectRow(r)
        self._table.blockSignals(False)
        self._on_select()

    def _selected(self) -> list[str]:
        return [self._table.item(i.row(), 0).data(Qt.ItemDataRole.UserRole)
                for i in self._table.selectionModel().selectedRows()]

    def _on_select(self) -> None:
        ids = self._selected()
        one = len(ids) == 1
        for w in (self._name, self._initials, self._role, self._email, self._active,
                  self._save):
            w.setEnabled(one)
        self._merge.setEnabled(len(ids) == 2)
        self._remove.setEnabled(bool(ids))
        if not one:
            for w in (self._name, self._initials, self._email):
                w.clear()
            self._role.setEditText("")
            self._primers.clear()
            return
        p = self._db.persons[ids[0]]
        self._name.setText(p.name)
        self._initials.setText(p.initials)
        self._role.setEditText(p.role)
        self._email.setText(p.email)
        self._active.setChecked(p.active)
        recs = sorted(r.name for r in self._db.primers_of(p.pid))
        self._primers.setText(
            f"<span style='color:palette(mid)'>Ordered {len(recs)} primer(s)"
            + (": " + ", ".join(recs[:12]) + (" …" if len(recs) > 12 else "") if recs else "")
            + "</span>")

    @_keeps_front_of_parent
    def _save_person(self) -> None:
        ids = self._selected()
        if len(ids) != 1:
            return
        try:
            self._db.update_person(ids[0], {
                "name": self._name.text(), "initials": self._initials.text().strip(),
                "role": self._role.currentText().strip(), "email": self._email.text().strip(),
                "active": self._active.isChecked()})
        except PrimerError as exc:
            QMessageBox.warning(self, "People", str(exc))
            return
        self._fill(select=ids)

    @_keeps_front_of_parent
    def _add(self) -> None:
        name, ok = QInputDialog.getText(self, "Add Person", "Name of the new lab member:")
        if not ok or not name.strip():
            return
        try:
            pid = self._db.add_person(name)
        except PrimerError as exc:
            QMessageBox.warning(self, "Add Person", str(exc))
            return
        self._fill(select=[pid])

    @_keeps_front_of_parent
    def _merge_selected(self) -> None:
        ids = self._selected()
        if len(ids) != 2:
            return
        a, b = (self._db.persons[i] for i in ids)
        keep_name, ok = QInputDialog.getItem(
            self, "Merge", f"{a.name} and {b.name} are one person. Keep which name?",
            [a.name, b.name], 0, False)
        if not ok:
            return
        keep, drop = (a, b) if keep_name == a.name else (b, a)
        n = self._db.merge_people(keep.pid, drop.pid)
        QMessageBox.information(self, "Merge", f"{drop.name} merged into {keep.name}; "
                                f"{n} primer(s) now show {keep.name}.")
        self._fill(select=[keep.pid])

    @_keeps_front_of_parent
    def _remove_selected(self) -> None:
        ids = self._selected()
        if not ids:
            return
        names = ", ".join(self._db.persons[i].name for i in ids)
        ans = QMessageBox.question(
            self, "Remove", f"Remove {names}?\n\nSomeone who ordered primers stays as a "
            "former member (shown in grey), so the record of who ordered what is kept.")
        if ans != QMessageBox.StandardButton.Yes:
            return
        for pid in ids:
            self._db.remove_person(pid)
        self._fill()


__all__ = ["AddPrimersDialog", "PeopleDialog", "PeopleField", "PrimerDatabaseDialog", "apply_lab_user_setting",
           "lab_user_setting"]
