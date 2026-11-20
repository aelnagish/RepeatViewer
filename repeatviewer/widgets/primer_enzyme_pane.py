"""Enzymes tab of the primer workbench (Primers rail): restriction map on the sequence.

Pick an enzyme set and how often an enzyme may cut; ticked enzymes are drawn on
the bases (a line at each cut). With a Bismark .cov or modkit bedMethyl file,
sites methylated in a context that blocks the enzyme are drawn "probably not
cut" (red, dashed). Hover a site for the position, the reads and the reason.
Digest opens the Digest page (tandem-array ladders, methylation, partial cutting).

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import json
import threading
from pathlib import Path

from PySide6.QtCore import QObject, QSettings, Qt, QTimer, Signal
from PySide6.QtGui import QColor
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import enzymes as ez
from repeatviewer.widgets.primer_design_pane import _row, _rule, _section, _spin
from repeatviewer.widgets.primer_pcr_pane import FIELD_W, PCR_LABEL_W

_SETTINGS_KEY = "primers/workbench_enzymes"
_RECENT_KEY = "primers/workbench_meth_files"
_BROWSE = "__browse__"

SET_HELP = (
    "<b>Set</b>: which enzymes to look at. <b>Methylation-sensitive</b> holds the enzymes "
    "with known CG or CHG sensitivity (HpaII and MspI as a pair).<br>"
    "<b>Show</b>: enzymes cutting this often are ticked and drawn on the sequence. Tick or "
    "untick any enzyme in the list below."
)
METH_HELP = (
    "<b>Data</b>: a Bismark <b>.cov</b> or modkit <b>bedMethyl</b> file. Rows of the "
    "chromosome named like this sequence are read (an extracted region named "
    "<i>chrom:start-end</i> is shifted to its own positions).<br>"
    "<b>Blocked when</b>: a cytosine of the site (or of a context overlapping it) is "
    "methylated at least this much, with at least this many reads, in a context that "
    "blocks the enzyme. Such sites are drawn red and dashed: <i>probably not cut</i>.<br>"
    "Sensitivity per enzyme (CG, CHG, Dam, Dcm) from the Promega and NEB tables and REBASE "
    "(McClelland et al. 1994). Dam and Dcm only matter for DNA grown in E. coli. Without "
    "data the map shows the labels only."
)
LIST_HELP = (
    "<b>Cuts</b>: sites on this sequence. <b>In your data</b>: sites probably not cut "
    "because of methylation in your file.<br>Click an enzyme to go to its first site. "
    "<b>Digest</b> shows the ticked enzymes as a gel: monomer ladders, sites lost "
    "in some copies, methylated sites from your data and partial cutting."
)
PAIR_NOTE = ("HpaII and MspI both cut CCGG; only HpaII is blocked by methylated CG, so the "
             "pair shows methylation on a gel.")


class _Relay(QObject):
    done = Signal(int, object, object, object)   # generation, rows, MethData, error


def _combo(items, parent: QWidget, current=None) -> QComboBox:
    box = QComboBox(parent)
    for label, data in items:
        box.addItem(label, data)
    if current is not None:
        box.setCurrentIndex(max(0, box.findData(current)))
    box.setFixedWidth(FIELD_W)
    return box


class EnzymePane(QWidget):
    """``panel`` is the PrimersPanel (host, open file, sequence)."""

    def __init__(self, panel, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._panel = panel
        self._rows: list[ez.EnzymeRow] = []
        self._listed: list[ez.EnzymeRow] = []
        self._meth: ez.MethData | None = None
        self._ticks: dict[str, bool] = {}      # the user's own ticks (this sequence)
        self._generation = 0
        self.active = False
        self._relay = _Relay(self)
        self._relay.done.connect(self._on_done)
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(300)
        self._timer.timeout.connect(self._run)
        s = self._load_settings()

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        self.stack = QStackedWidget(self)
        outer.addWidget(self.stack)
        self.map_page = QWidget(self.stack)
        self.stack.addWidget(self.map_page)
        from repeatviewer.widgets.primer_digest_pane import DigestPage

        self.digest_page = DigestPage(self, self.stack)
        self.stack.addWidget(self.digest_page)
        self.digest_page.backRequested.connect(self.show_map_page)
        lay = QVBoxLayout(self.map_page)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(10)

        sec, sl, _ = _section("Enzyme set", self, "Enzyme set", SET_HELP)
        self.set_combo = _combo([(ez.SET_LABELS[k], k) for k in ez.SETS], sec, s.enzyme_set)
        sl.addLayout(_row("Set", self.set_combo, parent=sec, label_w=PCR_LABEL_W))
        self.show_combo = _combo([(ez.SHOW_LABELS[k], k) for k in ez.SHOWS], sec, s.show)
        sl.addLayout(_row("Show", self.show_combo, parent=sec, label_w=PCR_LABEL_W))
        lay.addWidget(sec)
        lay.addWidget(_rule(self))

        sec, sl, _ = _section("Methylation", self, "Methylation", METH_HELP)
        self.data_combo = QComboBox(sec)
        self.data_combo.setFixedWidth(FIELD_W)
        sl.addLayout(_row("Data", self.data_combo, parent=sec, label_w=PCR_LABEL_W))
        self.min_pct = _spin(0, 100, s.min_percent, suffix=" %", parent=sec)
        sl.addLayout(_row("Methylated ≥", self.min_pct, parent=sec, label_w=PCR_LABEL_W))
        self.min_reads = _spin(1, 10000, s.min_reads, parent=sec)
        sl.addLayout(_row("Reads ≥", self.min_reads, parent=sec, label_w=PCR_LABEL_W))
        self.data_note = QLabel(sec)
        self.data_note.setWordWrap(True)
        self.data_note.setTextFormat(Qt.TextFormat.RichText)
        sl.addWidget(self.data_note)
        lay.addWidget(sec)
        lay.addWidget(_rule(self))

        sec, sl, self.list_title = _section("Enzymes", self, "Enzymes", LIST_HELP)
        self.table = QTableWidget(0, 4, sec)
        self.table.setHorizontalHeaderLabels(["Enzyme", "Site", "Cuts", "In your data"])
        self.table.verticalHeader().hide()
        self.table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        head = self.table.horizontalHeader()
        head.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        head.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        head.setStretchLastSection(True)
        for c, w in enumerate((84, 78, 44)):
            self.table.setColumnWidth(c, w)
        self.table.setMinimumHeight(180)
        self.table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        from repeatviewer.widgets.primer_design_pane import _card

        sl.addWidget(_card(self.table, sec), 1)
        self.pair_note = QLabel(PAIR_NOTE, sec)
        self.pair_note.setWordWrap(True)
        self.pair_note.setEnabled(False)
        sl.addWidget(self.pair_note)
        self.select_all_button = QPushButton("Select All", sec)
        self.select_all_button.setToolTip(
            "Tick every enzyme in the list; press again to untick them all")
        self.digest_button = QPushButton("Digest", sec)
        self.digest_button.setToolTip("Digest with the ticked enzymes")
        br = QHBoxLayout()
        br.addStretch(1)
        br.addWidget(self.select_all_button)
        br.addWidget(self.digest_button)
        sl.addLayout(br)
        lay.addWidget(sec, 1)

        self._fill_data_combo(s.meth_path)
        for box in (self.set_combo, self.show_combo):
            box.currentIndexChanged.connect(self._on_set_changed)
        self.data_combo.currentIndexChanged.connect(self._on_data_changed)
        self.min_pct.valueChanged.connect(self.schedule)
        self.min_reads.valueChanged.connect(self.schedule)
        self.table.itemChanged.connect(self._on_item_changed)
        self.table.itemSelectionChanged.connect(self._on_row_selected)
        self.digest_button.clicked.connect(self._on_digest)
        self.select_all_button.clicked.connect(self._on_select_all)
        self._set_note("", None)
        self.pair_note.hide()
        self.digest_button.setEnabled(False)

    # ------------------------------------------------------------------ settings

    def _load_settings(self) -> ez.MapSettings:
        raw = QSettings().value(_SETTINGS_KEY, "")
        try:
            return ez.MapSettings.from_dict(json.loads(raw) if raw else None)
        except (TypeError, ValueError):
            return ez.MapSettings()

    def settings(self) -> ez.MapSettings:
        data = self.data_combo.currentData()
        return ez.MapSettings(enzyme_set=str(self.set_combo.currentData()),
                              show=str(self.show_combo.currentData()),
                              meth_path="" if data in (None, _BROWSE) else str(data),
                              min_percent=float(self.min_pct.value()),
                              min_reads=int(self.min_reads.value()))

    def _recent(self) -> list[str]:
        raw = QSettings().value(_RECENT_KEY, "")
        try:
            out = json.loads(raw) if raw else []
        except (TypeError, ValueError):
            out = []
        return [p for p in out if isinstance(p, str) and Path(p).is_file()]

    def _remember(self, path: str) -> None:
        if not path:
            return
        recent = [path] + [p for p in self._recent() if p != path]
        QSettings().setValue(_RECENT_KEY, json.dumps(recent[:10]))

    def _fill_data_combo(self, current: str = "") -> None:
        box = self.data_combo
        box.blockSignals(True)
        box.clear()
        box.addItem("None", "")
        paths: list[str] = []
        for p in [str(x) for x in ez.methylation_files_near(self._panel.fasta_path())] + \
                self._recent() + ([current] if current and Path(current).is_file() else []):
            if p not in paths:
                paths.append(p)
        for p in paths:
            box.addItem(Path(p).name, p)
            box.setItemData(box.count() - 1, p, Qt.ItemDataRole.ToolTipRole)
        box.addItem("Browse…", _BROWSE)
        i = box.findData(current) if current else 0
        box.setCurrentIndex(max(0, i))
        box.blockSignals(False)

    def reset_for_sequence(self) -> None:
        self._ticks = {}
        self.show_map_page()
        self._fill_data_combo(self.settings().meth_path)
        self.schedule()

    # ------------------------------------------------------------------ live run

    def _on_set_changed(self, *_a) -> None:
        self._ticks = {}
        self.schedule()

    def _on_data_changed(self, *_a) -> None:
        if self.data_combo.currentData() == _BROWSE:
            start = str(self._panel.fasta_path().parent) if self._panel.fasta_path() else ""
            path, _ = QFileDialog.getOpenFileName(
                self, "Methylation File", start,
                "Methylation (*.cov *.cov.gz *.bed *.bed.gz *.bedmethyl *.bedmethyl.gz);;"
                "All files (*)")
            self._fill_data_combo(path or "")
            if path:
                self._remember(path)
        self.schedule()

    def schedule(self, *_a) -> None:
        if self.active:
            self._timer.start()

    def _run(self) -> None:
        self._generation += 1
        gen = self._generation
        s = self.settings()
        QSettings().setValue(_SETTINGS_KEY, json.dumps(s.to_dict()))
        seq = self._panel._template
        seq_id = self._panel._seq_id
        if not seq:
            self._rows, self._meth = [], None
            self._fill()
            return
        if s.meth_path:
            self._remember(s.meth_path)
            self._set_note(f"Reading {Path(s.meth_path).name}…", None)
        relay = self._relay

        def work() -> None:
            try:
                meth = ez.load_methylation(s.meth_path, seq_id, len(seq)) if s.meth_path else None
                rows = ez.restriction_map(seq, s, meth)
            except Exception as exc:  # noqa: BLE001 - shown in the pane
                relay.done.emit(gen, None, None, str(exc))
                return
            relay.done.emit(gen, rows, meth, None)

        threading.Thread(target=work, daemon=True).start()

    def _on_done(self, gen: int, rows, meth, error) -> None:
        if gen != self._generation:
            return
        if error:
            self._rows, self._meth = [], None
            self._fill()
            self._set_note(f"Restriction map failed: {error}", "amber")
            return
        self._rows, self._meth = list(rows or []), meth
        self._fill()
        if self.stack.currentWidget() is self.digest_page:
            self.digest_page.schedule()

    # ------------------------------------------------------------------ table

    def _ticked_default(self, row: ez.EnzymeRow) -> bool:
        return ez.in_range(row.cuts, str(self.show_combo.currentData())) and row.cuts > 0

    def is_ticked(self, name: str) -> bool:
        row = next((r for r in self._rows if r.name == name), None)
        if row is None:
            return False
        return self._ticks.get(name, self._ticked_default(row))

    def ticked_rows(self) -> list[ez.EnzymeRow]:
        return [r for r in self._listed if self.is_ticked(r.name) and r.cuts > 0]

    def _fill(self) -> None:
        show = str(self.show_combo.currentData())
        self._listed = ez.listed_rows(self._rows, show)
        has_data = self._meth is not None and bool(self._meth.calls)
        self.table.blockSignals(True)
        self.table.setRowCount(len(self._listed))
        grey = QColor(140, 140, 146)
        amber = QColor(180, 83, 9)
        for r, row in enumerate(self._listed):
            name = QTableWidgetItem(row.name)
            name.setFlags(name.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            name.setCheckState(Qt.CheckState.Checked if self.is_ticked(row.name)
                               else Qt.CheckState.Unchecked)
            name.setData(Qt.ItemDataRole.UserRole, row.name)
            name.setToolTip(ez.sensitivity_text(row.name))
            data = "–"
            if has_data and row.cuts:
                data = f"{row.blocked} blocked"
            items = [name, QTableWidgetItem(row.site), QTableWidgetItem(str(row.cuts)),
                     QTableWidgetItem(data)]
            inside = ez.in_range(row.cuts, show)
            for c, it in enumerate(items):
                if not inside:
                    it.setForeground(grey)
                elif c == 3 and row.blocked:
                    it.setForeground(amber)
                self.table.setItem(r, c, it)
        self.table.blockSignals(False)
        names = {r.name for r in self._listed}
        self.pair_note.setVisible({"HpaII", "MspI"} <= names)
        self._update_title()
        self._update_note()
        self.digest_button.setEnabled(bool(self.ticked_rows()))
        self.push_overlay()

    def _update_title(self) -> None:
        n = len(self.ticked_rows())
        self.list_title.setText(f"<b>Enzymes</b> &nbsp;<span style='color:gray'>{n} shown, "
                                "site · cuts · in your data</span>")

    def _update_note(self) -> None:
        ticked = self.ticked_rows()
        sites = sum(r.cuts for r in ticked)
        shown = f"{len(ticked)} enzyme{'s' if len(ticked) != 1 else ''} shown, {sites:,} sites"
        meth = self._meth
        if meth is None:
            self._set_note(f"{shown}. Sensitivity per enzyme from REBASE and the Promega / NEB "
                           "tables (CG, CHG, Dam, Dcm). Without data: labels only.", None)
        elif not meth.calls:
            self._set_note(meth.note or "No methylation calls on this sequence.", "amber")
        else:
            blocked = sum(r.blocked for r in ticked)
            self._set_note(f"{shown}; methylation from {Path(meth.path).name} "
                           f"({len(meth.calls):,} calls on this sequence): {blocked:,} "
                           "sites probably not cut.", None)

    def _set_note(self, text: str, tone) -> None:
        colour = "#8a4b00" if tone == "amber" else "gray"
        self.data_note.setText(f"<span style='color:{colour}'>{text}</span>" if text else "")
        self.data_note.setVisible(bool(text))

    def _on_item_changed(self, item: QTableWidgetItem) -> None:
        if item.column() != 0:
            return
        name = item.data(Qt.ItemDataRole.UserRole)
        if name:
            self._ticks[str(name)] = item.checkState() == Qt.CheckState.Checked
            self._update_title()
            self._update_note()
            self.digest_button.setEnabled(bool(self.ticked_rows()))
            self.push_overlay()
            if self.stack.currentWidget() is self.digest_page:
                self.digest_page.schedule()

    def _on_row_selected(self) -> None:
        rows = self.table.selectionModel().selectedRows() if self.table.selectionModel() else []
        if not rows or rows[0].row() >= len(self._listed):
            return
        row = self._listed[rows[0].row()]
        if not row.sites:
            return
        host = self._panel._host
        show = getattr(host, "primers_host_show_range", None) if host is not None else None
        if show is not None:
            site = row.sites[0]
            n = len(self._panel._template)
            try:
                show(max(1, site.start - 60), min(n, site.end + 60))
            except Exception:  # noqa: BLE001
                pass

    # ------------------------------------------------------------------ drawing

    def overlay_sites(self) -> list[tuple]:
        """(start, end, cut, label, blocked, tooltip) of the ticked enzymes' sites."""
        out = []
        for r in self.ticked_rows():
            for s in r.sites:
                out.append((s.start, s.end, s.cut, r.name, s.blocked, s.tooltip(r.site)))
                if len(out) >= ez.MAX_DRAWN_SITES:
                    return out
        return out

    def push_overlay(self) -> None:
        host = self._panel._host
        if host is None or not hasattr(host, "primers_host_set_enzymes"):
            return
        try:
            host.primers_host_set_enzymes(self.overlay_sites() if self.active else [])
        except Exception:  # noqa: BLE001 - drawing is a convenience
            pass

    # ------------------------------------------------------------------ digest

    def show_map_page(self) -> None:
        self.stack.setCurrentWidget(self.map_page)

    def _on_select_all(self) -> None:
        """Tick every enzyme in the list, or untick them all when all are ticked."""
        items = [self.table.item(r, 0) for r in range(self.table.rowCount())]
        items = [i for i in items if i is not None
                 and i.flags() & Qt.ItemFlag.ItemIsUserCheckable]
        if not items:
            return
        every = all(i.checkState() == Qt.CheckState.Checked for i in items)
        state = Qt.CheckState.Unchecked if every else Qt.CheckState.Checked
        for item in items:
            item.setCheckState(state)

    def _on_digest(self) -> None:
        if not self.ticked_rows():
            return
        self.stack.setCurrentWidget(self.digest_page)
        self.digest_page.recompute()
