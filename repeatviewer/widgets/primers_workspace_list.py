"""The list under the map of the Primers and Probes workspace.

One table for what Design found (candidates) and for the primers already on the
sequence. It holds no data of its own: the Design pane's candidate table and the
panel's primer table stay the source (they are kept, not shown), and this list
shows their rows, passes ticks and clicks on to them, and follows their
selection. So every rule of those two tables (what a tick means, what a click
shows on the map, the right-click menu of a primer) still applies.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QHeaderView,
    QTableWidget,
    QTableWidgetItem,
    QWidget,
)

KIND_CANDIDATE = "candidate"
KIND_PRIMER = "primer"
KIND_LABELS = {KIND_CANDIDATE: "candidate", KIND_PRIMER: "on the sequence"}

COL_TICK, COL_NAME, COL_KIND, COL_WHERE, COL_TM = range(5)
HEADERS = ["", "Name", "Kind", "Product or Site", "Tm"]


def _text(table: QTableWidget, row: int, column: int) -> str:
    item = table.item(row, column)
    return item.text() if item is not None else ""


def _state(table: QTableWidget, row: int, column: int = 0) -> Qt.CheckState:
    item = table.item(row, column)
    return item.checkState() if item is not None else Qt.CheckState.Unchecked


def _selected_row(table: QTableWidget) -> int:
    model = table.selectionModel()
    rows = sorted({i.row() for i in model.selectedIndexes()}) if model is not None else []
    return rows[0] if rows else -1


class WorkspaceList(QTableWidget):
    """Candidates and primers on the sequence, as one list."""

    picked = Signal(str)     # KIND_CANDIDATE or KIND_PRIMER: the row the user is on
    changed = Signal()       # rows were rebuilt (counts may differ)

    def __init__(self, candidates: QTableWidget, primers: QTableWidget,
                 parent: QWidget | None = None, *, primer_columns=(0, 1, 2, 3)) -> None:
        super().__init__(0, len(HEADERS), parent)
        self._candidates = candidates
        self._primers = primers
        # tick, name, Tm, site of the primer table
        self._p_tick, self._p_name, self._p_tm, self._p_site = primer_columns
        self._rows: list[tuple[str, int]] = []
        self._busy = False
        self._pending = ""                       # kind selected before its row was here
        self.context_menu_for_primer = None      # callable(pos in the primer table)

        self.setHorizontalHeaderLabels(HEADERS)
        self.verticalHeader().hide()
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.setShowGrid(False)
        self.setAlternatingRowColors(True)
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        head = self.horizontalHeader()
        head.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        head.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        head.setStretchLastSection(True)
        for column, width in ((COL_TICK, 28), (COL_NAME, 150), (COL_KIND, 130), (COL_WHERE, 220)):
            self.setColumnWidth(column, width)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(0)
        self._timer.timeout.connect(self.rebuild)
        for table in (candidates, primers):
            model = table.model()
            model.rowsInserted.connect(self._schedule)
            model.rowsRemoved.connect(self._schedule)
            model.modelReset.connect(self._schedule)
            model.dataChanged.connect(self._schedule)
        candidates.itemSelectionChanged.connect(lambda: self._follow(KIND_CANDIDATE))
        primers.itemSelectionChanged.connect(lambda: self._follow(KIND_PRIMER))
        self.itemChanged.connect(self._on_item_changed)
        self.itemSelectionChanged.connect(self._on_selection)
        self.customContextMenuRequested.connect(self._on_context_menu)
        self.rebuild()

    # ------------------------------------------------------------------ rows

    def _schedule(self, *_args) -> None:
        self._timer.start()

    def _source(self, kind: str) -> QTableWidget:
        return self._candidates if kind == KIND_CANDIDATE else self._primers

    def counts(self) -> tuple[int, int]:
        """(candidates, primers on the sequence)."""
        return (sum(1 for kind, _ in self._rows if kind == KIND_CANDIDATE),
                sum(1 for kind, _ in self._rows if kind == KIND_PRIMER))

    def row_of(self, kind: str, source_row: int) -> int:
        try:
            return self._rows.index((kind, int(source_row)))
        except ValueError:
            return -1

    def kind_of_row(self, row: int) -> str:
        return self._rows[row][0] if 0 <= row < len(self._rows) else ""

    def rebuild(self) -> None:
        """Read both source tables again."""
        self._timer.stop()
        cand, prim = self._candidates, self._primers
        rows: list[tuple[str, int]] = []
        cells: list[tuple[Qt.CheckState, str, str, str, str]] = []
        for r in range(cand.rowCount()):
            rows.append((KIND_CANDIDATE, r))
            cells.append((_state(cand, r), _text(cand, r, 1), KIND_LABELS[KIND_CANDIDATE],
                          _text(cand, r, 2), _text(cand, r, 3)))
        for r in range(prim.rowCount()):
            rows.append((KIND_PRIMER, r))
            cells.append((_state(prim, r, self._p_tick), _text(prim, r, self._p_name),
                          KIND_LABELS[KIND_PRIMER], _text(prim, r, self._p_site),
                          _text(prim, r, self._p_tm)))
        self._busy = True
        try:
            self._rows = rows
            self.setRowCount(len(rows))
            for row, (state, *texts) in enumerate(cells):
                tick = QTableWidgetItem()
                tick.setFlags(Qt.ItemFlag.ItemIsUserCheckable | Qt.ItemFlag.ItemIsEnabled
                              | Qt.ItemFlag.ItemIsSelectable)
                tick.setCheckState(state)
                self.setItem(row, COL_TICK, tick)
                for column, text in enumerate(texts, start=1):
                    self.setItem(row, column, QTableWidgetItem(text))
        finally:
            self._busy = False
        self._show_source_selection()
        self.changed.emit()
        # a source selected a row that was not in the list yet (a primer just
        # added, the best pair of a new run): now it is, so say so
        kind, self._pending = self._pending, ""
        if kind and _selected_row(self._source(kind)) >= 0:
            self.picked.emit(kind)

    # ------------------------------------------------------------- selection

    def _show_source_selection(self) -> None:
        """Select the row the sources have selected (a primer wins over a candidate)."""
        target = -1
        for kind in (KIND_PRIMER, KIND_CANDIDATE):
            source_row = _selected_row(self._source(kind))
            if source_row >= 0:
                target = self.row_of(kind, source_row)
                if target >= 0:
                    break
        self._busy = True
        try:
            if target >= 0:
                self.selectRow(target)
            else:
                self.clearSelection()
        finally:
            self._busy = False

    def _follow(self, kind: str) -> None:
        """A source table changed its selection (a click on the map, a new run)."""
        if self._busy:
            return
        source_row = _selected_row(self._source(kind))
        if source_row < 0:
            if self.kind_of_row(self.currentRow()) == kind and self.selectedIndexes():
                self._show_source_selection()
            return
        row = self.row_of(kind, source_row)
        if row < 0:
            self._pending = kind    # the rows are about to be rebuilt
            self._schedule()
            return
        self._busy = True
        try:
            # one row is "the" row: leaving a primer for a candidate drops the primer
            other = self._source(KIND_CANDIDATE if kind == KIND_PRIMER else KIND_PRIMER)
            if kind == KIND_CANDIDATE and _selected_row(other) >= 0:
                other.clearSelection()
            self.selectRow(row)
        finally:
            self._busy = False
        self.picked.emit(kind)

    def _on_selection(self) -> None:
        """The user clicked a row here."""
        if self._busy:
            return
        rows = sorted({i.row() for i in self.selectedIndexes()})
        if not rows or not (0 <= rows[0] < len(self._rows)):
            return
        kind, source_row = self._rows[rows[0]]
        source = self._source(kind)
        other = self._source(KIND_CANDIDATE if kind == KIND_PRIMER else KIND_PRIMER)
        self._busy = True
        try:
            if _selected_row(other) >= 0:
                other.clearSelection()
            if _selected_row(source) != source_row:
                source.selectRow(source_row)
        finally:
            self._busy = False
        self.picked.emit(kind)

    # ----------------------------------------------------------------- ticks

    def _on_item_changed(self, item: QTableWidgetItem) -> None:
        if self._busy or item.column() != COL_TICK:
            return
        row = item.row()
        if not (0 <= row < len(self._rows)):
            return
        kind, source_row = self._rows[row]
        column = 0 if kind == KIND_CANDIDATE else self._p_tick
        source_item = self._source(kind).item(source_row, column)
        if source_item is not None and source_item.checkState() != item.checkState():
            source_item.setCheckState(item.checkState())

    # ------------------------------------------------------------------ menu

    def _on_context_menu(self, pos) -> None:
        row = self.rowAt(pos.y())
        if self.kind_of_row(row) != KIND_PRIMER or self.context_menu_for_primer is None:
            return
        if row not in {i.row() for i in self.selectedIndexes()}:
            self.selectRow(row)
        # the same point of the screen, in the primer table's own coordinates
        there = self._primers.viewport().mapFromGlobal(self.viewport().mapToGlobal(pos))
        self.context_menu_for_primer(there)
