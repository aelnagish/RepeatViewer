"""Edit > Find in this document: the Find window and its background worker.

The window only collects options and shows results; the main window runs the
search (``repeatviewer.backends.find_in_file``) and moves the Annotation view
from hit to hit. Buttons are Find All / Find Previous / Find Next
(agreed exception to the Cancel / Save / Run rule).
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QObject, QSettings, Qt, Signal, Slot
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import (
    QAbstractItemView,
    QButtonGroup,
    QCheckBox,
    QDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.find_in_file import FindRequest, FoundItem, find_in_fasta
from repeatviewer.backends.motif_search import FindCancelled, FindOptions
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.spin_boxes import NoWheelSpinBox as QSpinBox

_SETTINGS = "find"

QUERY_HELP = (
    "Type a <b>DNA motif</b> (IUPAC codes allowed), a <b>peptide</b> (with Protein), "
    "or part of a <b>sequence name</b> or <b>annotation name</b>.<br>"
    "Spaces, line breaks and numbers in pasted sequence are ignored.<br>"
    "<b>Numbers only</b> (114590001, 1000-2000, scaffold_5:1000-2000) go to that position."
)
STRAND_HELP = (
    "<b>Include reverse complement</b>: report hits on the minus strand too "
    "(shown as strand -).<br>A palindrome such as GAATTC is reported once."
)
IUPAC_HELP = (
    "<b>On</b>: IUPAC ambiguity codes are treated as wildcards.<br><br>"
    "<b>Query</b> — e.g. R matches A or G, Y matches C or T, N matches any base.<br>"
    "<b>Sequence</b> — ambiguous letters in the sequence (N, R, Y, …) also count "
    "as a match for the corresponding motif base.<br><br>"
    "<b>Off</b>: each letter matches only itself (an R in the motif only matches "
    "an R in the sequence)."
)
GAPS_HELP = (
    "Ignore gap characters (- and .) in the sequence, so a motif is found even when an "
    "alignment gap splits it.<br>Positions are those of the gapped sequence."
)
MISMATCH_HELP = (
    "How many positions may differ from the motif (substitutions only, no indels).<br>"
    "Useful for FISH probes, primers and diverged satellite monomers."
)
PROTEIN_HELP = (
    "<b>Nucleotides</b>: search the DNA.<br>"
    "<b>Protein</b>: the motif is a peptide, searched in all six reading frames "
    "(standard genetic code). X in the motif matches any amino acid."
)
NAMES_HELP = (
    "<b>Sequence names</b>: sequences whose name contains the text.<br>"
    "<b>Annotation names</b>: annotations whose name or type contains the text.<br>"
    "<b>Case sensitive</b> applies to these two only."
)


class FindWorker(QObject):
    """Runs find_in_fasta off the UI thread (progress per sequence, Stop)."""

    total = Signal(int)  # number of sequences, once counted
    progress = Signal(int, str)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, path: Path, request: FindRequest, features, cancel=None):
        super().__init__()
        self._path = Path(path)
        self._request = request
        self._features = features
        self._cancel = cancel
        self._step = 1

    @Slot()
    def run(self) -> None:
        from repeatviewer.backends.job_cancel import CANCELLED_MESSAGE

        try:
            from repeatviewer.browser_file_metadata import scan_fasta_headers

            count = int(scan_fasta_headers(self._path)[0])
        except Exception:  # noqa: BLE001 - progress is optional
            count = 0
        self.total.emit(count)
        # About 200 updates at most (rule 15: big files stay smooth).
        self._step = max(1, count // 200)

        def on_progress(done: int, name: str) -> None:
            if done % self._step == 0 or done == count:
                self.progress.emit(done, name)

        try:
            result = find_in_fasta(
                self._path,
                self._request,
                features_by_id=self._features,
                cancelled=self._cancel,
                progress=on_progress,
            )
        except FindCancelled:
            self.failed.emit(CANCELLED_MESSAGE)
            return
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(str(exc))
            return
        self.finished.emit(result)


def _check(text: str, parent: QWidget) -> QCheckBox:
    return QCheckBox(text, parent)


def _with_help(widget: QWidget, title: str, text: str, parent: QWidget) -> QWidget:
    wrap = QWidget(parent)
    row = QHBoxLayout(wrap)
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(6)
    row.addWidget(widget)
    row.addWidget(make_help_button(wrap, title, text))
    row.addStretch(1)
    return wrap


class FindDialog(QDialog):
    """Non-modal Find window."""

    findRequested = Signal(int)  # +1 Find Next, -1 Find Previous, 0 Find All
    itemActivated = Signal(int)  # a row of the results list was clicked
    saveRequested = Signal()
    closed = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Find in this document")
        self.setModal(False)
        self.setWindowModality(Qt.WindowModality.NonModal)
        self._items: list[FoundItem] = []
        self._filling = False

        self._query = QLineEdit(self)
        self._query.setClearButtonEnabled(True)
        self._query.setPlaceholderText("Motif, sequence name, annotation name or position")
        self._query.setMinimumWidth(320)
        self._query.returnPressed.connect(lambda: self.findRequested.emit(1))

        # Labels / hierarchy match the Annotation Preview Find tab.
        self._in_sequence = _check("Sequence", self)
        self._nucleotides = QRadioButton("Nucleotides", self)
        self._protein = QRadioButton("Protein (6 frames)", self)
        self._kind_group = QButtonGroup(self)
        self._kind_group.addButton(self._nucleotides)
        self._kind_group.addButton(self._protein)
        self._both_strands = _check("Include reverse complement", self)
        self._iupac = _check("Interpret IUPAC codes", self)
        self._allow_gaps = _check("Allow gaps in sequence", self)
        self._mismatches = QSpinBox(self)
        self._mismatches.setRange(0, 10)
        self._in_names = _check("Sequence names", self)
        self._in_annotations = _check("Annotation names", self)
        self._case = _check("Case sensitive", self)

        form = QGridLayout()
        form.setHorizontalSpacing(10)
        form.setVerticalSpacing(6)
        form.addWidget(QLabel("Find:", self), 0, 0)
        form.addWidget(
            _with_help(self._query, "Find in this document", QUERY_HELP, self),
            0,
            1,
        )
        form.addWidget(
            QLabel("Search in:", self),
            1,
            0,
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop,
        )

        search_col = QWidget(self)
        search_layout = QVBoxLayout(search_col)
        search_layout.setContentsMargins(0, 0, 0, 0)
        search_layout.setSpacing(4)
        search_layout.addWidget(self._in_sequence)

        # Indented under Sequence (keep ? help on each option).
        seq_opts = QWidget(search_col)
        seq_layout = QVBoxLayout(seq_opts)
        seq_layout.setContentsMargins(22, 0, 0, 0)
        seq_layout.setSpacing(4)
        kind_row = QWidget(seq_opts)
        kind_layout = QHBoxLayout(kind_row)
        kind_layout.setContentsMargins(0, 0, 0, 0)
        kind_layout.setSpacing(6)
        kind_layout.addWidget(self._nucleotides)
        kind_layout.addWidget(
            make_help_button(kind_row, "Nucleotides or Protein", PROTEIN_HELP)
        )
        kind_layout.addStretch(1)
        seq_layout.addWidget(kind_row)
        seq_layout.addWidget(self._protein)
        seq_layout.addWidget(
            _with_help(
                self._both_strands, "Include reverse complement", STRAND_HELP, self
            )
        )
        seq_layout.addWidget(
            _with_help(self._iupac, "Interpret IUPAC codes", IUPAC_HELP, self)
        )
        seq_layout.addWidget(
            _with_help(self._allow_gaps, "Allow gaps in sequence", GAPS_HELP, self)
        )
        mm_row = QWidget(seq_opts)
        mm_layout = QHBoxLayout(mm_row)
        mm_layout.setContentsMargins(0, 0, 0, 0)
        mm_layout.addWidget(QLabel("Max mismatches:", mm_row))
        mm_layout.addWidget(self._mismatches)
        mm_layout.addWidget(make_help_button(mm_row, "Max mismatches", MISMATCH_HELP))
        mm_layout.addStretch(1)
        seq_layout.addWidget(mm_row)
        search_layout.addWidget(seq_opts)

        # Same level as Sequence.
        search_layout.addWidget(self._in_names)
        search_layout.addWidget(self._in_annotations)
        search_layout.addWidget(
            _with_help(self._case, "Names", NAMES_HELP, self)
        )
        form.addWidget(search_col, 1, 1)

        self._status = QLabel("", self)
        self._status.setWordWrap(True)
        self._status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        self._table = QTableWidget(0, 5, self)
        self._table.setHorizontalHeaderLabels(["Sequence", "Start", "End", "Strand", "Found In"])
        self._table.verticalHeader().setVisible(False)
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        header = self._table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        for col in range(1, 5):
            header.setSectionResizeMode(col, QHeaderView.ResizeMode.ResizeToContents)
        self._table.setMinimumHeight(140)
        self._table.hide()
        self._table.itemSelectionChanged.connect(self._on_table_selection)

        results_line = QFrame(self)
        results_line.setFrameShape(QFrame.Shape.HLine)
        results_line.setFrameShadow(QFrame.Shadow.Sunken)

        self._save_btn = QPushButton("Save Hits as Annotations…", self)
        self._save_btn.setEnabled(False)
        self._save_btn.clicked.connect(self.saveRequested.emit)
        self._all_btn = QPushButton("Find All", self)
        self._all_btn.clicked.connect(lambda: self.findRequested.emit(0))
        self._prev_btn = QPushButton("Find Previous", self)
        self._prev_btn.clicked.connect(lambda: self.findRequested.emit(-1))
        self._next_btn = QPushButton("Find Next", self)
        self._next_btn.setDefault(True)
        self._next_btn.clicked.connect(lambda: self.findRequested.emit(1))
        buttons = QHBoxLayout()
        buttons.addWidget(self._save_btn)
        buttons.addStretch(1)
        for btn in (self._all_btn, self._prev_btn, self._next_btn):
            buttons.addWidget(btn)

        root = QVBoxLayout(self)
        root.addLayout(form)
        root.addWidget(results_line)
        root.addWidget(self._status)
        root.addWidget(self._table, stretch=1)
        root.addLayout(buttons)

        # Same keys as the Edit menu while this window is in front.
        for keys, direction in (((Qt.Key.Key_F3, "Ctrl+G"), 1), ((QKeySequence("Shift+F3"),), -1)):
            for key in keys:
                shortcut = QShortcut(QKeySequence(key), self)
                shortcut.activated.connect(lambda d=direction: self.findRequested.emit(d))

        self._in_sequence.toggled.connect(self._sync_enabled)
        self._protein.toggled.connect(self._sync_enabled)
        self._load_settings()
        self._sync_enabled()

    # ------------------------------------------------------------------ options

    def _load_settings(self) -> None:
        s = QSettings()
        s.beginGroup(_SETTINGS)

        def flag(key: str, default: bool) -> bool:
            return str(s.value(key, default)).lower() in {"1", "true", "yes"}

        self._query.setText(str(s.value("query", "") or ""))
        self._in_sequence.setChecked(flag("in_sequence", True))
        (self._protein if flag("protein", False) else self._nucleotides).setChecked(True)
        self._both_strands.setChecked(flag("both_strands", True))
        # One IUPAC control (legacy split query/sequence settings → OR).
        if s.contains("iupac"):
            iupac_on = flag("iupac", True)
        else:
            iupac_on = flag("iupac_query", True) or flag("iupac_sequence", False)
        self._iupac.setChecked(iupac_on)
        self._allow_gaps.setChecked(flag("allow_gaps", False))
        try:
            self._mismatches.setValue(int(s.value("mismatches", 0)))
        except (TypeError, ValueError):
            self._mismatches.setValue(0)
        self._in_names.setChecked(flag("in_names", True))
        self._in_annotations.setChecked(flag("in_annotations", True))
        self._case.setChecked(flag("case_sensitive", False))
        s.endGroup()

    def save_settings(self) -> None:
        s = QSettings()
        s.beginGroup(_SETTINGS)
        s.setValue("query", self._query.text())
        s.setValue("in_sequence", self._in_sequence.isChecked())
        s.setValue("protein", self._protein.isChecked())
        s.setValue("both_strands", self._both_strands.isChecked())
        s.setValue("iupac", self._iupac.isChecked())
        s.setValue("iupac_query", self._iupac.isChecked())
        s.setValue("iupac_sequence", self._iupac.isChecked())
        s.setValue("allow_gaps", self._allow_gaps.isChecked())
        s.setValue("mismatches", self._mismatches.value())
        s.setValue("in_names", self._in_names.isChecked())
        s.setValue("in_annotations", self._in_annotations.isChecked())
        s.setValue("case_sensitive", self._case.isChecked())
        s.endGroup()

    def _sync_enabled(self) -> None:
        on = self._in_sequence.isChecked()
        for widget in (
            self._nucleotides,
            self._protein,
            self._both_strands,
            self._iupac,
            self._allow_gaps,
            self._mismatches,
        ):
            widget.setEnabled(on)

    def query_text(self) -> str:
        return self._query.text().strip()

    def request(self) -> FindRequest:
        iupac = self._iupac.isChecked()
        options = FindOptions(
            protein=self._protein.isChecked(),
            both_strands=self._both_strands.isChecked(),
            iupac_in_query=iupac,
            iupac_in_sequence=iupac,
            allow_gaps=self._allow_gaps.isChecked(),
            max_mismatches=self._mismatches.value(),
        )
        return FindRequest(
            query=self.query_text(),
            options=options,
            in_sequence=self._in_sequence.isChecked(),
            in_names=self._in_names.isChecked(),
            in_annotations=self._in_annotations.isChecked(),
            case_sensitive=self._case.isChecked(),
        )

    def focus_query(self) -> None:
        self._query.setFocus(Qt.FocusReason.OtherFocusReason)
        self._query.selectAll()

    # ------------------------------------------------------------------ results

    def set_busy(self, busy: bool) -> None:
        for btn in (self._next_btn, self._prev_btn, self._all_btn):
            btn.setEnabled(not busy)
        self._save_btn.setEnabled(not busy and any(i.kind == "sequence" for i in self._items))

    def set_status(self, text: str) -> None:
        self._status.setText(text)

    def set_results(self, items: list[FoundItem]) -> None:
        self._items = list(items)
        self._filling = True
        try:
            shown = self._items[:2000]  # the list is for browsing; navigation covers all
            self._table.setRowCount(len(shown))
            kinds = {"sequence": "Sequence", "name": "Sequence name", "annotation": "Annotation"}
            for row, item in enumerate(shown):
                found_in = kinds.get(item.kind, item.kind)
                if item.kind == "annotation" and item.label:
                    found_in = f"Annotation: {item.label}"
                elif item.kind == "sequence" and item.mismatches:
                    found_in = f"Sequence, {item.mismatches} mismatch" + ("es" if item.mismatches > 1 else "")
                values = [item.seq_id, f"{item.start:,}", f"{item.end:,}", item.strand, found_in]
                for col, value in enumerate(values):
                    cell = QTableWidgetItem(value)
                    if col in (1, 2):
                        cell.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                    self._table.setItem(row, col, cell)
        finally:
            self._filling = False
        self._table.setVisible(bool(self._items))
        self._save_btn.setEnabled(any(i.kind == "sequence" for i in self._items))

    def set_current(self, index: int) -> None:
        self._filling = True
        try:
            if 0 <= index < self._table.rowCount():
                self._table.selectRow(index)
                self._table.scrollToItem(self._table.item(index, 0))
            else:
                self._table.clearSelection()
        finally:
            self._filling = False

    def _on_table_selection(self) -> None:
        if self._filling:
            return
        rows = self._table.selectionModel().selectedRows()
        if rows:
            self.itemActivated.emit(rows[0].row())

    def closeEvent(self, event) -> None:  # noqa: N802
        self.save_settings()
        self.closed.emit()
        super().closeEvent(event)
