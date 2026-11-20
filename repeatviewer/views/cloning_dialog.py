"""Cloning: put a PCR product into a vector and see what comes out.

Pick the insert, the method (TA, blunt, restriction) and the vector; the window
shows both orientations, the finished plasmid, the colony PCR sizes and what to
watch out for, including a tandem array that will not survive an ordinary
strain.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import html
import json
from pathlib import Path

from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import cloning as cl
from repeatviewer.backends.primers import pcr
from repeatviewer.backends.primers import tails as rtail
from repeatviewer.views.sanger_dialog import _read_fasta
from repeatviewer.widgets.help_button import make_help_button

_SETTINGS_KEY = "primers/cloning"
READ_LENGTH = 700

CLONING_HELP = (
    "<b>Insert</b>: the PCR product, with the primer tails already on it. Use the open "
    "file, or Browse.<br>"
    "<b>Method</b>: <i>TA</i> needs the 3′ A that Taq leaves; <i>Blunt</i> needs a "
    "proofreading product; <i>Restriction</i> uses the sites in the primer tails.<br>"
    "<b>Vector</b>: pGEM-T Easy, pJET1.2/blunt or your own.<br>"
    "<b>Empty vector product</b>: what your colony-PCR pair gives on the vector with no "
    "insert. Your vector's map or manual says it. With it, the window gives the size a "
    "positive colony shows, and how far a sequencing read gets into the insert.<br>"
    "No vector sequence ships with RepeatViewer, so every size here comes from the "
    "numbers you give, never from a guessed map."
)
RESULT_HELP = (
    "Both orientations of the insert. A TA or blunt ligation gives both; two different "
    "enzymes give only the first. <b>Colony PCR</b> is the empty-vector product plus the "
    "insert, so a colony with the insert runs that much higher than one without."
)


def _row(label: str, *widgets, parent: QWidget) -> QHBoxLayout:
    row = QHBoxLayout()
    lab = QLabel(label, parent)
    lab.setFixedWidth(150)
    lab.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    row.addWidget(lab)
    for w in widgets:
        row.addWidget(w)
    row.addStretch(1)
    return row


class CloningDialog(QDialog):
    """``insert`` is the product to clone; everything else is chosen here."""

    def __init__(self, insert: str = "", *, insert_name: str = "", folder: Path | None = None,
                 features=(), parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Cloning")
        self.setObjectName("CloningDialog")
        self._insert = str(insert or "").upper()
        self._insert_name = insert_name or "insert"
        self._folder = Path(folder) if folder else Path.cwd()
        self._features = list(features or ())
        self._result = cl.CloneResult()
        saved = self._load()

        lay = QVBoxLayout(self)
        lay.setSpacing(8)
        head = QHBoxLayout()
        self.insert_label = QLabel(self)
        self.insert_label.setWordWrap(True)
        head.addWidget(self.insert_label, 1)
        self.insert_button = QPushButton("Insert Sequence", self)
        head.addWidget(self.insert_button)
        head.addWidget(make_help_button(self, "Cloning", CLONING_HELP))
        lay.addLayout(head)

        self.record_combo = QComboBox(self)
        self.record_combo.hide()
        lay.addWidget(self.record_combo)

        self.vector = QComboBox(self)
        for v in cl.PRESETS:
            self.vector.addItem(v.name, v.name)
        self.vector.setCurrentIndex(max(0, self.vector.findData(saved.get("vector"))))
        lay.addLayout(_row("Vector", self.vector, parent=self))

        self.method = QComboBox(self)
        lay.addLayout(_row("Method", self.method, parent=self))

        self.polymerase = QComboBox(self)
        for name in pcr.POLYMERASES:
            self.polymerase.addItem(name, name)
        self.polymerase.setCurrentIndex(max(0, self.polymerase.findData(saved.get("polymerase"))))
        lay.addLayout(_row("Polymerase", self.polymerase, parent=self))

        self.left = QComboBox(self)
        self.right = QComboBox(self)
        for box in (self.left, self.right):
            box.addItem("none", "")
            for name in rtail.TAIL_ENZYMES:
                box.addItem(name, name)
        self.enzyme_row = _row("Enzymes", self.left, self.right, parent=self)
        lay.addLayout(self.enzyme_row)

        self.vector_length = QSpinBox(self)
        self.vector_length.setRange(0, 1_000_000)
        self.vector_length.setSuffix(" bp")
        self.vector_length.setGroupSeparatorShown(True)
        self.empty_product = QSpinBox(self)
        self.empty_product.setRange(0, 100_000)
        self.empty_product.setSuffix(" bp")
        self.empty_product.setGroupSeparatorShown(True)
        self.empty_product.setValue(int(saved.get("empty_product", 0)))
        lay.addLayout(_row("Vector length", self.vector_length, parent=self))
        lay.addLayout(_row("Empty vector product", self.empty_product, parent=self))

        self.verdict = QLabel(self)
        self.verdict.setWordWrap(True)
        self.verdict.setObjectName("cloningVerdict")
        lay.addWidget(self.verdict)

        table_head = QHBoxLayout()
        table_head.addWidget(QLabel("<b>Clone</b>", self), 1)
        table_head.addWidget(make_help_button(self, "Clone", RESULT_HELP))
        lay.addLayout(table_head)
        self.table = QTableWidget(0, 4, self)
        self.table.setHorizontalHeaderLabels(["Orientation", "Insert", "Plasmid", "Colony PCR"])
        self.table.verticalHeader().hide()
        self.table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self.table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self.table.setShowGrid(False)
        self.table.setAlternatingRowColors(True)
        self.table.horizontalHeader().setDefaultAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self.table.horizontalHeader().setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        self.table.setMaximumHeight(110)
        lay.addWidget(self.table)

        self.notes = QLabel(self)
        self.notes.setWordWrap(True)
        self.notes.setTextFormat(Qt.TextFormat.RichText)
        lay.addWidget(self.notes, 1)

        buttons = QHBoxLayout()
        self.gel_button = QPushButton("Colony PCR Gel", self)
        self.report_button = QPushButton("Save Report", self)
        close = QPushButton("Close", self)
        buttons.addStretch(1)
        for b in (self.gel_button, self.report_button, close):
            buttons.addWidget(b)
        lay.addLayout(buttons)

        for box in (self.vector, self.method, self.polymerase):
            box.setMinimumWidth(220)
        self.insert_button.clicked.connect(self._on_pick_insert)
        self.record_combo.currentIndexChanged.connect(self._on_record)
        self.vector.currentIndexChanged.connect(self._on_vector)
        for box in (self.method, self.polymerase, self.left, self.right):
            box.currentIndexChanged.connect(self.recompute)
        for spin in (self.vector_length, self.empty_product):
            spin.valueChanged.connect(self.recompute)
        self.gel_button.clicked.connect(self._on_gel)
        self.report_button.clicked.connect(self._on_report)
        close.clicked.connect(self.accept)
        self.resize(760, 640)
        self._on_vector()

    # ------------------------------------------------------------------ settings

    def _load(self) -> dict:
        raw = QSettings().value(_SETTINGS_KEY, "")
        try:
            d = json.loads(raw) if raw else {}
        except (TypeError, ValueError):
            d = {}
        return d if isinstance(d, dict) else {}

    def _save(self) -> None:
        QSettings().setValue(_SETTINGS_KEY, json.dumps(
            {"vector": self.vector.currentData(), "method": self.method.currentData(),
             "polymerase": self.polymerase.currentData(),
             "empty_product": int(self.empty_product.value())}))

    # ------------------------------------------------------------------ inputs

    def vector_preset(self) -> cl.Vector:
        return cl.preset(str(self.vector.currentData()))

    def _on_vector(self, *_a) -> None:
        v = self.vector_preset()
        saved = self._load().get("method")
        self.method.blockSignals(True)
        self.method.clear()
        for m in v.methods:
            self.method.addItem(cl.METHOD_LABELS[m], m)
        i = self.method.findData(saved)
        self.method.setCurrentIndex(i if i >= 0 else 0)
        self.method.blockSignals(False)
        own = v.length == 0
        self.vector_length.setEnabled(own)
        if not own:
            self.vector_length.setValue(v.length)
        self.recompute()

    def set_insert(self, sequence: str, name: str = "", features=()) -> None:
        self._insert = str(sequence or "").upper()
        if name:
            self._insert_name = name
        self._features = list(features or ())
        self.recompute()

    def _on_pick_insert(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Insert Sequence", str(self._folder),
            "Sequences (*.fa *.fasta *.fas *.fna *.txt *.seq);;All files (*)")
        if not path:
            return
        records = _read_fasta(Path(path))
        if not records:
            self.verdict.setText("That file holds no sequence.")
            return
        self._folder = Path(path).parent
        self.record_combo.blockSignals(True)
        self.record_combo.clear()
        for name, seq in records:
            self.record_combo.addItem(f"{name} ({len(seq):,} bp)", (name, seq))
        self.record_combo.setVisible(len(records) > 1)
        self.record_combo.setCurrentIndex(0)
        self.record_combo.blockSignals(False)
        self.set_insert(records[0][1], records[0][0])

    def _on_record(self, index: int) -> None:
        data = self.record_combo.itemData(index)
        if data:
            self.set_insert(data[1], data[0])

    # ------------------------------------------------------------------ result

    def recompute(self, *_a) -> None:
        v = self.vector_preset()
        method = str(self.method.currentData() or cl.METHOD_TA)
        restriction = method == cl.METHOD_RE
        for i in range(self.enzyme_row.count()):
            w = self.enzyme_row.itemAt(i).widget()
            if w is not None:
                w.setVisible(restriction)
        if v.length == 0:
            v = cl.Vector(v.name, int(self.vector_length.value()), v.methods, v.primers, v.note)
        self._result = cl.plan(
            self._insert, v, method,
            enzymes=(str(self.left.currentData() or ""), str(self.right.currentData() or "")),
            polymerase=str(self.polymerase.currentData() or ""),
            empty_product=int(self.empty_product.value()), features=self._features)
        self._save()
        self._show(v, method)

    def _show(self, vector: cl.Vector, method: str) -> None:
        n = len(self._insert)
        self.insert_label.setText(
            f"<b>{html.escape(self._insert_name)}</b>, {n:,} bp" if n
            else "No insert yet: press Insert Sequence.")
        tone, line = cl.summary(self._result, vector, method)
        bg, border, fg = (("#E7F6EC", "#9BD3AE", "#14532D") if tone == "green"
                          else ("#FFF4E5", "#E8A33D", "#5A3700"))
        self.verdict.setStyleSheet(
            f"#cloningVerdict {{ background: {bg}; border: 1px solid {border};"
            f" border-radius: 6px; padding: 6px; color: {fg}; }}")
        self.verdict.setText(line)

        clones = self._result.clones
        if method == cl.METHOD_RE:
            left, right = self.left.currentData(), self.right.currentData()
            if left and right and left != right:
                clones = clones[:1]          # two different ends: one way in only
        self.table.setRowCount(len(clones))
        for r, c in enumerate(clones):
            cells = [c.label(), f"{c.insert_length:,} bp",
                     f"{c.total_length:,} bp" if c.total_length else "give the vector length",
                     f"{c.colony_product:,} bp" if c.colony_product else "give the empty product"]
            for col, text in enumerate(cells):
                self.table.setItem(r, col, QTableWidgetItem(text))
        head = self.table.horizontalHeader().height()
        rows = sum(self.table.rowHeight(r) for r in range(self.table.rowCount()))
        self.table.setFixedHeight(head + rows + 4)

        # the first warning is already in the box above
        bits = [f"<span style='color:#8a4b00'>{html.escape(w)}.</span>"
                for w in self._result.warnings[1:]]
        bits += [html.escape(t) + "." for t in self._result.notes]
        if vector.note:
            bits.append(html.escape(vector.note))
        if self._result.empty_product and self._insert:
            reach = cl.reads_into_insert(self._result.empty_product, READ_LENGTH)
            names = ", ".join(p.name for p in vector.primers)
            covered = ("the whole insert" if reach >= len(self._insert)
                       else f"about {reach:,} bp of it from each side")
            bits.append(f"Sequencing primers: {html.escape(names)}. A {READ_LENGTH} base read "
                        f"covers {covered}.")
        self.notes.setText("<br>".join(bits))
        has = bool(self._result.clones and self._result.empty_product)
        self.gel_button.setEnabled(has)
        self.report_button.setEnabled(bool(self._insert))

    # ------------------------------------------------------------------ actions

    def _safe_name(self) -> str:
        return "".join(c if c.isalnum() or c in "-_." else "_" for c in self._insert_name)[:40]

    def _free_path(self, suffix: str) -> Path:
        base = self._safe_name() or "clone"
        out = self._folder / f"{base}{suffix}"
        n = 2
        while out.exists():
            out = self._folder / f"{base}_{n}{suffix}"
            n += 1
        return out

    def _show_output(self, path: Path) -> None:
        main = self.parent().window() if self.parent() else None
        show = getattr(main, "_show_new_output", None)
        if show is not None:
            show(path)

    def _on_gel(self) -> None:
        doc = cl.colony_gel(self._result, title=f"{self._safe_name()} colony PCR")
        if doc is None:
            return
        from repeatviewer.backends.primers import gel_doc as gd

        out = self._free_path(gd.GEL_SUFFIX)
        try:
            doc.save(out)
        except OSError as exc:
            self.verdict.setText(f"Could not write the gel: {exc}")
            return
        self.gel_path = out
        self._show_output(out)

    def _on_report(self) -> None:
        out = self._free_path("_cloning.html")
        try:
            cl.write_report(out, self._insert_name, self._result, self.vector_preset(),
                            str(self.method.currentData()))
        except OSError as exc:
            self.verdict.setText(f"Could not write the report: {exc}")
            return
        self.report_path = out
        self._show_output(out)
