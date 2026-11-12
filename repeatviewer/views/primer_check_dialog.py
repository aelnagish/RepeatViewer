"""Primer Database › Check a Sequence: test the selected primers on one sequence.

A small window for one job. Everything Map Primers asks that does not matter here
(which inputs, which primer mode, where primers come from) is left out; "All
Settings…" opens the full Map Primers window filled in the same way.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QFormLayout,
    QFrame,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers.oligos import KIND_LABELS, Primer, PrimerDatabase
from repeatviewer.backends.primers.thermo import ThermoSettings, analyze_oligo
from repeatviewer.widgets.path_browse import PathBrowseRow
from repeatviewer.widgets.spin_boxes import NoWheelSpinBox as QSpinBox

_FASTA_FILTER = "FASTA (*.fa *.fasta *.fna *.fas *.ffn *.fa.gz *.fasta.gz);;All files (*)"
_FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas", ".ffn", ".fa.gz", ".fasta.gz")
_SETTINGS = "primers/check_window"

NO_SEQUENCE = ("Choose the sequence to test the primers on: click a FASTA in the Files "
               "panel, or use Browse.")


def _is_fasta(path: Path) -> bool:
    name = path.name.lower()
    return path.is_file() and name.endswith(_FASTA_SUFFIXES)


def _card(parent: QWidget, title: str) -> tuple[QFrame, QVBoxLayout]:
    card = QFrame(parent)
    card.setObjectName("CheckCard")
    card.setFrameShape(QFrame.Shape.StyledPanel)
    card.setStyleSheet("QFrame#CheckCard { border: 1px solid palette(mid); "
                       "border-radius: 8px; }")
    lay = QVBoxLayout(card)
    lay.setContentsMargins(16, 12, 16, 14)
    lay.setSpacing(8)
    lay.addWidget(QLabel(f"<b>{title}</b>", card))
    return card, lay


class CheckPrimersDialog(QDialog):
    """Which primers, which sequence, how strict. Run."""

    runRequested = Signal(dict)          # {"sequence", "product_min", ...}
    allSettingsRequested = Signal(dict)

    def __init__(self, db: PrimerDatabase | None, names: list[str], *,
                 thermo: ThermoSettings | None = None, sequence: Path | None = None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Check Primers on a Sequence")
        self._names = list(names)
        primers = list(db.primers) if db is not None else []
        if self._names:
            wanted = set(self._names)
            primers = [p for p in primers if p.name in wanted]
        self._primers = primers

        # ---- what is tested
        tested, tl = _card(self, "Primers")
        if self._names:
            self._table = QTableWidget(len(primers), 5, tested)
            self._table.setHorizontalHeaderLabels(
                ["Name", "Direction", "Sequence 5′→3′", "Length", "Tm Primer3 °C"])
            self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
            self._table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
            self._table.setFocusPolicy(Qt.FocusPolicy.NoFocus)
            self._table.verticalHeader().setVisible(False)
            self._table.verticalHeader().setDefaultSectionSize(26)
            self._table.setShowGrid(False)
            self._table.setAlternatingRowColors(True)
            hh = self._table.horizontalHeader()
            hh.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
            hh.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
            hh.setStretchLastSection(True)
            for i, p in enumerate(primers):
                self._fill_row(i, p, thermo)
            self._table.resizeColumnsToContents()
            for c in range(self._table.columnCount()):
                self._table.setColumnWidth(c, self._table.columnWidth(c) + 12)
            shown = min(6, max(1, len(primers)))
            frame = 2 * self._table.frameWidth()
            self._table.setFixedHeight(shown * 26 + hh.sizeHint().height() + frame + 2)
            if len(primers) <= shown:
                self._table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            self._table_width = sum(self._table.columnWidth(c)
                                    for c in range(self._table.columnCount())) + frame
            tl.addWidget(self._table)
        else:
            self._table = None
            n = len(primers)
            tl.addWidget(QLabel(
                f"All {n} saved primer{'s' if n != 1 else ''} in the Primer Database. "
                "Select primers in the Database first to test only those.", tested))

        # ---- the sequence
        seq_card, sl = _card(self, "Sequence")
        self._seq = PathBrowseRow(seq_card, placeholder="Click a FASTA in the Files panel, "
                                  "or Browse", browse_title="Choose the sequence FASTA",
                                  mode="file", file_filter=_FASTA_FILTER)
        self._seq.set_default_tooltip("Click a FASTA in the Files panel to fill this, or use "
                                      "Browse.")
        if sequence is not None:
            self._seq.set_path(sequence, emit=False)
        sl.addWidget(self._seq)

        # ---- how strict
        st = QSettings()
        strict, kl = _card(self, "Products and mismatches")
        self._pmin = QSpinBox(strict)
        self._pmax = QSpinBox(strict)
        for box in (self._pmin, self._pmax):
            box.setRange(1, 1_000_000)
            box.setSuffix(" bp")
            box.setMinimumWidth(110)
        self._pmin.setValue(int(st.value(f"{_SETTINGS}/product_min", 40)))
        self._pmax.setValue(int(st.value(f"{_SETTINGS}/product_max", 5000)))
        prod = QWidget(strict)
        pl = QHBoxLayout(prod)
        pl.setContentsMargins(0, 0, 0, 0)
        pl.addWidget(self._pmin)
        pl.addWidget(QLabel("to", prod))
        pl.addWidget(self._pmax)
        pl.addStretch(1)
        self._mm = QSpinBox(strict)
        self._mm.setRange(0, 10)
        self._mm.setValue(int(st.value(f"{_SETTINGS}/max_mismatches", 2)))
        self._end = QSpinBox(strict)
        self._end.setRange(0, 15)
        self._end.setValue(int(st.value(f"{_SETTINGS}/three_prime_exact", 3)))
        mm = QWidget(strict)
        ml = QHBoxLayout(mm)
        ml.setContentsMargins(0, 0, 0, 0)
        ml.addWidget(self._mm)
        ml.addWidget(QLabel("last", mm))
        ml.addWidget(self._end)
        ml.addWidget(QLabel("3′ bases must match", mm))
        ml.addStretch(1)
        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.addRow("Product size:", prod)
        form.addRow("Maximum mismatches:", mm)
        kl.addLayout(form)
        strict.setToolTip("Mismatches: a primer from another species may not match "
                          "perfectly. The last 3′ bases must match for the polymerase to "
                          "start, so a product is only counted when they do.")

        # ---- buttons
        more = QPushButton("All Settings…", self)
        more.setToolTip("Open the full Map Primers window with these primers")
        more.clicked.connect(self._all_settings)
        cancel = QPushButton("Cancel", self)
        cancel.clicked.connect(self.reject)
        self._run_btn = QPushButton("Run", self)
        self._run_btn.setDefault(True)
        self._run_btn.clicked.connect(self._run)
        bar = QHBoxLayout()
        bar.addWidget(more)
        bar.addStretch(1)
        bar.addWidget(cancel)
        bar.addWidget(self._run_btn)

        lay = QVBoxLayout(self)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.setSpacing(12)
        lay.addWidget(tested)
        lay.addWidget(seq_card)
        lay.addWidget(strict)
        lay.addStretch(1)
        lay.addLayout(bar)
        # Wide enough to show the whole table (primers are rarely over 40 bases).
        width = max(640, min(1000, getattr(self, "_table_width", 0) + 2 * 16 + 2 * 16 + 4))
        self.resize(width, self.sizeHint().height())

    def _fill_row(self, i: int, p: Primer, thermo: ThermoSettings | None) -> None:
        from repeatviewer.backends.primers.inventory import kind_of

        try:
            tm = f"{analyze_oligo(p.sequence, thermo).tm:.1f}"
        except Exception:  # noqa: BLE001 - a primer Primer3 cannot read
            tm = ""
        values = (p.name, KIND_LABELS.get(kind_of(p.name, p.kind), "Oligo"), p.sequence,
                  str(len(p.sequence)), tm)
        for j, v in enumerate(values):
            item = QTableWidgetItem(v)
            if j == 2:
                item.setFont(self._mono())
            if j >= 3:
                item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            self._table.setItem(i, j, item)

    @staticmethod
    def _mono():
        from PySide6.QtGui import QFontDatabase

        return QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)

    # ------------------------------------------------------------ inputs

    def apply_path_from_panel(self, path: Path) -> None:
        """A click in the Files panel: a FASTA becomes the sequence."""
        path = Path(path)
        if _is_fasta(path):
            self._seq.set_path(path)
            QTimer.singleShot(0, self.raise_)

    def sequence(self) -> Path | None:
        return self._seq.path()

    def values(self) -> dict:
        lo, hi = sorted((self._pmin.value(), self._pmax.value()))
        seq = self.sequence()
        return {"sequence": str(seq) if seq else "", "primer_names": list(self._names),
                "product_min": lo, "product_max": hi,
                "max_mismatches": self._mm.value(),
                "three_prime_exact": self._end.value()}

    def _remember(self) -> None:
        v = self.values()
        st = QSettings()
        for key in ("product_min", "product_max", "max_mismatches", "three_prime_exact"):
            st.setValue(f"{_SETTINGS}/{key}", v[key])

    # ------------------------------------------------------------ buttons

    def _run(self) -> None:
        seq = self.sequence()
        if seq is None:
            QMessageBox.information(self, "Check Primers", NO_SEQUENCE)
            return
        if not seq.is_file():
            QMessageBox.warning(self, "Check Primers", f"The file does not exist:\n{seq}")
            return
        if not self._primers:
            QMessageBox.information(self, "Check Primers",
                                    "There are no primers to test in the Primer Database.")
            return
        self._remember()
        self.runRequested.emit(self.values())
        self.accept()

    def _all_settings(self) -> None:
        self._remember()
        self.allSettingsRequested.emit(self.values())
        self.accept()


__all__ = ["CheckPrimersDialog", "NO_SEQUENCE"]
