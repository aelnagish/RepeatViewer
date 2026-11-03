"""Tree viewer → Export → Distance matrix…

Pairwise distances between the tips: along the tree (patristic, the sum of the
branch lengths between two tips) or from the alignment (p-distance, share of
differing residues over the columns both sequences have). Shown as a table and
saved as a tab-separated file next to the tree.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QApplication,
    QDialog,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
)

from repeatviewer.backends import tree_links as tk
from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox
from repeatviewer.widgets.help_button import make_help_button

MAX_SHOWN = 400          # larger matrices are saved, not shown
SOURCE_TREE = "tree"
SOURCE_ALIGNMENT = "alignment"
# Default dialog size: half of the previous 760×520 so long tip names
# do not dominate the Preview; the window and columns stay drag-resizable.
_DIALOG_WIDTH = 380
_DIALOG_HEIGHT = 520
_MIN_COL_WIDTH = 56
_MIN_ROW_LABEL_WIDTH = 72

HELP = (
    "<b>Tree (patristic)</b>: the sum of the branch lengths on the path between two tips, "
    "in substitutions per site, for the tree as shown (rerooting does not change it).<br>"
    "<b>Alignment (p-distance)</b>: the share of positions where two sequences differ, "
    "counted over the columns both have (gaps and N skipped). NA = no shared column.<br>"
    "<b>Save as TSV</b> writes the full matrix next to the tree."
)


class TreeDistancesDialog(QDialog):
    fileWritten = Signal(object)

    def __init__(self, tree_path: Path, model, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle(f"Distance matrix: {Path(tree_path).name}")
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, False)
        self.setModal(False)
        self.setWindowModality(Qt.WindowModality.NonModal)
        self.resize(_DIALOG_WIDTH, _DIALOG_HEIGHT)
        self.setMinimumSize(300, 240)
        self.setSizeGripEnabled(True)
        self._tree_path = Path(tree_path)
        self._model = model
        self._source = tk.sequences_for_tree(self._tree_path)
        self._names: list[str] = []
        self._matrix = np.zeros((0, 0))

        self._source_box = QComboBox(self)
        self._source_box.addItem("Tree (patristic)", SOURCE_TREE)
        if self._source is not None and self._source.aligned:
            self._source_box.addItem(f"Alignment (p-distance): {self._source.path.name}",
                                     SOURCE_ALIGNMENT)
        self._source_box.currentIndexChanged.connect(self._compute)
        self._note = QLabel(self)
        self._note.setStyleSheet("color: palette(mid);")
        self._table = QTableWidget(self)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        hh.setStretchLastSection(False)
        hh.setMinimumSectionSize(_MIN_COL_WIDTH)
        vh = self._table.verticalHeader()
        vh.setDefaultAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        )
        copy_btn = QPushButton("Copy", self)
        copy_btn.clicked.connect(self._copy)
        save_btn = QPushButton("Save as TSV…", self)
        save_btn.clicked.connect(self._save)
        close_btn = QPushButton("Close", self)
        close_btn.clicked.connect(self.close)

        top = QHBoxLayout()
        top.addWidget(QLabel("Distances:", self))
        top.addWidget(self._source_box, 1)
        top.addWidget(make_help_button(self, "Distance matrix", HELP))
        buttons = QHBoxLayout()
        buttons.addWidget(self._note, 1)
        buttons.addWidget(copy_btn)
        buttons.addWidget(save_btn)
        buttons.addWidget(close_btn)
        lay = QVBoxLayout(self)
        lay.addLayout(top)
        lay.addWidget(self._table, 1)
        lay.addLayout(buttons)
        self._compute()

    @property
    def source(self) -> str:
        return str(self._source_box.currentData())

    def matrix(self) -> tuple[list[str], np.ndarray]:
        return self._names, self._matrix

    def _apply_half_widths(self) -> None:
        """Data columns at half content width; tip-name column full width (figure 1)."""
        n = self._table.columnCount()
        if n <= 0:
            return
        self._table.resizeColumnsToContents()
        for c in range(n):
            self._table.setColumnWidth(
                c, max(_MIN_COL_WIDTH, self._table.columnWidth(c) // 2)
            )
        vh = self._table.verticalHeader()
        fm = vh.fontMetrics()
        label_w = max(
            (fm.horizontalAdvance(name) for name in self._names),
            default=_MIN_ROW_LABEL_WIDTH,
        ) + 16
        # Row labels like the original matrix: wide enough for the tip names.
        # Cap so a long name does not consume the whole half-width dialog.
        cap = max(200, int(self.width() * 0.65))
        vh.setFixedWidth(max(_MIN_ROW_LABEL_WIDTH, min(label_w, cap)))

    def _compute(self, *_args) -> None:
        if self.source == SOURCE_ALIGNMENT and self._source is not None:
            aln = tk.read_alignment(self._source.path)
            tips = {n.name for n in self._model.nodes if n.is_tip}
            rows = [i for i, name in enumerate(aln.names)
                    if name in tips or tk.tip_key(name) in tips]
            aln = tk.Alignment([aln.names[i] for i in rows], [aln.seqs[i] for i in rows], True)
            self._names, self._matrix = tk.p_distance_matrix(aln)
        else:
            self._names, self._matrix = tk.patristic_matrix(self._model)
        n = len(self._names)
        if n > MAX_SHOWN:
            self._table.clear()
            self._table.setRowCount(0)
            self._table.setColumnCount(0)
            self._note.setText(f"{n} sequences: too many to show here; use Save as TSV.")
            return
        self._note.setText(f"{n} sequences")
        self._table.setRowCount(n)
        self._table.setColumnCount(n)
        self._table.setHorizontalHeaderLabels(self._names)
        self._table.setVerticalHeaderLabels(self._names)
        for i, name in enumerate(self._names):
            for header_item in (
                self._table.horizontalHeaderItem(i),
                self._table.verticalHeaderItem(i),
            ):
                if header_item is not None:
                    header_item.setToolTip(name)
        for a in range(n):
            for b in range(n):
                v = self._matrix[a, b]
                item = QTableWidgetItem("NA" if not np.isfinite(v) else f"{v:.4g}")
                item.setTextAlignment(
                    Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
                )
                self._table.setItem(a, b, item)
        self._apply_half_widths()

    def default_path(self) -> Path:
        kind = "patristic" if self.source == SOURCE_TREE else "pdistance"
        return self._tree_path.parent / f"{self._tree_path.name}_distances_{kind}.tsv"

    def _copy(self) -> None:
        QApplication.clipboard().setText(tk.matrix_tsv(self._names, self._matrix))
        self._note.setText("Copied: paste into a spreadsheet.")

    def save_to(self, path: Path) -> Path:
        path = Path(path)
        path.write_text(tk.matrix_tsv(self._names, self._matrix), encoding="utf-8")
        self.fileWritten.emit(path)
        return path

    def _save(self) -> None:
        chosen, _ = QFileDialog.getSaveFileName(self, "Save distance matrix",
                                                str(self.default_path()),
                                                "Tab-separated (*.tsv);;All files (*)")
        if chosen:
            path = self.save_to(Path(chosen))
            self._note.setText(f"Saved {path.name}")
