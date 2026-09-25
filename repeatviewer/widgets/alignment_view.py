"""Small alignment window opened from the tree viewer ("Show in alignment").

Paints only the visible cells, so alignments of thousands of sequences and
columns scroll smoothly. Names stay fixed on the left, a ruler on top; the
chosen sequence is highlighted and scrolled into view. "Differences only"
shows a dot where a residue equals the consensus, so an odd copy stands out.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QPoint, QRect, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontDatabase, QFontMetrics, QPainter, QPalette
from PySide6.QtWidgets import (
    QAbstractScrollArea,
    QCheckBox,
    QHBoxLayout,
    QLabel,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.tree_links import Alignment, consensus, first_residue
from repeatviewer.widgets.help_button import make_help_button

# Pale backgrounds: nucleotides (A green, C blue, G amber, T/U red), protein by
# physico-chemical group (Clustal-like). Text stays black on all of them.
_NUC = {"A": "#b7e4b0", "C": "#b3cdf2", "G": "#f6d8a0", "T": "#f2b3b3", "U": "#f2b3b3"}
_PROT_GROUPS = {
    "AVLIMFWC": "#b3cdf2",   # hydrophobic
    "KR": "#f2b3b3",         # positive
    "DE": "#e7b3f2",         # negative
    "NQST": "#b7e4b0",       # polar
    "G": "#f6d8a0",
    "P": "#f9f1a5",
    "HY": "#b3e8e6",
}
_PROT = {aa: color for group, color in _PROT_GROUPS.items() for aa in group}
_HIGHLIGHT = QColor(255, 196, 0, 90)

HELP = (
    "All sequences of the tree, aligned. The sequence you chose in the tree is "
    "highlighted and scrolled into view; choose another tip to move the highlight.<br>"
    "<b>Differences only</b>: a dot where a residue is the same as the consensus (the most "
    "common residue in that column), so the changes of an odd copy stand out.<br>"
    "Click a name to highlight that sequence here."
)


class AlignmentCanvas(QAbstractScrollArea):
    rowClicked = Signal(int)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
        font.setPointSize(max(9, font.pointSize()))
        self.viewport().setFont(font)
        self._aln: Alignment | None = None
        self._consensus = ""
        self._protein = False
        self._diff = False
        self._highlight: int | None = None
        self._set_metrics(font)

    def _set_metrics(self, font: QFont) -> None:
        fm = QFontMetrics(font)
        self._cw = max(7, fm.horizontalAdvance("W"))
        self._rh = fm.height() + 2
        self._ascent = fm.ascent() + 1
        self._ruler_h = self._rh + 4
        self._name_w = 120

    # -------------------------------------------------------------- data

    def set_alignment(self, aln: Alignment) -> None:
        self._aln = aln
        self._consensus = consensus(aln.seqs) if aln.aligned else ""
        text = "".join(aln.seqs[:50])
        letters = [c for c in text if c.isalpha()]
        nuc = sum(c in "ACGTUN" for c in letters)
        self._protein = bool(letters) and nuc / len(letters) < 0.9
        fm = QFontMetrics(self.viewport().font())
        longest = max((fm.horizontalAdvance(n) for n in aln.names), default=60)
        self._name_w = min(320, max(80, longest + 12))
        self._highlight = None
        self._update_scrollbars()
        self.viewport().update()

    def set_differences_only(self, on: bool) -> None:
        self._diff = bool(on)
        self.viewport().update()

    def highlight_row(self, row: int | None) -> None:
        self._highlight = row
        if row is not None and self._aln is not None:
            visible_rows = max(1, (self.viewport().height() - self._ruler_h) // self._rh)
            self.verticalScrollBar().setValue(max(0, row - visible_rows // 2))
            col = first_residue(self._aln.seqs[row])
            visible_cols = max(1, (self.viewport().width() - self._name_w) // self._cw)
            h = self.horizontalScrollBar()
            if not (h.value() <= col < h.value() + visible_cols):
                h.setValue(max(0, col - 5))
        self.viewport().update()

    @property
    def highlighted(self) -> int | None:
        return self._highlight

    # -------------------------------------------------------------- scrolling

    def _update_scrollbars(self) -> None:
        if self._aln is None:
            return
        vis_rows = max(1, (self.viewport().height() - self._ruler_h) // self._rh)
        vis_cols = max(1, (self.viewport().width() - self._name_w) // self._cw)
        v, h = self.verticalScrollBar(), self.horizontalScrollBar()
        v.setRange(0, max(0, len(self._aln.seqs) - vis_rows))
        v.setPageStep(vis_rows)
        h.setRange(0, max(0, self._aln.width - vis_cols))
        h.setPageStep(vis_cols)
        h.setSingleStep(5)

    def resizeEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().resizeEvent(event)
        self._update_scrollbars()

    def scrollContentsBy(self, dx: int, dy: int) -> None:  # noqa: N802 (Qt API)
        self.viewport().update()

    def mousePressEvent(self, event) -> None:  # noqa: N802 (Qt API)
        if self._aln is None:
            return
        pos = event.position().toPoint()
        if pos.y() < self._ruler_h:
            return
        row = self.verticalScrollBar().value() + (pos.y() - self._ruler_h) // self._rh
        if 0 <= row < len(self._aln.seqs):
            self._highlight = row
            self.viewport().update()
            self.rowClicked.emit(row)

    # -------------------------------------------------------------- painting

    def _color(self, c: str) -> str | None:
        return (_PROT if self._protein else _NUC).get(c)

    def paintEvent(self, event) -> None:  # noqa: N802 (Qt API)
        p = QPainter(self.viewport())
        pal = self.palette()
        base = pal.color(QPalette.ColorRole.Base)
        text = pal.color(QPalette.ColorRole.Text)
        mid = pal.color(QPalette.ColorRole.Mid)
        p.fillRect(self.viewport().rect(), base)
        aln = self._aln
        if aln is None:
            return
        vw, vh = self.viewport().width(), self.viewport().height()
        r0 = self.verticalScrollBar().value()
        c0 = self.horizontalScrollBar().value()
        n_rows = min(len(aln.seqs) - r0, (vh - self._ruler_h) // self._rh + 1)
        n_cols = (vw - self._name_w) // self._cw + 1
        x0 = self._name_w
        # ruler
        p.setPen(mid)
        for j in range(c0, min(aln.width, c0 + n_cols)):
            if (j + 1) % 10 == 0:
                x = x0 + (j - c0) * self._cw
                p.drawText(QRect(x - 40, 0, 40 + self._cw, self._ruler_h - 2),
                           int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignBottom),
                           str(j + 1))
        p.drawLine(0, self._ruler_h - 1, vw, self._ruler_h - 1)
        black = QColor("black")
        for k in range(max(0, n_rows)):
            row = r0 + k
            y = self._ruler_h + k * self._rh
            seq = aln.seqs[row]
            for j in range(c0, min(len(seq), c0 + n_cols)):
                c = seq[j]
                x = x0 + (j - c0) * self._cw
                same = self._diff and self._consensus and c == self._consensus[j]
                color = None if same else self._color(c)
                if color:
                    p.fillRect(x, y, self._cw, self._rh, QColor(color))
                    p.setPen(black)
                else:
                    p.setPen(mid if (same or c in "-.") else text)
                glyph = "·" if same and c not in "-." else c
                p.drawText(QPoint(x + 1, y + self._ascent), glyph)
            if row == self._highlight:
                p.fillRect(0, y, vw, self._rh, _HIGHLIGHT)
        # names on top of the residues that scrolled under them
        p.fillRect(0, self._ruler_h, self._name_w, vh, base)
        for k in range(max(0, n_rows)):
            row = r0 + k
            y = self._ruler_h + k * self._rh
            if row == self._highlight:
                p.fillRect(0, y, self._name_w, self._rh, _HIGHLIGHT)
            p.setPen(text)
            name = QFontMetrics(p.font()).elidedText(aln.names[row], Qt.TextElideMode.ElideRight,
                                                     self._name_w - 8)
            p.drawText(QPoint(4, y + self._ascent), name)
        p.setPen(mid)
        p.drawLine(self._name_w - 2, 0, self._name_w - 2, vh)


class AlignmentWindow(QWidget):
    """Top-level window: one alignment, reused while tips are chosen."""

    def __init__(self, path: Path, aln: Alignment, *, note: str = "",
                 parent: QWidget | None = None) -> None:
        super().__init__(parent, Qt.WindowType.Window)
        self.path = Path(path)
        self.setWindowTitle(f"Alignment: {self.path.name}" + (f" ({note})" if note else ""))
        self.resize(900, 420)
        self.canvas = AlignmentCanvas(self)
        self.canvas.set_alignment(aln)
        self._aln = aln
        self._diff = QCheckBox("Differences only", self)
        self._diff.setEnabled(aln.aligned)
        self._diff.toggled.connect(self.canvas.set_differences_only)
        self._status = QLabel(self)
        self._status.setStyleSheet("color: palette(mid);")
        self.canvas.rowClicked.connect(self._show_row)
        top = QHBoxLayout()
        top.addWidget(self._diff)
        top.addWidget(make_help_button(self, "Alignment", HELP))
        top.addSpacing(12)
        top.addWidget(self._status, 1)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.addLayout(top)
        lay.addWidget(self.canvas, 1)
        self._show_row(None)

    def _show_row(self, row: int | None) -> None:
        aln = self._aln
        info = f"{len(aln.seqs)} sequences, {aln.width} columns"
        if not aln.aligned:
            info += " (not aligned)"
        if row is not None:
            residues = sum(c not in "-." for c in aln.seqs[row])
            info = f"{aln.names[row]}: {residues} residues · " + info
        self._status.setText(info)

    def show_tip(self, tip: str) -> bool:
        """Highlight the row of a tree tip; False when it is not in this file."""
        row = self._aln.row_of(tip)
        if row is None:
            return False
        self.canvas.highlight_row(row)
        self._show_row(row)
        return True
