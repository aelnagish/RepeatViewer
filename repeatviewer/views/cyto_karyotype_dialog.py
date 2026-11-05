"""Edit karyotype: chromosomes of a species (number, lengths, centromere index, copies, bands)."""

from __future__ import annotations

import copy
from pathlib import Path

from PySide6.QtCore import QEvent, QObject, QPointF, QRectF, QSettings, QSize, Qt, QThread, QTimer, Signal
from PySide6.QtGui import QBrush, QColor, QFont, QIcon, QKeySequence, QPainter, QPainterPath, QPen, QPixmap
from PySide6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QCheckBox,
    QColorDialog,
    QDialog,
    QFrame,
    QGridLayout,
    QMenu,
    QScrollArea,
    QSizePolicy,
    QToolButton,
    QFileDialog,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QSpinBox,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.cytogenetics import karyotype as kt
from repeatviewer.backends.cytogenetics import lists as lists_mod
from repeatviewer.backends.cytogenetics.lists import load_lists
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cyto_panel_targets import PanelPathFields
from repeatviewer.widgets.cytogenetics_style import rounded_box
from repeatviewer.widgets.help_button import make_help_button

HELP = (
    "<b>Karyotype</b><br>"
    "One set of chromosomes for a species. It is saved in the folder and used by Signals, Idiogram and Figure.<br>"
    "<b>Species</b>: choose one or type a new name; it is added to the species list on Save.<br>"
    "<b>Basic number x and Ploidy</b> apply as you type; chromosomes are added or removed at the end of the list.<br>"
    "<b>Homologues</b>: how the copies of a chromosome are named (1a, 1b or 1.1, 1.2 or your own list). "
    "Signals marked under the old names keep them, so choose this before marking.<br>"
    "<b>Unit</b>: µm, % of the haploid genome, or Mb from an assembly.<br>"
    "<b>Table</b>: double-click a cell to edit it. p and q are the arm lengths; CI (centromere index) is the short arm "
    "as percent of the chromosome, and the type (m, sm, st, t) follows Levan et al. (1964). Copies can differ from the "
    "ploidy (aneuploidy). Drag ⋮⋮ to change the order; a blue dot marks a chromosome whose homologues differ.<br>"
    "<b>Import</b>: Measurements (TSV) reads chromosome and length (and CI, or p and q), for example from KaryoScope. "
    "Lengths from Assembly reads a sizes file, .fai, BED, GFF or FASTA and matches Chr1, Chr01 … to chromosome 1.<br>"
    "<b>Preview</b>: click a homologue to edit it below. Drag the round centromere handle, the bottom end (length), "
    "a band (move) or its top or bottom edge (resize). <b>Draw band</b>: drag along a chromosome; the colour boxes "
    "choose the kind (NOR, DAPI+, CMA+, FISH probes) or, with a band selected, change it; + picks an own colour. "
    "Backspace removes the selected band.<br>"
    "<b>Same as the other homologues</b> (ticked by default): edits apply to every copy. Untick it to give one "
    "homologue its own length, centromere and bands, for example a heteromorphic pair.<br>"
    "<b>Right-click</b> a homologue: Mirror Bands (an inversion), Copy Bands to the Other Homologues, "
    "Duplicate as B Chromosome, Mark as Missing (aneuploidy). Right-click a NOR band: Draw as Satellite shows the "
    "secondary constriction with the satellite beyond it.<br>"
    "<b>Presets</b>: save the bands of this karyotype as a set for the species, and apply it to another karyotype.<br>"
    "<b>Undo</b>: Cmd+Z and Cmd+Shift+Z undo and redo changes inside this window. Nothing is written before Save."
)

HELP_ASSEMBLY = (
    "<b>Source</b>: a sizes file (name and length per line, for example <i>genome.fa_sizes.txt</i>), a .fai index, "
    "a BED or GFF file (the last end per sequence) or a FASTA file.<br>"
    "<b>FASTA</b> is read in the background and a sizes file is saved beside it, so the next time is instant. "
    "Tools > Cytogenetics > Make Chromosome Sizes File makes one for any FASTA.<br>"
    "<b>From R</b>: <i>rtracklayer::export(gr, \"x.bed\")</i>, or "
    "<i>write.table(data.frame(names(seqlengths(gr)), seqlengths(gr)), \"x.sizes\", sep = \"\\t\", "
    "quote = FALSE, row.names = FALSE, col.names = FALSE)</i>.<br>"
    "<b>Matching</b>: the current link, the same name, or Chr1, chr01, Csat_chr1 … for chromosome 1. "
    "Choose another sequence in the menu to change it. Save sets the lengths in Mb."
)
UNITS = (("µm", "µm"), ("% of genome (relative)", "%"), ("Mb (assembly)", "Mb"))
PLOIDIES = (1, 2, 3, 4, 5, 6, 8)
BAND_COLOURS = {
    "DAPI+ heterochromatin": "#3B6FD8",
    "CMA+ band": "#E8C23D",
}
CUSTOM_KIND = "Custom…"
OWN_LIST = "Own list…"
LINK_CSS = "QPushButton { color:#1F4E9A; border:none; } QPushButton:hover { text-decoration: underline; }"


class _Edit(QLineEdit):
    focused = Signal()

    def focusInEvent(self, event) -> None:  # noqa: N802
        self.focused.emit()
        super().focusInEvent(event)


class _ReorderTable(QTableWidget):
    rowMoved = Signal(int, int)

    def __init__(self, *args) -> None:
        super().__init__(*args)
        self.setDragEnabled(True)
        self.setAcceptDrops(True)
        self.viewport().setAcceptDrops(True)
        self.setDragDropMode(QTableWidget.DragDropMode.InternalMove)
        self.setDropIndicatorShown(True)
        self.setSelectionBehavior(QTableWidget.SelectionBehavior.SelectRows)
        self.setSelectionMode(QTableWidget.SelectionMode.SingleSelection)

    def dropEvent(self, event) -> None:  # noqa: N802
        src = self.currentRow()
        dst = self.rowAt(int(event.position().y()))
        if dst < 0:
            dst = self.rowCount() - 1
        event.setDropAction(Qt.DropAction.IgnoreAction)
        event.ignore()
        if src >= 0 and dst != src:
            self.rowMoved.emit(src, dst)


class KaryoPreview(QWidget):
    """Chromosome bars side by side, one per type or all homologues in groups. Drag a centromere or a band edge."""

    selected = Signal(int)
    edited = Signal()

    BAR_W = 22
    MAX_H = 300

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumSize(300, 340)
        self.setMouseTracking(True)
        self.karyotype: kt.Karyotype | None = None
        self.index = -1
        self.band_index = -1
        self.groups = False
        self._drag: tuple[str, int] | None = None
        self._hover = -1

    def _bars(self, i: int, chrom: kt.KaryoChromosome) -> list[QRectF]:
        """One rectangle per drawn bar of chromosome ``i`` (several in homologue groups)."""
        k = self.karyotype
        n = max(1, len(k.chromosomes))
        slot = max(30.0, (self.width() - 16) / n)
        longest = max((c.length for c in k.chromosomes), default=0.0)
        bottom = 44 if self.groups else 26
        h = self.MAX_H * (chrom.length / longest) if longest > 0 and chrom.length > 0 else self.MAX_H * 0.5
        h = max(30.0, min(h, self.height() - bottom - 14.0))
        copies = max(1, k.copies_of(chrom)) if self.groups else 1
        gap = 4.0
        w = self.BAR_W if copies == 1 else max(6.0, min(16.0, (slot - 10 - gap * (copies - 1)) / copies))
        total = copies * w + gap * (copies - 1)
        x0 = 8 + slot * (i + 0.5) - total / 2
        top = self.height() - bottom - h
        return [QRectF(x0 + j * (w + gap), top, w, h) for j in range(copies)]

    def band_colour(self, band: kt.Band) -> QColor:
        return QColor(band.color or BAND_COLOURS.get(band.kind, "#8A8D93"))

    def _group_rect(self, i: int, chrom: kt.KaryoChromosome) -> QRectF:
        bars = self._bars(i, chrom)
        return QRectF(bars[0].topLeft(), bars[-1].bottomRight())

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.fillRect(self.rect(), QColor("#FFFFFF"))
        k = self.karyotype
        if k is None or not k.chromosomes:
            p.setPen(QColor("#5A5D63"))
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "No chromosomes")
            return
        accent = self.palette().color(self.palette().ColorRole.Highlight)
        bold = QFont(self.font())
        bold.setBold(True)
        small = QFont(self.font())
        small.setPointSizeF(max(7.0, small.pointSizeF() - 2))
        for i, c in enumerate(k.chromosomes):
            bars = self._bars(i, c)
            group = self._group_rect(i, c)
            if i == self.index:
                p.setPen(QPen(accent, 1.2, Qt.PenStyle.DashLine))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRoundedRect(group.adjusted(-8, -8, 8, 22 if not self.groups else 40), 8, 8)
            for r in bars:
                radius = r.width() / 2
                for b_i, b in enumerate(c.bands):
                    top = r.top() + b.start * r.height()
                    rect = QRectF(r.left(), top, r.width(), max(2.0, (b.end - b.start) * r.height()))
                    if b.kind.startswith("Secondary"):
                        brush = QBrush(QColor("#8A8D93"), Qt.BrushStyle.BDiagPattern)
                    else:
                        brush = QBrush(self.band_colour(b))
                    p.setPen(Qt.PenStyle.NoPen)
                    p.setBrush(brush)
                    if kt.manual_probe(b.kind):
                        cy_b = (rect.top() + rect.bottom()) / 2
                        p.setPen(QPen(self.band_colour(b), 3))
                        p.setBrush(QColor("#FFFFFF"))
                        p.drawEllipse(QPointF(r.center().x(), cy_b), 6, 6)
                    else:
                        p.drawRect(rect)
                p.setPen(QPen(QColor("#6C6F75"), 1.5))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRoundedRect(r, radius, radius)
                cy = r.top() + c.centromere_fraction * r.height()
                if i == self.index or i == self._hover:
                    # the handle you drag: wider, in the accent colour
                    p.setPen(QPen(accent, 1.5))
                    p.setBrush(accent if i == self.index else QColor("#FFFFFF"))
                    p.drawRoundedRect(QRectF(r.left() - 6, cy - 5, r.width() + 12, 10), 5, 5)
                else:
                    p.setBrush(QColor("#FFFFFF"))
                    p.drawRoundedRect(QRectF(r.left() - 3, cy - 3, r.width() + 6, 6), 3, 3)
            label_y = group.bottom() + 4
            p.setPen(QColor("#1C1C1E"))
            p.setFont(bold)
            p.drawText(QRectF(group.center().x() - 40, label_y, 80, 18), Qt.AlignmentFlag.AlignCenter, c.name)
            if self.groups:
                p.setFont(small)
                p.setPen(QColor("#6E6E73"))
                p.drawText(QRectF(group.center().x() - 45, label_y + 17, 90, 16), Qt.AlignmentFlag.AlignCenter,
                           kt.group_name(k.copies_of(c)))
            if i == self.index and 0 <= self.band_index < len(c.bands):
                b = c.bands[self.band_index]
                p.setPen(QPen(accent, 2))
                for f in (b.start, b.end):
                    y = group.top() + f * group.height()
                    p.drawLine(QPointF(group.left() - 6, y), QPointF(group.right() + 6, y))

    def _hit(self, pos: QPointF) -> tuple[str, int] | None:
        k = self.karyotype
        if k is None:
            return None
        if 0 <= self.index < len(k.chromosomes):
            c = k.chromosomes[self.index]
            g = self._group_rect(self.index, c)
            if 0 <= self.band_index < len(c.bands) and g.left() - 10 <= pos.x() <= g.right() + 10:
                b = c.bands[self.band_index]
                for edge, f in (("start", b.start), ("end", b.end)):
                    if abs(pos.y() - (g.top() + f * g.height())) <= 4:
                        return edge, self.index
        for i, c in enumerate(k.chromosomes):
            g = self._group_rect(i, c)
            cy = g.top() + c.centromere_fraction * g.height()
            if g.left() - 8 <= pos.x() <= g.right() + 8 and abs(pos.y() - cy) <= 7:
                return "centromere", i
        for i, c in enumerate(k.chromosomes):
            if self._group_rect(i, c).adjusted(-8, -8, 8, 8).contains(pos):
                return "select", i
        return None

    def mousePressEvent(self, event) -> None:  # noqa: N802
        hit = self._hit(event.position())
        if hit is None:
            return
        kind, i = hit
        if i != self.index:
            self.index = i
            self.band_index = 0 if self.karyotype.chromosomes[i].bands else -1
            self.selected.emit(i)
        self._drag = None if kind == "select" else (kind, i)
        self.update()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        k = self.karyotype
        if self._drag is None or k is None:
            hit = self._hit(event.position())
            hover = hit[1] if hit and hit[0] == "centromere" else -1
            if hover != self._hover:
                self._hover = hover
                self.update()
            self.setCursor(Qt.CursorShape.SizeVerCursor if hit and hit[0] != "select" else Qt.CursorShape.ArrowCursor)
            return
        kind, i = self._drag
        c = k.chromosomes[i]
        g = self._group_rect(i, c)
        f = max(0.0, min(1.0, (event.position().y() - g.top()) / max(g.height(), 1.0)))
        if kind == "centromere":
            c.centromere_index = round(max(5.0, min(95.0, f * 100.0)), 1)
        elif 0 <= self.band_index < len(c.bands):
            b = c.bands[self.band_index]
            if kind == "start":
                b.start = min(f, b.end)
            else:
                b.end = max(f, b.start)
        self.edited.emit()
        self.update()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        self._drag = None

    def leaveEvent(self, event) -> None:  # noqa: N802
        if self._hover != -1:
            self._hover = -1
            self.update()


class KaryotypeCanvas(QWidget):
    """Preview of Edit karyotype: every homologue (or one per chromosome), editable by dragging.

    Select: click a homologue, a band or a chromosome name; drag a band to move it, its top or bottom edge to resize
    it, the round centromere handle or the bottom end of a chromosome. Draw band: drag along a chromosome.
    """

    selectionChanged = Signal()
    aboutToEdit = Signal()
    edited = Signal()
    menuRequested = Signal(object)  # global QPoint

    BAR_W = 22.0
    HOM_GAP = 8.0
    GROUP_GAP = 36.0
    LABEL_W = 26.0
    TOP = 30.0

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumSize(420, 300)
        self.setMouseTracking(True)
        self.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
        self.karyotype: kt.Karyotype | None = None
        self.every = True
        self.tool = "select"
        self.draw_kind = kt.BAND_KINDS[0]
        self.draw_colour = ""
        self.probe_colours: dict[str, str] = {}
        self.ci = -1  # chromosome
        self.hi = -1  # homologue (-1 = the whole chromosome)
        self.bi = -1  # band of the selected homologue
        self._items: list[tuple[int, int, QRectF, kt.KaryoChromosome, bool]] = []
        self._labels: list[tuple[int, QRectF]] = []
        self._scale = 1.0
        self._drag: dict | None = None
        self._ghost: tuple[int, int, float, float] | None = None

    # ------------------------------------------------------------------ geometry

    def _layout(self) -> None:
        self._items, self._labels = [], []
        k = self.karyotype
        if k is None or not k.chromosomes:
            return
        groups = []
        for ci, c in enumerate(k.chromosomes):
            n = k.copies_of(c) if self.every else 1
            groups.append((ci, c, max(1, n)))
        bar, gap, ggap, lab = self.BAR_W, self.HOM_GAP, self.GROUP_GAP, self.LABEL_W
        need = sum(lab + n * bar + (n - 1) * gap for _ci, _c, n in groups) + ggap * (len(groups) - 1)
        avail = max(100.0, self.width() - 40.0)
        f = min(1.0, avail / need) if need > 0 else 1.0
        bar, gap, ggap, lab = max(6.0, bar * f), gap * f, ggap * f, max(14.0, lab * f)
        need = sum(lab + n * bar + (n - 1) * gap for _ci, _c, n in groups) + ggap * (len(groups) - 1)
        lengths = [kt.homologue_view(c, i).length for _ci, c, n in groups for i in range(n)]
        longest = max([v for v in lengths if v > 0] or [1.0])
        self._scale = max(0.1, (self.height() - self.TOP - 56.0) / longest)
        x = max(20.0, (self.width() - need) / 2)
        for ci, c, n in groups:
            heights = []
            for i in range(n):
                view = kt.homologue_view(c, i if self.every else -1)
                heights.append((view, (view.length if view.length > 0 else longest) * self._scale))
            gh = max(h for _v, h in heights)
            self._labels.append((ci, QRectF(x, self.TOP, lab - 4, gh)))
            x += lab
            for i, (view, h) in enumerate(heights):
                hi = i if self.every else -1
                missing = self.every and kt.is_missing(c, i)
                self._items.append((ci, hi, QRectF(x, self.TOP, bar, h), view, missing))
                x += bar + gap
            x += ggap - gap

    def resizeEvent(self, event) -> None:  # noqa: N802
        self._layout()
        super().resizeEvent(event)

    def refresh(self) -> None:
        self._layout()
        self.update()

    def _item(self, ci: int, hi: int):
        return next((it for it in self._items if it[0] == ci and it[1] == hi), None)

    # ------------------------------------------------------------------ painting

    def band_brush(self, band: kt.Band) -> QBrush:
        probe = kt.manual_probe(band.kind)
        if band.color:
            return QBrush(QColor(band.color))
        if probe:
            return QBrush(QColor(self.probe_colours.get(probe, "#2FB457")))
        if band.kind.startswith("Secondary"):
            return QBrush(QColor("#8A8D93"), Qt.BrushStyle.BDiagPattern)
        return QBrush(QColor(BAND_COLOURS.get(band.kind, "#8A8D93")))

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.fillRect(self.rect(), QColor("#FFFFFF"))
        k = self.karyotype
        if k is None or not k.chromosomes:
            p.setPen(QColor("#5A5D63"))
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Add a chromosome with + under the table")
            return
        outline = QColor("#6C6F75")
        blue = QColor("#1F4E9A")
        small = QFont(self.font())
        small.setPointSizeF(max(8.0, small.pointSizeF() - 2))
        bold = QFont(self.font())
        bold.setBold(True)
        bold.setPointSizeF(bold.pointSizeF() + 1)
        for ci, hi, r, view, missing in self._items:
            path = chromosome_path(r, view)
            if missing:
                p.setPen(QPen(QColor("#B0B2B8"), 1.2, Qt.PenStyle.DashLine))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawPath(path)
            else:
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(QColor("#FFFFFF"))
                p.drawPath(path)
                p.save()
                p.setClipPath(path)
                for b in view.bands:
                    p.setBrush(self.band_brush(b))
                    p.drawRect(QRectF(r.left(), r.top() + b.start * r.height(), r.width(),
                                      max(1.5, (b.end - b.start) * r.height())))
                if self._ghost and self._ghost[:2] == (ci, hi):
                    a, z = sorted(self._ghost[2:])
                    p.setBrush(QColor(255, 176, 60, 140))
                    p.drawRect(QRectF(r.left(), r.top() + a * r.height(), r.width(), max(2.0, (z - a) * r.height())))
                p.restore()
                p.setPen(QPen(outline, 1.4))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawPath(path)
                cy = r.top() + view.centromere_fraction * r.height()
                p.setBrush(QColor("#FFFFFF"))
                p.drawEllipse(QPointF(r.center().x(), cy), 5, 5)
            if self.every:
                p.setFont(small)
                p.setPen(QColor("#5A5D63"))
                suffix = k.homologue_label(k.chromosomes[ci], hi)[len(k.chromosomes[ci].name):] or k.chromosomes[ci].name
                p.drawText(QRectF(r.left() - 10, r.bottom() + 10, r.width() + 20, 16), Qt.AlignmentFlag.AlignHCenter, suffix)
            if ci == self.ci and (hi == self.hi or not self.every):
                p.setPen(QPen(QColor("#2F6FDB"), 1.6, Qt.PenStyle.DashLine))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRoundedRect(r.adjusted(-7, -8, 7, 8), 10, 10)
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(blue)
                p.drawRoundedRect(QRectF(r.center().x() - 6, r.bottom() - 3, 12, 7), 3, 3)
                if 0 <= self.bi < len(view.bands):
                    b = view.bands[self.bi]
                    y0, y1 = r.top() + b.start * r.height(), r.top() + b.end * r.height()
                    p.setPen(QPen(blue, 1.6))
                    p.drawRect(QRectF(r.left() - 4, y0 - 3, r.width() + 8, y1 - y0 + 6))
                    p.setPen(Qt.PenStyle.NoPen)
                    for y in (y0 - 3, y1 + 3):
                        p.drawRoundedRect(QRectF(r.center().x() - 8, y - 2.5, 16, 5), 2, 2)
        p.setFont(bold)
        p.setPen(QColor("#1D1D1F"))
        for ci, rect in self._labels:
            if self.every and ci == self.ci and self.hi < 0:
                p.setPen(QColor("#1F4E9A"))
            p.drawText(rect, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, k.chromosomes[ci].name)
            p.setPen(QColor("#1D1D1F"))
        self._paint_scale(p, small)

    def _paint_scale(self, p: QPainter, font: QFont) -> None:
        k = self.karyotype
        if k is None or not any(c.length > 0 for c in k.chromosomes):
            return
        value = next((v for v in (0.1, 0.2, 0.5, 1, 2, 4, 5, 10, 20, 50, 100, 200, 500) if v * self._scale >= 50), 500)
        w = value * self._scale
        x1, y = self.width() - 20.0, self.height() - 30.0
        p.setPen(QPen(QColor("#1D1D1F"), 2))
        p.drawLine(QPointF(x1 - w, y), QPointF(x1, y))
        p.setFont(font)
        p.setPen(QColor("#5A5D63"))
        p.drawText(QRectF(x1 - w - 30, y + 3, w + 60, 16), Qt.AlignmentFlag.AlignHCenter, f"{value:g} {k.unit}")

    # ------------------------------------------------------------------ mouse

    def _hit(self, pos: QPointF) -> tuple[str, int, int, int]:
        """(what, chromosome, homologue, band): what is cen, end, edge0, edge1, band, hom, label or ''."""
        for ci, rect in self._labels:
            if rect.adjusted(-4, -4, 4, 4).contains(pos) and self.every:
                return "label", ci, -1, -1
        for ci, hi, r, view, missing in self._items:
            if not r.adjusted(-8, -10, 8, 10).contains(pos):
                continue
            selected = ci == self.ci and (hi == self.hi or not self.every)
            if missing:
                return "hom", ci, hi, -1
            cy = r.top() + view.centromere_fraction * r.height()
            if abs(pos.y() - cy) <= 7:
                return "cen", ci, hi, -1
            if selected and abs(pos.y() - r.bottom()) <= 6:
                return "end", ci, hi, -1
            f = (pos.y() - r.top()) / max(1.0, r.height())
            if selected and 0 <= self.bi < len(view.bands):
                b = view.bands[self.bi]
                for what, edge in (("edge0", b.start), ("edge1", b.end)):
                    if abs(pos.y() - (r.top() + edge * r.height())) <= 5:
                        return what, ci, hi, self.bi
            for bi in range(len(view.bands) - 1, -1, -1):
                b = view.bands[bi]
                if b.start - 0.004 <= f <= b.end + 0.004:
                    return "band", ci, hi, bi
            return "hom", ci, hi, -1
        return "", -1, -1, -1

    def _owner(self):
        c = self.karyotype.chromosomes[self.ci]
        return kt.target(c, self.hi)

    def _select(self, ci: int, hi: int, bi: int) -> None:
        if (ci, hi, bi) != (self.ci, self.hi, self.bi):
            self.ci, self.hi, self.bi = ci, hi, bi
            self.selectionChanged.emit()
        self.update()

    def mousePressEvent(self, event) -> None:  # noqa: N802
        pos = event.position()
        what, ci, hi, bi = self._hit(pos)
        if event.button() == Qt.MouseButton.RightButton:
            if what:
                self._select(ci, hi, bi if what in ("band", "edge0", "edge1") else -1)
                self.menuRequested.emit(event.globalPosition().toPoint())
            return
        if event.button() != Qt.MouseButton.LeftButton:
            return
        if not what:
            self._select(-1, -1, -1)
            return
        if what == "label":
            self._select(ci, -1, -1)
            return
        if self.tool == "draw" and what != "label" and not (self.every and kt.is_missing(self.karyotype.chromosomes[ci], hi)):
            self._select(ci, hi, -1)
            r = self._item(ci, hi)[2]
            f = max(0.0, min(1.0, (pos.y() - r.top()) / r.height()))
            self._ghost = (ci, hi, f, f)
            self._drag = {"what": "draw"}
            return
        self._select(ci, hi, bi if what in ("band", "edge0", "edge1") else -1)
        if what in ("cen", "end", "band", "edge0", "edge1"):
            r = self._item(ci, hi)[2]
            self._drag = {"what": what, "y0": pos.y(), "r": QRectF(r), "started": False}
            if what == "band":
                b = self._owner().bands[self.bi]
                self._drag["band"] = (b.start, b.end)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        pos = event.position()
        d = self._drag
        if d is None:
            what = self._hit(pos)[0]
            shape = {"cen": Qt.CursorShape.SizeVerCursor, "end": Qt.CursorShape.SizeVerCursor,
                     "edge0": Qt.CursorShape.SizeVerCursor, "edge1": Qt.CursorShape.SizeVerCursor,
                     "band": Qt.CursorShape.OpenHandCursor}.get(what)
            if self.tool == "draw" and what in ("hom", "band", "cen"):
                shape = Qt.CursorShape.CrossCursor
            self.setCursor(shape or Qt.CursorShape.ArrowCursor)
            return
        if d["what"] == "draw":
            r = self._item(self._ghost[0], self._ghost[1])[2]
            f = max(0.0, min(1.0, (pos.y() - r.top()) / r.height()))
            self._ghost = (*self._ghost[:3], f)
            self.update()
            return
        if not d["started"]:
            if abs(pos.y() - d["y0"]) < 2:
                return
            d["started"] = True
            self.aboutToEdit.emit()
        r: QRectF = d["r"]
        owner = self._owner()
        f = (pos.y() - r.top()) / max(1.0, r.height())
        if d["what"] == "cen":
            owner.centromere_index = round(max(0.5, min(50.0, 100.0 * f)), 2)
        elif d["what"] == "end":
            owner.length = round(max(0.1, (pos.y() - r.top()) / self._scale), 2)
        elif d["what"] == "band":
            s0, e0 = d["band"]
            delta = (pos.y() - d["y0"]) / max(1.0, r.height())
            delta = max(-s0, min(1.0 - e0, delta))
            b = owner.bands[self.bi]
            b.start, b.end = round(s0 + delta, 4), round(e0 + delta, 4)
        else:
            b = owner.bands[self.bi]
            f = max(0.0, min(1.0, f))
            if d["what"] == "edge0":
                b.start = round(min(f, b.end - 0.005), 4)
            else:
                b.end = round(max(f, b.start + 0.005), 4)
        self._layout()
        self.update()
        self.edited.emit()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        d, self._drag = self._drag, None
        if d and d["what"] == "draw" and self._ghost:
            ci, hi, a, z = self._ghost
            self._ghost = None
            a, z = sorted((a, z))
            if z - a < 0.01:
                a, z = max(0.0, a - 0.02), min(1.0, z + 0.02)
            self.aboutToEdit.emit()
            owner = self._owner()
            satellite = self.draw_kind.startswith("Secondary") and (a < 0.25 or z > 0.75)
            owner.bands.append(kt.Band(self.draw_kind, round(a, 4), round(z, 4), self.draw_colour, satellite))
            self.bi = len(owner.bands) - 1
            self.update()
            self.edited.emit()
            self.selectionChanged.emit()

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        what, ci, hi, _bi = self._hit(event.position())
        if what in ("hom", "band") and self.every:
            self._select(ci, -1, -1)

    def keyPressEvent(self, event) -> None:  # noqa: N802
        if event.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace) and self.ci >= 0 and self.bi >= 0:
            self.aboutToEdit.emit()
            del self._owner().bands[self.bi]
            self.bi = -1
            self.update()
            self.edited.emit()
            self.selectionChanged.emit()
            return
        super().keyPressEvent(event)


def chromosome_path(r: QRectF, chrom: kt.KaryoChromosome, centromere_gap: float = 3.0) -> QPainterPath:
    """Outline of a chromosome: arms split at the centromere and, for NOR bands drawn as satellite, a narrow neck
    (secondary constriction) between the satellite and the rest of the arm."""
    h = r.height()
    cuts: list[tuple[float, float, bool]] = []
    if centromere_gap > 0:
        cy = r.top() + chrom.centromere_fraction * h
        cuts.append((cy - centromere_gap, cy + centromere_gap, False))
    for b in chrom.bands:
        if not b.satellite or b.end - b.start <= 0:
            continue
        span = b.end - b.start
        if b.start < 1.0 - b.end:  # satellite at the p end: neck on the side of the centromere
            n0, n1 = b.start + 0.55 * span, b.end
        else:
            n0, n1 = b.start, b.start + 0.45 * span
        cuts.append((r.top() + n0 * h, r.top() + n1 * h, True))
    cuts.sort()
    path = QPainterPath()
    y = r.top()
    for y0, y1, neck in cuts + [(r.bottom(), r.bottom(), False)]:
        y0 = max(y, y0)
        if y0 - y > 1:
            rad = min(r.width() / 2, (y0 - y) / 2)
            path.addRoundedRect(QRectF(r.left(), y, r.width(), y0 - y), rad, rad)
        if neck and y1 > y0:
            nw = max(2.0, r.width() * 0.34)
            path.addRect(QRectF(r.center().x() - nw / 2, y0 - 0.5, nw, y1 - y0 + 1))
        y = max(y, y1)
    return path


def _swatch_icon(brush: QBrush, size: int = 16) -> QIcon:
    pm = QPixmap(size, size)
    pm.fill(Qt.GlobalColor.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(QColor("#FFFFFF"))
    p.drawRoundedRect(QRectF(0.5, 0.5, size - 1, size - 1), 3, 3)
    p.setBrush(brush)
    p.setPen(QPen(QColor("#B5B5BA"), 1))
    p.drawRoundedRect(QRectF(0.5, 0.5, size - 1, size - 1), 3, 3)
    p.end()
    return QIcon(pm)


def _dot_icon(colour: str = "#2F6FDB") -> QIcon:
    pm = QPixmap(12, 12)
    pm.fill(Qt.GlobalColor.transparent)
    p = QPainter(pm)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(QColor(colour))
    p.drawEllipse(QRectF(3, 3, 6, 6))
    p.end()
    return QIcon(pm)


COL_HANDLE, COL_NAME, COL_LEN, COL_P, COL_Q, COL_CI, COL_TYPE, COL_COPIES = range(8)
SHOW_KEY = "cytogenetics/karyotype_show"


class KaryotypeDialog(QDialog):
    """Create or edit one karyotype. Changes are written on Save; Cmd+Z undoes inside the window."""

    def __init__(self, folder: Path, name: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Edit karyotype")
        self.setMinimumSize(1220, 720)
        self._folder = Path(folder).expanduser().resolve()
        existing = kt.load_karyotype(self._folder, name) if name else None
        self._old_name = existing.name if existing else ""
        self._k: kt.Karyotype = copy.deepcopy(existing) if existing else kt.new_karyotype("New karyotype", 8, 2)
        self._listed = [c.name for c in self._k.chromosomes]
        self._saved = ""
        self._seq_names: list[str] = []
        self._undo: list[dict] = []
        self._redo: list[dict] = []
        self._filling = False
        self._lists = load_lists(self._folder)

        # ---- top row
        self._name = QLineEdit(self._k.name, self)
        self._name.setFixedWidth(120)
        self._species = NoWheelComboBox(self)
        self._species.setEditable(True)
        self._species.addItems([""] + lists_mod.names(self._lists, "species"))
        self._species.setCurrentText(self._k.species)
        self._species.setFixedWidth(160)
        self._species.setToolTip("Choose a species or type a new one: it is added to the species list on Save")
        self._basic = QSpinBox(self)
        self._basic.setRange(1, 60)
        self._basic.setValue(max(1, len(self._k.a_chromosomes()) or self._k.basic_number))
        self._basic.setFixedWidth(60)
        self._ploidy = NoWheelComboBox(self)
        for n in PLOIDIES:
            self._ploidy.addItem(f"{n}x", n)
        self._ploidy.setCurrentIndex(max(0, self._ploidy.findData(self._k.ploidy)))
        self._ploidy.setFixedWidth(60)
        self._formula = QLabel(self)
        self._formula.setStyleSheet("QLabel { background:#E3EBF8; color:#1F4E9A; font-weight:600; "
                                    "border-radius:10px; padding:3px 10px; }")
        self._hom_names = NoWheelComboBox(self)
        for key, (text, _stored) in kt.HOMOLOGUE_STYLES.items():
            self._hom_names.addItem(text, key)
        self._hom_names.addItem(OWN_LIST, "custom")
        self._hom_names.setFixedWidth(100)
        self._unit = NoWheelComboBox(self)
        for text, value in UNITS:
            self._unit.addItem(text, value)
        self._unit.setCurrentIndex(max(0, self._unit.findData(self._k.unit)))
        self._unit.setFixedWidth(100)
        self._basic_timer = QTimer(self)
        self._basic_timer.setSingleShot(True)
        self._basic_timer.setInterval(500)
        head = QHBoxLayout()
        head.setSpacing(6)
        for i, (text, w) in enumerate((("Name", self._name), ("Species", self._species), ("Basic number x", self._basic),
                                       ("Ploidy", self._ploidy))):
            if i:
                head.addSpacing(8)
            head.addWidget(QLabel(text, self))
            head.addWidget(w)
        head.addSpacing(10)
        head.addWidget(self._formula)
        head.addSpacing(10)
        for text, w in (("Homologues", self._hom_names), ("Unit", self._unit)):
            head.addWidget(QLabel(text, self))
            head.addWidget(w)
            head.addSpacing(8)
        head.addStretch(1)
        head.addWidget(make_help_button(self, "Karyotype", HELP))

        # ---- table
        self._table = _ReorderTable(0, 8, self)
        self._table.setHorizontalHeaderLabels(["", "Name", "Length", "p", "q", "CI", "Type", "Copies"])
        self._table.verticalHeader().setVisible(False)
        self._table.setShowGrid(False)
        self._table.setWordWrap(False)
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        hh.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        hh.setStretchLastSection(False)
        for col, w in ((COL_HANDLE, 24), (COL_NAME, 56), (COL_LEN, 82), (COL_P, 50), (COL_Q, 50), (COL_CI, 50),
                       (COL_TYPE, 44), (COL_COPIES, 54)):
            self._table.setColumnWidth(col, w)
        self._table.setEditTriggers(QTableWidget.EditTrigger.DoubleClicked | QTableWidget.EditTrigger.EditKeyPressed
                                    | QTableWidget.EditTrigger.AnyKeyPressed)
        self._table.rowMoved.connect(self._move_row)
        self._table.itemChanged.connect(self._on_item_changed)
        self._table.currentCellChanged.connect(lambda r, _c, _pr, _pc: self._on_row(r))
        plus = QPushButton("+", self)
        minus = QPushButton("−", self)
        for b, tip in ((plus, "Add a chromosome"), (minus, "Remove the selected chromosome")):
            b.setFixedSize(30, 24)
            b.setToolTip(tip)
        plus.clicked.connect(self._add_chromosome)
        minus.clicked.connect(self._remove_chromosome)
        sort_btn = QPushButton("Sort", self)
        sort_menu = QMenu(sort_btn)
        for key, text in (("listed", "As Listed"), ("length", "Longest First"), ("ci", "By Centromere Index")):
            sort_menu.addAction(text, lambda k=key: self._order(k))
        sort_btn.setMenu(sort_menu)
        imp_btn = QPushButton("Import", self)
        imp_menu = QMenu(imp_btn)
        imp_menu.addAction("Measurements (TSV)…", self._import_measurements)
        imp_menu.addAction("Lengths from Assembly…", self._take_from_assembly)
        imp_btn.setMenu(imp_menu)
        below = QHBoxLayout()
        below.setSpacing(4)
        below.addWidget(plus)
        below.addWidget(minus)
        below.addStretch(1)
        below.addSpacing(4)
        below.addWidget(sort_btn)
        below.addWidget(imp_btn)
        left = QVBoxLayout()
        left.setSpacing(8)
        left.addWidget(rounded_box(self._table, self), 1)
        left.addLayout(below)
        left_host = QWidget(self)
        left_host.setLayout(left)
        left.setContentsMargins(0, 0, 0, 0)
        left_host.setFixedWidth(440)

        # ---- preview toolbar
        self._canvas = KaryotypeCanvas(self)
        self._canvas.karyotype = self._k
        self._canvas.probe_colours = self._probe_colours()
        self._canvas.every = str(QSettings().value(SHOW_KEY, "every")) != "one"
        self._canvas.selectionChanged.connect(self._on_canvas_selection)
        self._canvas.aboutToEdit.connect(self._checkpoint)
        self._canvas.edited.connect(self._on_canvas_edited)
        self._canvas.menuRequested.connect(self._canvas_menu)
        self._tool_group = QButtonGroup(self)
        tools = QHBoxLayout()
        tools.setSpacing(0)
        for key, text in (("select", "Select"), ("draw", "Draw band")):
            b = QPushButton(text, self)
            b.setCheckable(True)
            b.setChecked(key == "select")
            b.clicked.connect(lambda _c=False, k=key: self._set_tool(k))
            self._tool_group.addButton(b)
            tools.addWidget(b)
        bar = QHBoxLayout()
        bar.setSpacing(5)
        bar.addLayout(tools)
        bar.addSpacing(10)
        self._swatch_group = QButtonGroup(self)
        self._swatches: list[tuple[QToolButton, str, str]] = []
        for kind, colour, tip in self._swatch_kinds():
            sb = QToolButton(self)
            sb.setCheckable(True)
            sb.setAutoRaise(True)
            sb.setIconSize(QSize(18, 18))
            sb.setIcon(_swatch_icon(self._canvas.band_brush(kt.Band(kind, 0, 1, colour)), 18))
            sb.setToolTip(tip)
            sb.clicked.connect(lambda _c=False, k=kind, c=colour: self._on_swatch(k, c))
            self._swatch_group.addButton(sb)
            self._swatches.append((sb, kind, colour))
            bar.addWidget(sb)
        if self._swatches:
            self._swatches[0][0].setChecked(True)
            self._canvas.draw_kind, self._canvas.draw_colour = self._swatches[0][1], self._swatches[0][2]
        own = QToolButton(self)
        own.setText("+")
        own.setAutoRaise(True)
        own.setToolTip("Own colour for the selected band or the next band you draw")
        own.clicked.connect(self._own_colour)
        bar.addWidget(own)
        self._presets_btn = QPushButton("Presets", self)
        self._presets_menu = QMenu(self._presets_btn)
        self._presets_menu.aboutToShow.connect(self._fill_presets_menu)
        self._presets_btn.setMenu(self._presets_menu)
        self._presets_btn.setToolTip("Band sets saved for this species")
        bar.addSpacing(6)
        bar.addWidget(self._presets_btn)
        bar.addStretch(1)
        bar.addWidget(QLabel("Show", self))
        self._show = NoWheelComboBox(self)
        self._show.addItem("Every homologue", True)
        self._show.addItem("One per chromosome", False)
        self._show.setCurrentIndex(0 if self._canvas.every else 1)
        self._show.setFixedWidth(160)
        self._show.activated.connect(self._on_show)
        bar.addWidget(self._show)

        # ---- inspector
        self._title = QLabel(self)
        tf = self._title.font()
        tf.setBold(True)
        tf.setPointSizeF(tf.pointSizeF() + 1)
        self._title.setFont(tf)
        self._kind_note = QLabel(self)
        self._kind_note.setStyleSheet("color: palette(dark);")
        self._len_edit, self._p_edit, self._q_edit = (QLineEdit(self) for _ in range(3))
        for e in (self._len_edit, self._p_edit, self._q_edit):
            e.setFixedWidth(90)
        self._len_edit.editingFinished.connect(self._on_len_field)
        self._p_edit.editingFinished.connect(self._on_arm_fields)
        self._q_edit.editingFinished.connect(self._on_arm_fields)
        form = QGridLayout()
        form.setHorizontalSpacing(10)
        form.setVerticalSpacing(6)
        title_row = QHBoxLayout()
        title_row.setSpacing(8)
        title_row.addWidget(self._title)
        title_row.addWidget(self._kind_note)
        title_row.addStretch(1)
        form.addLayout(title_row, 0, 0, 1, 2)
        for row, (text, e) in enumerate((("Length", self._len_edit), ("Short arm p", self._p_edit),
                                         ("Long arm q", self._q_edit)), start=1):
            form.addWidget(QLabel(text, self), row, 0)
            form.addWidget(e, row, 1)
        form.setRowStretch(4, 1)
        self._same = QCheckBox("Same as the other homologues", self)
        self._same.setToolTip("Untick to give this homologue its own length, centromere and bands")
        self._same.toggled.connect(self._on_same)
        self._bands_box = QVBoxLayout()
        self._bands_box.setSpacing(2)
        bands_host = QWidget(self)
        bands_host.setLayout(self._bands_box)
        self._bands_box.setContentsMargins(0, 0, 0, 0)
        self._bands_scroll = QScrollArea(self)
        self._bands_scroll.setWidget(bands_host)
        self._bands_scroll.setWidgetResizable(True)
        self._bands_scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._bands_scroll.setStyleSheet("QScrollArea, QScrollArea > QWidget > QWidget { background: transparent; }")
        self._bands_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._copy_link = QPushButton(self)
        self._copy_link.setStyleSheet(LINK_CSS)
        self._copy_link.setCursor(Qt.CursorShape.PointingHandCursor)
        self._copy_link.clicked.connect(self._copy_bands)
        right = QVBoxLayout()
        right.setSpacing(6)
        right.addWidget(self._same)
        right.addWidget(self._bands_scroll, 1)
        link_row = QHBoxLayout()
        link_row.addWidget(self._copy_link)
        link_row.addStretch(1)
        right.addLayout(link_row)
        inspector = QHBoxLayout()
        inspector.setContentsMargins(12, 10, 12, 10)
        inspector.setSpacing(22)
        form_host = QWidget(self)
        form_host.setLayout(form)
        form.setContentsMargins(0, 0, 0, 0)
        form_host.setFixedWidth(210)
        inspector.addWidget(form_host)
        inspector.addLayout(right, 1)
        insp_host = QWidget(self)
        insp_host.setLayout(inspector)
        insp_host.setFixedHeight(150)

        right_col = QVBoxLayout()
        right_col.setSpacing(8)
        right_col.setContentsMargins(0, 0, 0, 0)
        right_col.addLayout(bar)
        right_col.addWidget(rounded_box(self._canvas, self), 1)
        right_col.addWidget(rounded_box(insp_host, self))

        body = QHBoxLayout()
        body.setSpacing(14)
        body.addWidget(left_host)
        body.addLayout(right_col, 1)

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
        root.setContentsMargins(18, 14, 18, 14)
        root.setSpacing(10)
        root.addLayout(head)
        root.addLayout(body, 1)
        root.addLayout(buttons)

        # ---- signals
        self._name.textEdited.connect(self._on_name)
        self._species.editTextChanged.connect(lambda t: setattr(self._k, "species", t.strip()))
        self._basic.valueChanged.connect(lambda _v: self._basic_timer.start())
        self._basic.editingFinished.connect(self._on_basic)
        self._basic_timer.timeout.connect(self._on_basic)
        self._ploidy.activated.connect(self._on_ploidy)
        self._hom_names.activated.connect(self._on_hom_names)
        self._unit.activated.connect(self._on_unit)

        self._fill_table()
        self._sync_head()
        if self._k.chromosomes:
            self._table.setCurrentCell(0, COL_NAME)
        self._sync_inspector()

    # ------------------------------------------------------------------ helpers

    def saved_name(self) -> str:
        return self._saved

    def _probe_colours(self) -> dict[str, str]:
        out = {}
        for e in self._lists.get("probes", []):
            dye = str(e.get("dye") or "")
            if dye:
                out[str(e["name"])] = lists_mod.dye_color(self._lists, dye, "#2FB457")
        return out

    def _swatch_kinds(self) -> list[tuple[str, str, str]]:
        kinds = [("Secondary constriction (NOR)", "", "NOR (secondary constriction), hatched"),
                 ("DAPI+ heterochromatin", "", "DAPI+ heterochromatin"),
                 ("CMA+ band", "", "CMA+ band")]
        for probe in list(self._canvas.probe_colours)[:4]:
            kinds.append((kt.SIGNAL_KIND_PREFIX + probe, "", f"FISH signal of {probe}"))
        return kinds

    def _sel(self) -> tuple[kt.KaryoChromosome | None, int]:
        c = self._canvas
        if 0 <= c.ci < len(self._k.chromosomes):
            return self._k.chromosomes[c.ci], c.hi
        return None, -1

    def _refresh(self) -> None:
        self._fill_table()
        self._sync_head()
        self._canvas.refresh()
        self._sync_inspector()

    # ------------------------------------------------------------------ undo

    def _checkpoint(self) -> None:
        self._undo.append(self._k.to_dict())
        del self._undo[:-100]
        self._redo.clear()

    def _restore(self, state: dict) -> None:
        sel = (self._canvas.ci, self._canvas.hi, self._canvas.bi)
        self._k = kt.Karyotype.from_dict(state)
        self._canvas.karyotype = self._k
        n = len(self._k.chromosomes)
        self._canvas.ci, self._canvas.hi, self._canvas.bi = (sel if sel[0] < n else (-1, -1, -1))
        for w, value in ((self._basic, len(self._k.a_chromosomes())),):
            w.blockSignals(True)
            w.setValue(max(1, value))
            w.blockSignals(False)
        self._ploidy.setCurrentIndex(max(0, self._ploidy.findData(self._k.ploidy)))
        self._unit.setCurrentIndex(max(0, self._unit.findData(self._k.unit)))
        self._refresh()

    def _undo_step(self) -> None:
        if self._undo:
            self._redo.append(self._k.to_dict())
            self._restore(self._undo.pop())

    def _redo_step(self) -> None:
        if self._redo:
            self._undo.append(self._k.to_dict())
            self._restore(self._redo.pop())

    def _text_focus(self) -> bool:
        w = QApplication.focusWidget()
        return isinstance(w, (QLineEdit, QSpinBox)) or (w is not None and w.parent() is self._table.viewport())

    def event(self, e) -> bool:  # noqa: D401
        if e.type() == QEvent.Type.ShortcutOverride and not self._text_focus():
            if e.matches(QKeySequence.StandardKey.Undo) or e.matches(QKeySequence.StandardKey.Redo):
                e.accept()
                return True
        return super().event(e)

    def keyPressEvent(self, e) -> None:  # noqa: N802
        if not self._text_focus():
            if e.matches(QKeySequence.StandardKey.Undo):
                self._undo_step()
                return
            if e.matches(QKeySequence.StandardKey.Redo):
                self._redo_step()
                return
        super().keyPressEvent(e)

    # ------------------------------------------------------------------ top row

    def _sync_head(self) -> None:
        k = self._k
        self._formula.setText(k.formula)
        style = k.homologue_style
        i = self._hom_names.findData(style)
        if style == "custom":
            text = ", ".join(k.homologue_suffixes())
            self._hom_names.setItemText(i, text if len(text) <= 14 else text[:13] + "…")
        else:
            self._hom_names.setItemText(self._hom_names.findData("custom"), OWN_LIST)
        self._hom_names.setCurrentIndex(i)
        example = ", ".join(k.homologue_labels(k.chromosomes[0])) if k.chromosomes else ""
        self._hom_names.setToolTip(f"How the homologues are named: {example}" if example else "")
        header = self._table.horizontalHeaderItem(COL_LEN)
        if header is not None:
            header.setText(f"Length ({k.unit})")

    def _on_name(self, text: str) -> None:
        self._k.name = text.strip()

    def _on_hom_names(self, index: int) -> None:
        key = self._hom_names.itemData(index)
        self._checkpoint()
        if key == "custom":
            current = ", ".join(self._k.homologue_suffixes())
            text, ok = QInputDialog.getText(self, "Homologue names",
                                            "Names added to the chromosome name, separated by commas\n"
                                            "(for example: a, b, c or -I, -II, -III):", text=current)
            names = [t.strip() for t in text.split(",") if t.strip()] if ok else []
            if len(names) != len(set(names)):
                QMessageBox.information(self, "Homologue names", "Every homologue needs its own name.")
                names = []
            if names:
                self._k.homologue_letters = ",".join(names) + ("," if len(names) == 1 else "")
        else:
            self._k.homologue_letters = kt.HOMOLOGUE_STYLES[key][1]
        self._sync_head()
        self._canvas.refresh()
        self._sync_inspector()

    def _on_show(self, index: int) -> None:
        self._canvas.every = bool(self._show.itemData(index))
        QSettings().setValue(SHOW_KEY, "every" if self._canvas.every else "one")
        self._canvas.hi, self._canvas.bi = -1, -1
        self._canvas.refresh()
        self._sync_inspector()

    def _on_basic(self) -> None:
        self._basic_timer.stop()
        n = self._basic.value()
        k = self._k
        a = k.a_chromosomes()
        if n == len(a):
            return
        dropped = a[n:]
        if any(c.length > 0 or c.bands for c in dropped):
            answer = QMessageBox.question(self, "Basic number",
                                          f"Remove {len(dropped)} chromosome(s) with measurements from the end of the list?")
            if answer != QMessageBox.StandardButton.Yes:
                self._basic.blockSignals(True)
                self._basic.setValue(len(a))
                self._basic.blockSignals(False)
                return
        self._checkpoint()
        k.set_basic_number(n)
        self._listed = [c.name for c in k.chromosomes]
        self._canvas.ci = min(self._canvas.ci, len(k.chromosomes) - 1)
        self._refresh()

    def _on_ploidy(self, index: int) -> None:
        self._checkpoint()
        self._k.ploidy = int(self._ploidy.itemData(index))
        for c in self._k.chromosomes:
            kt.drop_extra_homologues(self._k, c)
        self._refresh()

    def _on_unit(self, index: int) -> None:
        self._k.unit = str(self._unit.itemData(index))
        self._sync_head()
        self._canvas.refresh()

    # ------------------------------------------------------------------ table

    @staticmethod
    def _fmt(v: float) -> str:
        return f"{v:g}" if v else ""

    def _fill_table(self) -> None:
        self._filling = True
        k = self._k
        t = self._table
        current = self._canvas.ci
        t.setRowCount(len(k.chromosomes))
        for r, c in enumerate(k.chromosomes):
            p, q = c.arms
            values = ["⋮⋮", c.name, self._fmt(c.length), self._fmt(round(p, 2)), self._fmt(round(q, 2)),
                      self._fmt(c.centromere_index) if c.length else "", c.type if c.length else "",
                      str(k.copies_of(c))]
            for col, text in enumerate(values):
                item = QTableWidgetItem(text)
                item.setTextAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
                flags = Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsDragEnabled
                if col not in (COL_HANDLE, COL_TYPE):
                    flags |= Qt.ItemFlag.ItemIsEditable
                item.setFlags(flags)
                if col == COL_HANDLE:
                    item.setForeground(QColor("#9A9AA0"))
                    item.setToolTip("Drag to change the order")
                if col == COL_NAME and kt.differs(c):
                    item.setIcon(_dot_icon())
                    item.setToolTip("Its homologues differ")
                if col == COL_COPIES and c.copies is None:
                    item.setToolTip("Same as the ploidy")
                t.setItem(r, col, item)
        if 0 <= current < t.rowCount():
            t.setCurrentCell(current, max(1, t.currentColumn()))
        self._filling = False

    def _on_row(self, r: int) -> None:
        if self._filling or r < 0 or r == self._canvas.ci:
            return
        self._canvas.ci, self._canvas.hi, self._canvas.bi = r, -1, -1
        self._canvas.update()
        self._sync_inspector()

    def _num(self, text: str) -> float | None:
        try:
            return float(text.replace(",", ".").strip())
        except ValueError:
            return None

    def _on_item_changed(self, item: QTableWidgetItem) -> None:
        if self._filling:
            return
        r, col = item.row(), item.column()
        if not 0 <= r < len(self._k.chromosomes):
            return
        c = self._k.chromosomes[r]
        text = item.text().strip()
        self._checkpoint()
        if col == COL_NAME:
            if text and text not in {x.name for x in self._k.chromosomes if x is not c}:
                self._listed = [text if n == c.name else n for n in self._listed]
                c.name = text
        elif col == COL_LEN:
            v = self._num(text)
            if v is not None and v >= 0:
                c.length = v
        elif col in (COL_P, COL_Q):
            p, q = c.arms
            v = self._num(text)
            if v is not None and v >= 0:
                c.set_arms(v if col == COL_P else p, v if col == COL_Q else q)
        elif col == COL_CI:
            v = self._num(text)
            if v is not None and 0 <= v <= 50:
                c.centromere_index = v
        elif col == COL_COPIES:
            v = self._num(text)
            if text == "":
                c.copies = None
            elif v is not None and v >= 1:
                c.copies = None if int(v) == self._k.ploidy else int(v)
            kt.drop_extra_homologues(self._k, c)
        QTimer.singleShot(0, self._refresh)

    def _move_row(self, src: int, dst: int) -> None:
        self._checkpoint()
        chroms = self._k.chromosomes
        chroms.insert(dst, chroms.pop(src))
        self._listed = [c.name for c in chroms]
        self._canvas.ci, self._canvas.hi, self._canvas.bi = dst, -1, -1
        self._refresh()

    def _order(self, key: str) -> None:
        self._checkpoint()
        k = self._k
        if key == "listed":
            rank = {name: i for i, name in enumerate(self._listed)}
            k.chromosomes.sort(key=lambda c: (c.b_chromosome, rank.get(c.name, len(rank))))
        else:
            k.chromosomes[:] = sorted(k.ordered(key), key=lambda c: c.b_chromosome)
        self._canvas.ci, self._canvas.hi, self._canvas.bi = (0 if k.chromosomes else -1), -1, -1
        self._refresh()

    def _add_chromosome(self) -> None:
        self._checkpoint()
        k = self._k
        used = {c.name for c in k.chromosomes}
        n = len(k.a_chromosomes()) + 1
        while str(n) in used:
            n += 1
        new = kt.KaryoChromosome(name=str(n))
        pos = len(k.a_chromosomes())
        k.chromosomes.insert(pos, new)
        self._listed.append(new.name)
        self._set_basic_silently()
        self._canvas.ci, self._canvas.hi, self._canvas.bi = pos, -1, -1
        self._refresh()

    def _remove_chromosome(self) -> None:
        r = self._canvas.ci if self._canvas.ci >= 0 else self._table.currentRow()
        if not 0 <= r < len(self._k.chromosomes):
            return
        c = self._k.chromosomes[r]
        if (c.length > 0 or c.bands) and QMessageBox.question(
                self, "Remove chromosome", f"Remove chromosome {c.name} with its measurements?") != QMessageBox.StandardButton.Yes:
            return
        self._checkpoint()
        del self._k.chromosomes[r]
        self._listed = [n for n in self._listed if n != c.name]
        self._set_basic_silently()
        self._canvas.ci = min(r, len(self._k.chromosomes) - 1)
        self._canvas.hi = self._canvas.bi = -1
        self._refresh()

    def _set_basic_silently(self) -> None:
        n = max(1, len(self._k.a_chromosomes()))
        self._k.basic_number = n
        self._basic.blockSignals(True)
        self._basic.setValue(n)
        self._basic.blockSignals(False)

    # ------------------------------------------------------------------ preview

    def _set_tool(self, tool: str) -> None:
        self._canvas.tool = tool
        self._canvas.setCursor(Qt.CursorShape.CrossCursor if tool == "draw" else Qt.CursorShape.ArrowCursor)

    def _on_swatch(self, kind: str, colour: str) -> None:
        self._canvas.draw_kind, self._canvas.draw_colour = kind, colour
        band = self._selected_band()
        if band is not None:
            self._checkpoint()
            band.kind, band.color = kind, colour
            band.satellite = band.satellite and kind.startswith("Secondary")
            self._after_band_change()

    def _own_colour(self) -> None:
        band = self._selected_band()
        start = QColor(band.color) if band is not None and band.color else QColor("#E8A33D")
        colour = QColorDialog.getColor(start, self, "Band colour")
        if not colour.isValid():
            return
        self._canvas.draw_colour = colour.name()
        if band is not None:
            self._checkpoint()
            band.color = colour.name()
            self._after_band_change()

    def _selected_band(self) -> kt.Band | None:
        c, hi = self._sel()
        if c is None:
            return None
        bands = kt.target(c, hi).bands
        return bands[self._canvas.bi] if 0 <= self._canvas.bi < len(bands) else None

    def _after_band_change(self) -> None:
        self._canvas.update()
        self._sync_inspector()

    def _on_canvas_selection(self) -> None:
        ci = self._canvas.ci
        if 0 <= ci < self._table.rowCount() and self._table.currentRow() != ci:
            self._filling = True
            self._table.setCurrentCell(ci, max(1, self._table.currentColumn()))
            self._filling = False
        self._sync_inspector()

    def _on_canvas_edited(self) -> None:
        r = self._canvas.ci
        self._fill_table()
        self._sync_head()
        self._sync_inspector(fields_only=self._canvas._drag is not None)
        if 0 <= r < self._table.rowCount():
            self._table.setCurrentCell(r, max(1, self._table.currentColumn()))

    def _canvas_menu(self, global_pos) -> None:
        c, hi = self._sel()
        if c is None:
            return
        menu = QMenu(self)
        band = self._selected_band()
        if band is not None:
            menu.addAction("Colour…", self._own_colour)
            if band.kind.startswith("Secondary"):
                a = menu.addAction("Draw as Satellite")
                a.setCheckable(True)
                a.setChecked(band.satellite)
                a.toggled.connect(lambda on: self._band_flag("satellite", on))
            menu.addAction("Remove Band", self._remove_band)
            menu.addSeparator()
        if self._canvas.every and hi >= 0 and self._k.copies_of(c) > 1:
            a = menu.addAction("Same as the Other Homologues")
            a.setCheckable(True)
            a.setChecked(not kt.is_separate(c, hi))
            a.toggled.connect(lambda on: self._same.setChecked(on))
            menu.addAction("Copy Bands to the Other Homologues", self._copy_bands)
        menu.addAction("Mirror Bands", self._mirror)
        menu.addSeparator()
        menu.addAction("Duplicate as B Chromosome", self._duplicate_b)
        if self._canvas.every and hi >= 0:
            a = menu.addAction("Mark as Missing")
            a.setCheckable(True)
            a.setChecked(kt.is_missing(c, hi))
            a.toggled.connect(self._set_missing)
        menu.exec(global_pos)

    def _band_flag(self, attr: str, on: bool) -> None:
        band = self._selected_band()
        if band is not None:
            self._checkpoint()
            setattr(band, attr, on)
            self._after_band_change()

    def _remove_band(self, index: int | None = None) -> None:
        c, hi = self._sel()
        if c is None:
            return
        i = self._canvas.bi if index is None else index
        bands = kt.target(c, hi).bands
        if 0 <= i < len(bands):
            self._checkpoint()
            del bands[i]
            self._canvas.bi = -1
            self._canvas.update()
            self._sync_inspector()

    def _mirror(self) -> None:
        c, hi = self._sel()
        if c is not None:
            self._checkpoint()
            kt.mirror_bands(c, hi if self._canvas.every else -1)
            self._refresh()

    def _copy_bands(self) -> None:
        c, hi = self._sel()
        if c is not None and hi >= 0:
            self._checkpoint()
            kt.copy_bands_to_others(self._k, c, hi)
            self._refresh()

    def _duplicate_b(self) -> None:
        c, hi = self._sel()
        if c is None:
            return
        self._checkpoint()
        new = kt.add_b_chromosome(self._k, c, hi)
        self._listed.append(new.name)
        self._canvas.ci, self._canvas.hi, self._canvas.bi = len(self._k.chromosomes) - 1, -1, -1
        self._refresh()

    def _set_missing(self, on: bool) -> None:
        c, hi = self._sel()
        if c is not None and hi >= 0:
            self._checkpoint()
            kt.set_missing(c, hi, on)
            self._refresh()

    def _fill_presets_menu(self) -> None:
        m = self._presets_menu
        m.clear()
        species = self._k.species or "Any species"
        m.addAction("Save Bands as Preset…", self._save_preset)
        m.addSeparator()
        names = kt.list_presets(self._folder, self._k.species)
        if not names:
            a = m.addAction(f"No presets for {species}")
            a.setEnabled(False)
        for n in names:
            m.addAction(n, lambda name=n: self._apply_preset(name))

    def _save_preset(self) -> None:
        if not any(c.bands for c in self._k.chromosomes):
            QMessageBox.information(self, "Presets", "There are no bands to save yet.")
            return
        name, ok = QInputDialog.getText(self, "Save bands as preset",
                                        f"Name of the band set for {self._k.species or 'any species'}:")
        if ok and name.strip():
            if name.strip() in kt.list_presets(self._folder, self._k.species) and QMessageBox.question(
                    self, "Presets", f"\"{name.strip()}\" exists. Replace it?") != QMessageBox.StandardButton.Yes:
                return
            kt.save_preset(self._folder, self._k.species, name.strip(), self._k)

    def _apply_preset(self, name: str) -> None:
        self._checkpoint()
        missing = kt.apply_preset(self._folder, self._k.species, name, self._k)
        self._refresh()
        if missing:
            QMessageBox.information(self, "Presets", "Not in this karyotype: " + ", ".join(missing))

    # ------------------------------------------------------------------ inspector

    def _sync_inspector(self, fields_only: bool = False) -> None:
        c, hi = self._sel()
        enabled = c is not None
        for w in (self._len_edit, self._p_edit, self._q_edit):
            w.setEnabled(enabled)
        if c is None:
            self._title.setText("")
            self._kind_note.setText("Click a chromosome")
            for e in (self._len_edit, self._p_edit, self._q_edit):
                e.clear()
            self._same.setVisible(False)
            self._copy_link.setVisible(False)
            self._fill_bands([], None)
            return
        k = self._k
        every_hom = self._canvas.every and hi >= 0
        view = kt.homologue_view(c, hi if every_hom else -1)
        self._title.setText(k.homologue_label(c, hi) if every_hom else c.name)
        self._kind_note.setText(f"{view.type}, CI {view.centromere_index:g}" if view.length else "no length yet")
        p, q = view.arms
        for e, v in ((self._len_edit, view.length), (self._p_edit, round(p, 2)), (self._q_edit, round(q, 2))):
            if not e.hasFocus():
                e.setText(self._fmt(v))
        if fields_only:
            return
        show_same = every_hom and k.copies_of(c) > 1 and not kt.is_missing(c, hi)
        self._same.setVisible(show_same)
        self._same.blockSignals(True)
        self._same.setChecked(not kt.is_separate(c, hi) if show_same else True)
        self._same.blockSignals(False)
        others = [k.homologue_label(c, i) for i in range(k.copies_of(c)) if i != hi]
        self._copy_link.setText("Copy bands to " + ", ".join(others))
        self._copy_link.setVisible(show_same and kt.is_separate(c, hi))
        self._fill_bands(view.bands, view)

    def _fill_bands(self, bands: list[kt.Band], view: kt.KaryoChromosome | None) -> None:
        while self._bands_box.count():
            item = self._bands_box.takeAt(0)
            if item.widget() is not None:
                item.widget().deleteLater()
        if view is None:
            self._bands_box.addStretch(1)
            return
        if not bands:
            hint = QLabel("No bands. Choose Draw band and drag along a chromosome.", self)
            hint.setStyleSheet("color: palette(dark);")
            hint.setWordWrap(True)
            self._bands_box.addWidget(hint)
        cf = view.centromere_fraction
        for i, b in enumerate(bands):
            row = QFrame(self)
            row.setObjectName("bandRow")
            if i == self._canvas.bi:
                row.setStyleSheet("QFrame#bandRow { background:#E3EBF8; border-radius:6px; }")
            lay = QHBoxLayout(row)
            lay.setContentsMargins(6, 2, 2, 2)
            lay.setSpacing(8)
            sw = QLabel(row)
            sw.setPixmap(_swatch_icon(self._canvas.band_brush(b), 14).pixmap(14, 14))
            name = QLabel(kt.manual_probe(b.kind) and f"FISH {kt.manual_probe(b.kind)}" or b.kind, row)
            name.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
            rng = QLineEdit(f"{kt.fraction_to_arm_text(b.start, cf)} to {kt.fraction_to_arm_text(b.end, cf)}", row)
            rng.setFixedWidth(120)
            rng.setFrame(False)
            rng.setToolTip("Where the band sits, for example q 0.41 to 0.59 (0 = centromere, 1 = telomere)")
            rng.editingFinished.connect(lambda e=rng, n=i: self._on_band_range(n, e.text()))
            x = QToolButton(row)
            x.setText("×")
            x.setAutoRaise(True)
            x.setToolTip("Remove band")
            x.clicked.connect(lambda _c=False, n=i: self._remove_band(n))
            lay.addWidget(sw)
            lay.addWidget(name, 1)
            lay.addWidget(rng)
            lay.addWidget(x)
            row.mousePressEvent = lambda _e, n=i: self._pick_band(n)  # noqa: B023
            self._bands_box.addWidget(row)
        self._bands_box.addStretch(1)

    def _pick_band(self, i: int) -> None:
        self._canvas.bi = i
        self._canvas.update()
        self._sync_inspector()

    def _on_band_range(self, i: int, text: str) -> None:
        c, hi = self._sel()
        if c is None:
            return
        owner = kt.target(c, hi if self._canvas.every else -1)
        if not 0 <= i < len(owner.bands):
            return
        cf = kt.homologue_view(c, hi if self._canvas.every else -1).centromere_fraction
        parts = [s.strip() for s in text.replace("–", " to ").split(" to ")]
        if len(parts) != 2:
            return
        a, z = (kt.arm_text_to_fraction(s, cf) for s in parts)
        if a is None or z is None:
            return
        b = owner.bands[i]
        if (round(b.start, 4), round(b.end, 4)) == tuple(sorted((round(a, 4), round(z, 4)))):
            return
        self._checkpoint()
        b.start, b.end = sorted((a, z))
        self._canvas.update()

    def _owner_for_fields(self):
        c, hi = self._sel()
        if c is None:
            return None, None
        return c, kt.target(c, hi if self._canvas.every else -1)

    def _on_len_field(self) -> None:
        c, owner = self._owner_for_fields()
        v = self._num(self._len_edit.text())
        if owner is None or v is None or v < 0 or v == (owner.length or 0):
            return
        self._checkpoint()
        owner.length = v
        self._refresh()

    def _on_arm_fields(self) -> None:
        c, owner = self._owner_for_fields()
        p, q = self._num(self._p_edit.text()), self._num(self._q_edit.text())
        if owner is None or p is None or q is None or p + q <= 0:
            return
        p, q = sorted((max(0.0, p), max(0.0, q)))
        length, ci = round(p + q, 4), round(100.0 * p / (p + q), 2)
        if (length, ci) == (round(owner.length or 0, 4), round(owner.centromere_index or 0, 2)):
            return
        self._checkpoint()
        owner.length, owner.centromere_index = length, ci
        self._refresh()

    def _on_same(self, on: bool) -> None:
        c, hi = self._sel()
        if c is None or hi < 0:
            return
        if on and kt.is_separate(c, hi) and QMessageBox.question(
                self, "Homologues", f"Make {self._k.homologue_label(c, hi)} the same as the others again? "
                                    "Its own length, centromere and bands are dropped.") != QMessageBox.StandardButton.Yes:
            self._same.blockSignals(True)
            self._same.setChecked(False)
            self._same.blockSignals(False)
            return
        self._checkpoint()
        kt.set_separate(c, hi, not on)
        self._refresh()

    # ------------------------------------------------------------------ import

    def _import_measurements(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Import measurements", str(self._folder), "Tables (*.tsv *.csv *.txt);;All files (*)")
        if not path:
            return
        try:
            rows = kt.parse_measurement_table(Path(path).read_text(encoding="utf-8", errors="replace"))
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Import measurements", str(exc))
            return
        if not rows:
            QMessageBox.information(self, "Import measurements",
                                    "No chromosome rows found. The table needs a chromosome column and a length "
                                    "(or p and q) column.")
            return
        self._checkpoint()
        unit = self._k.unit if self._k.unit in ("µm", "%") else "µm"
        missing = kt.apply_measurements(self._k, rows, unit)
        self._unit.setCurrentIndex(max(0, self._unit.findData(unit)))
        self._refresh()
        text = f"{len(rows) - len(missing)} chromosome(s) updated."
        if missing:
            text += "\nNot in this karyotype: " + ", ".join(missing)
        QMessageBox.information(self, "Import measurements", text)

    def _take_from_assembly(self) -> None:
        dialog = AssemblyLengthsDialog(self._folder, self._k, self)
        dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        dialog.accepted.connect(lambda d=dialog: self._apply_assembly(d))
        dialog.show()  # not modal: a click on a sizes file or FASTA in the Files panel fills Source

    def _apply_assembly(self, dialog: "AssemblyLengthsDialog") -> None:
        self._checkpoint()
        self._seq_names = dialog.sequence_names()
        missing = kt.apply_sizes(self._k, dialog.links(), dialog.sizes())
        self._unit.setCurrentIndex(max(0, self._unit.findData("Mb")))
        self._refresh()
        if missing:
            QMessageBox.information(self, "Lengths from assembly",
                                    "No sequence chosen for chromosome " + ", ".join(missing) + ". Their lengths are unchanged.")

    # ------------------------------------------------------------------ save

    def _save(self) -> None:
        k = self._k
        k.name = self._name.text().strip()
        k.species = self._species.currentText().strip()
        if not k.name:
            QMessageBox.information(self, "Karyotype", "Give the karyotype a name.")
            return
        names = [c.name for c in k.chromosomes]
        if len(set(names)) != len(names) or not all(names):
            QMessageBox.information(self, "Karyotype", "Every chromosome needs its own name.")
            return
        if k.name != self._old_name and kt.karyotype_path(self._folder, k.name).exists():
            answer = QMessageBox.question(self, "Karyotype", f"A karyotype named \"{k.name}\" exists. Replace it?")
            if answer != QMessageBox.StandardButton.Yes:
                return
        k.basic_number = max(1, len(k.a_chromosomes()))
        try:
            if k.species and k.species.lower() not in {n.lower() for n in lists_mod.names(self._lists, "species")}:
                lists = load_lists(self._folder)
                lists_mod.add_entry(lists, "species", k.species)
                lists_mod.save_lists(self._folder, lists)
            kt.save_karyotype(self._folder, k, old_name=self._old_name)
        except OSError as exc:
            QMessageBox.warning(self, "Karyotype", f"Could not save:\n{exc}")
            return
        from repeatviewer.widgets.cyto_images_ui import notify_changed

        self._saved = k.name
        notify_changed(self._folder)
        self.accept()


class _SizesReader(QObject):
    """Reads sequence lengths on a worker thread (a genome FASTA takes a while)."""

    progress = Signal(int)
    done = Signal(dict, str)

    def __init__(self, path: Path) -> None:
        super().__init__()
        self.path = path
        self.stopped = False

    def run(self) -> None:
        try:
            sizes = kt.read_sizes(self.path, self.progress.emit, lambda: self.stopped)
            if sizes and kt.is_fasta(self.path) and not self.stopped:
                side = kt.sizes_path_for(self.path)
                if not side.exists() or side.stat().st_mtime < self.path.stat().st_mtime:
                    try:
                        kt.write_sizes(side, sizes)
                    except OSError:
                        pass
            self.done.emit(sizes, "")
        except Exception as exc:  # noqa: BLE001
            self.done.emit({}, str(exc))


class AssemblyLengthsDialog(PanelPathFields, QDialog):
    """Pick a sizes file, index, BED, GFF or FASTA and match its sequences to the chromosomes."""

    def __init__(self, folder: Path, karyotype: kt.Karyotype, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        from repeatviewer.widgets.path_browse import PathBrowseRow

        self.setWindowTitle("Lengths from assembly")
        self.setMinimumSize(620, 520)
        self._folder = folder
        self._k = karyotype
        self._sizes: dict[str, int] = {}
        self._thread: QThread | None = None
        self._reader: _SizesReader | None = None
        self._loaded: Path | None = None
        sources = kt.sizes_sources(folder)
        self._path = PathBrowseRow(self, initial=sources[0] if sources else None, browse_title="Lengths from assembly",
                                   file_filter=kt.SIZES_FILTER, start_dir=str(folder),
                                   placeholder="Sizes file, .fai, BED, GFF or FASTA", elide_middle=True)
        self._path.path_changed.connect(self._load)
        self.watch_path_fields(self._path)
        self._status = QLabel(self)
        self._status.setStyleSheet("color: gray;")
        self._bar = QProgressBar(self)
        self._bar.setRange(0, 100)
        self._bar.setTextVisible(False)
        self._bar.setFixedHeight(6)
        self._bar.hide()
        self._table = QTableWidget(0, 3, self)
        self._table.setHorizontalHeaderLabels(["Chromosome", "Sequence", "Length (Mb)"])
        self._table.verticalHeader().setVisible(False)
        self._table.verticalHeader().setDefaultSectionSize(30)
        self._table.setShowGrid(False)
        self._table.setAlternatingRowColors(True)
        self._table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)

        source = QHBoxLayout()
        source.addWidget(QLabel("Source", self))
        source.addWidget(self._path, 1)
        source.addWidget(make_help_button(self, "Lengths from assembly", HELP_ASSEMBLY))
        cancel = QPushButton("Cancel", self)
        self._save = QPushButton("Save", self)
        self._save.setDefault(True)
        cancel.clicked.connect(self.reject)
        self._save.clicked.connect(self.accept)
        buttons = QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(cancel)
        buttons.addWidget(self._save)
        root = QVBoxLayout(self)
        root.setContentsMargins(20, 20, 20, 16)
        root.setSpacing(10)
        root.addLayout(source)
        root.addWidget(self._status)
        root.addWidget(self._bar)
        root.addWidget(rounded_box(self._table, self), 1)
        root.addLayout(buttons)
        self._fill()
        QTimer.singleShot(0, self._load)

    # ---------------------------------------------------------------- reading

    def _load(self) -> None:
        path = self._path.path()
        if path is not None and path == self._loaded:
            return
        self._loaded = path
        self._stop()
        self._sizes = {}
        if not path or not Path(path).is_file():
            self._status.setText("Choose a sizes file, .fai, BED, GFF or FASTA.")
            self._fill()
            return
        path = Path(path)
        self._status.setText(f"Reading {path.name}…")
        self._bar.setValue(0)
        self._bar.show()
        self._save.setEnabled(False)
        self._thread = QThread(self)
        self._reader = _SizesReader(path)
        self._reader.moveToThread(self._thread)
        self._thread.started.connect(self._reader.run)
        self._reader.progress.connect(self._bar.setValue, Qt.ConnectionType.QueuedConnection)
        self._reader.done.connect(self._on_read, Qt.ConnectionType.QueuedConnection)
        self._reader.done.connect(self._thread.quit)
        self._thread.start()

    def _on_read(self, sizes: dict, error: str) -> None:
        self._bar.hide()
        self._save.setEnabled(True)
        path = Path(self._path.path())
        if error or not sizes:
            self._status.setText(f"Could not read lengths from {path.name}. {error}".strip())
        else:
            note = ""
            if kt.is_fasta(path) and kt.sizes_path_for(path).exists():
                note = f" Saved as {kt.sizes_path_for(path).name} beside it."
            self._status.setText(f"{len(sizes)} sequences in {path.name}.{note}")
        self._sizes = dict(sizes)
        self._fill()

    def _stop(self) -> None:
        if self._reader is not None:
            self._reader.stopped = True
        if self._thread is not None:
            self._thread.quit()
            self._thread.wait(3000)
        self._thread = None
        self._reader = None

    def done(self, result: int) -> None:  # noqa: D401
        self._stop()
        super().done(result)

    # ---------------------------------------------------------------- matching

    def _fill(self) -> None:
        names = self.sequence_names()
        links = kt.match_sequences(self._k, names)
        self._table.setRowCount(len(self._k.chromosomes))
        for r, c in enumerate(self._k.chromosomes):
            item = QTableWidgetItem(c.name)
            item.setFlags(Qt.ItemFlag.ItemIsEnabled)
            self._table.setItem(r, 0, item)
            box = NoWheelComboBox(self)
            box.addItem("")
            box.addItems(names)
            box.setCurrentText(links.get(c.name, ""))
            box.currentTextChanged.connect(lambda _t, r=r: self._show_length(r))
            holder = QWidget(self)
            lay = QHBoxLayout(holder)
            lay.setContentsMargins(4, 2, 4, 2)
            lay.addWidget(box)
            self._table.setCellWidget(r, 1, holder)
            length = QTableWidgetItem()
            length.setFlags(Qt.ItemFlag.ItemIsEnabled)
            length.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            self._table.setItem(r, 2, length)
            self._show_length(r)

    def _box(self, r: int) -> NoWheelComboBox | None:
        holder = self._table.cellWidget(r, 1)
        return holder.findChild(NoWheelComboBox) if holder is not None else None

    def _show_length(self, r: int) -> None:
        box, item = self._box(r), self._table.item(r, 2)
        if box is None or item is None:
            return
        bp = self._sizes.get(box.currentText())
        item.setText(f"{bp / 1e6:.3f}" if bp else "")

    # ---------------------------------------------------------------- results

    def sizes(self) -> dict[str, int]:
        return dict(self._sizes)

    def sequence_names(self) -> list[str]:
        return list(self._sizes)

    def links(self) -> dict[str, str]:
        out = {}
        for r, c in enumerate(self._k.chromosomes):
            box = self._box(r)
            if box is not None and box.currentText():
                out[c.name] = box.currentText()
        return out
