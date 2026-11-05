"""Idiogram tab of the Cytogenetics workspace: in situ (FISH) against in silico (assembly)."""

from __future__ import annotations

from pathlib import Path

import re

from PySide6.QtCore import QPoint, QPointF, QRectF, QSize, Qt, Signal
from PySide6.QtGui import QBrush, QColor, QFont, QIcon, QImage, QPainter, QPainterPath, QPen, QPixmap, QPolygonF
from PySide6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QCheckBox,
    QColorDialog,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMenu,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QStackedWidget,
    QToolButton,
    QToolTip,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import idiogram as idg
from repeatviewer.backends.cytogenetics import karyotype as kt
from repeatviewer.backends.cytogenetics import lists as lists_mod
from repeatviewer.backends.cytogenetics.model import image_details
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cytogenetics_style import clear_layout
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.views.cyto_karyotype_dialog import BAND_COLOURS, KaryoPreview, chromosome_path
from repeatviewer.widgets.sequence_annotation_rail import _CollapsibleSection

SIDE_W = 330
FALLBACK_COLOURS = ("#2FB457", "#E8405F", "#E8A33D", "#5B74FF", "#B05BD6", "#2BB5C4")
STATUS_COLOURS = {
    idg.STATUS_AGREES: ("#E6F4EA", "#1E6B35"),
    idg.STATUS_NONE: ("#F1F1F3", "#45484D"),
    idg.STATUS_FISH: ("#FDECEF", "#A3213B"),
    idg.STATUS_ASM: ("#FFF3DD", "#8A5300"),
    idg.STATUS_DIFFERS: ("#FDECEF", "#A3213B"),
    idg.STATUS_UNLINKED: ("#F1F1F3", "#45484D"),
}
LINK_CSS = "QPushButton { color:#1F4E9A; border:none; text-align:left; } QPushButton:hover { text-decoration: underline; }"
GENOME_NONE, GENOME_BROWSE = "(no genome)", "Browse…"

HELP_TEXT = (
    "<b>Idiogram</b><br>"
    "The karyotype drawn as in Edit karyotype, with what you measured in the Signals tab on top.<br>"
    "<b>Karyotype</b>: which karyotype to show; each one keeps its own view (probes, assembly, layout).<br>"
    "<b>Compare with</b>: a second karyotype drawn under the first, for example a related species.<br>"
    "<b>Show</b>: every homologue (1a, 1b, 1c) or one bar per chromosome.<br>"
    "<b>Layers</b>: ■ karyotype bands (drawn in Edit karyotype), ▶ measured signals (pooled over the marked cells; "
    "hover for how often each was seen), dashed bar = assembly (choose a genome under Assembly).<br>"
    "<b>Probes</b>: untick to hide a probe; <i>Seen in at least</i> hides rare signals.<br>"
    "<b>Measured, not in the karyotype</b>: signals with no band of that probe nearby. "
    "<b>Add to karyotype</b> makes them bands.<br>"
    "<b>Drag</b> a triangle to correct a signal (Layout > Drag signals). Click it to open the cells it comes from.<br>"
    "<b>Add to figure</b>: Colour, or Journal (black and white with a table)."
)


SHOW_ITEMS = (("Every homologue", True), ("One per chromosome", False))
PAINT_KINDS = (("NOR (hatched)", kt.BAND_KINDS[0]), ("DAPI+ (blue)", kt.BAND_KINDS[1]), ("CMA+ (yellow)", kt.BAND_KINDS[2]))
OWN_COLOUR = "Own colour…"
TOOLS = (("select", "Select"), ("paint", "Paint band"), ("place", "Place signal"))


def _section_label(text: str) -> QLabel:
    label = QLabel(text)
    label.setWordWrap(True)
    label.setStyleSheet("color: palette(dark); font-size: 12px;")
    return label


class IdiogramCanvas(QWidget):
    """Chromosomes drawn like the Edit karyotype window, grouped by chromosome (one bar per homologue or one per
    chromosome). Layers: karyotype bands (boxes), measured signals (triangles on the left), assembly (thin bar)."""

    bandMoved = Signal(object, float)  # FishBand, new fraction
    orderChanged = Signal(list)  # chromosome names
    bandClicked = Signal(int, int, QPoint)  # row, band, global position
    arrayClicked = Signal(int, int, QPoint)  # row, array, global position

    BAR_W, HOM_GAP, GROUP_GAP, LABEL_W, ASM_W, ASM_GAP = 26.0, 14.0, 46.0, 28.0, 9.0, 5.0
    TOP, BOTTOM = 30, 44

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.rows: list[idg.ChromosomeRow] = []
        self.heights: list[float] = []
        self.colours: dict[str, str] = {}
        self.karyotype: kt.Karyotype | None = None
        self.title = ""
        self.show_bands = True
        self.show_measured = True
        self.show_assembly = True
        self.assembly_active = False
        self.drag_enabled = True
        self.reorder_enabled = False
        self.fit = True  # shrink to the width (Figure); the Idiogram tab keeps full size and scrolls
        self._drag_band: tuple[int, int] | None = None
        self._drag_frac = 0.0
        self._drag_slot: int | None = None
        self._drag_x = 0.0
        self._press: tuple[QPointF, tuple] | None = None
        self._rects: list[tuple[QRectF, QRectF]] = []
        self._labels: list[tuple[str, QRectF]] = []
        self._size = (0, 0)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setMinimumHeight(420)

    def set_rows(self, rows: list[idg.ChromosomeRow], heights: list[float], colours: dict[str, str]) -> None:
        self.rows, self.heights, self.colours = rows, heights, colours
        self._size = (0, 0)
        n_groups = len(self._groups())
        f_min = 10.0 / self.BAR_W if self.fit else 1.0  # smallest scale used by _layout
        asm = (self.ASM_W + self.ASM_GAP) * f_min if self.show_assembly and self.assembly_active else 0.0
        smallest = (n_groups * self.LABEL_W + len(rows) * (self.BAR_W * f_min + asm)
                    + (len(rows) - n_groups) * max(4.0, self.HOM_GAP * f_min)
                    + max(0, n_groups - 1) * max(14.0, self.GROUP_GAP * f_min) + 50)
        self.setMinimumWidth(max(300, int(smallest) + 1))  # a scroll bar only when even the smallest bars do not fit
        self.update()

    def sizeHint(self) -> QSize:  # noqa: D102
        return QSize(max(600, len(self.rows) * 48), 560)

    # ---------------------------------------------------------------- geometry

    def _groups(self) -> list[list[int]]:
        groups: list[list[int]] = []
        for i, row in enumerate(self.rows):
            if groups and self.rows[groups[-1][0]].chromosome.name == row.chromosome.name:
                groups[-1].append(i)
            else:
                groups.append([i])
        return groups

    def _has_asm(self) -> bool:
        return self.show_assembly and self.assembly_active and any(r.sequence_length for r in self.rows)

    def _layout(self, width: float, height: float) -> None:
        self._rects, self._labels = [], []
        groups = self._groups()
        if not groups:
            return
        has_asm = self._has_asm()
        top = self.TOP + (22 if self.title else 0)
        bottom = self.BOTTOM + (34 if self.assembly_active and self.show_assembly else 0)
        lab = self.LABEL_W
        n_bars, n_groups = len(self.rows), len(groups)
        asm_full = (self.ASM_W + self.ASM_GAP) if has_asm else 0.0
        fixed = lab * n_groups  # chromosome names do not shrink
        scalable = n_bars * (self.BAR_W + asm_full) + (n_bars - n_groups) * self.HOM_GAP + (n_groups - 1) * self.GROUP_GAP
        f = min(1.0, max(0.05, (width - 40 - fixed) / scalable)) if scalable > 0 and self.fit else 1.0
        bar = max(10.0, self.BAR_W * f)
        gap = max(4.0, self.HOM_GAP * f)
        ggap = max(14.0, self.GROUP_GAP * f)
        asm_w = asm_full * f
        need = fixed + n_bars * (bar + asm_w) + (n_bars - n_groups) * gap + (n_groups - 1) * ggap
        usable = max(80.0, height - top - bottom)
        x = max(20.0, (width - need) / 2)
        self._rects = [(QRectF(), QRectF())] * len(self.rows)
        for g in groups:
            gh = max(usable * (self.heights[i] if i < len(self.heights) else 1.0) for i in g)
            self._labels.append((self.rows[g[0]].chromosome.name, QRectF(x, top, lab - 6, gh)))
            x += lab
            for i in g:
                h = usable * (self.heights[i] if i < len(self.heights) else 1.0)
                fish = QRectF(x, top, bar, h)
                asm = QRectF(x + bar + self.ASM_GAP * f, top, self.ASM_W * f, h) if has_asm else QRectF()
                self._rects[i] = (fish, asm)
                x += bar + asm_w + gap
            x += ggap - gap

    def _bars(self, i: int) -> tuple[QRectF, QRectF]:
        if self._size != (self.width(), self.height()) or len(self._rects) != len(self.rows):
            self._size = (self.width(), self.height())
            self._layout(self.width(), self.height())
        return self._rects[i] if i < len(self._rects) else (QRectF(), QRectF())

    def _band_rect(self, fish: QRectF, row: idg.ChromosomeRow, band: idg.FishBand, fraction: float | None = None) -> QRectF:
        f = band.fraction if fraction is None else fraction
        y = fish.top() + f * fish.height()
        return QRectF(fish.left() - 13, y - 7, 13, 14)

    def colour(self, probe: str) -> QColor:
        return QColor(self.colours.get(probe, "#8A8D93"))

    def _view(self, row: idg.ChromosomeRow) -> kt.KaryoChromosome:
        k = self.karyotype
        if not row.homologue or k is None:
            return row.chromosome
        suffixes = k.homologue_suffixes()
        i = suffixes.index(row.homologue) if row.homologue in suffixes else -1
        return kt.homologue_view(row.chromosome, i)

    def _band_brush(self, b: kt.Band) -> QBrush:
        probe = kt.manual_probe(b.kind)
        if b.color:
            return QBrush(QColor(b.color))
        if probe:
            return QBrush(self.colour(probe))
        if b.kind.startswith("Secondary"):
            return QBrush(QColor("#8A8D93"), Qt.BrushStyle.BDiagPattern)
        return QBrush(QColor(BAND_COLOURS.get(b.kind, "#8A8D93")))

    # ---------------------------------------------------------------- painting

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        self.paint_on(p, self.width(), self.height())

    def paint_on(self, p: QPainter, width: int, height: int) -> None:
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.fillRect(QRectF(0, 0, width, height), QColor("#FFFFFF"))
        if not self.rows:
            p.setPen(QColor("#5A5D63"))
            p.drawText(QRectF(0, 0, width, height), Qt.AlignmentFlag.AlignCenter,
                       "No karyotype yet: Signals > 1 Karyotype, or Edit karyotype on the right.")
            return
        self._size = (width, height)
        self._layout(width, height)
        small = QFont(self.font())
        small.setPixelSize(11)
        bold = QFont(self.font())
        bold.setPixelSize(15)
        bold.setBold(True)
        if self.title:
            p.setFont(bold)
            p.setPen(QColor("#1C1C1E"))
            p.drawText(QRectF(20, 6, width - 40, 22), Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, self.title)
        outline = QColor("#6C6F75")
        for i, row in enumerate(self.rows):
            fish, asm = self._rects[i]
            if self._drag_slot == i:
                dx = self._drag_x - fish.center().x()
                fish, asm = fish.translated(dx, 0), asm.translated(dx, 0)
            view = self._view(row)
            path = chromosome_path(fish, view)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor("#FFFFFF"))
            p.drawPath(path)
            if self.show_bands:
                p.save()
                p.setClipPath(path)
                for b in view.bands:
                    p.setBrush(self._band_brush(b))
                    p.drawRect(QRectF(fish.left(), fish.top() + b.start * fish.height(), fish.width(),
                                      max(1.5, (b.end - b.start) * fish.height())))
                p.restore()
            p.setPen(QPen(outline, 1.4))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawPath(path)
            cy = fish.top() + view.centromere_fraction * fish.height()
            p.setBrush(QColor("#FFFFFF"))
            p.drawEllipse(QPointF(fish.center().x(), cy), 5, 5)
            if self.show_measured:
                for j, band in enumerate(row.fish):
                    frac = self._drag_frac if self._drag_band == (i, j) else band.fraction
                    y = fish.top() + frac * fish.height()
                    tri = QPolygonF([QPointF(fish.left() - 11, y - 5.5), QPointF(fish.left() - 11, y + 5.5),
                                     QPointF(fish.left() - 2, y)])
                    p.setPen(QPen(QColor("#1D1D1F"), 0.8))
                    p.setBrush(self.colour(band.probe))
                    p.drawPolygon(tri)
                    if self._has_asm() and row.sequence_length and j not in row.matched_fish:
                        p.setPen(QPen(QColor("#A3213B"), 1.2, Qt.PenStyle.DashLine))
                        p.setBrush(Qt.BrushStyle.NoBrush)
                        p.drawRect(QRectF(fish.left() - 14, y - 8, 16, 16))
            if not asm.isNull() and row.sequence_length:
                p.setPen(QPen(QColor("#8A8D93"), 1, Qt.PenStyle.DashLine))
                p.setBrush(QColor("#EEF0F4"))
                p.drawRoundedRect(asm, 3, 3)
                p.setPen(Qt.PenStyle.NoPen)
                for array in row.silico:
                    top = asm.top() + array.start * asm.height()
                    p.setBrush(self.colour(array.probe))
                    p.drawRect(QRectF(asm.left() + 1, top, asm.width() - 2, max(3.0, (array.end - array.start) * asm.height())))
            p.setFont(small)
            p.setPen(QColor("#5A5D63"))
            if row.homologue:
                p.drawText(QRectF(fish.left() - 12, fish.bottom() + 8, fish.width() + 24, 14),
                           Qt.AlignmentFlag.AlignHCenter, row.homologue.lstrip("."))
            if self.assembly_active and self.show_assembly:
                text = row.sequence or "not linked"
                y0 = fish.bottom() + (24 if row.homologue else 8)
                p.drawText(QRectF(fish.left() - 20, y0, fish.width() + 40, 14), Qt.AlignmentFlag.AlignHCenter,
                           p.fontMetrics().elidedText(text, Qt.TextElideMode.ElideMiddle, int(fish.width() + 38)))
        p.setFont(bold)
        p.setPen(QColor("#1C1C1E"))
        for name, rect in self._labels:
            p.drawText(rect, Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, name)

    # ---------------------------------------------------------------- mouse

    def _hit(self, pos: QPointF):
        for i, row in enumerate(self.rows):
            fish, asm = self._bars(i)
            if self.show_measured:
                for j, band in enumerate(row.fish):
                    if self._band_rect(fish, row, band).adjusted(-3, -2, 3, 2).contains(pos):
                        return ("band", i, j)
            if not asm.isNull() and row.sequence_length:
                for j, array in enumerate(row.silico):
                    top = asm.top() + array.start * asm.height()
                    h = max(5.0, (array.end - array.start) * asm.height())
                    if QRectF(asm.left() - 3, top - 2, asm.width() + 6, h + 4).contains(pos):
                        return ("array", i, j)
            if fish.adjusted(-6, -6, 6, 40).contains(pos) or (not asm.isNull() and asm.contains(pos)):
                return ("slot", i, -1)
        return None

    def _row_at_x(self, x: float) -> int:
        best, dist = 0, None
        for i in range(len(self.rows)):
            d = abs(self._bars(i)[0].center().x() - x)
            if dist is None or d < dist:
                best, dist = i, d
        return best

    def mousePressEvent(self, event) -> None:  # noqa: N802
        hit = self._hit(event.position())
        if hit is None or event.button() != Qt.MouseButton.LeftButton:
            return
        kind, i, j = hit
        self._press = (event.position(), hit)
        if kind == "band" and self.drag_enabled:
            self._drag_band = (i, j)
            self._drag_frac = self.rows[i].fish[j].fraction
        elif self.reorder_enabled and kind in ("slot", "band", "array"):
            self._drag_slot = i
            self._drag_x = event.position().x()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        pos = event.position()
        if self._drag_band is not None:
            i, _j = self._drag_band
            fish, _asm = self._bars(i)
            self._drag_frac = max(0.0, min(1.0, (pos.y() - fish.top()) / max(fish.height(), 1.0)))
            self.update()
            return
        if self._drag_slot is not None:
            self._drag_x = pos.x()
            self.update()
            return
        hit = self._hit(pos)
        if hit is None:
            QToolTip.hideText()
            self.unsetCursor()
            return
        kind, i, j = hit
        row = self.rows[i]
        if kind == "band":
            b = row.fish[j]
            self.setCursor(Qt.CursorShape.SizeVerCursor if self.drag_enabled else Qt.CursorShape.ArrowCursor)
            where = f"{b.arm} {b.position}" if b.arm in ("p", "q") else b.position
            text = (f"{b.probe}, {where}\nseen on {b.seen} of {b.total} marked homologues in {b.cells} "
                    f"cell{'s' if b.cells != 1 else ''}")
            if self._has_asm() and j not in row.matched_fish and row.sequence_length:
                text += "\nno array of this probe nearby in the assembly"
            QToolTip.showText(event.globalPosition().toPoint(), text, self)
        elif kind == "array":
            a = row.silico[j]
            mb = row.sequence_length / 1e6
            text = (f"{a.probe} array on {row.sequence}\n{a.start * mb:.2f} to {a.end * mb:.2f} Mb, "
                    f"{a.bp / 1000:.1f} kb covered by {a.hits} hit{'s' if a.hits != 1 else ''}")
            QToolTip.showText(event.globalPosition().toPoint(), text, self)
            self.unsetCursor()
        else:
            self.setCursor(Qt.CursorShape.OpenHandCursor if self.reorder_enabled else Qt.CursorShape.ArrowCursor)
            QToolTip.hideText()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        press, self._press = self._press, None
        if press is not None and (event.position() - press[0]).manhattanLength() < 4:
            kind, i, j = press[1]
            self._drag_band = None
            self._drag_slot = None
            self.update()
            if kind == "band":
                self.bandClicked.emit(i, j, event.globalPosition().toPoint())
            elif kind == "array":
                self.arrayClicked.emit(i, j, event.globalPosition().toPoint())
            return
        if self._drag_band is not None:
            i, j = self._drag_band
            self._drag_band = None
            band = self.rows[i].fish[j]
            if abs(self._drag_frac - band.fraction) > 0.005:
                self.bandMoved.emit(band, self._drag_frac)
            self.update()
        elif self._drag_slot is not None:
            src = self._drag_slot
            self._drag_slot = None
            dst = self._row_at_x(event.position().x())
            names: list[str] = []
            for r in self.rows:
                if r.chromosome.name not in names:
                    names.append(r.chromosome.name)
            moving = self.rows[src].chromosome.name
            target = self.rows[dst].chromosome.name
            if moving != target:
                names.remove(moving)
                names.insert(names.index(target) + (1 if dst > src else 0), moving)
                self.orderChanged.emit(names)
            self.update()

    # ---------------------------------------------------------------- export

    def to_image(self, scale: float = 2.0) -> QImage:
        size = QSize(max(self.width(), 600), max(self.height(), 420))
        image = QImage(int(size.width() * scale), int(size.height() * scale), QImage.Format.Format_ARGB32)
        image.fill(QColor("#FFFFFF"))
        p = QPainter(image)
        p.scale(scale, scale)
        self.paint_on(p, size.width(), size.height())
        p.end()
        self._size = (0, 0)
        return image

    def to_svg(self, path: Path) -> None:
        from PySide6.QtSvg import QSvgGenerator

        size = QSize(max(self.width(), 600), max(self.height(), 420))
        gen = QSvgGenerator()
        gen.setFileName(str(path))
        gen.setSize(size)
        gen.setViewBox(QRectF(0, 0, size.width(), size.height()))
        gen.setTitle("Idiogram")
        p = QPainter(gen)
        self.paint_on(p, size.width(), size.height())
        p.end()
        self._size = (0, 0)


def _font(p: QPainter, family: str, size: float, *, bold: bool = False, italic: bool = False) -> QFont:
    """A font whose size is ``size`` painter units (1 unit = 1 pt in the paper drawing) on any device."""
    f = QFont(family)
    dpi = p.device().logicalDpiY() if p.device() is not None else 72
    f.setPointSizeF(max(1.0, size * 72.0 / max(dpi, 1)))
    f.setBold(bold)
    f.setItalic(italic)
    return f


def _bar_path(r: QRectF) -> QPainterPath:
    path = QPainterPath()
    path.addRoundedRect(r, r.width() / 2, r.width() / 2)
    return path


def _draw_chromosome(p: QPainter, r: QRectF, chrom: kt.KaryoChromosome, bands: list[tuple[float, str, bool]],
                     colours: dict[str, str], *, outline: QColor, line: float) -> None:
    """One chromosome: white bar, karyotype bands, FISH signals (fraction, probe, centromeric) and the centromere."""
    path = chromosome_path(r, chrom, 0.0) if any(b.satellite for b in chrom.bands) else _bar_path(r)
    p.setPen(Qt.PenStyle.NoPen)
    p.setBrush(QColor("#FFFFFF"))
    p.drawPath(path)
    p.save()
    p.setClipPath(path)
    for b in chrom.bands:
        rect = QRectF(r.left(), r.top() + b.start * r.height(), r.width(), max(1.5, (b.end - b.start) * r.height()))
        probe = kt.manual_probe(b.kind)
        if probe:
            p.setBrush(QColor(colours.get(probe, "#8A8D93")))
        elif b.kind.startswith("Secondary"):
            p.setBrush(QBrush(outline, Qt.BrushStyle.BDiagPattern))
        else:
            p.setBrush(QColor(b.color or BAND_COLOURS.get(b.kind, "#8A8D93")))
        p.drawRect(rect)
    for f, probe, cen in bands:
        h = r.width() * (0.6 if cen else 0.45)
        y = r.top() + f * r.height()
        p.setBrush(QColor(colours.get(probe, "#8A8D93")))
        p.drawRect(QRectF(r.left(), y - h / 2, r.width(), h))
    p.restore()
    p.setPen(QPen(outline, line))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawPath(path)
    cy = r.top() + chrom.centromere_fraction * r.height()
    ch = max(3.0, r.width() * 0.28)
    p.setBrush(QColor("#FFFFFF"))
    p.drawRoundedRect(QRectF(r.left() - r.width() * 0.2, cy - ch / 2, r.width() * 1.4, ch), ch / 2, ch / 2)


def _fish_by_chrom(bands: list[idg.FishBand]) -> dict[str, list[tuple[float, str, bool]]]:
    out: dict[str, list] = {}
    for b in bands:
        out.setdefault(b.chromosome, []).append((b.fraction, b.probe, b.arm == "cen"))
    return out


def paper_height(k: kt.Karyotype, table_rows: int, opts: dict) -> float:
    lh = float(opts.get("size", 8.0)) * 1.6
    h = 20 + (lh * 1.6 if opts.get("formula") else 0) + 200 + lh * 2.6 + lh * 1.8
    if opts.get("table") and table_rows:
        h += (table_rows + 1) * lh + 10
    return h


def paint_paper(p: QPainter, w: float, h: float, k: kt.Karyotype | None, bands: list[idg.FishBand],
                colours: dict[str, str], opts: dict) -> None:
    """The idiogram as a journal figure, in points: black and white chromosomes, coloured signals, legend, table."""
    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    p.fillRect(QRectF(0, 0, w, h), QColor("#FFFFFF"))
    if k is None or not k.chromosomes:
        return
    family, size = str(opts.get("font") or "Arial"), float(opts.get("size") or 8.0)
    lh = size * 1.6
    black = QColor("#000000") if opts.get("outlines", True) else QColor("#6C6F75")
    line = max(0.6, size * 0.12)
    x0, y = 16.0, 16.0
    if opts.get("formula", True):
        parts = [("Idiogram of ", True, False)]
        if k.species:
            parts.append((k.species, True, True))
        else:
            parts.append((k.name, True, False))
        base = f"   {k.formula}   (x = {k.basic_number})"
        parts.append((base, False, False))
        x = x0
        for text, bold, italic in parts:
            p.setFont(_font(p, family, size * 1.15, bold=bold, italic=italic))
            p.setPen(QColor("#000000"))
            p.drawText(QPointF(x, y + size * 1.15), text)
            x += p.fontMetrics().horizontalAdvance(text)
        y += lh * 1.6
    area = 200.0
    n = len(k.chromosomes)
    right = 56.0 if opts.get("scale", True) else 16.0
    slot = (w - x0 - right) / n
    bar_w = max(6.0, min(22.0, slot * 0.32))
    longest = max((c.length for c in k.chromosomes), default=0.0)
    fish = _fish_by_chrom(bands)
    bottom = y + area
    for i, c in enumerate(k.chromosomes):
        hh = area * (c.length / longest) if longest > 0 and c.length > 0 else area * 0.6
        cx = x0 + slot * (i + 0.5)
        r = QRectF(cx - bar_w / 2, bottom - hh, bar_w, hh)
        _draw_chromosome(p, r, c, fish.get(c.name, []), colours, outline=black, line=line)
        p.setPen(QColor("#000000"))
        p.setFont(_font(p, family, size * 1.15, bold=True))
        p.drawText(QRectF(cx - slot / 2, bottom + 4, slot, lh), Qt.AlignmentFlag.AlignCenter, c.name)
        if opts.get("types", True) and c.type:
            p.setFont(_font(p, family, size * 0.9))
            p.setPen(QColor("#45484D"))
            p.drawText(QRectF(cx - slot / 2, bottom + 4 + lh, slot, lh), Qt.AlignmentFlag.AlignCenter, c.type)
    if opts.get("scale", True) and longest > 0:
        value = idg.nice_scale(longest)
        bar_h = area * value / longest
        sx = w - right + 14
        p.setPen(QPen(QColor("#000000"), line * 1.5))
        p.drawLine(QPointF(sx, bottom), QPointF(sx, bottom - bar_h))
        p.setFont(_font(p, family, size))
        p.drawText(QPointF(sx + 5, bottom - bar_h / 2 + size / 3), f"{value:g} {k.unit}")
    y = bottom + lh * 2.6
    # legend
    probes = list(dict.fromkeys([b.probe for b in bands] + [kt.manual_probe(b.kind) for c in k.chromosomes
                                                             for b in c.bands if kt.manual_probe(b.kind)]))
    kinds = list(dict.fromkeys(b.kind for c in k.chromosomes for b in c.bands if not kt.manual_probe(b.kind)))
    x = x0
    p.setFont(_font(p, family, size))
    for label, brush in [(pr, QBrush(QColor(colours.get(pr, "#8A8D93")))) for pr in probes] + [
            ({"Secondary constriction (NOR)": "NOR", "DAPI+ heterochromatin": "DAPI+", "CMA+ band": "CMA+"}.get(kd, kd),
             QBrush(black, Qt.BrushStyle.BDiagPattern) if kd.startswith("Secondary")
             else QBrush(QColor(BAND_COLOURS.get(kd, "#8A8D93")))) for kd in kinds]:
        sw = QRectF(x, y, size * 1.4, size)
        p.setPen(QPen(black, line * 0.6) if brush.style() != Qt.BrushStyle.SolidPattern else Qt.PenStyle.NoPen)
        p.setBrush(brush)
        p.drawRect(sw)
        p.setPen(QColor("#000000"))
        p.drawText(QPointF(x + size * 1.9, y + size * 0.85), label)
        x += size * 1.9 + p.fontMetrics().horizontalAdvance(label) + size * 2
    y += lh * 1.8
    # table of positions
    rows = idg.paper_table(k, bands) if opts.get("table", True) else []
    if rows:
        p.setPen(QPen(QColor("#000000"), line))
        p.drawLine(QPointF(x0, y - lh * 0.3), QPointF(w - 16, y - lh * 0.3))
        widths = [size * 5, size * 7.5, size * 4.5, size * 5]
        heads = ["Chr", "Length", "CI", "Type", "Signals (seen)"]
        for n_row, cells in enumerate([heads] + rows):
            p.setFont(_font(p, family, size * 0.95, bold=n_row == 0))
            x = x0
            for col, text in enumerate(cells):
                width = widths[col] if col < len(widths) else w - 16 - x
                elided = p.fontMetrics().elidedText(text, Qt.TextElideMode.ElideRight, int(width))
                p.drawText(QRectF(x, y, width, lh), Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, elided)
                x += width
            y += lh


def paint_pair(p: QPainter, w: float, h: float, top: tuple[kt.Karyotype, list[idg.FishBand]] | None,
               bottom: tuple[kt.Karyotype, list[idg.FishBand]] | None, colours: dict[str, str], link: bool,
               font: QFont) -> None:
    """Two karyotypes, one above the other; dashed lines join the same probe on the same chromosome."""
    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    p.fillRect(QRectF(0, 0, w, h), QColor("#FFFFFF"))
    if top is None or bottom is None:
        p.setPen(QColor("#5A5D63"))
        p.setFont(font)
        p.drawText(QRectF(0, 0, w, h), Qt.AlignmentFlag.AlignCenter,
                   "Choose two karyotypes (Top and Bottom) to compare. Make a second one under Layout > Edit karyotype.")
        return
    (tk, tb), (bk, bb) = top, bottom
    n = max(len(tk.chromosomes), len(bk.chromosomes), 1)
    left = 150.0
    slot = (w - left - 20) / n
    bar_w = max(8.0, min(24.0, slot * 0.3))
    half = (h - 40) / 2
    area = half - 46
    longest = max([c.length for c in tk.chromosomes + bk.chromosomes] or [0.0])
    grey = QColor("#6C6F75")
    rects: dict[tuple[int, int], QRectF] = {}
    bold = QFont(font)
    bold.setBold(True)
    small = QFont(font)
    small.setPointSizeF(max(7.0, font.pointSizeF() - 1))
    for row, (k, bands) in enumerate(((tk, tb), (bk, bb))):
        base = 20 + row * (half + 20) + area
        label = QFont(bold)
        label.setItalic(bool(k.species))
        p.setFont(label)
        p.setPen(QColor("#1C1C1E"))
        p.drawText(QRectF(8, base - area, left - 16, 40), Qt.AlignmentFlag.AlignLeft | Qt.TextFlag.TextWordWrap,
                   k.species or k.name)
        p.setFont(small)
        p.setPen(QColor("#6E6E73"))
        p.drawText(QRectF(8, base - area + 40, left - 16, 18), Qt.AlignmentFlag.AlignLeft, k.formula)
        fish = _fish_by_chrom(bands)
        for i, c in enumerate(k.chromosomes):
            hh = area * (c.length / longest) if longest > 0 and c.length > 0 else area * 0.6
            cx = left + slot * (i + 0.5)
            r = QRectF(cx - bar_w / 2, base - hh, bar_w, hh)
            rects[(row, i)] = r
            _draw_chromosome(p, r, c, fish.get(c.name, []), colours, outline=grey, line=1.4)
            p.setFont(bold)
            p.setPen(QColor("#1C1C1E"))
            p.drawText(QRectF(cx - slot / 2, base + 4, slot, 18), Qt.AlignmentFlag.AlignCenter, c.name)
    links, only_top, only_bottom = idg.pair_compare(tk, tb, bk, bb)
    if link:
        for ln in links:
            a, b = rects.get((0, ln.index)), rects.get((1, ln.bottom_index))
            if a is None or b is None:
                continue
            y1 = a.top() + ln.top * a.height()
            y2 = b.top() + ln.bottom * b.height()
            path = QPainterPath(QPointF(a.right() + 4, y1))
            mid = (a.bottom() + b.top()) / 2
            path.cubicTo(QPointF(a.right() + 14, mid), QPointF(b.right() + 14, mid), QPointF(b.right() + 4, y2))
            p.setPen(QPen(QColor(colours.get(ln.probe, "#8A8D93")), 1.5, Qt.PenStyle.DashLine))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawPath(path)
    p.setFont(small)
    for row, items, k in ((0, only_top, tk), (1, only_bottom, bk)):
        for i, band in items:
            r = rects.get((row, i))
            if r is None:
                continue
            y = r.top() + band.fraction * r.height()
            colour = QColor(colours.get(band.probe, "#8A8D93")).darker(150)
            p.setPen(colour)
            p.drawText(QPointF(r.right() + 8, y + 4), f"{band.probe} only in {k.species or k.name}")


class PaperCanvas(QWidget):
    """Screen preview of the paper drawing, scaled to the width."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.karyotype: kt.Karyotype | None = None
        self.bands: list[idg.FishBand] = []
        self.colours: dict[str, str] = {}
        self.opts: dict = dict(idg.PAPER_DEFAULTS)
        self.setMinimumSize(300, 300)

    def logical_size(self) -> tuple[float, float]:
        rows = len(idg.paper_table(self.karyotype, self.bands)) if self.karyotype else 0
        return idg.PAPER_WIDTH_PT, paper_height(self.karyotype, rows, self.opts) if self.karyotype else 300.0

    def paint_into(self, p: QPainter) -> None:
        w, h = self.logical_size()
        paint_paper(p, w, h, self.karyotype, self.bands, self.colours, self.opts)

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.fillRect(self.rect(), self.palette().color(self.palette().ColorRole.Window))
        w, h = self.logical_size()
        scale = min((self.width() - 32) / w, (self.height() - 32) / h)
        ox, oy = (self.width() - w * scale) / 2, 16
        p.fillRect(QRectF(ox + 2, oy + 3, w * scale, h * scale), QColor(0, 0, 0, 30))
        p.translate(ox, oy)
        p.scale(scale, scale)
        self.paint_into(p)


class PairCanvas(QWidget):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.top = None
        self.bottom = None
        self.colours: dict[str, str] = {}
        self.link = True
        self.setMinimumSize(300, 420)

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        paint_pair(p, self.width(), self.height(), self.top, self.bottom, self.colours, self.link, self.font())


class IdiogramEditCanvas(KaryoPreview):
    """Edit view: the karyotype drawn large; drag length, centromere and band edges, paint bands, place signals."""

    changed = Signal()  # an edit is finished: save the karyotype

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumSize(300, 420)
        self.tool = "select"
        self.kind = kt.BAND_KINDS[0]
        self.kind_colour = ""
        self.probe_colours: dict[str, str] = {}
        self.fish: dict[str, list[tuple[float, str, bool]]] = {}  # chromosome -> (fraction, probe, centromeric)
        self._paint: tuple[int, float, float] | None = None
        self._len_drag: tuple[int, float, float] | None = None

    @property
    def MAX_H(self) -> float:  # noqa: N802 - overrides the class constant of KaryoPreview
        return max(120.0, self.height() - 110.0)

    def band_colour(self, band: kt.Band) -> QColor:
        probe = kt.manual_probe(band.kind)
        if probe:
            return QColor(self.probe_colours.get(probe, "#8A8D93"))
        return super().band_colour(band)

    def _fraction(self, i: int, y: float) -> float:
        g = self._group_rect(i, self.karyotype.chromosomes[i])
        return max(0.0, min(1.0, (y - g.top()) / max(g.height(), 1.0)))

    def _on_top(self, pos: QPointF) -> bool:
        k = self.karyotype
        if k is None or not 0 <= self.index < len(k.chromosomes):
            return False
        g = self._group_rect(self.index, k.chromosomes[self.index])
        return g.left() - 6 <= pos.x() <= g.right() + 6 and abs(pos.y() - (g.top() - 4)) <= 6

    def mousePressEvent(self, event) -> None:  # noqa: N802
        k = self.karyotype
        if k is None or event.button() != Qt.MouseButton.LeftButton:
            return
        pos = event.position()
        if self.tool == "select":
            if self._on_top(pos):
                c = k.chromosomes[self.index]
                g = self._group_rect(self.index, c)
                longest = max((x.length for x in k.chromosomes), default=0.0) or 1.0
                self._len_drag = (self.index, longest, g.bottom())
                return
            super().mousePressEvent(event)
            return
        hit = self._hit(pos)
        if hit is None:
            return
        i = hit[1]
        if i != self.index:
            self.index = i
            self.band_index = -1
            self.selected.emit(i)
        f = self._fraction(i, pos.y())
        if self.tool == "paint":
            self._paint = (i, f, f)
        else:
            c = k.chromosomes[i]
            c.bands.append(kt.Band(kt.SIGNAL_KIND_PREFIX + self.kind, max(0.0, f - 0.012), min(1.0, f + 0.012)))
            self.band_index = len(c.bands) - 1
            self.selected.emit(i)
            self.changed.emit()
        self.update()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        pos = event.position()
        k = self.karyotype
        if self._len_drag is not None and k is not None:
            i, longest, bottom = self._len_drag
            h = max(20.0, bottom - pos.y())
            k.chromosomes[i].length = round(max(0.02 * longest, h / self.MAX_H * longest), 3)
            self.edited.emit()
            self.update()
            return
        if self._paint is not None:
            i, f0, _f = self._paint
            self._paint = (i, f0, self._fraction(i, pos.y()))
            self.update()
            return
        if self.tool == "select":
            super().mouseMoveEvent(event)
            if self._drag is None and self._on_top(pos):
                self.setCursor(Qt.CursorShape.SizeVerCursor)
        else:
            self.setCursor(Qt.CursorShape.CrossCursor if self._hit(pos) else Qt.CursorShape.ArrowCursor)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        k = self.karyotype
        if self._len_drag is not None:
            self._len_drag = None
            self.changed.emit()
        elif self._paint is not None and k is not None:
            i, a, b = self._paint
            self._paint = None
            if abs(b - a) > 0.01:
                c = k.chromosomes[i]
                c.bands.append(kt.Band(self.kind, min(a, b), max(a, b), self.kind_colour))
                self.band_index = len(c.bands) - 1
                self.selected.emit(i)
                self.changed.emit()
            self.update()
        else:
            dragged = self._drag is not None
            super().mouseReleaseEvent(event)
            if dragged:
                self.changed.emit()

    def paintEvent(self, event) -> None:  # noqa: N802
        super().paintEvent(event)
        k = self.karyotype
        if k is None or not k.chromosomes:
            return
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        accent = self.palette().color(self.palette().ColorRole.Highlight)
        for i, c in enumerate(k.chromosomes):  # the FISH signals marked in Signals, pooled
            for f, probe, cen in self.fish.get(c.name, []):
                for r in self._bars(i, c):
                    h = 12.0 if cen else 8.0
                    y = r.top() + f * r.height()
                    clip = QPainterPath()
                    clip.addRoundedRect(r, r.width() / 2, r.width() / 2)
                    p.save()
                    p.setClipPath(clip)
                    p.setPen(Qt.PenStyle.NoPen)
                    p.setBrush(QColor(self.probe_colours.get(probe, "#8A8D93")))
                    p.drawRect(QRectF(r.left() + 1, y - h / 2, r.width() - 2, h))
                    p.restore()
        if self.tool == "select" and 0 <= self.index < len(k.chromosomes):
            g = self._group_rect(self.index, k.chromosomes[self.index])
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(accent)
            p.drawRoundedRect(QRectF(g.center().x() - 16, g.top() - 9, 32, 7), 3, 3)
        if self._paint is not None:
            i, a, b = self._paint
            g = self._group_rect(i, k.chromosomes[i])
            top = g.top() + min(a, b) * g.height()
            colour = QColor(self.kind_colour or BAND_COLOURS.get(self.kind, "#8A8D93"))
            colour.setAlpha(140)
            p.setBrush(colour)
            p.drawRect(QRectF(g.left(), top, g.width(), abs(b - a) * g.height()))


class _Popover(QFrame):
    """Small floating panel shown where the user clicked (closes on a click outside)."""

    def __init__(self, parent: QWidget) -> None:
        super().__init__(parent, Qt.WindowType.Popup)
        self.setFrameShape(QFrame.Shape.StyledPanel)
        self.setAutoFillBackground(True)
        self.lay = QVBoxLayout(self)
        self.lay.setContentsMargins(14, 12, 14, 12)
        self.lay.setSpacing(8)

    def title(self, text: str) -> None:
        label = QLabel(f"<b>{text}</b>", self)
        label.setWordWrap(True)
        self.lay.addWidget(label)

    def note(self, text: str) -> QLabel:
        label = QLabel(text, self)
        label.setWordWrap(True)
        label.setStyleSheet("color: gray;")
        self.lay.addWidget(label)
        return label

    def popup_at(self, pos: QPoint) -> None:
        self.setFixedWidth(520)
        self.adjustSize()
        screen = QApplication.screenAt(pos) or QApplication.primaryScreen()
        area = screen.availableGeometry()
        x = min(pos.x() + 12, area.right() - self.width() - 8)
        y = min(pos.y() + 12, area.bottom() - self.height() - 8)
        self.move(max(area.left(), x), max(area.top(), y))
        self.show()


def _cell_name(image: ci.CytoImage) -> str:
    stem = Path(image.original_name or image.file).stem
    m = re.search(r"(\d+)\D*$", stem)
    return f"cell {m.group(1).zfill(2)}" if m else stem[:14]


class IdiogramTab(QWidget):
    name = "Idiogram"
    karyotypesChanged = Signal()
    figureRequested = Signal()
    openSignalRequested = Signal(str, str)  # image id, signal id
    openRequested = Signal(object)  # Path of a file to open in the viewer

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._folder: Path | None = None
        self._settings = idg.IdiogramSettings()
        self._karyotypes: list[kt.Karyotype] = []
        self._karyotype: kt.Karyotype | None = None
        self._cache: dict = {}
        self._lengths: dict[str, int] = {}
        self._fish: idg.FishData | None = None
        self._rows: list[idg.ChromosomeRow] = []
        self._colours: dict[str, str] = {}
        self._dirty = True
        self._busy = False

        # ---- main column
        self._title = QLabel(self)
        self._title.setTextFormat(Qt.TextFormat.RichText)
        f = self._title.font()
        f.setPointSizeF(f.pointSizeF() + 2)
        f.setBold(True)
        self._title.setFont(f)
        self._title.setWordWrap(True)
        self._title.setMinimumWidth(120)
        self._title.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self._legend = QHBoxLayout()
        self._legend.setSpacing(14)
        self._title.hide()
        head = QHBoxLayout()
        head.setSpacing(8)

        self._canvas = IdiogramCanvas(self)
        self._canvas.bandMoved.connect(self._on_band_moved)
        self._canvas.orderChanged.connect(self._on_order_dragged)
        self._canvas.bandClicked.connect(self._on_band_clicked)
        self._canvas.arrayClicked.connect(self._on_array_clicked)
        self._mode = "compare"
        self._edit = IdiogramEditCanvas(self)
        self._edit.selected.connect(self._on_edit_selected)
        self._edit.edited.connect(self._sync_edit_side)
        self._edit.changed.connect(self._on_edit_changed)
        frame = QFrame(self)
        frame.setStyleSheet("QFrame#idioFrame { border:1px solid #E6E6EA; border-radius:8px; background:#FFFFFF; }")
        frame.setObjectName("idioFrame")
        fl = QVBoxLayout(frame)
        fl.setContentsMargins(8, 12, 8, 8)
        scroll = QScrollArea(frame)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._canvas2 = IdiogramCanvas(self)
        self._canvas2.hide()
        self._canvas2.bandClicked.connect(lambda i, j, pos: None)
        canvases = QWidget(self)
        canvases.setObjectName("idioCanvases")
        canvases.setStyleSheet("QWidget#idioCanvases { background:#FFFFFF; }")
        cl = QVBoxLayout(canvases)
        cl.setContentsMargins(0, 0, 0, 0)
        cl.setSpacing(10)
        cl.addWidget(self._canvas, 1)
        cl.addWidget(self._canvas2, 1)
        scroll.setWidget(canvases)
        scroll.setStyleSheet("QScrollArea { background:#FFFFFF; }")
        self._stack = QStackedWidget(frame)
        self._stack.addWidget(scroll)
        self._stack.addWidget(self._edit)
        fl.addWidget(self._stack)

        # view switch and the Edit tools
        self._mode_btns: dict[str, QPushButton] = {}
        modes = QHBoxLayout()
        modes.setSpacing(4)
        group = QButtonGroup(self)
        for key, text in (("compare", "Compare"), ("edit", "Edit"), ("pair", "Two karyotypes"), ("paper", "Paper")):
            b = QPushButton(text, self)
            b.setCheckable(True)
            b.setChecked(key == "compare")
            b.clicked.connect(lambda _c=False, k=key: self._set_mode(k))
            group.addButton(b)
            self._mode_btns[key] = b
            modes.addWidget(b)
        modes.addStretch(1)
        self._tools_host = QWidget(self)
        tools = QHBoxLayout(self._tools_host)
        tools.setContentsMargins(0, 0, 0, 0)
        tools.setSpacing(4)
        tools.addWidget(QLabel("Tool", self._tools_host))
        tool_group = QButtonGroup(self)
        self._tool_btns: dict[str, QPushButton] = {}
        for key, text in TOOLS:
            b = QPushButton(text, self._tools_host)
            b.setCheckable(True)
            b.setChecked(key == "select")
            b.clicked.connect(lambda _c=False, k=key: self._set_tool(k))
            tool_group.addButton(b)
            self._tool_btns[key] = b
            tools.addWidget(b)
        tools.addSpacing(8)
        self._kind = NoWheelComboBox(self._tools_host)
        self._kind.setMinimumWidth(170)
        self._kind.activated.connect(self._on_kind)
        tools.addWidget(self._kind)
        self._tools_host.hide()
        modes.addWidget(self._tools_host)
        self._pair_host = QWidget(self)
        pair = QHBoxLayout(self._pair_host)
        pair.setContentsMargins(0, 0, 0, 0)
        pair.setSpacing(6)
        self._top_box = NoWheelComboBox(self._pair_host)
        self._bottom_box = NoWheelComboBox(self._pair_host)
        for text, box in (("Top", self._top_box), ("Bottom", self._bottom_box)):
            box.setMinimumWidth(190)
            box.activated.connect(self._on_pair_choice)
            pair.addWidget(QLabel(text, self._pair_host))
            pair.addWidget(box)
            pair.addSpacing(6)
        self._link = QCheckBox("Link the same probe", self._pair_host)
        self._link.toggled.connect(self._on_pair_choice)
        pair.addWidget(self._link)
        self._pair_host.hide()
        modes.addWidget(self._pair_host)
        self._pair = PairCanvas(frame)
        self._paper = PaperCanvas(frame)
        self._stack.addWidget(self._pair)
        self._stack.addWidget(self._paper)

        self._note_fish = _section_label("")
        self._note_asm = _section_label("")
        self._note_warn = QLabel(self)
        self._note_warn.setWordWrap(True)
        self._note_warn.setStyleSheet("color:#A3213B; font-size:12px;")
        notes = QHBoxLayout()
        notes.setSpacing(20)
        notes.addWidget(self._note_fish, 1)
        notes.addWidget(self._note_asm, 1)

        main = QVBoxLayout()
        main.setContentsMargins(16, 12, 16, 12)
        main.setSpacing(14)
        main.addLayout(head)
        main.addWidget(frame, 1)
        main.addWidget(self._note_warn)
        for b in self._mode_btns.values():
            b.hide()
        self._note_fish.hide()
        self._note_asm.hide()
        main_host = QWidget(self)
        main_host.setObjectName("cytoMain")
        main_host.setLayout(main)

        # ---- side panel
        body = QWidget(self)
        body.setObjectName("cytoAside")
        side = QVBoxLayout(body)
        side.setContentsMargins(14, 12, 14, 12)
        side.setSpacing(14)

        self._layers_section = _CollapsibleSection("Layers", body, expanded=True)
        self._layer_boxes: dict[str, QCheckBox] = {}
        for key, text, mark in (("show_bands", "Karyotype bands", "■"), ("show_measured", "Measured signals", "▶"),
                                ("show_assembly", "Assembly", "┆")):
            row = QHBoxLayout()
            box = QCheckBox(text, self)
            box.toggled.connect(lambda on, k=key: self._on_layer(k, on))
            self._layer_boxes[key] = box
            sym = QLabel(mark, self)
            sym.setStyleSheet("color:#8A8D93;")
            row.addWidget(box, 1)
            row.addWidget(sym)
            self._layers_section.body_layout().addLayout(row)
        side.addWidget(self._layers_section)

        self._probe_section = _CollapsibleSection("Probes", body, expanded=True)
        self._probe_host = QWidget(self._probe_section)
        self._probe_lay = QVBoxLayout(self._probe_host)
        self._probe_lay.setContentsMargins(0, 0, 0, 0)
        self._probe_section.add_widget(self._probe_host)
        self._share = NoWheelComboBox(self)
        for value, text in idg.THRESHOLDS:
            self._share.addItem(text, value)
        self._share.activated.connect(self._on_share)
        r = QHBoxLayout()
        lab = QLabel("Seen in at least", self)
        lab.setFixedWidth(110)
        r.addWidget(lab)
        r.addWidget(self._share, 1)
        self._probe_section.body_layout().addLayout(r)
        side.addWidget(self._probe_section)

        self._measured_section = _CollapsibleSection("Measured, not in the karyotype", body, expanded=True)
        self._measured_grid = QGridLayout()
        self._measured_grid.setHorizontalSpacing(8)
        self._measured_grid.setVerticalSpacing(4)
        self._measured_section.body_layout().addLayout(self._measured_grid)
        self._add_measured = QPushButton("Add to karyotype", self)
        self._add_measured.setToolTip("Make these measured signals bands of the karyotype")
        self._add_measured.clicked.connect(self._add_measured_to_karyotype)
        self._measured_section.body_layout().addWidget(self._add_measured, 0, Qt.AlignmentFlag.AlignLeft)
        side.addWidget(self._measured_section)

        self._silico_section = _CollapsibleSection("Assembly", body, expanded=False)
        self._genome = NoWheelComboBox(self)
        self._genome.activated.connect(self._on_genome)
        g = QHBoxLayout()
        lab = QLabel("Genome", self)
        lab.setFixedWidth(90)
        g.addWidget(lab)
        g.addWidget(self._genome, 1)
        self._silico_section.body_layout().addLayout(g)
        self._source_host = QWidget(self._silico_section)
        self._source_lay = QVBoxLayout(self._source_host)
        self._source_lay.setContentsMargins(0, 0, 0, 0)
        self._source_lay.setSpacing(6)
        self._silico_section.add_widget(self._source_host)
        self._min_kb = QDoubleSpinBox(self)
        self._min_kb.setRange(0, 100000)
        self._min_kb.setDecimals(1)
        self._min_kb.setSuffix(" kb")
        self._min_kb.editingFinished.connect(self._on_numbers)
        self._tol = QDoubleSpinBox(self)
        self._tol.setRange(1, 50)
        self._tol.setDecimals(0)
        self._tol.setSuffix(" % of length")
        self._tol.editingFinished.connect(self._on_numbers)
        for text, w in (("Hide arrays below", self._min_kb), ("Agree within", self._tol)):
            r = QHBoxLayout()
            lab = QLabel(text, self)
            lab.setFixedWidth(110)
            r.addWidget(lab)
            r.addWidget(w, 1)
            self._silico_section.body_layout().addLayout(r)
        side.addWidget(self._silico_section)

        self._layout_section = _CollapsibleSection("Layout", body, expanded=False)
        self._karyo = NoWheelComboBox(self)
        self._karyo.activated.connect(self._on_karyotype)
        self._order = NoWheelComboBox(self)
        for key in idg.ORDERS:
            self._order.addItem(idg.ORDER_LABELS[key], key)
        self._order.activated.connect(self._on_layout)
        self._lengths_box = NoWheelComboBox(self)
        for key in idg.LENGTHS:
            self._lengths_box.addItem(idg.LENGTH_LABELS[key], key)
        self._lengths_box.activated.connect(self._on_layout)
        for text, w in (("Order", self._order), ("Lengths", self._lengths_box)):
            r = QHBoxLayout()
            lab = QLabel(text, self)
            lab.setFixedWidth(90)
            r.addWidget(lab)
            r.addWidget(w, 1)
            self._layout_section.body_layout().addLayout(r)
        self._show = NoWheelComboBox(self)
        for text, value in SHOW_ITEMS:
            self._show.addItem(text, value)
        self._show.setToolTip("Homologue groups draw each homologue: bivalents, trivalents … up to hexavalents")
        self._show.activated.connect(self._on_layout)
        self._drag = QCheckBox("Drag signals to correct their position", self)
        self._drag.toggled.connect(self._on_layout)
        self._layout_section.add_widget(self._drag)

        self._names_section = self._silico_section  # chromosome names belong to the assembly comparison
        self._silico_section.add_widget(_section_label("Chromosome names: FISH number to assembly sequence"))
        self._names_grid = QGridLayout()
        self._names_grid.setHorizontalSpacing(8)
        self._names_grid.setVerticalSpacing(4)
        self._silico_section.body_layout().addLayout(self._names_grid)
        side.addWidget(self._silico_section)
        side.addWidget(self._layout_section)
        side.addStretch(1)

        save = QPushButton("Save", self)
        menu = QMenu(save)
        menu.addAction("Picture (PNG)…", self._save_png)
        menu.addAction("Vector (SVG)…", self._save_svg)
        menu.addAction("Comparison table (TSV)…", self._save_tsv)
        save.setMenu(menu)
        fig = QPushButton("Add to figure", self)
        fig.setToolTip("Add this idiogram as a panel of the current figure")
        fig_menu = QMenu(fig)
        fig_menu.addAction("Colour", lambda: self._add_to_figure("compare"))
        fig_menu.addAction("Journal (Black and White, with Table)", lambda: self._add_to_figure("paper"))
        fig.setMenu(fig_menu)
        edit = QPushButton("Edit karyotype", self)
        edit.clicked.connect(self._edit_karyotype)
        bottom = QHBoxLayout()
        bottom.setSpacing(6)
        bottom.addWidget(edit, 1)
        bottom.addWidget(fig, 1)
        bottom.addWidget(save)
        side.addLayout(bottom)

        aside = QScrollArea(self)
        aside.setWidgetResizable(True)
        aside.setFrameShape(QFrame.Shape.NoFrame)
        aside.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        aside.setWidget(body)
        from repeatviewer.widgets.cytogenetics_style import fit_side_scroll

        fit_side_scroll(aside)

        self._side_stack = QStackedWidget(self)
        self._side_stack.addWidget(aside)
        self._side_stack.addWidget(self._build_edit_side())
        self._side_stack.addWidget(self._build_paper_side())
        self._side_stack.setMinimumWidth(320)  # the right column never clips its content
        self._side_stack.setMaximumWidth(SIDE_W)

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addWidget(main_host, 1)
        line = QFrame(self)
        line.setFrameShape(QFrame.Shape.VLine)
        line.setStyleSheet("color: palette(mid);")
        root.addWidget(line)
        root.addWidget(self._side_stack)

        # ---- top bar: Karyotype, Compare with, Show, "?"
        self._compare = NoWheelComboBox(self)
        self._compare.activated.connect(self._on_compare)
        for text, w, width in (("Karyotype", self._karyo, 200), ("Compare with", self._compare, 160),
                               ("Show", self._show, 160)):
            w.setFixedWidth(width)
            head.addWidget(QLabel(text, self))
            head.addWidget(w)
            head.addSpacing(14)
        self._show.setToolTip("Every homologue (1a, 1b, 1c) or one bar per chromosome")
        head.addStretch(1)
        head.addWidget(make_help_button(self, "Idiogram", HELP_TEXT))

    # ------------------------------------------------------------------ loading

    def refresh(self, folder: Path | None) -> None:
        self._folder = Path(folder) if folder else None
        self._dirty = True
        if self.isVisible():
            self._rebuild()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        if self._dirty:
            self._rebuild()

    def _rebuild(self) -> None:
        self._dirty = False
        if self._folder is None:
            return
        self._busy = True
        try:
            self._settings = idg.load_settings(self._folder)
            self._karyotypes = kt.list_karyotypes(self._folder)
            self._karyotype = next((k for k in self._karyotypes if k.name == self._settings.karyotype), None) \
                or (self._karyotypes[0] if self._karyotypes else None)
            self._colours = self._probe_colours()
            self._fill_controls()
            self._recompute()
        finally:
            self._busy = False

    def _probe_colours(self) -> dict[str, str]:
        colours: dict[str, str] = {}
        for image in ci.load_images(self._folder):
            for c in image_details(image).channels:
                if c.probe and c.probe not in colours and c.display.color.lower() not in ("#ffffff", "#5b74ff"):
                    colours[c.probe] = c.display.color
        lists = lists_mod.load_lists(self._folder)
        for e in lists.get("probes", []):
            name, dye = str(e.get("name") or ""), str(e.get("dye") or "")
            if name and name not in colours and dye:
                colours[name] = lists_mod.dye_color(lists, dye, "#8A8D93")
        return colours

    def _colour(self, probe: str, n: int) -> str:
        if probe not in self._colours:
            self._colours[probe] = FALLBACK_COLOURS[n % len(FALLBACK_COLOURS)]
        return self._colours[probe]

    def _genome_path(self) -> Path | None:
        return idg.resolve_genome(self._folder, self._settings.genome) if self._folder else None

    def _load_lengths(self) -> dict[str, int]:
        genome = self._genome_path()
        if genome is None:
            return {}
        key = ("lengths", str(genome), genome.stat().st_mtime)
        if key not in self._cache:
            QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
            try:
                self._cache[key] = idg.sequence_lengths(genome)
            except Exception:  # noqa: BLE001
                self._cache[key] = {}
            finally:
                QApplication.restoreOverrideCursor()
        return self._cache[key]

    def _all_probes(self) -> list[str]:
        probes: list[str] = []
        if self._fish is not None:
            probes += sorted(self._fish.cells_per_probe)
        for image in ci.load_images(self._folder) if self._folder else []:
            for c in image_details(image).channels:
                if c.probe and c.probe not in probes:
                    probes.append(c.probe)
        return probes

    def _recompute(self) -> None:
        k = self._karyotype
        s = self._settings
        if self._folder is None or k is None:
            self._rows = []
            self._canvas.set_rows([], [], {})
            self._canvas2.hide()
            self._fill_measured([])
            self._edit.karyotype = None
            self._edit.update()
            self._title.setText("Idiogram")
            self._note_fish.setText("Make a karyotype in the Signals tab (Edit karyotype) to see the idiogram.")
            self._note_asm.setText("")
            self._note_warn.setText("")
            return
        self._fish = idg.collect_fish(self._folder, k, per_homologue=s.homologues)
        self._lengths = self._load_lengths()
        genome = self._genome_path()
        hits: dict[str, dict] = {}
        probes = [p for p in self._all_probes() if p not in s.hidden_probes]
        if genome is not None and self._lengths:
            QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
            try:
                for probe in probes:
                    try:
                        hits[probe] = idg.features_for_probe(genome, list(self._lengths), probe, s.source(probe), self._cache)
                    except Exception:  # noqa: BLE001 - a broken annotation file must not stop the view
                        hits[probe] = {}
            finally:
                QApplication.restoreOverrideCursor()
        silico = idg.collect_silico(k, self._lengths, hits, min_array_kb=s.min_array_kb, flips=s.flips)
        self._rows = idg.build_rows(k, s, self._fish, silico, self._lengths)
        colours = {p: self._colour(p, n) for n, p in enumerate(self._all_probes())}
        for c in k.chromosomes:  # probes drawn as karyotype bands
            for b in c.bands:
                pr = kt.manual_probe(b.kind)
                if pr and pr not in colours:
                    colours[pr] = self._colour(pr, len(colours))
        other = next((x for x in self._karyotypes if x.name == s.compare_with and x is not k), None)
        for cv, kk in ((self._canvas, k), (self._canvas2, other)):
            cv.fit = False  # full-size chromosomes; a scroll bar when they do not fit
            cv.karyotype = kk
            cv.show_bands, cv.show_measured, cv.show_assembly = s.show_bands, s.show_measured, s.show_assembly
            cv.assembly_active = genome is not None
            cv.drag_enabled = s.drag
            cv.reorder_enabled = s.order == "custom" and cv is self._canvas
            cv.title = f"{kk.name}  {kk.formula}" if (kk is not None and other is not None) else ""
        self._canvas.set_rows(self._rows, idg.relative_heights(self._rows, s.lengths), colours)
        if other is not None:
            fish2 = idg.collect_fish(self._folder, other, per_homologue=s.homologues)
            silico2 = idg.collect_silico(other, self._lengths, hits, min_array_kb=s.min_array_kb, flips=s.flips)
            rows2 = idg.build_rows(other, s, fish2, silico2, self._lengths)
            self._canvas2.set_rows(rows2, idg.relative_heights(rows2, s.lengths), colours)
            self._canvas2.show()
        else:
            self._canvas2.hide()
        self._fill_measured(self._unmatched())
        self._edit.karyotype = k
        self._edit.groups = s.homologues
        self._edit.probe_colours = colours
        fish: dict[str, list] = {}
        for row in self._rows:
            for band in row.fish:
                fish.setdefault(row.chromosome.name, []).append((band.fraction, band.probe, band.arm == "cen"))
        self._edit.fish = fish
        if not 0 <= self._edit.index < len(k.chromosomes):
            self._edit.index = 0 if k.chromosomes else -1
        self._edit.update()
        self._fill_kinds()
        self._sync_edit_side()
        self._refresh_pair()
        self._refresh_paper()
        self._probe_colour_map = colours
        f = self._fish
        self._note_fish.setText(f"Left bar: signals from {f.homologues} marked homologues in {f.cells} "
                                f"cell{'s' if f.cells != 1 else ''}; the number is how often the signal was seen.")
        if genome is None:
            self._note_asm.setText("Right bar: choose a genome under In silico to compare with the assembly.")
        else:
            kinds = sorted({idg.SOURCE_LABELS[s.source(p).kind].lower() for p in probes})
            self._note_asm.setText(f"Right bar: scaled to the assembly chromosome ({genome.name}); marks are arrays of at "
                                   f"least {s.min_array_kb:g} kb from {', '.join(kinds) or 'no source'}.")
        warn = []
        if f.unplaced:
            warn.append(f"{f.unplaced} signal{'s' if f.unplaced != 1 else ''} not shown: the chromosome is not measured "
                        "(use Measure in the Signals tab) or the arm is not set.")
        unlinked = [r.chromosome.name for r in self._rows if not r.sequence_length and genome is not None]
        if unlinked:
            warn.append("No assembly sequence for chromosome " + ", ".join(dict.fromkeys(unlinked)) + " (see Chromosome names).")
        self._note_warn.setText(" ".join(warn))
        self._fill_names()

    # ------------------------------------------------------------------ controls

    def _clear(self, layout) -> None:
        clear_layout(layout)

    def _fill_legend(self, probes: list[str], colours: dict[str, str]) -> None:
        return  # the Probes box is the legend
        self._clear(self._legend)
        for text, css in [(p, f"background:{colours.get(p, '#8A8D93')}; border-radius:2px;") for p in probes] + [
            ("In situ (FISH)", "background:#FFFFFF; border:1px solid #8A8D93; border-radius:2px;"),
            ("In silico (assembly)", "background:#EEF0F4; border:1px dashed #8A8D93; border-radius:2px;"),
        ]:
            chip = QLabel(self)
            chip.setFixedSize(14, 10)
            chip.setStyleSheet(css)
            row = QHBoxLayout()
            row.setSpacing(6)
            row.addWidget(chip)
            row.addWidget(QLabel(text, self))
            self._legend.addLayout(row)

    def _fill_controls(self) -> None:
        s = self._settings
        # probes
        self._clear(self._probe_lay)
        fish_cells = self._fish.cells_per_probe if self._fish else {}
        if self._karyotype is not None and self._folder is not None:
            fish_cells = idg.collect_fish(self._folder, self._karyotype).cells_per_probe
        lists = lists_mod.load_lists(self._folder) if self._folder else {}
        for n_, probe in enumerate(self._all_probes_with(fish_cells)):
            n = fish_cells.get(probe, 0)
            row = QHBoxLayout()
            row.setSpacing(8)
            box = QCheckBox(probe, self._probe_host)
            pm = QPixmap(14, 10)
            pm.fill(QColor(self._colour(probe, n_)))
            box.setIcon(QIcon(pm))
            box.setChecked(probe not in s.hidden_probes)
            box.toggled.connect(lambda on, p=probe: self._on_probe(p, on))
            dye = str((lists_mod.entry(lists, "probes", probe) or {}).get("dye") or "") if lists else ""
            info = QLabel(f"{dye}   {n} cell{'s' if n != 1 else ''}".strip(), self._probe_host)
            info.setStyleSheet("color: palette(dark);")
            row.addWidget(box, 1)
            row.addWidget(info)
            self._probe_lay.addLayout(row)
        if not self._probe_lay.count():
            self._probe_lay.addWidget(_section_label("No probes yet. Set the probe of each channel in the Image tab."))
        self._share.setCurrentIndex(max(0, self._share.findData(s.min_share)))
        # genome
        self._genome.clear()
        self._genome.addItem(GENOME_NONE, "")
        for path in idg.genome_candidates(self._folder):
            self._genome.addItem(path.name, path.name)
        if s.genome and self._genome.findData(s.genome) < 0:
            self._genome.addItem(Path(s.genome).name, s.genome)
        self._genome.addItem(GENOME_BROWSE, "\x00browse")
        self._genome.setCurrentIndex(max(0, self._genome.findData(s.genome)))
        self._fill_sources()
        self._min_kb.setValue(s.min_array_kb)
        self._tol.setValue(s.tolerance * 100)
        # layout
        self._karyo.clear()
        for k in self._karyotypes:
            self._karyo.addItem(f"{k.name} {k.formula}", k.name)
        if self._karyotype is not None:
            self._karyo.setCurrentIndex(max(0, self._karyo.findData(self._karyotype.name)))
        self._order.setCurrentIndex(max(0, self._order.findData(s.order)))
        self._lengths_box.setCurrentIndex(max(0, self._lengths_box.findData(s.lengths)))
        self._show.setCurrentIndex(max(0, self._show.findData(bool(s.homologues))))
        self._compare.blockSignals(True)
        self._compare.clear()
        self._compare.addItem("None", "")
        for k in self._karyotypes:
            if self._karyotype is None or k.name != self._karyotype.name:
                self._compare.addItem(f"{k.name} {k.formula}", k.name)
        self._compare.setCurrentIndex(max(0, self._compare.findData(s.compare_with)))
        self._compare.blockSignals(False)
        for key, box in self._layer_boxes.items():
            box.blockSignals(True)
            box.setChecked(bool(getattr(s, key)))
            box.blockSignals(False)
        self._layer_boxes["show_assembly"].setEnabled(bool(s.genome))
        self._drag.blockSignals(True)
        self._drag.setChecked(s.drag)
        self._drag.blockSignals(False)

    def _all_probes_with(self, fish_cells: dict[str, int]) -> list[str]:
        probes = sorted(fish_cells)
        for image in ci.load_images(self._folder) if self._folder else []:
            for c in image_details(image).channels:
                if c.probe and c.probe not in probes:
                    probes.append(c.probe)
        return probes

    def _fill_sources(self) -> None:
        self._clear(self._source_lay)
        s = self._settings
        enabled = bool(s.genome)
        for probe in [p for p in self._all_probes() if p not in s.hidden_probes]:
            src = s.source(probe)
            kind = NoWheelComboBox(self._source_host)
            for key in idg.SOURCE_KINDS:
                kind.addItem(idg.SOURCE_LABELS[key], key)
            kind.setCurrentIndex(max(0, kind.findData(src.kind)))
            kind.setEnabled(enabled)
            kind.activated.connect(lambda _i, p=probe, cb=kind: self._on_source_kind(p, cb.currentData()))
            r = QHBoxLayout()
            lab = QLabel(f"{probe} from", self._source_host)
            lab.setFixedWidth(90)
            lab.setWordWrap(True)
            r.addWidget(lab)
            r.addWidget(kind, 1)
            self._source_lay.addLayout(r)
            match = QLineEdit(src.match, self._source_host)
            match.setPlaceholderText(f"name contains: {probe}")
            match.setToolTip("Features whose name or type contain this text count as hits of the probe")
            match.setEnabled(enabled and src.kind != "none")
            match.editingFinished.connect(lambda p=probe, w=match: self._on_source_match(p, w.text()))
            r2 = QHBoxLayout()
            r2.addSpacing(98)
            r2.addWidget(match, 1)
            self._source_lay.addLayout(r2)
            if src.kind == "file":
                pick = QPushButton(Path(src.file).name if src.file else "Choose file…", self._source_host)
                pick.setToolTip(src.file or "GFF or BED file with the arrays")
                pick.clicked.connect(lambda _c=False, p=probe: self._choose_source_file(p))
                r3 = QHBoxLayout()
                r3.addSpacing(98)
                r3.addWidget(pick, 1)
                self._source_lay.addLayout(r3)
        if not enabled:
            self._source_lay.addWidget(_section_label("Choose a genome FASTA of this folder to compare with the assembly."))

    def _fill_names(self) -> None:
        self._clear(self._names_grid)
        k = self._karyotype
        if k is None:
            return
        seqs = sorted(self._lengths)
        self._names_grid.addWidget(_section_label("FISH"), 0, 0)
        self._names_grid.addWidget(_section_label("Assembly"), 0, 1)
        self._names_grid.addWidget(_section_label("Flip"), 0, 2)
        for r, chrom in enumerate(k.chromosomes, start=1):
            self._names_grid.addWidget(QLabel(chrom.name, self), r, 0)
            combo = NoWheelComboBox(self)
            combo.setEditable(True)
            combo.addItems([""] + seqs)
            combo.setCurrentText(chrom.assembly_seq)
            combo.activated.connect(lambda _i, c=chrom.name, cb=combo: self._on_link(c, cb.currentText()))
            combo.lineEdit().editingFinished.connect(lambda c=chrom.name, cb=combo: self._on_link(c, cb.currentText()))
            self._names_grid.addWidget(combo, r, 1)
            flip = QCheckBox(self)
            flip.setToolTip("The assembly sequence starts at the q telomere")
            flip.setChecked(self._settings.flips.get(chrom.name, False))
            flip.toggled.connect(lambda on, c=chrom.name: self._on_flip(c, on))
            self._names_grid.addWidget(flip, r, 2)
        self._names_grid.setColumnStretch(1, 1)

    # ------------------------------------------------------------------ changes

    def _store(self, *, controls: bool = False) -> None:
        if self._folder is None or self._busy:
            return
        idg.save_settings(self._folder, self._settings)
        self._busy = True
        try:
            if controls:
                self._fill_controls()
            self._recompute()
        finally:
            self._busy = False

    def _on_probe(self, probe: str, on: bool) -> None:
        hidden = set(self._settings.hidden_probes)
        (hidden.discard if on else hidden.add)(probe)
        self._settings.hidden_probes = sorted(hidden)
        self._store()
        self._fill_sources()

    def _on_share(self, index: int) -> None:
        self._settings.min_share = float(self._share.itemData(index))
        self._store()

    def _on_genome(self, index: int) -> None:
        data = self._genome.itemData(index)
        if data == "\x00browse":
            path, _ = QFileDialog.getOpenFileName(self, "Genome FASTA", str(self._folder or ""),
                                                  "FASTA (*.fa *.fasta *.fna *.fas);;All files (*)")
            if not path:
                self._genome.setCurrentIndex(max(0, self._genome.findData(self._settings.genome)))
                return
            data = path
        self._settings.genome = data or ""
        self._store(controls=True)

    def _on_source_kind(self, probe: str, kind: str) -> None:
        self._settings.source(probe).kind = kind
        if kind == "file" and not self._settings.source(probe).file:
            self._choose_source_file(probe)
            return
        self._store()
        self._fill_sources()

    def _on_source_match(self, probe: str, text: str) -> None:
        if self._settings.source(probe).match != text.strip():
            self._settings.source(probe).match = text.strip()
            self._store()

    def _choose_source_file(self, probe: str) -> None:
        path, _ = QFileDialog.getOpenFileName(self, f"Arrays of {probe}", str(self._folder or ""),
                                              "Annotations (*.gff *.gff3 *.gtf *.bed);;All files (*)")
        if path:
            self._settings.source(probe).file = path
            self._settings.source(probe).kind = "file"
        self._store()
        self._fill_sources()

    def _on_numbers(self) -> None:
        self._settings.min_array_kb = float(self._min_kb.value())
        self._settings.tolerance = float(self._tol.value()) / 100.0
        self._store()

    def _on_karyotype(self, index: int) -> None:
        name = str(self._karyo.itemData(index) or "")
        if not name or name == self._settings.karyotype:
            return
        self._settings = idg.switch_karyotype(self._settings, name)  # each karyotype keeps its own view
        if self._settings.compare_with == name:
            self._settings.compare_with = ""
        self._karyotype = next((k for k in self._karyotypes if k.name == name), self._karyotype)
        self._store(controls=True)

    def _on_compare(self, index: int) -> None:
        self._settings.compare_with = str(self._compare.itemData(index) or "")
        self._store()

    def _on_layer(self, key: str, on: bool) -> None:
        setattr(self._settings, key, on)
        for cv in (self._canvas, self._canvas2):
            setattr(cv, key, on)
            cv.update()
        if self._folder is not None and not self._busy:
            idg.save_settings(self._folder, self._settings)

    # ------------------------------------------------------------------ measured, not in the karyotype

    def _unmatched(self) -> list[tuple[idg.ChromosomeRow, idg.FishBand]]:
        out = []
        for row in self._rows:
            view = self._canvas._view(row)
            for band in row.fish:
                near = any(kt.manual_probe(b.kind) == band.probe and b.start - 0.03 <= band.fraction <= b.end + 0.03
                           for b in view.bands)
                if not near:
                    out.append((row, band))
        return out

    def _fill_measured(self, items: list) -> None:
        self._clear(self._measured_grid)
        self._unmatched_items = items
        self._measured_section.setVisible(bool(items))
        for r, (row, band) in enumerate(items[:12]):
            dot = QLabel(self)
            pm = QPixmap(10, 10)
            pm.fill(Qt.GlobalColor.transparent)
            pt = QPainter(pm)
            pt.setRenderHint(QPainter.RenderHint.Antialiasing)
            pt.setPen(Qt.PenStyle.NoPen)
            pt.setBrush(self._canvas.colour(band.probe))
            pt.drawEllipse(QRectF(0, 0, 10, 10))
            pt.end()
            dot.setPixmap(pm)
            where = f"{band.arm} {band.position}" if band.arm in ("p", "q") else band.position
            for c, w in enumerate((QLabel(row.name, self), dot, QLabel(band.probe, self), QLabel(where, self),
                                   QLabel(f"{band.seen}/{band.total}", self))):
                if c >= 3:
                    w.setStyleSheet("color: palette(dark);")
                self._measured_grid.addWidget(w, r, c)
        self._measured_grid.setColumnStretch(2, 1)
        if len(items) > 12:
            more = QLabel(f"+ {len(items) - 12} more", self)
            more.setStyleSheet("color: palette(dark);")
            self._measured_grid.addWidget(more, 12, 0, 1, 5)

    def _add_measured_to_karyotype(self) -> None:
        k = self._karyotype
        items = getattr(self, "_unmatched_items", [])
        if k is None or not items or self._folder is None:
            return
        suffixes = k.homologue_suffixes()
        for row, band in items:
            chrom = k.chromosome(row.chromosome.name)
            if chrom is None:
                continue
            i = suffixes.index(row.homologue) if row.homologue in suffixes else -1
            owner = kt.target(chrom, i)
            f = band.fraction
            owner.bands.append(kt.Band(kt.SIGNAL_KIND_PREFIX + band.probe, max(0.0, f - 0.015), min(1.0, f + 0.015)))
        kt.save_karyotype(self._folder, k, old_name=k.name)
        from repeatviewer.widgets.cyto_images_ui import notify_changed

        notify_changed(self._folder)
        self.karyotypesChanged.emit()
        self._rebuild()

    def _on_layout(self, *_args) -> None:
        s = self._settings
        s.order = str(self._order.currentData() or "karyotype")
        s.lengths = str(self._lengths_box.currentData() or "karyotype")
        s.homologues = bool(self._show.currentData())
        s.drag = self._drag.isChecked()
        if s.order == "custom" and not s.custom_order and self._karyotype is not None:
            s.custom_order = [c.name for c in self._karyotype.chromosomes]
        self._store()

    def _on_order_dragged(self, names: list) -> None:
        self._settings.custom_order = [str(n) for n in names]
        self._store()

    def _on_link(self, chrom: str, seq: str) -> None:
        k = self._karyotype
        if k is None or self._folder is None:
            return
        c = k.chromosome(chrom)
        seq = seq.strip()
        if c is None or c.assembly_seq == seq:
            return
        c.assembly_seq = seq
        kt.save_karyotype(self._folder, k)
        self.karyotypesChanged.emit()
        self._store()

    def _on_flip(self, chrom: str, on: bool) -> None:
        self._settings.flips[chrom] = bool(on)
        self._store()

    def _on_band_moved(self, band: idg.FishBand, fraction: float) -> None:
        k = self._karyotype
        if k is None or self._folder is None:
            return
        chrom = k.chromosome(band.chromosome)
        if chrom is None:
            return
        n = idg.move_band(self._folder, band, fraction, chrom.centromere_fraction)
        if n:
            from repeatviewer.widgets.cyto_images_ui import notify_changed

            notify_changed(self._folder)  # refresh() rebuilds this tab and the Signals tab

    def _edit_karyotype(self) -> None:
        if self._folder is None:
            return
        from repeatviewer.views.cyto_karyotype_dialog import KaryotypeDialog

        dialog = KaryotypeDialog(self._folder, self._karyotype.name if self._karyotype else "", self.window())
        dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)

        def done() -> None:
            self._settings.karyotype = dialog.saved_name() or self._settings.karyotype
            idg.save_settings(self._folder, self._settings)
            self.karyotypesChanged.emit()
            self._rebuild()

        dialog.accepted.connect(done)
        dialog.show()  # not modal: the Files panel stays clickable (Lengths from assembly)

    # ------------------------------------------------------------------ saving

    def _default_name(self, suffix: str) -> str:
        base = self._karyotype.name if self._karyotype else "idiogram"
        return str((self._folder or Path.home()) / f"{kt._safe_name(base)}_idiogram{suffix}")

    def _save_png(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save idiogram picture", self._default_name(".png"), "PNG (*.png)")
        if path and not self._canvas.to_image(2.0).save(path, "PNG"):
            QMessageBox.warning(self, "Save", f"Could not write {path}")

    def _save_svg(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save idiogram as vector", self._default_name(".svg"), "SVG (*.svg)")
        if path:
            try:
                self._canvas.to_svg(Path(path))
            except Exception as exc:  # noqa: BLE001
                QMessageBox.warning(self, "Save", str(exc))

    def _save_tsv(self) -> None:
        path, _ = QFileDialog.getSaveFileName(self, "Save comparison table", self._default_name(".tsv"), "TSV (*.tsv)")
        if path:
            try:
                Path(path).write_text(idg.comparison_tsv(self._rows), encoding="utf-8")
            except OSError as exc:
                QMessageBox.warning(self, "Save", str(exc))

    # ------------------------------------------------------------------ Edit view

    def _build_edit_side(self) -> QWidget:
        body = QWidget(self)
        body.setObjectName("cytoAside")
        side = QVBoxLayout(body)
        side.setContentsMargins(16, 16, 16, 16)
        side.setSpacing(10)
        self._e_title = QLabel(body)
        f = self._e_title.font()
        f.setBold(True)
        self._e_title.setFont(f)
        side.addWidget(self._e_title)
        self._e_len, self._e_p, self._e_q = QLineEdit(body), QLineEdit(body), QLineEdit(body)
        self._e_units: list[QLabel] = []
        for text, w in (("Length", self._e_len), ("Short arm (p)", self._e_p), ("Long arm (q)", self._e_q)):
            w.setFixedWidth(110)
            w.setAlignment(Qt.AlignmentFlag.AlignRight)
            w.editingFinished.connect(self._on_edit_numbers)
            r = QHBoxLayout()
            r.addWidget(QLabel(text, body), 1)
            r.addWidget(w)
            side.addLayout(r)
        self._e_info = QLabel(body)
        self._e_info.setWordWrap(True)
        self._e_info.setStyleSheet("color: gray;")
        side.addWidget(self._e_info)
        side.addSpacing(8)
        self._e_bands_title = QLabel(body)
        self._e_bands_title.setFont(f)
        side.addWidget(self._e_bands_title)
        self._e_bands = QListWidget(body)
        self._e_bands.setFrameShape(QFrame.Shape.NoFrame)
        self._e_bands.setStyleSheet("QListWidget { background: transparent; }")
        self._e_bands.currentRowChanged.connect(self._on_edit_band_row)
        side.addWidget(self._e_bands, 1)
        self._e_remove = QPushButton("Remove band", body)
        self._e_remove.clicked.connect(self._remove_edit_band)
        r = QHBoxLayout()
        r.addWidget(self._e_remove)
        r.addStretch(1)
        side.addLayout(r)
        return body

    def _set_mode(self, mode: str) -> None:
        self._mode = mode
        for key, b in self._mode_btns.items():
            b.setChecked(key == mode)
        edit = mode == "edit"
        self._stack.setCurrentIndex({"compare": 0, "edit": 1, "pair": 2, "paper": 3}[mode])
        self._side_stack.setCurrentIndex({"edit": 1, "paper": 2}.get(mode, 0))
        self._side_stack.setVisible(mode != "pair")
        self._tools_host.setVisible(edit)
        self._pair_host.setVisible(mode == "pair")
        for w in (self._note_fish, self._note_asm, self._note_warn):
            w.setVisible(mode == "compare")
        self._sync_edit_side()
        self._refresh_pair()
        self._refresh_paper()

    def _set_tool(self, tool: str) -> None:
        self._edit.tool = tool
        for key, b in self._tool_btns.items():
            b.setChecked(key == tool)
        self._fill_kinds()
        self._edit.update()

    def _fill_kinds(self) -> None:
        """Band kinds for Paint band, probes for Place signal."""
        box = self._kind
        box.blockSignals(True)
        box.clear()
        if self._edit.tool == "place":
            probes = self._all_probes() if self._folder else []
            for p in probes:
                box.addItem(p, p)
            if not probes:
                box.addItem("No probes yet", "")
            box.setCurrentIndex(max(0, box.findData(self._edit.kind)))
            self._edit.kind = str(box.currentData() or "")
            self._edit.kind_colour = ""
        else:
            for text, kind in PAINT_KINDS:
                box.addItem(text, kind)
            k = self._karyotype
            own = {(b.kind, b.color) for c in (k.chromosomes if k else []) for b in c.bands
                   if b.color and not kt.manual_probe(b.kind)}
            for kind, colour in sorted(own):
                box.addItem(kind, (kind, colour))
            box.addItem(OWN_COLOUR, "\x00own")
            current = (self._edit.kind, self._edit.kind_colour) if self._edit.kind_colour else self._edit.kind
            i = box.findData(current)
            box.setCurrentIndex(max(0, i))
            if i < 0:
                self._edit.kind, self._edit.kind_colour = PAINT_KINDS[0][1], ""
        box.setEnabled(self._edit.tool != "select")
        box.blockSignals(False)

    def _on_kind(self, index: int) -> None:
        data = self._kind.itemData(index)
        if data == "\x00own":
            name, ok = QInputDialog.getText(self, "Own band kind", "Name of the band kind:")
            if not ok or not name.strip():
                self._fill_kinds()
                return
            colour = QColorDialog.getColor(QColor("#E8A33D"), self, f"Colour of {name.strip()}")
            if not colour.isValid():
                self._fill_kinds()
                return
            self._edit.kind, self._edit.kind_colour = name.strip(), colour.name()
            self._fill_kinds()
            return
        if isinstance(data, (tuple, list)):
            self._edit.kind, self._edit.kind_colour = str(data[0]), str(data[1])
        else:
            self._edit.kind, self._edit.kind_colour = str(data or ""), ""

    def _edit_chrom(self) -> kt.KaryoChromosome | None:
        k = self._karyotype
        i = self._edit.index
        return k.chromosomes[i] if k is not None and 0 <= i < len(k.chromosomes) else None

    def _on_edit_selected(self, _i: int) -> None:
        self._sync_edit_side()

    def _sync_edit_side(self) -> None:
        c = self._edit_chrom()
        k = self._karyotype
        on = c is not None
        for w in (self._e_len, self._e_p, self._e_q, self._e_remove):
            w.setEnabled(on)
        self._e_bands.blockSignals(True)
        self._e_bands.clear()
        if c is None:
            self._e_title.setText("No chromosome")
            self._e_info.setText("Make a karyotype first (Edit karyotype under Layout in Compare).")
            self._e_bands_title.clear()
            for w in (self._e_len, self._e_p, self._e_q):
                w.clear()
            self._e_bands.blockSignals(False)
            return
        unit = k.unit
        p, q = c.arms
        self._e_title.setText(f"Chromosome {c.name}")
        self._e_len.setText(f"{c.length:g} {unit}")
        self._e_p.setText(f"{round(p, 3):g} {unit}")
        self._e_q.setText(f"{round(q, 3):g} {unit}")
        ci_text = f"Centromere index {c.centromere_index:g}, type {c.type}. " if c.length else ""
        self._e_info.setText(ci_text + "Every change is saved to the karyotype and shows in Signals and Figure.")
        self._e_bands_title.setText(f"Bands and signals on {c.name}")
        colours = self._edit.probe_colours
        for row in self._rows:
            if row.chromosome.name != c.name:
                continue
            for band in row.fish:
                hom = f" {row.name}" if row.homologue else ""
                item = QListWidgetItem(self._swatch(colours.get(band.probe, "#8A8D93")),
                                       f"{band.probe}{hom} (FISH, {band.label})")
                item.setData(Qt.ItemDataRole.UserRole, -1)
                item.setToolTip("Marked in the Signals tab; move it there or drag it in Compare")
                self._e_bands.addItem(item)
        f = c.centromere_fraction
        for n, b in enumerate(c.bands):
            probe = kt.manual_probe(b.kind)
            if probe:
                text = f"{probe} (placed by hand, {kt.fraction_to_arm_text((b.start + b.end) / 2, f)})"
            else:
                text = f"{b.kind}, {kt.fraction_to_arm_text(b.start, f)} to {kt.fraction_to_arm_text(b.end, f)}"
            item = QListWidgetItem(self._swatch(self._edit.band_colour(b).name(), hatched=b.kind.startswith("Secondary")), text)
            item.setData(Qt.ItemDataRole.UserRole, n)
            self._e_bands.addItem(item)
        for r in range(self._e_bands.count()):
            if self._edit.band_index >= 0 and self._e_bands.item(r).data(Qt.ItemDataRole.UserRole) == self._edit.band_index:
                self._e_bands.setCurrentRow(r)
        if not self._e_bands.count():
            self._e_bands.addItem("None yet. Use Paint band or Place signal.")
        self._e_bands.blockSignals(False)
        self._e_remove.setEnabled(0 <= self._edit.band_index < len(c.bands))

    @staticmethod
    def _swatch(colour: str, hatched: bool = False):
        pix = QPixmap(14, 10)
        pix.fill(QColor("#FFFFFF") if hatched else QColor(colour))
        if hatched:
            p = QPainter(pix)
            p.fillRect(pix.rect(), QBrush(QColor("#45484D"), Qt.BrushStyle.BDiagPattern))
            p.end()
        return QIcon(pix)

    def _on_edit_band_row(self, row: int) -> None:
        item = self._e_bands.item(row)
        n = item.data(Qt.ItemDataRole.UserRole) if item is not None else -1
        self._edit.band_index = n if isinstance(n, int) else -1
        c = self._edit_chrom()
        self._e_remove.setEnabled(c is not None and 0 <= self._edit.band_index < len(c.bands))
        self._edit.update()

    def _on_edit_numbers(self) -> None:
        c = self._edit_chrom()
        if c is None:
            return

        def num(w: QLineEdit) -> float | None:
            try:
                return float(w.text().split()[0].replace(",", ".")) if w.text().split() else None
            except ValueError:
                return None

        length, p, q = num(self._e_len), num(self._e_p), num(self._e_q)
        old_p, old_q = c.arms
        if p is not None and q is not None and (abs(p - old_p) > 1e-6 or abs(q - old_q) > 1e-6) and p + q > 0:
            c.set_arms(p, q)
        elif length is not None and length > 0 and abs(length - c.length) > 1e-9:
            c.length = round(length, 4)
        else:
            self._sync_edit_side()
            return
        self._on_edit_changed()

    def _remove_edit_band(self) -> None:
        c = self._edit_chrom()
        if c is None or not 0 <= self._edit.band_index < len(c.bands):
            return
        del c.bands[self._edit.band_index]
        self._edit.band_index = -1
        self._on_edit_changed()

    def _on_edit_changed(self) -> None:
        """Save the karyotype after each finished edit; the other tabs pick it up."""
        k = self._karyotype
        if k is None or self._folder is None:
            return
        try:
            kt.save_karyotype(self._folder, k)
        except OSError as exc:
            QMessageBox.warning(self, "Idiogram", f"Could not save the karyotype:\n{exc}")
            return
        self._edit.update()
        self._canvas.update()
        self._sync_edit_side()
        self.karyotypesChanged.emit()

    # ------------------------------------------------------------------ click to jump

    def _on_band_clicked(self, i: int, j: int, pos: QPoint) -> None:
        if not 0 <= i < len(self._rows) or not 0 <= j < len(self._rows[i].fish) or self._folder is None:
            return
        row, band = self._rows[i], self._rows[i].fish[j]
        where = f"{band.arm} {band.position}" if band.arm in ("p", "q") else band.position
        pop = _Popover(self)
        pop.title(f"{band.probe} on chromosome {row.name}, {where}")
        pop.note(f"Seen on {band.seen} of {band.total} marked homologues in {band.cells} "
                 f"cell{'s' if band.cells != 1 else ''}. Click a cell to open it in Signals with the signal selected.")
        first: dict[str, str] = {}
        for image_id, signal_id in band.refs:
            first.setdefault(image_id, signal_id)
        images = {im.id: im for im in ci.load_images(self._folder)}
        grid = QGridLayout()
        grid.setSpacing(6)
        for n, (image_id, signal_id) in enumerate(list(first.items())[:21]):
            image = images.get(image_id)
            if image is None:
                continue
            b = QToolButton(pop)
            b.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextUnderIcon)
            b.setAutoRaise(True)
            b.setText(_cell_name(image))
            b.setToolTip(image.original_name or image.file)
            pix = QPixmap(str(ci.image_file_path(self._folder, image)))
            if not pix.isNull():
                b.setIcon(QIcon(pix.scaled(64, 52, Qt.AspectRatioMode.KeepAspectRatioByExpanding,
                                           Qt.TransformationMode.SmoothTransformation).copy(0, 0, 64, 52)))
                b.setIconSize(QSize(64, 52))
            b.clicked.connect(lambda _c=False, im=image_id, sg=signal_id: (pop.close(), self.openSignalRequested.emit(im, sg)))
            grid.addWidget(b, n // 7, n % 7)
        pop.lay.addLayout(grid)
        if len(first) > 21:
            pop.note(f"and {len(first) - 21} more cells.")
        if j not in row.matched_fish and row.sequence_length:
            pop.note("No array of this probe nearby in the assembly.")
        pop.popup_at(pos)

    def _on_array_clicked(self, i: int, j: int, pos: QPoint) -> None:
        if not 0 <= i < len(self._rows) or not 0 <= j < len(self._rows[i].silico):
            return
        row, array = self._rows[i], self._rows[i].silico[j]
        genome = self._genome_path()
        start, end = idg.array_span_bp(array, row.sequence_length, self._settings.flips.get(row.chromosome.name, False))
        pop = _Popover(self)
        pop.title(f"{array.probe} array on {row.sequence}, {start / 1e6:.1f} to {end / 1e6:.1f} Mb")
        size = (end - start + 1) / 1e6
        text = f"{size:.2f} Mb, {array.hits} hit{'s' if array.hits != 1 else ''}."
        if j not in row.matched_silico:
            text += " No FISH signal here: maybe too short or too diverged to see."
        pop.note(text)
        buttons = QHBoxLayout()
        open_btn = QPushButton("Open in sequence viewer", pop)
        open_btn.setEnabled(genome is not None)
        open_btn.clicked.connect(lambda: (pop.close(), self.openRequested.emit(genome)))
        extract = QPushButton("Extract array as FASTA", pop)
        extract.setEnabled(genome is not None)
        extract.clicked.connect(lambda: (pop.close(), self._extract_array(genome, row.sequence, start, end)))
        buttons.addWidget(open_btn)
        buttons.addWidget(extract)
        buttons.addStretch(1)
        pop.lay.addLayout(buttons)
        pop.popup_at(pos)

    def _extract_array(self, genome: Path | None, seq: str, start: int, end: int) -> None:
        if genome is None:
            return
        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        out = confirm_existing_file_path(self, idg.region_output_path(genome, seq, start, end),
                                         title="Extract array as FASTA")
        if out is None:
            return
        QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
        try:
            idg.extract_region(genome, seq, start, end, out)
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, "Extract array as FASTA", str(exc))
            return
        finally:
            QApplication.restoreOverrideCursor()
        self.openRequested.emit(out)

    # ------------------------------------------------------------------ Two karyotypes

    def _pooled(self, k: kt.Karyotype) -> list[idg.FishBand]:
        s = self._settings
        return idg.visible_bands(idg.collect_fish(self._folder, k), set(s.hidden_probes), s.min_share)

    def _refresh_pair(self) -> None:
        if self._mode != "pair" or self._folder is None:
            return
        s = self._settings
        names = [k.name for k in self._karyotypes]
        top = self._karyotype
        bottom = next((k for k in self._karyotypes if k.name == s.compare_with and k is not top), None) \
            or next((k for k in self._karyotypes if k is not top), None)
        for box, chosen in ((self._top_box, top), (self._bottom_box, bottom)):
            box.blockSignals(True)
            box.clear()
            for k in self._karyotypes:
                box.addItem(f"{k.species or k.name} ({k.formula})", k.name)
            box.setCurrentIndex(names.index(chosen.name) if chosen is not None else -1)
            box.blockSignals(False)
        self._link.blockSignals(True)
        self._link.setChecked(s.link_probes)
        self._link.blockSignals(False)
        self._pair.colours = self._edit.probe_colours
        self._pair.link = s.link_probes
        if top is not None and bottom is not None:
            self._pair.top = (top, self._pooled(top))
            self._pair.bottom = (bottom, self._pooled(bottom))
        else:
            self._pair.top = self._pair.bottom = None
        self._pair.update()

    def _on_pair_choice(self, *_args) -> None:
        s = self._settings
        top = str(self._top_box.currentData() or "")
        bottom = str(self._bottom_box.currentData() or "")
        s.link_probes = self._link.isChecked()
        s.compare_with = bottom
        if top and top != s.karyotype:
            s.karyotype = top
            self._karyotype = next((k for k in self._karyotypes if k.name == top), self._karyotype)
            self._store(controls=True)
        else:
            self._store()
        self._refresh_pair()

    # ------------------------------------------------------------------ Paper

    def _build_paper_side(self) -> QWidget:
        from PySide6.QtWidgets import QFontComboBox

        body = QWidget(self)
        body.setObjectName("cytoAside")
        side = QVBoxLayout(body)
        side.setContentsMargins(16, 16, 16, 16)
        side.setSpacing(10)
        title = QLabel("Paper style", body)
        f = title.font()
        f.setBold(True)
        title.setFont(f)
        side.addWidget(title)
        self._paper_boxes: dict[str, QCheckBox] = {}
        for key, text in (("outlines", "Black outlines, white chromosomes"), ("formula", "Formula line"),
                          ("types", "Type under each chromosome"), ("table", "Table of positions"),
                          ("scale", "Scale bar")):
            box = QCheckBox(text, body)
            box.toggled.connect(self._on_paper_option)
            self._paper_boxes[key] = box
            side.addWidget(box)
        self._paper_font = QFontComboBox(body)
        self._paper_font.currentFontChanged.connect(self._on_paper_option)
        self._paper_size = QDoubleSpinBox(body)
        self._paper_size.setRange(5, 14)
        self._paper_size.setDecimals(1)
        self._paper_size.setSingleStep(0.5)
        self._paper_size.setSuffix(" pt")
        self._paper_size.valueChanged.connect(self._on_paper_option)
        for text, w in (("Font", self._paper_font), ("Size", self._paper_size)):
            r = QHBoxLayout()
            lab = QLabel(text, body)
            lab.setFixedWidth(50)
            r.addWidget(lab)
            r.addWidget(w, 1)
            side.addLayout(r)
        side.addSpacing(6)
        save = QPushButton("Save SVG · PDF · TIFF", body)
        menu = QMenu(save)
        menu.addAction("Vector (SVG)…", lambda: self._save_paper("svg"))
        menu.addAction("PDF…", lambda: self._save_paper("pdf"))
        menu.addAction("TIFF, 600 dpi…", lambda: self._save_paper("tif"))
        save.setMenu(menu)
        fig = QPushButton("Add to figure", body)
        fig.setToolTip("Add the paper style idiogram as a panel of the current figure")
        fig.clicked.connect(lambda: self._add_to_figure("paper"))
        for b in (save, fig):
            r = QHBoxLayout()
            r.addWidget(b)
            r.addStretch(1)
            side.addLayout(r)
        note = QLabel("The figure is 180 mm wide (two journal columns); font sizes are the printed sizes.", body)
        note.setWordWrap(True)
        note.setStyleSheet("color: gray;")
        side.addWidget(note)
        side.addStretch(1)
        return body

    def _sync_paper_side(self) -> None:
        opts = self._settings.paper
        for key, box in self._paper_boxes.items():
            box.blockSignals(True)
            box.setChecked(bool(opts.get(key, True)))
            box.blockSignals(False)
        self._paper_font.blockSignals(True)
        self._paper_font.setCurrentFont(QFont(str(opts.get("font") or "Arial")))
        self._paper_font.blockSignals(False)
        self._paper_size.blockSignals(True)
        self._paper_size.setValue(float(opts.get("size") or 8.0))
        self._paper_size.blockSignals(False)

    def _refresh_paper(self) -> None:
        if self._mode != "paper" or self._folder is None:
            return
        self._sync_paper_side()
        k = self._karyotype
        self._paper.karyotype = k
        self._paper.bands = self._pooled(k) if k is not None else []
        self._paper.colours = self._edit.probe_colours
        self._paper.opts = self._settings.paper
        self._paper.update()

    def _on_paper_option(self, *_args) -> None:
        opts = self._settings.paper
        for key, box in self._paper_boxes.items():
            opts[key] = box.isChecked()
        opts["font"] = self._paper_font.currentFont().family()
        opts["size"] = float(self._paper_size.value())
        if self._folder is not None:
            idg.save_settings(self._folder, self._settings)
        self._paper.opts = opts
        self._paper.update()

    def _save_paper(self, kind: str) -> None:
        from PySide6.QtCore import QMarginsF, QSizeF
        from PySide6.QtGui import QPageLayout, QPageSize, QPdfWriter

        filters = {"svg": "SVG (*.svg)", "pdf": "PDF (*.pdf)", "tif": "TIFF (*.tif *.tiff)"}
        path, _ = QFileDialog.getSaveFileName(self, "Save idiogram", self._default_name(f"_paper.{kind}"), filters[kind])
        if not path:
            return
        w, h = self._paper.logical_size()
        try:
            if kind == "svg":
                from PySide6.QtSvg import QSvgGenerator

                gen = QSvgGenerator()
                gen.setFileName(path)
                gen.setResolution(72)
                gen.setSize(QSize(int(w), int(h)))
                gen.setViewBox(QRectF(0, 0, w, h))
                gen.setTitle("Idiogram")
                p = QPainter(gen)
                self._paper.paint_into(p)
                p.end()
            elif kind == "pdf":
                writer = QPdfWriter(path)
                writer.setResolution(1200)
                writer.setPageLayout(QPageLayout(QPageSize(QSizeF(w, h), QPageSize.Unit.Point),
                                                 QPageLayout.Orientation.Portrait, QMarginsF(0, 0, 0, 0)))
                p = QPainter(writer)
                p.scale(1200 / 72, 1200 / 72)
                self._paper.paint_into(p)
                p.end()
            else:
                dpi = 600
                image = QImage(int(w * dpi / 72), int(h * dpi / 72), QImage.Format.Format_RGB32)
                image.setDotsPerMeterX(int(dpi / 0.0254))
                image.setDotsPerMeterY(int(dpi / 0.0254))
                image.fill(QColor("#FFFFFF"))
                p = QPainter(image)
                p.scale(dpi / 72, dpi / 72)
                self._paper.paint_into(p)
                p.end()
                from PySide6.QtGui import QImageWriter

                writer = QImageWriter(path, b"tiff")
                writer.setCompression(1)  # LZW, lossless
                if not writer.write(image):
                    raise OSError(f"Could not write {path}: {writer.errorString()}")
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Save", str(exc))
            return

    def _add_to_figure(self, style: str) -> None:
        self._settings.figure_style = style
        if self._folder is not None:
            idg.save_settings(self._folder, self._settings)
        self.figureRequested.emit()
