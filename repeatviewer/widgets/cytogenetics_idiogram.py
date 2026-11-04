"""Idiogram tab of the Cytogenetics workspace: in situ (FISH) against in silico (assembly)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QPointF, QRectF, QSize, Qt, Signal
from PySide6.QtGui import QColor, QFont, QImage, QPainter, QPen
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QToolTip,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import idiogram as idg
from repeatviewer.backends.cytogenetics import karyotype as kt
from repeatviewer.backends.cytogenetics.model import image_details
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.help_button import make_help_button
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
    "Each chromosome has two bars. <b>Left (in situ)</b>: the FISH signals you marked in the Signals tab, pooled over all cells; "
    "the number is on how many of the marked homologues the signal was seen.<br>"
    "<b>Right (in silico)</b>: the assembly sequence linked to the chromosome, with the arrays of the probe "
    "(hits closer than 50 kb are merged; short arrays are hidden because FISH cannot see them).<br>"
    "<b>Status</b>: Agrees when every FISH band has an array within the tolerance and the other way round; "
    "FISH only, Assembly only or Differs otherwise; No assembly when the chromosome has no linked sequence.<br>"
    "<b>Drag a band</b> (when allowed in Layout) to correct where the signal sits; all signals pooled in that band are changed "
    "and kept as set by hand.<br>"
    "<b>Chromosome names</b>: which assembly sequence belongs to which FISH chromosome. Tick Flip when the sequence starts at the q end.<br>"
    "Signals without an arm (chromosome not measured) cannot be placed; measure the chromosome in the Signals tab."
)


def _section_label(text: str) -> QLabel:
    label = QLabel(text)
    label.setWordWrap(True)
    label.setStyleSheet("color: palette(dark); font-size: 12px;")
    return label


class IdiogramCanvas(QWidget):
    """Pairs of bars, one pair per chromosome (or homologue)."""

    bandMoved = Signal(object, float)  # FishBand, new fraction
    orderChanged = Signal(list)  # chromosome names

    SLOT_MIN = 74
    FISH_W, ASM_W, GAP = 26, 16, 8
    TOP, BOTTOM = 26, 78

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.rows: list[idg.ChromosomeRow] = []
        self.heights: list[float] = []
        self.colours: dict[str, str] = {}
        self.drag_enabled = True
        self.reorder_enabled = False
        self._drag_band: tuple[int, int] | None = None
        self._drag_frac = 0.0
        self._drag_slot: int | None = None
        self._drag_x = 0.0
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.setMinimumHeight(420)

    def set_rows(self, rows: list[idg.ChromosomeRow], heights: list[float], colours: dict[str, str]) -> None:
        self.rows, self.heights, self.colours = rows, heights, colours
        self.setMinimumWidth(max(300, len(rows) * self.SLOT_MIN + 20))
        self.update()

    def sizeHint(self) -> QSize:  # noqa: D102
        return QSize(max(600, len(self.rows) * 96), 560)

    # ---------------------------------------------------------------- geometry

    def _slot(self) -> float:
        return max(float(self.SLOT_MIN), (self.width() - 20) / max(1, len(self.rows)))

    def _bars(self, i: int) -> tuple[QRectF, QRectF]:
        slot = self._slot()
        cx = 10 + slot * (i + 0.5)
        usable = max(80.0, self.height() - self.TOP - self.BOTTOM)
        h = usable * (self.heights[i] if i < len(self.heights) else 1.0)
        bottom = self.TOP + usable
        pair = self.FISH_W + self.GAP + self.ASM_W
        x0 = cx - pair / 2 + 8  # leave room on the left for the "seen" numbers
        fish = QRectF(x0, bottom - h, self.FISH_W, h)
        asm = QRectF(x0 + self.FISH_W + self.GAP, bottom - h, self.ASM_W, h)
        return fish, asm

    def _band_rect(self, fish: QRectF, row: idg.ChromosomeRow, band: idg.FishBand, fraction: float | None = None) -> QRectF:
        f = band.fraction if fraction is None else fraction
        tall = band.arm == "cen"
        h = 20.0 if tall else 12.0
        y = fish.top() + f * fish.height()
        return QRectF(fish.left(), y - h / 2, fish.width(), h)

    def colour(self, probe: str) -> QColor:
        return QColor(self.colours.get(probe, "#8A8D93"))

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
                       "No karyotype yet. Make one in the Signals tab (Edit karyotype).")
            return
        small = QFont(self.font())
        small.setPointSizeF(max(7.5, small.pointSizeF() - 2))
        bold = QFont(self.font())
        bold.setBold(True)
        for i, row in enumerate(self.rows):
            fish, asm = self._bars(i)
            if self._drag_slot == i:
                shift = self._drag_x - (fish.left() + fish.width())
                fish.translate(shift, 0)
                asm.translate(shift, 0)
            # FISH bar
            p.setPen(QPen(QColor("#6C6F75"), 1.5))
            p.setBrush(QColor("#FFFFFF"))
            p.drawRoundedRect(fish, self.FISH_W / 2, self.FISH_W / 2)
            for j, band in enumerate(row.fish):
                frac = self._drag_frac if self._drag_band == (i, j) else None
                r = self._band_rect(fish, row, band, frac)
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(self.colour(band.probe))
                p.drawRect(r)
                if j not in row.matched_fish and row.sequence_length:
                    p.setPen(QPen(QColor("#A3213B"), 1.2, Qt.PenStyle.DashLine))
                    p.setBrush(Qt.BrushStyle.NoBrush)
                    p.drawRect(r.adjusted(-3, -2, 3, 2))
                p.setFont(small)
                p.setPen(QColor("#45484D"))
                p.drawText(QRectF(fish.left() - 46, r.center().y() - 8, 42, 16),
                           Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter, band.label)
            cy = fish.top() + row.chromosome.centromere_fraction * fish.height()
            p.setPen(QPen(QColor("#6C6F75"), 1.5))
            p.setBrush(QColor("#FFFFFF"))
            p.drawRoundedRect(QRectF(fish.left() - 2, cy - 3, fish.width() + 4, 6), 3, 3)
            # assembly bar
            if row.sequence_length:
                p.setPen(QPen(QColor("#8A8D93"), 1, Qt.PenStyle.DashLine))
                p.setBrush(QColor("#EEF0F4"))
                p.drawRoundedRect(asm, 4, 4)
                p.setPen(Qt.PenStyle.NoPen)
                for array in row.silico:
                    top = asm.top() + array.start * asm.height()
                    h = max(3.0, (array.end - array.start) * asm.height())
                    p.setBrush(self.colour(array.probe))
                    p.drawRect(QRectF(asm.left() + 1, top, asm.width() - 2, h))
            # labels
            bottom = max(fish.bottom(), asm.bottom())
            slot = self._slot()
            text_rect = QRectF(fish.center().x() + (self.ASM_W + self.GAP) / 2 - slot / 2, bottom + 6, slot, 18)
            p.setFont(bold)
            p.setPen(QColor("#1C1C1E"))
            p.drawText(text_rect, Qt.AlignmentFlag.AlignCenter, row.name)
            p.setFont(small)
            p.setPen(QColor("#5A5D63"))
            seq = row.sequence or "not linked"
            p.drawText(text_rect.translated(0, 18), Qt.AlignmentFlag.AlignCenter,
                       p.fontMetrics().elidedText(seq, Qt.TextElideMode.ElideMiddle, int(slot) - 4))
            bg, fg = STATUS_COLOURS.get(row.status, ("#F1F1F3", "#45484D"))
            w = p.fontMetrics().horizontalAdvance(row.status) + 16
            pill = QRectF(text_rect.center().x() - w / 2, text_rect.top() + 38, w, 18)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(bg))
            p.drawRoundedRect(pill, 9, 9)
            p.setPen(QColor(fg))
            p.drawText(pill, Qt.AlignmentFlag.AlignCenter, row.status)

    # ---------------------------------------------------------------- mouse

    def _hit(self, pos: QPointF):
        for i, row in enumerate(self.rows):
            fish, asm = self._bars(i)
            for j, band in enumerate(row.fish):
                if self._band_rect(fish, row, band).adjusted(-4, -3, 4, 3).contains(pos):
                    return ("band", i, j)
            if row.sequence_length:
                for j, array in enumerate(row.silico):
                    top = asm.top() + array.start * asm.height()
                    h = max(5.0, (array.end - array.start) * asm.height())
                    if QRectF(asm.left() - 3, top - 2, asm.width() + 6, h + 4).contains(pos):
                        return ("array", i, j)
            name_rect = QRectF(fish.left() - 20, max(fish.bottom(), asm.bottom()), 80, 60)
            if name_rect.contains(pos) or fish.contains(pos) or asm.contains(pos):
                return ("slot", i, -1)
        return None

    def mousePressEvent(self, event) -> None:  # noqa: N802
        hit = self._hit(event.position())
        if hit is None or event.button() != Qt.MouseButton.LeftButton:
            return
        kind, i, j = hit
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
            if j not in row.matched_fish and row.sequence_length:
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
            dst = int(max(0, min(len(self.rows) - 1, (event.position().x() - 10) // self._slot())))
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


class IdiogramTab(QWidget):
    name = "Idiogram"
    karyotypesChanged = Signal()

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
        self._legend = QHBoxLayout()
        self._legend.setSpacing(14)
        head = QHBoxLayout()
        head.addWidget(self._title)
        head.addStretch(1)
        head.addLayout(self._legend)
        head.addWidget(make_help_button(self, "Idiogram", HELP_TEXT))

        self._canvas = IdiogramCanvas(self)
        self._canvas.bandMoved.connect(self._on_band_moved)
        self._canvas.orderChanged.connect(self._on_order_dragged)
        frame = QFrame(self)
        frame.setStyleSheet("QFrame#idioFrame { border:1px solid #E6E6EA; border-radius:8px; background:#FFFFFF; }")
        frame.setObjectName("idioFrame")
        fl = QVBoxLayout(frame)
        fl.setContentsMargins(8, 12, 8, 8)
        scroll = QScrollArea(frame)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(self._canvas)
        scroll.setStyleSheet("QScrollArea { background:#FFFFFF; }")
        fl.addWidget(scroll)

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
        main.setContentsMargins(24, 16, 24, 16)
        main.setSpacing(14)
        main.addLayout(head)
        main.addWidget(frame, 1)
        main.addLayout(notes)
        main.addWidget(self._note_warn)
        main_host = QWidget(self)
        main_host.setObjectName("cytoMain")
        main_host.setLayout(main)

        # ---- side panel
        body = QWidget(self)
        body.setObjectName("cytoAside")
        side = QVBoxLayout(body)
        side.setContentsMargins(14, 12, 14, 12)
        side.setSpacing(14)

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
        lab = QLabel("Show signals seen in at least", self)
        lab.setWordWrap(True)
        lab.setFixedWidth(110)
        r.addWidget(lab)
        r.addWidget(self._share, 1)
        self._probe_section.body_layout().addLayout(r)
        side.addWidget(self._probe_section)

        self._silico_section = _CollapsibleSection("In silico", body, expanded=True)
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

        self._layout_section = _CollapsibleSection("Layout", body, expanded=True)
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
        for text, w in (("Karyotype", self._karyo), ("Order", self._order), ("Lengths", self._lengths_box)):
            r = QHBoxLayout()
            lab = QLabel(text, self)
            lab.setFixedWidth(90)
            r.addWidget(lab)
            r.addWidget(w, 1)
            self._layout_section.body_layout().addLayout(r)
        self._homs = QCheckBox("Draw each homologue (a, b, c) separately", self)
        self._homs.toggled.connect(self._on_layout)
        self._drag = QCheckBox("Drag signals to correct their position", self)
        self._drag.toggled.connect(self._on_layout)
        self._layout_section.add_widget(self._homs)
        self._layout_section.add_widget(self._drag)
        edit = QPushButton("Edit karyotype (lengths, centromeres, bands)", self)
        edit.setFlat(True)
        edit.setStyleSheet(LINK_CSS)
        edit.setCursor(Qt.CursorShape.PointingHandCursor)
        edit.clicked.connect(self._edit_karyotype)
        self._layout_section.add_widget(edit)
        side.addWidget(self._layout_section)

        self._names_section = _CollapsibleSection("Chromosome names", body, expanded=True)
        self._names_section.add_widget(_section_label("FISH number to assembly sequence (edit when the order differs)"))
        self._names_grid = QGridLayout()
        self._names_grid.setHorizontalSpacing(8)
        self._names_grid.setVerticalSpacing(4)
        self._names_section.body_layout().addLayout(self._names_grid)
        side.addWidget(self._names_section)
        side.addStretch(1)

        save = QPushButton("Save ▾", self)
        menu = QMenu(save)
        menu.addAction("Picture (PNG)…", self._save_png)
        menu.addAction("Vector (SVG)…", self._save_svg)
        menu.addAction("Comparison table (TSV)…", self._save_tsv)
        save.setMenu(menu)
        fig = QPushButton("Add to figure", self)
        fig.setEnabled(False)
        fig.setToolTip("The Figure tab comes in the next step")
        bottom = QHBoxLayout()
        bottom.addWidget(save)
        bottom.addWidget(fig)
        bottom.addStretch(1)
        side.addLayout(bottom)

        aside = QScrollArea(self)
        aside.setWidgetResizable(True)
        aside.setFrameShape(QFrame.Shape.NoFrame)
        aside.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        aside.setWidget(body)
        aside.setFixedWidth(SIDE_W)

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addWidget(main_host, 1)
        line = QFrame(self)
        line.setFrameShape(QFrame.Shape.VLine)
        line.setStyleSheet("color: palette(mid);")
        root.addWidget(line)
        root.addWidget(aside)
        self.setStyleSheet("#cytoAside { background: palette(alternate-base); } #cytoMain { background: palette(base); }")

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
        self._canvas.drag_enabled = s.drag
        self._canvas.reorder_enabled = s.order == "custom"
        self._canvas.set_rows(self._rows, idg.relative_heights(self._rows, s.lengths), colours)
        head, species = idg.title_text(k, probes)
        self._title.setText(f"{head} <i>{species}</i> karyotype")
        self._fill_legend(probes, colours)
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
        while layout.count():
            item = layout.takeAt(0)
            if item.widget() is not None:
                item.widget().deleteLater()
            elif item.layout() is not None:
                self._clear(item.layout())

    def _fill_legend(self, probes: list[str], colours: dict[str, str]) -> None:
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
        for probe in self._all_probes_with(fish_cells):
            n = fish_cells.get(probe, 0)
            box = QCheckBox(f"{probe} ({n} cell{'s' if n != 1 else ''})", self._probe_host)
            box.setChecked(probe not in s.hidden_probes)
            box.toggled.connect(lambda on, p=probe: self._on_probe(p, on))
            self._probe_lay.addWidget(box)
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
        for box, value in ((self._homs, s.homologues), (self._drag, s.drag)):
            box.blockSignals(True)
            box.setChecked(value)
            box.blockSignals(False)

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
        self._settings.karyotype = str(self._karyo.itemData(index) or "")
        self._karyotype = next((k for k in self._karyotypes if k.name == self._settings.karyotype), self._karyotype)
        self._store(controls=True)

    def _on_layout(self, *_args) -> None:
        s = self._settings
        s.order = str(self._order.currentData() or "karyotype")
        s.lengths = str(self._lengths_box.currentData() or "karyotype")
        s.homologues = self._homs.isChecked()
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

        dialog = KaryotypeDialog(self._folder, self._karyotype.name if self._karyotype else "", self)
        if dialog.exec():
            self._settings.karyotype = dialog.saved_name() or self._settings.karyotype
            idg.save_settings(self._folder, self._settings)
            self.karyotypesChanged.emit()
            self._rebuild()

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
