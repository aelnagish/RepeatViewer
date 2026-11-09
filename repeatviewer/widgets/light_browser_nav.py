"""Navigation bar of the light genome browser: back / forward, the genome chip, the
sequence chip (a filterable list in genome order) and the position box.

Vendored from GenomeFlow ``genomeflow/widgets/genome_navigator.py`` (commit 5c2221b:
ChromPicker, BreadcrumbBar and their chips) so the Preview looks like the genome
browser. Refresh by copying again rather than editing here.
"""

from __future__ import annotations

from PySide6.QtCore import QPointF, QRectF, Qt, Signal
from PySide6.QtGui import QColor, QFont, QFontMetrics, QPainter, QPen
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QStyle,
    QStyleOptionToolButton,
    QStylePainter,
    QToolButton,
    QVBoxLayout,
    QWidget,
)


class ChromPicker(QFrame):
    """The sequence list under the chromosome chip: a filter box and one scrolling
    column, never taller than about 14 rows."""

    picked = Signal(str)

    def __init__(self, chroms: list[tuple[str, int]], current: str,
                 parent: QWidget | None = None) -> None:
        from PySide6.QtWidgets import QListWidget, QListWidgetItem

        super().__init__(parent, Qt.WindowType.Popup | Qt.WindowType.FramelessWindowHint)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
        # a rounded card like the menus: only the card is painted, the window
        # around its corners stays transparent
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground)
        self.setFrameShape(QFrame.Shape.NoFrame)
        pal = self.palette()
        base = pal.color(pal.ColorRole.Base)
        dark = base.lightness() < 128
        self._card = (QColor(base), QColor("#3a3f47" if dark else "#d5d9e0"))
        self.setStyleSheet(
            "QListWidget { border: none; background: transparent; }"
            "QListWidget::item { padding: 2px 6px; border-radius: 5px; }"
            "QListWidget::item:selected { background: #3d7be0; color: white; }"
            "QListWidget::item:hover:!selected { background: rgba(61, 123, 224, 40); }")
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 8, 8, 8)
        lay.setSpacing(4)
        self.filter = QLineEdit(self)
        self.filter.setPlaceholderText(f"Filter {len(chroms):,} sequences…")
        self.filter.setClearButtonEnabled(True)
        lay.addWidget(self.filter)
        self.list = QListWidget(self)
        self.list.setUniformItemSizes(True)
        mono = QFont(self.font())
        for name, n in chroms:
            it = QListWidgetItem(f"{name}   {n:,} bp")
            it.setData(Qt.ItemDataRole.UserRole, name)
            it.setFont(mono)
            self.list.addItem(it)
            if name == current:
                self.list.setCurrentItem(it)
        lay.addWidget(self.list)
        rows = min(14, max(3, len(chroms)))
        self.list.setFixedHeight(rows * (self.list.sizeHintForRow(0) or 20) + 6)
        self.setFixedWidth(max(240, self.list.sizeHintForColumn(0) + 40))
        self.filter.textChanged.connect(self._filter)
        self.filter.returnPressed.connect(self._pick_first)
        self.list.itemClicked.connect(self._pick)
        self.list.itemActivated.connect(self._pick)
        if self.list.currentItem() is not None:
            self.list.scrollToItem(self.list.currentItem(),
                                   self.list.ScrollHint.PositionAtCenter)
        self.filter.setFocus()

    def paintEvent(self, event) -> None:
        fill, edge = self._card
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        p.setPen(QPen(edge, 1))
        p.setBrush(fill)
        p.drawRoundedRect(QRectF(self.rect()).adjusted(0.5, 0.5, -0.5, -0.5), 10, 10)
        p.end()

    def _filter(self, text: str) -> None:
        q = text.strip().lower()
        for i in range(self.list.count()):
            it = self.list.item(i)
            it.setHidden(bool(q) and q not in str(it.data(Qt.ItemDataRole.UserRole)).lower())

    def _pick_first(self) -> None:
        for i in range(self.list.count()):
            it = self.list.item(i)
            if not it.isHidden():
                self._pick(it)
                return

    def _pick(self, item) -> None:
        self.picked.emit(str(item.data(Qt.ItemDataRole.UserRole)))
        self.close()


# ---------------------------------------------------------------- breadcrumb bar


_CHIP_CSS = """
QToolButton#navChip {{ border: 1px solid {line}; border-radius: 7px; padding: 3px 10px;
                       background: palette(base); }}
QToolButton#navChip:hover {{ border-color: #9aa3ad; }}
QToolButton#navChip[current="true"] {{ border: 1px solid #3d7be0; color: #2a62c9;
                                        background: {sel}; font-weight: 600; }}
QToolButton#navRound {{ border: 1px solid {line}; border-radius: 13px; min-width: 24px;
                        min-height: 24px; max-width: 24px; max-height: 24px;
                        background: palette(base); }}
QToolButton#navRound:disabled {{ color: #b5bcc4; }}
QLabel#navSep {{ color: #9aa3ad; }}
QLineEdit#navPos {{ border: 1px solid {line}; border-radius: 7px; padding: 0px 8px;
                   background: palette(base); }}
QLineEdit#navPos:focus {{ border-color: #3d7be0; }}
"""


class _LeftChip(QToolButton):
    """A chip whose text starts at the left (a tool button always centres it), at the
    same inset as the text in the position box beside it."""

    TEXT_INSET = 9          # the position box: 1 px border + 8 px padding

    def paintEvent(self, event) -> None:
        p = QStylePainter(self)
        opt = QStyleOptionToolButton()
        self.initStyleOption(opt)
        text = opt.text
        opt.text = ""
        p.drawComplexControl(QStyle.ComplexControl.CC_ToolButton, opt)   # the frame
        p.setPen(self.palette().color(self.foregroundRole()))
        p.drawText(self.rect().adjusted(self.TEXT_INSET, 0, -self.TEXT_INSET, 0),
                   Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft, text)
        p.end()


class _ArrowButton(QToolButton):
    """A round Back / Forward button whose chevron is drawn at the exact centre
    (a font's ‹ › glyphs sit off-centre)."""

    def __init__(self, direction: int, parent=None) -> None:
        super().__init__(parent)
        self.direction = direction

    def paintEvent(self, event) -> None:
        super().paintEvent(event)       # the circle, from the style sheet
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        col = self.palette().color(self.foregroundRole())
        if not self.isEnabled():
            col = QColor("#b5bcc4")
        pen = QPen(col, 1.8)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        pen.setJoinStyle(Qt.PenJoinStyle.RoundJoin)
        p.setPen(pen)
        cx, cy = self.width() / 2, self.height() / 2
        w, h = 2.6, 5.0                 # half width and half height of the chevron
        tip = cx + self.direction * w
        back = cx - self.direction * w
        p.drawPolyline([QPointF(back, cy - h), QPointF(tip, cy), QPointF(back, cy + h)])
        p.end()


class BreadcrumbBar(QWidget):
    """‹ › chr4 › 2.10–2.18 Mb › BLGH_03307: back / forward and where the view is.
    The chromosome chip lists the sequences, the range chip edits the position, the
    gene chip zooms to the gene in the middle of the view."""

    backRequested = Signal()
    forwardRequested = Signal()
    chromRequested = Signal(str)
    rangeClicked = Signal()
    positionEntered = Signal(str)       # the position box: Enter on chrom:start-end
    geneClicked = Signal()

    GENOME_W = 170
    CHROM_W = 110
    POS_W = 190                         # until the font is known (then POS_SAMPLE fits)
    POS_SAMPLE = "1,000,000,000-1,000,000,000"      # a region on a 1 Gbp chromosome

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(6)
        self.back = self._round("‹", "Back (⌘[)", self.backRequested)
        self.forward = self._round("›", "Forward (⌘])", self.forwardRequested)
        lay.addWidget(self.back)
        lay.addWidget(self.forward)
        lay.addSpacing(4)
        # the genome (a drop-down the browser window fills), then the sequence: both
        # at a fixed width, long names shortened in the middle with …
        self.genome = self._chip("Genome: open another one, or one used before", left=True)
        self.genome.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self.genome.setFixedWidth(self.GENOME_W)
        self.genome.setStyleSheet("QToolButton::menu-indicator { image: none; width: 0px; }")
        self.genome.setVisible(False)
        self.genome_sep = self._sep()
        self.genome_sep.setVisible(False)       # no › between the parts
        lay.addWidget(self.genome)
        lay.addWidget(self.genome_sep)
        self.chrom = self._chip("Sequence: pick another one", left=True)
        self.chrom.setFixedWidth(self.CHROM_W)
        self._chrom_name = ""
        self._genome_name = ""
        self.chrom.clicked.connect(self._open_picker)
        self.picker: ChromPicker | None = None
        # the region in view as chrom:start-end, typed over to go elsewhere; a fixed
        # width, and a long one keeps its end in sight (the start scrolls out)
        self.range = QLineEdit(self)
        self.range.setObjectName("navPos")
        self.range.setFixedWidth(self.POS_W)
        self.range.setToolTip("The region in view on this sequence. Type a position, a "
                              "region (1,710,785-1,710,959, or scaffold_9:… for another "
                              "sequence) or a gene and press Enter.")
        self.range.returnPressed.connect(
            lambda: self.positionEntered.emit(self.range.text().strip()))
        self.range.editingFinished.connect(self._show_where)
        # the same font and height as the genome and sequence chips beside it
        self.range.setFont(self.chrom.font())
        self._where = ""
        self.gene = self._chip("The feature in the middle of the view: click to zoom to it")
        self.gene.clicked.connect(self.geneClicked)
        self.gene.setProperty("current", True)
        self.seps = [self._sep(), self._sep()]
        lay.addWidget(self.chrom)
        lay.addWidget(self.seps[0])
        lay.addWidget(self.range)
        lay.addWidget(self.seps[1])
        lay.addWidget(self.gene)
        self._chroms: list[tuple[str, int]] = []
        line = "#3a3f47" if self.palette().color(
            self.palette().ColorRole.Base).lightness() < 128 else "#d5d9e0"
        self.setStyleSheet(_CHIP_CSS.format(line=line, sel="rgba(61,123,224,0.10)"))
        self.set_where("", 0, 0, "")
        self.set_history(False, False)

    def _round(self, text: str, tip: str, signal) -> QToolButton:
        b = _ArrowButton(-1 if text == "‹" else 1, self)
        b.setObjectName("navRound")
        b.setToolTip(tip)
        b.clicked.connect(signal)
        return b

    def _chip(self, tip: str, left: bool = False) -> QToolButton:
        b = _LeftChip(self) if left else QToolButton(self)
        b.setObjectName("navChip")
        b.setToolTip(tip)
        b.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        return b

    def _sep(self) -> QLabel:
        s = QLabel("›", self)
        s.setObjectName("navSep")
        return s

    def set_chroms(self, chroms: list[tuple[str, int]]) -> None:
        self._chroms = list(chroms)

    def _fit(self, button: QToolButton, text: str, width: int, suffix: str = "") -> None:
        """``text`` (+ suffix) shortened in the middle to fit a chip ``width`` wide."""
        room = width - 2 * _LeftChip.TEXT_INSET - 4 - \
            button.fontMetrics().horizontalAdvance(suffix)
        short = button.fontMetrics().elidedText(text, Qt.TextElideMode.ElideMiddle, room)
        button.setText(short + suffix)

    def set_genome(self, name: str, menu: QMenu | None = None) -> None:
        """Show the genome chip (``name``, a drop-down ``menu``) before the sequence."""
        self._genome_name = name
        if menu is not None:
            self.genome.setMenu(menu)
        self._fit(self.genome, name, self.GENOME_W)
        self.genome.setToolTip(f"Genome: {name}\nClick to open another one, or one used "
                               "before")
        self.genome.setVisible(bool(name))

    def _open_picker(self) -> None:
        self.picker = ChromPicker(self._chroms, self._chrom_name, self)
        self.picker.picked.connect(self.chromRequested)
        self.picker.move(self.chrom.mapToGlobal(self.chrom.rect().bottomLeft()))
        self.picker.show()

    def _match_chip_font(self) -> None:
        """The position box uses the very font of the genome / sequence chips (on macOS a
        line edit's own font is larger than a tool button's)."""
        self.chrom.ensurePolished()
        f = QFont(self.chrom.font())
        f.setPointSizeF(f.pointSizeF())     # explicit, so the class font cannot win
        for w in (self.genome, self.chrom, self.range):
            if w.font() != f or not w.testAttribute(Qt.WidgetAttribute.WA_SetFont):
                w.setFont(f)
        # just wide enough for the longest likely region (1 Gbp chromosome), plus the
        # box's padding (8 px each side) and border; a longer one keeps its end in view
        width = QFontMetrics(f).horizontalAdvance(self.POS_SAMPLE) + 2 * 8 + 2 + 6
        if self.range.width() != width:
            self.range.setFixedWidth(width)

    def showEvent(self, event) -> None:
        super().showEvent(event)
        self._match_chip_font()

    def set_where(self, chrom: str, a: int, b: int, gene: str) -> None:
        self._chrom_name = chrom
        self._fit(self.chrom, chrom or "—", self.CHROM_W)
        self.chrom.setToolTip(f"{chrom}\nClick to pick another sequence" if chrom
                              else "Sequence: pick another one")
        self.chrom.setVisible(bool(chrom))
        has_range = bool(chrom) and b > a
        self._where = f"{a + 1:,}-{b:,}" if has_range else ""
        self._match_chip_font()
        self.range.setFixedHeight(self.chrom.sizeHint().height())   # as tall as the chips
        self._show_where()
        self.range.setVisible(has_range)
        for sep in self.seps:
            sep.setVisible(False)               # no › between the parts
        self.gene.setText(gene)
        self.gene.setVisible(False)             # no gene name after the position

    def _show_where(self) -> None:
        """The region in view, unless it is being typed over; the end stays in sight."""
        if self.range.hasFocus():
            return
        self.range.setText(self._where)
        self.range.setCursorPosition(len(self._where))

    def set_history(self, can_back: bool, can_forward: bool) -> None:
        self.back.setEnabled(can_back)
        self.forward.setEnabled(can_forward)
