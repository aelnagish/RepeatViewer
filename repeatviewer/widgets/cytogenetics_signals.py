"""Signals tab of the Cytogenetics workspace: number chromosomes, mark signals, measure axes."""

from __future__ import annotations

import math
from pathlib import Path

from PySide6.QtCore import QPointF, QRectF, QSettings, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QPainter, QPen
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import signals as sg
from repeatviewer.backends.cytogenetics import signals_session as ss
from repeatviewer.backends.cytogenetics.karyotype import list_karyotypes
from repeatviewer.backends.cytogenetics.model import load_slides
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cytogenetics_style import clear_layout
from repeatviewer.widgets.cytogenetics_image import ImageCanvas, rgb_to_qimage
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.sequence_annotation_rail import _CollapsibleSection

_SETTINGS = "cytogenetics_workspace"
SIDE_W = 340
SAVE_DELAY_MS = 800
GOLD = QColor("#FFD35C")
LABEL_COLOUR = QColor("#FFE08A")
NEW_KARYOTYPE = "New karyotype…"

HELP_TEXT = (
    "<b>Signals</b><br>"
    "Name the chromosomes of one cell, then mark where each probe signal sits.<br>"
    "<b>1 Number chromosomes</b>: click a chromosome to give it the next name (1a, 1b, 1c, 2a …). "
    "Its axis is found from the DAPI image: square handles are the telomeres (p, q), the round one the centromere. "
    "Drag them when they are wrong. Click a name on the right to choose another one.<br>"
    "<b>2 Mark signals</b>: click a spot to add a signal of the chosen probe. With <b>Snap</b> on it moves to the centre "
    "of the spot. It is given to the nearest chromosome, and its arm and position are worked out from the axis.<br>"
    "<b>In both steps</b>: drag a name or a signal to move it, right-click for more (change probe, move to another "
    "chromosome, strength, delete), Backspace deletes the selected one.<br>"
    "<b>Moving the picture</b>: mouse wheel or trackpad zooms, Option + drag or the middle button moves it.<br>"
    "Arm and position that you set by hand are never recomputed."
)


def _read(key: str, default: str) -> str:
    try:
        return str(QSettings().value(f"{_SETTINGS}/{key}", default))
    except Exception:  # noqa: BLE001
        return default


def _write(key: str, value) -> None:
    try:
        QSettings().setValue(f"{_SETTINGS}/{key}", value)
    except Exception:  # noqa: BLE001
        pass


class MarkCanvas(ImageCanvas):
    """Image with chromosome labels, axes and signal rings. Emits clicks in image pixels."""

    pressedAt = Signal(float, float)
    movedTo = Signal(float, float)
    releasedAt = Signal(float, float)
    deleteRequested = Signal()
    menuAt = Signal(float, float, object)  # image x, y, global position

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.session: ss.SignalsSession | None = None
        self.selected: tuple[str, str] | None = None
        self.rubber: tuple[QPointF, QPointF] | None = None
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self._tool_cursor = Qt.CursorShape.CrossCursor
        self._panning = False

    def set_tool_cursor(self, shape: Qt.CursorShape) -> None:
        self._tool_cursor = shape
        self.setCursor(shape)

    def mousePressEvent(self, event) -> None:  # noqa: N802
        self.setFocus()
        alt = bool(event.modifiers() & Qt.KeyboardModifier.AltModifier)
        if event.button() == Qt.MouseButton.MiddleButton or (event.button() == Qt.MouseButton.LeftButton and alt):
            self._panning = True
            self.start_pan(event.position())
        elif event.button() == Qt.MouseButton.LeftButton:
            p = self.image_point(event.position())
            self.pressedAt.emit(p.x(), p.y())
        elif event.button() == Qt.MouseButton.RightButton:
            p = self.image_point(event.position())
            self.menuAt.emit(p.x(), p.y(), event.globalPosition().toPoint())

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if self._panning:
            self.move_pan(event.position())
        elif event.buttons() & Qt.MouseButton.LeftButton:
            p = self.image_point(event.position())
            self.movedTo.emit(p.x(), p.y())

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if self._panning:
            self._panning = False
            self.end_pan()
            self.setCursor(self._tool_cursor)
        elif event.button() == Qt.MouseButton.LeftButton:
            p = self.image_point(event.position())
            self.releasedAt.emit(p.x(), p.y())

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        pass

    def keyPressEvent(self, event) -> None:  # noqa: N802
        if event.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            self.deleteRequested.emit()
        else:
            super().keyPressEvent(event)

    def paintEvent(self, event) -> None:  # noqa: N802
        super().paintEvent(event)
        s = self.session
        if s is None:
            return
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        font = QFont(self.font())
        font.setBold(True)
        font.setPointSizeF(max(10.0, font.pointSizeF() + 1))
        p.setFont(font)
        for m in s.details.chromosomes:
            chosen = self.selected == ("mark", m.id)
            pts = [q for q in (m.telomere_p, m.centromere, m.telomere_q) if q is not None]
            if len(pts) >= 2:
                p.setPen(QPen(QColor(255, 255, 255, 150), 1, Qt.PenStyle.DashLine))
                for a, b in zip(pts, pts[1:]):
                    p.drawLine(self.screen_point(*a), self.screen_point(*b))
                if m.centromere is not None:
                    p.setBrush(QColor(255, 255, 255, 200))
                    p.setPen(Qt.PenStyle.NoPen)
                    p.drawEllipse(self.screen_point(*m.centromere), 3.0, 3.0)
                    p.setBrush(Qt.BrushStyle.NoBrush)
            anchor = self.screen_point(m.x, m.y)
            rect = QRectF(anchor.x() - 30, anchor.y() - 26, 60, 18)
            p.setPen(QColor(0, 0, 0))
            p.drawText(rect.translated(1, 1), Qt.AlignmentFlag.AlignCenter, m.label)
            p.setPen(LABEL_COLOUR)
            p.drawText(rect, Qt.AlignmentFlag.AlignCenter, m.label)
            if chosen:
                p.setPen(QPen(GOLD, 2, Qt.PenStyle.DashLine))
                p.drawEllipse(anchor, 16.0, 16.0)
                # handles: square telomeres, round centromere (drag them to correct)
                p.setPen(QPen(QColor("#000000"), 1))
                for which in ("telomere_p", "telomere_q"):
                    pt = getattr(m, which)
                    if pt is not None:
                        q = self.screen_point(*pt)
                        p.setBrush(QColor("#FFFFFF"))
                        p.drawRect(QRectF(q.x() - 5, q.y() - 5, 10, 10))
                        p.setPen(QColor("#FFFFFF"))
                        p.drawText(QRectF(q.x() + 7, q.y() - 9, 20, 18), Qt.AlignmentFlag.AlignVCenter, which[-1])
                        p.setPen(QPen(QColor("#000000"), 1))
                if m.centromere is not None:
                    p.setBrush(QColor("#FFFFFF"))
                    p.drawEllipse(self.screen_point(*m.centromere), 6.0, 6.0)
                p.setBrush(Qt.BrushStyle.NoBrush)
        for sig in s.details.signals:
            chosen = self.selected == ("signal", sig.id)
            c = self.screen_point(sig.x, sig.y)
            p.setPen(QPen(GOLD if chosen else QColor("#FFFFFF"), 2.5 if chosen else 1.5))
            p.drawEllipse(c, 9.0, 9.0)
        if self.rubber is not None:
            a, b = self.screen_point(self.rubber[0].x(), self.rubber[0].y()), self.screen_point(self.rubber[1].x(), self.rubber[1].y())
            p.setPen(QPen(QColor("#5CE1FF"), 2))
            p.drawLine(a, b)
            p.setBrush(QColor("#5CE1FF"))
            p.drawEllipse(a, 4.0, 4.0)
            p.setBrush(Qt.BrushStyle.NoBrush)
            if math.hypot(b.x() - a.x(), b.y() - a.y()) > 12:
                p.drawEllipse(b, 4.0, 4.0)


def _qimage_rgb(image):
    """QImage -> numpy RGB array (for snapping signals and finding chromosome axes on plain pictures)."""
    import numpy as np
    from PySide6.QtGui import QImage

    img = image.convertToFormat(QImage.Format.Format_RGB888)
    w, h, bpl = img.width(), img.height(), img.bytesPerLine()
    buf = np.frombuffer(img.constBits(), dtype=np.uint8, count=bpl * h).reshape(h, bpl)
    return buf[:, : w * 3].reshape(h, w, 3).copy()


def _combo(items: list[tuple[str, str]], current: str) -> NoWheelComboBox:
    box = NoWheelComboBox()
    for text, data in items:
        box.addItem(text, data)
    idx = box.findData(current)
    box.setCurrentIndex(idx if idx >= 0 else 0)
    return box


class SignalsTab(QWidget):
    name = "Signals"
    imageChanged = Signal(str)
    idiogramRequested = Signal()
    karyotypesChanged = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._folder: Path | None = None
        self._session: ss.SignalsSession | None = None
        self._tool = "number"
        self._pending_label = ""
        self._drag: tuple[str, str] | None = None
        self._moved = False
        self._saving = False
        self._pending = ""
        self._label_buttons: dict[str, QPushButton] = {}

        # ---- toolbar: 1 Number chromosomes | 2 Mark signals, probe, snap, zoom, cell navigation
        self._canvas = MarkCanvas(self)
        self._canvas.pressedAt.connect(self._on_press)
        self._canvas.movedTo.connect(self._on_move)
        self._canvas.releasedAt.connect(self._on_release)
        self._canvas.deleteRequested.connect(self._delete_selected)
        self._canvas.menuAt.connect(self._on_menu)
        self._tool_group = QButtonGroup(self)
        self._tool_group.setExclusive(True)
        bar = QHBoxLayout()
        bar.setContentsMargins(10, 8, 10, 6)
        bar.setSpacing(8)
        for key, text in ss.TOOLS:
            b = QPushButton(text, self)
            b.setCheckable(True)
            b.setChecked(key == self._tool)
            b.clicked.connect(lambda _c=False, k=key: self._set_tool(k))
            self._tool_group.addButton(b)
            bar.addWidget(b)
        bar.addSpacing(12)
        bar.addWidget(QLabel("Probe", self))
        self._probe = NoWheelComboBox(self)
        self._probe.setMinimumWidth(140)
        self._probe.currentIndexChanged.connect(lambda _i: self._update_hint())
        bar.addWidget(self._probe)
        self._auto = QCheckBox("Snap to the centre of the spot", self)
        self._auto.setChecked(_read("signals_auto", "true") != "false")
        self._auto.toggled.connect(lambda on: _write("signals_auto", "true" if on else "false"))
        bar.addWidget(self._auto)
        bar.addStretch(1)
        for text, slot in (("−", lambda: self._canvas.set_zoom(self._canvas.zoom() / 1.25)),
                           ("+", lambda: self._canvas.set_zoom(self._canvas.zoom() * 1.25)),
                           ("Fit", lambda: self._canvas.fit())):
            zb = QPushButton(text, self)
            zb.clicked.connect(slot)
            bar.addWidget(zb)
        bar.addSpacing(10)
        self._prev = QPushButton("‹", self)
        self._next = QPushButton("›", self)
        self._prev.setToolTip("Previous cell of this slide")
        self._next.setToolTip("Next cell of this slide")
        self._prev.clicked.connect(lambda: self._step_cell(-1))
        self._next.clicked.connect(lambda: self._step_cell(1))
        self._nav = QLabel(self)
        bar.addWidget(self._prev)
        bar.addWidget(self._nav)
        bar.addWidget(self._next)
        bar.addWidget(make_help_button(self, "Signals", HELP_TEXT))

        # ---- bottom row: what the next click does, then the cells of the slide
        self._hint = QLabel(self)
        self._hint.setObjectName("cytoMuted")
        from PySide6.QtWidgets import QSizePolicy

        self._hint.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self._cells_row = QHBoxLayout()
        self._cells_row.setContentsMargins(0, 0, 0, 0)
        self._cells_host = QWidget(self)
        self._cells_host.setLayout(self._cells_row)
        bottom = QHBoxLayout()
        bottom.setContentsMargins(10, 6, 10, 6)
        bottom.addWidget(self._hint, 1)
        bottom.addWidget(self._cells_host)

        left = QVBoxLayout()
        left.setContentsMargins(0, 0, 0, 0)
        left.setSpacing(0)
        left.addLayout(bar)
        left.addWidget(self._canvas, 1)
        left.addLayout(bottom)

        # ---- side panel
        body = QWidget(self)
        body.setObjectName("cytoAside")
        side = QVBoxLayout(body)
        side.setContentsMargins(10, 4, 10, 8)
        side.setSpacing(12)

        self._chrom_section = _CollapsibleSection("Chromosomes", body, expanded=True)
        self._chrom_count = QLabel(self)
        self._chrom_count.setObjectName("cytoMuted")
        self._karyo = NoWheelComboBox(self)
        self._karyo.activated.connect(self._on_karyotype)
        self._edit_karyo = QPushButton("Edit karyotype", self)
        self._edit_karyo.setFlat(True)
        self._edit_karyo.setObjectName("cytoLink")
        self._edit_karyo.clicked.connect(lambda: self._open_karyotype(edit=True))
        kr = QHBoxLayout()
        kr.addWidget(QLabel("Karyotype", self))
        kr.addWidget(self._karyo, 1)
        self._labels_grid = QGridLayout()
        self._labels_grid.setSpacing(4)
        self._chrom_section.body_layout().addWidget(self._chrom_count)
        self._chrom_section.body_layout().addLayout(kr)
        self._chrom_section.body_layout().addWidget(self._edit_karyo, 0, Qt.AlignmentFlag.AlignLeft)
        self._chrom_section.body_layout().addLayout(self._labels_grid)
        side.addWidget(self._chrom_section)

        self._sel_section = _CollapsibleSection("Selected", body, expanded=True)
        self._sel_host = QWidget(self._sel_section)
        self._sel_lay = QVBoxLayout(self._sel_host)
        self._sel_lay.setContentsMargins(0, 0, 0, 0)
        self._sel_lay.setSpacing(6)
        self._sel_section.add_widget(self._sel_host)
        side.addWidget(self._sel_section)

        self._sum_section = _CollapsibleSection("Summary", body, expanded=True)
        self._sum_title = QLabel(self)
        self._save_tsv = QPushButton("Save TSV", self)
        self._save_tsv.clicked.connect(self._save_summary)
        sum_head = QHBoxLayout()
        sum_head.addWidget(self._sum_title, 1)
        sum_head.addWidget(self._save_tsv)
        self._table = QTableWidget(0, 4, self)
        self._table.setHorizontalHeaderLabels(["Chr", "Probe", "Position", "Seen"])
        self._table.verticalHeader().setVisible(False)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        hh.setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self._table.setShowGrid(False)
        self._show_idio = QPushButton("Show on idiogram", self)
        self._show_idio.setFlat(True)
        self._show_idio.setObjectName("cytoLink")
        self._show_idio.clicked.connect(self.idiogramRequested)
        self._sum_section.body_layout().addLayout(sum_head)
        self._sum_section.body_layout().addWidget(self._table)
        self._sum_section.body_layout().addWidget(self._show_idio, 0, Qt.AlignmentFlag.AlignLeft)
        side.addWidget(self._sum_section)
        side.addStretch(1)

        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(body)
        scroll.setFixedWidth(SIDE_W)

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addLayout(left, 1)
        root.addWidget(scroll)

        self._selected: tuple[str, str] | None = None
        self._save_timer = QTimer(self)
        self._save_timer.setSingleShot(True)
        self._save_timer.setInterval(SAVE_DELAY_MS)
        self._save_timer.timeout.connect(self.flush)
        self._set_tool("number")
        self._show_empty()

    # ------------------------------------------------------------------ loading

    def refresh(self, folder: Path | None) -> None:
        self._folder = Path(folder) if folder else None
        if self._saving or self._folder is None:
            return
        s = self._session
        if s is None:
            return
        if not any(i.id == s.image.id for i in ci.load_images(self._folder)):
            self._session = None
            self._show_empty()
        elif not (s.dirty or self._save_timer.isActive()):
            self.show_image(s.image.id, keep_view=True)

    def image_id(self) -> str:
        return self._session.image.id if self._session else ""

    def follow(self, image_id: str) -> None:
        if not image_id or image_id == self.image_id():
            self._pending = ""
            return
        if self.isVisible():
            self.show_image(image_id)
        else:
            self._pending = image_id

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        pending, self._pending = self._pending, ""
        if pending and pending != self.image_id():
            self.show_image(pending)

    def hideEvent(self, event) -> None:  # noqa: N802
        self.flush()
        super().hideEvent(event)

    def show_image(self, image_id: str, *, keep_view: bool = False) -> None:
        self.flush()
        if self._folder is None:
            return
        image = next((i for i in ci.load_images(self._folder) if i.id == image_id), None)
        if image is None:
            self._session = None
            self._show_empty()
            return
        previous = self._session.image.id if self._session else ""
        self._session = ss.SignalsSession(self._folder, image)
        self._canvas.session = self._session
        if previous != image_id:
            self._selected = None
            self._pending_label = ""
        self._canvas.selected = self._selected
        message = ""
        qimage = None
        try:
            qimage = rgb_to_qimage(self._session.view.merged()) if self._session.view.is_multichannel else None
            if qimage is None:
                from PySide6.QtGui import QImage

                qimage = QImage(str(ci.image_file_path(self._folder, image)))
                if qimage.isNull():
                    qimage, message = None, "The picture file is missing."
                else:
                    self._session.set_picture(_qimage_rgb(qimage))
        except Exception as exc:  # noqa: BLE001
            message = f"Cannot show this image:\n{exc}"
        self._canvas.set_image(qimage, keep_view=keep_view and previous == image_id, message=message)
        self._fill_probes()
        self._fill_karyotypes()
        self._sync_all()
        if previous != image_id:
            self.imageChanged.emit(image_id)

    def _show_empty(self) -> None:
        self._canvas.session = None
        self._canvas.set_image(None, keep_view=False, message="No image selected. Choose Mark signals in the Library.")
        self._clear_layout(self._sel_lay)
        self._table.setRowCount(0)
        self._sum_title.setText("")
        self._chrom_count.setText("")
        self._clear_cells()

    # ------------------------------------------------------------------ tools and canvas events

    def _set_tool(self, tool: str) -> None:
        self._tool = tool if tool in dict(ss.TOOLS) else "number"
        for b, (key, _text) in zip(self._tool_group.buttons(), ss.TOOLS):
            b.setChecked(key == self._tool)
        self._canvas.set_tool_cursor(Qt.CursorShape.CrossCursor)
        self._update_hint()
        self._sync_labels()

    def _update_hint(self) -> None:
        s = self._session
        if s is None:
            self._hint.setText("")
            return
        if self._tool == "number":
            if s.karyotype is None:
                self._hint.setText("Make a karyotype first: Chromosomes > Karyotype > New karyotype…")
                return
            nxt = self._pending_label if self._pending_label and self._pending_label not in s.used_labels() else s.next_label()
            self._hint.setText(ss.TOOL_HINTS["number"].format(next=nxt) if nxt else
                               "All chromosomes are named. Drag a name, a telomere or the centromere to correct it.")
        else:
            probe = self._probe.currentData()
            self._hint.setText(ss.TOOL_HINTS["signal"].format(probe=probe) if probe else
                               "Add a probe first: Library > Edit lists > Probes.")

    def _tolerance(self, pixels: float) -> float:
        return pixels / max(self._canvas.zoom(), 0.05)

    def _handle_at(self, x: float, y: float) -> tuple[str, str] | None:
        """(mark id, which) of a handle of the selected chromosome under the pointer."""
        s = self._session
        if s is None or not self._selected or self._selected[0] != "mark":
            return None
        mark = s.details.chromosome(self._selected[1])
        if mark is None:
            return None
        tol = self._tolerance(9)
        for which in ("centromere", "telomere_p", "telomere_q"):
            pt = getattr(mark, which)
            if pt is not None and math.hypot(pt[0] - x, pt[1] - y) <= tol:
                return mark.id, which
        return None

    def _on_press(self, x: float, y: float) -> None:
        s = self._session
        if s is None:
            return
        self._moved = False
        handle = self._handle_at(x, y)
        if handle is not None:
            self._drag = ("handle", handle[0], handle[1])
            return
        hit = s.hit(x, y, max(self._tolerance(ss.HIT_RADIUS * 0.8), 6.0))
        if hit is not None:
            self._select(hit)
            self._drag = hit
            return
        self._drag = None
        if self._tool == "number":
            self._number_at(x, y)
        else:
            self._signal_at(x, y)

    def _on_move(self, x: float, y: float) -> None:
        s = self._session
        if s is None or not self._drag:
            return
        self._moved = True
        if self._drag[0] == "handle":
            s.set_point(self._drag[1], self._drag[2], (x, y))
        elif self._drag[0] == "mark":
            s.move_mark(self._drag[1], x, y)
        else:
            s.move_signal(self._drag[1], x, y)
        self._canvas.update()

    def _on_release(self, x: float, y: float) -> None:
        if self._drag:
            self._drag = None
            if getattr(self, "_moved", False):
                self._changed()

    def _on_menu(self, x: float, y: float, global_pos) -> None:
        from PySide6.QtWidgets import QMenu

        s = self._session
        if s is None:
            return
        hit = s.hit(x, y, max(self._tolerance(ss.HIT_RADIUS * 0.8), 6.0))
        if hit is None:
            return
        self._select(hit)
        menu = QMenu(self)
        if hit[0] == "signal":
            sig = s.signal(hit[1])
            mark = s.details.chromosome(sig.chromosome_id) if sig.chromosome_id else None
            head = menu.addAction(f"Signal on {mark.label if mark else 'no chromosome'} · {sig.probe}")
            head.setEnabled(False)
            sub = menu.addMenu("Change probe")
            for probe, dye in s.probes():
                act = sub.addAction(f"{probe} · {dye}" if dye else probe)
                act.setCheckable(True)
                act.setChecked(probe == sig.probe)
                act.triggered.connect(lambda _c=False, pr=probe: self._edit_signal(sig.id, probe=pr))
            sub = menu.addMenu("Move to chromosome")
            for m in sorted(s.details.chromosomes, key=lambda m: ss.natural_key(m.label)):
                act = sub.addAction(m.label)
                act.setCheckable(True)
                act.setChecked(m.id == sig.chromosome_id)
                act.triggered.connect(lambda _c=False, mid=m.id: self._edit_signal(sig.id, chromosome_id=mid))
            sub = menu.addMenu("Strength")
            for value in ss.STRENGTHS:
                act = sub.addAction(value.capitalize())
                act.setCheckable(True)
                act.setChecked(value == sig.strength)
                act.triggered.connect(lambda _c=False, v=value: self._edit_signal(sig.id, strength=v))
            menu.addSeparator()
            menu.addAction("Delete signal", self._delete_selected)
        else:
            mark = s.details.chromosome(hit[1])
            head = menu.addAction(f"Chromosome {mark.label}")
            head.setEnabled(False)
            sub = menu.addMenu("Rename")
            for label in s.labels():
                if label == mark.label or s.mark_for_label(label) is None:
                    act = sub.addAction(label)
                    act.setCheckable(True)
                    act.setChecked(label == mark.label)
                    act.triggered.connect(lambda _c=False, lb=label: self._relabel(mark.id, lb))
            menu.addAction("Find the axis again", lambda: (s.auto_axis(mark.id), self._changed()))
            menu.addAction("Swap p and q", lambda: (s.swap_arms(mark.id), self._changed()))
            menu.addSeparator()
            menu.addAction("Remove chromosome", self._delete_selected)
        menu.exec(global_pos)

    def _step_cell(self, step: int) -> None:
        s = self._session
        if s is None:
            return
        cells = s.slide_cells()
        ids = [c.image_id for c in cells]
        if s.image.id in ids:
            k = ids.index(s.image.id) + step
            if 0 <= k < len(ids):
                self.show_image(ids[k])

    def _number_at(self, x: float, y: float) -> None:
        s = self._session
        if s is None:
            return
        if s.karyotype is None:
            self._hint.setText("Make a karyotype first (Chromosomes > Karyotype).")
            return
        label = self._pending_label if self._pending_label and self._pending_label not in s.used_labels() else s.next_label()
        if not label:
            self._hint.setText("All chromosomes of this karyotype are numbered.")
            return
        mark = s.add_mark(label, x, y)
        s.auto_axis(mark.id)
        self._pending_label = ""
        self._select(("mark", mark.id))
        self._changed()

    def _signal_at(self, x: float, y: float) -> None:
        s = self._session
        if s is None:
            return
        probe = self._probe.currentData()
        if not probe:
            self._update_hint()
            return
        sig = s.add_signal(probe, x, y, refine=self._auto.isChecked())
        self._select(("signal", sig.id))
        self._changed()

    def _delete_selected(self) -> None:
        s = self._session
        if s is None or not self._selected:
            return
        kind, ident = self._selected
        (s.remove_mark if kind == "mark" else s.remove_signal)(ident)
        self._select(None)
        self._changed()

    # ------------------------------------------------------------------ side panel

    def _select(self, hit: tuple[str, str] | None) -> None:
        self._selected = hit
        self._canvas.selected = hit
        self._canvas.update()
        self._sync_selected()
        self._sync_labels()

    def _changed(self) -> None:
        self._canvas.update()
        self._sync_all()
        self._save_timer.start()

    def _sync_all(self) -> None:
        self._update_hint()
        self._sync_labels()
        self._sync_selected()
        self._sync_summary()
        self._sync_cells()

    def _fill_probes(self) -> None:
        s = self._session
        self._probe.clear()
        for probe, dye in (s.probes() if s else []):
            self._probe.addItem(f"{probe} · {dye}" if dye else probe, probe)
        if self._probe.count() == 0:
            self._probe.addItem("no probe set", "")

    def _fill_karyotypes(self) -> None:
        s = self._session
        self._karyo.blockSignals(True)
        self._karyo.clear()
        if s is not None:
            for k in s.karyotypes:
                self._karyo.addItem(f"{k.name} {k.formula}", k.name)
            self._karyo.addItem(NEW_KARYOTYPE, "")
            if s.karyotype is not None:
                idx = self._karyo.findData(s.karyotype.name)
                self._karyo.setCurrentIndex(max(idx, 0))
            else:
                self._karyo.setCurrentIndex(-1)
        self._karyo.blockSignals(False)
        self._edit_karyo.setEnabled(s is not None and s.karyotype is not None)

    def _on_karyotype(self, index: int) -> None:
        s = self._session
        if s is None:
            return
        name = self._karyo.itemData(index)
        if not name:
            self._open_karyotype(edit=False)
            return
        s.set_karyotype(next((k for k in s.karyotypes if k.name == name), None))
        self._changed()

    def _open_karyotype(self, *, edit: bool) -> None:
        s = self._session
        if s is None or self._folder is None:
            return
        from repeatviewer.views.cyto_karyotype_dialog import KaryotypeDialog

        dialog = KaryotypeDialog(self._folder, s.karyotype.name if edit and s.karyotype else "", self)
        accepted = dialog.exec()
        s.karyotypes = list_karyotypes(self._folder)
        if accepted and dialog.saved_name():
            s.karyotype = next((k for k in s.karyotypes if k.name == dialog.saved_name()), s.karyotype)
            s.details.karyotype = s.karyotype.name if s.karyotype else ""
            s.dirty = True
            self._save_timer.start()
            self.karyotypesChanged.emit()
        self._fill_karyotypes()
        self._sync_all()
        self._canvas.update()

    def _sync_labels(self) -> None:
        s = self._session
        clear_layout(self._labels_grid)
        self._label_buttons = {}
        if s is None:
            return
        labels = s.labels()
        used = s.used_labels()
        selected_label = ""
        if self._selected and self._selected[0] == "mark":
            m = s.details.chromosome(self._selected[1])
            selected_label = m.label if m else ""
        nxt = self._pending_label if self._pending_label and self._pending_label not in used else s.next_label()
        for n, label in enumerate(labels):
            b = QPushButton(label, self)
            on = label == selected_label
            if on:
                css = "background:#1F4E9A; color:#FFFFFF; border:1px solid #1F4E9A;"
            elif label == nxt and self._tool == "number":
                css = "background:#FFFFFF; color:#1C1C1E; border:2px solid #1F4E9A;"
            elif label in used:
                css = "background:#EEF2FA; color:#1C1C1E; border:1px solid #C9D3E6;"
            else:
                css = "background:#FFFFFF; color:#8A8D93; border:1px dashed #C9D3E6;"
            b.setStyleSheet(f"QPushButton {{ {css} border-radius:5px; padding:4px 0; font-size:12px; }}")
            b.clicked.connect(lambda _c=False, l=label: self._on_label_clicked(l))
            self._labels_grid.addWidget(b, n // 8, n % 8)
            self._label_buttons[label] = b
        self._chrom_count.setText(f"{len(used & set(labels))} of {len(labels)}")

    def _on_label_clicked(self, label: str) -> None:
        s = self._session
        if s is None:
            return
        mark = s.mark_for_label(label)
        if mark is not None:
            self._select(("mark", mark.id))
        else:
            self._pending_label = label
            self._set_tool("number")
            self._hint.setText(f"Click chromosome {label} in the picture")
            self._sync_labels()

    def _clear_layout(self, layout) -> None:
        clear_layout(layout)

    def _sync_selected(self) -> None:
        self._clear_layout(self._sel_lay)
        s = self._session
        if s is None:
            return
        sel = self._selected
        if sel and sel[0] == "signal" and s.signal(sel[1]) is None:
            sel = self._selected = None
        if sel and sel[0] == "mark" and s.details.chromosome(sel[1]) is None:
            sel = self._selected = None
        if sel is None:
            self._sel_section._title.setText("Selected")
            note = QLabel("Nothing selected. Use Select and click a chromosome or a signal.", self._sel_host)
            note.setWordWrap(True)
            self._sel_lay.addWidget(note)
            return
        if sel[0] == "signal":
            sig = s.signal(sel[1])
            mark = s.details.chromosome(sig.chromosome_id) if sig.chromosome_id else None
            self._sel_section._title.setText(f"Selected: signal on {mark.label}" if mark else "Selected: signal")
            form = QGridLayout()
            form.setColumnStretch(1, 1)
            probe = _combo([(f"{p} · {d}" if d else p, p) for p, d in s.probes()] or [(sig.probe, sig.probe)], sig.probe)
            chrom = _combo([("(none)", "")] + [(m.label, m.id) for m in s.details.chromosomes], sig.chromosome_id)
            arm = _combo([("(not set)", "")] + [(v, k) for k, v in ss.ARM_LABELS.items()], sig.arm)
            pos = _combo([("(not set)", "")] + [(v.capitalize(), v) for v in ss.POSITION_CHOICES], sig.position)
            strength = _combo([("(not set)", "")] + [(v.capitalize(), v) for v in ss.STRENGTHS], sig.strength)
            rows = (("Probe", probe), ("Chromosome", chrom), ("Arm", arm), ("Position", pos), ("Strength", strength))
            for r, (text, widget) in enumerate(rows):
                form.addWidget(QLabel(text, self._sel_host), r, 0)
                form.addWidget(widget, r, 1)
            dist = QLabel(s.distance_text(sig), self._sel_host)
            dist.setWordWrap(True)
            dist.setStyleSheet("color: palette(dark);")
            form.addWidget(QLabel("Distance", self._sel_host), len(rows), 0, Qt.AlignmentFlag.AlignTop)
            form.addWidget(dist, len(rows), 1)
            self._sel_lay.addLayout(form)
            probe.activated.connect(lambda _i: self._edit_signal(sig.id, probe=probe.currentData()))
            chrom.activated.connect(lambda _i: self._edit_signal(sig.id, chromosome_id=chrom.currentData()))
            arm.activated.connect(lambda _i: self._edit_signal(sig.id, arm=arm.currentData()))
            pos.activated.connect(lambda _i: self._edit_signal(sig.id, position=pos.currentData()))
            strength.activated.connect(lambda _i: self._edit_signal(sig.id, strength=strength.currentData()))
            remove = QPushButton("Remove signal", self._sel_host)
            remove.clicked.connect(self._delete_selected)
            self._sel_lay.addWidget(remove, 0, Qt.AlignmentFlag.AlignLeft)
        else:
            mark = s.details.chromosome(sel[1])
            self._sel_section._title.setText(f"Selected: chromosome {mark.label}")
            free = [l for l in s.labels() if l == mark.label or s.mark_for_label(l) is None]
            label = _combo([(l, l) for l in free] or [(mark.label, mark.label)], mark.label)
            form = QGridLayout()
            form.setColumnStretch(1, 1)
            form.addWidget(QLabel("Name", self._sel_host), 0, 0)
            form.addWidget(label, 0, 1)
            measured = "axis measured" if mark.centromere is not None else "not measured: use Measure to set the telomeres"
            note = QLabel(measured, self._sel_host)
            note.setStyleSheet("color: palette(dark);")
            form.addWidget(QLabel("Axis", self._sel_host), 1, 0)
            form.addWidget(note, 1, 1)
            self._sel_lay.addLayout(form)
            label.activated.connect(lambda _i: self._relabel(mark.id, label.currentData()))
            remove = QPushButton("Remove chromosome", self._sel_host)
            remove.clicked.connect(self._delete_selected)
            self._sel_lay.addWidget(remove, 0, Qt.AlignmentFlag.AlignLeft)

    def _edit_signal(self, signal_id: str, **fields) -> None:
        s = self._session
        if s is None:
            return
        s.update_signal(signal_id, **fields)
        self._changed()

    def _relabel(self, mark_id: str, label: str) -> None:
        s = self._session
        if s is None:
            return
        try:
            s.relabel_mark(mark_id, label)
        except ValueError as exc:
            QMessageBox.information(self, "Chromosome name", str(exc))
        self._changed()

    # ------------------------------------------------------------------ summary and cells

    def _sync_summary(self) -> None:
        s = self._session
        if s is None:
            return
        rows, cells = s.summary()
        self._sum_title.setText(f"<b>Summary over {cells} cell{'s' if cells != 1 else ''}</b>")
        unsig = s.unsignalled()
        self._table.setRowCount(len(rows) + (1 if unsig else 0))
        for r, row in enumerate(rows):
            where = f"{row.arm} {row.position}" if row.arm in ("p", "q") else row.position
            for c, text in enumerate((row.chromosome, row.probe, where, f"{row.seen} / {row.total}")):
                item = QTableWidgetItem(text)
                if c == 3:
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                if c == 1:
                    item.setForeground(QColor(s.probe_colour(row.probe)))
                self._table.setItem(r, c, item)
        if unsig:
            r = len(rows)
            for c, text in enumerate((", ".join(unsig), "none", "", "")):
                self._table.setItem(r, c, QTableWidgetItem(text))
        self._table.resizeRowsToContents()
        rows_h = sum(self._table.rowHeight(r) for r in range(self._table.rowCount())) or 26
        self._table.setFixedHeight(self._table.horizontalHeader().sizeHint().height() + rows_h + 4)
        self._table.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)

    def _save_summary(self) -> None:
        s = self._session
        if s is None or self._folder is None:
            return
        rows, _n = s.summary()
        path, _ = QFileDialog.getSaveFileName(self, "Save summary", str(self._folder / "signal_summary.tsv"), "TSV (*.tsv)")
        if not path:
            return
        try:
            Path(path).write_text(sg.summary_tsv(rows), encoding="utf-8")
        except OSError as exc:
            QMessageBox.warning(self, "Save summary", str(exc))

    def _clear_cells(self) -> None:
        self._clear_layout(self._cells_row)

    def _sync_cells(self) -> None:
        s = self._session
        self._clear_cells()
        if s is None:
            return
        cells = s.slide_cells()
        slide = next((sl.name for sl in load_slides(self._folder) if sl.id == s.details.slide_id), "")
        lead = QLabel(f"Cells in {slide}" if slide else "Cells", self._cells_host)
        lead.setObjectName("cytoMuted")
        ids = [c.image_id for c in cells]
        pos = ids.index(s.image.id) + 1 if s.image.id in ids else 0
        self._nav.setText(f"cell {cells[pos - 1].label} of {len(cells)}" if pos else "")
        self._prev.setEnabled(pos > 1)
        self._next.setEnabled(0 < pos < len(cells))
        self._cells_row.addWidget(lead)
        for cell in cells:
            b = QPushButton(cell.label, self._cells_host)
            b.setCheckable(True)
            b.setChecked(cell.image_id == s.image.id)
            b.setToolTip("marked" if cell.marked else "not marked yet")
            if not cell.marked:
                f = b.font()
                f.setItalic(True)
                b.setFont(f)
            b.clicked.connect(lambda _c=False, i=cell.image_id: self.show_image(i))
            self._cells_row.addWidget(b)
        done = sum(1 for c in cells if c.marked)
        self._cells_row.addWidget(QLabel(f"{done} of {len(cells)} marked", self._cells_host))

    # ------------------------------------------------------------------ saving

    def flush(self) -> None:
        self._save_timer.stop()
        s = self._session
        if s is None or self._folder is None or not s.dirty:
            return
        try:
            s.save()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Signals", f"Could not save the marks:\n{exc}")
            return
        from repeatviewer.widgets.cyto_images_ui import notify_changed

        self._saving = True
        try:
            notify_changed(self._folder)
        finally:
            self._saving = False
