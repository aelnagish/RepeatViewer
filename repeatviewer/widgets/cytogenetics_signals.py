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
    "Number the chromosomes of one cell and mark where each probe signal sits.<br>"
    "<b>Number chromosomes</b>: click the chromosomes in order (1a, 1b, 1c, 2a …). The next name is taken from the karyotype; "
    "click a name on the right to choose another one.<br>"
    "<b>Add signal</b>: click a signal. It is moved to the centre of the bright spot (if \"Find signals automatically\" is on) "
    "and given to the nearest chromosome.<br>"
    "<b>Measure</b>: drag from the p (short arm) telomere to the q telomere. The centromere is placed from the centromere index "
    "of the karyotype. Arm and position of signals on that chromosome are then worked out.<br>"
    "<b>Select</b>: click a chromosome or a signal to edit it, drag to move it, press Delete to remove it.<br>"
    "Middle or right mouse button moves the picture, the wheel zooms.<br>"
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
        if event.button() == Qt.MouseButton.LeftButton:
            p = self.image_point(event.position())
            self.pressedAt.emit(p.x(), p.y())
        elif event.button() in (Qt.MouseButton.MiddleButton, Qt.MouseButton.RightButton):
            self._panning = True
            self.start_pan(event.position())

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
        self._tool = "select"
        self._pending_label = ""
        self._drag: tuple[str, str] | None = None
        self._measure_start: QPointF | None = None
        self._measure_mark = ""
        self._saving = False
        self._pending = ""
        self._label_buttons: dict[str, QPushButton] = {}

        # ---- toolbar
        self._canvas = MarkCanvas(self)
        self._canvas.pressedAt.connect(self._on_press)
        self._canvas.movedTo.connect(self._on_move)
        self._canvas.releasedAt.connect(self._on_release)
        self._canvas.deleteRequested.connect(self._delete_selected)
        self._tool_group = QButtonGroup(self)
        self._tool_group.setExclusive(True)
        bar = QHBoxLayout()
        bar.setContentsMargins(10, 8, 10, 6)
        for key, text in ss.TOOLS:
            b = QPushButton(text, self)
            b.setCheckable(True)
            b.setChecked(key == self._tool)
            b.clicked.connect(lambda _c=False, k=key: self._set_tool(k))
            self._tool_group.addButton(b)
            bar.addWidget(b)
        bar.addSpacing(10)
        bar.addWidget(QLabel("Probe", self))
        self._probe = NoWheelComboBox(self)
        self._probe.setMinimumWidth(130)
        bar.addWidget(self._probe)
        self._auto = QCheckBox("Find signals automatically", self)
        self._auto.setChecked(_read("signals_auto", "true") != "false")
        self._auto.toggled.connect(lambda on: _write("signals_auto", "true" if on else "false"))
        bar.addWidget(self._auto)
        bar.addStretch(1)
        self._hint = QLabel(self)
        self._hint.setStyleSheet("color: palette(dark);")
        bar.addWidget(self._hint)
        bar.addWidget(make_help_button(self, "Signals", HELP_TEXT))

        # ---- cells row
        self._cells_row = QHBoxLayout()
        self._cells_row.setContentsMargins(10, 6, 10, 6)
        self._cells_host = QWidget(self)
        self._cells_host.setLayout(self._cells_row)

        left = QVBoxLayout()
        left.setContentsMargins(0, 0, 0, 0)
        left.setSpacing(0)
        left.addLayout(bar)
        left.addWidget(self._canvas, 1)
        left.addWidget(self._cells_host)

        # ---- side panel
        body = QWidget(self)
        side = QVBoxLayout(body)
        side.setContentsMargins(10, 4, 10, 8)
        side.setSpacing(12)

        self._chrom_section = _CollapsibleSection("Chromosomes", body, expanded=True)
        self._chrom_count = QLabel(self)
        self._chrom_count.setStyleSheet("color: palette(dark);")
        self._karyo = NoWheelComboBox(self)
        self._karyo.activated.connect(self._on_karyotype)
        self._edit_karyo = QPushButton("Edit karyotype", self)
        self._edit_karyo.setFlat(True)
        self._edit_karyo.setStyleSheet("QPushButton { color:#1F4E9A; border:none; } QPushButton:hover { text-decoration: underline; }")
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
        self._table.horizontalHeader().setSectionResizeMode(QHeaderView.ResizeMode.Stretch)
        self._table.setShowGrid(False)
        self._show_idio = QPushButton("Show on idiogram", self)
        self._show_idio.setFlat(True)
        self._show_idio.setStyleSheet("QPushButton { color:#1F4E9A; border:none; } QPushButton:hover { text-decoration: underline; }")
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
        self._set_tool("select")
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
        self._tool = tool
        for b, (key, _text) in zip(self._tool_group.buttons(), ss.TOOLS):
            b.setChecked(key == tool)
        self._hint.setText(ss.TOOL_HINTS[tool] + (" (p end first)" if tool == "measure" else ""))
        self._canvas.set_tool_cursor(Qt.CursorShape.ArrowCursor if tool == "select" else Qt.CursorShape.CrossCursor)

    def _on_press(self, x: float, y: float) -> None:
        s = self._session
        if s is None:
            return
        tol = ss.HIT_RADIUS / max(self._canvas.zoom(), 0.05) * 0.8
        if self._tool == "select":
            hit = s.hit(x, y, max(tol, 6.0))
            self._select(hit)
            self._drag = hit
        elif self._tool == "number":
            self._number_at(x, y)
        elif self._tool == "signal":
            self._signal_at(x, y)
        elif self._tool == "measure":
            hit = s.hit(x, y, max(tol * 2, 10.0))
            mark_id = hit[1] if hit and hit[0] == "mark" else (self._selected[1] if self._selected and self._selected[0] == "mark" else "")
            if not mark_id:
                self._hint.setText("Number the chromosome first (Number chromosomes), then measure it.")
                return
            self._measure_mark = mark_id
            self._measure_start = QPointF(x, y)
            self._canvas.rubber = (QPointF(x, y), QPointF(x, y))

    def _on_move(self, x: float, y: float) -> None:
        s = self._session
        if s is None:
            return
        if self._tool == "select" and self._drag:
            kind, ident = self._drag
            if kind == "mark":
                s.move_mark(ident, x, y)
            else:
                s.move_signal(ident, x, y)
            self._canvas.update()
        elif self._tool == "measure" and self._measure_start is not None:
            self._canvas.rubber = (self._measure_start, QPointF(x, y))
            self._canvas.update()

    def _on_release(self, x: float, y: float) -> None:
        s = self._session
        if s is None:
            return
        if self._tool == "select" and self._drag:
            self._drag = None
            self._changed()
        elif self._tool == "measure" and self._measure_start is not None:
            a, b = self._measure_start, QPointF(x, y)
            self._measure_start = None
            self._canvas.rubber = None
            if math.hypot(b.x() - a.x(), b.y() - a.y()) >= 5:
                s.set_axis(self._measure_mark, (a.x(), a.y()), (b.x(), b.y()))
                self._select(("mark", self._measure_mark))
                self._changed()
            else:
                self._canvas.update()

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
        self._pending_label = ""
        self._select(("mark", mark.id))
        self._changed()
        nxt = s.next_label()
        self._hint.setText(f"Next: {nxt}" if nxt else "All chromosomes numbered.")

    def _signal_at(self, x: float, y: float) -> None:
        s = self._session
        if s is None:
            return
        probe = self._probe.currentData()
        if not probe:
            self._hint.setText("Set a probe on a channel in the Image tab first.")
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
        while self._labels_grid.count():
            w = self._labels_grid.takeAt(0).widget()
            if w is not None:
                w.deleteLater()
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
            elif label in used:
                css = "background:#CFE0D4; color:#1C1C1E; border:1px solid #9DBBA5;"
            elif label == nxt and self._tool == "number":
                css = "background:#EEF2FA; color:#1C1C1E; border:2px solid #1F4E9A;"
            else:
                css = "background:#EEF2FA; color:#1C1C1E; border:1px solid #C9D3E6;"
            b.setStyleSheet(f"QPushButton {{ {css} border-radius:5px; padding:4px 0; font-size:12px; }}")
            b.clicked.connect(lambda _c=False, l=label: self._on_label_clicked(l))
            self._labels_grid.addWidget(b, n // 8, n % 8)
            self._label_buttons[label] = b
        self._chrom_count.setText(f"{len(used & set(labels))} of {len(labels)} numbered")

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
        while layout.count():
            item = layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
            elif item.layout() is not None:
                self._clear_layout(item.layout())

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
        self._table.setFixedHeight(28 + 24 * max(1, self._table.rowCount()))

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
        lead.setStyleSheet("color: palette(dark);")
        self._cells_row.addWidget(lead)
        for cell in cells:
            b = QPushButton(cell.label, self._cells_host)
            current = cell.image_id == s.image.id
            if current:
                css = "background:#E3EBF8; color:#1C1C1E; border:1px solid #1F4E9A;"
            elif cell.marked:
                css = "background:#FFFFFF; color:#1C1C1E; border:1px solid #C4C4CA;"
            else:
                css = "background:#F1F1F3; color:#5A5D63; border:1px solid #C4C4CA;"
            b.setStyleSheet(f"QPushButton {{ {css} border-radius:5px; padding:3px 8px; font-size:12px; }}")
            b.clicked.connect(lambda _c=False, i=cell.image_id: self.show_image(i))
            self._cells_row.addWidget(b)
        self._cells_row.addStretch(1)
        done = sum(1 for c in cells if c.marked)
        self._cells_row.addWidget(QLabel(f"{done} of {len(cells)} cells marked", self._cells_host))

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
