"""Image tab of the Cytogenetics workspace: view one cell, set channels, levels and the scale bar."""

from __future__ import annotations

from pathlib import Path

import numpy as np
from PySide6.QtCore import QPointF, QRectF, QSettings, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QImage, QPainter, QPen
from PySide6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSlider,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import image_session as isn
from repeatviewer.backends.cytogenetics.lists import load_lists
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cytogenetics_style import clear_layout
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.sequence_annotation_rail import _CollapsibleSection

_SETTINGS = "cytogenetics_workspace"
SIDE_W = 320
SAVE_DELAY_MS = 800
RENDER_DELAY_MS = 30

HELP_TEXT = (
    "<b>Image</b><br>"
    "One cell. The raw channel TIFFs are never changed: levels, gamma and colours are display settings.<br>"
    "<b>Mouse wheel</b> zooms, <b>drag</b> moves the picture.<br>"
    "<b>Channels</b>: tick to show or hide, click the colour square to change the colour, "
    "set the dye and the probe, move the sliders for levels (black point and white point) and gamma.<br>"
    "<b>Gamma</b> other than 1.0 is a non-linear change: say so in the figure legend (the log below tells you).<br>"
    "<b>Scale bar</b> needs the pixel size; it is read from the TIFF when the file has it.<br>"
    "<b>Adjustment log</b> lists every change for the figure legend. "
    "<b>Apply to whole slide</b> gives all other cells of the slide the same colours, levels and gamma."
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


def _cell_name(image: ci.CytoImage, number: int) -> str:
    import re

    m = re.search(r"(\d+)\D*$", image.original_name or "")
    return m.group(1).zfill(2) if m else str(number).zfill(2)


def rgb_to_qimage(rgb: np.ndarray) -> QImage:
    rgb = np.ascontiguousarray(rgb, dtype=np.uint8)
    h, w = rgb.shape[:2]
    return QImage(rgb.data, w, h, 3 * w, QImage.Format.Format_RGB888).copy()


class ImageCanvas(QWidget):
    """Zoom and pan view with an optional scale bar."""

    zoomChanged = Signal(float)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumSize(320, 240)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self._image: QImage | None = None
        self._message = ""
        self._zoom = 1.0
        self._fit = True
        self._off = QPointF(0, 0)
        self._drag: QPointF | None = None
        self._bar: tuple[float, str] | None = None  # (length in image pixels, label)

    def set_image(self, image: QImage | None, *, keep_view: bool, message: str = "") -> None:
        self._message = message
        new_size = image is not None and (self._image is None or image.size() != self._image.size())
        self._image = image
        if image is not None and (new_size or not keep_view):
            self.fit()
        self.update()

    def set_scale_bar(self, bar: tuple[float, str] | None) -> None:
        self._bar = bar
        self.update()

    def zoom(self) -> float:
        return self._zoom

    def fit(self) -> None:
        self._fit = True
        self._apply_fit()

    def _apply_fit(self) -> None:
        if self._image is None or self._image.isNull():
            return
        z = min(self.width() / self._image.width(), self.height() / self._image.height())
        self._zoom = max(0.02, z)
        self._off = QPointF(0, 0)
        self.zoomChanged.emit(self._zoom)
        self.update()

    def set_zoom(self, zoom: float, anchor: QPointF | None = None) -> None:
        if self._image is None:
            return
        zoom = min(max(zoom, 0.02), 32.0)
        anchor = anchor if anchor is not None else QPointF(self.width() / 2, self.height() / 2)
        origin = self._origin()
        ix, iy = (anchor.x() - origin.x()) / self._zoom, (anchor.y() - origin.y()) / self._zoom
        self._zoom = zoom
        self._fit = False
        self._off = QPointF(anchor.x() - ix * zoom, anchor.y() - iy * zoom)
        self.zoomChanged.emit(self._zoom)
        self.update()

    def _origin(self) -> QPointF:
        """Top-left of the picture on screen, centred or clamped."""
        if self._image is None:
            return QPointF(0, 0)
        w, h = self._image.width() * self._zoom, self._image.height() * self._zoom
        if self._fit:
            return QPointF((self.width() - w) / 2, (self.height() - h) / 2)
        x = (self.width() - w) / 2 if w <= self.width() else min(0.0, max(self.width() - w, self._off.x()))
        y = (self.height() - h) / 2 if h <= self.height() else min(0.0, max(self.height() - h, self._off.y()))
        return QPointF(x, y)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        if self._fit:
            self._apply_fit()

    def wheelEvent(self, event) -> None:  # noqa: N802
        steps = event.angleDelta().y() / 120.0
        if steps:
            self.set_zoom(self._zoom * (1.18 ** steps), event.position())
        event.accept()

    def image_point(self, pos: QPointF) -> QPointF:
        """Screen position -> image pixel."""
        origin = self._origin()
        return QPointF((pos.x() - origin.x()) / self._zoom, (pos.y() - origin.y()) / self._zoom)

    def screen_point(self, x: float, y: float) -> QPointF:
        origin = self._origin()
        return QPointF(origin.x() + x * self._zoom, origin.y() + y * self._zoom)

    def start_pan(self, pos: QPointF) -> None:
        self._drag = pos
        self.setCursor(Qt.CursorShape.ClosedHandCursor)

    def move_pan(self, pos: QPointF) -> None:
        if self._drag is not None and self._image is not None:
            delta = pos - self._drag
            self._drag = pos
            origin = self._origin()
            self._fit = False
            self._off = QPointF(origin.x() + delta.x(), origin.y() + delta.y())
            self.update()

    def end_pan(self) -> None:
        self._drag = None
        self.unsetCursor()

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self.start_pan(event.position())

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        self.move_pan(event.position())

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        self.end_pan()

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        self.fit()

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.fillRect(self.rect(), QColor("#0b0c10"))
        if self._image is None or self._image.isNull():
            p.setPen(QColor("#aaaaaa"))
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter | Qt.TextFlag.TextWordWrap,
                       self._message or "No image selected. Double-click an image in the Library.")
            return
        origin = self._origin()
        target = QRectF(origin.x(), origin.y(), self._image.width() * self._zoom, self._image.height() * self._zoom)
        p.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, self._zoom < 1.0)
        p.drawImage(target, self._image)
        if self._message:
            p.setPen(QColor("#ffb4b4"))
            p.drawText(self.rect().adjusted(12, 12, -12, -12), Qt.AlignmentFlag.AlignTop | Qt.TextFlag.TextWordWrap, self._message)
        if self._bar is not None:
            length = self._bar[0] * self._zoom
            if 4 <= length <= self.width() * 0.7:
                x1, y = self.width() - 20.0, self.height() - 22.0
                x0 = x1 - length
                p.fillRect(QRectF(x0 - 1, y - 1, length + 2, 6), QColor("#000000"))
                p.fillRect(QRectF(x0, y, length, 4), QColor("#ffffff"))
                font = QFont(self.font())
                font.setPointSizeF(max(8.0, font.pointSizeF()))
                p.setFont(font)
                label = self._bar[1]
                rect = QRectF(x0 - 40, y - 24, length + 80, 20)
                p.setPen(QPen(QColor("#000000")))
                p.drawText(rect.translated(1, 1), Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignBottom, label)
                p.setPen(QPen(QColor("#ffffff")))
                p.drawText(rect, Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignBottom, label)


_Canvas = ImageCanvas  # old name


def _slider(lo: int, hi: int, value: int) -> QSlider:
    s = QSlider(Qt.Orientation.Horizontal)
    s.setRange(lo, hi)
    s.setValue(max(lo, min(hi, value)))
    return s


class _ChannelCard(QFrame):
    """Controls of one channel."""

    changed = Signal()
    removeRequested = Signal(int)
    colourRequested = Signal(int)

    def __init__(self, session: isn.ImageSession, index: int, dyes: list[str], probes: list[str], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._s = session
        self._i = index
        c = session.channels[index]
        self.setFrameShape(QFrame.Shape.StyledPanel)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(8, 6, 8, 6)
        lay.setSpacing(4)

        top = QHBoxLayout()
        self._show = QCheckBox(self)
        self._show.setChecked(c.display.visible)
        self._show.setToolTip("Show this channel")
        self._show.toggled.connect(self._on_visible)
        self._dye = NoWheelComboBox(self)
        self._dye.setEditable(True)
        self._dye.addItems(dyes)
        self._dye.setCurrentText(c.dye)
        self._dye.setFixedWidth(78)
        self._dye.setToolTip("Dye of this channel")
        self._dye.lineEdit().editingFinished.connect(self._on_dye)
        name = QLabel(self)
        name.setStyleSheet("color: palette(dark); font-size: 11px;")
        name.setToolTip(c.original_name or c.file)
        name.setText(name.fontMetrics().elidedText(c.original_name or Path(c.file).name, Qt.TextElideMode.ElideMiddle, 96))
        name.setMinimumWidth(0)
        self._colour = QPushButton(self)
        self._colour.setFixedSize(24, 24)
        self._colour.setToolTip("Channel colour")
        self._colour.clicked.connect(lambda: self.colourRequested.emit(self._i))
        self.paint_colour()
        remove = QPushButton("×", self)
        remove.setStyleSheet("QPushButton { padding: 0px; }")
        remove.setFixedSize(24, 24)
        remove.setToolTip("Remove this channel from the image")
        remove.clicked.connect(lambda: self.removeRequested.emit(self._i))
        top.addWidget(self._show)
        top.addWidget(self._dye)
        top.addWidget(name, 1)
        top.addWidget(self._colour)
        top.addWidget(remove)
        lay.addLayout(top)

        self._probe = NoWheelComboBox(self)
        self._probe.setEditable(True)
        self._probe.addItems([""] + probes)
        self._probe.setCurrentText(c.probe)
        self._probe.setToolTip("Probe detected in this channel")
        self._probe.lineEdit().setPlaceholderText("probe in this channel")
        self._probe.lineEdit().editingFinished.connect(self._on_probe)
        lay.addWidget(self._probe)

        lo, hi = session.level_range(index)
        self._low = _slider(int(lo), int(hi), int(c.display.low))
        self._high = _slider(int(lo), int(hi), int(c.display.high or hi))
        self._gamma = _slider(5, 20, int(round(c.display.gamma * 10)))
        self._levels_val, self._low_val, self._gamma_val = QLabel(self), QLabel(self), QLabel(self)
        for text, slider, value in (("Levels", self._high, self._levels_val), ("Black", self._low, self._low_val),
                                    ("Gamma", self._gamma, self._gamma_val)):
            row = QHBoxLayout()
            label = QLabel(text, self)
            label.setFixedWidth(46)
            value.setFixedWidth(84 if text == "Levels" else 34)
            value.setStyleSheet("font-size: 11px;")
            value.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            row.addWidget(label)
            row.addWidget(slider, 1)
            row.addWidget(value)
            lay.addLayout(row)
        self._low.valueChanged.connect(self._on_levels)
        self._high.valueChanged.connect(self._on_levels)
        self._gamma.valueChanged.connect(self._on_gamma)
        self._sync_labels()

    def paint_colour(self) -> None:
        self._colour.setStyleSheet(
            f"background:{self._s.channels[self._i].display.color}; border:1px solid #888; border-radius:4px;"
        )

    def refresh_levels(self) -> None:
        d = self._s.channels[self._i].display
        for slider, value in ((self._low, d.low), (self._high, d.high), (self._gamma, d.gamma * 10)):
            slider.blockSignals(True)
            slider.setValue(int(round(value)))
            slider.blockSignals(False)
        self._sync_labels()

    def _sync_labels(self) -> None:
        self._low_val.setText(f"{self._low.value()}")
        self._levels_val.setText(f"{self._low.value()} to {self._high.value()}")
        self._gamma_val.setText(f"{self._gamma.value() / 10:.1f}")

    def _on_visible(self, on: bool) -> None:
        self._s.set_visible(self._i, on)
        self.changed.emit()

    def _on_dye(self) -> None:
        self._s.channels[self._i].dye = self._dye.currentText().strip()
        self.changed.emit()

    def _on_probe(self) -> None:
        self._s.channels[self._i].probe = self._probe.currentText().strip()
        self.changed.emit()

    def _on_levels(self) -> None:
        low, high = self._low.value(), self._high.value()
        if high <= low:  # keep white above black
            sender = self.sender()
            if sender is self._low:
                high = low + 1
                self._high.blockSignals(True)
                self._high.setValue(high)
                self._high.blockSignals(False)
            else:
                low = max(self._low.minimum(), high - 1)
                self._low.blockSignals(True)
                self._low.setValue(low)
                self._low.blockSignals(False)
        self._s.set_levels(self._i, low, high)
        self._sync_labels()
        self.changed.emit()

    def _on_gamma(self, value: int) -> None:
        self._s.set_gamma(self._i, value / 10.0)
        self._sync_labels()
        self.changed.emit()


class ImageTab(QWidget):
    """Image: canvas on the left, channels / scale bar / adjustment log on the right."""

    name = "Image"
    imageChanged = Signal(str)  # image id after a new image is shown

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._folder: Path | None = None
        self._session: isn.ImageSession | None = None
        self._cards: list[_ChannelCard] = []
        self._merge_mode = _read("image_mode", "merge") != "split"
        self._saving = False
        self._neighbours = (0, 0, "", "")
        self._image_w = 0
        self._pending = ""

        # toolbar
        self._canvas = ImageCanvas(self)
        self._canvas.zoomChanged.connect(self._on_zoom)
        out_btn, in_btn, fit_btn, one_btn = (QPushButton(t, self) for t in ("−", "+", "Fit", "100 %"))
        out_btn.clicked.connect(lambda: self._canvas.set_zoom(self._canvas.zoom() / 1.25))
        in_btn.clicked.connect(lambda: self._canvas.set_zoom(self._canvas.zoom() * 1.25))
        fit_btn.clicked.connect(self._canvas.fit)
        one_btn.clicked.connect(lambda: self._canvas.set_zoom(1.0))
        self._zoom_label = QLabel(self)
        self._zoom_label.setMinimumWidth(48)
        self._merge_btn = QPushButton("Merge", self)
        self._split_btn = QPushButton("Channels side by side", self)
        for b in (self._merge_btn, self._split_btn):
            b.setCheckable(True)
        self._merge_btn.clicked.connect(lambda: self._set_mode(True))
        self._split_btn.clicked.connect(lambda: self._set_mode(False))
        self._prev = QPushButton("◀", self)
        self._next = QPushButton("▶", self)
        self._prev.clicked.connect(lambda: self._go(2))
        self._next.clicked.connect(lambda: self._go(3))
        self._nav = QLabel(self)
        bar = QHBoxLayout()
        bar.setContentsMargins(10, 8, 10, 6)
        for w in (out_btn, in_btn, fit_btn, one_btn, self._zoom_label, self._merge_btn, self._split_btn):
            bar.addWidget(w)
        bar.addStretch(1)
        bar.addWidget(self._prev)
        bar.addWidget(self._nav)
        bar.addWidget(self._next)
        bar.addWidget(make_help_button(self, "Image", HELP_TEXT))

        left = QVBoxLayout()
        left.setContentsMargins(0, 0, 0, 0)
        left.setSpacing(0)
        left.addLayout(bar)
        left.addWidget(self._canvas, 1)

        # side panel
        self._side_body = QWidget(self)
        self._side_body.setObjectName("cytoAside")
        self._side_body.setMaximumWidth(SIDE_W)
        side = QVBoxLayout(self._side_body)
        side.setContentsMargins(10, 4, 10, 8)
        side.setSpacing(12)

        self._chan_section = _CollapsibleSection("Channels", self._side_body, expanded=True)
        self._chan_host = QWidget(self._chan_section)
        self._chan_lay = QVBoxLayout(self._chan_host)
        self._chan_lay.setContentsMargins(0, 0, 0, 0)
        self._chan_lay.setSpacing(6)
        self._chan_section.add_widget(self._chan_host)
        self._add_btn = QPushButton("+ Add channel (TIFF)", self)
        self._auto_btn = QPushButton("Auto levels (0.1 % saturated)", self)
        self._add_btn.clicked.connect(self._add_channel)
        self._auto_btn.clicked.connect(self._auto_all)
        self._chan_section.add_widget(self._add_btn)
        self._chan_section.add_widget(self._auto_btn)
        side.addWidget(self._chan_section)

        self._scale_section = _CollapsibleSection("Scale bar", self._side_body, expanded=True)
        self._scale_show = QCheckBox("Show", self)
        self._scale_show.toggled.connect(self._on_scale)
        self._pixel = QLineEdit(self)
        self._pixel.setPlaceholderText("µm per pixel")
        self._pixel.editingFinished.connect(self._on_pixel)
        self._length = NoWheelComboBox(self)
        self._length.addItems(isn.LENGTH_CHOICES)
        self._length.currentTextChanged.connect(self._on_scale)
        form = QVBoxLayout()
        form.addWidget(self._scale_show)
        for text, widget in (("Pixel size (µm)", self._pixel), ("Length", self._length)):
            r = QHBoxLayout()
            lab = QLabel(text, self)
            lab.setFixedWidth(100)
            r.addWidget(lab)
            r.addWidget(widget, 1)
            form.addLayout(r)
        self._scale_section.body_layout().addLayout(form)
        side.addWidget(self._scale_section)

        self._log_section = _CollapsibleSection("Adjustment log", self._side_body, expanded=True)
        self._log_label = QLabel(self)
        self._log_label.setWordWrap(True)
        self._log_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._slide_btn = QPushButton("Apply to whole slide", self)
        self._slide_btn.clicked.connect(self._apply_to_slide)
        self._log_section.add_widget(self._log_label)
        self._log_section.add_widget(self._slide_btn)
        side.addWidget(self._log_section)
        side.addStretch(1)

        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(self._side_body)
        scroll.setFixedWidth(SIDE_W)

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addLayout(left, 1)
        root.addWidget(scroll)

        self._render_timer = QTimer(self)
        self._render_timer.setSingleShot(True)
        self._render_timer.setInterval(RENDER_DELAY_MS)
        self._render_timer.timeout.connect(self._render)
        self._save_timer = QTimer(self)
        self._save_timer.setSingleShot(True)
        self._save_timer.setInterval(SAVE_DELAY_MS)
        self._save_timer.timeout.connect(self.flush)
        self._sync_mode_buttons()
        self._show_empty()

    # ------------------------------------------------------------------ loading

    def refresh(self, folder: Path | None) -> None:
        """Called by the workspace on every store change. Ignores changes this tab wrote itself."""
        self._folder = Path(folder) if folder else None
        if self._saving:
            return
        if self._session is not None and self._folder is not None:
            current = next((i for i in ci.load_images(self._folder) if i.id == self._session.image.id), None)
            if current is None:
                self._session = None
                self._show_empty()
            elif not self._save_timer.isActive():
                self.show_image(current.id, keep_view=True)

    def follow(self, image_id: str) -> None:
        """Another tab moved to this cell: show it now if visible, else when this tab is shown."""
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

    def image_id(self) -> str:
        return self._session.image.id if self._session else ""

    def show_image(self, image_id: str, *, keep_view: bool = False) -> None:
        self.flush()
        if self._folder is None:
            return
        image = next((i for i in ci.load_images(self._folder) if i.id == image_id), None)
        if image is None:
            self._session = None
            self._show_empty()
            return
        self._session = isn.ImageSession(self._folder, image)
        self._neighbours = isn.neighbours(self._folder, image_id)
        self._build_panel()
        self._render(keep_view=keep_view)
        self.imageChanged.emit(image_id)

    def _show_empty(self) -> None:
        self._canvas.set_image(None, keep_view=False)
        self._canvas.set_scale_bar(None)
        self._nav.setText("")
        self._prev.setEnabled(False)
        self._next.setEnabled(False)
        for w in (self._add_btn, self._auto_btn, self._slide_btn):
            w.setEnabled(False)
        self._clear_cards()
        self._log_label.setText("")

    def _clear_cards(self) -> None:
        self._cards = []
        clear_layout(self._chan_lay)

    # ------------------------------------------------------------------ side panel

    def _build_panel(self) -> None:
        s = self._session
        assert s is not None
        self._clear_cards()
        lists = load_lists(self._folder) if self._folder else None
        dyes = [d.get("name", "") for d in lists.get("dyes", [])] if lists else []
        probes = [p.get("name", "") for p in lists.get("probes", [])] if lists else []
        if s.is_multichannel:
            from repeatviewer.backends.cytogenetics.figure import ordered_channels

            order = [s.channels.index(c) for c in ordered_channels(s.details)]
            self._cards = [None] * len(s.channels)
            for i in order:
                card = _ChannelCard(s, i, dyes, probes, self._chan_host)
                card.changed.connect(self._on_card_changed)
                card.removeRequested.connect(self._remove_channel)
                card.colourRequested.connect(self._pick_colour)
                self._chan_lay.addWidget(card)
                self._cards[i] = card
        else:
            note = QLabel("A single picture without separate channels. Add a TIFF to turn it into a channel image.", self._chan_host)
            note.setWordWrap(True)
            self._chan_lay.addWidget(note)
        self._add_btn.setEnabled(True)
        self._auto_btn.setEnabled(s.is_multichannel)
        self._slide_btn.setEnabled(s.is_multichannel and s.slide_cell_count() > 1)
        self._slide_btn.setText(
            f"Apply to whole slide ({s.slide_cell_count()} cells)" if s.slide_cell_count() > 1 else "Apply to whole slide"
        )
        for w in (self._scale_show, self._pixel, self._length):
            w.blockSignals(True)
        self._scale_show.setChecked(s.scale_show)
        self._pixel.setText(f"{s.pixel_size_um:g}" if s.pixel_size_um else "")
        self._length.setCurrentText("Auto" if not s.scale_length_um else f"{s.scale_length_um:g} µm")
        for w in (self._scale_show, self._pixel, self._length):
            w.blockSignals(False)
        self._scale_section.setEnabled(s.is_multichannel)
        self._sync_log()
        pos, total, prev_id, next_id = self._neighbours
        self._nav.setText(f"cell {_cell_name(s.image, pos)} of {total}" if total else "")
        self._prev.setEnabled(bool(prev_id))
        self._next.setEnabled(bool(next_id))

    def _sync_log(self) -> None:
        s = self._session
        if s is None:
            return
        from repeatviewer.backends.cytogenetics.importer import refresh_log

        refresh_log(s.details)
        self._log_label.setText("\n".join(f"• {line}" for line in s.details.adjustment_log) or "No adjustments.")

    def _on_card_changed(self) -> None:
        self._sync_log()
        self._schedule()

    def _schedule(self) -> None:
        self._render_timer.start()
        self._save_timer.start()

    def _pick_colour(self, index: int) -> None:
        s = self._session
        if s is None:
            return
        colour = QColorDialog.getColor(QColor(s.channels[index].display.color), self, "Channel colour")
        if colour.isValid():
            s.set_color(index, colour.name())
            self._cards[index].paint_colour()
            self._schedule()

    def _auto_all(self) -> None:
        if self._session is None:
            return
        try:
            self._session.auto_levels()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Auto levels", str(exc))
            return
        for card in self._cards:
            card.refresh_levels()
        self._sync_log()
        self._schedule()

    def _add_channel(self) -> None:
        s = self._session
        if s is None:
            return
        path, _ = QFileDialog.getOpenFileName(self, "Add channel", str(self._folder or ""), "TIFF (*.tif *.tiff)")
        if not path:
            return
        try:
            s.add_channel(Path(path))
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Add channel", str(exc))
            return
        self._build_panel()
        self._schedule()

    def _remove_channel(self, index: int) -> None:
        s = self._session
        if s is None:
            return
        answer = QMessageBox.question(
            self, "Remove channel",
            f"Remove {s.label(index)} from this image?\nA copy in the project is deleted; the original microscope file is not touched.",
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        s.remove_channel(index)
        self._build_panel()
        self._render(keep_view=True)
        self._save_timer.start()

    def _on_pixel(self) -> None:
        s = self._session
        if s is None:
            return
        text = self._pixel.text().strip().replace(",", ".")
        try:
            s.set_pixel_size(float(text) if text else None)
        except ValueError:
            QMessageBox.information(self, "Pixel size", "Enter the pixel size in µm, for example 0.0645.")
            self._pixel.setText(f"{s.pixel_size_um:g}" if s.pixel_size_um else "")
            return
        self._on_scale()

    def _on_scale(self, *_args) -> None:
        s = self._session
        if s is None:
            return
        s.scale_show = self._scale_show.isChecked()
        s.scale_length_um = isn.parse_length(self._length.currentText())
        self._update_bar()
        self._save_timer.start()

    def _apply_to_slide(self) -> None:
        s = self._session
        if s is None:
            return
        n = s.slide_cell_count() - 1
        answer = QMessageBox.question(
            self, "Apply to whole slide",
            f"Give the other {n} cell{'s' if n != 1 else ''} of this slide the same colours, levels and gamma "
            "(matched by dye)? Raw files are not changed.",
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        self.flush()
        changed = s.apply_to_slide()
        self._notify()
        QMessageBox.information(self, "Apply to whole slide", f"{changed} other cell{'s' if changed != 1 else ''} updated.")

    # ------------------------------------------------------------------ view

    def _set_mode(self, merge: bool) -> None:
        self._merge_mode = merge
        _write("image_mode", "merge" if merge else "split")
        self._sync_mode_buttons()
        self._render()

    def _sync_mode_buttons(self) -> None:
        self._merge_btn.setChecked(self._merge_mode)
        self._split_btn.setChecked(not self._merge_mode)

    def _on_zoom(self, zoom: float) -> None:
        self._zoom_label.setText(f"{zoom * 100:.0f} %")

    def _go(self, which: int) -> None:
        target = self._neighbours[which]
        if target:
            self.show_image(target)

    def _render(self, keep_view: bool = True) -> None:
        s = self._session
        if s is None:
            return
        try:
            if not s.is_multichannel:
                image = QImage(str(ci.image_file_path(s.folder, s.image)))
                self._canvas.set_image(image if not image.isNull() else None, keep_view=keep_view,
                                       message="" if not image.isNull() else "The picture file is missing.")
                self._canvas.set_scale_bar(None)
                return
            rgb = s.merged() if self._merge_mode else s.side_by_side()
        except Exception as exc:  # noqa: BLE001 - missing TIFF, size mismatch
            self._canvas.set_image(None, keep_view=False, message=f"Cannot show this image:\n{exc}")
            self._canvas.set_scale_bar(None)
            return
        self._canvas.set_image(rgb_to_qimage(rgb), keep_view=keep_view)
        self._image_w = rgb.shape[1]
        self._update_bar()

    def _update_bar(self) -> None:
        s = self._session
        if s is None or not s.is_multichannel or not s.scale_show or not self._merge_mode or not self._image_w:
            self._canvas.set_scale_bar(None)
            return
        bar = s.scale_bar(self._image_w)
        self._canvas.set_scale_bar(None if bar is None else (bar[1], f"{bar[0]:g} µm"))

    # ------------------------------------------------------------------ saving

    def flush(self) -> None:
        """Write pending display changes (log, record, merged preview)."""
        self._save_timer.stop()
        s = self._session
        if s is None or self._folder is None:
            return
        try:
            s.save()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Image", f"Could not save the display settings:\n{exc}")
            return
        self._notify()

    def _notify(self) -> None:
        if self._folder is None:
            return
        from repeatviewer.widgets.cyto_images_ui import notify_changed

        self._saving = True
        try:
            notify_changed(self._folder)
        finally:
            self._saving = False

    def hideEvent(self, event) -> None:  # noqa: N802
        self.flush()
        super().hideEvent(event)
