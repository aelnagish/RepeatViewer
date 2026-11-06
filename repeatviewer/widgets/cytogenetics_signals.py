"""Signals tab of the Cytogenetics workspace: find chromosomes and signals, review them, pair homologues.

Step 1 Chromosomes: the counterstain is segmented (outlines, midlines, centromeres); Split, Merge,
Brush and Erase correct it. Step 2 Signals: every probe channel is searched for spots, drawn as
short ticks across the chromosome in the probe colour. Step 3 Pair: chromosomes are straightened,
paired into homologues (named after the karyotype) and summarised over the cells; a karyotype can be
made from the measured cells. Items the detection is unsure about are marked for checking.
"""

from __future__ import annotations

import copy
import math
import threading
from pathlib import Path

import numpy as np
from PySide6.QtCore import QObject, QPointF, QRectF, QSettings, Qt, QTimer, Signal
from PySide6.QtGui import QBrush, QColor, QFont, QImage, QKeySequence, QPainter, QPainterPath, QPen
from PySide6.QtWidgets import (
    QAbstractItemView,
    QButtonGroup,
    QCheckBox,
    QColorDialog,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QToolButton,
    QToolTip,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import channels as cyto_channels
from repeatviewer.backends.cytogenetics import lists as lists_mod
from repeatviewer.backends.cytogenetics import pairing as pr
from repeatviewer.backends.cytogenetics import segment as sm
from repeatviewer.backends.cytogenetics import signals as sg
from repeatviewer.backends.cytogenetics import signals_session as ss
from repeatviewer.backends.cytogenetics.signals_detect import PLANES, ProbeRow
from repeatviewer.backends.cytogenetics import spots as sp
from repeatviewer.backends.cytogenetics.karyotype import list_karyotypes
from repeatviewer.backends.cytogenetics.model import STRENGTHS, UNNAMED, ChromosomeMark, ImageDetails, Signal as CytoSignal, image_details, load_slides
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cytogenetics_image import ImageCanvas, _slider, rgb_to_qimage
from repeatviewer.widgets.cytogenetics_style import clear_layout
from repeatviewer.widgets.sequence_annotation_rail import _CollapsibleSection

_SETTINGS = "cytogenetics_workspace"
SIDE_W = 330
FIELD_W = 132
SAVE_DELAY_MS = 800
BLUE = QColor("#4C8DFF")
AMBER = QColor("#E8A33D")
LABEL_COLOUR = QColor("#FFFFFF")
NONE_KARYOTYPE = "__none__"
NEW_KARYOTYPE = "__new__"
EDIT_KARYOTYPE = "__edit__"
STEPS = (("chromosomes", "1 Chromosomes"), ("signals", "2 Signals"), ("pair", "3 Pair"))
TOOLS = (("select", "Select"), ("split", "Split"), ("merge", "Merge"), ("brush", "Brush"), ("erase", "Erase"),
         ("add", "Add signal"))
STEP_TOOLS = {"chromosomes": {"select", "split", "merge", "brush", "erase"}, "signals": {"select", "add"}}
SHOW = (("all", "Outlines and ticks"), ("ticks", "Ticks only"), ("image", "Image only"))
UNDO_STEPS = 60

HELP_DETECT = (
    "<b>Chromosomes from</b>: the counterstain channel the chromosomes are found in.<br>"
    "<b>Sensitivity</b>: higher finds fainter chromosomes; lower gives tighter outlines.<br>"
    "<b>Separate touching</b>: higher splits touching chromosomes harder (and may cut a chromosome at its centromere).<br>"
    "<b>Smallest</b>: objects smaller than this are left out (0 = automatic).<br>"
    "<b>Detect chromosomes</b> replaces the chromosomes of this cell, with their names. "
    "<b>Detect signals</b> searches every probe channel again; signals you added by hand stay."
)
HELP_PROBES = (
    "<b>Probes</b>: one row per channel of the cell. Untick a probe to hide its ticks.<br>"
    "<b>+</b> puts a probe from your lists on the first channel without one, or makes a <b>New Probe</b>. "
    "<b>−</b> takes the probe off the selected channel (its signals go too).<br>"
    "<b>Sensitivity</b>: higher finds weaker spots; the count n changes as you move it.<br>"
    "<b>Smallest spot</b>: diameter of the smallest spot. <b>Join spots within</b>: spots closer than this count as one.<br>"
    "<b>Same settings for the whole slide</b>: probes and settings go to every cell of the slide when saved."
)
HELP_TOOLS = (
    "<b>Select</b>: click a chromosome or signal; drag the white centromere tick along the midline, or a midline end "
    "to correct the length; drag a signal to move it. Backspace removes the selected one.<br>"
    "<b>Split</b>: draw a line across touching chromosomes to cut them apart.<br>"
    "<b>Merge</b>: click two pieces one after the other to join them.<br>"
    "<b>Brush</b> / <b>Erase</b>: paint onto or erase from a chromosome; painting on empty space adds a new one.<br>"
    "<b>Add signal</b>: click a spot to add a signal of the probe selected in the Probes table.<br>"
    "Right-click for more. Cmd+Z undoes. Mouse wheel zooms, Option + drag moves the picture."
)
HELP_PAIR = (
    "<b>Pairing</b>: chromosomes are straightened and paired by length, centromere index and signal pattern. "
    "With a karyotype they are named after it.<br>"
    "Drag a chromosome onto another to swap their names. Double-click one to open Edit karyotype.<br>"
    "<b>Across cells</b>: the table counts the cells with each signal on each chromosome.<br>"
    "<b>New from cells…</b> makes a karyotype from the mean lengths, centromere indices and the signals seen in at "
    "least the chosen share of cells; <b>Update</b> puts these values into the chosen karyotype. Both open Edit "
    "karyotype first, so you can check before saving."
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


def _qimage_rgb(image):
    """QImage -> numpy RGB array (for plain pictures)."""
    img = image.convertToFormat(QImage.Format.Format_RGB888)
    w, h, bpl = img.width(), img.height(), img.bytesPerLine()
    buf = np.frombuffer(img.constBits(), dtype=np.uint8, count=bpl * h).reshape(h, bpl)
    return buf[:, : w * 3].reshape(h, w, 3).copy()


def _rgba_image(rgba: np.ndarray) -> QImage:
    h, w = rgba.shape[:2]
    data = np.ascontiguousarray(rgba, dtype=np.uint8)
    img = QImage(data.data, w, h, 4 * w, QImage.Format.Format_RGBA8888)
    return img.copy()


def _swatch(colour: str, size: int = 12):
    from repeatviewer.views.cyto_new_probe_dialog import swatch_icon

    return swatch_icon(colour, size)


class _Relay(QObject):
    """Delivers a worker thread's result on the GUI thread (one per tab, never at module level)."""

    done = Signal(object, object)  # callback, result


def _segmented(parent: QWidget, items, on_click) -> tuple[QWidget, QButtonGroup, dict]:
    """A row of checkable buttons that touch, used for steps and tools (app buttons, no own style)."""
    host = QWidget(parent)
    row = QHBoxLayout(host)
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(0)
    group = QButtonGroup(host)
    group.setExclusive(True)
    buttons = {}
    for key, text in items:
        b = QPushButton(text, host)
        b.setCheckable(True)
        b.setSizePolicy(QSizePolicy.Policy.Minimum, QSizePolicy.Policy.Fixed)  # never squeezed below its text
        b.clicked.connect(lambda _c=False, k=key: on_click(k))
        group.addButton(b)
        row.addWidget(b)
        buttons[key] = b
    return host, group, buttons


# ---------------------------------------------------------------------------
# canvas
# ---------------------------------------------------------------------------


class SignalsCanvas(ImageCanvas):
    """Image with chromosome outlines, midlines, centromere ticks and signal ticks."""

    pressedAt = Signal(float, float, object)  # x, y, modifiers
    movedTo = Signal(float, float)
    releasedAt = Signal(float, float)
    hoverAt = Signal(float, float, object)  # x, y, global position
    deleteRequested = Signal()
    menuAt = Signal(float, float, object)
    doubleClickedAt = Signal(float, float)
    undoRequested = Signal()
    redoRequested = Signal()

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.session: ss.SignalsSession | None = None
        self.selected: tuple[str, str] | None = None
        self.picked: list[str] = []  # Merge: first piece
        self.overlay: QImage | None = None
        self.show_mode = "all"
        self.line: tuple[QPointF, QPointF] | None = None  # Split line in image pixels
        self.brush_at: tuple[float, float, float] | None = None  # x, y, radius in image pixels
        self.status = ""
        self.hidden: set[str] = set()
        self.colours: dict[str, str] = {}
        self._edge_cache: dict = {}
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
        self.setMouseTracking(True)
        self._panning = False

    # ---- events

    def event(self, ev) -> bool:  # noqa: D401
        from PySide6.QtCore import QEvent

        if ev.type() == QEvent.Type.ShortcutOverride:
            if ev.matches(QKeySequence.StandardKey.Undo) or ev.matches(QKeySequence.StandardKey.Redo):
                ev.accept()
                return True
        return super().event(ev)

    def mousePressEvent(self, event) -> None:  # noqa: N802
        self.setFocus()
        alt = bool(event.modifiers() & Qt.KeyboardModifier.AltModifier)
        if event.button() == Qt.MouseButton.MiddleButton or (event.button() == Qt.MouseButton.LeftButton and alt):
            self._panning = True
            self.start_pan(event.position())
        elif event.button() == Qt.MouseButton.LeftButton:
            p = self.image_point(event.position())
            self.pressedAt.emit(p.x(), p.y(), event.modifiers())
        elif event.button() == Qt.MouseButton.RightButton:
            p = self.image_point(event.position())
            self.menuAt.emit(p.x(), p.y(), event.globalPosition().toPoint())

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        p = self.image_point(event.position())
        if self._panning:
            self.move_pan(event.position())
        elif event.buttons() & Qt.MouseButton.LeftButton:
            self.movedTo.emit(p.x(), p.y())
        else:
            self.hoverAt.emit(p.x(), p.y(), event.globalPosition().toPoint())

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if self._panning:
            self._panning = False
            self.end_pan()
        elif event.button() == Qt.MouseButton.LeftButton:
            p = self.image_point(event.position())
            self.releasedAt.emit(p.x(), p.y())

    def end_pan(self) -> None:
        cursor = self.cursor()
        super().end_pan()
        self.setCursor(cursor if cursor.shape() != Qt.CursorShape.ClosedHandCursor else self.tool_cursor)

    tool_cursor = Qt.CursorShape.ArrowCursor

    def center_on(self, x: float, y: float) -> None:
        """Move the view so an image point is in the middle (zooming in when the whole picture is shown)."""
        if self._image is None:
            return
        if self._fit:
            self._zoom = min(self._zoom * 2.0, 8.0)
            self._fit = False
        self._off = QPointF(self.width() / 2 - x * self._zoom, self.height() / 2 - y * self._zoom)
        self.update()

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        p = self.image_point(event.position())
        self.doubleClickedAt.emit(p.x(), p.y())

    def keyPressEvent(self, event) -> None:  # noqa: N802
        if event.matches(QKeySequence.StandardKey.Undo):
            self.undoRequested.emit()
        elif event.matches(QKeySequence.StandardKey.Redo):
            self.redoRequested.emit()
        elif event.key() in (Qt.Key.Key_Delete, Qt.Key.Key_Backspace):
            self.deleteRequested.emit()
        else:
            super().keyPressEvent(event)

    def leaveEvent(self, event) -> None:  # noqa: N802
        if self.brush_at is not None:
            self.brush_at = None
            self.update()
        super().leaveEvent(event)

    # ---- drawing

    def _image_rect(self) -> QRectF | None:
        if self._image is None or self._image.isNull():
            return None
        o = self._origin()
        return QRectF(o.x(), o.y(), self._image.width() * self._zoom, self._image.height() * self._zoom)

    def _tick(self, p: QPainter, line, at: tuple[float, float], half: float, colour: QColor, width: float,
              dashed: bool = False) -> None:
        """A short line across the chromosome at the midline point nearest to ``at``."""
        pts = np.asarray(line, dtype=np.float64)
        k = int(np.argmin(np.hypot(pts[:, 0] - at[0], pts[:, 1] - at[1])))
        a, b = pts[max(0, k - 2)], pts[min(len(pts) - 1, k + 2)]
        d = b - a
        n = math.hypot(*d) or 1.0
        nx, ny = -d[1] / n, d[0] / n
        c = self.screen_point(float(pts[k, 0]), float(pts[k, 1]))
        pen = QPen(colour, width)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        if dashed:
            pen.setStyle(Qt.PenStyle.DashLine)
        p.setPen(pen)
        p.drawLine(QPointF(c.x() - nx * half, c.y() - ny * half), QPointF(c.x() + nx * half, c.y() + ny * half))

    def paintEvent(self, event) -> None:  # noqa: N802
        super().paintEvent(event)
        s = self.session
        if s is None:
            return
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        rect = self._image_rect()
        if rect is not None and self.overlay is not None and self.show_mode == "all":
            p.drawImage(rect, self.overlay)
        if rect is not None and self.show_mode == "all":
            self._draw_outlines(p)
        z = self._zoom
        font = QFont(self.font())
        font.setBold(True)
        p.setFont(font)
        marks = {m.id: m for m in s.details.chromosomes}
        lines = {m.id: _line_of(m) for m in s.details.chromosomes}
        halves = {m.id: min(16.0, max(6.0, 0.9 * pr.half_width(m) * z + 2.0)) for m in s.details.chromosomes}
        if self.show_mode != "image":
            for m in s.details.chromosomes:
                line = lines[m.id]
                if not line or len(line) < 2:
                    continue
                if self.show_mode == "all":
                    pen = QPen(QColor(255, 255, 255, 130), 1, Qt.PenStyle.DashLine)
                    p.setPen(pen)
                    path = QPainterPath(self.screen_point(*line[0]))
                    for q in line[1:]:
                        path.lineTo(self.screen_point(*q))
                    p.drawPath(path)
                if m.centromere is not None:
                    self._tick(p, line, m.centromere, halves[m.id], QColor("#FFFFFF"), 2.5)
            for sig in s.details.signals:
                if sig.probe in self.hidden:
                    continue
                colour = QColor(self.colours.get(sig.probe, "#FFFFFF"))
                chosen = self.selected == ("signal", sig.id)
                mark = marks.get(sig.chromosome_id)
                line = lines.get(sig.chromosome_id)
                if mark is not None and line and len(line) >= 2:
                    if chosen:
                        self._tick(p, line, (sig.x, sig.y), halves[mark.id] + 3, BLUE, 6.5)
                    self._tick(p, line, (sig.x, sig.y), halves[mark.id] + 1, colour, 3.0, dashed=bool(sig.check))
                else:
                    c = self.screen_point(sig.x, sig.y)
                    p.setPen(QPen(BLUE if chosen else colour, 2.5))
                    r = 6.0
                    p.drawLine(QPointF(c.x() - r, c.y() - r), QPointF(c.x() + r, c.y() + r))
                    p.drawLine(QPointF(c.x() - r, c.y() + r), QPointF(c.x() + r, c.y() - r))
            # labels and "?" badges
            for m in s.details.chromosomes:
                anchor = self.screen_point(m.x, m.y)
                label = m.label if m.named else m.homologue
                box = QRectF(anchor.x() + halves.get(m.id, 10) + 2, anchor.y() - 10, 60, 20)
                p.setPen(QColor(0, 0, 0, 200))
                p.drawText(box.translated(1, 1), Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, label)
                p.setPen(LABEL_COLOUR)
                p.drawText(box, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, label)
                if m.check:
                    w = p.fontMetrics().horizontalAdvance(label)
                    self._badge(p, QPointF(box.x() + w + 11, anchor.y()))
            for sig in s.details.signals:
                if sig.check and sig.probe not in self.hidden:
                    c = self.screen_point(sig.x, sig.y)
                    self._badge(p, QPointF(c.x() + 12, c.y() - 12))
            # selected chromosome: handles on the midline ends
            if self.selected and self.selected[0] == "mark":
                m = marks.get(self.selected[1])
                if m is not None and lines.get(m.id):
                    for end in (lines[m.id][0], lines[m.id][-1]):
                        q = self.screen_point(*end)
                        p.setPen(QPen(QColor("#000000"), 1))
                        p.setBrush(QColor("#FFFFFF"))
                        p.drawRect(QRectF(q.x() - 4, q.y() - 4, 8, 8))
                    p.setBrush(Qt.BrushStyle.NoBrush)
        if self.line is not None:
            a, b = self.screen_point(self.line[0].x(), self.line[0].y()), self.screen_point(self.line[1].x(), self.line[1].y())
            pen = QPen(BLUE, 2, Qt.PenStyle.DashLine)
            p.setPen(pen)
            p.drawLine(a, b)
        if self.brush_at is not None:
            x, y, r = self.brush_at
            c = self.screen_point(x, y)
            p.setPen(QPen(QColor("#FFFFFF"), 1.5))
            p.drawEllipse(c, r * z, r * z)
        if self.status:
            fm = p.fontMetrics()
            p.setFont(self.font())
            fm = p.fontMetrics()
            w = fm.horizontalAdvance(self.status) + 20
            box = QRectF(14, self.height() - 14 - 26, w, 26)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(255, 255, 255, 235))
            p.drawRoundedRect(box, 8, 8)
            p.setPen(QColor("#1D1D1F"))
            p.drawText(box, Qt.AlignmentFlag.AlignCenter, self.status)
        p.end()

    def _edges(self, m: ChromosomeMark) -> list:
        key = (m.id, id(m.mask))
        hit = self._edge_cache.get(m.id)
        if hit is not None and hit[0] == key:
            return hit[1]
        dec = sm.decode_mask(m.mask) if m.mask else None
        lines = mask_edges(*dec) if dec is not None else []
        self._edge_cache[m.id] = (key, lines)
        return lines

    def _draw_outlines(self, p: QPainter) -> None:
        s = self.session
        o = self._origin()
        p.save()
        p.translate(o)
        p.scale(self._zoom, self._zoom)
        selected = self.selected[1] if self.selected and self.selected[0] == "mark" else ""
        for m in s.details.chromosomes:
            lines = self._edges(m)
            if not lines:
                continue
            chosen = m.id == selected or m.id in self.picked
            shadow = QPen(QColor(0, 0, 0, 150), 3.5)
            shadow.setCosmetic(True)
            p.setPen(shadow)
            p.drawLines(lines)
            colour = BLUE if chosen else (AMBER if m.check else QColor(255, 255, 255, 230))
            pen = QPen(colour, 2.5 if chosen else 1.4)
            pen.setCosmetic(True)
            if m.check and not chosen:
                pen.setStyle(Qt.PenStyle.DashLine)
            p.setPen(pen)
            p.drawLines(lines)
        p.restore()

    def _badge(self, p: QPainter, centre: QPointF) -> None:
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(AMBER)
        p.drawEllipse(centre, 8.0, 8.0)
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.setPen(QColor("#1D1D1F"))
        p.drawText(QRectF(centre.x() - 8, centre.y() - 8, 16, 16), Qt.AlignmentFlag.AlignCenter, "?")


def _line_of(m: ChromosomeMark) -> list | None:
    """Midline of a detected chromosome, or telomere-centromere-telomere of one marked by hand."""
    if m.midline and len(m.midline) >= 2:
        return m.midline
    pts = [q for q in (m.telomere_p, m.centromere, m.telomere_q) if q is not None]
    if len(pts) < 2:
        return None
    out = []
    for a, b in zip(pts, pts[1:]):  # densify so ticks find a nearby point
        for t in np.linspace(0, 1, 12, endpoint=False):
            out.append((a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t))
    out.append(pts[-1])
    return out


def mask_edges(mask: np.ndarray, x0: int, y0: int) -> list:
    """Outline of a mask as horizontal and vertical QLineF runs in image pixels (drawn with a cosmetic pen,
    so the outline is equally thin and visible at every zoom, also on large images)."""
    from PySide6.QtCore import QLineF

    m = np.pad(np.asarray(mask, dtype=bool), 1)
    out = []
    # horizontal edges: between row r-1 and r
    diff = m[1:, :] != m[:-1, :]
    for r in np.flatnonzero(diff.any(axis=1)):
        row = diff[r].astype(np.int8)
        d = np.diff(np.concatenate(([0], row, [0])))
        for a, b in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)):
            y = y0 + r - 0.0
            out.append(QLineF(x0 + a - 1, y, x0 + b - 1, y))
    diff = m[:, 1:] != m[:, :-1]
    for c in np.flatnonzero(diff.any(axis=0)):
        col = diff[:, c].astype(np.int8)
        d = np.diff(np.concatenate(([0], col, [0])))
        for a, b in zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1)):
            x = x0 + c
            out.append(QLineF(x, y0 + a - 1, x, y0 + b - 1))
    return out


def build_overlay(session: ss.SignalsSession, selected: str = "", picked: list[str] | None = None) -> QImage | None:
    """Pale fill of every chromosome (outlines are drawn as lines by the canvas)."""
    labels = session.labels_map()
    if labels is None or not labels.any():
        return None
    h, w = labels.shape
    rgba = np.zeros((h, w, 4), dtype=np.uint8)
    inside = labels > 0
    rgba[inside] = (90, 140, 255, 46)
    return _rgba_image(rgba)


# ---------------------------------------------------------------------------
# step 3: straightened chromosomes and probe profiles
# ---------------------------------------------------------------------------


class StripsView(QWidget):
    """Karyogram of straightened chromosomes, homologues side by side, p telomere at the top."""

    selected = Signal(str)  # mark id
    swapRequested = Signal(str, str)
    openRequested = Signal(str)

    STRIP_W = 20

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumHeight(260)
        self.groups: list[tuple[str, list[dict]]] = []  # (chromosome name, [{id, image, length, cen, bands, label}])
        self.current = ""
        self.checks: dict[str, str] = {}
        self._rects: list[tuple[QRectF, str]] = []
        self._drag = ""
        self._drag_pos: QPointF | None = None

    def set_groups(self, groups, checks: dict[str, str]) -> None:
        self.groups = groups
        self.checks = checks
        self.update()

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.fillRect(self.rect(), QColor("#0B0D12"))
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        self._rects = []
        if not self.groups:
            p.setPen(QColor("#B9C3D6"))
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Find the chromosomes first (step 1).")
            return
        p.setPen(QColor("#B9C3D6"))
        p.drawText(QRectF(16, 6, self.width() - 32, 18), Qt.AlignmentFlag.AlignLeft,
                   "p telomere at the top. Drag a chromosome onto another to swap their names.")
        longest = max((it["length"] for _n, items in self.groups for it in items), default=1.0) or 1.0
        top, bottom_pad = 34.0, 66.0
        max_h = max(60.0, self.height() - top - bottom_pad)
        slot = max(40.0, (self.width() - 20) / max(1, len(self.groups)))
        font = QFont(self.font())
        for gi, (name, items) in enumerate(self.groups):
            n = len(items)
            gap = 8.0
            total = n * self.STRIP_W + (n - 1) * gap
            x0 = 10 + slot * (gi + 0.5) - total / 2
            lowest = top
            for j, it in enumerate(items):
                h = max(20.0, max_h * it["length"] / longest)
                r = QRectF(x0 + j * (self.STRIP_W + gap), top, self.STRIP_W, h)
                lowest = max(lowest, r.bottom())
                self._rects.append((r, it["id"]))
                path = QPainterPath()
                path.addRoundedRect(r, self.STRIP_W / 2, self.STRIP_W / 2)
                p.save()
                p.setClipPath(path)
                if it.get("image") is not None:
                    p.drawImage(r, it["image"])
                else:
                    p.fillRect(r, QColor("#4C7BFF"))
                for frac, colour in it["bands"]:
                    p.fillRect(QRectF(r.x(), r.y() + frac * h - 3, r.width(), 6), QColor(colour))
                p.restore()
                if it.get("cen") is not None:
                    p.fillRect(QRectF(r.x() - 2, r.y() + it["cen"] * h - 1.5, r.width() + 4, 3), QColor("#FFFFFF"))
                if it["id"] == self.current:
                    p.setPen(QPen(BLUE, 2))
                    p.drawRoundedRect(r.adjusted(-3, -3, 3, 3), 12, 12)
                p.setPen(QColor("#B9C3D6"))
                font.setBold(False)
                p.setFont(font)
                p.drawText(QRectF(r.x() - 10, r.bottom() + 2, r.width() + 20, 16), Qt.AlignmentFlag.AlignHCenter, it["label"])
            font.setBold(True)
            p.setFont(font)
            p.setPen(QColor("#FFFFFF"))
            p.drawText(QRectF(x0 - 20, lowest + 18, total + 40, 18), Qt.AlignmentFlag.AlignHCenter,
                       name if name != UNNAMED else "not paired")
            if name in self.checks:
                badge = QRectF(x0 + total / 2 - 26, lowest + 38, 52, 16)
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(AMBER)
                p.drawRoundedRect(badge, 8, 8)
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.setPen(QColor("#1D1D1F"))
                f2 = QFont(font)
                f2.setPointSizeF(max(7.0, font.pointSizeF() - 2))
                p.setFont(f2)
                p.drawText(badge, Qt.AlignmentFlag.AlignCenter, "check")
        if self._drag and self._drag_pos is not None:
            p.setPen(QPen(BLUE, 2, Qt.PenStyle.DashLine))
            p.drawEllipse(self._drag_pos, 10, 10)
        p.end()

    def _hit(self, pos: QPointF) -> str:
        for r, ident in self._rects:
            if r.adjusted(-4, -4, 4, 4).contains(pos):
                return ident
        return ""

    def mousePressEvent(self, event) -> None:  # noqa: N802
        ident = self._hit(event.position())
        if ident:
            self.current = ident
            self._drag = ident
            self.selected.emit(ident)
            self.update()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if self._drag:
            self._drag_pos = event.position()
            self.update()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        target = self._hit(event.position())
        if self._drag and target and target != self._drag:
            self.swapRequested.emit(self._drag, target)
        self._drag = ""
        self._drag_pos = None
        self.update()

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        ident = self._hit(event.position())
        if ident:
            self.openRequested.emit(ident)


class ProfileView(QWidget):
    """One track per probe along the selected chromosome: homologue a solid, b dashed; grey band = spread of the
    signal over the cells in scope."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setMinimumHeight(220)
        self.title = ""
        self.cen: float | None = None
        self.tracks: list[dict] = []  # {name, colour, lines: [(values, dashed)], spread: (lo, hi) | None, note}

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.fillRect(self.rect(), self.palette().base())
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        if not self.tracks:
            p.setPen(self.palette().placeholderText().color())
            p.drawText(self.rect(), Qt.AlignmentFlag.AlignCenter, "Click a chromosome above to see its probe tracks.")
            return
        left, right = 110.0, 20.0
        x0, x1 = left, self.width() - right
        row_h = max(22.0, min(40.0, (self.height() - 40) / max(1, len(self.tracks))))
        text = self.palette().text().color()
        muted = self.palette().placeholderText().color()
        if self.cen is not None:
            xc = x0 + (x1 - x0) * self.cen
            p.setPen(QPen(muted, 1, Qt.PenStyle.DashLine))
            p.drawLine(QPointF(xc, 4), QPointF(xc, 10 + row_h * len(self.tracks)))
        for i, tr in enumerate(self.tracks):
            base = 10 + row_h * (i + 1) - 4
            p.setPen(text)
            p.drawText(QRectF(0, base - row_h + 4, left - 10, row_h), Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter, tr["name"])
            if tr.get("spread"):
                lo, hi = tr["spread"]
                p.fillRect(QRectF(x0 + (x1 - x0) * lo, base - row_h + 8, max(3.0, (x1 - x0) * (hi - lo)), row_h - 8),
                           QColor(140, 140, 150, 60))
            for values, dashed in tr["lines"]:
                if values is None or len(values) < 2:
                    continue
                v = np.asarray(values, dtype=np.float64)
                top = float(v.max()) or 1.0
                path = QPainterPath()
                for k, val in enumerate(v):
                    pt = QPointF(x0 + (x1 - x0) * k / (len(v) - 1), base - (row_h - 10) * val / top)
                    path.moveTo(pt) if k == 0 else path.lineTo(pt)
                pen = QPen(QColor(tr["colour"]), 2)
                if dashed:
                    pen.setStyle(Qt.PenStyle.DashLine)
                p.setPen(pen)
                p.drawPath(path)
            if tr.get("note"):
                p.setPen(muted)
                p.drawText(QRectF(x0, base - row_h + 2, x1 - x0, 14), Qt.AlignmentFlag.AlignRight, tr["note"])
        p.setPen(muted)
        y = 10 + row_h * len(self.tracks) + 14
        p.drawText(QPointF(x0, y), "p-ter")
        p.drawText(QPointF(x1 - 30, y), "q-ter")
        if self.cen is not None:
            p.drawText(QPointF(x0 + (x1 - x0) * self.cen - 10, y), "cen")
        p.end()


# ---------------------------------------------------------------------------
# the tab
# ---------------------------------------------------------------------------


FORM_LABELS = ("Chromosomes from", "Separate touching", "Join spots within", "Centromere index", "Seen in at least",
               "Smallest spot", "Sensitivity")


def _form_row(parent: QWidget, label: str, widget: QWidget) -> QHBoxLayout:
    """Label of one common width, then the field taking the rest of the panel: every field of the side panel
    starts at the same place and has the same width, whatever the font."""
    row = QHBoxLayout()
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(8)
    lab = QLabel(label, parent)
    fm = lab.fontMetrics()
    lab.setFixedWidth(max(fm.horizontalAdvance(t) for t in FORM_LABELS + (label,)) + 6)
    row.addWidget(lab)
    widget.setMinimumWidth(60)
    widget.setMaximumWidth(16777215)
    widget.setSizePolicy(QSizePolicy.Policy.Expanding, widget.sizePolicy().verticalPolicy())
    row.addWidget(widget, 1)
    return row


def _number(text: str, default: float) -> float:
    """First number in a field such as "2 px"."""
    import re

    m = re.search(r"\d+(?:[.,]\d+)?", text or "")
    return float(m.group(0).replace(",", ".")) if m else default


def _muted(text: str, parent: QWidget) -> QLabel:
    lab = QLabel(text, parent)
    lab.setWordWrap(True)
    lab.setStyleSheet("color: gray;")
    lab.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
    return lab


class SignalsTab(QWidget):
    name = "Signals"
    imageChanged = Signal(str)
    idiogramRequested = Signal()
    karyotypesChanged = Signal()
    figureRequested = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        from repeatviewer.widgets.help_button import make_help_button

        self._folder: Path | None = None
        self._session: ss.SignalsSession | None = None
        self._step = _read("signals_step", "chromosomes")
        if self._step not in dict(STEPS):
            self._step = "chromosomes"
        self._tool = "select"
        self._selected: tuple[str, str] | None = None
        self._drag: tuple[str, str] | None = None
        self._saving = False
        self._pending = ""
        self._busy = False
        self._undo: list[dict] = []
        self._redo: list[dict] = []
        self._probe_key = ""
        self._overlay_key = None
        self._relay = _Relay(self)
        self._relay.done.connect(lambda cb, result: cb(result))
        self._pairing_checks: dict[str, str] = {}
        self._strip_cache: dict[str, QImage] = {}

        # ---- toolbar
        # two rows, so nothing is squeezed when the Files panel takes room: karyotype, steps, zoom and cells
        # on top; the tools of the step (or the step 3 view) below
        bar = QHBoxLayout()
        bar.setContentsMargins(10, 8, 10, 4)
        bar.setSpacing(8)
        bar2 = QHBoxLayout()
        bar2.setContentsMargins(10, 0, 10, 6)
        bar2.setSpacing(8)
        self._karyo = NoWheelComboBox(self)
        self._karyo.setMinimumWidth(130)
        self._karyo.setMaximumWidth(200)
        self._karyo.setSizeAdjustPolicy(NoWheelComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self._karyo.setMinimumContentsLength(14)
        self._karyo.setToolTip("Karyotype of this cell: names the chromosomes when they are paired")
        self._karyo.activated.connect(self._on_karyotype)
        bar.addWidget(self._karyo)
        steps_host, _g, self._step_buttons = _segmented(self, STEPS, self._set_step)
        bar.addWidget(steps_host)
        self._sep = QFrame(self)
        self._sep.setFrameShape(QFrame.Shape.VLine)
        self._sep.setFrameShadow(QFrame.Shadow.Sunken)
        bar2.addWidget(self._sep)
        self._tools_host, _g2, self._tool_buttons = _segmented(self, TOOLS, self._set_tool)
        bar2.addWidget(self._tools_host)
        self._show = NoWheelComboBox(self)
        for key, text in SHOW:
            self._show.addItem(text, key)
        self._show.setMinimumWidth(120)
        self._show.setCurrentIndex(max(0, self._show.findData(_read("signals_show", "all"))))
        self._show.currentIndexChanged.connect(self._on_show)
        bar2.addWidget(self._show)
        # step 3 controls
        self._view_host, _g3, self._view_buttons = _segmented(self, (("image", "Image"), ("strips", "Straightened")), self._set_view)
        bar2.addWidget(self._view_host)
        self._sort_label = QLabel("Sort by", self)
        bar2.addWidget(self._sort_label)
        self._sort = NoWheelComboBox(self)
        for key, text in (("karyotype", "Karyotype order"), ("length", "Length")):
            self._sort.addItem(text, key)
        self._sort.setMinimumWidth(120)
        self._sort.currentIndexChanged.connect(lambda _i: self._sync_pair())
        bar2.addWidget(self._sort)
        bar2.addStretch(1)
        bar.addStretch(1)
        zoom_row = bar2  # zoom sits with the tools, the cell navigation with the steps
        for text, tip, slot in (("−", "Zoom out", lambda: self._canvas.set_zoom(self._canvas.zoom() / 1.25)),
                                ("+", "Zoom in", lambda: self._canvas.set_zoom(self._canvas.zoom() * 1.25)),
                                ("Fit", "Fit the picture", lambda: self._canvas.fit())):
            zb = QToolButton(self)  # small buttons: a push button has a wide minimum on macOS
            zb.setText(text)
            zb.setToolTip(tip)
            zb.clicked.connect(slot)
            zoom_row.addWidget(zb)
        self._prev = QToolButton(self)
        self._next = QToolButton(self)
        for b, text, tip in ((self._prev, "‹", "Previous cell of this slide"), (self._next, "›", "Next cell of this slide")):
            b.setText(text)
            b.setToolTip(tip)
        self._prev.clicked.connect(lambda: self._step_cell(-1))
        self._next.clicked.connect(lambda: self._step_cell(1))
        self._nav = QLabel(self)
        bar.addWidget(self._prev)
        bar.addWidget(self._nav)
        bar.addWidget(self._next)
        bar_host = QWidget(self)
        bars = QVBoxLayout(bar_host)
        bars.setContentsMargins(0, 0, 0, 0)
        bars.setSpacing(0)
        bars.addLayout(bar)
        bars.addLayout(bar2)

        # ---- canvas and step 3 views
        self._canvas = SignalsCanvas(self)
        self._canvas.pressedAt.connect(self._on_press)
        self._canvas.movedTo.connect(self._on_move)
        self._canvas.releasedAt.connect(self._on_release)
        self._canvas.hoverAt.connect(self._on_hover)
        self._canvas.deleteRequested.connect(self._delete_selected)
        self._canvas.menuAt.connect(self._on_menu)
        self._canvas.doubleClickedAt.connect(self._on_double)
        self._canvas.undoRequested.connect(self.undo)
        self._canvas.redoRequested.connect(self.redo)
        self._strips = StripsView(self)
        self._strips.selected.connect(self._on_strip_selected)
        self._strips.swapRequested.connect(self._swap_names)
        self._strips.openRequested.connect(lambda _i: self._open_karyotype(edit=True))
        self._profile = ProfileView(self)
        self._profile_title = QLabel(self)
        self._profile_avg = QCheckBox("Average over cells", self)
        self._profile_avg.setChecked(True)
        self._profile_avg.toggled.connect(lambda _on: self._sync_profile())
        strips_page = QWidget(self)
        sl = QVBoxLayout(strips_page)
        sl.setContentsMargins(0, 0, 0, 0)
        sl.setSpacing(0)
        sl.addWidget(self._strips, 3)
        prof_head = QHBoxLayout()
        prof_head.setContentsMargins(16, 8, 16, 0)
        prof_head.addWidget(self._profile_title, 1)
        prof_head.addWidget(self._profile_avg)
        sl.addLayout(prof_head)
        sl.addWidget(self._profile, 2)
        self._views = QStackedWidget(self)
        self._views.addWidget(self._canvas)
        self._views.addWidget(strips_page)

        left = QVBoxLayout()
        left.setContentsMargins(0, 0, 0, 0)
        left.setSpacing(0)
        left.addWidget(bar_host)
        left.addWidget(self._views, 1)

        # ---- side panel (steps 1 and 2)
        self._side_review = self._build_review_side(make_help_button)
        self._side_pair = self._build_pair_side(make_help_button)
        self._sides = QStackedWidget(self)
        for body in (self._side_review, self._side_pair):
            scroll = QScrollArea(self)
            scroll.setWidgetResizable(True)
            scroll.setFrameShape(QFrame.Shape.NoFrame)
            scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            scroll.setWidget(body)
            from repeatviewer.widgets.cytogenetics_style import fit_side_scroll

            fit_side_scroll(scroll)
            self._sides.addWidget(scroll)
        self._sides.setMinimumWidth(280)
        self._sides.setMaximumWidth(SIDE_W)

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        left_host = QWidget(self)
        left_host.setLayout(left)
        root.addWidget(left_host, 1)
        line = QFrame(self)
        line.setFrameShape(QFrame.Shape.VLine)
        line.setFrameShadow(QFrame.Shadow.Plain)
        line.setStyleSheet("color: palette(mid);")
        root.addWidget(line)
        root.addWidget(self._sides)

        self._save_timer = QTimer(self)
        self._save_timer.setSingleShot(True)
        self._save_timer.setInterval(SAVE_DELAY_MS)
        self._save_timer.timeout.connect(self.flush)
        self._spot_timer = QTimer(self)
        self._spot_timer.setSingleShot(True)
        self._spot_timer.setInterval(300)
        self._spot_timer.timeout.connect(self._detect_selected_probe)
        self._set_view(_read("signals_view", "strips"))
        self._set_step(self._step)
        self._show_empty()

    # ------------------------------------------------------------------ side panels

    def _build_review_side(self, make_help_button) -> QWidget:
        body = QWidget(self)
        side = QVBoxLayout(body)
        side.setContentsMargins(12, 10, 12, 10)
        side.setSpacing(12)

        # Detect
        sec = _CollapsibleSection("Detect", body, expanded=True)
        self._det_from = NoWheelComboBox(self)
        self._det_from.setToolTip("Channel the chromosomes are found in")
        self._det_from.activated.connect(self._on_det_from)
        self._det_sens = _slider(0, 100, 50)
        self._det_sep = _slider(0, 100, 50)
        self._det_small = QLineEdit(self)
        self._det_small.setPlaceholderText("auto")
        self._det_small.setToolTip("Smallest chromosome area in pixels; empty = automatic")
        for w in (self._det_sens, self._det_sep):
            w.sliderReleased.connect(self._remember_detect)
        sec.add_layout(_form_row(body, "Chromosomes from", self._det_from))
        sec.add_layout(_form_row(body, "Sensitivity", self._det_sens))
        sec.add_layout(_form_row(body, "Separate touching", self._det_sep))
        sec.add_layout(_form_row(body, "Smallest", self._det_small))
        brow = QHBoxLayout()
        self._btn_detect = QPushButton("Detect Chromosomes", body)
        self._btn_detect.clicked.connect(self._detect_chromosomes)
        self._btn_signals = QPushButton("Detect Signals", body)
        self._btn_signals.clicked.connect(self._detect_all_signals)
        brow.addWidget(self._btn_detect, 1)
        brow.addWidget(self._btn_signals, 1)
        brow.addWidget(make_help_button(body, "Detect", HELP_DETECT))
        sec.add_layout(brow)
        side.addWidget(sec)

        # Probes
        sec = _CollapsibleSection("Probes", body, expanded=True)
        self._probes = QTableWidget(0, 4, body)
        self._probes.setHorizontalHeaderLabels(["", "Probe", "Channel", "n"])
        self._probes.verticalHeader().setVisible(False)
        self._probes.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._probes.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self._probes.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._probes.setShowGrid(False)
        self._probes.setWordWrap(False)
        self._probes.setTextElideMode(Qt.TextElideMode.ElideMiddle)
        self._probes.installEventFilter(self)
        hh = self._probes.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)  # columns stay resizable
        hh.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        hh.setStretchLastSection(False)
        hh.sectionResized.connect(self._probe_col_dragged)
        self._probes.itemChanged.connect(self._on_probe_item)
        self._probes.itemSelectionChanged.connect(self._on_probe_selected)
        sec.add_widget(self._probes)
        prow = QHBoxLayout()
        prow.setSpacing(0)
        self._probe_add = QPushButton("+", body)
        self._probe_add.setToolTip("Add a probe")
        self._probe_add.setFixedWidth(28)
        self._probe_add.clicked.connect(self._probe_add_menu)
        self._probe_remove = QPushButton("−", body)
        self._probe_remove.setToolTip("Take the probe off the selected channel")
        self._probe_remove.setFixedWidth(28)
        self._probe_remove.clicked.connect(self._remove_probe)
        prow.addWidget(self._probe_add)
        prow.addWidget(self._probe_remove)
        prow.addStretch(1)
        self._btn_detect_all = QPushButton("Detect All", body)
        self._btn_detect_all.clicked.connect(self._detect_all_signals)
        prow.addWidget(self._btn_detect_all)
        prow.addSpacing(6)
        prow.addWidget(make_help_button(body, "Probes", HELP_PROBES))
        sec.add_layout(prow)
        side.addWidget(sec)

        # selected probe
        self._probe_sec = _CollapsibleSection("Probe", body, expanded=True)
        self._p_channel = NoWheelComboBox(self)
        self._p_channel.activated.connect(self._on_probe_channel)
        self._p_dye = NoWheelComboBox(self)
        self._p_dye.activated.connect(self._on_probe_dye)
        self._p_colour = QPushButton(self)
        self._p_colour.clicked.connect(self._on_probe_colour)
        sens_host = QWidget(self)
        sh = QHBoxLayout(sens_host)
        sh.setContentsMargins(0, 0, 0, 0)
        self._p_sens = _slider(0, 100, 50)
        self._p_sens_value = QLabel("50", sens_host)
        self._p_sens_value.setFixedWidth(26)
        sh.addWidget(self._p_sens, 1)
        sh.addWidget(self._p_sens_value)
        self._p_sens.valueChanged.connect(self._on_probe_settings)
        self._p_small = QLineEdit(self)
        self._p_join = QLineEdit(self)
        for w in (self._p_small, self._p_join):
            w.editingFinished.connect(self._on_probe_settings)
        self._p_slide = QCheckBox("Same settings for the whole slide", self)
        self._p_slide.setChecked(_read("signals_same_slide", "true") != "false")
        self._p_slide.toggled.connect(self._on_same_slide)
        self._probe_sec.add_layout(_form_row(body, "Channel", self._p_channel))
        self._probe_sec.add_layout(_form_row(body, "Dye", self._p_dye))
        self._probe_sec.add_layout(_form_row(body, "Colour", self._p_colour))
        self._probe_sec.add_layout(_form_row(body, "Sensitivity", sens_host))
        self._probe_sec.add_layout(_form_row(body, "Smallest spot", self._p_small))
        self._probe_sec.add_layout(_form_row(body, "Join spots within", self._p_join))
        self._probe_sec.add_widget(self._p_slide)
        self._probe_sec.add_widget(_muted("The count updates as you move the slider. Untick a probe to hide its ticks.", body))
        side.addWidget(self._probe_sec)

        # Check
        sec = _CollapsibleSection("Check", body, expanded=True)
        self._check_text = _muted("", body)
        self._check_text.setStyleSheet("")
        sec.add_widget(self._check_text)
        crow = QHBoxLayout()
        self._check_prev = QPushButton("‹ Previous", body)
        self._check_next = QPushButton("Next to Check ›", body)
        self._check_prev.clicked.connect(lambda: self._step_check(-1))
        self._check_next.clicked.connect(lambda: self._step_check(1))
        crow.addWidget(self._check_prev, 1)
        crow.addWidget(self._check_next, 1)
        sec.add_layout(crow)
        arow = QHBoxLayout()
        self._accept_one = QPushButton("Accept", body)
        self._accept_one.clicked.connect(self._accept_selected)
        self._accept_all = QPushButton("Accept All in This Cell", body)
        self._accept_all.clicked.connect(self._accept_everything)
        arow.addWidget(self._accept_one, 1)
        arow.addWidget(self._accept_all, 1)
        arow.addWidget(make_help_button(body, "Tools", HELP_TOOLS))
        sec.add_layout(arow)
        side.addWidget(sec)

        # Selected
        self._sel_section = _CollapsibleSection("Selected", body, expanded=True)
        self._sel_host = QWidget(self._sel_section)
        self._sel_lay = QVBoxLayout(self._sel_host)
        self._sel_lay.setContentsMargins(0, 0, 0, 0)
        self._sel_lay.setSpacing(6)
        self._sel_section.add_widget(self._sel_host)
        side.addWidget(self._sel_section)
        side.addStretch(1)
        return body

    def _build_pair_side(self, make_help_button) -> QWidget:
        body = QWidget(self)
        side = QVBoxLayout(body)
        side.setContentsMargins(12, 10, 12, 10)
        side.setSpacing(12)

        sec = _CollapsibleSection("Pairing", body, expanded=True)
        sec.add_widget(_muted("Paired by length, centromere index and signal pattern.", body))
        self._pair_mode = NoWheelComboBox(self)
        self._pair_mode.addItem("Automatic", "auto")
        self._pair_mode.addItem("By hand", "hand")
        self._pair_mode.setCurrentIndex(max(0, self._pair_mode.findData(_read("signals_pair_mode", "auto"))))
        self._pair_mode.currentIndexChanged.connect(lambda _i: _write("signals_pair_mode", self._pair_mode.currentData()))
        sec.add_layout(_form_row(body, "Pair", self._pair_mode))
        self._pair_status = QLabel(body)
        self._pair_status.setWordWrap(True)
        self._pair_status.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        sec.add_widget(self._pair_status)
        row = QHBoxLayout()
        self._pair_again = QPushButton("Pair Again", body)
        self._pair_again.clicked.connect(lambda: self._pair(force=True))
        self._pair_accept = QPushButton("Accept Pairs", body)
        self._pair_accept.clicked.connect(self._accept_pairs)
        row.addWidget(self._pair_again, 1)
        row.addWidget(self._pair_accept, 1)
        row.addWidget(make_help_button(body, "Pairing", HELP_PAIR))
        sec.add_layout(row)
        side.addWidget(sec)

        sec = _CollapsibleSection("Across cells", body, expanded=True)
        self._scope = NoWheelComboBox(self)
        for key, text in (("slide", "This slide"), ("cell", "This cell"), ("all", "All slides")):
            self._scope.addItem(text, key)
        self._scope.setCurrentIndex(max(0, self._scope.findData(_read("signals_scope", "slide"))))
        self._scope.currentIndexChanged.connect(self._on_scope)
        self._seen = NoWheelComboBox(self)
        for pct in (50, 25, 75, 100):
            self._seen.addItem(f"{pct} %", pct / 100.0)
        self._seen.currentIndexChanged.connect(lambda _i: self._sync_pair_table())
        sec.add_layout(_form_row(body, "Cells", self._scope))
        sec.add_layout(_form_row(body, "Seen in at least", self._seen))
        side.addWidget(sec)

        sec = _CollapsibleSection("Signals per chromosome", body, expanded=True)
        self._per_note = _muted("", body)
        sec.add_widget(self._per_note)
        self._per = QTableWidget(0, 1, body)
        self._per.verticalHeader().setVisible(False)
        self._per.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._per.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self._per.setShowGrid(False)
        ph = self._per.horizontalHeader()
        ph.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        ph.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self._per.cellClicked.connect(self._on_per_clicked)
        sec.add_widget(self._per)
        self._per_warn = QLabel(body)
        self._per_warn.setWordWrap(True)
        self._per_warn.setStyleSheet("color: #9A5B00;")
        self._per_warn.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        sec.add_widget(self._per_warn)
        side.addWidget(sec)

        sec = _CollapsibleSection("Karyotype", body, expanded=True)
        sec.add_widget(_muted("Mean lengths, arms and the signals seen in at least the chosen share of cells become the "
                              "karyotype.", body))
        row = QHBoxLayout()
        self._k_new = QPushButton("New from Cells…", body)
        self._k_new.clicked.connect(lambda: self._karyotype_from_cells(update=False))
        self._k_update = QPushButton("Update", body)
        self._k_update.clicked.connect(lambda: self._karyotype_from_cells(update=True))
        row.addWidget(self._k_new, 1)
        row.addWidget(self._k_update, 1)
        sec.add_layout(row)
        side.addWidget(sec)
        side.addStretch(1)
        row = QHBoxLayout()
        b = QPushButton("Show on Idiogram", body)
        b.clicked.connect(self.idiogramRequested)
        f = QPushButton("Add to Figure", body)
        f.clicked.connect(lambda: self.figureRequested.emit(self.image_id()) if self.image_id() else None)
        row.addWidget(b, 1)
        row.addWidget(f, 1)
        side.addLayout(row)
        save = QPushButton("Save ▾", body)
        menu = QMenu(save)
        menu.addAction("Summary (TSV)…", self._save_summary)
        menu.addAction("Chromosome Measurements (TSV)…", self._save_measurements)
        save.setMenu(menu)
        side.addWidget(save)
        return body

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
        elif not (s.dirty or self._save_timer.isActive() or self._busy):
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
        if self._folder is None or self._busy:
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
            self._undo.clear()
            self._redo.clear()
            self._strip_cache.clear()
            self._pairing_checks = {}
        message = ""
        qimage = None
        try:
            qimage = rgb_to_qimage(self._session.view.merged()) if self._session.view.is_multichannel else None
            if qimage is None:
                qimage = QImage(str(ci.image_file_path(self._folder, image)))
                if qimage.isNull():
                    qimage, message = None, "The picture file is missing."
                else:
                    self._session.set_picture(_qimage_rgb(qimage))
        except Exception as exc:  # noqa: BLE001
            message = f"Cannot show this image:\n{exc}"
        self._canvas.set_image(qimage, keep_view=keep_view and previous == image_id, message=message)
        self._overlay_key = None
        self._fill_karyotypes()
        self._fill_detect()
        self._sync_all()
        self._refresh_overlay(force=True)
        if previous != image_id:
            self.imageChanged.emit(image_id)

    def _show_empty(self) -> None:
        self._canvas.session = None
        self._canvas.overlay = None
        self._canvas.status = ""
        self._canvas.set_image(None, keep_view=False, message="No image selected. Choose Mark signals in the Library.")
        clear_layout(self._sel_lay)
        self._probes.setRowCount(0)
        self._nav.setText("")
        self._strips.set_groups([], {})
        self._profile.tracks = []
        self._profile.update()

    # ------------------------------------------------------------------ steps, tools, view

    def _set_step(self, step: str) -> None:
        self._step = step if step in dict(STEPS) else "chromosomes"
        _write("signals_step", self._step)
        for key, b in self._step_buttons.items():
            b.setChecked(key == self._step)
        pair = self._step == "pair"
        self._tools_host.setVisible(not pair)
        self._show.setVisible(not pair)
        self._sep.setVisible(False)
        self._view_host.setVisible(pair)
        self._sort_label.setVisible(pair)
        self._sort.setVisible(pair)
        self._sides.setCurrentIndex(1 if pair else 0)
        if pair:
            self._views.setCurrentIndex(1 if self._view == "strips" else 0)
            self._pair()
        else:
            self._views.setCurrentIndex(0)
        allowed = STEP_TOOLS.get(self._step, {"select"})
        for key, b in self._tool_buttons.items():
            b.setVisible(key in allowed)  # each step shows only its own tools
        if self._tool not in allowed:
            self._set_tool("select")
        self._fill_probe_form()
        self._sync_status()

    def _set_tool(self, tool: str) -> None:
        self._tool = tool if tool in dict(TOOLS) else "select"
        for key, b in self._tool_buttons.items():
            b.setChecked(key == self._tool)
        self._canvas.picked = []
        cursor = {"split": Qt.CursorShape.CrossCursor, "brush": Qt.CursorShape.BlankCursor,
                  "erase": Qt.CursorShape.BlankCursor, "add": Qt.CursorShape.CrossCursor}.get(self._tool, Qt.CursorShape.ArrowCursor)
        self._canvas.setCursor(cursor)
        self._canvas.tool_cursor = cursor
        self._canvas.brush_at = None
        self._refresh_overlay(force=True)

    def _set_view(self, view: str) -> None:
        self._view = view if view in ("image", "strips") else "strips"
        _write("signals_view", self._view)
        for key, b in self._view_buttons.items():
            b.setChecked(key == self._view)
        if getattr(self, "_step", "") == "pair":
            self._views.setCurrentIndex(1 if self._view == "strips" else 0)

    def _on_show(self, _index: int) -> None:
        self._canvas.show_mode = str(self._show.currentData())
        _write("signals_show", self._canvas.show_mode)
        self._canvas.update()

    # ------------------------------------------------------------------ undo

    def _snapshot(self) -> dict:
        d = self._session.details
        return {"chromosomes": [m.to_dict() for m in d.chromosomes], "signals": [s.to_dict() for s in d.signals],
                "detection": copy.deepcopy(d.detection), "channels": [c.to_dict() for c in d.channels]}

    def _restore(self, snap: dict) -> None:
        from repeatviewer.backends.cytogenetics.model import Channel

        d = self._session.details
        d.chromosomes = [m for m in (ChromosomeMark.from_dict(x) for x in snap["chromosomes"]) if m]
        d.signals = [s for s in (CytoSignal.from_dict(x) for x in snap["signals"]) if s]
        d.detection = copy.deepcopy(snap["detection"])
        d.channels = [c for c in (Channel.from_dict(x) for x in snap["channels"]) if c]
        self._session.dirty = True

    def _remember(self) -> None:
        if self._session is None:
            return
        self._undo.append(self._snapshot())
        del self._undo[:-UNDO_STEPS]
        self._redo.clear()

    def undo(self) -> None:
        if self._session is None or not self._undo or self._busy:
            return
        self._redo.append(self._snapshot())
        self._restore(self._undo.pop())
        self._select(None)
        self._changed(structure=True)

    def redo(self) -> None:
        if self._session is None or not self._redo or self._busy:
            return
        self._undo.append(self._snapshot())
        self._restore(self._redo.pop())
        self._select(None)
        self._changed(structure=True)

    # ------------------------------------------------------------------ background work

    def _run(self, label: str, work, done) -> None:
        """Run ``work()`` on a worker thread (never the GUI thread), then ``done(result)`` here."""
        if self._busy:
            return
        self._busy = True
        self._canvas.status = f"{label}…"
        self._canvas.update()
        for b in (self._btn_detect, self._btn_signals, self._btn_detect_all):
            b.setEnabled(False)

        def finish(result) -> None:
            self._busy = False
            for b in (self._btn_detect, self._btn_signals, self._btn_detect_all):
                b.setEnabled(True)
            if isinstance(result, Exception):
                QMessageBox.warning(self, "Signals", f"{label} failed:\n{result}")
                self._sync_status()
                return
            done(result)

        def target() -> None:
            try:
                result = work()
            except Exception as exc:  # noqa: BLE001
                result = exc
            self._relay.done.emit(finish, result)

        threading.Thread(target=target, daemon=True).start()

    # ------------------------------------------------------------------ detection

    def _fill_detect(self) -> None:
        s = self._session
        self._det_from.blockSignals(True)
        self._det_from.clear()
        if s is not None:
            if s.details.channels:
                from repeatviewer.backends.cytogenetics.signals_detect import short_channel_name

                self._det_from.addItem("Automatic", "")
                for i, c in enumerate(s.details.channels):
                    self._det_from.addItem(f"{i + 1} · {c.dye or short_channel_name(c.original_name or c.file)}", c.id)
                chosen = s.details.detection.get("counterstain", "")
                self._det_from.setCurrentIndex(max(0, self._det_from.findData(chosen)))
            else:
                self._det_from.addItem("Blue of the picture", "")
            st = s.detect_settings()
            self._det_sens.setValue(int(round(st.sensitivity * 100)))
            self._det_sep.setValue(int(round(st.separate * 100)))
            self._det_small.setText(f"{st.smallest_px} px" if st.smallest_px else "")
        self._det_from.blockSignals(False)

    def _on_det_from(self, _index: int) -> None:
        s = self._session
        if s is None or not s.details.channels:
            return
        self._remember()
        s.choose_counterstain(str(self._det_from.currentData() or ""))
        self._changed(structure=False)

    def _detect_values(self) -> sm.DetectSettings:
        small = int(_number(self._det_small.text(), 0))
        return sm.DetectSettings(self._det_sens.value() / 100.0, self._det_sep.value() / 100.0, small)

    def _remember_detect(self) -> None:
        if self._session is not None:
            self._session.set_detect_settings(self._detect_values())
            self._save_timer.start()

    def _detect_chromosomes(self) -> None:
        s = self._session
        if s is None:
            return
        if any(m.named for m in s.details.chromosomes):
            answer = QMessageBox.question(self, "Detect Chromosomes",
                                          "Detect again? The chromosomes of this cell are replaced, with their names.")
            if answer != QMessageBox.StandardButton.Yes:
                return
        if s.counterstain_index() is None and s.picture is None:
            QMessageBox.information(self, "Detect Chromosomes", "This cell has no counterstain image.")
            return
        settings = self._detect_values()
        guess = bool(s.details.channels) and not s.counterstain_chosen()

        def work():
            idx = s.guess_counterstain() if guess else None
            if idx is not None:
                grey = np.asarray(cyto_channels.to_grey(s._layers()[idx][0]), dtype=float)
            else:
                grey = s.counterstain()
            return idx, sm.detect_chromosomes(grey, settings) if grey is not None else []

        def done(result) -> None:
            if self._session is not s:
                return
            idx, found = result
            if not found:
                name = self._det_from.currentText() if not guess or idx is None else \
                    f"{idx + 1} · {s.details.channels[idx].dye or s.details.channels[idx].original_name}"
                QMessageBox.information(self, "Detect Chromosomes",
                                        f"No chromosomes found in channel {name}.\nChoose the channel that shows the "
                                        "chromosomes in Chromosomes from, or move Sensitivity up.")
                self._sync_status()
                return
            self._remember()
            if idx is not None:
                s.choose_counterstain(s.details.channels[idx].id)  # remembered; change it in Chromosomes from
                self._fill_detect()
            n = s.apply_chromosomes(found, settings)
            self._select(None)
            self._changed(structure=True)
            self._canvas.status = f"{n} chromosomes found"
            self._canvas.update()

        self._run("Finding chromosomes", work, done)

    def _detect_all_signals(self) -> None:
        s = self._session
        if s is None:
            return
        rows = [r for r in s.probe_rows() if r.probe]
        if not rows:
            QMessageBox.information(self, "Detect Signals", "Put a probe on a channel first (Probes, +).")
            return
        jobs = []
        labels = s.labels_map()
        lab = labels if labels is not None and labels.any() else None
        for r in rows:
            arr = s.probe_array(r)
            if arr is not None:
                jobs.append((r, arr, s.spot_settings(r.probe)))

        def work():
            return [(r, sp.find_spots(arr, st, lab), st) for r, arr, st in jobs]

        def done(results) -> None:
            if self._session is not s:
                return
            self._remember()
            for r, found, st in results:
                s.apply_spots(r, found, st)
            self._changed(structure=False)

        self._run("Finding signals", work, done)

    def _detect_selected_probe(self) -> None:
        s = self._session
        row = self._current_probe_row()
        if s is None or row is None or not row.probe or self._busy:
            return
        arr = s.probe_array(row)
        if arr is None:
            return
        st = self._probe_values()
        labels = s.labels_map()
        lab = labels if labels is not None and labels.any() else None

        def done(found) -> None:
            if self._session is not s:
                return
            self._remember()
            s.apply_spots(row, found, st)
            self._changed(structure=False)

        self._run(f"Finding {row.probe} signals", lambda: sp.find_spots(arr, st, lab), done)

    # ------------------------------------------------------------------ probes table

    def _current_probe_row(self) -> "ProbeRow | None":
        s = self._session
        if s is None:
            return None
        return next((r for r in s.probe_rows() if r.key == self._probe_key), None)

    def _fill_probes(self) -> None:
        s = self._session
        self._probes.blockSignals(True)
        rows = s.probe_rows() if s else []
        self._probes.setRowCount(len(rows))
        self._canvas.colours = {r.probe: r.colour for r in rows if r.probe}
        self._canvas.hidden = {r.probe for r in rows if r.probe and not r.visible}
        selected_row = -1
        for i, r in enumerate(rows):
            show = QTableWidgetItem()
            show.setData(Qt.ItemDataRole.UserRole, r.key)
            if r.probe:
                show.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable | Qt.ItemFlag.ItemIsUserCheckable)
                show.setCheckState(Qt.CheckState.Checked if r.visible else Qt.CheckState.Unchecked)
                show.setIcon(_swatch(r.colour))
            else:
                show.setFlags(Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable)
            probe = QTableWidgetItem(r.probe or "No probe")
            if not r.probe:
                f = probe.font()
                f.setItalic(True)
                probe.setFont(f)
                probe.setForeground(QBrush(QColor("#8A8D93")))
            chan = QTableWidgetItem(r.channel_text)
            count = QTableWidgetItem(str(r.count) if r.probe else "")
            for c, item in enumerate((show, probe, chan, count)):
                item.setTextAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
                self._probes.setItem(i, c, item)
            if r.key == self._probe_key:
                selected_row = i
        if selected_row < 0 and rows:
            first = next((i for i, r in enumerate(rows) if r.probe), 0)
            selected_row = first
            self._probe_key = rows[first].key
        if selected_row >= 0:
            self._probes.selectRow(selected_row)
        self._probes.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._probes.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        QTimer.singleShot(0, self._size_probe_columns)
        self._probes.resizeRowsToContents()
        h = self._probes.horizontalHeader().sizeHint().height() + sum(self._probes.rowHeight(r) for r in range(len(rows))) + 4
        self._probes.setFixedHeight(max(60, h))
        self._probes.blockSignals(False)
        self._fill_probe_form()

    def _size_probe_columns(self) -> None:
        """First sizes share the panel width; afterwards the user's own column widths stay."""
        if getattr(self, "_probe_cols_user", False):
            return
        w = self._probes.viewport().width()
        if w < 120:
            return
        self._sizing_cols = True
        rest = w - 38 - 30
        for col, width in enumerate((38, int(rest * 0.45), rest - int(rest * 0.45), 30)):
            self._probes.setColumnWidth(col, width)
        self._sizing_cols = False

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        from PySide6.QtCore import QEvent

        if obj is self._probes and event.type() == QEvent.Type.Resize:
            QTimer.singleShot(0, self._size_probe_columns)
        return super().eventFilter(obj, event)

    def _probe_col_dragged(self, *_args) -> None:
        if not getattr(self, "_sizing_cols", False) and self._probes.isVisible():
            self._probe_cols_user = True  # the user's own widths stay

    def _on_probe_item(self, item: QTableWidgetItem) -> None:
        s = self._session
        if s is None or item.column() != 0:
            return
        key = item.data(Qt.ItemDataRole.UserRole)
        row = next((r for r in s.probe_rows() if r.key == key), None)
        if row is None or not row.probe:
            return
        s.set_probe_visible(row.probe, item.checkState() == Qt.CheckState.Checked)
        self._canvas.hidden = {r.probe for r in s.probe_rows() if r.probe and not r.visible}
        self._canvas.update()
        self._save_timer.start()

    def _on_probe_selected(self) -> None:
        items = self._probes.selectedItems()
        if not items:
            return
        key = self._probes.item(items[0].row(), 0).data(Qt.ItemDataRole.UserRole)
        if key != self._probe_key:
            self._probe_key = key
            self._fill_probe_form()

    def _fill_probe_form(self) -> None:
        s = self._session
        row = self._current_probe_row()
        self._probe_sec.setVisible(row is not None and bool(row.probe) and self._step != "pair")
        if s is None or row is None or not row.probe:
            return
        self._probe_sec._title.setText(row.probe)
        lists = lists_mod.load_lists(self._folder)
        for w in (self._p_channel, self._p_dye, self._p_sens):
            w.blockSignals(True)
        self._p_channel.clear()
        for r in s.probe_rows():
            if r.key == row.key or not r.probe:
                self._p_channel.addItem(f"{r.number} · {r.file}", r.key)
        self._p_channel.setToolTip(row.file)
        self._p_channel.setCurrentIndex(max(0, self._p_channel.findData(row.key)))
        self._p_dye.clear()
        for name in lists_mod.names(lists, "dyes"):
            self._p_dye.addItem(_swatch(lists_mod.dye_color(lists, name, "#ffffff")), name, name)
        if row.dye and self._p_dye.findData(row.dye) < 0:
            self._p_dye.addItem(row.dye, row.dye)
        self._p_dye.addItem("Add New Dye…", "__add__")
        self._p_dye.setCurrentIndex(max(0, self._p_dye.findData(row.dye)))
        self._p_colour.setIcon(_swatch(row.colour, 14))
        self._p_colour.setText(QColor(row.colour).name().upper())
        st = s.spot_settings(row.probe)
        self._p_sens.setValue(int(round(st.sensitivity * 100)))
        self._p_sens_value.setText(str(self._p_sens.value()))
        self._p_small.setText(f"{st.smallest_px:g} px")
        self._p_join.setText(f"{st.join_px:g} px")
        for w in (self._p_channel, self._p_dye, self._p_sens):
            w.blockSignals(False)

    def _probe_values(self) -> sp.SpotSettings:
        return sp.SpotSettings.from_dict({"sensitivity": self._p_sens.value() / 100.0,
                                          "smallest_px": _number(self._p_small.text(), 2.0),
                                          "join_px": _number(self._p_join.text(), 3.0)})

    def _on_probe_settings(self, *_args) -> None:
        self._p_sens_value.setText(str(self._p_sens.value()))
        row = self._current_probe_row()
        if self._session is None or row is None or not row.probe:
            return
        self._session.set_spot_settings(row.probe, self._probe_values())
        if self._p_slide.isChecked():
            self._session.copy_settings_to_slide()
        self._spot_timer.start()  # the count follows the slider
        self._save_timer.start()

    def _on_same_slide(self, on: bool) -> None:
        _write("signals_same_slide", "true" if on else "false")
        if on and self._session is not None:
            self._session.copy_settings_to_slide()
            self._save_timer.start()

    def _on_probe_channel(self, _index: int) -> None:
        s = self._session
        row = self._current_probe_row()
        key = self._p_channel.currentData()
        if s is None or row is None or not key or key == row.key:
            return
        self._remember()
        s.remove_probe(row.key)
        s.assign_probe(key, row.probe, row.dye, row.colour)
        self._probe_key = key
        self._changed(structure=False)

    def _on_probe_dye(self, _index: int) -> None:
        s = self._session
        row = self._current_probe_row()
        dye = self._p_dye.currentData()
        if s is None or row is None:
            return
        if dye == "__add__":
            from PySide6.QtWidgets import QInputDialog

            name, ok = QInputDialog.getText(self, "Add New Dye", "Name of the dye:")
            name = name.strip() if ok else ""
            if not name:
                self._fill_probe_form()
                return
            lists = lists_mod.load_lists(self._folder)
            lists_mod.add_entry(lists, "dyes", name)
            lists_mod.save_lists(self._folder, lists)
            dye = name
        self._remember()
        lists = lists_mod.load_lists(self._folder)
        s.assign_probe(row.key, row.probe, dye, lists_mod.dye_color(lists, dye, row.colour))
        self._changed(structure=False)

    def _on_probe_colour(self) -> None:
        s = self._session
        row = self._current_probe_row()
        if s is None or row is None:
            return
        c = QColorDialog.getColor(QColor(row.colour), self, f"Colour of {row.probe}")
        if c.isValid():
            self._remember()
            s.assign_probe(row.key, row.probe, row.dye, c.name())
            self._changed(structure=False)

    def _probe_add_menu(self) -> None:
        s = self._session
        if s is None:
            return
        menu = QMenu(self)
        lists = lists_mod.load_lists(self._folder)
        used = {r.probe for r in s.probe_rows() if r.probe}
        free = s.free_rows() if s.details.channels else [None] if s.free_planes() else []
        head = menu.addAction("From your lists")
        head.setEnabled(False)
        any_listed = False
        for e in lists.get("probes", []):
            name = str(e.get("name") or "")
            if not name or name in used:
                continue
            dye = str(e.get("dye") or "")
            act = menu.addAction(_swatch(lists_mod.dye_color(lists, dye, "#ffffff")), f"{name} · {dye}" if dye else name)
            act.setEnabled(bool(free))
            act.triggered.connect(lambda _c=False, n=name, d=dye: self._add_listed_probe(n, d))
            any_listed = True
        if not any_listed:
            none = menu.addAction("No other probes in the lists")
            none.setEnabled(False)
        menu.addSeparator()
        menu.addAction("New Probe…", self._new_probe)
        menu.addAction("Edit Lists…", self._edit_lists)
        menu.exec(self._probe_add.mapToGlobal(self._probe_add.rect().bottomLeft()))

    def _first_free_key(self) -> str:
        s = self._session
        if s.details.channels:
            free = s.free_rows()
            return free[0].key if free else ""
        planes = s.free_planes()
        return f"plane:{planes[0]}" if planes else ""

    def _add_listed_probe(self, name: str, dye: str) -> None:
        s = self._session
        key = self._first_free_key() if s else ""
        if not key:
            QMessageBox.information(self, "Probes", "Every channel of this cell already has a probe.")
            return
        self._remember()
        s.assign_probe(key, name, dye)
        self._probe_key = key
        if self._p_slide.isChecked():
            s.copy_settings_to_slide()
        self._changed(structure=False)

    def _new_probe(self) -> None:
        s = self._session
        if s is None or self._folder is None:
            return
        from repeatviewer.views.cyto_new_probe_dialog import NewProbeDialog

        if s.details.channels:
            channels = [(r.key, f"{r.number} · {r.file}") for r in s.free_rows()]
        else:
            channels = [(f"plane:{p}", PLANES[p]) for p in s.free_planes()]
        slide = next((sl.name for sl in load_slides(self._folder) if sl.id == s.details.slide_id), "")
        dialog = NewProbeDialog(self._folder, channels, slide, self.window())
        dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)

        def done(accepted: int) -> None:
            if not accepted or self._session is not s:
                return
            self._remember()
            if dialog.save_to_list():
                s.add_probe_to_lists(dialog.name(), dialog.dye(), dialog.sequence())
            s.assign_probe(dialog.channel_key(), dialog.name(), dialog.dye(), dialog.colour())
            self._probe_key = dialog.channel_key()
            if dialog.every_cell():
                s.copy_settings_to_slide()
            self._changed(structure=False)

        dialog.finished.connect(done)
        dialog.show()  # not modal: a Files click fills the sequence

    def _edit_lists(self) -> None:
        if self._folder is None:
            return
        try:
            from repeatviewer.views.cyto_lists_dialog import ListsDialog

            dialog = ListsDialog(self._folder, self.window())
            dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
            dialog.finished.connect(lambda _r: self._fill_probes())
            dialog.show()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.information(self, "Edit Lists", f"Open Tools > Cytogenetics > Edit Lists.\n({exc})")

    def _remove_probe(self) -> None:
        s = self._session
        row = self._current_probe_row()
        if s is None or row is None or not row.probe:
            return
        n = sum(1 for x in s.details.signals if x.probe == row.probe)
        if n:
            answer = QMessageBox.question(self, "Remove Probe", f"Take {row.probe} off channel {row.number}? "
                                          f"Its {n} signal{'s' if n != 1 else ''} in this cell are removed too.")
            if answer != QMessageBox.StandardButton.Yes:
                return
        self._remember()
        s.remove_probe(row.key)
        self._changed(structure=False)

    # ------------------------------------------------------------------ canvas events

    def _tolerance(self, pixels: float) -> float:
        return pixels / max(self._canvas.zoom(), 1e-6)

    def _brush_radius(self) -> float:
        s = self._session
        widths = [pr.half_width(m) for m in s.details.chromosomes if m.mask] if s else []
        return max(2.5, 0.8 * float(np.median(widths))) if widths else 5.0

    def _handle_at(self, x: float, y: float) -> tuple[str, str] | None:
        """("centromere" | "telomere_p" | "telomere_q", mark id) of the selected chromosome near a point."""
        s = self._session
        if s is None or not self._selected or self._selected[0] != "mark":
            return None
        m = s.details.chromosome(self._selected[1])
        if m is None:
            return None
        tol = self._tolerance(9)
        for which in ("centromere", "telomere_p", "telomere_q"):
            pt = getattr(m, which)
            if pt is not None and math.hypot(pt[0] - x, pt[1] - y) <= tol:
                return which, m.id
        return None

    def _on_press(self, x: float, y: float, modifiers) -> None:
        s = self._session
        if s is None or self._busy:
            return
        self._drag = None
        self._moved = False
        tool = self._tool
        if self._step == "pair":
            tool = "select"
        if tool == "select":
            handle = self._handle_at(x, y)
            if handle:
                self._remember()
                self._drag = handle
                return
            hit = s.hit(x, y, self._tolerance(10))
            if hit and hit[0] == "signal":
                self._remember()
                self._drag = hit
            if self._step == "pair" and hit and hit[0] == "mark" and self._selected and self._selected[0] == "mark" \
                    and bool(modifiers & Qt.KeyboardModifier.ShiftModifier):
                self._swap_names(self._selected[1], hit[1])
                return
            self._select(hit)
        elif tool == "split":
            self._canvas.line = (QPointF(x, y), QPointF(x, y))
        elif tool == "merge":
            m = s.mark_at(x, y)
            if m is None:
                return
            if self._canvas.picked and self._canvas.picked[0] != m.id:
                self._remember()
                keep = s.merge_marks([self._canvas.picked[0], m.id])
                self._canvas.picked = [keep]
                self._select(("mark", keep))
                self._changed(structure=True)
            else:
                self._canvas.picked = [m.id]
                self._select(("mark", m.id))
        elif tool in ("brush", "erase"):
            m = s.mark_at(x, y)
            if tool == "erase" and m is None:
                return
            self._remember()
            r = self._brush_radius()
            ident = s.brush_mark(m.id if m else "", x, y, r, tool == "brush", measure_now=False)
            self._drag = ("brush", ident)
            self._refresh_overlay(force=True)
        elif tool == "add":
            row = self._current_probe_row()
            if row is None or not row.probe:
                QMessageBox.information(self, "Add signal", "Select a probe in the Probes table first.")
                return
            self._remember()
            sig = s.add_signal(row.probe, x, y, refine=True)
            self._select(("signal", sig.id))
            self._changed(structure=False)

    def _on_move(self, x: float, y: float) -> None:
        s = self._session
        if s is None:
            return
        if self._tool in ("brush", "erase") and self._step != "pair":
            self._canvas.brush_at = (x, y, self._brush_radius())
        if self._canvas.line is not None:
            self._canvas.line = (self._canvas.line[0], QPointF(x, y))
            self._canvas.update()
            return
        if self._drag is None:
            self._canvas.update()
            return
        kind, ident = self._drag
        self._moved = True
        if kind == "centromere":
            s.move_centromere(ident, x, y)
        elif kind in ("telomere_p", "telomere_q"):
            s.move_end(ident, kind, x, y)
        elif kind == "signal":
            s.move_signal(ident, x, y)
        elif kind == "brush" and ident:
            s.brush_mark(ident, x, y, self._brush_radius(), self._tool == "brush", measure_now=False)
            self._refresh_overlay(force=True)
        self._canvas.update()

    def _on_release(self, x: float, y: float) -> None:
        s = self._session
        if s is None:
            return
        if self._canvas.line is not None:
            a, b = self._canvas.line
            self._canvas.line = None
            if math.hypot(b.x() - a.x(), b.y() - a.y()) > 3:
                self._remember()
                ids = s.split_at_line((a.x(), a.y()), (b.x(), b.y()))
                if not ids:
                    self._undo.pop()
                self._changed(structure=True)
            self._canvas.update()
            return
        drag, self._drag = self._drag, None
        if drag is None:
            return
        if drag[0] == "brush":
            if drag[1] and s.details.chromosome(drag[1]) is not None:
                s.remeasure(drag[1])
                s._renumber_unnamed()
                self._select(("mark", drag[1]))
            self._changed(structure=True)
        elif self._moved:
            self._changed(structure=False)
        elif drag[0] in ("centromere", "telomere_p", "telomere_q", "signal"):
            self._undo.pop() if self._undo else None

    def _on_hover(self, x: float, y: float, global_pos) -> None:
        s = self._session
        if s is None:
            return
        if self._tool in ("brush", "erase") and self._step != "pair":
            self._canvas.brush_at = (x, y, self._brush_radius())
            self._canvas.update()
        hit = s.hit(x, y, self._tolerance(8))
        text = self._tooltip(hit) if hit else ""
        if text:
            QToolTip.showText(global_pos, text, self._canvas)
        else:
            QToolTip.hideText()

    def _mark_title(self, m: ChromosomeMark) -> str:
        s = self._session
        name = m.label if m.named else f"Chromosome {m.homologue} (not paired)"
        parts = [name]
        um = s.length_um(m)
        if um is not None:
            parts.append(f"{um:.1f} µm")
        elif m.midline:
            parts.append(f"{sm.polyline_length(m.midline):.0f} px")
        ci_val = s.centromere_index(m)
        if ci_val is not None:
            from repeatviewer.backends.cytogenetics.karyotype import levan_type

            parts.append(f"CI {ci_val:.0f} ({levan_type(ci_val)})")
        return " · ".join(parts)

    def _position_text(self, sig) -> str:
        if sig.rel_distance is None or not sig.arm:
            return "not on a chromosome" if not sig.chromosome_id else "position not measured"
        return "at the centromere" if sig.arm == "cen" else f"{sig.arm} {sig.rel_distance:.2f}"

    def _tooltip(self, hit) -> str:
        s = self._session
        kind, ident = hit
        if kind == "mark":
            m = s.details.chromosome(ident)
            if m is None:
                return ""
            lines = [f"<b>{self._mark_title(m)}</b>"]
            for sig in s.details.signals:
                if sig.chromosome_id == m.id:
                    lines.append(f"{sig.probe}: {self._position_text(sig)}")
            if m.check:
                lines.append(f"<span style='color:#9A5B00'>Check: {m.check}</span>")
            return "<br>".join(lines)
        sig = s.signal(ident)
        if sig is None:
            return ""
        m = s.details.chromosome(sig.chromosome_id)
        lines = [f"<b>{sig.probe}</b>", f"{m.label if m and m.named else 'chromosome ' + (m.homologue if m else '?')}: "
                 f"{self._position_text(sig)}"]
        if sig.strength:
            lines.append(sig.strength)
        if sig.check:
            lines.append(f"<span style='color:#9A5B00'>Check: {sig.check}</span>")
        return "<br>".join(lines)

    def _on_double(self, x: float, y: float) -> None:
        s = self._session
        if s is None:
            return
        m = s.mark_at(x, y)
        if m is None:
            self._canvas.fit()
        elif self._step == "pair":
            self._open_karyotype(edit=True)

    def _on_menu(self, x: float, y: float, global_pos) -> None:
        s = self._session
        if s is None:
            return
        hit = s.hit(x, y, self._tolerance(10))
        menu = QMenu(self)
        if hit is None:
            menu.addAction("Add Chromosome Here", lambda: self._add_chromosome(x, y))
            if self._step == "signals":
                row = self._current_probe_row()
                if row is not None and row.probe:
                    menu.addAction(f"Add {row.probe} Signal Here", lambda: self._add_signal_at(row.probe, x, y))
        elif hit[0] == "mark":
            self._select(hit)
            m = s.details.chromosome(hit[1])
            if m.check:
                menu.addAction("Accept", self._accept_selected)
            menu.addAction("Swap p and q Arms", lambda: self._edit(lambda: s.swap_mark_arms(m.id)))
            if self._step == "pair" and s.karyotype is not None:
                sub = menu.addMenu("Name")
                for label in s.labels():
                    sub.addAction(label, lambda lab=label: self._rename_mark(m.id, lab))
            menu.addSeparator()
            menu.addAction("Remove Chromosome", self._delete_selected)
        else:
            self._select(hit)
            sig = s.signal(hit[1])
            if sig.check:
                menu.addAction("Accept", self._accept_selected)
            sub = menu.addMenu("Probe")
            for r in s.probe_rows():
                if r.probe:
                    act = sub.addAction(_swatch(r.colour), r.probe, lambda p=r.probe: self._edit(lambda: s.update_signal(sig.id, probe=p)))
                    act.setCheckable(True)
                    act.setChecked(r.probe == sig.probe)
            sub = menu.addMenu("Strength")
            for st in STRENGTHS:
                act = sub.addAction(st.capitalize(), lambda v=st: self._edit(lambda: s.update_signal(sig.id, strength=v)))
                act.setCheckable(True)
                act.setChecked(sig.strength == st)
            menu.addSeparator()
            menu.addAction("Remove Signal", self._delete_selected)
        menu.exec(global_pos)

    def _edit(self, fn) -> None:
        self._remember()
        fn()
        self._changed(structure=False)

    def _add_chromosome(self, x: float, y: float) -> None:
        s = self._session
        self._remember()
        ident = s.add_chromosome_at(x, y)
        if not ident:
            self._undo.pop()
            self._canvas.status = "No chromosome found there. Use Brush to paint one."
            self._canvas.update()
            return
        self._select(("mark", ident))
        self._changed(structure=True)

    def _add_signal_at(self, probe: str, x: float, y: float) -> None:
        self._remember()
        sig = self._session.add_signal(probe, x, y, refine=True)
        self._select(("signal", sig.id))
        self._changed(structure=False)

    def _delete_selected(self) -> None:
        s = self._session
        if s is None or not self._selected:
            return
        kind, ident = self._selected
        self._remember()
        (s.remove_mark if kind == "mark" else s.remove_signal)(ident)
        self._select(None)
        self._changed(structure=kind == "mark")

    def _rename_mark(self, mark_id: str, label: str) -> None:
        s = self._session
        other = s.mark_for_label(label)
        self._remember()
        if other is not None and other.id != mark_id:
            pr.swap_names(s.details, mark_id, other.id)
        else:
            parts = ss.split_label(s.karyotype, label)
            m = s.details.chromosome(mark_id)
            if parts and m:
                m.chromosome, m.homologue = parts
        s.dirty = True
        self._changed(structure=False)

    # ------------------------------------------------------------------ selection and review

    def select_signal(self, signal_id: str) -> None:
        s = self._session
        if s is not None and signal_id and s.signal(signal_id) is not None:
            if self._step == "pair":
                self._set_step("signals")
            self._select(("signal", signal_id))

    def _select(self, hit) -> None:
        self._selected = hit
        self._canvas.selected = hit
        self._refresh_overlay(force=False)
        self._canvas.update()
        self._sync_selected()
        self._sync_check()

    def _refresh_overlay(self, force: bool = False) -> None:
        s = self._session
        if s is None:
            self._canvas.overlay = None
            return
        sel = self._selected[1] if self._selected and self._selected[0] == "mark" else ""
        key = (id(s), tuple((m.id, id(m.mask)) for m in s.details.chromosomes))  # outlines are drawn apart
        if force or key != self._overlay_key:
            if force:
                s._labels_cache = None
            self._canvas.overlay = build_overlay(s, sel, self._canvas.picked)
            self._overlay_key = key
        self._canvas.update()

    def _step_check(self, step: int) -> None:
        s = self._session
        if s is None:
            return
        items = s.to_check()
        if not items:
            return
        keys = [(k, i) for k, i, _r in items]
        cur = tuple(self._selected) if self._selected else None
        k = keys.index(cur) + step if cur in keys else (0 if step > 0 else len(keys) - 1)
        k %= len(keys)
        kind, ident = keys[k]
        self._select((kind, ident))
        item = s.details.chromosome(ident) if kind == "mark" else s.signal(ident)
        if item is not None:
            x, y = (item.x, item.y)
            self._canvas.center_on(x, y)

    def _accept_selected(self) -> None:
        s = self._session
        if s is None or not self._selected:
            return
        self._remember()
        s.accept(*self._selected)
        self._changed(structure=False)
        self._step_check(1)

    def _accept_everything(self) -> None:
        s = self._session
        if s is None:
            return
        self._remember()
        s.accept_all()
        self._changed(structure=False)

    def _sync_check(self) -> None:
        s = self._session
        items = s.to_check() if s else []
        for b in (self._check_prev, self._check_next, self._accept_all):
            b.setEnabled(bool(items))
        self._accept_one.setEnabled(bool(self._selected) and bool(items) and
                                    any((k, i) == tuple(self._selected) for k, i, _r in items))
        if s is None:
            self._check_text.setText("")
            return
        if not items:
            self._check_text.setText("Nothing to check in this cell.")
            self._fit_check_text()
            return
        lines = []
        for kind, ident, reason in items[:4]:
            if kind == "mark":
                m = s.details.chromosome(ident)
                lines.append(f"{m.label if m.named else m.homologue}: {reason}")
            else:
                sig = s.signal(ident)
                m = s.details.chromosome(sig.chromosome_id) if sig else None
                where = (m.label if m.named else m.homologue) if m else "no chromosome"
                lines.append(f"{sig.probe} on {where}: {reason}")
        more = len(items) - 4
        self._check_text.setText("\n".join(lines) + (f"\nand {more} more" if more > 0 else ""))
        self._fit_check_text()

    def _fit_check_text(self) -> None:
        """Wrapped label in the side panel: give it the height its text needs at the panel's width."""
        width = max(200, self._sides.width() - 44)
        self._check_text.setMinimumHeight(self._check_text.heightForWidth(width))

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        QTimer.singleShot(0, self._fit_check_text)

    def _sync_selected(self) -> None:
        clear_layout(self._sel_lay)
        s = self._session
        if s is None or not self._selected:
            self._sel_section._title.setText("Selected")
            self._sel_lay.addWidget(_muted("Click a chromosome or a signal.", self._sel_host))
            return
        kind, ident = self._selected
        host = self._sel_host
        if kind == "mark":
            m = s.details.chromosome(ident)
            if m is None:
                return
            self._sel_section._title.setText(f"Selected: {m.label if m.named else 'chromosome ' + m.homologue}")
            name = NoWheelComboBox(host)
            name.addItem("Not paired", "")
            for label in s.labels():
                name.addItem(label, label)
            name.setCurrentIndex(max(0, name.findData(m.label if m.named else "")))
            name.activated.connect(lambda _i, box=name: self._rename_mark(m.id, box.currentData()) if box.currentData() else None)
            name.setEnabled(s.karyotype is not None)
            self._sel_lay.addLayout(_form_row(host, "Name", name))
            um = s.length_um(m)
            length = QLineEdit(f"{um:.2f} µm" if um is not None else (f"{sm.polyline_length(m.midline):.0f} px" if m.midline else ""), host)
            length.setReadOnly(True)
            self._sel_lay.addLayout(_form_row(host, "Length", length))
            ci_val = s.centromere_index(m)
            cie = QLineEdit(f"{ci_val:.1f}" if ci_val is not None else "", host)
            cie.setReadOnly(True)
            self._sel_lay.addLayout(_form_row(host, "Centromere index", cie))
            hint = _muted("Drag the white centromere tick along the midline, or a midline end.", host)
            hint.setMinimumHeight(hint.heightForWidth(250))
            self._sel_lay.addWidget(hint)
            rm = QPushButton("Remove Chromosome", host)
            rm.clicked.connect(self._delete_selected)
            self._sel_lay.addWidget(rm, 0, Qt.AlignmentFlag.AlignLeft)
        else:
            sig = s.signal(ident)
            if sig is None:
                return
            self._sel_section._title.setText(f"Selected: {sig.probe} signal")
            probe = NoWheelComboBox(host)
            for r in s.probe_rows():
                if r.probe:
                    probe.addItem(_swatch(r.colour), r.probe, r.probe)
            probe.setCurrentIndex(max(0, probe.findData(sig.probe)))
            probe.activated.connect(lambda _i, box=probe: self._edit(lambda: s.update_signal(sig.id, probe=box.currentData())))
            self._sel_lay.addLayout(_form_row(host, "Probe", probe))
            m = s.details.chromosome(sig.chromosome_id)
            on = QLineEdit((m.label if m.named else f"chromosome {m.homologue}") if m else "none", host)
            on.setReadOnly(True)
            self._sel_lay.addLayout(_form_row(host, "Chromosome", on))
            pos = QLineEdit(self._position_text(sig), host)
            pos.setReadOnly(True)
            self._sel_lay.addLayout(_form_row(host, "Position", pos))
            strength = NoWheelComboBox(host)
            for st in STRENGTHS:
                strength.addItem(st.capitalize(), st)
            strength.setCurrentIndex(max(0, strength.findData(sig.strength)))
            strength.activated.connect(lambda _i, box=strength: self._edit(lambda: s.update_signal(sig.id, strength=box.currentData())))
            self._sel_lay.addLayout(_form_row(host, "Strength", strength))
            rm = QPushButton("Remove Signal", host)
            rm.clicked.connect(self._delete_selected)
            self._sel_lay.addWidget(rm, 0, Qt.AlignmentFlag.AlignLeft)

    def _changed(self, *, structure: bool) -> None:
        if structure:
            self._strip_cache.clear()
        self._refresh_overlay(force=structure)
        self._sync_all()
        self._save_timer.start()

    def _sync_all(self) -> None:
        self._fill_probes()
        self._sync_selected()
        self._sync_check()
        self._sync_status()
        self._sync_cells()
        if self._step == "pair":
            self._sync_pair()

    def _sync_status(self) -> None:
        s = self._session
        if s is None or self._busy:
            return
        n = len(s.details.chromosomes)
        k = len([x for x in s.details.signals if x.probe not in self._canvas.hidden])
        checks = s.to_check()
        cm = sum(1 for c in checks if c[0] == "mark")
        cs = sum(1 for c in checks if c[0] == "signal")
        parts = [f"{n} chromosome{'s' if n != 1 else ''}", f"{k} signal{'s' if k != 1 else ''}"]
        todo = []
        if cm:
            todo.append(f"{cm} chromosome{'s' if cm != 1 else ''}")
        if cs:
            todo.append(f"{cs} signal{'s' if cs != 1 else ''}")
        if todo:
            parts.append(" and ".join(todo) + " to check")
        if not n:
            parts = ["Press Detect Chromosomes to find the chromosomes"]
        self._canvas.status = " · ".join(parts)
        self._canvas.update()

    # ------------------------------------------------------------------ karyotype

    def _fill_karyotypes(self) -> None:
        s = self._session
        self._karyo.blockSignals(True)
        self._karyo.clear()
        if s is not None:
            for k in s.karyotypes:
                self._karyo.addItem(f"{k.name} {k.formula}", k.name)
            self._karyo.addItem("None (build from cells)", NONE_KARYOTYPE)
            self._karyo.insertSeparator(self._karyo.count())
            self._karyo.addItem("New Karyotype…", NEW_KARYOTYPE)
            self._karyo.addItem("Edit Karyotype…", EDIT_KARYOTYPE)
            idx = self._karyo.findData(s.karyotype.name) if s.karyotype else self._karyo.findData(NONE_KARYOTYPE)
            self._karyo.setCurrentIndex(max(0, idx))
        self._karyo.blockSignals(False)
        self._k_update.setText(f"Update {s.karyotype.name}" if s and s.karyotype else "Update")
        self._k_update.setEnabled(bool(s and s.karyotype))

    def _on_karyotype(self, index: int) -> None:
        s = self._session
        if s is None:
            return
        data = self._karyo.itemData(index)
        if data == NEW_KARYOTYPE:
            self._fill_karyotypes()
            self._open_karyotype(edit=False)
            return
        if data == EDIT_KARYOTYPE:
            self._fill_karyotypes()
            self._open_karyotype(edit=True)
            return
        self._remember()
        s.set_karyotype(None if data == NONE_KARYOTYPE else next((k for k in s.karyotypes if k.name == data), None))
        self._fill_karyotypes()
        self._changed(structure=False)

    def _open_karyotype(self, *, edit: bool, prefilled=None, replaces: str = "") -> None:
        s = self._session
        if s is None or self._folder is None:
            return
        from repeatviewer.views.cyto_karyotype_dialog import KaryotypeDialog

        name = replaces if prefilled is not None else (s.karyotype.name if edit and s.karyotype else "")
        dialog = KaryotypeDialog(self._folder, name, self.window(), karyotype=prefilled)
        dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)

        def done(accepted: int) -> None:
            s2 = self._session
            if s2 is None or self._folder is None:
                return
            s2.karyotypes = list_karyotypes(self._folder)
            if accepted and dialog.saved_name():
                s2.karyotype = next((k for k in s2.karyotypes if k.name == dialog.saved_name()), s2.karyotype)
                s2.details.karyotype = s2.karyotype.name if s2.karyotype else ""
                s2.dirty = True
                self._save_timer.start()
                self.karyotypesChanged.emit()
            self._fill_karyotypes()
            self._sync_all()

        dialog.finished.connect(done)
        dialog.show()  # not modal: the Files panel stays clickable (Lengths from assembly)

    # ------------------------------------------------------------------ step 3

    def _scope_details(self) -> list[ImageDetails]:
        s = self._session
        scope = self._scope.currentData()
        if scope == "cell" or self._folder is None:
            return [s.details]
        out = []
        for image in ci.load_images(self._folder):
            d = s.details if image.id == s.image.id else image_details(image)
            if scope == "slide" and d.slide_id != s.details.slide_id:
                continue
            out.append(d)
        return out

    def _on_scope(self, _index: int) -> None:
        _write("signals_scope", self._scope.currentData())
        self._sync_pair_table()
        self._sync_profile()

    def _pair(self, force: bool = False) -> None:
        s = self._session
        if s is None or not s.details.chromosomes:
            self._sync_pair()
            return
        unnamed = any(not m.named for m in s.details.chromosomes if m.midline)
        if force or (unnamed and self._pair_mode.currentData() == "auto"):
            if not force and any(m.named for m in s.details.chromosomes):
                pass
            self._remember()
            if s.karyotype is not None:
                pairing = pr.pair_with_karyotype(s.details, s.karyotype)
            else:
                pairing = pr.pair_without_karyotype(s.details, 2)
            pr.apply_pairing(s.details, pairing)
            self._pairing_checks = pairing.checks
            s.dirty = True
            self._save_timer.start()
        self._sync_pair()

    def _accept_pairs(self) -> None:
        self._pairing_checks = {}
        self._sync_pair()

    def _swap_names(self, a: str, b: str) -> None:
        s = self._session
        if s is None:
            return
        self._remember()
        pr.swap_names(s.details, a, b)
        s.dirty = True
        for m in (s.details.chromosome(a), s.details.chromosome(b)):
            if m is not None:
                self._pairing_checks.pop(m.chromosome, None)
        self._save_timer.start()
        self._sync_pair()
        self._refresh_overlay()

    def _strip_image(self, m: ChromosomeMark, grey: np.ndarray) -> QImage | None:
        key = f"{m.id}:{id(m.mask)}:{len(m.midline or [])}"
        if key in self._strip_cache:
            return self._strip_cache[key]
        if not m.midline or len(m.midline) < 2:
            return None
        strip = pr.straighten(grey, m.midline, pr.half_width(m) + 1)
        if strip.size == 0:
            return None
        lo, hi = np.percentile(grey[:: max(1, grey.shape[0] // 300), :: max(1, grey.shape[1] // 300)], (1, 99.8))
        v = np.clip((strip - lo) / max(hi - lo, 1e-6), 0, 1)
        rgba = np.zeros(strip.shape + (4,), dtype=np.uint8)
        rgba[..., 0] = (v * 76).astype(np.uint8)
        rgba[..., 1] = (v * 123).astype(np.uint8)
        rgba[..., 2] = (v * 255).astype(np.uint8)
        rgba[..., 3] = 255
        img = _rgba_image(rgba)
        self._strip_cache[key] = img
        return img

    def _sync_pair(self) -> None:
        s = self._session
        if s is None:
            return
        grey = s.counterstain()
        colours = {r.probe: r.colour for r in s.probe_rows() if r.probe}
        groups: dict[str, list[dict]] = {}
        for m in s.details.chromosomes:
            if not m.midline or len(m.midline) < 2:
                continue
            bands = []
            for sig in s.details.signals:
                if sig.chromosome_id == m.id and sig.probe not in self._canvas.hidden:
                    f = pr.signal_fraction(m, sig.x, sig.y)
                    if f is not None:
                        bands.append((f, colours.get(sig.probe, "#FFFFFF")))
            total = sm.polyline_length(m.midline)
            cen = sg.arc_position(m.midline, m.centromere) / total if m.centromere and total else None
            groups.setdefault(m.chromosome if m.named else UNNAMED, []).append({
                "id": m.id, "image": self._strip_image(m, grey) if grey is not None else None, "length": total,
                "cen": cen, "bands": bands, "label": m.homologue.lstrip(".") if m.named else m.homologue})
        order = [c.name for c in s.karyotype.chromosomes] if s.karyotype else []
        rank = {n: i for i, n in enumerate(order)}
        from repeatviewer.backends.cytogenetics.importer import natural_key

        if self._sort.currentData() == "length":
            keyf = lambda kv: (kv[0] == UNNAMED, -max(it["length"] for it in kv[1]))  # noqa: E731
        else:
            keyf = lambda kv: (kv[0] == UNNAMED, rank.get(kv[0], len(rank)), natural_key(kv[0]))  # noqa: E731
        items = sorted(groups.items(), key=keyf)
        for _name, its in items:
            its.sort(key=lambda it: natural_key(it["label"]))
        self._strips.set_groups(items, self._pairing_checks)
        paired = sum(1 for n, _ in items if n != UNNAMED)
        sure = paired - len([c for c in self._pairing_checks if c in groups])
        text = f"{sure} of {paired} pairs sure." if paired else "Not paired yet."
        if self._pairing_checks:
            text += " Check " + ", ".join(sorted(self._pairing_checks, key=natural_key)) + "."
        unp = len(groups.get(UNNAMED, []))
        if unp:
            text += f" {unp} not paired."
        self._pair_status.setText(text)
        self._sync_pair_table()
        self._sync_profile()

    def _sync_pair_table(self) -> None:
        s = self._session
        if s is None:
            return
        cells = self._scope_details()
        order = [c.name for c in s.karyotype.chromosomes] if s.karyotype else None
        self._stats = pr.stats_over_cells(cells, order)
        rows = s.probe_rows()
        probes = [r for r in rows if r.probe]
        n_cells = sum(1 for d in cells if any(m.named for m in d.chromosomes))
        self._per_note.setText(f"Cells with the signal, of {n_cells}. Click a row to see its tracks.")
        self._per.clear()
        self._per.setColumnCount(1 + len(probes))
        self._per.setRowCount(len(self._stats))
        header = [QTableWidgetItem("Chr")] + [QTableWidgetItem(_swatch(r.colour), r.probe) for r in probes]
        for c, item in enumerate(header):
            self._per.setHorizontalHeaderItem(c, item)
        need = float(self._seen.currentData() or 0.5)
        warn = []
        for i, st in enumerate(self._stats):
            self._per.setItem(i, 0, QTableWidgetItem(st.name))
            for j, r in enumerate(probes, start=1):
                n = st.seen.get(r.probe, 0)
                item = QTableWidgetItem(str(n) if n else "·")
                if n and st.cells:
                    share = n / st.cells
                    item.setBackground(QColor("#DCE8FB") if share >= need else QColor("#FFF1D6"))
                    if share < need and share >= need / 2:
                        warn.append(f"{st.name}: {r.probe} in {n} of {st.cells} cells")
                self._per.setItem(i, j, item)
        self._per.setColumnWidth(0, 36)
        width = max(40, (self._per.viewport().width() - 36) // max(1, len(probes))) if self._per.viewport().width() > 100 else 52
        for j, r in enumerate(probes, start=1):
            self._per.setColumnWidth(j, width)
            self._per.horizontalHeaderItem(j).setToolTip(r.probe)
        self._per.resizeRowsToContents()
        h = self._per.horizontalHeader().sizeHint().height() + sum(self._per.rowHeight(r) for r in range(self._per.rowCount())) + 4
        self._per.setFixedHeight(max(60, min(420, h)))
        self._per_warn.setText("\n".join(warn[:3]))

    def _on_per_clicked(self, row: int, _col: int) -> None:
        if 0 <= row < len(getattr(self, "_stats", [])):
            name = self._stats[row].name
            s = self._session
            m = next((m for m in s.details.chromosomes if m.named and m.chromosome == name), None)
            if m is not None:
                self._on_strip_selected(m.id)

    def _on_strip_selected(self, mark_id: str) -> None:
        self._strips.current = mark_id
        self._strips.update()
        self._selected = ("mark", mark_id)
        self._canvas.selected = self._selected
        self._sync_profile()

    def _sync_profile(self) -> None:
        s = self._session
        if s is None:
            return
        mark = s.details.chromosome(self._strips.current) if self._strips.current else None
        if mark is None or not mark.midline:
            self._profile.tracks = []
            self._profile_title.setText("")
            self._profile.update()
            return
        same = [m for m in s.details.chromosomes if m.midline and
                ((m.named and mark.named and m.chromosome == mark.chromosome) or m.id == mark.id)]
        same.sort(key=lambda m: m.homologue)
        title = mark.chromosome if mark.named else f"chromosome {mark.homologue}"
        self._profile_title.setText(f"<b>Chromosome {title} along its length</b>"
                                    f"<span style='color:gray'>  {', '.join(m.label if m.named else m.homologue for m in same[:2])}"
                                    f"{' (solid, dashed)' if len(same) > 1 else ''}</span>")
        stats = next((x for x in getattr(self, "_stats", []) if x.name == mark.chromosome), None)
        tracks = []
        grey = s.counterstain()
        if grey is not None:
            tracks.append({"name": "DAPI", "colour": "#4C7BFF",
                           "lines": [(pr.profile(grey, m.midline, pr.half_width(m)), k > 0) for k, m in enumerate(same[:2])]})
        for r in s.probe_rows():
            if not r.probe:
                continue
            arr = s.probe_array(r)
            if arr is None:
                continue
            lines = [(pr.profile(arr, m.midline, pr.half_width(m)), k > 0) for k, m in enumerate(same[:2])]
            spread, note = None, ""
            if stats is not None and self._profile_avg.isChecked() and stats.positions.get(r.probe):
                pos = stats.positions[r.probe]
                spread = (min(pos), max(pos))
                f = float(np.median(pos))
                cen = stats.ci / 100.0
                arm = "p" if f < cen else "q"
                rel = (cen - f) / cen if arm == "p" and cen else (f - cen) / (1 - cen) if cen < 1 else 0
                note = f"{arm} {rel:.2f} · {stats.seen.get(r.probe, 0)} of {stats.cells} cells"
            tracks.append({"name": r.probe, "colour": r.colour, "lines": lines, "spread": spread, "note": note})
        total = sm.polyline_length(mark.midline)
        self._profile.cen = sg.arc_position(mark.midline, mark.centromere) / total if mark.centromere and total else None
        self._profile.tracks = tracks
        self._profile.update()

    def _karyotype_from_cells(self, *, update: bool) -> None:
        s = self._session
        if s is None or self._folder is None:
            return
        stats = pr.stats_over_cells(self._scope_details(), [c.name for c in s.karyotype.chromosomes] if s.karyotype else None)
        if not stats:
            QMessageBox.information(self, "Karyotype from cells", "Pair the chromosomes of at least one cell first.")
            return
        colours = {r.probe: r.colour for r in s.probe_rows() if r.probe}
        need = float(self._seen.currentData() or 0.5)
        n_cells = max((st.cells for st in stats), default=0)
        if update and s.karyotype is not None:
            k = pr.karyotype_from_cells(stats, s.karyotype.name, seen_at_least=need, colours=colours, base=s.karyotype)
            self._open_karyotype(edit=True, prefilled=k, replaces=s.karyotype.name)
        else:
            slide = next((sl for sl in load_slides(self._folder) if sl.id == s.details.slide_id), None)
            species = s.image.metadata.get("species") or (slide.species if slide else "")
            base_name = species or (slide.name if slide else "Karyotype")
            k = pr.karyotype_from_cells(stats, f"{base_name} from {n_cells} cell{'s' if n_cells != 1 else ''}",
                                        species=species, ploidy=s.karyotype.ploidy if s.karyotype else 2,
                                        seen_at_least=need, colours=colours)
            self._open_karyotype(edit=False, prefilled=k, replaces="")

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

    def _save_measurements(self) -> None:
        s = self._session
        if s is None or self._folder is None:
            return
        stats = pr.stats_over_cells(self._scope_details(), [c.name for c in s.karyotype.chromosomes] if s.karyotype else None)
        path, _ = QFileDialog.getSaveFileName(self, "Save measurements", str(self._folder / "chromosome_measurements.tsv"), "TSV (*.tsv)")
        if not path:
            return
        lines = ["chromosome\tcells\thomologues\tlength_um\tlength_px\trelative_length_pct\tcentromere_index"]
        for st in stats:
            um = f"{st.length_um:.3f}" if st.length_um is not None else ""
            lines.append(f"{st.name}\t{st.cells}\t{st.homologues}\t{um}\t{st.length_px:.1f}\t{st.rel_length:.2f}\t{st.ci:.2f}")
        try:
            Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")
        except OSError as exc:
            QMessageBox.warning(self, "Save measurements", str(exc))

    # ------------------------------------------------------------------ cells

    def _step_cell(self, step: int) -> None:
        s = self._session
        if s is None:
            return
        ids = [c.image_id for c in s.slide_cells()]
        if s.image.id in ids:
            k = ids.index(s.image.id) + step
            if 0 <= k < len(ids):
                self.show_image(ids[k])

    def _sync_cells(self) -> None:
        s = self._session
        if s is None:
            self._nav.setText("")
            return
        cells = s.slide_cells()
        ids = [c.image_id for c in cells]
        pos = ids.index(s.image.id) + 1 if s.image.id in ids else 0
        self._nav.setText(f"cell {cells[pos - 1].label} of {len(cells)}" if pos else "")
        self._prev.setEnabled(pos > 1)
        self._next.setEnabled(0 < pos < len(cells))

    # ------------------------------------------------------------------ saving

    def flush(self) -> None:
        self._save_timer.stop()
        s = self._session
        if s is None or self._folder is None or not s.dirty or self._busy:
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
