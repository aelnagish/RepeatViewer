"""Figure tab of the Cytogenetics workspace: panels, letters, labels, arrowheads, scale bar, legend."""

from __future__ import annotations

import math
from pathlib import Path

from PySide6.QtCore import QPointF, QRectF, QSettings, QSizeF, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QFont, QGuiApplication, QImage, QPageSize, QPainter, QPainterPath, QPen, QTextDocument
from PySide6.QtWidgets import (
    QCheckBox,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QGraphicsDropShadowEffect,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSlider,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import figure as fg
from repeatviewer.backends.cytogenetics import idiogram as idg
from repeatviewer.backends.cytogenetics.model import image_details, load_slides
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cytogenetics_image import rgb_to_qimage
from repeatviewer.widgets.cytogenetics_style import clear_layout, link_button
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.sequence_annotation_rail import _CollapsibleSection

_SETTINGS = "cytogenetics_workspace"
SIDE_W = 330
PREVIEW_W = 680  # page width on screen (as in the mockup)
PAGE_PAD = 16
SAVE_DELAY_MS = 400
NEW_FIGURE = "New figure…"
IDIOGRAM_LOGICAL_W = 760.0

HELP_TEXT = (
    "<b>Figure</b><br>"
    "Panels of one or more cells (each channel and the merge) and the idiogram, with letters, channel labels, "
    "arrowheads, a scale bar and a legend written from the image details.<br>"
    "<b>On the page</b>: click a panel to select it; with Arrowheads on, a click places an arrowhead "
    "(right-click removes it, the mouse wheel over it turns it). Drag across a panel to crop it; double-click to show the whole image again.<br>"
    "<b>Same crop for all panels</b> keeps the channels of one cell aligned.<br>"
    "<b>Legend</b>: written from probes, dyes, stage, species, karyotype and the adjustment log. Edit it on the page; "
    "<i>Write again</i> replaces your edits.<br>"
    "<b>Save</b> writes the panels (without the legend) as SVG, PDF or TIFF at 600 dpi; copy the legend and the methods "
    "paragraph into your manuscript."
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


# ---------------------------------------------------------------------------
# Painting (shared by the page preview and every export)
# ---------------------------------------------------------------------------


class FigurePainter:
    """Holds the panel pictures of one figure and paints them at any scale."""

    def __init__(self) -> None:
        self.folder: Path | None = None
        self.fig: fg.Figure | None = None
        self.images: dict[str, ci.CytoImage] = {}
        self.pictures: dict[str, QImage] = {}  # panel id -> cropped picture
        self.full_width: dict[str, int] = {}  # panel id -> width of the uncropped source
        self.errors: dict[str, str] = {}
        self.idiogram: idg.IdiogramView | None = None
        self._canvas = None
        self._cache: dict = {}
        self._key: dict[str, tuple] = {}

    def load(self, folder: Path, fig: fg.Figure) -> None:
        self.folder, self.fig = folder, fig
        self.images = {i.id: i for i in ci.load_images(folder)}
        alive = set()
        for p in fig.panels:
            alive.add(p.id)
            if p.is_idiogram:
                continue
            image = self.images.get(p.image_id)
            crop = fig.crop_of(p)
            stamp = self._stamp(image)
            key = (p.image_id, p.view, crop, stamp)
            if self._key.get(p.id) == key and p.id in self.pictures:
                continue
            self._key[p.id] = key
            self.errors.pop(p.id, None)
            if image is None:
                self.pictures.pop(p.id, None)
                self.errors[p.id] = "image removed"
                continue
            try:
                if image_details(image).channels:
                    rgb = fg.render_panel(folder, image, p, (0.0, 0.0, 1.0, 1.0), self._cache)
                    self.full_width[p.id] = rgb.shape[1]
                    self.pictures[p.id] = rgb_to_qimage(fg.crop_array(rgb, crop))
                else:
                    pic = QImage(str(ci.image_file_path(folder, image)))
                    if pic.isNull():
                        raise ValueError("picture file missing")
                    self.full_width[p.id] = pic.width()
                    x0, y0, x1, y1 = crop
                    self.pictures[p.id] = pic.copy(int(x0 * pic.width()), int(y0 * pic.height()),
                                                   max(1, int((x1 - x0) * pic.width())), max(1, int((y1 - y0) * pic.height())))
            except Exception as exc:  # noqa: BLE001 - missing TIFF, removed channel
                self.pictures.pop(p.id, None)
                self.errors[p.id] = str(exc)
        for gone in set(self.pictures) - alive:
            self.pictures.pop(gone, None)
        if any(p.is_idiogram for p in fig.panels):
            self.idiogram = idg.compute_view(folder, cache=self._cache)

    def forget(self) -> None:
        """Drop cached pictures (after the images changed)."""
        self._key.clear()
        self._cache.clear()

    @staticmethod
    def _stamp(image: ci.CytoImage | None) -> tuple:
        if image is None:
            return ()
        d = image_details(image)
        return tuple((c.id, c.display.color, c.display.low, c.display.high, c.display.gamma, c.display.visible) for c in d.channels)

    def aspects(self) -> dict[str, float]:
        out = {}
        for pid, pic in self.pictures.items():
            out[pid] = pic.height() / max(1, pic.width())
        for p in self.fig.panels if self.fig else []:
            out.setdefault(p.id, 0.75)
        return out

    def idiogram_aspect(self) -> float:
        n = len(self.idiogram.rows) if self.idiogram else 8
        return 0.42 if n <= 12 else 0.32

    def layout(self) -> tuple[list[tuple[fg.Panel, fg.Box]], float]:
        return fg.layout(self.fig, self.aspects(), self.idiogram_aspect())

    def _idiogram_canvas(self):
        if self._canvas is None:
            from repeatviewer.widgets.cytogenetics_idiogram import IdiogramCanvas

            self._canvas = IdiogramCanvas()
        return self._canvas

    def paint(self, p: QPainter, px_per_mm: float, *, selected: str = "") -> None:
        fig = self.fig
        if fig is None:
            return
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.setRenderHint(QPainter.RenderHint.SmoothPixmapTransform, True)
        boxes, _h = self.layout()
        last_image = next((pp.id for pp, _b in reversed(boxes) if not pp.is_idiogram), "")
        pt = 0.3528 * px_per_mm  # one point in device pixels
        for index, (panel, box) in enumerate(boxes):
            r = QRectF(box.x * px_per_mm, box.y * px_per_mm, box.w * px_per_mm, box.h * px_per_mm)
            if panel.is_idiogram:
                self._paint_idiogram(p, r)
            else:
                self._paint_image_panel(p, r, panel, pt, show_bar=fig.scale_bar and panel.id == last_image)
            if fig.letters:
                font = QFont(fig.font_family or "Arial")
                font.setPixelSize(max(6, int(round(fig.font_pt * 1.5 * pt))))
                font.setBold(True)
                p.setFont(font)
                p.setPen(QColor("#1C1C1E") if panel.is_idiogram else QColor("#FFFFFF"))
                mm = pt / 0.3528
                p.drawText(r.adjusted(1.6 * mm, 1.0 * mm, 0, 0), Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop,
                           fg.letter(index))
            if selected == panel.id:
                p.setPen(QPen(QColor("#1F4E9A"), 2, Qt.PenStyle.DashLine))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRect(r.adjusted(1, 1, -1, -1))

    def _paint_idiogram(self, p: QPainter, r: QRectF) -> None:
        view = self.idiogram
        p.fillRect(r, QColor("#FFFFFF"))
        if view is None or not view.rows:
            p.setPen(QColor("#5A5D63"))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter, "No idiogram yet (make a karyotype and mark signals)")
            return
        if view.settings.figure_style == "paper":
            from repeatviewer.widgets.cytogenetics_idiogram import paint_paper, paper_height

            w = idg.PAPER_WIDTH_PT
            h = paper_height(view.karyotype, len(idg.paper_table(view.karyotype, view.pooled)), view.settings.paper)
            scale = min(r.width() / w, r.height() / h)
            p.save()
            p.translate(r.left() + (r.width() - w * scale) / 2, r.top() + (r.height() - h * scale) / 2)
            p.scale(scale, scale)
            paint_paper(p, w, h, view.karyotype, view.pooled, view.colours, view.settings.paper)
            p.restore()
            return
        canvas = self._idiogram_canvas()
        logical_h = IDIOGRAM_LOGICAL_W * r.height() / max(r.width(), 1.0)
        canvas.resize(int(IDIOGRAM_LOGICAL_W), int(max(logical_h, 200)))
        canvas.set_rows(view.rows, view.heights, view.colours)
        canvas.resize(int(IDIOGRAM_LOGICAL_W), int(max(logical_h, 200)))
        p.save()
        p.translate(r.topLeft())
        p.scale(r.width() / canvas.width(), r.height() / canvas.height())
        canvas.paint_on(p, canvas.width(), canvas.height())
        p.restore()

    def _paint_image_panel(self, p: QPainter, r: QRectF, panel: fg.Panel, pt: float, *, show_bar: bool) -> None:
        fig = self.fig
        mm = pt / 0.3528  # one millimetre in device pixels
        p.fillRect(r, QColor("#030405"))
        pic = self.pictures.get(panel.id)
        image = self.images.get(panel.image_id)
        details = image_details(image) if image else None
        if pic is None:
            p.setPen(QColor("#FF9DAE"))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter | Qt.TextFlag.TextWordWrap, self.errors.get(panel.id, "missing"))
            return
        p.drawImage(r, pic)
        crop = fig.crop_of(panel)
        cw, chh = crop[2] - crop[0], crop[3] - crop[1]

        def to_panel(fx: float, fy: float) -> QPointF | None:
            u, v = (fx - crop[0]) / cw, (fy - crop[1]) / chh
            if not (0 <= u <= 1 and 0 <= v <= 1):
                return None
            return QPointF(r.left() + u * r.width(), r.top() + v * r.height())

        if fig.chromosome_numbers and details is not None and pic is not None and image is not None:
            full_w = self.full_width.get(panel.id, 0)
            full_h = full_w * (pic.height() / max(pic.width(), 1)) * cw / chh if full_w else 0
            font = QFont(fig.font_family or "Arial")
            font.setPixelSize(max(5, int(round(fig.font_pt * 0.9 * pt))))
            font.setBold(True)
            p.setFont(font)
            for m in details.chromosomes:
                if not full_w or not full_h:
                    break
                q = to_panel(m.x / full_w, m.y / full_h)
                if q is not None:
                    p.setPen(QColor("#000000"))
                    p.drawText(QRectF(q.x() - 20 * pt + pt * 0.4, q.y() - 4 * pt + pt * 0.4, 40 * pt, 4 * pt),
                               Qt.AlignmentFlag.AlignCenter, m.label)
                    p.setPen(QColor("#FFE08A"))
                    p.drawText(QRectF(q.x() - 20 * pt, q.y() - 4 * pt, 40 * pt, 4 * pt), Qt.AlignmentFlag.AlignCenter, m.label)
        if fig.channel_labels:
            label = fg.auto_label(panel, details)
            if label:
                font = QFont(fig.font_family or "Arial")
                font.setPixelSize(max(5, int(round(fig.font_pt * pt))))
                font.setBold(True)
                p.setFont(font)
                p.setPen(QColor(fg.label_colour(panel, details)))
                p.drawText(r.adjusted(0, 1.2 * mm, -1.6 * mm, 0), Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignTop, label)
        if fig.arrowheads:
            size = 2.8 * mm
            for a in panel.arrows:
                q = to_panel(a.x, a.y)
                if q is not None:
                    paint_arrow(p, q, a.angle, size)
        if show_bar and details is not None:
            frac = fg.scale_bar_fraction(details, crop, self.full_width.get(panel.id, 0), fig.scale_um)
            if frac and frac < 0.9:
                length = frac * r.width()
                margin = 2.2 * mm
                bar_h = max(1.0, 0.8 * mm)
                p.fillRect(QRectF(r.right() - margin - length, r.bottom() - margin - bar_h, length, bar_h), QColor("#FFFFFF"))


def paint_arrow(p: QPainter, tip: QPointF, angle_deg: float, size: float) -> None:
    """Filled arrowhead with its tip at ``tip``, pointing in ``angle_deg`` (screen degrees)."""
    a = math.radians(angle_deg)
    dx, dy = math.cos(a), math.sin(a)
    back = QPointF(tip.x() - dx * size, tip.y() - dy * size)
    nx, ny = -dy, dx
    half = size * 0.42
    notch = QPointF(tip.x() - dx * size * 0.72, tip.y() - dy * size * 0.72)
    path = QPainterPath(tip)
    path.lineTo(back.x() + nx * half, back.y() + ny * half)
    path.lineTo(notch)
    path.lineTo(back.x() - nx * half, back.y() - ny * half)
    path.closeSubpath()
    p.setPen(QPen(QColor(0, 0, 0, 160), max(0.6, size * 0.05)))
    p.setBrush(QColor("#FFFFFF"))
    p.drawPath(path)


# ---------------------------------------------------------------------------
# The page on screen
# ---------------------------------------------------------------------------


class FigurePage(QWidget):
    """White page: panels painted by FigurePainter, the legend editable below them."""

    panelClicked = Signal(str)
    changed = Signal()  # arrows or crop changed
    legendEdited = Signal()

    def __init__(self, painter: FigurePainter, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.fp = painter
        self.selected = ""
        self.setFixedWidth(PREVIEW_W)
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.setStyleSheet("FigurePage { background: #FFFFFF; }")
        self.legend = QTextEdit(self)
        self.legend.setFrameShape(QFrame.Shape.NoFrame)
        self.legend.setStyleSheet("QTextEdit { border: none; background: transparent; color: #1C1C1E; }")
        self.legend.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.legend.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self.legend.textChanged.connect(self._fit_legend)
        self.legend.textChanged.connect(self.legendEdited)
        shadow = QGraphicsDropShadowEffect(self)
        shadow.setBlurRadius(20)
        shadow.setOffset(0, 2)
        shadow.setColor(QColor(0, 0, 0, 40))
        self.setGraphicsEffect(shadow)
        self._press: QPointF | None = None
        self._drag_rect: QRectF | None = None
        self._panels_h = 0.0

    def px_per_mm(self) -> float:
        fig = self.fp.fig
        return (PREVIEW_W - 2 * PAGE_PAD) / (fig.page_width if fig else 174.0)

    def relayout(self) -> None:
        fig = self.fp.fig
        if fig is None:
            return
        _boxes, h_mm = self.fp.layout()
        self._panels_h = h_mm * self.px_per_mm()
        font = QFont(fig.font_family or "Arial")
        font.setPixelSize(max(9, int(round(fig.font_pt * 0.3528 * self.px_per_mm() * 1.25))))
        self.legend.document().setDefaultFont(font)
        self._fit_legend()
        self.update()

    def _fit_legend(self) -> None:
        width = PREVIEW_W - 2 * PAGE_PAD
        self.legend.document().setTextWidth(width)
        h = int(self.legend.document().size().height()) + 6
        top = int(PAGE_PAD + self._panels_h + 12)
        self.legend.setGeometry(PAGE_PAD - 4, top, width + 8, h)
        self.setFixedHeight(top + h + PAGE_PAD)

    def set_legend_html(self, html: str) -> None:
        self.legend.blockSignals(True)
        self.legend.setHtml(html)
        self.legend.blockSignals(False)
        self._fit_legend()

    def paintEvent(self, event) -> None:  # noqa: N802
        p = QPainter(self)
        p.fillRect(self.rect(), QColor("#FFFFFF"))
        if self.fp.fig is None or not self.fp.fig.panels:
            p.setPen(QColor("#5A5D63"))
            p.drawText(QRectF(0, 0, self.width(), max(self._panels_h, 120) + PAGE_PAD * 2), Qt.AlignmentFlag.AlignCenter,
                       "No panels yet. Choose Add to figure on an image in the Library,\nor + Add under Panels.")
            return
        p.translate(PAGE_PAD, PAGE_PAD)
        self.fp.paint(p, self.px_per_mm(), selected=self.selected)
        if self._drag_rect is not None:
            p.setPen(QPen(QColor("#FFD35C"), 1.5, Qt.PenStyle.DashLine))
            p.setBrush(QColor(255, 211, 92, 40))
            p.drawRect(self._drag_rect)

    # ------------------------------------------------------------------ mouse

    def _panel_at(self, pos: QPointF) -> tuple[fg.Panel, QRectF] | None:
        if self.fp.fig is None:
            return None
        k = self.px_per_mm()
        local = QPointF(pos.x() - PAGE_PAD, pos.y() - PAGE_PAD)
        for panel, box in self.fp.layout()[0]:
            r = QRectF(box.x * k, box.y * k, box.w * k, box.h * k)
            if r.contains(local):
                return panel, r
        return None

    def _source_fraction(self, panel: fg.Panel, r: QRectF, pos: QPointF) -> tuple[float, float]:
        crop = self.fp.fig.crop_of(panel)
        u = (pos.x() - PAGE_PAD - r.left()) / max(r.width(), 1.0)
        v = (pos.y() - PAGE_PAD - r.top()) / max(r.height(), 1.0)
        u, v = max(0.0, min(1.0, u)), max(0.0, min(1.0, v))
        return crop[0] + u * (crop[2] - crop[0]), crop[1] + v * (crop[3] - crop[1])

    def _arrow_near(self, panel: fg.Panel, fx: float, fy: float) -> fg.Arrow | None:
        crop = self.fp.fig.crop_of(panel)
        tol = 0.04 * (crop[2] - crop[0])
        return min((a for a in panel.arrows if math.hypot(a.x - fx, a.y - fy) <= tol),
                   key=lambda a: math.hypot(a.x - fx, a.y - fy), default=None)

    def mousePressEvent(self, event) -> None:  # noqa: N802
        hit = self._panel_at(event.position())
        if hit is None:
            return
        panel, r = hit
        self.selected = panel.id
        self.panelClicked.emit(panel.id)
        if panel.is_idiogram:
            self.update()
            return
        if event.button() == Qt.MouseButton.RightButton:
            fx, fy = self._source_fraction(panel, r, event.position())
            arrow = self._arrow_near(panel, fx, fy)
            if arrow is not None:
                panel.arrows.remove(arrow)
                self.changed.emit()
            self.update()
            return
        self._press = event.position()
        self.update()

    def mouseMoveEvent(self, event) -> None:  # noqa: N802
        if self._press is None:
            return
        a, b = self._press, event.position()
        if math.hypot(b.x() - a.x(), b.y() - a.y()) > 6:
            self._drag_rect = QRectF(QPointF(min(a.x(), b.x()) - PAGE_PAD, min(a.y(), b.y()) - PAGE_PAD),
                                     QPointF(max(a.x(), b.x()) - PAGE_PAD, max(a.y(), b.y()) - PAGE_PAD))
            self.update()

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802
        if self._press is None:
            return
        start, self._press = self._press, None
        hit = self._panel_at(start)
        self._drag_rect, dragged = None, self._drag_rect is not None
        if hit is None:
            self.update()
            return
        panel, r = hit
        fig = self.fp.fig
        if dragged:
            x0, y0 = self._source_fraction(panel, r, start)
            x1, y1 = self._source_fraction(panel, r, event.position())
            crop = fg.normalise_crop(x0, y0, x1, y1)
            if crop is not None:
                if fig.same_crop:
                    fig.crop = crop
                else:
                    panel.crop = crop
                self.changed.emit()
        elif fig.arrowheads:
            fx, fy = self._source_fraction(panel, r, start)
            panel.arrows.append(fg.Arrow(fx, fy))
            self.changed.emit()
        self.update()

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        hit = self._panel_at(event.position())
        if hit is None or hit[0].is_idiogram:
            return
        panel, _r = hit
        fig = self.fp.fig
        # the first click of the double-click placed an arrowhead: take it back
        if fig.arrowheads and panel.arrows:
            panel.arrows.pop()
        if fig.same_crop:
            fig.crop = None
        else:
            panel.crop = None
        self.changed.emit()

    def wheelEvent(self, event) -> None:  # noqa: N802
        hit = self._panel_at(event.position())
        if hit is None or hit[0].is_idiogram:
            event.ignore()
            return
        panel, r = hit
        fx, fy = self._source_fraction(panel, r, event.position())
        arrow = self._arrow_near(panel, fx, fy)
        if arrow is None:
            event.ignore()
            return
        arrow.angle = (arrow.angle + (15 if event.angleDelta().y() > 0 else -15)) % 360
        self.changed.emit()
        self.update()
        event.accept()


# ---------------------------------------------------------------------------
# The tab
# ---------------------------------------------------------------------------


class FigureTab(QWidget):
    name = "Figure"
    crumbChanged = Signal(str)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._folder: Path | None = None
        self._fig: fg.Figure | None = None
        self._old_name = ""
        self._selected = ""
        self._busy = False
        self._dirty = True
        self._painter = FigurePainter()

        # ---- page area
        self._page = FigurePage(self._painter, self)
        self._page.panelClicked.connect(self._select_panel)
        self._page.changed.connect(self._on_page_changed)
        self._page.legendEdited.connect(self._on_legend_edited)
        holder = QWidget(self)
        holder.setObjectName("cytoPage")
        holder.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        hl = QVBoxLayout(holder)
        hl.setContentsMargins(28, 28, 28, 28)
        hl.addWidget(self._page, 0, Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
        hl.addStretch(1)
        page_scroll = QScrollArea(self)
        page_scroll.setWidgetResizable(True)
        page_scroll.setFrameShape(QFrame.Shape.NoFrame)
        page_scroll.setWidget(holder)

        # ---- side panel
        body = QWidget(self)
        body.setObjectName("cytoAside")
        body.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        side = QVBoxLayout(body)
        side.setContentsMargins(14, 12, 14, 12)
        side.setSpacing(14)

        page_sec = _CollapsibleSection("Page", body, expanded=True)
        self._which = NoWheelComboBox(self)
        self._which.activated.connect(self._on_which)
        self._name = QLineEdit(self)
        self._name.setPlaceholderText("Figure name")
        self._name.editingFinished.connect(self._on_name)
        self._width = NoWheelComboBox(self)
        for text, mm in fg.WIDTHS:
            self._width.addItem(text, mm)
        self._width.activated.connect(self._on_page_setting)
        self._custom = QDoubleSpinBox(self)
        self._custom.setRange(30, 400)
        self._custom.setSuffix(" mm")
        self._custom.setDecimals(0)
        self._custom.editingFinished.connect(self._on_page_setting)
        self._layout = NoWheelComboBox(self)
        for text, key in fg.LAYOUTS:
            self._layout.addItem(text, key)
        self._layout.activated.connect(self._on_page_setting)
        self._spacing = QSlider(Qt.Orientation.Horizontal, self)
        self._spacing.setRange(0, 10)
        from repeatviewer.widgets.cytogenetics_style import app_slider

        app_slider(self._spacing)
        self._spacing.sliderReleased.connect(self._on_page_setting)
        self._spacing.valueChanged.connect(lambda _v: None if self._spacing.isSliderDown() else self._on_page_setting())
        from PySide6.QtWidgets import QFontComboBox, QSpinBox

        self._font = QFontComboBox(self)
        self._font.setMinimumWidth(60)
        self._font.setSizeAdjustPolicy(QFontComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self._font.currentFontChanged.connect(self._on_page_setting)
        self._font_size = QDoubleSpinBox(self)
        self._font_size.setRange(5, 16)
        self._font_size.setDecimals(1)
        self._font_size.setSingleStep(0.5)
        self._font_size.setSuffix(" pt")
        self._font_size.setFixedWidth(84)
        self._font_size.editingFinished.connect(self._on_page_setting)
        self._rows = QSpinBox(self)
        self._cols = QSpinBox(self)
        for box in (self._rows, self._cols):
            box.setRange(1, 12)
            box.editingFinished.connect(self._on_page_setting)
        self._grid_host = QWidget(self)
        gl = QHBoxLayout(self._grid_host)
        gl.setContentsMargins(0, 0, 0, 0)
        gl.addWidget(self._rows, 1)
        gl.addWidget(QLabel("rows ×", self._grid_host))
        gl.addWidget(self._cols, 1)
        gl.addWidget(QLabel("columns", self._grid_host))
        self._custom.setFixedWidth(84)
        for text, w in (("Figure", self._which), ("Name", self._name), ("Width", self._width),
                        ("Layout", self._layout), ("", self._grid_host), ("Spacing", self._spacing), ("Font", self._font)):
            r = QHBoxLayout()
            lab = QLabel(text, self)
            lab.setFixedWidth(90)
            r.addWidget(lab)
            r.addWidget(w, 1)
            if w is self._width:
                r.addWidget(self._custom)
            if w is self._font:
                r.addWidget(self._font_size)
            if w is self._grid_host:
                self._grid_row = (lab, w)
            page_sec.body_layout().addLayout(r)
        side.addWidget(page_sec)

        panels_sec = _CollapsibleSection("Panels", body, expanded=True)
        add = QPushButton("+ Add", self)
        add_menu = QMenu(add)
        add_menu.addAction("Channels and merge of the current cell", lambda: self._add_current(all_channels=True))
        add_menu.addAction("Merge of the current cell", lambda: self._add_current(all_channels=False))
        add_menu.addAction("Idiogram", self.add_idiogram)
        add.setMenu(add_menu)
        head = QHBoxLayout()
        head.addStretch(1)
        head.addWidget(add)
        panels_sec.body_layout().addLayout(head)
        self._panel_host = QWidget(panels_sec)
        self._panel_lay = QVBoxLayout(self._panel_host)
        self._panel_lay.setContentsMargins(0, 0, 0, 0)
        self._panel_lay.setSpacing(6)
        panels_sec.add_widget(self._panel_host)
        self._edit_host = QWidget(panels_sec)
        self._edit_lay = QVBoxLayout(self._edit_host)
        self._edit_lay.setContentsMargins(0, 4, 0, 0)
        self._edit_lay.setSpacing(6)
        panels_sec.add_widget(self._edit_host)
        self._same_crop = QCheckBox("Same crop for all panels", self)
        self._same_crop.toggled.connect(self._on_marks)
        panels_sec.add_widget(self._same_crop)
        hint = QLabel("Drag across a panel to crop it, double-click to see the whole image.", self)
        hint.setObjectName("cytoMuted")
        hint.setWordWrap(True)
        panels_sec.add_widget(hint)
        side.addWidget(panels_sec)

        marks_sec = _CollapsibleSection("Marks", body, expanded=True)
        self._letters = QCheckBox("Panel letters", self)
        self._labels = QCheckBox("Channel labels (colour of the dye)", self)
        self._bar = QCheckBox("Scale bar on last panel", self)
        self._bar.setMinimumWidth(0)
        self._bar_len = QDoubleSpinBox(self)
        self._bar_len.setRange(0.5, 500)
        self._bar_len.setDecimals(1)
        self._bar_len.setSuffix(" µm")
        self._bar_len.setFixedWidth(86)
        self._bar_len.editingFinished.connect(self._on_marks)
        self._arrows = QCheckBox("Arrowheads (click to place)", self)
        self._numbers = QCheckBox("Chromosome numbers", self)
        for box in (self._letters, self._labels, self._bar, self._arrows, self._numbers):
            box.toggled.connect(self._on_marks)
        marks_sec.add_widget(self._letters)
        marks_sec.add_widget(self._labels)
        br = QHBoxLayout()
        br.addWidget(self._bar, 1)
        br.addWidget(self._bar_len)
        marks_sec.body_layout().addLayout(br)
        marks_sec.add_widget(self._arrows)
        marks_sec.add_widget(self._numbers)
        side.addWidget(marks_sec)

        legend_sec = _CollapsibleSection("Legend text", body, expanded=True)
        note = QLabel("Written from probes, dyes, stage and the adjustment log. Edit it on the page.", self)
        note.setObjectName("cytoMuted")
        note.setWordWrap(True)
        legend_sec.add_widget(note)
        self._number = QLineEdit(self)
        self._number.setPlaceholderText("Figure 1.")
        self._number.editingFinished.connect(self._on_number)
        nr = QHBoxLayout()
        lab = QLabel("Starts with", self)
        lab.setFixedWidth(90)
        nr.addWidget(lab)
        nr.addWidget(self._number, 1)
        legend_sec.body_layout().addLayout(nr)
        lr = QHBoxLayout()
        copy_methods = QPushButton("Copy methods paragraph", self)
        copy_methods.clicked.connect(self._copy_methods)
        copy_legend = QPushButton("Copy legend", self)
        copy_legend.clicked.connect(self._copy_legend)
        lr.addWidget(copy_methods, 1)
        lr.addWidget(copy_legend, 1)
        legend_sec.body_layout().addLayout(lr)
        again = link_button("Write the legend again (drops your edits)", self)
        again.clicked.connect(self._rewrite_legend)
        legend_sec.add_widget(again)
        side.addWidget(legend_sec)
        side.addStretch(1)

        save = QPushButton("Save ▾", self)
        save_menu = QMenu(save)
        save_menu.addAction("SVG…", lambda: self._export("svg"))
        save_menu.addAction("PDF…", lambda: self._export("pdf"))
        save_menu.addAction("TIFF 600 dpi…", lambda: self._export("tiff"))
        save.setMenu(save_menu)
        fmt = QLabel("SVG · PDF · TIFF 600 dpi", self)
        fmt.setObjectName("cytoMuted")
        bottom = QHBoxLayout()
        bottom.addWidget(save)
        bottom.addWidget(fmt)
        bottom.addStretch(1)
        bottom.addWidget(make_help_button(self, "Figure", HELP_TEXT))
        side.addLayout(bottom)

        aside = QScrollArea(self)
        aside.setWidgetResizable(True)
        aside.setFrameShape(QFrame.Shape.NoFrame)
        aside.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        aside.setWidget(body)
        aside.setMinimumWidth(240)
        aside.setMaximumWidth(SIDE_W)

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addWidget(page_scroll, 1)
        line = QFrame(self)
        line.setFrameShape(QFrame.Shape.VLine)
        line.setStyleSheet("color: #E0E0E4;")
        root.addWidget(line)
        root.addWidget(aside)

        self._save_timer = QTimer(self)
        self._save_timer.setSingleShot(True)
        self._save_timer.setInterval(SAVE_DELAY_MS)
        self._save_timer.timeout.connect(self.flush)
        self.current_image = ""  # set by the workspace (cell shown in Image / Signals)

    # ------------------------------------------------------------------ loading

    def refresh(self, folder: Path | None) -> None:
        self._folder = Path(folder) if folder else None
        self._painter.forget()
        self._dirty = True
        if self.isVisible():
            self._load()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        if self._dirty:
            self._load()
        self.crumbChanged.emit(self.crumb())

    def hideEvent(self, event) -> None:  # noqa: N802
        self.flush()
        super().hideEvent(event)

    def crumb(self) -> str:
        if self._folder is None:
            return ""
        name = fg.figure_path(self._folder, self._fig.name).name if self._fig else ""
        return f"{self._folder.name} › {name}" if name else self._folder.name

    def _load(self) -> None:
        self._dirty = False
        if self._folder is None:
            return
        figures = fg.list_figures(self._folder)
        if self._fig is None or not any(f.name == self._fig.name for f in figures):
            wanted = _read("figure_name", "")
            self._fig = next((f for f in figures if f.name == wanted), None) or (figures[0] if figures else None)
            if self._fig is None:
                self._fig = fg.Figure(name=fg.new_figure_name(self._folder))
            self._old_name = self._fig.name if any(f.name == self._fig.name for f in figures) else ""
        self._redraw(controls=True)

    def _ensure(self) -> fg.Figure:
        if self._fig is None:
            self._load()
        if self._fig is None:  # no folder
            self._fig = fg.Figure()
        return self._fig

    def _redraw(self, *, controls: bool = False, legend: bool = True) -> None:
        fig = self._fig
        if fig is None or self._folder is None:
            return
        self._painter.load(self._folder, fig)
        self._page.selected = self._selected
        self._page.relayout()
        if legend:
            self._page.set_legend_html(fig.legend or fg.legend_html(self._folder, fig))
        if controls:
            self._fill_controls()
        self._fill_panels()
        self._page.update()

    # ------------------------------------------------------------------ controls

    def _fill_controls(self) -> None:
        fig = self._fig
        self._busy = True
        try:
            self._which.clear()
            names = [f.name for f in fg.list_figures(self._folder)]
            if fig.name not in names:
                names.append(fig.name)
            for n in names:
                self._which.addItem(n, n)
            self._which.addItem(NEW_FIGURE, "\x00new")
            self._which.setCurrentIndex(max(0, self._which.findData(fig.name)))
            self._name.setText(fig.name)
            idx = next((i for i, (_t, mm) in enumerate(fg.WIDTHS) if mm == fig.width_mm), len(fg.WIDTHS) - 1)
            self._width.setCurrentIndex(idx)
            self._custom.setValue(fig.custom_mm)
            self._custom.setVisible(fig.width_mm <= 0)
            self._layout.setCurrentIndex(max(0, self._layout.findData(fig.layout)))
            self._spacing.setValue(fig.spacing)
            self._font.blockSignals(True)
            self._font.setCurrentFont(QFont(fig.font_family or "Arial"))
            self._font.blockSignals(False)
            self._font_size.setValue(fig.font_pt)
            self._rows.setValue(fig.rows)
            self._cols.setValue(fig.cols)
            for w in self._grid_row:
                w.setVisible(fig.layout == "custom")
            for box, value in ((self._same_crop, fig.same_crop), (self._letters, fig.letters), (self._labels, fig.channel_labels),
                               (self._bar, fig.scale_bar), (self._arrows, fig.arrowheads), (self._numbers, fig.chromosome_numbers)):
                box.setChecked(value)
            self._bar_len.setValue(fig.scale_um)
            self._number.setText(fig.number)
        finally:
            self._busy = False

    def _clear(self, layout) -> None:
        clear_layout(layout)

    def _fill_panels(self) -> None:
        self._clear(self._panel_lay)
        fig = self._fig
        slides = {s.id: s.name for s in load_slides(self._folder)} if self._folder else {}
        if fig.panels:
            from PySide6.QtWidgets import QListWidget, QListWidgetItem

            lst = QListWidget(self._panel_host)
            lst.setTextElideMode(Qt.TextElideMode.ElideRight)
            lst.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
            for i, panel in enumerate(fig.panels):
                image = self._painter.images.get(panel.image_id)
                slide = slides.get(image_details(image).slide_id, "") if image else ""
                text = fg.source_text(panel, image, slide)
                item = QListWidgetItem(f"{fg.letter(i)}    {text}")
                item.setToolTip(text)
                item.setData(Qt.ItemDataRole.UserRole, panel.id)
                lst.addItem(item)
                if panel.id == self._selected:
                    lst.setCurrentItem(item)
            lst.setFixedHeight(min(8, lst.count()) * max(lst.sizeHintForRow(0), 16) + 2 * lst.frameWidth() + 2)
            lst.itemClicked.connect(lambda it: self._select_panel(it.data(Qt.ItemDataRole.UserRole)))
            self._panel_lay.addWidget(lst)
        if not fig.panels:
            empty = QLabel("No panels yet. Use Add to figure in the Library or + Add.", self._panel_host)
            empty.setObjectName("cytoMuted")
            empty.setWordWrap(True)
            self._panel_lay.addWidget(empty)
        self._fill_editor()

    def _fill_editor(self) -> None:
        self._clear(self._edit_lay)
        fig = self._fig
        panel = fig.panel(self._selected) if fig else None
        if panel is None:
            return
        image = self._painter.images.get(panel.image_id)
        if not panel.is_idiogram:
            show = NoWheelComboBox(self._edit_host)
            details = image_details(image) if image else None
            if details is not None and details.channels:
                show.addItem("Merge", fg.VIEW_MERGE)
                for c in fg.ordered_channels(details):
                    show.addItem(f"{c.dye or 'channel'}" + (f" · {c.probe}" if c.probe else ""), c.id)
            else:
                show.addItem("Picture", fg.VIEW_MERGE)
            show.setCurrentIndex(max(0, show.findData(panel.view)))
            show.activated.connect(lambda _i, p=panel, cb=show: self._set_panel(p, view=cb.currentData()))
            label = QLineEdit(panel.label, self._edit_host)
            label.setPlaceholderText(fg.auto_label(fg.Panel(id="x", view=panel.view), details) or "no label")
            label.editingFinished.connect(lambda p=panel, w=label: self._set_panel(p, label=w.text().strip()))
            for text, w in (("Show", show), ("Label", label)):
                r = QHBoxLayout()
                lab = QLabel(text, self._edit_host)
                lab.setFixedWidth(60)
                r.addWidget(lab)
                r.addWidget(w, 1)
                self._edit_lay.addLayout(r)
        r = QHBoxLayout()
        for text, slot in (("Move up", lambda: self._move(-1)), ("Move down", lambda: self._move(1)), ("Remove", self._remove)):
            b = QPushButton(text, self._edit_host)
            b.clicked.connect(slot)
            r.addWidget(b, 1)
        self._edit_lay.addLayout(r)
        if panel.arrows:
            b = QPushButton("Clear arrowheads", self._edit_host)
            b.clicked.connect(lambda: self._set_panel(panel, arrows=[]))
            r2 = QHBoxLayout()
            r2.addWidget(b)
            r2.addStretch(1)
            self._edit_lay.addLayout(r2)

    # ------------------------------------------------------------------ actions

    def add_cell(self, image_id: str, *, all_channels: bool = True) -> None:
        if self._folder is None:
            return
        fig = self._ensure()
        image = next((i for i in ci.load_images(self._folder) if i.id == image_id), None)
        if image is None:
            return
        if all_channels:
            new = fg.add_cell(fig, image)
        else:
            new = [fg.Panel(id=fg.new_id(), image_id=image.id)]
            fig.panels.extend(new)
        self._selected = new[-1].id if new else self._selected
        self._changed(controls=True)

    def add_idiogram(self) -> None:
        if self._folder is None:
            return
        fig = self._ensure()
        p = fg.add_idiogram(fig)
        self._selected = p.id
        self._changed(controls=True)

    def _add_current(self, *, all_channels: bool) -> None:
        if not self.current_image:
            QMessageBox.information(self, "Figure", "Open a cell first (double-click it in the Library), "
                                                    "or use Add to figure on an image in the Library.")
            return
        self.add_cell(self.current_image, all_channels=all_channels)

    def _select_panel(self, panel_id: str) -> None:
        self._selected = panel_id
        self._page.selected = panel_id
        self._page.update()
        self._fill_panels()

    def _set_panel(self, panel: fg.Panel, **values) -> None:
        for key, value in values.items():
            setattr(panel, key, value)
        self._changed()

    def _move(self, step: int) -> None:
        fg.move_panel(self._fig, self._selected, step)
        self._changed()

    def _remove(self) -> None:
        fg.remove_panel(self._fig, self._selected)
        self._selected = ""
        self._changed()

    def _on_page_changed(self) -> None:
        self._changed()

    def _changed(self, *, controls: bool = False) -> None:
        self._redraw(controls=controls, legend=not (self._fig and self._fig.legend))
        self._save_timer.start()

    def _on_which(self, index: int) -> None:
        data = self._which.itemData(index)
        self.flush()
        if data == "\x00new":
            name, ok = QInputDialog.getText(self, "New figure", "Name:", text=fg.new_figure_name(self._folder))
            if not ok or not name.strip():
                self._fill_controls()
                return
            self._fig = fg.Figure(name=name.strip())
            self._old_name = ""
            self._save_timer.start()
        else:
            self._fig = next((f for f in fg.list_figures(self._folder) if f.name == data), self._fig)
            self._old_name = self._fig.name
        self._selected = ""
        _write("figure_name", self._fig.name)
        self._redraw(controls=True)
        self.crumbChanged.emit(self.crumb())

    def _on_name(self) -> None:
        fig = self._fig
        text = self._name.text().strip()
        if fig is None or not text or text == fig.name:
            return
        if any(f.name == text for f in fg.list_figures(self._folder)):
            QMessageBox.information(self, "Figure", f"A figure named \"{text}\" exists already.")
            self._name.setText(fig.name)
            return
        fig.name = text
        _write("figure_name", text)
        self.flush()
        self._fill_controls()
        self.crumbChanged.emit(self.crumb())

    def _on_page_setting(self, *_args) -> None:
        fig = self._fig
        if fig is None or self._busy:
            return
        fig.width_mm = float(self._width.currentData())
        fig.custom_mm = float(self._custom.value())
        self._custom.setVisible(fig.width_mm <= 0)
        fig.layout = str(self._layout.currentData())
        fig.spacing = int(self._spacing.value())
        fig.font_pt = float(self._font_size.value())
        fig.font_family = self._font.currentFont().family()
        fig.rows, fig.cols = int(self._rows.value()), int(self._cols.value())
        for w in self._grid_row:
            w.setVisible(fig.layout == "custom")
        self._changed()

    def _on_marks(self, *_args) -> None:
        fig = self._fig
        if fig is None or self._busy:
            return
        fig.same_crop = self._same_crop.isChecked()
        fig.letters = self._letters.isChecked()
        fig.channel_labels = self._labels.isChecked()
        fig.scale_bar = self._bar.isChecked()
        fig.scale_um = float(self._bar_len.value())
        fig.arrowheads = self._arrows.isChecked()
        fig.chromosome_numbers = self._numbers.isChecked()
        self._changed()

    def _on_number(self) -> None:
        fig = self._fig
        if fig is not None and self._number.text() != fig.number:
            fig.number = self._number.text()
            self._changed()

    def _on_legend_edited(self) -> None:
        if self._fig is not None:
            self._fig.legend = self._page.legend.toHtml()
            self._save_timer.start()

    def _rewrite_legend(self) -> None:
        if self._fig is not None:
            self._fig.legend = ""
            self._changed()

    def _copy_methods(self) -> None:
        if self._fig is not None and self._folder is not None:
            QGuiApplication.clipboard().setText(fg.methods_text(self._folder, self._fig))

    def _copy_legend(self) -> None:
        QGuiApplication.clipboard().setText(self._page.legend.toPlainText())

    def flush(self) -> None:
        self._save_timer.stop()
        if self._fig is None or self._folder is None:
            return
        try:
            fg.save_figure(self._folder, self._fig, old_name=self._old_name)
            self._old_name = self._fig.name
        except OSError as exc:
            QMessageBox.warning(self, "Figure", f"Could not save the figure:\n{exc}")

    # ------------------------------------------------------------------ export

    def _export(self, kind: str) -> None:
        fig = self._fig
        if fig is None or not fig.panels or self._folder is None:
            QMessageBox.information(self, "Save", "The figure has no panels yet.")
            return
        suffix = {"svg": ".svg", "pdf": ".pdf", "tiff": ".tif"}[kind]
        filt = {"svg": "SVG (*.svg)", "pdf": "PDF (*.pdf)", "tiff": "TIFF (*.tif *.tiff)"}[kind]
        default = str(self._folder / f"{fg._safe(fig.name)}{suffix}")
        path, _ = QFileDialog.getSaveFileName(self, f"Save {fig.name}", default, filt)
        if not path:
            return
        _boxes, h_mm = self._painter.layout()
        w_mm = fig.page_width
        try:
            if kind == "svg":
                from PySide6.QtSvg import QSvgGenerator

                k = 96 / 25.4
                gen = QSvgGenerator()
                gen.setFileName(path)
                gen.setSize(QSizeF(w_mm * k, h_mm * k).toSize())
                gen.setViewBox(QRectF(0, 0, w_mm * k, h_mm * k))
                gen.setResolution(96)
                gen.setTitle(fig.name)
                p = QPainter(gen)
                self._painter.paint(p, k)
                p.end()
            elif kind == "pdf":
                from PySide6.QtGui import QPdfWriter

                writer = QPdfWriter(path)
                writer.setResolution(600)
                writer.setPageSize(QPageSize(QSizeF(w_mm, h_mm), QPageSize.Unit.Millimeter, fig.name))
                from PySide6.QtCore import QMarginsF
                from PySide6.QtGui import QPageLayout

                writer.setPageMargins(QMarginsF(0, 0, 0, 0), QPageLayout.Unit.Millimeter)
                p = QPainter(writer)
                self._painter.paint(p, 600 / 25.4)
                p.end()
            else:
                k = 600 / 25.4
                image = QImage(fg.mm_to_px(w_mm, 600), fg.mm_to_px(h_mm, 600), QImage.Format.Format_RGB32)
                image.fill(QColor("#FFFFFF"))
                dpm = int(round(600 / 0.0254))
                image.setDotsPerMeterX(dpm)
                image.setDotsPerMeterY(dpm)
                p = QPainter(image)
                self._painter.paint(p, k)
                p.end()
                if not image.save(path, "TIFF"):
                    alt = str(Path(path).with_suffix(".png"))
                    if not image.save(alt, "PNG"):
                        raise OSError(f"Could not write {path}")
                    QMessageBox.information(self, "Save", f"TIFF is not available in this Qt; saved a 600 dpi PNG instead:\n{alt}")
                    return
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Save", f"Could not save the figure:\n{exc}")
            return
        doc = QTextDocument()
        doc.setHtml(self._page.legend.toHtml())
        QGuiApplication.clipboard().setText(doc.toPlainText())
        QMessageBox.information(self, "Save", f"Saved {Path(path).name}.\nThe legend text was copied to the clipboard.")
