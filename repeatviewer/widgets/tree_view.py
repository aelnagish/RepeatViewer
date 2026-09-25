"""Tree tab of the Preview panel: an interactive view of a Newick tree.

Opens ``.treefile`` (IQ-TREE), ``.nwk``, ``.newick``, ``.tre`` and ``.tree``.
Rectangular, circular and radial layouts; tip names, support or branch-length
labels, font size, scale bar, zoom and pan; export to SVG, PDF and PNG.

The geometry comes from :mod:`repeatviewer.backends.tree_layout` (pure,
unit-tested); this module only turns it into a ``QGraphicsScene``. Exported
figures are always drawn in black on white, whatever the app theme.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from pathlib import Path

from PySide6.QtCore import QPointF, QRectF, QSizeF, Qt, Signal
from PySide6.QtGui import (
    QBrush,
    QColor,
    QFont,
    QFontMetricsF,
    QImage,
    QPageSize,
    QPainter,
    QPainterPath,
    QPalette,
    QPen,
)
from PySide6.QtWidgets import (
    QCheckBox,
    QFileDialog,
    QGraphicsScene,
    QGraphicsSimpleTextItem,
    QGraphicsView,
    QHBoxLayout,
    QLabel,
    QMenu,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import tree_layout as tl
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.spin_boxes import NoWheelSpinBox

BRANCH_LABELS_NONE = "none"
BRANCH_LABELS_SUPPORT = "support"
BRANCH_LABELS_LENGTH = "length"

EXPORT_FORMATS = (("SVG", "svg"), ("PDF", "pdf"), ("PNG", "png"))
PNG_SCALE = 4.0            # 4 × 72 dpi ≈ 300 dpi
_ZOOM_STEP = 1.12
_TRACKPAD_ZOOM_BASE = 1.0022
_MIN_TREE_PX = 260.0       # smallest drawing size (rectangular width / radius)
_RECT_TREE_WIDTH_PX = 520.0
_MARGIN_PX = 24.0

SUPPORT_HELP = (
    "<b>Support</b>: the label on each internal node.<br>"
    "IQ-TREE shows <b>SH-aLRT/UFBoot</b>, e.g. 92/97. A branch is well supported "
    "at SH-aLRT 80 or more and UFBoot 95 or more.<br>"
    "FastME bootstrap: 70 or more. FastTree: local support from 0 to 1, close to 1 is strong.<br>"
    "<b>Branch length</b>: substitutions per site."
)


@dataclass
class TreeViewOptions:
    layout: str = tl.LAYOUT_RECTANGULAR
    order: str = tl.ORDER_FILE
    tip_labels: bool = True
    branch_labels: str = BRANCH_LABELS_SUPPORT
    font_pt: int = 10
    scale_bar: bool = True


# The last settings used, so the next tree opens the same way (this session).
_last_options = TreeViewOptions()


def export_default_path(tree_path: Path, layout: str, ext: str) -> Path:
    """Rule 1: keep the full input name, then a suffix (``x.treefile_circular.svg``)."""
    return tree_path.parent / f"{tree_path.name}_{layout}.{ext}"


# --------------------------------------------------------------------------
# Scene building
# --------------------------------------------------------------------------


def _pixel_scale(layout: tl.Layout, n_tips: int, row_h: float) -> tuple[float, float]:
    """(sx, sy) from layout units to pixels."""
    extent = layout.extent if layout.extent > 0 else 1.0
    if layout.kind == tl.LAYOUT_RECTANGULAR:
        return _RECT_TREE_WIDTH_PX / extent, row_h
    if layout.kind == tl.LAYOUT_CIRCULAR:
        # Leave one row of height per tip around the circle so labels do not overlap.
        radius = max(_MIN_TREE_PX, n_tips * row_h / (2 * math.pi))
    else:
        radius = max(_MIN_TREE_PX, math.sqrt(n_tips) * row_h * 2.0)
    s = radius / extent
    return s, s


def build_tree_scene(
    scene: QGraphicsScene,
    model: tl.TreeModel,
    opts: TreeViewOptions,
    *,
    ink: QColor,
    muted: QColor,
    for_export: bool = False,
) -> None:
    """Fill ``scene`` with the tree drawn per ``opts``."""
    scene.clear()
    shown = tl.ordered(model, opts.order)
    layout = tl.compute_layout(shown, opts.layout)
    font = QFont()
    font.setPointSizeF(float(opts.font_pt))
    small = QFont(font)
    small.setPointSizeF(max(5.0, opts.font_pt * 0.8))
    fm = QFontMetricsF(font)
    fm_small = QFontMetricsF(small)
    row_h = fm.height() * 1.3
    sx, sy = _pixel_scale(layout, shown.n_tips, row_h)

    def px(point: tuple[float, float]) -> QPointF:
        return QPointF(point[0] * sx, point[1] * sy)

    pen = QPen(ink, 1.2)
    pen.setCosmetic(not for_export)
    pen.setCapStyle(Qt.PenCapStyle.FlatCap)
    path = QPainterPath()
    for x0, y0, x1, y1 in layout.lines:
        path.moveTo(x0 * sx, y0 * sy)
        path.lineTo(x1 * sx, y1 * sy)
    for arc in layout.arcs:
        r = arc.radius * sx
        rect = QRectF(-r, -r, 2 * r, 2 * r)
        path.arcMoveTo(rect, arc.start_deg)
        path.arcTo(rect, arc.start_deg, arc.span_deg)
    scene.addPath(path, pen)

    brush = QBrush(ink)
    gap = fm.averageCharWidth() * 0.6
    if opts.tip_labels:
        for tip in shown.tips():
            node = shown.nodes[tip]
            item = QGraphicsSimpleTextItem(node.name)
            item.setFont(font)
            item.setBrush(brush)
            if not for_export:
                item.setToolTip(f"{node.name}\nBranch length: {tl.branch_length_text(node.length)}")
            w = item.boundingRect().width()
            h = item.boundingRect().height()
            p = px(layout.xy[tip])
            if layout.kind == tl.LAYOUT_RECTANGULAR:
                item.setPos(p.x() + gap, p.y() - h / 2.0)
            else:
                angle = layout.tip_angle[tip]
                rad = math.radians(angle)
                anchor = QPointF(p.x() + gap * math.cos(rad), p.y() - gap * math.sin(rad))
                rotation, flipped = tl.label_rotation(angle)
                if flipped:
                    item.setTransformOriginPoint(w, h / 2.0)
                    item.setPos(anchor.x() - w, anchor.y() - h / 2.0)
                else:
                    item.setTransformOriginPoint(0.0, h / 2.0)
                    item.setPos(anchor.x(), anchor.y() - h / 2.0)
                item.setRotation(rotation)
            scene.addItem(item)

    if opts.branch_labels != BRANCH_LABELS_NONE:
        label_brush = QBrush(muted)
        for node in shown.nodes:
            if node.parent is None:
                continue
            if opts.branch_labels == BRANCH_LABELS_SUPPORT:
                text = node.support
            else:
                text = tl.branch_length_text(node.length)
            if not text:
                continue
            item = QGraphicsSimpleTextItem(text)
            item.setFont(small)
            item.setBrush(label_brush)
            w = item.boundingRect().width()
            h = fm_small.height()
            if opts.branch_labels == BRANCH_LABELS_SUPPORT:
                p = px(layout.xy[node.index])
                if layout.kind == tl.LAYOUT_RECTANGULAR:
                    item.setPos(p.x() - w - 2.0, p.y() - h)
                else:
                    item.setPos(p.x() + 2.0, p.y() - h)
            else:
                p = px(layout.branch_mid[node.index])
                item.setPos(p.x() - w / 2.0, p.y() - h)
            scene.addItem(item)

    if opts.scale_bar and shown.has_lengths:
        value = tl.nice_scale_length(layout.extent)
        if value > 0:
            bounds = scene.itemsBoundingRect()
            length_px = value * sx
            y = bounds.bottom() + row_h
            x = bounds.left()
            bar = QPainterPath()
            bar.moveTo(x, y)
            bar.lineTo(x + length_px, y)
            scene.addPath(bar, pen)
            label = QGraphicsSimpleTextItem(tl.scale_text(value))
            label.setFont(small)
            label.setBrush(brush)
            label.setPos(x + (length_px - label.boundingRect().width()) / 2.0, y + 2.0)
            scene.addItem(label)

    rect = scene.itemsBoundingRect().adjusted(-_MARGIN_PX, -_MARGIN_PX, _MARGIN_PX, _MARGIN_PX)
    scene.setSceneRect(rect)


# --------------------------------------------------------------------------
# Export
# --------------------------------------------------------------------------


def export_tree(model: tl.TreeModel, opts: TreeViewOptions, out: Path, *, title: str = "") -> Path:
    """Write the tree to ``out`` (.svg, .pdf or .png), black on white."""
    scene = QGraphicsScene()
    build_tree_scene(scene, model, opts, ink=QColor(0, 0, 0), muted=QColor(70, 70, 70),
                     for_export=True)
    rect = scene.sceneRect()
    ext = out.suffix.lower()
    if ext == ".svg":
        from PySide6.QtSvg import QSvgGenerator

        gen = QSvgGenerator()
        gen.setFileName(str(out))
        gen.setSize(rect.size().toSize())
        gen.setViewBox(QRectF(0, 0, rect.width(), rect.height()))
        gen.setTitle(title or out.stem)
        painter = QPainter(gen)
        try:
            painter.fillRect(QRectF(0, 0, rect.width(), rect.height()), QColor(255, 255, 255))
            scene.render(painter, QRectF(0, 0, rect.width(), rect.height()), rect)
        finally:
            painter.end()
    elif ext == ".pdf":
        from PySide6.QtGui import QPageLayout, QPdfWriter
        from PySide6.QtCore import QMarginsF

        writer = QPdfWriter(str(out))
        writer.setResolution(300)
        writer.setTitle(title or out.stem)
        writer.setPageLayout(QPageLayout(
            QPageSize(QSizeF(rect.width(), rect.height()), QPageSize.Unit.Point),
            QPageLayout.Orientation.Portrait, QMarginsF(0, 0, 0, 0),
        ))
        painter = QPainter(writer)
        try:
            target = QRectF(painter.viewport())
            painter.fillRect(target, QColor(255, 255, 255))
            scene.render(painter, target, rect)
        finally:
            painter.end()
    elif ext == ".png":
        image = QImage(int(math.ceil(rect.width() * PNG_SCALE)),
                       int(math.ceil(rect.height() * PNG_SCALE)), QImage.Format.Format_ARGB32)
        image.fill(QColor(255, 255, 255))
        image.setDotsPerMeterX(int(72 * PNG_SCALE / 0.0254))
        image.setDotsPerMeterY(int(72 * PNG_SCALE / 0.0254))
        painter = QPainter(image)
        try:
            painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
            scene.render(painter, QRectF(0, 0, image.width(), image.height()), rect)
        finally:
            painter.end()
        if not image.save(str(out), "PNG"):
            raise OSError(f"Could not write {out}")
    else:
        raise ValueError(f"Unsupported export format: {ext}")
    return out


# --------------------------------------------------------------------------
# Widgets
# --------------------------------------------------------------------------


class TreeCanvas(QGraphicsView):
    """Zoom with the wheel or trackpad, drag to pan, double-click to fit."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setScene(QGraphicsScene(self))
        self.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        self.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
        self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)

    def zoom_by(self, factor: float) -> None:
        current = self.transform().m11()
        target = min(40.0, max(0.02, current * factor))
        self.scale(target / current, target / current)

    def fit(self, *, max_scale: float = 1.0) -> None:
        rect = self.scene().sceneRect()
        if rect.isEmpty():
            return
        self.resetTransform()
        view = self.viewport().rect()
        if view.width() <= 0 or view.height() <= 0:
            return
        s = min(view.width() / rect.width(), view.height() / rect.height(), max_scale)
        self.scale(s, s)
        self.centerOn(rect.center())

    def wheelEvent(self, event) -> None:  # noqa: N802 (Qt API)
        pixel = event.pixelDelta()
        angle = event.angleDelta()
        if not pixel.isNull() and (pixel.y() != 0 or pixel.x() != 0):
            dy = float(pixel.y() if pixel.y() != 0 else pixel.x())
            factor = _TRACKPAD_ZOOM_BASE ** dy
        elif angle.y() != 0:
            factor = _ZOOM_STEP ** (angle.y() / 120.0)
        else:
            return
        self.zoom_by(factor)
        event.accept()

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802 (Qt API)
        self.fit(max_scale=4.0)
        event.accept()


class TreePanel(QWidget):
    """Controls + canvas for one tree file."""

    fileWritten = Signal(object)  # Path of an exported figure

    def __init__(self, path: Path | None = None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("TreePanel")
        self._path: Path | None = None
        self._model: tl.TreeModel | None = None
        self._opts = replace(_last_options)
        self._block = False

        self._layout_box = NoWheelComboBox(self)
        for label, value in (("Rectangular", tl.LAYOUT_RECTANGULAR),
                             ("Circular", tl.LAYOUT_CIRCULAR),
                             ("Radial (unrooted)", tl.LAYOUT_RADIAL)):
            self._layout_box.addItem(label, value)
        self._order_box = NoWheelComboBox(self)
        for label, value in (("As in file", tl.ORDER_FILE),
                             ("Increasing", tl.ORDER_INCREASING),
                             ("Decreasing", tl.ORDER_DECREASING)):
            self._order_box.addItem(label, value)
        self._order_box.setToolTip("Sort each node's branches by clade size (ladderize).")
        self._branch_box = NoWheelComboBox(self)
        for label, value in (("Support", BRANCH_LABELS_SUPPORT),
                             ("Branch length", BRANCH_LABELS_LENGTH),
                             ("None", BRANCH_LABELS_NONE)):
            self._branch_box.addItem(label, value)
        self._tips_check = QCheckBox("Tip names", self)
        self._scale_check = QCheckBox("Scale bar", self)
        self._font_spin = NoWheelSpinBox(self)
        self._font_spin.setRange(5, 36)
        self._font_spin.setSuffix(" pt")
        self._font_spin.setFixedWidth(70)
        fit_btn = QPushButton("Fit", self)
        fit_btn.setToolTip("Show the whole tree (or double-click the tree).")
        fit_btn.clicked.connect(lambda: self._canvas.fit(max_scale=4.0))
        export_btn = QPushButton("Export", self)
        menu = QMenu(export_btn)
        for label, ext in EXPORT_FORMATS:
            menu.addAction(f"{label}…", lambda checked=False, e=ext: self._export(e))
        export_btn.setMenu(menu)

        row1 = QHBoxLayout()
        row1.setSpacing(6)
        row1.addWidget(QLabel("Layout:", self))
        row1.addWidget(self._layout_box)
        row1.addSpacing(8)
        row1.addWidget(QLabel("Order:", self))
        row1.addWidget(self._order_box)
        row1.addSpacing(8)
        row1.addWidget(QLabel("Branch labels:", self))
        row1.addWidget(self._branch_box)
        row1.addWidget(make_help_button(self, "Branch labels", SUPPORT_HELP))
        row1.addStretch(1)
        row2 = QHBoxLayout()
        row2.setSpacing(6)
        row2.addWidget(self._tips_check)
        row2.addWidget(self._scale_check)
        row2.addSpacing(8)
        row2.addWidget(QLabel("Font:", self))
        row2.addWidget(self._font_spin)
        row2.addStretch(1)
        row2.addWidget(fit_btn)
        row2.addWidget(export_btn)

        self._canvas = TreeCanvas(self)
        self._status = QLabel(self)
        self._status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 6, 8, 4)
        layout.setSpacing(4)
        layout.addLayout(row1)
        layout.addLayout(row2)
        layout.addWidget(self._canvas, stretch=1)
        layout.addWidget(self._status)

        self._apply_options_to_controls()
        for box in (self._layout_box, self._order_box, self._branch_box):
            box.currentIndexChanged.connect(self._on_controls_changed)
        for check in (self._tips_check, self._scale_check):
            check.toggled.connect(self._on_controls_changed)
        self._font_spin.valueChanged.connect(self._on_controls_changed)

        if path is not None:
            self.load_tree(path)

    # ------------------------------------------------------------ state

    def options(self) -> TreeViewOptions:
        return replace(self._opts)

    def model(self) -> tl.TreeModel | None:
        return self._model

    def path(self) -> Path | None:
        return self._path

    def scene(self) -> QGraphicsScene:
        return self._canvas.scene()

    def _apply_options_to_controls(self) -> None:
        self._block = True
        try:
            for box, value in ((self._layout_box, self._opts.layout),
                               (self._order_box, self._opts.order),
                               (self._branch_box, self._opts.branch_labels)):
                idx = box.findData(value)
                box.setCurrentIndex(max(0, idx))
            self._tips_check.setChecked(self._opts.tip_labels)
            self._scale_check.setChecked(self._opts.scale_bar)
            self._font_spin.setValue(self._opts.font_pt)
        finally:
            self._block = False

    def set_options(self, opts: TreeViewOptions) -> None:
        self._opts = replace(opts)
        self._apply_options_to_controls()
        self._redraw(refit=True)

    def _on_controls_changed(self, *_args) -> None:
        if self._block:
            return
        global _last_options
        old_layout = self._opts.layout
        self._opts = TreeViewOptions(
            layout=str(self._layout_box.currentData()),
            order=str(self._order_box.currentData()),
            tip_labels=self._tips_check.isChecked(),
            branch_labels=str(self._branch_box.currentData()),
            font_pt=self._font_spin.value(),
            scale_bar=self._scale_check.isChecked(),
        )
        _last_options = replace(self._opts)
        self._redraw(refit=self._opts.layout != old_layout)

    # ------------------------------------------------------------ loading

    def load_tree(self, path: Path) -> bool:
        self._path = Path(path).expanduser().resolve()
        try:
            self._model = tl.read_tree_file(self._path)
        except (OSError, tl.TreeReadError) as exc:
            self._model = None
            self._canvas.scene().clear()
            self._status.setText(f"Could not read {self._path.name}: {exc}")
            return False
        if not self._model.has_support and self._opts.branch_labels == BRANCH_LABELS_SUPPORT:
            self._block = True
            self._branch_box.setCurrentIndex(self._branch_box.findData(BRANCH_LABELS_NONE))
            self._block = False
            self._opts = replace(self._opts, branch_labels=BRANCH_LABELS_NONE)
        self._redraw(refit=True)
        return True

    def _palette_colors(self) -> tuple[QColor, QColor]:
        pal = self.palette()
        ink = pal.color(QPalette.ColorRole.WindowText)
        muted = QColor(ink)
        muted.setAlphaF(0.7)
        return ink, muted

    def _redraw(self, *, refit: bool) -> None:
        if self._model is None:
            return
        ink, muted = self._palette_colors()
        self._canvas.setBackgroundBrush(self.palette().color(QPalette.ColorRole.Base))
        build_tree_scene(self._canvas.scene(), self._model, self._opts, ink=ink, muted=muted)
        if refit:
            self._canvas.fit(max_scale=1.0)
        m = self._model
        parts = [f"{m.n_tips} tips", f"{m.n_internal} internal nodes"]
        if not m.has_lengths:
            parts.append("no branch lengths in the file")
        if m.has_support and self._opts.branch_labels == BRANCH_LABELS_SUPPORT:
            if any("/" in n.support for n in m.nodes):
                parts.append("support: SH-aLRT/UFBoot")
        self._status.setText(" · ".join(parts))

    def showEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().showEvent(event)
        if self._model is not None:
            self._canvas.fit(max_scale=1.0)

    # ------------------------------------------------------------ export

    def _export(self, ext: str) -> None:
        if self._model is None or self._path is None:
            QMessageBox.information(self, "Export tree", "There is no tree to export.")
            return
        default = export_default_path(self._path, self._opts.layout, ext)
        label = dict((e, lab) for lab, e in EXPORT_FORMATS)[ext]
        chosen, _ = QFileDialog.getSaveFileName(
            self, f"Export tree as {label}", str(default), f"{label} (*.{ext})"
        )
        if not chosen:
            return
        out = Path(chosen)
        if out.suffix.lower() != f".{ext}":
            out = out.with_name(out.name + f".{ext}")
        try:
            export_tree(self._model, self._opts, out, title=self._path.name)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Export tree", f"Could not write {out.name}:\n{exc}")
            return
        self._status.setText(f"Saved {out.name}")
        self.fileWritten.emit(out)
