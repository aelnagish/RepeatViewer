"""Tree tab of the Preview panel: an interactive view of a Newick tree.

Opens ``.treefile`` (IQ-TREE), ``.nwk``, ``.newick``, ``.tre`` and ``.tree``.

* View: rectangular, circular and radial layouts; tip names, support or
  branch-length labels, font size, scale bar, zoom and pan.
* Edit (right-click a branch, or the Edit menu for the selected branch):
  reroot here, midpoint root, swap branches, collapse / expand a clade,
  colour a clade, show one clade only; collapse all branches below a support
  value; undo (Cmd+Z / Ctrl+Z).
* Keep: colours, collapsed clades and view settings are remembered in a hidden
  file next to the tree (``.<tree name>.view.json``). Rerooting and swapping
  change the tree itself, which is never overwritten: "Save tree as…" writes a
  new file (``<tree name>_rerooted.treefile``, rule 1).
* Export: SVG, PDF and PNG of exactly what is shown, black on white.

The geometry and the edits come from the pure, unit-tested modules
:mod:`repeatviewer.backends.tree_layout` and :mod:`repeatviewer.backends.tree_edit`;
this module only turns them into a ``QGraphicsScene``.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

from PySide6.QtCore import QPointF, QRectF, QSizeF, Qt, QTimer, QUrl, Signal
from PySide6.QtGui import (
    QBrush,
    QDesktopServices,
    QLinearGradient,
    QColor,
    QFont,
    QFontMetricsF,
    QImage,
    QKeySequence,
    QPageSize,
    QPainter,
    QPainterPath,
    QPainterPathStroker,
    QPalette,
    QPen,
    QPolygonF,
    QShortcut,
)
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QColorDialog,
    QFileDialog,
    QGraphicsPathItem,
    QGraphicsPolygonItem,
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

from repeatviewer.backends import tree_edit as te
from repeatviewer.backends import tree_layout as tl
from repeatviewer.backends import tree_traits as tt
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
_HIT_WIDTH = 7.0           # clickable width around a branch, in pixels
_SELECTED = QColor(0, 110, 230)
ROLE_NODE = 0              # QGraphicsItem.data key: display node index

SUPPORT_HELP = (
    "<b>Support</b>: the label on each internal node.<br>"
    "IQ-TREE shows <b>SH-aLRT/UFBoot</b>, e.g. 92/97. A branch is well supported "
    "at SH-aLRT 80 or more and UFBoot 95 or more.<br>"
    "FastME bootstrap: 70 or more. FastTree: local support from 0 to 1, close to 1 is strong.<br>"
    "<b>Branch length</b>: substitutions per site."
)
COLLAPSE_HELP = (
    "<b>Collapse below support</b>: branches with support under this value are merged "
    "into their parent, so only well-supported groupings remain.<br>"
    "IQ-TREE: the UFBoot value is used (the second number, e.g. 97 in 92/97). "
    "FastTree values (0 to 1) are read as percentages.<br>"
    "0 = off. The tree file is not changed; use <b>Save tree as…</b> to keep the result."
)
EDIT_HELP = (
    "<b>Right-click a branch</b> (or select it and use <b>Edit</b>) to reroot there, swap "
    "its branches, collapse or colour its clade, or show that clade only.<br>"
    "Undo: Cmd+Z / Ctrl+Z.<br>"
    "Colours and collapsed clades are remembered for this tree. Rerooting and swapping "
    "change the tree: the file is never overwritten, use <b>Save tree as…</b>."
)


TRAITS_HELP = (
    "<b>Colour tips by</b> a column of the trait table next to the tree "
    "(<i>tree name</i>_traits.tsv): chromosome, tandem array, methylation, or any column "
    "you add yourself in a spreadsheet.<br>"
    "Text columns get distinct colour-blind-safe colours; number columns a colour scale "
    "(low = purple, high = yellow). A legend is drawn next to the tree and is part of "
    "every export.<br>"
    "<b>Trait table → Create or update</b> fills the table from the sequence headers "
    "(and a methylation file, if given)."
)


@dataclass
class TreeViewOptions:
    layout: str = tl.LAYOUT_RECTANGULAR
    order: str = tl.ORDER_FILE
    tip_labels: bool = True
    branch_labels: str = BRANCH_LABELS_SUPPORT
    font_pt: int = 10
    scale_bar: bool = True


@dataclass
class TreeViewState:
    """What the user did to this tree view (clades by their tip names)."""

    collapsed: set[tuple[str, ...]] = field(default_factory=set)
    colors: dict[tuple[str, ...], str] = field(default_factory=dict)
    subtree: tuple[str, ...] | None = None
    collapse_below: int = 0          # support threshold, 0 = off
    edited: str = ""                 # "", "rerooted" or "edited": the tree differs from the file
    color_by: str = ""               # trait-table column colouring the tips ("" = none)
    trait_settings: dict = field(default_factory=dict)  # methylation_file, array_gap

    def copy(self) -> "TreeViewState":
        return TreeViewState(set(self.collapsed), dict(self.colors), self.subtree,
                             self.collapse_below, self.edited, self.color_by,
                             dict(self.trait_settings))


@dataclass
class Display:
    """The tree as drawn: display model plus per-node clade keys and colours."""

    model: tl.TreeModel
    info: dict[int, te.CollapsedInfo]
    keys: list[tuple[str, ...]]
    colors: list[str | None]
    n_contracted: int = 0
    tip_colors: dict[int, str] = field(default_factory=dict)   # trait colour per tip
    tip_notes: dict[int, str] = field(default_factory=dict)    # trait values for tooltips
    legend: tt.ColourScale | None = None


# The last settings used, so the next tree opens the same way (this session).
_last_options = TreeViewOptions()


def export_default_path(tree_path: Path, layout: str, ext: str) -> Path:
    """Rule 1: keep the full input name, then a suffix (``x.treefile_circular.svg``)."""
    return tree_path.parent / f"{tree_path.name}_{layout}.{ext}"


def saved_tree_default_path(tree_path: Path, edited: str) -> Path:
    """Rule 1: ``x.treefile_rerooted.treefile`` (never the file itself)."""
    return tree_path.parent / f"{tree_path.name}_{edited or 'edited'}.treefile"


def view_state_path(tree_path: Path) -> Path:
    """Hidden settings file next to the tree (not listed in the file panel)."""
    return tree_path.parent / f".{tree_path.name}.view.json"


# --------------------------------------------------------------------------
# Display pipeline (no drawing)
# --------------------------------------------------------------------------


def tree_structure(model: tl.TreeModel, state: TreeViewState) -> tuple[tl.TreeModel, int]:
    """The tree after weak branches are collapsed and a clade is shown alone."""
    n = 0
    if state.collapse_below > 0:
        model, n = te.collapse_weak(model, state.collapse_below)
    if state.subtree:
        idx = te.find_clade(model, state.subtree)
        if idx is not None and not model.nodes[idx].is_tip:
            model = te.subtree(model, idx)
    return model, n


def make_display(model: tl.TreeModel, opts: TreeViewOptions, state: TreeViewState) -> Display:
    m, n_contracted = tree_structure(model, state)
    collapsed = set()
    for key in state.collapsed:
        idx = te.find_clade(m, key)
        if idx is not None and not m.nodes[idx].is_tip and idx != 0:
            collapsed.add(idx)
    d, info = te.display_model(m, collapsed)
    d = tl.ordered(d, opts.order)
    keys: list[tuple[str, ...]] = [()] * len(d.nodes)
    for node in reversed(d.nodes):
        i = node.index
        if i in info:
            keys[i] = info[i].key
        elif node.is_tip:
            keys[i] = (node.name,)
        else:
            acc: set[str] = set()
            for c in node.children:
                acc.update(keys[c])
            keys[i] = tuple(sorted(acc))
    colors: list[str | None] = [None] * len(d.nodes)
    for node in d.nodes:
        own = state.colors.get(keys[node.index])
        colors[node.index] = own or (colors[node.parent] if node.parent is not None else None)
    return Display(d, info, keys, colors, n_contracted)


# --------------------------------------------------------------------------
# Scene building
# --------------------------------------------------------------------------


class _BranchItem(QGraphicsPathItem):
    """One node's branch (and its connector); easy to hit with the mouse."""

    def shape(self) -> QPainterPath:  # noqa: D401 (Qt API)
        stroker = QPainterPathStroker()
        stroker.setWidth(_HIT_WIDTH)
        return stroker.createStroke(self.path())


def _pixel_scale(layout: tl.Layout, n_tips: int, row_h: float) -> tuple[float, float]:
    """(sx, sy) from layout units to pixels."""
    extent = layout.extent if layout.extent > 0 else 1.0
    if layout.kind == tl.LAYOUT_RECTANGULAR:
        return _RECT_TREE_WIDTH_PX / extent, row_h
    if layout.kind == tl.LAYOUT_CIRCULAR:
        radius = max(_MIN_TREE_PX, n_tips * row_h / (2 * math.pi))
    else:
        radius = max(_MIN_TREE_PX, math.sqrt(n_tips) * row_h * 2.0)
    s = radius / extent
    return s, s


def build_tree_scene(
    scene: QGraphicsScene,
    display: Display | tl.TreeModel,
    opts: TreeViewOptions,
    *,
    ink: QColor,
    muted: QColor,
    for_export: bool = False,
    selected: tuple[str, ...] | None = None,
) -> None:
    """Fill ``scene`` with the tree drawn per ``opts``."""
    if isinstance(display, tl.TreeModel):
        display = make_display(display, opts, TreeViewState())
    scene.clear()
    shown = display.model
    info = display.info
    layout = tl.compute_layout(shown, opts.layout, {i: c.depth for i, c in info.items()})
    font = QFont()
    font.setPointSizeF(float(opts.font_pt))
    small = QFont(font)
    small.setPointSizeF(max(5.0, opts.font_pt * 0.8))
    fm = QFontMetricsF(font)
    fm_small = QFontMetricsF(small)
    row_h = fm.height() * 1.3
    sx, sy = _pixel_scale(layout, shown.n_tips, row_h)
    radial = layout.kind != tl.LAYOUT_RECTANGULAR

    def px(point: tuple[float, float]) -> QPointF:
        return QPointF(point[0] * sx, point[1] * sy)

    def color_for(i: int) -> QColor:
        c = display.colors[i]
        return QColor(c) if c else QColor(ink)

    def pen_for(i: int, width: float = 1.2) -> QPen:
        is_sel = selected is not None and display.keys[i] == selected and not for_export
        pen = QPen(_SELECTED if is_sel else color_for(i), 2.6 if is_sel else width)
        pen.setCosmetic(not for_export)
        pen.setCapStyle(Qt.PenCapStyle.FlatCap)
        return pen

    # --- branches: one item per node (its own branch + its connector) -------
    depth = tl.depths(shown)
    for node in shown.nodes:
        i = node.index
        path = QPainterPath()
        if opts.layout == tl.LAYOUT_CIRCULAR:
            a = math.radians(layout.node_angle[i])
            if node.parent is not None:
                r0 = depth[node.parent] * sx
                path.moveTo(r0 * math.cos(a), -r0 * math.sin(a))
                path.lineTo(px(layout.xy[i]))
            if not node.is_tip and depth[i] > 0:
                a0 = layout.node_angle[node.children[0]]
                a1 = layout.node_angle[node.children[-1]]
                if a1 > a0:
                    rr = depth[i] * sx
                    rect = QRectF(-rr, -rr, 2 * rr, 2 * rr)
                    path.arcMoveTo(rect, a0)
                    path.arcTo(rect, a0, a1 - a0)
        elif opts.layout == tl.LAYOUT_RECTANGULAR:
            if node.parent is not None:
                p1 = px(layout.xy[i])
                path.moveTo(px(layout.xy[node.parent]).x(), p1.y())
                path.lineTo(p1)
            if not node.is_tip:
                x = px(layout.xy[i]).x()
                path.moveTo(x, layout.xy[node.children[0]][1] * sy)
                path.lineTo(x, layout.xy[node.children[-1]][1] * sy)
        elif node.parent is not None:
            path.moveTo(px(layout.xy[node.parent]))
            path.lineTo(px(layout.xy[i]))
        if path.isEmpty():
            continue
        item = _BranchItem(path)
        item.setPen(pen_for(i))
        item.setData(ROLE_NODE, i)
        scene.addItem(item)

    # --- collapsed clades: triangles ------------------------------------------
    for i, ci in info.items():
        base = px(layout.xy[i])
        length = ci.depth * sx
        half = row_h * 0.42
        if radial:
            a = math.radians(layout.node_angle.get(i, 0.0))
            dx, dy = math.cos(a), -math.sin(a)
        else:
            dx, dy = 1.0, 0.0
        far = QPointF(base.x() + dx * length, base.y() + dy * length)
        nx, ny = -dy, dx
        tri = QPolygonF([base, QPointF(far.x() + nx * half, far.y() + ny * half),
                         QPointF(far.x() - nx * half, far.y() - ny * half)])
        item = QGraphicsPolygonItem(tri)
        item.setPen(pen_for(i))
        fill = QColor(color_for(i))
        fill.setAlphaF(0.18)
        item.setBrush(QBrush(fill))
        item.setData(ROLE_NODE, i)
        if not for_export:
            item.setToolTip(f"Collapsed clade: {ci.n_tips} sequences\nRight-click to expand")
        scene.addItem(item)

    # --- trait dots at the tips ---------------------------------------------------
    gap = fm.averageCharWidth() * 0.6
    dot = fm.height() * 0.62 if display.tip_colors else 0.0
    for tip, color in display.tip_colors.items():
        p = px(layout.xy[tip])
        extra = info[tip].depth * sx if tip in info else 0.0
        if radial:
            a = math.radians(layout.tip_angle[tip])
            dx, dy = math.cos(a), -math.sin(a)
        else:
            dx, dy = 1.0, 0.0
        reach = extra + gap + dot / 2.0
        cx, cy = p.x() + dx * reach, p.y() + dy * reach
        circle = scene.addEllipse(QRectF(cx - dot / 2.0, cy - dot / 2.0, dot, dot),
                                  QPen(QColor(color).darker(130), 0.6), QBrush(QColor(color)))
        circle.setData(ROLE_NODE, tip)
        if not for_export and tip in display.tip_notes:
            circle.setToolTip(display.tip_notes[tip])
    if dot:
        gap += dot + gap

    # --- tip labels -------------------------------------------------------------
    if opts.tip_labels:
        for tip in shown.tips():
            node = shown.nodes[tip]
            text = node.name
            item = QGraphicsSimpleTextItem(text)
            item.setFont(font)
            item.setBrush(QBrush(color_for(tip)))
            item.setData(ROLE_NODE, tip)
            if not for_export:
                tip_text = f"{text}\nBranch length: {tl.branch_length_text(node.length)}"
                if tip in display.tip_notes:
                    tip_text += "\n" + display.tip_notes[tip]
                item.setToolTip(tip_text)
            w = item.boundingRect().width()
            h = item.boundingRect().height()
            p = px(layout.xy[tip])
            extra = info[tip].depth * sx if tip in info else 0.0
            if layout.kind == tl.LAYOUT_RECTANGULAR:
                item.setPos(p.x() + extra + gap, p.y() - h / 2.0)
            else:
                angle = layout.tip_angle[tip]
                rad = math.radians(angle)
                reach = extra + gap
                anchor = QPointF(p.x() + reach * math.cos(rad), p.y() - reach * math.sin(rad))
                rotation, flipped = tl.label_rotation(angle)
                if flipped:
                    item.setTransformOriginPoint(w, h / 2.0)
                    item.setPos(anchor.x() - w, anchor.y() - h / 2.0)
                else:
                    item.setTransformOriginPoint(0.0, h / 2.0)
                    item.setPos(anchor.x(), anchor.y() - h / 2.0)
                item.setRotation(rotation)
            scene.addItem(item)

    # --- branch labels ------------------------------------------------------------
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

    # --- legend for the trait colours ------------------------------------------------
    if display.legend is not None:
        _draw_legend(scene, display.legend, font, small, ink, muted, for_export)

    # --- scale bar ------------------------------------------------------------------
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
            pen = QPen(ink, 1.2)
            pen.setCosmetic(not for_export)
            scene.addPath(bar, pen)
            label = QGraphicsSimpleTextItem(tl.scale_text(value))
            label.setFont(small)
            label.setBrush(QBrush(ink))
            label.setPos(x + (length_px - label.boundingRect().width()) / 2.0, y + 2.0)
            scene.addItem(label)

    rect = scene.itemsBoundingRect().adjusted(-_MARGIN_PX, -_MARGIN_PX, _MARGIN_PX, _MARGIN_PX)
    scene.setSceneRect(rect)


def _draw_legend(scene: QGraphicsScene, scale: tt.ColourScale, font: QFont, small: QFont,
                 ink: QColor, muted: QColor, for_export: bool) -> None:
    """Legend to the right of the tree: categories, or a colour bar for numbers."""
    bounds = scene.itemsBoundingRect()
    fm = QFontMetricsF(font)
    x = bounds.right() + fm.averageCharWidth() * 4
    y = bounds.top()
    title_font = QFont(font)
    title_font.setBold(True)
    title = QGraphicsSimpleTextItem(scale.column)
    title.setFont(title_font)
    title.setBrush(QBrush(ink))
    title.setPos(x, y)
    scene.addItem(title)
    y += fm.height() * 1.4
    size = fm.height() * 0.62
    if scale.numeric:
        bar_h, bar_w = fm.height() * 6, fm.height() * 0.9
        grad = QLinearGradient(0, y + bar_h, 0, y)   # low at the bottom
        for k in range(11):
            grad.setColorAt(k / 10, QColor(tt.viridis(k / 10)))
        scene.addRect(QRectF(x, y, bar_w, bar_h), QPen(QColor(ink), 0.6), QBrush(grad))
        for value, yy in ((scale.high, y), (scale.low, y + bar_h)):
            label = QGraphicsSimpleTextItem(f"{value:g}")
            label.setFont(small)
            label.setBrush(QBrush(ink))
            label.setPos(x + bar_w + 6, yy - label.boundingRect().height() / 2)
            scene.addItem(label)
        y += bar_h + fm.height() * 0.6
    else:
        for value, color in scale.categories:
            scene.addEllipse(QRectF(x, y + (fm.height() - size) / 2, size, size),
                             QPen(QColor(color).darker(130), 0.6), QBrush(QColor(color)))
            label = QGraphicsSimpleTextItem(value)
            label.setFont(font)
            label.setBrush(QBrush(ink))
            label.setPos(x + size + 6, y)
            scene.addItem(label)
            y += fm.height() * 1.15
    if scale.n_missing:
        note = QGraphicsSimpleTextItem(f"no value: {scale.n_missing}")
        note.setFont(small)
        note.setBrush(QBrush(muted))
        note.setPos(x, y + 2)
        scene.addItem(note)


# --------------------------------------------------------------------------
# Export
# --------------------------------------------------------------------------


def export_tree(model: tl.TreeModel | Display, opts: TreeViewOptions, out: Path, *,
                title: str = "") -> Path:
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
        from PySide6.QtCore import QMarginsF
        from PySide6.QtGui import QPageLayout, QPdfWriter

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


def _fish_files(folder: Path) -> set[str] | None:
    """FASTA files in ``folder`` that have FISH images attached (None: no images at all)."""
    try:
        from repeatviewer.backends.cyto_images import load_images

        images = load_images(folder)
    except Exception:  # noqa: BLE001 - optional information
        return None
    if not images:
        return None
    return {target.sequence for img in images for target in img.targets}


# --------------------------------------------------------------------------
# Widgets
# --------------------------------------------------------------------------


class TreeCanvas(QGraphicsView):
    """Zoom with the wheel or trackpad, drag to pan, double-click to fit.

    Clicking a branch, tip or collapsed clade reports its display node index."""

    nodeClicked = Signal(int)                 # -1 = empty space
    nodeContextMenu = Signal(int, object)     # index (-1 = empty), global QPoint

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setScene(QGraphicsScene(self))
        self.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        self.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
        self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
        self._press = None

    def node_at(self, view_pos) -> int:
        for item in self.items(view_pos):
            value = item.data(ROLE_NODE)
            if value is not None:
                return int(value)
        return -1

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

    def mousePressEvent(self, event) -> None:  # noqa: N802 (Qt API)
        self._press = event.position().toPoint()
        super().mousePressEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().mouseReleaseEvent(event)
        pos = event.position().toPoint()
        if (event.button() == Qt.MouseButton.LeftButton and self._press is not None
                and (pos - self._press).manhattanLength() < 4):
            self.nodeClicked.emit(self.node_at(pos))
        self._press = None

    def contextMenuEvent(self, event) -> None:  # noqa: N802 (Qt API)
        self.nodeContextMenu.emit(self.node_at(event.pos()), event.globalPos())
        event.accept()

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802 (Qt API)
        self.fit(max_scale=4.0)
        event.accept()


class TreePanel(QWidget):
    """Controls + canvas for one tree file."""

    fileWritten = Signal(object)  # Path of an exported figure or saved tree

    def __init__(self, path: Path | None = None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("TreePanel")
        self._path: Path | None = None
        self._file_model: tl.TreeModel | None = None
        self._model: tl.TreeModel | None = None
        self._state = TreeViewState()
        self._undo: list[tuple[tl.TreeModel, TreeViewState, TreeViewOptions]] = []
        self._display: Display | None = None
        self._selected: tuple[str, ...] | None = None
        self._opts = replace(_last_options)
        self._block = False
        self._loaded_mtime: float | None = None
        self._traits: tt.TraitTable | None = None
        self._trait_dialog = None
        self._trait_job = None
        self._trait_timer = QTimer(self)
        self._trait_timer.setInterval(200)
        self._trait_timer.timeout.connect(self._poll_trait_job)

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
        export_menu = QMenu(export_btn)
        for label, ext in EXPORT_FORMATS:
            export_menu.addAction(f"{label}…", lambda checked=False, e=ext: self._export(e))
        export_btn.setMenu(export_menu)

        self._edit_btn = QPushButton("Edit", self)
        self._edit_menu = QMenu(self._edit_btn)
        self._edit_menu.aboutToShow.connect(self._fill_edit_menu)
        self._edit_btn.setMenu(self._edit_menu)
        self._undo_btn = QPushButton("Undo", self)
        self._undo_btn.setToolTip("Undo the last edit (Cmd+Z / Ctrl+Z).")
        self._undo_btn.clicked.connect(self.undo)
        self._collapse_spin = NoWheelSpinBox(self)
        self._collapse_spin.setRange(0, 100)
        self._collapse_spin.setSpecialValueText("Off")
        self._collapse_spin.setFixedWidth(70)
        self._collapse_spin.editingFinished.connect(self._on_collapse_threshold)
        save_btn = QPushButton("Save tree as…", self)
        save_btn.setToolTip("Write the tree as shown (rerooted, collapsed…) to a new file.")
        save_btn.clicked.connect(self._save_tree_as)

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
        row3 = QHBoxLayout()
        row3.setSpacing(6)
        row3.addWidget(self._edit_btn)
        row3.addWidget(make_help_button(self, "Editing the tree", EDIT_HELP))
        row3.addWidget(self._undo_btn)
        row3.addSpacing(8)
        row3.addWidget(QLabel("Collapse below support:", self))
        row3.addWidget(self._collapse_spin)
        row3.addWidget(make_help_button(self, "Collapse below support", COLLAPSE_HELP))
        row3.addStretch(1)
        row3.addWidget(save_btn)

        self._color_box = NoWheelComboBox(self)
        self._color_box.setMinimumWidth(160)
        self._color_box.currentIndexChanged.connect(self._on_color_by_changed)
        self._traits_btn = QPushButton("Trait table", self)
        traits_menu = QMenu(self._traits_btn)
        traits_menu.addAction("Create or update…", self._open_trait_dialog)
        traits_menu.addAction("Reload", self._reload_traits)
        self._open_traits_action = traits_menu.addAction("Open in spreadsheet app", self._open_traits_file)
        self._traits_btn.setMenu(traits_menu)
        row4 = QHBoxLayout()
        row4.setSpacing(6)
        row4.addWidget(QLabel("Colour tips by:", self))
        row4.addWidget(self._color_box)
        row4.addWidget(make_help_button(self, "Colour tips by", TRAITS_HELP))
        row4.addWidget(self._traits_btn)
        row4.addStretch(1)

        self._canvas = TreeCanvas(self)
        self._canvas.nodeClicked.connect(self._on_node_clicked)
        self._canvas.nodeContextMenu.connect(self._on_context_menu)
        self._status = QLabel(self)
        self._status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        undo_shortcut = QShortcut(QKeySequence(QKeySequence.StandardKey.Undo), self)
        undo_shortcut.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        undo_shortcut.activated.connect(self.undo)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 6, 8, 4)
        layout.setSpacing(4)
        layout.addLayout(row1)
        layout.addLayout(row2)
        layout.addLayout(row3)
        layout.addLayout(row4)
        layout.addWidget(self._canvas, stretch=1)
        layout.addWidget(self._status)

        self._apply_options_to_controls()
        for box in (self._layout_box, self._order_box, self._branch_box):
            box.currentIndexChanged.connect(self._on_controls_changed)
        for check in (self._tips_check, self._scale_check):
            check.toggled.connect(self._on_controls_changed)
        self._font_spin.valueChanged.connect(self._on_controls_changed)
        self._update_undo_button()

        if path is not None:
            self.load_tree(path)

    # ------------------------------------------------------------ state

    def options(self) -> TreeViewOptions:
        return replace(self._opts)

    def state(self) -> TreeViewState:
        return self._state.copy()

    def model(self) -> tl.TreeModel | None:
        return self._model

    def display(self) -> Display | None:
        return self._display

    def path(self) -> Path | None:
        return self._path

    def scene(self) -> QGraphicsScene:
        return self._canvas.scene()

    def canvas(self) -> TreeCanvas:
        return self._canvas

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
            self._collapse_spin.setValue(self._state.collapse_below)
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
        self._save_view_state()

    # ------------------------------------------------------------ loading

    def load_tree(self, path: Path) -> bool:
        new_path = Path(path).expanduser().resolve()
        try:
            mtime = new_path.stat().st_mtime
        except OSError:
            mtime = None
        if (self._model is not None and new_path == self._path and mtime is not None
                and mtime == self._loaded_mtime):
            # Same unchanged file opened again (e.g. the file list refreshed):
            # keep unsaved edits and the undo history.
            self._redraw(refit=False)
            return True
        self._path = new_path
        self._loaded_mtime = mtime
        self._undo.clear()
        self._selected = None
        self._state = TreeViewState()
        try:
            self._file_model = tl.read_tree_file(self._path)
        except (OSError, tl.TreeReadError) as exc:
            self._file_model = self._model = self._display = None
            self._canvas.scene().clear()
            self._status.setText(f"Could not read {self._path.name}: {exc}")
            self._update_undo_button()
            return False
        self._model = self._file_model
        self._load_view_state()
        self._load_traits()
        if not self._model.has_support and self._opts.branch_labels == BRANCH_LABELS_SUPPORT:
            self._opts = replace(self._opts, branch_labels=BRANCH_LABELS_NONE)
        self._apply_options_to_controls()
        self._update_undo_button()
        self._redraw(refit=True)
        return True

    def _load_view_state(self) -> None:
        if self._path is None:
            return
        try:
            data = json.loads(view_state_path(self._path).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            return
        opts = data.get("options") or {}
        known = {k: v for k, v in opts.items() if k in TreeViewOptions.__dataclass_fields__}
        try:
            self._opts = replace(TreeViewOptions(), **known)
        except TypeError:
            pass
        self._state.collapsed = {tuple(k) for k in data.get("collapsed", []) if k}
        self._state.colors = {tuple(c["tips"]): str(c["color"])
                              for c in data.get("colors", []) if c.get("tips") and c.get("color")}
        sub = data.get("subtree")
        self._state.subtree = tuple(sub) if sub else None
        try:
            self._state.collapse_below = max(0, min(100, int(data.get("collapse_below", 0))))
        except (TypeError, ValueError):
            self._state.collapse_below = 0
        self._state.color_by = str(data.get("colour_by") or "")
        ts = data.get("trait_settings")
        self._state.trait_settings = dict(ts) if isinstance(ts, dict) else {}

    def _save_view_state(self) -> None:
        """Remember colours, collapsed clades and view settings for this tree."""
        if self._path is None or self._model is None:
            return
        data = {
            "version": 1,
            "options": asdict(self._opts),
            "collapsed": [list(k) for k in sorted(self._state.collapsed)],
            "colors": [{"tips": list(k), "color": c} for k, c in sorted(self._state.colors.items())],
            "subtree": list(self._state.subtree) if self._state.subtree else None,
            "collapse_below": self._state.collapse_below,
            "colour_by": self._state.color_by,
            "trait_settings": self._state.trait_settings,
        }
        try:
            view_state_path(self._path).write_text(json.dumps(data, indent=1), encoding="utf-8")
        except OSError:
            pass  # read-only folder: the view still works, it is just not remembered

    # ------------------------------------------------------------ drawing

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
        self._display = make_display(self._model, self._opts, self._state)
        self._apply_traits(self._display)
        build_tree_scene(self._canvas.scene(), self._display, self._opts, ink=ink, muted=muted,
                         selected=self._selected)
        if refit:
            self._canvas.fit(max_scale=1.0)
        self._update_status()

    def _update_status(self, message: str = "") -> None:
        m, d = self._model, self._display
        if m is None or d is None:
            return
        parts = [f"{m.n_tips} tips", f"{m.n_internal} internal nodes"]
        if not m.has_lengths:
            parts.append("no branch lengths in the file")
        if m.has_support and self._opts.branch_labels == BRANCH_LABELS_SUPPORT:
            if any("/" in n.support for n in m.nodes):
                parts.append("support: SH-aLRT/UFBoot")
        if d.n_contracted:
            parts.append(f"{d.n_contracted} weak branches collapsed")
        if self._state.subtree:
            parts.append(f"showing one clade ({len(self._state.subtree)} tips)")
        if self._state.edited:
            parts.append(f"tree {self._state.edited}: use Save tree as… to keep it")
        if self._selected:
            parts.append(f"selected: {len(self._selected)} tips")
        if message:
            parts.append(message)
        self._status.setText(" · ".join(parts))

    def showEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().showEvent(event)
        if self._model is not None:
            self._canvas.fit(max_scale=1.0)

    # ------------------------------------------------------------ selection and menus

    def _key_for(self, index: int) -> tuple[str, ...] | None:
        if self._display is None or index < 0 or index >= len(self._display.keys):
            return None
        return self._display.keys[index]

    def _on_node_clicked(self, index: int) -> None:
        key = self._key_for(index)
        self._selected = None if key == self._selected else key
        self._redraw(refit=False)

    def select_clade(self, key: tuple[str, ...] | None) -> None:
        self._selected = tuple(sorted(key)) if key else None
        self._redraw(refit=False)

    def _is_collapsed(self, key) -> bool:
        return key in self._state.collapsed

    def _add_clade_actions(self, menu: QMenu, key: tuple[str, ...] | None) -> None:
        d = self._display
        idx = d.keys.index(key) if (d is not None and key in d.keys) else -1
        node = d.model.nodes[idx] if idx >= 0 else None
        is_root = node is not None and node.parent is None
        is_clade = node is not None and (not node.is_tip or idx in d.info)
        is_internal = node is not None and not node.is_tip
        label = "branch" if node is None or not node.is_tip or idx in d.info else f"“{node.name}”"
        act = menu.addAction(f"Reroot on this {label}", lambda: self.reroot_at(key))
        act.setEnabled(node is not None and not is_root and not self._state.subtree)
        act = menu.addAction("Swap branches", lambda: self.swap_at(key))
        act.setEnabled(is_internal)
        if key is not None and self._is_collapsed(key):
            menu.addAction("Expand clade", lambda: self.toggle_collapse(key))
        else:
            act = menu.addAction("Collapse clade", lambda: self.toggle_collapse(key))
            act.setEnabled(is_internal and not is_root)
        act = menu.addAction("Show only this clade", lambda: self.show_only(key))
        act.setEnabled(is_clade and not is_root)
        act = menu.addAction("Colour clade…", lambda: self._choose_color(key))
        act.setEnabled(node is not None)
        if key is not None and key in self._state.colors:
            menu.addAction("Remove colour", lambda: self.set_clade_color(key, None))
        if node is not None and node.is_tip and idx not in d.info:
            menu.addSeparator()
            menu.addAction("Copy name", lambda: QApplication.clipboard().setText(node.name))
            header = (self._traits.rows.get(node.name, {}).get("header", "")
                      if self._traits else "")
            if header and header != node.name:
                menu.addAction("Copy original header",
                               lambda: QApplication.clipboard().setText(header))

    def _add_tree_actions(self, menu: QMenu) -> None:
        act = menu.addAction("Midpoint root", self.midpoint_root)
        act.setEnabled(self._model is not None and not self._state.subtree)
        act = menu.addAction("Back to full tree", self.full_tree)
        act.setEnabled(bool(self._state.subtree))
        act = menu.addAction("Expand all", self.expand_all)
        act.setEnabled(bool(self._state.collapsed))
        act = menu.addAction("Remove all colours", self.clear_colors)
        act.setEnabled(bool(self._state.colors))
        menu.addSeparator()
        act = menu.addAction("Undo", self.undo)
        act.setEnabled(bool(self._undo))
        menu.addAction("Save tree as…", self._save_tree_as)

    def _fill_edit_menu(self) -> None:
        menu = self._edit_menu
        menu.clear()
        if self._selected:
            title = menu.addAction(f"Selected clade ({len(self._selected)} tips)")
            title.setEnabled(False)
            self._add_clade_actions(menu, self._selected)
        else:
            hint = menu.addAction("Click a branch to select it")
            hint.setEnabled(False)
        menu.addSeparator()
        self._add_tree_actions(menu)

    def _on_context_menu(self, index: int, global_pos) -> None:
        if self._model is None:
            return
        menu = QMenu(self)
        key = self._key_for(index)
        if key is not None:
            self._selected = key
            self._redraw(refit=False)
            self._add_clade_actions(menu, key)
            menu.addSeparator()
        self._add_tree_actions(menu)
        menu.exec(global_pos)

    # ------------------------------------------------------------ edits

    def _push_undo(self) -> None:
        if self._model is not None:
            self._undo.append((self._model, self._state.copy(), replace(self._opts)))
            del self._undo[:-50]
        self._update_undo_button()

    def _update_undo_button(self) -> None:
        self._undo_btn.setEnabled(bool(self._undo))

    def _after_edit(self, message: str = "", *, refit: bool = False) -> None:
        self._update_undo_button()
        self._redraw(refit=refit)
        if message:
            self._update_status(message)
        self._save_view_state()

    def undo(self) -> None:
        if not self._undo:
            return
        self._model, self._state, self._opts = self._undo.pop()
        self._apply_options_to_controls()
        self._after_edit("undone")

    def _node_in_model(self, key) -> int | None:
        if self._model is None or key is None:
            return None
        return te.find_clade(self._model, key)

    def reroot_at(self, key) -> None:
        idx = self._node_in_model(key)
        if idx is None or idx == 0:
            return
        self._push_undo()
        self._model = te.reroot(self._model, idx)
        self._state.edited = "rerooted"
        self._after_edit("rerooted", refit=True)

    def midpoint_root(self) -> None:
        if self._model is None:
            return
        self._push_undo()
        self._model = te.midpoint_root(self._model)
        self._state.edited = "rerooted"
        self._after_edit("midpoint-rooted", refit=True)

    def swap_at(self, key) -> None:
        idx = self._node_in_model(key)
        if idx is None or self._model.nodes[idx].is_tip:
            return
        self._push_undo()
        self._model = te.swap_children(self._model, idx)
        if self._opts.order != tl.ORDER_FILE:
            # a sorted order would undo the swap right away
            self._opts = replace(self._opts, order=tl.ORDER_FILE)
            self._apply_options_to_controls()
        self._state.edited = self._state.edited or "edited"
        self._after_edit("branches swapped")

    def toggle_collapse(self, key) -> None:
        if key is None:
            return
        self._push_undo()
        if key in self._state.collapsed:
            self._state.collapsed.discard(key)
            self._after_edit("expanded")
        else:
            self._state.collapsed.add(key)
            self._after_edit("collapsed")

    def expand_all(self) -> None:
        if not self._state.collapsed:
            return
        self._push_undo()
        self._state.collapsed.clear()
        self._after_edit("all clades expanded")

    def show_only(self, key) -> None:
        if key is None:
            return
        self._push_undo()
        self._state.subtree = key
        self._selected = None
        self._after_edit("showing one clade", refit=True)

    def full_tree(self) -> None:
        if not self._state.subtree:
            return
        self._push_undo()
        self._state.subtree = None
        self._after_edit("full tree", refit=True)

    def set_clade_color(self, key, color: str | None) -> None:
        if key is None:
            return
        self._push_undo()
        if color:
            self._state.colors[key] = QColor(color).name()
        else:
            self._state.colors.pop(key, None)
        self._after_edit("colour set" if color else "colour removed")

    def clear_colors(self) -> None:
        if not self._state.colors:
            return
        self._push_undo()
        self._state.colors.clear()
        self._after_edit("colours removed")

    def _choose_color(self, key) -> None:
        start = QColor(self._state.colors.get(key, "#d62728"))
        color = QColorDialog.getColor(start, self, "Colour clade")
        if color.isValid():
            self.set_clade_color(key, color.name())

    def set_collapse_below(self, value: int) -> None:
        value = max(0, min(100, int(value)))
        if value == self._state.collapse_below:
            return
        self._push_undo()
        self._state.collapse_below = value
        self._block = True
        self._collapse_spin.setValue(value)
        self._block = False
        self._after_edit()

    def _on_collapse_threshold(self) -> None:
        if not self._block:
            self.set_collapse_below(self._collapse_spin.value())

    # ------------------------------------------------------------ traits

    def traits(self) -> tt.TraitTable | None:
        return self._traits

    def _load_traits(self) -> None:
        self._traits = None
        if self._path is not None:
            path = tt.traits_path(self._path)
            if path.is_file():
                try:
                    self._traits = tt.read_trait_table(path)
                except (OSError, ValueError, UnicodeDecodeError):
                    self._traits = None
            # RT-domain trees come with a lineage column: colour by it on first opening.
            if (not self._state.color_by and self._traits is not None
                    and "lineage" in self._traits.columns
                    and not view_state_path(self._path).exists()):
                self._state.color_by = "lineage"
        self._refresh_color_box()

    def _refresh_color_box(self) -> None:
        self._block = True
        try:
            self._color_box.clear()
            self._color_box.addItem("None", "")
            for col in (self._traits.columns if self._traits else []):
                if col not in ("header", "start", "end"):
                    self._color_box.addItem(col, col)
            idx = self._color_box.findData(self._state.color_by)
            if idx < 0:
                self._state.color_by = ""
                idx = 0
            self._color_box.setCurrentIndex(idx)
            self._color_box.setEnabled(self._traits is not None)
            self._open_traits_action.setEnabled(self._traits is not None)
        finally:
            self._block = False

    def set_color_by(self, column: str) -> None:
        idx = self._color_box.findData(column)
        if idx >= 0:
            self._color_box.setCurrentIndex(idx)

    def _on_color_by_changed(self, *_args) -> None:
        if self._block:
            return
        self._state.color_by = str(self._color_box.currentData() or "")
        self._redraw(refit=False)
        self._save_view_state()

    def _apply_traits(self, display: Display) -> None:
        if self._traits is None:
            return
        m = display.model
        tips = [t for t in m.tips() if t not in display.info]
        names = {t: m.nodes[t].name for t in tips}
        cols = [c for c in self._traits.columns if c != "header"]
        for t, name in names.items():
            row = self._traits.rows.get(name)
            if row:
                values = [f"{c}: {row.get(c, '')}" for c in cols if row.get(c, "")]
                if values:
                    display.tip_notes[t] = "\n".join(values)
        column = self._state.color_by
        if column and column in self._traits.columns:
            scale = tt.colour_scale(self._traits, column, list(names.values()))
            display.tip_colors = {t: scale.colors[n] for t, n in names.items() if n in scale.colors}
            display.legend = scale

    def _reload_traits(self) -> None:
        self._load_traits()
        self._redraw(refit=False)
        self._update_status("trait table reloaded" if self._traits else "no trait table yet")

    def _open_traits_file(self) -> None:
        if self._path is not None and tt.traits_path(self._path).is_file():
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(tt.traits_path(self._path))))

    def _open_trait_dialog(self) -> None:
        from repeatviewer.views.tree_traits_dialog import TraitTableDialog

        if self._model is None or self._path is None:
            return
        if self._trait_dialog is not None and self._trait_dialog.isVisible():
            self._trait_dialog.raise_()
            return
        dialog = TraitTableDialog(settings=self._state.trait_settings, tree_path=self._path,
                                  parent=self)
        dialog.setModal(False)
        dialog.finished.connect(lambda code, d=dialog: self._on_trait_dialog_finished(d, code))
        self._trait_dialog = dialog
        dialog.show()

    def trait_dialog_open(self) -> bool:
        return self._trait_dialog is not None and self._trait_dialog.isVisible()

    def apply_path_from_panel(self, path: Path) -> bool:
        """A file clicked in the Files panel fills the open trait dialog."""
        if self.trait_dialog_open():
            if self._trait_dialog.apply_path_from_panel(path):
                self._trait_dialog.raise_()
                return True
        return False

    def _on_trait_dialog_finished(self, dialog, code: int) -> None:
        self._trait_dialog = None
        if code != dialog.DialogCode.Accepted:
            return
        self._state.trait_settings = dialog.settings()
        self._save_view_state()
        if dialog.should_run():
            self.build_traits()

    def build_traits(self, *, wait: bool = False) -> None:
        """Create or update the trait table in the background."""
        import concurrent.futures

        if self._model is None or self._path is None or self._trait_job is not None:
            return
        tips = [self._model.nodes[t].name for t in self._model.tips()]
        settings = dict(self._state.trait_settings)
        meth = settings.get("methylation_file") or None
        gap = int(settings.get("array_gap") or tt.DEFAULT_ARRAY_GAP)
        tree_path = self._path
        existing = self._traits
        fish = _fish_files(tree_path.parent)

        def work():
            table = tt.build_trait_table(tree_path, tips, methylation_file=Path(meth) if meth else None,
                                         array_gap=gap, fish_files=fish, existing=existing)
            return tt.write_trait_table(tt.traits_path(tree_path), table), bool(meth)

        pool = concurrent.futures.ThreadPoolExecutor(max_workers=1)
        self._trait_job = (pool, pool.submit(work))
        self._update_status("building the trait table…")
        if wait:
            self._trait_job[1].result()
            self._poll_trait_job()
        else:
            self._trait_timer.start()

    def _poll_trait_job(self) -> None:
        if self._trait_job is None:
            self._trait_timer.stop()
            return
        pool, future = self._trait_job
        if not future.done():
            return
        self._trait_timer.stop()
        self._trait_job = None
        pool.shutdown(wait=False)
        try:
            path, with_meth = future.result()
        except Exception as exc:  # noqa: BLE001 - shown to the user
            QMessageBox.warning(self, "Trait table", f"Could not build the trait table:\n{exc}")
            self._update_status("trait table failed")
            return
        self._load_traits()
        if not self._state.color_by and self._traits:
            for col in (("methylation",) if with_meth else ()) + ("array", "chromosome", "source_file"):
                if col in self._traits.columns:
                    self._state.color_by = col
                    break
            self._refresh_color_box()
        self._redraw(refit=False)
        self._save_view_state()
        self._update_status(f"saved {path.name}")
        self.fileWritten.emit(path)

    # ------------------------------------------------------------ saving

    def tree_to_save(self) -> tl.TreeModel | None:
        """The tree as shown: rerooted / swapped, weak branches collapsed, one
        clade if only one is shown, in the chosen order (collapsed triangles are
        a drawing aid and keep all their tips)."""
        if self._model is None:
            return None
        m, _n = tree_structure(self._model, self._state)
        return tl.ordered(m, self._opts.order)

    def save_tree(self, out: Path) -> Path:
        tree = self.tree_to_save()
        if tree is None:
            raise ValueError("There is no tree to save.")
        out.write_text(te.to_newick(tree) + "\n", encoding="utf-8")
        return out

    def _save_tree_as(self) -> None:
        if self._model is None or self._path is None:
            QMessageBox.information(self, "Save tree", "There is no tree to save.")
            return
        edited = self._state.edited
        if self._state.subtree:
            edited = f"{edited}_clade" if edited else "clade"
        if self._state.collapse_below and not self._state.subtree:
            edited = f"{edited}_collapsed" if edited else "collapsed"
        default = saved_tree_default_path(self._path, edited)
        chosen, _ = QFileDialog.getSaveFileName(
            self, "Save tree as", str(default), "Newick tree (*.treefile *.nwk *.tree)"
        )
        if not chosen:
            return
        out = Path(chosen)
        if out.resolve() == self._path:
            QMessageBox.warning(self, "Save tree",
                                "Choose a new file name: the original tree is never overwritten.")
            return
        if not tl.is_tree_file(out):
            out = out.with_name(out.name + ".treefile")
        try:
            self.save_tree(out)
        except (OSError, ValueError) as exc:
            QMessageBox.warning(self, "Save tree", f"Could not write {out.name}:\n{exc}")
            return
        self._update_status(f"saved {out.name}")
        self.fileWritten.emit(out)

    # ------------------------------------------------------------ export

    def _export(self, ext: str) -> None:
        if self._display is None or self._path is None:
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
            export_tree(self._display, self._opts, out, title=self._path.name)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Export tree", f"Could not write {out.name}:\n{exc}")
            return
        self._update_status(f"saved {out.name}")
        self.fileWritten.emit(out)
