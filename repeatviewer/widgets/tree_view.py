"""Tree tab of the Preview panel: an interactive view of a Newick tree.

Opens ``.treefile`` (IQ-TREE), ``.nwk``, ``.newick``, ``.tre`` and ``.tree``.

* View: rectangular, circular and radial layouts; tip names, support or
  branch-length labels, font size, scale bar, zoom and pan. Display settings
  live in the right Display rail (same chrome as the sequence Annotation view).
* Edit (right-click a branch, or the Edit menu for the selected branch):
  reroot here, midpoint root, swap branches, collapse / expand a clade,
  colour a clade, show one clade only; collapse all branches below a support
  value; undo (Cmd+Z / Ctrl+Z) and redo (Cmd+Shift+Z / Ctrl+Shift+Z).
* Keep: colours, collapsed clades and view settings are remembered in a hidden
  file next to the tree (``.<tree name>.view.json``). Rerooting and swapping
  change the tree itself, which is never overwritten: Save → "Tree"
  writes a new file (``<tree name>_rerooted.treefile``, rule 1).
* Save ▾: figure as shown (SVG, PDF, PNG), tree (.treefile), distance matrix,
  copy Newick. Tree runs no longer write figures; they come from here.

The geometry and the edits come from the pure, unit-tested modules
:mod:`repeatviewer.backends.tree_layout` and :mod:`repeatviewer.backends.tree_edit`;
this module only turns them into a ``QGraphicsScene``.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field, replace
from pathlib import Path

from PySide6.QtCore import (
    QEvent,
    QLineF,
    QPoint,
    QPointF,
    QRect,
    QRectF,
    QSettings,
    QSizeF,
    Qt,
    QTimer,
    QUrl,
    Signal,
)
from PySide6.QtGui import (
    QBrush,
    QDesktopServices,
    QGuiApplication,
    QLinearGradient,
    QColor,
    QFont,
    QFontMetricsF,
    QImage,
    QKeySequence,
    QMouseEvent,
    QPageSize,
    QPainter,
    QPainterPath,
    QPainterPathStroker,
    QPalette,
    QPen,
    QPolygonF,
    QTransform,
)
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QColorDialog,
    QComboBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QGraphicsEllipseItem,
    QGraphicsPathItem,
    QGraphicsPolygonItem,
    QGraphicsRectItem,
    QGraphicsScene,
    QGraphicsSimpleTextItem,
    QGraphicsView,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMenu,
    QMessageBox,
    QPushButton,
    QRubberBand,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import tree_edit as te
from repeatviewer.backends import tree_find as tf
from repeatviewer.backends import tree_layout as tl
from repeatviewer.backends import tree_traits as tt
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.preview_controls import (
    ANNOTATION_FIT_BUTTON_RENDERED_WIDTH,
    ANNOTATION_ZOOM_SLIDER_WIDTH,
    OVERVIEW_SELECT_FIT_HEIGHT,
    OVERVIEW_TOOLBAR_HEIGHT,
    ZOOM_SLIDER_HEIGHT,
    apply_annotation_zoom_slider_style,
    make_overview_zoom_button,
)
from repeatviewer.widgets.sequence_annotation_rail import (
    _h_divider,
    _page_title_with_rule,
    _section_label,
)
from repeatviewer.widgets.spin_boxes import NoWheelDoubleSpinBox, NoWheelSlider, NoWheelSpinBox
from repeatviewer.widgets.tree_display_rail import TreeDisplayRail

BRANCH_LABELS_NONE = "none"
BRANCH_LABELS_SUPPORT = "support"
BRANCH_LABELS_LENGTH = "length"
TIP_LABEL_NAME = "name"

EXPORT_FORMATS = (("SVG", "svg"), ("PDF", "pdf"), ("PNG", "png"))
PNG_SCALE = 4.0            # 4 × 72 dpi ≈ 300 dpi
_ZOOM_STEP = 1.12
_TRACKPAD_ZOOM_BASE = 1.0022
_MIN_TREE_PX = 260.0       # smallest drawing size (rectangular width / radius)
_RECT_TREE_WIDTH_PX = 520.0
SPREAD_MIN, SPREAD_MAX = 100, 500   # % (Spread)
_MARGIN_PX = 24.0
_HIT_WIDTH = 7.0           # clickable width around a branch, in pixels
_PICK_RADIUS_PX = 8        # a click reaches the nearest branch / tip within this many
                           # screen pixels, whatever the zoom
_SELECTED = QColor(0, 110, 230)
ROLE_NODE = 0              # QGraphicsItem.data key: display node index

TIPS_HELP = (
    "<b>Align</b>: all tip names on one line (rectangular) or one circle (circular), "
    "joined to their branches by dotted lines, so they never sit on top of each other. "
    "Not for the radial layout.<br>"
    "<b>Text</b>: the name, or the name with a column of the trait table next to it "
    "(chromosome, array, lineage…), or that column instead of the name. The trait "
    "table is made in Colour tips → Trait table.<br>"
    "<b>Max length</b>: shorten long names to this many characters, keeping their end "
    "(where the coordinates are) or their start. Off = full names. Hover over a tip for "
    "the full name."
)
SCALE_HELP = (
    "<b>Length</b>: the bar's length in substitutions per site. Auto picks a round value "
    "(1, 2 or 5 × 10ⁿ) about a fifth of the tree's depth.<br>"
    "The scale bar is drawn for proportional branch lengths only."
)
SUPPORT_HELP = (
    "<b>Show: Support</b>: the label on each internal node.<br>"
    "IQ-TREE shows <b>SH-aLRT/UFBoot</b>, e.g. 92/97. A branch is well supported "
    "at SH-aLRT 80 or more and UFBoot 95 or more.<br>"
    "FastME bootstrap: 70 or more. FastTree: local support from 0 to 1, close to 1 is strong.<br>"
    "<b>Branch length</b>: substitutions per site.<br>"
    "<b>Hide crowded</b>: a label is shown only if it fits on its branch and touches no "
    "line, tip name or other label (labels on longer branches are kept first). Turn it "
    "off to see every label. Hover over a branch to read its length and support.<br>"
    "<b>Hide zeros</b>: no label on branches whose length is 0 or rounds to 0.<br>"
    "<b>Decimals</b>: places after the point for branch lengths (0.0277 shows as 0.03 "
    "with 2). Hover over a tip name for the length with more digits. "
    "<b>Font</b>: Auto = 80 % of the tip-name font."
)
DRAW_HELP = (
    "<b>Order</b>: sort each node's branches by clade size (ladderize).<br>"
    "<b>Branches</b>: <b>Proportional</b> draws the lengths from the file. "
    "<b>Equal</b> ignores lengths and shows the branching order only, with all tips "
    "lined up (a cladogram). <b>Square root</b> keeps longer branches longer but shrinks "
    "one very long branch, so it no longer squeezes the rest of the tree. "
    "The scale bar is shown for Proportional only.<br>"
    "<b>Spread</b>: makes the tree larger while the font stays the same; the circular "
    "layout also opens a space in the centre. Use it when names or branches crowd.<br>"
    "<b>Root angle</b>: turns a circular or radial tree, e.g. to put a long branch "
    "where it fits the figure. <b>Flip</b>: mirrors the tree left to right."
)
COLLAPSE_HELP = (
    "<b>Collapse below support</b>: branches with support under this value are merged "
    "into their parent, so only well-supported groupings remain.<br>"
    "IQ-TREE: the UFBoot value is used (the second number, e.g. 97 in 92/97). "
    "FastTree values (0 to 1) are read as percentages.<br>"
    "0 = off. The tree file is not changed; use <b>Save → Tree</b> to keep the result.<br>"
    "<b>Distance below</b>: clades whose tips all lie closer than this to the clade's own "
    "node (substitutions per site) are drawn as triangles, e.g. near-identical copies. "
    "Right-click a triangle to open it. Off = 0."
)
EDIT_HELP = (
    "<b>Right-click a branch</b> (or select it and use <b>Edit</b>) to reroot there, swap "
    "its branches, collapse or colour its clade, or show that clade only.<br>"
    "<b>Extract sequences</b> (same menu) saves the copies of that clade to a new FASTA "
    "next to the tree: <b>Without gaps</b> (your original records, full headers) or "
    "<b>As aligned</b> (the tree's alignment, gaps kept).<br>"
    "<b>Double-click a tip</b> (or right-click it → Show in alignment) to open the tree's "
    "alignment in a Preview tab with that sequence selected; right-click a clade → "
    "Show clade in alignment selects all of its sequences. There, right-click → "
    "Extract Sequence or Copy Sequence.<br>"
    "Undo: Cmd+Z / Ctrl+Z. Redo: Cmd+Shift+Z / Ctrl+Shift+Z.<br>"
    "Colours and collapsed clades are remembered for this tree. Rerooting and swapping "
    "change the tree: the file is never overwritten, use <b>Save → Tree</b>."
)


FIND_HELP = (
    "<b>Find</b> tips by name or by the values in the trait table (chromosome, array, "
    "lineage, any column you added).<br>"
    "Type one or several words (separated by spaces or commas): a tip matches when it "
    "contains <b>any</b> of them.<br>"
    "<b>Find next / previous</b> select one tip at a time and scroll to it; a tip inside a "
    "collapsed clade selects that clade.<br>"
    "<b>Find all</b> selects every match, so Edit, right-click and Save work on them "
    "(for example Extract sequences)."
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
    align_labels: bool = False       # tip names on one line / circle, dotted lines to the tips
    branch_lengths: str = tl.BRANCHES_PROPORTIONAL
    spread: int = 100                # % : tree size relative to the font (circular: open centre)
    hide_crowded: bool = True        # skip branch labels that would touch lines, names, labels
    hide_zeros: bool = True          # no length label when it is 0 or rounds to 0
    root_angle: int = 0              # circular / radial: turn the tree (degrees)
    flip: bool = False               # mirror left-right
    line_weight: float = 1.2
    tip_label: str = TIP_LABEL_NAME  # "name", "name+col:<column>" or "col:<column>"
    max_chars: int = 0               # shorten tip names to this many characters (0 = off)
    keep_end: bool = True            # shortened names keep their end (coordinates)
    branch_decimals: int = 2         # places after the point for branch-length labels
    branch_font_pt: int = 0          # 0 = 80 % of the tip font
    scale_length: float = 0.0        # 0 = automatic
    scale_font_pt: int = 0           # 0 = like the branch labels


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
    collapse_distance: float = 0.0   # collapse clades whose tips lie this close (0 = off)
    expanded: set[tuple[str, ...]] = field(default_factory=set)  # opened although within distance

    def copy(self) -> "TreeViewState":
        return TreeViewState(set(self.collapsed), dict(self.colors), self.subtree,
                             self.collapse_below, self.edited, self.color_by,
                             dict(self.trait_settings), self.collapse_distance,
                             set(self.expanded))


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
    tip_text: dict[int, str] = field(default_factory=dict)     # label when not the plain name
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
    if state.collapse_distance > 0:
        sets = te.tip_sets(m)
        for idx in te.clades_within(m, state.collapse_distance):
            if tuple(sorted(sets[idx])) not in state.expanded:
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


def spread_factor(opts: TreeViewOptions) -> float:
    return max(SPREAD_MIN, min(SPREAD_MAX, int(opts.spread))) / 100.0


def inner_fraction(opts: TreeViewOptions) -> float:
    """Circular layout: empty centre (share of the tree depth) growing with Spread."""
    return min(0.5, 0.15 * (spread_factor(opts) - 1.0))


def _pixel_scale(layout: tl.Layout, n_tips: int, row_h: float,
                 spread: float = 1.0) -> tuple[float, float]:
    """(sx, sy) from layout units to pixels."""
    extent = layout.extent if layout.extent > 0 else 1.0
    if layout.kind == tl.LAYOUT_RECTANGULAR:
        return _RECT_TREE_WIDTH_PX * spread / extent, row_h
    if layout.kind == tl.LAYOUT_CIRCULAR:
        radius = max(_MIN_TREE_PX, n_tips * row_h / (2 * math.pi))
    else:
        radius = max(_MIN_TREE_PX, math.sqrt(n_tips) * row_h * 2.0)
    s = radius * spread / extent
    return s, s


def display_leaves(display: Display) -> list[int]:
    """Display nodes drawn as a leaf (a tip or a collapsed clade), top to bottom."""
    nodes = display.model.nodes
    out: list[int] = []
    stack = [0]
    while stack:
        i = stack.pop()
        node = nodes[i]
        if node.is_tip or i in display.info:
            out.append(i)
        else:
            stack.extend(reversed(node.children))
    return out


def selection_branch_nodes(
    display: Display, selected: tuple[str, ...] | None
) -> frozenset[int]:
    """Display nodes whose branch should be blue for ``selected`` tips.

    Tips that sit next to each other in the drawing form one group, drawn as the
    spanning subtree (every branch on a path between them, including the LCA
    connector). Groups separated by unselected tips are drawn on their own, so
    two distant picks do not light up the branches between them.
    """
    if not selected:
        return frozenset()
    sel = set(selected)
    keys = display.keys
    runs: list[list[int]] = []
    current: list[int] = []
    for i in display_leaves(display):
        if any(tip in sel for tip in keys[i]):
            current.append(i)
        elif current:
            runs.append(current)
            current = []
    if current:
        runs.append(current)
    out: set[int] = set()
    for run in runs:
        out |= _run_branch_nodes(display, run, sel)
    return frozenset(out)


def _run_branch_nodes(display: Display, run: list[int], sel: set[str]) -> set[int]:
    nodes = display.model.nodes
    keys = display.keys
    run_sel = {tip for i in run for tip in keys[i] if tip in sel}
    n_sel = len(run_sel)
    counts = [0] * len(nodes)
    for i in run:
        counts[i] = sum(1 for tip in keys[i] if tip in run_sel)
    for node in reversed(nodes):
        if node.children:
            counts[node.index] = sum(counts[c] for c in node.children)
    out: set[int] = set()
    for i, node in enumerate(nodes):
        under = counts[i]
        if under <= 0:
            continue
        if under < n_sel:
            # Proper edge on a path between selected tips (includes tip stubs).
            out.add(i)
            continue
        # under == n_sel: highlight only the LCA (no child holds every selected tip).
        if node.is_tip:
            if node.name in run_sel or i in display.info:
                out.add(i)
        elif node.children and all(counts[c] < n_sel for c in node.children):
            out.add(i)
    return out


def shorten(text: str, max_chars: int, keep_end: bool = True) -> str:
    """Tip name cut to ``max_chars`` (0 = unchanged), keeping its end or start."""
    if max_chars <= 0 or len(text) <= max_chars:
        return text
    return "…" + text[-max_chars:] if keep_end else text[:max_chars] + "…"


class _Obstacles:
    """Lines and boxes branch labels must not touch; a coarse grid keeps lookups fast."""

    def __init__(self, cell: float) -> None:
        self._cell = max(4.0, cell)
        self._lines: dict[tuple[int, int], list[QLineF]] = {}
        self._shapes: dict[tuple[int, int], list[QPolygonF]] = {}

    def _cells(self, r: QRectF):
        c = self._cell
        for gx in range(int(math.floor(r.left() / c)), int(math.floor(r.right() / c)) + 1):
            for gy in range(int(math.floor(r.top() / c)), int(math.floor(r.bottom() / c)) + 1):
                yield gx, gy

    def add_line(self, a: QPointF, b: QPointF) -> None:
        # Long lines are cut into pieces about one cell long, so each sits in few cells.
        n = max(1, int(math.ceil(QLineF(a, b).length() / self._cell)))
        for k in range(n):
            p = a + (b - a) * (k / n)
            q = a + (b - a) * ((k + 1) / n)
            piece = QLineF(p, q)
            for key in self._cells(QRectF(p, q).normalized()):
                self._lines.setdefault(key, []).append(piece)

    def add_path(self, path: QPainterPath) -> None:
        for poly in path.toSubpathPolygons():
            for k in range(len(poly) - 1):
                if poly[k] != poly[k + 1]:
                    self.add_line(poly[k], poly[k + 1])

    def add_shape(self, poly: QPolygonF) -> None:
        for key in self._cells(poly.boundingRect()):
            self._shapes.setdefault(key, []).append(poly)

    def hits(self, poly: QPolygonF) -> bool:
        edges = [QLineF(poly[k], poly[(k + 1) % len(poly)]) for k in range(len(poly))]
        for key in self._cells(poly.boundingRect()):
            for line in self._lines.get(key, ()):
                if poly.containsPoint(line.p1(), Qt.FillRule.OddEvenFill):
                    return True
                if any(e.intersects(line)[0] == QLineF.IntersectionType.BoundedIntersection
                       for e in edges):
                    return True
            for other in self._shapes.get(key, ()):
                if other.intersects(poly):
                    return True
        return False


def _is_zero_text(text: str) -> bool:
    """True for "0", "0.00", "-0.00" and the like."""
    try:
        return float(text) == 0.0
    except ValueError:
        return False


def _branch_tooltip(node: tl.Node, n_seqs: int, has_lengths: bool) -> str:
    lines = []
    if has_lengths:
        lines.append(f"Branch length: {tl.branch_length_text(node.length, 4)}")
    if node.support:
        lines.append(f"Support: {node.support}")
    if n_seqs > 1:
        lines.append(f"{n_seqs} sequences below")
    return "\n".join(lines)


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
    mode = opts.branch_lengths if opts.branch_lengths in tl.BRANCH_TRANSFORMS else tl.BRANCHES_PROPORTIONAL
    geom = tl.transform_lengths(shown, mode)     # drawn lengths; labels keep the file's
    extra_depth = {i: tl.transform_length(c.depth, mode) for i, c in info.items()}
    layout = tl.compute_layout(geom, opts.layout, extra_depth,
                               inner_fraction=inner_fraction(opts))
    layout = tl.orient(layout, opts.root_angle, opts.flip)
    rect_dir = -1.0 if layout.mirrored else 1.0     # rectangular: which way the tips point
    font = QFont()
    font.setPointSizeF(float(opts.font_pt))
    small = QFont(font)
    small.setPointSizeF(float(opts.branch_font_pt) if opts.branch_font_pt > 0
                        else max(5.0, opts.font_pt * 0.8))
    scale_font = QFont(font)
    scale_font.setPointSizeF(float(opts.scale_font_pt) if opts.scale_font_pt > 0
                             else small.pointSizeF())
    decimals = max(1, min(8, int(opts.branch_decimals)))
    fm = QFontMetricsF(font)
    fm_small = QFontMetricsF(small)
    row_h = fm.height() * 1.3
    sx, sy = _pixel_scale(layout, shown.n_tips, row_h, spread_factor(opts))
    radial = layout.kind != tl.LAYOUT_RECTANGULAR
    align = opts.align_labels and layout.kind != tl.LAYOUT_RADIAL

    def px(point: tuple[float, float]) -> QPointF:
        return QPointF(point[0] * sx, point[1] * sy)

    def color_for(i: int) -> QColor:
        c = display.colors[i]
        return QColor(c) if c else QColor(ink)

    selected_branches = (
        frozenset()
        if for_export or selected is None
        else selection_branch_nodes(display, selected)
    )

    def pen_for(i: int, width: float | None = None) -> QPen:
        width = opts.line_weight if width is None else width
        is_sel = i in selected_branches
        pen = QPen(_SELECTED if is_sel else color_for(i), max(2.6, width + 1.4) if is_sel else width)
        pen.setCosmetic(not for_export)
        pen.setCapStyle(Qt.PenCapStyle.FlatCap)
        return pen

    # Lines, names and labels already drawn: branch labels must not touch them.
    obstacles = _Obstacles(fm_small.height() * 2.0)
    n_below = [0] * len(shown.nodes)        # sequences under each node (tooltips)
    for node in reversed(shown.nodes):
        if node.is_tip:
            n_below[node.index] = info[node.index].n_tips if node.index in info else 1
        if node.parent is not None:
            n_below[node.parent] += n_below[node.index]

    # --- branches: one item per node (its own branch + its connector) -------
    depth = [d + layout.inner for d in tl.depths(geom)]
    for node in shown.nodes:
        i = node.index
        path = QPainterPath()
        if opts.layout == tl.LAYOUT_CIRCULAR:
            a = math.radians(layout.node_angle[i])
            if node.parent is not None:
                r0 = depth[node.parent] * sx
                path.moveTo(r0 * math.cos(a), -r0 * math.sin(a))
                path.lineTo(px(layout.xy[i]))
            if not node.is_tip and depth[i] > 0 and node.children:
                ends = (layout.node_angle[node.children[0]], layout.node_angle[node.children[-1]])
                a0, a1 = min(ends), max(ends)
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
        if not for_export:
            tip_text = _branch_tooltip(node, n_below[i], shown.has_lengths)
            if tip_text:
                item.setToolTip(tip_text)
        scene.addItem(item)
        obstacles.add_path(path)

    # --- collapsed clades: triangles ------------------------------------------
    for i, ci in info.items():
        base = px(layout.xy[i])
        length = extra_depth[i] * sx
        half = row_h * 0.42
        if radial:
            a = math.radians(layout.node_angle.get(i, 0.0))
            dx, dy = math.cos(a), -math.sin(a)
        else:
            dx, dy = rect_dir, 0.0
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
        obstacles.add_shape(tri)

    # --- where each tip's dot and name start ----------------------------------------
    # At the end of its branch (or collapsed triangle), or, with aligned names, on
    # one line / circle for all tips, joined to the branch end by a dotted line.
    def tip_dir(tip: int) -> tuple[float, float]:
        if radial:
            a = math.radians(layout.tip_angle[tip])
            return math.cos(a), -math.sin(a)
        return rect_dir, 0.0

    anchors: dict[int, QPointF] = {}
    leader_pen = QPen(muted, 0.8, Qt.PenStyle.DotLine)
    leader_pen.setCosmetic(not for_export)
    for tip in shown.tips():
        p = px(layout.xy[tip])
        dx, dy = tip_dir(tip)
        extra = extra_depth.get(tip, 0.0) * sx
        end = QPointF(p.x() + dx * extra, p.y() + dy * extra)
        if align:
            if layout.kind == tl.LAYOUT_RECTANGULAR:
                target = QPointF(rect_dir * layout.extent * sx, p.y())
            else:
                r = layout.extent * sx
                target = QPointF(dx * r, dy * r)
            if math.hypot(target.x() - end.x(), target.y() - end.y()) > 2.0:
                leader = QPainterPath(end)
                leader.lineTo(target)
                scene.addPath(leader, leader_pen)
                obstacles.add_path(leader)
            end = target
        anchors[tip] = end

    # --- trait dots at the tips ---------------------------------------------------
    gap = fm.averageCharWidth() * 0.6
    dot = fm.height() * 0.62 if display.tip_colors else 0.0
    for tip, color in display.tip_colors.items():
        p = anchors[tip]
        dx, dy = tip_dir(tip)
        reach = gap + dot / 2.0
        cx, cy = p.x() + dx * reach, p.y() + dy * reach
        circle = scene.addEllipse(QRectF(cx - dot / 2.0, cy - dot / 2.0, dot, dot),
                                  QPen(QColor(color).darker(130), 0.6), QBrush(QColor(color)))
        circle.setData(ROLE_NODE, tip)
        obstacles.add_shape(QPolygonF(circle.rect()))
        if not for_export and tip in display.tip_notes:
            circle.setToolTip(display.tip_notes[tip])
    if dot:
        gap += dot + gap

    # --- tip labels -------------------------------------------------------------
    if opts.tip_labels:
        for tip in shown.tips():
            node = shown.nodes[tip]
            if tip in info:        # collapsed clade: shorten the name, keep "+ n more"
                ci = info[tip]
                text = (shorten(ci.first_tip, opts.max_chars, opts.keep_end)
                        + f" + {ci.n_tips - 1} more")
            else:
                text = shorten(display.tip_text.get(tip, node.name), opts.max_chars,
                               opts.keep_end)
            tip_sel = (
                selected is not None
                and not for_export
                and (display.keys[tip] == selected or node.name in selected)
            )
            item = QGraphicsSimpleTextItem(text)
            tip_font = QFont(font)
            if tip_sel:
                tip_font.setBold(True)
            item.setFont(tip_font)
            item.setBrush(QBrush(_SELECTED if tip_sel else color_for(tip)))
            item.setData(ROLE_NODE, tip)
            if not for_export:
                tip_text = f"{node.name}\nBranch length: {tl.branch_length_text(node.length, 4)}"
                if tip in display.tip_notes:
                    tip_text += "\n" + display.tip_notes[tip]
                item.setToolTip(tip_text)
            w = item.boundingRect().width()
            h = item.boundingRect().height()
            p = anchors[tip]
            if layout.kind == tl.LAYOUT_RECTANGULAR:
                x = p.x() + gap if rect_dir > 0 else p.x() - gap - w
                item.setPos(x, p.y() - h / 2.0)
            else:
                angle = layout.tip_angle[tip]
                rad = math.radians(angle)
                reach = gap
                anchor = QPointF(p.x() + reach * math.cos(rad), p.y() - reach * math.sin(rad))
                rotation, flipped = tl.label_rotation(angle)
                if flipped:
                    item.setTransformOriginPoint(w, h / 2.0)
                    item.setPos(anchor.x() - w, anchor.y() - h / 2.0)
                else:
                    item.setTransformOriginPoint(0.0, h / 2.0)
                    item.setPos(anchor.x(), anchor.y() - h / 2.0)
                item.setRotation(rotation)
            if tip_sel:
                # Soft blue wash behind the name (follows rotation as a child).
                pad = 2.5
                wash = QColor(_SELECTED)
                wash.setAlphaF(0.18)
                bg = QGraphicsRectItem(
                    item.boundingRect().adjusted(-pad, -pad, pad, pad), item
                )
                bg.setBrush(QBrush(wash))
                bg.setPen(QPen(Qt.PenStyle.NoPen))
                bg.setZValue(-1)
            scene.addItem(item)
            obstacles.add_shape(item.mapToScene(item.boundingRect()))

    # --- branch labels ------------------------------------------------------------
    # Each label runs along its branch (turned so it never reads upside down), 2 px
    # off the line: lengths in the middle, support values at the node end. With
    # "Hide crowded" a label must fit on its branch and touch no line, tip name or
    # other label; it tries the other side of the branch before giving up.
    if opts.branch_labels != BRANCH_LABELS_NONE:
        label_brush = QBrush(muted)
        support_mode = opts.branch_labels == BRANCH_LABELS_SUPPORT
        ascent, cap = fm_small.ascent(), fm_small.capHeight()
        label_gap = 2.0
        candidates: list[tuple[float, int, str]] = []
        for node in shown.nodes:
            if node.parent is None:
                continue
            if support_mode:
                text = node.support
            else:
                text = tl.branch_length_text(node.length, decimals=decimals)
                if opts.hide_zeros and _is_zero_text(text):
                    continue
            if not text:
                continue
            branch_px = max(0.0, geom.nodes[node.index].length) * sx
            candidates.append((branch_px, node.index, text))
        candidates.sort(key=lambda c: -c[0])     # longest branches first

        for branch_px, i, text in candidates:
            item = QGraphicsSimpleTextItem(text)
            item.setFont(small)
            item.setBrush(label_brush)
            item.setData(ROLE_NODE, i)
            w = item.boundingRect().width()
            if opts.hide_crowded and w + 4.0 > branch_px:
                continue                       # wider than its branch
            if layout.kind == tl.LAYOUT_RECTANGULAR:
                u = QPointF(rect_dir, 0.0)
                rotation = 0.0
            else:
                a = math.radians(layout.node_angle[i])
                u = QPointF(math.cos(a), -math.sin(a))   # parent -> node on screen
                rotation, _flipped = tl.label_rotation(layout.node_angle[i])
            if support_mode:
                centre = px(layout.xy[i]) - u * (w / 2.0 + 3.0)
            else:
                centre = px(layout.branch_mid[i])
            for side in (1, -1):               # above the line, then below
                t = QTransform()
                t.translate(centre.x(), centre.y())
                t.rotate(rotation)
                if side > 0:
                    t.translate(-w / 2.0, -label_gap - ascent)
                else:
                    t.translate(-w / 2.0, label_gap - (ascent - cap))
                ink_box = QRectF(0.0, ascent - cap, w, cap).adjusted(-1.0, -1.0, 1.0, 1.0)
                poly = t.map(QPolygonF(ink_box))
                if opts.hide_crowded and obstacles.hits(poly):
                    continue
                item.setTransform(t)
                if opts.hide_crowded:
                    obstacles.add_shape(poly)
                scene.addItem(item)
                break

    # --- legend for the trait colours ------------------------------------------------
    if display.legend is not None:
        _draw_legend(scene, display.legend, font, small, ink, muted, for_export)

    # --- scale bar ------------------------------------------------------------------
    if opts.scale_bar and shown.has_lengths and mode == tl.BRANCHES_PROPORTIONAL:
        value = (opts.scale_length if opts.scale_length > 0
                 else tl.nice_scale_length(layout.extent - layout.inner))
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
            label.setFont(scale_font)
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
                title: str = "", png_dpi: int | None = None) -> Path:
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
        dpi = int(png_dpi) if png_dpi and png_dpi > 0 else int(round(72 * PNG_SCALE))
        scale = max(1.0, float(dpi) / 72.0)
        image = QImage(int(math.ceil(rect.width() * scale)),
                       int(math.ceil(rect.height() * scale)), QImage.Format.Format_ARGB32)
        image.fill(QColor(255, 255, 255))
        image.setDotsPerMeterX(int(dpi / 0.0254))
        image.setDotsPerMeterY(int(dpi / 0.0254))
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


class _Section(QWidget):
    """Collapsible Display section (MSA Tracks font; arrow left, ? right-aligned)."""

    def __init__(self, title: str, parent: QWidget, *, help_text: str = "",
                 key: str = "", expanded: bool = False) -> None:
        super().__init__(parent)
        self._key = key
        self._title_text = title
        self._expanded = bool(expanded)

        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        # Tight header ↔ body gap (collapsed rows stay compact).
        root.setSpacing(3)

        self._header = QWidget(self)
        self._header.setCursor(Qt.CursorShape.PointingHandCursor)
        self._header.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        head = QHBoxLayout(self._header)
        # Half the previous vertical pad so Layout / Tip names / … sit closer
        # to their dividers.
        head.setContentsMargins(0, 0, 0, 0)
        head.setSpacing(6)

        # Arrow left of title; "?" flush right when help is provided.
        self._arrow = QLabel(self._header)
        self._arrow.setFixedWidth(12)
        self._arrow.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        )
        self._title = QLabel(title, self._header)
        # Regular weight (page titles like Display stay bold elsewhere).
        self._title.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Preferred
        )
        arrow_font = self._arrow.font()
        arrow_font.setPointSizeF(max(6.0, self._title.font().pointSizeF() * 0.75))
        self._arrow.setFont(arrow_font)
        head.addWidget(self._arrow, 0)
        head.addWidget(self._title, 0)
        head.addStretch(1)
        if help_text:
            help_btn = make_help_button(self._header, title, help_text)
            help_btn.setCursor(Qt.CursorShape.ArrowCursor)
            head.addWidget(help_btn, 0, Qt.AlignmentFlag.AlignVCenter)
        self._header.mousePressEvent = self._on_header_pressed  # type: ignore[method-assign]
        root.addWidget(self._header)

        # Compatibility for older tests that looked for a checkable toggle.
        self.toggle = self._title

        self.labels: list[QLabel] = []
        self.body = QWidget(self)
        self.form = QFormLayout(self.body)
        self.form.setContentsMargins(0, 0, 0, 2)
        self.form.setHorizontalSpacing(6)
        self.form.setVerticalSpacing(2)
        self.form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
        self.form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        root.addWidget(self.body)
        self.set_expanded(self._expanded)

    def is_expanded(self) -> bool:
        return bool(self._expanded)

    def set_expanded(self, expanded: bool) -> None:
        self._expanded = bool(expanded)
        self.body.setVisible(self._expanded)
        self._arrow.setText("▼" if self._expanded else "▶")

    def _on_header_pressed(self, event) -> None:  # noqa: ANN001
        if event.button() == Qt.MouseButton.LeftButton:
            self.set_expanded(not self._expanded)
        event.accept()

    def add(self, label: str, widget: QWidget) -> None:
        text = QLabel(label, self.body)
        text.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.labels.append(text)
        self.form.addRow(text, widget)


_COL1_PX = 64      # options panel: width of the first column (number boxes, first tick box)
_COL_GAP_PX = 6    # space before the second column


def _grid_row(parent: QWidget, first: QWidget, second: QWidget) -> QWidget:
    """Two-column row of the options panel: ``first`` as wide as the number boxes,
    ``second`` starting at the same x in every row and, if it is a box, running to
    the right edge like the drop-down boxes above it."""
    row = QWidget(parent)
    lay = QHBoxLayout(row)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(_COL_GAP_PX)
    first.setFixedWidth(_COL1_PX)
    lay.addWidget(first, 0)
    if isinstance(second, QComboBox):
        second.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        lay.addWidget(second, 1)
    else:
        lay.addWidget(second, 0)
        lay.addStretch(1)
    return row


def _pair(parent: QWidget, *widgets: QWidget) -> QWidget:
    row = QWidget(parent)
    lay = QHBoxLayout(row)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(6)
    for w in widgets:
        lay.addWidget(w)
    lay.addStretch(1)
    return row


def _slider_value_row(
    parent: QWidget, slider: QWidget, value: QLabel, *, value_width: int = 44
) -> QWidget:
    """Slider fills the row; value label stays fixed on the right (aligned lengths)."""
    row = QWidget(parent)
    lay = QHBoxLayout(row)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(6)
    slider.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    value.setFixedWidth(value_width)
    value.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
    value.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    lay.addWidget(slider, 1)
    lay.addWidget(value, 0)
    return row


def _point_segment_distance(p: QPointF, a: QPointF, b: QPointF) -> float:
    dx, dy = b.x() - a.x(), b.y() - a.y()
    length2 = dx * dx + dy * dy
    t = 0.0 if length2 == 0 else max(0.0, min(1.0, ((p.x() - a.x()) * dx + (p.y() - a.y()) * dy) / length2))
    return math.hypot(p.x() - (a.x() + t * dx), p.y() - (a.y() + t * dy))


class TreeCanvas(QGraphicsView):
    """Shift+scroll to zoom, scroll to pan, left-drag to select tips (Cmd/Ctrl adds to the
    selection), middle/Alt-drag to pan.

    Clicking a branch, tip or collapsed clade reports its display node index."""

    nodeClicked = Signal(int)                 # -1 = empty space
    nodeContextMenu = Signal(int, object)     # index (-1 = empty), global QPoint
    nodeDoubleClicked = Signal(int)           # -1 = empty space
    zoomChanged = Signal()                    # after wheel / zoom_by / fit
    rubberBandNodes = Signal(object)          # list[int] tip/collapsed indices in the band
    nodeToggled = Signal(int)                 # Cmd/Ctrl+click: add or remove from the selection
    rubberBandAdd = Signal(object)            # Cmd/Ctrl+drag: add the band's tips

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setScene(QGraphicsScene(self))
        self.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        self.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
        # Left-drag = rubber-band tip select; pan with middle button or Alt+left-drag.
        self.setDragMode(QGraphicsView.DragMode.NoDrag)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setResizeAnchor(QGraphicsView.ViewportAnchor.AnchorViewCenter)
        self._press: QPoint | None = None
        self._rubber_origin: QPoint | None = None
        self._rubber: QRubberBand | None = None
        self._panning = False
        self._pan_last: QPoint | None = None

    def node_at(self, view_pos) -> int:
        """The branch, tip or collapsed clade nearest to ``view_pos`` on screen, within
        ``_PICK_RADIUS_PX`` whatever the zoom; a name or label under the mouse wins.
        -1 when nothing is that close (empty space)."""
        r = _PICK_RADIUS_PX
        pos = QPointF(view_pos)
        best, best_d = -1, float(r) + 0.001
        area = QRect(int(pos.x()) - r, int(pos.y()) - r, 2 * r + 1, 2 * r + 1)
        for item in self.items(area):
            value = item.data(ROLE_NODE)
            if value is None:
                continue
            d = self._screen_distance(item, pos)
            if d < best_d:
                best, best_d = int(value), d
                if d == 0.0 and not isinstance(item, _BranchItem):
                    break          # right on a tip name, label or triangle
        return best

    def _screen_distance(self, item, pos: QPointF) -> float:
        """Distance in screen pixels from ``pos`` to ``item`` (0 inside a name or shape)."""
        to_view = item.sceneTransform() * self.viewportTransform()
        if isinstance(item, _BranchItem):
            best = float("inf")
            for poly in item.path().toSubpathPolygons():
                pts = [to_view.map(p) for p in poly]
                for a, b in zip(pts, pts[1:]):
                    best = min(best, _point_segment_distance(pos, a, b))
            return best
        box = to_view.map(QPolygonF(item.boundingRect()))
        if box.containsPoint(pos, Qt.FillRule.OddEvenFill):
            return 0.0
        n = len(box)
        return min(_point_segment_distance(pos, box[k], box[(k + 1) % n]) for k in range(n))

    def show_node(self, index: int) -> None:
        """Scroll so the tip, label or collapsed clade ``index`` is in view."""
        best = None
        for item in self.scene().items():
            if item.data(ROLE_NODE) != index or isinstance(item, _BranchItem):
                continue
            if isinstance(item, QGraphicsSimpleTextItem):
                best = item
                break
            best = best or item
        if best is not None:
            self.ensureVisible(best.sceneBoundingRect(), 80, 80)

    def tip_nodes_in_view_rect(self, view_rect: QRect) -> list[int]:
        """Display indices of tip labels / trait dots intersecting ``view_rect``."""
        if view_rect.isEmpty():
            return []
        scene_rect = self.mapToScene(view_rect).boundingRect()
        found: list[int] = []
        seen: set[int] = set()
        for item in self.scene().items(
            scene_rect, Qt.ItemSelectionMode.IntersectsItemShape
        ):
            if not isinstance(item, (QGraphicsSimpleTextItem, QGraphicsEllipseItem)):
                continue
            value = item.data(ROLE_NODE)
            if value is None:
                continue
            idx = int(value)
            if idx in seen:
                continue
            seen.add(idx)
            found.append(idx)
        return found

    def current_scale(self) -> float:
        return abs(float(self.transform().m11()))

    def zoom_by(self, factor: float) -> None:
        current = self.transform().m11()
        target = min(40.0, max(0.02, current * factor))
        self.scale(target / current, target / current)
        self.zoomChanged.emit()

    def set_scale(self, scale: float) -> None:
        target = min(40.0, max(0.02, float(scale)))
        current = self.current_scale()
        if current <= 1e-12:
            return
        self.scale(target / current, target / current)
        self.zoomChanged.emit()

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
        self.zoomChanged.emit()

    def wheelEvent(self, event) -> None:  # noqa: N802 (Qt API)
        # Match MSA: Shift+scroll zooms; plain scroll pans the view.
        mods = event.modifiers() | QGuiApplication.keyboardModifiers()
        if not (mods & Qt.KeyboardModifier.ShiftModifier):
            super().wheelEvent(event)
            return
        pixel = event.pixelDelta()
        angle = event.angleDelta()
        if not pixel.isNull() and (pixel.y() != 0 or pixel.x() != 0):
            dy = float(pixel.y() if pixel.y() != 0 else pixel.x())
            factor = _TRACKPAD_ZOOM_BASE ** dy
        elif angle.y() != 0 or angle.x() != 0:
            delta = float(angle.y() if angle.y() != 0 else angle.x())
            factor = _ZOOM_STEP ** (delta / 120.0)
        else:
            return
        self.zoom_by(factor)
        event.accept()

    def _want_pan(self, event: QMouseEvent) -> bool:
        if event.button() == Qt.MouseButton.MiddleButton:
            return True
        if event.button() != Qt.MouseButton.LeftButton:
            return False
        mods = event.modifiers() | QGuiApplication.keyboardModifiers()
        return bool(mods & Qt.KeyboardModifier.AltModifier)

    def _clear_rubber(self) -> None:
        if self._rubber is not None:
            self._rubber.hide()
            self._rubber.deleteLater()
            self._rubber = None
        self._rubber_origin = None

    def mousePressEvent(self, event) -> None:  # noqa: N802 (Qt API)
        if self._want_pan(event):
            self._panning = True
            self._pan_last = event.position().toPoint()
            self._clear_rubber()
            self.setCursor(Qt.CursorShape.ClosedHandCursor)
            event.accept()
            return
        if event.button() == Qt.MouseButton.LeftButton:
            self._press = event.position().toPoint()
            self._rubber_origin = self._press
            event.accept()
            return
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event) -> None:  # noqa: N802 (Qt API)
        pos = event.position().toPoint()
        if self._panning and self._pan_last is not None:
            delta = pos - self._pan_last
            self._pan_last = pos
            hbar = self.horizontalScrollBar()
            vbar = self.verticalScrollBar()
            hbar.setValue(hbar.value() - delta.x())
            vbar.setValue(vbar.value() - delta.y())
            event.accept()
            return
        if (
            self._rubber_origin is not None
            and event.buttons() & Qt.MouseButton.LeftButton
        ):
            if (pos - self._rubber_origin).manhattanLength() >= 4:
                if self._rubber is None:
                    self._rubber = QRubberBand(
                        QRubberBand.Shape.Rectangle, self.viewport()
                    )
                self._rubber.setGeometry(QRect(self._rubber_origin, pos).normalized())
                self._rubber.show()
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802 (Qt API)
        pos = event.position().toPoint()
        if event.button() in (Qt.MouseButton.MiddleButton, Qt.MouseButton.LeftButton) and self._panning:
            self._panning = False
            self._pan_last = None
            self.unsetCursor()
            event.accept()
            return
        if event.button() == Qt.MouseButton.LeftButton:
            additive = bool(event.modifiers() & Qt.KeyboardModifier.ControlModifier)
            if self._rubber is not None and self._rubber.isVisible():
                rect = self._rubber.geometry()
                self._clear_rubber()
                nodes = self.tip_nodes_in_view_rect(rect)
                if additive:
                    self.rubberBandAdd.emit(nodes)
                else:
                    self.rubberBandNodes.emit(nodes)
            elif (
                self._press is not None
                and (pos - self._press).manhattanLength() < 4
            ):
                if additive:
                    self.nodeToggled.emit(self.node_at(pos))
                else:
                    self.nodeClicked.emit(self.node_at(pos))
            self._press = None
            self._clear_rubber()
            event.accept()
            return
        super().mouseReleaseEvent(event)

    def contextMenuEvent(self, event) -> None:  # noqa: N802 (Qt API)
        self.nodeContextMenu.emit(self.node_at(event.pos()), event.globalPos())
        event.accept()

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802 (Qt API)
        # The panel decides: a tip opens its sequence, elsewhere the tree is fitted.
        self.nodeDoubleClicked.emit(self.node_at(event.position().toPoint()))
        event.accept()


class TreePanel(QWidget):
    """Controls + canvas for one tree file."""

    fileWritten = Signal(object)  # Path of an exported figure or saved tree
    # (alignment path, record ids, folder for extracts): open the alignment in
    # Preview with these rows selected; its extracts are saved next to the tree
    alignmentRequested = Signal(object, object, object)

    def __init__(self, path: Path | None = None, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("TreePanel")
        self._path: Path | None = None
        self._file_model: tl.TreeModel | None = None
        self._model: tl.TreeModel | None = None
        self._state = TreeViewState()
        self._undo: list[tuple[tl.TreeModel, TreeViewState, TreeViewOptions]] = []
        self._redo: list[tuple[tl.TreeModel, TreeViewState, TreeViewOptions]] = []
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
        self._align_check = QCheckBox("Align names", self)
        self._crowded_check = QCheckBox("Hide crowded", self)
        self._zeros_check = QCheckBox("Hide zeros", self)
        self._zeros_check.setToolTip("No label on branches whose length is 0 or rounds to 0.")
        self._lengths_box = NoWheelComboBox(self)
        for label, value in (("Proportional", tl.BRANCHES_PROPORTIONAL),
                             ("Equal", tl.BRANCHES_EQUAL),
                             ("Square root", tl.BRANCHES_SQRT)):
            self._lengths_box.addItem(label, value)
        self._spread_slider = NoWheelSlider(Qt.Orientation.Horizontal, self)
        self._spread_slider.setObjectName("annotationZoomSlider")
        self._spread_slider.setRange(SPREAD_MIN, SPREAD_MAX)
        self._spread_slider.setSingleStep(25)
        self._spread_slider.setPageStep(50)
        self._spread_slider.setValue(100)
        self._spread_slider.setMinimumWidth(72)
        self._spread_slider.setFixedHeight(18)
        self._spread_slider.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        apply_annotation_zoom_slider_style(self._spread_slider)
        self._spread_slider.setToolTip(
            "Tree size relative to the font (circular: also opens the centre)."
        )
        self._spread_label = QLabel("100%", self)
        self._spread_label.setMinimumWidth(40)
        self._spread_label.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        self._font_spin = NoWheelSpinBox(self)
        self._font_spin.setRange(5, 36)
        self._font_spin.setSuffix(" pt")
        self._font_spin.setFixedWidth(70)
        self._fit_btn = QPushButton("Fit", self)
        self._fit_btn.setToolTip("Show the whole tree (or double-click an empty part of the tree).")
        self._fit_btn.setFixedWidth(ANNOTATION_FIT_BUTTON_RENDERED_WIDTH)
        self._fit_btn.setFixedHeight(OVERVIEW_SELECT_FIT_HEIGHT)
        self._fit_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._fit_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._fit_btn.clicked.connect(self._on_fit_clicked)
        # Export opens an MSA-style Type + File dialog (SVG / PDF / PNG).
        self._save_btn = QPushButton("Export", self)
        self._save_btn.setFixedHeight(OVERVIEW_SELECT_FIT_HEIGHT)
        self._save_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._save_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._save_btn.setToolTip(
            "Export the tree figure (SVG, PDF or PNG) — choose type and file."
        )
        self._save_btn.clicked.connect(self._open_export_dialog)

        # Edit lives on the Display rail (native chrome like MSA Delete/Save).
        # Undo is Ctrl/Cmd+Z (and in the Edit menu) — no separate Undo button.
        self._edit_btn = QPushButton("Edit", self)
        self._edit_btn.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._edit_btn.setMinimumWidth(0)
        self._edit_btn.setToolTip("Edit the selected tip or clade. Undo: Ctrl+Z / Cmd+Z.")
        self._edit_menu = QMenu(self._edit_btn)
        self._edit_menu.aboutToShow.connect(self._fill_edit_menu)
        self._edit_btn.setMenu(self._edit_menu)

        zoom_out = make_overview_zoom_button("−", self)
        zoom_out.setToolTip("Zoom out (Shift+scroll)")
        zoom_out.clicked.connect(lambda: self._nudge_zoom(-3))
        zoom_in = make_overview_zoom_button("+", self)
        zoom_in.setToolTip("Zoom in (Shift+scroll)")
        zoom_in.clicked.connect(lambda: self._nudge_zoom(3))
        self._zoom_out_btn = zoom_out
        self._zoom_in_btn = zoom_in
        self._zoom_slider = NoWheelSlider(Qt.Orientation.Horizontal, self)
        self._zoom_slider.setObjectName("annotationZoomSlider")
        self._zoom_slider.setRange(0, 100)
        self._zoom_slider.setValue(0)
        self._zoom_slider.setFixedWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
        self._zoom_slider.setMinimumWidth(ANNOTATION_ZOOM_SLIDER_WIDTH)
        self._zoom_slider.setFixedHeight(ZOOM_SLIDER_HEIGHT)
        self._zoom_slider.setToolTip(
            "Tree size: left = fit, right = zoom in. Drag, − / +, or Shift+scroll."
        )
        apply_annotation_zoom_slider_style(self._zoom_slider)
        self._zoom_slider.valueChanged.connect(self._on_zoom_slider)
        self._fit_scale = 1.0
        self._zoom_block = False
        self._collapse_spin = NoWheelSpinBox(self)
        self._collapse_spin.setRange(0, 100)
        self._collapse_spin.setSpecialValueText("Off")
        self._collapse_spin.setToolTip("Merge branches with support below this value (Off = 0).")
        self._collapse_spin.editingFinished.connect(self._on_collapse_threshold)

        self._color_box = NoWheelComboBox(self)
        self._color_box.currentIndexChanged.connect(self._on_color_by_changed)
        self._traits_btn = QPushButton("Trait table", self)
        traits_menu = QMenu(self._traits_btn)
        traits_menu.addAction("Create or update…", self._open_trait_dialog)
        traits_menu.addAction("Reload", self._reload_traits)
        self._open_traits_action = traits_menu.addAction("Open in spreadsheet app",
                                                         self._open_traits_file)
        self._traits_btn.setMenu(traits_menu)

        # --- controls added with the options panel -----------------------------
        # Same slider as the Annotation and Dotplot panels (accent track, round
        # handle, no wheel; a click on the track jumps there).
        self._angle_slider = NoWheelSlider(Qt.Orientation.Horizontal, self)
        self._angle_slider.setObjectName("annotationZoomSlider")
        self._angle_slider.setRange(0, 71)            # steps of 5°
        self._angle_slider.setPageStep(9)
        self._angle_slider.setMinimumWidth(72)
        self._angle_slider.setFixedHeight(18)
        self._angle_slider.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        apply_annotation_zoom_slider_style(self._angle_slider)
        self._angle_label = QLabel("0°", self)
        self._angle_label.setMinimumWidth(36)
        self._angle_label.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        self._flip_check = QCheckBox("Flip", self)
        self._flip_check.setToolTip("Mirror the tree left to right.")
        self._tiplabel_box = NoWheelComboBox(self)
        self._maxchars_spin = NoWheelSpinBox(self)
        self._maxchars_spin.setRange(0, 200)
        self._maxchars_spin.setSpecialValueText("Off")
        self._maxchars_spin.setToolTip("Shorten tip names to this many characters (Off = full names).")
        self._keep_box = NoWheelComboBox(self)
        self._keep_box.addItem("end", True)
        self._keep_box.addItem("start", False)
        self._keep_box.setToolTip("Which part of a shortened name to keep.")
        self._digits_spin = NoWheelSpinBox(self)
        self._digits_spin.setRange(1, 8)
        self._bfont_spin = NoWheelSpinBox(self)
        self._bfont_spin.setRange(0, 36)
        self._bfont_spin.setSpecialValueText("Auto")
        self._bfont_spin.setSuffix(" pt")
        self._lw_spin = NoWheelDoubleSpinBox(self)
        self._lw_spin.setRange(0.5, 6.0)
        self._lw_spin.setSingleStep(0.25)
        self._lw_spin.setDecimals(2)
        self._collapse_dist = NoWheelDoubleSpinBox(self)
        self._collapse_dist.setRange(0.0, 100.0)
        self._collapse_dist.setDecimals(4)
        self._collapse_dist.setSingleStep(0.005)
        self._collapse_dist.setSpecialValueText("Off")
        self._collapse_dist.editingFinished.connect(self._on_collapse_distance)
        self._collapse_dist.setToolTip("Collapse clades whose tips lie closer than this (Off = 0).")
        self._scale_len = NoWheelDoubleSpinBox(self)
        self._scale_len.setRange(0.0, 1000.0)
        self._scale_len.setDecimals(4)
        self._scale_len.setSingleStep(0.01)
        self._scale_len.setSpecialValueText("Auto")
        self._sfont_spin = NoWheelSpinBox(self)
        self._sfont_spin.setRange(0, 36)
        self._sfont_spin.setSpecialValueText("Auto")
        self._sfont_spin.setSuffix(" pt")
        self._stats_label = QLabel(self)
        self._stats_label.setWordWrap(True)
        self._stats_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._tips_check.setText("Show")
        self._align_check.setText("Align")
        self._scale_check.setText("Show")
        # Panel fields: combos fill the width, number boxes share one small width.
        for combo in (self._layout_box, self._order_box, self._lengths_box, self._branch_box,
                      self._tiplabel_box, self._color_box, self._keep_box):
            combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
            combo.setMinimumContentsLength(4)
            combo.setMaximumWidth(16777215)
            combo.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        for spin in (self._maxchars_spin, self._font_spin, self._digits_spin,
                     self._bfont_spin, self._lw_spin, self._collapse_spin, self._collapse_dist,
                     self._scale_len, self._sfont_spin):
            spin.setFixedWidth(_COL1_PX)
        self._keep_box.setMinimumContentsLength(3)

        # --- top row: − / + / slider · Fit · Export (MSA order)
        self._status = QLabel(self)
        self._status.setObjectName("TreeToolbarStatus")
        self._status.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        self._status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._status.setStyleSheet("color: palette(mid);")

        top = QHBoxLayout()
        top.setSpacing(6)
        top.setAlignment(Qt.AlignmentFlag.AlignVCenter)
        top.addWidget(self._zoom_out_btn, 0, Qt.AlignmentFlag.AlignVCenter)
        top.addWidget(self._zoom_in_btn, 0, Qt.AlignmentFlag.AlignVCenter)
        top.addWidget(self._zoom_slider, 0, Qt.AlignmentFlag.AlignVCenter)
        top.addWidget(self._fit_btn, 0, Qt.AlignmentFlag.AlignVCenter)
        top.addWidget(self._save_btn, 0, Qt.AlignmentFlag.AlignVCenter)
        top.addStretch(1)
        top.addWidget(self._status, 0, Qt.AlignmentFlag.AlignVCenter)

        # --- options panel (right) — MSA Display section + divider pattern -----
        panel = QWidget(self)
        pl = QVBoxLayout(panel)
        pl.setContentsMargins(10, 6, 10, 6)
        # Half the old gap between section headers / dividers.
        pl.setSpacing(4)

        def _add_section(sec: _Section) -> None:
            if pl.count() > 0:
                pl.addWidget(_h_divider(panel))
            pl.addWidget(sec)

        sec = _Section("Layout", panel, help_text=DRAW_HELP, key="layout")
        sec.add("Layout", self._layout_box)
        sec.add("Order", self._order_box)
        sec.add("Branches", self._lengths_box)
        # Same value-column width so Spread / Angle tracks end flush.
        self._spread_row = _slider_value_row(
            panel, self._spread_slider, self._spread_label, value_width=44
        )
        sec.add("Spread", self._spread_row)
        self._angle_row = _slider_value_row(
            panel, self._angle_slider, self._angle_label, value_width=44
        )
        sec.add("Angle", self._angle_row)
        sec.add("", self._flip_check)
        self._layout_section = sec
        _add_section(sec)
        sec = _Section("Tip names", panel, help_text=TIPS_HELP, key="tips")
        sec.add("", _grid_row(panel, self._tips_check, self._align_check))
        sec.add("Text", self._tiplabel_box)
        sec.add("Max", _grid_row(panel, self._maxchars_spin, self._keep_box))
        sec.add("Font", self._font_spin)
        _add_section(sec)
        sec = _Section("Branches", panel, help_text=SUPPORT_HELP, key="branches")
        sec.add("Labels", self._branch_box)
        sec.add("Digits", self._digits_spin)
        sec.add("Font", self._bfont_spin)
        sec.add("", self._crowded_check)
        sec.add("", self._zeros_check)
        sec.add("Width", self._lw_spin)
        _add_section(sec)
        sec = _Section("Collapse", panel, help_text=COLLAPSE_HELP, key="collapse")
        sec.add("Support", self._collapse_spin)
        sec.add("Distance", self._collapse_dist)
        _add_section(sec)
        sec = _Section("Colour tips", panel, help_text=TRAITS_HELP, key="colour")
        sec.add("By", self._color_box)
        sec.add("", _pair(panel, self._traits_btn))
        _add_section(sec)
        sec = _Section("Scale bar", panel, help_text=SCALE_HELP, key="scale")
        sec.add("", self._scale_check)
        sec.add("Length", self._scale_len)
        sec.add("Font", self._sfont_spin)
        _add_section(sec)
        self._report_btn = QPushButton("Open run report", self)
        self._report_btn.setToolTip("The tree run's report (in its _details folder).")
        self._report_btn.clicked.connect(self._open_report)
        sec = _Section("Statistics", panel, key="stats")
        sec.form.addRow(self._stats_label)
        sec.form.addRow(_pair(panel, self._report_btn))
        _add_section(sec)
        pl.addStretch(1)
        # one label column for all sections, so every field starts at the same x
        labels = [lab for sec in panel.findChildren(_Section) for lab in sec.labels]
        width = max((lab.sizeHint().width() for lab in labels), default=0)
        for lab in labels:
            lab.setFixedWidth(width)
        # Keep section titles at MSA Tracks size (do not shrink with MacSmallSize).
        for sec in panel.findChildren(_Section):
            sec._title.setAttribute(Qt.WidgetAttribute.WA_MacSmallSize, False)

        # Edit footer (MSA Annotation Delete/Save button chrome; Undo = Ctrl/Cmd+Z).
        edit_footer = QWidget(self)
        edit_footer.setObjectName("treeDisplayFooter")
        efl = QVBoxLayout(edit_footer)
        efl.setContentsMargins(10, 8, 10, 10)
        efl.setSpacing(8)
        efl.addWidget(_h_divider(edit_footer))
        edit_row = QHBoxLayout()
        edit_row.setContentsMargins(0, 0, 0, 0)
        edit_row.setSpacing(8)
        edit_help = make_help_button(edit_footer, "Editing the tree", EDIT_HELP)
        edit_row.addWidget(self._edit_btn, 1)
        edit_row.addWidget(edit_help, 0)
        efl.addLayout(edit_row)

        # MSA-style right rail: Window surface + Display tab + section dividers.
        self._find_hits: list[str] = []
        self._find_index = -1
        self._find_signature: tuple | None = None
        self._rail = TreeDisplayRail(panel, self, footer=edit_footer,
                                     find_page=self._build_find_page())
        self._options_panel = self._rail.options_scroll()

        self._canvas = TreeCanvas(self)
        # No box and no own background: the tree sits on whatever the Preview
        # panel paints (white in RepeatViewer's light theme, dark in dark mode).
        self._canvas.setFrameShape(QGraphicsView.Shape.NoFrame)
        self._canvas.setBackgroundBrush(Qt.BrushStyle.NoBrush)
        self._canvas.viewport().setAutoFillBackground(False)
        self._canvas.setStyleSheet("QGraphicsView { background: transparent; border: none; }")
        self._canvas.nodeClicked.connect(self._on_node_clicked)
        self._canvas.nodeContextMenu.connect(self._on_context_menu)
        self._canvas.nodeDoubleClicked.connect(self._on_node_double_clicked)
        self._canvas.rubberBandNodes.connect(self._on_rubber_band_nodes)
        self._canvas.nodeToggled.connect(self._on_node_toggled)
        self._canvas.rubberBandAdd.connect(self._on_rubber_band_add)
        self._canvas.zoomChanged.connect(self._sync_zoom_slider_from_canvas)

        # Toolbar only over the tree column — Display rail is full height (MSA-like).
        # Status (tips / nodes / support legend) sits on the right, like MSA sequence count.
        toolbar = QWidget(self)
        toolbar.setObjectName("TreeToolbar")
        toolbar.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        toolbar.setFixedHeight(OVERVIEW_TOOLBAR_HEIGHT)
        toolbar.setStyleSheet("QWidget#TreeToolbar { background: palette(midlight);"
                              " border: none; border-bottom: 1px solid palette(mid); }")
        top.setContentsMargins(6, 2, 8, 2)
        toolbar.setLayout(top)
        left = QVBoxLayout()
        left.setContentsMargins(0, 0, 0, 0)
        left.setSpacing(0)
        left.addWidget(toolbar, 0)
        left.addWidget(self._canvas, stretch=1)
        middle = QHBoxLayout()
        middle.setContentsMargins(0, 0, 0, 0)
        middle.setSpacing(0)
        middle.addLayout(left, stretch=1)
        middle.addWidget(self._rail, stretch=0)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addLayout(middle, stretch=1)

        self._apply_options_to_controls()
        for box in (self._layout_box, self._order_box, self._branch_box, self._lengths_box,
                    self._tiplabel_box, self._keep_box):
            box.currentIndexChanged.connect(self._on_controls_changed)
        for check in (self._tips_check, self._scale_check, self._align_check,
                      self._crowded_check, self._zeros_check, self._flip_check):
            check.toggled.connect(self._on_controls_changed)
        for spin in (self._font_spin, self._maxchars_spin, self._digits_spin,
                     self._bfont_spin, self._lw_spin, self._scale_len, self._sfont_spin):
            spin.valueChanged.connect(self._on_controls_changed)
        self._spread_slider.valueChanged.connect(self._on_controls_changed)
        self._angle_slider.valueChanged.connect(self._on_controls_changed)

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
                               (self._branch_box, self._opts.branch_labels),
                               (self._lengths_box, self._opts.branch_lengths)):
                idx = box.findData(value)
                box.setCurrentIndex(max(0, idx))
            self._tips_check.setChecked(self._opts.tip_labels)
            self._scale_check.setChecked(self._opts.scale_bar)
            self._align_check.setChecked(self._opts.align_labels)
            self._crowded_check.setChecked(self._opts.hide_crowded)
            self._zeros_check.setChecked(self._opts.hide_zeros)
            self._font_spin.setValue(self._opts.font_pt)
            spread_pct = int(spread_factor(self._opts) * 100)
            self._spread_slider.setValue(spread_pct)
            self._spread_label.setText(f"{spread_pct}%")
            angle = int(self._opts.root_angle) % 360
            self._angle_slider.setValue(angle // 5)
            self._angle_label.setText(f"{angle}°")
            self._flip_check.setChecked(self._opts.flip)
            self._maxchars_spin.setValue(max(0, int(self._opts.max_chars)))
            self._keep_box.setCurrentIndex(0 if self._opts.keep_end else 1)
            self._digits_spin.setValue(max(1, min(8, int(self._opts.branch_decimals))))
            self._bfont_spin.setValue(max(0, int(self._opts.branch_font_pt)))
            self._lw_spin.setValue(float(self._opts.line_weight))
            self._scale_len.setValue(max(0.0, float(self._opts.scale_length)))
            self._sfont_spin.setValue(max(0, int(self._opts.scale_font_pt)))
            self._collapse_dist.setValue(self._state.collapse_distance)
            self._fill_tiplabel_box()
            self._sync_enabled_controls()
            self._collapse_spin.setValue(self._state.collapse_below)
        finally:
            self._block = False

    def _sync_enabled_controls(self) -> None:
        kind = self._layout_box.currentData()
        radial = kind == tl.LAYOUT_RADIAL
        tips = self._tips_check.isChecked()
        self._align_check.setEnabled(not radial and tips)
        self._align_check.setToolTip("Not for the radial layout." if radial else
                                     "All names on one line or circle, dotted lines to the tips.")
        for w in (self._tiplabel_box, self._maxchars_spin, self._keep_box):
            w.setEnabled(tips)
        self._tiplabel_box.setEnabled(tips and self._tiplabel_box.count() > 1)
        self._keep_box.setEnabled(tips and self._maxchars_spin.value() > 0)
        labels = self._branch_box.currentData()
        self._crowded_check.setEnabled(labels != BRANCH_LABELS_NONE)
        self._zeros_check.setEnabled(labels == BRANCH_LABELS_LENGTH)
        self._bfont_spin.setEnabled(labels != BRANCH_LABELS_NONE)
        self._digits_spin.setEnabled(labels == BRANCH_LABELS_LENGTH)
        proportional = self._lengths_box.currentData() == tl.BRANCHES_PROPORTIONAL
        self._scale_check.setEnabled(proportional)
        for w in (self._scale_len, self._sfont_spin):
            w.setEnabled(proportional and self._scale_check.isChecked())
        # Angle turns circular / radial trees; keep the row visible for all layouts.
        rect = kind == tl.LAYOUT_RECTANGULAR
        self._layout_section.form.setRowVisible(self._angle_row, True)
        self._angle_slider.setEnabled(not rect)
        self._angle_label.setEnabled(not rect)
        self._angle_slider.setToolTip(
            "Rectangular trees have no angle — switch Layout to Circular or Radial."
            if rect
            else "Turn the tree (steps of 5°)."
        )

    def _fill_tiplabel_box(self) -> None:
        """Name, and name + / only each trait-table column (when there is a table)."""
        box = self._tiplabel_box
        wanted = self._opts.tip_label
        was = self._block
        self._block = True
        try:
            box.clear()
            box.addItem("Name", TIP_LABEL_NAME)
            for col in (self._traits.columns if self._traits else []):
                if col == "header":
                    continue
                box.addItem(f"Name + {col}", f"name+col:{col}")
                box.addItem(f"{col} only", f"col:{col}")
            idx = box.findData(wanted)
            box.setCurrentIndex(max(0, idx))
            if idx < 0:
                self._opts.tip_label = TIP_LABEL_NAME
        finally:
            self._block = was

    # ------------------------------------------------------------ find
    def _build_find_page(self) -> QWidget:
        page = QWidget(self)
        layout = QVBoxLayout(page)
        layout.setContentsMargins(10, 10, 10, 10)
        layout.setSpacing(8)
        title, _rule = _page_title_with_rule("Find", page)
        layout.addWidget(title)
        self._find_query = QLineEdit(page)
        self._find_query.setPlaceholderText("Tip name or trait value…")
        self._find_query.setClearButtonEnabled(True)
        self._find_query.returnPressed.connect(lambda: self.find_step(1))
        form = QFormLayout()
        form.setContentsMargins(0, 0, 0, 0)
        form.setHorizontalSpacing(8)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        find_row = QHBoxLayout()
        find_row.setContentsMargins(0, 0, 0, 0)
        find_row.setSpacing(6)
        find_row.addWidget(self._find_query, 1)
        find_row.addWidget(make_help_button(page, "Find in the tree", FIND_HELP), 0)
        form.addRow("Find:", find_row)
        layout.addLayout(form)
        layout.addWidget(_section_label("Search in", page))
        self._find_names_check = QCheckBox("Tip names", page)
        self._find_names_check.setChecked(True)
        self._find_traits_check = QCheckBox("Trait values", page)
        self._find_traits_check.setChecked(True)
        self._find_case_check = QCheckBox("Case sensitive", page)
        for box in (self._find_names_check, self._find_traits_check, self._find_case_check):
            layout.addWidget(box)
        self._find_prev_btn = QPushButton("Find previous", page)
        self._find_next_btn = QPushButton("Find next", page)
        self._find_all_btn = QPushButton("Find all", page)
        self._find_prev_btn.clicked.connect(lambda: self.find_step(-1))
        self._find_next_btn.clicked.connect(lambda: self.find_step(1))
        self._find_all_btn.clicked.connect(self.find_all)
        nav = QHBoxLayout()
        nav.setContentsMargins(0, 0, 0, 0)
        nav.setSpacing(6)
        nav.addWidget(self._find_prev_btn, 1)
        nav.addWidget(self._find_next_btn, 1)
        layout.addLayout(nav)
        layout.addWidget(self._find_all_btn)
        layout.addStretch(1)
        return page

    def focus_find(self) -> None:
        """Edit > Find: open the Find tab and put the cursor in the search box."""
        self._rail.show_tab(TreeDisplayRail.TAB_FIND)
        self._find_query.setFocus()
        self._find_query.selectAll()

    def _find_order(self) -> list[tuple[str, int]]:
        """(tip name, display node that shows it) from top to bottom of the drawing;
        tips inside a collapsed clade point at the clade."""
        d = self._display
        if d is None:
            return []
        out: list[tuple[str, int]] = []
        stack = [0]
        while stack:
            i = stack.pop()
            node = d.model.nodes[i]
            if i in d.info:
                out.extend((name, i) for name in d.keys[i])
            elif node.is_tip:
                out.append((node.name, i))
            else:
                stack.extend(reversed(node.children))
        return out

    def _compute_find(self) -> tuple[list[str], dict[str, int]]:
        order = self._find_order()
        table = self._traits.rows if self._traits is not None else None
        hits = tf.find_tips(
            [name for name, _i in order], self._find_query.text(), traits=table,
            in_names=self._find_names_check.isChecked(),
            in_traits=self._find_traits_check.isChecked(),
            case_sensitive=self._find_case_check.isChecked())
        return hits, dict(order)

    def _find_state(self) -> tuple:
        return (self._find_query.text(), self._find_names_check.isChecked(),
                self._find_traits_check.isChecked(), self._find_case_check.isChecked(),
                id(self._model))

    def find_step(self, direction: int) -> None:
        """Select the next (+1) or previous (-1) match and scroll to it."""
        if self._display is None:
            return
        signature = self._find_state()
        if signature != self._find_signature:
            self._find_signature = signature
            self._find_index = -1
        hits, where = self._compute_find()
        self._find_hits = hits
        if not hits:
            self._update_status("Find: no match" if self._find_query.text().strip() else "")
            return
        self._find_index = tf.step_index(self._find_index, len(hits), direction)
        tip = hits[self._find_index]
        node = where[tip]
        d = self._display
        self._selected = d.keys[node] if node in d.info else (tip,)
        self._redraw(refit=False)
        self._canvas.show_node(node)
        note = f"Find: {self._find_index + 1} of {len(hits)}"
        if node in d.info:
            note += " (inside a collapsed clade)"
        self._update_status(note)

    def find_all(self) -> None:
        """Select every match."""
        if self._display is None:
            return
        self._find_signature = self._find_state()
        hits, _where = self._compute_find()
        self._find_hits = hits
        self._find_index = -1
        self._selected = tuple(sorted(hits)) if hits else None
        self._redraw(refit=False)
        self._update_status(f"Find: {len(hits)} selected" if hits else "Find: no match")

    # ------------------------------------------------------------ state
    def _show_options_panel(self, on: bool, *, save: bool = True) -> None:
        """Show or collapse the Display rail (keeps the vertical Display tab)."""
        self._rail.set_panel_collapsed(not bool(on), persist=save)

    def _on_fit_clicked(self) -> None:
        self._canvas.fit(max_scale=4.0)
        self._fit_scale = max(1e-6, self._canvas.current_scale())
        self._set_zoom_slider_value(0)

    def _nudge_zoom(self, delta: int) -> None:
        self._zoom_slider.setValue(
            max(0, min(100, int(self._zoom_slider.value()) + int(delta)))
        )

    def _set_zoom_slider_value(self, value: int) -> None:
        self._zoom_block = True
        try:
            self._zoom_slider.setValue(max(0, min(100, int(value))))
        finally:
            self._zoom_block = False

    def _on_zoom_slider(self, value: int) -> None:
        if self._zoom_block:
            return
        base = self._fit_scale if self._fit_scale > 1e-9 else 1.0
        # 0 → fit scale; 100 → 16× fit (comfortable for tip labels).
        target = base * (16.0 ** (max(0, min(100, int(value))) / 100.0))
        self._zoom_block = True
        try:
            self._canvas.set_scale(target)
        finally:
            self._zoom_block = False

    def _sync_zoom_slider_from_canvas(self) -> None:
        if self._zoom_block:
            return
        base = self._fit_scale if self._fit_scale > 1e-9 else self._canvas.current_scale()
        if base <= 1e-12:
            return
        ratio = max(1e-12, self._canvas.current_scale() / base)
        # Invert scale = base * 16^(v/100)
        value = int(round(100.0 * math.log(ratio) / math.log(16.0)))
        self._set_zoom_slider_value(value)

    def set_options(self, opts: TreeViewOptions) -> None:
        self._opts = replace(opts)
        self._apply_options_to_controls()
        self._redraw(refit=True)

    def _on_controls_changed(self, *_args) -> None:
        if self._block:
            return
        global _last_options
        old = self._opts
        spread_pct = int(self._spread_slider.value())
        self._spread_label.setText(f"{spread_pct}%")
        self._angle_label.setText(f"{self._angle_slider.value() * 5}°")
        self._opts = replace(
            old,
            layout=str(self._layout_box.currentData()),
            order=str(self._order_box.currentData()),
            tip_labels=self._tips_check.isChecked(),
            branch_labels=str(self._branch_box.currentData()),
            font_pt=self._font_spin.value(),
            scale_bar=self._scale_check.isChecked(),
            align_labels=self._align_check.isChecked(),
            branch_lengths=str(self._lengths_box.currentData()),
            spread=spread_pct,
            hide_crowded=self._crowded_check.isChecked(),
            hide_zeros=self._zeros_check.isChecked(),
            root_angle=self._angle_slider.value() * 5,
            flip=self._flip_check.isChecked(),
            tip_label=str(self._tiplabel_box.currentData() or TIP_LABEL_NAME),
            max_chars=self._maxchars_spin.value(),
            keep_end=bool(self._keep_box.currentData()),
            branch_decimals=self._digits_spin.value(),
            branch_font_pt=self._bfont_spin.value(),
            line_weight=self._lw_spin.value(),
            scale_length=self._scale_len.value(),
            scale_font_pt=self._sfont_spin.value(),
        )
        self._sync_enabled_controls()
        _last_options = replace(self._opts)
        refit = (self._opts.layout != old.layout or self._opts.spread != old.spread
                 or self._opts.branch_lengths != old.branch_lengths
                 or self._opts.align_labels != old.align_labels
                 or self._opts.flip != old.flip or self._opts.max_chars != old.max_chars
                 or self._opts.root_angle != old.root_angle)
        self._redraw(refit=refit)
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
        self._redo.clear()
        self._selected = None
        self._state = TreeViewState()
        try:
            self._file_model = tl.read_tree_file(self._path)
        except (OSError, tl.TreeReadError) as exc:
            self._file_model = self._model = self._display = None
            self._canvas.scene().clear()
            self._status.setText(f"Could not read {self._path.name}: {exc}")
            return False
        self._model = self._file_model
        self._load_view_state()
        self._load_traits()
        if not self._model.has_support and self._opts.branch_labels == BRANCH_LABELS_SUPPORT:
            self._opts = replace(self._opts, branch_labels=BRANCH_LABELS_NONE)
        self._apply_options_to_controls()
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
        try:
            self._state.collapse_distance = max(0.0, float(data.get("collapse_distance", 0.0)))
        except (TypeError, ValueError):
            self._state.collapse_distance = 0.0
        self._state.expanded = {tuple(k) for k in data.get("expanded", []) if k}
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
            "collapse_distance": self._state.collapse_distance,
            "expanded": [list(k) for k in sorted(self._state.expanded)],
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
        # transparent: the Preview panel's own background shows through
        # (figures are still black on white)
        self._canvas.setBackgroundBrush(Qt.BrushStyle.NoBrush)
        self._display = make_display(self._model, self._opts, self._state)
        self._apply_traits(self._display)
        build_tree_scene(self._canvas.scene(), self._display, self._opts, ink=ink, muted=muted,
                         selected=self._selected)
        if refit:
            self._canvas.fit(max_scale=1.0)
            self._fit_scale = max(1e-6, self._canvas.current_scale())
            self._set_zoom_slider_value(0)
        self._update_status()

    def _update_status(self, message: str = "") -> None:
        m, d = self._model, self._display
        if m is None or d is None:
            return
        parts = [f"{m.n_tips} tips", f"{m.n_internal} internal nodes"]
        depth = max(tl.depths(m)) if m.has_lengths else 0.0
        stats = [f"{m.n_tips} tips", f"{m.n_internal} internal nodes"]
        if m.has_lengths:
            stats.append(f"Tree depth: {tl.branch_length_text(depth, 4)}")
        if len(d.info):
            stats.append(f"Collapsed clades shown: {len(d.info)}")
        self._stats_label.setText("\n".join(stats))
        self._report_btn.setEnabled(self.report_path() is not None)
        if not m.has_lengths:
            parts.append("no branch lengths in the file")
        support_tip = ""
        if m.has_support and self._opts.branch_labels == BRANCH_LABELS_SUPPORT:
            if any("/" in n.support for n in m.nodes):
                parts.append("support: SH-aLRT/UFBoot")
                support_tip = (
                    "Branch labels are SH-aLRT / UFBoot (e.g. 92/97).\n"
                    "Well supported at SH-aLRT ≥ 80 and UFBoot ≥ 95."
                )
            else:
                support_tip = (
                    "Branch labels are support values from the tree file.\n"
                    "FastME bootstrap: trust ≥ 70. FastTree local support: 0–1, near 1 is strong."
                )
        if d.n_contracted:
            parts.append(f"{d.n_contracted} weak branches collapsed")
        if self._state.subtree:
            parts.append(f"showing one clade ({len(self._state.subtree)} tips)")
        if self._state.edited:
            parts.append(f"tree {self._state.edited}: Save → Tree keeps it")
        if self._selected:
            parts.append(f"selected: {len(self._selected)} tips")
        if message:
            parts.append(message)
        self._status.setText(" · ".join(parts))
        self._status.setToolTip(support_tip)

    def showEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().showEvent(event)
        if self._model is not None:
            self._canvas.fit(max_scale=1.0)
            self._fit_scale = max(1e-6, self._canvas.current_scale())
            self._set_zoom_slider_value(0)

    # ------------------------------------------------------------ selection and menus

    def _key_for(self, index: int) -> tuple[str, ...] | None:
        if self._display is None or index < 0 or index >= len(self._display.keys):
            return None
        return self._display.keys[index]

    def _on_node_clicked(self, index: int) -> None:
        key = self._key_for(index)
        self._selected = None if key == self._selected else key
        self._redraw(refit=False)

    def _band_names(self, indices) -> set[str]:
        """Tip names (a collapsed clade counts all its tips) for display indices."""
        d = self._display
        names: set[str] = set()
        if d is None:
            return names
        for raw in indices or ():
            try:
                i = int(raw)
            except (TypeError, ValueError):
                continue
            if i < 0 or i >= len(d.keys):
                continue
            node = d.model.nodes[i]
            if node.is_tip or i in d.info:
                names.update(d.keys[i])
        return names

    def _on_rubber_band_nodes(self, indices) -> None:
        """Left-drag rectangle: select every tip (or collapsed clade) whose label is hit."""
        if self._display is None:
            return
        names = self._band_names(indices)
        self._selected = tuple(sorted(names)) if names else None
        self._redraw(refit=False)

    def _on_rubber_band_add(self, indices) -> None:
        """Cmd/Ctrl + drag: add the tips in the rectangle to the selection."""
        names = self._band_names(indices)
        if not names:
            return
        self._selected = tuple(sorted(set(self._selected or ()) | names))
        self._redraw(refit=False)

    def _on_node_toggled(self, index: int) -> None:
        """Cmd/Ctrl + click: add the tip, branch or collapsed clade to the selection,
        or take it out when all its tips are already selected."""
        key = self._key_for(index)
        if not key:
            return
        current = set(self._selected or ())
        if all(name in current for name in key):
            current -= set(key)
        else:
            current |= set(key)
        self._selected = tuple(sorted(current)) if current else None
        self._redraw(refit=False)

    def _tip_at(self, index: int) -> str | None:
        """Name of the tip at a display index (not a collapsed clade), else None."""
        d = self._display
        if d is None or index < 0 or index >= len(d.keys) or index in d.info:
            return None
        node = d.model.nodes[index]
        return node.name if node.is_tip else None

    def _on_node_double_clicked(self, index: int) -> None:
        tip = self._tip_at(index)
        if tip is None:
            self._canvas.fit(max_scale=4.0)
            return
        self._selected = self._key_for(index)
        self._redraw(refit=False)
        self.show_in_alignment(tip)

    def _sequence_source(self):
        from repeatviewer.backends import tree_links

        return tree_links.sequences_for_tree(self._path) if self._path else None

    def show_in_alignment(self, tips) -> bool:
        """Open the tree's alignment in a Preview tab with ``tips`` (one name or a
        clade) selected, so its right-click menu (Extract, Copy…) works on them."""
        from repeatviewer.backends import tree_links

        tips = [tips] if isinstance(tips, str) else list(tips)
        source = self._sequence_source()
        if source is None:
            QMessageBox.information(
                self, "Show in alignment",
                "No alignment or FASTA was found for this tree. RepeatViewer looks for "
                "the tree run's .aln.fa, or a FASTA with the tree's name next to it.")
            return False
        try:
            aln = tree_links.read_alignment(source.path, aligned=source.aligned)
        except OSError as exc:
            QMessageBox.warning(self, "Show in alignment", f"Could not read {source.path.name}:\n{exc}")
            return False
        ids, missing = [], []
        for tip in tips:
            row = aln.row_of(tip)
            if row is None:
                missing.append(tip)
            else:
                name = aln.names[row]
                ids.append(name.split()[0] if name.split() else name)
        if not ids:
            self._update_status(f"Not found in {source.path.name}: {', '.join(missing[:3])}")
            return False
        self.alignmentRequested.emit(source.path, ids, self._path.parent)
        message = f"{source.path.name}: {len(ids)} sequence{'s' if len(ids) != 1 else ''} selected"
        if missing:
            message += f", {len(missing)} not found"
        if source.note:
            message += f" ({source.note})"
        self._update_status(message)
        return True

    def extract_sequences(self, key, *, aligned: bool, path: Path | None = None) -> Path | None:
        """Write the sequences of a clade (or one tip) to a new FASTA next to the tree."""
        from repeatviewer.backends import tree_links

        if not key or self._path is None:
            return None
        tips = list(key)
        got = tree_links.clade_sequences(self._path, tips, aligned=aligned)
        if not got.records:
            QMessageBox.information(
                self, "Extract sequences",
                "No sequences were found for this selection"
                + (f" ({got.note})." if got.note else ".")
                + " RepeatViewer looks for the tree run's input FASTA and alignment.")
            return None
        if path is None:
            default = tree_links.extract_default_path(self._path, len(got.records), aligned=aligned)
            chosen, _ = QFileDialog.getSaveFileName(self, "Extract sequences", str(default),
                                                    "FASTA (*.fa *.fasta)")
            if not chosen:
                return None
            path = Path(chosen)
        try:
            path.write_text(tree_links.fasta_text(got.records), encoding="utf-8")
        except OSError as exc:
            QMessageBox.warning(self, "Extract sequences", f"Could not write {path.name}:\n{exc}")
            return None
        message = f"saved {len(got.records)} sequence(s) to {path.name}"
        if got.note:
            message += f" ({got.note})"
        if got.missing:
            shown = ", ".join(got.missing[:5]) + ("…" if len(got.missing) > 5 else "")
            QMessageBox.warning(self, "Extract sequences",
                                f"{len(got.missing)} sequence(s) were not found and are not "
                                f"in {path.name}: {shown}")
        self._update_status(message)
        self.fileWritten.emit(path)
        return path

    def select_clade(self, key: tuple[str, ...] | None) -> None:
        self._selected = tuple(sorted(key)) if key else None
        self._redraw(refit=False)

    def _is_collapsed(self, key) -> bool:
        """Drawn as a triangle: collapsed by hand or by distance."""
        if key in self._state.collapsed:
            return True
        d = self._display
        return d is not None and any(d.keys[i] == key for i in d.info)

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
        if key:
            n = len(key)
            sub = QMenu("Extract sequence" if n == 1 else f"Extract sequences ({n})", menu)
            sub.setToolTipsVisible(True)
            menu.addMenu(sub)
            act = sub.addAction("Without gaps…", lambda: self.extract_sequences(key, aligned=False))
            act.setToolTip("Your original records with their full headers, gaps removed, "
                           "not rotated or masked.")
            act = sub.addAction("As aligned…", lambda: self.extract_sequences(key, aligned=True))
            source = self._sequence_source()
            act.setToolTip("The rows of the tree's alignment"
                           + (f" ({source.path.name})" if source is not None else "")
                           + ": gaps kept, rotated or masked as in this run.")
            sub.setEnabled(self._path is not None)
            if n > 1:
                act = menu.addAction(f"Show clade in alignment ({n})",
                                     lambda: self.show_in_alignment(list(key)))
                act.setEnabled(source is not None)
        if node is not None and node.is_tip and idx not in d.info:
            menu.addSeparator()
            act = menu.addAction("Show in alignment", lambda: self.show_in_alignment(node.name))
            act.setEnabled(self._sequence_source() is not None)
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
        act.setEnabled(bool(self._state.collapsed or self._state.collapse_distance))
        act = menu.addAction("Remove all colours", self.clear_colors)
        act.setEnabled(bool(self._state.colors))
        menu.addSeparator()
        act = menu.addAction("Undo", self.undo)
        act.setEnabled(bool(self._undo))
        act = menu.addAction("Redo", self.redo)
        act.setEnabled(bool(self._redo))
        menu.addAction("Save tree (.treefile)", self._save_tree_as)
        menu.addAction("Distance matrix…", self._open_distances)
        menu.addAction("Copy Newick", self._copy_newick)

    def _fill_edit_menu(self) -> None:
        menu = self._edit_menu
        menu.clear()
        if self._selected:
            d = self._display
            is_clade = d is not None and self._selected in d.keys
            label = (
                f"Selected clade ({len(self._selected)} tips)"
                if is_clade
                else f"Selected tips ({len(self._selected)})"
            )
            title = menu.addAction(label)
            title.setEnabled(False)
            self._add_clade_actions(menu, self._selected)
        else:
            hint = menu.addAction("Click a branch, or drag over tip names (Cmd/Ctrl adds)")
            hint.setEnabled(False)
        menu.addSeparator()
        self._add_tree_actions(menu)

    def _click_in_selection(self, key: tuple[str, ...] | None) -> bool:
        """True when ``key`` is already covered by the current multi/clade selection."""
        if not key or not self._selected:
            return False
        if key == self._selected:
            return True
        # Tip or collapsed clade inside a rubber-band / clade selection.
        return all(name in self._selected for name in key)

    def _select_for_context_menu(self, index: int) -> tuple[str, ...] | None:
        """Update selection for a right-click; keep multi-select when clicking inside it."""
        key = self._key_for(index)
        if key is None:
            return self._selected
        if not self._click_in_selection(key):
            self._selected = key
            self._redraw(refit=False)
        return self._selected

    def _on_context_menu(self, index: int, global_pos) -> None:
        if self._model is None:
            return
        menu = QMenu(self)
        selected = self._select_for_context_menu(index)
        if selected is not None:
            self._add_clade_actions(menu, selected)
            menu.addSeparator()
        self._add_tree_actions(menu)
        menu.exec(global_pos)

    # ------------------------------------------------------------ edits

    def _push_undo(self) -> None:
        # A new edit: what was undone can no longer be redone.
        if self._model is not None:
            self._undo.append((self._model, self._state.copy(), replace(self._opts)))
            del self._undo[:-50]
        self._redo.clear()

    def _after_edit(self, message: str = "", *, refit: bool = False) -> None:
        self._redraw(refit=refit)
        if message:
            self._update_status(message)
        self._save_view_state()

    def event(self, event) -> bool:  # noqa: N802 (Qt API)
        # Cmd/Ctrl+Z inside the tree is the tree's own undo. Claiming the key
        # here (instead of a QShortcut) keeps it working in every window and
        # stops it competing with the main window's Edit > Undo shortcut, which
        # would make the key do nothing. Same for Redo (Cmd+Shift+Z / Ctrl+Shift+Z).
        if event.type() == QEvent.Type.ShortcutOverride and (
            event.matches(QKeySequence.StandardKey.Undo)
            or event.matches(QKeySequence.StandardKey.Redo)
        ):
            event.accept()
            return True
        return super().event(event)

    def keyPressEvent(self, event) -> None:  # noqa: N802 (Qt API)
        if event.matches(QKeySequence.StandardKey.Undo):
            self.undo()
            event.accept()
            return
        if event.matches(QKeySequence.StandardKey.Redo):
            self.redo()
            event.accept()
            return
        super().keyPressEvent(event)

    def _snapshot(self) -> tuple[tl.TreeModel, TreeViewState, TreeViewOptions]:
        return (self._model, self._state.copy(), replace(self._opts))

    def undo(self) -> None:
        if not self._undo:
            return
        if self._model is not None:
            self._redo.append(self._snapshot())
        self._model, self._state, self._opts = self._undo.pop()
        self._apply_options_to_controls()
        self._after_edit("undone")

    def can_redo(self) -> bool:
        return bool(self._redo)

    def redo(self) -> None:
        if not self._redo:
            return
        if self._model is not None:
            self._undo.append(self._snapshot())
        self._model, self._state, self._opts = self._redo.pop()
        self._apply_options_to_controls()
        self._after_edit("redone")

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
        if self._is_collapsed(key):
            self._state.collapsed.discard(key)
            self._state.expanded.add(key)      # also opens a clade collapsed by distance
            self._after_edit("expanded")
        else:
            self._state.collapsed.add(key)
            self._state.expanded.discard(key)
            self._after_edit("collapsed")

    def expand_all(self) -> None:
        if not self._state.collapsed and not self._state.collapse_distance:
            return
        self._push_undo()
        self._state.collapsed.clear()
        self._state.expanded.clear()
        self._state.collapse_distance = 0.0
        self._apply_options_to_controls()
        self._after_edit("all clades expanded")

    def set_collapse_distance(self, value: float) -> None:
        value = max(0.0, float(value))
        if abs(value - self._state.collapse_distance) < 1e-12:
            return
        self._push_undo()
        self._state.collapse_distance = value
        self._state.expanded.clear()
        self._block = True
        self._collapse_dist.setValue(value)
        self._block = False
        self._after_edit()

    def _on_collapse_distance(self) -> None:
        if not self._block:
            self.set_collapse_distance(self._collapse_dist.value())

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
            self._fill_tiplabel_box()
            self._sync_enabled_controls()
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
        mode = self._opts.tip_label
        if ":" in mode:
            col = mode.split(":", 1)[1]
            if col in self._traits.columns:
                for t, name in names.items():
                    value = self._traits.rows.get(name, {}).get(col, "")
                    if mode.startswith("col:"):
                        display.tip_text[t] = value or name
                    elif value:
                        display.tip_text[t] = f"{name}  {value}"

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

    def report_path(self) -> Path | None:
        """The tree run's report, if this tree came from one."""
        from repeatviewer.backends.tree_links import run_base

        if self._path is None:
            return None
        rb = run_base(self._path)
        if rb is None:
            return None
        base = rb[0]
        for candidate in (self._path.parent / f"{base}_details" / f"{base}_report.html",
                          self._path.parent / f"{base}_report.html"):
            if candidate.is_file():
                return candidate
        return None

    def _open_report(self) -> None:
        report = self.report_path()
        if report is not None:
            QDesktopServices.openUrl(QUrl.fromLocalFile(str(report)))

    def _current_structure(self) -> tl.TreeModel | None:
        """The tree as shown: rerooted, weak branches merged, one clade if chosen."""
        if self._model is None:
            return None
        return tree_structure(self._model, self._state)[0]

    def _copy_newick(self) -> None:
        m = self._current_structure()
        if m is None:
            return
        QApplication.clipboard().setText(te.to_newick(m))
        self._update_status("Newick copied")

    def _open_distances(self) -> None:
        from repeatviewer.views.tree_distances_dialog import TreeDistancesDialog

        m = self._current_structure()
        if m is None or self._path is None:
            return
        dlg = TreeDistancesDialog(self._path, m, self)
        dlg.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        dlg.fileWritten.connect(self.fileWritten)
        self._dist_dialog = dlg
        dlg.show()

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
            self, "Save tree", str(default), "Newick tree (*.treefile *.nwk *.tree)"
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

    def _open_export_dialog(self) -> None:
        """MSA-style Export window: choose Type (SVG/PDF/PNG) and File."""
        if self._display is None or self._path is None:
            QMessageBox.information(self, "Export", "There is no tree to export.")
            return
        from repeatviewer.views.export_tree_dialog import ExportTreeDialog

        dialog = ExportTreeDialog(
            source_path=self._path,
            layout_tag=str(self._opts.layout or "tree"),
            default_format="svg",
            parent=self.window() if self.window() is not None else self,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        out = dialog.chosen_path()
        if out is None:
            return
        try:
            export_tree(
                self._display,
                self._opts,
                out,
                title=self._path.name,
                png_dpi=dialog.dpi(),
            )
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Export", f"Could not write {out.name}:\n{exc}")
            return
        self._update_status(f"exported {out.name}")
        self.fileWritten.emit(out)

    def _export(self, ext: str) -> None:
        """Legacy helper used by tests: write a figure next to the tree."""
        if self._display is None or self._path is None:
            QMessageBox.information(self, "Export", "There is no tree to export.")
            return
        out = export_default_path(self._path, self._opts.layout, ext)
        try:
            export_tree(self._display, self._opts, out, title=self._path.name)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Export", f"Could not write {out.name}:\n{exc}")
            return
        self._update_status(f"exported {out.name}")
        self.fileWritten.emit(out)
