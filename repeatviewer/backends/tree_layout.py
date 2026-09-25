"""Tree layouts for the Preview panel's Tree tab (pure: no Qt).

Reads a Newick (or NEXUS) tree into a small node model and computes drawing
coordinates for three layouts, in the spirit of the Geneious tree viewer:

* **rectangular**: tips in rows, x = distance from the root;
* **circular**: tips around a circle, radius = distance from the root;
* **radial**: unrooted "equal-angle" layout (Felsenstein), each subtree gets
  a wedge proportional to its number of tips.

Coordinates are in layout units: for the rectangular layout x is in branch
length units and y in rows (one per tip); circular and radial use branch
length units on both axes. The Qt widget scales them to pixels.

``.treefile`` is IQ-TREE's extension for a Newick file; ``.nwk``, ``.newick``,
``.tre`` and ``.tree`` are the same format under other names.
"""

from __future__ import annotations

import io
import math
from dataclasses import dataclass, field
from pathlib import Path

TREE_SUFFIXES = frozenset({".treefile", ".nwk", ".newick", ".tre", ".tree"})

LAYOUT_RECTANGULAR = "rectangular"
LAYOUT_CIRCULAR = "circular"
LAYOUT_RADIAL = "radial"
LAYOUTS = (LAYOUT_RECTANGULAR, LAYOUT_CIRCULAR, LAYOUT_RADIAL)

ORDER_FILE = "file"
ORDER_INCREASING = "increasing"
ORDER_DECREASING = "decreasing"
ORDERS = (ORDER_FILE, ORDER_INCREASING, ORDER_DECREASING)

# How branch lengths are drawn (Geneious: "Transform branches")
BRANCHES_PROPORTIONAL = "proportional"   # as in the file
BRANCHES_EQUAL = "equal"                 # cladogram: topology only, tips aligned
BRANCHES_SQRT = "sqrt"                   # square root: one very long branch no longer squeezes the rest
BRANCH_TRANSFORMS = (BRANCHES_PROPORTIONAL, BRANCHES_EQUAL, BRANCHES_SQRT)


class TreeReadError(ValueError):
    """The file is not a readable Newick / NEXUS tree."""


def is_tree_file(path: Path | str) -> bool:
    return Path(path).suffix.lower() in TREE_SUFFIXES


# --------------------------------------------------------------------------
# Model
# --------------------------------------------------------------------------


@dataclass
class Node:
    index: int
    parent: int | None
    children: list[int] = field(default_factory=list)
    name: str = ""        # tip name; internal nodes usually ""
    length: float = 0.0   # branch length to the parent
    support: str = ""     # internal label, e.g. "92/97" (SH-aLRT/UFBoot) or "100"

    @property
    def is_tip(self) -> bool:
        return not self.children


@dataclass
class TreeModel:
    nodes: list[Node]      # preorder; nodes[0] is the root
    has_lengths: bool      # False when the file has no branch lengths
    has_support: bool

    @property
    def root(self) -> Node:
        return self.nodes[0]

    def tips(self) -> list[int]:
        """Tip indices in drawing order (depth-first, children in list order)."""
        out: list[int] = []
        stack = [0]
        while stack:
            node = self.nodes[stack.pop()]
            if node.is_tip:
                out.append(node.index)
            else:
                stack.extend(reversed(node.children))
        return out

    @property
    def n_tips(self) -> int:
        return sum(1 for n in self.nodes if n.is_tip)

    @property
    def n_internal(self) -> int:
        return len(self.nodes) - self.n_tips


def _format_support(value: float) -> str:
    if value == int(value):
        return str(int(value))
    return f"{value:.3g}"


def model_from_phylo(tree) -> TreeModel:
    """Convert a Bio.Phylo tree into a :class:`TreeModel` (preorder)."""
    nodes: list[Node] = []
    has_lengths = False
    has_support = False
    stack: list[tuple[object, int | None]] = [(tree.root, None)]
    while stack:
        clade, parent = stack.pop()
        index = len(nodes)
        is_tip = not clade.clades
        length = clade.branch_length
        if length is not None and parent is not None:
            has_lengths = True
        support = ""
        name = clade.name or ""
        if not is_tip:
            if getattr(clade, "confidence", None) is not None:
                support = _format_support(float(clade.confidence))
            elif name:
                # IQ-TREE writes "SH-aLRT/UFBoot" (e.g. "92/97"), which Newick
                # parsers keep as the internal node's name.
                support, name = name, ""
        if support:
            has_support = True
        nodes.append(Node(index=index, parent=parent, name=name,
                          length=float(length) if length is not None else 0.0,
                          support=support))
        if parent is not None:
            nodes[parent].children.append(index)
        for child in reversed(clade.clades):
            stack.append((child, index))
    if not has_lengths:
        for node in nodes[1:]:
            node.length = 1.0
    return TreeModel(nodes=nodes, has_lengths=has_lengths, has_support=has_support)


def read_tree_text(text: str) -> TreeModel:
    """Parse the first tree in Newick or NEXUS text."""
    from Bio import Phylo

    stripped = text.strip()
    if not stripped:
        raise TreeReadError("The file is empty.")
    fmt = "nexus" if stripped[:6].upper() == "#NEXUS" else "newick"
    try:
        trees = list(Phylo.parse(io.StringIO(stripped), fmt))
    except Exception as exc:  # noqa: BLE001 - Bio raises several types
        raise TreeReadError(f"Not a readable {fmt.upper()} tree: {exc}") from exc
    if not trees:
        raise TreeReadError("No tree found in the file.")
    model = model_from_phylo(trees[0])
    if model.n_tips < 2:
        raise TreeReadError("The tree has fewer than two tips.")
    return model


def read_tree_file(path: Path | str) -> TreeModel:
    return read_tree_text(Path(path).read_text(encoding="utf-8", errors="replace"))


# --------------------------------------------------------------------------
# Order (ladderize) and depths
# --------------------------------------------------------------------------


def _tip_counts(model: TreeModel) -> list[int]:
    counts = [0] * len(model.nodes)
    for node in reversed(model.nodes):  # children come after parents in preorder
        counts[node.index] = 1 if node.is_tip else sum(counts[c] for c in node.children)
    return counts


def ordered(model: TreeModel, order: str) -> TreeModel:
    """Copy of ``model`` with children sorted by clade size (ladderize).

    ``increasing``: small clades first (drawn at the top); ``decreasing``:
    large clades first; ``file``: as written in the file. Ties keep file order.
    """
    nodes = [Node(n.index, n.parent, list(n.children), n.name, n.length, n.support)
             for n in model.nodes]
    if order in (ORDER_INCREASING, ORDER_DECREASING):
        counts = _tip_counts(model)
        reverse = order == ORDER_DECREASING
        for node in nodes:
            node.children.sort(key=lambda c: counts[c], reverse=reverse)
    return TreeModel(nodes=nodes, has_lengths=model.has_lengths, has_support=model.has_support)


def depths(model: TreeModel) -> list[float]:
    """Distance of every node from the root (sum of branch lengths)."""
    out = [0.0] * len(model.nodes)
    for node in model.nodes[1:]:
        out[node.index] = out[node.parent] + max(0.0, node.length)
    return out


def transform_lengths(model: TreeModel, mode: str) -> TreeModel:
    """Copy of ``model`` with branch lengths as drawn for ``mode``.

    ``equal``: every node sits one step left of its tallest child, so all tips
    line up (a cladogram); ``sqrt``: square root of each length. Labels keep
    showing the file's lengths; only the drawing changes.
    """
    nodes = [Node(n.index, n.parent, list(n.children), n.name, n.length, n.support)
             for n in model.nodes]
    if mode == BRANCHES_SQRT:
        for node in nodes[1:]:
            node.length = math.sqrt(max(0.0, node.length))
    elif mode == BRANCHES_EQUAL:
        height = [0] * len(nodes)
        for node in reversed(nodes):
            if node.children:
                height[node.index] = 1 + max(height[c] for c in node.children)
        for node in nodes[1:]:
            node.length = float(height[node.parent] - height[node.index])
    return TreeModel(nodes=nodes, has_lengths=model.has_lengths, has_support=model.has_support)


def transform_length(value: float, mode: str) -> float:
    """One length (e.g. a collapsed clade's depth) as drawn for ``mode``."""
    if mode == BRANCHES_SQRT:
        return math.sqrt(max(0.0, value))
    if mode == BRANCHES_EQUAL:
        return 1.0
    return value


# --------------------------------------------------------------------------
# Layouts
# --------------------------------------------------------------------------


@dataclass
class Arc:
    radius: float
    start_deg: float   # counter-clockwise from +x, as in QPainterPath.arcTo
    span_deg: float


@dataclass
class Layout:
    kind: str
    xy: list[tuple[float, float]]          # node positions
    lines: list[tuple[float, float, float, float]]
    arcs: list[Arc]                        # circular layout only, centred on (0, 0)
    tip_angle: dict[int, float]            # label direction per tip, degrees
    branch_mid: dict[int, tuple[float, float]]  # middle of each node's own branch
    extent: float                          # max root-to-tip distance
    n_rows: int                            # rectangular: number of tip rows
    node_angle: dict[int, float] = field(default_factory=dict)  # circular / radial, degrees
    inner: float = 0.0                     # circular: empty radius left at the centre


def _rectangular(model: TreeModel, depth: list[float]) -> Layout:
    y = [0.0] * len(model.nodes)
    for row, tip in enumerate(model.tips()):
        y[tip] = float(row)
    for node in reversed(model.nodes):
        if not node.is_tip:
            y[node.index] = (y[node.children[0]] + y[node.children[-1]]) / 2.0
    xy = [(depth[i], y[i]) for i in range(len(model.nodes))]
    lines = []
    mid = {}
    for node in model.nodes:
        if node.parent is not None:
            px = depth[node.parent]
            lines.append((px, y[node.index], depth[node.index], y[node.index]))
            mid[node.index] = ((px + depth[node.index]) / 2.0, y[node.index])
        if not node.is_tip:
            x = depth[node.index]
            lines.append((x, y[node.children[0]], x, y[node.children[-1]]))
    tips = model.tips()
    return Layout(LAYOUT_RECTANGULAR, xy, lines, [], {t: 0.0 for t in tips}, mid,
                  max(depth), len(tips))


def _polar(r: float, deg: float) -> tuple[float, float]:
    rad = math.radians(deg)
    # Screen y grows downwards: negate so angles run counter-clockwise on screen.
    return (r * math.cos(rad), -r * math.sin(rad))


def _circular(model: TreeModel, depth: list[float], inner: float = 0.0) -> Layout:
    tips = model.tips()
    n = len(tips)
    angle = [0.0] * len(model.nodes)
    for i, tip in enumerate(tips):
        angle[tip] = 360.0 * i / n
    for node in reversed(model.nodes):
        if not node.is_tip:
            angle[node.index] = (angle[node.children[0]] + angle[node.children[-1]]) / 2.0
    # ``inner`` moves the root off the centre, so the branches near it open up.
    depth = [d + inner for d in depth]
    xy = [_polar(depth[i], angle[i]) for i in range(len(model.nodes))]
    lines = []
    arcs = []
    mid = {}
    for node in model.nodes:
        i = node.index
        if node.parent is not None:
            r0 = depth[node.parent]
            x0, y0 = _polar(r0, angle[i])
            lines.append((x0, y0, xy[i][0], xy[i][1]))
            mid[i] = _polar((r0 + depth[i]) / 2.0, angle[i])
        if not node.is_tip:
            a0 = angle[node.children[0]]
            a1 = angle[node.children[-1]]
            if depth[i] > 0 and a1 > a0:
                arcs.append(Arc(depth[i], a0, a1 - a0))
    return Layout(LAYOUT_CIRCULAR, xy, lines, arcs, {t: angle[t] for t in tips}, mid,
                  max(depth), n, dict(enumerate(angle)), inner)


def _radial(model: TreeModel) -> Layout:
    """Equal-angle unrooted layout: wedges proportional to tip counts."""
    counts = _tip_counts(model)
    xy = [(0.0, 0.0)] * len(model.nodes)
    direction = [0.0] * len(model.nodes)
    wedge = {0: (0.0, 360.0)}
    for node in model.nodes:  # preorder: a parent's wedge is set before its children
        start, span = wedge[node.index]
        total = counts[node.index]
        cursor = start
        for child in node.children:
            share = span * counts[child] / total
            deg = cursor + share / 2.0
            direction[child] = deg
            dx, dy = _polar(max(0.0, model.nodes[child].length), deg)
            px, py = xy[node.index]
            xy[child] = (px + dx, py + dy)
            wedge[child] = (cursor, share)
            cursor += share
    lines = []
    mid = {}
    for node in model.nodes[1:]:
        px, py = xy[node.parent]
        x, y = xy[node.index]
        lines.append((px, py, x, y))
        mid[node.index] = ((px + x) / 2.0, (py + y) / 2.0)
    tips = model.tips()
    extent = max(math.hypot(x, y) for x, y in xy)
    return Layout(LAYOUT_RADIAL, xy, lines, [], {t: direction[t] for t in tips}, mid,
                  extent, len(tips), dict(enumerate(direction)))


def compute_layout(model: TreeModel, kind: str = LAYOUT_RECTANGULAR,
                   extra_depth: dict[int, float] | None = None, *,
                   inner_fraction: float = 0.0) -> Layout:
    """Coordinates for ``kind``. ``extra_depth`` reserves room beyond some tips
    (collapsed clades drawn as triangles) so the drawing scale includes them.
    ``inner_fraction`` (circular): empty centre as a share of the tree depth."""
    if kind == LAYOUT_RADIAL:
        layout = _radial(model)
        if extra_depth:
            layout.extent = max(layout.extent, max(
                math.hypot(*layout.xy[i]) + d for i, d in extra_depth.items()))
        return layout
    depth = depths(model)
    if kind == LAYOUT_CIRCULAR:
        inner = max(0.0, inner_fraction) * (max(depth) or 1.0)
        layout = _circular(model, depth, inner)
    else:
        layout = _rectangular(model, depth)
    if extra_depth:
        layout.extent = max(layout.extent, max(depth[i] + layout.inner + d
                                               for i, d in extra_depth.items()))
    return layout


def label_rotation(angle_deg: float) -> tuple[float, bool]:
    """Screen rotation for a tip label pointing along ``angle_deg``.

    Returns (rotation in degrees, clockwise as QGraphicsItem.setRotation expects;
    flipped). Labels on the left half are turned 180° so they never read
    upside down; ``flipped`` then means the text ends at the tip instead of
    starting there.
    """
    a = angle_deg % 360.0
    flipped = 90.0 < a < 270.0
    rotation = -a + (180.0 if flipped else 0.0)
    return rotation, flipped


# --------------------------------------------------------------------------
# Labels and scale bar
# --------------------------------------------------------------------------


def branch_length_text(length: float, digits: int = 4) -> str:
    """Branch length with ``digits`` significant digits (Geneious style)."""
    if length == 0:
        return "0"
    return f"{length:.{digits}g}"


def nice_scale_length(extent: float, fraction: float = 0.2) -> float:
    """A round scale-bar length (1, 2 or 5 × 10^k) near ``fraction`` of ``extent``."""
    if extent <= 0 or not math.isfinite(extent):
        return 0.0
    target = extent * fraction
    power = 10.0 ** math.floor(math.log10(target))
    for step in (5.0, 2.0, 1.0):
        if step * power <= target:
            return step * power
    return power


def scale_text(value: float) -> str:
    return f"{value:g}"
