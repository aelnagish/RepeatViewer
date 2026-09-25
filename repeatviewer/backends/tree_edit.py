"""Tree edits for the Tree tab (pure: no Qt).

Operations on :class:`repeatviewer.backends.tree_layout.TreeModel`, each
returning a new model (the input is never changed, so the viewer can undo):

* :func:`reroot` on a branch, :func:`midpoint_root`;
* :func:`swap_children` (rotate a node; the tree itself is unchanged);
* :func:`collapse_weak`: contract branches with support below a threshold
  into polytomies;
* :func:`subtree`: one clade on its own;
* :func:`display_model`: collapsed clades drawn as single "tips".

Support values belong to branches (splits), but Newick stores them on nodes.
Rerooting therefore goes through an unrooted edge list, so every support stays
on its own branch (a common mistake when rerooting Newick trees).

Clades are identified by their tip names (:func:`clade_key`), which survive
reordering and reloading, so view settings can refer to them.
"""

from __future__ import annotations

from dataclasses import dataclass

from repeatviewer.backends.tree_layout import Node, TreeModel

# --------------------------------------------------------------------------
# Helpers
# --------------------------------------------------------------------------


def _copy(model: TreeModel, nodes: list[Node]) -> TreeModel:
    return TreeModel(nodes=nodes, has_lengths=model.has_lengths,
                     has_support=any(n.support for n in nodes))


def tip_sets(model: TreeModel) -> list[frozenset[str]]:
    """Tip names below every node."""
    out: list[frozenset[str]] = [frozenset()] * len(model.nodes)
    for node in reversed(model.nodes):
        if node.is_tip:
            out[node.index] = frozenset([node.name])
        else:
            acc: set[str] = set()
            for c in node.children:
                acc |= out[c]
            out[node.index] = frozenset(acc)
    return out


def clade_key(model: TreeModel, node: int) -> tuple[str, ...]:
    """Sorted tip names of the clade below ``node`` (stable across edits)."""
    return tuple(sorted(tip_sets(model)[node]))


def find_clade(model: TreeModel, key) -> int | None:
    wanted = frozenset(key)
    for i, tips in enumerate(tip_sets(model)):
        if tips == wanted:
            return i
    return None


def _rebuild(root, children_of, make) -> list[Node]:
    """Build a preorder node list. ``children_of(key)`` lists child keys,
    ``make(key, index, parent_index)`` returns the Node without children."""
    nodes: list[Node] = []
    stack = [(root, None)]
    while stack:
        key, parent = stack.pop()
        index = len(nodes)
        node = make(key, index, parent)
        nodes.append(node)
        if parent is not None:
            nodes[parent].children.append(index)
        for child in reversed(children_of(key)):
            stack.append((child, index))
    return nodes


# --------------------------------------------------------------------------
# Simple edits
# --------------------------------------------------------------------------


def swap_children(model: TreeModel, node: int) -> TreeModel:
    """Reverse the order of ``node``'s branches (same tree, other drawing)."""
    nodes = [Node(n.index, n.parent, list(n.children), n.name, n.length, n.support)
             for n in model.nodes]
    nodes[node].children.reverse()
    # keep the list in preorder so indices stay meaningful for drawing
    keyed = {n.index: n for n in nodes}
    return _copy(model, _rebuild(
        0, lambda k: keyed[k].children,
        lambda k, i, p: Node(i, p, [], keyed[k].name, keyed[k].length, keyed[k].support)))


def subtree(model: TreeModel, node: int) -> TreeModel:
    """The clade below ``node`` as a tree of its own."""
    src = model.nodes
    return _copy(model, _rebuild(
        node, lambda k: src[k].children,
        lambda k, i, p: Node(i, p, [], src[k].name, 0.0 if p is None else src[k].length,
                             src[k].support)))


# --------------------------------------------------------------------------
# Support values and collapsing
# --------------------------------------------------------------------------


def _numbers(label: str) -> list[float]:
    out = []
    for part in str(label or "").split("/"):
        try:
            out.append(float(part))
        except ValueError:
            return []
    return out


def support_scale(model: TreeModel) -> float:
    """100 when supports are fractions (FastTree 0–1), else 1."""
    values = [v for n in model.nodes for v in _numbers(n.support)]
    return 100.0 if values and max(values) <= 1.0 else 1.0


def support_value(label: str, scale: float = 1.0) -> float | None:
    """The value compared with a collapse threshold, on a 0–100 scale.

    IQ-TREE "SH-aLRT/UFBoot" uses the UFBoot value (the second one)."""
    vals = _numbers(label)
    if not vals:
        return None
    return vals[-1] * scale


def collapse_weak(model: TreeModel, threshold: float) -> tuple[TreeModel, int]:
    """Contract every internal branch whose support is below ``threshold``.

    The node's branches join its parent (a polytomy); their lengths grow by the
    contracted branch, so root-to-tip distances stay the same. Branches without
    a support value are kept. Returns (new model, number contracted).
    """
    scale = support_scale(model)
    src = model.nodes
    weak = [False] * len(src)
    for n in src:
        if n.parent is not None and not n.is_tip:
            value = support_value(n.support, scale)
            weak[n.index] = value is not None and value < threshold
    new_parent: dict[int, int | None] = {}
    extra = [0.0] * len(src)
    kept: list[int] = []
    for n in src:  # preorder: parents first
        if n.parent is None:
            new_parent[n.index] = None
        else:
            p = n.parent
            if weak[p]:
                new_parent[n.index] = new_parent[p]
                extra[n.index] = extra[p] + src[p].length
            else:
                new_parent[n.index] = p
        if not weak[n.index]:
            kept.append(n.index)
    kids: dict[int, list[int]] = {k: [] for k in kept}
    for k in kept:
        p = new_parent[k]
        if p is not None:
            kids[p].append(k)
    nodes = _rebuild(0, lambda k: kids[k],
                     lambda k, i, p: Node(i, p, [], src[k].name, src[k].length + extra[k],
                                          src[k].support))
    return _copy(model, nodes), sum(weak)


# --------------------------------------------------------------------------
# Rerooting
# --------------------------------------------------------------------------


@dataclass
class _Edge:
    length: float
    support: str


def _unrooted(model: TreeModel) -> tuple[dict[int, list[int]], dict[frozenset, _Edge]]:
    """Adjacency and edges; a two-child root is dissolved into one branch."""
    adj: dict[int, list[int]] = {n.index: [] for n in model.nodes}
    edges: dict[frozenset, _Edge] = {}

    def link(a: int, b: int, length: float, support: str) -> None:
        adj[a].append(b)
        adj[b].append(a)
        edges[frozenset((a, b))] = _Edge(length, support)

    root = model.root
    for n in model.nodes[1:]:
        if n.parent == 0 and len(root.children) == 2:
            continue
        link(n.index, n.parent, n.length, "" if n.is_tip else n.support)
    if len(root.children) == 2:
        a, b = (model.nodes[c] for c in root.children)
        label = (a.support if not a.is_tip else "") or (b.support if not b.is_tip else "")
        link(a.index, b.index, a.length + b.length, label)
        del adj[0]
    return adj, edges


def _root_on_edge(model: TreeModel, adj, edges, u: int, v: int, from_u: float) -> TreeModel:
    """New root on branch u–v, ``from_u`` along it from u."""
    src = model.nodes
    edge = edges[frozenset((u, v))]
    from_u = min(max(0.0, from_u), edge.length)

    def count_side(start: int, avoid: int) -> int:
        """Tips reachable from ``start`` without crossing ``avoid``."""
        seen, stack, total = {avoid}, [start], 0
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            total += 1 if src[x].is_tip else 0
            stack.extend(adj[x])
        return total

    u_tips, v_tips = count_side(u, v), count_side(v, u)
    # The branch's support describes one split: show it once, on the larger side.
    u_label = edge.support if (not src[u].is_tip and u_tips >= v_tips) else ""
    v_label = edge.support if (not src[v].is_tip and not u_label) else ""
    ROOT = -1

    def children_of(key):
        node, came_from = key
        if node == ROOT:
            return [(u, ROOT), (v, ROOT)]
        return [(x, node) for x in adj[node] if x != came_from
                and not (came_from == ROOT and x in (u, v) and {node, x} == {u, v})]

    def make(key, index, parent):
        node, came_from = key
        if node == ROOT:
            return Node(index, None, [], "", 0.0, "")
        if came_from == ROOT:
            length = from_u if node == u else edge.length - from_u
            label = u_label if node == u else v_label
        else:
            e = edges[frozenset((node, came_from))]
            length, label = e.length, e.support
        return Node(index, parent, [], src[node].name, length, label)

    return _copy(model, _rebuild((ROOT, None), children_of, make))


def reroot(model: TreeModel, node: int, position: float = 0.5) -> TreeModel:
    """Root on the branch above ``node`` (``position`` = fraction from ``node``)."""
    if node == 0:
        return model
    adj, edges = _unrooted(model)
    parent = model.nodes[node].parent
    if parent == 0 and 0 not in adj:  # old two-child root was dissolved
        other = [c for c in model.root.children if c != node][0]
        u, v = node, other
    else:
        u, v = node, parent
    length = edges[frozenset((u, v))].length
    return _root_on_edge(model, adj, edges, u, v, length * position)


def _farthest(adj, edges, start: int) -> tuple[int, float, dict[int, int]]:
    dist, prev, stack = {start: 0.0}, {start: -1}, [start]
    while stack:
        x = stack.pop()
        for y in adj[x]:
            if y not in dist:
                dist[y] = dist[x] + edges[frozenset((x, y))].length
                prev[y] = x
                stack.append(y)
    far = max(dist, key=dist.get)
    return far, dist[far], prev


def midpoint_root(model: TreeModel) -> TreeModel:
    """Root halfway along the longest tip-to-tip path."""
    adj, edges = _unrooted(model)
    any_tip = next(n.index for n in model.nodes if n.is_tip)
    a, _d, _p = _farthest(adj, edges, any_tip)
    b, total, prev = _farthest(adj, edges, a)
    half = total / 2.0
    # walk from b back to a until half the path is covered
    walked, x = 0.0, b
    while prev[x] != -1:
        y = prev[x]
        length = edges[frozenset((x, y))].length
        if walked + length >= half:
            return _root_on_edge(model, adj, edges, x, y, half - walked)
        walked += length
        x = y
    return model


# --------------------------------------------------------------------------
# Collapsed clades (display only)
# --------------------------------------------------------------------------


@dataclass
class CollapsedInfo:
    n_tips: int
    depth: float       # deepest tip below the collapsed node, from that node
    first_tip: str
    key: tuple[str, ...] = ()   # clade key of the collapsed node


def display_model(model: TreeModel, collapsed: set[int]) -> tuple[TreeModel, dict[int, CollapsedInfo]]:
    """Model where each collapsed node is a single tip (the viewer draws a
    triangle). Returns the model and, per new node index, what it hides."""
    src = model.nodes
    if not collapsed:
        return model, {}
    depth_below = [0.0] * len(src)
    counts = [0] * len(src)
    first = [""] * len(src)
    for n in reversed(src):
        if n.is_tip:
            counts[n.index], first[n.index] = 1, n.name
        else:
            counts[n.index] = sum(counts[c] for c in n.children)
            depth_below[n.index] = max(depth_below[c] + src[c].length for c in n.children)
            first[n.index] = first[n.children[0]]
    info: dict[int, CollapsedInfo] = {}
    sets = tip_sets(model)

    def children_of(k):
        return [] if k in collapsed else src[k].children

    def make(k, i, p):
        if k in collapsed and not src[k].is_tip:
            info[i] = CollapsedInfo(counts[k], depth_below[k], first[k], tuple(sorted(sets[k])))
            return Node(i, p, [], f"{first[k]} + {counts[k] - 1} more", src[k].length, src[k].support)
        return Node(i, p, [], src[k].name, src[k].length, src[k].support)

    return _copy(model, _rebuild(0, children_of, make)), info


# --------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------

_NEEDS_QUOTES = set(" \t():;,[]'")


def _quote(name: str) -> str:
    if name and _NEEDS_QUOTES & set(name):
        return "'" + name.replace("'", "''") + "'"
    return name


def to_newick(model: TreeModel) -> str:
    """Newick text with branch lengths and support labels."""
    parts: list[str] = [""] * len(model.nodes)
    for node in reversed(model.nodes):
        if node.is_tip:
            text = _quote(node.name)
        else:
            text = "(" + ",".join(parts[c] for c in node.children) + ")" + _quote(node.support)
        if node.parent is not None and model.has_lengths:
            text += f":{node.length:.10g}"
        parts[node.index] = text
    return parts[0] + ";"
