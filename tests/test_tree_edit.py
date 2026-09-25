"""Tree tab edits (pure): reroot, midpoint, swap, collapse, subtree, collapsed display."""

from __future__ import annotations

import pytest

from repeatviewer.backends import tree_edit as te
from repeatviewer.backends import tree_layout as tl

IQTREE = "((A:0.1,B:0.2)92/97:0.05,(C:0.3,(D:0.1,E:0.1)100/100:0.2)88/99:0.1,F:0.4);"


def _m(text=IQTREE):
    return tl.read_tree_text(text)


def _names(m):
    return [m.nodes[t].name for t in m.tips()]


def _by_name(m):
    return {n.name: n.index for n in m.nodes if n.name}


def _labelled_splits(m):
    """Unrooted split → support label, the way a Newick reader should see it."""
    sets = te.tip_sets(m)
    taxa = sets[0]
    first = min(taxa)
    out = {}
    for n in m.nodes[1:]:
        s = sets[n.index]
        if first in s:
            s = taxa - s
        if 1 < len(s) < len(taxa) - 1:
            if n.support or s not in out:
                out[s] = n.support or out.get(s, "")
    return out


def _tip_distances(m):
    d = tl.depths(m)
    sets = te.tip_sets(m)
    names = {n.name: n.index for n in m.nodes if n.is_tip}
    out = {}
    for a in names:
        for b in names:
            if a < b:
                ia, ib = names[a], names[b]
                # deepest common ancestor
                anc = max((n.index for n in m.nodes if a in sets[n.index] and b in sets[n.index]),
                          key=lambda i: d[i])
                out[(a, b)] = round(d[ia] + d[ib] - 2 * d[anc], 9)
    return out


def test_reroot_keeps_supports_on_their_branches():
    m = _m()
    before = _labelled_splits(m)
    for target in ("C", "D", "B", "F"):
        r = te.reroot(m, _by_name(m)[target])
        assert _labelled_splits(r) == before, target
        assert sorted(_names(r)) == sorted(_names(m))


def test_reroot_keeps_tip_to_tip_distances():
    m = _m()
    r = te.reroot(m, _by_name(m)["D"])
    assert _tip_distances(r) == _tip_distances(m)


def test_reroot_puts_the_target_next_to_the_root():
    m = _m()
    r = te.reroot(m, _by_name(m)["F"])
    root_children = [r.nodes[c] for c in r.root.children]
    assert "F" in [c.name for c in root_children]
    assert sum(c.length for c in root_children) == pytest.approx(0.4)


def test_reroot_on_a_two_child_root_tree():
    # A two-child root splits one branch in two; IQ-TREE labels it once.
    m = _m("((A:1,B:1)90:1,(C:1,D:1):1,);".replace(",)", ")"))
    r = te.reroot(m, _by_name(m)["A"])
    assert _labelled_splits(r) == _labelled_splits(m)
    assert _tip_distances(r) == _tip_distances(m)


def test_midpoint_root_balances_the_longest_path():
    m = _m("((A:1,B:1):1,(C:1,D:9):1);")
    r = te.midpoint_root(m)
    d = tl.depths(r)
    tips = {r.nodes[t].name: d[t] for t in r.tips()}
    assert tips["D"] == pytest.approx(6.0)          # longest path A–D = 12, half = 6
    assert max(tips["A"], tips["B"]) == pytest.approx(6.0)
    assert _tip_distances(r) == _tip_distances(m)


def test_swap_children_changes_drawing_order_only():
    m = _m("((A,B),(C,D));")
    root_first_child = m.root.children[0]
    s = te.swap_children(m, 0)
    assert _names(s) == ["C", "D", "A", "B"]
    s2 = te.swap_children(m, root_first_child)
    assert _names(s2) == ["B", "A", "C", "D"]
    assert _labelled_splits(s) == _labelled_splits(m)


def test_collapse_weak_uses_ufboot_and_keeps_distances():
    m = _m()
    c, n = te.collapse_weak(m, 98)       # 92/97 goes, 88/99 and 100/100 stay
    assert n == 1
    root_kids = {c.nodes[k].name for k in c.root.children}
    assert {"A", "B", "F"} <= root_kids   # A and B now hang from the root
    # the contracted branch is added to its children: root-to-tip distances stay
    old = {m.nodes[t].name: d for t, d in zip(m.tips(), [tl.depths(m)[t] for t in m.tips()])}
    new = {c.nodes[t].name: d for t, d in zip(c.tips(), [tl.depths(c)[t] for t in c.tips()])}
    assert new == pytest.approx(old)
    assert te.collapse_weak(m, 0)[1] == 0


def test_collapse_weak_scales_fasttree_support():
    m = _m("((A:1,B:1)0.95:1,(C:1,D:1)0.5:1,E:1);")
    _c, n = te.collapse_weak(m, 90)
    assert n == 1


def test_subtree_and_clade_keys():
    m = _m()
    cd = m.nodes[_by_name(m)["C"]].parent
    key = te.clade_key(m, cd)
    assert key == ("C", "D", "E")
    s = te.subtree(m, cd)
    assert sorted(_names(s)) == ["C", "D", "E"] and s.root.length == 0.0
    assert te.find_clade(m, key) == cd
    assert te.find_clade(te.reroot(m, _by_name(m)["A"]), ("D", "E")) is not None


def test_display_model_collapses_to_one_tip():
    m = _m()
    cd = m.nodes[_by_name(m)["C"]].parent
    d, info = te.display_model(m, {cd})
    assert d.n_tips == 4                              # A, B, F + one collapsed clade
    (idx, ci), = info.items()
    assert ci.n_tips == 3 and ci.first_tip == "C"
    assert ci.depth == pytest.approx(0.3)             # C (0.3) vs D/E (0.2 + 0.1)
    assert d.nodes[idx].is_tip and d.nodes[idx].name == "C + 2 more"


def _random_newick(n_tips, seed):
    import random
    rng = random.Random(seed)
    nodes = [f"t{i}:{rng.uniform(0.01, 0.5):.4f}" for i in range(n_tips)]
    k = 0
    while len(nodes) > 3:
        a = nodes.pop(rng.randrange(len(nodes)))
        b = nodes.pop(rng.randrange(len(nodes)))
        k += 1
        nodes.append(f"({a},{b}){k}/{rng.randint(50, 100)}:{rng.uniform(0.01, 0.3):.4f}")
    return "(" + ",".join(nodes) + ");"


@pytest.mark.parametrize("seed", [1, 2, 3])
def test_reroot_anywhere_on_random_trees(seed):
    m = _m(_random_newick(60, seed))
    splits, dists = _labelled_splits(m), _tip_distances(m)
    for node in range(1, len(m.nodes)):
        r = te.reroot(m, node, 0.3)
        assert _labelled_splits(r) == splits
        assert _tip_distances(r) == pytest.approx(dists)
        # rerooting the result again keeps everything too
        r2 = te.reroot(r, len(r.nodes) - 1)
        assert _labelled_splits(r2) == splits
    mid = te.midpoint_root(m)
    assert _labelled_splits(mid) == splits


def test_newick_round_trip_keeps_names_support_and_lengths():
    m = _m("((\'CroSat 1 mono\':0.1,B:0.2)92/97:0.05,(C:0.3,D:0.1)88/99:0.1,F:0.4);")
    back = tl.read_tree_text(te.to_newick(m))
    assert sorted(_names(back)) == sorted(_names(m))
    assert "CroSat 1 mono" in _names(back)
    assert _labelled_splits(back) == _labelled_splits(m)
    assert _tip_distances(back) == pytest.approx(_tip_distances(m))


def test_layout_extent_includes_collapsed_triangles():
    m = _m(IQTREE.replace("F:0.4", "F:0.1"))   # the collapsed clade is the deepest part
    cd = m.nodes[_by_name(m)["C"]].parent
    d, info = te.display_model(m, {cd})
    (idx, ci), = info.items()
    plain = tl.compute_layout(d, tl.LAYOUT_RECTANGULAR)
    room = tl.compute_layout(d, tl.LAYOUT_RECTANGULAR, {idx: ci.depth})
    assert plain.extent == pytest.approx(0.25)   # B; the triangle is not counted
    assert room.extent == pytest.approx(0.4)    # 0.1 to the node + 0.3 inside it
    assert ci.key == ("C", "D", "E")
