"""Tree tab layout maths (no Qt): parsing, ordering, the three layouts."""

from __future__ import annotations

import math

import pytest

from repeatviewer.backends import tree_layout as tl

IQTREE = "((A:0.1,B:0.2)92/97:0.05,(C:0.3,(D:0.1,E:0.1)100:0.2)88/99:0.1,F:0.4);"


def _model(text=IQTREE):
    return tl.read_tree_text(text)


def _by_name(model):
    return {n.name: n for n in model.nodes if n.name}


def test_tree_suffixes():
    for name in ("x.fa_tree.treefile", "x.nwk", "x.newick", "x.tre", "x.TREE"):
        assert tl.is_tree_file(name)
    assert not tl.is_tree_file("x.fa")


def test_reads_iqtree_support_labels():
    m = _model()
    assert m.n_tips == 6
    assert m.has_lengths and m.has_support
    supports = sorted(n.support for n in m.nodes if n.support)
    assert supports == ["100", "88/99", "92/97"]
    # internal support labels never become tip names
    assert set(_by_name(m)) == set("ABCDEF")


def test_no_branch_lengths_become_unit_lengths():
    m = _model("((A,B),(C,D));")
    assert not m.has_lengths
    assert all(n.length == 1.0 for n in m.nodes[1:])


def test_nexus_is_read():
    m = _model("#NEXUS\nbegin trees;\ntree t1 = ((A:1,B:1):1,C:2);\nend;\n")
    assert sorted(_by_name(m)) == ["A", "B", "C"]


@pytest.mark.parametrize("text", ["", "not a tree", "(A);"])
def test_bad_input_raises(text):
    with pytest.raises(tl.TreeReadError):
        tl.read_tree_text(text)


def test_order_increasing_and_decreasing():
    m = _model("((A,(B,C)),D);")
    names = lambda mm: [mm.nodes[t].name for t in mm.tips()]  # noqa: E731
    assert names(tl.ordered(m, tl.ORDER_FILE)) == ["A", "B", "C", "D"]
    assert names(tl.ordered(m, tl.ORDER_INCREASING)) == ["D", "A", "B", "C"]
    assert names(tl.ordered(m, tl.ORDER_DECREASING)) == ["B", "C", "A", "D"]
    # ordering copies; the original is unchanged
    assert names(m) == ["A", "B", "C", "D"]


def test_rectangular_positions():
    m = _model()
    lay = tl.compute_layout(m, tl.LAYOUT_RECTANGULAR)
    n = _by_name(m)
    ax, ay = lay.xy[n["A"].index]
    assert ax == pytest.approx(0.15)          # 0.05 + 0.1
    assert ay == 0.0                          # first row
    assert lay.xy[n["F"].index] == (pytest.approx(0.4), 5.0)
    assert lay.n_rows == 6
    assert lay.extent == pytest.approx(0.4)   # root → C, D and F are all 0.4
    # parent sits halfway between its first and last child
    parent = m.nodes[n["A"].parent]
    assert lay.xy[parent.index][1] == pytest.approx(0.5)
    # one horizontal per non-root node + one vertical per internal node
    assert len(lay.lines) == (len(m.nodes) - 1) + m.n_internal


def test_circular_tips_on_their_radius():
    m = _model()
    lay = tl.compute_layout(m, tl.LAYOUT_CIRCULAR)
    d = tl.depths(m)
    for tip in m.tips():
        x, y = lay.xy[tip]
        assert math.hypot(x, y) == pytest.approx(d[tip])
    angles = [lay.tip_angle[t] for t in m.tips()]
    assert angles == sorted(angles) and angles[0] == 0.0 and angles[-1] < 360.0
    assert all(a.span_deg > 0 for a in lay.arcs)


def test_radial_branch_lengths_are_true_distances():
    m = _model()
    lay = tl.compute_layout(m, tl.LAYOUT_RADIAL)
    for node in m.nodes[1:]:
        (px, py), (x, y) = lay.xy[node.parent], lay.xy[node.index]
        assert math.hypot(x - px, y - py) == pytest.approx(node.length)


def test_label_rotation_never_upside_down():
    assert tl.label_rotation(0.0) == (0.0, False)
    rot, flipped = tl.label_rotation(180.0)
    assert flipped and rot == pytest.approx(0.0)
    rot, flipped = tl.label_rotation(45.0)
    assert not flipped and rot == pytest.approx(-45.0)


@pytest.mark.parametrize("extent,expected", [(0.6, 0.1), (1.0, 0.2), (0.02, 0.002), (37, 5)])
def test_nice_scale_length(extent, expected):
    assert tl.nice_scale_length(extent) == pytest.approx(expected)


def test_branch_length_text():
    assert tl.branch_length_text(0.013912) == "0.01391"
    assert tl.branch_length_text(0.0) == "0"
