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
    assert tl.branch_length_text(0.02767, decimals=2) == "0.03"
    assert tl.branch_length_text(0.2088, decimals=2) == "0.21"
    assert tl.branch_length_text(0.00154, decimals=2) == "0.00"
    assert tl.branch_length_text(0.0, decimals=2) == "0"


def test_branch_transforms_and_inner_space():
    m = tl.read_tree_text("((A:0.1,B:0.2):0.05,(C:0.3,(D:0.1,E:0.1):0.2):0.1,F:4.0);")
    eq = tl.transform_lengths(m, tl.BRANCHES_EQUAL)
    d = tl.depths(eq)
    assert len({round(d[t], 9) for t in eq.tips()}) == 1       # cladogram: tips aligned
    sq = tl.transform_lengths(m, tl.BRANCHES_SQRT)
    f = next(n for n in sq.nodes if n.name == "F")
    assert abs(f.length - 2.0) < 1e-9
    assert next(n for n in m.nodes if n.name == "F").length == 4.0   # original untouched
    assert tl.transform_length(4.0, tl.BRANCHES_SQRT) == 2.0
    plain = tl.compute_layout(m, tl.LAYOUT_CIRCULAR)
    holed = tl.compute_layout(m, tl.LAYOUT_CIRCULAR, inner_fraction=0.25)
    assert plain.inner == 0 and abs(holed.inner - 1.0) < 1e-9       # 25 % of depth 4.0
    assert abs(holed.extent - plain.extent - 1.0) < 1e-9
    root_r = (holed.xy[0][0] ** 2 + holed.xy[0][1] ** 2) ** 0.5
    assert abs(root_r - 1.0) < 1e-9


def test_orient_rotates_and_mirrors():
    m = tl.read_tree_text("((A:1,B:1):1,C:2);")
    circ = tl.compute_layout(m, tl.LAYOUT_CIRCULAR)
    a = next(n.index for n in m.nodes if n.name == "A")
    turned = tl.orient(circ, 90)
    assert abs(turned.tip_angle[a] - (circ.tip_angle[a] + 90)) < 1e-9
    x, y = circ.xy[a]
    tx, ty = turned.xy[a]
    assert abs(tx ** 2 + ty ** 2 - (x ** 2 + y ** 2)) < 1e-9          # same radius
    rect = tl.compute_layout(m, tl.LAYOUT_RECTANGULAR)
    flipped = tl.orient(rect, 45, mirror=True)                         # angle ignored
    assert flipped.mirrored and flipped.xy[a] == (-rect.xy[a][0], rect.xy[a][1])
    assert tl.orient(rect) is rect


def test_other_programs_tree_names_are_trees():
    from pathlib import Path
    from repeatviewer.backends.tree_layout import is_tree_file

    for name in ("x.contree", "RAxML_bestTree.run1", "RAxML_bipartitions.run1",
                 "RAxML_result.run1", "a.treefile"):
        assert is_tree_file(Path(name)), name
    assert not is_tree_file(Path("RAxML_info.run1"))
    assert not is_tree_file(Path("notes.txt"))


def test_bom_and_crlf_do_not_leak_into_labels(tmp_path):
    from repeatviewer.backends.tree_layout import read_tree_file

    path = tmp_path / "t.nwk"
    path.write_bytes(b"\xef\xbb\xbf(A:0.1,B:0.2,(C:0.3,D:0.4)95:0.05);\r\n")
    model = read_tree_file(path)
    assert [n.support for n in model.nodes if n.children and n.support] == ["95"]


def test_nexus_quotes_removed_from_tip_names():
    from repeatviewer.backends.tree_layout import read_tree_text

    text = ("#NEXUS\nbegin trees;\ntranslate 1 'Homo sapiens', 2 Pan, 3 'it''s';\n"
            "tree T = ((1:1,2:1):1,3:2);\nend;\n")
    names = sorted(n.name for n in read_tree_text(text).nodes if n.is_tip)
    assert names == ["Homo sapiens", "Pan", "it's"]


def test_multi_tree_and_duplicate_tips_are_reported():
    from repeatviewer.backends.tree_layout import read_tree_text

    notes = read_tree_text("(A:1,A:1,B:1);\n(A:1,B:1);\n").notes
    assert any("2 trees" in n for n in notes)
    assert any("more than once" in n for n in notes)
    assert read_tree_text("(A:1,B:1);").notes == []
