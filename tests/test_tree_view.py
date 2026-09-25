"""Tree tab (Preview panel): loads Newick files, redraws, exports SVG/PDF/PNG."""

from __future__ import annotations

from pathlib import Path

import pytest
from PySide6.QtWidgets import QGraphicsSimpleTextItem

from repeatviewer.backends import tree_edit as te
from repeatviewer.backends import tree_layout as tl

IQTREE = "((A:0.1,B:0.2)92/97:0.05,(C:0.3,(D:0.1,E:0.1)100:0.2)88/99:0.1,F:0.4);\n"


@pytest.fixture
def tree_file(tmp_path: Path) -> Path:
    path = tmp_path / "sat.fa_tree.treefile"
    path.write_text(IQTREE, encoding="utf-8")
    return path


def _texts(panel) -> list[str]:
    return sorted(i.text() for i in panel.scene().items() if isinstance(i, QGraphicsSimpleTextItem))


def test_loads_tree_with_tip_names_and_support(qapp, tree_file):
    from repeatviewer.widgets.tree_view import TreePanel, TreeViewOptions

    panel = TreePanel(tree_file)
    panel.set_options(TreeViewOptions())
    texts = _texts(panel)
    for name in "ABCDEF":
        assert name in texts
    assert "92/97" in texts and "88/99" in texts
    assert "0.05" in texts         # scale bar: 1/5 of the 0.4 extent, rounded down


@pytest.mark.parametrize("layout", tl.LAYOUTS)
def test_every_layout_draws(qapp, tree_file, layout):
    from repeatviewer.widgets.tree_view import TreePanel, TreeViewOptions

    panel = TreePanel(tree_file)
    panel.set_options(TreeViewOptions(layout=layout, branch_labels="length"))
    assert not panel.scene().sceneRect().isEmpty()
    assert "0.4" in _texts(panel)  # F's branch length


def test_tip_names_and_labels_can_be_hidden(qapp, tree_file):
    from repeatviewer.widgets.tree_view import TreePanel, TreeViewOptions

    panel = TreePanel(tree_file)
    panel.set_options(TreeViewOptions(tip_labels=False, branch_labels="none", scale_bar=False))
    assert _texts(panel) == []


def test_tree_without_support_starts_without_support_labels(qapp, tmp_path):
    from repeatviewer.widgets.tree_view import BRANCH_LABELS_NONE, TreePanel

    path = tmp_path / "plain.nwk"
    path.write_text("((A:1,B:1):1,(C:1,D:1):1);", encoding="utf-8")
    panel = TreePanel(path)
    assert panel.options().branch_labels == BRANCH_LABELS_NONE


def test_unreadable_file_shows_a_message(qapp, tmp_path):
    from repeatviewer.widgets.tree_view import TreePanel

    path = tmp_path / "broken.tre"
    path.write_text("this is not newick", encoding="utf-8")
    panel = TreePanel()
    assert panel.load_tree(path) is False
    assert panel.model() is None


@pytest.mark.parametrize("ext,magic", [("svg", b"<?xml"), ("pdf", b"%PDF"), ("png", b"\x89PNG")])
def test_export_formats(qapp, tree_file, ext, magic):
    from repeatviewer.widgets.tree_view import TreeViewOptions, export_default_path, export_tree

    out = export_default_path(tree_file, "circular", ext)
    assert out.name == f"sat.fa_tree.treefile_circular.{ext}"   # rule 1 naming
    model = tl.read_tree_file(tree_file)
    export_tree(model, TreeViewOptions(layout="circular"), out)
    assert out.read_bytes()[: len(magic)] == magic


def test_preview_panel_opens_trees_in_the_tree_view(qapp, tree_file):
    from repeatviewer.widgets.preview_tabs import FilePreviewPage, is_previewable_file

    assert is_previewable_file(tree_file)
    page = FilePreviewPage(tree_file)
    assert page.is_tree_preview_active()
    assert page._tree_panel.model().n_tips == 6


# --------------------------------------------------------------------------
# A2: editing in the viewer
# --------------------------------------------------------------------------


def _panel(qapp, tree_file):
    from repeatviewer.widgets.tree_view import TreePanel, TreeViewOptions

    panel = TreePanel(tree_file)
    panel.set_options(TreeViewOptions())
    return panel


def _tip_order(panel):
    d = panel.display().model
    return [d.nodes[t].name for t in d.tips()]


def _branch_items(panel):
    from repeatviewer.widgets.tree_view import _BranchItem

    return [i for i in panel.scene().items() if isinstance(i, _BranchItem)]


def test_clicking_a_branch_selects_its_clade(qapp, tree_file):
    panel = _panel(qapp, tree_file)
    canvas = panel.canvas()
    canvas.resize(800, 600)
    item = next(i for i in panel.scene().items()
                if isinstance(i, QGraphicsSimpleTextItem) and i.text() == "D")
    view_pos = canvas.mapFromScene(item.sceneBoundingRect().center())
    idx = canvas.node_at(view_pos)
    assert panel.display().keys[idx] == ("D",)
    panel._on_node_clicked(idx)
    assert panel._selected == ("D",)
    panel._on_node_clicked(idx)            # clicking again clears
    assert panel._selected is None


def test_reroot_midpoint_and_undo(qapp, tree_file):
    panel = _panel(qapp, tree_file)
    before = te.to_newick(panel.model())
    panel.reroot_at(("F",))
    root_kids = [panel.model().nodes[k].name for k in panel.model().root.children]
    assert "F" in root_kids
    assert panel.state().edited == "rerooted"
    panel.midpoint_root()
    panel.undo()
    panel.undo()
    assert te.to_newick(panel.model()) == before
    assert not panel._undo


def test_swap_collapse_expand(qapp, tree_file):
    panel = _panel(qapp, tree_file)
    order = _tip_order(panel)
    panel.swap_at(("A", "B"))
    assert _tip_order(panel)[:2] == ["B", "A"] and order[:2] == ["A", "B"]
    panel.toggle_collapse(("C", "D", "E"))
    d = panel.display()
    assert d.model.n_tips == 4 and len(d.info) == 1
    assert "C + 2 more" in _texts(panel)
    panel.toggle_collapse(("C", "D", "E"))
    assert panel.display().model.n_tips == 6


def test_colour_and_show_only_this_clade(qapp, tree_file):
    panel = _panel(qapp, tree_file)
    panel.set_clade_color(("C", "D", "E"), "#d62728")
    red = [i for i in _branch_items(panel) if i.pen().color().name() == "#d62728"]
    assert len(red) >= 5                   # the clade's branches, not the rest
    panel.show_only(("C", "D", "E"))
    assert sorted(_tip_order(panel)) == ["C", "D", "E"]
    panel.full_tree()
    assert panel.display().model.n_tips == 6


def test_collapse_below_support(qapp, tree_file):
    panel = _panel(qapp, tree_file)
    panel.set_collapse_below(98)           # 92/97 goes (UFBoot 97)
    assert panel.display().n_contracted == 1
    panel.set_collapse_below(0)
    assert panel.display().n_contracted == 0


def test_view_is_remembered_next_to_the_tree(qapp, tree_file):
    from repeatviewer.widgets.tree_view import TreePanel, view_state_path

    panel = _panel(qapp, tree_file)
    panel.set_clade_color(("A", "B"), "#1f77b4")
    panel.toggle_collapse(("C", "D", "E"))
    sidecar = view_state_path(tree_file)
    assert sidecar.name == ".sat.fa_tree.treefile.view.json"   # hidden: not in the file list
    again = TreePanel(tree_file)
    assert again.state().colors == {("A", "B"): "#1f77b4"}
    assert again.state().collapsed == {("C", "D", "E")}


def test_save_tree_as_writes_a_new_file(qapp, tree_file):
    from repeatviewer.widgets.tree_view import saved_tree_default_path

    panel = _panel(qapp, tree_file)
    original = tree_file.read_text()
    panel.reroot_at(("F",))
    out = saved_tree_default_path(tree_file, panel.state().edited)
    assert out.name == "sat.fa_tree.treefile_rerooted.treefile"   # rule 1
    panel.save_tree(out)
    saved = tl.read_tree_text(out.read_text())
    assert "F" in [saved.nodes[k].name for k in saved.root.children]
    assert tree_file.read_text() == original                       # never overwritten


def test_export_shows_collapsed_clade(qapp, tree_file, tmp_path):
    from repeatviewer.widgets.tree_view import export_tree

    panel = _panel(qapp, tree_file)
    panel.toggle_collapse(("C", "D", "E"))
    out = tmp_path / "x.svg"
    export_tree(panel.display(), panel.options(), out)
    svg = out.read_text()
    assert "C + 2 more" in svg and ">D<" not in svg


def test_context_menu_actions(qapp, tree_file):
    from PySide6.QtWidgets import QMenu

    panel = _panel(qapp, tree_file)
    menu = QMenu()
    panel._add_clade_actions(menu, ("C", "D", "E"))
    enabled = {a.text(): a.isEnabled() for a in menu.actions() if a.text()}
    assert enabled["Reroot on this branch"] and enabled["Swap branches"]
    assert enabled["Collapse clade"] and enabled["Show only this clade"]
    tip_menu = QMenu()
    panel._add_clade_actions(tip_menu, ("F",))
    tip = {a.text(): a.isEnabled() for a in tip_menu.actions() if a.text()}
    assert tip["Reroot on this “F”"] and not tip["Swap branches"] and not tip["Collapse clade"]


def test_reopening_the_same_unchanged_tree_keeps_unsaved_edits(qapp, tree_file):
    import os

    panel = _panel(qapp, tree_file)
    panel.reroot_at(("F",))
    panel.load_tree(tree_file)             # e.g. the file list refreshed
    assert panel.state().edited == "rerooted" and panel._undo
    st = tree_file.stat()
    os.utime(tree_file, (st.st_atime, st.st_mtime + 5))   # the file changed on disk
    panel.load_tree(tree_file)
    assert panel.state().edited == "" and not panel._undo


# --------------------------------------------------------------------------
# A3: trait table colours
# --------------------------------------------------------------------------


@pytest.fixture
def located_tree(tmp_path: Path) -> Path:
    """A tree run with copies from two arrays on Chr1 and one on Chr2."""
    tips = ["Chr1_100-256_+", "Chr1_260-416_-", "Chr1_420-576_+",
            "Chr1_90000-90156_+", "Chr2_5000-5156_-", "Chr2_5160-5316_+"]
    tree = tmp_path / "sat.fa_tree.treefile"
    tree.write_text(f"((({tips[0]}:0.1,{tips[1]}:0.1)95/99:0.05,{tips[2]}:0.2)90/97:0.05,"
                    f"({tips[3]}:0.2,({tips[4]}:0.1,{tips[5]}:0.1)99/100:0.1)80/96:0.05);\n")
    details = tmp_path / "sat.fa_tree_details"
    details.mkdir()
    lines = ["tree_name\toriginal_header\tsource_file"]
    for t in tips:
        chrom, rest = t.split("_", 1)
        span, strand = rest.rsplit("_", 1)
        lines.append(f"{t}\t{chrom}:{span}({strand})\tsat.fa")
    (details / "sat.fa_tree.names.tsv").write_text("\n".join(lines) + "\n")
    return tree


def test_trait_table_colours_tips_by_array(qapp, located_tree):
    from repeatviewer.backends import tree_traits as tt

    panel = _panel(qapp, located_tree)
    assert not panel._color_box.isEnabled()           # no table yet
    panel.build_traits(wait=True)
    assert tt.traits_path(located_tree).is_file()
    assert panel.state().color_by == "array"          # sensible default
    colours = set(panel.display().tip_colors.values())
    assert len(colours) == 3                          # two arrays on Chr1, one on Chr2
    texts = _texts(panel)
    assert "array" in texts and any(t.startswith("Chr2:") for t in texts)   # legend
    notes = " ".join(panel.display().tip_notes.values())
    assert "chromosome: Chr1" in notes


def test_colour_by_methylation_gives_a_colour_bar(qapp, located_tree, tmp_path):
    cov = tmp_path / "leaf_bismark.cov"
    cov.write_text("Chr1\t150\t150\t100\t10\t0\nChr2\t5100\t5100\t0\t0\t10\n")
    panel = _panel(qapp, located_tree)
    panel._state.trait_settings = {"methylation_file": str(cov), "array_gap": 10_000}
    panel.build_traits(wait=True)
    assert panel.state().color_by == "methylation"
    legend = panel.display().legend
    assert legend.numeric and (legend.low, legend.high) == (0.0, 100.0)
    assert legend.n_missing == 4                      # only two copies have sites
    assert "no value: 4" in _texts(panel)


def test_colour_choice_is_remembered_and_exported(qapp, located_tree, tmp_path):
    from repeatviewer.widgets.tree_view import TreePanel, export_tree

    panel = _panel(qapp, located_tree)
    panel.build_traits(wait=True)
    panel.set_color_by("chromosome")
    again = TreePanel(located_tree)
    assert again.state().color_by == "chromosome"
    out = tmp_path / "t.svg"
    export_tree(again.display(), again.options(), out)
    svg = out.read_text()
    assert ">chromosome<" in svg and ">Chr1<" in svg


def test_user_columns_appear_after_reload(qapp, located_tree):
    from repeatviewer.backends import tree_traits as tt

    panel = _panel(qapp, located_tree)
    panel.build_traits(wait=True)
    table = tt.read_trait_table(tt.traits_path(located_tree))
    table.columns.append("karyotype")
    for i, row in enumerate(table.rows.values()):
        row["karyotype"] = "B" if i % 2 else "A"
    tt.write_trait_table(tt.traits_path(located_tree), table)
    panel._reload_traits()
    panel.set_color_by("karyotype")
    assert set(panel.display().tip_colors.values()) == {tt.CATEGORICAL[0], tt.CATEGORICAL[1]}


def test_trait_dialog_rules(qapp, located_tree, tmp_path):
    from PySide6.QtWidgets import QPushButton

    from repeatviewer.views.tree_traits_dialog import TraitTableDialog

    dlg = TraitTableDialog(tree_path=located_tree)
    names = [b.text() for b in dlg.findChildren(QPushButton) if b.text() in ("Cancel", "Save", "Run")]
    assert names == ["Cancel", "Save", "Run"]                         # rule 2
    cov = tmp_path / "x.cov"
    cov.write_text("Chr1\t1\t1\t0\t0\t1\n")
    fasta = tmp_path / "x.fa"
    fasta.write_text(">a\nACGT\n")
    assert dlg.apply_path_from_panel(cov) and not dlg.apply_path_from_panel(fasta)   # rule 7
    assert dlg.settings()["methylation_file"].endswith("x.cov")
    assert dlg.settings()["array_gap"] == 10_000
