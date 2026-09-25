"""Tree tab (Preview panel): loads Newick files, redraws, exports SVG/PDF/PNG."""

from __future__ import annotations

from pathlib import Path

import pytest
from PySide6.QtWidgets import QGraphicsSimpleTextItem

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
