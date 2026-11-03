"""Tree tip → sequence: finding a tree's alignment and the alignment window."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QGraphicsSimpleTextItem, QMenu

from repeatviewer.backends import tree_links as tk

IQTREE = "((A:0.1,B:0.2)92/97:0.05,(C:0.3,(D:0.1,E:0.1)100:0.2)88/99:0.1,F:0.4);\n"
ALN = ">A\nACGTACGT--\n>B\nACGTACGA--\n>C\nACGTTCGTAA\n>D\n--GTACGTAA\n>E\nACGTACGTAA\n>F\nTTTTACGTAA\n"


def _run(tmp_path: Path) -> Path:
    tree = tmp_path / "sat.fa_tree.treefile"
    tree.write_text(IQTREE)
    (tmp_path / "sat.fa_tree.aln.fa").write_text(ALN)
    details = tmp_path / "sat.fa_tree_details"
    details.mkdir()
    (details / "sat.fa_tree.names.tsv").write_text(
        "tree_name\toriginal_header\tsource_file\n"
        + "".join(f"{n}\t{n} copy\tsat.fa\n" for n in "ABCDEF"))
    (tmp_path / "sat.fa").write_text(ALN.replace("-", ""))
    return tree


def test_run_tree_and_base():
    assert tk.run_tree(Path("/x/a.fa_tree.treefile_rerooted.treefile")) == Path("/x/a.fa_tree.treefile")
    assert tk.run_base(Path("a.fa_tree_2_masked.treefile")) == ("a.fa_tree_2", True)
    assert tk.run_base(Path("a.nwk")) is None


def test_sources_in_order(tmp_path):
    tree = _run(tmp_path)
    src = tk.sequences_for_tree(tree)
    assert src.path.name == "sat.fa_tree.aln.fa" and src.aligned
    saved = tmp_path / "sat.fa_tree.treefile_rerooted.treefile"
    saved.write_text(IQTREE)
    assert tk.sequences_for_tree(saved).path == src.path
    masked = tmp_path / "sat.fa_tree_masked.treefile"
    masked.write_text(IQTREE)
    (tmp_path / "sat.fa_tree_details" / "sat.fa_tree_masked.aln.fa").write_text(ALN)
    assert tk.sequences_for_tree(masked).path.parent.name == "sat.fa_tree_details"
    # no alignment left: the input FASTA from names.tsv, marked as not aligned
    (tmp_path / "sat.fa_tree.aln.fa").unlink()
    fallback = tk.sequences_for_tree(tree)
    assert fallback.path.name == "sat.fa" and not fallback.aligned and fallback.note


def test_iqtree_and_sibling_naming(tmp_path):
    (tmp_path / "mono.aln.fasta").write_text(ALN)
    t = tmp_path / "mono.aln.fasta.treefile"
    t.write_text(IQTREE)
    assert tk.sequences_for_tree(t).path.name == "mono.aln.fasta"
    (tmp_path / "other.fa").write_text(ALN)
    n = tmp_path / "other.nwk"
    n.write_text(IQTREE)
    assert tk.sequences_for_tree(n).path.name == "other.fa"
    lone = tmp_path / "lone.nwk"
    lone.write_text(IQTREE)
    assert tk.sequences_for_tree(lone) is None


def test_rows_consensus_and_unaligned(tmp_path):
    f = tmp_path / "x.fa"
    f.write_text(">chr1:10-20(+) sat\nACG\n>chr1:10-20(+) again\nACT\n>Seq B\nAC-\n")
    aln = tk.read_alignment(f)
    assert aln.aligned and aln.width == 3
    assert aln.row_of("chr1_10-20") == 0          # first word, made tree-safe
    assert aln.row_of("chr1_10-20_2") == 1        # duplicate numbered like the tree builder
    assert aln.row_of("Seq") == 2 and aln.row_of("nope") is None
    assert tk.consensus(aln.seqs) == "ACG"
    assert tk.first_residue("--AC") == 2
    g = tmp_path / "y.fa"
    g.write_text(">a\nACGT\n>b\nAC\n")
    assert not tk.read_alignment(g).aligned


def test_double_click_tip_opens_alignment_highlighted(qapp, tmp_path):
    from repeatviewer.widgets.tree_view import TreePanel

    panel = TreePanel(_run(tmp_path))
    canvas = panel.canvas()
    canvas.resize(800, 600)
    item = next(i for i in panel.scene().items()
                if isinstance(i, QGraphicsSimpleTextItem) and i.text() == "D")
    idx = canvas.node_at(canvas.mapFromScene(item.sceneBoundingRect().center()))
    panel._on_node_double_clicked(idx)
    win = panel._aln_window
    assert win is not None and win.path.name == "sat.fa_tree.aln.fa"
    qapp.processEvents()
    assert win.canvas.highlighted == 3
    assert "D:" in win._status.text()
    # another tip reuses the window
    assert panel.show_in_alignment("F") and panel._aln_window is win
    assert win.canvas.highlighted == 5
    win.canvas.set_differences_only(True)
    win.canvas.grab()                                  # paints without error
    # right-click menu offers it for a tip
    menu = QMenu()
    panel._add_clade_actions(menu, ("F",))
    assert {a.text(): a.isEnabled() for a in menu.actions()}["Show in alignment"]
    win.close()


def test_double_click_empty_space_still_fits(qapp, tmp_path):
    from repeatviewer.widgets.tree_view import TreePanel

    panel = TreePanel(_run(tmp_path))
    panel._on_node_double_clicked(-1)
    assert panel._aln_window is None


def test_distance_matrices():
    from repeatviewer.backends import tree_layout as tl

    m = tl.read_tree_text("((A:0.1,B:0.2):0.05,C:0.3);")
    names, d = tk.patristic_matrix(m)
    ia, ib, ic = names.index("A"), names.index("B"), names.index("C")
    assert abs(d[ia, ib] - 0.3) < 1e-9 and abs(d[ia, ic] - 0.45) < 1e-9 and d[ia, ia] == 0
    aln = tk.Alignment(["x", "y", "z"], ["ACGT", "ACGA", "--NN"], True)
    names, p = tk.p_distance_matrix(aln)
    assert abs(p[0, 1] - 0.25) < 1e-9
    import math
    assert math.isnan(p[0, 2])
    tsv = tk.matrix_tsv(names, p)
    assert tsv.splitlines()[0] == "\tx\ty\tz" and "NA" in tsv


def test_distance_dialog_saves_next_to_tree(qapp, tmp_path):
    from repeatviewer.backends import tree_layout as tl
    from repeatviewer.views.tree_distances_dialog import SOURCE_ALIGNMENT, TreeDistancesDialog

    tree = _run(tmp_path)
    dlg = TreeDistancesDialog(tree, tl.read_tree_file(tree))
    names, matrix = dlg.matrix()
    assert sorted(names) == list("ABCDEF") and matrix.shape == (6, 6)
    out = dlg.save_to(dlg.default_path())
    assert out.name == "sat.fa_tree.treefile_distances_patristic.tsv" and out.is_file()
    dlg._source_box.setCurrentIndex(dlg._source_box.findData(SOURCE_ALIGNMENT))
    assert dlg.default_path().name.endswith("_pdistance.tsv")
    assert dlg.matrix()[1].shape == (6, 6)


def test_alignment_in_details_folder_is_found_first(tmp_path):
    tree = _run(tmp_path)
    (tmp_path / "sat.fa_tree.aln.fa").rename(tmp_path / "sat.fa_tree_details" / "sat.fa_tree.aln.fa")
    src = tk.sequences_for_tree(tree)
    assert src.path == tmp_path / "sat.fa_tree_details" / "sat.fa_tree.aln.fa" and src.aligned


def test_export_button_and_report_button(qapp, tmp_path):
    from repeatviewer.widgets.tree_view import TreePanel

    tree = _run(tmp_path)
    panel = TreePanel(tree)
    # Export opens a Type+File dialog (no save submenu). Distance / Newick live under Edit.
    assert panel._save_btn.text() == "Export"
    assert panel._save_btn.menu() is None
    panel._fill_edit_menu()
    edit_texts = [a.text() for a in panel._edit_menu.actions() if a.text()]
    assert "Distance matrix…" in edit_texts
    assert "Copy Newick" in edit_texts
    assert "Save tree (.treefile)" in edit_texts
    assert panel.report_path() is None and not panel._report_btn.isEnabled()
    report = tmp_path / "sat.fa_tree_details" / "sat.fa_tree_report.html"
    report.write_text("<html></html>")
    panel._update_status()
    assert panel.report_path() == report and panel._report_btn.isEnabled()


def test_clade_sequences_original_and_aligned(tmp_path):
    tree = _run(tmp_path)
    (tmp_path / "sat.fa").write_text(">A copy\nACGTACGT\n>B copy\nACGTACGA\n>C copy\nACGTTCGTAA\n"
                                     ">D copy\nGTACGTAA\n>E copy\nACGTACGTAA\n>F copy\nTTTTACGTAA\n")
    got = tk.clade_sequences(tree, ["D", "E"], aligned=False)
    assert got.records == [("D copy", "GTACGTAA"), ("E copy", "ACGTACGTAA")] and not got.missing
    got = tk.clade_sequences(tree, ["D", "nope"], aligned=True)
    assert got.records == [("D", "--GTACGTAA")] and got.missing == ["nope"]
    # a tree from elsewhere: the alignment without gaps, and said so
    other = tmp_path / "other.nwk"
    other.write_text(IQTREE)
    (tmp_path / "other.fa").write_text(ALN)
    got = tk.clade_sequences(other, ["A"], aligned=False)
    assert got.records == [("A", "ACGTACGT")] and "gaps removed" in got.note
    assert tk.extract_default_path(tree, 2, aligned=False).name == "sat.fa_tree.treefile_2seqs.fa"
    (tmp_path / "sat.fa_tree.treefile_2seqs.fa").write_text("")
    assert tk.extract_default_path(tree, 2, aligned=False).name == "sat.fa_tree.treefile_2seqs_2.fa"
    assert tk.extract_default_path(tree, 1, aligned=True).name == "sat.fa_tree.treefile_1seq_aligned.fa"


def test_extract_from_the_clade_menu(qapp, tmp_path):
    from repeatviewer.widgets.tree_view import TreePanel

    tree = _run(tmp_path)
    panel = TreePanel(tree)
    menu = QMenu()
    panel._add_clade_actions(menu, ("C", "D", "E"))
    sub = next(a.menu() for a in menu.actions() if a.menu() is not None)
    assert sub.title() == "Extract sequences (3)"
    assert [a.text() for a in sub.actions()] == ["Original sequences…", "As aligned…"]
    written = []
    panel.fileWritten.connect(written.append)
    out = panel.extract_sequences(("D", "E"), aligned=True, path=tmp_path / "de.fa")
    assert out.read_text().startswith(">D\n--GTACGTAA\n>E\n") and written == [out]
