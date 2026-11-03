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


def test_double_click_tip_asks_preview_to_select_it(qapp, tmp_path):
    from repeatviewer.widgets.tree_view import TreePanel

    panel = TreePanel(_run(tmp_path))
    asked = []
    panel.alignmentRequested.connect(lambda path, ids, folder: asked.append((Path(path).name, list(ids))))
    idx = next(n.index for n in panel._display.model.nodes if n.name == "D")
    panel._on_node_double_clicked(idx)
    assert asked == [("sat.fa_tree.aln.fa", ["D"])]
    # a clade selects all of its rows; unknown names are skipped
    assert panel.show_in_alignment(["C", "D", "nope"])
    assert asked[-1] == ("sat.fa_tree.aln.fa", ["C", "D"])
    menu = QMenu()
    panel._add_clade_actions(menu, ("F",))
    assert {a.text(): a.isEnabled() for a in menu.actions()}["Show in alignment"]
    menu = QMenu()
    panel._add_clade_actions(menu, ("C", "D", "E"))
    assert {a.text(): a.isEnabled() for a in menu.actions()}["Show clade in alignment (3)"]


def test_preview_opens_alignment_with_rows_selected(qapp, tmp_path):
    from repeatviewer.widgets.preview_tabs import PreviewTabWidget

    tree = _run(tmp_path)
    tabs = PreviewTabWidget()
    tabs.resize(1000, 700)
    tabs.open_file(tree, replace=False)
    tree_page = tabs.currentWidget()
    tree_page.tree_panel().show_in_alignment(["C", "E"])
    page = tabs.currentWidget()
    assert page is not tree_page and page.path.name == "sat.fa_tree.aln.fa"
    overview = page._annotation_overview
    assert [overview._record_ids[r] for r in overview._bars.selected_rows()] == ["C", "E"]
    # asking again reuses that tab
    tree_page.tree_panel().show_in_alignment("A")
    assert tabs.currentWidget() is page and tabs.count() == 2
    assert [overview._record_ids[r] for r in overview._bars.selected_rows()] == ["A"]


def test_double_click_empty_space_still_fits(qapp, tmp_path):
    from repeatviewer.widgets.tree_view import TreePanel

    panel = TreePanel(_run(tmp_path))
    asked = []
    panel.alignmentRequested.connect(lambda *a: asked.append(a))
    panel._on_node_double_clicked(-1)
    assert asked == []


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
    assert [a.text() for a in sub.actions()] == ["Without gaps…", "As aligned…"]
    written = []
    panel.fileWritten.connect(written.append)
    out = panel.extract_sequences(("D", "E"), aligned=True, path=tmp_path / "de.fa")
    assert out.read_text().startswith(">D\n--GTACGTAA\n>E\n") and written == [out]


def test_without_gaps_strips_gaps_when_the_input_was_an_alignment(tmp_path):
    tree = _run(tmp_path)
    (tmp_path / "sat.fa").write_text(ALN)            # the tree was built from an MSA
    got = tk.clade_sequences(tree, ["A", "D"], aligned=False)
    assert got.records == [("A copy", "ACGTACGT"), ("D copy", "GTACGTAA")]
    aligned = tk.clade_sequences(tree, ["D"], aligned=True)
    assert aligned.records == [("D", "--GTACGTAA")]


def test_extract_from_tree_alignment_saves_next_to_tree_without_gaps(qapp, tmp_path):
    from repeatviewer.widgets.preview_tabs import PreviewTabWidget

    tree = _run(tmp_path)
    details = tmp_path / "sat.fa_tree_details"
    aln = details / "sat.fa_tree.aln.fa"            # where tree runs keep it
    (tmp_path / "sat.fa_tree.aln.fa").rename(aln)
    aln.write_text(ALN.replace(">D\n", ">D copy 7 chr2:10-20\n"))
    tabs = PreviewTabWidget()
    tabs.resize(1000, 700)
    tabs.open_file(tree, replace=False)
    written = []
    tabs.sequenceExtracted.connect(written.append)
    tabs.currentWidget().tree_panel().show_in_alignment("D")
    overview = tabs.currentWidget()._annotation_overview
    overview._extract_selected(without_gaps=True)
    overview._extract_selected(without_gaps=False)
    assert [p.parent for p in written] == [tmp_path, tmp_path]        # next to the tree
    nogap, gapped = written
    assert nogap.name.endswith("_D_nogap.fa") and gapped.name.endswith("_D.fa")
    assert nogap.read_text() == ">D copy 7 chr2:10-20\nGTACGTAA\n"    # full header, no gaps
    assert gapped.read_text().splitlines() == [">D copy 7 chr2:10-20", "--GTACGTAA"]


def test_self_dotplot_opens_in_a_new_tab_and_details_extracts_go_next_to_tree(qapp, tmp_path):
    from repeatviewer.widgets.preview_tabs import PreviewTabWidget

    _run(tmp_path)
    aln = tmp_path / "sat.fa_tree_details" / "sat.fa_tree.aln.fa"
    (tmp_path / "sat.fa_tree.aln.fa").rename(aln)
    tabs = PreviewTabWidget()
    tabs.resize(1000, 700)
    tabs.open_file(aln, replace=False)                 # opened from Files, not the tree
    page = tabs.currentWidget()
    page.select_records(["D"])
    overview = page._annotation_overview
    assert overview._output_folder() == tmp_path       # next to the tree anyway
    overview.showSelfDotplotRequested.emit("D")
    assert tabs.count() == 2 and tabs.tabText(1) == "Dotplot: D"
    assert tabs.widget(0) is page and tabs.currentIndex() == 1
    for _ in range(5):
        qapp.processEvents()
    assert tabs.widget(1)._fasta_stack.currentIndex() == 2          # Dotplot (Self), not Overview
    tabs.open_file(aln, replace=False)                 # the file's own tab is still the alignment
    assert tabs.currentWidget() is page and tabs.count() == 2


def test_deleting_sequences_asks_first_and_names_the_tree(qapp, tmp_path, monkeypatch):
    from repeatviewer.widgets import multi_sequence_overview as mso
    from repeatviewer.widgets.preview_tabs import PreviewTabWidget

    tree = _run(tmp_path)
    aln = tmp_path / "sat.fa_tree_details" / "sat.fa_tree.aln.fa"
    (tmp_path / "sat.fa_tree.aln.fa").rename(aln)
    assert mso.tree_built_from(aln) == tree
    assert mso.tree_built_from(tmp_path / "sat.fa") is None
    asked = []
    answer = [False]
    monkeypatch.setattr(mso, "confirm_delete_sequences",
                        lambda parent, n, path, t: asked.append((n, t)) or answer[0])
    tabs = PreviewTabWidget()
    tabs.open_file(aln, replace=False)
    page = tabs.currentWidget()
    page.select_records(["B", "C"])
    overview = page._annotation_overview
    overview._delete_selected_msa_rows(hard_delete=True)          # Cancel
    assert asked == [(2, tree)] and aln.read_text().count(">") == 6
    answer[0] = True
    overview._delete_selected_msa_rows(hard_delete=True)          # Delete
    assert aln.read_text().count(">") == 4
