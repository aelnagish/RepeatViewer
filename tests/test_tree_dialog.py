"""Tree dialog: shared rows above the tabs, IQ-TREE / FastME / FastTree / Options tabs."""

from __future__ import annotations

from PySide6.QtCore import QElapsedTimer
from PySide6.QtWidgets import QPushButton, QRadioButton, QToolButton

from repeatviewer.job_options import (
    TREE_ALIGN_ALWAYS,
    TREE_ALIGN_AUTO,
    TREE_ALIGN_NEVER,
    TREE_ENGINE_FASTME,
    TREE_ENGINE_FASTTREE,
    TREE_ENGINE_IQTREE,
    TREE_ROOT_MIDPOINT,
    TREE_ROOT_OUTGROUP,
    TREE_SEQ_AS_IS,
    TreeOptions,
)


def _dialog(qapp, tmp_path, **kw):
    from repeatviewer.views.tree_dialog import TreeParamsDialog

    fa = tmp_path / "sat.fa"
    fa.write_text(">m1\nACGTACGTAC\n>m2\nACGTACGTAA\n>m3\nACGAACGTAA\n>m4\nACGAAC\n")
    return TreeParamsDialog(files=[fa], **kw)


def _wait_summary(dlg, qapp, timeout_ms: int = 5000) -> None:
    """Sequence lengths load in a background thread after the dialog opens."""
    if not dlg._files:
        return
    timer = QElapsedTimer()
    timer.start()
    while dlg._summary is None and timer.elapsed() < timeout_ms:
        qapp.processEvents()
    assert dlg._summary is not None, "input summary did not load"
    dlg._sync()


def _shown(dlg, widget) -> bool:
    return widget.isVisibleTo(dlg)


def _normalised(opts: TreeOptions) -> TreeOptions:
    """What the window writes back: MAFFT is chosen automatically, per-program PATH boxes."""
    from dataclasses import replace

    from repeatviewer.backends import tree_builder as tb

    legacy = opts.use_system_path
    return replace(opts, align_mode=TREE_ALIGN_ALWAYS if opts.align_mode == TREE_ALIGN_AUTO
                   else opts.align_mode, mafft_strategy=tb.MAFFT_CHOOSE, adjust_direction=True,
                   mafft_gap_open=1.53, mafft_offset=0.0, bnni=True,
                   right_click_action="run",
                   sequence_mode=TREE_SEQ_AS_IS,
                   use_system_iqtree=opts.use_system_iqtree or legacy,
                   use_system_fastme=opts.use_system_fastme or legacy,
                   use_system_fasttree=opts.use_system_fasttree or legacy,
                   use_system_path=False)


def _tab_texts(dlg):
    return [dlg._tabs.tabText(i) for i in range(dlg._tabs.count())]


def test_default_options_round_trip(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path, initial=TreeOptions(threads=2))
    _wait_summary(dlg, qapp)
    assert dlg.options() == _normalised(TreeOptions(threads=2))
    assert dlg._support_iq.currentText() == "SH-aLRT + UFBoot, 1000 each"   # not a custom entry
    assert dlg._support_iq.count() == 4


def test_custom_options_round_trip_keeps_custom_support(qapp, tmp_path):
    opts = TreeOptions(
        engine=TREE_ENGINE_FASTME, model="GTR+F+G4", runs=3, seed=7, threads=3,
        min_coverage=0.5, root_mode=TREE_ROOT_OUTGROUP, outgroup=["m2"], alrt=0, ufboot=2000,
        distance_bootstrap=250, right_click_engine=TREE_ENGINE_FASTTREE,
        show_process_window=True, use_system_mafft=True, use_system_fastme=True,
        coverage_reference="longest", align_mode=TREE_ALIGN_NEVER,
    )
    dlg = _dialog(qapp, tmp_path, initial=opts)
    _wait_summary(dlg, qapp)
    assert dlg.options() == _normalised(opts)
    assert dlg._support_iq.currentText() == "SH-aLRT 0, UFBoot 2000"
    assert dlg._support_fme.currentText() == "Bootstrap, 250 replicates"
    assert dlg._tabs.currentIndex() == 1                     # opens on the saved engine's tab


def test_buttons_and_tabs(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    labels = [b.text() for b in dlg._button_box.findChildren(QPushButton)]
    assert labels == ["Restore Defaults", "Cancel", "Save", "Run"] and "OK" not in labels
    assert _tab_texts(dlg) == ["IQ-TREE", "FastME", "FastTree", "Options"]
    assert not hasattr(dlg, "_more_toggle")
    assert dlg._show_process_window is not None
    assert set(dlg._system_boxes) == {"iqtree", "fastme", "fasttree", "mafft"}
    assert set(dlg._found) == {"iqtree", "fastme", "fasttree", "mafft"}
    assert not hasattr(dlg, "_download_btn")
    assert not hasattr(dlg, "_seq_mode")


def test_shared_rows_always_visible_above_tabs(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    _wait_summary(dlg, qapp)
    shared = (dlg._align_mode, dlg._rotate, dlg._root_mode, dlg._min_cov, dlg._cov_ref)
    for i in range(dlg._tabs.count()):
        dlg._tabs.setCurrentIndex(i)
        for w in shared:
            assert _shown(dlg, w)
    assert dlg._shared.parentWidget() is dlg and dlg._shared is not dlg._tabs.currentWidget()
    assert "1 of 4 sequences" in dlg._frag_hint.text()    # m4 is 6 bp vs median 10 bp
    dlg._min_cov.setValue(0)
    assert "None of the 4" in dlg._frag_hint.text()


def test_engine_specific_options_sit_on_their_tabs(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    assert dlg._model.parentWidget() is not None
    for w in (dlg._support_iq, dlg._model, dlg._mask, dlg._runs, dlg._seed, dlg._threads):
        assert dlg._iq_tab.isAncestorOf(w)
    assert dlg._fme_tab.isAncestorOf(dlg._support_fme) and dlg._fme_tab.isAncestorOf(dlg._threads_fme)
    assert dlg._ft_tab.isAncestorOf(dlg._support_ft) and dlg._ft_tab.isAncestorOf(dlg._threads_ft)


def test_engine_comes_from_the_open_tab(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    assert dlg.options().engine == TREE_ENGINE_IQTREE
    dlg._tabs.setCurrentIndex(1)
    assert dlg.options().engine == TREE_ENGINE_FASTME
    dlg._tabs.setCurrentIndex(2)
    assert dlg.options().engine == TREE_ENGINE_FASTTREE
    dlg._tabs.setCurrentIndex(3)                             # Options: the last engine tab
    assert dlg.options().engine == TREE_ENGINE_FASTTREE


def test_threads_are_shared_between_tabs(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    dlg._threads_fme.setValue(7)
    assert dlg._threads.value() == 7 and dlg._threads_ft.value() == 7
    dlg._tabs.setCurrentIndex(2)
    assert dlg.options().threads == 7


def test_right_click_engine_radios(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path, initial=TreeOptions(engine=TREE_ENGINE_FASTME))
    assert list(dlg._right_click_buttons) == [TREE_ENGINE_IQTREE, TREE_ENGINE_FASTME,
                                              TREE_ENGINE_FASTTREE]
    dlg._right_click_buttons[TREE_ENGINE_FASTTREE].setChecked(True)
    o = dlg.options()
    assert o.right_click_engine == TREE_ENGINE_FASTTREE and o.engine == TREE_ENGINE_FASTME


def test_per_tool_system_path_boxes(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    dlg._system_boxes["fastme"].setChecked(True)
    dlg._system_boxes["mafft"].setChecked(True)
    o = dlg.options()
    assert (o.use_system_iqtree, o.use_system_fastme, o.use_system_fasttree,
            o.use_system_mafft) == (False, True, False, True)
    legacy = _dialog(qapp, tmp_path, initial=TreeOptions(use_system_path=True))
    o = legacy.options()
    assert o.use_system_iqtree and o.use_system_fastme and o.use_system_fasttree
    assert not o.use_system_mafft and not o.use_system_path


def test_rows_that_do_not_apply_are_hidden(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    _wait_summary(dlg, qapp)
    assert _shown(dlg, dlg._mask)
    dlg._align_mode.setCurrentIndex(dlg._align_mode.findData(TREE_ALIGN_NEVER))
    assert not _shown(dlg, dlg._rotate)                   # rotation needs a new alignment


def test_support_choices_write_the_replicates(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    dlg._support_iq.setCurrentIndex(dlg._support_iq.findData("0,1000"))
    o = dlg.options()
    assert (o.alrt, o.ufboot, o.bnni) == (0, 1000, True)
    dlg._support_iq.setCurrentIndex(dlg._support_iq.findData("0,0"))
    assert (dlg.options().alrt, dlg.options().ufboot) == (0, 0)
    dlg._support_fme.setCurrentIndex(dlg._support_fme.findData(1000))
    assert dlg.options().distance_bootstrap == 1000


def test_root_and_outgroup(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    _wait_summary(dlg, qapp)
    assert not dlg._outgroup_btn.isEnabled()
    dlg._root_mode.setCurrentIndex(dlg._root_mode.findData(TREE_ROOT_OUTGROUP))
    assert dlg._outgroup_btn.isEnabled() and dlg._outgroup_btn.text() == "Outgroup…"
    dlg._outgroup = ["m3"]
    dlg._sync()
    assert dlg._outgroup_btn.text() == "Outgroup (1)…" and dlg._outgroup_btn.toolTip() == "m3"
    assert dlg.options().outgroup == ["m3"] and dlg.options().root_mode == TREE_ROOT_OUTGROUP
    dlg._root_mode.setCurrentIndex(dlg._root_mode.findData(TREE_ROOT_MIDPOINT))
    assert dlg.options().root_mode == TREE_ROOT_MIDPOINT


def test_outgroup_picker(qapp):
    from PySide6.QtCore import Qt

    from repeatviewer.views.tree_dialog import OutgroupDialog

    dlg = OutgroupDialog(["a", "b", "c"], ["b"])
    assert dlg.chosen() == ["b"]
    dlg._list.item(2).setCheckState(Qt.CheckState.Checked)
    dlg._filter.setText("c")
    assert dlg._list.item(0).isHidden() and not dlg._list.item(2).isHidden()
    assert dlg.chosen() == ["b", "c"]
    labels = [b.text() for b in dlg.findChildren(QPushButton)]
    assert labels == ["Cancel", "Save"]


def test_question_mark_help(qapp, tmp_path):
    from PySide6.QtWidgets import QLabel

    dlg = _dialog(qapp, tmp_path)
    radios = [r.text() for r in dlg.findChildren(QRadioButton)]
    assert [r.split(" (")[0] for r in radios] == ["IQ-TREE", "FastME", "FastTree"]
    assert not [lab for lab in dlg.findChildren(QLabel) if lab.text() == "ⓘ"]
    helps = [b for b in dlg.findChildren(QToolButton) if b.text() == "?"]
    # shared: Alignment, Monomers, Root, Skip fragments;
    # IQ-TREE: Branch support, Model, Methylation check, Searches;
    # FastME and FastTree: support  (no REXdb on Options)
    assert len(helps) == 10
    joined = " ".join(b.toolTip() for b in helps)
    for piece in ("L-INS-i", "FFT-NS-2", "Align tool", "ModelFinder", "UFBoot", "Midpoint", "CHG"):
        assert piece in joined


def test_protein_input_hides_dna_options(qapp, tmp_path):
    from repeatviewer.views.tree_dialog import TreeParamsDialog

    fa = tmp_path / "p.fa"
    fa.write_text("".join(f">p{i}\nMKVLLAQEEFWRSTLK\n" for i in range(4)))
    dlg = TreeParamsDialog(files=[fa], initial=TreeOptions(methylation_mask="CG"))
    _wait_summary(dlg, qapp)
    assert not _shown(dlg, dlg._rotate) and not _shown(dlg, dlg._mask)
    assert dlg.options().methylation_mask == ""


def test_restore_defaults(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path, initial=TreeOptions(engine=TREE_ENGINE_FASTME, runs=9,
                                                      rotate_monomers=True, methylation_mask="CG",
                                                      right_click_engine=TREE_ENGINE_FASTTREE))
    dlg._restore_defaults()
    o = dlg.options()
    assert (o.engine, o.runs, o.model, o.rotate_monomers, o.methylation_mask) == (
        "iqtree", 5, "MFP", False, "")
    assert o.right_click_engine == TREE_ENGINE_IQTREE


def test_opens_without_selection(qapp):
    from repeatviewer.views.tree_dialog import TreeParamsDialog

    dlg = TreeParamsDialog(files=[])
    assert not _shown(dlg, dlg._frag_hint)
    assert dlg.options().outgroup == []


def test_rexdb_not_on_options_tab(qapp, tmp_path):
    """REXdb is picked when running RT-domain tree; not shown on Tree Options."""
    rex = tmp_path / "Viridiplantae_v4.0.fasta"
    rex.write_text(">Ty1-RT__REXdb_ID1\nMKV\n")
    (tmp_path / "Viridiplantae_v4.0.classification").write_text("REXdb_ID1\tClass_I\tAle\n")
    opts = TreeOptions(threads=2, rexdb_fasta=str(rex), rexdb_references=5, rotate_monomers=True)
    dlg = _dialog(qapp, tmp_path, initial=opts)
    _wait_summary(dlg, qapp)
    assert dlg.options() == _normalised(opts)
    assert not hasattr(dlg, "_rexdb")
    assert not hasattr(dlg, "_rexdb_refs")
    other = tmp_path / "copies.fa"
    other.write_text(">a\nACGT\n")
    assert not dlg.apply_path_from_panel(other)
    assert dlg.apply_path_from_panel(rex)
    assert dlg.options().rexdb_fasta == str(rex)
    assert dlg.options().sequence_mode == TREE_SEQ_AS_IS
    assert dlg.width() == dlg._DIALOG_WIDTH


def test_options_tab_executable_status(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    dlg._tabs.setCurrentIndex(dlg._options_index)
    dlg._update_found()
    # Unchecked PATH boxes hide the status line.
    assert dlg._found["iqtree"].isHidden()
    assert dlg._found["fastme"].isHidden()
    dlg._system_boxes["iqtree"].setChecked(True)
    dlg._update_found()
    status = dlg._found["iqtree"].text()
    assert status
    assert status.startswith("/") or "built-in" in status
    assert "IQ-TREE" in dlg._system_boxes["iqtree"].toolTip()
    assert "bundled" in dlg._system_boxes["iqtree"].toolTip().lower() or "3.1.4" in dlg._system_boxes["iqtree"].toolTip()
