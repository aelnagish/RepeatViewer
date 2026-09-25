"""Tree dialog: essentials up front, More options, Options tab, round trips."""

from __future__ import annotations

from PySide6.QtWidgets import QPushButton, QRadioButton, QToolButton

from repeatviewer.job_options import (
    TREE_ALIGN_ALWAYS,
    TREE_ALIGN_AUTO,
    TREE_ALIGN_NEVER,
    TREE_ENGINE_FASTME,
    TREE_ENGINE_FASTTREE,
    TREE_RIGHT_CLICK_DIALOG,
    TREE_ROOT_MIDPOINT,
    TREE_ROOT_OUTGROUP,
    TREE_SEQ_RT,
    TreeOptions,
)


def _dialog(qapp, tmp_path, **kw):
    from repeatviewer.views.tree_dialog import TreeParamsDialog

    fa = tmp_path / "sat.fa"
    fa.write_text(">m1\nACGTACGTAC\n>m2\nACGTACGTAA\n>m3\nACGAACGTAA\n>m4\nACGAAC\n")
    return TreeParamsDialog(files=[fa], **kw)


def _shown(dlg, widget) -> bool:
    return widget.isVisibleTo(dlg)


def _normalised(opts: TreeOptions) -> TreeOptions:
    """What the window writes back: MAFFT is chosen automatically, one PATH switch."""
    from dataclasses import replace

    from repeatviewer.backends import tree_builder as tb

    system = opts.use_system_path or opts.use_system_mafft
    return replace(opts, align_mode=TREE_ALIGN_ALWAYS if opts.align_mode == TREE_ALIGN_AUTO
                   else opts.align_mode, mafft_strategy=tb.MAFFT_CHOOSE, adjust_direction=True,
                   mafft_gap_open=1.53, mafft_offset=0.0, bnni=True,
                   use_system_mafft=system, use_system_path=system)


def test_default_options_round_trip(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path, initial=TreeOptions(threads=2))
    assert dlg.options() == _normalised(TreeOptions(threads=2))
    assert dlg._support_iq.currentText() == "SH-aLRT + UFBoot, 1000 each"   # not a custom entry
    assert dlg._support_iq.count() == 4


def test_custom_options_round_trip_keeps_custom_support(qapp, tmp_path):
    opts = TreeOptions(
        engine=TREE_ENGINE_FASTME, model="GTR+F+G4", runs=3, seed=7, threads=3,
        min_coverage=0.5, root_mode=TREE_ROOT_OUTGROUP, outgroup=["m2"], alrt=0, ufboot=2000,
        distance_bootstrap=250, right_click_action=TREE_RIGHT_CLICK_DIALOG,
        show_process_window=True, use_system_mafft=True, coverage_reference="longest",
        align_mode=TREE_ALIGN_NEVER,
    )
    dlg = _dialog(qapp, tmp_path, initial=opts)
    assert dlg.options() == _normalised(opts)
    assert dlg._support_iq.currentText() == "SH-aLRT 0, UFBoot 2000"
    assert dlg._support_fme.currentText() == "Bootstrap, 250 replicates"


def test_buttons_tabs_and_options_tab(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    labels = [b.text() for b in dlg._button_box.findChildren(QPushButton)]
    assert labels == ["Restore Defaults", "Cancel", "Save", "Run"] and "OK" not in labels
    assert [dlg._tabs.tabText(i) for i in range(dlg._tabs.count())] == ["Tree", "Options"]
    assert dlg._show_process_window is not None and dlg._use_system_path is not None
    assert "IQ-TREE:" in dlg._found.text() and "MAFFT:" in dlg._found.text()


def test_essentials_up_front_rest_under_more_options(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    for w in (dlg._seq_mode, dlg._align_mode, dlg._rotate, dlg._support_stack, dlg._root_mode):
        assert _shown(dlg, w)
    for w in (dlg._min_cov, dlg._model, dlg._mask, dlg._runs, dlg._seed, dlg._threads):
        assert not _shown(dlg, w)
    dlg._more_toggle.setChecked(True)
    for w in (dlg._min_cov, dlg._model, dlg._mask, dlg._runs, dlg._seed, dlg._threads):
        assert _shown(dlg, w)
    assert "1 of 4 sequences" in dlg._frag_hint.text()    # m4 is 6 bp vs median 10 bp
    dlg._min_cov.setValue(0)
    assert "None of the 4" in dlg._frag_hint.text()


def test_rows_that_do_not_apply_are_hidden(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    dlg._more_toggle.setChecked(True)
    assert dlg._support_stack.currentWidget() is dlg._support_iq
    dlg._engine_buttons[TREE_ENGINE_FASTME].setChecked(True)
    assert dlg._support_stack.currentWidget() is dlg._support_fme
    assert not _shown(dlg, dlg._model) and not _shown(dlg, dlg._runs)
    dlg._engine_buttons[TREE_ENGINE_FASTTREE].setChecked(True)
    assert dlg._support_stack.currentWidget() is dlg._support_ft
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
    assert not dlg._outgroup_btn.isEnabled() and not _shown(dlg, dlg._outgroup_note)
    dlg._root_mode.setCurrentIndex(dlg._root_mode.findData(TREE_ROOT_OUTGROUP))
    assert dlg._outgroup_btn.isEnabled() and "None chosen" in dlg._outgroup_note.text()
    dlg._outgroup = ["m3"]
    dlg._sync()
    assert "1 chosen: m3" in dlg._outgroup_note.text()
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
    assert "Maximum likelihood (IQ-TREE)" in radios and "Minimum evolution (FastME)" in radios
    assert not [lab for lab in dlg.findChildren(QLabel) if lab.text() == "ⓘ"]
    helps = [b for b in dlg.findChildren(QToolButton) if b.text() == "?"]
    # Sequences, REXdb, Alignment, Monomers, Method, Branch support, Root,
    # Skip fragments, Model, Methylation check, Searches
    assert len(helps) == 11
    joined = " ".join(b.toolTip() for b in helps)
    for piece in ("L-INS-i", "FFT-NS-2", "Align tool", "ModelFinder", "UFBoot", "Midpoint", "CHG"):
        assert piece in joined


def test_protein_input_hides_dna_options(qapp, tmp_path):
    from repeatviewer.views.tree_dialog import TreeParamsDialog

    fa = tmp_path / "p.fa"
    fa.write_text("".join(f">p{i}\nMKVLLAQEEFWRSTLK\n" for i in range(4)))
    dlg = TreeParamsDialog(files=[fa], initial=TreeOptions(methylation_mask="CG"))
    dlg._more_toggle.setChecked(True)
    assert not _shown(dlg, dlg._rotate) and not _shown(dlg, dlg._mask)
    assert not dlg._seq_mode.isEnabled()
    assert dlg.options().methylation_mask == ""


def test_restore_defaults(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path, initial=TreeOptions(engine=TREE_ENGINE_FASTME, runs=9,
                                                      rotate_monomers=True, methylation_mask="CG"))
    dlg._restore_defaults()
    o = dlg.options()
    assert (o.engine, o.runs, o.model, o.rotate_monomers, o.methylation_mask) == (
        "iqtree", 5, "MFP", False, "")


def test_opens_without_selection(qapp):
    from repeatviewer.views.tree_dialog import TreeParamsDialog

    dlg = TreeParamsDialog(files=[])
    dlg._more_toggle.setChecked(True)
    assert not _shown(dlg, dlg._frag_hint)
    assert dlg.options().outgroup == []


def test_rt_domain_mode(qapp, tmp_path):
    rex = tmp_path / "Viridiplantae_v4.0.fasta"
    rex.write_text(">Ty1-RT__REXdb_ID1\nMKV\n")
    (tmp_path / "Viridiplantae_v4.0.classification").write_text("REXdb_ID1\tClass_I\tAle\n")
    opts = TreeOptions(threads=2, sequence_mode=TREE_SEQ_RT, rexdb_fasta=str(rex),
                       rexdb_references=5, rotate_monomers=True)
    dlg = _dialog(qapp, tmp_path, initial=opts)
    assert dlg.options() == _normalised(opts)
    assert _shown(dlg, dlg._rexdb_refs) and not _shown(dlg, dlg._rotate)
    assert "Found" in dlg._rexdb_status.text()
    other = tmp_path / "copies.fa"
    other.write_text(">a\nACGT\n")
    assert not dlg.apply_path_from_panel(other)
    assert dlg.apply_path_from_panel(rex)
    dlg._seq_mode.setCurrentIndex(0)
    assert _shown(dlg, dlg._rotate) and not _shown(dlg, dlg._rexdb_refs)
