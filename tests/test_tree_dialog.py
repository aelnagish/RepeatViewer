"""Tree dialog: options round-trip, UI rules (Cancel / Save / Run, Options tab)."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QPushButton

from repeatviewer.job_options import (
    TREE_ENGINE_FASTME,
    TREE_RIGHT_CLICK_DIALOG,
    TREE_ROOT_OUTGROUP,
    TreeOptions,
)


def _dialog(qapp, tmp_path, **kw):
    from repeatviewer.views.tree_dialog import TreeParamsDialog

    fa = tmp_path / "sat.fa"
    fa.write_text(">m1\nACGTACGTAC\n>m2\nACGTACGTAA\n>m3\nACGAACGTAA\n>m4\nACGAAC\n")
    return TreeParamsDialog(files=[fa], **kw)


def test_default_options_round_trip(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path, initial=TreeOptions(threads=2))
    assert dlg.options() == TreeOptions(threads=2)


def test_custom_options_round_trip(qapp, tmp_path):
    opts = TreeOptions(
        engine=TREE_ENGINE_FASTME, model="GTR+F+G4", runs=3, seed=7, threads=3,
        min_coverage=0.5, adjust_direction=False, mafft_strategy="ginsi",
        root_mode=TREE_ROOT_OUTGROUP, outgroup=["m2"], alrt=0, ufboot=2000,
        distance_bootstrap=250, right_click_action=TREE_RIGHT_CLICK_DIALOG,
        show_process_window=True, use_system_mafft=True, use_system_path=True,
        mafft_gap_open=2.0, mafft_offset=0.1,
    )
    dlg = _dialog(qapp, tmp_path, initial=opts)
    assert dlg.options() == opts


def test_buttons_are_cancel_save_run_in_order(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    labels = [b.text() for b in dlg._button_box.findChildren(QPushButton)]
    assert labels == ["Restore Defaults", "Cancel", "Save", "Run"]
    assert "OK" not in labels


def test_tabs_and_options_tab(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path)
    tabs = [dlg._tabs.tabText(i) for i in range(dlg._tabs.count())]
    assert tabs == ["Alignment", "Tree", "Support and rooting", "Options"]
    assert dlg._show_process_window is not None and dlg._use_system_path is not None


def test_input_check_and_outgroup_list(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path, initial=TreeOptions(root_mode=TREE_ROOT_OUTGROUP))
    assert "4</b> sequences" in dlg._check.text()
    assert "1</b> fragment" in dlg._check.text()  # m4 is 6 bp vs median 10 bp
    names = [dlg._outgroup_list.item(i).text() for i in range(dlg._outgroup_list.count())]
    assert names == ["m1", "m2", "m3", "m4"]
    dlg._outgroup_list.item(2).setCheckState(Qt.CheckState.Checked)
    assert dlg.options().outgroup == ["m3"]
    dlg._min_cov.setValue(0)
    assert "fragment" not in dlg._check.text()


def test_restore_defaults(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path, initial=TreeOptions(engine=TREE_ENGINE_FASTME, runs=9))
    dlg._restore_defaults()
    o = dlg.options()
    assert (o.engine, o.runs, o.model) == ("iqtree", 5, "MFP")


def test_no_selection_message(qapp):
    from repeatviewer.views.tree_dialog import TreeParamsDialog

    dlg = TreeParamsDialog(files=[])
    assert "No FASTA file selected" in dlg._check.text()
    assert dlg.options().outgroup == []
