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
        mafft_gap_open=2.0, mafft_offset=0.1, coverage_reference="longest",
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


def test_fragment_hint_and_outgroup_list(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path, initial=TreeOptions(root_mode=TREE_ROOT_OUTGROUP))
    assert not hasattr(dlg, "_check")  # no input strip above the tabs
    assert "1 of 4 sequences" in dlg._frag_hint.text()  # m4 is 6 bp vs median 10 bp
    names = [dlg._outgroup_list.item(i).text() for i in range(dlg._outgroup_list.count())]
    assert names == ["m1", "m2", "m3", "m4"]
    dlg._outgroup_list.item(2).setCheckState(Qt.CheckState.Checked)
    assert dlg.options().outgroup == ["m3"]
    dlg._min_cov.setValue(0)
    assert "None of the 4" in dlg._frag_hint.text()


def test_professional_labels_and_info_icons(qapp, tmp_path):
    from PySide6.QtWidgets import QLabel, QRadioButton

    dlg = _dialog(qapp, tmp_path)
    radios = [r.text() for r in dlg.findChildren(QRadioButton)]
    assert "Maximum likelihood (IQ-TREE)" in radios
    assert "Minimum evolution (FastME)" in radios
    assert not any("recommended" in r.lower() and "IQ-TREE" in r for r in radios)
    tips = [l.toolTip() for l in dlg.findChildren(QLabel) if l.text() == "ⓘ"]
    assert len(tips) == 3
    joined = " ".join(tips)
    for piece in ("L-INS-i", "G-INS-i", "E-INS-i", "FFT-NS-2", "ModelFinder", "TN93", "HKY", "K2P"):
        assert piece in joined


def test_inapplicable_options_are_greyed(qapp, tmp_path):
    from repeatviewer.job_options import TREE_ALIGN_NEVER, TREE_ENGINE_FASTTREE

    dlg = _dialog(qapp, tmp_path)
    assert dlg._model.isEnabled() and not dlg._dist_boot.isEnabled()
    dlg._engine_buttons[TREE_ENGINE_FASTME].setChecked(True)
    assert not dlg._model.isEnabled() and not dlg._ufboot.isEnabled()
    assert not dlg._bnni.isEnabled() and dlg._dist_boot.isEnabled()
    dlg._engine_buttons[TREE_ENGINE_FASTTREE].setChecked(True)
    assert not dlg._dist_boot.isEnabled() and not dlg._runs.isEnabled()
    dlg._engine_buttons["iqtree"].setChecked(True)
    dlg._ufboot.setValue(0)
    assert not dlg._bnni.isEnabled()
    dlg._set_combo(dlg._align_mode, TREE_ALIGN_NEVER)
    assert not dlg._strategy.isEnabled() and not dlg._more_box.isEnabled()
    assert not dlg._adjust_direction.isEnabled()
    assert dlg._min_cov.isEnabled()


def test_protein_input_greys_adjust_direction(qapp, tmp_path):
    from repeatviewer.views.tree_dialog import TreeParamsDialog

    fa = tmp_path / "p.fa"
    fa.write_text("".join(f">p{i}\nMKVLLAQEEFWRSTLK\n" for i in range(4)))
    dlg = TreeParamsDialog(files=[fa])
    assert not dlg._adjust_direction.isEnabled()


def test_restore_defaults(qapp, tmp_path):
    dlg = _dialog(qapp, tmp_path, initial=TreeOptions(engine=TREE_ENGINE_FASTME, runs=9))
    dlg._restore_defaults()
    o = dlg.options()
    assert (o.engine, o.runs, o.model) == ("iqtree", 5, "MFP")


def test_opens_without_selection(qapp):
    from repeatviewer.views.tree_dialog import TreeParamsDialog

    dlg = TreeParamsDialog(files=[])
    assert not dlg._frag_hint.isVisible()
    assert dlg.options().outgroup == []
