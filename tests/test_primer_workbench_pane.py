"""Primer workbench Design tab (offscreen)."""

import random
import time

import pytest
from PySide6.QtCore import QCoreApplication, QSettings


@pytest.fixture()
def isolated_settings(qapp, tmp_path):
    org, app = QCoreApplication.organizationName(), QCoreApplication.applicationName()
    QCoreApplication.setOrganizationName("RVWorkbenchTests")
    QCoreApplication.setApplicationName(f"t{id(tmp_path)}")
    QSettings().clear()
    yield
    QSettings().clear()
    QCoreApplication.setOrganizationName(org)
    QCoreApplication.setApplicationName(app)


class _Host:
    def __init__(self, path, seq):
        self.path, self.seq, self.shown = path, seq, None
        self.features = []

    def primers_host_path(self):
        return self.path

    def primers_host_sequence(self):
        return "t", self.seq

    def primers_host_features(self):
        return list(self.features)

    def primers_host_apply_features(self, features):
        self.features = list(features)

    def primers_host_show(self, start, end):
        self.shown = (start, end)

    def primers_host_set_marks(self, positions):
        pass


def _wait(app, pane, seconds=20):
    end = time.time() + seconds
    while pane._busy and time.time() < end:
        app.processEvents()
        time.sleep(0.02)


def test_design_tab_region_advice_and_keep(qapp, isolated_settings, tmp_path):
    from repeatviewer.widgets.primers_panel import PrimersPanel

    r = random.Random(3)
    seq = "".join(r.choice("ACGT") for _ in range(3000))
    fa = tmp_path / "t.fa"
    fa.write_text(">t\n" + seq + "\n")
    panel = PrimersPanel(width=330)
    host = _Host(fa, seq)
    panel.set_host(host)
    panel.reload()
    panel.show_design_tab()
    pane = panel.design_pane
    assert panel.pages.currentWidget() is pane
    panel.on_sequence_selection(500, 1999)          # 1,500 bp, product max 1,000
    assert pane.region().length == 1500
    assert not pane.advice.isHidden() and pane.fit_button.isVisibleTo(pane)
    pane.fit_button.click()
    assert pane.advice.isHidden()
    panel.on_sequence_selection(1400, 1600)
    pane.prod_min.setValue(250)
    pane.prod_max.setValue(600)
    pane.run_design()
    _wait(qapp, pane)
    assert pane.table.rowCount() > 0
    assert host.shown is not None and host.shown[0] <= 1400 and host.shown[1] >= 1600
    pane.keep_button.click()
    assert len(panel.primers()) == 2 and panel.pages.currentIndex() == 1
    panel.on_sequence_selection(10, 10)              # a click is not a region
    assert pane.region().start == 1400


def test_overlay_follows_region_candidates_and_tab(qapp, isolated_settings, tmp_path):
    from repeatviewer.widgets.primers_panel import PrimersPanel

    r = random.Random(4)
    seq = "".join(r.choice("ACGT") for _ in range(2000))
    fa = tmp_path / "t.fa"
    fa.write_text(">t\n" + seq + "\n")
    calls = []

    class Host(_Host):
        def primers_host_set_workbench(self, region, cands, selected=-1):
            calls.append((region, list(cands), selected))

        def primers_host_show_range(self, start, end):
            self.shown = (start, end)

    host = Host(fa, seq)
    panel = PrimersPanel(width=330)
    panel.set_host(host)
    panel.reload()
    panel.show_design_tab()
    pane = panel.design_pane
    panel.on_sequence_selection(900, 1000)
    assert calls[-1][0] == (900, 1000, "amplify") and calls[-1][1] == []
    pane.prod_min.setValue(150)
    pane.prod_max.setValue(500)
    pane.run_design()
    _wait(qapp, pane)
    region, cands, selected = calls[-1]
    assert cands and selected == 0
    assert all(f_s < 900 and r_e > 1000 for f_s, _fe, _rs, r_e in cands)
    assert pane.region().start == 900          # showing a candidate keeps the region
    panel.tabs.setCurrentIndex(1)
    assert calls[-1] == (None, [], -1)
