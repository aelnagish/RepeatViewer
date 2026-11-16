"""Primer workbench Design tab (offscreen)."""

import random
import time

import pytest
from PySide6.QtCore import QCoreApplication, QSettings, Qt


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


def test_drag_primer_preview_apply_and_revert(qapp, isolated_settings, tmp_path):
    from repeatviewer.backends.primers import on_sequence as osq
    from repeatviewer.backends.primers.oligos import Primer
    from repeatviewer.widgets.primers_panel import PrimersPanel

    r = random.Random(6)
    seq = "".join(r.choice("ACGT") for _ in range(800))
    fa = tmp_path / "t.fa"
    fa.write_text(">t\n" + seq + "\n")
    edits = []

    class Host(_Host):
        def primers_host_edit_primer(self, span, strand="+", flagged=False):
            edits.append((span, strand))

    panel = PrimersPanel(width=330)
    panel.set_host(Host(fa, seq))
    panel.reload()
    panel.tabs.setCurrentIndex(1)
    panel.add_primers([Primer("P_F", seq[199:221])])
    panel._select_row_of("P_F")
    assert edits[-1] == ((200, 221), "+")
    panel.on_primer_dragged(210, 231, False)         # move by 10
    panel.on_primer_dragged(210, 231, True)
    assert panel._primer("P_F").sequence == seq[209:231]
    assert not panel.apply_button.isHidden()
    panel._revert_drag()
    assert panel._primer("P_F").sequence == seq[199:221]
    panel.on_primer_dragged(200, 226, True)          # 3' end longer
    panel.apply_button.click()
    p = panel._primer("P_F")
    assert p.sequence == seq[199:226] and p.extra.get("edited") is True
    assert [x.sequence for x in osq.load_sequence_primers(fa)] == [seq[199:226]]


def test_site_navigator_and_amber_flags(qapp, isolated_settings, tmp_path):
    from repeatviewer.backends.primers.oligos import Primer
    from repeatviewer.widgets.primers_panel import PrimersPanel

    r = random.Random(8)
    mono = "".join(r.choice("ACGT") for _ in range(170))
    seq = "".join(r.choice("ACGT") for _ in range(100)) + mono * 4
    fa = tmp_path / "t.fa"
    fa.write_text(">t\n" + seq + "\n")
    edits = []

    class Host(_Host):
        def primers_host_edit_primer(self, span, strand="+", flagged=False):
            edits.append((span, flagged))

        def primers_host_show_range(self, start, end):
            self.shown = (start, end)

    host = Host(fa, seq)
    panel = PrimersPanel(width=330)
    panel.set_host(host)
    panel.reload()
    panel.tabs.setCurrentIndex(1)
    panel.add_primers([Primer("S_F", mono[30:44])])      # 14 nt: short, low Tm
    panel._select_row_of("S_F")
    panel._refresh_inspector()
    assert panel.site_count.text() == "<b>1</b> of 4"
    panel.site_next.click()
    assert panel.site_count.text() == "<b>2</b> of 4" and host.shown[0] == 100 + 170 + 31
    assert not panel.flags_label.isHidden() and "shorter than 18" in panel.flags_label.text()
    panel._push_edit_primer()
    assert edits[-1][1] is True


def test_live_pcr_tab(qapp, isolated_settings, tmp_path):
    from repeatviewer.backends.primers import thermo as th
    from repeatviewer.backends.primers.oligos import Primer
    from repeatviewer.widgets.primers_panel import PrimersPanel

    r = random.Random(12)
    seq = "".join(r.choice("ACGT") for _ in range(2500))
    fa = tmp_path / "t.fa"
    fa.write_text(">t\n" + seq + "\n")
    calls = []

    class Host(_Host):
        def primers_host_set_workbench(self, region, cands, selected=-1):
            calls.append((region, list(cands), selected))

        def primers_host_show_range(self, start, end):
            self.shown = (start, end)

    panel = PrimersPanel(width=330)
    panel.set_host(Host(fa, seq))
    panel.reload()
    panel.add_primers([Primer("A_F", seq[299:321]), Primer("A_R", th.reverse_complement(seq[1678:1700]))])
    panel.tabs.setCurrentIndex(2)
    pane = panel.pcr_pane
    assert pane.forward.currentData() == "A_F" and pane.reverse.currentData() == "A_R"
    end = time.time() + 20
    while (pane._result is None or pane._timer.isActive()) and time.time() < end:
        qapp.processEvents()
        time.sleep(0.02)
    prods = pane.products()
    assert prods and prods[0].size == 1401 and (prods[0].start, prods[0].end) == (300, 1700)
    assert calls[-1][1] and calls[-1][2] == 0
    assert "longer" in pane.status_label.text()          # Design range 100 to 1,000
    pane.fit_button.click()
    assert "Ready for PCR" in pane.status_label.text()
    assert panel.design_pane.prod_max.value() >= 1401


def test_enzymes_tab_draws_sites_and_methylation(qapp, isolated_settings, tmp_path):
    from repeatviewer.widgets.primers_panel import PrimersPanel

    r = random.Random(5)
    seq = "".join(r.choice("AT") for _ in range(600))
    seq = seq[:100] + "GAATTC" + seq[106:300] + "AACCGGTT" + seq[308:]
    fa = tmp_path / "t.fa"
    fa.write_text(">t\n" + seq + "\n")
    (tmp_path / "t.bismark.cov").write_text("t\t304\t304\t95\t19\t1\n")   # C of CG in CCGG
    drawn = []

    class Host(_Host):
        def primers_host_set_enzymes(self, sites):
            drawn.append(list(sites))

        def primers_host_show_range(self, start, end):
            self.shown = (start, end)

    QSettings().setValue("primers/workbench_enzymes",
                         '{"enzyme_set": "methyl", "show": "many"}')
    panel = PrimersPanel(width=330)
    panel.set_host(Host(fa, seq))
    panel.reload()
    pane = panel.enzyme_pane
    assert pane.data_combo.findText("t.bismark.cov") > 0
    pane.data_combo.setCurrentIndex(pane.data_combo.findText("t.bismark.cov"))
    panel.tabs.setCurrentIndex(3)

    def settle():
        end = time.time() + 20
        while (pane._timer.isActive() or not pane._rows) and time.time() < end:
            qapp.processEvents()
            time.sleep(0.02)
        for _ in range(5):
            qapp.processEvents()
            time.sleep(0.02)

    settle()
    names = {s[3]: s for s in drawn[-1]}
    assert names["HpaII"][4] is True and names["MspI"][4] is False
    assert "probably not cut" in names["HpaII"][5].lower()
    assert names["EcoRI"][2] == 102
    assert not pane.pair_note.isHidden()
    # untick EcoRI: gone from the map
    row = [pane.table.item(i, 0).text() for i in range(pane.table.rowCount())].index("EcoRI")
    pane.table.item(row, 0).setCheckState(Qt.CheckState.Unchecked)
    assert "EcoRI" not in {s[3] for s in drawn[-1]}
    panel.tabs.setCurrentIndex(1)
    assert drawn[-1] == []
    pane.digest_button.click()
    assert list(tmp_path.glob("t_t_digest*.rvgel"))
