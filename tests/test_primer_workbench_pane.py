"""Primer workbench Design tab (offscreen)."""

import json
import random
import time
from pathlib import Path

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
    assert panel.pages.currentWidget().widget() is pane
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
    page = pane.digest_page
    assert pane.stack.currentWidget() is page
    labels = [ln.label for ln in page.lanes()]
    assert labels == ["HpaII (your data)", "MspI", "HpaII (no data)"] or \
        labels[0].startswith("HpaII (your data)")
    assert page.gel._doc is not None and page.table.rowCount() > 0
    assert "methylated in your data" in page.why.text()
    page.open_button.click()
    assert list(tmp_path.glob("t_t_digest*.rvgel"))


def test_primer_chips_warn_about_a_restriction_tail(qapp, isolated_settings, tmp_path):
    from repeatviewer.backends.primers import tails as tl
    from repeatviewer.backends.primers import thermo as th
    from repeatviewer.backends.primers.oligos import Primer
    from repeatviewer.widgets.primers_panel import PrimersPanel

    r = random.Random(21)
    seq = "".join(r.choice("ACGT") for _ in range(1200))
    seq = seq[:600] + "GAATTC" + seq[606:]          # EcoRI inside the open sequence
    fa = tmp_path / "t.fa"
    fa.write_text(">t\n" + seq + "\n")
    panel = PrimersPanel(width=330)
    panel.set_host(_Host(fa, seq))
    panel.reload()
    fwd = "atgaattc" + seq[99:121]                  # 2 bases of padding only
    rev = tl.tail_for("EcoRI", binding="x") + th.reverse_complement(seq[1000:1022])
    panel.add_primers([Primer("A_F", fwd), Primer("A_R", rev)], tick=True)
    panel._select_row_of("A_F")
    panel._refresh_inspector()
    chips = panel.flags_label.text()
    assert "Only 2 bases" in chips
    assert "also cuts the open sequence" in chips
    assert "carries EcoRI too" in chips and "A_R" in chips


def test_design_tab_carries_the_window_settings(qapp, isolated_settings, tmp_path):
    from repeatviewer.backends.primers import design as dz
    from repeatviewer.widgets.primer_design_pane import RUN_FILES, RUN_THIS
    from repeatviewer.widgets.primers_panel import PrimersPanel

    r = random.Random(9)
    seq = "".join(r.choice("ACGT") for _ in range(2000))
    fa = tmp_path / "t.fa"
    fa.write_text(">t\n" + seq + "\n")
    panel = PrimersPanel(width=330)
    panel.set_host(_Host(fa, seq))
    panel.reload()
    pane = panel.design_pane

    # every condition of the old window is here and lands in the settings
    pane.gc_min.setValue(35.0)
    pane.gc_max.setValue(65.0)
    pane.tm_diff.setValue(1.5)
    pane.num_return.setValue(7)
    pane.probe_check.setChecked(True)
    pane.poly_x.setValue(3)
    pane.gc_clamp.setValue(1)
    pane.end_stab.setValue(8.0)
    pane.probe_size_opt.setValue(24)
    pane.probe_tm_opt.setValue(68.0)
    pane.na.setValue(60.0)
    pane.mg.setValue(2.0)
    pane.dntp.setValue(0.8)
    s = pane.settings()
    assert (s.gc_min, s.gc_max, s.max_tm_diff, s.num_return) == (35.0, 65.0, 1.5, 7)
    assert s.pick_probe and (s.max_poly_x, s.gc_clamp, s.max_end_stability) == (3, 1, 8.0)
    assert (s.probe_size_opt, s.probe_tm_opt) == (24, 68.0)
    assert (s.thermo.monovalent_mM, s.thermo.divalent_mM, s.thermo.dntp_mM) == (60.0, 2.0, 0.8)

    # Advanced is closed until asked for
    assert pane.adv_body.isHidden() and "Poly-X 3" in pane.adv_summary.text()
    pane.adv_toggle.setChecked(True)
    assert not pane.adv_body.isHidden() and pane.adv_toggle.text() == "Hide"

    # Save writes the settings that right-click running reads
    pane.save_button.click()
    saved = dz.DesignSettings.from_dict(
        json.loads(QSettings().value("primers/workbench_design")))
    assert saved.gc_max == 65.0 and saved.num_return == 7
    assert "Saved" in pane.saved_note.text()

    # Restore Defaults puts the mode's own values back
    pane.defaults_button.click()
    assert pane.settings().gc_max == dz.DesignSettings().gc_max
    assert not pane.probe_check.isChecked()


def test_design_on_selected_files_asks_for_the_job(qapp, isolated_settings, tmp_path):
    from repeatviewer.widgets.primer_design_pane import RUN_FILES
    from repeatviewer.widgets.primers_panel import PrimersPanel

    r = random.Random(10)
    seq = "".join(r.choice("ACGT") for _ in range(1500))
    fa = tmp_path / "t.fa"
    fa.write_text(">t\n" + seq + "\n")
    panel = PrimersPanel(width=330)
    panel.set_host(_Host(fa, seq))
    panel.reload()
    pane = panel.design_pane
    asked = []
    pane.jobRequested.connect(asked.append)

    panel.on_sequence_selection(100, 600)
    assert pane.region() is not None
    for b in pane.run_group.buttons():
        if b.property("run") == RUN_FILES:
            b.click()
    assert pane.region_section.isHidden()        # the job runs on whole files
    assert pane.region() is None
    assert pane.design_button.text() == "Run"
    pane.ask_check.setChecked(True)
    pane.db_check.setChecked(True)
    pane.design_button.click()
    assert len(asked) == 1
    assert asked[0]["annotate_design"] and asked[0]["save_design_to_database"]
    assert asked[0]["settings"]["size_opt"] == pane.size_opt.value()
    assert not pane._busy                        # the panel itself ran nothing


def test_probe_mode_shows_the_fish_rows_and_designs_probes(qapp, isolated_settings, tmp_path):
    from repeatviewer.backends.primers import design as dz
    from repeatviewer.backends.primers import workbench as wb
    from repeatviewer.widgets.primers_panel import PrimersPanel

    r = random.Random(5)
    mono = "".join(r.choice("ACGT") for _ in range(170))
    copies = []
    for _ in range(12):
        m = list(mono)
        for _ in range(4):
            m[r.randrange(170)] = r.choice("ACGT")
        copies.append("".join(m))
    fa = tmp_path / "CroSat3.fa"
    fa.write_text("".join(f">m{i}\n{c}\n" for i, c in enumerate(copies)))

    panel = PrimersPanel(width=330)
    panel.set_host(_Host(fa, copies[0]))
    panel.reload()
    pane = panel.design_pane
    assert not pane.is_probe_mode() and pane.product_row.isVisibleTo(pane)

    modes = [pane.mode.itemData(i) for i in range(pane.mode.count())]
    assert dz.MODE_FISH in modes
    pane.mode.setCurrentIndex(modes.index(dz.MODE_FISH))
    assert pane.mode.currentText() == "Repetitive Probes"
    # PCR rows go, probe rows come
    assert not pane.product_row.isVisibleTo(pane)
    assert not pane.tm_diff_row.isVisibleTo(pane)
    assert not pane.probe_row.isVisibleTo(pane)
    assert pane.cross_section.isVisibleTo(pane) and pane.hyb_section.isVisibleTo(pane)
    assert pane.wash_section.isVisibleTo(pane)
    assert pane.num_label.text() == "Probes"
    shown = {str(b.property("use")) for b in pane.use_group.buttons() if b.isVisibleTo(pane)}
    assert shown == set(wb.PROBE_USES)          # no Amplify for a probe

    # the lab protocol is the default, and it reaches the settings
    s = pane.settings()
    assert (s.hyb.formamide_pct, s.hyb.ssc, s.hyb.temp_c) == (50.0, 2.0, 37.0)
    assert (s.hyb.wash_formamide_pct, s.hyb.wash_ssc, s.hyb.wash_temp_c) == (20.0, 0.1, 42.0)
    pane.min_cov.setValue(90)
    pane.cov_mm.setValue(2)
    s = pane.settings()
    assert s.mode == dz.MODE_FISH and s.min_coverage == 0.9 and s.coverage_mismatches == 2

    pane.run_design()
    end = time.time() + 40
    while pane._busy and time.time() < end:
        qapp.processEvents()
        time.sleep(0.02)
    assert pane.table.rowCount() >= 1
    assert pane.table.horizontalHeaderItem(1).text() == "Probe"
    assert pane.table.item(0, 1).text() == "P1"
    assert "nt" in pane.table.item(0, 2).text() and "°C" in pane.table.item(0, 3).text()
    assert "probes" in pane.cand_title.text().lower()

    # keeping a probe adds one primer, not a pair
    kept = []
    pane.keepRequested.connect(kept.append)
    pane._on_keep()
    assert len(kept) == 1 and len(kept[0]) == 1
    assert "probe" in kept[0][0].note.lower() or kept[0][0].kind == "probe"


def test_pcr_tab_runs_on_several_ticked_files(qapp, isolated_settings, tmp_path, monkeypatch):
    """Run on: Selected files tests one pair on every ticked FASTA, file column and all."""
    from repeatviewer.backends.primers import thermo as th
    from repeatviewer.backends.primers.oligos import Primer
    from repeatviewer.widgets import primer_pcr_pane as pp
    from repeatviewer.widgets.primers_panel import PrimersPanel

    r = random.Random(21)
    seq = "".join(r.choice("ACGT") for _ in range(2000))
    files = []
    for i in range(2):
        fa = tmp_path / f"fam{i}.fa"
        fa.write_text(f">fam{i}\n" + seq + "\n")
        files.append(fa)

    class MW:
        def _pcr_template_candidates(self):
            return list(files)

    monkeypatch.setattr("repeatviewer.widgets.primers_panel.find_main_window",
                        lambda: MW(), raising=False)

    panel = PrimersPanel(width=330)
    panel.set_host(_Host(files[0], seq))
    panel.reload()
    panel.add_primers([Primer("F", seq[199:221]),
                       Primer("R", th.reverse_complement(seq[1098:1120]))])
    panel.tabs.setCurrentIndex(2)
    pane = panel.pcr_pane
    assert pane.run_where() == pp.RUN_THIS

    for b in pane.run_group.buttons():
        if b.property("run") == pp.RUN_FILES:
            b.click()
    assert pane.run_where() == pp.RUN_FILES
    assert pane.templates() == files
    assert pane.table.columnCount() == 5
    assert pane.table.horizontalHeaderItem(1).text() == "File"
    assert "2 files" in pane.run_note.text()

    end = time.time() + 20
    while (pane._result is None or pane._timer.isActive()) and time.time() < end:
        qapp.processEvents()
        time.sleep(0.02)
    prods = pane.products()
    assert len({Path(p.template).name for p in prods}) == 2
    names = {pane.table.item(r, 1).text() for r in range(pane.table.rowCount())}
    assert names == {"fam0.fa", "fam1.fa"}


def test_the_rail_has_no_primers_button(qapp, isolated_settings, tmp_path):
    """Primers and Probes is reached from the menu, not from a rail button."""
    from repeatviewer.widgets.sequence_annotation_rail import (
        PRIMERS_TAB_INDEX,
        SequenceAnnotationRail,
    )

    rail = SequenceAnnotationRail()
    buttons = rail._tab_buttons
    assert [b.isVisible() for b in buttons[:4]] == [True] * 4 or not rail.isVisible()
    assert buttons[PRIMERS_TAB_INDEX].isHidden()
    rail.show_primers_tab()                      # the menu still opens the page
    assert rail._active_tab == PRIMERS_TAB_INDEX
