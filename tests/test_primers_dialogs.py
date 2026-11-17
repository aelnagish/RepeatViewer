"""Primers and Probes: the menu, the Database window, saved options and the job worker.

Design, Map and PCR are the panel and its tabs now; their tests live in
test_primer_workbench_pane.py, test_map_page.py and test_primers_workbench.py.
"""

from __future__ import annotations

from pathlib import Path

import pytest

pytest.importorskip("primer3")

from PySide6.QtCore import QCoreApplication, QSettings  # noqa: E402

from primer_testdata import make_data  # noqa: E402

from repeatviewer.backends.primers import design as dz  # noqa: E402
from repeatviewer.backends.primers import jobs as pj  # noqa: E402
from repeatviewer.job_options import (  # noqa: E402
    PRIMER_RIGHT_CLICK_RUN,
    PrimerOptions,
    load_primer_options,
    save_primer_options,
)


@pytest.fixture()
def isolated_settings(qapp, tmp_path):
    org, app = QCoreApplication.organizationName(), QCoreApplication.applicationName()
    QCoreApplication.setOrganizationName("RVPrimerTests")
    QCoreApplication.setApplicationName(f"t{id(tmp_path)}")
    QSettings().clear()
    yield
    QSettings().clear()
    QCoreApplication.setOrganizationName(org)
    QCoreApplication.setApplicationName(app)


@pytest.fixture()
def data(tmp_path):
    return make_data(tmp_path / "d", genome_flank=5_000)


def test_database_dialog_lists_and_analyses(qapp, isolated_settings, data):
    from repeatviewer.views.primer_dialogs import PrimerDatabaseDialog

    db_path = data.primers
    d = PrimerDatabaseDialog(db_path)
    assert d._table.rowCount() == 4
    d._table.selectRow(0)
    qapp.processEvents()
    assert d._form_ids and not hasattr(d, "_details")   # no check panel any more
    d._filter.setText("Sat")
    assert d._table.rowCount() == 2


def test_primer_options_persist(qapp, isolated_settings):
    opts = PrimerOptions(design=dz.DesignSettings(mode=dz.MODE_FISH, fish_label="Cy5").to_dict(),
                         database_path="/x/db.fa", save_design_to_database=True,
                         right_click_action=PRIMER_RIGHT_CLICK_RUN)
    save_primer_options(opts)
    back = load_primer_options()
    assert back.design_settings().mode == dz.MODE_FISH
    assert back.design_settings().fish_label == "Cy5"
    assert back.database_path == "/x/db.fa" and back.save_design_to_database
    assert back.right_click_action == PRIMER_RIGHT_CLICK_RUN
    # Unknown / old keys never break loading.
    QSettings().setValue("primers/design_json", "{not json")
    assert load_primer_options().design_settings().mode == dz.MODE_PCR


def test_worker_runs_and_reports(qapp, isolated_settings, data):
    from repeatviewer.primer_actions import PRIMER_JOB_DESIGN, PRIMER_JOB_MAP, PrimerJobWorker

    got, failed, steps = [], [], []
    w = PrimerJobWorker(PRIMER_JOB_DESIGN, [data.gene], dz.DesignSettings(), report_path=None,
                        save_to_database=True, database_path=str(data.root / "db.fa"))
    w.finished.connect(got.append)
    w.failed.connect(failed.append)
    w.progress.connect(lambda v, t: steps.append(v))
    w.run()
    assert not failed and got[0]["pairs"] == 5 and got[0]["saved_to_database"] == 10
    assert Path(got[0]["report"]).name == "gene.fa_primers.html"
    assert max(steps) <= 3 and steps[0] == 0
    w2 = PrimerJobWorker(PRIMER_JOB_MAP, [data.gene], pj.MapSettings(forward=""), report_path=None)
    w2.failed.connect(failed.append)
    w2.run()
    assert failed and "choose a saved primer" in failed[0]


def test_worker_cancel(qapp, isolated_settings, data):
    from repeatviewer.backends.job_cancel import CANCELLED_MESSAGE, JobCancellation
    from repeatviewer.primer_actions import PRIMER_JOB_DESIGN, PrimerJobWorker

    cancel = JobCancellation()
    cancel.request()
    failed = []
    w = PrimerJobWorker(PRIMER_JOB_DESIGN, [data.gene], dz.DesignSettings(), report_path=None,
                        cancel=cancel)
    w.failed.connect(failed.append)
    w.run()
    assert failed == [CANCELLED_MESSAGE]


def test_menu_texts_and_file_family(qapp):
    from PySide6.QtWidgets import QMenu

    from repeatviewer.browser_file_groups import file_family_key
    from repeatviewer.primer_actions import (
        DESIGN_SUBMENU_ITEMS,
        MAP_SUBMENU_ITEMS,
        PRIMER_MANUAL,
        PRIMER_MENU_ITEMS,
        PrimerActionsMixin,
    )

    class _Host(PrimerActionsMixin):
        pass

    menu = QMenu()
    _Host()._populate_primers_menu(menu)
    texts = [a.text() for a in menu.actions() if not a.isSeparator()]
    # the submenus come first, then the items; the PCR mixin adds the rest
    assert texts == ["Design", "Map", "Database"]
    submenus = [a.menu() for a in menu.actions() if a.menu() is not None]
    assert [a.text() for a in menu.actions() if a.menu() is not None] == ["Design", "Map"]
    assert [x.text() for x in submenus[0].actions()] == list(DESIGN_SUBMENU_ITEMS)
    assert [x.text() for x in submenus[1].actions()] == list(MAP_SUBMENU_ITEMS)
    assert list(PRIMER_MENU_ITEMS)[:3] == ["Design", "Map", "Database"]
    assert all("." not in t and "…" not in t for t in texts)
    assert PRIMER_MANUAL.is_file()
    html = PRIMER_MANUAL.read_text(encoding="utf-8")
    assert "—" not in html
    for label in dz.MODE_LABELS.values():
        assert label in html, label
    assert file_family_key("gene.fa_primers.html") == "gene.fa"
    assert file_family_key("gene.fa_primer_map_2.html") == "gene.fa"


def test_designed_dialog_lists_fish_probes(qapp, tmp_path):
    from repeatviewer.backends.primers.on_sequence import load_sequence_primers
    from repeatviewer.views.primer_edit_dialog import DesignedPairsDialog

    fasta = tmp_path / "mono.fa"
    seq = "ACGTTGCAAGCTTAGCCTAGGATCCAAGTCGATCGGATCCTTAGCAAGTCCGATTACGATCGGCTA"
    fasta.write_text(">m1\n" + seq + "\n")
    probe = {"name": "consensus_FISH1", "sequence": seq[5:45], "start": 6, "end": 45,
             "length": 40, "tm": 52.3}
    d = DesignedPairsDialog([{"name": "consensus_FISH1", "template": "consensus",
                              "probe": probe}], [fasta])
    assert d.windowTitle() == "Designed Probes" and d.table.rowCount() == 1
    assert d.table.horizontalHeaderItem(1).text() == "Probe"
    assert d.add_ticked() == [fasta]
    assert [p.name for p in load_sequence_primers(fasta)] == ["consensus_FISH1"]
    d.close()


def test_tail_box_pads_six_bases_and_warns(qapp, isolated_settings):
    from repeatviewer.backends.primers import tails as tl
    from repeatviewer.backends.primers.oligos import Primer
    from repeatviewer.views.primer_edit_dialog import PAD_AUTO, EditPrimerDialog

    bind = "ATGCTAGCTAGGTTCAGATCC"
    template = "TTTT" + bind + "AAAA" + "GAATTC" + "CCCC"
    dlg = EditPrimerDialog(Primer("A_F", bind), template=template)
    assert dlg.pad_combo.currentData() == PAD_AUTO      # 6 bases by default
    dlg.enzyme_combo.setCurrentIndex(
        [dlg.enzyme_combo.itemData(i) for i in range(dlg.enzyme_combo.count())].index("BamHI"))
    tail = dlg._chosen_tail()
    assert len(tail) == tl.PAD_BASES + 6 and tail.endswith("ggatcc")
    assert "does not cut the open sequence" in dlg.tail_note.text()

    i = [dlg.enzyme_combo.itemData(x) for x in range(dlg.enzyme_combo.count())].index("EcoRI")
    dlg.enzyme_combo.setCurrentIndex(i)
    assert "also cuts the open sequence" in dlg.tail_note.text()
    assert "Not cutting here" in dlg.tail_note.text()

    j = [dlg.pad_combo.itemData(x) for x in range(dlg.pad_combo.count())].index("AT")
    dlg.pad_combo.setCurrentIndex(j)
    assert "Only 2 bases" in dlg.tail_note.text()
    dlg.apply_tail.click()
    assert dlg._clean().startswith("atgaattc")
    dlg.deleteLater()


def test_menu_items_open_the_panel_not_a_window(qapp, tmp_path):
    """Design and Restriction Digest show the sequence with the right panel page."""
    from PySide6.QtWidgets import QMenu

    from repeatviewer.backends.primers import design as dz
    from repeatviewer.pcr_actions import PcrActionsMixin
    from repeatviewer.primer_actions import PrimerActionsMixin

    opened: list[tuple] = []

    class _Host(PrimerActionsMixin, PcrActionsMixin):
        def _primers_ready(self):
            return True

        def _open_primers_panel(self, *, tab="design", mode=None, digest=False):
            opened.append((tab, mode, digest))
            return None

    host = _Host()
    menu = QMenu()
    host._populate_primers_menu(menu)
    design = [a for a in menu.actions() if a.text() == "Design"][0].menu()
    [a for a in design.actions() if a.text() == "Primers"][0].trigger()
    [a for a in design.actions() if a.text() == "Repetitive Probes"][0].trigger()
    [a for a in menu.actions() if a.text() == "Restriction Digest"][0].trigger()
    assert opened == [("design", dz.MODE_PCR, False),
                      ("design", dz.MODE_FISH, False),
                      ("enzymes", None, True)]


def test_map_for_probes_drops_the_pcr_rules(qapp):
    """Map > Repetitive Probes searches probes only, with no 3' end rule."""
    from dataclasses import replace

    from repeatviewer.job_options import load_primer_options

    opts = load_primer_options()
    ms = dict(opts.map)
    ms.update(search_forward=False, search_reverse=False, search_probe=True,
              pairs_only=False, three_prime_exact=0, use_max_mismatches=True)
    probe_opts = replace(opts, map=ms)
    assert probe_opts.map["search_probe"] is True
    assert probe_opts.map["search_forward"] is False
    assert probe_opts.map["three_prime_exact"] == 0

    from repeatviewer.backends.primers import jobs as pj

    settings = pj.MapSettings.from_dict(probe_opts.map)
    assert settings.search_probe and not settings.search_forward
    assert settings.three_prime_exact == 0 and not settings.pairs_only

def test_mode_list_is_grouped(qapp, isolated_settings, tmp_path):
    """Every mode sits in one family, and the Design pane offers them all."""
    from repeatviewer.widgets.primers_panel import PrimersPanel

    grouped = [m for _f, groups in dz.MODE_FAMILIES for _g, modes in groups for m in modes]
    assert sorted(grouped) == sorted(dz.ALL_MODES) and len(grouped) == len(set(grouped))
    panel = PrimersPanel(width=330)
    pane = panel.design_pane
    offered = [pane.mode.itemData(i) for i in range(pane.mode.count())]
    assert [m for m in offered if m] and set(m for m in offered if m) <= set(grouped)
    assert pane.mode.currentData() in grouped


def test_toolbar_ids_renamed(qapp, isolated_settings):
    from repeatviewer.views.toolbar_dialog import (
        TOOLBAR_LABELS,
        TOOLS_TOOLBAR_IDS,
        _rename_old_toolbar_ids,
    )

    assert "primers" in TOOLS_TOOLBAR_IDS and TOOLBAR_LABELS["primers"] == "Primers and Probes"
    assert "primer_design" not in TOOLBAR_LABELS and "primer_map" not in TOOLBAR_LABELS
    assert _rename_old_toolbar_ids(["blast", "primer_design", "tree", "primer_map"]) == [
        "blast", "primers", "tree"]
    assert _rename_old_toolbar_ids(["primer_map"]) == ["primers"]


def test_the_old_windows_are_gone():
    """Design, Map and PCR live in the panel and in tabs now, not in windows."""
    import importlib

    from repeatviewer.views import primer_dialogs

    assert not hasattr(primer_dialogs, "PrimerDesignDialog")
    assert not hasattr(primer_dialogs, "MapPrimersDialog")
    with pytest.raises(ModuleNotFoundError):
        importlib.import_module("repeatviewer.views.pcr_dialog")


def test_design_opens_the_workspace_behind_a_tool_tab(qapp, tmp_path):
    """A tool tab that comes and goes must not send Design to the wrong tab."""
    import random

    pytest.importorskip("PySide6.QtPdf")   # the preview tabs need it

    from PySide6.QtWidgets import QLabel, QWidget

    from repeatviewer.pcr_actions import PcrActionsMixin
    from repeatviewer.primer_actions import PrimerActionsMixin
    from repeatviewer.views.primers_workspace import PrimersWorkspacePage
    from repeatviewer.widgets.preview_tabs import PreviewTabWidget

    r = random.Random(3)
    fa = tmp_path / "CroSat4.fasta"
    fa.write_text(">CroSat4\n" + "".join(r.choice("ACGT") for _ in range(209)) + "\n")
    tabs = PreviewTabWidget()

    class Host(PrimerActionsMixin, PcrActionsMixin, QWidget):
        def __init__(self):
            super().__init__()
            self._preview_tabs = tabs
            self._primer_options = load_primer_options()

        def _selected_files(self):
            return [fa]

        def _ensure_preview_visible(self):
            pass

        def _set_status_message(self, *a, **k):
            pass

    host = Host()
    tabs.open_tool_tab(QLabel("one"), "Sanger Check")
    host._design_in_panel(dz.MODE_PCR)
    assert isinstance(tabs.currentWidget(), PrimersWorkspacePage)

    # The same tool runs again: its old tab goes, every tab after it moves up.
    tabs.open_tool_tab(QLabel("two"), "Sanger Check")
    host._design_in_panel(dz.MODE_PCR)
    page = tabs.currentWidget()
    assert isinstance(page, PrimersWorkspacePage) and page.sequence_path() == fa.resolve()
    assert page.primers_panel().tabs.currentIndex() == 0
    # A window left behind would be found by find_main_window in later tests.
    from PySide6.QtCore import QCoreApplication, QEvent
    from PySide6.QtTest import QTest

    QTest.qWait(200)             # queued refreshes finish before their widgets go
    for widget in (host, tabs):
        widget.close()
        widget.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)


def test_toolbar_popup_lists_what_is_inside_design_and_map(qapp):
    """The toolbar popup is a flat list: a submenu must show its items, not a dead row."""
    pytest.importorskip("PySide6.QtPdf")   # the main window needs it
    from PySide6.QtWidgets import QLabel, QMenu, QPushButton

    from repeatviewer.main_window import _ToolbarPopup

    fired: list[str] = []
    menu = QMenu()
    design = menu.addMenu("Design")
    design.addAction("Primers", lambda: fired.append("design primers"))
    design.addAction("Repetitive Probes", lambda: fired.append("design probes"))
    map_menu = menu.addMenu("Map")
    map_menu.addAction("Primers", lambda: fired.append("map primers"))
    menu.addAction("Database", lambda: fired.append("database"))

    popup = _ToolbarPopup(list(menu.actions()))
    rows = [w for w in popup._panel.children() if isinstance(w, (QLabel, QPushButton))]
    assert [(type(w).__name__, w.text()) for w in rows] == [
        ("QLabel", "Design"), ("QPushButton", "Primers"), ("QPushButton", "Repetitive Probes"),
        ("QLabel", "Map"), ("QPushButton", "Primers"), ("QPushButton", "Database")]
    rows[4].click()                              # Map > Primers
    assert fired == ["map primers"]
