"""The Primers and Probes workspace tab (offscreen)."""

import random

import pytest
from PySide6.QtCore import QCoreApplication, QSettings
from PySide6.QtWidgets import QLabel, QWidget

from repeatviewer.backends.primers.oligos import Primer


@pytest.fixture()
def isolated_settings(qapp, tmp_path):
    org, app = QCoreApplication.organizationName(), QCoreApplication.applicationName()
    QCoreApplication.setOrganizationName("RVWorkspaceTests")
    QCoreApplication.setApplicationName(f"t{id(tmp_path)}")
    QSettings().clear()
    yield
    QSettings().clear()
    QCoreApplication.setOrganizationName(org)
    QCoreApplication.setApplicationName(app)


def _fasta(path, records):
    path.write_text("".join(f">{name}\n{seq}\n" for name, seq in records))
    return path


def _seq(seed, n=600):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


@pytest.fixture()
def two_files(tmp_path):
    return (_fasta(tmp_path / "knob.fa", [("knob", _seq(1))]),
            _fasta(tmp_path / "other.fasta", [("first", _seq(2)), ("second", _seq(3))]))


def _page(path, **kw):
    from repeatviewer.views.primers_workspace import PrimersWorkspacePage

    return PrimersWorkspacePage(path, **kw)


def test_sections_drive_the_panel(qapp, isolated_settings, two_files):
    from repeatviewer.views import primers_workspace as pw

    page = _page(two_files[0])
    assert page.tab_name() == "Primers and Probes" and page.path is None
    labels = [page.sections.tabText(i) for i in range(page.sections.count())]
    assert labels == ["Design", "PCR", "Digest"]
    tabs = page.primers.tabs

    assert page.section() == pw.SECTION_DESIGN and tabs.currentIndex() == pw.PANEL_DESIGN
    assert tabs.isHidden()                           # the sections do its job
    assert not page.bottom.isHidden() and page.bottom_stack.currentIndex() == 0

    page.set_section("pcr")
    assert page.section() == pw.SECTION_PCR and tabs.currentIndex() == pw.PANEL_PCR
    assert page.primers.pcr_pane.active
    assert not page.bottom.isHidden() and page.bottom_stack.currentIndex() == 1
    assert page.bottom.isAncestorOf(page.primers.pcr_pane.table)   # products under the map

    page.set_section("enzymes")                      # the old name of Digest
    assert page.section() == pw.SECTION_DIGEST and tabs.currentIndex() == pw.PANEL_ENZYMES
    assert page.primers.enzyme_pane.active and not page.primers.pcr_pane.active
    assert page.bottom.isHidden()

    page.sections.setCurrentIndex(0)                 # a click on Design
    assert tabs.currentIndex() == pw.PANEL_DESIGN and tabs.isHidden()
    assert page.primers.design_pane.overlay_enabled
    assert not page.bottom.isHidden() and page.bottom_stack.currentIndex() == 0


def test_the_panel_can_move_the_section(qapp, isolated_settings, two_files):
    from repeatviewer.views import primers_workspace as pw

    page = _page(two_files[0], section="pcr")
    assert page.section() == pw.SECTION_PCR
    page.primers.tabs.setCurrentIndex(pw.PANEL_DESIGN)   # the panel went to Design by itself
    assert page.section() == pw.SECTION_DESIGN
    assert page.bottom_stack.currentIndex() == 0
    # the primer page needs a primer: with none selected the settings stay
    page.primers.tabs.setCurrentIndex(pw.PANEL_PRIMER)
    qapp.processEvents()
    assert page.primers.tabs.currentIndex() == pw.PANEL_DESIGN
    assert pw.section_of_panel_page(3) == pw.SECTION_DIGEST


def test_nothing_can_switch_the_primers_page_away(qapp, isolated_settings, two_files):
    from repeatviewer.widgets.sequence_annotation_rail import PRIMERS_TAB_INDEX

    page = _page(two_files[0])
    rail = page.panel._rail
    assert rail._tab_strip.isHidden()
    assert rail._active_tab == PRIMERS_TAB_INDEX and not rail._panel_collapsed
    assert page.set_sequence(two_files[1])
    assert rail._tab_strip.isHidden()
    assert rail._active_tab == PRIMERS_TAB_INDEX and not rail._panel_collapsed
    # the Annotation view of a file keeps its own remembered rail tab
    from repeatviewer.widgets.sequence_annotation_rail import SequenceAnnotationRail

    assert SequenceAnnotationRail()._active_tab != PRIMERS_TAB_INDEX


def test_sequence_and_record(qapp, isolated_settings, two_files, tmp_path):
    knob, other = two_files
    page = _page(knob)
    seen = []
    page.sequenceChanged.connect(seen.append)
    assert page.sequence_path() == knob.resolve()
    assert page.record_combo.isHidden() and page.record_label.isHidden()   # one sequence

    assert page.set_sequence(other)
    assert page.sequence_path() == other.resolve() and seen == [other.resolve()]
    assert page.sequence_row.path() == other.resolve()
    assert not page.record_combo.isHidden()
    assert [page.record_combo.itemData(i) for i in range(2)] == ["first", "second"]
    assert page.panel.primers_host_sequence()[0] == "first"
    page.record_combo.setCurrentIndex(1)
    page._on_record_chosen(1)
    assert page.panel.primers_host_sequence()[0] == "second"

    assert page.set_sequence(other) and seen == [other.resolve()]          # same file: nothing
    assert not page.set_sequence(tmp_path / "missing.fa")
    assert page.sequence_path() == other.resolve()


def test_the_sequence_field_follows_browse_and_refuses_nonsense(qapp, isolated_settings,
                                                               two_files, tmp_path):
    knob, other = two_files
    page = _page(knob)
    page.sequence_row.set_path(other)                 # as Browse does
    assert page.sequence_path() == other.resolve()
    page.sequence_row.set_path(tmp_path / "missing.fa")
    assert page.sequence_path() == other.resolve()
    assert page.sequence_row.path() == other.resolve()
    page.sequence_row.clear()
    assert page.sequence_row.path() == other.resolve()


def test_primers_kept_here_are_saved_and_announced(qapp, isolated_settings, tmp_path):
    from repeatviewer.backends.primers import on_sequence as osq

    template = ("GGCATCGATCGATTACGGA" * 3 + "CCTGCTAGTCGTAGACCTTGAA" + _seq(5, 200)
                + "GGCATACGATTTGCAGCAAGT" + "TTGACCGGTA" * 4)
    fasta = _fasta(tmp_path / "t.fa", [("t", template)])
    page = _page(fasta)
    saved = []
    page.annotationsSaved.connect(saved.append)
    page.primers.add_primers([Primer("F", "CCTGCTAGTCGTAGACCTTGAA"),
                              Primer("R", "ACTTGCTGCAAATCGTATGCC")], tick=True)
    qapp.processEvents()
    assert osq.primers_sidecar_path(fasta).is_file()
    names = {osq.feature_primer_name(f) for f in page.panel.primers_host_features()
             if osq.is_primer_feature(f)}
    assert names == {"F", "R"} and saved
    assert not page.has_unsaved_changes() and page.confirm_leave()

    # another view of the same file reads them, and this one follows a save made there
    again = _page(fasta)
    names = {osq.feature_primer_name(f) for f in again.panel.primers_host_features()
             if osq.is_primer_feature(f)}
    assert names == {"F", "R"}
    again.primers.add_primers([Primer("F2", "GGCATCGATCGATTACGGA")], tick=True)
    qapp.processEvents()
    page.reload_annotations()
    names = {osq.feature_primer_name(f) for f in page.panel.primers_host_features()
             if osq.is_primer_feature(f)}
    assert "F2" in names
    from repeatviewer.widgets.sequence_annotation_rail import PRIMERS_TAB_INDEX

    assert page.panel._rail._active_tab == PRIMERS_TAB_INDEX


@pytest.fixture()
def make_host(qapp):
    """A stand-in main window with the primer menu actions; removed after the test.

    A window left behind would be found by ``find_main_window`` in later tests.
    """
    from PySide6.QtCore import QEvent

    from repeatviewer.job_options import load_primer_options
    from repeatviewer.pcr_actions import PcrActionsMixin
    from repeatviewer.primer_actions import PrimerActionsMixin

    made = []

    class Host(PrimerActionsMixin, PcrActionsMixin, QWidget):
        def __init__(self, tabs, selected):
            super().__init__()
            self._preview_tabs = tabs
            self._primer_options = load_primer_options()
            self.selected = list(selected)
            self.status = []

        def _selected_files(self):
            return list(self.selected)

        def _ensure_preview_visible(self):
            pass

        def _set_status_message(self, text, *a, **k):
            self.status.append(text)

    def build(tabs, selected):
        host = Host(tabs, selected)
        made.extend((host, tabs))
        return host

    yield build
    from PySide6.QtTest import QTest

    QTest.qWait(200)             # queued refreshes finish before their widgets go
    for widget in made:
        widget.close()
        widget.deleteLater()
    QCoreApplication.sendPostedEvents(None, QEvent.Type.DeferredDelete)


def test_every_menu_item_opens_the_one_workspace(qapp, isolated_settings, two_files, make_host):
    pytest.importorskip("PySide6.QtPdf")   # the preview tabs need it
    from repeatviewer.backends.primers import design as dz
    from repeatviewer.views import primers_workspace as pw
    from repeatviewer.widgets.preview_tabs import PreviewTabWidget

    knob, other = two_files
    tabs = PreviewTabWidget()
    host = make_host(tabs, [knob])
    tabs.open_file(knob)
    tabs.open_tool_tab(QLabel("one"), "Sanger Check")

    primers = host._design_in_panel(dz.MODE_PCR) or host._primers_workspace().primers
    page = tabs.currentWidget()
    assert isinstance(page, pw.PrimersWorkspacePage) and page is host._primers_workspace()
    assert page.sequence_path() == knob.resolve() and page.section() == pw.SECTION_DESIGN
    assert primers.design_pane.mode.currentData() == dz.MODE_PCR

    # a tool tab comes and goes, another item is clicked: still the same one workspace
    tabs.open_tool_tab(QLabel("two"), "Sanger Check")
    host._open_pcr_in_panel()
    assert tabs.currentWidget() is page and page.section() == pw.SECTION_PCR
    host._open_digest_in_panel()
    assert tabs.currentWidget() is page and page.section() == pw.SECTION_DIGEST
    host._design_in_panel(dz.MODE_FISH)
    assert page.section() == pw.SECTION_DESIGN
    assert page.primers.design_pane.mode.currentData() == dz.MODE_FISH
    count = sum(isinstance(tabs.widget(i), pw.PrimersWorkspacePage) for i in range(tabs.count()))
    assert count == 1
    # the file's own Annotation view never shows the Primers page
    from repeatviewer.widgets.sequence_annotation_rail import PRIMERS_TAB_INDEX

    file_page = tabs.page_for_path(knob)
    panel = getattr(file_page, "_annotation_panel", None)
    assert panel is None or panel._rail._active_tab != PRIMERS_TAB_INDEX

    # another file selected in Files, then the menu: the workspace moves to it
    host.selected = [other]
    host._design_in_panel(dz.MODE_PCR)
    assert tabs.currentWidget() is page and page.sequence_path() == other.resolve()
    # nothing selected: it stays on what it has
    host.selected = []
    host._open_pcr_in_panel()
    assert page.sequence_path() == other.resolve() and page.section() == pw.SECTION_PCR


def test_a_click_in_files_sets_the_sequence_only_when_the_workspace_is_in_front(
        qapp, isolated_settings, two_files, tmp_path, make_host):
    pytest.importorskip("PySide6.QtPdf")
    from repeatviewer.backends.primers import design as dz
    from repeatviewer.widgets.preview_tabs import PreviewTabWidget

    knob, other = two_files
    note = tmp_path / "notes.txt"
    note.write_text("hello\n")
    tabs = PreviewTabWidget()
    host = make_host(tabs, [knob])
    assert not host._primers_workspace_takes_file(other)          # no workspace yet
    host._design_in_panel(dz.MODE_PCR)
    page = host._primers_workspace()

    assert host._primers_workspace_takes_file(other)
    assert page.sequence_path() == other.resolve()
    assert not host._primers_workspace_takes_file(note)           # not a sequence: Preview opens it
    assert not host._primers_workspace_takes_file(tmp_path / "gone.fa")

    tabs.open_file(knob)                                          # a file tab in front
    assert tabs.currentWidget() is not page
    assert not host._primers_workspace_takes_file(knob)
    assert page.sequence_path() == other.resolve()


def test_closing_the_tab(qapp, isolated_settings, two_files, make_host):
    pytest.importorskip("PySide6.QtPdf")
    from repeatviewer.backends.primers import design as dz
    from repeatviewer.widgets.preview_tabs import PreviewTabWidget

    tabs = PreviewTabWidget()
    host = make_host(tabs, [two_files[0]])
    host._design_in_panel(dz.MODE_PCR)
    page = host._primers_workspace()
    assert tabs.confirm_leave_page(page)             # nothing unsaved: no question
    tabs._close_tab(tabs.indexOf(page))
    assert host._primers_workspace() is None
    host._open_pcr_in_panel()                        # and it opens again
    assert host._primers_workspace() is not None


def test_a_file_without_a_sequence_is_refused(qapp, isolated_settings, tmp_path, monkeypatch,
                                              make_host):
    pytest.importorskip("PySide6.QtPdf")
    from PySide6.QtWidgets import QMessageBox

    from repeatviewer.backends.primers import design as dz
    from repeatviewer.widgets.preview_tabs import PreviewTabWidget

    empty = tmp_path / "empty.fa"
    empty.write_text("")
    said = []
    monkeypatch.setattr(QMessageBox, "information",
                        staticmethod(lambda *a, **k: said.append(a[2])))
    tabs = PreviewTabWidget()
    host = make_host(tabs, [empty])
    assert host._design_in_panel(dz.MODE_PCR) is None
    assert host._primers_workspace() is None and said and "empty.fa" in said[0]


TEMPLATE = ("GGCATCGATCGATTACGGA" * 3 + "CCTGCTAGTCGTAGACCTTGAA" + _seq(5, 200)
            + "GGCATACGATTTGCAGCAAGT" + "TTGACCGGTA" * 4)


def _names(page, kind_label):
    from repeatviewer.widgets import primers_workspace_list as wl

    return [page.list.item(r, wl.COL_NAME).text() for r in range(page.list.rowCount())
            if page.list.item(r, wl.COL_KIND).text() == kind_label]


def _wait(app, pane, seconds=20):
    import time

    end = time.time() + seconds
    while pane._busy and time.time() < end:
        app.processEvents()
        time.sleep(0.02)
    app.processEvents()


def test_the_list_under_the_map_shows_the_primers(qapp, isolated_settings, tmp_path):
    from repeatviewer.views import primers_workspace as pw
    from repeatviewer.widgets import primers_workspace_list as wl

    fasta = _fasta(tmp_path / "t.fa", [("t", TEMPLATE)])
    page = _page(fasta)
    primers = page.primers
    # the two source tables are kept but never shown; their buttons sit under the map
    assert page._kept.isHidden()
    assert primers.table.parentWidget() is page._kept
    assert primers.design_pane.table.parentWidget() is page._kept
    for button in (primers.add_button, primers.database_button,
                   primers.design_pane.report_button, primers.design_pane.keep_button):
        assert page.bottom.isAncestorOf(button)
    assert primers.add_button.text() == "Add Primer"
    assert primers.database_button.text() == "From Database"
    assert page.splitter.indexOf(page.panel._map) == 0 and page.splitter.indexOf(page.bottom) == 1
    assert page.list.rowCount() == 0

    primers.add_primers([Primer("F", "CCTGCTAGTCGTAGACCTTGAA"),
                         Primer("R", "ACTTGCTGCAAATCGTATGCC")], tick=False)
    qapp.processEvents()
    assert _names(page, "on the sequence") == ["F", "R"] and page.list.counts() == (0, 2)

    # a tick in the list is a tick of the primer
    row = page.list.row_of(wl.KIND_PRIMER, 0)
    from PySide6.QtCore import Qt

    page.list.item(row, wl.COL_TICK).setCheckState(Qt.CheckState.Checked)
    qapp.processEvents()
    assert primers._ticked == {"F"}
    page.list.item(row, wl.COL_TICK).setCheckState(Qt.CheckState.Unchecked)
    qapp.processEvents()
    assert primers._ticked == set()

    # the primer just added is the selected one: its details show, with the
    # settings folded to one line
    assert page.list.currentRow() == page.list.row_of(wl.KIND_PRIMER, 1)
    assert primers._selected == "R" and primers.tabs.currentIndex() == pw.PANEL_PRIMER
    assert "Standard PCR" in page.settings_summary.text() and "Tm" in page.settings_summary.text()
    page.show_settings_button.click()
    assert primers.tabs.currentIndex() == pw.PANEL_DESIGN
    page.list.cellClicked.emit(page.list.row_of(wl.KIND_PRIMER, 1), wl.COL_NAME)   # same row again
    assert primers.tabs.currentIndex() == pw.PANEL_PRIMER

    # a click on another primer in the list
    page.show_settings_button.click()
    page.list.selectRow(page.list.row_of(wl.KIND_PRIMER, 0))
    assert primers._selected == "F" and primers.tabs.currentIndex() == pw.PANEL_PRIMER
    # a primer picked somewhere else (a click on the map) is followed
    page.show_settings_button.click()
    primers._select_row_of("R")
    qapp.processEvents()
    assert page.list.currentRow() == page.list.row_of(wl.KIND_PRIMER, 1)
    assert primers._selected == "R" and primers.tabs.currentIndex() == pw.PANEL_PRIMER
    primers._select_row_of("F")
    qapp.processEvents()
    assert page.list.currentRow() == page.list.row_of(wl.KIND_PRIMER, 0)
    assert primers._selected == "F" and primers.tabs.currentIndex() == pw.PANEL_PRIMER

    # in PCR a picked primer does not pull the workspace back to Design
    page.set_section("pcr")
    primers._select_row_of("R")
    qapp.processEvents()
    assert page.section() == pw.SECTION_PCR
    # back in Design the selected primer's details are there
    page.set_section("design")
    assert primers.tabs.currentIndex() in (pw.PANEL_DESIGN, pw.PANEL_PRIMER)


def test_candidates_join_the_list_and_keep_ticked_works(qapp, isolated_settings, tmp_path):
    from PySide6.QtCore import Qt

    from repeatviewer.views import primers_workspace as pw
    from repeatviewer.widgets import primers_workspace_list as wl

    fasta = _fasta(tmp_path / "long.fa", [("long", _seq(10, 3000))])
    page = _page(fasta)
    primers, design = page.primers, page.primers.design_pane
    primers.on_sequence_selection(1400, 1600)
    design.prod_min.setValue(250)
    design.prod_max.setValue(600)
    design.run_design()
    _wait(qapp, design)
    n = design.table.rowCount()
    assert n > 0 and page.list.counts() == (n, 0)
    assert _names(page, "candidate")[0] == "Pair 1"
    assert page.list.item(0, wl.COL_TICK).checkState() == Qt.CheckState.Checked   # best pair
    assert page.list.currentRow() == 0                                            # and shown
    assert "region 1,400 to 1,600" in page.settings_summary_text()

    if n > 1:
        page.list.item(1, wl.COL_TICK).setCheckState(Qt.CheckState.Checked)
        qapp.processEvents()
        assert len(design.ticked_pairs()) == 2
        page.list.item(1, wl.COL_TICK).setCheckState(Qt.CheckState.Unchecked)
        qapp.processEvents()
    assert len(design.ticked_pairs()) == 1

    design.keep_button.click()                       # the real button, now under the map
    qapp.processEvents()
    assert len(primers.primers()) == 2
    assert page.list.counts() == (n, 2) and len(_names(page, "on the sequence")) == 2

    # a click on a candidate goes back to the settings and drops the primer selection
    page.list.selectRow(page.list.row_of(wl.KIND_PRIMER, 0))
    assert primers.tabs.currentIndex() == pw.PANEL_PRIMER
    page.list.selectRow(0)
    assert primers.tabs.currentIndex() == pw.PANEL_DESIGN and primers._selected == ""

    # another sequence: the list empties with its sources
    other = _fasta(tmp_path / "o.fa", [("o", _seq(11, 800))])
    assert page.set_sequence(other)
    qapp.processEvents()
    assert page.list.counts() == (0, 0)


def test_a_path_field_that_had_the_cursor_takes_the_clicked_file(
        qapp, isolated_settings, two_files, tmp_path, make_host):
    pytest.importorskip("PySide6.QtPdf")
    from repeatviewer.backends.primers import design as dz
    from repeatviewer.widgets.preview_tabs import PreviewTabWidget

    knob, other = two_files
    families = _fasta(tmp_path / "families.fa", [("x", _seq(7, 300))])
    tabs = PreviewTabWidget()
    host = make_host(tabs, [knob])
    host._design_in_panel(dz.MODE_FISH)
    page = host._primers_workspace()
    tabs.resize(1300, 800)
    tabs.show()
    qapp.processEvents()
    field = page.primers.design_pane.avoid_path
    assert field.isVisible()                          # Repetitive Probes shows it

    assert not page.fill_armed_path(families)         # no field had the cursor
    page._on_focus_changed(None, page.sequence_row.line_edit())
    assert not page.fill_armed_path(families)         # the Sequence field is not one of them
    page._on_focus_changed(None, field.line_edit())   # the cursor goes into Avoid file
    page._on_focus_changed(field.line_edit(), tabs)   # then to Files: still armed
    assert host._primers_workspace_takes_file(families)
    assert field.path() == families.resolve()
    assert page.sequence_path() == knob.resolve()     # the sequence did not change
    # the next click is an ordinary one again
    assert host._primers_workspace_takes_file(other)
    assert page.sequence_path() == other.resolve()
