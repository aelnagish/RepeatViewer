"""Sanger Check as a centre tab: table, coverage bar, chromatogram, settings panel."""

import random

import pytest
from PySide6.QtCore import QCoreApplication, QSettings

from repeatviewer.backends.primers import sanger as sg
from repeatviewer.backends.primers import thermo as th
from tests.test_sanger_dialog import _ab1, _insert


@pytest.fixture()
def isolated_settings(qapp, tmp_path):
    org, app = QCoreApplication.organizationName(), QCoreApplication.applicationName()
    QCoreApplication.setOrganizationName("RVSangerPage")
    QCoreApplication.setApplicationName(f"t{id(tmp_path)}")
    QSettings().clear()
    yield
    QSettings().clear()
    QCoreApplication.setOrganizationName(org)
    QCoreApplication.setApplicationName(app)


def _page(tmp_path, seq, reads):
    from repeatviewer.views.sanger_page import SangerPage

    page = SangerPage(seq, expected_name="pGEM insert", folder=tmp_path)
    page.add_reads(reads)
    return page


def test_page_shows_reads_differences_and_coverage(qapp, isolated_settings, tmp_path):
    seq = _insert(600, seed=21)
    read = list(seq[20:360])
    read[150] = "A" if read[150] != "A" else "G"
    fwd = _ab1(tmp_path / "cloneA_T7.ab1", "".join(read), [42] * len(read))
    rev = _ab1(tmp_path / "cloneA_SP6.ab1", th.reverse_complement(seq[250:]), [40] * 350)
    page = _page(tmp_path, seq, [fwd, rev])

    assert page.tab_name() == "Sanger Check · pGEM insert"
    assert "1 difference" in page.verdict.text()
    assert page.table.rowCount() == 2
    assert page.table.item(0, 0).text() == "cloneA_T7"
    assert page.table.item(0, 1).text() == "forward"
    assert page.table.item(0, 2).text() == "21 to 360"
    assert page.table.item(0, 4).text().startswith("171:")
    assert page.table.item(0, 5).text() == "Q42"
    assert page.table.item(1, 1).text() == "reverse" and page.table.item(1, 4).text() == "none"

    # the coverage bar knows the length, the reads and the gaps
    assert page.coverage._length == len(seq)
    assert [c.name for c in page.coverage._checks] == ["cloneA_T7", "cloneA_SP6"]
    assert sg.gaps_in_coverage(len(seq), page.coverage._checks)[0] == (1, 20)
    page.coverage.resize(800, 120)
    page.coverage.grab()                       # paints without raising

    # picking a difference moves the chromatogram to it
    page.table.selectRow(0)
    assert page.trace._centre == 150 and page.trace._expected == seq[170]
    assert "expected 171" in page.trace_title.text()
    page.deleteLater()


def test_quality_threshold_and_saving(qapp, isolated_settings, tmp_path):
    import json

    seq = _insert(400, seed=22)
    read = list(seq[10:300])
    read[100] = "A" if read[100] != "A" else "G"
    quality = [40] * len(read)
    quality[100] = 18                          # a poor call, not the clone
    poor = _ab1(tmp_path / "cloneB.ab1", "".join(read), quality)
    page = _page(tmp_path, seq, [poor])

    diff = page._rows[0][1]
    assert diff is not None and page._doubtful(diff)
    page.good_q.setValue(10)                   # trust calls from Q10 up
    page.recheck()
    assert not page._doubtful(page._rows[0][1])

    page.cutoff.setValue(0.02)
    page.report_check.setChecked(True)
    page.save_button.click()
    saved = json.loads(QSettings().value("primers/sanger_page"))
    assert saved["good_quality"] == 10 and saved["cutoff"] == 0.02 and saved["report"]
    page.defaults_button.click()
    assert page.good_q.value() == sg.GOOD_QUALITY and not page.report_check.isChecked()
    page.deleteLater()


def test_report_is_written_when_asked(qapp, isolated_settings, tmp_path):
    seq = _insert(300, seed=23)
    read = _ab1(tmp_path / "c.ab1", seq[:200], [40] * 200)
    page = _page(tmp_path, seq, [read])
    page.report_check.setChecked(True)
    page.check_button.click()
    out = tmp_path / "pGEM_insert_sanger.html"
    assert out.is_file() and "pGEM insert" in out.read_text()
    assert "Report" in page.saved_note.text()
    page.deleteLater()


def test_empty_page_asks_for_input(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.sanger_page import SangerPage

    page = SangerPage(folder=tmp_path)
    assert "Click an .ab1 read" in page.verdict.text()
    assert page.table.rowCount() == 0 and page.coverage._length == 0
    page.set_expected(_insert(200, seed=24), "insert_1")
    assert "Add the .ab1 reads" in page.verdict.text()
    assert page.tab_name() == "Sanger Check · insert_1"
    page.deleteLater()
