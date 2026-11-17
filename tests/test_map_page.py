"""Map as a centre tab: per-sequence table, binding picture, settings panel."""

import json
import random

import pytest
from PySide6.QtCore import QCoreApplication, QSettings

from repeatviewer.backends.primers import jobs as pj
from repeatviewer.backends.primers import thermo as th


@pytest.fixture()
def isolated_settings(qapp, tmp_path):
    org, app = QCoreApplication.organizationName(), QCoreApplication.applicationName()
    QCoreApplication.setOrganizationName("RVMapPage")
    QCoreApplication.setApplicationName(f"t{id(tmp_path)}")
    QSettings().clear()
    yield
    QSettings().clear()
    QCoreApplication.setOrganizationName(org)
    QCoreApplication.setApplicationName(app)


def _files(tmp_path):
    r = random.Random(4)
    mono = "".join(r.choice("ACGT") for _ in range(170))
    other = "".join(r.choice("ACGT") for _ in range(182))
    other = other[:44] + mono[10:32] + other[66:]          # only the forward primer binds
    (tmp_path / "CroSat3_monomer.fa").write_text(">CroSat3_monomer\n" + mono + "\n")
    (tmp_path / "CroSat1_monomer.fa").write_text(">CroSat1_monomer\n" + other + "\n")
    return ([tmp_path / "CroSat3_monomer.fa", tmp_path / "CroSat1_monomer.fa"],
            mono[10:32], th.reverse_complement(mono[120:142]))


def _run(qapp, page):
    import time

    page.run_map()
    end = time.time() + 40
    while page._busy and time.time() < end:
        qapp.processEvents()
        time.sleep(0.02)


def test_map_page_lists_every_sequence_and_draws_the_sites(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.map_page import MapPage

    files, fwd, rev = _files(tmp_path)
    page = MapPage(files, database_path="")
    page.source.setCurrentIndex(
        [page.source.itemData(i) for i in range(page.source.count())].index(pj.MAP_SPECIFIC))
    page.forward.setCurrentText(fwd)
    page.reverse.setCurrentText(rev)
    assert page.tab_name().startswith("Map · ")

    _run(qapp, page)
    assert page.table.rowCount() == 2
    first = [page.table.item(0, c).text() for c in range(6)]
    assert first[0] == "CroSat3_monomer" and first[1] == "11" and first[2] == "121"
    assert first[4] == "132 bp"
    second = [page.table.item(1, c).text() for c in range(6)]
    assert second[0] == "CroSat1_monomer" and second[2] == "no site" and second[4] == "none"
    assert "2 of 3" not in page.verdict.text() and "1 of 2" in page.verdict.text()

    # the picture knows both rows and paints
    assert [row["id"] for row in page.picture._rows] == ["CroSat3_monomer", "CroSat1_monomer"]
    page.picture.resize(700, 140)
    page.picture.grab()
    page.deleteLater()


def test_map_page_settings_and_saving(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.map_page import MapPage

    files, fwd, rev = _files(tmp_path)
    page = MapPage(files, database_path="")
    page.source.setCurrentIndex(
        [page.source.itemData(i) for i in range(page.source.count())].index(pj.MAP_SPECIFIC))
    page.forward.setCurrentText(fwd)
    page.mismatches.setValue(1)
    page.exact.setValue(0)
    page.size_min.setValue(50)
    page.size_max.setValue(900)
    page.pairs_only.setChecked(True)
    s = page.settings()
    assert s.mode == pj.MAP_SPECIFIC and s.forward == fwd and not s.use_reverse
    assert (s.max_mismatches, s.three_prime_exact) == (1, 0)
    assert (s.product_min, s.product_max) == (50, 900) and s.pairs_only

    page.save_button.click()
    saved = json.loads(QSettings().value("primers/map_page"))
    assert saved["max_mismatches"] == 1 and saved["product_max"] == 900
    page.defaults_button.click()
    assert page.mismatches.value() == pj.MapSettings().max_mismatches
    page.deleteLater()


def test_map_page_asks_for_files_and_primers(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.map_page import MapPage

    page = MapPage([], database_path="")
    page.run_map()
    assert "Tick the FASTA files" in page.verdict.text()
    files, _fwd, _rev = _files(tmp_path)
    page.set_files(files)
    page.run_map()
    assert "at least one primer" in page.verdict.text()
    assert "2 files" in page.files_note.text()
    page.deleteLater()


def test_report_is_written_only_when_asked(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.map_page import MapPage

    files, fwd, rev = _files(tmp_path)
    page = MapPage(files, database_path="")
    page.source.setCurrentIndex(
        [page.source.itemData(i) for i in range(page.source.count())].index(pj.MAP_SPECIFIC))
    page.forward.setCurrentText(fwd)
    page.reverse.setCurrentText(rev)
    _run(qapp, page)
    assert not list(tmp_path.glob("*primer_map*"))        # nothing written yet

    page.report_check.setChecked(True)
    _run(qapp, page)
    assert page.report_path.is_file() and "Report" in page.saved_note.text()
    page.deleteLater()
