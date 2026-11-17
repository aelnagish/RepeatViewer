"""Cloning as a centre tab: clone table, warnings, colony-PCR gel, settings panel."""

import json
import random

import pytest
from PySide6.QtCore import QCoreApplication, QSettings

from repeatviewer.backends.primers import cloning as cl


@pytest.fixture()
def isolated_settings(qapp, tmp_path):
    org, app = QCoreApplication.organizationName(), QCoreApplication.applicationName()
    QCoreApplication.setOrganizationName("RVCloningPage")
    QCoreApplication.setApplicationName(f"t{id(tmp_path)}")
    QSettings().clear()
    yield
    QSettings().clear()
    QCoreApplication.setOrganizationName(org)
    QCoreApplication.setApplicationName(app)


def _seq(n=900, seed=31):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def test_page_plans_both_orientations_and_draws_the_gel(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.cloning_page import CloningPage

    seq = _seq()
    page = CloningPage(seq, insert_name="CroSat3 product", folder=tmp_path)
    assert page.tab_name() == "Cloning · CroSat3 product"
    assert page.vector.currentText() == "pGEM-T Easy"
    assert not page.enzyme_row.isVisibleTo(page)          # TA needs no enzymes
    assert page.vector_length.value() == 3015 and not page.vector_length.isEnabled()

    names = [page.polymerase.itemData(i) for i in range(page.polymerase.count())]
    page.polymerase.setCurrentIndex(names.index("Taq"))
    page.empty_product.setValue(250)
    assert page.table.rowCount() == 2
    assert page.table.item(0, 2).text() == f"{3015 + len(seq):,} bp"
    assert page.table.item(0, 3).text() == f"{250 + len(seq):,} bp"
    assert page.gel.doc() is not None and page.gel_button.isEnabled()
    assert "M13 forward" in page.notes.text()

    # a proofreading polymerase cannot go into a TA vector
    page.polymerase.setCurrentIndex(names.index("Phusion / Q5"))
    assert "blunt ends" in page.verdict.text()
    page.deleteLater()


def test_page_warns_about_a_tandem_array_and_writes_the_files(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.cloning_page import CloningPage

    r = random.Random(33)
    mono = "".join(r.choice("ACGT") for _ in range(170))
    page = CloningPage(mono * 10, insert_name="CroSat3 array", folder=tmp_path)
    names = [page.polymerase.itemData(i) for i in range(page.polymerase.count())]
    page.polymerase.setCurrentIndex(names.index("Taq"))
    page.empty_product.setValue(250)
    assert "tandem array" in page.verdict.text() and "recA" in page.verdict.text()

    page.gel_button.click()
    assert page.gel_path.is_file() and page.gel_path.suffix == ".rvgel"
    page.report_button.click()
    assert "CroSat3 array" in page.report_path.read_text()
    assert "Report" in page.saved_note.text()
    page.deleteLater()


def test_restriction_shows_the_enzymes_and_one_orientation(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.cloning_page import CloningPage

    clean = "".join(random.Random(34).choice("AT") for _ in range(500))
    page = CloningPage(clean, insert_name="insert", folder=tmp_path)
    vectors = [page.vector.itemData(i) for i in range(page.vector.count())]
    page.vector.setCurrentIndex(vectors.index("Your own vector"))
    methods = [page.method.itemData(i) for i in range(page.method.count())]
    page.method.setCurrentIndex(methods.index(cl.METHOD_RE))
    assert page.enzyme_row.isVisibleTo(page) and page.vector_length.isEnabled()
    assert "Choose the enzyme" in page.verdict.text()

    for box, enzyme in ((page.left, "EcoRI"), (page.right, "BamHI")):
        box.setCurrentIndex([box.itemData(i) for i in range(box.count())].index(enzyme))
    page.vector_length.setValue(2700)
    assert page.table.rowCount() == 1                     # two ends, one way in
    saved = json.loads(QSettings().value("primers/cloning"))
    assert saved["vector"] == "Your own vector" and saved["method"] == cl.METHOD_RE
    page.deleteLater()
