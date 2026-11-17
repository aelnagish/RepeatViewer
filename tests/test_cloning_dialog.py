"""Cloning window: vector presets, warnings, colony sizes, gel and report."""

import random

import pytest
from PySide6.QtCore import QCoreApplication, QSettings

from repeatviewer.backends.primers import cloning as cl


@pytest.fixture()
def isolated_settings(qapp, tmp_path):
    org, app = QCoreApplication.organizationName(), QCoreApplication.applicationName()
    QCoreApplication.setOrganizationName("RVCloningTests")
    QCoreApplication.setApplicationName(f"t{id(tmp_path)}")
    QSettings().clear()
    yield
    QSettings().clear()
    QCoreApplication.setOrganizationName(org)
    QCoreApplication.setApplicationName(app)


def _insert(n=900, seed=12):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def test_vector_decides_the_methods_and_sizes(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.cloning_dialog import CloningDialog

    seq = _insert()
    dlg = CloningDialog(seq, insert_name="CroSat3 product", folder=tmp_path)
    assert dlg.vector.currentText() == "pGEM-T Easy"
    assert [dlg.method.itemData(i) for i in range(dlg.method.count())] == [cl.METHOD_TA]
    assert dlg.vector_length.value() == 3015 and not dlg.vector_length.isEnabled()
    assert not dlg.left.isVisible()                    # enzymes only for restriction

    dlg.polymerase.setCurrentIndex(
        [dlg.polymerase.itemData(i) for i in range(dlg.polymerase.count())].index("Taq"))
    dlg.empty_product.setValue(250)
    assert dlg.table.rowCount() == 2
    assert dlg.table.item(0, 2).text() == f"{3015 + len(seq):,} bp"
    assert dlg.table.item(0, 3).text() == f"{250 + len(seq):,} bp"
    assert "colony PCR" in dlg.verdict.text()
    assert "M13 forward" in dlg.notes.text() and "575 bp" in dlg.notes.text()

    dlg.polymerase.setCurrentIndex(
        [dlg.polymerase.itemData(i) for i in range(dlg.polymerase.count())].index("Phusion / Q5"))
    assert "blunt ends" in dlg.verdict.text()
    dlg.deleteLater()


def test_own_vector_restriction_and_one_orientation(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.cloning_dialog import CloningDialog

    clean = "".join(random.Random(4).choice("AT") for _ in range(500))
    dlg = CloningDialog(clean, insert_name="insert", folder=tmp_path)
    dlg.vector.setCurrentIndex(
        [dlg.vector.itemData(i) for i in range(dlg.vector.count())].index("Your own vector"))
    dlg.method.setCurrentIndex(
        [dlg.method.itemData(i) for i in range(dlg.method.count())].index(cl.METHOD_RE))
    assert dlg.left.isVisibleTo(dlg) and dlg.vector_length.isEnabled()
    assert "Choose the enzyme" in dlg.verdict.text()

    for box, enzyme in ((dlg.left, "EcoRI"), (dlg.right, "BamHI")):
        box.setCurrentIndex([box.itemData(i) for i in range(box.count())].index(enzyme))
    dlg.vector_length.setValue(2700)
    dlg.empty_product.setValue(180)
    assert dlg.table.rowCount() == 1                   # two ends, one way in
    assert dlg.table.item(0, 2).text() == f"{2700 + len(clean):,} bp"
    dlg.deleteLater()


def test_array_warning_gel_and_report(qapp, isolated_settings, tmp_path):
    from repeatviewer.views.cloning_dialog import CloningDialog

    r = random.Random(6)
    mono = "".join(r.choice("ACGT") for _ in range(170))
    dlg = CloningDialog(mono * 10, insert_name="CroSat3 array", folder=tmp_path)
    dlg.polymerase.setCurrentIndex(
        [dlg.polymerase.itemData(i) for i in range(dlg.polymerase.count())].index("Taq"))
    dlg.empty_product.setValue(250)
    assert "tandem array" in dlg.verdict.text()
    assert "30 °C" in dlg.verdict.text()

    dlg.gel_button.click()
    assert dlg.gel_path.is_file() and dlg.gel_path.suffix == ".rvgel"
    dlg.report_button.click()
    text = dlg.report_path.read_text()
    assert "CroSat3 array" in text and "tandem array" in text
    dlg.deleteLater()


def test_menu_offers_cloning(qapp):
    from PySide6.QtWidgets import QMenu

    from repeatviewer.pcr_actions import PcrActionsMixin
    from repeatviewer.primer_actions import PrimerActionsMixin

    class _Host(PrimerActionsMixin, PcrActionsMixin):
        pass

    menu = QMenu()
    _Host()._populate_primers_menu(menu)
    texts = [a.text() for a in menu.actions() if not a.isSeparator()]
    assert texts == ["Design", "Map", "Database", "In-silico PCR", "Restriction Digest",
                     "Gel", "Sanger Check", "Cloning"]
