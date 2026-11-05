"""Edit karyotype window: one homologue edited on its own, drawing a band, undo, B chromosome, species list."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")

from PySide6.QtCore import QPoint, Qt  # noqa: E402
from PySide6.QtTest import QTest  # noqa: E402
from PySide6.QtWidgets import QApplication  # noqa: E402

from repeatviewer.backends.cytogenetics import karyotype as kt  # noqa: E402
from repeatviewer.backends.cytogenetics.lists import load_lists, names  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def _dialog(tmp_path):
    from repeatviewer.views.cyto_karyotype_dialog import KaryotypeDialog

    k = kt.new_karyotype("test", 2, 3, "Hapazlem bazaza")
    k.chromosomes[0].length, k.chromosomes[0].centromere_index = 20, 37.5
    k.chromosomes[1].length, k.chromosomes[1].centromere_index = 17.5, 26
    kt.save_karyotype(tmp_path, k)
    d = KaryotypeDialog(tmp_path, "test")
    d.resize(1220, 720)
    d.show()
    QApplication.processEvents()
    return d


def test_homologue_draw_undo_and_save(app, tmp_path):
    d = _dialog(tmp_path)
    cv = d._canvas
    ci, hi, r, _view, _m = cv._items[1]  # chromosome 1, homologue b
    assert (ci, hi) == (0, 1)
    cv._select(0, 1, -1)
    assert d._same.isVisible() and d._same.isChecked()
    d._same.setChecked(False)
    d._len_edit.setText("21")
    d._on_len_field()
    c = d._k.chromosomes[0]
    assert kt.homologue_view(c, 1).length == 21 and c.length == 20
    d._set_tool("draw")
    ci, hi, r, _view, _m = next(it for it in cv._items if it[:2] == (0, 1))
    x = int(r.center().x())
    QTest.mousePress(cv, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, QPoint(x, int(r.top() + 0.6 * r.height())))
    QTest.mouseMove(cv, QPoint(x, int(r.top() + 0.7 * r.height())))
    QTest.mouseRelease(cv, Qt.MouseButton.LeftButton, Qt.KeyboardModifier.NoModifier, QPoint(x, int(r.top() + 0.7 * r.height())))
    c = d._k.chromosomes[0]
    assert len(kt.homologue_view(c, 1).bands) == 1 and c.bands == []
    d._undo_step()
    assert kt.homologue_view(d._k.chromosomes[0], 1).bands == []
    d._redo_step()
    assert len(kt.homologue_view(d._k.chromosomes[0], 1).bands) == 1
    d._duplicate_b()
    assert d._formula.text() == "2n = 3x + 1B = 7"
    d._save()
    saved = kt.load_karyotype(tmp_path, "test")
    assert saved.chromosomes[-1].b_chromosome
    assert kt.homologue_view(saved.chromosomes[0], 1).length == 21
    assert "Hapazlem bazaza" in names(load_lists(tmp_path), "species")


def test_band_edge_without_selection_and_species(app, tmp_path, monkeypatch):
    from PySide6.QtCore import QPointF
    from PySide6.QtWidgets import QInputDialog

    d = _dialog(tmp_path)
    c = d._k.chromosomes[0]
    c.bands.append(kt.Band("CMA+ band", 0.6, 0.8))
    d._canvas.refresh()
    ci, hi, r, _v, _m = d._canvas._items[0]
    what, *_rest, bi = d._canvas._hit(QPointF(r.center().x(), r.top() + 0.8 * r.height()))
    assert (what, bi) == ("edge1", 0)
    monkeypatch.setattr(QInputDialog, "getText", lambda *a, **k: ("Crocus sativus", True))
    d._on_species(d._species.findData("Add New Species…"))
    assert d._k.species == "Crocus sativus"
    assert "Crocus sativus" in names(load_lists(tmp_path), "species")
