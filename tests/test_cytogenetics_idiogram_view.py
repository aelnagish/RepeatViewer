"""Idiogram tab redesign: per-karyotype views, layers, measured signals that are not in the karyotype."""

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
pytest.importorskip("PySide6")

from PySide6.QtWidgets import QApplication  # noqa: E402

from repeatviewer.backends.cytogenetics import idiogram as idg  # noqa: E402
from repeatviewer.backends.cytogenetics import karyotype as kt  # noqa: E402
from test_cytogenetics_idiogram import _folder, _mark_cells  # noqa: E402


@pytest.fixture(scope="module")
def app():
    return QApplication.instance() or QApplication([])


def test_each_karyotype_keeps_its_view():
    s = idg.IdiogramSettings(karyotype="A", min_share=0.3, show_bands=False)
    s = idg.switch_karyotype(s, "B")
    assert s.karyotype == "B" and s.views["A"]["min_share"] == 0.3
    s.min_share = 0.9
    s = idg.switch_karyotype(s, "A")
    assert s.min_share == 0.3 and s.show_bands is False
    again = idg.IdiogramSettings.from_dict(s.to_dict())
    assert again.views["B"]["min_share"] == 0.9


def test_measured_added_to_karyotype(app, tmp_path):
    from repeatviewer.widgets.cytogenetics_idiogram import IdiogramTab

    work, _k = _folder(tmp_path)
    _mark_cells(work)
    t = IdiogramTab()
    t.refresh(work)
    t.show()
    QApplication.processEvents()
    assert len(t._unmatched_items) == 1  # CroSat2 at the q end of chromosome 1, no band there yet
    t._add_measured_to_karyotype()
    k = kt.load_karyotype(work, "Cs")
    assert any(kt.manual_probe(b.kind) == "CroSat2" for b in k.chromosomes[0].bands)
    assert t._unmatched_items == []
