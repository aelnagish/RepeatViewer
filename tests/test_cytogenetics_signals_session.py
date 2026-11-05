"""Signals tab logic: numbering, signals, measuring, summary, karyotype helpers."""

from __future__ import annotations

import numpy as np
import pytest

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import image_session as isn
from repeatviewer.backends.cytogenetics import karyotype as kt
from repeatviewer.backends.cytogenetics import model as md
from repeatviewer.backends.cytogenetics import signals_session as ss
from repeatviewer.backends.cytogenetics.importer import ImportRequest, run_request

tifffile = pytest.importorskip("tifffile")


def _tiff(path, spot=None, shape=(80, 100)):
    arr = np.full(shape, 100, dtype=np.uint16)
    if spot:
        x, y = spot
        arr[y - 2:y + 3, x - 2:x + 3] = 3000
    tifffile.imwrite(str(path), arr)


def _folder(tmp_path, cells=2):
    src = tmp_path / "scope"
    src.mkdir()
    for n in range(1, cells + 1):
        _tiff(src / f"s1_c0{n}_DAPI.tif")
        _tiff(src / f"s1_c0{n}_FITC.tif", spot=(40, 30))
    work = tmp_path / "data"
    work.mkdir()
    run_request(ImportRequest(folder=work, source=src, slide_name="S1", probes_by_dye={"FITC": "CroSat2"}))
    k = kt.new_karyotype("C. sativus", 2, 2, species="Crocus sativus")
    k.chromosomes[0].centromere_index = 40
    kt.save_karyotype(work, k)
    return work


def _session(work, index=0):
    return ss.SignalsSession(work, ci.load_images(work)[index])


def test_picks_karyotype_and_labels(tmp_path):
    s = _session(_folder(tmp_path))
    assert s.karyotype.name == "C. sativus"
    assert s.labels() == ["1a", "1b", "2a", "2b"]
    assert s.next_label() == "1a"
    s.add_mark("1a", 10, 10)
    s.add_mark("1b", 30, 10)
    assert s.next_label() == "2a"
    with pytest.raises(ValueError):
        s.add_mark("1a", 5, 5)
    with pytest.raises(ValueError):
        s.add_mark("9z", 5, 5)
    s.relabel_mark(s.mark_for_label("1b").id, "2b")
    assert s.mark_for_label("2b") is not None and s.mark_for_label("1b") is None


def test_signal_refine_and_assignment(tmp_path):
    s = _session(_folder(tmp_path))
    a = s.add_mark("1a", 20, 60)
    b = s.add_mark("1b", 80, 60)
    sig = s.add_signal("CroSat2", 42, 31)
    assert (round(sig.x), round(sig.y)) == (40, 30)  # moved onto the bright spot
    assert sig.chromosome_id in (a.id, b.id)
    far = s.add_signal("CroSat2", 90, 70, refine=False)
    assert far.chromosome_id == b.id


def test_measure_classifies_signal(tmp_path):
    s = _session(_folder(tmp_path))
    mark = s.add_mark("1a", 40, 10)
    s.set_axis(mark.id, (40, 0), (40, 100))  # ci 40 -> centromere at y = 40
    assert mark.centromere == (40.0, 40.0)
    sig = s.add_signal("CroSat2", 40, 100, refine=False)
    assert sig.arm == "q" and sig.position == "terminal" and sig.rel_distance == 1.0
    near = s.add_signal("CroSat2", 40, 42, refine=False)
    assert near.arm == "cen" and near.position == "pericentromeric"
    s.update_signal(near.id, arm="p", position="interstitial")
    s.set_axis(mark.id, (40, 0), (40, 100))
    assert near.arm == "p" and near.manual
    assert "of arm length" in s.distance_text(sig)


def test_remove_mark_unassigns(tmp_path):
    s = _session(_folder(tmp_path))
    mark = s.add_mark("1a", 40, 10)
    sig = s.add_signal("CroSat2", 40, 12, refine=False)
    assert sig.chromosome_id == mark.id
    s.remove_mark(mark.id)
    assert sig.chromosome_id == ""
    s.remove_signal(sig.id)
    assert not s.details.signals


def test_hit(tmp_path):
    s = _session(_folder(tmp_path))
    s.add_mark("1a", 20, 20)
    sig = s.add_signal("CroSat2", 60, 60, refine=False)
    assert s.hit(61, 59) == ("signal", sig.id)
    assert s.hit(21, 21)[0] == "mark"
    assert s.hit(95, 5) is None


def test_save_keeps_display_changes(tmp_path):
    work = _folder(tmp_path)
    marks = _session(work)
    marks.add_mark("1a", 10, 10)
    shown = isn.ImageSession(work, ci.load_images(work)[0])
    shown.set_levels(0, 10, 500)
    shown.save()
    marks.save()  # started before the display change
    again = md.image_details(ci.load_images(work)[0])
    assert again.chromosomes and (again.channels[0].display.low, again.channels[0].display.high) == (10, 500)
    shown.set_gamma(0, 1.4)
    shown.save()  # saved after the marks
    assert md.image_details(ci.load_images(work)[0]).chromosomes


def test_slide_cells_and_summary(tmp_path):
    work = _folder(tmp_path)
    for index in (0, 1):
        s = _session(work, index)
        mark = s.add_mark("1a", 40, 10)
        s.set_axis(mark.id, (40, 0), (40, 100))
        s.add_signal("CroSat2", 40, 90, refine=False)
        s.save()
    s = _session(work, 0)
    cells = s.slide_cells()
    assert [c.label for c in cells] == ["01", "02"] and all(c.marked for c in cells)
    rows, marked = s.summary()
    assert marked == 2 and rows[0].chromosome == "1" and (rows[0].seen, rows[0].total) == (2, 2)
    assert s.unsignalled() == []


def test_arm_text_and_bands():
    assert kt.fraction_to_arm_text(0.92, 0.4) == "q 0.87"
    f = kt.arm_text_to_fraction("q 0.87", 0.4)
    assert f == pytest.approx(0.922)
    assert kt.arm_text_to_fraction("p 0.5", 0.4) == pytest.approx(0.2)
    assert kt.arm_text_to_fraction("x 1", 0.4) is None
    c = kt.KaryoChromosome("1", centromere_index=40, bands=[kt.Band("Secondary constriction (NOR)", 0.9, 1.0)])
    assert kt.band_summary(c) == "NOR (q)" and kt.band_summary(kt.KaryoChromosome("2")) == "none"


def test_auto_axis_finds_telomeres_and_centromere(tmp_path):
    src = tmp_path / "scope2"
    src.mkdir()
    dapi = np.full((120, 120), 50, dtype=np.uint16)
    dapi[20:100, 55:65] = 3000          # vertical chromosome, y 20..99
    dapi[44:48, 55:58] = 50              # constriction at y ~46 (p arm is the short top part)
    dapi[44:48, 62:65] = 50
    tifffile.imwrite(str(src / "c01_DAPI.tif"), dapi)
    tifffile.imwrite(str(src / "c01_FITC.tif"), np.full((120, 120), 100, dtype=np.uint16))
    work = tmp_path / "w2"
    work.mkdir()
    run_request(ImportRequest(folder=work, source=src, slide_name="S"))
    kt.save_karyotype(work, kt.new_karyotype("K", 1, 2))
    s = ss.SignalsSession(work, ci.load_images(work)[0])
    mark = s.add_mark("1a", 60, 70)
    assert s.auto_axis(mark.id)
    assert mark.centromere is not None and abs(mark.centromere[1] - 46) < 6
    assert mark.telomere_p[1] < mark.telomere_q[1]  # p (short arm) is at the top
    assert s.estimate_axis(5, 5) is None
    s.set_point(mark.id, "centromere", (60, 60))
    assert mark.centromere == (60.0, 60.0)
    s.swap_arms(mark.id)
    assert mark.telomere_p[1] > mark.telomere_q[1]


def test_probes_include_lists_and_picture_snap(tmp_path):
    work = _folder(tmp_path)
    from repeatviewer.backends.cytogenetics import lists as ls
    lists = ls.load_lists(work)
    ls.add_entry(lists, "probes", "45S rDNA", dye="Cy3")
    ls.save_lists(work, lists)
    s = _session(work)
    assert ("CroSat2", "FITC") in s.probes() and ("45S rDNA", "Cy3") in s.probes()
    s2 = _session(work)
    s2.details.channels = []
    pic = np.zeros((80, 100, 3), dtype=np.uint8)
    pic[28:33, 58:63] = (255, 30, 40)  # red spot
    s2.set_picture(pic)
    x, y = s2.refine_position("45S rDNA", 57, 27)
    assert (round(x), round(y)) == (60, 30)
