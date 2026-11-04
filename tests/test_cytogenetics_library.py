"""Library logic, the lists editor and the slide-session import request."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import library as lb
from repeatviewer.backends.cytogenetics import lists as ls
from repeatviewer.backends.cytogenetics import model as md
from repeatviewer.backends.cytogenetics.importer import ImportRequest, run_request, scan_tiffs, targets_for_probes

tifffile = pytest.importorskip("tifffile")


def _image(i, *, slide="", probes=(), species="", stage="", signals=0, seq=""):
    img = ci.CytoImage(id=i, file=f"{i}.png", original_name=f"cell {i}",
                       targets=[ci.ImageTarget(seq, probe=probes[0])] if seq and probes else [],
                       metadata={k: v for k, v in (("species", species), ("stage", stage)) if v})
    md.set_image_details(img, md.ImageDetails(
        slide_id=slide, channels=[md.Channel(id="c" + i, file="c.tif", probe=p) for p in probes],
        signals=[md.Signal(id=f"s{n}", probe="CroSat2", x=0, y=0) for n in range(signals)]))
    return img


def _folder(tmp_path):
    s1 = md.add_slide(tmp_path, "Slide 12", date="2026-10-14", species="Crocus sativus")
    s2 = md.add_slide(tmp_path, "Slide 9")
    ci.save_images(tmp_path, [
        _image("1", slide=s1.id, probes=("CroSat2", "45S"), species="Crocus sativus", stage="Mitotic metaphase", signals=2, seq="CroSat2.fasta"),
        _image("2", slide=s1.id, probes=("CroSat2",), species="Crocus sativus", stage="Interphase"),
        _image("10", slide=s2.id, probes=("CroSat1",), species="C. cartwrightianus", stage="Mitotic metaphase"),
        _image("3"),
    ])
    return lb.build_items(tmp_path)


def test_filter_and_facets(tmp_path):
    items = _folder(tmp_path)
    assert len(items) == 4 and items[0].marked and items[0].channels == 2
    assert [i.image.id for i in lb.apply_filter(items, lb.LibraryFilter(text="interphase"))] == ["2"]
    assert [i.image.id for i in lb.apply_filter(items, lb.LibraryFilter(text="slide 12 45s"))] == ["1"]
    assert {i.image.id for i in lb.apply_filter(items, lb.LibraryFilter(probe={"CroSat2"}))} == {"1", "2"}
    assert [i.image.id for i in lb.apply_filter(items, lb.LibraryFilter(sequence="CroSat2.fasta"))] == ["1"]
    assert [i.image.id for i in lb.apply_filter(items, lb.LibraryFilter(signals={lb.SIGNALS_MARKED}))] == ["1"]
    flt = lb.LibraryFilter(species={"Crocus sativus"})
    counts = dict(lb.facet_counts(items, flt, "probe"))
    assert counts == {"45S": 1, "CroSat1": 0, "CroSat2": 2}
    assert dict(lb.facet_counts(items, flt, "species")) == {"C. cartwrightianus": 1, "Crocus sativus": 2}
    assert lb.facet_counts(items, lb.LibraryFilter(), "signals") == [(lb.SIGNALS_MARKED, 1), (lb.SIGNALS_NOT_MARKED, 3)]
    assert lb.LibraryFilter().is_empty() and not flt.is_empty()


def test_grouping(tmp_path):
    items = _folder(tmp_path)
    groups = lb.group_items(items, "Slide")
    assert [g.title for g in groups] == ["Slide 9", "Slide 12", lb.NO_SLIDE]
    assert groups[1].subtitle == "2026-10-14 · Crocus sativus · CroSat2, 45S · 2 images"
    by_probe = {g.title: [i.image.id for i in g.items] for g in lb.group_items(items, "Probe")}
    assert by_probe["CroSat2"] == ["1", "2"] and by_probe["No probe"] == ["3"]
    assert [g.title for g in lb.group_items(items, "Probe")][-1] == "No probe"
    assert [i.image.id for i in lb.group_items(items, "Nothing")[0].items] == ["1", "2", "3", "10"]
    assert lb.group_items([], "Slide") == []


def test_lists_editor(tmp_path):
    ci.save_images(tmp_path, [_image("1", probes=("CroSat2",), species="C. sat", stage="Interphase")])
    ed = ls.ListsEditor(tmp_path)
    ed.add("probes", "CroSat2", sequence="CroSat2.fasta")
    ed.add("species", "C. sat")
    with pytest.raises(ValueError):
        ed.add("species", "c. sat")
    ed.rename("species", "C. sat", "Crocus sativus")
    with pytest.raises(ValueError):
        ed.rename("species", "Crocus sativus", "")
    assert ed.used_by("species", "Crocus sativus") == 1 and not ed.can_remove("species", "Crocus sativus")
    with pytest.raises(ValueError):
        ed.remove("species", "Crocus sativus")
    ed.add("species", "Crocus cartwrightianus")
    ed.remove("species", "Crocus cartwrightianus")
    ed.add("probes", "CroSat2 old")
    ed.merge("probes", ["CroSat2 old"], "CroSat2")
    ed.set_value("probes", "CroSat2", "sequence", "")
    assert ls.entry(ed.lists, "probes", "CroSat2").get("sequence") is None
    # nothing is written before apply
    assert ci.load_images(tmp_path)[0].metadata["species"] == "C. sat"
    assert ed.apply() == 1
    assert ci.load_images(tmp_path)[0].metadata["species"] == "Crocus sativus"
    assert ls.names(ls.load_lists(tmp_path), "species") == ["Crocus sativus"]


def _tif(path, value=1000, shape=(30, 40)):
    arr = np.zeros(shape, dtype=np.uint16)
    arr[5:15, 5:20] = value
    tifffile.imwrite(path, arr, imagej=True, resolution=(1 / 0.0645, 1 / 0.0645), metadata={"unit": "um"})


def test_run_request(tmp_path):
    data = tmp_path / "data"
    data.mkdir()
    src = tmp_path / "scope"
    (src / "sub").mkdir(parents=True)
    for c in (1, 2):
        for dye in ("DAPI", "FITC"):
            _tif(src / f"s12_c0{c}_{dye}.tif")
    _tif(src / "sub" / "x_Cy3.tif")
    _tif(src / "mystery.tif")
    assert [p.name for p in scan_tiffs(src)] == ["mystery.tif", "s12_c01_DAPI.tif", "s12_c01_FITC.tif", "s12_c02_DAPI.tif", "s12_c02_FITC.tif"]
    assert len(scan_tiffs(src, recursive=True)) == 6
    lists = ls.load_lists(data)
    ls.add_entry(lists, "probes", "CroSat2", sequence="CroSat2.fasta")
    ls.save_lists(data, lists)
    req = ImportRequest(folder=data, source=src, slide_name="Slide 12", slide_values={"species": "Crocus sativus", "date": "2026-10-14"},
                        stage="Mitotic metaphase", probes_by_dye={"FITC": "CroSat2"})
    slide, images, errors = run_request(req)
    assert slide.name == "Slide 12" and len(images) == 2
    assert errors == ["mystery.tif: no known dye in the file name"]
    first = ci.load_images(data)[0]
    assert first.metadata == {"species": "Crocus sativus", "date": "2026-10-14", "stage": "Mitotic metaphase"}
    assert [(t.sequence, t.probe, t.dye) for t in first.targets] == [("CroSat2.fasta", "CroSat2", "FITC")]
    assert ci.images_for_sequence(data / "CroSat2.fasta") == ci.load_images(data)
    assert "Crocus sativus" in ls.names(ls.load_lists(data), "species")
    assert len(md.load_slides(data)) == 1
    empty = tmp_path / "empty"
    empty.mkdir()
    s2, imgs2, errs2 = run_request(ImportRequest(folder=data, source=empty))
    assert s2 is None and imgs2 == [] and errs2 == ["No channel TIFFs found."]
    assert len(md.load_slides(data)) == 1
    assert targets_for_probes(ls.load_lists(data), {"FITC": "Unknown"}) == []


def test_details_rows(tmp_path):
    items = _folder(tmp_path)
    rows = dict(lb.details_rows(items[0]))
    assert rows["Slide"] == "Slide 12" and rows["Channels"] == "channel (CroSat2) · channel (45S)"
    assert rows["Signals"] == "2 on 0 chromosomes" and rows["Sequences"] == "CroSat2.fasta"
    assert "Pixel size" not in rows and "Notes" not in rows


def test_probe_info_and_texts(tmp_path):
    items = _folder(tmp_path)
    first = items[0]
    assert [row[0] for row in first.probe_info] == ["CroSat2", "45S"]
    assert first.probe_info[0][3] == "CroSat2.fasta"
    assert lb.signals_text(first) == "0 chromosomes numbered, 2 signals marked"
    assert lb.signals_text(items[3]) == "Not marked yet"
    assert lb.probe_colours(items) == {"CroSat2": "#ffffff", "45S": "#ffffff", "CroSat1": "#ffffff"}
