"""Cytogenetics backend: model v2, karyotypes, lists, channels, import, signals."""

from __future__ import annotations

import json
import struct
import zlib
from pathlib import Path

import numpy as np
import pytest

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import channels as ch
from repeatviewer.backends.cytogenetics import karyotype as ky
from repeatviewer.backends.cytogenetics import lists as ls
from repeatviewer.backends.cytogenetics import model as md
from repeatviewer.backends.cytogenetics import signals as sg
from repeatviewer.backends.cytogenetics.importer import group_channel_files, import_cells, render_preview

tifffile = pytest.importorskip("tifffile")


def _tif(path: Path, value: int, *, pixel_um: float | None = 0.0645, shape=(40, 60)) -> Path:
    arr = np.zeros(shape, dtype=np.uint16)
    arr[10:20, 10:30] = value
    if pixel_um:
        tifffile.imwrite(path, arr, imagej=True, resolution=(1 / pixel_um, 1 / pixel_um), metadata={"unit": "um"})
    else:
        tifffile.imwrite(path, arr)
    return path


# --- store compatibility ---------------------------------------------------

def test_v1_index_loads_and_extra_fields_survive(tmp_path):
    store = ci.store_dir(tmp_path)
    store.mkdir()
    (store / "index.json").write_text(json.dumps({"version": 1, "images": [
        {"id": "a", "file": "a.png", "targets": [{"sequence": "CroSat2.fasta"}], "metadata": {"stage": "Interphase"}}]}))
    (img,) = ci.load_images(tmp_path)
    assert md.image_details(img).channels == []
    details = md.image_details(img)
    details.signals.append(md.Signal(id="s1", probe="CroSat2", x=1, y=2))
    md.set_image_details(img, details)
    ci.save_images(tmp_path, [img])
    (again,) = ci.load_images(tmp_path)
    assert md.image_details(again).signals[0].probe == "CroSat2"
    assert json.loads((store / "index.json").read_text())["version"] == 2
    ci.rename_sequence(tmp_path, "CroSat2.fasta", "CroSat2b.fasta")
    (renamed,) = ci.load_images(tmp_path)
    assert md.image_details(renamed).signals and renamed.targets[0].sequence == "CroSat2b.fasta"


def test_remove_image_deletes_copied_channels_only(tmp_path):
    store = ci.store_dir(tmp_path)
    store.mkdir()
    outside = _tif(tmp_path / "keep_DAPI.tif", 5)
    (store / "x.png").write_bytes(b"x")
    (store / "x_c1.tif").write_bytes(b"x")
    img = ci.CytoImage(id="x", file="x.png", extra={"channels": [
        {"id": "c1", "file": "x_c1.tif"}, {"id": "c2", "file": str(outside)}]})
    ci.save_images(tmp_path, [img])
    assert ci.remove_image(tmp_path, "x")
    assert not (store / "x_c1.tif").exists() and not (store / "x.png").exists()
    assert outside.exists()


# --- karyotype -------------------------------------------------------------

def test_karyotype_formula_labels_and_order(tmp_path):
    k = ky.new_karyotype("C. sativus", 8, 3, species="Crocus sativus")
    assert k.formula == "2n = 3x = 24"
    assert k.homologue_labels(k.chromosomes[1]) == ["2a", "2b", "2c"]
    k.chromosomes[0].copies = 4
    assert k.formula == "2n = 3x + 1 = 25"
    for c, (length, ci_) in zip(k.chromosomes, [(5, 40), (9, 30), (7, 48)]):
        c.length, c.centromere_index = length, ci_
    assert [c.name for c in k.ordered("length")][:3] == ["2", "3", "1"]
    assert [c.name for c in k.ordered("ci") if c.name in ("1", "2", "3")] == ["3", "1", "2"]
    k.set_basic_number(10)
    assert [c.name for c in k.chromosomes][-2:] == ["9", "10"]
    k.set_basic_number(8)
    path = ky.save_karyotype(tmp_path, k)
    assert path.name == "C._sativus.karyotype.json"
    k2 = ky.load_karyotype(tmp_path, "C. sativus")
    assert k2.formula == k.formula and k2.chromosomes[0].copies == 4
    k2.name = "Crocus sativus"
    ky.save_karyotype(tmp_path, k2, old_name="C. sativus")
    assert [x.name for x in ky.list_karyotypes(tmp_path)] == ["Crocus sativus"]


@pytest.mark.parametrize("ci_, kind", [(48, "m"), (37.5, "m"), (30, "sm"), (20, "st"), (5, "t"), (70, "sm")])
def test_levan(ci_, kind):
    assert ky.levan_type(ci_) == kind


def test_measurement_table_averages_homologues():
    text = "Chromosome\tp\tq\n1a\t4\t6\n1b\t4\t8\n2a\t3\t3\n"
    rows = ky.parse_measurement_table(text)
    assert rows[0][0] == "1" and rows[0][1] == pytest.approx(11.0)
    assert rows[0][2] == pytest.approx((40 + 100 * 4 / 12) / 2)
    k = ky.new_karyotype("k", 2, 2)
    assert ky.apply_measurements(k, rows + [("9", 1.0, None)]) == ["9"]
    assert k.chromosomes[1].centromere_index == 50.0
    csv_rows = ky.parse_measurement_table("name,length (µm),CI\n1,10,0.4\n")
    assert csv_rows == [("1", 10.0, 40.0)]
    with pytest.raises(ValueError):
        ky.parse_measurement_table("a\tb\n1\t2\n")


def test_assembly_links_and_lengths():
    k = ky.new_karyotype("k", 3, 2)
    assert ky.guess_assembly_links(k, ["Csat_chr01", "Csat_chr02", "scaffold_9"]) == 2
    assert ky.apply_assembly_lengths(k, {"Csat_chr01": 2_500_000, "Csat_chr02": 1_000_000}) == ["3"]
    assert k.chromosomes[0].length == 2.5 and k.unit == "Mb"


# --- lists -----------------------------------------------------------------

def test_lists_defaults_rename_merge_and_dyes(tmp_path):
    lists = ls.load_lists(tmp_path)
    assert "FITC" in ls.names(lists, "dyes") and "Interphase" in ls.names(lists, "stages")
    ls.add_entry(lists, "probes", "CroSat2")
    ls.add_entry(lists, "probes", "crosat2")
    assert ls.names(lists, "probes") == ["CroSat2"]
    img = ci.CytoImage(id="a", file="a.png", targets=[ci.ImageTarget("x.fasta", probe="CroSat2", dye="FITC")],
                       metadata={"species": "C. sat"})
    md.set_image_details(img, md.ImageDetails(
        channels=[md.Channel(id="c", file="c.tif", dye="FITC", probe="CroSat2")],
        signals=[md.Signal(id="s", probe="CroSat2", x=0, y=0)]))
    ci.save_images(tmp_path, [img])
    md.save_slides(tmp_path, [md.Slide(id="s1", name="Slide 1", species="C. sat")])
    ls.add_entry(lists, "species", "C. sat")
    assert ls.rename_entry(tmp_path, lists, "probes", "CroSat2", "CroSat2_v2") == 1
    (img,) = ci.load_images(tmp_path)
    d = md.image_details(img)
    assert img.targets[0].probe == d.channels[0].probe == d.signals[0].probe == "CroSat2_v2"
    ls.add_entry(lists, "species", "Crocus sativus")
    ls.merge_entries(tmp_path, lists, "species", ["C. sat"], "Crocus sativus")
    assert ci.load_images(tmp_path)[0].metadata["species"] == "Crocus sativus"
    assert md.load_slides(tmp_path)[0].species == "Crocus sativus"
    assert ls.names(lists, "species") == ["Crocus sativus"]
    with pytest.raises(ValueError):
        ls.remove_entry(tmp_path, lists, "species", "Crocus sativus")
    ls.save_lists(tmp_path, lists)
    assert ls.names(ls.load_lists(tmp_path), "probes") == ["CroSat2_v2"]
    assert ls.dye_for_filename("s12_c03_DAPI.tif", lists) == ("DAPI", 2)
    assert ls.dye_for_filename("cell-4 488.TIFF", lists) == ("FITC", 2)
    assert ls.dye_for_filename("notes.tif", lists) == ("", -1)


# --- channels --------------------------------------------------------------

def test_read_channel_pixel_size(tmp_path):
    d = ch.read_channel(_tif(tmp_path / "a.tif", 1000))
    assert d.bits == 16 and d.pixel_size_um == pytest.approx(0.0645, rel=1e-3)
    assert ch.read_channel(_tif(tmp_path / "b.tif", 1000, pixel_um=None)).pixel_size_um is None
    arr = np.zeros((8, 8), dtype=np.uint16)
    tifffile.imwrite(tmp_path / "c.tif", arr, resolution=(1e4 / 0.1, 1e4 / 0.1), resolutionunit="CENTIMETER")
    assert ch.read_channel(tmp_path / "c.tif").pixel_size_um == pytest.approx(0.1, rel=1e-3)


def test_composite_colours_levels_and_log():
    a = np.zeros((4, 4), dtype=np.uint16)
    a[0, 0] = 1000
    blue = md.ChannelDisplay(color="#0000ff", low=0, high=1000)
    green = md.ChannelDisplay(color="#00ff00", low=0, high=500)
    rgb = ch.composite([(a, blue), (a, green)])
    assert rgb.dtype == np.uint8 and tuple(rgb[0, 0]) == (0, 255, 255) and tuple(rgb[1, 1]) == (0, 0, 0)
    green.visible = False
    assert tuple(ch.composite([(a, blue), (a, green)])[0, 0]) == (0, 0, 255)
    with pytest.raises(ValueError):
        ch.composite([(a, blue), (np.zeros((3, 3)), green)])
    big = np.arange(10000, dtype=np.uint16).reshape(100, 100)
    lo, hi = ch.auto_levels(big, 1.0)
    assert 90 < lo < 110 and 9890 < hi < 9910
    assert ch.composite([(big, md.ChannelDisplay())], max_side=25).shape == (25, 25, 3)
    log = ch.adjustment_log([("DAPI", 16, blue), ("FITC", 16, md.ChannelDisplay(low=0, high=9, gamma=1.5))])
    assert log[0] == "Raw TIFFs unchanged (16 bit)" and log[1] == "DAPI: levels 0 to 1000, linear"
    assert "non-linear" in log[2]


def test_write_png_is_valid(tmp_path):
    rgb = np.zeros((3, 5, 3), dtype=np.uint8)
    rgb[1, 2] = (10, 20, 30)
    path = tmp_path / "p.png"
    ch.write_png(path, rgb)
    data = path.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n"
    w, h = struct.unpack(">II", data[16:24])
    assert (w, h) == (5, 3)
    idat = data[data.index(b"IDAT") + 4:data.index(b"IEND") - 8]
    raw = zlib.decompress(idat)
    assert raw[1 * 16 + 1 + 2 * 3: 1 * 16 + 1 + 2 * 3 + 3] == bytes((10, 20, 30))


# --- import ----------------------------------------------------------------

def test_group_and_import_session(tmp_path):
    src = tmp_path / "scope"
    src.mkdir()
    files = [_tif(src / f"s12_c{c:02d}_{dye}.tif", v) for c in (1, 2, 10) for dye, v in (("DAPI", 3000), ("FITC", 800), ("Cy3", 400))]
    (src / "notes.txt").write_text("x")
    _tif(src / "s12_c01_DAPI_2.tif", 5)
    lists = ls.load_lists(tmp_path)
    groups, unmatched = group_channel_files(files + [src / "notes.txt"], lists)
    assert [g.key for g in groups] == ["s12_c01", "s12_c02", "s12_c10"]
    assert groups[0].dyes == ["Cy3", "DAPI", "FITC"]
    assert [p.name for p in unmatched] == ["notes.txt"]
    data = tmp_path / "data"
    data.mkdir()
    slide = md.add_slide(data, "Slide 12", species="Crocus sativus")
    images, errors = import_cells(data, groups, lists=lists, slide_id=slide.id,
                                  metadata={"stage": "Mitotic metaphase"}, probes_by_dye={"FITC": "CroSat2"})
    assert errors == [] and len(images) == 3
    loaded = ci.load_images(data)
    d = md.image_details(loaded[0])
    assert d.slide_id == slide.id and d.pixel_size_um == pytest.approx(0.0645, rel=1e-3)
    assert {c.dye: c.probe for c in d.channels}["FITC"] == "CroSat2"
    assert all(md.channel_path(data, c).is_file() for c in d.channels)
    assert d.adjustment_log[0] == "Raw TIFFs unchanged (16 bit)"
    preview = ci.image_file_path(data, loaded[0])
    assert preview.read_bytes()[:4] == b"\x89PNG"
    d.channels[0].display.visible = False
    md.set_image_details(loaded[0], d)
    assert render_preview(data, loaded[0]) == preview
    ci.remove_image(data, loaded[0].id)
    assert not any(md.channel_path(data, c).exists() for c in d.channels)


def test_import_without_copy_and_bad_cell(tmp_path):
    src = tmp_path / "scope"
    src.mkdir()
    _tif(src / "a_DAPI.tif", 10)
    _tif(src / "a_FITC.tif", 10, shape=(10, 10))
    _tif(src / "b_DAPI.tif", 10)
    lists = ls.load_lists(tmp_path)
    groups, _ = group_channel_files(sorted(src.iterdir()), lists)
    images, errors = import_cells(tmp_path, groups, lists=lists, copy=False)
    assert len(images) == 1 and errors and errors[0].startswith("a:")
    channel = md.image_details(images[0]).channels[0]
    assert Path(channel.file).is_absolute() and md.channel_path(tmp_path, channel) == (src / "b_DAPI.tif").resolve()
    assert len(list(ci.store_dir(tmp_path).glob("*.tif"))) == 0


# --- signals ---------------------------------------------------------------

def _mark(mid="m1", chrom="2", hom="a"):
    return md.ChromosomeMark(id=mid, chromosome=chrom, homologue=hom, centromere=(0, 40), telomere_p=(0, 0), telomere_q=(0, 100))


@pytest.mark.parametrize("y, arm, position", [(38, "cen", "pericentromeric"), (20, "p", "interstitial"),
                                              (1, "p", "terminal"), (95, "q", "subterminal")])
def test_classify(y, arm, position):
    s = md.Signal(id="s", probe="P", x=1, y=y)
    sg.classify(s, _mark())
    assert (s.arm, s.position) == (arm, position)


def test_nearest_summary_and_manual():
    d1 = md.ImageDetails(chromosomes=[_mark("a", "1", "a"), md.ChromosomeMark(id="b", chromosome="2", homologue="a", x=200, y=0)])
    s = md.Signal(id="s", probe="CroSat2", x=198, y=3)
    sg.assign_and_classify(d1, s)
    assert s.chromosome_id == "b" and s.arm == ""
    far = md.Signal(id="f", probe="CroSat2", x=500, y=500)
    sg.assign_and_classify(d1, far, max_distance=20)
    assert far.chromosome_id == ""
    manual = md.Signal(id="m", probe="45S", x=0, y=99, arm="p", position="terminal", manual=True)
    sg.classify(manual, d1.chromosomes[0])
    assert manual.arm == "p"
    d1.signals = [s, md.Signal(id="t", probe="45S", x=0, y=99, chromosome_id="a")]
    sg.classify(d1.signals[1], d1.chromosomes[0])
    d2 = md.ImageDetails(chromosomes=[md.ChromosomeMark(id="c", chromosome="2", homologue="b"),
                                      md.ChromosomeMark(id="d", chromosome="3", homologue="a")],
                         signals=[md.Signal(id="u", probe="CroSat2", x=0, y=0, chromosome_id="c")])
    rows = sg.summarize([d1, d2], order=["2", "1", "3"])
    assert [(r.chromosome, r.probe, r.seen, r.total, r.cells) for r in rows] == [("2", "CroSat2", 2, 2, 2), ("1", "45S", 1, 1, 1)]
    assert rows[1].position == "terminal"
    assert sg.unsignalled_chromosomes([d1, d2], rows) == ["3"]
    assert sg.summary_tsv(rows).splitlines()[1] == "2\tCroSat2\t\t\t2\t2\t2"
