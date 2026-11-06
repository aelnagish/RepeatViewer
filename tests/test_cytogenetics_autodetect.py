"""Signals auto-detection: chromosome segmentation, spot detection, review edits, pairing, karyotype from cells."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("scipy")

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import karyotype as kt
from repeatviewer.backends.cytogenetics import pairing as pr
from repeatviewer.backends.cytogenetics import segment as sm
from repeatviewer.backends.cytogenetics import signals_session as ss
from repeatviewer.backends.cytogenetics import spots as sp
from repeatviewer.backends.cytogenetics.importer import ImportRequest, run_request
from repeatviewer.backends.cytogenetics.model import UNNAMED

tifffile = pytest.importorskip("tifffile")

# three rod chromosomes (x0, y0, x1, y1, width) and a signal on the first and third
RODS = [(20, 30, 20, 150, 10), (70, 30, 70, 120, 10), (120, 40, 120, 100, 10)]
SPOTS = [(20, 50), (120, 90)]


def _rods(shape=(180, 160), touching=False):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    img = np.full(shape, 10.0)
    rods = list(RODS)
    if touching:
        rods.append((26, 30, 26, 150, 10))  # a second rod right next to the first
    for x0, y0, x1, y1, w in rods:
        d = np.where((yy >= y0) & (yy <= y1), np.abs(xx - x0), 1e9)
        img += 200.0 * np.exp(-(d / (w / 2.2)) ** 4)
    # a narrowing (centromere) on rod 1 at y = 70
    img[66:74, 14:27] *= 0.55
    return img


def _spots(shape=(180, 160)):
    yy, xx = np.mgrid[0:shape[0], 0:shape[1]]
    img = np.full(shape, 5.0)
    for x, y in SPOTS:
        img += 400.0 * np.exp(-((xx - x) ** 2 + (yy - y) ** 2) / 4.0)
    return img


def test_mask_round_trip():
    m = np.zeros((7, 9), dtype=bool)
    m[2:5, 3:8] = True
    m[0, 0] = True
    enc = sm.encode_mask(m, 4, 5)
    back, x0, y0 = sm.decode_mask(enc)
    assert (x0, y0) == (4, 5) and np.array_equal(back, m)
    assert sm.decode_mask({"x": 0}) is None


def test_detects_rods_with_midline_and_centromere():
    found = sm.detect_chromosomes(_rods())
    assert len(found) == 3
    longest = found[0]
    assert 105 <= longest.length <= 130  # rod 1 is 120 px
    xs = [p[0] for p in longest.midline]
    assert max(xs) - min(xs) < 6  # straight, vertical
    cy = longest.midline[longest.centromere][1]
    assert 60 <= cy <= 80  # the narrowing


def test_split_and_merge():
    found = sm.detect_chromosomes(_rods())
    rod = found[0]
    pieces = sm.split_along(rod.mask, rod.x0, rod.y0, (5, 100), (35, 100))
    assert len(pieces) == 2
    joined, _x, _y = sm.union(pieces)
    assert abs(int(joined.sum()) - rod.area) < rod.area * 0.1


def test_spots_on_chromosomes():
    labels = sm.label_image([(o.mask, o.x0, o.y0) for o in sm.detect_chromosomes(_rods())], (180, 160))
    found = sp.find_spots(_spots(), sp.SpotSettings(), labels)
    assert len(found) == 2
    assert all(sp.assign((s.x, s.y), labels) >= 0 for s in found)


def _folder(tmp_path):
    src = tmp_path / "scope"
    src.mkdir()
    tifffile.imwrite(str(src / "c01_DAPI.tif"), _rods().astype(np.uint16))
    tifffile.imwrite(str(src / "c01_FITC.tif"), _spots().astype(np.uint16))
    work = tmp_path / "data"
    work.mkdir()
    run_request(ImportRequest(folder=work, source=src, slide_name="S1", probes_by_dye={"FITC": "CroSat2"}))
    return work


def test_session_detect_review_and_save(tmp_path):
    work = _folder(tmp_path)
    s = ss.SignalsSession(work, ci.load_images(work)[0])
    rows = s.probe_rows()
    assert [r.probe for r in rows] == ["CroSat2"]
    assert s.detect_chromosomes() == 3
    assert all(m.chromosome == UNNAMED and m.midline for m in s.details.chromosomes)
    assert s.detect_signals(rows[0]) == 2
    assert all(sig.chromosome_id and sig.arm for sig in s.details.signals)
    first = s.details.chromosomes[0]
    ids = s.split_mark(first.id, (5, 100), (35, 100))
    assert len(ids) == 2 and len(s.details.chromosomes) == 4
    keep = s.merge_marks(ids)
    assert len(s.details.chromosomes) == 3 and s.details.chromosome(keep) is not None
    s.move_centromere(keep, 20, 110)
    assert s.details.chromosome(keep).centromere is not None
    s.accept_all()
    assert s.to_check() == []
    s.save()
    again = ss.SignalsSession(work, ci.load_images(work)[0])
    assert len(again.details.chromosomes) == 3 and len(again.details.signals) == 2
    assert again.details.chromosomes[0].mask is not None


def test_pairing_and_karyotype_from_cells(tmp_path):
    work = _folder(tmp_path)
    s = ss.SignalsSession(work, ci.load_images(work)[0])
    s.detect_chromosomes()
    s.detect_all_signals()
    p = pr.pair_without_karyotype(s.details, 1)
    pr.apply_pairing(s.details, p)
    assert sorted(m.chromosome for m in s.details.chromosomes) == ["1", "2", "3"]
    k = kt.new_karyotype("Test", 3, 1)
    for c, length in zip(k.chromosomes, (12, 9, 6)):
        c.length = length
    p2 = pr.pair_with_karyotype(s.details, k)
    pr.apply_pairing(s.details, p2)
    longest = max(s.details.chromosomes, key=lambda m: sm.polyline_length(m.midline))
    assert longest.chromosome == "1"
    stats = pr.stats_over_cells([s.details], ["1", "2", "3"])
    assert [x.name for x in stats] == ["1", "2", "3"]
    made = pr.karyotype_from_cells(stats, "From cells", ploidy=1)
    assert [c.name for c in made.chromosomes] == ["1", "2", "3"]
    assert any(b.kind == "FISH CroSat2" for c in made.chromosomes for b in c.bands)
    strip = pr.straighten(s.counterstain(), longest.midline, pr.half_width(longest))
    assert strip.shape[0] >= 100
