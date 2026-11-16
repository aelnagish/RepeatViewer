"""In-silico PCR, gel, photo comparison, genome view and bench card."""

from __future__ import annotations

import os
import random

import numpy as np
import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from repeatviewer.backends.primers import bench, gel_doc as gd, gel_photo as gp, pcr  # noqa: E402
from repeatviewer.backends.primers import thermo as th  # noqa: E402
from repeatviewer.backends.primers.oligos import Primer  # noqa: E402

pytest.importorskip("primer3")


def _rnd(n: int, seed: int) -> str:
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


MONO = _rnd(306, 1)
ARRAY = _rnd(2000, 2) + MONO * 6 + _rnd(2000, 3)
F = Primer("CroSat2_F", "atggtacc" + MONO[20:42])
R = Primer("CroSat2_R", th.reverse_complement(MONO[260:281]))


@pytest.fixture()
def template(tmp_path):
    path = tmp_path / "chr3_region.fa"
    path.write_text(">chr3\n" + ARRAY + "\n")
    return path


def test_ladder_pcr_and_outputs(template):
    res = pcr.run_pcr([template], [F, R], pcr.PcrSettings(mode=pcr.MODE_LADDER))
    assert res.monomer == 306
    labels = [p.label for p in res.products if p.kind in ("expected", "ladder")]
    assert labels[:3] == ["1×", "2×", "3×"]
    first = res.products[0]
    assert first.sequence.startswith("ATGGTACC") and first.size == len(first.sequence)
    assert first.sites == 6 and first.yield_ == 1.0
    out = pcr.write_outputs(res)
    assert out["fasta"].name == "chr3_region.fa_PCR_CroSat2_F+R.fa"
    assert (out["details"] / f"{out['base']}_history.json").is_file()
    back = pcr.load_result(out["result_json"])
    assert len(back.products) == len(res.products)
    from repeatviewer.backends.primers.on_sequence import load_sequence_primers
    assert {p.name for p in load_sequence_primers(out["fasta"])} == {"CroSat2_F", "CroSat2_R"}
    again = pcr.write_outputs(res)
    assert again["base"].endswith("_2")


def test_removed_modes_open_as_standard():
    for old in ("bisulfite", "presence"):
        assert pcr.PcrSettings(mode=old).mode == pcr.MODE_STANDARD
        assert old not in pcr.MODES


def test_gel_document(tmp_path, template):
    res = pcr.run_pcr([template], [F, R], pcr.PcrSettings(mode=pcr.MODE_LADDER))
    doc = gd.gel_from_pcr(res)
    assert doc.lanes[0].kind == gd.LANE_LADDER and len(doc.lanes) == 2
    drawn = gd.layout(doc)
    assert drawn and all(0 <= d.y <= 1 for d in drawn)
    # bigger bands run less far
    a = gd.migration(300, 1.5)
    b = gd.migration(3000, 1.5)
    assert a > b
    assert abs(gd.size_at(doc, gd.migration(1000, doc.agarose, doc.run_factor)) - 1000) < 1
    doc.lanes[1].cut_with = "KpnI"
    cut = gd.lane_bands(doc.lanes[1])
    assert any(b.size == 1 or b.size < 20 for b in cut) or len(cut) > len(doc.lanes[1].bands)
    path = doc.save(tmp_path / "x.rvgel")
    assert gd.GelDoc.load(path).lanes[1].cut_with == "KpnI"
    assert "<svg" in gd.to_svg(doc)


def test_photo_reading():
    H, W = 600, 300
    img = np.full((H, W), 30.0)

    def put(x, size, amp):
        y = int(40 + gd.migration(size, 1.5, 0.88) * 520)
        for dy in range(-3, 4):
            img[y + dy, x - 25:x + 25] += amp * np.exp(-dy * dy / 4)

    for s, m in gd.LADDERS["GeneRuler 1 kb Plus"]:
        put(60, s, 80 * m)
    for s in (284, 590, 896, 450, 765):
        put(150, s, 120)
    img += np.random.default_rng(0).normal(0, 3, img.shape)
    lanes = gp.find_lanes(img)
    assert len(lanes) == 2
    fit = gp.fit_ladder([b.y for b in gp.find_bands(img, lanes[0])],
                        [s for s, _ in gd.LADDERS["GeneRuler 1 kb Plus"]])
    sizes = sorted(round(fit.size_at(b.y)) for b in gp.find_bands(img, lanes[1]))
    assert all(abs(a / b - 1) < 0.07 for a, b in zip(sizes, [284, 450, 590, 765, 896]))
    rows = gp.compare(sizes, [(284, "1×", 1), (590, "2×", 1), (896, "3×", 1), (447, "off-target", .1)])
    status = {s for _o, _e, s, _l in rows}
    assert {"matched", "off-target", "unexplained"} <= status


def test_bench_numbers():
    mix = bench.mastermix([F, R], reactions=10, volume_ul=25)
    assert abs(sum(r.per_reaction for r in mix) - 24) < 0.05
    steps, cycles, why = bench.program([F, R], longest_bp=1500)
    assert cycles == 30 and steps[3].seconds == 90 and "Annealing" in why


def test_genome_sites(template):
    recs = pcr.genome_sites(template, [F, R])
    assert recs[0].counts().get(0) == 12 and recs[0].products
    assert "primer_bind" in pcr.sites_gff3(recs)


def test_windows_build(tmp_path, template):
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    from repeatviewer.views.bench_card_dialog import BenchCardDialog
    from repeatviewer.views.gel_photo_dialog import GelPhotoDialog
    from repeatviewer.views.pcr_dialog import PcrWindow
    from repeatviewer.views.primer_genome_dialog import PrimerGenomeDialog
    from repeatviewer.widgets.gel_view import GelDocumentPanel

    res = pcr.run_pcr([template], [F, R], pcr.PcrSettings(mode=pcr.MODE_LADDER))
    out = pcr.write_outputs(res)
    win = PcrWindow(pcr.PcrSettings(), templates=[template], primers=[F, R])
    win.show_result(res, out)
    assert win.table.rowCount() >= 6
    doc = gd.gel_from_pcr(res)
    path = doc.save(tmp_path / "g.rvgel")
    panel = GelDocumentPanel()
    panel.load(path)
    panel.lanes_table.selectRow(1)
    panel._duplicate_lane()
    assert len(gd.GelDoc.load(path).lanes) in (2, 3)
    BenchCardDialog([F, R], result=res)
    dlg = GelPhotoDialog(doc)
    from PySide6.QtGui import QImage
    img = QImage(200, 300, QImage.Format.Format_Grayscale8)
    img.fill(20)
    dlg.load_image(img)
    gdlg = PrimerGenomeDialog([F, R], genome=template)
    gdlg._on_done(pcr.genome_sites(template, [F, R]))
    app.processEvents()
