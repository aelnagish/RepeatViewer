"""Primer workbench backend (region, advice, candidates)."""

import random

from repeatviewer.backends.primers import design as dz
from repeatviewer.backends.primers import workbench as wb


def _seq(n=3000, seed=5):
    r = random.Random(seed)
    return "".join(r.choice("ACGT") for _ in range(n))


def test_parse_region_and_primer3_args():
    r = wb.parse_region("CroSat3:1,850..2,250", 3000)
    assert (r.start, r.end, r.length) == (1850, 2250, 401)
    assert wb.parse_region("5000-6000", 3000) is None
    assert wb.primer3_region_args(r) == {"target": (1849, 2250)}
    assert wb.primer3_region_args(wb.Region(10, 20, "inside")) == {"included": (9, 20)}
    assert wb.primer3_region_args(wb.Region(10, 20, "avoid")) == {"excluded": [(9, 20)]}


def test_region_advice_and_fit():
    s = dz.DesignSettings()                     # product 100 to 1,000
    long = wb.Region(500, 1999)                 # 1,500 bp
    adv = wb.region_advice(long, s, 3000)
    assert adv and adv.can_fix and adv.product_min >= 1500
    fixed = wb.fit_settings(s, adv)
    assert wb.region_advice(long, fixed, 3000) is None
    assert wb.region_advice(wb.Region(1, 100), s, 3000).can_fix is False   # no room outside
    assert wb.region_advice(wb.Region(1000, 1300), s, 3000) is None


def test_candidates_cover_the_region(tmp_path):
    seq = _seq()
    region = wb.Region(1400, 1600)
    s = dz.DesignSettings(product_min=250, product_max=600)
    run = wb.design_candidates("t", seq, s, region)
    assert run.pairs
    assert all(wb.covers(p, region) for p in run.pairs)
    prims = wb.candidate_primers(run.pairs[0])
    assert len(prims) == 2
    fa = tmp_path / "t.fa"
    fa.write_text(">t\n" + seq + "\n")
    report = wb.write_candidates_report(fa, "t", seq, s, region, run)
    assert report.is_file() and "Region 1,400" in report.read_text(encoding="utf-8")


def test_primer_for_span_move_and_resize():
    from repeatviewer.backends.primers import on_sequence as osq
    from repeatviewer.backends.primers import thermo as th

    seq = _seq(600, 9)
    fwd = "gaattc" + seq[99:121]                       # tail + 22 nt at 100..121
    place = osq.place_primer(fwd, seq, name="F")
    assert (place.start, place.end, place.strand) == (100, 121, "+")
    moved = wb.primer_for_span(fwd, place, seq, 110, 131)
    assert moved == "gaattc" + seq[109:131]
    longer = wb.primer_for_span(fwd, place, seq, 100, 125)   # 3' end + 4
    assert longer == "gaattc" + seq[99:125]
    assert wb.primer_for_span(fwd, place, seq, 100, 105) == fwd   # below 10 nt: unchanged
    rev = th.reverse_complement(seq[399:420])
    rplace = osq.place_primer(rev, seq, name="R")
    assert rplace.strand == "-"
    assert wb.primer_for_span(rev, rplace, seq, 405, 425) == th.reverse_complement(seq[404:425])
