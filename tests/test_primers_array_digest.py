"""Tandem-array digest: ladders from site loss, methylation and partial cutting."""

import random

from repeatviewer.backends.primers import array_digest as ad
from repeatviewer.backends.primers import enzymes as ez


def _array(copies=30, lose=(3, 10, 17, 24)):
    r = random.Random(4)
    mono = "".join(r.choice("AT") for _ in range(170))
    mono = mono[:40] + "CCGG" + mono[44:]
    out = []
    for k in range(copies):
        m = list(mono)
        if k in lose:
            m[41] = "T"
        out.append("".join(m))
    return "".join(out)


def test_period_estimate_and_annotation():
    seq = _array()
    info = ad.array_info(seq)
    assert (info.period, info.copies, info.source) == (170, 30, "estimated")

    class F:
        def __init__(self, s, e):
            self.start, self.end, self.name = s, e, "CroSat3 monomer"

    feats = [F(1 + 170 * k, 170 * (k + 1)) for k in range(30)]
    info = ad.array_info(seq, feats)
    assert (info.period, info.copies, info.source) == (170, 30, "annotation")
    assert ad.array_info("ACGT" * 3).period == 0 or ad.array_info("ACGT").period == 0


def test_partial_fragments_math():
    frags = ad.partial_fragments(300, [101, 201], [1.0, 1.0])
    assert frags == {100: 3.0}           # three 100 bp pieces
    frags = ad.partial_fragments(300, [101, 201], [0.5, 0.5])
    assert abs(frags[100] - (0.5 + 0.25 + 0.5)) < 1e-9    # 1-100, 101-200, 201-300
    assert abs(frags[200] - 0.5) < 1e-9 and abs(frags[300] - 0.25) < 1e-9


def test_ladder_from_lost_sites_and_methylation():
    seq = _array()
    info = ad.array_info(seq)
    rows = ez.restriction_map(seq, ez.MapSettings(), None, names=["MspI"])
    s = ad.DigestSettings(completeness=1.0)
    lane = ad.digest_lane(seq, rows[0], info, s, use_data=False)
    by = {b.label: b for b in lane.bands}
    assert "1×" in by and "2×" in by and lane.sites == 26
    assert "Only 26 of about 30 copies" in lane.why
    every = ad.digest_lane(seq, rows[0], info, ad.DigestSettings(completeness=1.0,
                                                                 site_loss=ad.LOSS_EVERY),
                           use_data=False)
    assert every.added == 4 and "2×" not in {b.label for b in every.bands}

    calls = {}
    for i in range(len(seq) - 1):
        if seq[i:i + 2] == "CG":
            calls[i + 1] = (100.0, 10)
    meth = ez.MethData(calls=calls)
    rows = ez.restriction_map(seq, ez.MapSettings(), meth, names=["HpaII", "MspI"])
    lanes = ad.digest_lanes(seq, rows, info, s, has_data=True)
    assert [ln.label for ln in lanes] == ["HpaII (your data)", "MspI", "HpaII (no data)"]
    hpa = lanes[0]
    assert [b.label for b in hpa.bands] == ["uncut"] and hpa.blocked == 26
    doc = ad.gel_doc(lanes, len(seq), s)
    assert len(doc.lanes) == 4 and doc.lanes[0].kind == "ladder"
