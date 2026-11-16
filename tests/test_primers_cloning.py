"""Cloning a PCR product into a vector: ends, orientations, colony PCR, warnings."""

import random

from repeatviewer.backends.primers import cloning as cl
from repeatviewer.backends.primers import tails as tl
from repeatviewer.backends.primers import thermo as th

_R = random.Random(8)
INSERT = "".join(_R.choice("ACGT") for _ in range(900))
# the same, with no site of any tail enzyme in it (checked below)
CLEAN = "".join(_R.choice("AT") for _ in range(300)) + "GGCCAATTGGCCAATT" * 20


def test_presets_and_methods():
    assert [v.name for v in cl.PRESETS][:2] == ["pGEM-T Easy", "pJET1.2/blunt"]
    assert cl.preset("pGEM-T Easy").takes(cl.METHOD_TA)
    assert not cl.preset("pGEM-T Easy").takes(cl.METHOD_BLUNT)
    assert cl.preset("pJET1.2/blunt").takes(cl.METHOD_BLUNT)
    assert cl.preset("not a vector").name == "Your own vector"
    assert cl.preset("Your own vector").takes(cl.METHOD_RE)


def test_ta_with_a_proofreading_polymerase_is_flagged():
    v = cl.preset("pGEM-T Easy")
    ok = cl.plan(INSERT, v, cl.METHOD_TA, polymerase="Taq", empty_product=250)
    assert ok.warnings == [] and "3′ A" in ok.ends
    assert any("A-tail" in n or "Taq" in n for n in ok.notes)
    bad = cl.plan(INSERT, v, cl.METHOD_TA, polymerase="Phusion / Q5")
    assert any("blunt ends" in w for w in bad.warnings)
    blunt = cl.plan(INSERT, cl.preset("pJET1.2/blunt"), cl.METHOD_BLUNT, polymerase="Taq")
    assert any("3′ A" in w for w in blunt.warnings)
    wrong = cl.plan(INSERT, cl.preset("pGEM-T Easy"), cl.METHOD_BLUNT)
    assert any("does not take" in w for w in wrong.warnings)


def test_both_orientations_plasmid_and_colony_sizes():
    v = cl.preset("pGEM-T Easy")
    res = cl.plan(INSERT, v, cl.METHOD_TA, polymerase="Taq", empty_product=250)
    assert [c.orientation for c in res.clones] == ["forward", "reverse"]
    assert all(c.total_length == 3015 + len(INSERT) for c in res.clones)
    assert all(c.colony_product == 250 + len(INSERT) for c in res.clones)
    tone, line = cl.summary(res, v, cl.METHOD_TA)
    assert tone == "green" and "colony PCR 1,150 bp" in line

    unknown = cl.plan(INSERT, cl.preset("Your own vector"), cl.METHOD_TA, polymerase="Taq")
    assert unknown.clones[0].total_length == 0 and unknown.clones[0].colony_product == 0
    assert any("empty vector" in n for n in unknown.notes)
    assert any("vector's length" in n for n in unknown.notes)


def test_restriction_ends_and_an_enzyme_cutting_the_insert():
    v = cl.preset("Your own vector")
    assert not tl.cuts_in("EcoRI", CLEAN) and not tl.cuts_in("BamHI", CLEAN)
    fwd = tl.tail_for("EcoRI", binding=CLEAN[:20]) + CLEAN[:20]
    rev = tl.tail_for("BamHI", binding="x") + th.reverse_complement(CLEAN[-20:])
    clean = cl.plan(CLEAN, v, cl.METHOD_RE, enzymes=("EcoRI", "BamHI"),
                    forward_primer=fwd, reverse_primer=rev)
    assert clean.warnings == [] and "EcoRI on the left" in clean.ends

    missing = cl.plan(CLEAN, v, cl.METHOD_RE, enzymes=("EcoRI", ""))
    assert any("Choose the enzyme" in w for w in missing.warnings)
    no_site = cl.plan(CLEAN, v, cl.METHOD_RE, enzymes=("EcoRI", "BamHI"),
                      forward_primer="ATGC" + CLEAN[:20])
    assert any("no restriction site" in w for w in no_site.warnings)

    cut = CLEAN[:400] + "GAATTC" + CLEAN[400:]
    inside = cl.plan(cut, v, cl.METHOD_RE, enzymes=("EcoRI", "BamHI"))
    assert any("also cuts inside the insert" in w for w in inside.warnings)
    same = cl.plan(CLEAN, v, cl.METHOD_RE, enzymes=("EcoRI", "EcoRI"))
    assert any("either way" in n for n in same.notes)


def test_tandem_array_warning():
    r = random.Random(3)
    mono = "".join(r.choice("ACGT") for _ in range(170))
    res = cl.plan(mono * 8, cl.preset("pGEM-T Easy"), cl.METHOD_TA, polymerase="Taq")
    assert any("tandem array" in w and "recA" in w for w in res.warnings)
    assert cl.plan(mono, cl.preset("pGEM-T Easy"), cl.METHOD_TA,
                   polymerase="Taq").warnings == []


def test_reads_gel_and_report(tmp_path):
    v = cl.preset("pGEM-T Easy")
    res = cl.plan(INSERT, v, cl.METHOD_TA, polymerase="Taq", empty_product=250)
    assert cl.reads_into_insert(250, 700) == 575
    assert cl.reads_into_insert(0, 700) == 700

    doc = cl.colony_gel(res, title="x")
    assert [ln.label for ln in doc.lanes][1:] == ["No insert", "With insert"]
    assert [b.size for b in doc.lanes[2].bands] == [1150]
    assert cl.colony_gel(cl.plan(INSERT, v, cl.METHOD_TA)) is None

    out = cl.write_report(tmp_path / "c.html", "CroSat3 product", res, v, cl.METHOD_TA)
    text = out.read_text()
    assert "CroSat3 product" in text and "pGEM-T Easy" in text and "M13 forward" in text
    assert "1,150 bp" in text and "\u2014" not in text
