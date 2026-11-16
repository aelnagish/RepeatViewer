"""Restriction map with methylation from data (primer workbench Enzymes tab)."""

from repeatviewer.backends.primers import enzymes as ez

# HpaII/MspI site CCGG at 7..10: outer C 7 (CHG), inner C 8 (CG), G 9 / G 10 bottom Cs.
SEQ = "AAAATTCCGGTTAAAGAATTCAAACTGCAGAAA"


def _site(name, start):
    return next(s for s in ez.find_sites(name, SEQ) if s.start == start)


def test_sites_and_cuts():
    assert [(s.start, s.end, s.cut) for s in ez.find_sites("EcoRI", SEQ)] == [(16, 21, 17)]
    assert [(s.start, s.cut) for s in ez.find_sites("PstI", SEQ)] == [(25, 30)]
    assert ez.fragments(100, [11, 51]) == [10, 40, 50]
    assert "EcoRI" in ez.set_names(ez.SET_COMMON6)
    assert {"HpaII", "MspI"} <= set(ez.set_names(ez.SET_METHYL))


def test_contexts_both_strands():
    assert ez.c_context(SEQ, 6) == ("+", "CHG")    # C7 of CCGG
    assert ez.c_context(SEQ, 7) == ("+", "CG")     # C8
    assert ez.c_context(SEQ, 8) == ("-", "CG")     # G9 = bottom C of the CG
    assert ez.c_context(SEQ, 9) == ("-", "CHG")    # G10 = bottom outer C


def test_hpaii_mspi_pair():
    cg = ez.MethData(calls={8: (90.0, 12)})
    chg = ez.MethData(calls={7: (90.0, 12)})
    for meth, hpa, msp in ((cg, True, False), (chg, True, True)):
        a, b = _site("HpaII", 7), _site("MspI", 7)
        ez.judge_site(a, SEQ, meth, min_percent=50, min_reads=5)
        ez.judge_site(b, SEQ, meth, min_percent=50, min_reads=5)
        assert (a.blocked, b.blocked) == (hpa, msp)
    low = _site("HpaII", 7)
    ez.judge_site(low, SEQ, ez.MethData(calls={8: (90.0, 3)}), min_percent=50, min_reads=5)
    assert not low.blocked and low.calls


def test_merged_cpg_row_counts_for_both_strands():
    s = _site("HpaII", 7)
    ez.judge_site(s, SEQ, ez.MethData(calls={8: (80.0, 10)}), min_percent=50, min_reads=5)
    assert [c.pos for c in s.calls] == [8]


def test_pst_chg_and_no_data():
    s = _site("PstI", 25)
    ez.judge_site(s, SEQ, ez.MethData(calls={28: (70.0, 9)}), min_percent=50, min_reads=5)
    assert s.blocked and "mCHG" in s.reason
    e = _site("EcoRI", 16)
    ez.judge_site(e, SEQ, None, min_percent=50, min_reads=5)
    assert not e.blocked and "mCG" in e.reason


def test_load_bismark_and_bedmethyl(tmp_path):
    cov = tmp_path / "x.bismark.cov"
    cov.write_text("t\t8\t8\t90\t9\t1\nother\t8\t8\t0\t0\t5\n")
    m = ez.load_methylation(cov, "t", len(SEQ))
    assert m.calls == {8: (90.0, 10)}
    bed = tmp_path / "x.bedmethyl"
    cols = ["chr1", "1006", "1007", "m", "10", "+", "1006", "1007", "255,0,0", "10", "90.0",
            "9", "1", "0", "0", "0", "0", "0"]
    bed.write_text("\t".join(cols) + "\n")
    m = ez.load_methylation(bed, "chr1:1000-1032", len(SEQ))
    assert m.calls == {8: (90.0, 10)}
    assert "No methylation rows" in ez.load_methylation(cov, "zz", len(SEQ)).note


def test_map_and_listing():
    s = ez.MapSettings(enzyme_set=ez.SET_METHYL, show=ez.SHOW_MANY)
    rows = ez.restriction_map(SEQ, s, ez.MethData(calls={8: (90.0, 12)}))
    by = {r.name: r for r in rows}
    assert by["HpaII"].blocked == 1 and by["MspI"].blocked == 0
    listed = ez.listed_rows(rows, s.show)
    assert listed[0].cuts >= 1 and listed[-1].cuts == 0
