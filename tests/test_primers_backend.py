"""Tools > Primers back end: thermodynamics, oligos, binding search, design modes."""

from __future__ import annotations

import random

import pytest

primer3 = pytest.importorskip("primer3")

from primer_testdata import make_data, revcomp  # noqa: E402

from repeatviewer.backends.primers import binding as bd  # noqa: E402
from repeatviewer.backends.primers import design as dz  # noqa: E402
from repeatviewer.backends.primers import gel  # noqa: E402
from repeatviewer.backends.primers import thermo as th  # noqa: E402
from repeatviewer.backends.primers.oligos import (  # noqa: E402
    KIND_FORWARD,
    KIND_OLIGO,
    KIND_PROBE,
    KIND_REVERSE,
    Primer,
    PrimerDatabase,
    kind_from_name,
    load_primers,
    parse_primer_text,
    write_primers,
)


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    return make_data(tmp_path_factory.mktemp("primer_data"))


def _records(path):
    return list(bd.iter_fasta(path))


# --------------------------------------------------------------------------
# Thermodynamics
# --------------------------------------------------------------------------


def test_tm_is_primer3_on_the_binding_part():
    seq = "AGCTAGCTAGCATCGATCGAT"
    assert th.calc_tm(seq) == pytest.approx(primer3.calc_tm(seq), abs=0.01)
    # A lower-case 5' tail does not change the template Tm.
    assert th.calc_tm("gtaaaacgacggccagt" + seq) == pytest.approx(th.calc_tm(seq), abs=0.01)


@pytest.mark.parametrize("tm_method", th.TM_METHODS)
@pytest.mark.parametrize("salt", th.SALT_METHODS)
def test_tm_settings_reach_primer3(tm_method, salt):
    s = th.ThermoSettings(tm_method=tm_method, salt_method=salt, divalent_mM=3.0)
    seq = "ACGTTGCAAGTCCGATTAGC"
    expected = primer3.calc_tm(seq, mv_conc=50, dv_conc=3.0, dntp_conc=0.6, dna_conc=50,
                               tm_method=tm_method, salt_corrections_method=salt)
    assert th.calc_tm(seq, s) == pytest.approx(expected, abs=0.01)


def test_degenerate_tm_is_mean_of_variants():
    lo, hi = th.calc_tm_range("ACGTRGCAAGTCCGATTAGC")
    assert lo < hi
    assert lo <= th.calc_tm("ACGTRGCAAGTCCGATTAGC") <= hi


def test_clean_oligo_and_tail():
    assert th.clean_oligo("5'-acg TTU-3'") == "acgTTT"
    assert th.split_tail("gtaaACGT") == ("gtaa", "ACGT")
    assert th.split_tail("acgt") == ("", "ACGT")  # all lower case: all of it binds
    with pytest.raises(th.PrimerError):
        th.clean_oligo("ACGTX")


def test_hairpin_is_dot_bracket_and_dimer_has_pairs():
    hp = th.hairpin("GGGGCCCCAAAATTTTGGGGCCCC")
    assert hp.found and hp.dg < -5
    assert hp.lines[0].startswith("5′ GGGG")
    assert set(hp.lines[1].strip()) <= set("().")
    dimer = th.homodimer("TTGCAGCTAGCTAGCTGCAA")
    assert dimer.found and any("|" in ln for ln in dimer.lines)


def test_three_prime_dg_and_analysis_warnings():
    assert th.three_prime_dg("AAAAAGCGCG") < th.three_prime_dg("GCGCGAAAAA")
    st = th.analyze_oligo("AAAAAAAAAAAAAAAAAAAA")
    assert any("run of" in w for w in st.warnings)
    assert st.binding_length == 20


def test_no_em_dash_in_user_text():
    """Rule: no em dash in UI / help / report text (reads as AI-written)."""
    import inspect

    from repeatviewer import primer_actions
    from repeatviewer.backends.primers import jobs, report
    from repeatviewer.views import primer_dialogs

    for module in (th, bd, dz, gel, jobs, report, primer_dialogs, primer_actions):
        assert "—" not in inspect.getsource(module), module.__name__


# --------------------------------------------------------------------------
# Oligos and the Primer Database
# --------------------------------------------------------------------------


def test_kind_from_name():
    assert kind_from_name("CroSat2_F") == KIND_FORWARD
    assert kind_from_name("CroSat2-R2") == KIND_REVERSE
    assert kind_from_name("Sat probe") == KIND_PROBE
    assert kind_from_name("random") == KIND_OLIGO


def test_parse_primer_text_shapes():
    tsv = parse_primer_text("A_F\tACGTACGTACGT\nA_R  TTTTGGGGCCCC\nGGGGAAAACCCCTTTT")
    assert [p.name for p in tsv] == ["A_F", "A_R", "primer_3"]
    assert [p.kind for p in tsv][:2] == [KIND_FORWARD, KIND_REVERSE]
    fa = parse_primer_text('>X_F kind=probe note="ladder"\nacgtACGTACGTAC\n')
    assert fa[0].kind == KIND_PROBE and fa[0].note == "ladder"
    assert fa[0].tail == "acgt" and fa[0].binding == "ACGTACGTAC"


def test_database_roundtrip_dedupe_and_names(tmp_path):
    path = tmp_path / "db.fa"
    db = PrimerDatabase(path)
    added = db.add([Primer("P_F", "ACGTACGTACGTAAA"), Primer("P_F", "GGGGCCCCAAAATTT"),
                    Primer("Other", "ACGTACGTACGTAAA")])
    # Same sequence is not added twice; a clashing name gets _2.
    assert [p.name for p in added] == ["P_F", "P_F_2"]
    db.save()
    again = PrimerDatabase(path)
    assert again.names() == ["P_F", "P_F_2"]
    assert again.get("P_F").kind == KIND_FORWARD
    assert again.remove(["P_F_2"]) == 1
    write_primers(tmp_path / "x.fa", again.primers)
    assert [p.name for p in load_primers(tmp_path / "x.fa")] == ["P_F"]


# --------------------------------------------------------------------------
# Binding search and in-silico PCR
# --------------------------------------------------------------------------


def _plant(seq: str, at: int, site: str) -> str:
    return seq[:at] + site + seq[at + len(site):]


def test_sites_both_strands_mismatches_and_three_prime_rule():
    rng = random.Random(3)
    tpl = "".join(rng.choice("ACGT") for _ in range(3000))
    primer = "ACGTTGCAAGTCCGATTAGC"
    one_mm = primer[:5] + ("A" if primer[5] != "A" else "C") + primer[6:]
    end_mm = primer[:-1] + ("A" if primer[-1] != "A" else "C")
    tpl = _plant(tpl, 100, primer)
    tpl = _plant(tpl, 900, revcomp(one_mm))
    tpl = _plant(tpl, 2000, end_mm)
    st = bd.SearchSettings(max_mismatches=1, three_prime_exact=3)
    sites = bd.find_sites_in_sequence("p", primer, "t", tpl, st)
    got = {(s.strand, s.start, s.mismatches) for s in sites}
    assert ("+", 101, 0) in got
    assert ("-", 901, 1) in got
    # A mismatch in the last three 3' bases stops the polymerase: not a site.
    assert not any(s.start == 2001 for s in sites)
    minus = next(s for s in sites if s.strand == "-")
    assert minus.mismatch_positions == [6]
    assert minus.site == one_mm or revcomp(tpl[900:920]) == minus.site


def test_seeded_and_exhaustive_agree():
    rng = random.Random(5)
    tpl = "".join(rng.choice("ACGT") for _ in range(20_000))
    primer = tpl[5000:5020]
    tpl = _plant(tpl, 12_000, revcomp(primer[:3] + "T" + primer[4:]))
    st = bd.SearchSettings(max_mismatches=2, three_prime_exact=3)
    a = bd.find_sites_in_sequence("p", primer, "t", tpl, st, exhaustive=True)
    b = bd.find_sites_in_sequence("p", primer, "t", tpl, st, exhaustive=False)
    assert [(s.strand, s.start, s.mismatches) for s in a] == \
           [(s.strand, s.start, s.mismatches) for s in b]


@pytest.mark.parametrize("from_end", [4, 5, 6, 7, 8, 10, 15])
def test_genome_search_finds_mismatches_inside_the_seed(from_end):
    """The seeded genome search must find every site the exhaustive one finds,
    including sites whose mismatch sits in the 3' k-mer used as the seed."""
    rng = random.Random(from_end)
    tpl = "".join(rng.choice("ACGT") for _ in range(30_000))
    primer = "ACGTTGCAGGTCAATGCCTA"
    i = len(primer) - from_end
    site = primer[:i] + ("A" if primer[i] != "A" else "C") + primer[i + 1:]
    tpl = _plant(tpl, 10_000, site)
    tpl = _plant(tpl, 20_000, revcomp(site))
    st = bd.SearchSettings(max_mismatches=2, three_prime_exact=3)
    a = bd.find_sites_in_sequence("p", primer, "t", tpl, st, exhaustive=True)
    b = bd.GenomeScanner([("p", primer)], st).scan("t", tpl)
    c = bd.find_sites_in_sequence("p", primer, "t", tpl, st, exhaustive=False)
    key = [(s.strand, s.start, s.mismatches) for s in a]
    assert len(key) == 2
    assert key == [(s.strand, s.start, s.mismatches) for s in b] == \
           [(s.strand, s.start, s.mismatches) for s in c]


def test_seed_versions():
    assert bd.seed_variants("ACGT", 2, 0, at_end=True) == ["ACGT"]
    one = bd.seed_variants("ACGT", 2, 1, at_end=True)
    assert len(one) == 1 + 2 * 3 and all(v.endswith("GT") for v in one)
    rc = bd.seed_variants("ACGT", 2, 1, at_end=False)
    assert all(v.startswith("AC") for v in rc)
    assert bd.seed_variant_count("ACGTACGTACGT", 3, 2) == 1 + 9 * 3 + 36 * 9


def test_genome_scanner_finds_site_across_chunk_border(monkeypatch):
    monkeypatch.setattr(bd, "CHUNK_BP", 1000)
    rng = random.Random(9)
    tpl = "".join(rng.choice("ACGT") for _ in range(5000))
    primer = "ACGTTGCAAGTCCGATTAGC"
    tpl = _plant(tpl, 990, primer)  # crosses the 1000 bp chunk border
    sc = bd.GenomeScanner([("p", primer)], bd.SearchSettings(max_mismatches=0))
    sites = sc.scan("chr", tpl)
    assert [(s.start, s.strand) for s in sites if s.mismatches == 0] == [(991, "+")]


def test_template_n_never_matches_and_degenerate_primer_does():
    tpl = "A" * 50 + "ACGTTGCAAGTCCGATTAGC" + "N" * 30 + "A" * 50
    st = bd.SearchSettings(max_mismatches=0, three_prime_exact=0)
    assert bd.find_sites_in_sequence("deg", "ACGTTGCARGTCCGATTAGC", "t", tpl, st)
    assert not bd.find_sites_in_sequence("n", "NNNNNNNNNNNNNNNNNNNN", "t", "N" * 40, st)


def test_products_and_single_primer_products():
    rng = random.Random(11)
    tpl = "".join(rng.choice("ACGT") for _ in range(2000))
    f, r = tpl[100:120], revcomp(tpl[480:500])
    st = bd.SearchSettings(max_mismatches=0)
    sites = bd.find_sites_in_sequence("F", f, "t", tpl, st) + \
        bd.find_sites_in_sequence("R", r, "t", tpl, st)
    prods = bd.products_for_pair(sites, "F", "R", max_size=1000)
    assert [(p.start, p.end, p.size) for p in prods] == [(101, 500, 400)]
    # An inverted copy of F lets F prime on its own.
    tpl2 = tpl[:1500] + revcomp(f) + tpl[1520:]
    sites2 = bd.find_sites_in_sequence("F", f, "t", tpl2, st)
    own = bd.products_for_pair(sites2, "F", "R", max_size=2000)
    assert any(p.forward == "F" and p.reverse == "F" for p in own)


def test_consensus_threshold_and_gaps():
    rows = ["ACGT-A", "ACGTTA", "ATGTTA", "ACGTTA"]
    assert bd.consensus(rows) == "ACGTTA"
    # 100 % needs C and T together at column 2: Y.
    assert bd.consensus(rows, threshold=1.0)[1] == "Y"


# --------------------------------------------------------------------------
# Design modes (coordinates must point at the primer on the template)
# --------------------------------------------------------------------------


def _check_on_template(pair: dz.DesignedPair, template: str) -> None:
    t = template.upper()
    for o in (pair.forward, pair.reverse, pair.probe):
        if o is None:
            continue
        bind = th.binding_part(o.sequence)
        site = t[o.start - 1:o.end]
        if o.strand == "+":
            assert bind == site, (o.name, bind, site)
        else:
            assert bind == revcomp(site), (o.name, bind, revcomp(site))


def _settings(mode: str, **kw) -> dz.DesignSettings:
    return dz.DesignSettings(mode=mode, **{**dz.mode_preset(mode), **kw})


def test_pcr_pairs_sit_on_the_template(data):
    (rid, seq), = _records(data.gene)
    s = _settings(dz.MODE_PCR)
    run = dz.run_primer3(rid, seq, s)
    assert len(run.pairs) == s.num_return
    for p in run.pairs:
        _check_on_template(p, seq)
        assert p.product_size == p.reverse.end - p.forward.start + 1
        assert abs(p.forward.tm - p.reverse.tm) <= s.max_tm_diff + 0.01


def test_only_the_kept_modes_are_offered():
    assert dz.DESIGN_MODES == (dz.MODE_PCR, dz.MODE_QPCR, dz.MODE_LADDER, dz.MODE_FISH)
    assert dz.MODE_LABELS[dz.MODE_QPCR] == "qPCR"
    listed = [m for _f, groups in dz.MODE_FAMILIES for _g, modes in groups for m in modes]
    assert listed == [dz.MODE_PCR, dz.MODE_QPCR, dz.MODE_LADDER, dz.MODE_FISH, dz.MODE_PAINTING]
    # settings saved with a removed mode open as Standard PCR
    for old in ("bsp", "msp", "conserved", "specific", "insertion"):
        assert dz.DesignSettings.from_dict({"mode": old, "mask_tandem": True}).mode == dz.MODE_PCR


def test_target_and_tails(data):
    (rid, seq), = _records(data.gene)
    s = _settings(dz.MODE_PCR, forward_tail="gtaaaacgacggccagt", reverse_tail="caggaaacagctatgac")
    run = dz.run_primer3(rid, seq, s, target=(2000, 2050))
    assert run.pairs
    for p in run.pairs:
        assert p.forward.start <= 2001 and p.reverse.end >= 2050
        assert p.forward.sequence.startswith("gtaaaacgacggccagt")
        assert p.reverse.sequence.startswith("caggaaacagctatgac")
        _check_on_template(p, seq)


def test_qpcr_probe_inside_product(data):
    (rid, seq), = _records(data.gene)
    s = _settings(dz.MODE_QPCR, pick_probe=True)
    run = dz.run_primer3(rid, seq, s, pick_probe=True)
    assert run.pairs
    for p in run.pairs:
        assert 70 <= p.product_size <= 180
        assert p.probe is not None
        assert p.forward.end < p.probe.start and p.probe.end < p.reverse.start
        assert p.probe.tm > p.forward.tm + 4
        _check_on_template(p, seq)


def test_ladder_primers_face_outward(data):
    n = len(data.mono)
    s = _settings(dz.MODE_LADDER)
    run = dz.design_ladder("mono", data.mono, s)
    assert run.pairs
    dimer = data.mono + data.mono
    for p in run.pairs:
        f, r = p.forward, p.reverse
        assert 1 <= f.start <= n and 1 <= r.end <= n
        # Reverse sits upstream of forward in the monomer: they point apart.
        assert r.end < f.start
        assert p.product_size == (n - f.start + 1) + r.end
        assert th.binding_part(f.sequence) == dimer[f.start - 1:f.end]
        assert th.binding_part(r.sequence) == revcomp(dimer[n + r.start - 1:n + r.end])
    bands = [size for size, _w in dz.ladder_bands(run.pairs[0].product_size, n)]
    assert bands[1] - bands[0] == n


def test_fish_probes_are_in_every_monomer(data):
    members = _records(data.monomers)
    s = _settings(dz.MODE_FISH)
    rid, tseq = dz.representative_record(members)
    probes = dz.design_fish(rid, tseq, members, s)
    assert probes
    for p in probes:
        assert p.forward is None and p.reverse is None and p.probe is not None
        assert s.size_min <= p.probe.length <= s.size_max
        assert p.coverage["target_amplified"] / p.coverage["target_total"] >= s.min_coverage
        assert p.extra["order"].startswith("5'-Cy3-")
    spans = sorted((p.probe.start, p.probe.end) for p in probes)
    assert all(a[1] < b[0] for a, b in zip(spans, spans[1:])), "probes overlap"


# --------------------------------------------------------------------------
# Gel
# --------------------------------------------------------------------------


def test_gel_ladder_choice_and_band_order():
    assert gel.pick_ladder([120, 180]) == "25 bp (small)"
    assert gel.pick_ladder([300, 900]) == "100 bp"
    assert gel.pick_ladder([500, 8000]) == "1 kb Plus"
    lo, hi = 50.0, 5000.0
    assert gel.band_y(100, lo, hi, 0, 100) > gel.band_y(1000, lo, hi, 0, 100)
    bands = gel.merge_bands([200, 201, 202, 400])
    assert [b.size for b in bands] == [400, 201]
    svg = gel.gel_svg([gel.Lane("pair", bands)], ladder="100 bp")
    assert svg.startswith("<svg") and "pair" in svg


# ---------------------------------------------------------------- FISH hybridisation


def test_hyb_tm_formamide_and_salt():
    from repeatviewer.backends.primers import thermo as tt
    probe = "ATGCGTACGTTAGCATCGATCGTAGCTAGCTAGGCTAGCTAGCATCGATCGATCGTAGCTAG"[:50]
    c = tt.HybConditions()                      # 50 % formamide, 2x SSC
    assert c.sodium_mM() == 390.0 and abs(c.sodium_mM(wash=True) - 19.5) < 1e-9
    no_fa = tt.hyb_tm(probe, tt.HybConditions(formamide_pct=0))
    assert abs(no_fa - tt.hyb_tm(probe, c) - 30.0) < 0.01     # 0.60 °C per %
    assert tt.hyb_tm(probe, c, wash=True) < no_fa
    assert tt.HybConditions.from_dict(c.to_dict()) == c


def test_unsuitable_regions_found():
    import random
    from repeatviewer.backends.primers import low_complexity as lc
    rnd = random.Random(3)
    r = lambda n: "".join(rnd.choice("ACGT") for _ in range(n))  # noqa: E731
    stem = r(15)
    rc = stem[::-1].translate(str.maketrans("ACGT", "TGCA"))
    seq = r(100) + "AT" * 15 + r(100) + stem + r(30) + rc + r(100) + "A" * 25 + r(200)
    kinds = {(x.start, x.kind) for x in lc.unsuitable_regions(seq)}
    assert (99, "simple repeat") in kinds and (99, "low complexity") in kinds
    assert (230, "inverted repeat") in kinds and (275, "inverted repeat") in kinds
    assert (390, "simple repeat") in kinds
    assert not [x for x in lc.unsuitable_regions(r(3000)) if x.kind == "low complexity"]


def test_dimer_consensus_becomes_monomer():
    from repeatviewer.backends.primers import design as dz
    mono = "ACGTTGCAAGCTTAGCCTAGGATCCAAGTCGATCGGATCCTTAGCAAGTCC"
    assert dz.monomer_from_dimer(mono + mono) == (mono, True)
    assert dz.monomer_from_dimer(mono + mono[::-1])[1] is False


def test_ladder_turns_template_mispriming_off():
    from repeatviewer.backends.primers import design as dz
    assert all(v == -1 for v in dz.LADDER_NO_TEMPLATE_MISPRIMING.values())


def test_fish_probes_use_hybridisation_tm_and_skip_regions():
    import random
    from repeatviewer.backends.primers import design as dz
    from repeatviewer.backends.primers import low_complexity as lc
    from repeatviewer.backends.primers import thermo as tt
    rnd = random.Random(7)
    mono = "".join(rnd.choice("ACGT") for _ in range(300))
    copies = [(f"c{i}", mono) for i in range(5)]
    s = dz.DesignSettings(mode=dz.MODE_FISH, **dz.mode_preset(dz.MODE_FISH))
    skip = [lc.Region(0, 150, "low complexity")]
    pairs = dz.design_fish("m", mono, copies, s, avoid=[("other", mono[200:280])], skip=skip)
    assert pairs
    for p in pairs:
        pr = p.probe
        assert pr.start > 150
        assert s.tm_min <= pr.tm <= s.tm_max
        assert abs(pr.tm - tt.hyb_tm(pr.sequence, s.hyb)) < 0.01
        assert pr.extra["wash_tm"] >= s.hyb.wash_temp_c
    hit = [p for p in pairs if p.coverage["cross_hits"]]
    assert all(p.coverage["cross_names"] == ["other"] for p in hit)
    back = dz.DesignSettings.from_dict(s.to_dict())
    assert back.hyb == s.hyb
