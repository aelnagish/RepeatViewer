"""RT-domain mode (Tools → Tree, "RT domain from DNA") with a small synthetic REXdb."""

from __future__ import annotations

import random
import shutil
from itertools import product

import pytest

from repeatviewer.backends import rt_domains as rd
from repeatviewer.backends import tree_builder as tb
from repeatviewer.backends import tree_traits as tt
from repeatviewer.job_options import TREE_SEQ_RT, TreeOptions

AAS = "ACDEFGHIKLMNPQRSTVWY"
_CODONS: dict[str, list[str]] = {}
for _c, _a in zip(("".join(p) for p in product("TCAG", repeat=3)),
                  "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG"):
    _CODONS.setdefault(_a, []).append(_c)

LINEAGES = {"Ty1-RT": ("Ale", "SIRE"), "Ty3-RT": ("Tekay", "Athila")}
SUPERFAMILY = {"Ty1-RT": "Ty1/copia", "Ty3-RT": "Ty3/gypsy"}

needs_mafft = pytest.mark.skipif(shutil.which("mafft") is None, reason="MAFFT not installed")


def _mutate(rng, seq, rate):
    return "".join(rng.choice(AAS) if rng.random() < rate else c for c in seq)


def _rc(s):
    return s[::-1].translate(str.maketrans("ACGT", "TGCA"))


@pytest.fixture(scope="module")
def world(tmp_path_factory):
    rng = random.Random(3)
    root = tmp_path_factory.mktemp("rexdb")
    fasta, table = root / "Mini_v1.fasta", root / "Mini_v1.classification"
    ancestors = {}
    with fasta.open("w") as fa, table.open("w") as cl:
        n = 0
        for rt_type, lineages in LINEAGES.items():
            base = "".join(rng.choice(AAS) for _ in range(220))
            for lin in lineages:
                anc = _mutate(rng, base, 0.35)
                ancestors[lin] = anc
                for _k in range(8):
                    n += 1
                    rid = f"REXdb_ID{n}"
                    fa.write(f">{rt_type}__{rid}\n{_mutate(rng, anc, 0.08)}\n")
                    cl.write(f"{rid}\tClass_I\tLTR\t{SUPERFAMILY[rt_type]}\t{lin}\n")
    records, truth = [], {}
    for lin, anc in ancestors.items():
        for k in range(3):
            prot = _mutate(rng, anc, 0.12)
            dna = "".join(rng.choice(_CODONS[a]) for a in prot)
            flank = lambda m: "".join(rng.choice("ACGT") for _ in range(m))  # noqa: E731
            dna = flank(300) + dna + flank(250)
            if k == 1:
                dna = _rc(dna)
            records.append((f"{lin}_{k}", dna))
            truth[f"{lin}_{k}"] = lin
    records.append(("decoy", "".join(rng.choice("ACGT") for _ in range(1500))))
    return fasta, records, truth, root


def test_load_rexdb_and_errors(world, tmp_path):
    fasta, *_ = world
    db = rd.load_rexdb(fasta)
    assert len(db.rt) == 32
    rid = next(iter(db.rt))
    assert db.rt[rid][0] == "Ty1-RT" and db.lineage(rid) == "Ale"
    assert db.superfamily(rid) == "Ty1/copia"
    assert rd.looks_like_rexdb(fasta)
    lone = tmp_path / "x.fasta"
    lone.write_text(">a\nMK\n")
    assert not rd.looks_like_rexdb(lone)
    with pytest.raises(rd.RexdbError, match="classification"):
        rd.load_rexdb(lone)


def test_six_frames_and_dna_span():
    dna = "ATGAAACCCTAA" + "GG"
    frames = {(s, f): p for s, f, p in rd.six_frames(dna)}
    assert frames[("+", 0)].startswith("MKP")
    assert "*" not in "".join(frames.values())
    start, end = rd.dna_span(len(dna), "+", 0, 1, 3)
    assert (start, end) == (1, 9)


@needs_mafft
def test_rt_domains_found_and_classified(world, tmp_path):
    fasta, records, truth, _root = world
    mafft = (["mafft"], None)
    res = rd.rt_domains_for(records, fasta, mafft=mafft, cache_dir=tmp_path / "cache",
                            per_lineage=2)
    assert set(res.hits) == set(truth)
    assert res.missing == ["decoy"]
    assert all(res.hits[n].lineage == lin for n, lin in truth.items())
    assert {h.strand for h in res.hits.values()} == {"+", "-"}
    assert len(res.references) == 2 * 4
    assert "no RT domain found" in res.table_tsv()
    # cached profiles: a second call builds nothing
    built = []
    rd.rt_domains_for(records, fasta, mafft=mafft, cache_dir=tmp_path / "cache",
                      per_lineage=0, say=built.append)
    assert not any("Building" in m for m in built)


@needs_mafft
def test_tree_builder_rt_records_and_traits(world, tmp_path):
    fasta, records, truth, _root = world
    recs = [tb.SeqRecordLite(n, f"chr1:{i}-{i + 9}(+) {n}", "copies.fa", s)
            for i, (n, s) in enumerate(records, 1)]
    opts = TreeOptions(sequence_mode=TREE_SEQ_RT, rexdb_fasta=str(fasta), rexdb_references=1,
                       use_system_mafft=True, methylation_mask="CG")
    warnings: list[str] = []
    out, res = tb.rt_domain_records(recs, opts, lambda _m: None, None,
                                    warnings)
    names = [r.name for r in out]
    assert "decoy" not in names and len(names) == len(truth) + 4
    assert sum(r.source == "REXdb" for r in out) == 4
    assert all(set(r.seq) <= set(AAS + "X") for r in out)
    assert any("No RT domain" in w for w in warnings)
    assert any("CG/CHG" in w for w in warnings)
    assert tb.run_tag(opts, dna=True) == "_rt"
    tree = tmp_path / "copies.fa_tree_rt.treefile"
    tree.write_text("();\n")
    table = tt.read_trait_table(tb._write_rt_traits(tree, res, out))
    assert table.rows["Ale_0"]["lineage"] == "Ale"
    ref = next(r.name for r in out if r.source == "REXdb")
    assert table.rows[ref]["kind"] == "REXdb" and table.rows[ref]["rexdb_lineage"]
    row, methods = tb._rt_texts(res, opts, len(records))
    assert "12 of 13" in row and "Mini_v1" in methods


def test_rt_mode_needs_rexdb(tmp_path):
    recs = [tb.SeqRecordLite(f"s{i}", f"s{i}", "a.fa", "ACGT" * 30) for i in range(4)]
    with pytest.raises(tb.TreeBuildError, match="needs REXdb"):
        tb.rt_domain_records(recs, TreeOptions(sequence_mode=TREE_SEQ_RT), lambda _m: None,
                             None, [])
