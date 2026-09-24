"""Tree builder: input checks, split handling, rooting, reports, run folders."""

from __future__ import annotations

import shutil
import textwrap
from pathlib import Path

import pytest

from repeatviewer.backends import tree_builder as tb
from repeatviewer.job_options import (
    TREE_ENGINE_FASTME,
    TREE_ENGINE_IQTREE,
    TREE_ROOT_MIDPOINT,
    TREE_ROOT_NONE,
    TREE_ROOT_OUTGROUP,
    TreeOptions,
)


def _rec(name, seq):
    return tb.SeqRecordLite(name=name, original=name, source="x.fa", seq=seq)


# --- input handling ---------------------------------------------------------


def test_newick_safe_keeps_readable_names():
    assert tb.newick_safe("CroSat3 monomer (chr1):12-180") == "CroSat3_monomer_chr1_12-180"
    assert tb.newick_safe("  ") == "seq"


def test_read_sequences_makes_names_unique(tmp_path):
    fa = tmp_path / "a.fa"
    fa.write_text(">x one\nACGT\n>x two\nACGA\n>y\nAAAA\n")
    recs = tb.read_sequences([fa])
    assert [r.name for r in recs] == ["x", "x_2", "y"]
    assert recs[1].original == "x two"


def test_looks_like_dna():
    assert tb.looks_like_dna([_rec("a", "ACGTNNACGT-")])
    assert not tb.looks_like_dna([_rec("a", "MKVLLAQEEFWRS")])


def test_coverage_filter_excludes_fragments_relative_to_median():
    recs = [_rec("a", "A" * 100), _rec("b", "A" * 98), _rec("c", "A" * 60), _rec("d", "A" * 102)]
    kept, excluded = tb.coverage_filter(recs, 0.7)
    assert [r.name for r in kept] == ["a", "b", "d"]
    assert excluded[0][0].name == "c" and excluded[0][1] == pytest.approx(60 / 99)
    kept, excluded = tb.coverage_filter(recs, 0.0)
    assert len(kept) == 4 and not excluded


def test_coverage_filter_ignores_gaps():
    recs = [_rec("a", "AC--GT"), _rec("b", "ACGTGT"), _rec("c", "A-----")]
    kept, excluded = tb.coverage_filter(recs, 0.5)
    assert [r.name for r, _ in excluded] == ["c"]


def test_alignment_stats_counts_informative_sites():
    seqs = [
        "AAGT-",
        "AAGT-",
        "ACCTA",
        "ACCTA",
    ]
    st = tb.alignment_stats(seqs)
    # col1 constant, col2 A/A/C/C informative, col3 G/G/C/C informative,
    # col4 constant, col5 A/A + gaps -> constant
    assert (st.length, st.informative_sites, st.constant_sites) == (5, 2, 3)
    assert st.gap_fraction == pytest.approx(2 / 20)


def test_alignment_stats_singletons_are_not_informative():
    st = tb.alignment_stats(["AAAA", "AAAA", "AAAC", "AAAG"])
    assert st.informative_sites == 0


# --- splits, support, rooting -----------------------------------------------

NWK = "((A:1,B:1)90/95:1,(C:1,D:1)80/85:1,(E:1,F:1)70/75:1);"


def test_splits_are_orientation_independent():
    t1 = tb.read_newick(NWK)
    t2 = tb.read_newick("((C:1,D:1):1,((A:1,B:1):1,(E:1,F:1):1):1);")
    assert tb.rf_distance(t1, t2) == 0
    t3 = tb.read_newick("((A:1,C:1):1,(B:1,D:1):1,(E:1,F:1):1);")
    assert tb.rf_distance(t1, t3) == 4


def _label_by_split(tree):
    taxa = sorted(t.name for t in tree.get_terminals())
    ref, all_t = taxa[0], frozenset(taxa)
    out = {}
    for c in tree.get_nonterminals():
        side = frozenset(t.name for t in c.get_terminals())
        if ref in side:
            side = all_t - side
        if 1 < len(side) < len(taxa) - 1 and c.name:
            out[side] = c.name
    return out


@pytest.mark.parametrize("mode,outgroup", [
    (TREE_ROOT_OUTGROUP, ["E"]),
    (TREE_ROOT_OUTGROUP, ["E", "F"]),
    (TREE_ROOT_OUTGROUP, ["A"]),
    (TREE_ROOT_MIDPOINT, []),
])
def test_rooting_keeps_supports_on_their_splits(mode, outgroup):
    tree = tb.read_newick(NWK)
    before = _label_by_split(tree)
    sentence = tb.apply_rooting(tree, mode, outgroup)
    after = _label_by_split(tree)
    assert before == after
    assert sentence


def test_rooting_with_outgroup_puts_it_at_the_root():
    tree = tb.read_newick(NWK)
    tb.apply_rooting(tree, TREE_ROOT_OUTGROUP, ["E", "F"])
    root_children = [sorted(t.name for t in c.get_terminals()) for c in tree.root.clades]
    assert ["E", "F"] in root_children


def test_rooting_missing_outgroup_leaves_tree_unrooted():
    tree = tb.read_newick(NWK)
    assert "Unrooted" in tb.apply_rooting(tree, TREE_ROOT_OUTGROUP, ["Z"])
    assert "Unrooted" in tb.apply_rooting(tree, TREE_ROOT_NONE, [])


def test_support_from_replicates():
    main = tb.read_newick("((A,B),(C,D),(E,F));")
    reps = [tb.read_newick(x) for x in (
        "((A,B),(C,D),(E,F));",
        "((A,B),(C,E),(D,F));",
        "((A,C),(B,D),(E,F));",
        "((A,B),(C,D),(E,F));",
    )]
    tb.support_from_replicates(main, reps)
    got = _label_by_split(main)
    assert got[frozenset("AB") if "A" not in "AB" else frozenset("CDEF")] == "75"
    assert got[frozenset("CD")] == "50"
    assert got[frozenset("EF")] == "75"


# --- IQ-TREE output parsing -------------------------------------------------


def test_parse_iqtree_report_and_runtrees():
    report = textwrap.dedent("""
        Best-fit model according to BIC: TIM2e+G4
        Log-likelihood of the tree: -850.1348 (s.e. 41.2)
    """)
    assert tb.iqtree_model_from_report(report) == "TIM2e+G4"
    assert tb.iqtree_loglik_from_report(report) == pytest.approx(-850.1348)
    runs = tb.parse_runtrees("[ lh=-850.1348313 ](A,B,(C,D));\n[ lh=-851.2 ](A,C,(B,D));\n")
    assert [lh for lh, _ in runs] == [pytest.approx(-850.1348313), pytest.approx(-851.2)]


# --- output naming (Jiangzhao rule 1: keep the full input name) --------------


def test_outputs_next_to_input_keep_full_name_and_never_overwrite(tmp_path):
    fa = tmp_path / "CroSat3_monomers.fa"
    fa.write_text(">a\nA\n")
    l1 = tb.plan_outputs([fa])
    assert l1.base == "CroSat3_monomers.fa_tree"
    assert l1.path(".treefile") == tmp_path / "CroSat3_monomers.fa_tree.treefile"
    assert l1.path("_report.html").parent == tmp_path
    assert l1.path(".iqtree") == tmp_path / "CroSat3_monomers.fa_tree_details" / "CroSat3_monomers.fa_tree.iqtree"
    l1.path(".treefile").write_text("(a);")
    l2 = tb.plan_outputs([fa])
    assert l2.base == "CroSat3_monomers.fa_tree_2"
    assert tb.output_base_name([fa, tmp_path / "b.fa"]) == "CroSat3_monomers.fa_multi"


# --- methods text -------------------------------------------------------------


def _report(**kw):
    base = dict(
        input_files=["a.fa"], n_input=20, n_used=18, excluded=[("x", 0.4), ("y", 0.5)],
        dna=True, aligned_by="L-INS-i", reversed_by_mafft=["z"],
        stats=tb.AlignmentStats(18, 250, 60, 150, 0.1), engine=TREE_ENGINE_IQTREE,
        engine_version="2.3.6", mafft_version="v7.505", model="TN+F+G4", loglik=-1.0,
        runs=5, runs_same_topology=5, support="Branch support: x.",
        rooting="Rooted with outgroup o.", seed=12345,
    )
    base.update(kw)
    return tb.TreeRunReport(**base)


def test_methods_paragraph_iqtree():
    text = tb.methods_paragraph(_report(), TreeOptions())
    for piece in ("18 nucleotide sequences were aligned with MAFFT v7.505 using the L-INS-i strategy",
                  "--adjustdirectionaccurately",
                  "2 sequences shorter than 70% of the median sequence length were excluded",
                  "60 of which were parsimony-informative", "IQ-TREE 2.3.6",
                  "TN+F+G4, selected by ModelFinder according to the Bayesian information criterion",
                  "from 5 independent searches (random seeds 12345–12349) was retained",
                  "1000 SH-aLRT replicates and 1000 ultrafast bootstrap replicates "
                  "(UFBoot2, with the -bnni correction)",
                  "The tree was rooted with o as outgroup."):
        assert piece in text, piece
    assert "(s)" not in text and ")(" not in text and ") (" not in text


def test_methods_paragraph_fastme_and_citations():
    r = _report(engine=TREE_ENGINE_FASTME, aligned_by="", excluded=[])
    opts = TreeOptions(engine=TREE_ENGINE_FASTME)
    text = tb.methods_paragraph(r, opts)
    assert "provided alignment" in text and "FastME 2.3.6" in text and "TN93" in text
    assert "100 bootstrap replicates" in text and "The tree was rooted with o as outgroup." in text
    cites = tb.citations_for(r, opts)
    assert any("FastME 2.0" in c for c in cites) and not any("MAFFT" in c for c in cites)


# --- full run with a fake engine (no external tree program needed) ------------


FAKE_FASTME = """#!/usr/bin/env python3
import sys, re
args = dict(a.split("=", 1) for a in sys.argv[1:] if a.startswith("--") and "=" in a)
if "--version" in sys.argv:
    print("FastME 2.1.6.1"); sys.exit(0)
names = [l.split()[0] for l in open(args["--input_data"]).read().splitlines()[1:]]
def nwk(ns):
    return "(" + ",".join(f"{n}:0.1" for n in ns[:2]) + ",(" + ",".join(f"{n}:0.1" for n in ns[2:]) + "):0.1);"
open(args["--output_tree"], "w").write(nwk(names) + "\\n")
if "-B" in sys.argv:
    out = sys.argv[sys.argv.index("-B") + 1]
    open(out, "w").write("\\n".join(nwk(names) for _ in range(4)) + "\\n")
"""


def test_build_tree_run_end_to_end_with_fake_fastme(tmp_path, monkeypatch):
    fake = tmp_path / "bin" / "fastme"
    fake.parent.mkdir()
    fake.write_text(FAKE_FASTME)
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", f"{fake.parent}:" + __import__("os").environ["PATH"])
    fa = tmp_path / "aln.fa"
    fa.write_text(">s1\nACGTACGTAA\n>s2\nACGTACGTTA\n>s3\nACGAACGTTT\n>s4\nACGAACCTTT\n>s5\nTCGAACCTTT\n")
    opts = TreeOptions(engine=TREE_ENGINE_FASTME, use_system_path=True, threads=1,
                       distance_bootstrap=4, root_mode=TREE_ROOT_MIDPOINT)
    details, files, report = tb.build_tree_run([fa], opts, stamp="t")
    assert details == tmp_path / "aln.fa_tree_details"
    for suffix in (".treefile", ".aln.fa", ".svg", ".pdf", ".png", "_report.html"):
        assert (tmp_path / f"aln.fa_tree{suffix}") in files  # next to the input
    for suffix in ("_methods.txt", "_commands.txt", ".names.tsv"):
        assert (details / f"aln.fa_tree{suffix}") in files
    tree_text = (tmp_path / "aln.fa_tree.treefile").read_text()
    assert "s1" in tree_text and "S00001" not in tree_text  # safe IDs mapped back
    assert "100" in tree_text  # bootstrap support computed from the replicate file
    html = report.read_text()
    assert "provided alignment" in (details / "aln.fa_tree_methods.txt").read_text()
    assert "<svg" in html and "Midpoint-rooted" in html
    assert not (details / ".work").exists()


def test_missing_engine_is_an_error_not_a_downgrade(tmp_path, monkeypatch):
    monkeypatch.setenv("PATH", str(tmp_path))
    fa = tmp_path / "a.fa"
    fa.write_text("".join(f">s{i}\nACGT\n" for i in range(5)))
    with pytest.raises(tb.TreeBuildError, match="conda install -c bioconda iqtree"):
        tb.build_tree_run([fa], TreeOptions(use_system_path=True), stamp="t")
    assert not any(tmp_path.glob("a.fa_tree*"))


def test_too_few_sequences(tmp_path, monkeypatch):
    fake = tmp_path / "iqtree2"
    fake.write_text("#!/bin/sh\n")
    fake.chmod(0o755)
    monkeypatch.setenv("PATH", str(tmp_path))
    fa = tmp_path / "a.fa"
    fa.write_text(">a\nACGT\n>b\nACGT\n>c\nACGT\n")
    with pytest.raises(tb.TreeBuildError, match="at least 4"):
        tb.build_tree_run([fa], TreeOptions(use_system_path=True), stamp="t")


@pytest.mark.skipif(not (shutil.which("iqtree2") and shutil.which("mafft")),
                    reason="needs iqtree2 and mafft on PATH")
def test_build_tree_run_real_iqtree(tmp_path):
    import random
    rng = random.Random(1)
    base = "".join(rng.choice("ACGT") for _ in range(160))
    seqs = []
    for i in range(8):
        s = list(base)
        for _ in range(8 + i * 3):
            j = rng.randrange(len(s))
            s[j] = rng.choice("ACGT")
        seqs.append("".join(s))
    fa = tmp_path / "sat.fa"
    fa.write_text("".join(f">m{i}\n{s}\n" for i, s in enumerate(seqs)))
    opts = TreeOptions(use_system_path=True, use_system_mafft=True, threads=1, runs=2,
                       root_mode=TREE_ROOT_OUTGROUP, outgroup=["m7"])
    details, files, report = tb.build_tree_run([fa], opts, stamp="t")
    assert (details / "sat.fa_tree.iqtree").is_file()
    assert "of 2 found the same topology" in report.read_text()


def test_input_summary(tmp_path):
    fa = tmp_path / "a.fa"
    fa.write_text(">a\nAAGT\n>b\nAAGT\n>c\nACCT\n>d\nACCA\n")
    s = tb.input_summary([fa], 0.7)
    assert (s.n_sequences, s.aligned, s.length, s.informative_sites, s.fragments) == (4, True, 4, 2, 0)
    fb = tmp_path / "b.fa"
    fb.write_text(">a\nAAGTAAGT\n>b\nAAGTAAG\n>c\nAAG\n")
    s = tb.input_summary([fb], 0.7)
    assert (s.aligned, s.length, s.fragments) == (False, None, 1)
    assert tb.input_summary([tmp_path / "missing.fa"], 0.7).error


def test_fragment_count_matches_coverage_filter():
    lengths = [100, 98, 60, 102, 30]
    recs = [_rec(str(i), "A" * n) for i, n in enumerate(lengths)]
    for cov in (0.0, 0.5, 0.7, 0.99):
        assert tb.fragment_count(lengths, cov) == len(tb.coverage_filter(recs, cov)[1])


def test_root_branch_support_shown_once():
    tree = tb.read_newick("((A:1,B:1)90/95:1,(C:1,D:1)80/85:1,(E:1,F:1)70/75:1);")
    tb.apply_rooting(tree, TREE_ROOT_MIDPOINT, [])
    labels = [c.name for c in tree.get_nonterminals() if c.name]
    assert len(labels) == len(set(labels)) == 3



# --- round 2: crash handling, separate searches, length reference -------------


def test_coverage_filter_longest_reference():
    recs = [_rec("a", "A" * 800), _rec("b", "A" * 200), _rec("c", "A" * 150),
            _rec("d", "A" * 120), _rec("e", "A" * 790)]
    kept, excluded = tb.coverage_filter(recs, 0.7, "median")
    # The median (200 bp) is itself a fragment length, so only d (120 bp) goes.
    assert [r.name for r, _ in excluded] == ["d"]
    kept, excluded = tb.coverage_filter(recs, 0.7, "longest")
    assert [r.name for r, _ in excluded] == ["b", "c", "d"]
    assert tb.fragment_count([r.seq.__len__() for r in recs], 0.7, "longest") == 3


def test_iqtree_thread_args():
    assert tb.iqtree_thread_args(1) == ["-T", "1"]
    assert tb.iqtree_thread_args(4) == ["-T", "AUTO", "--threads-max", "4"]


def test_summarize_tool_error_openmp_and_stack_traces():
    omp = ("Run 5 gave best log-likelihood\nOMP: Error #34: System unable to allocate "
           "necessary resources for OMP thread:\nERROR: STACK TRACE FOR DEBUGGING:\n"
           "ERROR: 2   _sigtramp()\nERROR: *** IQ-TREE CRASHES WITH SIGNAL ABORTED")
    msg = tb.summarize_tool_error("IQ-TREE", -6, omp)
    assert "could not create more threads" in msg and "_sigtramp" not in msg
    msg = tb.summarize_tool_error("IQ-TREE", 2, "ERROR: Sequence S00003 contains too many characters\n")
    assert msg.endswith("ERROR: Sequence S00003 contains too many characters")


FAKE_IQTREE = """#!/usr/bin/env python3
import sys, os
a = sys.argv
prefix = a[a.index("--prefix") + 1]
seed = int(a[a.index("--seed") + 1])
model = a[a.index("-m") + 1]
open(os.environ["FAKE_LOG"], "a").write(" ".join(a[1:]) + "\\n")
if "--version" in a:
    print("IQ-TREE version 3.1.3"); sys.exit(0)
if os.environ.get("FAKE_CRASH") == str(seed):
    print("OMP: Error #34: System unable to allocate necessary resources for OMP thread")
    sys.exit(134)
names = [l[1:].strip() for l in open(a[a.index("-s") + 1]) if l.startswith(">")]
tree = "(" + ",".join(names[:2]) + ",(" + ",".join(names[2:]) + "));"
if seed % 2 == 0:
    tree = "(" + names[0] + "," + names[2] + ",(" + ",".join([names[1]] + names[3:]) + "));"
open(prefix + ".treefile", "w").write(tree + "\\n")
best = "TN+F+G4" if model == "MFP" else model
open(prefix + ".iqtree", "w").write(
    f"Best-fit model according to BIC: {best}\\nLog-likelihood of the tree: {-100 - seed % 7}.5\\n")
open(prefix + ".log", "w").write("log " + str(seed))
"""


def _fake_iqtree(tmp_path, monkeypatch, crash_seed=None):
    fake = tmp_path / "bin" / "iqtree3"
    fake.parent.mkdir(exist_ok=True)
    fake.write_text(FAKE_IQTREE.replace("\\\\n", "\\n"))
    fake.chmod(0o755)
    log = tmp_path / "calls.txt"
    monkeypatch.setenv("FAKE_LOG", str(log))
    if crash_seed is not None:
        monkeypatch.setenv("FAKE_CRASH", str(crash_seed))
    monkeypatch.setenv("PATH", f"{fake.parent}:" + __import__("os").environ["PATH"])
    fa = tmp_path / "aln.fa"
    fa.write_text(">s1\nACGTACGTAA\n>s2\nACGTACGTTA\n>s3\nACGAACGTTT\n"
                  ">s4\nACGAACCTTT\n>s5\nTCGAACCTTT\n")
    return fa, log


def test_independent_searches_are_separate_processes(tmp_path, monkeypatch):
    fa, log = _fake_iqtree(tmp_path, monkeypatch)
    opts = TreeOptions(use_system_path=True, threads=4, runs=3, seed=10)
    run_dir, files, report = tb.build_tree_run([fa], opts, stamp="t")
    calls = [c for c in log.read_text().splitlines() if "--prefix" in c]
    assert len(calls) == 3
    assert "-m MFP" in calls[0] and "-m TN+F+G4" in calls[1] and "-m TN+F+G4" in calls[2]
    assert all("--runs" not in c for c in calls)
    assert all("-T AUTO --threads-max 4" in c for c in calls)
    assert [c.split("--seed ")[1].split()[0] for c in calls] == ["10", "11", "12"]
    html = report.read_text()
    assert "of 3 found the same topology" in html
    assert (run_dir / "aln.fa_tree.iqtree").is_file()
    assert report.parent == tmp_path


def test_failed_run_writes_log_next_to_input(tmp_path, monkeypatch):
    fa, _log = _fake_iqtree(tmp_path, monkeypatch, crash_seed=12345)
    opts = TreeOptions(use_system_path=True, threads=2, runs=1)
    with pytest.raises(tb.TreeBuildError) as err:
        tb.build_tree_run([fa], opts, stamp="t")
    msg = str(err.value)
    assert "could not create more threads" in msg
    failed = tmp_path / "aln.fa_tree_failed.log"
    assert failed.is_file() and f"Full log: {failed}" in msg
    text = failed.read_text()
    assert "iqtree3" in text and "Reason:" in text
    assert not any(tmp_path.glob("aln.fa_tree.*")) and not any(tmp_path.glob("aln.fa_tree_*[!g]"))


def test_one_failed_search_is_skipped_with_warning(tmp_path, monkeypatch):
    fa, _log = _fake_iqtree(tmp_path, monkeypatch, crash_seed=11)
    opts = TreeOptions(use_system_path=True, threads=1, runs=3, seed=10)
    run_dir, files, report = tb.build_tree_run([fa], opts, stamp="t")
    html = report.read_text()
    assert "1 of 3 IQ-TREE searches failed" in html



def test_iqtree_citation_follows_installed_version():
    opts = TreeOptions()
    c3 = tb.citations_for(_report(engine_version="3.1.3"), opts)
    c2 = tb.citations_for(_report(engine_version="2.3.6"), opts)
    assert any("IQ-TREE 3" in c for c in c3) and not any("IQ-TREE 2" in c for c in c3)
    assert any("IQ-TREE 2" in c for c in c2)


def test_methods_single_search_and_unrooted():
    text = tb.methods_paragraph(_report(runs=1, rooting="Unrooted (x)."), TreeOptions(runs=1))
    assert "A single tree search was run (random seed 12345)." in text
    assert text.endswith("The tree is unrooted.")



def test_outgroup_from_another_file_is_reported(tmp_path, monkeypatch):
    fa, _log = _fake_iqtree(tmp_path, monkeypatch)
    opts = TreeOptions(use_system_path=True, threads=1, runs=1,
                       root_mode=TREE_ROOT_OUTGROUP, outgroup=["CroSat9"])
    _details, _files, report = tb.build_tree_run([fa], opts, stamp="t")
    html = report.read_text()
    assert "not in this input: CroSat9" in html and "Unrooted" in html
