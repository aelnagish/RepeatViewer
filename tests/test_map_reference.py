"""Map to Reference backend: parsers, projection, consensus, variants, lift-over,
multiple placements, output layout, and real runs when the engines are present."""

from __future__ import annotations

import json
import random
import shutil
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from repeatviewer.backends import map_reference as mr
from repeatviewer.backends.sequence_features import SequenceFeature, sidecar_path
from repeatviewer.job_options import (
    MAP_ENGINE_BOWTIE2,
    MAP_ENGINE_MAFFT,
    MAP_ENGINE_MINIMAP2,
    MAP_MULTI_ALL,
    MAP_MULTI_NONE,
    MapOptions,
)


def _read(row, start=0, ref_row="", ins=(), **kw):
    r = mr.MappedRead(name=kw.pop("name", "r"), ref_id="ref", start=start, end=start + len(row),
                      strand=kw.pop("strand", "+"), row=row, insertions=list(ins),
                      ref_row=ref_row, **kw)
    return mr.compute_stats(r)


# --- parsers ---------------------------------------------------------------


def test_project_gapped_pair_insertion_deletion_and_clipping():
    ref = "AC--GTACGT-"
    qry = "TCAAG-ACGTA"
    start, end, row, ref_row, ins, clipped, clip_left = mr.project_gapped_pair(ref, qry)
    assert (start, end) == (0, 8)
    assert row == "TCG-ACGT"
    assert ref_row == "ACGTACGT"
    assert ins == [(2, "AA")]
    assert clipped == 1 and clip_left == 0          # the trailing A has no reference column


def test_project_gapped_pair_no_shared_column():
    assert mr.project_gapped_pair("ACGT----", "----ACGT") is None


def test_parse_cs_long_form():
    row, ref, ins = mr.parse_cs("=ACG*ct+tt=GG-aa=T", 10)
    assert row == "ACGTGG--T"
    assert ref == "ACGCGGAAT"
    assert ins == [(14, "TT")]


def test_parse_paf_line_reverse_strand_and_stats():
    line = ("q1\t12\t1\t11\t-\tchr\t100\t20\t30\t9\t10\t60\ttp:A:P\tcs:Z:=ACGT*ag=ACGTA\n")
    name, qlen, read = mr.parse_paf_line(line)
    assert (name, qlen) == ("q1", 12)
    assert (read.start, read.end, read.strand) == (20, 30, "-")
    assert read.mismatches == 1 and read.matches == 9
    assert read.clip_left == 1                       # qlen - qend on the minus strand
    assert mr.parse_paf_line("q2\t50\t0\t0\t*\t*\t0\t0\t0\t0\t0\t0\n")[2] is None


def test_parse_sam_line_cigar_md_and_softclip():
    # read: 2 soft-clipped, 4 M (one mismatch), 1 I, 3 M, 2 D, 2 M
    seq = "TTACGTAGCAGG"
    line = "\t".join(["r1", "16", "chr", "101", "42", "2S4M1I3M2D2M", "*", "0", "0", seq, "*",
                      "AS:i:10", "MD:Z:2C4^TT2"]) + "\n"
    name, qlen, read = mr.parse_sam_line(line)
    assert read.strand == "-" and read.start == 100
    assert read.row == "ACGTGCA--GG"
    assert read.insertions == [(104, "A")]
    assert read.ref_row == "ACCTGCATTGG"
    assert read.mismatches == 1 and read.deleted == 2 and read.inserted == 1
    assert read.clip_left == 2
    unmapped = "\t".join(["r2", "4", "*", "0", "0", "*", "*", "0", "0", "ACGT", "*"]) + "\n"
    assert mr.parse_sam_line(unmapped)[2] is None


def test_parse_mafft_alignment_reversed_names():
    pairs = [("REF", "ACGTACGTAC"), ("Q000001", "--GTAC----"), ("_R_Q000002", "ACGTTCG---"),
             ("Q000003", "----------")]
    reads, unmapped = mr.parse_mafft_alignment(pairs, "REF", "ref",
                                               {"Q000001": "a", "Q000002": "b", "Q000003": "c"})
    assert [(r.name, r.strand, r.start, r.end) for r in reads] == [("a", "+", 2, 6),
                                                                   ("b", "-", 0, 7)]
    assert reads[1].mismatches == 1
    assert unmapped == ["c"]


# --- filtering and multiple placements --------------------------------------


def test_placement_verdicts():
    opts = MapOptions(min_identity=0.8, min_aligned_bp=5)
    good = _read("ACGTACGTAC", ref_row="ACGTACGTAC", mapq=40)
    assert mr.placement_verdict(good, opts, MAP_ENGINE_MINIMAP2) is None
    low = _read("ACGTTTTTTT", ref_row="ACGTACGTAC")
    assert "identity" in mr.placement_verdict(low, opts, MAP_ENGINE_MINIMAP2)
    short = _read("ACG", ref_row="ACG")
    assert "shorter" in mr.placement_verdict(short, opts, MAP_ENGINE_MINIMAP2)
    sec = _read("ACGTACGTAC", ref_row="ACGTACGTAC", kind="secondary")
    assert mr.placement_verdict(sec, opts, MAP_ENGINE_MINIMAP2) == ""
    assert mr.placement_verdict(sec, MapOptions(multiple_hits=MAP_MULTI_ALL, min_aligned_bp=5),
                                MAP_ENGINE_MINIMAP2) is None
    amb = _read("ACGTACGTAC", ref_row="ACGTACGTAC", mapq=0)
    none = MapOptions(multiple_hits=MAP_MULTI_NONE, min_aligned_bp=5)
    assert "ambiguous" in mr.placement_verdict(amb, none, MAP_ENGINE_MINIMAP2)
    assert mr.placement_verdict(amb, none, MAP_ENGINE_MAFFT) is None


def test_filter_placements_collects_unused():
    opts = MapOptions(min_aligned_bp=5)
    kept = _read("ACGTACGTAC", ref_row="ACGTACGTAC", name="a")
    res = mr.filter_placements([kept], ["a", "b"], opts, MAP_ENGINE_MINIMAP2, unmapped=["c"])
    assert [r.name for r in res.kept] == ["a"]
    assert set(res.unused) == {"b", "c"}


# --- consensus and variants -------------------------------------------------


def _counts(ref_len, reads):
    counts = np.zeros((ref_len, 6), dtype=np.int32)
    ins = {}
    for r in reads:
        mr.add_to_counts(counts, r)
        for p, s in r.insertions:
            ins.setdefault(p, Counter())[s] += 1
    return counts, ins


def test_consensus_majority_iupac_deletion_and_no_coverage():
    ref = "ACGTACGTAA"
    reads = [_read("ACGTA-GT", ref_row=ref[:8]), _read("ACCTA-GT", ref_row=ref[:8]),
             _read("ACCTACGT", ref_row=ref[:8], ins=[(4, "GG")]),
             _read("ACGTA-GT", ref_row=ref[:8], ins=[(4, "GG")])]
    counts, ins = _counts(len(ref), reads)
    cons, m = mr.build_consensus(ref, counts, ins, threshold=0.0, no_coverage="N")
    # position 2: G/C tie -> S; position 5 deleted by 3 of 4; GG inserted by 2 of 4 (50 %)
    assert cons == "ACSTGGAGTNN"
    assert m[5] == -1 and m[4] == 6
    cons_ref, _ = mr.build_consensus(ref, counts, ins, no_coverage="reference")
    assert cons_ref.endswith("AA")
    cons_skip, _ = mr.build_consensus(ref, counts, ins, no_coverage="skip")
    assert cons_skip == "ACSTGGAGT"


def test_consensus_threshold_gives_iupac():
    counts = np.zeros((1, 6), dtype=np.int32)
    counts[0, 0], counts[0, 2], counts[0, 1] = 6, 3, 1          # A 60 %, G 30 %, C 10 %
    assert mr.consensus_calls(counts, 0.0) == ["A"]
    assert mr.consensus_calls(counts, 0.75) == ["R"]            # A + G reach 90 %
    assert mr.consensus_calls(counts, 0.95) == ["V"]            # A + C + G


def test_call_variants_with_context_and_features():
    ref = "ACGTTCGAAA"
    reads = [_read("ATGTTCGAAA", ref_row=ref), _read("ACGTTCGAAA", ref_row=ref),
             _read("ACGTT-GAAA", ref_row=ref, ins=[(8, "C")])]
    counts, ins = _counts(len(ref), reads)
    feats = [SequenceFeature(start=5, end=8, name="motif")]
    vs = mr.call_variants(ref, counts, ins, feats)
    got = {(v.pos, v.alt): v for v in vs}
    assert got[(2, "T")].context == "CG" and got[(2, "T")].count == 1
    assert got[(6, "-")].features == "motif" and got[(6, "-")].context == "CG"
    assert got[(9, "+C")].depth == 3
    assert mr.call_variants(ref, counts, ins, feats, min_freq=0.5) == []


# --- lift-over ---------------------------------------------------------------


def test_lift_features_status_and_coordinates():
    ref = "AAAACCCCGGGGTTTT"
    feats = [SequenceFeature(start=1, end=4, name="A"), SequenceFeature(start=5, end=8, name="C"),
             SequenceFeature(start=9, end=12, name="G"), SequenceFeature(start=13, end=16, name="T")]
    # covers reference 2..15: A partial, C one mismatch + 2 inserted bases,
    # G three of four bases deleted, T partial
    row = "AAA" + "CCTC" + "G---" + "TTT"
    read = _read(row, start=1, ref_row=ref[1:15], ins=[(6, "XX")], clip_left=3)
    got = {feats[lf.index].name: lf for lf in mr.lift_features(read, feats)}
    assert got["A"].status == mr.STATUS_PARTIAL
    assert got["C"].status == mr.STATUS_CHANGED and got["C"].differences == 3
    assert got["G"].status == mr.STATUS_LOST
    assert got["T"].status == mr.STATUS_PARTIAL
    # coordinates on the whole sequence: 3 clipped bases, then AAA, then C…
    assert (got["C"].start, got["C"].end) == (3 + 4, 3 + 9)
    intact = _read(ref[4:8], start=4, ref_row=ref[4:8])
    only = mr.lift_features(intact, feats)
    assert [(feats[x.index].name, x.status, x.start, x.end) for x in only] == [
        ("C", mr.STATUS_INTACT, 1, 4)]


def test_feature_index_finds_long_overlapping_features():
    feats = [SequenceFeature(start=1, end=1000, name="array"),
             SequenceFeature(start=500, end=510, name="motif"),
             SequenceFeature(start=2000, end=2010, name="far")]
    idx = mr.FeatureIndex(feats)
    assert [feats[i].name for i in idx.overlapping(600, 700)] == ["array", "motif"]


def test_rle_roundtrip():
    v = np.array([0, 0, 3, 3, 3, 1])
    assert mr.rle(v) == [[0, 2], [3, 3], [1, 1]]
    assert mr.unrle(mr.rle(v)).tolist() == v.tolist()


# --- engine choice and naming ------------------------------------------------


def test_choose_engine():
    small = mr.InputInfo(n=6, longest=400, lengths_sample=[150, 300, 400])
    assert mr.choose_engine({"sat": 828}, small) == MAP_ENGINE_MAFFT
    short = mr.InputInfo(n=10_000, longest=150, lengths_sample=[150] * 100)
    assert mr.choose_engine({"chr1": 10**7}, short) == MAP_ENGINE_BOWTIE2
    long = mr.InputInfo(n=1000, longest=30_000, lengths_sample=[15_000] * 100)
    assert mr.choose_engine({"chr1": 10**7, "chr2": 10**7}, long) == MAP_ENGINE_MINIMAP2
    assert mr.auto_minimap2_preset(short) == "sr"
    assert mr.auto_minimap2_preset(long) == "map-ont"


def test_plan_outputs_keeps_full_name_and_never_overwrites(tmp_path):
    q = tmp_path / "variant2.fa"
    q.write_text(">a\nACGT\n")
    first = mr.plan_outputs([q])
    assert first.rvmap.name == "variant2.fa_map.rvmap"
    assert first.details.name == "variant2.fa_map_details"
    first.rvmap.write_text("{}")
    second = mr.plan_outputs([q])
    assert second.base == "variant2.fa_map_2"


def test_minimap2_and_bowtie2_arguments():
    opts = MapOptions(threads=3, multiple_hits=MAP_MULTI_ALL)
    args = mr.minimap2_arguments(opts, "map-ont")
    assert "--cs=long" in args and "-N" in args and "--paf-no-hit" in args
    assert "--secondary=no" in mr.minimap2_arguments(MapOptions(), "sr")
    b = mr.bowtie2_arguments(opts)
    assert b[0] == "--sensitive-local" and "-k" in b


# --- real runs ----------------------------------------------------------------


def _crosat(tmp_path: Path):
    rng = random.Random(3)

    def rnd(n):
        return "".join(rng.choice("ACGT") for _ in range(n))

    def mut(s, rate):
        return "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < rate else c
                       for c in s)

    unit = rnd(118)
    ref = "".join(mut(unit, 0.08) for _ in range(7))
    (tmp_path / "ref.fa").write_text(">CroSat1\n" + ref + "\n")
    (tmp_path / "ref.RV.gff3").write_text(
        "##gff-version 3\n"
        "CroSat1\tRepeatViewer\tmotif\t395\t410\t.\t+\t.\tID=Motif_F;Name=Motif_F;type=motif\n"
        "CroSat1\tRepeatViewer\tmotif\t150\t165\t.\t+\t.\tID=Motif_E;Name=Motif_E;type=motif\n")
    variants = {
        "v1": mut(ref[0:420], 0.02),
        "v2": mr.revcomp(mut(ref[50:200], 0.02)),
        "v3": mut(ref[230:400] + ref[430:620], 0.02),
        "junk": rnd(200),
    }
    (tmp_path / "variants.fa").write_text("".join(f">{k}\n{v}\n" for k, v in variants.items()))
    return tmp_path / "ref.fa", tmp_path / "variants.fa"


def _mafft_options(ref: Path) -> MapOptions | None:
    from repeatviewer.backends.sequence_aligner import mafft_command

    if mafft_command(use_system_mafft=False) is not None:
        return MapOptions(reference_fasta=str(ref), engine=MAP_ENGINE_MAFFT, threads=2)
    if mafft_command(use_system_mafft=True) is not None:
        return MapOptions(reference_fasta=str(ref), engine=MAP_ENGINE_MAFFT, threads=2,
                          use_system_path=True)
    return None


def test_mafft_run_end_to_end(tmp_path):
    ref, variants = _crosat(tmp_path)
    opts = _mafft_options(ref)
    if opts is None:
        pytest.skip("MAFFT not available")
    rvmap, files, report = mr.run_mapping([variants], opts)
    assert rvmap.name == "variants.fa_map.rvmap" and rvmap.parent == tmp_path
    data = mr.load_rvmap(rvmap)
    assert data["summary"]["n_placed"] == 3 and data["summary"]["n_reversed"] == 1
    contig = data["contigs"][0]
    reads = {r["name"]: r for r in contig["reads"]}
    assert reads["v2"]["strand"] == "-"
    motif_f = next(i for i, f in enumerate(contig["features"]) if f["name"] == "Motif_F")
    status_v3 = {x[0]: x[1] for x in reads["v3"]["lift"]}
    assert status_v3[motif_f] == mr.STATUS_LOST                # the 30-bp deletion
    assert [u[0] for u in data["unused"]] == ["junk"]
    details = tmp_path / data["details_dir"]
    for key in ("consensus", "mapped", "unused", "variants", "features", "reads", "report"):
        assert (details / data["details"][key]).is_file(), key
    mapped = (details / data["details"]["mapped"]).read_text()
    assert ">v2_reversed " in mapped
    # the placed sequences carry the lifted annotations (seqmap JSON sidecar)
    mapped_sidecar = sidecar_path(details / data["details"]["mapped"])
    assert mapped_sidecar.is_file()
    consensus_sidecar = sidecar_path(details / data["details"]["consensus"])
    assert "Motif_E" in consensus_sidecar.read_text(encoding="utf-8")
    assert report.is_file() and "Methods text" in report.read_text()


def test_missing_reference_is_a_clear_error(tmp_path):
    q = tmp_path / "q.fa"
    q.write_text(">a\nACGT\n")
    with pytest.raises(mr.MapReferenceError, match="reference"):
        mr.run_mapping([q], MapOptions(reference_fasta=""))
    with pytest.raises(mr.MapReferenceError, match="not found"):
        mr.run_mapping([q], MapOptions(reference_fasta=str(tmp_path / "none.fa")))


def test_mafft_refuses_multi_sequence_reference(tmp_path):
    ref = tmp_path / "ref.fa"
    ref.write_text(">a\nACGTACGTAC\n>b\nACGTACGTAC\n")
    q = tmp_path / "q.fa"
    q.write_text(">x\nACGTACG\n")
    with pytest.raises(mr.MapReferenceError, match="one reference sequence"):
        mr.run_mapping([q], MapOptions(reference_fasta=str(ref), engine=MAP_ENGINE_MAFFT))
    assert not list(tmp_path.glob("*_map_details"))           # nothing half-written


@pytest.mark.parametrize("engine", [MAP_ENGINE_MINIMAP2, MAP_ENGINE_BOWTIE2])
def test_read_mappers_end_to_end(tmp_path, engine):
    exe = {"minimap2": "minimap2", "bowtie2": "bowtie2"}[engine]
    if shutil.which(exe) is None:
        pytest.skip(f"{exe} not on PATH")
    rng = random.Random(5)
    genome = "".join(rng.choice("ACGT") for _ in range(20_000))
    (tmp_path / "g.fa").write_text(">chr1\n" + genome + "\n")
    reads = []
    for i in range(800):
        p = rng.randrange(0, len(genome) - 150)
        s = genome[p:p + 150]
        reads.append((f"r{i}", mr.revcomp(s) if i % 2 else s))
    reads.append(("nohit", "".join(rng.choice("ACGT") for _ in range(150))))
    (tmp_path / "reads.fq").write_text("".join(f"@{n}\n{s}\n+\n{'I' * len(s)}\n"
                                               for n, s in reads))
    opts = MapOptions(reference_fasta=str(tmp_path / "g.fa"), engine=engine,
                      use_system_path=True, threads=2)
    rvmap, _files, _report = mr.run_mapping([tmp_path / "reads.fq"], opts)
    data = json.loads(rvmap.read_text())
    assert data["engine"] == engine
    assert data["summary"]["n_placed"] >= 790
    assert data["summary"]["n_reversed"] >= 390
    assert "nohit" in [u[0] for u in data["unused"]]
    contig = data["contigs"][0]
    assert contig["consensus"].count("N") < len(genome) * 0.05


def test_reference_lengths_line_by_line_and_fai(tmp_path):
    ref = tmp_path / "g.fa"
    ref.write_text(">chr1 desc\nACGT\nACG\n>chr2\nAA\n")
    assert mr.reference_lengths(ref) == {"chr1": 7, "chr2": 2}
    (tmp_path / "g.fa.fai").write_text("chr1\t700\t6\t4\t5\nchr2\t200\t20\t2\t3\n")
    assert mr.reference_lengths(ref) == {"chr1": 700, "chr2": 200}


def test_scan_inputs_cap_marks_incomplete(tmp_path):
    q = tmp_path / "reads.fa"
    q.write_text("".join(f">r{i}\nACGTACGT\n" for i in range(30)))
    info = mr.scan_inputs([q], max_records=10)
    assert info.n == 10 and not info.complete
    assert mr.choose_engine({"sat": 800}, info) != MAP_ENGINE_MAFFT


def test_long_reference_without_base_detail(tmp_path, monkeypatch):
    if shutil.which("minimap2") is None:
        pytest.skip("minimap2 not on PATH")
    from repeatviewer.widgets.contig_view import ContigModel

    monkeypatch.setattr(mr, "DETAIL_LIMIT", 1000)
    rng = random.Random(9)
    genome = "".join(rng.choice("ACGT") for _ in range(30_000))
    (tmp_path / "g.fa").write_text(">chr1\n" + genome + "\n")
    reads = []
    for i in range(20):
        p = rng.randrange(0, len(genome) - 3000)
        s = list(genome[p:p + 3000])
        s[100] = "A" if s[100] != "A" else "C"
        reads.append((f"L{i}", "".join(s)))
    (tmp_path / "long.fa").write_text("".join(f">{n}\n{s}\n" for n, s in reads))
    opts = MapOptions(reference_fasta=str(tmp_path / "g.fa"), engine=MAP_ENGINE_MINIMAP2,
                      use_system_path=True, threads=2)
    rvmap, _f, _r = mr.run_mapping([tmp_path / "long.fa"], opts)
    data = mr.load_rvmap(rvmap)
    contig = data["contigs"][0]
    assert contig["detailed"] is False and contig["reference"] == ""
    assert contig["bin_size"] >= 1 and contig["coverage"]
    model = ContigModel.from_json(contig)
    assert all(r.mism_pos for r in model.reads)          # marks come from the stored list
    assert any("coverage and the mapped reads" in w for w in data["warnings"])
