"""Genome browser tracks: parsers, interval index, coverage, BAM readers."""

from __future__ import annotations

import itertools

import numpy as np
import pytest

from repeatviewer.backends import genome_tracks as gt


def _item(s, e, **kw):
    return gt.TrackItem(s, e, **kw)


def test_chrom_index_query_and_rows():
    ix = gt.ChromIndex([_item(0, 100), _item(50, 60), _item(200, 300), _item(10, 500)])
    assert sorted(ix.items[i].start for i in ix.query(150, 210)) == [10, 200]
    assert ix.query(600, 700) == []
    assert sorted(ix.items[i].start for i in ix.query(55, 56)) == [0, 10, 50]
    rows, n = ix.lanes()
    assert n == 3
    # overlapping items never share a row
    for a in range(len(ix.items)):
        for b in range(a + 1, len(ix.items)):
            if rows[a] == rows[b]:
                assert ix.items[a].end + 2 <= ix.items[b].start or \
                    ix.items[b].end + 2 <= ix.items[a].start


def test_coverage_bins_exact():
    starts = np.array([0, 10, 60], dtype=np.int64)
    ends = np.array([100, 20, 70], dtype=np.int64)
    cov = gt.coverage_bins(starts, ends, 10, length=100)
    assert cov.tolist() == [1, 2, 1, 1, 1, 1, 2, 1, 1, 1]
    cov64 = gt.coverage_bins(np.array([5]), np.array([200]), 64)
    assert cov64[0] == pytest.approx(59 / 64)
    assert cov64[1] == pytest.approx(1.0)


def test_items_coverage_respects_blocks():
    it = _item(0, 100, blocks=[(0, 10, 0), (90, 100, 10)])
    cov = gt.items_coverage([it], 0, 100, 100)
    assert cov[5] == 1 and cov[50] == 0 and cov[95] == 1


def test_match_chrom():
    names = ["Chr1", "chr2.1", "scaffold_9"]
    assert gt.match_chrom(names, "Chr1") == "Chr1"
    assert gt.match_chrom(names, "1") == "Chr1"
    assert gt.match_chrom(names, "chr2") == "chr2.1"
    assert gt.match_chrom(names, "9") == "scaffold_9"
    assert gt.match_chrom(names, "chr3") is None


def test_parse_locus():
    assert gt.parse_locus("chr5:35,618,000-35,630,000") == ("chr5", 35618000, 35630000)
    assert gt.parse_locus("chr5:100") == ("chr5", 100, 100)
    assert gt.parse_locus("hello") is None


def test_gff_gene_models(tmp_path):
    gff = tmp_path / "a.gff3"
    gff.write_text(
        "##gff-version 3\n"
        "c1\t.\tgene\t101\t1000\t.\t+\t.\tID=g1;Name=GeneA\n"
        "c1\t.\tmRNA\t101\t1000\t.\t+\t.\tID=t1;Parent=g1;Name=GeneA.1\n"
        "c1\t.\texon\t101\t300\t.\t+\t.\tParent=t1\n"
        "c1\t.\texon\t701\t1000\t.\t+\t.\tParent=t1\n"
        "c1\t.\tCDS\t201\t300\t.\t+\t0\tParent=t1\n"
        "c1\t.\tCDS\t701\t800\t.\t+\t2\tParent=t1\n"
        "c1\t.\tLTR_retrotransposon\t2001\t7000\t.\t-\t.\tID=te1;Name=Gypsy-1\n")
    t = gt.open_track(gff)
    assert t.kind == gt.KIND_ANNOTATION
    got = {it.name: it for it, _row in t.fetch("c1", 0, 10_000)}
    assert set(got) == {"GeneA", "GeneA.1", "Gypsy-1"}
    tx = got["GeneA.1"]
    assert tx.blocks == [(100, 300, 0), (700, 1000, 0)]
    assert tx.thick == [(200, 300), (700, 800)]
    assert got["Gypsy-1"].strand == "-"
    assert t.find("gypsy")[0][1].name == "Gypsy-1"


def test_gtf_without_transcript_lines(tmp_path):
    gtf = tmp_path / "a.gtf"
    gtf.write_text(
        'c1\tx\texon\t1\t50\t.\t-\t.\tgene_id "g"; transcript_id "t1"; gene_name "Foo";\n'
        'c1\tx\texon\t81\t100\t.\t-\t.\tgene_id "g"; transcript_id "t1"; gene_name "Foo";\n')
    t = gt.open_track(gtf)
    (it, _row), = t.fetch("c1", 0, 200)
    assert (it.start, it.end, it.name) == (0, 100, "Foo")
    assert len(it.blocks) == 2


def test_bed_track(tmp_path):
    bed = tmp_path / "a.bed"
    bed.write_text("c1\t10\t20\tTE1\t0\t+\nc1\t30\t40\tTE2\t0\t-\n")
    t = gt.open_track(bed)
    assert [(i.start, i.end) for i, _r in t.fetch("c1", 0, 100)] == [(10, 20), (30, 40)]


def test_paf_track_and_cs(tmp_path):
    paf = tmp_path / "a.paf"
    paf.write_text(
        "r1\t30\t0\t30\t+\tc1\t1000\t100\t130\t27\t30\t60\ttp:A:P\tcs:Z::10*ag:5-tt:5+c:8\n"
        "r2\t50\t0\t50\t-\tc1\t1000\t500\t550\t50\t50\t0\ttp:A:S\n")
    t = gt.open_track(paf)
    items = [i for i, _r in t.fetch("c1", 0, 1000)]
    r1 = items[0]
    assert r1.score == pytest.approx(0.9) and r1.mapq == 60 and not r1.secondary
    assert items[1].secondary and items[1].strand == "-"
    gt.cs_details(r1)
    assert r1.mism == [(110, "G")]
    assert r1.ins == [123]
    assert r1.blocks == [(100, 116, 0), (118, 131, 0)]
    cov = t.coverage("c1", 0, 1000, 1000)
    assert cov[105] == 1 and cov[300] == 0 and cov[520] == 1


def test_blast_track_sides(tmp_path):
    f = tmp_path / "hits.outfmt6"
    f.write_text("TE1\tchr1\t95.5\t100\t3\t1\t1\t100\t5100\t5001\t1e-40\t180\n"
                 "TE2\tchr2\t80.0\t50\t5\t0\t11\t60\t201\t250\t1e-10\t60\n")
    assert gt.track_kind_for(f) == "blast"
    t = gt.open_track(f, reference_names=["chr1", "chr2"])
    assert t.genome_side == "subject"
    (it, _r), = t.fetch("chr1", 0, 10_000)
    assert (it.start, it.end, it.strand, it.name) == (5000, 5100, "-", "TE1")
    assert it.score == pytest.approx(0.955)
    q = gt.BlastTrack(f, genome_side="query")
    assert q.chrom_names() == ["TE1", "TE2"]


def _sam_text():
    return ("@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:c1\tLN:200\n"
            "r1\t0\tc1\t11\t60\t5M2D5M1I4M\t*\t0\t0\tACGTAACGTATTTTT\t*\tNM:i:4\n"
            "r2\t16\tc1\t41\t30\t3S10M\t*\t0\t0\tNNNACGTACGTAC\t*\n"
            "r3\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\t*\n")


def test_sam_track(tmp_path):
    sam = tmp_path / "a.sam"
    sam.write_text(_sam_text())
    t = gt.open_track(sam)
    items = [i for i, _r in t.fetch("c1", 0, 200)]
    assert [i.name for i in items] == ["r1", "r2"]
    r1, r2 = items
    assert (r1.start, r1.end) == (10, 26)
    assert r1.blocks == [(10, 15, 0), (17, 22, 5), (22, 26, 11)]
    assert r1.ins == [22]
    assert r1.score == pytest.approx(1 - 4 / 17)
    assert r2.strand == "-" and r2.blocks == [(40, 50, 3)]
    ref = "A" * 200
    mism = gt.mismatches_against(r2, ref, 0, 200)
    assert mism[0] == (40, "C") or mism[0][0] == 41   # ACGT… vs AAAA…


def _write_bam(tmp_path, index=True):
    pysam = pytest.importorskip("pysam")
    sam = tmp_path / "a.sam"
    sam.write_text(_sam_text())
    bam = tmp_path / "a.bam"
    pysam.sort("-o", str(bam), str(sam))
    if index:
        pysam.index(str(bam))
    return bam


def test_memory_bam_matches_sam(tmp_path):
    bam = _write_bam(tmp_path, index=False)
    t = gt.MemoryBamTrack(bam)
    items = [i for i, _r in t.fetch("c1", 0, 200)]
    assert [i.name for i in items] == ["r1", "r2"]
    assert items[0].blocks == [(10, 15, 0), (17, 22, 5), (22, 26, 11)]
    assert items[0].seq == "ACGTAACGTATTTTT"
    assert items[0].score == pytest.approx(1 - 4 / 17)
    with pytest.raises(gt.TrackError):
        gt.MemoryBamTrack(bam, limit=10)


def test_indexed_bam_window_and_cache(tmp_path):
    bam = _write_bam(tmp_path)
    t = gt.open_track(bam)
    assert isinstance(t, gt.IndexedBamTrack)
    assert [i.name for i, _r in t.fetch("c1", 0, 200)] == ["r1", "r2"]
    assert [i.name for i, _r in t.fetch("c1", 45, 46)] == ["r2"]
    # wide window: no reads, coverage only once the cache is built
    assert t.fetch("c1", 0, gt.BAM_READ_WINDOW + 1) == []
    assert t.coverage("c1", 0, gt.BAM_READ_WINDOW + 1, 10) is None
    out = gt.build_bam_coverage(bam, bin_size=10)
    assert out.is_file()
    t.reload_coverage()
    cov = t.coverage("c1", 0, gt.BAM_READ_WINDOW + 1, 10)
    assert cov is not None and cov.max() > 0
    assert gt.load_coverage_cache(bam)["c1"][1] == pytest.approx(0.8)


def test_build_coverage_without_pysam(tmp_path, monkeypatch):
    bam = _write_bam(tmp_path, index=False)
    monkeypatch.setattr(gt, "pysam_available", lambda: False)
    out = gt.build_bam_coverage(bam, bin_size=10)
    data = gt.load_coverage_cache(bam)
    assert out.is_file() and data["__bin__"] == 10
    assert data["c1"][1] == pytest.approx(0.8)   # r1 blocks 10-15 and 17-20 in bin 10-19


def test_session_roundtrip(tmp_path):
    ref = tmp_path / "g.fa"
    ref.write_text(">c1\nACGT\n")
    bed = tmp_path / "a.bed"
    bed.write_text("c1\t0\t2\tx\n")
    path = gt.session_path_for(ref)
    assert path.name == "g.fa_browser.rvgb"
    gt.write_session(path, ref, [gt.open_track(bed)], "c1:1-4")
    data = gt.read_session(path)
    assert data["reference"] == str(ref) and data["tracks"][0]["path"] == str(bed)
    assert gt.session_path_for(ref).name == "g.fa_browser_2.rvgb"


def test_fasta_index_and_fetch(tmp_path):
    fa = tmp_path / "g.fa"
    fa.write_text(">c1 desc\nACGTACGTAC\nGTACG\n>c2\nTTTT\nGG\n")
    idx = gt.fasta_index(fa)
    assert {k: v.length for k, v in idx.items()} == {"c1": 15, "c2": 6}
    assert (tmp_path / "g.fa.fai").is_file()
    assert gt.fetch_sequence(fa, idx["c1"]) == "ACGTACGTACGTACG"
    assert gt.fetch_sequence(fa, idx["c1"], 8, 12) == "ACGT"
    assert gt.fetch_sequence(fa, gt.fasta_index(fa)["c2"]) == "TTTTGG"


def test_pack_right_and_left_keep_existing_rows():
    last, first = [], []
    a = gt._pack_right(np.array([0, 5, 20]), np.array([10, 15, 30]), last, first)
    assert a == [0, 1, 0]
    b = gt._pack_right(np.array([31, 33]), np.array([40, 50]), last, first)
    assert b == [1, 0]           # row 1 frees first (ended at 15)
    c = gt._pack_left(np.array([-20, -10]), np.array([-2, -3]), last, first)
    assert sorted(c) == [0, 1] and len(last) == 2
    # nothing overlaps within a row
    spans = {}
    for st, en, rr in zip([0, 5, 20, 31, 33, -20, -10], [10, 15, 30, 40, 50, -2, -3], a + b + c):
        spans.setdefault(rr, []).append((st, en))
    for row in spans.values():
        row.sort()
        assert all(x[1] <= y[0] for x, y in itertools.pairwise(row))


def test_bam_rows_stay_put_while_panning(tmp_path):
    pysam = pytest.importorskip("pysam")
    import random

    rng = random.Random(3)
    lines = ["@SQ\tSN:c1\tLN:20000"]
    for i, s in enumerate(sorted(rng.randrange(0, 19800) for _ in range(800))):
        lines.append(f"r{i}\t0\tc1\t{s + 1}\t60\t100M\t*\t0\t0\t*\t*")
    sam = tmp_path / "p.sam"
    sam.write_text("\n".join(lines) + "\n")
    bam = tmp_path / "p.bam"
    pysam.sort("-o", str(bam), str(sam))
    pysam.index(str(bam))
    t = gt.IndexedBamTrack(bam)
    t.TILE = 1000
    first = {it.name: row for it, row in t.fetch("c1", 5000, 6000)}
    serial = t.layout_version("c1", 5000, 6000)
    for s in range(5000, 9000, 250):          # pan right, then back left
        t.fetch("c1", s, s + 1000)
    for s in range(9000, 3000, -250):
        t.fetch("c1", s, s + 1000)
    again = {it.name: row for it, row in t.fetch("c1", 5000, 6000)}
    assert again == first and t.layout_version("c1", 5000, 6000) == serial
    # every read in the window, none twice, and no overlaps in a row
    names = [it.name for it, _r in t.fetch("c1", 3000, 4000)]
    assert len(names) == len(set(names))
    expect = {f"r{i}" for i, ln in enumerate(lines[1:])
              if int(ln.split("\t")[3]) - 1 < 4000 and int(ln.split("\t")[3]) + 99 > 3000}
    assert set(names) == expect
    by_row: dict = {}
    for it, r in t.fetch("c1", 2000, 10000):
        by_row.setdefault(r, []).append((it.start, it.end))
    for spans in by_row.values():
        spans.sort()
        assert all(a[1] <= b[0] for a, b in itertools.pairwise(spans))


def test_chrom_index_orders_items_without_a_python_sort():
    items = [_item(30, 40), _item(10, 50), _item(10, 20), _item(5, 9)]
    ix = gt.ChromIndex(items)
    assert [(it.start, it.end) for it in ix.items] == [(5, 9), (10, 20), (10, 50), (30, 40)]
    assert ix.starts.tolist() == [5, 10, 10, 30] and ix.ends.tolist() == [9, 20, 50, 40]
    in_order = [_item(1, 5), _item(2, 3), _item(2, 4)]
    same = list(in_order)
    gt.ChromIndex(in_order)
    assert all(a is b for a, b in zip(in_order, same))       # left as it was


def _many_reads_bam(tmp_path, n=3000, length=60000):
    pysam = pytest.importorskip("pysam")
    import random

    rng = random.Random(5)
    lines = [f"@SQ\tSN:c1\tLN:{length}"]
    for i, s in enumerate(sorted(rng.randrange(0, length - 200) for _ in range(n))):
        lines.append(f"r{i}\t0\tc1\t{s + 1}\t60\t100M\t*\t0\t0\t*\t*")
    sam = tmp_path / "m.sam"
    sam.write_text("\n".join(lines) + "\n")
    bam = tmp_path / "m.bam"
    pysam.sort("-o", str(bam), str(sam))
    pysam.index(str(bam))
    return bam


def test_bam_never_blocks_while_drawing_and_says_when_a_tile_is_ready(tmp_path, monkeypatch):
    import threading
    import time

    t = gt.IndexedBamTrack(_many_reads_bam(tmp_path))
    real = t._decode_tile
    gate = threading.Event()

    def slow(bam, key, tile):
        gate.wait(5)                     # a slow disk
        return real(bam, key, tile)

    monkeypatch.setattr(t, "_decode_tile", slow)
    ready = threading.Event()
    t.on_tile_ready = ready.set
    t0 = time.monotonic()
    with t.drawing():
        assert t.fetch("c1", 30000, 32000) == []           # nothing yet, no waiting
        assert not t.loaded_near("c1", 30000, 32000)
    assert time.monotonic() - t0 < 0.5
    gate.set()
    assert ready.wait(5)                                    # the view is told
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline:
        with t.drawing():
            got = t.fetch("c1", 30000, 32000)
        if got and t.tile_complete("c1", 30000, 32000):
            break
        time.sleep(0.02)
    assert got and t.loaded_near("c1", 30000, 32000)
    # outside drawing, a fetch waits for its reads as before
    assert t.fetch("c1", 1000, 3000)


def test_coverage_peak_is_remembered_per_region(tmp_path):
    t = gt.IndexedBamTrack(_many_reads_bam(tmp_path))
    t.fetch("c1", 10000, 12000)
    r = t._region
    calls = []
    real = r.depth
    r.depth = lambda a, b: calls.append((a, b)) or real(a, b)
    p1 = t.coverage_peak("c1", 9000, 13000)
    p2 = t.coverage_peak("c1", 9000, 13000)
    assert p1 == p2 > 0 and len(calls) == 1


def test_worker_prepares_tile_blocks_and_depth(tmp_path):
    """A tile read in the background arrives sorted, with its depth worked out, so
    adding it to the view while panning is cheap; the filter is applied there too."""
    import time

    t = gt.IndexedBamTrack(_many_reads_bam(tmp_path))
    t._prefetch("c1", 1)
    deadline = time.monotonic() + 5
    while time.monotonic() < deadline and t._cache.get(("c1", 1)) is None:
        time.sleep(0.01)
    items, _cut = t._cache.get(("c1", 1))
    assert items.prepared is not None
    _blk, (s0, depth) = items.prepared
    ix = gt.ChromIndex(list(items))
    assert s0 == t.TILE and (depth == ix.window_depth(s0, s0 + depth.size)).all()
    t.fetch("c1", t.TILE + 100, t.TILE + 2000)
    assert t._region.chunks[1]._depth[1] is depth          # used as is
    # with a read filter, the kept reads are prepared instead
    keep = lambda it: it.strand == "+"                     # noqa: E731
    pre = gt._prepare_tile(list(items), t.TILE, t.TILE, keep)
    assert pre.kept_for is keep and all(it.strand == "+" for it in pre.kept)
    assert pre.kept.prepared is not None


def test_bedgraph_signal_is_the_highest_value_per_column(tmp_path):
    p = tmp_path / "s.bedgraph"
    p.write_text("track type=bedGraph\nchr1\t0\t100\t1\nchr1\t100\t200\t5\n"
                 "chr1\t1000\t5000\t-2\n")
    t = gt.open_track(p)
    assert t.kind == gt.KIND_SIGNAL and t.format == "bedgraph"
    got = t.signal("chr1", 0, 6000, 6)
    assert got[0] == 5 and (got[1:5] == -2).all() and np.isnan(got[5])
    assert t.value_at("chr1", 150) == 5 and t.value_at("chr1", 600) is None
    assert t.fetch("chr1", 0, 100) == [] and t.row_count("chr1") == 0


def test_bigwig_signal_uses_the_file_zoom_levels(tmp_path):
    pyBigWig = pytest.importorskip("pyBigWig")
    p = tmp_path / "s.bw"
    bw = pyBigWig.open(str(p), "w")
    bw.addHeader([("chr1", 10000)])
    bw.addEntries(["chr1", "chr1"], [100, 5000], ends=[200, 5100], values=[3.0, 8.0])
    bw.close()
    t = gt.open_track(p)
    assert t.kind == gt.KIND_SIGNAL and t.format == "bigwig"
    got = t.signal("chr1", 0, 10000, 10)
    assert got[0] == 3 and got[5] == 8 and np.isnan(got[2])
    assert t.value_at("chr1", 150) == 3
    # zoomed past one base per column, and a view past the end of the sequence
    assert (t.signal("chr1", 150, 154, 8) == 3).all()
    tail = t.signal("chr1", 9000, 11000, 4)
    assert np.isnan(tail[2:]).all()


_VCF = """##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##INFO=<ID=AF,Number=A,Type=Float,Description="AF">
##FORMAT=<ID=GT,Number=1,Type=String,Description="GT">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\ts1\ts2\ts3
chr1\t120\t.\tAT\tA\t20\t.\tAF=0.25\tGT\t./.\t0|0\t0|1
chr1\t150\trs1\tA\tG\t50\tPASS\t.\tGT\t0/1\t1/1\t0/0
chr1\t500\trs3\tC\tT,G\t.\tq10\t.\tGT\t1/2\t0/0\t.
"""


def test_vcf_variants_and_genotypes(tmp_path):
    p = tmp_path / "v.vcf"
    p.write_text(_VCF)
    t = gt.open_track(p)
    assert t.kind == gt.KIND_VARIANT and t.samples == ["s1", "s2", "s3"]
    items, gts = t.variants("chr1", 0, 1000)
    assert [(it.start, it.end, it.name, it.kind) for it in items] == [
        (119, 121, "AT>A", "INDEL"), (149, 150, "rs1", "SNP"), (499, 500, "rs3", "SNP")]
    assert gts.tolist() == [[gt.GT_NOCALL, gt.GT_HOMREF, gt.GT_HET],
                            [gt.GT_HET, gt.GT_HOMALT, gt.GT_HOMREF],
                            [gt.GT_HET, gt.GT_HOMREF, gt.GT_NOCALL]]
    assert t.allele_freq(items[0], gts[0]) == 0.25            # from INFO AF
    assert t.allele_freq(items[1], gts[1]) == 0.5             # from the genotypes
    assert items[2].info["FILTER"] == "q10" and t.count("chr1") == 3
    assert t.search("rs3")[1] == 1
    assert gt.genotype_code("1") == gt.GT_HOMALT and gt.genotype_code(".") == gt.GT_NOCALL


def test_indexed_vcf_is_read_by_window(tmp_path):
    pysam = pytest.importorskip("pysam")
    p = tmp_path / "v.vcf"
    p.write_text(_VCF)
    pysam.tabix_index(str(p), preset="vcf", force=True)
    t = gt.open_track(tmp_path / "v.vcf.gz")
    assert t._vf is not None and t.chrom_names() == ["chr1"]
    items, gts = t.variants("chr1", 100, 200)
    assert [it.name for it in items] == ["AT>A", "rs1"]
    assert gts[1].tolist() == [gt.GT_HET, gt.GT_HOMALT, gt.GT_HOMREF]
    assert t.variants("chr1", 0, gt.VARIANT_WINDOW + 10)[0] == []    # too wide


def _spliced_bam(tmp_path):
    pysam = pytest.importorskip("pysam")
    hdr = {"HD": {"VN": "1.6", "SO": "coordinate"}, "SQ": [{"SN": "chr1", "LN": 10000}]}
    out = tmp_path / "rna.bam"
    rows = [(100, [(0, 50), (3, 400), (0, 50)], "+")] * 3 + \
           [(120, [(0, 30), (3, 400), (0, 30), (3, 1000), (0, 40)], "-")] * 2 + \
           [(2000, [(0, 60)], None)]
    with pysam.AlignmentFile(str(out), "wb", header=hdr) as f:
        for k, (pos, cig, xs) in enumerate(sorted(rows, key=lambda r: r[0])):
            a = pysam.AlignedSegment()
            a.query_name, a.reference_id, a.reference_start = f"r{k}", 0, pos
            a.cigartuples = cig
            a.query_sequence = "A" * sum(n for op, n in cig if op == 0)
            a.mapping_quality = 60
            if xs:
                a.set_tag("XS", xs)
            f.write(a)
    pysam.index(str(out))
    return out


def test_splice_junctions_are_counted_per_strand(tmp_path):
    bam = _spliced_bam(tmp_path)
    t = gt.open_track(bam)
    j = t.junctions("chr1", 0, 5000).tolist()
    assert j == [[150, 550, 3, 2, 0], [580, 1580, 0, 2, 0]]
    assert t.junctions("chr1", 1600, 5000).tolist() == []
    # the junction track of the same file, kept in a session as such
    jt = gt.open_track(bam, settings={"view": "junctions"})
    assert jt.kind == gt.KIND_JUNCTION and jt.to_json()["settings"]["view"] == "junctions"
    assert jt.junctions("chr1", 0, 5000).shape == t.junctions("chr1", 0, 5000).shape
    # a SAM read's blocks give the same junctions (gap >= MIN_INTRON)
    it = gt.TrackItem(0, 900, blocks=[(0, 50, 0), (450, 500, 50), (505, 900, 100)])
    assert gt.item_junctions(it) == [(50, 450, 2)]
