"""Read colouring and read details (IGV-style 'Colour by' and the read popup)."""

from __future__ import annotations

import itertools
import random

import pytest

from repeatviewer.backends import alignment_style as style
from repeatviewer.backends import genome_tracks as gt


def _paired_sam(n_pairs: int = 60) -> str:
    """Proper pairs with ~300 bp inserts, plus one far pair, one outward pair, one
    pair with the mate on another sequence and one with the mate unmapped."""
    rnd = random.Random(1)
    lines = ["@SQ\tSN:chr1\tLN:20000", "@SQ\tSN:chr2\tLN:5000",
             "@RG\tID:rgA\tSM:sample1\tLB:lib1", "@RG\tID:rgB\tSM:sample2\tLB:lib1"]
    seq = "A" * 50

    def pair(name, p1, p2, tlen, f1=99, f2=147, rg="rgA", rnext="="):
        lines.append(f"{name}\t{f1}\tchr1\t{p1}\t60\t50M\t{rnext}\t{p2}\t{tlen}\t{seq}\t*\t"
                     f"RG:Z:{rg}\tHP:i:{1 + (p1 % 2)}")
        if rnext == "=":
            lines.append(f"{name}\t{f2}\tchr1\t{p2}\t60\t50M\t=\t{p1}\t{-tlen}\t{seq}\t*\t"
                         f"RG:Z:{rg}\tHP:i:{1 + (p1 % 2)}")

    for i in range(n_pairs):
        p1 = 100 + i * 150
        size = 300 + rnd.randint(-20, 20)
        pair(f"p{i}", p1, p1 + size - 50, size, rg="rgA" if i % 2 else "rgB")
    pair("far", 1000, 6000, 5050, f1=97, f2=145)           # not proper, too far apart
    pair("out", 2800, 3000, 250, f1=81, f2=161)             # left read reverse, mate on +
    lines.append("other\t97\tchr1\t4000\t60\t50M\tchr2\t100\t0\t" + seq + "\t*\tRG:Z:rgA")
    lines.append("lone\t105\tchr1\t4500\t0\t50M\t=\t4500\t0\t" + seq + "\t*\tRG:Z:rgA")
    return "\n".join(lines) + "\n"


@pytest.fixture
def sam_track(tmp_path):
    path = tmp_path / "p.sam"
    path.write_text(_paired_sam())
    return gt.open_track(path)


def _by_name(track, name):
    return [it for it, _r in track.fetch("chr1", 0, 20000) if it.name == name]


def test_sam_read_groups_and_tags(sam_track):
    assert sam_track.read_groups == {"rgA": {"SM": "sample1", "LB": "lib1"},
                                     "rgB": {"SM": "sample2", "LB": "lib1"}}
    it = _by_name(sam_track, "p1")[0]
    assert style.tag(it, "RG") == "rgA"
    assert style.tag(it, "HP") in (1, 2)
    assert style.mate(it)[0] == "chr1"


def test_insert_size_and_orientation_colours(sam_track):
    lo, hi = sam_track.insert_range("chr1")
    assert 270 <= lo <= 300 <= hi <= 330
    kw = {"chrom": "chr1", "insert_range": (lo, hi)}
    normal = _by_name(sam_track, "p3")[0]
    assert style.pair_orientation(normal) == "LR"
    assert style.read_color(normal, style.COLOR_INSERT, **kw) == style.NEUTRAL
    far = _by_name(sam_track, "far")[0]
    assert style.read_color(far, style.COLOR_INSERT, **kw) == style.LARGER
    out = _by_name(sam_track, "out")[0]
    assert style.pair_orientation(out) == "RL"
    assert style.read_color(out, style.COLOR_INSERT, **kw) == style.OUTWARD
    other = _by_name(sam_track, "other")[0]
    assert style.read_color(other, style.COLOR_INSERT, **kw) == style.category_color("chr2")
    lone = _by_name(sam_track, "lone")[0]
    assert style.read_color(lone, style.COLOR_INSERT, **kw) == style.MATE_UNMAPPED


def test_strand_group_and_tag_colours(sam_track):
    first, second = _by_name(sam_track, "p2")
    # first-of-pair strand: both mates get the colour of the first read's strand
    assert style.read_color(first, style.COLOR_FIRST_STRAND) == \
        style.read_color(second, style.COLOR_FIRST_STRAND) == style.FWD
    assert style.read_color(second, style.COLOR_STRAND) == style.REV
    groups = sam_track.read_groups
    a, b = _by_name(sam_track, "p1")[0], _by_name(sam_track, "p2")[0]
    assert style.read_color(a, style.COLOR_READ_GROUP) != style.read_color(b, style.COLOR_READ_GROUP)
    assert style.read_color(a, style.COLOR_SAMPLE, read_groups=groups) == \
        style.category_color("sample1")
    assert style.read_color(a, style.COLOR_LIBRARY, read_groups=groups) == \
        style.read_color(b, style.COLOR_LIBRARY, read_groups=groups)
    assert style.read_color(a, style.COLOR_TAG, tag_name="HP") == \
        style.category_color(style.tag(a, "HP"))
    assert style.read_color(a, style.COLOR_TAG, tag_name="ZZ") == style.NEUTRAL
    assert style.read_color(a, style.COLOR_MAPQ).startswith("#")


def test_details_rows(sam_track):
    it = _by_name(sam_track, "p1")[0]
    rows = dict(style.details(it, chrom="chr1"))
    assert rows["Read"] == "p1"
    assert rows["MAPQ"] == "60"
    assert "proper pair" in rows["Flags"]
    assert rows["Mate"].startswith("chr1:")
    assert rows["Pair orientation"] == "LR"
    assert rows["RG"] == "rgA"


def test_feature_details(tmp_path):
    bed = tmp_path / "a.bed"
    bed.write_text("chr1\t10\t50\tTE1\t0\t-\n")
    t = gt.open_track(bed)
    it = t.fetch("chr1", 0, 100)[0][0]
    rows = dict(style.details(it, chrom="chr1"))
    assert rows["Name"] == "TE1" and rows["Position"] == "chr1:11–50 (-)"


def _bam(tmp_path):
    pysam = pytest.importorskip("pysam")
    sam = tmp_path / "p.sam"
    sam.write_text(_paired_sam())
    bam = tmp_path / "p.bam"
    pysam.sort("-o", str(bam), str(sam))
    return pysam, bam


def test_memory_bam_has_mates_tags_and_groups(tmp_path):
    _pysam, bam = _bam(tmp_path)
    t = gt.MemoryBamTrack(bam)
    assert set(t.read_groups) == {"rgA", "rgB"}
    it = _by_name(t, "p1")[0]
    assert style.tag(it, "RG") == "rgA" and style.tag(it, "HP") in (1, 2)
    assert style.mate(it)[0] == "chr1" and style.mate(it)[2] != 0
    lo, hi = t.insert_range("chr1")
    assert lo <= 300 <= hi


def test_indexed_bam_insert_range_and_groups(tmp_path):
    pysam, bam = _bam(tmp_path)
    pysam.index(str(bam))
    t = gt.open_track(bam)
    assert isinstance(t, gt.IndexedBamTrack)
    assert t.read_groups["rgB"]["SM"] == "sample2"
    lo, hi = t.insert_range("chr1")
    assert lo <= 300 <= hi
    far = _by_name(t, "far")[0]
    assert style.read_color(far, style.COLOR_INSERT, chrom="chr1",
                            insert_range=(lo, hi)) == style.LARGER
    rows = dict(style.details(far, chrom="chr1"))
    assert rows["Insert size"] == "5,050 bp" and rows["HP"] in ("1", "2")


def test_bam_tags_parser():
    buf = (b"RGZrgA\x00" + b"NMC\x03" + b"XSi" + (-5).to_bytes(4, "little", signed=True)
           + b"AAAx")
    assert gt.bam_tags(buf) == {"RG": "rgA", "NM": 3, "XS": -5, "AA": "x"}


def test_group_colours_follow_header_order(sam_track):
    pal = style.group_colors(sam_track.read_groups, style.COLOR_SAMPLE)
    assert list(pal) == ["sample1", "sample2"] and len(set(pal.values())) == 2
    a = _by_name(sam_track, "p1")[0]
    assert style.read_color(a, style.COLOR_SAMPLE, read_groups=sam_track.read_groups,
                            palette=pal) == pal["sample1"]
    assert style.group_colors(sam_track.read_groups, style.COLOR_LIBRARY) == \
        {"lib1": style.CATEGORY[0]}


def _clip_sam() -> str:
    ref = "ACGT" * 100
    read = ref[100:130] + "GGGGG" + ref[130:160] + "TTTTTTTTTTTTTTTTTTTTTTTTT"
    qual = "I" * 10 + "#" + "I" * (len(read) - 11)      # '#' = Q2 at offset 10
    read = read[:10] + ("A" if read[10] != "A" else "C") + read[11:]
    return ("@SQ\tSN:chr1\tLN:400\n"
            f"c1\t0\tchr1\t101\t60\t30M5I30M25S\t*\t0\t0\t{read}\t{qual}\n"), ref


def test_clips_insertions_and_qualities_from_sam(tmp_path):
    text, ref = _clip_sam()
    path = tmp_path / "c.sam"
    path.write_text(text)
    it = gt.open_track(path).fetch("chr1", 0, 400)[0][0]
    assert style.clips(it) == (0, 25)
    assert style.insertions(it) == [(130, 5)]
    q = style.qualities(it)
    assert q[10] == 2 and q[0] == 40
    mm = gt.mismatches_against(it, ref, 0, 400, offsets=True)
    assert mm == [(110, mm[0][1], 10)]
    assert style.quality_alpha(q[10]) < 0.2 and style.quality_alpha(30) == 1.0
    rows = dict(style.details(it))
    assert rows["Clipped"] == "0 bp left · 25 bp right"
    assert rows["Insertions"].startswith("1 (largest 5 bp)")


def test_clips_insertions_and_qualities_from_bam(tmp_path):
    pysam = pytest.importorskip("pysam")
    text, _ref = _clip_sam()
    sam = tmp_path / "c.sam"
    sam.write_text(text)
    bam = tmp_path / "c.bam"
    pysam.sort("-o", str(bam), str(sam))
    for t in (gt.MemoryBamTrack(bam), None):
        if t is None:
            pysam.index(str(bam))
            t = gt.open_track(bam)
        it = t.fetch("chr1", 0, 400)[0][0]
        assert style.clips(it) == (0, 25)
        assert style.insertions(it) == [(130, 5)]
        assert style.qualities(it)[10] == 2


def test_paf_unaligned_ends_count_as_clips(tmp_path):
    paf = tmp_path / "a.paf"
    paf.write_text("q1\t1000\t100\t900\t-\tchr1\t5000\t200\t1000\t780\t800\t60\n")
    it = gt.open_track(paf).fetch("chr1", 0, 5000)[0][0]
    # reverse strand: the query's tail is on the reference's left
    assert style.clips(it) == (100, 100)


def test_group_rows_per_group(sam_track):
    fn = style.group_function(style.COLOR_SAMPLE, read_groups=sam_track.read_groups)
    sam_track.set_grouping(("sample",), fn)
    counts = sam_track.group_rows("chr1")
    assert set(counts) == {"sample1", "sample2"}
    seen: dict = {}
    for it, g, row in sam_track.fetch_groups("chr1", 0, 20000):
        assert g == fn(it)
        seen.setdefault((g, row), []).append((it.start, it.end))
    for spans in seen.values():          # no two reads of a group overlap in one row
        spans.sort()
        assert all(b[0] >= a[1] for a, b in itertools.pairwise(spans))
    sam_track.set_grouping(None, None)
    assert not sam_track.grouped and sam_track.fetch_groups("chr1", 0, 100) == []


def test_indexed_bam_grouping_rebuilds_rows(tmp_path):
    pysam, bam = _bam(tmp_path)
    pysam.index(str(bam))
    t = gt.open_track(bam)
    fn = style.group_function("pair_orientation")
    t.set_grouping(("pair_orientation",), fn)
    got = t.fetch_groups("chr1", 0, 20000)
    labels = {g for _it, g, _r in got}
    assert {"LR", "RL"} <= labels
    assert set(t.group_rows("chr1")) == labels
    v1 = t.layout_version("chr1", 0, 20000)
    assert isinstance(v1, tuple)


def test_sort_keys():
    class It:
        def __init__(self, start, seq, strand="+", mapq=60):
            self.start, self.end, self.seq, self.strand, self.mapq = start, start + len(seq), seq, strand, mapq
            self.blocks = [(start, self.end, 0)]
            self.name = f"r{start}"
            self.info = {}
    ref = "AAAAAAAAAA"
    items = [It(0, "AAAAA"), It(1, "AAACA"), It(2, "AACAA"), It(3, "AGAAA")]
    # base at 4: r0 A, r1 C, r2 C, r3 G -> C (2 reads) first, then G, the reference A last
    fn = style.sort_key_function("base", 4, items, reference=ref)
    assert [it.name for it in sorted(items, key=fn)] == ["r1", "r2", "r3", "r0"]
    assert style.base_at(items[0], 9) == "-"
    fn = style.sort_key_function("mapq", 4, items)
    items[2].mapq = 70
    assert min(items, key=fn).name == "r2"


def _flagged_sam(tmp_path):
    seq = "A" * 50
    lines = ["@SQ\tSN:chr1\tLN:5000"]
    for name, fl, mq in (("ok", 0, 60), ("dup", 1024, 60), ("qc", 512, 60), ("sec", 256, 60),
                         ("supp", 2048, 60), ("low", 0, 3)):
        lines.append(f"{name}\t{fl}\tchr1\t101\t{mq}\t50M\t*\t0\t0\t{seq}\t*")
    path = tmp_path / "f.sam"
    path.write_text("\n".join(lines) + "\n")
    return path


def test_read_filters(tmp_path):
    t = gt.open_track(_flagged_sam(tmp_path))
    def names():
        return sorted(it.name for it, _r in t.fetch("chr1", 0, 5000))
    assert len(names()) == 6
    t.set_filter(*style.read_filter({}))                  # IGV defaults
    assert names() == ["low", "ok", "sec", "supp"]
    t.set_filter(*style.read_filter({"hide_secondary": True, "hide_supplementary": True,
                                     "min_mapq": 10, "hide_duplicates": False}))
    assert names() == ["dup", "ok"]
    assert t.row_count("chr1") == 2
    t.set_filter(*style.read_filter({"hide_duplicates": False, "hide_qcfail": False}))
    assert len(names()) == 6


def test_view_as_pairs_puts_mates_in_one_row(sam_track):
    sam_track.set_pairs(True)
    rows: dict = {}
    for it, row in sam_track.fetch("chr1", 0, 20000):
        rows.setdefault(it.name, set()).add(row)
    for name in ("p1", "p2", "p10", "far"):
        assert len(rows[name]) == 1, name
    sam_track.set_pairs(False)


def test_bam_view_as_pairs_and_filter(tmp_path):
    pysam, bam = _bam(tmp_path)
    pysam.index(str(bam))
    t = gt.open_track(bam)
    t.set_pairs(True)
    rows: dict = {}
    for it, row in t.fetch("chr1", 0, 20000):
        rows.setdefault(it.name, set()).add(row)
    assert all(len(r) == 1 for n, r in rows.items() if n.startswith("p"))
    t.set_filter(*style.read_filter({"min_mapq": 10}))
    assert "lone" not in {it.name for it, _r in t.fetch("chr1", 0, 20000)}   # MAPQ 0


def test_feature_sequences_from_the_reference(tmp_path):
    ref = "AAAACCCCGGGGTTTT" * 10
    gff = tmp_path / "m.gff3"
    gff.write_text(
        "chr1\t.\tmRNA\t1\t16\t.\t-\t.\tID=m1;Name=M1\n"
        "chr1\t.\texon\t1\t4\t.\t-\t.\tParent=m1\n"
        "chr1\t.\texon\t9\t16\t.\t-\t.\tParent=m1\n"
        "chr1\t.\tCDS\t3\t4\t.\t-\t0\tParent=m1\n"
        "chr1\t.\tCDS\t9\t12\t.\t-\t0\tParent=m1\n")
    it = next(i for i, _r in gt.open_track(gff).fetch("chr1", 0, 100) if i.name == "M1")
    def fetch(s, e):
        return ref[s:e]
    kinds = [k for k, _l, _n in style.feature_sequences(it, fetch)]
    assert kinds == ["genomic", "spliced", "cds"]
    assert style.feature_sequence(it, fetch) == style.reverse_complement(ref[:16])
    assert style.feature_sequence(it, fetch, "spliced") == \
        style.reverse_complement("AAAA" + "GGGGTTTT")
    assert style.feature_sequence(it, fetch, "cds") == style.reverse_complement("AA" + "GGGG")
