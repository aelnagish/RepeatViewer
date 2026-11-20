"""Large and compressed FASTA: index, regions, bgzip (repeatviewer.backends.large_fasta)."""

from __future__ import annotations

import gzip
import random

import pytest

from repeatviewer.backends import large_fasta as lf


def _records(seed: int = 7):
    rnd = random.Random(seed)
    return [
        ("chr1", "first one", "".join(rnd.choice("ACGTacgtN") for _ in range(1003))),
        ("chr2", "", "".join(rnd.choice("ACGT") for _ in range(240))),
        ("scaffold_5", "long", "".join(rnd.choice("ACGT") for _ in range(70_001))),
        ("empty", "", ""),
    ]


def _fasta_text(records, width=60, crlf=False) -> str:
    nl = "\r\n" if crlf else "\n"
    out = []
    for name, desc, seq in records:
        out.append(f">{name} {desc}".rstrip() + nl)
        for i in range(0, len(seq), width):
            out.append(seq[i: i + width] + nl)
    return "".join(out)


def _check_all(index, records):
    by = {n: s for n, _, s in records}
    assert index.names() == [n for n, _, _ in records]
    for name, seq in by.items():
        assert index.entry(name).length == len(seq)
        if not seq:
            continue
        for start, end in [(1, len(seq)), (1, 1), (len(seq), len(seq)), (59, 62),
                           (60, 121), (max(1, len(seq) // 2), len(seq) - 1)]:
            if start <= end:
                assert lf.fetch(index, name, start, end) == seq[start - 1: end], (name, start, end)


@pytest.mark.parametrize("width,crlf", [(60, False), (80, True), (10_000_000, False)])
def test_plain_index_and_fetch(tmp_path, width, crlf):
    recs = _records()
    path = tmp_path / "genome.fna"
    path.write_bytes(_fasta_text(recs, width, crlf).encode())
    index = lf.get_index(path)
    assert index.kind == lf.KIND_PLAIN
    _check_all(index, recs)
    assert (tmp_path / "genome.fna.fai").is_file()
    again = lf.load_index(path)
    assert again is not None and again.names() == index.names()
    _check_all(again, recs)


def test_gzip_and_bgzip(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    recs = _records()
    text = _fasta_text(recs).encode()
    gz = tmp_path / "genome.fa.gz"
    gz.write_bytes(gzip.compress(text))
    assert lf.is_fasta_name(gz) and lf.storage_kind(gz) == lf.KIND_GZIP
    index = lf.get_index(gz)
    _check_all(index, recs)
    assert not (tmp_path / "genome.fa.gz.fai").exists()  # plain gzip: index in the cache

    bgz = lf.convert_to_bgzip(gz)
    assert bgz.name == "genome_bgzip.fa.gz"
    assert lf.storage_kind(bgz) == lf.KIND_BGZF
    assert gzip.decompress(bgz.read_bytes()) == text
    bindex = lf.get_index(bgz)
    assert len(bindex.blocks) >= 2
    _check_all(bindex, recs)
    assert (tmp_path / "genome_bgzip.fa.gz.fai").is_file()
    assert (tmp_path / "genome_bgzip.fa.gz.gzi").is_file()
    _check_all(lf.load_index(bgz), recs)


def test_irregular_lines_still_read(tmp_path, monkeypatch):
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    path = tmp_path / "odd.fa"
    path.write_text(">a\nACGT\nAC\nGTACGT\n>b\nTTTT\n")
    index = lf.get_index(path)
    assert index.entry("a").line_bases == 0
    assert lf.fetch(index, "a", 3, 9) == "GTACGTA"
    assert lf.fetch(index, "b") == "TTTT"
    assert not (tmp_path / "odd.fa.fai").exists()


def test_errors(tmp_path):
    dup = tmp_path / "dup.fa"
    dup.write_text(">x\nAC\n>x\nGT\n")
    with pytest.raises(lf.LargeFastaError, match="more than one"):
        lf.build_index(dup)
    bad = tmp_path / "bad.fa"
    bad.write_text("ACGT\n")
    with pytest.raises(lf.LargeFastaError, match="header"):
        lf.build_index(bad)


def test_regions_and_output(tmp_path):
    recs = _records()
    path = tmp_path / "Csat.fa"
    path.write_text(_fasta_text(recs))
    index = lf.get_index(path)
    r = lf.parse_region("scaffold_5:1,000-2,000", index)
    assert (r.name, r.start, r.end) == ("scaffold_5", 1000, 2000)
    assert lf.parse_region("chr2", index) == lf.Region("chr2", 1, 240)
    assert lf.parse_region("100-90000", index, current="scaffold_5").end == 70_001
    with pytest.raises(lf.LargeFastaError, match="past the end"):
        lf.parse_region("chr2:500-600", index)
    with pytest.raises(lf.LargeFastaError, match="name too"):
        lf.parse_region("10-20", index)

    assert lf.opening_region(index, "scaffold_5") == lf.Region("scaffold_5", 1, 70_001)
    monkeypatch_limit = lf.WHOLE_SEQUENCE_LIMIT
    try:
        lf.WHOLE_SEQUENCE_LIMIT = 50_000
        assert lf.opening_region(index, "scaffold_5", window=5_000) == lf.Region("scaffold_5", 1, 5_000)
    finally:
        lf.WHOLE_SEQUENCE_LIMIT = monkeypatch_limit

    out = lf.region_output_path(index, r)
    assert out.name == "Csat.fa_scaffold_5_1000-2000.fa"
    lf.write_region(index, r, out)
    lines = out.read_text().splitlines()
    assert lines[0] == ">scaffold_5:1000-2000 from Csat.fa"
    assert "".join(lines[1:]) == recs[2][2][999:2000]
    whole = lf.Region("chr1", 1, 1003)
    assert lf.region_output_path(index, whole).name == "Csat.fa_chr1.fa"
    out.write_text("x")
    assert lf.numbered_path(out).name == "Csat.fa_scaffold_5_1000-2000_2.fa"
    assert lf.format_bp(438_200_000) == "438.2 Mb" and lf.format_bp(3_550_000_000) == "3.55 Gb"


def test_small_chunks_and_odd_endings(tmp_path, monkeypatch):
    """Chunk edges anywhere, no final newline, blank lines, irregular records."""
    monkeypatch.setenv("XDG_CACHE_HOME", str(tmp_path / "cache"))
    monkeypatch.setattr(lf, "_CHUNK", 7)
    rnd = random.Random(3)
    for trial in range(60):
        recs = []
        text = []
        for r in range(rnd.randint(1, 4)):
            seq = "".join(rnd.choice("ACGT") for _ in range(rnd.randint(0, 50)))
            width = rnd.randint(1, 12)
            lines = [seq[i: i + width] for i in range(0, len(seq), width)]
            if lines and rnd.random() < 0.3:  # make it irregular
                k = rnd.randrange(len(lines))
                lines[k:k + 1] = [lines[k][:1], lines[k][1:]] if len(lines[k]) > 1 else [lines[k]]
            body = "".join(line + "\n" for line in lines)
            if rnd.random() < 0.3:
                body += "\n"
            text.append(f">s{r}\n{body}")
            recs.append((f"s{r}", "", seq))
        data = "".join(text)
        if rnd.random() < 0.5:
            data = data.rstrip("\n")
        path = tmp_path / f"t{trial}.fa"
        path.write_text(data)
        index = lf.build_index(path)
        for name, _, seq in recs:
            assert index.entry(name).length == len(seq)
            if seq:
                a = rnd.randint(1, len(seq))
                b = rnd.randint(a, len(seq))
                assert lf.fetch(index, name, a, b) == seq[a - 1: b], (data, name, a, b)
