"""Chromosome painting (Chorus2): BED parsing, pool selection, outputs, and
running a stand-in Chorus2 from a fake conda env."""

from __future__ import annotations

import json
import os
import random
import stat
import sys
from pathlib import Path

import pytest

pytest.importorskip("primer3")

from repeatviewer.backends.primers import painting as pt  # noqa: E402
from repeatviewer.backends.primers.thermo import PrimerError  # noqa: E402


def _genome(tmp_path: Path) -> tuple[Path, dict[str, int]]:
    rng = random.Random(3)
    lengths = {"chr1": 300_000, "chr2": 200_000}
    path = tmp_path / "genome.fa"
    with open(path, "w") as fh:
        for name, n in lengths.items():
            seq = "".join(rng.choice("ACGT") for _ in range(n))
            fh.write(f">{name} test\n")
            for i in range(0, n, 80):
                fh.write(seq[i:i + 80] + "\n")
    return path, lengths


def _chorus_bed(path: Path, lengths: dict[str, int], *, every: int = 60, scores=False) -> Path:
    rng = random.Random(5)
    with open(path, "w") as fh:
        for chrom, n in lengths.items():
            for start in range(1, n - 45, every):
                seq = "".join(rng.choice("ACGT") for _ in range(45))
                row = [chrom, str(start), str(start + 44), seq]
                if scores:
                    row += [str(rng.randint(0, 2000)), "+"]
                fh.write("\t".join(row) + "\n")
    return path


def test_parse_bed_and_errors(tmp_path):
    _g, lengths = _genome(tmp_path)
    bed = _chorus_bed(tmp_path / "p.bed", lengths)
    probes = pt.parse_probe_bed(bed)
    assert probes[0].chrom == "chr1" and probes[0].start == 1 and probes[0].length == 45
    assert probes[0].score is None and probes[0].strand == "+"
    bad = tmp_path / "bad.bed"
    bad.write_text("chr1\t1\t45\n")
    with pytest.raises(PrimerError, match="Chorus2 BED"):
        pt.parse_probe_bed(bad)
    bad.write_text("chr1\t1\t45\tHELLO\n")
    with pytest.raises(PrimerError, match="not a DNA"):
        pt.parse_probe_bed(bad)


def test_regions_parse():
    assert pt.parse_regions("chr1, chr2:1-5,000,000; chr3:10-20") == [
        ("chr1", 1, 0), ("chr2", 1, 5_000_000), ("chr3", 10, 20)]
    assert pt.parse_regions("") == []
    with pytest.raises(PrimerError):
        pt.parse_regions("chr1:50-10")


def test_selection_spacing_density_pool_strand(tmp_path):
    _g, lengths = _genome(tmp_path)
    probes = pt.parse_probe_bed(_chorus_bed(tmp_path / "p.bed", lengths, every=30))
    # Spacing: probes 30 bp apart overlap (45 nt), so every other one goes at gap 25.
    s = pt.PaintingSettings(spacing=25, strand=pt.STRAND_PLUS)
    chosen, stats = pt.select_probes(probes, s, lengths)
    for a, b in zip(chosen, chosen[1:]):
        if a.chrom == b.chrom:
            assert b.start > a.end + 25
    assert stats["after_spacing"] == len(chosen) < len(probes)
    # Density 1 per kb -> about 500 probes on 500 kb.
    s = pt.PaintingSettings(target=pt.TARGET_DENSITY, density_per_kb=1.0, strand=pt.STRAND_PLUS)
    chosen, _ = pt.select_probes(probes, s, lengths)
    assert 450 <= len(chosen) <= 510
    # Pool size: exactly the number asked for (there are enough).
    s = pt.PaintingSettings(target=pt.TARGET_POOL, pool_size=300)
    chosen, _ = pt.select_probes(probes, s, lengths)
    assert len(chosen) == 300
    per = {c: sum(1 for p in chosen if p.chrom == c) for c in lengths}
    assert abs(per["chr1"] / per["chr2"] - 1.5) < 0.2  # even along the genome
    # Both strands: alternate, the minus ones reverse-complemented.
    by_pos = {(p.chrom, p.start): p for p in probes}
    chr1 = [p for p in chosen if p.chrom == "chr1"]
    assert [p.strand for p in chr1[:4]] == ["+", "-", "+", "-"]
    assert chr1[1].sequence == pt.revcomp(by_pos[("chr1", chr1[1].start)].sequence)
    # Regions.
    s = pt.PaintingSettings(regions="chr2:1-100000", strand=pt.STRAND_PLUS)
    chosen, _ = pt.select_probes(probes, s, lengths)
    assert chosen and all(p.chrom == "chr2" and p.end <= 100_000 for p in chosen)


def test_selection_ngs_score_window(tmp_path):
    _g, lengths = _genome(tmp_path)
    probes = pt.parse_probe_bed(_chorus_bed(tmp_path / "p.bed", lengths, scores=True))
    s = pt.PaintingSettings(use_ngs=True, ngs_low_quantile=0.1, ngs_high_quantile=0.9,
                            spacing=0, strand=pt.STRAND_PLUS)
    chosen, stats = pt.select_probes(probes, s, lengths)
    lo, hi = stats["score_window"]
    assert all(lo < p.score < hi for p in chosen)
    assert 0.7 * len(probes) < len(chosen) < 0.85 * len(probes)


def test_import_writes_report_bed_fasta_order(tmp_path):
    genome, lengths = _genome(tmp_path)
    bed = _chorus_bed(tmp_path / "genome.fa.bed", lengths)
    s = pt.PaintingSettings(source=pt.SOURCE_IMPORT, import_bed=str(bed), genome_fasta=str(genome),
                            target=pt.TARGET_POOL, pool_size=200, primer5="TAATACGACTCACTATAGGG",
                            pool_name="Chr1Paint")
    main, files, payload = pt.run_painting(s)
    assert main == tmp_path / "genome.fa.bed_painting.html"
    names = {f.name for f in files}
    assert {"genome.fa.bed_painting_probes.bed", "genome.fa.bed_painting_probes.fa",
            "genome.fa.bed_painting_order.tsv", "genome.fa.bed_painting.json"} <= names
    html = main.read_text()
    assert "Chromosome painting probes" in html and "<svg" in html and "Methods text" in html
    assert "—" not in html and "Zhang" in html
    out_bed = next(f for f in files if f.name.endswith("_probes.bed"))
    rows = [ln.split("\t") for ln in out_bed.read_text().splitlines()[1:]]
    assert len(rows) == 200 == payload["selected"]
    src = {(p.chrom, p.start): p for p in pt.parse_probe_bed(bed)}
    first = rows[0]
    assert (first[0], int(first[1]) + 1) in src  # 0-based BED start
    order = next(f for f in files if f.name.endswith("_order.tsv")).read_text().splitlines()
    assert order[0].startswith("pool\tname\tsequence")
    assert order[1].split("\t")[2].startswith("TAATACGACTCACTATAGGG")
    assert order[1].split("\t")[1].startswith("Chr1Paint_chr1_")
    data = json.loads(next(f for f in files if f.suffix == ".json").read_text())
    assert data["kind"] == "painting" and sum(r["probes"] for r in data["per_chromosome"]) == 200
    # A second run never overwrites: _2.
    main2, _f, _p = pt.run_painting(s)
    assert main2.name == "genome.fa.bed_painting_2.html"


def test_nothing_selected_explains(tmp_path):
    genome, lengths = _genome(tmp_path)
    bed = _chorus_bed(tmp_path / "x.bed", lengths)
    s = pt.PaintingSettings(source=pt.SOURCE_IMPORT, import_bed=str(bed), regions="chrX")
    with pytest.raises(PrimerError, match="No probe is left"):
        pt.run_painting(s)


def test_server_command(tmp_path):
    genome, _ = _genome(tmp_path)
    s = pt.PaintingSettings(genome_fasta=str(genome), threads=16, use_ngs=True,
                            ngs_fastq="/data/r1.fq.gz, /data/r2.fq.gz")
    text = pt.server_command_text(s)
    lines = text.splitlines()
    assert lines[0] == "conda activate chorus2"
    assert lines[1].startswith("Chorus2 -g genome.fa -i genome.fa -s chorus2_out -t 16 -l 45")
    assert "--docker True" in lines[1]
    assert lines[2].startswith("ChorusNGSfilter -g genome.fa -i r1.fq.gz,r2.fq.gz -z gz")
    assert "-p chorus2_out/genome.fa.bed -o chorus2_out/probes_kmer_score.bed" in lines[2]
    assert str(tmp_path) not in text


def test_missing_env_is_explained(tmp_path):
    inst = pt.find_chorus2(str(tmp_path / "no_such_env"))
    assert not inst.ok
    genome, _ = _genome(tmp_path)
    s = pt.PaintingSettings(genome_fasta=str(genome), conda_env=str(tmp_path / "no_such_env"))
    with pytest.raises(PrimerError, match="conda create -n chorus2"):
        pt.run_painting(s)


FAKE_CHORUS = r'''
import argparse, os, random
ap = argparse.ArgumentParser()
for a in ("-g", "-i", "-s", "-t", "-l", "--homology", "-d", "--step", "--ploidy", "--docker",
          "-j", "-b"):
    ap.add_argument(a)
a = ap.parse_args()
os.makedirs(a.s, exist_ok=True)
assert os.path.isfile(a.j) and os.path.isfile(a.b), "jellyfish / bwa paths"
L = int(a.l)
seqs, name = {}, None
for line in open(a.i):
    line = line.strip()
    if line.startswith(">"):
        name = line[1:].split()[0]; seqs[name] = []
    else:
        seqs[name].append(line)
with open(os.path.join(a.s, os.path.basename(a.i) + ".bed"), "w") as out:
    for n, parts in seqs.items():
        s = "".join(parts)
        for st in range(1, len(s) - L, 200):
            out.write(f"{n}\t{st}\t{st + L - 1}\t{s[st - 1:st - 1 + L]}\n")
print("Job finshed!!")
'''


@pytest.mark.skipif(sys.platform == "win32", reason="Chorus2 runs on Linux / macOS")
def test_run_here_with_a_stand_in_chorus2(tmp_path):
    env = tmp_path / "envs" / "chorus2"
    bin_dir = env / "bin"
    bin_dir.mkdir(parents=True)
    os.symlink(sys.executable, bin_dir / "python")
    (bin_dir / "Chorus2").write_text(FAKE_CHORUS)
    for tool in ("bwa", "jellyfish"):
        p = bin_dir / tool
        p.write_text("#!/bin/sh\nexit 0\n")
        p.chmod(p.stat().st_mode | stat.S_IEXEC)
    (bin_dir / "Chorus2").chmod(0o755)
    inst = pt.find_chorus2(str(env))
    assert inst.ok, inst.message
    genome, lengths = _genome(tmp_path)
    s = pt.PaintingSettings(genome_fasta=str(genome), conda_env=str(env), length=40,
                            strand=pt.STRAND_PLUS)
    lines = []
    main, files, payload = pt.run_painting(s, progress=lines.append)
    assert main.name == "genome.fa_painting.html"
    assert (main.parent / "genome.fa_painting_details" / "chorus2" / "genome.fa.bed").is_file()
    expected = sum(len(range(1, n - 40, 200)) for n in lengths.values())
    assert payload["input"] == expected and payload["selected"] == expected
    assert any("Chorus2" in ln for ln in lines)
    fa = next(f for f in files if f.name.endswith("_probes.fa")).read_text().splitlines()
    seq = "".join(ln.strip() for ln in genome.read_text().splitlines()[1:3001] if not
                  ln.startswith(">"))
    assert fa[1] == seq[:40]
