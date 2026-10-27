"""Synthetic, deterministic data for the primer tests (Tools > Primers).

Not a test module: imported by ``test_primers_*.py`` (pytest puts ``tests/``
on ``sys.path``).
"""

from __future__ import annotations

import random
from dataclasses import dataclass
from pathlib import Path

_COMP = str.maketrans("ACGTacgt", "TGCAtgca")


def revcomp(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


@dataclass
class PrimerData:
    root: Path
    gene: Path
    monomers: Path
    array: Path
    te_aln: Path
    locus: Path
    promoter: Path
    genome: Path
    refgene: Path
    primers: Path
    gene_seq: str
    mono: str
    locus_seq: str
    te_start: int  # 1-based TE start in locus.fa
    te_end: int
    gene_f: str
    gene_r: str
    sat_f: str
    sat_r: str


def _write(path: Path, records: list[tuple[str, str]]) -> Path:
    with open(path, "w", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), 70):
                fh.write(seq[i:i + 70] + "\n")
    return path


def make_data(out: Path, *, genome_flank: int = 60_000) -> PrimerData:
    out.mkdir(parents=True, exist_ok=True)
    rng = random.Random(7)

    def rnd(n: int, gc: float = 0.45) -> str:
        w = [(1 - gc) / 2, gc / 2, gc / 2, (1 - gc) / 2]
        return "".join(rng.choices("ACGT", weights=w, k=n))

    def mutate(s: str, rate: float) -> str:
        return "".join(rng.choice("ACGT".replace(c, "")) if rng.random() < rate else c for c in s)

    left, rep = rnd(800), rnd(400).lower()
    mid, right = rnd(700) + "AG" * 30 + rnd(300), rnd(800)
    gene_seq = left + rep + mid + right
    gene = _write(out / "gene.fa", [("gene1", gene_seq)])

    mono = rnd(180, 0.38)
    monomers = _write(out / "monomers.fa",
                      [(f"CroSatX_{i + 1}", mutate(mono, 0.03)) for i in range(40)])
    array = _write(out / "array.fa", [("array1", "".join(mutate(mono, 0.02) for _ in range(12)))])

    te = rnd(1200, 0.42)
    sub_a, sub_b = mutate(te, 0.04), mutate(te, 0.06)
    rows = []
    for i in range(12):
        s = list(mutate(sub_a if i < 7 else sub_b, 0.01))
        if i % 4 == 0:
            k = rng.randrange(100, 1000)
            s[k:k + 5] = "-----"
        rows.append((f"{'TEfamA' if i < 7 else 'TEfamB'}_{i + 1}", "".join(s)))
    te_aln = _write(out / "te_family_aln.fa", rows)

    tsd = rnd(5)
    flank_l, flank_r = rnd(900), rnd(900)
    locus_seq = flank_l + tsd + te + tsd + flank_r
    te_start = len(flank_l) + len(tsd) + 1
    te_end = te_start + len(te) - 1
    locus = _write(out / "locus.fa", [("locus1", locus_seq)])

    cpg = "".join(rng.choice(["CG", "A", "T", "C", "G", "AT", "TA", "CA", "TG"]) for _ in range(900))
    promoter = _write(out / "promoter.fa", [("prom1", rnd(300) + cpg + rnd(300))])

    c1 = rnd(genome_flank) + gene_seq.upper() + rnd(genome_flank)
    c2 = rnd(genome_flank) + "".join(mutate(mono, 0.02) for _ in range(60)) + rnd(genome_flank)
    c3 = rnd(genome_flank) + locus_seq + rnd(genome_flank)
    genome = _write(out / "genome.fa", [("chr1", c1), ("chr2", c2), ("chr3", c3)])
    refgene = _write(out / "refgene.fa", [("actin_like", rnd(1500, 0.5))])

    gene_f = left[100:121]
    gene_r = revcomp(mid[200:221])
    sat_f, sat_r = mono[20:40], revcomp(mono[120:140])
    primers = _write(out / "primers.fa", [("Gene_F", gene_f), ("Gene_R", gene_r),
                                          ("Sat_F", sat_f), ("Sat_R", sat_r)])
    return PrimerData(out, gene, monomers, array, te_aln, locus, promoter, genome, refgene,
                      primers, gene_seq, mono, locus_seq, te_start, te_end,
                      gene_f, gene_r, sat_f, sat_r)
