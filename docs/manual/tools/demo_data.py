"""Build the demo folder used for the manual's screenshots.

Everything is synthetic and made from a fixed random seed, so every run gives
the same sequences. Real tools (MAFFT, IQ-TREE, minimap2, Primer3) are then run
through RepeatViewer's own backends to make realistic results.
"""

from __future__ import annotations

import random
from pathlib import Path

RNG = random.Random(20261110)
BASES = "ACGT"


def _rand(n: int) -> str:
    return "".join(RNG.choice(BASES) for _ in range(n))


def _mutate(seq: str, rate: float) -> str:
    out = []
    for b in seq:
        r = RNG.random()
        if r < rate * 0.8:
            out.append(RNG.choice([x for x in BASES if x != b]))
        elif r < rate * 0.9:
            continue
        elif r < rate:
            out.append(b + RNG.choice(BASES))
        else:
            out.append(b)
    return "".join(out)


def _revcomp(seq: str) -> str:
    return seq[::-1].translate(str.maketrans("ACGTacgt", "TGCAtgca"))


def _write_fasta(path: Path, records: list[tuple[str, str]], width: int = 70) -> Path:
    with path.open("w", encoding="utf-8") as fh:
        for name, seq in records:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")
    return path


def make_ltr_element(folder: Path) -> Path:
    """A 5.6 kb LTR retrotransposon with annotations in its sidecar."""
    from repeatviewer.backends.sequence_features import SequenceFeature, save_features

    tsd = "GATCA"
    ltr = "TGTTG" + _rand(330) + "CAACA"
    pbs = "TGGTATCAGAGCCGGATTCGAGG"
    internal = _rand(4500)
    ppt = "AAGGAGAAGGGGGA"
    body = ltr + _rand(8) + pbs + internal + ppt + _rand(6) + _mutate(ltr, 0.004)
    seq = _rand(300) + tsd + body + tsd + _rand(300)
    path = _write_fasta(folder / "Ty3_gypsy_demo.fa", [("Ty3_gypsy_demo", seq)])

    off = 300 + len(tsd)
    l5 = (off + 1, off + len(ltr))
    pbs_s = l5[1] + 9
    int_s = pbs_s + len(pbs)
    ppt_s = int_s + len(internal)
    l3_s = ppt_s + len(ppt) + 6
    l3 = (l3_s, off + len(body))
    feats = [
        SequenceFeature(301, 305, "TSD5 GATCA", "TSD"),
        SequenceFeature(*l5, "LTR5", "LTR"),
        SequenceFeature(pbs_s, pbs_s + len(pbs) - 1, "PBS (tRNA-Met)", "PBS"),
        SequenceFeature(int_s + 120, int_s + 1110, "GAG", "domain"),
        SequenceFeature(int_s + 1300, int_s + 1600, "PROT", "domain"),
        SequenceFeature(int_s + 1700, int_s + 2400, "RT", "domain"),
        SequenceFeature(int_s + 2450, int_s + 2850, "RH", "domain"),
        SequenceFeature(int_s + 2900, int_s + 3800, "INT", "domain"),
        SequenceFeature(int_s + 3900, int_s + 4300, "CHD", "domain"),
        SequenceFeature(int_s + 60, int_s + 4400, "ORF1 (gag-pol)", "ORF", track=1),
        SequenceFeature(ppt_s, ppt_s + len(ppt) - 1, "PPT", "PPT"),
        SequenceFeature(*l3, "LTR3", "LTR"),
        SequenceFeature(l3[1] + 1, l3[1] + 5, "TSD3 GATCA", "TSD"),
    ]
    save_features(path, sequence_id="Ty3_gypsy_demo", features=feats)
    return path


def make_satellites(folder: Path) -> tuple[Path, Path, Path]:
    """Monomers of three subfamilies, an array reference, and ONT-like reads."""
    from repeatviewer.backends.sequence_features import SequenceFeature, save_features

    base = _rand(178)
    subfams = {"SatA": _mutate(base, 0.06), "SatB": _mutate(base, 0.12),
               "SatC": _mutate(base, 0.18)}
    monomers = []
    for fam, cons in subfams.items():
        for i in range(1, 9):
            chrom = RNG.choice(["chr1", "chr2", "chr3", "chr5"])
            start = RNG.randrange(1_000_000, 40_000_000)
            name = f"{fam}_{i:02d} {chrom}:{start}-{start + len(cons) - 1}"
            monomers.append((name, _mutate(cons, 0.03)))
    mono = _write_fasta(folder / "satellite_monomers.fa", monomers)

    unit = subfams["SatA"]
    array = _rand(400) + "".join(_mutate(unit, 0.01) for _ in range(12)) + _rand(400)
    ref = _write_fasta(folder / "satellite_array.fa", [("SatA_array", array)])
    feats = []
    pos = 401
    step = len(array[400:-400]) // 12
    for i in range(12):
        feats.append(SequenceFeature(pos, pos + step - 1, f"SatA monomer {i + 1}",
                                     "tandem repeat", track=i % 2))
        pos += step
    feats.append(SequenceFeature(401 + step * 3 + 40, 401 + step * 3 + 64,
                                 "CENP-B box", "others"))
    save_features(ref, sequence_id="SatA_array", features=feats)

    reads = []
    for i in range(1, 41):
        a = RNG.randrange(0, len(array) - 700)
        b = a + RNG.randrange(500, 1400)
        piece = _mutate(array[a:b], 0.035)
        if i % 3 == 0:
            piece = _revcomp(piece)
        reads.append((f"read_{i:03d}", piece))
    reads.append(("read_unrelated", _rand(600)))
    reads_path = _write_fasta(folder / "ont_reads.fa", reads)
    return mono, ref, reads_path


def paths(folder: Path) -> dict:
    """The demo files of an existing demo folder."""
    return {"ltr": folder / "Ty3_gypsy_demo.fa", "monomers": folder / "satellite_monomers.fa",
            "array": folder / "satellite_array.fa", "reads": folder / "ont_reads.fa",
            "tree": next(folder.glob("satellite_monomers.fa_*.treefile"), None),
            "paf": next(folder.glob("ont_reads.fa_*.paf"), None),
            "sam": next(folder.glob("ont_reads.fa_*.sam"), None),
            "pcr": next(folder.glob("satellite_array.fa_PCR_*.fa"), None),
            "gel": next(folder.glob("*.rvgel"), None),
            "aln": next(folder.glob("satellite_monomers.fa_*_details/*.aln.fa"), None),
            "tree_report": next(folder.glob("satellite_monomers.fa_*_details/*_report.html"), None),
            "map_report": next(folder.glob("ont_reads.fa_*_details/*_report.html"), None)}


def build(folder: Path, *, run_tools: bool = True, say=print) -> dict:
    folder.mkdir(parents=True, exist_ok=True)
    out = {"ltr": make_ltr_element(folder)}
    out["monomers"], out["array"], out["reads"] = make_satellites(folder)
    sub = folder / "candidates"
    sub.mkdir(exist_ok=True)
    for i, fam in enumerate(["DNA_TIR_01", "LINE_02", "Copia_03", "Helitron_04"]):
        _write_fasta(sub / f"{fam}.fa", [(fam, _rand(1200 + 400 * i))])
    if not run_tools:
        return out
    from repeatviewer.job_options import MapOptions, TreeOptions
    from repeatviewer.backends.tree_builder import build_tree_run
    from repeatviewer.backends.map_reference import run_mapping

    say("Building the demo tree (MAFFT + IQ-TREE)…")
    topts = TreeOptions(runs=1, ufboot=1000, alrt=1000, threads=2)
    out["tree_details"], out["tree_files"], out["tree_report"] = build_tree_run(
        [out["monomers"]], topts)
    say("Mapping the demo reads (minimap2)…")
    mopts = MapOptions(reference_fasta=str(out["array"]), engine="minimap2", threads=2)
    out["map_record"], out["map_files"], out["map_report"] = run_mapping(
        [out["reads"]], mopts)
    say("Mapping the demo reads again, as SAM for the browser view…")
    sopts = MapOptions(reference_fasta=str(out["array"]), engine="minimap2", threads=2,
                       minimap2_output="sam")
    run_mapping([out["reads"]], sopts)
    say("Placing demo primers and running in-silico PCR…")
    make_primers_and_pcr(out)
    out.update({k: v for k, v in paths(folder).items() if k not in out or out[k] is None})
    return out


def _read_seq(path: Path) -> str:
    return "".join(l.strip() for l in path.read_text().splitlines() if not l.startswith(">"))


def make_primers_and_pcr(out: dict) -> None:
    from repeatviewer.backends.primers import gel_doc as gd
    from repeatviewer.backends.primers import pcr
    from repeatviewer.backends.primers.oligos import Primer
    from repeatviewer.backends.primers.on_sequence import save_sequence_primers

    ltr = _read_seq(out["ltr"])
    rt_f = Primer("Ty3_RT_F", ltr[2900:2921], source="demo")
    rt_r = Primer("Ty3_RT_R", _revcomp(ltr[3380:3401]), source="demo")
    ltr_f = Primer("Ty3_LTR_F", ltr[400:420], source="demo")
    save_sequence_primers(out["ltr"], [rt_f, rt_r, ltr_f])

    arr = _read_seq(out["array"])
    sat_f = Primer("SatA_F", arr[430:450], source="demo")
    sat_r = Primer("SatA_R", _revcomp(arr[530:550]), source="demo")
    save_sequence_primers(out["array"], [sat_f, sat_r])
    settings = pcr.PcrSettings(mode=pcr.MODE_LADDER, max_product=2500)
    result = pcr.run_pcr([out["array"]], [sat_f, sat_r], settings)
    written = pcr.write_outputs(result)
    doc = gd.gel_from_pcr(result)
    doc.title = Path(written["fasta"]).stem
    doc.save(out["array"].parent / f"{doc.title}{gd.GEL_SUFFIX}")


if __name__ == "__main__":
    import sys
    root = Path(__file__).resolve().parents[3]
    sys.path.insert(0, str(root))
    target = Path(sys.argv[1]) if len(sys.argv) > 1 else Path("RepeatViewer demo")
    print(build(target))
