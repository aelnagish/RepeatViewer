"""RT-domain trees for retrotransposons (Tools → Tree, "RT domain from DNA").

Plant TE papers classify LTR retrotransposons by placing their reverse
transcriptase (RT) domain in a tree with reference RTs of known lineages
(Ty1/copia: Ale, Angela, SIRE…; Ty3/gypsy: Tekay, Athila, CRM…). This module:

1. reads **REXdb** (Neumann et al. 2019), which the user downloads once
   (``Viridiplantae_v4.0.fasta`` with ``Viridiplantae_v4.0.classification``
   next to it; not shipped with RepeatViewer);
2. builds one profile HMM per RT type (Ty1-RT, Ty3-RT, LINE-RT, …) from REXdb
   sequences (MAFFT + pyhmmer), cached so this happens once;
3. finds the RT domain in the six-frame translation of every TE copy and cuts
   it out (pyhmmer ``hmmsearch``);
4. names each copy's closest REXdb lineage: one profile per lineage, as in
   TEsorter and DANTE, searched with the copies' RT domains;
5. picks up to *n* REXdb reference RTs per lineage of the RT types found.

The tree builder then aligns the RT proteins and builds a protein tree.
"""

from __future__ import annotations

import hashlib
import random
import subprocess
import tempfile
from collections import Counter
from dataclasses import dataclass, field
from pathlib import Path

RT_TYPES = ("Ty1-RT", "Ty3-RT", "LINE-RT", "PARA-RT", "DIRS-RT", "Penelope-RT")
MIN_PROFILE_SEQS = 10          # an RT type needs this many REXdb sequences for a profile
PROFILE_SAMPLE = 200           # sequences per RT-type profile (seeded sample)
MIN_LINEAGE_SEQS = 5           # a lineage needs this many for its own profile
LINEAGE_SAMPLE = 60            # sequences per lineage profile
DEFAULT_EVALUE = 1e-5
REF_MIN_LENGTH = 0.9           # references: at least this share of their RT type's median
PROFILE_VERSION = 2


class RexdbError(ValueError):
    """REXdb is missing or not in the expected format."""


# --------------------------------------------------------------------------
# REXdb
# --------------------------------------------------------------------------


@dataclass
class Rexdb:
    fasta: Path
    rt: dict[str, tuple[str, str]]          # REXdb id → (RT type, protein sequence)
    lineage_path: dict[str, tuple[str, ...]]  # REXdb id → classification path

    def lineage(self, rexdb_id: str) -> str:
        path = self.lineage_path.get(rexdb_id, ())
        return path[-1] if path else "unclassified"

    def superfamily(self, rexdb_id: str) -> str:
        path = self.lineage_path.get(rexdb_id, ())
        for level in path:
            if "/" in level:            # Ty1/copia, Ty3/gypsy, …
                return level
        return path[1] if len(path) > 1 else (path[0] if path else "")

    @property
    def name(self) -> str:
        return self.fasta.stem


def classification_path(fasta: Path) -> Path:
    fasta = Path(fasta)
    return fasta.with_suffix(".classification")


def looks_like_rexdb(path: Path) -> bool:
    """A protein FASTA with its REXdb ``.classification`` table next to it."""
    path = Path(path)
    return path.is_file() and classification_path(path).is_file()


def load_rexdb(fasta: Path | str) -> Rexdb:
    fasta = Path(fasta).expanduser()
    if not fasta.is_file():
        raise RexdbError(f"REXdb file not found: {fasta}")
    table = classification_path(fasta)
    if not table.is_file():
        raise RexdbError(
            f"The REXdb classification table is missing: {table.name} must be next to "
            f"{fasta.name} (both are in the REXdb download)."
        )
    lineage_path: dict[str, tuple[str, ...]] = {}
    for line in table.read_text(encoding="utf-8", errors="replace").splitlines():
        parts = [p for p in line.rstrip("\n").split("\t") if p]
        if len(parts) >= 2:
            lineage_path[parts[0]] = tuple(parts[1:])
    rt: dict[str, tuple[str, str]] = {}
    name, chunks = None, []

    def flush() -> None:
        if name is None or "__" not in name:
            return
        domain, rid = name.split("__", 1)
        if domain in RT_TYPES:
            rt[rid] = (domain, "".join(chunks).replace("*", "").upper())

    with fasta.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                flush()
                name, chunks = line[1:].split()[0] if line[1:].strip() else "", []
            else:
                chunks.append(line.strip())
    flush()
    if not rt:
        raise RexdbError(f"{fasta.name} has no RT domains (expected headers like Ty1-RT__REXdb_ID1).")
    return Rexdb(fasta, rt, lineage_path)


# --------------------------------------------------------------------------
# Profiles (built from REXdb, cached)
# --------------------------------------------------------------------------


def _cache_key(fasta: Path) -> str:
    st = fasta.stat()
    raw = f"{fasta.resolve()}|{st.st_size}|{int(st.st_mtime)}|v{PROFILE_VERSION}"
    return hashlib.sha1(raw.encode()).hexdigest()[:16]


def _align_protein(seqs: list[tuple[str, str]], mafft: tuple[list[str], dict], work: Path,
                   threads: int) -> Path:
    raw = work / "in.fa"
    raw.write_text("".join(f">{n}\n{s}\n" for n, s in seqs), encoding="utf-8")
    command, env = mafft
    out = work / "aln.fa"
    with out.open("w", encoding="utf-8") as fh:
        subprocess.run([*command, "--auto", "--amino", "--quiet", "--thread", str(max(1, threads)),
                        str(raw)], stdout=fh, stderr=subprocess.PIPE, env=env, check=True)
    return out


def _safe(text: str) -> str:
    return "".join(c if c.isalnum() or c in "-_." else "_" for c in text)


def _build_one(name: str, seqs: list[tuple[str, str]], cached: Path, *,
               mafft: tuple[list[str], dict], threads: int):
    from pyhmmer import easel, plan7

    alphabet = easel.Alphabet.amino()
    if cached.is_file():
        with plan7.HMMFile(str(cached)) as fh:
            return fh.read()
    with tempfile.TemporaryDirectory() as tmp:
        aln = _align_protein(seqs, mafft, Path(tmp), threads)
        with easel.MSAFile(str(aln), digital=True, alphabet=alphabet, format="afa") as fh:
            msa = fh.read()
    msa.name = name.encode()
    hmm, _profile, _opt = plan7.Builder(alphabet).build_msa(msa, plan7.Background(alphabet))
    with cached.open("wb") as fh:
        hmm.write(fh)
    return hmm


def _sample(ids: list[str], n: int) -> list[str]:
    ids = sorted(ids)
    return ids if len(ids) <= n else sorted(random.Random(1).sample(ids, n))


def build_profiles(rexdb: Rexdb, *, mafft: tuple[list[str], dict], cache_dir: Path,
                   threads: int = 1, say=None) -> dict[str, object]:
    """One HMM per RT type with enough REXdb sequences (cached next time)."""
    folder = Path(cache_dir) / f"{rexdb.name}_{_cache_key(rexdb.fasta)}"
    folder.mkdir(parents=True, exist_ok=True)
    by_type: dict[str, list[str]] = {}
    for rid, (rt_type, _s) in rexdb.rt.items():
        by_type.setdefault(rt_type, []).append(rid)
    hmms: dict[str, object] = {}
    for rt_type in RT_TYPES:
        ids = by_type.get(rt_type, [])
        if len(ids) < MIN_PROFILE_SEQS:
            continue
        cached = folder / f"{rt_type}.hmm"
        if say and not cached.is_file():
            say(f"Building the {rt_type} profile from REXdb (only the first time)…")
        hmms[rt_type] = _build_one(rt_type, [(i, rexdb.rt[i][1]) for i in _sample(ids, PROFILE_SAMPLE)],
                                   cached, mafft=mafft, threads=threads)
    return hmms


def build_lineage_profiles(rexdb: Rexdb, rt_types: set[str], *, mafft: tuple[list[str], dict],
                           cache_dir: Path, threads: int = 1, say=None) -> dict[tuple[str, str], object]:
    """One HMM per (RT type, lineage) with enough REXdb sequences (cached)."""
    folder = Path(cache_dir) / f"{rexdb.name}_{_cache_key(rexdb.fasta)}" / "lineages"
    folder.mkdir(parents=True, exist_ok=True)
    groups: dict[tuple[str, str], list[str]] = {}
    for rid, (rt_type, _s) in rexdb.rt.items():
        if rt_type in rt_types:
            groups.setdefault((rt_type, rexdb.lineage(rid)), []).append(rid)
    todo = [k for k, ids in groups.items() if len(ids) >= MIN_LINEAGE_SEQS]
    missing = [k for k in todo if not (folder / f"{_safe(k[0])}__{_safe(k[1])}.hmm").is_file()]
    if say and missing:
        say(f"Building {len(missing)} lineage profiles from REXdb (only the first time)…")
    out = {}
    for key in sorted(todo):
        cached = folder / f"{_safe(key[0])}__{_safe(key[1])}.hmm"
        seqs = [(i, rexdb.rt[i][1]) for i in _sample(groups[key], LINEAGE_SAMPLE)]
        out[key] = _build_one(f"{key[0]}|{key[1]}", seqs, cached, mafft=mafft, threads=threads)
    return out


# --------------------------------------------------------------------------
# Finding RT domains in DNA
# --------------------------------------------------------------------------


def six_frames(dna: str) -> list[tuple[str, int, str]]:
    """(strand, frame offset, protein) for the six reading frames; stops → X."""
    from Bio.Seq import Seq

    seq = dna.upper().replace("-", "").replace(".", "")
    rc = str(Seq(seq).reverse_complement())
    out = []
    for strand, s in (("+", seq), ("-", rc)):
        for f in range(3):
            part = s[f:]
            part = part[: len(part) - len(part) % 3]
            protein = str(Seq(part).translate()).replace("*", "X") if part else ""
            out.append((strand, f, protein))
    return out


def dna_span(length: int, strand: str, frame: int, aa_from: int, aa_to: int) -> tuple[int, int]:
    """1-based DNA coordinates (on the input strand) of amino acids aa_from..aa_to."""
    start = frame + 3 * (aa_from - 1) + 1
    end = frame + 3 * aa_to
    if strand == "+":
        return start, end
    return length - end + 1, length - start + 1


@dataclass
class RTHit:
    name: str
    rt_type: str
    strand: str
    dna_start: int
    dna_end: int
    protein: str
    score: float
    evalue: float
    lineage: str = ""
    superfamily: str = ""
    lineage_score: float | None = None


def _name(value) -> str:
    return value.decode() if isinstance(value, (bytes, bytearray)) else str(value)


def find_rt_domains(records: list[tuple[str, str]], hmms: dict[str, object], *,
                    evalue: float = DEFAULT_EVALUE, cpus: int = 1) -> dict[str, RTHit]:
    """Best RT domain per record (name → hit) over all RT profiles."""
    import pyhmmer
    from pyhmmer import easel

    alphabet = easel.Alphabet.amino()
    frames: dict[str, tuple[str, str, int, str, int]] = {}
    seqs = []
    for i, (name, dna) in enumerate(records):
        for strand, f, protein in six_frames(dna):
            if len(protein) < 30:
                continue
            key = f"{i}|{strand}{f}"
            frames[key] = (name, strand, f, protein, len(dna.replace("-", "").replace(".", "")))
            seqs.append(easel.TextSequence(name=key.encode(), sequence=protein).digitize(alphabet))
    if not seqs:
        return {}
    block = easel.DigitalSequenceBlock(alphabet, seqs)
    best: dict[str, RTHit] = {}
    types = list(hmms)
    for rt_type, hits in zip(types, pyhmmer.hmmsearch([hmms[t] for t in types], block,
                                                       cpus=cpus, E=evalue, domE=evalue)):
        for hit in hits:
            if not hit.included:
                continue
            dom = hit.best_domain
            if dom.i_evalue > evalue:
                continue
            name, strand, f, protein, length = frames[_name(hit.name)]
            if name in best and best[name].score >= dom.score:
                continue
            a, b = dom.env_from, dom.env_to
            start, end = dna_span(length, strand, f, a, b)
            best[name] = RTHit(name, rt_type, strand, start, end, protein[a - 1:b],
                               float(dom.score), float(dom.i_evalue))
    return best


def classify(hits: dict[str, RTHit], rexdb: Rexdb, lineage_hmms: dict[tuple[str, str], object],
             *, cpus: int = 1) -> None:
    """Closest REXdb lineage for every RT domain (best-scoring lineage profile), in place."""
    import pyhmmer
    from pyhmmer import easel

    if not hits or not lineage_hmms:
        return
    alphabet = easel.Alphabet.amino()
    names = list(hits)
    block = easel.DigitalSequenceBlock(alphabet, [
        easel.TextSequence(name=str(i).encode(), sequence=hits[n].protein.replace("X", ""))
        .digitize(alphabet) for i, n in enumerate(names)])
    keys = list(lineage_hmms)
    superfamily_of = {}
    for rid in rexdb.rt:
        superfamily_of.setdefault(rexdb.lineage(rid), rexdb.superfamily(rid))
    # best lineage per (copy, RT type); the RT type the domain was found with wins
    best: dict[tuple[int, str], tuple[float, str]] = {}
    for key, top in zip(keys, pyhmmer.hmmsearch([lineage_hmms[k] for k in keys], block,
                                                 cpus=cpus, E=1.0)):
        rt_type, lineage = key
        for hit in top:
            i = int(_name(hit.name))
            if hit.score > best.get((i, rt_type), (float("-inf"), ""))[0]:
                best[(i, rt_type)] = (float(hit.score), lineage)
    for i, n in enumerate(names):
        h = hits[n]
        choice = best.get((i, h.rt_type))
        if choice is None:
            others = [v for (j, _t), v in best.items() if j == i]
            if not others:
                continue
            choice = max(others)
        score, lineage = choice
        h.lineage = lineage
        h.superfamily = superfamily_of.get(lineage, "")
        h.lineage_score = score


def pick_references(rexdb: Rexdb, rt_types: set[str], per_lineage: int, *,
                    seed: int = 1) -> list[tuple[str, str, str]]:
    """Up to ``per_lineage`` REXdb RTs per lineage of the RT types found:
    (REXdb id, lineage, protein)."""
    if per_lineage <= 0 or not rt_types:
        return []
    by_lineage: dict[tuple[str, str], list[str]] = {}
    for rid, (rt_type, _seq) in rexdb.rt.items():
        if rt_type in rt_types:
            by_lineage.setdefault((rt_type, rexdb.lineage(rid)), []).append(rid)
    # partial REXdb RTs make poor references: keep those at least REF_MIN_LENGTH
    # of the median length of their RT type
    lengths: dict[str, list[int]] = {}
    for rid, (rt_type, seq) in rexdb.rt.items():
        if rt_type in rt_types:
            lengths.setdefault(rt_type, []).append(len(seq))
    floor = {t: REF_MIN_LENGTH * sorted(v)[len(v) // 2] for t, v in lengths.items()}
    rng = random.Random(seed)
    out = []
    for key in sorted(by_lineage):
        ids = sorted(by_lineage[key])
        whole = [i for i in ids if len(rexdb.rt[i][1]) >= floor[key[0]]]
        ids = whole or ids
        chosen = ids if len(ids) <= per_lineage else sorted(rng.sample(ids, per_lineage))
        out += [(rid, key[1], rexdb.rt[rid][1]) for rid in chosen]
    return out


@dataclass
class RTResult:
    hits: dict[str, RTHit]
    missing: list[str]
    references: list[tuple[str, str, str]]
    by_type: Counter = field(default_factory=Counter)
    rexdb: Rexdb | None = None

    def table_tsv(self) -> str:
        lines = ["name\trt_type\tstrand\tdna_start\tdna_end\tlength_aa\tscore\tevalue"
                 "\tclosest_lineage\tsuperfamily\tlineage_score"]
        for h in self.hits.values():
            ls = "" if h.lineage_score is None else f"{h.lineage_score:.1f}"
            lines.append(f"{h.name}\t{h.rt_type}\t{h.strand}\t{h.dna_start}\t{h.dna_end}"
                         f"\t{len(h.protein)}\t{h.score:.1f}\t{h.evalue:.2g}\t{h.lineage}"
                         f"\t{h.superfamily}\t{ls}")
        for name in self.missing:
            lines.append(f"{name}\t\t\t\t\t\t\t\tno RT domain found\t\t")
        return "\n".join(lines) + "\n"


def reference_name(rexdb_id: str, lineage: str) -> str:
    return f"REXdb_{lineage}_{rexdb_id.replace('REXdb_', '')}"


def rt_domains_for(records: list[tuple[str, str]], rexdb_fasta: Path | str, *,
                   mafft: tuple[list[str], dict], cache_dir: Path, per_lineage: int = 3,
                   threads: int = 1, evalue: float = DEFAULT_EVALUE, say=None) -> RTResult:
    """Everything above in one call."""
    rexdb = load_rexdb(rexdb_fasta)
    hmms = build_profiles(rexdb, mafft=mafft, cache_dir=cache_dir, threads=threads, say=say)
    if say:
        say("Finding RT domains in the six reading frames…")
    hits = find_rt_domains(records, hmms, evalue=evalue, cpus=threads)
    types = {h.rt_type for h in hits.values()}
    lineage_hmms = build_lineage_profiles(rexdb, types, mafft=mafft, cache_dir=cache_dir,
                                          threads=threads, say=say)
    if say:
        say("Finding the closest REXdb lineage…")
    classify(hits, rexdb, lineage_hmms, cpus=threads)
    missing = [n for n, _s in records if n not in hits]
    refs = pick_references(rexdb, types, per_lineage)
    return RTResult(hits, missing, refs, Counter(h.rt_type for h in hits.values()), rexdb)
