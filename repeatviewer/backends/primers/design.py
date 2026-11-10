"""Primer design: Primer3 plus the satellite modes.

Modes (``DesignSettings.mode``):

* ``pcr``        – standard PCR pair (Primer3).
* ``qpcr``       – qPCR pair (+ optional hydrolysis probe), and a matching pair
                   on a single-copy reference gene.
* ``ladder``     – satellite ladder PCR: outward-facing primers in one monomer,
                   giving 1×, 2×, 3× … monomer bands on a tandem array.
* ``fish``       – FISH oligo probes conserved across monomers / copies.

Everything here is Qt-free; the dialogs only fill :class:`DesignSettings`.
"""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import asdict, dataclass, field

import numpy as np

from repeatviewer.backends.primers import binding as bd
from repeatviewer.backends.primers.thermo import (
    PrimerError,
    ThermoSettings,
    analyze_oligo,
    calc_tm,
    gc_percent,
    heterodimer,
    require_primer3,
    reverse_complement,
)

MODE_PCR = "pcr"
MODE_QPCR = "qpcr"
MODE_LADDER = "ladder"
MODE_FISH = "fish"
DESIGN_MODES = (MODE_PCR, MODE_QPCR, MODE_LADDER, MODE_FISH)
MODE_LABELS = {
    MODE_PCR: "Standard PCR",
    MODE_QPCR: "qPCR",
    MODE_LADDER: "Ladder PCR",
    MODE_FISH: "Satellite repeat oligos",
    "painting": "Chromosome painting (Chorus2)",
}
# Chromosome painting is a FISH-probe mode run by Chorus2 (see ``painting.py``),
# not by the primer engine: it is in the Mode list but not in DESIGN_MODES.
MODE_PAINTING = "painting"
ALL_MODES = DESIGN_MODES + (MODE_PAINTING,)
# Mode dropdown: family > group > modes (headers are not selectable).
MODE_FAMILIES: tuple[tuple[str, tuple[tuple[str, tuple[str, ...]], ...]], ...] = (
    ("PCR primers", (
        ("General", (MODE_PCR, MODE_QPCR)),
        ("Satellites", (MODE_LADDER,)),
    )),
    ("FISH probes", (
        ("Satellites", (MODE_FISH,)),
        ("Chromosomes", (MODE_PAINTING,)),
    )),
)
FISH_MODES = (MODE_FISH, MODE_PAINTING)
MAX_TEMPLATE_BP = 2_000_000

FISH_LABELS = ("None", "FAM", "Cy3", "Cy5", "TAMRA", "Alexa Fluor 488", "Alexa Fluor 594",
               "ATTO 550", "Biotin", "Digoxigenin")


@dataclass
class DesignSettings:
    mode: str = MODE_PCR
    # Primer conditions
    size_min: int = 18
    size_opt: int = 20
    size_max: int = 25
    tm_min: float = 57.0
    tm_opt: float = 60.0
    tm_max: float = 63.0
    gc_min: float = 40.0
    gc_max: float = 60.0
    max_tm_diff: float = 2.0
    product_min: int = 100
    product_max: int = 1000
    num_return: int = 5
    min_three_prime_distance: int = 3
    max_poly_x: int = 4
    gc_clamp: int = 0
    max_end_stability: float = 9.0
    pick_probe: bool = False
    probe_size_min: int = 18
    probe_size_opt: int = 22
    probe_size_max: int = 30
    probe_tm_min: float = 67.0
    probe_tm_opt: float = 69.0
    probe_tm_max: float = 72.0
    thermo: ThermoSettings = field(default_factory=ThermoSettings)
    # Tails
    forward_tail: str = ""
    reverse_tail: str = ""
    # qPCR
    reference_fasta: str = ""
    # Ladder / FISH: how many input copies a primer or probe must bind
    min_coverage: float = 0.90
    coverage_mismatches: int = 1
    coverage_three_prime_exact: int = 3
    # FISH
    fish_label: str = "Cy3"
    fish_avoid_fasta: str = ""

    def __post_init__(self) -> None:
        if isinstance(self.thermo, dict):
            self.thermo = ThermoSettings.from_dict(self.thermo)
        if self.mode not in ALL_MODES:
            self.mode = MODE_PCR
        self.size_min = max(5, int(self.size_min))
        self.size_max = max(self.size_min, int(self.size_max))
        self.size_opt = min(max(int(self.size_opt), self.size_min), self.size_max)
        self.product_min = max(1, int(self.product_min))
        self.product_max = max(self.product_min, int(self.product_max))
        self.num_return = max(1, min(100, int(self.num_return)))

    def to_dict(self) -> dict:
        d = asdict(self)
        d["thermo"] = self.thermo.to_dict()
        return d

    @classmethod
    def from_dict(cls, data: dict | None) -> DesignSettings:
        data = dict(data or {})
        known = {k: data[k] for k in cls.__dataclass_fields__ if k in data}
        return cls(**known)


def mode_preset(mode: str) -> dict:
    """Condition values the dialog fills in when a mode is chosen."""
    if mode == MODE_QPCR:
        return dict(size_min=18, size_opt=20, size_max=24, tm_min=58.0, tm_opt=60.0,
                    tm_max=62.0, gc_min=40.0, gc_max=60.0, max_tm_diff=1.0,
                    product_min=70, product_max=180, gc_clamp=0, max_poly_x=3)
    if mode == MODE_LADDER:
        return dict(coverage_mismatches=1, size_min=18, size_opt=20, size_max=25, tm_min=56.0, tm_opt=60.0,
                    tm_max=64.0, gc_min=30.0, gc_max=70.0, max_tm_diff=2.0,
                    product_min=40, product_max=2000)
    if mode == MODE_FISH:
        return dict(size_min=30, size_opt=36, size_max=45, tm_min=60.0, tm_opt=72.0,
                    tm_max=85.0, gc_min=35.0, gc_max=65.0, coverage_mismatches=3,
                    min_coverage=0.9)
    return dict(size_min=18, size_opt=20, size_max=25, tm_min=57.0, tm_opt=60.0,
                tm_max=63.0, gc_min=40.0, gc_max=60.0, max_tm_diff=2.0,
                product_min=100, product_max=1000, gc_clamp=0, max_poly_x=4)


# --------------------------------------------------------------------------
# Results
# --------------------------------------------------------------------------


@dataclass
class DesignedOligo:
    name: str
    sequence: str                 # full oligo (lower-case tail + binding part)
    kind: str                     # forward / reverse / probe
    start: int = 0                # 1-based on the template (binding part)
    end: int = 0
    strand: str = "+"
    tm: float = 0.0
    gc: float = 0.0
    length: int = 0
    end_dg: float = 0.0
    hairpin_dg: float = 0.0
    hairpin_tm: float = 0.0
    self_dimer_dg: float = 0.0
    degeneracy: int = 1
    penalty: float = 0.0
    warnings: list[str] = field(default_factory=list)
    hairpin_lines: list[str] = field(default_factory=list)
    self_dimer_lines: list[str] = field(default_factory=list)
    extra: dict = field(default_factory=dict)

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class DesignedPair:
    rank: int
    name: str
    template: str
    forward: DesignedOligo | None = None
    reverse: DesignedOligo | None = None
    probe: DesignedOligo | None = None
    product_size: int = 0
    penalty: float = 0.0
    tm_diff: float = 0.0
    dimer_dg: float = 0.0
    dimer_lines: list[str] = field(default_factory=list)
    group: str = ""               # e.g. "M" / "U", "target" / "reference"
    notes: list[str] = field(default_factory=list)
    specificity: dict = field(default_factory=dict)
    coverage: dict = field(default_factory=dict)
    bands: list[dict] = field(default_factory=list)   # expected gel bands
    extra: dict = field(default_factory=dict)

    def oligos(self) -> list[DesignedOligo]:
        return [o for o in (self.forward, self.reverse, self.probe) if o is not None]

    def to_dict(self) -> dict:
        d = asdict(self)
        for key in ("forward", "reverse", "probe"):
            obj = getattr(self, key)
            d[key] = obj.to_dict() if obj is not None else None
        return d


def make_oligo(name: str, seq: str, kind: str, *, start: int = 0, end: int = 0,
               strand: str = "+", thermo: ThermoSettings | None = None,
               penalty: float = 0.0, tail: str = "") -> DesignedOligo:
    full = (tail.lower() + seq.upper()) if tail else seq
    st = analyze_oligo(full, thermo)
    return DesignedOligo(
        name=name, sequence=st.sequence, kind=kind, start=start, end=end, strand=strand,
        tm=st.tm, gc=st.gc, length=st.binding_length, end_dg=st.end_dg,
        hairpin_dg=st.hairpin.dg, hairpin_tm=st.hairpin.tm, self_dimer_dg=st.self_dimer.dg,
        degeneracy=st.degeneracy, penalty=round(penalty, 3), warnings=list(st.warnings),
        hairpin_lines=list(st.hairpin.lines), self_dimer_lines=list(st.self_dimer.lines),
    )


def finish_pair(pair: DesignedPair, thermo: ThermoSettings) -> DesignedPair:
    """Tm difference and cross-dimer of a pair."""
    if pair.forward is not None and pair.reverse is not None:
        pair.tm_diff = round(abs(pair.forward.tm - pair.reverse.tm), 2)
        dim = heterodimer(pair.forward.sequence, pair.reverse.sequence, thermo)
        pair.dimer_dg = dim.dg
        pair.dimer_lines = list(dim.lines)
        if dim.dg < -9.0:
            pair.notes.append(f"primer-dimer ΔG {dim.dg:.1f} kcal/mol")
    return pair


# --------------------------------------------------------------------------
# Primer3
# --------------------------------------------------------------------------


def primer3_global_args(s: DesignSettings, *, num_return: int | None = None,
                        product: tuple[int, int] | list | None = None,
                        pick_probe: bool | None = None) -> dict:
    t = s.thermo
    ranges = product if isinstance(product, list) else [list(product or (s.product_min, s.product_max))]
    args = {
        "PRIMER_TASK": "generic",
        "PRIMER_PICK_LEFT_PRIMER": 1,
        "PRIMER_PICK_RIGHT_PRIMER": 1,
        "PRIMER_PICK_INTERNAL_OLIGO": int(bool(s.pick_probe if pick_probe is None else pick_probe)),
        "PRIMER_OPT_SIZE": s.size_opt,
        "PRIMER_MIN_SIZE": s.size_min,
        "PRIMER_MAX_SIZE": min(35, s.size_max),
        "PRIMER_OPT_TM": s.tm_opt,
        "PRIMER_MIN_TM": s.tm_min,
        "PRIMER_MAX_TM": s.tm_max,
        "PRIMER_MIN_GC": s.gc_min,
        "PRIMER_MAX_GC": s.gc_max,
        "PRIMER_PAIR_MAX_DIFF_TM": s.max_tm_diff,
        "PRIMER_PRODUCT_SIZE_RANGE": [[int(a), int(b)] for a, b in ranges if b >= a],
        "PRIMER_NUM_RETURN": int(num_return or s.num_return),
        "PRIMER_MAX_POLY_X": s.max_poly_x,
        "PRIMER_GC_CLAMP": s.gc_clamp,
        "PRIMER_MAX_END_STABILITY": s.max_end_stability,
        "PRIMER_SALT_MONOVALENT": t.monovalent_mM,
        "PRIMER_SALT_DIVALENT": t.divalent_mM,
        "PRIMER_DNTP_CONC": t.dntp_mM,
        "PRIMER_DNA_CONC": t.oligo_nM,
        "PRIMER_TM_FORMULA": 1 if t.tm_method == "santalucia" else 0,
        "PRIMER_SALT_CORRECTIONS": {"schildkraut": 0, "santalucia": 1, "owczarzy": 2}[t.salt_method],
        "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT": 1,
        "PRIMER_EXPLAIN_FLAG": 1,
        "PRIMER_LIBERAL_BASE": 1,
        "PRIMER_MAX_NS_ACCEPTED": 0,
        # Pairs that differ by more than a base or two (Primer3 otherwise
        # returns the same primer shifted by one base several times).
        "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE": s.min_three_prime_distance,
        "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE": s.min_three_prime_distance,
    }
    if args["PRIMER_PICK_INTERNAL_OLIGO"]:
        args.update({
            "PRIMER_INTERNAL_OPT_SIZE": s.probe_size_opt,
            "PRIMER_INTERNAL_MIN_SIZE": s.probe_size_min,
            "PRIMER_INTERNAL_MAX_SIZE": min(35, s.probe_size_max),
            "PRIMER_INTERNAL_OPT_TM": s.probe_tm_opt,
            "PRIMER_INTERNAL_MIN_TM": s.probe_tm_min,
            "PRIMER_INTERNAL_MAX_TM": s.probe_tm_max,
            "PRIMER_INTERNAL_MIN_GC": 30.0,
            "PRIMER_INTERNAL_MAX_GC": 80.0,
            "PRIMER_INTERNAL_SALT_MONOVALENT": t.monovalent_mM,
            "PRIMER_INTERNAL_SALT_DIVALENT": t.divalent_mM,
            "PRIMER_INTERNAL_DNTP_CONC": t.dntp_mM,
            "PRIMER_INTERNAL_DNA_CONC": t.oligo_nM,
        })
    return args


def _explain(result: dict) -> list[str]:
    out = []
    for key, label in (("PRIMER_LEFT_EXPLAIN", "Left"), ("PRIMER_RIGHT_EXPLAIN", "Right"),
                       ("PRIMER_INTERNAL_EXPLAIN", "Probe"), ("PRIMER_PAIR_EXPLAIN", "Pairs")):
        text = result.get(key)
        if text and not str(text).startswith("considered 0"):
            out.append(f"{label}: {text}")
    return out


@dataclass
class Primer3Run:
    pairs: list[DesignedPair]
    explain: list[str]


def run_primer3(template_id: str, template: str, s: DesignSettings, *,
                excluded: list[tuple[int, int]] | None = None,
                target: tuple[int, int] | None = None,
                included: tuple[int, int] | None = None,
                ok_regions: list[list[int]] | None = None,
                left_primer: str | None = None,
                num_return: int | None = None,
                product: tuple[int, int] | list | None = None,
                pick_probe: bool | None = None,
                misprime_lib: dict | None = None,
                name_prefix: str = "",
                coord_map: Callable[[int, int, str], tuple[int, int, str]] | None = None,
                ) -> Primer3Run:
    """One Primer3 call. Intervals are 0-based half-open on ``template``.

    ``coord_map(start1, end1, strand)`` turns template coordinates into the
    coordinates reported to the user (for converted / excised templates).
    """
    p3 = require_primer3()
    if len(template) > MAX_TEMPLATE_BP:
        raise PrimerError(
            f"{template_id} is {len(template):,} bp. Primer design works on templates up to "
            f"{MAX_TEMPLATE_BP:,} bp: extract the region first (Sequence → Extract)."
        )
    seq_args: dict = {"SEQUENCE_ID": template_id, "SEQUENCE_TEMPLATE": template.upper()}
    if excluded:
        seq_args["SEQUENCE_EXCLUDED_REGION"] = [[a, b - a] for a, b in excluded if b > a]
    if target:
        seq_args["SEQUENCE_TARGET"] = [[target[0], max(1, target[1] - target[0])]]
    if included:
        seq_args["SEQUENCE_INCLUDED_REGION"] = [included[0], max(1, included[1] - included[0])]
    if ok_regions:
        seq_args["SEQUENCE_PRIMER_PAIR_OK_REGION_LIST"] = ok_regions
    if left_primer:
        seq_args["SEQUENCE_PRIMER"] = left_primer
    gargs = primer3_global_args(s, num_return=num_return, product=product, pick_probe=pick_probe)
    if left_primer:
        gargs["PRIMER_PICK_LEFT_PRIMER"] = 1
    try:
        res = p3.bindings.design_primers(seq_args, gargs, misprime_lib=misprime_lib or None)
    except (OSError, ValueError, RuntimeError) as exc:
        raise PrimerError(f"Primer3: {exc}") from exc
    pairs: list[DesignedPair] = []
    n = int(res.get("PRIMER_PAIR_NUM_RETURNED", 0))
    for i in range(n):
        lpos, llen = res[f"PRIMER_LEFT_{i}"]
        rpos, rlen = res[f"PRIMER_RIGHT_{i}"]
        lseq = res[f"PRIMER_LEFT_{i}_SEQUENCE"]
        rseq = res[f"PRIMER_RIGHT_{i}_SEQUENCE"]
        f_start, f_end, f_strand = lpos + 1, lpos + llen, "+"
        r_start, r_end, r_strand = rpos - rlen + 2, rpos + 1, "-"
        if coord_map is not None:
            f_start, f_end, f_strand = coord_map(f_start, f_end, "+")
            r_start, r_end, r_strand = coord_map(r_start, r_end, "-")
        tag = f"{name_prefix}{i + 1}"
        fwd = make_oligo(f"{template_id}_{tag}_F", lseq, "forward", start=f_start, end=f_end,
                         strand=f_strand, thermo=s.thermo,
                         penalty=float(res.get(f"PRIMER_LEFT_{i}_PENALTY", 0.0)),
                         tail=s.forward_tail)
        rev = make_oligo(f"{template_id}_{tag}_R", rseq, "reverse", start=r_start, end=r_end,
                         strand=r_strand, thermo=s.thermo,
                         penalty=float(res.get(f"PRIMER_RIGHT_{i}_PENALTY", 0.0)),
                         tail=s.reverse_tail)
        probe = None
        if f"PRIMER_INTERNAL_{i}_SEQUENCE" in res:
            ipos, ilen = res[f"PRIMER_INTERNAL_{i}"]
            p_start, p_end, p_strand = ipos + 1, ipos + ilen, "+"
            if coord_map is not None:
                p_start, p_end, p_strand = coord_map(p_start, p_end, "+")
            probe = make_oligo(f"{template_id}_{tag}_P", res[f"PRIMER_INTERNAL_{i}_SEQUENCE"],
                               "probe", start=p_start, end=p_end, strand=p_strand,
                               thermo=s.thermo,
                               penalty=float(res.get(f"PRIMER_INTERNAL_{i}_PENALTY", 0.0)))
        pair = DesignedPair(
            rank=i + 1, name=f"{template_id}_{tag}", template=template_id,
            forward=fwd, reverse=rev, probe=probe,
            product_size=int(res.get(f"PRIMER_PAIR_{i}_PRODUCT_SIZE", 0)),
            penalty=round(float(res.get(f"PRIMER_PAIR_{i}_PENALTY", 0.0)), 3),
        )
        pairs.append(finish_pair(pair, s.thermo))
    return Primer3Run(pairs=pairs, explain=_explain(res))


# --------------------------------------------------------------------------
# Ranking
# --------------------------------------------------------------------------


def rerank(pairs: list[DesignedPair], *, rename: str | None = None) -> list[DesignedPair]:
    """Number pairs 1, 2, 3 …; with ``rename`` (e.g. "LAD") names follow the rank."""
    for i, p in enumerate(pairs, start=1):
        p.rank = i
        if rename is not None:
            p.name = f"{p.template}_{rename}{i}"
            for o, tag in ((p.forward, "F"), (p.reverse, "R"), (p.probe, "P")):
                if o is not None:
                    o.name = f"{p.name}_{tag}"
    return pairs


# --------------------------------------------------------------------------
# Expected gel bands for a pair on its own template(s)
# --------------------------------------------------------------------------


def pair_bands_on(pair: DesignedPair, records: list[tuple[str, str]],
                  settings: bd.SearchSettings | None = None, max_size: int = 20000) -> list[int]:
    if pair.forward is None or pair.reverse is None:
        return []
    st = settings or bd.SearchSettings(max_mismatches=0, three_prime_exact=3)
    primers = [(pair.forward.name, pair.forward.sequence), (pair.reverse.name, pair.reverse.sequence)]
    sizes = []
    for rid, seq in records:
        sites = []
        for name, pseq in primers:
            sites.extend(bd.find_sites_in_sequence(name, pseq, rid, seq, st, exhaustive=True))
        for pr in bd.products_for_pair(sites, pair.forward.name, pair.reverse.name,
                                       max_size=max_size):
            sizes.append(pr.size)
    return sizes


# --------------------------------------------------------------------------
# Satellite ladder PCR and FISH probes
# --------------------------------------------------------------------------


def design_ladder(rid: str, monomer: str, s: DesignSettings, misprime: dict | None = None
                  ) -> Primer3Run:
    """Outward-facing primers in one monomer: the forward primer in copy 1 and the
    reverse in copy 2 of a dimer, so a tandem array gives d, d + n, d + 2n … bands."""
    n = len(monomer)
    if n < s.size_min * 2 + 10:
        raise PrimerError(f"{rid} is only {n} bp: too short for ladder primers.")
    dimer = (monomer + monomer).upper()
    lo = max(s.product_min, s.size_min * 2)
    # ≤ one monomer: the two primers then never overlap each other
    hi = min(s.product_max, n)
    if hi < lo:
        lo, hi = s.size_min * 2, n

    def mono_coords(a: int, b: int, strand: str) -> tuple[int, int, str]:
        a0 = (a - 1) % n + 1
        return a0, a0 + (b - a), strand

    run = run_primer3(rid, dimer, s, ok_regions=[[0, n, n, n]], product=(lo, hi),
                      misprime_lib=misprime, name_prefix="LAD", coord_map=mono_coords,
                      num_return=s.num_return)
    run.pairs = [p for p in run.pairs if p.dimer_dg >= -9.0] or run.pairs
    for p in run.pairs:
        p.extra["monomer_length"] = n
        p.notes.append(f"smallest band {p.product_size} bp, then +{n} bp per extra monomer")
    return run


def ladder_bands(product: int, monomer: int, copies: int = 8) -> list[tuple[int, float]]:
    """Bands from an array of ``copies`` monomers: (size, relative amount)."""
    out = []
    for k in range(copies - 1):
        size = product + k * monomer
        out.append((size, (copies - 1 - k) / (copies - 1)))
    return out


def design_fish(rid: str, template: str, members: list[tuple[str, str]], s: DesignSettings,
                avoid: list[tuple[str, str]] | None = None,
                progress: Callable[[str], None] | None = None) -> list[DesignedPair]:
    """Oligo probes conserved across all monomers/copies, without hairpins."""
    th = s.thermo
    t = template.upper()
    lengths = sorted({s.size_min, s.size_opt, s.size_max})
    sep = "N" * (s.size_max + 2)
    pool = sep.join(seq.upper() for _n, seq in members)
    pool_mask = bd.encode_template(pool)
    starts = np.cumsum([0] + [len(seq) + len(sep) for _n, seq in members])
    avoid_pool = sep.join(seq.upper() for _n, seq in (avoid or []))
    avoid_mask = bd.encode_template(avoid_pool) if avoid_pool else None
    st = bd.SearchSettings(max_mismatches=s.coverage_mismatches, three_prime_exact=0)
    cands = []
    for L in lengths:
        if progress is not None:
            progress(f"Scoring {L}-nt probes")
        for i in range(0, len(t) - L + 1):
            seq = t[i:i + L]
            if "N" in seq:
                continue
            gc = gc_percent(seq)
            if not (s.gc_min <= gc <= s.gc_max):
                continue
            cands.append((i, L, seq, gc))
    scored = []
    for i, L, seq, gc in cands:
        hit_rows = set()
        for strand_seq in (seq, reverse_complement(seq)):
            idx, _m = bd._score_positions(pool_mask, bd.encode_primer(strand_seq), st, True)
            if len(idx):
                rows = np.searchsorted(starts, idx, side="right") - 1
                hit_rows.update(rows.tolist())
        cov = len(hit_rows) / max(1, len(members))
        if cov < s.min_coverage:
            continue
        cross = 0
        if avoid_mask is not None:
            for strand_seq in (seq, reverse_complement(seq)):
                idx, _m = bd._score_positions(avoid_mask, bd.encode_primer(strand_seq), st, True)
                cross += len(idx)
        scored.append((i, L, seq, gc, cov, cross))
    scored.sort(key=lambda x: (x[5] > 0, -x[4], abs(x[3] - 48)))
    pairs: list[DesignedPair] = []
    taken: list[tuple[int, int]] = []
    for i, L, seq, gc, cov, cross in scored:
        if any(not (i + L <= a or i >= b) for a, b in taken):
            continue
        tm = calc_tm(seq, th)
        if not (s.tm_min <= tm <= s.tm_max):
            continue
        oligo = make_oligo(f"{rid}_FISH{len(pairs) + 1}", seq, "probe", start=i + 1,
                           end=i + L, strand="+", thermo=th)
        if oligo.hairpin_dg < -6.0 or oligo.self_dimer_dg < -12.0:
            continue
        k = len(pairs) + 1
        label = "" if s.fish_label == "None" else s.fish_label
        pair = DesignedPair(rank=k, name=f"{rid}_FISH{k}", template=rid, probe=oligo,
                            penalty=round((1 - cov) * 20 + (5 if cross else 0), 2))
        pair.coverage = {"target_total": len(members), "target_amplified": int(round(cov * len(members))),
                         "cross_hits": cross}
        pair.extra["label"] = label
        pair.extra["order"] = (f"5'-{label}-" if label else "5'-") + seq + "-3'"
        pair.notes.append(f"in {pair.coverage['target_amplified']}/{len(members)} copies"
                          + (f"; {cross} hits in the repeats to avoid" if cross else ""))
        pairs.append(pair)
        taken.append((i, i + L))
        if len(pairs) >= s.num_return:
            break
    return pairs


def representative_record(records: list[tuple[str, str]], k: int = 8,
                          sample: int = 60) -> tuple[str, str]:
    """The copy most like the others (shares the most k-mers): a stand-in
    consensus for unaligned monomers / copies."""
    if len(records) <= 2:
        return records[0]
    pool = records[:sample]

    def kmers(seq: str) -> set[str]:
        s = seq.upper()
        return {s[i:i + k] for i in range(len(s) - k + 1)}

    sets = [kmers(q) for _n, q in pool]
    best, best_score = 0, -1.0
    for i, a in enumerate(sets):
        score = sum(len(a & b) / max(1, len(a | b)) for j, b in enumerate(sets) if j != i)
        if score > best_score:
            best, best_score = i, score
    return pool[best]
