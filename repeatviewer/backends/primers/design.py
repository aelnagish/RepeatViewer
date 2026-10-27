"""Primer design: Primer3 plus the repeat-aware modes.

Modes (``DesignSettings.mode``):

* ``pcr``        – standard PCR pair (Primer3), repeat-aware masking.
* ``qpcr``       – qPCR pair (+ optional hydrolysis probe), and a matching pair
                   on a single-copy reference gene for copy-number ΔCt.
* ``bsp``        – bisulfite PCR: primers on the converted strand, no CpG inside.
* ``msp``        – methylation-specific PCR: M and U pairs whose 3' ends sit on CpGs.
* ``conserved``  – primers that catch a whole TE / repeat family (from an MSA),
                   degenerate where needed.
* ``specific``   – primers that catch one subfamily of an MSA and not the rest.
* ``ladder``     – satellite ladder PCR: outward-facing primers in one monomer,
                   giving 1×, 2×, 3× … monomer bands on a tandem array.
* ``fish``       – FISH oligo probes conserved across monomers / copies.
* ``insertion``  – TE insertion genotyping: flanking pair + a primer inside the TE.

Everything here is Qt-free; the dialogs only fill :class:`DesignSettings`.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import asdict, dataclass, field

import numpy as np

from repeatviewer.backends.primers import binding as bd
from repeatviewer.backends.primers.thermo import (
    PrimerError,
    ThermoSettings,
    analyze_oligo,
    binding_part,
    calc_tm,
    degeneracy,
    expand_iupac,
    gc_percent,
    heterodimer,
    iupac_code,
    require_primer3,
    reverse_complement,
)

MODE_PCR = "pcr"
MODE_QPCR = "qpcr"
MODE_BSP = "bsp"
MODE_MSP = "msp"
MODE_CONSERVED = "conserved"
MODE_SPECIFIC = "specific"
MODE_LADDER = "ladder"
MODE_FISH = "fish"
MODE_INSERTION = "insertion"
DESIGN_MODES = (MODE_PCR, MODE_QPCR, MODE_BSP, MODE_MSP, MODE_CONSERVED, MODE_SPECIFIC,
                MODE_LADDER, MODE_FISH, MODE_INSERTION)
MODE_LABELS = {
    MODE_PCR: "Standard PCR",
    MODE_QPCR: "qPCR (copy number)",
    MODE_BSP: "Bisulfite PCR (BSP)",
    MODE_MSP: "Methylation-specific PCR (MSP)",
    MODE_CONSERVED: "Conserved in a family (MSA)",
    MODE_SPECIFIC: "Subfamily-specific (MSA)",
    MODE_LADDER: "Satellite ladder PCR",
    MODE_FISH: "FISH oligo probes",
    MODE_INSERTION: "TE insertion (presence / absence)",
}
MSA_MODES = (MODE_CONSERVED, MODE_SPECIFIC)
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
    # Regions (1-based inclusive; 0 = not set)
    target_start: int = 0
    target_end: int = 0
    included_start: int = 0
    included_end: int = 0
    # Repeat-aware masking
    mask_lowercase: bool = True
    mask_features: bool = True
    mask_feature_types: str = ""          # comma list; empty = every type
    mask_tandem: bool = True
    mask_file: str = ""                   # GFF / BED / RepeatMasker .out
    mispriming_library: str = ""          # FASTA of repeats (Primer3 mispriming)
    # Genome specificity
    check_genome: bool = False
    genome_fasta: str = ""
    genome_max_mismatches: int = 2
    genome_three_prime_exact: int = 3
    genome_seed: int = bd.DEFAULT_SEED
    genome_max_product: int = 3000
    genome_max_products: int = 1          # drop pairs making more products (0 = keep all)
    drop_unspecific: bool = True
    # Tails
    forward_tail: str = ""
    reverse_tail: str = ""
    # qPCR
    reference_fasta: str = ""
    # Bisulfite
    bisulfite_strand: str = "top"          # top / bottom
    min_converted: int = 3
    # MSA modes
    min_coverage: float = 0.90
    max_degeneracy: int = 16
    degenerate_min_freq: float = 0.10
    coverage_mismatches: int = 1
    coverage_three_prime_exact: int = 3
    target_pattern: str = ""
    max_off_target: float = 0.10
    # FISH
    fish_label: str = "Cy3"
    fish_avoid_fasta: str = ""
    # TE insertion
    te_start: int = 0
    te_end: int = 0
    te_feature: str = ""
    flank_max: int = 400

    def __post_init__(self) -> None:
        if isinstance(self.thermo, dict):
            self.thermo = ThermoSettings.from_dict(self.thermo)
        if self.mode not in DESIGN_MODES:
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

    def search_settings(self) -> bd.SearchSettings:
        return bd.SearchSettings(max_mismatches=self.genome_max_mismatches,
                                 three_prime_exact=self.genome_three_prime_exact)


def mode_preset(mode: str) -> dict:
    """Condition values the dialog fills in when a mode is chosen."""
    if mode == MODE_QPCR:
        return dict(size_min=18, size_opt=20, size_max=24, tm_min=58.0, tm_opt=60.0,
                    tm_max=62.0, gc_min=40.0, gc_max=60.0, max_tm_diff=1.0,
                    product_min=70, product_max=180, gc_clamp=0, max_poly_x=3)
    if mode in (MODE_BSP, MODE_MSP):
        return dict(size_min=20, size_opt=25, size_max=30, tm_min=52.0, tm_opt=58.0,
                    tm_max=64.0, gc_min=15.0, gc_max=60.0, max_tm_diff=3.0,
                    product_min=100, product_max=350, gc_clamp=0, max_poly_x=6)
    if mode in MSA_MODES:
        return dict(coverage_mismatches=1, min_coverage=0.9, size_min=18, size_opt=21, size_max=26, tm_min=54.0, tm_opt=60.0,
                    tm_max=66.0, gc_min=30.0, gc_max=70.0, max_tm_diff=4.0,
                    product_min=100, product_max=1500)
    if mode == MODE_LADDER:
        return dict(coverage_mismatches=1, size_min=18, size_opt=20, size_max=25, tm_min=56.0, tm_opt=60.0,
                    tm_max=64.0, gc_min=30.0, gc_max=70.0, max_tm_diff=2.0,
                    product_min=40, product_max=2000)
    if mode == MODE_FISH:
        return dict(size_min=30, size_opt=36, size_max=45, tm_min=60.0, tm_opt=72.0,
                    tm_max=85.0, gc_min=35.0, gc_max=65.0, coverage_mismatches=3,
                    min_coverage=0.9)
    if mode == MODE_INSERTION:
        return dict(size_min=18, size_opt=22, size_max=27, tm_min=57.0, tm_opt=60.0,
                    tm_max=64.0, gc_min=35.0, gc_max=65.0, max_tm_diff=2.0,
                    product_min=150, product_max=700)
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
# Repeat-aware masking
# --------------------------------------------------------------------------


def merge_intervals(intervals: list[tuple[int, int]], gap: int = 0) -> list[tuple[int, int]]:
    """0-based half-open intervals, merged."""
    out: list[list[int]] = []
    for s, e in sorted((int(a), int(b)) for a, b in intervals if b > a):
        if out and s <= out[-1][1] + gap:
            out[-1][1] = max(out[-1][1], e)
        else:
            out.append([s, e])
    return [(a, b) for a, b in out]


def lowercase_intervals(seq: str) -> list[tuple[int, int]]:
    out = []
    start = None
    for i, ch in enumerate(seq):
        low = ch.islower()
        if low and start is None:
            start = i
        elif not low and start is not None:
            out.append((start, i))
            start = None
    if start is not None:
        out.append((start, len(seq)))
    return out


def tandem_intervals(record_id: str, seq: str) -> list[tuple[int, int]]:
    try:
        from repeatviewer.backends.tandem_repeat_detect import detect_tandem_repeat_features

        res = detect_tandem_repeat_features(record_id, seq)
    except Exception:  # noqa: BLE001 - pytrf missing or failed: no tandem mask
        return []
    return merge_intervals([(f.start - 1, f.end) for f in res.features], gap=5)


def file_intervals(path: str, record_id: str) -> list[tuple[int, int]]:
    """Repeat intervals from a GFF / BED / RepeatMasker .out file for one record."""
    from pathlib import Path

    p = Path(path).expanduser()
    if not p.is_file():
        return []
    out = []
    name = p.name.lower()
    with open(p, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.split()
            try:
                if name.endswith(".out") and len(cols) >= 7 and cols[0].isdigit():
                    rid, s, e = cols[4], int(cols[5]) - 1, int(cols[6])
                elif name.endswith(".bed") and len(cols) >= 3:
                    rid, s, e = cols[0], int(cols[1]), int(cols[2])
                else:
                    cols = line.rstrip("\n").split("\t")
                    if len(cols) < 5:
                        continue
                    rid, s, e = cols[0], int(cols[3]) - 1, int(cols[4])
            except ValueError:
                continue
            if rid == record_id or rid.split()[0] == record_id:
                out.append((s, e))
    return merge_intervals(out)


def mask_intervals(record_id: str, seq: str, settings: DesignSettings,
                   features: list | None = None) -> tuple[list[tuple[int, int]], dict]:
    """Regions primers must not overlap (0-based half-open) and what they came from."""
    sources: dict[str, int] = {}
    parts: list[tuple[int, int]] = []

    def add(label: str, ivs: list[tuple[int, int]]) -> None:
        if ivs:
            sources[label] = sum(e - s for s, e in merge_intervals(ivs))
            parts.extend(ivs)

    if settings.mask_lowercase:
        add("soft-masked (lower case)", lowercase_intervals(seq))
    if settings.mask_features and features:
        wanted = {t.strip().lower() for t in settings.mask_feature_types.split(",") if t.strip()}
        skip = {"primer_bind", "primer", "pcr product", "probe", "orf", "cds", "gene", "exon",
                "mrna", "domain"}
        ivs = []
        for f in features:
            ftype = str(getattr(f, "feature_type", "") or "").lower()
            if wanted:
                if ftype not in wanted:
                    continue
            elif ftype in skip:
                continue
            ivs.append((int(f.start) - 1, int(f.end)))
        add("annotated repeats", ivs)
    if settings.mask_tandem:
        add("tandem repeats", tandem_intervals(record_id, seq.upper()))
    if settings.mask_file:
        add("repeat annotation file", file_intervals(settings.mask_file, record_id))
    return merge_intervals(parts), sources


def load_mispriming_library(path: str, limit_bp: int = 5_000_000) -> dict[str, str]:
    from pathlib import Path

    p = Path(path).expanduser() if path else None
    if p is None or not p.is_file():
        return {}
    lib: dict[str, str] = {}
    total = 0
    for rid, seq in bd.iter_fasta(p):
        s = "".join(ch for ch in seq.upper() if ch in "ACGTN")
        if not s:
            continue
        lib[f"{rid}_{len(lib)}"] = s
        total += len(s)
        if total > limit_bp:
            break
    return lib


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
# Genome specificity
# --------------------------------------------------------------------------


def check_genome(pairs: list[DesignedPair], s: DesignSettings, *,
                 progress: Callable[[str], None] | None = None,
                 cancel: Callable[[], bool] | None = None,
                 keep_products: int = 50) -> dict:
    """Search every primer of ``pairs`` in the genome; count sites and products.

    Fills ``pair.specificity`` = {sites_forward, sites_reverse, products,
    product_sizes (first ``keep_products``), unique}. Returns a summary.
    """
    from pathlib import Path

    gpath = Path(s.genome_fasta).expanduser()
    if not gpath.is_file():
        raise PrimerError(f"Genome FASTA not found: {s.genome_fasta}")
    primers: dict[str, str] = {}
    for p in pairs:
        for o in p.oligos():
            if o.kind in ("forward", "reverse") or (o.kind == "probe" and p.forward is None):
                primers[o.name] = binding_part(o.sequence)
    if not primers:
        return {}
    settings = s.search_settings()
    large = gpath.stat().st_size > bd.EXHAUSTIVE_MAX_BP
    site_counts: dict[str, int] = {n: 0 for n in primers}
    products: dict[str, list[tuple[str, int, int, int]]] = {p.name: [] for p in pairs}
    product_counts: dict[str, int] = {p.name: 0 for p in pairs}
    total_bp = 0
    n_records = 0

    def on_record(rid: str, seq: str, sites: list[bd.BindingSite]) -> None:
        nonlocal total_bp, n_records
        total_bp += len(seq)
        n_records += 1
        for site in sites:
            site_counts[site.primer] = site_counts.get(site.primer, 0) + 1
        for p in pairs:
            if p.forward is None or p.reverse is None:
                continue
            prods = bd.products_for_pair(sites, p.forward.name, p.reverse.name,
                                         max_size=s.genome_max_product)
            product_counts[p.name] += len(prods)
            room = keep_products - len(products[p.name])
            for pr in prods[:max(0, room)]:
                products[p.name].append((rid, pr.start, pr.end, pr.size))

    bd.search_records(list(primers.items()), bd.iter_fasta(gpath), settings, large=large,
                      seed=s.genome_seed, progress=progress, cancel=cancel,
                      on_record=on_record)
    for p in pairs:
        spec = {
            "genome": gpath.name,
            "sites_forward": site_counts.get(p.forward.name, 0) if p.forward else 0,
            "sites_reverse": site_counts.get(p.reverse.name, 0) if p.reverse else 0,
            "sites_probe": site_counts.get(p.probe.name, 0) if p.probe else 0,
            "products": product_counts.get(p.name, 0),
            "product_list": [list(x) for x in products.get(p.name, [])],
        }
        if p.forward is None and p.probe is not None:
            spec["sites_probe"] = site_counts.get(p.probe.name, 0)
        spec["unique"] = spec["products"] == 1
        p.specificity = spec
        if p.forward is not None:
            if spec["products"] == 0:
                p.notes.append(f"no product in {gpath.name}")
            elif spec["products"] > 1:
                p.notes.append(f"{spec['products']:,} products in {gpath.name}")
    return {"genome": str(gpath), "records": n_records, "bp": total_bp,
            "seeded": large, "seed": s.genome_seed if large else 0}


def drop_unspecific(pairs: list[DesignedPair], s: DesignSettings) -> tuple[list[DesignedPair], int]:
    if not (s.check_genome and s.drop_unspecific and s.genome_max_products > 0):
        return pairs, 0
    kept = [p for p in pairs if not p.specificity
            or p.specificity.get("products", 0) <= s.genome_max_products]
    return kept, len(pairs) - len(kept)


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
# Bisulfite (BSP / MSP)
# --------------------------------------------------------------------------


def bisulfite_convert(seq: str, *, keep_cpg: bool) -> str:
    """Top-strand bisulfite conversion: C → T, CpG C kept when ``keep_cpg``."""
    s = seq.upper()
    out = []
    for i, ch in enumerate(s):
        if ch == "C":
            if keep_cpg and i + 1 < len(s) and s[i + 1] == "G":
                out.append("C")
            else:
                out.append("T")
        else:
            out.append(ch)
    return "".join(out)


def cpg_positions(seq: str) -> list[int]:
    s = seq.upper()
    return [i for i in range(len(s) - 1) if s[i] == "C" and s[i + 1] == "G"]


def _strand_template(seq: str, strand: str) -> tuple[str, Callable]:
    """The strand to convert, and a map back to original top-strand coordinates."""
    if strand == "bottom":
        rc = reverse_complement(seq.upper())
        n = len(seq)

        def back(a: int, b: int, primer_strand: str) -> tuple[int, int, str]:
            return n - b + 1, n - a + 1, "-" if primer_strand == "+" else "+"
        return rc, back

    def same(a: int, b: int, primer_strand: str) -> tuple[int, int, str]:
        return a, b, primer_strand
    return seq.upper(), same


def converted_c_count(original: str, start0: int, end0: int) -> int:
    """Non-CpG C in the original strand inside [start0, end0): they become T."""
    s = original.upper()
    return sum(1 for i in range(start0, end0)
               if s[i] == "C" and not (i + 1 < len(s) and s[i + 1] == "G"))


def design_bsp(rid: str, seq: str, s: DesignSettings, excluded: list[tuple[int, int]],
               misprime: dict | None = None) -> Primer3Run:
    strand_seq, back = _strand_template(seq, s.bisulfite_strand)
    conv = bisulfite_convert(strand_seq, keep_cpg=True)
    ex = list(excluded) if s.bisulfite_strand == "top" else [
        (len(seq) - b, len(seq) - a) for a, b in excluded]
    ex += [(i, i + 2) for i in cpg_positions(strand_seq)]
    run = run_primer3(rid, conv, s, excluded=merge_intervals(ex), num_return=max(40, s.num_return * 8),
                      misprime_lib=misprime, name_prefix="BSP", coord_map=back)
    kept = []
    for p in run.pairs:
        f0 = p.forward.start - 1 if s.bisulfite_strand == "top" else len(seq) - p.forward.end
        r0 = p.reverse.start - 1 if s.bisulfite_strand == "top" else len(seq) - p.reverse.end
        fc = converted_c_count(strand_seq, f0, f0 + p.forward.length)
        # a reverse primer gains specificity from converted C on the strand it copies
        rc = converted_c_count(strand_seq, r0, r0 + p.reverse.length)
        p.extra.update({"converted_c_forward": fc, "converted_c_reverse": rc,
                        "strand": s.bisulfite_strand})
        if fc < s.min_converted or rc < s.min_converted:
            continue
        p.notes.append(f"{fc}/{rc} converted C (F/R); no CpG in primers")
        kept.append(p)
    kept = rerank(kept[: s.num_return])
    for p in kept:
        p.name = f"{rid}_BSP{p.rank}"
        p.forward.name, p.reverse.name = f"{p.name}_F", f"{p.name}_R"
    explain = run.explain
    if run.pairs and not kept:
        explain = explain + [f"Primer3 found {len(run.pairs)} pairs, but none had ≥ "
                             f"{s.min_converted} converted C in both primers."]
    return Primer3Run(pairs=kept, explain=explain)


def design_msp(rid: str, seq: str, s: DesignSettings) -> list[DesignedPair]:
    """MSP: an M pair and a U pair at the same place, 3' ends on CpG C."""
    strand_seq, back = _strand_template(seq, s.bisulfite_strand)
    tm_seq = bisulfite_convert(strand_seq, keep_cpg=True)
    tu_seq = bisulfite_convert(strand_seq, keep_cpg=False)
    cpgs = cpg_positions(strand_seq)
    cpg_set = set(cpgs)
    th = s.thermo
    fwd_c: list[tuple[float, int, int, str, float]] = []
    rev_c: list[tuple[float, int, int, str, float]] = []
    n = len(strand_seq)
    for c in cpgs:
        for L in range(s.size_min, s.size_max + 1):
            a = c - L + 1
            if a >= 0:
                f = tm_seq[a:c + 1]
                tm = calc_tm(f, th)
                if s.tm_min <= tm <= s.tm_max and s.gc_min <= gc_percent(f) <= s.gc_max:
                    k = sum(1 for i in range(a, c + 1) if i in cpg_set)
                    conv = converted_c_count(strand_seq, a, c + 1)
                    pen = abs(tm - s.tm_opt) + abs(L - s.size_opt) * 0.3 - min(k, 3) * 0.8 \
                        - min(conv, 5) * 0.4
                    fwd_c.append((pen, a, c + 1, f, tm))
            b = c + L
            if b <= n:
                r = reverse_complement(tm_seq[c:b])
                tm = calc_tm(r, th)
                if s.tm_min <= tm <= s.tm_max and s.gc_min <= gc_percent(r) <= s.gc_max:
                    k = sum(1 for i in range(c, b) if i in cpg_set)
                    conv = converted_c_count(strand_seq, c, b)
                    pen = abs(tm - s.tm_opt) + abs(L - s.size_opt) * 0.3 - min(k, 3) * 0.8 \
                        - min(conv, 5) * 0.4
                    rev_c.append((pen, c, b, r, tm))
    fwd_c.sort()
    rev_c.sort()
    fwd_c, rev_c = fwd_c[:300], rev_c[:300]
    combos = []
    for fp in fwd_c:
        for rp in rev_c:
            size = rp[2] - fp[1]
            if rp[1] < fp[2] or not (s.product_min <= size <= s.product_max):
                continue
            if abs(fp[4] - rp[4]) > s.max_tm_diff:
                continue
            combos.append((fp[0] + rp[0] + abs(fp[4] - rp[4]) * 0.5, fp, rp, size))
    combos.sort(key=lambda x: x[0])
    pairs: list[DesignedPair] = []
    used: set[tuple[int, int]] = set()
    for pen, fp, rp, size in combos:
        if (fp[1], rp[2]) in used or len(pairs) >= 2 * s.num_return:
            continue
        used.add((fp[1], rp[2]))
        k = len(pairs) // 2 + 1
        for group, conv in (("M", tm_seq), ("U", tu_seq)):
            fseq = conv[fp[1]:fp[2]]
            rseq = reverse_complement(conv[rp[1]:rp[2]])
            fs, fe, fstr = back(fp[1] + 1, fp[2], "+")
            rs, re_, rstr = back(rp[1] + 1, rp[2], "-")
            name = f"{rid}_MSP{k}{group}"
            pair = DesignedPair(
                rank=k, name=name, template=rid,
                forward=make_oligo(f"{name}_F", fseq, "forward", start=fs, end=fe, strand=fstr,
                                   thermo=th, tail=s.forward_tail),
                reverse=make_oligo(f"{name}_R", rseq, "reverse", start=rs, end=re_, strand=rstr,
                                   thermo=th, tail=s.reverse_tail),
                product_size=size, penalty=round(pen, 2), group=group,
            )
            n_cpg = sum(1 for i in list(range(fp[1], fp[2])) + list(range(rp[1], rp[2]))
                        if i in cpg_set)
            pair.notes.append(("methylated" if group == "M" else "unmethylated")
                              + f" template; {n_cpg} CpG under the primers")
            pair.extra["strand"] = s.bisulfite_strand
            pairs.append(finish_pair(pair, th))
    return pairs


# --------------------------------------------------------------------------
# MSA modes: conserved family primers and subfamily-specific primers
# --------------------------------------------------------------------------


@dataclass
class _Window:
    strand: str            # "+" forward, "-" reverse
    start: int             # 0-based consensus start
    length: int
    seq: str               # primer sequence (5'→3')
    tm: float
    gc: float
    deg: int
    hit: np.ndarray        # bool per alignment row: primer binds
    penalty: float

    @property
    def end(self) -> int:
        return self.start + self.length


def _match_matrix(M: np.ndarray, masks: np.ndarray) -> np.ndarray:
    return (M & masks[None, :]) != 0


def _window_ok(match: np.ndarray, L: int, max_mm: int, exact: int, at_end: bool) -> np.ndarray:
    """(n_rows, n_windows) True where the window binds with ≤ max_mm, exact 3' end."""
    n_rows, P = match.shape
    nw = P - L + 1
    if nw <= 0:
        return np.zeros((n_rows, 0), dtype=bool)
    cs = np.concatenate([np.zeros((n_rows, 1), dtype=np.int32),
                         np.cumsum(match, axis=1, dtype=np.int32)], axis=1)
    hits = cs[:, L:L + nw] - cs[:, :nw]
    ok = (L - hits) <= max_mm
    k = min(exact, L)
    if k:
        if at_end:
            seg = cs[:, L:L + nw] - cs[:, L - k:L - k + nw]
        else:
            seg = cs[:, k:k + nw] - cs[:, :nw]
        ok &= seg == k
    return ok


def _bits(mask: int) -> str:
    return "".join(b for b, bit in zip("ACGT", (1, 2, 4, 8)) if mask & bit)


def design_msa(names: list[str], rows: list[str], s: DesignSettings,
               progress: Callable[[str], None] | None = None) -> tuple[list[DesignedPair], dict]:
    """Conserved (whole family) or subfamily-specific primers from an alignment."""
    if len(rows) < 2:
        raise PrimerError("These modes need an alignment of at least two sequences.")
    masks = bd.alignment_masks([r.upper() for r in rows])
    n_rows = masks.shape[0]
    specific = s.mode == MODE_SPECIFIC
    if specific:
        pat = s.target_pattern.strip()
        if not pat:
            raise PrimerError("Type the target sequences' names (or part of them) in "
                              "\"Target sequences\".")
        import re

        try:
            rx = re.compile(pat, re.I)
            target = np.array([bool(rx.search(n)) for n in names])
        except re.error:
            target = np.array([pat.lower() in n.lower() for n in names])
        if target.sum() == 0:
            raise PrimerError(f"No sequence name matches \"{pat}\".")
        if target.sum() == n_rows:
            raise PrimerError("Every sequence matches the target names: nothing to tell apart.")
    else:
        target = np.ones(n_rows, dtype=bool)
    other = ~target
    gap_frac = (masks == 0).mean(axis=0)
    cols = np.nonzero(gap_frac <= 0.5)[0]
    if len(cols) < s.size_min * 2:
        raise PrimerError("The alignment has too few well-covered columns for primers.")
    M = masks[:, cols]
    Mt = M[target]
    # per-column base frequencies in the target rows
    counts = bd.base_counts(Mt)[:, :4]
    totals = np.maximum(counts.sum(axis=1), 1e-9)
    freq = counts / totals[:, None]
    major = np.argmax(freq, axis=1)
    cons_mask = np.array([(1, 2, 4, 8)[i] for i in major], dtype=np.uint8)
    deg_mask = np.zeros(len(cols), dtype=np.uint8)
    for j, bit in enumerate((1, 2, 4, 8)):
        deg_mask |= np.where(freq[:, j] >= s.degenerate_min_freq, bit, 0).astype(np.uint8)
    deg_mask = np.where(deg_mask == 0, cons_mask, deg_mask)
    deg_bits = np.array([bin(int(m)).count("1") for m in deg_mask], dtype=np.float64)
    log_deg = np.concatenate([[0.0], np.cumsum(np.log2(deg_bits))])
    cons_seq = "".join("ACGT"[i] for i in major)
    deg_seq = "".join(iupac_code(_bits(int(m))) for m in deg_mask)
    match_c = _match_matrix(M, cons_mask)
    match_d = _match_matrix(M, deg_mask)
    mm = s.coverage_mismatches
    ex = s.coverage_three_prime_exact
    th = s.thermo
    windows: dict[str, list[_Window]] = {"+": [], "-": []}
    P = len(cols)
    max_deg_log = math.log2(max(1, s.max_degeneracy))
    for L in range(s.size_min, min(s.size_max, P) + 1):
        if progress is not None:
            progress(f"Scoring {L}-mers across {n_rows} sequences")
        for strand, at_end in (("+", True), ("-", False)):
            ok_c = _window_ok(match_c, L, mm, ex, at_end)
            ok_d = _window_ok(match_d, L, mm, ex, at_end)
            nw = ok_c.shape[1]
            cov_c = ok_c[target].mean(axis=0)
            cov_d = ok_d[target].mean(axis=0)
            wdeg = log_deg[L:L + nw] - log_deg[:nw]
            use_d = (cov_d > cov_c + 1e-9) & (wdeg <= max_deg_log)
            cov = np.where(use_d, cov_d, cov_c)
            cand = np.nonzero(cov >= s.min_coverage)[0]
            if specific and len(cand):
                off = np.where(use_d[cand], ok_d[other][:, cand].mean(axis=0),
                               ok_c[other][:, cand].mean(axis=0))
                cand = cand[off <= max(s.max_off_target, 0.5)]
            for i in cand.tolist():
                d = bool(use_d[i])
                raw = (deg_seq if d else cons_seq)[i:i + L]
                seq = raw if strand == "+" else reverse_complement(raw)
                gc = gc_percent(seq)
                if not (s.gc_min <= gc <= s.gc_max):
                    continue
                variants = expand_iupac(seq, 8)
                tm = sum(calc_tm(v, th) for v in variants) / len(variants)
                if not (s.tm_min <= tm <= s.tm_max):
                    continue
                hit = (ok_d if d else ok_c)[:, i].copy()
                c = float(hit[target].mean())
                pen = abs(tm - s.tm_opt) + abs(L - s.size_opt) * 0.3 + (1 - c) * 25 \
                    + (2 ** float(wdeg[i]) - 1 if d else 0) * 0.25
                if specific:
                    pen += float(hit[other].mean()) * 25
                windows[strand].append(_Window(strand, i, L, seq, round(tm, 2), round(gc, 1),
                                               degeneracy(seq), hit, pen))
    for strand in windows:
        windows[strand].sort(key=lambda w: w.penalty)
        # keep the best window per start position, then the best few hundred
        seen: set[int] = set()
        best = []
        for w in windows[strand]:
            key = w.start if strand == "+" else w.end
            if key in seen:
                continue
            seen.add(key)
            best.append(w)
        windows[strand] = best[:250]
    combos = []
    for f in windows["+"]:
        for r in windows["-"]:
            size = r.end - f.start
            if r.start < f.end or not (s.product_min <= size <= s.product_max):
                continue
            if abs(f.tm - r.tm) > s.max_tm_diff:
                continue
            both = f.hit & r.hit
            cov_t = float(both[target].mean())
            if cov_t < s.min_coverage * 0.95:
                continue
            pen = f.penalty + r.penalty + (1 - cov_t) * 40 + abs(f.tm - r.tm) * 0.5
            cov_o = float(both[other].mean()) if other.any() else 0.0
            if specific:
                if cov_o > s.max_off_target:
                    continue
                pen += cov_o * 60
            combos.append((pen, f, r, size, cov_t, cov_o, both))
    combos.sort(key=lambda x: x[0])
    pairs: list[DesignedPair] = []
    used_f: dict[int, int] = {}
    used_r: dict[int, int] = {}
    label = "SPEC" if specific else "CONS"
    for pen, f, r, size, cov_t, cov_o, both in combos:
        if used_f.get(f.start, 0) >= 2 or used_r.get(r.end, 0) >= 2:
            continue
        used_f[f.start] = used_f.get(f.start, 0) + 1
        used_r[r.end] = used_r.get(r.end, 0) + 1
        k = len(pairs) + 1
        name = f"{label}{k}"
        pair = DesignedPair(
            rank=k, name=name, template="consensus",
            forward=make_oligo(f"{name}_F", f.seq, "forward", start=f.start + 1, end=f.end,
                               strand="+", thermo=th, tail=s.forward_tail),
            reverse=make_oligo(f"{name}_R", r.seq, "reverse", start=r.start + 1, end=r.end,
                               strand="-", thermo=th, tail=s.reverse_tail),
            product_size=size, penalty=round(pen, 2),
        )
        n_t = int(target.sum())
        n_o = int(other.sum())
        pair.coverage = {
            "target_total": n_t, "target_amplified": int(both[target].sum()),
            "other_total": n_o, "other_amplified": int(both[other].sum()) if n_o else 0,
            "forward_binds": int(f.hit[target].sum()), "reverse_binds": int(r.hit[target].sum()),
            "missed": [names[i] for i in np.nonzero(target & ~both)[0][:20].tolist()],
            "off_target": [names[i] for i in np.nonzero(other & both)[0][:20].tolist()],
        }
        msg = f"amplifies {pair.coverage['target_amplified']}/{n_t}"
        msg += " target sequences" if specific else " sequences"
        if specific:
            msg += f", {pair.coverage['other_amplified']}/{n_o} others"
        pair.notes.append(msg)
        if f.deg > 1 or r.deg > 1:
            pair.notes.append(f"degeneracy F ×{f.deg}, R ×{r.deg}")
        pair.extra["alignment_columns"] = [int(cols[f.start]) + 1, int(cols[r.end - 1]) + 1]
        pairs.append(finish_pair(pair, th))
        if len(pairs) >= s.num_return:
            break
    info = {"consensus": cons_seq, "degenerate_consensus": deg_seq, "rows": n_rows,
            "targets": int(target.sum()), "columns": int(P),
            "target_names": [names[i] for i in np.nonzero(target)[0].tolist()][:200]}
    return pairs, info


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


# --------------------------------------------------------------------------
# TE insertion genotyping
# --------------------------------------------------------------------------


def design_insertion(rid: str, seq: str, s: DesignSettings, te: tuple[int, int],
                     excluded: list[tuple[int, int]] | None = None,
                     misprime: dict | None = None) -> tuple[list[DesignedPair], list[str]]:
    """Flanking pair (empty site) + a primer inside the TE (occupied site)."""
    a, b = te  # 1-based inclusive
    n = len(seq)
    if not (1 <= a < b <= n):
        raise PrimerError("Set the TE start and end (1-based) inside the sequence.")
    te_len = b - a + 1
    left = seq[: a - 1]
    right = seq[b:]
    if len(left) < s.size_min + 20 or len(right) < s.size_min + 20:
        raise PrimerError("The TE needs at least ~40 bp of flank on both sides.")
    empty = left + right
    j = len(left)
    inc_lo = max(0, j - s.flank_max)
    inc_hi = min(len(empty), j + s.flank_max)
    ex = []
    for xa, xb in excluded or []:
        # keep only masks in the flanks, moved onto the excised template
        if xb <= a - 1:
            ex.append((xa, xb))
        elif xa >= b:
            ex.append((xa - te_len, xb - te_len))

    def to_genome(p: int, q: int, strand: str) -> tuple[int, int, str]:
        if q <= j:
            return p, q, strand
        return p + te_len, q + te_len, strand

    run = run_primer3(rid, empty, s, excluded=merge_intervals(ex), target=(j - 1, j + 1),
                      included=(inc_lo, inc_hi), misprime_lib=misprime, name_prefix="INS",
                      coord_map=to_genome, num_return=s.num_return)
    explain = list(run.explain)
    pairs = []
    for p in run.pairs:
        empty_size = p.product_size
        full_size = empty_size + te_len
        # inside primer: the forward primer fixed, a reverse primer in the first 300 bp of the TE
        f0 = p.forward.start - 1
        te_win = min(300, te_len)
        # The TE band must run clearly apart from the empty-site band.
        gap = max(40, int(empty_size * 0.2))
        ranges = [[max(60, s.product_min // 2), empty_size - gap],
                  [empty_size + gap, max(empty_size + gap, s.product_max)]]
        sub = run_primer3(rid, seq, s, left_primer=binding_part(p.forward.sequence),
                          ok_regions=[[f0, p.forward.length, a - 1, te_win]],
                          product=[r for r in ranges if r[1] >= r[0]],
                          num_return=1, name_prefix="TE")
        te_primer = None
        te_size = 0
        if sub.pairs:
            te_primer = sub.pairs[0].reverse
            te_primer.name = f"{p.name}_TE_R"
            te_primer.kind = "reverse"
            te_size = sub.pairs[0].product_size
        else:
            explain.extend(f"TE primer for {p.name}: {x}" for x in sub.explain)
        p.extra.update({"te_start": a, "te_end": b, "te_length": te_len,
                        "empty_product": empty_size, "occupied_flank_product": full_size,
                        "te_product": te_size})
        if te_primer is not None:
            p.extra["te_primer"] = te_primer.to_dict()
            dim = heterodimer(te_primer.sequence, p.forward.sequence, s.thermo)
            p.extra["te_primer_dimer_dg"] = dim.dg
        p.product_size = empty_size
        p.notes.append(f"empty site {empty_size} bp; with the TE {full_size:,} bp (F+R)"
                       + (f", {te_size} bp (F+TE)" if te_size else "; no TE primer found"))
        p.bands = [
            {"lane": "empty/empty", "sizes": [empty_size]},
            {"lane": "empty/TE", "sizes": [empty_size] + ([te_size] if te_size else [])
             + ([full_size] if full_size <= 5000 else [])},
            {"lane": "TE/TE", "sizes": ([te_size] if te_size else [])
             + ([full_size] if full_size <= 5000 else [])},
        ]
        pairs.append(p)
    return pairs, explain


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
