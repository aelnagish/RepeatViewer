"""Primer thermodynamics: Tm, GC, hairpins, dimers and 3' end stability.

Thin layer over primer3-py (the same Primer3 engine Geneious uses), with the
settings shown in the Tm Calculation section of the primer windows. Every
function takes plain strings, so the rest of the primer code never touches
primer3 directly and can be tested without Qt.

A primer may carry a lower-case 5' tail (adapter, restriction site, M13…).
The tail does not bind the template, so Tm and binding use only the upper-case
3' part (see :func:`split_tail`); dimers and hairpins use the whole oligo.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import itertools
import math
import re
from dataclasses import asdict, dataclass, field

# --------------------------------------------------------------------------
# Alphabet
# --------------------------------------------------------------------------

IUPAC_BASES: dict[str, str] = {
    "A": "A", "C": "C", "G": "G", "T": "T",
    "R": "AG", "Y": "CT", "S": "CG", "W": "AT", "K": "GT", "M": "AC",
    "B": "CGT", "D": "AGT", "H": "ACT", "V": "ACG", "N": "ACGT",
}
_COMPLEMENT = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn-",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn-",
)
_VALID = set(IUPAC_BASES) | {c.lower() for c in IUPAC_BASES}
# Bases written back from a set of bases (inverse of IUPAC_BASES).
_CODE_FOR = {frozenset(v): k for k, v in IUPAC_BASES.items()}


def reverse_complement(seq: str) -> str:
    return str(seq).translate(_COMPLEMENT)[::-1]


def iupac_code(bases) -> str:
    """IUPAC letter for a set of bases (e.g. {'A', 'G'} → 'R')."""
    key = frozenset(b for b in bases if b in "ACGT")
    return _CODE_FOR.get(key, "N") if key else "N"


class PrimerError(ValueError):
    """A primer job cannot go on (bad sequence, missing primer3, bad settings)."""


PRIMER3_MISSING = (
    "The primer tools need primer3-py.\n\n"
    "Install it in the RepeatViewer environment and restart:\n"
    "    pip install primer3-py\n"
    "or\n"
    "    conda install -c bioconda primer3-py"
)


def primer3_available() -> bool:
    try:
        import primer3  # noqa: F401
    except ImportError:
        return False
    return True


def require_primer3():
    try:
        import primer3
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise PrimerError(PRIMER3_MISSING) from exc
    return primer3


def primer3_version() -> str:
    try:
        import primer3

        return str(getattr(primer3, "__version__", ""))
    except ImportError:
        return ""


def clean_oligo(text: str) -> str:
    """Oligo as typed → IUPAC letters only, case kept (lower-case = 5' tail).

    Spaces, digits, dashes and 5'/3' marks are dropped; U becomes T. Raises
    :class:`PrimerError` for any other letter.
    """
    raw = str(text or "").strip()
    raw = re.sub(r"^5['’′]?-?|-?3['’′]?$", "", raw)
    raw = raw.replace("U", "T").replace("u", "t")
    out = []
    for ch in raw:
        if ch in _VALID:
            out.append(ch)
        elif ch.isspace() or ch.isdigit() or ch in "-_.,;:'’′":
            continue
        else:
            raise PrimerError(f"'{ch}' is not a DNA base (use A C G T or IUPAC codes).")
    return "".join(out)


def split_tail(seq: str) -> tuple[str, str]:
    """Return (5' tail, binding part).

    A lower-case 5' prefix followed by upper case is a tail. An oligo typed all
    in lower case (or all upper case) has no tail: all of it binds.
    """
    s = str(seq or "")
    m = re.match(r"^([a-z]+)([A-Z].*)$", s)
    if not m:
        return "", s.upper()
    return m.group(1), m.group(2).upper()


def binding_part(seq: str) -> str:
    return split_tail(seq)[1]


def degeneracy(seq: str) -> int:
    n = 1
    for ch in str(seq).upper():
        n *= len(IUPAC_BASES.get(ch, "N"))
    return n


def expand_iupac(seq: str, limit: int = 256) -> list[str]:
    """All plain-ACGT versions of a degenerate oligo (at most ``limit``, evenly picked)."""
    choices = [IUPAC_BASES.get(ch, "ACGT") for ch in str(seq).upper()]
    total = 1
    for c in choices:
        total *= len(c)
    if total <= limit:
        return ["".join(p) for p in itertools.product(*choices)]
    step = total / limit
    out = []
    for k in range(limit):
        idx = int(k * step)
        letters = []
        for c in reversed(choices):
            idx, r = divmod(idx, len(c))
            letters.append(c[r])
        out.append("".join(reversed(letters)))
    return out


def gc_percent(seq: str) -> float:
    """GC % of an oligo; degenerate bases count by their share of G/C."""
    s = str(seq).upper()
    if not s:
        return 0.0
    total = 0.0
    for ch in s:
        bases = IUPAC_BASES.get(ch, "ACGT")
        total += sum(1 for b in bases if b in "GC") / len(bases)
    return 100.0 * total / len(s)


def longest_run(seq: str) -> tuple[str, int]:
    """Longest mononucleotide run (base, length)."""
    best = ("", 0)
    for base, grp in itertools.groupby(str(seq).upper()):
        n = sum(1 for _ in grp)
        if n > best[1]:
            best = (base, n)
    return best


def longest_dinucleotide_repeat(seq: str) -> int:
    """Copies of the longest dinucleotide repeat (ATATAT → 3)."""
    s = str(seq).upper()
    best = 1
    for start in range(min(2, len(s))):
        i = start
        while i + 1 < len(s):
            unit = s[i:i + 2]
            if unit[0] == unit[1]:
                i += 1
                continue
            n = 1
            j = i + 2
            while s[j:j + 2] == unit:
                n += 1
                j += 2
            best = max(best, n)
            i = j if n > 1 else i + 1
    return best


def gc_clamp_count(seq: str, last: int = 5) -> int:
    return sum(1 for ch in str(seq).upper()[-last:] if ch in "GC")


# SantaLucia (1998) unified nearest-neighbour ΔG°37 (kcal/mol), for the
# 3' end stability (ΔG of the last five bases, as Primer3 reports it).
_NN_DG = {
    "AA": -1.00, "TT": -1.00, "AT": -0.88, "TA": -0.58,
    "CA": -1.45, "TG": -1.45, "GT": -1.44, "AC": -1.44,
    "CT": -1.28, "AG": -1.28, "GA": -1.30, "TC": -1.30,
    "CG": -2.17, "GC": -2.24, "GG": -1.84, "CC": -1.84,
}


def three_prime_dg(seq: str, last: int = 5) -> float:
    """ΔG (kcal/mol) of the duplex formed by the last ``last`` bases.

    More negative = a more stable 3' end, which primes more readily on partial
    matches (Primer3 rejects < −9 by default). Degenerate bases use the mean.
    """
    s = binding_part(seq)[-last:]
    if len(s) < 2:
        return 0.0
    variants = expand_iupac(s, 32)
    values = []
    for v in variants:
        dg = sum(_NN_DG.get(v[i:i + 2], -1.3) for i in range(len(v) - 1))
        for end in (v[0], v[-1]):
            dg += 0.98 if end in "GC" else 1.03
        values.append(dg)
    return round(sum(values) / len(values), 2)


# --------------------------------------------------------------------------
# Settings
# --------------------------------------------------------------------------

TM_METHODS = ("santalucia", "breslauer")
SALT_METHODS = ("santalucia", "owczarzy", "schildkraut")
TM_METHOD_LABELS = {"santalucia": "SantaLucia 1998", "breslauer": "Breslauer 1986"}
SALT_METHOD_LABELS = {
    "santalucia": "SantaLucia 1998",
    "owczarzy": "Owczarzy 2004",
    "schildkraut": "Schildkraut 1965",
}


@dataclass
class ThermoSettings:
    """Reaction conditions for Tm / ΔG (Geneious' Tm Calculation section)."""

    tm_method: str = "santalucia"
    salt_method: str = "santalucia"
    monovalent_mM: float = 50.0     # Na+ / K+
    divalent_mM: float = 1.5        # Mg2+
    dntp_mM: float = 0.6
    oligo_nM: float = 50.0
    temperature_c: float = 37.0     # temperature for hairpin / dimer ΔG

    def __post_init__(self) -> None:
        if self.tm_method not in TM_METHODS:
            self.tm_method = "santalucia"
        if self.salt_method not in SALT_METHODS:
            self.salt_method = "santalucia"
        self.monovalent_mM = max(0.0, float(self.monovalent_mM))
        self.divalent_mM = max(0.0, float(self.divalent_mM))
        self.dntp_mM = max(0.0, float(self.dntp_mM))
        self.oligo_nM = max(0.01, float(self.oligo_nM))
        self.temperature_c = float(self.temperature_c)

    def conc_kwargs(self) -> dict:
        return {
            "mv_conc": self.monovalent_mM,
            "dv_conc": self.divalent_mM,
            "dntp_conc": self.dntp_mM,
            "dna_conc": self.oligo_nM,
        }

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict | None) -> ThermoSettings:
        data = dict(data or {})
        known = {k: data[k] for k in cls.__dataclass_fields__ if k in data}
        return cls(**known)

    def describe(self) -> str:
        return (
            f"{TM_METHOD_LABELS[self.tm_method]} Tm, {SALT_METHOD_LABELS[self.salt_method]} "
            f"salt correction; {self.monovalent_mM:g} mM monovalent, "
            f"{self.divalent_mM:g} mM Mg²⁺, {self.dntp_mM:g} mM dNTP, "
            f"{self.oligo_nM:g} nM oligo"
        )


# primer3 limits the thermodynamic oligo length (THAL_MAX_ALIGN = 60).
_THAL_MAX = 60


def _thal_seq(seq: str) -> str:
    """primer3 thal accepts ACGT only (≤ 60 nt): degenerate bases → first choice."""
    s = "".join(IUPAC_BASES.get(ch, "A")[0] for ch in str(seq).upper())
    return s[-_THAL_MAX:]


def calc_tm(seq: str, settings: ThermoSettings | None = None) -> float:
    """Melting temperature (°C) of the binding part; mean over degenerate versions."""
    s = binding_part(seq)
    if not s:
        return float("nan")
    p3 = require_primer3()
    st = settings or ThermoSettings()
    values = []
    for v in expand_iupac(s, 64):
        values.append(p3.bindings.calc_tm(
            v, **st.conc_kwargs(), tm_method=st.tm_method,
            salt_corrections_method=st.salt_method, max_nn_length=60,
        ))
    return round(sum(values) / len(values), 2)


def calc_tm_range(seq: str, settings: ThermoSettings | None = None) -> tuple[float, float]:
    s = binding_part(seq)
    p3 = require_primer3()
    st = settings or ThermoSettings()
    values = [
        p3.bindings.calc_tm(v, **st.conc_kwargs(), tm_method=st.tm_method,
                            salt_corrections_method=st.salt_method, max_nn_length=60)
        for v in expand_iupac(s, 64)
    ]
    return round(min(values), 2), round(max(values), 2)


@dataclass
class Structure:
    """A hairpin / dimer from Primer3's thal: Tm, ΔG and a text drawing."""

    kind: str
    tm: float = 0.0
    dg: float = 0.0          # kcal/mol
    lines: list[str] = field(default_factory=list)

    @property
    def found(self) -> bool:
        return bool(self.lines) or self.dg < 0

    def to_dict(self) -> dict:
        return asdict(self)


def _structure(kind: str, result, *, label_a: str = "5'", label_b: str = "3'") -> Structure:
    found = bool(getattr(result, "structure_found", False))
    if not found:
        return Structure(kind=kind)
    lines = list(getattr(result, "ascii_structure_lines", None) or [])
    dg = round(float(result.dg) / 1000.0, 2)
    tm = round(float(result.tm), 1)
    tidy = _tidy_hairpin(lines) if kind == "hairpin" else _tidy_structure(lines)
    return Structure(kind=kind, tm=tm, dg=dg, lines=tidy)


def _tidy_hairpin(lines: list[str]) -> list[str]:
    """primer3's hairpin lines (``--///---\\\\\\---`` over the sequence) → dot-bracket::

        5′ AAACGAACCGTGGGGATTTTC 3′
           ..(((...)))..........
    """
    bodies = [ln.partition("\t")[2] for ln in lines]
    if len(bodies) < 2:
        return [b for b in bodies if b.strip()]
    marks, seq = bodies[0], bodies[1]
    table = str.maketrans({"/": "(", "\\": ")", "-": "."})
    return ["5′ " + seq + " 3′", "   " + marks.translate(table)]


def _tidy_structure(lines: list[str]) -> list[str]:
    """primer3's four SEQ/STR lines → a duplex drawing with ends and pairs marked::

        5′ ACACGGAGGAG   A   GA          3′
                      GCT AGC
                      ||| |||
                      CGA TCG
        3′ AG   A   GAGGAGGCACA          5′
    """
    bodies = [ln.partition("\t")[2] for ln in lines]
    if len(bodies) < 4:
        return [b for b in bodies if b.strip()]
    top_free, top_pair, bot_pair, bot_free = bodies[:4]
    width = max(len(b) for b in bodies[:4])
    bars = "".join("|" if i < len(top_pair) and i < len(bot_pair) and top_pair[i].strip()
                   and bot_pair[i].strip() else " " for i in range(width))
    out = [
        "5′ " + top_free.ljust(width) + " 3′",
        "   " + top_pair.ljust(width),
        "   " + bars,
        "   " + bot_pair.ljust(width),
        "3′ " + bot_free.ljust(width) + " 5′",
    ]
    return [x.rstrip() for x in out]


def hairpin(seq: str, settings: ThermoSettings | None = None) -> Structure:
    p3 = require_primer3()
    st = settings or ThermoSettings()
    s = _thal_seq(seq)
    if len(s) < 4:
        return Structure(kind="hairpin")
    r = p3.bindings.calc_hairpin(s, **st.conc_kwargs(), temp_c=st.temperature_c,
                                 output_structure=True)
    return _structure("hairpin", r)


def homodimer(seq: str, settings: ThermoSettings | None = None) -> Structure:
    p3 = require_primer3()
    st = settings or ThermoSettings()
    s = _thal_seq(seq)
    if len(s) < 2:
        return Structure(kind="self-dimer")
    r = p3.bindings.calc_homodimer(s, **st.conc_kwargs(), temp_c=st.temperature_c,
                                   output_structure=True)
    return _structure("self-dimer", r)


def heterodimer(a: str, b: str, settings: ThermoSettings | None = None) -> Structure:
    p3 = require_primer3()
    st = settings or ThermoSettings()
    sa, sb = _thal_seq(a), _thal_seq(b)
    if len(sa) < 2 or len(sb) < 2:
        return Structure(kind="cross-dimer")
    r = p3.bindings.calc_heterodimer(sa, sb, **st.conc_kwargs(), temp_c=st.temperature_c,
                                     output_structure=True)
    return _structure("cross-dimer", r)


def duplex_tm(primer: str, site: str, settings: ThermoSettings | None = None) -> float:
    """Tm of a primer on its (possibly mismatched) binding site.

    ``site`` is the template stretch the primer copies (same orientation as
    the primer), so the strand it anneals to is its reverse complement.
    """
    p3 = require_primer3()
    st = settings or ThermoSettings()
    a = _thal_seq(binding_part(primer))
    b = _thal_seq(reverse_complement(site.upper()))
    if len(a) < 2 or len(b) < 2:
        return float("nan")
    r = p3.bindings.calc_heterodimer(a, b, **st.conc_kwargs(), temp_c=st.temperature_c)
    if not r.structure_found:
        return 0.0
    return round(float(r.tm), 1)


# --------------------------------------------------------------------------
# One oligo, analysed
# --------------------------------------------------------------------------

# Warning thresholds (the usual primer-design rules of thumb).
WARN_GC = (40.0, 60.0)
WARN_RUN = 5
WARN_DINUC = 4
WARN_DIMER_DG = -9.0       # kcal/mol
WARN_HAIRPIN_DG = -3.0
WARN_END_DG = -9.0
WARN_CLAMP_MAX = 3


@dataclass
class OligoStats:
    sequence: str
    length: int
    binding_length: int
    tm: float
    tm_min: float
    tm_max: float
    gc: float
    degeneracy: int
    end_dg: float
    hairpin: Structure
    self_dimer: Structure
    warnings: list[str] = field(default_factory=list)

    def to_dict(self) -> dict:
        d = asdict(self)
        d["hairpin"] = self.hairpin.to_dict()
        d["self_dimer"] = self.self_dimer.to_dict()
        return d


def analyze_oligo(seq: str, settings: ThermoSettings | None = None) -> OligoStats:
    st = settings or ThermoSettings()
    full = clean_oligo(seq)
    tail, bind = split_tail(full)
    if not bind:
        raise PrimerError("Empty oligo.")
    tm_lo, tm_hi = calc_tm_range(bind, st)
    tm = calc_tm(bind, st)
    gc = gc_percent(bind)
    hp = hairpin(full, st)
    sd = homodimer(full, st)
    end_dg = three_prime_dg(bind)
    warnings: list[str] = []
    if not (WARN_GC[0] <= gc <= WARN_GC[1]):
        warnings.append(f"GC {gc:.0f} % (aim for 40–60 %)")
    base, run = longest_run(bind)
    if run >= WARN_RUN:
        warnings.append(f"run of {run} {base}")
    di = longest_dinucleotide_repeat(bind)
    if di >= WARN_DINUC:
        warnings.append(f"dinucleotide repeat ×{di}")
    if gc_clamp_count(bind) > WARN_CLAMP_MAX:
        warnings.append(f"{gc_clamp_count(bind)} G/C in the last 5 bases")
    if gc_clamp_count(bind, 2) == 0:
        warnings.append("no G/C at the 3' end")
    if hp.dg < WARN_HAIRPIN_DG and hp.tm > st.temperature_c:
        warnings.append(f"hairpin ΔG {hp.dg:.1f}")
    if sd.dg < WARN_DIMER_DG:
        warnings.append(f"self-dimer ΔG {sd.dg:.1f}")
    if end_dg < WARN_END_DG:
        warnings.append(f"very stable 3' end (ΔG {end_dg:.1f})")
    deg = degeneracy(bind)
    if deg > 1:
        warnings.append(f"degenerate ×{deg}")
    return OligoStats(
        sequence=full, length=len(full), binding_length=len(bind), tm=tm,
        tm_min=tm_lo, tm_max=tm_hi, gc=round(gc, 1), degeneracy=deg,
        end_dg=end_dg, hairpin=hp, self_dimer=sd, warnings=warnings,
    )


def pair_dimer_warning(dimer: Structure) -> str:
    return f"primer-dimer ΔG {dimer.dg:.1f}" if dimer.dg < WARN_DIMER_DG else ""


def molecular_weight(seq: str) -> float:
    """Approximate MW (g/mol) of a single-stranded DNA oligo (anhydrous, 5'-OH)."""
    w = {"A": 313.21, "C": 289.18, "G": 329.21, "T": 304.2}
    s = binding_part(seq) if seq else ""
    full = str(seq).upper()
    total = 0.0
    for ch in full:
        bases = IUPAC_BASES.get(ch, "ACGT")
        total += sum(w[b] for b in bases) / len(bases)
    del s
    return round(total - 61.96, 1) if full else 0.0


def extinction_260(seq: str) -> float:
    """ε260 (L/(mol·cm)) by nearest neighbours, for ordering / OD to µM."""
    eps_nn = {
        "AA": 27400, "AC": 21200, "AG": 25000, "AT": 22800, "CA": 21200,
        "CC": 14600, "CG": 18000, "CT": 15200, "GA": 25200, "GC": 17600,
        "GG": 21600, "GT": 20000, "TA": 23400, "TC": 16200, "TG": 19000,
        "TT": 16800,
    }
    eps_n = {"A": 15400, "C": 7400, "G": 11500, "T": 8700}
    s = "".join(IUPAC_BASES.get(ch, "A")[0] for ch in str(seq).upper())
    if not s:
        return 0.0
    if len(s) == 1:
        return float(eps_n[s])
    total = sum(eps_nn[s[i:i + 2]] for i in range(len(s) - 1))
    total -= sum(eps_n[s[i]] for i in range(1, len(s) - 1))
    return float(total)


def fmt_tm(value: float) -> str:
    return "–" if value is None or (isinstance(value, float) and math.isnan(value)) else f"{value:.1f}"
