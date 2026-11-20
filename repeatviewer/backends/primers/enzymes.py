"""Restriction map of one sequence, with methylation from the user's own data.

Sites come from Biopython's copy of REBASE. Methylation sensitivity (CG, CHG,
Dam, Dcm) comes from published tables (Promega "Effect of site-specific
methylation on Promega restriction enzymes"; McClelland, Nelson & Raschke
1994, Nucleic Acids Res 22:3640; REBASE). With a Bismark ``.cov`` or a modkit
bedMethyl file, a site is drawn "probably not cut" when a cytosine of the site
(or of a context that overlaps it) is methylated in a context that blocks the
enzyme, at or above the chosen level and read depth.

Context of each cytosine is read on the sequence, on both strands: a C at
position i is CG when i+1 is G, CHG when i+1 is A/C/T and i+2 is G, else CHH;
a G at position i is the C of the other strand, read the same way leftwards.
Methylation files report the C of either strand at its own position (the G
position for the bottom strand); files merged to one CpG row (Bismark
``--merge_CpG``) are read through the partner position.

Dam and Dcm are labels only: they matter for DNA grown in E. coli, not plant
genomic DNA.

Pure Python, no Qt; unit-tested in ``tests/test_primers_enzymes.py``.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import os
import re
import threading
from dataclasses import dataclass, field
from pathlib import Path

# --------------------------------------------------------------------------
# Enzyme sets and sensitivity
# --------------------------------------------------------------------------

SENSITIVE, INSENSITIVE, OVERLAP = "s", "i", "ol"

# name: (CG, CHG, Dam, Dcm); s = blocked, ol = blocked when the methylated
# context overlaps the site, i = not blocked, None = not known.
# Promega table, MspI and HpaII from McClelland et al. 1994 (MspI cuts CmCGG,
# not mCCGG; HpaII is blocked by either).
SENSITIVITY: dict[str, tuple] = {
    "AatII": ("s", "i", "i", "i"),
    "ApaI": ("ol", "i", "i", "ol"),
    "AvaI": ("s", "i", "i", "i"),
    "BamHI": ("i", "ol", "i", "i"),
    "BclI": ("i", "i", "s", "i"),
    "BglI": ("ol", "ol", "i", "i"),
    "BglII": ("i", "ol", "i", "i"),
    "ClaI": ("s", "i", "ol", "i"),
    "DdeI": ("i", "ol", "i", "i"),
    "EcoRI": ("ol", "i", "i", "i"),
    "HaeIII": ("i", "ol", "i", "i"),
    "HhaI": ("s", "ol", "i", "i"),
    "HindIII": ("i", "i", "i", "i"),
    "HpaII": ("s", "s", "i", "i"),
    "KpnI": ("i", "i", "i", "i"),
    "MboI": ("i", "i", "s", "i"),
    "MluI": ("s", "i", "i", "i"),
    "MspI": ("i", "s", "i", "i"),
    "NheI": ("ol", "ol", "i", "i"),
    "NotI": ("s", "s", "i", "i"),
    "PstI": ("i", "s", "i", "i"),
    "PvuI": ("s", "ol", "i", "i"),
    "SacI": ("i", None, "i", "i"),
    "SacII": ("s", "s", "i", "i"),
    "SalI": ("s", "i", "i", "i"),
    "ScaI": ("i", "i", "i", "i"),
    "SmaI": ("s", "s", "i", "i"),
    "SphI": ("i", "i", "i", "i"),
    "TaqI": ("i", "i", "ol", "i"),
    "XbaI": ("i", "i", "ol", "i"),
    "XhoI": ("s", "i", "i", "i"),
}

SET_COMMON6 = "common6"
SET_FOUR = "four"
SET_METHYL = "methyl"
SET_ALL = "all"
SET_LABELS = {
    SET_COMMON6: "Common 6-cutters",
    SET_FOUR: "4- and 5-cutters",
    SET_METHYL: "Methylation-sensitive",
    SET_ALL: "All commercial",
}
SETS = (SET_COMMON6, SET_FOUR, SET_METHYL, SET_ALL)
_SET_NAMES = {
    SET_COMMON6: ("EcoRI", "BamHI", "HindIII", "KpnI", "SacI", "XhoI", "XbaI", "SalI",
                  "PstI", "SmaI", "SphI", "ScaI", "NheI", "BglII", "ClaI", "ApaI", "MluI",
                  "SacII", "EcoRV", "NcoI", "NdeI", "SpeI", "NotI"),
    SET_FOUR: ("MspI", "HpaII", "HaeIII", "HhaI", "MboI", "Sau3AI", "TaqI", "AluI", "RsaI",
               "MseI", "HinfI", "DdeI"),
    SET_METHYL: ("HpaII", "MspI", "PstI", "NotI", "SmaI", "SacII", "HhaI", "SalI", "XhoI",
                 "MluI", "AatII", "AvaI", "ClaI", "PvuI", "EcoRI", "BamHI", "HaeIII"),
}

SHOW_ONCE = "once"
SHOW_FEW = "few"
SHOW_MANY = "many"
SHOW_ANY = "any"
SHOW_NONE = "none"
SHOW_LABELS = {
    SHOW_ONCE: "Cut once",
    SHOW_FEW: "Cut 1 to 3 times",
    SHOW_MANY: "Cut 1 to 30 times",
    SHOW_ANY: "Cut any number",
    SHOW_NONE: "Do not cut",
}
SHOWS = (SHOW_ONCE, SHOW_FEW, SHOW_MANY, SHOW_ANY, SHOW_NONE)
_SHOW_RANGE = {SHOW_ONCE: (1, 1), SHOW_FEW: (1, 3), SHOW_MANY: (1, 30), SHOW_ANY: (1, None),
               SHOW_NONE: (0, 0)}

# Sets up to this size list every enzyme (out of range ones unticked, grey).
LIST_ALL_UP_TO = 60
MAX_DRAWN_SITES = 4000


def show_range(show: str) -> tuple[int, int | None]:
    return _SHOW_RANGE.get(show, _SHOW_RANGE[SHOW_MANY])


def in_range(cuts: int, show: str) -> bool:
    lo, hi = show_range(show)
    return cuts >= lo and (hi is None or cuts <= hi)


def _restriction():
    from Bio import Restriction

    return Restriction


def enzyme(name: str):
    enz = getattr(_restriction(), name, None)
    if enz is None or not hasattr(enz, "site"):
        raise ValueError(f"Unknown enzyme '{name}'.")
    return enz


def set_names(key: str) -> list[str]:
    """Enzyme names of a set, unknown names dropped, defined cutters only."""
    R = _restriction()
    if key == SET_ALL:
        names = sorted(str(e) for e in R.CommOnly)
    else:
        names = list(_SET_NAMES.get(key, _SET_NAMES[SET_COMMON6]))
    out = []
    for n in names:
        enz = getattr(R, n, None)
        if enz is not None and enz.is_defined():
            out.append(n)
    return out


def blocking_contexts(name: str) -> frozenset[str]:
    s = SENSITIVITY.get(name)
    if s is None:
        return frozenset()
    return frozenset(c for c, v in zip(("CG", "CHG"), s[:2]) if v in (SENSITIVE, OVERLAP))


def sensitivity_text(name: str) -> str:
    """One line: what blocks the enzyme ("not known" when no table has it)."""
    s = SENSITIVITY.get(name)
    if s is None:
        return "Methylation sensitivity not recorded here; drawn as cut."
    words = {SENSITIVE: "blocked", OVERLAP: "blocked when overlapping", INSENSITIVE: "not blocked",
             None: "not known"}
    parts = [f"m{c}: {words[v]}" for c, v in zip(("CG", "CHG"), s[:2])]
    host = [f"{c}: {words[v]}" for c, v in zip(("Dam", "Dcm"), s[2:])
            if v in (SENSITIVE, OVERLAP)]
    text = "; ".join(parts)
    if host:
        text += " (E. coli DNA only: " + "; ".join(host) + ")"
    return text


def is_defined(name: str) -> bool:
    try:
        return bool(enzyme(name).is_defined())
    except ValueError:
        return False


# --------------------------------------------------------------------------
# Methylation data
# --------------------------------------------------------------------------


@dataclass
class MethData:
    """Methylation of one sequence: 1-based position → (percent, reads)."""

    path: str = ""
    chrom: str = ""
    offset: int = 0            # sequence position 1 = chrom position offset + 1
    calls: dict[int, tuple[float, int]] = field(default_factory=dict)
    note: str = ""

    def at(self, pos: int) -> tuple[float, int] | None:
        return self.calls.get(pos)


_REGION_ID = re.compile(r"^(?P<chrom>.+?)[:_](?P<start>\d[\d,]*)[-_.]+(?P<end>\d[\d,]*)$")


def chrom_candidates(seq_id: str) -> list[tuple[str, int]]:
    """(chromosome name, offset) to look for in a methylation file.

    The id itself, its first word, and ``chrom:start-end`` for an extracted
    region (positions then shift by start - 1).
    """
    out: list[tuple[str, int]] = []
    sid = (seq_id or "").strip()
    first = sid.split()[0] if sid else ""
    for name in (sid, first):
        if name and (name, 0) not in out:
            out.append((name, 0))
    m = _REGION_ID.match(first)
    if m:
        start = int(m.group("start").replace(",", ""))
        out.append((m.group("chrom"), max(0, start - 1)))
    return out


_CACHE: dict[tuple, MethData] = {}
_CACHE_LOCK = threading.Lock()


def methylation_format(path: Path) -> str:
    from repeatviewer.backends.tree_traits import methylation_format as fmt

    return fmt(Path(path))


def load_methylation(path: str | Path, seq_id: str, length: int, *,
                     chunk: int = 2_000_000) -> MethData:
    """5mC calls of the sequence from a Bismark .cov or a modkit bedMethyl file.

    Only rows of the matching chromosome inside the sequence are kept. Cached
    by file, modification time and sequence.
    """
    import pandas as pd

    path = Path(path)
    try:
        mtime = os.path.getmtime(path)
    except OSError as exc:
        return MethData(path=str(path), note=f"Cannot read {path.name}: {exc}")
    key = (str(path), mtime, seq_id, int(length))
    with _CACHE_LOCK:
        hit = _CACHE.get(key)
    if hit is not None:
        return hit
    try:
        fmt = methylation_format(path)
    except (OSError, ValueError) as exc:
        return MethData(path=str(path), note=str(exc))
    cands = chrom_candidates(seq_id)
    wanted = {c for c, _o in cands}
    if fmt == "bismark":
        cols, names = [0, 1, 4, 5], ["chrom", "pos", "m", "u"]
    else:
        cols, names = [0, 1, 3, 11, 12], ["chrom", "pos", "code", "m", "u"]
    found: dict[str, list] = {}
    reader = pd.read_csv(path, sep=r"\s+", header=None, comment="#", usecols=cols, names=names,
                         chunksize=chunk, compression="infer", dtype={"chrom": str},
                         on_bad_lines="skip", engine="c")
    for df in reader:
        df = df[df["chrom"].isin(wanted)]
        if "code" in df:
            df = df[df["code"].astype(str) == "m"]
        if not df.empty:
            for chrom, sub in df.groupby("chrom", sort=False):
                found.setdefault(chrom, []).append(sub)
    data = MethData(path=str(path))
    for chrom, offset in cands:
        if chrom not in found:
            continue
        sub = pd.concat(found[chrom])
        pos = sub["pos"].to_numpy().astype("int64")
        if fmt == "bedmethyl":
            pos = pos + 1          # 0-based start
        pos = pos - offset
        m = sub["m"].to_numpy().astype(float)
        u = sub["u"].to_numpy().astype(float)
        keep = (pos >= 1) & (pos <= length)
        if not keep.any():
            continue
        calls = {}
        for p, a, b in zip(pos[keep], m[keep], u[keep]):
            depth = int(a + b)
            if depth > 0:
                calls[int(p)] = (float(100.0 * a / depth), depth)
        data.chrom, data.offset, data.calls = chrom, offset, calls
        break
    if not data.calls:
        tried = ", ".join(sorted(wanted)) or "(no name)"
        data.note = (f"No methylation rows for {tried} in {path.name}. The sequence name has "
                     "to match the chromosome name in the file.")
    with _CACHE_LOCK:
        if len(_CACHE) > 16:
            _CACHE.clear()
        _CACHE[key] = data
    return data


# --------------------------------------------------------------------------
# Sites
# --------------------------------------------------------------------------

_COMP = {"A": "T", "C": "G", "G": "C", "T": "A"}


def c_context(seq: str, i: int) -> tuple[str, str] | None:
    """(strand, context) of a cytosine at 0-based ``i`` (top C or bottom C under a G)."""
    n = len(seq)
    b = seq[i]
    if b == "C":
        n1 = seq[i + 1] if i + 1 < n else ""
        n2 = seq[i + 2] if i + 2 < n else ""
        if n1 == "G":
            return "+", "CG"
        if n1 in "ACT" and n1 and n2 == "G":
            return "+", "CHG"
        return ("+", "CHH") if n1 and n2 else None
    if b == "G":
        n1 = _COMP.get(seq[i - 1], "") if i >= 1 else ""
        n2 = _COMP.get(seq[i - 2], "") if i >= 2 else ""
        if n1 == "G":
            return "-", "CG"
        if n1 in "ACT" and n1 and n2 == "G":
            return "-", "CHG"
        return ("-", "CHH") if n1 and n2 else None
    return None


@dataclass
class MethCall:
    pos: int          # 1-based
    strand: str
    context: str
    percent: float
    reads: int
    blocks: bool = False      # the context blocks this enzyme (any level)


@dataclass
class Site:
    enzyme: str
    start: int        # 1-based, recognition site
    end: int
    cut: int          # 1-based base right of the top-strand cut
    blocked: bool = False
    calls: list[MethCall] = field(default_factory=list)
    reason: str = ""

    def tooltip(self, site_text: str) -> str:
        head = f"<b>{self.enzyme}</b> ({site_text}) at {self.start:,}–{self.end:,}, cuts after {self.cut - 1:,}"
        lines = [head]
        for c in self.calls:
            lines.append(f"m{c.context} {c.strand} at {c.pos:,}: {c.percent:.0f} % of {c.reads} reads")
        lines.append(self.reason or ("Cut." if not self.blocked else "Probably not cut."))
        return "<br>".join(lines)


@dataclass
class EnzymeRow:
    name: str
    site: str
    sites: list[Site]

    @property
    def cuts(self) -> int:
        return len(self.sites)

    @property
    def blocked(self) -> int:
        return sum(1 for s in self.sites if s.blocked)


@dataclass
class MapSettings:
    enzyme_set: str = SET_COMMON6
    show: str = SHOW_MANY
    meth_path: str = ""
    min_percent: float = 50.0
    min_reads: int = 5

    def to_dict(self) -> dict:
        return {"enzyme_set": self.enzyme_set, "show": self.show, "meth_path": self.meth_path,
                "min_percent": self.min_percent, "min_reads": self.min_reads}

    @classmethod
    def from_dict(cls, d: dict | None) -> MapSettings:
        d = d or {}
        s = cls()
        if d.get("enzyme_set") in SETS:
            s.enzyme_set = d["enzyme_set"]
        if d.get("show") in SHOWS:
            s.show = d["show"]
        s.meth_path = str(d.get("meth_path") or "")
        try:
            s.min_percent = float(d.get("min_percent", s.min_percent))
            s.min_reads = int(d.get("min_reads", s.min_reads))
        except (TypeError, ValueError):
            pass
        return s


def _cut_positions(enz, seq) -> list[int]:
    try:
        return sorted(int(p) for p in enz.search(seq, linear=True))
    except Exception:  # noqa: BLE001 - odd enzymes (no defined cut)
        return []


def find_sites(name: str, sequence: str) -> list[Site]:
    """Recognition sites (both strands) and their top-strand cut, 1-based."""
    from Bio.Seq import Seq

    enz = enzyme(name)
    seq = sequence.upper()
    size = int(enz.size)
    cuts = _cut_positions(enz, Seq(seq))
    if not cuts:
        return []
    out = []
    seen = set()
    for m in enz.compsite.finditer(seq):
        start0 = m.start()
        if start0 in seen or start0 + size > len(seq):
            continue
        seen.add(start0)
        centre = start0 + size / 2.0
        cut = min(cuts, key=lambda c: abs(c - 1 - centre))
        out.append(Site(enzyme=name, start=start0 + 1, end=start0 + size, cut=cut))
    return out


def judge_site(site: Site, sequence: str, meth: MethData | None, *, min_percent: float,
               min_reads: int) -> None:
    """Fill ``site.calls`` / ``blocked`` / ``reason`` from the data and the enzyme's sensitivity."""
    seq = sequence.upper()
    sens = dict(zip(("CG", "CHG"), (SENSITIVITY.get(site.enzyme) or (None, None))[:2]))
    if meth is None or not meth.calls:
        site.reason = sensitivity_text(site.enzyme)
        return
    blocked_by = []
    used: set[int] = set()
    # Cytosines of the site, plus the overlapping context one base either side
    # (an overlapping CG needs the flanking base).
    lo = max(0, site.start - 2)
    hi = min(len(seq), site.end + 1)
    for i in range(lo, hi):
        ctx = c_context(seq, i)
        if ctx is None:
            continue
        strand, context = ctx
        inside = site.start - 1 <= i <= site.end - 1
        row = i + 1
        call = meth.at(row)
        if call is None and context == "CG":
            # merged CpG files keep one row for both Cs of the CG
            row = i if strand == "-" else i + 2
            call = meth.at(row)
        if call is None or row in used:
            continue
        used.add(row)
        pct, reads = call
        if not inside and context != "CG" and context != "CHG":
            continue
        rule = sens.get(context)
        counts = rule == OVERLAP or (rule == SENSITIVE and inside)
        site.calls.append(MethCall(pos=i + 1, strand=strand, context=context, percent=pct,
                                   reads=reads, blocks=counts))
        if counts and pct >= min_percent and reads >= min_reads:
            blocked_by.append(f"m{context} at {i + 1:,} ({pct:.0f} %, {reads} reads)")
    if blocked_by:
        site.blocked = True
        site.reason = (f"Probably not cut: {site.enzyme} is blocked by "
                       + ", ".join(sorted({b.split()[0] for b in blocked_by})) + "; "
                       + "; ".join(blocked_by) + ".")
    elif site.calls:
        site.reason = "Cut: " + sensitivity_text(site.enzyme)
    else:
        site.reason = "No reads at this site in your data. " + sensitivity_text(site.enzyme)


def blocked_fraction(site: Site, min_reads: int) -> float | None:
    """Share of molecules methylated in a blocking context (highest C), None without reads."""
    vals = [c.percent / 100.0 for c in site.calls if c.blocks and c.reads >= min_reads]
    return max(vals) if vals else None


def restriction_map(sequence: str, settings: MapSettings, meth: MethData | None = None,
                    names: list[str] | None = None) -> list[EnzymeRow]:
    """One row per enzyme of the set (or ``names``), sites judged against the data."""
    rows = []
    for name in names if names is not None else set_names(settings.enzyme_set):
        try:
            enz = enzyme(name)
            sites = find_sites(name, sequence)
        except ValueError:
            continue
        for s in sites:
            judge_site(s, sequence, meth, min_percent=settings.min_percent,
                       min_reads=settings.min_reads)
        rows.append(EnzymeRow(name=name, site=str(enz.site), sites=sites))
    return rows


def listed_rows(rows: list[EnzymeRow], show: str) -> list[EnzymeRow]:
    """Rows for the table: every enzyme of a small set, the in-range ones of a big set."""
    if len(rows) <= LIST_ALL_UP_TO:
        inside = [r for r in rows if in_range(r.cuts, show)]
        outside = [r for r in rows if not in_range(r.cuts, show)]
        return sorted(inside, key=lambda r: (r.cuts, r.name)) + \
            sorted(outside, key=lambda r: (r.cuts == 0, r.cuts, r.name))
    return sorted((r for r in rows if in_range(r.cuts, show)), key=lambda r: (r.cuts, r.name))


def fragments(length: int, cuts: list[int]) -> list[int]:
    """Fragment sizes of a linear sequence cut at 1-based positions (base right of the cut)."""
    edges = sorted({c for c in cuts if 1 < c <= length})
    out, prev = [], 1
    for c in edges:
        out.append(c - prev)
        prev = c
    out.append(length - prev + 1)
    return [f for f in out if f > 0]


METH_SUFFIXES = (".cov", ".cov.gz", ".bed", ".bed.gz", ".bedmethyl", ".bedmethyl.gz")


def methylation_files_near(path: str | Path | None, limit: int = 30) -> list[Path]:
    """Methylation files in the folder of the open file (newest first)."""
    if not path:
        return []
    folder = Path(path).parent
    try:
        files = [p for p in folder.iterdir()
                 if p.is_file() and p.name.lower().endswith(METH_SUFFIXES)]
    except OSError:
        return []
    files.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    return files[:limit]
