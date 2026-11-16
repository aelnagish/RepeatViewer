"""In-silico PCR (Tools > Primers and Probes > PCR). Qt-free.

Kinds of PCR:

* **Standard**: every product the primers can make on each template (all
  primer combinations, as in a real tube), with the mismatches per primer.
* **Satellite ladder**: the same products grouped by size; when the primers
  sit in a tandem repeat, sizes step by the monomer and are labelled 1×, 2×, …

Each product carries the primers' 5' tails (copied in from cycle 3 on) and a
relative yield used by the gel: fewer mismatches, an exact 3' end and a
shorter product make more DNA; identical products from repeated sites add up.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import datetime as _dt
import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Callable

from repeatviewer.backends.primers import binding as bd
from repeatviewer.backends.primers import thermo as th
from repeatviewer.backends.primers.oligos import Primer

ProgressFn = Callable[[str], None] | None
CancelFn = Callable[[], bool] | None

MODE_STANDARD = "standard"
MODE_LADDER = "ladder"
MODES = (MODE_STANDARD, MODE_LADDER)
MODE_LABELS = {
    MODE_STANDARD: "Standard",
    MODE_LADDER: "Satellite ladder",
}

POLYMERASES = {
    # name: (bases per second, length where yield has dropped to 1/e)
    "Taq": (1000 / 60, 3000),
    "Phusion / Q5": (1000 / 20, 10000),
    "Long-range": (1000 / 60, 20000),
}

MAX_PRODUCTS = 5000
MAX_WRITTEN = 500
PCR_TAG = "_PCR"
RESULT_FORMAT = "repeatviewer-pcr-1"


@dataclass
class PcrSettings:
    mode: str = MODE_STANDARD
    max_mismatches: int = 2
    three_prime_exact: int = 3
    min_product: int = 30
    max_product: int = 5000
    polymerase: str = "Taq"
    annealing_c: float = 58.0
    cycles: int = 30
    dimers: bool = True

    def __post_init__(self) -> None:
        if self.mode not in MODES:   # also the removed "bisulfite" and "presence" modes
            self.mode = MODE_STANDARD
        if self.polymerase not in POLYMERASES:
            self.polymerase = "Taq"
        self.max_mismatches = max(0, int(self.max_mismatches))
        self.three_prime_exact = max(0, int(self.three_prime_exact))
        self.min_product = max(1, int(self.min_product))
        self.max_product = max(self.min_product, int(self.max_product))

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict | None) -> PcrSettings:
        data = dict(data or {})
        return cls(**{k: data[k] for k in cls.__dataclass_fields__ if k in data})

    def search(self) -> bd.SearchSettings:
        return bd.SearchSettings(max_mismatches=self.max_mismatches,
                                 three_prime_exact=self.three_prime_exact, both_strands=True)


@dataclass
class PcrProduct:
    lane: str                 # template label (file name, accession, "methylated" …)
    template: str             # template FASTA path
    record: str
    start: int                # 1-based, top strand, binding parts only
    end: int
    forward: str
    reverse: str
    mm_forward: int = 0
    mm_reverse: int = 0
    size: int = 0             # incl. tails
    sites: int = 1            # identical products (same size) merged
    kind: str = "product"     # expected / off-target / ladder n× / dimer
    label: str = ""
    yield_: float = 0.0       # relative amount, 0..1 within the lane
    sequence: str = ""        # incl. tails ("" when not kept)
    locations: list[str] = field(default_factory=list)  # first few record:start-end

    def to_dict(self) -> dict:
        d = asdict(self)
        d["yield"] = d.pop("yield_")
        return d

    @classmethod
    def from_dict(cls, d: dict) -> PcrProduct:
        d = dict(d)
        d["yield_"] = d.pop("yield", d.get("yield_", 0.0))
        return cls(**{k: d[k] for k in cls.__dataclass_fields__ if k in d})


@dataclass
class PcrResult:
    settings: PcrSettings
    primers: list[Primer]
    templates: list[str]
    products: list[PcrProduct] = field(default_factory=list)
    lanes: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)
    monomer: int = 0          # ladder step found (0 = none)

    def to_dict(self) -> dict:
        return {
            "format": RESULT_FORMAT,
            "settings": self.settings.to_dict(),
            "primers": [p.to_dict() for p in self.primers],
            "templates": list(self.templates),
            "products": [p.to_dict() for p in self.products],
            "lanes": list(self.lanes),
            "warnings": list(self.warnings),
            "monomer": self.monomer,
        }

    @classmethod
    def from_dict(cls, d: dict) -> PcrResult:
        return cls(settings=PcrSettings.from_dict(d.get("settings")),
                   primers=[Primer.from_dict(x) for x in d.get("primers") or []],
                   templates=list(d.get("templates") or []),
                   products=[PcrProduct.from_dict(x) for x in d.get("products") or []],
                   lanes=list(d.get("lanes") or []), warnings=list(d.get("warnings") or []),
                   monomer=int(d.get("monomer") or 0))

    def lane_products(self, lane: str) -> list[PcrProduct]:
        return [p for p in self.products if p.lane == lane]


# --------------------------------------------------------------------------
# Yield model
# --------------------------------------------------------------------------


def product_efficiency(size: int, mm_f: int, mm_r: int, polymerase: str = "Taq") -> float:
    """Relative amplification of one product (no units)."""
    _rate, decay = POLYMERASES.get(polymerase, POLYMERASES["Taq"])
    mm = 0.35 ** (mm_f + mm_r)
    return mm * math.exp(-max(0, size) / decay)


def extension_seconds(longest: int, polymerase: str = "Taq") -> int:
    rate, _ = POLYMERASES.get(polymerase, POLYMERASES["Taq"])
    secs = max(15, int(math.ceil(longest / rate / 5.0) * 5))
    return secs


# --------------------------------------------------------------------------
# Running
# --------------------------------------------------------------------------


def _tails(primers: list[Primer]) -> dict[str, str]:
    return {p.name: th.split_tail(p.sequence)[0].upper() for p in primers}


def _products_on_record(rid: str, seq: str, sites: list[bd.BindingSite], s: PcrSettings,
                        tails: dict[str, str], lane: str, template: str,
                        keep_seq: bool) -> list[PcrProduct]:
    key = {(x.primer, x.strand, x.start): x for x in sites}
    prods = bd.predict_products(sites, min_size=s.min_product, max_size=s.max_product,
                                max_products=MAX_PRODUCTS)
    out = []
    by_end = {(x.primer, x.end): x for x in sites if x.strand == "-"}
    for p in prods:
        f = key.get((p.forward, "+", p.start))
        r = by_end.get((p.reverse, p.end))
        tf = tails.get(p.forward, "")
        tr = tails.get(p.reverse, "")
        size = p.size + len(tf) + len(tr)
        sequence = ""
        if keep_seq:
            sequence = tf + seq[p.start - 1:p.end].upper() + th.reverse_complement(tr)
        out.append(PcrProduct(
            lane=lane, template=template, record=rid, start=p.start, end=p.end,
            forward=p.forward, reverse=p.reverse,
            mm_forward=f.mismatches if f else 0, mm_reverse=r.mismatches if r else 0,
            size=size, sequence=sequence, locations=[f"{rid}:{p.start:,}-{p.end:,}"]))
    return out


def _merge_same(products: list[PcrProduct]) -> list[PcrProduct]:
    """One row per (lane, size, primer pair, mismatches); ``sites`` counts copies."""
    merged: dict[tuple, PcrProduct] = {}
    for p in products:
        k = (p.lane, p.size, p.forward, p.reverse, p.mm_forward + p.mm_reverse)
        if k in merged:
            m = merged[k]
            m.sites += 1
            if len(m.locations) < 5:
                m.locations.extend(p.locations)
        else:
            merged[k] = p
    return list(merged.values())


def find_monomer(sizes: list[int], *, tolerance: float = 0.04) -> int:
    """Common step between product sizes (a tandem repeat's monomer), or 0."""
    uniq = sorted(set(sizes))
    if len(uniq) < 3:
        return 0
    diffs = [b - a for a, b in zip(uniq, uniq[1:]) if b - a > 20]
    if not diffs:
        return 0
    best, best_n = 0, 0
    for d in diffs:
        n = sum(1 for x in diffs if abs(x - d) <= tolerance * d)
        if n > best_n or (n == best_n and d < best):
            best, best_n = d, n
    return best if best_n >= 2 else 0


def _label(products: list[PcrProduct], s: PcrSettings, monomer: int) -> None:
    """Kinds and labels per lane, then yields normalised per lane."""
    lanes: dict[str, list[PcrProduct]] = {}
    for p in products:
        lanes.setdefault(p.lane, []).append(p)
    for lane, prods in lanes.items():
        real = [p for p in prods if p.kind != "dimer"]
        if not real:
            continue
        best = min(real, key=lambda p: (p.mm_forward + p.mm_reverse,
                                        0 if p.forward != p.reverse else 1, p.size))
        base = best.size
        for p in real:
            mm = p.mm_forward + p.mm_reverse
            if s.mode == MODE_LADDER and monomer and p.forward == best.forward \
                    and p.reverse == best.reverse:
                k = round((p.size - base) / monomer) + 1
                if k >= 1 and abs(p.size - (base + (k - 1) * monomer)) <= max(6, 0.04 * monomer * k):
                    p.kind = "expected" if k == 1 else "ladder"
                    p.label = f"{k}×"
                    continue
            if p is best or (p.size == best.size and mm == best.mm_forward + best.mm_reverse
                             and p.forward == best.forward):
                p.kind = "expected"
            else:
                p.kind = "off-target"
        for p in real:
            p.yield_ = product_efficiency(p.size, p.mm_forward, p.mm_reverse, s.polymerase) * p.sites
        top = max((p.yield_ for p in prods), default=0.0) or 1.0
        for p in prods:
            if p.kind == "dimer":
                p.yield_ = 0.12
            else:
                p.yield_ = round(p.yield_ / top, 4)


def _dimer_products(primers: list[Primer], lanes: list[str], templates: dict[str, str],
                    s: PcrSettings) -> list[PcrProduct]:
    if not s.dimers:
        return []
    worst, worst_dg = None, 0.0
    for i, a in enumerate(primers):
        for b in primers[i:]:
            try:
                d = th.homodimer(a.sequence) if a is b else th.heterodimer(a.sequence, b.sequence)
            except Exception:  # noqa: BLE001 - primer3 missing
                return []
            if d.dg < worst_dg:
                worst, worst_dg = (a, b), d.dg
    if worst is None or worst_dg > -5.0:
        return []
    a, b = worst
    size = max(20, len(a.sequence) + len(b.sequence) - 8)
    out = []
    for lane in lanes:
        out.append(PcrProduct(lane=lane, template=templates.get(lane, ""), record="",
                              start=0, end=0, forward=a.name, reverse=b.name, size=size,
                              kind="dimer", label="primer-dimer",
                              locations=[f"{a.name} × {b.name}, ΔG {worst_dg:.1f}"]))
    return out


def lane_label(path: Path) -> str:
    name = Path(path).name
    for ext in (".fa.gz", ".fasta.gz", ".fa", ".fasta", ".fna", ".fas"):
        if name.lower().endswith(ext):
            return name[: -len(ext)]
    return Path(path).stem


def run_pcr(templates: list[Path], primers: list[Primer], settings: PcrSettings, *,
            progress: ProgressFn = None, cancel: CancelFn = None) -> PcrResult:
    """All products of ``primers`` on each template FASTA."""
    s = settings
    if not primers:
        raise th.PrimerError("Choose at least one primer.")
    if not templates:
        raise th.PrimerError("Choose at least one template FASTA.")
    result = PcrResult(settings=s, primers=list(primers),
                       templates=[str(Path(t)) for t in templates])
    tails = _tails(primers)
    pseqs = [(p.name, p.sequence) for p in primers]
    all_products: list[PcrProduct] = []
    lane_template: dict[str, str] = {}
    for tpl in templates:
        tpl = Path(tpl).expanduser()
        if not tpl.is_file():
            raise th.PrimerError(f"Template not found: {tpl}")
        base = lane_label(tpl)
        large = tpl.stat().st_size > bd.EXHAUSTIVE_MAX_BP
        keep = True

        for lane in (base,):
            result.lanes.append(lane)
            lane_template[lane] = str(tpl)

            def records():
                yield from bd.iter_fasta(tpl)

            found: list[PcrProduct] = []

            def on_record(rid, seq, sites, _lane=lane, _found=found):
                if len(_found) >= MAX_PRODUCTS:
                    return
                _found.extend(_products_on_record(rid, seq, sites, s, tails, _lane, str(tpl),
                                                  keep_seq=keep))

            bd.search_records(pseqs, records(), s.search(), large=large, progress=progress,
                              cancel=cancel, on_record=on_record)
            if len(found) >= MAX_PRODUCTS:
                result.warnings.append(f"{lane}: more than {MAX_PRODUCTS:,} products, the rest "
                                       "was not counted.")
            all_products.extend(found)
    merged = _merge_same(all_products)
    result.monomer = find_monomer([p.size for p in merged]) if s.mode == MODE_LADDER else 0
    if s.mode == MODE_LADDER and not result.monomer and merged:
        result.warnings.append("No ladder found: the products do not step by a common size.")
    merged += _dimer_products(primers, result.lanes, lane_template, s)
    _label(merged, s, result.monomer)
    merged.sort(key=lambda p: (result.lanes.index(p.lane) if p.lane in result.lanes else 0,
                               {"expected": 0, "ladder": 1, "off-target": 2}.get(p.kind, 3),
                               p.size))
    result.products = merged
    for lane in result.lanes:
        if not [p for p in result.lane_products(lane) if p.kind != "dimer"]:
            result.warnings.append(f"{lane}: no product.")
    return result


# --------------------------------------------------------------------------
# Writing
# --------------------------------------------------------------------------


def output_base(first_template: Path, primers: list[Primer]) -> str:
    names = [p.name for p in primers[:2]]
    tag = "+".join(names) if names else "primers"
    if len(names) == 2:
        a, b = names
        stem = _common_stem(a, b)
        if stem and len(stem) >= 3:
            tag = f"{a}+{b[len(stem):]}" if b[len(stem):] else a
    safe = "".join(ch if ch.isalnum() or ch in "+-_." else "_" for ch in tag)
    return f"{Path(first_template).name}{PCR_TAG}_{safe}"


def _common_stem(a: str, b: str) -> str:
    n = 0
    for x, y in zip(a, b):
        if x != y:
            break
        n += 1
    stem = a[:n]
    # keep "CroSat2_" (cut back to a separator)
    for sep in ("_", "-", "."):
        if sep in stem:
            return stem[: stem.rfind(sep) + 1]
    return ""


def unique_base(folder: Path, base: str) -> str:
    cand, n = base, 2
    while (folder / f"{cand}.fa").exists() or (folder / f"{cand}_details").exists():
        cand = f"{base}_{n}"
        n += 1
    return cand


def product_name(p: PcrProduct, i: int) -> str:
    tag = p.label.replace("×", "x") if p.label and p.kind != "dimer" else p.kind.replace(" ", "_")
    return f"{p.forward}+{p.reverse}_{tag}_{p.size}bp_{i}"


def _lift_features(template: Path, record: str, start: int, end: int, offset: int) -> list:
    """Template annotations inside [start, end] moved onto the product."""
    from repeatviewer.backends.sequence_features import (
        SequenceFeature,
        load_features_by_id_from_json,
    )

    try:
        feats = load_features_by_id_from_json(template, [record]).get(record) or []
    except Exception:  # noqa: BLE001 - no annotations
        return []
    out = []
    for f in feats:
        if f.feature_type == "primer_bind" or f.end < start or f.start > end:
            continue
        a = max(f.start, start) - start + 1 + offset
        b = min(f.end, end) - start + 1 + offset
        partial = f.start < start or f.end > end
        out.append(SequenceFeature(start=a, end=b, name=f.name + (" (partial)" if partial else ""),
                                   feature_type=f.feature_type, color=f.color, strand=f.strand,
                                   source="PCR"))
    return out


def write_outputs(result: PcrResult, *, folder: Path | None = None,
                  base: str | None = None) -> dict:
    """Products FASTA (+ annotations, primers) next to the first template, the rest in
    ``<base>_details/``: products.tsv, history.json, result.json.

    Returns ``{"fasta": path, "details": folder, "files": [...]}``.
    """
    from repeatviewer.backends.primers import on_sequence as osq
    from repeatviewer.backends.sequence_features import save_annotations_multi

    first = Path(result.templates[0])
    folder = Path(folder) if folder else first.parent
    base = unique_base(folder, base or output_base(first, result.primers))
    details = folder / f"{base}_details"
    details.mkdir(parents=True, exist_ok=True)
    fasta = folder / f"{base}.fa"
    written: list[tuple[str, PcrProduct]] = []
    with open(fasta, "w", encoding="utf-8") as fh:
        for i, p in enumerate(result.products, start=1):
            if p.kind == "dimer" or not p.sequence:
                continue
            if len(written) >= MAX_WRITTEN:
                break
            name = product_name(p, i)
            desc = f"{p.locations[0] if p.locations else p.record} lane={p.lane}"
            fh.write(f">{name} {desc}\n")
            for k in range(0, len(p.sequence), 80):
                fh.write(p.sequence[k:k + 80] + "\n")
            written.append((name, p))
    by_name = {pr.name: pr for pr in result.primers}
    feats: dict[str, list] = {}
    lengths: dict[str, int] = {}
    for name, p in written:
        lengths[name] = len(p.sequence)
        fl = []
        for pname, strand in ((p.forward, "+"), (p.reverse, "-")):
            pr = by_name.get(pname)
            if pr is None:
                continue
            hit = osq.place_primer(pr.sequence.upper(), p.sequence, name=pname)
            if hit is not None:
                fl.append(osq.primer_feature(pr, hit))
        tail_f = len(th.split_tail(by_name[p.forward].sequence)[0]) if p.forward in by_name else 0
        if p.record and "(" not in p.record:
            fl.extend(_lift_features(Path(p.template), p.record, p.start, p.end, tail_f))
        feats[name] = sorted(fl, key=lambda f: (f.start, f.end))
    if written:
        save_annotations_multi(fasta, features_by_id=feats, lengths_by_id=lengths,
                               record_ids=[n for n, _ in written])
        osq.save_sequence_primers(fasta, [Primer(name=pr.name, sequence=pr.sequence.upper(),
                                                 kind=pr.kind, note=pr.note, source="PCR")
                                          for pr in result.primers])
    tsv = details / f"{base}_products.tsv"
    with open(tsv, "w", encoding="utf-8") as fh:
        fh.write("lane\tkind\tlabel\tsize\tforward\treverse\tmismatches_forward\t"
                 "mismatches_reverse\tsites\trelative_yield\twhere\n")
        for p in result.products:
            fh.write("\t".join(str(x) for x in (
                p.lane, p.kind, p.label, p.size, p.forward, p.reverse, p.mm_forward,
                p.mm_reverse, p.sites, f"{p.yield_:.3f}", "; ".join(p.locations))) + "\n")
    history = {
        "made": _dt.datetime.now().isoformat(timespec="seconds"),
        "tool": "RepeatViewer in-silico PCR",
        "kind": MODE_LABELS.get(result.settings.mode, result.settings.mode),
        "templates": result.templates,
        "primers": [{"name": p.name, "sequence": p.sequence} for p in result.primers],
        "settings": result.settings.to_dict(),
        "products_written": len(written),
    }
    hist = details / f"{base}_history.json"
    hist.write_text(json.dumps(history, indent=1) + "\n", encoding="utf-8")
    res = details / f"{base}_result.json"
    res.write_text(json.dumps(result.to_dict(), indent=1) + "\n", encoding="utf-8")
    files = [p for p in (fasta if written else None, tsv, hist, res) if p is not None]
    if not written:
        try:
            fasta.unlink()
        except OSError:
            pass
    return {"fasta": fasta if written else None, "details": details, "base": base,
            "files": files, "result_json": res}


def load_result(path: Path) -> PcrResult:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("format") != RESULT_FORMAT:
        raise th.PrimerError(f"{Path(path).name} is not a RepeatViewer PCR result.")
    return PcrResult.from_dict(data)


def describe(result: PcrResult) -> str:
    """'7 products: 5 ladder steps, 1 off-target, 1 primer-dimer' style summary."""
    real = [p for p in result.products if p.kind != "dimer"]
    parts = []
    exp = sum(1 for p in real if p.kind == "expected")
    lad = sum(1 for p in real if p.kind == "ladder")
    off = sum(1 for p in real if p.kind == "off-target")
    dim = sum(1 for p in result.products if p.kind == "dimer")
    if exp:
        parts.append(f"{exp} expected")
    if lad:
        parts.append(f"{lad} ladder step{'s' if lad > 1 else ''}"
                     + (f" ({result.monomer} bp)" if result.monomer else ""))
    if off:
        parts.append(f"{off} off-target")
    if dim:
        parts.append("primer-dimer" if dim == len(result.lanes) else f"{dim} primer-dimer")
    n = len(real)
    head = f"{n} product{'s' if n != 1 else ''}"
    return head + (": " + ", ".join(parts) if parts else "")


# --------------------------------------------------------------------------
# Binding sites across a genome (Binding Sites on Chromosomes)
# --------------------------------------------------------------------------


@dataclass
class RecordSites:
    record: str
    length: int
    # Sites kept for drawing (all of them up to MAX_SITES_PER_RECORD).
    sites: list[tuple[str, str, int, int]] = field(default_factory=list)  # primer, strand, pos, mm
    # Products kept for drawing: start, end, size, left primer, right primer.
    products: list[tuple[int, int, int, str, str]] = field(default_factory=list)
    # Full counts (not limited by what is kept for drawing).
    site_counts: dict[str, dict[int, int]] = field(default_factory=dict)   # primer -> mm -> n
    pair_counts: dict[tuple[str, str], int] = field(default_factory=dict)  # (left, right) -> n
    sites_truncated: bool = False
    products_truncated: bool = False

    def counts(self, primers: set[str] | None = None) -> dict[int, int]:
        out: dict[int, int] = {}
        for name, by_mm in self.site_counts.items():
            if primers is not None and name not in primers:
                continue
            for mm, n in by_mm.items():
                out[mm] = out.get(mm, 0) + n
        return out

    def total_sites(self, primers: set[str] | None = None) -> int:
        return sum(self.counts(primers).values())

    def pair_products(self, left: str, right: str) -> int:
        """Products with ``left`` on the forward site and ``right`` on the reverse site
        (either order, since both primers go into the tube)."""
        n = self.pair_counts.get((left, right), 0)
        if left != right:
            n += self.pair_counts.get((right, left), 0)
        return n


MAX_SITES_PER_RECORD = 20000
MAX_PRODUCTS_PER_RECORD = 200000
MAX_PRODUCTS_DRAWN = 5000


def genome_sites(genome: Path, primers: list[Primer], *, max_mismatches: int = 2,
                 three_prime_exact: int = 3, max_product: int = 5000,
                 progress: ProgressFn = None, cancel: CancelFn = None) -> list[RecordSites]:
    """Sites of each primer on every record of a FASTA, plus products per record.

    Counts are complete; only what is kept for drawing is limited (the
    ``*_truncated`` flags say when).
    """
    genome = Path(genome).expanduser()
    if not genome.is_file():
        raise th.PrimerError(f"Genome FASTA not found: {genome}")
    s = bd.SearchSettings(max_mismatches=max_mismatches, three_prime_exact=three_prime_exact,
                          both_strands=True)
    out: list[RecordSites] = []

    def on_record(rid, seq, sites):
        rs = RecordSites(record=rid, length=len(seq))
        for x in sites:
            by_mm = rs.site_counts.setdefault(x.primer, {})
            by_mm[x.mismatches] = by_mm.get(x.mismatches, 0) + 1
        for x in sites[:MAX_SITES_PER_RECORD]:
            pos = x.start if x.strand == "+" else x.end
            rs.sites.append((x.primer, x.strand, pos, x.mismatches))
        rs.sites_truncated = len(sites) > MAX_SITES_PER_RECORD
        prods = bd.predict_products(sites, min_size=30, max_size=max_product,
                                    max_products=MAX_PRODUCTS_PER_RECORD)
        rs.products_truncated = len(prods) >= MAX_PRODUCTS_PER_RECORD
        for p in prods:
            key = (p.forward, p.reverse)
            rs.pair_counts[key] = rs.pair_counts.get(key, 0) + 1
        rs.products = [(p.start, p.end, p.size, p.forward, p.reverse)
                       for p in prods[:MAX_PRODUCTS_DRAWN]]
        out.append(rs)

    bd.search_records([(p.name, p.sequence) for p in primers], bd.iter_fasta(genome), s,
                      large=genome.stat().st_size > bd.EXHAUSTIVE_MAX_BP,
                      progress=progress, cancel=cancel, on_record=on_record)
    return out


def sites_gff3(records: list[RecordSites]) -> str:
    lines = ["##gff-version 3"]
    for rs in records:
        for primer, strand, pos, mm in rs.sites:
            lines.append(f"{rs.record}\tRepeatViewer\tprimer_bind\t{max(1, pos)}\t{max(1, pos)}\t"
                         f"{mm}\t{strand}\t.\tName={primer};mismatches={mm}")
    return "\n".join(lines) + "\n"
