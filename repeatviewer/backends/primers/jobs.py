"""Run a primer job end to end and write its files.

* :func:`run_design` – Tools → Primers → Design Primers
* :func:`run_map`    – Tools → Primers → Map Primers to Sequences

Each run writes ``<input>_primers.html`` (or ``_primer_map.html``) next to
the first input: the report the preview opens. Its ``<base>_details/``
folder holds the primers as FASTA, the tables, the gel picture and the full
result as JSON. The main window asks before an existing report is
overwritten (Overwrite / Cancel / Rename), so the name can also be given.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import csv
import datetime as _dt
import json
from collections.abc import Callable
from dataclasses import asdict, dataclass, field
from pathlib import Path

from repeatviewer.backends.primers import binding as bd
from repeatviewer.backends.primers import design as dz
from repeatviewer.backends.primers import low_complexity as lc
from repeatviewer.backends.primers.gel import Lane, gel_svg, merge_bands, pick_ladder
from repeatviewer.backends.primers.oligos import (
    KIND_FORWARD,
    KIND_OLIGO,
    KIND_PROBE,
    KIND_REVERSE,
    Primer,
    PrimerDatabase,
    load_primers,
)
from repeatviewer.backends.primers.thermo import (
    PrimerError,
    ThermoSettings,
    analyze_oligo,
    binding_part,
    clean_oligo,
    duplex_tm,
    primer3_version,
    require_primer3,
)

REPORT_SUFFIX = ".html"
RESULT_SUFFIX = ".json"
RESULT_FORMAT = "repeatviewer-primers"
ANNOTATION_SOURCE = "primers"
FEATURE_PRIMER = "primer_bind"
FEATURE_PROBE = "probe"
FEATURE_PRODUCT = "PCR product"
COLOR_FORWARD = "#2e7d32"
COLOR_REVERSE = "#7b1fa2"
COLOR_PROBE = "#ef6c00"
COLOR_PRODUCT = "#90a4ae"
MAX_TEMPLATE_IN_FILE = 300_000
MAX_SITES_IN_FILE = 5000
MAX_PRODUCTS_IN_FILE = 5000
MAX_SITE_TM = 1500

ProgressFn = Callable[[str], None] | None
CancelFn = Callable[[], bool] | None


def _check(cancel: CancelFn) -> None:
    if cancel is not None and cancel():
        raise bd.PrimerCancelled()


def _say(progress: ProgressFn, text: str) -> None:
    if progress is not None:
        progress(text)


# --------------------------------------------------------------------------
# Output layout
# --------------------------------------------------------------------------


@dataclass
class OutputLayout:
    parent: Path
    base: str
    details: Path

    @property
    def main(self) -> Path:
        """The HTML report next to the input (opened in the preview)."""
        return self.parent / f"{self.base}{REPORT_SUFFIX}"

    @property
    def result(self) -> Path:
        return self.path(RESULT_SUFFIX)

    def path(self, suffix: str) -> Path:
        return self.details / f"{self.base}{suffix}"


DESIGN_TAG = "primers"
MAP_TAG = "primer_map"


def default_report_path(first_input: Path, tag: str, out_dir: Path | None = None) -> Path:
    """``<input name>_<tag>.html`` (full input name kept, e.g. ``gene.fa_primers.html``)."""
    first = Path(first_input).expanduser().resolve()
    parent = Path(out_dir).expanduser().resolve() if out_dir else first.parent
    return parent / f"{first.name}_{tag}{REPORT_SUFFIX}"


def plan_outputs(first_input: Path, tag: str, out_dir: Path | None = None,
                 *, report_path: Path | None = None) -> OutputLayout:
    """Where the result goes.

    ``report_path`` (already confirmed by the user: overwrite or a new name)
    is used as given. Without it the first free ``_2``, ``_3`` … name is taken.
    """
    if report_path is not None:
        rp = Path(report_path).expanduser().resolve()
        base = rp.name[: -len(REPORT_SUFFIX)] if rp.name.lower().endswith(REPORT_SUFFIX) else rp.name
        details = rp.parent / f"{base}_details"
        details.mkdir(parents=True, exist_ok=True)
        return OutputLayout(rp.parent, base, details)
    first = Path(first_input).expanduser().resolve()
    parent = Path(out_dir).expanduser().resolve() if out_dir else first.parent
    stem = f"{first.name}_{tag}"
    n = 1
    while True:
        base = stem if n == 1 else f"{stem}_{n}"
        details = parent / f"{base}_details"
        if not (parent / f"{base}{REPORT_SUFFIX}").exists():
            try:
                details.mkdir(parents=True, exist_ok=False)
                return OutputLayout(parent, base, details)
            except FileExistsError:
                pass
        n += 1


# --------------------------------------------------------------------------
# Inputs
# --------------------------------------------------------------------------


@dataclass
class TemplateSet:
    path: Path
    names: list[str]
    rows: list[str]               # as read (gaps kept)
    aligned: bool

    def ungapped(self) -> list[tuple[str, str]]:
        return [(n, r.replace("-", "").replace(".", "")) for n, r in zip(self.names, self.rows)]


def read_templates(path: Path, *, max_records: int | None = None) -> TemplateSet:
    names, rows = [], []
    for rid, seq in bd.iter_fasta(Path(path)):
        names.append(rid)
        rows.append(seq)
        if max_records is not None and len(names) >= max_records:
            break
    if not names:
        raise PrimerError(f"No sequences in {Path(path).name}.")
    return TemplateSet(Path(path), names, rows, bd.is_aligned(rows))


def record_features(path: Path) -> dict[str, list]:
    """Annotations RepeatViewer keeps for this FASTA."""
    try:
        from repeatviewer.backends.sequence_features import load_sequence_map_data

        data = load_sequence_map_data(Path(path))
        return {rid: list(f) for rid, f in data.features_by_id.items()}
    except Exception:  # noqa: BLE001 - no sidecar / unreadable: no features
        return {}


# --------------------------------------------------------------------------
# Annotations written back to the sequences
# --------------------------------------------------------------------------


def annotate_fasta(fasta: Path, features_by_record: dict[str, list]) -> list[Path]:
    """Replace earlier primer annotations of this FASTA with ``features_by_record``."""
    from repeatviewer.backends.sequence_features import (
        load_features_by_id_from_json,
        load_sequence_map_data,
        save_annotations_multi,
    )

    fasta = Path(fasta).expanduser().resolve()
    data = load_sequence_map_data(fasta)
    record_ids = list(data.record_ids)
    lengths = dict(data.lengths_by_id)
    json_feats = load_features_by_id_from_json(fasta, record_ids)
    merged = {}
    for rid in record_ids:
        current = list(json_feats.get(rid) or data.features_by_id.get(rid) or [])
        current = [f for f in current if str(f.source or "") != ANNOTATION_SOURCE]
        merged[rid] = current + list(features_by_record.get(rid, []))
    json_path, gff_path = save_annotations_multi(
        fasta, features_by_id=merged, lengths_by_id=lengths, record_ids=record_ids)
    return [p for p in (json_path, gff_path) if p is not None]


def _feature(start: int, end: int, name: str, ftype: str, color: str, strand: str):
    from repeatviewer.backends.sequence_features import SequenceFeature

    return SequenceFeature(start=start, end=end, name=name, feature_type=ftype, color=color,
                           strand=strand, source=ANNOTATION_SOURCE)


def pair_features(pair: dz.DesignedPair) -> list:
    feats = []
    for o in pair.oligos():
        if not o.start or not o.end:
            continue
        if o.kind == "probe":
            feats.append(_feature(o.start, o.end, o.name, FEATURE_PROBE, COLOR_PROBE, o.strand))
        else:
            color = COLOR_FORWARD if o.strand == "+" else COLOR_REVERSE
            feats.append(_feature(o.start, o.end, o.name, FEATURE_PRIMER, color, o.strand))
    te = pair.extra.get("te_primer")
    if te and te.get("start"):
        feats.append(_feature(te["start"], te["end"], te["name"], FEATURE_PRIMER,
                              COLOR_REVERSE, te.get("strand", "-")))
    if pair.forward is not None and pair.reverse is not None and pair.forward.start \
            and pair.reverse.end and pair.reverse.end > pair.forward.start \
            and pair.forward.strand == "+":
        feats.append(_feature(pair.forward.start, pair.reverse.end,
                              f"{pair.name} product ({pair.reverse.end - pair.forward.start + 1} bp)",
                              FEATURE_PRODUCT, COLOR_PRODUCT, "+"))
    return feats


# --------------------------------------------------------------------------
# Writing results
# --------------------------------------------------------------------------


def _template_entry(rid: str, seq: str, **extra) -> dict:
    entry = {"id": rid, "length": len(seq)}
    if len(seq) <= MAX_TEMPLATE_IN_FILE:
        entry["sequence"] = seq
    entry.update(extra)
    return entry


def _write_tsv(path: Path, header: list[str], rows: list[list]) -> Path:
    with open(path, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh, delimiter="\t", lineterminator="\n")
        w.writerow(header)
        w.writerows(rows)
    return path


def _primers_from_pairs(pairs: list[dz.DesignedPair]) -> list[Primer]:
    out = []
    seen = set()
    for p in pairs:
        oligos = list(p.oligos())
        te = p.extra.get("te_primer")
        if te:
            oligos.append(dz.DesignedOligo(**{k: te[k] for k in dz.DesignedOligo.__dataclass_fields__
                                              if k in te}))
        for o in oligos:
            if o.name in seen:
                continue
            seen.add(o.name)
            note = f"Tm {o.tm:.1f} °C, {o.length} nt"
            if p.product_size:
                note += f", product {p.product_size} bp"
            if p.extra.get("order"):
                note += f", order {p.extra['order']}"
            partners_ = [x.name for x in oligos if x is not o and x.kind != o.kind
                         and x.kind in ("forward", "reverse") and o.kind in ("forward", "reverse")]
            out.append(Primer(name=o.name, sequence=o.sequence, kind=o.kind, note=note,
                              source=p.template,
                              extra={"status": "designed", "target": p.template,
                                     **({"pair": ", ".join(partners_)} if partners_ else {})}))
    return out


def write_result(layout: OutputLayout, payload: dict, *, pairs: list[dz.DesignedPair] | None = None,
                 sites: list[bd.BindingSite] | None = None,
                 products: list[bd.Product] | None = None) -> list[Path]:
    from repeatviewer.backends.primers.report import write_report

    files: list[Path] = []
    if pairs:
        prim = _primers_from_pairs(pairs)
        from repeatviewer.backends.primers.oligos import write_primers

        files.append(write_primers(layout.path("_primers.fa"), prim))
        rows = []
        for p in pairs:
            for o in p.oligos():
                rows.append([p.rank, p.name, p.group, o.name, o.kind, o.sequence, o.start, o.end,
                             o.strand, o.length, f"{o.tm:.1f}", f"{o.gc:.1f}", f"{o.end_dg:.2f}",
                             f"{o.hairpin_dg:.2f}", f"{o.self_dimer_dg:.2f}", p.product_size,
                             f"{p.dimer_dg:.2f}", p.specificity.get("products", ""),
                             "; ".join(p.notes + o.warnings)])
        files.append(_write_tsv(layout.path("_pairs.tsv"), [
            "rank", "pair", "group", "oligo", "kind", "sequence", "start", "end", "strand",
            "length", "tm", "gc", "end_dG", "hairpin_dG", "self_dimer_dG", "product",
            "pair_dimer_dG", "genome_products", "notes"], rows))
    if sites:
        files.append(_write_tsv(layout.path("_binding_sites.tsv"), [
            "primer", "sequence", "strand", "start", "end", "mismatches",
            "mismatch_positions", "binding_tm", "site"],
            [[s.primer, s.record, s.strand, s.start, s.end, s.mismatches,
              ",".join(map(str, s.mismatch_positions)), s.tm or "", s.site] for s in sites]))
    if products:
        files.append(_write_tsv(layout.path("_products.tsv"), [
            "sequence", "forward", "reverse", "start", "end", "size", "mismatches"],
            [[p.record, p.forward, p.reverse, p.start, p.end, p.size, p.mismatches]
             for p in products]))
        with_seq = [p for p in products if p.sequence]
        if with_seq:
            fa = layout.path("_products.fa")
            with open(fa, "w", encoding="utf-8") as fh:
                for p in with_seq[:MAX_PRODUCTS_IN_FILE]:
                    fh.write(f">{p.record}:{p.start}-{p.end} {p.forward}+{p.reverse} "
                             f"{p.size}bp\n{p.sequence}\n")
            files.append(fa)
    gel = payload.get("gel") or {}
    if gel.get("lanes"):
        svg = gel_svg([Lane.from_dict(x) for x in gel["lanes"]], ladder=gel.get("ladder", ""),
                      title="Virtual gel")
        gpath = layout.path("_gel.svg")
        gpath.write_text(svg, encoding="utf-8")
        files.append(gpath)
    files.append(layout.result)
    payload["files"] = {p.name: str(p) for p in files}
    payload["report"] = str(layout.main)
    layout.result.write_text(json.dumps(payload, indent=1, default=_json_default) + "\n",
                             encoding="utf-8")
    write_report(layout.main, payload)
    return [layout.main] + files


def _json_default(obj):
    if hasattr(obj, "to_dict"):
        return obj.to_dict()
    if hasattr(obj, "item"):
        return obj.item()
    return str(obj)


def load_result(path: Path) -> dict:
    """The JSON result written into the details folder."""
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    if data.get("format") != RESULT_FORMAT:
        raise PrimerError(f"{Path(path).name} is not a RepeatViewer primer result.")
    return data


def _base_payload(kind: str, title: str, inputs: list[Path]) -> dict:
    return {
        "format": RESULT_FORMAT,
        "version": 1,
        "kind": kind,
        "title": title,
        "created": _dt.datetime.now().isoformat(timespec="seconds"),
        "inputs": [str(p) for p in inputs],
        "engine": f"Primer3 (primer3-py {primer3_version()})",
        "warnings": [],
        "explain": [],
    }


def _gel_from_pairs(pairs: list[dz.DesignedPair], lane_sizes: dict[str, list[int]],
                    *, max_lanes: int = 10) -> dict:
    lanes: list[Lane] = []
    all_sizes: list[int] = []
    for p in pairs[:max_lanes]:
        sizes = lane_sizes.get(p.name) or ([p.product_size] if p.product_size else [])
        if not sizes:
            continue
        bands = merge_bands(sizes)
        if p.product_size:
            for b in bands:
                if abs(b.size - p.product_size) <= max(3, p.product_size * 0.03):
                    b.highlight = True
        lanes.append(Lane(label=p.name, bands=bands))
        all_sizes.extend(sizes)
    if not lanes:
        return {}
    return {"ladder": pick_ladder(all_sizes), "lanes": [lane.to_dict() for lane in lanes]}


# --------------------------------------------------------------------------
# Design
# --------------------------------------------------------------------------


def run_design(inputs: list[Path], s: dz.DesignSettings, *, progress: ProgressFn = None,
               cancel: CancelFn = None, out_dir: Path | None = None,
               annotate: bool = True, max_templates: int = 20,
               report_path: Path | None = None,
               ) -> tuple[Path, list[Path], dict]:
    """Design primers on the selected FASTA files; returns (report, files, payload)."""
    require_primer3()
    inputs = [Path(p).expanduser().resolve() for p in inputs]
    if not inputs:
        raise PrimerError("Select a FASTA file in the Files panel first.")
    _say(progress, "Reading templates")
    whole = s.mode in (dz.MODE_LADDER, dz.MODE_FISH)
    sets = [read_templates(p, max_records=None if whole else max_templates) for p in inputs]
    payload = _base_payload("design", f"{dz.MODE_LABELS[s.mode]}: {inputs[0].name}", inputs)
    payload["mode"] = s.mode
    payload["mode_label"] = dz.MODE_LABELS[s.mode]
    payload["settings"] = s.to_dict()
    payload["thermo"] = s.hyb.describe() if s.mode == dz.MODE_FISH else s.thermo.describe()
    payload["templates"] = []
    pairs: list[dz.DesignedPair] = []
    lane_sizes: dict[str, list[int]] = {}
    explain: list[str] = []
    annotations: dict[Path, dict[str, list]] = {}

    for ts in sets:
        records = ts.ungapped()
        if s.mode in (dz.MODE_LADDER, dz.MODE_FISH):
            if ts.aligned:
                rid, tseq = "consensus", bd.consensus(ts.rows)
            else:
                # the copy most like the others stands in for a consensus
                short = [r for r in records if len(r[1]) <= 20_000] or records
                rid, tseq = dz.representative_record(short)
            mono, was_dimer = dz.monomer_from_dimer(tseq)
            if was_dimer:
                explain.append(f"{rid} is a dimer consensus (the monomer written twice, as in "
                               f"TAREAN's consensus_dimer.fasta): its first {len(mono)} bp "
                               "monomer was used.")
                tseq = mono
                if len(records) == 1:
                    records = [(rid, tseq)]
            payload["templates"].append(_template_entry(rid, tseq, file=str(ts.path)))
            _check(cancel)
            if s.mode == dz.MODE_LADDER:
                pairs.extend(_design_ladder(rid, tseq, records, s, lane_sizes, explain))
            else:
                avoid = []
                if s.fish_avoid_fasta and Path(s.fish_avoid_fasta).expanduser().is_file():
                    avoid = list(bd.iter_fasta(Path(s.fish_avoid_fasta).expanduser()))
                if s.fish_avoid_other_inputs:
                    for other in sets[1:]:
                        avoid.extend(other.ungapped())
                    if len(sets) > 1:
                        payload["avoided_files"] = [str(o.path) for o in sets[1:]]
                skip: list = []
                if s.fish_skip_unsuitable:
                    skip = lc.unsuitable_regions(tseq)
                    payload["unsuitable"] = [r.to_dict() for r in skip]
                members = records if len(records) > 1 else [(rid, tseq)]
                pairs.extend(dz.design_fish(rid, tseq, members, s, avoid, progress, skip=skip))
                if not pairs:
                    explain.append("No probe passed the coverage / Tm / GC / hairpin "
                                   "filters. Lower \"Min coverage\" or widen the length, "
                                   "Tm or GC range.")
            break   # one template set for these modes
        for rid, seq in records:
            _check(cancel)
            _say(progress, f"Designing on {rid} ({len(seq):,} bp)")
            payload["templates"].append(_template_entry(rid, seq, file=str(ts.path)))
            run = dz.run_primer3(rid, seq, s, pick_probe=s.pick_probe or None)
            found = run.pairs
            explain.extend(f"{rid}: {x}" for x in run.explain)
            for p in found:
                lane_sizes[p.name] = dz.pair_bands_on(p, [(rid, seq)])
            pairs.extend(found)
            if annotate and found:
                feats_out = []
                for p in found[:1]:
                    feats_out.extend(pair_features(p))
                annotations.setdefault(ts.path, {})[rid] = feats_out
        if s.mode == dz.MODE_QPCR and s.reference_fasta:
            pairs.extend(_design_reference(s, payload, lane_sizes, explain))
    pairs = _trim_per_template(pairs, s.num_return)
    payload["pairs"] = [p.to_dict() for p in pairs]
    payload["explain"] = explain
    payload["gel"] = _gel_from_pairs(pairs, lane_sizes)
    if not pairs:
        payload["warnings"].append("No primers met the conditions. See \"Why\" below.")
    layout = plan_outputs(inputs[0], DESIGN_TAG, out_dir, report_path=report_path)
    for fasta, feats in annotations.items():
        try:
            kept = {p.name for p in pairs}
            feats = {rid: [f for f in fl if any(f.name.startswith(k) for k in kept)]
                     for rid, fl in feats.items()}
            for path in annotate_fasta(fasta, feats):
                payload.setdefault("annotated", []).append(str(path))
        except Exception as exc:  # noqa: BLE001 - annotations are a bonus, never fail the run
            payload["warnings"].append(f"Could not add annotations to {fasta.name}: {exc}")
    _say(progress, "Writing results")
    files = write_result(layout, payload, pairs=pairs)
    return layout.main, files, payload


def _trim_per_template(pairs: list[dz.DesignedPair], n: int) -> list[dz.DesignedPair]:
    count: dict[tuple[str, str], int] = {}
    out = []
    for p in pairs:
        key = (p.template, p.group)
        count[key] = count.get(key, 0) + 1
        if count[key] <= n:
            out.append(p)
    return out


def _design_ladder(rid: str, tseq: str, records: list[tuple[str, str]], s: dz.DesignSettings,
                   lane_sizes: dict[str, list[int]], explain: list[str]) -> list[dz.DesignedPair]:
    """Ladder primers on the monomer ``tseq``, ranked by how many input copies they bind."""
    want = s.num_return
    s.num_return = min(100, want * 4)
    try:
        run = dz.design_ladder(rid, tseq, s)
    finally:
        s.num_return = want
    explain.extend(run.explain)
    mono = len(tseq)
    long_arrays = [(n, q) for n, q in records if len(q) >= 3 * mono]
    st = bd.SearchSettings(max_mismatches=s.coverage_mismatches,
                           three_prime_exact=s.coverage_three_prime_exact)
    for p in run.pairs:
        cov = 0
        for n_, q in records:
            f = bd.find_sites_in_sequence("F", p.forward.sequence, n_, q + q, st)
            r = bd.find_sites_in_sequence("R", p.reverse.sequence, n_, q + q, st)
            cov += bool(f) and bool(r)
        p.coverage = {"target_total": len(records), "target_amplified": cov}
        if len(records) > 1:
            p.notes.append(f"both primers bind {cov}/{len(records)} input copies")
        if long_arrays:
            p.extra["_lane"] = dz.pair_bands_on(p, long_arrays[:5], st)
        else:
            sizes = []
            for size, w in dz.ladder_bands(p.product_size, mono):
                sizes.extend([size] * max(1, int(round(w * 7))))
            p.extra["_lane"] = sizes
            p.extra["simulated_array"] = 8
    # most conserved first: primers that sit in every monomer
    ranked = sorted(run.pairs, key=lambda q: (
        -q.coverage.get("target_amplified", 0), q.penalty))[:want]
    out = []
    for p in dz.rerank(ranked, rename="LAD"):
        lane_sizes[p.name] = p.extra.pop("_lane", [])
        out.append(p)
    return out


def _design_reference(s: dz.DesignSettings, payload: dict,
                      lane_sizes: dict[str, list[int]], explain: list[str]) -> list[dz.DesignedPair]:
    ref = Path(s.reference_fasta).expanduser()
    if not ref.is_file():
        explain.append(f"Reference gene FASTA not found: {s.reference_fasta}")
        return []
    out = []
    for rid, seq in list(bd.iter_fasta(ref))[:1]:
        payload["templates"].append(_template_entry(rid, seq, file=str(ref), role="reference"))
        run = dz.run_primer3(rid, seq.replace("-", ""), s,
                             name_prefix="REF", pick_probe=s.pick_probe or None)
        explain.extend(f"{rid} (reference): {x}" for x in run.explain)
        for p in run.pairs:
            p.group = "reference"
            p.notes.append("single-copy reference gene pair")
            lane_sizes[p.name] = [p.product_size]
        out.extend(run.pairs[: s.num_return])
    return out


# --------------------------------------------------------------------------
# Map primers
# --------------------------------------------------------------------------

MAP_SPECIFIC = "specific"
MAP_DATABASE = "database"
TEST_ON_SEQUENCES = "sequences"
TEST_ON_CONSENSUS = "consensus"


@dataclass
class MapSettings:
    mode: str = MAP_SPECIFIC
    use_forward: bool = True
    forward: str = ""
    use_reverse: bool = False
    reverse: str = ""
    use_probe: bool = False
    probe: str = ""
    primer_file: str = ""           # database mode: FASTA; empty = Primer Database
    # Database mode: only these saved primers (empty = all), e.g. the ones selected
    # in the Database window for "Check a Sequence".
    primer_names: list = field(default_factory=list)
    search_forward: bool = True
    search_reverse: bool = True
    search_probe: bool = True
    pairs_only: bool = False
    use_included: bool = False
    included_start: int = 1
    included_end: int = 0
    use_target: bool = False
    target_start: int = 1
    target_end: int = 1
    use_product_range: bool = True
    product_min: int = 1
    product_max: int = 5000
    use_optimal: bool = False
    optimal_size: int = 0
    use_max_mismatches: bool = True
    max_mismatches: int = 2
    three_prime_exact: int = 3
    thermo: ThermoSettings = field(default_factory=ThermoSettings)
    test_on: str = TEST_ON_SEQUENCES
    consensus_threshold: float = 0.0
    consensus_ignore_gaps: bool = True
    check_genome: bool = False
    genome_fasta: str = ""
    genome_seed: int = bd.DEFAULT_SEED
    annotate: bool = True

    def __post_init__(self) -> None:
        if isinstance(self.thermo, dict):
            self.thermo = ThermoSettings.from_dict(self.thermo)
        if self.mode not in (MAP_SPECIFIC, MAP_DATABASE):
            self.mode = MAP_SPECIFIC
        if self.test_on not in (TEST_ON_SEQUENCES, TEST_ON_CONSENSUS):
            self.test_on = TEST_ON_SEQUENCES

    def to_dict(self) -> dict:
        d = asdict(self)
        d["thermo"] = self.thermo.to_dict()
        return d

    @classmethod
    def from_dict(cls, data: dict | None) -> MapSettings:
        data = dict(data or {})
        return cls(**{k: data[k] for k in cls.__dataclass_fields__ if k in data})

    def search(self) -> bd.SearchSettings:
        return bd.SearchSettings(
            max_mismatches=self.max_mismatches if self.use_max_mismatches else 0,
            three_prime_exact=self.three_prime_exact)


def resolve_primer(text: str, kind: str, db: PrimerDatabase | None, label: str) -> Primer:
    """A primer field: a Primer Database name, or a sequence typed in."""
    value = str(text or "").strip()
    if not value:
        raise PrimerError(f"{label}: choose a saved primer or type its sequence.")
    if db is not None:
        hit = db.get(value)
        if hit is not None:
            return Primer(name=hit.name, sequence=hit.sequence, kind=kind, note=hit.note,
                          source="Primer Database")
    try:
        seq = clean_oligo(value)
    except PrimerError as exc:
        raise PrimerError(f"{label}: \"{value}\" is neither a saved primer nor a sequence "
                          f"({exc})") from exc
    if len(binding_part(seq)) < 6:
        raise PrimerError(f"{label}: at least 6 bases are needed.")
    name = {"forward": "Forward primer", "reverse": "Reverse primer", "probe": "Probe"}[kind]
    return Primer(name=name, sequence=seq, kind=kind, source="typed")


def selection_slots(primers: list[Primer]) -> dict[str, str] | None:
    """Selected Database primers as Forward / Reverse / Probe boxes of Map Primers.

    A pair (or one primer, or a pair plus a probe) fills the boxes; anything that
    does not fit three boxes (two forwards, four primers …) returns None and is
    tested as "the selected primers" instead."""
    from repeatviewer.backends.primers.inventory import kind_of

    if not 1 <= len(primers) <= 3:
        return None
    slots: dict[str, str] = {}
    unknown: list[str] = []
    for p in primers:
        k = kind_of(p.name, p.kind)
        if k in (KIND_FORWARD, KIND_REVERSE, KIND_PROBE):
            if k in slots:
                return None
            slots[k] = p.name
        else:
            unknown.append(p.name)
    for name in unknown:
        free = next((k for k in (KIND_FORWARD, KIND_REVERSE) if k not in slots), None)
        if free is None:
            return None
        slots[free] = name
    return slots


def map_primer_list(ms: MapSettings, db: PrimerDatabase | None) -> list[Primer]:
    if ms.mode == MAP_SPECIFIC:
        out = []
        if ms.use_forward:
            out.append(resolve_primer(ms.forward, KIND_FORWARD, db, "Forward primer"))
        if ms.use_reverse:
            out.append(resolve_primer(ms.reverse, KIND_REVERSE, db, "Reverse primer"))
        if ms.use_probe:
            out.append(resolve_primer(ms.probe, KIND_PROBE, db, "Probe"))
        if not out:
            raise PrimerError("Tick at least one of Forward primer, Reverse primer or Probe.")
        names = [p.name for p in out]
        if len(set(names)) < len(names):
            for i, p in enumerate(out):
                p.name = f"{p.name}_{i + 1}"
        return out
    if ms.primer_file:
        primers = load_primers(Path(ms.primer_file))
        where = Path(ms.primer_file).name
    else:
        primers = list(db.primers) if db is not None else []
        where = "the Primer Database"
    if not primers:
        raise PrimerError(f"No primers in {where}.")
    if ms.primer_names and not ms.primer_file:
        names = set(ms.primer_names)
        chosen = [p for p in primers if p.name in names]
        if not chosen:
            raise PrimerError("The selected primers are no longer in the Primer Database.")
        return chosen
    wanted = set()
    if ms.search_forward:
        wanted |= {KIND_FORWARD, KIND_OLIGO}
    if ms.search_reverse:
        wanted |= {KIND_REVERSE, KIND_OLIGO}
    if ms.search_probe:
        wanted |= {KIND_PROBE}
    chosen = [p for p in primers if p.kind in wanted]
    if not chosen:
        raise PrimerError("No saved primer of the ticked kinds (Forward / Reverse / Probe).")
    return chosen


def run_map(inputs: list[Path], ms: MapSettings, *, database_path: Path | None = None,
            progress: ProgressFn = None, cancel: CancelFn = None,
            out_dir: Path | None = None,
            report_path: Path | None = None) -> tuple[Path, list[Path], dict]:
    """Test primers against the selected sequences (Geneious' Map Primers)."""
    inputs = [Path(p).expanduser().resolve() for p in inputs]
    if not inputs:
        raise PrimerError("Select the sequences to test in the Files panel first.")
    from repeatviewer.backends.primers.inventory import open_database

    db = open_database(database_path)
    primers = map_primer_list(ms, db)
    _say(progress, f"Testing {len(primers)} primer(s)")
    payload = _base_payload("map", f"Map Primers: {inputs[0].name}", inputs)
    payload["engine"] = "RepeatViewer primer search; Tm by Primer3"
    payload["settings"] = ms.to_dict()
    payload["thermo"] = ms.thermo.describe()
    stats = []
    for p in primers:
        st = analyze_oligo(p.sequence, ms.thermo)
        d = p.to_dict()
        d.update({"tm": st.tm, "gc": st.gc, "length": st.binding_length, "end_dg": st.end_dg,
                  "hairpin_dg": st.hairpin.dg, "self_dimer_dg": st.self_dimer.dg,
                  "warnings": st.warnings, "hairpin_lines": st.hairpin.lines,
                  "self_dimer_lines": st.self_dimer.lines, "degeneracy": st.degeneracy})
        stats.append(d)
    payload["primers"] = stats
    templates: list[tuple[str, str, Path]] = []
    for path in inputs:
        ts = read_templates(path)
        if ms.test_on == TEST_ON_CONSENSUS and ts.aligned:
            cons = bd.consensus(ts.rows, threshold=ms.consensus_threshold,
                                ignore_gaps=ms.consensus_ignore_gaps)
            templates.append((f"{path.name} consensus", cons, path))
        else:
            if ms.test_on == TEST_ON_CONSENSUS:
                payload["warnings"].append(f"{path.name} is not an alignment: tested on its "
                                           "sequences instead of a consensus.")
            templates.extend((rid, seq, path) for rid, seq in ts.ungapped())
    total = sum(len(t[1]) for t in templates)
    large = total > bd.EXHAUSTIVE_MAX_BP
    settings = ms.search()
    sites: list[bd.BindingSite] = []
    seqs: dict[str, str] = {}
    offsets: dict[str, int] = {}
    for rid, seq, _path in templates:
        _check(cancel)
        lo, hi = 0, len(seq)
        if ms.use_included:
            lo = max(0, ms.included_start - 1)
            hi = min(len(seq), ms.included_end) if ms.included_end > 0 else len(seq)
        sub = seq[lo:hi]
        seqs[rid] = seq
        offsets[rid] = lo
        found = bd.search_records([(p.name, p.sequence) for p in primers], [(rid, sub)],
                                  settings, large=large, seed=ms.genome_seed, progress=progress,
                                  cancel=cancel)
        for s_ in found:
            s_.start += lo
            s_.end += lo
        sites.extend(found)
    for s_ in sites[:MAX_SITE_TM]:
        pr = next(p for p in primers if p.name == s_.primer)
        s_.tm = duplex_tm(pr.sequence, s_.site, ms.thermo)
    lo_size = ms.product_min if ms.use_product_range else 1
    hi_size = ms.product_max if ms.use_product_range else 50_000
    allowed = None
    if ms.mode == MAP_SPECIFIC:
        names = [p.name for p in primers if p.kind != KIND_PROBE]
        allowed = {(a, b) for a in names for b in names}
    products = bd.predict_products([s for s in sites if _kind_of(primers, s.primer) != KIND_PROBE],
                                   min_size=lo_size, max_size=hi_size, sequences=seqs,
                                   allowed_pairs=allowed)
    if ms.use_target:
        products = [p for p in products if p.start <= ms.target_start and p.end >= ms.target_end]
    if ms.use_optimal and ms.optimal_size > 0:
        products.sort(key=lambda p: (abs(p.size - ms.optimal_size), p.mismatches))
    if ms.mode == MAP_DATABASE and ms.pairs_only:
        used = {p.forward for p in products} | {p.reverse for p in products}
        sites = [s for s in sites if s.primer in used]
    if ms.mode == MAP_DATABASE:
        payload["stock"] = stock_summary(primers, sites, products, len(templates))
    payload["templates"] = [_template_entry(rid, seq, file=str(path)) for rid, seq, path
                            in templates[:200]]
    payload["sites"] = [s.to_dict() for s in sites[:MAX_SITES_IN_FILE]]
    payload["products"] = [p.to_dict() for p in products[:MAX_PRODUCTS_IN_FILE]]
    payload["counts"] = {"primers": len(primers), "sequences": len(templates),
                         "sites": len(sites), "products": len(products),
                         "sequences_with_product": len({p.record for p in products})}
    if len(sites) > MAX_SITES_IN_FILE:
        payload["warnings"].append(f"{len(sites):,} binding sites: the first "
                                   f"{MAX_SITES_IN_FILE:,} are listed (all in the TSV).")
    lanes = []
    per_record: dict[str, list[int]] = {}
    for p in products:
        per_record.setdefault(p.record, []).append(p.size)
    for rid, sizes in list(per_record.items())[:12]:
        lanes.append(Lane(label=rid, bands=merge_bands(sizes)))
    if lanes:
        all_sizes = [p.size for p in products]
        payload["gel"] = {"ladder": pick_ladder(all_sizes), "lanes": [x.to_dict() for x in lanes]}
    if ms.check_genome and ms.genome_fasta:
        _say(progress, f"Searching {Path(ms.genome_fasta).name}")
        payload["genome"] = _map_genome(primers, ms, progress, cancel)
    layout = plan_outputs(inputs[0], MAP_TAG, out_dir, report_path=report_path)
    if ms.annotate and sites:
        by_file: dict[Path, dict[str, list]] = {}
        rid_file = {rid: path for rid, _seq, path in templates}
        cons_ids = {rid for rid, _s, _p in templates if rid.endswith(" consensus")}
        for s_ in sites:
            if s_.record in cons_ids:
                continue
            kind = _kind_of(primers, s_.primer)
            color = COLOR_PROBE if kind == KIND_PROBE else (
                COLOR_FORWARD if s_.strand == "+" else COLOR_REVERSE)
            label = s_.primer + (f" ({s_.mismatches} mm)" if s_.mismatches else "")
            by_file.setdefault(rid_file[s_.record], {}).setdefault(s_.record, []).append(
                _feature(s_.start, s_.end, label,
                         FEATURE_PROBE if kind == KIND_PROBE else FEATURE_PRIMER, color, s_.strand))
        for p in products[:500]:
            if p.record in cons_ids:
                continue
            by_file.setdefault(rid_file[p.record], {}).setdefault(p.record, []).append(
                _feature(p.start, p.end, f"{p.forward}+{p.reverse} ({p.size} bp)",
                         FEATURE_PRODUCT, COLOR_PRODUCT, "+"))
        for fasta, feats in by_file.items():
            try:
                for path in annotate_fasta(fasta, feats):
                    payload.setdefault("annotated", []).append(str(path))
            except Exception as exc:  # noqa: BLE001
                payload["warnings"].append(f"Could not add annotations to {fasta.name}: {exc}")
    _say(progress, "Writing results")
    files = write_result(layout, payload, sites=sites, products=products)
    return layout.main, files, payload


STOCK_FIELDS = ("box", "freezer", "room", "status", "tested", "organism", "target",
                "ordered_by", "ordered_on", "anneal_temp")


def primer_location(p: Primer) -> str:
    """``box BI-8 · -20, 'general stuff' shelf · Sterilbankraum`` (what is known)."""
    e_ = p.extra or {}
    parts = []
    if e_.get("box"):
        parts.append(f"box {e_['box']}")
    for k in ("freezer", "room"):
        if e_.get(k):
            parts.append(str(e_[k]))
    return " · ".join(parts)


def stock_summary(primers: list[Primer], sites: list[bd.BindingSite],
                  products: list[bd.Product], n_sequences: int) -> dict:
    """What the saved primers offer for these sequences.

    * ready pairs: two saved primers that give a product (a saved pair first)
    * single primers: bind, but no saved partner gives a product
    """
    from repeatviewer.backends.primers.inventory import is_saved_pair, partners

    by_name = {p.name: p for p in primers}
    pairs_of = partners(primers)
    site_mm: dict[tuple, int] = {}
    for s_ in sites:
        key = (s_.primer, s_.record, s_.strand, s_.start if s_.strand == "+" else s_.end)
        site_mm[key] = min(s_.mismatches, site_mm.get(key, 99))
    ready: dict[tuple[str, str], dict] = {}
    for pr in products:
        mm_f = site_mm.get((pr.forward, pr.record, "+", pr.start), 0)
        mm_r = site_mm.get((pr.reverse, pr.record, "-", pr.end), 0)
        key = tuple(sorted((pr.forward, pr.reverse)))
        entry = ready.setdefault(key, {
            "forward": pr.forward, "reverse": pr.reverse,
            "saved_pair": is_saved_pair(pr.forward, pr.reverse, pairs_of),
            "mm_forward": mm_f, "mm_reverse": mm_r, "sizes": [], "records": set(),
            "products": 0})
        if mm_f + mm_r < entry["mm_forward"] + entry["mm_reverse"]:
            entry.update(forward=pr.forward, reverse=pr.reverse, mm_forward=mm_f, mm_reverse=mm_r)
        entry["products"] += 1
        entry["records"].add(pr.record)
        if pr.size not in entry["sizes"] and len(entry["sizes"]) < 6:
            entry["sizes"].append(pr.size)
    ready_list = []
    for entry in ready.values():
        entry["records"] = len(entry["records"])
        entry["sizes"].sort()
        for side in ("forward", "reverse"):
            p = by_name.get(entry[side])
            entry[f"{side}_where"] = primer_location(p) if p else ""
            entry[f"{side}_tested"] = (p.extra or {}).get("tested", "") if p else ""
        ready_list.append(entry)
    ready_list.sort(key=lambda d: (not d["saved_pair"], d["mm_forward"] + d["mm_reverse"],
                                   -d["records"], d["products"], d["sizes"][0] if d["sizes"] else 0))
    in_product = {d["forward"] for d in ready_list} | {d["reverse"] for d in ready_list}
    single: dict[str, dict] = {}
    for s_ in sites:
        if s_.primer in in_product:
            continue
        d = single.setdefault(s_.primer, {"name": s_.primer, "mismatches": s_.mismatches,
                                          "sites": 0, "record": s_.record,
                                          "strand": s_.strand, "start": s_.start})
        d["sites"] += 1
        if s_.mismatches < d["mismatches"]:
            d.update(mismatches=s_.mismatches, record=s_.record, strand=s_.strand,
                     start=s_.start)
    single_list = sorted(single.values(), key=lambda d: (d["mismatches"], -d["sites"]))
    for d in single_list:
        p = by_name.get(d["name"])
        d["where"] = primer_location(p) if p else ""
        d["partners"] = pairs_of.get(d["name"], [])
    return {"ready": ready_list[:200], "ready_total": len(ready_list),
            "single": single_list[:200], "single_total": len(single_list),
            "primers": len(primers), "sequences": n_sequences}


def _kind_of(primers: list[Primer], name: str) -> str:
    return next((p.kind for p in primers if p.name == name), KIND_OLIGO)


def _map_genome(primers: list[Primer], ms: MapSettings, progress: ProgressFn,
                cancel: CancelFn) -> dict:
    gpath = Path(ms.genome_fasta).expanduser()
    if not gpath.is_file():
        raise PrimerError(f"Genome FASTA not found: {ms.genome_fasta}")
    counts = {p.name: 0 for p in primers}
    n_products = 0
    sample: list[list] = []
    names = [p.name for p in primers if p.kind != KIND_PROBE]
    allowed = {(a, b) for a in names for b in names}
    hi = ms.product_max if ms.use_product_range else 5000
    lo = ms.product_min if ms.use_product_range else 1

    def on_record(rid, seq, sites):
        nonlocal n_products
        for s_ in sites:
            counts[s_.primer] = counts.get(s_.primer, 0) + 1
        prods = bd.predict_products(sites, min_size=lo, max_size=hi, allowed_pairs=allowed)
        n_products += len(prods)
        for p in prods[: max(0, 50 - len(sample))]:
            sample.append([rid, p.forward, p.reverse, p.start, p.end, p.size])

    bd.search_records([(p.name, p.sequence) for p in primers], bd.iter_fasta(gpath),
                      ms.search(), large=gpath.stat().st_size > bd.EXHAUSTIVE_MAX_BP,
                      seed=ms.genome_seed, progress=progress, cancel=cancel, on_record=on_record)
    return {"genome": str(gpath), "sites": counts, "products": n_products,
            "product_list": sample}


def save_pairs_to_database(payload: dict, db_path: Path | None = None,
                           pair_names: list[str] | None = None) -> list[Primer]:
    """Primers of a design result → Primer Database (identical sequences skipped)."""
    from repeatviewer.backends.primers.inventory import open_database

    db = open_database(db_path)
    pairs = []
    for d in payload.get("pairs", []):
        if pair_names and d.get("name") not in pair_names:
            continue
        p = dz.DesignedPair(rank=d.get("rank", 0), name=d.get("name", ""),
                            template=d.get("template", ""), product_size=d.get("product_size", 0),
                            extra=dict(d.get("extra") or {}))
        for key in ("forward", "reverse", "probe"):
            o = d.get(key)
            if o:
                setattr(p, key, dz.DesignedOligo(**{k: o[k] for k in
                                                   dz.DesignedOligo.__dataclass_fields__ if k in o}))
        pairs.append(p)
    added = db.add(_primers_from_pairs(pairs))
    db.save()
    return added
