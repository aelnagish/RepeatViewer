"""Primers that live on a sequence (Preview > Primers tab, Edit Primer window).

Qt-free. A FASTA keeps its primers in a hidden sidecar next to it,
``{stem}.primers_RV.json`` (name, full oligo with a lower-case 5' tail, note).
Where a primer binds is not stored there: it is worked out from the sequence
(:func:`place_primer`) and drawn as a ``primer_bind`` annotation, the same
feature type Design and Map Primers already write.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass, field
from pathlib import Path

from repeatviewer.backends.primers import binding as bd
from repeatviewer.backends.primers import thermo as th
from repeatviewer.backends.primers.oligos import Primer, unique_name

PRIMERS_SUFFIX = ".primers_RV.json"
FORMAT_VERSION = 1

# Same family as Design / Map Primers (jobs.ANNOTATION_SOURCE etc.), kept here
# so this module does not import the job runner.
ANNOTATION_SOURCE = "primers"
FEATURE_PRIMER = "primer_bind"
COLOR_FORWARD = "#2f6fdb"
COLOR_REVERSE = "#d9730d"

MIN_ANNEALED = 12          # shortest 3' part accepted when a tail is guessed
PLACE_MAX_MISMATCHES = 4


# --------------------------------------------------------------------------
# Sidecar
# --------------------------------------------------------------------------


def primers_sidecar_path(fasta_path: Path | str) -> Path:
    """``{stem}.primers_RV.json`` beside the FASTA (same stem rule as the annotations)."""
    from repeatviewer.backends.sequence_features import sidecar_path

    seqmap = sidecar_path(Path(fasta_path))
    return seqmap.with_name(seqmap.name[: -len(".seqmap_RV.json")] + PRIMERS_SUFFIX)


def load_sequence_primers(fasta_path: Path | str) -> list[Primer]:
    path = primers_sidecar_path(fasta_path)
    if not path.is_file():
        return []
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    out: list[Primer] = []
    for item in data.get("primers") or []:
        try:
            out.append(Primer.from_dict(item))
        except (th.PrimerError, TypeError, ValueError):
            continue
    return out


def save_sequence_primers(fasta_path: Path | str, primers: list[Primer]) -> Path | None:
    """Write the sidecar; an empty list removes it."""
    path = primers_sidecar_path(fasta_path)
    if not primers:
        try:
            path.unlink()
        except OSError:
            pass
        return None
    payload = {"format_version": FORMAT_VERSION,
               "primers": [p.to_dict() for p in primers]}
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return path


# --------------------------------------------------------------------------
# Where a primer binds
# --------------------------------------------------------------------------


@dataclass
class Placement:
    """Best binding site of a primer on one sequence (annealed part only)."""

    name: str
    strand: str                 # "+" points right, "-" points left
    start: int                  # 1-based inclusive, top strand
    end: int
    mismatch_positions: list[int] = field(default_factory=list)  # template, 1-based
    annealed: int = 0           # length of the binding part
    tail: int = 0               # 5' tail length (does not bind)
    site: str = ""              # template stretch in the primer's orientation
    three_exact: int = 0        # exact bases counted from the 3' end
    other_sites: int = 0        # further sites with <= max mismatches

    @property
    def mismatches(self) -> int:
        return len(self.mismatch_positions)

    def to_dict(self) -> dict:
        return {k: getattr(self, k) for k in self.__dataclass_fields__}


def _site_to_placement(site: bd.BindingSite, annealed: int, tail: int) -> Placement:
    mm_t = []
    for i in site.mismatch_positions:
        mm_t.append(site.start + i - 1 if site.strand == "+" else site.end - i + 1)
    if site.mismatch_positions:
        three = annealed - max(site.mismatch_positions)
    else:
        three = annealed
    return Placement(name=site.primer, strand=site.strand, start=site.start, end=site.end,
                     mismatch_positions=sorted(mm_t), annealed=annealed, tail=tail,
                     site=site.site, three_exact=three)


def _rank(site: bd.BindingSite, annealed: int) -> tuple:
    three = annealed - max(site.mismatch_positions) if site.mismatch_positions else annealed
    return (site.mismatches, -min(three, 8), site.start)


def find_all_sites(primer_seq: str, template: str, *, max_mismatches: int = 2,
                   three_prime_exact: int = 3, name: str = "primer",
                   record: str = "seq") -> list[bd.BindingSite]:
    settings = bd.SearchSettings(max_mismatches=max_mismatches,
                                 three_prime_exact=three_prime_exact, both_strands=True)
    return bd.find_sites_in_sequence(name, primer_seq, record, template.upper(), settings)


def place_primer(primer_seq: str, template: str, *, max_mismatches: int = PLACE_MAX_MISMATCHES,
                 name: str = "primer") -> Placement | None:
    """Best site of the binding part (lower-case tail ignored), or None."""
    bind = th.binding_part(primer_seq)
    tail = len(primer_seq) - len(bind)
    if not bind or not template:
        return None
    sites = find_all_sites(primer_seq, template, max_mismatches=max_mismatches,
                           three_prime_exact=0, name=name)
    if not sites:
        return None
    sites.sort(key=lambda s: _rank(s, len(bind)))
    best = _site_to_placement(sites[0], len(bind), tail)
    best.other_sites = sum(1 for s in sites[1:] if s.mismatches <= 2)
    return best


def guess_tail(primer_seq: str, template: str, *, max_mismatches: int = 2,
               min_annealed: int = MIN_ANNEALED) -> str:
    """Primer with a lower-case 5' tail when only its 3' part binds the template.

    Typed oligos are usually all upper case; if the whole oligo does not bind
    but a 3' part of at least ``min_annealed`` bases does, the 5' rest is the
    tail. Returns the oligo unchanged when it binds as typed (or nowhere).
    """
    seq = th.clean_oligo(primer_seq)
    if th.split_tail(seq)[0]:
        return seq
    full = seq.upper()
    if place_primer(full, template, max_mismatches=max_mismatches) is not None:
        return full
    # Each mismatch costs as much as 3 bases of tail, so a lone mismatch near
    # the tail does not push the tail boundary past it.
    best, best_score = full, None
    for cut in range(1, len(full) - min_annealed + 1):
        part = full[cut:]
        hit = place_primer(part, template, max_mismatches=max_mismatches)
        if hit is None or hit.mismatches > max_mismatches:
            continue
        score = hit.mismatches * 3 + cut
        if best_score is None or score < best_score:
            best, best_score = full[:cut].lower() + part, score
    return best


# --------------------------------------------------------------------------
# Numbers shown in the side panel and the Edit Primer window
# --------------------------------------------------------------------------


@dataclass
class PrimerReport:
    length: int
    annealed: int
    tail: int
    gc: float
    tm_annealed: float          # on its site (mismatches included) or of the binding part
    tm_full: float              # whole oligo, from cycle 3 on
    hairpin: th.Structure
    self_dimer: th.Structure
    cross: dict[str, th.Structure] = field(default_factory=dict)  # other primer name -> dimer
    placement: Placement | None = None

    def status(self, dg: float, *, limit: float) -> str:
        """"ok" / "check" / "bad" for a ΔG (kcal/mol)."""
        if dg > limit / 2:
            return "ok"
        if dg > limit:
            return "check"
        return "bad"


def _safe(fn, *a, **k):
    try:
        return fn(*a, **k)
    except Exception:  # noqa: BLE001 - primer3 missing or odd input: show nothing
        return None


def analyse_primer(primer: Primer, template: str = "", *, others: list[Primer] | None = None,
                   settings: th.ThermoSettings | None = None) -> PrimerReport:
    seq = primer.sequence
    bind = th.binding_part(seq)
    place = place_primer(seq, template, name=primer.name) if template else None
    tm_ann = float("nan")
    if place is not None and place.mismatches and place.site:
        tm_ann = _safe(th.duplex_tm, bind, place.site, settings) or float("nan")
    if not tm_ann or math.isnan(tm_ann):
        tm_ann = _safe(th.calc_tm, bind, settings) or float("nan")
    tm_full = _safe(th.calc_tm, seq.upper(), settings) or float("nan")
    hp = _safe(th.hairpin, seq, settings) or th.Structure(kind="hairpin")
    sd = _safe(th.homodimer, seq, settings) or th.Structure(kind="self-dimer")
    cross = {}
    for o in others or []:
        if o.name == primer.name:
            continue
        d = _safe(th.heterodimer, seq, o.sequence, settings)
        if d is not None:
            cross[o.name] = d
    return PrimerReport(length=len(seq), annealed=len(bind), tail=len(seq) - len(bind),
                        gc=th.gc_percent(seq), tm_annealed=tm_ann, tm_full=tm_full,
                        hairpin=hp, self_dimer=sd, cross=cross, placement=place)


def dg_status(dg: float, *, limit: float = th.WARN_DIMER_DG) -> str:
    if dg >= limit / 2:
        return "ok"
    if dg >= limit:
        return "check"
    return "bad"


# --------------------------------------------------------------------------
# Edits
# --------------------------------------------------------------------------


def fix_mismatches(primer_seq: str, placement: Placement) -> str:
    """Binding part replaced by the template stretch (the tail is kept)."""
    tail, _bind = th.split_tail(primer_seq)
    if not placement.site:
        return primer_seq
    return tail + placement.site.upper()


def resize_on_template(primer_seq: str, placement: Placement, template: str,
                       five_delta: int, three_delta: int) -> str:
    """Extend (+) or trim (-) the annealed part at either end, following the template.

    ``five_delta`` / ``three_delta`` count bases at the primer's 5' / 3' end.
    The tail and the bases still inside the old site are kept as they were.
    """
    tail, bind = th.split_tail(primer_seq)
    tpl = template.upper()
    n = len(tpl)
    if placement.strand == "+":
        start = placement.start - five_delta
        end = placement.end + three_delta
    else:
        start = placement.start - three_delta
        end = placement.end + five_delta
    start = max(1, start)
    end = min(n, end)
    if end - start + 1 < 8:
        return primer_seq
    stretch = tpl[start - 1:end]
    new = stretch if placement.strand == "+" else th.reverse_complement(stretch)
    # Keep the user's own bases (mismatches) where the old site overlaps.
    old_off = (placement.start - start) if placement.strand == "+" else (end - placement.end)
    out = list(new)
    for i, ch in enumerate(bind):
        j = i + old_off
        if 0 <= j < len(out):
            out[j] = ch
    return tail + "".join(out)


def snap_to_tm(primer_seq: str, placement: Placement, template: str, target: float,
               *, settings: th.ThermoSettings | None = None, min_len: int = 15,
               max_len: int = 36) -> str:
    """Move the 5' end of the annealed part until its Tm is closest to ``target``.

    The 3' end stays where it is (it decides what gets copied).
    """
    best, best_gap = primer_seq, float("inf")
    cur = len(th.binding_part(primer_seq))
    for length in range(min_len, max_len + 1):
        cand = resize_on_template(primer_seq, placement, template, length - cur, 0)
        tm = _safe(th.calc_tm, th.binding_part(cand), settings)
        if tm is None or math.isnan(tm):
            continue
        gap = abs(tm - target)
        if gap < best_gap:
            best, best_gap = cand, gap
    return best


# --------------------------------------------------------------------------
# 5' tails
# --------------------------------------------------------------------------

TAIL_PRESETS = {
    "M13 forward (-20)": "GTAAAACGACGGCCAGT",
    "M13 reverse": "CAGGAAACAGCTATGAC",
    "T7 promoter": "TAATACGACTCACTATAGGG",
    "SP6 promoter": "ATTTAGGTGACACTATAG",
}

COMMON_ENZYMES = ("EcoRI", "BamHI", "HindIII", "KpnI", "SacI", "XhoI", "NotI", "XbaI",
                  "SalI", "PstI", "NcoI", "NdeI", "SpeI", "SmaI", "EcoRV", "BglII",
                  "NheI", "ApaI", "ClaI", "SphI", "MluI", "AscI", "PacI", "BsaI", "BsmBI")


def enzyme_site(name: str) -> str:
    """Recognition site of a restriction enzyme (Biopython's REBASE copy)."""
    from Bio import Restriction

    enz = getattr(Restriction, name, None)
    if enz is None:
        raise th.PrimerError(f"Unknown enzyme '{name}'.")
    return str(enz.site).upper()


def enzyme_label(name: str) -> str:
    """'KpnI  GGTAC^C' style label for the drop-down."""
    from Bio import Restriction

    enz = getattr(Restriction, name, None)
    if enz is None:
        return name
    site = str(enz.site)
    cut = enz.fst5
    if cut is not None and 0 <= cut <= len(site):
        site = site[:cut] + "^" + site[cut:]
    return f"{name}  {site}"


def enzyme_cuts(name: str, sequence: str, *, linear: bool = True) -> list[int]:
    """1-based cut positions of the enzyme on the sequence (top strand)."""
    from Bio import Restriction
    from Bio.Seq import Seq

    enz = getattr(Restriction, name, None)
    if enz is None or not sequence:
        return []
    return list(enz.search(Seq(sequence.upper()), linear=linear))


def enzymes_not_cutting(sequence: str, candidates=COMMON_ENZYMES) -> list[str]:
    return [e for e in candidates if not enzyme_cuts(e, sequence)]


def with_tail(primer_seq: str, tail: str) -> str:
    """Replace the 5' tail (lower case in the result)."""
    bind = th.binding_part(primer_seq)
    return th.clean_oligo(tail).lower() + bind


def restriction_tail(enzyme: str, padding: str = "AT") -> str:
    return (padding + enzyme_site(enzyme)).lower()


# --------------------------------------------------------------------------
# Features for the annotation view
# --------------------------------------------------------------------------


def primer_feature(primer: Primer, place: Placement):
    from repeatviewer.backends.sequence_features import SequenceFeature

    label = primer.name
    if place.tail:
        label += f" (+{place.tail} tail)"
    if place.mismatches:
        label += f" · {place.mismatches} mm"
    color = COLOR_FORWARD if place.strand == "+" else COLOR_REVERSE
    return SequenceFeature(start=place.start, end=place.end, name=label,
                           feature_type=FEATURE_PRIMER, color=color, strand=place.strand,
                           source=ANNOTATION_SOURCE)


def feature_primer_name(feat) -> str:
    """Primer name from a primer_bind feature label written by :func:`primer_feature`."""
    name = str(getattr(feat, "name", "") or "")
    for sep in (" (+", " · "):
        if sep in name:
            name = name.split(sep, 1)[0]
    return name.strip()


def is_primer_feature(feat) -> bool:
    return str(getattr(feat, "feature_type", "")) == FEATURE_PRIMER


def merge_primer_features(features: list, primers: list[Primer], placements: dict[str, Placement | None]) -> list:
    """Annotation list with the arrows of ``primers`` replaced (other features kept).

    Arrows of primers no longer in the list are removed when they came from
    this panel (source "primers"); arrows from other sources stay.
    """
    names = {p.name for p in primers}
    kept = []
    for f in features:
        if is_primer_feature(f) and (feature_primer_name(f) in names
                                     or str(getattr(f, "source", "")) == ANNOTATION_SOURCE):
            continue
        kept.append(f)
    for p in primers:
        place = placements.get(p.name)
        if place is not None:
            kept.append(primer_feature(p, place))
    return sorted(kept, key=lambda f: (f.start, f.end))


def primers_from_features(features: list, template: str, known: list[Primer]) -> list[Primer]:
    """Primers for arrows that have no stored sequence (older Design / Map runs).

    The sequence is read off the template at the arrow; such primers have no tail.
    """
    have = {p.name for p in known}
    out: list[Primer] = []
    for f in features:
        if not is_primer_feature(f):
            continue
        name = feature_primer_name(f)
        if not name or name in have:
            continue
        stretch = template[f.start - 1:f.end].upper()
        if not stretch:
            continue
        seq = stretch if f.is_forward() else th.reverse_complement(stretch)
        try:
            out.append(Primer(name=name, sequence=seq, source="annotation"))
        except th.PrimerError:
            continue
        have.add(name)
    return out


def new_primer_name(base: str, primers: list[Primer]) -> str:
    return unique_name(base, {p.name for p in primers})


# --------------------------------------------------------------------------
# Writing primers onto a FASTA on disk (Design results, PCR products)
# --------------------------------------------------------------------------


def add_primers_to_fasta(fasta: Path | str, primers: list[Primer], *,
                         record_id: str | None = None) -> dict[str, Placement | None]:
    """Keep ``primers`` with the FASTA (same name: replaced) and draw their arrows.

    Each primer is placed on ``record_id`` (or, when None, on the first record
    it binds). Returns the placements by primer name.
    """
    from repeatviewer.backends.sequence_features import (
        load_features_by_id_from_json,
        load_sequence_map_data,
        save_annotations_multi,
    )

    fasta = Path(fasta).expanduser().resolve()
    stored = load_sequence_primers(fasta)
    by_name = {p.name: p for p in stored}
    for p in primers:
        by_name[p.name] = p
    save_sequence_primers(fasta, list(by_name.values()))

    data = load_sequence_map_data(fasta)
    record_ids = list(data.record_ids)
    lengths = dict(data.lengths_by_id)
    feats = load_features_by_id_from_json(fasta, record_ids)
    seqs = {rid: seq for rid, seq in bd.iter_fasta(fasta)}
    # iter_fasta ids may carry the description; match on the first token too.
    by_token = {rid.split()[0]: seq for rid, seq in seqs.items()}
    placements: dict[str, Placement | None] = {}
    per_record: dict[str, list[Primer]] = {}
    for p in primers:
        targets = [record_id] if record_id else record_ids
        placements[p.name] = None
        for rid in targets:
            seq = seqs.get(rid) or by_token.get(str(rid).split()[0], "")
            hit = place_primer(p.sequence, seq, name=p.name) if seq else None
            if hit is not None:
                placements[p.name] = hit
                per_record.setdefault(rid, []).append(p)
                break
    merged = {}
    for rid in record_ids:
        current = list(feats.get(rid) or data.features_by_id.get(rid) or [])
        prims = per_record.get(rid, [])
        if prims:
            names = {p.name for p in prims}
            current = [f for f in current
                       if not (is_primer_feature(f) and feature_primer_name(f) in names)]
            current += [primer_feature(p, placements[p.name]) for p in prims]
            current.sort(key=lambda f: (f.start, f.end))
        merged[rid] = current
    save_annotations_multi(fasta, features_by_id=merged, lengths_by_id=lengths,
                           record_ids=record_ids)
    return placements


def primers_of_pair(pair: dict) -> list[Primer]:
    """Primer objects of a designed pair (``DesignedPair.to_dict()``)."""
    out = []
    for key in ("forward", "reverse", "probe"):
        o = pair.get(key)
        if not o or not o.get("sequence"):
            continue
        try:
            out.append(Primer(name=o.get("name") or f"{pair.get('name', 'pair')}_{key[0].upper()}",
                              sequence=o["sequence"], source="Primer Design",
                              note=f"{pair.get('name', '')}, product {pair.get('product_size', '')} bp"))
        except th.PrimerError:
            continue
    return out
