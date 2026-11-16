"""Primer workbench logic (Design tab beside the sequence). Qt-free.

The workbench designs on the sequence that is open, around a region the user
drags on the map, clicks as an annotation or types. Rule (agreed 16 Nov 2026):
the settings drive the automatic design, the user's hand always wins. Nothing
here blocks a choice; :func:`region_advice` explains a conflict and offers a fix.
"""

from __future__ import annotations

from dataclasses import dataclass

from repeatviewer.backends.motif_search import parse_position
from repeatviewer.backends.primers import design as dz
from repeatviewer.backends.primers.oligos import Primer
from repeatviewer.backends.primers.on_sequence import primers_of_pair

USE_AMPLIFY = "amplify"
USE_INSIDE = "inside"
USE_AVOID = "avoid"
USES = (USE_AMPLIFY, USE_INSIDE, USE_AVOID)
USE_LABELS = {USE_AMPLIFY: "Amplify", USE_INSIDE: "Inside", USE_AVOID: "Avoid"}

#: Modes the workbench designs itself; Ladder and FISH stay in the Design window.
WORKBENCH_MODES = (dz.MODE_PCR, dz.MODE_QPCR)


@dataclass
class Region:
    start: int          # 1-based, inclusive
    end: int            # 1-based, inclusive
    use: str = USE_AMPLIFY

    def __post_init__(self) -> None:
        self.start, self.end = sorted((int(self.start), int(self.end)))
        if self.use not in USES:
            self.use = USE_AMPLIFY

    @property
    def length(self) -> int:
        return self.end - self.start + 1

    def label(self) -> str:
        return f"{self.start:,} – {self.end:,} · {USE_LABELS[self.use]} ({self.length:,} bp)"


def parse_region(text: str, seq_len: int, use: str = USE_AMPLIFY) -> Region | None:
    """ "120-380", "1,850..2,250" or "CroSat3:120-380" → Region (clipped to the sequence)."""
    pos = parse_position(text)
    if pos is None or seq_len < 1:
        return None
    _name, a, b = pos
    a, b = max(1, a), min(seq_len, b)
    if b < a:
        return None
    return Region(a, b, use)


def primer3_region_args(region: Region | None) -> dict:
    """Keyword arguments for :func:`design.run_primer3` (0-based half-open)."""
    if region is None:
        return {}
    a, b = region.start - 1, region.end
    if region.use == USE_AMPLIFY:
        return {"target": (a, b)}
    if region.use == USE_INSIDE:
        return {"included": (a, b)}
    return {"excluded": [(a, b)]}


@dataclass
class RegionAdvice:
    """Why Design would find nothing, and the settings that would fit."""

    message: str
    product_min: int | None = None
    product_max: int | None = None

    @property
    def can_fix(self) -> bool:
        return self.product_min is not None and self.product_max is not None


def region_advice(region: Region | None, s: dz.DesignSettings, seq_len: int) -> RegionAdvice | None:
    """A conflict between the region and the settings, or None when they fit."""
    if region is None:
        return None
    if region.use == USE_AMPLIFY:
        room_left = region.start - 1
        room_right = seq_len - region.end
        if room_left < s.size_min or room_right < s.size_min:
            return RegionAdvice(
                f"Primers must sit outside the region, but there are only {room_left} bp "
                f"before it and {room_right} bp after it (primers need {s.size_min} nt). "
                "Make the region smaller or use Inside.")
        need = region.length + 2 * s.size_min
        if s.product_max < need:
            return RegionAdvice(
                f"This {region.length:,} bp region needs products of at least {need:,} bp, "
                f"but the product range is {s.product_min:,} to {s.product_max:,} bp. "
                "Design would find nothing.",
                product_min=need, product_max=max(need + 200, need * 13 // 10))
        return None
    if region.use == USE_INSIDE:
        if region.length < s.product_min:
            return RegionAdvice(
                f"Primers and product must fit inside {region.length:,} bp, but the product "
                f"range starts at {s.product_min:,} bp.",
                product_min=max(2 * s.size_min + 10, region.length // 3),
                product_max=region.length)
        if region.length < 2 * s.size_min + 10:
            return RegionAdvice(f"{region.length} bp is too short for two primers.")
        return None
    free = seq_len - region.length
    if free < s.product_min:
        return RegionAdvice(
            f"Only {free:,} bp are left outside the avoided region, less than the product "
            f"range minimum ({s.product_min:,} bp).")
    return None


def fit_settings(s: dz.DesignSettings, advice: RegionAdvice) -> dz.DesignSettings:
    """A copy of ``s`` with the product range the advice suggests."""
    d = s.to_dict()
    if advice.can_fix:
        d["product_min"], d["product_max"] = advice.product_min, advice.product_max
    return dz.DesignSettings.from_dict(d)


def design_candidates(seq_id: str, seq: str, s: dz.DesignSettings,
                      region: Region | None) -> dz.Primer3Run:
    """Primer3 on the open sequence around the region; pairs come back as candidates."""
    return dz.run_primer3(seq_id or "sequence", seq, s, name_prefix="C",
                          **primer3_region_args(region))


def candidate_primers(pair: dz.DesignedPair) -> list[Primer]:
    """The primers to keep for a ticked candidate pair."""
    return primers_of_pair(pair.to_dict())


def covers(pair: dz.DesignedPair, region: Region | None) -> bool:
    """True when the pair's product spans the whole region."""
    if region is None or pair.forward is None or pair.reverse is None:
        return False
    return pair.forward.start <= region.start and pair.reverse.end >= region.end


def write_candidates_report(fasta, seq_id: str, seq: str, s: dz.DesignSettings,
                            region: Region | None, run: dz.Primer3Run):
    """The usual primer report for a workbench run (written only when asked)."""
    from pathlib import Path

    from repeatviewer.backends.primers import jobs as pj

    fasta = Path(fasta)
    payload = pj._base_payload("design", f"{dz.MODE_LABELS[s.mode]}: {fasta.name}", [fasta])
    payload.update(mode=s.mode, mode_label=dz.MODE_LABELS[s.mode], settings=s.to_dict(),
                   thermo=s.thermo.describe(),
                   templates=[pj._template_entry(seq_id, seq, file=str(fasta))],
                   pairs=[p.to_dict() for p in run.pairs], explain=list(run.explain))
    if region is not None:
        payload["explain"].insert(0, f"Region {region.label()}.")
    lane_sizes = {p.name: dz.pair_bands_on(p, [(seq_id, seq)]) for p in run.pairs}
    payload["gel"] = pj._gel_from_pairs(run.pairs, lane_sizes)
    if not run.pairs:
        payload["warnings"].append("No primers met the conditions. See \"Why\" below.")
    layout = pj.plan_outputs(fasta, pj.DESIGN_TAG)
    pj.write_result(layout, payload, pairs=run.pairs)
    return layout.main


#: Shortest binding part a primer can be dragged to (hard limit: no meaningful Tm below).
MIN_DRAG_NT = 10


def primer_for_span(primer_seq: str, placement, template: str, start: int, end: int) -> str:
    """The primer after dragging its binding part to ``start``..``end`` (1-based).

    Moving (same length, new place) takes the template bases there; resizing keeps
    the user's own bases (mismatches) where the old site overlaps. The 5' tail stays.
    """
    from repeatviewer.backends.primers import on_sequence as osq
    from repeatviewer.backends.primers import thermo as th

    tpl = template.upper()
    start, end = max(1, int(start)), min(len(tpl), int(end))
    if end - start + 1 < MIN_DRAG_NT:
        return primer_seq
    same_len = (end - start) == (placement.end - placement.start)
    if same_len and start == placement.start:
        return primer_seq
    if same_len:
        tail, _bind = th.split_tail(primer_seq)
        stretch = tpl[start - 1:end]
        return tail + (stretch if placement.strand == "+" else th.reverse_complement(stretch))
    if placement.strand == "+":
        five, three = placement.start - start, end - placement.end
    else:
        three, five = placement.start - start, end - placement.end
    return osq.resize_on_template(primer_seq, placement, template, five, three)


@dataclass
class SiteInfo:
    """One place a primer binds on the open sequence (for the site navigator)."""

    start: int
    end: int
    strand: str
    mismatches: int
    annealed: int
    tm: float
    three_exact: int


def binding_sites(primer_seq: str, template: str, *, max_mismatches: int = 2,
                  three_prime_exact: int = 3, settings=None) -> list[SiteInfo]:
    """Every site of the primer on the sequence, in sequence order, with its own Tm."""
    from repeatviewer.backends.primers import on_sequence as osq
    from repeatviewer.backends.primers import thermo as th

    bind = th.binding_part(primer_seq)
    n = len(bind)
    out = []
    for s in osq.find_all_sites(bind, template, max_mismatches=max_mismatches,
                                three_prime_exact=three_prime_exact):
        exact = n - max(s.mismatch_positions) if s.mismatch_positions else n
        try:
            tm = th.duplex_tm(bind, s.site, settings) if s.site else 0.0
        except Exception:  # noqa: BLE001 - Tm is informative
            tm = 0.0
        out.append(SiteInfo(s.start, s.end, s.strand, s.mismatches, s.end - s.start + 1,
                            float(tm), int(exact)))
    out.sort(key=lambda x: (x.start, x.strand))
    return out


#: A mismatch this close to the 3' end (1 = the last base) is flagged.
NEAR_THREE_PRIME = 5


def primer_flags(primer_seq: str, placement, s: dz.DesignSettings,
                 region: Region | None = None) -> list[str]:
    """Where a primer is outside the Design settings (shown in amber, never blocked)."""
    from repeatviewer.backends.primers import thermo as th

    flags = []
    bind = th.binding_part(primer_seq)
    if not bind:
        return flags
    tm = th.calc_tm(bind, s.thermo)
    if tm == tm:  # not NaN
        if tm < s.tm_min:
            flags.append(f"Tm {tm:.1f} °C, below {s.tm_min:g} to {s.tm_max:g}")
        elif tm > s.tm_max:
            flags.append(f"Tm {tm:.1f} °C, above {s.tm_min:g} to {s.tm_max:g}")
    n = len(bind)
    if n < s.size_min:
        flags.append(f"{n} nt, shorter than {s.size_min} to {s.size_max}")
    elif n > s.size_max:
        flags.append(f"{n} nt, longer than {s.size_min} to {s.size_max}")
    if placement is not None and placement.mismatch_positions:
        three = placement.end if placement.strand == "+" else placement.start
        nearest = min(abs(three - p) + 1 for p in placement.mismatch_positions)
        if nearest <= NEAR_THREE_PRIME:
            flags.append(f"Mismatch {nearest} nt from the 3′ end")
    if placement is not None and region is not None:
        overlaps = placement.start <= region.end and region.start <= placement.end
        if region.use == USE_AVOID and overlaps:
            flags.append(f"In the avoided region ({region.start:,} to {region.end:,})")
        if region.use == USE_INSIDE and not (region.start <= placement.start
                                             and placement.end <= region.end):
            flags.append(f"Outside the Inside region ({region.start:,} to {region.end:,})")
    return flags
