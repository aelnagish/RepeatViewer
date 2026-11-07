"""Bench card numbers: mastermix and PCR program. Qt-free."""

from __future__ import annotations

import math
from dataclasses import dataclass

from repeatviewer.backends.primers import thermo as th
from repeatviewer.backends.primers.oligos import Primer

# Polymerase recipes (final concentrations), per the usual supplier protocols.
RECIPES = {
    "Taq": {"buffer_x": 10, "mg_mM": 1.5, "mg_stock_mM": 25.0, "dntp_mM": 0.2,
            "primer_uM": 0.2, "enzyme_u_per_ul": 0.025, "enzyme_stock": 5.0,
            "denature": (95, 180), "cycle_denature": (95, 30), "final": (72, 300),
            "rate_s_per_kb": 60, "anneal_offset": -3.0, "enzyme": "Taq 5 U/µl"},
    "Phusion / Q5": {"buffer_x": 5, "mg_mM": 0.0, "mg_stock_mM": 0.0, "dntp_mM": 0.2,
                     "primer_uM": 0.5, "enzyme_u_per_ul": 0.02, "enzyme_stock": 2.0,
                     "denature": (98, 30), "cycle_denature": (98, 10), "final": (72, 300),
                     "rate_s_per_kb": 20, "anneal_offset": 3.0, "enzyme": "Phusion 2 U/µl"},
    "Long-range": {"buffer_x": 10, "mg_mM": 0.0, "mg_stock_mM": 0.0, "dntp_mM": 0.35,
                   "primer_uM": 0.3, "enzyme_u_per_ul": 0.04, "enzyme_stock": 5.0,
                   "denature": (94, 120), "cycle_denature": (94, 15), "final": (68, 600),
                   "rate_s_per_kb": 60, "anneal_offset": -3.0, "enzyme": "Long-range 5 U/µl"},
}


@dataclass
class MixRow:
    reagent: str
    per_reaction: float   # µl
    total: float          # µl


def mastermix(primers: list[Primer], *, reactions: int = 10, extra_percent: float = 10.0,
              volume_ul: float = 25.0, template_ul: float = 1.0, polymerase: str = "Taq",
              primer_stock_uM: float = 10.0, dntp_stock_mM: float = 10.0) -> list[MixRow]:
    """Volumes per reaction and for ``reactions`` (+ extra %), template added separately."""
    r = RECIPES.get(polymerase, RECIPES["Taq"])
    n = reactions * (1 + extra_percent / 100.0)
    rows: list[tuple[str, float]] = []
    rows.append((f"{r['buffer_x']}× buffer", volume_ul / r["buffer_x"]))
    if r["mg_mM"] and r["mg_stock_mM"]:
        rows.append((f"MgCl₂ {r['mg_stock_mM']:g} mM", volume_ul * r["mg_mM"] / r["mg_stock_mM"]))
    rows.append((f"dNTP {dntp_stock_mM:g} mM each", volume_ul * r["dntp_mM"] / dntp_stock_mM))
    for p in primers[:4]:
        rows.append((f"{p.name} {primer_stock_uM:g} µM",
                     volume_ul * r["primer_uM"] / primer_stock_uM))
    rows.append((r["enzyme"], volume_ul * r["enzyme_u_per_ul"] / r["enzyme_stock"]))
    used = sum(v for _n, v in rows) + template_ul
    water = max(0.0, volume_ul - used)
    out = [MixRow("Water", round(water, 2), round(water * n, 1))]
    out += [MixRow(name, round(v, 2), round(v * n, 1)) for name, v in rows]
    return out


@dataclass
class Step:
    name: str
    temp_c: float
    seconds: int
    in_cycle: bool = False


def suggested_annealing(primers: list[Primer], polymerase: str = "Taq",
                        settings: th.ThermoSettings | None = None) -> tuple[float, list[float]]:
    tms = []
    for p in primers[:2]:
        try:
            tms.append(th.calc_tm(p.binding, settings))
        except Exception:  # noqa: BLE001 - primer3 missing
            continue
    tms = [t for t in tms if t and not math.isnan(t)]
    if not tms:
        return 58.0, []
    r = RECIPES.get(polymerase, RECIPES["Taq"])
    ta = min(tms) + r["anneal_offset"]
    ta = max(45.0, min(72.0, round(ta)))
    return ta, tms


def extension_seconds(longest_bp: int, polymerase: str = "Taq") -> int:
    r = RECIPES.get(polymerase, RECIPES["Taq"])
    secs = longest_bp / 1000.0 * r["rate_s_per_kb"]
    return int(max(15, math.ceil(secs / 5.0) * 5))


def program(primers: list[Primer], *, longest_bp: int = 1000, polymerase: str = "Taq",
            cycles: int = 30, annealing_c: float | None = None) -> tuple[list[Step], int, str]:
    """Steps, cycles and a sentence saying where the numbers come from."""
    r = RECIPES.get(polymerase, RECIPES["Taq"])
    ta, tms = suggested_annealing(primers, polymerase)
    if annealing_c is not None:
        ta = annealing_c
    ext = extension_seconds(longest_bp, polymerase)
    steps = [
        Step("Initial denaturation", *r["denature"]),
        Step("Denaturation", *r["cycle_denature"], in_cycle=True),
        Step("Annealing", ta, 30 if polymerase == "Taq" else 20, in_cycle=True),
        Step("Extension", 72 if polymerase != "Long-range" else 68, ext, in_cycle=True),
        Step("Final extension", *r["final"]),
        Step("Hold", 10, 0),
    ]
    why = []
    if tms:
        sign = "−" if r["anneal_offset"] < 0 else "+"
        why.append(f"Annealing {ta:g} °C = lower Tm {min(tms):.1f} °C {sign} "
                   f"{abs(r['anneal_offset']):g}")
        if len(tms) == 2 and abs(tms[0] - tms[1]) > 3:
            lo = round(min(tms) - 4)
            hi = round(max(tms))
            why.append(f"the Tm values differ by {abs(tms[0] - tms[1]):.1f} °C: try a gradient "
                       f"{lo} to {hi} °C")
    why.append(f"extension {fmt_seconds(ext)} for {longest_bp:,} bp "
               f"({r['rate_s_per_kb']} s per kb, {polymerase})")
    return steps, cycles, "; ".join(why) + "."


def fmt_seconds(s: int) -> str:
    if s <= 0:
        return "hold"
    if s < 60:
        return f"{s} s"
    m, sec = divmod(s, 60)
    return f"{m} min" if not sec else f"{m} min {sec} s"
