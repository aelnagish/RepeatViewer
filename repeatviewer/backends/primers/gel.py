"""Virtual agarose gel: where bands of a given size run.

Migration is linear in log(size) between the ladder's smallest and largest
band, which is close to what a real gel shows over a 2–3 log range. The same
geometry is used by the Qt gel view and by the SVG in the HTML report, so
both always agree.
"""

from __future__ import annotations

import math
from dataclasses import asdict, dataclass, field
from html import escape

LADDERS: dict[str, tuple[int, ...]] = {
    "100 bp": (100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1200, 1500, 2000, 3000),
    "1 kb Plus": (75, 200, 300, 400, 500, 700, 1000, 1500, 2000, 3000, 4000, 5000, 7000,
                  10000, 20000),
    "50 bp": (50, 100, 150, 200, 250, 300, 350, 400, 450, 500, 600, 700, 800, 900, 1000),
    "25 bp (small)": (25, 50, 75, 100, 125, 150, 175, 200, 250, 300, 350, 400, 450, 500),
}
DEFAULT_LADDER = "1 kb Plus"
# Bands heavier on the ladder (as most commercial ladders print them).
_BRIGHT = {500, 1000, 3000}


def pick_ladder(sizes: list[int]) -> str:
    """Ladder whose range suits the products best."""
    if not sizes:
        return DEFAULT_LADDER
    lo, hi = min(sizes), max(sizes)
    if hi <= 450 and lo >= 20:
        return "25 bp (small)" if hi < 200 else "50 bp"
    if hi <= 2500 and lo >= 80:
        return "100 bp"
    return DEFAULT_LADDER


@dataclass
class Band:
    size: int
    label: str = ""
    weight: float = 1.0     # relative brightness 0–1
    highlight: bool = False  # the expected (on-target) product

    def to_dict(self) -> dict:
        return asdict(self)


@dataclass
class Lane:
    label: str
    bands: list[Band] = field(default_factory=list)
    note: str = ""

    def to_dict(self) -> dict:
        return {"label": self.label, "note": self.note, "bands": [b.to_dict() for b in self.bands]}

    @classmethod
    def from_dict(cls, d: dict) -> Lane:
        return cls(label=d.get("label", ""), note=d.get("note", ""),
                   bands=[Band(**{k: b[k] for k in ("size", "label", "weight", "highlight") if k in b})
                          for b in d.get("bands", [])])


def size_range(ladder: str, lanes: list[Lane]) -> tuple[float, float]:
    sizes = list(LADDERS.get(ladder, LADDERS[DEFAULT_LADDER]))
    for lane in lanes:
        sizes.extend(b.size for b in lane.bands if b.size > 0)
    lo = max(10.0, min(sizes) * 0.8)
    hi = max(sizes) * 1.15
    return lo, hi


def band_y(size: float, lo: float, hi: float, top: float, bottom: float) -> float:
    """Vertical position of a band: big near the wells (top), small at the bottom."""
    size = min(max(size, lo), hi)
    t = (math.log(hi) - math.log(size)) / max(1e-9, math.log(hi) - math.log(lo))
    return top + t * (bottom - top)


def merge_bands(sizes: list[int], *, resolution: float = 0.03, max_bands: int = 40) -> list[Band]:
    """Products whose sizes differ by < ``resolution`` (fraction) run as one band.

    A tandem array gives many products of each ladder size; the band weight
    grows with their number, so the gel looks like a real ladder PCR.
    """
    if not sizes:
        return []
    sizes = sorted(sizes)
    groups: list[list[int]] = [[sizes[0]]]
    for s in sizes[1:]:
        if (s - groups[-1][0]) / max(1, groups[-1][0]) <= resolution:
            groups[-1].append(s)
        else:
            groups.append([s])
    top = max(len(g) for g in groups)
    bands = []
    for g in groups:
        size = int(round(sum(g) / len(g)))
        label = f"{size:,} bp" + (f" ×{len(g)}" if len(g) > 1 else "")
        bands.append(Band(size=size, label=label, weight=0.35 + 0.65 * len(g) / top))
    bands.sort(key=lambda b: -b.weight)
    bands = sorted(bands[:max_bands], key=lambda b: -b.size)
    return bands


def gel_svg(lanes: list[Lane], *, ladder: str = DEFAULT_LADDER, width: int | None = None,
            height: int = 360, title: str = "") -> str:
    """Self-contained SVG of the gel (for the HTML report and Save image)."""
    lane_w, gap, left, right = 64, 14, 58, 16
    n = len(lanes) + 1
    w = width or (left + right + n * lane_w + (n - 1) * gap)
    top, bottom = 52, height - 18
    lo, hi = size_range(ladder, lanes)
    out = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{height}" '
        f'viewBox="0 0 {w} {height}" font-family="Helvetica, Arial, sans-serif">',
        f'<rect x="0" y="0" width="{w}" height="{height}" rx="8" fill="#16181d"/>',
    ]
    if title:
        out.append(f'<text x="{w / 2:.0f}" y="16" fill="#c9cdd6" font-size="11" '
                   f'text-anchor="middle">{escape(title)}</text>')

    def lane_x(i: int) -> float:
        return left + i * (lane_w + gap)

    labels = [ladder] + [lane.label for lane in lanes]
    for i, lab in enumerate(labels):
        x = lane_x(i)
        out.append(f'<rect x="{x:.0f}" y="{top - 14}" width="{lane_w}" height="6" rx="1.5" '
                   f'fill="#2c3039"/>')
        short = lab if len(lab) <= 12 else lab[:11] + "…"
        out.append(f'<text x="{x + lane_w / 2:.0f}" y="{top - 22}" fill="#c9cdd6" '
                   f'font-size="10" text-anchor="middle"><title>{escape(lab)}</title>'
                   f'{escape(short)}</text>')
    for size in LADDERS.get(ladder, LADDERS[DEFAULT_LADDER]):
        y = band_y(size, lo, hi, top, bottom)
        alpha = 0.95 if size in _BRIGHT else 0.6
        out.append(f'<rect x="{lane_x(0) + 6:.0f}" y="{y - 1.5:.1f}" width="{lane_w - 12}" '
                   f'height="3" rx="1.5" fill="#f2f4f8" fill-opacity="{alpha}"/>')
        out.append(f'<text x="{lane_x(0) - 6:.0f}" y="{y + 3:.1f}" fill="#8b919c" '
                   f'font-size="9" text-anchor="end">{size:,}</text>')
    for i, lane in enumerate(lanes, start=1):
        x = lane_x(i)
        for b in lane.bands:
            y = band_y(b.size, lo, hi, top, bottom)
            color = "#7dd3a0" if b.highlight else "#f2f4f8"
            opacity = 0.3 + 0.7 * max(0.0, min(1.0, b.weight))
            out.append(f'<rect x="{x + 6:.0f}" y="{y - 2:.1f}" width="{lane_w - 12}" '
                       f'height="4" rx="2" fill="{color}" fill-opacity="{opacity:.2f}">'
                       f'<title>{escape(b.label or f"{b.size:,} bp")}</title></rect>')
    out.append("</svg>")
    return "\n".join(out)
