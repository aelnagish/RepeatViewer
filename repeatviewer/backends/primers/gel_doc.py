"""In-silico agarose gel document (``.rvgel``). Qt-free.

A gel holds lanes of bands (size, relative amount, label, optional
sequence so a lane can be cut with a restriction enzyme). The migration model
puts bands on a log-size scale whose resolving range follows the agarose
percentage (Sambrook-style table, interpolated), with sizes outside that range
squeezed towards the wells or the front, so a 0.7 % gel separates 10 kb bands
and a 3 % gel separates 50 bp steps. Brightness follows the mass of DNA
(amount × length) unless switched off: short products look fainter, as with
SYBR / GelRed staining.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass, field
from pathlib import Path

GEL_SUFFIX = ".rvgel"
GEL_FORMAT = "repeatviewer-gel-1"

# (size bp, relative mass) ; heavier reference bands as on the product sheets.
LADDERS: dict[str, list[tuple[int, float]]] = {
    "GeneRuler 1 kb Plus": [(20000, 1), (10000, 1), (7000, 1), (5000, 3), (4000, 1), (3000, 1),
                            (2000, 1), (1500, 3), (1000, 1), (700, 1), (500, 3), (400, 1),
                            (300, 1), (200, 1), (75, 1)],
    "GeneRuler 1 kb": [(10000, 1), (8000, 1), (6000, 1), (5000, 1), (4000, 1), (3500, 1),
                       (3000, 3), (2500, 1), (2000, 1), (1500, 1), (1000, 3), (750, 1),
                       (500, 1), (250, 1)],
    "GeneRuler 100 bp Plus": [(3000, 1), (2000, 1), (1500, 1), (1200, 1), (1000, 1), (900, 1),
                              (800, 1), (700, 1), (600, 1), (500, 3), (400, 1), (300, 1),
                              (200, 1), (100, 1)],
    "GeneRuler 50 bp": [(1000, 1), (900, 1), (800, 1), (700, 1), (600, 1), (500, 1), (400, 1),
                        (300, 1), (250, 3), (200, 1), (150, 1), (100, 1), (50, 1)],
    "Lambda DNA / HindIII": [(23130, 1), (9416, 1), (6557, 1), (4361, 1), (2322, 1),
                             (2027, 1), (564, 1), (125, 1)],
}
DEFAULT_LADDER = "GeneRuler 1 kb Plus"

# Agarose % -> (smallest, largest) size separated well (bp).
_RANGES = [(0.5, 1000, 30000), (0.7, 800, 12000), (1.0, 500, 10000), (1.2, 400, 7000),
           (1.5, 200, 3000), (2.0, 100, 2000), (3.0, 50, 1000), (4.0, 20, 500)]

RUNS = {"Short (30 min)": 0.72, "Standard (45 min)": 0.88, "Long (75 min)": 1.0}
DEFAULT_RUN = "Standard (45 min)"
LOOKS = ("Dark (SYBR / GelRed)", "Light (inverted, for print)")

LANE_LADDER = "ladder"
LANE_PCR = "pcr"
LANE_CUSTOM = "custom"


def resolving_range(agarose: float) -> tuple[float, float]:
    a = min(max(float(agarose), _RANGES[0][0]), _RANGES[-1][0])
    for (p0, lo0, hi0), (p1, lo1, hi1) in zip(_RANGES, _RANGES[1:]):
        if p0 <= a <= p1:
            t = (a - p0) / (p1 - p0) if p1 > p0 else 0.0
            lo = math.exp(math.log(lo0) + t * (math.log(lo1) - math.log(lo0)))
            hi = math.exp(math.log(hi0) + t * (math.log(hi1) - math.log(hi0)))
            return lo, hi
    return float(_RANGES[-1][1]), float(_RANGES[-1][2])


def migration(size: float, agarose: float = 1.0, run: float = RUNS[DEFAULT_RUN]) -> float:
    """Distance from the wells as a fraction of the gel (0 = well; > 1 ran off)."""
    lo, hi = resolving_range(agarose)
    mid = (math.log(lo) + math.log(hi)) / 2
    k = 2 * math.log(9) / (math.log(hi) - math.log(lo))
    x = math.log(max(1.0, float(size)))
    y = 1.0 / (1.0 + math.exp(k * (x - mid)))
    return y * run * 1.02


@dataclass
class GelBand:
    size: int
    amount: float = 1.0           # relative molar amount
    label: str = ""
    sequence: str = ""            # kept for digests (empty: not cuttable)
    ref: str = ""                 # e.g. product name in a PCR products FASTA

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, d: dict) -> GelBand:
        return cls(**{k: d[k] for k in cls.__dataclass_fields__ if k in d})


@dataclass
class GelLane:
    label: str
    kind: str = LANE_CUSTOM
    bands: list[GelBand] = field(default_factory=list)
    cut_with: str = ""            # enzyme name, "" = uncut
    source: str = ""              # file the lane came from
    note: str = ""

    def to_dict(self) -> dict:
        d = asdict(self)
        d["bands"] = [b.to_dict() for b in self.bands]
        return d

    @classmethod
    def from_dict(cls, d: dict) -> GelLane:
        return cls(label=d.get("label", ""), kind=d.get("kind", LANE_CUSTOM),
                   bands=[GelBand.from_dict(b) for b in d.get("bands") or []],
                   cut_with=d.get("cut_with", ""), source=d.get("source", ""),
                   note=d.get("note", ""))


@dataclass
class GelDoc:
    ladder: str = DEFAULT_LADDER
    agarose: float = 1.5
    run: str = DEFAULT_RUN
    look: str = LOOKS[0]
    mass_brightness: bool = True
    show_dimers: bool = True
    lanes: list[GelLane] = field(default_factory=list)
    title: str = ""

    def to_dict(self) -> dict:
        return {"format": GEL_FORMAT, "ladder": self.ladder, "agarose": self.agarose,
                "run": self.run, "look": self.look, "mass_brightness": self.mass_brightness,
                "show_dimers": self.show_dimers, "title": self.title,
                "lanes": [ln.to_dict() for ln in self.lanes]}

    @classmethod
    def from_dict(cls, d: dict) -> GelDoc:
        if d.get("format") != GEL_FORMAT:
            raise ValueError("Not a RepeatViewer gel file.")
        doc = cls(ladder=d.get("ladder", DEFAULT_LADDER), agarose=float(d.get("agarose", 1.5)),
                  run=d.get("run", DEFAULT_RUN), look=d.get("look", LOOKS[0]),
                  mass_brightness=bool(d.get("mass_brightness", True)),
                  show_dimers=bool(d.get("show_dimers", True)), title=d.get("title", ""))
        doc.lanes = [GelLane.from_dict(x) for x in d.get("lanes") or []]
        if doc.ladder not in LADDERS:
            doc.ladder = DEFAULT_LADDER
        if doc.run not in RUNS:
            doc.run = DEFAULT_RUN
        return doc

    def save(self, path: Path) -> Path:
        path = Path(path)
        path.write_text(json.dumps(self.to_dict(), indent=1) + "\n", encoding="utf-8")
        return path

    @classmethod
    def load(cls, path: Path) -> GelDoc:
        return cls.from_dict(json.loads(Path(path).read_text(encoding="utf-8")))

    @property
    def run_factor(self) -> float:
        return RUNS.get(self.run, RUNS[DEFAULT_RUN])


def ladder_lane(name: str = DEFAULT_LADDER) -> GelLane:
    bands = [GelBand(size=s, amount=m / max(1, s) * 1000.0, label=f"{s:,}")
             for s, m in LADDERS.get(name, LADDERS[DEFAULT_LADDER])]
    return GelLane(label="L", kind=LANE_LADDER, bands=bands, note=name)


def set_ladder(doc: GelDoc, name: str) -> None:
    doc.ladder = name if name in LADDERS else DEFAULT_LADDER
    for i, lane in enumerate(doc.lanes):
        if lane.kind == LANE_LADDER:
            doc.lanes[i] = ladder_lane(doc.ladder)


# --------------------------------------------------------------------------
# Digest
# --------------------------------------------------------------------------


def digest_sizes(sequence: str, enzyme: str) -> list[int]:
    """Fragment lengths of a linear sequence cut with ``enzyme`` (uncut: [len])."""
    if not sequence:
        return []
    if not enzyme:
        return [len(sequence)]
    from Bio import Restriction
    from Bio.Seq import Seq

    enz = getattr(Restriction, enzyme, None)
    if enz is None:
        return [len(sequence)]
    frags = enz.catalyse(Seq(sequence.upper()), linear=True)
    return [len(f) for f in frags]


def lane_bands(lane: GelLane, *, show_dimers: bool = True) -> list[GelBand]:
    """Bands as they run: digested when the lane is cut, dimers optional."""
    out = []
    for b in lane.bands:
        if not show_dimers and "dimer" in b.label.lower():
            continue
        if lane.cut_with and b.sequence:
            for n in digest_sizes(b.sequence, lane.cut_with):
                out.append(GelBand(size=n, amount=b.amount, label=f"{b.label} / {lane.cut_with}"
                                   if b.label else lane.cut_with, ref=b.ref))
        else:
            out.append(b)
    return out


# --------------------------------------------------------------------------
# Layout (shared by the Qt view and the SVG export)
# --------------------------------------------------------------------------


@dataclass
class DrawnBand:
    lane: int
    band: GelBand
    y: float          # fraction of the gel height from the wells
    intensity: float  # 0..1
    thickness: float  # fraction of a normal band


def layout(doc: GelDoc) -> list[DrawnBand]:
    drawn: list[tuple[int, GelBand, float, float]] = []
    for i, lane in enumerate(doc.lanes):
        for b in lane_bands(lane, show_dimers=doc.show_dimers):
            y = migration(b.size, doc.agarose, doc.run_factor)
            if y > 1.0:
                continue
            mass = b.amount * (b.size if doc.mass_brightness else 1.0)
            drawn.append((i, b, y, mass))
    if not drawn:
        return []
    # Normalise ladders and samples separately so a faint sample is still visible.
    top_lad = max((m for i, _b, _y, m in drawn if doc.lanes[i].kind == LANE_LADDER), default=0) or 1
    top_smp = max((m for i, _b, _y, m in drawn if doc.lanes[i].kind != LANE_LADDER), default=0) or 1
    out = []
    for i, b, y, m in drawn:
        top = top_lad if doc.lanes[i].kind == LANE_LADDER else top_smp
        rel = max(0.0, min(1.0, m / top))
        inten = 0.12 + 0.88 * rel ** 0.6 if rel > 0 else 0.0
        thick = 0.8 + 0.7 * rel
        if "dimer" in b.label.lower():
            thick = 2.2
            inten = min(inten, 0.3)
        out.append(DrawnBand(lane=i, band=b, y=y, intensity=inten, thickness=thick))
    return out


def ladder_marks(doc: GelDoc) -> list[tuple[int, float, bool]]:
    """(size, y, heavy) of the ladder bands for the size axis."""
    out = []
    for s, m in LADDERS.get(doc.ladder, []):
        y = migration(s, doc.agarose, doc.run_factor)
        if y <= 1.0:
            out.append((s, y, m > 1))
    return out


def size_at(doc: GelDoc, y: float) -> float:
    """Inverse of :func:`migration` (for the photo comparison and hover)."""
    lo, hi = resolving_range(doc.agarose)
    mid = (math.log(lo) + math.log(hi)) / 2
    k = 2 * math.log(9) / (math.log(hi) - math.log(lo))
    f = y / (doc.run_factor * 1.02)
    f = min(max(f, 1e-6), 1 - 1e-6)
    return math.exp(mid + math.log(1 / f - 1) / k)


def to_svg(doc: GelDoc, *, lane_w: int = 56, gap: int = 20, height: int = 600) -> str:
    """Stand-alone SVG of the gel (sizes on the left, lane labels below)."""
    n = max(1, len(doc.lanes))
    left = 70
    top = 30
    w = left + n * (lane_w + gap) + gap
    h = top + height + 60
    dark = doc.look == LOOKS[0]
    bg = "#151a22" if dark else "#f4f4f2"
    band_rgb = (232, 241, 255) if dark else (20, 20, 20)
    parts = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" '
             f'viewBox="0 0 {w} {h}" font-family="Helvetica, Arial, sans-serif">',
             '<defs><filter id="b" x="-10%" y="-300%" width="120%" height="700%">'
             '<feGaussianBlur stdDeviation="1.4"/></filter></defs>',
             f'<rect x="{left - 6}" y="{top - 16}" width="{w - left}" height="{height + 24}" '
             f'rx="6" fill="{bg}"/>']
    for i in range(n):
        x = left + gap + i * (lane_w + gap)
        parts.append(f'<rect x="{x}" y="{top - 8}" width="{lane_w}" height="5" '
                     f'fill="{"#39414d" if dark else "#c9c9c4"}"/>')
    for s, y, heavy in ladder_marks(doc):
        yy = top + y * height
        weight = ' font-weight="bold"' if heavy else ""
        parts.append(f'<text x="{left - 12}" y="{yy + 4:.1f}" font-size="11" text-anchor="end" '
                     f'fill="#5f6368"{weight}>{s:,}</text>')
    parts.append('<g filter="url(#b)">')
    for d in layout(doc):
        x = left + gap + d.lane * (lane_w + gap)
        yy = top + d.y * height
        t = 3.0 * d.thickness
        r, g, b = band_rgb
        parts.append(f'<rect x="{x}" y="{yy - t / 2:.1f}" width="{lane_w}" height="{t:.1f}" '
                     f'fill="rgb({r},{g},{b})" opacity="{d.intensity:.2f}"/>')
    parts.append("</g>")
    for i, lane in enumerate(doc.lanes):
        x = left + gap + i * (lane_w + gap) + lane_w / 2
        label = lane.label.replace("&", "&amp;").replace("<", "&lt;")
        num = "L" if lane.kind == LANE_LADDER else str(i + 1)
        parts.append(f'<text x="{x:.0f}" y="{top + height + 26}" font-size="12" '
                     f'text-anchor="middle" fill="#1d1d1f" font-weight="bold">{num}</text>')
        parts.append(f'<text x="{x:.0f}" y="{top + height + 42}" font-size="10" '
                     f'text-anchor="middle" fill="#5f6368">{label[:14]}</text>')
    parts.append("</svg>")
    return "\n".join(parts)


# --------------------------------------------------------------------------
# From a PCR result
# --------------------------------------------------------------------------


def lanes_from_pcr(result, *, min_yield: float = 0.01, max_bands: int = 60) -> list[GelLane]:
    """One lane per PCR lane (template / methylation state)."""
    lanes = []
    for lane in result.lanes:
        bands = []
        for p in result.lane_products(lane):
            if p.yield_ < min_yield and p.kind != "dimer":
                continue
            label = p.label or p.kind
            if p.kind == "off-target":
                label = "off-target"
            bands.append(GelBand(size=p.size, amount=max(p.yield_, 0.001), label=label,
                                 sequence=p.sequence if len(p.sequence) <= 50000 else "",
                                 ref=f"{p.forward}+{p.reverse} {p.size} bp"))
            if len(bands) >= max_bands:
                break
        lanes.append(GelLane(label=lane, kind=LANE_PCR, bands=bands,
                             source=result.templates[0] if result.templates else ""))
    return lanes


def gel_from_pcr(result, *, ladder: str = "") -> GelDoc:
    sizes = [p.size for p in result.products if p.kind != "dimer"] or [500]
    if not ladder:
        ladder = pick_ladder_name(sizes)
    doc = GelDoc(ladder=ladder, agarose=pick_agarose(sizes))
    doc.lanes = [ladder_lane(ladder)] + lanes_from_pcr(result)
    return doc


def pick_ladder_name(sizes: list[int]) -> str:
    hi = max(sizes) if sizes else 1000
    lo = min(sizes) if sizes else 100
    if hi <= 900 and lo >= 40:
        return "GeneRuler 100 bp Plus" if lo >= 90 else "GeneRuler 50 bp"
    return "GeneRuler 1 kb Plus"


def pick_agarose(sizes: list[int]) -> float:
    if not sizes:
        return 1.0
    mid = sorted(sizes)[len(sizes) // 2]
    if mid < 200:
        return 3.0
    if mid < 600:
        return 2.0
    if mid < 1500:
        return 1.5
    if mid < 5000:
        return 1.0
    return 0.7


def unique_gel_path(folder: Path, base: str) -> Path:
    cand = folder / f"{base}{GEL_SUFFIX}"
    n = 2
    while cand.exists():
        cand = folder / f"{base}_{n}{GEL_SUFFIX}"
        n += 1
    return cand
