"""Idiogram: FISH signals (in situ) against repeat positions in the assembly (in silico). No Qt.

In situ: every marked signal is put on its chromosome as a fraction of the length from the
p telomere (0) to the q telomere (1), using its arm and relative distance from the
centromere. Signals of the same probe at the same place are pooled into one band; the
band says on how many of the marked homologues it was seen.

In silico: hits of the probe on the assembly sequence linked to each chromosome are merged
into arrays; arrays shorter than ``min_array_kb`` are dropped (FISH cannot see them).
Positions are fractions of the sequence length; ``flip`` turns a sequence written q to p.

Comparison per chromosome: a FISH band and an assembly array agree when they are closer
than ``tolerance`` (fraction of the length).
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

from ..cyto_images import CytoImage, load_images, save_images, store_dir
from .importer import natural_key
from .karyotype import KaryoChromosome, Karyotype
from .model import Signal, image_details, set_image_details
from .signals import position_class

SETTINGS_FILE = "idiogram.json"
SOURCE_KINDS = ("annotations", "file", "none")
SOURCE_LABELS = {"annotations": "Saved annotations", "file": "Annotation file (GFF / BED)", "none": "None"}
ORDERS = ("karyotype", "length", "custom")
ORDER_LABELS = {"karyotype": "As in karyotype", "length": "Longest first", "custom": "Custom (drag a chromosome)"}
LENGTHS = ("karyotype", "assembly", "same")
LENGTH_LABELS = {"karyotype": "From karyotype", "assembly": "From assembly (Mb)", "same": "All the same"}
THRESHOLDS = ((0.5, "50 % of homologues"), (0.0, "any"))
STATUS_AGREES, STATUS_FISH, STATUS_ASM, STATUS_DIFFERS = "Agrees", "FISH only", "Assembly only", "Differs"
STATUS_NONE, STATUS_UNLINKED = "No signal", "No assembly"
# fraction of the arm from the centromere used when a signal has a position class but no measured distance
POSITION_DEFAULTS = {"pericentromeric": 0.05, "interstitial": 0.5, "subterminal": 0.87, "terminal": 0.97}
FISH_FASTA_SUFFIXES = (".fa", ".fasta", ".fna", ".fas")


# ---------------------------------------------------------------------------
# Settings (one file per folder)
# ---------------------------------------------------------------------------


@dataclass
class ProbeSource:
    kind: str = "annotations"  # see SOURCE_KINDS
    file: str = ""  # GFF / BED path for kind "file"
    match: str = ""  # text a feature name or type must contain ("" = the probe name)


PAPER_DEFAULTS = {"outlines": True, "formula": True, "types": True, "table": True, "scale": True,
                  "font": "Arial", "size": 8.0}
PAPER_WIDTH_PT = 510.0  # 180 mm, a two-column journal figure


@dataclass
class IdiogramSettings:
    karyotype: str = ""
    hidden_probes: list[str] = field(default_factory=list)
    min_share: float = 0.5
    genome: str = ""  # FASTA path (absolute, or a name in the folder)
    sources: dict[str, ProbeSource] = field(default_factory=dict)
    min_array_kb: float = 5.0
    tolerance: float = 0.1
    order: str = "karyotype"
    custom_order: list[str] = field(default_factory=list)
    lengths: str = "karyotype"
    homologues: bool = False
    drag: bool = True
    flips: dict[str, bool] = field(default_factory=dict)
    compare_with: str = ""  # second karyotype of the "Two karyotypes" view
    link_probes: bool = True
    paper: dict = field(default_factory=lambda: dict(PAPER_DEFAULTS))
    figure_style: str = "compare"  # how the idiogram is drawn in a figure: "compare" or "paper"

    def source(self, probe: str) -> ProbeSource:
        return self.sources.setdefault(probe, ProbeSource())

    def to_dict(self) -> dict:
        d = asdict(self)
        d["sources"] = {k: asdict(v) for k, v in self.sources.items()}
        return d

    @classmethod
    def from_dict(cls, raw) -> "IdiogramSettings":
        raw = raw if isinstance(raw, dict) else {}
        s = cls()
        for key in ("karyotype", "genome", "order", "lengths", "compare_with", "figure_style"):
            if isinstance(raw.get(key), str):
                setattr(s, key, raw[key])
        for key in ("min_share", "min_array_kb", "tolerance"):
            try:
                setattr(s, key, float(raw.get(key, getattr(s, key))))
            except (TypeError, ValueError):
                pass
        if isinstance(raw.get("paper"), dict):
            s.paper.update({k: v for k, v in raw["paper"].items() if k in PAPER_DEFAULTS})
        for key in ("homologues", "drag", "link_probes"):
            if key in raw:
                setattr(s, key, bool(raw[key]))
        s.hidden_probes = [str(x) for x in raw.get("hidden_probes") or []]
        s.custom_order = [str(x) for x in raw.get("custom_order") or []]
        s.flips = {str(k): bool(v) for k, v in (raw.get("flips") or {}).items()}
        for probe, src in (raw.get("sources") or {}).items():
            if isinstance(src, dict):
                kind = src.get("kind") if src.get("kind") in SOURCE_KINDS else "annotations"
                s.sources[str(probe)] = ProbeSource(kind, str(src.get("file") or ""), str(src.get("match") or ""))
        if s.order not in ORDERS:
            s.order = "karyotype"
        if s.lengths not in LENGTHS:
            s.lengths = "karyotype"
        return s


def settings_path(folder: Path) -> Path:
    return store_dir(Path(folder).expanduser().resolve()) / SETTINGS_FILE


def load_settings(folder: Path) -> IdiogramSettings:
    try:
        return IdiogramSettings.from_dict(json.loads(settings_path(folder).read_text(encoding="utf-8")))
    except (OSError, ValueError):
        return IdiogramSettings()


def save_settings(folder: Path, settings: IdiogramSettings) -> None:
    path = settings_path(folder)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(settings.to_dict(), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)


# ---------------------------------------------------------------------------
# In situ
# ---------------------------------------------------------------------------


def signal_fraction(signal: Signal, centromere_fraction: float) -> float | None:
    """Place of a signal from the p telomere (0) to the q telomere (1); None without an arm."""
    c = max(0.0, min(1.0, centromere_fraction))
    if signal.arm == "cen":
        return c
    if signal.arm not in ("p", "q"):
        return None
    rel = signal.rel_distance
    if rel is None:
        rel = POSITION_DEFAULTS.get(signal.position, 0.5)
    rel = max(0.0, min(1.0, float(rel)))
    return c - rel * c if signal.arm == "p" else c + rel * (1.0 - c)


def fraction_to_arm(fraction: float, centromere_fraction: float, cen_zone: float = 0.04) -> tuple[str, float]:
    """Inverse of ``signal_fraction``: (arm, relative distance from the centromere)."""
    c = max(0.0, min(1.0, centromere_fraction))
    f = max(0.0, min(1.0, fraction))
    if abs(f - c) <= cen_zone:
        return "cen", 0.0
    if f < c:
        return "p", (c - f) / c if c > 0 else 0.0
    return "q", (f - c) / (1.0 - c) if c < 1 else 0.0


@dataclass
class FishBand:
    chromosome: str
    probe: str
    arm: str
    position: str
    fraction: float
    seen: int  # homologues with this signal
    total: int  # homologues marked
    cells: int
    homologue: str = ""  # set when homologues are drawn separately
    refs: list[tuple[str, str]] = field(default_factory=list)  # (image id, signal id)

    @property
    def share(self) -> float:
        return self.seen / self.total if self.total else 0.0

    @property
    def label(self) -> str:
        return f"{self.seen}/{self.total}"


@dataclass
class FishData:
    bands: list[FishBand]
    cells: int  # cells with numbered chromosomes
    homologues: int  # marked homologues in those cells
    unplaced: int  # signals on a chromosome but without an arm (not measured)
    cells_per_probe: dict[str, int]


def _uses(image: CytoImage, karyotype: Karyotype) -> bool:
    d = image_details(image)
    if not d.chromosomes:
        return False
    if d.karyotype:
        return d.karyotype == karyotype.name
    names = {c.name for c in karyotype.chromosomes}
    return all(m.chromosome in names for m in d.chromosomes)


def collect_fish(folder: Path, karyotype: Karyotype, *, per_homologue: bool = False,
                 images: list[CytoImage] | None = None) -> FishData:
    images = [i for i in (images if images is not None else load_images(folder)) if _uses(i, karyotype)]
    cfrac = {c.name: c.centromere_fraction for c in karyotype.chromosomes}
    totals: dict[tuple[str, str], int] = {}
    groups: dict[tuple, dict] = {}
    unplaced = 0
    probe_cells: dict[str, set[str]] = {}
    homologues = 0
    for image in images:
        d = image_details(image)
        marks = {m.id: m for m in d.chromosomes}
        for m in marks.values():
            if m.chromosome not in cfrac:
                continue
            homologues += 1
            key = (m.chromosome, m.homologue if per_homologue else "")
            totals[key] = totals.get(key, 0) + 1
        for s in d.signals:
            mark = marks.get(s.chromosome_id)
            if mark is None or not s.probe or mark.chromosome not in cfrac:
                continue
            probe_cells.setdefault(s.probe, set()).add(image.id)
            f = signal_fraction(s, cfrac[mark.chromosome])
            if f is None:
                unplaced += 1
                continue
            hom = mark.homologue if per_homologue else ""
            position = s.position or position_class(s.arm, s.rel_distance)
            gkey = (mark.chromosome, hom, s.probe, s.arm, position)
            g = groups.setdefault(gkey, {"fractions": [], "seen": set(), "cells": set(), "refs": []})
            g["fractions"].append(f)
            g["seen"].add((image.id, mark.id))
            g["cells"].add(image.id)
            g["refs"].append((image.id, s.id))
    bands = []
    for (chrom, hom, probe, arm, position), g in groups.items():
        bands.append(FishBand(
            chromosome=chrom, probe=probe, arm=arm, position=position,
            fraction=sum(g["fractions"]) / len(g["fractions"]),
            seen=len(g["seen"]), total=totals.get((chrom, hom), len(g["seen"])), cells=len(g["cells"]),
            homologue=hom, refs=g["refs"],
        ))
    bands.sort(key=lambda b: (natural_key(b.chromosome), b.homologue, b.fraction, b.probe))
    return FishData(bands, len(images), homologues, unplaced, {p: len(c) for p, c in probe_cells.items()})


def visible_bands(data: FishData, hidden: set[str], min_share: float) -> list[FishBand]:
    return [b for b in data.bands if b.probe not in hidden and b.share + 1e-9 >= min_share]


def move_band(folder: Path, band: FishBand, fraction: float, centromere_fraction: float) -> int:
    """Correct the place of a pooled band: every signal in it gets the new arm and distance (set by hand)."""
    arm, rel = fraction_to_arm(fraction, centromere_fraction)
    wanted: dict[str, set[str]] = {}
    for image_id, signal_id in band.refs:
        wanted.setdefault(image_id, set()).add(signal_id)
    images = list(load_images(folder))
    changed = 0
    for image in images:
        ids = wanted.get(image.id)
        if not ids:
            continue
        d = image_details(image)
        for s in d.signals:
            if s.id in ids:
                s.arm, s.rel_distance, s.position, s.manual = arm, round(rel, 3), position_class(arm, rel), True
                changed += 1
        set_image_details(image, d)
    if changed:
        save_images(folder, images)
    band.arm, band.fraction = arm, fraction
    return changed


# ---------------------------------------------------------------------------
# In silico
# ---------------------------------------------------------------------------


@dataclass
class SilicoArray:
    chromosome: str
    probe: str
    start: float  # fraction of the sequence (after flip), start <= end
    end: float
    bp: int  # bases covered by hits
    hits: int

    @property
    def middle(self) -> float:
        return (self.start + self.end) / 2


def genome_candidates(folder: Path) -> list[Path]:
    folder = Path(folder).expanduser().resolve()
    try:
        return sorted((p for p in folder.iterdir() if p.is_file() and p.suffix.lower() in FISH_FASTA_SUFFIXES),
                      key=lambda p: natural_key(p.name))
    except OSError:
        return []


def resolve_genome(folder: Path, genome: str) -> Path | None:
    if not genome:
        return None
    path = Path(genome)
    if not path.is_absolute():
        path = Path(folder).expanduser().resolve() / genome
    return path if path.is_file() else None


def sequence_lengths(genome: Path) -> dict[str, int]:
    """Lengths from a .fai or sizes file beside the genome when there is one, else read once and saved as a sizes file."""
    from .karyotype import _newer, read_sizes, sizes_path_for, write_sizes

    side = sizes_path_for(genome)
    fresh = _newer(side, genome) or _newer(Path(str(genome) + ".fai"), genome)
    sizes = read_sizes(genome)
    if not fresh and sizes:
        try:
            write_sizes(side, sizes)
        except OSError:
            pass
    return sizes


def array_span_bp(array: "SilicoArray", length: int, flipped: bool) -> tuple[int, int]:
    """Start and end (1-based, inclusive) of an array on the assembly sequence, undoing a flip."""
    a, b = array.start, array.end
    if flipped:
        a, b = 1.0 - b, 1.0 - a
    return max(1, int(round(a * length)) + 1), max(1, int(round(b * length)))


def region_output_path(genome: Path, seq: str, start: int, end: int) -> Path:
    """'genome.fa' -> 'genome.fa_Chr2_191400001-193000000.fa' (full input name kept)."""
    safe = re.sub(r"[^A-Za-z0-9._-]+", "_", seq)
    return Path(genome).with_name(f"{Path(genome).name}_{safe}_{start}-{end}.fa")


def extract_region(genome: Path, seq: str, start: int, end: int, out: Path) -> Path:
    """Write bases ``start``..``end`` (1-based, inclusive) of one sequence as FASTA, reading the genome line by line."""
    parts: list[str] = []
    pos, inside = 0, False
    with open(genome, "r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith(">"):
                if inside:
                    break
                inside = line[1:].split(None, 1)[0] == seq if line[1:].strip() else False
                pos = 0
                continue
            if not inside:
                continue
            text = line.strip()
            lo, hi = pos + 1, pos + len(text)
            pos = hi
            if hi < start:
                continue
            parts.append(text[max(0, start - lo): end - lo + 1])
            if hi >= end:
                break
    bases = "".join(parts)
    if not bases:
        raise ValueError(f"{seq}:{start}-{end} not found in {Path(genome).name}.")
    body = "\n".join(bases[i:i + 60] for i in range(0, len(bases), 60))
    Path(out).write_text(f">{seq}:{start}-{end}\n{body}\n", encoding="utf-8")
    return Path(out)


def _matches(text: str, needle: str) -> bool:
    return needle.lower() in text.lower()


def features_for_probe(genome: Path, seq_names: list[str], probe: str, source: ProbeSource,
                       cache: dict | None = None) -> dict[str, list[tuple[int, int]]]:
    """(start, end) hits per sequence for one probe, 1-based inclusive."""
    needle = (source.match or probe).strip()
    if source.kind == "none" or not needle:
        return {}
    out: dict[str, list[tuple[int, int]]] = {}
    if source.kind == "annotations":
        key = ("annotations", str(genome), genome.stat().st_mtime if genome.exists() else 0)
        features = cache.get(key) if cache is not None else None
        if features is None:
            from ..sequence_features import load_features_by_id_from_json

            features = load_features_by_id_from_json(genome, list(seq_names))
            if cache is not None:
                cache[key] = features
        for seq, feats in features.items():
            for f in feats:
                if _matches(f.name, needle) or _matches(f.feature_type, needle):
                    out.setdefault(seq, []).append((int(f.start), int(f.end)))
    elif source.kind == "file" and source.file:
        path = Path(source.file)
        if not path.is_file():
            return {}
        key = ("file", str(path), path.stat().st_mtime)
        rows = cache.get(key) if cache is not None else None
        if rows is None:
            from ..gff_features import parse_annotation

            rows = parse_annotation(path)
            if cache is not None:
                cache[key] = rows
        for r in rows:
            text = " ".join([r.label(), r.feature_type] + [str(v) for v in r.attributes.values()])
            if _matches(text, needle):
                out.setdefault(r.seqid, []).append((int(r.start), int(r.end)))
    return out


def merge_hits(hits: list[tuple[int, int]], length: int, *, gap: int | None = None) -> list[tuple[int, int, int, int]]:
    """Merge hits closer than ``gap`` into arrays: (start, end, bases covered, number of hits)."""
    if not hits:
        return []
    gap = gap if gap is not None else max(50_000, int(length * 0.01))
    hits = sorted((min(a, b), max(a, b)) for a, b in hits)
    arrays = []
    s, e, covered, n = hits[0][0], hits[0][1], hits[0][1] - hits[0][0] + 1, 1
    for a, b in hits[1:]:
        if a - e <= gap:
            covered += max(0, b - max(a, e + 1) + 1)
            e = max(e, b)
            n += 1
        else:
            arrays.append((s, e, covered, n))
            s, e, covered, n = a, b, b - a + 1, 1
    arrays.append((s, e, covered, n))
    return arrays


def collect_silico(karyotype: Karyotype, lengths: dict[str, int], hits_by_probe: dict[str, dict[str, list[tuple[int, int]]]],
                   *, min_array_kb: float = 5.0, flips: dict[str, bool] | None = None) -> list[SilicoArray]:
    flips = flips or {}
    out = []
    for chrom in karyotype.chromosomes:
        seq = chrom.assembly_seq
        length = lengths.get(seq or "")
        if not seq or not length:
            continue
        for probe, by_seq in hits_by_probe.items():
            for start, end, covered, n in merge_hits(by_seq.get(seq, []), length):
                if covered < min_array_kb * 1000:
                    continue
                a, b = (start - 1) / length, end / length
                if flips.get(chrom.name):
                    a, b = 1.0 - b, 1.0 - a
                out.append(SilicoArray(chrom.name, probe, max(0.0, a), min(1.0, b), covered, n))
    return out


# ---------------------------------------------------------------------------
# Comparison and layout
# ---------------------------------------------------------------------------


@dataclass
class ChromosomeRow:
    chromosome: KaryoChromosome
    homologue: str  # "" when pooled
    sequence: str
    sequence_length: int
    fish: list[FishBand]
    silico: list[SilicoArray]
    status: str
    matched_fish: set[int] = field(default_factory=set)  # indices into fish
    matched_silico: set[int] = field(default_factory=set)

    @property
    def name(self) -> str:
        return f"{self.chromosome.name}{self.homologue}"


def _distance(band: FishBand, array: SilicoArray) -> float:
    if array.start <= band.fraction <= array.end:
        return 0.0
    return min(abs(band.fraction - array.start), abs(band.fraction - array.end))


def compare_row(row: ChromosomeRow, tolerance: float, linked: bool) -> None:
    row.matched_fish, row.matched_silico = set(), set()
    for i, band in enumerate(row.fish):
        for j, array in enumerate(row.silico):
            if array.probe == band.probe and _distance(band, array) <= tolerance:
                row.matched_fish.add(i)
                row.matched_silico.add(j)
    if not linked:
        row.status = STATUS_UNLINKED
        return
    fish_only = len(row.fish) - len(row.matched_fish)
    asm_only = len(row.silico) - len(row.matched_silico)
    if not row.fish and not row.silico:
        row.status = STATUS_NONE
    elif fish_only and asm_only:
        row.status = STATUS_DIFFERS
    elif fish_only:
        row.status = STATUS_FISH if not row.matched_fish else STATUS_DIFFERS
    elif asm_only:
        row.status = STATUS_ASM if not row.matched_silico else STATUS_DIFFERS
    else:
        row.status = STATUS_AGREES


def ordered_chromosomes(karyotype: Karyotype, settings: IdiogramSettings, lengths: dict[str, int]) -> list[KaryoChromosome]:
    chroms = list(karyotype.chromosomes)
    if settings.order == "length":
        if settings.lengths == "assembly":
            chroms.sort(key=lambda c: -lengths.get(c.assembly_seq or "", 0))
        else:
            chroms.sort(key=lambda c: -c.length)
    elif settings.order == "custom" and settings.custom_order:
        rank = {n: i for i, n in enumerate(settings.custom_order)}
        chroms.sort(key=lambda c: rank.get(c.name, len(rank) + karyotype.chromosomes.index(c)))
    return chroms


def build_rows(karyotype: Karyotype, settings: IdiogramSettings, fish: FishData, silico: list[SilicoArray],
               lengths: dict[str, int]) -> list[ChromosomeRow]:
    hidden = set(settings.hidden_probes)
    bands = visible_bands(fish, hidden, settings.min_share)
    arrays = [a for a in silico if a.probe not in hidden]
    rows = []
    for chrom in ordered_chromosomes(karyotype, settings, lengths):
        homs = [""]
        if settings.homologues:
            homs = [label[len(chrom.name):] for label in karyotype.homologue_labels(chrom)] or [""]
        seq = chrom.assembly_seq if chrom.assembly_seq in lengths else ""
        for hom in homs:
            row = ChromosomeRow(
                chromosome=chrom, homologue=hom, sequence=chrom.assembly_seq, sequence_length=lengths.get(seq, 0),
                fish=[b for b in bands if b.chromosome == chrom.name and b.homologue == hom],
                silico=[a for a in arrays if a.chromosome == chrom.name], status="",
            )
            compare_row(row, settings.tolerance, bool(seq))
            rows.append(row)
    return rows


def relative_heights(rows: list[ChromosomeRow], mode: str) -> list[float]:
    """Height of each row from 0 to 1 (1 = the longest)."""
    if mode == "same":
        return [1.0] * len(rows)
    if mode == "assembly":
        values = [float(r.sequence_length) for r in rows]
    else:
        values = [float(r.chromosome.length) for r in rows]
    top = max(values, default=0.0)
    if top <= 0:
        return [1.0] * len(rows)
    return [max(0.15, v / top) if v > 0 else 0.5 for v in values]


def comparison_tsv(rows: list[ChromosomeRow]) -> str:
    lines = ["chromosome\tsequence\tstatus\tsource\tprobe\tarm\tposition\tfrom\tto\tseen\tmarked\tbases\tmatched"]
    for r in rows:
        if not r.fish and not r.silico:
            lines.append(f"{r.name}\t{r.sequence}\t{r.status}\t\t\t\t\t\t\t\t\t\t")
        for i, b in enumerate(r.fish):
            lines.append(f"{r.name}\t{r.sequence}\t{r.status}\tin situ\t{b.probe}\t{b.arm}\t{b.position}\t{b.fraction:.3f}\t"
                         f"{b.fraction:.3f}\t{b.seen}\t{b.total}\t\t{'yes' if i in r.matched_fish else 'no'}")
        for j, a in enumerate(r.silico):
            lines.append(f"{r.name}\t{r.sequence}\t{r.status}\tin silico\t{a.probe}\t\t\t{a.start:.3f}\t{a.end:.3f}\t\t\t"
                         f"{a.bp}\t{'yes' if j in r.matched_silico else 'no'}")
    return "\n".join(lines) + "\n"


def title_text(karyotype: Karyotype, probes: list[str]) -> tuple[str, str]:
    """("CroSat2 and 45S rDNA on the", "C. sativus") plus " karyotype" is added by the view."""
    if not probes:
        return "Signals on the", karyotype.species or karyotype.name
    head = probes[0] if len(probes) == 1 else ", ".join(probes[:-1]) + " and " + probes[-1]
    return f"{head} on the", karyotype.species or karyotype.name


# ---------------------------------------------------------------------------
# Everything the Idiogram view (and the Figure) needs, in one call
# ---------------------------------------------------------------------------

FALLBACK_COLOURS = ("#2FB457", "#E8405F", "#E8A33D", "#5B74FF", "#B05BD6", "#2BB5C4")


@dataclass
class IdiogramView:
    karyotype: Karyotype | None
    settings: IdiogramSettings
    fish: FishData | None
    rows: list[ChromosomeRow]
    heights: list[float]
    colours: dict[str, str]
    probes: list[str]  # shown probes
    lengths: dict[str, int]
    genome: Path | None
    pooled: list[FishBand] = field(default_factory=list)  # signals over all homologues (paper style)


def all_probes(folder: Path, fish: FishData | None) -> list[str]:
    probes: list[str] = sorted(fish.cells_per_probe) if fish else []
    for image in load_images(folder):
        for c in image_details(image).channels:
            if c.probe and c.probe not in probes:
                probes.append(c.probe)
    return probes


def probe_colours(folder: Path, probes: list[str]) -> dict[str, str]:
    """Colour of each probe: its channel colour (not white or the DAPI blue), else a fixed palette."""
    colours: dict[str, str] = {}
    for image in load_images(folder):
        for c in image_details(image).channels:
            if c.probe and c.probe not in colours and c.display.color.lower() not in ("#ffffff", "#5b74ff"):
                colours[c.probe] = c.display.color
    for n, p in enumerate(probes):
        colours.setdefault(p, FALLBACK_COLOURS[n % len(FALLBACK_COLOURS)])
    return colours


def compute_view(folder: Path, settings: IdiogramSettings | None = None, cache: dict | None = None) -> IdiogramView:
    from .karyotype import list_karyotypes

    folder = Path(folder).expanduser().resolve()
    settings = settings or load_settings(folder)
    cache = cache if cache is not None else {}
    karyotypes = list_karyotypes(folder)
    k = next((x for x in karyotypes if x.name == settings.karyotype), None) or (karyotypes[0] if karyotypes else None)
    if k is None:
        return IdiogramView(None, settings, None, [], [], {}, [], {}, None)
    fish = collect_fish(folder, k, per_homologue=settings.homologues)
    genome = resolve_genome(folder, settings.genome)
    lengths: dict[str, int] = {}
    if genome is not None:
        key = ("lengths", str(genome), genome.stat().st_mtime)
        if key not in cache:
            try:
                cache[key] = sequence_lengths(genome)
            except Exception:  # noqa: BLE001 - unreadable FASTA: no in silico bar
                cache[key] = {}
        lengths = cache[key]
    every = all_probes(folder, fish)
    probes = [p for p in every if p not in settings.hidden_probes]
    hits: dict[str, dict] = {}
    if genome is not None and lengths:
        for probe in probes:
            try:
                hits[probe] = features_for_probe(genome, list(lengths), probe, settings.source(probe), cache)
            except Exception:  # noqa: BLE001 - a broken annotation file must not stop the view
                hits[probe] = {}
    silico = collect_silico(k, lengths, hits, min_array_kb=settings.min_array_kb, flips=settings.flips)
    rows = build_rows(k, settings, fish, silico, lengths)
    pooled = fish if not settings.homologues else collect_fish(folder, k)
    return IdiogramView(k, settings, fish, rows, relative_heights(rows, settings.lengths),
                        probe_colours(folder, every), probes, lengths, genome,
                        visible_bands(pooled, set(settings.hidden_probes), settings.min_share))


# ---------------------------------------------------------------------------
# Paper style and two karyotypes
# ---------------------------------------------------------------------------

POSITION_SHORT = {"terminal": "ter", "subterminal": "subter", "interstitial": "int", "pericentromeric": "peri",
                  "centromeric": "cen"}


def nice_scale(longest: float) -> float:
    """A round scale bar length near a third of the longest chromosome (1, 2, 5, 10, 20 …)."""
    if longest <= 0:
        return 0.0
    target = longest / 3.0
    step = 10 ** int(f"{target:e}".split("e")[1])
    for m in (1, 2, 5, 10):
        if m * step >= target * 0.7:
            return float(m * step)
    return float(10 * step)


def band_text(band: FishBand) -> str:
    """'45S rDNA q ter (21/21)'."""
    where = POSITION_SHORT.get(band.position, band.position)
    arm = f"{band.arm} " if band.arm in ("p", "q") else ""
    return f"{band.probe} {arm}{where} ({band.label})".replace("  ", " ")


def paper_table(karyotype: Karyotype, bands: list[FishBand]) -> list[list[str]]:
    """Rows Chr, Length, CI, Type, Signals for chromosomes with a signal or a band."""
    rows = []
    for c in karyotype.chromosomes:
        texts = [band_text(b) for b in bands if b.chromosome == c.name]
        for b in c.bands:
            mid = (b.start + b.end) / 2
            arm = "p" if mid < c.centromere_fraction else "q"
            probe = b.kind[len("FISH "):] if b.kind.startswith("FISH ") else ""
            name = probe or {"Secondary constriction (NOR)": "NOR", "DAPI+ heterochromatin": "DAPI+",
                             "CMA+ band": "CMA+"}.get(b.kind, b.kind)
            texts.append(f"{name} {arm}")
        if texts:
            length = f"{c.length:g} {karyotype.unit}" if c.length else ""
            rows.append([c.name, length, f"{c.centromere_index:.0f}" if c.length else "", c.type, ", ".join(texts)])
    return rows


@dataclass
class PairLink:
    probe: str
    index: int  # chromosome position in the top karyotype
    top: float  # fraction on the top chromosome
    bottom: float
    bottom_index: int = -1  # chromosome position in the bottom karyotype


def pair_compare(top: Karyotype, top_bands: list[FishBand], bottom: Karyotype,
                 bottom_bands: list[FishBand]) -> tuple[list[PairLink], list[tuple[int, FishBand]], list[tuple[int, FishBand]]]:
    """Link the same probe on the same chromosome (by name, else by position). Returns links and the unmatched bands."""
    bottom_index = {c.name: i for i, c in enumerate(bottom.chromosomes)}
    links: list[PairLink] = []
    only_top: list[tuple[int, FishBand]] = []
    used: set[int] = set()
    for i, c in enumerate(top.chromosomes):
        j = bottom_index.get(c.name, i if i < len(bottom.chromosomes) else -1)
        other = bottom.chromosomes[j].name if j >= 0 else None
        for b in (x for x in top_bands if x.chromosome == c.name):
            match = next((n for n, x in enumerate(bottom_bands)
                          if n not in used and x.chromosome == other and x.probe == b.probe), None)
            if match is None:
                only_top.append((i, b))
            else:
                used.add(match)
                links.append(PairLink(b.probe, i, b.fraction, bottom_bands[match].fraction, j))
    names = {c.name: i for i, c in enumerate(bottom.chromosomes)}
    only_bottom = [(names.get(b.chromosome, 0), b) for n, b in enumerate(bottom_bands) if n not in used]
    return links, only_top, only_bottom
