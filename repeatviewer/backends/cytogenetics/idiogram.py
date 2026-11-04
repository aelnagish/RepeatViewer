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
        for key in ("karyotype", "genome", "order", "lengths"):
            if isinstance(raw.get(key), str):
                setattr(s, key, raw[key])
        for key in ("min_share", "min_array_kb", "tolerance"):
            try:
                setattr(s, key, float(raw.get(key, getattr(s, key))))
            except (TypeError, ValueError):
                pass
        for key in ("homologues", "drag"):
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
    from ..sequence_features import fasta_record_lengths

    return {str(k): int(v) for k, v in fasta_record_lengths(genome).items()}


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
