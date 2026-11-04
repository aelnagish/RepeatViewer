"""Marking one cell: number chromosomes, place signals, measure axes. No Qt.

``SignalsSession`` edits the chromosome marks and signals of one image against
a karyotype. Saving only replaces those parts of the stored record, so display
changes made in the Image tab are never overwritten.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from ..cyto_images import CytoImage, load_images, update_image
from . import channels as ch
from . import signals as sg
from .image_session import ImageSession
from .importer import channel_layers, natural_key
from .karyotype import Karyotype, list_karyotypes
from .lists import dye_color, load_lists
from .model import POSITIONS, STRENGTHS, ChromosomeMark, ImageDetails, Signal, image_details, load_slides, new_id, set_image_details

ARM_LABELS = {"cen": "Centromere", "p": "p", "q": "q"}
POSITION_CHOICES = POSITIONS
TOOLS = (("select", "Select"), ("number", "Number chromosomes"), ("signal", "Add signal"), ("measure", "Measure"))
TOOL_HINTS = {
    "select": "Click a chromosome or a signal to edit it",
    "number": "Click chromosomes in order: 1a, 1b, 1c, 2a…",
    "signal": "Click a signal; it is assigned to the nearest chromosome",
    "measure": "Drag from telomere to telomere through the centromere",
}
HIT_RADIUS = 18.0  # image pixels


def pick_karyotype(karyotypes: list[Karyotype], name: str = "", species: str = "") -> Karyotype | None:
    """The image's own karyotype, else one for its species, else the first."""
    if name:
        for k in karyotypes:
            if k.name == name:
                return k
    if species:
        for k in karyotypes:
            if k.species.lower() == species.lower() or k.name.lower() == species.lower():
                return k
    return karyotypes[0] if karyotypes else None


def split_label(karyotype: Karyotype, label: str) -> tuple[str, str] | None:
    """"2b" -> ("2", "b") using the karyotype's chromosome names; None if the label does not exist."""
    for chrom in karyotype.chromosomes:
        for hom in karyotype.homologue_labels(chrom):
            if hom == label:
                return chrom.name, label[len(chrom.name):]
    return None


@dataclass
class CellStatus:
    image_id: str
    label: str
    marked: bool


class SignalsSession:
    def __init__(self, folder: Path, image: CytoImage, karyotype: Karyotype | None = None) -> None:
        self.folder = Path(folder).expanduser().resolve()
        self.image = image
        self.view = ImageSession(self.folder, image)
        self.details: ImageDetails = self.view.details
        self.karyotypes = list_karyotypes(self.folder)
        slide_species = ""
        for s in load_slides(self.folder):
            if s.id == self.details.slide_id:
                slide_species = s.species
        species = image.metadata.get("species") or slide_species
        self.karyotype: Karyotype | None = karyotype or pick_karyotype(self.karyotypes, self.details.karyotype, species)
        self.dirty = False

    # ---------------------------------------------------------------- karyotype and labels

    def set_karyotype(self, karyotype: Karyotype | None) -> None:
        self.karyotype = karyotype
        self.details.karyotype = karyotype.name if karyotype else ""
        self.dirty = True

    def labels(self) -> list[str]:
        return self.karyotype.all_labels() if self.karyotype else []

    def used_labels(self) -> set[str]:
        return {m.label for m in self.details.chromosomes}

    def next_label(self) -> str:
        """First karyotype label not yet numbered in this cell ("" when all are used)."""
        used = self.used_labels()
        return next((l for l in self.labels() if l not in used), "")

    def mark_for_label(self, label: str) -> ChromosomeMark | None:
        return next((m for m in self.details.chromosomes if m.label == label), None)

    # ---------------------------------------------------------------- marks

    def add_mark(self, label: str, x: float, y: float) -> ChromosomeMark:
        parts = split_label(self.karyotype, label) if self.karyotype else None
        if parts is None:
            raise ValueError(f"{label} is not in the karyotype")
        if self.mark_for_label(label) is not None:
            raise ValueError(f"{label} is already numbered in this cell")
        mark = ChromosomeMark(id=new_id(), chromosome=parts[0], homologue=parts[1], x=float(x), y=float(y))
        self.details.chromosomes.append(mark)
        self.reclassify()
        self.dirty = True
        return mark

    def move_mark(self, mark_id: str, x: float, y: float) -> None:
        mark = self.details.chromosome(mark_id)
        if mark is not None:
            mark.x, mark.y = float(x), float(y)
            self.dirty = True

    def remove_mark(self, mark_id: str) -> None:
        self.details.chromosomes = [m for m in self.details.chromosomes if m.id != mark_id]
        for s in self.details.signals:
            if s.chromosome_id == mark_id:
                s.chromosome_id = ""
                if not s.manual:
                    s.arm, s.position, s.rel_distance = "", "", None
        self.reclassify()
        self.dirty = True

    def relabel_mark(self, mark_id: str, label: str) -> None:
        mark = self.details.chromosome(mark_id)
        parts = split_label(self.karyotype, label) if self.karyotype else None
        if mark is None or parts is None:
            raise ValueError(f"{label} is not in the karyotype")
        other = self.mark_for_label(label)
        if other is not None and other.id != mark_id:
            raise ValueError(f"{label} is already numbered in this cell")
        mark.chromosome, mark.homologue = parts
        self.dirty = True

    def set_axis(self, mark_id: str, p_tel: tuple[float, float], q_tel: tuple[float, float]) -> None:
        """Telomere to telomere; the centromere is put where the karyotype's centromere index says."""
        mark = self.details.chromosome(mark_id)
        if mark is None:
            return
        chrom = self.karyotype.chromosome(mark.chromosome) if self.karyotype else None
        f = chrom.centromere_fraction if chrom is not None and chrom.centromere_index > 0 else 0.5
        mark.telomere_p, mark.telomere_q = (float(p_tel[0]), float(p_tel[1])), (float(q_tel[0]), float(q_tel[1]))
        mark.centromere = (p_tel[0] + (q_tel[0] - p_tel[0]) * f, p_tel[1] + (q_tel[1] - p_tel[1]) * f)
        self.reclassify()
        self.dirty = True

    def set_centromere(self, mark_id: str, point: tuple[float, float]) -> None:
        mark = self.details.chromosome(mark_id)
        if mark is not None:
            mark.centromere = (float(point[0]), float(point[1]))
            self.reclassify()
            self.dirty = True

    # ---------------------------------------------------------------- signals

    def probes(self) -> list[tuple[str, str]]:
        """(probe, dye) of the channels that have a probe, for the Probe menu."""
        out = []
        for c in self.details.channels:
            if c.probe and (c.probe, c.dye) not in out:
                out.append((c.probe, c.dye))
        return out

    def probe_colour(self, probe: str) -> str:
        lists = load_lists(self.folder)
        for c in self.details.channels:
            if c.probe == probe:
                return c.display.color
        return dye_color(lists, "") if probe else "#ffffff"

    def refine_position(self, probe: str, x: float, y: float, radius: float = 14.0) -> tuple[float, float]:
        """Move a click to the centre of the bright spot of the probe's channel nearby."""
        channel_index = next((i for i, c in enumerate(self.details.channels) if c.probe == probe), None)
        if channel_index is None:
            return float(x), float(y)
        array = ch.to_grey(channel_layers(self.folder, self.details, self.view.cache)[channel_index][0]).astype(np.float64)
        h, w = array.shape
        r = int(radius)
        x0, x1 = max(0, int(x) - r), min(w, int(x) + r + 1)
        y0, y1 = max(0, int(y) - r), min(h, int(y) + r + 1)
        if x1 <= x0 or y1 <= y0:
            return float(x), float(y)
        window = array[y0:y1, x0:x1]
        lo, hi = float(window.min()), float(window.max())
        if hi <= lo:
            return float(x), float(y)
        weights = np.where(window >= lo + 0.5 * (hi - lo), window - lo, 0.0)
        total = weights.sum()
        if total <= 0:
            return float(x), float(y)
        ys, xs = np.mgrid[y0:y1, x0:x1]
        return float((weights * xs).sum() / total), float((weights * ys).sum() / total)

    def add_signal(self, probe: str, x: float, y: float, *, refine: bool = True, strength: str = "") -> Signal:
        if refine:
            x, y = self.refine_position(probe, x, y)
        signal = Signal(id=new_id(), probe=probe, x=float(x), y=float(y), strength=strength)
        sg.assign_and_classify(self.details, signal, max_distance=None)
        self.details.signals.append(signal)
        self.dirty = True
        return signal

    def signal(self, signal_id: str) -> Signal | None:
        return next((s for s in self.details.signals if s.id == signal_id), None)

    def update_signal(self, signal_id: str, **fields) -> None:
        """Set probe, arm, position, strength or chromosome_id. Arm and position by hand are kept ("manual")."""
        s = self.signal(signal_id)
        if s is None:
            return
        for key, value in fields.items():
            if key in ("probe", "strength", "arm", "position", "chromosome_id"):
                setattr(s, key, value)
        if "arm" in fields or "position" in fields:
            s.manual = True
        if "chromosome_id" in fields and not s.manual:
            sg.classify(s, self.details.chromosome(s.chromosome_id))
        self.dirty = True

    def move_signal(self, signal_id: str, x: float, y: float) -> None:
        s = self.signal(signal_id)
        if s is None:
            return
        s.x, s.y = float(x), float(y)
        if not s.manual:
            s.chromosome_id = ""
        self.reclassify()
        self.dirty = True

    def remove_signal(self, signal_id: str) -> None:
        self.details.signals = [s for s in self.details.signals if s.id != signal_id]
        self.dirty = True

    def reclassify(self) -> None:
        for s in self.details.signals:
            if s.chromosome_id and self.details.chromosome(s.chromosome_id) is None:
                s.chromosome_id = ""
            if not s.manual:
                if not s.chromosome_id:
                    mark = sg.nearest_chromosome(s.x, s.y, self.details.chromosomes)
                    s.chromosome_id = mark.id if mark else ""
                sg.classify(s, self.details.chromosome(s.chromosome_id))

    def distance_text(self, signal: Signal) -> str:
        if signal.rel_distance is None:
            return "not measured (set the centromere and telomeres with Measure)"
        arm = signal.arm
        if arm == "cen":
            return f"{signal.rel_distance:.2f} of the arm length, at the centromere"
        return f"{signal.rel_distance:.2f} of arm length from centromere"

    # ---------------------------------------------------------------- hit testing

    def hit(self, x: float, y: float, tolerance: float = HIT_RADIUS) -> tuple[str, str] | None:
        """("signal", id) or ("mark", id) at a point; signals win."""
        best = None
        for s in self.details.signals:
            d = math.hypot(s.x - x, s.y - y)
            if d <= tolerance and (best is None or d < best[0]):
                best = (d, ("signal", s.id))
        if best:
            return best[1]
        for m in self.details.chromosomes:
            d = sg.distance_to_mark(x, y, m)
            if d <= tolerance * 1.5 and (best is None or d < best[0]):
                best = (d, ("mark", m.id))
        return best[1] if best else None

    # ---------------------------------------------------------------- the slide

    def slide_cells(self) -> list[CellStatus]:
        slide_id = self.details.slide_id
        images = [i for i in load_images(self.folder) if image_details(i).slide_id == slide_id]
        images.sort(key=lambda i: natural_key(i.original_name or i.id))
        out = []
        for n, i in enumerate(images, start=1):
            d = self.details if i.id == self.image.id else image_details(i)
            out.append(CellStatus(i.id, _short_cell_name(i, n), bool(d.signals or d.chromosomes)))
        return out

    def summary(self) -> tuple[list[sg.SummaryRow], int]:
        """Rows over all cells of the slide (this cell with unsaved changes) and the number of cells marked."""
        slide_id = self.details.slide_id
        details = []
        for i in load_images(self.folder):
            d = image_details(i)
            if d.slide_id != slide_id:
                continue
            details.append(self.details if i.id == self.image.id else d)
        marked = [d for d in details if d.chromosomes]
        order = [c.name for c in self.karyotype.chromosomes] if self.karyotype else None
        return sg.summarize(marked, order), len(marked)

    def unsignalled(self) -> list[str]:
        slide_id = self.details.slide_id
        details = [self.details if i.id == self.image.id else image_details(i)
                   for i in load_images(self.folder) if image_details(i).slide_id == slide_id]
        marked = [d for d in details if d.chromosomes]
        rows = sg.summarize(marked)
        return sg.unsignalled_chromosomes(marked, rows)

    # ---------------------------------------------------------------- saving

    def save(self) -> None:
        """Replace only karyotype, chromosomes and signals of the stored record."""
        fresh = next((i for i in load_images(self.folder) if i.id == self.image.id), self.image)
        stored = image_details(fresh)
        stored.karyotype = self.details.karyotype
        stored.chromosomes = self.details.chromosomes
        stored.signals = self.details.signals
        set_image_details(fresh, stored)
        update_image(self.folder, fresh)
        self.dirty = False


def _short_cell_name(image: CytoImage, number: int) -> str:
    """"s12_c03" -> "03"; falls back to the running number."""
    import re

    m = re.search(r"(\d+)\D*$", image.original_name or "")
    return m.group(1).zfill(2) if m else str(number).zfill(2)
