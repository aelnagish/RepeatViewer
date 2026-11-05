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
TOOLS = (("number", "2  Number chromosomes"), ("signal", "3  Mark signals"))
TOOL_HINTS = {
    "number": "Click a chromosome to name it {next}. Drag a name, a telomere or the centromere to move it. Right-click for more.",
    "signal": "Click a spot to add a {probe} signal. Drag a signal to move it. Right-click for more.",
}
COUNTERSTAINS = ("dapi", "hoechst", "pi")
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

    # ---------------------------------------------------------------- images used for snapping

    picture: np.ndarray | None = None  # RGB of a plain picture (set by the view, which can read any format)

    def set_picture(self, rgb: np.ndarray | None) -> None:
        self.picture = rgb

    def _counterstain(self) -> np.ndarray | None:
        """Grey image of the chromosomes: the DAPI-like channel, else the blue of a picture."""
        if self.details.channels:
            layers = channel_layers(self.folder, self.details, self.view.cache)
            idx = next((i for i, c in enumerate(self.details.channels) if c.dye.lower() in COUNTERSTAINS), 0)
            return ch.to_grey(layers[idx][0]).astype(np.float64)
        if self.picture is not None and self.picture.ndim == 3:
            rgb = self.picture[..., :3].astype(np.float64)
            blue = rgb[..., 2]
            return blue if blue.mean() > 0.6 * rgb.mean() else rgb.mean(axis=2)
        return None

    def estimate_axis(self, x: float, y: float, radius: int = 160):
        """(telomere, centromere or None, telomere) of the chromosome under (x, y), from the counterstain.

        The chromosome is the bright patch around the click; its long axis is its main direction and the
        centromere the narrowest place between 15 and 85 % of its length. None when nothing bright is there.
        """
        grey = self._counterstain()
        if grey is None:
            return None
        h, w = grey.shape
        x0, x1 = max(0, int(x) - radius), min(w, int(x) + radius + 1)
        y0, y1 = max(0, int(y) - radius), min(h, int(y) + radius + 1)
        win = grey[y0:y1, x0:x1]
        if win.size < 50:
            return None
        thr = _otsu(win)
        mask = win > thr
        cx, cy = int(x) - x0, int(y) - y0
        if not (0 <= cy < mask.shape[0] and 0 <= cx < mask.shape[1]):
            return None
        if not mask[cy, cx]:
            ys, xs = np.nonzero(mask[max(0, cy - 8):cy + 9, max(0, cx - 8):cx + 9])
            if not len(xs):
                return None
            k = int(np.argmin((xs - min(cx, 8)) ** 2 + (ys - min(cy, 8)) ** 2))
            cx, cy = max(0, cx - 8) + int(xs[k]), max(0, cy - 8) + int(ys[k])
        comp = _component(mask, cy, cx)
        if comp.shape[0] < 30:
            return None
        pts = comp[:, ::-1].astype(np.float64)  # (x, y)
        mean = pts.mean(axis=0)
        cov = np.cov((pts - mean).T)
        vals, vecs = np.linalg.eigh(cov)
        v = vecs[:, int(np.argmax(vals))]
        t = (pts - mean) @ v
        tmin, tmax = float(t.min()), float(t.max())
        a = mean + v * tmin + np.array([x0, y0])
        b = mean + v * tmax + np.array([x0, y0])
        cen = None
        if tmax - tmin > 12:
            bins = 40
            hist, _edges = np.histogram(t, bins=bins, range=(tmin, tmax))
            smooth = np.convolve(hist, np.ones(3) / 3, mode="same")
            lo, hi = int(bins * 0.15), int(bins * 0.85)
            k = lo + int(np.argmin(smooth[lo:hi]))
            if smooth[k] < 0.8 * float(np.median(smooth[lo:hi])):
                tc = tmin + (k + 0.5) * (tmax - tmin) / bins
                cen = mean + v * tc + np.array([x0, y0])
        return (float(a[0]), float(a[1])), (None if cen is None else (float(cen[0]), float(cen[1]))), (float(b[0]), float(b[1]))

    def auto_axis(self, mark_id: str) -> bool:
        """Set telomeres and centromere of a numbered chromosome from the image. p is the shorter arm."""
        mark = self.details.chromosome(mark_id)
        if mark is None:
            return False
        found = self.estimate_axis(mark.x, mark.y)
        if found is None:
            return False
        a, cen, b = found
        if cen is None:
            self.set_axis(mark_id, a, b)
            return True
        da = ((a[0] - cen[0]) ** 2 + (a[1] - cen[1]) ** 2) ** 0.5
        db = ((b[0] - cen[0]) ** 2 + (b[1] - cen[1]) ** 2) ** 0.5
        p_tel, q_tel = (a, b) if da <= db else (b, a)
        mark.telomere_p, mark.telomere_q, mark.centromere = p_tel, q_tel, cen
        self.reclassify()
        self.dirty = True
        return True

    def set_point(self, mark_id: str, which: str, point: tuple[float, float]) -> None:
        """Move one handle of a chromosome: "telomere_p", "telomere_q" or "centromere"."""
        mark = self.details.chromosome(mark_id)
        if mark is None or which not in ("telomere_p", "telomere_q", "centromere"):
            return
        setattr(mark, which, (float(point[0]), float(point[1])))
        self.reclassify()
        self.dirty = True

    def swap_arms(self, mark_id: str) -> None:
        mark = self.details.chromosome(mark_id)
        if mark is not None and mark.telomere_p is not None and mark.telomere_q is not None:
            mark.telomere_p, mark.telomere_q = mark.telomere_q, mark.telomere_p
            self.reclassify()
            self.dirty = True

    # ---------------------------------------------------------------- signals

    def probes(self) -> list[tuple[str, str]]:
        """(probe, dye) for the Probe menu: the channels' probes first, then the other probes of the lists."""
        out = []
        for c in self.details.channels:
            if c.probe and (c.probe, c.dye) not in out:
                out.append((c.probe, c.dye))
        known = {p for p, _d in out}
        for e in load_lists(self.folder).get("probes", []):
            name = str(e.get("name") or "")
            if name and name not in known:
                out.append((name, str(e.get("dye") or "")))
                known.add(name)
        return out

    def probe_colour(self, probe: str) -> str:
        for c in self.details.channels:
            if c.probe == probe:
                return c.display.color
        lists = load_lists(self.folder)
        e = next((e for e in lists.get("probes", []) if e.get("name") == probe), None)
        dye = str((e or {}).get("dye") or "")
        return dye_color(lists, dye, "#ffffff") if dye else "#ffffff"

    def refine_position(self, probe: str, x: float, y: float, radius: float = 14.0) -> tuple[float, float]:
        """Move a click to the centre of the bright spot of the probe's channel nearby."""
        channel_index = next((i for i, c in enumerate(self.details.channels) if c.probe == probe), None)
        if channel_index is not None:
            array = ch.to_grey(channel_layers(self.folder, self.details, self.view.cache)[channel_index][0]).astype(np.float64)
        elif self.picture is not None and self.picture.ndim == 3:
            colour = np.array(ch.hex_to_rgb(self.probe_colour(probe)))
            colour = colour / (np.linalg.norm(colour) or 1.0)
            rgb = self.picture[..., :3].astype(np.float64)
            array = rgb @ colour - 0.5 * (rgb.sum(axis=2) - rgb @ colour)
        else:
            return float(x), float(y)
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


def _otsu(values: np.ndarray) -> float:
    data = values.astype(np.float64).ravel()
    lo, hi = float(data.min()), float(data.max())
    if hi <= lo:
        return hi
    hist, edges = np.histogram(data, bins=256, range=(lo, hi))
    centers = (edges[:-1] + edges[1:]) / 2
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * centers)
    mean0 = m0 / np.maximum(w0, 1)
    mean1 = (m0[-1] - m0) / np.maximum(w1, 1)
    between = w0 * w1 * (mean0 - mean1) ** 2
    return float(centers[int(np.argmax(between))])


def _component(mask: np.ndarray, row: int, col: int) -> np.ndarray:
    """(row, col) of the 4-connected patch of True around one pixel."""
    try:
        from scipy import ndimage

        labels, _n = ndimage.label(mask)
        lab = labels[row, col]
        return np.argwhere(labels == lab) if lab else np.zeros((0, 2), dtype=int)
    except ImportError:
        pass
    seen = np.zeros(mask.shape, dtype=bool)
    stack = [(row, col)]
    seen[row, col] = True
    out = []
    h, w = mask.shape
    while stack:
        r, c = stack.pop()
        out.append((r, c))
        for rr, cc in ((r + 1, c), (r - 1, c), (r, c + 1), (r, c - 1)):
            if 0 <= rr < h and 0 <= cc < w and mask[rr, cc] and not seen[rr, cc]:
                seen[rr, cc] = True
                stack.append((rr, cc))
    return np.array(out, dtype=int)
