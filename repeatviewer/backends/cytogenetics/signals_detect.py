"""Detection and review for the Signals tab (mixed into ``SignalsSession``). No Qt.

Step 1 finds the chromosomes in the counterstain and lets the user split, merge, paint and
erase them; step 2 finds the signals of every probe channel. Everything found is stored in
the ordinary chromosome marks and signals of the image, so the Idiogram, Figure and the
summaries read it as before. Detected chromosomes stay unnamed ("?") until step 3 pairs them.
"""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

import numpy as np

from . import channels as ch
from . import segment as sm
from . import signals as sg
from . import spots as sp
from .importer import channel_layers
from .lists import add_entry, dye_color, load_lists, save_lists
from .model import UNNAMED, ChromosomeMark, Signal, new_id

COUNTERSTAINS = ("dapi", "hoechst", "pi")
PLANES = ("Red", "Green", "Blue")
CHECK_WEAK = "weak signal"
CHECK_OFF = "not on a chromosome"


@dataclass
class ProbeRow:
    """One row of the Probes table: a probe channel of the cell (or a colour of a plain picture)."""

    key: str  # channel id, or "plane:<0-2>" for pictures
    number: int  # channel number shown to the user (1-based)
    file: str  # short channel label ("ch2_FITC.tif", "Green")
    probe: str  # "" = no probe yet
    dye: str
    colour: str
    count: int = 0
    visible: bool = True

    @property
    def channel_text(self) -> str:
        return f"{self.number} · {self.dye}" if self.dye and self.probe else f"{self.number} · {self.file}"


def short_channel_name(name: str) -> str:
    """"cell05_FITC.tif" stays; "1571 - 158_2.tif · Layer 2" (a layer of a Photoshop file) becomes "Layer 2"."""
    name = str(name).replace("\\", "/").rsplit("/", 1)[-1]
    if " · " in name:
        return name.rsplit(" · ", 1)[-1]
    return name


class DetectMixin:
    """Needs ``folder``, ``details``, ``view`` (ImageSession), ``picture``, ``dirty`` and ``reclassify``."""

    _labels_cache: tuple | None = None

    # ---------------------------------------------------------------- settings

    def detect_settings(self) -> sm.DetectSettings:
        return sm.DetectSettings.from_dict(self.details.detection.get("chromosomes"))

    def set_detect_settings(self, settings: sm.DetectSettings) -> None:
        self.details.detection["chromosomes"] = settings.to_dict()
        self.dirty = True

    def spot_settings(self, probe: str) -> sp.SpotSettings:
        return sp.SpotSettings.from_dict((self.details.detection.get("probes") or {}).get(probe))

    def set_spot_settings(self, probe: str, settings: sp.SpotSettings) -> None:
        self.details.detection.setdefault("probes", {})[probe] = settings.to_dict()
        self.dirty = True

    def hidden_probes(self) -> set[str]:
        return set(self.details.detection.get("hidden") or [])

    def set_probe_visible(self, probe: str, visible: bool) -> None:
        hidden = self.hidden_probes()
        (hidden.discard if visible else hidden.add)(probe)
        self.details.detection["hidden"] = sorted(hidden)
        self.dirty = True

    # ---------------------------------------------------------------- arrays

    def _layers(self):
        return channel_layers(self.folder, self.details, self.view.cache)

    def counterstain_index(self) -> int | None:
        """The channel the chromosomes are found in: the one chosen in "Chromosomes from", else the DAPI-like
        dye, else the first channel that is not a Photoshop "Background" layer."""
        chans = self.details.channels
        if not chans:
            return None
        chosen = self.details.detection.get("counterstain")
        if chosen:
            for i, c in enumerate(chans):
                if c.id == chosen:
                    return i
        for i, c in enumerate(chans):
            if re.sub(r"[^a-z]", "", c.dye.lower()) in COUNTERSTAINS:
                return i
        for i, c in enumerate(chans):
            if short_channel_name(c.original_name or c.file).lower() != "background":
                return i
        return 0

    def counterstain_chosen(self) -> bool:
        """True when the user chose the channel or a channel has a counterstain dye (no guessing needed)."""
        if self.details.detection.get("counterstain"):
            return True
        return any(re.sub(r"[^a-z]", "", c.dye.lower()) in COUNTERSTAINS for c in self.details.channels)

    def choose_counterstain(self, channel_id: str) -> None:
        """"" = automatic."""
        if channel_id:
            self.details.detection["counterstain"] = channel_id
        else:
            self.details.detection.pop("counterstain", None)
        self._labels_cache = None
        self.dirty = True

    def guess_counterstain(self) -> int | None:
        """The channel that looks like chromosomes: not flat, and the largest bright area (probes are small
        spots, the counterstain covers whole chromosomes). Reads every channel: run it off the GUI thread."""
        layers = self._layers()
        best, best_score = None, -1.0
        for i, (array, _display) in enumerate(layers):
            g = ch.to_grey(array).astype(np.float64)
            g = g[:: max(1, g.shape[0] // 600), :: max(1, g.shape[1] // 600)]
            if g.std() < 1e-6:
                continue
            share = float((g > sm.otsu(g)).mean())
            if share > best_score:
                best, best_score = i, share
        return best

    def counterstain(self) -> np.ndarray | None:
        """Grey counterstain of the cell: the DAPI channel, else the blue of a picture."""
        if self.details.channels:
            idx = self.counterstain_index()
            return ch.to_grey(self._layers()[idx][0]).astype(np.float64)
        if self.picture is not None and self.picture.ndim == 3:
            rgb = self.picture[..., :3].astype(np.float64)
            blue = rgb[..., 2]
            return blue if blue.mean() > 0.6 * rgb.mean() else rgb.mean(axis=2)
        return None

    def _measure_grey(self) -> np.ndarray | None:
        g = self.counterstain()
        return None if g is None else sm.prepared(g)

    def image_shape(self) -> tuple[int, int] | None:
        g = self.counterstain()
        return None if g is None else g.shape[:2]

    def probe_array(self, row: ProbeRow) -> np.ndarray | None:
        if row.key.startswith("plane:"):
            if self.picture is None or self.picture.ndim != 3:
                return None
            plane = int(row.key.split(":", 1)[1])
            rgb = self.picture[..., :3].astype(np.float64)
            # the plane, minus what the other two planes share (white / blue counterstain stays out)
            others = np.delete(rgb, plane, axis=2).max(axis=2)
            return np.clip(rgb[..., plane] - 0.8 * others, 0, None)
        idx = next((i for i, c in enumerate(self.details.channels) if c.id == row.key), None)
        if idx is None:
            return None
        return ch.to_grey(self._layers()[idx][0]).astype(np.float64)

    # ---------------------------------------------------------------- probe table

    def probe_rows(self) -> list[ProbeRow]:
        hidden = self.hidden_probes()
        counts: dict[str, int] = {}
        for s in self.details.signals:
            counts[s.probe] = counts.get(s.probe, 0) + 1
        rows = []
        if self.details.channels:
            cs = self.counterstain_index()
            for i, c in enumerate(self.details.channels):
                if i == cs or (short_channel_name(c.original_name or c.file).lower() == "background" and not c.probe):
                    continue  # the counterstain, and the empty "Background" layer of a Photoshop file
                rows.append(ProbeRow(c.id, i + 1, short_channel_name(c.original_name or c.file), c.probe, c.dye, c.display.color or "#ffffff",
                                     counts.get(c.probe, 0) if c.probe else 0, c.probe not in hidden))
            return rows
        for e in self.details.detection.get("picture_probes") or []:
            plane = int(e.get("plane", 0)) % 3
            probe = str(e.get("probe") or "")
            rows.append(ProbeRow(f"plane:{plane}", plane + 1, PLANES[plane], probe, str(e.get("dye") or ""),
                                 str(e.get("colour") or "#ffffff"), counts.get(probe, 0), probe not in hidden))
        return rows

    def free_rows(self) -> list[ProbeRow]:
        return [r for r in self.probe_rows() if not r.probe]

    def free_planes(self) -> list[int]:
        used = {int(e.get("plane", 0)) for e in self.details.detection.get("picture_probes") or []}
        return [p for p in range(3) if p not in used and p != 2]

    def assign_probe(self, key: str, probe: str, dye: str = "", colour: str = "") -> ProbeRow | None:
        """Put a probe on a channel (or a picture colour). Returns its row."""
        lists = load_lists(self.folder)
        if not dye:
            e = next((e for e in lists.get("probes", []) if str(e.get("name", "")).lower() == probe.lower()), None)
            dye = str((e or {}).get("dye") or "")
        colour = colour or (dye_color(lists, dye, "") if dye else "") or "#ffffff"
        if key.startswith("plane:"):
            plane = int(key.split(":", 1)[1])
            entries = [e for e in self.details.detection.get("picture_probes") or [] if int(e.get("plane", -1)) != plane]
            entries.append({"probe": probe, "dye": dye, "colour": colour, "plane": plane})
            self.details.detection["picture_probes"] = entries
        else:
            channel = next((c for c in self.details.channels if c.id == key), None)
            if channel is None:
                return None
            old = channel.probe
            channel.probe, channel.dye = probe, dye or channel.dye
            channel.display.color = colour or channel.display.color
            if old and old != probe:
                for s in self.details.signals:
                    if s.probe == old:
                        s.probe = probe
        self.dirty = True
        return next((r for r in self.probe_rows() if r.key == key), None)

    def remove_probe(self, key: str) -> None:
        """Take the probe off a channel; its signals go too."""
        row = next((r for r in self.probe_rows() if r.key == key), None)
        if row is None:
            return
        if key.startswith("plane:"):
            plane = int(key.split(":", 1)[1])
            self.details.detection["picture_probes"] = [e for e in self.details.detection.get("picture_probes") or []
                                                        if int(e.get("plane", -1)) != plane]
        else:
            channel = next((c for c in self.details.channels if c.id == key), None)
            if channel is not None:
                channel.probe = ""
        if row.probe:
            self.details.signals = [s for s in self.details.signals if s.probe != row.probe]
        self.dirty = True

    def add_probe_to_lists(self, name: str, dye: str, sequence: str = "") -> None:
        lists = load_lists(self.folder)
        entry = add_entry(lists, "probes", name, dye=dye)
        if sequence:
            entry["sequence"] = sequence
        if dye:
            add_entry(lists, "dyes", dye)
        save_lists(self.folder, lists)

    # ---------------------------------------------------------------- chromosomes

    def labels_map(self) -> np.ndarray | None:
        """int32 map, 0 = background, i + 1 = ``details.chromosomes[i]`` (marks with a mask only)."""
        shape = self.image_shape()
        if shape is None:
            return None
        key = (shape, tuple((m.id, id(m.mask)) for m in self.details.chromosomes))
        if self._labels_cache is not None and self._labels_cache[0] == key:
            return self._labels_cache[1]
        out = np.zeros(shape, dtype=np.int32)
        for i, m in enumerate(self.details.chromosomes):
            dec = sm.decode_mask(m.mask) if m.mask else None
            if dec is not None:
                out[sm.paste(shape, *dec)] = i + 1
        self._labels_cache = (key, out)
        return out

    def mark_mask(self, mark: ChromosomeMark):
        return sm.decode_mask(mark.mask) if mark.mask else None

    def _apply_measure(self, mark: ChromosomeMark, measured: sm.Measured) -> None:
        mark.mask = sm.encode_mask(measured.mask, measured.x0, measured.y0)
        line = measured.midline
        if len(line) >= 2:
            mark.midline = [[round(x, 2), round(y, 2)] for x, y in line]
            mark.telomere_p = tuple(line[0])
            mark.telomere_q = tuple(line[-1])
            k = measured.centromere if measured.centromere >= 0 else len(line) // 2
            mark.centromere = tuple(line[k])
            mid = line[len(line) // 2]
            mark.x, mark.y = float(mid[0]), float(mid[1])
        else:
            cx, cy = measured.centre()
            mark.x, mark.y = cx, cy
            mark.midline = None
        mark.check = "; ".join(measured.checks)
        mark.auto = True

    def _renumber_unnamed(self) -> None:
        unnamed = [m for m in self.details.chromosomes if not m.named]
        unnamed.sort(key=lambda m: -sm.polyline_length(m.midline or []))
        for i, m in enumerate(unnamed, start=1):
            m.homologue = str(i)

    def detect_chromosomes(self, settings: sm.DetectSettings | None = None) -> int:
        """Replace the chromosomes of the cell by the ones found in the counterstain. Returns how many."""
        grey = self.counterstain()
        if grey is None:
            raise ValueError("This cell has no counterstain (DAPI) image")
        settings = settings or self.detect_settings()
        return self.apply_chromosomes(sm.detect_chromosomes(grey, settings), settings)

    def apply_chromosomes(self, found: list, settings: sm.DetectSettings) -> int:
        """Store chromosomes found by ``segment.detect_chromosomes`` (run in a worker thread by the view)."""
        self.set_detect_settings(settings)
        marks = []
        for i, obj in enumerate(found, start=1):
            mark = ChromosomeMark(id=new_id(), chromosome=UNNAMED, homologue=str(i))
            self._apply_measure(mark, obj)
            marks.append(mark)
        self.details.chromosomes = marks
        for s in self.details.signals:
            s.chromosome_id = ""
        self._reassign_signals()
        self.dirty = True
        return len(marks)

    def remeasure(self, mark_id: str) -> None:
        mark = self.details.chromosome(mark_id)
        dec = self.mark_mask(mark) if mark else None
        if dec is None:
            return
        measured = sm.measure(dec[0], dec[1], dec[2], self.image_shape(), self._measure_grey())
        others = [sm.Measured(*d) for d in (self.mark_mask(m) for m in self.details.chromosomes if m.id != mark_id) if d]
        others.append(measured)
        sm.flag_sizes(others) if len(others) >= 3 else None
        self._apply_measure(mark, measured)
        self._reassign_signals()
        self.dirty = True

    def _new_mark_from(self, crop, x0, y0, like: ChromosomeMark | None = None) -> ChromosomeMark:
        mark = ChromosomeMark(id=new_id(), chromosome=UNNAMED, homologue="0")
        self._apply_measure(mark, sm.measure(crop, x0, y0, self.image_shape(), self._measure_grey()))
        return mark

    def split_mark(self, mark_id: str, a: tuple[float, float], b: tuple[float, float]) -> list[str]:
        """Cut a chromosome along a line drawn across it. Returns the ids of the pieces."""
        mark = self.details.chromosome(mark_id)
        dec = self.mark_mask(mark) if mark else None
        if dec is None:
            return []
        pieces = sm.split_along(dec[0], dec[1], dec[2], a, b)
        if len(pieces) < 2:
            return [mark_id]
        idx = self.details.chromosomes.index(mark)
        first = pieces[0]
        self._apply_measure(mark, sm.measure(first[0], first[1], first[2], self.image_shape(), self._measure_grey()))
        ids = [mark.id]
        for crop, x0, y0 in pieces[1:]:
            if crop.sum() < 6:
                continue
            new = self._new_mark_from(crop, x0, y0)
            self.details.chromosomes.insert(idx + 1, new)
            ids.append(new.id)
        if mark.named:
            for i in ids[1:]:
                m = self.details.chromosome(i)
                m.chromosome, m.homologue = UNNAMED, "0"
        self._renumber_unnamed()
        self._reassign_signals()
        self.dirty = True
        return ids

    def split_at_line(self, a: tuple[float, float], b: tuple[float, float]) -> list[str]:
        """Split every chromosome the line crosses."""
        labels = self.labels_map()
        if labels is None:
            return []
        n = max(2, int(math.hypot(b[0] - a[0], b[1] - a[1])))
        hit = []
        h, w = labels.shape
        for t in np.linspace(0, 1, n):
            x, y = a[0] + (b[0] - a[0]) * t, a[1] + (b[1] - a[1]) * t
            xi, yi = int(round(x)), int(round(y))
            if 0 <= yi < h and 0 <= xi < w and labels[yi, xi]:
                mid = self.details.chromosomes[labels[yi, xi] - 1].id
                if mid not in hit:
                    hit.append(mid)
        out = []
        for mid in hit:
            out.extend(self.split_mark(mid, a, b))
        return out

    def merge_marks(self, ids: list[str]) -> str:
        """One chromosome from several; keeps the first one's name."""
        marks = [m for m in (self.details.chromosome(i) for i in ids) if m is not None]
        parts = [d for d in (self.mark_mask(m) for m in marks) if d is not None]
        if len(parts) < 2:
            return marks[0].id if marks else ""
        crop, x0, y0 = sm.union(parts)
        keep = marks[0]
        self._apply_measure(keep, sm.measure(crop, x0, y0, self.image_shape(), self._measure_grey()))
        drop = {m.id for m in marks[1:]}
        self.details.chromosomes = [m for m in self.details.chromosomes if m.id not in drop]
        for s in self.details.signals:
            if s.chromosome_id in drop:
                s.chromosome_id = ""
        self._renumber_unnamed()
        self._reassign_signals()
        self.dirty = True
        return keep.id

    def brush_mark(self, mark_id: str, x: float, y: float, radius: float, add: bool, *, measure_now: bool = True) -> str:
        """Paint onto (or erase from) a chromosome; painting on empty space with no id makes a new one.

        Returns the id that was changed ("" when an erase removed the whole chromosome).
        """
        mark = self.details.chromosome(mark_id) if mark_id else None
        if mark is None:
            if not add:
                return ""
            crop = np.ones((1, 1), dtype=bool)
            mark = ChromosomeMark(id=new_id(), chromosome=UNNAMED, homologue="0", x=x, y=y, auto=True)
            mark.mask = sm.encode_mask(crop, int(x), int(y))
            self.details.chromosomes.append(mark)
        dec = self.mark_mask(mark)
        if dec is None:
            return ""
        crop, x0, y0 = sm.brush(dec[0], dec[1], dec[2], x, y, radius, add)
        if not crop.any():
            self.remove_mark(mark.id)
            return ""
        mark.mask = sm.encode_mask(crop, x0, y0)
        if measure_now:
            self.remeasure(mark.id)
            self._renumber_unnamed()
        self.dirty = True
        return mark.id

    def add_chromosome_at(self, x: float, y: float) -> str:
        """Add the bright object under a click as a chromosome (for one the detection missed)."""
        grey = self.counterstain()
        if grey is None:
            return ""
        found = sm.object_at(grey, x, y, self.detect_settings().sensitivity)
        if found is None:
            return ""
        crop, x0, y0 = found
        labels = self.labels_map()
        if labels is not None:
            taken = labels[y0:y0 + crop.shape[0], x0:x0 + crop.shape[1]] > 0
            crop = crop & ~taken
            if crop.sum() < 6:
                return ""
        mark = self._new_mark_from(crop, x0, y0)
        self.details.chromosomes.append(mark)
        self._renumber_unnamed()
        self._reassign_signals()
        self.dirty = True
        return mark.id

    def mark_at(self, x: float, y: float) -> ChromosomeMark | None:
        labels = self.labels_map()
        if labels is not None:
            xi, yi = int(round(x)), int(round(y))
            if 0 <= yi < labels.shape[0] and 0 <= xi < labels.shape[1] and labels[yi, xi]:
                return self.details.chromosomes[labels[yi, xi] - 1]
        return None

    def move_centromere(self, mark_id: str, x: float, y: float) -> None:
        """Slide the centromere along the midline to the point nearest to (x, y)."""
        mark = self.details.chromosome(mark_id)
        if mark is None:
            return
        if mark.midline and len(mark.midline) >= 2:
            pts = np.asarray(mark.midline, dtype=np.float64)
            k = int(np.argmin(np.hypot(pts[:, 0] - x, pts[:, 1] - y)))
            k = max(1, min(len(pts) - 2, k))
            mark.centromere = (float(pts[k, 0]), float(pts[k, 1]))
            if k > len(pts) - 1 - k:  # p stays the shorter arm
                mark.midline = mark.midline[::-1]
                mark.telomere_p, mark.telomere_q = mark.telomere_q, mark.telomere_p
        else:
            mark.centromere = (float(x), float(y))
        mark.check = "; ".join(c for c in mark.check.split("; ") if c and c != sm.CHECK_CENTROMERE)
        self.reclassify()
        self.dirty = True

    def move_end(self, mark_id: str, which: str, x: float, y: float) -> None:
        """Move a telomere: the midline is cut at the nearest point and goes straight to the new end."""
        mark = self.details.chromosome(mark_id)
        if mark is None or not mark.midline or len(mark.midline) < 3:
            if mark is not None:
                setattr(mark, which, (float(x), float(y)))
                self.reclassify()
                self.dirty = True
            return
        pts = [tuple(p) for p in mark.midline]
        cen = mark.centromere
        if which == "telomere_q":
            pts = pts[::-1]
        arr = np.asarray(pts)
        ci = int(np.argmin(np.hypot(arr[:, 0] - cen[0], arr[:, 1] - cen[1]))) if cen else len(pts) // 2
        k = int(np.argmin(np.hypot(arr[:ci, 0] - x, arr[:ci, 1] - y))) if ci > 0 else 0
        new = [(float(x), float(y))] + pts[k + 1:] if k + 1 < len(pts) else [(float(x), float(y))] + pts[-2:]
        if which == "telomere_q":
            new = new[::-1]
        mark.midline = [[round(a, 2), round(b, 2)] for a, b in new]
        mark.telomere_p, mark.telomere_q = tuple(mark.midline[0]), tuple(mark.midline[-1])
        self.reclassify()
        self.dirty = True

    def swap_mark_arms(self, mark_id: str) -> None:
        """p and q the other way round (the midline is turned too)."""
        mark = self.details.chromosome(mark_id)
        if mark is None:
            return
        if mark.midline:
            mark.midline = mark.midline[::-1]
        mark.telomere_p, mark.telomere_q = mark.telomere_q, mark.telomere_p
        self.reclassify()
        self.dirty = True

    def length_um(self, mark: ChromosomeMark) -> float | None:
        if not mark.midline:
            return None
        px = sm.polyline_length(mark.midline)
        return px * self.details.pixel_size_um if self.details.pixel_size_um else None

    def arms_px(self, mark: ChromosomeMark) -> tuple[float, float] | None:
        if not mark.midline or mark.centromere is None:
            return None
        total = sm.polyline_length(mark.midline)
        at = sg.arc_position(mark.midline, mark.centromere)
        return at, total - at

    def centromere_index(self, mark: ChromosomeMark) -> float | None:
        arms = self.arms_px(mark)
        if not arms or sum(arms) <= 0:
            return None
        return round(100.0 * min(arms) / sum(arms), 1)

    # ---------------------------------------------------------------- signals

    def _reassign_signals(self) -> None:
        labels = self.labels_map()
        for s in self.details.signals:
            if s.manual:
                continue
            if labels is not None and labels.any():
                k = sp.assign((s.x, s.y), labels)
                s.chromosome_id = self.details.chromosomes[k].id if k >= 0 else ""
            else:
                s.chromosome_id = ""
            if s.auto:
                checks = [c for c in s.check.split("; ") if c and c != CHECK_OFF]
                if not s.chromosome_id:
                    checks.append(CHECK_OFF)
                s.check = "; ".join(checks)
        self.reclassify()

    def detect_signals(self, row: ProbeRow, settings: sp.SpotSettings | None = None) -> int:
        """Replace the detected signals of one probe. Signals added by hand stay. Returns how many were found."""
        if not row.probe:
            return 0
        array = self.probe_array(row)
        if array is None:
            return 0
        settings = settings or self.spot_settings(row.probe)
        labels = self.labels_map()
        found = sp.find_spots(array, settings, labels if labels is not None and labels.any() else None)
        return self.apply_spots(row, found, settings)

    def apply_spots(self, row: ProbeRow, found: list, settings: sp.SpotSettings) -> int:
        """Store spots found by ``spots.find_spots`` for one probe (replacing its detected ones)."""
        self.set_spot_settings(row.probe, settings)
        labels = self.labels_map()
        self.details.signals = [s for s in self.details.signals if not (s.auto and s.probe == row.probe)]
        for spot in found:
            k = sp.assign((spot.x, spot.y), labels) if labels is not None and labels.any() else -1
            checks = []
            if spot.peak < 1.3:
                checks.append(CHECK_WEAK)
            if k < 0:
                checks.append(CHECK_OFF)
            sig = Signal(id=new_id(), probe=row.probe, x=round(spot.x, 2), y=round(spot.y, 2),
                         chromosome_id=self.details.chromosomes[k].id if k >= 0 else "",
                         strength=sp.strength(spot), auto=True, check="; ".join(checks), peak=round(spot.peak, 2))
            sg.classify(sig, self.details.chromosome(sig.chromosome_id))
            self.details.signals.append(sig)
        self.dirty = True
        return len(found)

    def detect_all_signals(self) -> dict[str, int]:
        return {r.probe: self.detect_signals(r) for r in self.probe_rows() if r.probe}

    # ---------------------------------------------------------------- review

    def to_check(self) -> list[tuple[str, str, str]]:
        """("mark" | "signal", id, reason) of everything still to look at, chromosomes first."""
        out = [("mark", m.id, m.check) for m in self.details.chromosomes if m.check]
        out += [("signal", s.id, s.check) for s in self.details.signals if s.check]
        return out

    def accept(self, kind: str, ident: str) -> None:
        item = self.details.chromosome(ident) if kind == "mark" else self.signal(ident)
        if item is not None:
            item.check = ""
            self.dirty = True

    def accept_all(self) -> None:
        for m in self.details.chromosomes:
            m.check = ""
        for s in self.details.signals:
            s.check = ""
        self.dirty = True

    def copy_settings_to_slide(self) -> int:
        """Detection settings and probe assignment of this cell to the other cells of its slide (on save)."""
        self._copy_to_slide = True
        self.dirty = True
        return 1
