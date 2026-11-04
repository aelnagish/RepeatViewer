"""Figure builder: panels with letters, labels, arrowheads, scale bar and a written legend. No Qt.

A figure is kept in the store as ``<name>.rvfig`` (JSON). Positions on the page are in mm;
positions inside a panel (arrowheads, crop) are fractions of the source image (0..1), so
they survive any change of size or resolution.
"""

from __future__ import annotations

import json
import os
import re
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from ..cyto_images import CytoImage, _color_word, _join_and, load_images, store_dir
from . import channels as ch
from .importer import channel_layers, natural_key
from .karyotype import list_karyotypes
from .model import ImageDetails, image_details, load_slides, new_id, ordered_channels  # noqa: F401 (re-export)

SUFFIX = ".rvfig"
WIDTHS = (("Two columns (174 mm)", 174.0), ("One column (84 mm)", 84.0), ("Custom", 0.0))
LAYOUTS = (("2 × 2", "2x2"), ("4 in a row", "4row"), ("3 + idiogram", "3+idiogram"), ("One row per cell", "auto"))
FONTS = (("Arial 8 pt", 8.0), ("Arial 10 pt", 10.0))
VIEW_MERGE = "merge"
VIEW_IDIOGRAM = "idiogram"
LEGEND_GAP_MM = 3.0


@dataclass
class Arrow:
    x: float  # fraction of the source image width
    y: float
    angle: float = 135.0  # direction the arrow points to, degrees (screen: 0 = right, 90 = down; 135 = down-left)


@dataclass
class Panel:
    id: str
    image_id: str = ""
    view: str = VIEW_MERGE  # "merge", a channel id, or "idiogram"
    label: str = ""  # "" = automatic (dye, probe or "Merge")
    arrows: list[Arrow] = field(default_factory=list)
    crop: tuple[float, float, float, float] | None = None  # x0, y0, x1, y1 fractions; None = whole image

    @property
    def is_idiogram(self) -> bool:
        return self.view == VIEW_IDIOGRAM


@dataclass
class Figure:
    name: str = "Figure 1"
    width_mm: float = 174.0
    custom_mm: float = 120.0
    layout: str = "2x2"
    spacing: int = 2  # 0..10, half millimetres
    font_pt: float = 8.0
    panels: list[Panel] = field(default_factory=list)
    same_crop: bool = True
    crop: tuple[float, float, float, float] | None = None
    letters: bool = True
    channel_labels: bool = True
    scale_bar: bool = True
    scale_um: float = 10.0
    arrowheads: bool = True
    chromosome_numbers: bool = False
    number: str = "Figure 1."  # bold start of the legend
    legend: str = ""  # "" = written automatically; otherwise the user's text (simple HTML)

    @property
    def page_width(self) -> float:
        return self.custom_mm if self.width_mm <= 0 else self.width_mm

    @property
    def gap_mm(self) -> float:
        return self.spacing * 0.5

    def panel(self, panel_id: str) -> Panel | None:
        return next((p for p in self.panels if p.id == panel_id), None)

    def crop_of(self, panel: Panel) -> tuple[float, float, float, float]:
        crop = self.crop if self.same_crop else panel.crop
        return crop or (0.0, 0.0, 1.0, 1.0)

    def to_dict(self) -> dict:
        return asdict(self)

    @classmethod
    def from_dict(cls, raw) -> "Figure":
        raw = raw if isinstance(raw, dict) else {}
        fig = cls()
        for key, value in raw.items():
            if key == "panels" or not hasattr(fig, key):
                continue
            if key == "crop":
                fig.crop = _crop(value)
                continue
            try:
                current = getattr(fig, key)
                setattr(fig, key, type(current)(value) if current is not None and not isinstance(current, str) else str(value))
            except (TypeError, ValueError):
                continue
        for p in raw.get("panels") or []:
            if not isinstance(p, dict) or not p.get("id"):
                continue
            arrows = []
            for a in p.get("arrows") or []:
                try:
                    arrows.append(Arrow(float(a["x"]), float(a["y"]), float(a.get("angle", 135.0))))
                except (KeyError, TypeError, ValueError):
                    continue
            fig.panels.append(Panel(id=str(p["id"]), image_id=str(p.get("image_id") or ""),
                                    view=str(p.get("view") or VIEW_MERGE), label=str(p.get("label") or ""),
                                    arrows=arrows, crop=_crop(p.get("crop"))))
        return fig


def _crop(value) -> tuple[float, float, float, float] | None:
    try:
        x0, y0, x1, y1 = (max(0.0, min(1.0, float(v))) for v in value)
    except (TypeError, ValueError):
        return None
    x0, x1 = min(x0, x1), max(x0, x1)
    y0, y1 = min(y0, y1), max(y0, y1)
    if x1 - x0 < 0.02 or y1 - y0 < 0.02:
        return None
    return (x0, y0, x1, y1)


def normalise_crop(x0: float, y0: float, x1: float, y1: float) -> tuple[float, float, float, float] | None:
    return _crop((x0, y0, x1, y1))


# ---------------------------------------------------------------------------
# Store
# ---------------------------------------------------------------------------


def _safe(name: str) -> str:
    return re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip()).strip("._") or "Figure"


def figure_path(folder: Path, name: str) -> Path:
    return store_dir(Path(folder).expanduser().resolve()) / f"{_safe(name)}{SUFFIX}"


def list_figures(folder: Path) -> list[Figure]:
    directory = store_dir(Path(folder).expanduser().resolve())
    out = []
    for path in sorted(directory.glob(f"*{SUFFIX}"), key=lambda p: natural_key(p.name)) if directory.is_dir() else []:
        try:
            out.append(Figure.from_dict(json.loads(path.read_text(encoding="utf-8"))))
        except (OSError, ValueError):
            continue
    return out


def save_figure(folder: Path, fig: Figure, *, old_name: str = "") -> Path:
    path = figure_path(folder, fig.name)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(".tmp")
    tmp.write_text(json.dumps(fig.to_dict(), indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    os.replace(tmp, path)
    if old_name and _safe(old_name) != _safe(fig.name):
        try:
            figure_path(folder, old_name).unlink()
        except OSError:
            pass
    return path


def delete_figure(folder: Path, name: str) -> bool:
    try:
        figure_path(folder, name).unlink()
        return True
    except OSError:
        return False


def new_figure_name(folder: Path) -> str:
    used = {f.name for f in list_figures(folder)}
    n = 1
    while f"Figure {n}" in used:
        n += 1
    return f"Figure {n}"


# ---------------------------------------------------------------------------
# Panels
# ---------------------------------------------------------------------------


def letter(index: int) -> str:
    out = ""
    index += 1
    while index:
        index, rem = divmod(index - 1, 26)
        out = chr(65 + rem) + out
    return out


def panels_for_cell(image: CytoImage) -> list[Panel]:
    """One panel per channel plus the merge (as in the mockup); one panel for a plain picture."""
    details = image_details(image)
    if not details.channels:
        return [Panel(id=new_id(), image_id=image.id)]
    panels = [Panel(id=new_id(), image_id=image.id, view=c.id) for c in ordered_channels(details)]
    panels.append(Panel(id=new_id(), image_id=image.id, view=VIEW_MERGE))
    return panels


def add_cell(fig: Figure, image: CytoImage) -> list[Panel]:
    new = panels_for_cell(image)
    fig.panels.extend(new)
    if len(fig.panels) > 4 and fig.layout == "2x2":
        fig.layout = "auto"
    return new


def add_idiogram(fig: Figure) -> Panel:
    p = Panel(id=new_id(), view=VIEW_IDIOGRAM)
    fig.panels.append(p)
    if fig.layout == "2x2" and len(fig.panels) == 4:
        fig.layout = "3+idiogram"
    return p


def move_panel(fig: Figure, panel_id: str, step: int) -> None:
    i = next((k for k, p in enumerate(fig.panels) if p.id == panel_id), None)
    if i is None:
        return
    j = max(0, min(len(fig.panels) - 1, i + step))
    fig.panels.insert(j, fig.panels.pop(i))


def remove_panel(fig: Figure, panel_id: str) -> None:
    fig.panels = [p for p in fig.panels if p.id != panel_id]


def auto_label(panel: Panel, details: ImageDetails | None) -> str:
    if panel.label:
        return panel.label
    if panel.is_idiogram:
        return ""
    if panel.view == VIEW_MERGE or details is None:
        return "Merge" if details is not None and details.channels else ""
    c = next((c for c in details.channels if c.id == panel.view), None)
    if c is None:
        return ""
    return c.probe if c.probe and c.dye.upper() != "DAPI" else (c.dye or c.probe)


def label_colour(panel: Panel, details: ImageDetails | None) -> str:
    """Colour of the channel label: the dye colour, lightened for dark blues (mockup: DAPI #8EA0FF)."""
    if panel.view == VIEW_MERGE or details is None or panel.is_idiogram:
        return "#FFFFFF"
    c = next((c for c in details.channels if c.id == panel.view), None)
    if c is None:
        return "#FFFFFF"
    r, g, b = ch.hex_to_rgb(c.display.color)
    if 0.299 * r + 0.587 * g + 0.114 * b < 0.45:  # too dark on black: mix with white
        r, g, b = (0.55 * r + 0.45, 0.55 * g + 0.45, 0.55 * b + 0.45)
    return "#{:02X}{:02X}{:02X}".format(int(r * 255), int(g * 255), int(b * 255))


def source_text(panel: Panel, image: CytoImage | None, slide_name: str = "") -> str:
    """"Slide 12 cell 03 · FITC" for the panel list."""
    if panel.is_idiogram:
        return "Idiogram (in situ vs in silico)"
    if image is None:
        return "missing image"
    details = image_details(image)
    m = re.search(r"(\d+)\D*$", image.original_name or "")
    cell = f"cell {m.group(1).zfill(2)}" if m else (image.original_name or image.id)
    head = f"{slide_name} {cell}".strip()
    if panel.view == VIEW_MERGE:
        return f"{head} · merge" if details.channels else head
    c = next((c for c in details.channels if c.id == panel.view), None)
    return f"{head} · {c.dye or c.probe}" if c else f"{head} · missing channel"


def render_panel(folder: Path, image: CytoImage, panel: Panel, crop: tuple[float, float, float, float],
                 cache: dict | None = None, max_side: int = 0) -> np.ndarray:
    """RGB uint8 of a panel (cropped). Plain pictures need Qt to read: returns an empty array for them."""
    details = image_details(image)
    if not details.channels:
        return np.zeros((0, 0, 3), dtype=np.uint8)
    layers = channel_layers(folder, details, cache)
    if panel.view == VIEW_MERGE:
        rgb = ch.composite(layers, max_side=max_side)
    else:
        idx = next((i for i, c in enumerate(details.channels) if c.id == panel.view), None)
        if idx is None:
            raise ValueError("The channel of this panel was removed")
        rgb = ch.single_channel_view(layers[idx][0], layers[idx][1], max_side=max_side)
    return crop_array(rgb, crop)


def crop_array(rgb: np.ndarray, crop: tuple[float, float, float, float]) -> np.ndarray:
    h, w = rgb.shape[:2]
    x0, y0, x1, y1 = crop
    a, b = int(round(x0 * w)), int(round(x1 * w))
    c, d = int(round(y0 * h)), int(round(y1 * h))
    return rgb[c:max(d, c + 1), a:max(b, a + 1)]


def scale_bar_fraction(details: ImageDetails, crop: tuple[float, float, float, float], full_width_px: int,
                       length_um: float) -> float | None:
    """Scale bar length as a fraction of the panel width, or None without a pixel size."""
    if not details.pixel_size_um or full_width_px <= 0:
        return None
    width_um = full_width_px * (crop[2] - crop[0]) * details.pixel_size_um
    return length_um / width_um if width_um > 0 else None


# ---------------------------------------------------------------------------
# Layout (mm)
# ---------------------------------------------------------------------------


@dataclass
class Box:
    x: float
    y: float
    w: float
    h: float


def columns_for(fig: Figure) -> int:
    n = len([p for p in fig.panels if not p.is_idiogram])
    if fig.layout == "2x2":
        return 2
    if fig.layout == "4row":
        return 4
    if fig.layout == "3+idiogram":
        return 3
    # one row per cell: as many columns as panels of the first cell (at most 5)
    first = fig.panels[0].image_id if fig.panels else ""
    per_cell = sum(1 for p in fig.panels if p.image_id == first and not p.is_idiogram) or 2
    return max(1, min(5, per_cell, n or 1))


def layout(fig: Figure, aspects: dict[str, float], idiogram_aspect: float = 0.45) -> tuple[list[tuple[Panel, Box]], float]:
    """Boxes of all panels and the height used, both in mm. ``aspects`` = height / width per panel id."""
    width = fig.page_width
    gap = fig.gap_mm
    cols = columns_for(fig)
    cell_w = (width - gap * (cols - 1)) / cols
    boxes: list[tuple[Panel, Box]] = []
    y = 0.0
    row: list[Panel] = []

    def flush(row_panels: list[Panel]) -> None:
        nonlocal y
        if not row_panels:
            return
        h = max(cell_w * aspects.get(p.id, 0.75) for p in row_panels)
        for k, p in enumerate(row_panels):
            boxes.append((p, Box(k * (cell_w + gap), y, cell_w, h)))
        y += h + gap

    for p in fig.panels:
        if p.is_idiogram:
            flush(row)
            row = []
            h = width * idiogram_aspect
            boxes.append((p, Box(0.0, y, width, h)))
            y += h + gap
            continue
        row.append(p)
        if len(row) == cols:
            flush(row)
            row = []
    flush(row)
    return boxes, max(0.0, y - gap)


# ---------------------------------------------------------------------------
# Legend and methods
# ---------------------------------------------------------------------------


def _images_of(folder: Path, fig: Figure) -> list[CytoImage]:
    ids = []
    for p in fig.panels:
        if p.image_id and p.image_id not in ids:
            ids.append(p.image_id)
    by_id = {i.id: i for i in load_images(folder)}
    return [by_id[i] for i in ids if i in by_id]


def legend_html(folder: Path, fig: Figure) -> str:
    """Legend written from probes, dyes, stage, species, karyotype and the adjustment log."""
    images = _images_of(folder, fig)
    number = f"<b>{_esc(fig.number)}</b> " if fig.number.strip() else ""
    if not images:
        return number + ("Idiogram of the FISH signals compared with the assembly." if any(p.is_idiogram for p in fig.panels) else "")
    first = images[0]
    details = image_details(first)
    meta = first.metadata
    slides = {s.id: s for s in load_slides(folder)}
    slide = slides.get(details.slide_id)
    technique = meta.get("technique") or (slide.technique if slide else "") or "FISH"
    species = meta.get("species") or (slide.species if slide else "")
    stage = (meta.get("stage") or "").strip()
    counter = meta.get("counterstain") or (slide.counterstain if slide else "")
    probes, counter_colour = [], ""
    for c in ordered_channels(details):
        word = _color_word(c.display.color)
        if c.dye.upper() in ("DAPI", "HOECHST", "PI") or (counter and c.dye.lower() == counter.lower()):
            counter, counter_colour = counter or c.dye, word
            continue
        if c.probe:
            probes.append(f"{_esc(c.probe)} ({', '.join(x for x in (word, c.dye) if x)})")
    if not probes:
        for t in first.targets:
            probes.append(_esc(t.probe_label))
    parts = [f"{_esc(technique)} of {_join_and(probes)}" if probes else _esc(technique)]
    subject = ""
    if stage:
        subject = stage[0].lower() + stage[1:]
        if stage.lower() not in ("interphase", "extended fibre") and "chromosome" not in subject:
            subject += " chromosomes"
    if species:
        subject = f"{subject} of <i>{_esc(species)}</i>" if subject else f"<i>{_esc(species)}</i>"
        karyo = next((k for k in list_karyotypes(folder) if k.name == details.karyotype or k.species == species), None)
        if karyo is not None and karyo.chromosomes:
            subject += f" ({karyo.formula})"
    if subject:
        parts[0] += f" on {subject}"
    text = parts[0]
    if counter:
        text += f", counterstained with {_esc(counter)}" + (f" ({counter_colour})" if counter_colour else "")
    text += ". "
    # panel list
    names = []
    for i, p in enumerate(fig.panels):
        img = next((im for im in images if im.id == p.image_id), None)
        lab = auto_label(p, image_details(img) if img else None) if not p.is_idiogram else "idiogram of the signals compared with the assembly"
        if lab:
            names.append(f"({letter(i)}) {_esc('merge' if lab == 'Merge' else lab)}")
    if len(names) > 1:
        text += ", ".join(names)
        arrows = [letter(i) for i, p in enumerate(fig.panels) if p.arrows] if fig.arrowheads else []
        if arrows:
            text += "; arrowheads mark the signals discussed in the text"
        text += ". "
    # adjustments
    channel_sets = [image_details(im).channels for im in images if image_details(im).channels]
    if channel_sets:
        sentence = ch.legend_sentence([(c.dye, c.bits, c.display) for c in channel_sets[0]])
        same = all([(c.dye, c.display.low, c.display.high, c.display.gamma) for c in cs] ==
                   [(c.dye, c.display.low, c.display.high, c.display.gamma) for c in channel_sets[0]] for cs in channel_sets)
        if len(channel_sets) > 1 and same:
            sentence = sentence.rstrip(".") + ", identically for all cells."
        text += sentence + " "
    if fig.scale_bar and any(image_details(im).pixel_size_um for im in images):
        text += f"Scale bar {fig.scale_um:g} µm."
    return (number + text).strip()


def methods_text(folder: Path, fig: Figure) -> str:
    """Methods paragraph: capture, merge, every level and gamma setting."""
    images = _images_of(folder, fig)
    with_channels = [im for im in images if image_details(im).channels]
    if not with_channels:
        return "Images were prepared in RepeatViewer; no channel adjustments were made."
    d = image_details(with_channels[0])
    dyes = ", ".join(c.dye for c in ordered_channels(d) if c.dye)
    bits = sorted({c.bits for c in d.channels})
    text = (f"Each channel ({dyes}) was captured as a separate {'/'.join(str(b) for b in bits)}-bit image and the channels "
            f"were merged additively in RepeatViewer. Raw image files were not changed. ")
    lines = []
    for c in ordered_channels(d):
        dd = c.display
        level = f"display range {dd.low:g} to {dd.high:g}" if dd.is_set else "automatic display range"
        gamma = f", gamma {dd.gamma:g}" if abs(dd.gamma - 1.0) > 1e-6 else ""
        lines.append(f"{c.dye or 'channel'}: {level}{gamma}")
    text += "Display settings: " + "; ".join(lines) + ". "
    if len(with_channels) > 1:
        text += f"The same settings were used for all {len(with_channels)} cells shown. "
    if d.pixel_size_um:
        text += f"Pixel size {d.pixel_size_um:g} µm."
    return text.strip()


def _esc(text: str) -> str:
    return str(text).replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")


def html_to_text(html: str) -> str:
    text = re.sub(r"<br\s*/?>", "\n", html)
    text = re.sub(r"<[^>]+>", "", text)
    return text.replace("&lt;", "<").replace("&gt;", ">").replace("&amp;", "&").strip()


def mm_to_px(mm: float, dpi: float) -> int:
    return int(round(mm / 25.4 * dpi))
