"""One open multi-channel image: display changes, scale bar, saving. No Qt.

The Image tab calls this; it never touches the raw TIFFs. Every change goes to
``ChannelDisplay`` and is written by ``save()`` (log, record, merged preview).
"""

from __future__ import annotations

import shutil
from pathlib import Path

import numpy as np

from ..cyto_images import CytoImage, load_images, store_dir, update_image
from . import channels as ch
from .importer import channel_layers, natural_key, refresh_log, render_preview
from .lists import dye_color, dye_for_filename, load_lists
from .model import Channel, ChannelDisplay, ImageDetails, channel_path, image_details, new_id, set_image_details

NICE_LENGTHS_UM = (0.5, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000)
LENGTH_CHOICES = ("Auto", "1 µm", "2 µm", "5 µm", "10 µm", "20 µm", "50 µm", "100 µm")


def nice_length(target_um: float) -> float:
    """Largest round length (1, 2, 5, 10 ...) that is not longer than ``target_um``."""
    best = NICE_LENGTHS_UM[0]
    for value in NICE_LENGTHS_UM:
        if value <= target_um:
            best = value
    return float(best)


def parse_length(text: str) -> float:
    """"10 µm" -> 10.0; "Auto" or anything unreadable -> 0.0."""
    try:
        return float(str(text).lower().replace("µm", "").replace("um", "").strip())
    except ValueError:
        return 0.0


def tile(views: list[np.ndarray], columns: int = 2, gap: int = 6) -> np.ndarray:
    """Lay equal-sized (H, W, 3) uint8 views out on a dark grid."""
    if not views:
        raise ValueError("No views")
    h, w = views[0].shape[:2]
    columns = max(1, min(columns, len(views)))
    rows = (len(views) + columns - 1) // columns
    out = np.full((rows * h + (rows - 1) * gap, columns * w + (columns - 1) * gap, 3), 12, dtype=np.uint8)
    for n, view in enumerate(views):
        r, c = divmod(n, columns)
        out[r * (h + gap): r * (h + gap) + h, c * (w + gap): c * (w + gap) + w] = view
    return out


class ImageSession:
    """Display state of one image, with a per-session cache of the channel arrays."""

    def __init__(self, folder: Path, image: CytoImage) -> None:
        self.folder = Path(folder).expanduser().resolve()
        self.image = image
        self.details: ImageDetails = image_details(image)
        self.cache: dict = {}
        scale = image.extra.get("scale_bar") if isinstance(image.extra.get("scale_bar"), dict) else {}
        self.scale_show: bool = bool(scale.get("show", True))
        self.scale_length_um: float = float(scale.get("length_um") or 0.0)  # 0 = auto

    # ---------------------------------------------------------------- channels

    @property
    def channels(self) -> list[Channel]:
        return self.details.channels

    @property
    def is_multichannel(self) -> bool:
        return bool(self.details.channels)

    def label(self, index: int) -> str:
        c = self.channels[index]
        return f"{c.dye} · {c.probe}" if c.dye and c.probe else (c.dye or c.probe or f"Channel {index + 1}")

    def level_range(self, index: int) -> tuple[float, float]:
        """Slider range for the levels of one channel."""
        c = self.channels[index]
        top = float(2 ** min(c.bits, 16) - 1) if c.bits < 32 else 65535.0
        return 0.0, max(top, c.display.high)

    def set_visible(self, index: int, visible: bool) -> None:
        self.channels[index].display.visible = bool(visible)

    def set_color(self, index: int, color: str) -> None:
        self.channels[index].display.color = color

    def set_levels(self, index: int, low: float, high: float) -> None:
        d = self.channels[index].display
        d.low, d.high = float(low), float(max(high, low + 1))

    def set_gamma(self, index: int, gamma: float) -> None:
        self.channels[index].display.gamma = float(min(max(gamma, 0.1), 5.0))

    def auto_levels(self, index: int | None = None) -> None:
        """Auto levels for one channel, or for all when ``index`` is None."""
        layers = channel_layers(self.folder, self.details, self.cache)
        for i, (array, display) in enumerate(layers):
            if index is None or i == index:
                display.low, display.high = ch.auto_levels(array)
                display.gamma = 1.0

    def add_channel(self, path: Path, *, dye: str = "", probe: str = "", copy: bool = True) -> Channel:
        """Add a TIFF as one more channel. It must have the size of the others."""
        path = Path(path)
        data = ch.read_channel(path)
        grey = ch.to_grey(data.array)
        if self.channels:
            existing = ch.to_grey(channel_layers(self.folder, self.details, self.cache)[0][0])
            if existing.shape != grey.shape:
                raise ValueError(
                    f"{path.name} is {grey.shape[1]}×{grey.shape[0]}, the image is {existing.shape[1]}×{existing.shape[0]}"
                )
        lists = load_lists(self.folder)
        dye = dye or dye_for_filename(path.name, lists)[0] or ""
        cid = new_id()
        if copy:
            directory = store_dir(self.folder)
            directory.mkdir(parents=True, exist_ok=True)
            stored = directory / f"{self.image.id}_{cid[:8]}{path.suffix.lower()}"
            shutil.copy2(path, stored)
            file_ref = stored.name
        else:
            file_ref = str(path.resolve())
        low, high = ch.auto_levels(grey)
        channel = Channel(id=cid, file=file_ref, dye=dye, probe=probe, bits=data.bits, original_name=path.name,
                          display=ChannelDisplay(color=dye_color(lists, dye) if dye else "#ffffff", low=low, high=high))
        self.channels.append(channel)
        if not self.details.pixel_size_um and data.pixel_size_um:
            self.details.pixel_size_um = data.pixel_size_um
        self.cache.clear()
        return channel

    def remove_channel(self, index: int) -> None:
        """Remove a channel from the image. A copy in the store is deleted, other files are left."""
        channel = self.channels.pop(index)
        path = channel_path(self.folder, channel)
        if not Path(channel.file).is_absolute():
            try:
                path.unlink()
            except OSError:
                pass
        self.cache.clear()

    # ---------------------------------------------------------------- views

    def merged(self, max_side: int = 0) -> np.ndarray:
        return ch.composite(channel_layers(self.folder, self.details, self.cache), max_side=max_side)

    def side_by_side(self, max_side: int = 0, columns: int = 2) -> np.ndarray:
        """One panel per channel (in its colour) plus the merge, tiled."""
        layers = channel_layers(self.folder, self.details, self.cache)
        views = [ch.single_channel_view(a, d, max_side=max_side) for a, d in layers]
        if len(layers) > 1:
            views.append(ch.composite(layers, max_side=max_side))
        return tile(views, columns)

    def panel_labels(self) -> list[str]:
        labels = [self.label(i) for i in range(len(self.channels))]
        if len(labels) > 1:
            labels.append("Merge")
        return labels

    # ---------------------------------------------------------------- scale bar

    @property
    def pixel_size_um(self) -> float | None:
        return self.details.pixel_size_um

    def set_pixel_size(self, value: float | None) -> None:
        self.details.pixel_size_um = float(value) if value and value > 0 else None

    def scale_bar(self, image_width_px: int) -> tuple[float, float] | None:
        """(length in µm, length in image pixels) or None without a pixel size."""
        pixel = self.details.pixel_size_um
        if not pixel:
            return None
        length = self.scale_length_um or nice_length(image_width_px * pixel * 0.2)
        return length, length / pixel

    # ---------------------------------------------------------------- saving

    def save(self) -> None:
        """Write log, record and merged preview. Raw TIFFs stay untouched."""
        refresh_log(self.details)
        set_image_details(self.image, self.details)
        self.image.extra["scale_bar"] = {"show": self.scale_show, "length_um": self.scale_length_um}
        update_image(self.folder, self.image)
        if self.is_multichannel:
            render_preview(self.folder, self.image, cache=self.cache)

    def apply_to_slide(self) -> int:
        """Copy colours, levels and gamma of this image to the other cells of the slide (matched by dye).

        Returns how many other images were changed. Call ``save()`` first for this one.
        """
        slide_id = self.details.slide_id
        if not slide_id:
            return 0
        source = {c.dye: c.display for c in self.channels if c.dye}
        changed = 0
        for other in load_images(self.folder):
            if other.id == self.image.id:
                continue
            details = image_details(other)
            if details.slide_id != slide_id or not details.channels:
                continue
            touched = False
            for c in details.channels:
                d = source.get(c.dye)
                if d is not None:
                    c.display = ChannelDisplay(color=d.color, visible=d.visible, low=d.low, high=d.high, gamma=d.gamma)
                    touched = True
            if not touched:
                continue
            refresh_log(details)
            set_image_details(other, details)
            update_image(self.folder, other)
            try:
                render_preview(self.folder, other)
            except Exception:  # noqa: BLE001 - one unreadable cell must not stop the rest
                continue
            changed += 1
        return changed

    def slide_cell_count(self) -> int:
        slide_id = self.details.slide_id
        if not slide_id:
            return 0
        return sum(1 for o in load_images(self.folder) if image_details(o).slide_id == slide_id)


def neighbours(folder: Path, image_id: str) -> tuple[int, int, str, str]:
    """(position 1-based, total, previous id, next id) among the cells of the same slide.

    Images without a slide are compared with all images without a slide.
    """
    images = list(load_images(folder))
    current = next((i for i in images if i.id == image_id), None)
    if current is None:
        return 0, 0, "", ""
    slide_id = image_details(current).slide_id
    group = [i for i in images if image_details(i).slide_id == slide_id]
    group.sort(key=lambda i: natural_key(i.original_name or i.id))
    ids = [i.id for i in group]
    pos = ids.index(image_id)
    prev_id = ids[pos - 1] if pos > 0 else ""
    next_id = ids[pos + 1] if pos + 1 < len(ids) else ""
    return pos + 1, len(ids), prev_id, next_id
