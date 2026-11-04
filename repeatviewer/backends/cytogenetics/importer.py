"""Import a slide session: channel TIFFs grouped into cells, copied into the store."""

from __future__ import annotations

import re
import shutil
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path

from ..cyto_images import CytoImage, ImageTarget, image_file_path, load_images, save_images, store_dir
from . import channels as ch
from .lists import Lists, dye_color, dye_for_filename, filename_tokens
from .model import Channel, ChannelDisplay, ImageDetails, image_details, new_id, set_image_details, channel_path

PREVIEW_SIDE = 2048


def natural_key(text: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", text)]


@dataclass
class CellGroup:
    key: str  # file name without the channel word, e.g. "s12_c03"
    files: list[tuple[Path, str]] = field(default_factory=list)  # (path, dye)

    @property
    def dyes(self) -> list[str]:
        return [d for _, d in self.files]


def group_channel_files(paths: list[Path], lists: Lists) -> tuple[list[CellGroup], list[Path]]:
    """Cells from file names: the name without its dye word is the cell.

    ``s12_c03_DAPI.tif`` + ``s12_c03_FITC.tif`` -> cell ``s12_c03`` with two channels.
    Files without a known dye word are returned separately.
    """
    groups: dict[str, CellGroup] = {}
    unmatched: list[Path] = []
    for path in sorted((Path(p) for p in paths), key=lambda p: natural_key(p.name)):
        if not ch.is_tiff(path):
            unmatched.append(path)
            continue
        dye, index = dye_for_filename(path.name, lists)
        if not dye:
            unmatched.append(path)
            continue
        tokens = filename_tokens(path.name)
        key = "_".join(tokens[:index] + tokens[index + 1:]) or path.parent.name
        group = groups.setdefault(key, CellGroup(key=key))
        if dye in group.dyes:  # same dye twice in one cell: keep them apart
            unmatched.append(path)
            continue
        group.files.append((path, dye))
    ordered = sorted(groups.values(), key=lambda g: natural_key(g.key))
    return ordered, unmatched


def channel_layers(folder: Path, details: ImageDetails, cache: dict | None = None):
    """[(array, display)] for the image's channels, reading each file once per cache."""
    layers = []
    for channel in details.channels:
        key = str(channel_path(folder, channel))
        if cache is not None and key in cache:
            array = cache[key]
        else:
            array = ch.read_channel(channel_path(folder, channel)).array
            if cache is not None:
                cache[key] = array
        layers.append((array, channel.display))
    return layers


def refresh_log(details: ImageDetails) -> None:
    details.adjustment_log = ch.adjustment_log([(c.dye, c.bits, c.display) for c in details.channels])


def render_preview(folder: Path, image: CytoImage, *, cache: dict | None = None) -> Path:
    """Write the merged preview PNG (``image.file``) from the channels and their display."""
    details = image_details(image)
    rgb = ch.composite(channel_layers(folder, details, cache), max_side=PREVIEW_SIDE)
    path = image_file_path(folder, image)
    ch.write_png(path, rgb)
    return path


def import_cells(
    folder: Path,
    groups: list[CellGroup],
    *,
    lists: Lists,
    slide_id: str = "",
    metadata: dict[str, str] | None = None,
    targets: list[ImageTarget] | None = None,
    probes_by_dye: dict[str, str] | None = None,
    copy: bool = True,
    pixel_size_um: float | None = None,
    progress=None,
) -> tuple[list[CytoImage], list[str]]:
    """Create one image record per cell. Returns (new images, error lines).

    ``pixel_size_um`` overrides the value read from the files. ``progress(done, total)``
    is called after each cell; return False from it to stop.
    """
    folder = Path(folder).expanduser().resolve()
    directory = store_dir(folder)
    directory.mkdir(parents=True, exist_ok=True)
    probes_by_dye = probes_by_dye or {}
    created: list[CytoImage] = []
    errors: list[str] = []
    for done, group in enumerate(groups, start=1):
        image_id = new_id()
        copied: list[Path] = []
        try:
            channels: list[Channel] = []
            layers = []
            pixel = pixel_size_um
            for path, dye in group.files:
                data = ch.read_channel(path)
                grey = ch.to_grey(data.array)
                low, high = ch.auto_levels(grey)
                cid = new_id()
                if copy:
                    stored = directory / f"{image_id}_{cid[:8]}{path.suffix.lower()}"
                    shutil.copy2(path, stored)
                    copied.append(stored)
                    file_ref = stored.name
                else:
                    file_ref = str(path.resolve())
                display = ChannelDisplay(color=dye_color(lists, dye), low=low, high=high)
                channels.append(Channel(id=cid, file=file_ref, dye=dye, probe=probes_by_dye.get(dye, ""),
                                        bits=data.bits, original_name=path.name, display=display))
                layers.append((grey, display))
                if pixel is None and data.pixel_size_um:
                    pixel = data.pixel_size_um
            details = ImageDetails(slide_id=slide_id, pixel_size_um=pixel, channels=channels)
            refresh_log(details)
            image = CytoImage(
                id=image_id,
                file=f"{image_id}.png",
                original_name=group.key,
                added=datetime.now().isoformat(timespec="seconds"),
                targets=list(targets or []),
                metadata={k: str(v) for k, v in (metadata or {}).items() if str(v).strip()},
            )
            set_image_details(image, details)
            ch.write_png(image_file_path(folder, image), ch.composite(layers, max_side=PREVIEW_SIDE))
            created.append(image)
        except Exception as exc:  # one bad cell must not stop the session
            for p in copied:
                try:
                    p.unlink()
                except OSError:
                    pass
            errors.append(f"{group.key}: {exc}")
        if progress is not None and progress(done, len(groups)) is False:
            break
    if created:
        save_images(folder, list(load_images(folder)) + created)
    return created, errors
