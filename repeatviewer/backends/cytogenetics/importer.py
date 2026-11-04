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
from .model import Channel, ChannelDisplay, ImageDetails, channel_path, image_details, new_id, ordered_channels, set_image_details

PREVIEW_SIDE = 2048


def natural_key(text: str):
    return [int(t) if t.isdigit() else t.lower() for t in re.split(r"(\d+)", text)]


@dataclass
class CellGroup:
    key: str  # file name without the channel word ("s12_c03"), the cell folder ("cell1"), or the file stem
    files: list[tuple] = field(default_factory=list)  # (path, dye) or (path, dye, plane index)
    colours: list[str] = field(default_factory=list)  # colour per entry; "" = from the dye

    @property
    def dyes(self) -> list[str]:
        return [f[1] for f in self.files]

    def entries(self) -> list[tuple[Path, str, int]]:
        return [(f[0], f[1], f[2] if len(f) > 2 else -1) for f in self.files]

    @property
    def size(self) -> int:
        return len(self.files)


def group_channel_files(paths: list[Path], lists: Lists) -> tuple[list[CellGroup], list[Path]]:
    """Cells from file names: the name without its dye word is the cell.

    ``s12_c03_DAPI.tif`` + ``s12_c03_FITC.tif`` -> cell ``s12_c03`` with two channels.
    Files without a known dye word are returned separately.
    """
    groups: dict[str, CellGroup] = {}
    unmatched: list[Path] = []
    for path in sorted((Path(p) for p in paths), key=lambda p: natural_key(p.name)):
        if not ch.is_image(path):
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


def plan_cells(paths: list[Path], lists: Lists) -> tuple[list[CellGroup], list[str]]:
    """Every TIFF goes somewhere; nothing needs a dye name. Returns (cells, notes).

    - A multichannel TIFF (ImageJ, OME, several pages) is one cell, split into its channels.
    - Files named with a dye word (s12_c03_DAPI.tif) are grouped by the rest of the name, as before.
    - Other files (1.tif, 2.tif …) are grouped by folder: one folder = one cell, one file = one channel,
      named Channel 1, 2 … in the file order; dye and probe are set later in the Image tab.
    - An RGB merge on its own becomes one cell with its red, green and blue planes as channels.
    """
    notes: list[str] = []
    singles: list[Path] = []
    multi: list[CellGroup] = []
    for path in sorted((Path(p) for p in paths), key=lambda p: natural_key(str(p))):
        if not ch.is_image(path):
            continue
        try:
            info = ch.inspect_tiff(path)
        except Exception as exc:  # noqa: BLE001 - unreadable file: report, keep going
            notes.append(f"{path.name}: cannot be read ({exc})")
            continue
        if info.kind == "channels":
            group = CellGroup(key=path.stem)
            for i in range(info.count):
                dye = dye_for_filename(info.names[i], lists)[0] if info.names[i] else ""
                group.files.append((path, dye, i))
                group.colours.append("" if dye else info.colours[i])
            multi.append(group)
        else:
            singles.append(path)
    named, rest = group_channel_files(singles, lists)
    by_folder: dict[Path, list[Path]] = {}
    for path in rest:
        by_folder.setdefault(path.parent, []).append(path)
    for parent, files in by_folder.items():
        files.sort(key=lambda p: natural_key(p.name))
        group = CellGroup(key=parent.name)
        n = 0
        for path in files:
            try:
                info = ch.inspect_tiff(path)
            except Exception:  # noqa: BLE001 - reported above
                continue
            planes = range(3) if info.kind == "rgb" and len(files) == 1 else [-1]
            for k in planes:
                group.files.append((path, "", k))
                group.colours.append(info.colours[k] if k >= 0 else ch.UNKNOWN_COLOURS[n % len(ch.UNKNOWN_COLOURS)])
                n += 1
        if group.files:
            multi.append(group)
    cells = named + multi
    # cells with the same name from different folders: add the folder to tell them apart
    seen: dict[str, int] = {}
    for g in cells:
        seen[g.key] = seen.get(g.key, 0) + 1
    for g in cells:
        if seen[g.key] > 1:
            g.key = f"{g.entries()[0][0].parent.name}_{g.key}"
    return sorted(cells, key=lambda g: natural_key(g.key)), notes


def channel_layers(folder: Path, details: ImageDetails, cache: dict | None = None):
    """[(array, display)] for the image's channels, reading each file once per cache."""
    layers = []
    for channel in details.channels:
        key = f"{channel_path(folder, channel)}#{channel.index}"
        if cache is not None and key in cache:
            array = cache[key]
        else:
            array = ch.read_plane(channel_path(folder, channel), channel.index).array
            if cache is not None:
                cache[key] = array
        layers.append((array, channel.display))
    return layers


def refresh_log(details: ImageDetails) -> None:
    named = [(c.dye or f"Channel {details.channels.index(c) + 1}", c.bits, c.display) for c in ordered_channels(details)]
    details.adjustment_log = ch.adjustment_log(named)


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
            stored_as: dict[Path, str] = {}
            for n, (path, dye, index) in enumerate(group.entries()):
                data = ch.read_plane(path, index)
                grey = ch.to_grey(data.array)
                low, high = ch.auto_levels(grey)
                cid = new_id()
                if path not in stored_as:
                    if copy:
                        stored = directory / f"{image_id}_{cid[:8]}{path.suffix.lower()}"
                        shutil.copy2(path, stored)
                        copied.append(stored)
                        stored_as[path] = stored.name
                    else:
                        stored_as[path] = str(path.resolve())
                file_ref = stored_as[path]
                colour = group.colours[n] if n < len(group.colours) and group.colours[n] else ""
                colour = colour or (dye_color(lists, dye) if dye else ch.UNKNOWN_COLOURS[n % len(ch.UNKNOWN_COLOURS)])
                display = ChannelDisplay(color=colour, low=low, high=high)
                label = path.name if index < 0 else f"{path.name} · plane {index + 1}"
                channels.append(Channel(id=cid, file=file_ref, dye=dye, probe=probes_by_dye.get(dye, "") if dye else "",
                                        bits=data.bits, original_name=label, display=display, index=index))
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


# ---------------------------------------------------------------------------
# One slide session = one request (the Import dialog builds it)
# ---------------------------------------------------------------------------

from .lists import add_entry, ensure_used_values, entry, load_lists, save_lists  # noqa: E402
from .model import add_slide, load_slides, save_slides  # noqa: E402

SLIDE_FIELDS = ("date", "species", "accession", "technique", "counterstain")


@dataclass
class ImportRequest:
    folder: Path                 # the data folder (where the store lives)
    source: Path                 # folder with the microscope TIFFs
    recursive: bool = False
    slide_name: str = ""
    slide_values: dict = field(default_factory=dict)   # SLIDE_FIELDS
    stage: str = ""
    probes_by_dye: dict = field(default_factory=dict)  # dye -> probe name
    copy: bool = True
    pixel_size_um: float | None = None


def scan_tiffs(source: Path, recursive: bool = False) -> list[Path]:
    source = Path(source).expanduser()
    if not source.is_dir():
        return []
    walker = source.rglob("*") if recursive else source.iterdir()
    return sorted((p for p in walker if p.is_file() and ch.is_image(p) and not p.name.startswith(".")),
                  key=lambda p: natural_key(str(p)))


def targets_for_probes(lists: Lists, probes_by_dye: dict) -> list[ImageTarget]:
    """Sequence links for the probes that have a FASTA in the Probes list (keeps right-click on sequences working)."""
    out = []
    for dye, probe in probes_by_dye.items():
        e = entry(lists, "probes", probe) if probe else None
        seq = str((e or {}).get("sequence") or "").strip()
        if seq:
            out.append(ImageTarget(seq, probe=probe, dye=dye, color=dye_color(lists, dye, "")))
    return out


def run_request(req: ImportRequest, *, progress=None):
    """Create the slide and its images. Returns (slide or None, new images, error lines).

    The slide is dropped again when no image could be imported.
    """
    folder = Path(req.folder).expanduser().resolve()
    lists = load_lists(folder)
    files = scan_tiffs(req.source, req.recursive)
    groups, errors = plan_cells(files, lists)
    if not groups:
        return None, [], errors or ["No picture files found (TIFF, JPG, PNG or BMP)."]
    values = {k: str(req.slide_values.get(k, "")).strip() for k in SLIDE_FIELDS}
    slide = add_slide(folder, req.slide_name or Path(req.source).name, **values)
    metadata = {k: v for k, v in values.items() if v}
    if req.stage:
        metadata["stage"] = req.stage
    images, errs = import_cells(
        folder, groups, lists=lists, slide_id=slide.id, metadata=metadata,
        targets=targets_for_probes(lists, req.probes_by_dye), probes_by_dye=req.probes_by_dye,
        copy=req.copy, pixel_size_um=req.pixel_size_um, progress=progress,
    )
    errors += errs
    if not images:
        save_slides(folder, [s for s in load_slides(folder) if s.id != slide.id])
        return None, [], errors
    lists = load_lists(folder)
    for probe in req.probes_by_dye.values():
        if probe:
            add_entry(lists, "probes", probe)
    if ensure_used_values(folder, lists) or req.probes_by_dye:
        save_lists(folder, lists)
    return slide, images, errors


def add_channel_image(folder: Path, path: Path, *, copy: bool = True) -> CytoImage | None:
    """Add one multichannel TIFF as a channel image (no slide). Returns the new image, or None."""
    folder = Path(folder).expanduser().resolve()
    lists = load_lists(folder)
    cells, _notes = plan_cells([Path(path)], lists)
    if not cells:
        return None
    images, errors = import_cells(folder, cells[:1], lists=lists, copy=copy)
    if errors:
        raise ValueError(errors[0])
    return images[0] if images else None
