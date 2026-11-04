"""Channel images: read TIFFs, set levels, compose a colour view, write a preview.

Raw files are never changed: levels, gamma and colours are display settings,
applied here on a copy. Everything is numpy; no Qt.
"""

from __future__ import annotations

import re
import struct
import zlib
from dataclasses import dataclass
from pathlib import Path

import numpy as np

from .model import ChannelDisplay

TIFF_SUFFIXES = (".tif", ".tiff")
# Pictures that can be imported as channels or merges (JPG and the rest are read without tifffile)
IMAGE_SUFFIXES = TIFF_SUFFIXES + (".jpg", ".jpeg", ".png", ".bmp", ".gif", ".webp")
DEFAULT_SATURATED = 0.1  # percent of pixels clipped at each end by Auto levels


@dataclass
class ChannelData:
    array: np.ndarray  # 2D (grey) or 3D (H, W, 3/4 for an RGB file)
    bits: int
    pixel_size_um: float | None


def is_tiff(path: Path) -> bool:
    return Path(path).suffix.lower() in TIFF_SUFFIXES


def is_image(path: Path) -> bool:
    return Path(path).suffix.lower() in IMAGE_SUFFIXES


def _read_other(path: Path) -> np.ndarray:
    """Pixels of any picture without tifffile: Pillow when installed, else Qt (both read LZW TIFFs and JPG)."""
    try:
        from PIL import Image  # type: ignore[import-not-found]

        with Image.open(str(path)) as im:
            if im.mode in ("P", "CMYK", "LA", "YCbCr", "HSV", "1"):
                im = im.convert("RGBA" if "A" in im.mode or im.mode == "P" else "RGB")
            return np.array(im)
    except ImportError:
        pass
    from PySide6.QtGui import QImage

    img = QImage(str(path))
    if img.isNull():
        raise ValueError(f"{Path(path).name} is not a picture that can be read")
    if img.isGrayscale():
        deep = img.depth() > 8
        img = img.convertToFormat(QImage.Format.Format_Grayscale16 if deep else QImage.Format.Format_Grayscale8)
        dtype, channels = (np.uint16, 1) if deep else (np.uint8, 1)
    else:
        img = img.convertToFormat(QImage.Format.Format_RGB888)
        dtype, channels = np.uint8, 3
    width, height, stride = img.width(), img.height(), img.bytesPerLine()
    raw = np.frombuffer(img.constBits(), dtype=np.uint8, count=stride * height).reshape(height, stride)
    item = np.dtype(dtype).itemsize
    raw = raw[:, : width * channels * item].copy()
    data = raw.view(dtype).reshape(height, width, channels) if channels > 1 else raw.view(dtype).reshape(height, width)
    return data


def _as_rgb_or_grey(array: np.ndarray) -> np.ndarray:
    """A colour picture whose three planes are equal is a grey channel."""
    if array.ndim == 3 and array.shape[-1] >= 3:
        rgb = array[..., :3]
        if np.array_equal(rgb[..., 0], rgb[..., 1]) and np.array_equal(rgb[..., 1], rgb[..., 2]):
            return np.ascontiguousarray(rgb[..., 0])
        return np.ascontiguousarray(rgb)
    return array


def _fallback_error(path: Path, exc: Exception) -> Exception:
    text = str(exc)
    if "imagecodecs" in text:
        text = ("the file is compressed (LZW or similar) and could not be opened. "
                "Install the decoder with: pip install imagecodecs")
    return ValueError(f"{Path(path).name}: {text}")


def _tifffile():
    try:
        import tifffile
    except ImportError as exc:  # pragma: no cover - depends on the install
        raise RuntimeError("Reading channel TIFFs needs tifffile: pip install tifffile") from exc
    return tifffile


def _pixel_size(tif) -> float | None:
    """Micrometres per pixel from OME, ImageJ or TIFF resolution tags; None if not trustworthy."""
    ome = getattr(tif, "ome_metadata", None)
    if ome:
        m = re.search(r'PhysicalSizeX="([0-9.eE+-]+)"', ome)
        if m:
            unit = re.search(r'PhysicalSizeXUnit="([^"]+)"', ome)
            value = float(m.group(1))
            u = unit.group(1) if unit else "µm"
            factor = {"µm": 1.0, "um": 1.0, "nm": 1e-3, "mm": 1e3}.get(u)
            if factor and value > 0:
                return value * factor
    page = tif.pages[0]
    tag = page.tags.get("XResolution")
    if tag is None:
        return None
    try:
        num, den = tag.value
        per_unit = num / den if den else 0.0
    except (TypeError, ValueError):
        return None
    if per_unit <= 0:
        return None
    ij = getattr(tif, "imagej_metadata", None) or {}
    unit = str(ij.get("unit", "")).lower() if isinstance(ij, dict) else ""
    if unit in ("micron", "µm", "um", "\\u00B5m"):
        return 1.0 / per_unit
    res_unit = page.tags.get("ResolutionUnit")
    res = int(getattr(res_unit, "value", 2)) if res_unit is not None else 2
    if res == 3:  # centimetre
        return 1e4 / per_unit
    return None  # inch resolutions (72/96/300 dpi) say nothing about the sample


def read_channel(path: Path) -> ChannelData:
    """Read one picture (grey or RGB). Multi-page TIFFs: the first page. Other formats and TIFFs tifffile
    cannot decode (LZW without imagecodecs) are read with Pillow or Qt."""
    pixel = None
    try:
        if not is_tiff(path):
            raise ValueError("not a TIFF")
        tifffile = _tifffile()
        with tifffile.TiffFile(str(path)) as tif:
            pixel = _pixel_size(tif)
            array = tif.pages[0].asarray()
    except Exception as exc:  # noqa: BLE001 - try the other readers before giving up
        try:
            array = _as_rgb_or_grey(_read_other(path))
        except Exception:  # noqa: BLE001
            raise _fallback_error(path, exc) from exc
    array = np.asarray(array)
    while array.ndim > 2 and array.shape[0] == 1:
        array = array[0]
    if array.ndim == 3 and array.shape[0] in (3, 4) and array.shape[-1] not in (3, 4):
        array = np.moveaxis(array, 0, -1)
    bits = int(array.dtype.itemsize * 8) if array.dtype.kind in "ui" else 32
    return ChannelData(array=array, bits=bits, pixel_size_um=pixel)


UNKNOWN_COLOURS = ("#33ff33", "#ff3333", "#4d6bff", "#ffd400", "#ff33ff", "#33e5ff", "#ffffff", "#ff9933")


@dataclass
class TiffInfo:
    kind: str  # "grey", "colour" (one channel saved in colour), "rgb" (a merged picture) or "channels"
    count: int  # planes: 1 for grey and colour, 3 for rgb, n for channels
    names: list[str]
    colours: list[str]
    width: int = 0
    height: int = 0


def colour_of(rgb: np.ndarray) -> str | None:
    """The colour of an RGB picture that holds one channel (Photoshop or ImageJ save a channel tinted:
    blue DAPI, green FITC, yellow, white ...). None when the planes differ: then it is a merge of channels."""
    if rgb.ndim != 3 or rgb.shape[-1] < 3:
        return None
    planes = rgb[..., :3].astype(np.float32)
    step = max(1, int(np.ceil(np.sqrt(planes.shape[0] * planes.shape[1] / 250_000))))
    planes = planes[::step, ::step]
    peaks = planes.reshape(-1, 3).max(axis=0)
    top = float(peaks.max())
    if top <= 0:
        return "#ffffff"
    intensity = planes.max(axis=-1)
    full = float(np.iinfo(rgb.dtype).max) if rgb.dtype.kind in "ui" else top
    for k in range(3):
        expected = intensity * (peaks[k] / top)
        if float(np.abs(planes[..., k] - expected).mean()) > 0.01 * full:
            return None
    r, g, b = (int(round(255 * p / top)) for p in peaks)
    return f"#{r:02x}{g:02x}{b:02x}"


def _rgb_info(array: np.ndarray) -> TiffInfo:
    h, w = array.shape[:2]
    if array.ndim == 3 and array.shape[-1] >= 3:
        colour = colour_of(array)
        if colour is not None:
            return TiffInfo("colour", 1, [""], [colour], w, h)
        return TiffInfo("rgb", 3, ["red", "green", "blue"], ["#ff3333", "#33ff33", "#4d6bff"], w, h)
    return TiffInfo("grey", 1, [""], ["#ffffff"], w, h)


def _channel_axis(series) -> tuple[str, int] | None:
    axes = str(getattr(series, "axes", ""))
    shape = tuple(series.shape)
    for ax in ("C", "I", "Q", "P"):
        if ax in axes and shape[axes.index(ax)] > 1:
            return ax, shape[axes.index(ax)]
    return None


def _imagej_luts(tif) -> list[str]:
    ij = getattr(tif, "imagej_metadata", None) or {}
    out = []
    for lut in ij.get("LUTs", []) if isinstance(ij, dict) else []:
        try:
            arr = np.asarray(lut)
            r, g, b = (int(arr[k][-1]) for k in range(3))
            out.append(f"#{r:02x}{g:02x}{b:02x}")
        except (IndexError, TypeError, ValueError):
            out.append("")
    return out


def _channel_names(tif, count: int) -> list[str]:
    ome = getattr(tif, "ome_metadata", None) or ""
    names = re.findall(r'<Channel[^>]*\sName="([^"]*)"', ome) if ome else []
    if not names:
        ij = getattr(tif, "imagej_metadata", None) or {}
        labels = ij.get("Labels") if isinstance(ij, dict) else None
        if isinstance(labels, (list, tuple)):
            names = [str(x) for x in labels]
    return [names[i] if i < len(names) else "" for i in range(count)]


def inspect_tiff(path: Path) -> TiffInfo:
    """What one picture holds: one grey plane, an RGB merge, or several channels (ImageJ, OME, multi-page)."""
    if not is_tiff(path):
        return _inspect_other(path)
    try:
        return _inspect_tiff(path)
    except Exception as exc:  # noqa: BLE001
        try:
            return _inspect_other(path)
        except Exception:  # noqa: BLE001
            raise _fallback_error(path, exc) from exc


def _inspect_other(path: Path) -> TiffInfo:
    return _rgb_info(_as_rgb_or_grey(_read_other(path)))


def _inspect_tiff(path: Path) -> TiffInfo:
    tifffile = _tifffile()
    with tifffile.TiffFile(str(path)) as tif:
        series = tif.series[0]
        found = _channel_axis(series)
        page = tif.pages[0]
        if found is not None:
            count = min(found[1], 16)
            colours = _imagej_luts(tif)
            colours = [colours[i] if i < len(colours) and colours[i] else UNKNOWN_COLOURS[i % len(UNKNOWN_COLOURS)]
                       for i in range(count)]
            return TiffInfo("channels", count, _channel_names(tif, count), colours,
                            int(getattr(page, "imagewidth", 0)), int(getattr(page, "imagelength", 0)))
        w, h = int(getattr(page, "imagewidth", 0)), int(getattr(page, "imagelength", 0))
        if int(getattr(page, "samplesperpixel", 1)) >= 3:
            array = _as_rgb_or_grey(np.asarray(page.asarray()))
            if array.ndim == 3 and array.shape[0] in (3, 4) and array.shape[-1] not in (3, 4):
                array = np.moveaxis(array, 0, -1)
            return _rgb_info(array)
        return TiffInfo("grey", 1, [""], ["#ffffff"], w, h)


def read_plane(path: Path, index: int = -1) -> ChannelData:
    """One plane of a TIFF. ``index`` -1 reads the first page as before; otherwise the channel (or RGB plane).

    Other axes (z, time) are reduced: z by maximum projection, the rest by taking the first.
    """
    if index < 0:
        return read_channel(path)
    try:
        return _read_plane_tiff(path, index)
    except Exception as exc:  # noqa: BLE001
        data = read_channel(path).array  # raises a clear error when nothing can read it
        if data.ndim == 3:
            data = data[..., min(index, data.shape[-1] - 1)]
        elif not is_tiff(path):
            pass
        else:
            raise _fallback_error(path, exc) from exc
        bits = int(data.dtype.itemsize * 8) if data.dtype.kind in "ui" else 32
        return ChannelData(array=np.ascontiguousarray(data), bits=bits, pixel_size_um=None)


def _read_plane_tiff(path: Path, index: int) -> ChannelData:
    if not is_tiff(path):
        raise ValueError("not a TIFF")
    tifffile = _tifffile()
    with tifffile.TiffFile(str(path)) as tif:
        series = tif.series[0]
        pixel = _pixel_size(tif)
        found = _channel_axis(series)
        if found is not None:
            axes = str(series.axes)
            data = np.asarray(series.asarray())
            data = np.take(data, index, axis=axes.index(found[0]))
            axes = axes.replace(found[0], "", 1)
            for ax in [a for a in axes if a not in "YXS"]:
                k = axes.index(ax)
                data = data.max(axis=k) if ax == "Z" else np.take(data, 0, axis=k)
                axes = axes.replace(ax, "", 1)
            if data.ndim == 3:  # leftover samples
                data = data[..., 0]
        else:
            data = np.asarray(tif.pages[0].asarray())
            if data.ndim == 3 and data.shape[0] in (3, 4) and data.shape[-1] not in (3, 4):
                data = np.moveaxis(data, 0, -1)
            if data.ndim == 3:
                data = data[..., min(index, data.shape[-1] - 1)]
    bits = int(data.dtype.itemsize * 8) if data.dtype.kind in "ui" else 32
    return ChannelData(array=data, bits=bits, pixel_size_um=pixel)


def split_rgb(array: np.ndarray) -> list[np.ndarray]:
    """An RGB(A) image as three grey channels (red, green, blue)."""
    if array.ndim != 3 or array.shape[-1] < 3:
        raise ValueError("Not an RGB image")
    return [np.ascontiguousarray(array[..., i]) for i in range(3)]


def to_grey(array: np.ndarray) -> np.ndarray:
    """A channel must be one plane; RGB files are reduced to their brightest plane."""
    if array.ndim == 2:
        return array
    if array.ndim == 3 and array.shape[-1] >= 3:
        return array[..., :3].max(axis=-1)
    raise ValueError(f"Unsupported image shape {array.shape}")


def _sample(array: np.ndarray, max_pixels: int = 2_000_000) -> np.ndarray:
    flat = array.reshape(-1)
    step = max(1, flat.size // max_pixels)
    return flat[::step]


def auto_levels(array: np.ndarray, saturated: float = DEFAULT_SATURATED) -> tuple[float, float]:
    """(low, high) leaving ``saturated`` percent of pixels clipped at each end."""
    data = _sample(to_grey(array)).astype(np.float64)
    if data.size == 0:
        return 0.0, 1.0
    low, high = np.percentile(data, [saturated, 100.0 - saturated])
    if high <= low:
        low, high = float(data.min()), float(data.max())
        if high <= low:
            high = low + 1.0
    return float(low), float(high)


def apply_levels(array: np.ndarray, low: float, high: float, gamma: float = 1.0) -> np.ndarray:
    """Linear stretch to 0..1 (then gamma). Returns float32."""
    grey = to_grey(array).astype(np.float32)
    span = float(high - low) or 1.0
    out = np.clip((grey - float(low)) / span, 0.0, 1.0)
    if gamma and abs(gamma - 1.0) > 1e-6:
        out = np.power(out, 1.0 / float(gamma), dtype=np.float32)
    return out


def hex_to_rgb(color: str) -> tuple[float, float, float]:
    c = str(color or "").strip().lstrip("#")
    if len(c) == 3:
        c = "".join(ch * 2 for ch in c)
    try:
        return (int(c[0:2], 16) / 255.0, int(c[2:4], 16) / 255.0, int(c[4:6], 16) / 255.0)
    except (ValueError, IndexError):
        return (1.0, 1.0, 1.0)


def downsample(array: np.ndarray, max_side: int) -> np.ndarray:
    """Every n-th pixel so the longer side is at most ``max_side`` (fast, for previews)."""
    if max_side <= 0:
        return array
    step = int(np.ceil(max(array.shape[0], array.shape[1]) / float(max_side)))
    return array[::step, ::step] if step > 1 else array


def composite(layers: list[tuple[np.ndarray, ChannelDisplay]], *, max_side: int = 0) -> np.ndarray:
    """Additive colour merge of the visible channels. Returns uint8 (H, W, 3).

    Channels whose display has no levels yet use auto levels.
    """
    visible = [(a, d) for a, d in layers if d.visible]
    if not layers:
        raise ValueError("No channels")
    shapes = {to_grey(a).shape for a, _ in layers}
    if len(shapes) != 1:
        raise ValueError("The channel images do not have the same size: " + ", ".join(f"{w}×{h}" for h, w in sorted(shapes)))
    h, w = next(iter(shapes))
    template = downsample(np.zeros((h, w), dtype=np.uint8), max_side)
    rgb = np.zeros(template.shape + (3,), dtype=np.float32)
    for array, display in visible:
        grey = downsample(to_grey(array), max_side)
        low, high = (display.low, display.high) if display.is_set else auto_levels(grey)
        scaled = apply_levels(grey, low, high, display.gamma)
        r, g, b = hex_to_rgb(display.color)
        rgb[..., 0] += scaled * r
        rgb[..., 1] += scaled * g
        rgb[..., 2] += scaled * b
    return (np.clip(rgb, 0.0, 1.0) * 255.0 + 0.5).astype(np.uint8)


def single_channel_view(array: np.ndarray, display: ChannelDisplay, *, grey: bool = False, max_side: int = 0) -> np.ndarray:
    """One channel alone, in its colour (or grey). uint8 (H, W, 3)."""
    d = ChannelDisplay(color="#ffffff" if grey else display.color, visible=True,
                       low=display.low, high=display.high, gamma=display.gamma)
    return composite([(array, d)], max_side=max_side)


def adjustment_log(channels: list[tuple[str, int, ChannelDisplay]]) -> list[str]:
    """Plain lines for the log and the figure legend. Input: (dye, bits, display)."""
    if not channels:
        return []
    bits = sorted({b for _, b, _ in channels})
    lines = [f"Raw TIFFs unchanged ({', '.join(str(b) for b in bits)} bit)"]
    nonlinear = False
    for dye, _bits, d in channels:
        if not d.visible:
            lines.append(f"{dye}: hidden")
            continue
        levels = f"levels {d.low:g} to {d.high:g}" if d.is_set else "auto levels"
        if abs(d.gamma - 1.0) > 1e-6:
            nonlinear = True
            lines.append(f"{dye}: {levels}, gamma {d.gamma:g} (non-linear)")
        else:
            lines.append(f"{dye}: {levels}, linear")
    if nonlinear:
        lines.append("Gamma other than 1.0 is a non-linear change: state it in the figure legend.")
    return lines


def legend_sentence(channels: list[tuple[str, int, ChannelDisplay]]) -> str:
    """One sentence for a figure legend about the adjustments."""
    gammas = [d for _, _, d in channels if d.visible and abs(d.gamma - 1.0) > 1e-6]
    if gammas:
        return "Brightness and contrast were adjusted for the whole image; gamma was changed for " + \
            ", ".join(dye for dye, _, d in channels if d in gammas) + "."
    return "Brightness and contrast were adjusted linearly for the whole image."


# ---------------------------------------------------------------------------
# PNG writer (previews are read by Qt; no Pillow needed)
# ---------------------------------------------------------------------------


def write_png(path: Path, rgb: np.ndarray) -> None:
    rgb = np.ascontiguousarray(rgb, dtype=np.uint8)
    if rgb.ndim != 3 or rgb.shape[2] != 3:
        raise ValueError("write_png expects (H, W, 3) uint8")
    h, w, _ = rgb.shape
    raw = b"".join(b"\x00" + rgb[y].tobytes() for y in range(h))

    def chunk(kind: bytes, data: bytes) -> bytes:
        return struct.pack(">I", len(data)) + kind + data + struct.pack(">I", zlib.crc32(kind + data) & 0xFFFFFFFF)

    png = b"\x89PNG\r\n\x1a\n" + chunk(b"IHDR", struct.pack(">IIBBBBB", w, h, 8, 2, 0, 0, 0)) + \
        chunk(b"IDAT", zlib.compress(raw, 6)) + chunk(b"IEND", b"")
    Path(path).write_bytes(png)
