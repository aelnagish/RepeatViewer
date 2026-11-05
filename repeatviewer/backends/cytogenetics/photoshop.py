"""Photoshop layers saved inside a TIFF (tag 37724, "Save layers" in Photoshop's TIFF options).

Qt-free. Reads 8-bit, 16-bit and 32-bit layers in either byte order, raw, RLE or ZIP compressed.
"""

from __future__ import annotations

import struct
import zlib
from dataclasses import dataclass

import numpy as np

TAG = 37724


@dataclass
class PsLayer:
    name: str
    visible: bool
    text: bool  # a type layer (chromosome numbers, labels)
    rgb: np.ndarray  # H x W x 3 on the full canvas
    coverage: float = 1.0  # share of the canvas the layer's box covers (cut-out chromosomes are small)
    box: tuple[int, int, int, int] = (0, 0, 0, 0)  # top, left, bottom, right on the canvas

    @property
    def empty(self) -> bool:
        return not bool(self.rgb.any())


def picture_layers(layers: list[PsLayer]) -> list[PsLayer]:
    """The layers that are channels: visible pixel layers with something in them, no text.

    The channels are the largest layers; layers less than half their size are cut-out chromosomes, arrows or
    labels. When the channels fill only part of the canvas (a karyotype figure), they are cropped to that part.
    """
    shown = [layer for layer in layers if layer.visible and not layer.text and not layer.empty]
    if not shown:
        return []
    biggest = max(layer.coverage for layer in shown)
    kept = [layer for layer in shown if layer.coverage >= MIN_SHARE_OF_BIGGEST * biggest]
    top = min(layer.box[0] for layer in kept)
    left = min(layer.box[1] for layer in kept)
    bottom = max(layer.box[2] for layer in kept)
    right = max(layer.box[3] for layer in kept)
    h, w = kept[0].rgb.shape[:2]
    if bottom > top and right > left and (bottom - top) * (right - left) < 0.9 * h * w:
        kept = [PsLayer(layer.name, layer.visible, layer.text, layer.rgb[top:bottom, left:right].copy(),
                        layer.coverage, layer.box) for layer in kept]
    return kept


MIN_SHARE_OF_BIGGEST = 0.5


def _packbits(data, size):
    out = bytearray(); i = 0
    while i < len(data) and len(out) < size:
        n = data[i]; i += 1
        if n < 128:
            out += data[i:i + n + 1]; i += n + 1
        elif n > 128:
            out += bytes([data[i]]) * (257 - n); i += 1
    return bytes(out[:size]).ljust(size, b"\0")

def layers_from_block(block: bytes, width: int, height: int) -> list["PsLayer"]:
    """Every pixel layer in a Photoshop data block (TIFF tag 37724), placed on the full canvas."""
    head = b"Adobe Photoshop Document Data Block\0"
    if not block.startswith(head):
        return []
    pos = len(head)
    if block[pos:pos + 4] in (b"8BIM", b"8B64"):
        e = ">"
    elif block[pos:pos + 4] in (b"MIB8", b"46B8"):
        e = "<"
    else:
        return []
    def u(fmt, at): return struct.unpack_from(e + fmt, block, at)[0]
    def key(at):
        k = block[at:at + 4]
        return k if e == ">" else k[::-1]
    layers = []
    while pos + 12 <= len(block):
        if key(pos) not in (b"8BIM", b"8B64"):
            break
        k = key(pos + 4); length = u("I", pos + 8); body = pos + 12
        if k in (b"Layr", b"Lr16", b"Lr32"):
            depth = {b"Layr": 1, b"Lr16": 2, b"Lr32": 4}[k]
            layers = _parse_layers(block, body, e, u, depth, width, height)
            if layers:
                return layers
        pos = body + length + ((-length) % 4)
    return layers

def _parse_layers(block, at, e, u, depth, W, H):
    count = abs(u("h", at)); at += 2
    recs = []
    for _ in range(count):
        top, left, bottom, right = struct.unpack_from(e + "4i", block, at); at += 16
        nch = u("H", at); at += 2
        chans = []
        for _c in range(nch):
            chans.append((u("h", at), u("I", at + 2))); at += 6
        at += 4  # blend signature
        at += 4  # blend key
        opacity, clipping, flags = block[at], block[at + 1], block[at + 2]; at += 4
        extra = u("I", at); at += 4
        end = at + extra
        mask_len = u("I", at); at += 4 + mask_len
        blend_len = u("I", at); at += 4 + blend_len
        nlen = block[at]
        name = block[at + 1:at + 1 + nlen].decode("latin-1", "replace")
        info = block[at + 1 + nlen:end]
        text = b"TySh" in info or b"hSyT" in info
        at = end
        recs.append(dict(box=(top, left, bottom, right), chans=chans, name=name, visible=not (flags & 2),
                         opacity=opacity, text=text))
    dtype = {1: np.uint8, 2: np.dtype(e + "u2"), 4: np.dtype(e + "f4")}[depth]
    out = []
    for r in recs:
        top, left, bottom, right = r["box"]; h, w = bottom - top, right - left
        planes = {}
        for cid, clen in r["chans"]:
            comp = u("H", at); data = block[at + 2:at + clen]; at += clen
            if cid < 0 or h <= 0 or w <= 0:
                continue  # transparency and masks
            size = h * w * depth
            if comp == 0:
                raw = data[:size]
            elif comp == 1:
                counts = struct.unpack_from(e + f"{h}H", data, 0)
                p = 2 * h; rows = []
                for c in counts:
                    rows.append(_packbits(data[p:p + c], w * depth)); p += c
                raw = b"".join(rows)
            else:
                raw = zlib.decompress(data)
                if comp == 3 and depth == 1:
                    a = np.frombuffer(raw, np.uint8).reshape(h, w)
                    raw = np.cumsum(a, axis=1, dtype=np.uint8).tobytes()
            planes[cid] = np.frombuffer(raw[:size], dtype).reshape(h, w)
        if not planes or h <= 0 or w <= 0:
            continue
        img = np.zeros((H, W, 3), dtype=planes[next(iter(planes))].dtype)
        y0, x0 = max(top, 0), max(left, 0); y1, x1 = min(bottom, H), min(right, W)
        for cid in range(3):
            p = planes.get(cid, planes.get(0))
            img[y0:y1, x0:x1, cid] = p[y0 - top:y1 - top, x0 - left:x1 - left]
        cover = max(0, y1 - y0) * max(0, x1 - x0) / float(W * H) if W and H else 0.0
        out.append(PsLayer(name=r["name"], visible=r["visible"], text=r["text"], rgb=img, coverage=cover,
                           box=(y0, x0, max(y0, y1), max(x0, x1))))
    return out
