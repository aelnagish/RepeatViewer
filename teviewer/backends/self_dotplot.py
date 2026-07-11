"""Fast in-process self dotplot for single-sequence FASTA preview."""

from __future__ import annotations

import gzip
from pathlib import Path

import numpy as np
from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QImage, QPainter, QPen

_CANVAS_PX = 1024
_PLOT_MARGIN = 48
_DNA = set("ACGTUN")
_MAX_KMER_POSITIONS = 800  # cap before O(k^2) pair loops get too slow


def read_single_fasta_sequence(path: Path) -> tuple[str, str]:
    resolved = path.expanduser().resolve()
    header = ""
    parts: list[str] = []
    if resolved.name.lower().endswith(".gz"):
        opener = gzip.open(resolved, "rt", encoding="utf-8", errors="replace")
    else:
        opener = resolved.open(encoding="utf-8", errors="replace")
    with opener as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if parts:
                    break
                header = line[1:].strip()
                continue
            if not header:
                header = resolved.stem
            parts.append("".join(ch for ch in line.upper() if ch in _DNA))
    sequence = "".join(parts)
    if not sequence:
        raise ValueError(f"No sequence in {resolved.name}")
    return header, sequence


def _build_kmer_index(sequence: str, *, word_size: int) -> dict[str, list[int]]:
    seq = "".join(ch for ch in sequence.upper() if ch in _DNA)
    index: dict[str, list[int]] = {}
    for pos in range(len(seq) - word_size + 1):
        kmer = seq[pos : pos + word_size]
        if "N" in kmer:
            continue
        index.setdefault(kmer, []).append(pos)
    return index


def _positions_for_plotting(positions: list[int]) -> list[int]:
    """Thin very large k-mer hit lists so pair loops stay responsive."""
    sorted_pos = sorted(positions)
    if len(sorted_pos) <= _MAX_KMER_POSITIONS:
        return sorted_pos
    step = len(sorted_pos) // _MAX_KMER_POSITIONS + 1
    return sorted_pos[::step]


def _rasterize_dotplot(index: dict[str, list[int]], *, seq_len: int, plot_px: int) -> np.ndarray:
    """
    Standard dottup-style dotplot: one pixel per matching word pair (i, j).

    Lines in the image emerge from dense dots (tandem/direct repeats), not from
    drawing segments between hits.
    """
    grid = np.zeros((plot_px, plot_px), dtype=np.uint8)
    if seq_len <= 1:
        return grid

    scale = (plot_px - 1) / (seq_len - 1)
    xs: list[int] = []
    ys: list[int] = []

    for positions in index.values():
        plot_positions = _positions_for_plotting(positions)
        for a in plot_positions:
            for b in plot_positions:
                xs.append(int(round(a * scale)))
                ys.append(plot_px - 1 - int(round(b * scale)))

    if not xs:
        return grid

    xs_arr = np.asarray(xs, dtype=np.int32)
    ys_arr = np.asarray(ys, dtype=np.int32)
    valid = (xs_arr >= 0) & (xs_arr < plot_px) & (ys_arr >= 0) & (ys_arr < plot_px)
    grid[ys_arr[valid], xs_arr[valid]] = 1
    return grid


def _grid_to_qimage(grid: np.ndarray) -> QImage:
    height, width = grid.shape
    rgba = np.zeros((height, width, 4), dtype=np.uint8)
    rgba[..., 3] = 255
    rgba[grid > 0] = [26, 26, 26, 255]
    rgba[grid == 0] = [255, 255, 255, 255]
    image = QImage(rgba.data, width, height, width * 4, QImage.Format.Format_RGBA8888)
    return image.copy()


def generate_self_dotplot_png(
    fasta_path: Path,
    output_png: Path,
    *,
    word_size: int,
) -> Path:
    """Render a self dotplot PNG (dottup / Flexidot-style match points)."""
    _header, sequence = read_single_fasta_sequence(fasta_path)
    seq = "".join(ch for ch in sequence.upper() if ch in _DNA)
    seq_len = max(len(seq), 1)
    plot_px = _CANVAS_PX - _PLOT_MARGIN * 2

    image = QImage(_CANVAS_PX, _CANVAS_PX, QImage.Format.Format_ARGB32)
    image.fill(QColor("#ffffff"))

    index = _build_kmer_index(sequence, word_size=word_size)
    painter = QPainter(image)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)
    painter.setPen(QPen(QColor("#d8d8dc"), 1))
    painter.drawRect(_PLOT_MARGIN, _PLOT_MARGIN, plot_px, plot_px)

    if index:
        grid = _rasterize_dotplot(index, seq_len=seq_len, plot_px=plot_px)
        plot_image = _grid_to_qimage(grid)
        painter.drawImage(_PLOT_MARGIN, _PLOT_MARGIN, plot_image)
    else:
        painter.setPen(QColor("#666666"))
        painter.drawText(
            image.rect(),
            int(Qt.AlignmentFlag.AlignCenter),
            "No matches for this word size.",
        )

    painter.end()

    output_png = output_png.expanduser().resolve()
    output_png.parent.mkdir(parents=True, exist_ok=True)
    if not image.save(str(output_png), "PNG"):
        raise OSError(f"Could not write dotplot image: {output_png}")
    return output_png
