"""MSA overview plot (TEtrimmer / CIAlign-style residue heatmap).

Adapted from TEtrimmer ``cialign.drawMiniAlignment`` / ``MSAcluster.plot_msa``
without crop arrows or start/end annotations.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Literal

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt
import numpy as np
from Bio import SeqIO
from matplotlib.colors import ListedColormap

MSA_PLOT_SEQ_AUTO = "auto"
MSA_PLOT_SEQ_NUCLEOTIDE = "nucleotide"
MSA_PLOT_SEQ_PROTEIN = "protein"
MSA_PLOT_FORMAT_PNG = "png"
MSA_PLOT_FORMAT_PDF = "pdf"

SequenceType = Literal["auto", "nucleotide", "protein"]
OutputFormat = Literal["png", "pdf"]

# IUPAC nucleotide ambiguity codes (default white so users can recolour them).
_NUCLEOTIDE_AMBIGUITY_LETTERS = ("I", "R", "Y", "W", "S", "M", "K", "B", "H", "D", "V")

# Nucleotide defaults (A green, C purple-blue, G black, T/U pink); ambiguities white.
NUCLEOTIDE_COLORS_DEFAULT: dict[str, str] = {
    "A": "#00C800",
    "C": "#7040E8",
    "G": "#000000",
    "T": "#F06090",
    "U": "#F06090",
    "N": "#C8C8C8",
    "X": "#C8C8C8",
    "-": "#FFFFFF",
    **{letter: "#FFFFFF" for letter in _NUCLEOTIDE_AMBIGUITY_LETTERS},
}

# Same Wong / Paul Tol palette as Preferences → Color-blind nucleotide.
NUCLEOTIDE_COLORS_COLORBLIND: dict[str, str] = {
    "A": "#E69F00",
    "C": "#0072B2",
    "G": "#000000",
    "T": "#009E73",
    "U": "#009E73",
    "N": "#C8C8C8",
    "X": "#C8C8C8",
    "-": "#FFFFFF",
    **{letter: "#FFFFFF" for letter in _NUCLEOTIDE_AMBIGUITY_LETTERS},
}

# Categories for the nucleotide colour picker.
NUCLEOTIDE_CATEGORIES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("Bases", ("A", "T", "U", "G", "C")),
    ("Ambiguity", _NUCLEOTIDE_AMBIGUITY_LETTERS),
    ("Other / gap", ("N", "X", "-")),
)

# Jalview / ClustalX residue colours.
AMINO_ACID_COLORS_CLUSTAL: dict[str, str] = {
    "A": "#80A0F0",
    "R": "#F01505",
    "N": "#00FF00",
    "D": "#C048C0",
    "C": "#F08080",
    "Q": "#00FF00",
    "E": "#C048C0",
    "G": "#F09048",
    "H": "#15A4A4",
    "I": "#80A0F0",
    "L": "#80A0F0",
    "K": "#F01505",
    "M": "#80A0F0",
    "F": "#80A0F0",
    "P": "#FFFF00",
    "S": "#00FF00",
    "T": "#00FF00",
    "W": "#80A0F0",
    "Y": "#15A4A4",
    "V": "#80A0F0",
    "B": "#C8C8C8",
    "Z": "#C8C8C8",
    "X": "#C8C8C8",
    "*": "#C8C8C8",
    "-": "#FFFFFF",
}

# Jalview Zappo (physicochemical property groups).
AMINO_ACID_COLORS_ZAPPO: dict[str, str] = {
    "A": "#FF7F7F",
    "I": "#FF7F7F",
    "L": "#FF7F7F",
    "M": "#FF7F7F",
    "V": "#FF7F7F",
    "F": "#FFC44D",
    "W": "#FFC44D",
    "Y": "#FFC44D",
    "K": "#7C7CFF",
    "R": "#7C7CFF",
    "H": "#FF00FF",
    "D": "#FF0000",
    "E": "#FF0000",
    "N": "#00FF00",
    "Q": "#00FF00",
    "S": "#00FF00",
    "T": "#00FF00",
    "C": "#FFFF00",
    "G": "#FF7F00",
    "P": "#FF7F00",
    "B": "#FFFFFF",
    "Z": "#FFFFFF",
    "X": "#FFFFFF",
    "*": "#FFFFFF",
    "-": "#FFFFFF",
}

# Jalview Taylor (high-contrast per-residue colours).
AMINO_ACID_COLORS_TAYLOR: dict[str, str] = {
    "A": "#CCFF00",
    "R": "#0000FF",
    "N": "#CC00FF",
    "D": "#FF0000",
    "C": "#FFFF00",
    "Q": "#FF00CC",
    "E": "#FF0066",
    "G": "#FF9900",
    "H": "#0066FF",
    "I": "#66FF00",
    "L": "#33FF00",
    "K": "#6600FF",
    "M": "#00FF00",
    "F": "#00FF66",
    "P": "#FFCC00",
    "S": "#FF3300",
    "T": "#FF6600",
    "W": "#00CCFF",
    "Y": "#00FFCC",
    "V": "#99FF00",
    "B": "#FFFFFF",
    "Z": "#FFFFFF",
    "X": "#FFFFFF",
    "*": "#FFFFFF",
    "-": "#FFFFFF",
}

# Clustal X–style amino-acid categories for the colour picker layout.
# (C is listed under Cysteine rather than Hydrophobic to avoid duplicates.)
AMINO_ACID_CATEGORIES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("Hydrophobic", ("A", "I", "L", "M", "F", "W", "V")),
    ("Positive charge", ("K", "R")),
    ("Negative charge", ("D", "E")),
    ("Polar", ("N", "Q", "S", "T")),
    ("Cysteine", ("C",)),
    ("Glycine", ("G",)),
    ("Proline", ("P",)),
    ("Aromatic", ("H", "Y")),
    ("Other / gap", ("B", "Z", "X", "*", "-")),
)

DEFAULT_NUCLEOTIDE_COLORS = NUCLEOTIDE_COLORS_DEFAULT
DEFAULT_AMINO_ACID_COLORS = AMINO_ACID_COLORS_CLUSTAL

NUCLEOTIDE_COLOR_PRESETS: tuple[tuple[str, dict[str, str]], ...] = (
    ("Default", NUCLEOTIDE_COLORS_DEFAULT),
    ("Color-blind", NUCLEOTIDE_COLORS_COLORBLIND),
)

AMINO_ACID_COLOR_PRESETS: tuple[tuple[str, dict[str, str]], ...] = (
    ("Clustal", AMINO_ACID_COLORS_CLUSTAL),
    ("Zappo", AMINO_ACID_COLORS_ZAPPO),
    ("Taylor", AMINO_ACID_COLORS_TAYLOR),
)

_NUCLEOTIDE_CHARS = frozenset("ACGTUNacgtun-")
_AMINO_ACID_CHARS = frozenset("ACDEFGHIKLMNPQRSTVWYBXZ*acdefghiklmnpqrstvwybxz*-")


@dataclass
class MsaPlotOptions:
    """User settings for the MSA overview plot."""

    sequence_type: str = MSA_PLOT_SEQ_AUTO
    output_format: str = MSA_PLOT_FORMAT_PNG
    dpi: int = 100
    # Figure centimeters per sequence row (user-adjustable; converted to inches for matplotlib).
    sequence_height_cm: float = 0.5
    # White gap between consecutive sequence rows (cm); 0 disables gaps.
    sequence_distance_cm: float = 0.01
    nucleotide_colors: dict[str, str] = field(
        default_factory=lambda: dict(DEFAULT_NUCLEOTIDE_COLORS)
    )
    amino_acid_colors: dict[str, str] = field(
        default_factory=lambda: dict(DEFAULT_AMINO_ACID_COLORS)
    )


def default_msa_plot_options() -> MsaPlotOptions:
    return MsaPlotOptions()


def msa_plot_output_path(input_fasta: Path, *, output_format: str) -> Path:
    """Default output next to the input MSA (keeps the original filename)."""
    from repeatviewer.tool_output_paths import kept_input_name

    path = input_fasta.expanduser().resolve()
    ext = "pdf" if str(output_format).lower() == MSA_PLOT_FORMAT_PDF else "png"
    return path.parent / f"{kept_input_name(path)}_msa_plot.{ext}"


def fasta_to_alignment_array(path: Path) -> tuple[np.ndarray, list[str]]:
    """Load a FASTA alignment as an uppercase character array (rows = sequences)."""
    names: list[str] = []
    seqs: list[list[str]] = []
    for record in SeqIO.parse(str(path.expanduser().resolve()), "fasta"):
        names.append(str(record.id))
        seqs.append(list(str(record.seq).upper()))
    if not seqs:
        raise ValueError(f"No sequences found in {path}")
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(
            "Sequences have different lengths; provide an aligned FASTA (MSA)."
        )
    if lengths.pop() == 0:
        raise ValueError("Alignment has zero columns.")
    return np.asarray(seqs, dtype="<U1"), names


def detect_sequence_type(arr: np.ndarray) -> str:
    """Guess nucleotide vs protein from residue composition (gaps ignored)."""
    flat = arr.ravel()
    residues = [ch for ch in flat if ch not in {"-", ".", " "}]
    if not residues:
        return MSA_PLOT_SEQ_NUCLEOTIDE
    sample = residues[: min(len(residues), 50_000)]
    n_nuc = sum(1 for ch in sample if ch in _NUCLEOTIDE_CHARS)
    n_aa = sum(1 for ch in sample if ch in _AMINO_ACID_CHARS)
    # Prefer nucleotide when almost everything fits ACGTUN.
    if n_nuc / max(1, len(sample)) >= 0.95:
        return MSA_PLOT_SEQ_NUCLEOTIDE
    if n_aa / max(1, len(sample)) >= 0.8:
        return MSA_PLOT_SEQ_PROTEIN
    return MSA_PLOT_SEQ_NUCLEOTIDE if n_nuc >= n_aa else MSA_PLOT_SEQ_PROTEIN


def _normalize_color_map(
    colors: dict[str, str] | None,
    defaults: dict[str, str],
) -> dict[str, str]:
    merged = {key.upper(): value for key, value in defaults.items()}
    if colors:
        for key, value in colors.items():
            text = str(value).strip()
            if not text:
                continue
            if not text.startswith("#"):
                text = f"#{text}"
            merged[str(key).upper()] = text
    if "-" not in merged:
        merged["-"] = "#FFFFFF"
    return merged


def _colors_for_type(sequence_type: str, options: MsaPlotOptions) -> dict[str, str]:
    if sequence_type == MSA_PLOT_SEQ_PROTEIN:
        return _normalize_color_map(options.amino_acid_colors, DEFAULT_AMINO_ACID_COLORS)
    return _normalize_color_map(options.nucleotide_colors, DEFAULT_NUCLEOTIDE_COLORS)


def _numeric_alignment(
    arr: np.ndarray,
    color_map: dict[str, str],
) -> tuple[np.ndarray, ListedColormap, int]:
    """Map residues to integer codes and a matching colormap (unknown → gap color).

    Returns ``(numeric, cmap, gap_index)`` where ``gap_index`` is the code used for
    blank spacer rows between sequences (the ``-`` / white color).
    """
    gap_color = color_map.get("-", "#FFFFFF")
    # Stable order: known residues first, then anything else seen in the MSA.
    ordered: list[str] = []
    for key in color_map:
        if key not in ordered:
            ordered.append(key)
    for ch in np.unique(arr):
        token = str(ch)
        if token not in ordered:
            ordered.append(token)
            color_map[token] = gap_color if token in {"-", ".", " "} else color_map.get(
                "N", color_map.get("X", gap_color)
            )

    index = {ch: i for i, ch in enumerate(ordered)}
    numeric = np.empty(arr.shape, dtype=np.float64)
    for y in range(arr.shape[0]):
        for x in range(arr.shape[1]):
            numeric[y, x] = index[str(arr[y, x])]
    cmap = ListedColormap([color_map[ch] for ch in ordered])
    gap_index = int(index.get("-", 0))
    return numeric, cmap, gap_index


# Vertical resolution used when expanding each sequence row for white gaps.
_SEQ_ROW_PIXELS = 20


def _gap_pixels(
    sequence_height_cm: float,
    sequence_distance_cm: float,
    *,
    row_pixels: int = _SEQ_ROW_PIXELS,
) -> int:
    """Map distance cm → integer spacer rows relative to sequence row height."""
    height = max(1e-6, float(sequence_height_cm))
    gap = max(0.0, float(sequence_distance_cm))
    if gap <= 0:
        return 0
    return max(1, int(round(row_pixels * gap / height)))


def _insert_row_gaps(
    numeric: np.ndarray,
    gap_code: float,
    *,
    row_pixels: int = _SEQ_ROW_PIXELS,
    gap_pixels: int = 1,
) -> tuple[np.ndarray, list[float]]:
    """Expand rows so a thin white band separates each sequence; return y-centers."""
    n_seq, n_pos = numeric.shape
    if n_seq == 0:
        return numeric, []
    blocks: list[np.ndarray] = []
    centers: list[float] = []
    cursor = 0.0
    for i in range(n_seq):
        blocks.append(np.repeat(numeric[i : i + 1], row_pixels, axis=0))
        centers.append(cursor + (row_pixels - 1) / 2.0)
        cursor += row_pixels
        if i < n_seq - 1 and gap_pixels > 0:
            blocks.append(
                np.full((gap_pixels, n_pos), gap_code, dtype=numeric.dtype)
            )
            cursor += gap_pixels
    return np.vstack(blocks), centers


def _figure_size(
    n_seq: int,
    n_pos: int,
    *,
    sequence_height_cm: float,
    sequence_distance_cm: float = 0.01,
) -> tuple[float, float]:
    width = min(48.0, max(8.0, n_pos / 80.0))
    row_h = max(0.05, min(2.0, float(sequence_height_cm) / 2.54))
    gap_h = max(0.0, float(sequence_distance_cm) / 2.54)
    height = min(
        80.0,
        max(2.0, n_seq * row_h + max(0, n_seq - 1) * gap_h),
    )
    return width, height


def _x_tick_step(n_pos: int, *, target_ticks: int = 12) -> int:
    """Choose a round step so the x-axis has roughly ``target_ticks`` labels."""
    if n_pos <= 1:
        return 1
    if n_pos <= target_ticks:
        return 1
    raw = n_pos / float(target_ticks)
    magnitude = 10 ** int(np.floor(np.log10(raw)))
    for mult in (1, 2, 5, 10):
        step = int(mult * magnitude)
        if step <= 0:
            continue
        if n_pos / step <= target_ticks * 1.35:
            return step
    return max(1, int(10 * magnitude))


def write_msa_plot(
    input_fasta: Path,
    output_path: Path | None = None,
    *,
    options: MsaPlotOptions | None = None,
) -> Path:
    """Render an MSA overview image and return the written path."""
    opts = options or default_msa_plot_options()
    input_fasta = input_fasta.expanduser().resolve()
    if not input_fasta.is_file():
        raise FileNotFoundError(f"MSA FASTA not found: {input_fasta}")

    fmt = str(opts.output_format).lower().strip()
    if fmt not in {MSA_PLOT_FORMAT_PNG, MSA_PLOT_FORMAT_PDF}:
        raise ValueError(f"Unsupported output format: {opts.output_format}")
    out = (
        output_path.expanduser().resolve()
        if output_path is not None
        else msa_plot_output_path(input_fasta, output_format=fmt)
    )
    if out.suffix.lower() not in {".png", ".pdf"}:
        out = out.with_suffix(f".{fmt}")

    arr, names = fasta_to_alignment_array(input_fasta)
    seq_type = str(opts.sequence_type).lower().strip()
    if seq_type == MSA_PLOT_SEQ_AUTO:
        seq_type = detect_sequence_type(arr)
    elif seq_type not in {MSA_PLOT_SEQ_NUCLEOTIDE, MSA_PLOT_SEQ_PROTEIN}:
        raise ValueError(f"Unknown sequence type: {opts.sequence_type}")

    # Normalize ambiguity / gap synonyms before coloring.
    arr = np.where(np.isin(arr, [".", " "]), "-", arr)

    color_map = _colors_for_type(seq_type, opts)
    # Keep T and U identical for DNA/RNA plots.
    if seq_type == MSA_PLOT_SEQ_NUCLEOTIDE and "U" in color_map:
        color_map["T"] = color_map["U"]
    numeric, cmap, gap_index = _numeric_alignment(arr, color_map)
    n_seq, n_pos = numeric.shape
    dpi = max(50, min(600, int(opts.dpi)))

    display, y_centers = _insert_row_gaps(
        numeric,
        float(gap_index),
        gap_pixels=_gap_pixels(opts.sequence_height_cm, opts.sequence_distance_cm),
    )
    n_display = display.shape[0]

    fig, ax = plt.subplots(
        figsize=_figure_size(
            n_seq,
            n_pos,
            sequence_height_cm=opts.sequence_height_cm,
            sequence_distance_cm=opts.sequence_distance_cm,
        ),
        dpi=dpi,
        facecolor="white",
    )
    # Pin vmin/vmax so integer codes map 1:1 onto the ListedColormap
    # (otherwise matplotlib stretches the used codes across the whole palette).
    ax.imshow(
        display,
        cmap=cmap,
        aspect="auto",
        interpolation="nearest",
        vmin=-0.5,
        vmax=len(cmap.colors) - 0.5,
    )
    ax.set_xlim(-0.5, n_pos - 0.5)
    ax.set_ylim(n_display - 0.5, -0.5)  # first sequence at the top
    ax.set_facecolor("white")

    x_step = _x_tick_step(n_pos)
    x_ticks = list(range(0, n_pos, x_step))
    if n_pos > 1 and x_ticks[-1] != n_pos - 1 and (n_pos - 1 - x_ticks[-1]) >= max(1, x_step // 2):
        x_ticks.append(n_pos - 1)
    ax.set_xticks(x_ticks)
    ax.set_xticklabels([str(t) for t in x_ticks])

    if n_seq <= 11:
        y_step = 1
    elif n_seq <= 500:
        y_step = 10
    else:
        y_step = 100
    y_indices = list(range(0, n_seq, y_step))
    ax.set_yticks([y_centers[i] for i in y_indices])
    ax.set_yticklabels([str(i) for i in y_indices])
    ax.tick_params(axis="both", labelsize=10)
    ax.set_xlabel("Alignment position")
    ax.set_ylabel("Sequence")

    kind = "Protein" if seq_type == MSA_PLOT_SEQ_PROTEIN else "Nucleotide"
    ax.set_title(
        f"{input_fasta.name}  ({kind}, {n_seq} × {n_pos})",
        fontsize=12,
        pad=10,
    )

    fig.tight_layout()
    out.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out, dpi=dpi, bbox_inches="tight", facecolor="white")
    plt.close(fig)
    del numeric, arr, names
    if not out.is_file():
        raise RuntimeError(f"Failed to write MSA plot: {out}")
    return out
