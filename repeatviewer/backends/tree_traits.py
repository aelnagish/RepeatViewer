"""Trait table for the Tree tab: colour tips by chromosome, array, methylation…

One tab-separated table per tree, next to it (rule 1):
``sat.fa_tree.treefile`` → ``sat.fa_tree.treefile_traits.tsv``. The first
column is the tip name; every other column is a trait. It is the same idea as
iTOL datasets or ggtree's ``%<+%``, so the figure is reproducible and the table
can go into a supplement. Users can edit it or add columns in any spreadsheet.

:func:`build_trait_table` pre-fills what RepeatViewer already knows:

* **header / chromosome / start / end / strand**: from each tip's original
  FASTA header (``chrom:start-end(strand)``, as written when copies are
  extracted), found through the tree run's ``names.tsv``;
* **array**: copies on one chromosome closer than ``array_gap`` bp form one
  tandem array, labelled by its span (e.g. ``Chr1:12.40-12.91 Mb``);
* **methylation**: % methylated over each copy's interval, from a Bismark
  ``.cov`` or a modkit bedMethyl file (read once, streamed);
* **source_file / fish**: when a tree combines several FASTA files, which file
  each tip came from and whether FISH images are attached to that file.

Columns are numeric (colour scale) or categorical (distinct colours).
Pure Python + numpy, no Qt; unit-tested in ``tests/test_tree_traits.py``.
"""

from __future__ import annotations

import csv
import gzip
import math
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

DEFAULT_ARRAY_GAP = 10_000
NAME_COLUMN = "name"
MISSING = {"", "na", "n/a", "nan", "none", "-"}

# Okabe & Ito (2008) followed by Paul Tol's "muted" set: 17 colour-blind-safe hues.
CATEGORICAL = [
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7", "#000000",
    "#332288", "#88CCEE", "#44AA99", "#117733", "#999933", "#DDCC77", "#CC6677",
    "#882255", "#AA4499",
]
OTHER_COLOR = "#BBBBBB"
# viridis at 0, 0.25, 0.5, 0.75, 1
_VIRIDIS = [(68, 1, 84), (59, 82, 139), (33, 145, 140), (94, 201, 98), (253, 231, 37)]


def traits_path(tree_path: Path) -> Path:
    tree_path = Path(tree_path)
    return tree_path.parent / f"{tree_path.name}_traits.tsv"


def names_table_path(tree_path: Path) -> Path | None:
    """``names.tsv`` of the tree run that wrote ``tree_path``, if there is one."""
    tree_path = Path(tree_path)
    name = tree_path.name
    for suffix in ("_masked.treefile", ".treefile"):
        if name.endswith(suffix):
            base = name[: -len(suffix)]
            candidate = tree_path.parent / f"{base}_details" / f"{base}.names.tsv"
            return candidate if candidate.is_file() else None
    return None


# --------------------------------------------------------------------------
# Reading and writing the table
# --------------------------------------------------------------------------


@dataclass
class TraitTable:
    columns: list[str]                                   # trait columns (not "name")
    rows: dict[str, dict[str, str]] = field(default_factory=dict)  # tip → column → value

    def values(self, column: str) -> dict[str, str]:
        return {tip: row.get(column, "") for tip, row in self.rows.items()}


def read_trait_table(path: Path) -> TraitTable:
    with Path(path).open(encoding="utf-8-sig", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if not header:
            raise ValueError("The trait table is empty.")
        header = [h.strip() for h in header]
        rows: dict[str, dict[str, str]] = {}
        for line in reader:
            if not line or not line[0].strip():
                continue
            values = [v.strip() for v in line] + [""] * (len(header) - len(line))
            rows[values[0]] = dict(zip(header[1:], values[1:]))
    return TraitTable([h for h in header[1:] if h], rows)


def write_trait_table(path: Path, table: TraitTable) -> Path:
    with Path(path).open("w", encoding="utf-8", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t", lineterminator="\n")
        writer.writerow([NAME_COLUMN, *table.columns])
        for tip, row in table.rows.items():
            writer.writerow([tip, *(row.get(c, "") for c in table.columns)])
    return Path(path)


# --------------------------------------------------------------------------
# Colours
# --------------------------------------------------------------------------


def _number(value: str) -> float | None:
    if value is None or str(value).strip().lower() in MISSING:
        return None
    try:
        x = float(str(value).strip().rstrip("%"))
    except ValueError:
        return None
    return x if math.isfinite(x) else None


def is_numeric(values: list[str]) -> bool:
    present = [v for v in values if str(v).strip().lower() not in MISSING]
    return bool(present) and all(_number(v) is not None for v in present)


def viridis(t: float) -> str:
    t = min(1.0, max(0.0, t)) * (len(_VIRIDIS) - 1)
    i = min(int(t), len(_VIRIDIS) - 2)
    f = t - i
    a, b = _VIRIDIS[i], _VIRIDIS[i + 1]
    r, g, bl = (round(a[k] + (b[k] - a[k]) * f) for k in range(3))
    return f"#{r:02x}{g:02x}{bl:02x}"


@dataclass
class ColourScale:
    column: str
    numeric: bool
    colors: dict[str, str]                  # tip → colour (tips without a value absent)
    categories: list[tuple[str, str]] = field(default_factory=list)  # (value, colour)
    low: float = 0.0
    high: float = 0.0
    n_missing: int = 0


def colour_scale(table: TraitTable, column: str, tips: list[str]) -> ColourScale:
    """Colours for ``tips`` from ``column`` (numeric → viridis, else categories)."""
    raw = {t: table.rows.get(t, {}).get(column, "") for t in tips}
    present = {t: v for t, v in raw.items() if str(v).strip().lower() not in MISSING}
    missing = len(tips) - len(present)
    if is_numeric(list(present.values())):
        nums = {t: _number(v) for t, v in present.items()}
        low, high = min(nums.values()), max(nums.values())
        span = high - low
        colors = {t: viridis((x - low) / span if span else 0.5) for t, x in nums.items()}
        return ColourScale(column, True, colors, [], low, high, missing)
    counts: dict[str, int] = {}
    for v in present.values():
        counts[v] = counts.get(v, 0) + 1
    ordered = sorted(counts, key=lambda v: (-counts[v], v))
    palette = {v: CATEGORICAL[i] for i, v in enumerate(ordered[: len(CATEGORICAL)])}
    colors = {t: palette.get(v, OTHER_COLOR) for t, v in present.items()}
    cats = [(v, palette[v]) for v in ordered[: len(CATEGORICAL)]]
    if len(ordered) > len(CATEGORICAL):
        cats.append((f"{len(ordered) - len(CATEGORICAL)} more values", OTHER_COLOR))
    return ColourScale(column, False, colors, cats, 0.0, 0.0, missing)


# --------------------------------------------------------------------------
# Coordinates from headers, arrays
# --------------------------------------------------------------------------

_HEADER = re.compile(r"^(?P<chrom>\S+?):(?P<start>\d[\d,]*)-(?P<end>\d[\d,]*)(?:\((?P<strand>[+-])\))?")
_SAFE = re.compile(r"^(?P<chrom>.+?)_(?P<start>\d+)-(?P<end>\d+)(?:_(?P<strand>[+-]))?$")


@dataclass
class Location:
    chrom: str
    start: int
    end: int
    strand: str = ""


def parse_location(header: str) -> Location | None:
    """``chrom:start-end(strand)`` from a FASTA header, or its tree-safe form
    ``chrom_start-end_strand`` (tree programs do not allow ``:`` and brackets)."""
    text = str(header or "").strip()
    first = text.split()[0] if text else ""
    m = _HEADER.match(first) or _SAFE.match(first)
    if not m:
        return None
    start = int(m.group("start").replace(",", ""))
    end = int(m.group("end").replace(",", ""))
    if end < start:
        start, end = end, start
    return Location(m.group("chrom"), start, end, m.group("strand") or "")


def _span_label(chrom: str, start: int, end: int) -> str:
    """``Chr1:12.40-12.91 Mb`` for long arrays; ``Chr1:12.40 Mb (1.0 kb)`` when the
    array is short compared with its position; ``Chr1:1.0-1.3 kb`` near the start."""
    if start < 1_000_000:
        return f"{chrom}:{start / 1e3:.1f}-{end / 1e3:.1f} kb"
    if end - start >= 100_000:
        return f"{chrom}:{start / 1e6:.2f}-{end / 1e6:.2f} Mb"
    return f"{chrom}:{start / 1e6:.2f} Mb ({(end - start) / 1e3:.1f} kb)"


def assign_arrays(locations: dict[str, Location], gap: int = DEFAULT_ARRAY_GAP) -> dict[str, str]:
    """Group copies into tandem arrays: same chromosome, less than ``gap`` bp apart."""
    by_chrom: dict[str, list[tuple[str, Location]]] = {}
    for tip, loc in locations.items():
        by_chrom.setdefault(loc.chrom, []).append((tip, loc))
    out: dict[str, str] = {}
    for chrom, items in by_chrom.items():
        items.sort(key=lambda it: (it[1].start, it[1].end))
        group: list[tuple[str, Location]] = []
        end = -1

        def close() -> None:
            if group:
                label = _span_label(chrom, min(g[1].start for g in group),
                                    max(g[1].end for g in group))
                for tip, _loc in group:
                    out[tip] = label

        for tip, loc in items:
            if group and loc.start - end > gap:
                close()
                group = []
            group.append((tip, loc))
            end = max(end, loc.end)
        close()
    return out


# --------------------------------------------------------------------------
# Methylation over intervals
# --------------------------------------------------------------------------


def _open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open(encoding="utf-8", errors="replace")


def methylation_format(path: Path) -> str:
    """"bismark" (.cov: chrom pos end % meth unmeth) or "bedmethyl" (modkit)."""
    with _open_text(Path(path)) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track")):
                continue
            # modkit writes columns 10-18 separated by spaces in some versions
            parts = line.split()
            if len(parts) >= 13 and parts[5] in ("+", "-", "."):
                return "bedmethyl"
            if len(parts) == 6:
                return "bismark"
            break
    raise ValueError(f"{Path(path).name} is neither a Bismark .cov nor a bedMethyl file.")


def interval_methylation(path: Path, locations: dict[str, Location], *,
                         chunk: int = 2_000_000, progress=None) -> dict[str, float | None]:
    """% methylated (reads) over each tip's interval: sum(methylated) / sum(all).

    Bismark .cov: columns 5 and 6 (methylated, unmethylated) at position col 2.
    bedMethyl (modkit): 5mC rows ("m"), N_mod (col 12) and N_canonical (col 13)
    at the 0-based start (col 2). Only chromosomes in ``locations`` are kept.
    """
    import pandas as pd

    path = Path(path)
    fmt = methylation_format(path)
    by_chrom: dict[str, list[tuple[str, int, int]]] = {}
    for tip, loc in locations.items():
        by_chrom.setdefault(loc.chrom, []).append((tip, loc.start, loc.end))
    arrays = {c: (np.array([s for _t, s, _e in iv]), np.array([e for _t, _s, e in iv]),
                  [t for t, _s, _e in iv]) for c, iv in by_chrom.items()}
    meth = {t: 0.0 for t in locations}
    total = {t: 0.0 for t in locations}
    if fmt == "bismark":
        cols, names = [0, 1, 4, 5], ["chrom", "pos", "m", "u"]
    else:
        cols, names = [0, 1, 3, 11, 12], ["chrom", "pos", "code", "m", "u"]
    reader = pd.read_csv(path, sep=r"\s+", header=None, comment="#", usecols=cols, names=names,
                         chunksize=chunk, compression="infer", dtype={"chrom": str},
                         on_bad_lines="skip", engine="c")
    for n, df in enumerate(reader):
        if progress is not None:
            progress(n)
        if "code" in df:
            df = df[df["code"].astype(str) == "m"]
        df = df[df["chrom"].isin(arrays)]
        if df.empty:
            continue
        for chrom, sub in df.groupby("chrom", sort=False):
            starts, ends, tips = arrays[chrom]
            pos = sub["pos"].to_numpy(dtype=np.int64)
            order = np.argsort(pos, kind="stable")
            pos = pos[order]
            m = np.cumsum(np.concatenate([[0.0], sub["m"].to_numpy(dtype=float)[order]]))
            a = np.cumsum(np.concatenate([[0.0], (sub["m"].to_numpy(dtype=float)
                                                 + sub["u"].to_numpy(dtype=float))[order]]))
            lo = np.searchsorted(pos, starts, side="left")
            hi = np.searchsorted(pos, ends, side="right")
            for tip, i, j in zip(tips, lo, hi):
                meth[tip] += m[j] - m[i]
                total[tip] += a[j] - a[i]
    return {t: (100.0 * meth[t] / total[t] if total[t] > 0 else None) for t in locations}


# --------------------------------------------------------------------------
# Building the table
# --------------------------------------------------------------------------


def read_names_table(path: Path | None) -> dict[str, tuple[str, str]]:
    """tree name → (original header, source file) from a tree run's names.tsv."""
    if path is None or not Path(path).is_file():
        return {}
    out = {}
    with Path(path).open(encoding="utf-8") as fh:
        next(fh, None)
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                out[parts[0]] = (parts[1], parts[2])
    return out


def build_trait_table(
    tree_path: Path,
    tips: list[str],
    *,
    methylation_file: Path | None = None,
    array_gap: int = DEFAULT_ARRAY_GAP,
    fish_files: set[str] | None = None,
    existing: TraitTable | None = None,
    progress=None,
) -> TraitTable:
    """Pre-fill a trait table for ``tips``; columns the user added are kept."""
    names = read_names_table(names_table_path(tree_path))
    headers = {t: names.get(t, (t, ""))[0] for t in tips}
    sources = {t: names.get(t, ("", ""))[1] for t in tips}
    locations = {t: loc for t in tips if (loc := parse_location(headers[t]))}
    columns: list[str] = []
    rows: dict[str, dict[str, str]] = {t: {} for t in tips}
    if any(headers[t] != t for t in tips):
        columns.append("header")
        for t in tips:
            rows[t]["header"] = headers[t]
    if locations:
        columns += ["chromosome", "start", "end", "strand", "array"]
        arrays = assign_arrays(locations, array_gap)
        for t, loc in locations.items():
            rows[t].update(chromosome=loc.chrom, start=str(loc.start), end=str(loc.end),
                           strand=loc.strand, array=arrays[t])
    if methylation_file is not None and locations:
        levels = interval_methylation(Path(methylation_file), locations, progress=progress)
        columns.append("methylation")
        for t, v in levels.items():
            rows[t]["methylation"] = "" if v is None else f"{v:.1f}"
    distinct_sources = {s for s in sources.values() if s}
    if len(distinct_sources) > 1:
        columns.append("source_file")
        for t in tips:
            rows[t]["source_file"] = sources[t]
        if fish_files is not None:
            columns.append("fish")
            for t in tips:
                rows[t]["fish"] = "yes" if sources[t] in fish_files else "no"
    if existing is not None:
        for col in existing.columns:
            if col not in columns:
                columns.append(col)
                for t in tips:
                    rows[t][col] = existing.rows.get(t, {}).get(col, "")
    return TraitTable(columns, rows)
