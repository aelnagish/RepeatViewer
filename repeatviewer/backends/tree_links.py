"""From a tree tip to its sequence (tree viewer → "Show in alignment").

Pure Python, no Qt; unit-tested in ``tests/test_tree_links.py``.

Where the sequences of a tree are found, first hit wins:

1. a RepeatViewer tree run: ``<base>.treefile`` → ``<base>_details/<base>.aln.fa``
   (older runs: ``<base>.aln.fa`` next to it);
   ``<base>_masked.treefile`` → ``<base>_details/<base>_masked.aln.fa``; trees
   saved from the viewer (``<base>.treefile_rerooted.treefile`` …) use the run's;
2. IQ-TREE's own naming: ``aln.fa.treefile`` → ``aln.fa``;
3. a FASTA with the tree's stem next to it (``x.nwk`` → ``x.fa`` / ``x.fasta`` …);
4. the run's original input FASTA (from ``names.tsv``), unaligned.

Tip names are matched to FASTA headers as the tree builder names them (first
word, made Newick-safe); IQ-TREE-style trees of other origins usually match too.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

import numpy as np

FASTA_SUFFIXES = (".fa", ".fasta", ".fas", ".faa", ".fna", ".aln", ".afa")
_SAVED_TREE = re.compile(r"^(?P<base>.+?\.treefile)_[A-Za-z0-9_]+\.treefile$")


def run_tree(tree_path: Path) -> Path:
    """The tree-run tree a viewer-saved tree came from (``x.treefile_rerooted.treefile``
    → ``x.treefile``); any other path unchanged."""
    tree_path = Path(tree_path)
    m = _SAVED_TREE.match(tree_path.name)
    return tree_path.parent / m.group("base") if m else tree_path


def run_base(tree_path: Path) -> tuple[str, bool] | None:
    """(base, masked) of a tree-run tree, or None."""
    name = run_tree(tree_path).name
    if name.endswith("_masked.treefile"):
        return name[: -len("_masked.treefile")], True
    if name.endswith(".treefile"):
        return name[: -len(".treefile")], False
    return None


@dataclass
class SequenceSource:
    path: Path
    aligned: bool          # False: the original, unaligned input
    note: str = ""         # shown in the window title when not the tree's alignment


def sequences_for_tree(tree_path: Path) -> SequenceSource | None:
    tree_path = Path(tree_path)
    folder = tree_path.parent
    rb = run_base(tree_path)
    if rb is not None:
        base, masked = rb
        details = folder / f"{base}_details"
        name = f"{base}_masked.aln.fa" if masked else f"{base}.aln.fa"
        # runs since 25 Sep 2026 keep the alignment in the details folder
        for candidate in (details / name, folder / name):
            if candidate.is_file():
                return SequenceSource(candidate, True)
        # IQ-TREE run on an alignment: aln.fa.treefile → aln.fa
        plain = folder / base
        if plain.is_file() and plain.suffix.lower() in FASTA_SUFFIXES:
            return SequenceSource(plain, True)
    stem = tree_path.name.rsplit(".", 1)[0]
    for suffix in (".aln.fa", *FASTA_SUFFIXES):
        candidate = folder / f"{stem}{suffix}"
        if candidate.is_file():
            return SequenceSource(candidate, True)
    if rb is not None:
        from repeatviewer.backends.tree_traits import read_names_table

        base = rb[0]
        names = read_names_table(folder / f"{base}_details" / f"{base}.names.tsv")
        for _header, source in names.values():
            candidate = folder / source
            if source and candidate.is_file():
                return SequenceSource(candidate, False, "input FASTA, not aligned")
    return None


def tip_key(header: str) -> str:
    """How the tree builder turns a FASTA header into a tip name (first word, safe)."""
    word = header.strip().split()[0] if header.strip() else ""
    cleaned = re.sub(r"[^A-Za-z0-9_.|\-]+", "_", word)
    return cleaned.strip("_") or "seq"


@dataclass
class Alignment:
    names: list[str]       # FASTA ids as written
    seqs: list[str]        # upper case
    aligned: bool

    @property
    def width(self) -> int:
        return max((len(s) for s in self.seqs), default=0)

    def row_of(self, tip: str) -> int | None:
        """Row of a tree tip: exact id, then the tree builder's name for it
        (duplicates numbered ``_2``, ``_3`` as the builder does)."""
        for i, n in enumerate(self.names):
            if n == tip:
                return i
        used: dict[str, int] = {}
        for i, n in enumerate(self.names):
            key = tip_key(n)
            k = used.get(key, 0)
            used[key] = k + 1
            if (key if k == 0 else f"{key}_{k + 1}") == tip:
                return i
        return None


def read_alignment(path: Path, *, aligned: bool = True) -> Alignment:
    names: list[str] = []
    parts: list[list[str]] = []
    with Path(path).open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                names.append(line[1:].strip())
                parts.append([])
            elif parts:
                parts[-1].append(line.strip())
    seqs = ["".join(p).upper() for p in parts]
    if aligned and len({len(s) for s in seqs}) > 1:
        aligned = False
    return Alignment(names, seqs, aligned)


def consensus(seqs: list[str]) -> str:
    """Most common residue per column; a gap only when gaps outnumber every
    residue. Ties between residues → first in ASCII order."""
    if not seqs:
        return ""
    width = max(len(s) for s in seqs)
    arr = np.full((len(seqs), width), ord("-"), dtype=np.uint8)
    for i, s in enumerate(seqs):
        arr[i, : len(s)] = np.frombuffer(s.encode("ascii", "replace"), dtype=np.uint8)
    out = bytearray(width)
    for j in range(width):
        counts = np.bincount(arr[:, j], minlength=256)
        gaps = int(counts[ord("-")] + counts[ord(".")])
        counts[ord("-")] = counts[ord(".")] = 0
        best = int(counts.argmax())
        out[j] = ord("-") if gaps > counts[best] else best
    return out.decode("ascii")


def first_residue(seq: str) -> int:
    """Column of the first non-gap residue (0 if none)."""
    for j, c in enumerate(seq):
        if c not in "-.":
            return j
    return 0


# --------------------------------------------------------------------------
# Distance matrices (tree viewer → Export → Distance matrix…)
# --------------------------------------------------------------------------


def patristic_matrix(model) -> tuple[list[str], np.ndarray]:
    """Tip names and the path length between every two tips of the tree."""
    from repeatviewer.backends.tree_layout import depths

    tips = model.tips()
    depth = np.array(depths(model))
    ancestors: dict[int, list[int]] = {}
    for t in tips:
        path, i = [], t
        while i is not None:
            path.append(i)
            i = model.nodes[i].parent
        ancestors[t] = path
    n = len(tips)
    out = np.zeros((n, n))
    for a in range(n):
        seen = set(ancestors[tips[a]])
        for b in range(a + 1, n):
            mrca = next(i for i in ancestors[tips[b]] if i in seen)
            out[a, b] = out[b, a] = depth[tips[a]] + depth[tips[b]] - 2 * depth[mrca]
    return [model.nodes[t].name for t in tips], out


def p_distance_matrix(aln: Alignment) -> tuple[list[str], np.ndarray]:
    """Share of differing residues over the columns both sequences have
    (gaps skipped); NaN when two sequences share no column."""
    width = aln.width
    arr = np.full((len(aln.seqs), width), ord("-"), dtype=np.uint8)
    for i, seq in enumerate(aln.seqs):
        arr[i, : len(seq)] = np.frombuffer(seq.encode("ascii", "replace"), dtype=np.uint8)
    gap = (arr == ord("-")) | (arr == ord(".")) | (arr == ord("N")) | (arr == ord("X"))
    n = len(aln.seqs)
    out = np.zeros((n, n))
    for a in range(n):
        both = ~gap[a] & ~gap[a + 1:]
        diff = (arr[a] != arr[a + 1:]) & both
        shared = both.sum(axis=1)
        with np.errstate(invalid="ignore", divide="ignore"):
            d = np.where(shared > 0, diff.sum(axis=1) / shared, np.nan)
        out[a, a + 1:] = d
        out[a + 1:, a] = d
    return list(aln.names), out


def matrix_tsv(names: list[str], matrix: np.ndarray, digits: int = 5) -> str:
    lines = ["\t" + "\t".join(names)]
    for name, row in zip(names, matrix):
        cells = ("NA" if not np.isfinite(v) else f"{v:.{digits}g}" for v in row)
        lines.append(name + "\t" + "\t".join(cells))
    return "\n".join(lines) + "\n"


# --------------------------------------------------------------------------
# Extract the sequences of a clade (tree viewer → Extract sequences)
# --------------------------------------------------------------------------


@dataclass
class Extracted:
    records: list[tuple[str, str]]   # (FASTA header, sequence)
    missing: list[str]               # tips whose sequence was not found
    note: str = ""                   # where the sequences came from, if not as asked


def _read_fasta_records(path: Path) -> list[tuple[str, str]]:
    out: list[tuple[str, list[str]]] = []
    with Path(path).open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\r\n")
            if line.startswith(">"):
                out.append((line[1:].strip(), []))
            elif out:
                out[-1][1].append(line.strip())
    return [(h, "".join(parts)) for h, parts in out]


def clade_sequences(tree_path: Path, tips: list[str], *, aligned: bool) -> Extracted:
    """Sequences of ``tips``: from the tree run's input FASTA with their original
    headers (``aligned=False``) or from the tree's alignment, gaps kept."""
    tree_path = Path(tree_path)
    if not aligned:
        rb = run_base(tree_path)
        if rb is not None:
            from repeatviewer.backends.tree_traits import read_names_table

            base = rb[0]
            names = read_names_table(tree_path.parent / f"{base}_details" / f"{base}.names.tsv")
            if names:
                files: dict[str, dict[str, str]] = {}
                records, missing = [], []
                for tip in tips:
                    header, source = names.get(tip, ("", ""))
                    path = tree_path.parent / source if source else None
                    if path is None or not path.is_file():
                        missing.append(tip)
                        continue
                    if source not in files:
                        table: dict[str, str] = {}
                        for h, seq in _read_fasta_records(path):
                            table.setdefault(h, seq)
                            table.setdefault(h.split()[0] if h.split() else h, seq)
                        files[source] = table
                    table = files[source]
                    seq = table.get(header) or table.get(header.split()[0] if header.split() else "")
                    if seq is None:
                        missing.append(tip)
                    else:
                        records.append((header, seq))
                if records or not missing:
                    return Extracted(records, missing)
        # no input FASTA to go back to: the alignment without gaps
        got = clade_sequences(tree_path, tips, aligned=True)
        return Extracted([(h, s.replace("-", "").replace(".", "")) for h, s in got.records],
                         got.missing, "from the alignment, gaps removed (input FASTA not found)")
    src = sequences_for_tree(tree_path)
    if src is None:
        return Extracted([], list(tips), "no alignment or FASTA found for this tree")
    aln = read_alignment(src.path, aligned=src.aligned)
    records, missing = [], []
    for tip in tips:
        row = aln.row_of(tip)
        if row is None:
            missing.append(tip)
        else:
            records.append((aln.names[row], aln.seqs[row]))
    return Extracted(records, missing, "" if src.aligned else "input FASTA, not aligned")


def fasta_text(records: list[tuple[str, str]], width: int = 80) -> str:
    lines = []
    for header, seq in records:
        lines.append(f">{header}")
        lines += [seq[i:i + width] for i in range(0, len(seq), width)] or [""]
    return "\n".join(lines) + "\n"


def extract_default_path(tree_path: Path, n: int, *, aligned: bool) -> Path:
    """Rule 1: ``x.treefile_7seqs.fa`` / ``x.treefile_7seqs_aligned.fa``; never
    an existing file (``_2``, ``_3`` …)."""
    tree_path = Path(tree_path)
    stem = f"{tree_path.name}_{n}seq{'s' if n != 1 else ''}{'_aligned' if aligned else ''}"
    path = tree_path.parent / f"{stem}.fa"
    k = 2
    while path.exists():
        path = tree_path.parent / f"{stem}_{k}.fa"
        k += 1
    return path
