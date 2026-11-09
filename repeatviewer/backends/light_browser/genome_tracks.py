# Vendored, read-only copy of GenomeFlow's genomeflow/backends/genome_tracks.py
# (GenomeFlow commit 5c2221b). Used by RepeatViewer's light genome browser preview.
# Changes: imports point at RepeatViewer; the gene-editing track is left out.
# Refresh by copying the file again rather than editing it here.
"""Genome browser tracks: annotations and alignments drawn under a reference.

The contig view (``widgets/contig_view.py``) shows a reference with any number
of extra tracks, like a genome browser:

* **Annotation tracks**: GFF3 / GTF / BED / VCF (``gff_features``). Exons, CDS
  and UTRs are folded into their transcript so a gene draws as one model.
* **Alignment tracks**: minimap2 PAF, BLAST tabular (``-outfmt 6/7``), SAM and
  BAM / CRAM.

Performance model
-----------------
Text formats (PAF, BLAST, SAM, GFF…) are read once into per-sequence interval
indexes (sorted numpy arrays), so a window query is two binary searches,
whatever the file size. Coverage comes from a fixed 64 bp bin profile built
once per sequence with ``np.add.at``.

BAM / CRAM are never read whole when they are indexed: with pysam (macOS /
Linux) only the visible window is fetched, padded and cached, and reads are
only fetched below ``BAM_READ_WINDOW``. Wider views show coverage from a
``.rvcov.npz`` cache built once in the background (``build_bam_coverage``), or
a hint to zoom in. Without pysam (Windows wheels do not exist) a small BAM
(< ``BAM_MEMORY_LIMIT``) is decoded by a pure-Python reader into memory, and a
larger one gets a clear message. ``samtools`` on PATH is used for CRAM.

Pure Python + numpy, no Qt; tested in ``tests/test_genome_tracks.py``.
"""

from __future__ import annotations

import gzip
import json
import re
import shutil
import struct
import zlib
from collections import OrderedDict
from collections.abc import Callable, Iterable
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

KIND_ANNOTATION = "annotation"
KIND_ALIGNMENT = "alignment"
KIND_SIGNAL = "signal"           # a value along the genome (bigWig, bedGraph)
KIND_VARIANT = "variant"         # VCF: variant sites and each sample's genotype
KIND_JUNCTION = "junction"       # splice junctions of an RNA-seq BAM's reads

COV_BIN = 64                     # bp per bin of the in-memory coverage profile
BAM_READ_WINDOW = 100_000        # wider than this: coverage only for BAM / CRAM
BAM_MEMORY_LIMIT = 300_000_000   # bytes: pure-Python BAM reader loads files up to this
BAM_FETCH_LIMIT = 50_000         # reads kept per fetched window (downsampled beyond)
SESSION_SUFFIX = ".rvgb"
SESSION_FORMAT = "repeatviewer-browser"

ANNOTATION_SUFFIXES = (".gff", ".gff3", ".gtf", ".bed", ".gff.gz", ".gff3.gz",
                       ".gtf.gz", ".bed.gz")
VCF_SUFFIXES = (".vcf", ".vcf.gz", ".vcf.bgz", ".bcf")
PAF_SUFFIXES = (".paf", ".paf.gz")
SAM_SUFFIXES = (".sam", ".sam.gz")
BAM_SUFFIXES = (".bam", ".cram")
BIGWIG_SUFFIXES = (".bw", ".bigwig")
BEDGRAPH_SUFFIXES = (".bedgraph", ".bdg", ".bedgraph.gz", ".bdg.gz")
BLAST_SUFFIXES = (".blast", ".blastn", ".outfmt6", ".m8", ".tsv", ".txt", ".tab")

_CHILD_TYPES = {"exon", "cds", "five_prime_utr", "three_prime_utr", "utr", "start_codon",
                "stop_codon", "5utr", "3utr"}


class TrackError(RuntimeError):
    """A track file that cannot be read; the message is shown to the user."""


# --------------------------------------------------------------------------
# Items and interval index
# --------------------------------------------------------------------------


@dataclass
class TrackItem:
    """One drawn interval: a feature or an alignment. Coordinates are 0-based,
    end exclusive, on the reference."""

    start: int
    end: int
    strand: str = "."
    name: str = ""
    kind: str = ""                 # feature type, or "alignment"
    score: float | None = None     # alignment identity 0..1
    mapq: int | None = None
    secondary: bool = False
    info: dict = field(default_factory=dict)
    # aligned blocks / exons: (ref_start, ref_end, query_offset); empty = one block
    blocks: list[tuple[int, int, int]] = field(default_factory=list)
    thick: list[tuple[int, int]] = field(default_factory=list)   # CDS parts
    ins: list[int] = field(default_factory=list)                 # insertion positions
    mism: list[tuple[int, str]] | None = None                    # (pos, base), lazily
    seq: str = ""                                                # query bases (SAM/BAM)
    cs: str = ""                                                 # minimap2 cs tag (PAF)

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def ref_blocks(self):
        """Aligned blocks on the reference, for drawing and depth."""
        return self.blocks or ((self.start, self.end, 0),)


class ChromIndex:
    """Intervals of one sequence, sorted by start, for fast window queries."""

    def __init__(self, items: list[TrackItem]) -> None:
        starts = np.fromiter((it.start for it in items), dtype=np.int64, count=len(items))
        ends = np.fromiter((it.end for it in items), dtype=np.int64, count=len(items))
        if len(items) > 1:
            # by start, then end; reads from a sorted BAM are usually in order already,
            # so the list is only rebuilt when it is not
            order = np.lexsort((ends, starts))
            if (order[1:] < order[:-1]).any():
                items[:] = [items[i] for i in order.tolist()]
                starts, ends = starts[order], ends[order]
        self.items = items
        self.starts = starts
        self.ends = ends
        self._max_end = np.maximum.accumulate(self.ends) if len(items) else self.ends
        self._lanes: list[int] | None = None
        self.lane_count = 0
        self._cov: np.ndarray | None = None
        self._blk: tuple[np.ndarray, np.ndarray] | None = None

    def __len__(self) -> int:
        return len(self.items)

    def query_range(self, start: int, end: int) -> tuple[int, int]:
        """Index range [lo, hi) holding every item that may overlap [start, end)."""
        hi = int(np.searchsorted(self.starts, end, side="left"))
        lo = int(np.searchsorted(self._max_end, start, side="right"))
        return lo, max(lo, hi)

    def query(self, start: int, end: int) -> list[int]:
        lo, hi = self.query_range(start, end)
        if hi <= lo:
            return []
        idx = np.nonzero(self.ends[lo:hi] > start)[0] + lo
        return idx.tolist()

    def lanes(self) -> tuple[list[int], int]:
        """Packed rows for every item (computed once, so rows stay put while scrolling)."""
        if self._lanes is None:
            self._lanes, self.lane_count = pack_rows(self.starts, self.ends)
        return self._lanes, self.lane_count

    def block_arrays(self) -> tuple[np.ndarray, np.ndarray]:
        """Starts and ends of every aligned block / exon (introns and deletions
        excluded), built once."""
        if self._blk is None:
            bs: list[int] = []
            be: list[int] = []
            for it in self.items:
                for s, e, _q in it.ref_blocks:
                    bs.append(s)
                    be.append(e)
            self._blk = (np.array(bs, dtype=np.int64), np.array(be, dtype=np.int64))
        return self._blk

    def coverage_profile(self) -> np.ndarray:
        """Depth per ``COV_BIN`` bin (mean over the bin), built once."""
        if self._cov is None:
            if not len(self.items):
                self._cov = np.zeros(1, dtype=np.float32)
            else:
                bs, be = self.block_arrays()
                self._cov = coverage_bins(bs, be, COV_BIN)
        return self._cov

    def window_coverage(self, start: int, end: int, n: int) -> np.ndarray:
        """Exact depth (max per output point) across [start, end), vectorised."""
        return _max_per_point(self.window_depth(start, end), n)

    def window_depth(self, start: int, end: int) -> np.ndarray:
        """Depth at every base of [start, end)."""
        span = max(1, end - start)
        bs, be = self.block_arrays()
        sel = (be > start) & (bs < end)
        diff = np.zeros(span + 1, dtype=np.int32)
        if sel.any():
            a = np.clip(bs[sel], start, end) - start
            b = np.clip(be[sel], start, end) - start
            np.add.at(diff, a, 1)
            np.add.at(diff, b, -1)
        return np.cumsum(diff[:span]).astype(np.float32)


def _max_per_point(values: np.ndarray, n: int) -> np.ndarray:
    """``n`` maxima over equal slices of ``values`` (each slice at least one value)."""
    if n <= 0:
        return np.zeros(0, dtype=np.float32)
    size = values.size
    if size == 0:
        return np.zeros(n, dtype=np.float32)
    idx = np.minimum((np.arange(n) * size) // n, size - 1)
    return np.maximum.reduceat(values, idx).astype(np.float32) if n <= size else \
        values[idx].astype(np.float32)


def pack_rows(starts: np.ndarray, ends: np.ndarray, gap: int = 2) -> tuple[list[int], int]:
    """Greedy first-fit rows for intervals sorted by start (genome-browser 'expanded')."""
    import heapq

    rows: list[int] = []
    free: list[tuple[int, int]] = []      # (end + gap, row) of each row's last item
    spare: list[int] = []                 # rows whose last item ended (min-heap)
    count = 0
    for s, e in zip(starts.tolist(), ends.tolist()):
        while free and free[0][0] <= s:
            _end, r = heapq.heappop(free)
            heapq.heappush(spare, r)
        if spare:
            r = heapq.heappop(spare)
        else:
            r = count
            count += 1
        rows.append(r)
        heapq.heappush(free, (e + gap, r))
    return rows, count


def _pair_span(it, chrom: str | None):
    """(pair key, mate start, approximate mate end) when ``it`` should share a row
    with its mate in "view as pairs", else None."""
    from repeatviewer.backends.light_browser import alignment_style as style

    f = style.flag(it)
    if f is None or not f & style.FLAG_PAIRED or f & (style.FLAG_MATE_UNMAPPED | 0x900):
        return None
    mchrom, mpos, _tlen = style.mate(it)
    if mpos < 0 or (mchrom is not None and chrom is not None and mchrom != chrom):
        return None
    return it.name, mpos, mpos + (it.end - it.start)


def pack_pairs(items, starts: np.ndarray, ends: np.ndarray, chrom: str | None = None,
               gap: int = 2) -> tuple[list[int], int]:
    """Like ``pack_rows`` (items sorted by start), with the two reads of a pair in one
    row and the row kept free between them."""
    import heapq

    rows: list[int] = []
    free: list[tuple[int, int]] = []
    spare: list[int] = []
    placed: dict = {}
    count = 0
    for it, s_, e_ in zip(items, starts.tolist(), ends.tolist()):
        span = _pair_span(it, chrom)
        if span is not None and span[0] in placed:
            rows.append(placed.pop(span[0]))
            continue
        while free and free[0][0] <= s_:
            _end, r = heapq.heappop(free)
            heapq.heappush(spare, r)
        if spare:
            r = heapq.heappop(spare)
        else:
            r = count
            count += 1
        rows.append(r)
        end = e_
        if span is not None and span[1] >= s_:
            end = max(e_, span[2])
            placed[span[0]] = r
        heapq.heappush(free, (end + gap, r))
    return rows, count


def coverage_bins(starts: np.ndarray, ends: np.ndarray, bin_size: int,
                  length: int | None = None) -> np.ndarray:
    """Mean depth per bin from intervals (exact base counts, vectorised)."""
    starts = np.asarray(starts, dtype=np.int64)
    ends = np.asarray(ends, dtype=np.int64)
    top = int(length if length is not None else (ends.max() if ends.size else 0))
    n = max(1, -(-top // bin_size))
    diff = np.zeros(n + 1, dtype=np.float64)
    if starts.size:
        # bases each interval contributes to its first and last bin, whole bins between
        sb = starts // bin_size
        eb = (ends - 1) // bin_size
        same = sb == eb
        np.add.at(diff, sb[same], (ends - starts)[same])
        np.add.at(diff, sb[same] + 1, -(ends - starts)[same])
        o = ~same
        if o.any():
            s, e, sbo, ebo = starts[o], ends[o], sb[o], eb[o]
            first = (sbo + 1) * bin_size - s
            last = e - ebo * bin_size
            # partial first bin
            np.add.at(diff, sbo, first)
            np.add.at(diff, sbo + 1, -first)
            # whole bins
            np.add.at(diff, sbo + 1, bin_size)
            np.add.at(diff, ebo, -bin_size)
            # partial last bin
            np.add.at(diff, ebo, last)
            np.add.at(diff, ebo + 1, -last)
    cov = np.cumsum(diff[:n]) / bin_size
    return cov.astype(np.float32)


def window_profile(profile: np.ndarray, bin_size: int, start: int, end: int,
                   n: int) -> np.ndarray:
    """``n`` values (max of the bins under each) for [start, end) from a binned profile."""
    if n <= 0 or end <= start or profile.size == 0:
        return np.zeros(max(0, n), dtype=np.float32)
    lo = max(0, start // bin_size)
    hi = min(profile.size, max(lo + 1, -(-end // bin_size)))
    if lo >= profile.size:
        return np.zeros(n, dtype=np.float32)
    return _max_per_point(profile[lo:hi], n)


def items_coverage(items: Iterable[TrackItem], start: int, end: int, n: int) -> np.ndarray:
    """Depth at ``n`` points across [start, end) from a list of items (blocks honoured)."""
    span = max(1, end - start)
    diff = np.zeros(span + 1, dtype=np.int32)
    for it in items:
        parts = it.blocks or [(it.start, it.end, 0)]
        for s, e, _q in parts:
            s, e = max(s, start), min(e, end)
            if e > s:
                diff[s - start] += 1
                diff[e - start] -= 1
    depth = np.cumsum(diff[:span])
    if n <= 0:
        return np.zeros(0, dtype=np.float32)
    edges = np.linspace(0, span, n + 1).astype(np.int64)
    out = np.zeros(n, dtype=np.float32)
    for i in range(n):
        a, b = edges[i], max(edges[i] + 1, edges[i + 1])
        out[i] = depth[a:min(b, span)].max() if a < span else 0
    return out


# --------------------------------------------------------------------------
# Sequence-name matching
# --------------------------------------------------------------------------


def match_chrom(names: Iterable[str], wanted: str) -> str | None:
    """Track sequence name for reference ``wanted`` (exact, then chr/Chr/version-insensitive)."""
    names = list(names)
    if wanted in names:
        return wanted

    def norm(s: str) -> str:
        s = s.split()[0] if s.split() else s
        s = re.sub(r"^(chr|chromosome|scaffold)[_\-]?", "", s, flags=re.IGNORECASE)
        return re.sub(r"\.\d+$", "", s).lower()

    target = norm(wanted)
    hits = [n for n in names if norm(n) == target]
    return hits[0] if len(hits) == 1 else None


# --------------------------------------------------------------------------
# Tracks
# --------------------------------------------------------------------------


class Track:
    """A file shown as a track. Subclasses fill ``_chroms`` or override fetch."""

    kind = KIND_ALIGNMENT
    format = ""
    read_window: int | None = None     # reads are fetched only below this window

    def __init__(self, path: Path, name: str | None = None) -> None:
        self.path = Path(path)
        self.name = name or self.path.name
        self._chroms: dict[str, ChromIndex] = {}
        self._alias: dict[str, str | None] = {}
        self.warnings: list[str] = []
        self.settings: dict = {}       # display settings saved with a session
        self.read_groups: dict[str, dict[str, str]] = {}   # @RG ID -> {SM, LB, PL…}
        self._insert_ranges: dict[str, tuple[int, int] | None] = {}
        self._group_key = None          # what the rows are grouped by (None: no groups)
        self._group_fn = None
        self._group_lanes: dict = {}
        self._filter_key = None         # which reads are hidden (None: all shown)
        self._filter_fn = None
        self._filtered: dict = {}
        self.pairs = False              # mates in one row ("view as pairs")
        self._pair_lanes: dict = {}

    # names
    def chrom_names(self) -> list[str]:
        return list(self._chroms)

    def resolve(self, ref_id: str) -> str | None:
        if ref_id not in self._alias:
            self._alias[ref_id] = match_chrom(self.chrom_names(), ref_id)
        return self._alias[ref_id]

    def index(self, ref_id: str) -> ChromIndex | None:
        key = self.resolve(ref_id)
        ix = self._chroms.get(key) if key is not None else None
        if ix is None or self._filter_fn is None:
            return ix
        got = self._filtered.get(key)
        if got is None:
            keep = self._filter_fn
            got = self._filtered[key] = ChromIndex([it for it in ix.items if keep(it)])
        return got

    def set_filter(self, key, fn) -> None:
        """Hide the reads for which ``fn(item)`` is False (``key`` names the setting)."""
        key = key if fn is not None else None
        if key == self._filter_key:
            return
        self._filter_key, self._filter_fn = key, fn
        self._filtered = {}
        self._group_lanes = {}
        self._pair_lanes = {}
        self._layout_changed()

    def set_pairs(self, on: bool) -> None:
        if bool(on) == self.pairs:
            return
        self.pairs = bool(on)
        self._group_lanes = {}
        self._pair_lanes = {}
        self._layout_changed()

    def _layout_changed(self) -> None:
        self._grouping_changed()

    def _rows_of(self, ix: ChromIndex, ref_id: str) -> tuple[list[int], int]:
        if not self.pairs:
            return ix.lanes()
        got = self._pair_lanes.get(id(ix))
        if got is None:
            got = self._pair_lanes[id(ix)] = pack_pairs(ix.items, ix.starts, ix.ends,
                                                        self.resolve(ref_id))
        return got

    def count(self, ref_id: str) -> int:
        ix = self.index(ref_id)
        return len(ix) if ix is not None else 0

    # groups: reads split by a property (strand, read group, tag…), each group with
    # its own rows

    def set_grouping(self, key, fn) -> None:
        """Group by ``fn(item) -> label`` (``key`` names the setting); None: no groups."""
        key = key if fn is not None else None
        if key == self._group_key:
            return
        self._group_key, self._group_fn = key, fn
        self._group_lanes = {}
        self._grouping_changed()

    def _grouping_changed(self) -> None:
        pass

    @property
    def grouped(self) -> bool:
        return self._group_fn is not None

    def _lanes_for(self, ix: ChromIndex, chrom: str | None = None):
        got = self._group_lanes.get(id(ix))
        if got is None:
            fn = self._group_fn
            labels = [fn(it) for it in ix.items]
            rows = [0] * len(labels)
            by: dict = {}
            for i, g in enumerate(labels):
                by.setdefault(g, []).append(i)
            counts = {}
            for g, idx in by.items():
                arr = np.asarray(idx)
                if self.pairs:
                    lanes, n = pack_pairs([ix.items[i] for i in idx], ix.starts[arr],
                                          ix.ends[arr], chrom)
                else:
                    lanes, n = pack_rows(ix.starts[arr], ix.ends[arr])
                for i, r in zip(idx, lanes):
                    rows[i] = r
                counts[g] = n
            got = self._group_lanes[id(ix)] = (labels, rows, counts)
        return got

    def fetch_groups(self, ref_id: str, start: int, end: int) -> list[tuple]:
        """(item, group label, row within the group) overlapping [start, end)."""
        ix = self.index(ref_id)
        if ix is None or not self.grouped:
            return []
        labels, rows, _n = self._lanes_for(ix, self.resolve(ref_id))
        return [(ix.items[i], labels[i], rows[i]) for i in ix.query(start, end)]

    def group_rows(self, ref_id: str) -> dict:
        """Group label -> rows it uses on this sequence."""
        ix = self.index(ref_id)
        if ix is None or not self.grouped:
            return {}
        return self._lanes_for(ix, self.resolve(ref_id))[2]

    def insert_range(self, ref_id: str) -> tuple[int, int] | None:
        """Expected insert size (0.5–99.5 percentile of proper pairs) on this sequence,
        used to colour pairs that are too far apart or too close."""
        if ref_id not in self._insert_ranges:
            from repeatviewer.backends.light_browser.alignment_style import insert_bounds

            self._insert_ranges[ref_id] = insert_bounds(self._insert_sample(ref_id))
        return self._insert_ranges[ref_id]

    def _insert_sample(self, ref_id: str):
        ix = self.index(ref_id)
        if ix is None:
            return []
        step = max(1, len(ix.items) // 20000)
        return ix.items[::step]

    # queries
    def fetch(self, ref_id: str, start: int, end: int) -> list[tuple[TrackItem, int]]:
        """(item, row) overlapping [start, end); rows are stable for the sequence."""
        ix = self.index(ref_id)
        if ix is None:
            return []
        rows, _n = self._rows_of(ix, ref_id)
        return [(ix.items[i], rows[i]) for i in ix.query(start, end)]

    def row_count(self, ref_id: str) -> int:
        ix = self.index(ref_id)
        if ix is None:
            return 0
        return self._rows_of(ix, ref_id)[1]

    def coverage(self, ref_id: str, start: int, end: int, n: int) -> np.ndarray | None:
        ix = self.index(ref_id)
        if ix is None:
            return None
        if (end - start) / max(1, n) < COV_BIN:
            return ix.window_coverage(start, end, n)
        return window_profile(ix.coverage_profile(), COV_BIN, start, end, n)

    def count_in(self, ref_id: str, start: int, end: int) -> int:
        """Upper bound of the items overlapping [start, end) (no item is read)."""
        ix = self.index(ref_id)
        if ix is None:
            return 0
        lo, hi = ix.query_range(start, end)
        return hi - lo

    def layout_version(self, ref_id: str, start: int, end: int):
        """Changes when the rows of the items move (drawn tiles are then redrawn)."""
        return 0

    def tile_complete(self, ref_id: str, start: int, end: int) -> bool:
        """False when [start, end) may still miss items (a drawn tile is then not kept)."""
        return True

    def junctions(self, ref_id: str, start: int, end: int) -> np.ndarray:
        """Splice junctions (introns, N in CIGAR) of the reads that overlap
        [start, end): rows of (start, end, reads on +, reads on -, reads with no
        strand), the strand from the XS / TS tag."""
        ix = self.index(ref_id)
        return _window_junctions([ix] if ix is not None else [], start, end)

    def coverage_peak(self, ref_id: str, start: int, end: int) -> float:
        cov = self.coverage(ref_id, start, end, 256)
        return float(cov.max()) if cov is not None and len(cov) else 0.0

    def reads_visible(self, start: int, end: int) -> bool:
        return self.read_window is None or end - start <= self.read_window

    def find(self, text: str, limit: int = 50) -> list[tuple[str, TrackItem]]:
        """Items whose name contains ``text`` (case-insensitive), for Go to."""
        t = text.lower()
        out = []
        for chrom, ix in self._chroms.items():
            for it in ix.items:
                if t in it.name.lower():
                    out.append((chrom, it))
                    if len(out) >= limit:
                        return out
        return out

    def search(self, text: str, limit: int = 200) -> tuple[list[tuple[str, TrackItem]], int]:
        """Features whose name, type or any attribute (ID, gene_id, product, Note…)
        contains ``text``, case-insensitive: (up to ``limit`` of them, how many in
        all). Names that match come first. Only annotation tracks have features."""
        if self.kind != KIND_ANNOTATION:
            return [], 0
        q = text.strip().lower()
        if not q:
            return [], 0
        by_name: list = []
        other: list = []
        for chrom, ix in getattr(self, "_chroms", {}).items():
            for it in ix.items:
                if q in it.name.lower():
                    by_name.append((chrom, it))
                elif q in it.kind.lower() or any(
                        q in str(v).lower() for k, v in it.info.items()
                        if not str(k).startswith("_")):
                    other.append((chrom, it))
        # exact names first, then names that start with it, then the rest
        by_name.sort(key=lambda ci: (ci[1].name.lower() != q,
                                     not ci[1].name.lower().startswith(q)))
        hits = by_name + other
        return hits[:limit], len(hits)

    def to_json(self) -> dict:
        return {"path": str(self.path), "name": self.name, "settings": dict(self.settings)}


def _open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open(encoding="utf-8", errors="replace")


class AnnotationTrack(Track):
    kind = KIND_ANNOTATION
    format = "annotation"

    def __init__(self, path: Path, name: str | None = None) -> None:
        super().__init__(path, name)
        from repeatviewer.backends.gff_features import parse_annotation

        p = self.path
        if p.name.lower().endswith(".gz") and not p.name.lower().endswith(".vcf.gz"):
            feats = _parse_gz_annotation(p)
        else:
            feats = parse_annotation(p)
        if not feats:
            raise TrackError(f"No features found in {p.name}.")
        by_chrom: dict[str, list[TrackItem]] = {}
        for chrom, item in fold_gene_models(feats):
            by_chrom.setdefault(chrom, []).append(item)
        self._chroms = {c: ChromIndex(v) for c, v in by_chrom.items()}

    def feature_types(self) -> list[str]:
        seen: dict[str, None] = {}
        for ix in self._chroms.values():
            for it in ix.items:
                seen.setdefault(it.kind, None)
        return list(seen)


def _parse_gz_annotation(path: Path):
    """GFF / GTF / BED compressed with gzip (the plain parsers read text files)."""
    import tempfile

    from repeatviewer.backends.gff_features import parse_annotation

    inner = path.name[:-3]
    with tempfile.TemporaryDirectory() as tmp:
        target = Path(tmp) / inner
        with gzip.open(path, "rb") as src, target.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        return parse_annotation(target)


def fold_gene_models(feats) -> list[tuple[str, TrackItem]]:
    """GffFeatures → (seqid, item). Exons / CDS / UTRs become blocks of their parent
    transcript (GFF3 ``Parent`` or GTF ``transcript_id``); the rest stay as they are."""
    parents: dict[tuple[str, str], TrackItem] = {}
    out: list[tuple[str, TrackItem]] = []
    children = []
    for f in feats:
        ftype = str(f.feature_type or "feature")
        attrs = f.attributes or {}
        if ftype.lower() in _CHILD_TYPES:
            children.append(f)
            continue
        item = TrackItem(f.start - 1, f.end, f.strand if f.strand in "+-" else ".",
                         f.label(), ftype, info=dict(attrs))
        out.append((f.seqid, item))
        for key in (attrs.get("ID"), attrs.get("transcript_id")):
            if key and ftype.lower() not in ("gene",):
                parents[(f.seqid, key)] = item
    orphans: dict[tuple[str, str], TrackItem] = {}
    for f in children:
        attrs = f.attributes or {}
        keys = [k for k in (attrs.get("Parent", "").split(",") + [attrs.get("transcript_id", "")])
                if k]
        targets = [parents[(f.seqid, k)] for k in keys if (f.seqid, k) in parents]
        if not targets and attrs.get("transcript_id"):
            # GTF without transcript lines: build the transcript from its parts
            key = (f.seqid, attrs["transcript_id"])
            t = orphans.get(key)
            if t is None:
                t = TrackItem(f.start - 1, f.end, f.strand if f.strand in "+-" else ".",
                              attrs.get("gene_name") or attrs["transcript_id"], "transcript",
                              info={k: v for k, v in attrs.items() if k != "exon_number"})
                orphans[key] = t
                out.append((f.seqid, t))
            t.start, t.end = min(t.start, f.start - 1), max(t.end, f.end)
            targets = [t]
        if not targets:
            out.append((f.seqid, TrackItem(f.start - 1, f.end,
                                           f.strand if f.strand in "+-" else ".",
                                           f.label(), str(f.feature_type), info=dict(attrs))))
            continue
        for t in targets:
            if f.feature_type.lower() == "cds":
                t.thick.append((f.start - 1, f.end))
            else:
                t.blocks.append((f.start - 1, f.end, 0))
    for _c, t in out:
        if t.thick and not t.blocks:
            t.blocks = [(s, e, 0) for s, e in t.thick]
        if t.blocks:
            t.blocks = _merge_blocks(t.blocks)
        t.thick.sort()
    return out


def _merge_blocks(blocks: list[tuple[int, int, int]]) -> list[tuple[int, int, int]]:
    blocks = sorted(blocks)
    merged = [list(blocks[0])]
    for s, e, q in blocks[1:]:
        if s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e, q])
    return [tuple(b) for b in merged]


# ---------------------------------------------------------------- PAF


_CS_TOKEN = re.compile(r"(=[ACGTNacgtn]+|:[0-9]+|\*[a-z][a-z]|\+[a-z]+|-[a-z]+|~[a-z]{2}[0-9]+[a-z]{2})")


def parse_paf_record(line: str) -> tuple[str, TrackItem] | None:
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 12:
        return None
    try:
        qlen = int(cols[1])
        qs, qe = int(cols[2]), int(cols[3])
        ts, te = int(cols[7]), int(cols[8])
        nmatch, alen, mapq = int(cols[9]), int(cols[10]), int(cols[11])
    except ValueError:
        return None
    tags = {}
    for t in cols[12:]:
        if len(t) > 5 and t[2] == ":" and t[4] == ":":
            tags[t[:2]] = t[5:]
    item = TrackItem(ts, te, cols[4], cols[0], "alignment",
                     score=nmatch / alen if alen else None, mapq=mapq,
                     secondary=tags.get("tp", "P") not in ("P", "I"),
                     info={"query length": qlen, "query": f"{qs + 1:,}–{qe:,}",
                           "matches": nmatch, "aligned": alen},
                     cs=tags.get("cs", ""))
    if "NM" in tags:
        item.info["edit distance"] = int(tags["NM"])
    head, tail = qs, qlen - qe          # unaligned query ends (minimap2 clips them)
    clip = (head, tail) if cols[4] == "+" else (tail, head)
    if clip != (0, 0):
        item.info["_clip"] = clip
    return cols[5], item


def cs_details(item: TrackItem) -> None:
    """Mismatches, deletions (as blocks) and insertions from a cs tag, once."""
    if item.mism is not None:
        return
    item.mism = []
    if not item.cs:
        return
    pos = item.start
    blocks: list[tuple[int, int, int]] = []
    block_start = pos
    for tok in _CS_TOKEN.findall(item.cs):
        op, body = tok[0], tok[1:]
        if op == "=":
            pos += len(body)
        elif op == ":":
            pos += int(body)
        elif op == "*":
            item.mism.append((pos, body[1].upper()))
            pos += 1
        elif op == "+":
            item.ins.append(pos)
            item.info.setdefault("_ins", []).append((pos, len(body)))
        elif op in "-~":
            n = len(body) if op == "-" else int(body[2:-2])
            if pos > block_start:
                blocks.append((block_start, pos, 0))
            pos += n
            block_start = pos
    if pos > block_start:
        blocks.append((block_start, pos, 0))
    if len(blocks) > 1:
        item.blocks = blocks


class PafTrack(Track):
    format = "paf"

    def __init__(self, path: Path, name: str | None = None) -> None:
        super().__init__(path, name)
        by_chrom: dict[str, list[TrackItem]] = {}
        with _open_text(self.path) as fh:
            for line in fh:
                rec = parse_paf_record(line)
                if rec is not None:
                    by_chrom.setdefault(rec[0], []).append(rec[1])
        if not by_chrom:
            raise TrackError(f"No alignments found in {self.path.name} (not a PAF file?).")
        self._chroms = {c: ChromIndex(v) for c, v in by_chrom.items()}


# ---------------------------------------------------------------- BLAST

BLAST_DEFAULT_COLUMNS = ["qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
                         "qstart", "qend", "sstart", "send", "evalue", "bitscore"]


def looks_like_blast_table(path: Path, lines: int = 20) -> bool:
    try:
        with _open_text(path) as fh:
            for _ in range(lines):
                line = fh.readline()
                if not line:
                    break
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12:
                    return False
                float(cols[2])
                [int(c) for c in cols[6:10]]
                float(cols[10])
                return True
    except (OSError, ValueError):
        return False
    return False


class BlastTrack(Track):
    """BLAST tabular output. ``genome_side`` says which column holds the reference
    sequence: ``subject`` (TEs searched against a genome; the default), ``query``,
    or ``auto`` (whichever side's names match ``reference_names``)."""

    format = "blast"

    def __init__(self, path: Path, name: str | None = None, *, genome_side: str = "auto",
                 reference_names: Iterable[str] = ()) -> None:
        super().__init__(path, name)
        cols_names = list(BLAST_DEFAULT_COLUMNS)
        rows: list[list[str]] = []
        with _open_text(self.path) as fh:
            for line in fh:
                if line.startswith("#"):
                    if line.lower().startswith("# fields:"):
                        cols_names = _blast_fields(line)
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) >= 12:
                    rows.append(cols)
        if not rows:
            raise TrackError(f"No hits found in {self.path.name} (not BLAST tabular output?).")
        col = {n: i for i, n in enumerate(cols_names)}
        need = ("qseqid", "sseqid", "pident", "qstart", "qend", "sstart", "send")
        if any(k not in col for k in need):
            col = {n: i for i, n in enumerate(BLAST_DEFAULT_COLUMNS)}
        side = genome_side
        if side == "auto":
            refs = set(reference_names)
            subj = sum(1 for r in rows[:2000] if r[col["sseqid"]] in refs)
            qry = sum(1 for r in rows[:2000] if r[col["qseqid"]] in refs)
            side = "query" if qry > subj else "subject"
        self.genome_side = side
        on, other = ("sseqid", "qseqid") if side == "subject" else ("qseqid", "sseqid")
        s_key, e_key = ("sstart", "send") if side == "subject" else ("qstart", "qend")
        o_s, o_e = ("qstart", "qend") if side == "subject" else ("sstart", "send")
        by_chrom: dict[str, list[TrackItem]] = {}
        for r in rows:
            try:
                a, b = int(r[col[s_key]]), int(r[col[e_key]])
                oa, ob = int(r[col[o_s]]), int(r[col[o_e]])
                pident = float(r[col["pident"]])
            except (ValueError, IndexError):
                continue
            rev = (a > b) != (oa > ob)
            s, e = min(a, b) - 1, max(a, b)
            info = {"hit on": f"{min(oa, ob):,}–{max(oa, ob):,}",
                    "identity": f"{pident:.1f}%"}
            for k in ("length", "evalue", "bitscore", "mismatch", "gapopen"):
                if k in col and col[k] < len(r):
                    info[k] = r[col[k]]
            by_chrom.setdefault(r[col[on]], []).append(
                TrackItem(s, e, "-" if rev else "+", r[col[other]], "alignment",
                          score=pident / 100.0, info=info))
        self._chroms = {c: ChromIndex(v) for c, v in by_chrom.items()}


def _blast_fields(line: str) -> list[str]:
    names = {"query acc.ver": "qseqid", "query id": "qseqid", "subject acc.ver": "sseqid",
             "subject id": "sseqid", "% identity": "pident", "alignment length": "length",
             "mismatches": "mismatch", "gap opens": "gapopen", "q. start": "qstart",
             "q. end": "qend", "s. start": "sstart", "s. end": "send", "evalue": "evalue",
             "bit score": "bitscore"}
    fields = line.split(":", 1)[1].split(",")
    return [names.get(f.strip().lower(), f.strip().lower()) for f in fields]


# ---------------------------------------------------------------- SAM / BAM

_CIGAR = re.compile(r"(\d+)([MIDNSHP=X])")
_CIGAR_OPS = "MIDNSHP=X"


def cigar_layout(pos: int, ops: Iterable[tuple[str, int]]) -> tuple[int, list, list[int], int]:
    """(end, blocks (ref_start, ref_end, query_offset), insertion positions, matched bases)."""
    ref = pos
    q = 0
    blocks: list[tuple[int, int, int]] = []
    ins: list[int] = []
    aligned = 0
    for op, n in ops:
        if op in "M=X":
            if blocks and blocks[-1][1] == ref and blocks[-1][2] + (ref - blocks[-1][0]) == q:
                s0, _e0, q0 = blocks[-1]
                blocks[-1] = (s0, ref + n, q0)
            else:
                blocks.append((ref, ref + n, q))
            ref += n
            q += n
            aligned += n
        elif op == "I":
            ins.append(ref)
            q += n
        elif op == "S":
            q += n
        elif op in "DN":
            ref += n
    return ref, blocks, ins, aligned


def cigar_extras(pos: int, ops) -> tuple[list[tuple[int, int]], tuple[int, int]]:
    """(insertions as (reference position, length), (left, right) soft-clip lengths)."""
    ins: list[tuple[int, int]] = []
    ref = pos
    left = right = 0
    seen_aligned = False
    for op, n in ops:
        if op in "M=XDN":
            ref += n
            seen_aligned = True
        elif op == "I":
            ins.append((ref, n))
        elif op == "S":
            if seen_aligned:
                right += n
            else:
                left += n
    return ins, (left, right)


def _nm_identity(nm: int | None, aligned: int, ins_dels: int) -> float | None:
    if nm is None or aligned + ins_dels <= 0:
        return None
    return max(0.0, 1.0 - nm / (aligned + ins_dels))


def parse_read_groups(header_lines: Iterable[str]) -> dict[str, dict[str, str]]:
    """@RG lines of a SAM header -> {ID: {SM: …, LB: …}}."""
    groups: dict[str, dict[str, str]] = {}
    for line in header_lines:
        if not line.startswith("@RG"):
            continue
        fields = dict(f.split(":", 1) for f in line.rstrip("\n").split("\t")[1:] if ":" in f)
        if "ID" in fields:
            groups[fields.pop("ID")] = fields
    return groups


def parse_sam_record(line: str, keep_seq: bool = True) -> tuple[str, TrackItem] | None:
    if not line or line.startswith("@"):
        return None
    cols = line.rstrip("\n").split("\t")
    if len(cols) < 11:
        return None
    try:
        flag = int(cols[1])
        pos = int(cols[3]) - 1
        mapq = int(cols[4])
    except ValueError:
        return None
    if flag & 4 or cols[2] == "*" or cols[5] == "*":
        return None
    ops = [(op, int(n)) for n, op in _CIGAR.findall(cols[5])]
    end, blocks, ins, aligned = cigar_layout(pos, ops)
    tags = {t[:2]: t[5:] for t in cols[11:] if len(t) > 5}
    indel = sum(n for op, n in ops if op in "ID")
    nm = int(tags["NM"]) if "NM" in tags and tags["NM"].lstrip("-").isdigit() else None
    item = TrackItem(pos, end, "-" if flag & 16 else "+", cols[0], "alignment",
                     score=_nm_identity(nm, aligned, indel), mapq=mapq,
                     secondary=bool(flag & 0x900),
                     info={"flag": flag, "cigar": cols[5] if len(cols[5]) < 60 else
                           cols[5][:57] + "…",
                           "_rnext": cols[2] if cols[6] == "=" else cols[6],
                           "_pnext": int(cols[7]) - 1 if cols[7].isdigit() else -1,
                           "_tlen": int(cols[8]) if cols[8].lstrip("-").isdigit() else 0,
                           "_tags": {t[:2]: t[5:] for t in cols[11:] if len(t) > 5}},
                     blocks=blocks if len(blocks) > 1 or (blocks and blocks[0][2]) else [],
                     ins=ins, seq=cols[9] if keep_seq and cols[9] != "*" else "")
    if not item.blocks and blocks:
        item.blocks = blocks
    _put_extras(item, pos, ops, cols[10] if keep_seq and cols[10] != "*" else None)
    return cols[2], item


def _put_extras(item: TrackItem, pos: int, ops, qual) -> None:
    """Insertion lengths, soft clips and base qualities (only when present)."""
    ins, clip = cigar_extras(pos, ops)
    if any(n > 1 for _p, n in ins):
        item.info["_ins"] = ins
    if clip != (0, 0):
        item.info["_clip"] = clip
    if qual:
        item.info["_qual"] = qual


def mismatches_against(item: TrackItem, reference: str, start: int, end: int, *,
                       offsets: bool = False) -> list[tuple]:
    """Bases of ``item`` differing from ``reference`` within [start, end): (position,
    base), or (position, base, query offset or -1) with ``offsets``."""
    if item.mism is not None and not item.seq:
        hits = [m for m in item.mism if start <= m[0] < end]
        return [(m[0], m[1], -1) for m in hits] if offsets else hits
    if not item.seq or not reference:
        return []
    refb = _reference_codes(reference)
    seqb = np.frombuffer(item.seq.upper().encode("ascii", "replace"), dtype=np.uint8)
    n_ref = len(refb)
    out = []
    for s, e, q in item.blocks or [(item.start, item.end, 0)]:
        lo, hi = max(s, start), min(e, end, n_ref)
        if hi <= lo:
            continue
        k0 = q + lo - s
        k1 = min(k0 + hi - lo, len(seqb))
        if k1 <= k0:
            continue
        qa = seqb[k0:k1]
        ra = refb[lo:lo + (k1 - k0)]
        hit = np.nonzero((qa != ra) & (qa != 78) & (ra != 78))[0]
        for j in hit.tolist():
            base = chr(qa[j])
            out.append((lo + j, base, k0 + j) if offsets else (lo + j, base))
    return out


_REF_CODES: list = [None, None]      # (reference str, its upper-case bytes)


def _reference_codes(reference: str) -> np.ndarray:
    if _REF_CODES[0] is not reference:
        _REF_CODES[0] = reference
        _REF_CODES[1] = np.frombuffer(reference.upper().encode("ascii", "replace"),
                                      dtype=np.uint8)
    return _REF_CODES[1]


class SamTrack(Track):
    format = "sam"

    SEQ_BUDGET = 300_000_000   # query bases kept in memory for mismatch drawing

    def __init__(self, path: Path, name: str | None = None) -> None:
        super().__init__(path, name)
        by_chrom: dict[str, list[TrackItem]] = {}
        kept = 0
        header: list[str] = []
        with _open_text(self.path) as fh:
            for line in fh:
                if line.startswith("@"):
                    header.append(line)
                    continue
                rec = parse_sam_record(line, keep_seq=kept < self.SEQ_BUDGET)
                if rec is not None:
                    kept += len(rec[1].seq)
                    by_chrom.setdefault(rec[0], []).append(rec[1])
        if kept >= self.SEQ_BUDGET:
            self.warnings.append("Large file: differences are shown for the first reads only.")
        if not by_chrom:
            raise TrackError(f"No mapped reads found in {self.path.name}.")
        self._chroms = {c: ChromIndex(v) for c, v in by_chrom.items()}
        self.read_groups = parse_read_groups(header)


# BAM


def pysam_available() -> bool:
    try:
        import pysam  # noqa: F401
    except Exception:  # noqa: BLE001 - missing or broken build
        return False
    return True


def bam_index_path(path: Path) -> Path | None:
    for cand in (path.with_name(path.name + ".bai"), path.with_suffix(".bai"),
                 path.with_name(path.name + ".csi"), path.with_name(path.name + ".crai")):
        if cand.is_file():
            return cand
    return None


def coverage_cache_path(path: Path) -> Path:
    return path.with_name(path.name + ".rvcov.npz")


class _WindowCache:
    def __init__(self, size: int = 6) -> None:
        self._data: OrderedDict = OrderedDict()
        self.size = size

    def get(self, key):
        if key in self._data:
            self._data.move_to_end(key)
            return self._data[key]
        return None

    def put(self, key, value) -> None:
        self._data[key] = value
        self._data.move_to_end(key)
        while len(self._data) > self.size:
            self._data.popitem(last=False)


def _pack_right(starts, ends, last_end: list[int], first_start: list[int],
                gap: int = 2) -> list[int]:
    """Rows for reads sorted by start that all start right of the rows already used:
    the earlier rows keep their reads (so drawn tiles stay valid)."""
    import heapq

    free = [(le + gap, r) for r, le in enumerate(last_end)]
    heapq.heapify(free)
    spare: list[int] = []
    rows = []
    for s_, e_ in zip(starts.tolist(), ends.tolist()):
        while free and free[0][0] <= s_:
            heapq.heappush(spare, heapq.heappop(free)[1])
        if spare:
            r = heapq.heappop(spare)
        else:
            r = len(last_end)
            last_end.append(e_)
            first_start.append(s_)
        last_end[r] = max(last_end[r], e_)
        rows.append(r)
        heapq.heappush(free, (e_ + gap, r))
    return rows


def _pack_left(starts, ends, last_end: list[int], first_start: list[int],
               gap: int = 2) -> list[int]:
    """Mirror of ``_pack_right`` for reads added on the left of the used rows."""
    import heapq

    s_list, e_list = starts.tolist(), ends.tolist()
    free = [(-(fs - gap), r) for r, fs in enumerate(first_start)]
    heapq.heapify(free)
    spare: list[int] = []
    rows = [0] * len(s_list)
    for i in sorted(range(len(s_list)), key=lambda k: -e_list[k]):
        s_, e_ = s_list[i], e_list[i]
        while free and -free[0][0] >= e_:
            heapq.heappush(spare, heapq.heappop(free)[1])
        if spare:
            r = heapq.heappop(spare)
        else:
            r = len(first_start)
            last_end.append(e_)
            first_start.append(s_)
        first_start[r] = min(first_start[r], s_)
        rows[i] = r
        heapq.heappush(free, (-(s_ - gap), r))
    return rows


def _pack_right_pairs(items, starts, ends, last_end: list[int], first_start: list[int],
                      placed: dict, chrom: str | None, gap: int = 2) -> list[int]:
    """``_pack_right`` for "view as pairs": a read whose mate lies to its right keeps
    the row free up to the mate, and the mate then takes that row."""
    import heapq

    free = [(le + gap, r) for r, le in enumerate(last_end)]
    heapq.heapify(free)
    spare: list[int] = []
    rows = []
    for it, s_, e_ in zip(items, starts.tolist(), ends.tolist()):
        span = _pair_span(it, chrom)
        if span is not None and span[0] in placed:
            r = placed.pop(span[0])
            rows.append(r)
            last_end[r] = max(last_end[r], e_)
            continue
        while free and free[0][0] <= s_:
            heapq.heappush(spare, heapq.heappop(free)[1])
        if spare:
            r = heapq.heappop(spare)
        else:
            r = len(last_end)
            last_end.append(e_)
            first_start.append(s_)
        end = e_
        if span is not None and span[1] >= s_:
            end = max(e_, span[2])
            placed[span[0]] = r
        last_end[r] = max(last_end[r], end)
        rows.append(r)
        heapq.heappush(free, (end + gap, r))
    return rows


def _pack_left_pairs(items, starts, ends, last_end: list[int], first_start: list[int],
                     placed: dict, chrom: str | None, gap: int = 2) -> list[int]:
    """Mirror of ``_pack_right_pairs`` for tiles added on the left."""
    import heapq

    s_list, e_list = starts.tolist(), ends.tolist()
    free = [(-(fs - gap), r) for r, fs in enumerate(first_start)]
    heapq.heapify(free)
    spare: list[int] = []
    rows = [0] * len(s_list)
    for i in sorted(range(len(s_list)), key=lambda k: -e_list[k]):
        it, s_, e_ = items[i], s_list[i], e_list[i]
        span = _pair_span(it, chrom)
        if span is not None and span[0] in placed:
            r = placed.pop(span[0])
            if first_start[r] >= e_ + gap or first_start[r] >= s_:
                rows[i] = r
                first_start[r] = min(first_start[r], s_)
                continue
        while free and -free[0][0] >= e_:
            heapq.heappush(spare, heapq.heappop(free)[1])
        if spare:
            r = heapq.heappop(spare)
        else:
            r = len(first_start)
            last_end.append(e_)
            first_start.append(s_)
        start = s_
        if span is not None and span[1] < s_:
            start = min(s_, span[1])
            placed[span[0]] = r
        first_start[r] = min(first_start[r], start)
        rows[i] = r
        heapq.heappush(free, (-(start - gap), r))
    return rows


class _BamRegion:
    """Reads of consecutive BAM tiles with stable rows. It grows tile by tile at
    either end while panning, so rows already drawn never move."""

    MAX_TILES = 32

    def __init__(self, key: str, serial: int, tile: int, group_fn=None,
                 pairs: bool = False) -> None:
        self.key = key
        self.serial = serial
        self.tile = tile
        self.group_fn = group_fn
        self.pairs = pairs
        self.placed: dict = {}                # pair key -> row kept for the mate
        self.group_placed: dict = {}
        self.groups: dict = {}                # label -> (last_end, first_start)
        self.lo = self.hi = 0                 # tiles [lo, hi)
        self.chunks: dict[int, ChromIndex] = {}
        self.last_end: list[int] = []
        self.first_start: list[int] = []
        self.max_len = 0

    def _chunk(self, items: list[TrackItem]) -> ChromIndex:
        ix = ChromIndex(items)
        pre = getattr(items, "prepared", None)
        if pre is not None:
            # block arrays and depth were worked out by the worker that read the tile
            ix._blk, ix._depth = pre
        if len(ix):
            self.max_len = max(self.max_len, int((ix.ends - ix.starts).max()))
        return ix

    def add_right(self, t: int, items: list[TrackItem]) -> None:
        ix = self._chunk(items)
        if not self.chunks:
            self.lo = t
        if self.pairs:
            ix._lanes = _pack_right_pairs(ix.items, ix.starts, ix.ends, self.last_end,
                                          self.first_start, self.placed, self.key)
        else:
            ix._lanes = _pack_right(ix.starts, ix.ends, self.last_end, self.first_start)
        self._pack_groups(ix, _pack_right_pairs if self.pairs else _pack_right)
        self.chunks[t] = ix
        self.hi = t + 1

    def add_left(self, t: int, items: list[TrackItem]) -> None:
        ix = self._chunk(items)
        if self.pairs:
            ix._lanes = _pack_left_pairs(ix.items, ix.starts, ix.ends, self.last_end,
                                         self.first_start, self.placed, self.key)
        else:
            ix._lanes = _pack_left(ix.starts, ix.ends, self.last_end, self.first_start)
        self._pack_groups(ix, _pack_left_pairs if self.pairs else _pack_left)
        self.chunks[t] = ix
        self.lo = t

    def _pack_groups(self, ix: ChromIndex, pack) -> None:
        if self.group_fn is None:
            return
        labels = [self.group_fn(it) for it in ix.items]
        rows = [0] * len(labels)
        by: dict = {}
        for i, g in enumerate(labels):
            by.setdefault(g, []).append(i)
        for g, idx in by.items():
            last_end, first_start = self.groups.setdefault(g, ([], []))
            arr = np.asarray(idx)
            if self.pairs:
                got = pack([ix.items[i] for i in idx], ix.starts[arr], ix.ends[arr],
                           last_end, first_start, self.group_placed.setdefault(g, {}),
                           self.key)
            else:
                got = pack(ix.starts[arr], ix.ends[arr], last_end, first_start)
            for i, r in zip(idx, got):
                rows[i] = r
        ix._labels, ix._group_rows = labels, rows

    def query_groups(self, start: int, end: int) -> list[tuple]:
        out = []
        for ix in self._chunks_for(start, end):
            labels, rows = ix._labels, ix._group_rows
            out.extend((ix.items[i], labels[i], rows[i]) for i in ix.query(start, end))
        return out

    def complete_from(self) -> int:
        return 0 if self.lo == 0 else self.lo * self.tile + self.max_len

    def covers(self, start: int, end: int, length: int) -> bool:
        return start >= self.complete_from() and (end <= self.hi * self.tile or
                                                  self.hi * self.tile >= length)

    def _chunks_for(self, start: int, end: int):
        first = max(self.lo, (start - self.max_len) // self.tile)
        for t in range(first, min(self.hi, end // self.tile + 1)):
            yield self.chunks[t]

    def query(self, start: int, end: int) -> list[tuple[TrackItem, int]]:
        out = []
        for ix in self._chunks_for(start, end):
            rows = ix._lanes
            out.extend((ix.items[i], rows[i]) for i in ix.query(start, end))
        return out

    def _lanes_of(self, ix) -> np.ndarray:
        lanes = getattr(ix, "_lanes_np", None)
        if lanes is None or getattr(ix, "_lanes_of", None) is not ix._lanes:
            lanes = ix._lanes_np = np.asarray(ix._lanes, dtype=np.int64)
            ix._lanes_of = ix._lanes
        return lanes

    def query_below(self, start: int, end: int, rows: int) -> list[tuple[TrackItem, int]]:
        """Like ``query``, but only the reads in rows 0..rows-1 (the ones a track of
        that height shows): a deep place has far more reads than rows."""
        out = []
        for ix in self._chunks_for(start, end):
            lo, hi = ix.query_range(start, end)
            if hi <= lo:
                continue
            lanes = self._lanes_of(ix)
            sel = np.nonzero((ix.ends[lo:hi] > start) & (lanes[lo:hi] < rows))[0] + lo
            items = ix.items
            out.extend((items[i], int(lanes[i])) for i in sel.tolist())
        return out

    def max_row(self, start: int, end: int) -> int:
        """The lowest row used by the reads overlapping [start, end) (-1: none),
        without making a Python object per read."""
        top = -1
        for ix in self._chunks_for(start, end):
            lo, hi = ix.query_range(start, end)
            if hi <= lo:
                continue
            idx = np.nonzero(ix.ends[lo:hi] > start)[0]
            if not len(idx):
                continue
            top = max(top, int(self._lanes_of(ix)[idx + lo].max()))
        return top

    def count(self, start: int, end: int) -> int:
        n = 0
        for ix in self._chunks_for(start, end):
            lo, hi = ix.query_range(start, end)
            n += hi - lo
        return n

    def _chunk_depth(self, t: int) -> tuple[int, np.ndarray]:
        """(first base, depth per base) of the reads starting in tile t, built once."""
        ix = self.chunks[t]
        got = getattr(ix, "_depth", None)
        if got is None:
            s0 = t * self.tile
            e0 = max(s0 + self.tile, int(ix.ends.max()) if len(ix) else s0 + self.tile)
            got = ix._depth = (s0, ix.window_depth(s0, e0))
        return got

    def depth(self, start: int, end: int) -> np.ndarray:
        total = np.zeros(max(1, end - start), dtype=np.float32)
        for ix in self._chunks_for(start, end):
            t = ix.items[0].start // self.tile if len(ix) else None
            if t is None:
                continue
            s0, arr = self._chunk_depth(t)
            a, b = max(start, s0), min(end, s0 + arr.size)
            if b > a:
                total[a - start:b - start] += arr[a - s0:b - s0]
        return total


MIN_INTRON = 20      # without the CIGAR, a gap between blocks this long is a junction


def item_junctions(it) -> list[tuple[int, int, int]]:
    """(intron start, intron end, strand 0 + / 1 - / 2 unknown) of one alignment."""
    seg = getattr(it, "seg", None)
    out: list[tuple[int, int, int]] = []
    if seg is not None:
        cig = seg.cigartuples
        if not cig or not any(op == 3 for op, _n in cig):
            return out
        strand = 2
        for tag in ("XS", "TS"):
            if seg.has_tag(tag):
                v = seg.get_tag(tag)
                strand = 0 if v == "+" else 1 if v == "-" else 2
                if tag == "TS" and seg.is_reverse and strand < 2:
                    strand = 1 - strand      # TS is given relative to the read
                break
        ref = seg.reference_start
        for op, n in cig:
            if op == 3:
                out.append((ref, ref + n, strand))
            if op in (0, 2, 3, 7, 8):        # M D N = X move along the reference
                ref += n
        return out
    blocks = it.ref_blocks
    tags = it.info.get("_tags", {}) if isinstance(it.info, dict) else {}
    xs = str(tags.get("XS", tags.get("TS", "")))
    strand = 0 if xs.endswith("+") else 1 if xs.endswith("-") else 2
    for (s0, e0, _q0), (s1, _e1, _q1) in zip(blocks, blocks[1:]):
        if s1 - e0 >= MIN_INTRON:
            out.append((e0, s1, strand))
    return out


def _chunk_junctions(ix) -> np.ndarray:
    """Junction rows (start, end, strand) of every read in a chunk, worked out once."""
    got = getattr(ix, "_junctions", None)
    if got is None:
        rows = [j for it in ix.items for j in item_junctions(it)]
        got = np.array(rows, dtype=np.int64).reshape(-1, 3)
        ix._junctions = got
    return got


def _window_junctions(chunks, start: int, end: int) -> np.ndarray:
    """(start, end, + count, - count, unknown count) of the junctions overlapping
    [start, end), from the reads of ``chunks``."""
    parts = [_chunk_junctions(ix) for ix in chunks]
    parts = [a for a in parts if len(a)]
    if not parts:
        return np.zeros((0, 5), dtype=np.int64)
    j = np.concatenate(parts)
    j = j[(j[:, 1] > start) & (j[:, 0] < end)]
    if not len(j):
        return np.zeros((0, 5), dtype=np.int64)
    keys, inv = np.unique(j[:, :2], axis=0, return_inverse=True)
    counts = np.zeros((len(keys), 3), dtype=np.int64)
    np.add.at(counts, (inv.ravel(), j[:, 2]), 1)
    return np.concatenate([keys, counts], axis=1)


class JunctionTrack(Track):
    """The splice junctions of an alignment track (the same BAM / SAM file), drawn
    as arcs. ``source`` is that track; the view links it, or it is opened here."""

    kind = KIND_JUNCTION
    format = "junctions"

    def __init__(self, path: Path, name: str | None = None, source: Track | None = None,
                 reference_fasta: Path | None = None) -> None:
        super().__init__(path, name or f"{Path(path).name} junctions")
        self.source = source
        self._ref_fa = reference_fasta
        self.settings["view"] = "junctions"

    def _src(self) -> Track:
        if self.source is None:
            self.source = open_track(self.path, reference_fasta=self._ref_fa)
        return self.source

    def chrom_names(self) -> list[str]:
        return self._src().chrom_names()

    def resolve(self, ref_id: str) -> str | None:
        return self._src().resolve(ref_id)

    def reads_visible(self, start: int, end: int) -> bool:
        return self._src().reads_visible(start, end)

    def fetch(self, ref_id: str, start: int, end: int) -> list[tuple[TrackItem, int]]:
        return []

    def row_count(self, ref_id: str) -> int:
        return 0

    def count(self, ref_id: str) -> int:
        return 0

    def count_in(self, ref_id: str, start: int, end: int) -> int:
        return 0

    def junctions(self, ref_id: str, start: int, end: int) -> np.ndarray:
        return self._src().junctions(ref_id, start, end)

    def drawing(self):
        d = getattr(self._src(), "drawing", None)
        if d is None:
            from contextlib import nullcontext

            return nullcontext()
        return d()


class _TileItems(list):
    """A decoded tile's reads, sorted, with their aligned blocks and depth already
    worked out (in the worker thread, so the screen does not wait for it)."""

    prepared: tuple | None = None
    kept: list | None = None          # the reads the filter keeps ...
    kept_for = None                   # ... and that filter


def _prepare_tile(items: list, s0: int, tile: int, keep=None) -> list:
    out = _TileItems(items)
    try:
        shown = out
        if keep is not None:
            shown = out.kept = _TileItems(it for it in out if keep(it))
            out.kept_for = keep
        ix = ChromIndex(shown)                # sorts the list in place
        e0 = max(s0 + tile, int(ix.ends.max()) if len(ix) else s0 + tile)
        shown.prepared = (ix.block_arrays(), (s0, ix.window_depth(s0, e0)))
    except Exception:  # noqa: BLE001 - worked out on the main thread instead
        out.prepared = out.kept = out.kept_for = None
    return out


class IndexedBamTrack(Track):
    """BAM / CRAM with an index, read window by window through pysam."""

    format = "bam"
    read_window = BAM_READ_WINDOW
    TILE = 25_000               # smaller tiles: the first reads of a new place come sooner
    PAINT_WAIT = 0.015          # while drawing, wait at most this long for a tile (s)
    HTS_THREADS = 2             # htslib threads decompressing the BAM

    def __init__(self, path: Path, name: str | None = None, *,
                 reference_fasta: Path | None = None) -> None:
        super().__init__(path, name)
        import pysam

        kw = {}
        if self.path.suffix.lower() == ".cram" and reference_fasta:
            kw["reference_filename"] = str(reference_fasta)
        try:
            self._bam = pysam.AlignmentFile(str(self.path), "rc" if self.path.suffix.lower()
                                            == ".cram" else "rb", threads=self.HTS_THREADS,
                                            **kw)
        except (OSError, ValueError) as exc:
            raise TrackError(f"Cannot open {self.path.name}: {exc}") from exc
        if not self._bam.has_index():
            raise TrackError(f"{self.path.name} has no index (.bai / .csi). "
                             "Sort and index it first (samtools sort, samtools index).")
        self._lengths = dict(zip(self._bam.references, self._bam.lengths))
        # neighbouring tiles are decoded ahead in a worker (own pysam handle), so panning
        # rarely waits for the BAM
        import threading
        from concurrent.futures import ThreadPoolExecutor

        self._open_args = (str(self.path), "rc" if self.path.suffix.lower() == ".cram"
                           else "rb", kw)
        self._lock = threading.Lock()
        self._pending: dict = {}
        self._local = threading.local()
        # two workers: the tile under the view is decoded while a neighbour still is
        self._pool = ThreadPoolExecutor(max_workers=2, thread_name_prefix="bam-prefetch")
        self.on_tile_ready = None   # called (from a worker) when a tile is decoded
        self._nonblocking = 0       # > 0 while the view draws: never wait long for a tile
        self._deadline = 0.0
        self._cache = _WindowCache(size=16)
        self._region = None          # (name, start, end, ChromIndex, serial)
        self._region_serial = 0
        self._cov_cache = load_coverage_cache(self.path)
        try:
            self.read_groups = {str(g.get("ID")): {k: str(v) for k, v in g.items() if k != "ID"}
                                for g in self._bam.header.to_dict().get("RG", []) if "ID" in g}
        except (AttributeError, ValueError):
            self.read_groups = {}

    def chrom_names(self) -> list[str]:
        return list(self._lengths)

    def _insert_sample(self, ref_id: str):
        """Up to 20,000 proper pairs from the middle of the sequence."""
        key = self.resolve(ref_id)
        if key is None:
            return []
        mid = self._lengths.get(key, 0) // 2
        out = []
        with self._lock:
            for start in (mid, 0):
                for n, a in enumerate(self._bam.fetch(key, start)):
                    if a.is_proper_pair and not a.is_secondary and not a.is_supplementary:
                        out.append(BamItem(a))
                        if len(out) >= 20000:
                            return out
                    if n > 200_000:
                        break
                if len(out) >= 20:
                    break
        return out

    def count(self, ref_id: str) -> int:
        key = self.resolve(ref_id)
        if key is None:
            return 0
        try:
            for st in self._bam.get_index_statistics():
                if st.contig == key:
                    return int(st.mapped)
        except (AttributeError, ValueError):
            pass
        return -1

    def _decode_tile(self, bam, key: str, t: int) -> tuple[list[TrackItem], bool]:
        s, e = t * self.TILE, (t + 1) * self.TILE
        items: list[TrackItem] = []
        truncated = False
        for a in bam.fetch(key, s, e):
            # each read belongs to the tile it starts in
            if a.is_unmapped or a.reference_start < s:
                continue
            if len(items) >= BAM_FETCH_LIMIT:
                truncated = True
                break
            items.append(BamItem(a))
        return items, truncated

    def _tile(self, key: str, t: int, deadline: float | None = None):
        """Tile ``t``'s reads. With a ``deadline`` (while the view draws) a tile that
        is not decoded by then is left to a worker and None comes back: the view
        draws what it has and is told (``on_tile_ready``) when the tile is there."""
        with self._lock:
            hit = self._cache.get((key, t))
            fut = self._pending.get((key, t))
        if hit is not None:
            return hit
        if deadline is not None:
            import time
            from concurrent.futures import TimeoutError as _Timeout

            if fut is None:
                self._prefetch(key, t)
                with self._lock:
                    fut = self._pending.get((key, t))
                    hit = self._cache.get((key, t))
                if hit is not None:
                    return hit
            if fut is None:
                return None
            try:
                return fut.result(timeout=max(0.0, deadline - time.monotonic()))
            except _Timeout:
                return None
            except Exception:  # noqa: BLE001 - decode it here instead
                pass
        elif fut is not None:
            try:
                return fut.result()
            except Exception:  # noqa: BLE001, S110 - decode it here instead
                pass
        res = self._decode_tile(self._bam, key, t)
        with self._lock:
            self._cache.put((key, t), res)
        return res

    def drawing(self):
        """``with track.drawing():`` the view's paint never blocks on the BAM."""
        from contextlib import contextmanager

        @contextmanager
        def _ctx():
            import time

            if not self._nonblocking:
                # one wait for the whole frame, however many times it asks for reads
                self._deadline = time.monotonic() + self.PAINT_WAIT
            self._nonblocking += 1
            try:
                yield
            finally:
                self._nonblocking -= 1
        return _ctx()

    def _cancel_far(self, key: str, keep: set) -> None:
        """Drop queued (not started) decodes of tiles the view no longer needs."""
        with self._lock:
            for k, fut in list(self._pending.items()):
                if k not in keep and fut.cancel():
                    self._pending.pop(k, None)

    def _prefetch(self, key: str, t: int) -> None:
        if t < 0 or t * self.TILE >= self._lengths.get(key, 0):
            return
        with self._lock:
            if self._cache.get((key, t)) is not None or (key, t) in self._pending:
                return
            self._pending[(key, t)] = self._pool.submit(self._prefetch_job, key, t)

    def _prefetch_job(self, key: str, t: int):
        import pysam

        try:
            bam = getattr(self._local, "bam", None)
            if bam is None:
                path, mode, kw = self._open_args
                bam = self._local.bam = pysam.AlignmentFile(path, mode,
                                                            threads=self.HTS_THREADS, **kw)
            res = self._decode_tile(bam, key, t)
            res = (_prepare_tile(res[0], t * self.TILE, self.TILE, self._filter_fn), res[1])
            with self._lock:
                self._cache.put((key, t), res)
            return res
        finally:
            with self._lock:
                self._pending.pop((key, t), None)
            cb = self.on_tile_ready
            if cb is not None:
                try:
                    cb()
                except Exception:  # noqa: BLE001, S110 - the view may be gone
                    pass

    def _region_for(self, ref_id: str, start: int, end: int) -> _BamRegion | None:
        """Reads from one view width left to one right of [start, end), in tiles that
        are added as the view moves (decoded ahead in the background)."""
        key = self.resolve(ref_id)
        if key is None or not self.reads_visible(start, end):
            return None
        length = self._lengths.get(key, end)
        span = max(1, end - start)
        # what must be loaded: the view itself; half a view either side is read ahead
        # in the background (so the first reads of a new place come sooner)
        want_lo, want_hi = max(0, start), min(length, end)
        T = self.TILE
        ahead = range(max(0, (start - span // 2) // T),
                      min(length, end + span // 2) // T + 1)
        deadline = self._deadline if self._nonblocking else None
        r = self._region
        if r is not None and r.key == key and r.covers(want_lo, want_hi, length):
            for t in (*ahead, r.lo - 1, r.hi):
                self._prefetch(key, t)
            return r
        if (r is None or r.key != key or want_hi // T < r.lo - 1 or want_lo // T > r.hi
                or r.hi - r.lo >= _BamRegion.MAX_TILES):
            first = want_lo // T
            # a new place: the tiles of the view first, queued reads elsewhere dropped
            self._cancel_far(key, {(key, t) for t in ahead})
            for t in range(first, want_hi // T + 1):
                self._prefetch(key, t)
            got = self._tile(key, first, deadline)
            if got is None:
                for t in ahead:
                    self._prefetch(key, t)
                return r if r is not None and r.key == key else None
            self._region_serial += 1
            r = self._region = _BamRegion(key, self._region_serial, T, self._group_fn,
                                          self.pairs)
            r.add_right(first, self._kept(got[0]))
        self.truncated = False
        while r.hi * T < want_hi:
            got = self._tile(key, r.hi, deadline)
            if got is None:
                break
            items, cut = got
            self.truncated |= cut
            r.add_right(r.hi, self._kept(items))
        while r.lo > 0 and r.complete_from() > want_lo:
            got = self._tile(key, r.lo - 1, deadline)
            if got is None:
                break
            items, cut = got
            self.truncated |= cut
            r.add_left(r.lo - 1, self._kept(items))
        for t in (*ahead, r.lo - 1, r.hi):
            self._prefetch(key, t)
        return r

    def fetch(self, ref_id: str, start: int, end: int) -> list[tuple[TrackItem, int]]:
        r = self._region_for(ref_id, start, end)
        return r.query(start, end) if r is not None else []

    def fetch_shown(self, ref_id: str, start: int, end: int,
                    rows: int) -> list[tuple[TrackItem, int]]:
        """``fetch`` limited to the reads in the first ``rows`` rows."""
        r = self._region_for(ref_id, start, end)
        return r.query_below(start, end, rows) if r is not None else []

    def max_row_in(self, ref_id: str, start: int, end: int) -> int:
        """The lowest row of the reads overlapping [start, end) (-1: none)."""
        r = self._region_for(ref_id, start, end)
        return r.max_row(start, end) if r is not None else -1

    def row_count(self, ref_id: str) -> int:
        r = self._region
        return len(r.last_end) if r is not None else 0

    def loaded_near(self, ref_id: str, start: int, end: int) -> bool:
        """Whether any reads of [start, end) are loaded (a new place: not yet)."""
        r = self._region
        return (r is not None and r.key == self.resolve(ref_id) and bool(r.chunks)
                and r.lo * r.tile < end and r.hi * r.tile > start)

    def rows_known(self, ref_id: str) -> bool:
        """Whether reads of this sequence are loaded, so row_count means something
        (not while another sequence's reads are still in memory, or none yet)."""
        r = self._region
        return r is not None and r.key == self.resolve(ref_id)

    def count_in(self, ref_id: str, start: int, end: int) -> int:
        r = self._region_for(ref_id, start, end)
        return r.count(start, end) if r is not None else 0

    def layout_version(self, ref_id: str, start: int, end: int):
        r = self._region_for(ref_id, start, end)
        if r is None:
            return 0
        return (r.serial, tuple(sorted(r.groups))) if self.grouped else r.serial

    def _grouping_changed(self) -> None:
        self._region = None              # rows are packed again (groups, pairs, filter)

    def _kept(self, items: list) -> list:
        keep = self._filter_fn
        if keep is not None and getattr(items, "kept_for", None) is keep:
            return items.kept             # filtered by the worker already
        return items if keep is None else [it for it in items if keep(it)]

    def fetch_groups(self, ref_id: str, start: int, end: int) -> list[tuple]:
        if not self.grouped:
            return []
        r = self._region_for(ref_id, start, end)
        return r.query_groups(start, end) if r is not None else []

    def group_rows(self, ref_id: str) -> dict:
        r = self._region
        if r is None or r.key != self.resolve(ref_id) or not self.grouped:
            return {}
        return {g: len(le) for g, (le, _fs) in r.groups.items()}

    def coverage_peak(self, ref_id: str, start: int, end: int) -> float:
        """Highest depth in [start, end) among the reads already loaded (no new fetch)."""
        r = self._region
        if r is None or r.key != self.resolve(ref_id):
            return 0.0
        a, b = max(start, r.lo * r.tile), min(end, r.hi * r.tile)
        if b <= a:
            return 0.0
        # remembered per region and loaded tiles: the same view asks on every frame
        memo = r.__dict__.setdefault("_peaks", {})
        k = (a, b, r.lo, r.hi)
        if k not in memo:
            if len(memo) > 64:
                memo.clear()
            memo[k] = float(r.depth(a, b).max())
        return memo[k]

    def tile_complete(self, ref_id: str, start: int, end: int) -> bool:
        r = self._region
        key = self.resolve(ref_id)
        return r is not None and r.key == key and r.covers(start, end,
                                                          self._lengths.get(key, end))

    def coverage(self, ref_id: str, start: int, end: int, n: int) -> np.ndarray | None:
        key = self.resolve(ref_id)
        if key is None:
            return None
        if self.reads_visible(start, end):
            r = self._region_for(ref_id, start, end)
            return _max_per_point(r.depth(start, end), n) if r is not None else None
        if self._cov_cache is not None and key in self._cov_cache:
            prof = self._cov_cache[key]
            return window_profile(prof, self._cov_cache["__bin__"], start, end, n)
        return None

    def find(self, text: str, limit: int = 50) -> list[tuple[str, TrackItem]]:
        return []

    def junctions(self, ref_id: str, start: int, end: int) -> np.ndarray:
        if not self.reads_visible(start, end):
            return np.zeros((0, 5), dtype=np.int64)
        r = self._region_for(ref_id, start, end)
        if r is None:
            return np.zeros((0, 5), dtype=np.int64)
        chunks = list(r._chunks_for(start, end))
        # the junctions of a set of loaded chunks are merged once; a window of them
        # is then a cheap filter (panning asks for a new window every frame)
        key = (id(r), tuple(id(ix) for ix in chunks))
        cache = self.__dict__.setdefault("_junction_cache", {})
        hit = cache.get(key)
        # the chunks themselves are kept with it: an id is only theirs while they live
        got = hit[1] if hit is not None and all(
            a is b for a, b in zip(hit[0], chunks)) else None
        if got is None:
            lo = min((int(ix.starts.min()) for ix in chunks if len(ix)), default=start)
            hi = max((int(ix.ends.max()) for ix in chunks if len(ix)), default=end)
            got = _window_junctions(chunks, lo - 1, hi + 1)
            cache[key] = (tuple(chunks), got)
            while len(cache) > 6:
                cache.pop(next(iter(cache)))
        if not len(got):
            return got
        return got[(got[:, 1] > start) & (got[:, 0] < end)]

    def reload_coverage(self) -> None:
        self._cov_cache = load_coverage_cache(self.path)


class BamItem:
    """A BAM record shown as a track item. Only position, strand and flags are read
    when the tile is decoded; blocks, identity, bases and tags come from the record
    when first needed, so a tile of 10,000 reads decodes in a few milliseconds."""

    __slots__ = ("_layout", "end", "mapq", "secondary", "seg", "start", "strand")
    kind = "alignment"
    cs = ""
    mism = None
    thick = ()

    def __init__(self, a) -> None:
        self.seg = a
        self.start = a.reference_start
        self.end = a.reference_end or self.start + 1
        flag = a.flag
        self.strand = "-" if flag & 16 else "+"
        self.mapq = a.mapping_quality
        self.secondary = bool(flag & 0x900)
        self._layout = None

    def _lay(self):
        if self._layout is None:
            ops = [(_CIGAR_OPS[op], n) for op, n in (self.seg.cigartuples or [])]
            _end, blocks, ins, aligned = cigar_layout(self.start, ops)
            indel = sum(n for op, n in ops if op in "ID")
            self._layout = (blocks, ins, aligned, indel)
        return self._layout

    @property
    def name(self) -> str:
        return self.seg.query_name or ""

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def ref_blocks(self):
        blocks = self.seg.get_blocks()
        return [(s, e, 0) for s, e in blocks] if blocks else ((self.start, self.end, 0),)

    @property
    def blocks(self):
        return self._lay()[0]

    @property
    def ins(self):
        return self._lay()[1]

    @property
    def seq(self) -> str:
        return self.seg.query_sequence or ""

    @property
    def score(self):
        a = self.seg
        nm = a.get_tag("NM") if a.has_tag("NM") else None
        _b, _i, aligned, indel = self._lay()
        return _nm_identity(nm, aligned, indel)

    @property
    def info(self) -> dict:
        return {"flag": self.seg.flag, "cigar": self.seg.cigarstring or ""}


def _from_pysam(a) -> TrackItem:
    ops = [(_CIGAR_OPS[op], n) for op, n in (a.cigartuples or [])]
    end, blocks, ins, aligned = cigar_layout(a.reference_start, ops)
    indel = sum(n for op, n in ops if op in "ID")
    nm = a.get_tag("NM") if a.has_tag("NM") else None
    return TrackItem(a.reference_start, end, "-" if a.is_reverse else "+", a.query_name,
                     "alignment", score=_nm_identity(nm, aligned, indel),
                     mapq=a.mapping_quality, secondary=a.is_secondary or a.is_supplementary,
                     info={"flag": a.flag}, blocks=blocks, ins=ins,
                     seq=a.query_sequence or "")


class MemoryBamTrack(Track):
    """Pure-Python BAM reader (no pysam, e.g. on Windows): decodes the whole file
    into memory, so it is limited to ``BAM_MEMORY_LIMIT`` bytes."""

    format = "bam"

    def __init__(self, path: Path, name: str | None = None, *,
                 limit: int = BAM_MEMORY_LIMIT) -> None:
        super().__init__(path, name)
        size = self.path.stat().st_size
        if size > limit:
            raise TrackError(
                f"{self.path.name} is {size / 1e6:,.0f} MB. Without pysam, BAM files up to "
                f"{limit / 1e6:,.0f} MB are read into memory. Install pysam (macOS / Linux: "
                "pip install pysam) to browse large indexed BAM files, or convert a region to "
                "SAM/PAF (samtools view file.bam chr:start-end).")
        by_chrom: dict[str, list[TrackItem]] = {}
        header: list[str] = []
        for chrom, item in iter_bam(self.path, header_out=header):
            by_chrom.setdefault(chrom, []).append(item)
        if not by_chrom:
            raise TrackError(f"No mapped reads found in {self.path.name}.")
        self._chroms = {c: ChromIndex(v) for c, v in by_chrom.items()}
        self.read_groups = parse_read_groups(header)


def iter_bam(path: Path, *, keep_seq: bool = True,
             header_out: list[str] | None = None) -> Iterable[tuple[str, TrackItem]]:
    """Decode a BAM file sequentially (BGZF is plain multi-member gzip). The SAM
    header lines go to ``header_out`` when given."""
    with gzip.open(path, "rb") as fh:
        if fh.read(4) != b"BAM\x01":
            raise TrackError(f"{Path(path).name} is not a BAM file.")
        (l_text,) = struct.unpack("<i", fh.read(4))
        text = fh.read(l_text)
        if header_out is not None:
            header_out.extend(text.rstrip(b"\x00").decode(errors="replace").splitlines())
        (n_ref,) = struct.unpack("<i", fh.read(4))
        names = []
        for _ in range(n_ref):
            (l_name,) = struct.unpack("<i", fh.read(4))
            names.append(fh.read(l_name)[:-1].decode())
            fh.read(4)
        seq_codes = "=ACMGRSVTWYHKDBN"
        while True:
            head = fh.read(4)
            if len(head) < 4:
                break
            (block_size,) = struct.unpack("<i", head)
            rec = fh.read(block_size)
            (ref_id, pos, l_name, mapq, _bin, n_cig, flag, l_seq, nref, npos,
             tlen) = struct.unpack("<iiBBHHHiiii", rec[:32])
            if ref_id < 0 or flag & 4:
                continue
            o = 32
            qname = rec[o:o + l_name - 1].decode(errors="replace")
            o += l_name
            cig = struct.unpack(f"<{n_cig}I", rec[o:o + 4 * n_cig])
            o += 4 * n_cig
            ops = [(_CIGAR_OPS[c & 0xF], c >> 4) for c in cig]
            seq = ""
            if keep_seq and l_seq:
                raw = rec[o:o + (l_seq + 1) // 2]
                seq = "".join(seq_codes[b >> 4] + seq_codes[b & 0xF] for b in raw)[:l_seq]
            o += (l_seq + 1) // 2
            qual = rec[o:o + l_seq] if keep_seq and l_seq and rec[o] != 0xFF else None
            o += l_seq
            nm = _bam_int_tag(rec[o:], b"NM")
            end, blocks, ins, aligned = cigar_layout(pos, ops)
            indel = sum(n for op, n in ops if op in "ID")
            item = TrackItem(
                pos, end, "-" if flag & 16 else "+", qname, "alignment",
                score=_nm_identity(nm, aligned, indel), mapq=mapq,
                secondary=bool(flag & 0x900),
                info={"flag": flag, "_rnext": names[nref] if 0 <= nref < len(names) else "*",
                      "_pnext": npos, "_tlen": tlen, "_rawtags": rec[o:]},
                blocks=blocks, ins=ins, seq=seq)
            _put_extras(item, pos, ops, qual)
            yield names[ref_id], item


_TAG_SIZES = {b"A": 1, b"c": 1, b"C": 1, b"s": 2, b"S": 2, b"i": 4, b"I": 4, b"f": 4}
_TAG_FMT = {b"c": "<b", b"C": "<B", b"s": "<h", b"S": "<H", b"i": "<i", b"I": "<I"}


def _bam_int_tag(buf: bytes, want: bytes) -> int | None:
    o = 0
    n = len(buf)
    while o + 3 <= n:
        tag, typ = buf[o:o + 2], buf[o + 2:o + 3]
        o += 3
        if typ in _TAG_SIZES:
            size = _TAG_SIZES[typ]
            if tag == want and typ in _TAG_FMT:
                return struct.unpack(_TAG_FMT[typ], buf[o:o + size])[0]
            o += size
        elif typ in (b"Z", b"H"):
            end = buf.find(b"\x00", o)
            o = (end if end >= 0 else n) + 1
        elif typ == b"B":
            sub = buf[o:o + 1]
            (count,) = struct.unpack("<i", buf[o + 1:o + 5])
            o += 5 + count * _TAG_SIZES.get(sub, 1)
        else:
            return None
    return None


def bam_tags(buf: bytes) -> dict[str, object]:
    """All optional fields of a BAM record (the bytes after the qualities)."""
    out: dict[str, object] = {}
    o, n = 0, len(buf)
    while o + 3 <= n:
        name, typ = buf[o:o + 2].decode(errors="replace"), buf[o + 2:o + 3]
        o += 3
        if typ == b"A":
            out[name] = buf[o:o + 1].decode(errors="replace")
            o += 1
        elif typ in _TAG_FMT:
            size = _TAG_SIZES[typ]
            out[name] = struct.unpack(_TAG_FMT[typ], buf[o:o + size])[0]
            o += size
        elif typ == b"f":
            out[name] = round(struct.unpack("<f", buf[o:o + 4])[0], 6)
            o += 4
        elif typ in (b"Z", b"H"):
            end = buf.find(b"\x00", o)
            end = n if end < 0 else end
            out[name] = buf[o:end].decode(errors="replace")
            o = end + 1
        elif typ == b"B":
            sub = buf[o:o + 1]
            (count,) = struct.unpack("<i", buf[o + 1:o + 5])
            o += 5 + count * _TAG_SIZES.get(sub, 1)
            out[name] = f"array of {count}"
        else:
            break
    return out


# Coverage cache for BAM (built once, read at wide zoom)


def build_bam_coverage(path: Path, *, bin_size: int = 1000,
                       progress: Callable[[str, float], None] | None = None,
                       cancel: Callable[[], bool] | None = None) -> Path:
    """Mean depth per ``bin_size`` bin for every sequence of a BAM, saved as
    ``<bam>.rvcov.npz``. Uses pysam when present, else the pure-Python reader."""
    path = Path(path)
    profiles: dict[str, np.ndarray] = {}
    if pysam_available():
        import pysam

        with pysam.AlignmentFile(str(path), "rb") as bam:
            total = sum(bam.lengths) or 1
            done = 0
            for name, length in zip(bam.references, bam.lengths):
                n = -(-length // bin_size)
                diff = np.zeros(n * bin_size + 1, dtype=np.int32)
                for a in bam.fetch(name) if bam.has_index() else ():
                    if cancel is not None and cancel():
                        raise TrackError("Cancelled.")
                    if a.is_unmapped or a.is_secondary or a.is_supplementary:
                        continue
                    for s, e in a.get_blocks():
                        diff[s] += 1
                        diff[min(e, len(diff) - 1)] -= 1
                depth = np.cumsum(diff[:n * bin_size])
                profiles[name] = depth.reshape(n, bin_size).mean(axis=1).astype(np.float32)
                done += length
                if progress:
                    progress(name, done / total)
    else:
        by: dict[str, list[tuple[int, int]]] = {}
        for k, (chrom, it) in enumerate(iter_bam(path, keep_seq=False)):
            if cancel is not None and k % 10000 == 0 and cancel():
                raise TrackError("Cancelled.")
            if it.secondary:
                continue
            by.setdefault(chrom, []).extend((s, e) for s, e, _q in it.blocks or
                                            [(it.start, it.end, 0)])
        for chrom, spans in by.items():
            arr = np.array(spans, dtype=np.int64)
            profiles[chrom] = coverage_bins(arr[:, 0], arr[:, 1], bin_size)
    out = coverage_cache_path(path)
    np.savez_compressed(out, __bin__=np.array([bin_size]), __mtime__=np.array(
        [path.stat().st_mtime]), **profiles)
    return out


def load_coverage_cache(path: Path) -> dict | None:
    cache = coverage_cache_path(Path(path))
    if not cache.is_file():
        return None
    try:
        data = np.load(cache)
        if float(data["__mtime__"][0]) < Path(path).stat().st_mtime - 1:
            return None
        out = {k: data[k] for k in data.files if not k.startswith("__")}
        out["__bin__"] = int(data["__bin__"][0])
        return out
    except (OSError, ValueError, KeyError, zlib.error):
        return None


# --------------------------------------------------------------------------
# Signal tracks: a number along the genome (bigWig, bedGraph)
# --------------------------------------------------------------------------


class SignalTrack(Track):
    """A value per base range, drawn as bars. ``signal`` gives one value per pixel
    column: the highest in that column (NaN where the file has no data)."""

    kind = KIND_SIGNAL

    def fetch(self, ref_id: str, start: int, end: int) -> list[tuple[TrackItem, int]]:
        return []

    def row_count(self, ref_id: str) -> int:
        return 0

    def count(self, ref_id: str) -> int:
        return 0

    def count_in(self, ref_id: str, start: int, end: int) -> int:
        return 0

    def coverage(self, ref_id: str, start: int, end: int, n: int) -> np.ndarray | None:
        return self.signal(ref_id, start, end, n)

    def signal(self, ref_id: str, start: int, end: int, n: int) -> np.ndarray | None:
        raise NotImplementedError

    def value_at(self, ref_id: str, pos: int) -> float | None:
        got = self.signal(ref_id, pos, pos + 1, 1)
        if got is None or not len(got) or np.isnan(got[0]):
            return None
        return float(got[0])


class BigWigTrack(SignalTrack):
    format = "bigwig"

    def __init__(self, path: Path, name: str | None = None) -> None:
        super().__init__(path, name)
        try:
            import pyBigWig
        except ImportError as exc:
            raise TrackError("bigWig files need the pyBigWig library "
                             "(pip install pyBigWig).") from exc
        try:
            self._bw = pyBigWig.open(str(self.path))
        except (OSError, RuntimeError) as exc:
            raise TrackError(f"Cannot open {self.path.name}: {exc}") from exc
        if self._bw is None or not self._bw.isBigWig():
            raise TrackError(f"{self.path.name} is not a bigWig file.")
        self._lengths = dict(self._bw.chroms())
        self._memo: tuple | None = None

    def chrom_names(self) -> list[str]:
        return list(self._lengths)

    def signal(self, ref_id: str, start: int, end: int, n: int) -> np.ndarray | None:
        key = self.resolve(ref_id)
        if key is None or n <= 0:
            return None
        length = self._lengths[key]
        a, b = max(0, int(start)), min(length, int(end))
        if b <= a:
            return np.full(n, np.nan, dtype=np.float32)
        memo = self._memo
        if memo is not None and memo[0] == (key, a, b, n):
            return memo[1]                   # the same view asks on every frame
        try:
            if b - a <= n:
                # zoomed in past one base per column: each base's value, spread out
                vals = np.asarray(self._bw.values(key, a, b, numpy=True), dtype=np.float32)
                idx = np.minimum((np.arange(n) * (b - a)) // n, b - a - 1)
                out = vals[idx]
            else:
                # the file's zoom levels give the maximum per column directly
                out = np.array([np.nan if v is None else v for v in
                                self._bw.stats(key, a, b, nBins=n, type="max")],
                               dtype=np.float32)
        except (RuntimeError, OSError):
            return None
        if (a, b) != (start, end):
            # part of the view lies past the sequence end: no data there
            full = np.full(n, np.nan, dtype=np.float32)
            i0 = int((a - start) * n / max(1, end - start))
            i1 = max(i0 + 1, int((b - start) * n / max(1, end - start)))
            full[i0:i1] = _resample(out, i1 - i0)
            out = full
        self._memo = ((key, a, b, n), out)
        return out


def _resample(values: np.ndarray, n: int) -> np.ndarray:
    if n <= 0 or not len(values):
        return np.full(max(0, n), np.nan, dtype=np.float32)
    idx = np.minimum((np.arange(n) * len(values)) // n, len(values) - 1)
    return values[idx]


class BedGraphTrack(SignalTrack):
    """bedGraph (chrom, start, end, value), read into memory, sorted per sequence."""

    format = "bedgraph"

    def __init__(self, path: Path, name: str | None = None) -> None:
        super().__init__(path, name)
        rows: dict[str, tuple[list, list, list]] = {}
        try:
            with _open_text(self.path) as fh:
                for line in fh:
                    if not line.strip() or line.startswith(("#", "track", "browser")):
                        continue
                    f = line.split()
                    if len(f) < 4:
                        continue
                    try:
                        s0, e0, v = int(f[1]), int(f[2]), float(f[3])
                    except ValueError:
                        continue
                    r = rows.setdefault(f[0], ([], [], []))
                    r[0].append(s0)
                    r[1].append(e0)
                    r[2].append(v)
        except (OSError, EOFError) as exc:
            raise TrackError(f"Cannot read {self.path.name}: {exc}") from exc
        if not rows:
            raise TrackError(f"No values found in {self.path.name}.")
        self._data: dict[str, tuple[np.ndarray, np.ndarray, np.ndarray]] = {}
        for chrom, (st, en, va) in rows.items():
            st_a = np.asarray(st, dtype=np.int64)
            order = np.argsort(st_a, kind="stable")
            self._data[chrom] = (st_a[order], np.asarray(en, dtype=np.int64)[order],
                                 np.asarray(va, dtype=np.float32)[order])
        # ends in order too (intervals that do not overlap, the usual case): the
        # first interval of a window is found by a search
        self._ends_sorted = {c: _sorted(d[1]) for c, d in self._data.items()}

    def chrom_names(self) -> list[str]:
        return list(self._data)

    def signal(self, ref_id: str, start: int, end: int, n: int) -> np.ndarray | None:
        key = self.resolve(ref_id)
        if key is None or n <= 0:
            return None
        st, en, va = self._data[key]
        out = np.full(n, np.nan, dtype=np.float32)
        span = max(1, end - start)
        lo = int(np.searchsorted(en, start, side="right")) if self._ends_sorted[key] else 0
        hi = int(np.searchsorted(st, end, side="left"))
        if hi <= lo:
            return out
        s_, e_, v_ = st[lo:hi], en[lo:hi], va[lo:hi]
        keep = (e_ > start) & (s_ < end)
        s_, e_, v_ = s_[keep], e_[keep], v_[keep]
        if not len(v_):
            return out
        # every interval counts in the columns its two ends fall in ...
        b0 = np.clip(((s_ - start) * n) // span, 0, n - 1)
        b1 = np.clip(((e_ - 1 - start) * n) // span, 0, n - 1)
        np.fmax.at(out, b0, v_)
        np.fmax.at(out, b1, v_)
        # ... and in every column between, through the column centres it covers
        centre = start + (np.arange(n) + 0.5) * span / n
        j = np.searchsorted(s_, centre, side="right") - 1
        ok = j >= 0
        jj = np.where(ok, j, 0)
        inside = ok & (e_[jj] > centre)
        out = np.where(inside, np.fmax(out, v_[jj]), out)
        return out.astype(np.float32)


def _sorted(arr: np.ndarray) -> bool:
    return len(arr) < 2 or bool((arr[1:] >= arr[:-1]).all())


# --------------------------------------------------------------------------
# Variant tracks: VCF sites and per-sample genotypes
# --------------------------------------------------------------------------

GT_NONE, GT_NOCALL, GT_HOMREF, GT_HET, GT_HOMALT = 0, 1, 2, 3, 4
GT_NAMES = {GT_NOCALL: "no call", GT_HOMREF: "hom ref", GT_HET: "het",
            GT_HOMALT: "hom var"}
VARIANT_WINDOW = 2_000_000       # an indexed VCF is read for views up to this wide
VARIANT_FETCH_LIMIT = 200_000    # variants read per window at most


def genotype_code(gt: str) -> int:
    """VCF GT text ("0/1", "1|1", "./.", "0") -> GT_* code."""
    if not gt:
        return GT_NOCALL
    alleles = gt.replace("|", "/").split("/")
    if any(a in (".", "") for a in alleles):
        return GT_NOCALL
    if all(a == "0" for a in alleles):
        return GT_HOMREF
    if len(set(alleles)) == 1:
        return GT_HOMALT
    return GT_HET


def _variant_item(chrom: str, pos: int, vid: str, ref: str, alt: str, qual: str,
                  filt: str, info_raw: str) -> TrackItem:
    from repeatviewer.backends.gff_features import _parse_vcf_info, _vcf_variant_type

    start = pos - 1
    end = start + max(1, len(ref) if ref not in (".", "") else 1)
    info = _parse_vcf_info(info_raw)
    info["REF"], info["ALT"] = ref, alt
    if qual not in (".", ""):
        info["QUAL"] = qual
    if filt not in (".", ""):
        info["FILTER"] = filt
    name = vid if vid not in (".", "") else f"{ref}>{alt.split(',')[0]}"
    return TrackItem(start, end, ".", name, _vcf_variant_type(ref, alt), info=info)


class VariantTrack(Track):
    """VCF / BCF. Each variant is a TrackItem (ID or REF>ALT as its name, INFO,
    REF, ALT, QUAL, FILTER in ``info``); ``variants`` adds a genotype matrix
    (variants x samples, GT_* codes). A bgzipped file with a tabix / CSI index is
    read window by window (pysam); any other file is read whole."""

    kind = KIND_VARIANT
    format = "vcf"

    def __init__(self, path: Path, name: str | None = None) -> None:
        super().__init__(path, name)
        self.samples: list[str] = []
        self._gts: dict[str, np.ndarray] = {}
        self._vf = None
        self._win: tuple | None = None          # (key, lo, hi, ChromIndex, gts)
        low = self.path.name.lower()
        indexed = low.endswith((".gz", ".bgz", ".bcf")) and any(
            Path(str(self.path) + ext).is_file() for ext in (".tbi", ".csi"))
        if indexed and pysam_available():
            self._open_indexed()
        elif low.endswith(".bcf"):
            raise TrackError("BCF files need pysam and an index (bcftools index).")
        else:
            self._read_all()

    # ---- whole file

    def _read_all(self) -> None:
        rows: dict[str, list] = {}
        try:
            with _open_text(self.path) as fh:
                for line in fh:
                    if line.startswith("##"):
                        continue
                    if line.startswith("#"):
                        cols = line.rstrip("\n").split("\t")
                        self.samples = cols[9:] if len(cols) > 9 else []
                        continue
                    f = line.rstrip("\n").split("\t")
                    if len(f) < 8:
                        continue
                    try:
                        pos = int(f[1])
                    except ValueError:
                        continue
                    it = _variant_item(f[0], pos, f[2], f[3], f[4], f[5], f[6], f[7])
                    codes = b""
                    if len(f) > 9:
                        fmt = f[8].split(":")
                        gi = fmt.index("GT") if "GT" in fmt else -1
                        codes = bytes(genotype_code(c.split(":")[gi]) if gi >= 0
                                      else GT_NOCALL for c in f[9:])
                    rows.setdefault(f[0], []).append((it, codes))
        except (OSError, EOFError, UnicodeDecodeError) as exc:
            raise TrackError(f"Cannot read {self.path.name}: {exc}") from exc
        if not rows:
            raise TrackError(f"No variants found in {self.path.name}.")
        n = len(self.samples)
        for chrom, got in rows.items():
            got.sort(key=lambda r: (r[0].start, r[0].end))
            self._chroms[chrom] = ChromIndex([it for it, _c in got])
            gts = np.full((len(got), n), GT_NOCALL, dtype=np.uint8)
            for i, (_it, codes) in enumerate(got):
                if codes:
                    gts[i, :min(n, len(codes))] = np.frombuffer(codes, dtype=np.uint8)[:n]
            self._gts[chrom] = gts

    # ---- indexed

    def _open_indexed(self) -> None:
        import pysam

        try:
            self._vf = pysam.VariantFile(str(self.path))
        except (OSError, ValueError) as exc:
            raise TrackError(f"Cannot open {self.path.name}: {exc}") from exc
        self.samples = list(self._vf.header.samples)
        self.read_window = VARIANT_WINDOW
        self._names = list(self._vf.header.contigs) or list(self._vf.index.keys())

    def chrom_names(self) -> list[str]:
        return self._names if self._vf is not None else list(self._chroms)

    def _window(self, key: str, start: int, end: int):
        w = self._win
        if w is not None and w[0] == key and w[1] <= start and end <= w[2]:
            return w[3], w[4]
        span = max(1, end - start)
        lo, hi = max(0, start - span), end + span
        items: list[TrackItem] = []
        codes: list[bytes] = []
        try:
            for rec in self._vf.fetch(key, lo, hi):
                alts = ",".join(rec.alts or ())
                info = ";".join(f"{k}={v if not isinstance(v, tuple) else ','.join(map(str, v))}"
                                for k, v in rec.info.items())
                it = _variant_item(key, rec.pos, rec.id or ".", rec.ref or ".", alts or ".",
                                   "." if rec.qual is None else f"{rec.qual:g}",
                                   ";".join(rec.filter.keys()) or ".", info)
                items.append(it)
                row = []
                for smp in rec.samples.values():
                    g = smp.get("GT") if "GT" in rec.format else None
                    if not g or any(a is None for a in g):
                        row.append(GT_NOCALL)
                    elif all(a == 0 for a in g):
                        row.append(GT_HOMREF)
                    elif len(set(g)) == 1:
                        row.append(GT_HOMALT)
                    else:
                        row.append(GT_HET)
                codes.append(bytes(row))
                if len(items) >= VARIANT_FETCH_LIMIT:
                    hi = it.end            # what was read is complete up to here
                    break
        except (ValueError, OSError) as exc:
            raise TrackError(f"Cannot read {self.path.name}: {exc}") from exc
        ix = ChromIndex(items)
        n = len(self.samples)
        gts = np.full((len(items), n), GT_NOCALL, dtype=np.uint8)
        for i, c in enumerate(codes):
            if c:
                gts[i, :n] = np.frombuffer(c, dtype=np.uint8)[:n]
        self._win = (key, lo, hi, ix, gts)
        return ix, gts

    # ---- queries

    def index(self, ref_id: str) -> ChromIndex | None:
        key = self.resolve(ref_id)
        if key is None:
            return None
        if self._vf is not None:
            w = self._win
            return w[3] if w is not None and w[0] == key else None
        return self._chroms.get(key)

    def variants(self, ref_id: str, start: int, end: int):
        """(variant items overlapping [start, end), their genotypes: an array of
        len(items) x len(samples) GT_* codes)."""
        key = self.resolve(ref_id)
        empty = ([], np.zeros((0, len(self.samples)), dtype=np.uint8))
        if key is None:
            return empty
        if self._vf is not None:
            if not self.reads_visible(start, end):
                return empty
            ix, gts = self._window(key, start, end)
        else:
            ix, gts = self._chroms.get(key), self._gts.get(key)
            if ix is None:
                return empty
        idx = ix.query(start, end)
        return [ix.items[i] for i in idx], gts[idx] if len(idx) else empty[1]

    def fetch(self, ref_id: str, start: int, end: int) -> list[tuple[TrackItem, int]]:
        return [(it, 0) for it in self.variants(ref_id, start, end)[0]]

    def row_count(self, ref_id: str) -> int:
        return len(self.samples)

    def count(self, ref_id: str) -> int:
        if self._vf is not None:
            return 0                    # not known without reading the whole file
        ix = self._chroms.get(self.resolve(ref_id) or "")
        return len(ix) if ix is not None else 0

    def count_in(self, ref_id: str, start: int, end: int) -> int:
        return len(self.variants(ref_id, start, end)[0])

    def allele_freq(self, item: TrackItem, gts_row: np.ndarray | None) -> float | None:
        """Alternate allele frequency: INFO AF, or worked out from the genotypes."""
        af = item.info.get("AF")
        if af:
            try:
                return max(0.0, min(1.0, float(str(af).split(",")[0])))
            except ValueError:
                pass
        if gts_row is None or not len(gts_row):
            return None
        called = gts_row[gts_row >= GT_HOMREF]
        if not len(called):
            return None
        alt = (called == GT_HET).sum() + 2 * (called == GT_HOMALT).sum()
        return float(alt) / (2 * len(called))

    def search(self, text: str, limit: int = 200):
        """Variants whose ID or INFO contains ``text`` (whole-file VCFs only)."""
        if self._vf is not None:
            return [], 0
        q = text.strip().lower()
        if not q:
            return [], 0
        hits = [(c, it) for c, ix in self._chroms.items() for it in ix.items
                if q in it.name.lower() or any(q in str(v).lower() for v in it.info.values())]
        return hits[:limit], len(hits)


# --------------------------------------------------------------------------
# Opening files and sessions
# --------------------------------------------------------------------------


def track_kind_for(path: Path) -> str | None:
    """'annotation' / 'paf' / 'sam' / 'bam' / 'bigwig' / 'bedgraph' / 'blast', or None
    for an unknown file."""
    name = Path(path).name.lower()
    if name.endswith(VCF_SUFFIXES):
        return "vcf"
    if name.endswith(BIGWIG_SUFFIXES):
        return "bigwig"
    if name.endswith(BEDGRAPH_SUFFIXES):
        return "bedgraph"
    if name.endswith(ANNOTATION_SUFFIXES):
        return "annotation"
    if name.endswith(PAF_SUFFIXES):
        return "paf"
    if name.endswith(SAM_SUFFIXES):
        return "sam"
    if name.endswith(BAM_SUFFIXES):
        return "bam"
    if looks_like_blast_table(Path(path)):
        return "blast"
    return None


def open_track(path: Path, *, reference_names: Iterable[str] = (),
               reference_fasta: Path | None = None, settings: dict | None = None) -> Track:
    path = Path(path)
    if not path.is_file():
        raise TrackError(f"{path} does not exist.")
    kind = track_kind_for(path)
    if (settings or {}).get("view") == "junctions" and kind in ("bam", "sam"):
        track: Track = JunctionTrack(path, reference_fasta=reference_fasta)
        track.settings.update(settings or {})
        return track
    if kind == "annotation":  # read-only copy: no gene editing (GeneEditTrack)
        track = AnnotationTrack(path)
    elif kind == "paf":
        track = PafTrack(path)
    elif kind == "sam":
        track = SamTrack(path)
    elif kind == "bam":
        if pysam_available():
            track = IndexedBamTrack(path, reference_fasta=reference_fasta)
        elif path.suffix.lower() == ".cram":
            raise TrackError("CRAM files need pysam (macOS / Linux: pip install pysam).")
        else:
            track = MemoryBamTrack(path)
    elif kind == "vcf":
        track = VariantTrack(path)
    elif kind == "bigwig":
        track = BigWigTrack(path)
    elif kind == "bedgraph":
        track = BedGraphTrack(path)
    elif kind == "blast":
        track = BlastTrack(path, genome_side=(settings or {}).get("genome_side", "auto"),
                           reference_names=reference_names)
    else:
        raise TrackError(f"{path.name}: unknown track format. Use GFF/GTF/BED/VCF, PAF, SAM, "
                         "BAM, bigWig, bedGraph or BLAST tabular output (-outfmt 6).")
    if settings:
        track.settings.update(settings)
    return track


SESSION_EXTRAS = ("bookmarks", "recent", "density_track")     # kept when only the tracks are rewritten


def write_session(path: Path, reference: Path, tracks: list[Track] | list[dict],
                  locus: str = "", extras: dict | None = None) -> Path:
    """Write a session. ``extras`` (bookmarks, recent places) replace the file's own;
    those not given are kept from the file as it is."""
    data = {"format": SESSION_FORMAT, "version": 1, "reference": str(reference),
            "tracks": [t.to_json() if isinstance(t, Track) else t for t in tracks],
            "locus": locus}
    old: dict = {}
    if Path(path).is_file():
        try:
            old = json.loads(Path(path).read_text(encoding="utf-8"))
        except (OSError, ValueError):
            old = {}
    for key in SESSION_EXTRAS:
        if extras is not None and key in extras:
            data[key] = extras[key]
        elif isinstance(old, dict) and key in old:
            data[key] = old[key]
    Path(path).write_text(json.dumps(data, indent=1), encoding="utf-8")
    return Path(path)


def update_session_extras(path: Path, extras: dict) -> None:
    """Rewrite only the bookmarks / recent places of a session file."""
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return
    if not isinstance(data, dict) or data.get("format") != SESSION_FORMAT:
        return
    for key in SESSION_EXTRAS:
        if key in extras:
            data[key] = extras[key]
    Path(path).write_text(json.dumps(data, indent=1), encoding="utf-8")


def read_session(path: Path) -> dict:
    try:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        raise TrackError(f"Cannot read {Path(path).name}: {exc}") from exc
    if data.get("format") != SESSION_FORMAT:
        raise TrackError(f"{Path(path).name} is not a genome browser session.")
    base = Path(path).parent

    def resolve(p: str) -> str:
        q = Path(p)
        return str(q if q.is_absolute() else (base / q))

    data["reference"] = resolve(data.get("reference", ""))
    for t in data.get("tracks", []):
        t["path"] = resolve(t.get("path", ""))
    return data


def session_path_for(reference: Path) -> Path:
    """``<reference name>_browser.rvgb`` next to the reference (never overwrites)."""
    reference = Path(reference)
    base = reference.with_name(reference.name + "_browser")
    cand = Path(str(base) + SESSION_SUFFIX)
    k = 2
    while cand.exists():
        cand = Path(f"{base}_{k}{SESSION_SUFFIX}")
        k += 1
    return cand


def parse_locus(text: str) -> tuple[str, int, int] | None:
    """``chr5:35,618,000-35,630,000`` / ``chr5:35618000`` / ``chr5`` → (name, start1, end1)."""
    t = text.strip().replace(",", "")
    m = re.fullmatch(r"(.+?):(\d+)(?:\s*[-–.]+\s*(\d+))?", t)
    if m:
        a = int(m.group(2))
        b = int(m.group(3)) if m.group(3) else a
        return m.group(1), min(a, b), max(a, b)
    return None


# --------------------------------------------------------------------------
# Reference FASTA with random access (samtools-style .fai)
# --------------------------------------------------------------------------


@dataclass
class FaiEntry:
    name: str
    length: int
    offset: int
    line_bases: int
    line_bytes: int


def fasta_index(path: Path) -> dict[str, FaiEntry]:
    """Read ``<fasta>.fai`` when current, otherwise build it in one pass (and save it
    next to the FASTA when the folder is writable). Gzipped FASTA is not indexable."""
    path = Path(path)
    if path.name.lower().endswith(".gz"):
        raise TrackError("A gzipped reference cannot be browsed; decompress it first.")
    fai = path.with_name(path.name + ".fai")
    try:
        if fai.is_file() and fai.stat().st_mtime >= path.stat().st_mtime:
            out = {}
            for line in fai.read_text(encoding="utf-8").splitlines():
                c = line.split("\t")
                if len(c) >= 5:
                    out[c[0]] = FaiEntry(c[0], int(c[1]), int(c[2]), int(c[3]), int(c[4]))
            if out:
                return out
    except (OSError, ValueError):
        pass
    out: dict[str, FaiEntry] = {}
    cur: FaiEntry | None = None
    with path.open("rb") as fh:
        pos = 0
        for raw in fh:
            n = len(raw)
            if raw.startswith(b">"):
                words = raw[1:].split()
                name = words[0].decode(errors="replace") if words else ""
                cur = FaiEntry(name, 0, pos + n, 0, 0)
                out[name] = cur
            elif cur is not None:
                bases = len(raw.rstrip(b"\r\n"))
                if cur.line_bases == 0:
                    cur.line_bases, cur.line_bytes = bases, n
                cur.length += bases
            pos += n
    try:
        fai.write_text("".join(f"{e.name}\t{e.length}\t{e.offset}\t{e.line_bases}\t"
                               f"{e.line_bytes}\n" for e in out.values()), encoding="utf-8")
    except OSError:
        pass
    return out


def fetch_sequence(path: Path, entry: FaiEntry, start: int = 0, end: int | None = None) -> str:
    """Bases [start, end) of one record, read by offset (upper case)."""
    end = entry.length if end is None else min(end, entry.length)
    if end <= start or entry.line_bases <= 0:
        return ""
    lb, lB = entry.line_bases, entry.line_bytes
    first = entry.offset + (start // lb) * lB + start % lb
    last = entry.offset + ((end - 1) // lb) * lB + (end - 1) % lb
    with Path(path).open("rb") as fh:
        fh.seek(first)
        raw = fh.read(last - first + 1)
    return raw.replace(b"\n", b"").replace(b"\r", b"").decode("ascii", errors="replace").upper()
