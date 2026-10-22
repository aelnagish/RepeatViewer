"""Read alignment files (PAF, and .1aln via ALNtoPAF) as dotplots.

The dotplot canvas takes a ``DotplotResult`` full of segments, so anything that
can be turned into segments can be plotted. PAF is a natural fit: it is the
output of FastGA and many other aligners, so reading it means a long
alignment can be computed once and viewed instantly afterwards, as often as
wanted.

Two layouts are supported:

* **One sequence pair** -- query against target, coordinates as in the file.
  This is what a self-comparison of a single-sequence FASTA produces.
* **All sequences** -- every query and target laid end to end on their
  respective axes, with boundaries recorded so the viewer can draw separators
  and label each sequence. This is how a multi-scaffold genome comparison is
  normally shown.

``.1aln`` is a binary ONEcode format and cannot be parsed here; it is converted
with FastGA's ``ALNtoPAF`` first.
"""

from __future__ import annotations

import heapq
import re
import subprocess
import sys
from collections.abc import Callable
from dataclasses import dataclass, field
from pathlib import Path

from repeatviewer.external_process import run_logged_process

from repeatviewer.backends.dotplot import (
    STRAND_FORWARD,
    STRAND_REVERSE,
    DotplotResult,
    DotplotSegment,
)

PAF_SUFFIXES = frozenset({".paf"})
ALN_SUFFIXES = frozenset({".1aln"})

# Keep at most this many alignment records in memory after a full-file scan.
# Unlike "first N lines", we keep the *longest* hits (FastGA ALNplot style) so
# a multi-GB PAF still shows structure instead of only early short noise.
MAX_PAF_RECORDS = 500_000

# FastGA ALNplot (-n) default: plot only the longest N alignments so short
# noise does not bury the structural diagonals on whole-genome PAFs.
# https://github.com/thegenemyers/FASTGA#ALNplot
MAX_PLOT_ALIGNMENTS = 100_000

# Built-in Dotplot writes this comment so Preview can show the run settings.
# Example: ``# RepeatViewer-dotplot word_size=15 include_reverse=1``
_PAF_RV_HEADER_PREFIX = "# RepeatViewer-dotplot"

# Concatenated "All sequences" view: past this many scaffolds the axes and
# tick labels become unusable (and building them stalls the UI).
MAX_ALL_SEQUENCES = 128

# When All-sequences is unavailable, Compare offers an overview of this many
# longest scaffolds (same layout, smaller subset).
MAX_OVERVIEW_SEQUENCES = 64

# Compare dropdown: filling a QComboBox with tens of thousands of pairs freezes
# Qt. Cap the list; Filter searches the full pair map and refills up to this.
MAX_PAIR_COMBO_ITEMS = 500


def _natural_sort_key(name: str) -> list:
    """Scaffold_1 before Scaffold_9 before Scaffold_112."""
    return [
        int(part) if part.isdigit() else part.lower()
        for part in re.split(r"(\d+)", str(name))
    ]


class PafError(RuntimeError):
    """The alignment file could not be read or converted."""


@dataclass
class AxisEntry:
    """One sequence occupying part of an axis in the concatenated layout."""

    name: str
    offset: int          # start of this sequence along the axis
    length: int

    @property
    def end(self) -> int:
        return self.offset + self.length


def format_paf_rv_header(
    *,
    word_size: int | None = None,
    include_reverse: bool | None = None,
) -> str:
    """One comment line recording built-in Dotplot run settings."""
    parts = [_PAF_RV_HEADER_PREFIX]
    if word_size is not None:
        parts.append(f"word_size={max(0, int(word_size))}")
    if include_reverse is not None:
        parts.append(f"include_reverse={1 if include_reverse else 0}")
    return " ".join(parts)


def _parse_paf_rv_header(line: str) -> dict[str, str]:
    """Pull ``key=value`` tokens from a RepeatViewer-dotplot comment line."""
    text = line.strip()
    if not text.startswith("#"):
        return {}
    body = text.lstrip("#").strip()
    if not body.lower().startswith("repeatviewer"):
        return {}
    out: dict[str, str] = {}
    for token in body.split():
        if "=" not in token:
            continue
        key, _, value = token.partition("=")
        key = key.strip().lower()
        value = value.strip()
        if key and value:
            out[key] = value
    return out


@dataclass
class PafAlignments:
    """Parsed PAF, grouped so either layout can be produced without re-reading."""

    records: list[tuple[str, int, int, int, str, str, int, int, int]] = field(
        default_factory=list
    )
    query_lengths: dict[str, int] = field(default_factory=dict)
    target_lengths: dict[str, int] = field(default_factory=dict)
    truncated: bool = False
    # Total alignment lines seen while scanning (may exceed ``len(records)``).
    total_records_seen: int = 0
    # Best (longest) span per pair across the whole file — used for Compare
    # even when only the longest records are kept in ``records``.
    pair_spans: dict[tuple[str, str], int] = field(default_factory=dict)
    # From ``# RepeatViewer-dotplot word_size=…`` when the built-in writer
    # recorded it; None for FastGA / foreign PAFs.
    word_size: int | None = None
    include_reverse: bool | None = None
    # Source path — used to re-scan one Compare pair from a huge PAF.
    source_path: Path | None = None

    def unique_pairs(self) -> set[tuple[str, str]]:
        """Observed (query, target) pairs (full-file map when available)."""
        if self.pair_spans:
            return set(self.pair_spans)
        return {
            (qname, tname)
            for qname, _ql, _qs, _qe, _st, tname, _tl, _ts, _te in self.records
        }

    def pairs(self) -> list[tuple[str, str]]:
        """Sequence pairs for the Compare menu.

        Order: self-comparisons first (natural name order), then cross pairs
        (query then target, natural-sorted).

        Only sorts pairs that actually appear in the file. An earlier version
        walked every name × name combination, which hung forever on multi-
        FastGA PAFs with hundreds of thousands of contigs.
        """
        pair_set = self.unique_pairs()
        self_pairs = sorted(
            (p for p in pair_set if p[0] == p[1]),
            key=lambda p: _natural_sort_key(p[0]),
        )
        cross_pairs = sorted(
            (p for p in pair_set if p[0] != p[1]),
            key=lambda p: (_natural_sort_key(p[0]), _natural_sort_key(p[1])),
        )
        return self_pairs + cross_pairs

    def sequence_count(self) -> int:
        """Distinct sequence names on either axis."""
        return len(set(self.query_lengths) | set(self.target_lengths))

    def supports_all_sequences_view(self) -> bool:
        """Whether the concatenated multi-scaffold layout is practical."""
        return self.sequence_count() <= MAX_ALL_SEQUENCES

    def pairs_for_combo(
        self,
        *,
        max_items: int = MAX_PAIR_COMBO_ITEMS,
        filter_text: str = "",
    ) -> tuple[list[tuple[str, str]], int]:
        """Pairs to list in the Compare menu, plus the total unique count.

        Order: self pairs first (natural name), then cross pairs (natural
        query/target). Longest span / length are tie-breakers only. Caps the
        list so the dropdown stays responsive. ``filter_text`` matches against
        the full pair map (not only the capped default list).
        """
        pair_set = self.unique_pairs()
        total = len(pair_set)
        if total == 0:
            return [], 0
        matched = filter_alignment_pairs(pair_set, filter_text)
        if not matched and filter_text.strip():
            return [], total
        lengths = _unified_sequence_lengths(self)
        spans = self.pair_spans
        ranked = sorted(
            matched,
            key=lambda pair: _pair_combo_rank(pair, lengths=lengths, spans=spans),
        )
        limit = max(1, int(max_items))
        return ranked[:limit], total

    def is_single_pair(self) -> bool:
        """True when at most one (query, target) pair appears (no full sort)."""
        if self.pair_spans:
            return len(self.pair_spans) <= 1
        seen: set[tuple[str, str]] = set()
        for qname, _ql, _qs, _qe, _st, tname, _tl, _ts, _te in self.records:
            seen.add((qname, tname))
            if len(seen) > 1:
                return False
        return True

def _pair_combo_rank(
    pair: tuple[str, str],
    *,
    lengths: dict[str, int],
    spans: dict[tuple[str, str], int],
) -> tuple:
    """Self first, then natural names; span/length only break ties."""
    query, target = pair
    q_len = lengths.get(query, 0)
    t_len = lengths.get(target, 0)
    best = int(spans.get(pair, 0)) if spans else 0
    return (
        0 if query == target else 1,
        _natural_sort_key(query),
        _natural_sort_key(target),
        -best,
        -max(q_len, t_len),
        -min(q_len, t_len),
    )


def filter_alignment_pairs(
    pairs: set[tuple[str, str]] | list[tuple[str, str]],
    filter_text: str,
) -> list[tuple[str, str]]:
    """Keep pairs whose query, target, or ``q vs t`` label contains ``filter_text``.

    Case-insensitive substring match. Empty filter returns all pairs (unsorted).
    """
    items = list(pairs)
    needle = (filter_text or "").strip().lower()
    if not needle:
        return items
    out: list[tuple[str, str]] = []
    for query, target in items:
        label = query if query == target else f"{query} vs {target}"
        if (
            needle in query.lower()
            or needle in target.lower()
            or needle in label.lower()
        ):
            out.append((query, target))
    return out


def is_paf_file(path: Path) -> bool:
    return path.suffix.lower() in PAF_SUFFIXES


def is_aln_file(path: Path) -> bool:
    return path.suffix.lower() in ALN_SUFFIXES


def alntopaf_executable() -> str | None:
    from repeatviewer.conda_tools import which_tool

    return which_tool("ALNtoPAF")


def convert_aln_to_paf(path: Path, out_path: Path | None = None) -> Path:
    """Convert a binary .1aln to PAF using FastGA's ALNtoPAF."""
    exe = alntopaf_executable()
    if not exe:
        from repeatviewer.backends.fastga_dotplot import install_hint

        raise PafError(
            ".1aln is a binary format and needs FastGA's ALNtoPAF to read.\n\n"
            + install_hint()
        )
    source = path.expanduser().resolve()
    target = out_path or source.with_suffix(".paf")
    try:
        cmd = [exe, str(source)]
        proc = run_logged_process(
            cmd,
            tool="ALNtoPAF",
            inputs=[source],
            cwd=source.parent,
            stdout_is_result=True,
            output_path=target,
        )
        target.write_text(proc.stdout, encoding="utf-8")
    except OSError as exc:
        raise PafError(f"Could not run ALNtoPAF: {exc}") from exc
    if proc.returncode != 0:
        raise PafError(
            f"ALNtoPAF failed (exit {proc.returncode}).\n{proc.stderr.strip()[:400]}"
        )
    return target


def parse_paf(
    path: Path,
    *,
    max_records: int = MAX_PAF_RECORDS,
    progress: Callable[[str], None] | None = None,
    is_cancelled: Callable[[], bool] | None = None,
) -> PafAlignments:
    """Read a PAF file into grouped alignment records.

    Only the first 9 columns are needed; anything after them (CIGAR, tags) is
    ignored, so files from any aligner are accepted.

    The whole file is scanned (with progress). Sequence lengths and pair
    identity are recorded for every line. Alignment rows kept in memory are
    the **longest** ``max_records`` hits — the FastGA ALNplot approach — not
    merely the first lines of the file (which left late contig pairs empty).

    ``progress`` receives short status strings (percent when file size is
    known). ``is_cancelled`` is polled between lines so a Stop control can
    abort a multi-hundred-megabyte read without killing the process.
    """
    resolved = path.expanduser().resolve()
    result = PafAlignments(source_path=resolved)
    keep_limit = max(1, int(max_records))
    # Min-heap of (span, seq, record) — drop the shortest when over capacity.
    heap: list[tuple[int, int, tuple]] = []
    seq = 0
    try:
        try:
            file_size = resolved.stat().st_size
        except OSError:
            file_size = 0
        # Do not call handle.tell() while iterating lines — on text streams
        # next() disables telling ("telling position disabled by next() call").
        bytes_read = 0
        next_report_at = 0
        report_step = max(file_size // 20, 1_000_000) if file_size > 0 else 0
        with resolved.open(encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                if is_cancelled is not None and is_cancelled():
                    raise PafError("Cancelled.")
                # PAF is ASCII in practice; len(raw) is a good byte proxy.
                bytes_read += len(raw)
                line = raw.strip()
                if not line:
                    continue
                if line.startswith("#"):
                    meta = _parse_paf_rv_header(line)
                    if "word_size" in meta:
                        try:
                            result.word_size = int(meta["word_size"])
                        except ValueError:
                            pass
                    if "include_reverse" in meta:
                        result.include_reverse = meta["include_reverse"] not in {
                            "0",
                            "false",
                            "no",
                        }
                    continue
                parts = line.split("\t")
                if len(parts) < 9:
                    continue
                try:
                    qname = parts[0]
                    qlen = int(parts[1])
                    qstart = int(parts[2])
                    qend = int(parts[3])
                    strand = parts[4]
                    tname = parts[5]
                    tlen = int(parts[6])
                    tstart = int(parts[7])
                    tend = int(parts[8])
                except (ValueError, IndexError):
                    continue
                span = min(max(0, qend - qstart), abs(tend - tstart))
                rec = (qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend)
                result.total_records_seen += 1
                if qlen > result.query_lengths.get(qname, 0):
                    result.query_lengths[qname] = qlen
                if tlen > result.target_lengths.get(tname, 0):
                    result.target_lengths[tname] = tlen
                pair = (qname, tname)
                if span > result.pair_spans.get(pair, -1):
                    result.pair_spans[pair] = span
                seq += 1
                item = (span, seq, rec)
                if len(heap) < keep_limit:
                    heapq.heappush(heap, item)
                elif span > heap[0][0]:
                    heapq.heapreplace(heap, item)
                if (
                    progress is not None
                    and file_size > 0
                    and bytes_read >= next_report_at
                ):
                    pct = min(99, int(100 * bytes_read / file_size))
                    progress(f"Reading alignments… {pct}%")
                    next_report_at = bytes_read + report_step
    except OSError as exc:
        raise PafError(f"Could not read {resolved.name}: {exc}") from exc

    result.records = [rec for _span, _seq, rec in heap]
    if result.total_records_seen > len(result.records):
        result.truncated = True
    if not result.records:
        raise PafError(f"No alignments found in {resolved.name}.")
    return result


def parse_paf_pair(
    path: Path,
    query: str,
    target: str,
    *,
    max_records: int = MAX_PLOT_ALIGNMENTS,
    progress: Callable[[str], None] | None = None,
    is_cancelled: Callable[[], bool] | None = None,
) -> PafAlignments:
    """Re-scan a PAF for one Compare pair (FastGA selection-style).

    Whole-genome PAFs only keep the global longest hits in memory; picking a
    contig pair re-reads the file so that pair's alignments are not missing.
    """
    resolved = path.expanduser().resolve()
    result = PafAlignments(source_path=resolved)
    keep_limit = max(1, int(max_records))
    heap: list[tuple[int, int, tuple]] = []
    seq = 0
    try:
        try:
            file_size = resolved.stat().st_size
        except OSError:
            file_size = 0
        bytes_read = 0
        next_report_at = 0
        report_step = max(file_size // 20, 1_000_000) if file_size > 0 else 0
        with resolved.open(encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                if is_cancelled is not None and is_cancelled():
                    raise PafError("Cancelled.")
                bytes_read += len(raw)
                line = raw.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t", 9)
                if len(parts) < 9:
                    continue
                if parts[0] != query or parts[5] != target:
                    continue
                try:
                    qlen = int(parts[1])
                    qstart = int(parts[2])
                    qend = int(parts[3])
                    strand = parts[4]
                    tlen = int(parts[6])
                    tstart = int(parts[7])
                    tend = int(parts[8])
                except (ValueError, IndexError):
                    continue
                span = min(max(0, qend - qstart), abs(tend - tstart))
                rec = (query, qlen, qstart, qend, strand, target, tlen, tstart, tend)
                result.total_records_seen += 1
                if qlen > result.query_lengths.get(query, 0):
                    result.query_lengths[query] = qlen
                if tlen > result.target_lengths.get(target, 0):
                    result.target_lengths[target] = tlen
                seq += 1
                item = (span, seq, rec)
                if len(heap) < keep_limit:
                    heapq.heappush(heap, item)
                elif span > heap[0][0]:
                    heapq.heapreplace(heap, item)
                if (
                    progress is not None
                    and file_size > 0
                    and bytes_read >= next_report_at
                ):
                    pct = min(99, int(100 * bytes_read / file_size))
                    progress(f"Loading {query} vs {target}… {pct}%")
                    next_report_at = bytes_read + report_step
    except OSError as exc:
        raise PafError(f"Could not read {resolved.name}: {exc}") from exc

    result.records = [rec for _span, _seq, rec in heap]
    result.pair_spans[(query, target)] = max(
        (min(max(0, r[3] - r[2]), abs(r[8] - r[7])) for r in result.records),
        default=0,
    )
    if result.total_records_seen > len(result.records):
        result.truncated = True
    if not result.records:
        raise PafError(f"No alignments for {query} vs {target} in {resolved.name}.")
    return result

def _segment(qstart: int, qend: int, strand: str, tstart: int, tend: int,
             x_off: int, y_off: int) -> DotplotSegment | None:
    """One PAF record as a dotplot segment, with axis offsets applied.

    A '-' record means the query is reverse-complemented relative to the
    target, so it runs up-right: from (qstart, tend) down to (qend, tstart).
    """
    span = min(qend - qstart, abs(tend - tstart))
    if span <= 0:
        return None
    if strand == "-":
        return DotplotSegment(
            x=qstart + x_off, y=tend + y_off, length=span, strand=STRAND_REVERSE
        )
    return DotplotSegment(
        x=qstart + x_off, y=tstart + y_off, length=span, strand=STRAND_FORWARD
    )


def _transpose_segment(seg: DotplotSegment) -> DotplotSegment:
    """Mirror a segment across y=x (fills the other triangle of the matrix).

    Pairwise PAFs store each unordered pair once (A→B, not B→A). The
    All-sequences view is symmetric, so the missing half is drawn by
    reflecting every hit.
    """
    if seg.strand == STRAND_REVERSE:
        # Reverse runs (x,y)→(x+L,y−L); mirror is (y−L,x+L)→(y,x).
        return DotplotSegment(
            x=seg.end_y,
            y=seg.end_x,
            length=seg.length,
            strand=STRAND_REVERSE,
        )
    return DotplotSegment(
        x=seg.y, y=seg.x, length=seg.length, strand=STRAND_FORWARD
    )


def keep_longest_segments(
    segments: list[DotplotSegment],
    *,
    limit: int = MAX_PLOT_ALIGNMENTS,
) -> tuple[list[DotplotSegment], int, bool]:
    """Keep the longest alignments for plotting (FastGA ``ALNplot -n``).

    Short noisy hits drown whole-genome overviews; keeping the longest
    ``limit`` records (default 100,000) matches FastGA's static viewer.
    Returned segments are ordered short→long so longer diagonals paint on top.

    Returns ``(kept, total_before_cap, was_capped)``.
    """
    import heapq

    total = len(segments)
    cap = max(0, int(limit))
    if total <= cap or cap <= 0:
        kept = sorted(segments, key=lambda s: int(s.length))
        return kept, total, False
    kept = heapq.nlargest(cap, segments, key=lambda s: int(s.length))
    kept.sort(key=lambda s: int(s.length))
    return kept, total, True


def _plot_note(*, file_truncated: bool, plot_capped: bool, kept: int, total: int) -> str:
    bits: list[str] = []
    if plot_capped:
        bits.append(
            f"Showing longest {kept:,} of {total:,} hits (FastGA ALNplot-style)."
        )
    if file_truncated:
        bits.append("File truncated: too many records.")
    return " ".join(bits)


def _intersects(seg: DotplotSegment, lo: int, hi: int) -> bool:
    """True when the segment overlaps [lo, hi] on both axes."""
    x_lo, x_hi = seg.x, seg.end_x
    y_lo, y_hi = (seg.y, seg.end_y) if seg.end_y >= seg.y else (seg.end_y, seg.y)
    return not (x_hi < lo or x_lo > hi or y_hi < lo or y_lo > hi)


def result_for_pair(
    alignments: PafAlignments,
    query: str,
    target: str,
    *,
    min_length: int = 0,
    window: tuple[int, int] | None = None,
) -> DotplotResult:
    """Dotplot of one query against one target, in file coordinates.

    X is the query length and Y is the target length. When the two sequences
    differ in size the plot is rectangular (not a square padded to
    ``max(q, t)``).

    ``window`` restricts a *self* comparison to a base range (1-based
    inclusive) on both axes — the square block on the diagonal. For unequal
    query/target pairs the shared window is ignored and the full extents are
    used.
    """
    q_len = max(int(alignments.query_lengths.get(query, 0)), 1)
    t_len = max(int(alignments.target_lengths.get(target, 0)), 1)
    is_self = query == target

    # Shared square window only makes sense when both axes are the same sequence.
    use_window = window is not None and is_self
    if use_window:
        assert window is not None
        lo = max(1, min(window[0], window[1]))
        hi = min(q_len, max(window[0], window[1]))
        shift = lo - 1
        plot_q = plot_t = max(hi - shift, 1)
    else:
        lo, hi = 1, max(q_len, t_len)
        shift = 0
        plot_q, plot_t = q_len, t_len

    segments: list[DotplotSegment] = []
    for qname, _ql, qstart, qend, strand, tname, _tl, tstart, tend in alignments.records:
        if qname != query or tname != target:
            continue
        seg = _segment(qstart, qend, strand, tstart, tend, 0, 0)
        if seg is None or seg.length < min_length:
            continue
        if use_window:
            if not _intersects(seg, lo - 1, hi):
                continue
            seg = DotplotSegment(
                x=seg.x - shift, y=seg.y - shift, length=seg.length, strand=seg.strand
            )
        segments.append(seg)

    segments, total_hits, plot_capped = keep_longest_segments(segments)
    label = query if is_self else f"{query} vs {target}"
    recorded_word = int(alignments.word_size or 0)
    return DotplotResult(
        sequence_id=label,
        sequence_length=plot_q,
        word_size=recorded_word,
        segments=segments,
        include_reverse=True,
        raw_match_count=total_hits,
        truncated=alignments.truncated or plot_capped,
        note=_plot_note(
            file_truncated=alignments.truncated,
            plot_capped=plot_capped,
            kept=len(segments),
            total=total_hits,
        ),
        coordinate_offset=shift,
        full_sequence_length=q_len,
        target_id="" if is_self else target,
        target_length=0 if is_self else plot_t,
    )


def _axis_layout(lengths: dict[str, int]) -> tuple[list[AxisEntry], int]:
    """Lay sequences end to end, longest first, returning entries and total."""
    entries: list[AxisEntry] = []
    offset = 0
    for name, length in sorted(lengths.items(), key=lambda kv: (-kv[1], kv[0])):
        entries.append(AxisEntry(name=name, offset=offset, length=length))
        offset += length
    return entries, max(offset, 1)


def _unified_sequence_lengths(alignments: PafAlignments) -> dict[str, int]:
    """Merge query and target lengths so every sequence appears on both axes.

    Pairwise PAFs (e.g. ``*_multi_dot.paf``) only list a name as query or as
    target depending on pair order. Building X from queries and Y from targets
    alone puts the same sequence at different offsets, so the identity
    diagonal breaks into staggered blocks instead of one straight line.
    """
    lengths: dict[str, int] = {}
    for name, n in alignments.query_lengths.items():
        lengths[name] = max(lengths.get(name, 0), int(n))
    for name, n in alignments.target_lengths.items():
        lengths[name] = max(lengths.get(name, 0), int(n))
    return lengths


def _result_for_lengths(
    alignments: PafAlignments,
    lengths: dict[str, int],
    *,
    min_length: int = 0,
    title: str | None = None,
) -> tuple[DotplotResult, list[AxisEntry], list[AxisEntry]]:
    """Concatenated multi-sequence plot for a given name→length map."""
    if not lengths:
        raise PafError("No sequences to plot.")
    entries, total = _axis_layout(lengths)
    x_entries = list(entries)
    y_entries = list(entries)
    offset_of = {e.name: e.offset for e in entries}
    allowed = set(offset_of)
    n_seq = len(entries)
    # Two-sequence pairwise matrices are cross-only (no A×A / B×B clutter).
    skip_self = n_seq == 2

    primary: list[DotplotSegment] = []
    for qname, _ql, qstart, qend, strand, tname, _tl, tstart, tend in alignments.records:
        if qname not in allowed or tname not in allowed:
            continue
        if skip_self and qname == tname:
            continue
        seg = _segment(
            qstart, qend, strand, tstart, tend,
            offset_of[qname], offset_of[tname],
        )
        if seg is None or seg.length < min_length:
            continue
        primary.append(seg)

    primary, total_hits, plot_capped = keep_longest_segments(primary)
    segments: list[DotplotSegment] = []
    for seg in primary:
        segments.append(seg)
        mirrored = _transpose_segment(seg)
        if (
            mirrored.x != seg.x
            or mirrored.y != seg.y
            or mirrored.strand != seg.strand
        ):
            segments.append(mirrored)

    n = len(entries)
    recorded_word = int(alignments.word_size or 0)
    return (
        DotplotResult(
            sequence_id=title or f"{n} x {n} sequences",
            sequence_length=max(total, 1),
            word_size=recorded_word,
            segments=segments,
            include_reverse=True,
            raw_match_count=total_hits,
            truncated=alignments.truncated or plot_capped,
            note=_plot_note(
                file_truncated=alignments.truncated,
                plot_capped=plot_capped,
                kept=len(primary),
                total=total_hits,
            ),
        ),
        x_entries,
        y_entries,
    )


def result_for_all(
    alignments: PafAlignments,
    *,
    min_length: int = 0,
) -> tuple[DotplotResult, list[AxisEntry], list[AxisEntry]]:
    """Dotplot of every sequence, laid end to end on both axes.

    X and Y share one sequence order and the same offsets so self-comparison
    blocks sit on a single straight identity diagonal.

    Returns the result plus the axis layouts, so the viewer can draw separator
    lines and label each sequence along the edges.

    Raises ``PafError`` when there are too many sequences for a readable
    concatenated view (see ``MAX_ALL_SEQUENCES``).
    """
    lengths = _unified_sequence_lengths(alignments)
    n_seq = len(lengths)
    if n_seq > MAX_ALL_SEQUENCES:
        raise PafError(
            f"All-sequences view supports at most {MAX_ALL_SEQUENCES} "
            f"sequences (this file has {n_seq}). Pick one pair in Compare, "
            f"or use Overview (longest {MAX_OVERVIEW_SEQUENCES})."
        )
    return _result_for_lengths(alignments, lengths, min_length=min_length)


def result_for_overview(
    alignments: PafAlignments,
    *,
    min_length: int = 0,
    max_sequences: int = MAX_OVERVIEW_SEQUENCES,
) -> tuple[DotplotResult, list[AxisEntry], list[AxisEntry]]:
    """All-sequences layout limited to the longest ``max_sequences`` scaffolds.

    Used when the full genome has too many contigs for All-sequences.
    """
    lengths = _unified_sequence_lengths(alignments)
    limit = max(1, int(max_sequences))
    if len(lengths) <= limit:
        return _result_for_lengths(alignments, lengths, min_length=min_length)
    top = dict(
        sorted(lengths.items(), key=lambda kv: (-kv[1], kv[0]))[:limit]
    )
    n = len(top)
    return _result_for_lengths(
        alignments,
        top,
        min_length=min_length,
        title=f"Overview: {n} longest of {len(lengths)} sequences",
    )


def load_alignment_file(
    path: Path,
    *,
    progress: Callable[[str], None] | None = None,
    is_cancelled: Callable[[], bool] | None = None,
) -> tuple[PafAlignments, Path]:
    """Read a .paf, or convert a .1aln first. Returns (alignments, paf path)."""
    resolved = path.expanduser().resolve()
    if is_aln_file(resolved):
        if progress is not None:
            progress(f"Converting {resolved.name}…")
        converted = convert_aln_to_paf(resolved)
        return (
            parse_paf(
                converted, progress=progress, is_cancelled=is_cancelled
            ),
            converted,
        )
    return (
        parse_paf(resolved, progress=progress, is_cancelled=is_cancelled),
        resolved,
    )
