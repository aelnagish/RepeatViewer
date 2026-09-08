"""Read alignment files (PAF, and .1aln via ALNtoPAF) as dotplots.

The dotplot canvas takes a ``DotplotResult`` full of segments, so anything that
can be turned into segments can be plotted. PAF is a natural fit: it is the
output of FastGA, minimap2 and many other aligners, so reading it means a long
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

import shutil
import subprocess
from repeatviewer.external_process import run_logged_process
from dataclasses import dataclass, field
from pathlib import Path

from repeatviewer.backends.dotplot import (
    STRAND_FORWARD,
    STRAND_REVERSE,
    DotplotResult,
    DotplotSegment,
)

PAF_SUFFIXES = frozenset({".paf"})
ALN_SUFFIXES = frozenset({".1aln"})

# A PAF of a whole-genome comparison can hold millions of records. Past this
# the plot is unreadable and the memory pointless, so stop and say so.
MAX_PAF_RECORDS = 2_000_000


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


@dataclass
class PafAlignments:
    """Parsed PAF, grouped so either layout can be produced without re-reading."""

    records: list[tuple[str, int, int, int, str, str, int, int, int]] = field(
        default_factory=list
    )
    query_lengths: dict[str, int] = field(default_factory=dict)
    target_lengths: dict[str, int] = field(default_factory=dict)
    truncated: bool = False

    def pairs(self) -> list[tuple[str, str]]:
        """Sequence pairs present, most alignments first."""
        counts: dict[tuple[str, str], int] = {}
        for qname, _ql, _qs, _qe, _st, tname, _tl, _ts, _te in self.records:
            key = (qname, tname)
            counts[key] = counts.get(key, 0) + 1
        return [k for k, _ in sorted(counts.items(), key=lambda kv: -kv[1])]

    def is_single_pair(self) -> bool:
        return len(self.pairs()) <= 1


def is_paf_file(path: Path) -> bool:
    return path.suffix.lower() in PAF_SUFFIXES


def is_aln_file(path: Path) -> bool:
    return path.suffix.lower() in ALN_SUFFIXES


def alntopaf_executable() -> str | None:
    return shutil.which("ALNtoPAF")


def convert_aln_to_paf(path: Path, out_path: Path | None = None) -> Path:
    """Convert a binary .1aln to PAF using FastGA's ALNtoPAF."""
    exe = alntopaf_executable()
    if not exe:
        raise PafError(
            ".1aln is a binary format and needs FastGA's ALNtoPAF to read.\n\n"
            "Install FastGA with:  conda install -c bioconda fastga"
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


def parse_paf(path: Path, *, max_records: int = MAX_PAF_RECORDS) -> PafAlignments:
    """Read a PAF file into grouped alignment records.

    Only the first 9 columns are needed; anything after them (CIGAR, tags) is
    ignored, so files from any aligner are accepted.
    """
    resolved = path.expanduser().resolve()
    result = PafAlignments()
    try:
        with resolved.open(encoding="utf-8", errors="replace") as handle:
            for raw in handle:
                line = raw.strip()
                if not line or line.startswith("#"):
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
                result.records.append(
                    (qname, qlen, qstart, qend, strand, tname, tlen, tstart, tend)
                )
                if qlen > result.query_lengths.get(qname, 0):
                    result.query_lengths[qname] = qlen
                if tlen > result.target_lengths.get(tname, 0):
                    result.target_lengths[tname] = tlen
                if len(result.records) >= max_records:
                    result.truncated = True
                    break
    except OSError as exc:
        raise PafError(f"Could not read {resolved.name}: {exc}") from exc

    if not result.records:
        raise PafError(f"No alignments found in {resolved.name}.")
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

    ``window`` restricts the plot to a base range (1-based inclusive) on both
    axes -- the square block on the diagonal -- which is how a small region of
    a very large scaffold is inspected without redrawing everything.
    """
    q_len = alignments.query_lengths.get(query, 0)
    t_len = alignments.target_lengths.get(target, 0)
    full_length = max(q_len, t_len, 1)

    if window is not None:
        lo = max(1, min(window[0], window[1]))
        hi = min(full_length, max(window[0], window[1]))
    else:
        lo, hi = 1, full_length
    shift = lo - 1

    segments: list[DotplotSegment] = []
    for qname, _ql, qstart, qend, strand, tname, _tl, tstart, tend in alignments.records:
        if qname != query or tname != target:
            continue
        seg = _segment(qstart, qend, strand, tstart, tend, 0, 0)
        if seg is None or seg.length < min_length:
            continue
        if window is not None:
            if not _intersects(seg, lo - 1, hi):
                continue
            seg = DotplotSegment(
                x=seg.x - shift, y=seg.y - shift, length=seg.length, strand=seg.strand
            )
        segments.append(seg)

    label = query if query == target else f"{query} vs {target}"
    return DotplotResult(
        sequence_id=label,
        sequence_length=max(hi - shift, 1),
        word_size=0,
        segments=segments,
        include_reverse=True,
        raw_match_count=len(segments),
        truncated=alignments.truncated,
        note="Truncated: too many records." if alignments.truncated else "",
        coordinate_offset=shift,
        full_sequence_length=full_length,
    )


def _axis_layout(lengths: dict[str, int]) -> tuple[list[AxisEntry], int]:
    """Lay sequences end to end, longest first, returning entries and total."""
    entries: list[AxisEntry] = []
    offset = 0
    for name, length in sorted(lengths.items(), key=lambda kv: (-kv[1], kv[0])):
        entries.append(AxisEntry(name=name, offset=offset, length=length))
        offset += length
    return entries, max(offset, 1)


def result_for_all(
    alignments: PafAlignments,
    *,
    min_length: int = 0,
) -> tuple[DotplotResult, list[AxisEntry], list[AxisEntry]]:
    """Dotplot of every sequence, laid end to end on both axes.

    Returns the result plus the axis layouts, so the viewer can draw separator
    lines and label each sequence along the edges.
    """
    x_entries, x_total = _axis_layout(alignments.query_lengths)
    y_entries, y_total = _axis_layout(alignments.target_lengths)
    x_off = {e.name: e.offset for e in x_entries}
    y_off = {e.name: e.offset for e in y_entries}

    segments: list[DotplotSegment] = []
    for qname, _ql, qstart, qend, strand, tname, _tl, tstart, tend in alignments.records:
        seg = _segment(
            qstart, qend, strand, tstart, tend,
            x_off.get(qname, 0), y_off.get(tname, 0),
        )
        if seg is not None and seg.length >= min_length:
            segments.append(seg)

    # The canvas assumes a square plot, so both axes share one span.
    axis_length = max(x_total, y_total, 1)
    return (
        DotplotResult(
            sequence_id=f"{len(x_entries)} x {len(y_entries)} sequences",
            sequence_length=axis_length,
            word_size=0,
            segments=segments,
            include_reverse=True,
            raw_match_count=len(segments),
            truncated=alignments.truncated,
            note="Truncated: too many records." if alignments.truncated else "",
        ),
        x_entries,
        y_entries,
    )


def load_alignment_file(path: Path) -> tuple[PafAlignments, Path]:
    """Read a .paf, or convert a .1aln first. Returns (alignments, paf path)."""
    resolved = path.expanduser().resolve()
    if is_aln_file(resolved):
        converted = convert_aln_to_paf(resolved)
        return parse_paf(converted), converted
    return parse_paf(resolved), resolved
