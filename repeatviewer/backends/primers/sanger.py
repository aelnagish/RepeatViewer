"""Sanger verification: an .ab1 read against the sequence the clone should have.

Reads a trace file, trims the poor ends by quality (Mott's algorithm, the rule
phred and most assemblers use), aligns what is left to the expected sequence in
both directions, and lists every difference with the base-call quality behind
it. A difference with a high quality score is real; a difference with a low one
is usually the machine.

Plain text reads (.seq, .fasta) work too; they simply have no trace and no
quality, so every difference is reported without a quality.

The verdict is deliberately cautious: a clone counts as confirmed only when the
reads cover the whole expected sequence and no difference is left above the
quality cut.

Pure Python, no Qt; unit-tested in ``tests/test_primers_sanger.py``.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

from repeatviewer.backends.primers import thermo as th

TRACE_SUFFIXES = (".ab1", ".abi", ".fasta", ".fa", ".seq", ".txt")
AB1_SUFFIXES = (".ab1", ".abi")

# Mott trimming: bases are kept while the running sum of (cutoff - error) is
# positive. 0.05 is phred's default (about Q13).
TRIM_CUTOFF = 0.05
MIN_TRIMMED = 20
# A difference at or above this quality is called real, below it doubtful.
GOOD_QUALITY = 30
CHANNELS = ("G", "A", "T", "C")
_DATA_TAGS = ("DATA9", "DATA10", "DATA11", "DATA12")


@dataclass
class Trace:
    """One Sanger read: bases, quality and (for .ab1) the four dye channels."""

    name: str = ""
    path: str = ""
    bases: str = ""
    quality: list[int] = field(default_factory=list)
    locations: list[int] = field(default_factory=list)
    channels: dict[str, list[int]] = field(default_factory=dict)

    @property
    def has_trace(self) -> bool:
        return bool(self.channels and self.locations)

    def mean_quality(self, start: int = 0, end: int | None = None) -> float:
        q = self.quality[start:end if end is not None else len(self.quality)]
        return sum(q) / len(q) if q else 0.0


def read_trace(path: str | Path) -> Trace:
    """Read an .ab1 (with traces) or a plain sequence file."""
    path = Path(path)
    if path.suffix.lower() in AB1_SUFFIXES:
        return _read_ab1(path)
    return _read_plain(path)


def _read_ab1(path: Path) -> Trace:
    from Bio import SeqIO

    with path.open("rb") as fh:
        rec = SeqIO.read(fh, "abi")
    raw = rec.annotations.get("abif_raw", {}) if hasattr(rec, "annotations") else {}
    trace = Trace(name=path.stem, path=str(path), bases=str(rec.seq).upper(),
                  quality=[int(q) for q in rec.letter_annotations.get("phred_quality", [])])
    locs = raw.get("PLOC2") or raw.get("PLOC1") or []
    trace.locations = [int(x) for x in locs]
    order = raw.get("FWO_1") or b"GATC"
    if isinstance(order, bytes):
        order = order.decode("ascii", "replace")
    order = str(order).upper()
    if sorted(order) != sorted(CHANNELS):      # unreadable order: the usual one
        order = "GATC"
    for tag, base in zip(_DATA_TAGS, order):
        data = raw.get(tag)
        if data:
            trace.channels[base] = [int(v) for v in data]
    return trace


def _read_plain(path: Path) -> Trace:
    text = path.read_text(encoding="utf-8", errors="replace")
    lines = [ln.strip() for ln in text.splitlines()]
    name = path.stem
    if lines and lines[0].startswith(">"):
        name = lines[0][1:].split()[0] or name
        lines = lines[1:]
    seq = "".join(ln for ln in lines if not ln.startswith(">"))
    seq = "".join(c for c in seq.upper() if c.isalpha())
    return Trace(name=name, path=str(path), bases=seq)


def mott_trim(quality: list[int], *, cutoff: float = TRIM_CUTOFF,
              minimum: int = MIN_TRIMMED) -> tuple[int, int]:
    """Best stretch by Mott's rule; (0, len) when there is no quality."""
    n = len(quality)
    if n == 0:
        return 0, 0
    best_sum = best_start = best_end = 0
    run_sum = 0.0
    run_start = 0
    for i, q in enumerate(quality):
        score = cutoff - 10.0 ** (-max(0, int(q)) / 10.0)
        if run_sum + score <= 0:
            run_sum = 0.0
            run_start = i + 1
            continue
        run_sum += score
        if run_sum > best_sum:
            best_sum, best_start, best_end = run_sum, run_start, i + 1
    if best_end - best_start < minimum:
        return 0, n
    return best_start, best_end


@dataclass
class Difference:
    """One disagreement between the read and the expected sequence (1-based)."""

    position: int            # on the expected sequence
    expected: str
    found: str
    kind: str                # mismatch, insertion (extra base), deletion (missing base)
    quality: int = -1        # -1: no quality in the file
    read_position: int = 0   # 1-based on the aligned (trimmed, turned) read
    trace_index: int = -1    # 0-based on the read as the file holds it

    def __post_init__(self) -> None:
        # Bio.Align hands out numpy integers; keep plain ints for the UI and JSON.
        self.position = int(self.position)
        self.read_position = int(self.read_position)
        self.trace_index = int(self.trace_index)

    @property
    def doubtful(self) -> bool:
        return 0 <= self.quality < GOOD_QUALITY

    def text(self) -> str:
        q = "no quality" if self.quality < 0 else f"Q{self.quality}"
        if self.kind == "mismatch":
            what = f"{self.expected} → {self.found}"
        elif self.kind == "deletion":
            what = f"{self.expected} missing"
        else:
            what = f"extra {self.found}"
        return f"{self.position:,}: {what} ({q})"


@dataclass
class ReadCheck:
    """One read against the expected sequence."""

    name: str = ""
    path: str = ""
    strand: str = "+"
    length: int = 0
    trimmed: tuple[int, int] = (0, 0)
    start: int = 0           # 1-based span covered on the expected sequence
    end: int = 0
    identity: float = 0.0
    differences: list[Difference] = field(default_factory=list)
    trace: Trace | None = None
    note: str = ""

    @property
    def covered(self) -> int:
        return max(0, self.end - self.start + 1) if self.end else 0

    @property
    def real_differences(self) -> list[Difference]:
        return [d for d in self.differences if not d.doubtful]

    def summary(self) -> str:
        if not self.end:
            return f"{self.name}: does not match the expected sequence"
        way = "forward" if self.strand == "+" else "reverse"
        n = len(self.real_differences)
        doubt = len(self.differences) - n
        bits = [f"{self.name}: {way}, {self.start:,} to {self.end:,} "
                f"({self.identity * 100:.1f} % identity)"]
        if n:
            bits.append(f"{n} difference{'s' if n != 1 else ''}")
        if doubt:
            bits.append(f"{doubt} doubtful (low quality)")
        if not n and not doubt:
            bits.append("no differences")
        return ", ".join(bits)


def _aligner():
    from Bio import Align

    a = Align.PairwiseAligner(mode="local")
    a.match_score = 2
    a.mismatch_score = -3
    a.open_gap_score = -7
    a.extend_gap_score = -2
    return a


def _differences(expected: str, read: str, blocks, qual, offset: int, strand: str,
                 read_len: int) -> list[Difference]:
    out: list[Difference] = []
    tgt, qry = blocks

    def file_index(i: int) -> int:
        """Where an aligned-read position sits in the file's own base order."""
        j = offset + i
        if strand == "-":
            j = offset + (read_len - offset) - 1 - i if read_len else j
        return j

    def qual_at(i: int) -> int:
        j = file_index(i)
        return int(qual[j]) if qual and 0 <= j < len(qual) else -1

    for n, ((ts, te), (qs, qe)) in enumerate(zip(tgt, qry)):
        for k in range(te - ts):
            e, f = expected[ts + k], read[qs + k]
            if e != f:
                out.append(Difference(position=ts + k + 1, expected=e, found=f, kind="mismatch",
                                      quality=qual_at(qs + k), read_position=qs + k + 1,
                                      trace_index=file_index(qs + k)))
        if n + 1 < len(tgt):
            nts, nqs = tgt[n + 1][0], qry[n + 1][0]
            for k in range(nts - te):            # bases of the expected not in the read
                out.append(Difference(position=te + k + 1, expected=expected[te + k], found="-",
                                      kind="deletion", quality=qual_at(qe),
                                      read_position=qe, trace_index=file_index(qe)))
            for k in range(nqs - qe):            # bases in the read that should not be there
                out.append(Difference(position=te + 1, expected="-", found=read[qe + k],
                                      kind="insertion", quality=qual_at(qe + k),
                                      read_position=qe + k + 1,
                                      trace_index=file_index(qe + k)))
    return out


def check_read(expected: str, trace: Trace, *, cutoff: float = TRIM_CUTOFF) -> ReadCheck:
    """Trim the read, align it both ways, and list the differences."""
    expected = str(expected or "").upper()
    check = ReadCheck(name=trace.name, path=trace.path, length=len(trace.bases), trace=trace)
    if not expected or not trace.bases:
        check.note = "Nothing to compare."
        return check
    start, end = mott_trim(trace.quality, cutoff=cutoff) if trace.quality else (0, len(trace.bases))
    check.trimmed = (start, end)
    piece = trace.bases[start:end]
    if len(piece) < 10:
        check.note = "The read is too poor to use."
        return check
    aligner = _aligner()
    best = None
    for strand, seq in (("+", piece), ("-", th.reverse_complement(piece))):
        try:
            alignment = aligner.align(expected, seq)[0]
        except (ValueError, IndexError):
            continue
        if best is None or alignment.score > best[0].score:
            best = (alignment, strand, seq)
    if best is None:
        check.note = "The read does not align to the expected sequence."
        return check
    alignment, strand, seq = best
    check.strand = strand
    tgt, qry = alignment.aligned
    if len(tgt) == 0:
        check.note = "The read does not align to the expected sequence."
        return check
    check.start, check.end = int(tgt[0][0]) + 1, int(tgt[-1][1])
    check.differences = _differences(expected, seq, (tgt, qry), trace.quality, start, strand,
                                     end)
    aligned_bases = sum(int(b) - int(a) for a, b in tgt)
    wrong = sum(1 for d in check.differences if d.kind == "mismatch")
    check.identity = (aligned_bases - wrong) / aligned_bases if aligned_bases else 0.0
    if trace.quality:
        check.note = (f"Trimmed to {end - start:,} of {len(trace.bases):,} bases, "
                      f"mean Q{trace.mean_quality(start, end):.0f}.")
    else:
        check.note = f"{len(trace.bases):,} bases, no quality in this file."
    return check


def check_reads(expected: str, paths, *, cutoff: float = TRIM_CUTOFF) -> list[ReadCheck]:
    out = []
    for p in paths:
        try:
            trace = read_trace(p)
        except Exception as exc:  # noqa: BLE001 - a bad file is one bad read
            out.append(ReadCheck(name=Path(p).stem, path=str(p),
                                 note=f"Could not read this file: {exc}"))
            continue
        out.append(check_read(expected, trace, cutoff=cutoff))
    return out


def gaps_in_coverage(length: int, checks: list[ReadCheck]) -> list[tuple[int, int]]:
    """Stretches of the expected sequence no read covers (1-based, inclusive)."""
    spans = sorted((c.start, c.end) for c in checks if c.end)
    gaps, at = [], 1
    for s, e in spans:
        if s > at:
            gaps.append((at, s - 1))
        at = max(at, e + 1)
    if at <= length:
        gaps.append((at, length))
    return gaps


def verdict(expected: str, checks: list[ReadCheck]) -> tuple[str, str]:
    """(tone, sentence): green when every base is covered and nothing is wrong."""
    n = len(expected or "")
    usable = [c for c in checks if c.end]
    if not usable:
        return "amber", "No read matches the expected sequence."
    gaps = gaps_in_coverage(n, usable)
    missing = sum(b - a + 1 for a, b in gaps)
    real = sum(len(c.real_differences) for c in usable)
    doubtful = sum(len(c.differences) - len(c.real_differences) for c in usable)
    if real:
        where = "; ".join(d.text() for c in usable for d in c.real_differences[:4])
        return "amber", f"{real} difference{'s' if real != 1 else ''} from the expected " \
                        f"sequence ({where})."
    if missing:
        first = ", ".join(f"{a:,} to {b:,}" for a, b in gaps[:3])
        return "amber", (f"No difference in what the reads cover, but {missing:,} bp are not "
                         f"read ({first}).")
    extra = f" {doubtful} low-quality difference{'s' if doubtful != 1 else ''} ignored." \
        if doubtful else ""
    return "green", f"Insert confirmed: all {n:,} bp read, no difference.{extra}"


def write_report(path: str | Path, expected_name: str, expected: str,
                 checks: list[ReadCheck]) -> Path:
    """A small HTML report beside the reads."""
    import html as _html

    tone, line = verdict(expected, checks)
    colour = "#14532D" if tone == "green" else "#8a4b00"
    rows = []
    for c in checks:
        diffs = "<br>".join(_html.escape(d.text()) for d in c.differences) or "none"
        rows.append(f"<tr><td>{_html.escape(c.name)}</td><td>{_html.escape(c.summary())}</td>"
                    f"<td>{diffs}</td></tr>")
    gaps = gaps_in_coverage(len(expected), [c for c in checks if c.end])
    gap_text = ", ".join(f"{a:,} to {b:,}" for a, b in gaps) or "none"
    doc = f"""<!DOCTYPE html>
<html><head><meta charset="utf-8"><title>Sanger check: {_html.escape(expected_name)}</title>
<style>body{{font-family:-apple-system,Segoe UI,sans-serif;margin:2em;max-width:60em}}
table{{border-collapse:collapse;width:100%}}td,th{{border-bottom:1px solid #ddd;padding:6px;
text-align:left;vertical-align:top}}</style></head><body>
<h1>Sanger check: {_html.escape(expected_name)}</h1>
<p style="color:{colour}"><b>{_html.escape(line)}</b></p>
<p>Expected sequence: {len(expected):,} bp. Reads: {len(checks)}. Not read: {gap_text}.</p>
<table><tr><th>Read</th><th>Result</th><th>Differences</th></tr>{''.join(rows)}</table>
<p style="color:gray">Reads trimmed by Mott's rule (cutoff {TRIM_CUTOFF}); a difference below
Q{GOOD_QUALITY} is treated as doubtful. Made by RepeatViewer.</p>
</body></html>
"""
    out = Path(path)
    out.write_text(doc, encoding="utf-8")
    return out
