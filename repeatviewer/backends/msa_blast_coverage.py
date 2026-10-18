"""TE-Aid-style BLAST hit coverage for MSA views.

Primary Cov path: majority consensus (threshold 0) → BLAST vs the project
genome → one coverage array aligned to MSA columns (consensus has no gaps and
the same width as the alignment).

Helpers also support per-row ungapped queries with gaps projected to 0.
"""

from __future__ import annotations

import csv
import tempfile
from collections.abc import Iterable, Sequence
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from repeatviewer.backends.msa_coordinates import GAP_CHARS, ungap_sequence


def interpolate_coverage_across_gaps(
    gapped_sequence: str,
    coverage: Sequence[int | float],
) -> list[float]:
    """Fill gap columns by connecting flanking nucleotide coverage values.

    Stored coverage keeps gaps at 0; for plotting, bridge each gap run with a
    linear ramp between the left and right non-gap depths. Leading / trailing
    gaps hold the nearest nucleotide value (or 0 if none).
    """
    seq = str(gapped_sequence or "")
    n = len(seq)
    if n <= 0:
        return []
    raw = [float(coverage[i]) if i < len(coverage) else 0.0 for i in range(n)]
    is_gap = [ch in GAP_CHARS for ch in seq]
    # Indices of non-gap columns.
    anchors = [i for i, g in enumerate(is_gap) if not g]
    if not anchors:
        return [0.0] * n
    out = list(raw)
    # Leading gaps: hold first nucleotide coverage.
    first = anchors[0]
    first_v = raw[first]
    for i in range(first):
        out[i] = first_v
    # Trailing gaps: hold last nucleotide coverage.
    last = anchors[-1]
    last_v = raw[last]
    for i in range(last + 1, n):
        out[i] = last_v
    # Interior gap runs: linear between left and right anchors.
    for a, b in zip(anchors, anchors[1:]):
        if b <= a + 1:
            continue
        left_v = raw[a]
        right_v = raw[b]
        span = float(b - a)
        for i in range(a + 1, b):
            t = (i - a) / span
            out[i] = left_v + (right_v - left_v) * t
    return out


def blast_hit_coverage(
    query_len: int,
    hits: Iterable[tuple[int, int]],
) -> list[int]:
    """Per-site hit depth for an ungapped query of length ``query_len``.

    ``hits`` are 1-based inclusive ``(qstart, qend)`` pairs (BLAST ``V7``/``V8``
    or ``qstart``/``qend``). Overlapping hits add.
    """
    length = max(0, int(query_len))
    depth = [0] * length
    if length <= 0:
        return depth
    for raw_start, raw_end in hits:
        a = int(raw_start)
        b = int(raw_end)
        if b < a:
            a, b = b, a
        # 1-based inclusive → 0-based half-open, clipped to query.
        start0 = max(0, a - 1)
        end0 = min(length, b)
        if end0 <= start0:
            continue
        for i in range(start0, end0):
            depth[i] += 1
    return depth


def project_coverage_to_gapped(
    gapped_sequence: str,
    ungapped_depth: Sequence[int],
) -> list[int]:
    """Map ungapped coverage onto MSA columns; gap columns are ``0``."""
    seq = str(gapped_sequence or "")
    if not seq:
        return []
    _ungapped, keep = ungap_sequence(seq)
    out = [0] * len(seq)
    if keep is None:
        # No gaps: 1:1 with ungapped depth (truncate/pad safely).
        for i in range(len(seq)):
            out[i] = int(ungapped_depth[i]) if i < len(ungapped_depth) else 0
        return out
    for ungapped_i, gapped_i in enumerate(keep):
        if ungapped_i < len(ungapped_depth):
            out[int(gapped_i)] = int(ungapped_depth[ungapped_i])
    return out


def parse_teaid_blast_qspans(blast_path: Path | str) -> list[tuple[int, int]]:
    """Read ``(qstart, qend)`` from a TE-Aid CSV (``V7``/``V8``) or TSV table."""
    path = Path(blast_path)
    if not path.is_file():
        return []
    text = path.read_text(encoding="utf-8", errors="replace")
    if not text.strip():
        return []
    # TE-Aid writes comma-separated V1..V13; plain blast may be tab + header.
    sample = text[:4096]
    delimiter = "," if sample.count(",") >= sample.count("\t") else "\t"
    reader = csv.reader(text.splitlines(), delimiter=delimiter)
    rows = list(reader)
    if not rows:
        return []
    header = [str(c).strip().lower() for c in rows[0]]
    start_i: int | None = None
    end_i: int | None = None
    if "v7" in header and "v8" in header:
        start_i = header.index("v7")
        end_i = header.index("v8")
        data_rows = rows[1:]
    elif "qstart" in header and "qend" in header:
        start_i = header.index("qstart")
        end_i = header.index("qend")
        data_rows = rows[1:]
    else:
        # No header: TE-Aid style 13 columns, V7/V8 at indices 6/7.
        data_rows = rows
        start_i, end_i = 6, 7
    out: list[tuple[int, int]] = []
    for row in data_rows:
        if len(row) <= max(start_i, end_i):
            continue
        try:
            out.append((int(float(row[start_i])), int(float(row[end_i]))))
        except (TypeError, ValueError):
            continue
    return out


def gapped_coverage_from_blast_table(
    gapped_sequence: str,
    blast_path: Path | str | None,
) -> list[int]:
    """Build MSA-column coverage from a BLAST/TE-Aid hit table path."""
    ungapped, _keep = ungap_sequence(gapped_sequence)
    if not ungapped:
        return [0] * len(gapped_sequence)
    if blast_path is None:
        return project_coverage_to_gapped(gapped_sequence, [0] * len(ungapped))
    spans = parse_teaid_blast_qspans(blast_path)
    depth = blast_hit_coverage(len(ungapped), spans)
    return project_coverage_to_gapped(gapped_sequence, depth)


def compute_row_gapped_coverage(
    gapped_sequence: str,
    *,
    genome_path: Path,
    output_dir: Path,
    e_value: float,
    num_threads: int = 1,
    record_id: str = "query",
    use_system_blast: bool = False,
    cancel=None,
) -> list[int]:
    """Ungap → genome BLAST → TE-Aid coverage → gapped MSA array.

    Returns zeros when BLAST finds no hits or fails in a soft way
    (``blast_n_zero``). Raises ``RuntimeError`` on hard tool failures.
    """
    from repeatviewer.backends.teaid_plotter import teaid_genome_blast

    ungapped, _keep = ungap_sequence(gapped_sequence)
    width = len(gapped_sequence)
    if not ungapped:
        return [0] * width

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    safe = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in record_id)[
        :80
    ] or "query"
    query_fa = output_dir / f"{safe}.fa"
    SeqIO.write(
        [SeqRecord(Seq(ungapped), id=safe, name=safe, description=safe)],
        str(query_fa),
        "fasta",
    )
    result = teaid_genome_blast(
        str(query_fa),
        str(genome_path),
        str(output_dir),
        e_value,
        num_threads=max(1, int(num_threads)),
        use_system_blast=bool(use_system_blast),
        output_label=safe,
        cancel=cancel,
    )
    if result == "blast_n_zero" or result is False or not result:
        return project_coverage_to_gapped(gapped_sequence, [0] * len(ungapped))
    return gapped_coverage_from_blast_table(gapped_sequence, Path(result))


def compute_msa_gapped_coverages(
    record_ids: Sequence[str],
    sequences: Sequence[str],
    *,
    genome_path: Path,
    e_value: float,
    num_threads: int = 1,
    use_system_blast: bool = False,
    cancel=None,
    progress=None,
) -> dict[str, list[int]]:
    """Run coverage for every MSA row; ``progress(completed, detail)`` optional."""
    if len(record_ids) != len(sequences):
        raise ValueError("record_ids and sequences length mismatch")
    genome = Path(genome_path)
    if not genome.is_file():
        raise FileNotFoundError(f"Genome not found: {genome}")
    out: dict[str, list[int]] = {}
    with tempfile.TemporaryDirectory(prefix="rv_msa_cov_") as tmp:
        tmp_path = Path(tmp)
        total = len(record_ids)
        for i, (rid, seq) in enumerate(zip(record_ids, sequences)):
            if cancel is not None and bool(getattr(cancel, "requested", False)):
                break
            row_dir = tmp_path / f"row_{i}"
            row_dir.mkdir(parents=True, exist_ok=True)
            out[str(rid)] = compute_row_gapped_coverage(
                seq,
                genome_path=genome,
                output_dir=row_dir,
                e_value=e_value,
                num_threads=num_threads,
                record_id=str(rid),
                use_system_blast=use_system_blast,
                cancel=cancel,
            )
            if progress is not None:
                progress(i + 1, f"{rid} ({i + 1}/{total})")
    return out


def compute_msa_consensus_blast_coverage(
    sequences: Sequence[str],
    *,
    genome_path: Path,
    e_value: float,
    num_threads: int = 1,
    use_system_blast: bool = False,
    cancel=None,
    progress=None,
) -> list[int]:
    """Build consensus (threshold 0), BLAST vs genome, return MSA-column depths.

    Consensus length equals the MSA width (no gap characters), so coverage maps
    1:1 onto alignment columns for the top-of-MSA plot.
    """
    from repeatviewer.backends.teaid_plotter import (
        _alignment_dumb_consensus,
        teaid_genome_blast,
    )

    seqs = [str(s) for s in sequences]
    if not seqs:
        return []
    width = len(seqs[0])
    if width <= 0:
        return []
    if any(len(s) != width for s in seqs):
        raise ValueError("MSA rows must share one length")
    genome = Path(genome_path)
    if not genome.is_file():
        raise FileNotFoundError(f"Genome not found: {genome}")

    if progress is not None:
        progress(1, "consensus (threshold 0)")
    if cancel is not None and bool(getattr(cancel, "requested", False)):
        return [0] * width

    consensus = _alignment_dumb_consensus(seqs, threshold=0.0, ambiguous="N")
    if not consensus:
        return [0] * width

    if progress is not None:
        progress(2, "BLAST vs genome")
    if cancel is not None and bool(getattr(cancel, "requested", False)):
        return [0] * width

    with tempfile.TemporaryDirectory(prefix="rv_msa_cons_cov_") as tmp:
        tmp_path = Path(tmp)
        query_fa = tmp_path / "consensus.fa"
        SeqIO.write(
            [
                SeqRecord(
                    Seq(consensus),
                    id="consensus",
                    name="consensus",
                    description="threshold=0",
                )
            ],
            str(query_fa),
            "fasta",
        )
        result = teaid_genome_blast(
            str(query_fa),
            str(genome),
            str(tmp_path),
            e_value,
            num_threads=max(1, int(num_threads)),
            use_system_blast=bool(use_system_blast),
            output_label="consensus",
            cancel=cancel,
        )
        if result == "blast_n_zero" or result is False or not result:
            return [0] * len(consensus)
        spans = parse_teaid_blast_qspans(Path(result))
        depth = blast_hit_coverage(len(consensus), spans)
        # Consensus has no gaps and matches MSA width — return as-is.
        if len(depth) < width:
            depth = depth + [0] * (width - len(depth))
        elif len(depth) > width:
            depth = depth[:width]
        return depth
