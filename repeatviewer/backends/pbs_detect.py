"""LTR-guided, strand-aware PBS annotation for RepeatViewer.

Uses *existing* terminal-repeat features; never reruns terminal-repeat detection.
Returns SequenceFeature objects with 1-based inclusive coordinates, as used by
terminal_repeat_detect.py. 
Require PBS database gathered from Neumann et al., 2019
"""
from __future__ import annotations

import csv
import logging
import re
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, TYPE_CHECKING

from Bio import SeqIO
from repeatviewer.backends.sequence_features import SequenceFeature
from repeatviewer.backends.job_cancel import (
    JobCancellation, JobCancelledError, check_cancelled, run_subprocess_with_cancel,
)

if TYPE_CHECKING:
    from repeatviewer.backends.terminal_repeat_detect import TerminalRepeatResult

MOTIF_FILENAME = "PBS_unique_motifs_Neumann2019.fasta"
OUTFMT = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore"
FIELDS = OUTFMT.split()[1:]


@dataclass(frozen=True)
class PBSSettings:
    # Window is an upper bound; max_distance constrains PBS placement.
    window: int = 40
    min_length: int = 15
    min_identity: float = 95.0  # percent
    min_coverage: float = 0.90  # fraction of reference motif aligned
    max_distance_from_ltr5_nt: int = 20
    orientation: str = "auto"  # auto, +, -
    num_threads: int = 1

    def __post_init__(self):
        if (
            self.window < 1
            or self.min_length < 1
            or self.num_threads < 1
            or self.max_distance_from_ltr5_nt < 0
        ):
            raise ValueError(
                "window, min_length and num_threads must be >= 1; "
                "max_distance_from_ltr5_nt must be >= 0"
            )
        if not 0 <= self.min_identity <= 100 or not 0 <= self.min_coverage <= 1:
            raise ValueError("min_identity must be 0..100 and min_coverage must be 0..1")
        if self.orientation not in {"auto", "+", "-"}:
            raise ValueError("orientation must be auto, + or -")


@dataclass(frozen=True)
class PBSCandidate:
    motif_id: str
    trna: str
    start: int  # on original input FASTA, 1-based inclusive
    end: int
    strand: str
    identity_percent: float
    alignment_length: int
    motif_length: int
    motif_coverage: float
    distance_from_ltr5_nt: int
    bitscore: float
    evalue: float


@dataclass(frozen=True)
class PBSResult:
    features: tuple[SequenceFeature, ...]  # best PBS only (or empty)
    best_candidate: PBSCandidate | None
    candidates: tuple[PBSCandidate, ...]
    message: str
    orientation: str | None
    # Mapping of biologically labeled LTRs; the caller may use these to
    # update strand and LTR5/LTR3 labels after auto orientation selection.
    ltr5: SequenceFeature | None
    ltr3: SequenceFeature | None


def default_motif_fasta() -> Path:
    """Prefer ``repeatviewer/database/``, then repo ``database/``."""
    package = Path(__file__).resolve().parents[1]  # repeatviewer/
    candidates = (
        package / "database" / MOTIF_FILENAME,
        package.parent / "database" / MOTIF_FILENAME,
        Path(__file__).resolve().parent / MOTIF_FILENAME,
    )
    return next((p for p in candidates if p.is_file()), candidates[0])


def motif_database_label(path: Path | None = None) -> str:
    """Short display name for the PBS motif FASTA (hover / Source field)."""
    source = Path(path) if path is not None else default_motif_fasta()
    stem = source.stem or source.name
    if stem == Path(MOTIF_FILENAME).stem or "Neumann2019" in stem:
        return "Neumann2019"
    return stem or "PBS motifs"


def pbs_motif_source(path: Path | None = None) -> str:
    """``SequenceFeature.source`` stamp for PBS motif hits (not an import)."""
    return f"PBS motif: {motif_database_label(path)}"


def pbs_motif_database_available(path: Path | None = None) -> bool:
    """True when the Neumann PBS motif FASTA is present on disk."""
    source = Path(path) if path is not None else default_motif_fasta()
    try:
        return source.expanduser().resolve().is_file()
    except OSError:
        return False


def pbs_feature_from_candidate(
    hit: PBSCandidate,
    *,
    motif_database: Path | None = None,
) -> SequenceFeature:
    """Build a PBS ``SequenceFeature`` with identity, coverage, and DB source."""
    # Keep readable feature names: avoid a percent character (GFF serializers
    # may otherwise encode it as %25).
    safe_trna = re.sub(r"[^A-Za-z0-9_.-]", "_", hit.trna or "unknown")
    return SequenceFeature(
        start=hit.start,
        end=hit.end,
        name=f"PBS_{safe_trna}",
        feature_type="PBS",
        strand=hit.strand,
        source=pbs_motif_source(motif_database),
        pident=float(hit.identity_percent),
        subject_id=str(hit.motif_id or "").strip(),
        coverage=float(hit.motif_coverage),
    )


def _extract_ltrs(terminal_repeat_result: "TerminalRepeatResult | Iterable[SequenceFeature]"):
    features = (terminal_repeat_result.features
                if hasattr(terminal_repeat_result, "features") else terminal_repeat_result)
    # Prefer the primary LTR5/LTR3 pair from Feature detection; otherwise use
    # the outermost two LTR features (works for Transfer / multi-pair maps).
    named = [
        f for f in features
        if getattr(f, "feature_type", None) == "LTR"
        and re.match(r"^LTR[53](_\d+)?$", str(getattr(f, "name", "") or ""))
    ]
    primary = [
        f for f in named
        if re.match(r"^LTR[53]$", str(getattr(f, "name", "") or ""))
    ]
    pool = primary if len(primary) >= 2 else named
    if len(pool) < 2:
        pool = [
            f for f in features if getattr(f, "feature_type", None) == "LTR"
        ]
    if len(pool) < 2:
        return None
    left_to_right = sorted(pool, key=lambda f: (int(f.start), int(f.end)))
    return left_to_right[0], left_to_right[-1]


def _ranking(hit: PBSCandidate):
    # Prefer a stronger and longer PBS match over a short full-coverage motif.
    # Do not let motif coverage alone elevate short 12-nt hits.
    return (
        -hit.identity_percent,
        -hit.alignment_length,
        -hit.motif_coverage,
        hit.distance_from_ltr5_nt,
        -hit.bitscore,
        hit.strand != "+",
        hit.motif_id,
    )


def _candidate_rows(table: Path, motif_lengths: dict[str, int],
                    strand: str, region_start: int, region_end: int,
                    anchor: int, settings: PBSSettings) -> list[PBSCandidate]:
    hits = []
    with table.open(encoding="utf-8", newline="") as handle:
        for columns in csv.reader(handle, delimiter="\t"):
            if not columns:
                continue
            if len(columns) != len(FIELDS):
                raise ValueError(f"Expected {len(FIELDS)} BLAST fields, got {len(columns)}")
            row = dict(zip(FIELDS, columns))
            motif_id = row["sseqid"]
            motif_len = motif_lengths.get(motif_id)
            if not motif_len:
                continue
            qstart, qend = int(row["qstart"]), int(row["qend"])
            sstart, send = int(row["sstart"]), int(row["send"])
            # Both DNA alignments must be plus-relative in the oriented query.
            if (qend - qstart) * (send - sstart) < 0:
                continue
            qlo, qhi = sorted((qstart, qend))
            if strand == "+":
                start, end = region_start + qlo - 1, region_start + qhi - 1
                distance = start - anchor - 1
            else:
                start, end = region_end - qhi + 1, region_end - qlo + 1
                distance = anchor - end - 1
            length = int(row["length"])
            # Require ungapped matches; BLAST alignment length otherwise includes gaps.
            if (
                int(row["gapopen"]) != 0
                or qhi - qlo + 1 != length
                or abs(send - sstart) + 1 != length
            ):
                continue
            identity = float(row["pident"])
            coverage = (abs(send - sstart) + 1) / motif_len
            if (
                length < settings.min_length
                or identity < settings.min_identity
                or coverage < settings.min_coverage
                or distance < 0
                or distance > settings.max_distance_from_ltr5_nt
                or not (region_start <= start <= end <= region_end)
            ):
                continue
            trna = motif_id.split("tRNA=", 1)[1] if "tRNA=" in motif_id else ""
            hits.append(PBSCandidate(motif_id=motif_id, trna=trna,
                start=start, end=end, strand=strand, identity_percent=identity,
                alignment_length=length, motif_length=motif_len,
                motif_coverage=round(coverage, 4), distance_from_ltr5_nt=distance,
                bitscore=float(row["bitscore"]), evalue=float(row["evalue"])))
    return hits


def _run_command(command: list[str], cancel: JobCancellation | None):
    check_cancelled(cancel)
    run_subprocess_with_cancel(command, cancel=cancel)
    check_cancelled(cancel)


def auto_annotate_pbs(
    fasta_path: Path,
    terminal_repeat_result: "TerminalRepeatResult | Iterable[SequenceFeature]",
    *,
    motif_database: Path | None = None,
    settings: PBSSettings | None = None,
    cancel: JobCancellation | None = None,
) -> PBSResult:
    """Identify the best PBS using detected terminal repeats, without re-detecting them.

    LTR left/right coordinates are assumed 1-based inclusive (terminal_repeat_detect
    passes BLAST coordinates straight to SequenceFeature). Requires two LTRs;
    absence of LTRs leads to a clean no-result return. For reverse orientation,
    the right LTR becomes biological LTR5 and the left becomes LTR3.

    Matching is a computational *candidate* only; no experimental validation.
    """
    cfg = settings or PBSSettings()
    check_cancelled(cancel)
    fasta_path = Path(fasta_path).expanduser().resolve()
    record = SeqIO.read(str(fasta_path), "fasta")
    n = len(record.seq)
    pair = _extract_ltrs(terminal_repeat_result)
    if pair is None:
        return PBSResult((), None, (), "PBS skipped: two LTR features are required.", None, None, None)
    left, right = pair
    left_end, right_start = int(left.end), int(right.start)
    if not (1 <= int(left.start) <= left_end < right_start <= int(right.end) <= n):
        raise ValueError("Invalid or overlapping LTR coordinates")
    windows = []
    if cfg.orientation in {"auto", "+"}:
        a, b = left_end + 1, min(left_end + cfg.window, right_start - 1)
        if a <= b:
            windows.append(("+", left_end, a, b))
    if cfg.orientation in {"auto", "-"}:
        a, b = max(right_start - cfg.window, left_end + 1), right_start - 1
        if a <= b:
            windows.append(("-", right_start, a, b))
    if not windows:
        return PBSResult((), None, (), "PBS skipped: no internal region between LTRs.", None, left, right)

    source = Path(motif_database) if motif_database is not None else default_motif_fasta()
    source = source.expanduser().resolve()
    if not source.is_file():
        raise FileNotFoundError(f"Local PBS motif FASTA not found: {source}")
    from repeatviewer.backends.database_prepare import _bundled_or_system
    from repeatviewer.backends.similarity_annotation import (
        load_use_system_blast_from_settings,
    )

    use_system = load_use_system_blast_from_settings()
    makeblastdb = _bundled_or_system("makeblastdb", use_system_blast=use_system)
    blastn = _bundled_or_system("blastn", use_system_blast=use_system)
    motif_lengths = {entry.id: len(entry.seq) for entry in SeqIO.parse(str(source), "fasta")}
    if not motif_lengths:
        raise ValueError(f"PBS motif FASTA is empty: {source}")

    hits: list[PBSCandidate] = []
    try:
        # Build in a temporary directory: no permission requirement for the
        # installed motif FASTA; no stale indices; no permanent side effects.
        with tempfile.TemporaryDirectory(prefix="repeatviewer_pbs_") as tmp:
            work = Path(tmp)
            db = work / "pbs_motifs"
            _run_command([makeblastdb, "-in", str(source), "-dbtype", "nucl",
                          "-out", str(db)], cancel)
            for strand, anchor, a, b in windows:
                segment = record.seq[a - 1:b]
                oriented = segment if strand == "+" else segment.reverse_complement()
                query = work / f"window_{'plus' if strand == '+' else 'minus'}.fa"
                query.write_text(f">pbs_window\n{oriented}\n", encoding="utf-8")
                output = work / f"hits_{'plus' if strand == '+' else 'minus'}.tsv"
                _run_command([blastn, "-task", "blastn", "-query", str(query),
                              "-db", str(db), "-out", str(output), "-outfmt", OUTFMT,
                              "-word_size", "7", "-dust", "no", "-evalue", "1000",
                              "-max_target_seqs", "10000", "-num_threads", str(cfg.num_threads)], cancel)
                hits.extend(_candidate_rows(output, motif_lengths, strand, a, b, anchor, cfg))
    except JobCancelledError:
        raise
    except Exception:
        logging.exception("PBS auto-annotation failed for %s", fasta_path.name)
        raise

    hits.sort(key=_ranking)
    best = hits[0] if hits else None
    if best is None:
        return PBSResult((), None, (), "No PBS motif passed the candidate filters.", None, left, right)
    feature = pbs_feature_from_candidate(best, motif_database=source)
    db_label = motif_database_label(source)
    ltr5, ltr3 = (left, right) if best.strand == "+" else (right, left)
    return PBSResult(
        features=(feature,),
        best_candidate=best,
        candidates=tuple(hits),
        message=(
            f"Detected PBS candidate ({best.start}-{best.end}, strand {best.strand}, "
            f"tRNA {best.trna or 'unknown'}, identity {best.identity_percent:.1f}%, "
            f"coverage {best.motif_coverage * 100.0:.1f}%, DB {db_label})."
        ),
        orientation=best.strand,
        ltr5=ltr5,
        ltr3=ltr3,
    )


def merge_pbs_features(existing: list[SequenceFeature],
                       detected: Iterable[SequenceFeature]) -> list[SequenceFeature]:
    """Replace prior PBS features, preserving all other annotation types."""
    return [f for f in existing if f.feature_type != "PBS"] + list(detected)
