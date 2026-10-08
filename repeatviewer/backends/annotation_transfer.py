"""Alignment-based annotation transfer (Geneious-style lift via minimap2).

Source sequence owns the annotations; target receives lifted coordinates after
pairwise alignment (source = reference, target = query).
"""

from __future__ import annotations

import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from repeatviewer.backends.job_cancel import JobCancellation, JobCancelledError
from repeatviewer.backends.sequence_features import SequenceFeature


@dataclass
class AnnotationTransferOptions:
    """Thresholds for keeping a lifted feature."""

    min_identity: float = 0.5
    min_coverage: float = 0.5
    threads: int = 2
    allow_reverse: bool = True


@dataclass
class LiftedAnnotation:
    source: SequenceFeature
    status: str  # intact / changed / partial / lost
    identity: float
    feature: SequenceFeature | None = None  # None when lost / unplaced


@dataclass
class AnnotationTransferResult:
    lifted: list[LiftedAnnotation] = field(default_factory=list)
    message: str = ""
    alignment_identity: float = 0.0
    strand: str = "+"

    @property
    def features(self) -> list[SequenceFeature]:
        return [item.feature for item in self.lifted if item.feature is not None]

    def counts(self) -> dict[str, int]:
        out = {"intact": 0, "changed": 0, "partial": 0, "lost": 0, "transferred": 0}
        for item in self.lifted:
            out[item.status] = out.get(item.status, 0) + 1
            if item.feature is not None:
                out["transferred"] += 1
        return out


class AnnotationTransferError(RuntimeError):
    """User-facing transfer failure."""


def transfer_annotations(
    *,
    source_id: str,
    source_sequence: str,
    target_id: str,
    target_sequence: str,
    features: list[SequenceFeature],
    opts: AnnotationTransferOptions | None = None,
    cancel: JobCancellation | None = None,
) -> AnnotationTransferResult:
    """Align target to source and lift ``features`` onto target coordinates."""
    from repeatviewer.backends.map_reference import (
        MAP_ENGINE_MINIMAP2,
        STATUS_CHANGED,
        STATUS_INTACT,
        STATUS_LOST,
        STATUS_PARTIAL,
        engine_executable,
        lift_features,
        parse_paf_line,
    )
    from repeatviewer.conda_tools import prepare_command
    import subprocess

    cfg = opts or AnnotationTransferOptions()
    if not source_sequence or not target_sequence:
        raise AnnotationTransferError("Source and target sequences must be non-empty.")
    if not features:
        return AnnotationTransferResult(message="No features to transfer.")

    if cancel is not None and cancel.is_cancelled:
        raise JobCancelledError("Annotation transfer cancelled.")

    exe = engine_executable(MAP_ENGINE_MINIMAP2)
    if not exe:
        raise AnnotationTransferError(
            "minimap2 was not found. Install it or use the bundled copy."
        )

    with tempfile.TemporaryDirectory(prefix="rv_annot_xfer_") as tmp:
        tmp_path = Path(tmp)
        ref_fa = tmp_path / "source.fa"
        qry_fa = tmp_path / "target.fa"
        paf_path = tmp_path / "align.paf"
        SeqIO.write(
            [SeqRecord(Seq(source_sequence), id=source_id or "source", description="")],
            str(ref_fa),
            "fasta",
        )
        SeqIO.write(
            [SeqRecord(Seq(target_sequence), id=target_id or "target", description="")],
            str(qry_fa),
            "fasta",
        )
        # asm20 for TE / contig-scale; cs=long required by parse_paf_line.
        cmd = [
            exe,
            "-x",
            "asm20",
            "-t",
            str(max(1, int(cfg.threads))),
            "-c",
            "--cs=long",
            "--secondary=no",
            str(ref_fa),
            str(qry_fa),
        ]
        try:
            proc = subprocess.run(
                prepare_command(cmd),
                capture_output=True,
                text=True,
                timeout=600,
                check=False,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            raise AnnotationTransferError(f"Could not run minimap2: {exc}") from exc

        if cancel is not None and cancel.is_cancelled:
            raise JobCancelledError("Annotation transfer cancelled.")

        if proc.returncode != 0:
            detail = (proc.stderr or proc.stdout or "").strip()
            raise AnnotationTransferError(
                f"minimap2 failed (exit {proc.returncode}).\n{detail[:500]}"
            )
        paf_path.write_text(proc.stdout or "", encoding="utf-8")

        best = None
        best_score = -1.0
        for line in (proc.stdout or "").splitlines():
            if not line.strip() or line.startswith("@"):
                continue
            try:
                _qname, _qlen, read = parse_paf_line(line)
            except Exception:  # noqa: BLE001
                continue
            if read is None:
                continue
            if not cfg.allow_reverse and read.strand == "-":
                continue
            score = float(read.identity) * max(1, read.aligned_bp)
            if score > best_score:
                best_score = score
                best = read

        if best is None:
            return AnnotationTransferResult(
                message="No alignment between source and target.",
            )

        lifted_raw = lift_features(best, list(features), reference=source_sequence)
        results: list[LiftedAnnotation] = []
        src_len = len(source_sequence)
        tgt_len = len(target_sequence)

        for lf in lifted_raw:
            if not (0 <= lf.index < len(features)):
                continue
            src = features[lf.index]
            flen = max(1, int(src.end) - int(src.start) + 1)
            covered = 0
            if lf.start > 0 and lf.end >= lf.start:
                covered = lf.end - lf.start + 1
            coverage = covered / flen
            status = lf.status
            keep = (
                status in (STATUS_INTACT, STATUS_CHANGED, STATUS_PARTIAL)
                and lf.start > 0
                and lf.end >= lf.start
                and float(lf.identity) >= float(cfg.min_identity)
                and coverage >= float(cfg.min_coverage)
            )
            if status == STATUS_LOST or not keep:
                results.append(
                    LiftedAnnotation(
                        source=src,
                        status=STATUS_LOST if status == STATUS_LOST else status,
                        identity=float(lf.identity),
                        feature=None,
                    )
                )
                continue

            start = max(1, min(tgt_len, int(lf.start)))
            end = max(1, min(tgt_len, int(lf.end)))
            if end < start:
                start, end = end, start

            strand = src.strand
            if best.strand == "-":
                # Query is reverse-complemented vs reference; flip feature strand.
                strand = "-" if src.is_forward() else "+"

            results.append(
                LiftedAnnotation(
                    source=src,
                    status=status,
                    identity=float(lf.identity),
                    feature=SequenceFeature(
                        start=start,
                        end=end,
                        name=src.name,
                        feature_type=src.feature_type,
                        color=src.color,
                        strand=strand,
                        track=src.track,
                        parent=src.parent,
                        source=src.source or "Annotation transfer",
                    ),
                )
            )

        # Features that did not overlap the alignment at all → lost.
        seen_idx = {lf.index for lf in lifted_raw}
        for index, src in enumerate(features):
            if index in seen_idx:
                continue
            results.append(
                LiftedAnnotation(
                    source=src,
                    status=STATUS_LOST,
                    identity=0.0,
                    feature=None,
                )
            )

        counts = {"intact": 0, "changed": 0, "partial": 0, "lost": 0, "transferred": 0}
        for item in results:
            counts[item.status] = counts.get(item.status, 0) + 1
            if item.feature is not None:
                counts["transferred"] += 1

        msg = (
            f"Aligned {target_id} to {source_id} "
            f"(identity {best.identity:.1%}, strand {best.strand}). "
            f"Transferred {counts['transferred']} of {len(features)} feature(s) "
            f"(intact {counts.get('intact', 0)}, changed {counts.get('changed', 0)}, "
            f"partial {counts.get('partial', 0)}, lost {counts.get('lost', 0)})."
        )
        _ = src_len  # reserved for future span checks
        return AnnotationTransferResult(
            lifted=results,
            message=msg,
            alignment_identity=float(best.identity),
            strand=best.strand,
        )
