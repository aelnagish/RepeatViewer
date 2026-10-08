"""Hybrid annotation transfer (whole-seq lift + per-feature blastn).

Strategy (default ``hybrid``):
1. **Sequence→sequence** with minimap2, then lift features through the
   alignment (Geneious-style Transfer by Alignment / liftOver-like).
2. **Per-feature blastn** for any feature that was lost or failed thresholds —
   extract that feature’s DNA and search the target (Annotate from Database).

Per feature, the better placement (higher identity × coverage) wins when both
methods succeed. Requires a RepeatViewer ``*.seqmap_RV.json`` on the source.
"""

from __future__ import annotations

import subprocess
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

from repeatviewer.backends.job_cancel import JobCancellation, JobCancelledError
from repeatviewer.backends.sequence_features import SequenceFeature, resolve_sidecar_path

_DEFAULT_TRANSFER_SOURCE = "Annotation transfer"
# Whole-seq lift is preferred when alignment is at least this strong.
_LIFT_ALIGN_MIN_IDENTITY = 0.35
_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".ffn", ".faa"}


@dataclass(frozen=True)
class AnnotatedTransferSource:
    """One annotated sequence that can be transferred onto a target."""

    fasta_path: Path
    sequence_id: str
    label: str = ""


def _is_fasta_path(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".fa.gz") or name.endswith(".fasta.gz"):
        return True
    return path.suffix.lower() in _FASTA_SUFFIXES


def _iter_folder_fastas(folder: Path) -> list[Path]:
    """Non-recursive FASTA files in ``folder`` (sorted by name)."""
    found: list[Path] = []
    try:
        entries = sorted(folder.iterdir(), key=lambda p: p.name.lower())
    except OSError:
        return []
    for entry in entries:
        if entry.is_file() and _is_fasta_path(entry):
            found.append(entry)
    return found


def annotated_records_in_fasta(fasta_path: Path) -> list[AnnotatedTransferSource]:
    """Return records in ``fasta_path`` that already have saved annotations.

    Read-only: does not write or rename anything next to the FASTA.
    """
    path = Path(fasta_path).expanduser().resolve()
    if not path.is_file() or resolve_sidecar_path(path) is None:
        return []
    from repeatviewer.backends.sequence_features import load_features_by_id_from_json

    try:
        record_ids = [str(rec.id) for rec in SeqIO.parse(str(path), "fasta")]
    except (OSError, ValueError):
        return []
    if not record_ids:
        return []
    by_id = load_features_by_id_from_json(path, record_ids)
    out: list[AnnotatedTransferSource] = []
    for rid in record_ids:
        feats = by_id.get(rid) or []
        if not feats:
            continue
        out.append(
            AnnotatedTransferSource(
                fasta_path=path,
                sequence_id=str(rid),
                label=f"{path.stem}:{rid}" if len(record_ids) > 1 else path.stem,
            )
        )
    return out


def discover_annotated_transfer_sources(path: Path) -> list[AnnotatedTransferSource]:
    """Resolve a file or folder into annotated FASTA sequence sources.

    - File: sequences in that FASTA that have a ``*.seqmap_RV.json`` (or legacy)
      sidecar with features.
    - Folder: every FASTA in the folder that has annotations (non-recursive).
    """
    root = Path(path).expanduser().resolve()
    if root.is_dir():
        sources: list[AnnotatedTransferSource] = []
        for fasta in _iter_folder_fastas(root):
            sources.extend(annotated_records_in_fasta(fasta))
        return sources
    if root.is_file():
        return annotated_records_in_fasta(root)
    return []


@dataclass
class AnnotationTransferOptions:
    """Thresholds and strategy for annotation transfer."""

    min_identity: float = 0.5
    min_coverage: float = 0.5
    threads: int = 2
    allow_reverse: bool = True
    evalue: float = 10.0
    # hybrid | liftover | blastn
    strategy: str = "hybrid"
    annotation_source: str = _DEFAULT_TRANSFER_SOURCE
    use_system_blast: bool = False


@dataclass
class LiftedAnnotation:
    source: SequenceFeature
    status: str  # intact / changed / partial / lost
    identity: float
    feature: SequenceFeature | None = None
    method: str = ""  # liftover | blastn
    coverage: float = 0.0


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
        out = {
            "intact": 0,
            "changed": 0,
            "partial": 0,
            "lost": 0,
            "transferred": 0,
            "liftover": 0,
            "blastn": 0,
        }
        for item in self.lifted:
            out[item.status] = out.get(item.status, 0) + 1
            if item.feature is not None:
                out["transferred"] += 1
                key = str(item.method or "")
                if key in out:
                    out[key] += 1
        return out


class AnnotationTransferError(RuntimeError):
    """User-facing transfer failure."""


def require_annotation_json(source_fasta: Path) -> Path:
    """Return the RepeatViewer JSON sidecar path, or raise if missing."""
    path = Path(source_fasta).expanduser().resolve()
    sidecar = resolve_sidecar_path(path)
    if sidecar is None:
        raise AnnotationTransferError(
            "Source has no RepeatViewer annotation JSON "
            f"({path.stem}.seqmap_RV.json).\n"
            "Open the source FASTA in Annotation, add features, and Save first."
        )
    return sidecar


def _emit(progress, message: str) -> None:
    if callable(progress):
        try:
            progress(str(message))
        except Exception:  # noqa: BLE001
            pass


def _feature_dna(sequence: str, feature: SequenceFeature) -> str:
    start = max(1, int(feature.start))
    end = min(len(sequence), int(feature.end))
    if end < start or not sequence:
        return ""
    chunk = sequence[start - 1 : end]
    if feature.is_forward():
        return chunk
    return str(Seq(chunk).reverse_complement())


def _write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    lines: list[str] = []
    for seq_id, seq in records:
        lines.append(f">{seq_id}")
        for i in range(0, len(seq), 80):
            lines.append(seq[i : i + 80])
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def _parse_blast_rows(text: str) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    for line in text.splitlines():
        if not line.strip():
            continue
        fields = line.split("\t")
        if len(fields) < 10:
            continue
        rows.append(
            {
                "qseqid": fields[0],
                "sseqid": fields[1],
                "pident": fields[2],
                "length": fields[3],
                "qstart": fields[4],
                "qend": fields[5],
                "sstart": fields[6],
                "send": fields[7],
                "evalue": fields[8],
                "bitscore": fields[9],
            }
        )
    return rows


def _hit_status(
    *,
    identity: float,
    coverage: float,
    min_identity: float,
    min_coverage: float,
) -> str:
    if identity >= max(0.98, min_identity) and coverage >= max(0.98, min_coverage):
        return "intact"
    if identity >= min_identity and coverage >= min_coverage:
        if coverage < 0.95 or identity < 0.95:
            return "partial" if coverage < 0.95 else "changed"
        return "changed"
    return "lost"


def _score(item: LiftedAnnotation) -> float:
    if item.feature is None:
        return -1.0
    return float(item.identity) * max(0.0, float(item.coverage))


def _make_placed(
    *,
    src: SequenceFeature,
    start: int,
    end: int,
    strand: str,
    identity: float,
    coverage: float,
    status: str,
    method: str,
    ann_source: str,
    source_id: str,
    tgt_len: int,
) -> LiftedAnnotation:
    start = max(1, min(tgt_len, int(start)))
    end = max(1, min(tgt_len, int(end)))
    if end < start:
        start, end = end, start
    return LiftedAnnotation(
        source=src,
        status=status,
        identity=float(identity),
        coverage=float(coverage),
        method=method,
        feature=SequenceFeature(
            start=start,
            end=end,
            name=src.name,
            feature_type=src.feature_type,
            color=src.color,
            strand=strand,
            track=src.track,
            parent=src.parent,
            source=ann_source,
            pident=float(identity) * 100.0,
            subject_id=str(source_id or ""),
        ),
    )


def _transfer_via_liftover(
    *,
    source_id: str,
    source_sequence: str,
    target_id: str,
    target_sequence: str,
    features: list[SequenceFeature],
    cfg: AnnotationTransferOptions,
    cancel: JobCancellation | None,
    progress,
) -> tuple[dict[int, LiftedAnnotation], float, str]:
    """Whole-sequence minimap2 alignment → lift features. Returns by feature index."""
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

    _emit(
        progress,
        "Method A — sequence→sequence (minimap2), then lift features through "
        "the alignment…",
    )
    exe = engine_executable(MAP_ENGINE_MINIMAP2)
    if not exe:
        _emit(progress, "minimap2 not found — skipping whole-sequence lift.")
        return {}, 0.0, "+"

    if cancel is not None and cancel.requested:
        raise JobCancelledError("Annotation transfer cancelled.")

    with tempfile.TemporaryDirectory(prefix="rv_annot_xfer_lift_") as tmp:
        tmp_path = Path(tmp)
        ref_fa = tmp_path / "source.fa"
        qry_fa = tmp_path / "target.fa"
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
        _emit(progress, f"minimap2 aligning {target_id} → {source_id}…")
        try:
            proc = subprocess.run(
                prepare_command(cmd),
                capture_output=True,
                text=True,
                timeout=600,
                check=False,
            )
        except (OSError, subprocess.SubprocessError) as exc:
            _emit(progress, f"minimap2 failed ({exc}) — skipping lift.")
            return {}, 0.0, "+"

        if cancel is not None and cancel.requested:
            raise JobCancelledError("Annotation transfer cancelled.")
        if proc.returncode != 0:
            detail = (proc.stderr or proc.stdout or "").strip()[:200]
            _emit(progress, f"minimap2 exit {proc.returncode}: {detail}")
            return {}, 0.0, "+"

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
            _emit(progress, "No whole-sequence alignment — skipping lift.")
            return {}, 0.0, "+"

        align_id = float(best.identity)
        _emit(
            progress,
            f"Whole-seq alignment identity {align_id:.1%}, strand {best.strand}. "
            "Lifting features…",
        )
        if align_id < _LIFT_ALIGN_MIN_IDENTITY:
            _emit(
                progress,
                f"Alignment identity {align_id:.1%} is low for reliable lift "
                f"(<{_LIFT_ALIGN_MIN_IDENTITY:.0%}) — will rely more on blastn.",
            )

        lifted_raw = lift_features(best, list(features), reference=source_sequence)
        tgt_len = len(target_sequence)
        ann_source = str(cfg.annotation_source or _DEFAULT_TRANSFER_SOURCE).strip()
        by_index: dict[int, LiftedAnnotation] = {}
        seen: set[int] = set()

        for lf in lifted_raw:
            if not (0 <= lf.index < len(features)):
                continue
            seen.add(lf.index)
            src = features[lf.index]
            flen = max(1, int(src.end) - int(src.start) + 1)
            covered = 0
            if lf.start > 0 and lf.end >= lf.start:
                covered = lf.end - lf.start + 1
            coverage = covered / float(flen)
            identity = float(lf.identity)
            status = lf.status
            keep = (
                status in (STATUS_INTACT, STATUS_CHANGED, STATUS_PARTIAL)
                and lf.start > 0
                and lf.end >= lf.start
                and identity >= float(cfg.min_identity)
                and coverage >= float(cfg.min_coverage)
            )
            if status == STATUS_LOST or not keep:
                by_index[lf.index] = LiftedAnnotation(
                    source=src,
                    status=STATUS_LOST if status == STATUS_LOST else status,
                    identity=identity,
                    coverage=coverage,
                    method="liftover",
                    feature=None,
                )
                continue
            strand = src.strand
            if best.strand == "-":
                strand = "-" if src.is_forward() else "+"
            by_index[lf.index] = _make_placed(
                src=src,
                start=lf.start,
                end=lf.end,
                strand=strand,
                identity=identity,
                coverage=coverage,
                status=status if status in (STATUS_INTACT, STATUS_CHANGED, STATUS_PARTIAL) else "changed",
                method="liftover",
                ann_source=ann_source,
                source_id=source_id,
                tgt_len=tgt_len,
            )

        for index, src in enumerate(features):
            if index in seen:
                continue
            by_index[index] = LiftedAnnotation(
                source=src,
                status=STATUS_LOST,
                identity=0.0,
                coverage=0.0,
                method="liftover",
                feature=None,
            )

        placed = sum(1 for item in by_index.values() if item.feature is not None)
        _emit(
            progress,
            f"Lift placed {placed}/{len(features)} feature(s) "
            f"(alignment identity {align_id:.1%}).",
        )
        return by_index, align_id, best.strand


def _transfer_via_blastn(
    *,
    source_id: str,
    source_sequence: str,
    target_id: str,
    target_sequence: str,
    features: list[SequenceFeature],
    feature_indices: list[int] | None,
    cfg: AnnotationTransferOptions,
    cancel: JobCancellation | None,
    progress,
) -> dict[int, LiftedAnnotation]:
    """Per-feature DNA blastn against the target. Optional subset of indices."""
    from repeatviewer.backends.similarity_annotation import _blast_executable
    from repeatviewer.external_process import run_logged_process

    indices = (
        list(feature_indices)
        if feature_indices is not None
        else list(range(len(features)))
    )
    if not indices:
        return {}

    _emit(
        progress,
        "Method B — extract each feature’s DNA and blastn against the target…",
    )
    try:
        blastn = _blast_executable("blastn", use_system_blast=cfg.use_system_blast)
        makeblastdb = _blast_executable(
            "makeblastdb", use_system_blast=cfg.use_system_blast
        )
    except FileNotFoundError as exc:
        raise AnnotationTransferError(str(exc)) from exc

    query_records: list[tuple[str, str]] = []
    feature_by_qid: dict[str, SequenceFeature] = {}
    index_by_qid: dict[str, int] = {}
    total = len(indices)
    for step, index in enumerate(indices, start=1):
        if cancel is not None and cancel.requested:
            raise JobCancelledError("Annotation transfer cancelled.")
        feat = features[index]
        name = str(feat.name or feat.feature_type or f"feature {index + 1}")
        _emit(
            progress,
            f"Extracting feature DNA {step}/{total}: {name} "
            f"({feat.start}-{feat.end})",
        )
        dna = _feature_dna(source_sequence, feat)
        if len(dna) < 11:
            continue
        qid = f"f{index}"
        query_records.append((qid, dna))
        feature_by_qid[qid] = feat
        index_by_qid[qid] = index

    by_index: dict[int, LiftedAnnotation] = {
        index: LiftedAnnotation(
            source=features[index],
            status="lost",
            identity=0.0,
            coverage=0.0,
            method="blastn",
            feature=None,
        )
        for index in indices
    }
    if not query_records:
        _emit(progress, "No feature DNA long enough for blastn (≥11 bp).")
        return by_index

    ann_source = str(cfg.annotation_source or _DEFAULT_TRANSFER_SOURCE).strip()
    tgt_len = len(target_sequence)

    with tempfile.TemporaryDirectory(prefix="rv_annot_xfer_blast_") as tmp:
        tmp_path = Path(tmp)
        query_fa = tmp_path / "features.fa"
        subject_fa = tmp_path / "target.fa"
        db_prefix = tmp_path / "target_db"
        out_tsv = tmp_path / "hits.tsv"
        _write_fasta(query_fa, query_records)
        _write_fasta(subject_fa, [(target_id or "target", target_sequence)])

        _emit(
            progress,
            f"Indexing target ({len(target_sequence):,} bp) and running blastn "
            f"on {len(query_records)} feature(s)…",
        )
        try:
            run_logged_process(
                [
                    makeblastdb,
                    "-in",
                    str(subject_fa),
                    "-dbtype",
                    "nucl",
                    "-out",
                    str(db_prefix),
                ],
                tool="makeblastdb",
                inputs=[subject_fa],
                check=True,
                cancel=cancel,
            )
            run_logged_process(
                [
                    blastn,
                    "-query",
                    str(query_fa),
                    "-db",
                    str(db_prefix),
                    "-out",
                    str(out_tsv),
                    "-outfmt",
                    "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore",
                    "-evalue",
                    str(float(cfg.evalue)),
                    "-dust",
                    "no",
                    "-max_hsps",
                    "5",
                    "-max_target_seqs",
                    "5",
                    "-num_threads",
                    str(max(1, int(cfg.threads))),
                ],
                tool="blastn",
                inputs=[query_fa, subject_fa],
                check=True,
                output_path=out_tsv,
                cancel=cancel,
            )
        except JobCancelledError:
            raise
        except Exception as exc:  # noqa: BLE001
            raise AnnotationTransferError(
                f"blastn failed while transferring annotations.\n{str(exc).strip()[:500]}"
            ) from exc

        if cancel is not None and cancel.requested:
            raise JobCancelledError("Annotation transfer cancelled.")

        _emit(progress, "Parsing blastn hits…")
        rows = _parse_blast_rows(
            out_tsv.read_text(encoding="utf-8") if out_tsv.is_file() else ""
        )
        best_by_qid: dict[str, dict[str, str]] = {}
        for row in rows:
            qid = row["qseqid"]
            if qid not in feature_by_qid:
                continue
            try:
                sstart = int(row["sstart"])
                send = int(row["send"])
                qstart = int(row["qstart"])
                qend = int(row["qend"])
                pident = float(row["pident"]) / 100.0
                bits = float(row["bitscore"])
                aln_len = int(row["length"])
            except (TypeError, ValueError):
                continue
            reverse_hit = sstart > send
            if reverse_hit and not cfg.allow_reverse:
                continue
            feat = feature_by_qid[qid]
            qlen = max(1, abs(int(feat.end) - int(feat.start)) + 1)
            qcov = abs(qend - qstart) + 1
            coverage = min(1.0, max(aln_len, qcov) / float(qlen))
            if pident < float(cfg.min_identity) or coverage < float(cfg.min_coverage):
                continue
            prev = best_by_qid.get(qid)
            if prev is None or bits > float(prev["bitscore"]):
                best_by_qid[qid] = {
                    **row,
                    "_identity": f"{pident}",
                    "_coverage": f"{coverage}",
                }

        for qid, row in best_by_qid.items():
            index = index_by_qid[qid]
            feat = feature_by_qid[qid]
            identity = float(row["_identity"])
            coverage = float(row["_coverage"])
            sstart = int(row["sstart"])
            send = int(row["send"])
            reverse_hit = sstart > send
            status = _hit_status(
                identity=identity,
                coverage=coverage,
                min_identity=float(cfg.min_identity),
                min_coverage=float(cfg.min_coverage),
            )
            by_index[index] = _make_placed(
                src=feat,
                start=min(sstart, send),
                end=max(sstart, send),
                strand="-" if reverse_hit else "+",
                identity=identity,
                coverage=coverage,
                status=status,
                method="blastn",
                ann_source=ann_source,
                source_id=source_id,
                tgt_len=tgt_len,
            )

    placed = sum(1 for item in by_index.values() if item.feature is not None)
    _emit(progress, f"blastn placed {placed}/{len(indices)} feature(s).")
    return by_index


def transfer_annotations(
    *,
    source_id: str,
    source_sequence: str,
    target_id: str,
    target_sequence: str,
    features: list[SequenceFeature],
    opts: AnnotationTransferOptions | None = None,
    cancel: JobCancellation | None = None,
    source_fasta: Path | None = None,
    progress=None,
) -> AnnotationTransferResult:
    """Hybrid transfer: whole-seq lift first, then blastn for leftovers."""
    cfg = opts or AnnotationTransferOptions()
    strategy = str(cfg.strategy or "hybrid").strip().lower()
    if strategy not in {"hybrid", "liftover", "blastn"}:
        strategy = "hybrid"

    _emit(
        progress,
        "Transfer strategy: hybrid — (1) minimap2 sequence alignment + feature "
        "lift, then (2) blastn on remaining feature DNAs."
        if strategy == "hybrid"
        else (
            "Transfer strategy: liftover only (minimap2 sequence→sequence)."
            if strategy == "liftover"
            else "Transfer strategy: blastn only (feature DNA → target)."
        ),
    )
    if source_fasta is not None:
        _emit(progress, "Checking RepeatViewer annotation JSON…")
        require_annotation_json(source_fasta)
    if not source_sequence or not target_sequence:
        raise AnnotationTransferError("Source and target sequences must be non-empty.")
    if not features:
        return AnnotationTransferResult(message="No features to transfer.")
    if cancel is not None and cancel.requested:
        raise JobCancelledError("Annotation transfer cancelled.")

    lift_by_index: dict[int, LiftedAnnotation] = {}
    align_id = 0.0
    align_strand = "+"
    if strategy in {"hybrid", "liftover"}:
        lift_by_index, align_id, align_strand = _transfer_via_liftover(
            source_id=source_id,
            source_sequence=source_sequence,
            target_id=target_id,
            target_sequence=target_sequence,
            features=features,
            cfg=cfg,
            cancel=cancel,
            progress=progress,
        )

    need_blast: list[int] = []
    if strategy == "blastn":
        need_blast = list(range(len(features)))
    elif strategy == "hybrid":
        for index, src in enumerate(features):
            item = lift_by_index.get(index)
            if item is None or item.feature is None:
                need_blast.append(index)
        if need_blast:
            _emit(
                progress,
                f"{len(need_blast)} feature(s) not placed by lift — "
                "trying per-feature blastn…",
            )
    # liftover-only: no blastn

    blast_by_index: dict[int, LiftedAnnotation] = {}
    if need_blast and strategy in {"hybrid", "blastn"}:
        blast_by_index = _transfer_via_blastn(
            source_id=source_id,
            source_sequence=source_sequence,
            target_id=target_id,
            target_sequence=target_sequence,
            features=features,
            feature_indices=need_blast,
            cfg=cfg,
            cancel=cancel,
            progress=progress,
        )

    results: list[LiftedAnnotation] = []
    for index, src in enumerate(features):
        lift_item = lift_by_index.get(index)
        blast_item = blast_by_index.get(index)
        chosen: LiftedAnnotation | None = None
        if lift_item is not None and blast_item is not None:
            chosen = lift_item if _score(lift_item) >= _score(blast_item) else blast_item
        elif lift_item is not None:
            chosen = lift_item
        elif blast_item is not None:
            chosen = blast_item
        else:
            chosen = LiftedAnnotation(
                source=src,
                status="lost",
                identity=0.0,
                coverage=0.0,
                method="",
                feature=None,
            )
        results.append(chosen)

    counts = AnnotationTransferResult(lifted=results).counts()
    msg = (
        f"Hybrid transfer {source_id} → {target_id}: "
        f"placed {counts['transferred']} of {len(features)} "
        f"(liftover {counts.get('liftover', 0)}, blastn {counts.get('blastn', 0)}; "
        f"intact {counts.get('intact', 0)}, changed {counts.get('changed', 0)}, "
        f"partial {counts.get('partial', 0)}, lost {counts.get('lost', 0)}). "
        f"Whole-seq identity {align_id:.1%}."
        if strategy == "hybrid"
        else (
            f"Liftover transfer {source_id} → {target_id}: "
            f"placed {counts['transferred']} of {len(features)} "
            f"(whole-seq identity {align_id:.1%})."
            if strategy == "liftover"
            else (
                f"blastn transfer {source_id} → {target_id}: "
                f"placed {counts['transferred']} of {len(features)}."
            )
        )
    )
    ann = str(cfg.annotation_source or "").strip()
    if ann:
        msg = f"{msg} Source label: {ann}."
    _emit(progress, msg)
    return AnnotationTransferResult(
        lifted=results,
        message=msg,
        alignment_identity=align_id,
        strand=align_strand,
    )
