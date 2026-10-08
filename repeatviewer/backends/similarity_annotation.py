"""Similarity-based annotation against user-defined HMM / protein / DNA databases."""

from __future__ import annotations

import os
import platform
import subprocess
from repeatviewer.external_process import run_logged_process
import tempfile
from collections.abc import Callable
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass
from pathlib import Path

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from repeatviewer.backends.thread_budget import analysis_threads_per_job, load_batch_workers
from repeatviewer.backends.auto_annotation_settings import OrfSettings
from repeatviewer.backends.pyhmmer_search import collect_pyhmmer_domain_hits, hmm_uses_nucleotide_query
from repeatviewer.backends.similarity_annotation_settings import (
    SimilarityAnnotationSettings,
    similarity_feature_name,
)
from repeatviewer.backends.database_prepare import (
    database_search_readiness,
    prepared_cdd_prefix,
    pyhmmer_install_hint,
)
from repeatviewer.backends.gui_functions import check_cmd_in_path, get_original_file_path
from repeatviewer.backends.job_cancel import JobCancellation, JobCancelledError, check_cancelled
from repeatviewer.backends.orf_detect import OrfDetectResult, detect_orf_features
from repeatviewer.backends.sequence_features import (
    REPEATVIEWER_INTERNAL_SOURCE,
    SequenceFeature,
    load_sequence_map_data,
)
from repeatviewer.project import (
    DATABASE_TYPE_CDD,
    DATABASE_TYPE_DNA,
    DATABASE_TYPE_PROTEIN,
    UserDatabase,
)

SIMILARITY_DB_COLORS: tuple[str, ...] = (
    "#e41a1c",
    "#377eb8",
    "#4daf4a",
    "#984ea3",
    "#ff7f00",
    "#a65628",
    "#f781bf",
    "#1b9e77",
    "#7570b3",
    "#d95f02",
)
# Default E-value cutoffs for similarity annotation (user-configurable in settings).
DEFAULT_BLAST_EVALUE = 1e-4
DEFAULT_HMM_EVALUE = 1e-4
DEFAULT_CDD_EVALUE = 1e-4
DEFAULT_MAX_HITS_PER_DATABASE = 10
BLAST_SEARCH_MAX_TARGET_SEQS = 500


@dataclass(frozen=True)
class SimilarityDatabaseRun:
    database: UserDatabase
    color: str


@dataclass(frozen=True)
class SimilarityAnnotationResult:
    features: tuple[SequenceFeature, ...]
    message: str


def color_for_database_index(index: int, *, database_name: str | None = None) -> str:
    from repeatviewer.backends.annotation_colors import similarity_database_color

    return similarity_database_color(database_name or "", index)


def _os_type() -> str:
    return platform.system()


def _blast_executable(name: str, *, use_system_blast: bool) -> str:
    if use_system_blast and check_cmd_in_path(name):
        return name
    blast_dir = Path(get_original_file_path()) / "blast"
    suffix = "_linux" if _os_type() == "Linux" else "_mac" if _os_type() == "Darwin" else ".exe"
    bundled = blast_dir / f"{name}{suffix}"
    if bundled.is_file():
        return str(bundled)
    if check_cmd_in_path(name):
        return name
    raise FileNotFoundError(
        f"'{name}' was not found. Install BLAST+ or enable system BLAST in Preferences."
    )


def similarity_search_needs_orfs(databases: tuple[SimilarityDatabaseRun, ...]) -> bool:
    """Whether any selected database searches translated ORF products."""
    for entry in databases:
        db = entry.database
        if db.db_type == DATABASE_TYPE_CDD:
            continue
        if db.is_hmm:
            if hmm_uses_nucleotide_query(db.path):
                continue
            return True
        if db.db_type == DATABASE_TYPE_PROTEIN:
            return True
    return False


def _orfs_from_features(features: tuple[SequenceFeature, ...]) -> tuple[SequenceFeature, ...]:
    return tuple(
        feat
        for feat in features
        if str(feat.feature_type).upper() == "ORF" or str(feat.name or "").upper() == "ORF"
    )


def _resolve_orfs_for_similarity(
    *,
    sequence_id: str,
    sequence: str,
    features: list[SequenceFeature],
    orf_settings: OrfSettings | None,
) -> OrfDetectResult:
    existing = _orfs_from_features(tuple(features))
    if existing:
        return OrfDetectResult(
            existing,
            f"Using {len(existing)} ORF(s) from the current annotation.",
        )
    return detect_orf_features(sequence_id, sequence, settings=orf_settings)


def _orf_id(orf: SequenceFeature) -> str:
    return f"orf_{orf.start}_{orf.end}_{orf.strand}"


def _parse_orf_id(query_id: str) -> tuple[int, int, str] | None:
    parts = str(query_id).split("_")
    if len(parts) < 4 or parts[0] != "orf":
        return None
    try:
        return int(parts[1]), int(parts[2]), parts[3]
    except ValueError:
        return None


def _map_protein_coords_to_nucleotide(
    orf: SequenceFeature,
    qstart: int,
    qend: int,
) -> tuple[int, int] | None:
    qstart = int(qstart)
    qend = int(qend)
    if qstart <= 0 or qend <= 0:
        return None
    if qstart > qend:
        qstart, qend = qend, qstart
    if orf.strand == "-":
        nt_end = orf.end - (qstart - 1) * 3
        nt_start = orf.end - qend * 3 + 1
    else:
        nt_start = orf.start + (qstart - 1) * 3
        nt_end = orf.start + qend * 3 - 1
    nt_start = max(orf.start, min(orf.end, nt_start))
    nt_end = max(orf.start, min(orf.end, nt_end))
    if nt_start > nt_end:
        nt_start, nt_end = nt_end, nt_start
    return nt_start, nt_end


def _protein_records_from_orfs(
    sequence: str,
    orfs: tuple[SequenceFeature, ...],
    *,
    min_protein_length: int = 10,
) -> list[SeqRecord]:
    records: list[SeqRecord] = []
    for orf in orfs:
        subseq = sequence[orf.start - 1 : orf.end]
        if orf.strand == "-":
            subseq = str(Seq(subseq).reverse_complement())
        protein = str(Seq(subseq).translate(to_stop=False)).rstrip("*")
        if len(protein) < min_protein_length:
            continue
        records.append(
            SeqRecord(
                Seq(protein),
                id=_orf_id(orf),
                description=orf.name or "ORF",
            )
        )
    return records


def _cdd_domain_label(stitle: str, sseqid: str = "") -> str:
    """Extract the CDD domain name from rpstblastn stitle (same rule as TE-Aid)."""
    description = (stitle or sseqid or "hit").strip()
    if not description:
        return "hit"
    parts = [part.strip() for part in description.split(",") if part.strip()]
    if len(parts) > 1:
        return parts[1]
    return parts[0]


def _clamp_evalue(value: float, *, default: float) -> float:
    if value <= 0.0:
        return default
    return max(value, 1e-200)


def _max_hits(value: int) -> int:
    return max(1, int(value))


def _top_blast_hits(hits: list[dict[str, str]], *, max_hits: int) -> list[dict[str, str]]:
    """Return up to max_hits BLAST rows ranked by bit score."""
    ranked: list[tuple[float, dict[str, str]]] = []
    for hit in hits:
        try:
            bitscore = float(hit["bitscore"])
        except (KeyError, TypeError, ValueError):
            continue
        ranked.append((bitscore, hit))
    ranked.sort(key=lambda row: row[0], reverse=True)
    limit = _max_hits(max_hits)
    return [hit for _score, hit in ranked[:limit]]


def _top_cdd_hits(
    hits: list[dict[str, str]],
    *,
    max_evalue: float,
    max_hits: int,
) -> list[dict[str, str]]:
    """Return up to max_hits CDD rows ranked by bit score within the E-value cutoff."""
    ranked: list[tuple[float, dict[str, str]]] = []
    for hit in hits:
        try:
            evalue = float(hit["evalue"])
            score = float(hit["score"])
        except (KeyError, TypeError, ValueError):
            continue
        if evalue > max_evalue:
            continue
        ranked.append((score, hit))
    ranked.sort(key=lambda row: row[0], reverse=True)
    limit = _max_hits(max_hits)
    return [hit for _score, hit in ranked[:limit]]


def _top_pyhmmer_domain_hits(
    raw_hits: list[tuple[float, float | None, str, str, int, int]],
    *,
    max_hits: int,
) -> list[tuple[str, str, int, int, float | None]]:
    """Return up to max_hits pyhmmer domains ranked by bit score."""
    sorted_hits = sorted(raw_hits, key=lambda row: row[0], reverse=True)
    limit = _max_hits(max_hits)
    return [
        (qname, hmm_name, ali_from, ali_to, evalue)
        for _score, evalue, qname, hmm_name, ali_from, ali_to in sorted_hits[:limit]
    ]


def _blast_hit_bitscore(hit: dict[str, str]) -> float | None:
    try:
        return float(hit["bitscore"])
    except (KeyError, TypeError, ValueError):
        return None


def _cdd_hit_score(hit: dict[str, str]) -> float | None:
    try:
        return float(hit["score"])
    except (KeyError, TypeError, ValueError):
        return None


def _cdd_hit_passes_evalue(hit: dict[str, str], *, max_evalue: float) -> bool:
    try:
        return float(hit["evalue"]) <= max_evalue
    except (KeyError, TypeError, ValueError):
        return False


def _interval_contains(outer: tuple[int, int], inner: tuple[int, int]) -> bool:
    o_start, o_end = sorted(outer)
    i_start, i_end = sorted(inner)
    return o_start <= i_start and o_end >= i_end


def _intervals_overlap(left: tuple[int, int], right: tuple[int, int]) -> bool:
    l_start, l_end = sorted(left)
    r_start, r_end = sorted(right)
    return not (l_end < r_start or r_end < l_start)


def _intervals_partially_overlap(left: tuple[int, int], right: tuple[int, int]) -> bool:
    return (
        _intervals_overlap(left, right)
        and not _interval_contains(left, right)
        and not _interval_contains(right, left)
    )


def _normalize_interval(start: int, end: int) -> tuple[int, int]:
    start_i, end_i = int(start), int(end)
    if start_i <= end_i:
        return start_i, end_i
    return end_i, start_i


def _interval_length(interval: tuple[int, int]) -> int:
    start, end = sorted(interval)
    return max(0, end - start + 1)


def _payload_evalue(payload: object) -> float:
    """Extract E-value from a BLAST hit dict or pyhmmer payload tuple."""
    if isinstance(payload, dict):
        try:
            return float(payload["evalue"])
        except (KeyError, TypeError, ValueError):
            return float("inf")
    if isinstance(payload, tuple) and len(payload) >= 5:
        evalue = payload[4]
        if evalue is None:
            return float("inf")
        try:
            return float(evalue)
        except (TypeError, ValueError):
            return float("inf")
    return float("inf")


def _dedupe_identical_intervals(
    scored: list[tuple[int, int, float, object]],
) -> list[tuple[int, int, float, object]]:
    """Keep one hit per exact genomic interval; highest bit score wins."""
    best: dict[tuple[int, int], tuple[int, int, float, object]] = {}
    for start, end, score, payload in scored:
        interval = _normalize_interval(start, end)
        row = (interval[0], interval[1], float(score), payload)
        previous = best.get(interval)
        if previous is None or row[2] > previous[2]:
            best[interval] = row
    return list(best.values())


def _filter_full_overlaps(
    scored: list[tuple[int, int, float, object]],
) -> list[tuple[int, int, float, object]]:
    """Remove shorter hits fully nested inside a longer hit when E-value is not better.

    When a longer hit fully covers a shorter one and the shorter hit's E-value is
    greater than or equal to the covering hit's E-value, the shorter hit is dropped.
    This applies even when both E-values are 0.
    """
    if not scored:
        return []

    ranked = sorted(
        scored,
        key=lambda row: (_interval_length((row[0], row[1])), row[2]),
        reverse=True,
    )
    kept: list[tuple[int, int, float, object]] = []
    for start, end, score, payload in ranked:
        interval = (start, end)
        length = _interval_length(interval)
        inner_evalue = _payload_evalue(payload)
        skip = False
        remove_indices: list[int] = []
        for idx, (kept_start, kept_end, kept_score, kept_payload) in enumerate(kept):
            kept_interval = (kept_start, kept_end)
            kept_length = _interval_length(kept_interval)
            kept_evalue = _payload_evalue(kept_payload)
            if (
                _interval_contains(kept_interval, interval)
                and kept_length > length
                and inner_evalue >= kept_evalue
            ):
                skip = True
                break
            if (
                _interval_contains(interval, kept_interval)
                and length > kept_length
                and kept_evalue >= inner_evalue
            ):
                remove_indices.append(idx)
        for idx in reversed(remove_indices):
            kept.pop(idx)
        if not skip:
            kept.append((start, end, score, payload))
    return kept


def _cluster_partial_overlaps_keep_best(
    scored: list[tuple[int, int, float, object]],
) -> list[tuple[int, int, float, object]]:
    """Merge partially overlapping hits; keep the highest-scoring hit in each group."""
    if len(scored) <= 1:
        return list(scored)

    count = len(scored)
    parent = list(range(count))

    def find(index: int) -> int:
        while parent[index] != index:
            parent[index] = parent[parent[index]]
            index = parent[index]
        return index

    def union(left: int, right: int) -> None:
        root_left = find(left)
        root_right = find(right)
        if root_left != root_right:
            parent[root_right] = root_left

    for left in range(count):
        left_interval = (scored[left][0], scored[left][1])
        for right in range(left + 1, count):
            right_interval = (scored[right][0], scored[right][1])
            if _intervals_partially_overlap(left_interval, right_interval):
                union(left, right)

    groups: dict[int, list[tuple[int, int, float, object]]] = {}
    for index, row in enumerate(scored):
        groups.setdefault(find(index), []).append(row)

    return [max(group, key=lambda item: item[2]) for group in groups.values()]


def _cap_scored_intervals(
    scored: list[tuple[int, int, float, object]],
    *,
    max_hits: int,
) -> list[tuple[int, int, float, object]]:
    limit = _max_hits(max_hits)
    ranked = sorted(scored, key=lambda row: row[2], reverse=True)
    return ranked[:limit]


def _select_scored_intervals(
    scored: list[tuple[int, int, float, object]],
    *,
    max_hits: int,
) -> list[tuple[int, int, float, object]]:
    """Hit selection: identical intervals, full-overlap filter, partial clustering, cap."""
    if not scored:
        return []
    limit = _max_hits(max_hits)
    scored = _dedupe_identical_intervals(scored)
    after_full = _filter_full_overlaps(scored)
    if len(after_full) <= limit:
        return after_full
    after_partial = _cluster_partial_overlaps_keep_best(after_full)
    if len(after_partial) <= limit:
        return after_partial
    return _cap_scored_intervals(after_partial, max_hits=limit)


def _blast_hit_interval(hit: dict[str, str]) -> tuple[int, int] | None:
    try:
        start = int(hit["qstart"])
        end = int(hit["qend"])
    except (KeyError, TypeError, ValueError):
        return None
    return start, end


def _blast_hit_nt_interval(
    hit: dict[str, str],
    lookup: dict[str, SequenceFeature],
) -> tuple[int, int] | None:
    orf_key = str(hit.get("qseqid", "")).strip()
    orf = lookup.get(orf_key)
    if orf is None:
        parsed = _parse_orf_id(orf_key)
        if parsed is None:
            return None
        start, end, strand = parsed
        orf = SequenceFeature(start=start, end=end, name="ORF", feature_type="ORF", strand=strand)
    try:
        mapped = _map_protein_coords_to_nucleotide(
            orf,
            int(hit["qstart"]),
            int(hit["qend"]),
        )
    except (KeyError, TypeError, ValueError):
        return None
    return mapped


def _blast_hits_to_scored(
    hits: list[dict[str, str]],
    *,
    score_getter=_blast_hit_bitscore,
    interval_getter=_blast_hit_interval,
) -> list[tuple[int, int, float, dict[str, str]]]:
    scored: list[tuple[int, int, float, dict[str, str]]] = []
    for hit in hits:
        score = score_getter(hit)
        interval = interval_getter(hit)
        if score is None or interval is None:
            continue
        start, end = interval
        scored.append((start, end, score, hit))
    return scored


def _dedupe_blast_hits_by_interval(
    hits: list[dict[str, str]],
    *,
    interval_getter,
    score_getter=_blast_hit_bitscore,
) -> list[dict[str, str]]:
    """Drop duplicate BLAST rows that cover the same interval; keep top bit score."""
    best: dict[tuple[int, int], tuple[float, dict[str, str]]] = {}
    for hit in hits:
        interval = interval_getter(hit)
        score = score_getter(hit)
        if interval is None or score is None:
            continue
        key = _normalize_interval(*interval)
        current = best.get(key)
        if current is None or score > current[0]:
            best[key] = (score, hit)
    ranked = sorted(best.values(), key=lambda row: row[0], reverse=True)
    return [hit for _score, hit in ranked]


def _dedupe_domain_features_by_interval(
    features: list[SequenceFeature],
    *,
    scores: list[float],
) -> list[SequenceFeature]:
    """Keep one domain feature per genomic interval; highest score wins."""
    if len(features) != len(scores):
        return features
    best: dict[tuple[int, int, str], tuple[float, SequenceFeature]] = {}
    for feat, score in zip(features, scores):
        key = (feat.start, feat.end, feat.strand)
        current = best.get(key)
        if current is None or score > current[0]:
            best[key] = (score, feat)
    return [feat for _score, feat in sorted(best.values(), key=lambda row: row[1].start)]


def _select_blast_hits(
    hits: list[dict[str, str]],
    *,
    max_hits: int,
    score_getter=_blast_hit_bitscore,
    interval_getter=_blast_hit_interval,
) -> list[dict[str, str]]:
    selected = _select_scored_intervals(
        _blast_hits_to_scored(
            hits,
            score_getter=score_getter,
            interval_getter=interval_getter,
        ),
        max_hits=max_hits,
    )
    return _dedupe_blast_hits_by_interval(
        [hit for _start, _end, _score, hit in selected],
        interval_getter=interval_getter,
        score_getter=score_getter,
    )


def _select_blast_protein_hits_on_genome(
    hits: list[dict[str, str]],
    lookup: dict[str, SequenceFeature],
    *,
    max_hits: int,
) -> list[dict[str, str]]:
    scored: list[tuple[int, int, float, dict[str, str]]] = []
    for hit in hits:
        score = _blast_hit_bitscore(hit)
        interval = _blast_hit_nt_interval(hit, lookup)
        if score is None or interval is None:
            continue
        start, end = interval
        scored.append((start, end, score, hit))
    selected = _select_scored_intervals(scored, max_hits=max_hits)
    return _dedupe_blast_hits_by_interval(
        [hit for _start, _end, _score, hit in selected],
        interval_getter=lambda row: _blast_hit_nt_interval(row, lookup),
        score_getter=_blast_hit_bitscore,
    )


def _select_pyhmmer_protein_hits_on_genome(
    raw_hits: list[tuple[float, float | None, str, str, int, int]],
    lookup: dict[str, SequenceFeature],
    *,
    max_hits: int,
) -> list[tuple[str, str, int, int, float | None]]:
    scored: list[tuple[int, int, float, tuple[str, str, int, int, float | None]]] = []
    for score, evalue, qname, hmm_name, ali_from, ali_to in raw_hits:
        orf = lookup.get(qname)
        if orf is None:
            parsed = _parse_orf_id(qname)
            if parsed is None:
                continue
            start, end, strand = parsed
            orf = SequenceFeature(start=start, end=end, name="ORF", feature_type="ORF", strand=strand)
        mapped = _map_protein_coords_to_nucleotide(orf, ali_from, ali_to)
        if mapped is None:
            continue
        nt_start, nt_end = mapped
        payload = (qname, hmm_name, ali_from, ali_to, evalue)
        scored.append((nt_start, nt_end, float(score), payload))
    selected = _select_scored_intervals(scored, max_hits=max_hits)
    return [payload for _start, _end, _score, payload in selected]


def _select_blast_hits_by_coverage(
    hits: list[dict[str, str]],
    *,
    max_hits: int,
    score_getter=_blast_hit_bitscore,
) -> list[dict[str, str]]:
    return _select_blast_hits(
        hits,
        max_hits=max_hits,
        score_getter=score_getter,
    )


def _pyhmmer_rows_to_scored(
    raw_hits: list[tuple[float, float | None, str, str, int, int]],
) -> list[tuple[int, int, float, tuple[str, str, int, int, float | None]]]:
    return [
        (
            int(ali_from),
            int(ali_to),
            float(score),
            (qname, hmm_name, ali_from, ali_to, evalue),
        )
        for score, evalue, qname, hmm_name, ali_from, ali_to in raw_hits
    ]


def _select_pyhmmer_hits(
    raw_hits: list[tuple[float, float | None, str, str, int, int]],
    *,
    max_hits: int,
) -> list[tuple[str, str, int, int, float | None]]:
    selected = _select_scored_intervals(_pyhmmer_rows_to_scored(raw_hits), max_hits=max_hits)
    return [payload for _start, _end, _score, payload in selected]


def _select_pyhmmer_hits_by_coverage(
    raw_hits: list[tuple[float, float | None, str, str, int, int]],
    *,
    max_hits: int,
) -> list[tuple[str, str, int, int, float | None]]:
    return _select_pyhmmer_hits(raw_hits, max_hits=max_hits)


def _search_max_target_seqs(*, display_limit: int) -> int:
    return max(BLAST_SEARCH_MAX_TARGET_SEQS, _max_hits(display_limit))


def _run_blast_tabular(
    *,
    program: str,
    query_path: Path,
    db_path: Path,
    use_system_blast: bool,
    report_evalue: float | str = DEFAULT_BLAST_EVALUE,
    max_target_seqs: int = DEFAULT_MAX_HITS_PER_DATABASE,
    num_threads: int = 1,
    cancel: JobCancellation | None = None,
) -> list[dict[str, str]]:
    exe = _blast_executable(program, use_system_blast=use_system_blast)
    with tempfile.TemporaryDirectory(prefix="repeatviewer_simann_") as tmp:
        out_path = Path(tmp) / "hits.tsv"
        cmd = [
            exe,
            "-query",
            str(query_path),
            "-db",
            str(db_path),
            "-out",
            str(out_path),
            "-outfmt",
            "6 qseqid sseqid pident length qstart qend sstart send evalue bitscore stitle",
            "-evalue",
            str(report_evalue),
            "-max_target_seqs",
            str(_max_hits(max_target_seqs)),
            "-num_threads",
            str(max(1, num_threads)),
        ]
        run_logged_process(
            cmd,
            tool=program,
            inputs=[query_path, db_path],
            check=True,
            output_path=out_path,
            cancel=cancel,
        )
        if not out_path.is_file() or out_path.stat().st_size == 0:
            return []
        rows: list[dict[str, str]] = []
        for line in out_path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 11:
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
                    "stitle": fields[10] if len(fields) == 11 else " ".join(fields[10:]),
                }
            )
        return rows


def _orf_lookup(orfs: tuple[SequenceFeature, ...]) -> dict[str, SequenceFeature]:
    return {_orf_id(orf): orf for orf in orfs}


def _domain_features_from_protein_hits(
    hits: list[tuple[str, str, int, int, float | None]],
    lookup: dict[str, SequenceFeature],
    *,
    database_name: str,
    show_database_name: bool = True,
) -> list[SequenceFeature]:
    features: list[SequenceFeature] = []
    seen: set[tuple[int, int, str]] = set()
    for orf_key, target_name, ali_from, ali_to, evalue in hits:
        orf = lookup.get(orf_key)
        if orf is None:
            parsed = _parse_orf_id(orf_key)
            if parsed is None:
                continue
            start, end, strand = parsed
            orf = SequenceFeature(
                start=start,
                end=end,
                name="ORF",
                feature_type="ORF",
                strand=strand,
            )
        mapped = _map_protein_coords_to_nucleotide(orf, ali_from, ali_to)
        if mapped is None:
            continue
        nt_start, nt_end = mapped
        label = similarity_feature_name(
            database_name,
            target_name,
            show_database_name=show_database_name,
            evalue=evalue,
        )
        interval_key = (nt_start, nt_end, orf.strand)
        if interval_key in seen:
            continue
        seen.add(interval_key)
        features.append(
            SequenceFeature(
                start=nt_start,
                end=nt_end,
                name=label,
                feature_type="domain",
                strand=orf.strand,
            )
        )
    return features


def _domain_features_from_nucleotide_hmm_hits(
    hits: list[tuple[str, str, int, int, float | None]],
    *,
    database_name: str,
    show_database_name: bool = True,
) -> list[SequenceFeature]:
    features: list[SequenceFeature] = []
    seen: set[tuple[int, int, str]] = set()
    for _query_id, hmm_name, ali_from, ali_to, evalue in hits:
        start, end = sorted((int(ali_from), int(ali_to)))
        label = similarity_feature_name(
            database_name,
            hmm_name,
            show_database_name=show_database_name,
            evalue=evalue,
        )
        interval_key = (start, end, "+")
        if interval_key in seen:
            continue
        seen.add(interval_key)
        features.append(
            SequenceFeature(
                start=start,
                end=end,
                name=label,
                feature_type="domain",
                strand="+",
            )
        )
    return features


def _annotate_protein_database(
    sequence: str,
    orfs: tuple[SequenceFeature, ...],
    database: UserDatabase,
    *,
    settings: SimilarityAnnotationSettings,
    use_system_blast: bool,
    blast_threads: int = 1,
    cancel: JobCancellation | None = None,
) -> tuple[tuple[SequenceFeature, ...], str]:
    from repeatviewer.backends.database_prepare import blast_search_db_path

    db_path = blast_search_db_path(
        database.path.expanduser().resolve(), protein=True
    )
    readiness = database_search_readiness(database)
    if readiness is not None:
        return (), f"{database.name}: {readiness}"
    proteins = _protein_records_from_orfs(sequence, orfs)
    if not proteins:
        return (), f"{database.name}: no ORFs long enough for BLASTP."

    lookup = _orf_lookup(orfs)
    report_evalue = _clamp_evalue(settings.blastp_evalue, default=DEFAULT_BLAST_EVALUE)
    max_hits = _max_hits(settings.max_hits_per_database)
    with tempfile.TemporaryDirectory(prefix="repeatviewer_simann_protein_") as tmp:
        query_faa = Path(tmp) / "orfs.faa"
        SeqIO.write(proteins, str(query_faa), "fasta")
        try:
            hits = _run_blast_tabular(
                program="blastp",
                query_path=query_faa,
                db_path=db_path,
                use_system_blast=use_system_blast,
                report_evalue=report_evalue,
                max_target_seqs=_search_max_target_seqs(display_limit=max_hits),
                num_threads=blast_threads,
                cancel=cancel,
            )
        except JobCancelledError:
            raise
        except FileNotFoundError as exc:
            return (), str(exc)
        except subprocess.CalledProcessError as exc:
            log_process_output("rpstblastn", stdout=exc.stdout, stderr=exc.stderr)
            detail = (exc.stderr or exc.stdout or "").strip()
            return (), f"{database.name}: BLASTP failed.{chr(10)}{detail}"

    top_hits = _select_blast_protein_hits_on_genome(hits, lookup, max_hits=max_hits)
    if not top_hits:
        return (), f"{database.name}: no BLASTP hits (E-value ≤ {report_evalue:g})."

    features: list[SequenceFeature] = []
    feature_scores: list[float] = []
    for hit in top_hits:
        orf_key = hit.get("qseqid", "")
        orf = lookup.get(orf_key)
        if orf is None:
            parsed = _parse_orf_id(orf_key)
            if parsed is None:
                continue
            start, end, strand = parsed
            orf = SequenceFeature(start=start, end=end, name="ORF", feature_type="ORF", strand=strand)
        try:
            mapped = _map_protein_coords_to_nucleotide(
                orf,
                int(hit["qstart"]),
                int(hit["qend"]),
            )
        except (KeyError, TypeError, ValueError):
            continue
        if mapped is None:
            continue
        nt_start, nt_end = mapped
        target = hit.get("stitle") or hit.get("sseqid") or "hit"
        target = target.split()[0] if target else "hit"
        bitscore = _blast_hit_bitscore(hit)
        if bitscore is None:
            continue
        features.append(
            SequenceFeature(
                start=nt_start,
                end=nt_end,
                name=similarity_feature_name(
                    database.name,
                    target,
                    show_database_name=settings.show_database_name_on_track,
                    evalue=hit.get("evalue"),
                ),
                feature_type="domain",
                strand=orf.strand,
            )
        )
        feature_scores.append(bitscore)
    features = _dedupe_domain_features_by_interval(features, scores=feature_scores)
    if not features:
        return (), f"{database.name}: no BLASTP hits (E-value ≤ {report_evalue:g})."
    return (
        tuple(features),
        f"{database.name}: {len(features)} BLASTP hit(s) (E-value ≤ {report_evalue:g}).",
    )


def _annotate_dna_database(
    sequence: str,
    seq_id: str,
    database: UserDatabase,
    *,
    settings: SimilarityAnnotationSettings,
    use_system_blast: bool,
    blast_threads: int = 1,
    cancel: JobCancellation | None = None,
) -> tuple[tuple[SequenceFeature, ...], str]:
    from repeatviewer.backends.database_prepare import blast_search_db_path

    db_path = blast_search_db_path(
        database.path.expanduser().resolve(), protein=False
    )
    readiness = database_search_readiness(database)
    if readiness is not None:
        return (), f"{database.name}: {readiness}"
    report_evalue = _clamp_evalue(settings.blastn_evalue, default=DEFAULT_BLAST_EVALUE)
    max_hits = _max_hits(settings.max_hits_per_database)
    with tempfile.TemporaryDirectory(prefix="repeatviewer_simann_dna_") as tmp:
        query_fa = Path(tmp) / "query.fa"
        SeqIO.write([SeqRecord(Seq(sequence), id=seq_id)], str(query_fa), "fasta")
        try:
            hits = _run_blast_tabular(
                program="blastn",
                query_path=query_fa,
                db_path=db_path,
                use_system_blast=use_system_blast,
                report_evalue=report_evalue,
                max_target_seqs=_search_max_target_seqs(display_limit=max_hits),
                num_threads=blast_threads,
                cancel=cancel,
            )
        except JobCancelledError:
            raise
        except FileNotFoundError as exc:
            return (), str(exc)
        except subprocess.CalledProcessError as exc:
            detail = (exc.stderr or exc.stdout or "").strip()
            return (), f"{database.name}: BLASTN failed.{chr(10)}{detail}"

    top_hits = _select_blast_hits_by_coverage(hits, max_hits=max_hits)
    if not top_hits:
        return (), f"{database.name}: no BLASTN hits (E-value ≤ {report_evalue:g})."

    features: list[SequenceFeature] = []
    feature_scores: list[float] = []
    for hit in top_hits:
        try:
            qstart = int(hit["qstart"])
            qend = int(hit["qend"])
        except (KeyError, TypeError, ValueError):
            continue
        if qstart > qend:
            qstart, qend = qend, qstart
        bitscore = _blast_hit_bitscore(hit)
        if bitscore is None:
            continue
        target = hit.get("stitle") or hit.get("sseqid") or "hit"
        target = target.split()[0] if target else "hit"
        features.append(
            SequenceFeature(
                start=qstart,
                end=qend,
                name=similarity_feature_name(
                    database.name,
                    target,
                    show_database_name=settings.show_database_name_on_track,
                    evalue=hit.get("evalue"),
                ),
                feature_type="domain",
                strand="+",
            )
        )
        feature_scores.append(bitscore)
    features = _dedupe_domain_features_by_interval(features, scores=feature_scores)
    if not features:
        return (), f"{database.name}: no BLASTN hits (E-value ≤ {report_evalue:g})."
    return (
        tuple(features),
        f"{database.name}: {len(features)} BLASTN hit(s) (E-value ≤ {report_evalue:g}).",
    )


def _annotate_cdd_database(
    sequence: str,
    seq_id: str,
    database: UserDatabase,
    *,
    settings: SimilarityAnnotationSettings,
    use_system_blast: bool,
    blast_threads: int = 1,
    cancel: JobCancellation | None = None,
) -> tuple[tuple[SequenceFeature, ...], str]:
    cdd_evalue = _clamp_evalue(settings.cdd_evalue, default=DEFAULT_CDD_EVALUE)
    max_hits = _max_hits(settings.max_hits_per_database)
    cdd_dir = database.path.expanduser().resolve()
    readiness = database_search_readiness(database)
    if readiness is not None:
        return (), f"{database.name}: {readiness}"
    cdd_profile = prepared_cdd_prefix(cdd_dir)

    with tempfile.TemporaryDirectory(prefix="repeatviewer_simann_cdd_") as tmp:
        query_fa = Path(tmp) / "query.fa"
        out_path = Path(tmp) / "hits.tsv"
        SeqIO.write([SeqRecord(Seq(sequence), id=seq_id)], str(query_fa), "fasta")
        try:
            exe = _blast_executable("rpstblastn", use_system_blast=use_system_blast)
            cmd = [
                exe,
                "-query",
                str(query_fa),
                "-db",
                str(cdd_profile),
                "-out",
                str(out_path),
                "-outfmt",
                "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue score stitle",
                "-evalue",
                str(cdd_evalue),
                "-max_target_seqs",
                str(_search_max_target_seqs(display_limit=max_hits)),
                "-num_threads",
                str(max(1, blast_threads)),
            ]
            run_logged_process(
                cmd,
                tool="rpstblastn",
                inputs=[query_fa, cdd_profile],
                check=True,
                output_path=out_path,
                cancel=cancel,
            )
        except JobCancelledError:
            raise
        except FileNotFoundError as exc:
            return (), str(exc)
        except subprocess.CalledProcessError as exc:
            detail = (exc.stderr or exc.stdout or "").strip()
            return (), f"{database.name}: rpstblastn failed.{chr(10)}{detail}"

        if not out_path.is_file() or out_path.stat().st_size == 0:
            return (), f"{database.name}: no CDD hits (E-value ≤ {cdd_evalue:g})."

        hits: list[dict[str, str]] = []
        for line in out_path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 13:
                continue
            hits.append(
                {
                    "qseqid": fields[0],
                    "sseqid": fields[1],
                    "qstart": fields[6],
                    "qend": fields[7],
                    "evalue": fields[10],
                    "score": fields[11],
                    "stitle": fields[12] if len(fields) == 13 else " ".join(fields[12:]),
                }
            )

    top_hits = _select_blast_hits_by_coverage(
        [
            hit
            for hit in hits
            if _cdd_hit_passes_evalue(hit, max_evalue=cdd_evalue)
        ],
        max_hits=max_hits,
        score_getter=_cdd_hit_score,
    )
    if not top_hits:
        return (), f"{database.name}: no CDD hits (E-value ≤ {cdd_evalue:g})."

    features: list[SequenceFeature] = []
    feature_scores: list[float] = []
    for hit in top_hits:
        try:
            qstart = int(hit["qstart"])
            qend = int(hit["qend"])
        except (KeyError, TypeError, ValueError):
            continue
        if qstart > qend:
            qstart, qend = qend, qstart
            strand = "-"
        else:
            strand = "+"
        score = _cdd_hit_score(hit)
        if score is None:
            continue
        domain_name = _cdd_domain_label(
            str(hit.get("stitle", "")),
            str(hit.get("sseqid", "")),
        )
        features.append(
            SequenceFeature(
                start=qstart,
                end=qend,
                name=similarity_feature_name(
                    database.name,
                    domain_name,
                    show_database_name=settings.show_database_name_on_track,
                    evalue=hit.get("evalue"),
                ),
                feature_type="domain",
                strand=strand,
            )
        )
        feature_scores.append(score)
    features = _dedupe_domain_features_by_interval(features, scores=feature_scores)
    if not features:
        return (), f"{database.name}: no CDD hits (E-value ≤ {cdd_evalue:g})."
    return (
        tuple(features),
        f"{database.name}: {len(features)} CDD hit(s) (E-value ≤ {cdd_evalue:g}).",
    )


def _annotate_protein_hmm_with_pyhmmer(
    sequence: str,
    orfs: tuple[SequenceFeature, ...],
    database: UserDatabase,
    *,
    settings: SimilarityAnnotationSettings,
) -> tuple[tuple[SequenceFeature, ...], str]:
    hmm_path = database.path.expanduser().resolve()
    proteins = _protein_records_from_orfs(sequence, orfs)
    if not proteins:
        return (), f"{database.name}: no ORFs long enough for HMM search."

    lookup = _orf_lookup(orfs)
    max_hits = _max_hits(settings.max_hits_per_database)
    with tempfile.TemporaryDirectory(prefix="repeatviewer_simann_hmm_") as tmp:
        query_faa = Path(tmp) / "orfs.faa"
        SeqIO.write(proteins, str(query_faa), "fasta")
        raw_hits = collect_pyhmmer_domain_hits(
            hmm_path,
            query_faa,
            report_e=_clamp_evalue(settings.hmm_evalue, default=DEFAULT_HMM_EVALUE),
        )
    selected = _select_pyhmmer_protein_hits_on_genome(raw_hits, lookup, max_hits=max_hits)
    if not selected:
        report_e = _clamp_evalue(settings.hmm_evalue, default=DEFAULT_HMM_EVALUE)
        return (), f"{database.name}: no HMM hits (E-value ≤ {report_e:g})."
    features = _domain_features_from_protein_hits(
        selected,
        lookup,
        database_name=database.name,
        show_database_name=settings.show_database_name_on_track,
    )
    if not features:
        return (), f"{database.name}: no HMM hits found."
    report_e = _clamp_evalue(settings.hmm_evalue, default=DEFAULT_HMM_EVALUE)
    return (
        tuple(features),
        f"{database.name}: {len(features)} pyhmmer domain hit(s) (E-value ≤ {report_e:g}).",
    )


def _annotate_dna_hmm_with_pyhmmer(
    sequence: str,
    seq_id: str,
    database: UserDatabase,
    *,
    settings: SimilarityAnnotationSettings,
) -> tuple[tuple[SequenceFeature, ...], str]:
    hmm_path = database.path.expanduser().resolve()
    max_hits = _max_hits(settings.max_hits_per_database)
    with tempfile.TemporaryDirectory(prefix="repeatviewer_simann_nhmm_") as tmp:
        query_fa = Path(tmp) / "query.fa"
        SeqIO.write([SeqRecord(Seq(sequence), id=seq_id or "query")], str(query_fa), "fasta")
        raw_hits = collect_pyhmmer_domain_hits(
            hmm_path,
            query_fa,
            report_e=_clamp_evalue(settings.hmm_evalue, default=DEFAULT_HMM_EVALUE),
        )
    selected = _select_pyhmmer_hits_by_coverage(raw_hits, max_hits=max_hits)
    if not selected:
        report_e = _clamp_evalue(settings.hmm_evalue, default=DEFAULT_HMM_EVALUE)
        return (), f"{database.name}: no HMM hits (E-value ≤ {report_e:g})."
    features = _domain_features_from_nucleotide_hmm_hits(
        selected,
        database_name=database.name,
        show_database_name=settings.show_database_name_on_track,
    )
    if not features:
        return (), f"{database.name}: no HMM hits found."
    report_e = _clamp_evalue(settings.hmm_evalue, default=DEFAULT_HMM_EVALUE)
    return (
        tuple(features),
        f"{database.name}: {len(features)} pyhmmer domain hit(s) (E-value ≤ {report_e:g}).",
    )


def _annotate_protein_hmm_database(
    sequence: str,
    orfs: tuple[SequenceFeature, ...],
    database: UserDatabase,
    *,
    settings: SimilarityAnnotationSettings,
) -> tuple[tuple[SequenceFeature, ...], str]:
    try:
        import pyhmmer  # type: ignore[import-untyped]  # noqa: F401
    except ImportError:
        return (), f"{database.name}: pyhmmer is not installed ({pyhmmer_install_hint()})."

    try:
        return _annotate_protein_hmm_with_pyhmmer(sequence, orfs, database, settings=settings)
    except Exception as exc:  # noqa: BLE001
        return (), f"{database.name}: pyhmmer search failed ({exc})."


def _annotate_dna_hmm_database(
    sequence: str,
    seq_id: str,
    database: UserDatabase,
    *,
    settings: SimilarityAnnotationSettings,
) -> tuple[tuple[SequenceFeature, ...], str]:
    try:
        import pyhmmer  # type: ignore[import-untyped]  # noqa: F401
    except ImportError:
        return (), f"{database.name}: pyhmmer is not installed ({pyhmmer_install_hint()})."

    try:
        return _annotate_dna_hmm_with_pyhmmer(sequence, seq_id, database, settings=settings)
    except Exception as exc:  # noqa: BLE001
        return (), f"{database.name}: pyhmmer search failed ({exc})."


def _annotate_hmm_database(
    sequence: str,
    seq_id: str,
    orfs: tuple[SequenceFeature, ...],
    database: UserDatabase,
    *,
    settings: SimilarityAnnotationSettings,
) -> tuple[tuple[SequenceFeature, ...], str]:
    hmm_path = database.path.expanduser().resolve()
    readiness = database_search_readiness(database)
    if readiness is not None:
        return (), f"{database.name}: {readiness}"

    if hmm_uses_nucleotide_query(hmm_path):
        return _annotate_dna_hmm_database(sequence, seq_id, database, settings=settings)
    return _annotate_protein_hmm_database(sequence, orfs, database, settings=settings)


def _next_track(features: list[SequenceFeature]) -> int:
    if not features:
        return 0
    return max(feat.track for feat in features) + 1


def _domain_feature_signature(feat: SequenceFeature) -> tuple[int, int, str, str]:
    return (feat.start, feat.end, feat.name, feat.strand)


def _domain_features_for_color(
    features: list[SequenceFeature],
    color: str,
) -> list[SequenceFeature]:
    return [
        feat
        for feat in features
        if feat.feature_type == "domain" and feat.color == color
    ]


def _remove_domain_features_for_color(
    features: list[SequenceFeature],
    color: str,
) -> list[SequenceFeature]:
    return [
        feat
        for feat in features
        if not (feat.feature_type == "domain" and feat.color == color)
    ]


def _annotate_single_similarity_database(
    *,
    entry: SimilarityDatabaseRun,
    sequence: str,
    seq_id: str,
    orfs: tuple[SequenceFeature, ...],
    settings: SimilarityAnnotationSettings,
    use_system_blast: bool,
    blast_threads: int = 1,
    cancel: JobCancellation | None = None,
) -> tuple[tuple[SequenceFeature, ...], str]:
    db = entry.database
    if db.db_type == DATABASE_TYPE_CDD:
        return _annotate_cdd_database(
            sequence,
            seq_id,
            db,
            settings=settings,
            use_system_blast=use_system_blast,
            blast_threads=blast_threads,
            cancel=cancel,
        )
    if db.is_hmm:
        return _annotate_hmm_database(
            sequence,
            seq_id,
            orfs,
            db,
            settings=settings,
        )
    if db.db_type == DATABASE_TYPE_PROTEIN:
        return _annotate_protein_database(
            sequence,
            orfs,
            db,
            settings=settings,
            use_system_blast=use_system_blast,
            blast_threads=blast_threads,
            cancel=cancel,
        )
    if db.db_type == DATABASE_TYPE_DNA:
        return _annotate_dna_database(
            sequence,
            seq_id,
            db,
            settings=settings,
            use_system_blast=use_system_blast,
            blast_threads=blast_threads,
            cancel=cancel,
        )
    return (), f"{db.name}: unsupported database type {db.display_type_label!r}."


def _similarity_parallel_workers(
    database_count: int,
    *,
    batch_workers: int | None = None,
) -> int:
    workers = max(1, int(batch_workers if batch_workers is not None else load_batch_workers()))
    return max(1, min(workers, database_count))


def _similarity_blast_threads(
    database_count: int,
    *,
    batch_workers: int | None = None,
) -> int:
    parallel = _similarity_parallel_workers(database_count, batch_workers=batch_workers)
    return analysis_threads_per_job(batch_workers, parallel_jobs=parallel)


def _stamp_features(
    features: tuple[SequenceFeature, ...],
    *,
    track: int,
    color: str,
) -> list[SequenceFeature]:
    stamped: list[SequenceFeature] = []
    for feat in features:
        stamped.append(
            SequenceFeature(
                start=feat.start,
                end=feat.end,
                name=feat.name,
                feature_type=feat.feature_type,
                color=color,
                strand=feat.strand,
                track=track,
                parent=feat.parent,
                source=REPEATVIEWER_INTERNAL_SOURCE,
            )
        )
    return stamped


def run_similarity_annotations(
    fasta_path: Path,
    databases: tuple[SimilarityDatabaseRun, ...],
    *,
    settings: SimilarityAnnotationSettings | None = None,
    use_system_blast: bool = False,
    orf_settings: OrfSettings | None = None,
    cancel: JobCancellation | None = None,
    base_features: tuple[SequenceFeature, ...] | None = None,
    progress: Callable[[str], None] | None = None,
) -> SimilarityAnnotationResult:
    """Run similarity annotation for each selected user database."""
    if not databases:
        return SimilarityAnnotationResult((), "No databases selected.")

    cfg = settings or SimilarityAnnotationSettings()

    data = load_sequence_map_data(fasta_path)
    sequence = data.sequence
    seq_id = data.sequence_id
    features = list(base_features if base_features is not None else data.features)
    messages: list[str] = []
    track = _next_track(features)

    needs_orfs = similarity_search_needs_orfs(databases)
    orfs: tuple[SequenceFeature, ...] = ()
    if needs_orfs:
        try:
            orf_result = _resolve_orfs_for_similarity(
                sequence_id=seq_id,
                sequence=sequence,
                features=features,
                orf_settings=orf_settings,
            )
            orfs = orf_result.features
            if orf_result.message:
                messages.append(orf_result.message)
        except RuntimeError as exc:
            orfs = ()
            messages.append(str(exc))

    total = len(databases)
    workers = _similarity_parallel_workers(total)
    blast_threads = _similarity_blast_threads(total)
    db_results: list[tuple[SimilarityDatabaseRun, tuple[SequenceFeature, ...], str] | None] = [
        None
    ] * total
    completed = 0

    def _run_one(index: int, entry: SimilarityDatabaseRun) -> tuple[int, SimilarityDatabaseRun, tuple[SequenceFeature, ...], str]:
        check_cancelled(cancel)
        detected, msg = _annotate_single_similarity_database(
            entry=entry,
            sequence=sequence,
            seq_id=seq_id,
            orfs=orfs,
            settings=cfg,
            use_system_blast=use_system_blast,
            blast_threads=blast_threads,
            cancel=cancel,
        )
        return index, entry, detected, msg

    if progress is not None:
        progress(f"0/{total} starting")

    if workers > 1 and total > 1:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [
                pool.submit(_run_one, index, entry)
                for index, entry in enumerate(databases)
            ]
            for future in as_completed(futures):
                check_cancelled(cancel)
                index, entry, detected, msg = future.result()
                db_results[index] = (entry, detected, msg)
                completed += 1
                if progress is not None:
                    # Completed count after each DB finishes.
                    progress(f"{completed}/{total} {entry.database.name}")
    else:
        for index, entry in enumerate(databases):
            check_cancelled(cancel)
            if progress is not None:
                progress(f"{index}/{total} {entry.database.name}")
            _, _, detected, msg = _run_one(index, entry)
            db_results[index] = (entry, detected, msg)
            completed = index + 1
            if progress is not None:
                progress(f"{completed}/{total} {entry.database.name}")

    for index, result in enumerate(db_results):
        if result is None:
            continue
        entry, detected, msg = result
        db = entry.database
        if msg:
            messages.append(msg)
        if not detected:
            continue
        color = entry.color or color_for_database_index(index, database_name=db.name)
        existing = _domain_features_for_color(features, color)
        new_sigs = {_domain_feature_signature(feat) for feat in detected}
        old_sigs = {_domain_feature_signature(feat) for feat in existing}
        if existing and new_sigs == old_sigs:
            messages.append(f"{db.name}: unchanged ({len(detected)} hit(s)).")
            continue
        target_track = min(feat.track for feat in existing) if existing else track
        features = _remove_domain_features_for_color(features, color)
        features.extend(_stamp_features(detected, track=target_track, color=color))
        if not existing:
            track += 1

    summary = "\n".join(messages) if messages else "Similarity annotation complete."
    return SimilarityAnnotationResult(features=tuple(features), message=summary)


def load_use_system_blast_from_settings() -> bool:
    from PySide6.QtCore import QSettings

    s = QSettings()
    s.beginGroup("project")
    raw = s.value("use_system_blast", False)
    s.endGroup()
    if isinstance(raw, str):
        return raw.strip().lower() in {"1", "true", "yes", "on"}
    return bool(raw)
