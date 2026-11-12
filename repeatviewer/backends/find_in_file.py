"""Search one FASTA file for Edit > Find: sequence motif, sequence names, annotation names.

No Qt here (unit tested). Big files are searched in a worker thread: the file is
streamed record by record, ``progress`` is called after each sequence (rule 13:
"N done, M left") and ``cancelled`` is checked so Stop works.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import gzip
from collections.abc import Callable, Iterator
from dataclasses import dataclass, field
from pathlib import Path

from repeatviewer.backends.motif_search import (
    FindCancelled,
    FindOptions,
    clean_query,
    find_in_sequence,
    query_problem,
)


@dataclass(frozen=True)
class FindRequest:
    query: str
    options: FindOptions = FindOptions()
    in_sequence: bool = True
    in_names: bool = True
    in_annotations: bool = True
    case_sensitive: bool = False  # sequence names and annotation names only


@dataclass(frozen=True)
class FoundItem:
    seq_id: str
    start: int  # 1-based, inclusive
    end: int
    strand: str = "+"
    kind: str = "sequence"  # "sequence" | "name" | "annotation"
    mismatches: int = 0
    frame: int = 0
    label: str = ""  # the annotation name for kind == "annotation"


@dataclass
class FileFindResult:
    items: list[FoundItem] = field(default_factory=list)
    truncated: bool = False
    record_ids: list[str] = field(default_factory=list)
    lengths: dict[str, int] = field(default_factory=dict)
    sequence_skipped: str | None = None  # why the motif search did not run


def _open_text(path: Path):
    if path.name.lower().endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return open(path, encoding="utf-8", errors="replace")


def iter_fasta_records(
    path: Path, *, cancelled: Callable[[], bool] | None = None
) -> Iterator[tuple[str, str]]:
    """(id, sequence) for each record, streamed.

    Ids are named exactly like the Annotation views name them: the first word of
    the header, duplicates renamed ``name#2``, ``name#3``, ... (see
    sequence_features.fasta_record_lengths / uniquify_record_ids).
    """
    seen: dict[str, int] = {}
    header: str | None = None
    chunks: list[str] = []

    def emit(raw_header: str) -> tuple[str, str] | None:
        words = raw_header[1:].strip().split()
        if not words:
            return None
        base = words[0]
        count = seen.get(base, 0) + 1
        seen[base] = count
        return (base if count == 1 else f"{base}#{count}"), "".join(chunks)

    with _open_text(Path(path)) as handle:
        for line_no, line in enumerate(handle):
            if cancelled is not None and line_no % 200_000 == 0 and cancelled():
                raise FindCancelled()
            if line.startswith(">"):
                if header is not None:
                    record = emit(header)
                    if record is not None:
                        yield record
                header = line.rstrip("\r\n")
                chunks = []
            elif header is not None:
                chunks.append(line.strip())
        if header is not None:
            record = emit(header)
            if record is not None:
                yield record


def text_matches(text: str, needle: str, case_sensitive: bool) -> bool:
    if not needle:
        return False
    return needle in text if case_sensitive else needle.lower() in text.lower()


def _load_saved_features(path: Path, record_ids: list[str], lengths: dict[str, int]):
    """The file's saved annotations, as the multi-FASTA overview loads them."""
    from repeatviewer.backends.sequence_features import (
        load_all_features_by_id,
        load_features_by_id_from_json,
        resolve_features_gff_path,
    )

    json_by_id = load_features_by_id_from_json(path, record_ids)
    if json_by_id and any(json_by_id.values()):
        return json_by_id
    gff = resolve_features_gff_path(path)
    if gff is not None:
        return load_all_features_by_id(gff, record_ids, lengths_by_id=lengths)
    return {}


def find_in_fasta(
    path: Path,
    request: FindRequest,
    *,
    features_by_id: dict | None = None,
    cancelled: Callable[[], bool] | None = None,
    progress: Callable[[int, str], None] | None = None,
) -> FileFindResult:
    """Everything ``request`` asks for in one FASTA file, in file order.

    ``features_by_id`` (id -> list of SequenceFeature-like objects with start,
    end, name, feature_type, strand) lets the caller pass annotations that are
    open in the viewer; otherwise the saved ones are read from disk.
    """
    path = Path(path).expanduser().resolve()
    result = FileFindResult()
    needle = request.query.strip()
    options = request.options
    motif_ok = False
    if request.in_sequence:
        problem = query_problem(needle, options)
        if problem is None:
            motif_ok = True
        else:
            result.sequence_skipped = problem
    limit = options.max_hits
    items: list[FoundItem] = []

    for done, (seq_id, sequence) in enumerate(iter_fasta_records(path, cancelled=cancelled), 1):
        if cancelled is not None and cancelled():
            raise FindCancelled()
        result.record_ids.append(seq_id)
        result.lengths[seq_id] = len(sequence)
        if request.in_names and len(items) < limit and text_matches(
            seq_id, needle, request.case_sensitive
        ):
            items.append(FoundItem(seq_id, 1, max(1, len(sequence)), "+", "name"))
        if motif_ok and len(items) < limit and sequence:
            found = find_in_sequence(
                sequence,
                clean_query(needle),
                options,
                seq_id=seq_id,
                cancelled=cancelled,
                limit=limit - len(items),
            )
            result.truncated = result.truncated or found.truncated
            items.extend(
                FoundItem(h.seq_id, h.start, h.end, h.strand, "sequence", h.mismatches, h.frame)
                for h in found.hits
            )
        if progress is not None:
            progress(done, seq_id)
        if len(items) >= limit:
            result.truncated = True

    if request.in_annotations and needle and len(items) < limit:
        features = features_by_id
        if features is None:
            try:
                features = _load_saved_features(path, result.record_ids, result.lengths)
            except (OSError, ValueError):
                features = {}
        for seq_id in result.record_ids:
            for feat in features.get(seq_id) or []:
                name = str(getattr(feat, "name", "") or getattr(feat, "feature_type", ""))
                if not (
                    text_matches(name, needle, request.case_sensitive)
                    or text_matches(str(getattr(feat, "feature_type", "")), needle, request.case_sensitive)
                ):
                    continue
                if len(items) >= limit:
                    result.truncated = True
                    break
                items.append(
                    FoundItem(
                        seq_id,
                        int(feat.start),
                        int(feat.end),
                        str(getattr(feat, "strand", "+") or "+"),
                        "annotation",
                        label=name,
                    )
                )

    order = {seq_id: i for i, seq_id in enumerate(result.record_ids)}
    kind_rank = {"name": 0, "annotation": 1, "sequence": 2}
    items.sort(key=lambda it: (order.get(it.seq_id, 0), it.start, kind_rank.get(it.kind, 3), it.end))
    result.items = items[:limit]
    return result


def next_index(items: list[FoundItem], current: int, direction: int) -> int:
    """Index of the next (+1) or previous (-1) item, wrapping around."""
    if not items:
        return -1
    if current < 0:
        return 0 if direction >= 0 else len(items) - 1
    return (current + (1 if direction >= 0 else -1)) % len(items)


def first_index_after(items: list[FoundItem], seq_id: str, position: int, record_ids: list[str]) -> int:
    """First item at or after ``position`` in ``seq_id`` (in file order), else 0."""
    if not items:
        return -1
    order = {rid: i for i, rid in enumerate(record_ids)}
    here = (order.get(seq_id, 0), position)
    for index, item in enumerate(items):
        if (order.get(item.seq_id, 0), item.start) >= here:
            return index
    return 0
