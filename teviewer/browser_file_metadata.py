"""Metadata columns for the browser file table."""

from __future__ import annotations

import gzip
import json
import threading
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path

from teviewer.app_preferences import config_dir
from teviewer.backends.sequence_features import summarize_single_sequence_fasta

_FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".fna"}
# Summarising a single-record FASTA loads the whole sequence into memory
# (Biopython parse + an uppercased copy) and scans it for GC content. Past this
# size that is slow enough to freeze the UI, so callers fall back to cheaper
# output. Shared with the preview panel, which performs the same load.
MAX_FASTA_SUMMARY_BYTES = 50 * 1024 * 1024
_MAX_FASTA_SUMMARY_BYTES = MAX_FASTA_SUMMARY_BYTES  # backwards-compatible alias
_FASTA_SCAN_BLOCK_BYTES = 1 << 20  # 1 MiB

_METADATA_CACHE: dict[str, tuple[float, int, "BrowserFileMetadata"]] = {}
_METADATA_CACHE_LOCK = threading.Lock()
_METADATA_CACHE_DIRTY = False
_METADATA_CACHE_SAVE_TIMER: threading.Timer | None = None


def _metadata_cache_path() -> Path:
    return config_dir() / "browser_file_metadata_cache.json"


def _load_metadata_cache_from_disk() -> None:
    path = _metadata_cache_path()
    if not path.is_file():
        return
    try:
        raw = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return
    if not isinstance(raw, dict):
        return
    loaded: dict[str, tuple[float, int, BrowserFileMetadata]] = {}
    for key, entry in raw.items():
        if not isinstance(entry, dict):
            continue
        try:
            mtime = float(entry["mtime"])
            size = int(entry["size"])
            meta = BrowserFileMetadata(**entry["meta"])
        except (KeyError, TypeError, ValueError):
            continue
        loaded[str(key)] = (mtime, size, meta)
    with _METADATA_CACHE_LOCK:
        _METADATA_CACHE.update(loaded)


def _schedule_metadata_cache_save() -> None:
    global _METADATA_CACHE_SAVE_TIMER
    with _METADATA_CACHE_LOCK:
        if _METADATA_CACHE_SAVE_TIMER is not None:
            _METADATA_CACHE_SAVE_TIMER.cancel()

        def _flush() -> None:
            global _METADATA_CACHE_SAVE_TIMER, _METADATA_CACHE_DIRTY
            with _METADATA_CACHE_LOCK:
                if not _METADATA_CACHE_DIRTY:
                    _METADATA_CACHE_SAVE_TIMER = None
                    return
                snapshot = dict(_METADATA_CACHE)
                _METADATA_CACHE_DIRTY = False
                _METADATA_CACHE_SAVE_TIMER = None
            path = _metadata_cache_path()
            try:
                path.parent.mkdir(parents=True, exist_ok=True)
                payload = {
                    key: {
                        "mtime": mtime,
                        "size": size,
                        "meta": asdict(meta),
                    }
                    for key, (mtime, size, meta) in snapshot.items()
                }
                path.write_text(json.dumps(payload), encoding="utf-8")
            except OSError:
                pass

        _METADATA_CACHE_SAVE_TIMER = threading.Timer(2.0, _flush)
        _METADATA_CACHE_SAVE_TIMER.daemon = True
        _METADATA_CACHE_SAVE_TIMER.start()


def _cache_get(path: Path) -> BrowserFileMetadata | None:
    key = str(path.expanduser().resolve())
    with _METADATA_CACHE_LOCK:
        entry = _METADATA_CACHE.get(key)
    if entry is None:
        return None
    mtime, size, meta = entry
    try:
        st = path.stat()
    except OSError:
        return None
    if float(st.st_mtime) != mtime or int(st.st_size) != size:
        return None
    return meta


def _cache_put(path: Path, meta: BrowserFileMetadata) -> None:
    global _METADATA_CACHE_DIRTY
    key = str(path.expanduser().resolve())
    try:
        st = path.stat()
        mtime = float(st.st_mtime)
        size = int(st.st_size)
    except OSError:
        return
    with _METADATA_CACHE_LOCK:
        _METADATA_CACHE[key] = (mtime, size, meta)
        _METADATA_CACHE_DIRTY = True
    _schedule_metadata_cache_save()


def _iter_fasta_blocks(path: Path, block_size: int = _FASTA_SCAN_BLOCK_BYTES):
    """Yield raw byte blocks of a FASTA file, decompressing .gz transparently."""
    opener = gzip.open(path, "rb") if path.name.lower().endswith(".gz") else path.open("rb")
    with opener as handle:
        while True:
            block = handle.read(block_size)
            if not block:
                return
            yield block


def scan_fasta_headers(
    path: Path,
    *,
    stop_at: int | None = None,
    want_description: bool = False,
    max_desc_len: int = 120,
) -> tuple[int, str]:
    """
    Count FASTA header lines (and optionally read the first description).

    Scans raw bytes and counts ``b"\\n>"`` with ``bytes.count`` rather than
    iterating decoded text lines. That avoids UTF-8 decoding the whole file and
    avoids a Python-level loop per line -- on a 500 MB FASTA it is the
    difference between tens of seconds and well under one.

    ``stop_at`` stops the scan once that many headers have been seen, and the
    returned count then *saturates* at ``stop_at`` -- it is a lower bound, only
    meaningful compared against ``stop_at``. Callers that only need "one
    sequence or many?" should pass ``stop_at=2``; callers that *display* the
    real total (the preview info bar, the Sequences column) must leave it
    ``None``.
    """
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        return 0, "—"

    count = 0
    description = "—"
    tail = b""  # last byte of the previous block: catches a "\n" / ">" split
    first_block = True

    try:
        for block in _iter_fasta_blocks(resolved):
            if first_block:
                first_block = False
                if want_description:
                    start = 0 if block.startswith(b">") else -1
                    if start == -1:
                        found = block.find(b"\n>")
                        start = found + 1 if found != -1 else -1
                    if start != -1:
                        end = block.find(b"\n", start)
                        raw = block[start + 1 : end if end != -1 else len(block)]
                        text = raw.decode("utf-8", errors="replace").strip()
                        if text:
                            description = (
                                text
                                if len(text) <= max_desc_len
                                else text[: max_desc_len - 1] + "…"
                            )
                if block.startswith(b">"):
                    count += 1
            elif tail == b"\n" and block.startswith(b">"):
                count += 1

            count += block.count(b"\n>")
            if stop_at is not None and count >= stop_at:
                return min(count, stop_at), description
            tail = block[-1:]
    except OSError:
        return 0, "—"

    return count, description


@dataclass(frozen=True)
class BrowserFileMetadata:
    modified_text: str
    description_text: str
    length_text: str
    size_text: str
    molecule_type: str
    sequence_count_text: str
    gc_text: str
    sort_name: str
    sort_description: str
    sort_modified: float
    sort_length: float
    sort_size: float
    sort_sequence_count: int
    sort_gc: float


def _format_modified(mtime: float) -> str:
    if mtime <= 0:
        return "—"
    dt = datetime.fromtimestamp(mtime)
    time_text = dt.strftime("%I:%M %p").lstrip("0")
    return f"{dt.strftime('%d %b %Y')} {time_text}"


def _format_file_size(num_bytes: int) -> str:
    if num_bytes < 1024:
        return f"{num_bytes} B"
    if num_bytes < 1024 * 1024:
        return f"{num_bytes / 1024:.1f} KB"
    if num_bytes < 1024 * 1024 * 1024:
        return f"{num_bytes / (1024 * 1024):.1f} MB"
    return f"{num_bytes / (1024 * 1024 * 1024):.1f} GB"


def is_fasta_path(path: Path) -> bool:
    name = path.name.lower()
    return path.suffix.lower() in _FASTA_SUFFIXES or name.endswith(".fa.gz")


def _molecule_type(path: Path) -> str:
    suffix = path.suffix.lower()
    name = path.name.lower()
    if is_fasta_path(path):
        return "DNA"
    if suffix in {".gff", ".gff3", ".gtf", ".bed"}:
        return "Annotation"
    if suffix in {".html", ".htm"}:
        return "HTML"
    if suffix in {".png", ".jpg", ".jpeg", ".gif", ".svg", ".pdf"}:
        return "Image" if suffix != ".pdf" else "PDF"
    if suffix in {".bam", ".sam", ".cram"}:
        return "Alignment"
    if suffix in {".out"}:
        return "Table"
    if suffix in {".txt", ".log", ".md"} or name.endswith("_blast.txt"):
        return "Text"
    if suffix in {".vcf"} or name.endswith(".vcf.gz"):
        return "Variant"
    if suffix:
        return suffix.lstrip(".").upper()
    return "File"


def minimal_name_metadata(path: Path) -> BrowserFileMetadata:
    """Cheapest row metadata: filename only (no stat or file reads)."""
    resolved = path.expanduser()
    try:
        resolved = resolved.resolve()
    except OSError:
        pass
    name = resolved.name
    return BrowserFileMetadata(
        modified_text="—",
        description_text="—",
        length_text="—",
        size_text="—",
        molecule_type="—",
        sequence_count_text="—",
        gc_text="—",
        sort_name=name.lower(),
        sort_description="",
        sort_modified=0.0,
        sort_length=-1.0,
        sort_size=0.0,
        sort_sequence_count=-1,
        sort_gc=-1.0,
    )


def quick_metadata_for_file(path: Path) -> BrowserFileMetadata:
    """Cheap metadata from ``stat`` + suffix only (no file reads).

    Used to populate large folders quickly; FASTA columns are filled later by
    :func:`enrich_fasta_metadata`.
    """
    resolved = path.expanduser()
    try:
        resolved = resolved.resolve()
    except OSError:
        pass
    try:
        st = resolved.stat()
        size = int(st.st_size)
        mtime = float(st.st_mtime)
    except OSError:
        size = 0
        mtime = 0.0

    molecule = _molecule_type(resolved)
    pending = is_fasta_path(resolved)
    return BrowserFileMetadata(
        modified_text=_format_modified(mtime),
        description_text="…" if pending else "—",
        length_text="…" if pending else "—",
        size_text=_format_file_size(size),
        molecule_type=molecule,
        sequence_count_text="…" if pending else "—",
        gc_text="…" if pending else "—",
        sort_name=resolved.name.lower(),
        sort_description="",
        sort_modified=mtime,
        sort_length=-1.0,
        sort_size=float(size),
        sort_sequence_count=-1 if pending else -1,
        sort_gc=-1.0,
    )


def _fasta_scan(path: Path, *, max_desc_len: int = 120) -> tuple[int, str]:
    """One pass: sequence count and first header description.

    The Sequences column shows the real total, so this must not use ``stop_at``.
    """
    return scan_fasta_headers(path, want_description=True, max_desc_len=max_desc_len)


def enrich_fasta_metadata(path: Path) -> BrowserFileMetadata:
    """Full FASTA metadata (reads the file). Safe to call off the UI thread."""
    path = path.expanduser().resolve()
    cached = _cache_get(path)
    if cached is not None:
        return cached

    try:
        st = path.stat()
        size = int(st.st_size)
        mtime = float(st.st_mtime)
    except OSError:
        size = 0
        mtime = 0.0

    seq_count, description = _fasta_scan(path)
    seq_text = str(seq_count) if seq_count else "—"
    length_text = "—"
    gc_text = "—"
    length_key = -1.0
    gc_key = -1.0
    if seq_count == 1 and size <= _MAX_FASTA_SUMMARY_BYTES:
        try:
            summary = summarize_single_sequence_fasta(path)
            length_text = f"{summary.length:,}"
            length_key = float(summary.length)
            if summary.gc_percent is not None:
                gc_text = f"{summary.gc_percent:.1f}"
                gc_key = float(summary.gc_percent)
        except (OSError, ValueError):
            pass

    meta = BrowserFileMetadata(
        modified_text=_format_modified(mtime),
        description_text=description,
        length_text=length_text,
        size_text=_format_file_size(size),
        molecule_type="DNA",
        sequence_count_text=seq_text,
        gc_text=gc_text,
        sort_name=path.name.lower(),
        sort_description=description.lower(),
        sort_modified=mtime,
        sort_length=length_key,
        sort_size=float(size),
        sort_sequence_count=seq_count,
        sort_gc=gc_key,
    )
    _cache_put(path, meta)
    return meta


def metadata_for_file(path: Path) -> BrowserFileMetadata:
    """Full metadata (may read FASTA). Prefer quick + enrich for large folders."""
    path = path.expanduser().resolve()
    if is_fasta_path(path):
        return enrich_fasta_metadata(path)
    return quick_metadata_for_file(path)


_load_metadata_cache_from_disk()
