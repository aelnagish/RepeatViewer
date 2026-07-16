"""Metadata columns for the browser file table."""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from teviewer.backends.sequence_features import summarize_single_sequence_fasta

_FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".fna"}
_MAX_FASTA_SUMMARY_BYTES = 50 * 1024 * 1024


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
    """One pass: sequence count and first header description."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        return 0, "—"
    count = 0
    description = "—"
    try:
        if resolved.name.lower().endswith(".gz"):
            opener = gzip.open(resolved, "rt", encoding="utf-8", errors="replace")
        else:
            opener = resolved.open(encoding="utf-8", errors="replace")
        with opener as handle:
            for line in handle:
                if not line.startswith(">"):
                    continue
                count += 1
                if count == 1:
                    text = line[1:].strip()
                    description = (
                        text if len(text) <= max_desc_len else text[: max_desc_len - 1] + "…"
                    )
    except OSError:
        return 0, "—"
    return count, description


def enrich_fasta_metadata(path: Path) -> BrowserFileMetadata:
    """Full FASTA metadata (reads the file). Safe to call off the UI thread."""
    path = path.expanduser().resolve()
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

    return BrowserFileMetadata(
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


def metadata_for_file(path: Path) -> BrowserFileMetadata:
    """Full metadata (may read FASTA). Prefer quick + enrich for large folders."""
    path = path.expanduser().resolve()
    if is_fasta_path(path):
        return enrich_fasta_metadata(path)
    return quick_metadata_for_file(path)
