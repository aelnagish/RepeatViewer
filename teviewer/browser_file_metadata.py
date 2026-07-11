"""Metadata columns for the browser file table."""

from __future__ import annotations

import gzip
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path

from teviewer.backends.sequence_features import summarize_single_sequence_fasta

_FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".fna"}
_DNA_TYPES = _FASTA_SUFFIXES | {".seq"}
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


def _format_modified(path: Path) -> str:
    try:
        stamp = path.stat().st_mtime
    except OSError:
        return "—"
    dt = datetime.fromtimestamp(stamp)
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


def _is_fasta(path: Path) -> bool:
    name = path.name.lower()
    return path.suffix.lower() in _FASTA_SUFFIXES or name.endswith(".fa.gz")


def _fasta_header_count(path: Path) -> int:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        return 0
    count = 0
    try:
        if resolved.name.lower().endswith(".gz"):
            with gzip.open(resolved, "rt", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.startswith(">"):
                        count += 1
        else:
            with resolved.open(encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.startswith(">"):
                        count += 1
    except OSError:
        return 0
    return count


def _fasta_first_description(path: Path, *, max_len: int = 120) -> str:
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        return "—"
    try:
        if resolved.name.lower().endswith(".gz"):
            with gzip.open(resolved, "rt", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.startswith(">"):
                        text = line[1:].strip()
                        return text if len(text) <= max_len else text[: max_len - 1] + "…"
        else:
            with resolved.open(encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.startswith(">"):
                        text = line[1:].strip()
                        return text if len(text) <= max_len else text[: max_len - 1] + "…"
    except OSError:
        return "—"
    return "—"


def _molecule_type(path: Path) -> str:
    suffix = path.suffix.lower()
    name = path.name.lower()
    if _is_fasta(path):
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


def metadata_for_file(path: Path) -> BrowserFileMetadata:
    path = path.expanduser().resolve()
    modified = _format_modified(path)
    molecule = _molecule_type(path)
    try:
        mtime = path.stat().st_mtime
    except OSError:
        mtime = 0.0

    if _is_fasta(path):
        try:
            size = path.stat().st_size
        except OSError:
            size = 0
        seq_count = _fasta_header_count(path)
        seq_text = str(seq_count) if seq_count else "—"
        length_text = "—"
        gc_text = "—"
        length_key = -1.0
        gc_key = -1.0
        description = _fasta_first_description(path)
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
            modified_text=modified,
            description_text=description,
            length_text=length_text,
            size_text=_format_file_size(size),
            molecule_type=molecule,
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

    try:
        size = path.stat().st_size
        size_text = _format_file_size(size)
        size_key = float(size)
    except OSError:
        size_text = "—"
        size_key = -1.0

    return BrowserFileMetadata(
        modified_text=modified,
        description_text="—",
        length_text="—",
        size_text=size_text,
        molecule_type=molecule,
        sequence_count_text="—",
        gc_text="—",
        sort_name=path.name.lower(),
        sort_description="",
        sort_modified=mtime,
        sort_length=-1.0,
        sort_size=size_key,
        sort_sequence_count=-1,
        sort_gc=-1.0,
    )
