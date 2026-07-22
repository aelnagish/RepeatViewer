"""Persisted file-list sort order for the current folder view."""

from __future__ import annotations

from enum import Enum
from pathlib import Path

from PySide6.QtCore import QSettings

FILES_SETTINGS_GROUP = "files"
DEFAULT_FILE_SORT_MODE = "name"


class FileSortMode(str, Enum):
    NAME = "name"
    TYPE = "type"
    TIME_DESC = "time_desc"
    TIME_ASC = "time_asc"


_FILE_SORT_LABELS: dict[FileSortMode, str] = {
    FileSortMode.NAME: "Name",
    FileSortMode.TYPE: "Type",
    FileSortMode.TIME_DESC: "Newest first",
    FileSortMode.TIME_ASC: "Oldest first",
}


def file_sort_mode_labels() -> list[tuple[FileSortMode, str]]:
    return [(mode, _FILE_SORT_LABELS[mode]) for mode in FileSortMode]


def load_file_sort_mode() -> FileSortMode:
    s = QSettings()
    s.beginGroup(FILES_SETTINGS_GROUP)
    raw = str(s.value("sort_mode", DEFAULT_FILE_SORT_MODE))
    s.endGroup()
    try:
        return FileSortMode(raw)
    except ValueError:
        return FileSortMode.NAME


def save_file_sort_mode(mode: FileSortMode) -> None:
    s = QSettings()
    s.beginGroup(FILES_SETTINGS_GROUP)
    s.setValue("sort_mode", mode.value)
    s.endGroup()


def _mtime(path: Path) -> float:
    try:
        return path.stat().st_mtime
    except OSError:
        return 0.0


def _type_group_key(path: Path) -> tuple[str, str]:
    """Group related extensions (e.g. .gff / .gff3 / .gtf) for Type sort."""
    name = path.name.lower()
    # Longer suffixes first so ``file.gff3`` is not treated as plain ``.gff``.
    grouped_suffixes: tuple[tuple[str, str], ...] = (
        ("annotation", ".gff3"),
        ("annotation", ".gff"),
        ("annotation", ".gtf"),
        ("annotation", ".bed"),
        ("annotation", ".vcf.gz"),
        ("annotation", ".vcf"),
        ("fasta", ".fasta"),
        ("fasta", ".fa"),
        ("fasta", ".fas"),
        ("fasta", ".fna"),
        ("fasta", ".faa"),
        ("html", ".html"),
        ("html", ".htm"),
        ("methylation", ".cov"),
        ("table", ".out"),
        ("table", ".tsv"),
        ("table", ".csv"),
        ("text", ".txt"),
        ("text", ".log"),
        ("text", ".md"),
        ("json", ".json"),
        ("xml", ".xml"),
        ("pdf", ".pdf"),
    )
    for group, suffix in grouped_suffixes:
        if name.endswith(suffix):
            return (group, name)
    suffix = path.suffix.lower() or "(no extension)"
    return (suffix, name)


def sort_file_paths(paths: list[Path], mode: FileSortMode) -> list[Path]:
    if mode == FileSortMode.NAME:
        return sorted(paths, key=lambda p: p.name.lower())
    if mode == FileSortMode.TYPE:
        return sorted(paths, key=_type_group_key)
    if mode == FileSortMode.TIME_DESC:
        return sorted(paths, key=_mtime, reverse=True)
    if mode == FileSortMode.TIME_ASC:
        return sorted(paths, key=_mtime)
    return sorted(paths, key=lambda p: p.name.lower())
