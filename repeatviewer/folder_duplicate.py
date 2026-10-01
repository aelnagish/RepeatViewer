"""Duplicate folders: naming, size check and the copy itself (no Qt, unit tested)."""

from __future__ import annotations

import os
import shutil
from collections.abc import Callable, Iterable
from dataclasses import dataclass
from pathlib import Path

# Ask before copying more than this.
LARGE_COPY_BYTES = 1024**3  # 1 GB
# Stop counting after this many entries so a huge folder never freezes the UI.
MEASURE_ENTRY_LIMIT = 200_000


class CopyCancelled(Exception):
    """The user pressed Stop; the partial copy has been removed."""


@dataclass(frozen=True)
class FolderSize:
    bytes: int
    files: int
    complete: bool  # False: counting stopped early (limit reached)


def unique_folder_copy_path(folder: Path) -> Path:
    """Next unused sibling: ``name_copy``, ``name_copy2``, … (same rule as file Duplicate).

    The whole folder name is kept, dots included (``run2.fasta_out`` becomes
    ``run2.fasta_out_copy``).
    """
    folder = Path(folder)
    candidate = folder.parent / f"{folder.name}_copy"
    n = 2
    while candidate.exists() or candidate.is_symlink():
        candidate = folder.parent / f"{folder.name}_copy{n}"
        n += 1
    return candidate


def outermost_folders(folders: Iterable[Path]) -> list[Path]:
    """Drop folders that sit inside another selected folder (they are copied with it)."""
    resolved = []
    for folder in folders:
        path = Path(folder).expanduser().resolve()
        if path not in resolved:
            resolved.append(path)
    return [
        path
        for path in resolved
        if not any(other != path and other in path.parents for other in resolved)
    ]


def measure_folders(
    folders: Iterable[Path],
    *,
    stop_above_bytes: int | None = None,
    entry_limit: int = MEASURE_ENTRY_LIMIT,
) -> FolderSize:
    """Total size of the regular files under ``folders``. Symlinks are not followed.

    Counting stops early (``complete=False``) once ``stop_above_bytes`` is passed
    or ``entry_limit`` entries were seen.
    """
    total = 0
    files = 0
    entries = 0
    stack = [Path(folder) for folder in folders]
    while stack:
        current = stack.pop()
        try:
            with os.scandir(current) as it:
                for entry in it:
                    entries += 1
                    if entries > entry_limit:
                        return FolderSize(total, files, False)
                    try:
                        if entry.is_dir(follow_symlinks=False):
                            stack.append(Path(entry.path))
                        elif entry.is_file(follow_symlinks=False):
                            total += entry.stat(follow_symlinks=False).st_size
                            files += 1
                    except OSError:
                        continue
                    if stop_above_bytes is not None and total > stop_above_bytes:
                        return FolderSize(total, files, False)
        except OSError:
            continue
    return FolderSize(total, files, True)


def format_size(num_bytes: int) -> str:
    size = float(num_bytes)
    for unit in ("bytes", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            if unit == "bytes":
                return f"{int(size)} bytes"
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


def copy_folder(
    src: Path,
    dest: Path,
    *,
    progress: Callable[[int], None] | None = None,
    cancelled: Callable[[], bool] | None = None,
) -> None:
    """Copy ``src`` to the new folder ``dest`` (which must not exist yet).

    Symlinks are copied as links, file dates are kept. ``progress`` gets the
    size of every file after it is copied; ``cancelled`` is checked before each
    file. On Stop or on any error the partial ``dest`` is removed, so no
    half-finished copy is left behind.
    """
    src = Path(src)
    dest = Path(dest)
    if dest.exists() or dest.is_symlink():
        raise FileExistsError(f'"{dest.name}" already exists')

    def copy_one(source: str, target: str) -> str:
        if cancelled is not None and cancelled():
            raise CopyCancelled()
        result = shutil.copy2(source, target)
        if progress is not None:
            try:
                progress(os.path.getsize(target))
            except OSError:
                progress(0)
        return result

    try:
        # CopyCancelled is not an OSError, so copytree lets it through at once.
        shutil.copytree(src, dest, symlinks=True, copy_function=copy_one)
    except BaseException:
        shutil.rmtree(dest, ignore_errors=True)
        raise
