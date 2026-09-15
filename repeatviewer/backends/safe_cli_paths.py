"""Rename paths that contain spaces so CLI tools (BLAST, MAFFT, …) can run safely."""

from __future__ import annotations

import os
from dataclasses import dataclass
from pathlib import Path


def basename_has_spaces(path: Path | str) -> bool:
    return " " in Path(path).name


def spaced_basename_target(path: Path) -> Path:
    """Return ``path`` with each space in the basename replaced by ``_``.

    One space → ``_``, two consecutive spaces → ``__``, and so on.
    """
    path = path.expanduser().resolve()
    safe_name = path.name.replace(" ", "_")
    if safe_name == path.name:
        return path
    candidate = path.with_name(safe_name)
    if not candidate.exists():
        return candidate
    stem = Path(safe_name).stem
    suffix = Path(safe_name).suffix
    index = 1
    while True:
        candidate = path.with_name(f"{stem}_{index}{suffix}")
        if not candidate.exists():
            return candidate
        index += 1


def rename_spaces_in_basename(path: Path) -> tuple[Path, bool]:
    """Rename ``path`` on disk when its basename contains spaces.

    Returns ``(final_path, renamed)``. If no rename is needed, returns the
    resolved original path and ``False``.
    """
    path = path.expanduser().resolve()
    if not basename_has_spaces(path):
        return path, False
    if not path.exists():
        raise FileNotFoundError(path)
    target = spaced_basename_target(path)
    path.rename(target)
    return target.resolve(), True


@dataclass(frozen=True)
class SpacedNameHit:
    """A file or folder whose basename contains one or more spaces."""

    path: Path
    kind: str  # "folder" | "file"
    display: str


def find_names_with_spaces(
    root: Path | str,
    *,
    limit: int = 500,
) -> tuple[list[SpacedNameHit], bool]:
    """Find folders and files with spaces in their basename under ``root``.

    Also includes ancestor directories of ``root`` whose names contain spaces
    (needed when the import path itself has spaces).
    """
    root = Path(root).expanduser().resolve()
    hits: list[SpacedNameHit] = []
    truncated = False

    # Build ancestor paths part by part (e.g. ``…/TEtrimmer output/…``).
    parts = root.parts
    if parts:
        current = Path(parts[0])
        for part in parts[1:]:
            current = current / part
            if " " not in part:
                continue
            if len(hits) >= limit:
                truncated = True
                break
            display = (
                f"[folder] {part}  (in import path)"
                if current != root
                else f"[folder] {part}"
            )
            hits.append(SpacedNameHit(path=current, kind="folder", display=display))

    if truncated or not root.is_dir():
        return hits, truncated

    root_s = str(root)
    for dirpath, dirnames, filenames in os.walk(root_s, followlinks=False):
        dirnames.sort()
        filenames.sort()
        parent = Path(dirpath)

        for name in dirnames:
            if " " not in name:
                continue
            if len(hits) >= limit:
                truncated = True
                dirnames.clear()
                break
            child = parent / name
            # Skip if this is ``root`` itself (already listed via ancestors).
            if child.resolve() == root and any(h.path == root for h in hits):
                continue
            try:
                rel = child.relative_to(root)
                display = f"[folder] {rel}"
            except ValueError:
                display = f"[folder] {child}"
            hits.append(SpacedNameHit(path=child, kind="folder", display=display))

        if truncated:
            break

        for name in filenames:
            if " " not in name:
                continue
            if len(hits) >= limit:
                truncated = True
                dirnames.clear()
                break
            child = parent / name
            try:
                rel = child.relative_to(root)
                display = f"[file] {rel}"
            except ValueError:
                display = f"[file] {child}"
            hits.append(SpacedNameHit(path=child, kind="file", display=display))

        if truncated:
            break

    return hits, truncated


def find_folders_with_spaces(
    root: Path | str,
    *,
    limit: int = 80,
) -> tuple[list[str], bool]:
    """Backward-compatible wrapper: folder display strings only."""
    hits, truncated = find_names_with_spaces(root, limit=limit)
    entries = [h.display.removeprefix("[folder] ").strip() for h in hits if h.kind == "folder"]
    return entries, truncated


def apply_path_renames(path: Path, renames: list[tuple[Path, Path]]) -> Path:
    """Rewrite ``path`` after a sequence of directory/file renames (old → new)."""
    result = path.expanduser().resolve()
    for old, new in renames:
        if result == old:
            result = new
            continue
        try:
            result = new / result.relative_to(old)
        except ValueError:
            continue
    return result


def rename_spaced_paths(paths: list[Path]) -> tuple[list[tuple[Path, Path]], list[str]]:
    """Rename each path's basename (spaces → ``_``), deepest paths first.

    Returns ``(renames, errors)`` where ``renames`` is ``[(old, new), ...]``.
    """
    unique: dict[str, Path] = {}
    for raw in paths:
        try:
            path = Path(raw).expanduser().resolve()
        except OSError:
            continue
        if basename_has_spaces(path) and path.exists():
            unique[str(path)] = path

    ordered = sorted(unique.values(), key=lambda p: len(p.parts), reverse=True)
    renames: list[tuple[Path, Path]] = []
    errors: list[str] = []
    for path in ordered:
        # Parent may already have been renamed; resolve via prior renames.
        live = apply_path_renames(path, renames)
        if not live.exists():
            errors.append(f"Missing: {path}")
            continue
        if not basename_has_spaces(live):
            continue
        try:
            new_path, renamed = rename_spaces_in_basename(live)
        except OSError as exc:
            errors.append(f"{live.name}: {exc}")
            continue
        if renamed:
            renames.append((live, new_path))
    return renames, errors
