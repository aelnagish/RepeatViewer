"""Rename paths that contain spaces so CLI tools (BLAST, MAFFT, …) can run safely."""

from __future__ import annotations

from pathlib import Path


def basename_has_spaces(path: Path | str) -> bool:
    return " " in Path(path).name


def spaced_basename_target(path: Path) -> Path:
    """Return ``path`` with spaces in the basename replaced by underscores."""
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
