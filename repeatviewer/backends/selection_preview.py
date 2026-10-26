"""Temporary combined FASTA for previewing several selected files together.

When two or more FASTA files are ticked (or highlighted) in the Files panel,
the Preview shows them as one multi-sequence file. The combined copy lives in
the app cache, never next to the user's files, and is rebuilt only when a
source file or its annotation sidecar changes.
"""

from __future__ import annotations

import hashlib
import json
import shutil
from pathlib import Path

SELECTION_PREVIEW_MIN_FILES = 2
SELECTION_PREVIEW_MAX_FILES = 100
SELECTION_PREVIEW_MAX_BYTES = 200 * 1024 * 1024
# Newest combined copies kept on disk (at most about 1 GB at the size limit).
_KEEP_CACHED = 5
_MANIFEST_NAME = "selection_preview.json"


def selection_preview_root() -> Path:
    from repeatviewer.app_preferences import cache_dir

    return cache_dir() / "selection_preview"


def _manifest_for(path: Path) -> Path:
    return path.parent / _MANIFEST_NAME


def is_selection_preview_path(path: object) -> bool:
    """True when ``path`` is a combined file built by :func:`build_selection_preview`."""
    if path is None:
        return False
    try:
        resolved = Path(str(path)).expanduser().resolve()
        root = selection_preview_root().expanduser().resolve()
    except (OSError, TypeError, ValueError):
        return False
    return resolved.parent.parent == root and _manifest_for(resolved).is_file()


def selection_preview_sources(path: Path) -> list[Path]:
    """Source FASTA files of a combined preview, in display order."""
    try:
        data = json.loads(_manifest_for(Path(path)).read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return []
    return [Path(p) for p in data.get("sources", []) if str(p).strip()]


def selection_preview_limit_message(files: list[Path]) -> str | None:
    """Why ``files`` cannot be previewed together, or ``None`` when they can."""
    if len(files) > SELECTION_PREVIEW_MAX_FILES:
        return (
            f"{len(files)} files are selected. The Preview shows at most "
            f"{SELECTION_PREVIEW_MAX_FILES} files together."
        )
    total = 0
    for path in files:
        try:
            total += Path(path).stat().st_size
        except OSError:
            continue
    if total > SELECTION_PREVIEW_MAX_BYTES:
        limit_mb = SELECTION_PREVIEW_MAX_BYTES // (1024 * 1024)
        return (
            f"The selected files are too large to preview together "
            f"({total / (1024 * 1024):.0f} MB, limit {limit_mb} MB)."
        )
    return None


def _fingerprint(files: list[Path]) -> str:
    from repeatviewer.backends.sequence_features import resolve_sidecar_path

    digest = hashlib.sha1()
    for path in files:
        for item in (path, resolve_sidecar_path(path)):
            if item is None:
                digest.update(b"-\n")
                continue
            try:
                st = item.stat()
            except OSError:
                digest.update(f"{item}|missing\n".encode())
                continue
            digest.update(f"{item}|{st.st_mtime_ns}|{st.st_size}\n".encode())
    return digest.hexdigest()[:16]


def _prune_cache(root: Path, keep: Path) -> None:
    try:
        folders = [p for p in root.iterdir() if p.is_dir() and p != keep]
    except OSError:
        return
    folders.sort(key=lambda p: p.stat().st_mtime, reverse=True)
    for old in folders[_KEEP_CACHED - 1:]:
        shutil.rmtree(old, ignore_errors=True)


def build_selection_preview(files: list[Path]) -> Path:
    """Write (or reuse) the combined FASTA for ``files`` and return its path.

    Records keep the order of ``files``. Annotations from each source sidecar
    are copied in the same way as Combine Sequences to One File.
    """
    from repeatviewer.backends.sequence_compose import (
        combine_annotation_json,
        group_sequences,
    )

    sources = [Path(p).expanduser().resolve() for p in files]
    if len(sources) < SELECTION_PREVIEW_MIN_FILES:
        raise ValueError("Select at least two FASTA files to preview together.")
    root = selection_preview_root()
    folder = root / _fingerprint(sources)
    out_path = folder / f"selection_{len(sources)}_files.fa"
    manifest = folder / _MANIFEST_NAME
    if out_path.is_file() and manifest.is_file():
        return out_path
    if folder.exists():
        shutil.rmtree(folder, ignore_errors=True)
    folder.mkdir(parents=True, exist_ok=True)
    group_sequences(sources, out_path)
    try:
        combine_annotation_json(sources, out_path)
    except Exception:  # noqa: BLE001 - annotations are optional for a preview
        pass
    manifest.write_text(
        json.dumps({"sources": [str(p) for p in sources]}, indent=1),
        encoding="utf-8",
    )
    _prune_cache(root, keep=folder)
    return out_path
