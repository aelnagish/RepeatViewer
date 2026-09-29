"""Pure helpers for the Edit menu (no Qt, so they can be unit tested anywhere)."""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
from typing import Any

# Finder writes these into any folder it opens; they do not count as a change.
_IGNORED_NAMES = {".DS_Store"}


def _count(n: int, noun: str) -> str:
    return f"{n} {noun}" if n == 1 else f"{n} {noun}s"


def browser_undo_label(entry: dict[str, Any] | None) -> str:
    """Menu text for Edit > Undo when the file or folder panel has focus.

    ``entry`` is the newest item of the main window's browser undo stack.
    """
    if not entry:
        return "Undo"
    kind = str(entry.get("kind") or "")
    if kind == "files":
        return f"Undo Delete {_count(len(entry.get('items') or []), 'File')}"
    if kind == "delete_folders":
        return f"Undo Delete {_count(len(entry.get('items') or []), 'Folder')}"
    if kind == "hide_folders":
        n = len(entry.get("folders") or [])
        return f"Undo Remove {_count(n, 'Folder')} from List"
    if kind == "transfer":
        action = str(entry.get("action") or "Copy")
        return f"Undo {action} {_count(int(entry.get('count') or 0), 'File')}"
    if kind == "created_folders":
        return f"Undo Duplicate {_count(len(entry.get('items') or []), 'Folder')}"
    if kind == "rename_file":
        return "Undo Rename File"
    if kind == "rename_folder":
        return "Undo Rename Folder"
    return "Undo"


def path_fingerprint(path: Path) -> str | None:
    """Content signature of a file or folder: names, sizes and modification times.

    Undo only removes a copy it made when the fingerprint still matches, so a
    copy you changed afterwards is never deleted. None if the path is missing.
    """
    path = Path(path)
    digest = hashlib.sha1()
    try:
        if path.is_file() and not path.is_symlink():
            st = path.stat()
            digest.update(f"file\0{st.st_size}\0{st.st_mtime_ns}".encode())
            return digest.hexdigest()
        if not path.is_dir() or path.is_symlink():
            return None
    except OSError:
        return None
    entries: list[str] = []
    for root, dirs, files in os.walk(path):
        dirs.sort()
        rel_root = os.path.relpath(root, path)
        for name in dirs:
            entries.append(f"d\0{os.path.join(rel_root, name)}")
        for name in files:
            if name in _IGNORED_NAMES:
                continue
            full = os.path.join(root, name)
            try:
                st = os.lstat(full)
            except OSError:
                continue
            entries.append(f"f\0{os.path.join(rel_root, name)}\0{st.st_size}\0{st.st_mtime_ns}")
    for line in sorted(entries):
        digest.update(line.encode("utf-8", "surrogateescape"))
        digest.update(b"\n")
    return digest.hexdigest()
