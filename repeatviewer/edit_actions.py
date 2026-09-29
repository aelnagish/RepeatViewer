"""Pure helpers for the Edit menu (no Qt, so they can be unit tested anywhere)."""

from __future__ import annotations

from typing import Any


def _count(n: int, noun: str) -> str:
    return f"{n} {noun}" if n == 1 else f"{n} {noun}s"


def browser_undo_label(entry: dict[str, Any] | None) -> str:
    """Menu text for Edit > Undo when the file or folder panel has focus.

    ``entry`` is the newest item of the main window's browser undo stack
    (``{"kind": "files" | "delete_folders" | "hide_folders", ...}``).
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
    return "Undo"
