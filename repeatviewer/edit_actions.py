"""Pure helpers for the Edit menu (no Qt, so they can be unit tested anywhere).

Developer: Abdullah Yusuf
"""

from __future__ import annotations

import hashlib
import os
from pathlib import Path
from typing import Any

# Finder writes these into any folder it opens; they do not count as a change.
_IGNORED_NAMES = {".DS_Store"}


def _count(n: int, noun: str) -> str:
    return f"{n} {noun}" if n == 1 else f"{n} {noun}s"


def _action_phrase(kind: str, count: int, action: str = "") -> str:
    """What a browser undo/redo step did, e.g. "Delete 3 Files" ("" if unknown)."""
    if kind == "files":
        return f"Delete {_count(count, 'File')}"
    if kind == "delete_folders":
        return f"Delete {_count(count, 'Folder')}"
    if kind == "hide_folders":
        return f"Remove {_count(count, 'Folder')} from List"
    if kind == "transfer":
        return f"{action or 'Copy'} {_count(count, 'File')}"
    if kind == "created_folders":
        return f"Duplicate {_count(count, 'Folder')}"
    if kind == "rename_file":
        return "Rename File"
    if kind == "rename_folder":
        return "Rename Folder"
    if kind == "batch_rename":
        return f"Rename {_count(count, 'File')}"
    return ""


def browser_undo_label(entry: dict[str, Any] | None) -> str:
    """Menu text for Edit > Undo when the file or folder panel has focus.

    ``entry`` is the newest item of the main window's browser undo stack.
    """
    if not entry:
        return "Undo"
    kind = str(entry.get("kind") or "")
    if kind in {"files", "delete_folders", "created_folders"}:
        count = len(entry.get("items") or [])
    elif kind == "hide_folders":
        count = len(entry.get("folders") or [])
    else:
        count = int(entry.get("count") or 0)
    phrase = _action_phrase(kind, count, str(entry.get("action") or ""))
    return f"Undo {phrase}" if phrase else "Undo"


def browser_redo_label(entry: dict[str, Any] | None) -> str:
    """Menu text for Edit > Redo; ``entry`` is the newest item of the redo stack."""
    if not entry:
        return "Redo"
    phrase = _action_phrase(
        str(entry.get("kind") or ""),
        int(entry.get("count") or 0),
        str(entry.get("action") or ""),
    )
    return f"Redo {phrase}" if phrase else "Redo"


# IUPAC complement; U pairs with A. Gaps and N stay as they are.
_DNA_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "TGCAAYRSWMKVHDBNtgcaayrswmkvhdbn",
)
_RNA_COMPLEMENT = str.maketrans(
    "ACGTURYSWKMBDHVNacgturyswkmbdhvn",
    "UGCAAYRSWMKVHDBNugcaayrswmkvhdbn",
)
_SEQUENCE_LETTERS = frozenset("ACGTURYSWKMBDHVNacgturyswkmbdhvn-.")
_IGNORED_IN_SEQUENCE = frozenset(" \t\r\n0123456789")


def _reverse_complement_block(lines: list[str]) -> tuple[list[str], int] | None:
    """Sequence lines -> reverse-complemented lines (same line width), base count."""
    letters = [c for line in lines for c in line if c not in _IGNORED_IN_SEQUENCE]
    if any(c not in _SEQUENCE_LETTERS for c in letters):
        return None
    sequence = "".join(letters)
    if not sequence:
        return [], 0
    upper = sequence.upper()
    table = _RNA_COMPLEMENT if "U" in upper and "T" not in upper else _DNA_COMPLEMENT
    rc = sequence.translate(table)[::-1]
    widths = [
        len([c for c in line if c not in _IGNORED_IN_SEQUENCE])
        for line in lines
    ]
    widths = [w for w in widths if w]
    width = widths[0] if len(widths) > 1 else len(rc)
    width = max(1, width)
    out = [rc[i : i + width] for i in range(0, len(rc), width)]
    bases = sum(1 for c in rc if c not in "-.")
    return out, bases


def reverse_complement_text(text: str) -> tuple[str, int] | None:
    """Reverse complement of pasted sequence text, or None if it is not DNA/RNA.

    Plain sequence: spaces, line breaks and numbers are ignored, and the result
    keeps the line width of the original (one line stays one line). FASTA:
    every record is reverse-complemented and its header gets
    " reverse_complement", as Sequence > Reverse Complement does. IUPAC codes,
    lower case and gaps are kept; RNA (U, no T) stays RNA. Returns the text and
    the number of bases.
    """
    if not text or not text.strip():
        return None
    lines = text.replace("\r\n", "\n").replace("\r", "\n").split("\n")
    trailing_newline = text.endswith(("\n", "\r"))
    if trailing_newline and lines and lines[-1] == "":
        lines = lines[:-1]
    if not any(line.lstrip().startswith(">") for line in lines):
        result = _reverse_complement_block(lines)
        if result is None or not result[0]:
            return None
        out, bases = result
        return "\n".join(out) + ("\n" if trailing_newline else ""), bases

    blocks: list[tuple[str | None, list[str]]] = []
    for line in lines:
        if line.lstrip().startswith(">"):
            blocks.append((line.strip(), []))
        elif blocks:
            blocks[-1][1].append(line)
        elif line.strip():
            blocks.append((None, [line]))
    out_lines: list[str] = []
    total = 0
    for header, seq_lines in blocks:
        result = _reverse_complement_block(seq_lines)
        if result is None:
            return None
        rc_lines, bases = result
        total += bases
        if header is not None:
            out_lines.append(f"{header} reverse_complement")
        out_lines.extend(rc_lines)
    if total == 0:
        return None
    return "\n".join(out_lines) + ("\n" if trailing_newline else ""), total


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
