"""Shared helpers for tool output paths and conflict numbering."""

from __future__ import annotations

from pathlib import Path


def next_numbered_path(path: Path) -> Path:
    """Return ``path`` if free, else ``stem_2.ext``, ``stem_3.ext``, …

    Numbering is inserted before the final suffix::

        reads.fa_FastGA.paf → reads.fa_FastGA_2.paf
        reads.fa_TEAid.html → reads.fa_TEAid_2.html
        reads.fa_0.8_co.fa  → reads.fa_0.8_co_2.fa
    """
    path = Path(path)
    if not path.exists():
        return path
    suffix = path.suffix
    base = path.name[: -len(suffix)] if suffix else path.name
    parent = path.parent
    n = 2
    while True:
        candidate = parent / (f"{base}_{n}{suffix}" if suffix else f"{base}_{n}")
        if not candidate.exists():
            return candidate
        n += 1
