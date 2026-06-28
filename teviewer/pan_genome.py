"""Single vs pan-genome configuration helpers."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

GENOME_MODE_SINGLE = "single"
GENOME_MODE_PAN = "pan"
MAX_PAN_GENOME_NAME_LEN = 10
_PAN_GENOME_NAME_RE = re.compile(r"^[A-Za-z0-9_\-]+$")


@dataclass
class PanGenomeEntry:
    """One genome in pan-genome mode (short label + FASTA path)."""

    name: str
    path: Path

    def __post_init__(self) -> None:
        self.name = normalize_pan_genome_name(self.name)
        self.path = Path(self.path).expanduser()


def normalize_pan_genome_name(name: str) -> str:
    text = str(name).strip()
    if not text:
        raise ValueError("Genome name is required.")
    if len(text) > MAX_PAN_GENOME_NAME_LEN:
        raise ValueError(
            f"Genome name must be at most {MAX_PAN_GENOME_NAME_LEN} characters."
        )
    if not _PAN_GENOME_NAME_RE.match(text):
        raise ValueError(
            "Genome name may only contain letters, digits, underscore, and hyphen."
        )
    return text


def genome_output_suffix(label: str) -> str:
    """Filename suffix for pan-genome outputs, e.g. ``_Ref``."""
    if not label:
        return ""
    safe = normalize_pan_genome_name(label)
    return f"_{safe}"


def tagged_output_stem(base: str, label: str) -> str:
    if not label:
        return base
    return f"{base}{genome_output_suffix(label)}"


def parse_pan_genomes_json(raw) -> list[PanGenomeEntry]:
    if not raw:
        return []
    if isinstance(raw, str):
        import json

        try:
            raw = json.loads(raw)
        except json.JSONDecodeError:
            return []
    if not isinstance(raw, list):
        return []
    entries: list[PanGenomeEntry] = []
    seen: set[str] = set()
    for item in raw:
        if not isinstance(item, dict):
            continue
        name = str(item.get("name", "")).strip()
        path_raw = item.get("path")
        if not name or not path_raw:
            continue
        try:
            entry = PanGenomeEntry(name=name, path=Path(str(path_raw)))
        except ValueError:
            continue
        if not entry.path.is_file():
            continue
        if entry.name in seen:
            continue
        seen.add(entry.name)
        entries.append(entry)
    return entries


def pan_genomes_to_json(entries: list[PanGenomeEntry]) -> str:
    import json

    payload = [{"name": e.name, "path": str(e.path.expanduser().resolve())} for e in entries]
    return json.dumps(payload)
