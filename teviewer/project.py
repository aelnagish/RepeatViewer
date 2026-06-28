"""TEviewer project paths and validation."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from teviewer.pan_genome import (
    GENOME_MODE_PAN,
    GENOME_MODE_SINGLE,
    PanGenomeEntry,
    normalize_pan_genome_name,
)

DATABASE_TYPE_HMM = "HMM"  # legacy persisted value; migrated to Protein/DNA + is_hmm
DATABASE_TYPE_PROTEIN = "Protein"
DATABASE_TYPE_DNA = "DNA"
DATABASE_TYPE_CDD = "CDD"
USER_DATABASE_TYPES = (
    DATABASE_TYPE_PROTEIN,
    DATABASE_TYPE_DNA,
)
SIMILARITY_DATABASE_TYPES = USER_DATABASE_TYPES + (DATABASE_TYPE_CDD,)


@dataclass
class UserDatabase:
    """User-defined sequence/HMM database from Preferences."""

    name: str
    path: Path
    db_type: str = DATABASE_TYPE_PROTEIN
    is_hmm: bool = False

    def __post_init__(self) -> None:
        self.name = str(self.name).strip()
        self.is_hmm = bool(self.is_hmm)
        self.db_type = str(self.db_type).strip() or DATABASE_TYPE_PROTEIN
        if self.db_type == DATABASE_TYPE_HMM:
            self.db_type = DATABASE_TYPE_PROTEIN
            self.is_hmm = True
        if self.db_type == DATABASE_TYPE_CDD:
            self.is_hmm = False
        elif self.db_type not in USER_DATABASE_TYPES:
            self.db_type = DATABASE_TYPE_PROTEIN
        resolved = _to_path(self.path)
        self.path = resolved if resolved is not None else Path(str(self.path))

    @property
    def display_type_label(self) -> str:
        if self.db_type == DATABASE_TYPE_CDD:
            return DATABASE_TYPE_CDD
        if self.is_hmm:
            return f"{self.db_type} HMM"
        return self.db_type

    @property
    def prepare_kind(self) -> str:
        if self.db_type == DATABASE_TYPE_CDD:
            return "cdd"
        if self.is_hmm:
            return "hmm"
        if self.db_type == DATABASE_TYPE_DNA:
            return "dna"
        return "protein"


@dataclass
class TEviewerProject:
    """Paths the user configures once per analysis session."""

    input_dir: Optional[Path] = None
    genome_file: Optional[Path] = None
    genome_mode: str = GENOME_MODE_SINGLE
    pan_genomes: list[PanGenomeEntry] = field(default_factory=list)
    cdd_dirs: list[Path] = field(default_factory=list)
    user_databases: list[UserDatabase] = field(default_factory=list)

    # Optional tool flags (mirror TEtrimmerGUI CLI)
    use_system_blast: bool = False
    use_system_sequence_viewer: bool = False
    # How many sequences to process at once in batch jobs (Blast, TEAid, etc.)
    batch_workers: int = 4

    def __post_init__(self) -> None:
        self.input_dir = _to_path(self.input_dir)
        self.genome_file = _to_path(self.genome_file)
        if self.genome_mode not in {GENOME_MODE_SINGLE, GENOME_MODE_PAN}:
            self.genome_mode = GENOME_MODE_SINGLE
        normalized_pan: list[PanGenomeEntry] = []
        seen_names: set[str] = set()
        for entry in self.pan_genomes:
            if isinstance(entry, PanGenomeEntry):
                item = entry
            elif isinstance(entry, dict):
                name = str(entry.get("name", "")).strip()
                raw_path = entry.get("path")
                if not name or not raw_path:
                    continue
                item = PanGenomeEntry(name=name, path=Path(str(raw_path)))
            else:
                continue
            if item.name in seen_names:
                continue
            seen_names.add(item.name)
            normalized_pan.append(
                PanGenomeEntry(name=item.name, path=_to_path(item.path) or item.path)
            )
        self.pan_genomes = normalized_pan
        if self.genome_mode == GENOME_MODE_PAN and self.pan_genomes:
            self.genome_file = _to_path(self.pan_genomes[0].path)
        self.cdd_dirs = [p for p in (_to_path(path) for path in self.cdd_dirs) if p is not None]
        normalized_dbs: list[UserDatabase] = []
        for entry in self.user_databases:
            if isinstance(entry, UserDatabase):
                if entry.name and entry.path:
                    normalized_dbs.append(entry)
            elif isinstance(entry, dict):
                name = str(entry.get("name", "")).strip()
                raw_path = entry.get("path")
                if not name or not raw_path:
                    continue
                normalized_dbs.append(
                    UserDatabase(
                        name=name,
                        path=Path(str(raw_path)),
                        db_type=str(entry.get("db_type", DATABASE_TYPE_PROTEIN)),
                        is_hmm=bool(entry.get("is_hmm", False)),
                    )
                )
        self.user_databases = normalized_dbs

    @property
    def cdd_dir(self) -> Optional[Path]:
        return self.cdd_dirs[0] if self.cdd_dirs else None

    @cdd_dir.setter
    def cdd_dir(self, value: Optional[Path]) -> None:
        resolved = _to_path(value)
        if resolved is None:
            self.cdd_dirs = []
        elif self.cdd_dirs:
            self.cdd_dirs[0] = resolved
        else:
            self.cdd_dirs = [resolved]

    def validate_input_dir(self) -> tuple[bool, str]:
        if self.input_dir is None:
            return False, "Input directory is not set."
        if not self.input_dir.is_dir():
            return False, f"Not a directory: {self.input_dir}"
        return True, str(self.input_dir)

    @property
    def is_pan_genome_mode(self) -> bool:
        return self.genome_mode == GENOME_MODE_PAN and bool(self.pan_genomes)

    def genome_targets(self) -> list[tuple[str, Path]]:
        """Return ``(label, fasta_path)`` pairs for genome-dependent tools."""
        if self.is_pan_genome_mode:
            return [(entry.name, entry.path.expanduser().resolve()) for entry in self.pan_genomes]
        if self.genome_file is not None and self.genome_file.is_file():
            return [("", self.genome_file.expanduser().resolve())]
        return []

    def pan_genome_by_name(self, name: str) -> PanGenomeEntry | None:
        try:
            key = normalize_pan_genome_name(name)
        except ValueError:
            return None
        for entry in self.pan_genomes:
            if entry.name == key:
                return entry
        return None

    def has_genome_configured(self) -> bool:
        if self.is_pan_genome_mode:
            return bool(self.pan_genomes) and all(
                entry.path.is_file() for entry in self.pan_genomes
            )
        return self.genome_file is not None and self.genome_file.is_file()

    def validate_genome(self) -> tuple[bool, str]:
        if self.is_pan_genome_mode:
            if not self.pan_genomes:
                return False, "Pan-genome mode requires at least one genome entry."
            for entry in self.pan_genomes:
                if not entry.path.is_file():
                    return False, f"Genome file not found for '{entry.name}': {entry.path}"
            return True, f"{len(self.pan_genomes)} pan-genome(s) configured"
        if self.genome_file is None:
            return False, "Genome FASTA is not set (required for Blast / Extend)."
        if not self.genome_file.is_file():
            return False, f"Genome file not found: {self.genome_file}"
        return True, str(self.genome_file)


def _to_path(value: Optional[Path | str]) -> Optional[Path]:
    if value is None or value == "":
        return None
    return Path(value).expanduser().resolve()
