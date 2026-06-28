"""TEviewer project paths and validation."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

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

    def validate_genome(self) -> tuple[bool, str]:
        if self.genome_file is None:
            return False, "Genome FASTA is not set (required for Blast / Extend)."
        if not self.genome_file.is_file():
            return False, f"Genome file not found: {self.genome_file}"
        return True, str(self.genome_file)


def _to_path(value: Optional[Path | str]) -> Optional[Path]:
    if value is None or value == "":
        return None
    return Path(value).expanduser().resolve()
