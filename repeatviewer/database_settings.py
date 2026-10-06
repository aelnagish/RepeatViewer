"""Persist user-defined databases and CDD paths (no UI dependencies)."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QSettings

from repeatviewer.project import (
    BLAST_MOLECULE_NUCL,
    DATABASE_SCOPE_PUBLIC,
    DATABASE_SCOPE_USER,
    DATABASE_TYPE_CDD,
    DATABASE_TYPE_DNA,
    DATABASE_TYPE_HMM,
    DATABASE_TYPE_PROTEIN,
    BlastDatabase,
    ManagedGenome,
    RepeatViewerProject,
    UserDatabase,
    normalize_database_scope,
)

SETTINGS_GROUP = "project"
CDD_DIRS_KEY = "cdd_dirs_json"
USER_DATABASES_KEY = "user_databases_json"
BLAST_DATABASES_KEY = "blast_databases_json"
MANAGED_GENOMES_KEY = "managed_genomes_json"
GENOME_FILE_HISTORY_KEY = "genome_file_history_json"
MAX_MANAGED_GENOMES = 50


def _path_from_settings(value) -> Path | None:
    if not value:
        return None
    text = str(value).strip()
    if not text:
        return None
    return Path(text).expanduser()


def load_cdd_dirs_from_settings(s: QSettings | None = None) -> list[Path]:
    owns_settings = s is None
    if owns_settings:
        s = QSettings()
        s.beginGroup(SETTINGS_GROUP)
    raw_json = s.value(CDD_DIRS_KEY)
    paths: list[Path] = []
    if raw_json:
        try:
            decoded = json.loads(str(raw_json))
            if isinstance(decoded, list):
                paths = [
                    p for p in (_path_from_settings(item) for item in decoded) if p is not None
                ]
        except (json.JSONDecodeError, TypeError):
            paths = []
    if not paths:
        single = _path_from_settings(s.value("cdd_dir"))
        if single is not None:
            paths = [single]
    if owns_settings:
        s.endGroup()
    return paths


def load_user_databases_from_settings() -> list[UserDatabase]:
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    raw = s.value(USER_DATABASES_KEY)
    s.endGroup()
    if not raw:
        return []
    try:
        decoded = json.loads(str(raw))
    except (json.JSONDecodeError, TypeError):
        return []
    if not isinstance(decoded, list):
        return []
    databases: list[UserDatabase] = []
    for item in decoded:
        if not isinstance(item, dict):
            continue
        name = str(item.get("name", "")).strip()
        path = _path_from_settings(item.get("path"))
        if not name or path is None:
            continue
        db_type = str(item.get("db_type", DATABASE_TYPE_PROTEIN))
        is_hmm = bool(item.get("is_hmm", False))
        if db_type == DATABASE_TYPE_HMM:
            db_type = DATABASE_TYPE_PROTEIN
            is_hmm = True
        # CDD migrated from Preferences gets Public by default; others User.
        default_scope = (
            DATABASE_SCOPE_PUBLIC
            if db_type == DATABASE_TYPE_CDD
            else DATABASE_SCOPE_USER
        )
        databases.append(
            UserDatabase(
                name=name,
                path=path,
                db_type=db_type,
                is_hmm=is_hmm,
                scope=str(item.get("scope", default_scope)),
                created_at=str(item.get("created_at", "")),
                show_in_annotation=bool(item.get("show_in_annotation", True)),
                annotation_name=str(item.get("annotation_name", "")),
            )
        )
    return databases


def load_similarity_annotation_databases() -> list[UserDatabase]:
    """Annotation-tab databases: user DBs + BLAST rows marked show_in_annotation.

    Display name prefers ``annotation_name`` (falls back to file name).
    """
    databases: list[UserDatabase] = []
    seen_names: set[str] = set()
    for db in load_user_databases_from_settings():
        if not db.show_in_annotation:
            continue
        label = db.annotation_label
        if label.lower() in seen_names:
            continue
        databases.append(
            UserDatabase(
                name=label,
                path=db.path,
                db_type=db.db_type,
                is_hmm=db.is_hmm,
                scope=db.scope,
                created_at=db.created_at,
                show_in_annotation=True,
                annotation_name=db.annotation_name,
            )
        )
        seen_names.add(label.lower())
    for blast in load_blast_databases_from_settings():
        if not blast.show_in_annotation:
            continue
        label = blast.annotation_label
        if label.lower() in seen_names:
            continue
        databases.append(
            UserDatabase(
                name=label,
                path=blast.source_fasta,
                db_type=(
                    DATABASE_TYPE_DNA
                    if blast.molecule_type == BLAST_MOLECULE_NUCL
                    else DATABASE_TYPE_PROTEIN
                ),
                created_at=blast.created_at,
                show_in_annotation=True,
                annotation_name=blast.annotation_name,
            )
        )
        seen_names.add(label.lower())
    cdd_dirs = load_cdd_dirs_from_settings()
    if not cdd_dirs:
        return databases
    cdd_dir = cdd_dirs[0].expanduser().resolve()
    if not str(cdd_dir):
        return databases
    # Prefer an explicit CDD row in user_databases if present.
    for db in databases:
        if db.db_type == DATABASE_TYPE_CDD:
            return databases
    # Also skip synthetic CDD if a hidden CDD row exists in settings.
    for db in load_user_databases_from_settings():
        if db.db_type == DATABASE_TYPE_CDD:
            return databases
    return [
        *databases,
        UserDatabase(
            name="CDD",
            path=cdd_dir,
            db_type=DATABASE_TYPE_CDD,
            scope=DATABASE_SCOPE_PUBLIC,
            show_in_annotation=True,
            annotation_name="CDD",
        ),
    ]


def load_blast_databases_from_settings() -> list[BlastDatabase]:
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    raw = s.value(BLAST_DATABASES_KEY)
    s.endGroup()
    if not raw:
        return []
    try:
        decoded = json.loads(str(raw))
    except (json.JSONDecodeError, TypeError):
        return []
    if not isinstance(decoded, list):
        return []
    databases: list[BlastDatabase] = []
    for item in decoded:
        if not isinstance(item, dict):
            continue
        name = str(item.get("name", "")).strip()
        source_fasta = _path_from_settings(item.get("source_fasta"))
        db_prefix_path = _path_from_settings(item.get("db_prefix_path"))
        if not name or source_fasta is None or db_prefix_path is None:
            continue
        databases.append(
            BlastDatabase(
                name=name,
                source_fasta=source_fasta,
                db_prefix_path=db_prefix_path,
                molecule_type=str(item.get("molecule_type", BLAST_MOLECULE_NUCL)),
                created_at=str(item.get("created_at", "")),
                scope=str(item.get("scope", DATABASE_SCOPE_USER)),
                show_in_annotation=bool(item.get("show_in_annotation", False)),
                annotation_name=str(item.get("annotation_name", "")),
            )
        )
    return databases


def reload_blast_databases_into_project(project: RepeatViewerProject) -> list[BlastDatabase]:
    """Replace ``project.blast_databases`` with the latest saved list.

    Call before showing the BLAST Database dropdown so Manage / Add / Remove
    (and other windows) always appear in the Run dialog.
    """
    databases = load_blast_databases_from_settings()
    project.blast_databases = databases
    return databases


def save_blast_databases_to_settings(
    blast_databases: list[BlastDatabase],
    *,
    s: QSettings | None = None,
) -> None:
    owns_settings = s is None
    if owns_settings:
        s = QSettings()
        s.beginGroup(SETTINGS_GROUP)
    s.setValue(
        BLAST_DATABASES_KEY,
        json.dumps(
            [
                {
                    "name": db.name,
                    "source_fasta": str(db.source_fasta),
                    "db_prefix_path": str(db.db_prefix_path),
                    "molecule_type": db.molecule_type,
                    "created_at": db.created_at,
                    "scope": normalize_database_scope(db.scope),
                    "show_in_annotation": bool(db.show_in_annotation),
                    "annotation_name": db.annotation_name,
                }
                for db in blast_databases
            ],
            indent=2,
        ),
    )
    if owns_settings:
        s.endGroup()


def _genome_path_key(path: Path) -> str:
    try:
        return str(path.expanduser().resolve())
    except OSError:
        return str(path.expanduser())


def load_managed_genomes_from_settings() -> list[ManagedGenome]:
    """Load managed genomes; migrate from genome_file_history_json when empty."""
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    raw = s.value(MANAGED_GENOMES_KEY)
    history_raw = s.value(GENOME_FILE_HISTORY_KEY)
    current = _path_from_settings(s.value("genome_file"))
    s.endGroup()

    genomes: list[ManagedGenome] = []
    if raw:
        try:
            decoded = json.loads(str(raw))
        except (json.JSONDecodeError, TypeError):
            decoded = None
        if isinstance(decoded, list):
            for item in decoded:
                if not isinstance(item, dict):
                    continue
                name = str(item.get("name", "")).strip()
                path = _path_from_settings(item.get("path"))
                if path is None:
                    continue
                genomes.append(
                    ManagedGenome(
                        name=name or path.stem,
                        path=path,
                        scope=str(item.get("scope", DATABASE_SCOPE_USER)),
                        created_at=str(item.get("created_at", "")),
                    )
                )

    if not genomes and history_raw:
        try:
            decoded = json.loads(str(history_raw))
        except (json.JSONDecodeError, TypeError):
            decoded = None
        if isinstance(decoded, list):
            for item in decoded:
                path = _path_from_settings(item)
                if path is None:
                    continue
                genomes.append(
                    ManagedGenome(
                        name=path.stem,
                        path=path,
                        scope=DATABASE_SCOPE_USER,
                    )
                )

    if current is not None:
        key = _genome_path_key(current)
        if not any(_genome_path_key(g.path) == key for g in genomes):
            genomes.insert(
                0,
                ManagedGenome(
                    name=current.stem,
                    path=current,
                    scope=DATABASE_SCOPE_USER,
                ),
            )

    # Deduplicate by resolved path, prefer first.
    seen: set[str] = set()
    unique: list[ManagedGenome] = []
    for g in genomes:
        key = _genome_path_key(g.path)
        if key in seen:
            continue
        seen.add(key)
        unique.append(g)
    return unique[:MAX_MANAGED_GENOMES]


def save_managed_genomes_to_settings(
    genomes: list[ManagedGenome],
    *,
    s: QSettings | None = None,
) -> None:
    owns_settings = s is None
    if owns_settings:
        s = QSettings()
        s.beginGroup(SETTINGS_GROUP)
    unique: list[ManagedGenome] = []
    seen: set[str] = set()
    for g in genomes:
        key = _genome_path_key(g.path)
        if key in seen:
            continue
        seen.add(key)
        unique.append(g)
    unique = unique[:MAX_MANAGED_GENOMES]
    s.setValue(
        MANAGED_GENOMES_KEY,
        json.dumps(
            [
                {
                    "name": g.name,
                    "path": str(g.path),
                    "scope": normalize_database_scope(g.scope),
                    "created_at": g.created_at,
                }
                for g in unique
            ],
            indent=2,
        ),
    )
    # Keep legacy history in sync for GenomeDialog / other consumers.
    s.setValue(
        GENOME_FILE_HISTORY_KEY,
        json.dumps([str(g.path) for g in unique], indent=2),
    )
    if owns_settings:
        s.endGroup()


def reload_managed_genomes_into_project(project: RepeatViewerProject) -> list[ManagedGenome]:
    genomes = load_managed_genomes_from_settings()
    project.managed_genomes = genomes
    return genomes


def save_database_settings(
    *,
    cdd_dirs: list[Path],
    user_databases: list[UserDatabase],
    s: QSettings | None = None,
) -> None:
    owns_settings = s is None
    if owns_settings:
        s = QSettings()
        s.beginGroup(SETTINGS_GROUP)
    s.setValue("cdd_dir", str(cdd_dirs[0]) if cdd_dirs else "")
    s.setValue(CDD_DIRS_KEY, json.dumps([str(p) for p in cdd_dirs], indent=2))
    s.setValue(
        USER_DATABASES_KEY,
        json.dumps(
            [
                {
                    "name": db.name,
                    "path": str(db.path),
                    "db_type": db.db_type,
                    "is_hmm": db.is_hmm,
                    "scope": normalize_database_scope(db.scope),
                    "created_at": db.created_at,
                    "show_in_annotation": bool(db.show_in_annotation),
                    "annotation_name": db.annotation_name,
                }
                for db in user_databases
            ],
            indent=2,
        ),
    )
    if owns_settings:
        s.endGroup()
