"""Persist user-defined databases and CDD paths (no UI dependencies)."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QSettings

from repeatviewer.project import (
    BLAST_MOLECULE_NUCL,
    DATABASE_TYPE_CDD,
    DATABASE_TYPE_HMM,
    DATABASE_TYPE_PROTEIN,
    BlastDatabase,
    RepeatViewerProject,
    UserDatabase,
)

SETTINGS_GROUP = "project"
CDD_DIRS_KEY = "cdd_dirs_json"
USER_DATABASES_KEY = "user_databases_json"
BLAST_DATABASES_KEY = "blast_databases_json"


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
        databases.append(
            UserDatabase(
                name=name,
                path=path,
                db_type=db_type,
                is_hmm=is_hmm,
            )
        )
    return databases


def load_similarity_annotation_databases() -> list[UserDatabase]:
    """User databases plus CDD when a CDD directory is configured in Preferences."""
    databases = load_user_databases_from_settings()
    cdd_dirs = load_cdd_dirs_from_settings()
    if not cdd_dirs:
        return databases
    cdd_dir = cdd_dirs[0].expanduser().resolve()
    if not str(cdd_dir):
        return databases
    return [
        *databases,
        UserDatabase(name="CDD", path=cdd_dir, db_type=DATABASE_TYPE_CDD),
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
                }
                for db in blast_databases
            ],
            indent=2,
        ),
    )
    if owns_settings:
        s.endGroup()


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
                }
                for db in user_databases
            ],
            indent=2,
        ),
    )
    if owns_settings:
        s.endGroup()
