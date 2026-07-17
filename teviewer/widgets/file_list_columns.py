"""Persisted visible columns for the Files table."""

from __future__ import annotations

from PySide6.QtCore import QSettings

from teviewer.widgets.file_list_table import (
    COL_CHECK,
    COL_DESCRIPTION,
    COL_GC,
    COL_LENGTH,
    COL_MODIFIED,
    COL_NAME,
    COL_SEQUENCES,
    COL_SIZE,
    COL_TYPE,
)

FILES_SETTINGS_GROUP = "files"

# Checkbox + Name are always shown.
ALWAYS_VISIBLE_COLUMNS = frozenset({COL_CHECK, COL_NAME})

TOGGLEABLE_COLUMNS: dict[int, str] = {
    COL_LENGTH: "Length",
    COL_SEQUENCES: "# Sequences",
    COL_GC: "% GC",
    COL_TYPE: "Molecule Type",
    COL_DESCRIPTION: "Description",
    COL_MODIFIED: "Modified",
    COL_SIZE: "Size",
}

DEFAULT_VISIBLE_COLUMNS = frozenset(
    ALWAYS_VISIBLE_COLUMNS | set(TOGGLEABLE_COLUMNS.keys())
)


def load_visible_file_columns() -> set[int]:
    s = QSettings()
    s.beginGroup(FILES_SETTINGS_GROUP)
    raw = s.value("visible_columns")
    s.endGroup()
    if not raw:
        return set(DEFAULT_VISIBLE_COLUMNS)
    try:
        if isinstance(raw, (list, tuple)):
            values = [int(v) for v in raw]
        else:
            values = [int(part) for part in str(raw).split(",") if part.strip()]
    except (TypeError, ValueError):
        return set(DEFAULT_VISIBLE_COLUMNS)
    visible = set(ALWAYS_VISIBLE_COLUMNS) | {
        column for column in values if column in TOGGLEABLE_COLUMNS
    }
    return visible


def save_visible_file_columns(columns: set[int]) -> None:
    visible = sorted(
        column
        for column in columns
        if column in TOGGLEABLE_COLUMNS or column in ALWAYS_VISIBLE_COLUMNS
    )
    s = QSettings()
    s.beginGroup(FILES_SETTINGS_GROUP)
    s.setValue("visible_columns", ",".join(str(column) for column in visible))
    s.endGroup()
