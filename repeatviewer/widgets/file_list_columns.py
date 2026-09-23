"""Persisted visible columns for the Files table."""

from __future__ import annotations

from PySide6.QtCore import QSettings

from repeatviewer.widgets.file_list_table import (
    COL_CHECK,
    COL_COLOR,
    COL_DESCRIPTION,
    COL_GC,
    COL_IMAGES,
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
    COL_TYPE: "Type",
    COL_DESCRIPTION: "Description",
    COL_MODIFIED: "Modified",
    COL_SIZE: "Size",
    COL_COLOR: "Color",
    COL_IMAGES: "Cytogenetic images",
}

DEFAULT_VISIBLE_COLUMNS = frozenset(
    ALWAYS_VISIBLE_COLUMNS | set(TOGGLEABLE_COLUMNS.keys())
)


def metadata_columns_visible(visible: set[int] | None = None) -> bool:
    """True when any optional metadata column is shown."""
    if visible is None:
        visible = load_visible_file_columns()
    # The images tick column is filled from the image store, not file metadata.
    return bool(set(visible) & (set(TOGGLEABLE_COLUMNS.keys()) - {COL_IMAGES}))


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
    if not _images_column_migration_done():
        # One-shot: the images column did not exist when this list was saved,
        # so show it once. Later hides by the user are respected.
        visible.add(COL_IMAGES)
        save_visible_file_columns(visible)
        _mark_images_column_migration_done()
    return visible


_IMAGES_COLUMN_MIGRATION_KEY = "images_column_migrated_v1"


def _images_column_migration_done() -> bool:
    s = QSettings()
    s.beginGroup(FILES_SETTINGS_GROUP)
    done = s.value(_IMAGES_COLUMN_MIGRATION_KEY, False)
    s.endGroup()
    return str(done).lower() in {"1", "true"}


def _mark_images_column_migration_done() -> None:
    s = QSettings()
    s.beginGroup(FILES_SETTINGS_GROUP)
    s.setValue(_IMAGES_COLUMN_MIGRATION_KEY, True)
    s.endGroup()


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


def load_column_widths() -> dict[int, int]:
    """Persisted section widths keyed by logical column index."""
    s = QSettings()
    s.beginGroup(FILES_SETTINGS_GROUP)
    raw = s.value("column_widths")
    s.endGroup()
    if not raw:
        return {}
    widths: dict[int, int] = {}
    try:
        text = str(raw)
        for part in text.split(","):
            part = part.strip()
            if not part or ":" not in part:
                continue
            col_text, width_text = part.split(":", 1)
            column = int(col_text)
            width = int(width_text)
            if width >= 22:
                widths[column] = width
    except (TypeError, ValueError):
        return {}
    return widths


def save_column_widths(widths: dict[int, int]) -> None:
    s = QSettings()
    s.beginGroup(FILES_SETTINGS_GROUP)
    payload = ",".join(f"{column}:{width}" for column, width in sorted(widths.items()))
    s.setValue("column_widths", payload)
    s.endGroup()


def load_column_order() -> list[int] | None:
    """Logical column indices in left-to-right visual order."""
    s = QSettings()
    s.beginGroup(FILES_SETTINGS_GROUP)
    raw = s.value("column_order")
    s.endGroup()
    if not raw:
        return None
    try:
        order = [int(part) for part in str(raw).split(",") if part.strip()]
    except (TypeError, ValueError):
        return None
    if not order:
        return None
    return order


def save_column_order(order: list[int]) -> None:
    s = QSettings()
    s.beginGroup(FILES_SETTINGS_GROUP)
    s.setValue("column_order", ",".join(str(column) for column in order))
    s.endGroup()
