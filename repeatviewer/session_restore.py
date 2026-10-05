"""Persist browser file selection and annotation viewport across app restarts."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QSettings

from repeatviewer.database_settings import SETTINGS_GROUP

_MAX_FOLDER_FILE_ENTRIES = 120
_MAX_ANNOTATION_VIEW_ENTRIES = 40


def _load_json_dict(raw: object) -> dict[str, object]:
    if not raw:
        return {}
    try:
        decoded = json.loads(str(raw))
    except (json.JSONDecodeError, TypeError):
        return {}
    if not isinstance(decoded, dict):
        return {}
    return decoded


def _trim_dict(d: dict[str, object], limit: int) -> dict[str, object]:
    if len(d) <= limit:
        return d
    # Keep the most recently written entries (dict preserves insertion order in Py3.7+).
    items = list(d.items())[-limit:]
    return dict(items)


def load_folder_file_memory() -> tuple[dict[str, str], dict[str, int]]:
    """Return (last_file_path_by_folder, last_scroll_by_folder) from QSettings."""
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    files_raw = _load_json_dict(s.value("folder_last_files_json"))
    scroll_raw = _load_json_dict(s.value("folder_last_scroll_json"))
    s.endGroup()
    files: dict[str, str] = {}
    for key, value in files_raw.items():
        text = str(key).strip()
        path = str(value or "").strip()
        if text and path:
            files[text] = path
    scrolls: dict[str, int] = {}
    for key, value in scroll_raw.items():
        text = str(key).strip()
        if not text:
            continue
        try:
            scrolls[text] = int(value)
        except (TypeError, ValueError):
            continue
    return files, scrolls


def save_folder_file_memory(
    last_file_by_folder: dict[str, str],
    last_scroll_by_folder: dict[str, int],
) -> None:
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    s.setValue(
        "folder_last_files_json",
        json.dumps(_trim_dict(dict(last_file_by_folder), _MAX_FOLDER_FILE_ENTRIES)),
    )
    scroll_json = {
        k: int(v) for k, v in last_scroll_by_folder.items() if str(k).strip()
    }
    s.setValue(
        "folder_last_scroll_json",
        json.dumps(_trim_dict(scroll_json, _MAX_FOLDER_FILE_ENTRIES)),
    )
    s.endGroup()


def _view_storage_key(path: Path, sequence_id: str = "") -> str:
    try:
        base = str(path.expanduser().resolve())
    except OSError:
        base = str(path)
    sid = str(sequence_id or "").strip()
    return f"{base}\0{sid}" if sid else base


def load_annotation_views() -> dict[str, dict]:
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    raw = _load_json_dict(s.value("annotation_views_by_path_json"))
    s.endGroup()
    views: dict[str, dict] = {}
    for key, value in raw.items():
        if not isinstance(value, dict):
            continue
        views[str(key)] = value
    return views


def save_annotation_view(path: Path, state: dict, *, sequence_id: str = "") -> None:
    if not state:
        return
    key = _view_storage_key(path, sequence_id)
    views = load_annotation_views()
    views[key] = dict(state)
    views = _trim_dict(views, _MAX_ANNOTATION_VIEW_ENTRIES)
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    s.setValue("annotation_views_by_path_json", json.dumps(views))
    s.endGroup()


def lookup_annotation_view(path: Path, *, sequence_id: str = "") -> dict | None:
    key = _view_storage_key(path, sequence_id)
    views = load_annotation_views()
    hit = views.get(key)
    if hit is not None:
        return hit
    if sequence_id:
        try:
            fallback = _view_storage_key(path, "")
        except OSError:
            fallback = str(path)
        return views.get(fallback)
    return None
