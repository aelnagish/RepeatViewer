"""User-configurable colors for annotation feature types and similarity databases."""

from __future__ import annotations

import json

from PySide6.QtCore import QSettings

from teviewer.backends.sequence_features import FEATURE_TYPE_COLORS
from teviewer.backends.similarity_annotation import SIMILARITY_DB_COLORS

_SETTINGS_GROUP = "annotation_colors"
_FEATURE_TYPES_KEY = "feature_types_json"
_SIMILARITY_DBS_KEY = "similarity_databases_json"

ANNOTATION_FEATURE_TYPE_LABELS: tuple[tuple[str, str], ...] = (
    ("LTR", "LTR"),
    ("TIR", "TIR"),
    ("ORF", "ORF"),
    ("Tandem", "tandem repeat"),
    ("PBS", "PBS"),
    ("PPT", "PPT"),
    ("Domain", "domain"),
    ("Others", "others"),
)


def _load_json_dict(key: str) -> dict[str, str]:
    s = QSettings()
    s.beginGroup(_SETTINGS_GROUP)
    raw = s.value(key, "{}")
    s.endGroup()
    if not raw:
        return {}
    try:
        decoded = json.loads(str(raw))
    except (json.JSONDecodeError, TypeError):
        return {}
    if not isinstance(decoded, dict):
        return {}
    colors: dict[str, str] = {}
    for item_key, value in decoded.items():
        color = str(value).strip()
        if item_key and color.startswith("#") and len(color) >= 4:
            colors[str(item_key)] = color
    return colors


def _save_json_dict(key: str, colors: dict[str, str]) -> None:
    s = QSettings()
    s.beginGroup(_SETTINGS_GROUP)
    s.setValue(key, json.dumps(colors, indent=2, sort_keys=True))
    s.endGroup()


def default_feature_type_colors() -> dict[str, str]:
    return dict(FEATURE_TYPE_COLORS)


def load_feature_type_colors() -> dict[str, str]:
    colors = default_feature_type_colors()
    colors.update(_load_json_dict(_FEATURE_TYPES_KEY))
    return colors


def save_feature_type_colors(colors: dict[str, str]) -> None:
    defaults = default_feature_type_colors()
    overrides = {
        key: value
        for key, value in colors.items()
        if defaults.get(key, "") != value and value.startswith("#")
    }
    _save_json_dict(_FEATURE_TYPES_KEY, overrides)


def feature_type_color(feature_type: str) -> str:
    colors = load_feature_type_colors()
    return colors.get(feature_type, colors.get("others", FEATURE_TYPE_COLORS["others"]))


def load_similarity_database_colors() -> dict[str, str]:
    return _load_json_dict(_SIMILARITY_DBS_KEY)


def save_similarity_database_colors(colors: dict[str, str]) -> None:
    cleaned = {
        name: color
        for name, color in colors.items()
        if name and str(color).startswith("#")
    }
    _save_json_dict(_SIMILARITY_DBS_KEY, cleaned)


def similarity_database_color(database_name: str, index: int = 0) -> str:
    stored = load_similarity_database_colors()
    if database_name in stored:
        return stored[database_name]
    return SIMILARITY_DB_COLORS[index % len(SIMILARITY_DB_COLORS)]


def snapshot_similarity_database_colors(database_names: list[str]) -> dict[str, str]:
    return {
        name: similarity_database_color(name, index)
        for index, name in enumerate(database_names)
    }
