"""User-configurable colors for annotation feature types and similarity databases."""

from __future__ import annotations

import json
import re

from PySide6.QtCore import QSettings

from repeatviewer.backends.sequence_features import FEATURE_TYPE_COLORS
from repeatviewer.backends.similarity_annotation import SIMILARITY_DB_COLORS

_SETTINGS_GROUP = "annotation_colors"
_FEATURE_TYPES_KEY = "feature_types_json"
_SIMILARITY_DBS_KEY = "similarity_databases_json"
_ADD_FEATURE_BUTTONS_KEY = "add_feature_buttons_json"

ANNOTATION_FEATURE_TYPE_LABELS: tuple[tuple[str, str], ...] = (
    ("LTR", "LTR"),
    ("TIR", "TIR"),
    ("ORF", "ORF"),
    ("TSD", "TSD"),
    ("Tandem", "tandem repeat"),
    ("PBS", "PBS"),
    ("PPT", "PPT"),
    ("Domain", "domain"),
    ("Others", "others"),
)

# Buttons shown on the sequence-map "Add feature" row (Domain is similarity-only).
DEFAULT_ADD_FEATURE_BUTTONS: tuple[tuple[str, str], ...] = (
    ("LTR", "LTR"),
    ("TIR", "TIR"),
    ("ORF", "ORF"),
    ("Tandem", "tandem repeat"),
    ("PBS", "PBS"),
    ("PPT", "PPT"),
    ("Others", "others"),
)


def slug_feature_type(label: str) -> str:
    cleaned = re.sub(r"\s+", " ", str(label or "").strip())
    if not cleaned:
        return "others"
    lower = cleaned.lower()
    for display, feature_type in DEFAULT_ADD_FEATURE_BUTTONS:
        if lower == display.lower() or lower == feature_type.lower():
            return feature_type
    if lower in FEATURE_TYPE_COLORS:
        return next(
            (key for key in FEATURE_TYPE_COLORS if key.lower() == lower),
            cleaned,
        )
    return cleaned


def _normalize_add_feature_buttons(
    raw: object,
) -> list[tuple[str, str]]:
    if not isinstance(raw, list):
        return list(DEFAULT_ADD_FEATURE_BUTTONS)
    buttons: list[tuple[str, str]] = []
    seen: set[str] = set()
    for item in raw:
        label = ""
        feature_type = ""
        if isinstance(item, dict):
            label = str(item.get("label") or item.get("name") or "").strip()
            feature_type = str(item.get("feature_type") or item.get("type") or "").strip()
        elif isinstance(item, (list, tuple)) and len(item) >= 2:
            label = str(item[0]).strip()
            feature_type = str(item[1]).strip()
        elif isinstance(item, str):
            label = item.strip()
        if not label:
            continue
        if not feature_type:
            feature_type = slug_feature_type(label)
        key = feature_type.lower()
        if key in seen:
            continue
        seen.add(key)
        buttons.append((label, feature_type))
    return buttons or list(DEFAULT_ADD_FEATURE_BUTTONS)


def load_add_feature_buttons() -> list[tuple[str, str]]:
    s = QSettings()
    s.beginGroup(_SETTINGS_GROUP)
    raw = s.value(_ADD_FEATURE_BUTTONS_KEY, "")
    s.endGroup()
    if not raw:
        return list(DEFAULT_ADD_FEATURE_BUTTONS)
    try:
        decoded = json.loads(str(raw))
    except (json.JSONDecodeError, TypeError):
        return list(DEFAULT_ADD_FEATURE_BUTTONS)
    return _normalize_add_feature_buttons(decoded)


def save_add_feature_buttons(buttons: list[tuple[str, str]]) -> None:
    cleaned = _normalize_add_feature_buttons(
        [{"label": label, "feature_type": feature_type} for label, feature_type in buttons]
    )
    payload = [
        {"label": label, "feature_type": feature_type}
        for label, feature_type in cleaned
    ]
    s = QSettings()
    s.beginGroup(_SETTINGS_GROUP)
    s.setValue(_ADD_FEATURE_BUTTONS_KEY, json.dumps(payload, indent=2))
    s.endGroup()


def reset_add_feature_buttons() -> list[tuple[str, str]]:
    defaults = list(DEFAULT_ADD_FEATURE_BUTTONS)
    save_add_feature_buttons(defaults)
    return defaults


def annotation_feature_type_labels_for_colors() -> list[tuple[str, str]]:
    """Feature types shown in the colors dialog (defaults + user buttons)."""
    ordered: list[tuple[str, str]] = []
    seen: set[str] = set()
    for label, feature_type in ANNOTATION_FEATURE_TYPE_LABELS:
        key = feature_type.lower()
        if key in seen:
            continue
        seen.add(key)
        ordered.append((label, feature_type))
    for label, feature_type in load_add_feature_buttons():
        key = feature_type.lower()
        if key in seen:
            continue
        seen.add(key)
        ordered.append((label, feature_type))
    return ordered


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


def pident_color_hex(pident: float, zmin: float = 50.0, zmax: float = 100.0) -> str:
    """Map percent identity to a NCBI-style low→high hex color.

    Low identity: cooler/redder; high identity: greener (matches genome-map
    ``_pident_colorscale``).
    """
    try:
        value = float(pident)
    except (TypeError, ValueError):
        value = zmin
    span = max(float(zmax) - float(zmin), 1.0)
    t = (value - float(zmin)) / span
    t = max(0.0, min(1.0, t))
    r = int(230 - 90 * t)
    g = int(90 + 140 * t)
    b = int(60 + 120 * t)
    return f"#{r:02x}{g:02x}{b:02x}"


def pident_scale_bounds(pidents: list[float]) -> tuple[float, float]:
    """Return (zmin, zmax) for coloring; ignore zeros (unset). Fallback 50–100."""
    values = [float(v) for v in pidents if v and float(v) > 0]
    if not values:
        return 50.0, 100.0
    zmin = min(values)
    zmax = max(values)
    if zmax - zmin < 1.0:
        zmin = max(0.0, zmin - 5.0)
        zmax = min(100.0, zmax + 5.0)
    return zmin, zmax
