"""Persisted font size for the Folders and Files dock panels."""

from __future__ import annotations

from PySide6.QtCore import QSettings

BROWSER_PANEL_SETTINGS_GROUP = "browser_panel"
DEFAULT_BROWSER_PANEL_FONT_SIZE = 12
MIN_BROWSER_PANEL_FONT_SIZE = 8
MAX_BROWSER_PANEL_FONT_SIZE = 16


def load_browser_panel_font_size() -> int:
    s = QSettings()
    s.beginGroup(BROWSER_PANEL_SETTINGS_GROUP)
    value = int(s.value("font_size", DEFAULT_BROWSER_PANEL_FONT_SIZE))
    s.endGroup()
    return max(
        MIN_BROWSER_PANEL_FONT_SIZE,
        min(MAX_BROWSER_PANEL_FONT_SIZE, value),
    )


def save_browser_panel_font_size(size: int) -> None:
    s = QSettings()
    s.beginGroup(BROWSER_PANEL_SETTINGS_GROUP)
    s.setValue(
        "font_size",
        max(
            MIN_BROWSER_PANEL_FONT_SIZE,
            min(MAX_BROWSER_PANEL_FONT_SIZE, int(size)),
        ),
    )
    s.endGroup()
