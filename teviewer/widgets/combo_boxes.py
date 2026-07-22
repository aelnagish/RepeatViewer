"""Combo boxes with application-specific input behaviour."""

from __future__ import annotations

from PySide6.QtGui import QWheelEvent
from PySide6.QtWidgets import QComboBox


class NoWheelComboBox(QComboBox):
    """Ignore wheel input so a containing scroll area can handle it."""

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802 (Qt API)
        event.ignore()
