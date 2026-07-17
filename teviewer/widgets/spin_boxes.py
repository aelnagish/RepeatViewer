"""Spin boxes that ignore mouse-wheel changes unless the user uses the arrows."""

from __future__ import annotations

from PySide6.QtGui import QWheelEvent
from PySide6.QtWidgets import QDoubleSpinBox, QSpinBox


class NoWheelSpinBox(QSpinBox):
    """Do not change value when the scroll wheel passes over the widget."""

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802 (Qt API)
        event.ignore()


class NoWheelDoubleSpinBox(QDoubleSpinBox):
    """Do not change value when the scroll wheel passes over the widget."""

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802 (Qt API)
        event.ignore()
