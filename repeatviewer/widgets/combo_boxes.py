"""Combo boxes with application-specific input behaviour."""

from __future__ import annotations

from PySide6.QtGui import QWheelEvent
from PySide6.QtWidgets import QComboBox


class NoWheelComboBox(QComboBox):
    """Ignore wheel input so a containing scroll area can handle it."""

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802 (Qt API)
        event.ignore()


def size_export_format_combo(combo: QComboBox, *, extra: int = 44) -> None:
    """Widen Type so PDF / PNG / SVG / GFF3 labels are fully visible (closed + popup).

    QComboBox defaults to a Preferred width near the shortest label (~55 px), so
    the macOS/Fusion popup truncates longer items to "…".
    """
    combo.setSizeAdjustPolicy(QComboBox.SizeAdjustPolicy.AdjustToContents)
    fm = combo.fontMetrics()
    text_w = max(
        (fm.horizontalAdvance(combo.itemText(i)) for i in range(combo.count())),
        default=fm.horizontalAdvance("PDF"),
    )
    width = max(100, int(text_w) + int(extra))
    combo.setMinimumWidth(width)
    view = combo.view()
    if view is not None:
        view.setMinimumWidth(width)
