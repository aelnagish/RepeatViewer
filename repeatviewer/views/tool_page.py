"""A tool window shown as a tab in the middle instead of floating on top.

Some tools (the Primer Database) are whole windows of their own. Rather than
rewriting them, this page holds one inside a tab: the same widget, the same
behaviour, no second window. Its own Close button is hidden, because a tab is
closed by the tab bar.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QDialog, QPushButton, QVBoxLayout, QWidget


class ToolDialogPage(QWidget):
    """``inner`` is shown as the page's whole content."""

    def __init__(self, inner: QDialog, name: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.path = None
        self._name = name
        self.inner = inner
        inner.setParent(self)
        inner.setWindowFlags(Qt.WindowType.Widget)
        inner.setSizeGripEnabled(False) if hasattr(inner, "setSizeGripEnabled") else None
        for button in inner.findChildren(QPushButton):
            if button.text().strip() in ("Close", "Cancel"):
                button.hide()
        lay = QVBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(0)
        lay.addWidget(inner)
        inner.show()

    def tab_name(self) -> str:
        return self._name

    def closeEvent(self, event) -> None:  # noqa: N802
        """Closing the tab closes the tool, so it can save or ask as it would."""
        self.inner.close()
        super().closeEvent(event)
