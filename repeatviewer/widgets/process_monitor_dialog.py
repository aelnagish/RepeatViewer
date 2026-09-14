"""Non-modal window that streams a long-running tool's command and progress."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtGui import QFont, QTextCursor
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
)


class ProcessMonitorDialog(QDialog):
    """Show a command line first, then append intermediate process output."""

    def __init__(self, title: str = "Process", parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setModal(False)
        self.setWindowFlag(Qt.WindowType.Window, True)
        self.resize(720, 420)

        self._status = QLabel("Running…", self)
        self._status.setStyleSheet("color: palette(mid);")

        self._log = QPlainTextEdit(self)
        self._log.setReadOnly(True)
        self._log.setLineWrapMode(QPlainTextEdit.LineWrapMode.WidgetWidth)
        font = QFont("Menlo")
        if not font.exactMatch():
            font = QFont("Consolas")
        if not font.exactMatch():
            font = QFont("Courier New")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(11)
        self._log.setFont(font)

        self._close_btn = QPushButton("Close", self)
        self._close_btn.setEnabled(False)
        self._close_btn.clicked.connect(self.close)

        footer = QHBoxLayout()
        footer.addWidget(self._status, 1)
        footer.addWidget(self._close_btn)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(8)
        layout.addWidget(self._log, 1)
        layout.addLayout(footer)

    def prepare(self, title: str) -> None:
        self.setWindowTitle(title)
        self._log.clear()
        self._status.setText("Running…")
        self._close_btn.setEnabled(False)
        self.show()
        self.raise_()
        self.activateWindow()

    def append_line(self, text: str) -> None:
        if text is None:
            return
        line = str(text).rstrip("\r\n")
        self._log.appendPlainText(line)
        cursor = self._log.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.End)
        self._log.setTextCursor(cursor)

    def mark_finished(self, ok: bool, detail: str = "") -> None:
        if ok:
            self._status.setText("Finished.")
            self.append_line("")
            self.append_line("— finished —")
        else:
            message = detail.strip() or "failed"
            self._status.setText("Failed.")
            self.append_line("")
            self.append_line(f"— failed: {message} —")
        self._close_btn.setEnabled(True)

    def show_failure_log(
        self,
        title: str,
        lines: list[str],
        *,
        error: str = "",
    ) -> None:
        """Open (or reveal) this window with a captured command log after a failure."""
        self.setWindowTitle(title)
        self._log.clear()
        for line in lines:
            self.append_line(line)
        self.mark_finished(False, error)
        self.show()
        self.raise_()
        self.activateWindow()
