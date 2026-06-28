"""Compact job progress for the status bar (right side)."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QHBoxLayout, QLabel, QProgressBar, QWidget


def truncate_filename(name: str, max_len: int = 20) -> str:
    if len(name) <= max_len:
        return name
    return name[:max_len] + "…"


class TaskProgressPanel(QWidget):
    """Task label (always) and optional batch progress bar (2+ sequences)."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(8, 0, 4, 0)
        layout.setSpacing(8)

        self._label = QLabel("")
        self._label.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)

        self._bar = QProgressBar()
        self._bar.setTextVisible(True)
        self._bar.setMinimum(0)
        self._bar.setMaximum(1)
        self._bar.setValue(0)
        self._bar.setFixedWidth(200)
        self._bar.setFixedHeight(16)

        layout.addWidget(self._label)
        layout.addWidget(self._bar)

        self._job_title = ""
        self._total = 1
        self.hide()

    def start(self, job_title: str, total: int) -> None:
        self._job_title = job_title
        self._total = max(1, total)
        show_bar = self._total > 1
        self._bar.setVisible(show_bar)
        if show_bar:
            self._bar.setMaximum(self._total)
            self._bar.setValue(0)
            self._bar.setFormat(f"0 / {self._total}")
        self._set_running_label("")
        self.show()

    def update(self, completed: int, current_file: str = "") -> None:
        completed = max(0, min(completed, self._total))
        if self._total > 1:
            self._bar.setValue(completed)
            remaining = self._total - completed
            self._bar.setFormat(f"{completed} done, {remaining} left")
        self._set_running_label(current_file)

    def reset(self) -> None:
        self._bar.setValue(0)
        self._label.clear()
        self.hide()

    def _set_running_label(self, filename: str) -> None:
        short = truncate_filename(filename, 20) if filename else ""
        if short:
            self._label.setText(f"{self._job_title} running {short}")
        else:
            self._label.setText(f"{self._job_title} running")
