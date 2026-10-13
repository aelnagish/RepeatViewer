"""Compact job progress for the status bar (right side)."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QHBoxLayout, QLabel, QProgressBar, QWidget


def iter_job_progress_hosts(widget: QWidget) -> list[QWidget]:
    """Windows that should show panel-job progress (local preview + MainWindow).

    Duplicated Overview/Annotation windows are Qt Tool windows: ``widget.window()``
    is the floating preview, which historically had no progress API, so Similarity
    annotation updates were dropped. Always include the floating window (when it
    supports progress) and the MainWindow (for the status bar + tool locking).
    """
    hosts: list[QWidget] = []
    seen: set[int] = set()

    def _add(candidate: QWidget | None) -> None:
        if candidate is None:
            return
        key = id(candidate)
        if key in seen:
            return
        if not hasattr(candidate, "_show_job_progress"):
            return
        seen.add(key)
        hosts.append(candidate)

    _add(widget.window())

    parent = widget.parentWidget()
    while parent is not None:
        if hasattr(parent, "_set_job_actions_enabled") and hasattr(
            parent, "_show_job_progress"
        ):
            _add(parent)
            break
        parent = parent.parentWidget()

    if not any(hasattr(host, "_set_job_actions_enabled") for host in hosts):
        app = QApplication.instance()
        if app is not None:
            for top in app.topLevelWidgets():
                if hasattr(top, "_set_job_actions_enabled") and hasattr(
                    top, "_show_job_progress"
                ):
                    _add(top)
                    break
    return hosts


def truncate_filename(name: str, max_len: int = 20) -> str:
    if len(name) <= max_len:
        return name
    return name[:max_len] + "…"


def _strip_leading_running(text: str) -> str:
    lowered = text.lower()
    if lowered.startswith("running "):
        return text[8:].lstrip()
    return text


def format_job_status_label(job_title: str, detail: str = "") -> str:
    """Build one status line without repeating the word 'running'."""
    detail = _strip_leading_running(detail.strip())
    if not detail:
        return job_title
    return f"{job_title} {truncate_filename(detail, 20)}"


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
        self._total = max(0, int(total))
        self._bar.setVisible(True)
        if self._total <= 0:
            # Indeterminate busy bar (opening a large file, unknown steps).
            self._bar.setRange(0, 0)
            self._bar.setFormat("")
        else:
            self._bar.setRange(0, self._total)
            self._bar.setValue(0)
            if self._total > 1:
                self._bar.setFormat("%v / %m")
            else:
                self._bar.setFormat("")
        self._set_running_label("")
        self.show()

    def update(self, completed: int, current_file: str = "") -> None:
        if self._total <= 0:
            # Busy mode — keep the bar spinning; only refresh the label.
            self._set_running_label(current_file)
            return
        completed = max(0, min(completed, self._total))
        self._bar.setValue(completed)
        if self._total > 1:
            self._bar.setFormat("%v / %m")
        self._set_running_label(current_file)

    def reset(self) -> None:
        self._bar.setValue(0)
        self._label.clear()
        self.hide()

    def _set_running_label(self, filename: str) -> None:
        self._label.setText(format_job_status_label(self._job_title, filename))
