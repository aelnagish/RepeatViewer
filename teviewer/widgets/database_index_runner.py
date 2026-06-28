"""Background database indexing with a non-modal progress window."""

from __future__ import annotations

import threading
from typing import TYPE_CHECKING

from PySide6.QtCore import QObject, Qt, QThread, Signal, Slot
from PySide6.QtWidgets import QMessageBox, QProgressDialog, QWidget

from teviewer.backends.database_prepare import (
    DatabaseTarget,
    format_prepare_error,
    prepare_database_target,
)

if TYPE_CHECKING:
    pass


class DatabaseIndexWorker(QObject):
    progress = Signal(int, int, str)
    finished = Signal(int, bool)
    failed = Signal(str, str)

    def __init__(
        self,
        targets: list[DatabaseTarget],
        *,
        use_system_blast: bool,
        cancel_event: threading.Event,
    ) -> None:
        super().__init__()
        self._targets = list(targets)
        self._use_system_blast = use_system_blast
        self._cancel_event = cancel_event

    @Slot()
    def run(self) -> None:
        total = len(self._targets)
        success_count = 0
        for index, target in enumerate(self._targets):
            if self._cancel_event.is_set():
                self.finished.emit(success_count, True)
                return
            self.progress.emit(index, total, target.label)
            try:
                prepare_database_target(
                    target,
                    use_system_blast=self._use_system_blast,
                    validate_search=True,
                )
            except Exception as exc:  # noqa: BLE001
                self.failed.emit(target.label, format_prepare_error(exc))
                self.finished.emit(success_count, False)
                return
            success_count += 1
        self.progress.emit(total, total, "")
        self.finished.emit(success_count, False)


class DatabaseIndexRunner(QObject):
    """Run database indexing on one background thread at a time."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._host = parent
        self._thread: QThread | None = None
        self._worker: DatabaseIndexWorker | None = None
        self._progress: QProgressDialog | None = None
        self._cancel_event = threading.Event()
        self._show_success_message = True
        self._target_count = 0

    def is_running(self) -> bool:
        return self._thread is not None and self._thread.isRunning()

    def start(
        self,
        targets: list[DatabaseTarget],
        *,
        use_system_blast: bool,
        show_success: bool = True,
    ) -> bool:
        if not targets:
            return False
        if self.is_running():
            QMessageBox.information(
                self._host,
                "Database indexing",
                "Database indexing is already running in the background.",
            )
            return False

        self._show_success_message = show_success
        self._target_count = len(targets)
        self._cancel_event.clear()

        progress = QProgressDialog(
            "Preparing databases…",
            "Cancel",
            0,
            len(targets),
            self._host,
        )
        progress.setWindowTitle("Database indexing")
        progress.setWindowModality(Qt.WindowModality.NonModal)
        progress.setMinimumDuration(0)
        progress.setAutoClose(False)
        progress.setAutoReset(False)
        progress.setValue(0)
        progress.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, False)

        worker = DatabaseIndexWorker(
            targets,
            use_system_blast=use_system_blast,
            cancel_event=self._cancel_event,
        )
        thread = QThread(self._host)
        worker.moveToThread(thread)

        thread.started.connect(worker.run)
        worker.progress.connect(self._on_progress)
        worker.failed.connect(self._on_failed)
        worker.finished.connect(self._on_finished)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._cleanup)
        progress.canceled.connect(self._request_cancel)

        self._thread = thread
        self._worker = worker
        self._progress = progress

        progress.show()
        thread.start()
        return True

    def _request_cancel(self) -> None:
        self._cancel_event.set()
        if self._progress is not None:
            self._progress.setLabelText("Cancelling after the current database…")

    @Slot(int, int, str)
    def _on_progress(self, current: int, total: int, label: str) -> None:
        progress = self._progress
        if progress is None:
            return
        progress.setMaximum(max(1, total))
        progress.setValue(min(current, total))
        if label:
            step = min(current + 1, total)
            progress.setLabelText(f"Preparing {label}… ({step} of {total})")

    @Slot(str, str)
    def _on_failed(self, label: str, message: str) -> None:
        QMessageBox.critical(self._host, label, message)

    @Slot(int, bool)
    def _on_finished(self, success_count: int, cancelled: bool) -> None:
        progress = self._progress
        if progress is not None:
            progress.setValue(progress.maximum())
            progress.close()
        if not self._show_success_message or success_count <= 0:
            return
        if cancelled:
            QMessageBox.information(
                self._host,
                "Database indexing",
                f"Prepared {success_count} database(s) before cancellation.",
            )
            return
        if success_count >= self._target_count:
            QMessageBox.information(
                self._host,
                "Database indexing",
                f"Prepared {success_count} database(s).",
            )

    def _cleanup(self) -> None:
        self._thread = None
        self._worker = None
        self._progress = None


def find_database_index_runner(widget: QWidget | None) -> DatabaseIndexRunner:
    current: QWidget | None = widget
    while current is not None:
        runner = getattr(current, "database_index_runner", None)
        if isinstance(runner, DatabaseIndexRunner):
            return runner
        current = current.parent()
    return DatabaseIndexRunner(widget)
