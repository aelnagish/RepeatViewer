"""Background workers for long-running TE tools."""

from __future__ import annotations

from typing import Any, Callable

from PySide6.QtCore import QObject, QThread, Signal


class Worker(QObject):
    """
    Run a callable on a background thread.

    Usage:
        worker = Worker(blast_fn, arg1, kw=val)
        worker.signals.finished.connect(on_done)
        worker.signals.error.connect(on_err)
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        thread.start()
    """

    class Signals(QObject):
        finished = Signal(object)
        error = Signal(str)
        progress = Signal(str)

    def __init__(self, func: Callable[..., Any], *args: Any, **kwargs: Any) -> None:
        super().__init__()
        self._func = func
        self._args = args
        self._kwargs = kwargs
        self.signals = Worker.Signals()

    def run(self) -> None:
        try:
            result = self._func(*self._args, **self._kwargs)
            self.signals.finished.emit(result)
        except Exception as exc:  # noqa: BLE001 — surface to UI
            self.signals.error.emit(str(exc))
