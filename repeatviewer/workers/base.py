"""Background workers for long-running TE tools."""

from __future__ import annotations

from typing import Any, Callable

from PySide6.QtCore import QObject, Signal, Slot


class Worker(QObject):
    """
    Run a callable on a background thread.

    Usage:
        worker = Worker(blast_fn, arg1, kw=val)
        worker.signals.finished.connect(on_done)  # prefer a QObject @Slot
        worker.signals.error.connect(on_err)
        thread = QThread()
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        thread.start()

    Signals live on this QObject. After ``moveToThread``, emits from ``run``
    are queued to main-thread QObject slots. Do **not** show dialogs from a
    plain Python ``lambda`` connected to these signals — use a ``@Slot`` on a
    main-thread QObject instead.
    """

    finished = Signal(object)
    error = Signal(str)
    progress = Signal(str)
    # Command lines / process stdout from run_logged_process → process monitor.
    log_line = Signal(str)

    def __init__(self, func: Callable[..., Any], *args: Any, **kwargs: Any) -> None:
        super().__init__()
        self._func = func
        self._args = args
        self._kwargs = kwargs
        # Back-compat: existing call sites use ``worker.signals.finished``.
        self.signals = self

    @Slot()
    def run(self) -> None:
        from repeatviewer.external_process import process_log_sink

        kwargs = dict(self._kwargs)
        # Route progress callbacks through this worker's signal so UI updates
        # are queued to the main thread (never touch widgets here).
        if "progress" in kwargs and kwargs["progress"] is not None:
            kwargs["progress"] = self.progress.emit

        def _log(message: str) -> None:
            self.log_line.emit(str(message))

        try:
            with process_log_sink(_log):
                result = self._func(*self._args, **kwargs)
            self.finished.emit(result)
        except Exception as exc:  # noqa: BLE001 — surface to UI
            self.error.emit(str(exc))
