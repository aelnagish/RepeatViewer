"""CPU thread budget for BLAST and other analysis backends."""

from __future__ import annotations

import os

from PySide6.QtCore import QSettings

# Threads reserved for the RepeatViewer UI, file browser, and Qt event loop.
_APP_THREAD_RESERVE = 2


DEFAULT_BATCH_WORKERS = 5


def load_batch_workers() -> int:
    return max(1, int(QSettings().value("batch_workers", DEFAULT_BATCH_WORKERS)))


def reserved_app_threads() -> int:
    cpus = os.cpu_count() or 4
    # User override from Preferences (0 == auto).
    try:
        from repeatviewer.app_preferences import load_ui_thread_reserve

        override = load_ui_thread_reserve()
    except Exception:
        override = 0
    if override > 0:
        return max(1, min(override, cpus - 1))
    if cpus <= 2:
        return 1
    return min(_APP_THREAD_RESERVE, cpus - 1)


def analysis_thread_pool() -> int:
    cpus = os.cpu_count() or 4
    return max(1, cpus - reserved_app_threads())


def analysis_threads_per_job(
    batch_workers: int | None = None,
    *,
    parallel_jobs: int = 1,
) -> int:
    """
    Threads for one analysis process when auto-budgeting from CPUs.

    ``batch_workers`` is Max. parallel jobs from Preferences.
    ``parallel_jobs`` is how many analysis processes run at once (e.g. similarity DBs).
    """
    workers = max(1, int(batch_workers if batch_workers is not None else load_batch_workers()))
    concurrent = max(1, min(workers, int(parallel_jobs)))
    return max(1, analysis_thread_pool() // concurrent)
