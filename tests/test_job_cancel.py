"""Job cancellation: Stop must kill child tools quickly and cleanly."""

from __future__ import annotations

import sys
import time

import pytest

from repeatviewer.backends.job_cancel import (
    JobCancellation,
    JobCancelledError,
    run_subprocess_with_cancel,
)
from repeatviewer.external_process import run_logged_process


def test_job_cancellation_request_sets_flag() -> None:
    cancel = JobCancellation()
    assert cancel.requested is False
    cancel.request()
    assert cancel.requested is True
    with pytest.raises(JobCancelledError):
        cancel.check()


def test_job_cancellation_reset_clears_flag() -> None:
    cancel = JobCancellation()
    cancel.request()
    cancel.reset()
    assert cancel.requested is False
    cancel.check()  # must not raise


def test_run_logged_process_raises_on_cancel() -> None:
    cancel = JobCancellation()
    command = [
        sys.executable,
        "-c",
        "import time; time.sleep(30)",
    ]

    def _stop_soon() -> None:
        time.sleep(0.15)
        cancel.request()

    import threading

    threading.Thread(target=_stop_soon, daemon=True).start()
    started = time.monotonic()
    with pytest.raises(JobCancelledError):
        run_logged_process(command, tool="sleeper", check=True, cancel=cancel)
    elapsed = time.monotonic() - started
    assert elapsed < 5.0, f"cancel took too long: {elapsed:.2f}s"


def test_run_subprocess_with_cancel_raises_on_cancel() -> None:
    cancel = JobCancellation()
    command = [
        sys.executable,
        "-c",
        "import time; time.sleep(30)",
    ]

    def _stop_soon() -> None:
        time.sleep(0.15)
        cancel.request()

    import threading

    threading.Thread(target=_stop_soon, daemon=True).start()
    started = time.monotonic()
    with pytest.raises(JobCancelledError):
        run_subprocess_with_cancel(command, cancel=cancel, check=True)
    elapsed = time.monotonic() - started
    assert elapsed < 5.0, f"cancel took too long: {elapsed:.2f}s"


def test_run_logged_process_succeeds_without_cancel() -> None:
    result = run_logged_process(
        [sys.executable, "-c", "print('ok', flush=True)"],
        tool="demo",
        check=True,
        cancel=None,
    )
    assert result.returncode == 0
    assert "ok" in (result.stdout or "")


def test_register_process_kills_if_already_cancelled(tmp_path) -> None:
    """Late register after Stop must still terminate the child."""
    import subprocess

    cancel = JobCancellation()
    cancel.request()
    proc = subprocess.Popen(
        [sys.executable, "-c", "import time; time.sleep(30)"],
    )
    try:
        cancel.register_process(proc)
        deadline = time.monotonic() + 3.0
        while proc.poll() is None and time.monotonic() < deadline:
            time.sleep(0.05)
        assert proc.poll() is not None, "process should be killed after late register"
    finally:
        if proc.poll() is None:
            proc.kill()
            proc.wait(timeout=2)
