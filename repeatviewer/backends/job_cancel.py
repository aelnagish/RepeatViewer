"""Cooperative cancellation for long-running background jobs."""

from __future__ import annotations

import logging
import os
import signal
import subprocess
import sys
import threading
import time
from pathlib import Path
from typing import Iterable

from repeatviewer.external_process import log_command, windows_hidden_process_kwargs

process_logger = logging.getLogger("repeatviewer.commands")


class JobCancelledError(Exception):
    """Raised when the user stops a running job."""


class JobCancellation:
    """Thread-safe cancel flag shared between UI and worker threads."""

    def __init__(self) -> None:
        self._event = threading.Event()
        self._lock = threading.Lock()
        self._processes: list[subprocess.Popen] = []

    def reset(self) -> None:
        self._event.clear()
        with self._lock:
            self._processes.clear()

    def request(self) -> None:
        """Mark cancelled and kill registered subprocesses without blocking the UI."""
        self._event.set()
        self.kill_processes()

    @property
    def requested(self) -> bool:
        return self._event.is_set()

    def check(self) -> None:
        if self._event.is_set():
            raise JobCancelledError("Job stopped by user")

    def register_process(self, proc: subprocess.Popen) -> None:
        with self._lock:
            self._processes.append(proc)
        # Stop may have been clicked just before register — kill immediately.
        if self._event.is_set():
            terminate_process(proc, wait_s=0.35)

    def unregister_process(self, proc: subprocess.Popen) -> None:
        with self._lock:
            try:
                self._processes.remove(proc)
            except ValueError:
                pass

    def kill_processes(self) -> None:
        """Signal all registered processes; finish waiting off the caller thread."""
        with self._lock:
            procs = list(self._processes)
            self._processes.clear()
        if not procs:
            return

        def _finish() -> None:
            for proc in procs:
                terminate_process(proc, wait_s=0.35)

        # SIGTERM/SIGKILL wait must not freeze the GUI thread when Stop is clicked.
        if threading.current_thread() is threading.main_thread():
            for proc in procs:
                _signal_stop(proc)
            threading.Thread(target=_finish, daemon=True, name="job-cancel-kill").start()
        else:
            _finish()


def check_cancelled(cancel: JobCancellation | None) -> None:
    if cancel is not None:
        cancel.check()


def _signal_stop(proc: subprocess.Popen) -> None:
    """Send a stop signal without waiting for the process to exit."""
    if proc.poll() is not None:
        return
    try:
        if sys.platform != "win32":
            try:
                os.killpg(proc.pid, signal.SIGTERM)
            except (ProcessLookupError, PermissionError, OSError):
                proc.terminate()
        else:
            proc.terminate()
    except OSError:
        pass


def terminate_process(proc: subprocess.Popen, *, wait_s: float = 0.35) -> None:
    """Stop a subprocess, including child workers when started in a new session."""
    if proc.poll() is not None:
        return
    _signal_stop(proc)
    deadline = time.monotonic() + max(0.05, float(wait_s))
    while proc.poll() is None and time.monotonic() < deadline:
        time.sleep(0.02)
    if proc.poll() is None:
        try:
            if sys.platform != "win32":
                try:
                    os.killpg(proc.pid, signal.SIGKILL)
                except (ProcessLookupError, PermissionError, OSError):
                    proc.kill()
            else:
                proc.kill()
        except OSError:
            pass
        try:
            proc.wait(timeout=0.5)
        except subprocess.TimeoutExpired:
            pass


def run_subprocess_with_cancel(
    cmd: Iterable[str],
    *,
    cancel: JobCancellation | None = None,
    env: dict[str, str] | None = None,
    cwd: str | os.PathLike[str] | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run a subprocess, polling for user cancellation and terminating it on stop."""
    from repeatviewer.conda_tools import prepare_command

    command = prepare_command(list(cmd))
    tool = Path(command[0]).name if command else "external"
    # Prefer the real tool name when wrapped as wsl … bash -lc '… exec Tool …'
    if tool in {"bash", "wsl"} and len(command) >= 2:
        tool = Path(str(list(cmd)[0])).name if list(cmd) else tool
    # Capture before spawn: drain threads do not inherit ContextVars.
    from repeatviewer.external_process import current_process_log_sink, emit_process_log

    sink = current_process_log_sink()
    log_command(tool, command, cwd=cwd)
    popen_kwargs = dict(windows_hidden_process_kwargs())
    if sys.platform != "win32":
        # Own process group so Stop can kill RepeatMasker workers (-pa), etc.
        popen_kwargs["start_new_session"] = True
    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
        cwd=str(cwd) if cwd is not None else None,
        **popen_kwargs,
    )
    started = time.monotonic()
    process_logger.info("START tool=%s pid=%s", tool, proc.pid)
    emit_process_log(f"START pid={proc.pid}", sink=sink)
    captured: dict[str, list[str]] = {"stdout": [], "stderr": []}

    def _drain(name: str, stream) -> None:
        if stream is None:
            return
        for line in iter(stream.readline, ""):
            captured[name].append(line)
            text = line.rstrip("\r\n")
            process_logger.info("%s %s | %s", tool, name, text)
            if text:
                emit_process_log(text, sink=sink)

    readers = [
        threading.Thread(target=_drain, args=("stdout", proc.stdout), daemon=True),
        threading.Thread(target=_drain, args=("stderr", proc.stderr), daemon=True),
    ]
    for reader in readers:
        reader.start()
    if cancel is not None:
        cancel.register_process(proc)
    try:
        while True:
            if cancel is not None and cancel.requested:
                terminate_process(proc)
                for reader in readers:
                    reader.join(timeout=0.5)
                process_logger.warning(
                    "CANCELLED tool=%s pid=%s exit=%s elapsed=%.3fs",
                    tool,
                    proc.pid,
                    proc.returncode,
                    time.monotonic() - started,
                )
                emit_process_log(
                    f"CANCELLED exit={proc.returncode} "
                    f"elapsed={time.monotonic() - started:.3f}s",
                    sink=sink,
                )
                raise JobCancelledError("Job stopped by user")
            if proc.poll() is not None:
                # Kill may have finished the process before we observed cancel.
                if cancel is not None and cancel.requested:
                    for reader in readers:
                        reader.join(timeout=0.5)
                    raise JobCancelledError("Job stopped by user")
                break
            time.sleep(0.05)
        for reader in readers:
            reader.join()
        stdout = "".join(captured["stdout"])
        stderr = "".join(captured["stderr"])
        result = subprocess.CompletedProcess(
            proc.args,
            proc.returncode,
            stdout,
            stderr,
        )
        elapsed = time.monotonic() - started
        process_logger.log(
            logging.INFO if proc.returncode == 0 else logging.ERROR,
            "FINISH tool=%s pid=%s exit=%s elapsed=%.3fs",
            tool,
            proc.pid,
            proc.returncode,
            elapsed,
        )
        emit_process_log(
            f"FINISH exit={proc.returncode} elapsed={elapsed:.3f}s",
            sink=sink,
        )
        if check and result.returncode != 0:
            from repeatviewer.external_process import command_text

            detail_parts = [
                f"{tool} exited with code {result.returncode}",
                f"Command: {command_text(result.args if isinstance(result.args, (list, tuple)) else [result.args])}",
            ]
            err_text = (result.stderr or "").strip()
            out_text = (result.stdout or "").strip()
            if err_text:
                detail_parts.append(err_text)
            elif out_text:
                detail_parts.append(out_text)
            raise subprocess.CalledProcessError(
                result.returncode,
                result.args,
                output=result.stdout,
                stderr="\n\n".join(detail_parts),
            )
        return result
    finally:
        if cancel is not None:
            cancel.unregister_process(proc)


CANCELLED_MESSAGE = "__CANCELLED__"
