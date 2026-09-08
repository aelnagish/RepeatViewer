"""Cooperative cancellation for long-running background jobs."""

from __future__ import annotations

import subprocess
import logging
import threading
import time
from typing import Iterable

from repeatviewer.external_process import log_command

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

    def unregister_process(self, proc: subprocess.Popen) -> None:
        with self._lock:
            try:
                self._processes.remove(proc)
            except ValueError:
                pass

    def kill_processes(self) -> None:
        with self._lock:
            procs = list(self._processes)
        for proc in procs:
            terminate_process(proc)
        with self._lock:
            self._processes.clear()


def check_cancelled(cancel: JobCancellation | None) -> None:
    if cancel is not None:
        cancel.check()


def terminate_process(proc: subprocess.Popen, *, wait_s: float = 3.0) -> None:
    if proc.poll() is not None:
        return
    try:
        proc.terminate()
    except OSError:
        pass
    deadline = time.monotonic() + wait_s
    while proc.poll() is None and time.monotonic() < deadline:
        time.sleep(0.05)
    if proc.poll() is None:
        try:
            proc.kill()
        except OSError:
            pass
        try:
            proc.wait(timeout=2)
        except subprocess.TimeoutExpired:
            pass


def run_subprocess_with_cancel(
    cmd: Iterable[str],
    *,
    cancel: JobCancellation | None = None,
    env: dict[str, str] | None = None,
    check: bool = True,
) -> subprocess.CompletedProcess[str]:
    """Run a subprocess, polling for user cancellation and terminating it on stop."""
    command = list(cmd)
    tool = command[0].rsplit("/", 1)[-1] if command else "external"
    log_command(tool, command)
    proc = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        env=env,
    )
    started = time.monotonic()
    process_logger.info("START tool=%s pid=%s", tool, proc.pid)
    captured: dict[str, list[str]] = {"stdout": [], "stderr": []}

    def _drain(name: str, stream) -> None:
        if stream is None:
            return
        for line in iter(stream.readline, ""):
            captured[name].append(line)
            process_logger.info("%s %s | %s", tool, name, line.rstrip("\r\n"))

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
                    reader.join()
                process_logger.warning(
                    "CANCELLED tool=%s pid=%s exit=%s elapsed=%.3fs",
                    tool,
                    proc.pid,
                    proc.returncode,
                    time.monotonic() - started,
                )
                raise JobCancelledError("Job stopped by user")
            if proc.poll() is not None:
                break
            time.sleep(0.25)
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
        process_logger.log(
            logging.INFO if proc.returncode == 0 else logging.ERROR,
            "FINISH tool=%s pid=%s exit=%s elapsed=%.3fs",
            tool,
            proc.pid,
            proc.returncode,
            time.monotonic() - started,
        )
        if check and result.returncode != 0:
            raise subprocess.CalledProcessError(
                result.returncode,
                result.args,
                output=result.stdout,
                stderr=result.stderr,
            )
        return result
    finally:
        if cancel is not None:
            cancel.unregister_process(proc)


CANCELLED_MESSAGE = "__CANCELLED__"
