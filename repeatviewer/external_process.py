"""Consistent logging for commands launched outside the Python process."""

from __future__ import annotations

import contextvars
import logging
import os
import shlex
import subprocess
import sys
import threading
import time
from collections.abc import Callable, Iterable, Iterator, Sequence
from contextlib import contextmanager
from pathlib import Path

logger = logging.getLogger("repeatviewer.commands")

# Packed Windows builds use console=False (no parent console). Child console
# tools such as blastn.exe / mafft.exe / wsl.exe then flash a new cmd window.
# Unpacked `python` already has a console, so children inherit it and no extra
# window appears. CREATE_NO_WINDOW keeps both cases quiet.
_CREATE_NO_WINDOW = getattr(subprocess, "CREATE_NO_WINDOW", 0x08000000)


def windows_hidden_process_kwargs() -> dict:
    """Kwargs so Windows GUI/frozen launches do not flash a console."""
    if sys.platform != "win32":
        return {}
    startupinfo = subprocess.STARTUPINFO()
    startupinfo.dwFlags |= subprocess.STARTF_USESHOWWINDOW
    startupinfo.wShowWindow = getattr(subprocess, "SW_HIDE", 0)
    return {"startupinfo": startupinfo, "creationflags": _CREATE_NO_WINDOW}


_ProcessLogCallback = Callable[[str], None]
_process_log_sink: contextvars.ContextVar[_ProcessLogCallback | None] = contextvars.ContextVar(
    "repeatviewer_process_log_sink",
    default=None,
)


@contextmanager
def process_log_sink(callback: _ProcessLogCallback | None) -> Iterator[None]:
    """Temporarily forward command/process lines to ``callback`` (UI monitor)."""
    token = _process_log_sink.set(callback)
    try:
        yield
    finally:
        _process_log_sink.reset(token)


def emit_process_log(message: str, *, sink: _ProcessLogCallback | None = None) -> None:
    callback = sink if sink is not None else _process_log_sink.get()
    if callback is None:
        return
    try:
        callback(str(message))
    except Exception:  # noqa: BLE001 - never break a running tool
        pass


def current_process_log_sink() -> _ProcessLogCallback | None:
    return _process_log_sink.get()


def command_text(command: Sequence[object]) -> str:
    """Return a shell-readable command without invoking a shell."""
    return shlex.join(str(part) for part in command)


def log_command(
    tool: str,
    command: Sequence[object],
    *,
    inputs: Iterable[str | os.PathLike[str]] = (),
    cwd: str | os.PathLike[str] | None = None,
) -> None:
    input_text = ", ".join(str(Path(item)) for item in inputs)
    cmd = command_text(command)
    logger.info(
        "RUN tool=%s inputs=[%s] cwd=%s command=%s",
        tool,
        input_text,
        str(cwd) if cwd is not None else str(Path.cwd()),
        cmd,
    )
    if _process_log_sink.get() is not None:
        emit_process_log("Command:")
        emit_process_log(cmd)
        if cwd is not None:
            emit_process_log(f"Working directory: {cwd}")
        if input_text:
            emit_process_log(f"Inputs: {input_text}")
        emit_process_log("")
        emit_process_log("—— program output ——")


def log_process_output(
    tool: str,
    *,
    stdout: str | bytes | None = None,
    stderr: str | bytes | None = None,
    stdout_is_result: bool = False,
    output_path: str | os.PathLike[str] | None = None,
) -> None:
    """Copy terminal-style process output into the application log."""
    sink = current_process_log_sink()
    for stream_name, content in (("stdout", stdout), ("stderr", stderr)):
        if not content:
            continue
        if stream_name == "stdout" and stdout_is_result:
            continue
        if isinstance(content, bytes):
            content = content.decode("utf-8", errors="replace")
        for line in content.rstrip().splitlines():
            logger.info("%s %s | %s", tool, stream_name, line)
            emit_process_log(line, sink=sink)
    if output_path is not None:
        logger.info("%s output | %s", tool, output_path)
        emit_process_log(f"output: {output_path}", sink=sink)


def run_logged_process(
    command: Sequence[object] | str,
    *,
    tool: str,
    inputs: Iterable[str | os.PathLike[str]] = (),
    cwd: str | os.PathLike[str] | None = None,
    env: dict[str, str] | None = None,
    check: bool = False,
    timeout: float | None = None,
    shell: bool = False,
    stdout_is_result: bool = False,
    output_path: str | os.PathLike[str] | None = None,
    progress_file: str | os.PathLike[str] | None = None,
) -> subprocess.CompletedProcess[str]:
    """Run a command while streaming its terminal output into the log."""
    if not isinstance(command, str):
        from repeatviewer.conda_tools import prepare_command

        command = prepare_command(command)
    display_command = (
        ["shell", "-c", command] if isinstance(command, str) else list(command)
    )
    # Capture before spawn: drain threads do not inherit ContextVars.
    sink = current_process_log_sink()
    log_command(tool, display_command, inputs=inputs, cwd=cwd)
    started = time.monotonic()
    try:
        proc = subprocess.Popen(
            command,
            cwd=str(cwd) if cwd is not None else None,
            env=env,
            shell=shell,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
            bufsize=1,
            **windows_hidden_process_kwargs(),
        )
    except OSError:
        logger.exception("LAUNCH FAILED tool=%s", tool)
        emit_process_log(f"LAUNCH FAILED tool={tool}", sink=sink)
        raise
    logger.info("START tool=%s pid=%s", tool, proc.pid)
    emit_process_log(f"START pid={proc.pid}", sink=sink)
    captured: dict[str, list[str]] = {"stdout": [], "stderr": []}

    def _drain(name: str, stream, *, write_lines: bool) -> None:
        if stream is None:
            return
        for line in iter(stream.readline, ""):
            captured[name].append(line)
            if write_lines:
                text = line.rstrip("\r\n")
                logger.info("%s %s | %s", tool, name, text)
                if text:
                    emit_process_log(text, sink=sink)
        stream.close()

    stdout_thread = threading.Thread(
        target=_drain,
        args=("stdout", proc.stdout),
        kwargs={"write_lines": not stdout_is_result},
        daemon=True,
    )
    stderr_thread = threading.Thread(
        target=_drain,
        args=("stderr", proc.stderr),
        kwargs={"write_lines": True},
        daemon=True,
    )
    stdout_thread.start()
    stderr_thread.start()
    progress_stop = threading.Event()

    def _tail_progress() -> None:
        if progress_file is None:
            return
        path = Path(progress_file)
        position = 0
        while True:
            try:
                has_unread = path.exists() and path.stat().st_size > position
            except OSError:
                has_unread = False
            if progress_stop.is_set() and not has_unread:
                break
            if path.exists():
                try:
                    with path.open(encoding="utf-8", errors="replace") as handle:
                        handle.seek(position)
                        for line in handle:
                            text = line.rstrip("\r\n")
                            logger.info("%s progress | %s", tool, text)
                            if text:
                                emit_process_log(text, sink=sink)
                        position = handle.tell()
                except OSError:
                    pass
            progress_stop.wait(0.1)

    progress_thread = threading.Thread(target=_tail_progress, daemon=True)
    progress_thread.start()
    try:
        returncode = proc.wait(timeout=timeout)
    except subprocess.TimeoutExpired:
        proc.kill()
        proc.wait()
        logger.error("TIMEOUT tool=%s pid=%s timeout=%ss", tool, proc.pid, timeout)
        emit_process_log(
            f"TIMEOUT tool={tool} pid={proc.pid} timeout={timeout}s",
            sink=sink,
        )
        raise
    finally:
        stdout_thread.join()
        stderr_thread.join()
        progress_stop.set()
        progress_thread.join()

    stdout = "".join(captured["stdout"])
    stderr = "".join(captured["stderr"])
    elapsed = time.monotonic() - started
    level = logging.INFO if returncode == 0 else logging.ERROR
    logger.log(
        level,
        "FINISH tool=%s pid=%s exit=%s elapsed=%.3fs",
        tool,
        proc.pid,
        returncode,
        elapsed,
    )
    emit_process_log(
        f"FINISH exit={returncode} elapsed={elapsed:.3f}s",
        sink=sink,
    )
    if stdout_is_result and stdout.strip():
        emit_process_log(
            f"(stdout captured as result, {len(stdout)} characters not shown)",
            sink=sink,
        )
    if output_path is not None:
        logger.info("%s output | %s", tool, output_path)
        emit_process_log(f"output: {output_path}", sink=sink)
    result = subprocess.CompletedProcess(display_command, returncode, stdout, stderr)
    if check and returncode != 0:
        raise subprocess.CalledProcessError(
            returncode, display_command, output=stdout, stderr=stderr
        )
    return result
