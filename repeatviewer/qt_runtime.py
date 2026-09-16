"""Qt / WebEngine environment setup — call before importing PySide6 widgets."""

from __future__ import annotations

import os
import sys
import threading
from ctypes import c_char_p, c_void_p, cdll
from pathlib import Path

from repeatviewer.branding import APP_DISPLAY_NAME

_MACOS_STDERR_SKIP = (
    b"IMKCFRunLoopWakeUpReliable",
    b"IMKClient subclass",
    b"IMKInputSession subclass",
    b"error messaging the mach port",
    b"TSMSendMessageToUIServer",
)


def configure_qt_runtime() -> None:
    """
    Reduce noisy stderr from Qt WebEngine and macOS pointer synthesis.

    Set REPEATVIEWER_QUIET_QT=0 to disable these defaults when debugging.
    On macOS, harmless TSM/IMK stderr lines are filtered by default; set
    REPEATVIEWER_FILTER_STDERR_FD=0 if the HTML preview (TEAid / genome map) misbehaves.
    """
    # Windows/Python 3.12: pandas → python-dateutil → six.moves can trip
    # Shiboken's import hook if PySide6 loaded first
    # (AttributeError: '_SixMetaPathImporter' object has no attribute '_path').
    # Preload pandas before any PySide6 import.
    if sys.platform == "win32":
        try:
            import pandas  # noqa: F401
        except ImportError:
            pass

    _configure_macos_app_name()
    # Must run even when Qt quieting is disabled -- it prevents hard crashes.
    _configure_matplotlib_backend()
    install_crash_diagnostics()
    if os.environ.get("REPEATVIEWER_QUIET_QT", "1").lower() in {"0", "false", "no"}:
        return

    _configure_qt_logging()
    _configure_macos_pointer()
    _configure_webengine()
    _install_stderr_filter()


def install_crash_diagnostics() -> None:
    """Opt-in diagnostics for segfaults and cross-thread Qt warnings.

    Enable with REPEATVIEWER_DEBUG_CRASH=1. Two things are installed:

    * ``faulthandler`` -- on a segfault, dumps the Python stack of every thread
      to stderr instead of the process vanishing silently.
    * A Qt message handler that, on "Timers cannot be started from another
      thread", prints the offending thread name and Python stack. If the stack
      shows no repeatviewer frames the call came from inside Qt (most likely the
      macOS accessibility bridge) and is not our bug.
    """
    if os.environ.get("REPEATVIEWER_DEBUG_CRASH", "").strip().lower() not in {
        "1",
        "true",
        "yes",
        "on",
    }:
        return

    import faulthandler
    import threading
    import traceback

    faulthandler.enable()

    try:
        from PySide6.QtCore import qInstallMessageHandler
    except ImportError:
        return

    def _handler(mode, context, message) -> None:
        text = str(message)
        sys.stderr.write(f"[qt] {text}\n")
        if "Timers cannot be started" in text or "another thread" in text:
            sys.stderr.write(
                f"[qt]   current thread: {threading.current_thread().name}\n"
            )
            sys.stderr.write("[qt]   Python stack:\n")
            for line in traceback.format_stack():
                sys.stderr.write("[qt]     " + line.rstrip() + "\n")
            sys.stderr.flush()
        del mode, context

    qInstallMessageHandler(_handler)
    sys.stderr.write(
        "[repeatviewer] crash diagnostics on (faulthandler + Qt message handler)\n"
    )


def _configure_matplotlib_backend() -> None:
    """
    Force matplotlib's headless Agg backend process-wide.

    Plot rendering runs on background worker threads. Because PySide6 is
    installed, matplotlib would otherwise auto-select the interactive QtAgg
    backend, which builds QObjects/QTimers on a non-GUI thread and segfaults
    ("QObject::startTimer: Timers cannot be started from another thread").

    Setting MPLBACKEND here covers every ``import matplotlib.pyplot`` in the
    codebase, so individual call sites cannot reintroduce the crash.
    Set REPEATVIEWER_MPL_BACKEND to override when debugging plots interactively.
    """
    os.environ.setdefault(
        "MPLBACKEND", os.environ.get("REPEATVIEWER_MPL_BACKEND", "Agg")
    )


def install_macos_stderr_fd_filter() -> None:
    """Filter TSM/IMK noise written directly to stderr FD (call after QApplication exists)."""
    if sys.platform != "darwin":
        return
    if os.environ.get("REPEATVIEWER_FILTER_STDERR_FD", "1").lower() in {"1", "true", "yes"}:
        _install_macos_stderr_fd_filter()


def _configure_macos_app_name() -> None:
    """Show RepeatViewer in the macOS menu bar instead of Python when launched via interpreter."""
    if sys.platform != "darwin":
        return
    apply_macos_app_name()


def apply_macos_app_name() -> None:
    """Best-effort rename of the macOS application menu for interpreter launches."""
    if sys.platform != "darwin":
        return
    try:
        libc = cdll.LoadLibrary("/usr/lib/libc.dylib")
        if hasattr(libc, "setprogname"):
            libc.setprogname(c_char_p(APP_DISPLAY_NAME.encode("utf-8")))
    except Exception:
        pass
    try:
        import Foundation  # type: ignore[import-not-found]

        Foundation.NSProcessInfo.processInfo().setProcessName_(APP_DISPLAY_NAME)
    except Exception:
        pass
    try:
        objc = cdll.LoadLibrary("/usr/lib/libobjc.A.dylib")
        objc.objc_getClass.restype = c_void_p
        objc.sel_registerName.restype = c_void_p
        objc.objc_msgSend.restype = c_void_p

        ns_process_info = objc.objc_getClass(b"NSProcessInfo")
        ns_string = objc.objc_getClass(b"NSString")
        process_info = objc.objc_msgSend(
            c_void_p(ns_process_info),
            c_void_p(objc.sel_registerName(b"processInfo")),
        )
        app_name = objc.objc_msgSend(
            c_void_p(ns_string),
            c_void_p(objc.sel_registerName(b"stringWithUTF8String:")),
            c_char_p(APP_DISPLAY_NAME.encode("utf-8")),
        )
        objc.objc_msgSend(
            c_void_p(process_info),
            c_void_p(objc.sel_registerName(b"setProcessName:")),
            c_void_p(app_name),
        )
    except Exception:
        pass

    # The bold application-menu title on macOS comes from CFBundleName, not the
    # process name. For an interpreter launch (no .app bundle) it defaults to
    # "python3.x"; override the main bundle's info dictionary so it reads as the
    # app name instead. Must run before the Cocoa menu bar is built.
    _set_macos_bundle_name()


def _set_macos_bundle_name() -> None:
    """Rewrite CFBundleName in the running process's main bundle (best effort)."""
    if sys.platform != "darwin":
        return
    try:
        objc = cdll.LoadLibrary("/usr/lib/libobjc.A.dylib")
        objc.objc_getClass.restype = c_void_p
        objc.objc_getClass.argtypes = [c_char_p]
        objc.sel_registerName.restype = c_void_p
        objc.sel_registerName.argtypes = [c_char_p]
        objc.objc_msgSend.restype = c_void_p

        def msg(receiver, selector, *args, argtypes=()):
            objc.objc_msgSend.argtypes = [c_void_p, c_void_p, *argtypes]
            return objc.objc_msgSend(
                c_void_p(receiver),
                c_void_p(objc.sel_registerName(selector)),
                *args,
            )

        ns_bundle = objc.objc_getClass(b"NSBundle")
        ns_string = objc.objc_getClass(b"NSString")
        if not ns_bundle or not ns_string:
            return
        main_bundle = msg(ns_bundle, b"mainBundle")
        if not main_bundle:
            return
        info = msg(main_bundle, b"infoDictionary")
        if not info:
            return

        def ns_str(text: str):
            return msg(
                ns_string,
                b"stringWithUTF8String:",
                c_char_p(text.encode("utf-8")),
                argtypes=(c_char_p,),
            )

        value = ns_str(APP_DISPLAY_NAME)
        key = ns_str("CFBundleName")
        if value and key:
            msg(
                info,
                b"setObject:forKey:",
                c_void_p(value),
                c_void_p(key),
                argtypes=(c_void_p, c_void_p),
            )
    except Exception:
        pass


def _configure_qt_logging() -> None:
    rules = os.environ.get("QT_LOGGING_RULES", "").strip()
    for rule in ("qt.pointer.dispatch=false",):
        if rule.split("=")[0] not in rules:
            rules = f"{rules};{rule}" if rules else rule
    os.environ["QT_LOGGING_RULES"] = rules


def _configure_macos_pointer() -> None:
    if sys.platform == "darwin":
        os.environ.setdefault("QT_MAC_DISABLE_TOUCH_SYNTHESIS", "1")


def _configure_webengine() -> None:
    flags = os.environ.get("QTWEBENGINE_CHROMIUM_FLAGS", "").strip()
    for flag in ("--log-level=3",):
        if flag not in flags:
            flags = f"{flags} {flag}".strip()
    if sys.platform == "darwin":
        for flag in (
            "--disable-gpu",
            "--disable-gpu-compositing",
            "--disable-gpu-rasterization",
            "--disable-features=Metal,Vulkan,WebGPU,VizDisplayCompositor",
            "--no-sandbox",
        ):
            if flag not in flags:
                flags = f"{flags} {flag}".strip()
    os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = flags

    if os.environ.get("QTWEBENGINE_DICTIONARIES_PATH"):
        return

    try:
        import PySide6
    except ImportError:
        return

    pyside_root = Path(PySide6.__file__).resolve().parent
    dict_dir = pyside_root / "qtwebengine_dictionaries"
    dict_dir.mkdir(exist_ok=True)
    os.environ["QTWEBENGINE_DICTIONARIES_PATH"] = str(dict_dir)


class _QtNoiseStderrFilter:
    """Drop known harmless Qt/WebEngine stderr lines."""

    _SKIP_FRAGMENTS = (
        "qt.pointer.dispatch:",
        "Path override failed for key base::DIR_APP_DICTIONARIES",
        "Path override failed for key",
        "IMKCFRunLoopWakeUpReliable",
        "IMKClient subclass",
        "IMKInputSession subclass",
        "error messaging the mach port",
        "TSMSendMessageToUIServer",
    )

    def __init__(self, stream):
        self._stream = stream

    def write(self, text: str) -> int:
        if text and not any(fragment in text for fragment in self._SKIP_FRAGMENTS):
            return self._stream.write(text)
        return len(text)

    def flush(self) -> None:
        self._stream.flush()

    def fileno(self):
        return self._stream.fileno()

    def isatty(self):
        return self._stream.isatty()


def _install_macos_stderr_fd_filter() -> None:
    """
    Filter harmless macOS Input Method Kit (IMK) stderr noise.

    Those messages are written by Apple frameworks directly to file descriptor 2,
    so wrapping sys.stderr alone is not enough.
    """
    if sys.platform != "darwin":
        return
    if os.environ.get("REPEATVIEWER_FILTER_STDERR", "1").lower() in {"0", "false", "no"}:
        return
    if getattr(_install_macos_stderr_fd_filter, "_installed", False):
        return

    def _should_drop(line: bytes) -> bool:
        return any(fragment in line for fragment in _MACOS_STDERR_SKIP)

    def _drain_stderr(read_fd: int, saved_stderr_fd: int) -> None:
        buffer = b""
        while True:
            chunk = os.read(read_fd, 4096)
            if not chunk:
                break
            buffer += chunk
            while b"\n" in buffer:
                line, buffer = buffer.split(b"\n", 1)
                payload = line + b"\n"
                if not _should_drop(payload):
                    os.write(saved_stderr_fd, payload)
        if buffer and not _should_drop(buffer):
            os.write(saved_stderr_fd, buffer)

    try:
        read_fd, write_fd = os.pipe()
        saved_stderr_fd = os.dup(2)
        os.dup2(write_fd, 2)
        os.close(write_fd)
        thread = threading.Thread(
            target=_drain_stderr,
            args=(read_fd, saved_stderr_fd),
            name="repeatviewer-stderr-filter",
            daemon=True,
        )
        thread.start()
        _install_macos_stderr_fd_filter._installed = True
    except OSError:
        pass


def _install_stderr_filter() -> None:
    if os.environ.get("REPEATVIEWER_FILTER_STDERR", "1").lower() in {"0", "false", "no"}:
        return
    if isinstance(sys.stderr, _QtNoiseStderrFilter):
        return
    sys.stderr = _QtNoiseStderrFilter(sys.stderr)
