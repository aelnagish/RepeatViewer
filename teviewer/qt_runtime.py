"""Qt / WebEngine environment setup — call before importing PySide6 widgets."""

from __future__ import annotations

import os
import sys
import threading
from ctypes import c_char_p, c_void_p, cdll
from pathlib import Path

from teviewer.branding import APP_DISPLAY_NAME

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

    Set TEVIEWER_QUIET_QT=0 to disable these defaults when debugging.
    On macOS, harmless TSM/IMK stderr lines are filtered by default; set
    TEVIEWER_FILTER_STDERR_FD=0 if the HTML preview (TEAid / genome map) misbehaves.
    """
    _configure_macos_app_name()
    if os.environ.get("TEVIEWER_QUIET_QT", "1").lower() in {"0", "false", "no"}:
        return

    _configure_qt_logging()
    _configure_macos_pointer()
    _configure_webengine()
    if sys.platform == "darwin":
        os.environ.setdefault("TEVIEWER_FILTER_STDERR_FD", "1")
    if os.environ.get("TEVIEWER_FILTER_STDERR_FD", "0").lower() in {"1", "true", "yes"}:
        _install_macos_stderr_fd_filter()
    _install_stderr_filter()


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
        for flag in ("--disable-gpu", "--no-sandbox"):
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
    if os.environ.get("TEVIEWER_FILTER_STDERR", "1").lower() in {"0", "false", "no"}:
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
            name="teviewer-stderr-filter",
            daemon=True,
        )
        thread.start()
        _install_macos_stderr_fd_filter._installed = True
    except OSError:
        pass


def _install_stderr_filter() -> None:
    if os.environ.get("TEVIEWER_FILTER_STDERR", "1").lower() in {"0", "false", "no"}:
        return
    if isinstance(sys.stderr, _QtNoiseStderrFilter):
        return
    sys.stderr = _QtNoiseStderrFilter(sys.stderr)
