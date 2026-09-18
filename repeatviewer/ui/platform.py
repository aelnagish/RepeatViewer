"""OS integration helpers."""

from __future__ import annotations

import platform
import subprocess
import sys
from pathlib import Path

from repeatviewer.external_process import log_command


def hide_from_macos_dock() -> None:
    """Hide this process from the macOS Dock (best effort).

    Used for frozen-app helper re-entries such as Flexidot, which re-launch the
    same ``.app`` binary and would otherwise show a second Dock icon until exit.
    """
    if sys.platform != "darwin":
        return
    # NSApplicationActivationPolicyProhibited = 2
    try:
        from ctypes import c_bool, c_char_p, c_long, c_void_p, cdll

        # Ensure AppKit classes are registered before objc_getClass.
        cdll.LoadLibrary("/System/Library/Frameworks/AppKit.framework/AppKit")
        objc = cdll.LoadLibrary("/usr/lib/libobjc.A.dylib")
        objc.objc_getClass.restype = c_void_p
        objc.objc_getClass.argtypes = [c_char_p]
        objc.sel_registerName.restype = c_void_p
        objc.sel_registerName.argtypes = [c_char_p]

        ns_app = objc.objc_getClass(b"NSApplication")
        if not ns_app:
            raise RuntimeError("NSApplication unavailable")
        shared_sel = objc.sel_registerName(b"sharedApplication")
        policy_sel = objc.sel_registerName(b"setActivationPolicy:")

        objc.objc_msgSend.restype = c_void_p
        objc.objc_msgSend.argtypes = [c_void_p, c_void_p]
        app = objc.objc_msgSend(c_void_p(ns_app), c_void_p(shared_sel))
        if not app:
            raise RuntimeError("sharedApplication unavailable")

        objc.objc_msgSend.restype = c_bool
        objc.objc_msgSend.argtypes = [c_void_p, c_void_p, c_long]
        if not objc.objc_msgSend(c_void_p(app), c_void_p(policy_sel), c_long(2)):
            raise RuntimeError("setActivationPolicy failed")
        return
    except Exception:
        pass

    # Fallback for older macOS / if AppKit messaging fails.
    try:
        from ctypes import POINTER, Structure, byref, c_uint32, cdll

        class ProcessSerialNumber(Structure):
            _fields_ = [("highLongOfPSN", c_uint32), ("lowLongOfPSN", c_uint32)]

        lib = cdll.LoadLibrary(
            "/System/Library/Frameworks/ApplicationServices.framework/ApplicationServices"
        )
        psn = ProcessSerialNumber(0, 2)  # kCurrentProcess
        lib.TransformProcessType.argtypes = [POINTER(ProcessSerialNumber), c_uint32]
        # kProcessTransformToBackgroundApplication = 2
        lib.TransformProcessType(byref(psn), c_uint32(2))
    except Exception:
        pass


def reveal_in_file_manager(path: Path) -> None:
    """Show a file or folder in the system file manager."""
    path = path.expanduser().resolve()
    system = platform.system()
    try:
        if system == "Darwin":
            command = ["open", "-R", str(path)]
        elif system == "Windows":
            command = ["explorer", "/select,", str(path)]
        elif path.is_file():
            command = ["xdg-open", str(path.parent)]
        else:
            command = ["xdg-open", str(path)]
        log_command("file-manager", command, inputs=[path])
        subprocess.run(command, check=False)
    except OSError:
        pass
