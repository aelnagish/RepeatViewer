"""Detect whether Qt WebEngine can be used without crashing the process.

On some macOS + PySide6 setups, creating ``QWebEngineView`` (or loading Plotly
HTML) prints ``GPUInfo not initialized on GpuInfoUpdate`` and then segfaults.
That cannot be caught inside the main process, so we probe in a short-lived
child first and fall back to the system browser when the probe dies.
"""

from __future__ import annotations

import logging
import os
import subprocess
import sys
import tempfile
import threading
from pathlib import Path

logger = logging.getLogger(__name__)

_PROBE_RESULT: bool | None = None
_PROBE_LOCK = threading.Lock()

_PROBE_SCRIPT = r"""
import os
import sys

# Match RepeatViewer macOS Chromium flags before any Qt import.
flags = os.environ.get("QTWEBENGINE_CHROMIUM_FLAGS", "").strip()
extra = (
    "--log-level=3 "
    "--disable-gpu "
    "--disable-gpu-compositing "
    "--disable-gpu-rasterization "
    "--disable-features=Metal,Vulkan,WebGPU,VizDisplayCompositor "
    "--in-process-gpu "
    "--no-sandbox"
)
for part in extra.split():
    key = part.split("=")[0]
    if key not in flags:
        flags = f"{flags} {part}".strip()
os.environ["QTWEBENGINE_CHROMIUM_FLAGS"] = flags

from PySide6.QtCore import Qt, QTimer, QUrl
from PySide6.QtWidgets import QApplication
from PySide6.QtWebEngineCore import QWebEngineSettings
from PySide6.QtWebEngineWidgets import QWebEngineView

QApplication.setAttribute(Qt.ApplicationAttribute.AA_ShareOpenGLContexts, True)
app = QApplication(sys.argv)
view = QWebEngineView()
settings = view.settings()
settings.setAttribute(QWebEngineSettings.WebAttribute.JavascriptEnabled, True)
settings.setAttribute(QWebEngineSettings.WebAttribute.Accelerated2dCanvasEnabled, False)
if hasattr(QWebEngineSettings.WebAttribute, "WebGLEnabled"):
    settings.setAttribute(QWebEngineSettings.WebAttribute.WebGLEnabled, False)
view.setHtml("<!DOCTYPE html><html><body>repeatviewer-webengine-ok</body></html>")

def _done(ok: bool = True) -> None:
    app.quit()

view.loadFinished.connect(lambda ok: QTimer.singleShot(50, _done))
QTimer.singleShot(4000, _done)
code = app.exec()
sys.exit(0 if code == 0 else 1)
"""


def _cache_path() -> Path:
    root = Path(tempfile.gettempdir()) / "repeatviewer_webengine_probe"
    root.mkdir(parents=True, exist_ok=True)
    try:
        import PySide6

        version = getattr(PySide6, "__version__", "unknown")
    except Exception:
        version = "unknown"
    py = f"{sys.version_info.major}.{sys.version_info.minor}"
    return root / f"safe-{sys.platform}-py{py}-pyside{version}.txt"


def _read_cached() -> bool | None:
    path = _cache_path()
    try:
        text = path.read_text(encoding="utf-8").strip()
    except OSError:
        return None
    if text == "1":
        return True
    if text == "0":
        return False
    return None


def _write_cached(ok: bool) -> None:
    try:
        _cache_path().write_text("1" if ok else "0", encoding="utf-8")
    except OSError:
        pass


def clear_webengine_probe_cache() -> None:
    """Forget a previous probe result (used after Qt upgrades / forced retest)."""
    global _PROBE_RESULT
    _PROBE_RESULT = None
    try:
        _cache_path().unlink(missing_ok=True)
    except OSError:
        pass


def mark_webengine_safe() -> None:
    """Record that WebEngine loaded successfully in the main GUI process."""
    global _PROBE_RESULT
    with _PROBE_LOCK:
        _PROBE_RESULT = True
        _write_cached(True)


def webengine_is_known_unsafe() -> bool:
    """True when a probe or prior run marked WebEngine unsafe on this machine."""
    with _PROBE_LOCK:
        if _PROBE_RESULT is False:
            return True
    cached = _read_cached()
    return cached is False


def webengine_is_safe(*, force_probe: bool = False) -> bool:
    """Return False when embedded HTML preview would likely segfault this machine."""
    global _PROBE_RESULT
    with _PROBE_LOCK:
        if _PROBE_RESULT is not None and not force_probe:
            return _PROBE_RESULT

        flag = os.environ.get("REPEATVIEWER_EMBED_HTML", "").strip().lower()
        if flag in {"0", "false", "no", "off"}:
            _PROBE_RESULT = False
            return False
        if flag in {"1", "true", "yes", "on"}:
            # Explicit opt-in: skip probe (user accepts crash risk).
            _PROBE_RESULT = True
            return True

        if sys.platform != "darwin":
            _PROBE_RESULT = True
            return True

        if not force_probe:
            cached = _read_cached()
            if cached is not None:
                _PROBE_RESULT = cached
                return cached

        ok = _run_webengine_probe()
        _PROBE_RESULT = ok
        _write_cached(ok)
        if not ok:
            logger.warning(
                "Qt WebEngine probe failed on this Mac; HTML files will open in the "
                "system browser. Set REPEATVIEWER_EMBED_HTML=1 to force in-app preview, "
                "or REPEATVIEWER_EMBED_HTML=0 to silence this and always use the browser."
            )
        return ok


def warmup_webengine_probe() -> None:
    """Start the macOS safety probe in the background after the UI is up."""
    if sys.platform != "darwin":
        return
    flag = os.environ.get("REPEATVIEWER_EMBED_HTML", "").strip().lower()
    if flag in {"0", "false", "no", "off", "1", "true", "yes", "on"}:
        return
    if _PROBE_RESULT is not None or _read_cached() is not None:
        return
    threading.Thread(
        target=webengine_is_safe,
        name="repeatviewer-webengine-probe",
        daemon=True,
    ).start()


def _run_webengine_probe() -> bool:
    """Spawn a child that creates a WebEngine view; treat crashes as unsafe."""
    env = os.environ.copy()
    # Keep the probe quiet and headless-ish.
    env.setdefault("QT_LOGGING_RULES", "*=false")
    env["REPEATVIEWER_QUIET_QT"] = "1"
    try:
        result = subprocess.run(
            [sys.executable, "-c", _PROBE_SCRIPT],
            env=env,
            capture_output=True,
            text=True,
            timeout=12,
            check=False,
        )
    except subprocess.TimeoutExpired:
        logger.warning("Qt WebEngine probe timed out")
        return False
    except OSError as exc:
        logger.warning("Qt WebEngine probe could not start: %s", exc)
        return False

    if result.returncode != 0:
        # Negative codes / 139 / 245 etc. usually mean SIGSEGV in the child.
        logger.warning(
            "Qt WebEngine probe exited with code %s (stderr: %s)",
            result.returncode,
            (result.stderr or "").strip()[:300],
        )
        return False
    return True


def webengine_browser_fallback_reason() -> str:
    return (
        "In-app HTML preview is disabled on this Mac because Qt WebEngine crashes "
        "(GPUInfo / segmentation fault).\n\n"
        "HTML files open in your system browser instead.\n\n"
        "To retry embedded preview after upgrading PySide6:\n"
        "  delete the probe cache or set REPEATVIEWER_EMBED_HTML=1\n"
        "To always use the browser:\n"
        "  REPEATVIEWER_EMBED_HTML=0"
    )
