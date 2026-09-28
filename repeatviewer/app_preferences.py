"""General application preferences: startup, performance, diagnostics.

These are small, self-contained QSettings-backed toggles used by the
Preferences dialog and read by the app at startup / by the thread budget.
"""

from __future__ import annotations

import logging
import logging.handlers
import os
import shutil
import sys
import tempfile
import threading
import traceback
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import QSettings

APP_SETTINGS_GROUP = "app"

DEFAULT_REOPEN_LAST_FOLDER = True
DEFAULT_RESTORE_GEOMETRY = True
DEFAULT_SHOW_SPLASH = True
DEFAULT_VERBOSE_LOGGING = False
# 0 means "auto" (let thread_budget pick its default reserve).
DEFAULT_UI_THREAD_RESERVE = 0
MAX_UI_THREAD_RESERVE = 8
_EXCEPTION_HOOKS_INSTALLED = False


def _get_bool(key: str, default: bool) -> bool:
    s = QSettings()
    s.beginGroup(APP_SETTINGS_GROUP)
    value = s.value(key, default)
    s.endGroup()
    if isinstance(value, str):
        return value.lower() in ("1", "true", "yes", "on")
    return bool(value)


def _set_value(key: str, value) -> None:
    s = QSettings()
    s.beginGroup(APP_SETTINGS_GROUP)
    s.setValue(key, value)
    s.endGroup()


def load_reopen_last_folder() -> bool:
    return _get_bool("reopen_last_folder", DEFAULT_REOPEN_LAST_FOLDER)


def save_reopen_last_folder(enabled: bool) -> None:
    _set_value("reopen_last_folder", bool(enabled))


def load_restore_geometry() -> bool:
    return _get_bool("restore_geometry", DEFAULT_RESTORE_GEOMETRY)


def save_restore_geometry(enabled: bool) -> None:
    _set_value("restore_geometry", bool(enabled))


def load_show_splash() -> bool:
    return _get_bool("show_splash", DEFAULT_SHOW_SPLASH)


def save_show_splash(enabled: bool) -> None:
    _set_value("show_splash", bool(enabled))


def load_cyto_experimenter() -> str:
    """Default experimenter name for attached cytogenetic images."""
    s = QSettings()
    s.beginGroup(APP_SETTINGS_GROUP)
    value = s.value("cyto_experimenter", "")
    s.endGroup()
    return str(value or "").strip()


def save_cyto_experimenter(name: str) -> None:
    _set_value("cyto_experimenter", str(name or "").strip())


def load_verbose_logging() -> bool:
    return _get_bool("verbose_logging", DEFAULT_VERBOSE_LOGGING)


def save_verbose_logging(enabled: bool) -> None:
    _set_value("verbose_logging", bool(enabled))


def load_ui_thread_reserve() -> int:
    s = QSettings()
    s.beginGroup(APP_SETTINGS_GROUP)
    try:
        value = int(s.value("ui_thread_reserve", DEFAULT_UI_THREAD_RESERVE))
    except (TypeError, ValueError):
        value = DEFAULT_UI_THREAD_RESERVE
    s.endGroup()
    return max(0, min(MAX_UI_THREAD_RESERVE, value))


def save_ui_thread_reserve(value: int) -> None:
    _set_value("ui_thread_reserve", max(0, min(MAX_UI_THREAD_RESERVE, int(value))))


# --- Filesystem locations -------------------------------------------------

def config_dir() -> Path:
    base = os.environ.get("XDG_CONFIG_HOME", str(Path.home() / ".config"))
    return Path(base) / "RepeatViewer"


def _legacy_teviewer_blast_db_dir() -> Path:
    """Pre-rename app stored named BLAST indexes here (not used as the default)."""
    base = os.environ.get("XDG_CONFIG_HOME", str(Path.home() / ".config"))
    return Path(base) / "TEviewer" / "blast_databases"


def _obsolete_hidden_blast_db_dirs() -> frozenset[Path]:
    """Old defaults under ``.config`` — still readable, but not the default."""
    paths: list[Path] = []
    for raw in (config_dir() / "blast_databases", _legacy_teviewer_blast_db_dir()):
        try:
            paths.append(raw.expanduser().resolve())
        except OSError:
            paths.append(raw.expanduser())
    return frozenset(paths)


def default_blast_db_base_dir() -> Path:
    """Named BLAST indexes: ``Documents/RepeatViewer/blast_databases`` on all OSes.

    Uses the platform Documents folder (macOS / Linux / Windows), not hidden
    ``.config`` or ``AppData``.
    """
    try:
        from PySide6.QtCore import QStandardPaths

        documents = QStandardPaths.writableLocation(
            QStandardPaths.StandardLocation.DocumentsLocation
        )
        if documents:
            return Path(documents) / "RepeatViewer" / "blast_databases"
    except Exception:
        pass
    return Path.home() / "Documents" / "RepeatViewer" / "blast_databases"


def _dir_has_entries(path: Path) -> bool:
    try:
        return path.is_dir() and any(path.iterdir())
    except OSError:
        return False


def load_blast_db_base_dir() -> Path:
    """Return the folder for new named BLAST databases.

    Preference order:
    1. Explicit Preferences path (unless it is an obsolete empty ``.config`` default)
    2. ``Documents/RepeatViewer/blast_databases`` on every platform
    """
    s = QSettings()
    s.beginGroup(APP_SETTINGS_GROUP)
    value = s.value("blast_db_base_dir")
    s.endGroup()
    text = str(value).strip() if value else ""
    documents_default = default_blast_db_base_dir()

    if text:
        path = Path(text).expanduser()
        try:
            resolved = path.resolve()
        except OSError:
            resolved = path
        # Migrate away from the old hidden default when that folder is unused.
        if resolved in _obsolete_hidden_blast_db_dirs() and not _dir_has_entries(path):
            return documents_default
        return path

    return documents_default


def save_blast_db_base_dir(path: Path) -> None:
    _set_value("blast_db_base_dir", str(path.expanduser()))


def diagnostics_dir() -> Path:
    """User-visible location for logs and crash reports."""
    if sys.platform == "darwin":
        return Path.home() / "Library" / "Logs" / "RepeatViewer"
    return config_dir() / "logs"


def cache_dir() -> Path:
    return Path(tempfile.gettempdir()) / "repeatviewer"


def cache_dirs() -> list[Path]:
    """All on-disk cache locations the app writes to."""
    return [cache_dir()]


def log_file_path() -> Path:
    return diagnostics_dir() / "repeatviewer.log"


def crash_reports_dir() -> Path:
    return diagnostics_dir() / "crash-reports"


def write_exception_report(
    *,
    title: str,
    exc_type,
    exc_value,
    exc_traceback,
    context: str = "",
) -> Path | None:
    """Write a human-readable exception report; best effort."""
    try:
        out_dir = crash_reports_dir()
        out_dir.mkdir(parents=True, exist_ok=True)
        stamp = datetime.now().strftime("%Y%m%d-%H%M%S")
        path = out_dir / f"{stamp}-{title}.log"
        with path.open("w", encoding="utf-8") as fh:
            fh.write(f"title: {title}\n")
            fh.write(f"time: {datetime.now().isoformat(timespec='seconds')}\n")
            fh.write(f"python: {sys.version}\n")
            if context:
                fh.write(f"context: {context}\n")
            fh.write("\n")
            traceback.print_exception(exc_type, exc_value, exc_traceback, file=fh)
        return path
    except OSError:
        return None


def _dir_size_bytes(path: Path) -> int:
    total = 0
    if not path.exists():
        return 0
    for root, _dirs, files in os.walk(path):
        for name in files:
            fp = Path(root) / name
            try:
                total += fp.stat().st_size
            except OSError:
                continue
    return total


def cache_size_bytes() -> int:
    return sum(_dir_size_bytes(p) for p in cache_dirs())


def human_size(num_bytes: int) -> str:
    size = float(num_bytes)
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            return f"{size:.0f} {unit}" if unit == "B" else f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


def clear_cache() -> int:
    """Delete the on-disk cache directories. Returns bytes freed (best effort)."""
    freed = cache_size_bytes()
    for path in cache_dirs():
        if path.exists():
            shutil.rmtree(path, ignore_errors=True)
    return freed


def reset_all_preferences() -> None:
    """Clear UI preferences in QSettings, keeping user workspace content.

    Preserves folders/genomes/databases, command wrappers, workflows, and
    folder colors so Restore defaults only resets Preference-dialog settings.
    """
    s = QSettings()
    # QSettings.allKeys() returns "group/key" paths.
    preserve_prefixes = ("project/", "workflows/", "folders/")
    preserve_keys = {"plugins/custom_commands_json"}
    saved: dict[str, object] = {}
    for key in s.allKeys():
        if key in preserve_keys or any(key.startswith(prefix) for prefix in preserve_prefixes):
            saved[key] = s.value(key)
    s.clear()
    for key, value in saved.items():
        s.setValue(key, value)


def apply_log_level(verbose: bool | None = None) -> None:
    """Set the root logging level from the verbose-logging preference."""
    if verbose is None:
        verbose = load_verbose_logging()
    logging.getLogger().setLevel(logging.DEBUG if verbose else logging.INFO)


def setup_file_logging() -> None:
    """Attach a rotating file handler so the log file (and the Diagnostics
    'Open log file' button) actually contain something. Safe to call once."""
    path = log_file_path()
    try:
        path.parent.mkdir(parents=True, exist_ok=True)
    except OSError:
        return
    root = logging.getLogger()
    for handler in root.handlers:
        if getattr(handler, "_repeatviewer_file_handler", False):
            return
    try:
        handler = logging.handlers.RotatingFileHandler(
            path, maxBytes=10_000_000, backupCount=5, encoding="utf-8"
        )
    except OSError:
        return
    handler._repeatviewer_file_handler = True  # type: ignore[attr-defined]
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
    )
    root.addHandler(handler)


def install_exception_logging_hooks() -> None:
    """Log uncaught exceptions and write crash reports automatically."""
    global _EXCEPTION_HOOKS_INSTALLED
    if _EXCEPTION_HOOKS_INSTALLED:
        return
    _EXCEPTION_HOOKS_INSTALLED = True

    prev_excepthook = sys.excepthook

    def _log_uncaught(exc_type, exc_value, exc_traceback) -> None:
        logger = logging.getLogger(__name__)
        logger.exception(
            "Uncaught exception", exc_info=(exc_type, exc_value, exc_traceback)
        )
        report = write_exception_report(
            title="uncaught-exception",
            exc_type=exc_type,
            exc_value=exc_value,
            exc_traceback=exc_traceback,
            context="main-thread",
        )
        if report is not None:
            logger.error("Exception report written to %s", report)
        prev_excepthook(exc_type, exc_value, exc_traceback)

    sys.excepthook = _log_uncaught

    if hasattr(threading, "excepthook"):
        prev_threading_hook = threading.excepthook

        def _thread_excepthook(args) -> None:
            logger = logging.getLogger(__name__)
            thread_name = getattr(args.thread, "name", "unknown")
            logger.exception(
                "Uncaught thread exception (%s)",
                thread_name,
                exc_info=(args.exc_type, args.exc_value, args.exc_traceback),
            )
            report = write_exception_report(
                title="thread-exception",
                exc_type=args.exc_type,
                exc_value=args.exc_value,
                exc_traceback=args.exc_traceback,
                context=f"thread:{thread_name}",
            )
            if report is not None:
                logger.error("Thread exception report written to %s", report)
            prev_threading_hook(args)

        threading.excepthook = _thread_excepthook
