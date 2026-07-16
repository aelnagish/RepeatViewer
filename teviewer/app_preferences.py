"""General application preferences: startup, performance, diagnostics.

These are small, self-contained QSettings-backed toggles used by the
Preferences dialog and read by the app at startup / by the thread budget.
"""

from __future__ import annotations

import logging
import logging.handlers
import os
import shutil
import tempfile
from pathlib import Path

from PySide6.QtCore import QSettings

APP_SETTINGS_GROUP = "app"

DEFAULT_REOPEN_LAST_FOLDER = True
DEFAULT_RESTORE_GEOMETRY = True
DEFAULT_VERBOSE_LOGGING = False
# 0 means "auto" (let thread_budget pick its default reserve).
DEFAULT_UI_THREAD_RESERVE = 0
MAX_UI_THREAD_RESERVE = 8


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
    return Path(base) / "TEviewer"


def cache_dir() -> Path:
    return Path(tempfile.gettempdir()) / "teviewer"


def cache_dirs() -> list[Path]:
    """All on-disk cache locations the app writes to."""
    tmp = Path(tempfile.gettempdir())
    return [tmp / "teviewer", tmp / "teviewer-igv-sessions"]


def log_file_path() -> Path:
    return config_dir() / "repeatviewer.log"


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
    """Clear every persisted preference (QSettings). Files/projects are untouched."""
    QSettings().clear()


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
        if getattr(handler, "_teviewer_file_handler", False):
            return
    try:
        handler = logging.handlers.RotatingFileHandler(
            path, maxBytes=1_000_000, backupCount=3, encoding="utf-8"
        )
    except OSError:
        return
    handler._teviewer_file_handler = True  # type: ignore[attr-defined]
    handler.setFormatter(
        logging.Formatter("%(asctime)s %(levelname)s %(name)s: %(message)s")
    )
    root.addHandler(handler)
