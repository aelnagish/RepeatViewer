"""Files / Folders panel clicks fill the path fields of open Cytogenetics windows, as in every tool window.

A window registers itself and lists its path fields (``PathBrowseRow``). A click in the Files or Folders panel
goes to the field the user clicked last in the newest visible window; the window decides if the path fits.
"""

from __future__ import annotations

import weakref
from pathlib import Path

_windows: list[weakref.ref] = []


def register(window) -> None:
    """Call once per window; the window needs ``accept_panel_path(path) -> bool``."""
    _windows[:] = [w for w in _windows if w() is not None and w() is not window]
    _windows.append(weakref.ref(window))


def _visible() -> list:
    out = []
    for ref in reversed(_windows):
        w = ref()
        try:
            if w is not None and w.isVisible():
                out.append(w)
        except RuntimeError:  # the C++ object is gone
            continue
    return out


def any_visible() -> bool:
    return bool(_visible())


def offer(path: Path) -> bool:
    """Give a clicked file or folder to the newest open window that takes it."""
    for w in _visible():
        try:
            if w.accept_panel_path(Path(path)):
                return True
        except RuntimeError:
            continue
    return False


class PanelPathFields:
    """Mixin for a window with path fields: remembers the field clicked last and fills it from the panels."""

    def watch_path_fields(self, *rows) -> None:
        self._panel_rows = [r for r in rows if r is not None]
        self._panel_row = self._panel_rows[0] if self._panel_rows else None
        for row in self._panel_rows:
            row.activated.connect(lambda r=row: setattr(self, "_panel_row", r))
        register(self)

    def accept_panel_path(self, path: Path) -> bool:
        row = getattr(self, "_panel_row", None)
        if row is None or not row.isVisible() or not row.isEnabled():
            return False
        mode = getattr(row, "_mode", "file")
        if mode == "directory" and path.is_file():
            path = path.parent
        if mode == "file" and path.is_dir():
            return False
        suffixes = getattr(row, "panel_suffixes", None)
        if suffixes and path.is_file() and not path.name.lower().endswith(tuple(suffixes)):
            return False
        row.set_text(str(path))
        return True
