"""UI widgets package.

Keep this module import-light: ``from repeatviewer.widgets.app_icon`` (and
similar) must not pull HTML preview / matplotlib / WebEngine at process start.
Import concrete modules directly, e.g. ``repeatviewer.widgets.preview_tabs``.
"""

from __future__ import annotations

__all__: list[str] = []
