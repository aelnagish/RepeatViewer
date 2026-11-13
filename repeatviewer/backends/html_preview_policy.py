"""When Plotly HTML may load inside Qt WebEngine (macOS crash avoidance)."""

from __future__ import annotations

import os
from pathlib import Path

from repeatviewer.backends.html_kinds import is_flexidot_html_file


def may_attempt_embedded_html_preview(*, force: bool = False) -> bool:
    """True when RepeatViewer should try loading HTML inside Qt WebEngine.

    Embedded preview is **on by default**. Set ``REPEATVIEWER_EMBED_HTML=0`` to
    always use the system browser. A failed background safety probe does not
    block preview — the main GUI process may still render Plotly HTML fine.
    """
    flag = os.environ.get("REPEATVIEWER_EMBED_HTML", "").strip().lower()
    if flag in {"0", "false", "no", "off"}:
        return False
    if force or flag in {"1", "true", "yes", "on"}:
        return True
    return True


def use_embedded_html_preview() -> bool:
    """Return False only when embedded HTML preview is turned off."""
    return may_attempt_embedded_html_preview()


def embedded_html_disabled_reason() -> str:
    """Human-readable reason embedded HTML preview is unavailable right now."""
    flag = os.environ.get("REPEATVIEWER_EMBED_HTML", "").strip().lower()
    if flag in {"0", "false", "no", "off"}:
        return "Embedded HTML preview disabled (REPEATVIEWER_EMBED_HTML=0)."
    return (
        "Embedded HTML preview could not be started.\n\n"
        "Install PySide6-Addons (pip install 'PySide6-Addons>=6.6') "
        "or use Open in browser."
    )


def use_flexidot_plotly_webengine() -> bool:
    """
    Flexidot Plotly inside Qt WebEngine enables Go-to → Annotation in RepeatViewer.

    Default: embedded preview (same as other HTML). Opt out on macOS if WebEngine
    crashes: ``REPEATVIEWER_FLEXIDOT_WEBENGINE=0`` opens Flexidot in the browser instead.
    """
    flag = os.environ.get("REPEATVIEWER_FLEXIDOT_WEBENGINE", "").strip().lower()
    if flag in {"0", "false", "no", "off"}:
        return False
    if flag in {"1", "true", "yes", "on"}:
        return use_embedded_html_preview()
    return use_embedded_html_preview()


def should_use_webengine_for_html(path: Path) -> bool:
    """True when this HTML file should load in the in-app WebEngine viewer."""
    if not may_attempt_embedded_html_preview():
        return False
    resolved = path.expanduser().resolve()
    if is_flexidot_html_file(resolved):
        return use_flexidot_plotly_webengine()
    return True


def flexidot_browser_fallback_reason(path: Path) -> str:
    resolved = path.expanduser().resolve()
    return (
        "Flexidot HTML is opened in your browser (embedded preview disabled).\n\n"
        f"File: {resolved.name}\n\n"
        "Go-to → Annotation only works when the dotplot is shown inside RepeatViewer.\n"
        "Remove REPEATVIEWER_FLEXIDOT_WEBENGINE=0 and restart to use embedded preview.\n\n"
        "Use the button below to reopen the plot in your browser."
    )
