"""When Plotly HTML may load inside Qt WebEngine (macOS crash avoidance)."""

from __future__ import annotations

import os
from pathlib import Path

from teviewer.backends.curation import is_flexidot_html_file


def use_embedded_html_preview() -> bool:
    """Return False when embedded HTML preview should not be used.

    Honors ``TEVIEWER_EMBED_HTML`` overrides, and otherwise defers to the
    Qt WebEngine crash probe (macOS + certain PySide6 builds segfault on
    ``QWebEngineView`` creation) so we never build a WebEngine view in the
    main process unless it is known to be safe on this machine.
    """
    from teviewer.backends.webengine_probe import webengine_is_safe

    return webengine_is_safe()


def embedded_html_disabled_reason() -> str:
    """Human-readable reason embedded HTML preview is unavailable right now."""
    flag = os.environ.get("TEVIEWER_EMBED_HTML", "").strip().lower()
    if flag in {"0", "false", "no", "off"}:
        return "Embedded HTML preview disabled (TEVIEWER_EMBED_HTML=0)."
    from teviewer.backends.webengine_probe import webengine_browser_fallback_reason

    return webengine_browser_fallback_reason()


def use_flexidot_plotly_webengine() -> bool:
    """
    Flexidot Plotly inside Qt WebEngine enables Go-to → Annotation in TEviewer.

    Default: embedded preview (same as other HTML). Opt out on macOS if WebEngine
    crashes: ``TEVIEWER_FLEXIDOT_WEBENGINE=0`` opens Flexidot in the browser instead.
    """
    flag = os.environ.get("TEVIEWER_FLEXIDOT_WEBENGINE", "").strip().lower()
    if flag in {"0", "false", "no", "off"}:
        return False
    if flag in {"1", "true", "yes", "on"}:
        return use_embedded_html_preview()
    return use_embedded_html_preview()


def should_use_webengine_for_html(path: Path) -> bool:
    """True when this HTML file should load in the in-app WebEngine viewer."""
    if not use_embedded_html_preview():
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
        "Go-to → Annotation only works when the dotplot is shown inside TEviewer.\n"
        "Remove TEVIEWER_FLEXIDOT_WEBENGINE=0 and restart to use embedded preview.\n\n"
        "Use the button below to reopen the plot in your browser."
    )
