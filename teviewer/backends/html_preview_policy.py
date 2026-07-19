"""When Plotly HTML may load inside Qt WebEngine (macOS crash avoidance)."""

from __future__ import annotations

import os
from pathlib import Path

from teviewer.backends.curation import is_flexidot_html_file


def use_embedded_html_preview() -> bool:
    """Return False when ``TEVIEWER_EMBED_HTML=0`` (browser-only for all HTML)."""
    flag = os.environ.get("TEVIEWER_EMBED_HTML", "").strip().lower()
    return flag not in {"0", "false", "no", "off"}


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
