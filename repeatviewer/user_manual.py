"""Where the bundled user manual lives, and opening it from the Help menu.

The manual is written in docs/manual/*.md and built into
assets/manuals/user_manual/ by docs/manual/tools/build_manual.py.
"""

from __future__ import annotations

from pathlib import Path

MANUAL_DIR = Path(__file__).resolve().parent / "assets" / "manuals" / "user_manual"
MANUAL_INDEX = MANUAL_DIR / "index.html"
MANUAL_PDF = MANUAL_DIR / "RepeatViewer_User_Manual.pdf"
ONLINE_MANUAL = "https://github.com/qjiangzhao/RepeatViewer/tree/main/docs/manual"


def manual_page(chapter: str = "") -> Path | None:
    """The bundled page of ``chapter`` (for example "12-primers"), or the index.

    Returns ``None`` when the manual was not built into this copy of the app.
    """
    page = MANUAL_DIR / f"{chapter}.html" if chapter else MANUAL_INDEX
    if not page.is_file():
        page = MANUAL_INDEX
    return page if page.is_file() else None


def manual_url(chapter: str = "", anchor: str = "") -> str:
    """A URL for the manual: the bundled page, or the online copy."""
    page = manual_page(chapter)
    if page is None:
        return ONLINE_MANUAL
    url = page.as_uri()
    return f"{url}#{anchor}" if anchor else url


def open_manual(chapter: str = "", anchor: str = "", *, pdf: bool = False) -> bool:
    """Open the manual in the system browser (or the PDF viewer)."""
    from PySide6.QtCore import QUrl
    from PySide6.QtGui import QDesktopServices

    if pdf and MANUAL_PDF.is_file():
        return QDesktopServices.openUrl(QUrl.fromLocalFile(str(MANUAL_PDF)))
    return QDesktopServices.openUrl(QUrl(manual_url(chapter, anchor)))
