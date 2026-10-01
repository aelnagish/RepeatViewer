"""Shared round "?" help button (rule 8: one help style everywhere).

The look follows the MSA clean dialog ("Clean rows with divergence"): a 22 x 22
round button next to the item. Hover shows the help as a tooltip; click opens
the same text in a message box.
"""

from __future__ import annotations

import re

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QHBoxLayout,
    QMessageBox,
    QToolButton,
    QWidget,
)

HELP_BUTTON_SIZE = 22
TOOLTIP_MAX_WIDTH = 420
_HTML_TAG = re.compile(r"<[A-Za-z/][^>]*>")

HELP_BUTTON_STYLE = (
    "QToolButton {"
    "  border: 1px solid palette(mid);"
    "  border-radius: 10px;"
    "  background: palette(base);"
    "  font-weight: 600;"
    "  padding: 0px;"
    "}"
    "QToolButton:hover { background: palette(alternate-base); }"
)


def help_tooltip_text(help_text: str) -> str:
    """Rich text is width-capped so long help does not span the screen."""
    if _HTML_TAG.search(help_text):
        return f"<div style='max-width:{TOOLTIP_MAX_WIDTH}px'>{help_text}</div>"
    return help_text


def show_help_message(parent: QWidget | None, title: str, help_text: str) -> None:
    """Show help with HTML bold labels; native macOS alerts ignore that markup."""
    box = QMessageBox(parent)
    box.setWindowTitle(title)
    box.setIcon(QMessageBox.Icon.Information)
    box.setTextFormat(Qt.TextFormat.RichText)
    # Keep body at normal weight so only explicit <b>…</b> labels stand out.
    body = (
        f"<div style='font-weight:400'>{help_text}</div>"
        if _HTML_TAG.search(help_text)
        else help_text
    )
    box.setText(body)
    box.setStandardButtons(QMessageBox.StandardButton.Ok)
    # Native alerts strip/ignore HTML and often render the whole body bold.
    box.setOption(QMessageBox.Option.DontUseNativeDialog, True)
    box.exec()


def make_help_button(parent: QWidget | None, title: str, help_text: str) -> QToolButton:
    """Round "?" button: hover shows ``help_text``, click opens it titled ``title``."""
    button = QToolButton(parent)
    button.setText("?")
    button.setToolTip(help_tooltip_text(help_text))
    button.setAutoRaise(False)
    button.setFixedSize(HELP_BUTTON_SIZE, HELP_BUTTON_SIZE)
    button.setStyleSheet(HELP_BUTTON_STYLE)
    button.setAccessibleName(f"Help: {title}")
    button.clicked.connect(
        lambda checked=False: show_help_message(button.window(), title, help_text)
    )
    return button


def field_with_help(
    field: QWidget,
    title: str,
    help_text: str,
    *,
    parent: QWidget | None = None,
    stretch: bool = False,
) -> QWidget:
    """``field`` with the "?" button to its right, for a form row."""
    wrap = QWidget(parent)
    row = QHBoxLayout(wrap)
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(6)
    row.addWidget(field)
    row.addWidget(make_help_button(wrap, title, help_text))
    if stretch:
        row.addStretch(1)
    return wrap
