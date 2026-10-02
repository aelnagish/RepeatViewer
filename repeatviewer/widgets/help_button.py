"""Shared round "?" help button (rule 8: one help style everywhere).

The look follows the MSA clean dialog ("Clean rows with divergence"): a 22 x 22
round button next to the item. Hover shows the help as a tooltip; click opens
the same text in a message box.
"""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QToolButton,
    QWidget,
)

HELP_BUTTON_SIZE = 22
TOOLTIP_MAX_WIDTH = 420

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
    return f"<div style='max-width:{TOOLTIP_MAX_WIDTH}px'>{help_text}</div>"


def show_help_message(parent: QWidget | None, title: str, help_text: str) -> None:
    """Show help; body stays normal weight (only explicit <b>…</b> is bold)."""
    box = QMessageBox(parent)
    box.setWindowTitle(title)
    box.setIcon(QMessageBox.Icon.Information)
    box.setTextFormat(Qt.TextFormat.RichText)
    # Keep body at normal weight so only explicit <b>…</b> labels stand out.
    # Native macOS alerts often force the whole setText() body bold.
    box.setText(
        f"<div style='font-weight:400; max-width:{TOOLTIP_MAX_WIDTH}px'>"
        f"{help_text}</div>"
    )
    box.setStandardButtons(QMessageBox.StandardButton.Ok)
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


def label_with_help(
    text: str,
    title: str,
    help_text: str,
    *,
    parent: QWidget | None = None,
) -> QWidget:
    """Right-aligned ``Foo:`` + ``?`` for a form label column.

    Use for path rows so Browse lines up with Define-inputs Browse (help is not
    after the Browse button).
    """
    wrap = QWidget(parent)
    row = QHBoxLayout(wrap)
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(4)
    lab = QLabel(text, wrap)
    lab.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    row.addWidget(lab, stretch=1)
    row.addWidget(make_help_button(wrap, title, help_text))
    return wrap