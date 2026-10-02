"""Rule 8: every "?" help button uses the one shared MSA-clean style."""

from __future__ import annotations

import re
from pathlib import Path

PACKAGE = Path(__file__).resolve().parents[1] / "repeatviewer"
SHARED = PACKAGE / "widgets" / "help_button.py"


def test_no_hand_built_help_buttons():
    """Only widgets/help_button.py may build a "?" button or its stylesheet."""
    pattern = re.compile(r'setText\(\s*"\?"\s*\)|border-radius:\s*10px')
    offenders = [
        f"{path.relative_to(PACKAGE.parent)}:{n}"
        for path in PACKAGE.rglob("*.py")
        if path != SHARED
        for n, line in enumerate(path.read_text(encoding="utf-8").splitlines(), 1)
        if pattern.search(line)
    ]
    assert not offenders, "use widgets.help_button.make_help_button: " + ", ".join(offenders)


def test_help_button_look_and_tooltip(qapp):
    from repeatviewer.widgets.help_button import (
        HELP_BUTTON_SIZE,
        HELP_BUTTON_STYLE,
        make_help_button,
    )

    button = make_help_button(None, "Window size", "Sliding window in bp.")
    assert button.text() == "?"
    assert button.width() == button.height() == HELP_BUTTON_SIZE
    assert button.styleSheet() == HELP_BUTTON_STYLE
    assert button.toolTip().startswith("<div style='max-width:")
    assert "Sliding window in bp." in button.toolTip()


def test_rich_help_tooltip_is_width_capped(qapp):
    from repeatviewer.widgets.help_button import make_help_button

    button = make_help_button(None, "Model", "<b>HKY</b> — unequal base frequencies.")
    assert button.toolTip().startswith("<div style='max-width:")
    assert "<b>HKY</b>" in button.toolTip()


def test_plain_help_message_is_not_bold(qapp, monkeypatch):
    from repeatviewer.widgets import help_button

    seen: dict[str, object] = {}

    class _FakeBox:
        Option = help_button.QMessageBox.Option
        StandardButton = help_button.QMessageBox.StandardButton
        Icon = help_button.QMessageBox.Icon

        def __init__(self, parent=None):
            seen["parent"] = parent

        def setWindowTitle(self, title):
            seen["title"] = title

        def setIcon(self, icon):
            seen["icon"] = icon

        def setTextFormat(self, fmt):
            seen["format"] = fmt

        def setText(self, text):
            seen["text"] = text

        def setStandardButtons(self, buttons):
            seen["buttons"] = buttons

        def setOption(self, option, on=True):
            seen["option"] = (option, on)

        def exec(self):
            seen["exec"] = True

    monkeypatch.setattr(help_button, "QMessageBox", _FakeBox)
    help_button.show_help_message(None, "Show GFF features", "Short plain help.")
    assert "font-weight:400" in str(seen["text"])
    assert "Short plain help." in str(seen["text"])
    assert seen["option"] == (_FakeBox.Option.DontUseNativeDialog, True)


def test_field_with_help_puts_button_right_of_field(qapp):
    from PySide6.QtWidgets import QSpinBox, QToolButton

    from repeatviewer.widgets.help_button import field_with_help

    spin = QSpinBox()
    row = field_with_help(spin, "Window size", "Sliding window in bp.")
    layout = row.layout()
    assert layout.itemAt(0).widget() is spin
    assert isinstance(layout.itemAt(1).widget(), QToolButton)

