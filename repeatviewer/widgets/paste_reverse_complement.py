"""Edit > Paste Reverse Complement: paste the reverse complement of the clipboard.

Works in any editable text (the Preview editor in Edit text mode, the Find box,
name fields). Shared by the Edit menu and the Preview editor's right-click menu.
"""

from __future__ import annotations

from PySide6.QtWidgets import (
    QApplication,
    QLineEdit,
    QMenu,
    QPlainTextEdit,
    QTextEdit,
    QWidget,
)

from repeatviewer.edit_actions import reverse_complement_text

MENU_TEXT = "Paste Reverse Complement"
NOT_A_SEQUENCE = "The clipboard does not hold a DNA sequence"


def editable_text_widget(widget: QWidget | None) -> QWidget | None:
    """``widget`` if it is a text box the user can type in, else None."""
    if isinstance(widget, QLineEdit) and not widget.isReadOnly() and widget.isEnabled():
        return widget
    if isinstance(widget, (QTextEdit, QPlainTextEdit)) and not widget.isReadOnly():
        return widget
    return None


def clipboard_has_text() -> bool:
    clipboard = QApplication.clipboard()
    return bool(clipboard is not None and clipboard.text().strip())


def paste_reverse_complement(widget: QWidget | None) -> tuple[bool, str]:
    """Insert the reverse complement of the clipboard at the cursor of ``widget``.

    Replaces the selection, like Paste, and is one Undo step. Returns
    (pasted, message for the status bar).
    """
    target = editable_text_widget(widget)
    if target is None:
        return False, "Click in an editable text first (for example Edit text in Preview)"
    clipboard = QApplication.clipboard()
    result = reverse_complement_text(clipboard.text() if clipboard is not None else "")
    if result is None:
        return False, NOT_A_SEQUENCE
    text, bases = result
    if isinstance(target, QLineEdit):
        # One line: a FASTA header has no place here, only the bases.
        lines = [line for line in text.splitlines() if not line.startswith(">")]
        target.insert("".join(lines))
    else:
        target.insertPlainText(text)
        target.ensureCursorVisible()
    return True, f"Pasted the reverse complement ({bases:,} bp)"


def add_to_context_menu(menu: QMenu, widget: QWidget, on_message=None) -> None:
    """Put Paste Reverse Complement right after Paste in a text box's right-click menu."""
    action = menu.addAction(MENU_TEXT)
    action.setEnabled(editable_text_widget(widget) is not None and clipboard_has_text())

    def run() -> None:
        _ok, message = paste_reverse_complement(widget)
        if on_message is not None:
            on_message(message)

    action.triggered.connect(run)
    paste = next((a for a in menu.actions() if a.objectName() == "edit-paste"), None)
    if paste is None:
        return
    actions = menu.actions()
    index = actions.index(paste)
    if index + 1 < len(actions) and actions[index + 1] is not action:
        menu.removeAction(action)
        menu.insertAction(actions[index + 1], action)
