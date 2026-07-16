"""Automatic, app-wide polish for the tool parameter dialogs.

An application event filter catches the first show of any dialog defined in
``teviewer.views.tool_dialogs`` and gives every combo/spin box one uniform
width and every QFormLayout label a uniform right-aligned column, so the
parameter dialogs stay symmetric without per-dialog code.

IMPORTANT: the event filter must not mutate widgets (or call setProperty)
during event delivery — that re-enters the application event filters and can
crash PySide. So we only *schedule* the polish with a zero-delay timer and
track already-seen dialogs in a WeakSet (never setProperty).
"""

from __future__ import annotations

from weakref import WeakSet

from PySide6.QtCore import QEvent, QObject, Qt, QTimer
from PySide6.QtWidgets import (
    QAbstractSpinBox,
    QComboBox,
    QDialog,
    QFormLayout,
    QLabel,
)

_MIN_FIELD_WIDTH = 130
_MAX_FIELD_WIDTH = 300


def polish_form_dialog(dialog: QDialog) -> None:
    fields = list(dialog.findChildren(QComboBox)) + list(dialog.findChildren(QAbstractSpinBox))
    if not fields:
        return

    field_width = max([_MIN_FIELD_WIDTH] + [f.sizeHint().width() for f in fields])
    field_width = min(field_width, _MAX_FIELD_WIDTH)
    for field in fields:
        field.setFixedWidth(field_width)

    labels: list[QLabel] = []
    for form in dialog.findChildren(QFormLayout):
        for row in range(form.rowCount()):
            item = form.itemAt(row, QFormLayout.ItemRole.LabelRole)
            widget = item.widget() if item is not None else None
            if isinstance(widget, QLabel):
                labels.append(widget)
    if labels:
        label_width = max(lbl.sizeHint().width() for lbl in labels)
        for lbl in labels:
            lbl.setFixedWidth(label_width)
            lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)

    dialog.adjustSize()


class _DialogPolisher(QObject):
    def __init__(self, parent: QObject | None = None) -> None:
        super().__init__(parent)
        self._seen: WeakSet = WeakSet()

    def eventFilter(self, obj: QObject, event: QEvent) -> bool:  # noqa: N802
        # Keep this trivial: never mutate widgets or call setProperty here.
        if event.type() == QEvent.Type.Show and isinstance(obj, QDialog):
            module = type(obj).__module__ or ""
            if module.startswith("teviewer.views.tool_dialogs") and obj not in self._seen:
                self._seen.add(obj)
                QTimer.singleShot(0, lambda o=obj: _deferred_polish(o))
        return False


def _deferred_polish(dialog: QDialog) -> None:
    try:
        polish_form_dialog(dialog)
    except Exception:
        pass


_polisher: _DialogPolisher | None = None


def install_dialog_polisher(app) -> None:
    """Install the polisher once on the QApplication."""
    global _polisher
    if _polisher is None:
        _polisher = _DialogPolisher(app)
        app.installEventFilter(_polisher)
