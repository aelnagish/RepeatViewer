"""Automatic, app-wide polish for the tool parameter dialogs.

When a dialog from ``repeatviewer.views.tool_dialogs`` is first shown, give every
combo/spin box one uniform width and every QFormLayout label a uniform
right-aligned column, so the parameter dialogs stay symmetric without
per-dialog code.

IMPORTANT: never install a Python event filter on QApplication. On macOS,
PySide6 can segfault when an application-level filter runs during
``QObject::setProperty`` (nested ``getWrapperForQObject``). Polish uses a
per-dialog object filter instead, and only *schedules* work with a
zero-delay timer (never mutate widgets during event delivery).
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

_original_qdialog_init = None
_seen: WeakSet = WeakSet()
_installed = False


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


class _DialogShowWatcher(QObject):
    """Object-level filter: only receives events targeted at its dialog."""

    def eventFilter(self, obj: QObject, event: QEvent) -> bool:  # noqa: N802
        if obj is self.parent() and event.type() == QEvent.Type.Show:
            dialog = obj
            if isinstance(dialog, QDialog) and dialog not in _seen:
                _seen.add(dialog)
                QTimer.singleShot(0, lambda d=dialog: _deferred_polish(d))
        return False


def _deferred_polish(dialog: QDialog) -> None:
    try:
        polish_form_dialog(dialog)
    except Exception:
        pass


def _should_polish_dialog(dialog: QDialog) -> bool:
    module = type(dialog).__module__ or ""
    return module.startswith("repeatviewer.views.tool_dialogs")


def _patched_qdialog_init(self, *args, **kwargs) -> None:
    _original_qdialog_init(self, *args, **kwargs)
    if not _should_polish_dialog(self):
        return
    if getattr(self, "_repeatviewer_dialog_show_watcher", None) is not None:
        return
    watcher = _DialogShowWatcher(self)
    self.installEventFilter(watcher)
    self._repeatviewer_dialog_show_watcher = watcher


def install_dialog_polisher(app) -> None:
    """Install dialog polish once (per-dialog filters; never on QApplication)."""
    global _installed, _original_qdialog_init
    del app  # kept for call-site compatibility
    if _installed:
        return
    _original_qdialog_init = QDialog.__init__
    QDialog.__init__ = _patched_qdialog_init  # type: ignore[method-assign]
    _installed = True
