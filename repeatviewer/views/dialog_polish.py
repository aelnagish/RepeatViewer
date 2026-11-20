"""Automatic, app-wide polish for the tool parameter dialogs.

When a dialog from ``repeatviewer.views`` is first shown, give every
combo/spin box one uniform width and every QFormLayout label a uniform
right-aligned column, so the parameter dialogs stay symmetric without
per-dialog code.

IMPORTANT: never install a Python event filter on QApplication. On macOS,
PySide6 can segfault when an application-level filter runs during
``QObject::setProperty`` (nested ``getWrapperForQObject``). Polish uses a
per-dialog object filter instead.

Layout jumps after open: polish and ``prepare_for_show`` must run *before*
``show()`` (see ``prepare_params_dialog_for_show``). The Show filter only
covers dialogs that skipped that helper; it still defers with a zero-delay
timer (never mutate widgets during event delivery).
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


def polish_form_dialog(dialog: QDialog, *, adjust_size: bool = True) -> None:
    """Align combo/spin widths and form labels; optionally ``adjustSize``."""

    def _explicitly_fixed(widget) -> bool:
        # Honor widths already pinned by the dialog (e.g. Blast spin boxes at 150).
        return widget.minimumWidth() > 0 and widget.minimumWidth() == widget.maximumWidth()

    fields = [
        widget
        for widget in (
            list(dialog.findChildren(QComboBox))
            + list(dialog.findChildren(QAbstractSpinBox))
        )
        if not _explicitly_fixed(widget)
    ]
    if fields:
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

    if adjust_size:
        dialog.adjustSize()
    setattr(dialog, "_repeatviewer_form_polished", True)


def prepare_params_dialog_for_show(dialog: QDialog) -> None:
    """Polish + dialog-specific layout *before* the first paint.

    Call this immediately before ``dialog.show()`` so the first frame already
    has aligned columns and the final size (no post-show jump).
    """
    if not isinstance(dialog, QDialog):
        return
    try:
        prepare = getattr(dialog, "prepare_for_show", None)
        has_prepare = callable(prepare)
        # tool_dialogs + any dialog that implements prepare_for_show (Tree,
        # Map, Samtools, Primers, …).
        if _should_polish_dialog(dialog) or has_prepare:
            if not getattr(dialog, "_repeatviewer_form_polished", False):
                # Skip adjustSize when prepare_for_show will set the real height —
                # adjustSize uses the tallest tab and leaves a blank band.
                polish_form_dialog(dialog, adjust_size=not has_prepare)
        if has_prepare:
            prepare()
        _seen.add(dialog)
    except Exception:
        pass


class _DialogShowWatcher(QObject):
    """Object-level filter: only receives events targeted at its dialog."""

    def eventFilter(self, obj: QObject, event: QEvent) -> bool:  # noqa: N802
        if obj is self.parent() and event.type() == QEvent.Type.Show:
            dialog = obj
            if isinstance(dialog, QDialog) and dialog not in _seen:
                _seen.add(dialog)
                # Fallback when open paths forget prepare_params_dialog_for_show.
                QTimer.singleShot(0, lambda d=dialog: _deferred_polish(d))
        return False


def _deferred_polish(dialog: QDialog) -> None:
    try:
        if getattr(dialog, "_repeatviewer_form_polished", False):
            # Already laid out before show — do not polish or prepare again
            # (a second adjustSize / fit is the post-open jump).
            return
        prepare = getattr(dialog, "prepare_for_show", None)
        has_prepare = callable(prepare)
        polish_form_dialog(dialog, adjust_size=not has_prepare)
        if has_prepare:
            prepare()
    except Exception:
        pass


_POLISH_MODULE_PREFIXES = (
    "repeatviewer.views.tool_dialogs",
    "repeatviewer.views.tree_dialog",
    "repeatviewer.views.map_reference_dialog",
    "repeatviewer.views.samtools_dialog",
    "repeatviewer.views.primer_database_dialog",
    "repeatviewer.views.ncbi_download_dialog",
)


def _should_polish_dialog(dialog: QDialog) -> bool:
    module = type(dialog).__module__ or ""
    return module.startswith(_POLISH_MODULE_PREFIXES)


def _patched_qdialog_init(self, *args, **kwargs) -> None:
    _original_qdialog_init(self, *args, **kwargs)
    # Watch tool params dialogs and any dialog that implements prepare_for_show.
    if not (_should_polish_dialog(self) or hasattr(self, "prepare_for_show")):
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
