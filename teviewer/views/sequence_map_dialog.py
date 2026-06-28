"""Dialog hosting the interactive sequence map viewer."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QKeySequence, QShortcut
from PySide6.QtWidgets import QDialog, QHBoxLayout, QMessageBox, QPushButton, QVBoxLayout, QWidget

from teviewer.widgets.sequence_map_viewer import SequenceMapPanel


class SequenceMapDialog(QDialog):
    """Standalone sequence map window (legacy; prefer preview Annotation tab)."""

    def __init__(self, fasta_path: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._path = fasta_path.expanduser().resolve()
        self.setWindowTitle(f"Sequence map — {self._path.name}")
        self.setMinimumSize(920, 520)
        self.resize(1100, 640)

        self._panel = SequenceMapPanel(self._path, self, autoload=True)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_row = QHBoxLayout()
        btn_row.addStretch()
        btn_row.addWidget(close_btn)

        layout = QVBoxLayout(self)
        layout.addWidget(self._panel, stretch=1)
        layout.addLayout(btn_row)

        QShortcut(QKeySequence.StandardKey.Save, self, lambda: self._panel.save_annotations(silent=False))
        for key in (QKeySequence.StandardKey.Delete, QKeySequence(Qt.Key.Key_Backspace)):
            shortcut = QShortcut(key, self)
            shortcut.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            shortcut.activated.connect(
                lambda: self._panel._delete_feature(confirm_empty=False)
            )

    def fasta_path(self) -> Path:
        return self._path

    def closeEvent(self, event) -> None:  # noqa: N802
        if self._panel.is_dirty():
            answer = QMessageBox.question(
                self,
                "Unsaved annotations",
                "Save feature annotations before closing?",
                QMessageBox.StandardButton.Save
                | QMessageBox.StandardButton.Discard
                | QMessageBox.StandardButton.Cancel,
            )
            if answer == QMessageBox.StandardButton.Save:
                if not self._panel.save_annotations(silent=False):
                    event.ignore()
                    return
            elif answer == QMessageBox.StandardButton.Cancel:
                event.ignore()
                return
        super().closeEvent(event)


def open_sequence_map_dialog(fasta_path: Path, parent: QWidget | None = None) -> None:
    """Open or raise the sequence map for a single-sequence FASTA."""
    try:
        dialog = SequenceMapDialog(fasta_path, parent=parent)
    except ValueError as exc:
        QMessageBox.warning(parent, "Sequence map", str(exc))
        return
    dialog.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose)
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()
