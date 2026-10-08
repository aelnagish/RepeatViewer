"""Back / forward / up controls for the file browser."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import Signal
from PySide6.QtWidgets import QHBoxLayout, QSizePolicy, QToolButton, QWidget


class FolderNavigationBar(QWidget):
    """Compact navigation: file back/forward + up to parent folder.

    ◀ / ▶ step through recently opened files (including in the same folder).
    ↑ goes to the parent of the current folder.
    """

    navigate_requested = Signal(object, bool)  # folder Path, record_history
    file_navigate_requested = Signal(object, bool)  # file Path, record_history

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("FolderNavigationBar")
        self._roots_provider: Callable[[], list[Path]] = list
        self._current_folder: Path | None = None
        self._current_file: Path | None = None

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)

        self._back_btn = QToolButton(self)
        self._back_btn.setText("◀")
        self._back_btn.setToolTip("Back to previous file")
        self._back_btn.setAutoRaise(True)
        self._back_btn.clicked.connect(self._on_back)

        self._forward_btn = QToolButton(self)
        self._forward_btn.setText("▶")
        self._forward_btn.setToolTip("Forward to next file")
        self._forward_btn.setAutoRaise(True)
        self._forward_btn.clicked.connect(self._on_forward)

        self._up_btn = QToolButton(self)
        self._up_btn.setText("↑")
        self._up_btn.setToolTip("Up to parent folder")
        self._up_btn.setAutoRaise(True)
        self._up_btn.clicked.connect(self._on_up)

        for button in (self._back_btn, self._forward_btn, self._up_btn):
            button.setFixedSize(18, 18)
            layout.addWidget(button)

        self.setFixedHeight(18)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

        self._file_history: list[Path] = []
        self._file_history_index = -1
        self._update_buttons()

    def set_roots_provider(self, provider: Callable[[], list[Path]]) -> None:
        self._roots_provider = provider

    def set_folder(self, folder: Path | None, *, record_history: bool = True) -> None:
        """Update the current folder (for ↑). ``record_history`` is unused."""
        del record_history
        if folder is None:
            self._current_folder = None
            self._update_buttons()
            return
        try:
            self._current_folder = folder.expanduser().resolve()
        except OSError:
            self._current_folder = folder.expanduser()
        self._update_buttons()

    def set_file(self, path: Path | None, *, record_history: bool = True) -> None:
        """Update the current file; optionally append to last-used file history."""
        if path is None:
            self._current_file = None
            self._update_buttons()
            return
        try:
            resolved = path.expanduser().resolve()
        except OSError:
            resolved = path.expanduser()
        if record_history:
            if (
                self._file_history_index >= 0
                and self._file_history[self._file_history_index] == resolved
            ):
                pass
            else:
                if self._file_history_index >= 0:
                    self._file_history = self._file_history[
                        : self._file_history_index + 1
                    ]
                self._file_history.append(resolved)
                self._file_history_index = len(self._file_history) - 1
        self._current_file = resolved
        self._update_buttons()

    def _update_buttons(self) -> None:
        self._back_btn.setEnabled(self._file_history_index > 0)
        self._forward_btn.setEnabled(
            self._file_history_index >= 0
            and self._file_history_index < len(self._file_history) - 1
        )
        parent = (
            self._current_folder.parent if self._current_folder is not None else None
        )
        self._up_btn.setEnabled(parent is not None and parent.is_dir())

    def _on_back(self) -> None:
        if self._file_history_index <= 0:
            return
        self._file_history_index -= 1
        self._go_to_file_history_entry(self._file_history[self._file_history_index])

    def _on_forward(self) -> None:
        if (
            self._file_history_index < 0
            or self._file_history_index >= len(self._file_history) - 1
        ):
            return
        self._file_history_index += 1
        self._go_to_file_history_entry(self._file_history[self._file_history_index])

    def _on_up(self) -> None:
        if self._current_folder is None:
            return
        parent = self._current_folder.parent
        if parent.is_dir():
            self.navigate_requested.emit(parent, True)

    def _go_to_file_history_entry(self, path: Path) -> None:
        try:
            resolved = path.expanduser().resolve()
        except OSError:
            resolved = path.expanduser()
        self._current_file = resolved
        self._update_buttons()
        self.file_navigate_requested.emit(resolved, False)
