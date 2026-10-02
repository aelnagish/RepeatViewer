"""Shared file/folder path field: PathField line edit + Browse (MSA style)."""

from __future__ import annotations

from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import QEvent, Qt, Signal
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QLineEdit,
    QMenu,
    QPushButton,
    QSizePolicy,
    QWidget,
)

PATH_FIELD_MIN_WIDTH = 280


class PathBrowseRow(QWidget):
    """Line edit + Browse matching the MSA FASTA path-filling style.

    Form labels belong on the caller (``form.addRow("MSA FASTA:", row)``).
    """

    path_changed = Signal()
    activated = Signal(object)
    text_changed = Signal(str)

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        initial: Path | str | None = None,
        placeholder: str = "",
        browse_title: str = "Select file",
        file_filter: str = "All files (*)",
        start_dir: str = "",
        mode: str = "file",
        clear_button: bool = True,
        minimum_width: int = PATH_FIELD_MIN_WIDTH,
        line_edit: QLineEdit | None = None,
        browse_handler: Callable[[], None] | None = None,
    ) -> None:
        super().__init__(parent)
        if mode not in ("file", "directory", "file_or_directory"):
            raise ValueError(f"Unsupported PathBrowseRow mode: {mode!r}")
        self._mode = mode
        self._browse_title = browse_title
        self._file_filter = file_filter
        self._start_dir = start_dir
        self._default_tooltip = ""
        self._browse_handler = browse_handler

        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)

        self._path = line_edit if line_edit is not None else QLineEdit(self)
        if self._path.parent() is None:
            self._path.setParent(self)
        self._path.setObjectName("PathField")
        if placeholder:
            self._path.setPlaceholderText(placeholder)
        self._path.setMinimumWidth(minimum_width)
        self._path.setClearButtonEnabled(clear_button)
        self._path.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        self._path.installEventFilter(self)
        self._path.textChanged.connect(self._on_text_changed)
        self._path.editingFinished.connect(self.path_changed.emit)

        if initial is not None and not self._path.text().strip():
            self.set_path(initial)

        self._browse_btn = QPushButton("Browse", self)
        self._browse_btn.setSizePolicy(
            QSizePolicy.Policy.Fixed,
            QSizePolicy.Policy.Fixed,
        )
        self._browse_btn.clicked.connect(self._on_browse_clicked)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        layout.addWidget(self._path, stretch=1)
        layout.addWidget(self._browse_btn)

        self._default_tooltip = self._path.toolTip()
        self._sync_path_tooltip(self._path.text())

    def line_edit(self) -> QLineEdit:
        return self._path

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if obj is self._path and event.type() == QEvent.Type.FocusIn:
            self.activated.emit(self)
        return super().eventFilter(obj, event)

    def _on_text_changed(self, text: str) -> None:
        self._sync_path_tooltip(text)
        self.text_changed.emit(text)

    def _sync_path_tooltip(self, text: str) -> None:
        tip = text.strip()
        self._path.setToolTip(tip if tip else self._default_tooltip)

    def set_placeholder(self, text: str) -> None:
        self._path.setPlaceholderText(text)

    def set_default_tooltip(self, text: str) -> None:
        self._default_tooltip = text
        self._sync_path_tooltip(self._path.text())

    def set_browse_title(self, title: str) -> None:
        self._browse_title = title

    def set_file_filter(self, file_filter: str) -> None:
        self._file_filter = file_filter

    def set_start_dir(self, start_dir: str) -> None:
        self._start_dir = start_dir

    def set_browse_handler(self, handler: Callable[[], None] | None) -> None:
        self._browse_handler = handler

    def set_path(self, path: Path | str | None, *, emit: bool = True) -> None:
        if path is None:
            self._path.clear()
        else:
            text = str(path)
            if isinstance(path, Path) or ("/" in text or "\\" in text):
                try:
                    text = str(Path(text).expanduser().resolve())
                except OSError:
                    text = str(Path(text).expanduser())
            self._path.setText(text)
        if emit:
            self.path_changed.emit()

    def path(self) -> Path | None:
        text = self._path.text().strip()
        if not text:
            return None
        return Path(text).expanduser()

    def text(self) -> str:
        return self._path.text()

    def set_text(self, text: str, *, emit: bool = True) -> None:
        self._path.setText(text)
        if emit:
            self.path_changed.emit()

    def clear(self) -> None:
        self._path.clear()
        self.path_changed.emit()

    def focus_field(self) -> None:
        self._path.setFocus(Qt.FocusReason.OtherFocusReason)

    def _on_browse_clicked(self) -> None:
        if self._browse_handler is not None:
            self._browse_handler()
            return
        self._browse()

    def _browse(self) -> None:
        start = self._path.text().strip() or self._start_dir or str(Path.home())
        if self._mode == "directory":
            chosen = QFileDialog.getExistingDirectory(
                self,
                self._browse_title,
                start,
            )
            if chosen:
                self._path.setText(chosen)
                self.path_changed.emit()
            return
        if self._mode == "file_or_directory":
            menu = QMenu(self)
            act_file = menu.addAction("FASTA file…")
            act_folder = menu.addAction("Folder…")
            chosen_act = menu.exec(
                self._browse_btn.mapToGlobal(self._browse_btn.rect().bottomLeft())
            )
            if chosen_act is act_folder:
                chosen = QFileDialog.getExistingDirectory(
                    self,
                    "Choose a folder of FASTA files",
                    start,
                )
                if chosen:
                    self._path.setText(chosen)
                    self.path_changed.emit()
                return
            if chosen_act is not act_file:
                return
        chosen, _ = QFileDialog.getOpenFileName(
            self,
            self._browse_title,
            start,
            self._file_filter,
        )
        if chosen:
            self._path.setText(chosen)
            self.path_changed.emit()