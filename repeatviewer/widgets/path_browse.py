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
    QVBoxLayout,
    QWidget,
)

PATH_FIELD_MIN_WIDTH = 280


class PathBrowseRow(QWidget):
    """Line edit + Browse matching the MSA FASTA path-filling style.

    Form labels belong on the caller (``form.addRow("MSA FASTA:", row)``).

    ``browse_below``: put Browse under the path field (narrow side panels).
    ``elide_middle``: when unfocused, show only the file/folder name (last path
    component); hover tooltip and focus show the full path. ``path()`` always
    returns the full path.
    """

    path_changed = Signal()
    activated = Signal(object)
    text_changed = Signal(str)
    extra_clicked = Signal()

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
        browse_below: bool = False,
        elide_middle: bool = False,
        extra_button_label: str = "",
        extra_button_tooltip: str = "",
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
        self._elide_middle = bool(elide_middle)
        self._full_text = ""
        self._suppress_text_sync = False
        # True while the field shows the real path for editing (focused).
        self._editing_full = False

        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        # Container stays out of the tab chain; path edits only on click so
        # opening a tool dialog does not park the cursor in Reference / MSA.
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)

        self._path = line_edit if line_edit is not None else QLineEdit(self)
        if self._path.parent() is None:
            self._path.setParent(self)
        self._path.setObjectName("PathField")
        self._path.setFocusPolicy(Qt.FocusPolicy.ClickFocus)
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
        self._path.editingFinished.connect(self._on_editing_finished)

        if initial is not None and not self._path.text().strip():
            self.set_path(initial)

        self._extra_btn: QPushButton | None = None
        label = str(extra_button_label or "").strip()
        if label:
            self._extra_btn = QPushButton(label, self)
            self._extra_btn.setSizePolicy(
                QSizePolicy.Policy.Fixed,
                QSizePolicy.Policy.Fixed,
            )
            if extra_button_tooltip:
                self._extra_btn.setToolTip(extra_button_tooltip)
            self._extra_btn.clicked.connect(self.extra_clicked.emit)

        self._browse_btn = QPushButton("Browse", self)
        self._browse_btn.setSizePolicy(
            QSizePolicy.Policy.Fixed,
            QSizePolicy.Policy.Fixed,
        )
        self._browse_btn.clicked.connect(self._on_browse_clicked)

        if browse_below:
            layout = QVBoxLayout(self)
            layout.setContentsMargins(0, 0, 0, 0)
            layout.setSpacing(6)
            layout.addWidget(self._path)
            browse_row = QHBoxLayout()
            browse_row.setContentsMargins(0, 0, 0, 0)
            browse_row.setSpacing(6)
            if self._extra_btn is not None:
                browse_row.addWidget(self._extra_btn, 0)
            browse_row.addWidget(self._browse_btn, 0)
            browse_row.addStretch(1)
            layout.addLayout(browse_row)
        else:
            layout = QHBoxLayout(self)
            layout.setContentsMargins(0, 0, 0, 0)
            layout.setSpacing(6)
            layout.addWidget(self._path, stretch=1)
            if self._extra_btn is not None:
                layout.addWidget(self._extra_btn)
            layout.addWidget(self._browse_btn)

        self._default_tooltip = self._path.toolTip()
        self._sync_path_tooltip(self._full_text if self._elide_middle else self._path.text())
        if self._elide_middle:
            self._apply_display()

    def line_edit(self) -> QLineEdit:
        return self._path

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if obj is self._path:
            et = event.type()
            if et == QEvent.Type.FocusIn:
                self.activated.emit(self)
                if self._elide_middle:
                    self._editing_full = True
                    self._apply_display()
            elif et == QEvent.Type.FocusOut:
                if self._elide_middle:
                    # Only commit when the field was showing the real path. Never
                    # replace ``_full_text`` with an already-elided display string.
                    if self._editing_full:
                        self._full_text = self._path.text()
                    self._editing_full = False
                    self._apply_display()
                else:
                    # Prefer the filename end of long absolute paths when the
                    # field is not focused (Mapping Inputs / Reference, etc.).
                    self._reveal_path_end()
            elif et == QEvent.Type.Show:
                if not self._path.hasFocus():
                    self._reveal_path_end()
            elif et == QEvent.Type.Resize and self._elide_middle:
                if not self._editing_full and not self._path.hasFocus():
                    self._apply_display()
            elif et == QEvent.Type.Resize and not self._elide_middle:
                if not self._path.hasFocus():
                    self._reveal_path_end()
        return super().eventFilter(obj, event)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        if self._elide_middle and not self._path.hasFocus():
            self._apply_display()
        elif not self._elide_middle and not self._path.hasFocus():
            self._reveal_path_end()

    def _reveal_path_end(self) -> None:
        """Keep the end of a long path visible when the field is not edited."""
        if self._elide_middle:
            if not self._editing_full and not self._path.hasFocus():
                self._apply_display()
            return
        if self._path.hasFocus():
            return
        text = self._path.text()
        if text:
            self._path.setCursorPosition(len(text))

    def _on_editing_finished(self) -> None:
        if self._elide_middle and self._editing_full:
            self._full_text = self._path.text()
        self.path_changed.emit()

    def _on_text_changed(self, text: str) -> None:
        if self._suppress_text_sync:
            self._sync_path_tooltip(
                self._full_text if self._elide_middle else text
            )
            return
        if self._elide_middle:
            # Clear-button / typing while editing updates the real path.
            # Ignore display-only elided setText (handled via suppress).
            if self._editing_full or not text or text == self._full_text:
                self._full_text = text
            self._sync_path_tooltip(self._full_text)
            self.text_changed.emit(self._full_text)
        else:
            self._sync_path_tooltip(text)
            self.text_changed.emit(text)

    def _elided_display(self, text: str) -> str:
        """Show only the last path component; middle-elide that name if needed."""
        if not text:
            return ""
        try:
            name = Path(text).name
        except (OSError, ValueError, TypeError):
            name = ""
        if not name:
            stripped = text.rstrip("/\\")
            cut = max(stripped.rfind("/"), stripped.rfind("\\"))
            name = stripped[cut + 1 :] if cut >= 0 else stripped
        if not name:
            name = text
        fm = self._path.fontMetrics()
        # PathField uses padding 1px 6px; clear button sits in the right pad.
        # Underestimating here makes the name slightly too wide and QLineEdit
        # scrolls, clipping the first character.
        side_pad = 12  # left + right stylesheet padding
        clear_pad = 24 if self._path.isClearButtonEnabled() else 0
        width = max(24, self._path.width() - side_pad - clear_pad - 4)
        if fm.horizontalAdvance(name) <= width:
            return name
        return fm.elidedText(name, Qt.TextElideMode.ElideMiddle, width)

    def _apply_display(self) -> None:
        if not self._elide_middle:
            return
        editing = self._editing_full or self._path.hasFocus()
        if editing:
            shown = self._full_text
        else:
            shown = self._elided_display(self._full_text)
        if self._path.text() != shown:
            self._suppress_text_sync = True
            try:
                self._path.setText(shown)
            finally:
                self._suppress_text_sync = False
        # After shrinking from a long full path, reset scroll so the first
        # character of the filename is not clipped off the left edge.
        if not editing:
            self._path.setCursorPosition(0)
        self._sync_path_tooltip(self._full_text)

    def _sync_path_tooltip(self, text: str) -> None:
        tip = text.strip()
        self._path.setToolTip(tip if tip else self._default_tooltip)

    def set_placeholder(self, text: str) -> None:
        self._path.setPlaceholderText(text)

    def set_default_tooltip(self, text: str) -> None:
        self._default_tooltip = text
        tip_src = self._full_text if self._elide_middle else self._path.text()
        self._sync_path_tooltip(tip_src)

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
            self._full_text = ""
            self._suppress_text_sync = True
            try:
                self._path.clear()
            finally:
                self._suppress_text_sync = False
        else:
            text = str(path)
            if isinstance(path, Path) or ("/" in text or "\\" in text):
                try:
                    text = str(Path(text).expanduser().resolve())
                except OSError:
                    text = str(Path(text).expanduser())
            self._full_text = text
            if self._elide_middle:
                self._apply_display()
            else:
                self._path.setText(text)
                self._reveal_path_end()
        if emit:
            self.path_changed.emit()

    def path(self) -> Path | None:
        text = (
            self._full_text.strip()
            if self._elide_middle
            else self._path.text().strip()
        )
        if not text:
            return None
        path = Path(text).expanduser()
        try:
            return path.resolve()
        except OSError:
            return path

    def text(self) -> str:
        if self._elide_middle:
            return self._full_text
        return self._path.text()

    def set_text(self, text: str, *, emit: bool = True) -> None:
        self._full_text = text
        if self._elide_middle:
            self._apply_display()
        else:
            self._path.setText(text)
            self._reveal_path_end()
        if emit:
            self.path_changed.emit()

    def clear(self) -> None:
        self._full_text = ""
        self._suppress_text_sync = True
        try:
            self._path.clear()
        finally:
            self._suppress_text_sync = False
        self.path_changed.emit()

    def focus_field(self) -> None:
        self._path.setFocus(Qt.FocusReason.OtherFocusReason)

    def _on_browse_clicked(self) -> None:
        if self._browse_handler is not None:
            self._browse_handler()
            return
        self._browse()

    def _browse_start_dir(self) -> str:
        raw = (self.text().strip() if self._elide_middle else self._path.text().strip())
        raw = raw or self._start_dir or str(Path.home())
        try:
            path = Path(raw).expanduser()
            if path.is_file():
                return str(path.parent)
            if path.is_dir():
                return str(path)
            parent = path.parent
            if parent.is_dir():
                return str(parent)
        except OSError:
            pass
        return raw

    def _browse_file(self, start: str) -> None:
        chosen, _ = QFileDialog.getOpenFileName(
            self,
            self._browse_title,
            start,
            self._file_filter,
        )
        if chosen:
            self.set_path(chosen)

    def _browse_directory(self, start: str, *, title: str | None = None) -> None:
        chosen = QFileDialog.getExistingDirectory(
            self,
            title or self._browse_title,
            start,
        )
        if chosen:
            self.set_path(chosen)

    def _browse(self) -> None:
        start = self._browse_start_dir()
        if self._mode == "directory":
            self._browse_directory(start)
            return
        if self._mode == "file_or_directory":
            # Menu under Browse: pick a FASTA file or a folder of FASTAs.
            menu = QMenu(self)
            act_file = menu.addAction("FASTA file...")
            act_file.setData("file")
            act_folder = menu.addAction("Folder...")
            act_folder.setData("directory")
            chosen_act = menu.exec(
                self._browse_btn.mapToGlobal(self._browse_btn.rect().bottomLeft())
            )
            if chosen_act is None:
                return
            kind = str(chosen_act.data() or "")
            if kind == "directory":
                self._browse_directory(
                    start, title="Choose a folder of FASTA files"
                )
                return
            if kind != "file":
                return
        self._browse_file(start)
