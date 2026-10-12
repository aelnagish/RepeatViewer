"""Create-folder dialog: pick the destination, name it, validate before mkdir."""

from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import QSize, Qt
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

DEFAULT_FOLDER_NAME = "untitled folder"

# ``:`` is illegal in Finder, NUL is illegal everywhere, and a backslash in a
# segment is almost always a Windows path pasted by mistake.
_ILLEGAL_SEGMENT_CHARS = frozenset(':\\\0')
_RESERVED_SEGMENTS = frozenset({".", ".."})


def expand_name_tokens(name: str, *, now: datetime | None = None) -> str:
    """Expand ``{date}`` / ``{time}`` in a folder name (run folders are a habit here)."""
    stamp = now or datetime.now()
    return (
        name.replace("{date}", stamp.strftime("%Y-%m-%d"))
        .replace("{time}", stamp.strftime("%H%M"))
    )


def unique_folder_name(parent: Path, name: str) -> str:
    """``name`` if free, else ``name 2``, ``name 3``, … — mirrors the duplicate-file style."""
    if not (parent / name).exists():
        return name
    n = 2
    while (parent / f"{name} {n}").exists():
        n += 1
    return f"{name} {n}"


def common_prefix_name(paths: list[Path], *, minimum: int = 3) -> str:
    """Shared leading text of several file names, e.g. the four ``CroSat2_2*`` outputs.

    Returns ``""`` when the selection has no meaningful prefix, so the caller can
    fall back to the generic default rather than proposing a one-letter folder.
    """
    stems = [p.stem for p in paths if p.name]
    if len(stems) < 2:
        return stems[0].strip(" ._-") if stems else ""
    prefix = os.path.commonprefix(stems).strip(" ._-")
    return prefix if len(prefix) >= minimum else ""


def validate_folder_name(parent: Path | None, raw_name: str) -> tuple[bool, str]:
    """Check a destination + name pair before ``mkdir``.

    Returns ``(ok, message)``. Collision is reported so the accept handler can
    show the folder-exists reminder (no silent overwrite).
    """
    if parent is None:
        return False, "Choose a folder to create under."
    if not parent.is_dir():
        return False, f"Not a folder: {parent}"
    if not os.access(parent, os.W_OK):
        return False, f"No write permission in {parent.name}."

    raw = raw_name.strip()
    if not raw:
        return False, "Enter a folder name."

    name = expand_name_tokens(raw)
    if name.startswith(("/", "~")):
        return False, "Enter a name, not an absolute path."

    segments = [seg.strip() for seg in name.split("/")]
    if any(not seg for seg in segments):
        return False, "Remove the empty path segment."
    for seg in segments:
        if seg in _RESERVED_SEGMENTS:
            return False, f"“{seg}” is not a usable folder name."
        bad = sorted(_ILLEGAL_SEGMENT_CHARS.intersection(seg))
        if bad:
            shown = " ".join(repr(c) for c in bad)
            return False, f"Remove {shown} from the name."

    target = parent / name
    if target.exists():
        return False, "exists"
    if len(segments) > 1:
        return True, f"Creates {len(segments)} nested folders."
    if name != raw:
        return True, f"Creates “{name}”."
    return True, ""


class CreateFolderDialog(QDialog):
    """Ask for a destination and a folder name, with live validation."""

    # Preferred open size (width from screenshot; height fits content).
    _DEFAULT_WIDTH = 920

    def __init__(
        self,
        parent: QWidget | None,
        *,
        initial_destination: Path | None = None,
        destinations: list[Path] | None = None,
        default_name: str = DEFAULT_FOLDER_NAME,
        copy_count: int = 0,
        move_count: int = 0,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(
            "New Folder with Selected Items" if (copy_count or move_count) else "Create folder"
        )
        self.setMinimumWidth(640)
        self.setSizeGripEnabled(True)
        # Non-modal so Folder-panel clicks can fill Create under.
        self.setModal(False)
        self.setWindowModality(Qt.WindowModality.NonModal)
        # Parent stacking keeps us above RepeatViewer; never global always-on-top.
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, False)

        start = initial_destination
        if start is None and destinations:
            start = destinations[0]
        self._destination: Path | None = None
        self._did_center = False

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 14, 18, 14)
        layout.setSpacing(10)
        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setVerticalSpacing(8)

        self._destination_edit = QLineEdit(self)
        self._destination_edit.setMinimumWidth(420)
        self._destination_edit.setPlaceholderText("Folder path")
        self._destination_edit.setToolTip(
            "Click to edit the path, or use Browse… / click a folder in the Folders panel."
        )
        self._destination_edit.textChanged.connect(self._on_destination_text_changed)
        self._destination_edit.editingFinished.connect(self._on_destination_edited)
        browse_btn = QPushButton("Browse…")
        browse_btn.setAutoDefault(False)
        browse_btn.clicked.connect(self._browse_for_destination)
        destination_row = QHBoxLayout()
        destination_row.setContentsMargins(0, 0, 0, 0)
        destination_row.addWidget(self._destination_edit, 1)
        destination_row.addWidget(browse_btn, 0)
        destination_host = QWidget()
        destination_host.setLayout(destination_row)
        form.addRow("Create under:", destination_host)

        self._name_edit = QLineEdit(default_name)
        self._name_edit.setMinimumWidth(260)
        self._name_edit.textChanged.connect(self._revalidate)
        form.addRow("New folder name:", self._name_edit)
        layout.addLayout(form)

        count = copy_count or move_count
        if count:
            noun = "file" if count == 1 else "files"
            verb = "copied" if copy_count else "moved"
            note = QLabel(
                f"{count} selected {noun} will be {verb} into the new folder."
            )
            note.setWordWrap(True)
            layout.addWidget(note)

        self._status = QLabel("")
        self._status.setWordWrap(True)
        self._status.setTextFormat(Qt.TextFormat.PlainText)
        # Keep status from reserving a tall empty band when unused.
        self._status.setVisible(False)
        layout.addWidget(self._status)

        self._buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self._buttons.accepted.connect(self._accept)
        self._buttons.rejected.connect(self.reject)
        layout.addWidget(self._buttons)

        if start is not None:
            self.set_destination(start)
        self._revalidate()
        self.adjustSize()
        self.resize(self._DEFAULT_WIDTH, self.sizeHint().height())
        self._name_edit.setFocus()
        self._name_edit.selectAll()

    def sizeHint(self) -> QSize:  # noqa: N802
        hint = super().sizeHint()
        return QSize(max(hint.width(), self._DEFAULT_WIDTH), hint.height())

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        if self._did_center:
            return
        self._did_center = True
        # Fit height to content; keep the preferred width.
        self.adjustSize()
        self.resize(self._DEFAULT_WIDTH, self.sizeHint().height())
        parent = self.parentWidget()
        if parent is not None:
            geo = parent.frameGeometry()
            x = geo.x() + max(0, (geo.width() - self.width()) // 2)
            y = geo.y() + max(0, (geo.height() - self.height()) // 2)
            self.move(x, y)

    # ----- destinations -------------------------------------------------

    def set_destination(self, folder: Path, *, quiet: bool = False) -> None:
        """Set Create under from code (Browse… or Folders panel)."""
        try:
            folder = folder.expanduser().resolve()
        except OSError:
            folder = folder.expanduser()
        new_dest = folder if folder.is_dir() else None
        text = str(folder) if new_dest is not None else str(folder)
        if self._destination == new_dest and self._destination_edit.text() == text:
            return
        self._destination = new_dest
        if self._destination_edit.text() != text:
            self._destination_edit.blockSignals(True)
            self._destination_edit.setText(text)
            self._destination_edit.blockSignals(False)
        self._destination_edit.setToolTip(text)
        self._revalidate()

    def apply_folder_from_panel(self, folder: Path) -> bool:
        """Fill Create under when the user clicks a folder in the Folders panel."""
        try:
            resolved = folder.expanduser().resolve()
        except OSError:
            return False
        if not resolved.is_dir():
            return False
        self.set_destination(resolved, quiet=True)
        self.raise_()
        self.activateWindow()
        return True

    def _browse_for_destination(self) -> None:
        start = self.destination()
        picked = QFileDialog.getExistingDirectory(
            self,
            "Create folder under",
            str(start) if start is not None else str(Path.home()),
        )
        if not picked:
            return
        self.set_destination(Path(picked))

    def _on_destination_text_changed(self, _text: str = "") -> None:
        # Live-update while typing; full resolve on editingFinished.
        raw = self._destination_edit.text().strip()
        if not raw:
            self._destination = None
            self._revalidate()
            return
        candidate = Path(raw).expanduser()
        self._destination = candidate if candidate.is_dir() else None
        self._revalidate()

    def _on_destination_edited(self) -> None:
        text = self._destination_edit.text().strip()
        if not text:
            self._destination = None
            self._revalidate()
            return
        self.set_destination(Path(text))

    def destination(self) -> Path | None:
        return self._destination

    # ----- name and validation -----------------------------------------

    def folder_name(self) -> str:
        return expand_name_tokens(self._name_edit.text().strip())

    def target_path(self) -> Path | None:
        parent = self.destination()
        name = self.folder_name()
        if parent is None or not name:
            return None
        return parent / name

    def _revalidate(self) -> None:
        ok, message = self._validate()
        if message and message != "exists":
            self._status.setText(message)
            self._status.setVisible(True)
        else:
            self._status.clear()
            self._status.setVisible(False)
        button = self._buttons.button(QDialogButtonBox.StandardButton.Ok)
        if button is not None:
            button.setEnabled(ok or message == "exists")

    def _validate(self) -> tuple[bool, str]:
        return validate_folder_name(self.destination(), self._name_edit.text())

    def _accept(self) -> None:
        # Commit any in-progress path edit first.
        self._on_destination_edited()
        ok, message = self._validate()
        if message == "exists":
            name = self.folder_name()
            leaf = name.split("/")[-1] if name else self._name_edit.text().strip()
            QMessageBox.warning(
                self,
                self.windowTitle(),
                f"“{leaf}” folder name already exist.\n\n"
                "Please choose a different folder name.",
            )
            self._name_edit.selectAll()
            self._name_edit.setFocus()
            return
        if not ok:
            if message:
                QMessageBox.warning(self, self.windowTitle(), message)
            return
        self.accept()
