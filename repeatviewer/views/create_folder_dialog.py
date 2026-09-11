"""Create-folder dialog: pick the destination, name it, validate before mkdir.

The old flow was a bare ``QInputDialog.getText`` whose parent folder came from
whatever tree row happened to be under the mouse — invisible and uncorrectable —
and which reported collisions as a raw ``[Errno 17]`` after the fact. This dialog
shows the destination, lets you change it, and validates while you type.
"""

from __future__ import annotations

import os
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
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


def validate_folder_name(parent: Path | None, raw_name: str) -> tuple[bool, str, str]:
    """Check a destination + name pair before ``mkdir``.

    Returns ``(ok, message, suggestion)``. ``suggestion`` is non-empty only for a
    name collision, where the caller can offer the next free name instead.
    """
    if parent is None:
        return False, "Choose a folder to create in.", ""
    if not parent.is_dir():
        return False, f"Not a folder: {parent}", ""
    if not os.access(parent, os.W_OK):
        return False, f"No write permission in {parent.name}.", ""

    raw = raw_name.strip()
    if not raw:
        return False, "Enter a folder name.", ""

    name = expand_name_tokens(raw)
    if name.startswith(("/", "~")):
        return False, "Enter a name, not an absolute path.", ""

    segments = [seg.strip() for seg in name.split("/")]
    if any(not seg for seg in segments):
        return False, "Remove the empty path segment.", ""
    for seg in segments:
        if seg in _RESERVED_SEGMENTS:
            return False, f"“{seg}” is not a usable folder name.", ""
        bad = sorted(_ILLEGAL_SEGMENT_CHARS.intersection(seg))
        if bad:
            shown = " ".join(repr(c) for c in bad)
            return False, f"Remove {shown} from the name.", ""

    target = parent / name
    if target.exists():
        leaf = segments[-1]
        suggestion = "/".join([*segments[:-1], unique_folder_name(target.parent, leaf)])
        return False, f"“{leaf}” already exists here.", suggestion

    if len(segments) > 1:
        return True, f"Creates {len(segments)} nested folders.", ""
    if name != raw:
        return True, f"Creates “{name}”.", ""
    return True, "", ""


class CreateFolderDialog(QDialog):
    """Ask for a destination and a folder name, with live validation."""

    def __init__(
        self,
        parent: QWidget | None,
        *,
        destinations: list[Path],
        default_name: str = DEFAULT_FOLDER_NAME,
        move_count: int = 0,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Create folder")
        self._destinations: list[Path] = []

        layout = QVBoxLayout(self)
        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)

        self._destination_combo = QComboBox()
        self._destination_combo.setMinimumWidth(260)
        self._destination_combo.currentIndexChanged.connect(self._revalidate)
        browse_btn = QPushButton("Choose…")
        browse_btn.setAutoDefault(False)
        browse_btn.clicked.connect(self._browse_for_destination)
        destination_row = QHBoxLayout()
        destination_row.setContentsMargins(0, 0, 0, 0)
        destination_row.addWidget(self._destination_combo, 1)
        destination_row.addWidget(browse_btn, 0)
        destination_host = QWidget()
        destination_host.setLayout(destination_row)
        form.addRow("Create in:", destination_host)

        self._name_edit = QLineEdit(default_name)
        self._name_edit.setMinimumWidth(260)
        self._name_edit.textChanged.connect(self._revalidate)
        form.addRow("Name:", self._name_edit)
        layout.addLayout(form)

        self._hint = QLabel("{date} and {time} expand when the folder is created.")
        self._hint.setEnabled(False)
        layout.addWidget(self._hint)

        if move_count:
            noun = "file" if move_count == 1 else "files"
            moved = QLabel(f"{move_count} selected {noun} will be moved into the new folder.")
            moved.setWordWrap(True)
            layout.addWidget(moved)

        self._status = QLabel("")
        self._status.setWordWrap(True)
        self._status.setTextFormat(Qt.TextFormat.PlainText)
        layout.addWidget(self._status)

        self._suggest_btn = QPushButton("Use suggested name")
        self._suggest_btn.setAutoDefault(False)
        self._suggest_btn.setVisible(False)
        self._suggest_btn.clicked.connect(self._apply_suggestion)
        layout.addWidget(self._suggest_btn, 0, Qt.AlignmentFlag.AlignLeft)

        self._buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self._buttons.accepted.connect(self.accept)
        self._buttons.rejected.connect(self.reject)
        layout.addWidget(self._buttons)

        self.set_destinations(destinations)
        self._suggestion = ""
        self._revalidate()
        self._name_edit.setFocus()
        self._name_edit.selectAll()

    # ----- destinations -------------------------------------------------

    def set_destinations(self, destinations: list[Path]) -> None:
        seen: set[str] = set()
        ordered: list[Path] = []
        for raw in destinations:
            if raw is None:
                continue
            folder = Path(raw).expanduser()
            try:
                folder = folder.resolve()
            except OSError:
                continue
            if not folder.is_dir():
                continue
            key = str(folder)
            if key in seen:
                continue
            seen.add(key)
            ordered.append(folder)

        self._destinations = ordered
        self._destination_combo.blockSignals(True)
        self._destination_combo.clear()
        for folder in ordered:
            self._destination_combo.addItem(self._destination_label(folder), str(folder))
            self._destination_combo.setItemData(
                self._destination_combo.count() - 1,
                str(folder),
                Qt.ItemDataRole.ToolTipRole,
            )
        self._destination_combo.blockSignals(False)

    @staticmethod
    def _destination_label(folder: Path) -> str:
        try:
            shown = folder.relative_to(Path.home())
            return f"{folder.name}  —  ~/{shown.parent}" if shown.parent != Path(".") else folder.name
        except ValueError:
            return f"{folder.name}  —  {folder.parent}"

    def _browse_for_destination(self) -> None:
        start = self.destination()
        picked = QFileDialog.getExistingDirectory(
            self,
            "Create folder in",
            str(start) if start is not None else str(Path.home()),
        )
        if not picked:
            return
        folder = Path(picked).expanduser().resolve()
        existing = [Path(self._destination_combo.itemData(i)) for i in range(self._destination_combo.count())]
        self.set_destinations([folder, *existing])
        self._destination_combo.setCurrentIndex(0)
        self._revalidate()

    def destination(self) -> Path | None:
        raw = self._destination_combo.currentData()
        return Path(str(raw)) if raw else None

    # ----- name and validation -----------------------------------------

    def folder_name(self) -> str:
        return expand_name_tokens(self._name_edit.text().strip())

    def target_path(self) -> Path | None:
        parent = self.destination()
        name = self.folder_name()
        if parent is None or not name:
            return None
        return parent / name

    def _apply_suggestion(self) -> None:
        if self._suggestion:
            self._name_edit.setText(self._suggestion)
            self._name_edit.setFocus()
            self._name_edit.selectAll()

    def _revalidate(self) -> None:
        ok, message, suggestion = self._validate()
        self._suggestion = suggestion
        self._status.setText(message)
        self._suggest_btn.setVisible(bool(suggestion))
        if suggestion:
            self._suggest_btn.setText(f"Use “{suggestion}”")
        button = self._buttons.button(QDialogButtonBox.StandardButton.Ok)
        if button is not None:
            button.setEnabled(ok)

    def _validate(self) -> tuple[bool, str, str]:
        return validate_folder_name(self.destination(), self._name_edit.text())
