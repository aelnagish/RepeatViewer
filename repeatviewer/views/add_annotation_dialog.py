"""Dialog to pick an annotation file (GFF / BED / GTF) for a FASTA."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QLineEdit,
    QMessageBox,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.path_browse import PathBrowseRow

_ANNOTATION_TYPES: tuple[tuple[str, tuple[str, ...]], ...] = (
    ("GFF", (".gff", ".gff3")),
    ("BED", (".bed",)),
    ("GTF", (".gtf",)),
)


def annotation_type_for_path(path: Path) -> str | None:
    """Return GFF / BED / GTF for a path suffix, or None if unsupported."""
    suffix = path.suffix.lower()
    for label, suffixes in _ANNOTATION_TYPES:
        if suffix in suffixes:
            return label
    return None


class AddAnnotationDialog(QDialog):
    """Choose annotation type + file path (Browse, paste, or Files-panel click)."""

    def __init__(
        self,
        *,
        fasta_path: Path,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Add annotation")
        self.setMinimumWidth(480)
        self.setModal(False)
        self.setWindowModality(Qt.WindowModality.NonModal)
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, True)
        self._fasta_path = fasta_path.expanduser().resolve()
        self._start_dir = (
            self._fasta_path.parent if self._fasta_path.parent.is_dir() else Path.home()
        )
        self._result_path: Path | None = None

        self._type = QComboBox(self)
        for label, _suffixes in _ANNOTATION_TYPES:
            self._type.addItem(label, label)
        self._type.currentIndexChanged.connect(self._on_type_changed)

        self._file_row = PathBrowseRow(
            self,
            placeholder="Click a file in Files, or Browse",
            browse_title="Choose annotation file",
            start_dir=str(self._start_dir),
            mode="file",
            browse_handler=self._browse,
        )
        self._file = self._file_row.line_edit()
        self._file_row.set_default_tooltip(
            "Click a file in the Files panel to fill this path, or use Browse."
        )

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setSpacing(8)
        form.addRow("Type:", self._type)
        form.addRow("File:", self._file_row)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self._on_ok)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)
        layout.addLayout(form)
        layout.addWidget(buttons)

        self._suggest_sibling_annotation()

    def fasta_path(self) -> Path:
        return self._fasta_path

    def result_path(self) -> Path | None:
        return self._result_path

    def apply_path_from_panel(self, path: Path) -> bool:
        """Fill File from a Files-panel click; returns True if applied."""
        resolved = path.expanduser().resolve()
        if not resolved.is_file():
            return False
        if resolved == self._fasta_path:
            return False
        label = annotation_type_for_path(resolved)
        if label is None:
            return False
        self._set_type_label(label)
        self._file.setText(str(resolved))
        return True

    def _set_type_label(self, label: str) -> None:
        index = self._type.findData(label)
        if index < 0:
            return
        self._type.blockSignals(True)
        self._type.setCurrentIndex(index)
        self._type.blockSignals(False)

    def _selected_suffixes(self) -> tuple[str, ...]:
        label = str(self._type.currentData() or "GFF")
        for name, suffixes in _ANNOTATION_TYPES:
            if name == label:
                return suffixes
        return _ANNOTATION_TYPES[0][1]

    def _file_filter(self) -> str:
        label = str(self._type.currentData() or "GFF")
        suffixes = self._selected_suffixes()
        patterns = " ".join(f"*{suffix}" for suffix in suffixes)
        return f"{label} ({patterns});;All files (*)"

    def _on_type_changed(self, _index: int = 0) -> None:
        text = self._file.text().strip()
        if not text:
            return
        path = Path(text)
        suffixes = self._selected_suffixes()
        if path.suffix.lower() not in suffixes:
            self._file.clear()

    def _suggest_sibling_annotation(self) -> None:
        stem = self._fasta_path.stem
        parent = self._start_dir
        for label, suffixes in _ANNOTATION_TYPES:
            for suffix in suffixes:
                candidate = parent / f"{stem}{suffix}"
                if candidate.is_file():
                    self._set_type_label(label)
                    self._file.setText(str(candidate))
                    return

    def _browse(self) -> None:
        start = self._file.text().strip()
        start_dir = str(Path(start).expanduser().parent) if start else str(self._start_dir)
        chosen, _ = QFileDialog.getOpenFileName(
            self,
            "Choose annotation file",
            start_dir,
            self._file_filter(),
        )
        if not chosen:
            return
        path = Path(chosen)
        label = annotation_type_for_path(path)
        if label is not None:
            self._set_type_label(label)
        self._file.setText(str(path))
        self._file_row.set_file_filter(self._file_filter())

    def _on_ok(self) -> None:
        text = self._file.text().strip()
        if not text:
            QMessageBox.warning(self, "Add annotation", "Choose an annotation file.")
            return
        path = Path(text).expanduser().resolve()
        if not path.is_file():
            QMessageBox.warning(self, "Add annotation", f"File not found:\n{path}")
            return
        suffixes = self._selected_suffixes()
        if path.suffix.lower() not in suffixes:
            # If the path is a known annotation type, sync Type instead of blocking.
            label = annotation_type_for_path(path)
            if label is not None:
                self._set_type_label(label)
            else:
                selected = str(self._type.currentData() or "GFF")
                expected = ", ".join(suffixes)
                QMessageBox.warning(
                    self,
                    "Add annotation",
                    f"Selected type is {selected} ({expected}), but the file is "
                    f"“{path.name}”.\n\nChoose a matching file or change Type.",
                )
                return
        self._result_path = path
        self.accept()
