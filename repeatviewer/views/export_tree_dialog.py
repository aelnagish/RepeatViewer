"""Export dialog for the tree preview (MSA Annotation Export style).

Type + File (+ DPI for PNG). Cancel / OK.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings, Qt
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QMessageBox,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox
from repeatviewer.widgets.file_conflict import confirm_existing_file_path
from repeatviewer.widgets.help_button import field_with_help
from repeatviewer.widgets.path_browse import PathBrowseRow
from repeatviewer.widgets.spin_boxes import NoWheelSpinBox as QSpinBox

_FORMAT_SUFFIX = {
    "svg": ".svg",
    "png": ".png",
    "pdf": ".pdf",
}
_FORMAT_FILTER = {
    "svg": "SVG (*.svg);;All files (*)",
    "png": "PNG (*.png);;All files (*)",
    "pdf": "PDF (*.pdf);;All files (*)",
}
DPI_FORMATS = frozenset({"png"})

_TYPE_HELP = (
    "<b>SVG</b> and <b>PDF</b>: vector graphics of the tree as shown "
    "(sharp at any zoom).<br><br>"
    "<b>PNG</b>: a bitmap of the tree; uses the DPI setting."
)

_SETTINGS_GROUP = "tree_view_export"
_DEFAULT_DPI = 300
_MIN_DPI = 50
_MAX_DPI = 600


class ExportTreeDialog(QDialog):
    """Type + File (+ DPI for PNG), same chrome as Annotation Export."""

    def __init__(
        self,
        *,
        source_path: Path,
        layout_tag: str = "rectangular",
        default_format: str = "svg",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Export")
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, False)
        self.setMinimumWidth(480)
        self._source_path = source_path.expanduser()
        tag = str(layout_tag or "tree").strip() or "tree"
        self._layout_tag = "".join(c if c.isalnum() or c in "-_" else "_" for c in tag)
        self._start_dir = (
            self._source_path.parent
            if self._source_path.parent.is_dir()
            else Path.home()
        )
        self._chosen: Path | None = None
        self._dpi_value = _DEFAULT_DPI
        self._syncing_name = False

        settings = QSettings()
        settings.beginGroup(_SETTINGS_GROUP)
        saved_fmt = str(
            settings.value("format", default_format) or default_format
        ).lower()
        try:
            saved_dpi = int(settings.value("dpi", _DEFAULT_DPI))
        except (TypeError, ValueError):
            saved_dpi = _DEFAULT_DPI
        saved_dpi = max(_MIN_DPI, min(_MAX_DPI, saved_dpi))
        settings.endGroup()

        self._type = QComboBox(self)
        self._type.addItem("SVG", "svg")
        self._type.addItem("PDF", "pdf")
        self._type.addItem("PNG", "png")
        idx = self._type.findData(saved_fmt)
        if idx < 0:
            idx = self._type.findData(default_format)
        self._type.setCurrentIndex(max(0, idx))
        self._last_fmt = self._current_fmt()

        self._dpi = QSpinBox(self)
        self._dpi.setRange(_MIN_DPI, _MAX_DPI)
        self._dpi.setSingleStep(25)
        self._dpi.setValue(saved_dpi)
        self._dpi.setSuffix(" DPI")

        self._file_row = PathBrowseRow(
            self,
            placeholder="Choose export path, or Browse",
            browse_title="Export",
            start_dir=str(self._start_dir),
            mode="file",
            browse_handler=self._browse,
        )
        self._file = self._file_row.line_edit()
        self._file_row.set_default_tooltip(
            "Enter an export path, or use Browse to choose the file."
        )

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setSpacing(8)
        form.addRow(
            "Type:",
            field_with_help(self._type, "Export type", _TYPE_HELP, parent=self),
        )
        form.addRow("DPI:", self._dpi)
        form.addRow("File:", self._file_row)
        self._form = form

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=self,
        )
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)
        layout.addLayout(form)
        layout.addWidget(buttons)

        self._sync_filename_from_format()
        self._file_row.set_file_filter(self._file_filter())
        self._type.currentIndexChanged.connect(self._on_format_changed)
        self._sync_rows()

    def chosen_path(self) -> Path | None:
        return self._chosen

    def dpi(self) -> int:
        return int(self._dpi_value)

    def export_format(self) -> str:
        return self._current_fmt()

    def _current_fmt(self) -> str:
        fmt = str(self._type.currentData() or "svg")
        return fmt if fmt in _FORMAT_SUFFIX else "svg"

    def _file_filter(self) -> str:
        return _FORMAT_FILTER[self._current_fmt()]

    def _default_out_for(self, fmt: str) -> Path:
        return (
            self._source_path.parent
            / f"{self._source_path.name}_{self._layout_tag}{_FORMAT_SUFFIX[fmt]}"
        )

    def _default_out(self) -> Path:
        return self._default_out_for(self._current_fmt())

    def _sync_filename_from_format(self) -> None:
        self._syncing_name = True
        try:
            self._file.setText(str(self._default_out()))
        finally:
            self._syncing_name = False

    def _sync_rows(self) -> None:
        self._form.setRowVisible(self._dpi, self._current_fmt() in DPI_FORMATS)
        self.adjustSize()

    def _on_format_changed(self, _index: int = 0) -> None:
        previous = self._last_fmt
        fmt = self._current_fmt()
        self._last_fmt = fmt
        self._file_row.set_file_filter(self._file_filter())
        self._sync_rows()
        if self._syncing_name:
            return
        current = self._file.text().strip()
        if not current:
            self._sync_filename_from_format()
            return
        path = Path(current).expanduser()
        try:
            was_default = path.resolve() == self._default_out_for(previous).resolve()
        except OSError:
            was_default = False
        updated = (
            self._default_out_for(fmt)
            if was_default
            else path.with_suffix(_FORMAT_SUFFIX[fmt])
        )
        self._syncing_name = True
        try:
            self._file.setText(str(updated))
        finally:
            self._syncing_name = False

    def _browse(self) -> None:
        current = self._file.text().strip()
        start = str(Path(current).expanduser()) if current else str(self._default_out())
        chosen, _ = QFileDialog.getSaveFileName(
            self,
            self.windowTitle(),
            start,
            self._file_filter(),
        )
        if not chosen:
            return
        path = Path(chosen)
        suffix = _FORMAT_SUFFIX[self._current_fmt()]
        if path.suffix.lower() != suffix:
            path = path.with_suffix(suffix)
        self._file.setText(str(path))
        self._file_row.set_file_filter(self._file_filter())

    def _proposed_path(self) -> Path | None:
        text = self._file.text().strip()
        if not text:
            QMessageBox.warning(self, self.windowTitle(), "Choose an export file.")
            return None
        path = Path(text).expanduser()
        if path.name == text and "/" not in text and "\\" not in text:
            path = self._start_dir / path.name
        path = path.resolve()
        suffix = _FORMAT_SUFFIX[self._current_fmt()]
        if path.suffix.lower() != suffix:
            path = path.with_suffix(suffix)
        return path

    def _accept(self) -> None:
        fmt = self._current_fmt()
        out = self._proposed_path()
        if out is None:
            return
        confirmed = confirm_existing_file_path(
            self,
            out,
            title=self.windowTitle(),
            informative=(
                "Overwrite the existing file, or save under a different name?"
            ),
            required_suffix=_FORMAT_SUFFIX[fmt],
        )
        if confirmed is None:
            return
        self._chosen = confirmed
        self._dpi_value = int(self._dpi.value())
        settings = QSettings()
        settings.beginGroup(_SETTINGS_GROUP)
        settings.setValue("format", fmt)
        settings.setValue("dpi", self._dpi_value)
        settings.endGroup()
        self.accept()
