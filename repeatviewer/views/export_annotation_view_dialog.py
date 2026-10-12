"""Export annotation / overview view as SVG, PNG, or PDF (type + DPI + file)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings
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

_SETTINGS_GROUP = "annotation_view_export"
_DEFAULT_DPI = 400
_MIN_DPI = 50
_MAX_DPI = 600


def write_view_image_pdf(
    out_path: Path,
    img,
    *,
    dpi: int,
    content_width: int,
    content_height: int,
) -> None:
    """Write ``img`` to a PDF page matching the image (no large blank margins)."""
    from PySide6.QtCore import QMarginsF, QRect, QSizeF
    from PySide6.QtGui import QColor, QPageLayout, QPageSize, QPainter, QPdfWriter

    del content_width, content_height  # page follows the rendered image pixels
    writer = QPdfWriter(str(out_path))
    dpi = max(_MIN_DPI, min(_MAX_DPI, int(dpi)))
    writer.setResolution(dpi)
    # Size the page so device pixels match the image (no letterboxing / stretch).
    img_w = max(1, int(img.width()))
    img_h = max(1, int(img.height()))
    page_w_pt = img_w * 72.0 / float(dpi)
    page_h_pt = img_h * 72.0 / float(dpi)
    writer.setPageLayout(
        QPageLayout(
            QPageSize(QSizeF(page_w_pt, page_h_pt), QPageSize.Unit.Point),
            QPageLayout.Orientation.Portrait,
            QMarginsF(0, 0, 0, 0),
        )
    )
    painter = QPainter(writer)
    try:
        painter.fillRect(painter.viewport(), QColor(255, 255, 255))
        painter.drawImage(QRect(0, 0, img_w, img_h), img)
    finally:
        painter.end()


class ExportAnnotationViewDialog(QDialog):
    """Type + DPI + File layout matching Export annotation."""

    def __init__(
        self,
        *,
        source_path: Path,
        title: str = "Export annotation view",
        name_tag: str = "overview",
        default_format: str = "pdf",
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumWidth(480)
        self._source_path = source_path.expanduser()
        self._name_tag = str(name_tag or "overview").strip() or "overview"
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
        if saved_fmt not in _FORMAT_SUFFIX:
            saved_fmt = default_format if default_format in _FORMAT_SUFFIX else "pdf"
        try:
            saved_dpi = int(settings.value("dpi", _DEFAULT_DPI))
        except (TypeError, ValueError):
            saved_dpi = _DEFAULT_DPI
        saved_dpi = max(_MIN_DPI, min(_MAX_DPI, saved_dpi))
        settings.endGroup()

        self._type = QComboBox(self)
        self._type.addItem("PDF", "pdf")
        self._type.addItem("PNG", "png")
        self._type.addItem("SVG", "svg")
        idx = self._type.findData(saved_fmt)
        self._type.setCurrentIndex(max(0, idx))

        self._dpi = QSpinBox(self)
        self._dpi.setRange(_MIN_DPI, _MAX_DPI)
        self._dpi.setSingleStep(25)
        self._dpi.setValue(saved_dpi)
        self._dpi.setSuffix(" DPI")

        self._file_row = PathBrowseRow(
            self,
            placeholder="Choose export path, or Browse",
            browse_title=title,
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
        form.addRow("Type:", self._type)
        form.addRow("DPI:", self._dpi)
        form.addRow("File:", self._file_row)

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
        self._sync_dpi_enabled()

    def chosen_path(self) -> Path | None:
        return self._chosen

    def dpi(self) -> int:
        return int(self._dpi_value)

    def export_format(self) -> str:
        return str(self._type.currentData() or "svg")

    def _current_fmt(self) -> str:
        return str(self._type.currentData() or "svg")

    def _file_filter(self) -> str:
        return _FORMAT_FILTER[self._current_fmt()]

    def _default_out(self) -> Path:
        """``{source.name}_{tag}.pdf`` — keeps the original FASTA extension."""
        suffix = _FORMAT_SUFFIX[self._current_fmt()]
        return (
            self._source_path.parent
            / f"{self._source_path.name}_{self._name_tag}{suffix}"
        )

    def _sync_filename_from_format(self) -> None:
        self._syncing_name = True
        try:
            self._file.setText(str(self._default_out()))
        finally:
            self._syncing_name = False

    def _sync_dpi_enabled(self) -> None:
        # DPI applies to all three (PNG pixels, PDF/SVG resolution metadata).
        self._dpi.setEnabled(True)

    def _on_format_changed(self, _index: int = 0) -> None:
        self._file_row.set_file_filter(self._file_filter())
        self._sync_dpi_enabled()
        if self._syncing_name:
            return
        current = self._file.text().strip()
        suffix = _FORMAT_SUFFIX[self._current_fmt()]
        if current:
            path = Path(current).expanduser()
            updated = path.with_suffix(suffix)
            self._syncing_name = True
            try:
                self._file.setText(str(updated))
            finally:
                self._syncing_name = False
        else:
            self._sync_filename_from_format()

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
        out = self._proposed_path()
        if out is None:
            return
        suffix = _FORMAT_SUFFIX[self._current_fmt()]
        confirmed = confirm_existing_file_path(
            self,
            out,
            title=self.windowTitle(),
            informative=(
                "Overwrite the existing file, or save under a different name?"
            ),
            required_suffix=suffix,
        )
        if confirmed is None:
            return
        self._chosen = confirmed
        self._dpi_value = int(self._dpi.value())
        settings = QSettings()
        settings.beginGroup(_SETTINGS_GROUP)
        settings.setValue("format", self._current_fmt())
        settings.setValue("dpi", self._dpi_value)
        settings.endGroup()
        self.accept()
