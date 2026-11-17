"""Export dialog for the Annotation / Overview / Dotplot views.

Image types: PDF and SVG (vector) and PNG (bitmap, the only one that uses DPI).
When the caller passes ``annotation_provider`` the same dialog also offers the
annotation file types GFF3 / GTF / BED, so one Export button covers both.
"""

from __future__ import annotations

from collections.abc import Callable, Iterator
from contextlib import contextmanager
from pathlib import Path

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QMessageBox,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox
from repeatviewer.widgets.combo_boxes import size_export_format_combo
from repeatviewer.widgets.file_conflict import confirm_existing_file_path
from repeatviewer.widgets.path_browse import PathBrowseRow
from repeatviewer.widgets.spin_boxes import NoWheelSpinBox as QSpinBox

_FORMAT_SUFFIX = {
    "svg": ".svg",
    "png": ".png",
    "pdf": ".pdf",
    "gff3": ".gff3",
    "gtf": ".gtf",
    "bed": ".bed",
}

_FORMAT_FILTER = {
    "svg": "SVG (*.svg);;All files (*)",
    "png": "PNG (*.png);;All files (*)",
    "pdf": "PDF (*.pdf);;All files (*)",
    "gff3": "GFF3 (*.gff3 *.gff);;All files (*)",
    "gtf": "GTF (*.gtf);;All files (*)",
    "bed": "BED (*.bed);;All files (*)",
}

IMAGE_FORMATS = ("pdf", "png", "svg")
ANNOTATION_FORMATS = ("gff3", "gtf", "bed")
# Only bitmaps have a resolution; PDF and SVG are vector.
DPI_FORMATS = frozenset({"png"})

AnnotationProvider = Callable[[], "tuple[list[str], dict, dict]"]

_TYPE_HELP = (
    "<b>Annotation file</b> (GFF3, GTF, BED): the features of this file, for "
    "IGV, genome browsers or other tools. Unsaved edits shown in the view are "
    "included.<br><br>"
    "<b>Image</b> (PDF, PNG, SVG): a picture of the current view. PDF and SVG "
    "are vector graphics (sharp at any zoom, editable in Illustrator or "
    "Inkscape); PNG is a bitmap, so it uses the DPI setting."
)
_AUTO_GFF3_HELP = (
    "Annotations are always saved in the JSON file next to the FASTA "
    "(<i>name</i>.seqmap_RV.json). Tick this to also rewrite "
    "<i>name</i>.fasta_RV.gff3 every time annotations are saved, so the GFF3 "
    "stays up to date for other tools. Applies to all files."
)

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


def write_view_vector_pdf(
    out_path: Path,
    paint: Callable,
    *,
    width: int,
    height: int,
    logical_dpi: float = 72.0,
) -> None:
    """Paint a view straight onto a PDF page, so lines and text stay vectors.

    ``paint(painter)`` draws in screen pixels (0..width, 0..height). The PDF
    resolution is set to the screen's logical DPI so fonts keep the same size
    relative to the drawing as on screen; the page is sized to fit exactly.
    """
    from PySide6.QtCore import QMarginsF, QRect, QSizeF
    from PySide6.QtGui import QColor, QPageLayout, QPageSize, QPainter, QPdfWriter

    dpi = float(logical_dpi) if logical_dpi and logical_dpi > 0 else 72.0
    width = max(1, int(width))
    height = max(1, int(height))
    writer = QPdfWriter(str(out_path))
    writer.setResolution(int(round(dpi)))
    writer.setPageLayout(
        QPageLayout(
            QPageSize(
                QSizeF(width * 72.0 / dpi, height * 72.0 / dpi),
                QPageSize.Unit.Point,
            ),
            QPageLayout.Orientation.Portrait,
            QMarginsF(0, 0, 0, 0),
        )
    )
    painter = QPainter(writer)
    try:
        painter.fillRect(QRect(0, 0, width, height), QColor(255, 255, 255))
        paint(painter)
    finally:
        painter.end()


@contextmanager
def vector_export(root) -> Iterator[None]:
    """Make cached canvases under ``root`` paint live (vectors) while exporting.

    Canvases that blit a pixmap cache on screen expose a ``_vector_export``
    attribute; with it set they draw their content directly instead, so a
    PDF / SVG gets real shapes and text rather than an embedded bitmap.
    """
    from PySide6.QtWidgets import QWidget as _QWidget

    widgets = [root, *root.findChildren(_QWidget)]
    flagged = [w for w in widgets if hasattr(w, "_vector_export")]
    for w in flagged:
        w._vector_export = True
    try:
        yield
    finally:
        for w in flagged:
            try:
                w._vector_export = False
            except RuntimeError:
                pass


class ExportAnnotationViewDialog(QDialog):
    """Type + DPI + File. Image types always; annotation types when available."""

    def __init__(
        self,
        *,
        source_path: Path,
        title: str = "Export annotation view",
        name_tag: str = "overview",
        default_format: str = "pdf",
        annotation_provider: AnnotationProvider | None = None,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumWidth(480)
        self._source_path = source_path.expanduser()
        self._name_tag = str(name_tag or "overview").strip() or "overview"
        self._annotation_provider = annotation_provider
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
        if annotation_provider is not None:
            self._type.addItem("GFF3", "gff3")
            self._type.addItem("GTF", "gtf")
            self._type.addItem("BED", "bed")
            self._type.insertSeparator(self._type.count())
            self._type.addItem("PDF", "pdf")
            self._type.addItem("PNG", "png")
            self._type.addItem("SVG", "svg")
        else:
            self._type.addItem("PDF", "pdf")
            self._type.addItem("PNG", "png")
            self._type.addItem("SVG", "svg")
        idx = self._type.findData(saved_fmt)
        if idx < 0:
            idx = self._type.findData(default_format)
        self._type.setCurrentIndex(max(0, idx))
        size_export_format_combo(self._type)
        self._last_fmt = self._current_fmt()

        self._dpi = QSpinBox(self)
        self._dpi.setRange(_MIN_DPI, _MAX_DPI)
        self._dpi.setSingleStep(25)
        self._dpi.setValue(saved_dpi)
        self._dpi.setSuffix(" DPI")
        self._dpi.setMinimumWidth(max(110, self._dpi.sizeHint().width()))

        self._auto_gff3 = QCheckBox("Also write GFF3 on every Save", self)
        self._auto_gff3.setChecked(_auto_gff3_enabled())

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

        from repeatviewer.widgets.help_button import field_with_help

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setSpacing(8)
        type_field = (
            field_with_help(self._type, "Export type", _TYPE_HELP, parent=self)
            if annotation_provider is not None
            else self._type
        )
        form.addRow("Type:", type_field)
        form.addRow("DPI:", self._dpi)
        self._auto_gff3_row = field_with_help(
            self._auto_gff3, "Also write GFF3 on every Save", _AUTO_GFF3_HELP,
            parent=self,
        )
        form.addRow("", self._auto_gff3_row)
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

    def wrote_annotation_file(self) -> bool:
        """True when OK wrote a GFF3 / GTF / BED (nothing left for the caller)."""
        return self._chosen is not None and self._current_fmt() in ANNOTATION_FORMATS

    def _current_fmt(self) -> str:
        fmt = str(self._type.currentData() or "pdf")
        return fmt if fmt in _FORMAT_SUFFIX else "pdf"

    def _file_filter(self) -> str:
        return _FORMAT_FILTER[self._current_fmt()]

    def _default_out_for(self, fmt: str) -> Path:
        """Image: ``{source.name}_{tag}.pdf``; annotation: ``{source.name}_RV.gff3``."""
        if fmt in ANNOTATION_FORMATS:
            from repeatviewer.backends.sequence_features import export_annotation_path

            return export_annotation_path(self._source_path, fmt)
        return (
            self._source_path.parent
            / f"{self._source_path.name}_{self._name_tag}{_FORMAT_SUFFIX[fmt]}"
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
        fmt = self._current_fmt()
        self._form.setRowVisible(self._dpi, fmt in DPI_FORMATS)
        self._form.setRowVisible(
            self._auto_gff3_row,
            self._annotation_provider is not None and fmt in ANNOTATION_FORMATS,
        )
        self.adjustSize()

    def _sync_dpi_enabled(self) -> None:
        # Kept for older call sites; DPI row visibility now follows the type.
        self._sync_rows()

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
        if was_default:
            # Untouched suggestion: switch to the new type's default name
            # (image ``_annotation.pdf`` <-> annotation ``_RV.gff3``).
            updated = self._default_out_for(fmt)
        else:
            updated = path.with_suffix(_FORMAT_SUFFIX[fmt])
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
        annotation = fmt in ANNOTATION_FORMATS and self._annotation_provider is not None
        if annotation:
            record_ids, features_by_id, lengths_by_id = self._annotation_provider()
            if not features_by_id or not any(features_by_id.values()):
                QMessageBox.information(
                    self, self.windowTitle(), "No features to export."
                )
                return
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
        if annotation:
            from repeatviewer.views.export_annotation_dialog import (
                write_annotation_export,
            )

            if not write_annotation_export(
                self,
                self._source_path.expanduser().resolve(),
                fmt,
                confirmed,
                record_ids=list(record_ids),
                features_by_id=features_by_id,
                lengths_by_id=lengths_by_id,
                title=self.windowTitle(),
            ):
                return
            _set_auto_gff3_enabled(self._auto_gff3.isChecked())
        self._chosen = confirmed
        self._dpi_value = int(self._dpi.value())
        settings = QSettings()
        settings.beginGroup(_SETTINGS_GROUP)
        settings.setValue("format", fmt)
        settings.setValue("dpi", self._dpi_value)
        settings.endGroup()
        self.accept()


def _auto_gff3_enabled() -> bool:
    from repeatviewer.backends.sequence_features import auto_gff3_on_save_enabled

    return auto_gff3_on_save_enabled()


def _set_auto_gff3_enabled(enabled: bool) -> None:
    from repeatviewer.backends.sequence_features import set_auto_gff3_on_save

    set_auto_gff3_on_save(enabled)
