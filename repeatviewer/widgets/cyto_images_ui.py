"""Cytogenetic images: attach dialog, viewer, picker and the Cytogenetics panel.

Storage lives in :mod:`repeatviewer.backends.cyto_images`; this module is the
Qt side. Every change is announced through :func:`cyto_notifier` so the file
list tick column, annotation markers and open panels refresh together.
"""

from __future__ import annotations

import logging
import platform
from pathlib import Path

from PySide6.QtCore import QObject, QRectF, QSettings, QSize, Qt, QTimer, QUrl, Signal
from PySide6.QtGui import (
    QColor,
    QDesktopServices,
    QFont,
    QFontMetrics,
    QIcon,
    QImage,
    QImageReader,
    QPainter,
    QPixmap,
    QTextOption,
)
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QColorDialog,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGraphicsPixmapItem,
    QGraphicsScene,
    QGraphicsView,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSplitter,
    QTableWidget,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.cyto_images import (
    COUNTERSTAIN_CHOICES,
    DYE_CHOICES,
    IMAGE_SUFFIXES,
    METADATA_FIELDS,
    METADATA_KEYS,
    STAGE_CHOICES,
    TECHNIQUE_CHOICES,
    CytoImage,
    FeatureRef,
    ImageTarget,
    add_image,
    build_caption,
    default_dye_color,
    image_file_date,
    image_file_path,
    images_for_feature,
    images_for_sequence,
    is_image_path,
    unlink_target,
    update_image,
)

logger = logging.getLogger(__name__)

_FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".faa", ".fna"}
_SETTINGS_GROUP = "cytogenetics"
# Fields whose last value is offered again on the next attach.
_REMEMBERED_FIELDS = ("species", "accession", "technique", "counterstain", "stage", "microscope")
IMAGE_FILE_FILTER = "Images ({})".format(" ".join(f"*{s}" for s in IMAGE_SUFFIXES))


# ---------------------------------------------------------------------------
# Change notification
# ---------------------------------------------------------------------------


class CytoNotifier(QObject):
    """App-wide signal: images of ``folder`` changed."""

    changed = Signal(str)


_NOTIFIER: CytoNotifier | None = None


def cyto_notifier() -> CytoNotifier:
    global _NOTIFIER
    if _NOTIFIER is None:
        _NOTIFIER = CytoNotifier()
    return _NOTIFIER


def notify_changed(folder: Path) -> None:
    cyto_notifier().changed.emit(str(Path(folder).expanduser().resolve()))


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def image_paths_from_mime(mime) -> list[Path]:
    """Local image files carried by a drag (Finder / Explorer)."""
    if mime is None or not mime.hasUrls():
        return []
    out: list[Path] = []
    for url in mime.urls():
        if not url.isLocalFile():
            continue
        path = Path(url.toLocalFile())
        if path.is_file() and is_image_path(path):
            out.append(path)
    return out


def load_image(path: Path) -> QImage:
    reader = QImageReader(str(path))
    reader.setAutoTransform(True)
    image = reader.read()
    return image


_THUMB_CACHE: dict[tuple[str, float, int], QPixmap] = {}


def thumbnail(path: Path, side: int) -> QPixmap:
    try:
        mtime = path.stat().st_mtime
    except OSError:
        return QPixmap()
    key = (str(path), mtime, side)
    cached = _THUMB_CACHE.get(key)
    if cached is not None:
        return cached
    reader = QImageReader(str(path))
    reader.setAutoTransform(True)
    size = reader.size()
    if size.isValid() and (size.width() > side or size.height() > side):
        reader.setScaledSize(size.scaled(side, side, Qt.AspectRatioMode.KeepAspectRatio))
    pix = QPixmap.fromImage(reader.read())
    if len(_THUMB_CACHE) > 400:
        _THUMB_CACHE.clear()
    _THUMB_CACHE[key] = pix
    return pix


def _last_values() -> dict[str, str]:
    s = QSettings()
    s.beginGroup(_SETTINGS_GROUP)
    out = {key: str(s.value(f"last_{key}", "") or "") for key in _REMEMBERED_FIELDS}
    s.endGroup()
    return out


def _save_last_values(metadata: dict[str, str]) -> None:
    s = QSettings()
    s.beginGroup(_SETTINGS_GROUP)
    for key in _REMEMBERED_FIELDS:
        value = str(metadata.get(key, "") or "").strip()
        if value:
            s.setValue(f"last_{key}", value)
    s.endGroup()


def _fasta_files_in(folder: Path) -> list[str]:
    try:
        return sorted(
            p.name
            for p in folder.iterdir()
            if p.is_file() and p.suffix.lower() in _FASTA_SUFFIXES and not p.name.startswith(".")
        )
    except OSError:
        return []


def _features_of(fasta: Path) -> list[FeatureRef]:
    try:
        from repeatviewer.backends.sequence_features import load_sequence_map_data

        features = load_sequence_map_data(fasta).features
    except Exception:  # noqa: BLE001 - multi-sequence or unreadable FASTA
        return []
    return [FeatureRef(f.name or f.feature_type, f.start, f.end) for f in features]


def _feature_text(ref: FeatureRef) -> str:
    return f"{ref.name}  ({ref.start:,}–{ref.end:,})"


def _swatch_icon(color: str, size: int = 12) -> QIcon:
    pix = QPixmap(size, size)
    pix.fill(Qt.GlobalColor.transparent)
    if color:
        painter = QPainter(pix)
        painter.setRenderHint(QPainter.RenderHint.Antialiasing)
        painter.setBrush(QColor(color))
        painter.setPen(QColor(0, 0, 0, 90))
        painter.drawEllipse(1, 1, size - 2, size - 2)
        painter.end()
    return QIcon(pix)


# ---------------------------------------------------------------------------
# Attach / edit dialog
# ---------------------------------------------------------------------------


class _ColorButton(QPushButton):
    def __init__(self, color: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedWidth(64)
        self.clicked.connect(self._pick)
        self.set_color(color)

    def color(self) -> str:
        return self._color

    def set_color(self, color: str) -> None:
        self._color = str(color or "")
        self.setIcon(_swatch_icon(self._color, 14))
        self.setText("" if self._color else "None")

    def _pick(self) -> None:
        start = QColor(self._color) if self._color else QColor("#ff3b30")
        chosen = QColorDialog.getColor(start, self, "Signal colour")
        if chosen.isValid():
            self.set_color(chosen.name())


class _ProbeTable(QTableWidget):
    """Probes shown in the image: sequence, annotation, dye, colour."""

    COL_SEQ, COL_FEATURE, COL_PROBE, COL_DYE, COL_COLOR = range(5)

    def __init__(self, folder: Path, parent: QWidget | None = None) -> None:
        super().__init__(0, 5, parent)
        self._folder = folder
        self._fasta_names = _fasta_files_in(folder)
        self.setHorizontalHeaderLabels(["Sequence", "Annotation", "Probe name", "Dye", "Colour"])
        self.verticalHeader().setVisible(False)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        header = self.horizontalHeader()
        header.setSectionResizeMode(self.COL_SEQ, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(self.COL_FEATURE, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(self.COL_PROBE, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(self.COL_DYE, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(self.COL_COLOR, QHeaderView.ResizeMode.ResizeToContents)
        self.setMinimumHeight(84)
        self.setMaximumHeight(150)

    def add_target(self, target: ImageTarget | None = None) -> None:
        row = self.rowCount()
        self.insertRow(row)
        seq = QComboBox(self)
        names = list(self._fasta_names)
        if target is not None and target.sequence and target.sequence not in names:
            names.insert(0, target.sequence)
        seq.addItems(names)
        if target is not None and target.sequence:
            seq.setCurrentText(target.sequence)
        feature = QComboBox(self)
        probe = QLineEdit(self)
        dye = QComboBox(self)
        dye.setEditable(True)
        dye.addItem("")
        for name, color in DYE_CHOICES:
            dye.addItem(_swatch_icon(color), name)
        color_btn = _ColorButton("", self)
        self.setCellWidget(row, self.COL_SEQ, seq)
        self.setCellWidget(row, self.COL_FEATURE, feature)
        self.setCellWidget(row, self.COL_PROBE, probe)
        self.setCellWidget(row, self.COL_DYE, dye)
        self.setCellWidget(row, self.COL_COLOR, color_btn)

        def _refill_features(_index: int = 0, *, keep: FeatureRef | None = None) -> None:
            feature.clear()
            feature.addItem("(whole sequence)", None)
            refs = _features_of(self._folder / seq.currentText()) if seq.currentText() else []
            if keep is not None and not any(
                r.name == keep.name and r.start == keep.start and r.end == keep.end for r in refs
            ):
                refs.insert(0, keep)
            for ref in refs:
                feature.addItem(_feature_text(ref), ref)
            if keep is not None:
                for i in range(feature.count()):
                    ref = feature.itemData(i)
                    if ref is not None and (ref.name, ref.start, ref.end) == (keep.name, keep.start, keep.end):
                        feature.setCurrentIndex(i)
                        break
            _sync_probe_placeholder()

        def _sync_probe_placeholder(_index: int = 0) -> None:
            ref = feature.currentData()
            if ref is not None and ref.name:
                probe.setPlaceholderText(ref.name)
            else:
                probe.setPlaceholderText(Path(seq.currentText()).stem)

        def _on_dye(text: str) -> None:
            color = default_dye_color(text)
            if color:
                color_btn.set_color(color)

        seq.currentIndexChanged.connect(_refill_features)
        feature.currentIndexChanged.connect(_sync_probe_placeholder)
        dye.currentTextChanged.connect(_on_dye)
        _refill_features(keep=target.feature if target is not None else None)
        if target is not None:
            probe.setText(target.probe)
            dye.setCurrentText(target.dye)
            if target.color:
                color_btn.set_color(target.color)

    def targets(self) -> list[ImageTarget]:
        out: list[ImageTarget] = []
        for row in range(self.rowCount()):
            seq = self.cellWidget(row, self.COL_SEQ)
            feature = self.cellWidget(row, self.COL_FEATURE)
            probe = self.cellWidget(row, self.COL_PROBE)
            dye = self.cellWidget(row, self.COL_DYE)
            color = self.cellWidget(row, self.COL_COLOR)
            name = seq.currentText().strip() if seq is not None else ""
            if not name:
                continue
            out.append(
                ImageTarget(
                    sequence=name,
                    feature=feature.currentData() if feature is not None else None,
                    probe=probe.text() if probe is not None else "",
                    dye=dye.currentText() if dye is not None else "",
                    color=color.color() if color is not None else "",
                )
            )
        return out

    def remove_selected(self) -> None:
        rows = sorted({i.row() for i in self.selectedIndexes()}, reverse=True)
        if not rows and self.rowCount() > 1:
            rows = [self.rowCount() - 1]
        for row in rows:
            if self.rowCount() > 1:
                self.removeRow(row)


class AttachImageDialog(QDialog):
    """Preview on the left, optional details on the right. Save stores as is."""

    def __init__(
        self,
        folder: Path,
        *,
        source: Path | None = None,
        image: CytoImage | None = None,
        targets: list[ImageTarget] | None = None,
        allow_no_targets: bool = False,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self._allow_no_targets = bool(allow_no_targets)
        self._folder = Path(folder).expanduser().resolve()
        self._source = source
        self._image = image
        self._result: CytoImage | None = None
        editing = image is not None
        self.setWindowTitle("Image details" if editing else "Add cytogenetic image")
        self.setMinimumSize(880, 600)

        # Left: the picture, at its own proportions, with its file name
        preview_path = image_file_path(self._folder, image) if editing else source
        preview = QLabel(self)
        preview.setAlignment(Qt.AlignmentFlag.AlignCenter)
        preview.setFixedWidth(320)
        preview.setMinimumHeight(320)
        preview.setObjectName("CytoImagePreview")
        preview.setStyleSheet("#CytoImagePreview { background: #111214; border-radius: 6px; color: #bbbbbb; }")
        pix = thumbnail(Path(preview_path), 640) if preview_path else QPixmap()
        if pix.isNull():
            preview.setText("No preview")
        else:
            preview.setPixmap(pix.scaled(320, 420, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation))
            preview.setFixedHeight(min(420, max(200, int(320 * pix.height() / max(1, pix.width())))))
        file_name = (image.original_name if editing else Path(source).name) if (editing or source) else ""
        name = QLabel(self)
        name.setAlignment(Qt.AlignmentFlag.AlignCenter)
        name.setStyleSheet("color: gray;")
        name.setText(name.fontMetrics().elidedText(file_name, Qt.TextElideMode.ElideMiddle, 310))
        name.setToolTip(file_name)
        left_v = QVBoxLayout()
        left_v.setSpacing(6)
        left_v.addWidget(preview)
        left_v.addWidget(name)
        left_v.addStretch(1)

        # Right: probes, then the details in two groups; every field the same width
        self._probes = _ProbeTable(self._folder, self)
        for target in (image.targets if editing else targets or []):
            self._probes.add_target(target)
        if self._probes.rowCount() == 0:
            self._probes.add_target(None)
        self._probes.verticalHeader().setVisible(False)
        self._probes.verticalHeader().setDefaultSectionSize(30)
        self._probes.setMinimumHeight(30 * 3 + 30)
        self._probes.setMaximumHeight(30 * 5 + 30)
        add_btn = QPushButton("+", self)
        add_btn.setFixedWidth(32)
        add_btn.setToolTip("Add a probe (a multicolour image shows several)")
        add_btn.clicked.connect(lambda: self._probes.add_target(None))
        remove_btn = QPushButton("−", self)
        remove_btn.setFixedWidth(32)
        remove_btn.setToolTip("Remove the selected probe")
        remove_btn.clicked.connect(self._probes.remove_selected)
        probe_buttons = QHBoxLayout()
        probe_buttons.setSpacing(4)
        probe_buttons.addWidget(add_btn)
        probe_buttons.addWidget(remove_btn)
        probe_buttons.addStretch()

        existing = dict(image.metadata) if editing else {}
        defaults = {} if editing else self._defaults()
        labels = dict(METADATA_FIELDS)
        self._fields: dict[str, QWidget] = {}

        def form() -> QFormLayout:
            f = QFormLayout()
            f.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            f.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.AllNonFixedFieldsGrow)
            f.setRowWrapPolicy(QFormLayout.RowWrapPolicy.DontWrapRows)
            f.setHorizontalSpacing(12)
            f.setVerticalSpacing(8)
            return f

        def field(key: str) -> QWidget:
            value = existing.get(key, defaults.get(key, ""))
            if key in {"technique", "counterstain", "stage"}:
                widget = QComboBox(self)
                widget.setEditable(True)
                choices = {"technique": TECHNIQUE_CHOICES, "counterstain": COUNTERSTAIN_CHOICES, "stage": STAGE_CHOICES}[key]
                widget.addItem("")
                widget.addItems(list(choices))
                widget.setCurrentText(value)
            elif key == "notes":
                widget = QPlainTextEdit(self)
                widget.setPlainText(value)
                widget.setFixedHeight(64)
            else:
                widget = QLineEdit(self)
                widget.setText(value)
                if key == "date":
                    widget.setPlaceholderText("YYYY-MM-DD")
                elif key == "chromosomes":
                    widget.setPlaceholderText("for example 1, 3, 5 (subtelomeric)")
            widget.setMinimumWidth(260)
            self._fields[key] = widget
            return widget

        sample = form()
        for key in ("species", "accession", "stage", "chromosomes"):
            sample.addRow(labels[key], field(key))
        method = form()
        for key in ("technique", "counterstain", "microscope", "experimenter", "date"):
            method.addRow(labels[key], field(key))
        notes = form()
        notes.addRow(labels["notes"], field("notes"))
        # make the three label columns the same width so the fields line up
        width = max(self.fontMetrics().horizontalAdvance(labels[k]) for k in METADATA_KEYS) + 4
        for f in (sample, method, notes):
            for row in range(f.rowCount()):
                item = f.itemAt(row, QFormLayout.ItemRole.LabelRole)
                if item is not None and item.widget() is not None:
                    item.widget().setFixedWidth(width)
                    item.widget().setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)

        def heading(text: str) -> QLabel:
            label = QLabel(text, self)
            font = label.font()
            font.setBold(True)
            label.setFont(font)
            return label

        right_v = QVBoxLayout()
        right_v.setSpacing(8)
        right_v.addWidget(heading("Probes in this image"))
        right_v.addWidget(self._probes)
        right_v.addLayout(probe_buttons)
        right_v.addSpacing(10)
        right_v.addWidget(heading("Sample"))
        right_v.addLayout(sample)
        right_v.addSpacing(10)
        right_v.addWidget(heading("Method"))
        right_v.addLayout(method)
        right_v.addSpacing(10)
        right_v.addLayout(notes)
        right_v.addStretch(1)
        right = QWidget(self)
        right.setLayout(right_v)
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setWidget(right)

        body = QHBoxLayout()
        body.setSpacing(24)
        body.addLayout(left_v)
        body.addWidget(scroll, 1)

        hint = QLabel("Only the image is needed; everything else can be filled in later.", self)
        hint.setStyleSheet("color: gray;")
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Cancel | QDialogButtonBox.StandardButton.Save,
            parent=self,
        )
        save_btn = buttons.button(QDialogButtonBox.StandardButton.Save)
        save_btn.setDefault(True)
        save_btn.setAutoDefault(True)
        buttons.accepted.connect(self._save)
        buttons.rejected.connect(self.reject)
        bottom = QHBoxLayout()
        bottom.addWidget(hint, 1)
        bottom.addWidget(buttons)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 20, 20, 16)
        layout.setSpacing(16)
        layout.addLayout(body, 1)
        layout.addLayout(bottom)
        save_btn.setFocus()

    def _defaults(self) -> dict[str, str]:
        values = _last_values()
        if self._source is not None:
            values["date"] = image_file_date(self._source)
        try:
            from repeatviewer.app_preferences import load_cyto_experimenter

            values["experimenter"] = load_cyto_experimenter()
        except Exception:  # noqa: BLE001
            pass
        return values

    def metadata(self) -> dict[str, str]:
        out: dict[str, str] = {}
        for key, widget in self._fields.items():
            if isinstance(widget, QComboBox):
                text = widget.currentText()
            elif isinstance(widget, QPlainTextEdit):
                text = widget.toPlainText()
            else:
                text = widget.text()
            if text.strip():
                out[key] = text.strip()
        return out

    def result_image(self) -> CytoImage | None:
        return self._result

    def _save(self) -> None:
        targets = self._probes.targets()
        if not targets and not self._allow_no_targets:
            QMessageBox.warning(self, "Attach image", "Pick at least one sequence for this image.")
            return
        metadata = self.metadata()
        try:
            if self._image is None:
                self._result = add_image(self._folder, Path(self._source), targets=targets, metadata=metadata)
            else:
                self._image.targets = targets
                self._image.metadata = metadata
                update_image(self._folder, self._image)
                self._result = self._image
        except OSError as exc:
            QMessageBox.critical(self, "Attach image", f"Could not save the image:\n{exc}")
            return
        _save_last_values(metadata)
        notify_changed(self._folder)
        self.accept()


def attach_images(
    parent: QWidget | None,
    fasta_path: Path,
    images: list[Path] | None = None,
    *,
    feature: FeatureRef | None = None,
) -> int:
    """Attach flow shared by drag-and-drop and right-click. Returns # saved."""
    fasta_path = Path(fasta_path).expanduser().resolve()
    if not images:
        chosen, _ = QFileDialog.getOpenFileNames(
            parent,
            "Attach cytogenetic image",
            str(fasta_path.parent),
            IMAGE_FILE_FILTER,
        )
        images = [Path(p) for p in chosen]
    saved = 0
    for source in images:
        if not is_image_path(source):
            continue
        dialog = AttachImageDialog(
            fasta_path.parent,
            source=source,
            targets=[ImageTarget(fasta_path.name, feature=feature)],
            parent=parent,
        )
        if dialog.exec() == QDialog.DialogCode.Accepted:
            saved += 1
    return saved


def edit_image(parent: QWidget | None, folder: Path, image: CytoImage, *, allow_no_targets: bool = False) -> bool:
    dialog = AttachImageDialog(folder, image=image, allow_no_targets=allow_no_targets, parent=parent)
    return dialog.exec() == QDialog.DialogCode.Accepted


# ---------------------------------------------------------------------------
# Copy / export
# ---------------------------------------------------------------------------


def copy_image_to_clipboard(folder: Path, image: CytoImage) -> bool:
    qimage = load_image(image_file_path(folder, image))
    if qimage.isNull():
        return False
    QApplication.clipboard().setImage(qimage)
    return True


def render_with_caption(qimage: QImage, caption: str) -> QImage:
    """Image with a white caption strip underneath."""
    width = max(qimage.width(), 400)
    font = QFont()
    font.setPixelSize(max(14, width // 60))
    fm = QFontMetrics(font)
    margin = max(10, font.pixelSize())
    text_rect = fm.boundingRect(
        0, 0, width - 2 * margin, 10_000, int(Qt.TextFlag.TextWordWrap), caption
    )
    strip = text_rect.height() + 2 * margin
    out = QImage(width, qimage.height() + strip, QImage.Format.Format_RGB32)
    out.fill(QColor("white"))
    painter = QPainter(out)
    painter.setRenderHint(QPainter.RenderHint.Antialiasing)
    painter.drawImage((width - qimage.width()) // 2, 0, qimage)
    painter.setFont(font)
    painter.setPen(QColor("#111111"))
    option = QTextOption()
    option.setWrapMode(QTextOption.WrapMode.WordWrap)
    painter.drawText(
        QRectF(margin, qimage.height() + margin, width - 2 * margin, strip - margin),
        caption,
        option,
    )
    painter.end()
    return out


def export_with_caption(parent: QWidget | None, folder: Path, image: CytoImage) -> Path | None:
    source = image_file_path(folder, image)
    qimage = load_image(source)
    if qimage.isNull():
        QMessageBox.warning(parent, "Export image", "The image could not be read.")
        return None
    suggested = Path(folder) / f"{Path(image.title()).stem.replace(' + ', '_')}_captioned.png"
    target, _ = QFileDialog.getSaveFileName(
        parent,
        "Export image with caption",
        str(suggested),
        "PNG (*.png);;JPEG (*.jpg);;TIFF (*.tif)",
    )
    if not target:
        return None
    out = render_with_caption(qimage, build_caption(image))
    if not out.save(target):
        QMessageBox.critical(parent, "Export image", f"Could not write:\n{target}")
        return None
    return Path(target)


def reveal_image_file(folder: Path, image: CytoImage) -> None:
    path = image_file_path(folder, image)
    if platform.system() == "Darwin":
        import subprocess

        subprocess.Popen(["open", "-R", str(path)])  # noqa: S603,S607
    else:
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(path.parent)))


# ---------------------------------------------------------------------------
# Viewer
# ---------------------------------------------------------------------------


class _ZoomView(QGraphicsView):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setScene(QGraphicsScene(self))
        self._item = QGraphicsPixmapItem()
        self._item.setTransformationMode(Qt.TransformationMode.SmoothTransformation)
        self.scene().addItem(self._item)
        self.setDragMode(QGraphicsView.DragMode.ScrollHandDrag)
        self.setTransformationAnchor(QGraphicsView.ViewportAnchor.AnchorUnderMouse)
        self.setBackgroundBrush(QColor("#111111"))
        self.setFrameShape(QFrame.Shape.NoFrame)
        self._fit = True

    def set_image(self, qimage: QImage) -> None:
        self._item.setPixmap(QPixmap.fromImage(qimage))
        self.scene().setSceneRect(QRectF(self._item.pixmap().rect()))
        self.fit()

    def fit(self) -> None:
        self._fit = True
        self.resetTransform()
        if not self._item.pixmap().isNull():
            self.fitInView(self._item, Qt.AspectRatioMode.KeepAspectRatio)

    def actual_size(self) -> None:
        self._fit = False
        self.resetTransform()

    def zoom(self, factor: float) -> None:
        self._fit = False
        self.scale(factor, factor)

    def wheelEvent(self, event) -> None:  # noqa: N802
        steps = event.angleDelta().y() / 120.0
        if steps:
            self.zoom(1.15**steps)
            event.accept()

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        if self._fit:
            self.fit()


class CytoImageViewer(QDialog):
    """Zoomable image with a probe/dye legend and the generated caption."""

    showSourceRequested = Signal(str, object)  # sequence name, FeatureRef | None

    def __init__(self, folder: Path, images: list[CytoImage], index: int = 0, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self.setModal(False)
        self.resize(900, 720)
        self._folder = Path(folder)
        self._images = list(images)
        self._index = max(0, min(index, len(self._images) - 1))

        self._view = _ZoomView(self)
        self._legend = QLabel(self)
        self._legend.setTextFormat(Qt.TextFormat.RichText)
        self._legend.setWordWrap(True)
        self._caption = QLabel(self)
        self._caption.setWordWrap(True)
        self._caption.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        caption_font = self._caption.font()
        caption_font.setItalic(True)
        self._caption.setFont(caption_font)

        fit = QPushButton("Fit", self)
        fit.clicked.connect(self._view.fit)
        one = QPushButton("100%", self)
        one.clicked.connect(self._view.actual_size)
        zoom_in = QPushButton("+", self)
        zoom_in.setFixedWidth(32)
        zoom_in.clicked.connect(lambda: self._view.zoom(1.25))
        zoom_out = QPushButton("−", self)
        zoom_out.setFixedWidth(32)
        zoom_out.clicked.connect(lambda: self._view.zoom(0.8))
        self._prev = QPushButton("◀", self)
        self._prev.setFixedWidth(32)
        self._prev.clicked.connect(lambda: self._show(self._index - 1))
        self._next = QPushButton("▶", self)
        self._next.setFixedWidth(32)
        self._next.clicked.connect(lambda: self._show(self._index + 1))
        self._counter = QLabel(self)
        copy_btn = QPushButton("Copy", self)
        copy_btn.setToolTip("Copy the image to the clipboard")
        copy_btn.clicked.connect(self._copy)
        export_btn = QPushButton("Export with caption…", self)
        export_btn.clicked.connect(lambda: export_with_caption(self, self._folder, self._current()))

        bar = QHBoxLayout()
        for w in (self._prev, self._counter, self._next):
            bar.addWidget(w)
        bar.addSpacing(12)
        for w in (zoom_out, zoom_in, fit, one):
            bar.addWidget(w)
        bar.addStretch()
        bar.addWidget(copy_btn)
        bar.addWidget(export_btn)

        layout = QVBoxLayout(self)
        layout.addWidget(self._legend)
        layout.addWidget(self._view, 1)
        layout.addWidget(self._caption)
        layout.addLayout(bar)
        self._show(self._index)

    def _current(self) -> CytoImage:
        return self._images[self._index]

    def _show(self, index: int) -> None:
        if not self._images:
            return
        self._index = index % len(self._images)
        image = self._current()
        qimage = load_image(image_file_path(self._folder, image))
        self._view.set_image(qimage)
        self.setWindowTitle(f"{image.title()} — cytogenetic image")
        chips = []
        for target in image.targets:
            dot = (
                f"<span style='color:{target.color}; font-size:14pt;'>●</span> "
                if target.color
                else ""
            )
            dye = f" <span style='color:gray;'>({target.dye})</span>" if target.dye else ""
            chips.append(f"{dot}<b>{target.probe_label}</b>{dye}")
        counter = image.metadata.get("counterstain", "")
        if counter:
            chips.append(f"<span style='color:#4a7dff; font-size:14pt;'>●</span> {counter}")
        self._legend.setText("&nbsp;&nbsp;&nbsp;".join(chips))
        self._caption.setText(build_caption(image))
        many = len(self._images) > 1
        for w in (self._prev, self._next, self._counter):
            w.setVisible(many)
        self._counter.setText(f"{self._index + 1} / {len(self._images)}")

    def _copy(self) -> None:
        if copy_image_to_clipboard(self._folder, self._current()):
            self._caption.setText("Image copied to the clipboard.   " + build_caption(self._current()))


class CytoImagePicker(QDialog):
    """Thumbnails of several images; click one to open it."""

    def __init__(self, folder: Path, images: list[CytoImage], title: str, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(640, 420)
        self._folder = Path(folder)
        self._images = list(images)
        self._list = QListWidget(self)
        self._list.setViewMode(QListWidget.ViewMode.IconMode)
        self._list.setIconSize(QSize(160, 160))
        self._list.setResizeMode(QListWidget.ResizeMode.Adjust)
        self._list.setMovement(QListWidget.Movement.Static)
        self._list.setSpacing(8)
        self._list.setWordWrap(True)
        for i, image in enumerate(self._images):
            item = QListWidgetItem(QIcon(thumbnail(image_file_path(self._folder, image), 160)), image.title())
            item.setData(Qt.ItemDataRole.UserRole, i)
            item.setToolTip(build_caption(image))
            self._list.addItem(item)
        self._list.itemActivated.connect(self._open)
        self._list.itemClicked.connect(self._open)
        layout = QVBoxLayout(self)
        layout.addWidget(QLabel("Click an image to open it.", self))
        layout.addWidget(self._list, 1)

    def _open(self, item: QListWidgetItem) -> None:
        index = int(item.data(Qt.ItemDataRole.UserRole))
        viewer = CytoImageViewer(self._folder, self._images, index, self.parent())
        viewer.show()
        self.accept()


def open_images(parent: QWidget | None, folder: Path, images: list[CytoImage], *, title: str = "") -> None:
    """One image → viewer; several → picker."""
    if not images:
        return
    if len(images) == 1:
        CytoImageViewer(folder, images, 0, parent).show()
        return
    picker = CytoImagePicker(folder, images, title or f"{len(images)} cytogenetic images", parent)
    picker.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
    picker.show()


def _in_workspace(parent: QWidget | None, fasta_path: Path, images: list[CytoImage]) -> bool:
    """Show the images in the Cytogenetics workspace (Library filtered to the sequence; one image opens directly)."""
    window = parent.window() if parent is not None else None
    opener = getattr(window, "_open_cytogenetics_workspace", None)
    if opener is None or not images:
        return False
    opener(tab="Library", sequence=fasta_path, image_id=images[0].id if len(images) == 1 else "")
    return True


def open_images_for_sequence(parent: QWidget | None, fasta_path: Path) -> None:
    fasta_path = Path(fasta_path).expanduser().resolve()
    images = images_for_sequence(fasta_path)
    if _in_workspace(parent, fasta_path, images):
        return
    open_images(parent, fasta_path.parent, images, title=fasta_path.name)


def open_images_for_feature(parent: QWidget | None, fasta_path: Path, name: str, start: int, end: int) -> None:
    fasta_path = Path(fasta_path).expanduser().resolve()
    if _in_workspace(parent, fasta_path, images_for_feature(fasta_path, name, start, end)):
        return
    open_images(
        parent,
        fasta_path.parent,
        images_for_feature(fasta_path, name, start, end),
        title=f"{name} — {fasta_path.name}",
    )


# ---------------------------------------------------------------------------
# Cytogenetics panel (FASTA preview sub-tab)
# ---------------------------------------------------------------------------


class CytogeneticsPanel(QWidget):
    """Thumbnails on top; the selected image and its details below."""

    showSourceRequested = Signal(object)  # FeatureRef | None

    def __init__(self, fasta_path: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setAcceptDrops(True)
        self._path = Path(fasta_path).expanduser().resolve()
        self._images: list[CytoImage] = []

        self._thumbs = QListWidget(self)
        self._thumbs.setViewMode(QListWidget.ViewMode.IconMode)
        self._thumbs.setFlow(QListWidget.Flow.LeftToRight)
        self._thumbs.setWrapping(False)
        self._thumbs.setIconSize(QSize(88, 88))
        self._thumbs.setFixedHeight(128)
        self._thumbs.setMovement(QListWidget.Movement.Static)
        self._thumbs.setSpacing(6)
        self._thumbs.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._thumbs.currentRowChanged.connect(self._on_current_changed)
        self._thumbs.itemDoubleClicked.connect(lambda _item: self._open())

        attach_btn = QPushButton("Attach image…", self)
        attach_btn.clicked.connect(lambda: attach_images(self, self._path))
        top = QHBoxLayout()
        top.addWidget(self._thumbs, 1)
        top.addWidget(attach_btn, 0, Qt.AlignmentFlag.AlignTop)

        self._big = QLabel(self)
        self._big.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._big.setMinimumSize(200, 200)
        self._big.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Ignored)
        self._big.setObjectName("CytoBigImage")
        self._big.setStyleSheet("#CytoBigImage { background: #1e1e1e; border-radius: 4px; }")
        self._big_pixmap = QPixmap()

        self._details = QLabel(self)
        self._details.setTextFormat(Qt.TextFormat.RichText)
        self._details.setWordWrap(True)
        self._details.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
        self._details.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        details_scroll = QScrollArea(self)
        details_scroll.setWidgetResizable(True)
        details_scroll.setFrameShape(QFrame.Shape.NoFrame)
        details_scroll.setWidget(self._details)
        details_scroll.setMinimumWidth(240)

        splitter = QSplitter(Qt.Orientation.Horizontal, self)
        splitter.addWidget(self._big)
        splitter.addWidget(details_scroll)
        splitter.setStretchFactor(0, 3)
        splitter.setStretchFactor(1, 2)
        splitter.setSizes([600, 320])

        self._buttons: list[QPushButton] = []

        def _btn(label: str, tip: str, slot) -> QPushButton:
            b = QPushButton(label, self)
            b.setToolTip(tip)
            b.clicked.connect(slot)
            self._buttons.append(b)
            return b

        bar = QHBoxLayout()
        bar.addWidget(_btn("Open", "Zoom in, with a dye / probe legend", self._open))
        bar.addWidget(_btn("Copy", "Copy the image to the clipboard", self._copy))
        bar.addWidget(_btn("Export with caption…", "Save the image with a caption built from its details", self._export))
        bar.addWidget(_btn("Show source", "Jump to the probe in the Annotation view", self._show_source))
        bar.addStretch()
        bar.addWidget(_btn("Edit…", "Change probes or details", self._edit))
        bar.addWidget(_btn("Remove", "Remove this image from this sequence", self._remove))

        self._content = QWidget(self)
        content_v = QVBoxLayout(self._content)
        content_v.setContentsMargins(0, 0, 0, 0)
        content_v.addWidget(splitter, 1)
        content_v.addLayout(bar)

        self._empty = QLabel(
            "No cytogenetic images for this sequence yet.\n\n"
            "Drop an image here, onto the file in Files, or onto an annotation —\n"
            "or click “Attach image…”.",
            self,
        )
        self._empty.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._empty.setStyleSheet("color: palette(mid);")

        layout = QVBoxLayout(self)
        layout.setContentsMargins(8, 8, 8, 8)
        layout.addLayout(top)
        layout.addWidget(self._content, 1)
        layout.addWidget(self._empty, 1)

        cyto_notifier().changed.connect(self._on_store_changed)
        self.reload()

    # -- data -------------------------------------------------------------

    def set_path(self, fasta_path: Path) -> None:
        self._path = Path(fasta_path).expanduser().resolve()
        self.reload()

    def path(self) -> Path:
        return self._path

    def _on_store_changed(self, folder: str) -> None:
        if folder == str(self._path.parent):
            self.reload()

    def reload(self) -> None:
        current_id = self._current().id if self._current() is not None else None
        self._images = images_for_sequence(self._path)
        self._thumbs.blockSignals(True)
        self._thumbs.clear()
        select_row = 0
        for row, image in enumerate(self._images):
            item = QListWidgetItem(QIcon(thumbnail(image_file_path(self._path.parent, image), 88)), image.title())
            item.setToolTip(build_caption(image))
            self._thumbs.addItem(item)
            if image.id == current_id:
                select_row = row
        self._thumbs.blockSignals(False)
        has = bool(self._images)
        self._content.setVisible(has)
        self._thumbs.setVisible(has)
        self._empty.setVisible(not has)
        if has:
            self._thumbs.setCurrentRow(select_row)
            self._on_current_changed(select_row)

    def _current(self) -> CytoImage | None:
        row = self._thumbs.currentRow() if hasattr(self, "_thumbs") else -1
        if 0 <= row < len(getattr(self, "_images", [])):
            return self._images[row]
        return None

    def _on_current_changed(self, row: int) -> None:
        image = self._current()
        if image is None:
            self._big_pixmap = QPixmap()
            self._big.clear()
            self._details.clear()
            return
        self._big_pixmap = thumbnail(image_file_path(self._path.parent, image), 1400)
        self._rescale_big()
        self._details.setText(self._details_html(image))

    def _details_html(self, image: CytoImage) -> str:
        rows = []
        probes = []
        for target in image.targets:
            dot = f"<span style='color:{target.color};'>●</span> " if target.color else ""
            where = target.sequence
            if target.feature is not None:
                where += f", {target.feature.name} {target.feature.start:,}–{target.feature.end:,}"
            dye = f" — {target.dye}" if target.dye else ""
            probes.append(f"{dot}<b>{target.probe_label}</b>{dye}<br><span style='color:gray;'>{where}</span>")
        rows.append("<p>" + "<br>".join(probes) + "</p>")
        table = []
        for key, label in METADATA_FIELDS:
            value = image.metadata.get(key, "")
            if value:
                safe = value.replace("&", "&amp;").replace("<", "&lt;").replace("\n", "<br>")
                table.append(f"<tr><td style='color:gray; padding-right:8px;'>{label}</td><td>{safe}</td></tr>")
        if table:
            rows.append("<table>" + "".join(table) + "</table>")
        else:
            rows.append("<p style='color:gray;'>No details yet — use Edit… to add species, dye, stage …</p>")
        rows.append(f"<p style='color:gray;'><i>{build_caption(image)}</i></p>")
        rows.append(f"<p style='color:gray; font-size:small;'>File: {image.original_name}</p>")
        return "".join(rows)

    def _rescale_big(self) -> None:
        if self._big_pixmap.isNull():
            self._big.setText("Preview not available")
            return
        self._big.setPixmap(
            self._big_pixmap.scaled(
                self._big.size(),
                Qt.AspectRatioMode.KeepAspectRatio,
                Qt.TransformationMode.SmoothTransformation,
            )
        )

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._rescale_big()

    # -- actions ----------------------------------------------------------

    def _open(self) -> None:
        if self._current() is None:
            return
        viewer = CytoImageViewer(self._path.parent, self._images, self._thumbs.currentRow(), self)
        viewer.show()

    def _copy(self) -> None:
        image = self._current()
        if image is not None:
            copy_image_to_clipboard(self._path.parent, image)

    def _export(self) -> None:
        image = self._current()
        if image is not None:
            export_with_caption(self, self._path.parent, image)

    def _show_source(self) -> None:
        image = self._current()
        if image is None:
            return
        targets = image.targets_for_sequence(self._path.name)
        feature = next((t.feature for t in targets if t.feature is not None), None)
        self.showSourceRequested.emit(feature)

    def _edit(self) -> None:
        image = self._current()
        if image is not None:
            edit_image(self, self._path.parent, image)

    def _remove(self) -> None:
        image = self._current()
        if image is None:
            return
        others = [t.sequence for t in image.targets if t.sequence != self._path.name]
        if others:
            text = (
                f"Remove “{image.title()}” from {self._path.name}?\n\n"
                f"It stays attached to: {', '.join(sorted(set(others)))}."
            )
        else:
            text = f"Remove “{image.title()}”?\n\nThe stored copy of the image is deleted."
        answer = QMessageBox.question(self, "Remove image", text)
        if answer != QMessageBox.StandardButton.Yes:
            return
        unlink_target(self._path.parent, image.id, self._path.name)
        notify_changed(self._path.parent)

    # -- drag and drop ----------------------------------------------------

    def dragEnterEvent(self, event) -> None:  # noqa: N802
        if image_paths_from_mime(event.mimeData()):
            event.acceptProposedAction()
        else:
            event.ignore()

    def dragMoveEvent(self, event) -> None:  # noqa: N802
        self.dragEnterEvent(event)

    def dropEvent(self, event) -> None:  # noqa: N802
        paths = image_paths_from_mime(event.mimeData())
        if not paths:
            event.ignore()
            return
        event.acceptProposedAction()
        # Leave the drop handler before opening a modal dialog (macOS Finder waits otherwise).
        QTimer.singleShot(0, lambda: attach_images(self, self._path, paths))


__all__ = [
    "CytoImagePicker",
    "CytoImageViewer",
    "CytogeneticsPanel",
    "AttachImageDialog",
    "attach_images",
    "cyto_notifier",
    "edit_image",
    "image_paths_from_mime",
    "notify_changed",
    "open_images",
    "open_images_for_feature",
    "open_images_for_sequence",
]
