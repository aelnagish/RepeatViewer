"""Library tab of the Cytogenetics workspace: all images of a folder, grouped and filtered."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings, Qt, QThread, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSlider,
    QSplitter,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import library as lb
from repeatviewer.backends.cytogenetics.importer import ImportRequest, run_request
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.help_button import make_help_button

_SETTINGS = "cytogenetics_workspace"
THUMB_MIN, THUMB_MAX, THUMB_DEFAULT = 90, 240, 150
FILTER_W, DETAILS_W = 220, 290

HELP_TEXT = (
    "<b>Library</b><br>"
    "Every image of this folder, grouped by slide (or by probe, species or cell stage).<br>"
    "<b>Click</b> an image to see its details, <b>double-click</b> to open it.<br>"
    "<b>Filters</b> on the left narrow the list; the numbers show how many images each value would give.<br>"
    "<b>Import slide session</b> reads channel TIFFs from the microscope and sorts them into cells.<br>"
    "<b>Add image</b> stores one finished picture (merged RGB or a single file).<br>"
    "<b>Edit lists</b> changes the probes, species, stages, dyes and techniques offered."
)


def _read(key: str, default: str) -> str:
    try:
        return str(QSettings().value(f"{_SETTINGS}/{key}", default))
    except Exception:  # noqa: BLE001
        return default


def _write(key: str, value) -> None:
    try:
        QSettings().setValue(f"{_SETTINGS}/{key}", value)
    except Exception:  # noqa: BLE001
        pass


class _ImportWorker(QThread):
    """Runs one slide-session import off the GUI thread."""

    progressed = Signal(int, int)
    finishedWith = Signal(object)  # (slide, images, errors) or an Exception

    def __init__(self, request: ImportRequest, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._request = request
        self._stop = False

    def request_stop(self) -> None:
        self._stop = True

    def run(self) -> None:  # noqa: D102 - QThread entry point
        def progress(done: int, total: int):
            self.progressed.emit(done, total)
            return not self._stop

        try:
            self.finishedWith.emit(run_request(self._request, progress=progress))
        except Exception as exc:  # noqa: BLE001
            self.finishedWith.emit(exc)


class _Thumb(QFrame):
    clicked = Signal(str)
    doubleClicked = Signal(str)

    def __init__(self, item: lb.LibraryItem, folder: Path, side: int, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        from repeatviewer.widgets.cyto_images_ui import thumbnail

        self.image_id = item.image.id
        self.setFrameShape(QFrame.Shape.NoFrame)
        self.setCursor(Qt.CursorShape.PointingHandCursor)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(3)
        pic = QLabel(self)
        pic.setFixedSize(side, side)
        pic.setAlignment(Qt.AlignmentFlag.AlignCenter)
        pic.setStyleSheet("background:#000; color:#aaa; border-radius:4px;")
        pix = thumbnail(ci.image_file_path(folder, item.image), side)
        if pix.isNull():
            pic.setText("file missing")
        else:
            pic.setPixmap(pix.scaled(side, side, Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation))
        lay.addWidget(pic)
        name = QLabel(self)
        name.setFixedWidth(side)
        name.setText(name.fontMetrics().elidedText(item.title, Qt.TextElideMode.ElideMiddle, side))
        lay.addWidget(name)
        bits = []
        if item.channels:
            bits.append(f"{item.channels} ch")
        bits.append("marked" if item.marked else "not marked")
        sub = QLabel(" · ".join(bits), self)
        sub.setStyleSheet("color:palette(dark); font-size:11px;")
        lay.addWidget(sub)
        self.setToolTip("\n".join(f"{k}: {v}" for k, v in lb.details_rows(item)[:8]) or item.title)
        self.set_selected(False)

    def set_selected(self, on: bool) -> None:
        self.setStyleSheet(
            "_Thumb { background: palette(highlight); border-radius: 8px; }" if on
            else "_Thumb { background: transparent; border-radius: 8px; }"
        )

    def mousePressEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self.clicked.emit(self.image_id)
        super().mousePressEvent(event)

    def mouseDoubleClickEvent(self, event) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            self.doubleClicked.emit(self.image_id)
        super().mouseDoubleClickEvent(event)


class LibraryTab(QWidget):
    """Library: thumbnails grouped by slide, filters on the left, details on the right."""

    name = "Library"
    imageOpened = Signal(str)  # image id (double-click)

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._folder: Path | None = None
        self._items: list[lb.LibraryItem] = []
        self._flt = lb.LibraryFilter()
        self._selected = ""
        self._thumbs: dict[str, _Thumb] = {}
        self._worker: _ImportWorker | None = None
        self._cols_for: list[tuple[QGridLayout, list[_Thumb]]] = []

        # ---- top bar
        self._import_btn = QPushButton("Import slide session", self)
        self._add_btn = QPushButton("Add image", self)
        self._lists_btn = QPushButton("Edit lists", self)
        self._import_btn.clicked.connect(self.start_import)
        self._add_btn.clicked.connect(self.add_image)
        self._lists_btn.clicked.connect(self.edit_lists)
        self._count = QLabel(self)
        self._group = NoWheelComboBox(self)
        self._group.addItems(lb.GROUP_BYS)
        self._group.setCurrentText(_read("library_group", "Slide"))
        self._group.currentTextChanged.connect(self._on_group)
        self._size = QSlider(Qt.Orientation.Horizontal, self)
        self._size.setRange(THUMB_MIN, THUMB_MAX)
        try:
            self._size.setValue(int(_read("library_thumb", str(THUMB_DEFAULT))))
        except ValueError:
            self._size.setValue(THUMB_DEFAULT)
        self._size.setFixedWidth(110)
        self._size.sliderReleased.connect(self._render)
        self._size.valueChanged.connect(lambda v: (_write("library_thumb", int(v)), None) and None)
        top = QHBoxLayout()
        top.setContentsMargins(10, 8, 10, 4)
        for w in (self._import_btn, self._add_btn, self._lists_btn):
            top.addWidget(w)
        top.addStretch(1)
        top.addWidget(self._count)
        top.addWidget(QLabel("Group by", self))
        top.addWidget(self._group)
        top.addWidget(QLabel("Size", self))
        top.addWidget(self._size)
        top.addWidget(make_help_button(self, "Library", HELP_TEXT))

        # ---- progress row (import)
        self._progress_row = QWidget(self)
        pr = QHBoxLayout(self._progress_row)
        pr.setContentsMargins(10, 0, 10, 4)
        self._progress_label = QLabel(self._progress_row)
        self._progress = QProgressBar(self._progress_row)
        self._stop = QPushButton("Stop", self._progress_row)
        self._stop.clicked.connect(self._on_stop)
        pr.addWidget(self._progress_label)
        pr.addWidget(self._progress, 1)
        pr.addWidget(self._stop)
        self._progress_row.hide()

        # ---- sequence chip
        self._chip = QWidget(self)
        cr = QHBoxLayout(self._chip)
        cr.setContentsMargins(10, 0, 10, 4)
        self._chip_label = QLabel(self._chip)
        show_all = QPushButton("Show all images", self._chip)
        show_all.clicked.connect(lambda: self.set_sequence(""))
        cr.addWidget(self._chip_label)
        cr.addWidget(show_all)
        cr.addStretch(1)
        self._chip.hide()

        # ---- filters
        self._search = QLineEdit(self)
        self._search.setPlaceholderText("Search: probe, slide, note")
        self._search.setClearButtonEnabled(True)
        self._search.textChanged.connect(self._on_search)
        self._filter_body = QWidget(self)
        self._filter_lay = QVBoxLayout(self._filter_body)
        self._filter_lay.setContentsMargins(0, 0, 0, 0)
        filters = QWidget(self)
        fl = QVBoxLayout(filters)
        fl.setContentsMargins(10, 4, 6, 4)
        fl.addWidget(self._search)
        fl.addWidget(self._filter_body)
        fl.addStretch(1)
        filters_scroll = QScrollArea(self)
        filters_scroll.setWidgetResizable(True)
        filters_scroll.setFrameShape(QFrame.Shape.NoFrame)
        filters_scroll.setWidget(filters)
        filters_scroll.setFixedWidth(FILTER_W)

        # ---- grid
        self._grid_body = QWidget(self)
        self._grid_lay = QVBoxLayout(self._grid_body)
        self._grid_lay.setContentsMargins(8, 4, 8, 8)
        self._grid_lay.setSpacing(14)
        self._grid_scroll = QScrollArea(self)
        self._grid_scroll.setWidgetResizable(True)
        self._grid_scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._grid_scroll.setWidget(self._grid_body)
        self._relayout_timer = QTimer(self)
        self._relayout_timer.setSingleShot(True)
        self._relayout_timer.setInterval(80)
        self._relayout_timer.timeout.connect(self._relayout)

        # ---- details
        self._details_body = QWidget(self)
        self._details_lay = QVBoxLayout(self._details_body)
        self._details_lay.setContentsMargins(10, 4, 10, 8)
        self._details_lay.setSpacing(8)
        details_scroll = QScrollArea(self)
        details_scroll.setWidgetResizable(True)
        details_scroll.setFrameShape(QFrame.Shape.NoFrame)
        details_scroll.setWidget(self._details_body)
        details_scroll.setFixedWidth(DETAILS_W)

        split = QSplitter(Qt.Orientation.Horizontal, self)
        split.setChildrenCollapsible(False)
        split.addWidget(filters_scroll)
        split.addWidget(self._grid_scroll)
        split.addWidget(details_scroll)
        split.setStretchFactor(1, 1)

        root = QVBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addLayout(top)
        root.addWidget(self._progress_row)
        root.addWidget(self._chip)
        root.addWidget(split, 1)
        self._show_details(None)

    # ------------------------------------------------------------------ data

    def folder(self) -> Path | None:
        return self._folder

    def refresh(self, folder: Path | None) -> None:
        self._folder = Path(folder) if folder else None
        self._items = lb.build_items(self._folder) if self._folder else []
        known = {i.image.id for i in self._items}
        if self._selected not in known:
            self._selected = ""
        self._apply()

    def set_sequence(self, fasta_name: str) -> None:
        self._flt.sequence = str(fasta_name or "")
        self._chip.setVisible(bool(self._flt.sequence))
        self._chip_label.setText(f"Showing images linked to {self._flt.sequence}")
        self._apply()

    def _visible_items(self) -> list[lb.LibraryItem]:
        return lb.apply_filter(self._items, self._flt)

    def _apply(self) -> None:
        self._render()
        self._rebuild_filters()
        self._show_details(self._item(self._selected))

    def _item(self, image_id: str) -> lb.LibraryItem | None:
        return next((i for i in self._items if i.image.id == image_id), None)

    # ------------------------------------------------------------------ filters

    def _on_search(self, text: str) -> None:
        self._flt.text = text.strip()
        self._render()
        self._rebuild_filters()

    def _rebuild_filters(self) -> None:
        while self._filter_lay.count():
            w = self._filter_lay.takeAt(0).widget()
            if w is not None:
                w.deleteLater()
        for facet, title in (("probe", "Probe"), ("species", "Species"), ("stage", "Cell stage"), ("signals", "Signals")):
            counts = lb.facet_counts(self._items, self._flt, facet)
            if not counts:
                continue
            head = QLabel(f"<b>{title}</b>", self._filter_body)
            self._filter_lay.addWidget(head)
            chosen: set[str] = getattr(self._flt, facet)
            for value, n in counts:
                box = QCheckBox(f"{value}  ({n})", self._filter_body)
                box.setChecked(value in chosen)
                box.toggled.connect(lambda on, f=facet, v=value: self._toggle(f, v, on))
                self._filter_lay.addWidget(box)
        if not self._flt.is_empty():
            clear = QPushButton("Clear filters", self._filter_body)
            clear.clicked.connect(self._clear_filters)
            self._filter_lay.addWidget(clear)

    def _toggle(self, facet: str, value: str, on: bool) -> None:
        chosen: set[str] = getattr(self._flt, facet)
        (chosen.add if on else chosen.discard)(value)
        # the checkbox that sent this signal is replaced: do it after the signal ends
        QTimer.singleShot(0, self._apply)

    def _clear_filters(self) -> None:
        seq = self._flt.sequence
        self._flt = lb.LibraryFilter(sequence=seq)
        self._search.blockSignals(True)
        self._search.clear()
        self._search.blockSignals(False)
        QTimer.singleShot(0, self._apply)

    # ------------------------------------------------------------------ grid

    def _on_group(self, text: str) -> None:
        _write("library_group", text)
        self._render()

    def _clear_grid(self) -> None:
        self._thumbs = {}
        self._cols_for = []
        while self._grid_lay.count():
            w = self._grid_lay.takeAt(0).widget()
            if w is not None:
                w.deleteLater()

    def _render(self) -> None:
        self._clear_grid()
        items = self._visible_items()
        total = len(self._items)
        self._count.setText(f"{len(items)} of {total} images" if len(items) != total else f"{total} image{'s' if total != 1 else ''}")
        if self._folder is None or not items:
            msg = ("No images yet. Use Import slide session for channel TIFFs, or Add image for a finished picture."
                   if not total else "No image matches the filters.")
            label = QLabel(msg, self._grid_body)
            label.setWordWrap(True)
            label.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignLeft)
            self._grid_lay.addWidget(label)
            self._grid_lay.addStretch(1)
            return
        side = self._size.value()
        for group in lb.group_items(items, self._group.currentText()):
            section = QWidget(self._grid_body)
            sl = QVBoxLayout(section)
            sl.setContentsMargins(0, 0, 0, 0)
            sl.setSpacing(4)
            if group.title:
                head = QLabel(f"<b>{group.title}</b>  <span style='color:gray'>{group.subtitle}</span>", section)
                head.setTextFormat(Qt.TextFormat.RichText)
                sl.addWidget(head)
            grid = QGridLayout()
            grid.setSpacing(6)
            grid.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
            thumbs = []
            for item in group.items:
                t = _Thumb(item, self._folder, side, section)
                t.clicked.connect(self._on_clicked)
                t.doubleClicked.connect(self._on_double_clicked)
                self._thumbs[item.image.id] = t
                thumbs.append(t)
            self._cols_for.append((grid, thumbs))
            sl.addLayout(grid)
            self._grid_lay.addWidget(section)
        self._grid_lay.addStretch(1)
        if self._selected in self._thumbs:
            self._thumbs[self._selected].set_selected(True)
        self._relayout()

    def _relayout(self) -> None:
        if not self._cols_for:
            return
        cell = self._size.value() + 24
        cols = max(1, (self._grid_scroll.viewport().width() - 24) // cell)
        for grid, thumbs in self._cols_for:
            for t in thumbs:
                grid.removeWidget(t)
            for n, t in enumerate(thumbs):
                grid.addWidget(t, n // cols, n % cols)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._relayout_timer.start()

    # ------------------------------------------------------------------ selection and details

    def selected_id(self) -> str:
        return self._selected

    def select(self, image_id: str) -> None:
        self._on_clicked(image_id)

    def _on_clicked(self, image_id: str) -> None:
        if self._selected in self._thumbs:
            self._thumbs[self._selected].set_selected(False)
        self._selected = image_id
        if image_id in self._thumbs:
            self._thumbs[image_id].set_selected(True)
        self._show_details(self._item(image_id))

    def _on_double_clicked(self, image_id: str) -> None:
        self._on_clicked(image_id)
        self._open_selected()

    def _show_details(self, item: lb.LibraryItem | None) -> None:
        from repeatviewer.widgets.cyto_images_ui import thumbnail

        while self._details_lay.count():
            w = self._details_lay.takeAt(0).widget()
            if w is not None:
                w.deleteLater()
        if item is None or self._folder is None:
            hint = QLabel("Select an image to see its details.", self._details_body)
            hint.setWordWrap(True)
            self._details_lay.addWidget(hint)
            self._details_lay.addStretch(1)
            return
        side = DETAILS_W - 36
        pic = QLabel(self._details_body)
        pic.setFixedSize(side, int(side * 0.75))
        pic.setAlignment(Qt.AlignmentFlag.AlignCenter)
        pic.setStyleSheet("background:#000; color:#aaa; border-radius:4px;")
        pix = thumbnail(ci.image_file_path(self._folder, item.image), side)
        if pix.isNull():
            pic.setText("file missing")
        else:
            pic.setPixmap(pix.scaled(pic.size(), Qt.AspectRatioMode.KeepAspectRatio, Qt.TransformationMode.SmoothTransformation))
        self._details_lay.addWidget(pic)
        title = QLabel(f"<b>{item.title}</b>", self._details_body)
        title.setWordWrap(True)
        title.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._details_lay.addWidget(title)
        grid = QGridLayout()
        grid.setColumnStretch(1, 1)
        grid.setVerticalSpacing(4)
        for r, (k, v) in enumerate(lb.details_rows(item)):
            key = QLabel(k, self._details_body)
            key.setStyleSheet("color:palette(dark);")
            key.setAlignment(Qt.AlignmentFlag.AlignTop)
            val = QLabel(v, self._details_body)
            val.setWordWrap(True)
            val.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
            grid.addWidget(key, r, 0)
            grid.addWidget(val, r, 1)
        self._details_lay.addLayout(grid)
        row = QHBoxLayout()
        for text, slot in (("Open", self._open_selected), ("Edit", self._edit_selected), ("Remove", self._remove_selected)):
            b = QPushButton(text, self._details_body)
            b.clicked.connect(slot)
            row.addWidget(b)
        self._details_lay.addLayout(row)
        self._details_lay.addStretch(1)

    # ------------------------------------------------------------------ image actions

    def _open_selected(self) -> None:
        item = self._item(self._selected)
        if item is not None:
            self.imageOpened.emit(item.image.id)

    def _edit_selected(self) -> None:
        item = self._item(self._selected)
        if item is None or self._folder is None:
            return
        from repeatviewer.widgets.cyto_images_ui import edit_image

        edit_image(self, self._folder, item.image, allow_no_targets=True)

    def _remove_selected(self) -> None:
        item = self._item(self._selected)
        if item is None or self._folder is None:
            return
        answer = QMessageBox.question(
            self, "Remove image",
            f"Remove \"{item.title}\" from the Library?\nThe picture and its channel copies in the project are deleted; "
            "files outside the project are not touched.",
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        ci.remove_image(self._folder, item.image.id)
        self._selected = ""
        from repeatviewer.widgets.cyto_images_ui import notify_changed

        notify_changed(self._folder)

    def add_image(self) -> None:
        if self._folder is None:
            return
        from repeatviewer.widgets.cyto_images_ui import IMAGE_FILE_FILTER, AttachImageDialog

        chosen, _ = QFileDialog.getOpenFileNames(self, "Add cytogenetic image", str(self._folder), IMAGE_FILE_FILTER)
        for path in chosen:
            if ci.is_image_path(Path(path)):
                AttachImageDialog(self._folder, source=Path(path), targets=[], allow_no_targets=True, parent=self).exec()

    def edit_lists(self) -> None:
        if self._folder is None:
            return
        from repeatviewer.views.cyto_lists_dialog import CytoListsDialog

        CytoListsDialog(self._folder, self).exec()

    # ------------------------------------------------------------------ import

    def start_import(self) -> None:
        if self._folder is None or self._worker is not None:
            return
        from repeatviewer.views.cyto_import_dialog import CytoImportDialog

        dialog = CytoImportDialog(self._folder, self)
        dialog.exec()
        if not dialog.should_run():
            return
        self._worker = _ImportWorker(dialog.request(), self)
        self._worker.progressed.connect(self._on_progress)
        self._worker.finishedWith.connect(self._on_import_done)
        self._import_btn.setEnabled(False)
        self._progress.setRange(0, 0)
        self._progress_label.setText("Importing slide session…")
        self._progress_row.show()
        self._worker.start()

    def _on_progress(self, done: int, total: int) -> None:
        self._progress.setRange(0, max(1, total))
        self._progress.setValue(done)
        self._progress_label.setText(f"Importing slide session: {done} done, {max(0, total - done)} left")

    def _on_stop(self) -> None:
        if self._worker is not None:
            self._worker.request_stop()
            self._progress_label.setText("Stopping after the current cell…")

    def _on_import_done(self, result: object) -> None:
        worker, self._worker = self._worker, None
        if worker is not None:
            worker.wait()
            worker.deleteLater()
        self._progress_row.hide()
        self._import_btn.setEnabled(True)
        if isinstance(result, Exception):
            QMessageBox.critical(self, "Import slide session", f"The import failed:\n{result}")
            return
        slide, images, errors = result
        if self._folder is not None:
            from repeatviewer.widgets.cyto_images_ui import notify_changed

            notify_changed(self._folder)
        if images:
            self.select(images[0].id)
        text = (f"Imported {len(images)} cell{'s' if len(images) != 1 else ''} into slide \"{slide.name}\"."
                if images else "Nothing was imported.")
        if errors:
            text += "\n\n" + "\n".join(errors[:12]) + ("\n..." if len(errors) > 12 else "")
        (QMessageBox.warning if errors or not images else QMessageBox.information)(self, "Import slide session", text)

    def busy(self) -> bool:
        return self._worker is not None
