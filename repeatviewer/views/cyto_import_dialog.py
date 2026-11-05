"""Import: microscope pictures (TIFF, JPG, PNG …) of one slide, grouped into cells."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QSettings, Qt, QTimer
from PySide6.QtWidgets import (
    QGridLayout,
    QInputDialog,
    QMessageBox,
    QCheckBox,
    QDialog,
    QDoubleSpinBox,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QSizePolicy,
    QStackedWidget,
    QTabBar,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import lists as ls
from repeatviewer.backends.cytogenetics.importer import (
    ImportRequest,
    channel_rows,
    plan_cells,
    plan_probes_by_dye,
    scan_tiffs,
)
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cyto_panel_targets import PanelPathFields
from repeatviewer.widgets.cytogenetics_style import form_layout, rounded_box, section_label
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.path_browse import PathBrowseRow

_GROUP = "cytogenetics_import"

HELP_FILES = (
    "<b>Folder</b>: the folder with the microscope pictures of one slide (TIFF in any compression, JPG, PNG, BMP). "
    "Click a folder in the Folders panel or a file in the Files panel to fill it. Tick <b>Include subfolders</b> "
    "when each cell has its own folder.<br>"
    "<b>Cells</b>, decided in this order:<br>"
    "1. A multichannel TIFF (ImageJ, OME, several pages) is one cell, split into its channels.<br>"
    "2. Files named with a dye word (<i>s12_c03_DAPI.tif</i>, <i>s12_c03_FITC.tif</i>) are grouped by the rest of the name.<br>"
    "3. Any other files (<i>1.tif</i>, <i>2.tif</i> ...) are grouped by folder: one folder is one cell, "
    "each file one channel (Channel 1, 2 ...).<br>"
    "4. An RGB merge on its own is split into its red, green and blue planes.<br>"
    "Photoshop files (.psd) are skipped. Dyes and probes that are not in the file names are set later in the "
    "Image tab, once for the whole slide.<br>"
    "<b>Probe per dye</b> appears only for dyes found in file names."
)
HELP_SLIDE = (
    "Details that apply to every cell of this slide. They are kept on the slide and copied to each image, "
    "so the filters in the Library work. Anything left empty is simply not stored."
)
HELP_OPTIONS = (
    "<b>Copy the TIFFs</b>: keeps a copy in the project's cytogenetics folder. The originals on the "
    "microscope drive are never changed. Without copying, the images point to the original files, which "
    "must stay where they are.<br>"
    "<b>Pixel size</b>: read from the TIFF when it is stored there (needed for scale bars). "
    "Set it by hand when the files carry none."
)


def _settings() -> QSettings:
    return QSettings()


def _get(key: str, default: str = "") -> str:
    try:
        return str(_settings().value(f"{_GROUP}/{key}", default))
    except Exception:  # noqa: BLE001
        return default


def _put(key: str, value) -> None:
    try:
        _settings().setValue(f"{_GROUP}/{key}", value)
    except Exception:  # noqa: BLE001
        pass


ADD_ITEM = {"dyes": "Add a dye…", "species": "Add a species…", "techniques": "Add a technique…", "counterstains": "Add a counterstain…",
            "stages": "Add a cell stage…", "probes": "Add a probe…"}


def _combo(parent: QWidget, values: list[str], current: str = "", editable: bool = True, *,
           key: str = "", folder: Path | None = None, lists: dict | None = None, dye: str = "",
           on_change=None) -> NoWheelComboBox:
    """A menu of list values that ends with "Add a …": a new entry is saved to the folder's lists right away."""
    box = NoWheelComboBox(parent)
    box.setEditable(editable)
    box.addItems([v for v in values if v])
    extra = ADD_ITEM.get(key, "")
    if extra:
        box.addItem(extra)
    box.setCurrentText(current or "")
    box.setProperty("previous", current or "")

    def _activated(_i: int) -> None:
        text = box.currentText()
        if not extra or text != extra:
            box.setProperty("previous", text)
            return
        box.setCurrentText(str(box.property("previous") or ""))
        what = extra[len("Add a "):].rstrip("…")
        name, ok = QInputDialog.getText(box.window(), extra.rstrip("…"), f"Name of the {what}:")
        name = name.strip()
        if not ok or not name:
            return
        values: dict = {"dye": dye} if key == "probes" and dye else {}
        if key == "dyes":
            from PySide6.QtGui import QColor
            from PySide6.QtWidgets import QColorDialog

            colour = QColorDialog.getColor(QColor("#33ff33"), box.window(), f"Colour of {name}")
            values = {"color": colour.name() if colour.isValid() else "#ffffff", "words": name.lower()}
        if folder is not None and lists is not None:
            try:
                ls.add_entry(lists, key, name, **values)
                ls.save_lists(folder, lists)
            except (OSError, ValueError) as exc:
                QMessageBox.information(box.window(), extra.rstrip("…"), str(exc))
                return
        if box.findText(name) < 0:
            box.insertItem(box.count() - 1, name)
        box.setCurrentText(name)
        box.setProperty("previous", name)

    box.activated.connect(_activated)
    if on_change is not None:
        box.activated.connect(lambda _i: on_change(box.currentText()))
        if box.isEditable():
            box.lineEdit().editingFinished.connect(lambda: on_change(box.currentText()))
    return box


def _dot(colour: str):
    from PySide6.QtGui import QBrush, QColor, QPainter, QPen, QPixmap

    pix = QPixmap(14, 14)
    pix.fill(Qt.GlobalColor.transparent)
    p = QPainter(pix)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    p.setPen(QPen(QColor("#b0b0b5"), 1))
    p.setBrush(QBrush(QColor(colour or "#ffffff")))
    p.drawEllipse(1, 1, 12, 12)
    p.end()
    return pix


def _cell(widget: QWidget, parent: QWidget) -> QWidget:
    """A menu or label inside a table cell, with a little air around it."""
    host = QWidget(parent)
    lay = QHBoxLayout(host)
    lay.setContentsMargins(4, 2, 8, 2)
    lay.addWidget(widget)
    return host


class CytoImportDialog(PanelPathFields, QDialog):
    """Slide import in the layout of the app's tool windows. Run imports, Save only remembers the settings.

    Not modal: a click on a folder in the Folders panel (or a file in the Files panel) fills the Folder field.
    """

    def __init__(self, folder: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Import")
        self.setWindowModality(Qt.WindowModality.NonModal)
        self._folder = Path(folder).expanduser().resolve()
        self._lists = ls.load_lists(self._folder)
        self._run = False
        self._groups = []
        self._notes: list[str] = []
        self._probe_boxes: dict[str, NoWheelComboBox] = {}
        self._saved_probes: dict[str, str] = {}
        try:
            self._saved_probes = dict(json.loads(_get("probes_by_dye", "{}")))
        except (ValueError, TypeError):
            pass

        self._tab_bar = QTabBar(self)
        self._tab_bar.setExpanding(False)
        self._tab_bar.setDocumentMode(True)
        self._tab_bar.addTab("Import")
        self._tab_bar.addTab("Options")
        self._stack = QStackedWidget(self)
        self._stack.addWidget(self._build_import_page())
        self._stack.addWidget(self._build_options_tab())
        self._tab_bar.currentChanged.connect(self._stack.setCurrentIndex)

        cancel = QPushButton("Cancel", self)
        save = QPushButton("Save", self)
        run = QPushButton("Run", self)
        run.setDefault(True)
        cancel.clicked.connect(self.reject)
        save.clicked.connect(self._on_save)
        run.clicked.connect(self._on_run)
        buttons = QHBoxLayout()
        buttons.setSpacing(8)
        buttons.addStretch(1)
        for b in (cancel, save, run):
            buttons.addWidget(b)
        root = QVBoxLayout(self)
        root.setContentsMargins(22, 18, 22, 16)
        root.setSpacing(8)
        root.addWidget(self._tab_bar)
        root.addWidget(self._stack, 1)
        root.addLayout(buttons)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(300)
        self._timer.timeout.connect(self._rescan)
        self._source.line_edit().textChanged.connect(lambda _t: self._timer.start())
        self._recursive.toggled.connect(lambda _on: self._timer.start())
        self.watch_path_fields(self._source)
        self._rescan()
        self.resize(640, self.sizeHint().height())

    # ------------------------------------------------------------------ pages

    def _build_import_page(self) -> QWidget:
        page = QWidget(self)
        lay = QVBoxLayout(page)
        lay.setContentsMargins(16, 10, 16, 6)
        lay.setSpacing(6)
        form = form_layout()
        lists = self._lists
        files_head = QHBoxLayout()
        files_head.addWidget(section_label("Files", page))
        files_head.addStretch(1)
        files_head.addWidget(make_help_button(page, "Import", HELP_FILES + "<br>" + HELP_SLIDE))
        form.addRow(files_head)
        self._source = PathBrowseRow(
            page,
            initial=_get("source") or None,
            placeholder="Click a folder in the Folders panel, or Browse",
            browse_title="Choose the folder of the slide",
            mode="directory",
        )
        form.addRow("Folder", self._source)
        self._recursive = QCheckBox("Include subfolders (one folder per cell)", page)
        self._recursive.setChecked(_get("recursive", "true") == "true")
        form.addRow("", self._recursive)
        self._table = QTableWidget(0, 4, page)
        self._table.setHorizontalHeaderLabels(["Cell", "Channels", "Size", "From"])
        self._table.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._table.setTextElideMode(Qt.TextElideMode.ElideRight)
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self._table.verticalHeader().setVisible(False)
        self._table.verticalHeader().setDefaultSectionSize(24)
        header = self._table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        header.setStretchLastSection(True)
        header.setMinimumSectionSize(60)
        header.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self._table.setShowGrid(False)
        self._table.setAlternatingRowColors(True)
        self._table_box = rounded_box(self._table, page)
        form.addRow("Cells", self._table_box)

        form.addRow(section_label("Slide", page))
        self._slide_name = QLineEdit(page)
        self._slide_date = QLineEdit(page)
        self._slide_date.setPlaceholderText("YYYY-MM-DD")
        self._species = _combo(page, ls.names(lists, "species"), _get("species"), key="species", folder=self._folder, lists=lists)
        self._accession = QLineEdit(_get("accession"), page)
        self._technique = _combo(page, ls.names(lists, "techniques"), _get("technique", "FISH"), key="techniques", folder=self._folder, lists=lists)
        self._counterstain = _combo(page, ls.names(lists, "counterstains"), _get("counterstain", "DAPI"), key="counterstains", folder=self._folder, lists=lists)
        self._stage = _combo(page, ls.names(lists, "stages"), _get("stage", "Mitotic metaphase"), key="stages", folder=self._folder, lists=lists)
        lay.addLayout(form)
        # Slide: a 3 x 3 grid, small label above each field, every cell the same width; probe menus fill the next cells
        self._slide_grid = QGridLayout()
        self._slide_grid.setContentsMargins(0, 2, 0, 0)
        self._slide_grid.setHorizontalSpacing(16)
        self._slide_grid.setVerticalSpacing(10)
        for c in range(3):
            self._slide_grid.setColumnStretch(c, 1)
        self._slide_cells: list[QWidget] = []
        for text, w in (("Name", self._slide_name), ("Date", self._slide_date), ("Accession / line", self._accession),
                        ("Species", self._species), ("Technique", self._technique), ("Counterstain", self._counterstain),
                        ("Cell stage", self._stage)):
            self._add_slide_cell(text, w, page)
        self._slide_base = len(self._slide_cells)
        grid_host = QWidget(page)
        grid_host.setLayout(self._slide_grid)
        form.addRow("", grid_host)  # in the field column: same left and right edges as the Cells table
        # Channels: one row per channel position of the slide, dye and probe set once for the whole slide
        form.addRow(section_label("Channels", page))
        self._chan_table = QTableWidget(0, 4, page)
        self._chan_table.setHorizontalHeaderLabels(["", "Found as", "Dye", "Probe"])
        self._chan_table.verticalHeader().setVisible(False)
        self._chan_table.verticalHeader().setDefaultSectionSize(32)
        self._chan_table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self._chan_table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._chan_table.setShowGrid(False)
        self._chan_table.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        chan_header = self._chan_table.horizontalHeader()
        chan_header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)  # columns stay resizable
        chan_header.setStretchLastSection(True)
        chan_header.setDefaultAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        self._chan_table.setColumnWidth(0, 30)
        self._chan_box = rounded_box(self._chan_table, page)
        form.addRow("", self._chan_box)
        self._chan_choice: dict[int, dict] = {}  # position -> {"dye", "probe"} chosen by the user
        self._chan_rows: list = []
        return page

    def _add_slide_cell(self, text: str, field: QWidget, parent: QWidget) -> None:
        cell = QWidget(parent)
        v = QVBoxLayout(cell)
        v.setContentsMargins(0, 0, 0, 0)
        v.setSpacing(4)
        label = QLabel(text, cell)
        f = label.font()
        f.setPointSizeF(f.pointSizeF() - 2)
        label.setFont(f)
        label.setStyleSheet("color: gray;")
        field.setParent(cell)
        field.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        field.setMinimumWidth(0)
        field.setMaximumWidth(16777215)
        v.addWidget(label)
        v.addWidget(field)
        n = len(self._slide_cells)
        self._slide_grid.addWidget(cell, n // 3, n % 3)
        self._slide_cells.append(cell)

    def _build_options_tab(self) -> QWidget:
        page = QWidget(self)
        lay = QVBoxLayout(page)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.setSpacing(14)
        lay.setAlignment(Qt.AlignmentFlag.AlignTop)
        files_box = QGroupBox("Files", page)
        fl = QVBoxLayout(files_box)
        fl.setContentsMargins(10, 8, 10, 8)
        self._copy = QCheckBox("Copy the pictures into the project (originals are never changed)", files_box)
        self._copy.setChecked(_get("copy", "true") == "true")
        fl.addWidget(self._copy)
        self._merges = QCheckBox("Import merged and other pictures too (each as its own cell)", files_box)
        self._merges.setToolTip("A merge is a colour picture made from several channels (for example in Photoshop). "
                                "Other pictures are single files of another size than the channels next to them (a figure, a crop). "
                                "Both are left out by default: the channels are imported instead.")
        self._merges.setChecked(_get("merges", "false") == "true")
        self._merges.toggled.connect(lambda _on: self._timer.start())
        fl.addWidget(self._merges)
        pix_box = QGroupBox("Pixel size", page)
        pl = QVBoxLayout(pix_box)
        pl.setContentsMargins(10, 8, 10, 8)
        self._pix_auto = QRadioButton("From the TIFF", pix_box)
        self._pix_set = QRadioButton("Set by hand:", pix_box)
        self._pix_value = QDoubleSpinBox(pix_box)
        self._pix_value.setDecimals(4)
        self._pix_value.setRange(0.0001, 100.0)
        self._pix_value.setSuffix(" µm per pixel")
        try:
            self._pix_value.setValue(float(_get("pixel_um", "0.0645")))
        except ValueError:
            self._pix_value.setValue(0.0645)
        (self._pix_set if _get("pixel_set", "false") == "true" else self._pix_auto).setChecked(True)
        self._pix_value.setEnabled(self._pix_set.isChecked())
        self._pix_set.toggled.connect(self._pix_value.setEnabled)
        pl.addWidget(self._pix_auto)
        row = QHBoxLayout()
        row.addWidget(self._pix_set)
        row.addWidget(self._pix_value)
        row.addStretch(1)
        pl.addLayout(row)
        lay.addWidget(files_box)
        lay.addWidget(pix_box)
        help_row = QHBoxLayout()
        help_row.addWidget(make_help_button(page, "Import options", HELP_OPTIONS))
        help_row.addStretch(1)
        lay.addLayout(help_row)
        return page

    # ------------------------------------------------------------------ scanning

    def _source_path(self) -> Path | None:
        p = self._source.path()
        return Path(p) if p else None

    def _rescan(self) -> None:
        src = self._source_path()
        files = scan_tiffs(src, self._recursive.isChecked()) if src else []
        self._groups, notes = plan_cells(files, self._lists, merges=self._merges.isChecked())
        self._notes = notes
        self._table.setRowCount(len(self._groups))
        for r, g in enumerate(self._groups):
            sources = sorted({p.name for p, _d, _i in g.entries()}, key=lambda n: n.lower())
            origin = sources[0] if len(sources) == 1 else f"{len(sources)} files"
            named = [d for d in g.dyes if d]
            channels = ", ".join(named) if named and len(named) == g.size else f"{g.size} channel{'s' if g.size != 1 else ''}"
            size = ""
            try:
                from repeatviewer.backends.cytogenetics import channels as chn

                info = chn.inspect_tiff(g.entries()[0][0])
                size = f"{info.width} × {info.height}" if info.width else ""
            except Exception:  # noqa: BLE001
                pass
            for c, text in enumerate((g.key, channels, size, origin)):
                item = QTableWidgetItem(text)
                item.setToolTip(", ".join(sources) if c == 3 else text)
                self._table.setItem(r, c, item)
        rows = max(2, min(8, len(self._groups)))
        self._table_box.setFixedHeight(self._table.horizontalHeader().sizeHint().height() + rows * 24 + 8)
        header = self._table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)  # columns stay resizable
        header.setStretchLastSection(True)
        self._table.resizeColumnsToContents()
        rest = sum(self._table.columnWidth(c) + 12 for c in (1, 2, 3))
        for c in (1, 2, 3):
            self._table.setColumnWidth(c, self._table.columnWidth(c) + 12)
        self._table.setColumnWidth(0, max(self._table.columnWidth(0) + 12, self._table.viewport().width() - rest))
        if src is not None and not self._slide_name.text().strip():
            self._slide_name.setText(src.name)
        if src is not None and not self._slide_date.text().strip() and self._groups:
            self._slide_date.setText(ci.image_file_date(self._groups[0].entries()[0][0]))
        self._fill_probe_form()

    def _fill_probe_form(self) -> None:
        """The Channels table: what was found per channel position, with dye and probe menus."""
        rows = channel_rows(self._groups, self._lists)
        counter = {n.lower() for n in ls.names(self._lists, "counterstains")} | {"dapi", "hoechst"}
        dyes = ls.names(self._lists, "dyes")
        probes = ls.names(self._lists, "probes")
        t = self._chan_table
        t.setRowCount(len(rows))
        self._chan_rows = rows
        for r, row in enumerate(rows):
            chosen = self._chan_choice.get(row.position, {})
            row.dye = chosen.get("dye", row.dye)
            row.probe = chosen.get("probe", row.probe or (self._saved_probes.get(row.dye, "") if row.dye else ""))
            dot = QLabel(t)
            dot.setAlignment(Qt.AlignmentFlag.AlignCenter)
            dot.setPixmap(_dot(row.colour))
            t.setCellWidget(r, 0, dot)
            item = QTableWidgetItem(row.found_as)
            item.setToolTip(f"Channel {row.position}: {row.found_as}")
            t.setItem(r, 1, item)
            dye_box = _combo(t, dyes, row.dye, key="dyes", folder=self._folder, lists=self._lists,
                             on_change=lambda text, pos=row.position: self._on_channel(pos, "dye", text))
            t.setCellWidget(r, 2, _cell(dye_box, t))
            if row.dye.lower() in counter:
                note = QLabel("Counterstain", t)
                note.setStyleSheet("color: gray;")
                t.setCellWidget(r, 3, _cell(note, t))
            else:
                probe_box = _combo(t, probes, row.probe, key="probes", folder=self._folder, lists=self._lists,
                                   dye=row.dye,
                                   on_change=lambda text, pos=row.position: self._on_channel(pos, "probe", text))
                probe_box.lineEdit().setPlaceholderText("Choose or add…")
                t.setCellWidget(r, 3, _cell(probe_box, t))
        shown = max(1, min(6, len(rows)))
        self._chan_box.setFixedHeight(t.horizontalHeader().sizeHint().height() + shown * 32 + 8)
        self._chan_box.setVisible(bool(rows))
        if rows and t.columnWidth(1) < 60:
            width = max(t.viewport().width(), 420) - 30
            for c, w in ((1, width * 0.36), (2, width * 0.32)):
                t.setColumnWidth(c, int(w))

    def _on_channel(self, position: int, field: str, text: str) -> None:
        text = text.strip()
        if text.startswith("Add a "):
            return
        choice = self._chan_choice.setdefault(position, {})
        if choice.get(field) == text:
            return
        choice[field] = text
        if field == "dye":
            QTimer.singleShot(0, self._fill_probe_form)  # counterstain or probe menu, and the probe's dye

    # ------------------------------------------------------------------ result

    def _channel_plan(self) -> list[dict]:
        return [{"dye": row.dye, "probe": self._chan_choice.get(row.position, {}).get("probe", row.probe)}
                for row in self._chan_rows]

    def _probes_by_dye(self) -> dict[str, str]:
        return plan_probes_by_dye(self._channel_plan())

    def request(self) -> ImportRequest:
        return ImportRequest(
            folder=self._folder,
            source=self._source_path() or Path(""),
            recursive=self._recursive.isChecked(),
            slide_name=self._slide_name.text().strip(),
            slide_values={
                "date": self._slide_date.text().strip(),
                "species": self._species.currentText().strip(),
                "accession": self._accession.text().strip(),
                "technique": self._technique.currentText().strip(),
                "counterstain": self._counterstain.currentText().strip(),
            },
            stage=self._stage.currentText().strip(),
            probes_by_dye=self._probes_by_dye(),
            copy=self._copy.isChecked(),
            pixel_size_um=self._pix_value.value() if self._pix_set.isChecked() else None,
            merges=self._merges.isChecked(),
            channel_plan=self._channel_plan(),
        )

    def should_run(self) -> bool:
        return self._run

    def _remember(self) -> None:
        src = self._source_path()
        _put("source", str(src) if src else "")
        _put("recursive", "true" if self._recursive.isChecked() else "false")
        for key, box in (("species", self._species), ("technique", self._technique),
                         ("counterstain", self._counterstain), ("stage", self._stage)):
            _put(key, box.currentText().strip())
        _put("accession", self._accession.text().strip())
        _put("copy", "true" if self._copy.isChecked() else "false")
        _put("merges", "true" if self._merges.isChecked() else "false")
        _put("pixel_set", "true" if self._pix_set.isChecked() else "false")
        _put("pixel_um", self._pix_value.value())
        _put("probes_by_dye", json.dumps({**self._saved_probes, **self._probes_by_dye()}))

    def _on_save(self) -> None:
        self._remember()
        self.accept()

    def _on_run(self) -> None:
        from PySide6.QtWidgets import QMessageBox

        if not self._groups:
            text = "There is nothing to import in this folder."
            if self._notes:
                text += "\n\n" + "\n".join(self._notes[:8])
            QMessageBox.information(self, "Import", text)
            return
        self._remember()
        self._run = True
        self.accept()
