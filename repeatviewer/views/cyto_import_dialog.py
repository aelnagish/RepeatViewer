"""Import: microscope pictures (TIFF, JPG, PNG …) of one slide, grouped into cells."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QSettings, Qt, QTimer
from PySide6.QtWidgets import (
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
from repeatviewer.backends.cytogenetics.importer import ImportRequest, plan_cells, scan_tiffs
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cyto_panel_targets import PanelPathFields
from repeatviewer.widgets.cytogenetics_style import form_layout, rounded_box, section_label
from repeatviewer.widgets.help_button import field_with_help, make_help_button
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


def _combo(parent: QWidget, values: list[str], current: str = "", editable: bool = True) -> NoWheelComboBox:
    box = NoWheelComboBox(parent)
    box.setEditable(editable)
    box.addItem("")
    box.addItems(values)
    if current:
        box.setCurrentText(current)
    return box


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
        form.addRow(section_label("Files", page))
        self._source = PathBrowseRow(
            page,
            initial=_get("source") or None,
            placeholder="Click a folder in the Folders panel, or Browse",
            browse_title="Choose the folder of the slide",
            mode="directory",
        )
        form.addRow("Folder", field_with_help(self._source, "Import", HELP_FILES, stretch=True))
        self._recursive = QCheckBox("Include subfolders (one folder per cell)", page)
        self._recursive.setChecked(_get("recursive", "true") == "true")
        form.addRow("", self._recursive)
        self._table = QTableWidget(0, 3, page)
        self._table.setHorizontalHeaderLabels(["Cell", "Channels", "From"])
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self._table.verticalHeader().setVisible(False)
        self._table.verticalHeader().setDefaultSectionSize(24)
        header = self._table.horizontalHeader()
        header.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        header.setStretchLastSection(True)
        header.setMinimumSectionSize(60)
        self._table.setShowGrid(False)
        self._table.setAlternatingRowColors(True)
        self._table_box = rounded_box(self._table, page)
        form.addRow("Cells", self._table_box)
        self._info = QLabel(page)
        self._info.setWordWrap(True)
        self._info.setStyleSheet("color: gray;")
        form.addRow("", self._info)

        form.addRow(section_label("Slide", page))
        self._slide_name = QLineEdit(page)
        self._slide_date = QLineEdit(page)
        self._slide_date.setPlaceholderText("YYYY-MM-DD")
        self._species = _combo(page, ls.names(lists, "species"), _get("species"))
        self._accession = QLineEdit(_get("accession"), page)
        self._technique = _combo(page, ls.names(lists, "techniques"), _get("technique", "FISH"))
        self._counterstain = _combo(page, ls.names(lists, "counterstains"), _get("counterstain", "DAPI"))
        self._stage = _combo(page, ls.names(lists, "stages"), _get("stage", "Mitotic metaphase"))
        name_row = QHBoxLayout()
        name_row.setSpacing(6)
        name_row.addWidget(self._slide_name, 1)
        name_row.addWidget(make_help_button(page, "Slide", HELP_SLIDE))
        form.addRow("Name", name_row)
        for box in (self._species, self._technique, self._counterstain, self._stage):
            box.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        for label, w in (("Date", self._slide_date), ("Species", self._species), ("Accession / line", self._accession),
                         ("Technique", self._technique), ("Counterstain", self._counterstain), ("Cell stage", self._stage)):
            form.addRow(label, w)
        self._probe_form = form
        self._probe_rows_start = form.rowCount()
        lay.addLayout(form)
        return page

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
        top = QHBoxLayout()
        top.addStretch(1)
        top.addWidget(make_help_button(page, "Import options", HELP_OPTIONS))
        lay.addLayout(top)
        lay.addWidget(files_box)
        lay.addWidget(pix_box)
        return page

    # ------------------------------------------------------------------ scanning

    def _source_path(self) -> Path | None:
        p = self._source.path()
        return Path(p) if p else None

    def _rescan(self) -> None:
        src = self._source_path()
        files = scan_tiffs(src, self._recursive.isChecked()) if src else []
        self._groups, notes = plan_cells(files, self._lists)
        self._notes = notes
        self._table.setRowCount(len(self._groups))
        for r, g in enumerate(self._groups):
            names = [d or f"Channel {k + 1}" for k, d in enumerate(g.dyes)]
            sources = sorted({p.name for p, _d, _i in g.entries()}, key=lambda n: n.lower())
            origin = sources[0] if len(sources) == 1 else f"{len(sources)} files"
            for c, text in enumerate((g.key, ", ".join(names), origin)):
                self._table.setItem(r, c, QTableWidgetItem(text))
        if src is None:
            self._info.setText("Choose the folder of the slide.")
        elif not self._groups and notes:
            self._info.setText("These files could not be read:<br>" + "<br>".join(notes[:6]))
        elif not self._groups:
            self._info.setText("No pictures here (TIFF, JPG, PNG or BMP)." + ("" if self._recursive.isChecked() else " Tick Include subfolders when each cell has its own folder."))
        else:
            n_ch = sum(g.size for g in self._groups)
            text = (f"{len(self._groups)} cell{'s' if len(self._groups) != 1 else ''}, {n_ch} channel{'s' if n_ch != 1 else ''}. "
                    "Dyes and probes that are not in the file names are set later in the Image tab.")
            if notes:
                text += "<br>Not read: " + "; ".join(notes[:3])
            self._info.setText(text)
        rows = max(2, min(8, len(self._groups)))
        self._table_box.setFixedHeight(self._table.horizontalHeader().sizeHint().height() + rows * 24 + 8)
        self._table.resizeColumnsToContents()
        if src is not None and not self._slide_name.text().strip():
            self._slide_name.setText(src.name)
        if src is not None and not self._slide_date.text().strip() and self._groups:
            self._slide_date.setText(ci.image_file_date(self._groups[0].entries()[0][0]))
        self._fill_probe_form()

    def _fill_probe_form(self) -> None:
        dyes = sorted({d for g in self._groups for d in g.dyes if d})
        current = {d: b.currentText() for d, b in self._probe_boxes.items()}
        while self._probe_form.rowCount() > self._probe_rows_start:
            self._probe_form.removeRow(self._probe_rows_start)
        self._probe_boxes = {}
        probes = ls.names(self._lists, "probes")
        for dye in dyes:
            if dye == "DAPI":
                continue  # the counterstain carries no probe
            default = next((str(e["name"]) for e in self._lists.get("probes", []) if e.get("dye") == dye), "")
            box = _combo(self, probes, current.get(dye) or self._saved_probes.get(dye, "") or default)
            box.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            self._probe_boxes[dye] = box
            self._probe_form.addRow(f"Probe in {dye}", box)

    # ------------------------------------------------------------------ result

    def _probes_by_dye(self) -> dict[str, str]:
        return {d: b.currentText().strip() for d, b in self._probe_boxes.items() if b.currentText().strip()}

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
