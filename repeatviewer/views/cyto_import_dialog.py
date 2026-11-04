"""Import slide session: channel TIFFs from the microscope, grouped into cells."""

from __future__ import annotations

import json
from pathlib import Path

from PySide6.QtCore import QSettings, Qt, QTimer
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDoubleSpinBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QTableWidget,
    QTableWidgetItem,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import lists as ls
from repeatviewer.backends.cytogenetics.importer import ImportRequest, group_channel_files, scan_tiffs
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.path_browse import PathBrowseRow

_GROUP = "cytogenetics_import"

HELP_FILES = (
    "<b>Folder</b>: the folder with the microscope TIFFs of one slide.<br>"
    "<b>Cells</b>: files are grouped by their name without the dye word, so "
    "<i>s12_c03_DAPI.tif</i>, <i>s12_c03_FITC.tif</i> and <i>s12_c03_Cy3.tif</i> become cell "
    "<i>s12_c03</i> with three channels.<br>"
    "<b>Dye words</b> (DAPI, FITC, 488, Cy3 ...) come from <b>Edit lists</b> under Dyes and channels.<br>"
    "<b>Probe per dye</b>: the probe shown in that channel. A probe with a sequence in the Probes list "
    "also links the images to that FASTA file."
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


class CytoImportDialog(QDialog):
    """Slide session import. Run imports, Save only remembers the settings."""

    def __init__(self, folder: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Import slide session")
        self.setMinimumWidth(760)
        self._folder = Path(folder).expanduser().resolve()
        self._lists = ls.load_lists(self._folder)
        self._run = False
        self._groups = []
        self._probe_boxes: dict[str, NoWheelComboBox] = {}
        self._saved_probes: dict[str, str] = {}
        try:
            self._saved_probes = dict(json.loads(_get("probes_by_dye", "{}")))
        except (ValueError, TypeError):
            pass

        self._tabs = QTabWidget(self)
        self._tabs.addTab(self._build_files_tab(), "Files")
        self._tabs.addTab(self._build_slide_tab(), "Slide details")
        self._tabs.addTab(self._build_options_tab(), "Options")

        cancel = QPushButton("Cancel", self)
        save = QPushButton("Save", self)
        run = QPushButton("Run", self)
        run.setDefault(True)
        cancel.clicked.connect(self.reject)
        save.clicked.connect(self._on_save)
        run.clicked.connect(self._on_run)
        buttons = QHBoxLayout()
        buttons.addStretch(1)
        for b in (cancel, save, run):
            buttons.addWidget(b)
        root = QVBoxLayout(self)
        root.addWidget(self._tabs)
        root.addLayout(buttons)

        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(300)
        self._timer.timeout.connect(self._rescan)
        self._source.line_edit().textChanged.connect(lambda _t: self._timer.start())
        self._recursive.toggled.connect(lambda _on: self._timer.start())
        self._rescan()

    # ------------------------------------------------------------------ tabs

    def _build_files_tab(self) -> QWidget:
        page = QWidget(self)
        lay = QVBoxLayout(page)
        row = QHBoxLayout()
        row.addWidget(QLabel("Folder", page))
        self._source = PathBrowseRow(
            page,
            initial=_get("source") or None,
            placeholder="Folder with the channel TIFFs of one slide",
            browse_title="Choose the folder with the channel TIFFs",
            mode="directory",
        )
        row.addWidget(self._source, 1)
        row.addWidget(make_help_button(page, "Import slide session", HELP_FILES))
        lay.addLayout(row)
        self._recursive = QCheckBox("Include subfolders", page)
        self._recursive.setChecked(_get("recursive", "false") == "true")
        lay.addWidget(self._recursive)
        self._info = QLabel(page)
        self._info.setWordWrap(True)
        lay.addWidget(self._info)
        self._table = QTableWidget(0, 3, page)
        self._table.setHorizontalHeaderLabels(["Cell", "Channels", "Files"])
        self._table.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._table.setSelectionMode(QTableWidget.SelectionMode.NoSelection)
        self._table.verticalHeader().setVisible(False)
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.setColumnWidth(0, 200)
        self._table.setColumnWidth(1, 260)
        self._table.setMinimumHeight(200)
        lay.addWidget(self._table, 1)
        self._probe_form = QFormLayout()
        lay.addLayout(self._probe_form)
        return page

    def _build_slide_tab(self) -> QWidget:
        page = QWidget(self)
        form = QFormLayout(page)
        lists = self._lists
        self._slide_name = QLineEdit(page)
        self._slide_date = QLineEdit(page)
        self._slide_date.setPlaceholderText("YYYY-MM-DD")
        self._species = _combo(page, ls.names(lists, "species"), _get("species"))
        self._accession = QLineEdit(_get("accession"), page)
        self._technique = _combo(page, ls.names(lists, "techniques"), _get("technique", "FISH"))
        self._counterstain = _combo(page, ls.names(lists, "counterstains"), _get("counterstain", "DAPI"))
        self._stage = _combo(page, ls.names(lists, "stages"), _get("stage", "Mitotic metaphase"))
        top = QHBoxLayout()
        top.addWidget(QLabel("Details of the slide", page))
        top.addStretch(1)
        top.addWidget(make_help_button(page, "Slide details", HELP_SLIDE))
        form.addRow(top)
        for label, w in (("Slide name", self._slide_name), ("Date", self._slide_date), ("Species", self._species),
                         ("Accession / line", self._accession), ("Technique", self._technique),
                         ("Counterstain", self._counterstain), ("Cell stage", self._stage)):
            form.addRow(label, w)
        return page

    def _build_options_tab(self) -> QWidget:
        page = QWidget(self)
        lay = QVBoxLayout(page)
        top = QHBoxLayout()
        top.addStretch(1)
        top.addWidget(make_help_button(page, "Import options", HELP_OPTIONS))
        lay.addLayout(top)
        self._copy = QCheckBox("Copy the TIFFs into the project (originals are never changed)", page)
        self._copy.setChecked(_get("copy", "true") == "true")
        lay.addWidget(self._copy)
        self._pix_auto = QRadioButton("Pixel size from the TIFF", page)
        self._pix_set = QRadioButton("Set the pixel size:", page)
        self._pix_value = QDoubleSpinBox(page)
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
        lay.addWidget(self._pix_auto)
        row = QHBoxLayout()
        row.addWidget(self._pix_set)
        row.addWidget(self._pix_value)
        row.addStretch(1)
        lay.addLayout(row)
        lay.addStretch(1)
        return page

    # ------------------------------------------------------------------ scanning

    def _source_path(self) -> Path | None:
        p = self._source.path()
        return Path(p) if p else None

    def _rescan(self) -> None:
        src = self._source_path()
        files = scan_tiffs(src, self._recursive.isChecked()) if src else []
        self._groups, unmatched = group_channel_files(files, self._lists)
        tiffs = [p for p in unmatched if p.suffix.lower() in (".tif", ".tiff")]
        self._table.setRowCount(len(self._groups))
        for r, g in enumerate(self._groups):
            for c, text in enumerate((g.key, ", ".join(sorted(g.dyes)), str(len(g.files)))):
                self._table.setItem(r, c, QTableWidgetItem(text))
        if src is None:
            self._info.setText("Choose the folder with the channel TIFFs.")
        elif not self._groups:
            self._info.setText("No channel TIFFs with a known dye word in their name were found.")
        else:
            text = f"{len(self._groups)} cell{'s' if len(self._groups) != 1 else ''} with {sum(len(g.files) for g in self._groups)} files."
            if tiffs:
                names = ", ".join(p.name for p in tiffs[:3]) + (" ..." if len(tiffs) > 3 else "")
                text += f" {len(tiffs)} TIFF{'s' if len(tiffs) != 1 else ''} skipped (no known dye word, or the same dye twice in one cell): {names}"
            self._info.setText(text)
        if src is not None and not self._slide_name.text().strip():
            self._slide_name.setText(src.name)
        if src is not None and not self._slide_date.text().strip() and self._groups:
            self._slide_date.setText(ci.image_file_date(self._groups[0].files[0][0]))
        self._fill_probe_form()

    def _fill_probe_form(self) -> None:
        dyes = sorted({d for g in self._groups for d in g.dyes})
        current = {d: b.currentText() for d, b in self._probe_boxes.items()}
        while self._probe_form.rowCount():
            self._probe_form.removeRow(0)
        self._probe_boxes = {}
        probes = ls.names(self._lists, "probes")
        for dye in dyes:
            if dye == "DAPI":
                continue  # the counterstain carries no probe
            default = next((str(e["name"]) for e in self._lists.get("probes", []) if e.get("dye") == dye), "")
            box = _combo(self, probes, current.get(dye) or self._saved_probes.get(dye, "") or default)
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
            QMessageBox.information(self, "Import slide session", "There are no cells to import. Choose a folder with channel TIFFs first.")
            return
        self._remember()
        self._run = True
        self.accept()
