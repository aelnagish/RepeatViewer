"""New Probe window of the Signals tab: name, sequence, dye, colour and the channel it is on."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QIcon, QPixmap
from PySide6.QtWidgets import (
    QCheckBox,
    QColorDialog,
    QDialog,
    QGridLayout,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.cytogenetics import lists as lists_mod
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cyto_panel_targets import PanelPathFields

ADD_DYE = "Add New Dye…"
FIELD_W = 168


def swatch_icon(colour: str, size: int = 14) -> QIcon:
    pix = QPixmap(size, size)
    pix.fill(QColor(colour or "#ffffff"))
    return QIcon(pix)


class NewProbeDialog(PanelPathFields, QDialog):
    """Cancel / Save. ``channels`` are (key, text) of the channels the probe can go on, the first preselected."""

    def __init__(self, folder: Path, channels: list[tuple[str, str]], slide_name: str = "",
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        from repeatviewer.widgets.path_browse import PathBrowseRow

        self.setWindowTitle("New Probe")
        self.setWindowFlag(Qt.WindowType.Tool, True)  # stays above the RepeatViewer window
        self.setAttribute(Qt.WidgetAttribute.WA_MacAlwaysShowToolWindow, True)
        self._folder = Path(folder)
        self._lists = lists_mod.load_lists(self._folder)
        self._colour = ""

        self._name = QLineEdit(self)
        self._name.setFixedWidth(FIELD_W)
        self._seq = PathBrowseRow(self, placeholder="Click a FASTA in Files", browse_title="Sequence of the probe",
                                  start_dir=str(self._folder), elide_middle=True,
                                  file_filter="FASTA (*.fa *.fasta *.fna *.fas);;All files (*)")
        self._seq.panel_suffixes = (".fa", ".fasta", ".fna", ".fas")
        self.watch_path_fields(self._seq)
        hint = QLabel("Or click a FASTA in Files. Optional.", self)
        hint.setStyleSheet("color: gray;")
        self._dye = NoWheelComboBox(self)
        self._dye.setFixedWidth(FIELD_W)
        self._fill_dyes("")
        self._dye.activated.connect(self._on_dye)
        self._colour_btn = QPushButton(self)
        self._colour_btn.setFixedWidth(FIELD_W)
        self._colour_btn.clicked.connect(self._pick_colour)
        self._channel = NoWheelComboBox(self)
        self._channel.setFixedWidth(FIELD_W)
        for key, text in channels:
            self._channel.addItem(text, key)
        self._every = QCheckBox(f"Add to every cell of {slide_name}" if slide_name else "Add to every cell of the slide", self)
        self._every.setChecked(bool(slide_name))
        self._every.setEnabled(bool(slide_name))
        self._save_list = QCheckBox("Save to the probe list", self)
        self._save_list.setChecked(True)
        self._sync_colour()

        form = QGridLayout()
        form.setHorizontalSpacing(10)
        form.setVerticalSpacing(9)
        rows = (("Name", self._name), ("Sequence", self._seq), ("", hint), ("Dye", self._dye),
                ("Colour", self._colour_btn), ("Channel", self._channel), ("", self._every), ("", self._save_list))
        for r, (label, widget) in enumerate(rows):
            if label:
                form.addWidget(QLabel(label, self), r, 0)
            form.addWidget(widget, r, 1, alignment=Qt.AlignmentFlag.AlignLeft)
        form.setColumnMinimumWidth(0, 110)
        cancel = QPushButton("Cancel", self)
        save = QPushButton("Save", self)
        save.setDefault(True)
        cancel.clicked.connect(self.reject)
        save.clicked.connect(self._on_save)
        buttons = QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(cancel)
        buttons.addWidget(save)
        root = QVBoxLayout(self)
        root.setContentsMargins(18, 16, 18, 16)
        root.setSpacing(14)
        root.addLayout(form)
        root.addLayout(buttons)
        if not channels:
            self._channel.addItem("No free channel", "")
            self._channel.setEnabled(False)

    # ------------------------------------------------------------------ values

    def name(self) -> str:
        return self._name.text().strip()

    def dye(self) -> str:
        d = self._dye.currentData()
        return "" if d in (None, ADD_DYE) else str(d)

    def colour(self) -> str:
        return self._colour or lists_mod.dye_color(self._lists, self.dye(), "#ffffff")

    def channel_key(self) -> str:
        return str(self._channel.currentData() or "")

    def sequence(self) -> str:
        p = self._seq.path()
        return str(p) if p else ""

    def every_cell(self) -> bool:
        return self._every.isChecked() and self._every.isEnabled()

    def save_to_list(self) -> bool:
        return self._save_list.isChecked()

    # ------------------------------------------------------------------ helpers

    def _fill_dyes(self, current: str) -> None:
        self._dye.blockSignals(True)
        self._dye.clear()
        for name in lists_mod.names(self._lists, "dyes"):
            self._dye.addItem(swatch_icon(lists_mod.dye_color(self._lists, name, "#ffffff")), name, name)
        self._dye.addItem(ADD_DYE, ADD_DYE)
        idx = self._dye.findData(current) if current else 0
        self._dye.setCurrentIndex(max(0, idx))
        self._dye.blockSignals(False)

    def _on_dye(self, _index: int) -> None:
        if self._dye.currentData() == ADD_DYE:
            name, ok = QInputDialog.getText(self, "Add New Dye", "Name of the dye:")
            name = name.strip() if ok else ""
            if name:
                lists = lists_mod.load_lists(self._folder)
                lists_mod.add_entry(lists, "dyes", name)
                lists_mod.save_lists(self._folder, lists)  # saved at once, as everywhere in the app
                self._lists = lists
            self._fill_dyes(name)
        self._colour = ""
        self._sync_colour()

    def _pick_colour(self) -> None:
        c = QColorDialog.getColor(QColor(self.colour()), self, "Colour")
        if c.isValid():
            self._colour = c.name()
            self._sync_colour()

    def _sync_colour(self) -> None:
        self._colour_btn.setIcon(swatch_icon(self.colour()))
        self._colour_btn.setText("Own colour" if self._colour else "From dye")

    def _on_save(self) -> None:
        if not self.name():
            QMessageBox.information(self, "New Probe", "Give the probe a name.")
            return
        if not self.channel_key():
            QMessageBox.information(self, "New Probe", "Every channel of this cell already has a probe.")
            return
        self.accept()
