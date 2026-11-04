"""Edit slide: name, number, dates, species and the method of one slide. Its images follow unless set by hand."""

from __future__ import annotations

import copy
from pathlib import Path

from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.cytogenetics import lists as ls
from repeatviewer.backends.cytogenetics.model import Slide, update_slide
from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.cytogenetics_style import form_layout, section_label
from repeatviewer.widgets.help_button import make_help_button

HELP = (
    "<b>Slide</b>: the details of one microscope slide. Every image on the slide takes these values "
    "(species, accession, technique, counterstain, date, microscope, experimenter) unless you gave the image "
    "its own value in the Image tab.<br>"
    "<b>Assign images</b>: select images in the Library (Cmd-click or Shift-click for several), then right-click "
    "and choose Assign to slide, or drag them onto a slide heading."
)


class SlideDialog(QDialog):
    def __init__(self, folder: Path, slide: Slide, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Edit slide")
        self._folder = Path(folder)
        self._slide = copy.deepcopy(slide)
        lists = ls.load_lists(self._folder)

        def combo(key: str, value: str) -> NoWheelComboBox:
            box = NoWheelComboBox(self)
            box.setEditable(True)
            box.addItem("")
            box.addItems(ls.names(lists, key))
            box.setCurrentText(value)
            box.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            return box

        s = self._slide
        self._f = {
            "name": QLineEdit(s.name, self), "box": QLineEdit(s.box, self), "date": QLineEdit(s.date, self),
            "species": combo("species", s.species), "accession": QLineEdit(s.accession, self),
            "technique": combo("techniques", s.technique), "counterstain": combo("counterstains", s.counterstain),
            "probes": QLineEdit(s.probes, self), "hyb_temp": QLineEdit(s.hyb_temp, self),
            "fixation_date": QLineEdit(s.fixation_date, self), "microscope": QLineEdit(s.microscope, self),
            "experimenter": QLineEdit(s.experimenter, self), "protocol": QLineEdit(s.protocol, self),
        }
        for key in ("date", "fixation_date"):
            self._f[key].setPlaceholderText("YYYY-MM-DD")
        self._f["hyb_temp"].setPlaceholderText("e.g. 37 °C overnight")
        self._f["probes"].setPlaceholderText("e.g. CroSat2 (FITC), 45S rDNA (Cy3)")
        self._note = QPlainTextEdit(s.note, self)
        self._note.setFixedHeight(64)

        form = form_layout()
        head = QHBoxLayout()
        head.addWidget(section_label("Slide", self))
        head.addStretch(1)
        head.addWidget(make_help_button(self, "Edit slide", HELP))
        form.addRow(head)
        for label, key in (("Name", "name"), ("Slide / box number", "box"), ("Date", "date"), ("Species", "species"),
                           ("Accession / line", "accession")):
            form.addRow(label, self._f[key])
        form.addRow(section_label("Method", self))
        for label, key in (("Technique", "technique"), ("Counterstain", "counterstain"),
                           ("Probes / hybridisation", "probes"), ("Hybridisation temperature", "hyb_temp"),
                           ("Fixation date", "fixation_date"), ("Microscope", "microscope"),
                           ("Experimenter", "experimenter"), ("Protocol", "protocol")):
            form.addRow(label, self._f[key])
        form.addRow("Note", self._note)

        cancel = QPushButton("Cancel", self)
        save = QPushButton("Save", self)
        save.setDefault(True)
        cancel.clicked.connect(self.reject)
        save.clicked.connect(self._save)
        buttons = QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(cancel)
        buttons.addWidget(save)
        root = QVBoxLayout(self)
        root.setContentsMargins(22, 18, 22, 16)
        root.setSpacing(10)
        root.addLayout(form)
        root.addLayout(buttons)
        self.resize(520, self.sizeHint().height())

    def _save(self) -> None:
        s = self._slide
        for key, w in self._f.items():
            setattr(s, key, (w.currentText() if isinstance(w, NoWheelComboBox) else w.text()).strip())
        s.note = self._note.toPlainText().strip()
        if not s.name:
            QMessageBox.information(self, "Edit slide", "Give the slide a name.")
            return
        try:
            update_slide(self._folder, s)
        except OSError as exc:
            QMessageBox.warning(self, "Edit slide", f"Could not save:\n{exc}")
            return
        from repeatviewer.widgets.cyto_images_ui import notify_changed

        notify_changed(self._folder)
        self.accept()
