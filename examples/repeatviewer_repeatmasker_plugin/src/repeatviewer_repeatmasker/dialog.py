"""RepeatMasker parameter dialog for the example plugin."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

SETTINGS_GROUP = "plugin_repeatmasker_example"

# Common NCBI / RepeatMasker species examples (free text still allowed).
_SPECIES_EXAMPLES = (
    "human",
    "mouse",
    "rattus",
    "mammal",
    "drosophila",
    "arabidopsis",
    "rice",
    "maize",
    "danio",
    "worm",
)


class RepeatMaskerParamsDialog(QDialog):
    """Collect RepeatMasker options before the background job."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("RepeatMasker parameters")
        self.setMinimumWidth(560)

        s = QSettings()
        s.beginGroup(SETTINGS_GROUP)

        self._target_selected = QRadioButton("Mask selected FASTA file(s)", self)
        self._target_genome = QRadioButton(
            "Mask defined genome (selected FASTA = -lib)", self
        )
        if s.value("target", "selected") == "genome":
            self._target_genome.setChecked(True)
        else:
            self._target_selected.setChecked(True)

        self._lib_mode = QComboBox(self)
        self._lib_mode.addItem("Custom library (-lib)", "lib")
        self._lib_mode.addItem("Species / clade (-species)", "species")
        lib_mode = str(s.value("library_mode", "lib"))
        idx = self._lib_mode.findData(lib_mode)
        self._lib_mode.setCurrentIndex(idx if idx >= 0 else 0)

        self._library = QLineEdit(str(s.value("library_fasta", "")), self)
        self._library.setPlaceholderText("Custom repeat library FASTA")
        browse_lib = QPushButton("Browse…", self)
        browse_lib.clicked.connect(self._browse_library)
        lib_row = QWidget(self)
        lib_layout = QHBoxLayout(lib_row)
        lib_layout.setContentsMargins(0, 0, 0, 0)
        lib_layout.addWidget(self._library, 1)
        lib_layout.addWidget(browse_lib)

        self._species = QComboBox(self)
        self._species.setEditable(True)
        self._species.addItems(list(_SPECIES_EXAMPLES))
        self._species.setCurrentText(str(s.value("species", "human")))
        self._species.setToolTip(
            "NCBI Taxonomy name known to your RepeatMasker / FamDB install."
        )

        self._parallel = QSpinBox(self)
        self._parallel.setRange(1, 64)
        self._parallel.setValue(int(s.value("parallel", 4)))
        self._parallel.setToolTip(
            "RepeatMasker -pa: number of sequence batch jobs to run in parallel."
        )

        self._sensitivity = QComboBox(self)
        self._sensitivity.addItem("Default", "default")
        self._sensitivity.addItem("Sensitive (-s)", "sensitive")
        self._sensitivity.addItem("Quick (-q)", "quick")
        self._sensitivity.addItem("Rush (-qq)", "rush")
        sens = str(s.value("sensitivity", "default"))
        sidx = self._sensitivity.findData(sens)
        self._sensitivity.setCurrentIndex(sidx if sidx >= 0 else 0)

        self._engine = QComboBox(self)
        self._engine.addItem("Default engine", "")
        for eng in ("rmblast", "hmmer", "crossmatch", "abblast", "wublast"):
            self._engine.addItem(eng, eng)
        eng = str(s.value("engine", ""))
        eidx = self._engine.findData(eng)
        self._engine.setCurrentIndex(eidx if eidx >= 0 else 0)

        self._use_div = QCheckBox("Limit divergence (-div)", self)
        self._use_div.setChecked(s.value("use_div", True, type=bool))
        self._div = QSpinBox(self)
        self._div.setRange(1, 100)
        self._div.setSuffix(" %")
        self._div.setValue(int(s.value("div", 40)))
        self._div.setToolTip(
            "Only annotate copies less than this % diverged from the consensus."
        )

        self._cutoff = QSpinBox(self)
        self._cutoff.setRange(0, 1000)
        self._cutoff.setSpecialValueText("default (225)")
        self._cutoff.setValue(int(s.value("cutoff", 0)))
        self._cutoff.setToolTip("Score cutoff when using -lib (0 = RepeatMasker default).")

        self._gff = QCheckBox("Write GFF (-gff)", self)
        self._gff.setChecked(s.value("gff", True, type=bool))
        self._nolow = QCheckBox("Skip low-complexity / simple (-nolow)", self)
        self._nolow.setChecked(s.value("nolow", False, type=bool))
        self._norna = QCheckBox("Skip small RNA (-norna)", self)
        self._norna.setChecked(s.value("norna", True, type=bool))
        self._no_is = QCheckBox("Skip bacterial IS check (-no_is)", self)
        self._no_is.setChecked(s.value("no_is", True, type=bool))
        self._noint = QCheckBox("Only low-complexity / simple (-noint)", self)
        self._noint.setChecked(s.value("noint", False, type=bool))
        self._xsmall = QCheckBox("Soft-mask lowercase (-xsmall)", self)
        self._xsmall.setChecked(s.value("xsmall", False, type=bool))
        self._html = QCheckBox("HTML report (-html)", self)
        self._html.setChecked(s.value("html", False, type=bool))
        self._alignments = QCheckBox("Write alignments (-a)", self)
        self._alignments.setChecked(s.value("alignments", False, type=bool))
        self._uncurated = QCheckBox("Include uncurated FamDB (-uncurated)", self)
        self._uncurated.setChecked(s.value("uncurated", False, type=bool))

        s.endGroup()

        target_box = QGroupBox("Query", self)
        target_layout = QVBoxLayout(target_box)
        target_layout.addWidget(self._target_selected)
        target_layout.addWidget(self._target_genome)
        target_hint = QLabel(
            "Selected FASTA = sequences to mask, unless “defined genome” is chosen "
            "(then selection is the -lib library).",
            self,
        )
        target_hint.setWordWrap(True)
        target_hint.setStyleSheet("color: #666666;")
        target_layout.addWidget(target_hint)

        form = QFormLayout()
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(8)
        form.addRow("Repeat source:", self._lib_mode)
        form.addRow("Library FASTA:", lib_row)
        form.addRow("Species:", self._species)
        form.addRow("Parallel (-pa):", self._parallel)
        form.addRow("Sensitivity:", self._sensitivity)
        form.addRow("Engine (-e):", self._engine)
        div_row = QWidget(self)
        div_layout = QHBoxLayout(div_row)
        div_layout.setContentsMargins(0, 0, 0, 0)
        div_layout.addWidget(self._use_div)
        div_layout.addWidget(self._div)
        div_layout.addStretch(1)
        form.addRow("Divergence:", div_row)
        form.addRow("Lib cutoff (-cutoff):", self._cutoff)

        flags = QGroupBox("Output / filters", self)
        flags_layout = QVBoxLayout(flags)
        for w in (
            self._gff,
            self._nolow,
            self._norna,
            self._no_is,
            self._noint,
            self._xsmall,
            self._html,
            self._alignments,
            self._uncurated,
        ):
            flags_layout.addWidget(w)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addWidget(target_box)
        layout.addLayout(form)
        layout.addWidget(flags)
        layout.addWidget(buttons)

        self._lib_mode.currentIndexChanged.connect(self._sync_library_controls)
        self._target_selected.toggled.connect(self._sync_library_controls)
        self._target_genome.toggled.connect(self._sync_library_controls)
        self._use_div.toggled.connect(self._div.setEnabled)
        self._div.setEnabled(self._use_div.isChecked())
        self._sync_library_controls()

    def _browse_library(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select RepeatMasker library FASTA",
            self._library.text().strip() or "",
            "FASTA (*.fa *.fasta *.fas *.fna);;All files (*)",
        )
        if path:
            self._library.setText(path)

    def _sync_library_controls(self) -> None:
        genome_as_query = self._target_genome.isChecked()
        # When masking the genome, selection is always the -lib file.
        if genome_as_query:
            self._lib_mode.setCurrentIndex(self._lib_mode.findData("lib"))
            self._lib_mode.setEnabled(False)
            self._library.setEnabled(False)
            self._species.setEnabled(False)
        else:
            self._lib_mode.setEnabled(True)
            use_lib = self._lib_mode.currentData() == "lib"
            self._library.setEnabled(use_lib)
            self._species.setEnabled(not use_lib)
        self._cutoff.setEnabled(
            genome_as_query or self._lib_mode.currentData() == "lib"
        )

    def options(self) -> dict:
        genome_target = self._target_genome.isChecked()
        library_mode = "lib" if genome_target else str(self._lib_mode.currentData())
        return {
            "target": "genome" if genome_target else "selected",
            "library_mode": library_mode,
            "library_fasta": self._library.text().strip(),
            "species": self._species.currentText().strip(),
            "parallel": int(self._parallel.value()),
            "sensitivity": str(self._sensitivity.currentData()),
            "engine": str(self._engine.currentData() or ""),
            "use_div": self._use_div.isChecked(),
            "div": int(self._div.value()),
            "cutoff": int(self._cutoff.value()),
            "gff": self._gff.isChecked(),
            "nolow": self._nolow.isChecked(),
            "norna": self._norna.isChecked(),
            "no_is": self._no_is.isChecked(),
            "noint": self._noint.isChecked(),
            "xsmall": self._xsmall.isChecked(),
            "html": self._html.isChecked(),
            "alignments": self._alignments.isChecked(),
            "uncurated": self._uncurated.isChecked(),
        }

    def accept(self) -> None:
        opts = self.options()
        if opts["target"] == "selected" and opts["library_mode"] == "species":
            if not opts["species"]:
                from PySide6.QtWidgets import QMessageBox

                QMessageBox.warning(self, "RepeatMasker", "Enter a -species name.")
                return
        if opts["target"] == "selected" and opts["library_mode"] == "lib":
            lib = Path(opts["library_fasta"]).expanduser() if opts["library_fasta"] else None
            if lib is None or not lib.is_file():
                from PySide6.QtWidgets import QMessageBox

                QMessageBox.warning(
                    self,
                    "RepeatMasker",
                    "Choose a custom library FASTA for -lib, or switch to -species.",
                )
                return

        s = QSettings()
        s.beginGroup(SETTINGS_GROUP)
        for key, value in opts.items():
            s.setValue(key, value)
        s.endGroup()
        super().accept()
