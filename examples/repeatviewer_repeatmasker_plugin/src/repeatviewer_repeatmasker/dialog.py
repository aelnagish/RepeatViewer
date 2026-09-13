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
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.views.tool_dialogs import make_process_output_box

SETTINGS_GROUP = "plugin_repeatmasker_example"

# Common NCBI / RepeatMasker species examples (free text still allowed).
_SPECIES_EXAMPLES = (
    "other",
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

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        selected_files: list[Path] | None = None,
        genome_file: Path | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("RepeatMasker parameters")
        self.setMinimumWidth(560)

        selected = [p for p in (selected_files or []) if p is not None]
        default_lib = str(selected[0].expanduser().resolve()) if selected else ""
        default_genome = ""
        if genome_file is not None and Path(genome_file).is_file():
            default_genome = str(Path(genome_file).expanduser().resolve())

        s = QSettings()
        s.beginGroup(SETTINGS_GROUP)

        self._lib_mode = QComboBox(self)
        self._lib_mode.addItem("Custom library (-lib)", "lib")
        self._lib_mode.addItem("Species / clade (-species)", "species")
        lib_mode = str(s.value("library_mode", "lib"))
        idx = self._lib_mode.findData(lib_mode)
        self._lib_mode.setCurrentIndex(idx if idx >= 0 else 0)

        # TE library is the selected FASTA by default (used as -lib).
        saved_lib = str(s.value("library_fasta", "") or "").strip()
        self._library = QLineEdit(default_lib or saved_lib, self)
        self._library.setPlaceholderText("Selected TE consensus FASTA (-lib)")
        browse_lib = QPushButton("Browse", self)
        browse_lib.clicked.connect(self._browse_library)
        lib_row = QWidget(self)
        lib_layout = QHBoxLayout(lib_row)
        lib_layout.setContentsMargins(0, 0, 0, 0)
        lib_layout.addWidget(self._library, 1)
        lib_layout.addWidget(browse_lib)

        saved_genome = str(s.value("genome_fasta", "") or "").strip()
        self._genome = QLineEdit(default_genome or saved_genome, self)
        self._genome.setPlaceholderText("Genome FASTA to mask")
        browse_genome = QPushButton("Browse", self)
        browse_genome.clicked.connect(self._browse_genome)
        genome_row = QWidget(self)
        genome_layout = QHBoxLayout(genome_row)
        genome_layout.setContentsMargins(0, 0, 0, 0)
        genome_layout.addWidget(self._genome, 1)
        genome_layout.addWidget(browse_genome)

        self._species = QComboBox(self)
        self._species.setEditable(True)
        self._species.addItems(list(_SPECIES_EXAMPLES))
        self._species.setCurrentText(str(s.value("species", "other")))
        self._species.setToolTip(
            "NCBI Taxonomy name known to your RepeatMasker / FamDB install. "
            "Shown as “other” when using -lib."
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
        sens = str(s.value("sensitivity", "sensitive"))
        sidx = self._sensitivity.findData(sens)
        self._sensitivity.setCurrentIndex(sidx if sidx >= 0 else 0)

        self._engine = QComboBox(self)
        self._engine.addItem("Default engine", "")
        for eng in ("rmblast", "hmmer", "crossmatch", "abblast", "wublast"):
            self._engine.addItem(eng, eng)
        eng = str(s.value("engine", "rmblast"))
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

        process_box, self._show_process_window = make_process_output_box(
            self,
            tool_name="RepeatMasker",
            checked=s.value("show_process_window", False, type=bool),
        )

        s.endGroup()

        form = QFormLayout()
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(8)
        form.addRow("Repeat source:", self._lib_mode)
        form.addRow("TE library (-lib):", lib_row)
        form.addRow("Genome:", genome_row)
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

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Cancel)
        self._save_btn = buttons.addButton("Save", QDialogButtonBox.ButtonRole.ActionRole)
        self._run_btn = buttons.addButton("Run", QDialogButtonBox.ButtonRole.AcceptRole)
        self._save_btn.clicked.connect(self._save_settings)
        self._run_btn.clicked.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(flags)
        layout.addWidget(process_box)
        layout.addWidget(buttons)

        self._lib_mode.currentIndexChanged.connect(self._sync_library_controls)
        self._use_div.toggled.connect(self._div.setEnabled)
        self._div.setEnabled(self._use_div.isChecked())
        self._sync_library_controls()

    def _browse_library(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select TE library FASTA (-lib)",
            self._library.text().strip() or "",
            "FASTA (*.fa *.fasta *.fas *.fna);;All files (*)",
        )
        if path:
            self._library.setText(path)

    def _browse_genome(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Select genome FASTA",
            self._genome.text().strip() or "",
            "FASTA (*.fa *.fasta *.fas *.fna);;All files (*)",
        )
        if path:
            self._genome.setText(path)

    def _sync_library_controls(self) -> None:
        use_lib = self._lib_mode.currentData() == "lib"
        self._library.setEnabled(use_lib)
        self._cutoff.setEnabled(use_lib)
        self._species.setEnabled(not use_lib)
        if use_lib:
            self._species.setCurrentText("other")

    def options(self) -> dict:
        return {
            "target": "genome",
            "library_mode": str(self._lib_mode.currentData()),
            "library_fasta": self._library.text().strip(),
            "genome_fasta": self._genome.text().strip(),
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
            "show_process_window": self._show_process_window.isChecked(),
        }

    def _validate(self) -> bool:
        opts = self.options()
        genome = Path(opts["genome_fasta"]).expanduser() if opts["genome_fasta"] else None
        if genome is None or not genome.is_file():
            QMessageBox.warning(
                self,
                "RepeatMasker",
                "Choose a genome FASTA to mask.",
            )
            return False
        if opts["library_mode"] == "species":
            if not opts["species"] or opts["species"].lower() == "other":
                QMessageBox.warning(self, "RepeatMasker", "Enter a -species name.")
                return False
        else:
            lib = Path(opts["library_fasta"]).expanduser() if opts["library_fasta"] else None
            if lib is None or not lib.is_file():
                QMessageBox.warning(
                    self,
                    "RepeatMasker",
                    "Choose a TE library FASTA for -lib (usually the selected file).",
                )
                return False
        return True

    def _persist(self, opts: dict) -> None:
        s = QSettings()
        s.beginGroup(SETTINGS_GROUP)
        for key, value in opts.items():
            s.setValue(key, value)
        s.endGroup()

    def _save_settings(self) -> None:
        if not self._validate():
            return
        self._persist(self.options())
        QMessageBox.information(self, "RepeatMasker", "Settings saved.")

    def accept(self) -> None:
        if not self._validate():
            return
        self._persist(self.options())
        super().accept()
