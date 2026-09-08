from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.project import RepeatViewerProject
from repeatviewer.views.tool_dialogs import TeMultiInputFilesWidget, _run_cancel_button_box

from repeatviewer_minimap2.options import MINIMAP2_PRESETS, Minimap2Options

SETTINGS_GROUP = "plugin/minimap2"
_FASTA_FILTER = "FASTA (*.fa *.fasta *.fas *.fna);;All files (*)"


class Minimap2ParamsDialog(QDialog):
    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        project: RepeatViewerProject | None = None,
        initial_queries: list[Path] | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Minimap2")
        self.setMinimumWidth(520)
        self._settings = QSettings()
        self._project = project

        query_initial = list(initial_queries or [])
        self._query_input = TeMultiInputFilesWidget(
            self,
            initial=query_initial,
            start_empty=not query_initial,
            placeholder="Query FASTA to map",
            browse_title="Select query FASTA",
            file_filter=_FASTA_FILTER,
            add_button_label="Add another query FASTA",
            add_button_tooltip="Add another query sequence file. Each file is mapped separately.",
        )

        input_group = QGroupBox("Query input files", self)
        input_layout = QVBoxLayout(input_group)
        hint = QLabel(
            "FASTA files to map against the reference. Pre-filled from the Files panel when selected.",
            input_group,
        )
        hint.setWordWrap(True)
        input_layout.addWidget(hint)
        input_layout.addWidget(self._query_input)

        self._reference = QLineEdit(self)
        self._reference.setPlaceholderText("Reference genome FASTA")
        self._pan_combo: QComboBox | None = None
        if project is not None and project.is_pan_genome_mode and project.pan_genomes:
            self._pan_combo = QComboBox(self)
            for entry in project.pan_genomes:
                self._pan_combo.addItem(entry.name, str(entry.path))
            self._pan_combo.currentIndexChanged.connect(self._apply_pan_genome_selection)
            self._apply_pan_genome_selection()
        elif project is not None and project.genome_file is not None and project.genome_file.is_file():
            self._reference.setText(str(project.genome_file.expanduser().resolve()))

        ref_group = QGroupBox("Reference genome", self)
        ref_form = QFormLayout(ref_group)
        if self._pan_combo is not None:
            ref_form.addRow("Pan-genome:", self._pan_combo)
        ref_form.addRow("Reference FASTA:", self._reference_row())

        self._preset = QComboBox(self)
        for preset in MINIMAP2_PRESETS:
            self._preset.addItem(preset, preset)

        self._threads = QSpinBox(self)
        self._threads.setRange(1, 128)

        self._format = QComboBox(self)
        self._format.addItem("PAF", "paf")
        self._format.addItem("SAM", "sam")

        params_group = QGroupBox("Parameters", self)
        params_form = QFormLayout(params_group)
        params_form.addRow("Preset (-x):", self._preset)
        params_form.addRow("Threads:", self._threads)
        params_form.addRow("Output:", self._format)

        self._load_saved()

        buttons = _run_cancel_button_box(self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addWidget(input_group)
        layout.addWidget(ref_group)
        layout.addWidget(params_group)
        layout.addWidget(buttons)

    def _apply_pan_genome_selection(self) -> None:
        if self._pan_combo is None:
            return
        path = self._pan_combo.currentData()
        if path:
            self._reference.setText(str(path))

    def _reference_row(self) -> QWidget:
        row = QWidget(self)
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._reference)

        def browse() -> None:
            chosen, _ = QFileDialog.getOpenFileName(
                row,
                "Select reference FASTA",
                self._reference.text() or str(Path.home()),
                _FASTA_FILTER,
            )
            if chosen:
                self._reference.setText(chosen)

        browse_btn = QPushButton("Browse…", row)
        browse_btn.clicked.connect(browse)
        layout.addWidget(browse_btn)
        return row

    def _load_saved(self) -> None:
        self._settings.beginGroup(SETTINGS_GROUP)
        preset = str(self._settings.value("preset", "map-ont"))
        idx = self._preset.findData(preset)
        if idx >= 0:
            self._preset.setCurrentIndex(idx)
        self._threads.setValue(int(self._settings.value("threads", 4)))
        fmt = str(self._settings.value("output_format", "paf"))
        fmt_idx = self._format.findData(fmt)
        if fmt_idx >= 0:
            self._format.setCurrentIndex(fmt_idx)
        saved_ref = str(self._settings.value("reference_fasta", "")).strip()
        if saved_ref and not self._reference.text().strip():
            self._reference.setText(saved_ref)
        self._settings.endGroup()

    def _save(self) -> None:
        self._settings.beginGroup(SETTINGS_GROUP)
        self._settings.setValue("preset", str(self._preset.currentData()))
        self._settings.setValue("threads", int(self._threads.value()))
        self._settings.setValue("output_format", str(self._format.currentData()))
        self._settings.setValue("reference_fasta", self._reference.text().strip())
        self._settings.endGroup()

    def accept(self) -> None:
        queries = self.query_files()
        if not queries:
            QMessageBox.warning(
                self,
                "Minimap2",
                "Add at least one query FASTA file.",
            )
            return
        missing = [str(path) for path in queries if not path.is_file()]
        if missing:
            QMessageBox.warning(
                self,
                "Minimap2",
                "Query file not found:\n" + "\n".join(missing),
            )
            return
        reference = self._reference.text().strip()
        if not reference:
            QMessageBox.warning(
                self,
                "Minimap2",
                "Choose a reference genome FASTA.",
            )
            return
        ref_path = Path(reference).expanduser()
        if not ref_path.is_file():
            QMessageBox.warning(
                self,
                "Minimap2",
                f"Reference FASTA not found:\n{ref_path}",
            )
            return
        self._save()
        super().accept()

    def query_files(self) -> list[Path]:
        return self._query_input.paths()

    def options(self) -> Minimap2Options:
        return Minimap2Options(
            preset=str(self._preset.currentData()),
            threads=max(1, int(self._threads.value())),
            output_format=str(self._format.currentData()),
            reference_fasta=self._reference.text().strip(),
        )

