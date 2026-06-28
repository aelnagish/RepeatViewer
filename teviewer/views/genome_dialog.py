"""Choose or change the project genome FASTA (File → Define genome)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from teviewer.project import TEviewerProject
from teviewer.views.settings_dialog import (
    MAX_GENOME_HISTORY,
    _genome_history_label,
    add_genome_to_history,
    load_genome_history,
    remove_genome_from_history,
    save_workspace,
)


class GenomeDialog(QDialog):
    """Select the genome FASTA used for Blast, Extend, and TE genome map."""

    def __init__(self, project: TEviewerProject, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Define genome")
        self.setMinimumWidth(640)
        self._project = project

        self._genome = QComboBox()
        self._genome.setEditable(True)
        self._genome.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self._genome.setMinimumWidth(360)
        self._genome.lineEdit().setPlaceholderText("Select a recent genome or browse")
        self._genome.lineEdit().setClearButtonEnabled(True)
        self._genome.setToolTip(
            "Genome FASTA for Blast, Extend, and TE genome map. "
            f"Up to {MAX_GENOME_HISTORY} recent paths are remembered."
        )
        self._populate_genome_history(project.genome_file)
        self._genome.activated.connect(self._on_genome_history_selected)

        form = QFormLayout()
        form.addRow("Genome FASTA:", self._genome_row())

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        hint = QLabel(
            "Chromosome sizes and BLAST index are prepared automatically after you choose a genome."
        )
        hint.setWordWrap(True)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(hint)
        layout.addWidget(buttons)

    def _path_str(self, path: Path | None) -> str:
        return str(path) if path else ""

    def _populate_genome_history(self, current: Path | None) -> None:
        history = load_genome_history()
        if current is not None:
            try:
                resolved = current.expanduser().resolve()
            except OSError:
                resolved = None
            if resolved is not None and resolved.is_file() and resolved not in history:
                history.insert(0, resolved)
        self._genome.blockSignals(True)
        self._genome.clear()
        for path in history[:MAX_GENOME_HISTORY]:
            self._genome.addItem(_genome_history_label(path), str(path))
        self._genome.blockSignals(False)
        if current is not None:
            self._genome.lineEdit().setText(self._path_str(current))
        elif self._genome.count() > 0:
            self._genome.setCurrentIndex(0)
            self._on_genome_history_selected(0)

    def _on_genome_history_selected(self, index: int) -> None:
        if index < 0:
            return
        path = self._genome.itemData(index)
        if path:
            self._genome.lineEdit().setText(str(path))

    def _current_genome_path(self) -> Path | None:
        text = self._genome.lineEdit().text().strip()
        if not text:
            return None
        return Path(text).expanduser()

    def _genome_row(self) -> QWidget:
        container = QWidget()
        h = QHBoxLayout(container)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(self._genome, stretch=1)

        def browse() -> None:
            start = self._genome.lineEdit().text().strip() or str(Path.home())
            chosen, _ = QFileDialog.getOpenFileName(
                self,
                "Select genome FASTA",
                start,
                "FASTA (*.fa *.fasta *.fna);;All files (*)",
            )
            if chosen:
                self._genome.lineEdit().setText(chosen)

        browse_btn = QPushButton("Browse")
        browse_btn.clicked.connect(browse)
        h.addWidget(browse_btn)

        remove_btn = QPushButton("Remove")
        remove_btn.setToolTip("Remove the current path from recent genomes (file is not deleted).")
        remove_btn.clicked.connect(self._remove_genome_history_entry)
        h.addWidget(remove_btn)
        return container

    def _remove_genome_history_entry(self) -> None:
        path = self._current_genome_path()
        if path is None:
            return
        keep_text = self._genome.lineEdit().text()
        remove_genome_from_history(path)
        self._populate_genome_history(None)
        self._genome.lineEdit().setText(keep_text)

    def apply_to_project(self) -> TEviewerProject:
        genome = self._current_genome_path()
        self._project.genome_file = genome
        if genome is not None and genome.is_file():
            add_genome_to_history(genome)
        save_workspace(self._project)
        return self._project
