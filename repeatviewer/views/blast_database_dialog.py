"""Manage saved BLAST databases (Blast tool's "Database:" dropdown, Geneious-style).

The list UI now lives in :mod:`manage_references_dialog`. This module keeps
``AddBlastDatabaseDialog`` (used when auto-saving a genome index) and re-exports
``BlastDatabaseManagerDialog`` for existing call sites.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, QTimer, Signal, Slot
from PySide6.QtWidgets import (
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QLineEdit,
    QMessageBox,
    QProgressBar,
    QProgressDialog,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.app_preferences import load_blast_db_base_dir
from repeatviewer.backends.database_prepare import (
    build_named_blast_database,
    format_prepare_error,
)
from repeatviewer.backends.gui_functions import sanitize_name
from repeatviewer.project import (
    ANNOTATION_NAME_MAX_LENGTH,
    BLAST_MOLECULE_NUCL,
    BLAST_MOLECULE_PROT,
    BlastDatabase,
    normalize_annotation_name,
)
from repeatviewer.widgets.help_button import field_with_help
from repeatviewer.widgets.path_browse import PathBrowseRow

_SHORT_NAME_HELP = (
    "Optional short label (max "
    f"{ANNOTATION_NAME_MAX_LENGTH} characters) shown in BLAST database lists "
    "and on the Annotation tab. Leave blank to use the full Name."
)


def style_blast_db_progress(progress: QProgressDialog) -> None:
    """Make the makeblastdb busy bar wider (3/2) and shorter (1/2)."""
    progress.ensurePolished()
    base_w = max(progress.sizeHint().width(), 240)
    progress.setMinimumWidth(int(base_w * 3 / 2))
    bar = progress.findChild(QProgressBar)
    if bar is not None:
        hint_h = bar.sizeHint().height()
        if hint_h <= 0:
            hint_h = 16
        bar.setFixedHeight(max(6, hint_h // 2))
        hint_w = max(bar.sizeHint().width(), bar.minimumWidth(), 160)
        bar.setMinimumWidth(int(hint_w * 3 / 2))


# Re-export for callers that import from this module.
from repeatviewer.views.manage_references_dialog import (  # noqa: E402
    BlastDatabaseManagerDialog,
    ManageReferencesDialog,
)

_MOLECULE_LABELS = {
    BLAST_MOLECULE_NUCL: "Nucleotide",
    BLAST_MOLECULE_PROT: "Protein",
}

_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".faa"}


class _BuildDatabaseWorker(QObject):
    finished = Signal(object)
    failed = Signal(str)

    def __init__(
        self,
        source_fasta: Path,
        base_dir: Path,
        name: str,
        *,
        molecule_type: str,
        use_system_blast: bool,
    ) -> None:
        super().__init__()
        self._source_fasta = source_fasta
        self._base_dir = base_dir
        self._name = name
        self._molecule_type = molecule_type
        self._use_system_blast = use_system_blast

    @Slot()
    def run(self) -> None:
        try:
            db_prefix = build_named_blast_database(
                self._source_fasta,
                self._base_dir,
                self._name,
                molecule_type=self._molecule_type,
                use_system_blast=self._use_system_blast,
            )
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(format_prepare_error(exc))
            return
        self.finished.emit(db_prefix)


class AddBlastDatabaseDialog(QDialog):
    """Pick a FASTA + name, build the BLAST index in the background.

    Non-modal so a click in the Files panel can fill the FASTA path (same
    pattern as Add annotation / Define database).
    """

    def __init__(
        self,
        *,
        existing_names: set[str],
        use_system_blast: bool,
        initial_fasta: Path | None = None,
        allow_molecule_type: bool = True,
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Add BLAST database")
        self.setMinimumWidth(480)
        self.setModal(False)
        self.setWindowModality(Qt.WindowModality.NonModal)
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, True)
        self._existing_names = {n.strip().lower() for n in existing_names}
        self._use_system_blast = use_system_blast
        self._result_db: BlastDatabase | None = None
        self._thread: QThread | None = None
        self._worker: _BuildDatabaseWorker | None = None
        self._progress: QProgressDialog | None = None
        self._pending_name = ""
        self._pending_short_name = ""
        self._pending_fasta = Path()
        self._pending_molecule_type = BLAST_MOLECULE_NUCL

        self._name = QLineEdit(self)
        self._name.setPlaceholderText("Database name")
        self._short_name = QLineEdit(self)
        self._short_name.setMaxLength(ANNOTATION_NAME_MAX_LENGTH)
        self._short_name.setPlaceholderText(
            f"Optional — short label (max {ANNOTATION_NAME_MAX_LENGTH})"
        )
        self._fasta_row = PathBrowseRow(
            self,
            placeholder="Click a FASTA in Files, or Browse",
            browse_title="Choose FASTA file",
            file_filter="FASTA (*.fa *.fasta *.fna *.fas *.faa);;All files (*)",
            mode="file",
            browse_handler=self._browse,
        )
        self._fasta = self._fasta_row.line_edit()
        self._fasta_row.set_default_tooltip(
            "Click a FASTA file in the Files panel to fill this path, or use Browse."
        )

        # blastn uses nucleotide databases; blastp uses protein databases.
        self._molecule_type: QComboBox | None = None
        if allow_molecule_type:
            combo = QComboBox(self)
            for value, label in _MOLECULE_LABELS.items():
                combo.addItem(label, value)
            self._molecule_type = combo

        if initial_fasta is not None:
            self.apply_path_from_panel(initial_fasta)

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setSpacing(8)
        if self._molecule_type is not None:
            form.addRow("Type:", self._molecule_type)
        form.addRow("Name:", self._name)
        form.addRow(
            "Short name:",
            field_with_help(
                self._short_name,
                "Short name",
                _SHORT_NAME_HELP,
                parent=self,
            ),
        )
        form.addRow("FASTA:", self._fasta_row)

        self._buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self._buttons.accepted.connect(self._on_ok)
        self._buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(12)
        layout.addLayout(form)
        layout.addWidget(self._buttons)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        # Prefer FASTA for Files-panel / Browse fill; defer so layout is ready.
        QTimer.singleShot(0, self._focus_fasta_field)

    def _focus_fasta_field(self) -> None:
        if self.isVisible():
            self._fasta_row.focus_field()

    def _keep_in_front(self) -> None:
        if not self.isVisible():
            return
        self.raise_()
        self.activateWindow()

    def apply_path_from_panel(self, path: Path) -> bool:
        """Fill FASTA (and default Name) from a Files-panel click."""
        resolved = path.expanduser().resolve()
        if not resolved.is_file():
            return False
        if resolved.suffix.lower() not in _FASTA_SUFFIXES:
            return False
        current = self._fasta.text().strip()
        if current and Path(current).expanduser().resolve() == resolved:
            self._keep_in_front()
            QTimer.singleShot(0, self._keep_in_front)
            return True
        self._fasta.setText(str(resolved))
        if not self._name.text().strip():
            self._name.setText(resolved.stem)
        # Protein FASTA often uses .faa — switch type when obvious.
        if (
            self._molecule_type is not None
            and resolved.suffix.lower() == ".faa"
        ):
            index = self._molecule_type.findData(BLAST_MOLECULE_PROT)
            if index >= 0:
                self._molecule_type.setCurrentIndex(index)
        self._keep_in_front()
        QTimer.singleShot(0, self._keep_in_front)
        return True

    def _browse(self) -> None:
        start = self._fasta.text().strip()
        start_dir = str(Path(start).expanduser().parent) if start else str(Path.home())
        chosen, _ = QFileDialog.getOpenFileName(
            self,
            "Choose FASTA file",
            start_dir,
            "FASTA (*.fa *.fasta *.fna *.fas *.faa);;All files (*)",
        )
        if chosen:
            self.apply_path_from_panel(Path(chosen))

    def result_database(self) -> BlastDatabase | None:
        return self._result_db

    def _on_ok(self) -> None:
        name = self._name.text().strip()
        fasta_text = self._fasta.text().strip()
        if not name:
            QMessageBox.warning(self, "Add BLAST database", "Enter a database name.")
            return
        if name.lower() in self._existing_names:
            QMessageBox.warning(
                self, "Add BLAST database", f"A database named “{name}” already exists."
            )
            return
        if not fasta_text:
            QMessageBox.warning(self, "Add BLAST database", "Choose a FASTA file.")
            return
        fasta_path = Path(fasta_text).expanduser()
        if not fasta_path.is_file():
            QMessageBox.warning(self, "Add BLAST database", f"FASTA not found:\n{fasta_path}")
            return

        molecule_type = (
            str(self._molecule_type.currentData())
            if self._molecule_type is not None
            else BLAST_MOLECULE_NUCL
        )
        base_dir = load_blast_db_base_dir()
        safe_name = sanitize_name(name) or "database"
        target_dir = base_dir.expanduser().resolve() / safe_name
        if target_dir.exists():
            QMessageBox.warning(
                self,
                "Add BLAST database",
                f"A database folder already exists at:\n{target_dir}\n\n"
                "Choose a different name.",
            )
            return

        self._buttons.setEnabled(False)
        progress = QProgressDialog("Building BLAST database…", "", 0, 0, self)
        progress.setWindowTitle("Add BLAST database")
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.setCancelButton(None)
        progress.setMinimumDuration(0)
        progress.setAutoClose(False)
        progress.setAutoReset(False)
        style_blast_db_progress(progress)
        progress.show()
        self._progress = progress

        short_name = normalize_annotation_name(self._short_name.text())
        self._pending_name = name
        self._pending_short_name = short_name
        self._pending_fasta = fasta_path
        self._pending_molecule_type = molecule_type

        worker = _BuildDatabaseWorker(
            fasta_path,
            base_dir,
            name,
            molecule_type=molecule_type,
            use_system_blast=self._use_system_blast,
        )
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(self._on_build_finished)
        worker.failed.connect(self._on_build_failed)
        worker.finished.connect(thread.quit)
        worker.failed.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        self._thread = thread
        self._worker = worker
        thread.start()

    def _cleanup_progress(self) -> None:
        if self._progress is not None:
            self._progress.close()
            self._progress = None
        self._buttons.setEnabled(True)
        self._thread = None
        self._worker = None

    @Slot(object)
    def _on_build_finished(self, db_prefix: object) -> None:
        self._cleanup_progress()
        self._result_db = BlastDatabase(
            name=self._pending_name,
            source_fasta=self._pending_fasta,
            db_prefix_path=Path(str(db_prefix)),
            molecule_type=self._pending_molecule_type,
            created_at=datetime.now().isoformat(timespec="seconds"),
            annotation_name=self._pending_short_name,
        )
        self.accept()

    @Slot(str)
    def _on_build_failed(self, message: str) -> None:
        self._cleanup_progress()
        QMessageBox.critical(self, "Add BLAST database", message)


__all__ = [
    "AddBlastDatabaseDialog",
    "BlastDatabaseManagerDialog",
    "ManageReferencesDialog",
]
