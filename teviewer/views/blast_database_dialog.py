"""Manage saved BLAST databases (Blast tool's "Database:" dropdown, Geneious-style)."""

from __future__ import annotations

import shutil
from datetime import datetime
from pathlib import Path

from PySide6.QtCore import QObject, Qt, QThread, Signal, Slot
from PySide6.QtWidgets import (
    QAbstractItemView,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QLineEdit,
    QMessageBox,
    QProgressDialog,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from teviewer.app_preferences import load_blast_db_base_dir
from teviewer.backends.database_prepare import (
    blast_database_is_ready,
    build_named_blast_database,
    format_prepare_error,
)
from teviewer.backends.gui_functions import sanitize_name
from teviewer.database_settings import save_blast_databases_to_settings
from teviewer.project import (
    BLAST_MOLECULE_NUCL,
    BLAST_MOLECULE_PROT,
    BlastDatabase,
    TEviewerProject,
)
from teviewer.views.tool_dialogs import exec_modal_dialog

_MOLECULE_LABELS = {
    BLAST_MOLECULE_NUCL: "Nucleotide",
    BLAST_MOLECULE_PROT: "Protein (not usable by Blast yet)",
}


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
    """Pick a FASTA + name, build the BLAST index in the background."""

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
        self._existing_names = {n.strip().lower() for n in existing_names}
        self._use_system_blast = use_system_blast
        self._result_db: BlastDatabase | None = None
        self._thread: QThread | None = None
        self._worker: _BuildDatabaseWorker | None = None
        self._progress: QProgressDialog | None = None
        self._pending_name = ""
        self._pending_fasta = Path()
        self._pending_molecule_type = BLAST_MOLECULE_NUCL

        self._name = QLineEdit(self)
        self._name.setPlaceholderText("Database name")
        self._fasta = QLineEdit(self)
        self._fasta.setPlaceholderText("FASTA file")
        browse_btn = QPushButton("Browse", self)
        browse_btn.clicked.connect(self._browse)
        fasta_row = QWidget(self)
        fasta_h = QHBoxLayout(fasta_row)
        fasta_h.setContentsMargins(0, 0, 0, 0)
        fasta_h.addWidget(self._fasta, 1)
        fasta_h.addWidget(browse_btn)

        # The BLASTN pipeline can't search protein databases yet (see
        # BlastParamsDialog's dropdown filter), so the quick "Add as BLAST
        # database…" entry point from a FASTA's right-click menu skips the
        # type picker entirely and always builds nucleotide.
        self._molecule_type = QComboBox(self)
        for value, label in _MOLECULE_LABELS.items():
            self._molecule_type.addItem(label, value)

        if initial_fasta is not None:
            self._fasta.setText(str(initial_fasta))
            self._name.setText(initial_fasta.stem)

        form = QFormLayout()
        form.addRow("Name:", self._name)
        form.addRow("FASTA:", fasta_row)
        if allow_molecule_type:
            form.addRow("Type:", self._molecule_type)

        self._buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self._buttons.accepted.connect(self._on_ok)
        self._buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self._buttons)

    def _browse(self) -> None:
        chosen, _ = QFileDialog.getOpenFileName(
            self,
            "Choose FASTA file",
            str(Path.home()),
            "FASTA (*.fa *.fasta *.fna *.faa);;All files (*)",
        )
        if chosen:
            self._fasta.setText(chosen)
            if not self._name.text().strip():
                self._name.setText(Path(chosen).stem)

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

        molecule_type = str(self._molecule_type.currentData())
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
        progress.show()
        self._progress = progress

        self._pending_name = name
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
        )
        self.accept()

    @Slot(str)
    def _on_build_failed(self, message: str) -> None:
        self._cleanup_progress()
        QMessageBox.critical(self, "Add BLAST database", message)


class BlastDatabaseManagerDialog(QDialog):
    """File → manage saved BLAST databases (Blast tool's Database dropdown)."""

    def __init__(self, project: TEviewerProject, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Manage BLAST databases")
        self.setMinimumWidth(640)
        self.setMinimumHeight(360)
        self._project = project

        self._table = QTableWidget(0, 5, self)
        self._table.setHorizontalHeaderLabels(
            ["Name", "Type", "Source FASTA", "Created", "Status"]
        )
        self._table.horizontalHeader().setSectionResizeMode(2, QHeaderView.ResizeMode.Stretch)
        self._table.verticalHeader().setVisible(False)
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.itemSelectionChanged.connect(self._update_buttons)

        self._add_btn = QPushButton("Add Database…", self)
        self._remove_btn = QPushButton("Remove", self)
        self._add_btn.clicked.connect(self._add_database)
        self._remove_btn.clicked.connect(self._remove_database)

        btn_row = QHBoxLayout()
        btn_row.addWidget(self._add_btn)
        btn_row.addWidget(self._remove_btn)
        btn_row.addStretch(1)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.reject)
        close_btn = buttons.button(QDialogButtonBox.StandardButton.Close)
        if close_btn is not None:
            close_btn.clicked.connect(self.accept)

        layout = QVBoxLayout(self)
        layout.addWidget(self._table, 1)
        layout.addLayout(btn_row)
        layout.addWidget(buttons)

        self._populate()
        self._update_buttons()

    def _populate(self) -> None:
        self._table.setRowCount(0)
        for db in self._project.blast_databases:
            row = self._table.rowCount()
            self._table.insertRow(row)
            ready = blast_database_is_ready(db)
            self._table.setItem(row, 0, QTableWidgetItem(db.name))
            self._table.setItem(
                row, 1, QTableWidgetItem(_MOLECULE_LABELS.get(db.molecule_type, db.molecule_type))
            )
            self._table.setItem(row, 2, QTableWidgetItem(str(db.source_fasta)))
            self._table.setItem(row, 3, QTableWidgetItem(db.created_at))
            status_item = QTableWidgetItem("Ready" if ready else "Missing index files")
            if not ready:
                status_item.setForeground(Qt.GlobalColor.red)
            self._table.setItem(row, 4, status_item)

    def _selected_index(self) -> int | None:
        rows = {idx.row() for idx in self._table.selectedIndexes()}
        if not rows:
            return None
        return next(iter(rows))

    def _update_buttons(self) -> None:
        self._remove_btn.setEnabled(self._selected_index() is not None)

    def _persist(self) -> None:
        save_blast_databases_to_settings(self._project.blast_databases)

    def _add_database(self) -> None:
        existing_names = {db.name for db in self._project.blast_databases}
        dialog = AddBlastDatabaseDialog(
            existing_names=existing_names,
            use_system_blast=self._project.use_system_blast,
            parent=self,
        )
        if not exec_modal_dialog(dialog):
            return
        db = dialog.result_database()
        if db is None:
            return
        self._project.blast_databases.append(db)
        self._persist()
        self._populate()
        self._update_buttons()

    def _remove_database(self) -> None:
        index = self._selected_index()
        if index is None:
            return
        db = self._project.blast_databases[index]
        answer = QMessageBox.question(
            self,
            "Manage BLAST databases",
            f"Remove “{db.name}” and delete its index files?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        shutil.rmtree(db.db_prefix_path.parent, ignore_errors=True)
        del self._project.blast_databases[index]
        self._persist()
        self._populate()
        self._update_buttons()
