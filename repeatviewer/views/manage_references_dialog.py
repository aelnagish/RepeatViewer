"""Manage Database — BLAST, HMM, and Other databases in one window."""

from __future__ import annotations

import logging
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Literal

from PySide6.QtCore import QObject, Qt, QThread, QTimer, Signal, Slot
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFormLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QProgressDialog,
    QPushButton,
    QStyledItemDelegate,
    QStyleOptionViewItem,
    QTabBar,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.app_preferences import load_blast_db_base_dir
from repeatviewer.backends.database_prepare import (
    blast_database_is_ready,
    build_named_blast_database,
    collect_prepare_targets,
    describe_database_problem,
    format_prepare_error,
    hmm_indexing_available,
    prepare_index_is_ready,
    pyhmmer_install_hint,
)
from repeatviewer.backends.gui_functions import sanitize_name
from repeatviewer.database_settings import (
    reload_blast_databases_into_project,
    reload_managed_genomes_into_project,
    save_blast_databases_to_settings,
    save_database_settings,
    save_managed_genomes_to_settings,
)
from repeatviewer.pan_genome import GENOME_MODE_SINGLE
from repeatviewer.project import (
    ANNOTATION_NAME_MAX_LENGTH,
    BLAST_MOLECULE_NUCL,
    BLAST_MOLECULE_PROT,
    DATABASE_TYPE_CDD,
    DATABASE_TYPE_DNA,
    DB_KIND_CDD,
    DB_KIND_FASTA_DNA,
    DB_KIND_FASTA_PROTEIN,
    DB_KIND_HMM_DNA,
    DB_KIND_HMM_PROTEIN,
    DB_KIND_LABELS,
    DB_KINDS,
    BlastDatabase,
    ManagedGenome,
    RepeatViewerProject,
    UserDatabase,
    kind_label,
    normalize_annotation_name,
    user_fields_from_kind,
)
from repeatviewer.widgets.database_index_runner import find_database_index_runner
from repeatviewer.widgets.help_button import field_with_help, make_help_button
from repeatviewer.widgets.path_browse import PathBrowseRow

logger = logging.getLogger(__name__)

_NAME_HELP = (
    "<b>Name</b> identifies this entry in Manage Database and is used for "
    "BLAST index folders. It defaults to the file name; you can change it."
)
_ANNOTATION_NAME_HELP = (
    "<b>Annotation name</b> is the short label shown on the Annotation tab "
    f"(maximum {ANNOTATION_NAME_MAX_LENGTH} characters). "
    "Leave it empty if you do not use this database for annotation. "
    "When you check “Use for annotation”, it defaults to Name (truncated to "
    f"{ANNOTATION_NAME_MAX_LENGTH})."
)
_USE_FOR_ANNOTATION_HELP = (
    f"An Annotation name is required (max {ANNOTATION_NAME_MAX_LENGTH} characters)."
)
_SET_AS_MAIN_GENOME_HELP = (
    "Marks the selected FASTA DNA file as the reference genome. "
    "For example, this will be the default reference for BLASTn and TE-Aid."
)
_SET_MAIN_REFERENCE_HELP = (
    "Marks the selected FASTA DNA file as the reference genome. "
    "For example, this will be the default reference for BLASTn and TE-Aid."
)

# Legacy filter aliases used by callers that still pass role_filter=.
RoleFilter = Literal["all", "genome", "annotation", "blast"]
TypeFilter = Literal[
    "all",
    "fasta_dna",
    "fasta_protein",
    "hmm_dna",
    "hmm_protein",
    "cdd",
]
SectionFilter = Literal["blast", "hmm", "others"]
# What the caller can accept when using Manage Database → Select.
PickCompatible = Literal[
    "blast_nucl",
    "blast_prot",
    "blast_any",
    "hmm",
    "cdd",
]

_COL_NAME = 0
_COL_ANN_NAME = 1
_COL_TYPE = 2
_COL_PATH = 3
_COL_TIME = 4
_COL_ANNOTATION = 5
_COL_COUNT = 6

_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".faa"}

_BLAST_KINDS = {DB_KIND_FASTA_DNA, DB_KIND_FASTA_PROTEIN}
_HMM_KINDS = {DB_KIND_HMM_DNA, DB_KIND_HMM_PROTEIN}
_OTHER_KINDS = {DB_KIND_CDD}

_PICK_COMPATIBLE_LABELS: dict[str, str] = {
    "blast_nucl": "nucleotide BLAST / genome FASTA",
    "blast_prot": "protein BLAST FASTA",
    "blast_any": "BLAST FASTA",
    "hmm": "HMM database",
    "cdd": "CDD database",
}


def fasta_path_from_database_pick(payload: object) -> Path | None:
    """Resolve a Manage Database Select payload to a FASTA path when possible."""
    if isinstance(payload, BlastDatabase):
        return Path(payload.source_fasta).expanduser().resolve()
    if isinstance(payload, ManagedGenome):
        return Path(payload.path).expanduser().resolve()
    if isinstance(payload, UserDatabase):
        if payload.is_hmm or payload.db_type == DATABASE_TYPE_CDD:
            return None
        return Path(payload.path).expanduser().resolve()
    if payload is None:
        return None
    try:
        path = Path(payload).expanduser().resolve()
    except (TypeError, ValueError):
        return None
    return path


def database_name_from_pick(payload: object) -> str | None:
    """Display/registry name for combo widgets that store database names."""
    if isinstance(payload, BlastDatabase):
        return str(payload.name or "").strip() or None
    if isinstance(payload, ManagedGenome):
        return str(payload.name or "").strip() or None
    if isinstance(payload, UserDatabase):
        return str(payload.name or "").strip() or None
    return None


class _MiddleElideDelegate(QStyledItemDelegate):
    """Force middle-ellipsis like the Files panel (keep start + file name)."""

    def initStyleOption(self, option: QStyleOptionViewItem, index) -> None:  # noqa: N802
        super().initStyleOption(option, index)
        option.textElideMode = Qt.TextElideMode.ElideMiddle


class _EditableNameDelegate(_MiddleElideDelegate):
    """Inline editor for Name / Annotation name (double-click or Enter)."""

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        max_length: int | None = None,
        empty_display: str = "",
    ) -> None:
        super().__init__(parent)
        self._max_length = max_length
        self._empty_display = empty_display

    def createEditor(self, parent, option, index):  # noqa: N802
        editor = QLineEdit(parent)
        if self._max_length is not None:
            editor.setMaxLength(self._max_length)
        # Square corners — match the table cell (no rounded pill look).
        editor.setStyleSheet(
            "QLineEdit {"
            "  border: 1px solid palette(highlight);"
            "  border-radius: 0px;"
            "  padding: 0px 2px;"
            "  background: palette(base);"
            "  selection-background-color: palette(highlight);"
            "  selection-color: palette(highlighted-text);"
            "}"
        )
        return editor

    def setEditorData(self, editor, index) -> None:  # noqa: N802
        if not isinstance(editor, QLineEdit):
            super().setEditorData(editor, index)
            return
        real = index.data(Qt.ItemDataRole.UserRole)
        if real is None:
            text = str(index.data(Qt.ItemDataRole.DisplayRole) or "")
            if self._empty_display and text == self._empty_display:
                text = ""
        else:
            text = str(real)
        editor.setText(text)
        # Place caret at end — do not select the whole name on edit.
        # Defer past focus handling; macOS often selects-all when the editor gains focus.
        def _place_caret(ed: QLineEdit = editor) -> None:
            ed.deselect()
            ed.setCursorPosition(len(ed.text()))

        QTimer.singleShot(0, _place_caret)

    def setModelData(self, editor, model, index) -> None:  # noqa: N802
        if not isinstance(editor, QLineEdit):
            super().setModelData(editor, model, index)
            return
        text = editor.text().strip()
        if self._max_length is not None:
            text = text[: self._max_length]
        model.setData(index, text, Qt.ItemDataRole.UserRole)
        display = text if text else (self._empty_display or "")
        model.setData(index, display, Qt.ItemDataRole.DisplayRole)
        model.setData(index, display, Qt.ItemDataRole.EditRole)

def _now_iso() -> str:
    return datetime.now().isoformat(timespec="seconds")


def _existing_database_names(
    project: RepeatViewerProject,
    *,
    exclude_name: str | None = None,
) -> set[str]:
    """Lowercased Name values across BLAST, user/HMM/CDD, and managed genomes."""
    skip = (exclude_name or "").strip().lower()
    names: set[str] = set()
    for db in project.blast_databases:
        key = db.name.strip().lower()
        if key and key != skip:
            names.add(key)
    for db in project.user_databases:
        key = db.name.strip().lower()
        if key and key != skip:
            names.add(key)
    for genome in project.managed_genomes:
        key = genome.name.strip().lower()
        if key and key != skip:
            names.add(key)
    return names


def _existing_annotation_names(
    project: RepeatViewerProject,
    *,
    exclude_annotation_name: str | None = None,
) -> set[str]:
    """Lowercased non-empty Annotation names across BLAST and user databases."""
    skip = normalize_annotation_name(exclude_annotation_name).lower()
    names: set[str] = set()
    for db in project.blast_databases:
        key = normalize_annotation_name(db.annotation_name).lower()
        if key and key != skip:
            names.add(key)
    for db in project.user_databases:
        key = normalize_annotation_name(db.annotation_name).lower()
        if key and key != skip:
            names.add(key)
    return names


def _name_clash_message(name: str) -> str:
    return f"The name “{name}” is already used. Choose a unique Name."


def _annotation_name_clash_message(ann_name: str) -> str:
    return (
        f"The annotation name “{ann_name}” is already used. "
        "Choose a unique Annotation name."
    )


def _format_time(value: str) -> str:
    """Display Time as year-month-day only."""
    text = (value or "").strip()
    if not text:
        return ""
    if "T" in text:
        return text.split("T", 1)[0][:10]
    if " " in text:
        return text.split(" ", 1)[0][:10]
    return text[:10] if len(text) >= 10 else text


def _genome_path_key(path: Path) -> str:
    try:
        return str(path.expanduser().resolve())
    except OSError:
        return str(path.expanduser())


def _path_status_genome(path: Path) -> str:
    try:
        return "Ready" if path.expanduser().is_file() else "Missing"
    except OSError:
        return "Missing"


def _path_status_annotation(db: UserDatabase) -> str:
    path = db.path.expanduser()
    try:
        exists = path.is_file() or path.is_dir()
    except OSError:
        exists = False
    if not exists:
        return "Missing"
    kind = db.prepare_kind
    if describe_database_problem(kind, path) is not None:
        return "Missing"
    if not prepare_index_is_ready(kind, path):
        return "Needs index"
    return "Ready"


def _path_status_blast(db: BlastDatabase) -> str:
    try:
        if not db.source_fasta.expanduser().is_file():
            return "Missing"
    except OSError:
        return "Missing"
    if not blast_database_is_ready(db):
        return "Needs index"
    return "Ready"


def _ensure_chrom_sizes(path: Path) -> None:
    """Best-effort chromosome sizes sidecar for FASTA DNA."""
    try:
        from repeatviewer.backends.genome_prepare import ensure_genome_chrom_sizes_file

        ensure_genome_chrom_sizes_file(path)
    except Exception as exc:  # noqa: BLE001
        logger.warning("Could not write chrom.sizes for %s: %s", path, exc)


def _section_from_role(role: str) -> SectionFilter:
    if role in {"genome", "blast"}:
        return "blast"
    if role == "annotation":
        return "hmm"
    return "blast"


@dataclass
class _RefRow:
    name: str
    annotation_name: str
    kind: str
    path_display: str
    status: str
    created_at: str
    show_in_annotation: bool
    is_main: bool
    genome: ManagedGenome | None = None
    annotation: UserDatabase | None = None
    blast: BlastDatabase | None = None


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


class AddDatabaseDialog(QDialog):
    """Add a database of any type; destination tab is chosen from Type."""

    def __init__(
        self,
        project: RepeatViewerProject,
        *,
        initial_kind: str = DB_KIND_FASTA_PROTEIN,
        section: str = "blast",  # blast | hmm | others  (+ legacy aliases)
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        # Map legacy section names from older callers/tests.
        legacy = {
            "reference": "blast",
            "annotation": "hmm",
        }
        section = legacy.get(section, section)
        if section not in {"blast", "hmm", "others"}:
            section = "blast"
        self.setWindowTitle("Add database")
        self.setMinimumWidth(640)
        self.setModal(False)
        self.setWindowModality(Qt.WindowModality.NonModal)
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, True)
        self._project = project
        self._section = section  # opener tab (hint for initial type only)
        self._result_genome: ManagedGenome | None = None
        self._result_annotation: UserDatabase | None = None
        self._result_blast: BlastDatabase | None = None
        self._set_main = False
        self._thread = None
        self._worker = None
        self._progress: QProgressDialog | None = None

        self._type = QComboBox(self)
        for kind in DB_KINDS:
            self._type.addItem(DB_KIND_LABELS[kind], kind)
        idx = self._type.findData(initial_kind)
        if idx >= 0:
            self._type.setCurrentIndex(idx)
        elif self._type.count():
            self._type.setCurrentIndex(0)

        self._name = QLineEdit(self)
        self._name.setPlaceholderText("Defaults to file name")
        self._name.setMinimumWidth(400)

        self._ann_name = QLineEdit(self)
        self._ann_name.setMaxLength(ANNOTATION_NAME_MAX_LENGTH)
        self._ann_name.setPlaceholderText(
            f"Optional — Annotation tab label (max {ANNOTATION_NAME_MAX_LENGTH})"
        )
        self._ann_name.setMinimumWidth(400)

        self._path_row = PathBrowseRow(
            self,
            placeholder="Click a file in Files, or Browse",
            browse_title="Choose file",
            file_filter="All files (*)",
            mode="file",
            browse_handler=self._browse,
            minimum_width=400,
        )
        self._path = self._path_row.line_edit()

        self._show_annotation = QCheckBox("Use for annotation", self)
        self._as_main = QCheckBox("Set as reference genome", self)
        self._show_annotation_row = field_with_help(
            self._show_annotation,
            "Use for annotation",
            _USE_FOR_ANNOTATION_HELP,
            parent=self,
        )
        self._as_main_row = field_with_help(
            self._as_main,
            "Set as reference genome",
            _SET_AS_MAIN_GENOME_HELP,
            parent=self,
        )

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setSpacing(8)
        form.addRow("Type:", self._type)
        form.addRow(
            "Name:",
            field_with_help(self._name, "Name", _NAME_HELP, parent=self),
        )
        form.addRow(
            "Short name:",
            field_with_help(
                self._ann_name,
                "Short name",
                _ANNOTATION_NAME_HELP,
                parent=self,
            ),
        )
        form.addRow("Path:", self._path_row)
        form.addRow("", self._show_annotation_row)
        form.addRow("", self._as_main_row)

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

        self._type.currentIndexChanged.connect(self._sync_type_fields)
        self._path.textChanged.connect(self._maybe_fill_name_from_path)
        self._show_annotation.toggled.connect(self._sync_annotation_name_hint)
        self._sync_type_fields()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        QTimer.singleShot(0, self._path_row.focus_field)

    def _keep_in_front(self) -> None:
        if not self.isVisible():
            return
        self.raise_()
        self.activateWindow()

    def _current_kind(self) -> str:
        return str(self._type.currentData() or DB_KIND_FASTA_PROTEIN)

    def _destination_section(self) -> SectionFilter:
        kind = self._current_kind()
        if kind in _BLAST_KINDS:
            return "blast"
        if kind in _HMM_KINDS:
            return "hmm"
        return "others"

    def _sync_type_fields(self) -> None:
        kind = self._current_kind()
        is_dna = kind == DB_KIND_FASTA_DNA
        is_hmm = kind in _HMM_KINDS
        is_cdd = kind == DB_KIND_CDD

        self._as_main_row.setVisible(is_dna)
        if not is_dna:
            self._as_main.setChecked(False)

        if is_hmm or is_cdd:
            self._show_annotation.setChecked(True)
        if is_cdd and not self._name.text().strip():
            self._name.setText("CDD")
        self._sync_annotation_name_hint()

    def _sync_annotation_name_hint(self, _checked: bool = False) -> None:
        if self._show_annotation.isChecked():
            self._ann_name.setPlaceholderText(
                f"Used on Annotation tab (max {ANNOTATION_NAME_MAX_LENGTH})"
            )
        else:
            self._ann_name.setPlaceholderText(
                f"Optional — Annotation tab label (max {ANNOTATION_NAME_MAX_LENGTH})"
            )

    def _maybe_fill_name_from_path(self, _text: str = "") -> None:
        if self._name.text().strip():
            return
        path_text = self._path.text().strip()
        if not path_text:
            return
        path = Path(path_text)
        if self._current_kind() == DB_KIND_CDD:
            if not self._name.text().strip():
                self._name.setText(path.name or "CDD")
            return
        self._name.setText(path.stem)

    def apply_path_from_panel(self, path: Path) -> bool:
        resolved = path.expanduser().resolve()
        kind = self._current_kind()
        if kind == DB_KIND_CDD:
            if not resolved.is_dir():
                if resolved.is_file():
                    resolved = resolved.parent
                else:
                    return False
            self._path.setText(str(resolved))
            if not self._name.text().strip():
                self._name.setText("CDD")
            self._keep_in_front()
            return True
        if not resolved.is_file():
            return False
        if kind in _BLAST_KINDS:
            if resolved.suffix.lower() not in _FASTA_SUFFIXES:
                return False
        self._path.setText(str(resolved))
        if not self._name.text().strip():
            self._name.setText(resolved.stem)
        if resolved.suffix.lower() == ".faa" and kind == DB_KIND_FASTA_DNA:
            idx = self._type.findData(DB_KIND_FASTA_PROTEIN)
            if idx >= 0:
                self._type.setCurrentIndex(idx)
        self._keep_in_front()
        QTimer.singleShot(0, self._keep_in_front)
        return True

    def apply_folder_from_panel(self, folder: Path) -> bool:
        if self._current_kind() != DB_KIND_CDD:
            return False
        resolved = folder.expanduser().resolve()
        if not resolved.is_dir():
            return False
        self._path.setText(str(resolved))
        if not self._name.text().strip():
            self._name.setText("CDD")
        self._keep_in_front()
        return True

    def _browse(self) -> None:
        from PySide6.QtWidgets import QFileDialog

        kind = self._current_kind()
        start = self._path.text().strip()
        start_dir = str(Path(start).expanduser().parent) if start else str(Path.home())
        if kind == DB_KIND_CDD:
            chosen = QFileDialog.getExistingDirectory(self, "Choose CDD directory", start_dir)
            if chosen:
                self.apply_folder_from_panel(Path(chosen))
            return
        if kind in _BLAST_KINDS:
            filters = "FASTA (*.fa *.fasta *.fna *.fas *.faa);;All files (*)"
        elif kind in _HMM_KINDS:
            filters = "HMM (*.hmm);;All files (*)"
        else:
            filters = "All files (*)"
        chosen, _ = QFileDialog.getOpenFileName(self, "Choose file", start_dir, filters)
        if chosen:
            self.apply_path_from_panel(Path(chosen))

    def result_genome(self) -> ManagedGenome | None:
        return self._result_genome

    def result_annotation(self) -> UserDatabase | None:
        return self._result_annotation

    def result_blast(self) -> BlastDatabase | None:
        return self._result_blast

    def set_as_main(self) -> bool:
        return self._set_main

    def _on_ok(self) -> None:
        kind = self._current_kind()
        path_text = self._path.text().strip()
        if not path_text:
            QMessageBox.warning(self, self.windowTitle(), "Choose a path.")
            return
        path = Path(path_text).expanduser()

        name = self._name.text().strip() or path.stem
        if not name:
            QMessageBox.warning(self, self.windowTitle(), "Enter a name.")
            return

        show_ann = self._show_annotation.isChecked()
        ann_name = normalize_annotation_name(self._ann_name.text())
        # Annotation name is required when used for annotation; default to Name.
        if show_ann and not ann_name:
            ann_name = normalize_annotation_name(name)
            self._ann_name.setText(ann_name)

        if name.lower() in _existing_database_names(self._project):
            QMessageBox.warning(
                self, self.windowTitle(), _name_clash_message(name)
            )
            return
        if ann_name and ann_name.lower() in _existing_annotation_names(self._project):
            QMessageBox.warning(
                self,
                self.windowTitle(),
                _annotation_name_clash_message(ann_name),
            )
            return

        stamp = _now_iso()
        want_main = self._as_main.isChecked() and kind == DB_KIND_FASTA_DNA

        if kind == DB_KIND_CDD:
            if not path.is_dir():
                QMessageBox.warning(
                    self, self.windowTitle(), f"CDD directory not found:\n{path}"
                )
                return
        elif not path.is_file():
            QMessageBox.warning(self, self.windowTitle(), f"File not found:\n{path}")
            return

        if kind in _HMM_KINDS or kind == DB_KIND_CDD:
            db_type, is_hmm = user_fields_from_kind(kind)
            self._result_annotation = UserDatabase(
                name=name,
                path=path,
                db_type=db_type,
                is_hmm=is_hmm,
                created_at=stamp,
                show_in_annotation=show_ann,
                annotation_name=ann_name,
            )
            self.accept()
            return

        if kind not in _BLAST_KINDS:
            QMessageBox.warning(self, self.windowTitle(), f"Unsupported type: {kind}")
            return

        # FASTA DNA / Protein → BLAST tab (always build an index).
        molecule_type = (
            BLAST_MOLECULE_NUCL if kind == DB_KIND_FASTA_DNA else BLAST_MOLECULE_PROT
        )
        base_dir = load_blast_db_base_dir()
        safe_name = sanitize_name(name) or "database"
        target_dir = base_dir.expanduser().resolve() / safe_name
        if target_dir.exists():
            QMessageBox.warning(
                self,
                self.windowTitle(),
                f"A database folder already exists at:\n{target_dir}\n\n"
                "Choose a different name or file.",
            )
            return

        self._set_main = want_main
        if want_main:
            self._result_genome = ManagedGenome(name=name, path=path, created_at=stamp)

        self._buttons.setEnabled(False)
        progress = QProgressDialog("Building BLAST database…", "", 0, 0, self)
        progress.setWindowTitle(self.windowTitle())
        progress.setWindowModality(Qt.WindowModality.WindowModal)
        progress.setCancelButton(None)
        progress.setMinimumDuration(0)
        progress.setAutoClose(False)
        progress.setAutoReset(False)
        from repeatviewer.views.blast_database_dialog import style_blast_db_progress

        style_blast_db_progress(progress)
        progress.show()
        self._progress = progress
        self._pending_blast_name = name
        self._pending_blast_fasta = path
        self._pending_blast_molecule = molecule_type
        self._pending_stamp = stamp
        self._pending_show_annotation = show_ann
        self._pending_annotation_name = ann_name

        worker = _BuildDatabaseWorker(
            path,
            base_dir,
            name,
            molecule_type=molecule_type,
            use_system_blast=self._project.use_system_blast,
        )
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.finished.connect(self._on_blast_finished)
        worker.failed.connect(self._on_blast_failed)
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
    def _on_blast_finished(self, db_prefix: object) -> None:
        self._cleanup_progress()
        fasta = self._pending_blast_fasta
        molecule = self._pending_blast_molecule
        self._result_blast = BlastDatabase(
            name=self._pending_blast_name,
            source_fasta=fasta,
            db_prefix_path=Path(str(db_prefix)),
            molecule_type=molecule,
            created_at=self._pending_stamp,
            show_in_annotation=bool(getattr(self, "_pending_show_annotation", False)),
            annotation_name=str(getattr(self, "_pending_annotation_name", "")),
        )
        # Every FASTA DNA BLAST DB gets a chromosome sizes sidecar.
        if molecule == BLAST_MOLECULE_NUCL:
            _ensure_chrom_sizes(fasta)
        if self._set_main and self._result_genome is None:
            self._result_genome = ManagedGenome(
                name=self._pending_blast_name,
                path=fasta,
                created_at=self._pending_stamp,
            )
        self.accept()

    @Slot(str)
    def _on_blast_failed(self, message: str) -> None:
        self._cleanup_progress()
        QMessageBox.critical(self, "Add database", message)


# Back-compat alias for tests / callers.
AddReferenceDialog = AddDatabaseDialog


class ManageReferencesDialog(QDialog):
    """Manage Database — BLAST | HMM | Others."""

    references_changed = Signal()
    # Emitted when the user clicks Select (BlastDatabase | Path | UserDatabase).
    databaseSelected = Signal(object)

    def __init__(
        self,
        project: RepeatViewerProject,
        parent: QWidget | None = None,
        *,
        role_filter: RoleFilter | TypeFilter = "all",
        type_filter: TypeFilter | None = None,
        scope_filter: str = "all",  # ignored (legacy)
        pick_mode: bool = False,
    ) -> None:
        del scope_filter
        super().__init__(parent)
        self.setWindowTitle("Manage Database")
        self.setMinimumWidth(960)
        self.setMinimumHeight(420)
        # Independent floating window with a normal title bar (movable).
        # Avoid parenting under another modal dialog — on macOS that becomes a
        # non-draggable sheet sharing the host's title.
        self.setModal(False)
        self.setWindowModality(Qt.WindowModality.NonModal)
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, True)
        self._project = project
        self._add_dialog: AddDatabaseDialog | None = None
        self._rows: list[_RefRow] = []
        self._sort_column = _COL_NAME
        self._sort_order = Qt.SortOrder.AscendingOrder
        self._pick_mode = bool(pick_mode)
        self._pick_callback = None
        self._pick_compatible: frozenset[str] = frozenset()
        self._search_text = ""
        del type_filter  # type filter UI removed; tab already scopes kinds

        self._legacy_role = (
            role_filter if role_filter in {"genome", "annotation", "blast"} else "all"
        )

        # Compact left-aligned tabs (Alignment-style: no full-width base line).
        self._section_bar = QTabBar(self)
        self._section_bar.setExpanding(False)
        self._section_bar.setDocumentMode(True)
        self._section_bar.setDrawBase(False)
        self._section_bar.addTab("BLAST database")
        self._section_bar.addTab("HMM database")
        self._section_bar.addTab("Others")
        section = _section_from_role(self._legacy_role)
        self._section_bar.setCurrentIndex({"blast": 0, "hmm": 1, "others": 2}[section])

        # Legacy attribute kept so older callers do not break.
        self._role_filter = self._section_bar
        self._scope_filter = self._section_bar
        self._type_filter = self._section_bar

        tab_row = QHBoxLayout()
        tab_row.setContentsMargins(0, 0, 0, 0)
        tab_row.addWidget(self._section_bar, 0, Qt.AlignmentFlag.AlignLeft)
        tab_row.addStretch(1)

        self._search = QLineEdit(self)
        self._search.setPlaceholderText("Search name, annotation, type, or path…")
        self._search.setClearButtonEnabled(True)
        self._search.textChanged.connect(self._on_search_changed)

        self._table = QTableWidget(0, _COL_COUNT, self)
        self._table.setHorizontalHeaderLabels(
            [
                "Name",
                "Annotation name",
                "Type",
                "Path / source",
                "Time",
                "Annotation",
            ]
        )
        self._table.setTextElideMode(Qt.TextElideMode.ElideMiddle)
        # Explicit delegate — macOS styles sometimes ignore view elide mode.
        self._elide_delegate = _MiddleElideDelegate(self._table)
        self._name_delegate = _EditableNameDelegate(self._table)
        self._ann_name_delegate = _EditableNameDelegate(
            self._table,
            max_length=ANNOTATION_NAME_MAX_LENGTH,
            empty_display="-",
        )
        self._table.setItemDelegateForColumn(_COL_NAME, self._name_delegate)
        self._table.setItemDelegateForColumn(_COL_ANN_NAME, self._ann_name_delegate)
        self._table.setItemDelegateForColumn(_COL_PATH, self._elide_delegate)
        header = self._table.horizontalHeader()
        header.setSectionsClickable(True)
        header.setSortIndicatorShown(True)
        header.setSortIndicator(self._sort_column, self._sort_order)
        header.sectionClicked.connect(self._on_header_clicked)
        # Name gets leftover width so long database names stay readable.
        header.setSectionResizeMode(_COL_NAME, QHeaderView.ResizeMode.Stretch)
        header.setMinimumSectionSize(80)
        header.setSectionResizeMode(_COL_ANN_NAME, QHeaderView.ResizeMode.Interactive)
        header.resizeSection(_COL_ANN_NAME, 110)
        header.setSectionResizeMode(_COL_TYPE, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(_COL_PATH, QHeaderView.ResizeMode.Interactive)
        header.resizeSection(_COL_PATH, 200)
        header.setSectionResizeMode(_COL_TIME, QHeaderView.ResizeMode.ResizeToContents)
        header.setSectionResizeMode(_COL_ANNOTATION, QHeaderView.ResizeMode.ResizeToContents)
        self._table.verticalHeader().setVisible(False)
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        # Double-click or press Enter / type to edit Name and Annotation name.
        self._table.setEditTriggers(
            QAbstractItemView.EditTrigger.DoubleClicked
            | QAbstractItemView.EditTrigger.EditKeyPressed
            | QAbstractItemView.EditTrigger.AnyKeyPressed
        )
        self._table.itemSelectionChanged.connect(self._update_buttons)
        self._table.itemChanged.connect(self._on_item_changed)

        self._add_btn = QPushButton("Add…", self)
        self._remove_btn = QPushButton("Remove", self)
        self._select_btn = QPushButton("Select", self)
        self._select_btn.setToolTip(
            "Fill the database field that opened this window "
            "(only compatible types are accepted)."
        )
        self._select_btn.clicked.connect(self._select_current)
        self._table.doubleClicked.connect(self._on_table_double_clicked)
        self._set_main_btn = QPushButton("Set as reference genome", self)
        self._set_main_help = make_help_button(
            self, "Set as reference genome", _SET_MAIN_REFERENCE_HELP
        )
        self._add_btn.clicked.connect(self._add_database)
        self._remove_btn.clicked.connect(self._remove_selected)
        self._set_main_btn.clicked.connect(self._set_main_reference)

        self._save_btn = QPushButton("Save", self)
        self._save_btn.setDefault(True)
        self._save_btn.clicked.connect(self.accept)

        btn_row = QHBoxLayout()
        btn_row.addWidget(self._add_btn)
        btn_row.addWidget(self._remove_btn)
        btn_row.addWidget(self._select_btn)
        btn_row.addWidget(self._set_main_btn)
        btn_row.addWidget(self._set_main_help)
        btn_row.addStretch(1)
        btn_row.addWidget(self._save_btn)

        self._current_label = QLabel(self)
        self._current_label.setWordWrap(True)

        layout = QVBoxLayout(self)
        layout.addLayout(tab_row)
        layout.addWidget(self._search)
        layout.addWidget(self._table, 1)
        layout.addLayout(btn_row)
        layout.addWidget(self._current_label)

        self._section_bar.currentChanged.connect(self._on_section_changed)
        self._reload_from_settings()
        self._sync_section_chrome()
        self._populate()
        self._update_buttons()
        if self._pick_mode:
            self.set_pick_mode(True)

    def set_pick_mode(
        self,
        enabled: bool = True,
        *,
        callback=None,
        role_filter: str | None = None,
        compatible: PickCompatible | tuple[PickCompatible, ...] | list[str] | None = None,
    ) -> None:
        """Enable Select-to-fill mode for a tool database field.

        ``compatible`` limits Select to matching kinds, e.g. ``blast_nucl`` for
        TE-Aid / BLASTn, ``blast_prot`` for BLASTp/x, ``hmm`` for HMMER.
        """
        self._pick_mode = bool(enabled)
        self._pick_callback = callback if enabled else None
        if not enabled:
            self._pick_compatible = frozenset()
        elif compatible is None:
            self._pick_compatible = frozenset()
        elif isinstance(compatible, str):
            self._pick_compatible = frozenset({compatible})
        else:
            self._pick_compatible = frozenset(str(item) for item in compatible if item)
        if role_filter is not None:
            self.apply_role_filter(role_filter)
        elif self._pick_compatible:
            # Jump to the tab that matches the requested kind.
            if self._pick_compatible & {"blast_nucl", "blast_prot", "blast_any"}:
                self.apply_role_filter("blast")
            elif "hmm" in self._pick_compatible:
                self.apply_role_filter("annotation")
            elif "cdd" in self._pick_compatible:
                self.apply_role_filter("all")
                self._section_bar.setCurrentIndex(2)
        self._select_btn.setDefault(self._pick_mode)
        self._update_select_tooltip()
        self._update_buttons()
        if self._pick_mode:
            self.raise_()
            self.activateWindow()

    def _update_select_tooltip(self) -> None:
        if not self._pick_compatible:
            self._select_btn.setToolTip(
                "Fill the database field that opened this window."
            )
            return
        labels = [
            _PICK_COMPATIBLE_LABELS[key]
            for key in sorted(self._pick_compatible)
            if key in _PICK_COMPATIBLE_LABELS
        ]
        need = ", ".join(labels) if labels else "a compatible database"
        self._select_btn.setToolTip(f"Fill the opener’s field. Requires: {need}.")

    def apply_role_filter(self, role_filter: str) -> None:
        """Switch tab when an existing dialog is re-opened with a role hint."""
        section = _section_from_role(role_filter)
        idx = {"blast": 0, "hmm": 1, "others": 2}[section]
        self._legacy_role = (
            role_filter if role_filter in {"genome", "annotation", "blast"} else "all"
        )
        if self._section_bar.currentIndex() != idx:
            self._section_bar.setCurrentIndex(idx)
        else:
            self._sync_section_chrome()
            self._populate()

    def _on_section_changed(self, _index: int = 0) -> None:
        self._legacy_role = "all"
        self._sync_section_chrome()
        self._populate()

    def _current_section(self) -> SectionFilter:
        idx = self._section_bar.currentIndex() if hasattr(self, "_section_bar") else 0
        return {0: "blast", 1: "hmm", 2: "others"}.get(idx, "blast")  # type: ignore[return-value]

    def _sync_section_chrome(self) -> None:
        section = self._current_section()
        show_main = section == "blast"
        self._set_main_btn.setVisible(show_main)
        self._set_main_help.setVisible(show_main)

    def _on_header_clicked(self, logical_index: int) -> None:
        if logical_index < 0 or logical_index >= _COL_COUNT:
            return
        if logical_index == self._sort_column:
            self._sort_order = (
                Qt.SortOrder.DescendingOrder
                if self._sort_order == Qt.SortOrder.AscendingOrder
                else Qt.SortOrder.AscendingOrder
            )
        else:
            self._sort_column = logical_index
            self._sort_order = Qt.SortOrder.AscendingOrder
        self._table.horizontalHeader().setSortIndicator(
            self._sort_column, self._sort_order
        )
        self._populate()

    def _sort_rows(self, rows: list[_RefRow]) -> list[_RefRow]:
        col = self._sort_column
        reverse = self._sort_order == Qt.SortOrder.DescendingOrder

        def key(row: _RefRow) -> tuple:
            if col == _COL_NAME:
                return (row.name.lower(),)
            if col == _COL_ANN_NAME:
                return (row.annotation_name.lower(),)
            if col == _COL_TYPE:
                return (kind_label(row.kind).lower(),)
            if col == _COL_PATH:
                return (row.path_display.lower(),)
            if col == _COL_TIME:
                return (_format_time(row.created_at), row.created_at)
            if col == _COL_ANNOTATION:
                return (0 if row.show_in_annotation else 1, row.name.lower())
            return (row.name.lower(),)

        return sorted(rows, key=key, reverse=reverse)

    def project(self) -> RepeatViewerProject:
        return self._project

    def apply_path_from_panel(self, path: Path) -> bool:
        dialog = self._add_dialog
        if dialog is None or not dialog.isVisible():
            return False
        return dialog.apply_path_from_panel(path)

    def apply_folder_from_panel(self, folder: Path) -> bool:
        dialog = self._add_dialog
        if dialog is None or not dialog.isVisible():
            return False
        return dialog.apply_folder_from_panel(folder)

    def _reload_from_settings(self) -> None:
        from repeatviewer.backends.database_prepare import (
            blast_database_is_ready,
            blast_search_db_path,
        )
        from repeatviewer.database_settings import (
            load_cdd_dirs_from_settings,
            load_user_databases_from_settings,
        )

        self._project.user_databases = load_user_databases_from_settings()
        reload_blast_databases_into_project(self._project)
        reload_managed_genomes_into_project(self._project)
        changed = False
        # Promote legacy cdd_dirs-only CDD into a tracked user DB so Annotation
        # name survives close/reopen.
        if not any(db.db_type == DATABASE_TYPE_CDD for db in self._project.user_databases):
            cdd_dirs = load_cdd_dirs_from_settings()
            if cdd_dirs:
                self._project.user_databases.append(
                    UserDatabase(
                        name="CDD",
                        path=cdd_dirs[0],
                        db_type=DATABASE_TYPE_CDD,
                        show_in_annotation=True,
                        annotation_name="CDD",
                    )
                )
                self._project.cdd_dirs = [cdd_dirs[0]]
                changed = True
        # Promote annotation FASTA rows that already have a BLAST index into
        # blast_databases so BLASTn/p/x dropdowns list them.
        blast_paths = {
            _genome_path_key(db.source_fasta) for db in self._project.blast_databases
        }
        blast_names = {db.name.lower() for db in self._project.blast_databases}
        for db in self._project.user_databases:
            if db.kind not in _BLAST_KINDS:
                continue
            key = _genome_path_key(db.path)
            if key in blast_paths:
                continue
            protein = db.kind == DB_KIND_FASTA_PROTEIN
            prefix = blast_search_db_path(db.path, protein=protein)
            name = db.name
            if name.lower() in blast_names:
                name = f"{db.name}_{'prot' if protein else 'nucl'}"
            candidate = BlastDatabase(
                name=name,
                source_fasta=db.path,
                db_prefix_path=prefix,
                molecule_type=(
                    BLAST_MOLECULE_PROT if protein else BLAST_MOLECULE_NUCL
                ),
                created_at=db.created_at or _now_iso(),
                show_in_annotation=db.show_in_annotation,
                annotation_name=db.annotation_name,
            )
            if not blast_database_is_ready(candidate):
                continue
            self._project.blast_databases.append(candidate)
            blast_paths.add(key)
            blast_names.add(candidate.name.lower())
            changed = True
        # Re-apply annotation-name length cap; require a name when used for annotation.
        for db in self._project.user_databases:
            if db.show_in_annotation and not db.annotation_name.strip():
                db.annotation_name = normalize_annotation_name(db.name)
                changed = True
            capped = normalize_annotation_name(db.annotation_name)
            if capped != db.annotation_name:
                db.annotation_name = capped
                changed = True
        for db in self._project.blast_databases:
            if db.show_in_annotation and not db.annotation_name.strip():
                db.annotation_name = normalize_annotation_name(db.name)
                changed = True
            capped = normalize_annotation_name(db.annotation_name)
            if capped != db.annotation_name:
                db.annotation_name = capped
                changed = True
        if changed:
            save_database_settings(
                cdd_dirs=self._project.cdd_dirs,
                user_databases=self._project.user_databases,
            )
            save_blast_databases_to_settings(self._project.blast_databases)

    def _is_main_genome(self, path: Path) -> bool:
        current = self._project.genome_file
        if current is None:
            return False
        return _genome_path_key(current) == _genome_path_key(path)

    def _build_rows(self) -> list[_RefRow]:
        rows: list[_RefRow] = []
        blast_paths = {
            _genome_path_key(db.source_fasta) for db in self._project.blast_databases
        }
        # FASTA user DBs keyed by path so genome-only rows can share Annotation state.
        fasta_ann_by_path: dict[str, UserDatabase] = {}
        other_user_dbs: list[UserDatabase] = []
        for db in self._project.user_databases:
            if db.kind in _BLAST_KINDS:
                fasta_ann_by_path[_genome_path_key(db.path)] = db
            else:
                other_user_dbs.append(db)

        for db in self._project.blast_databases:
            path_key = _genome_path_key(db.source_fasta)
            linked = fasta_ann_by_path.pop(path_key, None)
            show_ann = bool(db.show_in_annotation) or (
                linked.show_in_annotation if linked is not None else False
            )
            ann_name = normalize_annotation_name(
                db.annotation_name
                or (linked.annotation_name if linked is not None else "")
            )
            rows.append(
                _RefRow(
                    name=db.name,
                    annotation_name=ann_name,
                    kind=db.kind,
                    path_display=str(db.source_fasta),
                    status=_path_status_blast(db),
                    created_at=db.created_at,
                    show_in_annotation=show_ann,
                    is_main=(
                        db.molecule_type == BLAST_MOLECULE_NUCL
                        and self._is_main_genome(db.source_fasta)
                    ),
                    blast=db,
                    annotation=linked,
                )
            )

        # Legacy managed genomes without a matching BLAST row still appear on BLAST tab.
        for g in self._project.managed_genomes:
            key = _genome_path_key(g.path)
            if key in blast_paths:
                continue
            linked = fasta_ann_by_path.pop(key, None)
            rows.append(
                _RefRow(
                    name=g.name,
                    annotation_name=normalize_annotation_name(
                        linked.annotation_name if linked is not None else ""
                    ),
                    kind=g.kind,
                    path_display=str(g.path),
                    status=_path_status_genome(g.path),
                    created_at=g.created_at,
                    show_in_annotation=(
                        linked.show_in_annotation if linked is not None else False
                    ),
                    is_main=self._is_main_genome(g.path),
                    genome=g,
                    annotation=linked,
                )
            )

        # Remaining FASTA user DBs that are neither genomes nor BLAST sources.
        for db in fasta_ann_by_path.values():
            rows.append(
                _RefRow(
                    name=db.name,
                    annotation_name=normalize_annotation_name(db.annotation_name),
                    kind=db.kind,
                    path_display=str(db.path),
                    status=_path_status_annotation(db),
                    created_at=db.created_at,
                    show_in_annotation=db.show_in_annotation,
                    is_main=False,
                    annotation=db,
                )
            )

        for db in other_user_dbs:
            rows.append(
                _RefRow(
                    name=db.name,
                    annotation_name=normalize_annotation_name(db.annotation_name),
                    kind=db.kind,
                    path_display=str(db.path),
                    status=_path_status_annotation(db),
                    created_at=db.created_at,
                    show_in_annotation=db.show_in_annotation,
                    is_main=False,
                    annotation=db,
                )
            )

        from repeatviewer.database_settings import load_cdd_dirs_from_settings

        if not any(db.db_type == DATABASE_TYPE_CDD for db in self._project.user_databases):
            cdd_dirs = load_cdd_dirs_from_settings()
            if cdd_dirs:
                cdd = UserDatabase(
                    name="CDD",
                    path=cdd_dirs[0],
                    db_type=DATABASE_TYPE_CDD,
                    show_in_annotation=True,
                    annotation_name="CDD",
                )
                rows.append(
                    _RefRow(
                        name=cdd.name,
                        annotation_name=cdd.annotation_name,
                        kind=cdd.kind,
                        path_display=str(cdd.path),
                        status=_path_status_annotation(cdd),
                        created_at=cdd.created_at,
                        show_in_annotation=True,
                        is_main=False,
                        annotation=cdd,
                    )
                )

        section = self._current_section()
        if section == "blast" or self._legacy_role in {"blast", "genome"}:
            rows = [
                r
                for r in rows
                if r.blast is not None
                or r.genome is not None
                or (r.annotation is not None and r.kind in _BLAST_KINDS)
            ]
        elif section == "hmm" or self._legacy_role == "annotation":
            rows = [r for r in rows if r.annotation is not None and r.kind in _HMM_KINDS]
        else:
            rows = [
                r for r in rows if r.annotation is not None and r.kind in _OTHER_KINDS
            ]
        # Hide entries that are not ready (no Status column).
        return [r for r in rows if r.status == "Ready"]

    def _on_search_changed(self, text: str) -> None:
        self._search_text = str(text or "").strip().lower()
        self._populate()

    def _row_matches_search(self, row: _RefRow) -> bool:
        needle = self._search_text
        if not needle:
            return True
        haystacks = (
            row.name,
            row.annotation_name,
            kind_label(row.kind),
            row.path_display,
        )
        return any(needle in str(part or "").lower() for part in haystacks)

    def _populate(self) -> None:
        filtered = [
            row for row in self._sort_rows(self._build_rows()) if self._row_matches_search(row)
        ]
        self._rows = filtered
        self._table.blockSignals(True)
        self._table.setRowCount(0)
        self._table.setRowCount(len(filtered))
        editable = (
            Qt.ItemFlag.ItemIsEnabled
            | Qt.ItemFlag.ItemIsSelectable
            | Qt.ItemFlag.ItemIsEditable
        )
        readonly = Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsSelectable
        for i, row in enumerate(filtered):
            ann_display = row.annotation_name if row.annotation_name.strip() else "-"
            name_item = QTableWidgetItem(row.name)
            name_item.setFlags(editable)
            name_item.setData(Qt.ItemDataRole.UserRole, row.name)
            name_item.setToolTip(f"{row.name}\nDouble-click to rename")
            self._table.setItem(i, _COL_NAME, name_item)

            ann_item = QTableWidgetItem(ann_display)
            ann_item.setFlags(editable)
            ann_item.setData(Qt.ItemDataRole.UserRole, row.annotation_name)
            tip = row.annotation_name or "(empty)"
            ann_item.setToolTip(
                f"{tip}\nDouble-click to edit (max {ANNOTATION_NAME_MAX_LENGTH})"
            )
            self._table.setItem(i, _COL_ANN_NAME, ann_item)

            type_item = QTableWidgetItem(kind_label(row.kind))
            type_item.setFlags(readonly)
            self._table.setItem(i, _COL_TYPE, type_item)

            path_item = QTableWidgetItem(row.path_display)
            path_item.setFlags(readonly)
            path_item.setToolTip(row.path_display)
            self._table.setItem(i, _COL_PATH, path_item)

            time_item = QTableWidgetItem(_format_time(row.created_at))
            time_item.setFlags(readonly)
            self._table.setItem(i, _COL_TIME, time_item)

            check_item = QTableWidgetItem()
            check_item.setFlags(
                Qt.ItemFlag.ItemIsEnabled
                | Qt.ItemFlag.ItemIsSelectable
                | Qt.ItemFlag.ItemIsUserCheckable
            )
            check_item.setCheckState(
                Qt.CheckState.Checked
                if row.show_in_annotation
                else Qt.CheckState.Unchecked
            )
            self._table.setItem(i, _COL_ANNOTATION, check_item)
        self._table.blockSignals(False)
        self._update_current_label()
        self._update_buttons()

    def _update_current_label(self) -> None:
        genome = self._project.genome_file
        if genome is None:
            self._current_label.setText("Reference genome: (none)")
            return
        # Prefer the managed / BLAST display name when available.
        label = genome.name
        key = _genome_path_key(genome)
        for g in self._project.managed_genomes:
            if _genome_path_key(g.path) == key and g.name:
                label = g.name
                break
        else:
            for db in self._project.blast_databases:
                if (
                    db.molecule_type == BLAST_MOLECULE_NUCL
                    and _genome_path_key(db.source_fasta) == key
                    and db.name
                ):
                    label = db.name
                    break
        self._current_label.setText(f"Reference genome: {label} — {genome}")
        self._current_label.setToolTip(str(genome))

    def _selected_rows(self) -> list[_RefRow]:
        indices = sorted({idx.row() for idx in self._table.selectedIndexes()})
        return [self._rows[i] for i in indices if 0 <= i < len(self._rows)]

    def _update_buttons(self) -> None:
        selected = self._selected_rows()
        self._remove_btn.setEnabled(bool(selected))
        can_select = (
            len(selected) == 1
            and self._selection_payload(selected[0]) is not None
            and self._row_is_compatible(selected[0])
        )
        self._select_btn.setEnabled(can_select)
        can_main = False
        if len(selected) == 1 and self._current_section() == "blast":
            row = selected[0]
            if row.kind == DB_KIND_FASTA_DNA and (
                row.blast is not None
                or row.genome is not None
                or (row.annotation is not None and row.kind == DB_KIND_FASTA_DNA)
            ):
                can_main = True
        self._set_main_btn.setEnabled(can_main)

    def _row_is_compatible(self, row: _RefRow) -> bool:
        """True when Select may fill the opener with this row."""
        needed = self._pick_compatible
        if not needed:
            # No restriction (or not in pick mode): any resolvable payload is OK.
            return self._selection_payload(row) is not None

        kind = row.kind
        blast = row.blast
        if "blast_nucl" in needed:
            if blast is not None and blast.molecule_type == BLAST_MOLECULE_NUCL:
                return True
            if row.genome is not None:
                return True
            if kind == DB_KIND_FASTA_DNA:
                return True
        if "blast_prot" in needed:
            if blast is not None and blast.molecule_type == BLAST_MOLECULE_PROT:
                return True
            if kind == DB_KIND_FASTA_PROTEIN:
                return True
        if "blast_any" in needed:
            if blast is not None or row.genome is not None:
                return True
            if kind in _BLAST_KINDS:
                return True
        if "hmm" in needed:
            if kind in _HMM_KINDS:
                return True
            if row.annotation is not None and row.annotation.is_hmm:
                return True
        if "cdd" in needed:
            if kind == DB_KIND_CDD:
                return True
            if (
                row.annotation is not None
                and row.annotation.db_type == DATABASE_TYPE_CDD
            ):
                return True
        return False

    def _compatible_requirement_text(self) -> str:
        if not self._pick_compatible:
            return "a database"
        labels = [
            _PICK_COMPATIBLE_LABELS[key]
            for key in sorted(self._pick_compatible)
            if key in _PICK_COMPATIBLE_LABELS
        ]
        return ", ".join(labels) if labels else "a compatible database"

    def _selection_payload(self, row: _RefRow):
        """Return the database object / path for Select, or None if unsupported."""
        if row.blast is not None:
            return row.blast
        if row.genome is not None:
            return row.genome.path
        if row.annotation is not None:
            return row.annotation
        return None

    def _on_table_double_clicked(self, _index) -> None:  # noqa: ANN001
        if self._pick_mode:
            self._select_current()

    def _select_current(self) -> None:
        selected = self._selected_rows()
        if len(selected) != 1:
            QMessageBox.information(
                self,
                "Select",
                "Select one database row first.",
            )
            return
        row = selected[0]
        if not self._row_is_compatible(row):
            QMessageBox.information(
                self,
                "Select",
                f"“{row.name}” is not compatible here.\n"
                f"This field needs: {self._compatible_requirement_text()}.",
            )
            return
        payload = self._selection_payload(row)
        if payload is None:
            QMessageBox.information(
                self,
                "Select",
                "That row cannot be used as a database source.",
            )
            return
        if not self._pick_mode and not callable(self._pick_callback):
            QMessageBox.information(
                self,
                "Select",
                "Open Manage Database from the tool field you want to fill "
                "(for example TE-Aid, BLAST, HMMER, or Transfer annotation), "
                "then click Select.",
            )
            return
        callback = self._pick_callback
        self.databaseSelected.emit(payload)
        if callable(callback):
            try:
                callback(payload)
            except Exception:  # noqa: BLE001
                logging.exception("Manage Database select callback failed")
        name = row.name or "database"
        if self._pick_mode:
            self._pick_callback = None
            self._pick_mode = False
            self._pick_compatible = frozenset()
            self._select_btn.setDefault(False)
            self._update_select_tooltip()
            self.accept()
            return
        self._current_label.setText(f"Selected: {name}")

    def _ensure_user_db_tracked(self, db: UserDatabase) -> UserDatabase:
        """Return the project-owned UserDatabase for ``db`` (append if missing).

        Synthetic CDD rows built from ``cdd_dirs`` are not in ``user_databases``
        until Annotation is toggled; without this, annotation_name never persists.
        """
        key = (_genome_path_key(db.path), db.db_type, bool(db.is_hmm))
        for existing in self._project.user_databases:
            ek = (_genome_path_key(existing.path), existing.db_type, bool(existing.is_hmm))
            if ek == key:
                existing.name = db.name or existing.name
                existing.annotation_name = db.annotation_name
                existing.show_in_annotation = db.show_in_annotation
                existing.created_at = existing.created_at or db.created_at
                return existing
        self._project.user_databases.append(db)
        if db.db_type == DATABASE_TYPE_CDD:
            self._project.cdd_dirs = [db.path]
        return db

    def _rebind_blast(self, row: _RefRow) -> BlastDatabase | None:
        """Return the live BlastDatabase for *row*, rebinding if the list was replaced."""
        blast = row.blast
        if blast is None:
            return None
        if any(db is blast for db in self._project.blast_databases):
            return blast
        key = _genome_path_key(blast.source_fasta)
        for db in self._project.blast_databases:
            if _genome_path_key(db.source_fasta) == key or db.name == row.name:
                row.blast = db
                return db
        return blast

    def _rebind_user_db(self, row: _RefRow) -> UserDatabase | None:
        """Return the live UserDatabase for *row*, rebinding if the list was replaced."""
        annotation = row.annotation
        if annotation is None:
            return None
        if any(db is annotation for db in self._project.user_databases):
            return annotation
        key = _genome_path_key(annotation.path)
        for db in self._project.user_databases:
            if _genome_path_key(db.path) == key or db.name == row.name:
                row.annotation = db
                return db
        return annotation

    def _on_item_changed(self, item: QTableWidgetItem) -> None:
        col = item.column()
        row_idx = item.row()
        if not (0 <= row_idx < len(self._rows)):
            return
        if col == _COL_NAME:
            self._on_name_edited(row_idx, item)
            return
        if col == _COL_ANN_NAME:
            self._on_annotation_name_edited(row_idx, item)
            return
        if col != _COL_ANNOTATION:
            return
        row = self._rows[row_idx]
        checked = item.checkState() == Qt.CheckState.Checked

        def _set_ann_name(raw: str) -> str:
            ann = normalize_annotation_name(raw)
            row.annotation_name = ann
            ann_cell = self._table.item(row_idx, _COL_ANN_NAME)
            if ann_cell is not None:
                self._table.blockSignals(True)
                ann_cell.setData(Qt.ItemDataRole.UserRole, ann)
                ann_cell.setText(ann if ann else "-")
                ann_cell.setToolTip(
                    f"{ann or '(empty)'}\nDouble-click to edit "
                    f"(max {ANNOTATION_NAME_MAX_LENGTH})"
                )
                self._table.blockSignals(False)
            return ann

        blast = self._rebind_blast(row)
        if blast is not None:
            if checked:
                raw = blast.annotation_name or blast.name
                blast.annotation_name = _set_ann_name(raw)
            blast.show_in_annotation = checked
            row.show_in_annotation = checked
            annotation = self._rebind_user_db(row)
            if annotation is not None:
                row.annotation = self._ensure_user_db_tracked(annotation)
                row.annotation.show_in_annotation = checked
                if checked:
                    row.annotation.annotation_name = blast.annotation_name
            self._persist_all()
            return

        annotation = self._rebind_user_db(row)
        if annotation is not None:
            if checked:
                raw = annotation.annotation_name or annotation.name
                annotation.annotation_name = _set_ann_name(raw)
            annotation.show_in_annotation = checked
            row.show_in_annotation = checked
            row.annotation = self._ensure_user_db_tracked(annotation)
            self._persist_all()
            return

        # Genome-only row: create a FASTA DNA user DB so Annotation can use it.
        if row.genome is not None and checked:
            ann_name = _set_ann_name(row.genome.name)
            user_db = UserDatabase(
                name=row.genome.name,
                path=row.genome.path,
                db_type=DATABASE_TYPE_DNA,
                created_at=row.genome.created_at or _now_iso(),
                show_in_annotation=True,
                annotation_name=ann_name,
            )
            row.annotation = self._ensure_user_db_tracked(user_db)
            row.show_in_annotation = True
            self._persist_all()
            self._populate()
            return
        if row.genome is not None and not checked:
            row.show_in_annotation = False
            self._persist_all()
            return

    def _on_name_edited(self, row_idx: int, item: QTableWidgetItem) -> None:
        row = self._rows[row_idx]
        new_name = str(item.data(Qt.ItemDataRole.UserRole) or item.text()).strip()
        if not new_name:
            self._revert_name_cell(item, row.name)
            QMessageBox.warning(self, "Manage Database", "Name cannot be empty.")
            return
        if new_name == row.name:
            return
        if new_name.lower() in _existing_database_names(
            self._project, exclude_name=row.name
        ):
            self._revert_name_cell(item, row.name)
            QMessageBox.warning(
                self, "Manage Database", _name_clash_message(new_name)
            )
            return
        if row.blast is not None:
            row.blast.name = new_name
        elif row.annotation is not None:
            row.annotation = self._ensure_user_db_tracked(row.annotation)
            row.annotation.name = new_name
        elif row.genome is not None:
            row.genome.name = new_name
        else:
            self._revert_name_cell(item, row.name)
            return
        row.name = new_name
        self._table.blockSignals(True)
        item.setText(new_name)
        item.setData(Qt.ItemDataRole.UserRole, new_name)
        item.setToolTip(f"{new_name}\nDouble-click to rename")
        self._table.blockSignals(False)
        self._persist_all()
        self._update_current_label()

    def _on_annotation_name_edited(self, row_idx: int, item: QTableWidgetItem) -> None:
        row = self._rows[row_idx]
        raw = str(item.data(Qt.ItemDataRole.UserRole) or "")
        if raw == "-":
            raw = ""
        new_ann = normalize_annotation_name(raw)
        if row.show_in_annotation and not new_ann:
            new_ann = normalize_annotation_name(row.name)
        if new_ann == row.annotation_name:
            self._table.blockSignals(True)
            item.setText(new_ann if new_ann else "-")
            item.setData(Qt.ItemDataRole.UserRole, new_ann)
            self._table.blockSignals(False)
            return
        if new_ann and new_ann.lower() in _existing_annotation_names(
            self._project, exclude_annotation_name=row.annotation_name
        ):
            self._revert_annotation_name_cell(item, row.annotation_name)
            QMessageBox.warning(
                self,
                "Manage Database",
                _annotation_name_clash_message(new_ann),
            )
            return
        row.annotation_name = new_ann
        if row.blast is not None:
            row.blast.annotation_name = new_ann
            if row.annotation is not None:
                row.annotation = self._ensure_user_db_tracked(row.annotation)
                row.annotation.annotation_name = new_ann
        elif row.annotation is not None:
            row.annotation = self._ensure_user_db_tracked(row.annotation)
            row.annotation.annotation_name = new_ann
        elif row.genome is not None and new_ann:
            user_db = UserDatabase(
                name=row.name,
                path=row.genome.path,
                db_type=DATABASE_TYPE_DNA,
                created_at=row.genome.created_at or _now_iso(),
                show_in_annotation=row.show_in_annotation,
                annotation_name=new_ann,
            )
            row.annotation = self._ensure_user_db_tracked(user_db)
        self._table.blockSignals(True)
        item.setText(new_ann if new_ann else "-")
        item.setData(Qt.ItemDataRole.UserRole, new_ann)
        item.setToolTip(
            f"{new_ann or '(empty)'}\nDouble-click to edit "
            f"(max {ANNOTATION_NAME_MAX_LENGTH})"
        )
        self._table.blockSignals(False)
        self._persist_all()

    def _revert_name_cell(self, item: QTableWidgetItem, name: str) -> None:
        self._table.blockSignals(True)
        item.setText(name)
        item.setData(Qt.ItemDataRole.UserRole, name)
        self._table.blockSignals(False)

    def _revert_annotation_name_cell(
        self, item: QTableWidgetItem, ann_name: str
    ) -> None:
        self._table.blockSignals(True)
        item.setText(ann_name if ann_name else "-")
        item.setData(Qt.ItemDataRole.UserRole, ann_name)
        self._table.blockSignals(False)

    def _persist_all(self) -> None:
        from repeatviewer.views.settings_dialog import save_workspace

        save_database_settings(
            cdd_dirs=self._project.cdd_dirs,
            user_databases=self._project.user_databases,
        )
        save_blast_databases_to_settings(self._project.blast_databases)
        save_managed_genomes_to_settings(self._project.managed_genomes)
        save_workspace(self._project)
        self.references_changed.emit()

    def _add_database(self) -> None:
        if self._add_dialog is not None and self._add_dialog.isVisible():
            self._add_dialog.raise_()
            self._add_dialog.activateWindow()
            return
        section = self._current_section()
        if section == "blast":
            initial = DB_KIND_FASTA_DNA
        elif section == "hmm":
            initial = DB_KIND_HMM_PROTEIN
        else:
            initial = DB_KIND_CDD
        dialog = AddDatabaseDialog(
            self._project, initial_kind=initial, section=section, parent=self
        )
        self._add_dialog = dialog
        dialog.finished.connect(self._on_add_finished)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _on_add_finished(self, result: int) -> None:
        dialog = self._add_dialog
        self._add_dialog = None
        if dialog is None or result != int(QDialog.DialogCode.Accepted):
            return
        dest_tab = 0
        if dialog.result_annotation() is not None:
            db = dialog.result_annotation()
            assert db is not None
            self._project.user_databases.append(db)
            if db.db_type == DATABASE_TYPE_CDD:
                self._project.cdd_dirs = [db.path]
                dest_tab = 2
            elif db.is_hmm:
                dest_tab = 1
            else:
                dest_tab = 0
            self._maybe_prepare_annotation(db)
        if dialog.result_blast() is not None:
            db = dialog.result_blast()
            assert db is not None
            self._project.blast_databases.append(db)
            dest_tab = 0
        if dialog.result_genome() is not None:
            g = dialog.result_genome()
            assert g is not None
            key = _genome_path_key(g.path)
            self._project.managed_genomes = [
                x for x in self._project.managed_genomes if _genome_path_key(x.path) != key
            ]
            self._project.managed_genomes.insert(0, g)
            if dialog.set_as_main():
                self._apply_main_genome(g.path, name=g.name)
            dest_tab = 0
        self._persist_all()
        # Show the tab that owns the new entry (Type decides destination).
        if self._section_bar.currentIndex() != dest_tab:
            self._section_bar.setCurrentIndex(dest_tab)
        else:
            self._populate()

    def _maybe_prepare_annotation(self, db: UserDatabase) -> None:
        targets = collect_prepare_targets(
            cdd_dir=db.path if db.db_type == DATABASE_TYPE_CDD else None,
            user_databases=[] if db.db_type == DATABASE_TYPE_CDD else [db],
        )
        pending = [t for t in targets if t.needs_index]
        if not pending:
            return
        if any(t.kind == "hmm" for t in pending) and not hmm_indexing_available():
            QMessageBox.warning(
                self,
                "HMM database",
                "An HMM database needs either:\n\n"
                "  • HMMER with hmmpress (conda install -c bioconda hmmer), or\n"
                f"  • pyhmmer ({pyhmmer_install_hint()})",
            )
            return
        host = self.parentWidget() or self
        runner = find_database_index_runner(host)
        use_system = self._project.use_system_blast
        QTimer.singleShot(
            0,
            lambda: runner.start(pending, use_system_blast=use_system, show_success=True),
        )

    def _remove_selected(self) -> None:
        selected = self._selected_rows()
        if not selected:
            return
        names = ", ".join(r.name for r in selected[:5])
        if len(selected) > 5:
            names += f", … (+{len(selected) - 5})"
        confirm = QMessageBox.question(
            self,
            "Manage Database",
            f"Remove {len(selected)} database(s)?\n\n{names}",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if confirm != QMessageBox.StandardButton.Yes:
            return

        genome_keys = {
            _genome_path_key(r.genome.path) for r in selected if r.genome is not None
        }
        # FASTA DNA BLAST set as main also tracks managed genomes by path.
        for r in selected:
            if r.blast is not None and r.blast.molecule_type == BLAST_MOLECULE_NUCL:
                genome_keys.add(_genome_path_key(r.blast.source_fasta))
        ann_names = {r.annotation.name for r in selected if r.annotation is not None}
        blast_names = {r.blast.name for r in selected if r.blast is not None}

        if genome_keys:
            self._project.managed_genomes = [
                g
                for g in self._project.managed_genomes
                if _genome_path_key(g.path) not in genome_keys
            ]
            if self._project.genome_file is not None:
                if _genome_path_key(self._project.genome_file) in genome_keys:
                    self._project.genome_file = None

        if ann_names:
            removed_cdd = any(
                r.annotation is not None
                and r.annotation.db_type == DATABASE_TYPE_CDD
                and r.annotation.name in ann_names
                for r in selected
            )
            self._project.user_databases = [
                db for db in self._project.user_databases if db.name not in ann_names
            ]
            if removed_cdd:
                self._project.cdd_dirs = []

        if blast_names:
            for r in selected:
                if r.blast is None or r.blast.name not in blast_names:
                    continue
                try:
                    prefix = r.blast.db_prefix_path.expanduser()
                    db_dir = prefix.parent
                    if db_dir.is_dir() and db_dir.name == (
                        sanitize_name(r.blast.name) or db_dir.name
                    ):
                        shutil.rmtree(db_dir, ignore_errors=True)
                except OSError:
                    pass
            self._project.blast_databases = [
                db for db in self._project.blast_databases if db.name not in blast_names
            ]

        self._persist_all()
        self._populate()

    def _fasta_path_for_row(self, row: _RefRow) -> Path | None:
        if row.genome is not None:
            return row.genome.path
        if row.blast is not None and row.kind == DB_KIND_FASTA_DNA:
            return row.blast.source_fasta
        if row.annotation is not None and row.kind == DB_KIND_FASTA_DNA:
            return row.annotation.path
        return None

    def _apply_main_genome(self, path: Path, *, name: str | None = None) -> None:
        from repeatviewer.views.settings_dialog import add_genome_to_history

        resolved = path.expanduser().resolve()
        if not resolved.is_file():
            QMessageBox.warning(
                self, "Manage Database", f"Genome file not found:\n{path}"
            )
            return
        _ensure_chrom_sizes(resolved)
        self._project.genome_file = resolved
        self._project.genome_mode = GENOME_MODE_SINGLE
        add_genome_to_history(resolved)
        key = _genome_path_key(resolved)
        others = [
            x for x in self._project.managed_genomes if _genome_path_key(x.path) != key
        ]
        g = ManagedGenome(
            name=name or resolved.stem,
            path=resolved,
            created_at=_now_iso(),
        )
        for existing in self._project.managed_genomes:
            if _genome_path_key(existing.path) == key and existing.created_at:
                g.created_at = existing.created_at
                g.name = existing.name or g.name
                break
        self._project.managed_genomes = [g, *others]

    def _set_main_reference(self) -> None:
        selected = self._selected_rows()
        if len(selected) != 1:
            return
        path = self._fasta_path_for_row(selected[0])
        if path is None:
            return
        self._apply_main_genome(path, name=selected[0].name)
        self._persist_all()
        self._populate()


def open_manage_references(
    project: RepeatViewerProject,
    parent: QWidget | None = None,
    *,
    role_filter: RoleFilter = "all",
    existing: ManageReferencesDialog | None = None,
) -> ManageReferencesDialog:
    """Show or raise Manage Database."""
    if existing is not None and existing.isVisible():
        existing.apply_role_filter(role_filter)
        existing.raise_()
        existing.activateWindow()
        return existing
    dialog = ManageReferencesDialog(project, parent, role_filter=role_filter)
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()
    return dialog


class BlastDatabaseManagerDialog(ManageReferencesDialog):
    """BLAST Manage button → Manage Database (BLAST tab)."""

    def __init__(
        self, project: RepeatViewerProject, parent: QWidget | None = None
    ) -> None:
        super().__init__(project, parent, role_filter="blast")
