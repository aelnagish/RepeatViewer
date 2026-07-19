"""Define CDD and user similarity-search databases (File → Define database)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QEvent, Qt, QTimer, Signal
from PySide6.QtWidgets import (
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from teviewer.backends.database_prepare import (
    CDD_MARKER,
    collect_prepare_targets,
    describe_database_problem,
    hmm_indexing_available,
    index_is_ready,
    pyhmmer_install_hint,
)
from teviewer.backends.pyhmmer_search import hmm_uses_nucleotide_query
from teviewer.project import (
    DATABASE_TYPE_DNA,
    DATABASE_TYPE_PROTEIN,
    TEviewerProject,
    USER_DATABASE_TYPES,
    UserDatabase,
)
from teviewer.views.settings_dialog import (
    _path_browse_row,
    _section_label,
    save_workspace,
)
from teviewer.widgets.database_index_runner import find_database_index_runner
_BLAST_DNA_INDEX_SUFFIXES = frozenset(
    {".nhr", ".nin", ".nsq", ".nog", ".nos", ".nod", ".not", ".ndb", ".nto", ".ntf", ".nni"}
)
_BLAST_PROT_INDEX_SUFFIXES = frozenset(
    {".phr", ".pin", ".psq", ".pog", ".pos", ".pod", ".pot", ".pdb", ".pto", ".ptf", ".pni"}
)
_HMM_INDEX_SUFFIXES = frozenset({".h3m", ".h3f", ".h3i", ".h3p"})
_DATABASE_FASTA_SUFFIXES = frozenset({".fa", ".fasta", ".faa", ".fna", ".fas", ".hmm"})
_MAX_USER_DATABASE_NAME_LENGTH = 50


def resolve_panel_database_path(path: Path) -> tuple[str, Path] | None:
    """Map a Files-panel entry to a CDD directory or user database file path."""
    path = path.expanduser().resolve()
    if not path.is_file():
        return None
    if path.name.lower() == CDD_MARKER.lower():
        return ("cdd", path.parent)
    suffix = path.suffix.lower()
    if suffix == ".hmm":
        return ("user", path)
    if suffix in _HMM_INDEX_SUFFIXES:
        base = path.with_suffix("")
        if base.suffix.lower() == ".hmm":
            return ("user", base)
        return None
    if suffix in _BLAST_DNA_INDEX_SUFFIXES | _BLAST_PROT_INDEX_SUFFIXES:
        return ("user", path.with_suffix(""))
    if suffix in _DATABASE_FASTA_SUFFIXES:
        return ("user", path)
    return None


def infer_user_database_defaults(path: Path) -> tuple[str, str, bool]:
    """Derive display name, database type, and HMM flag from a file path."""
    suffix = path.suffix.lower()
    name = path.stem[:_MAX_USER_DATABASE_NAME_LENGTH]
    if suffix == ".hmm":
        db_type = (
            DATABASE_TYPE_DNA
            if hmm_uses_nucleotide_query(path)
            else DATABASE_TYPE_PROTEIN
        )
        return name, db_type, True
    if suffix in _HMM_INDEX_SUFFIXES and path.with_suffix("").suffix.lower() == ".hmm":
        hmm_path = path.with_suffix("")
        db_type = (
            DATABASE_TYPE_DNA
            if hmm_uses_nucleotide_query(hmm_path)
            else DATABASE_TYPE_PROTEIN
        )
        return hmm_path.stem[:_MAX_USER_DATABASE_NAME_LENGTH], db_type, True
    if suffix in _BLAST_DNA_INDEX_SUFFIXES or suffix == ".fna":
        return name, DATABASE_TYPE_DNA, False
    if suffix in _BLAST_PROT_INDEX_SUFFIXES or suffix == ".faa":
        return name, DATABASE_TYPE_PROTEIN, False
    return name, DATABASE_TYPE_PROTEIN, False


class _UserDatabaseRow(QWidget):
    remove_requested = Signal()

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        initial: UserDatabase | None = None,
        removable: bool = True,
    ) -> None:
        super().__init__(parent)
        self._name = QLineEdit(self)
        self._name.setPlaceholderText("Database name")
        self._name.setMaxLength(_MAX_USER_DATABASE_NAME_LENGTH)
        self._name.setFixedWidth(96)
        self._name.setToolTip(f"Up to {_MAX_USER_DATABASE_NAME_LENGTH} characters.")
        self._path = QLineEdit(self)
        self._path.setPlaceholderText("Database file path")
        self._path.setMinimumWidth(140)
        self._type = QComboBox(self)
        self._type.setFixedWidth(88)
        for db_type in USER_DATABASE_TYPES:
            self._type.addItem(db_type)
        self._hmm = QCheckBox("HMM", self)
        self._hmm.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._hmm.setMinimumWidth(self._hmm.sizeHint().width() + 4)
        self._hmm.setToolTip(
            "Profile HMM model (hmmpress / pyhmmer). "
            "Protein HMM searches ORF amino acids; DNA HMM searches the nucleotide sequence."
        )
        if isinstance(initial, UserDatabase):
            self._name.setText(initial.name)
            self._path.setText(str(initial.path))
            type_index = self._type.findText(initial.db_type)
            if type_index >= 0:
                self._type.setCurrentIndex(type_index)
            self._hmm.setChecked(initial.is_hmm)

        browse_btn = QPushButton("Browse", self)
        browse_btn.clicked.connect(self._browse)
        self._remove_btn = QPushButton("Remove", self)
        self._remove_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._remove_btn.setVisible(removable)
        self._remove_btn.clicked.connect(self.remove_requested.emit)
        self._path.editingFinished.connect(self._apply_inferred_defaults)

        path_container = QWidget(self)
        path_layout = QHBoxLayout(path_container)
        path_layout.setContentsMargins(0, 0, 0, 0)
        path_layout.setSpacing(6)
        path_layout.addWidget(self._path, stretch=1)
        path_layout.addWidget(browse_btn)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        layout.addWidget(QLabel("Name", self))
        layout.addWidget(self._name)
        layout.addWidget(QLabel("Path", self))
        layout.addWidget(path_container, stretch=1)
        layout.addWidget(QLabel("Type", self))
        layout.addWidget(self._type)
        layout.addWidget(self._hmm, 0, Qt.AlignmentFlag.AlignVCenter)
        layout.addWidget(self._remove_btn, 0, Qt.AlignmentFlag.AlignVCenter)

    def _browse(self) -> None:
        start = self._path.text().strip() or str(Path.home())
        if self._hmm.isChecked():
            file_filter = "HMM profiles (*.hmm *.h3m);;All files (*)"
        elif self._type.currentText() == DATABASE_TYPE_DNA:
            file_filter = "DNA/FASTA (*.fa *.fasta *.fna);;All files (*)"
        else:
            file_filter = "Protein/FASTA (*.fa *.fasta *.faa);;All files (*)"
        picked, _ = QFileDialog.getOpenFileName(
            self,
            "Select database file",
            start,
            file_filter,
        )
        if picked:
            self._path.setText(picked)
            self._apply_inferred_defaults()

    def _apply_inferred_defaults(self) -> None:
        text = self._path.text().strip()
        if not text:
            return
        path = Path(text).expanduser()
        if not path.is_file():
            return
        name, db_type, is_hmm = infer_user_database_defaults(path)
        if not self._name.text().strip():
            self._name.setText(name)
        type_index = self._type.findText(db_type)
        if type_index >= 0:
            self._type.setCurrentIndex(type_index)
        if path.suffix.lower() == ".hmm" or any(
            path.name.endswith(suffix) for suffix in _HMM_INDEX_SUFFIXES
        ):
            self._hmm.setChecked(True)
        else:
            self._hmm.setChecked(is_hmm)

    def apply_path_from_panel(self, path: Path) -> None:
        path = path.expanduser().resolve()
        self._path.setText(str(path))
        name, db_type, is_hmm = infer_user_database_defaults(path)
        if not self._name.text().strip():
            self._name.setText(name)
        type_index = self._type.findText(db_type)
        if type_index >= 0:
            self._type.setCurrentIndex(type_index)
        self._hmm.setChecked(is_hmm)

    def database(self) -> UserDatabase | None:
        name = self._name.text().strip()
        path = self._path.text().strip()
        if not name or not path:
            return None
        return UserDatabase(
            name=name,
            path=Path(path).expanduser(),
            db_type=self._type.currentText(),
            is_hmm=self._hmm.isChecked(),
        )


class UserDatabasesWidget(QWidget):
    """Named user databases (Protein / DNA, optional HMM profile) with add/remove controls."""

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        initial: list[UserDatabase] | None = None,
    ) -> None:
        super().__init__(parent)
        self._rows: list[_UserDatabaseRow] = []
        self._active_row: _UserDatabaseRow | None = None

        self._rows_layout = QVBoxLayout()
        self._rows_layout.setContentsMargins(0, 0, 0, 0)
        self._rows_layout.setSpacing(6)

        self._add_btn = QPushButton("Add another database", self)
        self._add_btn.clicked.connect(self._on_add_row)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(6)
        outer.addLayout(self._rows_layout)
        outer.addWidget(self._add_btn, 0, Qt.AlignmentFlag.AlignLeft)

        if initial:
            for entry in initial:
                self.add_row(entry)
        if not self._rows:
            self.add_row()

    def _on_add_row(self) -> None:
        self.add_row()

    def add_row(self, entry: UserDatabase | None = None) -> None:
        row = _UserDatabaseRow(self, initial=entry, removable=len(self._rows) >= 1)
        row.remove_requested.connect(lambda r=row: self._remove_row(r))
        row._path.installEventFilter(self)
        self._rows.append(row)
        self._rows_layout.addWidget(row)
        self._refresh_remove_buttons()

    def eventFilter(self, watched, event) -> bool:  # noqa: N802 (Qt API)
        if event.type() == QEvent.Type.FocusIn:
            for row in self._rows:
                if watched is row._path:
                    self._active_row = row
                    break
        return super().eventFilter(watched, event)

    def apply_path_from_panel(self, path: Path) -> bool:
        path = path.expanduser().resolve()
        if not path.is_file():
            return False
        target = self._active_row
        if target is None:
            return False
        target.apply_path_from_panel(path)
        return True

    def _remove_row(self, row: _UserDatabaseRow) -> None:
        if row not in self._rows or len(self._rows) <= 1:
            return
        self._rows.remove(row)
        self._rows_layout.removeWidget(row)
        row.deleteLater()
        self._refresh_remove_buttons()

    def _refresh_remove_buttons(self) -> None:
        removable = len(self._rows) > 1
        for row in self._rows:
            row._remove_btn.setVisible(removable)

    def databases(self) -> list[UserDatabase]:
        out: list[UserDatabase] = []
        for row in self._rows:
            entry = row.database()
            if entry is not None:
                out.append(entry)
        return out


class DefineDatabaseDialog(QDialog):
    """Configure CDD and custom similarity annotation databases."""

    def __init__(self, project: TEviewerProject, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Define databases")
        self.setMinimumWidth(780)
        self.resize(820, self.sizeHint().height())
        self._project = project

        self._cdd_path = QLineEdit()
        self._cdd_path.setPlaceholderText("CDD database directory")
        if project.cdd_dir is not None:
            self._cdd_path.setText(str(project.cdd_dir))

        self._user_databases = UserDatabasesWidget(
            self,
            initial=list(project.user_databases),
        )
        self._focused_path_edit: QLineEdit | None = None
        self._cdd_path.installEventFilter(self)

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setFormAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setHorizontalSpacing(8)
        form.setVerticalSpacing(8)
        form.addRow(
            "CDD database dir:",
            _path_browse_row(
                self._cdd_path,
                parent=self,
                is_dir=True,
                browse_title="Select CDD database directory",
            ),
        )

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(10)
        layout.addLayout(form)
        layout.addWidget(_section_label("Other databases", self))
        layout.addWidget(self._user_databases)
        layout.addWidget(buttons)

    def showEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().showEvent(event)
        self._focused_path_edit = self._cdd_path
        self._cdd_path.setFocus(Qt.FocusReason.OtherFocusReason)

    def eventFilter(self, watched, event) -> bool:  # noqa: N802 (Qt API)
        if event.type() in {QEvent.Type.FocusIn, QEvent.Type.MouseButtonPress}:
            if watched is self._cdd_path:
                self._focused_path_edit = self._cdd_path
            else:
                for row in self._user_databases._rows:
                    if watched is row._path:
                        self._focused_path_edit = row._path
                        self._user_databases._active_row = row
                        break
        return super().eventFilter(watched, event)

    def apply_folder_from_panel(self, folder: Path) -> bool:
        """Fill the focused CDD directory field from a Folders-panel click."""
        if self._focused_path_edit is not self._cdd_path:
            return False
        if self.focusWidget() is not self._cdd_path:
            return False
        resolved = folder.expanduser().resolve()
        if not resolved.is_dir():
            return False
        self._cdd_path.setText(str(resolved))
        return True

    def apply_file_from_panel(self, path: Path) -> bool:
        """Fill the focused CDD or user-database path field from a Files panel click."""
        if self._focused_path_edit is None:
            return False
        resolved = resolve_panel_database_path(path)
        if resolved is None:
            return False
        kind, target_path = resolved
        if self._focused_path_edit is self._cdd_path:
            if kind != "cdd":
                return False
            self._cdd_path.setText(str(target_path))
            return True
        for row in self._user_databases._rows:
            if self._focused_path_edit is row._path:
                if kind != "user":
                    return False
                row.apply_path_from_panel(target_path)
                return True
        return False

    def accept(self) -> None:
        pending, blocked = self._pending_database_index_targets()
        if blocked:
            return
        use_system_blast = self._project.use_system_blast
        saved_targets = self._collect_database_targets()
        self._apply_to_project()
        host = self.parentWidget() or self.parent() or self
        super().accept()
        if pending:
            runner = find_database_index_runner(host)
            QTimer.singleShot(
                0,
                lambda: runner.start(
                    pending,
                    use_system_blast=use_system_blast,
                    show_success=True,
                ),
            )
        elif saved_targets:
            already_ready = [
                target.label
                for target in saved_targets
                if index_is_ready(target.kind, target.path)
            ]
            if already_ready:
                QMessageBox.information(
                    self,
                    "Define databases",
                    "Database settings saved.\n\n"
                    "No indexing was needed (index files already present):\n"
                    + "\n".join(f"  • {label}" for label in already_ready),
                )

    def _collect_database_targets(self):
        cdd_text = self._cdd_path.text().strip()
        cdd_dir = Path(cdd_text).expanduser().resolve() if cdd_text else None
        return collect_prepare_targets(
            cdd_dir=cdd_dir,
            user_databases=self._user_databases.databases(),
        )

    def _pending_database_index_targets(self) -> tuple[list, bool]:
        targets = self._collect_database_targets()
        if not targets:
            return [], False

        valid: list = []
        problems: list[str] = []
        for target in targets:
            problem = describe_database_problem(target.kind, target.path)
            if problem:
                problems.append(f"{target.label}\n{problem}")
            else:
                valid.append(target)

        if problems:
            extra = (
                "\n\nOther valid databases will still be saved and indexed."
                if valid
                else ""
            )
            QMessageBox.warning(
                self,
                "Database path",
                "\n\n".join(problems) + extra,
            )
            if not valid:
                return [], True

        pending = [target for target in valid if target.needs_index]
        if pending and any(target.kind == "hmm" for target in pending):
            if not hmm_indexing_available():
                QMessageBox.warning(
                    self,
                    "HMM database",
                    "An HMM database (e.g. Dfam) needs either:\n\n"
                    "  • HMMER with hmmpress (conda install -c bioconda hmmer), or\n"
                    f"  • pyhmmer ({pyhmmer_install_hint()})\n\n"
                    "Check the HMM box for .hmm files, then try again.",
                )
                return [], True
        return pending, False

    def _apply_to_project(self) -> TEviewerProject:
        cdd_text = self._cdd_path.text().strip()
        if cdd_text:
            cdd_dir = Path(cdd_text).expanduser().resolve()
            if describe_database_problem("cdd", cdd_dir) is None:
                self._project.cdd_dirs = [cdd_dir]
        else:
            self._project.cdd_dirs = []
        self._project.user_databases = self._user_databases.databases()
        save_workspace(self._project)
        return self._project

    def project(self) -> TEviewerProject:
        return self._project
