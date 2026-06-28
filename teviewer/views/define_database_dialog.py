"""Define CDD and user similarity-search databases (File → Define database)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, QTimer, Signal
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
    QVBoxLayout,
    QWidget,
)

from teviewer.backends.database_prepare import (
    collect_prepare_targets,
    describe_database_problem,
)
from teviewer.project import (
    DATABASE_TYPE_DNA,
    DATABASE_TYPE_PROTEIN,
    TEviewerProject,
    USER_DATABASE_TYPES,
    UserDatabase,
)
from teviewer.views.settings_dialog import _path_browse_row, _section_label, save_workspace
from teviewer.widgets.database_index_runner import find_database_index_runner

_MAX_USER_DATABASE_NAME_LENGTH = 10


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
        self._name.setFixedWidth(120)
        self._name.setToolTip(f"Up to {_MAX_USER_DATABASE_NAME_LENGTH} characters.")
        self._path = QLineEdit(self)
        self._path.setPlaceholderText("Database file path")
        self._type = QComboBox(self)
        self._type.setFixedWidth(88)
        for db_type in USER_DATABASE_TYPES:
            self._type.addItem(db_type)
        self._hmm = QCheckBox("HMM", self)
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
        self._remove_btn.setVisible(removable)
        self._remove_btn.clicked.connect(self.remove_requested.emit)

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
        layout.addWidget(self._hmm)
        layout.addWidget(self._remove_btn)

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
        self._rows.append(row)
        self._rows_layout.addWidget(row)
        self._refresh_remove_buttons()

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
        self.setMinimumWidth(680)
        self._project = project

        self._cdd_path = QLineEdit()
        self._cdd_path.setPlaceholderText("CDD database directory")
        if project.cdd_dir is not None:
            self._cdd_path.setText(str(project.cdd_dir))

        self._user_databases = UserDatabasesWidget(
            self,
            initial=list(project.user_databases),
        )

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

    def accept(self) -> None:
        pending = self._pending_database_index_targets()
        if pending is None:
            return
        use_system_blast = self._project.use_system_blast
        self._apply_to_project()
        host = self.parent() if self.parent() else self
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

    def _collect_database_targets(self):
        cdd_text = self._cdd_path.text().strip()
        cdd_dir = Path(cdd_text).expanduser().resolve() if cdd_text else None
        return collect_prepare_targets(
            cdd_dir=cdd_dir,
            user_databases=self._user_databases.databases(),
        )

    def _pending_database_index_targets(self) -> list | None:
        targets = self._collect_database_targets()
        if not targets:
            return []

        problems: list[str] = []
        for target in targets:
            problem = describe_database_problem(target.kind, target.path)
            if problem:
                problems.append(f"{target.label}\n{problem}")
        if problems:
            QMessageBox.warning(
                self,
                "Database path",
                "\n\n".join(problems),
            )
            return None

        return [target for target in targets if target.needs_index]

    def _apply_to_project(self) -> TEviewerProject:
        cdd_text = self._cdd_path.text().strip()
        self._project.cdd_dirs = (
            [Path(cdd_text).expanduser().resolve()] if cdd_text else []
        )
        self._project.user_databases = self._user_databases.databases()
        save_workspace(self._project)
        return self._project

    def project(self) -> TEviewerProject:
        return self._project
