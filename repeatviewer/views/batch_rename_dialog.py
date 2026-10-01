"""Edit > Batch Rename: rename many files at once, with a live preview.

Find / Replace With (plain text or a regular expression), Add Prefix, Add
Suffix and Keep File Extension. The preview lists every file with its new name
and marks names that cannot be used (taken, duplicated, invalid); Rename stays
off until every new name is fine. The main window performs the rename
(``MainWindow._batch_rename_files``), all at once and as one Undo step.
"""

from __future__ import annotations

import os
from pathlib import Path

from PySide6.QtCore import QSettings, Qt, QTimer
from PySide6.QtGui import QBrush, QColor, QFont, QPalette
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QDialog,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QPushButton,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.batch_rename import RenameOptions, RenamePlan, plan_batch_rename
from repeatviewer.widgets.help_button import HELP_BUTTON_SIZE, make_help_button

_SETTINGS = "batch_rename"
_PROBLEM_COLOR = QColor("#c62828")

FIND_HELP = (
    "<b>Find</b>: text to look for in each file name; every match is replaced.<br>"
    "<b>Replace With</b>: the text that takes its place. Leave it empty to remove "
    "the found text.<br>"
    "Example: Find <b>contig</b>, Replace With <b>Chr</b> turns contig_1.fa into Chr_1.fa."
)
MATCH_HELP = (
    "<b>Match Case</b>: Find only matches the same upper and lower case letters.<br>"
    "<b>Regular Expression</b>: Find is a Python regular expression, e.g. "
    "<b>^seq_(\\d+)$</b>. Use <b>\\1</b>, <b>\\2</b> in Replace With to reuse the "
    "parts in brackets."
)
PREFIX_HELP = "Text added at the start of every name, e.g. <b>CroSat_</b>."
SUFFIX_HELP = (
    "Text added at the end of every name, e.g. <b>_v2</b>.<br>"
    "With Keep File Extension it goes before the extension: sample_v2.fa."
)
EXTENSION_HELP = (
    "<b>On</b>: the extension (.fa, .fasta.gz, .gff3) is left out of Find and Replace "
    "and stays as it is.<br>"
    "<b>Off</b>: the whole name, extension included, can change."
)


def _with_help(widget: QWidget, title: str, text: str, parent: QWidget) -> QWidget:
    wrap = QWidget(parent)
    row = QHBoxLayout(wrap)
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(6)
    row.addWidget(widget, 1)
    row.addWidget(make_help_button(wrap, title, text))
    return wrap


def _aligned(widget: QWidget, parent: QWidget) -> QWidget:
    """``widget`` with room on the right where other rows have their "?" button."""
    wrap = QWidget(parent)
    row = QHBoxLayout(wrap)
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(6)
    row.addWidget(widget, 1)
    row.addSpacing(HELP_BUTTON_SIZE + row.spacing())  # no spacing next to a spacer
    return wrap


def _companions_for(paths: list[Path]) -> dict[str, list[Path]]:
    """Annotation sidecars (.seqmap_RV.json) that follow each FASTA."""
    from repeatviewer.backends.sequence_features import companion_annotation_paths_for_fasta

    companions: dict[str, list[Path]] = {}
    for path in paths:
        try:
            found = companion_annotation_paths_for_fasta(path)
        except OSError:
            found = []
        if found:
            companions[str(path)] = found
    return companions


class BatchRenameDialog(QDialog):
    """Options plus a preview table; Cancel / Rename."""

    def __init__(self, paths: list[Path], parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Batch Rename")
        self.setModal(True)
        self._paths = [Path(p) for p in paths]
        self._companions = _companions_for(self._paths)
        self._plan = RenamePlan([])

        folders = {p.parent for p in self._paths}
        count = f"{len(self._paths)} file{'s' if len(self._paths) != 1 else ''}"
        where = f" in {next(iter(folders)).name}" if len(folders) == 1 else ""
        heading = QLabel(f"Rename {count}{where}", self)
        font = QFont(heading.font())
        font.setBold(True)
        heading.setFont(font)

        self._find = QLineEdit(self)
        self._replace = QLineEdit(self)
        self._prefix = QLineEdit(self)
        self._suffix = QLineEdit(self)
        for field in (self._find, self._replace, self._prefix, self._suffix):
            field.setClearButtonEnabled(True)
            field.setMinimumWidth(300)
        self._match_case = QCheckBox("Match Case", self)
        self._regex = QCheckBox("Regular Expression", self)
        self._keep_extension = QCheckBox("Keep File Extension", self)

        form = QGridLayout()
        form.setHorizontalSpacing(10)
        form.setVerticalSpacing(6)
        form.addWidget(QLabel("Find:", self), 0, 0)
        form.addWidget(_with_help(self._find, "Find and Replace", FIND_HELP, self), 0, 1)
        form.addWidget(QLabel("Replace With:", self), 1, 0)
        form.addWidget(_aligned(self._replace, self), 1, 1)
        match_row = QWidget(self)
        match_layout = QHBoxLayout(match_row)
        match_layout.setContentsMargins(0, 0, 0, 0)
        match_layout.addWidget(self._match_case)
        match_layout.addSpacing(12)
        match_layout.addWidget(self._regex)
        match_layout.addWidget(make_help_button(match_row, "Match Case and Regular Expression", MATCH_HELP))
        match_layout.addStretch(1)
        form.addWidget(match_row, 2, 1)
        form.addWidget(QLabel("Add Prefix:", self), 3, 0)
        form.addWidget(_with_help(self._prefix, "Add Prefix", PREFIX_HELP, self), 3, 1)
        form.addWidget(QLabel("Add Suffix:", self), 4, 0)
        form.addWidget(_with_help(self._suffix, "Add Suffix", SUFFIX_HELP, self), 4, 1)
        ext_row = QWidget(self)
        ext_layout = QHBoxLayout(ext_row)
        ext_layout.setContentsMargins(0, 0, 0, 0)
        ext_layout.addWidget(self._keep_extension)
        ext_layout.addWidget(make_help_button(ext_row, "Keep File Extension", EXTENSION_HELP))
        ext_layout.addStretch(1)
        form.addWidget(ext_row, 5, 1)
        form.setColumnStretch(1, 1)

        self._table = QTableWidget(0, 2, self)
        self._table.setHorizontalHeaderLabels(["Current Name", "New Name"])
        self._table.verticalHeader().setVisible(False)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.setWordWrap(False)
        header = self._table.horizontalHeader()
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Stretch)
        header.setSectionResizeMode(1, QHeaderView.ResizeMode.Stretch)
        self._table.setMinimumHeight(180)

        self._status = QLabel("", self)
        self._status.setWordWrap(True)
        self._status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        hint = QLabel(
            "Annotation files (.seqmap_RV.json), file colors and cytogenetic image links "
            "follow their file. Edit > Undo takes the whole rename back.",
            self,
        )
        hint.setWordWrap(True)
        hint.setStyleSheet("color: palette(mid);")

        cancel_btn = QPushButton("Cancel", self)
        cancel_btn.clicked.connect(self.reject)
        self._rename_btn = QPushButton("Rename", self)
        self._rename_btn.setDefault(True)
        self._rename_btn.setEnabled(False)
        self._rename_btn.clicked.connect(self._accept_if_ok)
        buttons = QHBoxLayout()
        buttons.addStretch(1)
        buttons.addWidget(cancel_btn)
        buttons.addWidget(self._rename_btn)

        root = QVBoxLayout(self)
        root.setContentsMargins(22, 18, 22, 16)
        root.setSpacing(10)
        root.addWidget(heading)
        root.addLayout(form)
        root.addWidget(QLabel("Preview:", self))
        root.addWidget(self._table, stretch=1)
        root.addWidget(self._status)
        root.addWidget(hint)
        root.addLayout(buttons)
        self.resize(640, 560)

        # Typing updates the preview after a short pause (big selections stay smooth).
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(150)
        self._timer.timeout.connect(self._update_preview)
        for field in (self._find, self._replace, self._prefix, self._suffix):
            field.textChanged.connect(self._timer.start)
        for box in (self._match_case, self._regex, self._keep_extension):
            box.toggled.connect(self._update_preview)

        self._load_settings()
        self._update_preview()
        self._find.setFocus()

    # ------------------------------------------------------------------ options

    def _load_settings(self) -> None:
        s = QSettings()
        s.beginGroup(_SETTINGS)

        def flag(key: str, default: bool) -> bool:
            return str(s.value(key, default)).lower() in {"1", "true", "yes"}

        for box, key, default in (
            (self._match_case, "match_case", False),
            (self._regex, "regex", False),
            (self._keep_extension, "keep_extension", True),
        ):
            box.blockSignals(True)
            box.setChecked(flag(key, default))
            box.blockSignals(False)
        s.endGroup()

    def _save_settings(self) -> None:
        s = QSettings()
        s.beginGroup(_SETTINGS)
        s.setValue("match_case", self._match_case.isChecked())
        s.setValue("regex", self._regex.isChecked())
        s.setValue("keep_extension", self._keep_extension.isChecked())
        s.endGroup()

    def options(self) -> RenameOptions:
        return RenameOptions(
            find=self._find.text(),
            replace=self._replace.text(),
            match_case=self._match_case.isChecked(),
            regex=self._regex.isChecked(),
            prefix=self._prefix.text(),
            suffix=self._suffix.text(),
            keep_extension=self._keep_extension.isChecked(),
        )

    def plan(self) -> RenamePlan:
        """The rename to perform (as last shown in the preview)."""
        if self._timer.isActive():
            self._timer.stop()
            self._update_preview()
        return self._plan

    # ------------------------------------------------------------------ preview

    def _companion_files(self, path: Path) -> list[Path]:
        return self._companions.get(str(path), [])

    def _companion_renames(self, old: Path, new: Path) -> list[tuple[Path, Path]]:
        from repeatviewer.backends.sequence_features import sidecar_path

        found = self._companion_files(old)
        if not found:
            return []
        # The first one found is the current sidecar; older names stay put.
        source, target = found[0], sidecar_path(new)
        return [] if source == target else [(source, target)]

    def _update_preview(self) -> None:
        self._plan = plan_batch_rename(
            self._paths,
            self.options(),
            companion_files=self._companion_files,
            companion_renames=self._companion_renames,
            listdir=lambda folder: os.listdir(folder) if folder.is_dir() else [],
        )
        self._fill_table()
        self._update_status()

    def _fill_table(self) -> None:
        rows = self._plan.rows
        table = self._table
        table.setUpdatesEnabled(False)
        try:
            table.setRowCount(len(rows))
            muted = QBrush(self.palette().color(QPalette.ColorRole.PlaceholderText))
            problem = QBrush(_PROBLEM_COLOR)
            for index, row in enumerate(rows):
                old_item = QTableWidgetItem(row.source.name)
                old_item.setToolTip(str(row.source))
                new_item = QTableWidgetItem(row.new_name if row.changed else "(no change)")
                if row.problem:
                    new_item.setForeground(problem)
                    new_item.setToolTip(row.problem)
                elif not row.changed:
                    new_item.setForeground(muted)
                table.setItem(index, 0, old_item)
                table.setItem(index, 1, new_item)
        finally:
            table.setUpdatesEnabled(True)

    def _update_status(self) -> None:
        plan = self._plan
        if plan.error:
            self._status.setText(plan.error)
            self._status.setStyleSheet(f"color: {_PROBLEM_COLOR.name()};")
            self._rename_btn.setEnabled(False)
            return
        changed = len(plan.changed_rows)
        problems = plan.problem_rows
        total = len(plan.rows)
        if problems:
            first = problems[0]
            more = f" ({len(problems) - 1} more marked in red)" if len(problems) > 1 else ""
            self._status.setText(f"{first.source.name}: {first.problem}{more}.")
            self._status.setStyleSheet(f"color: {_PROBLEM_COLOR.name()};")
        elif changed:
            noun = "file" if total == 1 else "files"
            self._status.setText(f"{changed} of {total} {noun} will be renamed.")
            self._status.setStyleSheet("")
        else:
            self._status.setText("Type in Find, Add Prefix or Add Suffix to see the new names.")
            self._status.setStyleSheet("")
        self._rename_btn.setEnabled(plan.ok)

    def _accept_if_ok(self) -> None:
        if self.plan().ok:
            self._save_settings()
            self.accept()
