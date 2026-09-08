"""Parameters dialog for TE benchmarking."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSettings
from PySide6.QtWidgets import (
    QCheckBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

SETTINGS_GROUP = "plugin/te_benchmark"


class BenchmarkParamsDialog(QDialog):
    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("TE benchmarking")
        self.setMinimumWidth(520)
        self._settings = QSettings()

        self._tetrimmer_home = QLineEdit()
        self._reference = QLineEdit()
        self._less_class = QCheckBox("Simpler TE classes (merge LINE/SINE to Non-LTR)")
        self._debug = QCheckBox("Keep intermediate files (debug)")

        self._load_saved()

        form = QFormLayout()
        form.addRow("TEtrimmer folder:", self._dir_row(self._tetrimmer_home, "Select TEtrimmer repository root"))
        form.addRow(
            "Reference RM .out:",
            self._file_row(
                self._reference,
                "RepeatMasker .out from curated TE library",
                "RepeatMasker out (*.out);;All files (*)",
            ),
        )
        form.addRow(self._less_class)
        form.addRow(self._debug)

        hint = QLabel(
            "Select one or more <b>query</b> RepeatMasker .out files in the Files panel. "
            "The reference .out is the gold-standard library annotation (same for all queries). "
            "Genome FASTA comes from RepeatViewer Project settings."
        )
        hint.setWordWrap(True)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addWidget(hint)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def _load_saved(self) -> None:
        self._settings.beginGroup(SETTINGS_GROUP)
        self._tetrimmer_home.setText(self._settings.value("tetrimmer_home", ""))
        self._reference.setText(self._settings.value("reference_out", ""))
        self._less_class.setChecked(bool(self._settings.value("less_classification", False)))
        self._debug.setChecked(bool(self._settings.value("debug", False)))
        self._settings.endGroup()

    def _save(self) -> None:
        self._settings.beginGroup(SETTINGS_GROUP)
        self._settings.setValue("tetrimmer_home", self._tetrimmer_home.text().strip())
        self._settings.setValue("reference_out", self._reference.text().strip())
        self._settings.setValue("less_classification", self._less_class.isChecked())
        self._settings.setValue("debug", self._debug.isChecked())
        self._settings.endGroup()

    def accept(self) -> None:
        self._save()
        super().accept()

    def tetrimmer_home(self) -> Path | None:
        text = self._tetrimmer_home.text().strip()
        return Path(text).expanduser() if text else None

    def reference_out(self) -> Path | None:
        text = self._reference.text().strip()
        return Path(text).expanduser() if text else None

    def less_classification(self) -> bool:
        return self._less_class.isChecked()

    def debug(self) -> bool:
        return self._debug.isChecked()

    @staticmethod
    def _dir_row(line: QLineEdit, title: str) -> QWidget:
        row = QWidget()
        h = QHBoxLayout(row)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(line)

        def browse() -> None:
            chosen = QFileDialog.getExistingDirectory(
                row,
                title,
                line.text() or str(Path.home()),
            )
            if chosen:
                line.setText(chosen)

        btn = QPushButton("Browse…")
        btn.clicked.connect(browse)
        h.addWidget(btn)
        return row

    @staticmethod
    def _file_row(line: QLineEdit, title: str, filt: str) -> QWidget:
        row = QWidget()
        h = QHBoxLayout(row)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(line)

        def browse() -> None:
            chosen, _ = QFileDialog.getOpenFileName(
                row,
                title,
                line.text() or str(Path.home()),
                filt,
            )
            if chosen:
                line.setText(chosen)

        btn = QPushButton("Browse…")
        btn.clicked.connect(browse)
        h.addWidget(btn)
        return row
