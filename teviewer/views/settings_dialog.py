"""Project path configuration."""

from __future__ import annotations

import json
from collections.abc import Iterable
from pathlib import Path

from PySide6.QtCore import QSettings, Qt
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
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from teviewer.pan_genome import GENOME_MODE_SINGLE, parse_pan_genomes_json
from teviewer.project import TEviewerProject
from teviewer.views.toolbar_dialog import load_toolbar_font_size, save_toolbar_font_size
from teviewer.widgets.browser_panel_settings import (
    load_browser_panel_font_size,
    save_browser_panel_font_size,
)
from teviewer.widgets.preview_zoom import (
    MAX_PREVIEW_ZOOM_PERCENT,
    MIN_PREVIEW_ZOOM_PERCENT,
    PreviewFitMode,
    DefaultFastaPreviewTab,
    load_color_blind_nucleotide_scheme,
    load_default_fasta_preview_tab,
    load_preview_fit_mode,
    load_preview_zoom_percent,
    save_color_blind_nucleotide_scheme,
    save_default_fasta_preview_tab,
    save_preview_fit_mode,
    save_preview_zoom_percent,
)

from teviewer.database_settings import (
    SETTINGS_GROUP,
    load_cdd_dirs_from_settings,
    load_user_databases_from_settings,
    save_database_settings,
)

MAX_GENOME_HISTORY = 10


def _path_from_settings(value) -> Path | None:
    if not value:
        return None
    text = str(value).strip()
    if not text:
        return None
    return Path(text).expanduser()


def _path_list_from_settings(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value] if value.strip() else []
    return [str(v).strip() for v in value if v and str(v).strip()]


def save_workspace(
    project: TEviewerProject,
    *,
    root_folders: Iterable[Path] | None = None,
    last_folder: Path | None = None,
    expanded_folders: Iterable[Path] | None = None,
) -> None:
    """Persist project paths and folder tree workspace (survives app restart)."""
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    s.setValue("input_dir", str(project.input_dir) if project.input_dir else "")
    s.setValue("genome_file", str(project.genome_file) if project.genome_file else "")
    s.setValue("genome_mode", project.genome_mode)
    from teviewer.pan_genome import pan_genomes_to_json

    s.setValue("pan_genomes_json", pan_genomes_to_json(project.pan_genomes))
    save_database_settings(
        cdd_dirs=project.cdd_dirs,
        user_databases=project.user_databases,
        s=s,
    )
    s.setValue("use_system_blast", project.use_system_blast)
    s.setValue("use_system_sequence_viewer", project.use_system_sequence_viewer)
    s.setValue("batch_workers", project.batch_workers)
    if root_folders is not None:
        roots = [str(p.expanduser().resolve()) for p in root_folders if p.is_dir()]
        s.setValue("workspace_roots_json", json.dumps(roots))
        s.setValue("root_folders", roots)
    if last_folder is not None and last_folder.is_dir():
        s.setValue("last_folder", str(last_folder.expanduser().resolve()))
    if expanded_folders is not None:
        expanded = [
            str(p.expanduser().resolve()) for p in expanded_folders if p.is_dir()
        ]
        s.setValue("expanded_folders_json", json.dumps(expanded))
    s.endGroup()


def load_root_folders() -> list[Path]:
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    raw_json = s.value("workspace_roots_json")
    paths: list[Path | None] = []
    if raw_json:
        try:
            decoded = json.loads(str(raw_json))
            if isinstance(decoded, list):
                paths = [_path_from_settings(p) for p in decoded]
        except (json.JSONDecodeError, TypeError):
            paths = []
    if not paths:
        paths = [_path_from_settings(p) for p in _path_list_from_settings(s.value("root_folders"))]
    s.endGroup()
    return [p for p in paths if p is not None and p.is_dir()]


def load_last_folder() -> Path | None:
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    folder = _path_from_settings(s.value("last_folder"))
    s.endGroup()
    if folder is not None and folder.is_dir():
        return folder
    return None


def load_expanded_folders() -> list[Path]:
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    raw = s.value("expanded_folders_json")
    s.endGroup()
    if not raw:
        return []
    try:
        decoded = json.loads(str(raw))
    except (json.JSONDecodeError, TypeError):
        return []
    if not isinstance(decoded, list):
        return []
    paths: list[Path] = []
    seen: set[str] = set()
    for item in decoded:
        path = _path_from_settings(item)
        if path is None or not path.is_dir():
            continue
        key = str(path.expanduser().resolve())
        if key in seen:
            continue
        seen.add(key)
        paths.append(path)
    return paths


def load_genome_history() -> list[Path]:
    """Recently used genome FASTA paths (newest first, existing files only)."""
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    raw = s.value("genome_file_history_json")
    s.endGroup()
    if not raw:
        return []
    try:
        decoded = json.loads(str(raw))
    except (json.JSONDecodeError, TypeError):
        return []
    if not isinstance(decoded, list):
        return []
    paths: list[Path] = []
    seen: set[str] = set()
    for item in decoded:
        path = _path_from_settings(item)
        if path is None:
            continue
        key = str(path.expanduser().resolve())
        if key in seen or not path.is_file():
            continue
        seen.add(key)
        paths.append(path.expanduser().resolve())
    return paths


def save_genome_history(paths: list[Path]) -> None:
    unique: list[str] = []
    seen: set[str] = set()
    for path in paths:
        if not path:
            continue
        try:
            key = str(path.expanduser().resolve())
        except OSError:
            continue
        if key in seen:
            continue
        seen.add(key)
        unique.append(key)
    s = QSettings()
    s.beginGroup(SETTINGS_GROUP)
    s.setValue("genome_file_history_json", json.dumps(unique[:MAX_GENOME_HISTORY]))
    s.endGroup()


def add_genome_to_history(path: Path | None) -> list[Path]:
    if path is None:
        return load_genome_history()
    try:
        resolved = path.expanduser().resolve()
    except OSError:
        return load_genome_history()
    if not resolved.is_file():
        return load_genome_history()
    history = [p for p in load_genome_history() if p != resolved]
    history.insert(0, resolved)
    save_genome_history(history)
    return history


def remove_genome_from_history(path: Path | None) -> list[Path]:
    if path is None:
        return load_genome_history()
    try:
        target = path.expanduser().resolve()
    except OSError:
        return load_genome_history()
    history = [p for p in load_genome_history() if p != target]
    save_genome_history(history)
    return history


def _genome_history_label(path: Path) -> str:
    parent = path.parent.name or str(path.parent)
    return f"{path.name}  —  {parent}"


def _path_browse_row(
    line_edit: QLineEdit,
    *,
    parent: QWidget,
    is_dir: bool,
    browse_title: str,
    file_filter: str = "All files (*)",
    on_changed=None,
) -> QWidget:
    container = QWidget(parent)
    layout = QHBoxLayout(container)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(6)
    layout.addWidget(line_edit, stretch=1)

    def browse() -> None:
        start = line_edit.text().strip() or str(Path.home())
        if is_dir:
            chosen = QFileDialog.getExistingDirectory(parent, browse_title, start)
            if chosen:
                line_edit.setText(chosen)
                if on_changed is not None:
                    on_changed()
        else:
            chosen, _ = QFileDialog.getOpenFileName(
                parent,
                browse_title,
                start,
                file_filter,
            )
            if chosen:
                line_edit.setText(chosen)
                if on_changed is not None:
                    on_changed()

    browse_btn = QPushButton("Browse", container)
    browse_btn.clicked.connect(browse)
    layout.addWidget(browse_btn)
    return container


def _section_label(text: str, parent: QWidget) -> QLabel:
    label = QLabel(text, parent)
    label.setStyleSheet("font-weight: 600; margin-top: 10px;")
    return label


class SettingsDialog(QDialog):
    def __init__(self, project: TEviewerProject, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("TEviewer — Preferences")
        self.setMinimumWidth(520)
        self._project = project
        self._settings = QSettings()

        self._system_blast = QCheckBox()
        self._system_blast.setChecked(project.use_system_blast)
        self._system_viewer = QCheckBox()
        self._system_viewer.setChecked(project.use_system_sequence_viewer)
        self._color_blind_scheme = QCheckBox()
        self._color_blind_scheme.setChecked(load_color_blind_nucleotide_scheme())
        self._color_blind_scheme.setToolTip(
            "Use an orange/blue/teal nucleotide palette in the Annotation sequence view "
            "(AliView-style green/blue/red when off)."
        )
        self._batch_workers = QSpinBox()
        self._batch_workers.setRange(1, 32)
        self._batch_workers.setValue(project.batch_workers)
        self._batch_workers.setToolTip(
            "When several FASTA files are selected, run this many sequences in parallel. "
            "BLAST and TEAid thread counts are derived from CPU cores minus a small "
            "reserve for TEviewer, divided by this worker count."
        )
        self._toolbar_font = QSpinBox()
        self._toolbar_font.setRange(8, 16)
        self._toolbar_font.setValue(load_toolbar_font_size())
        self._toolbar_font.setSuffix(" pt")
        self._toolbar_font.setToolTip("Font size for function toolbar button labels.")
        self._browser_panel_font = QSpinBox()
        self._browser_panel_font.setRange(8, 16)
        self._browser_panel_font.setValue(load_browser_panel_font_size())
        self._browser_panel_font.setSuffix(" pt")
        self._browser_panel_font.setToolTip(
            "Font size for folder and file names in the Folders and Files panels."
        )
        self._preview_zoom = QSpinBox()
        self._preview_zoom.setRange(MIN_PREVIEW_ZOOM_PERCENT, MAX_PREVIEW_ZOOM_PERCENT)
        self._preview_zoom.setValue(load_preview_zoom_percent())
        self._preview_zoom.setSuffix(" %")
        self._preview_zoom.setToolTip(
            "Default scale for HTML plots, images, and text files in the Preview panel."
        )
        self._fasta_default_tab = QComboBox()
        for label, tab in (
            ("Text", DefaultFastaPreviewTab.TEXT),
            ("Annotation", DefaultFastaPreviewTab.ANNOTATION),
        ):
            self._fasta_default_tab.addItem(label, tab.value)
        fasta_tab_index = self._fasta_default_tab.findData(
            load_default_fasta_preview_tab().value
        )
        if fasta_tab_index >= 0:
            self._fasta_default_tab.setCurrentIndex(fasta_tab_index)
        self._fasta_default_tab.setToolTip(
            "Sub-tab shown when you single-click a single-sequence FASTA file."
        )
        self._preview_fit = QComboBox()
        for label, mode in (
            ("Fit window", PreviewFitMode.FIT_WINDOW),
            ("Fit width", PreviewFitMode.FIT_WIDTH),
            ("Original size", PreviewFitMode.ORIGINAL),
        ):
            self._preview_fit.addItem(label, mode.value)
        fit_index = self._preview_fit.findData(load_preview_fit_mode().value)
        if fit_index >= 0:
            self._preview_fit.setCurrentIndex(fit_index)

        for spin in (
            self._batch_workers,
            self._toolbar_font,
            self._browser_panel_font,
            self._preview_zoom,
        ):
            spin.setFixedWidth(88)
        self._preview_fit.setFixedWidth(160)
        self._fasta_default_tab.setFixedWidth(160)

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setFormAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setHorizontalSpacing(8)
        form.setVerticalSpacing(8)
        form.addRow("Use system BLAST:", self._system_blast)
        form.addRow("Use system sequence viewer:", self._system_viewer)
        form.addRow("Color-blind nucleotide:", self._color_blind_scheme)
        form.addRow("Parallel batch workers:", self._batch_workers)
        form.addRow("Toolbar button font:", self._toolbar_font)
        form.addRow("Folders / Files font:", self._browser_panel_font)
        form.addRow("Default preview zoom:", self._preview_zoom)
        form.addRow("Default preview fit:", self._preview_fit)
        form.addRow("Default FASTA preview tab:", self._fasta_default_tab)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(10)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def apply_to_project(self) -> TEviewerProject:
        self._project.use_system_blast = self._system_blast.isChecked()
        self._project.use_system_sequence_viewer = self._system_viewer.isChecked()
        self._project.batch_workers = self._batch_workers.value()
        save_toolbar_font_size(self._toolbar_font.value())
        save_browser_panel_font_size(self._browser_panel_font.value())
        save_preview_zoom_percent(self._preview_zoom.value())
        save_preview_fit_mode(PreviewFitMode(str(self._preview_fit.currentData())))
        save_default_fasta_preview_tab(
            DefaultFastaPreviewTab(str(self._fasta_default_tab.currentData()))
        )
        save_color_blind_nucleotide_scheme(self._color_blind_scheme.isChecked())
        self._save_settings()
        return self._project

    def _save_settings(self) -> None:
        s = self._settings
        s.beginGroup(SETTINGS_GROUP)
        s.setValue("preview_zoom_percent", self._preview_zoom.value())
        s.setValue("preview_fit_mode", str(self._preview_fit.currentData()))
        s.endGroup()
        save_workspace(self._project)

    @staticmethod
    def load_project_from_settings() -> TEviewerProject:
        s = QSettings()
        s.beginGroup(SETTINGS_GROUP)

        project = TEviewerProject(
            input_dir=_path_from_settings(s.value("input_dir"))
            or _path_from_settings(s.value("proof_curation_dir")),
            genome_file=_path_from_settings(s.value("genome_file")),
            genome_mode=str(s.value("genome_mode", GENOME_MODE_SINGLE)),
            pan_genomes=parse_pan_genomes_json(s.value("pan_genomes_json")),
            cdd_dirs=load_cdd_dirs_from_settings(s),
            user_databases=load_user_databases_from_settings(),
            use_system_blast=bool(s.value("use_system_blast", False)),
            use_system_sequence_viewer=bool(s.value("use_system_sequence_viewer", False)),
            batch_workers=int(s.value("batch_workers", 4)),
        )
        s.endGroup()
        return project
