"""Project path configuration."""

from __future__ import annotations

import json
import sys
from collections.abc import Iterable
from pathlib import Path
from typing import TYPE_CHECKING

from PySide6.QtCore import QRect, QSettings, QSize, Qt, QTimer, QUrl, Signal
from PySide6.QtGui import QColor, QDesktopServices, QPalette
from PySide6.QtWidgets import (
    QApplication,
    QCheckBox,
    QColorDialog,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.combo_boxes import NoWheelComboBox
from repeatviewer.widgets.spin_boxes import NoWheelSpinBox

if TYPE_CHECKING:
    from repeatviewer.plugins.registry import PluginRegistry

from repeatviewer.pan_genome import GENOME_MODE_SINGLE, parse_pan_genomes_json
from repeatviewer.project import RepeatViewerProject
from repeatviewer.views.toolbar_dialog import (
    ToolbarCustomizeWidget,
    build_toolbar_catalog,
    load_toolbar_order,
    save_toolbar_font_size,
    save_toolbar_order,
)
from repeatviewer.branding import APP_DISPLAY_NAME
from repeatviewer.app_preferences import (
    cache_size_bytes,
    crash_reports_dir,
    clear_cache,
    config_dir,
    human_size,
    load_blast_db_base_dir,
    load_cyto_experimenter,
    load_reopen_last_folder,
    load_restore_geometry,
    load_ui_thread_reserve,
    load_verbose_logging,
    log_file_path,
    reset_all_preferences,
    save_blast_db_base_dir,
    save_cyto_experimenter,
    save_reopen_last_folder,
    save_restore_geometry,
    save_ui_thread_reserve,
    save_verbose_logging,
)
from repeatviewer.conda_tools import (
    apply_conda_tools_path,
    describe_conda_tools_status,
    load_conda_tools_env,
    save_conda_tools_env,
)
from repeatviewer.ui.theme_mode import (
    DEFAULT_ACCENT_COLOR,
    THEME_LABELS,
    load_accent_color,
    load_theme_mode,
    save_accent_color,
    save_corner_shape,
    save_theme_mode,
)
from repeatviewer.widgets.browser_panel_settings import (
    load_browser_panel_font_size,
    save_browser_panel_font_size,
)
from repeatviewer.widgets.preview_zoom import (
    MAX_PREVIEW_ZOOM_PERCENT,
    MIN_PREVIEW_ZOOM_PERCENT,
    PreviewFitMode,
    DefaultFastaPreviewTab,
    load_color_blind_nucleotide_scheme,
    load_default_fasta_preview_tab,
    load_preview_fit_mode,
    load_preview_zoom_percent,
    load_show_sequence_summary,
    save_color_blind_nucleotide_scheme,
    save_default_fasta_preview_tab,
    save_preview_fit_mode,
    save_preview_zoom_percent,
    save_show_sequence_summary,
)

from repeatviewer.database_settings import (
    SETTINGS_GROUP,
    load_blast_databases_from_settings,
    load_cdd_dirs_from_settings,
    load_user_databases_from_settings,
    save_blast_databases_to_settings,
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


def _bool_from_settings(value, default: bool = False) -> bool:
    """Parse QSettings bools; string \"false\" must not become True via bool()."""
    if value is None:
        return default
    if isinstance(value, bool):
        return value
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"1", "true", "yes", "on"}
    return bool(value)


def _path_list_from_settings(value) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value] if value.strip() else []
    return [str(v).strip() for v in value if v and str(v).strip()]


def save_workspace(
    project: RepeatViewerProject,
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
    from repeatviewer.pan_genome import pan_genomes_to_json

    s.setValue("pan_genomes_json", pan_genomes_to_json(project.pan_genomes))
    save_database_settings(
        cdd_dirs=project.cdd_dirs,
        user_databases=project.user_databases,
        s=s,
    )
    save_blast_databases_to_settings(project.blast_databases, s=s)
    s.setValue("use_system_blast", project.use_system_blast)
    s.setValue("use_system_mafft", project.use_system_mafft)
    s.setValue("use_system_muscle", project.use_system_muscle)
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
    container.setSizePolicy(
        QSizePolicy.Policy.Expanding,
        QSizePolicy.Policy.Fixed,
    )
    layout = QHBoxLayout(container)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(6)
    line_edit.setSizePolicy(
        QSizePolicy.Policy.Expanding,
        QSizePolicy.Policy.Fixed,
    )
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
    browse_btn.setSizePolicy(
        QSizePolicy.Policy.Fixed,
        QSizePolicy.Policy.Fixed,
    )
    browse_btn.clicked.connect(browse)
    layout.addWidget(browse_btn)
    return container


def _section_label(text: str, parent: QWidget) -> QLabel:
    label = QLabel(text, parent)
    font = label.font()
    font.setBold(True)
    label.setFont(font)
    return label


# Uniform column sizing so every field lines up across all groups.
_FIELD_WIDTH = 170
_PATH_FIELD_MIN_WIDTH = 320
_LABEL_WIDTH = 200
_CONDA_STATUS_MAX_LINES = 4


class _CondaStatusLabel(QLabel):
    """Word-wrapped status that never expands the dialog to the unwrapped text width."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._wrap_width = _PATH_FIELD_MIN_WIDTH
        self.setWordWrap(True)
        self.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self.setAlignment(Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Preferred)

    def set_wrap_width(self, width: int) -> None:
        width = max(80, width)
        if width == self._wrap_width:
            return
        self._wrap_width = width
        self.updateGeometry()

    def hasHeightForWidth(self) -> bool:  # noqa: N802
        return True

    def heightForWidth(self, width: int) -> int:  # noqa: N802
        w = width if width > 0 else self._wrap_width
        fm = self.fontMetrics()
        bounds = fm.boundingRect(
            QRect(0, 0, w, 0),
            int(Qt.TextFlag.TextWordWrap),
            self.text() or " ",
        )
        return max(bounds.height() + 2, fm.lineSpacing())

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        # Critical: do NOT use the unwrapped string width (that forces a horizontal scrollbar).
        return QSize(80, self.fontMetrics().lineSpacing())

    def sizeHint(self) -> QSize:  # noqa: N802
        return QSize(self._wrap_width, self.heightForWidth(self._wrap_width))


class SettingsDialog(QDialog):
    colorBlindToggled = Signal(bool)  # emitted immediately when the checkbox is clicked

    def __init__(
        self,
        project: RepeatViewerProject,
        parent: QWidget | None = None,
        *,
        plugin_registry: PluginRegistry | None = None,
        toolbar_order: list[str] | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Preferences")
        self.setMinimumWidth(720)
        self.setMinimumHeight(420)
        self._project = project
        self._settings = QSettings()
        self._plugin_registry = plugin_registry

        self._system_viewer = QCheckBox()
        self._system_viewer.setChecked(project.use_system_sequence_viewer)
        self._color_blind_scheme = QCheckBox()
        self._color_blind_scheme.setChecked(load_color_blind_nucleotide_scheme())
        self._color_blind_scheme.setToolTip(
            "Use an orange/blue/teal nucleotide palette in the Annotation sequence view "
            "(AliView-style green/blue/red when off)."
        )
        self._color_blind_scheme.toggled.connect(self._on_color_blind_toggled)
        self._batch_workers = NoWheelSpinBox()
        self._batch_workers.setRange(1, 32)
        self._batch_workers.setValue(project.batch_workers)
        self._batch_workers.setToolTip(
            "When several FASTA files are selected, run this many sequences in parallel. "
            "BLAST and TEAid thread counts are derived from CPU cores minus a small "
            f"reserve for {APP_DISPLAY_NAME}, divided by this worker count."
        )
        self._browser_panel_font = NoWheelSpinBox()
        self._browser_panel_font.setRange(8, 16)
        self._browser_panel_font.setValue(load_browser_panel_font_size())
        self._browser_panel_font.setSuffix(" pt")
        self._browser_panel_font.setToolTip(
            "Font size for folder and file names in the Browser panel."
        )
        self._preview_zoom = NoWheelSpinBox()
        self._preview_zoom.setRange(MIN_PREVIEW_ZOOM_PERCENT, MAX_PREVIEW_ZOOM_PERCENT)
        self._preview_zoom.setValue(load_preview_zoom_percent())
        self._preview_zoom.setSuffix(" %")
        self._preview_zoom.setToolTip(
            "Default scale for HTML plots, images, and text files in the Preview panel."
        )
        self._fasta_default_tab = NoWheelComboBox()
        for label, tab in (
            ("Text", DefaultFastaPreviewTab.TEXT),
            ("Annotation", DefaultFastaPreviewTab.ANNOTATION),
            ("Dotplot (Self)", DefaultFastaPreviewTab.DOTPLOT),
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
        self._show_sequence_summary = QCheckBox("Show length and GC% in the preview")
        self._show_sequence_summary.setChecked(load_show_sequence_summary())
        self._show_sequence_summary.setToolTip(
            "GC% requires scanning every base. On very large single-sequence "
            "FASTA files (whole scaffolds or genomes) that costs a few seconds "
            "each time the file is opened. Turn this off to skip it."
        )
        self._preview_fit = NoWheelComboBox()
        for label, mode in (
            ("Fit window", PreviewFitMode.FIT_WINDOW),
            ("Fit width", PreviewFitMode.FIT_WIDTH),
            ("Original size", PreviewFitMode.ORIGINAL),
        ):
            self._preview_fit.addItem(label, mode.value)
        fit_index = self._preview_fit.findData(load_preview_fit_mode().value)
        if fit_index >= 0:
            self._preview_fit.setCurrentIndex(fit_index)

        self._theme_mode = NoWheelComboBox()
        for label, mode in THEME_LABELS:
            self._theme_mode.addItem(label, mode)
        theme_index = self._theme_mode.findData(load_theme_mode())
        if theme_index >= 0:
            self._theme_mode.setCurrentIndex(theme_index)
        self._theme_mode.setToolTip(
            "Application colour theme. Auto follows your operating system's "
            "light/dark setting."
        )

        # Accent colour (drives selection / highlight colour across the app).
        self._accent_color = load_accent_color()
        self._accent_swatch = QPushButton()
        self._accent_swatch.setFixedSize(_FIELD_WIDTH - 62, 28)  # swatch + Reset == field width
        self._accent_swatch.setCursor(Qt.CursorShape.PointingHandCursor)
        self._accent_swatch.setToolTip("Click to choose the highlight / selection colour.")
        self._accent_swatch.clicked.connect(self._pick_accent_color)
        self._accent_reset = QPushButton("Reset")
        self._accent_reset.setFixedWidth(56)
        self._accent_reset.clicked.connect(self._reset_accent_color)
        self._update_accent_swatch()

        # Performance
        self._ui_thread_reserve = NoWheelSpinBox()
        self._ui_thread_reserve.setRange(0, 8)
        self._ui_thread_reserve.setValue(load_ui_thread_reserve())
        self._ui_thread_reserve.setSpecialValueText("Auto")
        self._ui_thread_reserve.setToolTip(
            "CPU threads kept free for the interface while BLAST/TEAid jobs run. "
            "Auto reserves a small default. Raise it if the UI feels sluggish during "
            "heavy analyses; lower it to give jobs more cores."
        )

        # Startup
        self._reopen_last = QCheckBox()
        self._reopen_last.setChecked(load_reopen_last_folder())
        self._reopen_last.setToolTip("Reopen the last folder and genome when the app starts.")
        self._restore_geometry = QCheckBox()
        self._restore_geometry.setChecked(load_restore_geometry())
        self._restore_geometry.setToolTip("Restore the window size and position from the last session.")

        # Cytogenetic images
        self._cyto_experimenter = QLineEdit(load_cyto_experimenter())
        self._cyto_experimenter.setPlaceholderText("Your name")
        self._cyto_experimenter.setToolTip(
            "Filled in automatically when you attach a cytogenetic image to a sequence."
        )

        # BLAST databases
        self._blast_db_base_dir = QLineEdit(str(load_blast_db_base_dir()))
        self._blast_db_base_dir.setMinimumWidth(_PATH_FIELD_MIN_WIDTH)
        self._open_blast_db_btn = QPushButton("Open BLAST databases folder")
        self._open_blast_db_btn.clicked.connect(self._open_blast_db_folder)

        # External tools via conda (RepeatMasker, FastGA, bedtools, …)
        self._conda_tools_env = QLineEdit(load_conda_tools_env())
        self._conda_tools_env.setMinimumWidth(_PATH_FIELD_MIN_WIDTH)
        if sys.platform == "win32":
            self._conda_tools_env.setPlaceholderText(
                "WSL2 env name or path, e.g. tetrimmer"
            )
            self._conda_tools_env.setToolTip(
                "Conda/mamba environment inside WSL2 that provides external CLI tools "
                "(RepeatMasker, FastGA, bedtools, …).\n"
                "Enter the env name (e.g. tetrimmer), a Linux path "
                "(~/miniconda3/envs/tetrimmer), or Browse under \\\\wsl$\\…\n"
                "Leave empty to use only the Windows PATH. Bundled BLAST/MAFFT/MUSCLE "
                "stay on Windows; only Linux conda tools go through WSL."
            )
            conda_browse_title = "Choose WSL conda env (\\\\wsl$\\…) or Windows conda folder"
        else:
            self._conda_tools_env.setPlaceholderText("env name or path, e.g. tetrimmer")
            self._conda_tools_env.setToolTip(
                "Conda/mamba environment that provides external CLI tools "
                "(RepeatMasker, FastGA, bedtools, …).\n"
                "Enter the env name, or browse to the env folder "
                "(…/envs/your_env). Leave empty to use only the system PATH."
            )
            conda_browse_title = "Choose conda environment folder"
        # Do not resolve the env on every keystroke — that freezes typing.
        # Only update a cheap idle hint while typing; resolve on Test.
        self._conda_tools_env.textChanged.connect(self._show_conda_tools_idle_status)
        # Full status text; wraps within the path field; at most 4 lines (scroll if longer).
        self._conda_tools_status = _CondaStatusLabel()
        status_pal = self._conda_tools_status.palette()
        muted = status_pal.color(QPalette.ColorRole.PlaceholderText)
        if not muted.isValid() or muted.alpha() == 0:
            muted = status_pal.color(QPalette.ColorRole.WindowText)
            muted.setAlpha(180)
        status_pal.setColor(QPalette.ColorRole.WindowText, muted)
        self._conda_tools_status.setPalette(status_pal)

        self._conda_tools_status_scroll = QScrollArea()
        self._conda_tools_status_scroll.setWidget(self._conda_tools_status)
        self._conda_tools_status_scroll.setWidgetResizable(False)
        self._conda_tools_status_scroll.setFrameShape(QFrame.Shape.NoFrame)
        self._conda_tools_status_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self._conda_tools_status_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )
        self._conda_tools_status_scroll.setSizePolicy(
            QSizePolicy.Policy.Fixed,
            QSizePolicy.Policy.Fixed,
        )
        _line = max(self._conda_tools_status.fontMetrics().lineSpacing(), 1)
        self._conda_tools_status_scroll.setFixedSize(_PATH_FIELD_MIN_WIDTH, _line * 2 + 4)

        self._conda_test_btn = QPushButton("Test")
        self._conda_test_btn.setToolTip("Check whether the env resolves and which tools are visible.")
        self._conda_test_btn.clicked.connect(self._on_test_conda_tools)

        # Diagnostics
        self._verbose_logging = QCheckBox()
        self._verbose_logging.setChecked(load_verbose_logging())
        self._verbose_logging.setToolTip("Write detailed DEBUG-level messages to the log file.")
        self._clear_cache_btn = QPushButton()
        self._clear_cache_btn.clicked.connect(self._on_clear_cache)
        self._refresh_cache_button()
        self._open_log_btn = QPushButton("Open log file")
        self._open_log_btn.clicked.connect(self._open_log_file)
        self._open_config_btn = QPushButton("Open config folder")
        self._open_config_btn.clicked.connect(self._open_config_folder)
        self._copy_diag_btn = QPushButton("Copy diagnostics")
        self._copy_diag_btn.setToolTip(
            "Copy version, OS and folder locations to the clipboard for bug reports."
        )
        self._copy_diag_btn.clicked.connect(self._copy_diagnostics)

        # Every combo and spin box gets the same width so the columns are symmetric.
        for widget in (
            self._batch_workers,
            self._browser_panel_font,
            self._preview_zoom,
            self._ui_thread_reserve,
            self._preview_fit,
            self._fasta_default_tab,
            self._theme_mode,
        ):
            widget.setFixedWidth(_FIELD_WIDTH)

        def _group(title: str, rows: list) -> QGroupBox:
            box = QGroupBox(title)
            f = QFormLayout(box)
            f.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            f.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
            f.setHorizontalSpacing(10)
            f.setVerticalSpacing(8)
            for label, widget in rows:
                if label is None:
                    f.addRow(widget)
                else:
                    lbl = QLabel(label)
                    lbl.setFixedWidth(_LABEL_WIDTH)  # constant label column => fields align
                    lbl.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                    f.addRow(lbl, widget)
            return box

        # Accent swatch + Reset together fill exactly one field width, so Reset's
        # right edge lines up with the combos above/below.
        accent_field = QWidget()
        _accent_h = QHBoxLayout(accent_field)
        _accent_h.setContentsMargins(0, 0, 0, 0)
        _accent_h.setSpacing(6)
        _accent_h.addWidget(self._accent_swatch)
        _accent_h.addWidget(self._accent_reset)
        accent_field.setFixedWidth(_FIELD_WIDTH)

        appearance = _group("Appearance", [
            ("Theme:", self._theme_mode),
            ("Accent colour:", accent_field),
        ])
        seq_tools = _group("Sequence tools", [
            ("Use system sequence viewer:", self._system_viewer),
            ("Color-blind nucleotide:", self._color_blind_scheme),
        ])
        files_preview = _group("Files && preview", [
            ("Folders / Files font:", self._browser_panel_font),
            ("Default preview zoom:", self._preview_zoom),
            ("Default preview fit:", self._preview_fit),
            ("Default FASTA preview tab:", self._fasta_default_tab),
        ])
        performance = _group("Performance", [
            ("Threads for tools:", self._batch_workers),
            ("Threads for the UI:", self._ui_thread_reserve),
        ])
        startup = _group("Startup", [
            ("Reopen last folder:", self._reopen_last),
            ("Restore window size:", self._restore_geometry),
        ])
        cytogenetics = _group("Cytogenetic images", [
            ("Experimenter:", self._cyto_experimenter),
        ])
        blast_databases = _group("BLAST databases", [
            ("Base folder:", _path_browse_row(
                self._blast_db_base_dir,
                parent=self,
                is_dir=True,
                browse_title="Choose BLAST databases folder",
            )),
            ("", self._open_blast_db_btn),
        ])
        conda_tools_field = _path_browse_row(
            self._conda_tools_env,
            parent=self,
            is_dir=True,
            browse_title=conda_browse_title,
        )
        # Status sits under the path text box (form field column), wraps to that
        # field's right edge, max 4 lines.
        conda_tools = _group(
            "External tools (conda / WSL2)" if sys.platform == "win32" else "External tools (conda)",
            [
                ("Conda env name or path:", conda_tools_field),
                ("", self._conda_tools_status_scroll),
                ("", self._conda_test_btn),
            ],
        )
        diagnostics = _group("Diagnostics", [
            ("Verbose logging:", self._verbose_logging),
            ("", self._open_log_btn),
            ("", self._copy_diag_btn),
            ("", self._clear_cache_btn),
            ("", self._open_config_btn),
        ])

        content = QWidget()
        content_v = QVBoxLayout(content)
        content_v.setContentsMargins(10, 10, 10, 10)
        content_v.setSpacing(10)
        for box in (
            appearance,
            seq_tools,
            files_preview,
            performance,
            startup,
            cytogenetics,
            blast_databases,
            conda_tools,
            diagnostics,
        ):
            content_v.addWidget(box)
        content_v.addStretch()

        scroll = QScrollArea()
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(content)

        self._reset_btn = QPushButton("Restore defaults")
        self._reset_btn.clicked.connect(self._on_reset_defaults)
        config_lbl = QLabel(f"Config: {config_dir()}")
        config_lbl.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        config_lbl.setStyleSheet("color: palette(mid);")
        footer = QWidget()
        footer_h = QHBoxLayout(footer)
        footer_h.setContentsMargins(10, 0, 10, 4)
        footer_h.addWidget(self._reset_btn)
        footer_h.addStretch()
        footer_h.addWidget(config_lbl)

        general_tab = QWidget()
        general_layout = QVBoxLayout(general_tab)
        general_layout.setContentsMargins(0, 0, 0, 0)
        general_layout.setSpacing(8)
        general_layout.addWidget(scroll, stretch=1)
        general_layout.addWidget(footer)

        order = toolbar_order if toolbar_order is not None else load_toolbar_order(
            plugin_registry
        )
        catalog = build_toolbar_catalog(order, plugin_registry)
        self._toolbar_panel = ToolbarCustomizeWidget(
            catalog,
            registry=plugin_registry,
        )
        toolbar_tab = QWidget()
        toolbar_layout = QVBoxLayout(toolbar_tab)
        toolbar_layout.setContentsMargins(12, 12, 12, 16)
        toolbar_layout.addWidget(self._toolbar_panel)

        tabs = QTabWidget()
        tabs.addTab(general_tab, "General")
        tabs.addTab(toolbar_tab, "Toolbar")

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        ok_btn = buttons.button(QDialogButtonBox.StandardButton.Ok)
        if ok_btn is not None:
            ok_btn.setText("Save")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(10)
        layout.addWidget(tabs)
        layout.addWidget(buttons)

    def _on_color_blind_toggled(self, checked: bool) -> None:
        """Persist and broadcast immediately so open Annotation views recolor live,
        without waiting for OK - this is the kind of setting you want to preview."""
        save_color_blind_nucleotide_scheme(checked)
        self.colorBlindToggled.emit(checked)

    def accept(self) -> None:
        if not self._toolbar_panel.has_visible_tools():
            QMessageBox.warning(
                self,
                "Setting",
                "Show at least one toolbar button.",
            )
            return
        super().accept()

    def apply_toolbar_settings(self) -> tuple[list[str], int]:
        order = self._toolbar_panel.result_order()
        font_size = self._toolbar_panel.result_font_size()
        save_toolbar_order(order, self._plugin_registry)
        save_toolbar_font_size(font_size)
        return order, font_size

    def apply_to_project(self) -> RepeatViewerProject:
        self._project.use_system_sequence_viewer = self._system_viewer.isChecked()
        self._project.batch_workers = self._batch_workers.value()
        save_browser_panel_font_size(self._browser_panel_font.value())
        save_preview_zoom_percent(self._preview_zoom.value())
        save_preview_fit_mode(PreviewFitMode(str(self._preview_fit.currentData())))
        save_show_sequence_summary(self._show_sequence_summary.isChecked())
        save_default_fasta_preview_tab(
            DefaultFastaPreviewTab(str(self._fasta_default_tab.currentData()))
        )
        save_color_blind_nucleotide_scheme(self._color_blind_scheme.isChecked())
        save_theme_mode(str(self._theme_mode.currentData()))
        save_accent_color(self._accent_color)
        save_corner_shape("rounded")
        save_ui_thread_reserve(self._ui_thread_reserve.value())
        save_reopen_last_folder(self._reopen_last.isChecked())
        save_restore_geometry(self._restore_geometry.isChecked())
        save_cyto_experimenter(self._cyto_experimenter.text())
        save_verbose_logging(self._verbose_logging.isChecked())
        base_dir_text = self._blast_db_base_dir.text().strip()
        if base_dir_text:
            save_blast_db_base_dir(Path(base_dir_text))
        conda_spec = self._conda_tools_env.text().strip()
        save_conda_tools_env(conda_spec)
        apply_conda_tools_path(conda_spec)
        self._save_settings()
        return self._project

    def selected_theme_mode(self) -> str:
        return str(self._theme_mode.currentData())

    def selected_accent_color(self) -> str:
        return self._accent_color

    def appearance_signature(self) -> tuple[str, str]:
        """Values that require re-applying the theme when any of them changes."""
        return (
            self.selected_theme_mode(),
            self.selected_accent_color(),
        )

    # --- Appearance handlers ---------------------------------------------

    def _update_accent_swatch(self) -> None:
        self._accent_swatch.setStyleSheet(
            f"background-color: {self._accent_color}; border: 1px solid palette(mid); "
            "border-radius: 4px;"
        )

    def _pick_accent_color(self) -> None:
        color = QColorDialog.getColor(
            QColor(self._accent_color), self, "Choose accent colour"
        )
        if color.isValid():
            self._accent_color = color.name()
            self._update_accent_swatch()

    def _reset_accent_color(self) -> None:
        self._accent_color = DEFAULT_ACCENT_COLOR
        self._update_accent_swatch()

    # --- Diagnostics handlers --------------------------------------------

    def _refresh_cache_button(self) -> None:
        self._clear_cache_btn.setText(f"Clear cache ({human_size(cache_size_bytes())})")

    def _on_clear_cache(self) -> None:
        freed = clear_cache()
        self._refresh_cache_button()
        QMessageBox.information(
            self, "Cache cleared", f"Freed {human_size(freed)} of cached files."
        )

    def _open_log_file(self) -> None:
        path = log_file_path()
        if path.exists():
            target = path
        else:
            target = config_dir()
            target.mkdir(parents=True, exist_ok=True)
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(target)))

    def _open_config_folder(self) -> None:
        folder = config_dir()
        folder.mkdir(parents=True, exist_ok=True)
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))

    def _open_blast_db_folder(self) -> None:
        text = self._blast_db_base_dir.text().strip()
        folder = Path(text).expanduser() if text else load_blast_db_base_dir()
        folder.mkdir(parents=True, exist_ok=True)
        QDesktopServices.openUrl(QUrl.fromLocalFile(str(folder)))

    def _refresh_conda_tools_status(self) -> None:
        status = describe_conda_tools_status(self._conda_tools_env.text())
        self._conda_tools_status.setText(status)
        self._conda_tools_status.setToolTip(status)
        QTimer.singleShot(0, self._fit_conda_tools_status_height)

    def _show_conda_tools_idle_status(self) -> None:
        """Lightweight hint — no filesystem / conda lookup while typing."""
        text = self._conda_tools_env.text().strip()
        if text:
            status = f"Click Test to check “{text}”."
        else:
            status = "Enter a conda env name or path, then click Test."
        self._conda_tools_status.setText(status)
        self._conda_tools_status.setToolTip(status)
        QTimer.singleShot(0, self._fit_conda_tools_status_height)

    def _on_test_conda_tools(self) -> None:
        # Always show an immediate change so a repeat Test is visible even when
        # the resolved path is unchanged.
        self._conda_tools_status.setText("Testing…")
        self._conda_tools_status.setToolTip("Testing…")
        self._fit_conda_tools_status_height()
        app = QApplication.instance()
        if app is not None:
            app.processEvents()
        self._refresh_conda_tools_status()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._show_conda_tools_idle_status()
        QTimer.singleShot(0, self._fit_conda_tools_status_height)
        QTimer.singleShot(50, self._fit_conda_tools_status_height)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        if hasattr(self, "_conda_tools_status_scroll"):
            self._fit_conda_tools_status_height()

    def _conda_status_wrap_width(self) -> int:
        """Width of the path QLineEdit, capped so the dialog never grows horizontally."""
        path_w = self._conda_tools_env.width()
        if path_w < 80:
            parent = self._conda_tools_env.parentWidget()
            if parent is not None and parent.width() >= 120:
                browse = None
                for child in parent.findChildren(QPushButton):
                    if child.text() == "Browse":
                        browse = child
                        break
                if browse is not None:
                    path_w = parent.width() - browse.width() - 6
                else:
                    path_w = parent.width() - 80
            else:
                path_w = _PATH_FIELD_MIN_WIDTH
        # Cap to dialog client width minus label column / margins.
        dialog_cap = max(160, self.width() - _LABEL_WIDTH - 64)
        return max(120, min(path_w, dialog_cap))

    def _fit_conda_tools_status_height(self) -> None:
        """Wrap inside the path field width; show up to 4 lines (scroll after that)."""
        if not hasattr(self, "_conda_tools_status_scroll"):
            return
        label = self._conda_tools_status
        scroll = self._conda_tools_status_scroll
        wrap_w = self._conda_status_wrap_width()
        label.set_wrap_width(wrap_w)
        label.setFixedWidth(wrap_w)

        line_h = max(label.fontMetrics().lineSpacing(), 1)
        max_h = line_h * _CONDA_STATUS_MAX_LINES + 6
        needed = label.heightForWidth(wrap_w)
        height = min(max(needed, line_h), max_h)
        label.setFixedHeight(needed)  # full text height inside scroll area
        scroll.setFixedSize(wrap_w, height)
        scroll.updateGeometry()

    def _copy_diagnostics(self) -> None:
        import platform

        from PySide6 import __version__ as pyside_version
        from PySide6.QtCore import qVersion
        from PySide6.QtWidgets import QApplication

        from repeatviewer import __version__
        from repeatviewer.app_preferences import cache_dir, log_file_path

        try:
            from repeatviewer.backends.gui_functions import (
                flexidot_executable,
                flexidot_supports_html,
            )

            flex_exe = flexidot_executable() or "not found"
            flex_html = "yes" if flexidot_supports_html() else "no"
        except Exception as exc:  # noqa: BLE001
            flex_exe = f"error: {exc}"
            flex_html = "?"

        info = "\n".join([
            f"{APP_DISPLAY_NAME} {__version__}",
            f"Python {platform.python_version()} on {platform.system()} {platform.release()}",
            f"PySide6 {pyside_version} (Qt {qVersion()})",
            f"Config:   {config_dir()}",
            f"Cache:    {cache_dir()}",
            f"Log:      {log_file_path()}",
            f"Crashes:  {crash_reports_dir()}",
            f"Flexidot: {flex_exe} (html: {flex_html})",
        ])
        QApplication.clipboard().setText(info)
        QMessageBox.information(
            self, "Diagnostics copied", "System information was copied to the clipboard."
        )

    def _on_reset_defaults(self) -> None:
        confirm = QMessageBox.question(
            self,
            "Restore defaults",
            "Reset all preferences to their defaults?\n\n"
            "Your folders, genomes, files, workflows, and command wrappers "
            "are not affected. "
            f"{APP_DISPLAY_NAME} must be restarted for this to take full effect.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if confirm != QMessageBox.StandardButton.Yes:
            return
        reset_all_preferences()
        QMessageBox.information(
            self,
            "Defaults restored",
            f"Preferences were reset. Please restart {APP_DISPLAY_NAME}.",
        )
        self.reject()

    def _save_settings(self) -> None:
        s = self._settings
        s.beginGroup(SETTINGS_GROUP)
        s.setValue("preview_zoom_percent", self._preview_zoom.value())
        s.setValue("preview_fit_mode", str(self._preview_fit.currentData()))
        s.endGroup()
        save_workspace(self._project)

    @staticmethod
    def load_project_from_settings() -> RepeatViewerProject:
        s = QSettings()
        s.beginGroup(SETTINGS_GROUP)

        project = RepeatViewerProject(
            input_dir=_path_from_settings(s.value("input_dir"))
            or _path_from_settings(s.value("proof_curation_dir")),
            genome_file=_path_from_settings(s.value("genome_file")),
            genome_mode=str(s.value("genome_mode", GENOME_MODE_SINGLE)),
            pan_genomes=parse_pan_genomes_json(s.value("pan_genomes_json")),
            cdd_dirs=load_cdd_dirs_from_settings(s),
            user_databases=load_user_databases_from_settings(),
            blast_databases=load_blast_databases_from_settings(),
            use_system_blast=_bool_from_settings(s.value("use_system_blast", False)),
            use_system_mafft=_bool_from_settings(s.value("use_system_mafft", False)),
            use_system_muscle=_bool_from_settings(s.value("use_system_muscle", False)),
            use_system_sequence_viewer=_bool_from_settings(
                s.value("use_system_sequence_viewer", False)
            ),
            batch_workers=int(s.value("batch_workers", 4)),
        )
        s.endGroup()
        return project
