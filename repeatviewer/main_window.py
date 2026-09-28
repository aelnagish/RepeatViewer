"""RepeatViewer main window with folder navigator and file workspace."""

from __future__ import annotations

import logging
import os
import shutil
import platform
import tempfile
import time
import re
import uuid
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any, Literal

from PySide6.QtCore import (
    QEvent,
    QFileSystemWatcher,
    QItemSelectionModel,
    QMimeData,
    QObject,
    QPoint,
    QRect,
    QSettings,
    QSize,
    Qt,
    QThread,
    QTimer,
    Signal,
    Slot,
)
from PySide6.QtGui import (
    QAction,
    QBrush,
    QColor,
    QDrag,
    QFont,
    QFontMetrics,
    QKeySequence,
    QMouseEvent,
    QPainter,
    QPalette,
    QShortcut,
    QWheelEvent,
)
from PySide6.QtWidgets import (
    QApplication,
    QAbstractItemView,
    QCheckBox,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QDockWidget,
    QFrame,
    QGraphicsDropShadowEffect,
    QHeaderView,
    QFileDialog,
    QHBoxLayout,
    QInputDialog,
    QLineEdit,
    QMenu,
    QSizePolicy,
    QListWidget,
    QListWidgetItem,
    QMainWindow,
    QMessageBox,
    QPushButton,
    QStyle,
    QStyleOptionViewItem,
    QToolBar,
    QToolButton,
    QStatusBar,
    QTabWidget,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
    QLabel,
    QColorDialog,
)

_BROWSER_DOCK_WIDTH = 480
_BROWSER_FOLDER_TREE_WIDTH = 120
_BROWSER_FILES_PREVIEW_SPLIT = (170, 210)


class _ToolbarPopup(QWidget):
    """Toolbar dropdown that does not steal the mouse grab.

    ``QMenu.popup`` grabs the mouse, so the first click on another tool only
    dismisses the open menu. This lightweight popup stays below the toolbar and
    lets FastGA → MSA clean (and other) switches happen in one click.

    Visually approximates a macOS menu: rounded panel, soft hairline border,
    inset separators, and rounded blue item highlight.
    """

    closed = Signal()
    # Side/bottom room for the drop shadow. Top stays tight so the menu does
    # not sit on top of the toolbar button that opened it.
    _SHADOW = 12
    _TOP_PAD = 2

    def __init__(self, actions: list[QAction], parent: QWidget | None = None) -> None:
        # parent=None → true top-level window so move/setGeometry use global
        # screen coordinates (matching mapToGlobal). A MainWindow parent can
        # make macOS treat positions as parent-relative and shift the menu up.
        del parent
        super().__init__(
            None,
            Qt.WindowType.Tool
            | Qt.WindowType.FramelessWindowHint
            | Qt.WindowType.NoDropShadowWindowHint,
        )
        self.setAttribute(Qt.WidgetAttribute.WA_DeleteOnClose, True)
        self.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        self.setObjectName("ToolbarPopup")
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(self._SHADOW, self._TOP_PAD, self._SHADOW, self._SHADOW)
        outer.setSpacing(0)

        panel = QWidget(self)
        self._panel = panel
        panel.setObjectName("ToolbarPopupPanel")
        panel.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        layout = QVBoxLayout(panel)
        layout.setContentsMargins(5, 5, 5, 5)
        layout.setSpacing(1)

        panel.setStyleSheet(
            "QWidget#ToolbarPopupPanel {"
            "  background: #f6f6f6;"
            "  border: 1px solid #c8c8c8;"
            "  border-radius: 8px;"
            "}"
            "QWidget#ToolbarPopupPanel QPushButton {"
            "  border: 0;"
            "  border-radius: 5px;"
            "  padding: 4px 18px 4px 12px;"
            "  text-align: left;"
            "  background: transparent;"
            "  color: #1d1d1f;"
            "  font-size: 13px;"
            "}"
            "QWidget#ToolbarPopupPanel QPushButton:hover {"
            "  background: #0a84ff;"
            "  color: #ffffff;"
            "}"
            "QWidget#ToolbarPopupPanel QPushButton:disabled {"
            "  color: #aeaeb2;"
            "  background: transparent;"
            "}"
        )

        shadow = QGraphicsDropShadowEffect(panel)
        shadow.setBlurRadius(28)
        shadow.setOffset(0, 6)
        shadow.setColor(QColor(0, 0, 0, 48))
        panel.setGraphicsEffect(shadow)

        for action in actions:
            if action.isSeparator():
                row = QWidget(panel)
                row_layout = QHBoxLayout(row)
                row_layout.setContentsMargins(8, 1, 8, 1)
                row_layout.setSpacing(0)
                separator = QFrame(row)
                separator.setFixedHeight(1)
                separator.setStyleSheet("background: #d8d8d8; border: 0;")
                row_layout.addWidget(separator)
                layout.addWidget(row)
                continue
            button = QPushButton(action.text().replace("&", ""), panel)
            button.setCursor(Qt.CursorShape.ArrowCursor)
            button.setEnabled(action.isEnabled())
            button.setToolTip(action.toolTip())
            button.clicked.connect(
                lambda _checked=False, selected=action: self._run_action(selected)
            )
            layout.addWidget(button)

        outer.addWidget(panel)

    def _run_action(self, action: QAction) -> None:
        self.close()
        action.trigger()

    def show_below(self, widget: QWidget) -> None:
        """Place the visible menu panel fully below ``widget`` (global coords)."""
        self.adjustSize()
        gap = 3
        width = max(self.sizeHint().width(), 1)
        height = max(self.sizeHint().height(), 1)

        btn_top_left = widget.mapToGlobal(QPoint(0, 0))
        btn_bottom = widget.mapToGlobal(QPoint(0, widget.height())).y()

        # Window top-left in screen coordinates (top-level Tool window).
        x = btn_top_left.x() - self._SHADOW
        y = btn_bottom + gap - self._TOP_PAD

        screen = widget.screen()
        if screen is not None:
            available = screen.availableGeometry()
            x = max(available.left(), min(x, available.right() - width + 1))
            if y + height > available.bottom():
                y = max(available.top(), btn_top_left.y() - gap - height)
            self.setGeometry(QRect(x, y, width, height))
        else:
            self.setGeometry(QRect(x, y, width, height))

        self.show()
        self.raise_()
        # Snap again after show — macOS can adjust frameless tool windows.
        self._snap_panel_below(widget, gap=gap)
        self.setFocus(Qt.FocusReason.PopupFocusReason)

    def _snap_panel_below(self, widget: QWidget, *, gap: int = 3) -> None:
        """If the opaque panel still overlaps the button, nudge the window down."""
        panel = getattr(self, "_panel", None)
        if panel is None or not widget.isVisible():
            return
        btn_bottom = widget.mapToGlobal(QPoint(0, widget.height())).y()
        panel_top = panel.mapToGlobal(QPoint(0, 0)).y()
        if panel_top < btn_bottom + gap:
            self.move(self.x(), self.y() + (btn_bottom + gap - panel_top))

    def focusOutEvent(self, event) -> None:  # noqa: N802
        super().focusOutEvent(event)
        # Defer so a click on another toolbar tool can open its menu in the
        # same gesture after this popup closes.
        QTimer.singleShot(0, self._close_if_unfocused)

    def _close_if_unfocused(self) -> None:
        if self.isVisible() and not self.isActiveWindow():
            self.close()

    def hideEvent(self, event) -> None:  # noqa: N802
        super().hideEvent(event)
        self.closed.emit()


from repeatviewer import __version__
from repeatviewer.branding import APP_DISPLAY_NAME
from repeatviewer.backends.genome_prepare import ensure_genome_artifacts
from repeatviewer.backends.blast_genome_plot import (
    _map_chromosome_filter_to_genome,
    blast_genome_map_html_path_for_blast_table,
    divergence_plot_html_path_for_te_inputs,
    genome_map_html_path_for_te_inputs,
    GENOME_MAP_HTML_SUFFIX,
    blast_txt_path_for_query,
    chromosome_names_from_genome,
    genome_chrom_sizes_from_fasta,
    find_blast_table_for_query,
    invalidate_genome_chromosome_cache,
    is_blast_table_file,
    is_divergence_plot_html_file,
    is_genome_map_html_file,
    is_repeatmasker_out_file,
    is_te_hit_table_file,
    is_repeatviewer_hits_cache_file,
    te_hit_table_stem,
    write_blast_divergence_plot_html,
    GenomeMapNothingToPlotError,
    write_blast_genome_map_html,
    TE_INPUT_KIND_CONSENSUS,
    TE_INPUT_KIND_MASKING,
)
from repeatviewer.backends.te_abundance_plot import (
    te_abundance_plot_html_path_for_sources,
    write_te_abundance_plot_html,
)
from repeatviewer.backends.msa_plot import (
    MsaPlotOptions,
    msa_plot_output_path,
    write_msa_plot,
)
from repeatviewer.backends.job_cancel import CANCELLED_MESSAGE, JobCancellation, JobCancelledError

logger = logging.getLogger(__name__)
operation_logger = logging.getLogger("repeatviewer.operations")
from repeatviewer.backends.te_benchmark import (
    default_tetrimmer_home,
    require_bedtools,
    run_confusion_matrix_benchmark,
)
from repeatviewer.backends.te_library_benchmark import run_te_library_benchmark
from repeatviewer.backends.flexidot_go_to import (
    fasta_for_flexidot_html,
    parse_dotplot_go_to_targets,
)
from repeatviewer.backends.curation import (
    effective_batch_workers,
    is_flexidot_extra_plot_file,
    is_flexidot_html_file,
    run_batch,
)
from repeatviewer.backends.gui_functions import flexidot_executable
from repeatviewer.backends.gui_functions import sanitize_name as flexidot_sanitize_name
from repeatviewer.backends.sequence_aligner import (
    alignment_backends,
    run_alignment,
    run_alignment_combined,
)
from repeatviewer.backends.bismark_cov import is_bismark_cov_file
from repeatviewer.backends.gff_features import (
    build_gff_annotation_html,
    find_companion_gff,
    is_gene_annotation_file,
)
from repeatviewer.backends.sequence_features import (
    companion_annotation_files_for_delete,
    is_sequence_map_hidden_file,
)
from repeatviewer.plugins import get_registry, load_plugins
from repeatviewer.plugins.custom_command import (
    UserCommandPlugin,
    UserCommandRunDialog,
    delete_user_command_plugin_spec,
    existing_output_files,
    is_file_var_mode,
    is_folder_var_mode,
    is_path_var_mode,
    load_user_command_plugin_specs,
    normalize_var_mode,
    save_user_command_plugin_spec,
)
from repeatviewer.backends.sequence_viewer import open_fasta
from repeatviewer.job_options import (
    ALIGNER_FASTGA,
    ALIGNER_MAFFT,
    ALIGNER_MUSCLE,
    ALIGNERS,
    FASTGA_MODE_PAIRS,
    FASTGA_MODE_REFERENCE,
    ALIGN_MODE_COMBINED,
    BLAST_PROGRAM_BLASTN,
    BLAST_PROGRAM_BLASTP,
    BLAST_PROGRAM_BLASTX,
    EXTEND_RIGHT_CLICK_DIALOG,
    FASTGA_RIGHT_CLICK_DIALOG,
    FLEXIDOT_RIGHT_CLICK_DIALOG,
    CONSENSUS_RIGHT_CLICK_DIALOG,
    TEAID_RIGHT_CLICK_DIALOG,
    MSA_CLEAN_RIGHT_CLICK_DIALOG,
    TE_ANNOTATION_RIGHT_CLICK_DIALOG,
    TREE_RIGHT_CLICK_DIALOG,
    load_align_options,
    load_blast_options,
    load_consensus_options,
    load_extend_options,
    load_fastga_options,
    load_flexidot_options,
    load_msa_clean_options,
    load_te_annotation_options,
    load_teaid_options,
    load_te_genome_map_options,
    load_te_divergence_plot_options,
    load_msa_plot_options,
    load_repeat_quantify_options,
    load_te_confusion_matrix_benchmark_options,
    load_te_library_benchmark_options,
    load_tsd_search_options,
    load_tree_options,
    load_map_options,
    msa_clean_method_description,
    save_align_options,
    save_blast_options,
    save_consensus_options,
    save_extend_options,
    save_fastga_options,
    save_flexidot_options,
    save_msa_clean_options,
    save_te_annotation_options,
    save_teaid_options,
    save_te_genome_map_options,
    save_te_divergence_plot_options,
    save_msa_plot_options,
    save_repeat_quantify_options,
    save_te_confusion_matrix_benchmark_options,
    save_te_library_benchmark_options,
    save_tsd_search_options,
    save_tree_options,
    save_map_options,
    MAP_RIGHT_CLICK_DIALOG,
    TeGenomeMapOptions,
    TeDivergencePlotOptions,
    MsaPlotJobOptions,
    TeConfusionMatrixBenchmarkOptions,
    TeLibraryBenchmarkOptions,
    TsdSearchOptions,
)
from repeatviewer.views.tool_dialogs import (
    AlignParamsDialog,
    BlastParamsDialog,
    ConsensusParamsDialog,
    ExtendParamsDialog,
    FlexidotParamsDialog,
    TeAnnotationDialog,
    batch_options_from_te_annotation_job,
    OrfPredictionParamsDialog,
    MsaCleanParamsDialog,
    TeaAidParamsDialog,
    TsdSearchParamsDialog,
    TeDensityPlotDialog,
    TeDivergencePlotDialog,
    TeAbundancePlotDialog,
    MsaPlotDialog,
    ConfusionMatrixBenchmarkDialog,
    TeLibraryBenchmarkDialog,
    CommandPluginDesignerDialog,
    WrapperManageDialog,
    exec_modal_dialog,
)
from repeatviewer.views.plugin_manage_dialog import InstalledPluginManageDialog
from repeatviewer.backends.separate_fasta import (
    count_fasta_sequences,
    sanitize_name,
    separate_sequences,
    separation_output_dir,
)
from repeatviewer.views.genome_dialog import GenomeDialog, is_genome_fasta
from repeatviewer.views.create_folder_dialog import (
    DEFAULT_FOLDER_NAME,
    CreateFolderDialog,
    common_prefix_name,
)
from repeatviewer.views.define_database_dialog import DefineDatabaseDialog, resolve_panel_database_path
from repeatviewer.views.settings_dialog import (
    SettingsDialog,
    add_genome_to_history,
    load_genome_history,
    load_last_folder,
    load_expanded_folders,
    load_root_folders,
    save_workspace,
)
from repeatviewer.views.toolbar_dialog import (
    DEFAULT_TOOLBAR_ORDER,
    load_toolbar_font_size,
    load_toolbar_order,
)
from repeatviewer.widgets.browser_icons import (
    browser_file_icon,
    browser_folder_icon,
    browser_folder_spacer_icon,
)
from repeatviewer.widgets.browser_dock import PreviewDockWidget
from repeatviewer.widgets.browser_panel_settings import load_browser_panel_font_size
from repeatviewer.ui.platform import reveal_in_file_manager
from repeatviewer.ui.text_editor import open_in_text_editor
from repeatviewer.ui.theme import application_stylesheet
from repeatviewer.widgets.hidden_folders import (
    hide_root_folder,
    load_hidden_root_folders,
    remove_hidden_root,
    replace_hidden_root_path,
    unhide_root_folder,
)
from repeatviewer.widgets.folder_colors import (
    FOLDER_COLOR_PRESETS,
    apply_folder_item_color,
    apply_folder_colors_to_tree,
    get_folder_color,
    remove_folder_colors_under,
    replace_folder_color_paths,
    set_folder_color,
)
from repeatviewer.widgets.folder_navigation import FolderNavigationBar
from repeatviewer.widgets.folder_tree_delegate import FolderTreeItemDelegate
from repeatviewer.widgets.file_sort import (
    FileSortMode,
    load_file_sort_mode,
    sort_file_paths,
)
from repeatviewer.widgets.file_list_table import FileListRow, FileListTableWidget
from repeatviewer.widgets.preview_tabs import PreviewTabWidget, is_previewable_file, preview_path_priority
from repeatviewer.widgets.task_progress import TaskProgressPanel
from repeatviewer.widgets.database_index_runner import DatabaseIndexRunner

ROLE_PATH = Qt.ItemDataRole.UserRole
ROLE_LAZY_PLACEHOLDER = Qt.ItemDataRole.UserRole + 1

_BROWSER_ICON_SIZE = 16


def _unique_duplicate_path(path: Path) -> Path:
    """Next unused path in the same folder: ``name_copy.ext``, ``name_copy2.ext``, …"""
    parent = path.parent
    candidate = parent / f"{path.stem}_copy{path.suffix}"
    if not candidate.exists():
        return candidate
    n = 2
    while True:
        candidate = parent / f"{path.stem}_copy{n}{path.suffix}"
        if not candidate.exists():
            return candidate
        n += 1


@dataclass
class _FileClipboard:
    mode: Literal["copy", "cut"]
    paths: list[Path]


_NEW_OUTPUT_BG_LIGHT = QColor(214, 245, 214)
_NEW_OUTPUT_FG_LIGHT = QColor(0, 100, 0)
_NEW_OUTPUT_BG_DARK = QColor(36, 78, 48)
_NEW_OUTPUT_FG_DARK = QColor(180, 235, 190)
_SOURCE_BG_LIGHT = QColor(220, 235, 255)
_SOURCE_FG_LIGHT = QColor(0, 70, 140)
_SOURCE_BG_DARK = QColor(35, 55, 90)
_SOURCE_FG_DARK = QColor(170, 205, 255)

_FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".faa"}

# Extract Sequences also accepts the genome-style suffixes and gzipped input;
# it streams rather than parsing into memory, so size is not a reason to refuse.
_EXTRACT_SOURCE_SUFFIXES = _FASTA_SUFFIXES | {".fna", ".ffn"}
_GZIP_SUFFIXES = {".gz", ".bgz"}
# Tick-box list is fine for genomes (dozens of chroms); skip only for huge
# read sets where building thousands of QListWidget items is unusable.
_EXTRACT_PICK_LIST_MAX_RECORDS = 10_000
# Confirm before scattering this many single-sequence files into a folder.
_EXTRACT_MANY_FILES_WARNING = 200
# sseqid in a RepeatViewer BLAST table: the genome sequence that was hit.
_BLAST_SUBJECT_COLUMN = 2


def _is_extract_source(path: Path) -> bool:
    """True for a FASTA the Extract Sequences tool can read, gzipped or not."""
    suffixes = [s.lower() for s in path.suffixes]
    if suffixes and suffixes[-1] in _GZIP_SUFFIXES:
        suffixes = suffixes[:-1]
    return bool(suffixes) and suffixes[-1] in _EXTRACT_SOURCE_SUFFIXES

_TOOLBAR_TOOLTIPS = {
    "define_genome": "Set or change the project reference genome FASTA",
    "new_sequence": "Create a new empty sequence FASTA",
    "group_sequences": "Combine selected sequences into one FASTA file",
    "extract_sequences": "Extract selected records from a FASTA file",
    "reverse_complement": "Reverse-complement selected FASTA sequences",
    "orf_prediction": "Predict ORFs in selected FASTA sequences",
    "translate": "Translate selected DNA FASTA to protein",
    "aliview": "Open selected FASTA in AliView",
    "align": "Align FASTA using default parameters",
    "fastga": "FastGA whole-genome alignment (self-comparison one-by-one by default; optional all-pairs)",
    "tree": "Phylogenetic tree from the selected FASTA (MAFFT + IQ-TREE / FastME / FastTree)",
    "map_reference": "Map the selected sequences or reads to a reference (MAFFT / minimap2 / Bowtie2)",
    "blast": "BLAST selected FASTA against the project genome using default parameters",
    "te_genome_map": "TE density plot using default parameters",
    "te_divergence_plot": "TE divergence plot (divergence vs total TE length)",
    "extend": "Extend MSA ends using default parameters",
    "teaid": "TE-Aid interactive HTML report using default parameters",
    "msa_clean": "Trim or clean MSA using default parameters",
    "flexidot": "Dotplot comparison of selected sequences using default parameters",
    "consensus": "Click: set consensus parameters and run. Right-click: run with saved settings.",
    "te_annotation": "Run sequence annotation on selected single-sequence FASTA files",
}

_JOB_DISPLAY_TITLES = {
    "teaid": "TE-Aid",
    "blast": "BLAST",
    "extend": "MSA Extend",
    "msa_clean": "MSA Clean",
    "msa_clean_gappy": "MSA remove gappy columns",
    "msa_clean_div": "Clean row with divergence",
    "msa_clean_gap": "Clean row with gap",
    "consensus_generation": "Consensus",
    "flexidot": "Flexidot",
    "te_annotation": "Sequence Annotation",
    "reverse_complement": "Reverse complement",
    "orf_prediction": "ORF prediction",
    "dna_translation": "DNA translation",
    "te_abundance_plot": "TE abundance plot",
}


def _format_cli_failure(exc: BaseException) -> str:
    """Prefer stderr/command detail from failed external tools."""
    import subprocess

    if isinstance(exc, subprocess.CalledProcessError):
        detail = (exc.stderr or exc.output or "").strip()
        if detail:
            return detail
    return str(exc)


def _mafft_option_kwargs(options) -> dict[str, Any]:
    return {
        "num_threads": options.num_threads,
        "strategy": options.strategy,
        "sequence_type": options.sequence_type,
        "retree": options.retree,
        "max_iterate": options.max_iterate,
        "gap_open_penalty": options.gap_open_penalty,
        "offset": options.offset,
        "reorder": options.reorder,
        "adjust_direction": options.adjust_direction,
        "muscle_mode": options.muscle_mode,
        "muscle_num_threads": options.muscle_num_threads,
    }


class AlignWorker(QObject):
    progress = Signal(int, str)
    log_line = Signal(str)
    finished = Signal(list)
    failed = Signal(str)

    def __init__(
        self,
        seq_files: list[Path],
        *,
        options,
        mode: str = "one_by_one",
        max_workers: int = 4,
        use_system_mafft: bool = False,
        use_system_muscle: bool = False,
        engine: str = "mafft",
        stream_log: bool = False,
        cancel=None,
    ):
        super().__init__()
        self._seq_files = seq_files
        self._options = options
        self._mode = mode
        self._max_workers = max_workers
        self._use_system_mafft = use_system_mafft
        self._use_system_muscle = use_system_muscle
        self._engine = engine
        self._stream_log = bool(stream_log)
        self._cancel = cancel

    @Slot()
    def run(self) -> None:
        from concurrent.futures import ThreadPoolExecutor, as_completed

        from repeatviewer.external_process import process_log_sink

        def _log(message: str) -> None:
            self.log_line.emit(str(message))

        sink = _log  # always capture; UI shows live and/or on failure
        cancel = self._cancel

        if self._mode == ALIGN_MODE_COMBINED:
            try:
                with process_log_sink(sink):
                    if cancel is not None:
                        cancel.check()
                    self.progress.emit(0, "Aligning combined MSA…")
                    _log("Aligning combined MSA…")
                    output = run_alignment_combined(
                        self._seq_files,
                        **_mafft_option_kwargs(self._options),
                        use_system_mafft=self._use_system_mafft,
                        use_system_muscle=self._use_system_muscle,
                        engine=self._engine,
                        cancel=cancel,
                    )
                    self.progress.emit(1, "Alignment finished.")
                    _log("Alignment finished.")
                self.finished.emit([str(output)])
            except JobCancelledError:
                self.failed.emit(CANCELLED_MESSAGE)
            except Exception as exc:  # noqa: BLE001
                self.failed.emit(_format_cli_failure(exc))
            return

        total = len(self._seq_files)
        created: list[str] = []
        failures: list[str] = []
        workers = effective_batch_workers(self._max_workers, total)
        completed = 0

        def _align_one(seq_path: Path) -> tuple[Path, str | None, str | None]:
            try:
                if cancel is not None:
                    cancel.check()
                with process_log_sink(sink):
                    output = run_alignment(
                        seq_path,
                        **_mafft_option_kwargs(self._options),
                        use_system_mafft=self._use_system_mafft,
                        use_system_muscle=self._use_system_muscle,
                        engine=self._engine,
                        cancel=cancel,
                    )
                return seq_path, str(output), None
            except JobCancelledError:
                raise
            except Exception as exc:  # noqa: BLE001
                return seq_path, None, _format_cli_failure(exc)

        try:
            if workers == 1:
                for seq_path in self._seq_files:
                    if cancel is not None:
                        cancel.check()
                    self.progress.emit(completed, f"Aligning {seq_path.name}…")
                    _log(f"—— Aligning {seq_path.name} ——")
                    _, out_path, err = _align_one(seq_path)
                    completed += 1
                    if err:
                        failures.append(f"{seq_path.name}: {err}")
                    elif out_path:
                        created.append(out_path)
            else:
                with ThreadPoolExecutor(max_workers=workers) as pool:
                    futures = [pool.submit(_align_one, p) for p in self._seq_files]
                    for future in as_completed(futures):
                        if cancel is not None and cancel.requested:
                            for pending in futures:
                                pending.cancel()
                            raise JobCancelledError("Job stopped by user")
                        seq_path, out_path, err = future.result()
                        self.progress.emit(completed, f"Aligning {seq_path.name}…")
                        _log(f"Finished {seq_path.name}")
                        completed += 1
                        if err:
                            failures.append(f"{seq_path.name}: {err}")
                        elif out_path:
                            created.append(out_path)
            self.progress.emit(total, "Alignment finished.")
            _log("Alignment finished.")
            if failures and not created:
                self.failed.emit("\n".join(failures[:8]))
                return
            if failures:
                self.finished.emit(created)
                return
            self.finished.emit(created)
        except JobCancelledError:
            self.failed.emit(CANCELLED_MESSAGE)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(_format_cli_failure(exc))


class FastGaWorker(QObject):
    """Run FastGA alignment(s) off the UI thread.

    Modes (FastGaOptions.mode): each file against itself, all pairs into one
    file, or each file against a reference (indexed once). The loop lives in
    ``backends.fastga_dotplot.run_fastga_batch``. Indexing a large genome takes
    minutes, so progress is reported and cancellation is polled throughout.
    """

    progress = Signal(int, str)
    log_line = Signal(str)
    finished = Signal(list)
    failed = Signal(str)

    def __init__(
        self,
        seq_files: list[Path],
        options,
        *,
        cancel=None,
        stream_log: bool = False,
    ) -> None:
        super().__init__()
        self._seq_files = list(seq_files)
        self._options = options
        self._cancel = cancel
        self._stream_log = bool(stream_log)

    @Slot()
    def run(self) -> None:
        from repeatviewer.backends.fastga_dotplot import run_fastga_batch

        try:
            step = [1]

            def _on_step(n: int) -> None:
                step[0] = max(1, int(n))

            def _progress(message: str) -> None:
                self.progress.emit(step[0], message)

            def _log(message: str) -> None:
                self.log_line.emit(str(message))

            def _cancelled() -> bool:
                return bool(self._cancel is not None and self._cancel.requested)

            created = run_fastga_batch(
                self._seq_files,
                self._options,
                progress=_progress,
                on_log=_log,
                is_cancelled=_cancelled,
                on_step=_on_step,
            )
            self.finished.emit([str(p) for p in created])
        except Exception as exc:  # noqa: BLE001 - surface to the UI
            self.failed.emit(_format_cli_failure(exc))


class TreeWorker(QObject):
    """Build a phylogenetic tree off the UI thread (Tools → Tree)."""

    progress = Signal(int, str)
    log_line = Signal(str)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, seq_files: list[Path], options, *, cancel=None) -> None:
        super().__init__()
        self._seq_files = list(seq_files)
        self._options = options
        self._cancel = cancel

    @Slot()
    def run(self) -> None:
        from repeatviewer.backends.tree_builder import TreeBuildError, build_tree_run
        from repeatviewer.external_process import process_log_sink

        def _log(message: str) -> None:
            self.log_line.emit(str(message))

        steps = {"Aligning": 1, "Building": 2, "Drawing": 3}

        def _progress(message: str) -> None:
            step = next((v for k, v in steps.items() if message.startswith(k)), 1)
            self.progress.emit(step, message)
            self.log_line.emit(message)

        try:
            with process_log_sink(_log):
                run_dir, files, report = build_tree_run(
                    self._seq_files, self._options, progress=_progress, cancel=self._cancel
                )
            self.finished.emit({"details_dir": str(run_dir), "files": [str(f) for f in files],
                                "report": str(report)})
        except JobCancelledError:
            self.failed.emit(CANCELLED_MESSAGE)
        except TreeBuildError as exc:
            self.failed.emit(str(exc))
        except Exception as exc:  # noqa: BLE001 - surface to the UI
            self.failed.emit(_format_cli_failure(exc))


class MapReferenceWorker(QObject):
    """Map sequences or reads to a reference off the UI thread (Tools → Map to Reference)."""

    progress = Signal(int, str)
    log_line = Signal(str)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(self, seq_files: list[Path], options, *, cancel=None) -> None:
        super().__init__()
        self._seq_files = list(seq_files)
        self._options = options
        self._cancel = cancel

    @Slot()
    def run(self) -> None:
        from repeatviewer.backends.map_reference import MapReferenceError, run_mapping
        from repeatviewer.external_process import process_log_sink

        def _log(message: str) -> None:
            self.log_line.emit(str(message))

        steps = {"Reading": 1, "Indexing": 1, "Placing": 2, "Mapping": 2, "Collecting": 3,
                 "Consensus": 3, "Writing": 4}

        def _progress(message: str) -> None:
            step = next((v for k, v in steps.items() if message.startswith(k)), 1)
            self.progress.emit(step, message)
            self.log_line.emit(message)

        try:
            with process_log_sink(_log):
                rvmap, files, report = run_mapping(
                    self._seq_files, self._options, progress=_progress, cancel=self._cancel
                )
            self.finished.emit({"rvmap": str(rvmap), "files": [str(f) for f in files],
                                "report": str(report)})
        except JobCancelledError:
            self.failed.emit(CANCELLED_MESSAGE)
        except MapReferenceError as exc:
            self.failed.emit(str(exc))
        except Exception as exc:  # noqa: BLE001 - surface to the UI
            self.failed.emit(_format_cli_failure(exc))


class RepeatQuantifyWorker(QObject):
    """Quantify a repeat library across target sequences off the UI thread."""

    progress = Signal(int, str)
    finished = Signal(object)
    failed = Signal(str)
    # A new saved BLAST database (BlastDatabase) built for a genome; the main
    # window registers it so Add/Manage databases and other tools reuse it.
    database_created = Signal(object)

    def __init__(
        self, library: Path, targets: list[Path], *, engine: str, threads: int,
        params: dict | None = None, cancel=None,
        output_base: Path | None = None, stamp: str | None = None,
        blast_databases: list | None = None,
        save_blast_index: bool = False,
        blast_db_base_dir: Path | None = None,
        use_system_blast: bool = False,
    ) -> None:
        super().__init__()
        from repeatviewer.backends.repeat_quantify import run_stamp

        self._library = Path(library)
        self._targets = list(targets)
        self._engine = engine
        self._threads = threads
        self._params = params or {}
        self._cancel = cancel
        # Saved BLAST databases to reuse; new ones are added as they are built.
        self._blast_databases = list(blast_databases or [])
        self._save_blast_index = bool(save_blast_index)
        self._blast_db_base_dir = Path(blast_db_base_dir) if blast_db_base_dir else None
        self._use_system_blast = bool(use_system_blast)
        # None = a new folder next to each genome; otherwise new folders inside it.
        self._output_base = Path(output_base) if output_base else None
        # One stamp for the whole batch, taken when Run was clicked.
        self._stamp = stamp or run_stamp()

    @Slot()
    def run(self) -> None:
        from datetime import datetime

        from repeatviewer.backends.repeat_quantify import (
            RUN_PARAMETERS_FILENAME,
            create_repeat_quant_output_dir,
            library_duplicate_warnings,
            quantify_library,
            resolve_blast_database,
            write_landscape_html,
            write_run_parameters,
            write_summary_html,
        )

        out_dir: Path | None = None
        run_info: dict = {}

        def _record(status: str, error: str = "") -> None:
            if out_dir is None:
                return
            try:
                files = [
                    p for p in out_dir.rglob("*")
                    if p.is_file() and p.name != RUN_PARAMETERS_FILENAME
                ]
                write_run_parameters(
                    out_dir, **run_info, status=status,
                    finished=datetime.now(), outputs=files, error=error,
                )
            except OSError:
                logger.exception("Could not write %s", RUN_PARAMETERS_FILENAME)

        try:
            lib_warnings = library_duplicate_warnings(self._library)
            outputs: list[str] = []
            for index, target in enumerate(self._targets, 1):
                self.progress.emit(index, f"Quantifying {target.name}…")
                out_dir = create_repeat_quant_output_dir(
                    target, self._stamp, self._output_base
                )
                params = dict(self._params)
                if self._engine in ("blastn", "both"):
                    db_prefix = self._blast_database_for(target, index, resolve_blast_database)
                    if db_prefix is not None:
                        params["blast_db_prefix"] = str(db_prefix)
                run_info = dict(
                    library_fasta=self._library,
                    target_fasta=target,
                    engine=self._engine,
                    threads=self._threads,
                    params=params,
                    started=datetime.now(),
                )
                write_run_parameters(out_dir, **run_info)
                results = quantify_library(
                    self._library,
                    target,
                    out_dir,
                    engine=self._engine,
                    threads=self._threads,
                    use_system_blast=self._use_system_blast,
                    cancel=self._cancel,
                    **params,
                )
                summary = write_summary_html(
                    results,
                    out_dir / f"{target.name}_repeat_quant.html",
                    target_label=target.name,
                    warnings=lib_warnings,
                    library_fasta=self._library,
                    run_info=params,
                )
                outputs.append(str(summary))
                for result in results:
                    if result.engine == "repeatmasker" and result.landscape:
                        outputs.append(
                            str(
                                write_landscape_html(
                                    result,
                                    out_dir / f"{target.name}_repeat_landscape.html",
                                    target_label=target.name,
                                )
                            )
                        )
                _record("completed")
                out_dir = None
            self.finished.emit(outputs)
        except JobCancelledError:
            _record("cancelled")
            self.finished.emit([])  # stopped by user; end quietly
        except Exception as exc:  # noqa: BLE001 - surface to the UI
            message = _format_cli_failure(exc)
            _record("failed", message)
            self.failed.emit(message)

    def _blast_database_for(self, target: Path, index: int, resolve) -> Path | None:
        """Saved BLAST database prefix for ``target``, building one if asked."""
        from repeatviewer.backends.database_prepare import saved_blast_database_for

        if (
            self._save_blast_index
            and saved_blast_database_for(target, self._blast_databases) is None
        ):
            self.progress.emit(index, f"Indexing {target.name} as a BLAST database…")
        prefix, created = resolve(
            target,
            self._blast_databases,
            save_index=self._save_blast_index,
            base_dir=self._blast_db_base_dir,
            use_system_blast=self._use_system_blast,
        )
        if created is not None:
            self._blast_databases.append(created)
            self.database_created.emit(created)
        if self._cancel is not None:
            self._cancel.check()  # Stop pressed while makeblastdb ran
        return prefix


class GenomeMapWorker(QObject):
    """Build TE genome map HTML off the UI thread."""

    progress = Signal(int, str)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(
        self,
        te_input_files: list[Path],
        genome_fasta: Path,
        *,
        job_kwargs: dict | None = None,
    ):
        super().__init__()
        self._te_input_files = [p.expanduser().resolve() for p in te_input_files]
        self._genome_fasta = genome_fasta.expanduser().resolve()
        self._job_kwargs = job_kwargs or {}

    @Slot()
    def run(self) -> None:
        try:
            html_path = write_blast_genome_map_html(
                self._te_input_files,
                genome_fasta=self._genome_fasta,
                progress=self.progress.emit,
                **self._job_kwargs,
            )
            if html_path is None:
                self.failed.emit("Could not build the genome map.")
                return
            self.finished.emit(html_path)
        except JobCancelledError:
            self.failed.emit(CANCELLED_MESSAGE)
        except GenomeMapNothingToPlotError as exc:
            self.failed.emit(_format_cli_failure(exc))
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(_format_cli_failure(exc))


class ExtractSequencesWorker(QObject):
    """Extract (or just count) FASTA records named in an ID list, off the UI thread."""

    progress = Signal(int, str)
    finished = Signal(object)
    failed = Signal(str)

    def __init__(
        self,
        source: Path,
        ids: list[str],
        out_dir: Path,
        *,
        count_only: bool = False,
        single_file: bool = True,
        exclude: bool = False,
        output_path: Path | None = None,
        cancel=None,
    ):
        super().__init__()
        self._source = Path(source).expanduser().resolve()
        self._ids = list(ids)
        self._out_dir = Path(out_dir).expanduser().resolve()
        self._count_only = count_only
        self._single_file = single_file
        self._exclude = exclude
        self._output_path = (
            Path(output_path).expanduser().resolve() if output_path is not None else None
        )
        self._cancel = cancel

    def _emit_progress(self, done: int, total: int) -> None:
        if total > 0:
            percent = max(0, min(100, int(done * 100 / total)))
            self.progress.emit(percent, f"{percent}%")
        else:
            # Compressed input: the uncompressed size is unknown, so report
            # what has actually been read rather than a made-up percentage.
            self.progress.emit(0, f"{done / 1_000_000:.0f} MB read")

    @Slot()
    def run(self) -> None:
        from repeatviewer.backends.sequence_compose import (
            SequenceComposeError,
            count_matching_ids,
            extract_records_by_ids,
        )

        try:
            if self._count_only:
                summary = count_matching_ids(
                    self._source,
                    self._ids,
                    progress=self._emit_progress,
                    cancel=self._cancel,
                )
            else:
                summary = extract_records_by_ids(
                    self._source,
                    self._ids,
                    self._out_dir,
                    single_file=self._single_file,
                    exclude=self._exclude,
                    output_path=self._output_path,
                    progress=self._emit_progress,
                    cancel=self._cancel,
                )
            self.finished.emit(summary)
        except SequenceComposeError as exc:
            self.failed.emit(str(exc))
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(_format_cli_failure(exc))


class ConfusionMatrixBenchmarkWorker(QObject):
    """Run TEtrimmer confusion-matrix benchmark for one or more query .out files."""

    progress = Signal(int, str)
    finished = Signal(list)
    failed = Signal(str)

    def __init__(
        self,
        query_files: list[Path],
        *,
        reference_out: Path,
        genome_fasta: Path,
        tetrimmer_home: Path,
        less_classification: bool = False,
        debug: bool = False,
        cancel: JobCancellation | None = None,
    ) -> None:
        super().__init__()
        self._query_files = [p.expanduser().resolve() for p in query_files]
        self._reference_out = reference_out.expanduser().resolve()
        self._genome_fasta = genome_fasta.expanduser().resolve()
        self._tetrimmer_home = tetrimmer_home.expanduser().resolve()
        self._less_classification = less_classification
        self._debug = debug
        self._cancel = cancel

    @Slot()
    def run(self) -> None:
        try:
            created: list[Path] = []
            total = len(self._query_files)
            for index, query in enumerate(self._query_files, start=1):
                if self._cancel is not None and self._cancel.requested:
                    raise JobCancelledError("Job stopped by user")
                self.progress.emit(index - 1, f"Benchmarking: {query.name}…")
                outputs = run_confusion_matrix_benchmark(
                    query_out=query,
                    reference_out=self._reference_out,
                    genome_fasta=self._genome_fasta,
                    output_dir=query.parent,
                    tetrimmer_home=self._tetrimmer_home,
                    less_classification=self._less_classification,
                    debug=self._debug,
                )
                created.extend(outputs)
                self.progress.emit(index, f"Finished: {query.name}")
            self.finished.emit(created)
        except JobCancelledError:
            self.failed.emit(CANCELLED_MESSAGE)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(_format_cli_failure(exc))


class TeLibraryBenchmarkWorker(QObject):
    """Run TE consensus-library benchmark (Perfect / Good / Present)."""

    progress = Signal(int, str)
    finished = Signal(list)
    failed = Signal(str)

    def __init__(
        self,
        *,
        query_fasta: Path,
        standard_fasta: Path,
        rm_out: Path | None = None,
        run_repeatmasker: bool = True,
        threads: int = 4,
        cancel: JobCancellation | None = None,
    ) -> None:
        super().__init__()
        self._query_fasta = query_fasta.expanduser().resolve()
        self._standard_fasta = standard_fasta.expanduser().resolve()
        self._rm_out = rm_out.expanduser().resolve() if rm_out else None
        self._run_repeatmasker = run_repeatmasker
        self._threads = max(1, int(threads))
        self._cancel = cancel

    @Slot()
    def run(self) -> None:
        try:
            if self._cancel is not None and self._cancel.requested:
                raise JobCancelledError("Job stopped by user")
            self.progress.emit(0, f"Library benchmark: {self._query_fasta.name}…")
            outputs = run_te_library_benchmark(
                query_fasta=self._query_fasta,
                standard_fasta=self._standard_fasta,
                rm_out=self._rm_out,
                run_repeatmasker=self._run_repeatmasker,
                threads=self._threads,
            )
            self.progress.emit(1, "Finished")
            self.finished.emit(outputs)
        except JobCancelledError:
            self.failed.emit(CANCELLED_MESSAGE)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(_format_cli_failure(exc))


class FunctionWorker(QObject):
    progress = Signal(int, str)
    log_line = Signal(str)
    finished = Signal(dict)
    failed = Signal(str)

    def __init__(
        self,
        run_func,
        func_name: str,
        seq_files: list[Path],
        job_kwargs: dict | None = None,
        *,
        stream_log: bool = False,
    ):
        super().__init__()
        self._run_func = run_func
        self._func_name = func_name
        self._seq_files = seq_files
        self._job_kwargs = job_kwargs or {}
        self._stream_log = bool(stream_log)

    @Slot()
    def run(self) -> None:
        from repeatviewer.external_process import process_log_sink

        def _log(message: str) -> None:
            self.log_line.emit(str(message))

        try:
            with process_log_sink(_log):
                result = self._run_func(
                    self._func_name,
                    self._seq_files,
                    self.progress.emit,
                    **self._job_kwargs,
                )
            self.finished.emit(result)
        except JobCancelledError:
            self.failed.emit(CANCELLED_MESSAGE)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(_format_cli_failure(exc))


class WorkflowWorker(QObject):
    """Run a multi-step tool workflow on a background thread."""

    progress = Signal(int, str)
    log_line = Signal(str)
    finished = Signal(dict)
    failed = Signal(str)

    def __init__(self, run_func, job_kwargs: dict | None = None):
        super().__init__()
        self._run_func = run_func
        self._job_kwargs = job_kwargs or {}

    @Slot()
    def run(self) -> None:
        from repeatviewer.external_process import process_log_sink

        def _log(message: str) -> None:
            self.log_line.emit(str(message))

        try:
            with process_log_sink(_log):
                result = self._run_func(self.progress.emit, **self._job_kwargs)
            self.finished.emit(result)
        except JobCancelledError:
            self.failed.emit(CANCELLED_MESSAGE)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(_format_cli_failure(exc))


class FolderTreeWidget(QTreeWidget):
    """Folder-only navigator; accepts dropped folders/files."""

    folder_selected = Signal(Path)
    files_dropped_to_folder = Signal(Path, list)
    folder_dropped = Signal(Path)
    rename_folder_requested = Signal(Path)
    delete_folders_requested = Signal(list)
    clear_root_folders_requested = Signal(list)
    paste_to_folder_requested = Signal(Path)
    add_folder_requested = Signal()
    create_folder_requested = Signal(object)  # Path | None

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setHeaderHidden(True)
        # Cmd/Ctrl+click and Shift+click select multiple folders (delete, color, …).
        self.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self.setDragDropMode(QAbstractItemView.DragDropMode.DropOnly)
        self.setAcceptDrops(True)
        self.setDropIndicatorShown(True)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self.setHorizontalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.setVerticalScrollMode(QAbstractItemView.ScrollMode.ScrollPerPixel)
        self.horizontalScrollBar().setSingleStep(12)
        self.verticalScrollBar().setSingleStep(12)
        self.setTextElideMode(Qt.TextElideMode.ElideMiddle)
        self.setWordWrap(False)
        self.setUniformRowHeights(False)
        self.setAnimated(False)
        self.setIconSize(QSize(_BROWSER_ICON_SIZE, _BROWSER_ICON_SIZE))
        self.setIndentation(18)
        self.setRootIsDecorated(False)
        self.setMouseTracking(True)
        self.viewport().setMouseTracking(True)
        self.setItemDelegate(FolderTreeItemDelegate(self))
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
        self.customContextMenuRequested.connect(self._show_context_menu)
        header = self.header()
        # Size the column to content so deep/long names can scroll horizontally
        # instead of being clipped by a stretched viewport-width column.
        header.setStretchLastSection(False)
        header.setSectionResizeMode(0, QHeaderView.ResizeMode.Interactive)
        header.setMinimumSectionSize(100)
        self.currentItemChanged.connect(self._on_folder_current_changed)
        self.itemExpanded.connect(lambda *_: self.refresh_horizontal_extent())
        self.itemCollapsed.connect(lambda *_: self.refresh_horizontal_extent())
        self._suppress_current_changed = False
        # True while Qt is handling a click / key press: Qt applies the
        # selection itself then, so set_active_folder must not touch it.
        self._user_selecting = False
        self._active_folder_key: str | None = None
        self._drop_target_item: QTreeWidgetItem | None = None

    def drop_target_item(self) -> QTreeWidgetItem | None:
        """Folder row currently under an in-progress file drop, if any."""
        return self._drop_target_item

    def _set_drop_target_item(self, item: QTreeWidgetItem | None) -> None:
        if item is self._drop_target_item:
            return
        previous = self._drop_target_item
        self._drop_target_item = item
        viewport = self.viewport()
        if viewport is None:
            return
        for row in (previous, item):
            if row is None:
                continue
            rect = self.visualItemRect(row)
            if rect.isValid():
                viewport.update(rect.adjusted(-2, -2, 2, 2))

    def _folder_item_at(self, pos) -> QTreeWidgetItem | None:
        item = self.itemAt(pos)
        if item is None or self._is_lazy_placeholder(item):
            return None
        raw = item.data(0, ROLE_PATH)
        if not raw:
            return None
        return item

    @staticmethod
    def _mime_local_paths(mime) -> list[Path]:
        if mime is None or not mime.hasUrls():
            return []
        paths: list[Path] = []
        for url in mime.urls():
            if not url.isLocalFile():
                continue
            text = url.toLocalFile()
            if text:
                paths.append(Path(text))
        return paths

    def wheelEvent(self, event: QWheelEvent) -> None:  # noqa: N802
        """Scroll the folder tree more slowly than Qt's default trackpad/wheel speed."""
        # Keep horizontal and vertical motion; Shift+vertical becomes horizontal.
        pixel = event.pixelDelta()
        angle = event.angleDelta()
        if pixel.manhattanLength() > 0:
            dx, dy = float(pixel.x()), float(pixel.y())
        else:
            # angleDelta is 1/8 degree units; 120 ≈ one notch.
            dx, dy = float(angle.x()) / 8.0, float(angle.y()) / 8.0
        if (
            event.modifiers() & Qt.KeyboardModifier.ShiftModifier
            and abs(dx) < 0.5
            and abs(dy) >= 0.5
        ):
            dx, dy = dy, 0.0

        scale = 0.35
        moved = False
        if abs(dx) >= 0.5:
            hbar = self.horizontalScrollBar()
            hbar.setValue(hbar.value() - int(round(dx * scale)))
            moved = True
        if abs(dy) >= 0.5:
            vbar = self.verticalScrollBar()
            vbar.setValue(vbar.value() - int(round(dy * scale)))
            moved = True
        if moved:
            event.accept()
            return
        super().wheelEvent(event)

    def refresh_horizontal_extent(self) -> None:
        """Widen column to fit visible labels; keep at least the viewport width."""
        self.resizeColumnToContents(0)
        viewport_w = self.viewport().width() if self.viewport() is not None else 0
        if viewport_w > 0 and self.columnWidth(0) < viewport_w:
            self.setColumnWidth(0, viewport_w)

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self.refresh_horizontal_extent()

    def drawBranches(self, painter: QPainter, rect: QRect, index) -> None:  # noqa: N802
        """Skip stock branch artwork; indent guides are painted by the item delegate."""
        del painter, rect, index

    def set_active_folder(self, folder: Path | None) -> None:
        """Mark which folder owns the Files panel without wiping multi-selection.

        Called after every Files refresh. Must not ``ClearAndSelect``, or
        Cmd/Ctrl+click multi-select is immediately reduced to one row.
        """
        if folder is None:
            self._active_folder_key = None
            self.viewport().update()
            return
        try:
            self._active_folder_key = str(folder.expanduser().resolve())
        except OSError:
            self._active_folder_key = str(folder.expanduser())
        item = self.find_item_for_folder(folder)
        if item is not None and not self._user_selecting:
            # During a click/key press Qt sets the selection after this runs
            # (Cmd/Ctrl+click toggles, Shift extends). Selecting here with the
            # Current flag dropped the earlier rows and let the toggle unselect
            # the clicked one, so Cmd+click emptied the selection.
            index = self.indexFromItem(item)
            if index.isValid():
                self._suppress_current_changed = True
                try:
                    flags = QItemSelectionModel.SelectionFlag.Current
                    if not item.isSelected():
                        flags |= QItemSelectionModel.SelectionFlag.Select
                    self.selectionModel().setCurrentIndex(index, flags)
                finally:
                    QTimer.singleShot(
                        0,
                        lambda: setattr(self, "_suppress_current_changed", False),
                    )
        self.viewport().update()

    def find_item_for_folder(self, folder: Path) -> QTreeWidgetItem | None:
        try:
            key = str(folder.expanduser().resolve())
        except OSError:
            key = str(folder.expanduser())

        def walk(item: QTreeWidgetItem) -> QTreeWidgetItem | None:
            raw = item.data(0, ROLE_PATH)
            if raw:
                try:
                    item_key = str(Path(str(raw)).expanduser().resolve())
                except OSError:
                    item_key = str(raw)
                if item_key == key:
                    return item
            for index in range(item.childCount()):
                found = walk(item.child(index))
                if found is not None:
                    return found
            return None

        for index in range(self.topLevelItemCount()):
            found = walk(self.topLevelItem(index))
            if found is not None:
                return found
        return None

    @staticmethod
    def _is_lazy_placeholder(item: QTreeWidgetItem) -> bool:
        return bool(item.data(0, ROLE_LAZY_PLACEHOLDER))

    @staticmethod
    def _folder_name_for_item(item: QTreeWidgetItem) -> str:
        raw = item.data(0, ROLE_PATH)
        if raw:
            return Path(str(raw)).name.lower()
        label = item.text(0)
        if " (" in label:
            label = label.rsplit(" (", 1)[0]
        return label.lower()

    def apply_name_filter(self, query: str) -> None:
        """Show folder rows whose name matches ``query`` (case-insensitive substring)."""
        needle = query.strip().lower()

        if not needle:

            def reset(item: QTreeWidgetItem) -> None:
                item.setHidden(False)
                for index in range(item.childCount()):
                    child = item.child(index)
                    if self._is_lazy_placeholder(child):
                        child.setHidden(True)
                        continue
                    reset(child)

            for index in range(self.topLevelItemCount()):
                reset(self.topLevelItem(index))
            self.refresh_horizontal_extent()
            return

        def filter_item(item: QTreeWidgetItem) -> bool:
            if self._is_lazy_placeholder(item):
                item.setHidden(True)
                return False
            self_match = needle in self._folder_name_for_item(item)
            child_match = False
            for index in range(item.childCount()):
                if filter_item(item.child(index)):
                    child_match = True
            visible = self_match or child_match
            item.setHidden(not visible)
            if child_match:
                self.expandItem(item)
            return visible

        for index in range(self.topLevelItemCount()):
            filter_item(self.topLevelItem(index))
        self.refresh_horizontal_extent()

    def select_folder_item(self, item: QTreeWidgetItem) -> None:
        """Select exactly one folder row and make it current.

        Skipped while Qt is handling a click / key press in the tree: Qt sets
        the selection itself then, and clearing it here broke Cmd/Ctrl+click.
        """
        index = self.indexFromItem(item)
        if not index.isValid() or self._user_selecting:
            return
        self._suppress_current_changed = True
        try:
            self.selectionModel().select(
                index,
                QItemSelectionModel.SelectionFlag.ClearAndSelect
                | QItemSelectionModel.SelectionFlag.Current,
            )
            self.scrollToItem(item)
        finally:
            QTimer.singleShot(0, lambda: setattr(self, "_suppress_current_changed", False))

    def _emit_folder_for_item(self, item: QTreeWidgetItem | None) -> None:
        if item is None or bool(item.data(0, ROLE_LAZY_PLACEHOLDER)):
            return
        folder = self._folder_path_for_item(item)
        if folder is not None:
            self.folder_selected.emit(folder)

    def _on_folder_current_changed(
        self,
        current: QTreeWidgetItem | None,
        _previous: QTreeWidgetItem | None,
    ) -> None:
        if self._suppress_current_changed:
            return
        self._emit_folder_for_item(current)

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.position().toPoint()
            item = self.itemAt(pos)
            if item is not None and not bool(item.data(0, ROLE_LAZY_PLACEHOLDER)):
                index = self.indexFromItem(item)
                opt = QStyleOptionViewItem()
                opt.rect = self.visualRect(index)
                delegate = self.itemDelegate()
                if delegate is not None:
                    delegate.initStyleOption(opt, index)
                deco = self.style().subElementRect(
                    QStyle.SubElement.SE_ItemViewItemDecoration,
                    opt,
                    self,
                )
                chevron_right = (
                    deco.right()
                    if deco.isValid() and deco.width() > 0
                    else opt.rect.left() + self.iconSize().width()
                )
                # When the row is selected Qt's view-item decoration metrics can
                # shrink, making the chevron area feel "dead". Expand the hit
                # region by a few pixels to keep expand/collapse responsive.
                chevron_right = max(
                    chevron_right,
                    opt.rect.left() + self.iconSize().width() + 8,
                )
                if pos.x() <= chevron_right:
                    has_branch = item.childCount() > 0 or (
                        item.childIndicatorPolicy()
                        == QTreeWidgetItem.ChildIndicatorPolicy.ShowIndicator
                    )
                    if has_branch:
                        item.setExpanded(not item.isExpanded())
        self._user_selecting = True
        try:
            super().mousePressEvent(event)
        finally:
            self._user_selecting = False

    def keyPressEvent(self, event) -> None:  # noqa: N802
        self._user_selecting = True
        try:
            super().keyPressEvent(event)
        finally:
            self._user_selecting = False

    def _folder_path_for_item(self, item: QTreeWidgetItem | None) -> Path | None:
        if item is None:
            return None
        raw = item.data(0, ROLE_PATH)
        if not raw:
            return None
        path = Path(str(raw))
        return path if path.is_dir() else None

    def _selected_folder_items(self) -> list[tuple[Path, QTreeWidgetItem]]:
        pairs: list[tuple[Path, QTreeWidgetItem]] = []
        seen: set[str] = set()
        for item in self.selectedItems():
            folder = self._folder_path_for_item(item)
            if folder is None:
                continue
            key = str(folder.expanduser().resolve())
            if key in seen:
                continue
            seen.add(key)
            pairs.append((folder, item))
        return pairs

    def _is_top_level_item(self, item: QTreeWidgetItem) -> bool:
        return item.parent() is None

    def _show_context_menu(self, position) -> None:
        item = self.itemAt(position)
        if item is not None and item not in self.selectedItems():
            self.clearSelection()
            item.setSelected(True)

        folder_items = self._selected_folder_items() if item is not None else []
        self._show_folder_context_menu(position, folder_items)

    def _show_folder_context_menu(
        self,
        position,
        folder_items: list[tuple[Path, QTreeWidgetItem]],
    ) -> None:
        folders = [folder for folder, _item in folder_items]
        menu = QMenu(self)

        menu.addAction("Add Folder…", self.add_folder_requested.emit)
        if len(folders) == 1:
            suggested = folders[0]
        else:
            window = self.window()
            current = getattr(window, "_current_folder", None)
            suggested = current if isinstance(current, Path) and current.is_dir() else None
        if len(folders) <= 1:
            # The dialog picks the destination, so this is offered even with nothing
            # selected — it falls back to the workspace roots.
            menu.addAction(
                "Create Folder…",
                lambda f=suggested: self.create_folder_requested.emit(f),
            )

        if folders:
            menu.addSeparator()
            if len(folders) == 1:
                folder = folders[0]
                menu.addAction("Rename folder", lambda: self.rename_folder_requested.emit(folder))
                window = self.window()
                if hasattr(window, "_clipboard_can_paste") and window._clipboard_can_paste():
                    menu.addAction(
                        "Paste",
                        lambda f=folder: self.paste_to_folder_requested.emit(f),
                    )
            root_pairs = [
                (folder, tree_item)
                for folder, tree_item in folder_items
                if self._is_top_level_item(tree_item)
            ]
            if root_pairs:
                root_folders = [folder for folder, _tree_item in root_pairs]
                if len(root_folders) == 1:
                    remove_action = menu.addAction(
                        "Remove from list",
                        lambda: self.clear_root_folders_requested.emit(root_folders),
                    )
                else:
                    remove_action = menu.addAction(
                        f"Remove {len(root_folders)} from list",
                        lambda: self.clear_root_folders_requested.emit(root_folders),
                    )
                remove_action.setToolTip(
                    "Hide selected top-level folder(s) from RepeatViewer. Files on disk are not deleted."
                )
            menu.addSeparator()
            color_menu = menu.addMenu("Folder color")
            color_menu.addAction(
                "Choose color…",
                lambda: self._choose_folders_color(folder_items),
            )
            for label, hex_color in FOLDER_COLOR_PRESETS:
                color_menu.addAction(
                    label,
                    lambda checked=False, c=hex_color, pairs=folder_items: self._set_folders_color(pairs, c),
                )
            color_menu.addSeparator()
            color_menu.addAction(
                "Clear color",
                lambda: self._clear_folders_color(folder_items),
            )
            menu.addSeparator()
            if len(folders) == 1:
                menu.addAction(
                    "Delete folder",
                    lambda: self.delete_folders_requested.emit(folders),
                )
            else:
                menu.addAction(
                    f"Delete {len(folders)} folders",
                    lambda: self.delete_folders_requested.emit(folders),
                )
        menu.exec(self.mapToGlobal(position))
    def _set_folders_color(
        self,
        folder_items: list[tuple[Path, QTreeWidgetItem]],
        color: str,
    ) -> None:
        for folder, item in folder_items:
            set_folder_color(folder, color)
            apply_folder_item_color(item, folder)

    def _choose_folders_color(self, folder_items: list[tuple[Path, QTreeWidgetItem]]) -> None:
        first_folder, _item = folder_items[0]
        stored = get_folder_color(first_folder)
        initial = QColor(stored) if stored else QColor(FOLDER_COLOR_PRESETS[0][1])
        if not initial.isValid():
            initial = QColor(FOLDER_COLOR_PRESETS[0][1])

        dialog = QColorDialog(initial, self.window())
        dialog.setWindowTitle("Folder color")
        dialog.setStyleSheet("")
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        picked = dialog.selectedColor()
        if not picked.isValid():
            return
        self._set_folders_color(folder_items, picked.name())

    def _clear_folders_color(self, folder_items: list[tuple[Path, QTreeWidgetItem]]) -> None:
        for folder, item in folder_items:
            set_folder_color(folder, None)
            apply_folder_item_color(item, folder)

    def dragEnterEvent(self, event):  # noqa: N802 (Qt API)
        if event.mimeData().hasUrls():
            event.acceptProposedAction()
            return
        event.ignore()

    def dragMoveEvent(self, event):  # noqa: N802 (Qt API)
        paths = self._mime_local_paths(event.mimeData())
        if not paths:
            self._set_drop_target_item(None)
            event.ignore()
            return

        target_item = self._folder_item_at(event.position().toPoint())
        file_paths = [p for p in paths if p.is_file()]
        folder_paths = [p for p in paths if p.is_dir()]

        # File → folder copy: highlight the hovered folder as the drop target.
        if file_paths:
            if target_item is not None:
                self._set_drop_target_item(target_item)
                event.setDropAction(Qt.DropAction.CopyAction)
                event.accept()
                return
            self._set_drop_target_item(None)
            # Still allow folder-root drops in the same gesture.
            if folder_paths:
                event.acceptProposedAction()
                return
            event.ignore()
            return

        self._set_drop_target_item(None)
        if folder_paths:
            event.acceptProposedAction()
            return
        event.ignore()

    def dragLeaveEvent(self, event):  # noqa: N802 (Qt API)
        self._set_drop_target_item(None)
        super().dragLeaveEvent(event)

    def dropEvent(self, event):  # noqa: N802 (Qt API)
        paths = self._mime_local_paths(event.mimeData())
        self._set_drop_target_item(None)
        if not paths:
            event.ignore()
            return

        target_item = self._folder_item_at(event.position().toPoint())
        target_dir = None
        if target_item is not None:
            raw = target_item.data(0, ROLE_PATH)
            if raw:
                target_dir = Path(str(raw))
                if not target_dir.is_dir():
                    target_dir = None

        folder_paths = [p for p in paths if p.is_dir()]
        file_paths = [p for p in paths if p.is_file()]

        if folder_paths:
            for folder in folder_paths:
                self.folder_dropped.emit(folder)
        if target_dir and file_paths:
            self.files_dropped_to_folder.emit(target_dir, file_paths)
            event.setDropAction(Qt.DropAction.CopyAction)
        event.acceptProposedAction()


@dataclass
class _PendingTeVizBuild:
    """TE visualization run waiting for consensus FASTA BLAST to finish."""

    masking_paths: list[Path]
    consensus_fastas: list[Path]
    gene_paths: list[Path]
    methylation_paths: list[Path]
    map_options: TeGenomeMapOptions


_METADATA_UI_CHUNK = 24
_METADATA_STAT_CHUNK = 48
# Above this, avoid sync metadata/stat and O(n) preview/is_file scans on the GUI thread.
_LARGE_FILE_LIST = 400
_LARGE_JOB_OUTPUT_BATCH = 200


class _MetadataEnrichWorker(QObject):
    """Compute FASTA table metadata off the UI thread."""

    # The generation travels with the batch so this can connect straight to a
    # bound method on MainWindow. Connecting via a lambda gave Qt no receiver
    # QObject, so it used a direct connection and ran the slot on this worker
    # thread -- which then started a GUI-thread QTimer and crashed the app.
    batch_ready = Signal(list, int)  # list[tuple[str, BrowserFileMetadata]], generation
    finished = Signal()

    def __init__(self, paths: list[Path], generation: int) -> None:
        super().__init__()
        self._paths = paths
        self._generation = generation
        self._cancelled = False

    def cancel(self) -> None:
        self._cancelled = True

    @Slot()
    def run(self) -> None:
        from repeatviewer.browser_file_metadata import enrich_fasta_metadata

        batch: list[tuple[str, object]] = []
        batch_size = 48
        for path in self._paths:
            if self._cancelled:
                break
            try:
                meta = enrich_fasta_metadata(path)
                batch.append((str(path), meta))
            except Exception:  # noqa: BLE001
                continue
            if len(batch) >= batch_size:
                self.batch_ready.emit(list(batch), self._generation)
                batch.clear()
        if batch and not self._cancelled:
            self.batch_ready.emit(batch, self._generation)
        self.finished.emit()


class _ToolbarActionButton(QToolButton):
    """Toolbar button with a reliable right-click action.

    We do NOT use ``customContextMenuRequested`` for the right-click: on macOS
    that fires during the press while the button holds the mouse grab, so a slot
    that opens a modal dialog leaves the whole toolbar grabbed and frozen.
    Instead we swallow the right-button press (no grab, no context-menu event)
    and emit ``rightClicked`` on release, after Qt has finished with the event.
    """

    rightClicked = Signal()

    def mousePressEvent(self, event) -> None:  # noqa: N802 (Qt override)
        if event.button() == Qt.MouseButton.RightButton:
            event.accept()  # swallow: don't enter pressed state / grab
            return
        super().mousePressEvent(event)

    def mouseReleaseEvent(self, event) -> None:  # noqa: N802 (Qt override)
        if event.button() == Qt.MouseButton.RightButton:
            event.accept()
            if self.rect().contains(event.pos()):
                self.rightClicked.emit()
            return
        super().mouseReleaseEvent(event)


class MainWindow(QMainWindow):
    """Browser panel with folders, files, and preview."""

    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(APP_DISPLAY_NAME)
        self.resize(1680, 980)
        self.setMinimumSize(900, 700)
        # Cap default size / minimum to the usable screen so Windows title bar
        # (close) and status bar (Stop) are not pushed off-screen.
        from repeatviewer.ui.window_geometry import fit_widget_to_available_screen

        fit_widget_to_available_screen(self)
        self._toolbar_font_size = load_toolbar_font_size()
        self._browser_panel_font_size = load_browser_panel_font_size()
        self._functions_toolbar = None
        self._toolbar_container: QWidget | None = None
        self._toolbar_hbox = None
        self._toolbar_items: list[tuple] = []  # (widget, overflow_handler|None), ordered
        self._overflow_hidden: list[tuple] = []
        self._overflow_btn: QToolButton | None = None
        self._active_toolbar_menu: _ToolbarPopup | None = None
        self._active_toolbar_menu_button: QToolButton | None = None
        self._toolbar_font = QFont(self.font())
        self._browser_panel_font = QFont(self.font())
        self._apply_toolbar_theme()
        # Do not install a Python event filter on QApplication: PySide6 can
        # segfault on macOS when an app-level filter runs during
        # QObject::setProperty (nested getWrapperForQObject). See
        # repeatviewer.views.dialog_polish for the same constraint.

        self._project = SettingsDialog.load_project_from_settings()
        from repeatviewer.app_preferences import load_reopen_last_folder

        if not load_reopen_last_folder():
            # Keep global tool preferences, but start with no folder/genome open.
            self._project.input_dir = None
            self._project.genome_file = None
            self._project.pan_genomes = []
        if self._project.genome_file and self._project.genome_file.is_file():
            add_genome_to_history(self._project.genome_file)
        self._current_folder: Path | None = None
        self._saved_last_folder: Path | None = None
        self._worker_thread: QThread | None = None
        self._worker: FunctionWorker | AlignWorker | WorkflowWorker | None = None
        self._task_progress: TaskProgressPanel | None = None
        self._metadata_enrich_thread: QThread | None = None
        self._metadata_enrich_worker: _MetadataEnrichWorker | None = None
        self._metadata_enrich_generation = 0
        self._pending_metadata_ui: list[tuple[Path, object]] = []
        self._metadata_stat_queue: list[Path] = []
        self._metadata_ui_timer = QTimer(self)
        self._metadata_ui_timer.setSingleShot(True)
        self._metadata_ui_timer.timeout.connect(self._drain_metadata_ui_queue)
        self._job_running = False
        self._preview_open_progress = False
        self._job_finished_callback = None
        self._job_failed_callback = None
        self._last_job_source_files: list[Path] = []
        self._job_toolbar_actions: list[QAction | QToolButton] = []
        self._cached_folder_files: list[Path] = []
        self._highlight_output_keys: set[str] = set()
        self._highlight_source_keys: set[str] = set()
        # Soft-deleted files / folder actions for Cmd/Ctrl+Z restore.
        # Each entry is a dict with ``kind`` in {files, hide_folders, delete_folders}.
        self._browser_undo_stack: list[dict[str, Any]] = []
        self._deleted_files_trash_root = (
            Path(tempfile.gettempdir()) / "RepeatViewerDeletedFiles"
        )
        self._genome_dialog: GenomeDialog | None = None
        self._define_database_dialog: DefineDatabaseDialog | None = None
        self._add_annotation_dialog = None
        self._blast_params_dialog = None
        self._add_blast_database_dialog = None
        self._command_plugin_dialog: CommandPluginDesignerDialog | None = None
        self._plugin_run_dialog: UserCommandRunDialog | None = None
        self._msa_clean_options = load_msa_clean_options()
        self._blast_options = load_blast_options()
        self._fastga_options = load_fastga_options()
        self._tree_options = load_tree_options()
        self._map_options = load_map_options()
        self._extend_options = load_extend_options()
        self._tsd_search_options = load_tsd_search_options()
        self._teaid_options = load_teaid_options()
        self._align_options = load_align_options()
        self._te_genome_map_options = load_te_genome_map_options()
        # Superfamily filtering is controlled on the plot page, not in settings.
        self._te_genome_map_options.superfamilies_only = None
        self._te_divergence_plot_options = load_te_divergence_plot_options()
        self._msa_plot_options = load_msa_plot_options()
        self._consensus_options = load_consensus_options()
        self._flexidot_options = load_flexidot_options()
        self._te_annotation_options = load_te_annotation_options()
        from repeatviewer.backends.auto_annotation_settings import (
            load_orf_prediction_settings,
        )

        self._orf_prediction_settings = load_orf_prediction_settings()
        self._te_confusion_benchmark_options = load_te_confusion_matrix_benchmark_options()
        self._te_library_benchmark_options = load_te_library_benchmark_options()
        self.database_index_runner = DatabaseIndexRunner(self)
        self._plugin_registry = load_plugins()
        self._plugins_menu: QMenu | None = None
        self._others_menu: QMenu | None = None
        self._workflows_menu: QMenu | None = None
        self._toolbar_order = load_toolbar_order(self._plugin_registry)
        self._browser_panel: QWidget | None = None
        self._dock_folders: QDockWidget | None = None
        self._dock_files: QDockWidget | None = None
        self._dock_preview: QDockWidget | None = None
        self._preview_host: QWidget | None = None
        self._preview_file_tabs_host: QWidget | None = None
        self._preview_tabs: PreviewTabWidget | None = None
        self._pending_te_genome_map: TeGenomeMapOptions | None = None
        self._pending_te_genome_map_fastas: list[Path] = []
        self._pending_te_viz: _PendingTeVizBuild | None = None
        # Also holds RepeatQuantifyDialog: any non-modal dialog that implements
        # apply_path_from_panel is driven through this one slot.
        self._te_viz_input_dialog: (
            TeDensityPlotDialog
            | TeDivergencePlotDialog
            | TeAbundancePlotDialog
            | MsaPlotDialog
            | RepeatQuantifyDialog
            | None
        ) = None
        self._genome_chrom_names_cache: list[str] = []
        self._genome_chrom_sizes_cache: dict[str, int] = {}
        self._genome_chrom_cache_path: Path | None = None
        self._genome_map_thread: QThread | None = None
        self._genome_map_worker: GenomeMapWorker | None = None
        self._benchmark_thread: QThread | None = None
        self._benchmark_worker: ConfusionMatrixBenchmarkWorker | TeLibraryBenchmarkWorker | None = None
        self._extract_thread: QThread | None = None
        self._extract_worker: ExtractSequencesWorker | None = None
        self._extract_dialog = None
        self._extract_source: Path | None = None
        self._extract_count_only = False
        self._last_benchmark_source_files: list[Path] = []
        self._job_cancel = JobCancellation()
        self._genome_map_blast_txt: Path | None = None
        self._genome_map_te_paths: list[Path] = []
        self._queued_genome_map_jobs: list[tuple[Path, TeGenomeMapOptions | None]] = []
        self._stop_btn: QPushButton | None = None
        QTimer.singleShot(500, self._warm_genome_chromosome_cache)

        self.setDockOptions(
            QMainWindow.DockOption.AllowNestedDocks
            | QMainWindow.DockOption.AllowTabbedDocks
            | QMainWindow.DockOption.AnimatedDocks
        )
        self.setTabPosition(Qt.DockWidgetArea.LeftDockWidgetArea, QTabWidget.TabPosition.North)

        self._folder_tree = FolderTreeWidget()
        self._folder_tree.setObjectName("FolderTree")
        self._folder_tree.folder_selected.connect(self._on_folder_selected)
        self._folder_tree.files_dropped_to_folder.connect(self._copy_files_to_folder)
        self._folder_tree.folder_dropped.connect(
            lambda folder: self._add_root_folder(folder, unhide=True, warn_spaces=True),
        )

        self._file_list = FileListTableWidget()
        self._file_list.setObjectName("FileList")
        self._file_list.currentItemChanged.connect(self._on_file_current_item_changed)
        self._file_list.itemClicked.connect(self._on_file_item_clicked)
        self._file_list.itemDoubleClicked.connect(self._on_file_double_clicked)
        self._file_list.open_in_preview_requested.connect(self._open_files_in_preview)
        self._file_list.open_in_preview_tab_requested.connect(self._open_files_in_preview_tabs)
        self._file_list.open_in_text_editor_requested.connect(self._open_selected_in_text_editor)
        self._file_list.delete_files_requested.connect(self._delete_selected_files)
        self._file_list.cut_files_requested.connect(self._cut_files_to_clipboard)
        self._file_list.copy_files_requested.connect(self._copy_files_to_clipboard)
        self._file_list.paste_files_requested.connect(self._paste_files_shortcut)
        self._file_list.duplicate_files_requested.connect(self._duplicate_selected_files)
        self._file_list.rename_file_requested.connect(self._rename_file)
        self._file_list.extract_sequences_requested.connect(
            self._extract_sequences_from_list
        )
        self._file_list.extract_hit_subjects_requested.connect(
            self._extract_sequences_from_hit_table
        )
        self._file_list.group_sequences_requested.connect(self._group_sequences)
        self._file_list.new_folder_with_files_requested.connect(
            self._create_folder_with_files
        )
        self._file_list.separate_fasta_requested.connect(self._separate_fasta_file)
        self._file_list.separate_fasta_subfolders_requested.connect(
            self._separate_fasta_to_subfolders_file,
        )
        self._file_list.blast_genome_map_requested.connect(
            lambda p: self._open_te_density_plot_dialog(input_file=p),
        )
        self._file_list.blast_divergence_requested.connect(
            lambda p: self._open_te_divergence_plot_dialog(input_file=p),
        )
        self._file_list.add_blast_database_requested.connect(self._add_fasta_as_blast_database)
        self._file_list.add_annotation_requested.connect(self._add_annotation_for_fasta)
        self._file_list.attach_cyto_image_requested.connect(self._attach_cyto_images)
        self._file_list.cyto_images_requested.connect(self._open_cyto_images)
        from repeatviewer.widgets.cyto_images_ui import cyto_notifier

        cyto_notifier().changed.connect(self._on_cyto_images_changed)
        self._folder_tree.rename_folder_requested.connect(self._rename_folder)
        self._folder_tree.delete_folders_requested.connect(self._delete_folders)
        self._folder_tree.clear_root_folders_requested.connect(self._clear_root_folders)
        self._folder_tree.paste_to_folder_requested.connect(self._paste_files_to_folder)
        self._folder_tree.add_folder_requested.connect(self._pick_and_add_folder)
        self._folder_tree.create_folder_requested.connect(self._create_folder)
        self._folder_tree.itemExpanded.connect(self._on_folder_tree_expanded)
        self._folder_tree.itemCollapsed.connect(self._on_folder_tree_collapsed)
        self._apply_browser_panel_font()
        self._fs_watcher = QFileSystemWatcher(self)
        self._fs_watcher.directoryChanged.connect(self._on_watched_directory_changed)
        self._fs_refresh_timer = QTimer(self)
        self._fs_refresh_timer.setSingleShot(True)
        self._fs_refresh_timer.setInterval(400)
        self._fs_refresh_timer.timeout.connect(self._deferred_fs_folder_refresh)
        self._pending_fs_refresh_dir: str | None = None
        self._watcher_refresh_suppressed_until = 0.0
        self._folder_persist_timer = QTimer(self)
        self._folder_persist_timer.setSingleShot(True)
        self._folder_persist_timer.setInterval(250)
        self._folder_persist_timer.timeout.connect(self._persist_workspace)
        self._folder_tree_restoring = False
        self._hidden_root_folders = load_hidden_root_folders()
        self._file_clipboard: _FileClipboard | None = None

        self._build_dock_panels()

        status = QStatusBar()
        self.setStatusBar(status)
        self._status_linked_file: Path | None = None
        self._status_file_link = QLabel()
        self._status_file_link.setTextFormat(Qt.TextFormat.RichText)
        self._status_file_link.setOpenExternalLinks(False)
        self._status_file_link.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextBrowserInteraction
        )
        self._status_file_link.linkActivated.connect(self._on_status_file_link_activated)
        self._status_file_link.hide()
        # Left-side temporary widget: replaces status text when a file is selected.
        status.addWidget(self._status_file_link)
        self._status_link_restore_timer = QTimer(self)
        self._status_link_restore_timer.setSingleShot(True)
        self._status_link_restore_timer.timeout.connect(self._restore_status_file_link)
        self._task_progress = TaskProgressPanel()
        status.addPermanentWidget(self._task_progress)
        self._stop_btn = QPushButton("Stop")
        self._stop_btn.setToolTip("Stop the currently running job")
        self._stop_btn.setMaximumHeight(22)
        self._stop_btn.setEnabled(False)
        self._stop_btn.clicked.connect(self._stop_current_job)
        status.addPermanentWidget(self._stop_btn)
        self._build_menu()
        self._install_shortcuts()
        app = QApplication.instance()
        if app is not None:
            app.applicationStateChanged.connect(self._ensure_floating_docks_stay_visible)
        QTimer.singleShot(0, self._deferred_window_startup)

    def _deferred_window_startup(self) -> None:
        """Restore layout and title-bar tabs after the main window is shown."""
        self._restore_window_layout()
        self._mount_preview_title_file_tabs()
        self._complete_deferred_startup()

    def _complete_deferred_startup(self) -> None:
        """Load workspace folders after the window is shown (avoids blocking startup)."""
        self._load_initial_roots()
        self._apply_startup_folder()
        self._normalize_folder_tree_selection()
        self._update_directory_watchers()
        self._update_status()
        QTimer.singleShot(0, self._finalize_browser_layout)

    def _finalize_browser_layout(self) -> None:
        """Re-check dock tiling after the first paint (restoreState can lag)."""
        settings = QSettings()
        settings.beginGroup("mainwindow")
        has_saved_layout = bool(settings.value("toolbar_state_v5"))
        settings.endGroup()
        if self._browser_layout_needs_repair() or not has_saved_layout:
            self._apply_default_browser_layout(force=True)
        else:
            self._apply_default_browser_layout()
        from repeatviewer.ui.window_geometry import (
            ensure_floating_docks_on_screen,
            fit_widget_to_available_screen,
        )

        ensure_floating_docks_on_screen(self)
        fit_widget_to_available_screen(self)
        self._repair_browser_chrome()

    def _normalize_folder_tree_selection(self) -> None:
        """Ensure a current folder exists after restore; do not clear multi-select."""
        if self._folder_tree.currentItem() is not None:
            return
        if self._folder_tree.topLevelItemCount() <= 0:
            return
        item = self._folder_tree.topLevelItem(0)
        if item is None:
            return
        self._folder_tree.select_folder_item(item)

    def _update_directory_watchers(self) -> None:
        watched: set[str] = set()
        for root in self._collect_root_folders():
            if root.is_dir():
                watched.add(str(root.expanduser().resolve()))
        if self._current_folder and self._current_folder.is_dir():
            watched.add(str(self._current_folder.expanduser().resolve()))

        current = set(self._fs_watcher.directories())
        to_remove = sorted(current - watched)
        if to_remove:
            self._fs_watcher.removePaths(to_remove)
        to_add = sorted(watched - current)
        if to_add:
            self._fs_watcher.addPaths(to_add)

    def _on_watched_directory_changed(self, directory: str) -> None:
        if time.monotonic() < self._watcher_refresh_suppressed_until:
            return
        # Tool jobs write many files into the open folder. Rebuilding a large
        # file list on every change slows the job (disk + CPU) and freezes browsing.
        # Refresh once when the job finishes instead.
        if self._any_background_job_running():
            return
        self._pending_fs_refresh_dir = directory
        self._fs_refresh_timer.start()

    def _deferred_fs_folder_refresh(self) -> None:
        directory = self._pending_fs_refresh_dir
        self._pending_fs_refresh_dir = None
        if directory is None:
            return
        if self._any_background_job_running():
            return
        changed = Path(directory).expanduser().resolve()
        refreshed = False

        if changed.is_dir():
            refreshed = (
                self._refresh_subfolders_under(
                    changed,
                    keep_folder=self._current_folder,
                )
                is not None
            )

        if not refreshed:
            self._refresh_tree()

        if self._current_folder and self._current_folder.is_dir():
            self._populate_file_list(self._current_folder, open_preview=False)
        self._update_directory_watchers()

    def _build_file_columns_menu(self) -> None:
        from repeatviewer.widgets.file_list_columns import (
            TOGGLEABLE_COLUMNS,
            load_visible_file_columns,
        )

        self._file_columns_menu.clear()
        self._file_column_actions.clear()
        visible = load_visible_file_columns()
        for column, label in TOGGLEABLE_COLUMNS.items():
            action = self._file_columns_menu.addAction(label)
            action.setCheckable(True)
            action.setChecked(column in visible)
            action.toggled.connect(
                lambda checked, col=column: self._on_file_column_toggled(col, checked)
            )
            self._file_column_actions[column] = action

    def _on_file_column_toggled(self, column: int, checked: bool) -> None:
        from repeatviewer.widgets.file_list_columns import (
            load_visible_file_columns,
            metadata_columns_visible,
            save_visible_file_columns,
        )

        visible = load_visible_file_columns()
        if checked:
            visible.add(column)
        else:
            visible.discard(column)
        save_visible_file_columns(visible)
        self._file_list.apply_visible_columns(visible)
        if checked and metadata_columns_visible(visible):
            QTimer.singleShot(0, self._fill_visible_file_metadata)

    def _fill_visible_file_metadata(self) -> None:
        """Populate metadata columns after they are turned on or for stale rows."""
        from repeatviewer.browser_file_metadata import is_fasta_path

        paths = self._file_list.listed_paths()
        if not paths:
            return
        self._metadata_stat_queue = [path for path in paths if not is_fasta_path(path)]
        fasta_paths = [path for path in paths if is_fasta_path(path)]
        if fasta_paths:
            self._start_metadata_enrichment(fasta_paths)
        if self._metadata_stat_queue:
            self._schedule_metadata_ui_drain()

    def _schedule_metadata_ui_drain(self) -> None:
        if not self._metadata_ui_timer.isActive():
            self._metadata_ui_timer.start(0)

    def _collect_metadata_stat_chunk(self) -> list[tuple[Path, object]]:
        from repeatviewer.browser_file_metadata import quick_metadata_for_file

        chunk = self._metadata_stat_queue[:_METADATA_STAT_CHUNK]
        del self._metadata_stat_queue[:_METADATA_STAT_CHUNK]
        return [(path, quick_metadata_for_file(path)) for path in chunk]

    def _sync_folder_tree_label_for_folder(self, folder: Path) -> None:
        item = self._find_tree_item_for_path(folder.expanduser().resolve())
        if item is not None:
            self._set_folder_tree_item_label(item, folder)

    def _incremental_add_files_to_list(self, paths: list[Path]) -> None:
        """Append new rows without rebuilding the whole table."""
        if not paths:
            return
        query = self._file_filter.text().strip().lower()
        visible_paths: list[Path] = []
        for path in paths:
            if query and query not in path.name.lower():
                continue
            visible_paths.append(path)
        if not visible_paths:
            return

        was_sorting = self._file_list.isSortingEnabled()
        self._file_list.setSortingEnabled(False)
        self._file_list.setUpdatesEnabled(False)
        try:
            for path in visible_paths:
                row_item = self._file_list.add_file_row(
                    path,
                    fill_metadata_columns=False,
                )
                key = str(path.expanduser().resolve())
                if key in self._highlight_output_keys:
                    self._style_new_output_item(row_item)
                elif key in self._highlight_source_keys:
                    self._style_source_file_item(row_item)
        finally:
            self._file_list.setUpdatesEnabled(True)
            self._file_list.setSortingEnabled(was_sorting)
        self._file_list.apply_default_sort()

    def _restore_window_layout(self) -> None:
        import os

        self._purge_legacy_docks()
        settings = QSettings()
        settings.beginGroup("mainwindow")
        if os.environ.get("REPEATVIEWER_RESET_LAYOUT", "0").lower() in {"1", "true", "yes"}:
            settings.remove("geometry")
            settings.remove("toolbar_state_v5")
        # v5 = all browser docks share one area (no Left/Right central gutter).
        toolbar_state = settings.value("toolbar_state_v5")
        window_geometry = settings.value("geometry")
        settings.endGroup()
        from repeatviewer.app_preferences import load_restore_geometry

        if window_geometry and load_restore_geometry():
            try:
                if not self.restoreGeometry(window_geometry):
                    settings.beginGroup("mainwindow")
                    settings.remove("geometry")
                    settings.endGroup()
            except Exception:
                settings.beginGroup("mainwindow")
                settings.remove("geometry")
                settings.endGroup()
        if toolbar_state:
            try:
                if not self.restoreState(toolbar_state):
                    settings.beginGroup("mainwindow")
                    settings.remove("toolbar_state_v5")
                    settings.endGroup()
            except Exception:
                settings.beginGroup("mainwindow")
                settings.remove("toolbar_state_v5")
                settings.endGroup()
        from repeatviewer.ui.window_geometry import (
            ensure_floating_docks_on_screen,
            fit_widget_to_available_screen,
        )

        self._purge_legacy_docks()
        if self._browser_layout_needs_repair():
            self._apply_default_browser_layout(force=True)
        else:
            self._apply_default_browser_layout()
        self._ensure_preview_visible()
        fit_widget_to_available_screen(self)
        ensure_floating_docks_on_screen(self)
        QTimer.singleShot(0, lambda: fit_widget_to_available_screen(self))
        QTimer.singleShot(0, lambda: ensure_floating_docks_on_screen(self))
        # restoreState can hide the status bar; Stop lives there.
        status = self.statusBar()
        if status is not None:
            status.show()
        self._rebuild_top_toolbar()
        self._repair_browser_chrome()

    def _browser_layout_needs_repair(self) -> bool:
        """True when dock geometry leaves a visible central gutter or wrong areas."""
        if self._central_gap_width() > 12:
            return True
        for dock in (self._dock_folders, self._dock_files, self._dock_preview):
            if dock is None:
                continue
            if self.dockWidgetArea(dock) != Qt.DockWidgetArea.LeftDockWidgetArea:
                return True
        return False

    def _central_gap_width(self) -> int:
        central = self.centralWidget()
        if central is None:
            return 0
        return max(0, central.width())

    def _repair_browser_chrome(self) -> None:
        """Fix file-table header glitches and toolbar overlap after modal dialogs."""
        if hasattr(self, "_file_list"):
            self._file_list.repair_chrome()
        if self._functions_toolbar is not None:
            self._compact_functions_toolbar(self._functions_toolbar)

    def changeEvent(self, event: QEvent) -> None:  # noqa: N802 (Qt API)
        super().changeEvent(event)
        if (
            event.type() == QEvent.Type.ActivationChange
            and self.isActiveWindow()
            and hasattr(self, "_file_list")
        ):
            QTimer.singleShot(0, self._repair_browser_chrome)

    def _purge_legacy_docks(self) -> None:
        """Drop stale dock widgets from older layouts; keep the current three panels."""
        keep = {
            id(dock)
            for dock in (self._dock_folders, self._dock_files, self._dock_preview)
            if dock is not None
        }
        for dock in list(self.findChildren(QDockWidget)):
            if id(dock) in keep:
                continue
            name = dock.objectName().lower()
            if name in {"dock_browser", "dock_preview", "dock_folders", "dock_files"}:
                self.removeDockWidget(dock)
                dock.hide()
                dock.setWidget(QWidget())
                dock.deleteLater()

    def closeEvent(self, event) -> None:  # noqa: N802 (Qt API)
        self._stop_metadata_enrichment()
        if self._any_background_job_running():
            self._shutdown_background_threads()
        if self._preview_tabs is not None:
            self._preview_tabs.teardown_all_html_previews()
        self._persist_workspace()
        # Detached docks can outlive the main window (independent tool windows).
        # Redock/close them before saving layout so they vanish on quit.
        self._close_detached_panels()
        from repeatviewer.ui.window_geometry import ensure_floating_docks_on_screen

        ensure_floating_docks_on_screen(self)
        settings = QSettings()
        settings.beginGroup("mainwindow")
        settings.setValue("toolbar_state_v5", self.saveState())
        settings.setValue("geometry", self.saveGeometry())
        settings.endGroup()
        super().closeEvent(event)

    def _close_detached_panels(self) -> None:
        """Close every floating browser dock and standalone preview window."""
        for dock in self.findChildren(QDockWidget):
            if not dock.isFloating():
                continue
            dock.setAttribute(Qt.WidgetAttribute.WA_MacAlwaysShowToolWindow, False)
            dock.setFloating(False)
            if dock.isFloating():
                dock.hide()
                dock.close()
        from repeatviewer.widgets.preview_duplicate_window import PreviewDuplicateWindow

        app = QApplication.instance()
        if app is None:
            return
        for widget in list(app.topLevelWidgets()):
            if isinstance(widget, PreviewDuplicateWindow):
                widget.close()
            elif (
                isinstance(widget, QDockWidget)
                and widget.isFloating()
                and widget is not self
            ):
                widget.hide()
                widget.close()

    def _shutdown_background_threads(self, *, wait_ms: int = 30_000) -> None:
        """Stop and join background QThreads before teardown or window close."""
        self._job_cancel.request()
        self._queued_genome_map_jobs.clear()

        threads: list[QThread] = []
        for attr in (
            "_worker_thread",
            "_genome_map_thread",
            "_benchmark_thread",
            "_extract_thread",
        ):
            thread = getattr(self, attr, None)
            if isinstance(thread, QThread):
                threads.append(thread)

        for thread in threads:
            if thread.isRunning():
                thread.requestInterruption()
                thread.quit()

        for thread in threads:
            if not thread.isRunning():
                continue
            if thread.wait(wait_ms):
                continue
            logger.warning("Background thread did not stop in time; terminating.")
            thread.terminate()
            thread.wait(3000)

        self._worker = None
        self._worker_thread = None
        self._genome_map_worker = None
        self._genome_map_thread = None
        self._benchmark_worker = None
        self._benchmark_thread = None
        self._job_running = False
        if self._task_progress is not None:
            self._task_progress.reset()
        if self._stop_btn is not None:
            self._stop_btn.setEnabled(False)

    @staticmethod
    def _job_display_title(job_title: str) -> str:
        if job_title.startswith("plugin:"):
            plugin = get_registry().get(job_title[len("plugin:") :])
            if plugin is not None:
                return plugin.label
        if job_title.startswith("workflow:"):
            from repeatviewer.workflows import get_workflow_spec

            spec = get_workflow_spec(job_title[len("workflow:") :])
            if spec is not None:
                return f"Workflow {spec.name}"
            return "Workflow"
        key = job_title.strip().lower().replace(" ", "_")
        return _JOB_DISPLAY_TITLES.get(key, job_title)

    @staticmethod
    def _parse_progress_filename(label: str) -> str:
        if ": " in label:
            part = label.split(": ", 1)[1]
        elif label.lower().startswith("aligning "):
            part = label[9:]
        else:
            return ""
        part = part.removesuffix("…").strip()
        if not part or part.lower() in {"finished", "alignment finished."}:
            return ""
        return part

    def _running_job_title(self, job_title: str) -> str:
        display = self._job_display_title(job_title)
        if display.lower().startswith("running "):
            return display if display.endswith("…") else f"{display}…"
        return f"Running {display}…"

    def _show_job_progress(self, job_title: str, total_steps: int) -> None:
        self._job_running = True
        self._preview_open_progress = False
        self._job_progress_epoch = getattr(self, "_job_progress_epoch", 0) + 1
        self._active_job_log_title = job_title
        # Worker steps only; +1 reserved so the bar fills after results open.
        self._job_work_steps = max(1, int(total_steps))
        operation_logger.info(
            "JOB START name=%s total_steps=%s",
            self._job_display_title(job_title),
            self._job_work_steps,
        )
        self.statusBar().clearMessage()
        if self._task_progress is not None:
            self._task_progress.start(
                self._running_job_title(job_title),
                self._job_work_steps + 1,
            )
        if self._stop_btn is not None:
            self._stop_btn.setEnabled(True)

    def _hide_job_progress(self) -> None:
        self._job_running = False
        self._preview_open_progress = False
        if self._task_progress is not None:
            self._task_progress.reset()
        if self._stop_btn is not None and self._genome_map_thread is None and self._benchmark_thread is None:
            self._stop_btn.setEnabled(False)

    def _finish_background_job(self) -> None:
        """Worker thread done — keep progress visible until results are opened."""
        self._set_job_actions_enabled(True)
        if self._stop_btn is not None:
            self._stop_btn.setEnabled(False)
        self._job_open_epoch = getattr(self, "_job_progress_epoch", 0)
        if self._task_progress is not None and self._job_running:
            work = getattr(self, "_job_work_steps", 1)
            self._task_progress.update(work, "Opening results…")

    def _after_job_ui_ready(self) -> None:
        """Fill the progress bar only after refresh/preview has painted."""
        # A follow-on job (e.g. genome map after BLAST) may have taken the bar.
        if getattr(self, "_job_open_epoch", None) != getattr(
            self, "_job_progress_epoch", None
        ):
            return
        if self._task_progress is not None and self._job_running:
            total = max(1, int(self._task_progress._total))
            self._task_progress.update(total, "Done")
        self._hide_job_progress()
        # Large folders skipped sync metadata during the finish refresh.
        from repeatviewer.widgets.file_list_columns import metadata_columns_visible

        if (
            metadata_columns_visible()
            and self._file_list is not None
            and self._file_list.rowCount() > _LARGE_FILE_LIST
            and not self._any_background_job_running()
        ):
            QTimer.singleShot(0, self._fill_visible_file_metadata)

    def _schedule_after_job_ui_ready(self) -> None:
        """Defer fill/hide so preview paint is not interrupted by processEvents."""
        QTimer.singleShot(0, self._after_job_ui_ready)

    @Slot(dict)
    @Slot(list)
    @Slot(object)
    def _on_background_job_finished(self, result: object = None) -> None:
        """GUI-thread slot: never connect bare lambdas to worker.finished."""
        title = getattr(self, "_active_job_log_title", "job")
        self._log_job_finished(title, result)
        self._finish_background_job()
        callback = self._job_finished_callback
        self._job_finished_callback = None
        try:
            if callback is not None:
                callback(result)
        finally:
            self._schedule_after_job_ui_ready()

    @Slot(str)
    def _on_background_job_failed(self, error: str) -> None:
        """GUI-thread slot: never connect bare lambdas to worker.failed."""
        title = getattr(self, "_active_job_log_title", "job")
        self._log_job_failed(title, error)
        callback = self._job_failed_callback
        self._job_failed_callback = None
        try:
            if callback is not None:
                callback(error)
        finally:
            self._set_job_actions_enabled(True)
            if self._stop_btn is not None:
                self._stop_btn.setEnabled(False)
            self._hide_job_progress()

    def begin_preview_open_progress(self, filename: str, *, total: int = 0) -> bool:
        """Show status-bar progress while opening a large preview file.

        Does not mark a cancellable job (Stop stays disabled). Skipped when a
        real background job already owns the progress panel.

        Returns True when the bar is actually shown.
        """
        if self._job_running or self._task_progress is None:
            return False
        from repeatviewer.widgets.task_progress import truncate_filename

        self._preview_open_progress = True
        self.statusBar().clearMessage()
        short = truncate_filename(str(filename or ""), 28)
        title = f"Opening {short}…" if short else "Opening…"
        self._task_progress.start(title, int(total))
        # Preview open is sync on the GUI thread — pump paint so the bar is
        # visible during the read (safe here; not used on the job-finish path).
        from PySide6.QtCore import QEventLoop
        from PySide6.QtWidgets import QApplication

        QApplication.processEvents(
            QEventLoop.ProcessEventsFlag.ExcludeUserInputEvents
        )
        return True

    def update_preview_open_progress(self, completed: int, detail: str = "") -> None:
        if not self._preview_open_progress or self._task_progress is None:
            return
        self._task_progress.update(int(completed), detail)
        from PySide6.QtCore import QEventLoop
        from PySide6.QtWidgets import QApplication

        QApplication.processEvents(
            QEventLoop.ProcessEventsFlag.ExcludeUserInputEvents
        )

    def end_preview_open_progress(self) -> None:
        if not self._preview_open_progress:
            return
        self._preview_open_progress = False
        if self._job_running:
            return
        if self._task_progress is not None:
            self._task_progress.reset()

    def _set_status_file_link(self, path: Path | None) -> None:
        """Replace the left status text with a clickable filename (Finder link)."""
        self._status_link_restore_timer.stop()
        if path is None:
            self._status_linked_file = None
            self._status_file_link.hide()
            self._status_file_link.clear()
            self._update_status()
            return
        try:
            resolved = path.expanduser().resolve()
        except OSError:
            self._status_linked_file = None
            self._status_file_link.hide()
            self._status_file_link.clear()
            self._update_status()
            return
        self._status_linked_file = resolved
        safe_name = resolved.name.replace("&", "&amp;").replace("<", "&lt;")
        self._status_file_link.setText(f'<a href="reveal">{safe_name}</a>')
        self._status_file_link.setToolTip(f"Show in file manager\n{resolved}")
        # Clear status text first — showing both causes the stacked/overlap glitch.
        self.statusBar().clearMessage()
        self._status_file_link.show()

    def _restore_status_file_link(self) -> None:
        """After a timed status message, show the file link again if still set."""
        if self._status_linked_file is None:
            return
        if self.statusBar().currentMessage():
            return
        self._status_file_link.show()

    def _on_status_file_link_activated(self, _href: str) -> None:
        if self._status_linked_file is not None:
            reveal_in_file_manager(self._status_linked_file)

    def _set_status_message(
        self,
        message: str,
        timeout: int = 0,
        *,
        override_file_link: bool = False,
    ) -> None:
        """Show status text on the left.

        While a file link is active, routine messages are skipped so the link
        stays put. Job results pass ``override_file_link=True`` to temporarily
        replace the link; it returns after ``timeout``.
        """
        if self._job_running:
            return
        if self._status_linked_file is not None and not override_file_link:
            return
        self._status_file_link.hide()
        self.statusBar().showMessage(message, timeout)
        if self._status_linked_file is not None and timeout > 0:
            self._status_link_restore_timer.start(timeout + 50)
        elif self._status_linked_file is not None and timeout == 0:
            # Persistent override message — link stays hidden until cleared.
            self._status_link_restore_timer.stop()


    def _make_browser_dock(
        self,
        title: str,
        object_name: str,
        widget: QWidget,
        *,
        dock_cls: type[QDockWidget] = QDockWidget,
    ) -> QDockWidget:
        dock = dock_cls(title, self)
        dock.setObjectName(object_name)
        dock.setWidget(widget)
        dock.setAllowedAreas(Qt.DockWidgetArea.AllDockWidgetAreas)
        dock.setFeatures(
            QDockWidget.DockWidgetFeature.DockWidgetClosable
            | QDockWidget.DockWidgetFeature.DockWidgetMovable
            | QDockWidget.DockWidgetFeature.DockWidgetFloatable
        )
        # Floating docks default to Qt::Tool and hide when another app covers
        # the main window on macOS. Promote them to independent windows.
        dock.setAttribute(Qt.WidgetAttribute.WA_MacAlwaysShowToolWindow, True)
        dock.topLevelChanged.connect(self._on_browser_dock_floating_changed)
        return dock

    def _on_browser_dock_floating_changed(self, floating: bool) -> None:
        dock = self.sender()
        if not isinstance(dock, QDockWidget):
            return
        if floating:
            # Wait until Qt finishes its own float/show, then promote once.
            QTimer.singleShot(0, lambda d=dock: self._make_floating_dock_independent(d))
        else:
            dock.setProperty("_rv_float_independent", False)
            dock.setProperty("_rv_float_promoting", False)

    def _make_floating_dock_independent(self, dock: QDockWidget) -> None:
        """Turn a floated dock into a normal top-level window (not a Tool panel)."""
        if not dock.isFloating():
            return
        if bool(dock.property("_rv_float_independent")):
            self._detach_floating_dock_transient(dock)
            return
        if bool(dock.property("_rv_float_promoting")):
            return

        geo = QRect(dock.geometry())
        dock.setProperty("_rv_float_promoting", True)
        dock.setProperty("_rv_float_independent", True)
        dock.setAttribute(Qt.WidgetAttribute.WA_MacAlwaysShowToolWindow, True)

        # Hide before changing flags. Do not call QWindow.setFlags() in parallel —
        # that races with show() and logs ``Already setting window visible!``.
        dock.hide()
        dock.setWindowFlags(
            Qt.WindowType.Window | Qt.WindowType.FramelessWindowHint
        )
        dock.setGeometry(geo)

        def _show_independent() -> None:
            dock.setProperty("_rv_float_promoting", False)
            if not dock.isFloating():
                dock.setProperty("_rv_float_independent", False)
                return
            if not dock.isVisible():
                dock.show()
            self._detach_floating_dock_transient(dock)

        QTimer.singleShot(0, _show_independent)

    @staticmethod
    def _detach_floating_dock_transient(dock: QDockWidget) -> None:
        """Drop the transient-parent link so floated docks stay above other apps."""
        if not dock.isFloating():
            return
        handle = dock.windowHandle()
        if handle is None:
            win = dock.window()
            handle = win.windowHandle() if win is not None else None
        if handle is None:
            return
        # Only clear the transient parent. Changing QWindow flags here re-enters
        # the Cocoa visibility path and triggers ``Already setting window visible!``.
        handle.setTransientParent(None)

    def _ensure_floating_docks_stay_visible(self, state=None) -> None:
        """If macOS hid floated panels when the app deactivated, show them again."""
        for dock in (self._dock_folders, self._dock_files, self._dock_preview):
            if dock is None or not dock.isFloating():
                continue
            if bool(dock.property("_rv_float_promoting")):
                continue
            if not dock.isVisible():
                dock.setVisible(True)
            self._detach_floating_dock_transient(dock)

    def _build_dock_panels(self) -> None:
        # --- Folders (nav + filter live on the dock title bar) ---
        folder_side = QWidget()
        folder_side.setObjectName("FoldersPanel")
        folder_side_layout = QVBoxLayout(folder_side)
        folder_side_layout.setContentsMargins(4, 6, 4, 4)
        folder_side_layout.setSpacing(0)
        folder_side_layout.addWidget(self._folder_tree, stretch=1)
        self._folder_tree.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )

        self._folder_nav = FolderNavigationBar()
        self._folder_nav.setObjectName("FoldersDockNavHost")
        self._folder_nav.set_roots_provider(self._workspace_root_folders)
        self._folder_nav.navigate_requested.connect(self._on_folder_nav_requested)

        folders_trailing = QWidget()
        folders_trailing.setObjectName("FoldersDockTrailing")
        folders_trailing.setFixedHeight(18)
        folders_trailing_layout = QHBoxLayout(folders_trailing)
        folders_trailing_layout.setContentsMargins(0, 0, 0, 0)
        folders_trailing_layout.setSpacing(6)
        folders_trailing_layout.setAlignment(Qt.AlignmentFlag.AlignVCenter)
        self._folder_add_btn = QToolButton()
        self._folder_add_btn.setObjectName("FolderAddButton")
        self._folder_add_btn.setText("+")
        self._folder_add_btn.setToolTip("Add folder…")
        self._folder_add_btn.setAutoRaise(True)
        self._folder_add_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._folder_add_btn.setFixedSize(18, 18)
        self._folder_add_btn.clicked.connect(self._pick_and_add_folder)
        folders_trailing_layout.addWidget(
            self._folder_add_btn,
            0,
            Qt.AlignmentFlag.AlignVCenter,
        )
        folders_trailing_layout.addWidget(self._folder_nav, 0, Qt.AlignmentFlag.AlignVCenter)
        self._folder_filter = QLineEdit()
        self._folder_filter.setObjectName("FolderFilter")
        self._folder_filter.setPlaceholderText("Filter by name…")
        self._folder_filter.setClearButtonEnabled(True)
        self._folder_filter.setFrame(False)
        self._folder_filter.setAttribute(Qt.WidgetAttribute.WA_MacShowFocusRect, False)
        self._folder_filter.setFixedHeight(18)
        self._folder_filter.setMaximumHeight(18)
        self._folder_filter.setMinimumWidth(140)
        self._folder_filter.setMaximumWidth(280)
        self._folder_filter.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        self._folder_filter.textChanged.connect(self._apply_folder_filter)
        folders_trailing_layout.addWidget(
            self._folder_filter,
            0,
            Qt.AlignmentFlag.AlignVCenter,
        )

        # --- Files (list only; filter lives on the dock title bar) ---
        files_panel = QWidget()
        files_panel.setObjectName("FilesPanel")
        files_layout = QVBoxLayout(files_panel)
        files_layout.setContentsMargins(4, 6, 4, 4)
        files_layout.setSpacing(0)

        files_trailing = QWidget()
        files_trailing.setObjectName("FilesDockTrailing")
        files_trailing.setFixedHeight(18)
        files_trailing_layout = QHBoxLayout(files_trailing)
        files_trailing_layout.setContentsMargins(0, 0, 0, 0)
        files_trailing_layout.setSpacing(6)
        files_trailing_layout.setAlignment(Qt.AlignmentFlag.AlignVCenter)
        self._file_filter = QLineEdit()
        self._file_filter.setObjectName("FileFilter")
        self._file_filter.setPlaceholderText("Filter by name…")
        self._file_filter.setClearButtonEnabled(True)
        self._file_filter.setFrame(False)
        self._file_filter.setAttribute(Qt.WidgetAttribute.WA_MacShowFocusRect, False)
        self._file_filter.setFixedHeight(18)
        self._file_filter.setMaximumHeight(18)
        self._file_filter.setMinimumWidth(140)
        self._file_filter.setMaximumWidth(280)
        self._file_filter.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        self._file_filter.textChanged.connect(self._apply_file_filter)
        files_trailing_layout.addWidget(
            self._file_filter,
            0,
            Qt.AlignmentFlag.AlignVCenter,
        )

        self._file_columns_btn = QToolButton()
        self._file_columns_btn.setObjectName("FileColumnsButton")
        self._file_columns_btn.setText("Columns")
        self._file_columns_btn.setToolTip("Choose which file-table columns are visible")
        self._file_columns_btn.setFixedHeight(18)
        self._file_columns_btn.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self._file_columns_menu = QMenu(self._file_columns_btn)
        self._file_column_actions: dict[int, QAction] = {}
        self._build_file_columns_menu()
        self._file_columns_btn.setMenu(self._file_columns_menu)
        files_trailing_layout.addWidget(
            self._file_columns_btn,
            0,
            Qt.AlignmentFlag.AlignVCenter,
        )

        self._file_list.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        files_layout.addWidget(self._file_list, stretch=1)

        # --- Preview ---
        from repeatviewer.widgets.preview_controls import (
            AnnotationZoomBar,
            PreviewHeaderBar,
            PreviewZoomBar,
        )

        self._preview_host = QWidget()
        self._preview_host.setObjectName("PreviewHost")
        self._preview_host.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        preview_layout = QVBoxLayout(self._preview_host)
        preview_layout.setContentsMargins(0, 0, 0, 0)
        preview_layout.setSpacing(0)
        self._preview_zoom_bar = PreviewZoomBar()
        self._annotation_zoom_bar = AnnotationZoomBar()
        self._annotation_zoom_bar.hide()
        self._preview_header = PreviewHeaderBar(
            self._preview_zoom_bar,
            self._preview_host,
            annotation_zoom_bar=self._annotation_zoom_bar,
        )
        preview_layout.addWidget(self._preview_header)
        self._preview_tabs = PreviewTabWidget(
            zoom_bar=self._preview_zoom_bar,
            header_bar=self._preview_header,
            annotation_zoom_bar=self._annotation_zoom_bar,
        )
        self._preview_tabs.setObjectName("PreviewTabs")
        self._preview_tabs.annotationsSaved.connect(self._on_preview_annotations_saved)
        self._preview_tabs.sequenceExtracted.connect(self._on_preview_sequence_extracted)
        self._preview_tabs.fileWritten.connect(self._on_preview_file_written)
        self._preview_tabs.openFileRequested.connect(self._on_preview_open_requested)
        self._preview_tabs.dotplotRegionSelected.connect(self._on_dotplot_region_selected)
        preview_layout.addWidget(self._preview_tabs, stretch=1)
        self._preview_host.setMinimumHeight(120)
        self._preview_host.setMinimumWidth(200)

        from repeatviewer.widgets.dock_title_bar import BrowserDockTitleBar

        self._dock_folders = self._make_browser_dock("Folders", "dock_folders", folder_side)
        self._dock_files = self._make_browser_dock("Files", "dock_files", files_panel)
        self._dock_preview = self._make_browser_dock(
            "Preview", "dock_preview", self._preview_host, dock_cls=PreviewDockWidget
        )

        self._dock_files.setTitleBarWidget(
            BrowserDockTitleBar(self._dock_files, "Files", trailing=files_trailing)
        )
        # Open-file tabs beside "Preview"; Fit / zoom / Reset on the same title row.
        preview_file_tabs_host = QWidget()
        preview_file_tabs_host.setObjectName("PreviewDockFileTabs")
        preview_file_tabs_layout = QHBoxLayout(preview_file_tabs_host)
        preview_file_tabs_layout.setContentsMargins(0, 0, 0, 0)
        preview_file_tabs_layout.setSpacing(0)
        preview_file_tabs_host.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Preferred,
        )
        preview_file_tabs_host.setMinimumWidth(0)
        self._preview_file_tabs_host = preview_file_tabs_host
        preview_zoom_host = QWidget()
        preview_zoom_host.setObjectName("PreviewDockZoomHost")
        preview_zoom_host.setFixedHeight(18)
        preview_zoom_host.setSizePolicy(
            QSizePolicy.Policy.Fixed,
            QSizePolicy.Policy.Fixed,
        )
        preview_zoom_layout = QHBoxLayout(preview_zoom_host)
        preview_zoom_layout.setContentsMargins(0, 0, 0, 0)
        preview_zoom_layout.setSpacing(6)
        preview_zoom_layout.addWidget(self._preview_zoom_bar)
        # AnnotationZoomBar sits on PreviewHeaderBar next to Text/Annotation.
        preview_title_bar = BrowserDockTitleBar(
            self._dock_preview,
            "Preview",
            leading=preview_file_tabs_host,
            leading_stretch=1,
            trailing=preview_zoom_host,
        )
        # Use a real dock title bar (same as Folders/Files) so Qt can drag, float,
        # and redock Preview. A ZeroDockTitleBar left nothing to grab.
        self._dock_preview.setTitleBarWidget(preview_title_bar)
        self._preview_title_bar = preview_title_bar
        self._dock_folders.setTitleBarWidget(
            BrowserDockTitleBar(self._dock_folders, "Folders", trailing=folders_trailing)
        )

        # Central filler: when a saved layout leaves a dock gap, paint it with the
        # window color instead of showing uninitialized GPU memory (static/snow).
        central = QWidget(self)
        central.setObjectName("BrowserCentral")
        central.setAutoFillBackground(True)
        window_bg = self.palette().color(QPalette.ColorRole.Window)
        central_palette = central.palette()
        central_palette.setColor(QPalette.ColorRole.Window, window_bg)
        central.setPalette(central_palette)
        central.setSizePolicy(
            QSizePolicy.Policy.Fixed,
            QSizePolicy.Policy.Expanding,
        )
        central.setMinimumWidth(0)
        central.setMaximumWidth(0)
        self.setCentralWidget(central)
        self._browser_panel = central

        for dock in (self._dock_folders, self._dock_files, self._dock_preview):
            dock.setContentsMargins(0, 0, 0, 0)
            # Do not zero the panel widget's contents margins — that pulls the
            # first folder/file row up under the custom title bar.
            title = dock.titleBarWidget()
            if title is not None:
                title.raise_()

        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self._dock_folders)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self._dock_files)
        self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self._dock_preview)
        self.splitDockWidget(
            self._dock_folders,
            self._dock_files,
            Qt.Orientation.Horizontal,
        )
        self.splitDockWidget(
            self._dock_files,
            self._dock_preview,
            Qt.Orientation.Vertical,
        )
        self._dock_folders.setMinimumWidth(100)
        self._dock_files.setMinimumWidth(200)
        self._dock_preview.setMinimumWidth(200)

    def _mount_preview_title_file_tabs(self) -> None:
        """Attach mirrored file tabs to the Preview dock title (after the window exists)."""
        import os

        if os.environ.get("REPEATVIEWER_TITLE_FILE_TABS", "1").lower() in {"0", "false", "no"}:
            return
        host = getattr(self, "_preview_file_tabs_host", None)
        if host is not None and self._preview_tabs is not None:
            self._preview_tabs.mount_file_tab_bar(host)

    def _apply_default_browser_layout(self, *, force: bool = False) -> None:
        if self._dock_folders is None or self._dock_files is None or self._dock_preview is None:
            return
        width = max(self.width(), _BROWSER_DOCK_WIDTH)
        if force:
            folder_w = _BROWSER_FOLDER_TREE_WIDTH
        else:
            folder_w = min(
                max(_BROWSER_FOLDER_TREE_WIDTH, width // 6),
                _BROWSER_FOLDER_TREE_WIDTH + 40,
            )
        if force:
            for dock in (self._dock_folders, self._dock_files, self._dock_preview):
                dock.setFloating(False)
                dock.show()
            self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self._dock_folders)
            self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self._dock_files)
            self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self._dock_preview)
            self.splitDockWidget(
                self._dock_folders,
                self._dock_files,
                Qt.Orientation.Horizontal,
            )
            self.splitDockWidget(
                self._dock_files,
                self._dock_preview,
                Qt.Orientation.Vertical,
            )
            self.resizeDocks([self._dock_folders], [folder_w], Qt.Orientation.Horizontal)
            height = max(self.height(), 600)
            files_h, preview_h = _BROWSER_FILES_PREVIEW_SPLIT
            total_h = max(files_h + preview_h, height - 120)
            files_h = max(140, int(total_h * 0.32))
            preview_h = max(200, total_h - files_h)
            self.resizeDocks(
                [self._dock_files, self._dock_preview],
                [files_h, preview_h],
                Qt.Orientation.Vertical,
            )
        else:
            # Ensure panels remain visible after restore; do not override sizes.
            for dock in (self._dock_folders, self._dock_files, self._dock_preview):
                if not dock.isVisible() and not dock.isFloating():
                    dock.show()
        self._ensure_preview_visible()

    def _reset_browser_layout(self) -> None:
        settings = QSettings()
        settings.beginGroup("mainwindow")
        settings.remove("toolbar_state_v5")
        settings.remove("toolbar_state_v4")
        settings.remove("browser_main_split")
        settings.remove("browser_files_preview_split")
        settings.endGroup()
        self._purge_legacy_docks()
        self._apply_default_browser_layout(force=True)
        self._ensure_preview_visible()

    @staticmethod
    def _settings_menu_hint() -> str:
        if platform.system() == "Darwin":
            return f"{APP_DISPLAY_NAME} → Settings"
        return "File → Settings"

    @staticmethod
    def _genome_menu_hint() -> str:
        return "File → Define genome"

    def _populate_align_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        del skip_default
        menu.addAction("MAFFT", self._run_mafft_default)
        menu.addAction("MUSCLE", self._run_muscle_default)
        menu.addAction("FastGA", self._run_fastga_default)
        menu.addAction(
            "Align all selected files in one MSA",
            self._align_selected_combined,
        )
        menu.addSeparator()
        menu.addAction("Change alignment parameters", self._modify_align_parameters)

    def _new_sequence(self) -> None:
        """Type or paste a sequence and save it into the current folder."""
        from repeatviewer.backends.sequence_compose import (
            SequenceComposeError,
            write_new_sequence,
        )
        from repeatviewer.views.tool_dialogs import NewSequenceDialog

        folder = self._current_folder
        if folder is None:
            QMessageBox.information(
                self, "New sequence", "Open a folder first, so the file has somewhere to go."
            )
            return
        dialog = NewSequenceDialog(parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            path = write_new_sequence(folder, **dialog.values())
        except SequenceComposeError as exc:
            QMessageBox.warning(self, "New sequence", str(exc))
            return
        self._refresh_after_compose([path], f"Created {path.name}")

    def _group_sequences(self, files: list | None = None) -> None:
        """Combine FASTA files into one multi-record file.

        ``files`` comes from the right-click menu; without it the current
        selection is used, which is how the Sequence menu calls it.
        """
        from repeatviewer.backends.sequence_compose import (
            SequenceComposeError,
            default_group_name,
            group_sequences,
            safe_file_stem,
        )

        candidates = [Path(p) for p in files] if files else self._selected_files()
        files = [
            p for p in candidates
            if p.is_file() and p.suffix.lower() in _FASTA_SUFFIXES
        ]
        if len(files) < 2:
            QMessageBox.information(
                self,
                "Combine Sequences to One File",
                "Select two or more FASTA files to combine into one file.",
            )
            return
        default_name = default_group_name(files)
        dialog = QInputDialog(self)
        dialog.setWindowTitle("Combine Sequences to One File")
        dialog.setLabelText(
            f"Combine {len(files)} sequences to one file.\n\nName for the new file:"
        )
        dialog.setTextValue(default_name)
        dialog.setInputMode(QInputDialog.InputMode.TextInput)
        # QInputDialog otherwise shrinks; pin a width that fits long default stems.
        name_px = dialog.fontMetrics().horizontalAdvance(default_name) + 48
        dialog.resize(max(640, name_px + 120), dialog.sizeHint().height())
        dialog.setMinimumWidth(max(640, name_px + 120))
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        name = dialog.textValue()
        if not name.strip():
            return
        out_path = files[0].parent / f"{safe_file_stem(name)}.fasta"
        try:
            written, count = group_sequences(files, out_path)
        except SequenceComposeError as exc:
            QMessageBox.warning(self, "Combine Sequences to One File", str(exc))
            return
        self._refresh_after_compose(
            [written],
            f"Combined {count} sequences into {written.name}",
            source_files=files,
        )

    def _overview_extract_pick_records(
        self, source: Path
    ) -> list[tuple[str, str, int]] | None:
        """Reuse Overview IDs/lengths when that FASTA is already open in Preview."""
        tabs = self._preview_tabs
        if tabs is None:
            return None
        wanted = source.expanduser().resolve()
        for index in range(tabs.count()):
            widget = tabs.widget(index)
            overview = getattr(widget, "_annotation_overview", None)
            if overview is None:
                continue
            path = getattr(overview, "_path", None)
            if path is None:
                continue
            try:
                if Path(path).expanduser().resolve() != wanted:
                    continue
            except OSError:
                continue
            record_ids = list(getattr(overview, "_record_ids", None) or [])
            lengths = dict(getattr(overview, "_lengths", None) or {})
            if not record_ids:
                return None
            return [
                (rid, "", int(lengths.get(rid, 0)))
                for rid in record_ids
            ]
        return None

    def _extract_sequences_from_list(
        self,
        source_path: object = None,
        *,
        initial_ids: str = "",
        initial_column: int = 1,
    ) -> None:
        """Pull chosen records out of a multi-FASTA, leaving the original alone.

        Two ways to choose: tick boxes, or a list of sequence names. The pick
        list reuses Overview when that file is open, otherwise streams ids and
        lengths (same scan Overview uses). Very large read sets fall back to
        the ID tab.
        """
        from repeatviewer.backends.sequence_compose import (
            SequenceComposeError,
            list_records,
        )
        from repeatviewer.views.tool_dialogs import ExtractSequencesDialog

        if isinstance(source_path, (str, Path)):
            files = [Path(source_path)]
        else:
            files = [
                p for p in self._selected_files()
                if p.is_file() and _is_extract_source(p)
            ]
        if len(files) != 1:
            QMessageBox.information(
                self,
                "Extract Sequences from File",
                "Select one FASTA file containing several sequences.",
            )
            return
        source = files[0].expanduser().resolve()
        if not source.is_file():
            QMessageBox.warning(
                self, "Extract Sequences from File", f"{source.name} was not found."
            )
            return

        records: list[tuple[str, str, int]] | None = None
        records_note = ""
        records = self._overview_extract_pick_records(source)
        if records is None:
            try:
                QApplication.setOverrideCursor(Qt.CursorShape.WaitCursor)
                try:
                    records = list_records(source)
                finally:
                    QApplication.restoreOverrideCursor()
            except SequenceComposeError as exc:
                if not initial_ids:
                    QMessageBox.warning(self, "Extract Sequences from File", str(exc))
                    return
                records = None
                records_note = str(exc)

        if records is not None and len(records) > _EXTRACT_PICK_LIST_MAX_RECORDS:
            records_note = (
                f"{source.name} has {len(records):,} sequences — too many for the "
                "pick list. Use the ID list tab."
            )
            records = None

        dialog = ExtractSequencesDialog(
            records=records,
            source_name=source.name,
            initial_ids=initial_ids,
            initial_column=initial_column,
            records_note=records_note,
            parent=self,
        )
        self._extract_dialog = dialog
        self._extract_source = source
        dialog.check_requested.connect(self._check_extract_ids)
        try:
            accepted = dialog.exec() == QDialog.DialogCode.Accepted
        finally:
            self._extract_dialog = None
        if not accepted:
            return

        if dialog.use_id_list():
            self._run_id_list_extraction(source, dialog)
        else:
            self._run_pick_list_extraction(source, dialog)

    def _run_pick_list_extraction(self, source: Path, dialog) -> None:
        """The original tick-box path: small files, fast enough to do inline."""
        from repeatviewer.backends.sequence_compose import (
            SequenceComposeError,
            extract_records,
            proposed_extract_path,
            proposed_single_record_path,
        )
        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        selected = list(dialog.selected_ids())
        if not selected:
            return
        single_file = dialog.single_file()
        output_path: Path | None = None
        try:
            if single_file:
                proposed = proposed_extract_path(
                    source, selected, source.parent, single_file=True
                )
                confirmed = confirm_existing_file_path(
                    self,
                    proposed,
                    title="Extract Sequences from File",
                    informative=(
                        "Overwrite the existing file, or save under a different name?"
                    ),
                )
                if confirmed is None:
                    return
                output_path = confirmed
            elif len(selected) == 1:
                proposed = proposed_single_record_path(source, selected[0], source.parent)
                confirmed = confirm_existing_file_path(
                    self,
                    proposed,
                    title="Extract Sequences from File",
                    informative=(
                        "Overwrite the existing file, or save under a different name?"
                    ),
                )
                if confirmed is None:
                    return
                output_path = confirmed

            written = extract_records(
                source,
                set(selected),
                source.parent,
                single_file=single_file,
                output_path=output_path,
            )
        except SequenceComposeError as exc:
            QMessageBox.warning(self, "Extract Sequences from File", str(exc))
            return
        names = ", ".join(p.name for p in written[:4])
        if len(written) > 4:
            names += ", …"
        self._refresh_after_compose(
            written,
            f"Extracted {len(written)} file(s): {names}",
            source_files=[source],
        )

    def _run_id_list_extraction(self, source: Path, dialog) -> None:
        """ID-list path: one streaming pass in the background, with cancel."""
        from repeatviewer.backends.sequence_compose import proposed_extract_path
        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        parsed = dialog.parsed_ids()
        if not parsed.ids:
            return
        single_file = dialog.single_file()
        exclude = dialog.exclude()

        if not single_file and (exclude or len(parsed.ids) > _EXTRACT_MANY_FILES_WARNING):
            expected = (
                "one file per sequence that is NOT in the list"
                if exclude
                else f"up to {len(parsed.ids):,} files"
            )
            answer = QMessageBox.question(
                self,
                "Extract Sequences from File",
                f"This will write {expected} into:\n{source.parent}\n\n"
                "Tick \"Write all selected into one file\" to get a single "
                "multi-FASTA instead.\n\nContinue?",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
            if answer != QMessageBox.StandardButton.Yes:
                return

        output_path: Path | None = None
        if single_file and not exclude:
            proposed = proposed_extract_path(
                source, list(parsed.ids), source.parent, single_file=True
            )
            confirmed = confirm_existing_file_path(
                self,
                proposed,
                title="Extract Sequences from File",
                informative=(
                    "Overwrite the existing file, or save under a different name?"
                ),
            )
            if confirmed is None:
                return
            output_path = confirmed

        self._start_extract_job(
            source,
            parsed.ids,
            count_only=False,
            single_file=single_file,
            exclude=exclude,
            output_path=output_path,
        )

    @Slot(object)
    def _check_extract_ids(self, ids: object) -> None:
        """Check button: count matches without writing anything."""
        source = self._extract_source
        if source is None or self._extract_thread is not None:
            return
        self._start_extract_job(source, list(ids or []), count_only=True)

    def _start_extract_job(
        self,
        source: Path,
        ids: list[str],
        *,
        count_only: bool,
        single_file: bool = True,
        exclude: bool = False,
        output_path: Path | None = None,
    ) -> None:
        if not ids:
            return
        if self._extract_thread is not None:
            # Nearly always a Check still reading a large file when OK was pressed.
            if not count_only:
                QMessageBox.information(
                    self,
                    "Extract Sequences from File",
                    "The ID check is still running. Wait for it to finish, or press "
                    "Stop, then try again.",
                )
            return
        self._job_cancel.reset()
        self._extract_count_only = count_only
        self._extract_source = source
        title = "Checking sequence IDs" if count_only else "Extracting sequences"
        if not count_only:
            self._clear_job_file_highlights()
        self._set_status_message(f"{title} in {source.name}…", 0)
        if not count_only:
            self._job_running = True
            self._job_progress_epoch = getattr(self, "_job_progress_epoch", 0) + 1
            self._job_work_steps = 100
            self._set_job_actions_enabled(False)
            if self._task_progress is not None:
                # +1 reserved until extracted sequences are opened in preview.
                self._task_progress.start(
                    self._running_job_title("Extract sequences"), 101
                )
            if self._stop_btn is not None:
                self._stop_btn.setEnabled(True)

        self._extract_thread = QThread(self)
        self._extract_worker = ExtractSequencesWorker(
            source,
            ids,
            source.parent,
            count_only=count_only,
            single_file=single_file,
            exclude=exclude,
            output_path=output_path,
            cancel=lambda: self._job_cancel.requested,
        )
        self._extract_worker.moveToThread(self._extract_thread)
        self._extract_thread.started.connect(self._extract_worker.run)
        queued = Qt.ConnectionType.QueuedConnection
        self._extract_worker.progress.connect(self._on_extract_progress, queued)
        self._extract_worker.finished.connect(self._on_extract_finished, queued)
        self._extract_worker.failed.connect(self._on_extract_failed, queued)
        self._extract_worker.finished.connect(self._extract_thread.quit)
        self._extract_worker.failed.connect(self._extract_thread.quit)
        self._extract_thread.finished.connect(self._cleanup_extract_worker)
        self._extract_thread.start()

    @Slot(int, str)
    def _on_extract_progress(self, percent: int, text: str) -> None:
        if self._extract_count_only:
            if self._extract_dialog is not None:
                self._extract_dialog.set_check_progress(f"Checking… {text}")
            return
        if self._task_progress is not None:
            work = getattr(self, "_job_work_steps", 100)
            capped = min(max(0, int(percent)), work)
            self._task_progress.update(capped, text)

    @Slot(object)
    def _on_extract_finished(self, summary: object) -> None:
        source = self._extract_source
        name = source.name if source else "the file"
        if getattr(summary, "cancelled", False):
            if self._extract_count_only and self._extract_dialog is not None:
                self._extract_dialog.set_check_result("Check stopped.", ok=False)
            if not self._extract_count_only:
                self._hide_job_progress()
                self._set_job_actions_enabled(True)
            self._set_status_message("Extract sequences stopped", 6000)
            return

        if self._extract_count_only:
            found = summary.found_ids
            total = summary.requested_ids
            text = f"{found:,} of {total:,} IDs found in {name}."
            if summary.missing:
                shown = ", ".join(summary.missing[:5])
                if len(summary.missing) > 5:
                    shown += ", …"
                text += f"\nNot found: {shown}"
            if self._extract_dialog is not None:
                self._extract_dialog.set_check_result(text, ok=not summary.missing)
            self._set_status_message(f"{found:,} of {total:,} IDs found", 8000)
            return

        self._finish_background_job()
        detail = ""
        try:
            if not summary.written:
                self._set_status_message("No sequences extracted", 6000)
                return

            message = f"Wrote {summary.matched:,} sequence(s)"
            if len(summary.written) > 1:
                message += f" to {len(summary.written):,} file(s)"
            message += "."
            if summary.missing:
                message += f" {len(summary.missing):,} ID(s) not found."
            outputs = list(summary.written)
            if summary.missing_path is not None:
                outputs.append(summary.missing_path)
            sources = [self._extract_source] if getattr(self, "_extract_source", None) else []
            self._refresh_after_compose(outputs, message, source_files=sources)

            detail = message
            if summary.missing_path is not None:
                detail += f"\n\nThe missing names are listed in:\n{summary.missing_path.name}"
        finally:
            self._schedule_after_job_ui_ready()

        if not summary.written:
            QMessageBox.warning(
                self,
                "Extract Sequences from File",
                f"None of the {summary.requested_ids:,} IDs matched a sequence in "
                f"{name}.\n\nSequence names are matched up to the first space, "
                "the same text BLAST reports.",
            )
        elif detail:
            QMessageBox.information(self, "Extract Sequences from File", detail)

    @Slot(str)
    def _on_extract_failed(self, message: str) -> None:
        if not self._extract_count_only:
            self._hide_job_progress()
            self._set_job_actions_enabled(True)
        if message == CANCELLED_MESSAGE:
            self._set_status_message("Extract sequences stopped", 6000)
            return
        if self._extract_count_only and self._extract_dialog is not None:
            self._extract_dialog.set_check_result(message, ok=False)
            return
        QMessageBox.warning(self, "Extract Sequences from File", message)
        self._set_status_message("Extract sequences failed", 6000)

    def _cleanup_extract_worker(self) -> None:
        if self._job_cancel.requested:
            self._job_cancel.reset()
        if not self._extract_count_only:
            self._job_running = self._worker_thread is not None
            if self._worker_thread is None:
                self._set_job_actions_enabled(True)
            if self._stop_btn is not None and not self._job_running:
                self._stop_btn.setEnabled(False)
        thread = self._extract_thread
        if thread is not None and thread.isRunning():
            thread.quit()
            thread.wait(5000)
        if self._extract_worker is not None:
            self._extract_worker.deleteLater()
            self._extract_worker = None
        if self._extract_thread is not None:
            self._extract_thread.deleteLater()
            self._extract_thread = None

    def _extract_sequences_from_hit_table(self, blast_txt: object = None) -> None:
        """Right-click a *_b.txt: extract the genome sequences it hit.

        Reads the ``sseqid`` column itself, so no ID list has to be made by
        hand. The subject FASTA is the genome the BLAST ran against -- taken
        from the ``*.lnk.fa`` reference link written beside the results, then
        the project genome, and only then asked for.
        """
        from repeatviewer.backends.sequence_compose import (
            SequenceComposeError,
            read_id_file,
        )

        if isinstance(blast_txt, (str, Path)):
            table = Path(blast_txt).expanduser().resolve()
        else:
            selected = [p for p in self._selected_files() if is_blast_table_file(p)]
            if len(selected) != 1:
                QMessageBox.information(
                    self,
                    "Extract subject sequences",
                    "Select one BLAST hit table (*_b.txt).",
                )
                return
            table = selected[0].resolve()

        try:
            parsed = read_id_file(table, column=_BLAST_SUBJECT_COLUMN)
        except SequenceComposeError as exc:
            QMessageBox.warning(self, "Extract subject sequences", str(exc))
            return
        if not parsed.ids:
            QMessageBox.information(
                self,
                "Extract subject sequences",
                f"No subject names found in {table.name}.",
            )
            return

        subject = self._subject_fasta_for_hit_table(table)
        if subject is None:
            return
        self._extract_sequences_from_list(
            subject, initial_ids="\n".join(parsed.ids), initial_column=1
        )

    def _subject_fasta_for_hit_table(self, table: Path) -> Path | None:
        """The FASTA a hit table was searched against, asking only if needed."""
        links = sorted(table.parent.glob("*.lnk.fa"))
        for link in links:
            try:
                target = link.resolve()
            except OSError:
                continue
            if target.is_file():
                return target

        genome = getattr(self._project, "genome_file", None)
        if genome is not None and Path(genome).is_file():
            return Path(genome).expanduser().resolve()

        path, _ = QFileDialog.getOpenFileName(
            self,
            f"Genome or read FASTA searched by {table.name}",
            str(table.parent),
            "FASTA (*.fa *.fasta *.fas *.fna *.faa *.gz);;All files (*)",
        )
        if not path:
            return None
        return Path(path).expanduser().resolve()

    def _refresh_after_compose(
        self,
        created: list[Path],
        message: str,
        *,
        source_files: list[Path] | None = None,
    ) -> None:
        if self._current_folder:
            # Leave opening/selection to ``_auto_open_outputs`` so the new
            # file (not the input) becomes the deep-blue opened row.
            self._populate_file_list(
                self._current_folder,
                new_outputs=created,
                source_files=list(source_files or []),
                open_preview=False,
            )
        if created:
            self._auto_open_outputs(created)
        self._set_status_message(message, 8000)

    def _populate_fastga_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        if not skip_default:
            menu.addAction("Run FastGA with saved parameters", self._run_fastga_default)
        menu.addSeparator()
        menu.addAction("Change FastGA parameters", self._open_fastga_dialog)

    def _fastga_candidate_files(self) -> list[Path]:
        """Selected FASTA files suitable for FastGA (0, 1, or 2+). No UI."""
        selected = self._selected_files()
        return [
            p for p in selected
            if p.is_file() and p.suffix.lower() in _FASTA_SUFFIXES
        ]

    def _selected_files_for_fastga(self, options=None, files=None) -> list[Path]:
        from repeatviewer.backends.fastga_dotplot import fastga_available, install_hint

        if not fastga_available():
            QMessageBox.warning(
                self,
                "FastGA not found",
                f"FastGA is not installed or not on PATH.\n\n{install_hint()}",
            )
            return []
        opts = options if options is not None else self._fastga_options
        mode = getattr(opts, "mode", "self")
        files = list(files) if files else self._fastga_candidate_files()
        if mode == FASTGA_MODE_REFERENCE:
            from repeatviewer.backends.fastga_dotplot import reference_queries

            reference = Path(opts.reference_fasta) if opts.reference_fasta else None
            if reference is None or not reference.is_file():
                QMessageBox.information(
                    self,
                    "FastGA",
                    "Choose the reference first: Tools → Align → Change alignment "
                    "parameters → FastGA tab, then click the reference in Files."
                    if reference is None
                    else f"The reference was not found:\n{reference}",
                )
                return []
            if not reference_queries(files, reference):
                QMessageBox.information(
                    self,
                    "FastGA",
                    f"Select the FASTA files to compare with {reference.name} "
                    "(the reference itself is skipped).",
                )
                return []
            return files
        if not files:
            QMessageBox.information(
                self,
                "FastGA",
                "Select one or more FASTA files: each is compared against itself. "
                "For other modes see Tools → Align → Change alignment parameters "
                "→ FastGA tab.",
            )
            return []
        if mode == FASTGA_MODE_PAIRS and len(files) < 2:
            QMessageBox.information(
                self,
                "FastGA",
                "All pairs needs two or more selected FASTA files.\n"
                "Select more files, or choose another mode in the FastGA tab.",
            )
            return []
        return files

    def _run_fastga_default(self) -> None:
        opts = self._fastga_options
        if getattr(opts, "mode", "self") == FASTGA_MODE_REFERENCE and not (
            opts.reference_fasta and Path(opts.reference_fasta).is_file()
        ):
            # Nothing to run against yet: show the tab where it is chosen.
            self._open_align_dialog(initial_tab=ALIGNER_FASTGA)
            return
        files = self._selected_files_for_fastga(opts)
        if files:
            self._start_fastga_job(files, opts)

    def _quick_fastga_right_click(self) -> None:
        """Right-click FastGA: run or reopen dialog per Options tab."""
        if self._fastga_options.right_click_action == FASTGA_RIGHT_CLICK_DIALOG:
            self._open_fastga_dialog()
            return
        self._run_fastga_default()

    def _open_fastga_dialog(self) -> None:
        """FastGA button / Tools → FastGA: the Align window on its FastGA tab."""
        self._open_align_dialog(initial_tab=ALIGNER_FASTGA)

    # ------------------------------------------------------------------ Tree

    def _populate_tree_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        if not skip_default:
            menu.addAction("Tree with default parameters", self._run_tree_default)
        menu.addSeparator()
        # Settings can be changed without a selection (Run then asks for one).
        menu.addAction(
            "Change Tree parameters", lambda: self._open_tree_dialog(require_selection=False)
        )
        menu.addSeparator()
        menu.addAction(
            "Methylation-masked sensitivity run (CG and CHG)", self._run_tree_methylation_masked
        )
        menu.addAction("RT-domain tree (REXdb lineages)", self._run_tree_rt_domains)

    def _run_tree_methylation_masked(self) -> None:
        """Saved settings, plus a second tree without CG/CHG columns, compared in the report."""
        from dataclasses import replace

        from repeatviewer.backends.methylation_mask import MASK_CG_CHG

        files = self._selected_files_for_tree()
        if files:
            self._start_tree_job(files, replace(self._tree_options, methylation_mask=MASK_CG_CHG))

    def _tree_candidate_files(self) -> list[Path]:
        """Selected FASTA files (several files are combined into one tree). No UI."""
        return [
            p for p in self._selected_files()
            if p.is_file() and p.suffix.lower() in _FASTA_SUFFIXES
        ]

    def _selected_files_for_tree(self) -> list[Path]:
        files = self._tree_candidate_files()
        if not files:
            QMessageBox.information(
                self,
                "Tree",
                "Select one or more FASTA files in the file panel first.\n"
                "Several selected files are combined into one tree.",
            )
        return files

    def _run_tree_default(self) -> None:
        files = self._selected_files_for_tree()
        if files:
            self._start_tree_job(files, self._tree_options)

    def _quick_tree_right_click(self) -> None:
        """Right-click Tree: run with saved settings or reopen the dialog (Options tab)."""
        if self._tree_options.right_click_action == TREE_RIGHT_CLICK_DIALOG:
            self._open_tree_dialog()
            return
        self._run_tree_default()

    def _open_tree_dialog(self, *, require_selection: bool = True) -> None:
        """Click Tree → parameters dialog (Cancel / Save / Run).

        Without a FASTA selection the toolbar click only shows a message; the
        Tools → Tree → Change Tree parameters entry still opens the dialog."""
        from repeatviewer.views.tree_dialog import TreeParamsDialog

        files = self._tree_candidate_files()
        if require_selection and not files:
            self._selected_files_for_tree()  # shows the "select a FASTA file" message
            return
        dialog = TreeParamsDialog(
            initial=self._tree_options,
            files=files,
            parent=self,
        )

        def _accepted(d) -> None:
            self._tree_options = d.options()
            save_tree_options(self._tree_options)
            if not d.should_run():
                self._set_status_message("Tree settings saved", 4000)
                return
            # Non-modal: clicking REXdb in Files changes the selection, so the
            # files chosen when the dialog opened are used.
            run_files = files or [
                p for p in self._selected_files_for_tree()
                if str(p) != self._tree_options.rexdb_fasta
            ]
            if run_files:
                self._start_tree_job(run_files, self._tree_options)

        # Non-modal so a click in Files can fill the REXdb path (apply_path_from_panel).
        self._present_te_viz_dialog(dialog, _accepted)

    def _run_tree_rt_domains(self) -> None:
        """Saved settings, with the sequences replaced by their RT domains (REXdb)."""
        from repeatviewer.backends.rt_domains import looks_like_rexdb
        from repeatviewer.job_options import TREE_SEQ_RT

        rexdb = self._tree_options.rexdb_fasta
        if not rexdb or not looks_like_rexdb(Path(rexdb)):
            QMessageBox.information(
                self,
                "RT-domain tree",
                "Choose REXdb first: Tree parameters → Alignment tab → REXdb.\n\n"
                "Download Viridiplantae_v4.0.fasta and Viridiplantae_v4.0.classification "
                "from github.com/repeatexplorer/rexdb into one folder.",
            )
            self._open_tree_dialog(require_selection=False)
            return
        files = self._selected_files_for_tree()
        if files:
            self._start_tree_job(files, replace(self._tree_options, sequence_mode=TREE_SEQ_RT))

    def _start_tree_job(self, files: list[Path], options) -> None:
        title = (f"Tree ({files[0].name})" if len(files) == 1
                 else f"Tree ({len(files)} files combined)")
        worker = TreeWorker(files, options, cancel=self._job_cancel)
        self._attach_process_monitor(
            worker,
            title=title,
            enabled=bool(options.show_process_window),
            tool_key="tree",
            mirror_progress=False,
        )
        self._start_background_job(
            worker,
            job_title=title,
            total_steps=3,
            source_files=files,
            on_finished=self._on_tree_finished,
            on_failed=self._on_tree_failed,
        )

    @Slot(object)
    def _on_tree_finished(self, result) -> None:
        report = Path(result["report"])
        created = [Path(p) for p in result["files"] if Path(p).is_file()]
        if self._current_folder:
            # Only the tree sits next to the input (green) and opens in the tree
            # viewer (deep blue); the report and supporting files are in the
            # new <input>_tree_details folder.
            self._refresh_subfolders_under(self._current_folder)
            self._populate_file_list(
                self._current_folder,
                new_outputs=created,
                source_files=getattr(self, "_last_job_source_files", []),
                open_preview=False,
            )
        tree = next((p for p in created if p.name.endswith(".treefile")
                     and not p.name.endswith("_masked.treefile")), None)
        opened = tree if tree is not None else (report if report.is_file() else None)
        if opened is not None:
            self._auto_open_outputs([opened])
        self._set_status_message(
            f"Tree completed — {opened.name if opened else report.name} "
            f"(report in {report.parent.name})", 12000, override_file_link=True
        )

    @Slot(str)
    def _on_tree_failed(self, error: str) -> None:
        if error.strip() == CANCELLED_MESSAGE:
            return
        QMessageBox.warning(self, "Tree failed", error.strip() or "Unknown error.")

    # ------------------------------------------------------------------ Map to Reference

    _MAP_SUFFIXES = (".fa", ".fasta", ".fna", ".fas", ".fsa", ".fq", ".fastq")

    def _populate_map_menu(self, menu: QMenu) -> None:
        menu.addAction("Map with default parameters", self._run_map_default)
        menu.addSeparator()
        # Settings can be changed without a selection (Run then asks for one).
        menu.addAction(
            "Change Map to Reference parameters",
            lambda: self._open_map_dialog(require_selection=False),
        )

    def _map_candidate_files(self) -> list[Path]:
        """Selected FASTA / FASTQ files (gzipped too), without the reference. No UI."""
        reference = str(self._map_options.reference_fasta or "")
        out = []
        for p in self._selected_files():
            name = p.name.lower()
            if name.endswith(".gz"):
                name = name[:-3]
            if p.is_file() and name.endswith(self._MAP_SUFFIXES) and str(p) != reference:
                out.append(p)
        return out

    def _selected_files_for_map(self) -> list[Path]:
        files = self._map_candidate_files()
        if not files:
            QMessageBox.information(
                self,
                "Map to Reference",
                "Select the FASTA or FASTQ files to map in the file panel first.\n"
                "The reference is chosen in the Map to Reference window.",
            )
        return files

    def _run_map_default(self) -> None:
        if not self._map_options.reference_fasta:
            self._open_map_dialog()
            return
        files = self._selected_files_for_map()
        if files:
            self._start_map_job(files, self._map_options)

    def _quick_map_right_click(self) -> None:
        """Right-click: map with saved settings, or open the dialog (Options tab)."""
        if self._map_options.right_click_action == MAP_RIGHT_CLICK_DIALOG:
            self._open_map_dialog()
            return
        self._run_map_default()

    def _open_map_dialog(self, *, require_selection: bool = True) -> None:
        """Map to Reference → parameters dialog (Cancel / Save / Run), non-modal so a
        click in Files fills the Reference field."""
        from repeatviewer.views.map_reference_dialog import MapReferenceDialog

        files = self._map_candidate_files()
        if require_selection and not files:
            self._selected_files_for_map()  # shows the "select files" message
            return
        dialog = MapReferenceDialog(initial=self._map_options, files=files, parent=self)

        def _accepted(d) -> None:
            self._map_options = d.options()
            save_map_options(self._map_options)
            if not d.should_run():
                self._set_status_message("Map to Reference settings saved", 4000)
                return
            reference = self._map_options.reference_fasta
            run_files = [p for p in (files or self._selected_files_for_map())
                         if str(p) != reference]
            if not reference:
                QMessageBox.information(self, "Map to Reference",
                                        "Choose a reference first (Mapping tab).")
                return
            if run_files:
                self._start_map_job(run_files, self._map_options)

        self._present_te_viz_dialog(dialog, _accepted)

    def _start_map_job(self, files: list[Path], options) -> None:
        ref = Path(options.reference_fasta).name if options.reference_fasta else "reference"
        title = (f"Map {files[0].name} to {ref}" if len(files) == 1
                 else f"Map {len(files)} files to {ref}")
        worker = MapReferenceWorker(files, options, cancel=self._job_cancel)
        self._attach_process_monitor(
            worker,
            title=title,
            enabled=bool(options.show_process_window),
            tool_key="map_reference",
            mirror_progress=False,
        )
        self._start_background_job(
            worker,
            job_title=title,
            total_steps=4,
            source_files=files,
            on_finished=self._on_map_finished,
            on_failed=self._on_map_failed,
        )

    @Slot(object)
    def _on_map_finished(self, result) -> None:
        rvmap = Path(result["rvmap"])
        created = [Path(p) for p in result["files"] if Path(p).is_file()]
        if self._current_folder:
            # The .rvmap sits next to the input (green) and opens in the contig view
            # (deep blue); everything else is in <input>_map_details.
            self._refresh_subfolders_under(self._current_folder)
            self._populate_file_list(
                self._current_folder,
                new_outputs=created,
                source_files=getattr(self, "_last_job_source_files", []),
                open_preview=False,
            )
        if rvmap.is_file():
            self._auto_open_outputs([rvmap])
        self._set_status_message(
            f"Map to Reference completed — {rvmap.name} (report in {rvmap.stem}_details)",
            12000, override_file_link=True,
        )

    @Slot(str)
    def _on_map_failed(self, error: str) -> None:
        if error.strip() == CANCELLED_MESSAGE:
            return
        QMessageBox.warning(self, "Map to Reference failed", error.strip() or "Unknown error.")

    def _on_preview_open_requested(self, path: object) -> None:
        """A result file asked for from a preview (contig view → Results menu)."""
        if not isinstance(path, Path) or not path.is_file():
            return
        if path.suffix.lower() in _FASTA_SUFFIXES:
            self._open_in_sequence_viewer(path)
        else:
            self._open_preview_tab(path, replace=False)

    def _attach_process_monitor(
        self,
        worker,
        *,
        title: str,
        enabled: bool,
        tool_key: str,
        mirror_progress: bool = True,
    ) -> None:
        """Capture command output; show live when enabled, always show on failure."""
        from repeatviewer.widgets.process_monitor_dialog import (
            ProcessMonitorBridge,
            ProcessMonitorDialog,
        )

        attr = f"_{tool_key}_process_dialog"
        dialog = getattr(self, attr, None)
        if dialog is None:
            short = title.split("(")[0].strip() or tool_key
            dialog = ProcessMonitorDialog(short, parent=self)
            setattr(self, attr, dialog)

        buffer: list[str] = []
        setattr(self, f"_{tool_key}_process_log", buffer)

        if enabled:
            dialog.prepare(title)

        # Bound methods on a main-thread QObject + QueuedConnection — never
        # lambdas — so dialog.show()/append never run on the worker QThread
        # (macOS aborts when NSPanel is touched off the GUI thread).
        bridge_attr = f"_{tool_key}_process_bridge"
        old_bridge = getattr(self, bridge_attr, None)
        if old_bridge is not None:
            old_bridge.deleteLater()
        bridge = ProcessMonitorBridge(
            dialog,
            buffer,
            live=enabled,
            title=title,
            cancelled_message=CANCELLED_MESSAGE,
            parent=self,
        )
        setattr(self, bridge_attr, bridge)
        queued = Qt.ConnectionType.QueuedConnection
        if hasattr(worker, "log_line"):
            worker.log_line.connect(bridge.on_line, queued)
        if mirror_progress:
            worker.progress.connect(bridge.on_progress, queued)
        worker.finished.connect(bridge.on_finished, queued)
        worker.failed.connect(bridge.on_failed, queued)

    def _warn_folders_with_spaces(self, folders: list[Path]) -> None:
        """Remind the user when imported paths contain spaces; offer conversion."""
        from repeatviewer.backends.safe_cli_paths import (
            SpacedNameHit,
            find_names_with_spaces,
            rename_spaced_paths,
        )

        roots: list[Path] = []
        hits: list[SpacedNameHit] = []
        any_truncated = False
        seen: set[str] = set()
        for folder in folders:
            try:
                resolved = Path(folder).expanduser().resolve()
            except OSError:
                continue
            roots.append(resolved)
            found, truncated = find_names_with_spaces(resolved, limit=500)
            any_truncated = any_truncated or truncated
            for hit in found:
                key = str(hit.path)
                if key in seen:
                    continue
                seen.add(key)
                hits.append(hit)

        if not hits:
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Spaces in folder or file names")
        dialog.setMinimumWidth(520)
        dialog.setMinimumHeight(420)
        layout = QVBoxLayout(dialog)

        intro = QLabel(
            "Some folder or file names contain spaces. Import is allowed, but "
            "tools such as BLAST and TE-Aid can fail on paths with spaces."
        )
        intro.setWordWrap(True)
        layout.addWidget(intro)
        tip = QLabel(
            "Recommendation: replace each space with an underscore (_).\n"
            "Example: TEtrimmer output → TEtrimmer_output\n"
            "(one space → _ two spaces → __ and so on)"
        )
        tip.setWordWrap(True)
        layout.addWidget(tip)

        list_label = QLabel("Folders or files have space in their name:")
        list_label.setWordWrap(True)
        layout.addWidget(list_label)

        list_widget = QListWidget()
        list_widget.setSelectionMode(QListWidget.SelectionMode.NoSelection)
        for hit in hits:
            list_widget.addItem(hit.display)
        if any_truncated:
            list_widget.addItem("… (more names not shown)")
        layout.addWidget(list_widget, stretch=1)

        buttons = QDialogButtonBox()
        convert_btn = buttons.addButton(
            "Convert spaces to _", QDialogButtonBox.ButtonRole.AcceptRole
        )
        keep_btn = buttons.addButton(
            "Keep as is", QDialogButtonBox.ButtonRole.RejectRole
        )
        layout.addWidget(buttons)

        def _on_convert() -> None:
            renames, errors = rename_spaced_paths([hit.path for hit in hits])
            if errors and not renames:
                QMessageBox.warning(
                    dialog,
                    "Convert spaces",
                    "Could not rename any paths:\n\n" + "\n".join(errors[:20]),
                )
                return
            if renames:
                self._apply_space_renames_to_workspace(roots, renames)
            dialog.accept()
            if renames:
                lines = "\n".join(
                    f"• {old.name} → {new.name}" for old, new in renames[:40]
                )
                extra = ""
                if len(renames) > 40:
                    extra = f"\n… and {len(renames) - 40} more"
                err_text = ""
                if errors:
                    err_text = "\n\nSome renames failed:\n" + "\n".join(errors[:10])
                QMessageBox.information(
                    self,
                    "Spaces converted",
                    f"Renamed {len(renames)} folder(s)/file(s) "
                    f"(spaces → underscores):\n\n{lines}{extra}{err_text}",
                )
            elif errors:
                QMessageBox.warning(
                    dialog,
                    "Convert spaces",
                    "Could not rename:\n\n" + "\n".join(errors[:20]),
                )

        convert_btn.clicked.connect(_on_convert)
        keep_btn.clicked.connect(dialog.reject)
        dialog.exec()

    def _apply_space_renames_to_workspace(
        self,
        roots: list[Path],
        renames: list[tuple[Path, Path]],
    ) -> None:
        """Update folder-tree paths and project pointers after space conversion."""
        from repeatviewer.backends.safe_cli_paths import apply_path_renames

        def _map(path: Path | None) -> Path | None:
            if path is None:
                return None
            try:
                return apply_path_renames(Path(path), renames)
            except OSError:
                return path

        # Remap top-level roots (import path ancestors may have moved).
        for index in range(self._folder_tree.topLevelItemCount()):
            item = self._folder_tree.topLevelItem(index)
            raw = item.data(0, ROLE_PATH)
            if not raw:
                continue
            old = Path(str(raw))
            new = _map(old)
            if new is None or new == old:
                continue
            item.setData(0, ROLE_PATH, str(new))
            counts = self._folder_panel_counts(new) if new.is_dir() else (0, 0)
            item.setText(0, self._folder_tree_label(new, counts=counts))
            item.setToolTip(
                0, f"{new}\n{counts[0]} folders · {counts[1]} files"
            )
            # Drop stale children; rebuild if expanded.
            while item.childCount():
                item.removeChild(item.child(0))
            if new.is_dir():
                self._set_folder_expandable(item, new, counts=counts)
                if item.isExpanded() and counts[0] > 0:
                    self._populate_folder_tree_children(item)
            self._sync_folder_item_icon(item)

        if self._project.input_dir is not None:
            self._project.input_dir = _map(self._project.input_dir)
        if self._current_folder is not None:
            mapped = _map(self._current_folder)
            if mapped is not None and mapped.is_dir():
                self._on_folder_selected(mapped)
        if self._project.genome_file is not None:
            self._project.genome_file = _map(self._project.genome_file)
        if self._project.pan_genomes:
            updated = []
            for entry in self._project.pan_genomes:
                new_path = _map(Path(entry.path))
                updated.append(
                    type(entry)(name=entry.name, path=new_path or entry.path)
                )
            self._project.pan_genomes = updated

        self._sync_folder_tree_labels()
        self._persist_workspace()
        # roots unused except for documentation of call sites; keep API clear.
        _ = roots

    def _rename_paths_with_spaces(
        self,
        paths: list[Path],
        *,
        kind: str = "file",
    ) -> list[Path]:
        """Rename basenames that contain spaces to use ``_`` and remind the user."""
        from repeatviewer.backends.safe_cli_paths import (
            basename_has_spaces,
            rename_spaces_in_basename,
        )

        prepared: list[Path] = []
        renames: list[tuple[str, str]] = []
        for raw in paths:
            path = Path(raw).expanduser().resolve()
            if not basename_has_spaces(path):
                prepared.append(path)
                continue
            try:
                new_path, renamed = rename_spaces_in_basename(path)
            except OSError as exc:
                QMessageBox.warning(
                    self,
                    "Rename required",
                    f"Could not rename “{path.name}” (spaces break many "
                    f"command-line tools):\n{exc}",
                )
                prepared.append(path)
                continue
            prepared.append(new_path)
            if renamed:
                renames.append((path.name, new_path.name))
                if self._project.genome_file is not None:
                    try:
                        if path == self._project.genome_file.expanduser().resolve():
                            self._project.genome_file = new_path
                    except OSError:
                        pass
                if self._project.pan_genomes:
                    updated = []
                    for entry in self._project.pan_genomes:
                        try:
                            if path == Path(entry.path).expanduser().resolve():
                                updated.append(type(entry)(name=entry.name, path=new_path))
                            else:
                                updated.append(entry)
                        except OSError:
                            updated.append(entry)
                    self._project.pan_genomes = updated

        if renames:
            lines = "\n".join(f"• {old} → {new}" for old, new in renames)
            QMessageBox.information(
                self,
                "Spaces in file name",
                "Spaces in file names can break command-line tools.\n"
                f"Renamed {len(renames)} {kind}(s) (spaces → underscores):\n\n"
                f"{lines}\n\n"
                "Analysis will continue with the new name(s).",
            )
            if self._current_folder is not None:
                self._populate_file_list(self._current_folder)
            self._persist_workspace()
        return prepared

    def _ensure_genome_paths_without_spaces(self) -> None:
        """Rename configured genome FASTA(s) if their basenames contain spaces."""
        targets: list[Path] = []
        for _label, path in self._project.genome_targets():
            targets.append(Path(path))
        if not targets:
            return
        self._rename_paths_with_spaces(targets, kind="genome")

    def _start_fastga_job(self, files: list[Path], options) -> None:
        from repeatviewer.backends.fastga_dotplot import (
            fastga_step_count,
            reference_queries,
        )

        mode = getattr(options, "mode", "self")
        total_steps = fastga_step_count(files, options)
        source_files = list(files)
        if mode == FASTGA_MODE_PAIRS and len(files) >= 2:
            pair_count = len(files) * (len(files) - 1) // 2
            title = f"FastGA ({pair_count} pair(s) → multi)"
        elif mode == FASTGA_MODE_REFERENCE:
            reference = Path(options.reference_fasta)
            queries = reference_queries(files, reference)
            files = queries
            # Queries and reference are both inputs (blue highlight, rule 4).
            source_files = [*queries, reference]
            title = (
                f"FastGA ({queries[0].name} vs {reference.name})"
                if len(queries) == 1
                else f"FastGA ({len(queries)} files vs {reference.name})"
            )
        elif len(files) == 1:
            title = f"FastGA ({files[0].name} self)"
        else:
            title = f"FastGA ({len(files)} self-comparisons)"
        show_log = bool(getattr(options, "show_process_window", False))
        worker = FastGaWorker(
            files,
            options,
            cancel=self._job_cancel,
            stream_log=show_log,
        )
        self._attach_process_monitor(
            worker,
            title=title,
            enabled=show_log,
            tool_key="fastga",
            mirror_progress=False,
        )
        self._start_background_job(
            worker,
            job_title=title,
            total_steps=total_steps,
            source_files=source_files,
            on_finished=self._on_fastga_finished,
            on_failed=self._on_fastga_failed,
        )

    @Slot(list)
    def _on_fastga_finished(self, created: list) -> None:
        created_paths = [Path(p) for p in created if Path(p).is_file()]
        if self._current_folder:
            self._populate_file_list(
                self._current_folder,
                new_outputs=created_paths,
                source_files=getattr(self, "_last_job_source_files", []),
                open_preview=False,
            )
        if created_paths:
            self._auto_open_outputs(created_paths)
        names = ", ".join(p.name for p in created_paths)
        self._set_status_message(f"FastGA completed — {names}", 12000, override_file_link=True)

    @Slot(str)
    def _on_fastga_failed(self, error: str) -> None:
        from repeatviewer.backends.fastga_dotplot import install_hint

        body = error.strip()
        if "not installed" in body.lower() and "conda" not in body.lower():
            body += f"\n\n{install_hint()}"
        QMessageBox.warning(self, "FastGA failed", body)

    def _open_repeat_quantify(self) -> None:
        """Quantify a tandem-repeat / satellite library across a target FASTA."""
        from repeatviewer.views.tool_dialogs import RepeatQuantifyDialog

        start_dir = str(self._current_folder) if self._current_folder else ""

        def _on_accept(dialog: RepeatQuantifyDialog) -> None:
            save_repeat_quantify_options(dialog.options())
            if not dialog.should_run():
                self._set_status_message("Repeat quantify settings saved", 4000)
                return
            library = dialog.library()
            targets = dialog.targets()
            engine = dialog.engine()
            params = dialog.params()
            if library is None or not targets:
                return
            if engine in ("both", "repeatmasker"):
                from repeatviewer.backends.repeat_quantify import (
                    repeatmasker_executable_available,
                )

                if not repeatmasker_executable_available():
                    QMessageBox.warning(
                        self,
                        "RepeatMasker not found",
                        "RepeatMasker is not on PATH. Install it, or choose the "
                        "BLASTn engine instead.",
                    )
                    return
            from repeatviewer.app_preferences import load_blast_db_base_dir
            from repeatviewer.database_settings import (
                reload_blast_databases_into_project,
            )

            worker = RepeatQuantifyWorker(
                Path(library),
                targets,
                engine=engine,
                threads=self._project.batch_workers,
                params=params,
                cancel=self._job_cancel,
                output_base=dialog.output_dir(),
                stamp=dialog.run_stamp(),
                blast_databases=reload_blast_databases_into_project(self._project),
                save_blast_index=dialog.save_blast_index(),
                blast_db_base_dir=load_blast_db_base_dir(),
                use_system_blast=self._project.use_system_blast,
            )
            worker.database_created.connect(self._register_blast_database)
            self._start_background_job(
                worker,
                job_title="Repeat quantify",
                # 1 = the app's plain "Repeat quantify… <genome>" status line; a
                # count > 1 would switch the status bar to the batch progress bar.
                total_steps=1,
                source_files=targets,
                on_finished=self._on_repeat_quantify_finished,
                on_failed=self._on_repeat_quantify_failed,
            )

        from repeatviewer.database_settings import reload_blast_databases_into_project

        dialog = RepeatQuantifyDialog(
            start_dir=start_dir,
            initial=load_repeat_quantify_options(),
            blast_databases=reload_blast_databases_into_project(self._project),
            parent=self,
        )
        # Same non-modal presentation as the TE visualization dialogs, so the
        # Folders and Files panels stay clickable while it is open.
        self._present_te_viz_dialog(dialog, _on_accept)

    @Slot(object)
    def _register_blast_database(self, db) -> None:
        """Save a BLAST database a tool built, so Add/Manage and every tool see it."""
        from repeatviewer.database_settings import (
            load_blast_databases_from_settings,
            save_blast_databases_to_settings,
        )

        databases = load_blast_databases_from_settings()
        prefix = str(Path(db.db_prefix_path))
        if any(str(Path(d.db_prefix_path)) == prefix for d in databases):
            return
        databases.append(db)
        save_blast_databases_to_settings(databases)
        self._project.blast_databases = databases

    @Slot(object)
    def _on_repeat_quantify_finished(self, outputs) -> None:
        paths = [Path(p) for p in outputs if Path(p).is_file()]
        if self._current_folder:
            # Each run creates a new timestamped folder; show it in the tree.
            self._refresh_subfolders_under(self._current_folder)
            self._populate_file_list(
                self._current_folder,
                new_outputs=paths,
                source_files=getattr(self, "_last_job_source_files", []),
            )
        for index, path in enumerate(paths):
            self._open_preview_tab(path, replace=(index == 0))
        self._set_status_message(
            f"Repeat quantify — {len(paths)} report(s) written",
            12000,
            override_file_link=True,
        )

    @Slot(str)
    def _on_repeat_quantify_failed(self, error: str) -> None:
        QMessageBox.warning(
            self, "Repeat quantify failed", error.strip() or "Unknown error."
        )

    def _populate_blast_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        menu.addAction(
            "BLASTn",
            lambda checked=False: self._run_blast_program(BLAST_PROGRAM_BLASTN),
        )
        menu.addAction(
            "BLASTp",
            lambda checked=False: self._run_blast_program(BLAST_PROGRAM_BLASTP),
        )
        menu.addAction(
            "BLASTx",
            lambda checked=False: self._run_blast_program(BLAST_PROGRAM_BLASTX),
        )
        menu.addSeparator()
        menu.addAction("Change BLAST parameters", self._change_blast_parameters)

    def _populate_te_visualization_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        if not skip_default:
            menu.addAction("TE density plot", self._open_te_density_plot_dialog)
        menu.addAction("TE divergence plot", self._open_te_divergence_plot_dialog)
        menu.addAction("TE abundance plot", self._open_te_abundance_plot_dialog)
        menu.addAction("MSA plot", self._open_msa_plot_dialog)

    def _populate_benchmarking_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        if not skip_default:
            menu.addAction("Confusion matrix benchmarker", self._open_confusion_matrix_benchmark)
        menu.addAction("TE library benchmarker", self._open_te_library_benchmark)

    def _populate_extend_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        if not skip_default:
            menu.addAction("MSA Extend with default parameters", self._run_extend)
        menu.addSeparator()
        menu.addAction("Change MSA Extend parameters", self._modify_extend_parameters)

    def _populate_tsd_search_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        if not skip_default:
            menu.addAction("TSD search with default parameters", self._run_tsd_search)
        menu.addSeparator()
        menu.addAction("Change TSD search parameters", self._modify_tsd_search_parameters)

    def _populate_teaid_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        if not skip_default:
            menu.addAction("TE-Aid with default parameters", self._run_teaid)
        menu.addSeparator()
        menu.addAction("Change TE-Aid parameters", self._modify_teaid_parameters)

    def _populate_msa_clean_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        del skip_default
        menu.addAction("Combined MSA clean", self._run_default_msa_clean)
        menu.addSeparator()
        self._add_msa_clean_menu_action(
            menu,
            "Remove gappy columns",
            "msa_clean_gappy",
        )
        self._add_msa_clean_menu_action(
            menu,
            "Clean row with divergence",
            "msa_clean_div",
        )
        self._add_msa_clean_menu_action(
            menu,
            "Clean row with gap",
            "msa_clean_gap",
        )
        menu.addSeparator()
        menu.addAction("Change MSA clean parameters", self._modify_msa_clean_parameters)

    def _populate_consensus_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        if not skip_default:
            menu.addAction("Consensus with default parameters", self._run_consensus)
        menu.addSeparator()
        menu.addAction("Change consensus parameters", self._modify_consensus_parameters)

    def _populate_flexidot_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        if not skip_default:
            menu.addAction("Flexidot with default parameters", self._run_flexidot_combined)
        menu.addSeparator()
        menu.addAction("Change Flexidot parameters", self._modify_flexidot_parameters)

    def _populate_orf_prediction_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        if not skip_default:
            menu.addAction(
                "ORF Prediction with default parameters",
                self._run_orf_prediction,
            )
        menu.addSeparator()
        menu.addAction(
            "Change ORF Prediction parameters",
            self._modify_orf_prediction_parameters,
        )

    def _populate_te_annotation_menu(self, menu: QMenu, *, skip_default: bool = False) -> None:
        if not skip_default:
            menu.addAction(
                "Sequence annotation with default parameters",
                self._run_te_annotation_saved,
            )
        menu.addSeparator()
        menu.addAction(
            "Change Sequence annotation parameters",
            self._modify_te_annotation_parameters,
        )

    def _build_menu(self) -> None:
        file_menu = self.menuBar().addMenu("&File")
        add_folder_action = file_menu.addAction("Add &Folder…", self._pick_and_add_folder)
        add_folder_action.setShortcut(QKeySequence("Ctrl+Shift+O"))
        add_folder_action.setStatusTip("Add a folder as a workspace root in the Folders panel")
        new_folder_action = file_menu.addAction("New &Folder…", lambda: self._create_folder(None))
        new_folder_action.setShortcut(QKeySequence("Ctrl+Shift+N"))
        new_folder_action.setStatusTip("Create a folder inside the current workspace folder")
        file_menu.addSeparator()
        genome_action = file_menu.addAction("Define &genome", self._open_set_genome)
        genome_action.setShortcut(QKeySequence("Ctrl+Shift+G"))
        genome_action.setStatusTip("Set reference or pan-genome FASTA files")
        define_db_action = file_menu.addAction("Define &database", self._open_define_database)
        define_db_action.setStatusTip("Configure CDD and similarity annotation databases")
        file_menu.addAction("Restore deleted &files", self._restore_removed_folders)
        file_menu.addSeparator()
        settings_action = QAction("&Setting", self)
        settings_action.setShortcut(QKeySequence.StandardKey.Preferences)
        settings_action.setStatusTip("Preview defaults, BLAST options, and batch workers")
        settings_action.triggered.connect(self._open_settings)
        settings_action.setMenuRole(
            QAction.MenuRole.PreferencesRole
            if platform.system() == "Darwin"
            else QAction.MenuRole.ApplicationSpecificRole
        )
        file_menu.addAction(settings_action)
        file_menu.addSeparator()
        file_menu.addAction("&Quit", self.close, "Ctrl+Q")

        seq_menu = self.menuBar().addMenu("&Sequence")
        seq_menu.addAction("New Sequence", self._new_sequence)
        seq_menu.addAction("Combine Sequences to One File", self._group_sequences)
        seq_menu.addAction("Extract Sequences from File", self._extract_sequences_from_list)
        seq_menu.addAction(
            "Extract subject sequences from hit table",
            self._extract_sequences_from_hit_table,
        )
        seq_menu.addSeparator()
        seq_menu.addAction(
            "Reverse Complement",
            lambda: self._apply_function("reverse_complement"),
        )
        self._populate_orf_prediction_menu(seq_menu.addMenu("ORF &Prediction"))
        seq_menu.addAction(
            "Translate",
            lambda: self._apply_function(
                "dna_translation",
                orf_prediction_settings=self._orf_prediction_settings,
            ),
        )

        tools_menu = self.menuBar().addMenu("&Tools")
        tools_menu.addAction("Open in &AliView", self._open_selected_in_aliview)
        tools_menu.addSeparator()
        self._populate_blast_menu(tools_menu.addMenu("&BLAST"))
        self._populate_align_menu(tools_menu.addMenu("&Align"))
        self._populate_fastga_menu(tools_menu.addMenu("Fast&GA"))
        self._populate_tree_menu(tools_menu.addMenu("T&ree"))
        self._populate_map_menu(tools_menu.addMenu("&Map to Reference"))
        self._populate_extend_menu(tools_menu.addMenu("MSA &Extend"))
        self._populate_teaid_menu(tools_menu.addMenu("&TE-Aid"))
        self._populate_msa_clean_menu(tools_menu.addMenu("MSA &Clean"))
        self._populate_consensus_menu(tools_menu.addMenu("&Consensus"))
        self._populate_flexidot_menu(tools_menu.addMenu("&Flexidot"))
        self._populate_te_annotation_menu(tools_menu.addMenu("Sequence &Annotation"))
        tools_menu.addAction("Repeat &Quantify", self._open_repeat_quantify)
        self._populate_te_visualization_menu(tools_menu.addMenu("&Visualization"))

        workflows_menu = QMenu("&Workflows", self)
        self.menuBar().addMenu(workflows_menu)
        self._workflows_menu = workflows_menu
        self._rebuild_workflows_menu()

        others_menu = QMenu("&Others", self)
        self.menuBar().addMenu(others_menu)
        self._others_menu = others_menu
        self._plugins_menu = others_menu  # rebuild helpers refresh Others
        self._rebuild_others_menu()

        view_menu = self.menuBar().addMenu("&View")
        if self._dock_folders is not None:
            view_menu.addAction(self._dock_folders.toggleViewAction())
        if self._dock_files is not None:
            view_menu.addAction(self._dock_files.toggleViewAction())
        if self._dock_preview is not None:
            view_menu.addAction(self._dock_preview.toggleViewAction())
        view_menu.addSeparator()
        view_menu.addAction(
            "Reset &layout",
            self._reset_browser_layout,
        )

        help_menu = self.menuBar().addMenu("&Help")
        help_menu.addAction(f"&About {APP_DISPLAY_NAME}", self._show_about)
        self._connect_menu_operation_logging()
        self._rebuild_top_toolbar()

    @staticmethod
    def _clean_menu_label(text: str) -> str:
        return text.replace("&", "").replace("…", "...").strip()

    def _connect_menu_operation_logging(self) -> None:
        """Log menu-launched operations without logging ordinary file previews."""

        def connect_menu(menu: QMenu, parents: list[str]) -> None:
            menu_label = self._clean_menu_label(menu.title())
            path = [*parents, menu_label] if menu_label else parents
            for action in menu.actions():
                submenu = action.menu()
                if submenu is not None:
                    connect_menu(submenu, path)
                    continue
                if action.isSeparator() or action.property("operationLogConnected"):
                    continue
                label = self._clean_menu_label(action.text())
                menu_path = " > ".join([*path, label])
                action.triggered.connect(
                    lambda checked=False, operation=menu_path: self._log_menu_operation(
                        operation
                    )
                )
                action.setProperty("operationLogConnected", True)

        for action in self.menuBar().actions():
            menu = action.menu()
            if menu is not None:
                connect_menu(menu, [])

    def _log_menu_operation(self, operation: str) -> None:
        try:
            selected = self._selected_files()
        except Exception:  # noqa: BLE001
            selected = []
        operation_logger.info(
            "MENU operation=%s selected=[%s] folder=%s",
            operation,
            ", ".join(str(path) for path in selected),
            str(self._current_folder) if self._current_folder is not None else "",
        )

    def _install_shortcuts(self) -> None:
        QShortcut(QKeySequence("Ctrl+Shift+G"), self, self._open_set_genome)
        QShortcut(QKeySequence.StandardKey.Refresh, self, self._refresh_current_folder)
        QShortcut(QKeySequence(Qt.Key.Key_F5), self, self._refresh_current_folder)
        QShortcut(QKeySequence.StandardKey.Find, self, self._focus_browser_filter)
        for key in (QKeySequence.StandardKey.Delete, QKeySequence(Qt.Key.Key_Backspace)):
            shortcut = QShortcut(key, self._file_list)
            shortcut.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            shortcut.activated.connect(self._delete_selected_shortcut)
        for key, handler in (
            (QKeySequence.StandardKey.Copy, self._copy_files_shortcut),
            (QKeySequence.StandardKey.Cut, self._cut_files_shortcut),
            (QKeySequence.StandardKey.Paste, self._paste_files_shortcut),
        ):
            shortcut = QShortcut(key, self._file_list)
            shortcut.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            shortcut.activated.connect(handler)
        paste_tree = QShortcut(QKeySequence.StandardKey.Paste, self._folder_tree)
        paste_tree.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        paste_tree.activated.connect(self._paste_files_shortcut)
        find_tree = QShortcut(QKeySequence.StandardKey.Find, self._folder_tree)
        find_tree.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        find_tree.activated.connect(lambda: self._folder_filter.setFocus())
        undo_files = QShortcut(QKeySequence.StandardKey.Undo, self._file_list)
        undo_files.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        undo_files.activated.connect(self._undo_browser_action)
        undo_folders = QShortcut(QKeySequence.StandardKey.Undo, self._folder_tree)
        undo_folders.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
        undo_folders.activated.connect(self._undo_browser_action)

    def _focus_browser_filter(self) -> None:
        focus = QApplication.focusWidget()
        if focus is not None and (
            focus is self._folder_tree
            or focus is self._folder_filter
            or self._folder_tree.isAncestorOf(focus)
            or self._folder_filter.isAncestorOf(focus)
            or self._dock_folders.isAncestorOf(focus)
        ):
            self._folder_filter.setFocus()
            self._folder_filter.selectAll()
            return
        self._file_filter.setFocus()
        self._file_filter.selectAll()

    def _apply_folder_filter(self) -> None:
        if not hasattr(self, "_folder_filter"):
            return
        self._folder_tree.apply_name_filter(self._folder_filter.text())

    def _copy_files_shortcut(self) -> None:
        paths = self._selected_files()
        if paths:
            self._copy_files_to_clipboard(paths)

    def _cut_files_shortcut(self) -> None:
        paths = self._selected_files()
        if paths:
            self._cut_files_to_clipboard(paths)

    def _paste_files_shortcut(self) -> None:
        if self._current_folder is None:
            return
        self._paste_files_to_folder(self._current_folder)

    def _delete_selected_shortcut(self) -> None:
        paths = self._selected_files()
        if paths:
            self._delete_selected_files(paths)

    def _apply_toolbar_theme(self) -> None:
        size = self._toolbar_font_size
        self._toolbar_font = QFont(self.font())
        self._toolbar_font.setPointSize(size)
        self._toolbar_font.setBold(False)
        from repeatviewer.ui.theme_mode import load_accent_color

        self.setStyleSheet(
            application_stylesheet(
                toolbar_font_pt=size,
                browser_panel_font_pt=self._browser_panel_font_size,
                accent=load_accent_color(),
            )
        )
        if hasattr(self, "_preview_zoom_bar"):
            self._preview_zoom_bar.refresh_slider_style()
        if hasattr(self, "_annotation_zoom_bar"):
            self._annotation_zoom_bar.refresh_slider_style()
        if self._functions_toolbar is not None:
            self._compact_functions_toolbar(self._functions_toolbar)

    def _apply_browser_panel_font(self) -> None:
        self._browser_panel_font = QFont(self.font())
        self._browser_panel_font.setPointSize(self._browser_panel_font_size)
        self._browser_panel_font.setBold(False)
        if hasattr(self, "_folder_tree"):
            self._folder_tree.setFont(self._browser_panel_font)
            self._sync_folder_tree_placeholders()
        if hasattr(self, "_file_list"):
            self._file_list.setFont(self._browser_panel_font)

    def _rebuild_top_toolbar(self) -> None:
        if self._active_toolbar_menu is not None:
            self._active_toolbar_menu.close()
        self._active_toolbar_menu = None
        self._active_toolbar_menu_button = None
        if self._functions_toolbar is not None:
            self._functions_toolbar.removeEventFilter(self)
            self.removeToolBar(self._functions_toolbar)
            self._functions_toolbar.deleteLater()
            self._functions_toolbar = None
        self._job_toolbar_actions = []
        self._toolbar_items = []
        self._overflow_hidden = []
        self._overflow_btn = None
        self._toolbar_container = None
        self._toolbar_hbox = None
        self._build_top_toolbar()

    def _toolbar_button_height(self) -> int:
        # Compact single-row height with a little vertical padding so labels are
        # comfortably centred.
        return QFontMetrics(self._toolbar_font).height() + 12

    def _compact_functions_toolbar(self, bar: QToolBar) -> None:
        bar.setContentsMargins(0, 0, 0, 0)
        btn_h = self._toolbar_button_height()
        bar.setFixedHeight(btn_h + 8)  # padding above/below the single row
        for widget in bar.findChildren(QToolButton):
            # The overflow arrow gets the full row height so its glyph can be big.
            widget.setFixedHeight(btn_h + 8 if widget is self._overflow_btn else btn_h)
            widget.setSizePolicy(
                QSizePolicy.Policy.Fixed,
                QSizePolicy.Policy.Fixed,
            )

    def _clear_toolbar_button_highlights(
        self,
        except_btn: QToolButton | None = None,
    ) -> None:
        """Drop stale pressed/focus chrome when another toolbar control is used."""
        for control in self._job_toolbar_actions:
            if not isinstance(control, QToolButton):
                continue
            if control is except_btn:
                continue
            control.setDown(False)
            control.clearFocus()

    def _style_toolbar_button(self, widget) -> None:
        widget.setFont(self._toolbar_font)
        widget.setFixedHeight(self._toolbar_button_height())
        widget.setAutoRaise(True)
        widget.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        # Fixed size: keeps a stable width so the custom overflow logic can
        # measure which buttons fit (the native » chevron does not render with
        # the current theme, so we manage overflow ourselves).
        widget.setSizePolicy(
            QSizePolicy.Policy.Fixed,
            QSizePolicy.Policy.Fixed,
        )

    def _release_toolbar_button(self, btn: QToolButton) -> None:
        btn.setDown(False)
        btn.clearFocus()
        QTimer.singleShot(
            0,
            lambda button=btn: (button.setDown(False), button.clearFocus()),
        )

    def _register_toolbar_item(self, widget, overflow_handler) -> None:
        self._toolbar_items.append((widget, overflow_handler))

    def _add_toolbar_text_action(
        self,
        bar,
        text: str,
        slot,
        *,
        job_action: bool = True,
        tooltip: str = "",
    ) -> QToolButton:
        btn = QToolButton(self._toolbar_container)
        btn.setText(text)
        if tooltip:
            btn.setToolTip(tooltip)
        self._style_toolbar_button(btn)
        btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)

        def _invoke_slot() -> None:
            self._clear_toolbar_button_highlights(btn)
            try:
                slot()
            finally:
                self._release_toolbar_button(btn)

        btn.clicked.connect(_invoke_slot)
        self._toolbar_hbox.addWidget(btn)
        if job_action:
            self._job_toolbar_actions.append(btn)

        def _overflow(menu: QMenu) -> None:
            act = menu.addAction(text, slot)
            act.setEnabled(btn.isEnabled())

        self._register_toolbar_item(btn, _overflow)
        return btn

    def _add_toolbar_plain_button(
        self,
        bar,
        text: str,
        slot,
        *,
        job_action: bool = True,
        tooltip: str = "",
    ) -> QToolButton:
        """Text-only toolbar button added as a widget (avoids macOS QAction separator dots)."""
        btn = QToolButton(self._toolbar_container)
        btn.setText(text)
        if tooltip:
            btn.setToolTip(tooltip)
        self._style_toolbar_button(btn)
        btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)

        def _invoke_slot() -> None:
            self._clear_toolbar_button_highlights(btn)
            try:
                slot()
            finally:
                self._release_toolbar_button(btn)

        btn.clicked.connect(_invoke_slot)
        self._toolbar_hbox.addWidget(btn)
        if job_action:
            self._job_toolbar_actions.append(btn)

        def _overflow(menu: QMenu) -> None:
            act = menu.addAction(text, slot)
            act.setEnabled(btn.isEnabled())

        self._register_toolbar_item(btn, _overflow)
        return btn

    def _add_toolbar_spacer(self, bar: QToolBar, width: int = 10) -> None:
        """Invisible gap between toolbar groups."""
        spacer = QWidget(self._toolbar_container)
        spacer.setFixedWidth(width)
        spacer.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Preferred)
        spacer.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        self._toolbar_hbox.addWidget(spacer)
        self._register_toolbar_item(spacer, None)  # no menu entry for spacers

    def _add_toolbar_separator(self, bar: QToolBar) -> None:
        """Visible vertical rule between toolbar groups (also in overflow menu)."""
        del bar
        wrap = QWidget(self._toolbar_container)
        wrap.setObjectName("ToolbarSeparator")
        wrap.setFixedWidth(11)
        wrap.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Preferred)
        row = QHBoxLayout(wrap)
        row.setContentsMargins(5, 5, 5, 5)
        row.setSpacing(0)
        line = QFrame(wrap)
        line.setFrameShape(QFrame.Shape.VLine)
        line.setFixedWidth(1)
        line.setStyleSheet("background: #b0b0b0; border: 0; max-width: 1px;")
        row.addWidget(line)
        self._toolbar_hbox.addWidget(wrap)

        def _overflow(menu: QMenu) -> None:
            menu.addSeparator()

        self._register_toolbar_item(wrap, _overflow)
    def _add_toolbar_text_menu_button(
        self,
        bar,
        text: str,
        default_slot,
        menu_builder,
        *,
        tooltip: str = "",
        run_label: str | None = None,
        show_default_action: bool = True,
        right_click_slot=None,
    ) -> QToolButton:
        """Text-only toolbar button; click opens its actions and optional default action.

        ``right_click_slot`` runs a second action on right-click (same pattern as
        ``_add_toolbar_action_button``).
        """
        del bar
        btn = _ToolbarActionButton(self._toolbar_container)
        btn.setText(text)
        if tooltip:
            btn.setToolTip(tooltip)
        self._style_toolbar_button(btn)
        btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)

        def _show_menu_on_click() -> None:
            active_menu = self._active_toolbar_menu
            active_btn = self._active_toolbar_menu_button
            # Second click on the same tool closes the menu.
            if active_menu is not None and active_btn is btn:
                active_menu.close()
                return
            if active_menu is not None:
                active_menu.close()
            self._clear_toolbar_button_highlights(btn)

            # Build actions via a temporary QMenu (structure only; not shown).
            scratch = QMenu(btn)
            if show_default_action:
                scratch.addAction(run_label or text, default_slot)
            menu_builder(scratch, skip_default=True)
            actions = list(scratch.actions())
            if not actions:
                return

            popup = _ToolbarPopup(actions)
            self._active_toolbar_menu = popup
            self._active_toolbar_menu_button = btn
            btn.setDown(True)

            def _menu_closed() -> None:
                if self._active_toolbar_menu is popup:
                    self._active_toolbar_menu = None
                    self._active_toolbar_menu_button = None
                self._release_toolbar_button(btn)

            popup.closed.connect(_menu_closed)
            popup.show_below(btn)

        btn._repeatviewer_open_menu = _show_menu_on_click  # type: ignore[attr-defined]
        btn.clicked.connect(_show_menu_on_click)

        if right_click_slot is not None:
            btn.setContextMenuPolicy(Qt.ContextMenuPolicy.PreventContextMenu)
            btn.rightClicked.connect(
                lambda s=right_click_slot: QTimer.singleShot(0, s)
            )

        self._toolbar_hbox.addWidget(btn)
        self._job_toolbar_actions.append(btn)

        def _overflow(menu: QMenu) -> None:
            sub = menu.addMenu(text)
            sub.setEnabled(btn.isEnabled())
            if show_default_action:
                sub.addAction(run_label or text, default_slot)
            menu_builder(sub, skip_default=True)

        self._register_toolbar_item(btn, _overflow)
        return btn

    def _add_toolbar_action_button(
        self,
        bar,
        text: str,
        slot,
        *,
        tooltip: str = "",
        right_click_slot=None,
    ) -> QToolButton:
        """Text-only toolbar button whose click runs a single action.

        ``right_click_slot`` runs a second action on right-click -- used by
        BLAST, where left-click opens the dialog and right-click runs the
        program chosen in the dialog's Right Click tab with saved settings.
        """
        del bar
        btn = _ToolbarActionButton(self._toolbar_container)
        btn.setText(text)
        if tooltip:
            btn.setToolTip(tooltip)
        self._style_toolbar_button(btn)
        btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        btn.clicked.connect(slot)
        self._toolbar_hbox.addWidget(btn)
        self._job_toolbar_actions.append(btn)

        if right_click_slot is not None:
            btn.setContextMenuPolicy(Qt.ContextMenuPolicy.PreventContextMenu)
            # Run on the next event-loop tick so the release event fully
            # finishes before a modal dialog (e.g. TE-Aid params) opens.
            btn.rightClicked.connect(
                lambda s=right_click_slot: QTimer.singleShot(0, s)
            )

        def _overflow(menu: QMenu) -> None:
            action = menu.addAction(text, slot)
            action.setEnabled(btn.isEnabled())

        self._register_toolbar_item(btn, _overflow)
        return btn

    def _build_top_toolbar(self) -> None:
        bar = QToolBar("Functions", self)
        bar.setObjectName("FunctionsToolBar")
        bar.setMovable(True)
        bar.setFloatable(True)
        bar.setAllowedAreas(
            Qt.ToolBarArea.TopToolBarArea
            | Qt.ToolBarArea.BottomToolBarArea
            | Qt.ToolBarArea.LeftToolBarArea
            | Qt.ToolBarArea.RightToolBarArea
        )
        bar.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        bar.setMovable(False)
        bar.setFloatable(False)
        self._functions_toolbar = bar

        # All buttons live in ONE container widget with our own horizontal layout.
        # Because the toolbar holds a single (expanding) widget, its native »
        # extension never triggers — we manage overflow entirely ourselves.
        container = QWidget(bar)
        container.setObjectName("FunctionsToolBarInner")
        container.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        hbox = QHBoxLayout(container)
        hbox.setContentsMargins(6, 0, 6, 0)
        hbox.setSpacing(1)
        self._toolbar_container = container
        self._toolbar_hbox = hbox

        for tool_id in self._toolbar_order:
            self._add_toolbar_tool(bar, tool_id)

        # Stretch keeps buttons left and pushes the overflow arrow to the right,
        # leaving empty space on the right for future buttons.
        hbox.addStretch(1)

        # Custom overflow button — a plain text triangle so it never depends on
        # the theme's (broken) icon engine.
        self._overflow_btn = QToolButton(container)
        self._overflow_btn.setObjectName("ToolbarOverflowButton")
        self._overflow_btn.setText("▼")  # filled triangle reads much larger than ▾
        self._overflow_btn.setToolTip("More tools")
        self._overflow_btn.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextOnly)
        self._overflow_btn.setAutoRaise(True)
        self._overflow_btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self._overflow_btn.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self._overflow_btn.setFixedHeight(self._toolbar_button_height())
        overflow_menu = QMenu(self._overflow_btn)
        overflow_menu.aboutToShow.connect(self._rebuild_overflow_menu)
        self._overflow_btn.setMenu(overflow_menu)
        self._overflow_btn.hide()
        hbox.addWidget(self._overflow_btn)

        bar.addWidget(container)
        self.addToolBar(Qt.ToolBarArea.TopToolBarArea, bar)

        self._compact_functions_toolbar(bar)
        container.installEventFilter(self)
        QTimer.singleShot(0, self._relayout_toolbar_overflow)

    def eventFilter(self, obj, event):  # noqa: N802 (Qt API)
        if obj is self._toolbar_container and event.type() == QEvent.Type.Resize:
            self._relayout_toolbar_overflow()
        return super().eventFilter(obj, event)

    def _relayout_toolbar_overflow(self) -> None:
        if getattr(self, "_in_overflow_relayout", False):
            return  # guard against reentrancy from show/hide during event delivery
        container = self._toolbar_container
        if container is None or self._overflow_btn is None or not self._toolbar_items:
            return
        self._in_overflow_relayout = True
        try:
            self._do_relayout_toolbar_overflow(container)
        finally:
            self._in_overflow_relayout = False

    def _do_relayout_toolbar_overflow(self, container) -> None:
        layout = self._toolbar_hbox
        spacing = max(0, layout.spacing()) if layout is not None else 0
        left, _t, right, _b = layout.getContentsMargins() if layout is not None else (0, 0, 0, 0)
        available = container.width() - left - right
        if available <= 0:
            return  # not laid out yet; a later resize will re-run this

        widths = [w.sizeHint().width() for (w, _) in self._toolbar_items]
        total = sum(widths) + spacing * (len(self._toolbar_items) - 1)

        if total <= available:
            for widget, _ in self._toolbar_items:
                widget.setVisible(True)
            self._overflow_btn.setVisible(False)
            self._overflow_hidden = []
            return

        self._overflow_btn.setVisible(True)
        budget = available - self._overflow_btn.sizeHint().width() - spacing
        used = 0
        hidden: list[tuple] = []
        for index, (widget, handler) in enumerate(self._toolbar_items):
            need = widths[index] + (spacing if index > 0 else 0)
            if used + need <= budget:
                widget.setVisible(True)
                used += need
            else:
                widget.setVisible(False)
                hidden.append((widget, handler))
        self._overflow_hidden = hidden

    def _rebuild_overflow_menu(self) -> None:
        if self._overflow_btn is None:
            return
        menu = self._overflow_btn.menu()
        if menu is None:
            return
        menu.clear()
        for _widget, handler in self._overflow_hidden:
            if handler is not None:
                handler(menu)

    def _add_toolbar_tool(self, bar, tool_id: str) -> None:
        tip = _TOOLBAR_TOOLTIPS.get(tool_id, "")
        if tool_id.startswith("_sep"):
            self._add_toolbar_separator(bar)
        elif tool_id == "define_genome":
            self._add_toolbar_plain_button(
                bar,
                "Define genome",
                self._open_set_genome,
                job_action=False,
                tooltip=tip,
            )
        elif tool_id == "new_sequence":
            self._add_toolbar_text_action(
                bar,
                "New Sequence",
                self._new_sequence,
                job_action=False,
                tooltip=tip,
            )
        elif tool_id == "group_sequences":
            self._add_toolbar_text_action(
                bar,
                "Combine Sequences",
                self._group_sequences,
                job_action=False,
                tooltip=tip,
            )
        elif tool_id == "extract_sequences":
            self._add_toolbar_text_action(
                bar,
                "Extract Sequences",
                self._extract_sequences_from_list,
                job_action=False,
                tooltip=tip,
            )
        elif tool_id == "reverse_complement":
            self._add_toolbar_text_action(
                bar,
                "Reverse Complement",
                lambda: self._apply_function("reverse_complement"),
                tooltip=tip,
            )
        elif tool_id == "orf_prediction":
            self._add_toolbar_action_button(
                bar,
                "ORF Prediction",
                self._open_orf_prediction_dialog,
                tooltip=(
                    "Click: set parameters and run.\n"
                    "Right-click: run with saved parameters."
                ),
                right_click_slot=self._run_orf_prediction,
            )
        elif tool_id == "translate":
            self._add_toolbar_text_action(
                bar,
                "Translate",
                lambda: self._apply_function(
                    "dna_translation",
                    orf_prediction_settings=self._orf_prediction_settings,
                ),
                tooltip=tip,
            )
        elif tool_id == "aliview":
            self._add_toolbar_text_action(
                bar,
                "AliView",
                self._open_selected_in_aliview,
                job_action=False,
                tooltip=tip,
            )
        elif tool_id == "align":
            self._add_toolbar_action_button(
                bar,
                "Align",
                self._open_align_dialog,
                tooltip=(
                    "Click: set parameters and run (MAFFT, MUSCLE or FastGA by tab).\n"
                    "Right-click: run the aligner chosen in the Options tab."
                ),
                right_click_slot=self._quick_align_right_click,
            )
        elif tool_id == "fastga":
            self._add_toolbar_action_button(
                bar,
                "FastGA",
                self._open_fastga_dialog,
                tooltip=(
                    "Click: the Align window on its FastGA tab.\n"
                    "Right-click: action chosen in the Options tab."
                ),
                right_click_slot=self._quick_fastga_right_click,
            )
        elif tool_id == "tree":
            self._add_toolbar_action_button(
                bar,
                "Tree",
                lambda: self._open_tree_dialog(),
                tooltip=(
                    "Click: set parameters and run.\n"
                    "Right-click: action chosen in the Options tab "
                    "(default: build with saved parameters)."
                ),
                right_click_slot=self._quick_tree_right_click,
            )
        elif tool_id == "map_reference":
            self._add_toolbar_action_button(
                bar,
                "Map to Reference",
                lambda: self._open_map_dialog(),
                tooltip=(
                    "Click: set parameters and run.\n"
                    "Right-click: action chosen in the Options tab "
                    "(default: map with saved parameters)."
                ),
                right_click_slot=self._quick_map_right_click,
            )
        elif tool_id == "repeat_quantify":
            self._add_toolbar_text_action(
                bar,
                "Repeat Quantify",
                self._open_repeat_quantify,
                tooltip="Quantify a repeat / satellite library across a genome.",
            )
        elif tool_id == "blast":
            self._add_toolbar_action_button(
                bar,
                "BLAST",
                self._open_blast_dialog,
                tooltip=(
                    "Click: set parameters and run.\n"
                    "Right-click: run the program chosen in the Right Click tab, "
                    "with saved settings."
                ),
                right_click_slot=self._quick_blast_right_click,
            )
        elif tool_id == "te_genome_map":
            self._add_toolbar_text_menu_button(
                bar,
                "Visualization",
                self._open_te_density_plot_dialog,
                self._populate_te_visualization_menu,
                tooltip=tip,
                run_label="TE density plot",
            )
        elif tool_id == "te_divergence_plot":
            self._add_toolbar_text_action(
                bar,
                "TE divergence",
                self._open_te_divergence_plot_dialog,
                tooltip=tip,
            )
        elif tool_id == "extend":
            self._add_toolbar_action_button(
                bar,
                "MSA Extend",
                self._open_extend_dialog,
                tooltip=(
                    "Click: set parameters and run.\n"
                    "Right-click: run with saved settings (set in the Options tab)."
                ),
                right_click_slot=self._quick_extend_right_click,
            )
        elif tool_id == "teaid":
            self._add_toolbar_action_button(
                bar,
                "TE-Aid",
                self._open_teaid_dialog,
                tooltip=(
                    "Click: set parameters and run.\n"
                    "Right-click: run with saved settings (set in the Options tab)."
                ),
                right_click_slot=self._quick_teaid_right_click,
            )
        elif tool_id == "msa_clean":
            self._add_toolbar_text_menu_button(
                bar,
                "MSA Clean",
                self._run_default_msa_clean,
                self._populate_msa_clean_menu,
                tooltip=(
                    "Click: choose Combined or a single clean step.\n"
                    "Right-click: run the action chosen in the Options tab."
                ),
                show_default_action=False,
                right_click_slot=self._quick_msa_clean_right_click,
            )
        elif tool_id == "consensus":
            self._add_toolbar_action_button(
                bar,
                "Consensus",
                self._open_consensus_dialog,
                tooltip=(
                    "Click: set parameters and run.\n"
                    "Right-click: run with saved settings (set in the Options tab)."
                ),
                right_click_slot=self._quick_consensus_right_click,
            )
        elif tool_id == "flexidot":
            self._add_toolbar_action_button(
                bar,
                "Flexidot",
                self._open_flexidot_dialog,
                tooltip=(
                    "Click: set parameters and run.\n"
                    "Right-click: run with saved settings (set in the Options tab)."
                ),
                right_click_slot=self._quick_flexidot_right_click,
            )
        elif tool_id == "te_annotation":
            self._add_toolbar_action_button(
                bar,
                "Sequence Annotation",
                self._open_te_annotation_dialog,
                tooltip=(
                    "Click: set modules and run.\n"
                    "Right-click: run with saved selection (set in the Options tab)."
                ),
                right_click_slot=self._quick_te_annotation_right_click,
            )
        elif tool_id.startswith("plugin:"):
            plugin_id = tool_id[len("plugin:") :]
            plugin = self._plugin_registry.get(plugin_id)
            if plugin is None:
                return
            self._add_toolbar_text_action(
                bar,
                plugin.label,
                lambda checked=False, pid=plugin.id: self._apply_plugin(pid),
                tooltip=plugin.tooltip or tip,
            )
        elif tool_id.startswith("workflow:"):
            workflow_id = tool_id[len("workflow:") :]
            from repeatviewer.workflows import get_workflow_spec

            spec = get_workflow_spec(workflow_id)
            if spec is None:
                return
            self._add_toolbar_text_action(
                bar,
                spec.name,
                lambda checked=False, wid=spec.id: self._run_workflow(wid),
                tooltip=tip or f"Run workflow {spec.name}",
            )
        elif tool_id == "workflows_manage":
            self._add_toolbar_text_action(
                bar,
                "Design and manage workflows",
                self._open_workflow_designer,
                job_action=False,
                tooltip="Design and manage your workflows",
            )

    def _select_all_fasta_in_folder(self) -> None:
        if self._file_list.count() == 0:
            return
        self._file_list.clearSelection()
        selected = 0
        for index in range(self._file_list.count()):
            item = self._file_list.item(index)
            raw = item.data(ROLE_PATH)
            if not raw:
                continue
            path = Path(str(raw))
            if path.suffix.lower() in _FASTA_SUFFIXES:
                item.setChecked(True)
                selected += 1
        self._set_status_message(f"Selected {selected} FASTA file(s)", 4000)

    def _reveal_in_file_manager(self, paths: list[Path]) -> None:
        if paths:
            reveal_in_file_manager(paths[0])

    def _genome_path_resolved(self) -> Path | None:
        if not self._project.genome_file:
            return None
        try:
            path = self._project.genome_file.expanduser().resolve()
        except OSError:
            return None
        return path if path.is_file() else None

    def _pan_genome_names(self) -> list[str]:
        if self._project.is_pan_genome_mode:
            return [entry.name for entry in self._project.pan_genomes]
        return []

    def _job_total_steps(
        self,
        func_name: str,
        seq_files: list[Path],
        job_kwargs: dict[str, Any],
    ) -> int:
        n = len(seq_files)
        blast_opts = job_kwargs.get("blast_options")
        if func_name == "blast" and blast_opts is not None and blast_opts.active_database_name:
            return n
        if not self._project.is_pan_genome_mode:
            if func_name == "flexidot" and n > 1 and job_kwargs.get("flexidot_combined", True):
                return 1
            if func_name.startswith("plugin:"):
                plugin = self._plugin_registry.get(func_name[len("plugin:") :])
                if plugin is not None and not plugin.supports_parallel:
                    return n
            return n
        pan_n = len(self._project.pan_genomes)
        if func_name in {"blast", "teaid"}:
            return n * pan_n
        if func_name == "flexidot":
            if n > 1 and job_kwargs.get("flexidot_combined", True):
                return pan_n
            return n * pan_n
        if func_name.startswith("plugin:"):
            plugin = self._plugin_registry.get(func_name[len("plugin:") :])
            if plugin is not None and not plugin.supports_parallel:
                return n
            return n
        return n

    def _prepare_configured_genomes(self) -> None:
        targets = self._project.genome_targets()
        if not targets:
            return
        self._ensure_genome_paths_without_spaces()
        targets = self._project.genome_targets()
        if not targets:
            return
        self._set_status_message("Preparing genome index and chromosome sizes…", 0)
        QApplication.processEvents()
        from repeatviewer.external_process import process_log_sink
        from repeatviewer.widgets.process_monitor_dialog import ProcessMonitorDialog

        log_lines: list[str] = []

        def _log(message: str) -> None:
            log_lines.append(str(message))

        messages: list[str] = []
        ok_all = True
        with process_log_sink(_log):
            for label, path in targets:
                ok, msg = ensure_genome_artifacts(
                    path,
                    use_system_blast=self._project.use_system_blast,
                )
                ok_all = ok_all and ok
                messages.append(f"{label}: {msg}" if label else msg)
        if ok_all:
            preview = "; ".join(messages[:3])
            if len(messages) > 3:
                preview += f" (+{len(messages) - 3} more)"
            self._set_status_message(preview, 10000)
        else:
            QMessageBox.warning(self, "Genome preparation", "\n".join(messages))
            if log_lines:
                dialog = getattr(self, "_genome_prepare_process_dialog", None)
                if dialog is None:
                    dialog = ProcessMonitorDialog("Genome preparation", parent=self)
                    self._genome_prepare_process_dialog = dialog
                dialog.show_failure_log(
                    "Genome preparation",
                    log_lines,
                    error="\n".join(messages),
                )
            self._set_status_message("Genome preparation incomplete", 8000)

    def _on_genome_changed(self, genome_before: Path | None, genome_after: Path | None) -> None:
        invalidate_genome_chromosome_cache()
        self._genome_chrom_names_cache = []
        self._genome_chrom_sizes_cache = {}
        self._genome_chrom_cache_path = None
        if (
            not self._project.is_pan_genome_mode
            and genome_after is not None
            and genome_after.is_file()
        ):
            new_chroms = set(chromosome_names_from_genome(genome_after))
            old_filter = self._te_genome_map_options.chromosomes_only
            if new_chroms and old_filter and not _map_chromosome_filter_to_genome(
                old_filter, list(new_chroms)
            ):
                self._te_genome_map_options = replace(
                    self._te_genome_map_options,
                    chromosomes_only=None,
                )
                save_te_genome_map_options(self._te_genome_map_options)
                self._set_status_message(
                    "Cleared TE genome map chromosome filter (previous genome).",
                    8000,
                )
        if self._project.has_genome_configured():
            self._prepare_configured_genomes()
        self._update_status()

    def _open_set_genome(self) -> None:
        if self._genome_dialog is not None:
            self._genome_dialog.raise_()
            self._genome_dialog.activateWindow()
            return
        self._persist_workspace()
        genome_before = self._genome_path_resolved()
        dialog = GenomeDialog(
            self._project,
            self,
            selected_paths=self._selected_genome_fastas(),
        )
        self._genome_dialog = dialog
        dialog.setModal(False)
        dialog.setWindowModality(Qt.WindowModality.NonModal)
        dialog.finished.connect(self._on_genome_dialog_finished)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _selected_genome_fastas(self) -> list[Path]:
        return [
            p.expanduser().resolve()
            for p in self._selected_files()
            if is_genome_fasta(p)
        ]

    def _on_genome_dialog_finished(self, result: int) -> None:
        dialog = self._genome_dialog
        self._genome_dialog = None
        if dialog is None or result != int(QDialog.DialogCode.Accepted):
            return
        genome_before = self._genome_path_resolved()
        self._project = dialog.apply_to_project()
        self._ensure_genome_paths_without_spaces()
        genome_after = self._genome_path_resolved()
        self._on_genome_changed(genome_before, genome_after)
        self._persist_workspace()
        if self._project.is_pan_genome_mode:
            self._set_status_message(
                f"Pan-genome mode: {len(self._project.pan_genomes)} genome(s)",
                8000,
            )
        elif genome_after is not None:
            self._set_status_message(f"Genome: {genome_after.name}", 8000)
        else:
            self._set_status_message("Genome cleared", 5000)

    def _open_define_database(self) -> None:
        if self._define_database_dialog is not None:
            self._define_database_dialog.raise_()
            self._define_database_dialog.activateWindow()
            return
        dialog = DefineDatabaseDialog(self._project, self)
        self._define_database_dialog = dialog
        dialog.setModal(False)
        dialog.setWindowModality(Qt.WindowModality.NonModal)
        dialog.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, True)
        dialog.finished.connect(self._on_define_database_finished)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()
        # Do not auto-fill from the currently selected Files-panel item.
        # The dialog should only fill paths on an explicit Files-panel click.

    def _on_define_database_finished(self, result: int) -> None:
        dialog = self._define_database_dialog
        self._define_database_dialog = None
        if dialog is None:
            return
        if result == int(QDialog.DialogCode.Accepted):
            self._project = dialog.project()
            self._persist_workspace()
            if self._preview_tabs is not None:
                self._preview_tabs.refresh_similarity_databases()
        dialog.deleteLater()

    def _open_settings(self) -> None:
        dialog = SettingsDialog(
            self._project,
            self,
            plugin_registry=self._plugin_registry,
            toolbar_order=self._toolbar_order,
        )
        if dialog.exec():
            from repeatviewer.ui.theme_mode import (
                apply_theme,
                load_accent_color,
                load_theme_mode,
            )

            old_appearance = (
                load_theme_mode(),
                load_accent_color(),
            )
            self._project = dialog.apply_to_project()
            self._align_options = load_align_options()
            self._persist_workspace()
            if dialog.appearance_signature() != old_appearance:
                from repeatviewer.widgets.browser_icons import clear_browser_icon_cache

                apply_theme(dialog.selected_theme_mode())
                clear_browser_icon_cache()
                # Re-layer the window's own accents on top of the new palette.
                self._apply_toolbar_theme()
                # Re-tint file-list highlights for the new light/dark palette.
                if self._current_folder is not None:
                    self._apply_file_filter()
                self._sync_folder_tree_placeholders()
            from repeatviewer.app_preferences import apply_log_level

            apply_log_level()
            old_font_size = self._toolbar_font_size
            old_browser_size = self._browser_panel_font_size
            self._toolbar_order, new_font_size = dialog.apply_toolbar_settings()
            self._toolbar_font_size = new_font_size
            new_browser_size = load_browser_panel_font_size()
            self._browser_panel_font_size = new_browser_size
            if new_font_size != old_font_size or new_browser_size != old_browser_size:
                self._apply_toolbar_theme()
            if new_browser_size != old_browser_size:
                self._apply_browser_panel_font()
            self._rebuild_top_toolbar()
            self._update_status()
            if self._preview_tabs is not None:
                from repeatviewer.widgets.preview_zoom import (
                    load_preview_fit_mode,
                    load_preview_font_size,
                    load_preview_zoom_percent,
                )

                zoom_bar = self._preview_tabs.zoom_bar
                kind = getattr(zoom_bar, "current_zoom_kind", lambda: None)()
                if kind is not None and hasattr(zoom_bar, "set_zoom_kind"):
                    zoom_bar.set_zoom_kind(kind)
                else:
                    zoom_bar.set_values(
                        load_preview_zoom_percent(),
                        load_preview_fit_mode(),
                    )
                if hasattr(zoom_bar, "_font_spin"):
                    zoom_bar._block = True
                    zoom_bar._font_spin.setValue(load_preview_font_size())
                    zoom_bar._block = False
                self._preview_tabs._apply_zoom_to_current_tab()
                self._preview_tabs.apply_preview_font_size(load_preview_font_size())
                self._preview_tabs.refresh_nucleotide_color_scheme()

    def _collect_root_folders(self) -> list[Path]:
        roots: list[Path] = []
        for index in range(self._folder_tree.topLevelItemCount()):
            item = self._folder_tree.topLevelItem(index)
            raw = item.data(0, ROLE_PATH)
            if not raw:
                continue
            path = Path(str(raw)).expanduser().resolve()
            if path.is_dir():
                roots.append(path)
        return roots

    def _schedule_persist_workspace(self) -> None:
        if self._folder_tree_restoring:
            return
        self._folder_persist_timer.start()

    def _collect_expanded_folders(self) -> list[Path]:
        expanded: list[Path] = []
        seen: set[str] = set()

        def walk(item: QTreeWidgetItem) -> None:
            if item.isExpanded():
                raw = item.data(0, ROLE_PATH)
                if raw:
                    path = Path(str(raw)).expanduser().resolve()
                    if path.is_dir():
                        key = str(path)
                        if key not in seen:
                            seen.add(key)
                            expanded.append(path)
            for index in range(item.childCount()):
                child = item.child(index)
                if self._is_lazy_folder_placeholder(child):
                    continue
                walk(child)

        for index in range(self._folder_tree.topLevelItemCount()):
            walk(self._folder_tree.topLevelItem(index))
        return expanded

    def _apply_expanded_folders(self, folders: list[Path]) -> None:
        if not folders:
            for index in range(self._folder_tree.topLevelItemCount()):
                self._folder_tree.topLevelItem(index).setExpanded(True)
            return

        ordered = sorted(
            {path.expanduser().resolve() for path in folders if path.is_dir()},
            key=lambda path: (len(path.parts), str(path).lower()),
        )
        self._folder_tree.blockSignals(True)
        self._folder_tree_restoring = True
        try:
            for index in range(self._folder_tree.topLevelItemCount()):
                self._folder_tree.topLevelItem(index).setExpanded(False)

            def collapse_children(item: QTreeWidgetItem) -> None:
                for child_index in range(item.childCount()):
                    child = item.child(child_index)
                    if self._is_lazy_folder_placeholder(child):
                        continue
                    child.setExpanded(False)
                    collapse_children(child)

            for index in range(self._folder_tree.topLevelItemCount()):
                collapse_children(self._folder_tree.topLevelItem(index))

            for folder in ordered:
                item = self._ensure_folder_tree_path(folder)
                if item is not None:
                    item.setExpanded(True)
        finally:
            self._folder_tree_restoring = False
            self._folder_tree.blockSignals(False)
        self._populate_expanded_folder_tree_items(ordered if folders else None)

    def _populate_expanded_folder_tree_items(self, folders: list[Path] | None) -> None:
        def needs_populate(item: QTreeWidgetItem) -> bool:
            if item.childCount() == 0:
                return True
            return item.childCount() == 1 and self._is_lazy_folder_placeholder(item.child(0))

        def populate_if_needed(item: QTreeWidgetItem) -> None:
            if not item.isExpanded():
                return
            raw = item.data(0, ROLE_PATH)
            if not raw:
                return
            folder = Path(str(raw)).expanduser().resolve()
            if folder.is_dir() and self._folder_has_subfolders(folder) and needs_populate(item):
                self._populate_folder_tree_children(item)

        def walk(item: QTreeWidgetItem) -> None:
            populate_if_needed(item)
            for index in range(item.childCount()):
                child = item.child(index)
                if self._is_lazy_folder_placeholder(child):
                    continue
                walk(child)

        if folders:
            for folder in folders:
                item = self._find_tree_item_for_path(folder)
                if item is not None:
                    walk(item)
            return
        for index in range(self._folder_tree.topLevelItemCount()):
            walk(self._folder_tree.topLevelItem(index))

    def _persist_workspace(self) -> None:
        roots = self._collect_root_folders()
        if roots and self._project.input_dir is None:
            self._project.input_dir = roots[0]
        save_workspace(
            self._project,
            root_folders=roots,
            last_folder=self._current_folder,
            expanded_folders=self._collect_expanded_folders(),
        )

    def _load_initial_roots(self, reset: bool = False) -> None:
        if reset:
            self._folder_tree.clear()
        roots = load_root_folders()
        if not roots and self._project.input_dir and self._project.input_dir.is_dir():
            roots = [self._project.input_dir.resolve()]
        saved_last = load_last_folder()
        saved_expanded = load_expanded_folders()
        for folder in roots:
            self._add_root_folder(
                folder,
                unhide=True,
                persist=False,
                select=False,
                expand_root=False,
            )
        self._apply_expanded_folders(saved_expanded)
        self._sync_folder_tree_placeholders()
        if saved_last is not None:
            self._saved_last_folder = saved_last
        elif roots:
            self._saved_last_folder = roots[0].resolve()

    def _strip_lazy_folder_placeholders(self) -> None:
        def walk(item: QTreeWidgetItem) -> None:
            for index in range(item.childCount() - 1, -1, -1):
                child = item.child(index)
                if self._is_lazy_folder_placeholder(child):
                    item.removeChild(child)
                else:
                    walk(child)

        for index in range(self._folder_tree.topLevelItemCount()):
            walk(self._folder_tree.topLevelItem(index))

    def _add_lazy_folder_placeholder(self, item: QTreeWidgetItem) -> None:
        if item.isExpanded():
            return
        if item.childCount() > 0:
            if item.childCount() == 1 and self._is_lazy_folder_placeholder(item.child(0)):
                return
            return
        placeholder = QTreeWidgetItem([""])
        placeholder.setData(0, ROLE_LAZY_PLACEHOLDER, True)
        placeholder.setFlags(Qt.ItemFlag.NoItemFlags)
        item.addChild(placeholder)

    def _sync_folder_tree_placeholders(self) -> None:
        def walk(item: QTreeWidgetItem) -> None:
            raw = item.data(0, ROLE_PATH)
            if raw:
                folder = Path(str(raw)).expanduser().resolve()
                if folder.is_dir():
                    self._set_folder_expandable(item, folder)
            for index in range(item.childCount()):
                child = item.child(index)
                if self._is_lazy_folder_placeholder(child):
                    continue
                walk(child)

        for index in range(self._folder_tree.topLevelItemCount()):
            walk(self._folder_tree.topLevelItem(index))
        self._sync_folder_tree_icons()
        self._sync_folder_tree_labels()

    def _apply_startup_folder(self) -> None:
        folder = self._saved_last_folder
        self._saved_last_folder = None
        if folder is not None and folder.is_dir():
            self._navigate_to_folder(folder)
            self._persist_workspace()

    def _is_root_hidden(self, folder: Path) -> bool:
        return str(folder.expanduser().resolve()) in self._hidden_root_folders

    def _pick_and_add_folder(self) -> None:
        start = str(Path.home())
        if self._current_folder and self._current_folder.is_dir():
            start = str(self._current_folder)
        picked = QFileDialog.getExistingDirectory(self, "Select folder", start)
        if picked:
            self._add_root_folder(Path(picked), unhide=True, warn_spaces=True)

    def _add_root_folder(
        self,
        folder: Path,
        *,
        unhide: bool = False,
        persist: bool = True,
        select: bool = True,
        expand_root: bool = True,
        warn_spaces: bool = False,
    ) -> None:
        folder = folder.expanduser().resolve()
        if not folder.is_dir():
            return
        if unhide:
            self._hidden_root_folders = unhide_root_folder(self._hidden_root_folders, folder)
        elif self._is_root_hidden(folder):
            return
        already_present = False
        for i in range(self._folder_tree.topLevelItemCount()):
            item = self._folder_tree.topLevelItem(i)
            if Path(str(item.data(0, ROLE_PATH))) == folder:
                already_present = True
                if select:
                    self._folder_tree.select_folder_item(item)
                    self._on_folder_selected(folder)
                break
        if already_present:
            return
        counts = self._folder_panel_counts(folder)
        root = QTreeWidgetItem([self._folder_tree_label(folder, counts=counts)])
        root.setData(0, ROLE_PATH, str(folder))
        n_dirs, n_files = counts
        root.setToolTip(0, f"{folder}\n{n_dirs} folders · {n_files} files")
        apply_folder_item_color(root, folder)
        self._folder_tree.addTopLevelItem(root)
        self._set_folder_expandable(root, folder, counts=counts)
        if expand_root and n_dirs > 0:
            root.setExpanded(True)
            self._populate_folder_tree_children(root)
        self._sync_folder_item_icon(root)
        if select:
            self._folder_tree.select_folder_item(root)
            self._on_folder_selected(folder)
        if self._folder_tree.topLevelItemCount() == 1 and self._project.input_dir is None:
            self._project.input_dir = folder
        if persist:
            self._persist_workspace()
        self._apply_folder_filter()
        if warn_spaces:
            self._warn_folders_with_spaces([folder])

    def _sync_folder_item_icon(self, item: QTreeWidgetItem) -> None:
        if self._is_lazy_folder_placeholder(item):
            return
        # Always set a same-size icon so folder names align by tree level,
        # not by whether an expand arrow is present.
        if (
            item.childIndicatorPolicy()
            == QTreeWidgetItem.ChildIndicatorPolicy.ShowIndicator
        ):
            item.setIcon(
                0,
                browser_folder_icon(self._folder_tree, open_=item.isExpanded()),
            )
        else:
            item.setIcon(0, browser_folder_spacer_icon(self._folder_tree))

    def _sync_folder_tree_icons(self) -> None:
        def walk(item: QTreeWidgetItem) -> None:
            self._sync_folder_item_icon(item)
            for index in range(item.childCount()):
                child = item.child(index)
                if self._is_lazy_folder_placeholder(child):
                    continue
                walk(child)

        for index in range(self._folder_tree.topLevelItemCount()):
            walk(self._folder_tree.topLevelItem(index))

    def _on_folder_tree_expanded(self, item: QTreeWidgetItem) -> None:
        self._sync_folder_item_icon(item)
        self._populate_folder_tree_children(item)
        self._apply_folder_filter()
        self._schedule_persist_workspace()

    def _on_folder_tree_collapsed(self, item: QTreeWidgetItem) -> None:
        self._sync_folder_item_icon(item)
        self._schedule_persist_workspace()

    @staticmethod
    def _is_lazy_folder_placeholder(item: QTreeWidgetItem) -> bool:
        return bool(item.data(0, ROLE_LAZY_PLACEHOLDER))

    @staticmethod
    def _is_intermediate_teaid_folder(name: str) -> bool:
        return name.endswith("_self_blast")

    def _list_subfolders(self, folder: Path) -> list[Path]:
        try:
            with os.scandir(folder) as entries:
                children = [
                    Path(entry.path)
                    for entry in entries
                    if not entry.name.startswith(".")
                    and not self._is_intermediate_teaid_folder(entry.name)
                    and entry.is_dir(follow_symlinks=False)
                ]
            return sorted(
                children,
                key=lambda p: [
                    int(part) if part.isdigit() else part.lower()
                    for part in re.split(r"(\d+)", p.name)
                ],
            )
        except OSError:
            return []

    @staticmethod
    def _folder_has_subfolders(folder: Path) -> bool:
        try:
            with os.scandir(folder) as entries:
                for entry in entries:
                    if (
                        not entry.name.startswith(".")
                        and not MainWindow._is_intermediate_teaid_folder(entry.name)
                        and entry.is_dir(follow_symlinks=False)
                    ):
                        return True
        except OSError:
            return False
        return False

    def _files_for_folder(self, folder: Path) -> list[Path]:
        files: list[Path] = []
        try:
            with os.scandir(folder) as entries:
                for entry in entries:
                    if entry.name.startswith("."):
                        continue
                    try:
                        if not entry.is_file(follow_symlinks=False):
                            continue
                    except OSError:
                        continue
                    path = Path(entry.path)
                    if (
                        is_repeatviewer_hits_cache_file(path)
                        or is_sequence_map_hidden_file(path)
                        or is_flexidot_extra_plot_file(path)
                    ):
                        continue
                    files.append(path)
        except OSError:
            return []
        return files

    def _folder_panel_counts(self, folder: Path) -> tuple[int, int]:
        """Return ``(subdir_count, file_count)`` with one directory scan."""
        n_dirs = 0
        n_files = 0
        try:
            with os.scandir(folder) as entries:
                for entry in entries:
                    name = entry.name
                    if name.startswith("."):
                        continue
                    try:
                        is_dir = entry.is_dir(follow_symlinks=False)
                    except OSError:
                        continue
                    if is_dir:
                        if not self._is_intermediate_teaid_folder(name):
                            n_dirs += 1
                        continue
                    try:
                        if not entry.is_file(follow_symlinks=False):
                            continue
                    except OSError:
                        continue
                    path = Path(entry.path)
                    if (
                        is_repeatviewer_hits_cache_file(path)
                        or is_sequence_map_hidden_file(path)
                        or is_flexidot_extra_plot_file(path)
                    ):
                        continue
                    n_files += 1
        except OSError:
            return 0, 0
        return n_dirs, n_files

    def _folder_tree_label(self, folder: Path, *, counts: tuple[int, int] | None = None) -> str:
        folder = folder.expanduser().resolve()
        n_dirs, n_files = counts if counts is not None else self._folder_panel_counts(folder)
        return f"{folder.name} ({n_dirs}/{n_files})"

    def _set_folder_tree_item_label(self, item: QTreeWidgetItem, folder: Path) -> None:
        folder = folder.expanduser().resolve()
        counts = self._folder_panel_counts(folder)
        item.setText(0, self._folder_tree_label(folder, counts=counts))
        n_dirs, n_files = counts
        item.setToolTip(0, f"{folder}\n{n_dirs} folders · {n_files} files")
        self._set_folder_expandable(item, folder, counts=counts)
    def _sync_folder_tree_labels(self) -> None:
        def walk(item: QTreeWidgetItem) -> None:
            if self._is_lazy_folder_placeholder(item):
                return
            raw = item.data(0, ROLE_PATH)
            if raw:
                folder = Path(str(raw)).expanduser().resolve()
                if folder.is_dir():
                    self._set_folder_tree_item_label(item, folder)
            for index in range(item.childCount()):
                walk(item.child(index))

        for index in range(self._folder_tree.topLevelItemCount()):
            walk(self._folder_tree.topLevelItem(index))
        self._apply_folder_filter()

    def _populate_folder_tree_children(self, item: QTreeWidgetItem) -> None:
        if item.childCount() > 0:
            if item.childCount() == 1 and self._is_lazy_folder_placeholder(item.child(0)):
                item.removeChild(item.child(0))
            else:
                return
        raw = item.data(0, ROLE_PATH)
        if not raw:
            return
        folder = Path(str(raw)).expanduser().resolve()
        if not folder.is_dir():
            return
        self._populate_subfolders(folder, item)

    def _set_folder_expandable(
        self,
        item: QTreeWidgetItem,
        folder: Path,
        *,
        counts: tuple[int, int] | None = None,
    ) -> None:
        n_dirs = counts[0] if counts is not None else None
        has_subfolders = (
            n_dirs > 0 if n_dirs is not None else self._folder_has_subfolders(folder)
        )
        if has_subfolders:
            item.setChildIndicatorPolicy(
                QTreeWidgetItem.ChildIndicatorPolicy.ShowIndicator
            )
            if item.isExpanded():
                if item.childCount() == 0 or (
                    item.childCount() == 1
                    and self._is_lazy_folder_placeholder(item.child(0))
                ):
                    self._populate_folder_tree_children(item)
            else:
                self._add_lazy_folder_placeholder(item)
        else:
            if item.childCount() == 1 and self._is_lazy_folder_placeholder(item.child(0)):
                item.removeChild(item.child(0))
            item.setChildIndicatorPolicy(
                QTreeWidgetItem.ChildIndicatorPolicy.DontShowIndicatorWhenChildless
            )
        self._sync_folder_item_icon(item)

    def _populate_subfolders(self, folder: Path, parent: QTreeWidgetItem) -> None:
        for child in self._list_subfolders(folder):
            counts = self._folder_panel_counts(child)
            item = QTreeWidgetItem([self._folder_tree_label(child, counts=counts)])
            item.setData(0, ROLE_PATH, str(child))
            n_dirs, n_files = counts
            item.setToolTip(0, f"{child}\n{n_dirs} folders · {n_files} files")
            apply_folder_item_color(item, child)
            parent.addChild(item)
            self._set_folder_expandable(item, child, counts=counts)
        self._sync_folder_tree_placeholders()

    def _ensure_folder_tree_path(self, folder: Path) -> QTreeWidgetItem | None:
        folder = folder.expanduser().resolve()
        found = self._find_tree_item_for_path(folder)
        if found is not None:
            return found
        for index in range(self._folder_tree.topLevelItemCount()):
            root_item = self._folder_tree.topLevelItem(index)
            raw = root_item.data(0, ROLE_PATH)
            if not raw:
                continue
            root_path = Path(str(raw)).expanduser().resolve()
            try:
                rel_parts = folder.relative_to(root_path).parts
            except ValueError:
                continue
            current_item = root_item
            current_path = root_path
            for part in rel_parts:
                self._populate_folder_tree_children(current_item)
                current_path = current_path / part
                child = self._find_tree_item_for_path(current_path)
                if child is None:
                    break
                current_item = child
            else:
                return self._find_tree_item_for_path(folder)
        return None

    def _find_tree_item_for_path(self, folder: Path) -> QTreeWidgetItem | None:
        folder = folder.expanduser().resolve()

        def walk(item: QTreeWidgetItem) -> QTreeWidgetItem | None:
            if self._is_lazy_folder_placeholder(item):
                return None
            raw = item.data(0, ROLE_PATH)
            if not raw:
                return None
            if Path(str(raw)).expanduser().resolve() == folder:
                return item
            for index in range(item.childCount()):
                found = walk(item.child(index))
                if found is not None:
                    return found
            return None

        for index in range(self._folder_tree.topLevelItemCount()):
            found = walk(self._folder_tree.topLevelItem(index))
            if found is not None:
                return found
        return None

    def _refresh_subfolders_under(
        self,
        parent_folder: Path,
        *,
        keep_folder: Path | None = None,
    ) -> QTreeWidgetItem | None:
        parent_folder = parent_folder.expanduser().resolve()
        item = self._find_tree_item_for_path(parent_folder)
        if item is None:
            return None
        if keep_folder is None:
            keep_folder = (
                self._current_folder.expanduser().resolve()
                if self._current_folder is not None
                else None
            )
        else:
            keep_folder = keep_folder.expanduser().resolve()
        # Removing the selected child makes Qt select the next sibling and fire
        # folder_selected — suppress that so we do not jump to the folder below.
        self._folder_tree._suppress_current_changed = True
        try:
            while item.childCount():
                item.removeChild(item.child(0))
            self._populate_subfolders(parent_folder, item)
            item.setExpanded(True)
            self._sync_folder_tree_placeholders()
            if keep_folder is not None:
                keep_item = self._find_tree_item_for_path(keep_folder)
                if keep_item is not None:
                    self._folder_tree.selectionModel().select(
                        self._folder_tree.indexFromItem(keep_item),
                        QItemSelectionModel.SelectionFlag.ClearAndSelect
                        | QItemSelectionModel.SelectionFlag.Current,
                    )
                    self._folder_tree.scrollToItem(keep_item)
        finally:
            QTimer.singleShot(
                0,
                lambda: setattr(self._folder_tree, "_suppress_current_changed", False),
            )
        return item

    def _ensure_folder_visible_in_tree(self, folder: Path) -> QTreeWidgetItem | None:
        return self._ensure_folder_tree_path(folder)

    def _workspace_root_folders(self) -> list[Path]:
        roots: list[Path] = []
        for index in range(self._folder_tree.topLevelItemCount()):
            item = self._folder_tree.topLevelItem(index)
            raw = item.data(0, ROLE_PATH)
            if not raw:
                continue
            path = Path(str(raw)).expanduser().resolve()
            if path.is_dir():
                roots.append(path)
        return roots

    def _on_folder_nav_requested(self, folder: Path, record_history: bool) -> None:
        self._navigate_to_folder(folder, record_history=record_history)

    def _navigate_to_folder(
        self,
        folder: Path,
        *,
        new_outputs: list[Path] | None = None,
        source_files: list[Path] | None = None,
        record_history: bool = True,
        open_preview: bool = True,
    ) -> None:
        folder = folder.expanduser().resolve()
        self._ensure_folder_visible_in_tree(folder)
        item = self._find_tree_item_for_path(folder)
        if item is not None:
            parent = item.parent()
            while parent is not None:
                parent.setExpanded(True)
                parent = parent.parent()
            self._folder_tree.select_folder_item(item)
            self._folder_tree.scrollToItem(item)
        self._populate_file_list(
            folder,
            new_outputs=new_outputs,
            source_files=source_files,
            record_history=record_history,
            open_preview=open_preview,
        )

    def _on_folder_selected(self, folder: Path) -> None:
        self._sync_define_database_dialog_from_folder(folder)
        self._sync_te_viz_dialog_from_folder(folder)
        self._populate_file_list(folder, record_history=True)
        self._update_directory_watchers()
        self._persist_workspace()
        self._folder_tree.setFocus(Qt.FocusReason.OtherFocusReason)

    def _refresh_current_folder(self) -> None:
        if self._current_folder and self._current_folder.is_dir():
            self._populate_file_list(self._current_folder)
        else:
            self._set_status_message("Select a folder first", 4000)

    def _on_preview_annotations_saved(self, path: object) -> None:
        # Update/highlight the written GFF in Files, but keep the current FASTA
        # selection so the Annotation preview tab stays open after edits/deletes.
        self._sync_sidecar_in_file_list(path, select=False)

    def _on_preview_sequence_extracted(self, path: object) -> None:
        # Show the new FASTA in Files, but keep the current file selected/previewed.
        self._sync_sidecar_in_file_list(path, select=False)

    def _on_preview_file_written(self, path: object) -> None:
        # Exported tree figure etc.: show it in Files, keep the current preview.
        self._sync_sidecar_in_file_list(path, select=False)

    def _sync_sidecar_in_file_list(self, path: object, *, select: bool) -> None:
        if not isinstance(path, Path) or not self._current_folder:
            return
        resolved = path.expanduser().resolve()
        if resolved.parent != self._current_folder:
            return

        sidecar_key = str(resolved)
        prev_key: str | None = None
        current = self._file_list.currentItem()
        if current is not None:
            prev_path = self._preview_path_from_item(current)
            if prev_path is not None:
                prev_key = str(prev_path.expanduser().resolve())

        def _find_item_for_key(key: str) -> FileListRow | None:
            for index in range(self._file_list.count()):
                item = self._file_list.item(index)
                raw = item.data(ROLE_PATH)
                if raw and str(Path(str(raw)).expanduser().resolve()) == key:
                    return item
            return None

        self._file_list.blockSignals(True)
        try:
            # Keep the green “new output” tint even when we do not switch selection.
            self._highlight_output_keys.add(sidecar_key)

            sidecar_item = _find_item_for_key(sidecar_key)
            if sidecar_item is None:
                self._populate_file_list(self._current_folder, open_preview=False)
                # _populate_file_list clears blockSignals; re-block for the rest.
                self._file_list.blockSignals(True)
                # Populate without new_outputs preserves existing highlight keys;
                # re-assert this path so metadata refresh keeps the green style.
                self._highlight_output_keys.add(sidecar_key)
                sidecar_item = _find_item_for_key(sidecar_key)

            if sidecar_item is not None:
                self._style_new_output_item(sidecar_item)
                self._file_list.scrollToItem(
                    sidecar_item,
                    QAbstractItemView.ScrollHint.EnsureVisible,
                )

            if select:
                if sidecar_item is not None:
                    self._reveal_file_list_item(
                        sidecar_item,
                        scroll_hint=QAbstractItemView.ScrollHint.EnsureVisible,
                    )
            elif prev_key is not None:
                # Keep the previously selected FASTA / preview after a list refresh.
                prev_item = _find_item_for_key(prev_key)
                if prev_item is not None:
                    prev_item.setSelected(True)
                    self._file_list.setCurrentItem(prev_item)
        finally:
            self._file_list.blockSignals(False)

        if not self._file_list.hasFocus():
            QTimer.singleShot(0, self._focus_preview_widget)

    def _find_fasta_for_dotplot_go_to(self, seq_id: str) -> Path | None:
        found = self._find_fasta_for_sequence_id(seq_id)
        if found is not None:
            return found
        if self._preview_tabs is None:
            return None
        from repeatviewer.widgets.preview_tabs import FilePreviewPage

        current = self._preview_tabs.currentWidget()
        if isinstance(current, FilePreviewPage) and is_flexidot_html_file(current.path):
            return fasta_for_flexidot_html(current.path)
        return None

    def _on_dotplot_region_selected(self, payload: object) -> None:
        """Open (if needed) and jump the Annotation view to the region clicked in a Flexidot dotplot."""
        logger.debug("Go-to: received dotplot region-selected payload: %r", payload)
        if not isinstance(payload, dict):
            return

        targets = parse_dotplot_go_to_targets(payload)
        if not targets:
            logger.warning("Dotplot Go-to payload has no usable axis: %r", payload)
            return

        try:
            for seq, start, end in targets:
                if self._preview_tabs is None:
                    continue
                if self._preview_tabs.show_annotation_region_for_sequence(seq, start, end):
                    logger.debug("Go-to: jumped on open tab seq=%r %d-%d", seq, start, end)
                    self._ensure_preview_visible()
                    self._set_status_message(f"Jumped Annotation view to {seq}:{start}-{end}", 4000)
                    continue

                fasta_path = self._find_fasta_for_dotplot_go_to(seq)
                logger.debug("Go-to: seq=%r - folder search found: %s", seq, fasta_path)
                if fasta_path is None:
                    self._set_status_message(
                        f"Could not find a FASTA file for {seq} in this folder", 4000
                    )
                    continue

                self._open_preview_tab(fasta_path)

                def _jump_after_open(
                    *,
                    seq_id: str = seq,
                    region_start: int = start,
                    region_end: int = end,
                    fasta: Path = fasta_path,
                    attempt: int = 0,
                ) -> None:
                    if self._preview_tabs is None:
                        return
                    jumped = self._preview_tabs.show_annotation_region_for_sequence(
                        seq_id, region_start, region_end
                    )
                    if jumped:
                        logger.debug(
                            "Go-to: jump after opening %s seq=%r %d-%d",
                            fasta.name,
                            seq_id,
                            region_start,
                            region_end,
                        )
                        self._ensure_preview_visible()
                        self._set_status_message(
                            f"Jumped Annotation view to {seq_id}:{region_start}-{region_end}",
                            4000,
                        )
                        return
                    if attempt < 10:
                        QTimer.singleShot(50, lambda: _jump_after_open(attempt=attempt + 1))
                        return
                    self._set_status_message(
                        f"Opened {fasta.name} but could not jump to {seq_id}:{region_start}-{region_end}",
                        6000,
                    )

                QTimer.singleShot(0, _jump_after_open)
        except Exception:  # noqa: BLE001
            logger.exception("Failed to jump Annotation view for dotplot selection %r", payload)
            self._set_status_message(
                "Go-to-region failed - see terminal/log output for details", 6000
            )

    def _find_fasta_for_sequence_id(self, seq_id: str) -> Path | None:
        """Find a FASTA in the current folder that contains ``seq_id``."""
        from repeatviewer.backends.gff_features import _seqid_matches
        from repeatviewer.backends.sequence_features import list_fasta_record_ids

        if not self._current_folder:
            return None
        candidates = self._cached_folder_files or [
            p for p in self._current_folder.iterdir() if p.is_file()
        ]
        wanted = str(seq_id or "").strip()
        for path in candidates:
            if path.suffix.lower() not in _FASTA_SUFFIXES:
                continue
            try:
                record_ids = list_fasta_record_ids(path)
            except (OSError, ValueError):
                continue
            for record_id in record_ids:
                if (
                    record_id == wanted
                    or flexidot_sanitize_name(record_id) == wanted
                    or _seqid_matches(record_id, wanted)
                ):
                    return path
        return None

    def _focus_preview_widget(self) -> None:
        if self._preview_tabs is None:
            return
        current = self._preview_tabs.currentWidget()
        if current is not None:
            current.setFocus(Qt.FocusReason.OtherFocusReason)

    def _current_file_sort_mode(self) -> FileSortMode:
        return load_file_sort_mode()

    def _apply_file_filter(self) -> None:
        if not self._current_folder:
            return
        from repeatviewer.widgets.file_list_columns import metadata_columns_visible

        query = self._file_filter.text().strip().lower()
        self._file_list.clear_rows()
        visible: list[Path] = []
        for file_path in self._cached_folder_files:
            if query and query not in file_path.name.lower():
                continue
            visible.append(file_path)

        visible = sort_file_paths(visible, self._current_file_sort_mode())

        want_metadata = metadata_columns_visible()
        # Large folders: skip sync stat/metadata on the GUI thread; fill later.
        fill_metadata = want_metadata and len(visible) <= _LARGE_FILE_LIST
        defer_metadata = want_metadata and not fill_metadata
        was_sorting = self._file_list.isSortingEnabled()
        self._file_list.setSortingEnabled(False)
        self._file_list.setUpdatesEnabled(False)
        try:
            need_style = bool(self._highlight_output_keys or self._highlight_source_keys)
            for file_path in visible:
                row_item = self._file_list.add_file_row(
                    file_path,
                    fill_metadata_columns=fill_metadata,
                )
                if need_style:
                    key = str(file_path.expanduser().resolve())
                    if key in self._highlight_output_keys:
                        self._style_new_output_item(row_item)
                    elif key in self._highlight_source_keys:
                        self._style_source_file_item(row_item)
        finally:
            self._file_list.setUpdatesEnabled(True)
            self._file_list.setSortingEnabled(was_sorting)

        if not visible and self._cached_folder_files:
            self._file_list.add_placeholder_row("(no files match filter)")
        elif not self._cached_folder_files:
            self._file_list.add_placeholder_row("(empty folder)")

        self._file_list.apply_default_sort()
        self._file_list.repair_chrome()

        self._sync_files_status(
            len(visible),
            len(self._cached_folder_files),
        )
        if fill_metadata and not self._any_background_job_running():
            self._start_metadata_enrichment(visible)
        elif defer_metadata and not self._any_background_job_running():
            QTimer.singleShot(0, self._fill_visible_file_metadata)

    def _stop_metadata_enrichment(self, *, wait_ms: int = 1000) -> None:
        self._pending_metadata_ui.clear()
        self._metadata_stat_queue.clear()
        self._metadata_ui_timer.stop()
        worker = self._metadata_enrich_worker
        thread = self._metadata_enrich_thread
        self._metadata_enrich_worker = None
        self._metadata_enrich_thread = None
        if worker is not None:
            worker.cancel()
        if thread is not None:
            thread.quit()
            # Job start must not block the GUI for up to 1s waiting on enrichment.
            if wait_ms > 0:
                thread.wait(wait_ms)

    def _start_metadata_enrichment(self, paths: list[Path]) -> None:
        from repeatviewer.browser_file_metadata import is_fasta_path
        from repeatviewer.widgets.file_list_columns import metadata_columns_visible

        if not metadata_columns_visible():
            return
        self._stop_metadata_enrichment()
        fasta_paths = [p.expanduser().resolve() for p in paths if is_fasta_path(p)]
        if not fasta_paths:
            return
        self._metadata_enrich_generation += 1
        generation = self._metadata_enrich_generation
        thread = QThread(self)
        worker = _MetadataEnrichWorker(fasta_paths, generation)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        # Connect to the bound method, never a lambda: Qt needs a receiver
        # QObject to see that this crosses threads and queue the call. The
        # QueuedConnection is stated explicitly so this cannot regress.
        worker.batch_ready.connect(
            self._on_metadata_enrich_batch,
            Qt.ConnectionType.QueuedConnection,
        )
        worker.finished.connect(thread.quit)
        worker.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._on_metadata_enrich_thread_finished)
        self._metadata_enrich_thread = thread
        self._metadata_enrich_worker = worker
        thread.start()

    def _on_metadata_enrich_thread_finished(self) -> None:
        self._metadata_enrich_thread = None
        self._metadata_enrich_worker = None

    @Slot(list, int)
    def _on_metadata_enrich_batch(self, batch: list, generation: int) -> None:
        """Runs on the GUI thread: the connection is explicitly queued."""
        if generation != self._metadata_enrich_generation:
            return
        if self._file_list is None:
            return
        updates = [(Path(str(path_text)), meta) for path_text, meta in batch]
        self._enqueue_metadata_ui_updates(updates)

    def _enqueue_metadata_ui_updates(
        self,
        updates: list[tuple[Path, object]],
    ) -> None:
        if not updates or self._file_list is None:
            return
        self._pending_metadata_ui.extend(updates)
        self._schedule_metadata_ui_drain()

    def _drain_metadata_ui_queue(self) -> None:
        if self._file_list is None:
            self._pending_metadata_ui.clear()
            self._metadata_stat_queue.clear()
            return
        if self._metadata_stat_queue:
            self._pending_metadata_ui.extend(self._collect_metadata_stat_chunk())
        chunk = self._pending_metadata_ui[:_METADATA_UI_CHUNK]
        del self._pending_metadata_ui[:_METADATA_UI_CHUNK]
        if chunk:
            self._file_list.update_rows_metadata(chunk)
            self._refresh_metadata_row_highlights(chunk)
        if self._metadata_stat_queue or self._pending_metadata_ui:
            self._metadata_ui_timer.start(0)
        elif self._highlight_output_keys or self._highlight_source_keys:
            # Chunks already re-tint touched rows; only restyle when highlights
            # exist, and prefer the keyed path for large folders.
            self._refresh_file_list_highlights()

    def _refresh_metadata_row_highlights(
        self,
        updates: list[tuple[Path, object]],
    ) -> None:
        """Re-tint only rows touched by a metadata chunk (cheap during scroll)."""
        if self._file_list is None or not updates:
            return
        row_index = self._file_list._path_row_index()
        for path, _meta in updates:
            key = str(path.expanduser().resolve())
            row = row_index.get(key)
            if row is None:
                continue
            item = self._file_list.item(row)
            if item is None:
                continue
            if key in self._highlight_output_keys:
                self._style_new_output_item(item)
            elif key in self._highlight_source_keys:
                self._style_source_file_item(item)

    def _sync_files_status(
        self,
        visible_count: int,
        total_count: int,
    ) -> None:
        if self._job_running or not self._current_folder:
            return
        # File link owns the left status slot — do not paint folder text over it.
        if self._status_linked_file is not None:
            self.statusBar().clearMessage()
            self._status_file_link.show()
            return
        self._status_file_link.hide()
        if visible_count == total_count:
            count = f"{total_count} files"
        else:
            count = f"{visible_count} of {total_count} files"
        genome = (
            "genome set"
            if self._project.genome_file and self._project.genome_file.is_file()
            else "genome not set"
        )
        self.statusBar().showMessage(
            f"{self._current_folder.name} · {count} · {genome}",
            0,
        )

    def _populate_file_list(
        self,
        folder: Path,
        *,
        new_outputs: list[Path] | None = None,
        source_files: list[Path] | None = None,
        open_preview: bool = True,
        record_history: bool = True,
    ) -> None:
        target_folder = folder.expanduser().resolve()
        prev_folder = self._current_folder
        same_folder = False
        if prev_folder is not None:
            try:
                same_folder = prev_folder.expanduser().resolve() == target_folder
            except OSError:
                same_folder = False

        highlight_request = bool(new_outputs or source_files)
        large_folder = (
            same_folder
            and self._cached_folder_files is not None
            and len(self._cached_folder_files) > _LARGE_FILE_LIST
        )
        # Large folder + few known new outputs: skip full scandir/rebuild.
        # Membership is by basename in the current folder (unique on disk).
        few_new_outputs = len(new_outputs or []) <= _LARGE_JOB_OUTPUT_BATCH
        if (
            large_folder
            and highlight_request
            and few_new_outputs
            and self._file_list.rowCount() > 0
        ):
            self._current_folder = target_folder
            self._folder_tree.set_active_folder(target_folder)
            prev_highlight_keys = (
                self._highlight_output_keys | self._highlight_source_keys
            )
            if prev_highlight_keys:
                self._file_list.reapply_file_colors_for_keys(prev_highlight_keys)
            if new_outputs or source_files:
                self._highlight_output_keys = set()
                self._highlight_source_keys = set()
            highlight_paths: list[Path] = []

            def _track_fast(path: Path, *, as_output: bool) -> None:
                resolved = path.expanduser().resolve()
                if resolved.parent != self._current_folder:
                    return
                key = str(resolved)
                if as_output:
                    self._highlight_output_keys.add(key)
                    self._highlight_source_keys.discard(key)
                elif key not in self._highlight_output_keys:
                    self._highlight_source_keys.add(key)
                highlight_paths.append(resolved)

            if new_outputs:
                for path in new_outputs:
                    _track_fast(path, as_output=True)
            if source_files:
                for path in source_files:
                    _track_fast(path, as_output=False)

            existing_names = {p.name for p in self._cached_folder_files}
            added_paths: list[Path] = []
            for path in new_outputs or []:
                resolved = path.expanduser().resolve()
                if resolved.parent != self._current_folder:
                    continue
                if resolved.name in existing_names:
                    continue
                if not resolved.is_file() and not resolved.is_dir():
                    continue
                # Same absolute form as scandir (folder is already resolved).
                cached = self._current_folder / resolved.name
                self._cached_folder_files.append(cached)
                existing_names.add(cached.name)
                added_paths.append(cached)

            self._file_list.blockSignals(True)
            if added_paths:
                self._incremental_add_files_to_list(added_paths)
            else:
                # Outputs already in the list (rerun) — tint only those rows.
                self._apply_job_highlights_to_keys(
                    self._highlight_output_keys | self._highlight_source_keys
                )
            if highlight_paths:
                self._select_and_scroll_outputs(highlight_paths, open_preview=False)
                n_out = len(self._highlight_output_keys)
                n_src = len(self._highlight_source_keys)
                names = ", ".join(
                    Path(k).name for k in list(self._highlight_output_keys)[:3]
                )
                msg = f"Highlighted {n_out} new output(s) and {n_src} input(s)"
                if names:
                    msg += f" (e.g. {names})"
                self._set_status_message(msg, 12000, override_file_link=True)
            self._file_list.blockSignals(False)
            if highlight_paths and open_preview:
                self._open_best_preview_path(highlight_paths)
            if hasattr(self, "_folder_nav"):
                self._folder_nav.set_folder(
                    self._current_folder,
                    record_history=record_history,
                )
            self._sync_folder_tree_label_for_folder(self._current_folder)
            return

        prev_cached_keys = (
            {str(path) for path in self._cached_folder_files}
            if same_folder
            else set()
        )
        self._current_folder = target_folder
        self._folder_tree.set_active_folder(target_folder)
        # File name search is per-folder: clear it when navigating away so the
        # previous query does not hide files in the next directory.
        if not same_folder and self._file_filter.text():
            self._file_filter.blockSignals(True)
            self._file_filter.clear()
            self._file_filter.blockSignals(False)
        try:
            self._cached_folder_files = self._files_for_folder(self._current_folder)
        except OSError as exc:
            self._cached_folder_files = []
            self._file_list.clear_rows()
            self._set_status_message(f"Cannot read folder: {exc}", 6000)
            return

        # Keep highlight tint when the folder watcher refreshes the same directory
        # shortly after a job finishes; only reset when navigating away or when
        # new highlight targets are supplied explicitly.
        if new_outputs or source_files:
            self._highlight_output_keys = set()
            self._highlight_source_keys = set()
        elif not same_folder:
            self._highlight_output_keys = set()
            self._highlight_source_keys = set()
        highlight_paths = []

        def _track(path: Path, *, as_output: bool) -> None:
            resolved = path.expanduser().resolve()
            if resolved.parent != self._current_folder:
                return
            key = str(resolved)
            if as_output:
                self._highlight_output_keys.add(key)
                self._highlight_source_keys.discard(key)
            else:
                # Prefer output tint when a path is both source and result.
                if key not in self._highlight_output_keys:
                    self._highlight_source_keys.add(key)
            highlight_paths.append(resolved)

        if new_outputs:
            for path in new_outputs:
                _track(path, as_output=True)
        if source_files:
            for path in source_files:
                _track(path, as_output=False)

        # Apply green/blue tints while adding/rebuilding rows — skip a full-list
        # restyle of the *old* table first (that was O(n) resolve work for nothing).

        # Use unresolved absolute scandir paths for membership (no N× resolve).
        new_cache_keys = {str(path) for path in self._cached_folder_files}
        added_keys = new_cache_keys - prev_cached_keys
        removed_keys = prev_cached_keys - new_cache_keys

        # Same folder, highlights only, no membership change: do not rebuild.
        if (
            same_folder
            and highlight_paths
            and not added_keys
            and not removed_keys
        ):
            self._file_list.blockSignals(True)
            self._apply_job_highlights_to_keys(
                self._highlight_output_keys | self._highlight_source_keys
            )
            self._select_and_scroll_outputs(highlight_paths, open_preview=False)
            n_out = len(self._highlight_output_keys)
            n_src = len(self._highlight_source_keys)
            names = ", ".join(
                Path(k).name for k in list(self._highlight_output_keys)[:3]
            )
            msg = f"Highlighted {n_out} new output(s) and {n_src} input(s)"
            if names:
                msg += f" (e.g. {names})"
            self._set_status_message(msg, 12000, override_file_link=True)
            self._file_list.blockSignals(False)
            if open_preview:
                self._open_best_preview_path(highlight_paths)
            if hasattr(self, "_folder_nav"):
                self._folder_nav.set_folder(
                    self._current_folder,
                    record_history=record_history,
                )
            self._sync_folder_tree_label_for_folder(self._current_folder)
            return

        if (
            same_folder
            and highlight_paths
            and added_keys
            and not removed_keys
        ):
            added_paths = [
                path
                for path in self._cached_folder_files
                if str(path) in added_keys
            ]
            self._file_list.blockSignals(True)
            self._incremental_add_files_to_list(added_paths)
            self._select_and_scroll_outputs(highlight_paths, open_preview=False)
            n_out = len(self._highlight_output_keys)
            n_src = len(self._highlight_source_keys)
            names = ", ".join(
                Path(k).name for k in list(self._highlight_output_keys)[:3]
            )
            msg = f"Highlighted {n_out} new output(s) and {n_src} input(s)"
            if names:
                msg += f" (e.g. {names})"
            self._set_status_message(msg, 12000, override_file_link=True)
            self._file_list.blockSignals(False)
            if open_preview:
                self._open_best_preview_path(highlight_paths)
            if hasattr(self, "_folder_nav"):
                self._folder_nav.set_folder(
                    self._current_folder,
                    record_history=record_history,
                )
            self._sync_folder_tree_label_for_folder(self._current_folder)
            return

        # Full rebuild: restyle of the previous list is pointless; rows are recreated.
        preserve_ui = same_folder
        prev_key: str | None = None
        prev_scroll: int | None = None
        if preserve_ui:
            current = self._file_list.currentItem()
            if current is not None:
                prev_path = self._preview_path_from_item(current)
                if prev_path is not None:
                    prev_key = str(prev_path.expanduser().resolve())
            prev_scroll = int(self._file_list.verticalScrollBar().value())

        item = self._ensure_folder_visible_in_tree(target_folder)
        if item is not None and not same_folder:
            parent = item.parent()
            while parent is not None:
                if not parent.isExpanded():
                    parent.setExpanded(True)
                    self._populate_folder_tree_children(parent)
                parent = parent.parent()
            self._folder_tree.select_folder_item(item)
            self._folder_tree.scrollToItem(item)
        elif item is not None and same_folder:
            self._folder_tree.select_folder_item(item)

        self._file_list.blockSignals(True)
        self._apply_file_filter()

        if highlight_paths:
            self._select_and_scroll_outputs(highlight_paths, open_preview=False)
            n_out = len(self._highlight_output_keys)
            n_src = len(self._highlight_source_keys)
            names = ", ".join(
                Path(k).name for k in list(self._highlight_output_keys)[:3]
            )
            msg = f"Highlighted {n_out} new output(s) and {n_src} input(s)"
            if names:
                msg += f" (e.g. {names})"
            self._set_status_message(msg, 12000, override_file_link=True)
        else:
            self._sync_files_status(
                len(self._cached_folder_files),
                len(self._cached_folder_files),
            )

        if preserve_ui and not highlight_paths:
            if prev_key:
                match: FileListRow | None = None
                for index in range(self._file_list.count()):
                    item = self._file_list.item(index)
                    raw = item.data(ROLE_PATH)
                    if raw and str(raw) == prev_key:
                        match = item
                        break
                if match is not None:
                    match.setSelected(True)
                    self._file_list.setCurrentItem(match)
            if prev_scroll is not None:
                vsb = self._file_list.verticalScrollBar()
                vsb.setValue(min(prev_scroll, vsb.maximum()))
        self._file_list.blockSignals(False)
        if highlight_paths and open_preview:
            self._open_best_preview_path(highlight_paths)
        if hasattr(self, "_folder_nav"):
            self._folder_nav.set_folder(self._current_folder, record_history=record_history)
        if same_folder:
            self._sync_folder_tree_label_for_folder(self._current_folder)
        else:
            self._sync_folder_tree_labels()
        self._folder_tree.set_active_folder(self._current_folder)

    def _clear_job_file_highlights(self) -> None:
        """Drop previous job green/blue tints (called when a new job starts)."""
        prev_keys = self._highlight_output_keys | self._highlight_source_keys
        self._highlight_output_keys = set()
        self._highlight_source_keys = set()
        if not prev_keys:
            return
        # Only clear previously tinted rows — never restyle all 1700 files.
        self._file_list.reapply_file_colors_for_keys(prev_keys)

    def _apply_job_highlights_to_keys(self, keys: set[str]) -> None:
        """Tint only the given path keys (one path-index build)."""
        if not keys or self._file_list is None:
            return
        row_index = self._file_list._path_row_index()
        for key in keys:
            row = row_index.get(key)
            if row is None:
                continue
            item = self._file_list.item(row)
            if item is None:
                continue
            if key in self._highlight_output_keys:
                self._style_new_output_item(item)
            elif key in self._highlight_source_keys:
                self._style_source_file_item(item)

    def _refresh_file_list_highlights(self) -> None:
        """Re-apply new-output / source tint after metadata or sort updates."""
        n_rows = self._file_list.count()
        out_keys = self._highlight_output_keys
        src_keys = self._highlight_source_keys
        if not out_keys and not src_keys:
            # Nothing to tint. Callers that need a full clear should use
            # ``_clear_job_file_highlights`` (keyed) instead of walking all rows.
            return
        # Large folders: only touch highlighted rows (O(highlights), not O(folder)).
        if n_rows > _LARGE_FILE_LIST and len(out_keys) + len(src_keys) <= _LARGE_JOB_OUTPUT_BATCH:
            self._apply_job_highlights_to_keys(out_keys | src_keys)
            return
        for index in range(n_rows):
            item = self._file_list.item(index)
            raw = item.data(ROLE_PATH)
            if not raw:
                continue
            # ROLE_PATH is stored resolved; avoid another resolve() per row.
            key = str(raw)
            if key in out_keys:
                self._style_new_output_item(item)
            elif key in src_keys:
                self._style_source_file_item(item)
            else:
                self._file_list.reapply_file_color_at_row(index)

    def _palette_is_dark(self) -> bool:
        from PySide6.QtGui import QPalette

        return self.palette().color(QPalette.ColorRole.Window).lightness() < 128

    def _style_new_output_item(self, item: FileListRow) -> None:
        if self._palette_is_dark():
            bg, fg = _NEW_OUTPUT_BG_DARK, _NEW_OUTPUT_FG_DARK
        else:
            bg, fg = _NEW_OUTPUT_BG_LIGHT, _NEW_OUTPUT_FG_LIGHT
        item.setBackground(QBrush(bg))
        item.setForeground(QBrush(fg))
        font = item.font()
        font.setBold(True)
        item.setFont(font)
        item.setToolTip("Newly generated file from the last job")

    def _style_source_file_item(self, item: FileListRow) -> None:
        if self._palette_is_dark():
            bg, fg = _SOURCE_BG_DARK, _SOURCE_FG_DARK
        else:
            bg, fg = _SOURCE_BG_LIGHT, _SOURCE_FG_LIGHT
        item.setBackground(QBrush(bg))
        item.setForeground(QBrush(fg))
        font = item.font()
        font.setBold(True)
        item.setFont(font)
        item.setToolTip("Input file used for the last job")

    def _reveal_file_list_item(
        self,
        item: FileListRow,
        *,
        scroll_hint: QAbstractItemView.ScrollHint = QAbstractItemView.ScrollHint.EnsureVisible,
    ) -> None:
        item.setSelected(True)
        self._file_list.setCurrentItem(item)
        self._file_list.scrollToItem(item, scroll_hint)

    def _select_and_scroll_outputs(
        self,
        paths: list[Path],
        *,
        open_preview: bool = True,
    ) -> None:
        """Scroll job results into view; keep green/light-blue tints (no multi-select).

        The opened file is selected separately (deep blue) by ``_auto_open_outputs``.
        """
        path_set = {str(p.expanduser().resolve()) for p in paths}
        first_item: FileListRow | None = None
        preferred_item: FileListRow | None = None
        preferred_key: tuple[int, int] | None = None
        # Do not select every input/output — that would paint them all deep blue.
        self._file_list.clearSelection()
        self._file_list.setCurrentItem(None)
        # Tint only job keys — never restyle the whole large folder here.
        self._apply_job_highlights_to_keys(
            self._highlight_output_keys | self._highlight_source_keys
        )
        large_list = self._file_list.count() > _LARGE_FILE_LIST
        large_batch = len(path_set) > _LARGE_JOB_OUTPUT_BATCH

        def _consider(item: FileListRow, key: str) -> None:
            nonlocal first_item, preferred_item, preferred_key
            if first_item is None:
                first_item = item
            if large_batch:
                if preferred_item is None and key in self._highlight_output_keys:
                    preferred_item = item
                return
            candidate = Path(key)
            rank = (
                preview_path_priority(candidate),
                0 if key in self._highlight_output_keys else 1,
            )
            if preferred_key is None or rank < preferred_key:
                preferred_key = rank
                preferred_item = item

        if large_list and not large_batch:
            # Direct lookup for a few paths instead of scanning every row.
            row_index = self._file_list._path_row_index()
            for key in path_set:
                row = row_index.get(key)
                if row is None:
                    continue
                item = self._file_list.item(row)
                if item is not None:
                    _consider(item, key)
        else:
            for index in range(self._file_list.count()):
                item = self._file_list.item(index)
                raw = item.data(ROLE_PATH)
                if not raw:
                    continue
                key = str(raw)
                if key in path_set:
                    _consider(item, key)

        focus_item = preferred_item or first_item
        if focus_item is not None:
            self._file_list.scrollToItem(
                focus_item,
                QAbstractItemView.ScrollHint.EnsureVisible,
            )
        if open_preview and focus_item is not None:
            path = self._preview_path_from_item(focus_item)
            if path is not None and is_previewable_file(path):
                self._open_file_in_preview(path)
                self._reveal_file_list_item(focus_item)
                self._apply_job_highlights_to_keys(
                    self._highlight_output_keys | self._highlight_source_keys
                )

    def _open_best_preview_path(self, paths: list[Path]) -> None:
        if not paths:
            return
        # Large batches: pick without N× is_file() syscalls.
        if len(paths) > _LARGE_JOB_OUTPUT_BATCH:
            target = self._pick_preview_output(list(paths))
            if target is not None and is_previewable_file(target):
                self._open_file_in_preview(target)
            return
        candidates = [
            p.expanduser().resolve()
            for p in paths
            if p.is_file() and is_previewable_file(p)
        ]
        if not candidates:
            return

        def _rank(path: Path) -> tuple[int, int]:
            # Prefer new outputs over job inputs when priorities tie (e.g. both FASTA).
            is_output = str(path) in self._highlight_output_keys
            return (preview_path_priority(path), 0 if is_output else 1)

        best = min(candidates, key=_rank)
        self._open_file_in_preview(best)

    def _select_path_in_file_list(self, path: Path) -> None:
        """Select one file in the list (used after creating a genome map HTML)."""
        resolved = path.expanduser().resolve()
        key = str(resolved)
        row = self._file_list._path_row_index().get(key)
        if row is None:
            return
        item = self._file_list.item(row)
        if item is None:
            return
        self._reveal_file_list_item(
            item,
            scroll_hint=QAbstractItemView.ScrollHint.EnsureVisible,
        )
        # Re-assert green/blue on job rows only (selection paints deep blue).
        self._apply_job_highlights_to_keys(
            self._highlight_output_keys | self._highlight_source_keys
        )

    def _selected_files(self) -> list[Path]:
        files = []
        items = self._file_list.selectedItems()
        if not items:
            current = self._file_list.currentItem()
            if current is not None:
                items = [current]
        for item in items:
            raw = item.data(ROLE_PATH)
            if not raw:
                continue
            path = Path(str(raw))
            if path.is_file():
                files.append(path)
        return files

    def _preview_path_from_item(self, item: FileListRow | None) -> Path | None:
        if item is None:
            return None
        raw = item.data(ROLE_PATH)
        if not raw:
            return None
        path = Path(str(raw))
        return path if path.is_file() else None

    def _open_file_in_preview(self, path: Path) -> None:
        self._open_preview_tab(path, replace=True)
        self._set_status_file_link(path)

    def _on_file_current_item_changed(
        self,
        current: FileListRow | None,
        _previous: FileListRow | None,
    ) -> None:
        if self._file_list.suppress_selection_preview():
            return
        if self._define_database_dialog is not None and self._define_database_dialog.isVisible():
            return
        if self._te_viz_input_dialog is not None and self._te_viz_input_dialog.isVisible():
            return
        if self._plugin_run_dialog is not None and self._plugin_run_dialog.isVisible():
            return
        if self._tree_trait_dialog_open():
            return  # a click fills the trait dialog; keep the tree in Preview

        path = self._preview_path_from_item(current)
        if path is None:
            self._set_status_file_link(None)
            return
        QTimer.singleShot(0, lambda p=path: self._open_file_in_preview_if_allowed(p))

    def _open_file_in_preview_if_allowed(self, path: Path) -> None:
        if self._file_list.suppress_selection_preview():
            return
        self._open_file_in_preview(path)

    def _resolve_te_masking_input_path(self, path: Path) -> Path | None:
        if is_te_hit_table_file(path):
            return path.expanduser().resolve()
        return None

    def _resolve_te_consensus_input_path(self, path: Path) -> Path | None:
        if self._is_fasta_file(path):
            return path.expanduser().resolve()
        return None

    def _resolve_te_hit_input_path(self, path: Path) -> Path | None:
        resolved = self._resolve_te_masking_input_path(path)
        if resolved is not None:
            return resolved
        if path.suffix.lower() in _FASTA_SUFFIXES:
            blast_txt = find_blast_table_for_query(path)
            if blast_txt is not None:
                return blast_txt.expanduser().resolve()
        return None

    def _resolve_gene_input_path(self, path: Path) -> Path | None:
        if is_gene_annotation_file(path):
            return path.expanduser().resolve()
        return None

    def _resolve_methylation_input_path(self, path: Path) -> Path | None:
        if is_bismark_cov_file(path):
            return path.expanduser().resolve()
        return None

    def _guess_gene_input_file(self) -> Path | None:
        for path in self._selected_files():
            resolved = self._resolve_gene_input_path(path)
            if resolved is not None:
                return resolved
        current = self._file_list.currentItem()
        if current is not None:
            path = self._preview_path_from_item(current)
            if path is not None:
                return self._resolve_gene_input_path(path)
        return None

    def _guess_methylation_input_file(self) -> Path | None:
        for path in self._selected_files():
            resolved = self._resolve_methylation_input_path(path)
            if resolved is not None:
                return resolved
        current = self._file_list.currentItem()
        if current is not None:
            path = self._preview_path_from_item(current)
            if path is not None:
                return self._resolve_methylation_input_path(path)
        return None

    def _guess_te_masking_input_file(self) -> Path | None:
        for path in self._selected_files():
            resolved = self._resolve_te_masking_input_path(path)
            if resolved is not None:
                return resolved
        current = self._file_list.currentItem()
        if current is not None:
            path = self._preview_path_from_item(current)
            if path is not None:
                return self._resolve_te_masking_input_path(path)
        return None

    def _guess_te_consensus_input_file(self) -> Path | None:
        for path in self._selected_files():
            resolved = self._resolve_te_consensus_input_path(path)
            if resolved is not None:
                return resolved
        current = self._file_list.currentItem()
        if current is not None:
            path = self._preview_path_from_item(current)
            if path is not None:
                return self._resolve_te_consensus_input_path(path)
        return None

    def _guess_te_hit_input_file(self) -> Path | None:
        for path in self._selected_files():
            resolved = self._resolve_te_hit_input_path(path)
            if resolved is not None:
                return resolved
        current = self._file_list.currentItem()
        if current is not None:
            path = self._preview_path_from_item(current)
            if path is not None:
                return self._resolve_te_hit_input_path(path)
        return None

    def _sync_dialog_input_from_files(self) -> None:
        self._sync_te_viz_dialog_input_from_files()
        self._sync_define_database_dialog_from_files()
        self._sync_plugin_run_dialog_from_files()

    def _panel_path_for_plugin_sync(self) -> Path | None:
        for path in self._selected_files():
            if path.exists():
                return path.resolve()
        folder_item = self._folder_tree.currentItem()
        if folder_item is not None:
            raw = folder_item.data(0, ROLE_PATH)
            if raw:
                folder = Path(str(raw)).expanduser()
                if folder.is_dir():
                    return folder.resolve()
        return None

    def _sync_plugin_run_dialog_from_files(self, path: Path | None = None) -> None:
        if path is None:
            return
        dialog = self._plugin_run_dialog
        if dialog is not None and dialog.isVisible():
            dialog.apply_path_from_panel(path)
            dialog.raise_()
            return
        for plugin in self._plugin_registry.all_plugins():
            if not isinstance(plugin, UserCommandPlugin):
                continue
            active = plugin._active_run_dialog
            if isinstance(active, UserCommandRunDialog) and active.isVisible():
                active.apply_path_from_panel(path)
                active.raise_()
                return

    def _guess_database_input_file(self) -> Path | None:
        for path in self._selected_files():
            if resolve_panel_database_path(path) is not None:
                return path.expanduser().resolve()
        current = self._file_list.currentItem()
        if current is not None:
            path = self._preview_path_from_item(current)
            if path is not None and resolve_panel_database_path(path) is not None:
                return path.expanduser().resolve()
        return None

    def _sync_genome_dialog_from_files(self, clicked_path: Path | None = None) -> None:
        dialog = self._genome_dialog
        if dialog is None or not dialog.isVisible() or clicked_path is None:
            return
        dialog.apply_file_from_panel(clicked_path)

    def _sync_define_database_dialog_from_files(self, clicked_path: Path | None = None) -> None:
        dialog = self._define_database_dialog
        if dialog is None or not dialog.isVisible() or clicked_path is None:
            return
        if dialog.apply_file_from_panel(clicked_path):
            self._schedule_raise_define_database_dialog()

    def _sync_define_database_dialog_from_folder(self, folder: Path) -> None:
        dialog = self._define_database_dialog
        if dialog is None or not dialog.isVisible():
            return
        if dialog.apply_folder_from_panel(folder):
            self._schedule_raise_define_database_dialog()

    def _schedule_raise_define_database_dialog(self) -> None:
        QTimer.singleShot(0, self._raise_define_database_dialog)

    def _raise_define_database_dialog(self) -> None:
        dialog = self._define_database_dialog
        if dialog is None:
            return
        # Avoid stealing focus / flashing if already visible; the textbox change
        # from panel click is enough.
        if dialog.isVisible():
            return
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _sync_te_viz_dialog_input_from_files(self, clicked_path: Path | None = None) -> None:
        dialog = self._te_viz_input_dialog
        if dialog is None or not dialog.isVisible() or clicked_path is None:
            return
        if hasattr(dialog, "apply_path_from_panel"):
            if dialog.apply_path_from_panel(clicked_path):
                self._schedule_raise_te_viz_dialog()
            return
        filled = False
        masking_path = self._resolve_te_masking_input_path(clicked_path)
        if masking_path is not None and hasattr(dialog, "set_te_input_file"):
            filled = dialog.set_te_input_file(masking_path) or filled
        consensus_path = self._resolve_te_consensus_input_path(clicked_path)
        if consensus_path is not None and hasattr(dialog, "set_te_consensus_input_file"):
            filled = dialog.set_te_consensus_input_file(consensus_path) or filled
        gene_path = self._resolve_gene_input_path(clicked_path)
        if gene_path is not None and hasattr(dialog, "set_gene_input_file"):
            filled = dialog.set_gene_input_file(gene_path) or filled
        methylation_path = self._resolve_methylation_input_path(clicked_path)
        if methylation_path is not None and hasattr(dialog, "set_methylation_input_file"):
            filled = dialog.set_methylation_input_file(methylation_path) or filled
        if filled:
            self._schedule_raise_te_viz_dialog()

    def _sync_te_viz_dialog_from_folder(self, folder: Path) -> None:
        """Fill a folder field (e.g. Repeat quantify output) from a Folders-panel click."""
        dialog = self._te_viz_input_dialog
        if dialog is None or not dialog.isVisible():
            return
        apply = getattr(dialog, "apply_folder_from_panel", None)
        if apply is not None and apply(folder):
            self._schedule_raise_te_viz_dialog()

    def _schedule_raise_te_viz_dialog(self) -> None:
        QTimer.singleShot(0, self._raise_te_viz_dialog)

    def _raise_te_viz_dialog(self) -> None:
        dialog = self._te_viz_input_dialog
        if dialog is None:
            return
        # Avoid stealing focus / flashing if already visible.
        if dialog.isVisible():
            return
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _present_te_viz_dialog(self, dialog: QDialog, on_accept) -> None:
        """Show a non-modal TE settings dialog so the Files panel stays clickable."""
        if self._te_viz_input_dialog is not None:
            self._te_viz_input_dialog.close()
        self._te_viz_input_dialog = dialog
        dialog.setModal(False)
        dialog.setWindowModality(Qt.WindowModality.NonModal)
        dialog.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, True)

        def _finished(code: int) -> None:
            if self._te_viz_input_dialog is dialog:
                self._te_viz_input_dialog = None
            if code == int(QDialog.DialogCode.Accepted):
                on_accept(dialog)
            dialog.deleteLater()

        dialog.finished.connect(_finished)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()
        # Do not auto-fill from the currently selected Files-panel item.
        # The user should click a file in the panel to fill inputs.

    def _is_fasta_file(self, path: Path) -> bool:
        return path.suffix.lower() in {".fa", ".fasta", ".fna", ".fas", ".faa"}

    def _validate_te_consensus_fastas(self, fasta_paths: list[Path]) -> list[Path] | None:
        if not fasta_paths:
            return []
        resolved: list[Path] = []
        for fasta_path in fasta_paths:
            fasta_path = fasta_path.expanduser()
            if not fasta_path.is_file():
                QMessageBox.warning(
                    self,
                    "Visualization",
                    f"TE consensus file not found:\n{fasta_path}",
                )
                return None
            if not self._is_fasta_file(fasta_path):
                QMessageBox.warning(
                    self,
                    "Visualization",
                    "Each TE consensus input must be a FASTA file (.fasta / .fa / .fna).",
                )
                return None
            resolved.append(fasta_path.resolve())
        return resolved

    def _validate_gene_input_files(self, gene_paths: list[Path]) -> list[Path] | None:
        if not gene_paths:
            return []
        resolved: list[Path] = []
        for gene_path in gene_paths:
            validated = self._validate_gene_input_file(gene_path)
            if validated is None:
                return None
            resolved.append(validated)
        return resolved

    def _validate_methylation_input_files(
        self, methylation_paths: list[Path]
    ) -> list[Path] | None:
        if not methylation_paths:
            return []
        resolved: list[Path] = []
        for meth_path in methylation_paths:
            if meth_path is None or not meth_path.is_file():
                QMessageBox.warning(
                    self,
                    "Visualization",
                    f"DNA methylation file not found:\n{meth_path}",
                )
                return None
            if not is_bismark_cov_file(meth_path):
                QMessageBox.warning(
                    self,
                    "Visualization",
                    "Each DNA methylation input must be a Bismark coverage file "
                    "(*.cov or *.cov.gz).",
                )
                return None
            resolved.append(meth_path.expanduser().resolve())
        return resolved

    def _resolve_te_viz_hit_paths_and_kinds(
        self,
        pending: _PendingTeVizBuild,
        *,
        created_paths: list[Path] | None = None,
    ) -> tuple[list[Path], list[str]] | None:
        te_paths = list(pending.masking_paths)
        kinds = [TE_INPUT_KIND_MASKING] * len(pending.masking_paths)
        for fasta in pending.consensus_fastas:
            blast_txt = self._find_blast_table_after_job(
                fasta,
                created_paths=created_paths,
            )
            if blast_txt is None:
                return None
            te_paths.append(blast_txt.expanduser().resolve())
            kinds.append(TE_INPUT_KIND_CONSENSUS)
        return te_paths, kinds

    def _find_blast_table_after_job(
        self,
        fasta: Path,
        *,
        created_paths: list[Path] | None = None,
    ) -> Path | None:
        """Locate ``*_b.txt`` (or legacy ``*_blast.txt``) for a query after Blast."""
        from repeatviewer.backends.blast_genome_plot import strip_blast_table_suffix

        found = find_blast_table_for_query(fasta)
        if found is not None:
            return found
        if not created_paths:
            return None
        fasta = fasta.expanduser().resolve()
        expected_names = {
            f"{fasta.name}_b.txt",
            f"{fasta.name}_blast.txt",
            f"{fasta.stem}_b.txt",
            f"{fasta.stem}_blast.txt",
        }
        for path in created_paths:
            resolved = path.expanduser().resolve()
            if not resolved.is_file():
                continue
            if resolved.name in expected_names:
                return resolved
            stripped = strip_blast_table_suffix(resolved.name)
            if stripped in {fasta.name, fasta.stem}:
                return resolved
        return None

    def _resolve_te_viz_hit_paths(self, pending: _PendingTeVizBuild) -> list[Path] | None:
        resolved = self._resolve_te_viz_hit_paths_and_kinds(pending)
        if resolved is None:
            return None
        return resolved[0]

    def _start_te_viz_from_dialog(
        self,
        *,
        masking_paths: list[Path],
        consensus_fastas: list[Path],
        gene_paths: list[Path],
        methylation_paths: list[Path],
        map_options: TeGenomeMapOptions,
    ) -> None:
        need_blast = [
            fasta
            for fasta in consensus_fastas
            if find_blast_table_for_query(fasta) is None
        ]
        if need_blast:
            self._pending_te_viz = _PendingTeVizBuild(
                masking_paths=list(masking_paths),
                consensus_fastas=list(consensus_fastas),
                gene_paths=list(gene_paths),
                methylation_paths=list(methylation_paths),
                map_options=map_options,
            )
            self._pending_te_genome_map = map_options
            self._pending_te_genome_map_fastas = list(need_blast)
            self._set_status_message("Running Blast on TE consensus, then TE genome map…", 0)
            # Density plot needs the *_b.txt table even if the user turned off
            # "Export BLAST text file" for ordinary Blast runs.
            self._apply_function(
                "blast",
                preselected=need_blast,
                blast_options=replace(self._blast_options, write_text=True),
            )
            return
        te_paths = list(masking_paths)
        te_kinds = [TE_INPUT_KIND_MASKING] * len(masking_paths)
        for fasta in consensus_fastas:
            blast_txt = find_blast_table_for_query(fasta)
            if blast_txt is not None:
                te_paths.append(blast_txt.expanduser().resolve())
                te_kinds.append(TE_INPUT_KIND_CONSENSUS)
        self._open_blast_genome_map(
            te_paths,
            gene_files=gene_paths,
            methylation_files=methylation_paths,
            te_input_kinds=te_kinds,
            map_options=map_options,
        )

    def _validate_te_hit_inputs(self, hit_paths: list[Path]) -> list[Path] | None:
        if not hit_paths:
            QMessageBox.warning(
                self,
                "Visualization",
                "Select at least one RepeatMasker .out or BLAST hit table.",
            )
            return None
        resolved: list[Path] = []
        for hit_path in hit_paths:
            if hit_path is None or not hit_path.is_file():
                QMessageBox.warning(
                    self,
                    "Visualization",
                    f"TE input file not found:\n{hit_path}",
                )
                return None
            if not is_te_hit_table_file(hit_path):
                QMessageBox.warning(
                    self,
                    "Visualization",
                    "Each TE input must be a RepeatMasker .out file or BLAST hit table (*_b.txt).",
                )
                return None
            resolved.append(hit_path.expanduser().resolve())
        return resolved

    def _validate_gene_input_file(self, gene_path: Path | None) -> Path | None:
        if gene_path is None:
            return None
        gene_path = gene_path.expanduser()
        if not gene_path.is_file():
            QMessageBox.warning(
                self,
                "Visualization",
                f"Annotation file not found:\n{gene_path}",
            )
            return None
        if not is_gene_annotation_file(gene_path):
            QMessageBox.warning(
                self,
                "Visualization",
                "Annotation input must be a GFF, GFF3, GTF, BED, or VCF file.",
            )
            return None
        return gene_path.resolve()

    def _validate_te_hit_input(self, hit_path: Path | None) -> Path | None:
        if hit_path is None or not hit_path.is_file():
            QMessageBox.warning(
                self,
                "Visualization",
                "Select a RepeatMasker .out or BLAST hit table.",
            )
            return None
        if not is_te_hit_table_file(hit_path):
            QMessageBox.warning(
                self,
                "Visualization",
                "Input must be a RepeatMasker .out file or BLAST hit table (*_b.txt).",
            )
            return None
        return hit_path.expanduser().resolve()

    def _on_file_item_clicked(self, item: FileListRow) -> None:
        raw = item.data(ROLE_PATH)
        panel_path: Path | None = None
        if raw:
            candidate = Path(str(raw)).expanduser()
            if candidate.exists():
                panel_path = candidate.resolve()
        if panel_path is not None:
            self._sync_genome_dialog_from_files(panel_path)
            self._sync_plugin_run_dialog_from_files(panel_path)
            self._sync_te_viz_dialog_input_from_files(panel_path)
            self._sync_define_database_dialog_from_files(panel_path)
            self._sync_add_annotation_dialog_from_files(panel_path)
            self._sync_add_blast_database_dialog_from_files(panel_path)
            self._sync_tree_trait_dialog_from_files(panel_path)

    def _open_tree_panels(self):
        from repeatviewer.widgets.preview_tabs import FilePreviewPage

        for index in range(self._preview_tabs.count()):
            page = self._preview_tabs.widget(index)
            if isinstance(page, FilePreviewPage) and page.tree_panel() is not None:
                yield page.tree_panel()

    def _tree_trait_dialog_open(self) -> bool:
        return any(panel.trait_dialog_open() for panel in self._open_tree_panels())

    def _sync_tree_trait_dialog_from_files(self, clicked_path: Path) -> None:
        """Tree view → Trait table dialog: a clicked file fills its path field."""
        for panel in self._open_tree_panels():
            if panel.apply_path_from_panel(clicked_path):
                return

    def _on_file_double_clicked(self, item: FileListRow) -> None:
        raw = item.data(ROLE_PATH)
        if not raw:
            return
        path = Path(str(raw))
        if path.suffix.lower() in _FASTA_SUFFIXES:
            self._open_in_sequence_viewer(path)
            return
        self._open_preview_tab(path, replace=True)

    def _open_selected_in_aliview(self) -> None:
        selected = self._selected_files()
        fasta = [p for p in selected if p.suffix.lower() in _FASTA_SUFFIXES]
        if not fasta:
            QMessageBox.information(
                self,
                "No FASTA selected",
                "Select one or more .fa / .fasta files to open in AliView.",
            )
            return
        for path in fasta:
            self._open_in_sequence_viewer(path)

    def _open_selected_in_text_editor(self, paths: list[Path] | None = None) -> None:
        selected = list(paths) if paths else self._selected_files()
        files = [p.expanduser().resolve() for p in selected if p.is_file()]
        if not files:
            QMessageBox.information(
                self,
                "No file selected",
                "Select one or more files to open in a text editor.",
            )
            return
        opened = 0
        failures: list[str] = []
        for path in files:
            try:
                open_in_text_editor(path)
                opened += 1
            except (OSError, RuntimeError) as exc:
                failures.append(f"{path.name}: {exc}")
        if opened:
            label = files[0].name if opened == 1 else f"{opened} files"
            self._set_status_message(f"Opened in text editor: {label}", 5000)
        if failures:
            QMessageBox.warning(
                self,
                "Text editor",
                "Some files could not be opened:\n\n" + "\n".join(failures[:6]),
            )

    def _open_in_sequence_viewer(self, path: Path) -> None:
        try:
            open_fasta(path, use_system_viewer=self._project.use_system_sequence_viewer)
            self._set_status_message(f"Opened: {path.name}", 5000)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Open sequence failed", str(exc))

    def _run_align_default(self) -> None:
        """Toolbar click: align each selected FASTA separately."""
        self._align_sequences(engine=self._align_options.combined_aligner)

    def _run_mafft_default(self) -> None:
        self._align_sequences(engine=ALIGNER_MAFFT)

    def _run_muscle_default(self) -> None:
        self._align_sequences(engine=ALIGNER_MUSCLE)

    def _count_sequences_in_fastas(self, seq_files: list[Path]) -> int:
        total = 0
        for path in seq_files:
            try:
                total += count_fasta_sequences(path)
            except Exception:  # noqa: BLE001
                continue
        return total

    def _align_selected_combined(self) -> None:
        """One-off combined MSA from all sequences in the selection."""
        selected = self._selected_files()
        seq_files = [p for p in selected if p.suffix.lower() in _FASTA_SUFFIXES]
        if not seq_files:
            QMessageBox.information(
                self,
                "Align combined",
                "Select one or more FASTA files with at least two sequences total.",
            )
            return
        if self._count_sequences_in_fastas(seq_files) < 2:
            QMessageBox.information(
                self,
                "Align combined",
                "Need at least two sequences in the selected FASTA file(s) "
                "for a combined alignment.",
            )
            return
        self._align_sequences(
            seq_files,
            mode=ALIGN_MODE_COMBINED,
            engine=self._align_options.combined_aligner,
        )

    def _mode_for_aligner(self, aligner: str) -> str:
        """'combined' or 'one_by_one' from the aligner's saved combine flag."""
        combined = (
            self._align_options.muscle_combined
            if aligner == ALIGNER_MUSCLE
            else self._align_options.mafft_combined
        )
        return ALIGN_MODE_COMBINED if combined else "one_by_one"

    def _open_align_dialog(self, *, initial_tab: str | None = None) -> None:
        """Align parameters: Cancel / Save / Run (toolbar and menu use the same dialog).

        Non-modal, so a click in the Files panel can fill the FastGA reference.
        FastGA compares the FASTA files selected when the window opened (the
        reference click changes the selection); MAFFT / MUSCLE use the
        selection at Run, as before.
        """
        fastga_files = self._fastga_candidate_files()
        dialog = AlignParamsDialog(
            initial=self._align_options,
            use_system_mafft=self._project.use_system_mafft,
            use_system_muscle=self._project.use_system_muscle,
            fastga_initial=self._fastga_options,
            fastga_files=fastga_files,
            initial_tab=initial_tab,
            parent=self,
        )
        dialog.finished.connect(
            lambda _code: QTimer.singleShot(0, self._repair_browser_chrome)
        )

        def _accepted(d) -> None:
            self._align_options = d.options()
            save_align_options(self._align_options)
            self._fastga_options = d.fastga_options()
            save_fastga_options(self._fastga_options)
            self._project.use_system_mafft = d.use_system_mafft()
            self._project.use_system_muscle = d.use_system_muscle()
            self._persist_workspace()
            if not d.should_run():
                self._set_status_message("Alignment settings saved", 4000)
                return
            aligner = d.selected_aligner()
            if aligner == ALIGNER_FASTGA:
                files = self._selected_files_for_fastga(
                    self._fastga_options, files=fastga_files or None
                )
                if files:
                    self._start_fastga_job(files, self._fastga_options)
                return
            self._align_sequences(
                engine=aligner, mode=self._mode_for_aligner(aligner)
            )

        self._present_te_viz_dialog(dialog, _accepted)

    def _quick_align_right_click(self) -> None:
        """Right-click Align: run the aligner chosen in the Options tab."""
        aligner = self._align_options.right_click_aligner
        if aligner == ALIGNER_FASTGA:
            self._run_fastga_default()
            return
        self._align_sequences(engine=aligner, mode=self._mode_for_aligner(aligner))

    def _modify_align_parameters(self) -> None:
        """Menu: same Cancel / Save / Run dialog as the Align toolbar button."""
        self._open_align_dialog()

    def _align_sequences(
        self,
        paths: list[Path] | None = None,
        *,
        mode: str | None = None,
        engine: str | None = None,
    ) -> None:
        opts = self._align_options
        engine = engine or opts.combined_aligner
        if mode is None:
            mode = "one_by_one"

        selected = paths or self._selected_files()
        seq_files = [p for p in selected if p.suffix.lower() in _FASTA_SUFFIXES]
        if not seq_files:
            QMessageBox.information(
                self,
                "No FASTA selected",
                "Select one or more .fa / .fasta files for alignment.",
            )
            return
        if mode == ALIGN_MODE_COMBINED and self._count_sequences_in_fastas(seq_files) < 2:
            QMessageBox.information(
                self,
                "Align combined",
                "Need at least two sequences in the selected FASTA file(s).",
            )
            return
        available = alignment_backends(
            use_system_mafft=self._project.use_system_mafft,
            use_system_muscle=self._project.use_system_muscle,
        )
        if engine not in available:
            tool_name = engine.upper()
            use_system = (
                self._project.use_system_muscle
                if engine == ALIGNER_MUSCLE
                else self._project.use_system_mafft
            )
            if use_system:
                detail = (
                    f"System {tool_name} is enabled but no '{engine}' command was found "
                    f"on PATH. Disable 'Use {tool_name} from the system PATH' in "
                    "Alignment parameters to use the bundled copy."
                )
            else:
                detail = (
                    f"The bundled {tool_name} does not support this operating-system "
                    f"architecture. Install {tool_name} on PATH and enable "
                    f"'Use {tool_name} from the system PATH' in Alignment parameters."
                )
            QMessageBox.warning(
                self,
                f"{tool_name} not found",
                detail,
            )
            return
        title = f"Align ({engine.upper()})"
        if mode == ALIGN_MODE_COMBINED:
            title = f"{title} (combined)"
        total_steps = 1 if mode == ALIGN_MODE_COMBINED else len(seq_files)
        show_log = bool(getattr(opts, "show_process_window", False))
        worker = AlignWorker(
            seq_files,
            options=opts,
            mode=mode,
            max_workers=self._project.batch_workers,
            use_system_mafft=self._project.use_system_mafft,
            use_system_muscle=self._project.use_system_muscle,
            engine=engine,
            stream_log=show_log,
            cancel=self._job_cancel,
        )
        self._attach_process_monitor(worker, title=title, enabled=show_log, tool_key="align")
        self._start_background_job(
            worker,
            job_title=title,
            total_steps=total_steps,
            source_files=seq_files,
            on_finished=self._on_align_finished,
            on_failed=self._on_align_failed,
        )

    @Slot(list)
    def _on_align_finished(self, created: list) -> None:
        created_paths = [Path(p) for p in created if Path(p).is_file()]
        source_paths = getattr(self, "_last_job_source_files", [])
        if self._current_folder:
            self._populate_file_list(
                self._current_folder,
                new_outputs=created_paths,
                source_files=source_paths,
                open_preview=False,
            )
        self._auto_open_outputs(created_paths)
        names = ", ".join(p.name for p in created_paths[:6])
        if len(created_paths) > 6:
            names += ", …"
        self._set_status_message(
            f"Alignment completed — {len(created_paths)} file(s): {names}",
            12000,
            override_file_link=True,
        )

    @Slot(str)
    def _on_align_failed(self, error: str) -> None:
        if error == CANCELLED_MESSAGE:
            self._set_status_message("Stopped", 6000)
            return
        body = error.strip()
        QMessageBox.critical(self, "Alignment failed", body)

    def _run_blast(self) -> None:
        self._run_blast_program(BLAST_PROGRAM_BLASTN)

    def _run_blast_program(self, program: str) -> None:
        options = replace(self._blast_options, program=program)
        self._apply_function("blast", blast_options=options)

    def _quick_blast_right_click(self) -> None:
        """Right-click on BLAST: run the program chosen in the Right Click tab."""
        program = getattr(
            self._blast_options, "right_click_program", BLAST_PROGRAM_BLASTN
        )
        self._run_blast_program(program)

    def _open_blast_dialog(self) -> None:
        """Click BLAST -> parameters dialog. Run launches the active tab's program.

        Non-modal so Add/Manage BLAST can receive Files-panel clicks for FASTA
        paths. Cancel discards; Save/Run persist settings (and Run starts the job).
        """
        from repeatviewer.database_settings import reload_blast_databases_into_project

        if self._blast_params_dialog is not None and self._blast_params_dialog.isVisible():
            self._blast_params_dialog.raise_()
            self._blast_params_dialog.activateWindow()
            return

        reload_blast_databases_into_project(self._project)
        dialog = BlastParamsDialog(
            project=self._project, initial=self._blast_options, parent=self
        )
        dialog.setModal(False)
        dialog.setWindowModality(Qt.WindowModality.NonModal)
        dialog.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, True)
        self._blast_params_dialog = dialog
        dialog.finished.connect(self._on_blast_params_dialog_finished)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _on_blast_params_dialog_finished(self, result: int) -> None:
        dialog = self._blast_params_dialog
        self._blast_params_dialog = None
        QTimer.singleShot(0, self._repair_browser_chrome)
        if dialog is None:
            return
        if result != int(QDialog.DialogCode.Accepted):
            return
        self._blast_options = dialog.options()
        save_blast_options(self._blast_options)
        self._project.use_system_blast = dialog.use_system_blast()
        self._persist_workspace()
        if dialog.should_run():
            self._apply_function("blast", blast_options=self._blast_options)
        else:
            self._set_status_message("BLAST settings saved", 4000)

    def _change_blast_parameters(self) -> None:
        """Menu: same Cancel / Save / Run dialog as the BLAST toolbar button."""
        self._open_blast_dialog()

    def _sync_add_blast_database_dialog_from_files(
        self, clicked_path: Path | None = None
    ) -> None:
        """Fill Add BLAST database FASTA from a Files-panel click."""
        if clicked_path is None:
            return
        from repeatviewer.views.blast_database_dialog import BlastDatabaseManagerDialog

        direct = getattr(self, "_add_blast_database_dialog", None)
        if direct is not None and direct.isVisible():
            if direct.apply_path_from_panel(clicked_path):
                return

        managers: list[BlastDatabaseManagerDialog] = []
        params = self._blast_params_dialog
        if params is not None:
            mgr = getattr(params, "_blast_manager_dialog", None)
            if isinstance(mgr, BlastDatabaseManagerDialog):
                managers.append(mgr)
        for mgr in self.findChildren(BlastDatabaseManagerDialog):
            if mgr not in managers:
                managers.append(mgr)

        for manager in managers:
            if manager.isVisible() and manager.apply_path_from_panel(clicked_path):
                # Manager forwards to Add and keeps Add in front — do not raise Manage.
                return

    def _add_fasta_as_blast_database(self, path: Path) -> None:
        from repeatviewer.database_settings import save_blast_databases_to_settings
        from repeatviewer.views.blast_database_dialog import AddBlastDatabaseDialog

        existing = getattr(self, "_add_blast_database_dialog", None)
        if existing is not None and existing.isVisible():
            existing.apply_path_from_panel(path)
            existing.raise_()
            existing.activateWindow()
            return

        existing_names = {db.name for db in self._project.blast_databases}
        dialog = AddBlastDatabaseDialog(
            existing_names=existing_names,
            use_system_blast=self._project.use_system_blast,
            initial_fasta=path,
            parent=self,
        )
        self._add_blast_database_dialog = dialog
        dialog.finished.connect(self._on_add_blast_database_dialog_finished)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _on_add_blast_database_dialog_finished(self, result: int) -> None:
        from repeatviewer.database_settings import save_blast_databases_to_settings

        dialog = getattr(self, "_add_blast_database_dialog", None)
        self._add_blast_database_dialog = None
        if dialog is None:
            return
        if result != int(QDialog.DialogCode.Accepted):
            return
        db = dialog.result_database()
        if db is None:
            return
        self._project.blast_databases.append(db)
        save_blast_databases_to_settings(self._project.blast_databases)
        self._persist_workspace()
        self._set_status_message(f"Added BLAST database “{db.name}”", 4000)
        params = self._blast_params_dialog
        if params is not None and params.isVisible():
            refresh = getattr(params, "_refresh_database_combos", None)
            if callable(refresh):
                refresh()

    def _add_annotation_for_fasta(self, path: Path) -> None:
        """Open a non-modal Add annotation dialog for ``path``.

        Deferred one tick so the Files context menu can finish and release the
        click-to-preview suppress flag before the dialog opens.
        """
        self._file_list.clear_suppress_selection_preview()
        QTimer.singleShot(0, lambda p=path: self._run_add_annotation_for_fasta(p))

    def _sync_add_annotation_dialog_from_files(self, clicked_path: Path | None = None) -> None:
        dialog = self._add_annotation_dialog
        if dialog is None or not dialog.isVisible() or clicked_path is None:
            return
        if dialog.apply_path_from_panel(clicked_path):
            dialog.raise_()

    def _run_add_annotation_for_fasta(self, path: Path) -> None:
        from repeatviewer.views.add_annotation_dialog import AddAnnotationDialog

        self._file_list.clear_suppress_selection_preview()
        fasta = path.expanduser().resolve()
        if not fasta.is_file():
            QMessageBox.warning(self, "Add annotation", f"FASTA not found:\n{fasta}")
            return
        if self._add_annotation_dialog is not None:
            existing = self._add_annotation_dialog
            if existing.isVisible() and existing.fasta_path() == fasta:
                existing.raise_()
                existing.activateWindow()
                return
            existing.close()
            self._add_annotation_dialog = None

        dialog = AddAnnotationDialog(fasta_path=fasta, parent=self)
        self._add_annotation_dialog = dialog
        dialog.finished.connect(self._on_add_annotation_dialog_finished)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _on_add_annotation_dialog_finished(self, result: int) -> None:
        from repeatviewer.backends.sequence_features import (
            features_gff_path,
            import_annotation_as_rv_gff3,
        )

        dialog = self._add_annotation_dialog
        self._add_annotation_dialog = None
        self._file_list.clear_suppress_selection_preview()
        if dialog is None:
            return
        try:
            if result != int(QDialog.DialogCode.Accepted):
                return
            fasta = dialog.fasta_path()
            annotation = dialog.result_path()
            if annotation is None:
                return
            dest = features_gff_path(fasta)
            if dest.is_file() and dest.resolve() == annotation:
                QMessageBox.information(
                    self,
                    "Add annotation",
                    f"{dest.name} is already the RepeatViewer annotation for this FASTA.",
                )
                return
            try:
                result = import_annotation_as_rv_gff3(fasta, annotation)
            except Exception as exc:  # noqa: BLE001
                QMessageBox.warning(self, "Add annotation", str(exc))
                return

            gff_path = result.gff_path
            count = result.feature_count
            seqs = result.sequence_count
            self._sync_sidecar_in_file_list(gff_path, select=False)
            if self._preview_tabs is not None:
                self._preview_tabs.reload_paths([fasta])
            if seqs > 1:
                status = (
                    f"Added {count:,} feature(s) across {seqs:,} sequences → {gff_path.name}"
                )
            else:
                status = f"Added {count:,} feature(s) → {gff_path.name}"
            self._set_status_message(
                status,
                8000,
                override_file_link=True,
            )
            self._set_status_file_link(gff_path)
        finally:
            dialog.deleteLater()

    def _open_blast_genome_map_from_selection(
        self,
        *,
        map_options: TeGenomeMapOptions | None = None,
    ) -> None:
        selected = self._selected_files()
        if not selected:
            QMessageBox.information(
                self,
                "TE genome map",
                "Select one FASTA (with *_b.txt) or a hit table (*_b.txt / RepeatMasker .out).",
            )
            return
        blast_targets: list[Path] = []
        missing_blast_fastas: list[Path] = []
        invalid: list[str] = []
        for path in selected:
            if path.suffix.lower() in _FASTA_SUFFIXES:
                blast_txt = find_blast_table_for_query(path)
                if blast_txt is not None:
                    blast_targets.append(blast_txt)
                else:
                    missing_blast_fastas.append(path)
            elif is_te_hit_table_file(path):
                blast_targets.append(path)
            else:
                invalid.append(path.name)

        if invalid:
            QMessageBox.warning(
                self,
                "TE genome map",
                "Skipped unsupported selections:\n" + "\n".join(invalid[:8]),
            )

        if blast_targets:
            self._enqueue_genome_map_jobs(blast_targets, map_options)
        if missing_blast_fastas:
            self._start_blast_then_te_genome_map(missing_blast_fastas, map_options)
        if not blast_targets and not missing_blast_fastas:
            QMessageBox.warning(
                self,
                "TE genome map",
                "Select FASTA files, *_b.txt tables, or RepeatMasker .out files.",
            )

    def _enqueue_genome_map_jobs(
        self,
        blast_tables: list[Path],
        map_options: TeGenomeMapOptions | None,
    ) -> None:
        for bt in blast_tables:
            self._queued_genome_map_jobs.append((bt.expanduser().resolve(), map_options))
        self._start_next_genome_map_job()

    def _start_next_genome_map_job(self) -> None:
        if self._genome_map_thread is not None:
            return
        if not self._queued_genome_map_jobs:
            return
        blast_txt, map_options = self._queued_genome_map_jobs.pop(0)
        self._open_blast_genome_map(blast_txt, map_options=map_options)

    def _open_blast_divergence_from_selection(self) -> None:
        selected = self._selected_files()
        if not selected:
            QMessageBox.information(
                self,
                "TE divergence plot",
                "Select one FASTA (with *_b.txt) or a hit table (*_b.txt / RepeatMasker .out).",
            )
            return
        if len(selected) > 1:
            QMessageBox.information(
                self,
                "TE divergence plot",
                "Select a single FASTA, BLAST table, or RepeatMasker .out file.",
            )
            return
        path = selected[0]
        if path.suffix.lower() in _FASTA_SUFFIXES:
            blast_txt = find_blast_table_for_query(path)
            if blast_txt is not None:
                self._open_blast_divergence_plot(blast_txt)
            else:
                QMessageBox.warning(
                    self,
                    "TE divergence plot",
                    f"No BLAST table found for:\n{path.name}\n\n"
                    f"Expected {blast_txt_path_for_query(path).name} "
                    f"(or {path.name}_b.txt). Run Blast first.",
                )
            return
        if is_te_hit_table_file(path):
            self._open_blast_divergence_plot(path)
            return
        QMessageBox.warning(
            self,
            "TE divergence plot",
            "Select a FASTA with BLAST results, a *_b.txt table, or a RepeatMasker .out file.",
        )

    def _open_blast_divergence_plot(
        self,
        te_input_files: list[Path] | Path,
        *,
        bin_width: float | None = None,
    ) -> None:
        te_paths = (
            [te_input_files.expanduser().resolve()]
            if isinstance(te_input_files, Path)
            else [p.expanduser().resolve() for p in te_input_files]
        )
        te_paths = [p for p in te_paths if p.is_file()]
        if not te_paths:
            QMessageBox.warning(
                self,
                "TE divergence plot",
                "No TE input files found.\n\nRun Blast or select RepeatMasker .out / BLAST tables.",
            )
            return
        primary = te_paths[0]
        if len(te_paths) == 1:
            label = te_hit_table_stem(primary)
        else:
            label = f"{te_hit_table_stem(primary)} + {len(te_paths) - 1} more"
        proposed = divergence_plot_html_path_for_te_inputs(te_paths)
        output_html = self._confirm_plot_output_path(
            proposed,
            title="TE divergence plot",
        )
        if output_html is None:
            return
        bw = bin_width if bin_width is not None else self._te_divergence_plot_options.bin_width
        self._set_status_message("Building TE divergence plot…", 0)
        try:
            html_path = write_blast_divergence_plot_html(
                te_paths,
                query_label=label,
                bin_width=bw,
                output_html=output_html,
            )
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "TE divergence plot failed", str(exc))
            return
        if html_path is None:
            QMessageBox.warning(self, "TE divergence plot", "Could not build the divergence plot.")
            return
        self._open_divergence_in_preview(html_path)
        outputs = [html_path.resolve()]
        source = self._fasta_for_blast_table(primary) or primary
        if self._current_folder and self._current_folder.resolve() == html_path.parent.resolve():
            self._populate_file_list(html_path.parent.resolve(), new_outputs=outputs, source_files=[source])
        else:
            self._navigate_to_folder(html_path.parent.resolve(), new_outputs=outputs, source_files=[source])
        self._set_status_message(f"TE divergence plot: {html_path.name}", 8000, override_file_link=True)

    def _open_divergence_in_preview(self, html_path: Path) -> None:
        """Open divergence HTML in Preview with delayed reload for reliability."""
        from repeatviewer.widgets.preview_tabs import FilePreviewPage

        resolved = html_path.expanduser().resolve()
        if not resolved.is_file():
            return

        def _load() -> None:
            self._ensure_preview_visible()
            if self._preview_tabs is None:
                return
            self._preview_tabs.open_file(resolved, replace=True)
            page = self._preview_tabs.currentWidget()
            if isinstance(page, FilePreviewPage):
                page.reload_from_disk()

        _load()
        QTimer.singleShot(80, _load)
        QTimer.singleShot(450, _load)

    def _any_background_job_running(self) -> bool:
        return (
            self._worker_thread is not None
            or self._genome_map_thread is not None
            or self._benchmark_thread is not None
        )

    def _open_confusion_matrix_benchmark(self) -> None:
        selected = self._selected_files()
        query_files = [p for p in selected if is_repeatmasker_out_file(p)]
        if not query_files:
            QMessageBox.information(
                self,
                "Confusion matrix benchmarker",
                "Select one or more RepeatMasker .out files (query) in the Files panel.",
            )
            return
        if not self._project.genome_file or not self._project.genome_file.is_file():
            QMessageBox.warning(
                self,
                "Confusion matrix benchmarker",
                f"Set the genome FASTA under {self._genome_menu_hint()} first.",
            )
            return
        try:
            require_bedtools()
        except RuntimeError as exc:
            QMessageBox.warning(self, "Confusion matrix benchmarker", str(exc))
            return

        dialog = ConfusionMatrixBenchmarkDialog(
            initial=self._te_confusion_benchmark_options,
            parent=self,
        )
        if not exec_modal_dialog(dialog):
            return
        opts = dialog.options()
        self._te_confusion_benchmark_options = opts
        save_te_confusion_matrix_benchmark_options(opts)

        home_text = opts.tetrimmer_home.strip()
        if not home_text:
            detected = default_tetrimmer_home()
            if detected is None:
                QMessageBox.warning(
                    self,
                    "Confusion matrix benchmarker",
                    "Set the TEtrimmer repository folder in the dialog, or set "
                    "environment variable TETRIMMER_HOME.",
                )
                return
            home = detected
        else:
            home = Path(home_text).expanduser().resolve()

        reference = Path(opts.reference_out).expanduser().resolve()
        if not reference.is_file():
            QMessageBox.warning(
                self,
                "Confusion matrix benchmarker",
                "Choose a reference RepeatMasker .out file (curated TE library on the genome).",
            )
            return

        worker = ConfusionMatrixBenchmarkWorker(
            query_files,
            reference_out=reference,
            genome_fasta=self._project.genome_file,
            tetrimmer_home=home,
            less_classification=opts.less_classification,
            debug=opts.debug,
            cancel=self._job_cancel,
        )
        self._start_te_benchmark_job(
            worker,
            job_title="te_confusion_benchmark",
            total_steps=len(query_files),
            source_files=query_files,
        )

    def _open_te_library_benchmark(self) -> None:
        selected = self._selected_files()
        fasta_files = [p for p in selected if p.suffix.lower() in _FASTA_SUFFIXES]
        query_hint = fasta_files[0] if len(fasta_files) == 1 else None

        dialog = TeLibraryBenchmarkDialog(
            initial=self._te_library_benchmark_options,
            query_file=query_hint,
            parent=self,
        )
        if not exec_modal_dialog(dialog):
            return
        opts = dialog.options()
        self._te_library_benchmark_options = opts
        save_te_library_benchmark_options(opts)

        query = Path(opts.query_fasta).expanduser().resolve()
        standard = Path(opts.standard_fasta).expanduser().resolve()
        if not query.is_file():
            QMessageBox.warning(
                self,
                "TE library benchmarker",
                "Choose a query TE consensus library FASTA.",
            )
            return
        if not standard.is_file():
            QMessageBox.warning(
                self,
                "TE library benchmarker",
                "Choose a standard (curated) TE library FASTA.",
            )
            return

        rm_out: Path | None = None
        if opts.rm_out.strip():
            rm_out = Path(opts.rm_out).expanduser().resolve()
            if not rm_out.is_file():
                QMessageBox.warning(
                    self,
                    "TE library benchmarker",
                    f"RepeatMasker .out not found:\n{rm_out}",
                )
                return
        elif not opts.run_repeatmasker:
            QMessageBox.warning(
                self,
                "TE library benchmarker",
                "Provide an existing RepeatMasker .out file or enable Run RepeatMasker.",
            )
            return

        worker = TeLibraryBenchmarkWorker(
            query_fasta=query,
            standard_fasta=standard,
            rm_out=rm_out,
            run_repeatmasker=opts.run_repeatmasker,
            threads=opts.threads,
            cancel=self._job_cancel,
        )
        self._start_te_benchmark_job(
            worker,
            job_title="te_library_benchmark",
            total_steps=1,
            source_files=[query],
        )

    def _start_te_benchmark_job(
        self,
        worker: ConfusionMatrixBenchmarkWorker | TeLibraryBenchmarkWorker,
        *,
        job_title: str,
        total_steps: int,
        source_files: list[Path],
    ) -> None:
        if self._any_background_job_running():
            QMessageBox.information(
                self,
                "Job already running",
                "Wait for the current job to finish before starting another.",
            )
            return

        self._last_benchmark_source_files = list(source_files)
        self._clear_job_file_highlights()
        self._job_cancel.reset()
        self._show_job_progress(job_title, total_steps)
        self._set_job_actions_enabled(False)
        if self._stop_btn is not None:
            self._stop_btn.setEnabled(True)

        self._benchmark_thread = QThread(self)
        self._benchmark_worker = worker
        self._benchmark_worker.moveToThread(self._benchmark_thread)
        self._benchmark_thread.started.connect(self._benchmark_worker.run)
        queued = Qt.ConnectionType.QueuedConnection
        self._benchmark_worker.progress.connect(self._on_job_progress, queued)
        self._benchmark_worker.finished.connect(self._on_te_benchmark_finished, queued)
        self._benchmark_worker.failed.connect(self._on_te_benchmark_failed, queued)
        self._benchmark_worker.finished.connect(self._benchmark_thread.quit)
        self._benchmark_worker.failed.connect(self._benchmark_thread.quit)
        self._benchmark_thread.finished.connect(self._cleanup_benchmark_worker)
        self._benchmark_thread.finished.connect(self._update_status)
        self._benchmark_thread.start()

    @Slot(list)
    def _on_te_benchmark_finished(self, created_paths: list[Path]) -> None:
        self._log_job_finished(
            getattr(self, "_active_job_log_title", "benchmark"), created_paths
        )
        self._finish_background_job()
        try:
            paths = [
                Path(p).expanduser().resolve()
                for p in created_paths
                if Path(p).is_file()
            ]
            sources = self._last_benchmark_source_files
            self._refresh_folders_after_job(paths, sources)
            self._auto_open_outputs(paths)
            title = self._job_display_title("te_confusion_benchmark")
            if paths and any(p.name.endswith("_library_benchmark.html") for p in paths):
                title = self._job_display_title("te_library_benchmark")
            summary = f"{title} finished"
            if paths:
                names = ", ".join(p.name for p in paths[:4])
                if len(paths) > 4:
                    names += ", …"
                summary += f" — created: {names}"
            self._set_status_message(summary, 12000, override_file_link=True)
        finally:
            self._schedule_after_job_ui_ready()

    @Slot(str)
    def _on_te_benchmark_failed(self, error: str) -> None:
        self._log_job_failed(getattr(self, "_active_job_log_title", "benchmark"), error)
        self._hide_job_progress()
        self._set_job_actions_enabled(True)
        if error == CANCELLED_MESSAGE:
            self._set_status_message("Stopped benchmark", 6000, override_file_link=True)
            return
        self._set_status_message("Benchmark failed", 8000, override_file_link=True)
        QMessageBox.critical(self, "Benchmark failed", error)

    def _cleanup_benchmark_worker(self) -> None:
        if self._benchmark_worker:
            self._benchmark_worker.deleteLater()
            self._benchmark_worker = None
        if self._benchmark_thread:
            self._benchmark_thread.deleteLater()
            self._benchmark_thread = None

    def _open_te_density_plot_dialog(self, *, input_file: Path | None = None) -> None:
        if self._project.is_pan_genome_mode:
            QMessageBox.information(
                self,
                "Visualization",
                "Visualization is available only in single-genome mode.\n\n"
                f"Switch to single genome under {self._genome_menu_hint()}.",
            )
            return
        if not self._project.genome_file or not self._project.genome_file.is_file():
            QMessageBox.warning(
                self,
                "Visualization",
                f"Set the genome FASTA under {self._genome_menu_hint()} first.",
            )
            return
        chrom_names = self._genome_chromosome_names()
        chrom_sizes = self._genome_chromosome_sizes()
        if not chrom_names:
            QMessageBox.warning(
                self,
                "Visualization",
                "Could not read chromosome names from the genome FASTA.",
            )
            return
        masking_input: Path | None = None
        consensus_input: Path | None = None
        if input_file is not None:
            if is_te_hit_table_file(input_file):
                masking_input = input_file
            elif self._is_fasta_file(input_file):
                consensus_input = input_file
        else:
            masking_input = self._guess_te_masking_input_file()
            consensus_input = self._guess_te_consensus_input_file()

        def _on_accept(dialog: TeDensityPlotDialog) -> None:
            self._te_genome_map_options = dialog.options()
            save_te_genome_map_options(self._te_genome_map_options)
            masking_files: list[Path] = []
            if dialog.te_masking_files():
                validated_masking = self._validate_te_hit_inputs(dialog.te_masking_files())
                if validated_masking is None:
                    return
                masking_files = validated_masking
            consensus_files = self._validate_te_consensus_fastas(dialog.te_consensus_files())
            if consensus_files is None:
                return
            gene_files = self._validate_gene_input_files(dialog.gene_input_files())
            if dialog.gene_input_files() and gene_files is None:
                return
            methylation_files = self._validate_methylation_input_files(
                dialog.methylation_input_files()
            )
            if dialog.methylation_input_files() and methylation_files is None:
                return
            if not masking_files and not consensus_files and not gene_files and not methylation_files:
                QMessageBox.warning(
                    self,
                    "Visualization",
                    "Provide at least one TE masking file, TE consensus FASTA, gene annotation, "
                    "or DNA methylation file.",
                )
                return
            self._start_te_viz_from_dialog(
                masking_paths=masking_files,
                consensus_fastas=consensus_files,
                gene_paths=gene_files or [],
                methylation_paths=methylation_files or [],
                map_options=self._te_genome_map_options,
            )

        dialog = TeDensityPlotDialog(
            initial=self._te_genome_map_options,
            chromosome_names=chrom_names,
            chromosome_sizes=chrom_sizes,
            input_file=masking_input,
            consensus_input_file=consensus_input,
            parent=self,
        )
        self._present_te_viz_dialog(dialog, _on_accept)

    def _open_te_divergence_plot_dialog(self, *, input_file: Path | None = None) -> None:
        if input_file is None:
            input_file = self._guess_te_hit_input_file()

        def _on_accept(dialog: TeDivergencePlotDialog) -> None:
            self._te_divergence_plot_options = dialog.options()
            save_te_divergence_plot_options(self._te_divergence_plot_options)
            te_files = self._validate_te_hit_inputs(dialog.te_input_files())
            if te_files is None:
                return
            self._open_blast_divergence_plot(
                te_files,
                bin_width=self._te_divergence_plot_options.bin_width,
            )

        dialog = TeDivergencePlotDialog(
            initial=self._te_divergence_plot_options,
            input_file=input_file,
            parent=self,
        )
        self._present_te_viz_dialog(dialog, _on_accept)

    def _open_te_abundance_plot_dialog(self, *, input_file: Path | None = None) -> None:
        selected = [p for p in self._selected_files() if is_repeatmasker_out_file(p)]
        initial_files = [input_file] if input_file is not None else selected

        def _on_accept(dialog: TeAbundancePlotDialog) -> None:
            sources = dialog.input_files()
            if not sources:
                QMessageBox.warning(self, "TE abundance plot", "Select a RepeatMasker .out file.")
                return
            for source in sources:
                if not source.is_file() or not is_repeatmasker_out_file(source):
                    QMessageBox.warning(
                        self,
                        "TE abundance plot",
                        f"Select only valid RepeatMasker .out files.\n\nInvalid: {source}",
                    )
                    return
            self._start_te_abundance_plot(
                sources,
                metric=dialog.y_axis_metric(),
                top_n=dialog.top_n(),
            )

        dialog = TeAbundancePlotDialog(input_files=initial_files or None, parent=self)
        self._present_te_viz_dialog(dialog, _on_accept)

    def _open_msa_plot_dialog(self, *, input_file: Path | None = None) -> None:
        if input_file is None:
            selected = [
                p for p in self._selected_files()
                if p.is_file() and self._is_fasta_file(p)
            ]
            input_file = selected[0] if selected else None

        # Reload so palette migrations (Default colours) apply without restart.
        self._msa_plot_options = load_msa_plot_options()

        def _on_accept(dialog: MsaPlotDialog) -> None:
            self._msa_plot_options = dialog.options()
            save_msa_plot_options(self._msa_plot_options)
            if not dialog.should_run():
                self._set_status_message("MSA plot settings saved", 4000)
                return
            files = dialog.input_files()
            if not files:
                QMessageBox.warning(self, "MSA plot", "Select an aligned FASTA file.")
                return
            source = files[0].expanduser().resolve()
            if not source.is_file() or not self._is_fasta_file(source):
                QMessageBox.warning(
                    self,
                    "MSA plot",
                    f"Select a valid FASTA alignment.\n\nInvalid: {source}",
                )
                return
            self._start_msa_plot(source, options=self._msa_plot_options)

        dialog = MsaPlotDialog(
            initial=self._msa_plot_options,
            input_file=input_file,
            parent=self,
        )
        self._present_te_viz_dialog(dialog, _on_accept)

    def _start_msa_plot(self, source: Path, *, options: MsaPlotJobOptions) -> None:
        proposed = msa_plot_output_path(source, output_format=options.output_format)
        output_path = self._confirm_plot_output_path(proposed, title="MSA plot")
        if output_path is None:
            return
        self._start_background_job(
            FunctionWorker(
                self._run_msa_plot_job,
                "msa_plot",
                [source],
                {"options": options, "output_path": output_path},
            ),
            job_title="msa_plot",
            total_steps=1,
            source_files=[source],
            on_finished=self._on_msa_plot_finished,
            on_failed=self._on_worker_failed,
        )

    def _run_msa_plot_job(
        self,
        func_name: str,
        sources: list[Path],
        progress,
        *,
        options: MsaPlotJobOptions,
        output_path: Path | None = None,
    ) -> dict:
        progress(0, "Building MSA plot…")
        if self._job_cancel.requested:
            raise JobCancelledError()
        plot_opts = MsaPlotOptions(
            sequence_type=options.sequence_type,
            output_format=options.output_format,
            dpi=options.dpi,
            sequence_height_cm=options.sequence_height_cm,
            sequence_distance_cm=options.sequence_distance_cm,
            nucleotide_colors=dict(options.nucleotide_colors),
            amino_acid_colors=dict(options.amino_acid_colors),
        )
        written = write_msa_plot(sources[0], output_path, options=plot_opts)
        if self._job_cancel.requested:
            raise JobCancelledError()
        progress(1, f"Finished: {written.name}")
        return {
            "func_name": func_name,
            "created_outputs": [str(written.resolve())],
            "source_files": [str(p.resolve()) for p in sources],
            "total": 1,
            "failures": [],
        }

    @Slot(dict)
    def _on_msa_plot_finished(self, result: dict) -> None:
        created = [
            Path(str(raw)).expanduser().resolve()
            for raw in result.get("created_outputs", [])
        ]
        created = [p for p in created if p.is_file()]
        sources = [
            Path(str(raw)).expanduser().resolve()
            for raw in result.get("source_files", [])
        ]
        if not created:
            self._set_status_message("MSA plot finished", 8000, override_file_link=True)
            return
        out = created[0]
        self._refresh_folders_after_job(created, sources)
        if out.suffix.lower() == ".png":
            self._open_divergence_in_preview(out)
        self._set_status_message(f"MSA plot: {out.name}", 8000, override_file_link=True)

    def _start_te_abundance_plot(
        self,
        sources: list[Path],
        *,
        metric: str,
        top_n: int,
    ) -> None:
        proposed = te_abundance_plot_html_path_for_sources(sources)
        output_html = self._confirm_plot_output_path(
            proposed,
            title="TE abundance plot",
        )
        if output_html is None:
            return
        self._start_background_job(
            FunctionWorker(
                self._run_te_abundance_plot_job,
                "te_abundance_plot",
                sources,
                {"metric": metric, "top_n": top_n, "output_html": output_html},
            ),
            job_title="te_abundance_plot",
            total_steps=max(1, len(sources)),
            source_files=sources,
            on_finished=self._on_te_abundance_plot_finished,
            on_failed=self._on_worker_failed,
        )

    def _run_te_abundance_plot_job(
        self,
        func_name: str,
        sources: list[Path],
        progress,
        *,
        metric: str,
        top_n: int,
        output_html: Path | None = None,
    ) -> dict:
        progress(0, "Building TE abundance plot…")
        if self._job_cancel.requested:
            raise JobCancelledError()
        html_path = write_te_abundance_plot_html(
            sources,
            metric=metric,
            top_n=top_n,
            output_html=output_html,
        )
        if self._job_cancel.requested:
            raise JobCancelledError()
        progress(max(1, len(sources)), f"Finished: {html_path.name}")
        return {
            "func_name": func_name,
            "created_outputs": [str(html_path.resolve())],
            "source_files": [str(p.resolve()) for p in sources],
            "total": len(sources),
            "failures": [],
        }

    @Slot(dict)
    def _on_te_abundance_plot_finished(self, result: dict) -> None:
        created = [
            Path(str(raw)).expanduser().resolve()
            for raw in result.get("created_outputs", [])
        ]
        created = [p for p in created if p.is_file()]
        sources = [
            Path(str(raw)).expanduser().resolve()
            for raw in result.get("source_files", [])
        ]
        if not created:
            self._set_status_message("TE abundance plot finished", 8000, override_file_link=True)
            return
        html_path = created[0]
        self._open_divergence_in_preview(html_path)
        self._refresh_folders_after_job(created, sources)
        self._set_status_message(f"TE abundance plot: {html_path.name}", 8000, override_file_link=True)

    def _warm_genome_chromosome_cache(self) -> None:
        """Pre-read genome FASTA metadata so TE genome map dialogs open quickly."""
        if self._project.genome_file and self._project.genome_file.is_file():
            self._genome_chromosome_names()

    def _genome_chromosome_names(self) -> list[str]:
        """Chromosome list for dialogs (cached; avoids re-reading the whole genome FASTA)."""
        names, _ = self._genome_chromosome_meta()
        return names

    def _genome_chromosome_sizes(self) -> dict[str, int]:
        """Chromosome lengths (bp) for dialogs (cached with names)."""
        _, sizes = self._genome_chromosome_meta()
        return sizes

    def _genome_chromosome_meta(self) -> tuple[list[str], dict[str, int]]:
        if not self._project.genome_file or not self._project.genome_file.is_file():
            return [], {}
        path = self._project.genome_file.expanduser().resolve()
        if self._genome_chrom_cache_path == path:
            return self._genome_chrom_names_cache, self._genome_chrom_sizes_cache
        self._genome_chrom_cache_path = path
        self._genome_chrom_names_cache = chromosome_names_from_genome(path)
        self._genome_chrom_sizes_cache = genome_chrom_sizes_from_fasta(path)
        return self._genome_chrom_names_cache, self._genome_chrom_sizes_cache

    def _start_blast_then_te_genome_map(
        self,
        fasta_paths: list[Path],
        map_options: TeGenomeMapOptions | None,
    ) -> None:
        """Run Blast on FASTA(s), then build and open the TE genome map when Blast finishes."""
        if not fasta_paths:
            return
        if not self._project.genome_file or not self._project.genome_file.is_file():
            QMessageBox.warning(
                self,
                "Genome required",
                f"Set the genome FASTA under {self._genome_menu_hint()} before Blast / TE genome map.",
            )
            return
        self._pending_te_genome_map = map_options
        self._pending_te_genome_map_fastas = [p.expanduser().resolve() for p in fasta_paths]
        self._set_status_message("Running Blast, then TE genome map…", 0)
        self._apply_function(
            "blast",
            preselected=self._pending_te_genome_map_fastas,
            # Keep *_b.txt so the density plot can find the hit table after Blast.
            blast_options=replace(self._blast_options, write_text=True),
        )

    def _open_blast_genome_map(
        self,
        te_input_files: list[Path] | Path,
        *,
        gene_file: Path | None = None,
        gene_files: list[Path] | None = None,
        methylation_files: list[Path] | None = None,
        te_input_kinds: list[str] | None = None,
        map_options: TeGenomeMapOptions | None = None,
    ) -> None:
        te_paths = (
            [te_input_files.expanduser().resolve()]
            if isinstance(te_input_files, Path)
            else [p.expanduser().resolve() for p in te_input_files]
        )
        gene_paths: list[Path] = []
        if gene_files:
            gene_paths.extend(p.expanduser().resolve() for p in gene_files)
        if gene_file is not None:
            gene_paths.insert(0, gene_file.expanduser().resolve())
        methylation_paths: list[Path] = []
        if methylation_files:
            methylation_paths.extend(p.expanduser().resolve() for p in methylation_files)
        if te_paths and not any(p.is_file() for p in te_paths):
            query_fasta = self._fasta_for_blast_table(te_paths[0])
            if query_fasta is not None:
                self._start_blast_then_te_genome_map([query_fasta], map_options)
                return
            QMessageBox.warning(
                self,
                "TE genome map",
                f"TE input not found:\n{te_paths[0]}\n\nRun Blast on this sequence first.",
            )
            return
        te_paths = [p for p in te_paths if p.is_file()]
        if not te_paths and not gene_paths and not methylation_paths:
            QMessageBox.warning(
                self,
                "TE genome map",
                "Provide at least one TE hit table, gene annotation, or DNA methylation file.",
            )
            return
        if not self._project.genome_file or not self._project.genome_file.is_file():
            QMessageBox.warning(
                self,
                "Genome required",
                f"Set the genome FASTA under {self._genome_menu_hint()} before building the map.",
            )
            return
        if self._genome_map_thread is not None:
            QMessageBox.information(
                self,
                "TE genome map",
                "A genome map is already being built. Please wait for it to finish.",
            )
            return
        opts = map_options or self._te_genome_map_options
        if te_paths:
            if len(te_paths) == 1:
                query_label = te_hit_table_stem(te_paths[0])
            else:
                query_label = f"{te_hit_table_stem(te_paths[0])} + {len(te_paths) - 1} more"
        elif gene_paths:
            if len(gene_paths) == 1:
                query_label = gene_paths[0].stem
            else:
                query_label = f"{gene_paths[0].stem} + {len(gene_paths) - 1} more"
        elif methylation_paths:
            if len(methylation_paths) == 1:
                query_label = methylation_paths[0].stem
            else:
                query_label = (
                    f"{methylation_paths[0].stem} + {len(methylation_paths) - 1} more"
                )
        else:
            query_label = "genome map"
        if te_paths:
            proposed = genome_map_html_path_for_te_inputs(
                te_paths,
                gene_gff=gene_paths[0] if gene_paths else None,
            )
        elif gene_paths:
            proposed = genome_map_html_path_for_te_inputs(
                [],
                gene_gff=gene_paths[0],
            )
        else:
            proposed = (
                methylation_paths[0].parent
                / f"{methylation_paths[0].stem}{GENOME_MAP_HTML_SUFFIX}"
            )
        output_html = self._confirm_plot_output_path(
            proposed,
            title="TE density plot",
        )
        if output_html is None:
            return
        self._start_genome_map_build(
            te_paths,
            gene_files=gene_paths,
            methylation_files=methylation_paths,
            te_input_kinds=te_input_kinds,
            bin_size=opts.bin_size_bp,
            density_mode=opts.density_mode,
            color_by=opts.color_by,
            methylation_coverage_min=opts.methylation_coverage_min,
            split_strand_tracks=opts.split_strand_tracks,
            collapse_overlaps=opts.collapse_overlaps,
            show_te_diff_track=opts.show_te_diff_track,
            show_meth_diff_track=opts.show_meth_diff_track,
            show_window_ref_line=opts.show_window_ref_line,
            annotation_density_mode=opts.annotation_density_mode,
            annotation_feature_types=opts.annotation_feature_types,
            annotation_feature_types_by_file=opts.annotation_feature_types_by_file,
            query_label=query_label,
            chromosomes_only=opts.chromosomes_only or None,
            superfamilies_only=None,
            chromosome_sort=opts.chromosome_sort,
            output_html=output_html,
        )

    def _confirm_plot_output_path(self, proposed: Path, *, title: str) -> Path | None:
        """If ``proposed`` already exists, ask to overwrite, rename, or cancel."""
        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        return confirm_existing_file_path(
            self,
            proposed,
            title=title,
            informative="Overwrite the existing plot, or save under a different name?",
            force_html_suffix=True,
        )

    def _start_genome_map_build(
        self,
        te_input_files: list[Path],
        *,
        gene_files: list[Path] | None = None,
        methylation_files: list[Path] | None = None,
        te_input_kinds: list[str] | None = None,
        bin_size: int,
        density_mode: str,
        color_by: str,
        methylation_coverage_min: int,
        split_strand_tracks: bool,
        collapse_overlaps: bool,
        show_te_diff_track: bool,
        show_meth_diff_track: bool,
        show_window_ref_line: bool,
        annotation_density_mode: str,
        annotation_feature_types: list[str] | None,
        annotation_feature_types_by_file: dict[str, list[str]] | None,
        query_label: str,
        chromosomes_only: list[str] | None,
        superfamilies_only: list[str] | None,
        chromosome_sort: str,
        output_html: Path | None = None,
    ) -> None:
        genome = self._project.genome_file
        assert genome is not None
        te_paths = [p.expanduser().resolve() for p in te_input_files]
        gene_paths = [p.expanduser().resolve() for p in (gene_files or [])]
        methylation_paths = [p.expanduser().resolve() for p in (methylation_files or [])]
        primary = (
            te_paths[0]
            if te_paths
            else (gene_paths[0] if gene_paths else (methylation_paths[0] if methylation_paths else Path(".")))
        )
        self._genome_map_blast_txt = primary
        self._genome_map_te_paths = list(te_paths)
        self._job_cancel.reset()
        self._job_running = True
        self._set_job_actions_enabled(False)
        workers = max(1, min(self._project.batch_workers, 8))
        worker_note = f" ({workers} workers)" if workers > 1 else ""
        self._clear_job_file_highlights()
        self._set_status_message(f"Building TE genome map{worker_note}…", 0)
        self._job_progress_epoch = getattr(self, "_job_progress_epoch", 0) + 1
        self._job_work_steps = 6
        if self._task_progress is not None:
            # +1 reserved until the HTML preview actually opens.
            self._task_progress.start(self._running_job_title("TE genome map"), 7)
        if self._stop_btn is not None:
            self._stop_btn.setEnabled(True)
        self._genome_map_thread = QThread(self)
        self._genome_map_worker = GenomeMapWorker(
            te_paths,
            genome,
            job_kwargs={
                "output_html": output_html,
                "gene_gff": gene_paths or None,
                "methylation_files": methylation_paths or None,
                "te_input_kinds": te_input_kinds,
                "bin_size": bin_size,
                "density_mode": density_mode,
                "color_by": color_by,
                "methylation_coverage_min": methylation_coverage_min,
                "split_strand_tracks": split_strand_tracks,
                "collapse_overlaps": collapse_overlaps,
                "show_te_diff_track": show_te_diff_track,
                "show_meth_diff_track": show_meth_diff_track,
                "show_window_ref_line": show_window_ref_line,
                "annotation_density_mode": annotation_density_mode,
                "annotation_feature_types": annotation_feature_types,
                "annotation_feature_types_by_file": annotation_feature_types_by_file,
                "query_label": query_label,
                "chromosomes_only": chromosomes_only,
                "superfamilies_only": superfamilies_only,
                "chromosome_sort": chromosome_sort,
                "num_workers": self._project.batch_workers,
                "cancel": self._job_cancel,
            },
        )
        self._genome_map_worker.moveToThread(self._genome_map_thread)
        self._genome_map_thread.started.connect(self._genome_map_worker.run)
        queued = Qt.ConnectionType.QueuedConnection
        self._genome_map_worker.progress.connect(self._on_genome_map_progress, queued)
        self._genome_map_worker.finished.connect(self._on_genome_map_finished, queued)
        self._genome_map_worker.failed.connect(self._on_genome_map_failed, queued)
        self._genome_map_worker.finished.connect(self._genome_map_thread.quit)
        self._genome_map_worker.failed.connect(self._genome_map_thread.quit)
        self._genome_map_thread.finished.connect(self._cleanup_genome_map_worker)
        self._genome_map_thread.start()

    @Slot(object)
    def _on_genome_map_finished(self, html_path: object) -> None:
        blast_txt = self._genome_map_blast_txt
        self._finish_background_job()
        opened = False
        try:
            if html_path is None or not Path(html_path).is_file():
                self._set_status_message(
                    "TE genome map failed", 6000, override_file_link=True
                )
                return
            path = Path(html_path).resolve()
            if blast_txt is not None:
                self._highlight_genome_map_outputs(path, blast_txt)
            if (
                self._current_folder
                and path.parent.resolve() == self._current_folder.resolve()
            ):
                self._select_path_in_file_list(path)
            self._open_genome_map_in_preview(path)
            opened = True
            n_inputs = len(self._genome_map_te_paths)
            detail = (
                "all TE samples share one y-axis scale"
                if n_inputs > 1
                else "chromosome tracks, Mbp on X axis"
            )
            self._set_status_message(
                f"Genome map opened: {path.name} ({detail})",
                12000,
                override_file_link=True,
            )
        finally:
            self._schedule_after_job_ui_ready()

        if not opened:
            QMessageBox.warning(
                self, "TE genome map", "Could not build the genome map."
            )

    @Slot(int, str)
    def _on_genome_map_progress(self, value: int, label: str) -> None:
        if self._task_progress is not None:
            work = getattr(self, "_job_work_steps", 6)
            self._task_progress.update(min(max(0, int(value)), work), label)

    @Slot(str)
    def _on_genome_map_failed(self, error: str) -> None:
        self._hide_job_progress()
        self._set_job_actions_enabled(True)
        if error == CANCELLED_MESSAGE:
            self._set_status_message("Stopped TE genome map", 6000, override_file_link=True)
            return
        if error.startswith("There is nothing to plot"):
            QMessageBox.information(self, "TE density plot", error)
            self._set_status_message("TE density plot: nothing to plot", 8000, override_file_link=True)
            return
        QMessageBox.critical(self, "TE genome map failed", error)
        self._set_status_message("TE genome map failed", 6000, override_file_link=True)

    def _cleanup_genome_map_worker(self) -> None:
        cancelled = self._job_cancel.requested
        if cancelled:
            self._job_cancel.reset()
        self._job_running = self._worker_thread is not None
        thread = self._genome_map_thread
        if thread is not None and thread.isRunning():
            thread.quit()
            thread.wait(5000)
        if self._genome_map_worker:
            self._genome_map_worker.deleteLater()
            self._genome_map_worker = None
        if self._genome_map_thread:
            self._genome_map_thread.deleteLater()
            self._genome_map_thread = None
        self._genome_map_te_paths = []
        if self._worker_thread is None:
            self._set_job_actions_enabled(True)
        if self._stop_btn is not None and not self._job_running:
            self._stop_btn.setEnabled(False)
        if cancelled:
            self._queued_genome_map_jobs.clear()
        else:
            QTimer.singleShot(0, self._start_next_genome_map_job)

    def _highlight_genome_map_outputs(
        self,
        html_path: Path,
        blast_txt: Path,
    ) -> None:
        """Highlight the new *_genome_map.html (and query FASTA) in the Files panel."""
        folder = html_path.parent.resolve()
        outputs = [html_path.resolve()]
        source = self._fasta_for_blast_table(blast_txt) or blast_txt
        sources = [source.resolve()]
        if self._current_folder and self._current_folder.resolve() == folder:
            self._populate_file_list(folder, new_outputs=outputs, source_files=sources)
        else:
            self._navigate_to_folder(folder, new_outputs=outputs, source_files=sources)

    def _fasta_for_blast_table(self, blast_txt: Path) -> Path | None:
        """Infer query FASTA from a ``*_b.txt`` path when the table is missing."""
        from repeatviewer.backends.blast_genome_plot import strip_blast_table_suffix

        name = blast_txt.name
        stem = strip_blast_table_suffix(name)
        if stem is None:
            return None
        parent = blast_txt.parent
        candidates = [
            parent / stem,
            parent / f"{stem}.fa",
            parent / f"{stem}.fasta",
            parent / f"{stem}.fas",
        ]
        seen: set[str] = set()
        for query in candidates:
            key = str(query)
            if key in seen:
                continue
            seen.add(key)
            if query.suffix.lower() in _FASTA_SUFFIXES and query.is_file():
                return query.expanduser().resolve()
        return None

    def _open_genome_map_in_preview(self, html_path: Path) -> None:
        """Open the genome map HTML in Preview after the file is written."""
        from repeatviewer.widgets.preview_tabs import FilePreviewPage

        resolved = html_path.expanduser().resolve()
        if not resolved.is_file():
            return

        def _load() -> None:
            self._ensure_preview_visible()
            if self._preview_tabs is None:
                return
            self._preview_tabs.open_file(resolved, replace=True)
            page = self._preview_tabs.currentWidget()
            if isinstance(page, FilePreviewPage):
                page.reload_from_disk()

        _load()
        QTimer.singleShot(80, _load)
        QTimer.singleShot(450, _load)

    def _run_extend(self) -> None:
        self._apply_function("extend", extend=self._extend_options)

    def _open_extend_dialog(self) -> None:
        """Click MSA Extend: parameters dialog. Run runs; Save (Options tab) stores only."""
        dialog = ExtendParamsDialog(
            initial=self._extend_options,
            pan_genome_names=self._pan_genome_names(),
            project=self._project,
            project_genome=self._genome_path_resolved(),
            genome_history=load_genome_history(),
            parent=self,
        )
        try:
            if not dialog.exec():
                return
            self._extend_options = dialog.options()
            save_extend_options(self._extend_options)
            self._persist_workspace()
            run_after = dialog.should_run()
        finally:
            QTimer.singleShot(0, self._repair_browser_chrome)
        if run_after:
            self._run_extend()
        else:
            self._set_status_message("MSA Extend settings saved", 4000)

    def _quick_extend_right_click(self) -> None:
        """Right-click MSA Extend: run with saved settings, or open the dialog."""
        if self._extend_options.right_click_action == EXTEND_RIGHT_CLICK_DIALOG:
            self._open_extend_dialog()
        else:
            self._run_extend()

    def _run_tsd_search(self) -> None:
        self._apply_function("tsd_search", tsd_search_options=self._tsd_search_options)

    def _run_teaid(self) -> None:
        self._apply_function("teaid", teaid_options=self._teaid_options)

    def _open_teaid_dialog(self) -> None:
        """Click TE-Aid: parameters dialog. Run runs; Save (Options tab) stores only."""
        dialog = TeaAidParamsDialog(
            initial=self._teaid_options,
            project=self._project,
            project_genome=self._genome_path_resolved(),
            genome_history=load_genome_history(),
            parent=self,
        )
        try:
            if not dialog.exec():
                return
            self._teaid_options = dialog.options()
            save_teaid_options(self._teaid_options)
            self._persist_workspace()
            run_after = dialog.should_run()
        finally:
            QTimer.singleShot(0, self._repair_browser_chrome)
        if run_after:
            self._run_teaid()
        else:
            self._set_status_message("TE-Aid settings saved", 4000)

    def _quick_teaid_right_click(self) -> None:
        """Right-click TE-Aid: run with saved settings, or open the dialog."""
        if self._teaid_options.right_click_action == TEAID_RIGHT_CLICK_DIALOG:
            self._open_teaid_dialog()
        else:
            self._run_teaid()

    def _modify_teaid_parameters(self) -> None:
        """Menu: same Cancel / Save / Run dialog as the TE-Aid toolbar button."""
        self._open_teaid_dialog()

    def _modify_extend_parameters(self) -> None:
        """Menu: same Cancel / Save / Run dialog as the MSA Extend toolbar button."""
        self._open_extend_dialog()

    def _modify_tsd_search_parameters(self) -> None:
        """Menu / change parameters: Cancel / Save / Run."""
        dialog = TsdSearchParamsDialog(
            initial=self._tsd_search_options,
            parent=self,
        )
        try:
            if not dialog.exec():
                return
            self._tsd_search_options = dialog.options()
            save_tsd_search_options(self._tsd_search_options)
            run_after = dialog.should_run()
        finally:
            QTimer.singleShot(0, self._repair_browser_chrome)
        if run_after:
            self._run_tsd_search()
        else:
            self._set_status_message(
                "TSD search defaults: "
                f"left {self._tsd_search_options.left_boundary}, "
                f"right {self._tsd_search_options.right_boundary}, "
                f"radius {self._tsd_search_options.nt_radius} nt",
                5000,
            )

    def _show_gff_features_from_selection(self) -> None:
        selected = self._selected_files()
        fasta_files = [p for p in selected if p.suffix.lower() in _FASTA_SUFFIXES]
        if len(fasta_files) != 1:
            QMessageBox.information(
                self,
                "GFF annotations",
                "Select exactly one FASTA sequence file in the Files panel.",
            )
            return
        self._show_gff_features(fasta_files[0])

    def _show_gff_features(self, fasta_path: Path) -> None:
        if fasta_path.suffix.lower() not in _FASTA_SUFFIXES:
            QMessageBox.warning(
                self,
                "GFF annotations",
                "Choose a FASTA sequence file.",
            )
            return
        gff_path = find_companion_gff(fasta_path)
        if gff_path is None:
            start_dir = str(fasta_path.parent)
            picked, _ = QFileDialog.getOpenFileName(
                self,
                "Select GFF/GTF annotation file",
                start_dir,
                "Annotation (*.gff *.gff3 *.gtf);;All files (*)",
            )
            if not picked:
                return
            gff_path = Path(picked)
        if not gff_path.is_file():
            QMessageBox.warning(self, "GFF annotations", f"File not found:\n{gff_path}")
            return
        try:
            html_text = build_gff_annotation_html(fasta_path, gff_path)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "GFF annotations", str(exc))
            return
        temp_dir = Path(tempfile.gettempdir()) / "repeatviewer"
        temp_dir.mkdir(parents=True, exist_ok=True)
        html_path = temp_dir / f"gff_{fasta_path.stem}_{gff_path.stem}.html"
        html_path.write_text(html_text, encoding="utf-8")
        self._open_files_in_preview_tabs([html_path])
        self._set_status_message(
            f"GFF map: {fasta_path.name} + {gff_path.name}",
            6000,
        )

    def _extend_with_different_distances(self) -> None:
        if not self._project.has_genome_configured():
            QMessageBox.warning(
                self,
                "Genome required",
                f"Define genome FASTA under {self._genome_menu_hint()} before running MSA Extend.",
            )
            return
        seq_files = [
            p for p in self._selected_files() if p.suffix.lower() in _FASTA_SUFFIXES
        ]
        if not seq_files:
            QMessageBox.information(
                self,
                "No FASTA selected",
                "Select one or more sequence files in the Files panel.",
            )
            return
        dialog = ExtendParamsDialog(
            initial=self._extend_options,
            pan_genome_names=self._pan_genome_names(),
            project=self._project,
            project_genome=self._genome_path_resolved(),
            genome_history=load_genome_history(),
            parent=self,
        )
        if not dialog.exec() or not dialog.should_run():
            return
        opts = dialog.options()
        answer = QMessageBox.question(
            self,
            "Confirm MSA Extend",
            f"Extend {len(seq_files)} sequence file(s) by {opts.left} bp (left) "
            f"and {opts.right} bp (right)?\n\n"
            f"Output names will include _L{opts.left}_R{opts.right}.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        self._apply_function("extend", preselected=seq_files, extend=opts)

    def _run_consensus(self) -> None:
        self._apply_function(
            "consensus_generation",
            consensus_options=self._consensus_options,
        )

    def _open_consensus_dialog(self) -> None:
        """Click Consensus: parameters dialog. Run runs; Save (Options tab) stores only."""
        dialog = ConsensusParamsDialog(
            initial=self._consensus_options,
            parent=self,
        )
        try:
            if not dialog.exec():
                return
            self._consensus_options = dialog.options()
            save_consensus_options(self._consensus_options)
            self._persist_workspace()
            run_after = dialog.should_run()
        finally:
            QTimer.singleShot(0, self._repair_browser_chrome)
        if run_after:
            self._run_consensus()
        else:
            self._set_status_message("Consensus settings saved", 4000)

    def _quick_consensus_right_click(self) -> None:
        """Right-click Consensus: run with saved settings, or open the dialog."""
        if self._consensus_options.right_click_action == CONSENSUS_RIGHT_CLICK_DIALOG:
            self._open_consensus_dialog()
        else:
            self._run_consensus()

    def _modify_consensus_parameters(self) -> None:
        """Menu: same Cancel / Save / Run dialog as the Consensus toolbar button."""
        self._open_consensus_dialog()

    def _run_flexidot_combined(self) -> None:
        self._apply_function(
            "flexidot",
            flexidot_combined=True,
            flexidot_options=self._flexidot_options,
        )

    def _run_flexidot_saved(self) -> None:
        """Run Flexidot with saved settings (combined or per-file)."""
        self._apply_function(
            "flexidot",
            flexidot_combined=self._flexidot_options.combined,
            flexidot_options=self._flexidot_options,
        )

    def _open_flexidot_dialog(self) -> None:
        """Click Flexidot: parameters dialog. Run runs; Save (Options tab) stores only."""
        dialog = FlexidotParamsDialog(
            initial=self._flexidot_options,
            parent=self,
        )
        try:
            if not dialog.exec():
                return
            self._flexidot_options = dialog.options()
            save_flexidot_options(self._flexidot_options)
            self._persist_workspace()
            run_after = dialog.should_run()
        finally:
            QTimer.singleShot(0, self._repair_browser_chrome)
        if run_after:
            self._run_flexidot_saved()
        else:
            self._set_status_message("Flexidot settings saved", 4000)

    def _quick_flexidot_right_click(self) -> None:
        """Right-click Flexidot: run with saved settings, or open the dialog."""
        if self._flexidot_options.right_click_action == FLEXIDOT_RIGHT_CLICK_DIALOG:
            self._open_flexidot_dialog()
        else:
            self._run_flexidot_saved()

    def _reload_shared_annotation_settings(self) -> None:
        """Reload TE/annotation settings into every open sequence-map panel."""
        from repeatviewer.widgets.sequence_map_viewer import SequenceMapPanel

        for panel in self.findChildren(SequenceMapPanel):
            panel.reload_shared_annotation_settings()

    def _te_annotation_fasta_selection(self) -> list[Path] | None:
        """Validate selected single-sequence FASTAs for TE annotation, or None."""
        from repeatviewer.backends.te_annotation_batch import is_single_sequence_fasta

        selected = self._selected_files()
        if not selected:
            QMessageBox.information(
                self,
                "No sequence selected",
                "Select one or more single-sequence FASTA files in the Files panel.",
            )
            return None
        seq_files = [p for p in selected if p.suffix.lower() in _FASTA_SUFFIXES]
        if not seq_files:
            QMessageBox.warning(
                self,
                "Invalid selection",
                "Selected files do not include FASTA sequence files.",
            )
            return None
        single = [p for p in seq_files if is_single_sequence_fasta(p)]
        skipped = [p for p in seq_files if p not in single]
        if skipped:
            skipped_names = "\n".join(p.name for p in skipped[:12])
            if len(skipped) > 12:
                skipped_names += "\n..."
            QMessageBox.warning(
                self,
                "Skipped files",
                f"{len(skipped)} file(s) skipped — only single-sequence FASTA files "
                f"are supported:\n\n{skipped_names}",
            )
        if not single:
            QMessageBox.warning(
                self,
                "No valid files",
                "Select at least one single-sequence FASTA file for Sequence Annotation.",
            )
            return None
        return single

    def _confirm_te_annotation_overwrite(self, fastas: list[Path]) -> bool:
        """Ask before overwriting existing ``*.RV.gff3`` companions. True = continue."""
        from repeatviewer.backends.te_annotation_batch import existing_annotation_gff_paths

        existing = existing_annotation_gff_paths(fastas)
        if not existing:
            return True
        listing = "\n".join(f"  • {path.name}" for path in existing[:20])
        if len(existing) > 20:
            listing += f"\n  … and {len(existing) - 20} more"
        answer = QMessageBox.question(
            self,
            "Overwrite annotations?",
            "Annotation file(s) already exist for the selected sequence(s):\n\n"
            f"{listing}\n\n"
            "Running Sequence Annotation will overwrite them. Continue?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        return answer == QMessageBox.StandardButton.Yes

    def _run_te_annotation_saved(self) -> None:
        single = self._te_annotation_fasta_selection()
        if single is None:
            return
        opts = batch_options_from_te_annotation_job(self._te_annotation_options)
        if not opts.auto_modes and not opts.similarity_databases:
            QMessageBox.information(
                self,
                "Nothing selected",
                "No annotation modules are saved. Open Sequence Annotation and choose "
                "at least one auto-annotation module or similarity database.",
            )
            return
        if not self._confirm_te_annotation_overwrite(single):
            return
        self._apply_function(
            "te_annotation",
            preselected=single,
            te_annotation_options=opts,
        )

    def _open_te_annotation_dialog(self) -> None:
        """Click TE annotation: dialog. Run runs; Save (Options tab) stores only."""
        dialog = TeAnnotationDialog(
            initial=self._te_annotation_options,
            parent=self,
        )
        run_after = False
        opts = None
        try:
            if not exec_modal_dialog(dialog):
                return
            dialog.apply_settings_panels()
            self._te_annotation_options = dialog.job_options()
            save_te_annotation_options(self._te_annotation_options)
            self._persist_workspace()
            self._reload_shared_annotation_settings()
            run_after = dialog.should_run()
            opts = dialog.options()
        finally:
            QTimer.singleShot(0, self._repair_browser_chrome)
        if not run_after:
            self._set_status_message("Sequence Annotation settings saved", 4000)
            return
        if opts is None:
            return
        if not opts.auto_modes and not opts.similarity_databases:
            QMessageBox.information(
                self,
                "Nothing selected",
                "Select at least one auto-annotation module or similarity database.",
            )
            return
        single = self._te_annotation_fasta_selection()
        if single is None:
            return
        if not self._confirm_te_annotation_overwrite(single):
            return
        self._apply_function(
            "te_annotation",
            preselected=single,
            te_annotation_options=opts,
        )

    def _quick_te_annotation_right_click(self) -> None:
        """Right-click TE annotation: run with saved selection, or open the dialog."""
        if self._te_annotation_options.right_click_action == TE_ANNOTATION_RIGHT_CLICK_DIALOG:
            self._open_te_annotation_dialog()
        else:
            self._run_te_annotation_saved()

    def _run_orf_prediction(self) -> None:
        self._apply_function(
            "orf_prediction",
            orf_prediction_settings=self._orf_prediction_settings,
        )

    def _open_orf_prediction_dialog(self) -> None:
        """Toolbar left-click: settings with Run."""
        from repeatviewer.backends.auto_annotation_settings import (
            save_orf_prediction_settings,
        )

        dialog = OrfPredictionParamsDialog(
            initial=self._orf_prediction_settings,
            for_run=True,
            parent=self,
        )
        try:
            if not exec_modal_dialog(dialog):
                return
            self._orf_prediction_settings = dialog.options()
            save_orf_prediction_settings(self._orf_prediction_settings)
            self._persist_workspace()
            run_after = dialog.should_run()
        finally:
            QTimer.singleShot(0, self._repair_browser_chrome)
        if run_after:
            self._run_orf_prediction()
        else:
            self._set_status_message("ORF Prediction settings saved", 4000)

    def _modify_orf_prediction_parameters(self) -> None:
        """Menu: same Cancel / Save / Run dialog as the ORF Prediction toolbar button."""
        self._open_orf_prediction_dialog()

    def _run_te_annotation(self) -> None:
        """Tools menu default entry: run with saved selection."""
        self._run_te_annotation_saved()

    def _modify_te_annotation_parameters(self) -> None:
        """Menu: same Cancel / Save / Run dialog as the TE annotation toolbar button."""
        self._open_te_annotation_dialog()

    def _modify_flexidot_parameters(self) -> None:
        """Menu: same Cancel / Save / Run dialog as the Flexidot toolbar button."""
        self._open_flexidot_dialog()

    def _add_msa_clean_menu_action(
        self,
        menu: QMenu,
        label: str,
        method_id: str,
    ) -> None:
        action = menu.addAction(label, lambda: self._apply_msa_clean(method_id))
        desc = msa_clean_method_description(method_id)
        if desc:
            action.setToolTip(desc)
            action.setStatusTip(desc)

    def _modify_msa_clean_parameters(self) -> None:
        """Menu: Cancel / Save settings dialog (no Run — clean runs from the tool menu)."""
        self._open_msa_clean_dialog()

    def _open_msa_clean_dialog(self) -> None:
        """Open MSA Clean parameters (from Options right-click or Change menu)."""
        dialog = MsaCleanParamsDialog(
            initial=self._msa_clean_options,
            parent=self,
        )
        try:
            if not dialog.exec():
                return
            self._msa_clean_options = dialog.options()
            save_msa_clean_options(self._msa_clean_options)
            self._persist_workspace()
        finally:
            QTimer.singleShot(0, self._repair_browser_chrome)
        self._set_status_message("MSA clean settings saved", 4000)

    def _quick_msa_clean_right_click(self) -> None:
        """Right-click MSA Clean: run the action chosen in the Options tab."""
        action = self._msa_clean_options.right_click_action
        if action == MSA_CLEAN_RIGHT_CLICK_DIALOG:
            self._open_msa_clean_dialog()
        else:
            self._apply_msa_clean(action)

    def _run_default_msa_clean(self) -> None:
        self._apply_msa_clean("msa_clean")

    def _apply_msa_clean(self, func_name: str) -> None:
        self._apply_function(
            func_name,
            msa_clean_options=self._msa_clean_options,
        )

    def _ensure_preview_dock_attached(self) -> None:
        if self._dock_preview is None or self._preview_host is None:
            return
        if self._dock_preview.widget() is not self._preview_host:
            self._dock_preview.setWidget(self._preview_host)
        self._dock_preview.show()
        self._preview_host.show()

    def _ensure_preview_visible(self) -> None:
        self._ensure_preview_dock_attached()
        if self._dock_preview is not None:
            self._dock_preview.raise_()
        if self._preview_tabs is not None:
            self._preview_tabs.show()

    def _open_preview_tab(self, path: Path, *, replace: bool = False) -> None:
        if self._preview_tabs is None:
            return
        try:
            self._preview_tabs.open_file(path, replace=replace)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(
                self,
                "Preview failed",
                f"Could not open {path.name} in preview:\n{exc}",
            )
            return
        self._ensure_preview_visible()
        QTimer.singleShot(0, self._ensure_preview_visible)

    def _prompt_rename_text(
        self,
        *,
        title: str,
        label: str,
        initial: str,
    ) -> tuple[str, bool]:
        """Wide rename prompt so long filenames fit without horizontal scrolling."""
        from repeatviewer.widgets.file_conflict import prompt_rename_text

        return prompt_rename_text(
            self, title=title, label=label, initial=initial
        )

    # -- cytogenetic images -------------------------------------------------

    def _attach_cyto_images(self, fasta: Path, images: list) -> None:
        from repeatviewer.widgets.cyto_images_ui import attach_images

        saved = attach_images(self, Path(fasta), [Path(p) for p in images or []])
        if saved:
            noun = "image" if saved == 1 else "images"
            self._set_status_message(f"Attached {saved} cytogenetic {noun} to {Path(fasta).name}", 5000)

    def _open_cyto_images(self, fasta: Path) -> None:
        from repeatviewer.widgets.cyto_images_ui import open_images_for_sequence

        open_images_for_sequence(self, Path(fasta))

    def _on_cyto_images_changed(self, folder: str) -> None:
        current = getattr(self, "_current_folder", None)
        if current is not None and str(Path(current).expanduser().resolve()) == folder:
            self._file_list.refresh_cyto_images()

    def _rename_file(self, path: Path) -> None:
        from repeatviewer.backends.sequence_features import (
            companion_annotation_renames_for_fasta,
        )

        path = path.expanduser().resolve()
        if not path.is_file():
            return
        new_name, ok = self._prompt_rename_text(
            title="Rename file",
            label="New file name:",
            initial=path.name,
        )
        if not ok or not new_name.strip() or new_name.strip() == path.name:
            return
        target = (path.parent / new_name.strip()).expanduser().resolve()
        if target.exists():
            QMessageBox.warning(self, "Rename failed", f"Already exists:\n{target}")
            return

        companion_renames: list[tuple[Path, Path]] = []
        if path.suffix.lower() in _FASTA_SUFFIXES:
            companion_renames = companion_annotation_renames_for_fasta(path, target)
            blocked = [dst for _src, dst in companion_renames if dst.exists()]
            if blocked:
                QMessageBox.warning(
                    self,
                    "Rename failed",
                    "Cannot rename paired annotation file(s); already exists:\n"
                    + "\n".join(str(path) for path in blocked[:6]),
                )
                return

        try:
            path.rename(target)
            for src, dst in companion_renames:
                src.rename(dst)
        except OSError as exc:
            QMessageBox.critical(self, "Rename failed", str(exc))
            return
        if path.suffix.lower() in _FASTA_SUFFIXES and target.parent == path.parent:
            from repeatviewer.backends.cyto_images import rename_sequence

            try:
                rename_sequence(path.parent, path.name, target.name)
            except OSError:
                logger.exception("Could not move cytogenetic image links to %s", target.name)

        from repeatviewer.widgets.file_colors import replace_file_color_path

        replace_file_color_path(path, target)
        for src, dst in companion_renames:
            replace_file_color_path(src, dst)

        closed = [path, *[src for src, _dst in companion_renames]]
        if self._preview_tabs is not None:
            self._preview_tabs.close_paths(closed)
        if self._current_folder:
            self._populate_file_list(
                self._current_folder,
                new_outputs=[target],
                source_files=[dst for _src, dst in companion_renames],
                open_preview=True,
            )
        else:
            self._select_path_in_file_list(target)
        if companion_renames:
            extras = ", ".join(dst.name for _src, dst in companion_renames)
            self._set_status_message(
                f"Renamed to {target.name} (also {extras})",
                6000,
            )
        else:
            self._set_status_message(f"Renamed to {target.name}", 5000)

    def _rename_folder(self, folder: Path) -> None:
        folder = folder.expanduser().resolve()
        if not folder.is_dir():
            return
        new_name, ok = self._prompt_rename_text(
            title="Rename folder",
            label="New folder name:",
            initial=folder.name,
        )
        if not ok or not new_name.strip() or new_name.strip() == folder.name:
            return
        target = folder.parent / new_name.strip()
        if target.exists():
            QMessageBox.warning(self, "Rename failed", f"Already exists:\n{target}")
            return
        try:
            folder.rename(target)
        except OSError as exc:
            QMessageBox.critical(self, "Rename failed", str(exc))
            return
        self._replace_folder_path_in_tree(folder, target)
        replace_folder_color_paths(folder, target)
        apply_folder_colors_to_tree(self._folder_tree)
        self._hidden_root_folders = replace_hidden_root_path(
            self._hidden_root_folders, folder, target
        )
        if self._current_folder and self._current_folder.resolve() == folder:
            self._current_folder = target
            self._on_folder_selected(target)
        self._set_status_message(f"Renamed folder to {target.name}", 5000)

    def _clear_root_folders(self, folders: list[Path]) -> None:
        if not folders:
            return
        removed_folders: list[Path] = []
        for folder in folders:
            if self._clear_root_folder(folder, show_status=False, record_undo=False):
                removed_folders.append(folder.expanduser().resolve())
        if removed_folders:
            self._push_browser_undo(
                {"kind": "hide_folders", "folders": removed_folders}
            )
            label = "folder" if len(removed_folders) == 1 else "folders"
            self._set_status_message(
                f"Removed {len(removed_folders)} {label} from list "
                f"(⌘Z / Ctrl+Z to restore)",
                6000,
            )
            self._persist_workspace()

    def _clear_root_folder(
        self,
        folder: Path,
        *,
        show_status: bool = True,
        record_undo: bool = True,
    ) -> bool:
        """Hide a top-level root folder from the tree without deleting it on disk."""
        folder = folder.expanduser().resolve()
        is_top_level = False
        for index in range(self._folder_tree.topLevelItemCount()):
            item = self._folder_tree.topLevelItem(index)
            raw = item.data(0, ROLE_PATH)
            if raw and Path(str(raw)).resolve() == folder:
                is_top_level = True
                break
        if not is_top_level:
            return False

        self._hidden_root_folders = hide_root_folder(self._hidden_root_folders, folder)
        self._remove_folder_from_tree(folder)
        self._navigate_after_folder_removed(folder)
        if record_undo:
            self._push_browser_undo({"kind": "hide_folders", "folders": [folder]})

        if show_status:
            self._set_status_message(
                f"Removed from folder list (⌘Z / Ctrl+Z to restore): {folder.name}",
                6000,
            )
            self._persist_workspace()
        return True

    def _restore_removed_folders(self) -> None:
        hidden = sorted(self._hidden_root_folders)
        if not hidden:
            QMessageBox.information(
                self,
                "Restore deleted files",
                "No top-level folders were removed from the list.",
            )
            return

        dialog = QDialog(self)
        dialog.setWindowTitle("Restore deleted files")
        layout = QVBoxLayout(dialog)
        layout.addWidget(
            QLabel(
                "These top-level folders were hidden from the Folders panel. "
                "Restoring does not change files on disk."
            )
        )
        list_widget = QListWidget()
        for path_str in hidden:
            path = Path(path_str)
            item = QListWidgetItem(f"{path.name}  ({path_str})")
            item.setData(ROLE_PATH, path_str)
            if not path.is_dir():
                item.setForeground(QBrush(QColor("#888888")))
            list_widget.addItem(item)
        layout.addWidget(list_widget)

        buttons = QDialogButtonBox()
        restore_btn = buttons.addButton("Restore selected", QDialogButtonBox.ButtonRole.AcceptRole)
        restore_all_btn = buttons.addButton("Restore all", QDialogButtonBox.ButtonRole.ActionRole)
        close_btn = buttons.addButton(QDialogButtonBox.StandardButton.Close)
        layout.addWidget(buttons)

        def restore_paths(paths: list[Path]) -> None:
            restored = 0
            for path in paths:
                if not path.is_dir():
                    continue
                self._add_root_folder(path, unhide=True)
                restored += 1
            if restored:
                self._set_status_message(f"Restored {restored} folder(s) to the list", 5000)
            dialog.accept()

        def restore_selected() -> None:
            paths = []
            for item in list_widget.selectedItems():
                raw = item.data(ROLE_PATH)
                if raw:
                    paths.append(Path(str(raw)))
            if not paths:
                QMessageBox.information(dialog, "Restore", "Select one or more folders in the list.")
                return
            restore_paths(paths)

        def restore_all() -> None:
            restore_paths([Path(str(p)) for p in hidden])

        restore_btn.clicked.connect(restore_selected)
        restore_all_btn.clicked.connect(restore_all)
        close_btn.clicked.connect(dialog.reject)
        list_widget.itemDoubleClicked.connect(restore_selected)

        dialog.exec()

    def _dedupe_folder_selection(self, folders: list[Path]) -> list[Path]:
        """Drop nested folders when an ancestor is also selected."""
        resolved = [folder.expanduser().resolve() for folder in folders]
        unique = list(dict.fromkeys(resolved))
        return [
            folder
            for folder in unique
            if not any(
                folder != other and folder in other.parents for other in unique
            )
        ]

    def _navigate_after_folder_removed(self, removed: Path) -> None:
        if not self._current_folder:
            return
        removed = removed.expanduser().resolve()
        current = self._current_folder.expanduser().resolve()
        if current != removed and removed not in current.parents:
            return
        parent = removed.parent
        if parent.is_dir() and self._find_tree_item_for_path(parent) is not None:
            self._navigate_to_folder(parent)
        elif self._folder_tree.topLevelItemCount():
            raw = self._folder_tree.topLevelItem(0).data(0, ROLE_PATH)
            if raw:
                self._navigate_to_folder(Path(str(raw)))
        else:
            self._current_folder = None
            self._folder_tree.set_active_folder(None)
            self._cached_folder_files = []
            self._file_list.clear_rows()

    def _delete_folders(self, folders: list[Path]) -> None:
        targets = self._dedupe_folder_selection(folders)
        if not targets:
            return
        if len(targets) == 1:
            folder = targets[0]
            answer = QMessageBox.question(
                self,
                "Delete folder",
                f"Delete this folder and everything inside it?\n\n{folder}\n\n"
                "Press ⌘Z / Ctrl+Z afterward to restore.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
        else:
            names = "\n".join(str(folder) for folder in targets[:8])
            if len(targets) > 8:
                names += f"\n… and {len(targets) - 8} more"
            answer = QMessageBox.question(
                self,
                "Delete folders",
                f"Delete {len(targets)} folders and everything inside them?\n\n{names}\n\n"
                "Press ⌘Z / Ctrl+Z afterward to restore.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No,
            )
        if answer != QMessageBox.StandardButton.Yes:
            return

        batch_dir = self._new_trash_batch_dir()
        if batch_dir is None:
            QMessageBox.critical(
                self,
                "Delete folder failed",
                "Could not create a temporary trash folder for undo.",
            )
            return

        deleted = 0
        failed: list[str] = []
        batch: list[tuple[Path, Path]] = []
        for index, folder in enumerate(targets):
            try:
                resolved = folder.expanduser().resolve()
            except OSError as exc:
                failed.append(f"{folder.name}: {exc}")
                continue
            if not resolved.is_dir():
                continue
            dest = batch_dir / f"{index:04d}_{resolved.name}"
            try:
                shutil.move(str(resolved), str(dest))
            except OSError as exc:
                failed.append(f"{resolved.name}: {exc}")
                continue

            deleted += 1
            batch.append((resolved, dest))
            if self._preview_tabs is not None:
                self._preview_tabs.close_under_folder(resolved)
            self._remove_folder_from_tree(resolved)
            remove_folder_colors_under(resolved)
            self._hidden_root_folders = remove_hidden_root(
                self._hidden_root_folders, resolved
            )
            self._navigate_after_folder_removed(resolved)

        if batch:
            self._push_browser_undo({"kind": "delete_folders", "items": batch})
        else:
            try:
                batch_dir.rmdir()
            except OSError:
                pass

        if failed:
            QMessageBox.warning(
                self,
                "Delete partially failed",
                f"Deleted {deleted} folder(s).\n\nFailed:\n" + "\n".join(failed[:8]),
            )
        elif deleted:
            label = "folder" if deleted == 1 else "folders"
            self._set_status_message(
                f"Deleted {deleted} {label} (⌘Z / Ctrl+Z to restore)",
                6000,
            )
        elif failed:
            QMessageBox.critical(
                self,
                "Delete folder failed",
                "Failed:\n" + "\n".join(failed[:8]),
            )

    def _delete_folder(self, folder: Path) -> None:
        self._delete_folders([folder])

    def _remove_folder_from_tree(self, folder: Path) -> None:
        item = self._find_tree_item_for_path(folder.expanduser().resolve())
        if item is None:
            return
        parent = item.parent()
        if parent is not None:
            parent.removeChild(item)
        else:
            index = self._folder_tree.indexOfTopLevelItem(item)
            if index >= 0:
                self._folder_tree.takeTopLevelItem(index)
        self._folder_tree.refresh_horizontal_extent()

    def _replace_folder_path_in_tree(self, old: Path, new: Path) -> None:
        old_root = old.expanduser().resolve()
        new_root = new.expanduser().resolve()

        def walk(item: QTreeWidgetItem) -> None:
            raw = item.data(0, ROLE_PATH)
            if raw:
                path = Path(str(raw)).expanduser().resolve()
                if path == old_root or old_root in path.parents:
                    try:
                        rel = path.relative_to(old_root)
                        updated = new_root / rel if str(rel) != "." else new_root
                    except ValueError:
                        updated = path
                    item.setData(0, ROLE_PATH, str(updated))
                    self._set_folder_tree_item_label(item, updated)
                    item.setToolTip(0, str(updated))
            for index in range(item.childCount()):
                walk(item.child(index))

        for index in range(self._folder_tree.topLevelItemCount()):
            walk(self._folder_tree.topLevelItem(index))
        self._folder_tree.refresh_horizontal_extent()

    def _separate_fasta_file(self, path: Path, *, subfolder_per_sequence: bool = False) -> None:
        path = path.expanduser().resolve()
        if path.suffix.lower() not in _FASTA_SUFFIXES:
            QMessageBox.warning(self, "Not a FASTA file", "Select a .fa / .fasta file.")
            return
        try:
            seq_count = count_fasta_sequences(path)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Read FASTA failed", str(exc))
            return
        if seq_count <= 1:
            QMessageBox.information(
                self,
                "Separate FASTA",
                "This file already contains only one sequence.",
            )
            return

        output_dir = separation_output_dir(path)
        mapping_note = (
            f"Sequence_name_mapping_{sanitize_name(path.name)}.txt will be written next to the input file."
        )
        if subfolder_per_sequence:
            title = "Separate into subfolders"
            layout_note = (
                f"Each sequence will be written to its own subfolder under:\n{output_dir}\n"
                f"(e.g. {output_dir.name}/sequence_name/sequence_name{path.suffix.lower() or '.fa'})\n\n"
            )
            done_title = "Separate into subfolders completed"
        else:
            title = "Separate FASTA sequences"
            layout_note = f"into separate files under:\n{output_dir}\n\n"
            done_title = "Separate FASTA completed"
        answer = QMessageBox.question(
            self,
            title,
            f"Split {seq_count} sequences from\n{path.name}\n\n"
            f"{layout_note}"
            f"Special characters in headers will be sanitized. {mapping_note}\n\n"
            "Continue?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return

        try:
            created, written = separate_sequences(
                path,
                output_dir,
                subfolder_per_sequence=subfolder_per_sequence,
            )
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Separate FASTA failed", str(exc))
            return

        # New folder creation also triggers the FS watcher; suppress it briefly so
        # a rebuild cannot steal selection to the next sibling folder.
        self._watcher_refresh_suppressed_until = time.monotonic() + 3.0
        self._refresh_subfolders_under(path.parent, keep_folder=output_dir)
        self._navigate_to_folder(
            output_dir,
            new_outputs=created,
            source_files=[path],
        )
        if subfolder_per_sequence:
            status = f"Created {written} subfolder(s) in {output_dir.name}/"
            done_msg = f"Wrote {written} sequence(s), each in its own subfolder under:\n{output_dir}"
        else:
            status = f"Created {written} single-sequence file(s) in {output_dir.name}/"
            done_msg = f"Wrote {written} file(s) to:\n{output_dir}"
        self._set_status_message(status, 10000)
        QMessageBox.information(self, done_title, done_msg)

    def _separate_fasta_to_subfolders_file(self, path: Path) -> None:
        self._separate_fasta_file(path, subfolder_per_sequence=True)

    def _duplicate_selected_files(self, paths: list[Path]) -> None:
        if not paths:
            return
        created: list[Path] = []
        failed: list[str] = []
        for path in paths:
            src = path.expanduser().resolve()
            if not src.is_file():
                continue
            dest = _unique_duplicate_path(src)
            try:
                shutil.copy2(src, dest)
                created.append(dest)
            except OSError as exc:
                failed.append(f"{src.name}: {exc}")

        if self._current_folder:
            current = self._current_folder.resolve()
            self._populate_file_list(
                self._current_folder,
                new_outputs=[p for p in created if p.parent.resolve() == current],
                source_files=[
                    p.resolve() for p in paths if p.expanduser().resolve().parent == current
                ],
            )

        if failed:
            QMessageBox.warning(
                self,
                "Duplicate partially failed",
                f"Created {len(created)} copy/copies.\n\nFailed:\n" + "\n".join(failed[:8]),
            )
        elif created:
            self._set_status_message(
                f"Duplicated {len(created)} file(s)",
                5000,
            )

    def _push_browser_undo(self, entry: dict[str, Any]) -> None:
        self._browser_undo_stack.append(entry)
        while len(self._browser_undo_stack) > 20:
            old = self._browser_undo_stack.pop(0)
            self._discard_browser_undo_entry(old)

    def _discard_browser_undo_entry(self, entry: dict[str, Any]) -> None:
        kind = str(entry.get("kind") or "")
        if kind in {"files", "delete_folders"}:
            items = list(entry.get("items") or [])
            parents: set[Path] = set()
            for _original, trash in items:
                trash_path = Path(trash)
                parents.add(trash_path.parent)
                try:
                    if trash_path.is_file() or trash_path.is_symlink():
                        trash_path.unlink()
                    elif trash_path.is_dir():
                        shutil.rmtree(trash_path, ignore_errors=True)
                except OSError:
                    pass
            for parent in parents:
                try:
                    if parent.is_dir() and not any(parent.iterdir()):
                        parent.rmdir()
                except OSError:
                    pass

    def _new_trash_batch_dir(self) -> Path | None:
        try:
            self._deleted_files_trash_root.mkdir(parents=True, exist_ok=True)
            batch_dir = self._deleted_files_trash_root / uuid.uuid4().hex
            batch_dir.mkdir(parents=True, exist_ok=True)
            return batch_dir
        except OSError:
            return None

    def _expand_paths_for_delete(self, paths: list[Path]) -> list[Path]:
        """Include paired annotation sidecars (e.g. seqmap JSON with RV.gff3)."""
        expanded: list[Path] = []
        seen: set[str] = set()
        for path in paths:
            try:
                resolved = path.expanduser().resolve()
            except OSError:
                continue
            candidates = [resolved, *companion_annotation_files_for_delete(resolved)]
            for candidate in candidates:
                key = str(candidate)
                if key in seen:
                    continue
                seen.add(key)
                if candidate.is_file():
                    expanded.append(candidate)
        return expanded

    def _soft_delete_files(self, paths: list[Path]) -> tuple[int, list[str]]:
        """Move files into a session trash folder so Cmd+Z can restore them."""
        if not paths:
            return 0, []
        batch_dir = self._new_trash_batch_dir()
        if batch_dir is None:
            return 0, ["Could not create trash folder"]

        batch: list[tuple[Path, Path]] = []
        failed: list[str] = []
        for index, path in enumerate(paths):
            try:
                resolved = path.expanduser().resolve()
            except OSError as exc:
                failed.append(f"{path.name}: {exc}")
                continue
            if not resolved.is_file():
                failed.append(f"{path.name}: not a file")
                continue
            dest = batch_dir / f"{index:04d}_{resolved.name}"
            try:
                shutil.move(str(resolved), str(dest))
            except OSError as exc:
                failed.append(f"{resolved.name}: {exc}")
                continue
            batch.append((resolved, dest))

        if batch:
            self._push_browser_undo({"kind": "files", "items": batch})
        else:
            try:
                batch_dir.rmdir()
            except OSError:
                pass
        return len(batch), failed

    def _undo_browser_action(self) -> None:
        if not self._browser_undo_stack:
            self._set_status_message("Nothing to undo", 3000)
            return
        entry = self._browser_undo_stack.pop()
        kind = str(entry.get("kind") or "")
        if kind == "files":
            self._undo_deleted_files_batch(list(entry.get("items") or []))
        elif kind == "hide_folders":
            self._undo_hidden_folders(list(entry.get("folders") or []))
        elif kind == "delete_folders":
            self._undo_deleted_folders_batch(list(entry.get("items") or []))
        else:
            self._set_status_message("Nothing to undo", 3000)

    def _undo_deleted_files_batch(self, batch: list[tuple[Path, Path]]) -> None:
        restored = 0
        failed: list[str] = []
        for original, trash in batch:
            original = Path(original)
            trash = Path(trash)
            try:
                if not trash.is_file():
                    failed.append(f"{original.name}: trash copy missing")
                    continue
                original.parent.mkdir(parents=True, exist_ok=True)
                if original.exists():
                    failed.append(f"{original.name}: already exists")
                    continue
                shutil.move(str(trash), str(original))
                restored += 1
            except OSError as exc:
                failed.append(f"{original.name}: {exc}")
        try:
            parent = Path(batch[0][1]).parent if batch else None
            if parent is not None and parent.is_dir() and not any(parent.iterdir()):
                parent.rmdir()
        except OSError:
            pass

        if self._current_folder:
            self._populate_file_list(self._current_folder)

        if failed and restored:
            QMessageBox.warning(
                self,
                "Restore partially failed",
                f"Restored {restored} file(s).\n\nFailed:\n" + "\n".join(failed[:8]),
            )
        elif failed:
            QMessageBox.warning(
                self,
                "Restore failed",
                "Could not restore deleted file(s):\n" + "\n".join(failed[:8]),
            )
        elif restored:
            self._set_status_message(f"Restored {restored} deleted file(s)", 6000)

    def _undo_hidden_folders(self, folders: list[Path]) -> None:
        restored = 0
        for folder in folders:
            path = Path(folder)
            if not path.is_dir():
                continue
            self._add_root_folder(path, unhide=True, select=False)
            restored += 1
        if restored:
            self._persist_workspace()
            label = "folder" if restored == 1 else "folders"
            self._set_status_message(f"Restored {restored} {label} to the list", 6000)
        else:
            self._set_status_message("Could not restore folder(s) to the list", 5000)

    def _undo_deleted_folders_batch(self, batch: list[tuple[Path, Path]]) -> None:
        restored = 0
        failed: list[str] = []
        for original, trash in batch:
            original = Path(original)
            trash = Path(trash)
            try:
                if not trash.is_dir():
                    failed.append(f"{original.name}: trash copy missing")
                    continue
                if original.exists():
                    failed.append(f"{original.name}: already exists")
                    continue
                original.parent.mkdir(parents=True, exist_ok=True)
                shutil.move(str(trash), str(original))
                restored += 1
                parent_item = self._find_tree_item_for_path(original.parent)
                if parent_item is None:
                    self._add_root_folder(original, unhide=True, select=False)
                else:
                    self._refresh_subfolders_under(
                        original.parent, keep_folder=original
                    )
            except OSError as exc:
                failed.append(f"{original.name}: {exc}")
        try:
            parent = Path(batch[0][1]).parent if batch else None
            if parent is not None and parent.is_dir() and not any(parent.iterdir()):
                parent.rmdir()
        except OSError:
            pass

        if self._current_folder:
            self._populate_file_list(self._current_folder)

        if failed and restored:
            QMessageBox.warning(
                self,
                "Restore partially failed",
                f"Restored {restored} folder(s).\n\nFailed:\n" + "\n".join(failed[:8]),
            )
        elif failed:
            QMessageBox.warning(
                self,
                "Restore failed",
                "Could not restore deleted folder(s):\n" + "\n".join(failed[:8]),
            )
        elif restored:
            label = "folder" if restored == 1 else "folders"
            self._set_status_message(f"Restored {restored} deleted {label}", 6000)

    def _delete_selected_files(self, paths: list[Path]) -> None:
        if not paths:
            return
        paths = self._expand_paths_for_delete(paths)
        if not paths:
            return
        names = "\n".join(p.name for p in paths[:10])
        if len(paths) > 10:
            names += f"\n… and {len(paths) - 10} more"
        answer = QMessageBox.question(
            self,
            "Delete files",
            f"Delete {len(paths)} file(s)?\n\n{names}\n\n"
            "Press ⌘Z / Ctrl+Z afterward to restore.",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return

        deleted, failed = self._soft_delete_files(paths)

        if self._preview_tabs is not None:
            self._preview_tabs.close_paths(paths)

        if self._current_folder:
            self._populate_file_list(self._current_folder)

        if failed:
            QMessageBox.warning(
                self,
                "Delete partially failed",
                f"Deleted {deleted} file(s).\n\nFailed:\n" + "\n".join(failed[:8]),
            )
        else:
            self._set_status_message(
                f"Deleted {deleted} file(s) (⌘Z / Ctrl+Z to restore)",
                6000,
            )

    def _open_files_in_preview(self, paths: list[Path]) -> None:
        """Show the selected file (first only) in the active Preview tab.

        Other open tabs keep their content; use "Open in Additional Tab" to add
        a new one instead of reusing the active tab.
        """
        if self._preview_tabs is None or not paths:
            return
        path = paths[0].expanduser().resolve()
        self._open_preview_tab(path, replace=True)
        self._set_status_file_link(path)
        self._set_status_message(f"Opened in Preview: {path.name}", 5000)

    def _open_files_in_preview_tabs(self, paths: list[Path]) -> None:
        """Add one or more files as additional Preview tabs."""
        if self._preview_tabs is None or not paths:
            return
        self._preview_tabs.open_files_in_tabs(paths)
        self._ensure_preview_visible()
        if len(paths) == 1:
            self._set_status_message(
                f"Opened in Additional Tab: {paths[0].name}",
                5000,
            )
        else:
            self._set_status_message(f"Opened {len(paths)} files in Preview", 5000)

    def _folder_destination_candidates(self, preferred: Path | None) -> list[Path]:
        """Destinations offered by the create-folder dialog, best guess first."""
        candidates: list[Path | None] = [preferred, self._current_folder]
        candidates.extend(self._workspace_root_folders())
        out: list[Path] = []
        for raw in candidates:
            if not isinstance(raw, Path):
                continue
            folder = raw.expanduser()
            try:
                folder = folder.resolve()
            except OSError:
                continue
            if folder.is_dir():
                out.append(folder)
        return out

    def _create_folder(
        self,
        parent: Path | None = None,
        *,
        files: list[Path] | None = None,
    ) -> None:
        """Create a folder, optionally moving the selected files into it."""
        if not isinstance(parent, Path):
            parent = None
        move_files = [p for p in (files or []) if isinstance(p, Path) and p.is_file()]
        destinations = self._folder_destination_candidates(parent)
        if not destinations:
            QMessageBox.warning(
                self,
                "No folder available",
                "Add a folder to the Folders panel first (File → Add Folder).",
            )
            return

        default_name = common_prefix_name(move_files) if move_files else ""
        dialog = CreateFolderDialog(
            self,
            destinations=destinations,
            default_name=default_name or DEFAULT_FOLDER_NAME,
            move_count=len(move_files),
        )
        if not exec_modal_dialog(dialog):
            return
        new_dir = dialog.target_path()
        if new_dir is None:
            return

        try:
            new_dir.mkdir(parents=True, exist_ok=False)
        except OSError as exc:
            QMessageBox.critical(self, "Create folder failed", str(exc))
            return

        if move_files:
            self._transfer_files(move_files, new_dir, mode="cut")

        self._reveal_new_folder(new_dir)
        if move_files:
            noun = "file" if len(move_files) == 1 else "files"
            self._set_status_message(
                f"Created {new_dir.name} with {len(move_files)} {noun}",
                7000,
            )
        else:
            self._set_status_message(f"Created folder {new_dir.name}", 5000)

    def _create_folder_with_files(self, files: list[Path]) -> None:
        """Files-panel entry point: create a folder and move the selection into it."""
        self._create_folder(self._current_folder, files=list(files))

    def _reveal_new_folder(self, new_dir: Path) -> None:
        """Select the folder that was just created, not its parent."""
        # With a nested name (run3/monomers) the immediate parent is new too, so
        # refresh from the nearest ancestor the tree actually knows about.
        anchor: Path | None = None
        for ancestor in new_dir.parents:
            if self._find_tree_item_for_path(ancestor) is not None:
                anchor = ancestor
                break
        if anchor is None:
            # Created outside every workspace root — adopt it as one so it is reachable.
            if self._find_tree_item_for_path(new_dir) is None:
                self._add_root_folder(new_dir, unhide=True)
            return
        # Targeted refresh; a full _refresh_tree() rebuilds every root.
        self._refresh_subfolders_under(anchor, keep_folder=new_dir)
        self._navigate_to_folder(new_dir)

    def _copy_files_to_folder(self, target_dir: Path, files: list[Path]) -> None:
        self._transfer_files(files, target_dir, mode="copy")

    def _clipboard_can_paste(self) -> bool:
        if self._file_clipboard is None:
            return False
        return any(path.is_file() for path in self._file_clipboard.paths)

    def _set_file_clipboard(self, paths: list[Path], *, cut: bool) -> None:
        files = [p.expanduser().resolve() for p in paths if p.expanduser().resolve().is_file()]
        if not files:
            return
        self._file_clipboard = _FileClipboard(
            mode="cut" if cut else "copy",
            paths=files,
        )
        verb = "Cut" if cut else "Copied"
        self._set_status_message(f"{verb} {len(files)} file(s)", 4000)

    def _copy_files_to_clipboard(self, paths: list[Path]) -> None:
        self._set_file_clipboard(paths, cut=False)

    def _cut_files_to_clipboard(self, paths: list[Path]) -> None:
        self._set_file_clipboard(paths, cut=True)

    @staticmethod
    def _unique_path_in_folder(src: Path, dest_dir: Path) -> Path:
        dest_dir = dest_dir.expanduser().resolve()
        candidate = dest_dir / src.name
        if not candidate.exists() or candidate.resolve() == src.resolve():
            return candidate
        return _unique_duplicate_path(candidate)

    def _resolve_transfer_destination(
        self,
        src: Path,
        target_dir: Path,
        *,
        mode: Literal["copy", "cut"],
    ) -> Path | None:
        """Pick destination for a drag/paste; prompt if the name already exists."""
        from repeatviewer.widgets.file_conflict import confirm_existing_file_path

        target_dir = target_dir.expanduser().resolve()
        dest = (target_dir / src.name).expanduser().resolve()
        if src.resolve() == dest.resolve():
            return None
        if not dest.exists():
            return dest
        action = "move" if mode == "cut" else "copy"
        return confirm_existing_file_path(
            self,
            dest,
            title="File already exists",
            informative=(
                f"Overwrite the existing file, or {action} under a different name?"
            ),
        )

    def _transfer_files(
        self,
        files: list[Path],
        target_dir: Path,
        *,
        mode: Literal["copy", "cut"],
    ) -> None:
        target_dir = target_dir.expanduser().resolve()
        if not target_dir.is_dir():
            return

        created: list[Path] = []
        sources: list[Path] = []
        failed: list[str] = []

        for raw in files:
            src = raw.expanduser().resolve()
            if not src.is_file():
                continue
            dest = self._resolve_transfer_destination(src, target_dir, mode=mode)
            if dest is None:
                continue
            try:
                if mode == "cut":
                    if dest.exists():
                        dest.unlink()
                    shutil.move(src, dest)
                    from repeatviewer.widgets.file_colors import replace_file_color_path

                    replace_file_color_path(src, dest)
                else:
                    shutil.copy2(src, dest)
                created.append(dest)
                sources.append(src)
            except OSError as exc:
                failed.append(f"{src.name}: {exc}")

        if mode == "cut" and self._file_clipboard is not None:
            moved_sources = {str(p) for p in sources}
            remaining = [
                p for p in self._file_clipboard.paths if str(p.resolve()) not in moved_sources
            ]
            if remaining:
                self._file_clipboard = _FileClipboard(mode="cut", paths=remaining)
            else:
                self._file_clipboard = None

        if self._current_folder is not None:
            current = self._current_folder.resolve()
            self._populate_file_list(
                self._current_folder,
                new_outputs=[p for p in created if p.parent.resolve() == current],
                source_files=[p for p in sources if p.parent.resolve() == current],
                open_preview=False,
            )

        if created:
            self._sync_folder_tree_labels()

        verb = "Moved" if mode == "cut" else "Copied"
        if created:
            self._set_status_message(f"{verb} {len(created)} file(s) to {target_dir.name}", 7000)
        if failed:
            QMessageBox.warning(
                self,
                f"{'Move' if mode == 'cut' else 'Copy'} partially failed",
                "\n".join(failed[:10]),
            )

    def _paste_files_to_folder(self, target_dir: Path) -> None:
        if not self._clipboard_can_paste() or self._file_clipboard is None:
            return
        self._transfer_files(
            self._file_clipboard.paths,
            target_dir,
            mode=self._file_clipboard.mode,
        )

    def _refresh_tree(self) -> None:
        roots = []
        for i in range(self._folder_tree.topLevelItemCount()):
            item = self._folder_tree.topLevelItem(i)
            raw = item.data(0, ROLE_PATH)
            if raw:
                roots.append(Path(str(raw)))
        keep_folder = (
            self._current_folder.expanduser().resolve()
            if self._current_folder is not None
            else None
        )
        self._folder_tree.clear()
        for root in roots:
            if root.is_dir():
                self._add_root_folder(root, select=False)
        self._sync_folder_tree_placeholders()
        if keep_folder is not None:
            keep_item = self._find_tree_item_for_path(keep_folder)
            if keep_item is not None:
                self._folder_tree.select_folder_item(keep_item)
        self._update_directory_watchers()

    def _apply_function(self, func_name: str, *, preselected: list[Path] | None = None, **job_kwargs) -> None:
        selected = preselected or self._selected_files()
        if not selected:
            QMessageBox.information(self, "No sequence selected", "Select one or more sequence files in the Files panel.")
            return
        seq_files = [p for p in selected if p.suffix.lower() in _FASTA_SUFFIXES]
        if not seq_files:
            QMessageBox.warning(self, "Invalid selection", "Selected files do not include FASTA sequence files.")
            return
        blast_opts_for_gate = job_kwargs.get("blast_options")
        if func_name == "blast":
            from repeatviewer.database_settings import reload_blast_databases_into_project

            # Keep Run / right-click in sync with Manage BLAST add/remove.
            reload_blast_databases_into_project(self._project)
        uses_saved_blast_database = (
            func_name == "blast"
            and blast_opts_for_gate is not None
            and blast_opts_for_gate.active_database_name
        )
        if (
            func_name == "blast"
            and blast_opts_for_gate is not None
            and blast_opts_for_gate.program == BLAST_PROGRAM_BLASTP
            and not blast_opts_for_gate.blastp_database_name
        ):
            QMessageBox.warning(
                self,
                "Protein database required",
                "Select a protein database in the BLASTp tab under BLAST parameters.",
            )
            return
        if (
            func_name == "blast"
            and blast_opts_for_gate is not None
            and blast_opts_for_gate.program == BLAST_PROGRAM_BLASTX
            and not blast_opts_for_gate.blastx_database_name
        ):
            QMessageBox.warning(
                self,
                "Protein database required",
                "blastx searches a protein database. Select one in the BLASTx tab "
                "under BLAST parameters.",
            )
            return
        if (
            func_name in {"blast", "extend", "teaid", "flexidot"}
            and not uses_saved_blast_database
            and not self._project.has_genome_configured()
        ):
            QMessageBox.warning(
                self,
                "Genome required",
                f"Define genome under {self._genome_menu_hint()} before running this tool.",
            )
            return
        if func_name == "extend" and self._project.is_pan_genome_mode:
            names = self._pan_genome_names()
            extend_opts = job_kwargs.get("extend") or self._extend_options
            if names and not extend_opts.pan_genome_name:
                extend_opts = replace(extend_opts, pan_genome_name=names[0])
                job_kwargs["extend"] = extend_opts
        if func_name == "flexidot":
            if flexidot_executable() is None:
                QMessageBox.warning(
                    self,
                    "Flexidot not found",
                    "Flexidot is not on PATH.\n\nInstall with:\n  pip install flexidot\n\n"
                    "If it is installed, add your Python scripts directory to PATH "
                    "(e.g. ~/Library/Python/3.x/bin on macOS).",
                )
                return
        if func_name == "tsd_search":
            tsd_opts = job_kwargs.get("tsd_search_options") or self._tsd_search_options
            if tsd_opts.left_boundary <= 0 or tsd_opts.right_boundary <= 0:
                dialog = TsdSearchParamsDialog(
                    initial=tsd_opts,
                    title="TSD search parameters",
                    parent=self,
                )
                if not dialog.exec():
                    return
                tsd_opts = dialog.options()
                job_kwargs["tsd_search_options"] = tsd_opts
            if tsd_opts.left_boundary >= tsd_opts.right_boundary:
                QMessageBox.warning(
                    self,
                    "TSD search",
                    "Left boundary must be less than the right boundary.",
                )
                return
        show_log = False
        if func_name == "blast":
            blast_opts = job_kwargs.get("blast_options") or self._blast_options
            show_log = bool(getattr(blast_opts, "show_process_window", False))
        elif func_name == "flexidot":
            flex_opts = job_kwargs.get("flexidot_options") or self._flexidot_options
            show_log = bool(getattr(flex_opts, "show_process_window", False))
        worker = FunctionWorker(
            self._run_function_job,
            func_name,
            seq_files,
            job_kwargs,
            stream_log=show_log,
        )
        self._attach_process_monitor(
            worker,
            title=func_name,
            enabled=show_log,
            tool_key=func_name,
        )
        self._start_background_job(
            worker,
            job_title=func_name,
            total_steps=self._job_total_steps(func_name, seq_files, job_kwargs),
            source_files=seq_files,
            on_finished=self._on_worker_finished,
            on_failed=self._on_worker_failed,
        )

    def _command_plugin_work_dir(self, plugin: UserCommandPlugin) -> Path:
        runtime = getattr(plugin, "_job_options", {}).get("runtime_values") or {}
        outdir = str(runtime.get("RV_OUTDIR", "")).strip()
        if outdir:
            path = Path(outdir).expanduser()
            if path.is_dir():
                return path.resolve()
            if path.parent.exists():
                return path.parent.resolve()
        for key in ("RV_INPUT", "RV_GENOME"):
            raw = str(runtime.get(key, "")).strip()
            if not raw:
                continue
            path = Path(raw.split()[0]).expanduser()
            if path.is_file():
                return path.parent.resolve()
            if path.is_dir():
                return path.resolve()
            if path.parent.exists():
                return path.parent.resolve()
        for name, defn in plugin._defs.items():
            if not is_path_var_mode(defn.mode):
                continue
            raw = str(runtime.get(name, "")).strip()
            if not raw:
                continue
            path = Path(raw).expanduser()
            if path.is_file():
                return path.parent.resolve()
            if path.is_dir():
                return path.resolve()
            if path.parent.exists():
                return path.parent.resolve()
        if self._current_folder and self._current_folder.is_dir():
            return self._current_folder.expanduser().resolve()
        if self._project.input_dir and self._project.input_dir.is_dir():
            return self._project.input_dir.expanduser().resolve()
        return Path.cwd()

    def _plugin_job_source_files(self, plugin: UserCommandPlugin, work_dir: Path) -> list[Path]:
        runtime = getattr(plugin, "_job_options", {}).get("runtime_values") or {}
        sources: list[Path] = []
        for key in ("RV_INPUT", "RV_GENOME"):
            raw = str(runtime.get(key, "")).strip()
            if not raw:
                continue
            path = Path(raw.split()[0]).expanduser()
            if path.is_file():
                sources.append(path.resolve())
        for raw in str(runtime.get("RV_INPUTS", "")).split():
            path = Path(raw).expanduser()
            if path.is_file():
                sources.append(path.resolve())
        for name, defn in plugin._defs.items():
            if name in plugin.spec.output_vars:
                continue
            if not is_file_var_mode(defn.mode):
                continue
            raw = str(runtime.get(name, "")).strip()
            if not raw:
                continue
            path = Path(raw).expanduser()
            if path.is_file():
                sources.append(path.resolve())
        resolved_work_dir = work_dir.expanduser().resolve()
        if resolved_work_dir.is_dir():
            sources.append(resolved_work_dir)
        # Deduplicate while preserving order
        seen: set[str] = set()
        unique: list[Path] = []
        for path in sources:
            key = str(path)
            if key in seen:
                continue
            seen.add(key)
            unique.append(path)
        return unique

    def _apply_plugin(self, plugin_id: str, *, preselected: list[Path] | None = None) -> None:
        plugin = self._plugin_registry.get(plugin_id)
        if plugin is None:
            QMessageBox.warning(self, "Plugin", f"Plugin not found: {plugin_id}")
            return
        if isinstance(plugin, UserCommandPlugin):
            if not plugin.prepare_job(self, [], self._project):
                return
            work_dir = self._command_plugin_work_dir(plugin)
            files = [work_dir]
            error = plugin.validate(files, self._project)
            if error:
                QMessageBox.warning(self, plugin.label, error)
                return
            runtime = getattr(plugin, "_job_options", {}).get("runtime_values") or {}
            existing_outputs = [
                path
                for path in existing_output_files(
                    plugin.spec,
                    plugin._defs,
                    runtime_values=runtime,
                    work_dir=work_dir,
                )
                if path.is_file()
            ]
            if existing_outputs:
                listing = "\n".join(f"  • {path}" for path in existing_outputs)
                answer = QMessageBox.question(
                    self,
                    plugin.label,
                    "The following output path(s) already exist:\n\n"
                    f"{listing}\n\n"
                    "Running this plugin may overwrite them. Continue?",
                    QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                    QMessageBox.StandardButton.No,
                )
                if answer != QMessageBox.StandardButton.Yes:
                    return
            func_name = f"plugin:{plugin_id}"
            self._start_plugin_background_job(
                plugin,
                func_name=func_name,
                files=files,
                source_files=self._plugin_job_source_files(plugin, work_dir),
            )
            return
        selected = preselected or self._selected_files()
        initial: list[Path] = []
        if selected:
            initial = self._plugin_registry.filter_files(
                plugin,
                [path.expanduser().resolve() for path in selected],
            )
        if not initial and not plugin.picks_files_in_dialog and plugin.needs_selected_files:
            if plugin.fasta_only:
                QMessageBox.information(
                    self,
                    plugin.label,
                    "Select one or more FASTA files in the Files panel.",
                )
            else:
                QMessageBox.information(
                    self,
                    plugin.label,
                    "Select files in the Files panel.",
                )
            return
        if not plugin.prepare_job(self, initial, self._project):
            return
        job_options = getattr(plugin, "_job_options", {})
        raw_files = job_options.get("query_files") or [str(path) for path in initial]
        files = self._plugin_registry.filter_files(
            plugin,
            [Path(str(path)).expanduser().resolve() for path in raw_files],
        )
        if not files:
            if plugin.fasta_only:
                QMessageBox.warning(
                    self,
                    plugin.label,
                    "No valid query FASTA files were selected.",
                )
            else:
                QMessageBox.warning(self, plugin.label, "No files match this plugin.")
            return
        error = plugin.validate(files, self._project)
        if error:
            QMessageBox.warning(self, plugin.label, error)
            return
        func_name = f"plugin:{plugin_id}"
        self._start_plugin_background_job(
            plugin,
            func_name=func_name,
            files=files,
            source_files=files,
        )

    def _start_plugin_background_job(
        self,
        plugin,
        *,
        func_name: str,
        files: list[Path],
        source_files: list[Path],
    ) -> None:
        job_options = getattr(plugin, "_job_options", {}) or {}
        show_log = bool(job_options.get("show_process_window", False))
        worker = FunctionWorker(
            self._run_plugin_job,
            func_name,
            files,
            {},
            stream_log=show_log,
        )
        self._attach_process_monitor(
            worker,
            title=plugin.label,
            enabled=show_log,
            tool_key=f"plugin_{getattr(plugin, 'id', 'tool')}",
            mirror_progress=True,
        )
        self._start_background_job(
            worker,
            job_title=func_name,
            total_steps=self._job_total_steps(func_name, files, {}),
            source_files=source_files,
            on_finished=self._on_worker_finished,
            on_failed=self._on_worker_failed,
        )

    def _populate_plugins_menu_actions(self) -> None:
        """Deprecated path — use ``_rebuild_others_menu``."""
        self._rebuild_others_menu()

    def _ensure_others_menu(self) -> QMenu | None:
        """Return a live Others menu, recreating it if the cached one was deleted."""
        if self._menu_object_alive(self._others_menu):
            return self._others_menu

        self._others_menu = None
        self._plugins_menu = None
        for action in list(self.menuBar().actions()):
            candidate = action.menu()
            if not self._menu_object_alive(candidate):
                label = self._clean_menu_label(action.text()).lower()
                if label in {"others", "plugins"}:
                    self.menuBar().removeAction(action)
                continue
            title = self._clean_menu_label(candidate.title()).lower()
            if title in {"others", "plugins"}:
                # Migrate old Plugins title to Others.
                candidate.setTitle("&Others")
                self._others_menu = candidate
                self._plugins_menu = candidate
                return candidate

        self._remove_menubar_menus_named("Others")
        self._remove_menubar_menus_named("Plugins")
        menu = QMenu("&Others", self)
        self._others_menu = menu
        self._plugins_menu = menu
        insert_before = None
        for action in self.menuBar().actions():
            label = self._clean_menu_label(action.text())
            if label in {"View", "Help"}:
                insert_before = action
                break
        if insert_before is not None:
            self.menuBar().insertMenu(insert_before, menu)
        else:
            self.menuBar().addMenu(menu)
        return menu

    def _ensure_plugins_menu(self) -> QMenu | None:
        return self._ensure_others_menu()

    def _rebuild_others_menu(self) -> None:
        def _fill(menu: QMenu) -> None:
            menu.clear()
            plugin_menu = menu.addMenu("Plugin (Beta)")
            wrapper_menu = menu.addMenu("Command wrapper")

            plugin_menu.addAction("Manage plugin", self._open_plugin_manager)
            plugin_menu.addSeparator()
            installed = [
                plugin
                for plugin in self._plugin_registry.all_plugins()
                if not isinstance(plugin, UserCommandPlugin)
            ]
            if installed:
                for plugin in installed:
                    plugin_menu.addAction(
                        plugin.label,
                        lambda checked=False, pid=plugin.id: self._apply_plugin(pid),
                    )
            else:
                none_action = plugin_menu.addAction("(No installed plugins)")
                none_action.setEnabled(False)

            wrapper_menu.addAction("Manage wrapper", self._open_wrapper_manager)
            wrapper_menu.addSeparator()
            specs = load_user_command_plugin_specs()
            if specs:
                for spec in specs:
                    wrapper_menu.addAction(
                        spec.label,
                        lambda checked=False, pid=spec.id: self._apply_plugin(pid),
                    )
            else:
                empty = wrapper_menu.addAction("(No wrappers yet)")
                empty.setEnabled(False)

            self._others_menu = menu
            self._plugins_menu = menu

        menu = self._ensure_others_menu()
        if menu is None:
            return
        try:
            _fill(menu)
        except RuntimeError:
            self._others_menu = None
            self._plugins_menu = None
            self._remove_menubar_menus_named("Others")
            self._remove_menubar_menus_named("Plugins")
            menu = self._ensure_others_menu()
            if menu is None:
                return
            try:
                _fill(menu)
            except RuntimeError:
                return
        self._connect_menu_operation_logging()

    def _rebuild_plugins_menu(self) -> None:
        self._rebuild_others_menu()

    def _menu_object_alive(self, obj) -> bool:
        if obj is None:
            return False
        try:
            from shiboken6 import isValid

            return bool(isValid(obj))
        except Exception:  # noqa: BLE001
            try:
                _ = obj.objectName()
                return True
            except RuntimeError:
                return False

    def _remove_menubar_menus_named(self, title: str) -> None:
        """Drop menubar entries whose label matches ``title`` (including dead menus)."""
        wanted = title.strip().lower()
        for action in list(self.menuBar().actions()):
            label = self._clean_menu_label(action.text()).lower()
            menu = action.menu()
            menu_label = ""
            if self._menu_object_alive(menu):
                menu_label = self._clean_menu_label(menu.title()).lower()
            if label == wanted or menu_label == wanted:
                self.menuBar().removeAction(action)

    def _ensure_workflows_menu(self) -> QMenu | None:
        """Return a live Workflows menu, recreating it if the cached one was deleted."""
        if self._menu_object_alive(self._workflows_menu):
            return self._workflows_menu

        self._workflows_menu = None
        for action in list(self.menuBar().actions()):
            candidate = action.menu()
            if not self._menu_object_alive(candidate):
                if self._clean_menu_label(action.text()).lower() == "workflows":
                    self.menuBar().removeAction(action)
                continue
            if self._clean_menu_label(candidate.title()).lower() == "workflows":
                self._workflows_menu = candidate
                return candidate

        # Own the QMenu from the main window so macOS native menubar sync
        # cannot leave us holding a deleted temporary menu.
        menu = QMenu("&Workflows", self)
        self._workflows_menu = menu
        insert_before = None
        for action in self.menuBar().actions():
            if self._clean_menu_label(action.text()) in {
                "Others",
                "Plugins",
                "View",
                "Help",
            }:
                insert_before = action
                break
        if insert_before is not None:
            self.menuBar().insertMenu(insert_before, menu)
        else:
            self.menuBar().addMenu(menu)
        return menu

    def _rebuild_workflows_menu(self) -> None:
        from repeatviewer.workflows import ensure_default_workflows

        specs = ensure_default_workflows()

        def _fill(menu: QMenu) -> None:
            menu.clear()
            menu.addAction(
                "Design and manage your workflows",
                self._open_workflow_designer,
            )
            menu.addSeparator()
            if specs:
                for spec in specs:
                    menu.addAction(
                        spec.name,
                        lambda checked=False, wid=spec.id: self._run_workflow(wid),
                    )
            else:
                none_action = menu.addAction("(No workflows yet)")
                none_action.setEnabled(False)

        menu = self._ensure_workflows_menu()
        if menu is None:
            return
        try:
            _fill(menu)
        except RuntimeError:
            # Menu died mid-rebuild (common when triggered from the menu itself
            # on macOS). Tear down and recreate after that stack unwinds.
            self._workflows_menu = None
            self._remove_menubar_menus_named("Workflows")
            menu = self._ensure_workflows_menu()
            if menu is None:
                return
            try:
                _fill(menu)
            except RuntimeError:
                return
        self._connect_menu_operation_logging()

    def _schedule_workflows_menu_rebuild(self) -> None:
        """Defer rebuild until after the triggering menu action fully unwinds."""
        QTimer.singleShot(0, self._on_workflows_catalog_changed)

    def _on_workflows_catalog_changed(self) -> None:
        """Refresh Workflows menu and toolbar after save/rename/delete/reorder."""
        from repeatviewer.views.toolbar_dialog import workflow_tool_id
        from repeatviewer.workflows import load_workflow_specs

        valid_workflow_tools = {
            workflow_tool_id(spec.id) for spec in load_workflow_specs()
        }
        pruned = [
            tool_id
            for tool_id in self._toolbar_order
            if not str(tool_id).startswith("workflow:")
            or tool_id in valid_workflow_tools
        ]
        if pruned != list(self._toolbar_order):
            self._toolbar_order = pruned
            try:
                from repeatviewer.views.toolbar_dialog import save_toolbar_order

                save_toolbar_order(self._toolbar_order)
            except Exception:  # noqa: BLE001
                pass
        self._rebuild_workflows_menu()
        self._rebuild_top_toolbar()

    def _open_workflow_designer(self) -> None:
        from repeatviewer.views.workflow_dialog import WorkflowDesignerDialog

        dialog = WorkflowDesignerDialog(parent=self)
        # Live updates while the designer is open (deferred — safe on macOS).
        dialog.workflows_changed.connect(self._schedule_workflows_menu_rebuild)
        dialog.exec()
        # Final sync after close in case the last action did not emit.
        self._schedule_workflows_menu_rebuild()

    def _run_workflow(self, workflow_id: str) -> None:
        from repeatviewer.job_options import ALIGNER_MAFFT
        from repeatviewer.workflows import get_workflow_spec, run_workflow, validate_workflow_steps

        spec = get_workflow_spec(workflow_id)
        if spec is None:
            QMessageBox.warning(self, "Workflow", "That workflow was not found.")
            self._schedule_workflows_menu_rebuild()
            return
        error = validate_workflow_steps(spec.tool_ids())
        if error:
            QMessageBox.warning(self, "Workflow", error)
            return
        selected = self._selected_files()
        seq_files = [p for p in selected if p.suffix.lower() in _FASTA_SUFFIXES]
        if not seq_files:
            QMessageBox.information(
                self,
                "No sequence selected",
                "Select one or more FASTA files in the Files panel before running a workflow.",
            )
            return
        if spec.needs_genome() and not self._project.has_genome_configured():
            from repeatviewer.workflows import resolve_tool_info

            needing = []
            for s in spec.steps:
                info = resolve_tool_info(s.tool_id)
                if info is not None and info.needs_genome:
                    needing.append(s.tool_id)
            only_blast = set(needing) == {"blast"}
            if not (only_blast and self._blast_options.active_database_name):
                QMessageBox.warning(
                    self,
                    "Genome required",
                    f"Define genome under {self._genome_menu_hint()} before running this workflow.",
                )
                return
        if any(s.tool_id == "flexidot" for s in spec.steps) and flexidot_executable() is None:
            QMessageBox.warning(
                self,
                "Flexidot not found",
                "Flexidot is not on PATH.\n\nInstall with:\n  pip install flexidot",
            )
            return

        job_kwargs = {
            "spec": spec,
            "seq_files": seq_files,
            "project": self._project,
            "cancel": self._job_cancel,
            "max_workers": self._project.batch_workers,
            "blast_options": self._blast_options,
            "extend": self._extend_options,
            "teaid_options": self._teaid_options,
            "msa_clean_options": self._msa_clean_options,
            "consensus_options": self._consensus_options,
            "flexidot_options": self._flexidot_options,
            "orf_prediction_settings": self._orf_prediction_settings,
            "align_options": self._align_options,
            "use_system_mafft": self._project.use_system_mafft,
            "use_system_muscle": self._project.use_system_muscle,
            "align_engine": (
                self._align_options.right_click_aligner
                if self._align_options.right_click_aligner in ALIGNERS
                else self._align_options.combined_aligner
            )
            or ALIGNER_MAFFT,
        }

        def _run(progress_cb, **kwargs):
            return run_workflow(progress_cb=progress_cb, **kwargs)

        worker = WorkflowWorker(_run, job_kwargs)
        self._attach_process_monitor(
            worker,
            title=f"Workflow: {spec.name}",
            enabled=False,
            tool_key=f"workflow_{spec.id}",
        )
        self._start_background_job(
            worker,
            job_title=f"workflow:{spec.id}",
            total_steps=(
                len(spec.steps) * len(seq_files)
                if len(seq_files) > 1
                else len(spec.steps)
            ),
            source_files=seq_files,
            on_finished=self._on_worker_finished,
            on_failed=self._on_worker_failed,
        )

    def _open_wrapper_manager(self) -> None:
        dialog = WrapperManageDialog(
            on_save=self._save_designed_plugin,
            on_delete=self._delete_custom_plugin,
            parent=self,
        )
        exec_modal_dialog(dialog)

    def _open_plugin_manager(self) -> None:
        dialog = InstalledPluginManageDialog(
            parent=self,
            on_changed=self._reload_entry_point_plugins,
        )
        exec_modal_dialog(dialog)

    def _reload_entry_point_plugins(self) -> None:
        from repeatviewer.plugins.loader import reload_plugins
        from repeatviewer.views.toolbar_dialog import load_toolbar_order

        self._plugin_registry = reload_plugins()
        self._toolbar_order = load_toolbar_order(self._plugin_registry)
        QTimer.singleShot(0, self._rebuild_others_menu)
        QTimer.singleShot(0, self._rebuild_top_toolbar)

    def _open_command_plugin_designer(self) -> None:
        if self._command_plugin_dialog is not None and self._command_plugin_dialog.isVisible():
            self._command_plugin_dialog.raise_()
            self._command_plugin_dialog.activateWindow()
            return
        dialog = CommandPluginDesignerDialog(parent=self)
        self._command_plugin_dialog = dialog
        dialog.setModal(False)
        dialog.setWindowModality(Qt.WindowModality.NonModal)
        dialog.accepted.connect(lambda d=dialog: self._save_designed_plugin(d.spec()))
        dialog.finished.connect(self._on_command_plugin_dialog_finished)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _on_command_plugin_dialog_finished(self, _code: int) -> None:
        if self._command_plugin_dialog is not None:
            self._command_plugin_dialog.deleteLater()
            self._command_plugin_dialog = None

    def _ensure_wrapper_on_toolbar(self, plugin_id: str, *, enabled: bool) -> None:
        from repeatviewer.views.toolbar_dialog import plugin_tool_id, save_toolbar_order

        tool_id = plugin_tool_id(plugin_id)
        order = list(self._toolbar_order)
        if enabled:
            if tool_id not in order:
                order.append(tool_id)
                self._toolbar_order = order
                save_toolbar_order(self._toolbar_order, self._plugin_registry)
        else:
            if tool_id in order:
                order = [item for item in order if item != tool_id]
                self._toolbar_order = order
                save_toolbar_order(self._toolbar_order, self._plugin_registry)

    def _save_designed_plugin(self, spec, *, previous_id: str | None = None) -> None:
        plugin = UserCommandPlugin(spec)
        try:
            if previous_id and previous_id != spec.id:
                delete_user_command_plugin_spec(previous_id)
                self._plugin_registry.unregister(previous_id)
            was_existing = self._plugin_registry.get(spec.id) is not None
            save_user_command_plugin_spec(spec)
            self._plugin_registry.register_or_replace(plugin)
            # Defer menubar/toolbar refresh until after modal dialogs fully close
            # (macOS can otherwise hit a deleted QMenu during clear/rebuild).
            QTimer.singleShot(0, self._rebuild_others_menu)
            QTimer.singleShot(0, self._rebuild_top_toolbar)
            if was_existing or (previous_id and previous_id != spec.id):
                self._set_status_message(f"Updated wrapper {spec.label}", 5000)
            else:
                self._set_status_message(f"Added wrapper {spec.label}", 5000)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.warning(self, "Design command wrapper", str(exc))

    def _delete_custom_plugin(self, plugin_id: str) -> bool:
        if not delete_user_command_plugin_spec(plugin_id):
            return False
        self._plugin_registry.unregister(plugin_id)
        # Drop from toolbar if the user had added it via Setting → Toolbar.
        self._ensure_wrapper_on_toolbar(plugin_id, enabled=False)
        QTimer.singleShot(0, self._rebuild_others_menu)
        QTimer.singleShot(0, self._rebuild_top_toolbar)
        self._set_status_message(f"Deleted wrapper {plugin_id}", 5000)
        return True

    def _run_plugin_job(
        self,
        func_name: str,
        seq_files: list[Path],
        progress_cb,
        **job_kwargs: Any,
    ) -> dict:
        plugin_id = func_name[len("plugin:") :]
        return self._plugin_registry.run_plugin(
            plugin_id,
            seq_files,
            self._project,
            progress_cb,
            max_workers=self._project.batch_workers,
            cancel=self._job_cancel,
        )

    def _run_function_job(
        self,
        func_name: str,
        seq_files: list[Path],
        progress_cb,
        **job_kwargs: Any,
    ):
        return run_batch(
            func_name,
            seq_files,
            self._project,
            progress_cb,
            max_workers=self._project.batch_workers,
            cancel=self._job_cancel,
            **job_kwargs,
        )

    def _start_background_job(
        self,
        worker: FunctionWorker | AlignWorker | WorkflowWorker,
        *,
        job_title: str,
        total_steps: int,
        source_files: list[Path] | None = None,
        on_finished,
        on_failed,
    ) -> None:
        if self._worker_thread is not None:
            QMessageBox.information(
                self,
                "Job already running",
                "Wait for the current job to finish before starting another. "
                "You can still browse folders and preview files while it runs.",
            )
            return

        self._ensure_genome_paths_without_spaces()
        prepared_sources = list(source_files or [])
        if prepared_sources:
            prepared_sources = self._rename_paths_with_spaces(
                prepared_sources,
                kind="input file",
            )
            # Keep the worker in sync with renamed paths.
            if hasattr(worker, "_seq_files"):
                worker._seq_files = list(prepared_sources)
            if hasattr(worker, "_job_kwargs") and isinstance(
                getattr(worker, "_job_kwargs", None), dict
            ):
                if "seq_files" in worker._job_kwargs:
                    worker._job_kwargs["seq_files"] = list(prepared_sources)

        self._last_job_source_files = list(prepared_sources)
        self._clear_job_file_highlights()
        # Cancel enrichment without blocking the GUI (was wait(1000)).
        self._stop_metadata_enrichment(wait_ms=0)
        self._job_cancel.reset()
        self._job_finished_callback = on_finished
        self._job_failed_callback = on_failed
        self._show_job_progress(job_title, total_steps)
        self._set_job_actions_enabled(False)

        self._worker_thread = QThread(self)
        self._worker = worker
        self._worker.moveToThread(self._worker_thread)
        self._worker_thread.started.connect(self._worker.run)
        queued = Qt.ConnectionType.QueuedConnection
        self._worker.progress.connect(self._on_job_progress, queued)
        # Bound MainWindow slots only — bare lambdas/closures have no receiver
        # QObject and can run on the worker thread (QBasicTimer / segfault).
        self._worker.finished.connect(self._on_background_job_finished, queued)
        self._worker.failed.connect(self._on_background_job_failed, queued)
        self._worker.finished.connect(self._worker_thread.quit)
        self._worker.failed.connect(self._worker_thread.quit)
        self._worker_thread.finished.connect(self._cleanup_worker)
        self._worker_thread.finished.connect(self._update_status)
        self._worker_thread.start()

    @Slot(int, str)
    def _on_job_progress(self, value: int, label: str) -> None:
        operation_logger.info(
            "JOB STEP name=%s step=%s detail=%s",
            self._job_display_title(getattr(self, "_active_job_log_title", "job")),
            value,
            label,
        )
        if self._task_progress is not None:
            # Never fill the last slot here — reserved until results are shown.
            work = getattr(self, "_job_work_steps", self._task_progress._total)
            capped = min(max(0, int(value)), max(0, int(work)))
            # One-file jobs: keep the status as "Running BLAST…" without a name.
            detail = ""
            if int(work) > 1:
                detail = self._parse_progress_filename(label)
            self._task_progress.update(capped, detail)

    def _log_job_finished(self, job_title: str, result) -> None:
        if isinstance(result, (list, tuple)):
            outputs = ", ".join(str(item) for item in result)
        else:
            outputs = str(result) if result is not None else ""
        operation_logger.info(
            "JOB FINISH name=%s outputs=[%s]",
            self._job_display_title(job_title),
            outputs,
        )

    def _log_job_failed(self, job_title: str, error) -> None:
        operation_logger.error(
            "JOB ERROR name=%s error=%s",
            self._job_display_title(job_title),
            error,
        )

    def _set_job_actions_enabled(self, enabled: bool) -> None:
        for control in self._job_toolbar_actions:
            control.setEnabled(enabled)

    def _refresh_folders_after_job(
        self,
        created_paths: list[Path],
        source_paths: list[Path],
        *,
        extra_folders: set[Path] | None = None,
    ) -> None:
        """Refresh file list; show outputs in each input file's folder.

        Never auto-opens preview here — callers use ``_auto_open_outputs`` (or
        open a specific result) so large jobs do not open twice.
        """
        self._watcher_refresh_suppressed_until = time.monotonic() + 3.0
        touched = {p.parent.resolve() for p in created_paths} | {
            p.parent.resolve() for p in source_paths
        }
        if extra_folders:
            touched |= {folder.resolve() for folder in extra_folders}
        for folder in touched:
            ancestor = folder.parent
            while ancestor != ancestor.parent:
                item = self._find_tree_item_for_path(ancestor)
                if item is not None:
                    self._refresh_subfolders_under(ancestor)
                    break
                ancestor = ancestor.parent

        current = self._current_folder.resolve() if self._current_folder else None
        # Huge jobs: skip input tinting (every selected FASTA) — only mark outputs.
        # Keeps the finish path closer to O(new files) instead of O(inputs+outputs).
        highlight_sources = source_paths
        if len(created_paths) > _LARGE_JOB_OUTPUT_BATCH:
            highlight_sources = []
        if current is not None and current in touched:
            # created/source paths are already resolve()'d by the finish handler.
            self._populate_file_list(
                current,
                new_outputs=[p for p in created_paths if p.parent == current],
                source_files=[p for p in highlight_sources if p.parent == current],
                open_preview=False,
            )
            return

        if created_paths:
            target = created_paths[0].parent.resolve()
            self._navigate_to_folder(
                target,
                new_outputs=[p for p in created_paths if p.parent == target],
                source_files=[p for p in highlight_sources if p.parent == target],
                open_preview=False,
            )
            return

        if source_paths and self._current_folder:
            self._populate_file_list(self._current_folder, open_preview=False)

    @Slot(dict)
    def _on_worker_finished(self, result: dict) -> None:
        created_paths: list[Path] = []
        for raw in result.get("created_outputs", []):
            path = Path(str(raw)).expanduser()
            if path.is_file():
                created_paths.append(path.resolve())
            elif path.is_dir():
                created_paths.append(path.resolve())

        source_paths: list[Path] = []
        extra_folders: set[Path] = set()
        for raw in result.get("source_files", []):
            path = Path(str(raw)).expanduser()
            if path.is_file():
                source_paths.append(path.resolve())
            elif path.is_dir():
                extra_folders.add(path.resolve())

        self._refresh_folders_after_job(
            created_paths,
            source_paths,
            extra_folders=extra_folders,
        )
        func_name = result.get("func_name")
        if func_name == "te_annotation":
            fasta_sources = [
                p for p in source_paths if p.suffix.lower() in _FASTA_SUFFIXES
            ]
            if fasta_sources:
                self._auto_open_outputs(fasta_sources)
        elif func_name and str(func_name).startswith("plugin:"):
            preview_files = [p for p in created_paths if p.is_file()]
            self._auto_open_outputs(preview_files)
        else:
            self._auto_open_outputs(created_paths)

        pending_viz = self._pending_te_viz
        pending_fastas = self._pending_te_genome_map_fastas
        pending_map_opts = self._pending_te_genome_map
        if pending_viz and result.get("func_name") == "blast":
            self._pending_te_viz = None
            self._pending_te_genome_map_fastas = []
            self._pending_te_genome_map = None
            resolved = self._resolve_te_viz_hit_paths_and_kinds(
                pending_viz,
                created_paths=created_paths,
            )
            if resolved is None:
                example = (
                    f"{pending_viz.consensus_fastas[0].name}_b.txt"
                    if pending_viz.consensus_fastas
                    else "query.fa_b.txt"
                )
                QMessageBox.warning(
                    self,
                    "TE genome map",
                    "Blast finished but a BLAST table for the TE consensus was not found.\n\n"
                    "Expected a hit table named like:\n"
                    f"  {example}\n"
                    "next to the consensus FASTA (legacy *_blast.txt is also accepted).",
                )
            else:
                te_paths, te_kinds = resolved
                self._open_blast_genome_map(
                    te_paths,
                    gene_files=pending_viz.gene_paths,
                    methylation_files=pending_viz.methylation_paths,
                    te_input_kinds=te_kinds,
                    map_options=pending_viz.map_options,
                )
        elif pending_fastas and result.get("func_name") == "blast":
            self._pending_te_genome_map_fastas = []
            self._pending_te_genome_map = None
            self._open_te_genome_map_after_blast(pending_fastas, pending_map_opts, result)

        failures = result.get("failures", [])
        done_n = max(0, result["total"] - len(failures))
        if result.get("workflow_name"):
            title = f"Workflow {result['workflow_name']}"
            unit = "step(s)"
        else:
            title = MainWindow._job_display_title(func_name or "")
            unit = "sequence(s)"
        summary = f"{title} finished — {done_n}/{result['total']} {unit}"
        if created_paths:
            names = ", ".join(p.name for p in created_paths[:4])
            if len(created_paths) > 4:
                names += ", …"
            summary += f" — created: {names}"
        # TE-Aid opens its HTML in Preview and sets the file link; a timed
        # "finished" banner only flashes then covers that link.
        if func_name != "teaid":
            self._set_status_message(summary, 12000, override_file_link=True)

        if failures:
            msg = summary + "\n\nFailures:\n" + "\n".join(failures[:8])
            if len(failures) > 8:
                msg += "\n..."
            QMessageBox.warning(self, "Function completed with issues", msg)

    @Slot(str)
    def _on_worker_failed(self, error: str) -> None:
        self._pending_te_genome_map_fastas = []
        self._pending_te_genome_map = None
        self._pending_te_viz = None
        if error == CANCELLED_MESSAGE:
            self._set_status_message("Stopped", 6000)
            return
        self._set_status_message("Job failed", 8000)
        QMessageBox.critical(self, "Function failed", error)

    def _cleanup_worker(self) -> None:
        cancelled = self._job_cancel.requested
        if cancelled:
            self._job_cancel.reset()
        if self._worker:
            self._worker.deleteLater()
            self._worker = None
        if self._worker_thread:
            self._worker_thread.deleteLater()
            self._worker_thread = None
        if cancelled:
            self._set_status_message("Stopped", 6000)
        self._update_status()

    def _stop_current_job(self) -> None:
        """Stop currently running background job/map build."""
        if not (
            self._worker_thread is not None
            or self._genome_map_thread is not None
            or self._benchmark_thread is not None
            or self._extract_thread is not None
        ):
            self._set_status_message("No running task to stop", 3000)
            return

        # Kill child tools immediately (non-blocking on the GUI thread), then
        # ask worker QThreads to wind down. Do not wait here — that froze Stop.
        self._job_cancel.request()
        self._pending_te_genome_map_fastas = []
        self._pending_te_genome_map = None
        self._pending_te_viz = None
        self._queued_genome_map_jobs.clear()
        self._set_status_message("Stopping…", 8000)
        if self._stop_btn is not None:
            self._stop_btn.setEnabled(False)

        threads: list[QThread] = []
        if self._worker_thread is not None:
            threads.append(self._worker_thread)
        if self._genome_map_thread is not None:
            threads.append(self._genome_map_thread)
        if self._benchmark_thread is not None:
            threads.append(self._benchmark_thread)
        if self._extract_thread is not None:
            threads.append(self._extract_thread)

        for thread in threads:
            if thread.isRunning():
                thread.requestInterruption()
                thread.quit()

        # Re-enable toolbar promptly; cleanup finishes via failed/finished.
        self._finish_background_job()

    def _open_te_genome_map_after_blast(
        self,
        fasta_paths: list[Path],
        map_options: TeGenomeMapOptions | None,
        result: dict,
    ) -> None:
        failed_names: set[str] = set()
        for entry in result.get("failures", []):
            text = str(entry)
            if ": " in text:
                failed_names.add(text.split(": ", 1)[0].strip())
            else:
                failed_names.add(text.strip())
        queued: list[Path] = []
        for fasta in fasta_paths:
            if fasta.name in failed_names:
                continue
            blast_txt = self._find_blast_table_after_job(
                fasta,
                created_paths=[
                    Path(str(raw)).expanduser().resolve()
                    for raw in result.get("created_outputs", [])
                ],
            )
            if blast_txt is not None:
                queued.append(blast_txt)
        if queued:
            self._enqueue_genome_map_jobs(queued, map_options)
        else:
            QMessageBox.warning(
                self,
                "TE genome map",
                "Blast finished but the BLAST table for the TE genome map was not found.",
            )

    def _pick_preview_output(self, created_paths: list[Path]) -> Path | None:
        """Prefer HTML/PDF/images for preview after TE-Aid / Flexidot / genome map."""
        if not created_paths:
            return None

        def _as_path(raw: Path) -> Path:
            try:
                return raw.expanduser().resolve()
            except OSError:
                return raw.expanduser()

        # Large BLAST-style batches: avoid N× is_file() — trust worker paths.
        # Prefer a visual output if present; otherwise the last previewable file.
        if len(created_paths) > _LARGE_JOB_OUTPUT_BATCH:
            preferred_suffixes = {
                ".html",
                ".htm",
                ".pdf",
                ".png",
                ".jpg",
                ".jpeg",
                ".svg",
            }
            for raw in reversed(created_paths):
                path = _as_path(raw)
                if path.suffix.lower() in preferred_suffixes:
                    return path
            for raw in reversed(created_paths):
                path = _as_path(raw)
                if is_previewable_file(path):
                    return path
            return _as_path(created_paths[-1])

        files = [p for p in (_as_path(p) for p in created_paths) if p.is_file()]
        if not files:
            return None
        for path in files:
            if path.name.endswith("_sankey.html"):
                return path
        for path in files:
            if path.name.endswith("_library_benchmark.html"):
                return path
        for path in files:
            if path.name.endswith("_genome_map_region.html") or path.name.endswith(
                "_blast_genome_map_region.html"
            ):
                return path
        for path in files:
            if path.name.endswith("_genome_map.html") or path.name.endswith(
                "_blast_genome_map.html"
            ):
                return path
        in_folder = [
            p
            for p in files
            if self._current_folder and p.parent.resolve() == self._current_folder.resolve()
        ]
        search = in_folder if in_folder else files
        previewable = [p for p in search if is_previewable_file(p)]
        if not previewable:
            # Still return a created file so the Files list can select it.
            return search[-1] if search else files[-1]
        for suffix in (".html", ".htm", ".pdf", ".png", ".jpg", ".jpeg", ".svg"):
            for path in previewable:
                if path.suffix.lower() == suffix:
                    return path
        # Among equal priorities, prefer the last created path (newest output).
        order = {str(p): index for index, p in enumerate(files)}
        return min(
            previewable,
            key=lambda p: (preview_path_priority(p), -order.get(str(p), 0)),
        )

    def _auto_open_outputs(self, created_paths: list[Path]) -> None:
        """Open the best new output and select it (deep blue); others stay green."""
        target = self._pick_preview_output(created_paths)
        if target is None:
            return
        # Suppress selection→preview while we select, so we open exactly once.
        was_suppress = self._file_list.suppress_selection_preview()
        self._file_list.set_suppress_selection_preview(True)
        try:
            self._file_list.clearSelection()
            self._file_list.setCurrentItem(None)
            # Do not full-restyle the folder — select path re-tints job keys only.
            self._select_path_in_file_list(target)
            if is_previewable_file(target):
                self._open_preview_tab(target, replace=True)
                self._select_path_in_file_list(target)
        finally:
            if not was_suppress:
                # Keep suppress through deferred currentItemChanged timers.
                QTimer.singleShot(0, self._file_list.clear_suppress_selection_preview)

    def _update_status(self) -> None:
        if self._job_running:
            return
        if self._status_linked_file is not None:
            # Keep the clickable filename on the left; don't restore folder text.
            self.statusBar().clearMessage()
            self._status_file_link.show()
            return
        if self.statusBar().currentMessage():
            # A timed message (usually a job result) is on screen. Leave it alone.
            return
        if self._current_folder and self._cached_folder_files is not None:
            visible = self._file_list.rowCount()
            if visible == 1:
                only = self._file_list.item(0)
                if only is not None and only.data(ROLE_PATH) is None:
                    visible = 0
            self._sync_files_status(visible, len(self._cached_folder_files))
            return
        parts = [APP_DISPLAY_NAME]
        if self._project.input_dir:
            parts.append(self._project.input_dir.name)
        if self._project.is_pan_genome_mode:
            parts.append(f"pan ({len(self._project.pan_genomes)})")
        elif self._project.genome_file:
            parts.append(self._project.genome_file.name)
        self._set_status_message(" · ".join(parts))

    def _show_about(self) -> None:
        from repeatviewer.widgets.app_icon import app_icon

        dialog = QDialog(self)
        dialog.setWindowTitle(f"About {APP_DISPLAY_NAME}")
        dialog.setModal(True)

        icon_label = QLabel()
        icon_pix = app_icon().pixmap(QSize(64, 64))
        icon_label.setPixmap(icon_pix)
        icon_label.setFixedSize(64, 64)

        text_label = QLabel(
            f"<b>{APP_DISPLAY_NAME}</b> {__version__}<br/>"
            "<span style='color:#555555; font-weight:normal;'>"
            "Developed at the Chair of Molecular Botany,<br/>"
            "RWTH Aachen University"
            "</span>"
        )
        text_label.setTextFormat(Qt.TextFormat.RichText)
        text_label.setAlignment(
            Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft
        )
        text_label.setWordWrap(True)

        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(16)
        row.addWidget(icon_label, alignment=Qt.AlignmentFlag.AlignVCenter)
        row.addWidget(text_label, stretch=1, alignment=Qt.AlignmentFlag.AlignVCenter)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok)
        buttons.accepted.connect(dialog.accept)

        layout = QVBoxLayout(dialog)
        layout.setContentsMargins(20, 20, 20, 16)
        layout.setSpacing(16)
        layout.addLayout(row)
        layout.addWidget(buttons, alignment=Qt.AlignmentFlag.AlignRight)

        dialog.setFixedWidth(max(320, dialog.sizeHint().width()))
        dialog.exec()
