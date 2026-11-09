"""Dialogs for tool-specific parameters (right-click / custom runs)."""

from __future__ import annotations

import re
import sys
from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import QEvent, QItemSelectionModel, QPoint, QSize, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QMouseEvent, QStandardItemModel
from PySide6.QtWidgets import (
    QAbstractItemView,
    QAbstractScrollArea,
    QCheckBox,
    QColorDialog,
    QCompleter,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QButtonGroup,
    QPlainTextEdit,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QStyle,
    QStyleOptionViewItem,
    QStackedWidget,
    QTabBar,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
    QTextEdit,
)

# Use the wheel-ignoring spin boxes everywhere in these parameter dialogs, so
# scrolling over a field never changes its value by accident -- only the arrows
# or typing do. Aliased to the standard names so all existing construction
# sites pick them up unchanged.
from repeatviewer.widgets.spin_boxes import (  # noqa: E402
    NoWheelDoubleSpinBox as QDoubleSpinBox,
    NoWheelSpinBox as QSpinBox,
)
from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox  # noqa: E402
from repeatviewer.widgets.path_browse import PathBrowseRow  # noqa: E402
from repeatviewer.widgets.help_button import (  # noqa: E402
    HELP_BUTTON_SIZE,
    field_with_help,
    label_with_help,
    make_help_button,
)


def exec_modal_dialog(dialog: QDialog) -> bool:
    """Show a modal dialog above the main window (avoids WebEngine / dock click stealing)."""
    dialog.setModal(True)
    dialog.setWindowModality(Qt.WindowModality.ApplicationModal)
    dialog.adjustSize()
    dialog.raise_()
    dialog.activateWindow()
    return dialog.exec() == QDialog.DialogCode.Accepted


def _prepare_option_dialog(dialog: QDialog) -> None:
    dialog.setModal(True)
    dialog.setWindowModality(Qt.WindowModality.ApplicationModal)
    flags = dialog.windowFlags()
    dialog.setWindowFlags(
        flags
        | Qt.WindowType.Dialog
        | Qt.WindowType.WindowTitleHint
        | Qt.WindowType.WindowCloseButtonHint
    )


def system_path_executable_label(tool_name: str) -> str:
    """Checkbox text for Options → Executables (Align / BLAST / Map / Tree)."""
    return (
        f"Use {tool_name} from the system PATH or defined Conda environment"
    )


def resolve_system_executable(*names: str) -> str | None:
    """First match from the Settings conda env, then PATH."""
    import shutil

    from repeatviewer.conda_tools import which_tool

    for name in names:
        found = which_tool(name) or shutil.which(name)
        if found:
            return str(found)
    return None


def make_executable_status_label(parent: QWidget | None = None) -> QLabel:
    """Status line next to a PATH checkbox or under a named tool (FastGA / HISAT2)."""
    label = QLabel(parent)
    label.setWordWrap(False)
    label.setTextInteractionFlags(
        Qt.TextInteractionFlag.TextSelectableByMouse
    )
    label.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
    # Grey info line (path / not-found / not configured).
    label.setStyleSheet("color: palette(mid);")
    return label


def _executable_status_width(label: QLabel) -> int:
    """Usable width for eliding a long executable path."""
    width = label.width()
    if width >= 40:
        return width
    parent = label.parentWidget()
    while parent is not None:
        pw = parent.width()
        if pw >= 80:
            # Leave room for group-box / dialog chrome.
            return max(40, pw - 48)
        parent = parent.parentWidget()
    # Match TOOL_PARAMS_DIALOG_WIDTH when the widget is not laid out yet.
    return max(40, 608 - 80)


def _elide_middle(text: str, label: QLabel, *, reserve: int = 0) -> str:
    """Middle-elide ``text`` to fit the status label; full string stays in tooltip."""
    if not text:
        return text
    available = max(40, _executable_status_width(label) - max(0, reserve))
    return label.fontMetrics().elidedText(
        text, Qt.TextElideMode.ElideMiddle, available
    )


def set_path_checkbox_status(
    label: QLabel,
    path: str | None,
    *,
    checked: bool,
    has_builtin: bool = True,
) -> None:
    """Update status for a “Use … from PATH” checkbox.

    When checked: grey path (elided if long) or grey not-found message.
    When unchecked: hide the status line (no empty gap).
    """
    label.setStyleSheet("color: palette(mid);")
    if not checked:
        label.clear()
        label.setToolTip("")
        label.hide()
        return
    label.show()
    if path:
        elided = _elide_middle(path, label)
        label.setText(elided)
        label.setToolTip(path)
        return
    if has_builtin:
        text = "not found on PATH — using the built-in one"
    else:
        text = "not found on PATH or Conda environment"
    label.setText(text)
    label.setToolTip(text)


def set_named_executable_status(
    label: QLabel,
    name: str,
    path: str | None,
    *,
    missing_text: str = "not configured",
) -> None:
    """Grey ``**FastGA**: /path`` (elided) or ``**FastGA**: not configured``."""
    label.setStyleSheet("color: palette(mid);")
    label.show()
    prefix = f"<b>{name}</b>: "
    # Approximate bold-name width so the path elides inside the remaining space.
    reserve = label.fontMetrics().horizontalAdvance(f"{name}: ")
    if path:
        elided = _elide_middle(path, label, reserve=reserve)
        label.setText(f"{prefix}{elided}")
        label.setToolTip(path)
    else:
        label.setText(f"{prefix}{missing_text}")
        label.setToolTip(missing_text)


def executable_checkbox_row(
    parent: QWidget,
    checkbox: QCheckBox,
    status: QLabel,
) -> QWidget:
    """Checkbox; grey path / not-found on the next line, flush with the tick."""
    wrap = QWidget(parent)
    lay = QVBoxLayout(wrap)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(2)
    # Keep indicator tight to the label (macOS styles can add a wide gap).
    checkbox.setStyleSheet("QCheckBox { spacing: 6px; }")
    lay.addWidget(checkbox)
    # Same left edge as the checkbox — no extra indent before the status line.
    status.setContentsMargins(0, 0, 0, 0)
    status.hide()
    lay.addWidget(status)
    return wrap


def make_process_output_box(
    parent: QWidget,
    *,
    tool_name: str,
    checked: bool = True,
    close_when_finished: bool = True,
) -> tuple[QGroupBox, QCheckBox, QCheckBox]:
    """Options-tab group for the live process monitor window.

    Stays compact: never stretch the two checkboxes apart when the Options
    page is taller than its sizeHint (e.g. after visiting a Parameters tab).
    """
    checkbox = QCheckBox(f"Show process window while {tool_name} runs", parent)
    checkbox.setChecked(bool(checked))
    checkbox.setToolTip(
        f"Open a window that first shows the {tool_name} command line, then "
        "streams progress and tool messages. Useful for long runs."
    )
    close_box = QCheckBox("Close process window when finished", parent)
    close_box.setChecked(bool(close_when_finished))
    close_box.setEnabled(checkbox.isChecked())
    close_box.setToolTip(
        "Automatically close the process window after a successful run. "
        "On failure the window stays open so you can read the log."
    )
    checkbox.toggled.connect(close_box.setEnabled)
    box = QGroupBox("Process output", parent)
    box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
    layout = QVBoxLayout(box)
    layout.setContentsMargins(12, 10, 12, 10)
    layout.setSpacing(6)
    layout.setAlignment(Qt.AlignmentFlag.AlignTop)
    layout.addWidget(checkbox)
    layout.addWidget(close_box)
    return box, checkbox, close_box


def make_citation_box(parent: QWidget, html: str) -> QGroupBox:
    """Options-tab group listing tool paper citations (clickable DOIs)."""
    box = QGroupBox("Citation", parent)
    box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
    layout = QVBoxLayout(box)
    layout.setContentsMargins(12, 10, 12, 10)
    layout.setSpacing(4)
    layout.setAlignment(Qt.AlignmentFlag.AlignTop)
    label = QLabel(parent)
    label.setWordWrap(True)
    label.setTextFormat(Qt.TextFormat.RichText)
    label.setOpenExternalLinks(True)
    label.setTextInteractionFlags(Qt.TextInteractionFlag.TextBrowserInteraction)
    # RichText ignores QFont for HTML; set size in markup (smaller than body).
    label.setText(f"<div style='font-size:11px'>{html}</div>")
    label.setObjectName("ToolCitationBody")
    layout.addWidget(label)
    return box


def make_cancel_save_run_bar(
    dialog: QDialog,
    *,
    left_widgets: list[QWidget] | None = None,
    on_accept=None,
    include_run: bool = True,
) -> QWidget:
    """Footer with optional left controls plus Cancel / Save [/ Run].

    Sets ``dialog._close_action`` to ``\"save\"`` or ``\"run\"`` before closing.
    ``dialog.should_run()`` should return ``dialog._close_action == \"run\"``.
    Pass ``include_run=False`` for settings-only dialogs (Cancel / Save).
    """
    dialog._close_action = "cancel"
    bar = QWidget(dialog)
    row = QHBoxLayout(bar)
    row.setContentsMargins(0, 0, 0, 0)
    row.setSpacing(8)
    for widget in left_widgets or []:
        row.addWidget(widget)
    row.addStretch(1)
    cancel_btn = QPushButton("Cancel", bar)
    save_btn = QPushButton("Save", bar)
    save_btn.setDefault(True)
    dialog._save_button = save_btn
    dialog._ok_button = save_btn

    def _finish(action: str) -> None:
        check = getattr(dialog, "_run_check", None)
        if action == "run" and callable(check):
            problem = check(dialog)
            if problem:
                # Say what is missing and stay open, with everything still filled in.
                QMessageBox.information(dialog, dialog.windowTitle() or "Run", problem)
                return
        dialog._close_action = action
        if on_accept is not None:
            on_accept()
        else:
            dialog.accept()

    cancel_btn.clicked.connect(dialog.reject)
    save_btn.clicked.connect(lambda: _finish("save"))
    row.addWidget(cancel_btn)
    row.addWidget(save_btn)
    if include_run:
        run_btn = QPushButton("Run", bar)
        run_btn.setDefault(True)
        save_btn.setDefault(False)
        dialog._run_button = run_btn
        dialog._ok_button = run_btn
        run_btn.clicked.connect(lambda: _finish("run"))
        row.addWidget(run_btn)
    else:
        dialog._run_button = None
    return bar


def dialog_should_run(dialog: QDialog) -> bool:
    """True when the dialog was closed with Run."""
    return getattr(dialog, "_close_action", "") == "run"


# Outer width shared with AlignParamsDialog's default sizeHint (FastGA panel
# minimum 560 plus dialog margins / chrome). BLAST / Flexidot / TE-Aid match this.
TOOL_PARAMS_DIALOG_WIDTH = 608
# Wider dialog when Inputs need Browse/Remove path rows (e.g. Genome density plot).
# Use this when the user asks for a "longer" / "wider" window width.
TOOL_PARAMS_DIALOG_WIDTH_WIDE = 720
# Numeric/combo field width on Parameters forms (Align / BLAST / TE-Aid).
TOOL_PARAMS_FIELD_WIDTH = 260
# Longest Settings ``Foo:`` label across Align / BLAST / Mapping / TE-Aid / Flexidot /
# MSA / FastGA. Spacers and Settings labels share this width so every tool’s field
# column starts on the same vertical line. Keep this the longest real label — do
# not let per-dialog ``ref_label`` widen one tab past the others.
TOOL_PARAMS_LABEL_COLUMN_REF = "Minimum alignment length:"


def tool_params_dialog_width(*, define_inputs: bool = False) -> int:
    """Dialog width: 608 normally, 720 when Define-inputs path rows are shown."""
    if define_inputs:
        return TOOL_PARAMS_DIALOG_WIDTH_WIDE
    return TOOL_PARAMS_DIALOG_WIDTH


def pin_tool_tabs_height(tabs: QTabWidget, page_h: int | None) -> None:
    """Pin ``QTabWidget`` height to the active page so unused space cannot reappear.

    Pass ``None`` to clear pins **before** measuring. After ``resize``, call again
    with the measured active-page height. Required for Sequence Annotation,
    Flexidot, TE-Aid, HMMER, Genome density, MSA Clean, Consensus, etc.
    """
    if page_h is None:
        tabs.setMinimumHeight(0)
        tabs.setMaximumHeight(16777215)
        for i in range(tabs.count()):
            page = tabs.widget(i)
            if page is not None:
                page.setMinimumHeight(0)
                page.setMaximumHeight(16777215)
        return
    tab_bar_h = max(24, tabs.tabBar().sizeHint().height())
    tabs_h = tab_bar_h + max(1, int(page_h)) + 8
    tabs.setMinimumHeight(0)
    tabs.setMaximumHeight(tabs_h)
    current = tabs.currentWidget()
    if current is not None:
        current.setMinimumHeight(0)
        current.setMaximumHeight(max(1, int(page_h)))


def settings_form_label_column_width(
    host: QWidget,
    forms: list[QFormLayout | None] | None = None,
    *,
    ref_label: str | None = None,
) -> int:
    """Shared Settings label-column width (same across all tool dialogs)."""
    del forms, ref_label  # width is global; args kept for call-site compatibility
    return max(1, host.fontMetrics().horizontalAdvance(TOOL_PARAMS_LABEL_COLUMN_REF))


def sync_form_field_column_spacers(
    host: QWidget,
    spacer_labels: list[QLabel],
    forms: list[QFormLayout | None],
    *,
    ref_label: str | None = None,
) -> None:
    """Indent empty-label rows to the shared Settings field column.

    Spacer ``QLabel``s on Inputs / Mode / checkbox rows, and every Settings
    ``*:`` label in ``forms``, are fixed to ``TOOL_PARAMS_LABEL_COLUMN_REF``.
    Colon labels are right-aligned so short text (``Map mode:``, ``Format:``)
    sits next to the field instead of leaving a gap on the right of the label.
    """
    width = settings_form_label_column_width(host, forms, ref_label=ref_label)
    for lab in spacer_labels:
        lab.setFixedWidth(width)
    for form in forms:
        if form is None:
            continue
        for row in range(form.rowCount()):
            item = form.itemAt(row, QFormLayout.ItemRole.LabelRole)
            widget = item.widget() if item is not None else None
            if isinstance(widget, QLabel) and widget.text().endswith(":"):
                widget.setFixedWidth(width)
                widget.setAlignment(
                    Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
                )


def _run_cancel_button_box(parent: QWidget) -> QDialogButtonBox:
    buttons = QDialogButtonBox(
        QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
        parent=parent,
    )
    ok = buttons.button(QDialogButtonBox.StandardButton.Ok)
    if ok is not None:
        ok.setText("Run")
    return buttons


from repeatviewer.branding import APP_DISPLAY_NAME
from repeatviewer.backends.blast_genome_plot import (
    chromosome_names_from_genome,
    format_genome_length_bp,
    genome_chrom_sizes_from_fasta,
)
from repeatviewer.backends.database_prepare import blast_database_is_ready
from repeatviewer.project import BLAST_MOLECULE_NUCL, BLAST_MOLECULE_PROT, RepeatViewerProject
from repeatviewer.ui.theme import item_view_checkbox_stylesheet
from repeatviewer.job_options import (
    ALIGNER_MAFFT,
    ALIGNER_MUSCLE,
    TE_GENOME_MAP_COLOR_FAMILY,
    TE_GENOME_MAP_COLOR_IDENTITY,
    TE_GENOME_MAP_COLOR_ORDER,
    TE_GENOME_MAP_COLOR_SUPERFAMILY,
    TE_GENOME_MAP_DENSITY_HIT_COUNT,
    TE_GENOME_MAP_DENSITY_TE_LENGTH,
    TE_GENOME_MAP_DIALOG_TAB_ANNOTATION,
    TE_GENOME_MAP_DIALOG_TAB_METHYLATION,
    TE_GENOME_MAP_DIALOG_TAB_TE,
    TE_GENOME_MAP_DIALOG_TABS,
    DEFAULT_ANNOTATION_DENSITY,
    DEFAULT_ANNOTATION_FEATURE_TYPES,
    DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED,
    DEFAULT_SHOW_PROCESS_WINDOW,
    DEFAULT_TE_GENOME_MAP_BIN_BP,
    DEFAULT_TE_GENOME_MAP_COLLAPSE_OVERLAPS,
    DEFAULT_TE_GENOME_MAP_COLOR_BY,
    DEFAULT_TE_GENOME_MAP_DENSITY,
    DEFAULT_TE_GENOME_MAP_METHYLATION_COVERAGE_MIN,
    DEFAULT_TE_GENOME_MAP_SHOW_DIFF_TRACK,
    DEFAULT_TE_GENOME_MAP_SHOW_METH_DIFF_TRACK,
    DEFAULT_TE_GENOME_MAP_SPLIT_STRAND,
    CHROMOSOME_SORT_NATURAL,
    CHROMOSOME_SORT_SIZE_ASC,
    CHROMOSOME_SORT_SIZE_DESC,
    TeGenomeMapOptions,
    TeDivergencePlotOptions,
    MsaPlotJobOptions,
    RepeatQuantifyOptions,
    REPEAT_QUANTIFY_RIGHT_CLICK_DIALOG,
    REPEAT_QUANTIFY_RIGHT_CLICK_RERUN,
    TeConfusionMatrixBenchmarkOptions,
    TeLibraryBenchmarkOptions,
    DEFAULT_TE_BENCHMARK_THREADS,
    DEFAULT_TE_DIVERGENCE_BIN_WIDTH,
    DEFAULT_BLAST_EVALUE,
    DEFAULT_BLAST_MAX_MSA_LINES,
    DEFAULT_BLAST_MINIMUM_BIT_SCORE,
    DEFAULT_BLAST_THREADS,
    DEFAULT_BLAST_TOP_MSA_LINES,
    BLAST_PROGRAM_BLASTN,
    BLAST_PROGRAM_BLASTP,
    BLAST_PROGRAM_BLASTX,
    BLAST_PROGRAM_TBLASTN,
    BLAST_PROGRAM_TBLASTX,
    BLAST_PROGRAM_ONLINE,
    BLAST_PROGRAMS,
    BLAST_DIALOG_TABS,
    BLAST_ONLINE_ALGORITHMS,
    DIALOG_TAB_OPTIONS,
    opening_dialog_tab,
    persistable_dialog_tab,
    DEFAULT_BLAST_ONLINE_DATABASE,
    DEFAULT_BLAST_ONLINE_EVALUE,
    DEFAULT_BLAST_ONLINE_HITLIST,
    DEFAULT_BLAST_ONLINE_PROGRAM,
    DEFAULT_BLAST_ONLINE_WRITE_TEXT,
    DEFAULT_CONSENSUS_AMBIGUOUS,
    DEFAULT_CONSENSUS_THRESHOLD,
    CONSENSUS_RIGHT_CLICK_RUN,
    CONSENSUS_RIGHT_CLICK_DIALOG,
    DEFAULT_EXTEND_LEFT,
    DEFAULT_EXTEND_RIGHT,
    EXTEND_LENGTH_PRESETS,
    EXTEND_RIGHT_CLICK_RUN,
    EXTEND_RIGHT_CLICK_DIALOG,
    FASTGA_RIGHT_CLICK_RUN,
    FASTGA_MODE_SELF,
    FASTGA_MODE_PAIRS,
    FASTGA_MODE_REFERENCE,
    FASTGA_INPUT_SELECTED,
    FASTGA_INPUT_DEFINE,
    ALIGN_INPUT_SELECTED,
    ALIGN_INPUT_DEFINE,
    FASTGA_FORMAT_PAF,
    FASTGA_FORMAT_PAFX,
    FASTGA_FORMAT_PAFM,
    FASTGA_FORMAT_PAFS,
    FASTGA_FORMAT_PAFS_LONG,
    FASTGA_FORMAT_PSL,
    FASTGA_FORMAT_ALN,
    ALIGNER_FASTGA,
    DEFAULT_FLEXIDOT_LCS_SHADING,
    DEFAULT_FLEXIDOT_OUTPUT_FORMAT,
    DEFAULT_FLEXIDOT_WORDSIZE,
    DEFAULT_FLEXIDOT_MODES,
    DOTPLOT_DIALOG_TAB_BUILTIN,
    DOTPLOT_DIALOG_TAB_FLEXIDOT,
    FLEXIDOT_SEQUENCE_AA,
    FLEXIDOT_SEQUENCE_NUC,
    FLEXIDOT_RIGHT_CLICK_BUILTIN,
    FLEXIDOT_RIGHT_CLICK_FLEXIDOT,
    FLEXIDOT_RIGHT_CLICK_DIALOG,
    normalize_dotplot_dialog_tab,
    normalize_flexidot_right_click_action,
    TE_ANNOTATION_RIGHT_CLICK_RUN,
    TE_ANNOTATION_RIGHT_CLICK_DIALOG,
    TeAnnotationJobOptions,
    FLEXIDOT_OUTPUT_HTML,
    FLEXIDOT_OUTPUT_PNG,
    FLEXIDOT_OUTPUT_PDF,
    FLEXIDOT_OUTPUT_SVG,
    DEFAULT_DIV_THRESHOLD,
    DEFAULT_DIV_WINDOW_SIZE,
    DEFAULT_GAP_COL_GAP_THRESHOLD,
    DEFAULT_BLAST_EVALUE,
    DEFAULT_HMMER_EVALUE,
    DEFAULT_HMMER_DOMAIN_EVALUE,
    DEFAULT_HMMER_MAX_HITS,
    DEFAULT_HMMER_MINIMUM_BIT_SCORE,
    DEFAULT_HMMER_MAX_ORFS,
    DEFAULT_HMMER_WRITE_TEXT,
    HMMER_DIALOG_TAB_PARAMETERS,
    HMMER_DIALOG_TABS,
    HMMER_RIGHT_CLICK_RUN,
    HMMER_RIGHT_CLICK_DIALOG,
    DEFAULT_GAP_COL_MIN_NUCLEOTIDE,
    DEFAULT_GAP_COL_SIMI_GAP_THR,
    DEFAULT_GAP_COL_SIMILARITY_THR,
    DEFAULT_GAP_ROW_GAP_THRESHOLD,
    DEFAULT_GAP_ROW_WINDOW_SIZE,
    DEFAULT_MSA_CLEAN_ORDER,
    MSA_CLEAN_RIGHT_CLICK_COMBINED,
    MSA_CLEAN_RIGHT_CLICK_GAPPY,
    MSA_CLEAN_RIGHT_CLICK_DIV,
    MSA_CLEAN_RIGHT_CLICK_GAP,
    MSA_CLEAN_RIGHT_CLICK_DIALOG,
    msa_clean_method_description,
    DEFAULT_TEAID_EVALUE,
    DEFAULT_TEAID_X_AXIS_LABEL,
    DEFAULT_TEAID_SELF_BLAST_WORD_SIZE,
    DEFAULT_TEAID_HIT_LEN_HIGH_PCT,
    DEFAULT_TEAID_HIT_LEN_MID_PCT,
    DEFAULT_TEAID_HIT_COLOR_HIGH,
    DEFAULT_TEAID_HIT_COLOR_MID,
    DEFAULT_TEAID_HIT_LINE_TRANSPARENCY,
    MAX_TEAID_SELF_BLAST_WORD_SIZE,
    MIN_TEAID_SELF_BLAST_WORD_SIZE,
    TEAID_RIGHT_CLICK_RUN,
    TEAID_RIGHT_CLICK_DIALOG,
    clamp_teaid_hit_len_pct,
    clamp_teaid_hit_line_transparency,
    BlastOptions,
    HmmerOptions,
    ConsensusOptions,
    normalize_consensus_ambiguous,
    normalize_consensus_sequence_type,
    configure_report_evalue_spinbox,
    DivergenceRowOptions,
    AlignOptions,
    ALIGN_RUN_MODE_ADD_TO_MSA,
    ALIGN_RUN_MODE_COMBINED,
    ALIGN_RUN_MODE_SEPARATE,
    ALIGN_RUN_MODES,
    MAFFT_ADD_METHOD_ADD,
    MAFFT_ADD_METHOD_ADDFRAGMENTS,
    MAFFT_ADD_METHOD_ADDFULL,
    MAFFT_ADD_METHOD_ADDLONG,
    MAFFT_SEQUENCE_AUTO,
    MAFFT_SEQUENCE_NUCLEOTIDE,
    MAFFT_SEQUENCE_PROTEIN,
    MAFFT_STRATEGY_AUTO,
    MAFFT_STRATEGY_CUSTOM,
    MAFFT_STRATEGY_EINSI,
    MAFFT_STRATEGY_FFT_NS_1,
    MAFFT_STRATEGY_FFT_NS_2,
    MAFFT_STRATEGY_GINSI,
    MAFFT_STRATEGY_LINSI,
    MUSCLE_MODE_ALIGN,
    MUSCLE_MODE_SUPER5,
    MUSCLE_PERM_ABC,
    MUSCLE_PERM_ACB,
    MUSCLE_PERM_BCA,
    MUSCLE_PERM_NONE,
    ExtendOptions,
    FlexidotOptions,
    MAX_FLEXIDOT_WORDSIZE,
    MIN_FLEXIDOT_WORDSIZE,
    GapRowOptions,
    GappyColumnOptions,
    MsaCleanOptions,
    TeaAidOptions,
    TSD_BOUNDARY_REPEATVIEWER,
    TSD_BOUNDARY_SOURCES,
    TSD_BOUNDARY_USER,
    TsdSearchOptions,
    normalize_blast_evalue,
    normalize_flexidot_output_format,
    normalize_flexidot_sequence_type,
)
from repeatviewer.plugins.custom_command import (
    CommandVarDef,
    UserCommandPluginSpec,
    display_var_name,
    get_user_command_plugin_spec,
    is_affixed_variable,
    is_special_variable,
    load_user_command_plugin_specs,
    normalize_var_mode,
    normalize_var_name,
    plugin_id_from_label,
    special_variable_mode,
    special_variables_help_text,
    template_affixed_variables,
    template_variables,
    user_template_variables,
)


class _CommandVarRow(QWidget):
    """One variable row; fixed-width type / role / Remove keep columns aligned."""

    remove_requested = Signal(object)
    role_changed = Signal(object, str)

    _NAME_MIN_WIDTH = 140
    _DEFAULT_MIN_WIDTH = 220
    _MODE_WIDTH = 104
    _ROLE_WIDTH = 88
    _REMOVE_WIDTH = 80

    def __init__(
        self,
        initial: CommandVarDef | None = None,
        *,
        role: str = "",
        special: bool = False,
        parent=None,
    ) -> None:
        super().__init__(parent)
        var = initial or CommandVarDef(name="", mode="file", value="")
        self._special = bool(special) or (
            is_special_variable(var.name) and not is_affixed_variable(var.name)
        )
        name_key = normalize_var_name(var.name)
        self._affixed = is_affixed_variable(name_key)
        self._name = QLineEdit(display_var_name(name_key) if name_key else "", self)
        self._name.setPlaceholderText("${variable_name} or ${RV_INPUT}_b.text")
        self._name.setMinimumWidth(self._NAME_MIN_WIDTH)
        self._name.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._mode = QComboBox(self)
        self._mode.addItem("File", "file")
        self._mode.addItem("Folder", "folder")
        self._mode.addItem("Text", "text")
        self._mode.addItem("Integer", "integer")
        self._mode.addItem("Float", "float")
        self._mode.setFixedWidth(self._MODE_WIDTH)
        self._mode.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        if self._special and name_key:
            mode_key = special_variable_mode(name_key)
        elif self._affixed:
            mode_key = normalize_var_mode(var.mode) or special_variable_mode(name_key)
        else:
            mode_key = normalize_var_mode(var.mode)
        index = self._mode.findData(mode_key)
        if index >= 0:
            self._mode.setCurrentIndex(index)
        self._default = QLineEdit(var.value, self)
        self._default.setPlaceholderText("Default variable")
        self._default.setMinimumWidth(self._DEFAULT_MIN_WIDTH)
        self._default.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._role = QComboBox(self)
        self._role.addItem("—", "")
        self._role.addItem("Input", "input")
        self._role.addItem("Output", "output")
        self._role.setFixedWidth(self._ROLE_WIDTH)
        self._role.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._role.setToolTip(
            "Input: the FASTA (or path) passed into this step when used in a workflow.\n"
            "Output: a result file/folder produced by the command "
            "(shown in the file panel).\n"
            "Only one Input; multiple Outputs are allowed."
        )
        role_key = str(role or "").strip().lower()
        if role_key not in {"input", "output"}:
            role_key = ""
        role_index = self._role.findData(role_key)
        if role_index >= 0:
            self._role.setCurrentIndex(role_index)
        self._role.currentIndexChanged.connect(self._on_role_changed)
        remove_btn = QPushButton("Remove", self)
        remove_btn.clicked.connect(lambda: self.remove_requested.emit(self))
        self._remove_btn = remove_btn
        remove_btn.setFixedWidth(self._REMOVE_WIDTH)
        remove_btn.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(8)
        row.addWidget(self._name, 2)
        row.addWidget(self._mode, 0)
        row.addWidget(self._default, 3)
        row.addWidget(self._role, 0)
        row.addWidget(remove_btn, 0)
        if self._special:
            self._name.setReadOnly(True)
            self._mode.setEnabled(False)
            self._default.setEnabled(False)
            self._default.setPlaceholderText("Filled at run from selection / genome")
            self._name.setToolTip("Special variable filled at run.")
        elif self._affixed:
            self._name.setReadOnly(True)
            self._default.setPlaceholderText(
                "Leave empty to auto-fill as base path + suffix at run"
            )
            self._name.setToolTip(
                "Separate parameter from the glued suffix in the command "
                "(e.g. ${RV_INPUT}_b.text). Editable at run; default = base + suffix."
            )
        self._refresh_role_fields()

    def _on_role_changed(self, _index: int = 0) -> None:
        role = str(self._role.currentData() or "")
        self._refresh_role_fields()
        if role:
            self.role_changed.emit(self, role)

    def set_role(self, role: str) -> None:
        key = str(role or "").strip().lower()
        if key not in {"input", "output"}:
            key = ""
        index = self._role.findData(key)
        if index < 0:
            return
        self._role.blockSignals(True)
        self._role.setCurrentIndex(index)
        self._role.blockSignals(False)
        self._refresh_role_fields()

    def _refresh_role_fields(self) -> None:
        # Specials keep their built-in type; do not rewrite mode when I/O is set.
        if self._special:
            return
        role = str(self._role.currentData() or "")
        if role in {"input", "output"}:
            mode = normalize_var_mode(str(self._mode.currentData()))
            if mode not in {"file", "folder"}:
                index = self._mode.findData("file")
                if index >= 0:
                    self._mode.setCurrentIndex(index)

    def definition(self) -> tuple[CommandVarDef | None, str]:
        name = normalize_var_name(self._name.text())
        if not name:
            return None, ""
        shown = display_var_name(name)
        if self._name.text().strip() != shown:
            self._name.blockSignals(True)
            self._name.setText(shown)
            self._name.blockSignals(False)
        role = str(self._role.currentData() or "")
        if is_special_variable(name) and not is_affixed_variable(name):
            return (
                CommandVarDef(
                    name=name,
                    mode=special_variable_mode(name),
                    value="",
                ),
                role,
            )
        mode = str(self._mode.currentData())
        if is_affixed_variable(name) and not mode:
            mode = special_variable_mode(name)
        var = CommandVarDef(
            name=name,
            mode=mode,
            value=self._default.text().strip(),
        )
        return var, role

    def mark_removable(self, removable: bool) -> None:
        if self._special:
            self._remove_btn.setEnabled(False)
            self._remove_btn.setVisible(True)
            return
        self._remove_btn.setEnabled(removable)
        self._remove_btn.setVisible(True)

    @property
    def is_special(self) -> bool:
        return self._special


class CommandPluginDesignerDialog(QDialog):
    """Design a user-defined command plugin from a command template."""

    # Preferred open size. Must be reflected in sizeHint() because
    # exec_modal_dialog / dialog polish call adjustSize(), which otherwise
    # shrinks back to the content sizeHint and ignores resize().
    _DEFAULT_WIDTH = 1100
    _DEFAULT_HEIGHT = 640

    def __init__(
        self,
        *,
        initial: UserCommandPluginSpec | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Design command wrapper")
        self.setMinimumWidth(720)
        self.setMinimumHeight(520)
        self.resize(self._DEFAULT_WIDTH, self._DEFAULT_HEIGHT)
        self.setSizeGripEnabled(True)
        spec = initial or UserCommandPluginSpec(
            id="",
            label="",
            command_template="",
            output_vars=[],
            var_defs=[],
            fasta_only=True,
        )
        self._initial_fasta_only = spec.fasta_only
        self._initial_id = spec.id
        self._rebuilding_rows = False
        self._label = QLineEdit(spec.label, self)
        self._label.setPlaceholderText("Wrapper name (must be unique)")
        self._tooltip = QTextEdit(self)
        self._tooltip.setPlaceholderText(
            "Optional. Shown as the toolbar / menu tooltip, e.g. "
            "Run blastn of the selected FASTA against the defined genome."
        )
        self._tooltip.setPlainText(spec.tooltip)
        self._tooltip.setFixedHeight(52)
        self._tooltip.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._tooltip.setLineWrapMode(QTextEdit.LineWrapMode.WidgetWidth)
        self._command = QTextEdit(self)
        self._command.setFixedHeight(140)
        self._command.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._command.setPlaceholderText(
            "blastn -query ${RV_INPUT} -db ${RV_GENOME} "
            "-out ${RV_OUTDIR}/${RV_STEM}_b.txt "
            "-evalue 1e-5 -outfmt 6 -num_threads 4\n"
            "\n"
            "Glue a fixed suffix after a variable: ${RV_STEM}_b.txt\n"
            "Use ${RV_STEM} (name without extension), not ${RV_INPUT}_b.txt."
        )
        self._command.setPlainText(spec.command_template)

        special_hint = QLabel(
            "Special (auto-filled): ${RV_INPUT}, ${RV_INPUTS}, ${RV_GENOME}, "
            "${RV_FOLDER}, ${RV_OUTDIR}, ${RV_STEM}, ${RV_NAME}",
            self,
        )
        special_hint.setWordWrap(True)
        special_hint.setStyleSheet("color: #666666;")
        help_text = special_variables_help_text()
        help_btn = make_help_button(self, "Special variables", help_text)
        special_row = QWidget(self)
        special_layout = QHBoxLayout(special_row)
        special_layout.setContentsMargins(0, 0, 0, 0)
        special_layout.setSpacing(6)
        special_layout.addWidget(special_hint, 1)
        special_layout.addWidget(help_btn, 0)

        self._rows_host = QWidget(self)
        self._rows_layout = QVBoxLayout(self._rows_host)
        self._rows_layout.setContentsMargins(0, 0, 0, 0)
        self._rows_layout.setSpacing(6)
        self._rows_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        self._rows: list[_CommandVarRow] = []
        self._vars_scroll = QScrollArea(self)
        self._vars_scroll.setWidgetResizable(True)
        self._vars_scroll.setFrameShape(QFrame.Shape.StyledPanel)
        self._vars_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self._vars_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )
        self._vars_scroll.setWidget(self._rows_host)
        self._vars_scroll.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )
        # Scroll inside this area so upper fields stay put; dialog itself can resize.
        self._vars_scroll.setMinimumHeight(132)

        out_set = {normalize_var_name(v) for v in spec.output_vars}
        in_set = {normalize_var_name(v) for v in (spec.input_vars or [])}
        saved: dict[str, CommandVarDef] = {
            normalize_var_name(v.name): v
            for v in spec.var_defs
            if normalize_var_name(v.name) and not is_special_variable(v.name)
        }
        # Show specials used in the command first (read-only name/type; I/O editable).
        for name in template_variables(spec.command_template):
            if is_special_variable(name):
                role = "output" if name in out_set else ("input" if name in in_set else "")
                self._add_var_row(
                    CommandVarDef(name=name, mode=special_variable_mode(name), value=""),
                    role=role,
                    special=True,
                )
        # Affixed tokens like ${RV_INPUT}_b.text are separate editable parameters.
        for name in template_affixed_variables(spec.command_template):
            existing = saved.pop(name, None)
            role = "output" if name in out_set else ("input" if name in in_set else "")
            self._add_var_row(
                existing
                or CommandVarDef(name=name, mode=special_variable_mode(name), value=""),
                role=role,
            )
        for name, var in saved.items():
            role = "output" if name in out_set else ("input" if name in in_set else "")
            self._add_var_row(var, role=role)
        if not self._rows:
            self._add_var_row()

        detect_btn = QPushButton("Detect variables from command", self)
        detect_btn.clicked.connect(self._detect_variables)
        add_var_btn = QPushButton("Add variable", self)
        add_var_btn.clicked.connect(self._on_add_variable_clicked)

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        form.addRow("Wrapper name:", self._label)
        form.addRow("Description:", self._tooltip)
        form.addRow("Command template:", self._command)
        form.addRow("", special_row)
        btn_row = QWidget(self)
        btn_layout = QHBoxLayout(btn_row)
        btn_layout.setContentsMargins(0, 0, 0, 0)
        btn_layout.addWidget(detect_btn)
        btn_layout.addStretch(1)
        form.addRow("", btn_row)
        define_help = (
            "Input / Output mark variables for workflow chaining "
            "(sequence in → result files out).\n\n"
            "• Input: the file or folder this step receives from the previous "
            "step (or the selected FASTA). Only one Input can be defined.\n"
            "• Output: a file or folder this command produces. Multiple "
            "Outputs can be defined; they are also shown in the file panel.\n\n"
            "Input and Output variables must use type File or Folder."
        )
        define_label = QLabel("Define variables:", self)
        define_help_btn = make_help_button(self, "Input / Output", define_help)
        define_label_col = QWidget(self)
        define_label_layout = QVBoxLayout(define_label_col)
        define_label_layout.setContentsMargins(0, 0, 0, 0)
        define_label_layout.setSpacing(4)
        define_label_layout.setAlignment(Qt.AlignmentFlag.AlignTop | Qt.AlignmentFlag.AlignRight)
        define_label_layout.addWidget(define_label, 0, Qt.AlignmentFlag.AlignRight)
        define_label_layout.addWidget(define_help_btn, 0, Qt.AlignmentFlag.AlignRight)
        define_label_layout.addStretch(1)
        form.addRow(define_label_col, self._vars_scroll)
        add_row = QWidget(self)
        add_layout = QHBoxLayout(add_row)
        add_layout.setContentsMargins(0, 0, 0, 0)
        add_layout.addWidget(add_var_btn)
        add_layout.addStretch(1)
        form.addRow("", add_row)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        save_btn = buttons.button(QDialogButtonBox.StandardButton.Ok)
        if save_btn is not None:
            save_btn.setText("Save")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form, 1)
        layout.addWidget(buttons, 0)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        QTimer.singleShot(0, self._label.setFocus)

    def _on_add_variable_clicked(self) -> None:
        self._add_var_row()

    def _add_var_row(
        self,
        var: CommandVarDef | None = None,
        *,
        role: str = "",
        special: bool = False,
    ) -> None:
        if not hasattr(self, "_rows_host"):
            return
        row = _CommandVarRow(
            var,
            role=role,
            special=special,
            parent=self._rows_host,
        )
        row.remove_requested.connect(self._remove_var_row)
        row.role_changed.connect(self._on_role_changed)
        self._rows_layout.addWidget(row)
        self._rows.append(row)
        self._refresh_remove_buttons()
        if not self._rebuilding_rows:
            QTimer.singleShot(0, lambda r=row: self._scroll_var_row_into_view(r))

    def _remove_var_row(self, row: _CommandVarRow) -> None:
        if row not in self._rows or row.is_special:
            return
        user_rows = [r for r in self._rows if not r.is_special]
        if len(user_rows) <= 1 and row in user_rows:
            row._name.clear()
            row._default.clear()
            row.set_role("")
            row._mode.setCurrentIndex(row._mode.findData("file"))
            row._refresh_role_fields()
            return
        self._rows.remove(row)
        self._rows_layout.removeWidget(row)
        row.setParent(None)
        row.deleteLater()
        self._refresh_remove_buttons()

    def _scroll_var_row_into_view(self, row: _CommandVarRow) -> None:
        if not hasattr(self, "_vars_scroll"):
            return
        self._vars_scroll.ensureWidgetVisible(row, 0, 6)

    def _refresh_remove_buttons(self) -> None:
        user_rows = [r for r in self._rows if not r.is_special]
        removable = len(user_rows) > 1
        for row in self._rows:
            row.mark_removable(removable)

    def _on_role_changed(self, source_row: _CommandVarRow, role: str) -> None:
        # Only one Input; multiple Outputs are allowed.
        if role != "input":
            return
        for row in self._rows:
            if row is source_row:
                continue
            other, other_role = row.definition()
            del other
            if other_role == "input":
                row.set_role("")

    def _detect_variables(self) -> None:
        current_defs: dict[str, CommandVarDef] = {}
        roles: dict[str, str] = {}
        for row in self._rows:
            var, role = row.definition()
            if var is None:
                continue
            if role:
                roles[var.name] = role
            if is_special_variable(var.name) and not is_affixed_variable(var.name):
                continue
            current_defs[var.name] = var
        cmd = self._command.toPlainText()
        names_in_cmd = template_variables(cmd)
        specials = [n for n in names_in_cmd if is_special_variable(n)]
        affixed = template_affixed_variables(cmd)
        for name in affixed:
            if name not in current_defs:
                current_defs[name] = CommandVarDef(
                    name=name, mode=special_variable_mode(name), value=""
                )
        for name in names_in_cmd:
            if is_special_variable(name):
                continue
            if name not in current_defs:
                current_defs[name] = CommandVarDef(name=name, mode="file", value="")
        affixed_set = set(affixed)
        other_names = [n for n in sorted(current_defs) if n not in affixed_set]
        self._rebuilding_rows = True
        try:
            for row in list(self._rows):
                self._rows_layout.removeWidget(row)
                row.setParent(None)
                row.deleteLater()
            self._rows = []
            for name in specials:
                self._add_var_row(
                    CommandVarDef(name=name, mode=special_variable_mode(name), value=""),
                    role=roles.get(name, ""),
                    special=True,
                )
            for name in affixed:
                self._add_var_row(current_defs[name], role=roles.get(name, ""))
            for name in other_names:
                self._add_var_row(current_defs[name], role=roles.get(name, ""))
        finally:
            self._rebuilding_rows = False
        if not self._rows:
            self._add_var_row()
        if self._rows:
            QTimer.singleShot(0, lambda: self._scroll_var_row_into_view(self._rows[-1]))

    def sizeHint(self) -> QSize:  # noqa: N802
        hint = super().sizeHint()
        return QSize(
            max(hint.width(), self._DEFAULT_WIDTH),
            max(hint.height(), self._DEFAULT_HEIGHT),
        )

    def spec(self) -> UserCommandPluginSpec:
        inputs: list[str] = []
        outputs: list[str] = []
        var_defs: list[CommandVarDef] = []
        for row in self._rows:
            var, role = row.definition()
            if var is None:
                continue
            if is_special_variable(var.name) and not is_affixed_variable(var.name):
                if role == "input":
                    inputs.append(var.name)
                elif role == "output":
                    outputs.append(var.name)
                continue
            var_defs.append(var)
            if role == "input":
                inputs.append(var.name)
            elif role == "output":
                outputs.append(var.name)
        return UserCommandPluginSpec(
            id=plugin_id_from_label(self._label.text().strip()),
            label=self._label.text().strip(),
            command_template=self._command.toPlainText().strip(),
            input_vars=inputs[:1],
            output_vars=outputs,
            var_defs=var_defs,
            fasta_only=self._initial_fasta_only,
            tooltip=self._tooltip.toPlainText().strip(),
            show_on_toolbar=False,
        )

    def accept(self) -> None:
        spec = self.spec()
        if not spec.label:
            QMessageBox.warning(self, "Design command wrapper", "Wrapper name is required.")
            return
        label_key = spec.label.strip().lower()
        for existing in load_user_command_plugin_specs():
            if existing.id == self._initial_id:
                continue
            if existing.label.strip().lower() == label_key:
                QMessageBox.warning(
                    self,
                    "Design command wrapper",
                    f"A wrapper named “{spec.label}” already exists. Choose a different name.",
                )
                return
        if not spec.command_template:
            QMessageBox.warning(
                self, "Design command wrapper", "Command template is required."
            )
            return
        input_count = 0
        for row in self._rows:
            var, role = row.definition()
            if var is None:
                continue
            if role in {"input", "output"}:
                if normalize_var_mode(var.mode) not in {"file", "folder"}:
                    QMessageBox.warning(
                        self,
                        "Design command wrapper",
                        f"{role.title()} variable “{display_var_name(var.name)}” "
                        "must use type File or Folder.",
                    )
                    return
                if role == "input":
                    input_count += 1
        if input_count > 1:
            QMessageBox.warning(
                self,
                "Design command wrapper",
                "Only one variable can be marked as Input.",
            )
            return
        super().accept()


class WrapperManageDialog(QDialog):
    """List, create, edit, or delete command wrappers."""

    def __init__(
        self,
        *,
        on_save,
        on_delete,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Manage wrapper")
        self.setMinimumWidth(640)
        self.setMinimumHeight(400)
        self._on_save = on_save
        self._on_delete = on_delete

        self._list = QListWidget(self)
        self._list.currentItemChanged.connect(self._update_buttons)

        self._add_btn = QPushButton("New", self)
        self._edit_btn = QPushButton("Edit", self)
        self._duplicate_btn = QPushButton("Duplicate", self)
        self._delete_btn = QPushButton("Delete", self)
        self._add_btn.clicked.connect(self._add_wrapper)
        self._edit_btn.clicked.connect(self._edit_wrapper)
        self._duplicate_btn.clicked.connect(self._duplicate_wrapper)
        self._delete_btn.clicked.connect(self._delete_wrapper)
        self._list.itemDoubleClicked.connect(lambda _item: self._edit_wrapper())

        btn_row = QHBoxLayout()
        btn_row.addWidget(self._add_btn)
        btn_row.addWidget(self._edit_btn)
        btn_row.addWidget(self._duplicate_btn)
        btn_row.addWidget(self._delete_btn)
        btn_row.addStretch(1)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.reject)
        close_btn = buttons.button(QDialogButtonBox.StandardButton.Close)
        if close_btn is not None:
            close_btn.clicked.connect(self.accept)

        layout = QVBoxLayout(self)
        layout.addWidget(self._list, 1)
        layout.addLayout(btn_row)
        layout.addWidget(buttons)

        self._populate()
        self._update_buttons()

    def _populate(self) -> None:
        self._list.clear()
        specs = load_user_command_plugin_specs()
        if not specs:
            item = QListWidgetItem("(No command wrappers yet — click New)")
            item.setFlags(Qt.ItemFlag.NoItemFlags)
            self._list.addItem(item)
            return
        for index, spec in enumerate(specs, start=1):
            item = QListWidgetItem(f"{index}. {spec.label}")
            item.setData(Qt.ItemDataRole.UserRole, spec.id)
            if spec.tooltip:
                item.setToolTip(spec.tooltip)
            self._list.addItem(item)

    def _selected_id(self) -> str | None:
        item = self._list.currentItem()
        if item is None:
            return None
        plugin_id = item.data(Qt.ItemDataRole.UserRole)
        return str(plugin_id) if plugin_id else None

    def _update_buttons(self) -> None:
        has = self._selected_id() is not None
        self._edit_btn.setEnabled(has)
        self._duplicate_btn.setEnabled(has)
        self._delete_btn.setEnabled(has)

    def _add_wrapper(self) -> None:
        dialog = CommandPluginDesignerDialog(parent=self)
        if not exec_modal_dialog(dialog):
            return
        self._on_save(dialog.spec())
        self._populate()

    def _edit_wrapper(self) -> None:
        plugin_id = self._selected_id()
        if not plugin_id:
            return
        spec = get_user_command_plugin_spec(plugin_id)
        if spec is None:
            QMessageBox.warning(self, "Manage wrapper", f"Could not load wrapper {plugin_id!r}.")
            return
        dialog = CommandPluginDesignerDialog(initial=spec, parent=self)
        if not exec_modal_dialog(dialog):
            return
        self._on_save(dialog.spec(), previous_id=plugin_id)
        self._populate()

    @staticmethod
    def _unique_copy_label(base_label: str) -> str:
        existing = {
            spec.label.strip().lower() for spec in load_user_command_plugin_specs()
        }
        base = (base_label or "Wrapper").strip() or "Wrapper"
        candidate = f"{base} copy"
        if candidate.lower() not in existing:
            return candidate
        index = 2
        while f"{base} copy {index}".lower() in existing:
            index += 1
        return f"{base} copy {index}"

    def _duplicate_wrapper(self) -> None:
        plugin_id = self._selected_id()
        if not plugin_id:
            return
        spec = get_user_command_plugin_spec(plugin_id)
        if spec is None:
            QMessageBox.warning(
                self, "Manage wrapper", f"Could not load wrapper {plugin_id!r}."
            )
            return
        new_label = self._unique_copy_label(spec.label)
        new_id = plugin_id_from_label(new_label)
        existing_ids = {item.id for item in load_user_command_plugin_specs()}
        if new_id in existing_ids:
            stem = new_id
            index = 2
            while f"{stem}_{index}" in existing_ids:
                index += 1
            new_id = f"{stem}_{index}"
        duplicate = UserCommandPluginSpec(
            id=new_id,
            label=new_label,
            command_template=spec.command_template,
            output_vars=list(spec.output_vars),
            var_defs=[
                CommandVarDef(name=var.name, mode=var.mode, value=var.value)
                for var in spec.var_defs
            ],
            fasta_only=spec.fasta_only,
            tooltip=spec.tooltip,
            show_on_toolbar=False,
            input_vars=list(spec.input_vars or []),
        )
        self._on_save(duplicate)
        self._populate()
        for row in range(self._list.count()):
            item = self._list.item(row)
            if item is not None and item.data(Qt.ItemDataRole.UserRole) == new_id:
                self._list.setCurrentItem(item)
                break
        self._update_buttons()

    def _delete_wrapper(self) -> None:
        plugin_id = self._selected_id()
        if not plugin_id:
            return
        spec = get_user_command_plugin_spec(plugin_id)
        label = spec.label if spec is not None else plugin_id
        answer = QMessageBox.question(
            self,
            "Manage wrapper",
            f"Delete command wrapper “{label}”?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        if not self._on_delete(plugin_id):
            QMessageBox.warning(self, "Manage wrapper", f"Could not delete wrapper {plugin_id!r}.")
            return
        self._populate()
        self._update_buttons()


# Backward-compatible alias (older call sites / docs).
PluginManageDialog = WrapperManageDialog


class AlignParamsDialog(QDialog):
    def __init__(
        self,
        *,
        initial: AlignOptions | None = None,
        use_system_mafft: bool = False,
        use_system_muscle: bool = False,
        fastga_initial=None,
        fastga_files: list[Path] | None = None,
        fastga_default_reference: Path | None = None,
        default_msa: Path | None = None,
        initial_tab: str | None = None,
        title: str = "Alignment parameters",
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        from repeatviewer.job_options import FastGaOptions

        self.setWindowTitle(title)
        self._for_run = for_run
        opts = initial or AlignOptions()
        self._initial_align = opts
        self._default_msa = self._resolve_existing_align_path(default_msa)
        fastga_opts = fastga_initial or FastGaOptions()
        self._fastga_panel = FastGaPanel(
            initial=fastga_opts,
            files=fastga_files,
            default_reference=fastga_default_reference,
            parent=self,
        )

        self._strategy = QComboBox(self)
        for label, value in (
            ("Auto (MAFFT chooses)", MAFFT_STRATEGY_AUTO),
            ("FFT-NS-1 (very fast progressive)", MAFFT_STRATEGY_FFT_NS_1),
            ("FFT-NS-2 (fast progressive)", MAFFT_STRATEGY_FFT_NS_2),
            ("L-INS-i (local homology)", MAFFT_STRATEGY_LINSI),
            ("G-INS-i (global homology)", MAFFT_STRATEGY_GINSI),
            ("E-INS-i (large internal gaps)", MAFFT_STRATEGY_EINSI),
            ("Custom", MAFFT_STRATEGY_CUSTOM),
        ):
            self._strategy.addItem(label, value)
        idx = self._strategy.findData(opts.strategy)
        if idx >= 0:
            self._strategy.setCurrentIndex(idx)

        self._sequence_type = QComboBox(self)
        self._sequence_type.addItem("Auto-detect", MAFFT_SEQUENCE_AUTO)
        self._sequence_type.addItem("Nucleotide", MAFFT_SEQUENCE_NUCLEOTIDE)
        self._sequence_type.addItem("Protein", MAFFT_SEQUENCE_PROTEIN)
        idx = self._sequence_type.findData(opts.sequence_type)
        if idx >= 0:
            self._sequence_type.setCurrentIndex(idx)

        self._threads = QSpinBox(self)
        self._threads.setRange(1, 64)
        self._threads.setValue(max(1, opts.num_threads))

        self._retree = QSpinBox(self)
        self._retree.setRange(1, 100)
        self._retree.setValue(opts.retree)
        self._max_iterate = QSpinBox(self)
        self._max_iterate.setRange(0, 1000)
        self._max_iterate.setValue(opts.max_iterate)
        self._custom_retree = opts.retree
        self._custom_max_iterate = opts.max_iterate
        self._last_strategy: str | None = None

        self._gap_open = QDoubleSpinBox(self)
        self._gap_open.setRange(0.0, 20.0)
        self._gap_open.setDecimals(3)
        self._gap_open.setSingleStep(0.1)
        self._gap_open.setValue(opts.gap_open_penalty)
        self._offset = QDoubleSpinBox(self)
        self._offset.setRange(-10.0, 10.0)
        self._offset.setDecimals(3)
        self._offset.setSingleStep(0.05)
        self._offset.setValue(opts.offset)

        self._reorder = QCheckBox(
            "--reorder (group similar sequences in the output)", self
        )
        self._reorder.setChecked(opts.reorder)
        self._adjust_direction = QCheckBox(
            "--adjustdirectionaccurately (reverse-complement nucleotide sequences)",
            self,
        )
        self._adjust_direction.setChecked(opts.adjust_direction)

        for field in (
            self._strategy,
            self._sequence_type,
            self._threads,
            self._retree,
            self._max_iterate,
            self._gap_open,
            self._offset,
        ):
            field.setFixedWidth(260)

        self._mafft_combine = QRadioButton(
            "Combine selected sequences and align", self
        )
        self._mafft_separate = QRadioButton(
            "Align each selected file separately", self
        )
        self._mafft_add_to_msa = QRadioButton(
            "Add sequences to an existing MSA", self
        )
        self._mafft_mode_group = QButtonGroup(self)
        self._mafft_mode_group.addButton(self._mafft_combine)
        self._mafft_mode_group.addButton(self._mafft_separate)
        self._mafft_mode_group.addButton(self._mafft_add_to_msa)
        self._set_mafft_run_mode(
            getattr(opts, "mafft_run_mode", None)
            or (
                ALIGN_RUN_MODE_COMBINED
                if opts.mafft_combined
                else ALIGN_RUN_MODE_SEPARATE
            )
        )

        self._mafft_input_group = QButtonGroup(self)
        self._mafft_input_selected = QRadioButton("Use selected sequences", self)
        self._mafft_input_define = QRadioButton("Define inputs", self)
        self._mafft_input_group.addButton(self._mafft_input_selected)
        self._mafft_input_group.addButton(self._mafft_input_define)
        self._mafft_input_selected.setChecked(True)
        # Define-input paths are session-only; do not restore last run's paths.
        self._mafft_define_inputs = TeMultiInputFilesWidget(
            self,
            initial=None,
            start_empty=False,
            placeholder="Click in Files/Folders or Browse…",
            browse_title="Choose a FASTA file",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
            ),
            add_button_label="Add another input",
            add_button_tooltip=(
                "Add another FASTA file path. Click a file in the Files panel "
                "to fill a row, or use Browse."
            ),
            path_mode="file",
        )
        self._mafft_add_sequences = PathBrowseRow(
            self,
            placeholder="Click in Files, Browse, or Paste…",
            browse_title="Choose sequences to add (FASTA)",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
            ),
            minimum_width=280,
            extra_button_label="Paste…",
            extra_button_tooltip=(
                "Type or paste a sequence; saves a FASTA and fills this path."
            ),
        )
        self._mafft_add_sequences.extra_clicked.connect(
            lambda: self._paste_sequence_to_add(self._mafft_add_sequences)
        )
        self._mafft_add_msa = PathBrowseRow(
            self,
            initial=str(self._default_msa) if self._default_msa else None,
            placeholder="Click in Files or Browse…",
            browse_title="Choose the existing MSA (FASTA)",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
            ),
            minimum_width=280,
        )
        for row in (self._mafft_add_sequences, self._mafft_add_msa):
            row.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            row.line_edit().setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
            )

        self._mafft_add_method = QComboBox(self)
        for label, value in (
            ("--add (insert into MSA)", MAFFT_ADD_METHOD_ADD),
            ("--addfragments (short / partial)", MAFFT_ADD_METHOD_ADDFRAGMENTS),
            ("--addfull (full-length into fragmentary MSA)", MAFFT_ADD_METHOD_ADDFULL),
            ("--addlong (long sequences)", MAFFT_ADD_METHOD_ADDLONG),
        ):
            self._mafft_add_method.addItem(label, value)
        add_method = getattr(opts, "mafft_add_method", MAFFT_ADD_METHOD_ADD)
        add_idx = self._mafft_add_method.findData(add_method)
        if add_idx >= 0:
            self._mafft_add_method.setCurrentIndex(add_idx)
        self._mafft_add_method.setFixedWidth(260)
        self._mafft_add_keeplength = QCheckBox(
            "--keeplength (do not extend the MSA)", self
        )
        self._mafft_add_keeplength.setChecked(
            bool(getattr(opts, "mafft_add_keeplength", False))
        )

        mafft_tab = QWidget(self)
        mafft_root = QVBoxLayout(mafft_tab)
        mafft_root.setContentsMargins(16, 10, 16, 10)
        mafft_root.setSpacing(4)
        mafft_basic = QWidget(mafft_tab)
        self._mafft_basic = mafft_basic
        mafft_form = QFormLayout(mafft_basic)
        mafft_form.setContentsMargins(0, 0, 0, 0)
        mafft_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        # ExpandingFieldsGrow so Define-inputs path + Browse fill the dialog width
        # (Settings combos/spins stay fixed width and do not stretch).
        mafft_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        mafft_form.setHorizontalSpacing(12)
        mafft_form.setVerticalSpacing(4)

        def _mafft_section(title: str, parent: QWidget | None = None) -> QLabel:
            return QLabel(f"<b>{title}</b>", parent or mafft_tab)

        mafft_form.addRow(_mafft_section("Mode", mafft_basic))
        mafft_mode_box = QWidget(mafft_basic)
        mafft_mode_layout = QVBoxLayout(mafft_mode_box)
        mafft_mode_layout.setContentsMargins(0, 0, 0, 0)
        mafft_mode_layout.setSpacing(2)
        mafft_mode_layout.addWidget(
            field_with_help(
                self._mafft_combine,
                "MAFFT mode",
                _MAFFT_MODE_HELP,
                stretch=True,
            )
        )
        mafft_mode_layout.addWidget(self._mafft_separate)
        mafft_mode_layout.addWidget(
            field_with_help(
                self._mafft_add_to_msa,
                "Add sequences to an existing MSA",
                _MAFFT_ADD_TO_MSA_HELP,
                stretch=True,
            )
        )
        self._mafft_mode_label = QLabel("")
        mafft_form.addRow(self._mafft_mode_label, mafft_mode_box)

        mafft_form.addRow(_form_section_separator())
        mafft_form.addRow(_mafft_section("Inputs", mafft_basic))
        mafft_inputs_box = QWidget(mafft_basic)
        mafft_inputs_layout = QVBoxLayout(mafft_inputs_box)
        mafft_inputs_layout.setContentsMargins(0, 0, 0, 0)
        mafft_inputs_layout.setSpacing(2)
        mafft_inputs_layout.addWidget(
            field_with_help(
                self._mafft_input_selected,
                "MAFFT inputs",
                _MAFFT_INPUTS_HELP,
                stretch=True,
            )
        )
        mafft_inputs_layout.addWidget(self._mafft_input_define)
        self._mafft_inputs_radios_label = QLabel("")
        self._mafft_inputs_define_label = QLabel("")
        mafft_form.addRow(self._mafft_inputs_radios_label, mafft_inputs_box)
        mafft_form.addRow(self._mafft_inputs_define_label, self._mafft_define_inputs)
        self._mafft_add_sequences_label = QLabel("Sequences to add:", mafft_basic)
        self._mafft_add_msa_label = QLabel("Existing MSA:", mafft_basic)
        self._mafft_add_sequences_field = self._align_path_field_with_help_spacer(
            self._mafft_add_sequences, mafft_basic
        )
        self._mafft_add_msa_field = self._align_path_field_with_help_spacer(
            self._mafft_add_msa, mafft_basic
        )
        mafft_form.addRow(
            self._mafft_add_sequences_label, self._mafft_add_sequences_field
        )
        mafft_form.addRow(self._mafft_add_msa_label, self._mafft_add_msa_field)

        mafft_form.addRow(_form_section_separator())
        mafft_form.addRow(_mafft_section("Settings", mafft_basic))
        mafft_form.addRow(
            "Strategy:",
            field_with_help(self._strategy, "MAFFT strategy", _MAFFT_STRATEGY_HELP),
        )
        mafft_form.addRow(
            "Sequence type:",
            field_with_help(
                self._sequence_type, "Sequence type", _MAFFT_SEQUENCE_TYPE_HELP
            ),
        )
        mafft_form.addRow(
            "Threads:",
            field_with_help(self._threads, "Threads", _MAFFT_THREADS_HELP),
        )
        mafft_root.addWidget(mafft_basic)

        self._mafft_advanced_panel = QWidget(mafft_tab)
        self._mafft_advanced_panel.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        mafft_advanced_form = QFormLayout(self._mafft_advanced_panel)
        mafft_advanced_form.setContentsMargins(0, 4, 0, 0)
        mafft_advanced_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        mafft_advanced_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        mafft_advanced_form.setHorizontalSpacing(12)
        mafft_advanced_form.setVerticalSpacing(4)
        mafft_advanced_form.addRow(_form_section_separator())
        mafft_advanced_form.addRow(
            _mafft_section("Advanced Settings", self._mafft_advanced_panel)
        )
        mafft_advanced_form.addRow(
            "--retree:",
            field_with_help(self._retree, "Guide-tree rebuilds", _MAFFT_RETREE_HELP),
        )
        mafft_advanced_form.addRow(
            "--maxiterate:",
            field_with_help(
                self._max_iterate, "Refinement cycles", _MAFFT_MAXITERATE_HELP
            ),
        )
        mafft_advanced_form.addRow(
            "--op:",
            field_with_help(self._gap_open, "Gap-opening penalty", _MAFFT_OP_HELP),
        )
        mafft_advanced_form.addRow(
            "--ep:",
            field_with_help(self._offset, "Offset value", _MAFFT_EP_HELP),
        )
        mafft_advanced_form.addRow(
            "",
            field_with_help(self._reorder, "Reorder", _MAFFT_REORDER_HELP),
        )
        mafft_advanced_form.addRow(
            "",
            field_with_help(
                self._adjust_direction,
                "Adjust direction",
                _MAFFT_ADJUST_DIRECTION_HELP,
            ),
        )
        mafft_advanced_form.addRow(
            "Add method:",
            field_with_help(
                self._mafft_add_method, "MAFFT add method", _MAFFT_ADD_METHOD_HELP
            ),
        )
        mafft_advanced_form.addRow(
            "",
            field_with_help(
                self._mafft_add_keeplength,
                "Keep MSA length",
                _MAFFT_ADD_KEEPLENGTH_HELP,
            ),
        )
        mafft_root.addWidget(self._mafft_advanced_panel)
        self._mafft_tab = mafft_tab
        self._set_mafft_advanced_visible(False)
        self._mafft_combine.toggled.connect(self._on_mafft_run_mode_changed)
        self._mafft_separate.toggled.connect(self._on_mafft_run_mode_changed)
        self._mafft_add_to_msa.toggled.connect(self._on_mafft_run_mode_changed)
        self._mafft_input_selected.toggled.connect(self._sync_mafft_inputs)
        self._mafft_input_define.toggled.connect(self._sync_mafft_inputs)
        self._mafft_define_inputs.paths_changed.connect(
            self._on_mafft_define_paths_changed
        )
        self._mafft_add_sequences.path_changed.connect(
            self._on_mafft_define_paths_changed
        )
        self._mafft_add_msa.path_changed.connect(self._on_mafft_define_paths_changed)
        self._sync_mafft_inputs()
        mafft_tab.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )

        self._muscle_mode = QComboBox(self)
        self._muscle_mode.addItem("Align (high accuracy)", MUSCLE_MODE_ALIGN)
        self._muscle_mode.addItem("Super5 (large datasets)", MUSCLE_MODE_SUPER5)
        muscle_mode_index = self._muscle_mode.findData(opts.muscle_mode)
        if muscle_mode_index >= 0:
            self._muscle_mode.setCurrentIndex(muscle_mode_index)
        self._muscle_mode.setFixedWidth(260)
        self._muscle_threads = QSpinBox(self)
        self._muscle_threads.setRange(1, 64)
        self._muscle_threads.setValue(max(1, opts.muscle_num_threads))
        self._muscle_threads.setFixedWidth(260)

        self._muscle_perm = QComboBox(self)
        for label, value in (
            ("None (default)", MUSCLE_PERM_NONE),
            ("abc", MUSCLE_PERM_ABC),
            ("acb", MUSCLE_PERM_ACB),
            ("bca", MUSCLE_PERM_BCA),
        ):
            self._muscle_perm.addItem(label, value)
        muscle_perm_index = self._muscle_perm.findData(opts.muscle_perm)
        if muscle_perm_index >= 0:
            self._muscle_perm.setCurrentIndex(muscle_perm_index)
        self._muscle_perm.setFixedWidth(260)
        self._muscle_perturb = QSpinBox(self)
        self._muscle_perturb.setRange(0, 1_000_000)
        self._muscle_perturb.setValue(max(0, int(opts.muscle_perturb)))
        self._muscle_perturb.setSpecialValueText("Off (0)")
        self._muscle_perturb.setFixedWidth(260)

        self._muscle_combine = QRadioButton(
            "Combine selected sequences and align", self
        )
        self._muscle_separate = QRadioButton(
            "Align each selected file separately", self
        )
        self._muscle_add_to_msa = QRadioButton(
            "Add sequences to an existing MSA", self
        )
        self._muscle_mode_group = QButtonGroup(self)
        self._muscle_mode_group.addButton(self._muscle_combine)
        self._muscle_mode_group.addButton(self._muscle_separate)
        self._muscle_mode_group.addButton(self._muscle_add_to_msa)
        self._set_muscle_run_mode(
            getattr(opts, "muscle_run_mode", None)
            or (
                ALIGN_RUN_MODE_COMBINED
                if opts.muscle_combined
                else ALIGN_RUN_MODE_SEPARATE
            )
        )

        self._muscle_input_group = QButtonGroup(self)
        self._muscle_input_selected = QRadioButton("Use selected sequences", self)
        self._muscle_input_define = QRadioButton("Define inputs", self)
        self._muscle_input_group.addButton(self._muscle_input_selected)
        self._muscle_input_group.addButton(self._muscle_input_define)
        # Always open on Use selected (Define is for this run only).
        self._muscle_input_selected.setChecked(True)
        self._muscle_define_inputs = TeMultiInputFilesWidget(
            self,
            initial=None,
            start_empty=False,
            placeholder="Click in Files/Folders or Browse…",
            browse_title="Choose a FASTA file",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
            ),
            add_button_label="Add another input",
            add_button_tooltip=(
                "Add another FASTA file path. Click a file in the Files panel "
                "to fill a row, or use Browse."
            ),
            path_mode="file",
        )
        self._muscle_add_sequences = PathBrowseRow(
            self,
            placeholder="Click in Files, Browse, or Paste…",
            browse_title="Choose sequences to add (FASTA)",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
            ),
            minimum_width=280,
            extra_button_label="Paste…",
            extra_button_tooltip=(
                "Type or paste a sequence; saves a FASTA and fills this path."
            ),
        )
        self._muscle_add_sequences.extra_clicked.connect(
            lambda: self._paste_sequence_to_add(self._muscle_add_sequences)
        )
        self._muscle_add_msa = PathBrowseRow(
            self,
            initial=str(self._default_msa) if self._default_msa else None,
            placeholder="Click in Files or Browse…",
            browse_title="Choose the existing MSA (FASTA)",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
            ),
            minimum_width=280,
        )
        for row in (self._muscle_add_sequences, self._muscle_add_msa):
            row.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
            row.line_edit().setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
            )

        muscle_tab = QWidget(self)
        muscle_root = QVBoxLayout(muscle_tab)
        muscle_root.setContentsMargins(16, 10, 16, 10)
        muscle_root.setSpacing(4)
        muscle_basic = QWidget(muscle_tab)
        self._muscle_basic = muscle_basic
        muscle_form = QFormLayout(muscle_basic)
        muscle_form.setContentsMargins(0, 0, 0, 0)
        muscle_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        # Same as MAFFT: Define-inputs path rows grow; Settings stay fixed width.
        muscle_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        muscle_form.setHorizontalSpacing(12)
        muscle_form.setVerticalSpacing(4)

        def _muscle_section(title: str, parent: QWidget | None = None) -> QLabel:
            return QLabel(f"<b>{title}</b>", parent or muscle_tab)

        muscle_form.addRow(_muscle_section("Mode", muscle_basic))
        muscle_mode_box = QWidget(muscle_basic)
        muscle_mode_layout = QVBoxLayout(muscle_mode_box)
        muscle_mode_layout.setContentsMargins(0, 0, 0, 0)
        muscle_mode_layout.setSpacing(2)
        muscle_mode_layout.addWidget(
            field_with_help(
                self._muscle_combine,
                "MUSCLE mode",
                _MUSCLE_MODE_HELP,
                stretch=True,
            )
        )
        muscle_mode_layout.addWidget(self._muscle_separate)
        muscle_mode_layout.addWidget(
            field_with_help(
                self._muscle_add_to_msa,
                "Add sequences to an existing MSA",
                _MUSCLE_ADD_TO_MSA_HELP,
                stretch=True,
            )
        )
        self._muscle_mode_label = QLabel("")
        muscle_form.addRow(self._muscle_mode_label, muscle_mode_box)

        muscle_form.addRow(_form_section_separator())
        muscle_form.addRow(_muscle_section("Inputs", muscle_basic))
        muscle_inputs_box = QWidget(muscle_basic)
        muscle_inputs_layout = QVBoxLayout(muscle_inputs_box)
        muscle_inputs_layout.setContentsMargins(0, 0, 0, 0)
        muscle_inputs_layout.setSpacing(2)
        muscle_inputs_layout.addWidget(
            field_with_help(
                self._muscle_input_selected,
                "MUSCLE inputs",
                _MUSCLE_INPUTS_HELP,
                stretch=True,
            )
        )
        muscle_inputs_layout.addWidget(self._muscle_input_define)
        self._muscle_inputs_radios_label = QLabel("")
        self._muscle_inputs_define_label = QLabel("")
        muscle_form.addRow(self._muscle_inputs_radios_label, muscle_inputs_box)
        muscle_form.addRow(self._muscle_inputs_define_label, self._muscle_define_inputs)
        self._muscle_add_sequences_label = QLabel("Sequences to add:", muscle_basic)
        self._muscle_add_msa_label = QLabel("Existing MSA:", muscle_basic)
        self._muscle_add_sequences_field = self._align_path_field_with_help_spacer(
            self._muscle_add_sequences, muscle_basic
        )
        self._muscle_add_msa_field = self._align_path_field_with_help_spacer(
            self._muscle_add_msa, muscle_basic
        )
        muscle_form.addRow(
            self._muscle_add_sequences_label, self._muscle_add_sequences_field
        )
        muscle_form.addRow(self._muscle_add_msa_label, self._muscle_add_msa_field)

        muscle_form.addRow(_form_section_separator())
        muscle_form.addRow(_muscle_section("Settings", muscle_basic))
        muscle_form.addRow(
            "Algorithm:",
            field_with_help(
                self._muscle_mode, "MUSCLE algorithm", _MUSCLE_ALGORITHM_HELP
            ),
        )
        muscle_form.addRow(
            "Threads:",
            field_with_help(
                self._muscle_threads, "Threads", _MUSCLE_THREADS_HELP
            ),
        )
        muscle_root.addWidget(muscle_basic)

        self._muscle_advanced_panel = QWidget(muscle_tab)
        self._muscle_advanced_panel.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        muscle_advanced_form = QFormLayout(self._muscle_advanced_panel)
        muscle_advanced_form.setContentsMargins(0, 4, 0, 0)
        muscle_advanced_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        muscle_advanced_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        muscle_advanced_form.setHorizontalSpacing(12)
        muscle_advanced_form.setVerticalSpacing(4)
        muscle_advanced_form.addRow(_form_section_separator())
        muscle_advanced_form.addRow(
            _muscle_section("Advanced Settings", self._muscle_advanced_panel)
        )
        muscle_advanced_form.addRow(
            "Guide tree permutation:",
            field_with_help(
                self._muscle_perm, "Guide tree permutation", _MUSCLE_PERM_HELP
            ),
        )
        muscle_advanced_form.addRow(
            "Perturbation seed:",
            field_with_help(
                self._muscle_perturb, "Perturbation seed", _MUSCLE_PERTURB_HELP
            ),
        )
        muscle_root.addWidget(self._muscle_advanced_panel)
        self._muscle_tab = muscle_tab
        self._set_muscle_advanced_visible(False)
        self._muscle_combine.toggled.connect(self._on_muscle_run_mode_changed)
        self._muscle_separate.toggled.connect(self._on_muscle_run_mode_changed)
        self._muscle_add_to_msa.toggled.connect(self._on_muscle_run_mode_changed)
        self._muscle_input_selected.toggled.connect(self._sync_muscle_inputs)
        self._muscle_input_define.toggled.connect(self._sync_muscle_inputs)
        self._muscle_define_inputs.paths_changed.connect(
            self._on_muscle_define_paths_changed
        )
        self._muscle_add_sequences.path_changed.connect(
            self._on_muscle_define_paths_changed
        )
        self._muscle_add_msa.path_changed.connect(self._on_muscle_define_paths_changed)
        self._sync_muscle_inputs()
        muscle_tab.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )

        self._use_system_mafft = QCheckBox(
            system_path_executable_label("MAFFT"), self
        )
        self._use_system_mafft.setChecked(use_system_mafft)
        self._use_system_mafft.setToolTip(
            "Use mafft from PATH or the conda env in File → Settings → External tools "
            "instead of the bundled MAFFT."
        )
        self._use_system_muscle = QCheckBox(
            system_path_executable_label("MUSCLE"), self
        )
        self._use_system_muscle.setChecked(use_system_muscle)
        self._use_system_muscle.setToolTip(
            "Use muscle from PATH or the conda env in File → Settings → External tools "
            "instead of the bundled MUSCLE."
        )
        self._mafft_path_status = make_executable_status_label(self)
        self._muscle_path_status = make_executable_status_label(self)
        self._fastga_path_status = make_executable_status_label(self)
        self._use_system_mafft.toggled.connect(self._update_executable_status)
        self._use_system_muscle.toggled.connect(self._update_executable_status)

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)

        rc_box = QGroupBox("Right-click action", self)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(
            QLabel("Right-clicking the Align button runs, with saved settings:", self)
        )
        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("MAFFT", ALIGNER_MAFFT),
            ("MUSCLE", ALIGNER_MUSCLE),
            ("FastGA", ALIGNER_FASTGA),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_aligner == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
            rc_box_layout.addWidget(radio)
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[ALIGNER_MAFFT].setChecked(True)
        options_layout.addWidget(rc_box)

        exe_box = QGroupBox("Executables", self)
        exe_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        exe_box_layout = QVBoxLayout(exe_box)
        exe_box_layout.setContentsMargins(10, 8, 10, 8)
        exe_box_layout.setSpacing(8)
        exe_box_layout.addWidget(
            executable_checkbox_row(
                self, self._use_system_mafft, self._mafft_path_status
            )
        )
        exe_box_layout.addWidget(
            executable_checkbox_row(
                self, self._use_system_muscle, self._muscle_path_status
            )
        )
        exe_box_layout.addWidget(self._fastga_path_status)
        options_layout.addWidget(exe_box)

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="the aligner",
                checked=bool(
                    getattr(opts, "show_process_window", True)
                    or getattr(fastga_opts, "show_process_window", True)
                ),
                close_when_finished=bool(
                    getattr(opts, "close_process_window_when_finished", True)
                    or getattr(
                        fastga_opts, "close_process_window_when_finished", True
                    )
                ),
            )
        )
        options_layout.addWidget(process_box)
        from repeatviewer.tool_citations import CITATION_ALIGN

        options_layout.addWidget(make_citation_box(self, CITATION_ALIGN))
        options_layout.addStretch(1)

        self._aligner_tabs = QTabWidget(self)
        self._aligner_tabs.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._mafft_tab = mafft_tab
        self._muscle_tab = muscle_tab
        self._options_tab = options_tab
        # Use the panel directly (no scroll area) so Advanced Settings can grow
        # the dialog height instead of disappearing into an empty viewport.
        # Expanding horizontally so Define-inputs path rows fill like MAFFT/MUSCLE.
        self._fastga_panel.setMinimumWidth(560)
        self._fastga_panel.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )
        self._fastga_tab = self._fastga_panel
        self._fastga_compact_size: QSize | None = None
        self._muscle_compact_size: QSize | None = None
        self._mafft_compact_size: QSize | None = None
        self._fastga_panel.layout_changed.connect(self._on_fastga_panel_changed)
        self._aligner_tabs.addTab(mafft_tab, "MAFFT")
        self._aligner_tabs.addTab(muscle_tab, "MUSCLE")
        self._aligner_tabs.addTab(self._fastga_panel, "FastGA")
        self._aligner_tabs.addTab(options_tab, "Options")
        self._options_tab_index = self._aligner_tabs.indexOf(options_tab)

        self._close_action = "cancel"
        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._advanced_settings_btn = QPushButton("Advanced Settings ▸", self)
        self._advanced_settings_btn.setCheckable(True)
        self._advanced_settings_btn.setChecked(False)
        self._advanced_settings_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._advanced_settings_btn.toggled.connect(self._on_advanced_settings_toggled)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore, self._advanced_settings_btn],
        )

        # Connect after Advanced button exists — setCurrentWidget emits currentChanged.
        self._aligner_tabs.currentChanged.connect(self._on_aligner_tab_changed)
        tab_for = {
            ALIGNER_MAFFT: mafft_tab,
            ALIGNER_MUSCLE: muscle_tab,
            ALIGNER_FASTGA: self._fastga_panel,
            DIALOG_TAB_OPTIONS: options_tab,
        }
        open_tab = opening_dialog_tab(
            initial_tab
            if initial_tab in tab_for
            else getattr(opts, "last_dialog_tab", ALIGNER_MAFFT),
            allowed=tuple(tab_for.keys()),
            first_tab=ALIGNER_MAFFT,
        )
        self._aligner_tabs.setCurrentWidget(tab_for[open_tab])

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(14)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._aligner_tabs, stretch=0)
        layout.addWidget(self._button_box, stretch=0)

        self._strategy.currentIndexChanged.connect(self._sync_controls)
        self._sequence_type.currentIndexChanged.connect(self._sync_controls)
        self._sync_controls()
        self._update_advanced_settings_button()
        self._sync_align_field_columns()
        self._update_executable_status()
        QTimer.singleShot(0, self._fit_current_align_tab)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._sync_align_field_columns()
        if self._aligner_tabs.currentWidget() is self._options_tab:
            self._update_executable_status()
        QTimer.singleShot(0, self._fit_current_align_tab)

    def _fit_current_align_tab(self) -> None:
        """Shrink / grow the dialog to the active tab’s real content height."""
        if not self.isVisible():
            return
        current = self._aligner_tabs.currentWidget()
        if current is self._options_tab:
            self._fit_options_content()
        elif current is self._fastga_tab:
            self._refit_fastga_after_content_change()
        elif current is self._muscle_tab:
            if not self._muscle_advanced_visible():
                self._resize_for_muscle_advanced(False)
        elif current is self._mafft_tab:
            if not self._mafft_advanced_visible():
                self._resize_for_mafft_advanced(False)

    def _fit_options_content(self) -> None:
        page = self._options_tab
        # Clear any prior height pin so newly shown PATH status can grow sizeHint.
        self._pin_aligner_tab_height(None)
        page.ensurePolished()
        page.adjustSize()
        root = page.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        page_hint = page.minimumSizeHint().expandedTo(page.sizeHint())
        content_h = page_hint.height()
        if root is not None:
            content_h = max(
                content_h,
                root.totalMinimumSize().height(),
                root.totalSizeHint().height(),
            )
        chrome = self._dialog_chrome_height()
        height = chrome + content_h
        width = self._align_content_width()
        self.setMinimumHeight(0)
        self._pin_aligner_tab_height(content_h)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def _sync_align_field_columns(self) -> None:
        """Indent Mode / Inputs radios to the shared Settings field column."""
        mafft_form = (
            self._mafft_basic.layout()
            if hasattr(self, "_mafft_basic")
            else None
        )
        mafft_adv = (
            self._mafft_advanced_panel.layout()
            if hasattr(self, "_mafft_advanced_panel")
            else None
        )
        sync_form_field_column_spacers(
            self,
            [
                self._mafft_mode_label,
                self._mafft_inputs_radios_label,
                self._mafft_inputs_define_label,
                self._mafft_add_sequences_label,
                self._mafft_add_msa_label,
            ],
            [
                mafft_form if isinstance(mafft_form, QFormLayout) else None,
                mafft_adv if isinstance(mafft_adv, QFormLayout) else None,
            ],
        )
        muscle_form = (
            self._muscle_basic.layout()
            if hasattr(self, "_muscle_basic")
            else None
        )
        muscle_adv = (
            self._muscle_advanced_panel.layout()
            if hasattr(self, "_muscle_advanced_panel")
            else None
        )
        sync_form_field_column_spacers(
            self,
            [
                self._muscle_mode_label,
                self._muscle_inputs_radios_label,
                self._muscle_inputs_define_label,
                self._muscle_add_sequences_label,
                self._muscle_add_msa_label,
            ],
            [
                muscle_form if isinstance(muscle_form, QFormLayout) else None,
                muscle_adv if isinstance(muscle_adv, QFormLayout) else None,
            ],
        )
        if hasattr(self, "_fastga_panel") and hasattr(
            self._fastga_panel, "_sync_field_column"
        ):
            self._fastga_panel._sync_field_column()

    def _remind_fastga_if_needed(self) -> None:
        """Show install/bind instructions when FastGA is not configured."""
        from repeatviewer.backends.fastga_dotplot import (
            fastga_available,
            show_fastga_not_configured_dialog,
        )

        if fastga_available():
            return
        show_fastga_not_configured_dialog(self)

    def _on_aligner_tab_changed(self, _index: int) -> None:
        self._update_advanced_settings_button()
        self._sync_align_field_columns()
        current = self._aligner_tabs.currentWidget()
        if current is self._options_tab:
            self._update_executable_status()
        if current is self._fastga_tab:
            # Defer so the FastGA tab paints before any availability probe / reminder.
            QTimer.singleShot(0, self._remind_fastga_if_needed)
        if self.isVisible():
            QTimer.singleShot(0, self._fit_current_align_tab)

    def _on_advanced_settings_toggled(self, checked: bool) -> None:
        current = self._aligner_tabs.currentWidget()
        self._advanced_settings_btn.setText(
            "Advanced Settings ▾" if checked else "Advanced Settings ▸"
        )
        if current is self._fastga_tab:
            if checked:
                self._fastga_compact_size = self.size()
            self._fastga_panel.set_advanced_visible(checked)
            self._resize_for_fastga_advanced(checked)
        elif current is self._muscle_tab:
            if checked:
                self._muscle_compact_size = self.size()
            self._set_muscle_advanced_visible(checked)
            self._resize_for_muscle_advanced(checked)
        elif current is self._mafft_tab:
            if checked:
                self._mafft_compact_size = self.size()
            self._set_mafft_advanced_visible(checked)
            self._resize_for_mafft_advanced(checked)

    def _dialog_chrome_height(self) -> int:
        """Height of margins, tab bar, and footer buttons around the active tab."""
        lay = self.layout()
        margins = lay.contentsMargins() if lay is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 34
        spacing = lay.spacing() if lay is not None else 14
        tab_bar_h = self._aligner_tabs.tabBar().sizeHint().height() + 10
        buttons_h = self._button_box.sizeHint().height()
        return margin_h + spacing + tab_bar_h + buttons_h

    def _align_needs_wide_inputs(self) -> bool:
        """True when the active Align tab shows path rows that need width 720.

        Define inputs (MAFFT / MUSCLE / FastGA), add-to-MSA paths, or FastGA
        Compare-against-a reference (Reference genome path + Browse).
        """
        current = self._aligner_tabs.currentWidget()
        if current is self._mafft_tab:
            return (
                self.mafft_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA
                or self.mafft_input_source() == ALIGN_INPUT_DEFINE
            )
        if current is self._muscle_tab:
            return (
                self.muscle_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA
                or self.muscle_input_source() == ALIGN_INPUT_DEFINE
            )
        if current is self._fastga_tab:
            panel = self._fastga_panel
            return (
                panel.input_source() == FASTGA_INPUT_DEFINE
                or panel.mode() == FASTGA_MODE_REFERENCE
            )
        return False

    def _align_content_width(self) -> int:
        """608 normally; 720 when Define inputs or FastGA reference mode needs path fill."""
        width = tool_params_dialog_width(
            define_inputs=self._align_needs_wide_inputs()
        )
        self.setMinimumWidth(width)
        return width

    def _pin_aligner_tab_height(self, page_h: int | None) -> None:
        """Constrain the tab widget to the active page so other tabs cannot gap."""
        page = self._aligner_tabs.currentWidget()
        if page_h is None:
            self._aligner_tabs.setMinimumHeight(0)
            self._aligner_tabs.setMaximumHeight(16777215)
            if page is not None:
                page.setMinimumHeight(0)
                page.setMaximumHeight(16777215)
            return
        tab_bar_h = self._aligner_tabs.tabBar().sizeHint().height() + 10
        tabs_h = max(1, int(page_h) + tab_bar_h)
        self._aligner_tabs.setMinimumHeight(0)
        self._aligner_tabs.setMaximumHeight(tabs_h)
        if page is not None:
            page.setMinimumHeight(0)
            page.setMaximumHeight(max(1, int(page_h)))

    def _fastga_dialog_chrome_height(self) -> int:
        return self._dialog_chrome_height()

    def _on_fastga_panel_changed(self) -> None:
        """Grow/shrink when FastGA Mode or Inputs reveal extra rows."""
        if not self.isVisible():
            return
        if self._aligner_tabs.currentWidget() is not self._fastga_tab:
            return
        self._fastga_compact_size = None
        QTimer.singleShot(0, self._refit_fastga_after_content_change)

    def _refit_fastga_after_content_change(self) -> None:
        """Refit height after Mode/Inputs rows change.

        Measure the compact (no-Advanced) size first so a leftover Advanced
        height cannot leave a blank band under Settings. Restore Advanced only
        when the footer toggle is still checked.
        """
        if not self.isVisible():
            return
        if self._aligner_tabs.currentWidget() is not self._fastga_tab:
            return
        want_advanced = bool(self._advanced_settings_btn.isChecked())
        # Hide Advanced while measuring so its sizeHint cannot inflate the pin
        # or leave an empty gap when the panel is collapsed.
        self._fastga_panel.set_advanced_visible(False)
        self._fastga_compact_size = None
        self._resize_for_fastga_advanced(False)
        self._fastga_compact_size = self.size()
        if want_advanced:
            self._fastga_panel.set_advanced_visible(True)
            self._resize_for_fastga_advanced(True)
        self._update_advanced_settings_button()

    def _shrink_to_fastga_content(self) -> None:
        """Fit the dialog to the current FastGA form (no unused blank height)."""
        if self._aligner_tabs.currentWidget() is not self._fastga_tab:
            return
        if self._fastga_panel.advanced_visible():
            return
        self._resize_for_fastga_advanced(False)

    def _shrink_to_muscle_content(self) -> None:
        if self._aligner_tabs.currentWidget() is not self._muscle_tab:
            return
        if self._muscle_advanced_visible():
            return
        self._resize_for_muscle_advanced(False)

    def _muscle_advanced_visible(self) -> bool:
        return bool(self._muscle_advanced_panel.isVisible())

    def _set_muscle_advanced_visible(self, visible: bool) -> None:
        self._muscle_advanced_panel.setVisible(visible)
        if visible:
            self._muscle_advanced_panel.setMaximumHeight(16777215)
        else:
            self._muscle_advanced_panel.setMaximumHeight(0)
        self._muscle_advanced_panel.updateGeometry()
        root = self._muscle_tab.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        self._muscle_tab.updateGeometry()

    def _muscle_compact_content_height(self) -> int:
        root = self._muscle_tab.layout()
        margins = root.contentsMargins() if root is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 20
        spacing = root.spacing() if root is not None else 4
        self._muscle_basic.ensurePolished()
        basic_h = max(
            self._muscle_basic.minimumSizeHint().height(),
            self._muscle_basic.sizeHint().height(),
        )
        return margin_h + basic_h + spacing

    def _resize_for_fastga_advanced(self, expanded: bool) -> None:
        """Grow the dialog for Advanced Settings; shrink back when collapsed."""
        panel = self._fastga_panel
        adv = panel._advanced_panel
        # Keep panel visibility aligned with the requested state so a collapsed
        # Advanced block cannot keep dialog height while showing a blank gap.
        if expanded and not panel.advanced_visible():
            panel.set_advanced_visible(True)
        elif not expanded and panel.advanced_visible():
            # Caller may be measuring compact height; leave visibility alone when
            # the footer toggle still wants Advanced (handled by refit).
            pass
        # Unlock prior height pins before measuring (Add Define-input rows).
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        self._pin_aligner_tab_height(None)
        define = panel._define_inputs
        define.setMinimumHeight(0)
        define.setMaximumHeight(16777215)
        basic = panel._basic_panel
        basic.setMinimumHeight(0)
        basic.setMaximumHeight(16777215)
        basic.ensurePolished()
        basic.adjustSize()
        root = panel.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        panel.updateGeometry()
        self._aligner_tabs.updateGeometry()

        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        chrome = self._dialog_chrome_height()
        width = self._align_content_width()
        self.setMinimumWidth(width)

        if not expanded:
            content_h = panel.compact_content_height()
            height = chrome + content_h
            height = min(max(height, 220), max_h)
            self._pin_aligner_tab_height(content_h)
            # Cap then resize so macOS actually shrinks leftover Advanced height.
            self.setMaximumHeight(height)
            self.resize(width, height)
            self.setMaximumHeight(16777215)
            self._fastga_compact_size = QSize(width, height)
            return

        extra = max(0, adv.sizeHint().height())
        content_h = panel.compact_content_height() + extra
        self._pin_aligner_tab_height(content_h)
        base = self._fastga_compact_size or self.size()
        height = min(base.height() + extra + 12, max_h)
        height = max(
            height,
            min(chrome + content_h, max_h),
        )
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def _resize_for_muscle_advanced(self, expanded: bool) -> None:
        panel = self._muscle_tab
        adv = self._muscle_advanced_panel
        # Unlock prior height pins before measuring (Add Define-input rows).
        self._pin_aligner_tab_height(None)
        self._muscle_define_inputs.setMinimumHeight(0)
        self._muscle_define_inputs.setMaximumHeight(16777215)
        root = panel.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        self._muscle_basic.ensurePolished()
        self._muscle_basic.adjustSize()
        panel.updateGeometry()
        self._aligner_tabs.updateGeometry()

        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        chrome = self._dialog_chrome_height()
        width = self._align_content_width()

        if not expanded:
            content_h = self._muscle_compact_content_height()
            height = chrome + content_h
            height = min(max(height, 220), max_h)
            self.setMinimumHeight(0)
            self._aligner_tabs.setMinimumHeight(0)
            panel.setMinimumHeight(0)
            self._pin_aligner_tab_height(content_h)
            self.setMaximumHeight(height)
            self.resize(width, height)
            self.setMaximumHeight(16777215)
            self._muscle_compact_size = QSize(width, height)
            return

        extra = max(0, adv.sizeHint().height())
        content_h = self._muscle_compact_content_height() + extra
        self._pin_aligner_tab_height(content_h)
        base = self._muscle_compact_size or self.size()
        height = min(base.height() + extra + 12, max_h)
        height = max(
            height,
            min(chrome + content_h, max_h),
        )
        self.resize(width, height)

    def _shrink_to_mafft_content(self) -> None:
        if self._aligner_tabs.currentWidget() is not self._mafft_tab:
            return
        if self._mafft_advanced_visible():
            return
        self._resize_for_mafft_advanced(False)

    def _mafft_advanced_visible(self) -> bool:
        return bool(self._mafft_advanced_panel.isVisible())

    def _set_mafft_advanced_visible(self, visible: bool) -> None:
        self._mafft_advanced_panel.setVisible(visible)
        if visible:
            self._mafft_advanced_panel.setMaximumHeight(16777215)
        else:
            self._mafft_advanced_panel.setMaximumHeight(0)
        self._mafft_advanced_panel.updateGeometry()
        root = self._mafft_tab.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        self._mafft_tab.updateGeometry()

    def _mafft_compact_content_height(self) -> int:
        root = self._mafft_tab.layout()
        margins = root.contentsMargins() if root is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 20
        spacing = root.spacing() if root is not None else 4
        self._mafft_basic.ensurePolished()
        basic_h = max(
            self._mafft_basic.minimumSizeHint().height(),
            self._mafft_basic.sizeHint().height(),
        )
        return margin_h + basic_h + spacing

    def _resize_for_mafft_advanced(self, expanded: bool) -> None:
        panel = self._mafft_tab
        adv = self._mafft_advanced_panel
        # Unlock prior height pins before measuring (Add Define-input rows).
        self._pin_aligner_tab_height(None)
        self._mafft_define_inputs.setMinimumHeight(0)
        self._mafft_define_inputs.setMaximumHeight(16777215)
        root = panel.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        self._mafft_basic.ensurePolished()
        self._mafft_basic.adjustSize()
        panel.updateGeometry()
        self._aligner_tabs.updateGeometry()

        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        chrome = self._dialog_chrome_height()
        width = self._align_content_width()

        if not expanded:
            content_h = self._mafft_compact_content_height()
            height = chrome + content_h
            height = min(max(height, 220), max_h)
            self.setMinimumHeight(0)
            self._aligner_tabs.setMinimumHeight(0)
            panel.setMinimumHeight(0)
            self._pin_aligner_tab_height(content_h)
            self.setMaximumHeight(height)
            self.resize(width, height)
            self.setMaximumHeight(16777215)
            self._mafft_compact_size = QSize(width, height)
            return

        extra = max(0, adv.sizeHint().height())
        content_h = self._mafft_compact_content_height() + extra
        self._pin_aligner_tab_height(content_h)
        base = self._mafft_compact_size or self.size()
        height = min(base.height() + extra + 12, max_h)
        height = max(
            height,
            min(chrome + content_h, max_h),
        )
        self.resize(width, height)

    def _update_advanced_settings_button(self) -> None:
        if not hasattr(self, "_advanced_settings_btn"):
            return
        current = self._aligner_tabs.currentWidget()
        on_advanced_tab = current in (
            self._fastga_tab,
            self._muscle_tab,
            self._mafft_tab,
        )
        self._advanced_settings_btn.setVisible(on_advanced_tab)
        if not on_advanced_tab:
            return
        if current is self._fastga_tab:
            expanded = self._fastga_panel.advanced_visible()
        elif current is self._muscle_tab:
            expanded = self._muscle_advanced_visible()
        else:
            expanded = self._mafft_advanced_visible()
        self._advanced_settings_btn.blockSignals(True)
        self._advanced_settings_btn.setChecked(expanded)
        self._advanced_settings_btn.setText(
            "Advanced Settings ▾" if expanded else "Advanced Settings ▸"
        )
        self._advanced_settings_btn.blockSignals(False)

    def should_run(self) -> bool:
        """True when the dialog was closed with Run (not Save or Cancel)."""
        return dialog_should_run(self)

    def selected_aligner(self) -> str:
        """The aligner of the active tab; on the Options tab, the right-click one."""
        current = self._aligner_tabs.currentWidget()
        if current is self._mafft_tab:
            return ALIGNER_MAFFT
        if current is self._muscle_tab:
            return ALIGNER_MUSCLE
        if current is self._fastga_tab:
            return ALIGNER_FASTGA
        return self.right_click_aligner()

    def current_dialog_tab(self) -> str:
        """Active Align window tab key (mafft / muscle / fastga / options)."""
        current = self._aligner_tabs.currentWidget()
        if current is self._mafft_tab:
            return ALIGNER_MAFFT
        if current is self._muscle_tab:
            return ALIGNER_MUSCLE
        if current is self._fastga_tab:
            return ALIGNER_FASTGA
        return DIALOG_TAB_OPTIONS

    def fastga_panel(self) -> "FastGaPanel":
        return self._fastga_panel

    def apply_path_from_panel(self, path: Path) -> bool:
        """Files-panel click: fills Align Define inputs or FastGA reference."""
        current = self._aligner_tabs.currentWidget()
        if current is self._mafft_tab:
            return self._apply_mafft_path_from_panel(path)
        if current is self._muscle_tab:
            return self._apply_muscle_path_from_panel(path)
        if current is self._fastga_tab:
            return self._fastga_panel.apply_path_from_panel(path)
        return False

    def apply_folder_from_panel(self, folder: Path) -> bool:
        """Folders-panel click: fills Align/FastGA Define-inputs path."""
        current = self._aligner_tabs.currentWidget()
        if current is self._mafft_tab:
            return self._apply_mafft_folder_from_panel(folder)
        if current is self._muscle_tab:
            return self._apply_muscle_folder_from_panel(folder)
        if current is self._fastga_tab:
            return self._fastga_panel.apply_folder_from_panel(folder)
        return False

    @staticmethod
    def _resolve_existing_align_path(path: Path | str | None) -> Path | None:
        if path is None or str(path).strip() == "":
            return None
        try:
            candidate = Path(path).expanduser().resolve()
        except OSError:
            return None
        return candidate if candidate.is_file() else None

    @staticmethod
    def _align_path_field_with_help_spacer(
        row: PathBrowseRow, parent: QWidget
    ) -> QWidget:
        """PathBrowseRow + help-sized spacer so Browse lines up with Define inputs."""
        host = QWidget(parent)
        host.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        lay = QHBoxLayout(host)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(6)
        lay.addWidget(row, stretch=1)
        spacer = QWidget(host)
        spacer.setFixedSize(HELP_BUTTON_SIZE, HELP_BUTTON_SIZE)
        lay.addWidget(spacer)
        return host

    def mafft_run_mode(self) -> str:
        if self._mafft_add_to_msa.isChecked():
            return ALIGN_RUN_MODE_ADD_TO_MSA
        if self._mafft_combine.isChecked():
            return ALIGN_RUN_MODE_COMBINED
        return ALIGN_RUN_MODE_SEPARATE

    def _set_mafft_run_mode(self, mode: str | None) -> None:
        value = str(mode or ALIGN_RUN_MODE_SEPARATE)
        if value not in ALIGN_RUN_MODES:
            value = ALIGN_RUN_MODE_SEPARATE
        if value == ALIGN_RUN_MODE_ADD_TO_MSA:
            self._mafft_add_to_msa.setChecked(True)
        elif value == ALIGN_RUN_MODE_COMBINED:
            self._mafft_combine.setChecked(True)
        else:
            self._mafft_separate.setChecked(True)

    def _on_mafft_run_mode_changed(self, checked: bool = False) -> None:
        if not checked and self.sender() is not None:
            return
        self._sync_mafft_inputs()
        self._sync_mafft_add_controls()

    def muscle_run_mode(self) -> str:
        if self._muscle_add_to_msa.isChecked():
            return ALIGN_RUN_MODE_ADD_TO_MSA
        if self._muscle_combine.isChecked():
            return ALIGN_RUN_MODE_COMBINED
        return ALIGN_RUN_MODE_SEPARATE

    def _set_muscle_run_mode(self, mode: str | None) -> None:
        value = str(mode or ALIGN_RUN_MODE_SEPARATE)
        if value not in ALIGN_RUN_MODES:
            value = ALIGN_RUN_MODE_SEPARATE
        if value == ALIGN_RUN_MODE_ADD_TO_MSA:
            self._muscle_add_to_msa.setChecked(True)
        elif value == ALIGN_RUN_MODE_COMBINED:
            self._muscle_combine.setChecked(True)
        else:
            self._muscle_separate.setChecked(True)

    def _on_muscle_run_mode_changed(self, checked: bool = False) -> None:
        if not checked and self.sender() is not None:
            return
        self._sync_muscle_inputs()

    def mafft_input_source(self) -> str:
        if self.mafft_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA:
            return ALIGN_INPUT_DEFINE
        if self._mafft_input_define.isChecked():
            return ALIGN_INPUT_DEFINE
        return ALIGN_INPUT_SELECTED

    def mafft_run_files(self, *, panel_files: list[Path] | None = None) -> list[Path]:
        """FASTA inputs for Run: panel selection or expanded Define inputs."""
        from repeatviewer.backends.fastga_dotplot import expand_fastga_input_paths

        if self.mafft_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA:
            return []
        if self.mafft_input_source() == ALIGN_INPUT_DEFINE:
            return expand_fastga_input_paths(self._mafft_define_inputs.paths())
        return [Path(p) for p in (panel_files or [])]

    def mafft_add_paths(self) -> tuple[Path | None, Path | None]:
        """(sequences_to_add, existing_msa) for add-to-MSA mode."""
        seq = self._resolve_existing_align_path(self._mafft_add_sequences.path())
        msa = self._resolve_existing_align_path(self._mafft_add_msa.path())
        return seq, msa

    def _folder_for_pasted_add_sequence(self) -> Path | None:
        """Prefer the Existing MSA folder; otherwise the Sequences path folder."""
        for row in (
            self._mafft_add_msa,
            self._muscle_add_msa,
            self._mafft_add_sequences,
            self._muscle_add_sequences,
        ):
            existing = self._resolve_existing_align_path(row.path())
            if existing is not None:
                return existing.parent
        if self._default_msa is not None:
            return self._default_msa.parent
        return None

    def _paste_sequence_to_add(self, sequences_row: PathBrowseRow) -> None:
        """Open New Sequence, write a FASTA, and fill Sequences to add."""
        from repeatviewer.backends.sequence_compose import (
            SequenceComposeError,
            write_new_sequence,
        )

        folder = self._folder_for_pasted_add_sequence()
        if folder is None:
            QMessageBox.information(
                self,
                "Sequences to add",
                "Choose an Existing MSA first (or set Sequences to add to a "
                "folder), so the pasted sequence has somewhere to be saved.",
            )
            return
        dialog = NewSequenceDialog(
            parent=self,
            title="Sequences to add",
            show_metadata=False,
        )
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        try:
            path = write_new_sequence(folder, **dialog.values())
        except SequenceComposeError as exc:
            QMessageBox.warning(self, "Sequences to add", str(exc))
            return
        sequences_row.set_path(path)

    def _on_mafft_define_paths_changed(self) -> None:
        # Defer so new rows finish layout before height is measured/pinned.
        if not self.isVisible() or self._mafft_advanced_visible():
            return
        self._mafft_compact_size = None
        QTimer.singleShot(0, lambda: self._resize_for_mafft_advanced(False))

    def _sync_mafft_add_controls(self) -> None:
        add_mode = self.mafft_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA
        self._mafft_add_method.setEnabled(add_mode)
        self._mafft_add_keeplength.setEnabled(add_mode)

    def _sync_mafft_inputs(self, *, resize: bool = True) -> None:
        add_mode = self.mafft_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA
        define = (not add_mode) and self.mafft_input_source() == ALIGN_INPUT_DEFINE
        form = self._mafft_basic.layout()
        # Collapse the form row fully — a hidden TeMultiInputFilesWidget still
        # reports a large sizeHint and leaves a blank gap under Settings.
        if isinstance(form, QFormLayout):
            form.setRowVisible(self._mafft_inputs_radios_label, not add_mode)
            form.setRowVisible(self._mafft_define_inputs, define)
            form.setRowVisible(self._mafft_add_sequences_field, add_mode)
            form.setRowVisible(self._mafft_add_msa_field, add_mode)
        self._mafft_input_selected.setVisible(not add_mode)
        self._mafft_input_define.setVisible(not add_mode)
        self._mafft_define_inputs.setVisible(define)
        if define:
            self._mafft_define_inputs.setMinimumHeight(0)
            self._mafft_define_inputs.setMaximumHeight(16777215)
        else:
            self._mafft_define_inputs.setFixedHeight(0)
        for widget in (
            self._mafft_add_sequences,
            self._mafft_add_msa,
            self._mafft_add_sequences_field,
            self._mafft_add_msa_field,
            self._mafft_add_sequences_label,
            self._mafft_add_msa_label,
        ):
            widget.setVisible(add_mode)
        self._mafft_define_inputs.updateGeometry()
        self._mafft_basic.updateGeometry()
        self._sync_mafft_add_controls()
        if hasattr(self, "_muscle_mode_label"):
            self._sync_align_field_columns()
        if (
            resize
            and self.isVisible()
            and hasattr(self, "_advanced_settings_btn")
            and not self._mafft_advanced_visible()
        ):
            self._mafft_compact_size = None
            QTimer.singleShot(0, lambda: self._resize_for_mafft_advanced(False))

    def _apply_mafft_path_from_panel(self, path: Path) -> bool:
        try:
            candidate = Path(path).expanduser().resolve()
        except OSError:
            return False
        if self.mafft_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA:
            return self._fill_add_to_msa_path(
                candidate,
                sequences_row=self._mafft_add_sequences,
                msa_row=self._mafft_add_msa,
            )
        if self.mafft_input_source() != ALIGN_INPUT_DEFINE:
            return False
        return bool(self._mafft_define_inputs.accept_path(candidate))

    def _apply_mafft_folder_from_panel(self, folder: Path) -> bool:
        if self.mafft_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA:
            return False
        if self.mafft_input_source() != ALIGN_INPUT_DEFINE:
            return False
        try:
            candidate = Path(folder).expanduser().resolve()
        except OSError:
            return False
        if not candidate.is_dir():
            return False
        return bool(self._mafft_define_inputs.accept_path(candidate))

    def muscle_input_source(self) -> str:
        if self.muscle_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA:
            return ALIGN_INPUT_DEFINE
        if self._muscle_input_define.isChecked():
            return ALIGN_INPUT_DEFINE
        return ALIGN_INPUT_SELECTED

    def muscle_run_files(self, *, panel_files: list[Path] | None = None) -> list[Path]:
        """FASTA inputs for Run: panel selection or expanded Define inputs."""
        from repeatviewer.backends.fastga_dotplot import expand_fastga_input_paths

        if self.muscle_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA:
            return []
        if self.muscle_input_source() == ALIGN_INPUT_DEFINE:
            return expand_fastga_input_paths(self._muscle_define_inputs.paths())
        return [Path(p) for p in (panel_files or [])]

    def muscle_add_paths(self) -> tuple[Path | None, Path | None]:
        """(sequences_to_add, existing_msa) for add-to-MSA mode."""
        seq = self._resolve_existing_align_path(self._muscle_add_sequences.path())
        msa = self._resolve_existing_align_path(self._muscle_add_msa.path())
        return seq, msa

    def _on_muscle_define_paths_changed(self) -> None:
        # Defer so new rows finish layout before height is measured/pinned.
        if not self.isVisible() or self._muscle_advanced_visible():
            return
        self._muscle_compact_size = None
        QTimer.singleShot(0, lambda: self._resize_for_muscle_advanced(False))

    def _sync_muscle_inputs(self, *, resize: bool = True) -> None:
        add_mode = self.muscle_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA
        define = (not add_mode) and self.muscle_input_source() == ALIGN_INPUT_DEFINE
        form = self._muscle_basic.layout()
        if isinstance(form, QFormLayout):
            form.setRowVisible(self._muscle_inputs_radios_label, not add_mode)
            form.setRowVisible(self._muscle_define_inputs, define)
            form.setRowVisible(self._muscle_add_sequences_field, add_mode)
            form.setRowVisible(self._muscle_add_msa_field, add_mode)
        self._muscle_input_selected.setVisible(not add_mode)
        self._muscle_input_define.setVisible(not add_mode)
        self._muscle_define_inputs.setVisible(define)
        if define:
            self._muscle_define_inputs.setMinimumHeight(0)
            self._muscle_define_inputs.setMaximumHeight(16777215)
        else:
            self._muscle_define_inputs.setFixedHeight(0)
        for widget in (
            self._muscle_add_sequences,
            self._muscle_add_msa,
            self._muscle_add_sequences_field,
            self._muscle_add_msa_field,
            self._muscle_add_sequences_label,
            self._muscle_add_msa_label,
        ):
            widget.setVisible(add_mode)
        self._muscle_define_inputs.updateGeometry()
        self._muscle_basic.updateGeometry()
        if hasattr(self, "_muscle_mode_label"):
            self._sync_align_field_columns()
        if (
            resize
            and self.isVisible()
            and hasattr(self, "_advanced_settings_btn")
            and not self._muscle_advanced_visible()
        ):
            self._muscle_compact_size = None
            QTimer.singleShot(0, lambda: self._resize_for_muscle_advanced(False))

    def _apply_muscle_path_from_panel(self, path: Path) -> bool:
        try:
            candidate = Path(path).expanduser().resolve()
        except OSError:
            return False
        if self.muscle_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA:
            return self._fill_add_to_msa_path(
                candidate,
                sequences_row=self._muscle_add_sequences,
                msa_row=self._muscle_add_msa,
            )
        if self.muscle_input_source() != ALIGN_INPUT_DEFINE:
            return False
        return bool(self._muscle_define_inputs.accept_path(candidate))

    def _apply_muscle_folder_from_panel(self, folder: Path) -> bool:
        if self.muscle_run_mode() == ALIGN_RUN_MODE_ADD_TO_MSA:
            return False
        if self.muscle_input_source() != ALIGN_INPUT_DEFINE:
            return False
        try:
            candidate = Path(folder).expanduser().resolve()
        except OSError:
            return False
        if not candidate.is_dir():
            return False
        return bool(self._muscle_define_inputs.accept_path(candidate))

    @staticmethod
    def _fill_add_to_msa_path(
        candidate: Path,
        *,
        sequences_row: PathBrowseRow,
        msa_row: PathBrowseRow,
    ) -> bool:
        """Fill Sequences if empty, else Existing MSA if empty, else Sequences."""
        if not candidate.is_file():
            return False
        if sequences_row.path() is None:
            target = sequences_row
        elif msa_row.path() is None:
            target = msa_row
        else:
            target = sequences_row
        target.set_path(str(candidate))
        return True

    def fastga_run_files(self) -> list[Path]:
        """Query genomes for Run: panel selection or expanded Define inputs."""
        return self._fastga_panel.files()

    def fastga_options(self):
        return self._fastga_panel.options(
            right_click_action=FASTGA_RIGHT_CLICK_RUN,
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=(
                self._close_process_window.isChecked()
            ),
        )

    def right_click_aligner(self) -> str:
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                return value
        return ALIGNER_MAFFT

    def combined_for(self, aligner: str) -> bool:
        """Whether the given aligner should combine into one MSA."""
        mode = (
            self.muscle_run_mode()
            if aligner == ALIGNER_MUSCLE
            else self.mafft_run_mode()
        )
        return mode == ALIGN_RUN_MODE_COMBINED

    def run_mode_for(self, aligner: str) -> str:
        """Persisted Mode radio for MAFFT or MUSCLE."""
        if aligner == ALIGNER_MUSCLE:
            return self.muscle_run_mode()
        return self.mafft_run_mode()

    def _sync_controls(self) -> None:
        strategy = str(self._strategy.currentData())
        if self._last_strategy == MAFFT_STRATEGY_CUSTOM:
            self._custom_retree = self._retree.value()
            self._custom_max_iterate = self._max_iterate.value()

        if strategy == MAFFT_STRATEGY_AUTO:
            self._set_fixed_spinbox(self._retree, 0, special_text="MAFFT chooses")
            self._set_fixed_spinbox(self._max_iterate, -1, special_text="MAFFT chooses")
        elif strategy == MAFFT_STRATEGY_FFT_NS_1:
            self._set_fixed_spinbox(self._retree, 1)
            self._set_fixed_spinbox(self._max_iterate, 0)
        elif strategy == MAFFT_STRATEGY_FFT_NS_2:
            self._set_fixed_spinbox(self._retree, 2)
            self._set_fixed_spinbox(self._max_iterate, 0)
        elif strategy in {
            MAFFT_STRATEGY_LINSI,
            MAFFT_STRATEGY_GINSI,
            MAFFT_STRATEGY_EINSI,
        }:
            self._set_fixed_spinbox(self._retree, 2)
            self._set_fixed_spinbox(self._max_iterate, 1000)
        else:
            self._retree.setSpecialValueText("")
            self._retree.setRange(1, 100)
            self._retree.setValue(max(1, self._custom_retree))
            self._retree.setEnabled(True)
            self._max_iterate.setSpecialValueText("")
            self._max_iterate.setRange(0, 1000)
            self._max_iterate.setValue(max(0, self._custom_max_iterate))
            self._max_iterate.setEnabled(True)

        self._last_strategy = strategy
        protein = self._sequence_type.currentData() == MAFFT_SEQUENCE_PROTEIN
        self._adjust_direction.setEnabled(not protein)

    @staticmethod
    def _set_fixed_spinbox(
        spinbox: QSpinBox,
        value: int,
        *,
        special_text: str = "",
    ) -> None:
        spinbox.setSpecialValueText(special_text)
        spinbox.setRange(value, value)
        spinbox.setValue(value)
        spinbox.setEnabled(False)

    def options(self) -> AlignOptions:
        right_click = self.right_click_aligner()
        # FastGA is no MSA aligner: the combined / workflow aligner stays
        # MAFFT or MUSCLE even when FastGA is the right-click choice.
        combined = (
            right_click
            if right_click in (ALIGNER_MAFFT, ALIGNER_MUSCLE)
            else self._initial_align.combined_aligner
        )
        mafft_mode = self.mafft_run_mode()
        muscle_mode_run = self.muscle_run_mode()
        add_method = str(self._mafft_add_method.currentData() or MAFFT_ADD_METHOD_ADD)
        return AlignOptions(
            combined_aligner=combined,
            num_threads=max(1, int(self._threads.value())),
            strategy=str(self._strategy.currentData()),
            sequence_type=str(self._sequence_type.currentData()),
            retree=max(1, self._retree.value()),
            max_iterate=max(0, self._max_iterate.value()),
            gap_open_penalty=self._gap_open.value(),
            offset=self._offset.value(),
            reorder=self._reorder.isChecked(),
            adjust_direction=(
                self._adjust_direction.isChecked()
                and self._sequence_type.currentData() != MAFFT_SEQUENCE_PROTEIN
            ),
            muscle_mode=str(self._muscle_mode.currentData()),
            muscle_num_threads=max(1, int(self._muscle_threads.value())),
            muscle_perm=str(self._muscle_perm.currentData()),
            muscle_perturb=max(0, int(self._muscle_perturb.value())),
            muscle_input_source=ALIGN_INPUT_SELECTED,
            muscle_input_paths=[],
            mafft_input_source=ALIGN_INPUT_SELECTED,
            mafft_input_paths=[],
            mafft_run_mode=mafft_mode,
            muscle_run_mode=muscle_mode_run,
            mafft_combined=mafft_mode == ALIGN_RUN_MODE_COMBINED,
            muscle_combined=muscle_mode_run == ALIGN_RUN_MODE_COMBINED,
            mafft_add_method=add_method,
            mafft_add_keeplength=self._mafft_add_keeplength.isChecked(),
            last_dialog_tab=persistable_dialog_tab(
                self.current_dialog_tab(), first_tab=ALIGNER_MAFFT
            ),
            right_click_aligner=self.right_click_aligner(),
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=(
                self._close_process_window.isChecked()
            ),
        )

    def use_system_mafft(self) -> bool:
        return self._use_system_mafft.isChecked()

    def use_system_muscle(self) -> bool:
        return self._use_system_muscle.isChecked()

    def _update_executable_status(self, *_args) -> None:
        """Refresh PATH / FastGA lines on the Options → Executables box."""
        if not hasattr(self, "_mafft_path_status"):
            return
        set_path_checkbox_status(
            self._mafft_path_status,
            resolve_system_executable("mafft"),
            checked=self._use_system_mafft.isChecked(),
            has_builtin=True,
        )
        set_path_checkbox_status(
            self._muscle_path_status,
            resolve_system_executable("muscle"),
            checked=self._use_system_muscle.isChecked(),
            has_builtin=True,
        )
        from repeatviewer.backends.fastga_dotplot import fastga_executable

        set_named_executable_status(
            self._fastga_path_status,
            "FastGA",
            fastga_executable(),
            missing_text="not configured",
        )
        if (
            self.isVisible()
            and self._aligner_tabs.currentWidget() is self._options_tab
        ):
            QTimer.singleShot(0, self._fit_options_content)

    def _restore_defaults(self) -> None:
        """Reset MAFFT/MUSCLE options and use bundled binaries (not system PATH)."""
        defaults = AlignOptions()
        strategy_idx = self._strategy.findData(defaults.strategy)
        if strategy_idx >= 0:
            self._strategy.setCurrentIndex(strategy_idx)
        seq_idx = self._sequence_type.findData(defaults.sequence_type)
        if seq_idx >= 0:
            self._sequence_type.setCurrentIndex(seq_idx)
        self._threads.setValue(max(1, defaults.num_threads))
        self._custom_retree = defaults.retree
        self._custom_max_iterate = defaults.max_iterate
        self._retree.setValue(defaults.retree)
        self._max_iterate.setValue(defaults.max_iterate)
        self._gap_open.setValue(defaults.gap_open_penalty)
        self._offset.setValue(defaults.offset)
        self._reorder.setChecked(defaults.reorder)
        self._adjust_direction.setChecked(defaults.adjust_direction)
        self._set_mafft_run_mode(defaults.mafft_run_mode)
        add_idx = self._mafft_add_method.findData(defaults.mafft_add_method)
        if add_idx >= 0:
            self._mafft_add_method.setCurrentIndex(add_idx)
        self._mafft_add_keeplength.setChecked(defaults.mafft_add_keeplength)
        self._mafft_input_selected.setChecked(True)
        self._mafft_define_inputs.clear_paths()
        self._mafft_add_sequences.set_path(None)
        self._mafft_add_msa.set_path(
            str(self._default_msa) if self._default_msa else None
        )
        self._sync_mafft_inputs()

        muscle_idx = self._muscle_mode.findData(defaults.muscle_mode)
        if muscle_idx >= 0:
            self._muscle_mode.setCurrentIndex(muscle_idx)
        self._muscle_threads.setValue(max(1, defaults.muscle_num_threads))
        muscle_perm_idx = self._muscle_perm.findData(defaults.muscle_perm)
        if muscle_perm_idx >= 0:
            self._muscle_perm.setCurrentIndex(muscle_perm_idx)
        self._muscle_perturb.setValue(max(0, defaults.muscle_perturb))
        self._muscle_input_selected.setChecked(True)
        self._muscle_define_inputs.clear_paths()
        self._muscle_add_sequences.set_path(None)
        self._muscle_add_msa.set_path(
            str(self._default_msa) if self._default_msa else None
        )
        self._set_muscle_run_mode(defaults.muscle_run_mode)
        self._sync_muscle_inputs()

        for value, radio in self._right_click_buttons.items():
            radio.setChecked(value == defaults.right_click_aligner)
        self._show_process_window.setChecked(defaults.show_process_window)
        self._close_process_window.setChecked(
            defaults.close_process_window_when_finished
        )
        # Default: bundled MAFFT / MUSCLE binaries, not system PATH.
        self._use_system_mafft.setChecked(False)
        self._use_system_muscle.setChecked(False)
        self._fastga_panel.restore_defaults()
        self._sync_controls()
        self._update_executable_status()


class NewSequenceDialog(QDialog):
    """Type or paste a sequence and save it as a FASTA.

    ``show_metadata=False`` (Sequences to add): paste area only — name comes
    from a FASTA header when present, otherwise ``added_sequence``.
    """

    # Full New Sequence form.
    _DEFAULT_WIDTH = 720
    _DEFAULT_HEIGHT = 520
    # Sequences to add (paste only) — matches the usual resized paste window.
    _PASTE_DEFAULT_WIDTH = 870
    _PASTE_DEFAULT_HEIGHT = 440

    def __init__(
        self,
        *,
        parent=None,
        title: str = "New Sequence",
        show_metadata: bool = True,
    ) -> None:
        super().__init__(parent)
        from repeatviewer.backends.sequence_compose import (
            MOLECULE_DNA,
            MOLECULE_LABELS,
        )
        from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox

        self._show_metadata = bool(show_metadata)
        self.setWindowTitle(title)
        if self._show_metadata:
            self.resize(self._DEFAULT_WIDTH, self._DEFAULT_HEIGHT)
        else:
            self.setMinimumWidth(520)
            self.setMinimumHeight(280)
            self.resize(self._PASTE_DEFAULT_WIDTH, self._PASTE_DEFAULT_HEIGHT)
            self.setSizeGripEnabled(True)

        self._text = QPlainTextEdit(self)
        self._text.setPlaceholderText("Type or paste in sequence…")
        self._text.setTabChangesFocus(True)

        self._name = QLineEdit(self)
        self._name.setPlaceholderText("Required")
        self._description = QLineEdit(self)
        self._description.setPlaceholderText("Optional")
        self._organism = QLineEdit(self)
        self._organism.setPlaceholderText("Optional")

        self._molecule = QComboBox(self)
        for label, value in MOLECULE_LABELS:
            self._molecule.addItem(label, value)
        dna_idx = self._molecule.findData(MOLECULE_DNA)
        if dna_idx >= 0:
            self._molecule.setCurrentIndex(dna_idx)

        for field in (self._name, self._description, self._organism, self._molecule):
            field.setFixedWidth(320)

        self._summary = QLabel(self)
        self._summary.setStyleSheet("color: palette(mid);")

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(8)
        form.addRow("Name:", self._name)
        form.addRow("Description:", self._description)
        form.addRow("Organism:", self._organism)
        form.addRow("Type:", self._molecule)
        self._form_host = QWidget(self)
        self._form_host.setLayout(form)

        self._buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self._buttons.accepted.connect(self.accept)
        self._buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 14, 18, 14)
        layout.setSpacing(10)
        layout.addWidget(QLabel("Type or paste in sequence:", self))
        layout.addWidget(self._text, stretch=1)
        layout.addWidget(self._form_host)
        layout.addWidget(self._summary)
        layout.addWidget(self._buttons)

        if not self._show_metadata:
            self._form_host.hide()
            self._summary.hide()

        # Track whether the user has typed a name themselves: a pasted FASTA
        # header fills these in, but must never overwrite a typed value.
        self._name_touched = False
        self._description_touched = False
        self._name.textEdited.connect(self._on_name_edited)
        self._description.textEdited.connect(self._on_description_edited)

        self._text.textChanged.connect(self._on_text_changed)
        self._name.textChanged.connect(self._sync)
        self._molecule.currentIndexChanged.connect(self._sync)
        self._sync()

    def _on_name_edited(self, text: str) -> None:
        self._name_touched = bool(text.strip())

    def _on_description_edited(self, text: str) -> None:
        self._description_touched = bool(text.strip())

    def _on_text_changed(self) -> None:
        """Offer the name from a pasted FASTA header."""
        from repeatviewer.backends.sequence_compose import parse_pasted_header

        name, description = parse_pasted_header(self._text.toPlainText())
        if name and not self._name_touched:
            self._name.setText(name)
        if description and not self._description_touched:
            self._description.setText(description)
        self._sync()

    def _resolved_name(self) -> str:
        name = self._name.text().strip()
        if name:
            return name
        if not self._show_metadata:
            from repeatviewer.backends.sequence_compose import parse_pasted_header

            header_name, _ = parse_pasted_header(self._text.toPlainText())
            if header_name:
                return header_name
            return "added_sequence"
        return ""

    def _sync(self) -> None:
        from repeatviewer.backends.sequence_compose import check_sequence

        check = check_sequence(self._text.toPlainText(), str(self._molecule.currentData()))
        if self._show_metadata:
            parts = [f"{check.length:,} bp" if check.length else "empty"]
            if check.invalid:
                parts.append(f"unexpected characters: {check.invalid}")
            self._summary.setText("  ·  ".join(parts))
            ok = bool(check.length) and bool(self._name.text().strip())
        else:
            ok = bool(check.length)
        self._buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(ok)

    def values(self) -> dict:
        return {
            "name": self._resolved_name(),
            "sequence": self._text.toPlainText(),
            "description": (
                self._description.text().strip() if self._show_metadata else ""
            ),
            "organism": self._organism.text().strip() if self._show_metadata else "",
            "molecule": str(self._molecule.currentData()),
        }

    def sizeHint(self) -> QSize:  # noqa: N802
        # Keep resize() effective when dialog polish calls adjustSize().
        hint = super().sizeHint()
        if self._show_metadata:
            return QSize(
                max(hint.width(), self._DEFAULT_WIDTH),
                max(hint.height(), self._DEFAULT_HEIGHT),
            )
        return QSize(
            max(hint.width(), self._PASTE_DEFAULT_WIDTH),
            max(hint.height(), self._PASTE_DEFAULT_HEIGHT),
        )


class NameWithAnnotationDialog(QDialog):
    """Name field plus annotation copy checkbox (combine / separate)."""

    def __init__(
        self,
        *,
        title: str,
        intro: str,
        name_label: str,
        default_name: str,
        annotation_label: str,
        parent=None,
        name_is_folder: bool = False,
        parent_dir: Path | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self._name_is_folder = name_is_folder
        self._parent_dir = (
            parent_dir.expanduser().resolve() if parent_dir is not None else None
        )

        intro_label = QLabel(intro, self)
        intro_label.setWordWrap(True)

        self._name = QLineEdit(default_name, self)
        self._name.selectAll()
        name_px = self.fontMetrics().horizontalAdvance(default_name) + 48
        self.setMinimumWidth(max(520, min(name_px + 160, 900)))

        self._include_annotations = QCheckBox(annotation_label, self)
        self._include_annotations.setChecked(True)

        form = QFormLayout()
        form.addRow(name_label, self._name)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=self,
        )
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 14, 18, 14)
        layout.setSpacing(10)
        layout.addWidget(intro_label)
        layout.addLayout(form)
        layout.addWidget(self._include_annotations)
        layout.addWidget(buttons)

    def name_text(self) -> str:
        return self._name.text().strip()

    def include_annotations(self) -> bool:
        return self._include_annotations.isChecked()

    def resolved_folder(self) -> Path | None:
        """Folder path when ``name_is_folder``; otherwise ``None``."""
        if not self._name_is_folder or self._parent_dir is None:
            return None
        name = self.name_text()
        if not name:
            return None
        return (self._parent_dir / name).expanduser().resolve()

    def _accept(self) -> None:
        name = self.name_text()
        if not name:
            QMessageBox.warning(self, self.windowTitle(), "Enter a name.")
            return
        if Path(name).name != name or "/" in name or "\\" in name:
            QMessageBox.warning(
                self,
                self.windowTitle(),
                "Enter a name only (no folder path).",
            )
            return
        if self._name_is_folder and self._parent_dir is not None:
            target = (self._parent_dir / name).expanduser().resolve()
            if target.exists():
                QMessageBox.warning(
                    self,
                    self.windowTitle(),
                    f"“{name}” folder name already exist.\n\n"
                    "Please choose a different folder name.",
                )
                self._name.selectAll()
                self._name.setFocus()
                return
        self.accept()


class ModifyExtractFileNamesDialog(QDialog):
    """Edit per-sequence output file names when multi-extract paths conflict."""

    def __init__(
        self,
        *,
        record_ids: list[str],
        proposed_names: list[str],
        conflicting: set[str],
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Modify extract file names")
        self.setMinimumWidth(480)
        layout = QVBoxLayout(self)
        layout.addWidget(
            QLabel(
                "Some output file names already exist, or collide after special "
                "characters are sanitized. Edit the file names below "
                "(sequence names are used as file names)."
            )
        )
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        host = QWidget(scroll)
        host_layout = QVBoxLayout(host)
        self._edits: list[QLineEdit] = []
        for rid, name in zip(record_ids, proposed_names, strict=True):
            row = QHBoxLayout()
            label = QLabel(rid, host)
            label.setMinimumWidth(160)
            label.setToolTip(rid)
            row.addWidget(label)
            edit = QLineEdit(name, host)
            if name in conflicting:
                edit.setStyleSheet("QLineEdit { background: #fff3cd; }")
            self._edits.append(edit)
            row.addWidget(edit, stretch=1)
            wrap = QWidget(host)
            wrap.setLayout(row)
            host_layout.addWidget(wrap)
        host_layout.addStretch()
        scroll.setWidget(host)
        scroll.setMaximumHeight(320)
        layout.addWidget(scroll)
        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            parent=self,
        )
        buttons.accepted.connect(self._accept)
        buttons.rejected.connect(self.reject)
        layout.addWidget(buttons)
        self._record_ids = list(record_ids)

    def name_map(self) -> dict[str, str]:
        """Map record id → file name (with extension)."""
        mapping: dict[str, str] = {}
        used: set[str] = set()
        for rid, edit in zip(self._record_ids, self._edits, strict=True):
            name = edit.text().strip()
            if not name:
                raise ValueError("File names cannot be empty.")
            if Path(name).name != name or "/" in name or "\\" in name:
                raise ValueError("Enter a file name only (no folder path).")
            key = name.lower()
            if key in used:
                raise ValueError(f"File name “{name}” is still duplicated.")
            used.add(key)
            mapping[rid] = name
        return mapping

    def _accept(self) -> None:
        try:
            self.name_map()
        except ValueError as exc:
            QMessageBox.warning(self, self.windowTitle(), str(exc))
            return
        self.accept()


class _CheckboxToggleListWidget(QListWidget):
    """List where checkbox clicks toggle one row without range-select side effects.

    ``ExtendedSelection`` + selection↔check sync meant a normal click replaced
    the selection (wiping earlier ticks). Handling the indicator ourselves fixed
    that, but press/move/release could still apply a shift-style range and tick
    the gap (tick 1–3, then 5 → 4 also checked). Steal the full click gesture
    on the checkbox and let the dialog rebuild selection from checks only.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._checkbox_gesture = False

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if event.button() == Qt.MouseButton.LeftButton:
            pos = event.position().toPoint()
            item = self.itemAt(pos)
            if item is not None and self._hit_checkbox(item, pos):
                self._checkbox_gesture = True
                dialog = self._extract_dialog()
                if dialog is not None:
                    dialog._toggle_pick_checkbox(item)
                else:
                    checked = item.checkState() != Qt.CheckState.Checked
                    item.setCheckState(
                        Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
                    )
                event.accept()
                return
        self._checkbox_gesture = False
        super().mousePressEvent(event)

    def mouseMoveEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if self._checkbox_gesture:
            event.accept()
            return
        super().mouseMoveEvent(event)

    def mouseReleaseEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        if self._checkbox_gesture and event.button() == Qt.MouseButton.LeftButton:
            self._checkbox_gesture = False
            event.accept()
            return
        self._checkbox_gesture = False
        super().mouseReleaseEvent(event)

    def _extract_dialog(self) -> ExtractSequencesDialog | None:
        parent = self.parent()
        while parent is not None:
            if isinstance(parent, ExtractSequencesDialog):
                return parent
            parent = parent.parent()
        return None

    def _hit_checkbox(self, item: QListWidgetItem, pos: QPoint) -> bool:
        # Prefer a fixed left strip: stylesheet indicator images often make
        # SE_ItemViewItemCheckIndicator unreliable (empty or oversized).
        row = self.visualItemRect(item)
        if not row.contains(pos):
            return False
        if pos.x() <= row.left() + 24:
            return True
        opt = QStyleOptionViewItem()
        self.initViewItemOption(opt)
        opt.rect = row
        opt.features |= QStyleOptionViewItem.ViewItemFeature.HasCheckIndicator
        opt.checkState = item.checkState()
        if item.isSelected():
            opt.state |= QStyle.StateFlag.State_Selected
        check_rect = self.style().subElementRect(
            QStyle.SubElement.SE_ItemViewItemCheckIndicator, opt, self
        )
        return check_rect.isValid() and check_rect.adjusted(-1, -1, 1, 1).contains(pos)


class ExtractSequencesDialog(QDialog):
    """Pick records out of a multi-FASTA, by hand or from a list of names.

    Two ways in: tick boxes (fine for a genome or a handful of monomers) or an
    ID list (for huge read sets). ``records`` may be ``None`` when there are too
    many sequences for a tick list -- the pick tab then shows a note.
    """

    #: Emitted with the parsed id list when the user presses Check.
    check_requested = Signal(object)

    # Preferred open size. Must be reflected in sizeHint() because dialog
    # polish / adjustSize otherwise shrinks back to content and ignores resize().
    _DEFAULT_WIDTH = 870
    _DEFAULT_HEIGHT = 560

    def __init__(
        self,
        *,
        records: list[tuple[str, str, int]] | None,
        source_name: str,
        initial_ids: str = "",
        initial_column: int = 1,
        records_note: str = "",
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Extract Sequences from File")
        self.setMinimumWidth(640)
        self.setMinimumHeight(420)
        self.resize(self._DEFAULT_WIDTH, self._DEFAULT_HEIGHT)
        self.setSizeGripEnabled(True)
        self._source_name = source_name
        self._did_center = False

        self._tabs = QTabWidget(self)
        self._tabs.addTab(self._build_pick_tab(records, records_note), "Pick from list")
        self._tabs.addTab(self._build_id_tab(initial_ids, initial_column), "From ID list")

        self._single_file = QCheckBox("Write all selected into one file", self)
        self._single_file.setChecked(True)
        self._single_file.setToolTip(
            "On: one multi-FASTA holding every extracted sequence.\n"
            "Off: one file per sequence, which with a long ID list means "
            "one file per matching record."
        )
        self._file_name = QLineEdit(self)
        self._file_name.setPlaceholderText("Output file name")
        self._file_name_form = QFormLayout()
        self._file_name_form.addRow("File name:", self._file_name)
        self._include_annotations = QCheckBox(
            "Copy annotations for extracted sequences", self
        )
        self._include_annotations.setChecked(True)
        self._include_annotations.setToolTip(
            "When the source has a JSON annotation sidecar, write matching "
            "annotations beside the extracted FASTA file(s)."
        )

        self._summary = QLabel(self)
        self._summary.setStyleSheet("color: palette(mid);")

        self._buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        self._buttons.accepted.connect(self.accept)
        self._buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(18, 14, 18, 14)
        layout.setSpacing(10)
        layout.addWidget(self._tabs, stretch=1)
        layout.addWidget(self._single_file)
        layout.addLayout(self._file_name_form)
        layout.addWidget(self._include_annotations)
        layout.addWidget(self._summary)
        layout.addWidget(self._buttons)

        self._tabs.currentChanged.connect(self._sync)
        self._single_file.toggled.connect(self._on_single_file_toggled)
        # Opened with ids already in hand (from a hit table): start on that tab.
        if initial_ids.strip() or records is None:
            self._tabs.setCurrentIndex(1)
        self._sync()
        self._refresh_default_file_name()
        self._on_single_file_toggled(self._single_file.isChecked())

    # -- construction -------------------------------------------------------

    def _build_pick_tab(
        self, records: list[tuple[str, str, int]] | None, records_note: str
    ) -> QWidget:
        page = QWidget(self)
        layout = QVBoxLayout(page)
        layout.setContentsMargins(0, 8, 0, 0)
        layout.setSpacing(8)

        self._list = _CheckboxToggleListWidget(page)
        self._list.setStyleSheet(item_view_checkbox_stylesheet())
        self._pick_syncing = False
        self._filter: QLineEdit | None = None
        # Shift / Cmd (Ctrl) / click-drag rubber-band multi-select.
        self._list.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self._list.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectItems)
        self._list.setUniformItemSizes(True)
        self._list.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._list.setTextElideMode(Qt.TextElideMode.ElideMiddle)
        self._list.setMinimumWidth(0)
        self._list.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
        )

        if records is None:
            note = QLabel(
                records_note
                or "This file is too large to list sequence by sequence. "
                "Use the ID list tab instead.",
                page,
            )
            note.setWordWrap(True)
            note.setStyleSheet("color: palette(mid);")
            layout.addWidget(note)
            layout.addStretch(1)
            self._list.hide()
            return page

        for rec_id, description, length in records:
            label = f"{rec_id}    ({length:,} bp)"
            if description:
                label += f"   {description}"
            item = QListWidgetItem(label, self._list)
            item.setData(Qt.ItemDataRole.UserRole, rec_id)
            item.setFlags(
                item.flags()
                | Qt.ItemFlag.ItemIsUserCheckable
                | Qt.ItemFlag.ItemIsSelectable
                | Qt.ItemFlag.ItemIsEnabled
            )
            item.setCheckState(Qt.CheckState.Unchecked)

        self._filter = QLineEdit(page)
        self._filter.setPlaceholderText("Filter by name")
        self._filter.setClearButtonEnabled(True)
        self._filter.setToolTip(
            "Filter the list by sequence id or description. "
            "Select all applies to currently visible rows."
        )
        self._filter.textChanged.connect(self._apply_pick_filter)

        select_all = QPushButton("Select all", page)
        select_none = QPushButton("Select none", page)
        select_all.setToolTip("Select every visible (filtered) sequence")
        select_none.setToolTip("Clear the selection")
        select_all.clicked.connect(self._select_all_visible)
        select_none.clicked.connect(self._select_none)

        top = QHBoxLayout()
        top.setSpacing(8)
        top.addWidget(self._filter, stretch=1)
        top.addWidget(select_all)
        top.addWidget(select_none)

        layout.addLayout(top)
        layout.addWidget(self._list, stretch=1)
        self._list.itemSelectionChanged.connect(self._on_pick_selection_changed)
        self._list.itemChanged.connect(self._on_pick_item_changed)
        return page
    def _build_id_tab(self, initial_ids: str, initial_column: int) -> QWidget:
        page = QWidget(self)
        layout = QVBoxLayout(page)
        layout.setContentsMargins(0, 8, 0, 0)
        layout.setSpacing(8)

        intro = QLabel(
            "Paste one sequence name per line, or load a list. "
            "Blank lines, '#' comments and a leading '>' are ignored.",
            page,
        )
        intro.setWordWrap(True)
        intro.setStyleSheet("color: palette(mid);")

        self._id_text = QPlainTextEdit(page)
        self._id_text.setPlaceholderText("read_00123\nread_00456\n…")
        self._id_text.setPlainText(initial_ids)

        load = QPushButton("Load from file…", page)
        load.clicked.connect(self._load_id_file)
        load.setToolTip("Read the names from a text file or a hit table.")

        self._column = QSpinBox(page)
        self._column.setRange(1, 50)
        self._column.setValue(max(1, int(initial_column)))
        self._column.setToolTip(
            "Which column holds the names when the pasted text is a table.\n"
            "In a RepeatViewer BLAST table, column 1 is qseqid (the query) "
            "and column 2 is sseqid (the genome sequence hit)."
        )

        self._check_btn = QPushButton("Check", page)
        self._check_btn.clicked.connect(self._request_check)
        self._check_btn.setToolTip(
            "Count how many of these names exist in the file, without writing "
            "anything. Reads the file once."
        )

        self._exclude = QCheckBox("Exclude matching IDs (keep everything else)", page)
        self._exclude.setToolTip(
            "Write every sequence whose name is NOT in the list -- for dropping "
            "contaminant or organelle sequences."
        )

        self._check_result = QLabel(page)
        self._check_result.setWordWrap(True)

        row = QHBoxLayout()
        row.addWidget(load)
        row.addSpacing(12)
        row.addWidget(QLabel("IDs are in column:", page))
        row.addWidget(self._column)
        row.addStretch(1)
        row.addWidget(self._check_btn)

        layout.addWidget(intro)
        layout.addWidget(self._id_text, stretch=1)
        layout.addLayout(row)
        layout.addWidget(self._exclude)
        layout.addWidget(self._check_result)

        self._id_text.textChanged.connect(self._on_ids_changed)
        self._column.valueChanged.connect(self._on_ids_changed)
        self._exclude.toggled.connect(self._sync)
        return page

    # -- pick tab -----------------------------------------------------------

    def _apply_pick_filter(self, text: str = "") -> None:
        needle = str(text or "").strip().lower()
        for index in range(self._list.count()):
            item = self._list.item(index)
            if item is None:
                continue
            if not needle:
                item.setHidden(False)
                continue
            item.setHidden(needle not in item.text().lower())

    def _visible_items(self) -> list[QListWidgetItem]:
        items: list[QListWidgetItem] = []
        for index in range(self._list.count()):
            item = self._list.item(index)
            if item is not None and not item.isHidden():
                items.append(item)
        return items

    def _select_all_visible(self) -> None:
        self._pick_syncing = True
        try:
            self._list.clearSelection()
            for item in self._visible_items():
                item.setSelected(True)
                item.setCheckState(Qt.CheckState.Checked)
            # Clear checks on hidden rows so Select all means "what you see".
            for index in range(self._list.count()):
                item = self._list.item(index)
                if item is not None and item.isHidden():
                    item.setCheckState(Qt.CheckState.Unchecked)
                    item.setSelected(False)
        finally:
            self._pick_syncing = False
        self._sync()

    def _select_none(self) -> None:
        self._pick_syncing = True
        try:
            self._list.clearSelection()
            for index in range(self._list.count()):
                item = self._list.item(index)
                if item is not None:
                    item.setCheckState(Qt.CheckState.Unchecked)
        finally:
            self._pick_syncing = False
        self._sync()

    def _on_pick_selection_changed(self) -> None:
        if self._pick_syncing:
            return
        self._pick_syncing = True
        try:
            for index in range(self._list.count()):
                item = self._list.item(index)
                if item is None or item.isHidden():
                    continue
                want = (
                    Qt.CheckState.Checked
                    if item.isSelected()
                    else Qt.CheckState.Unchecked
                )
                if item.checkState() != want:
                    item.setCheckState(want)
        finally:
            self._pick_syncing = False
        self._sync()

    def _on_pick_item_changed(self, item: QListWidgetItem) -> None:
        """Keep row highlight in sync when a check changes (e.g. keyboard)."""
        if self._pick_syncing or item is None:
            return
        self._pick_syncing = True
        try:
            self._selection_from_checks()
        finally:
            self._pick_syncing = False
        self._sync()

    def _toggle_pick_checkbox(self, item: QListWidgetItem) -> None:
        """Toggle one checkbox and reselect exactly the checked rows (no gaps)."""
        if item is None:
            return
        self._pick_syncing = True
        try:
            checked = item.checkState() != Qt.CheckState.Checked
            item.setCheckState(
                Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
            )
            self._selection_from_checks()
            # Move focus without changing the selection set.
            self._list.setCurrentItem(item, QItemSelectionModel.SelectionFlag.NoUpdate)
        finally:
            self._pick_syncing = False
        self._sync()

    def _selection_from_checks(self) -> None:
        """Make the highlight match the checkboxes exactly (non-contiguous OK)."""
        model = self._list.selectionModel()
        if model is None:
            return
        model.clearSelection()
        for index in range(self._list.count()):
            item = self._list.item(index)
            if item is None or item.checkState() != Qt.CheckState.Checked:
                continue
            model.select(
                self._list.indexFromItem(item),
                QItemSelectionModel.SelectionFlag.Select,
            )

    def selected_ids(self) -> set[str]:
        return {
            self._list.item(i).data(Qt.ItemDataRole.UserRole)
            for i in range(self._list.count())
            if self._list.item(i) is not None
            and self._list.item(i).checkState() == Qt.CheckState.Checked
        }

    # -- id tab -------------------------------------------------------------

    def _load_id_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "Load sequence IDs",
            "",
            "ID lists and tables (*.txt *.tsv *.csv *.ids *.list);;All files (*)",
        )
        if not path:
            return
        try:
            text = Path(path).read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            QMessageBox.warning(self, "Load sequence IDs", str(exc))
            return
        self._id_text.setPlainText(text)

    def _on_ids_changed(self) -> None:
        self._check_result.clear()
        self._sync()

    def parsed_ids(self):
        from repeatviewer.backends.sequence_compose import parse_id_list

        return parse_id_list(self._id_text.toPlainText(), column=self._column.value())

    def _request_check(self) -> None:
        ids = self.parsed_ids().ids
        if not ids:
            self._check_result.setText("No sequence names found in the text above.")
            return
        self.set_check_running()
        self.check_requested.emit(list(ids))

    def set_check_running(self) -> None:
        self._check_btn.setEnabled(False)
        self._check_result.setStyleSheet("color: palette(mid);")
        self._check_result.setText("Checking…")

    def set_check_progress(self, text: str) -> None:
        """Progress while a check runs; leaves the button disabled."""
        self._check_result.setStyleSheet("color: palette(mid);")
        self._check_result.setText(text)

    def set_check_result(self, text: str, *, ok: bool = True) -> None:
        self._check_btn.setEnabled(True)
        self._check_result.setStyleSheet("" if ok else "color: palette(mid);")
        self._check_result.setText(text)

    def id_column(self) -> int:
        return self._column.value()

    def exclude(self) -> bool:
        return self._exclude.isChecked()

    # -- shared -------------------------------------------------------------

    def use_id_list(self) -> bool:
        return self._tabs.currentIndex() == 1

    def sizeHint(self) -> QSize:  # noqa: N802
        hint = super().sizeHint()
        return QSize(
            max(hint.width(), self._DEFAULT_WIDTH),
            max(hint.height(), self._DEFAULT_HEIGHT),
        )

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        if self._did_center:
            return
        self._did_center = True
        # Keep the preferred width even if polish tried to shrink us, then
        # center over the main window (figure-2 layout).
        self.resize(self._DEFAULT_WIDTH, self._DEFAULT_HEIGHT)
        parent = self.parentWidget()
        if parent is not None:
            geo = parent.frameGeometry()
            x = geo.x() + max(0, (geo.width() - self.width()) // 2)
            y = geo.y() + max(0, (geo.height() - self.height()) // 2)
            self.move(x, y)

    def single_file(self) -> bool:
        return self._single_file.isChecked()

    def include_annotations(self) -> bool:
        return self._include_annotations.isChecked()

    def output_file_name(self) -> str:
        """User-chosen name when writing one file; empty when one-per-sequence."""
        if not self.single_file():
            return ""
        return self._file_name.text().strip()

    def _on_single_file_toggled(self, checked: bool) -> None:
        # Show/hide every widget in the form row (label + field).
        for i in range(self._file_name_form.rowCount()):
            label_item = self._file_name_form.itemAt(i, QFormLayout.ItemRole.LabelRole)
            field_item = self._file_name_form.itemAt(i, QFormLayout.ItemRole.FieldRole)
            for item in (label_item, field_item):
                if item is None:
                    continue
                widget = item.widget()
                if widget is not None:
                    widget.setVisible(checked)
        if checked:
            self._refresh_default_file_name()

    def _refresh_default_file_name(self) -> None:
        if not self.single_file():
            return
        from repeatviewer.backends.separate_fasta import sanitize_name
        from repeatviewer.backends.sequence_compose import extract_source_suffix

        if self.use_id_list():
            ids = list(self.parsed_ids().ids)
        else:
            ids = list(self.selected_ids())
        suffix = extract_source_suffix(Path(self._source_name))
        if len(ids) == 1:
            default = f"{sanitize_name(ids[0])}{suffix}"
        else:
            default = f"{self._source_name}_extr{suffix}"
        current = self._file_name.text().strip()
        if not current or current == getattr(self, "_auto_file_name", ""):
            self._file_name.setText(default)
            self._auto_file_name = default

    def _sync(self) -> None:
        if self.use_id_list():
            parsed = self.parsed_ids()
            bits = []
            if parsed.ids:
                verb = "to exclude" if self.exclude() else "to extract"
                bits.append(f"{len(parsed.ids):,} name(s) {verb}")
            else:
                bits.append("No sequence names yet")
            if parsed.duplicates:
                bits.append(f"{parsed.duplicates:,} duplicate(s) collapsed")
            if parsed.ignored:
                bits.append(f"{parsed.ignored:,} line(s) ignored")
            self._summary.setText(" · ".join(bits))
            enabled = bool(parsed.ids)
        else:
            chosen = self.selected_ids()
            visible = sum(
                1
                for i in range(self._list.count())
                if self._list.item(i) is not None and not self._list.item(i).isHidden()
            )
            total = self._list.count()
            if chosen:
                if visible < total:
                    self._summary.setText(
                        f"{len(chosen)} selected · {visible} of {total} shown"
                    )
                else:
                    self._summary.setText(f"{len(chosen)} of {total} selected")
            else:
                self._summary.setText(
                    f"Nothing selected · {visible} of {total} shown"
                    if visible < total
                    else "Nothing selected"
                )
            enabled = bool(chosen)
        self._buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(enabled)
        self._refresh_default_file_name()


_MAFFT_MODE_HELP = (
    "<b>Combine selected sequences and align</b>: merge every sequence from the "
    "selected files into a single FASTA, then run MAFFT once on that combined "
    "file.<br><br>"
    "<b>Align each selected file separately</b>: align each selected FASTA file "
    "separately with MAFFT.<br><br>"
    "<b>Add sequences to an existing MSA</b>: provide unaligned sequences and an "
    "existing MSA; MAFFT inserts the new sequences with "
    "<code>--add</code> / related options (see Advanced Settings)."
)
_MAFFT_ADD_TO_MSA_HELP = (
    "Requires two files: <b>Sequences to add</b> (unaligned FASTA) and "
    "<b>Existing MSA</b> (aligned FASTA). The Files selection prefills the MSA "
    "path when you open this window. Click files in the Files panel to fill the "
    "paths, use Browse, or <b>Paste…</b> to type/paste a sequence (saved next to "
    "the MSA). Add-method flags live under Advanced Settings."
)
_MAFFT_ADD_METHOD_HELP = (
    "Used only in <b>Add sequences to an existing MSA</b> mode.<br><br>"
    "<b>--add</b>: insert new sequences into the existing alignment.<br><br>"
    "<b>--addfragments</b>: better for short / partial sequences.<br><br>"
    "<b>--addfull</b>: add full-length sequences into a fragmentary MSA.<br><br>"
    "<b>--addlong</b>: for sequences much longer than the MSA."
)
_MAFFT_ADD_KEEPLENGTH_HELP = (
    "MAFFT <b>--keeplength</b>: keep the existing MSA length (do not introduce "
    "new columns). Only used with Add sequences to an existing MSA."
)
_MAFFT_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): MAFFT uses the FASTA files that were "
    "selected in the Files panel when you click Run.<br><br>"
    "<b>Define inputs</b>: choose different FASTA files. Paths are for "
    "this run only and are not remembered next time."
)
_MAFFT_STRATEGY_HELP = (
    "MAFFT strategies trade speed for accuracy "
    "(see <a href='https://mafft.cbrc.jp/alignment/software/algorithms/algorithms.html'>"
    "MAFFT algorithms</a>).<br><br>"
    "<b>Auto</b>: MAFFT chooses a strategy.<br><br>"
    "<b>FFT-NS-1</b>: fastest progressive method "
    "(<code>--retree 1</code>).<br><br>"
    "<b>FFT-NS-2</b>: progressive with a rebuilt guide tree "
    "(<code>--retree 2</code>); a good default for large sets.<br><br>"
    "<b>L-INS-i</b>: iterative refinement with local pairwise alignment "
    "(<code>--localpair --maxiterate 1000</code>); one alignable domain with flanks.<br><br>"
    "<b>G-INS-i</b>: iterative refinement with global pairwise alignment "
    "(<code>--globalpair --maxiterate 1000</code>); assumes the whole region aligns.<br><br>"
    "<b>E-INS-i</b>: iterative refinement with generalized affine gaps "
    "(<code>--genafpair --maxiterate 1000</code>); conserved motifs in long "
    "unalignable regions.<br><br>"
    "<b>Custom</b>: set <code>--retree</code> and <code>--maxiterate</code> yourself "
    "in Advanced Settings."
)
_MAFFT_SEQUENCE_TYPE_HELP = (
    "<b>Auto-detect</b>: let MAFFT decide nucleotide vs protein.<br><br>"
    "<b>Nucleotide</b>: force <code>--nuc</code>.<br><br>"
    "<b>Protein</b>: force <code>--amino</code>."
)
_MAFFT_THREADS_HELP = "Number passed to MAFFT's <b>--thread</b> option."
_MAFFT_RETREE_HELP = (
    "MAFFT <b>--retree</b>: how many times to rebuild the guide tree. "
    "Editable only with the Custom strategy; presets fix this value."
)
_MAFFT_MAXITERATE_HELP = (
    "MAFFT <b>--maxiterate</b>: iterative refinement cycles. "
    "Editable only with the Custom strategy; INS-i presets use 1000."
)
_MAFFT_OP_HELP = (
    "MAFFT <b>--op</b>: gap-opening penalty. Default <b>1.53</b> "
    "(tuned for amino-acid alignments)."
)
_MAFFT_EP_HELP = (
    "MAFFT <b>--ep</b>: offset value (behaves like a gap-extension term). "
    "Default in this dialog is <b>0.0</b>."
)
_MAFFT_REORDER_HELP = (
    "MAFFT <b>--reorder</b>: order output sequences by similarity instead of "
    "keeping the input order."
)
_MAFFT_ADJUST_DIRECTION_HELP = (
    "MAFFT <b>--adjustdirectionaccurately</b>: reverse-complement nucleotide "
    "sequences as needed. Disabled for protein."
)

_MUSCLE_MODE_HELP = (
    "<b>Combine selected sequences and align</b>: merge every sequence from the "
    "selected files into a single FASTA, then run MUSCLE once on that combined "
    "file.<br><br>"
    "<b>Align each selected file separately</b>: align each selected FASTA file "
    "separately with MUSCLE.<br><br>"
    "<b>Add sequences to an existing MSA</b>: provide unaligned sequences and an "
    "existing MSA. MUSCLE 5 has no profile-add; it degaps the MSA, merges the "
    "new sequences, and re-aligns. Prefer MAFFT for true add-to-MSA."
)
_MUSCLE_ADD_TO_MSA_HELP = (
    "Requires two files: <b>Sequences to add</b> and <b>Existing MSA</b>. "
    "The Files selection prefills the MSA path when you open this window. "
    "Use Browse or <b>Paste…</b> to type/paste a sequence (saved next to the "
    "MSA). MUSCLE 5 re-aligns the merged sequences (it does not preserve MSA "
    "columns the way MAFFT <code>--add</code> does)."
)
_MUSCLE_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): MUSCLE uses the FASTA files that were "
    "selected in the Files panel when you click Run.<br><br>"
    "<b>Define inputs</b>: choose different FASTA files. Paths are for "
    "this run only and are not remembered next time."
)
_MUSCLE_ALGORITHM_HELP = (
    "<b>Align (high accuracy)</b>: MUSCLE <b>-align</b>, the normal high-accuracy "
    "method.<br><br>"
    "<b>Super5 (large datasets)</b>: MUSCLE <b>-super5</b>, intended for very large "
    "sequence sets."
)
_MUSCLE_THREADS_HELP = "Number passed to MUSCLE's <b>-threads</b> option."
_MUSCLE_PERM_HELP = (
    "MUSCLE <b>-perm</b>: guide-tree permutation for a single replicate alignment "
    "(<b>none</b>, <b>abc</b>, <b>acb</b>, or <b>bca</b>). Default is none."
)
_MUSCLE_PERTURB_HELP = (
    "MUSCLE <b>-perturb</b>: perturbation seed for a single replicate. "
    "<b>0</b> means do not perturb (MUSCLE default)."
)

_FASTGA_MODE_HELP = (
    "<b>Self comparison</b>: each selected FASTA file is compared with itself.<br><br>"
    "<b>Pairwise comparison</b>: each selected file is compared with every other "
    "selected file. Example: with <code>A.fa</code> and <code>B.fa</code>, every "
    "sequence in <code>A.fa</code> is compared with every sequence in "
    "<code>B.fa</code>. All pair results go into one combined file.<br><br>"
    "<b>Compare against a reference</b>: each selected query genome is aligned to one "
    "reference genome. The reference index is built once and reused."
)
_FASTGA_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): FastGA uses the FASTA files that were "
    "selected in the Files panel when you click Run.<br><br>"
    "<b>Define inputs</b>: choose different FASTA files. Paths are for "
    "this run only and are not remembered next time."
)
_FASTGA_REFERENCE_HELP = (
    "The genome every query file is aligned to (the target columns of the PAF).<br><br>"
    "Defaults to the project setting genome (Define Genome). Click another file in "
    "the Files panel while this window is open, or use Browse. If the reference is "
    "also among the query inputs it is skipped, so it is never compared with itself."
)
_FASTGA_COMBINED_HELP = (
    "Write the alignments of all selected files into one PAF/PSL file, "
    "&lt;first&gt;_multi_FastGA_vs_&lt;reference&gt;.paf, instead of one file per "
    "selected file.<br><br>"
    "Unavailable for <b>FastGA alignment (.1aln)</b> — that binary format cannot "
    "be concatenated, so each query gets its own .1aln file."
)
_FASTGA_FORMAT_HELP = (
    "<b>PAF</b> is the usual text alignment format. Optional <b>Include CIGAR</b> or "
    "<b>Include CS tag</b> add extra columns (only one at a time). Pairwise mode "
    "writes one combined <code>*_multi_FastGA.paf</code>.<br><br>"
    "<b>PSL</b> is UCSC PSL text; it does not open in the dotplot view.<br><br>"
    "<b>FastGA alignment (.1aln)</b> is FastGA's binary ONEcode output; it opens in "
    "the dotplot view. Pairwise and reference runs write <b>one .1aln per "
    "comparison</b> (binary files cannot be joined)."
)
_FASTGA_THREADS_HELP = "FastGA <b>-T</b>: number of threads."
_FASTGA_IDENTITY_HELP = (
    "FastGA <b>-i</b>: minimum alignment identity (0.70 = 70 %). Alignments below this "
    "are not reported."
)
_FASTGA_LENGTH_HELP = (
    "FastGA <b>-l</b>: shortest alignment to report, in bp. It cannot go below 100, "
    "because FastGA rejects limits shorter than its seed length.<br>"
    "Raise it for satellite arrays, where short matches can swamp the output."
)
_FASTGA_SEED_HELP = (
    "FastGA <b>-f</b>: seeds that occur more often than this are treated as repetitive "
    "and dropped before alignment.<br>"
    "FastGA's default is 10; RepeatViewer Align uses the same default. Higher values "
    "are slower and keep more repetitive seeds."
)
_FASTGA_CHAIN_COV_HELP = (
    "FastGA <b>-c</b>: a seed chain must cover at least this many bases in both genomes "
    "to become an alignment. FastGA default: 85."
)
_FASTGA_CHAIN_GAP_HELP = (
    "FastGA <b>-s</b>: largest spacing between two consecutive seeds in one chain. "
    "FastGA default: 1000."
)
_FASTGA_SYMMETRIC_HELP = (
    "FastGA <b>-S</b>: take seeds from both genomes, so the result does not depend on "
    "which one comes first. Finds more repeat copies in the second genome, but is "
    "slower. FastGA recommends it for repeat structure, not for synteny."
)
_FASTGA_SOFT_MASK_HELP = (
    "Alignments may only <b>start</b> in uppercase bases; lowercase (soft-masked) "
    "regions cannot start a hit. Default off. With soft masking on, a reference index "
    "is rebuilt every run instead of reused."
)
_FASTGA_KEEP_INDEX_HELP = (
    "FastGA <b>-k</b>: keep the genome index (.1gdb / .gix) so later runs skip building "
    "it. An index takes about 14 GB per gigabase of genome.<br><br>"
    "<b>Against itself</b> and <b>against each other</b>: kept next to each FASTA.<br><br>"
    "<b>Against a reference</b>: only the reference index is kept, in a hidden "
    ".repeatviewer_fastga folder next to it, and reused while the reference is unchanged. "
    "Indexes of the other files are always removed."
)


class FastGaPanel(QWidget):
    """The FastGA tab of the Alignment window.

    Holds the FastGA-only settings (mode, reference, output, alignment
    parameters, and index). Right-click action and the process window live
    in the Options tab of the owning dialog.
    """

    changed = Signal()
    # Emitted when Mode / Inputs rows show or hide (dialog should refit height).
    layout_changed = Signal()

    def __init__(
        self,
        *,
        initial=None,
        files=None,
        default_reference: Path | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        from repeatviewer.job_options import FastGaOptions

        opts = initial or FastGaOptions()
        self._initial = opts
        self._files = [Path(p) for p in (files or [])]
        self._default_reference = self._resolve_existing_path(default_reference)

        # --- Mode ---------------------------------------------------------
        self._mode_group = QButtonGroup(self)
        self._mode_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("Self comparison", FASTGA_MODE_SELF),
            ("Pairwise comparison", FASTGA_MODE_PAIRS),
            ("Compare against a reference", FASTGA_MODE_REFERENCE),
        ):
            radio = QRadioButton(label, self)
            self._mode_group.addButton(radio)
            self._mode_buttons[value] = radio
        self._set_mode(opts.mode)

        reference_initial = self._resolve_existing_path(opts.reference_fasta)
        if reference_initial is None:
            reference_initial = self._default_reference
        self._reference = PathBrowseRow(
            self,
            initial=str(reference_initial) if reference_initial else None,
            placeholder="Click in Files or Browse…",
            browse_title="Choose the reference FASTA",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fa.gz *.fasta.gz *.fna.gz);;"
                "All files (*)"
            ),
            minimum_width=280,
        )
        # Expand with the dialog (same fill as Define-inputs path rows).
        self._reference.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._reference.line_edit().setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._combined = QCheckBox("Write all into one file", self)
        self._combined.setChecked(bool(opts.reference_combined))

        # --- Output -------------------------------------------------------
        self._format = QComboBox(self)
        for label, value in (
            ("PAF", FASTGA_FORMAT_PAF),
            ("PSL", FASTGA_FORMAT_PSL),
            ("FastGA alignment (.1aln)", FASTGA_FORMAT_ALN),
        ):
            self._format.addItem(label, value)
        self._format.setFixedWidth(280)
        self._cigar = QCheckBox("Include CIGAR", self)
        self._cs = QCheckBox("Include CS tag", self)
        self._apply_saved_output_format(opts.output_format)

        # --- Settings (threads + advanced) --------------------------------
        self._threads = QSpinBox(self)
        self._threads.setRange(1, 128)
        self._threads.setValue(max(1, opts.threads))
        self._threads.setFixedWidth(280)

        self._identity = QDoubleSpinBox(self)
        self._identity.setRange(0.5, 1.0)
        self._identity.setSingleStep(0.05)
        self._identity.setDecimals(2)
        self._identity.setValue(opts.min_identity)
        self._min_length = QSpinBox(self)
        self._min_length.setRange(100, 1_000_000)
        self._min_length.setSingleStep(50)
        self._min_length.setValue(max(100, opts.min_length))
        self._seed_freq = QSpinBox(self)
        self._seed_freq.setRange(1, 100_000)
        self._seed_freq.setSingleStep(50)
        self._seed_freq.setValue(max(1, opts.seed_freq))
        self._chain_coverage = QSpinBox(self)
        self._chain_coverage.setRange(1, 1_000_000)
        self._chain_coverage.setValue(max(1, opts.chain_coverage))
        self._chain_gap = QSpinBox(self)
        self._chain_gap.setRange(1, 10_000_000)
        self._chain_gap.setSingleStep(100)
        self._chain_gap.setValue(max(1, opts.chain_gap))
        self._symmetric = QCheckBox("Symmetric alignment (-S)", self)
        self._symmetric.setChecked(opts.symmetric)
        self._soft_mask = QCheckBox("Soft masking (uppercase seeds only)", self)
        self._soft_mask.setChecked(bool(opts.soft_mask))
        for field in (
            self._identity, self._min_length, self._seed_freq,
            self._chain_coverage, self._chain_gap,
        ):
            field.setFixedWidth(280)

        self._keep_index = QCheckBox("Keep genome index (-k)", self)
        self._keep_index.setChecked(opts.keep_index)

        # --- Inputs -------------------------------------------------------
        self._input_group = QButtonGroup(self)
        self._input_selected = QRadioButton("Use selected sequences", self)
        self._input_define = QRadioButton("Define inputs", self)
        self._input_group.addButton(self._input_selected)
        self._input_group.addButton(self._input_define)
        # Always open on Use selected (Define is for this run only).
        self._input_selected.setChecked(True)
        # Define-input paths are session-only; do not restore last run's paths.
        self._define_inputs = TeMultiInputFilesWidget(
            self,
            initial=None,
            start_empty=False,
            placeholder="Click in Files/Folders or Browse…",
            browse_title="Choose a FASTA genome",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
            ),
            add_button_label="Add another input",
            add_button_tooltip=(
                "Add another FASTA file path. Click a file in the Files panel "
                "to fill a row, or use Browse."
            ),
            path_mode="file",
        )
        self._path_target = "input"  # or "reference"

        # --- Layout -------------------------------------------------------
        def _section(title: str, parent: QWidget | None = None) -> QLabel:
            return QLabel(f"<b>{title}</b>", parent or self)

        root = QVBoxLayout(self)
        root.setContentsMargins(16, 10, 16, 10)
        root.setSpacing(4)
        root.setAlignment(Qt.AlignmentFlag.AlignTop)

        basic = QWidget(self)
        self._basic_panel = basic
        basic.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )
        form = QFormLayout(basic)
        form.setContentsMargins(0, 0, 0, 0)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        # Match MAFFT/MUSCLE: Define-inputs path + Browse fill the row; fixed
        # Settings/Output fields (280) stay put.
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)

        form.addRow(_section("Mode", basic))
        mode_box = QWidget(basic)
        mode_layout = QVBoxLayout(mode_box)
        mode_layout.setContentsMargins(0, 0, 0, 0)
        mode_layout.setSpacing(2)
        for value, radio in self._mode_buttons.items():
            if value == FASTGA_MODE_SELF:
                mode_layout.addWidget(
                    field_with_help(
                        radio, "FastGA mode", _FASTGA_MODE_HELP, stretch=True
                    )
                )
            elif value == FASTGA_MODE_REFERENCE:
                mode_layout.addWidget(
                    field_with_help(
                        radio,
                        "Compare against a reference",
                        _FASTGA_REFERENCE_HELP,
                        stretch=True,
                    )
                )
            else:
                mode_layout.addWidget(radio)
        self._mode_label = QLabel("")
        form.addRow(self._mode_label, mode_box)
        self._reference_label = QLabel("Reference genome:", basic)
        # Trailing help-sized spacer matches Define-inputs rows (path/Browse align).
        self._reference_field = QWidget(basic)
        self._reference_field.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        ref_row = QHBoxLayout(self._reference_field)
        ref_row.setContentsMargins(0, 0, 0, 0)
        ref_row.setSpacing(6)
        ref_row.addWidget(self._reference, stretch=1)
        ref_help_spacer = QWidget(self._reference_field)
        ref_help_spacer.setFixedSize(HELP_BUTTON_SIZE, HELP_BUTTON_SIZE)
        ref_row.addWidget(ref_help_spacer)
        form.addRow(self._reference_label, self._reference_field)
        self._combined_field = field_with_help(
            self._combined, "Write all into one file", _FASTGA_COMBINED_HELP
        )
        self._combined_spacer = QLabel("")
        form.addRow(self._combined_spacer, self._combined_field)
        self._combined_label = form.labelForField(self._combined_field)

        form.addRow(_form_section_separator())
        form.addRow(_section("Inputs", basic))
        inputs_box = QWidget(basic)
        inputs_layout = QVBoxLayout(inputs_box)
        inputs_layout.setContentsMargins(0, 0, 0, 0)
        inputs_layout.setSpacing(2)
        inputs_layout.addWidget(
            field_with_help(
                self._input_selected, "FastGA inputs", _FASTGA_INPUTS_HELP, stretch=True
            )
        )
        inputs_layout.addWidget(self._input_define)
        self._inputs_radios_label = QLabel("")
        self._inputs_define_label = QLabel("")
        form.addRow(self._inputs_radios_label, inputs_box)
        form.addRow(self._inputs_define_label, self._define_inputs)

        form.addRow(_form_section_separator())
        form.addRow(_section("Output", basic))
        form.addRow(
            "Format:",
            field_with_help(self._format, "Output format", _FASTGA_FORMAT_HELP),
        )
        self._cigar_spacer = QLabel("")
        self._cs_spacer = QLabel("")
        form.addRow(self._cigar_spacer, self._cigar)
        form.addRow(self._cs_spacer, self._cs)

        form.addRow(_form_section_separator())
        form.addRow(_section("Settings", basic))
        form.addRow(
            "Threads:",
            field_with_help(self._threads, "Threads", _FASTGA_THREADS_HELP),
        )
        root.addWidget(basic)

        # Sibling of the basic form (not a QFormLayout row): hide/show is reliable.
        self._advanced_panel = QWidget(self)
        self._advanced_panel.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        advanced_form = QFormLayout(self._advanced_panel)
        advanced_form.setContentsMargins(0, 4, 0, 0)
        advanced_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        advanced_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        advanced_form.setHorizontalSpacing(12)
        advanced_form.setVerticalSpacing(4)
        advanced_form.addRow(_form_section_separator())
        self._advanced_header = _section("Advanced Settings", self._advanced_panel)
        advanced_form.addRow(self._advanced_header)
        advanced_form.addRow(
            "Minimum identity:",
            field_with_help(self._identity, "Minimum identity", _FASTGA_IDENTITY_HELP),
        )
        advanced_form.addRow(
            "Minimum alignment length:",
            field_with_help(self._min_length, "Minimum alignment length", _FASTGA_LENGTH_HELP),
        )
        advanced_form.addRow(
            "Seed cutoff:",
            field_with_help(self._seed_freq, "Seed cutoff", _FASTGA_SEED_HELP),
        )
        advanced_form.addRow(
            "Chain coverage:",
            field_with_help(self._chain_coverage, "Chain coverage", _FASTGA_CHAIN_COV_HELP),
        )
        advanced_form.addRow(
            "Chain gap:",
            field_with_help(self._chain_gap, "Chain gap", _FASTGA_CHAIN_GAP_HELP),
        )
        self._soft_mask_spacer = QLabel("")
        self._symmetric_spacer = QLabel("")
        self._keep_index_spacer = QLabel("")
        advanced_form.addRow(
            self._soft_mask_spacer,
            field_with_help(self._soft_mask, "Soft masking", _FASTGA_SOFT_MASK_HELP),
        )
        advanced_form.addRow(
            self._symmetric_spacer,
            field_with_help(self._symmetric, "Symmetric alignment", _FASTGA_SYMMETRIC_HELP),
        )
        advanced_form.addRow(
            self._keep_index_spacer,
            field_with_help(self._keep_index, "Keep genome index", _FASTGA_KEEP_INDEX_HELP),
        )
        self._basic_form = form
        self._advanced_form = advanced_form
        root.addWidget(self._advanced_panel)
        self.set_advanced_visible(False)
        self._sync_field_column()

        for radio in self._mode_buttons.values():
            radio.toggled.connect(self._sync)
        self._input_selected.toggled.connect(self._sync)
        self._input_define.toggled.connect(self._sync)
        self._define_inputs.paths_changed.connect(self._on_define_paths_changed)
        self._define_inputs.field_activated.connect(
            lambda: setattr(self, "_path_target", "input")
        )
        self._reference.activated.connect(
            lambda _w: setattr(self, "_path_target", "reference")
        )
        self._reference.text_changed.connect(lambda _t: self._sync())
        self._combined.toggled.connect(self._sync)
        self._format.currentIndexChanged.connect(self._on_format_changed)
        self._cigar.toggled.connect(self._on_cigar_toggled)
        self._cs.toggled.connect(self._on_cs_toggled)
        for w in (
            self._threads, self._identity, self._min_length, self._seed_freq,
            self._chain_coverage, self._chain_gap, self._symmetric, self._soft_mask,
            self._keep_index,
        ):
            if hasattr(w, "valueChanged"):
                w.valueChanged.connect(lambda _v: self.changed.emit())
            elif hasattr(w, "toggled"):
                w.toggled.connect(lambda _c: self.changed.emit())
        self._sync()

    # ----------------------------------------------------------------- state

    def set_advanced_visible(self, visible: bool) -> None:
        self._advanced_panel.setVisible(visible)
        # Force the layout to drop or reclaim height (setVisible alone can leave a gap).
        if visible:
            self._advanced_panel.setMaximumHeight(16777215)
        else:
            self._advanced_panel.setMaximumHeight(0)
        self._advanced_panel.updateGeometry()
        root = self.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        self.updateGeometry()

    def advanced_visible(self) -> bool:
        return bool(self._advanced_panel.isVisible())

    def compact_content_height(self) -> int:
        """Height of the FastGA form without Advanced Settings."""
        root = self.layout()
        margins = root.contentsMargins() if root is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 20
        spacing = root.spacing() if root is not None else 4
        form = getattr(self, "_basic_form", None)
        if isinstance(form, QFormLayout):
            form.activate()
            basic_h = max(form.totalMinimumSize().height(), form.totalSizeHint().height())
        else:
            basic_h = max(
                self._basic_panel.minimumSizeHint().height(),
                self._basic_panel.sizeHint().height(),
            )
        # Root keeps one spacing slot between basic and the collapsed Advanced panel.
        return margin_h + basic_h + spacing

    def _set_form_row_expanded(self, widget: QWidget, expanded: bool) -> None:
        """Show or fully collapse a QFormLayout field row (no leftover blank)."""
        form = getattr(self, "_basic_form", None)
        if isinstance(form, QFormLayout):
            form.setRowVisible(widget, expanded)
        widget.setVisible(expanded)
        if expanded:
            widget.setMinimumHeight(0)
            widget.setMaximumHeight(16777215)
        else:
            widget.setFixedHeight(0)
        widget.updateGeometry()

    def _apply_saved_output_format(self, fmt: str) -> None:
        fmt = str(fmt)
        self._cigar.blockSignals(True)
        self._cs.blockSignals(True)
        self._cigar.setChecked(False)
        self._cs.setChecked(False)
        if fmt in (FASTGA_FORMAT_PAFX, FASTGA_FORMAT_PAFM):
            idx = self._format.findData(FASTGA_FORMAT_PAF)
            self._format.setCurrentIndex(idx if idx >= 0 else 0)
            self._cigar.setChecked(True)
        elif fmt in (FASTGA_FORMAT_PAFS, FASTGA_FORMAT_PAFS_LONG):
            idx = self._format.findData(FASTGA_FORMAT_PAF)
            self._format.setCurrentIndex(idx if idx >= 0 else 0)
            self._cs.setChecked(True)
        elif fmt == FASTGA_FORMAT_PSL:
            idx = self._format.findData(FASTGA_FORMAT_PSL)
            self._format.setCurrentIndex(idx if idx >= 0 else 0)
        elif fmt == FASTGA_FORMAT_ALN:
            idx = self._format.findData(FASTGA_FORMAT_ALN)
            self._format.setCurrentIndex(idx if idx >= 0 else 0)
        else:
            idx = self._format.findData(FASTGA_FORMAT_PAF)
            self._format.setCurrentIndex(idx if idx >= 0 else 0)
        self._cigar.blockSignals(False)
        self._cs.blockSignals(False)
        self._update_paf_extras_enabled()

    def _output_format_value(self) -> str:
        fmt = str(self._format.currentData())
        if fmt == FASTGA_FORMAT_PAF:
            if self._cigar.isChecked():
                return FASTGA_FORMAT_PAFX
            if self._cs.isChecked():
                return FASTGA_FORMAT_PAFS
            return FASTGA_FORMAT_PAF
        return fmt

    def _update_paf_extras_enabled(self) -> None:
        paf = str(self._format.currentData()) == FASTGA_FORMAT_PAF
        self._cigar.setEnabled(paf)
        self._cs.setEnabled(paf)
        if not paf:
            self._cigar.blockSignals(True)
            self._cs.blockSignals(True)
            self._cigar.setChecked(False)
            self._cs.setChecked(False)
            self._cigar.blockSignals(False)
            self._cs.blockSignals(False)
        # .1aln cannot be concatenated — disable "Write all into one file".
        aln = str(self._format.currentData()) == FASTGA_FORMAT_ALN
        self._combined.setEnabled(not aln)
        if aln and self._combined.isChecked():
            self._combined.blockSignals(True)
            self._combined.setChecked(False)
            self._combined.blockSignals(False)

    def _on_format_changed(self, _index: int) -> None:
        self._update_paf_extras_enabled()
        self._sync()

    def _on_cigar_toggled(self, checked: bool) -> None:
        if checked:
            self._cs.blockSignals(True)
            self._cs.setChecked(False)
            self._cs.blockSignals(False)
        self._sync()

    def _on_cs_toggled(self, checked: bool) -> None:
        if checked:
            self._cigar.blockSignals(True)
            self._cigar.setChecked(False)
            self._cigar.blockSignals(False)
        self._sync()

    def _set_mode(self, mode: str) -> None:
        radio = self._mode_buttons.get(mode) or self._mode_buttons[FASTGA_MODE_SELF]
        radio.setChecked(True)

    def mode(self) -> str:
        for value, radio in self._mode_buttons.items():
            if radio.isChecked():
                return value
        return FASTGA_MODE_SELF

    def reference(self) -> Path | None:
        return self._reference.path()

    def input_source(self) -> str:
        if self._input_define.isChecked():
            return FASTGA_INPUT_DEFINE
        return FASTGA_INPUT_SELECTED

    def defined_input_paths(self) -> list[Path]:
        return self._define_inputs.paths()

    def files(self) -> list[Path]:
        """Query genomes FastGA will run on (selection or expanded Define inputs)."""
        if self.input_source() == FASTGA_INPUT_DEFINE:
            from repeatviewer.backends.fastga_dotplot import expand_fastga_input_paths

            return expand_fastga_input_paths(self.defined_input_paths())
        return list(self._files)

    def _on_define_paths_changed(self) -> None:
        """Row add/remove must grow the dialog — mode toggle alone is not enough."""
        self._sync()
        if self.input_source() == FASTGA_INPUT_DEFINE:
            # Always emit: _sync only fires layout_changed on Define on/off.
            self.layout_changed.emit()

    def _sync(self) -> None:
        mode = self.mode()
        reference_mode = mode == FASTGA_MODE_REFERENCE
        define = self.input_source() == FASTGA_INPUT_DEFINE
        prev_ref = bool(getattr(self, "_layout_reference_mode", False))
        prev_define = bool(getattr(self, "_layout_define_mode", False))
        # Collapse reference / Define rows fully — setVisible alone leaves blank
        # gaps in the Mode section and an oversized dialog height.
        self._set_form_row_expanded(self._reference_field, reference_mode)
        self._set_form_row_expanded(self._combined_field, reference_mode)
        for widget in (
            self._reference_label,
            self._combined_label,
            self._combined_spacer,
        ):
            if widget is not None:
                widget.setVisible(reference_mode)
        self._set_form_row_expanded(self._define_inputs, define)
        self._inputs_define_label.setVisible(define)

        if hasattr(self, "_basic_panel"):
            self._basic_panel.updateGeometry()
        root = self.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        self.updateGeometry()
        self._sync_field_column()

        self.changed.emit()
        self._layout_reference_mode = reference_mode
        self._layout_define_mode = define
        if reference_mode != prev_ref or define != prev_define:
            self.layout_changed.emit()

    def _sync_field_column(self) -> None:
        if not hasattr(self, "_inputs_radios_label"):
            return
        sync_form_field_column_spacers(
            self,
            [
                self._mode_label,
                self._reference_label,
                self._combined_spacer,
                self._inputs_radios_label,
                self._inputs_define_label,
                self._cigar_spacer,
                self._cs_spacer,
                self._soft_mask_spacer,
                self._symmetric_spacer,
                self._keep_index_spacer,
            ],
            [
                getattr(self, "_basic_form", None),
                getattr(self, "_advanced_form", None),
            ],
        )

    def apply_path_from_panel(self, path: Path) -> bool:
        """Fill Define-inputs or the reference from a Files-panel click."""
        from repeatviewer.backends.fastga_dotplot import is_fastga_fasta_path

        candidate = Path(path).expanduser()
        try:
            candidate = candidate.resolve()
        except OSError:
            return False
        if not candidate.is_file() or not is_fastga_fasta_path(candidate):
            return False

        define = self.input_source() == FASTGA_INPUT_DEFINE
        prefer_input = define and self._path_target == "input"
        prefer_reference = (
            self.mode() == FASTGA_MODE_REFERENCE and self._path_target == "reference"
        )
        if prefer_input or (define and not prefer_reference):
            if self._define_inputs.accept_path(candidate):
                self._path_target = "input"
                return True
        if self.mode() == FASTGA_MODE_REFERENCE:
            self._reference.set_path(candidate)
            self._path_target = "reference"
            return True
        if define and self._define_inputs.accept_path(candidate):
            self._path_target = "input"
            return True
        return False

    def apply_folder_from_panel(self, folder: Path) -> bool:
        """Fill a Define-inputs row from a Folders-panel click."""
        if self.input_source() != FASTGA_INPUT_DEFINE:
            return False
        candidate = Path(folder).expanduser()
        try:
            candidate = candidate.resolve()
        except OSError:
            return False
        if not candidate.is_dir():
            return False
        if self._define_inputs.accept_path(candidate):
            self._path_target = "input"
            return True
        return False

    @staticmethod
    def _resolve_existing_path(path) -> Path | None:
        if path is None or path == "":
            return None
        candidate = Path(path).expanduser()
        try:
            candidate = candidate.resolve()
        except OSError:
            return None
        return candidate if candidate.is_file() else None

    def restore_defaults(self) -> None:
        from repeatviewer.job_options import FastGaOptions

        d = FastGaOptions()
        self._set_mode(d.mode)
        self._input_selected.setChecked(True)
        self._define_inputs.clear_paths()
        self._combined.setChecked(d.reference_combined)
        self._apply_saved_output_format(d.output_format)
        self._threads.setValue(max(1, d.threads))
        self._identity.setValue(d.min_identity)
        self._min_length.setValue(max(100, d.min_length))
        self._seed_freq.setValue(d.seed_freq)
        self._chain_coverage.setValue(d.chain_coverage)
        self._chain_gap.setValue(d.chain_gap)
        self._symmetric.setChecked(d.symmetric)
        self._soft_mask.setChecked(d.soft_mask)
        self._keep_index.setChecked(d.keep_index)
        # Reference defaults to the project setting genome when available.
        if self._default_reference is not None:
            self._reference.set_path(self._default_reference)
        self._sync()

    def options(self, *, right_click_action: str, show_process_window: bool,
                close_process_window_when_finished: bool = True):
        from dataclasses import replace

        ref = self._reference.path()
        return replace(
            self._initial,
            output_format=self._output_format_value(),
            threads=int(self._threads.value()),
            min_identity=float(self._identity.value()),
            min_length=int(self._min_length.value()),
            seed_freq=int(self._seed_freq.value()),
            chain_coverage=int(self._chain_coverage.value()),
            chain_gap=int(self._chain_gap.value()),
            symmetric=self._symmetric.isChecked(),
            keep_index=self._keep_index.isChecked(),
            soft_mask=self._soft_mask.isChecked(),
            mode=self.mode(),
            reference_fasta=str(ref) if ref else "",
            reference_combined=self._combined.isChecked(),
            temp_dir="",
            input_source=self.input_source(),
            input_paths=[str(p) for p in self.defined_input_paths()],
            selected_pairs=None,
            right_click_action=right_click_action,
            show_process_window=bool(show_process_window),
            close_process_window_when_finished=bool(
                close_process_window_when_finished
            ),
        )


def reference_genome_combo_label(
    project: RepeatViewerProject | None,
    genome: Path | None,
) -> str:
    """Display name for the project reference genome in genome/database combos.

    Prefer the Manage Database display name when the FASTA is registered there;
    otherwise use the file name. Never wraps the label as ``Current genome (…)``.
    """
    if genome is None:
        return "Reference genome (not configured)"
    try:
        key = str(genome.expanduser().resolve())
    except OSError:
        key = str(genome)
    if project is not None:
        for managed in getattr(project, "managed_genomes", ()) or ():
            try:
                managed_path = Path(managed.path).expanduser().resolve()
            except OSError:
                continue
            if str(managed_path) == key and getattr(managed, "name", ""):
                return str(managed.name)
        for db in getattr(project, "blast_databases", ()) or ():
            if getattr(db, "molecule_type", None) != BLAST_MOLECULE_NUCL:
                continue
            try:
                source = Path(db.source_fasta).expanduser().resolve()
            except OSError:
                continue
            if str(source) == key and getattr(db, "name", ""):
                return str(db.name)
    return genome.name


class GenomeBlastDatabaseRow(QWidget):
    """Genome combo matching Blast Database: reference genome + nucleotide DBs + Manage."""

    selection_changed = Signal()

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        project: RepeatViewerProject | None = None,
        project_genome: Path | None = None,
        selected_path: str | None = None,
        manage_title: str = "Manage Database",
        help_title: str | None = None,
        help_text: str | None = None,
        show_manage: bool = True,
    ) -> None:
        super().__init__(parent)
        self._project = project
        self._project_genome = (
            project_genome.expanduser().resolve()
            if project_genome is not None and project_genome.is_file()
            else None
        )
        self._manage_title = manage_title
        self._blast_manager_dialog = None

        self._combo = QComboBox(self)
        self._combo.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        # Fixed height: TE-Aid pangenome rows must not squash when more genomes
        # are added and the dialog briefly measures under a tight max height.
        self._combo.ensurePolished()
        self._combo.setFixedHeight(max(22, self._combo.sizeHint().height()))
        self._combo.setSizePolicy(
            QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed
        )
        self._combo.currentIndexChanged.connect(lambda *_: self.selection_changed.emit())
        self.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        # Help sits on the combo row only (Align/BLAST), not centered on Manage.
        if help_title and help_text:
            layout.addWidget(field_with_help(self._combo, help_title, help_text))
        else:
            layout.addWidget(self._combo)
        if show_manage:
            manage = QPushButton("Manage Database", self)
            manage.clicked.connect(self._manage_databases)
            layout.addWidget(manage, 0, Qt.AlignmentFlag.AlignLeft)

        self.refresh(selected_path=selected_path)

    def selected_path(self) -> str | None:
        """Return the selected FASTA path, or None if reference genome is not configured."""
        data = self._combo.currentData()
        if data:
            return str(data)
        if self._project_genome is not None:
            return str(self._project_genome)
        return None

    def selected_label(self) -> str:
        """Short display name for the current combo choice (DB name or genome)."""
        data = self._combo.currentData()
        if not data:
            return reference_genome_combo_label(self._project, self._project_genome)
        text = self._combo.currentText().strip()
        if " (unavailable" in text:
            text = text.split(" (unavailable", 1)[0].strip()
        return text or Path(str(data)).stem

    def reset_to_default(self) -> None:
        if self._combo.count() > 0:
            self._combo.setCurrentIndex(0)

    def refresh(self, *, selected_path: str | None = None) -> None:
        """Fill combo like Blast Database (reference genome + nucleotide DBs)."""
        from repeatviewer.database_settings import reload_blast_databases_into_project

        # Preserve reference-genome selection (combo data is None) across refresh.
        prefer_current = False
        target = selected_path
        if self._combo.count() > 0 and selected_path is None:
            prefer_current = self._combo.currentData() is None
            target = self._combo.currentData()
        elif selected_path is not None and self._project_genome is not None:
            try:
                prefer_current = (
                    str(Path(str(selected_path)).expanduser().resolve())
                    == str(self._project_genome)
                )
            except OSError:
                prefer_current = False

        if self._project is not None:
            reload_blast_databases_into_project(self._project)

        self._combo.blockSignals(True)
        try:
            self._combo.clear()
            self._combo.addItem(
                reference_genome_combo_label(self._project, self._project_genome),
                None,
            )

            model = self._combo.model()
            databases = (
                list(self._project.blast_databases) if self._project is not None else []
            )
            for db in databases:
                if db.molecule_type != BLAST_MOLECULE_NUCL:
                    continue
                ready = blast_database_is_ready(db)
                label = db.name if ready else f"{db.name} (unavailable — index files missing)"
                data = str(Path(db.source_fasta).expanduser().resolve())
                self._combo.addItem(label, data)
                if not ready and isinstance(model, QStandardItemModel):
                    item = model.item(self._combo.count() - 1)
                    if item is not None:
                        item.setEnabled(False)

            if prefer_current or not target:
                self._combo.setCurrentIndex(0)
            else:
                resolved = str(Path(str(target)).expanduser().resolve())
                index = self._combo.findData(resolved)
                if index < 0 and self._project_genome is not None:
                    if resolved == str(self._project_genome):
                        index = 0
                if index < 0:
                    self._combo.addItem(Path(resolved).name, resolved)
                    index = self._combo.count() - 1
                self._combo.setCurrentIndex(index)
        finally:
            self._combo.blockSignals(False)
        self.selection_changed.emit()

    def _manage_databases(self) -> None:
        from repeatviewer.views.blast_database_dialog import BlastDatabaseManagerDialog

        if self._project is None:
            QMessageBox.information(
                self,
                self._manage_title,
                "No project is available to manage BLAST databases.",
            )
            return
        existing = self._blast_manager_dialog
        if existing is not None and existing.isVisible():
            existing.set_pick_mode(
                True,
                callback=self._on_database_picked,
                role_filter="blast",
                compatible="blast_nucl",
            )
            existing.raise_()
            existing.activateWindow()
            return
        # Parent to the owning dialog so Manage stays interactive under ApplicationModal.
        host: QWidget = self
        parent = self.parentWidget()
        while parent is not None:
            if isinstance(parent, QDialog):
                host = parent
                break
            parent = parent.parentWidget()
        dialog = BlastDatabaseManagerDialog(self._project, parent=host)
        self._blast_manager_dialog = dialog
        dialog.references_changed.connect(
            lambda: self.refresh(selected_path=self._combo.currentData())
        )
        dialog.finished.connect(self._on_manager_finished)
        dialog.set_pick_mode(
            True,
            callback=self._on_database_picked,
            role_filter="blast",
            compatible="blast_nucl",
        )
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _on_database_picked(self, payload) -> None:
        from repeatviewer.views.manage_references_dialog import (
            fasta_path_from_database_pick,
        )

        path = fasta_path_from_database_pick(payload)
        if path is None:
            QMessageBox.information(
                self,
                self._manage_title,
                "Select a nucleotide BLAST / genome FASTA database.",
            )
            return
        self.refresh(selected_path=str(path))

    def _on_manager_finished(self, _result: int = 0) -> None:
        self._blast_manager_dialog = None
        self.refresh(selected_path=self._combo.currentData())


class PanGenomeBlastListWidget(QWidget):
    """TE-Aid pangenome: stacked genome dropdowns + Manage Database + Add another."""

    size_changed = Signal()
    _MIN_GENOMES = 2

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        project: RepeatViewerProject | None = None,
        project_genome: Path | None = None,
        initial_paths: list[str] | None = None,
        help_title: str = "Genome",
        help_text: str | None = None,
        manage_title: str = "Manage Database",
    ) -> None:
        super().__init__(parent)
        self._project = project
        self._project_genome = (
            project_genome.expanduser().resolve()
            if project_genome is not None and project_genome.is_file()
            else None
        )
        self._help_title = help_title
        self._help_text = help_text
        self._manage_title = manage_title
        self._blast_manager_dialog = None
        self._rows: list[tuple[QWidget, GenomeBlastDatabaseRow, QPushButton]] = []
        # Grow with rows; never shrink genome combos to fit a short dialog.
        self.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Minimum
        )

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(6)

        self._rows_host = QWidget(self)
        self._rows_layout = QVBoxLayout(self._rows_host)
        self._rows_layout.setContentsMargins(0, 0, 0, 0)
        self._rows_layout.setSpacing(4)
        outer.addWidget(self._rows_host, 0)

        manage = QPushButton("Manage Database", self)
        manage.clicked.connect(self._manage_databases)
        outer.addWidget(manage, 0, Qt.AlignmentFlag.AlignLeft)

        add_btn = QPushButton("Add another genome", self)
        add_btn.setToolTip("Add another genome BLAST subject from Manage Database.")
        add_btn.clicked.connect(self._on_add_genome)
        outer.addWidget(add_btn, 0, Qt.AlignmentFlag.AlignLeft)

        paths = [p for p in (initial_paths or []) if str(p).strip()]
        while len(paths) < self._MIN_GENOMES:
            paths.append("")
        for index, path in enumerate(paths):
            selected = path or None
            if not selected and index > 0:
                selected = self._default_path_for_index(index)
            self._add_row(
                selected_path=selected,
                with_help=index == 0,
                emit_resize=False,
            )
        self._sync_remove_buttons()

    def _default_path_for_index(self, index: int) -> str | None:
        """Prefer a distinct nucleotide DB for row 2+ when available."""
        if self._project is None:
            return None
        databases = [
            db
            for db in self._project.blast_databases
            if db.molecule_type == BLAST_MOLECULE_NUCL and blast_database_is_ready(db)
        ]
        if not databases:
            return None
        pick = databases[min(index - 1, len(databases) - 1)]
        try:
            return str(Path(pick.source_fasta).expanduser().resolve())
        except OSError:
            return None

    def _add_row(
        self,
        *,
        selected_path: str | None = None,
        with_help: bool = False,
        emit_resize: bool = True,
    ) -> GenomeBlastDatabaseRow:
        wrap = QWidget(self._rows_host)
        wrap.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Fixed
        )
        row_layout = QHBoxLayout(wrap)
        row_layout.setContentsMargins(0, 0, 0, 0)
        row_layout.setSpacing(6)
        genome = GenomeBlastDatabaseRow(
            wrap,
            project=self._project,
            project_genome=self._project_genome,
            selected_path=selected_path,
            manage_title=self._manage_title,
            help_title=self._help_title if with_help else None,
            help_text=self._help_text if with_help else None,
            show_manage=False,
        )
        remove = QPushButton("Remove", wrap)
        remove.setToolTip("Remove this genome from the list")
        remove.setSizePolicy(
            QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed
        )
        remove.setFixedHeight(genome._combo.height())
        remove.clicked.connect(lambda: self._remove_row(wrap))
        row_layout.addWidget(genome, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        row_layout.addWidget(remove, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        row_layout.addStretch(1)
        self._rows.append((wrap, genome, remove))
        self._rows_layout.addWidget(wrap, 0)
        self._sync_remove_buttons()
        if emit_resize:
            self.size_changed.emit()
        return genome

    def _on_add_genome(self) -> None:
        self._add_row(with_help=False, emit_resize=True)

    def _remove_row(self, wrap: QWidget) -> None:
        if len(self._rows) <= self._MIN_GENOMES:
            return
        for index, (row_wrap, _genome, _remove) in enumerate(self._rows):
            if row_wrap is wrap:
                self._rows.pop(index)
                self._rows_layout.removeWidget(row_wrap)
                row_wrap.deleteLater()
                break
        self._sync_remove_buttons()
        self.size_changed.emit()

    def _sync_remove_buttons(self) -> None:
        can_remove = len(self._rows) > self._MIN_GENOMES
        for _wrap, _genome, remove in self._rows:
            remove.setEnabled(can_remove)
            remove.setVisible(can_remove)

    def refresh_all(self) -> None:
        for _wrap, genome, _remove in self._rows:
            genome.refresh(selected_path=genome.selected_path())

    def reset_to_defaults(self) -> None:
        while len(self._rows) > self._MIN_GENOMES:
            wrap, _genome, _remove = self._rows[-1]
            self._rows.pop()
            self._rows_layout.removeWidget(wrap)
            wrap.deleteLater()
        for _wrap, genome, _remove in self._rows:
            genome.reset_to_default()
        self._sync_remove_buttons()
        self.size_changed.emit()

    def snapshot_rows(self) -> list[tuple[str, str]]:
        """Return ``(name, path)`` pairs for filled genome rows."""
        from repeatviewer.pan_genome import normalize_pan_genome_name_loose

        out: list[tuple[str, str]] = []
        seen_names: set[str] = set()
        seen_paths: set[str] = set()
        for _wrap, genome, _remove in self._rows:
            path = genome.selected_path()
            if not path:
                continue
            try:
                resolved = str(Path(path).expanduser().resolve())
            except OSError:
                continue
            if resolved in seen_paths:
                continue
            base = normalize_pan_genome_name_loose(genome.selected_label())
            name = base
            suffix = 2
            while name in seen_names:
                # Keep within 10-char pan-genome name limit.
                trimmed = base[: max(1, 10 - len(str(suffix)))]
                name = normalize_pan_genome_name_loose(f"{trimmed}{suffix}")
                suffix += 1
            seen_names.add(name)
            seen_paths.add(resolved)
            out.append((name, resolved))
        return out

    def _manage_databases(self) -> None:
        from repeatviewer.views.blast_database_dialog import BlastDatabaseManagerDialog

        if self._project is None:
            QMessageBox.information(
                self,
                self._manage_title,
                "No project is available to manage BLAST databases.",
            )
            return
        existing = self._blast_manager_dialog
        if existing is not None and existing.isVisible():
            existing.set_pick_mode(
                True,
                callback=self._on_database_picked,
                role_filter="blast",
                compatible="blast_nucl",
            )
            existing.raise_()
            existing.activateWindow()
            return
        host: QWidget = self
        parent = self.parentWidget()
        while parent is not None:
            if isinstance(parent, QDialog):
                host = parent
                break
            parent = parent.parentWidget()
        dialog = BlastDatabaseManagerDialog(self._project, parent=host)
        self._blast_manager_dialog = dialog
        dialog.references_changed.connect(self.refresh_all)
        dialog.finished.connect(self._on_manager_finished)
        dialog.set_pick_mode(
            True,
            callback=self._on_database_picked,
            role_filter="blast",
            compatible="blast_nucl",
        )
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _on_database_picked(self, payload) -> None:
        from repeatviewer.views.manage_references_dialog import (
            fasta_path_from_database_pick,
        )

        path = fasta_path_from_database_pick(payload)
        if path is None:
            QMessageBox.information(
                self,
                self._manage_title,
                "Select a nucleotide BLAST / genome FASTA database.",
            )
            return
        self.refresh_all()
        target = None
        for _wrap, genome, _btn in self._rows:
            if genome._combo.currentData() is None:  # noqa: SLF001
                target = genome
                break
        if target is None and self._rows:
            target = self._rows[-1][1]
        if target is not None:
            target.refresh(selected_path=str(path))

    def _on_manager_finished(self, _result: int = 0) -> None:
        self._blast_manager_dialog = None
        self.refresh_all()


_HMMER_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): HMMER uses the FASTA files that "
    "were selected in the Files panel when you click Run.<br><br>"
    "<b>Define inputs</b>: choose different FASTA files. Paths are "
    "for this run only and are not remembered next time."
)
_HMMER_DATABASE_HELP = (
    "Profile HMM database from <b>File → Manage Database</b> (HMM DNA / HMM Protein). "
    "Protein HMMs search translated ORFs from nucleotide queries; "
    "DNA HMMs search nucleotide sequences directly."
)
_HMMER_EVALUE_HELP = (
    "Reporting sequence E-value threshold for pyhmmer <code>hmmscan</code> "
    "(<code>E</code> / <code>incE</code>)."
)
_HMMER_DOMAIN_EVALUE_HELP = (
    "Domain E-value threshold (<code>domE</code> / <code>incdomE</code>). "
    "Defaults to the same value as the sequence E-value."
)
_HMMER_THREADS_HELP = (
    "Number of CPU threads for pyhmmer <code>hmmscan</code> (modern API)."
)
_HMMER_MAX_HITS_HELP = (
    "Keep at most this many domain hits per query file, ranked by bit score."
)
_HMMER_WRITE_TEXT_HELP = (
    "Write a tab-separated hits table (<code>*_hmmer.txt</code>, "
    "<code>*_hmmer1.txt</code>, …) next to each query."
)
_HMMER_BIT_SCORE_HELP = (
    "Drop domain hits with bit score below this value (0 = no filter)."
)
_HMMER_MAX_ORFS_HELP = (
    "For protein HMMs on nucleotide queries: keep at most this many longest "
    "translated ORFs per sequence (0 = no limit)."
)


class HmmerParamsDialog(QDialog):
    """HMMER parameters (BLAST-style Inputs / Settings / Advanced / Options)."""

    _DIALOG_WIDTH = TOOL_PARAMS_DIALOG_WIDTH
    _SETTINGS_LABEL_REF = "Domain E-value:"

    def __init__(
        self,
        *,
        initial: HmmerOptions | None = None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("HMMER parameters")
        self.setMinimumWidth(self._DIALOG_WIDTH)
        opts = initial or HmmerOptions()
        self._for_run = for_run
        self._hmmer_compact_size: QSize | None = None
        self._define_database_dialog = None

        self._input_group = QButtonGroup(self)
        self._input_selected = QRadioButton("Use selected sequences", self)
        self._input_define = QRadioButton("Define inputs", self)
        self._input_group.addButton(self._input_selected)
        self._input_group.addButton(self._input_define)
        self._input_selected.setChecked(True)
        self._define_inputs = TeMultiInputFilesWidget(
            self,
            initial=None,
            start_empty=False,
            placeholder="Click in Files/Folders or Browse…",
            browse_title="Choose a FASTA file",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
            ),
            add_button_label="Add another input",
            add_button_tooltip=(
                "Add another FASTA file path. Click a file in the Files panel "
                "to fill a row, or use Browse."
            ),
            path_mode="file",
        )
        self._input_selected.toggled.connect(self._sync_hmmer_inputs)
        self._input_define.toggled.connect(self._sync_hmmer_inputs)
        self._define_inputs.paths_changed.connect(self._sync_hmmer_inputs)

        self._database = QComboBox(self)
        self._database.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        self._populate_database_combo(selected_name=opts.database_name)

        self._evalue = QDoubleSpinBox(self)
        configure_report_evalue_spinbox(self._evalue)
        self._evalue.setValue(normalize_blast_evalue(opts.e_value))
        self._evalue.setToolTip(_HMMER_EVALUE_HELP)

        self._threads = QSpinBox(self)
        self._threads.setRange(1, 128)
        self._threads.setValue(max(1, opts.num_threads))
        self._threads.setToolTip(_HMMER_THREADS_HELP)

        self._max_hits = QSpinBox(self)
        self._max_hits.setRange(1, 100_000)
        self._max_hits.setValue(max(1, opts.max_hits))
        self._max_hits.setToolTip(_HMMER_MAX_HITS_HELP)

        self._write_text = QCheckBox("Export hits table", self)
        self._write_text.setChecked(bool(opts.write_text))
        self._write_text.setToolTip(_HMMER_WRITE_TEXT_HELP)

        self._domain_evalue = QDoubleSpinBox(self)
        configure_report_evalue_spinbox(self._domain_evalue)
        self._domain_evalue.setValue(normalize_blast_evalue(opts.domain_e_value))
        self._domain_evalue.setToolTip(_HMMER_DOMAIN_EVALUE_HELP)

        self._min_bit = QDoubleSpinBox(self)
        self._min_bit.setRange(0.0, 10_000.0)
        self._min_bit.setDecimals(1)
        self._min_bit.setValue(max(0.0, opts.minimum_bit_score))
        self._min_bit.setToolTip(_HMMER_BIT_SCORE_HELP)

        self._max_orfs = QSpinBox(self)
        self._max_orfs.setRange(0, 100_000)
        self._max_orfs.setSpecialValueText("No limit")
        self._max_orfs.setValue(max(0, opts.max_orfs_per_sequence))
        self._max_orfs.setToolTip(_HMMER_MAX_ORFS_HELP)

        for field in (
            self._evalue,
            self._threads,
            self._max_hits,
            self._domain_evalue,
            self._min_bit,
            self._max_orfs,
        ):
            field.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)

        params_tab = QWidget(self)
        params_root = QVBoxLayout(params_tab)
        params_root.setContentsMargins(16, 10, 16, 10)
        params_root.setSpacing(4)
        params_root.setAlignment(Qt.AlignmentFlag.AlignTop)

        basic = QWidget(params_tab)
        self._hmmer_basic = basic
        basic.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        form = QFormLayout(basic)
        form.setContentsMargins(0, 0, 0, 0)
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        # Expanding so Define-inputs path + Browse fill the dialog width.
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)
        self._hmmer_basic_form = form

        def _section(title: str) -> QLabel:
            return QLabel(f"<b>{title}</b>", params_tab)

        form.addRow(_section("Inputs"))
        inputs_box = QWidget(basic)
        inputs_layout = QVBoxLayout(inputs_box)
        inputs_layout.setContentsMargins(0, 0, 0, 0)
        inputs_layout.setSpacing(2)
        inputs_layout.addWidget(
            field_with_help(
                self._input_selected,
                "HMMER inputs",
                _HMMER_INPUTS_HELP,
                stretch=True,
            )
        )
        inputs_layout.addWidget(self._input_define)
        self._inputs_radios_label = QLabel("")
        self._inputs_define_label = QLabel("")
        form.addRow(self._inputs_radios_label, inputs_box)
        form.addRow(self._inputs_define_label, self._define_inputs)

        form.addRow(_form_section_separator())
        form.addRow(_section("Settings"))
        db_field = QWidget(basic)
        db_layout = QVBoxLayout(db_field)
        db_layout.setContentsMargins(0, 0, 0, 0)
        db_layout.setSpacing(4)
        db_layout.addWidget(
            field_with_help(self._database, "HMM database", _HMMER_DATABASE_HELP)
        )
        manage = QPushButton("Manage Database…", db_field)
        manage.setCursor(Qt.CursorShape.PointingHandCursor)
        manage.clicked.connect(self._open_define_database)
        db_layout.addWidget(manage, 0, Qt.AlignmentFlag.AlignLeft)
        form.addRow("Database:", db_field)
        form.addRow(
            "E-value:",
            field_with_help(self._evalue, "E-value", _HMMER_EVALUE_HELP),
        )
        form.addRow(
            "Threads:",
            field_with_help(self._threads, "Threads", _HMMER_THREADS_HELP),
        )
        form.addRow(
            "Max hits:",
            field_with_help(self._max_hits, "Max hits", _HMMER_MAX_HITS_HELP),
        )
        form.addRow("", self._write_text)
        params_root.addWidget(basic)

        advanced = QWidget(params_tab)
        self._hmmer_advanced_panel = advanced
        advanced.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        advanced_form = QFormLayout(advanced)
        advanced_form.setContentsMargins(0, 4, 0, 0)
        advanced_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        advanced_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        advanced_form.setHorizontalSpacing(12)
        advanced_form.setVerticalSpacing(4)
        advanced_form.addRow(_form_section_separator())
        advanced_form.addRow(QLabel("<b>Advanced Settings</b>", advanced))
        advanced_form.addRow(
            "Domain E-value:",
            field_with_help(
                self._domain_evalue, "Domain E-value", _HMMER_DOMAIN_EVALUE_HELP
            ),
        )
        advanced_form.addRow(
            "Minimum bit score:",
            field_with_help(
                self._min_bit, "Minimum bit score", _HMMER_BIT_SCORE_HELP
            ),
        )
        advanced_form.addRow(
            "Max ORFs per sequence:",
            field_with_help(
                self._max_orfs, "Max ORFs per sequence", _HMMER_MAX_ORFS_HELP
            ),
        )
        self._hmmer_advanced_form = advanced_form
        params_root.addWidget(advanced)
        self._params_tab = params_tab

        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("Run HMMER with saved parameters", HMMER_RIGHT_CLICK_RUN),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[HMMER_RIGHT_CLICK_RUN].setChecked(True)

        rc_box = QGroupBox("Right-click action", self)
        rc_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(QLabel("Right-clicking the HMMER button does:", self))
        for radio in self._right_click_buttons.values():
            rc_box_layout.addWidget(radio)

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="HMMER",
                checked=bool(opts.show_process_window),
                close_when_finished=bool(opts.close_process_window_when_finished),
            )
        )

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)
        options_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        options_layout.addWidget(rc_box)
        options_layout.addWidget(process_box)
        options_layout.addStretch(1)
        self._options_tab = options_tab

        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        self._tabs.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._tabs.addTab(params_tab, "Parameters")
        self._tabs.addTab(options_tab, "Options")
        self._tabs.currentChanged.connect(self._on_hmmer_tab_changed)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._advanced_settings_btn = QPushButton("Advanced Settings ▸", self)
        self._advanced_settings_btn.setCheckable(True)
        self._advanced_settings_btn.setChecked(False)
        self._advanced_settings_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._advanced_settings_btn.toggled.connect(self._on_advanced_settings_toggled)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore, self._advanced_settings_btn],
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(14)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._tabs, 0)
        layout.addWidget(self._button_box, 0)

        open_tab = opening_dialog_tab(
            opts.last_dialog_tab,
            allowed=HMMER_DIALOG_TABS,
            first_tab=HMMER_DIALOG_TAB_PARAMETERS,
        )
        self._tabs.setCurrentWidget(
            options_tab if open_tab == DIALOG_TAB_OPTIONS else params_tab
        )

        self._sync_hmmer_inputs(resize=False)
        self._set_hmmer_advanced_visible(False)
        self._sync_field_column()
        self.resize(self._DIALOG_WIDTH, 420)
        QTimer.singleShot(0, lambda: self._resize_for_hmmer_advanced(False))

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._hmmer_compact_size = None
        self._populate_database_combo(
            selected_name=self._database.currentData()
            if self._database.currentIndex() >= 0
            else None
        )
        self._sync_field_column()
        QTimer.singleShot(0, lambda: self._resize_for_hmmer_advanced(False))

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def current_dialog_tab(self) -> str:
        if self._tabs.currentWidget() is self._options_tab:
            return DIALOG_TAB_OPTIONS
        return HMMER_DIALOG_TAB_PARAMETERS

    def input_source(self) -> str:
        if self._input_define.isChecked():
            return ALIGN_INPUT_DEFINE
        return ALIGN_INPUT_SELECTED

    def run_files(self, *, panel_files: list[Path] | None = None) -> list[Path]:
        from repeatviewer.backends.fastga_dotplot import expand_fastga_input_paths

        if self.input_source() == ALIGN_INPUT_DEFINE:
            return expand_fastga_input_paths(self._define_inputs.paths())
        return list(panel_files or [])

    def apply_path_from_panel(self, path: Path) -> bool:
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        return self._define_inputs.accept_path(Path(path))

    def apply_folder_from_panel(self, folder: Path) -> bool:
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        return self._define_inputs.accept_path(Path(folder))

    def _populate_database_combo(self, *, selected_name: str | None = None) -> None:
        from repeatviewer.database_settings import load_user_databases_from_settings

        current = selected_name
        if current is None and self._database.count():
            current = self._database.currentData()
        self._database.clear()
        self._database.addItem("Select an HMM database…", "")
        for db in load_user_databases_from_settings():
            if not db.is_hmm:
                continue
            label = f"{db.name} ({db.display_type_label})"
            self._database.addItem(label, db.name)
        if current:
            idx = self._database.findData(current)
            if idx >= 0:
                self._database.setCurrentIndex(idx)

    def _open_define_database(self) -> None:
        from repeatviewer.views.manage_references_dialog import ManageReferencesDialog

        if self._define_database_dialog is not None and self._define_database_dialog.isVisible():
            self._define_database_dialog.set_pick_mode(
                True,
                callback=self._on_hmm_database_picked,
                role_filter="annotation",
                compatible="hmm",
            )
            self._define_database_dialog.raise_()
            self._define_database_dialog.activateWindow()
            return
        host: QWidget = self
        parent = self.parentWidget()
        project = None
        while parent is not None:
            if isinstance(parent, QDialog):
                host = parent
            project = getattr(parent, "_project", project)
            parent = parent.parentWidget()
        if project is None:
            project = getattr(self.parent(), "_project", None)
        if project is None:
            from repeatviewer.project import RepeatViewerProject

            project = RepeatViewerProject()
        dialog = ManageReferencesDialog(project, host, role_filter="annotation")
        dialog.setModal(False)
        dialog.setWindowModality(Qt.WindowModality.NonModal)
        dialog.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, False)
        self._define_database_dialog = dialog
        dialog.references_changed.connect(
            lambda: self._populate_database_combo(
                selected_name=self._database.currentData()
                if self._database.currentIndex() >= 0
                else None
            )
        )
        dialog.finished.connect(self._on_define_database_finished)
        dialog.set_pick_mode(
            True,
            callback=self._on_hmm_database_picked,
            role_filter="annotation",
            compatible="hmm",
        )
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _on_hmm_database_picked(self, payload) -> None:
        from repeatviewer.views.manage_references_dialog import database_name_from_pick
        from repeatviewer.project import UserDatabase

        name = database_name_from_pick(payload)
        if not name and isinstance(payload, UserDatabase):
            name = str(payload.name or "").strip()
        if not name:
            QMessageBox.information(
                self,
                "HMMER",
                "Select an HMM DNA or HMM Protein database.",
            )
            return
        self._populate_database_combo(selected_name=name)

    def _on_define_database_finished(self, _result: int = 0) -> None:
        self._define_database_dialog = None
        self._populate_database_combo(
            selected_name=self._database.currentData()
            if self._database.currentIndex() >= 0
            else None
        )

    def _sync_hmmer_inputs(self, *, resize: bool = True) -> None:
        define = self.input_source() == ALIGN_INPUT_DEFINE
        form = self._hmmer_basic.layout()
        if isinstance(form, QFormLayout):
            form.setRowVisible(self._define_inputs, define)
        self._define_inputs.setVisible(define)
        self._inputs_define_label.setVisible(define)
        if define:
            self._define_inputs.setMinimumHeight(0)
            self._define_inputs.setMaximumHeight(16777215)
        else:
            self._define_inputs.setFixedHeight(0)
        self._define_inputs.updateGeometry()
        self._hmmer_basic.updateGeometry()
        self._sync_field_column()
        if resize and self.isVisible() and not self._advanced_settings_btn.isChecked():
            self._hmmer_compact_size = None
            self._resize_for_hmmer_advanced(False)

    def _sync_field_column(self) -> None:
        sync_form_field_column_spacers(
            self,
            [self._inputs_radios_label, self._inputs_define_label],
            [
                getattr(self, "_hmmer_basic_form", None),
                getattr(self, "_hmmer_advanced_form", None),
            ],
            ref_label=self._SETTINGS_LABEL_REF,
        )

    def _set_hmmer_advanced_visible(self, visible: bool) -> None:
        self._hmmer_advanced_panel.setVisible(visible)
        if visible:
            self._hmmer_advanced_panel.setMaximumHeight(16777215)
        else:
            self._hmmer_advanced_panel.setMaximumHeight(0)
        self._hmmer_advanced_panel.updateGeometry()
        root = self._params_tab.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        self._params_tab.updateGeometry()

    def _on_advanced_settings_toggled(self, checked: bool) -> None:
        self._advanced_settings_btn.setText(
            "Advanced Settings ▾" if checked else "Advanced Settings ▸"
        )
        if self._tabs.currentWidget() is not self._params_tab:
            return
        if checked:
            self._hmmer_compact_size = self.size()
        self._set_hmmer_advanced_visible(checked)
        self._resize_for_hmmer_advanced(checked)

    def _on_hmmer_tab_changed(self, _index: int) -> None:
        on_params = self._tabs.currentWidget() is self._params_tab
        self._advanced_settings_btn.setVisible(on_params)
        if not on_params:
            if self._advanced_settings_btn.isChecked():
                self._advanced_settings_btn.blockSignals(True)
                self._advanced_settings_btn.setChecked(False)
                self._advanced_settings_btn.setText("Advanced Settings ▸")
                self._advanced_settings_btn.blockSignals(False)
                self._set_hmmer_advanced_visible(False)
            self._fit_options_content()
        else:
            self._resize_for_hmmer_advanced(self._advanced_settings_btn.isChecked())

    def _dialog_chrome_height(self) -> int:
        lay = self.layout()
        margins = lay.contentsMargins() if lay is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 34
        spacing = lay.spacing() if lay is not None else 14
        tab_bar_h = self._tabs.tabBar().sizeHint().height() + 10
        buttons_h = self._button_box.sizeHint().height()
        return margin_h + spacing + tab_bar_h + buttons_h

    def _hmmer_compact_content_height(self) -> int:
        root = self._params_tab.layout()
        margins = root.contentsMargins() if root is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 20
        spacing = root.spacing() if root is not None else 4
        basic_h = max(
            self._hmmer_basic.minimumSizeHint().height(),
            self._hmmer_basic.sizeHint().height(),
        )
        return margin_h + basic_h + spacing

    def _resize_for_hmmer_advanced(self, expanded: bool) -> None:
        if self._tabs.currentWidget() is not self._params_tab:
            self._fit_options_content()
            return
        panel = self._params_tab
        adv = self._hmmer_advanced_panel
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        root = panel.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        panel.updateGeometry()
        self._tabs.updateGeometry()

        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        chrome = self._dialog_chrome_height()
        width = tool_params_dialog_width(
            define_inputs=self.input_source() == ALIGN_INPUT_DEFINE
        )
        self.setMinimumWidth(width)
        compact_h = self._hmmer_compact_content_height()

        if not expanded:
            height = chrome + compact_h
            height = min(max(height, 220), max_h)
            pin_tool_tabs_height(self._tabs, compact_h)
            self.setMaximumHeight(height)
            self.resize(width, height)
            self.setMaximumHeight(16777215)
            self._hmmer_compact_size = QSize(width, height)
            return

        extra = max(0, adv.sizeHint().height())
        page_h = compact_h + extra
        base = self._hmmer_compact_size or self.size()
        height = min(base.height() + extra + 12, max_h)
        height = max(height, min(chrome + page_h, max_h))
        pin_tool_tabs_height(self._tabs, page_h)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def _fit_options_content(self) -> None:
        page = self._options_tab
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        page.ensurePolished()
        page.adjustSize()
        root = page.layout()
        if root is not None:
            root.activate()
        page_hint = page.minimumSizeHint().expandedTo(page.sizeHint())
        content_h = page_hint.height()
        if root is not None:
            content_h = max(
                content_h,
                root.totalMinimumSize().height(),
                root.totalSizeHint().height(),
            )
        chrome = self._dialog_chrome_height()
        height = chrome + content_h
        pin_tool_tabs_height(self._tabs, content_h)
        self.setMaximumHeight(height)
        self.resize(self._DIALOG_WIDTH, height)
        self.setMaximumHeight(16777215)

    def _restore_defaults(self) -> None:
        defaults = HmmerOptions()
        self._input_selected.setChecked(True)
        self._define_inputs.clear_paths()
        self._sync_hmmer_inputs()
        self._populate_database_combo(selected_name="")
        self._evalue.setValue(DEFAULT_HMMER_EVALUE)
        self._threads.setValue(defaults.num_threads)
        self._max_hits.setValue(DEFAULT_HMMER_MAX_HITS)
        self._write_text.setChecked(DEFAULT_HMMER_WRITE_TEXT)
        self._domain_evalue.setValue(DEFAULT_HMMER_DOMAIN_EVALUE)
        self._min_bit.setValue(DEFAULT_HMMER_MINIMUM_BIT_SCORE)
        self._max_orfs.setValue(DEFAULT_HMMER_MAX_ORFS)
        self._show_process_window.setChecked(defaults.show_process_window)
        self._close_process_window.setChecked(
            defaults.close_process_window_when_finished
        )
        self._right_click_buttons[HMMER_RIGHT_CLICK_RUN].setChecked(True)

    def options(self) -> HmmerOptions:
        right_click_action = HMMER_RIGHT_CLICK_RUN
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click_action = value
                break
        db_name = self._database.currentData()
        return HmmerOptions(
            database_name=str(db_name or "").strip(),
            e_value=normalize_blast_evalue(self._evalue.value()),
            domain_e_value=normalize_blast_evalue(self._domain_evalue.value()),
            num_threads=max(1, self._threads.value()),
            max_hits=max(1, self._max_hits.value()),
            write_text=self._write_text.isChecked(),
            minimum_bit_score=max(0.0, float(self._min_bit.value())),
            max_orfs_per_sequence=max(0, self._max_orfs.value()),
            input_source=ALIGN_INPUT_SELECTED,
            input_paths=[],
            last_dialog_tab=persistable_dialog_tab(
                self.current_dialog_tab(),
                first_tab=HMMER_DIALOG_TAB_PARAMETERS,
            ),
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=self._close_process_window.isChecked(),
            right_click_action=right_click_action,
        )


def _make_extend_length_combo(
    initial: int, parent: QWidget | None = None
) -> QComboBox:
    """Editable dropdown: pick a preset length or type a custom bp value."""
    combo = QComboBox(parent)
    combo.setEditable(True)
    combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
    combo.setSizeAdjustPolicy(
        QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
    )
    combo.setMinimumContentsLength(8)
    for value in EXTEND_LENGTH_PRESETS:
        combo.addItem(f"{value} bp", int(value))
    value = max(0, min(500_000, int(initial)))
    index = combo.findData(value)
    if index >= 0:
        combo.setCurrentIndex(index)
    else:
        combo.setEditText(f"{value} bp")
    combo.setToolTip(
        "Choose a preset length or type a custom number of base pairs "
        "(digits, optionally with a bp suffix)."
    )
    return combo


def _extend_length_from_combo(combo: QComboBox) -> int:
    text = (combo.currentText() or "").strip().lower().replace(",", "")
    for suffix in (" bp", "bp"):
        if text.endswith(suffix):
            text = text[: -len(suffix)].strip()
            break
    if text.isdigit() or (text.startswith("-") and text[1:].isdigit()):
        return max(0, min(500_000, int(text)))
    data = combo.currentData()
    if isinstance(data, int):
        return max(0, min(500_000, data))
    return 0


class ExtendParamsDialog(QDialog):
    """Extend MSA ends.

    Left-click the toolbar button opens this dialog. The Parameters tab holds
    the distances and genome; the Options tab sets what right-clicking the
    toolbar button does.
    """

    def __init__(
        self,
        *,
        initial: ExtendOptions | None = None,
        title: str = "MSA Extend parameters",
        pan_genome_names: list[str] | None = None,
        project: RepeatViewerProject | None = None,
        project_genome: Path | None = None,
        genome_history: list[Path] | None = None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        opts = initial or ExtendOptions()
        self._for_run = for_run
        self._run_after = False
        # genome_history kept for call-site compatibility; Genome list uses BLAST DBs.
        _ = genome_history

        self._left = _make_extend_length_combo(opts.left, self)
        self._right = _make_extend_length_combo(opts.right, self)

        self._pan_genome = QComboBox()
        self._pan_genome_names = list(pan_genome_names or [])
        if self._pan_genome_names:
            self._pan_genome.addItems(self._pan_genome_names)
            selected = opts.pan_genome_name or self._pan_genome_names[0]
            index = self._pan_genome.findText(selected)
            if index >= 0:
                self._pan_genome.setCurrentIndex(index)

        self._genome_row = GenomeBlastDatabaseRow(
            self,
            project=project,
            project_genome=project_genome,
            selected_path=opts.genome_path,
            manage_title="MSA Extend parameters",
        )

        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("Run MSA Extend with saved parameters", EXTEND_RIGHT_CLICK_RUN),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[EXTEND_RIGHT_CLICK_RUN].setChecked(True)

        rc_box = QGroupBox("Right-click action", self)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(
            QLabel("Right-clicking the MSA Extend button does:", self)
        )
        for radio in self._right_click_buttons.values():
            rc_box_layout.addWidget(radio)

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="MSA Extend",
                checked=bool(opts.show_process_window),
                close_when_finished=bool(opts.close_process_window_when_finished),
            )
        )

        options_tab = QWidget(self)
        options_tab.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(12, 10, 12, 10)
        options_layout.setSpacing(10)
        options_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        options_layout.addWidget(rc_box)
        options_layout.addWidget(process_box)

        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        self._tabs.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        params_tab = QWidget(self)
        params_tab.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        params_form = QFormLayout(params_tab)
        params_form.setContentsMargins(12, 10, 12, 10)
        params_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        params_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        params_form.setHorizontalSpacing(12)
        params_form.setVerticalSpacing(8)
        if self._pan_genome_names:
            params_form.addRow("Pan-genome:", self._pan_genome)
        else:
            params_form.addRow("Reference genome:", self._genome_row)
        params_form.addRow("Extend left of MSA:", self._left)
        params_form.addRow("Extend right of MSA:", self._right)
        self._tabs.addTab(params_tab, "Parameters")
        self._tabs.addTab(options_tab, "Options")

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore],
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 10, 12, 10)
        layout.setSpacing(8)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._tabs, 0)
        layout.addWidget(self._button_box, 0)
        self._tabs.currentChanged.connect(self._on_tab_changed)
        self._on_tab_changed(self._tabs.currentIndex())

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        # Dialog polish may call adjustSize(); re-fit to the active tab after that.
        QTimer.singleShot(0, self._fit_to_current_tab)

    def _on_tab_changed(self, _tab_index: int) -> None:
        self._fit_to_current_tab()

    def _fit_to_current_tab(self) -> None:
        """Size the dialog to the active tab only (QTabWidget otherwise uses max of all)."""
        page = self._tabs.currentWidget()
        if page is None:
            return
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        page.ensurePolished()
        page.adjustSize()
        root = page.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        page_hint = page.minimumSizeHint().expandedTo(page.sizeHint())
        content_h = page_hint.height()
        if root is not None:
            content_h = max(
                content_h,
                root.totalMinimumSize().height(),
                root.totalSizeHint().height(),
            )
        pin_tool_tabs_height(self._tabs, content_h)
        tab_bar_h = max(24, self._tabs.tabBar().sizeHint().height())
        tabs_h = tab_bar_h + content_h + 8
        buttons_h = self._button_box.sizeHint().height()
        margins = self.layout().contentsMargins()
        spacing = self.layout().spacing()
        height = (
            tabs_h
            + buttons_h
            + margins.top()
            + margins.bottom()
            + spacing
            + 4
        )
        width = max(
            420,
            page_hint.width() + margins.left() + margins.right() + 24,
            self._button_box.sizeHint().width() + margins.left() + margins.right(),
        )
        self.setMinimumWidth(width)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _restore_defaults(self) -> None:
        defaults = ExtendOptions()
        left_index = self._left.findData(DEFAULT_EXTEND_LEFT)
        if left_index >= 0:
            self._left.setCurrentIndex(left_index)
        else:
            self._left.setEditText(f"{DEFAULT_EXTEND_LEFT} bp")
        right_index = self._right.findData(DEFAULT_EXTEND_RIGHT)
        if right_index >= 0:
            self._right.setCurrentIndex(right_index)
        else:
            self._right.setEditText(f"{DEFAULT_EXTEND_RIGHT} bp")
        self._genome_row.reset_to_default()
        self._right_click_buttons[EXTEND_RIGHT_CLICK_RUN].setChecked(True)
        self._show_process_window.setChecked(defaults.show_process_window)
        self._close_process_window.setChecked(
            defaults.close_process_window_when_finished
        )
        self._close_process_window.setEnabled(self._show_process_window.isChecked())

    def options(self) -> ExtendOptions:
        pan_name = None
        if self._pan_genome_names:
            pan_name = self._pan_genome.currentText().strip() or None
        genome_path = None
        if not self._pan_genome_names:
            genome_path = self._genome_row.selected_path()
        right_click_action = EXTEND_RIGHT_CLICK_RUN
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click_action = value
                break
        return ExtendOptions(
            left=_extend_length_from_combo(self._left),
            right=_extend_length_from_combo(self._right),
            pan_genome_name=pan_name,
            genome_path=genome_path,
            right_click_action=right_click_action,
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=self._close_process_window.isChecked(),
        )


class TsdSearchParamsDialog(QDialog):
    def __init__(
        self,
        *,
        initial: TsdSearchOptions | None = None,
        title: str = "TSD search parameters",
        parent=None,
        include_clustering: bool = True,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        opts = initial or TsdSearchOptions()
        self._include_clustering = bool(include_clustering)
        self._saved_left = max(0, int(opts.left_boundary))
        self._saved_right = max(0, int(opts.right_boundary))

        self._boundary_source = QComboBox()
        self._boundary_source.addItem("Auto detect", TSD_BOUNDARY_REPEATVIEWER)
        self._boundary_source.addItem("User define", TSD_BOUNDARY_USER)
        source = (
            opts.boundary_source
            if opts.boundary_source in TSD_BOUNDARY_SOURCES
            else TSD_BOUNDARY_REPEATVIEWER
        )
        idx = self._boundary_source.findData(source)
        self._boundary_source.setCurrentIndex(max(0, idx))
        from repeatviewer.backends.tsd_search import TSD_BOUNDARY_DEFINITION_HELP
        from repeatviewer.widgets.help_button import help_tooltip_text

        self._boundary_source.setToolTip(
            help_tooltip_text(TSD_BOUNDARY_DEFINITION_HELP)
        )

        self._radius = QSpinBox()
        self._radius.setRange(1, 10_000_000)
        self._radius.setValue(max(1, int(opts.nt_radius)))
        self._radius.setToolTip(
            "MSA columns to search on both sides of each TE edge "
            "(gaps count as positions)."
        )

        self._min_len = QSpinBox()
        self._min_len.setRange(2, 200)
        self._min_len.setValue(opts.min_tsd_length)
        self._min_len.setToolTip("Minimum TSD motif length.")

        self._max_len = QSpinBox()
        self._max_len.setRange(2, 200)
        self._max_len.setValue(opts.max_tsd_length)
        self._max_len.setToolTip("Maximum TSD motif length.")

        self._min_dist = QSpinBox()
        self._min_dist.setRange(0, 10_000_000)
        self._min_dist.setValue(opts.min_distance)
        self._min_dist.setToolTip(
            "Minimum distance between the left and right TSD copies "
            "(TE body length). Default 200 bp. 0 = no minimum."
        )
        self._saved_max_distance = max(0, int(opts.max_distance))

        self._write_text_report = QCheckBox("Write .tsd_search.txt report")
        self._write_text_report.setChecked(bool(opts.write_text_report))
        self._write_text_report.setToolTip(
            "When checked, write a clustered TSD text report beside the MSA "
            "(shows in Files). Off by default."
        )

        form = QFormLayout()
        form.addRow("Boundary definition:", self._boundary_source)
        form.addRow("Search range around boundary (bp):", self._radius)
        form.addRow("Min TSD length (bp):", self._min_len)
        form.addRow("Max TSD length (bp):", self._max_len)
        form.addRow("Min distance between TSDs (bp):", self._min_dist)
        form.addRow("", self._write_text_report)

        if self._include_clustering:
            self._max_sim = QDoubleSpinBox()
            self._max_sim.setRange(0.0, 1.0)
            self._max_sim.setDecimals(2)
            self._max_sim.setSingleStep(0.05)
            self._max_sim.setValue(opts.max_sim_threshold)
            self._max_sim.setToolTip(
                "Discard a TSD location when average motif similarity across "
                "copies exceeds this value (Jaccard on 2-mers). True genomic "
                "TSDs differ between insertions; identical flanks are rejected."
            )

            self._nt_threshold = QDoubleSpinBox()
            self._nt_threshold.setRange(0.0, 1.0)
            self._nt_threshold.setDecimals(2)
            self._nt_threshold.setSingleStep(0.05)
            self._nt_threshold.setValue(opts.nt_threshold)
            self._nt_threshold.setToolTip(
                "Clustering: poly-N / low-complexity filter threshold."
            )

            self._cluster_eps = QDoubleSpinBox()
            self._cluster_eps.setRange(0.1, 10.0)
            self._cluster_eps.setDecimals(1)
            self._cluster_eps.setSingleStep(0.1)
            self._cluster_eps.setValue(opts.cluster_eps)
            self._cluster_eps.setToolTip(
                "DBSCAN cluster distance (eps) for the Tools text report."
            )

            form.addRow("Max TSD similarity:", self._max_sim)
            form.addRow("Poly-N threshold:", self._nt_threshold)
            form.addRow("Cluster distance (eps):", self._cluster_eps)
        else:
            self._max_sim = None
            self._nt_threshold = None
            self._cluster_eps = None

        self._button_box = make_cancel_save_run_bar(self)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self._button_box)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def options(self) -> TsdSearchOptions:
        min_len = int(self._min_len.value())
        max_len = max(min_len, int(self._max_len.value()))
        defaults = TsdSearchOptions()
        source = self._boundary_source.currentData()
        if source not in TSD_BOUNDARY_SOURCES:
            source = TSD_BOUNDARY_REPEATVIEWER
        return TsdSearchOptions(
            left_boundary=self._saved_left,
            right_boundary=self._saved_right,
            nt_radius=int(self._radius.value()),
            min_tsd_length=min_len,
            max_tsd_length=max_len,
            min_distance=self._min_dist.value(),
            max_distance=self._saved_max_distance,
            max_sim_threshold=(
                self._max_sim.value()
                if self._max_sim is not None
                else defaults.max_sim_threshold
            ),
            nt_threshold=(
                self._nt_threshold.value()
                if self._nt_threshold is not None
                else defaults.nt_threshold
            ),
            cluster_eps=(
                self._cluster_eps.value()
                if self._cluster_eps is not None
                else defaults.cluster_eps
            ),
            boundary_source=str(source),
            write_text_report=self._write_text_report.isChecked(),
        )


class TsdBoundaryDialog(QDialog):
    """User-define left/right TE edges before a TSD scan.

    Drag a column range on the Overview MSA to fill Start / End. Prefills from
    RepeatViewer boundaries when available.
    """

    def __init__(
        self,
        *,
        left_boundary: int = 0,
        right_boundary: int = 0,
        nt_radius: int = 100,
        aln_len: int = 1,
        parent=None,
        title: str = "TSD search boundaries",
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        # Must stay non-modal and use show()+event-loop (not exec()) so the
        # Overview MSA can receive drag events while this dialog is open.
        # Never use WindowStaysOnTopHint (covers other apps); raise_ on MSA drag.
        self.setModal(False)
        self.setWindowModality(Qt.WindowModality.NonModal)
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, False)
        self._aln_len = max(1, int(aln_len))
        self._updating = False

        hint = QLabel("Drag a column range on the MSA to set Start and End", self)
        hint.setWordWrap(True)

        self._left = QSpinBox(self)
        self._left.setRange(0, self._aln_len)
        self._left.setSpecialValueText("unset")
        self._left.setValue(max(0, min(self._aln_len, int(left_boundary))))
        self._left.setToolTip("1-based left TE edge (0 = unset / whole MSA).")

        self._right = QSpinBox(self)
        self._right.setRange(0, self._aln_len)
        self._right.setSpecialValueText("unset")
        self._right.setValue(max(0, min(self._aln_len, int(right_boundary))))
        self._right.setToolTip("1-based right TE edge (0 = unset / whole MSA).")

        self._radius = QSpinBox(self)
        self._radius.setRange(1, 10_000_000)
        self._radius.setValue(max(1, int(nt_radius)))
        self._radius.setToolTip(
            "MSA columns to search on both sides of each TE edge "
            "(gaps count as positions)."
        )

        form = QFormLayout()
        form.addRow("Start (left boundary):", self._left)
        form.addRow("End (right boundary):", self._right)
        form.addRow("Search range around boundary (bp):", self._radius)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel,
            self,
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addWidget(hint)
        layout.addLayout(form)
        layout.addWidget(buttons)
        self.resize(420, 200)

    def set_from_msa_selection(self, start: int, end: int) -> None:
        """Fill Start/End from a 1-based inclusive MSA column selection."""
        if self._updating:
            return
        a = max(1, min(self._aln_len, int(start)))
        b = max(1, min(self._aln_len, int(end)))
        if b < a:
            a, b = b, a
        self._updating = True
        try:
            self._left.setValue(a)
            self._right.setValue(b)
        finally:
            self._updating = False
        # Keep this dialog in front of RepeatViewer while the user drags the MSA
        # (not above other applications).
        self.raise_()

    def left_boundary(self) -> int:
        return int(self._left.value())

    def right_boundary(self) -> int:
        return int(self._right.value())

    def nt_radius(self) -> int:
        return int(self._radius.value())

    def apply_to_options(self, opts: TsdSearchOptions) -> TsdSearchOptions:
        """Return ``opts`` with confirmed boundaries / radius."""
        from dataclasses import replace

        return replace(
            opts,
            left_boundary=self.left_boundary(),
            right_boundary=self.right_boundary(),
            nt_radius=self.nt_radius(),
        )


def _auto_resolve_tsd_boundaries(
    opts: TsdSearchOptions,
    *,
    file_path: Path | str,
    aln_len: int,
    reference_extent_1based: tuple[int, int] | None = None,
    reference_sequence: str | None = None,
    prediction=None,
) -> TsdSearchOptions:
    """Fill left/right from Auto detect / reference without a dialog."""
    from dataclasses import replace

    from repeatviewer.backends.tsd_search import (
        boundaries_from_sidecar_1based,
        preferred_tsd_boundaries_1based,
    )

    aln_len = max(1, int(aln_len))
    left, right = preferred_tsd_boundaries_1based(
        file_path,
        aln_len,
        boundary_source=opts.boundary_source,
        saved_left=opts.left_boundary,
        saved_right=opts.right_boundary,
        reference_extent_1based=reference_extent_1based,
        reference_sequence=reference_sequence,
        prediction=prediction,
    )
    if left <= 0 and right <= 0:
        rv = boundaries_from_sidecar_1based(file_path, aln_len)
        if rv is not None:
            left, right = rv
    return replace(opts, left_boundary=left, right_boundary=right)


def prompt_tsd_boundaries(
    *,
    opts: TsdSearchOptions,
    file_path: Path | str,
    aln_len: int,
    parent: QWidget | None = None,
    reference_extent_1based: tuple[int, int] | None = None,
    reference_sequence: str | None = None,
    on_dialog_ready=None,
    prediction=None,
) -> TsdSearchOptions | None:
    """Resolve TE edges for a TSD run.

    Auto detect / Reference: resolve silently (``prediction`` is the open
    MSA's Boundary prediction, used before ``*.boundaries_RV.json``).
    User define: show a non-modal dialog (MSA drag fills Start/End). Cancel →
    ``None``. Prefills from ``*.boundaries_RV.json`` when available.
    """
    from PySide6.QtCore import QEventLoop

    from repeatviewer.backends.tsd_search import boundaries_from_sidecar_1based

    if opts.boundary_source != TSD_BOUNDARY_USER:
        return _auto_resolve_tsd_boundaries(
            opts,
            file_path=file_path,
            aln_len=aln_len,
            reference_extent_1based=reference_extent_1based,
            reference_sequence=reference_sequence,
            prediction=prediction,
        )

    aln_len = max(1, int(aln_len))
    left, right = 0, 0
    from repeatviewer.backends.tsd_search import boundaries_from_prediction_1based

    rv = boundaries_from_prediction_1based(prediction) or boundaries_from_sidecar_1based(
        file_path, aln_len
    )
    if rv is not None:
        left, right = rv
    elif opts.left_boundary > 0 or opts.right_boundary > 0:
        left, right = opts.left_boundary, opts.right_boundary

    dialog = TsdBoundaryDialog(
        left_boundary=left,
        right_boundary=right,
        nt_radius=opts.nt_radius,
        aln_len=aln_len,
        parent=parent,
    )
    dialog.setModal(False)
    dialog.setWindowModality(Qt.WindowModality.NonModal)
    dialog.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, False)
    if callable(on_dialog_ready):
        on_dialog_ready(dialog)

    # QDialog.exec() forces modality even when setModal(False). Use show() +
    # a local event loop so the MSA stays interactive for drag-to-fill.
    loop = QEventLoop(dialog)
    finished_code = {"value": int(QDialog.DialogCode.Rejected)}

    def _on_finished(code: int) -> None:
        finished_code["value"] = int(code)
        loop.quit()

    dialog.finished.connect(_on_finished)
    dialog.show()
    dialog.raise_()
    dialog.activateWindow()
    loop.exec()
    if finished_code["value"] != int(QDialog.DialogCode.Accepted):
        return None
    return dialog.apply_to_options(opts)


_KEEP_CURRENT = object()

# Help text for BLAST dialog fields. BLAST+ option wording follows the NCBI
# BLAST Command Line Applications User Manual:
# https://www.ncbi.nlm.nih.gov/books/NBK1762/

_BLAST_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): BLAST uses the FASTA files that were "
    "selected in the Files panel when you click Run.<br><br>"
    "<b>Define inputs</b>: choose different FASTA files. Paths are for "
    "this run only and are not remembered next time."
)
_BLAST_DEFINE_INPUTS_HELP = (
    "Choose FASTA files for this BLAST run instead of the Files panel "
    "selection.<br><br>"
    "Click a file in the Files panel to fill a row, or use <b>Browse</b>. "
    "Paths are for this run only and are not remembered next time.<br><br>"
    "Each BLAST tab keeps its own Define-inputs list (BLASTn, BLASTp, BLASTx, "
    "and Online BLAST are independent)."
)
_BLAST_DATABASE_HELP = (
    "Subject database for the search.<br><br>"
    "BLASTn can use the project reference genome or a saved nucleotide BLAST "
    "database. BLASTp and BLASTx search a protein database chosen under "
    "<b>Manage Database</b>."
)
_BLAST_ONLINE_PROGRAM_HELP = (
    "NCBI BLAST suite:<br><br>"
    "<b>blastn</b> (DNA vs DNA)<br><br>"
    "<b>blastp</b> (protein vs protein)<br><br>"
    "<b>blastx</b> (translated DNA vs protein)<br><br>"
    "<b>tblastn</b> (protein vs translated DNA)<br><br>"
    "<b>tblastx</b> (translated DNA vs translated DNA)<br><br>"
    "Each query FASTA must contain <b>only one sequence</b>. "
    "Separate multi-sequence files first (right-click → Separate FASTA).<br><br>"
    "Results are written as an HTML report (<code>*_obn.html</code>, "
    "<code>*_obp.html</code>, <code>*_obx.html</code>, …) with "
    "Descriptions, Graphic Summary, and Lineage tabs."
)
_BLAST_ONLINE_WRITE_TEXT_HELP = (
    "Also keep the tabular hit table (<code>*_obn.txt</code>, …) next to the query. "
    "The HTML report is always written."
)
_BLAST_ONLINE_DATABASE_HELP = (
    "Public NCBI database searched on their servers (not a local RepeatViewer "
    "BLAST database). Choices follow NCBI web BLAST.<br><br>"
    "<b>Nucleotide (blastn / tblastn / tblastx)</b><br>"
    "<b>core_nt</b>: Core nucleotide: NCBI’s refined default; faster and more "
    "relevant than full <code>nt</code> (excludes some large eukaryotic "
    "chromosomes).<br>"
    "<b>nt</b>: Nucleotide collection (nr/nt): broad GenBank/EMBL/DDBJ set "
    "(excluding bulk divisions such as EST/GSS/STS/PAT/HTG/WGS).<br>"
    "<b>refseq_select_rna</b>: curated RefSeq Select transcripts (human/mouse "
    "and select prokaryotes); smaller and easier to interpret.<br>"
    "<b>refseq_rna</b>: NCBI RefSeq transcript (RNA) sequences.<br>"
    "<b>refseq_reference_genomes</b> / <b>refseq_genomes</b> / "
    "<b>refseq_genomic</b>: RefSeq genome assemblies at different scopes.<br>"
    "<b>wgs</b>: whole-genome shotgun contigs/scaffolds.<br>"
    "<b>est</b>: expressed sequence tags (cDNA fragments).<br>"
    "<b>tsa</b>: Transcriptome Shotgun Assembly.<br>"
    "<b>htgs</b>: unfinished high-throughput genomic sequences.<br>"
    "<b>gss</b>: genomic survey sequences.<br>"
    "<b>sts</b>: sequence tagged sites (markers).<br>"
    "<b>patnt</b>: patent nucleotide sequences.<br>"
    "<b>pdbnt</b>: nucleotide sequences from PDB structures.<br>"
    "<b>RefSeq_Gene</b>: human RefSeqGene loci.<br>"
    "<b>tls</b>: targeted loci (barcode / marker genes).<br>"
    "<b>sra</b>: Sequence Read Archive (very large; often slow or limited).<br><br>"
    "<b>Protein (blastp / blastx)</b><br>"
    "<b>nr</b>: non-redundant protein collection (broadest).<br>"
    "<b>refseq_select_prot</b>: RefSeq Select proteins.<br>"
    "<b>refseq_protein</b>: NCBI RefSeq proteins.<br>"
    "<b>swissprot</b>: UniProtKB/Swiss-Prot (highly curated).<br>"
    "<b>pdbaa</b>: PDB protein sequences.<br>"
    "<b>landmark</b>: compact model-organism proteomes.<br>"
    "<b>pataa</b>: patent protein sequences."
)
_BLAST_ONLINE_HITS_HELP = (
    "Maximum number of database hits to keep (NCBI <code>HITLIST_SIZE</code>)."
)
_BLAST_EVALUE_HELP = (
    "<b>-evalue</b>: Expectation value (E) threshold for saving hits.<br><br>"
    "Hits with an E-value above this threshold are discarded. Lower values keep "
    "fewer, more significant alignments. Typical default in BLAST+ is 10; "
    "RepeatViewer often uses a stricter value for TE work."
)
_BLAST_THREADS_HELP = (
    "<b>-num_threads</b>: Number of threads (CPUs) used to perform a search."
)
_BLAST_WORD_SIZE_HELP = (
    "<b>-word_size</b>: Length of the exact seed match (wordfinder) used to begin "
    "a search.<br><br>"
    "Smaller values are more sensitive but slower, especially on large or "
    "repetitive genomes. Typical BLASTn default is 11; protein searches "
    "(BLASTp / BLASTx) often use 3."
)
_BLAST_GAPOPEN_HELP = (
    "<b>-gapopen</b>: Cost to open a gap. Must be consistent with the scoring "
    "system used by the program."
)
_BLAST_GAPEXTEND_HELP = (
    "<b>-gapextend</b>: Cost to extend a gap."
)
_BLAST_PENALTY_HELP = (
    "<b>-penalty</b> (BLASTn): Penalty for a nucleotide mismatch. Must be a "
    "negative integer."
)
_BLAST_PERC_IDENTITY_HELP = (
    "<b>-perc_identity</b> (BLASTn): Percent identity cutoff.<br><br>"
    "Alignments below this percent identity are discarded. Use 0 to disable the "
    "filter."
)
_BLAST_QCOV_HELP = (
    "<b>-qcov_hsp_perc</b>: Percent query coverage per HSP.<br><br>"
    "Alignments that cover less than this percentage of the query length are "
    "discarded. Use 0 to disable the filter."
)
_BLAST_WRITE_TEXT_HELP = (
    "Keep the BLAST hits table (<code>*_bn.txt</code>, <code>*_bp.txt</code>, "
    "<code>*_bx.txt</code>) next to the query. "
    "When unchecked, the table is removed after BLAST."
)
_BLAST_MAX_SEQ_HELP = (
    "<b>Max sequence number</b><br><br>"
    "Maximum number of hit sequences included in the extracted FASTA used for "
    "multiple sequence alignment.<br><br>"
    "When a query has many BLAST hits, only this many sequences are kept so the "
    "downstream MSA stays manageable.<br><br>"
    "Default: <b>100</b>."
)
_BLAST_TOP_LONGEST_HELP = (
    "<b>Top longest sequences</b><br><br>"
    "Used with <b>Max sequence number</b>. If there are more BLAST hits than "
    "Max sequence number:<br><br>"
    "1. Sort hits by length and keep this many <b>longest</b> sequences.<br><br>"
    "2. Fill any remaining slots at random from the other hits until "
    "Max sequence number is reached.<br><br>"
    "Default: <b>100</b>."
)
class BlastParamsDialog(QDialog):
    # Same default outer width as AlignParamsDialog.
    _DIALOG_WIDTH = TOOL_PARAMS_DIALOG_WIDTH
    # Longest Settings label — keeps Inputs radios in the Settings field column.
    _SETTINGS_LABEL_REF = "Top longest sequences:"

    def __init__(
        self,
        *,
        project: RepeatViewerProject,
        initial: BlastOptions | None = None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("BLAST parameters")
        self._project = project
        self._for_run = for_run
        opts = initial or BlastOptions()
        self._initial_program = opts.program
        self.setMinimumWidth(self._DIALOG_WIDTH)

        # Always show the latest Manage-BLAST list (and current genome name).
        from repeatviewer.database_settings import reload_blast_databases_into_project

        reload_blast_databases_into_project(self._project)

        self._blastn_database = QComboBox()
        self._blastp_database = QComboBox()
        for combo in (self._blastn_database, self._blastp_database):
            combo.setFixedWidth(260)

        self._evalue = QDoubleSpinBox()
        configure_report_evalue_spinbox(self._evalue)
        self._evalue.setValue(normalize_blast_evalue(opts.e_value))
        self._evalue.setToolTip(_BLAST_EVALUE_HELP)

        self._threads = QSpinBox()
        self._threads.setRange(1, 128)
        self._threads.setValue(max(1, opts.num_threads))
        self._threads.setToolTip(_BLAST_THREADS_HELP)

        self._word_size = QSpinBox()
        self._word_size.setRange(4, 10_000)
        self._word_size.setValue(opts.word_size)
        self._word_size.setToolTip(_BLAST_WORD_SIZE_HELP)

        self._gap_open = QSpinBox()
        self._gap_open.setRange(0, 1000)
        self._gap_open.setValue(opts.gap_open)
        self._gap_open.setToolTip(_BLAST_GAPOPEN_HELP)

        self._gap_extend = QSpinBox()
        self._gap_extend.setRange(0, 1000)
        self._gap_extend.setValue(opts.gap_extend)
        self._gap_extend.setToolTip(_BLAST_GAPEXTEND_HELP)

        self._penalty = QSpinBox()
        self._penalty.setRange(-1000, 0)
        self._penalty.setValue(opts.penalty)
        self._penalty.setToolTip(_BLAST_PENALTY_HELP)

        self._percent_identity = QDoubleSpinBox()
        self._percent_identity.setRange(0.0, 100.0)
        self._percent_identity.setDecimals(1)
        self._percent_identity.setSuffix(" %")
        self._percent_identity.setValue(opts.percent_identity)
        self._percent_identity.setToolTip(_BLAST_PERC_IDENTITY_HELP)

        self._query_coverage = QDoubleSpinBox()
        self._query_coverage.setRange(0.0, 100.0)
        self._query_coverage.setDecimals(1)
        self._query_coverage.setSuffix(" %")
        self._query_coverage.setValue(opts.query_coverage_hsp)
        self._query_coverage.setToolTip(_BLAST_QCOV_HELP)

        self._blastn_write_text = QCheckBox("Export BLAST hits table")
        self._blastn_write_text.setChecked(opts.write_text)
        self._blastn_write_text.setToolTip(_BLAST_WRITE_TEXT_HELP)

        self._blastp_evalue = QDoubleSpinBox()
        configure_report_evalue_spinbox(self._blastp_evalue)
        self._blastp_evalue.setValue(normalize_blast_evalue(opts.blastp_e_value))
        self._blastp_evalue.setToolTip(_BLAST_EVALUE_HELP)

        self._blastp_threads = QSpinBox()
        self._blastp_threads.setRange(1, 128)
        self._blastp_threads.setValue(max(1, opts.blastp_num_threads))
        self._blastp_threads.setToolTip(_BLAST_THREADS_HELP)

        self._blastp_word_size = QSpinBox()
        self._blastp_word_size.setRange(2, 10_000)
        self._blastp_word_size.setValue(opts.blastp_word_size)
        self._blastp_word_size.setToolTip(_BLAST_WORD_SIZE_HELP)

        self._blastp_gap_open = QSpinBox()
        self._blastp_gap_open.setRange(0, 1000)
        self._blastp_gap_open.setValue(opts.blastp_gap_open)
        self._blastp_gap_open.setToolTip(_BLAST_GAPOPEN_HELP)

        self._blastp_gap_extend = QSpinBox()
        self._blastp_gap_extend.setRange(0, 1000)
        self._blastp_gap_extend.setValue(opts.blastp_gap_extend)
        self._blastp_gap_extend.setToolTip(_BLAST_GAPEXTEND_HELP)

        self._blastp_query_coverage = QDoubleSpinBox()
        self._blastp_query_coverage.setRange(0.0, 100.0)
        self._blastp_query_coverage.setDecimals(1)
        self._blastp_query_coverage.setSuffix(" %")
        self._blastp_query_coverage.setValue(opts.blastp_query_coverage_hsp)
        self._blastp_query_coverage.setToolTip(_BLAST_QCOV_HELP)
        self._blastp_write_text = QCheckBox("Export BLAST hits table")
        self._blastp_write_text.setChecked(opts.blastp_write_text)
        self._blastp_write_text.setToolTip(_BLAST_WRITE_TEXT_HELP)

        # --- blastx: nucleotide query translated against a protein database ---
        self._blastx_database = QComboBox()
        self._blastx_database.setFixedWidth(260)

        self._blastx_evalue = QDoubleSpinBox()
        configure_report_evalue_spinbox(self._blastx_evalue)
        self._blastx_evalue.setValue(normalize_blast_evalue(opts.blastx_e_value))
        self._blastx_evalue.setToolTip(_BLAST_EVALUE_HELP)

        self._blastx_threads = QSpinBox()
        self._blastx_threads.setRange(1, 128)
        self._blastx_threads.setValue(max(1, opts.blastx_num_threads))
        self._blastx_threads.setToolTip(_BLAST_THREADS_HELP)

        self._blastx_word_size = QSpinBox()
        self._blastx_word_size.setRange(2, 10_000)
        self._blastx_word_size.setValue(opts.blastx_word_size)
        self._blastx_word_size.setToolTip(_BLAST_WORD_SIZE_HELP)

        self._blastx_gap_open = QSpinBox()
        self._blastx_gap_open.setRange(0, 1000)
        self._blastx_gap_open.setValue(opts.blastx_gap_open)
        self._blastx_gap_open.setToolTip(_BLAST_GAPOPEN_HELP)

        self._blastx_gap_extend = QSpinBox()
        self._blastx_gap_extend.setRange(0, 1000)
        self._blastx_gap_extend.setValue(opts.blastx_gap_extend)
        self._blastx_gap_extend.setToolTip(_BLAST_GAPEXTEND_HELP)

        self._blastx_query_coverage = QDoubleSpinBox()
        self._blastx_query_coverage.setRange(0.0, 100.0)
        self._blastx_query_coverage.setDecimals(1)
        self._blastx_query_coverage.setSuffix(" %")
        self._blastx_query_coverage.setValue(opts.blastx_query_coverage_hsp)
        self._blastx_query_coverage.setToolTip(_BLAST_QCOV_HELP)
        self._blastx_write_text = QCheckBox("Export BLAST hits table")
        self._blastx_write_text.setChecked(opts.blastx_write_text)
        self._blastx_write_text.setToolTip(_BLAST_WRITE_TEXT_HELP)

        # --- Online BLAST (NCBI web) -----------------------------------------
        self._online_program = QComboBox()
        for label, value in (
            ("blastn", BLAST_PROGRAM_BLASTN),
            ("blastp", BLAST_PROGRAM_BLASTP),
            ("blastx", BLAST_PROGRAM_BLASTX),
            ("tblastn", BLAST_PROGRAM_TBLASTN),
            ("tblastx", BLAST_PROGRAM_TBLASTX),
        ):
            self._online_program.addItem(label, value)
        online_prog = (
            opts.online_program
            if opts.online_program in BLAST_ONLINE_ALGORITHMS
            else DEFAULT_BLAST_ONLINE_PROGRAM
        )
        idx = self._online_program.findData(online_prog)
        self._online_program.setCurrentIndex(idx if idx >= 0 else 0)
        self._online_program.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        self._online_program.setToolTip(_BLAST_ONLINE_PROGRAM_HELP)

        self._online_database = QComboBox()
        self._online_database.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        self._online_database.setToolTip(_BLAST_ONLINE_DATABASE_HELP)
        # Long NCBI labels: keep the field width, widen the popup list.
        self._online_database.view().setMinimumWidth(420)
        self._populate_online_database_combo(opts.online_database)

        self._online_evalue = QDoubleSpinBox()
        configure_report_evalue_spinbox(self._online_evalue)
        self._online_evalue.setValue(normalize_blast_evalue(opts.online_e_value))
        self._online_evalue.setToolTip(_BLAST_EVALUE_HELP)
        self._online_evalue.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)

        self._online_hitlist = QSpinBox()
        self._online_hitlist.setRange(1, 500)
        self._online_hitlist.setValue(max(1, min(500, int(opts.online_hitlist_size))))
        self._online_hitlist.setToolTip(_BLAST_ONLINE_HITS_HELP)
        self._online_hitlist.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)

        self._online_write_text = QCheckBox("Export BLAST hits table")
        self._online_write_text.setChecked(bool(opts.online_write_text))
        self._online_write_text.setToolTip(_BLAST_ONLINE_WRITE_TEXT_HELP)

        self._online_program.currentIndexChanged.connect(
            self._on_online_program_changed
        )

        # Shared post-processing limits: one spinbox set per program page, kept in sync.
        self._max_msa = QSpinBox()
        self._max_msa.setRange(10, 5000)
        self._max_msa.setValue(opts.max_msa_lines)
        self._max_msa.setToolTip(_BLAST_MAX_SEQ_HELP)
        self._top_msa = QSpinBox()
        self._top_msa.setRange(10, 5000)
        self._top_msa.setMaximum(self._max_msa.value())
        self._top_msa.setValue(opts.top_msa_lines)
        self._top_msa.setToolTip(_BLAST_TOP_LONGEST_HELP)

        self._blastp_max_msa = QSpinBox()
        self._blastp_max_msa.setRange(10, 5000)
        self._blastp_max_msa.setValue(opts.max_msa_lines)
        self._blastp_max_msa.setToolTip(_BLAST_MAX_SEQ_HELP)
        self._blastp_top_msa = QSpinBox()
        self._blastp_top_msa.setRange(10, 5000)
        self._blastp_top_msa.setMaximum(self._blastp_max_msa.value())
        self._blastp_top_msa.setValue(opts.top_msa_lines)
        self._blastp_top_msa.setToolTip(_BLAST_TOP_LONGEST_HELP)

        self._blastx_max_msa = QSpinBox()
        self._blastx_max_msa.setRange(10, 5000)
        self._blastx_max_msa.setValue(opts.max_msa_lines)
        self._blastx_max_msa.setToolTip(_BLAST_MAX_SEQ_HELP)
        self._blastx_top_msa = QSpinBox()
        self._blastx_top_msa.setRange(10, 5000)
        self._blastx_top_msa.setMaximum(self._blastx_max_msa.value())
        self._blastx_top_msa.setValue(opts.top_msa_lines)
        self._blastx_top_msa.setToolTip(_BLAST_TOP_LONGEST_HELP)

        self._max_msa_spins = (self._max_msa, self._blastp_max_msa, self._blastx_max_msa)
        self._top_msa_spins = (self._top_msa, self._blastp_top_msa, self._blastx_top_msa)
        for max_spin, top_spin in zip(self._max_msa_spins, self._top_msa_spins):
            max_spin.valueChanged.connect(top_spin.setMaximum)
        self._wire_synced_spins(self._max_msa_spins)
        self._wire_synced_spins(self._top_msa_spins)

        # Uniform field size so the inputs line up symmetrically.
        for field in (
            self._evalue,
            self._threads,
            self._word_size,
            self._gap_open,
            self._gap_extend,
            self._penalty,
            self._percent_identity,
            self._query_coverage,
            self._blastp_evalue,
            self._blastp_threads,
            self._blastp_word_size,
            self._blastp_gap_open,
            self._blastp_gap_extend,
            self._blastp_query_coverage,
            self._blastx_evalue,
            self._blastx_threads,
            self._blastx_word_size,
            self._blastx_gap_open,
            self._blastx_gap_extend,
            self._blastx_query_coverage,
            *self._max_msa_spins,
            *self._top_msa_spins,
        ):
            field.setFixedWidth(260)
            field.setAlignment(Qt.AlignmentFlag.AlignRight)

        # Shared Inputs radios; each program tab has its own Define-inputs list.
        self._input_group = QButtonGroup(self)
        self._input_selected = QRadioButton("Use selected sequences", self)
        self._input_define = QRadioButton("Define inputs", self)
        self._input_group.addButton(self._input_selected)
        self._input_group.addButton(self._input_define)
        self._input_selected.setChecked(True)

        def _make_blast_define_inputs() -> TeMultiInputFilesWidget:
            return TeMultiInputFilesWidget(
                self,
                initial=None,
                start_empty=False,
                placeholder="Click in Files/Folders or Browse…",
                browse_title="Choose a FASTA file",
                file_filter=(
                    "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
                ),
                add_button_label="Add another input",
                add_button_tooltip=(
                    "Add another FASTA file path. Click a file in the Files panel "
                    "to fill a row, or use Browse."
                ),
                path_mode="file",
            )

        self._blastn_define_inputs = _make_blast_define_inputs()
        self._blastp_define_inputs = _make_blast_define_inputs()
        self._blastx_define_inputs = _make_blast_define_inputs()
        self._online_define_inputs = _make_blast_define_inputs()
        self._define_inputs_by_tab = (
            self._blastn_define_inputs,
            self._blastp_define_inputs,
            self._blastx_define_inputs,
            self._online_define_inputs,
        )
        # Backward-compatible alias: active tab's widget (updated on tab change).
        self._define_inputs = self._blastn_define_inputs
        self._input_selected.toggled.connect(self._sync_blast_inputs)
        self._input_define.toggled.connect(self._sync_blast_inputs)
        for widget in self._define_inputs_by_tab:
            widget.paths_changed.connect(self._on_blast_define_paths_changed)

        self._inputs_panel = QWidget(self)
        self._inputs_panel.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )
        self._inputs_form = self._parameter_form(self._inputs_panel)
        self._inputs_form.setContentsMargins(16, 10, 16, 0)
        self._inputs_form.setVerticalSpacing(4)
        self._inputs_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        self._inputs_form.addRow(self._blast_section("Inputs", self._inputs_panel))
        self._inputs_radios_box = QWidget(self._inputs_panel)
        inputs_radios_layout = QVBoxLayout(self._inputs_radios_box)
        inputs_radios_layout.setContentsMargins(0, 0, 0, 0)
        inputs_radios_layout.setSpacing(2)
        inputs_radios_layout.addWidget(
            field_with_help(
                self._input_selected,
                "BLAST inputs",
                _BLAST_INPUTS_HELP,
                stretch=True,
            )
        )
        inputs_radios_layout.addWidget(
            field_with_help(
                self._input_define,
                "Define inputs",
                _BLAST_DEFINE_INPUTS_HELP,
                stretch=True,
            )
        )
        # Empty labels sized to Settings so radio dots line up with Database.
        self._inputs_radios_label = QLabel("")
        self._inputs_define_label = QLabel("")
        self._define_inputs_host = QWidget(self._inputs_panel)
        self._define_inputs_host.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )
        define_host_layout = QVBoxLayout(self._define_inputs_host)
        define_host_layout.setContentsMargins(0, 0, 0, 0)
        define_host_layout.setSpacing(0)
        for widget in self._define_inputs_by_tab:
            define_host_layout.addWidget(widget)
        self._inputs_form.addRow(self._inputs_radios_label, self._inputs_radios_box)
        self._inputs_form.addRow(self._inputs_define_label, self._define_inputs_host)
        self._inputs_form.addRow(_form_section_separator())
        self._sync_blast_inputs_field_column()
        self._sync_blast_inputs(resize=False)

        self._blast_compact_size: QSize | None = None

        self._use_system_blast = QCheckBox(system_path_executable_label("BLAST+"))
        self._use_system_blast.setChecked(project.use_system_blast)
        self._use_system_blast.setToolTip(
            "Use BLAST+ from PATH or the conda env in File → Settings → External tools "
            "instead of the bundled programs."
        )
        self._blast_path_status = make_executable_status_label(self)
        self._use_system_blast.toggled.connect(self._update_executable_status)

        (
            self._blastn_tab,
            self._blastn_basic,
            self._blastn_advanced,
        ) = self._build_program_tab(
            database_combo=self._blastn_database,
            evalue=self._evalue,
            threads=self._threads,
            word_size=self._word_size,
            max_msa=self._max_msa,
            top_msa=self._top_msa,
            write_text=self._blastn_write_text,
            advanced_rows=[
                (
                    "-gapopen:",
                    self._field_with_help(
                        self._gap_open, "-gapopen", _BLAST_GAPOPEN_HELP
                    ),
                ),
                (
                    "-gapextend:",
                    self._field_with_help(
                        self._gap_extend, "-gapextend", _BLAST_GAPEXTEND_HELP
                    ),
                ),
                (
                    "-penalty:",
                    self._field_with_help(
                        self._penalty, "-penalty", _BLAST_PENALTY_HELP
                    ),
                ),
                (
                    "-perc_identity:",
                    self._field_with_help(
                        self._percent_identity,
                        "-perc_identity",
                        _BLAST_PERC_IDENTITY_HELP,
                    ),
                ),
                (
                    "-qcov_hsp_perc:",
                    self._field_with_help(
                        self._query_coverage, "-qcov_hsp_perc", _BLAST_QCOV_HELP
                    ),
                ),
            ],
        )
        (
            self._blastp_tab,
            self._blastp_basic,
            self._blastp_advanced,
        ) = self._build_program_tab(
            database_combo=self._blastp_database,
            evalue=self._blastp_evalue,
            threads=self._blastp_threads,
            word_size=self._blastp_word_size,
            max_msa=self._blastp_max_msa,
            top_msa=self._blastp_top_msa,
            write_text=self._blastp_write_text,
            advanced_rows=[
                (
                    "-gapopen:",
                    self._field_with_help(
                        self._blastp_gap_open, "-gapopen", _BLAST_GAPOPEN_HELP
                    ),
                ),
                (
                    "-gapextend:",
                    self._field_with_help(
                        self._blastp_gap_extend, "-gapextend", _BLAST_GAPEXTEND_HELP
                    ),
                ),
                (
                    "-qcov_hsp_perc:",
                    self._field_with_help(
                        self._blastp_query_coverage,
                        "-qcov_hsp_perc",
                        _BLAST_QCOV_HELP,
                    ),
                ),
            ],
        )
        (
            self._blastx_tab,
            self._blastx_basic,
            self._blastx_advanced,
        ) = self._build_program_tab(
            database_combo=self._blastx_database,
            evalue=self._blastx_evalue,
            threads=self._blastx_threads,
            word_size=self._blastx_word_size,
            max_msa=self._blastx_max_msa,
            top_msa=self._blastx_top_msa,
            write_text=self._blastx_write_text,
            advanced_rows=[
                (
                    "-gapopen:",
                    self._field_with_help(
                        self._blastx_gap_open, "-gapopen", _BLAST_GAPOPEN_HELP
                    ),
                ),
                (
                    "-gapextend:",
                    self._field_with_help(
                        self._blastx_gap_extend, "-gapextend", _BLAST_GAPEXTEND_HELP
                    ),
                ),
                (
                    "-qcov_hsp_perc:",
                    self._field_with_help(
                        self._blastx_query_coverage,
                        "-qcov_hsp_perc",
                        _BLAST_QCOV_HELP,
                    ),
                ),
            ],
        )

        (
            self._online_tab,
            self._online_basic,
            self._online_advanced,
        ) = self._build_online_tab()

        # --- Options: right-click program + shared executable setting ---
        options_tab = QWidget()
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)

        rc_box = QGroupBox("Right-click action", self)
        rc_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_label = QLabel(
            "Right-clicking the BLAST button runs, with saved settings:", self
        )
        rc_label.setWordWrap(True)
        rc_box_layout.addWidget(rc_label)
        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("BLASTn", BLAST_PROGRAM_BLASTN),
            ("BLASTp", BLAST_PROGRAM_BLASTP),
            ("BLASTx", BLAST_PROGRAM_BLASTX),
            ("Online BLAST", BLAST_PROGRAM_ONLINE),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_program == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
            rc_box_layout.addWidget(radio)
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[BLAST_PROGRAM_BLASTN].setChecked(True)
        options_layout.addWidget(rc_box)

        exe_box = QGroupBox("Executables", self)
        exe_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        exe_box_layout = QVBoxLayout(exe_box)
        exe_box_layout.setContentsMargins(10, 8, 10, 8)
        exe_box_layout.setSpacing(8)
        exe_box_layout.addWidget(
            executable_checkbox_row(
                self, self._use_system_blast, self._blast_path_status
            )
        )
        options_layout.addWidget(exe_box)

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="BLAST",
                checked=bool(getattr(opts, "show_process_window", True)),
                close_when_finished=bool(
                    getattr(opts, "close_process_window_when_finished", True)
                ),
            )
        )
        options_layout.addWidget(process_box)
        from repeatviewer.tool_citations import CITATION_BLAST

        options_layout.addWidget(make_citation_box(self, CITATION_BLAST))
        # Absorb leftover height so group boxes stay tight (Align Options pattern).
        options_layout.addStretch(1)
        options_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        options_tab.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._options_tab = options_tab

        self._tabs = QTabBar(self)
        self._tabs.setExpanding(False)
        self._tabs.setDocumentMode(True)
        self._tabs.setDrawBase(False)
        self._tabs.addTab("BLASTn")
        self._tabs.addTab("BLASTp")
        self._tabs.addTab("BLASTx")
        self._tabs.addTab("Online BLAST")
        self._tabs.addTab("Options")
        self._options_tab_index = 4

        self._stack = QStackedWidget(self)
        self._stack.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._stack.addWidget(self._blastn_tab)
        self._stack.addWidget(self._blastp_tab)
        self._stack.addWidget(self._blastx_tab)
        self._stack.addWidget(self._online_tab)
        self._stack.addWidget(self._options_tab)

        self._populate_database_combos(
            blastn_name=opts.blast_database_name,
            blastp_name=opts.blastp_database_name,
            blastx_name=opts.blastx_database_name,
        )
        open_tab = opening_dialog_tab(
            getattr(opts, "last_dialog_tab", None),
            allowed=BLAST_DIALOG_TABS,
            first_tab=(
                opts.program
                if opts.program in BLAST_PROGRAMS
                else BLAST_PROGRAM_BLASTN
            ),
        )
        program_tab = {
            BLAST_PROGRAM_BLASTN: 0,
            BLAST_PROGRAM_BLASTP: 1,
            BLAST_PROGRAM_BLASTX: 2,
            BLAST_PROGRAM_ONLINE: 3,
        }.get(open_tab, 0)
        self._tabs.setCurrentIndex(program_tab)
        self._stack.setCurrentIndex(program_tab)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._advanced_settings_btn = QPushButton("Advanced Settings ▸", self)
        self._advanced_settings_btn.setCheckable(True)
        self._advanced_settings_btn.setChecked(False)
        self._advanced_settings_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._advanced_settings_btn.toggled.connect(self._on_advanced_settings_toggled)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore, self._advanced_settings_btn],
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        # Match Align (MAFFT→Mode): no extra gap under the tab underline; the
        # Inputs form's top margin (10) mirrors mafft_root's content padding.
        layout.setSpacing(0)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._tabs, 0)
        layout.addWidget(self._inputs_panel, 0)
        layout.addWidget(self._stack, 0)
        layout.addSpacing(8)
        layout.addWidget(self._button_box, 0)
        self._tabs.currentChanged.connect(self._on_blast_tab_changed)
        self._on_blast_tab_changed(self._tabs.currentIndex())
        self.resize(self._DIALOG_WIDTH, 420)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        # Refresh in case Manage / Add ran while this dialog was constructed,
        # and pick up a newly defined genome for the reference-genome row.
        self._refresh_database_combos()
        # Force a tight compact height after polish (no blank gap under Inputs).
        self._blast_compact_size = None
        self._sync_blast_inputs_field_column()
        if self._tabs.currentIndex() == self._options_tab_index:
            self._update_executable_status()
        QTimer.singleShot(0, lambda: self._resize_for_blast_advanced(False))

    def _refresh_database_combos(
        self,
        *,
        blastn_name=_KEEP_CURRENT,
        blastp_name=_KEEP_CURRENT,
        blastx_name=_KEEP_CURRENT,
    ) -> None:
        from repeatviewer.database_settings import reload_blast_databases_into_project

        reload_blast_databases_into_project(self._project)
        self._populate_database_combos(
            blastn_name=blastn_name,
            blastp_name=blastp_name,
            blastx_name=blastx_name,
        )

    def should_run(self) -> bool:
        return dialog_should_run(self)

    @staticmethod
    def _wire_synced_spins(spins: tuple) -> None:
        def _sync(source) -> None:
            value = source.value()

            def _apply() -> None:
                for spin in spins:
                    if spin is source:
                        continue
                    if spin.value() == value:
                        continue
                    spin.blockSignals(True)
                    spin.setValue(value)
                    spin.blockSignals(False)

            _apply()

        for spin in spins:
            spin.valueChanged.connect(lambda _v, s=spin: _sync(s))

    def _blast_section(self, title: str, parent: QWidget | None = None) -> QLabel:
        return QLabel(f"<b>{title}</b>", parent or self)

    def _build_program_tab(
        self,
        *,
        database_combo: QComboBox,
        evalue: QDoubleSpinBox,
        threads: QSpinBox,
        word_size: QSpinBox,
        max_msa: QSpinBox,
        top_msa: QSpinBox,
        write_text: QCheckBox,
        advanced_rows: list[tuple[str, QWidget]],
    ) -> tuple[QWidget, QWidget, QWidget]:
        tab = QWidget()
        tab.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        root = QVBoxLayout(tab)
        root.setContentsMargins(16, 0, 16, 10)
        root.setSpacing(4)
        basic = QWidget(tab)
        basic.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        form = self._parameter_form(basic)
        form.setContentsMargins(0, 0, 0, 0)
        form.setVerticalSpacing(4)

        form.addRow(self._blast_section("Settings", basic))
        form.addRow("Database:", self._database_field(database_combo))
        form.addRow(
            "-evalue:",
            self._field_with_help(evalue, "-evalue", _BLAST_EVALUE_HELP),
        )
        form.addRow(
            "Threads:",
            self._field_with_help(threads, "Threads", _BLAST_THREADS_HELP),
        )
        form.addRow(
            "-word_size:",
            self._field_with_help(word_size, "-word_size", _BLAST_WORD_SIZE_HELP),
        )
        form.addRow(
            "Max sequence number:",
            self._field_with_help(
                max_msa, "Max sequence number", _BLAST_MAX_SEQ_HELP
            ),
        )
        form.addRow(
            "Top longest sequences:",
            self._field_with_help(
                top_msa, "Top longest sequences", _BLAST_TOP_LONGEST_HELP
            ),
        )
        form.addRow(
            "",
            self._field_with_help(
                write_text, "Export BLAST hits table", _BLAST_WRITE_TEXT_HELP
            ),
        )
        root.addWidget(basic)

        advanced = QWidget(tab)
        advanced.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        advanced_form = QFormLayout(advanced)
        advanced_form.setContentsMargins(0, 4, 0, 0)
        advanced_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        advanced_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        advanced_form.setHorizontalSpacing(12)
        advanced_form.setVerticalSpacing(4)
        advanced_form.addRow(_form_section_separator())
        advanced_form.addRow(self._blast_section("Advanced Settings", advanced))
        for label, widget in advanced_rows:
            advanced_form.addRow(label, widget)
        root.addWidget(advanced)
        advanced.setVisible(False)
        advanced.setMaximumHeight(0)
        return (tab, basic, advanced)

    def _build_online_tab(self) -> tuple[QWidget, QWidget, QWidget]:
        """Online BLAST Settings page (NCBI web); no Advanced panel."""
        tab = QWidget()
        tab.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        root = QVBoxLayout(tab)
        root.setContentsMargins(16, 0, 16, 10)
        root.setSpacing(4)
        basic = QWidget(tab)
        basic.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        form = self._parameter_form(basic)
        form.setContentsMargins(0, 0, 0, 0)
        form.setVerticalSpacing(4)

        form.addRow(self._blast_section("Settings", basic))
        form.addRow(
            "Program:",
            self._field_with_help(
                self._online_program, "Program", _BLAST_ONLINE_PROGRAM_HELP
            ),
        )
        form.addRow(
            "Database:",
            self._field_with_help(
                self._online_database, "Database", _BLAST_ONLINE_DATABASE_HELP
            ),
        )
        form.addRow(
            "-evalue:",
            self._field_with_help(
                self._online_evalue, "-evalue", _BLAST_EVALUE_HELP
            ),
        )
        form.addRow(
            "Max hits:",
            self._field_with_help(
                self._online_hitlist, "Max hits", _BLAST_ONLINE_HITS_HELP
            ),
        )
        form.addRow(
            "",
            self._field_with_help(
                self._online_write_text,
                "Export BLAST hits table",
                _BLAST_ONLINE_WRITE_TEXT_HELP,
            ),
        )
        root.addWidget(basic)
        # Placeholder advanced widget (never shown) so tab shape matches others.
        advanced = QWidget(tab)
        advanced.setVisible(False)
        advanced.setMaximumHeight(0)
        root.addWidget(advanced)
        return (tab, basic, advanced)

    def _populate_online_database_combo(self, selected: str | None = None) -> None:
        from repeatviewer.backends.ncbi_blast_online import (
            ONLINE_BLAST_DATABASES,
            default_online_database,
        )

        program = str(self._online_program.currentData() or DEFAULT_BLAST_ONLINE_PROGRAM)
        choices = ONLINE_BLAST_DATABASES.get(
            program, ONLINE_BLAST_DATABASES[BLAST_PROGRAM_BLASTN]
        )
        wanted = str(selected or self._online_database.currentData() or "").strip()
        if not wanted or wanted not in {key for key, _ in choices}:
            wanted = default_online_database(program)
        self._online_database.blockSignals(True)
        self._online_database.clear()
        for key, label in choices:
            self._online_database.addItem(label, key)
        idx = self._online_database.findData(wanted)
        self._online_database.setCurrentIndex(idx if idx >= 0 else 0)
        self._online_database.blockSignals(False)

    def _on_online_program_changed(self, *_args) -> None:
        self._populate_online_database_combo()

    def input_source(self) -> str:
        if self._input_define.isChecked():
            return ALIGN_INPUT_DEFINE
        return ALIGN_INPUT_SELECTED

    def _active_define_inputs(self) -> TeMultiInputFilesWidget:
        index = self._tabs.currentIndex() if hasattr(self, "_tabs") else 0
        if 0 <= index < len(self._define_inputs_by_tab):
            return self._define_inputs_by_tab[index]
        return self._blastn_define_inputs

    def _on_blast_define_paths_changed(self) -> None:
        # Only the active tab's widget should drive a height refit.
        if self.sender() is not None and self.sender() is not self._active_define_inputs():
            return
        self._sync_blast_inputs()

    def run_files(self, *, panel_files: list[Path] | None = None) -> list[Path]:
        from repeatviewer.backends.fastga_dotplot import expand_fastga_input_paths

        if self.input_source() == ALIGN_INPUT_DEFINE:
            return expand_fastga_input_paths(self._active_define_inputs().paths())
        return [Path(p) for p in (panel_files or [])]

    def _sync_blast_inputs(self, *, resize: bool = True) -> None:
        define = self.input_source() == ALIGN_INPUT_DEFINE
        active = self._active_define_inputs()
        self._define_inputs = active
        # Collapse the form row fully — a hidden host still reports a large
        # sizeHint and leaves a blank gap above Settings.
        if hasattr(self, "_inputs_form"):
            self._inputs_form.setRowVisible(self._define_inputs_host, define)
        self._define_inputs_host.setVisible(define)
        for widget in self._define_inputs_by_tab:
            show = define and widget is active
            widget.setVisible(show)
            if show:
                widget.setMinimumHeight(0)
                widget.setMaximumHeight(16777215)
            else:
                widget.setFixedHeight(0)
            widget.updateGeometry()
        if define:
            self._define_inputs_host.setMinimumHeight(0)
            self._define_inputs_host.setMaximumHeight(16777215)
        else:
            self._define_inputs_host.setFixedHeight(0)
        self._define_inputs_host.updateGeometry()
        self._inputs_panel.updateGeometry()
        self._sync_blast_inputs_field_column()
        # Keep dialog height in sync when Define inputs opens/closes (Align-style).
        if (
            resize
            and self.isVisible()
            and hasattr(self, "_advanced_settings_btn")
            and not self._advanced_settings_btn.isChecked()
        ):
            self._blast_compact_size = None
            QTimer.singleShot(0, lambda: self._resize_for_blast_advanced(False))

    def _sync_blast_inputs_field_column(self) -> None:
        """Indent Inputs radios to the Settings field column (Align-style)."""
        if not hasattr(self, "_inputs_radios_label"):
            return
        basic = self._current_blast_basic() if hasattr(self, "_blastn_basic") else None
        form = basic.layout() if basic is not None else None
        sync_form_field_column_spacers(
            self,
            [self._inputs_radios_label, self._inputs_define_label],
            [form if isinstance(form, QFormLayout) else None],
            ref_label=self._SETTINGS_LABEL_REF,
        )

    def use_query_file(self, path: Path) -> bool:
        """Fill Define inputs with ``path`` (sequences sent from the contig view), on the
        BLASTn tab unless a nucleotide-query tab is already open."""
        self._input_define.setChecked(True)
        if self._tabs.currentIndex() in (1, self._options_tab_index):   # BLASTp / Options
            self._tabs.setCurrentIndex(0)
        return self.apply_path_from_panel(path)

    def apply_path_from_panel(self, path: Path) -> bool:
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        if self._tabs.currentIndex() == self._options_tab_index:
            return False
        try:
            candidate = Path(path).expanduser().resolve()
        except OSError:
            return False
        return bool(self._active_define_inputs().accept_path(candidate))

    def apply_folder_from_panel(self, folder: Path) -> bool:
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        if self._tabs.currentIndex() == self._options_tab_index:
            return False
        try:
            candidate = Path(folder).expanduser().resolve()
        except OSError:
            return False
        if not candidate.is_dir():
            return False
        return bool(self._active_define_inputs().accept_path(candidate))

    def _current_advanced_panel(self) -> QWidget | None:
        index = self._tabs.currentIndex()
        if index == 0:
            return self._blastn_advanced
        if index == 1:
            return self._blastp_advanced
        if index == 2:
            return self._blastx_advanced
        # Online BLAST has no Advanced Settings panel.
        return None

    def _on_blast_tab_changed(self, tab_index: int) -> None:
        """Switch program tabs like Align: no window flash, resize only when needed."""
        on_program = tab_index in (0, 1, 2, 3)
        self._stack.setCurrentIndex(tab_index)
        self._inputs_panel.setVisible(on_program)
        has_advanced = self._current_advanced_panel() is not None
        self._advanced_settings_btn.setVisible(on_program and has_advanced)
        if not on_program:
            self._update_executable_status()
            if self.isVisible():
                self._fit_to_content()
            return

        advanced_was_open = self._advanced_settings_btn.isChecked()
        if advanced_was_open:
            self._advanced_settings_btn.blockSignals(True)
            self._advanced_settings_btn.setChecked(False)
            self._advanced_settings_btn.setText("Advanced Settings ▸")
            self._advanced_settings_btn.blockSignals(False)
        for panel in (
            self._blastn_advanced,
            self._blastp_advanced,
            self._blastx_advanced,
        ):
            if panel.isVisible() or panel.maximumHeight() != 0:
                panel.setVisible(False)
                panel.setMaximumHeight(0)

        # Show this tab's Define-inputs list (independent of other BLAST tabs).
        self._sync_blast_inputs(resize=False)

        # Align keeps each tab's own size; always re-fit when switching back
        # from Options (or after Advanced) so height/pin stay consistent.
        if self.isVisible():
            self._blast_compact_size = None
            self._resize_for_blast_advanced(False)

    def _set_blast_advanced_visible(self, panel: QWidget, visible: bool) -> None:
        panel.setVisible(visible)
        if visible:
            panel.setMaximumHeight(16777215)
        else:
            panel.setMaximumHeight(0)
        panel.updateGeometry()
        page = self._stack.currentWidget()
        if page is not None:
            root = page.layout()
            if root is not None:
                root.invalidate()
                root.activate()
            page.updateGeometry()
        self._stack.updateGeometry()

    def _on_advanced_settings_toggled(self, checked: bool) -> None:
        self._advanced_settings_btn.setText(
            "Advanced Settings ▾" if checked else "Advanced Settings ▸"
        )
        panel = self._current_advanced_panel()
        if panel is None:
            return
        if checked:
            # Remember tight height before expanding (Align-style).
            self._blast_compact_size = self.size()
        self._set_blast_advanced_visible(panel, checked)
        # Resize immediately — no QTimer, avoids window flash.
        self._resize_for_blast_advanced(checked)

    def _blast_dialog_chrome_height(self) -> int:
        lay = self.layout()
        margins = lay.contentsMargins() if lay is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 34
        spacing = lay.spacing() if lay is not None else 8
        # tab bar + (inputs) + stack + buttons
        n_gaps = 3 if self._inputs_panel.isVisible() else 2
        tab_bar_h = self._tabs.sizeHint().height() + 4
        buttons_h = self._button_box.sizeHint().height()
        return margin_h + (spacing * n_gaps) + tab_bar_h + buttons_h

    def _current_blast_basic(self) -> QWidget | None:
        index = self._tabs.currentIndex()
        if index == 0:
            return self._blastn_basic
        if index == 1:
            return self._blastp_basic
        if index == 2:
            return self._blastx_basic
        if index == 3:
            return self._online_basic
        return None

    def _blast_stack_page_height(self, *, expanded: bool) -> int:
        """Height of the active stack page (Settings ± Advanced), not Inputs."""
        page = self._stack.currentWidget()
        basic = self._current_blast_basic()
        if page is None or basic is None:
            return 220
        root = page.layout()
        margins = root.contentsMargins() if root is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 10
        spacing = root.spacing() if root is not None else 4
        basic.ensurePolished()
        basic_h = max(basic.minimumSizeHint().height(), basic.sizeHint().height())
        # Keep the layout spacing before Advanced even when it is collapsed
        # (Align-style), so the pinned stack height matches the real layout.
        if not expanded:
            return margin_h + basic_h + spacing
        panel = self._current_advanced_panel()
        extra = max(0, panel.sizeHint().height()) if panel is not None else 0
        return margin_h + basic_h + spacing + extra

    def _blast_inputs_height(self) -> int:
        if not self._inputs_panel.isVisible():
            return 0
        self._inputs_panel.ensurePolished()
        return max(
            self._inputs_panel.minimumSizeHint().height(),
            self._inputs_panel.sizeHint().height(),
        )

    def _pin_stack_height(self, page_h: int | None) -> None:
        """Constrain the stack to the active page so max-page sizeHint cannot gap."""
        page = self._stack.currentWidget()
        if page_h is None:
            self._stack.setMinimumHeight(0)
            self._stack.setMaximumHeight(16777215)
            if page is not None:
                page.setMinimumHeight(0)
                page.setMaximumHeight(16777215)
            return
        self._stack.setMinimumHeight(0)
        self._stack.setMaximumHeight(page_h)
        if page is not None:
            page.setMinimumHeight(0)
            page.setMaximumHeight(page_h)

    def _resize_for_blast_advanced(self, expanded: bool) -> None:
        if self._tabs.currentIndex() == self._options_tab_index:
            self._pin_stack_height(None)
            self._fit_to_content()
            return
        page = self._stack.currentWidget()
        if page is not None:
            root = page.layout()
            if root is not None:
                root.invalidate()
                root.activate()
            page.updateGeometry()
        self._stack.updateGeometry()
        self._inputs_panel.updateGeometry()
        self._sync_blast_inputs_field_column()

        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        chrome = self._blast_dialog_chrome_height()
        inputs_h = self._blast_inputs_height()
        page_h = self._blast_stack_page_height(expanded=expanded)
        width = tool_params_dialog_width(
            define_inputs=self.input_source() == ALIGN_INPUT_DEFINE
        )
        self.setMinimumWidth(width)

        if not expanded:
            height = chrome + inputs_h + page_h
            height = min(max(height, 220), max_h)
            self.setMinimumHeight(0)
            self._pin_stack_height(page_h)
            self.setMaximumHeight(height)
            self.resize(width, height)
            self.setMaximumHeight(16777215)
            self._blast_compact_size = QSize(width, height)
            return

        self._pin_stack_height(page_h)
        height = chrome + inputs_h + page_h
        height = min(max(height, 220), max_h)
        if self._blast_compact_size is not None:
            # Grow from the remembered compact size when Advanced opens.
            height = max(height, self._blast_compact_size.height())
            height = min(height, max_h)
        self.resize(width, height)

    def _fit_to_blast_content(self) -> None:
        if self._tabs.currentIndex() == self._options_tab_index:
            self._fit_to_content()
            return
        self._resize_for_blast_advanced(self._advanced_settings_btn.isChecked())

    @staticmethod
    def _parameter_form(parent: QWidget) -> QFormLayout:
        form = QFormLayout(parent)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        return form

    def _database_field(self, combo: QComboBox) -> QWidget:
        """Database combo with Manage BLAST databases on the next line."""
        manage = QPushButton("Manage Database")
        manage.clicked.connect(self._manage_databases)
        box = QWidget()
        layout = QVBoxLayout(box)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        layout.addWidget(
            self._field_with_help(combo, "Database", _BLAST_DATABASE_HELP)
        )
        layout.addWidget(manage, 0, Qt.AlignmentFlag.AlignLeft)
        return box

    def _field_with_help(self, field: QWidget, title: str, help_text: str) -> QWidget:
        return field_with_help(field, title, help_text)

    def _fit_to_content(self) -> None:
        """Size width/height to the active Options page (not the widest stack page)."""
        page = self._stack.currentWidget()
        if page is None or self.layout() is None:
            return
        self._pin_stack_height(None)
        page.ensurePolished()
        page.adjustSize()
        page_hint = page.minimumSizeHint().expandedTo(page.sizeHint())
        # Ignore stretch-inflated hints: use content height from layout.
        root = page.layout()
        content_h = page_hint.height()
        if root is not None:
            root.activate()
            content_h = max(
                content_h,
                root.totalMinimumSize().height(),
                root.totalSizeHint().height(),
            )
        tab_chrome = max(24, self._tabs.sizeHint().height()) + 8
        buttons_h = self._button_box.sizeHint().height()
        buttons_w = self._button_box.sizeHint().width()
        margins = self.layout().contentsMargins()
        spacing = max(0, self.layout().spacing())
        width = max(
            self._DIALOG_WIDTH,
            page_hint.width() + margins.left() + margins.right() + 24,
            buttons_w + margins.left() + margins.right(),
        )
        height = (
            content_h
            + tab_chrome
            + buttons_h
            + margins.top()
            + margins.bottom()
            + spacing * 2
            + 4
        )
        self._stack.setMinimumHeight(0)
        self._stack.setMaximumHeight(max(1, content_h))
        page.setMinimumHeight(0)
        page.setMaximumHeight(max(1, content_h))
        self.setMinimumHeight(0)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def _populate_database_combo(
        self,
        combo: QComboBox,
        molecule_type: str,
        *,
        selected_name=_KEEP_CURRENT,
    ) -> None:
        target = combo.currentData() if selected_name is _KEEP_CURRENT else selected_name
        combo.clear()
        if molecule_type == BLAST_MOLECULE_NUCL:
            genome = (
                self._project.genome_file
                if self._project.has_genome_configured()
                else None
            )
            combo.addItem(reference_genome_combo_label(self._project, genome), None)
        else:
            combo.addItem("Select a protein database…", None)
        model = combo.model()
        for db in self._project.blast_databases:
            if db.molecule_type != molecule_type:
                continue
            ready = blast_database_is_ready(db)
            label = db.name if ready else f"{db.name} (unavailable — index files missing)"
            combo.addItem(label, db.name)
            if not ready and isinstance(model, QStandardItemModel):
                item = model.item(combo.count() - 1)
                if item is not None:
                    item.setEnabled(False)
        index = combo.findData(target) if target is not None else 0
        combo.setCurrentIndex(index if index >= 0 else 0)

    def _populate_database_combos(
        self,
        *,
        blastn_name=_KEEP_CURRENT,
        blastp_name=_KEEP_CURRENT,
        blastx_name=_KEEP_CURRENT,
    ) -> None:
        self._populate_database_combo(
            self._blastn_database,
            BLAST_MOLECULE_NUCL,
            selected_name=blastn_name,
        )
        self._populate_database_combo(
            self._blastp_database,
            BLAST_MOLECULE_PROT,
            selected_name=blastp_name,
        )
        # blastx also searches a protein database.
        self._populate_database_combo(
            self._blastx_database,
            BLAST_MOLECULE_PROT,
            selected_name=blastx_name,
        )

    def _manage_databases(self) -> None:
        from repeatviewer.views.blast_database_dialog import BlastDatabaseManagerDialog

        tab = self._tabs.currentIndex() if hasattr(self, "_tabs") else 0
        if tab == 0:
            compatible = "blast_nucl"
            role = "blast"
        elif tab in (1, 2):
            compatible = "blast_prot"
            role = "blast"
        else:
            # Online BLAST / Options — allow managing DBs, but Select needs a program tab.
            compatible = None
            role = "blast"

        # Non-modal manager (and Add dialog) so Files-panel clicks can fill FASTA.
        existing = getattr(self, "_blast_manager_dialog", None)
        if existing is not None and existing.isVisible():
            if compatible is not None:
                existing.set_pick_mode(
                    True,
                    callback=self._on_blast_database_picked,
                    role_filter=role,
                    compatible=compatible,
                )
            existing.raise_()
            existing.activateWindow()
            return
        dialog = BlastDatabaseManagerDialog(self._project, parent=self.window())
        self._blast_manager_dialog = dialog
        dialog.references_changed.connect(self._refresh_database_combos)
        dialog.finished.connect(self._on_blast_manager_finished)
        if compatible is not None:
            dialog.set_pick_mode(
                True,
                callback=self._on_blast_database_picked,
                role_filter=role,
                compatible=compatible,
            )
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _on_blast_database_picked(self, payload) -> None:
        from repeatviewer.views.manage_references_dialog import (
            database_name_from_pick,
            fasta_path_from_database_pick,
        )
        from repeatviewer.project import BlastDatabase

        tab = self._tabs.currentIndex() if hasattr(self, "_tabs") else 0
        if tab == 0:
            combo = self._blastn_database
            molecule = BLAST_MOLECULE_NUCL
        elif tab == 1:
            combo = self._blastp_database
            molecule = BLAST_MOLECULE_PROT
        elif tab == 2:
            combo = self._blastx_database
            molecule = BLAST_MOLECULE_PROT
        else:
            QMessageBox.information(
                self,
                "BLAST",
                "Switch to the BLASTn, BLASTp, or BLASTx tab, then Select a database.",
            )
            return

        name = database_name_from_pick(payload)
        if isinstance(payload, BlastDatabase):
            if payload.molecule_type != molecule:
                need = "nucleotide" if molecule == BLAST_MOLECULE_NUCL else "protein"
                QMessageBox.information(
                    self,
                    "BLAST",
                    f"This tab needs a {need} BLAST database.",
                )
                return
            name = str(payload.name or "").strip() or name
        elif tab == 0:
            # Reference-genome style path pick for BLASTn.
            path = fasta_path_from_database_pick(payload)
            if path is not None and self._project.has_genome_configured():
                genome = self._project.genome_file
                if genome is not None and path.resolve() == Path(genome).resolve():
                    self._populate_database_combo(
                        combo, molecule, selected_name=None
                    )
                    return

        if not name:
            QMessageBox.information(
                self,
                "BLAST",
                "Select a compatible BLAST database.",
            )
            return
        self._populate_database_combo(combo, molecule, selected_name=name)

    def _on_blast_manager_finished(self, _result: int = 0) -> None:
        self._blast_manager_dialog = None
        self._refresh_database_combos()

    def _restore_defaults(self) -> None:
        defaults = BlastOptions()
        self._evalue.setValue(defaults.e_value)
        self._threads.setValue(defaults.num_threads)
        self._word_size.setValue(defaults.word_size)
        self._gap_open.setValue(defaults.gap_open)
        self._gap_extend.setValue(defaults.gap_extend)
        self._penalty.setValue(defaults.penalty)
        self._percent_identity.setValue(defaults.percent_identity)
        self._query_coverage.setValue(defaults.query_coverage_hsp)
        self._blastn_write_text.setChecked(defaults.write_text)
        self._blastp_evalue.setValue(defaults.blastp_e_value)
        self._blastp_threads.setValue(defaults.blastp_num_threads)
        self._blastp_word_size.setValue(defaults.blastp_word_size)
        self._blastp_gap_open.setValue(defaults.blastp_gap_open)
        self._blastp_gap_extend.setValue(defaults.blastp_gap_extend)
        self._blastp_query_coverage.setValue(defaults.blastp_query_coverage_hsp)
        self._blastp_write_text.setChecked(defaults.blastp_write_text)
        self._blastx_evalue.setValue(defaults.blastx_e_value)
        self._blastx_threads.setValue(defaults.blastx_num_threads)
        self._blastx_word_size.setValue(defaults.blastx_word_size)
        self._blastx_gap_open.setValue(defaults.blastx_gap_open)
        self._blastx_gap_extend.setValue(defaults.blastx_gap_extend)
        self._blastx_query_coverage.setValue(defaults.blastx_query_coverage_hsp)
        self._blastx_write_text.setChecked(defaults.blastx_write_text)
        online_idx = self._online_program.findData(defaults.online_program)
        self._online_program.setCurrentIndex(online_idx if online_idx >= 0 else 0)
        self._populate_online_database_combo(defaults.online_database)
        self._online_evalue.setValue(defaults.online_e_value)
        self._online_hitlist.setValue(defaults.online_hitlist_size)
        self._online_write_text.setChecked(defaults.online_write_text)
        # Canonical copies drive synced values on the other program pages.
        self._max_msa.setValue(defaults.max_msa_lines)
        self._top_msa.setValue(defaults.top_msa_lines)
        self._blastn_database.setCurrentIndex(0)
        self._blastp_database.setCurrentIndex(0)
        self._blastx_database.setCurrentIndex(0)
        self._input_selected.setChecked(True)
        for widget in self._define_inputs_by_tab:
            widget.clear_paths()
        self._sync_blast_inputs()
        self._use_system_blast.setChecked(False)
        self._show_process_window.setChecked(defaults.show_process_window)
        self._close_process_window.setChecked(
            defaults.close_process_window_when_finished
        )
        self._right_click_buttons[BLAST_PROGRAM_BLASTN].setChecked(True)
        self._update_executable_status()

    def use_system_blast(self) -> bool:
        return self._use_system_blast.isChecked()

    def _update_executable_status(self, *_args) -> None:
        """Refresh PATH status on Options → Executables."""
        if not hasattr(self, "_blast_path_status"):
            return
        set_path_checkbox_status(
            self._blast_path_status,
            resolve_system_executable("blastn"),
            checked=self._use_system_blast.isChecked(),
            has_builtin=True,
        )
        if (
            self.isVisible()
            and self._tabs.currentIndex() == self._options_tab_index
        ):
            QTimer.singleShot(0, self._fit_to_content)

    def current_dialog_tab(self) -> str:
        """Active BLAST window tab key (blastn / blastp / blastx / online / options)."""
        return {
            0: BLAST_PROGRAM_BLASTN,
            1: BLAST_PROGRAM_BLASTP,
            2: BLAST_PROGRAM_BLASTX,
            3: BLAST_PROGRAM_ONLINE,
            4: DIALOG_TAB_OPTIONS,
        }.get(self._tabs.currentIndex(), BLAST_PROGRAM_BLASTN)

    def options(self) -> BlastOptions:
        tab = self.current_dialog_tab()
        program = tab if tab in (
            BLAST_PROGRAM_BLASTN,
            BLAST_PROGRAM_BLASTP,
            BLAST_PROGRAM_BLASTX,
            BLAST_PROGRAM_ONLINE,
        ) else self._initial_program
        # Define-input paths are session-only; always reopen on Use selected.
        # Shared limits: read from the blastn (canonical) copies; all pages stay synced.
        return BlastOptions(
            program=program,
            e_value=normalize_blast_evalue(self._evalue.value()),
            num_threads=self._threads.value(),
            word_size=self._word_size.value(),
            gap_open=self._gap_open.value(),
            gap_extend=self._gap_extend.value(),
            penalty=self._penalty.value(),
            percent_identity=self._percent_identity.value(),
            query_coverage_hsp=self._query_coverage.value(),
            max_msa_lines=self._max_msa.value(),
            top_msa_lines=self._top_msa.value(),
            write_text=self._blastn_write_text.isChecked(),
            write_gff=False,
            minimum_bit_score=DEFAULT_BLAST_MINIMUM_BIT_SCORE,
            blast_database_name=self._blastn_database.currentData(),
            blastp_e_value=normalize_blast_evalue(self._blastp_evalue.value()),
            blastp_num_threads=self._blastp_threads.value(),
            blastp_word_size=self._blastp_word_size.value(),
            blastp_gap_open=self._blastp_gap_open.value(),
            blastp_gap_extend=self._blastp_gap_extend.value(),
            blastp_query_coverage_hsp=self._blastp_query_coverage.value(),
            blastp_database_name=self._blastp_database.currentData(),
            blastp_write_text=self._blastp_write_text.isChecked(),
            blastp_write_gff=False,
            blastx_e_value=normalize_blast_evalue(self._blastx_evalue.value()),
            blastx_num_threads=self._blastx_threads.value(),
            blastx_word_size=self._blastx_word_size.value(),
            blastx_gap_open=self._blastx_gap_open.value(),
            blastx_gap_extend=self._blastx_gap_extend.value(),
            blastx_query_coverage_hsp=self._blastx_query_coverage.value(),
            blastx_database_name=self._blastx_database.currentData(),
            blastx_write_text=self._blastx_write_text.isChecked(),
            blastx_write_gff=False,
            online_program=str(
                self._online_program.currentData() or DEFAULT_BLAST_ONLINE_PROGRAM
            ),
            online_database=str(
                self._online_database.currentData() or DEFAULT_BLAST_ONLINE_DATABASE
            ),
            online_e_value=normalize_blast_evalue(self._online_evalue.value()),
            online_hitlist_size=int(self._online_hitlist.value()),
            online_write_text=self._online_write_text.isChecked(),
            input_source=ALIGN_INPUT_SELECTED,
            input_paths=[],
            last_dialog_tab=persistable_dialog_tab(
                tab, first_tab=BLAST_PROGRAM_BLASTN
            ),
            right_click_program=self._selected_right_click_program(),
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=(
                self._close_process_window.isChecked()
            ),
        )

    def _selected_right_click_program(self) -> str:
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                return value
        return BLAST_PROGRAM_BLASTN


_CONSENSUS_SEQUENCE_TYPE_HELP = (
    "<b>Auto-detect</b>: infer nucleotide vs protein from alignment "
    "residues (gaps ignored).<br><br>"
    "<b>Nucleotide</b>: treat columns as DNA/RNA (A, C, G, T/U, and IUPAC "
    "ambiguity codes).<br><br>"
    "<b>Protein</b>: treat columns as amino acids; use <b>X</b> (or "
    "<b>N</b>) for ambiguous positions."
)


class ConsensusParamsDialog(QDialog):
    """Consensus parameters.

    Left-click the toolbar button opens this dialog. The Parameters tab holds
    threshold and ambiguous base; the Options tab sets what right-clicking the
    toolbar button does.
    """

    def __init__(
        self,
        *,
        initial: ConsensusOptions | None = None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Consensus parameters")
        opts = initial or ConsensusOptions()
        self._for_run = for_run
        self._run_after = False

        self._sequence_type = QComboBox(self)
        self._sequence_type.addItem("Auto-detect", MAFFT_SEQUENCE_AUTO)
        self._sequence_type.addItem("Nucleotide", MAFFT_SEQUENCE_NUCLEOTIDE)
        self._sequence_type.addItem("Protein", MAFFT_SEQUENCE_PROTEIN)
        seq_idx = self._sequence_type.findData(
            normalize_consensus_sequence_type(opts.sequence_type)
        )
        self._sequence_type.setCurrentIndex(seq_idx if seq_idx >= 0 else 0)
        self._sequence_type.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        self._sequence_type.currentIndexChanged.connect(self._on_sequence_type_changed)

        self._threshold = QDoubleSpinBox()
        self._threshold.setRange(0.0, 1.0)
        self._threshold.setDecimals(2)
        self._threshold.setSingleStep(0.05)
        self._threshold.setValue(opts.threshold)
        self._threshold.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        self._threshold.setToolTip(
            "Minimum fraction of sequences that must share a residue for it to "
            "appear in the consensus (0–1)."
        )

        self._ambiguous = QComboBox()
        for base in ("N", "X"):
            self._ambiguous.addItem(base)
        amb_index = self._ambiguous.findText(opts.ambiguous.upper()[:1])
        if amb_index >= 0:
            self._ambiguous.setCurrentIndex(amb_index)
        else:
            self._ambiguous.setCurrentIndex(0)
        self._ambiguous.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        self._ambiguous.setToolTip(
            "Character used when no residue reaches the majority threshold."
        )
        self._ambiguous_label = QLabel("Ambiguous base:", self)

        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("Run Consensus with saved parameters", CONSENSUS_RIGHT_CLICK_RUN),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[CONSENSUS_RIGHT_CLICK_RUN].setChecked(True)

        rc_box = QGroupBox("Right-click action", self)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(
            QLabel("Right-clicking the Consensus button does:", self)
        )
        for radio in self._right_click_buttons.values():
            rc_box_layout.addWidget(radio)

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(12, 10, 12, 10)
        options_layout.setSpacing(8)
        options_layout.addWidget(rc_box)

        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        self._tabs.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        params_tab = QWidget(self)
        params_form = QFormLayout(params_tab)
        params_form.setContentsMargins(12, 10, 12, 10)
        params_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        params_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        params_form.setHorizontalSpacing(12)
        params_form.setVerticalSpacing(8)
        params_form.addRow(
            "Sequence type:",
            field_with_help(
                self._sequence_type,
                "Sequence type",
                _CONSENSUS_SEQUENCE_TYPE_HELP,
            ),
        )
        params_form.addRow("Majority threshold:", self._threshold)
        params_form.addRow(
            self._ambiguous_label,
            field_with_help(
                self._ambiguous,
                "Ambiguous character",
                "Character written at columns where no residue meets the "
                "majority threshold. Use <b>N</b> for nucleotide alignments "
                "and <b>X</b> for protein (Biopython consensus).",
            ),
        )
        self._tabs.addTab(params_tab, "Parameters")
        self._on_sequence_type_changed()
        self._tabs.addTab(options_tab, "Options")

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore],
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 10, 12, 10)
        layout.setSpacing(8)
        layout.addWidget(self._tabs)
        layout.addWidget(self._button_box)
        self._tabs.currentChanged.connect(self._on_tab_changed)
        self._on_tab_changed(self._tabs.currentIndex())

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        # Dialog polish calls adjustSize() after Show; re-shrink once that runs.
        QTimer.singleShot(0, self._fit_to_current_tab)

    def _on_tab_changed(self, tab_index: int) -> None:
        self._fit_to_current_tab()

    def _fit_to_current_tab(self) -> None:
        """Size the dialog to the active tab only (QTabWidget otherwise uses max of all)."""
        page = self._tabs.currentWidget()
        if page is None:
            return
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        page.ensurePolished()
        page.adjustSize()
        root = page.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        page_hint = page.minimumSizeHint().expandedTo(page.sizeHint())
        content_h = page_hint.height()
        if root is not None:
            content_h = max(
                content_h,
                root.totalMinimumSize().height(),
                root.totalSizeHint().height(),
            )
        pin_tool_tabs_height(self._tabs, content_h)
        tab_bar_h = max(24, self._tabs.tabBar().sizeHint().height())
        tabs_h = tab_bar_h + content_h + 8
        buttons_h = self._button_box.sizeHint().height()
        margins = self.layout().contentsMargins()
        spacing = self.layout().spacing()
        height = (
            tabs_h
            + buttons_h
            + margins.top()
            + margins.bottom()
            + spacing
            + 4
        )
        width = max(
            420,
            page_hint.width() + margins.left() + margins.right() + 24,
            self._button_box.sizeHint().width() + margins.left() + margins.right(),
        )
        self.setMinimumWidth(width)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _restore_defaults(self) -> None:
        self._sequence_type.setCurrentIndex(
            max(0, self._sequence_type.findData(MAFFT_SEQUENCE_AUTO))
        )
        self._threshold.setValue(DEFAULT_CONSENSUS_THRESHOLD)
        amb_index = self._ambiguous.findText(DEFAULT_CONSENSUS_AMBIGUOUS)
        self._ambiguous.setCurrentIndex(amb_index if amb_index >= 0 else 0)
        self._right_click_buttons[CONSENSUS_RIGHT_CLICK_RUN].setChecked(True)
        self._on_sequence_type_changed()

    def _on_sequence_type_changed(self, *_args) -> None:
        protein = self._sequence_type.currentData() == MAFFT_SEQUENCE_PROTEIN
        self._ambiguous_label.setText(
            "Ambiguous residue:" if protein else "Ambiguous base:"
        )
        if protein and self._ambiguous.currentText() == "N":
            x_idx = self._ambiguous.findText("X")
            if x_idx >= 0:
                self._ambiguous.setCurrentIndex(x_idx)
        elif not protein and self._sequence_type.currentData() != MAFFT_SEQUENCE_AUTO:
            if self._ambiguous.currentText() == "X":
                n_idx = self._ambiguous.findText("N")
                if n_idx >= 0:
                    self._ambiguous.setCurrentIndex(n_idx)

    def options(self) -> ConsensusOptions:
        right_click_action = CONSENSUS_RIGHT_CLICK_RUN
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click_action = value
                break
        return ConsensusOptions(
            threshold=self._threshold.value(),
            ambiguous=normalize_consensus_ambiguous(self._ambiguous.currentText()),
            sequence_type=str(self._sequence_type.currentData()),
            right_click_action=right_click_action,
        )


_BUILTIN_DOTPLOT_MODE_HELP = (
    "<b>Self comparison</b>: each selected FASTA is compared with itself "
    "(writes <code>*_self_dot.paf</code>, reopenable in Preview).<br><br>"
    "<b>Pairwise comparison</b>: each sequence is compared with every other "
    "(writes one <code>*_multi_dot.paf</code>). Needs two or more sequences. "
    "With three or more sequences, self-comparisons fill the diagonal blocks; "
    "with exactly two, only the cross-comparison is written. Preview mirrors "
    "each pair across the diagonal so both triangles of the matrix are filled."
)
_BUILTIN_DOTPLOT_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): Built-in Dotplot uses the FASTA "
    "files selected in the Files panel when you click Run.<br><br>"
    "<b>Define inputs</b>: choose different FASTA files. Paths are "
    "for this run only and are not remembered next time."
)
_BUILTIN_DOTPLOT_WORD_HELP = (
    "K-mer / word size for the built-in self-dotplot (same engine as Preview → "
    "Dotplot (Self)). Larger values are faster but less sensitive."
)
_BUILTIN_DOTPLOT_REVERSE_HELP = (
    "Include reverse-complement matches (inverted repeats, TIRs). Clear to "
    "plot forward-strand matches only."
)
_FLEXIDOT_MODE_HELP = (
    "<b>Self comparison</b>: each selected FASTA is plotted against itself "
    "(one Flexidot run per file).<br><br>"
    "<b>Pairwise comparison</b>: all selected sequences are compared with each "
    "other in one Flexidot run (pairwise for two sequences, all-vs-all for more)."
)
_FLEXIDOT_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): Flexidot uses the FASTA files that "
    "were selected in the Files panel when you click Run.<br><br>"
    "<b>Define inputs</b>: choose different FASTA files. Paths are "
    "for this run only and are not remembered next time."
)
_FLEXIDOT_OUTPUT_FORMAT_HELP = (
    f"Flexidot output format (<code>-f</code>). HTML gives interactive dotplots "
    f"with Go-to → Annotation inside {APP_DISPLAY_NAME} (forked Flexidot); it "
    "falls back to PNG automatically if that build isn't installed. "
    "PNG/SVG/PDF write a static image."
)
_FLEXIDOT_WORDSIZE_HELP = (
    f"K-mer length (<code>-k</code>), {MIN_FLEXIDOT_WORDSIZE}–"
    f"{MAX_FLEXIDOT_WORDSIZE}. Larger values are faster but less sensitive."
)
_FLEXIDOT_SEQUENCE_TYPE_HELP = (
    "Flexidot <code>-t</code>: sequence alphabet. DNA (nuc) is the default; "
    "choose Protein (aa) for amino-acid FASTA such as .faa files."
)
_FLEXIDOT_GFF_HELP = (
    "Shade the dotplot with a companion annotation GFF beside each FASTA "
    "(e.g. <code>Name.RV.gff3</code> from Annotation).<br><br>"
    "BLAST hit tables (<code>*_b.gff</code>) are skipped. "
    "If none is found, Flexidot runs without shading."
)
_FLEXIDOT_SUBSTITUTIONS_HELP = (
    "Flexidot <code>-S</code>: substitutions allowed inside each k-mer window. "
    "Zero uses the much faster exact-match engine."
)
_FLEXIDOT_WOBBLE_HELP = (
    "Treat ambiguity codes such as R, Y, and N as their possible bases "
    "(<code>-w</code>). This can substantially increase computation time."
)
_FLEXIDOT_MAX_N_HELP = (
    "Flexidot <code>--max_n</code>: maximum N percentage allowed in a k-mer when "
    "ambiguous-base expansion is enabled."
)
_FLEXIDOT_REVERSE_HELP = (
    "Enabled by default. Clear it to pass Flexidot <code>-r</code> / "
    "<code>--norev</code> and omit reverse-complement matches."
)
_FLEXIDOT_LINE_WIDTH_HELP = (
    "Flexidot <code>-A</code> / <code>--line_width</code>: thickness of forward "
    "and reverse-complement match lines in the generated dotplot."
)
_FLEXIDOT_LENGTH_SCALING_HELP = (
    "Scale each dotplot to the sequence lengths instead of a square panel "
    "(<code>-L</code>). Only affects pairwise and all-vs-all comparisons."
)


class FlexidotParamsDialog(QDialog):
    """Dotplot parameters: Built-in (self PAF) + Flexidot, BLAST-style tabs.

    Left-click the toolbar button opens this dialog. Built-in writes
    ``*_self_dot.paf`` / ``*_multi_dot.paf`` for Preview; Flexidot keeps
    HTML/PNG/SVG/PDF outputs. The Options tab sets right-click behaviour.
    """

    _DIALOG_WIDTH = TOOL_PARAMS_DIALOG_WIDTH
    _SETTINGS_LABEL_REF = "Output format:"
    _BUILTIN_SETTINGS_LABEL_REF = "Word size:"

    def __init__(
        self,
        *,
        initial: FlexidotOptions | None = None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        from repeatviewer.backends.dotplot import (
            DEFAULT_WORD_SIZE as _BUILTIN_DEFAULT_WORD,
            MAX_WORD_SIZE as _BUILTIN_MAX_WORD,
            MIN_WORD_SIZE as _BUILTIN_MIN_WORD,
        )

        self.setWindowTitle("Dotplot parameters")
        self.setMinimumWidth(self._DIALOG_WIDTH)
        opts = initial or FlexidotOptions()
        self._for_run = for_run
        self._run_after = False
        # Persisted but not exposed in the compact Advanced panel.
        self._lcs_shading = opts.lcs_shading
        self._only_vs_first_value = opts.only_vs_first
        self._flexidot_compact_size: QSize | None = None
        self._last_engine_tab = normalize_dotplot_dialog_tab(opts.last_dialog_tab)

        # ---- Built-in (*_self_dot.paf / *_multi_dot.paf) ----
        # builtin_combined True = Pairwise; False = Self (legacy key name).
        self._builtin_self = QRadioButton("Self comparison", self)
        self._builtin_pairs = QRadioButton("Pairwise comparison", self)
        self._builtin_mode_group = QButtonGroup(self)
        self._builtin_mode_group.addButton(self._builtin_self)
        self._builtin_mode_group.addButton(self._builtin_pairs)
        if opts.builtin_combined:
            self._builtin_pairs.setChecked(True)
        else:
            self._builtin_self.setChecked(True)
        # Compat aliases used by older tests.
        self._builtin_combine = self._builtin_pairs
        self._builtin_separate = self._builtin_self

        self._builtin_input_group = QButtonGroup(self)
        self._builtin_input_selected = QRadioButton("Use selected sequences", self)
        self._builtin_input_define = QRadioButton("Define inputs", self)
        self._builtin_input_group.addButton(self._builtin_input_selected)
        self._builtin_input_group.addButton(self._builtin_input_define)
        self._builtin_input_selected.setChecked(True)
        self._builtin_define_inputs = TeMultiInputFilesWidget(
            self,
            initial=None,
            start_empty=False,
            placeholder="Click in Files/Folders or Browse…",
            browse_title="Choose a FASTA file",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
            ),
            add_button_label="Add another input",
            add_button_tooltip=(
                "Add another FASTA file path. Click a file in the Files panel "
                "to fill a row, or use Browse."
            ),
            path_mode="file",
        )
        self._builtin_input_selected.toggled.connect(self._sync_builtin_inputs)
        self._builtin_input_define.toggled.connect(self._sync_builtin_inputs)
        self._builtin_define_inputs.paths_changed.connect(self._sync_builtin_inputs)

        self._builtin_word_size = QSpinBox()
        self._builtin_word_size.setRange(_BUILTIN_MIN_WORD, _BUILTIN_MAX_WORD)
        self._builtin_word_size.setValue(
            max(
                _BUILTIN_MIN_WORD,
                min(
                    _BUILTIN_MAX_WORD,
                    int(opts.builtin_word_size or _BUILTIN_DEFAULT_WORD),
                ),
            )
        )
        self._builtin_word_size.setFixedWidth(260)
        self._builtin_include_reverse = QCheckBox(
            "Include reverse-complement matches", self
        )
        self._builtin_include_reverse.setChecked(bool(opts.builtin_include_reverse))

        builtin_tab = QWidget(self)
        builtin_root = QVBoxLayout(builtin_tab)
        builtin_root.setContentsMargins(16, 10, 16, 10)
        builtin_root.setSpacing(4)
        builtin_basic = QWidget(builtin_tab)
        self._builtin_basic = builtin_basic
        builtin_basic.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )
        builtin_form = QFormLayout(builtin_basic)
        builtin_form.setContentsMargins(0, 0, 0, 0)
        builtin_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        builtin_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        builtin_form.setHorizontalSpacing(12)
        builtin_form.setVerticalSpacing(4)

        def _builtin_section(title: str) -> QLabel:
            return QLabel(f"<b>{title}</b>", builtin_tab)

        builtin_form.addRow(_builtin_section("Mode"))
        builtin_mode_box = QWidget(builtin_basic)
        builtin_mode_layout = QVBoxLayout(builtin_mode_box)
        builtin_mode_layout.setContentsMargins(0, 0, 0, 0)
        builtin_mode_layout.setSpacing(2)
        builtin_mode_layout.addWidget(
            field_with_help(
                self._builtin_self,
                "Built-in Dotplot mode",
                _BUILTIN_DOTPLOT_MODE_HELP,
                stretch=True,
            )
        )
        builtin_mode_layout.addWidget(self._builtin_pairs)
        self._builtin_mode_label = QLabel("")
        builtin_form.addRow(self._builtin_mode_label, builtin_mode_box)

        builtin_form.addRow(_form_section_separator())
        builtin_form.addRow(_builtin_section("Inputs"))
        builtin_inputs_box = QWidget(builtin_basic)
        builtin_inputs_layout = QVBoxLayout(builtin_inputs_box)
        builtin_inputs_layout.setContentsMargins(0, 0, 0, 0)
        builtin_inputs_layout.setSpacing(2)
        builtin_inputs_layout.addWidget(
            field_with_help(
                self._builtin_input_selected,
                "Built-in Dotplot inputs",
                _BUILTIN_DOTPLOT_INPUTS_HELP,
                stretch=True,
            )
        )
        builtin_inputs_layout.addWidget(self._builtin_input_define)
        self._builtin_inputs_radios_label = QLabel("")
        self._builtin_inputs_define_label = QLabel("")
        builtin_form.addRow(self._builtin_inputs_radios_label, builtin_inputs_box)
        builtin_form.addRow(
            self._builtin_inputs_define_label, self._builtin_define_inputs
        )

        builtin_form.addRow(_form_section_separator())
        builtin_form.addRow(_builtin_section("Settings"))
        builtin_form.addRow(
            "Word size:",
            field_with_help(
                self._builtin_word_size,
                "Word size",
                _BUILTIN_DOTPLOT_WORD_HELP,
            ),
        )
        self._builtin_reverse_spacer = QLabel("")
        builtin_form.addRow(
            self._builtin_reverse_spacer,
            field_with_help(
                self._builtin_include_reverse,
                "Reverse-complement matches",
                _BUILTIN_DOTPLOT_REVERSE_HELP,
            ),
        )
        self._builtin_basic_form = builtin_form
        builtin_root.addWidget(builtin_basic)
        builtin_root.addStretch(1)
        self._builtin_tab = builtin_tab

        # ---- Flexidot engine ----
        # combined True = Pairwise; False = Self (legacy key name).
        self._flexidot_self = QRadioButton("Self comparison", self)
        self._flexidot_pairs = QRadioButton("Pairwise comparison", self)
        self._mode_group = QButtonGroup(self)
        self._mode_group.addButton(self._flexidot_self)
        self._mode_group.addButton(self._flexidot_pairs)
        if opts.combined:
            self._flexidot_pairs.setChecked(True)
        else:
            self._flexidot_self.setChecked(True)
        # Compat aliases used by older tests.
        self._combine = self._flexidot_pairs
        self._separate = self._flexidot_self

        self._input_group = QButtonGroup(self)
        self._input_selected = QRadioButton("Use selected sequences", self)
        self._input_define = QRadioButton("Define inputs", self)
        self._input_group.addButton(self._input_selected)
        self._input_group.addButton(self._input_define)
        self._input_selected.setChecked(True)
        self._define_inputs = TeMultiInputFilesWidget(
            self,
            initial=None,
            start_empty=False,
            placeholder="Click in Files/Folders or Browse…",
            browse_title="Choose a FASTA file",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
            ),
            add_button_label="Add another input",
            add_button_tooltip=(
                "Add another FASTA file path. Click a file in the Files panel "
                "to fill a row, or use Browse."
            ),
            path_mode="file",
        )
        self._input_selected.toggled.connect(self._sync_flexidot_inputs)
        self._input_define.toggled.connect(self._sync_flexidot_inputs)
        self._define_inputs.paths_changed.connect(self._sync_flexidot_inputs)

        self._wordsize = QSpinBox()
        self._wordsize.setRange(MIN_FLEXIDOT_WORDSIZE, MAX_FLEXIDOT_WORDSIZE)
        self._wordsize.setValue(
            max(MIN_FLEXIDOT_WORDSIZE, min(MAX_FLEXIDOT_WORDSIZE, opts.wordsize))
        )
        self._wordsize.setToolTip(
            f"K-mer length (-k), {MIN_FLEXIDOT_WORDSIZE}–{MAX_FLEXIDOT_WORDSIZE}."
        )

        self._sequence_type = QComboBox()
        self._sequence_type.addItem("DNA / nucleotide (-t nuc)", FLEXIDOT_SEQUENCE_NUC)
        self._sequence_type.addItem("Protein (-t aa)", FLEXIDOT_SEQUENCE_AA)
        seq_type = normalize_flexidot_sequence_type(opts.sequence_type)
        seq_idx = self._sequence_type.findData(seq_type)
        self._sequence_type.setCurrentIndex(seq_idx if seq_idx >= 0 else 0)
        self._sequence_type.currentIndexChanged.connect(self._on_sequence_type_changed)

        self._use_gff = QCheckBox("Show GFF features on dotplot")
        self._use_gff.setChecked(opts.use_gff)

        self._output_format = QComboBox()
        self._output_format.addItem(
            f"HTML (interactive; Go-to in {APP_DISPLAY_NAME})",
            FLEXIDOT_OUTPUT_HTML,
        )
        self._output_format.addItem("PNG (static image)", FLEXIDOT_OUTPUT_PNG)
        self._output_format.addItem("SVG (vector)", FLEXIDOT_OUTPUT_SVG)
        self._output_format.addItem("PDF (vector)", FLEXIDOT_OUTPUT_PDF)
        fmt_index = self._output_format.findData(
            normalize_flexidot_output_format(opts.output_format)
        )
        self._output_format.setCurrentIndex(fmt_index if fmt_index >= 0 else 0)

        self._substitutions = QSpinBox()
        self._substitutions.setRange(0, 20)
        self._substitutions.setValue(opts.substitution_count)

        self._wobble = QCheckBox("Expand ambiguous IUPAC bases (-w)")
        self._wobble.setChecked(opts.wobble_conversion)

        self._max_n = QDoubleSpinBox()
        self._max_n.setRange(0.0, 100.0)
        self._max_n.setDecimals(1)
        self._max_n.setSuffix(" %")
        self._max_n.setValue(opts.max_n_percentage)
        self._wobble.toggled.connect(self._max_n.setEnabled)
        self._max_n.setEnabled(self._wobble.isChecked())

        self._include_reverse = QCheckBox(
            "Calculate reverse-complement matches (disable: -r / --norev)"
        )
        self._include_reverse.setChecked(opts.include_reverse)

        self._line_width = QDoubleSpinBox()
        self._line_width.setRange(0.1, 20.0)
        self._line_width.setDecimals(1)
        self._line_width.setSingleStep(0.2)
        self._line_width.setValue(opts.line_width)

        self._length_scaling = QCheckBox(
            "Scale plots by sequence length (-L, pairwise / all-vs-all)"
        )
        self._length_scaling.setChecked(opts.length_scaling)

        for field in (
            self._wordsize,
            self._sequence_type,
            self._output_format,
            self._substitutions,
            self._max_n,
            self._line_width,
        ):
            field.setFixedWidth(260)

        params_tab = QWidget(self)
        params_root = QVBoxLayout(params_tab)
        params_root.setContentsMargins(16, 10, 16, 10)
        params_root.setSpacing(4)

        basic = QWidget(params_tab)
        self._flexidot_basic = basic
        basic.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )
        form = QFormLayout(basic)
        form.setContentsMargins(0, 0, 0, 0)
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        # Expanding so Define-inputs path + Browse fill the dialog width.
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)

        def _section(title: str, parent: QWidget | None = None) -> QLabel:
            return QLabel(f"<b>{title}</b>", parent or params_tab)

        form.addRow(_section("Mode", basic))
        mode_box = QWidget(basic)
        mode_layout = QVBoxLayout(mode_box)
        mode_layout.setContentsMargins(0, 0, 0, 0)
        mode_layout.setSpacing(2)
        mode_layout.addWidget(
            field_with_help(
                self._flexidot_self,
                "Flexidot mode",
                _FLEXIDOT_MODE_HELP,
                stretch=True,
            )
        )
        mode_layout.addWidget(self._flexidot_pairs)
        self._mode_label = QLabel("")
        form.addRow(self._mode_label, mode_box)

        form.addRow(_form_section_separator())
        form.addRow(_section("Inputs", basic))
        inputs_box = QWidget(basic)
        inputs_layout = QVBoxLayout(inputs_box)
        inputs_layout.setContentsMargins(0, 0, 0, 0)
        inputs_layout.setSpacing(2)
        inputs_layout.addWidget(
            field_with_help(
                self._input_selected,
                "Flexidot inputs",
                _FLEXIDOT_INPUTS_HELP,
                stretch=True,
            )
        )
        inputs_layout.addWidget(self._input_define)
        self._inputs_radios_label = QLabel("")
        self._inputs_define_label = QLabel("")
        form.addRow(self._inputs_radios_label, inputs_box)
        form.addRow(self._inputs_define_label, self._define_inputs)

        form.addRow(_form_section_separator())
        form.addRow(_section("Outputs", basic))
        form.addRow(
            "Output format:",
            field_with_help(
                self._output_format, "Output format", _FLEXIDOT_OUTPUT_FORMAT_HELP
            ),
        )

        form.addRow(_form_section_separator())
        form.addRow(_section("Settings", basic))
        form.addRow(
            "k-mer:",
            field_with_help(self._wordsize, "k-mer", _FLEXIDOT_WORDSIZE_HELP),
        )
        form.addRow(
            "Sequence type:",
            field_with_help(
                self._sequence_type, "Sequence type", _FLEXIDOT_SEQUENCE_TYPE_HELP
            ),
        )
        self._use_gff_spacer = QLabel("")
        form.addRow(
            self._use_gff_spacer,
            field_with_help(self._use_gff, "Show GFF features", _FLEXIDOT_GFF_HELP),
        )
        self._flexidot_basic_form = form
        params_root.addWidget(basic)

        advanced = QWidget(params_tab)
        self._flexidot_advanced_panel = advanced
        advanced.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        advanced_form = QFormLayout(advanced)
        advanced_form.setContentsMargins(0, 4, 0, 0)
        advanced_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        advanced_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        advanced_form.setHorizontalSpacing(12)
        advanced_form.setVerticalSpacing(4)
        advanced_form.addRow(_form_section_separator())
        advanced_form.addRow(_section("Advanced Settings", advanced))
        advanced_form.addRow(
            "Allowed substitutions:",
            field_with_help(
                self._substitutions,
                "Allowed substitutions",
                _FLEXIDOT_SUBSTITUTIONS_HELP,
            ),
        )
        advanced_form.addRow(
            "Maximum N content:",
            field_with_help(self._max_n, "Maximum N content", _FLEXIDOT_MAX_N_HELP),
        )
        advanced_form.addRow(
            "Match line thickness:",
            field_with_help(
                self._line_width, "Match line thickness", _FLEXIDOT_LINE_WIDTH_HELP
            ),
        )
        self._wobble_spacer = QLabel("")
        self._include_reverse_spacer = QLabel("")
        self._length_scaling_spacer = QLabel("")
        advanced_form.addRow(
            self._wobble_spacer,
            field_with_help(
                self._wobble, "Expand ambiguous IUPAC bases", _FLEXIDOT_WOBBLE_HELP
            ),
        )
        advanced_form.addRow(
            self._include_reverse_spacer,
            field_with_help(
                self._include_reverse,
                "Reverse-complement matches",
                _FLEXIDOT_REVERSE_HELP,
            ),
        )
        advanced_form.addRow(
            self._length_scaling_spacer,
            field_with_help(
                self._length_scaling,
                "Length scaling",
                _FLEXIDOT_LENGTH_SCALING_HELP,
            ),
        )
        self._flexidot_advanced_form = advanced_form
        params_root.addWidget(advanced)
        advanced.setVisible(False)
        advanced.setMaximumHeight(0)
        self._params_tab = params_tab
        self._sync_field_column()

        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        saved_rc = normalize_flexidot_right_click_action(opts.right_click_action)
        for label, value in (
            ("Run Built-in Dotplot", FLEXIDOT_RIGHT_CLICK_BUILTIN),
            ("Run Flexidot", FLEXIDOT_RIGHT_CLICK_FLEXIDOT),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(saved_rc == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[FLEXIDOT_RIGHT_CLICK_BUILTIN].setChecked(True)

        rc_box = QGroupBox("Right-click action", self)
        rc_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(
            QLabel("Right-clicking the Dotplot button does:", self)
        )
        for radio in self._right_click_buttons.values():
            rc_box_layout.addWidget(radio)

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="Dotplot",
                checked=bool(getattr(opts, "show_process_window", True)),
                close_when_finished=bool(
                    getattr(opts, "close_process_window_when_finished", True)
                ),
            )
        )

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)
        options_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        options_layout.addWidget(rc_box)
        options_layout.addWidget(process_box)
        from repeatviewer.tool_citations import CITATION_DOTPLOT

        options_layout.addWidget(make_citation_box(self, CITATION_DOTPLOT))
        options_layout.addStretch(1)
        options_tab.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._options_tab = options_tab

        self._tabs = QTabWidget(self)
        self._tabs.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._tabs.addTab(builtin_tab, "Built-in")
        self._tabs.addTab(params_tab, "Flexidot")
        self._tabs.addTab(options_tab, "Options")
        open_tab = normalize_dotplot_dialog_tab(opts.last_dialog_tab)
        if open_tab == DOTPLOT_DIALOG_TAB_FLEXIDOT:
            self._tabs.setCurrentWidget(params_tab)
        else:
            self._tabs.setCurrentWidget(builtin_tab)
        self._tabs.currentChanged.connect(self._on_flexidot_tab_changed)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._advanced_settings_btn = QPushButton("Advanced Settings ▸", self)
        self._advanced_settings_btn.setCheckable(True)
        self._advanced_settings_btn.setChecked(False)
        self._advanced_settings_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._advanced_settings_btn.toggled.connect(self._on_advanced_settings_toggled)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore, self._advanced_settings_btn],
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(14)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._tabs, 0)
        layout.addWidget(self._button_box, 0)

        self._sync_builtin_inputs(resize=False)
        self._sync_flexidot_inputs(resize=False)
        self._on_sequence_type_changed()
        self._on_flexidot_tab_changed(self._tabs.currentIndex())
        self.resize(self._DIALOG_WIDTH, 420)
        QTimer.singleShot(0, self._fit_active_dotplot_tab)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._flexidot_compact_size = None
        self._sync_field_column()
        self._sync_builtin_field_column()
        QTimer.singleShot(0, self._fit_active_dotplot_tab)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def active_engine_tab(self) -> str:
        """Return ``builtin`` or ``flexidot`` for the active Parameters engine."""
        current = self._tabs.currentWidget()
        if current is self._params_tab:
            return DOTPLOT_DIALOG_TAB_FLEXIDOT
        if current is self._builtin_tab:
            return DOTPLOT_DIALOG_TAB_BUILTIN
        # Options tab: keep last engine selection from options().
        return normalize_dotplot_dialog_tab(
            getattr(self, "_last_engine_tab", DOTPLOT_DIALOG_TAB_BUILTIN)
        )

    def input_source(self) -> str:
        if self.active_engine_tab() == DOTPLOT_DIALOG_TAB_BUILTIN:
            if self._builtin_input_define.isChecked():
                return ALIGN_INPUT_DEFINE
            return ALIGN_INPUT_SELECTED
        if self._input_define.isChecked():
            return ALIGN_INPUT_DEFINE
        return ALIGN_INPUT_SELECTED

    def run_files(self, *, panel_files: list[Path] | None = None) -> list[Path]:
        from repeatviewer.backends.fastga_dotplot import expand_fastga_input_paths

        if self.input_source() == ALIGN_INPUT_DEFINE:
            widget = (
                self._builtin_define_inputs
                if self.active_engine_tab() == DOTPLOT_DIALOG_TAB_BUILTIN
                else self._define_inputs
            )
            return expand_fastga_input_paths(widget.paths())
        return [Path(p) for p in (panel_files or [])]

    def apply_path_from_panel(self, path: Path) -> bool:
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        if self._tabs.currentWidget() is self._options_tab:
            return False
        try:
            candidate = Path(path).expanduser().resolve()
        except OSError:
            return False
        widget = (
            self._builtin_define_inputs
            if self.active_engine_tab() == DOTPLOT_DIALOG_TAB_BUILTIN
            else self._define_inputs
        )
        return bool(widget.accept_path(candidate))

    def apply_folder_from_panel(self, folder: Path) -> bool:
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        if self._tabs.currentWidget() is self._options_tab:
            return False
        try:
            candidate = Path(folder).expanduser().resolve()
        except OSError:
            return False
        if not candidate.is_dir():
            return False
        widget = (
            self._builtin_define_inputs
            if self.active_engine_tab() == DOTPLOT_DIALOG_TAB_BUILTIN
            else self._define_inputs
        )
        return bool(widget.accept_path(candidate))

    def _sync_builtin_inputs(self, *, resize: bool = True) -> None:
        define = self._builtin_input_define.isChecked()
        form = self._builtin_basic.layout()
        if isinstance(form, QFormLayout):
            form.setRowVisible(self._builtin_define_inputs, define)
        self._builtin_define_inputs.setVisible(define)
        if define:
            self._builtin_define_inputs.setMinimumHeight(0)
            self._builtin_define_inputs.setMaximumHeight(16777215)
        else:
            self._builtin_define_inputs.setFixedHeight(0)
        self._builtin_define_inputs.updateGeometry()
        self._builtin_basic.updateGeometry()
        self._sync_builtin_field_column()
        if (
            resize
            and self.isVisible()
            and self._tabs.currentWidget() is self._builtin_tab
        ):
            self._fit_builtin_content()

    def _sync_flexidot_inputs(self, *, resize: bool = True) -> None:
        define = self._input_define.isChecked()
        form = self._flexidot_basic.layout()
        if isinstance(form, QFormLayout):
            form.setRowVisible(self._define_inputs, define)
        self._define_inputs.setVisible(define)
        if define:
            self._define_inputs.setMinimumHeight(0)
            self._define_inputs.setMaximumHeight(16777215)
        else:
            self._define_inputs.setFixedHeight(0)
        self._define_inputs.updateGeometry()
        self._flexidot_basic.updateGeometry()
        self._sync_field_column()
        if (
            resize
            and self.isVisible()
            and not self._advanced_settings_btn.isChecked()
            and self._tabs.currentWidget() is self._params_tab
        ):
            self._flexidot_compact_size = None
            self._resize_for_flexidot_advanced(False)

    def _sync_builtin_field_column(self) -> None:
        if not hasattr(self, "_builtin_inputs_radios_label"):
            return
        sync_form_field_column_spacers(
            self,
            [
                self._builtin_mode_label,
                self._builtin_inputs_radios_label,
                self._builtin_inputs_define_label,
                self._builtin_reverse_spacer,
            ],
            [getattr(self, "_builtin_basic_form", None)],
            ref_label=self._BUILTIN_SETTINGS_LABEL_REF,
        )

    def _sync_field_column(self) -> None:
        if not hasattr(self, "_inputs_radios_label"):
            return
        sync_form_field_column_spacers(
            self,
            [
                self._mode_label,
                self._inputs_radios_label,
                self._inputs_define_label,
                self._use_gff_spacer,
                self._wobble_spacer,
                self._include_reverse_spacer,
                self._length_scaling_spacer,
            ],
            [
                getattr(self, "_flexidot_basic_form", None),
                getattr(self, "_flexidot_advanced_form", None),
            ],
            ref_label=self._SETTINGS_LABEL_REF,
        )

    def _fit_active_dotplot_tab(self) -> None:
        current = self._tabs.currentWidget()
        if current is self._builtin_tab:
            self._fit_builtin_content()
        elif current is self._params_tab:
            self._resize_for_flexidot_advanced(
                self._advanced_settings_btn.isChecked()
            )
        else:
            self._fit_options_content()

    def _fit_builtin_content(self) -> None:
        page = self._builtin_tab
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        root = page.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        page.updateGeometry()
        basic_h = max(
            self._builtin_basic.minimumSizeHint().height(),
            self._builtin_basic.sizeHint().height(),
        )
        margins = root.contentsMargins() if root is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 20
        spacing = root.spacing() if root is not None else 4
        content_h = margin_h + basic_h + spacing
        chrome = self._dialog_chrome_height()
        width = tool_params_dialog_width(
            define_inputs=self._builtin_input_define.isChecked()
        )
        self.setMinimumWidth(width)
        height = chrome + content_h
        height = min(max(height, 220), 10_000)
        pin_tool_tabs_height(self._tabs, content_h)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def _on_flexidot_tab_changed(self, _index: int) -> None:
        current = self._tabs.currentWidget()
        on_flexidot = current is self._params_tab
        on_builtin = current is self._builtin_tab
        on_options = current is self._options_tab
        self._advanced_settings_btn.setVisible(on_flexidot)
        if on_flexidot or on_builtin:
            self._last_engine_tab = (
                DOTPLOT_DIALOG_TAB_FLEXIDOT if on_flexidot else DOTPLOT_DIALOG_TAB_BUILTIN
            )
        if not on_flexidot:
            if self._advanced_settings_btn.isChecked():
                self._advanced_settings_btn.blockSignals(True)
                self._advanced_settings_btn.setChecked(False)
                self._advanced_settings_btn.setText("Advanced Settings ▸")
                self._advanced_settings_btn.blockSignals(False)
                self._set_flexidot_advanced_visible(False)
        if not self.isVisible():
            return
        if on_options:
            self._fit_options_content()
        elif on_builtin:
            self._fit_builtin_content()
        else:
            self._resize_for_flexidot_advanced(False)

    def _set_flexidot_advanced_visible(self, visible: bool) -> None:
        panel = self._flexidot_advanced_panel
        panel.setVisible(visible)
        if visible:
            panel.setMaximumHeight(16777215)
        else:
            panel.setMaximumHeight(0)
        panel.updateGeometry()
        root = self._params_tab.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        self._params_tab.updateGeometry()
        self._tabs.updateGeometry()

    def _on_advanced_settings_toggled(self, checked: bool) -> None:
        self._advanced_settings_btn.setText(
            "Advanced Settings ▾" if checked else "Advanced Settings ▸"
        )
        if checked:
            self._flexidot_compact_size = self.size()
        self._set_flexidot_advanced_visible(checked)
        self._resize_for_flexidot_advanced(checked)

    def _dialog_chrome_height(self) -> int:
        lay = self.layout()
        margins = lay.contentsMargins() if lay is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 34
        spacing = lay.spacing() if lay is not None else 14
        tab_bar_h = self._tabs.tabBar().sizeHint().height() + 10
        buttons_h = self._button_box.sizeHint().height()
        return margin_h + spacing + tab_bar_h + buttons_h

    def _flexidot_compact_content_height(self) -> int:
        root = self._params_tab.layout()
        margins = root.contentsMargins() if root is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 20
        spacing = root.spacing() if root is not None else 4
        basic_h = max(
            self._flexidot_basic.minimumSizeHint().height(),
            self._flexidot_basic.sizeHint().height(),
        )
        return margin_h + basic_h + spacing

    def _resize_for_flexidot_advanced(self, expanded: bool) -> None:
        if self._tabs.currentWidget() is not self._params_tab:
            self._fit_options_content()
            return
        panel = self._params_tab
        adv = self._flexidot_advanced_panel
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        root = panel.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        panel.updateGeometry()
        self._tabs.updateGeometry()

        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        chrome = self._dialog_chrome_height()
        width = tool_params_dialog_width(
            define_inputs=self.input_source() == ALIGN_INPUT_DEFINE
        )
        self.setMinimumWidth(width)
        compact_h = self._flexidot_compact_content_height()

        if not expanded:
            height = chrome + compact_h
            height = min(max(height, 220), max_h)
            pin_tool_tabs_height(self._tabs, compact_h)
            self.setMaximumHeight(height)
            self.resize(width, height)
            self.setMaximumHeight(16777215)
            self._flexidot_compact_size = QSize(width, height)
            return

        extra = max(0, adv.sizeHint().height())
        page_h = compact_h + extra
        base = self._flexidot_compact_size or self.size()
        height = min(base.height() + extra + 12, max_h)
        height = max(height, min(chrome + page_h, max_h))
        pin_tool_tabs_height(self._tabs, page_h)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def _fit_options_content(self) -> None:
        page = self._options_tab
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        page.ensurePolished()
        page.adjustSize()
        root = page.layout()
        if root is not None:
            root.activate()
        page_hint = page.minimumSizeHint().expandedTo(page.sizeHint())
        content_h = page_hint.height()
        if root is not None:
            content_h = max(
                content_h,
                root.totalMinimumSize().height(),
                root.totalSizeHint().height(),
            )
        chrome = self._dialog_chrome_height()
        height = chrome + content_h
        pin_tool_tabs_height(self._tabs, content_h)
        self.setMaximumHeight(height)
        self.resize(self._DIALOG_WIDTH, height)
        self.setMaximumHeight(16777215)

    def _on_sequence_type_changed(self, *_args) -> None:
        is_protein = (
            str(self._sequence_type.currentData()) == FLEXIDOT_SEQUENCE_AA
        )
        self._include_reverse.setEnabled(not is_protein)
        if is_protein:
            self._include_reverse.setChecked(False)

    def _restore_defaults(self) -> None:
        defaults = FlexidotOptions()
        if defaults.builtin_combined:
            self._builtin_pairs.setChecked(True)
        else:
            self._builtin_self.setChecked(True)
        self._builtin_input_selected.setChecked(True)
        self._builtin_define_inputs.clear_paths()
        self._sync_builtin_inputs()
        self._builtin_word_size.setValue(defaults.builtin_word_size)
        self._builtin_include_reverse.setChecked(defaults.builtin_include_reverse)
        if defaults.combined:
            self._flexidot_pairs.setChecked(True)
        else:
            self._flexidot_self.setChecked(True)
        self._input_selected.setChecked(True)
        self._define_inputs.clear_paths()
        self._sync_flexidot_inputs()
        self._wordsize.setValue(defaults.wordsize)
        self._sequence_type.setCurrentIndex(
            max(0, self._sequence_type.findData(defaults.sequence_type))
        )
        self._output_format.setCurrentIndex(
            max(0, self._output_format.findData(defaults.output_format))
        )
        self._lcs_shading = defaults.lcs_shading
        self._use_gff.setChecked(defaults.use_gff)
        self._substitutions.setValue(defaults.substitution_count)
        self._wobble.setChecked(defaults.wobble_conversion)
        self._max_n.setValue(defaults.max_n_percentage)
        self._include_reverse.setChecked(defaults.include_reverse)
        self._only_vs_first_value = defaults.only_vs_first
        self._line_width.setValue(defaults.line_width)
        self._length_scaling.setChecked(defaults.length_scaling)
        self._show_process_window.setChecked(defaults.show_process_window)
        self._close_process_window.setChecked(
            defaults.close_process_window_when_finished
        )
        self._right_click_buttons[FLEXIDOT_RIGHT_CLICK_BUILTIN].setChecked(True)
        self._on_sequence_type_changed()

    def options(self) -> FlexidotOptions:
        fmt = self._output_format.currentData()
        right_click_action = FLEXIDOT_RIGHT_CLICK_BUILTIN
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click_action = value
                break
        right_click_action = normalize_flexidot_right_click_action(right_click_action)
        seq_type = normalize_flexidot_sequence_type(self._sequence_type.currentData())
        include_reverse = self._include_reverse.isChecked()
        if seq_type == FLEXIDOT_SEQUENCE_AA:
            include_reverse = False
        engine = self.active_engine_tab()
        if self._tabs.currentWidget() is self._options_tab:
            engine = normalize_dotplot_dialog_tab(
                getattr(self, "_last_engine_tab", DOTPLOT_DIALOG_TAB_BUILTIN)
            )
        return FlexidotOptions(
            wordsize=self._wordsize.value(),
            sequence_type=seq_type,
            lcs_shading=self._lcs_shading,
            use_gff=self._use_gff.isChecked(),
            output_format=normalize_flexidot_output_format(
                str(fmt or FLEXIDOT_OUTPUT_HTML)
            ),
            substitution_count=self._substitutions.value(),
            wobble_conversion=self._wobble.isChecked(),
            max_n_percentage=self._max_n.value(),
            include_reverse=include_reverse,
            only_vs_first=self._only_vs_first_value,
            line_width=self._line_width.value(),
            length_scaling=self._length_scaling.isChecked(),
            combined=self._flexidot_pairs.isChecked(),
            # Flexidot -m stays automatic from sequence count when Pairwise is on.
            modes=DEFAULT_FLEXIDOT_MODES,
            right_click_action=right_click_action,
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=(
                self._close_process_window.isChecked()
            ),
            builtin_word_size=int(self._builtin_word_size.value()),
            builtin_include_reverse=bool(self._builtin_include_reverse.isChecked()),
            builtin_combined=bool(self._builtin_pairs.isChecked()),
            last_dialog_tab=engine,
        )


# Public alias matching the Tools → Dotplot name.
DotplotParamsDialog = FlexidotParamsDialog


_TE_ANNOTATION_PARAMS_WIDTH = TOOL_PARAMS_DIALOG_WIDTH
# Cap Advanced Settings scroll height on short screens.
_TE_ANNOTATION_ADVANCED_SCROLL_MAX_HEIGHT = 420

_TE_ANN_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): Sequence Annotation uses the FASTA "
    "files that were selected in the Files panel when you click Run.<br><br>"
    "<b>Define inputs</b>: choose different FASTA files. Paths are "
    "for this run only and are not remembered next time."
)
_TE_ANN_FEATURE_DETECTION_HELP = (
    "Detect sequence features without an external database:<br><br>"
    "<b>LTR</b>: long terminal repeats (direct repeats typical of LTR "
    "retrotransposons).<br>"
    "<b>TIR</b>: terminal inverted repeats (typical of DNA transposons).<br>"
    "<b>Tandem</b>: tandem / satellite-style repeats.<br>"
    "<b>ORF</b>: open reading frames.<br>"
    "<b>TSD</b>: target site duplications near putative TE boundaries "
    "(MSA only; writes features to the annotation JSON for Overview).<br>"
    "<b>PBS</b>: primer binding site inside the 5′ LTR (Neumann 2019 motifs + "
    "blastn; requires two LTR features).<br>"
    "<b>All</b>: select or clear every detector (except TSD and PBS).<br><br>"
    "Thresholds (lengths, similarity, TSD / PBS settings, …) are under "
    "<b>Advanced Settings</b> → Feature detection settings. Threads are under "
    "<b>Settings</b> on this page (shared with Similarity annotation). The "
    "process window is under <b>Options</b>. Tools → TSD search also writes a "
    "clustered text report."
)
_TE_ANN_SIMILARITY_HELP = (
    "Search selected protein / DNA / HMM / CDD databases for similarity hits "
    "on each sequence. Coloured names match the feature-track colours.<br><br>"
    "Add or edit databases in <b>File → Manage Database</b>. "
    "<b>All</b> selects or clears every listed database.<br><br>"
    "E-value and hit limits are under <b>Advanced Settings</b> → "
    "Similarity annotation settings. Threads are under <b>Settings</b> on this "
    "page (shared with Feature detection). The process window is under "
    "<b>Options</b>."
)


def _section_title_label(text: str, parent: QWidget | None = None) -> QLabel:
    label = QLabel(text, parent)
    font = label.font()
    font.setBold(True)
    label.setFont(font)
    return label


def _section_title_with_help(
    text: str,
    help_title: str,
    help_text: str,
    parent: QWidget | None = None,
) -> QWidget:
    """Bold section title with a trailing ``?`` (Inputs / Feature detection style)."""
    row = QWidget(parent)
    layout = QHBoxLayout(row)
    layout.setContentsMargins(0, 0, 0, 0)
    layout.setSpacing(6)
    layout.addWidget(_section_title_label(text, row))
    layout.addWidget(make_help_button(row, help_title, help_text))
    layout.addStretch(1)
    return row


class _CheckboxFlowWidget(QWidget):
    """Lay out checkboxes in wrapping rows within a fixed content width."""

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        spacing: int = 10,
        row_spacing: int = 4,
        content_width: int = _TE_ANNOTATION_PARAMS_WIDTH - 40,
    ) -> None:
        super().__init__(parent)
        self._spacing = spacing
        self._row_spacing = row_spacing
        self._content_width = max(200, content_width)
        self._items: list[QWidget] = []
        self._last_flow_width = -1
        self._rows_layout = QVBoxLayout(self)
        self._rows_layout.setContentsMargins(0, 0, 0, 0)
        self._rows_layout.setSpacing(self._row_spacing)
        self.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)

    def set_items(self, items: list[QWidget]) -> None:
        self._items = list(items)
        for item in self._items:
            item.setParent(self)
        self._last_flow_width = -1
        self._reflow()

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._reflow()

    def sizeHint(self):  # noqa: N802
        from PySide6.QtCore import QSize

        return self._size_for_width(self._content_width)

    def minimumSizeHint(self):  # noqa: N802
        return self.sizeHint()

    def _available_width(self) -> int:
        width = self.width()
        if width > 40:
            return width
        return self._content_width

    def _pack_rows(self, available: int) -> list[list[QWidget]]:
        rows: list[list[QWidget]] = []
        current: list[QWidget] = []
        current_w = 0
        for item in self._items:
            item_w = max(1, item.sizeHint().width())
            need = item_w if not current else item_w + self._spacing
            if current and current_w + need > available:
                rows.append(current)
                current = [item]
                current_w = item_w
            else:
                current.append(item)
                current_w += need
        if current:
            rows.append(current)
        return rows

    def _size_for_width(self, available: int):
        from PySide6.QtCore import QSize

        if not self._items:
            return QSize(available, 0)
        rows = self._pack_rows(available)
        height = 0
        for index, row in enumerate(rows):
            height += max(chip.sizeHint().height() for chip in row)
            if index:
                height += self._row_spacing
        return QSize(available, height)

    def _reflow(self) -> None:
        available = self._available_width()
        if available == self._last_flow_width and self._rows_layout.count() > 0:
            return
        self._last_flow_width = available
        while self._rows_layout.count():
            item = self._rows_layout.takeAt(0)
            row_widget = item.widget()
            if row_widget is not None:
                row_layout = row_widget.layout()
                if row_layout is not None:
                    while row_layout.count():
                        child = row_layout.takeAt(0)
                        child_widget = child.widget()
                        if child_widget is not None:
                            child_widget.setParent(self)
                row_widget.deleteLater()

        if not self._items:
            self.updateGeometry()
            return

        available = self._available_width()
        for row_items in self._pack_rows(available):
            row_widget = QWidget(self)
            row_layout = QHBoxLayout(row_widget)
            row_layout.setContentsMargins(0, 0, 0, 0)
            row_layout.setSpacing(self._spacing)
            for chip in row_items:
                chip.setParent(row_widget)
                row_layout.addWidget(chip, 0, Qt.AlignmentFlag.AlignVCenter)
            row_layout.addStretch(1)
            self._rows_layout.addWidget(row_widget)
        self.updateGeometry()


class TeAnnotationDialog(QDialog):
    """Batch TE annotation: Parameters select modules; Advanced holds settings/colors."""

    def __init__(
        self,
        *,
        initial: TeAnnotationJobOptions | None = None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        from repeatviewer.backends.auto_annotation import AutoAnnotationMode
        from repeatviewer.backends.annotation_colors import similarity_database_color
        from repeatviewer.backends.auto_annotation_settings import (
            default_annotation_threads,
            load_annotation_run_options,
            load_auto_annotation_settings,
        )
        from repeatviewer.backends.similarity_annotation import SimilarityDatabaseRun
        from repeatviewer.backends.similarity_annotation_settings import (
            load_similarity_annotation_settings,
        )
        from repeatviewer.backends.te_annotation_batch import TeAnnotationOptions
        from repeatviewer.database_settings import load_similarity_annotation_databases
        from repeatviewer.widgets.annotation_colors_dialog import (
            AnnotationElementColorsPanel,
        )
        from repeatviewer.widgets.sequence_map_viewer import (
            AutoAnnotationSettingsPanel,
            SimilarityAnnotationSettingsPanel,
            TsdSearchSettingsPanel,
        )

        super().__init__(parent)
        self.setWindowTitle("Sequence Annotation")
        self._TsdSearchSettingsPanel = TsdSearchSettingsPanel
        self.setMinimumWidth(_TE_ANNOTATION_PARAMS_WIDTH)
        opts = initial or TeAnnotationJobOptions()
        self._for_run = for_run
        self._run_after = False
        self._te_ann_compact_size: QSize | None = None

        self._TeAnnotationOptions = TeAnnotationOptions
        self._SimilarityDatabaseRun = SimilarityDatabaseRun
        self._similarity_database_color = similarity_database_color
        self._AutoAnnotationMode = AutoAnnotationMode
        self._auto_all_updating = False
        self._similarity_all_updating = False

        selected_modes = set(opts.auto_modes)

        # --- Inputs (BLAST / Align figure-2 style) -------------------------
        self._input_group = QButtonGroup(self)
        self._input_selected = QRadioButton("Use selected sequences", self)
        self._input_define = QRadioButton("Define inputs", self)
        self._input_group.addButton(self._input_selected)
        self._input_group.addButton(self._input_define)
        self._input_selected.setChecked(True)
        self._define_inputs = TeMultiInputFilesWidget(
            self,
            initial=None,
            start_empty=False,
            placeholder="Click in Files/Folders or Browse…",
            browse_title="Choose a FASTA file",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
            ),
            add_button_label="Add another input",
            add_button_tooltip=(
                "Add another FASTA file path. Click a file in the Files panel "
                "to fill a row, or use Browse."
            ),
            path_mode="file",
        )

        # --- Feature detection (formerly Auto annotation) -----------------
        feature_section = QWidget(self)
        feature_outer = QVBoxLayout(feature_section)
        feature_outer.setContentsMargins(0, 0, 0, 0)
        feature_outer.setSpacing(6)
        feature_outer.addWidget(
            _section_title_with_help(
                "Feature detection",
                "Feature detection",
                _TE_ANN_FEATURE_DETECTION_HELP,
                feature_section,
            )
        )

        auto_checks_row = QHBoxLayout()
        auto_checks_row.setSpacing(10)
        self._auto_ann_all_check = QCheckBox("All", feature_section)
        self._auto_ann_all_check.setToolTip(
            "Select or clear all Feature detection modules"
        )
        self._auto_ann_all_check.toggled.connect(self._on_auto_all_toggled)
        auto_checks_row.addWidget(self._auto_ann_all_check)

        self._auto_ann_checks: dict = {}
        for label, mode in (
            ("LTR", AutoAnnotationMode.LTR),
            ("TIR", AutoAnnotationMode.TIR),
            ("Tandem", AutoAnnotationMode.TANDEM),
            ("ORF", AutoAnnotationMode.ORF),
            ("TSD", AutoAnnotationMode.TSD),
            ("PBS", AutoAnnotationMode.PBS),
        ):
            check = QCheckBox(label, feature_section)
            check.setChecked(mode.value in selected_modes)
            if mode == AutoAnnotationMode.TSD:
                check.setToolTip(
                    "Target site duplications near putative TE boundaries "
                    "(MSA only). Configure boundaries under Advanced → "
                    "Feature detection settings / Tools → TSD search."
                )
            elif mode == AutoAnnotationMode.PBS:
                check.setToolTip(
                    "Primer binding site inside the 5′ LTR (Neumann 2019 "
                    "motifs + blastn). Requires two LTR features. Configure "
                    "under Advanced → Feature detection settings."
                )
            check.toggled.connect(self._on_auto_check_toggled)
            self._auto_ann_checks[mode] = check
            auto_checks_row.addWidget(check)
        self._refresh_feature_detect_check_colors()

        auto_checks_row.addStretch(1)
        feature_outer.addLayout(auto_checks_row)
        self._sync_auto_all_check()

        # --- Similarity annotation ----------------------------------------
        sim_section = QWidget(self)
        sim_outer = QVBoxLayout(sim_section)
        sim_outer.setContentsMargins(0, 0, 0, 0)
        sim_outer.setSpacing(6)
        sim_outer.addWidget(
            _section_title_with_help(
                "Similarity annotation",
                "Similarity annotation",
                _TE_ANN_SIMILARITY_HELP,
                sim_section,
            )
        )

        self._similarity_ann_all_check = QCheckBox("All", sim_section)
        self._similarity_ann_all_check.setToolTip(
            "Select or clear all similarity databases"
        )
        self._similarity_ann_all_check.toggled.connect(self._on_similarity_all_toggled)

        self._similarity_db_rows: list[tuple] = []
        self._db_flow = _CheckboxFlowWidget(sim_section)

        databases = load_similarity_annotation_databases()
        select_all_dbs = opts.similarity_database_names is None
        selected_dbs = set(opts.similarity_database_names or ())
        if not databases:
            empty = QLabel(
                "No databases defined. Add them in File → Manage Database.",
                sim_section,
            )
            empty.setWordWrap(True)
            sim_outer.addWidget(empty)
            self._similarity_ann_all_check.hide()
            self._db_flow.hide()
        else:
            flow_items: list[QWidget] = [self._similarity_ann_all_check]
            for index, entry in enumerate(databases):
                check = self._make_similarity_db_check(
                    entry,
                    index,
                    self._db_flow,
                    checked=select_all_dbs or entry.name in selected_dbs,
                )
                self._similarity_db_rows.append((entry, check))
                flow_items.append(check)
            self._db_flow.set_items(flow_items)
            sim_outer.addWidget(self._db_flow)
        self._sync_similarity_all_check()

        # Thresholds/colors stay in Advanced; Threads + process window are on Settings.
        self._auto_settings_panel = AutoAnnotationSettingsPanel(
            load_auto_annotation_settings(), self, include_run_controls=False
        )
        self._similarity_settings_panel = SimilarityAnnotationSettingsPanel(
            load_similarity_annotation_settings(),
            self,
            include_run_controls=False,
        )
        # Feature types + databases only (no nucleotide / amino-acid palettes).
        self._colors_panel = AnnotationElementColorsPanel(
            self, show_residue_colors=False
        )

        auto_box = QGroupBox("Feature detection settings", self)
        auto_box_layout = QVBoxLayout(auto_box)
        auto_box_layout.addWidget(self._auto_settings_panel)
        # Sequence Annotation: RepeatViewer boundaries only (no User-define popup).
        self._tsd_settings_panel = self._TsdSearchSettingsPanel(
            auto_box, allow_user_boundary=False
        )
        auto_box_layout.addWidget(self._tsd_settings_panel)
        auto_box_layout.addWidget(self._auto_settings_panel.pbs_section_widget())

        sim_box = QGroupBox("Similarity annotation settings", self)
        sim_box_layout = QVBoxLayout(sim_box)
        sim_box_layout.addWidget(self._similarity_settings_panel)

        color_box = QGroupBox("Element colors", self)
        color_box_layout = QVBoxLayout(color_box)
        color_box_layout.addWidget(self._colors_panel)

        advanced_body = QWidget(self)
        advanced_body_layout = QVBoxLayout(advanced_body)
        advanced_body_layout.setContentsMargins(0, 0, 0, 0)
        advanced_body_layout.setSpacing(10)
        advanced_body_layout.addWidget(auto_box)
        advanced_body_layout.addWidget(sim_box)
        advanced_body_layout.addWidget(color_box)

        advanced_scroll = QScrollArea(self)
        advanced_scroll.setWidgetResizable(True)
        advanced_scroll.setFrameShape(QFrame.Shape.NoFrame)
        advanced_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        advanced_scroll.setWidget(advanced_body)
        advanced_scroll.setMaximumHeight(_TE_ANNOTATION_ADVANCED_SCROLL_MAX_HEIGHT)
        self._advanced_scroll = advanced_scroll

        advanced = QWidget(self)
        advanced_layout = QVBoxLayout(advanced)
        advanced_layout.setContentsMargins(0, 8, 0, 0)
        advanced_layout.setSpacing(8)
        advanced_layout.addWidget(_section_title_label("Advanced Settings", advanced))
        advanced_layout.addWidget(advanced_scroll)
        advanced.setVisible(False)
        advanced.setMaximumHeight(0)
        self._advanced_panel = advanced

        params_tab = QWidget(self)
        params_root = QVBoxLayout(params_tab)
        params_root.setContentsMargins(16, 10, 16, 10)
        params_root.setSpacing(4)
        params_root.setAlignment(Qt.AlignmentFlag.AlignTop)

        basic = QWidget(params_tab)
        self._te_ann_basic = basic
        basic.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )
        form = QFormLayout(basic)
        form.setContentsMargins(0, 0, 0, 0)
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)
        self._te_ann_basic_form = form

        # Inputs: title + radios flush-left (same rhythm as Feature detection /
        # BLAST title→control distance). Avoid the wide Settings label-column
        # indent — this page has no Settings colon fields to align to.
        inputs_section = QWidget(basic)
        inputs_outer = QVBoxLayout(inputs_section)
        inputs_outer.setContentsMargins(0, 0, 0, 0)
        inputs_outer.setSpacing(6)
        inputs_outer.addWidget(_section_title_label("Inputs", inputs_section))
        inputs_box = QWidget(inputs_section)
        inputs_layout = QVBoxLayout(inputs_box)
        inputs_layout.setContentsMargins(0, 0, 0, 0)
        inputs_layout.setSpacing(2)
        inputs_layout.addWidget(
            field_with_help(
                self._input_selected,
                "Sequence Annotation inputs",
                _TE_ANN_INPUTS_HELP,
                stretch=True,
            )
        )
        inputs_layout.addWidget(self._input_define)
        inputs_outer.addWidget(inputs_box)
        inputs_outer.addWidget(self._define_inputs)
        form.addRow(inputs_section)

        form.addRow(_form_section_separator())
        form.addRow(feature_section)
        form.addRow(_form_section_separator())
        form.addRow(sim_section)
        form.addRow(_form_section_separator())

        # Settings: Threads only (process window lives on Options, like BLAST).
        run_threads, _run_show, _run_close = load_annotation_run_options()
        settings_section = QWidget(basic)
        settings_outer = QVBoxLayout(settings_section)
        settings_outer.setContentsMargins(0, 0, 0, 0)
        settings_outer.setSpacing(6)
        settings_outer.addWidget(_section_title_label("Settings", settings_section))
        self._threads = QSpinBox(settings_section)
        self._threads.setRange(1, 256)
        self._threads.setValue(max(1, int(run_threads)))
        self._threads.setFixedWidth(260)
        threads_row = QWidget(settings_section)
        threads_row_layout = QHBoxLayout(threads_row)
        threads_row_layout.setContentsMargins(0, 0, 0, 0)
        threads_row_layout.setSpacing(8)
        threads_label = QLabel("Threads:", threads_row)
        threads_label.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter
        )
        threads_row_layout.addWidget(threads_label, 0, Qt.AlignmentFlag.AlignLeft)
        threads_row_layout.addWidget(
            field_with_help(
                self._threads,
                "Threads",
                "CPU threads inside each sequence job for Feature detection "
                "and Similarity annotation (like BLAST Threads). How many "
                "sequences run at once is Preferences → Max. parallel jobs. "
                f"Default is {default_annotation_threads()}.",
            ),
            0,
            Qt.AlignmentFlag.AlignLeft,
        )
        threads_row_layout.addStretch(1)
        settings_outer.addWidget(threads_row)
        form.addRow(settings_section)

        params_root.addWidget(basic)
        params_root.addWidget(advanced)
        params_tab.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._params_tab = params_tab

        self._input_selected.toggled.connect(self._sync_te_ann_inputs)
        self._input_define.toggled.connect(self._sync_te_ann_inputs)
        self._define_inputs.paths_changed.connect(self._on_te_ann_define_paths_changed)
        self._sync_te_ann_inputs(resize=False)

        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("Run Sequence annotation with saved selection", TE_ANNOTATION_RIGHT_CLICK_RUN),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[TE_ANNOTATION_RIGHT_CLICK_RUN].setChecked(True)

        rc_box = QGroupBox("Right-click action", self)
        rc_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(
            QLabel("Right-clicking the Sequence Annotation button does:", self)
        )
        for radio in self._right_click_buttons.values():
            rc_box_layout.addWidget(radio)

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="Sequence Annotation",
                checked=bool(_run_show),
                close_when_finished=bool(_run_close),
            )
        )

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)
        options_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        options_layout.addWidget(rc_box)
        options_layout.addWidget(process_box)
        from repeatviewer.tool_citations import CITATION_SEQUENCE_ANNOTATION

        options_layout.addWidget(
            make_citation_box(self, CITATION_SEQUENCE_ANNOTATION)
        )
        options_layout.addStretch(1)
        options_tab.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._options_tab = options_tab

        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        self._tabs.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._tabs.addTab(params_tab, "Parameters")
        self._tabs.addTab(options_tab, "Options")
        self._options_tab_index = 1

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._advanced_settings_btn = QPushButton("Advanced Settings ▸", self)
        self._advanced_settings_btn.setCheckable(True)
        self._advanced_settings_btn.setChecked(False)
        self._advanced_settings_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._advanced_settings_btn.toggled.connect(self._on_advanced_settings_toggled)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore, self._advanced_settings_btn],
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(8)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._tabs, 0)
        layout.addWidget(self._button_box, 0)
        self._tabs.currentChanged.connect(self._on_tab_changed)
        self._on_tab_changed(self._tabs.currentIndex())

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._te_ann_compact_size = None
        QTimer.singleShot(0, lambda: self._resize_for_advanced(False))

    def _on_tab_changed(self, tab_index: int) -> None:
        on_params = tab_index == 0
        self._advanced_settings_btn.setVisible(on_params)
        if not on_params:
            self._advanced_settings_btn.blockSignals(True)
            self._advanced_settings_btn.setChecked(False)
            self._advanced_settings_btn.blockSignals(False)
            self._advanced_settings_btn.setText("Advanced Settings ▸")
            self._advanced_panel.setVisible(False)
            self._advanced_panel.setMaximumHeight(0)
        QTimer.singleShot(0, lambda: self._resize_for_advanced(False))

    def _on_advanced_settings_toggled(self, checked: bool) -> None:
        self._advanced_settings_btn.setText(
            "Advanced Settings ▾" if checked else "Advanced Settings ▸"
        )
        self._resize_for_advanced(checked)

    def _te_ann_dialog_chrome_height(self) -> int:
        margins = self.layout().contentsMargins() if self.layout() else None
        chrome = (
            (margins.top() + margins.bottom()) if margins is not None else 34
        ) + 8
        tab_bar_h = max(24, self._tabs.tabBar().sizeHint().height()) + 8
        buttons_h = self._button_box.sizeHint().height()
        return chrome + tab_bar_h + buttons_h

    def _te_ann_dialog_width(self) -> int:
        margins = self.layout().contentsMargins() if self.layout() else None
        width = tool_params_dialog_width(
            define_inputs=self.input_source() == ALIGN_INPUT_DEFINE
        )
        return max(
            width,
            self._button_box.sizeHint().width()
            + (margins.left() + margins.right() if margins else 44),
        )

    def input_source(self) -> str:
        if self._input_define.isChecked():
            return ALIGN_INPUT_DEFINE
        return ALIGN_INPUT_SELECTED

    def run_files(self, *, panel_files: list[Path] | None = None) -> list[Path]:
        from repeatviewer.backends.fastga_dotplot import expand_fastga_input_paths
        from repeatviewer.backends.te_annotation_batch import is_annotatable_fasta

        if self.input_source() == ALIGN_INPUT_DEFINE:
            expanded = expand_fastga_input_paths(self._define_inputs.paths())
            return [p for p in expanded if is_annotatable_fasta(p)]
        return [Path(p) for p in (panel_files or [])]

    def apply_path_from_panel(self, path: Path) -> bool:
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        if self._tabs.currentWidget() is self._options_tab:
            return False
        try:
            candidate = Path(path).expanduser().resolve()
        except OSError:
            return False
        return bool(self._define_inputs.accept_path(candidate))

    def apply_folder_from_panel(self, folder: Path) -> bool:
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        if self._tabs.currentWidget() is self._options_tab:
            return False
        try:
            candidate = Path(folder).expanduser().resolve()
        except OSError:
            return False
        if not candidate.is_dir():
            return False
        return bool(self._define_inputs.accept_path(candidate))

    def _on_te_ann_define_paths_changed(self) -> None:
        self._sync_te_ann_inputs(resize=True)

    def _sync_te_ann_inputs(self, *, resize: bool = True) -> None:
        define = self.input_source() == ALIGN_INPUT_DEFINE
        self._define_inputs.setVisible(define)
        if define:
            self._define_inputs.setMinimumHeight(0)
            self._define_inputs.setMaximumHeight(16777215)
        else:
            self._define_inputs.setFixedHeight(0)
        self._define_inputs.updateGeometry()
        if hasattr(self, "_te_ann_basic"):
            self._te_ann_basic.updateGeometry()
        if (
            resize
            and self.isVisible()
            and not self._advanced_settings_btn.isChecked()
            and self._tabs.currentWidget() is self._params_tab
        ):
            self._te_ann_compact_size = None
            QTimer.singleShot(0, lambda: self._resize_for_advanced(False))

    def _te_ann_compact_page_height(self) -> int:
        """Parameters content height with Advanced collapsed."""
        was_visible = self._advanced_panel.isVisible()
        was_max = self._advanced_panel.maximumHeight()
        self._advanced_panel.setVisible(False)
        self._advanced_panel.setMaximumHeight(0)
        root = self._params_tab.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        height = max(
            self._params_tab.minimumSizeHint().height(),
            self._params_tab.sizeHint().height(),
        )
        self._advanced_panel.setVisible(was_visible)
        self._advanced_panel.setMaximumHeight(was_max)
        if root is not None:
            root.invalidate()
            root.activate()
        return height

    def _te_ann_advanced_extra_height(self) -> int:
        """Header + scroll viewport — not the full unwrapped settings/colors."""
        scroll_h = min(
            _TE_ANNOTATION_ADVANCED_SCROLL_MAX_HEIGHT,
            max(220, self._advanced_scroll.maximumHeight()),
        )
        return 36 + scroll_h + 12

    def _resize_for_advanced(self, expanded: bool) -> None:
        """Grow for Advanced Settings; shrink when collapsed (BLAST-style)."""
        if self._tabs.currentIndex() != 0:
            self._fit_options_tab()
            return
        width = self._te_ann_dialog_width()
        if hasattr(self, "_db_flow") and self._db_flow.isVisible():
            self._db_flow._last_flow_width = -1
            self._db_flow.setFixedWidth(max(200, width - 80))
            self._db_flow._reflow()

        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        chrome = self._te_ann_dialog_chrome_height()
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        self._define_inputs.setMinimumHeight(0)
        self._define_inputs.setMaximumHeight(16777215)
        if hasattr(self, "_te_ann_basic"):
            self._te_ann_basic.setMinimumHeight(0)
            self._te_ann_basic.setMaximumHeight(16777215)
        compact_page_h = self._te_ann_compact_page_height()

        self._advanced_panel.setVisible(expanded)
        if expanded:
            self._advanced_scroll.setMaximumHeight(
                _TE_ANNOTATION_ADVANCED_SCROLL_MAX_HEIGHT
            )
            self._advanced_panel.setMaximumHeight(16777215)
        else:
            self._advanced_panel.setMaximumHeight(0)

        root = self._params_tab.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        self._params_tab.updateGeometry()
        self._tabs.updateGeometry()
        self.setMinimumWidth(width)

        if not expanded:
            height = chrome + compact_page_h
            height = min(max(height, 220), max_h)
            pin_tool_tabs_height(self._tabs, compact_page_h)
            self.setMaximumHeight(height)
            self.resize(width, height)
            self.setMaximumHeight(16777215)
            self._te_ann_compact_size = QSize(width, height)
            return

        extra = self._te_ann_advanced_extra_height()
        page_h = compact_page_h + extra
        base = self._te_ann_compact_size or self.size()
        height = min(base.height() + extra + 12, max_h)
        height = max(
            height,
            min(chrome + page_h, max_h),
        )
        pin_tool_tabs_height(self._tabs, page_h)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def _fit_options_tab(self) -> None:
        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        chrome = self._te_ann_dialog_chrome_height()
        width = self._te_ann_dialog_width()
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        self._options_tab.ensurePolished()
        self._options_tab.adjustSize()
        root = self._options_tab.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        page_h = max(
            self._options_tab.minimumSizeHint().height(),
            self._options_tab.sizeHint().height(),
        )
        if root is not None:
            page_h = max(
                page_h,
                root.totalMinimumSize().height(),
                root.totalSizeHint().height(),
            )
        height = min(max(chrome + page_h, 220), max_h)
        pin_tool_tabs_height(self._tabs, page_h)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _make_similarity_db_check(
        self,
        entry,
        index: int,
        parent: QWidget,
        *,
        checked: bool = True,
    ) -> QCheckBox:
        color = self._similarity_database_color(entry.name, index)
        check = QCheckBox(entry.name, parent)
        check.setToolTip(f"{entry.display_type_label}\n{entry.path}")
        check.setChecked(checked)
        check.setStyleSheet(f"QCheckBox {{ color: {color}; }}")
        check.setSizePolicy(QSizePolicy.Policy.Maximum, QSizePolicy.Policy.Fixed)
        check.toggled.connect(self._on_similarity_check_toggled)
        return check

    def _refresh_similarity_db_checkbox_colors(self) -> None:
        for index, (entry, check) in enumerate(self._similarity_db_rows):
            color = self._similarity_database_color(entry.name, index)
            check.setStyleSheet(f"QCheckBox {{ color: {color}; }}")

    def _refresh_feature_detect_check_colors(self) -> None:
        from repeatviewer.backends.annotation_colors import feature_type_color

        for mode, check in self._auto_ann_checks.items():
            color = feature_type_color(mode.value)
            check.setStyleSheet(f"QCheckBox {{ color: {color}; }}")

    def _on_auto_check_toggled(self, _checked: bool) -> None:
        if self._auto_all_updating:
            return
        self._sync_auto_all_check()

    def _sync_auto_all_check(self) -> None:
        # "All" covers every shown detector, TSD and PBS included (they still
        # start unticked by default).
        from repeatviewer.backends.auto_annotation import AutoAnnotationMode

        checks = [
            check
            for _mode, check in self._auto_ann_checks.items()
            if not check.isHidden()
        ]
        if not checks:
            return
        all_checked = all(check.isChecked() for check in checks)
        self._auto_all_updating = True
        self._auto_ann_all_check.setChecked(all_checked)
        self._auto_all_updating = False

    def _on_auto_all_toggled(self, checked: bool) -> None:
        if self._auto_all_updating:
            return
        from repeatviewer.backends.auto_annotation import AutoAnnotationMode

        self._auto_all_updating = True
        for mode, check in self._auto_ann_checks.items():
            if check.isHidden():
                continue  # e.g. TSD on a single sequence
            check.setChecked(checked)
        self._auto_all_updating = False

    def _on_similarity_check_toggled(self, _checked: bool) -> None:
        if self._similarity_all_updating:
            return
        self._sync_similarity_all_check()

    def _sync_similarity_all_check(self) -> None:
        checks = [check for _entry, check in self._similarity_db_rows]
        if not checks:
            self._similarity_all_updating = True
            self._similarity_ann_all_check.setChecked(False)
            self._similarity_all_updating = False
            return
        all_checked = all(check.isChecked() for check in checks)
        self._similarity_all_updating = True
        self._similarity_ann_all_check.setChecked(all_checked)
        self._similarity_all_updating = False

    def _on_similarity_all_toggled(self, checked: bool) -> None:
        if self._similarity_all_updating:
            return
        self._similarity_all_updating = True
        for _entry, check in self._similarity_db_rows:
            check.setChecked(checked)
        self._similarity_all_updating = False

    def _restore_defaults(self) -> None:
        from repeatviewer.backends.auto_annotation_settings import AutoAnnotationSettings
        from repeatviewer.backends.similarity_annotation_settings import (
            SimilarityAnnotationSettings,
        )

        defaults = TeAnnotationJobOptions()
        self._input_selected.setChecked(True)
        self._define_inputs.clear_paths()
        self._sync_te_ann_inputs(resize=False)
        selected = set(defaults.auto_modes)
        self._auto_all_updating = True
        for mode, check in self._auto_ann_checks.items():
            check.setChecked(mode.value in selected)
        self._auto_all_updating = False
        self._sync_auto_all_check()

        self._similarity_all_updating = True
        for _entry, check in self._similarity_db_rows:
            check.setChecked(True)
        self._similarity_all_updating = False
        self._sync_similarity_all_check()

        # Rebuild settings panels from defaults by replacing child widgets.
        parent_auto = self._auto_settings_panel.parentWidget()
        parent_sim = self._similarity_settings_panel.parentWidget()
        from repeatviewer.backends.auto_annotation_settings import (
            default_annotation_threads,
        )
        from repeatviewer.job_options import TsdSearchOptions
        from repeatviewer.widgets.sequence_map_viewer import (
            AutoAnnotationSettingsPanel,
            SimilarityAnnotationSettingsPanel,
            TsdSearchSettingsPanel,
        )

        new_auto = AutoAnnotationSettingsPanel(
            AutoAnnotationSettings(), parent_auto, include_run_controls=False
        )
        new_sim = SimilarityAnnotationSettingsPanel(
            SimilarityAnnotationSettings(), parent_sim, include_run_controls=False
        )
        new_tsd = TsdSearchSettingsPanel(
            parent_auto, initial=TsdSearchOptions(), allow_user_boundary=False
        )
        new_pbs = new_auto.pbs_section_widget()
        if parent_auto is not None:
            layout = parent_auto.layout()
            if layout is not None:
                layout.replaceWidget(self._auto_settings_panel, new_auto)
                layout.replaceWidget(self._tsd_settings_panel, new_tsd)
                # PBS lives after TSD; drop the old section then append the new one.
                old_pbs = self._auto_settings_panel.pbs_section_widget()
                layout.removeWidget(old_pbs)
                old_pbs.deleteLater()
                layout.addWidget(new_pbs)
        self._auto_settings_panel.deleteLater()
        self._auto_settings_panel = new_auto
        self._tsd_settings_panel.deleteLater()
        self._tsd_settings_panel = new_tsd
        if parent_sim is not None:
            layout = parent_sim.layout()
            if layout is not None:
                layout.replaceWidget(self._similarity_settings_panel, new_sim)
        self._similarity_settings_panel.deleteLater()
        self._similarity_settings_panel = new_sim
        self._colors_panel.restore_defaults()
        self._right_click_buttons[TE_ANNOTATION_RIGHT_CLICK_RUN].setChecked(True)
        self._threads.setValue(default_annotation_threads())
        self._show_process_window.setChecked(defaults.show_process_window)
        self._close_process_window.setChecked(
            defaults.close_process_window_when_finished
        )
        self._close_process_window.setEnabled(self._show_process_window.isChecked())
        if self._advanced_settings_btn.isChecked():
            self._advanced_settings_btn.setChecked(False)
        else:
            self._resize_for_advanced(False)

    def apply_settings_panels(self) -> None:
        """Persist Settings + Advanced feature-detection, similarity, and colors."""
        from dataclasses import replace

        from repeatviewer.backends.auto_annotation_settings import (
            save_annotation_run_options,
            save_auto_annotation_settings,
        )
        from repeatviewer.backends.similarity_annotation_settings import (
            save_similarity_annotation_settings,
        )

        threads = max(1, int(self._threads.value()))
        show = bool(self._show_process_window.isChecked())
        close = bool(self._close_process_window.isChecked())
        save_annotation_run_options(
            num_threads=threads,
            show_process_window=show,
            close_process_window_when_finished=close,
        )
        auto = replace(
            self._auto_settings_panel.settings(),
            num_threads=threads,
            show_process_window=show,
            close_process_window_when_finished=close,
        )
        sim = replace(
            self._similarity_settings_panel.settings(),
            num_threads=threads,
            show_process_window=show,
            close_process_window_when_finished=close,
        )
        save_auto_annotation_settings(auto)
        save_similarity_annotation_settings(sim)
        from repeatviewer.job_options import save_tsd_search_options

        save_tsd_search_options(self._tsd_settings_panel.options())
        self._colors_panel.save()
        self._refresh_similarity_db_checkbox_colors()
        self._refresh_feature_detect_check_colors()

    def job_options(self) -> TeAnnotationJobOptions:
        auto_modes = tuple(
            mode.value
            for mode, check in self._auto_ann_checks.items()
            if check.isChecked()
        )
        selected_dbs = tuple(
            entry.name
            for entry, check in self._similarity_db_rows
            if check.isChecked()
        )
        all_dbs = tuple(entry.name for entry, _check in self._similarity_db_rows)
        similarity_database_names: tuple[str, ...] | None
        if not all_dbs or selected_dbs == all_dbs:
            similarity_database_names = None
        else:
            similarity_database_names = selected_dbs

        right_click_action = TE_ANNOTATION_RIGHT_CLICK_RUN
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click_action = value
                break
        return TeAnnotationJobOptions(
            auto_modes=auto_modes,
            similarity_database_names=similarity_database_names,
            right_click_action=right_click_action,
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=(
                self._close_process_window.isChecked()
            ),
        )

    def options(self):
        """Runtime batch options for ``run_te_annotation_on_fasta``."""
        auto_modes = tuple(
            mode for mode, check in self._auto_ann_checks.items() if check.isChecked()
        )
        similarity_databases: list = []
        for index, (entry, check) in enumerate(self._similarity_db_rows):
            if check.isChecked():
                similarity_databases.append(
                    self._SimilarityDatabaseRun(
                        database=entry,
                        color=self._similarity_database_color(entry.name, index),
                    )
                )
        return self._TeAnnotationOptions(
            auto_modes=auto_modes,
            similarity_databases=tuple(similarity_databases),
        )


def batch_options_from_te_annotation_job(
    job: TeAnnotationJobOptions,
):
    """Build runtime ``TeAnnotationOptions`` from saved job options."""
    from repeatviewer.backends.auto_annotation import AutoAnnotationMode
    from repeatviewer.backends.annotation_colors import similarity_database_color
    from repeatviewer.backends.similarity_annotation import SimilarityDatabaseRun
    from repeatviewer.backends.te_annotation_batch import TeAnnotationOptions
    from repeatviewer.database_settings import load_similarity_annotation_databases

    modes: list[AutoAnnotationMode] = []
    for raw in job.auto_modes:
        try:
            modes.append(AutoAnnotationMode(raw))
        except ValueError:
            continue
    databases = load_similarity_annotation_databases()
    select_all = job.similarity_database_names is None
    selected = set(job.similarity_database_names or ())
    runs: list[SimilarityDatabaseRun] = []
    for index, entry in enumerate(databases):
        if select_all or entry.name in selected:
            runs.append(
                SimilarityDatabaseRun(
                    database=entry,
                    color=similarity_database_color(entry.name, index),
                )
            )
    return TeAnnotationOptions(
        auto_modes=tuple(modes),
        similarity_databases=tuple(runs),
    )


class OrfPredictionParamsDialog(QDialog):
    """ORF Prediction settings for Sequence menu / toolbar (not TE annotation)."""

    def __init__(
        self,
        *,
        initial=None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        from repeatviewer.backends.auto_annotation_settings import (
            OrfSettings as _OrfSettings,
            normalize_orf_frames,
        )

        super().__init__(parent)
        self.setWindowTitle("ORF Prediction settings (orfipy)")
        self.setMinimumWidth(420)
        opts = initial or _OrfSettings()
        self._OrfSettings = _OrfSettings
        self._for_run = for_run
        self._run_after = False

        self._min_len = QSpinBox(self)
        self._min_len.setRange(3, 1_000_000)
        self._min_len.setValue(opts.min_length)
        self._max_len = QSpinBox(self)
        self._max_len.setRange(3, 1_000_000)
        self._max_len.setValue(opts.max_length)
        self._alt_starts = QLineEdit(self)
        self._alt_starts.setPlaceholderText("e.g. TTG, CTG")
        self._alt_starts.setText(opts.alternative_start_codons)
        self._alt_starts.setToolTip("Optional extra start codons (comma- or space-separated).")

        selected = set(normalize_orf_frames(opts.frames))
        frames_inner = QWidget(self)
        frames_inner.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        frames_grid = QGridLayout(frames_inner)
        frames_grid.setContentsMargins(0, 0, 0, 0)
        frames_grid.setHorizontalSpacing(4)
        frames_grid.setVerticalSpacing(2)
        self._frame_checks: dict[int, QCheckBox] = {}
        for col, frame in enumerate((1, 2, 3)):
            frames_grid.setColumnStretch(col, 0)
            fwd = QCheckBox(f"Frame {frame}", frames_inner)
            fwd.setToolTip(f"Forward reading frame {frame}")
            fwd.setChecked(frame in selected)
            self._frame_checks[frame] = fwd
            frames_grid.addWidget(fwd, 0, col, Qt.AlignmentFlag.AlignLeft)
            rev = QCheckBox(f"Reverse {frame}", frames_inner)
            rev.setToolTip(f"Reverse reading frame {frame}")
            rev.setChecked((-frame) in selected)
            self._frame_checks[-frame] = rev
            frames_grid.addWidget(rev, 1, col, Qt.AlignmentFlag.AlignLeft)
        frames_host = QWidget(self)
        frames_host.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        frames_host.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        frames_host_layout = QHBoxLayout(frames_host)
        frames_host_layout.setContentsMargins(0, 0, 0, 0)
        frames_host_layout.setSpacing(0)
        frames_host_layout.addWidget(
            frames_inner,
            0,
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop,
        )
        frames_host_layout.addStretch(1)

        def _bp_row(spin: QSpinBox) -> QWidget:
            row = QWidget(self)
            layout = QHBoxLayout(row)
            layout.setContentsMargins(0, 0, 0, 0)
            layout.setSpacing(6)
            layout.addWidget(spin)
            layout.addWidget(QLabel("bp", row))
            layout.addStretch()
            return row

        self._alt_starts.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        form.addRow("Min ORF length:", _bp_row(self._min_len))
        form.addRow("Max ORF length:", _bp_row(self._max_len))
        form.addRow(
            "Add other start codons:",
            field_with_help(
                self._alt_starts,
                "Other start codons",
                "Optional extra start codons (comma- or space-separated), "
                "e.g. <b>TTG, CTG</b>, in addition to ATG.",
                parent=self,
            ),
        )
        form.addRow(
            "Reading frames:",
            field_with_help(
                frames_host,
                "Reading frames",
                "Keep ORFs from these orfipy reading frames "
                "(Frame 1–3 forward, Reverse 1–3). "
                "Default is all six frames.",
                parent=self,
            ),
        )

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore],
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 10, 12, 10)
        layout.setSpacing(10)
        layout.addLayout(form)
        layout.addWidget(self._button_box)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _restore_defaults(self) -> None:
        from repeatviewer.backends.auto_annotation_settings import DEFAULT_ORF_FRAMES

        defaults = self._OrfSettings()
        self._min_len.setValue(defaults.min_length)
        self._max_len.setValue(defaults.max_length)
        self._alt_starts.setText(defaults.alternative_start_codons)
        selected = set(DEFAULT_ORF_FRAMES)
        for frame, check in self._frame_checks.items():
            check.setChecked(frame in selected)

    def options(self):
        from repeatviewer.backends.auto_annotation_settings import normalize_orf_frames

        frames = normalize_orf_frames(
            tuple(frame for frame, check in self._frame_checks.items() if check.isChecked())
        )
        return self._OrfSettings(
            min_length=self._min_len.value(),
            max_length=max(self._max_len.value(), self._min_len.value()),
            alternative_start_codons=self._alt_starts.text().strip(),
            frames=frames,
        )


_TEAID_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): TE-Aid uses the FASTA files that "
    "were selected in the Files panel when you click Run.<br><br>"
    "<b>Define inputs</b>: choose different FASTA files. Paths are "
    "for this run only and are not remembered next time."
)
_TEAID_MODE_HELP = (
    "<b>Single genome</b>: BLAST the consensus against one reference genome "
    "(project reference genome or a saved nucleotide BLAST database).<br><br>"
    "<b>Pangenome</b>: run TE-Aid on each provided genome and compare the "
    "results in one report (at least two genomes)."
)
_TEAID_PANGENOME_HELP = (
    "Run TE-Aid on each provided genome and compare them in one report "
    "(at least two genomes from the project genome or Manage Database)."
)
_TEAID_EVALUE_HELP = (
    "Reporting E-value for the genome BLAST search. TE-Aid BLAST-searches the "
    "genome first to find copies of the provided sequence, then draws those "
    "hits on the coverage track."
)
_TEAID_WORD_SIZE_HELP = (
    "BLASTN <code>word_size</code> for the self-BLAST step (range 4–25). TE-Aid "
    "self-BLASTs the input sequence to find repetitive regions and draw the "
    "self-alignment / terminal-repeat panels."
)
_TEAID_X_AXIS_LABEL_HELP = (
    "Title of the input-sequence axis in the BLAST hits, coverage and self "
    f"dot plots (default: \"{DEFAULT_TEAID_X_AXIS_LABEL}\"). Change it to "
    "describe the input, e.g. \"Consensus (bp)\" or \"satDNA monomer (bp)\". "
    "Leave empty to use the default."
)
_TEAID_TRANSPARENCY_HELP = (
    "Transparency of genome BLAST hit lines on the coverage track "
    "(0% = fully opaque, 100% = invisible)."
)
_TEAID_THREADS_HELP = (
    "Number of threads for TE-Aid genome BLAST and self-BLAST "
    "(<code>-num_threads</code>)."
)
_TEAID_HIT_HIGH_HELP = (
    "Hits covering at least this percent of the input / consensus length use "
    "the chosen color on the genome coverage track. Click the color swatch to "
    "change it."
)
_TEAID_HIT_MID_HELP = (
    "Hits at least this long but shorter than the full-length threshold use "
    "the mid color. Must stay below the full-length threshold. Click the "
    "color swatch to change it."
)
_TEAID_PROCESS_HELP = (
    "Open a window that shows the TE-Aid command line, then streams genome "
    "BLAST and self-BLAST progress. Useful for long runs."
)


class TeaAidParamsDialog(QDialog):
    """TE-Aid parameters (BLAST-style Inputs / Settings / Advanced / Options).

    Left-click the toolbar button opens this dialog. Parameters holds genome
    BLAST settings; Advanced Settings holds display tuning; Options sets
    right-click and process output.
    """

    _DIALOG_WIDTH = TOOL_PARAMS_DIALOG_WIDTH
    _SETTINGS_LABEL_REF = "Genome BLAST E-value:"

    def __init__(
        self,
        *,
        initial: TeaAidOptions | None = None,
        project: RepeatViewerProject | None = None,
        project_genome: Path | None = None,
        genome_history: list[Path] | None = None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("TE-Aid parameters")
        self.setMinimumWidth(self._DIALOG_WIDTH)
        opts = initial or TeaAidOptions()
        self._for_run = for_run
        self._run_after = False
        # Kept for call-site compatibility; Genome list uses BLAST DBs.
        _ = genome_history
        self._teaid_compact_size: QSize | None = None
        self._x_axis_label = QLineEdit(self)
        self._x_axis_label.setText(
            (opts.x_axis_label or "").strip() or DEFAULT_TEAID_X_AXIS_LABEL
        )
        self._x_axis_label.setPlaceholderText(DEFAULT_TEAID_X_AXIS_LABEL)
        self._x_axis_label.setMinimumWidth(220)

        high_pct, mid_pct = clamp_teaid_hit_len_pct(
            opts.hit_len_high_pct, opts.hit_len_mid_pct
        )

        self._input_group = QButtonGroup(self)
        self._input_selected = QRadioButton("Use selected sequences", self)
        self._input_define = QRadioButton("Define inputs", self)
        self._input_group.addButton(self._input_selected)
        self._input_group.addButton(self._input_define)
        self._input_selected.setChecked(True)
        self._define_inputs = TeMultiInputFilesWidget(
            self,
            initial=None,
            start_empty=False,
            placeholder="Click in Files/Folders or Browse…",
            browse_title="Choose a FASTA file",
            file_filter=(
                "FASTA (*.fa *.fasta *.fna *.fas *.faa *.fsa);;All files (*)"
            ),
            add_button_label="Add another input",
            add_button_tooltip=(
                "Add another FASTA file path. Click a file in the Files panel "
                "to fill a row, or use Browse."
            ),
            path_mode="file",
        )
        self._input_selected.toggled.connect(self._sync_teaid_inputs)
        self._input_define.toggled.connect(self._sync_teaid_inputs)
        self._define_inputs.paths_changed.connect(self._sync_teaid_inputs)

        self._mode_group = QButtonGroup(self)
        self._mode_single = QRadioButton("Single genome", self)
        self._mode_pan = QRadioButton("Pangenome", self)
        self._mode_group.addButton(self._mode_single)
        self._mode_group.addButton(self._mode_pan)
        if opts.use_pan_genome:
            self._mode_pan.setChecked(True)
        else:
            self._mode_single.setChecked(True)

        self._genome_row = GenomeBlastDatabaseRow(
            self,
            project=project,
            project_genome=project_genome,
            selected_path=opts.genome_path,
            manage_title="Manage Database",
        )
        pan_paths: list[str] = []
        if opts.pan_genomes:
            pan_paths = [
                str(path).strip()
                for _name, path in opts.pan_genomes
                if str(path).strip()
            ]
        elif project is not None and project.pan_genomes:
            pan_paths = [str(entry.path) for entry in project.pan_genomes]
        self._pan_list = PanGenomeBlastListWidget(
            self,
            project=project,
            project_genome=project_genome,
            initial_paths=pan_paths,
            manage_title="Manage Database",
        )
        self._pan_list.size_changed.connect(self._on_teaid_pan_list_resized)
        self._mode_single.toggled.connect(self._sync_teaid_mode)
        self._mode_pan.toggled.connect(self._sync_teaid_mode)

        self._evalue = QDoubleSpinBox()
        configure_report_evalue_spinbox(self._evalue)
        self._evalue.setValue(
            normalize_blast_evalue(opts.e_value, default=DEFAULT_TEAID_EVALUE)
        )

        self._word_size = QSpinBox()
        self._word_size.setRange(
            MIN_TEAID_SELF_BLAST_WORD_SIZE, MAX_TEAID_SELF_BLAST_WORD_SIZE
        )
        self._word_size.setValue(
            max(
                MIN_TEAID_SELF_BLAST_WORD_SIZE,
                min(MAX_TEAID_SELF_BLAST_WORD_SIZE, opts.self_blast_word_size),
            )
        )

        self._threads = QSpinBox()
        self._threads.setRange(1, 128)
        self._threads.setValue(max(1, int(opts.num_threads)))

        self._line_transparency = QDoubleSpinBox()
        self._line_transparency.setRange(0.0, 100.0)
        self._line_transparency.setDecimals(0)
        self._line_transparency.setSingleStep(5.0)
        self._line_transparency.setSuffix(" %")
        self._line_transparency.setValue(
            clamp_teaid_hit_line_transparency(opts.hit_line_transparency) * 100.0
        )

        self._hit_high_pct = QDoubleSpinBox()
        self._hit_high_pct.setRange(0.0, 100.0)
        self._hit_high_pct.setDecimals(1)
        self._hit_high_pct.setSingleStep(1.0)
        self._hit_high_pct.setSuffix(" %")
        self._hit_high_pct.setValue(high_pct)
        self._hit_color_high = opts.hit_color_high or DEFAULT_TEAID_HIT_COLOR_HIGH
        self._hit_color_high_btn = QPushButton()
        self._hit_color_high_btn.setFixedSize(56, 22)
        self._hit_color_high_btn.setToolTip("Color for full-length hits")
        self._hit_color_high_btn.clicked.connect(self._pick_hit_color_high)
        self._update_hit_color_button(self._hit_color_high_btn, self._hit_color_high)

        self._hit_mid_pct = QDoubleSpinBox()
        self._hit_mid_pct.setRange(0.0, 100.0)
        self._hit_mid_pct.setDecimals(1)
        self._hit_mid_pct.setSingleStep(1.0)
        self._hit_mid_pct.setSuffix(" %")
        self._hit_mid_pct.setValue(mid_pct)
        self._hit_color_mid = opts.hit_color_mid or DEFAULT_TEAID_HIT_COLOR_MID
        self._hit_color_mid_btn = QPushButton()
        self._hit_color_mid_btn.setFixedSize(56, 22)
        self._hit_color_mid_btn.setToolTip("Color for mid-length hits")
        self._hit_color_mid_btn.clicked.connect(self._pick_hit_color_mid)
        self._update_hit_color_button(self._hit_color_mid_btn, self._hit_color_mid)

        self._hit_high_pct.valueChanged.connect(self._on_hit_high_changed)
        self._hit_mid_pct.valueChanged.connect(self._on_hit_mid_changed)

        for field in (
            self._evalue,
            self._word_size,
            self._threads,
            self._line_transparency,
        ):
            field.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
            if hasattr(field, "setAlignment"):
                field.setAlignment(Qt.AlignmentFlag.AlignRight)

        # Hit rows: compact % + color + suffix (previous style; not full 260).
        for field in (self._hit_high_pct, self._hit_mid_pct):
            field.setFixedWidth(100)
            field.setAlignment(Qt.AlignmentFlag.AlignRight)

        high_row = QWidget(self)
        high_layout = QHBoxLayout(high_row)
        high_layout.setContentsMargins(0, 0, 0, 0)
        high_layout.setSpacing(8)
        high_layout.addWidget(self._hit_high_pct)
        high_layout.addWidget(self._hit_color_high_btn)
        high_layout.addWidget(QLabel("of sequence length", self))
        high_layout.addStretch(1)

        mid_row = QWidget(self)
        mid_layout = QHBoxLayout(mid_row)
        mid_layout.setContentsMargins(0, 0, 0, 0)
        mid_layout.setSpacing(8)
        mid_layout.addWidget(self._hit_mid_pct)
        mid_layout.addWidget(self._hit_color_mid_btn)
        mid_layout.addWidget(QLabel("but below full-length", self))
        mid_layout.addStretch(1)

        params_tab = QWidget(self)
        params_root = QVBoxLayout(params_tab)
        params_root.setContentsMargins(16, 10, 16, 10)
        params_root.setSpacing(4)
        basic = QWidget(params_tab)
        self._teaid_basic = basic
        basic.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )
        form = QFormLayout(basic)
        form.setContentsMargins(0, 0, 0, 0)
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        # Expanding so Define-inputs path + Browse fill the dialog width.
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)

        def _section(title: str) -> QLabel:
            return QLabel(f"<b>{title}</b>", params_tab)

        form.addRow(_section("Mode"))
        mode_box = QWidget(basic)
        mode_layout = QVBoxLayout(mode_box)
        mode_layout.setContentsMargins(0, 0, 0, 0)
        mode_layout.setSpacing(2)
        mode_layout.addWidget(
            field_with_help(
                self._mode_single,
                "TE-Aid mode",
                _TEAID_MODE_HELP,
                stretch=True,
            )
        )
        mode_layout.addWidget(
            field_with_help(
                self._mode_pan,
                "Pangenome",
                _TEAID_PANGENOME_HELP,
                stretch=True,
            )
        )
        self._mode_radios_label = QLabel("")
        self._mode_single_genome_label = QLabel("")
        self._mode_pan_list_label = QLabel("")
        form.addRow(self._mode_radios_label, mode_box)
        form.addRow(self._mode_single_genome_label, self._genome_row)
        form.addRow(self._mode_pan_list_label, self._pan_list)

        form.addRow(_form_section_separator())
        form.addRow(_section("Inputs"))
        inputs_box = QWidget(basic)
        inputs_layout = QVBoxLayout(inputs_box)
        inputs_layout.setContentsMargins(0, 0, 0, 0)
        inputs_layout.setSpacing(2)
        inputs_layout.addWidget(
            field_with_help(
                self._input_selected,
                "TE-Aid inputs",
                _TEAID_INPUTS_HELP,
                stretch=True,
            )
        )
        inputs_layout.addWidget(self._input_define)
        self._inputs_radios_label = QLabel("")
        self._inputs_define_label = QLabel("")
        form.addRow(self._inputs_radios_label, inputs_box)
        form.addRow(self._inputs_define_label, self._define_inputs)

        form.addRow(_form_section_separator())
        form.addRow(_section("Settings"))
        form.addRow(
            "Genome BLAST E-value:",
            field_with_help(
                self._evalue, "Genome BLAST E-value", _TEAID_EVALUE_HELP
            ),
        )
        form.addRow(
            "Self BLAST word size:",
            field_with_help(
                self._word_size, "Self BLAST word size", _TEAID_WORD_SIZE_HELP
            ),
        )
        form.addRow(
            "Threads:",
            field_with_help(self._threads, "Threads", _TEAID_THREADS_HELP),
        )
        self._teaid_basic_form = form
        params_root.addWidget(basic)

        advanced = QWidget(params_tab)
        self._teaid_advanced_panel = advanced
        advanced.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        advanced_form = QFormLayout(advanced)
        advanced_form.setContentsMargins(0, 4, 0, 0)
        advanced_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        advanced_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        advanced_form.setHorizontalSpacing(12)
        advanced_form.setVerticalSpacing(4)
        advanced_form.addRow(_form_section_separator())
        advanced_form.addRow(_section("Advanced Settings"))
        advanced_form.addRow(
            "Line transparency:",
            field_with_help(
                self._line_transparency,
                "Line transparency",
                _TEAID_TRANSPARENCY_HELP,
            ),
        )
        advanced_form.addRow(
            "Full-length hits (≥):",
            field_with_help(
                high_row, "Full-length hits", _TEAID_HIT_HIGH_HELP
            ),
        )
        advanced_form.addRow(
            "Mid-length hits (≥):",
            field_with_help(mid_row, "Mid-length hits", _TEAID_HIT_MID_HELP),
        )
        advanced_form.addRow(
            "X-axis label:",
            field_with_help(
                self._x_axis_label, "X-axis label", _TEAID_X_AXIS_LABEL_HELP
            ),
        )
        params_root.addWidget(advanced)
        params_root.addStretch(1)
        self._params_tab = params_tab

        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("Run TE-Aid with saved parameters", TEAID_RIGHT_CLICK_RUN),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[TEAID_RIGHT_CLICK_RUN].setChecked(True)

        rc_box = QGroupBox("Right-click action", self)
        rc_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(QLabel("Right-clicking the TE-Aid button does:", self))
        for radio in self._right_click_buttons.values():
            rc_box_layout.addWidget(radio)

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="TE-Aid",
                checked=bool(getattr(opts, "show_process_window", True)),
                close_when_finished=bool(
                    getattr(opts, "close_process_window_when_finished", True)
                ),
            )
        )
        self._show_process_window.setToolTip(_TEAID_PROCESS_HELP)

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)
        options_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        options_layout.addWidget(rc_box)
        options_layout.addWidget(process_box)
        from repeatviewer.tool_citations import CITATION_TEAID

        options_layout.addWidget(make_citation_box(self, CITATION_TEAID))
        options_layout.addStretch(1)
        options_tab.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._options_tab = options_tab

        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        self._tabs.tabBar().setDrawBase(False)
        self._tabs.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._tabs.addTab(params_tab, "Parameters")
        self._tabs.addTab(options_tab, "Options")
        self._tabs.currentChanged.connect(self._on_teaid_tab_changed)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._advanced_settings_btn = QPushButton("Advanced Settings ▸", self)
        self._advanced_settings_btn.setCheckable(True)
        self._advanced_settings_btn.setChecked(False)
        self._advanced_settings_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._advanced_settings_btn.toggled.connect(self._on_advanced_settings_toggled)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore, self._advanced_settings_btn],
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(14)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._tabs, 0)
        layout.addWidget(self._button_box, 0)

        self._sync_teaid_inputs(resize=False)
        self._sync_teaid_mode(resize=False)
        self._set_teaid_advanced_visible(False)
        self.resize(self._DIALOG_WIDTH, 420)
        QTimer.singleShot(0, lambda: self._resize_for_teaid_advanced(False))

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._teaid_compact_size = None
        self._sync_field_column()
        QTimer.singleShot(0, lambda: self._resize_for_teaid_advanced(False))

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _on_teaid_pan_list_resized(self) -> None:
        if self.isVisible() and self._mode_pan.isChecked():
            self._teaid_compact_size = None
            # Defer so the new genome row's sizeHint is ready before we pin height.
            QTimer.singleShot(0, self._resize_teaid_after_pan_change)

    def _resize_teaid_after_pan_change(self) -> None:
        if not self.isVisible() or not self._mode_pan.isChecked():
            return
        self._pan_list.updateGeometry()
        self._teaid_basic.updateGeometry()
        self._resize_for_teaid_advanced(self._advanced_settings_btn.isChecked())

    def _schedule_teaid_refit(self) -> None:
        """Refit the Parameters tab after rows appear/disappear.

        Also while Advanced Settings is open: the tab height is pinned, so
        without a refit the new rows (Pangenome list, Define inputs) get
        squeezed into the old height and their buttons overlap.
        Deferred so the form's new sizeHint is ready.
        """
        self._teaid_compact_size = None

        def _refit() -> None:
            if not self.isVisible() or self._tabs.currentWidget() is not self._params_tab:
                return
            self._teaid_basic.updateGeometry()
            self._resize_for_teaid_advanced(self._advanced_settings_btn.isChecked())

        QTimer.singleShot(0, _refit)

    def _sync_teaid_mode(self, *, resize: bool = True) -> None:
        pan = self._mode_pan.isChecked()
        form = self._teaid_basic.layout()
        if isinstance(form, QFormLayout):
            form.setRowVisible(self._genome_row, not pan)
            form.setRowVisible(self._pan_list, pan)
        self._genome_row.setVisible(not pan)
        self._pan_list.setVisible(pan)
        if pan:
            self._pan_list.setMinimumHeight(0)
            self._pan_list.setMaximumHeight(16777215)
            self._genome_row.setFixedHeight(0)
        else:
            self._genome_row.setMinimumHeight(0)
            self._genome_row.setMaximumHeight(16777215)
            self._pan_list.setFixedHeight(0)
        self._genome_row.updateGeometry()
        self._pan_list.updateGeometry()
        self._teaid_basic.updateGeometry()
        self._sync_field_column()
        if resize and self.isVisible():
            self._schedule_teaid_refit()

    def input_source(self) -> str:
        if self._input_define.isChecked():
            return ALIGN_INPUT_DEFINE
        return ALIGN_INPUT_SELECTED

    def run_files(self, *, panel_files: list[Path] | None = None) -> list[Path]:
        from repeatviewer.backends.fastga_dotplot import expand_fastga_input_paths

        if self.input_source() == ALIGN_INPUT_DEFINE:
            return expand_fastga_input_paths(self._define_inputs.paths())
        return [Path(p) for p in (panel_files or [])]

    def apply_path_from_panel(self, path: Path) -> bool:
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        if self._tabs.currentWidget() is self._options_tab:
            return False
        try:
            candidate = Path(path).expanduser().resolve()
        except OSError:
            return False
        return bool(self._define_inputs.accept_path(candidate))

    def apply_folder_from_panel(self, folder: Path) -> bool:
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        if self._tabs.currentWidget() is self._options_tab:
            return False
        try:
            candidate = Path(folder).expanduser().resolve()
        except OSError:
            return False
        if not candidate.is_dir():
            return False
        return bool(self._define_inputs.accept_path(candidate))

    def _sync_teaid_inputs(self, *, resize: bool = True) -> None:
        define = self.input_source() == ALIGN_INPUT_DEFINE
        form = self._teaid_basic.layout()
        if isinstance(form, QFormLayout):
            form.setRowVisible(self._define_inputs, define)
        self._define_inputs.setVisible(define)
        if define:
            self._define_inputs.setMinimumHeight(0)
            self._define_inputs.setMaximumHeight(16777215)
        else:
            self._define_inputs.setFixedHeight(0)
        self._define_inputs.updateGeometry()
        self._teaid_basic.updateGeometry()
        self._sync_field_column()
        if resize and self.isVisible():
            self._schedule_teaid_refit()

    def _sync_field_column(self) -> None:
        if not hasattr(self, "_inputs_radios_label"):
            return
        adv = (
            self._teaid_advanced_panel.layout()
            if hasattr(self, "_teaid_advanced_panel")
            else None
        )
        sync_form_field_column_spacers(
            self,
            [
                lab
                for lab in (
                    getattr(self, "_mode_radios_label", None),
                    getattr(self, "_mode_single_genome_label", None),
                    getattr(self, "_mode_pan_list_label", None),
                    self._inputs_radios_label,
                    self._inputs_define_label,
                )
                if lab is not None
            ],
            [
                getattr(self, "_teaid_basic_form", None),
                adv if isinstance(adv, QFormLayout) else None,
            ],
            ref_label=self._SETTINGS_LABEL_REF,
        )

    def _on_teaid_tab_changed(self, _index: int) -> None:
        on_params = self._tabs.currentWidget() is self._params_tab
        self._advanced_settings_btn.setVisible(on_params)
        if not on_params:
            if self._advanced_settings_btn.isChecked():
                self._advanced_settings_btn.blockSignals(True)
                self._advanced_settings_btn.setChecked(False)
                self._advanced_settings_btn.setText("Advanced Settings ▸")
                self._advanced_settings_btn.blockSignals(False)
                self._set_teaid_advanced_visible(False)
            if self.isVisible():
                self._fit_options_content()
            return
        if self.isVisible():
            self._resize_for_teaid_advanced(False)

    def _set_teaid_advanced_visible(self, visible: bool) -> None:
        panel = self._teaid_advanced_panel
        panel.setVisible(visible)
        if visible:
            panel.setMaximumHeight(16777215)
        else:
            panel.setMaximumHeight(0)
        panel.updateGeometry()
        root = self._params_tab.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        self._params_tab.updateGeometry()
        self._tabs.updateGeometry()

    def _on_advanced_settings_toggled(self, checked: bool) -> None:
        self._advanced_settings_btn.setText(
            "Advanced Settings ▾" if checked else "Advanced Settings ▸"
        )
        if checked:
            self._teaid_compact_size = self.size()
        self._set_teaid_advanced_visible(checked)
        self._resize_for_teaid_advanced(checked)

    def _dialog_chrome_height(self) -> int:
        lay = self.layout()
        margins = lay.contentsMargins() if lay is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 34
        spacing = lay.spacing() if lay is not None else 14
        tab_bar_h = self._tabs.tabBar().sizeHint().height() + 10
        buttons_h = self._button_box.sizeHint().height()
        return margin_h + spacing + tab_bar_h + buttons_h

    def _teaid_compact_content_height(self) -> int:
        root = self._params_tab.layout()
        margins = root.contentsMargins() if root is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 20
        spacing = root.spacing() if root is not None else 4
        basic_h = max(
            self._teaid_basic.minimumSizeHint().height(),
            self._teaid_basic.sizeHint().height(),
        )
        return margin_h + basic_h + spacing

    def _resize_for_teaid_advanced(self, expanded: bool) -> None:
        if self._tabs.currentWidget() is not self._params_tab:
            self._fit_options_content()
            return
        panel = self._params_tab
        adv = self._teaid_advanced_panel
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        root = panel.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        panel.updateGeometry()
        self._tabs.updateGeometry()

        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        chrome = self._dialog_chrome_height()
        width = tool_params_dialog_width(
            define_inputs=self.input_source() == ALIGN_INPUT_DEFINE
        )
        self.setMinimumWidth(width)
        compact_h = self._teaid_compact_content_height()

        if not expanded:
            height = chrome + compact_h
            height = min(max(height, 220), max_h)
            pin_tool_tabs_height(self._tabs, compact_h)
            self.setMaximumHeight(height)
            self.resize(width, height)
            self.setMaximumHeight(16777215)
            self._teaid_compact_size = QSize(width, height)
            return

        extra = max(0, adv.sizeHint().height())
        page_h = compact_h + extra
        base = self._teaid_compact_size or QSize(width, chrome + compact_h)
        height = min(base.height() + extra + 12, max_h)
        height = max(height, min(chrome + page_h, max_h))
        pin_tool_tabs_height(self._tabs, page_h)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def _fit_options_content(self) -> None:
        page = self._options_tab
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        page.ensurePolished()
        page.adjustSize()
        root = page.layout()
        if root is not None:
            root.activate()
        page_hint = page.minimumSizeHint().expandedTo(page.sizeHint())
        content_h = page_hint.height()
        if root is not None:
            content_h = max(
                content_h,
                root.totalMinimumSize().height(),
                root.totalSizeHint().height(),
            )
        chrome = self._dialog_chrome_height()
        height = chrome + content_h
        pin_tool_tabs_height(self._tabs, content_h)
        self.setMaximumHeight(height)
        self.resize(self._DIALOG_WIDTH, height)
        self.setMaximumHeight(16777215)

    def _fit_to_current_tab(self) -> None:
        if self._tabs.currentWidget() is self._options_tab:
            self._fit_options_content()
            return
        self._resize_for_teaid_advanced(self._advanced_settings_btn.isChecked())
    @staticmethod
    def _update_hit_color_button(button: QPushButton, color: str) -> None:
        hex_color = QColor(color).name(QColor.NameFormat.HexRgb)
        button.setStyleSheet(f"background-color: {hex_color};")

    def _pick_hit_color_high(self) -> None:
        chosen = QColorDialog.getColor(
            QColor(self._hit_color_high), self, "Hit color (full-length)"
        )
        if not chosen.isValid():
            return
        self._hit_color_high = chosen.name(QColor.NameFormat.HexRgb)
        self._update_hit_color_button(self._hit_color_high_btn, self._hit_color_high)

    def _pick_hit_color_mid(self) -> None:
        chosen = QColorDialog.getColor(
            QColor(self._hit_color_mid), self, "Hit color (mid-length)"
        )
        if not chosen.isValid():
            return
        self._hit_color_mid = chosen.name(QColor.NameFormat.HexRgb)
        self._update_hit_color_button(self._hit_color_mid_btn, self._hit_color_mid)

    def _on_hit_high_changed(self, value: float) -> None:
        high_pct, mid_pct = clamp_teaid_hit_len_pct(value, self._hit_mid_pct.value())
        if abs(self._hit_mid_pct.value() - mid_pct) > 1e-9:
            self._hit_mid_pct.blockSignals(True)
            self._hit_mid_pct.setValue(mid_pct)
            self._hit_mid_pct.blockSignals(False)
        if abs(self._hit_high_pct.value() - high_pct) > 1e-9:
            self._hit_high_pct.blockSignals(True)
            self._hit_high_pct.setValue(high_pct)
            self._hit_high_pct.blockSignals(False)

    def _on_hit_mid_changed(self, value: float) -> None:
        high_pct, mid_pct = clamp_teaid_hit_len_pct(self._hit_high_pct.value(), value)
        if abs(self._hit_mid_pct.value() - mid_pct) > 1e-9:
            self._hit_mid_pct.blockSignals(True)
            self._hit_mid_pct.setValue(mid_pct)
            self._hit_mid_pct.blockSignals(False)

    def _restore_defaults(self) -> None:
        defaults = TeaAidOptions()
        self._mode_single.setChecked(True)
        self._input_selected.setChecked(True)
        self._define_inputs.clear_paths()
        self._sync_teaid_inputs()
        self._sync_teaid_mode()
        self._genome_row.reset_to_default()
        self._pan_list.reset_to_defaults()
        self._evalue.setValue(DEFAULT_TEAID_EVALUE)
        self._word_size.setValue(DEFAULT_TEAID_SELF_BLAST_WORD_SIZE)
        self._threads.setValue(max(1, defaults.num_threads))
        self._hit_high_pct.blockSignals(True)
        self._hit_mid_pct.blockSignals(True)
        self._hit_high_pct.setValue(DEFAULT_TEAID_HIT_LEN_HIGH_PCT)
        self._hit_mid_pct.setValue(DEFAULT_TEAID_HIT_LEN_MID_PCT)
        self._hit_high_pct.blockSignals(False)
        self._hit_mid_pct.blockSignals(False)
        self._hit_color_high = DEFAULT_TEAID_HIT_COLOR_HIGH
        self._hit_color_mid = DEFAULT_TEAID_HIT_COLOR_MID
        self._update_hit_color_button(self._hit_color_high_btn, self._hit_color_high)
        self._update_hit_color_button(self._hit_color_mid_btn, self._hit_color_mid)
        self._line_transparency.setValue(DEFAULT_TEAID_HIT_LINE_TRANSPARENCY * 100.0)
        self._x_axis_label.setText(DEFAULT_TEAID_X_AXIS_LABEL)
        self._right_click_buttons[TEAID_RIGHT_CLICK_RUN].setChecked(True)
        self._show_process_window.setChecked(defaults.show_process_window)
        self._close_process_window.setChecked(
            defaults.close_process_window_when_finished
        )

    def options(self) -> TeaAidOptions:
        high_pct, mid_pct = clamp_teaid_hit_len_pct(
            self._hit_high_pct.value(), self._hit_mid_pct.value()
        )
        right_click_action = TEAID_RIGHT_CLICK_RUN
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click_action = value
                break
        return TeaAidOptions(
            e_value=normalize_blast_evalue(
                self._evalue.value(),
                default=DEFAULT_TEAID_EVALUE,
            ),
            self_blast_word_size=self._word_size.value(),
            num_threads=max(1, self._threads.value()),
            x_axis_label=self._x_axis_label.text().strip()
            or DEFAULT_TEAID_X_AXIS_LABEL,
            hit_len_high_pct=high_pct,
            hit_len_mid_pct=mid_pct,
            hit_color_high=self._hit_color_high or DEFAULT_TEAID_HIT_COLOR_HIGH,
            hit_color_mid=self._hit_color_mid or DEFAULT_TEAID_HIT_COLOR_MID,
            hit_line_transparency=clamp_teaid_hit_line_transparency(
                self._line_transparency.value() / 100.0
            ),
            right_click_action=right_click_action,
            genome_path=self._genome_row.selected_path(),
            use_pan_genome=self._mode_pan.isChecked(),
            pan_genomes=[
                (name.strip(), path.strip())
                for name, path in self._pan_list.snapshot_rows()
                if name.strip() and path.strip()
            ]
            if self._mode_pan.isChecked()
            else [],
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=(
                self._close_process_window.isChecked()
            ),
        )


class MsaCleanParamsDialog(QDialog):
    """MSA clean parameters (BLAST-style strategy tabs).

    Tabs: Remove gappy columns | Clean with divergence | Clean with gaps |
    Options. Each strategy tab has Inputs + Settings. Combined pipeline
    (which steps, in what order) lives under Options with right-click.
    """

    _DIALOG_WIDTH = TOOL_PARAMS_DIALOG_WIDTH
    _TAB_GAPPY = 0
    _TAB_DIV = 1
    _TAB_GAP = 2
    _TAB_OPTIONS = 3

    _INPUTS_HELP = (
        "<b>Use selected sequences</b> (default): clean the aligned FASTA files "
        "selected in the Files panel.<br><br>"
        "<b>Define inputs</b>: choose alignment files explicitly "
        "(click in Files or Browse). Each file is cleaned separately."
    )
    _GAP_THRESHOLD_HELP = (
        "Columns whose gap proportion is <b>above</b> this value are removed "
        "immediately.<br><br>"
        "Example: with threshold <b>0.80</b>, a column that is 85% gaps is dropped."
    )
    _SIMI_GAP_HELP = (
        "For columns that are <b>not</b> removed by the gap threshold, this is "
        "the maximum gap proportion allowed when applying the similarity keep "
        "rule.<br><br>"
        "Example: with <b>0.40</b>, only columns with fewer than 40% gaps are "
        "eligible to be kept based on similarity."
    )
    _SIMILARITY_HELP = (
        "Among columns that pass the similarity-check gap threshold, keep the "
        "column only if the most common nucleotide occupies at least this "
        "fraction of non-gap bases.<br><br>"
        "Example: with <b>0.70</b>, a column that is mostly T (≥70%) is kept; "
        "a mixed column is removed."
    )
    _MIN_NT_HELP = (
        "Columns with fewer than this many non-gap nucleotides are removed "
        "(too little data to trust).<br><br>"
        "Example: with <b>5</b>, a column that has only 3 real bases is dropped."
    )
    _DIV_THRESHOLD_HELP = (
        "End-cropping continues while the sum of per-base divergence scores "
        "in the sliding window is <b>below</b> this value. Cropped bases become "
        "gaps (<code>-</code>).<br><br>"
        "Higher values crop more aggressively from each sequence end."
    )
    _DIV_WINDOW_HELP = (
        "Length (bp) of the sliding window used at each sequence end when "
        "scoring divergence for cropping.<br><br>"
        "Example: <b>20</b> means each step looks at 20 positions from the end."
    )
    _ROW_GAP_THRESHOLD_HELP = (
        "End-cropping continues while the gap proportion inside the sliding "
        "window is <b>above</b> this value. Cropped bases become gaps "
        "(<code>-</code>).<br><br>"
        "Example: with <b>0.050</b>, cropping stops once the window is under "
        "5% gaps."
    )
    _ROW_GAP_WINDOW_HELP = (
        "Length (bp) of the sliding window used at each sequence end when "
        "measuring gap proportion for cropping.<br><br>"
        "Example: <b>300</b> means each step inspects 300 positions from the end."
    )
    _COMBINED_HELP = (
        "Combined MSA clean runs the <b>checked</b> steps below in list order "
        "on each input alignment.<br><br>"
        "• Check a step to include it; uncheck to skip it.<br>"
        "• Drag steps to change run order.<br>"
        "• At least one step must stay checked.<br><br>"
        "<b>Important:</b> each step cleans the <b>result of the previous "
        "step</b>, not the original file again. The pipeline is cumulative.<br><br>"
        "<b>Example:</b> with this order checked "
        "(1) <b>Remove gappy columns</b>, "
        "(2) <b>Clean rows with divergence</b>, "
        "(3) <b>Clean rows with gaps</b>:<br>"
        "first gappy columns are removed (using that tab’s Settings); "
        "then divergence cleaning runs on that cleaned alignment; "
        "then gap cleaning runs on the divergence cleaned result. "
        "If you uncheck divergence, gap cleaning still runs on the output of "
        "gappy column removal.<br><br>"
        "Each step uses the Settings from its own tab. Combined is used when "
        "you choose it under Right click action, or when you Run from the "
        "Options tab."
    )

    def __init__(
        self,
        *,
        initial: MsaCleanOptions | None = None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("MSA clean parameters")
        self.setMinimumWidth(self._DIALOG_WIDTH)
        opts = initial or MsaCleanOptions()
        self._for_run = for_run
        self._run_after = False
        self._msa_compact_size: QSize | None = None

        self._method_labels = {
            "msa_clean_div": "Clean rows with divergence",
            "msa_clean_gap": "Clean rows with gaps",
            "msa_clean_gappy": "Remove gappy columns",
        }
        enabled_methods = {
            "msa_clean_div": opts.use_divergence,
            "msa_clean_gap": opts.use_gap_row,
            "msa_clean_gappy": opts.use_gappy_columns,
        }

        def _make_form(parent: QWidget) -> QFormLayout:
            form = QFormLayout(parent)
            form.setContentsMargins(0, 0, 0, 0)
            form.setLabelAlignment(
                Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
            )
            form.setFieldGrowthPolicy(
                QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
            )
            form.setHorizontalSpacing(12)
            form.setVerticalSpacing(4)
            return form

        def _wrap_page(page: QWidget) -> QWidget:
            wrap = QWidget(self)
            wrap.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
            layout = QVBoxLayout(wrap)
            layout.setContentsMargins(16, 10, 16, 10)
            layout.setSpacing(4)
            layout.setAlignment(Qt.AlignmentFlag.AlignTop)
            layout.addWidget(page)
            return wrap

        def _spin_float(
            value: float, *, decimals: int, step: float, lo: float = 0.0, hi: float = 1.0
        ) -> QDoubleSpinBox:
            spin = QDoubleSpinBox(self)
            spin.setRange(lo, hi)
            spin.setDecimals(decimals)
            spin.setSingleStep(step)
            spin.setValue(value)
            spin.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
            return spin

        def _spin_int(value: int, *, lo: int, hi: int) -> QSpinBox:
            spin = QSpinBox(self)
            spin.setRange(lo, hi)
            spin.setValue(value)
            spin.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
            return spin

        def _section(title: str, parent: QWidget | None = None) -> QLabel:
            return QLabel(f"<b>{title}</b>", parent)

        def _make_define() -> TeMultiInputFilesWidget:
            return TeMultiInputFilesWidget(
                self,
                initial=None,
                allow_multiple=True,
                start_empty=False,
                placeholder="Click in Files/Folders or Browse…",
                browse_title="Select MSA FASTA",
                file_filter="FASTA (*.fa *.fasta *.fas *.fna *.faa);;All files (*)",
                add_button_label="Add another input",
                add_button_tooltip=(
                    "Add another aligned FASTA file path. Click a file in the "
                    "Files panel to fill a row, or use Browse."
                ),
                path_mode="file",
            )

        def _add_inputs_settings(
            form: QFormLayout,
            page: QWidget,
            *,
            method_id: str,
            selected: QRadioButton,
            define: QRadioButton,
            define_inputs: TeMultiInputFilesWidget,
            sel_spacer: QLabel,
            def_spacer: QLabel,
            settings_rows: list[tuple[str, QWidget, str]],
        ) -> None:
            del method_id  # strategy overview lives on per-parameter ? helps
            # BLAST / Align Inputs pattern: both radios in one field column;
            # Define path rows on the next form row with an empty label so Browse
            # lines up with Settings (never addRow(define_widget) spanning).
            form.addRow(_section("Inputs", page))
            inputs_box = QWidget(page)
            inputs_layout = QVBoxLayout(inputs_box)
            inputs_layout.setContentsMargins(0, 0, 0, 0)
            inputs_layout.setSpacing(2)
            inputs_layout.addWidget(
                field_with_help(
                    selected, "MSA clean inputs", self._INPUTS_HELP, stretch=True
                )
            )
            inputs_layout.addWidget(define)
            form.addRow(sel_spacer, inputs_box)
            form.addRow(def_spacer, define_inputs)
            form.addRow(_form_section_separator())
            form.addRow(_section("Settings", page))
            for label, widget, help_text in settings_rows:
                form.addRow(
                    label,
                    field_with_help(widget, label.rstrip(":"), help_text),
                )

        # --- Gappy tab ------------------------------------------------------
        self._gap_threshold = _spin_float(opts.gappy.gap_threshold, decimals=2, step=0.05)
        self._simi_gap = _spin_float(opts.gappy.simi_check_gap_thr, decimals=2, step=0.05)
        self._similarity = _spin_float(opts.gappy.similarity_thr, decimals=2, step=0.05)
        self._min_nt = _spin_int(opts.gappy.min_nucleotide, lo=1, hi=50)
        gappy_page = QWidget(self)
        gappy_page.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        self._gappy_form = _make_form(gappy_page)
        self._gappy_input_selected = QRadioButton("Use selected sequences", gappy_page)
        self._gappy_input_define = QRadioButton("Define inputs", gappy_page)
        gappy_group = QButtonGroup(gappy_page)
        gappy_group.addButton(self._gappy_input_selected)
        gappy_group.addButton(self._gappy_input_define)
        self._gappy_input_selected.setChecked(True)
        self._gappy_define_inputs = _make_define()
        self._gappy_sel_spacer = QLabel("")
        self._gappy_def_spacer = QLabel("")
        _add_inputs_settings(
            self._gappy_form,
            gappy_page,
            method_id="msa_clean_gappy",
            selected=self._gappy_input_selected,
            define=self._gappy_input_define,
            define_inputs=self._gappy_define_inputs,
            sel_spacer=self._gappy_sel_spacer,
            def_spacer=self._gappy_def_spacer,
            settings_rows=[
                ("Gap threshold:", self._gap_threshold, self._GAP_THRESHOLD_HELP),
                (
                    "Similarity check gap threshold:",
                    self._simi_gap,
                    self._SIMI_GAP_HELP,
                ),
                ("Similarity threshold:", self._similarity, self._SIMILARITY_HELP),
                ("Min nucleotides per column:", self._min_nt, self._MIN_NT_HELP),
            ],
        )
        self._gappy_tab = _wrap_page(gappy_page)

        # --- Divergence tab -------------------------------------------------
        self._div_threshold = _spin_float(
            opts.divergence.threshold, decimals=2, step=0.05
        )
        self._div_window = _spin_int(opts.divergence.window_size, lo=1, hi=10_000)
        div_page = QWidget(self)
        div_page.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        self._div_form = _make_form(div_page)
        self._div_input_selected = QRadioButton("Use selected sequences", div_page)
        self._div_input_define = QRadioButton("Define inputs", div_page)
        div_group = QButtonGroup(div_page)
        div_group.addButton(self._div_input_selected)
        div_group.addButton(self._div_input_define)
        self._div_input_selected.setChecked(True)
        self._div_define_inputs = _make_define()
        self._div_sel_spacer = QLabel("")
        self._div_def_spacer = QLabel("")
        _add_inputs_settings(
            self._div_form,
            div_page,
            method_id="msa_clean_div",
            selected=self._div_input_selected,
            define=self._div_input_define,
            define_inputs=self._div_define_inputs,
            sel_spacer=self._div_sel_spacer,
            def_spacer=self._div_def_spacer,
            settings_rows=[
                ("Divergence threshold:", self._div_threshold, self._DIV_THRESHOLD_HELP),
                ("Window size (bp):", self._div_window, self._DIV_WINDOW_HELP),
            ],
        )
        self._div_tab = _wrap_page(div_page)

        # --- Gaps tab -------------------------------------------------------
        self._row_gap_threshold = _spin_float(
            opts.gap_row.gap_threshold, decimals=3, step=0.01
        )
        self._row_gap_window = _spin_int(opts.gap_row.window_size, lo=1, hi=10_000)
        gap_page = QWidget(self)
        gap_page.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        self._gap_form = _make_form(gap_page)
        self._gap_input_selected = QRadioButton("Use selected sequences", gap_page)
        self._gap_input_define = QRadioButton("Define inputs", gap_page)
        gap_group = QButtonGroup(gap_page)
        gap_group.addButton(self._gap_input_selected)
        gap_group.addButton(self._gap_input_define)
        self._gap_input_selected.setChecked(True)
        self._gap_define_inputs = _make_define()
        self._gap_sel_spacer = QLabel("")
        self._gap_def_spacer = QLabel("")
        _add_inputs_settings(
            self._gap_form,
            gap_page,
            method_id="msa_clean_gap",
            selected=self._gap_input_selected,
            define=self._gap_input_define,
            define_inputs=self._gap_define_inputs,
            sel_spacer=self._gap_sel_spacer,
            def_spacer=self._gap_def_spacer,
            settings_rows=[
                (
                    "Gap proportion threshold:",
                    self._row_gap_threshold,
                    self._ROW_GAP_THRESHOLD_HELP,
                ),
                ("Window size (bp):", self._row_gap_window, self._ROW_GAP_WINDOW_HELP),
            ],
        )
        self._gap_tab = _wrap_page(gap_page)

        self._focused_file_input: TeMultiInputFilesWidget | None = None
        for widget in (
            self._gappy_define_inputs,
            self._div_define_inputs,
            self._gap_define_inputs,
        ):
            widget.field_activated.connect(self._on_file_input_field_activated)
            widget.paths_changed.connect(self._on_define_paths_changed)

        for selected, define in (
            (self._gappy_input_selected, self._gappy_input_define),
            (self._div_input_selected, self._div_input_define),
            (self._gap_input_selected, self._gap_input_define),
        ):
            selected.toggled.connect(self._sync_active_inputs)
            define.toggled.connect(self._sync_active_inputs)

        # --- Combined checklist (Options) -----------------------------------
        self._method_order = QListWidget()
        self._method_order.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self._method_order.setDefaultDropAction(Qt.DropAction.MoveAction)
        self._method_order.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self._method_order.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        self._method_order.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._method_order.setStyleSheet(
            "QListWidget { padding-top: 4px; padding-bottom: 4px; }\n"
            + item_view_checkbox_stylesheet()
        )
        self._method_order.setToolTip(
            "Check steps to include in Combined MSA clean. Drag to set execution order."
        )
        for method_id in opts.method_order:
            self._add_method_list_item(
                method_id,
                checked=enabled_methods.get(method_id, False),
            )
        listed = {
            str(self._method_order.item(i).data(Qt.ItemDataRole.UserRole))
            for i in range(self._method_order.count())
        }
        for method_id in DEFAULT_MSA_CLEAN_ORDER:
            if method_id in listed:
                continue
            self._add_method_list_item(
                method_id,
                checked=enabled_methods.get(method_id, False),
            )
        self._method_order.itemChanged.connect(self._on_method_item_changed)
        self._fit_method_list_height()

        combined_box = QGroupBox("Combined MSA clean", self)
        combined_layout = QVBoxLayout(combined_box)
        combined_layout.setContentsMargins(10, 8, 10, 8)
        combined_layout.setSpacing(6)
        combined_header = QWidget(combined_box)
        combined_header_layout = QHBoxLayout(combined_header)
        combined_header_layout.setContentsMargins(0, 0, 0, 0)
        combined_header_layout.setSpacing(6)
        combined_header_layout.addWidget(
            QLabel("Pipeline steps", combined_header), stretch=1
        )
        combined_header_layout.addWidget(
            make_help_button(combined_header, "Combined MSA clean", self._COMBINED_HELP)
        )
        combined_layout.addWidget(combined_header)
        hint = QLabel("Check steps to include. Drag to set order.", combined_box)
        hint.setWordWrap(True)
        combined_layout.addWidget(hint)
        combined_layout.addWidget(self._method_order)

        # --- Options: right-click + combined --------------------------------
        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        rc_box = QGroupBox("Right-click action", self)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(
            QLabel(
                "Right-clicking the MSA Clean button runs, with saved settings:",
                self,
            )
        )
        for label, value in (
            ("Combined MSA clean", MSA_CLEAN_RIGHT_CLICK_COMBINED),
            ("Remove gappy columns", MSA_CLEAN_RIGHT_CLICK_GAPPY),
            ("Clean rows with divergence", MSA_CLEAN_RIGHT_CLICK_DIV),
            ("Clean rows with gaps", MSA_CLEAN_RIGHT_CLICK_GAP),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
            rc_box_layout.addWidget(radio)
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[MSA_CLEAN_RIGHT_CLICK_COMBINED].setChecked(True)

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="MSA Clean",
                checked=bool(opts.show_process_window),
                close_when_finished=bool(opts.close_process_window_when_finished),
            )
        )

        options_tab = QWidget(self)
        options_tab.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)
        options_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        options_layout.addWidget(rc_box)
        options_layout.addWidget(combined_box)
        options_layout.addWidget(process_box)
        options_layout.addStretch(1)
        self._options_tab = options_tab

        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        self._tabs.tabBar().setExpanding(False)
        self._tabs.tabBar().setDrawBase(False)
        self._tabs.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._tabs.addTab(self._gappy_tab, "Remove gappy columns")
        self._tabs.addTab(self._div_tab, "Clean with divergence")
        self._tabs.addTab(self._gap_tab, "Clean with gaps")
        self._tabs.addTab(self._options_tab, "Options")

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore],
            on_accept=self._on_accept,
            include_run=True,
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(8)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._tabs, 0)
        layout.addWidget(self._button_box, 0)
        self._tabs.currentChanged.connect(self._on_tab_changed)
        self._sync_active_inputs()
        self._sync_field_column()
        self._on_tab_changed(self._tabs.currentIndex())

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._fit_method_list_height()
        self._sync_field_column()
        QTimer.singleShot(0, self._fit_to_current_tab)

    def _on_file_input_field_activated(self) -> None:
        sender = self.sender()
        if isinstance(sender, TeMultiInputFilesWidget):
            self._focused_file_input = sender

    def apply_path_from_panel(self, path: Path) -> bool:
        if self._tabs.currentWidget() is self._options_tab:
            return False
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        define = self._active_define_inputs()
        if define is None:
            return False
        try:
            candidate = Path(path).expanduser().resolve()
        except OSError:
            return False
        return bool(define.accept_path(candidate))

    def apply_folder_from_panel(self, folder: Path) -> bool:
        return self.apply_path_from_panel(folder)

    def _active_define_inputs(self) -> TeMultiInputFilesWidget | None:
        current = self._tabs.currentWidget()
        if current is self._gappy_tab:
            return self._gappy_define_inputs
        if current is self._div_tab:
            return self._div_define_inputs
        if current is self._gap_tab:
            return self._gap_define_inputs
        return None

    def _active_input_radios(self) -> tuple[QRadioButton, QRadioButton] | None:
        current = self._tabs.currentWidget()
        if current is self._gappy_tab:
            return self._gappy_input_selected, self._gappy_input_define
        if current is self._div_tab:
            return self._div_input_selected, self._div_input_define
        if current is self._gap_tab:
            return self._gap_input_selected, self._gap_input_define
        return None

    def input_source(self) -> str:
        radios = self._active_input_radios()
        if radios is not None and radios[1].isChecked():
            return ALIGN_INPUT_DEFINE
        return ALIGN_INPUT_SELECTED

    def run_files(self) -> list[Path]:
        define = self._active_define_inputs()
        if define is None:
            return []
        return define.paths()

    def selected_method(self) -> str:
        current = self._tabs.currentWidget()
        if current is self._gappy_tab:
            return "msa_clean_gappy"
        if current is self._div_tab:
            return "msa_clean_div"
        if current is self._gap_tab:
            return "msa_clean_gap"
        return "msa_clean"

    def _sync_active_inputs(self) -> None:
        for selected, define, widget, form in (
            (
                self._gappy_input_selected,
                self._gappy_input_define,
                self._gappy_define_inputs,
                self._gappy_form,
            ),
            (
                self._div_input_selected,
                self._div_input_define,
                self._div_define_inputs,
                self._div_form,
            ),
            (
                self._gap_input_selected,
                self._gap_input_define,
                self._gap_define_inputs,
                self._gap_form,
            ),
        ):
            use_define = define.isChecked()
            form.setRowVisible(widget, use_define)
            if use_define:
                widget.setMinimumHeight(0)
                widget.setMaximumHeight(16777215)
            else:
                widget.setFixedHeight(0)
            widget.setVisible(use_define)
            widget.updateGeometry()
        if self.isVisible() and self._tabs.currentWidget() is not self._options_tab:
            self._msa_compact_size = None
            QTimer.singleShot(0, self._fit_to_current_tab)

    def _on_define_paths_changed(self) -> None:
        if not self.isVisible() or self._tabs.currentWidget() is self._options_tab:
            return
        self._msa_compact_size = None
        QTimer.singleShot(0, self._fit_to_current_tab)

    def _sync_field_column(self) -> None:
        sync_form_field_column_spacers(
            self,
            [
                self._gappy_sel_spacer,
                self._gappy_def_spacer,
                self._div_sel_spacer,
                self._div_def_spacer,
                self._gap_sel_spacer,
                self._gap_def_spacer,
            ],
            [self._gappy_form, self._div_form, self._gap_form],
        )

    def _on_tab_changed(self, _tab_index: int) -> None:
        self._sync_active_inputs()
        self._fit_to_current_tab()

    def _pin_msa_tab_height(self, page_h: int | None) -> None:
        """Pin tab stack to current page height so unused tabs cannot leave a gap."""
        pin_tool_tabs_height(self._tabs, page_h)

    def _fit_to_current_tab(self) -> None:
        """Shrink/grow the dialog height to the active tab’s content (BLAST-style)."""
        page = self._tabs.currentWidget()
        if page is None:
            return
        # Unlock height pins before measuring (Define-input growth / tab switch).
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        self._pin_msa_tab_height(None)
        for widget in (
            self._gappy_define_inputs,
            self._div_define_inputs,
            self._gap_define_inputs,
        ):
            if widget.isVisible():
                widget.setMinimumHeight(0)
                widget.setMaximumHeight(16777215)

        # Measure the wrap’s inner form page when present (not the tallest tab).
        measure = page
        page_layout = page.layout()
        if page_layout is not None and page_layout.count() == 1:
            inner = page_layout.itemAt(0).widget()
            if inner is not None:
                measure = inner
        measure.ensurePolished()
        measure.adjustSize()
        measure_layout = measure.layout()
        if measure_layout is not None:
            measure_layout.invalidate()
            measure_layout.activate()
        if page_layout is not None:
            page_layout.invalidate()
            page_layout.activate()
        hint = measure.minimumSizeHint().expandedTo(measure.sizeHint())
        content_h = hint.height()
        if measure_layout is not None:
            content_h = max(
                content_h,
                measure_layout.totalMinimumSize().height(),
                measure_layout.totalSizeHint().height(),
            )
        if page_layout is not None:
            content_h = max(
                content_h,
                page_layout.totalMinimumSize().height(),
                page_layout.totalSizeHint().height(),
            )

        self._pin_msa_tab_height(content_h)
        tab_bar_h = max(24, self._tabs.tabBar().sizeHint().height())
        tabs_h = tab_bar_h + content_h + 8
        buttons_h = self._button_box.sizeHint().height()
        margins = self.layout().contentsMargins()
        spacing = self.layout().spacing()
        height = (
            tabs_h
            + buttons_h
            + margins.top()
            + margins.bottom()
            + spacing
            + 4
        )
        wide = self.input_source() == ALIGN_INPUT_DEFINE and page is not self._options_tab
        width = tool_params_dialog_width(define_inputs=wide)
        self.setMinimumWidth(width)
        self.setMaximumWidth(16777215)
        self.setMaximumHeight(height)
        target_w = max(width, self.width() if self.isVisible() else width)
        self.resize(target_w, height)
        # Keep width freely resizable; leave tab height pinned so content stays tight.
        self.setMaximumHeight(16777215)

    def _on_accept(self) -> None:
        if not self._checked_method_ids():
            QMessageBox.warning(
                self,
                "MSA clean",
                "Enable at least one step under Options → Combined MSA clean.",
            )
            self._tabs.setCurrentIndex(self._TAB_OPTIONS)
            return
        self.accept()

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _add_method_list_item(self, method_id: str, *, checked: bool) -> None:
        label_text = self._method_labels[method_id]
        item = QListWidgetItem(label_text)
        item.setData(Qt.ItemDataRole.UserRole, method_id)
        item.setFlags(
            item.flags()
            | Qt.ItemFlag.ItemIsUserCheckable
            | Qt.ItemFlag.ItemIsEnabled
            | Qt.ItemFlag.ItemIsSelectable
            | Qt.ItemFlag.ItemIsDragEnabled
        )
        item.setCheckState(
            Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        )
        self._method_order.addItem(item)

    def _fit_method_list_height(self) -> None:
        count = max(1, self._method_order.count())
        row_h = self._method_order.sizeHintForRow(0)
        if row_h <= 0:
            row_h = 22
        frame = 2 * self._method_order.frameWidth()
        padding = 16
        self._method_order.setFixedHeight(row_h * count + frame + padding)

    def _on_method_item_changed(self, item: QListWidgetItem) -> None:
        if item.checkState() == Qt.CheckState.Checked:
            return
        if self._checked_method_ids():
            return
        self._method_order.blockSignals(True)
        item.setCheckState(Qt.CheckState.Checked)
        self._method_order.blockSignals(False)

    def _checked_method_ids(self) -> list[str]:
        checked: list[str] = []
        for index in range(self._method_order.count()):
            item = self._method_order.item(index)
            if item.checkState() == Qt.CheckState.Checked:
                checked.append(str(item.data(Qt.ItemDataRole.UserRole)))
        return checked

    def _all_method_ids_in_order(self) -> list[str]:
        return [
            str(self._method_order.item(index).data(Qt.ItemDataRole.UserRole))
            for index in range(self._method_order.count())
        ]

    def _restore_defaults(self) -> None:
        defaults = MsaCleanOptions()
        self._div_threshold.setValue(defaults.divergence.threshold)
        self._div_window.setValue(defaults.divergence.window_size)
        self._row_gap_threshold.setValue(defaults.gap_row.gap_threshold)
        self._row_gap_window.setValue(defaults.gap_row.window_size)
        self._gap_threshold.setValue(defaults.gappy.gap_threshold)
        self._simi_gap.setValue(defaults.gappy.simi_check_gap_thr)
        self._similarity.setValue(defaults.gappy.similarity_thr)
        self._min_nt.setValue(defaults.gappy.min_nucleotide)
        enabled = {
            "msa_clean_div": defaults.use_divergence,
            "msa_clean_gap": defaults.use_gap_row,
            "msa_clean_gappy": defaults.use_gappy_columns,
        }
        self._method_order.blockSignals(True)
        self._method_order.clear()
        for method_id in defaults.method_order:
            self._add_method_list_item(
                method_id,
                checked=enabled.get(method_id, False),
            )
        self._method_order.blockSignals(False)
        self._fit_method_list_height()
        self._right_click_buttons[MSA_CLEAN_RIGHT_CLICK_COMBINED].setChecked(True)
        self._show_process_window.setChecked(defaults.show_process_window)
        self._close_process_window.setChecked(
            defaults.close_process_window_when_finished
        )
        self._close_process_window.setEnabled(self._show_process_window.isChecked())
        for selected in (
            self._gappy_input_selected,
            self._div_input_selected,
            self._gap_input_selected,
        ):
            selected.setChecked(True)
        self._sync_active_inputs()
        self._fit_to_current_tab()

    def options(self) -> MsaCleanOptions:
        method_order = tuple(self._all_method_ids_in_order())
        checked = set(self._checked_method_ids())
        right_click_action = MSA_CLEAN_RIGHT_CLICK_COMBINED
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click_action = value
                break
        return MsaCleanOptions(
            use_divergence="msa_clean_div" in checked,
            use_gap_row="msa_clean_gap" in checked,
            use_gappy_columns="msa_clean_gappy" in checked,
            method_order=method_order,
            gappy=GappyColumnOptions(
                gap_threshold=self._gap_threshold.value(),
                simi_check_gap_thr=self._simi_gap.value(),
                similarity_thr=self._similarity.value(),
                min_nucleotide=self._min_nt.value(),
            ),
            divergence=DivergenceRowOptions(
                threshold=self._div_threshold.value(),
                window_size=self._div_window.value(),
            ),
            gap_row=GapRowOptions(
                gap_threshold=self._row_gap_threshold.value(),
                window_size=self._row_gap_window.value(),
            ),
            right_click_action=right_click_action,
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=self._close_process_window.isChecked(),
        )



class ChromosomeMultiSelectWidget(QWidget):
    """Filterable multi-select list of chromosome names from the genome FASTA."""

    def __init__(
        self,
        chromosome_names: list[str],
        *,
        chromosome_sizes: dict[str, int] | None = None,
        preselected: list[str] | None = None,
        chromosome_sort: str = CHROMOSOME_SORT_NATURAL,
        defer_populate: bool = False,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self._all_chrom_names = list(chromosome_names)
        self._chrom_sizes = {
            str(k): int(v) for k, v in (chromosome_sizes or {}).items()
        }
        self._sort_mode = chromosome_sort
        self._preselected = None if preselected is None else set(preselected)
        self._check_boxes: list[QCheckBox] = []
        self._populated = False
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self._filter = QLineEdit(self)
        self._filter.setPlaceholderText("Type to filter (e.g. chr5)")
        self._filter.textChanged.connect(self._apply_filter)

        self._scroll = QScrollArea(self)
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QScrollArea.Shape.StyledPanel)
        self._scroll.setFixedHeight(150)
        self._scroll.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        self._checks_container = QWidget(self._scroll)
        self._checks_layout = QVBoxLayout(self._checks_container)
        self._checks_layout.setContentsMargins(6, 6, 6, 6)
        self._checks_layout.setSpacing(4)
        self._scroll.setWidget(self._checks_container)

        btn_row = QHBoxLayout()
        select_all_btn = QPushButton("Show all")
        select_all_btn.clicked.connect(self._check_all)
        clear_btn = QPushButton("Show none")
        clear_btn.clicked.connect(self._uncheck_all)
        if not self._all_chrom_names:
            select_all_btn.setEnabled(False)
            clear_btn.setEnabled(False)
        btn_row.addWidget(select_all_btn)
        btn_row.addWidget(clear_btn)
        min_size_lbl = QLabel("≥", self)
        min_size_lbl.setToolTip("Select chromosomes at least this large (skips small contigs)")
        self._min_size_mb = QDoubleSpinBox(self)
        self._min_size_mb.setRange(0.1, 10_000.0)
        self._min_size_mb.setDecimals(1)
        self._min_size_mb.setSingleStep(0.5)
        self._min_size_mb.setValue(2.0)
        self._min_size_mb.setFixedWidth(64)
        self._min_size_mb.setToolTip("Minimum chromosome size (Mb)")
        mb_unit_lbl = QLabel("Mb", self)
        self._min_size_confirm = QPushButton("Confirm", self)
        self._min_size_confirm.setToolTip("Select chromosomes ≥ the size above")
        self._min_size_confirm.clicked.connect(self._select_by_min_size_spinbox)
        if not self._all_chrom_names:
            self._min_size_mb.setEnabled(False)
            self._min_size_confirm.setEnabled(False)
        btn_row.addWidget(min_size_lbl)
        btn_row.addWidget(self._min_size_mb)
        btn_row.addWidget(mb_unit_lbl)
        btn_row.addWidget(self._min_size_confirm)
        btn_row.addStretch()

        layout.addWidget(self._filter)
        layout.addWidget(self._scroll)
        layout.addLayout(btn_row)
        self.setMaximumHeight(280)

        if defer_populate:
            self._checks_layout.addWidget(
                QLabel("Expand Advanced Settings to load chromosomes…", self._checks_container)
            )
            self._checks_layout.addStretch(1)
        else:
            self.ensure_populated()

    def ensure_populated(self) -> None:
        """Build the checkbox list once (expensive for large genomes)."""
        if self._populated:
            return
        self._populated = True
        # Pick initial selection before the first build — avoid a second rebuild.
        if self._preselected is not None:
            selected: set[str] | None = set(self._preselected)
        elif self._chrom_sizes:
            selected = self._names_at_least_mb(2.0) or None
        else:
            selected = None
        self._rebuild_check_list(selected)

    def _select_by_min_size_spinbox(self) -> None:
        self._select_min_size_mb(self._min_size_mb.value())

    def _names_at_least_mb(self, min_mb: float) -> set[str]:
        min_bp = int(float(min_mb) * 1_000_000)
        return {
            name
            for name in self._all_chrom_names
            if self._chrom_sizes.get(name, 0) >= min_bp
        }

    @staticmethod
    def _natural_sort_key(name: str) -> list:
        return [
            int(part) if part.isdigit() else part.lower()
            for part in re.split(r"(\d+)", str(name))
        ]

    def _sorted_names(self, mode: str) -> list[str]:
        names = list(self._all_chrom_names)
        if mode == CHROMOSOME_SORT_SIZE_DESC:
            return sorted(
                names,
                key=lambda c: (-self._chrom_sizes.get(c, 0), self._natural_sort_key(c)),
            )
        if mode == CHROMOSOME_SORT_SIZE_ASC:
            return sorted(
                names,
                key=lambda c: (self._chrom_sizes.get(c, 0), self._natural_sort_key(c)),
            )
        return sorted(names, key=self._natural_sort_key)

    def _display_label(self, name: str) -> str:
        size = self._chrom_sizes.get(name)
        if size is not None and size > 0:
            return f"{name} ({format_genome_length_bp(size)})"
        return name

    @staticmethod
    def _chrom_name_from_checkbox(cb: QCheckBox) -> str:
        raw = cb.property("chrom_name")
        if raw is not None:
            return str(raw)
        return cb.text()

    def _rebuild_check_list(self, selected: set[str] | None = None) -> None:
        if selected is None:
            selected = {
                self._chrom_name_from_checkbox(cb)
                for cb in self._check_boxes
                if cb.isChecked()
            }
            if not selected and self._preselected is not None:
                selected = set(self._preselected)

        self._checks_container.setUpdatesEnabled(False)
        self._scroll.setUpdatesEnabled(False)
        try:
            while self._checks_layout.count():
                item = self._checks_layout.takeAt(0)
                widget = item.widget()
                if widget is not None:
                    widget.setParent(None)
                    widget.deleteLater()
            self._check_boxes.clear()

            ordered = self._sorted_names(self._sort_mode)
            if not ordered:
                self._checks_layout.addWidget(
                    QLabel("No items available", self._checks_container)
                )
                self._checks_layout.addStretch(1)
                return

            # Batch-create checkboxes; defer show until the list is complete.
            for name in ordered:
                cb = QCheckBox(self._display_label(name), self._checks_container)
                cb.setProperty("chrom_name", name)
                cb.setChecked(bool(selected) and name in selected)
                self._check_boxes.append(cb)
                self._checks_layout.addWidget(cb)
            self._checks_layout.addStretch(1)
            self._apply_filter(self._filter.text())
        finally:
            self._checks_container.setUpdatesEnabled(True)
            self._scroll.setUpdatesEnabled(True)

    def set_chromosome_sort(self, mode: str) -> None:
        if mode not in (
            CHROMOSOME_SORT_NATURAL,
            CHROMOSOME_SORT_SIZE_DESC,
            CHROMOSOME_SORT_SIZE_ASC,
        ):
            mode = CHROMOSOME_SORT_NATURAL
        self._sort_mode = mode
        if self._populated:
            self._rebuild_check_list()

    def _apply_filter(self, text: str) -> None:
        if not self._populated:
            return
        needle = text.strip().lower()
        for cb in self._check_boxes:
            cb.setVisible(not needle or needle in cb.text().lower())

    def selected_names(self) -> list[str]:
        if not self._populated:
            if self._preselected is not None:
                return list(self._preselected)
            if self._chrom_sizes:
                return list(self._names_at_least_mb(2.0))
            return list(self._all_chrom_names)
        needle = self._filter.text().strip().lower()
        names: list[str] = []
        for cb in self._check_boxes:
            if not cb.isChecked():
                continue
            if needle and not cb.isVisible():
                continue
            names.append(self._chrom_name_from_checkbox(cb))
        return names

    def _check_all(self) -> None:
        self.ensure_populated()
        for cb in self._check_boxes:
            if cb.isVisible():
                cb.setChecked(True)

    def _uncheck_all(self) -> None:
        self.ensure_populated()
        for cb in self._check_boxes:
            if cb.isVisible():
                cb.setChecked(False)

    def _select_min_size_mb(self, min_mb: float) -> None:
        self.ensure_populated()
        selected = self._names_at_least_mb(min_mb)
        # Toggle checks in place — do not rebuild hundreds of widgets.
        for cb in self._check_boxes:
            name = self._chrom_name_from_checkbox(cb)
            cb.setChecked(name in selected)

    def name_count(self) -> int:
        return len(self._all_chrom_names)

    def set_chromosomes(
        self,
        chromosome_names: list[str],
        *,
        chromosome_sizes: dict[str, int] | None = None,
        preselected: list[str] | None = None,
    ) -> None:
        """Replace the genome chromosome list (e.g. after changing reference)."""
        self._all_chrom_names = list(chromosome_names)
        self._chrom_sizes = {
            str(k): int(v) for k, v in (chromosome_sizes or {}).items()
        }
        self._preselected = None if preselected is None else set(preselected)
        was_populated = self._populated
        self._populated = False
        self._check_boxes.clear()
        while self._checks_layout.count():
            item = self._checks_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()
        if was_populated:
            self.ensure_populated()
        else:
            self._checks_layout.addWidget(
                QLabel(
                    "Expand Advanced Settings to load chromosomes…",
                    self._checks_container,
                )
            )
            self._checks_layout.addStretch(1)


class _TeInputFileRow(QWidget):
    """Single file path row with optional remove button."""

    remove_requested = Signal()
    activated = Signal(object)
    path_changed = Signal()

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        initial: Path | None = None,
        removable: bool = True,
        placeholder: str = "RepeatMasker .out or BLAST hit table",
        browse_title: str = "Select RepeatMasker .out or BLAST hit table",
        file_filter: str = "TE hit tables (*.out *_b.txt *_blast.txt);;RepeatMasker (*.out);;All files (*)",
        start_dir: str = "",
        path_mode: str = "file",
        help_title: str = "",
        help_text: str = "",
        path_edit_width: int | None = None,
    ) -> None:
        super().__init__(parent)
        self.setSizePolicy(
            QSizePolicy.Policy.Preferred
            if path_edit_width is not None
            else QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        row_kwargs: dict = {
            "initial": initial,
            "placeholder": placeholder,
            "browse_title": browse_title,
            "file_filter": file_filter,
            "start_dir": start_dir,
            "mode": path_mode,
        }
        if path_edit_width is not None:
            row_kwargs["minimum_width"] = path_edit_width
        self._row = PathBrowseRow(self, **row_kwargs)
        self._path = self._row.line_edit()
        if path_edit_width is not None:
            self._path.setFixedWidth(path_edit_width)
            self._path.setSizePolicy(
                QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed
            )
            self._row.setSizePolicy(
                QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed
            )
        self._row.activated.connect(lambda _w: self.activated.emit(self))
        self._row.path_changed.connect(self.path_changed.emit)

        remove_btn = QPushButton("Remove", self)
        remove_btn.setVisible(removable)
        remove_btn.clicked.connect(self.remove_requested.emit)
        self._remove_btn = remove_btn

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        layout.addWidget(
            self._row,
            stretch=0 if path_edit_width is not None else 1,
        )
        layout.addWidget(remove_btn)
        # Keep path/Browse/Remove columns aligned: always reserve ? width.
        if help_text:
            layout.addWidget(make_help_button(self, help_title or "Help", help_text))
        else:
            help_spacer = QWidget(self)
            help_spacer.setFixedSize(HELP_BUTTON_SIZE, HELP_BUTTON_SIZE)
            layout.addWidget(help_spacer)
        if path_edit_width is not None:
            layout.addStretch(1)

    def set_removable(self, removable: bool) -> None:
        self._remove_btn.setVisible(removable)

    def set_path(self, path: Path | None) -> None:
        self._row.set_path(path)

    def path(self) -> Path | None:
        return self._row.path()


class TeMultiInputFilesWidget(QWidget):
    """One or more file paths with add/remove controls."""

    paths_changed = Signal()
    field_activated = Signal()

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        initial: list[Path] | None = None,
        start_empty: bool = False,
        placeholder: str = "RepeatMasker .out or BLAST hit table",
        browse_title: str = "Select RepeatMasker .out or BLAST hit table",
        file_filter: str = "TE hit tables (*.out *_b.txt *_blast.txt);;RepeatMasker (*.out);;All files (*)",
        add_button_label: str = "Add another input",
        add_button_tooltip: str = (
            "Add another file or folder path. Click a file or folder in the "
            "panel to fill a row, or use Browse."
        ),
        allow_multiple: bool = True,
        start_dir: str = "",
        path_mode: str = "file",
        help_title: str = "",
        help_text: str = "",
        path_edit_width: int | None = None,
    ) -> None:
        super().__init__(parent)
        self._path_edit_width = path_edit_width
        self.setSizePolicy(
            QSizePolicy.Policy.Preferred
            if path_edit_width is not None
            else QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Minimum,
        )
        self._placeholder = placeholder
        self._browse_title = browse_title
        self._file_filter = file_filter
        self._allow_multiple = allow_multiple
        self._start_dir = start_dir
        self._path_mode = path_mode
        self._help_title = help_title
        self._help_text = help_text

        self._rows_layout = QVBoxLayout()
        self._rows_layout.setContentsMargins(0, 0, 0, 0)
        self._rows_layout.setSpacing(2)
        self._rows: list[_TeInputFileRow] = []
        self._active_row: _TeInputFileRow | None = None

        self._add_btn = QPushButton(add_button_label, self)
        self._add_btn.setToolTip(add_button_tooltip)
        self._add_btn.clicked.connect(self._on_add_clicked)
        # Single-value fields (one TE library) hide the add control entirely.
        self._add_btn.setVisible(allow_multiple)

        footer = QHBoxLayout()
        footer.setContentsMargins(0, 0, 0, 0)
        footer.addWidget(self._add_btn)
        footer.addStretch(1)

        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)
        outer.addLayout(self._rows_layout)
        outer.addLayout(footer)

        if initial:
            for path in initial:
                self.add_row(path)
        if not self._rows and not start_empty:
            self.add_row()

    def _on_add_clicked(self) -> None:
        self.add_row(focus=True)
        self._add_btn.setEnabled(False)
        QTimer.singleShot(0, lambda: self._add_btn.setEnabled(True))

    def add_row(self, path: Path | None = None, *, focus: bool = False) -> None:
        # Put ? on the first row only (after Browse/Remove).
        first_row = not self._rows
        row = _TeInputFileRow(
            self,
            initial=path,
            removable=len(self._rows) >= 1,
            placeholder=self._placeholder,
            browse_title=self._browse_title,
            file_filter=self._file_filter,
            start_dir=self._start_dir,
            path_mode=self._path_mode,
            help_title=self._help_title if first_row else "",
            help_text=self._help_text if first_row else "",
            path_edit_width=self._path_edit_width,
        )
        row.remove_requested.connect(lambda r=row: self._remove_row(r))
        row.activated.connect(self._on_row_activated)
        row.path_changed.connect(self._emit_paths_changed)
        self._rows.append(row)
        self._rows_layout.addWidget(row)
        self._refresh_remove_buttons()
        # Unlock any height pin from a parent dialog so new rows are not crushed.
        self.setMinimumHeight(0)
        self.setMaximumHeight(16777215)
        self.updateGeometry()
        self._emit_paths_changed()
        if focus:
            self._on_row_activated(row)
            row._path.setFocus(Qt.FocusReason.OtherFocusReason)

    def _on_row_activated(self, row: _TeInputFileRow) -> None:
        self._active_row = row
        self.field_activated.emit()

    def _remove_row(self, row: _TeInputFileRow) -> None:
        if row not in self._rows or len(self._rows) <= 1:
            return
        if self._active_row is row:
            self._active_row = None
        self._rows.remove(row)
        self._rows_layout.removeWidget(row)
        row.deleteLater()
        self._refresh_remove_buttons()
        self._emit_paths_changed()

    def _emit_paths_changed(self) -> None:
        self.paths_changed.emit()

    def _refresh_remove_buttons(self) -> None:
        removable = self._allow_multiple and len(self._rows) > 1
        for row in self._rows:
            row.set_removable(removable)

    def set_path_at_focus(self, path: Path) -> bool:
        """Fill the path field that last had keyboard focus."""
        if self._active_row is None or self._active_row not in self._rows:
            return False
        self._active_row.set_path(path)
        return True

    def set_primary_path(self, path: Path) -> bool:
        return self.set_path_at_focus(path)

    def accept_path(self, path: Path) -> bool:
        """Fill the focused row, else the first empty row, else append / overwrite."""
        if self.set_path_at_focus(path):
            return True
        for row in self._rows:
            if row.path() is None:
                row.set_path(path)
                self._on_row_activated(row)
                return True
        if self._allow_multiple:
            self.add_row(path, focus=True)
            return True
        if self._rows:
            self._rows[0].set_path(path)
            self._on_row_activated(self._rows[0])
            return True
        self.add_row(path, focus=True)
        return True

    def clear_paths(self) -> None:
        while len(self._rows) > 1:
            self._remove_row(self._rows[-1])
        if self._rows:
            self._rows[0].set_path(None)
        else:
            self.add_row()

    def set_paths(self, paths: list[Path] | None) -> None:
        paths = list(paths or [])
        while len(self._rows) > 1:
            self._remove_row(self._rows[-1])
        if not paths:
            if self._rows:
                self._rows[0].set_path(None)
            else:
                self.add_row()
            return
        if not self._rows:
            self.add_row(paths[0])
        else:
            self._rows[0].set_path(paths[0])
        for path in paths[1:]:
            self.add_row(path)

    def paths(self) -> list[Path]:
        out: list[Path] = []
        for row in self._rows:
            picked = row.path()
            if picked is not None:
                out.append(picked.expanduser().resolve())
        return out

class TeHitTableInputWidget(QWidget):
    """RepeatMasker .out / BLAST hit table path entry with browse."""

    def __init__(self, parent: QWidget | None = None, *, initial: Path | None = None) -> None:
        super().__init__(parent)
        self._row = _TeInputFileRow(self, initial=initial, removable=False)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._row, stretch=1)
        self.setMinimumWidth(440)

    def set_path(self, path: Path) -> None:
        self._row.set_path(path)

    def path(self) -> Path | None:
        return self._row.path()


class GeneInputFileWidget(QWidget):
    """GFF/GTF/BED/VCF annotation path entry."""

    def __init__(self, parent: QWidget | None = None, *, initial: Path | None = None) -> None:
        super().__init__(parent)
        self._row = PathBrowseRow(
            self,
            initial=initial,
            placeholder="GFF/GTF/BED/VCF annotation",
            browse_title="Select annotation file",
            file_filter=(
                "Annotations (*.gff *.gff3 *.gtf *.bed *.vcf *.vcf.gz);;All files (*)"
            ),
            mode="file",
        )
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._row, stretch=1)

    def set_path(self, path: Path | None) -> None:
        self._row.set_path(path)

    def path(self) -> Path | None:
        return self._row.path()


def _form_section_separator() -> QWidget:
    """Thin grey horizontal rule between major dialog sections."""
    host = QWidget()
    host.setFixedHeight(7)
    layout = QVBoxLayout(host)
    layout.setContentsMargins(0, 2, 0, 2)
    line = QFrame()
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Plain)
    line.setFixedHeight(1)
    line.setStyleSheet("background-color: palette(mid); border: none; max-height: 1px;")
    layout.addWidget(line)
    return host


class AnnotationFeatureTypesWidget(QWidget):
    """Checklist of GFF/GTF feature types (or BED / VCF classes).

    Styled like ``ChromosomeMultiSelectWidget``: bordered scroll list of
    checkboxes plus Show all / Show none. Types appear only after annotation
    files are provided.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._check_boxes: list[QCheckBox] = []
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        self._scroll = QScrollArea(self)
        self._scroll.setWidgetResizable(True)
        self._scroll.setFrameShape(QScrollArea.Shape.StyledPanel)
        self._scroll.setFixedHeight(150)
        self._scroll.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Fixed,
        )
        self._scroll.setToolTip(
            "Feature types from column 3 of GFF/GTF files. BED files use a single "
            "'interval' type. VCF files use SNP, INDEL, MNP, SV, and other variant classes."
        )
        self._checks_container = QWidget(self._scroll)
        self._checks_layout = QVBoxLayout(self._checks_container)
        self._checks_layout.setContentsMargins(6, 6, 6, 6)
        self._checks_layout.setSpacing(4)
        self._scroll.setWidget(self._checks_container)

        btn_row = QHBoxLayout()
        btn_row.setContentsMargins(0, 0, 0, 0)
        btn_row.setSpacing(6)
        select_all_btn = QPushButton("Show all", self)
        select_all_btn.clicked.connect(self._check_all)
        clear_btn = QPushButton("Show none", self)
        clear_btn.clicked.connect(self._uncheck_all)
        btn_row.addWidget(select_all_btn)
        btn_row.addWidget(clear_btn)
        btn_row.addStretch(1)

        self._select_all_btn = select_all_btn
        self._clear_btn = clear_btn
        layout.addWidget(self._scroll)
        layout.addLayout(btn_row)
        self.setMaximumHeight(220)
        self._set_empty_state()

    def _clear_checks(self) -> None:
        while self._checks_layout.count():
            item = self._checks_layout.takeAt(0)
            widget = item.widget() if item is not None else None
            if widget is not None:
                widget.setParent(None)
                widget.deleteLater()
        self._check_boxes.clear()

    def _set_empty_state(self, message: str | None = None) -> None:
        self._clear_checks()
        self._checks_layout.addWidget(
            QLabel(
                message
                or "Add an annotation file to list feature types…",
                self._checks_container,
            )
        )
        self._checks_layout.addStretch(1)
        self._select_all_btn.setEnabled(False)
        self._clear_btn.setEnabled(False)

    def _set_types(self, types: list[str], *, checked: set[str] | None = None) -> None:
        prev = checked if checked is not None else {t.lower() for t in self.selected_types()}
        self._clear_checks()
        if not types:
            self._set_empty_state("No feature types found in the annotation file.")
            return
        for feature_type in types:
            cb = QCheckBox(feature_type, self._checks_container)
            cb.setChecked(feature_type.lower() in prev)
            self._check_boxes.append(cb)
            self._checks_layout.addWidget(cb)
        self._checks_layout.addStretch(1)
        self._select_all_btn.setEnabled(True)
        self._clear_btn.setEnabled(True)

    def _check_all(self) -> None:
        for cb in self._check_boxes:
            cb.setChecked(True)

    def _uncheck_all(self) -> None:
        for cb in self._check_boxes:
            cb.setChecked(False)

    def selected_types(self) -> list[str]:
        return [cb.text() for cb in self._check_boxes if cb.isChecked()]

    def set_checked_types(self, types: list[str]) -> None:
        want = {t.lower() for t in types if t.strip()}
        for cb in self._check_boxes:
            cb.setChecked(cb.text().lower() in want)

    @staticmethod
    def _default_checked_for_types(types: list[str]) -> set[str]:
        """Prefer saved defaults that exist in the file; otherwise check all."""
        available = {t.lower(): t for t in types}
        preferred = [
            available[t.lower()]
            for t in DEFAULT_ANNOTATION_FEATURE_TYPES
            if t.lower() in available
        ]
        if preferred:
            return {t.lower() for t in preferred}
        return {t.lower() for t in types}

    def refresh_from_paths(self, paths: list[Path]) -> None:
        from repeatviewer.backends.gff_features import list_feature_types_in_annotation

        valid = [raw.expanduser() for raw in paths if raw.expanduser().is_file()]
        if not valid:
            self._set_empty_state()
            return
        merged: set[str] = set()
        for path in valid:
            merged.update(list_feature_types_in_annotation(path))
        if not merged:
            self._set_empty_state("No feature types found in the annotation file.")
            return
        ordered = sorted(merged, key=str.lower)
        if merged == {"interval"}:
            self._set_types(["interval"], checked={"interval"})
            return
        merged_lower = {t.lower() for t in merged}
        if merged_lower.issubset({"snp", "indel", "mnp", "sv", "variant", "other"}):
            default = [
                t for t in ordered if t.lower() in {"snp", "indel", "mnp", "sv"}
            ] or ordered
            prev = (
                {t.lower() for t in self.selected_types()}
                or {t.lower() for t in default}
            )
            self._set_types(ordered, checked=prev)
            return
        prev = (
            {t.lower() for t in self.selected_types()}
            or self._default_checked_for_types(ordered)
        )
        self._set_types(ordered, checked=prev)


class PerFileAnnotationFeatureTypesWidget(QWidget):
    """Feature types per annotation file (one checklist when only one file is selected)."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._stack = QVBoxLayout(self)
        self._stack.setContentsMargins(0, 0, 0, 0)
        self._stack.setSpacing(6)
        self._single = AnnotationFeatureTypesWidget(self)
        self._per_file_host = QWidget(self)
        self._per_file_layout = QVBoxLayout(self._per_file_host)
        self._per_file_layout.setContentsMargins(0, 0, 0, 0)
        self._per_file_layout.setSpacing(6)
        self._widgets_by_key: dict[str, AnnotationFeatureTypesWidget] = {}
        self._saved_by_key: dict[str, list[str]] = {}
        self._global_fallback: list[str] = []
        self._stack.addWidget(self._single)
        self._stack.addWidget(self._per_file_host)
        self._per_file_host.hide()

    @staticmethod
    def _path_key(path: Path) -> str:
        return str(path.expanduser().resolve())

    def set_global_fallback_types(self, types: list[str] | None) -> None:
        self._global_fallback = [str(t).strip() for t in (types or []) if str(t).strip()]

    def set_saved_by_file(self, by_file: dict[str, list[str]] | None) -> None:
        self._saved_by_key = {
            self._path_key(Path(raw)): list(types)
            for raw, types in (by_file or {}).items()
            if types
        }

    def _clear_per_file_widgets(self) -> None:
        while self._per_file_layout.count():
            item = self._per_file_layout.takeAt(0)
            widget = item.widget() if item is not None else None
            if widget is not None:
                widget.deleteLater()
        self._widgets_by_key.clear()

    def refresh_from_paths(self, paths: list[Path]) -> None:
        valid = [p.expanduser() for p in paths if p.expanduser().is_file()]
        if len(valid) <= 1:
            self._single.show()
            self._per_file_host.hide()
            self._clear_per_file_widgets()
            if valid:
                self._single.refresh_from_paths(valid)
                key = self._path_key(valid[0])
                if key in self._saved_by_key:
                    self._single.set_checked_types(self._saved_by_key[key])
                elif self._global_fallback:
                    self._single.set_checked_types(self._global_fallback)
            else:
                self._single.refresh_from_paths([])
            return

        self._single.hide()
        self._per_file_host.show()
        self._clear_per_file_widgets()
        for path in valid:
            key = self._path_key(path)
            group = QGroupBox(path.name, self._per_file_host)
            group.setToolTip(str(path))
            group_layout = QVBoxLayout(group)
            widget = AnnotationFeatureTypesWidget(group)
            widget.refresh_from_paths([path])
            if key in self._saved_by_key:
                widget.set_checked_types(self._saved_by_key[key])
            elif self._saved_by_key:
                by_name = {
                    Path(raw).name: types for raw, types in self._saved_by_key.items()
                }
                if path.name in by_name:
                    widget.set_checked_types(by_name[path.name])
            elif self._global_fallback:
                widget.set_checked_types(self._global_fallback)
            group_layout.addWidget(widget)
            self._per_file_layout.addWidget(group)
            self._widgets_by_key[key] = widget

    def selected_types(self) -> list[str]:
        if self._single.isVisible():
            return self._single.selected_types()
        if len(self._widgets_by_key) == 1:
            return next(iter(self._widgets_by_key.values())).selected_types()
        return []

    def selected_types_by_file(self) -> dict[str, list[str]]:
        if self._single.isVisible():
            types = self._single.selected_types()
            return {}
        out: dict[str, list[str]] = {}
        for key, widget in self._widgets_by_key.items():
            picked = widget.selected_types()
            if picked:
                out[key] = picked
        return out


_GENOME_DENSITY_REFERENCE_HELP = (
    "Reference genome for chromosome coordinates. "
    "The TE masking files must be compatible with the reference genome. "
    "For example, the “query sequence” of a RepeatMasker <code>.out</code> "
    "file must have the same name as the chromosome name of the reference "
    "genome."
)
_GENOME_DENSITY_GENE_REFERENCE_HELP = (
    "Reference genome for chromosome coordinates. "
    "Annotation files must use the same chromosome names as this genome "
    "(GFF/GTF <code>seqid</code>, BED chrom, or VCF <code>CHROM</code>). "
    "Features whose chromosome is not in the reference are skipped."
)
_GENOME_DENSITY_METH_REFERENCE_HELP = (
    "Reference genome for chromosome coordinates. "
    "Bismark coverage files must use the same chromosome names as this genome "
    "(column 1). Rows whose chromosome is not in the reference are skipped."
)
_GENOME_DENSITY_TE_NAMING_HELP = (
    "TE classification is recognized by the traditional naming method "
    "<code>family#order/superfamily</code>. Example: in "
    "<code>rnd-4_family-1340#LTR/Copia</code>, "
    "<b>family</b> = <code>rnd-4_family-1340</code>, "
    "<b>order</b> = <code>LTR</code>, "
    "<b>superfamily</b> = <code>Copia</code>."
)
_GENOME_DENSITY_MASKING_HELP = (
    "RepeatMasker <code>.out</code> or a BLAST hit table "
    "(<code>*_b.txt</code>) for TE tracks.<br><br>"
    "The TE masking files must be compatible with the reference genome. "
    "For example, the “query sequence” of a RepeatMasker <code>.out</code> "
    "file must have the same name as the chromosome name of the reference "
    "genome.<br><br>"
    "BLAST tables are tab-separated with columns:<br>"
    "<code>qseqid sseqid pident length mismatch gapopen qstart qend "
    "sstart send sstrand evalue bitscore</code><br><br>"
    f"{_GENOME_DENSITY_TE_NAMING_HELP}<br><br>"
    "Add multiple files to show each in different tracks."
)
_GENOME_DENSITY_CONSENSUS_HELP = (
    "TE consensus FASTA files. Sequences from each file are BLASTed against "
    "the reference genome. Each file is plotted on its own track.<br><br>"
    f"{_GENOME_DENSITY_TE_NAMING_HELP}"
)
_GENOME_DENSITY_GENE_HELP = (
    "GFF, GTF, BED, or VCF annotation file. Each file is drawn on its own "
    "track below TE tracks.<br><br>"
    "<b>Chromosome names must match the reference genome</b> (seqid / chrom / "
    "VCF <code>CHROM</code>). Features whose chromosome is not in the reference "
    "are skipped and will not appear on the map.<br><br>"
    "<b>VCF</b>: each variant becomes an interval from <code>POS</code> through "
    "the REF allele length. Variants are classified as <b>SNP</b>, <b>INDEL</b>, "
    "<b>MNP</b>, <b>SV</b> (symbolic <code>&lt;…&gt;</code> alts), or other. Use "
    "<b>Display content</b> to choose which classes to plot (same density "
    "binning as GFF/BED features)."
)
_GENOME_DENSITY_METHYLATION_HELP = (
    "Bismark coverage (.cov / .cov.gz). Each file adds coverage and "
    "methylation tracks per chromosome.<br><br>"
    "<b>Chromosome names must match the reference genome</b> "
    "(column 1). Rows whose chromosome is not in the reference are skipped."
    "<br><br>"
    "Tab-separated columns: chrom, start, end, % methylation, count methylated, "
    "count unmethylated. Example:<br>"
    "<code>chr1&nbsp;&nbsp;1000&nbsp;&nbsp;1000&nbsp;&nbsp;75.0&nbsp;&nbsp;3&nbsp;&nbsp;1</code>"
)
_GENOME_DENSITY_METH_COVERAGE_MIN_HELP = (
    "Minimum read depth at a cytosine site for it to be used "
    "(methylated + unmethylated counts in the .cov file). Sites below this "
    "threshold are ignored when building methylation tracks.<br><br>"
    "Example: with min = 5, a site with 3 + 1 = 4 reads is dropped; "
    "3 + 2 = 5 reads is kept."
)
_GENOME_DENSITY_METH_DIFF_HELP = (
    "When exactly two Bismark coverage files are provided, add a per-chromosome "
    "difference track comparing methylation level "
    "(first file above the axis, second file below). "
    "File order follows the Inputs list top to bottom.<br><br>"
    "Enabled only when exactly two methylation files are provided."
)
_GENOME_DENSITY_Y_AXIS_HELP = (
    "<b>Total TE length (bp) per window</b>: sum of TE bases overlapping each "
    "window. Example: a 100 kb window covering 12 kb of TE sequence → "
    "Y = 12000.<br><br>"
    "<b>TE hits number per window</b>: count of TE hits overlapping each "
    "window. Example: the same window with 40 TE hits → Y = 40."
)
_GENOME_DENSITY_SPLIT_STRAND_HELP = (
    "Draw separate forward- and reverse-strand tracks per chromosome for "
    "RepeatMasker .out TE hits and for GFF/GTF/BED/VCF annotation features."
)
_GENOME_DENSITY_COLLAPSE_OVERLAPS_HELP = (
    "For <b>Total TE length (bp) per window</b>, count each base pair at most "
    "once when TEs overlap. When color is by order, superfamily, or family, "
    "overlaps are collapsed separately within each TE type; with % identity "
    "coloring, all TEs in a window are collapsed together. With separate "
    "forward/reverse tracks, collapse applies within each strand track only."
    "<br><br>"
    "Unavailable when Y axis is <b>TE hits number per window</b> "
    "(hit counts are not merged by overlap)."
)
_GENOME_DENSITY_TE_DIFF_HELP = (
    "When exactly two TE hit-table sources are selected (masking files and/or "
    "consensus BLAST tables), add a per-chromosome diff track "
    "(exclusive TE bp: above = first source only, below = second only). "
    "Track order follows input file order top to bottom. "
    "This step is optional and can take extra time on large genomes."
    "<br><br>"
    "Enabled only when exactly two TE sources are provided."
)
_GENOME_DENSITY_COLOR_BY_HELP = (
    "<b>% identity</b>: color by how similar each TE hit is to the consensus "
    "(continuous scale).<br><br>"
    "<b>TE order</b>: color by the broad class after <code>#</code> "
    "(before <code>/</code>). Example: in "
    "<code>rnd-4_family-1340#LTR/Copia</code>, order is <b>LTR</b>. "
    "Other orders include LINE, DNA, RC.<br><br>"
    "<b>TE superfamily</b>: color by the label after <code>/</code>. "
    "Example: in <code>rnd-4_family-1340#LTR/Copia</code>, superfamily is "
    "<b>Copia</b>. Other superfamilies include Gypsy, L1, Helitron.<br><br>"
    "<b>TE family</b>: color by the name before <code>#</code> "
    "(finest level; most colors). Example: in "
    "<code>rnd-4_family-1340#LTR/Copia</code>, family is "
    "<b>rnd-4_family-1340</b>."
)
_GENOME_DENSITY_WINDOW_HELP = (
    "Sliding window size for density bins. TE, Gene, and Methylation tracks "
    "use the same bin size; changing it on one tab updates the others."
)
_GENOME_DENSITY_CHROM_ORDER_HELP = (
    "Order of chromosomes in the plot. Shared by all track types."
)
_GENOME_DENSITY_CHROMOSOMES_HELP = (
    "Which chromosomes to include. Shared by TE, Gene, and Methylation — "
    "changing it here updates every track."
)
_GENOME_DENSITY_PLOT_ALL_HELP = (
    "When checked, <b>Run</b> from any tab (TE, Gene, or Methylation) builds "
    "all three track types in one plot, using each tab’s defined inputs. "
    "Tabs with no inputs are skipped. All included tracks share the same "
    "window size and chromosome selection, and appear as separate tracks on "
    "the map.<br><br>"
    "When unchecked, <b>Run</b> uses only the currently selected track tab."
)


class TeDensityPlotDialog(QDialog):
    """Genome density plot (TE / Gene / Methylation tracks)."""

    _DIALOG_WIDTH = TOOL_PARAMS_DIALOG_WIDTH_WIDE
    _ADVANCED_SCROLL_MAX_HEIGHT = 420

    def __init__(
        self,
        *,
        initial: TeGenomeMapOptions | None = None,
        chromosome_names: list[str] | None = None,
        chromosome_sizes: dict[str, int] | None = None,
        input_file: Path | None = None,
        consensus_input_file: Path | None = None,
        project: RepeatViewerProject | None = None,
        project_genome: Path | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Genome density plot")
        self._opts = initial or TeGenomeMapOptions()
        self._project = project
        self._project_genome = (
            project_genome.expanduser().resolve()
            if project_genome is not None and project_genome.is_file()
            else None
        )
        self.setMinimumWidth(self._DIALOG_WIDTH)
        self._te_density_compact_size: QSize | None = None

        genome_selected = (self._opts.genome_path or "").strip() or None

        def _make_genome_row(help_text: str) -> GenomeBlastDatabaseRow:
            return GenomeBlastDatabaseRow(
                self,
                project=project,
                project_genome=self._project_genome,
                selected_path=genome_selected,
                manage_title="Manage Database",
                help_title="Reference genome",
                help_text=help_text,
            )

        # Same reference on TE / Gene / Methylation — change one, update all.
        self._genome_row_te = _make_genome_row(_GENOME_DENSITY_REFERENCE_HELP)
        self._genome_row_ann = _make_genome_row(_GENOME_DENSITY_GENE_REFERENCE_HELP)
        self._genome_row_meth = _make_genome_row(_GENOME_DENSITY_METH_REFERENCE_HELP)
        self._genome_rows = (
            self._genome_row_te,
            self._genome_row_ann,
            self._genome_row_meth,
        )
        self._genome_row = self._genome_row_te
        self._syncing_genome_rows = False

        masking_initial: list[Path] | None = None
        if input_file is not None:
            masking_initial = [input_file]
        self._masking_input = TeMultiInputFilesWidget(
            self,
            initial=masking_initial,
            placeholder="RepeatMasker .out or BLAST hit table",
            add_button_label="Add another masking file",
            add_button_tooltip=(
                "Add another RepeatMasker .out or BLAST hit table. "
                "Each file is drawn on its own track (first row = top track)."
            ),
            help_title="TE masking file",
            help_text=_GENOME_DENSITY_MASKING_HELP,
        )
        self._consensus_input = TeMultiInputFilesWidget(
            self,
            initial=[consensus_input_file] if consensus_input_file is not None else None,
            placeholder="TE consensus FASTA (.fasta / .fa)",
            browse_title="Select TE consensus FASTA",
            file_filter="FASTA (*.fasta *.fa *.fna);;All files (*)",
            add_button_label="Add another TE consensus file",
            add_button_tooltip=(
                "Add another TE consensus FASTA. Sequences are BLASTed against the "
                "project genome; each file is plotted on its own track."
            ),
            help_title="TE consensus file",
            help_text=_GENOME_DENSITY_CONSENSUS_HELP,
        )
        self._gene_input = TeMultiInputFilesWidget(
            self,
            placeholder="GFF/GTF/BED/VCF annotation",
            browse_title="Select annotation file",
            file_filter="Annotations (*.gff *.gff3 *.gtf *.bed *.vcf *.vcf.gz);;All files (*)",
            add_button_label="Add another annotation file",
            add_button_tooltip=(
                "Add another GFF, GTF, BED, or VCF annotation file. "
                "Each file is drawn on its own track below TE tracks. "
                "Chromosome names must match the reference genome."
            ),
            help_title="Annotation files",
            help_text=_GENOME_DENSITY_GENE_HELP,
        )
        self._methylation_input = TeMultiInputFilesWidget(
            self,
            placeholder="Bismark coverage (.cov / .cov.gz)",
            browse_title="Select Bismark coverage file",
            file_filter=(
                "Bismark coverage (*.cov *.cov.gz *bismark*.cov*);;All files (*)"
            ),
            add_button_label="Add another methylation file",
            add_button_tooltip=(
                "Add another Bismark coverage file. Each file adds four tracks per "
                "chromosome: total (genome), covered (.cov lines), covered fraction, "
                "and methylation level (0–1). Chromosome names must match the "
                "reference genome."
            ),
            help_title="Methylation file",
            help_text=_GENOME_DENSITY_METHYLATION_HELP,
        )
        self._focused_file_input: TeMultiInputFilesWidget | None = None
        for widget in (
            self._masking_input,
            self._consensus_input,
            self._gene_input,
            self._methylation_input,
        ):
            widget.field_activated.connect(self._on_file_input_field_activated)
            widget.paths_changed.connect(self._on_file_rows_changed)

        def _make_bin_kb() -> QSpinBox:
            spin = QSpinBox(self)
            spin.setRange(5, 500)
            spin.setValue(max(5, self._opts.bin_size_bp // 1000))
            spin.setSuffix(" kb")
            return spin

        self._bin_kb = _make_bin_kb()
        self._bin_kb_gene = _make_bin_kb()
        self._bin_kb_meth = _make_bin_kb()
        BlastParamsDialog._wire_synced_spins(
            (self._bin_kb, self._bin_kb_gene, self._bin_kb_meth)
        )
        self._density = QComboBox(self)
        self._density.addItem("Total TE length (bp) per window", TE_GENOME_MAP_DENSITY_TE_LENGTH)
        self._density.addItem("TE hits number per window", TE_GENOME_MAP_DENSITY_HIT_COUNT)
        idx = self._density.findData(self._opts.density_mode)
        if idx >= 0:
            self._density.setCurrentIndex(idx)
        self._color_by = QComboBox(self)
        self._color_by.addItem("% identity", TE_GENOME_MAP_COLOR_IDENTITY)
        self._color_by.addItem("TE order", TE_GENOME_MAP_COLOR_ORDER)
        self._color_by.addItem("TE superfamily", TE_GENOME_MAP_COLOR_SUPERFAMILY)
        self._color_by.addItem("TE family", TE_GENOME_MAP_COLOR_FAMILY)
        idx_color = self._color_by.findData(self._opts.color_by)
        if idx_color >= 0:
            self._color_by.setCurrentIndex(idx_color)
        self._methylation_coverage_min = QSpinBox(self)
        self._methylation_coverage_min.setRange(1, 10_000)
        self._methylation_coverage_min.setValue(max(1, int(self._opts.methylation_coverage_min)))
        self._split_strand = QCheckBox("Separate forward / reverse strands", self)
        self._split_strand.setChecked(bool(self._opts.split_strand_tracks))
        self._collapse_overlaps = QCheckBox("Collapse overlapped TEs", self)
        self._collapse_overlaps.setChecked(bool(self._opts.collapse_overlaps))
        self._show_te_diff = QCheckBox("Show TE difference track (2 TE sources)", self)
        self._show_te_diff.setChecked(bool(self._opts.show_te_diff_track))
        self._show_meth_diff = QCheckBox(
            "Show methylation difference track (2 methylation files)", self
        )
        self._show_meth_diff.setChecked(bool(self._opts.show_meth_diff_track))

        self._annotation_density = QComboBox(self)
        self._annotation_density.addItem(
            "Total element length per window", TE_GENOME_MAP_DENSITY_TE_LENGTH
        )
        self._annotation_density.addItem(
            "Element hit number per window", TE_GENOME_MAP_DENSITY_HIT_COUNT
        )
        ann_idx = self._annotation_density.findData(self._opts.annotation_density_mode)
        if ann_idx >= 0:
            self._annotation_density.setCurrentIndex(ann_idx)
        self._annotation_features = PerFileAnnotationFeatureTypesWidget(self)
        self._annotation_features.set_global_fallback_types(self._opts.annotation_feature_types)
        if self._opts.annotation_feature_types_by_file:
            self._annotation_features.set_saved_by_file(
                self._opts.annotation_feature_types_by_file,
            )

        self._chromosome_sort = QComboBox(self)
        self._chromosome_sort.addItem("Natural (name order)", CHROMOSOME_SORT_NATURAL)
        self._chromosome_sort.addItem("Size (largest first)", CHROMOSOME_SORT_SIZE_DESC)
        self._chromosome_sort.addItem("Size (smallest first)", CHROMOSOME_SORT_SIZE_ASC)
        self._chromosome_sort.setToolTip(
            "Order of chromosomes in the list below (for selecting which to show)."
        )
        for field in (
            self._bin_kb,
            self._bin_kb_gene,
            self._bin_kb_meth,
            self._density,
            self._color_by,
            self._annotation_density,
            self._methylation_coverage_min,
            self._chromosome_sort,
        ):
            field.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        self._chrom_picker = ChromosomeMultiSelectWidget(
            chromosome_names or [],
            chromosome_sizes=chromosome_sizes,
            preselected=self._opts.chromosomes_only,
            chromosome_sort=self._opts.chromosome_sort,
            defer_populate=True,
            parent=self,
        )
        sort_idx = self._chromosome_sort.findData(self._opts.chromosome_sort)
        if sort_idx >= 0:
            self._chromosome_sort.blockSignals(True)
            self._chromosome_sort.setCurrentIndex(sort_idx)
            self._chromosome_sort.blockSignals(False)
        self._chromosome_sort.currentIndexChanged.connect(self._on_chromosome_sort_changed)
        self._chrom_picker.setToolTip(
            "Each chromosome is one row in the map. Many small scaffolds are slower "
            "than a few large ones — set a minimum size (Mb) and click Confirm."
        )
        for row in self._genome_rows:
            row.selection_changed.connect(self._on_genome_row_changed)

        def _section(title: str, parent: QWidget | None = None) -> QLabel:
            return QLabel(f"<b>{title}</b>", parent)

        def _make_form(parent: QWidget) -> QFormLayout:
            form = QFormLayout(parent)
            form.setContentsMargins(0, 0, 0, 0)
            form.setLabelAlignment(
                Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
            )
            # Expanding so path + Browse rows fill the dialog width.
            form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
            form.setHorizontalSpacing(12)
            form.setVerticalSpacing(4)
            return form

        # --- TE track page -------------------------------------------------
        te_page = QWidget(self)
        te_page.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        te_form = _make_form(te_page)
        te_form.addRow(_section("Inputs", te_page))
        te_form.addRow("Reference genome:", self._genome_row_te)
        te_form.addRow("TE masking file:", self._masking_input)
        te_form.addRow("TE consensus file:", self._consensus_input)
        te_form.addRow(_form_section_separator())
        te_form.addRow(_section("Settings", te_page))
        te_form.addRow(
            "Window size:",
            field_with_help(
                self._bin_kb, "Window size", _GENOME_DENSITY_WINDOW_HELP
            ),
        )
        te_form.addRow(
            "Y axis content:",
            field_with_help(
                self._density, "Y axis content", _GENOME_DENSITY_Y_AXIS_HELP
            ),
        )
        te_form.addRow(
            "Color by:",
            field_with_help(self._color_by, "Color by", _GENOME_DENSITY_COLOR_BY_HELP),
        )
        self._split_spacer = QLabel("")
        self._collapse_spacer = QLabel("")
        self._te_diff_spacer = QLabel("")
        te_form.addRow(
            self._split_spacer,
            field_with_help(
                self._split_strand,
                "Separate forward / reverse strands",
                _GENOME_DENSITY_SPLIT_STRAND_HELP,
            ),
        )
        te_form.addRow(
            self._collapse_spacer,
            field_with_help(
                self._collapse_overlaps,
                "Collapse overlapped TEs",
                _GENOME_DENSITY_COLLAPSE_OVERLAPS_HELP,
            ),
        )
        te_form.addRow(
            self._te_diff_spacer,
            field_with_help(
                self._show_te_diff,
                "Show TE difference track",
                _GENOME_DENSITY_TE_DIFF_HELP,
            ),
        )
        self._te_page = te_page
        self._te_form = te_form

        # --- Gene track page -----------------------------------------------
        ann_page = QWidget(self)
        ann_page.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        ann_form = _make_form(ann_page)
        ann_form.addRow(_section("Inputs", ann_page))
        ann_form.addRow("Reference genome:", self._genome_row_ann)
        ann_form.addRow("Annotation files:", self._gene_input)
        ann_form.addRow(_form_section_separator())
        ann_form.addRow(_section("Settings", ann_page))
        ann_form.addRow(
            "Window size:",
            field_with_help(
                self._bin_kb_gene, "Window size", _GENOME_DENSITY_WINDOW_HELP
            ),
        )
        ann_form.addRow("Display content:", self._annotation_features)
        ann_form.addRow("Y axis content:", self._annotation_density)
        self._ann_page = ann_page
        self._ann_form = ann_form

        # --- Methylation track page ----------------------------------------
        meth_page = QWidget(self)
        meth_page.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        meth_form = _make_form(meth_page)
        meth_form.addRow(_section("Inputs", meth_page))
        meth_form.addRow("Reference genome:", self._genome_row_meth)
        meth_form.addRow("Methylation file:", self._methylation_input)
        meth_form.addRow(_form_section_separator())
        meth_form.addRow(_section("Settings", meth_page))
        meth_form.addRow(
            "Window size:",
            field_with_help(
                self._bin_kb_meth, "Window size", _GENOME_DENSITY_WINDOW_HELP
            ),
        )
        self._meth_diff_spacer = QLabel("")
        meth_form.addRow(
            "Methylation coverage min:",
            field_with_help(
                self._methylation_coverage_min,
                "Methylation coverage min",
                _GENOME_DENSITY_METH_COVERAGE_MIN_HELP,
            ),
        )
        meth_form.addRow(
            self._meth_diff_spacer,
            field_with_help(
                self._show_meth_diff,
                "Show methylation difference track",
                _GENOME_DENSITY_METH_DIFF_HELP,
            ),
        )
        self._meth_page = meth_page
        self._meth_form = meth_form

        # Shared Advanced: chromosomes (applies to all tracks).
        advanced_body = QWidget(self)
        advanced_form = _make_form(advanced_body)
        advanced_form.setVerticalSpacing(6)
        advanced_form.addRow(
            "Chromosome order:",
            field_with_help(
                self._chromosome_sort,
                "Chromosome order",
                _GENOME_DENSITY_CHROM_ORDER_HELP,
            ),
        )
        chrom_field = QWidget(advanced_body)
        chrom_field_row = QHBoxLayout(chrom_field)
        chrom_field_row.setContentsMargins(0, 0, 0, 0)
        chrom_field_row.setSpacing(6)
        chrom_field_row.setAlignment(Qt.AlignmentFlag.AlignTop)
        chrom_field_row.addWidget(self._chrom_picker, stretch=1)
        chrom_field_row.addWidget(
            make_help_button(
                chrom_field,
                "Chromosomes to show",
                _GENOME_DENSITY_CHROMOSOMES_HELP,
            )
        )
        advanced_form.addRow("Chromosomes to show:", chrom_field)
        self._te_density_advanced_form = advanced_form

        advanced_scroll = QScrollArea(self)
        advanced_scroll.setWidgetResizable(True)
        advanced_scroll.setFrameShape(QFrame.Shape.NoFrame)
        advanced_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff
        )
        advanced_scroll.setWidget(advanced_body)
        advanced_scroll.setMaximumHeight(self._ADVANCED_SCROLL_MAX_HEIGHT)
        self._advanced_scroll = advanced_scroll

        advanced = QWidget(self)
        self._te_density_advanced_panel = advanced
        advanced.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        advanced_layout = QVBoxLayout(advanced)
        # Match tab-page wrap margins (16) so Advanced labels/fields line up
        # with Inputs / Settings inside the tabs.
        advanced_layout.setContentsMargins(16, 4, 16, 0)
        advanced_layout.setSpacing(6)
        advanced_layout.addWidget(_form_section_separator())
        advanced_layout.addWidget(_section("Advanced Settings", advanced))
        advanced_layout.addWidget(advanced_scroll)
        advanced.setVisible(False)
        advanced.setMaximumHeight(0)

        self._gene_input.paths_changed.connect(self._refresh_annotation_feature_types)
        self._refresh_annotation_feature_types()

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="Genome density plot",
                checked=bool(
                    getattr(self._opts, "show_process_window", DEFAULT_SHOW_PROCESS_WINDOW)
                ),
                close_when_finished=bool(
                    getattr(
                        self._opts,
                        "close_process_window_when_finished",
                        DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED,
                    )
                ),
            )
        )

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)
        options_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        self._plot_all_together = QCheckBox(
            "Plot all track types together", options_tab
        )
        self._plot_all_together.setChecked(
            bool(getattr(self._opts, "plot_all_tracks_together", False))
        )
        self._plot_all_together.setToolTip(_GENOME_DENSITY_PLOT_ALL_HELP)
        plot_box = QGroupBox("Tracks", options_tab)
        plot_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        plot_box_layout = QVBoxLayout(plot_box)
        plot_box_layout.addWidget(
            field_with_help(
                self._plot_all_together,
                "Plot all track types together",
                _GENOME_DENSITY_PLOT_ALL_HELP,
                stretch=True,
            )
        )
        options_layout.addWidget(plot_box)
        options_layout.addWidget(process_box)
        self._options_tab = options_tab

        def _wrap_page(page: QWidget) -> QWidget:
            wrap = QWidget(self)
            wrap.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
            wrap_layout = QVBoxLayout(wrap)
            wrap_layout.setContentsMargins(16, 10, 16, 10)
            wrap_layout.setSpacing(4)
            wrap_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
            wrap_layout.addWidget(page)
            return wrap

        self._te_wrap = _wrap_page(te_page)
        self._ann_wrap = _wrap_page(ann_page)
        self._meth_wrap = _wrap_page(meth_page)
        self._last_track_tab = TE_GENOME_MAP_DIALOG_TAB_TE

        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        self._tabs.tabBar().setDrawBase(False)
        self._tabs.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._tabs.addTab(self._te_wrap, "TE")
        self._tabs.addTab(self._ann_wrap, "Gene")
        self._tabs.addTab(self._meth_wrap, "Methylation")
        self._tabs.addTab(options_tab, "Options")
        self._options_tab_index = 3
        self._tabs.currentChanged.connect(self._on_tab_changed)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._advanced_settings_btn = QPushButton("Advanced Settings ▸", self)
        self._advanced_settings_btn.setCheckable(True)
        self._advanced_settings_btn.setChecked(False)
        self._advanced_settings_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._advanced_settings_btn.toggled.connect(self._on_advanced_settings_toggled)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore, self._advanced_settings_btn],
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(8)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._tabs, 0)
        layout.addWidget(advanced, 0)
        layout.addWidget(self._button_box, 0)

        open_tab = opening_dialog_tab(
            self._opts.last_dialog_tab,
            allowed=TE_GENOME_MAP_DIALOG_TABS,
            first_tab=TE_GENOME_MAP_DIALOG_TAB_TE,
        )
        tab_map = {
            TE_GENOME_MAP_DIALOG_TAB_TE: self._te_wrap,
            TE_GENOME_MAP_DIALOG_TAB_ANNOTATION: self._ann_wrap,
            TE_GENOME_MAP_DIALOG_TAB_METHYLATION: self._meth_wrap,
            DIALOG_TAB_OPTIONS: options_tab,
        }
        self._tabs.setCurrentWidget(tab_map.get(open_tab, self._te_wrap))
        if open_tab in (
            TE_GENOME_MAP_DIALOG_TAB_TE,
            TE_GENOME_MAP_DIALOG_TAB_ANNOTATION,
            TE_GENOME_MAP_DIALOG_TAB_METHYLATION,
        ):
            self._last_track_tab = open_tab
        self._masking_input.paths_changed.connect(self._sync_te_option_checkboxes)
        self._consensus_input.paths_changed.connect(self._sync_te_option_checkboxes)
        self._methylation_input.paths_changed.connect(self._sync_te_option_checkboxes)
        self._density.currentIndexChanged.connect(self._sync_te_option_checkboxes)
        self._sync_te_option_checkboxes()
        self._sync_field_column()
        self.resize(self._DIALOG_WIDTH, 480)
        QTimer.singleShot(0, lambda: self._resize_for_advanced(False))

    def _sync_te_option_checkboxes(self) -> None:
        """Grey TE / methylation options that do not apply to current inputs."""
        te_sources = len(self._masking_input.paths()) + len(
            self._consensus_input.paths()
        )
        self._show_te_diff.setEnabled(te_sources == 2)
        meth_files = len(self._methylation_input.paths())
        self._show_meth_diff.setEnabled(meth_files == 2)
        if meth_files != 2:
            self._show_meth_diff.setChecked(False)
        length_mode = (
            str(self._density.currentData()) == TE_GENOME_MAP_DENSITY_TE_LENGTH
        )
        self._collapse_overlaps.setEnabled(length_mode)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._te_density_compact_size = None
        self.setMinimumWidth(self._DIALOG_WIDTH)
        self._sync_field_column()
        QTimer.singleShot(0, lambda: self._resize_for_advanced(False))
        # Warm chromosome checkboxes after the dialog is visible so the first
        # Advanced Settings click does not stall on a large genome.
        QTimer.singleShot(50, self._chrom_picker.ensure_populated)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def current_dialog_tab(self) -> str:
        current = self._tabs.currentWidget()
        if current is self._options_tab:
            return DIALOG_TAB_OPTIONS
        if current is self._ann_wrap:
            return TE_GENOME_MAP_DIALOG_TAB_ANNOTATION
        if current is self._meth_wrap:
            return TE_GENOME_MAP_DIALOG_TAB_METHYLATION
        return TE_GENOME_MAP_DIALOG_TAB_TE

    def active_plot_track_tab(self) -> str:
        """Track tab used when not plotting all together."""
        tab = self.current_dialog_tab()
        if tab in (
            TE_GENOME_MAP_DIALOG_TAB_TE,
            TE_GENOME_MAP_DIALOG_TAB_ANNOTATION,
            TE_GENOME_MAP_DIALOG_TAB_METHYLATION,
        ):
            self._last_track_tab = tab
            return tab
        return self._last_track_tab

    def plot_all_tracks_together(self) -> bool:
        return self._plot_all_together.isChecked()

    def _on_track_page(self) -> bool:
        return self._tabs.currentWidget() is not self._options_tab

    def _active_track_wrap(self) -> QWidget:
        current = self._tabs.currentWidget()
        if current in (self._te_wrap, self._ann_wrap, self._meth_wrap):
            return current
        return self._te_wrap

    def _sync_field_column(self) -> None:
        sync_form_field_column_spacers(
            self,
            [
                self._split_spacer,
                self._collapse_spacer,
                self._te_diff_spacer,
                self._meth_diff_spacer,
            ],
            [self._te_form, self._ann_form, self._meth_form, self._te_density_advanced_form],
        )

    def _on_advanced_settings_toggled(self, checked: bool) -> None:
        self._advanced_settings_btn.setText(
            "Advanced Settings ▾" if checked else "Advanced Settings ▸"
        )
        # Expand the panel first so the click feels immediate; populate chromosomes
        # only if the idle warm-up has not finished yet.
        self._set_advanced_visible(checked)
        self._resize_for_advanced(checked)
        if checked:
            self._sync_field_column()
        if checked and not self._chrom_picker._populated:
            QTimer.singleShot(0, self._finish_advanced_chrom_populate)

    def _finish_advanced_chrom_populate(self) -> None:
        self._chrom_picker.ensure_populated()
        if self._advanced_settings_btn.isChecked():
            self._resize_for_advanced(True)

    def _fit_advanced_scroll_to_content(self) -> int:
        """Size Advanced scroll to its form body (chrom picker is already height-capped)."""
        body = self._advanced_scroll.widget()
        if body is None:
            need = 280
        else:
            body.updateGeometry()
            need = max(body.minimumSizeHint().height(), body.sizeHint().height())
            # ChromosomeMultiSelectWidget maxHeight is 280; keep this tight.
            need = min(self._ADVANCED_SCROLL_MAX_HEIGHT, max(200, need))
        self._advanced_scroll.setMaximumHeight(need)
        return need

    def _set_advanced_visible(self, visible: bool) -> None:
        panel = self._te_density_advanced_panel
        panel.setVisible(visible)
        if visible:
            self._fit_advanced_scroll_to_content()
            panel.setMaximumHeight(16777215)
        else:
            panel.setMaximumHeight(0)
        panel.updateGeometry()

    def _on_tab_changed(self, _index: int) -> None:
        tab = self.current_dialog_tab()
        if tab in (
            TE_GENOME_MAP_DIALOG_TAB_TE,
            TE_GENOME_MAP_DIALOG_TAB_ANNOTATION,
            TE_GENOME_MAP_DIALOG_TAB_METHYLATION,
        ):
            self._last_track_tab = tab
        on_track = self._on_track_page()
        self._advanced_settings_btn.setVisible(on_track)
        if not on_track:
            if self._advanced_settings_btn.isChecked():
                self._advanced_settings_btn.blockSignals(True)
                self._advanced_settings_btn.setChecked(False)
                self._advanced_settings_btn.setText("Advanced Settings ▸")
                self._advanced_settings_btn.blockSignals(False)
                self._set_advanced_visible(False)
            self._set_track_wraps_equal_height(False)
            self._te_density_advanced_panel.setVisible(False)
            self._fit_options_content()
        else:
            want = self._advanced_settings_btn.isChecked()
            self._set_advanced_visible(want)
            self._te_density_compact_size = None
            self._resize_for_advanced(want)

    def _track_wraps(self) -> tuple[QWidget, QWidget, QWidget]:
        return (self._te_wrap, self._ann_wrap, self._meth_wrap)

    def _set_track_wraps_equal_height(self, equal: bool) -> None:
        """While Advanced is open, all track pages share the tallest height."""
        wraps = self._track_wraps()
        if not equal:
            for wrap in wraps:
                wrap.setMinimumHeight(0)
            return
        for wrap in wraps:
            wrap.setMinimumHeight(0)
            root = wrap.layout()
            if root is not None:
                root.invalidate()
                root.activate()
            wrap.updateGeometry()
        target = max(
            max(wrap.minimumSizeHint().height(), wrap.sizeHint().height())
            for wrap in wraps
        )
        for wrap in wraps:
            wrap.setMinimumHeight(target)

    def _dialog_stack_height(self, *, expanded: bool) -> int:
        """Exact height of tabs + optional Advanced + footer (no estimate padding)."""
        self._tabs.updateGeometry()
        self._button_box.updateGeometry()
        tabs_h = max(self._tabs.minimumSizeHint().height(), self._tabs.sizeHint().height())
        buttons_h = max(
            self._button_box.minimumSizeHint().height(),
            self._button_box.sizeHint().height(),
        )
        adv_h = 0
        if expanded:
            panel = self._te_density_advanced_panel
            self._fit_advanced_scroll_to_content()
            panel.updateGeometry()
            adv_h = max(panel.minimumSizeHint().height(), panel.sizeHint().height())
        lay = self.layout()
        margins = lay.contentsMargins() if lay is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 34
        spacing = lay.spacing() if lay is not None else 8
        # layout always has tabs, advanced (0 when collapsed), buttons → 2 gaps
        return margin_h + tabs_h + adv_h + buttons_h + spacing * 2

    def _apply_dialog_height(self, height: int) -> None:
        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        width = self._DIALOG_WIDTH
        height = min(max(height, 220), max_h)
        self.setMinimumWidth(width)
        self.setMinimumHeight(0)
        # Cap height to the content stack so AlignTop does not leave a blank
        # under the buttons; callers unlock max before measuring when growing.
        self.setMaximumHeight(height)
        self.resize(width, height)

    def _on_file_rows_changed(self) -> None:
        """Grow/shrink the dialog when masking/consensus/gene/meth rows change."""
        if not self.isVisible() or not self._on_track_page():
            return
        # Let the new row finish laying out before measuring.
        QTimer.singleShot(0, self._reflow_after_file_rows_changed)

    def _reflow_after_file_rows_changed(self) -> None:
        if not self._on_track_page():
            return
        self._resize_for_advanced(self._advanced_settings_btn.isChecked())

    def _resize_for_advanced(self, expanded: bool) -> None:
        if not self._on_track_page():
            self._fit_options_content()
            return
        # Unlock height so sizeHints are not crushed by a previous maxHeight cap
        # (adding file rows otherwise overlaps Browse/Remove text).
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        wrap = self._active_track_wrap()
        root = wrap.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        wrap.updateGeometry()
        self._sync_field_column()

        if expanded:
            self._set_track_wraps_equal_height(True)
        else:
            self._set_track_wraps_equal_height(False)

        wrap_h = max(wrap.minimumSizeHint().height(), wrap.sizeHint().height())
        pin_tool_tabs_height(self._tabs, wrap_h)
        height = self._dialog_stack_height(expanded=expanded)
        self._apply_dialog_height(height)
        if not expanded:
            self._te_density_compact_size = QSize(self._DIALOG_WIDTH, self.height())

    def _fit_options_content(self) -> None:
        page = self._options_tab
        self.setMaximumHeight(16777215)
        pin_tool_tabs_height(self._tabs, None)
        page.ensurePolished()
        page.adjustSize()
        root = page.layout()
        if root is not None:
            root.activate()
        page_hint = page.minimumSizeHint().expandedTo(page.sizeHint())
        content_h = page_hint.height()
        if root is not None:
            content_h = max(
                content_h,
                root.totalMinimumSize().height(),
                root.totalSizeHint().height(),
            )
        pin_tool_tabs_height(self._tabs, content_h)
        tab_bar_h = self._tabs.tabBar().sizeHint().height() + 8
        buttons_h = self._button_box.sizeHint().height()
        lay = self.layout()
        margins = lay.contentsMargins() if lay is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 34
        spacing = lay.spacing() if lay is not None else 8
        height = margin_h + tab_bar_h + content_h + buttons_h + spacing * 2
        self._apply_dialog_height(height)

    def _restore_defaults(self) -> None:
        defaults = TeGenomeMapOptions()
        self._bin_kb.setValue(max(5, defaults.bin_size_bp // 1000))
        idx = self._density.findData(defaults.density_mode)
        if idx >= 0:
            self._density.setCurrentIndex(idx)
        idx_color = self._color_by.findData(defaults.color_by)
        if idx_color >= 0:
            self._color_by.setCurrentIndex(idx_color)
        self._split_strand.setChecked(defaults.split_strand_tracks)
        self._collapse_overlaps.setChecked(defaults.collapse_overlaps)
        self._show_te_diff.setChecked(defaults.show_te_diff_track)
        self._show_meth_diff.setChecked(defaults.show_meth_diff_track)
        self._sync_te_option_checkboxes()
        self._methylation_coverage_min.setValue(defaults.methylation_coverage_min)
        ann_idx = self._annotation_density.findData(defaults.annotation_density_mode)
        if ann_idx >= 0:
            self._annotation_density.setCurrentIndex(ann_idx)
        self._annotation_features.set_global_fallback_types(
            list(DEFAULT_ANNOTATION_FEATURE_TYPES)
        )
        self._annotation_features.set_saved_by_file({})
        self._refresh_annotation_feature_types()
        sort_idx = self._chromosome_sort.findData(defaults.chromosome_sort)
        if sort_idx >= 0:
            self._chromosome_sort.setCurrentIndex(sort_idx)
        self._syncing_genome_rows = True
        try:
            for row in self._genome_rows:
                row.reset_to_default()
        finally:
            self._syncing_genome_rows = False
        self._genome_row = self._genome_row_te
        self._chrom_picker.set_chromosome_sort(defaults.chromosome_sort)
        self._reload_chromosomes_for_reference()
        self._plot_all_together.setChecked(defaults.plot_all_tracks_together)
        self.setMinimumWidth(self._DIALOG_WIDTH)
        self._show_process_window.setChecked(defaults.show_process_window)
        self._close_process_window.setChecked(
            defaults.close_process_window_when_finished
        )
        if self._advanced_settings_btn.isChecked():
            self._advanced_settings_btn.setChecked(False)
        else:
            self._resize_for_advanced(False)

    def _on_file_input_field_activated(self) -> None:
        sender = self.sender()
        if isinstance(sender, TeMultiInputFilesWidget):
            self._focused_file_input = sender

    def apply_path_from_panel(self, path: Path) -> bool:
        if self._focused_file_input is None:
            return False
        filled = self._focused_file_input.set_path_at_focus(path)
        if filled and self._focused_file_input is self._gene_input:
            self._refresh_annotation_feature_types()
        return filled

    def set_te_input_file(self, path: Path) -> bool:
        return self._masking_input.set_path_at_focus(path)

    def set_te_consensus_input_file(self, path: Path) -> bool:
        return self._consensus_input.set_path_at_focus(path)

    def set_gene_input_file(self, path: Path) -> bool:
        filled = self._gene_input.set_path_at_focus(path)
        if filled:
            self._refresh_annotation_feature_types()
        return filled

    def _refresh_annotation_feature_types(self) -> None:
        self._annotation_features.refresh_from_paths(self.gene_input_files())

    def _on_chromosome_sort_changed(self) -> None:
        mode = str(self._chromosome_sort.currentData() or CHROMOSOME_SORT_NATURAL)
        self._chrom_picker.set_chromosome_sort(mode)

    def set_methylation_input_file(self, path: Path) -> bool:
        return self._methylation_input.set_path_at_focus(path)

    def te_masking_files(self) -> list[Path]:
        return self._masking_input.paths()

    def te_consensus_files(self) -> list[Path]:
        return self._consensus_input.paths()

    def te_input_files(self) -> list[Path]:
        """Masking hit-table paths (backward-compatible alias)."""
        return self.te_masking_files()

    def gene_input_files(self) -> list[Path]:
        return self._gene_input.paths()

    def methylation_input_files(self) -> list[Path]:
        return self._methylation_input.paths()

    def include_te_tracks(self) -> bool:
        if self.plot_all_tracks_together():
            return True
        return self.active_plot_track_tab() == TE_GENOME_MAP_DIALOG_TAB_TE

    def include_annotation_tracks(self) -> bool:
        if self.plot_all_tracks_together():
            return True
        return self.active_plot_track_tab() == TE_GENOME_MAP_DIALOG_TAB_ANNOTATION

    def include_methylation_tracks(self) -> bool:
        if self.plot_all_tracks_together():
            return True
        return self.active_plot_track_tab() == TE_GENOME_MAP_DIALOG_TAB_METHYLATION

    def gene_input_file(self) -> Path | None:
        paths = self.gene_input_files()
        return paths[0] if paths else None

    def reference_genome_path(self) -> Path | None:
        """FASTA used for chromosome coordinates (current genome or a nucleotide DB)."""
        raw = self._genome_row.selected_path()
        if not raw:
            return self._project_genome
        path = Path(raw).expanduser().resolve()
        return path if path.is_file() else self._project_genome

    def _on_genome_row_changed(self) -> None:
        """Keep TE / Gene / Methylation Reference genome combos in lockstep."""
        if self._syncing_genome_rows:
            return
        sender = self.sender()
        if not isinstance(sender, GenomeBlastDatabaseRow):
            return
        data = sender._combo.currentData()
        # refresh(None) means "keep current"; pass project genome / "" for reference.
        if data is None:
            sync_path: str | None = (
                str(self._project_genome) if self._project_genome is not None else ""
            )
        else:
            sync_path = str(data)
        self._syncing_genome_rows = True
        try:
            for row in self._genome_rows:
                if row is sender:
                    continue
                row.refresh(selected_path=sync_path)
            self._genome_row = sender
        finally:
            self._syncing_genome_rows = False
        self._on_reference_genome_changed()

    def _on_reference_genome_changed(self) -> None:
        if not self.isVisible():
            return
        self._reload_chromosomes_for_reference()
        QTimer.singleShot(
            0,
            lambda: self._resize_for_advanced(self._advanced_settings_btn.isChecked()),
        )

    def _reload_chromosomes_for_reference(self) -> None:
        genome = self.reference_genome_path()
        names = chromosome_names_from_genome(genome)
        sizes = genome_chrom_sizes_from_fasta(genome)
        self._chrom_picker.set_chromosomes(names, chromosome_sizes=sizes)

    def options(self) -> TeGenomeMapOptions:
        # Chromosome checkboxes are built lazily (open is fast); populate on Save/Run.
        self._chrom_picker.ensure_populated()
        mode = self._density.currentData()
        if mode not in (TE_GENOME_MAP_DENSITY_HIT_COUNT, TE_GENOME_MAP_DENSITY_TE_LENGTH):
            mode = TE_GENOME_MAP_DENSITY_HIT_COUNT
        ann_mode = self._annotation_density.currentData()
        if ann_mode not in (TE_GENOME_MAP_DENSITY_HIT_COUNT, TE_GENOME_MAP_DENSITY_TE_LENGTH):
            ann_mode = DEFAULT_ANNOTATION_DENSITY
        ann_types = self._annotation_features.selected_types()
        ann_by_file = self._annotation_features.selected_types_by_file()
        selected = self._chrom_picker.selected_names()
        if selected and len(selected) == self._chrom_picker.name_count():
            selected = []
        chrom_sort = str(self._chromosome_sort.currentData() or CHROMOSOME_SORT_NATURAL)
        if chrom_sort not in (
            CHROMOSOME_SORT_NATURAL,
            CHROMOSOME_SORT_SIZE_DESC,
            CHROMOSOME_SORT_SIZE_ASC,
        ):
            chrom_sort = CHROMOSOME_SORT_NATURAL
        genome = self.reference_genome_path()
        genome_path = ""
        if (
            genome is not None
            and self._project_genome is not None
            and genome != self._project_genome
        ):
            genome_path = str(genome)
        elif genome is not None and self._project_genome is None:
            genome_path = str(genome)
        return TeGenomeMapOptions(
            bin_size_bp=max(5_000, int(self._bin_kb.value()) * 1000),
            density_mode=str(mode),
            color_by=str(self._color_by.currentData()),
            methylation_coverage_min=max(1, int(self._methylation_coverage_min.value())),
            split_strand_tracks=self._split_strand.isChecked(),
            collapse_overlaps=self._collapse_overlaps.isChecked(),
            show_te_diff_track=self._show_te_diff.isChecked(),
            show_meth_diff_track=self._show_meth_diff.isChecked(),
            show_window_ref_line=self._opts.show_window_ref_line,
            annotation_density_mode=str(ann_mode),
            annotation_feature_types=ann_types or None,
            annotation_feature_types_by_file=ann_by_file or None,
            chromosomes_only=selected or None,
            superfamilies_only=None,
            chromosome_sort=chrom_sort,
            plot_all_tracks_together=self.plot_all_tracks_together(),
            genome_path=genome_path,
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=(
                self._close_process_window.isChecked()
            ),
            last_dialog_tab=persistable_dialog_tab(
                self.current_dialog_tab(),
                first_tab=TE_GENOME_MAP_DIALOG_TAB_TE,
            ),
        )


class TeDivergencePlotDialog(QDialog):
    """TE divergence histogram settings and input hit table."""

    def __init__(
        self,
        *,
        initial: TeDivergencePlotOptions | None = None,
        input_file: Path | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("TE divergence plot")
        self.setMinimumWidth(560)
        self._opts = initial or TeDivergencePlotOptions()

        self._input = TeMultiInputFilesWidget(
            self,
            initial=[input_file] if input_file is not None else None,
            allow_multiple=False,
        )
        self._focused_file_input: TeMultiInputFilesWidget | None = None
        self._input.field_activated.connect(self._on_file_input_field_activated)

        self._bin_width = QDoubleSpinBox(self)
        self._bin_width.setRange(0.001, 0.25)
        self._bin_width.setDecimals(3)
        self._bin_width.setSingleStep(0.001)
        self._bin_width.setValue(float(self._opts.bin_width))
        self._bin_width.setToolTip(
            "Width of each bin on the x-axis (divergence = 1 − identity, range 0–1)."
        )

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.addRow("TE input files:", self._input)
        form.addRow("Bin size:", self._bin_width)

        buttons = _run_cancel_button_box(self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 18, 20, 16)
        layout.setSpacing(12)
        layout.addLayout(form)
        layout.addWidget(buttons)
        self.setLayout(layout)

    def _on_file_input_field_activated(self) -> None:
        sender = self.sender()
        if isinstance(sender, TeMultiInputFilesWidget):
            self._focused_file_input = sender

    def apply_path_from_panel(self, path: Path) -> bool:
        if self._focused_file_input is None:
            return False
        return self._focused_file_input.set_path_at_focus(path)

    def set_te_input_file(self, path: Path) -> bool:
        return self._input.set_path_at_focus(path)

    def te_input_files(self) -> list[Path]:
        return self._input.paths()

    def options(self) -> TeDivergencePlotOptions:
        return TeDivergencePlotOptions(bin_width=float(self._bin_width.value()))


class TeAbundancePlotDialog(QDialog):
    """RepeatMasker family rank-abundance plot settings."""

    def __init__(
        self,
        *,
        input_files: list[Path] | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        from repeatviewer.backends.te_abundance_plot import TE_ABUNDANCE_BP, TE_ABUNDANCE_PERCENT

        _prepare_option_dialog(self)
        self.setWindowTitle("TE abundance plot")
        self.setMinimumWidth(560)

        self._input = TeMultiInputFilesWidget(
            self,
            initial=input_files,
            placeholder="RepeatMasker .out file",
            browse_title="Select RepeatMasker output",
            file_filter="RepeatMasker output (*.out);;All files (*)",
            add_button_label="Add another RepeatMasker file",
            add_button_tooltip="Add another RepeatMasker output as a separately coloured curve.",
        )
        self._metric = QComboBox(self)
        self._metric.addItem("Percentage occupied in genome", TE_ABUNDANCE_PERCENT)
        self._metric.addItem("Occupied nucleotides (bp)", TE_ABUNDANCE_BP)
        self._top_n = QSpinBox(self)
        self._top_n.setRange(1, 100_000)
        self._top_n.setValue(100)
        self._top_n.setToolTip("Show the most abundant TE families, ranked by occupied bp.")

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.addRow("RepeatMasker .out files:", self._input)
        form.addRow("Y axis:", self._metric)
        form.addRow("Top families:", self._top_n)
        buttons = _run_cancel_button_box(self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 18, 20, 16)
        layout.setSpacing(12)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def input_files(self) -> list[Path]:
        return self._input.paths()

    def y_axis_metric(self) -> str:
        return str(self._metric.currentData())

    def top_n(self) -> int:
        return int(self._top_n.value())

    def apply_path_from_panel(self, path: Path) -> bool:
        return self._input.set_path_at_focus(path)


from repeatviewer.widgets.residue_color_grid import ResidueColorGrid as _ResidueColorGrid


_MSA_PLOT_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): plot the aligned FASTA selected in "
    "the Files panel when you click Run.<br><br>"
    "<b>Define inputs</b>: choose a different aligned FASTA for this run only."
)
_MSA_PLOT_SEQ_TYPE_HELP = (
    "How residues are coloured. <b>Auto-detect</b> inspects the alignment; "
    "choose <b>Nucleotide</b> or <b>Protein</b> to force a palette."
)
_MSA_PLOT_FORMAT_HELP = "PNG for a bitmap image, or PDF for a vector figure."
_MSA_PLOT_DPI_HELP = (
    "Image resolution for PNG (higher = sharper and larger file). PDF ignores DPI."
)
_MSA_PLOT_HEIGHT_HELP = (
    "Figure height reserved for each sequence row (centimeters). "
    "Increase to make rows taller."
)
_MSA_PLOT_DISTANCE_HELP = (
    "White distance between consecutive sequence rows (centimeters). "
    "Set to 0 to remove spacing."
)


class MsaPlotDialog(QDialog):
    """MSA overview plot (Align-style Inputs / Settings / Advanced colors)."""

    _DIALOG_WIDTH = TOOL_PARAMS_DIALOG_WIDTH
    _COLORS_SCROLL_MAX_HEIGHT = 320
    _SETTINGS_LABEL_REF = "Sequence distance:"

    def __init__(
        self,
        *,
        initial: MsaPlotJobOptions | None = None,
        input_file: Path | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        from repeatviewer.backends.msa_plot import (
            AMINO_ACID_CATEGORIES,
            AMINO_ACID_COLOR_PRESETS,
            DEFAULT_AMINO_ACID_COLORS,
            DEFAULT_NUCLEOTIDE_COLORS,
            MSA_PLOT_FORMAT_PDF,
            MSA_PLOT_FORMAT_PNG,
            MSA_PLOT_SEQ_AUTO,
            MSA_PLOT_SEQ_NUCLEOTIDE,
            MSA_PLOT_SEQ_PROTEIN,
            NUCLEOTIDE_CATEGORIES,
            NUCLEOTIDE_COLOR_PRESETS,
        )

        _prepare_option_dialog(self)
        self.setWindowTitle("MSA plot")
        self.setMinimumWidth(self._DIALOG_WIDTH)
        self._opts = initial or MsaPlotJobOptions()
        self._default_nuc = dict(DEFAULT_NUCLEOTIDE_COLORS)
        self._default_aa = dict(DEFAULT_AMINO_ACID_COLORS)
        self._msa_compact_size: QSize | None = None
        self._msa_seq_auto = MSA_PLOT_SEQ_AUTO
        self._msa_seq_nuc = MSA_PLOT_SEQ_NUCLEOTIDE
        self._msa_seq_prot = MSA_PLOT_SEQ_PROTEIN
        self._msa_fmt_png = MSA_PLOT_FORMAT_PNG
        self._msa_fmt_pdf = MSA_PLOT_FORMAT_PDF
        self._msa_nuc_presets = list(NUCLEOTIDE_COLOR_PRESETS)
        self._msa_aa_presets = list(AMINO_ACID_COLOR_PRESETS)

        self._input_group = QButtonGroup(self)
        self._input_selected = QRadioButton("Use selected sequences", self)
        self._input_define = QRadioButton("Define inputs", self)
        self._input_group.addButton(self._input_selected)
        self._input_group.addButton(self._input_define)
        self._input_selected.setChecked(True)
        self._define_inputs = TeMultiInputFilesWidget(
            self,
            initial=None,
            allow_multiple=False,
            start_empty=False,
            placeholder="Aligned FASTA (.fa / .fasta / .faa)",
            browse_title="Select MSA FASTA",
            file_filter="FASTA (*.fa *.fasta *.fas *.fna *.faa);;All files (*)",
        )
        # Kept for call-site compatibility; Inputs uses panel selection by default.
        _ = input_file
        self._focused_file_input: TeMultiInputFilesWidget | None = None
        self._define_inputs.field_activated.connect(self._on_file_input_field_activated)
        self._input_selected.toggled.connect(self._sync_msa_inputs)
        self._input_define.toggled.connect(self._sync_msa_inputs)
        self._define_inputs.paths_changed.connect(self._sync_msa_inputs)

        self._seq_type = QComboBox(self)
        self._seq_type.addItem("Auto-detect", MSA_PLOT_SEQ_AUTO)
        self._seq_type.addItem("Nucleotide", MSA_PLOT_SEQ_NUCLEOTIDE)
        self._seq_type.addItem("Protein", MSA_PLOT_SEQ_PROTEIN)
        idx = self._seq_type.findData(self._opts.sequence_type)
        if idx >= 0:
            self._seq_type.setCurrentIndex(idx)

        self._format = QComboBox(self)
        self._format.addItem("PNG", MSA_PLOT_FORMAT_PNG)
        self._format.addItem("PDF", MSA_PLOT_FORMAT_PDF)
        idx_fmt = self._format.findData(self._opts.output_format)
        if idx_fmt >= 0:
            self._format.setCurrentIndex(idx_fmt)

        self._dpi = QSpinBox(self)
        self._dpi.setRange(50, 600)
        self._dpi.setSingleStep(25)
        self._dpi.setValue(int(self._opts.dpi))
        self._dpi.setSuffix(" DPI")

        self._seq_height = QDoubleSpinBox(self)
        self._seq_height.setRange(0.1, 5.0)
        self._seq_height.setDecimals(2)
        self._seq_height.setSingleStep(0.05)
        self._seq_height.setValue(float(self._opts.sequence_height_cm))
        self._seq_height.setSuffix(" cm")

        self._seq_distance = QDoubleSpinBox(self)
        self._seq_distance.setRange(0.0, 2.0)
        self._seq_distance.setDecimals(2)
        self._seq_distance.setSingleStep(0.01)
        self._seq_distance.setValue(float(self._opts.sequence_distance_cm))
        self._seq_distance.setSuffix(" cm")

        for field in (
            self._seq_type,
            self._format,
            self._dpi,
            self._seq_height,
            self._seq_distance,
        ):
            field.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)

        nuc_colors = dict(self._default_nuc)
        nuc_colors.update(self._opts.nucleotide_colors or {})
        if "U" in nuc_colors:
            nuc_colors["T"] = nuc_colors["U"]
        aa_colors = dict(self._default_aa)
        aa_colors.update(self._opts.amino_acid_colors or {})
        self._nuc_grid = _ResidueColorGrid(
            self,
            colors=nuc_colors,
            categories=list(NUCLEOTIDE_CATEGORIES),
        )
        self._aa_grid = _ResidueColorGrid(
            self,
            colors=aa_colors,
            categories=list(AMINO_ACID_CATEGORIES),
        )

        nuc_presets = QHBoxLayout()
        nuc_presets.setContentsMargins(0, 0, 0, 0)
        nuc_presets.setSpacing(6)
        nuc_preset_label = QLabel("Preset:", self)
        nuc_preset_label.setFixedWidth(_ResidueColorGrid.LABEL_WIDTH)
        nuc_preset_label.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        nuc_presets.addWidget(nuc_preset_label)
        nuc_tips = {
            "Default": "Standard A/C/G/T colours",
            "Color-blind": "Wong / Paul Tol palette (same as Annotation → nucleotide Color-blind)",
            "All black": "A, T, U, G, and C use a black background",
        }
        for label, palette in NUCLEOTIDE_COLOR_PRESETS:
            btn = QPushButton(label, self)
            btn.setToolTip(nuc_tips.get(label, f"Apply {label} nucleotide colours"))

            def _apply_nuc(colors: dict[str, str] = dict(palette)) -> None:
                synced = dict(colors)
                if "U" in synced:
                    synced["T"] = synced["U"]
                self._nuc_grid.set_colors(synced)

            btn.clicked.connect(lambda checked=False, apply=_apply_nuc: apply())
            nuc_presets.addWidget(btn)
        nuc_presets.addStretch(1)

        aa_presets = QHBoxLayout()
        aa_presets.setContentsMargins(0, 0, 0, 0)
        aa_presets.setSpacing(6)
        aa_preset_label = QLabel("Preset:", self)
        aa_preset_label.setFixedWidth(_ResidueColorGrid.LABEL_WIDTH)
        aa_preset_label.setAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        aa_presets.addWidget(aa_preset_label)
        aa_tips = {
            "Clustal": "Most common MSA colouring (ClustalX / Jalview)",
            "Zappo": "Clear physicochemical property groups",
            "Taylor": "High-contrast, almost one colour per residue",
        }
        for label, palette in AMINO_ACID_COLOR_PRESETS:
            btn = QPushButton(label, self)
            btn.setToolTip(aa_tips.get(label, f"Apply {label} amino-acid colours"))
            btn.clicked.connect(
                lambda checked=False, colors=dict(palette): self._aa_grid.set_colors(
                    colors
                )
            )
            aa_presets.addWidget(btn)
        aa_presets.addStretch(1)

        nuc_box = QGroupBox("Nucleotide colors", self)
        nuc_layout = QVBoxLayout(nuc_box)
        nuc_layout.addLayout(nuc_presets)
        nuc_layout.addWidget(self._nuc_grid)

        aa_box = QGroupBox("Amino-acid colors", self)
        aa_layout = QVBoxLayout(aa_box)
        aa_layout.addLayout(aa_presets)
        aa_layout.addWidget(self._aa_grid)

        root = QWidget(self)
        self._msa_root = root
        root.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        root_layout = QVBoxLayout(root)
        root_layout.setContentsMargins(16, 10, 16, 10)
        root_layout.setSpacing(4)

        basic = QWidget(root)
        self._msa_basic = basic
        basic.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        form = QFormLayout(basic)
        form.setContentsMargins(0, 0, 0, 0)
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        # Expanding so Define-inputs path + Browse fill the dialog width.
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)

        def _section(title: str) -> QLabel:
            return QLabel(f"<b>{title}</b>", root)

        form.addRow(_section("Inputs"))
        inputs_box = QWidget(basic)
        inputs_layout = QVBoxLayout(inputs_box)
        inputs_layout.setContentsMargins(0, 0, 0, 0)
        inputs_layout.setSpacing(2)
        inputs_layout.addWidget(
            field_with_help(
                self._input_selected,
                "MSA plot inputs",
                _MSA_PLOT_INPUTS_HELP,
                stretch=True,
            )
        )
        inputs_layout.addWidget(self._input_define)
        self._inputs_radios_label = QLabel("")
        self._inputs_define_label = QLabel("")
        form.addRow(self._inputs_radios_label, inputs_box)
        form.addRow(self._inputs_define_label, self._define_inputs)
        form.addRow(_form_section_separator())
        form.addRow(_section("Settings"))
        form.addRow(
            "Sequence type:",
            field_with_help(self._seq_type, "Sequence type", _MSA_PLOT_SEQ_TYPE_HELP),
        )
        form.addRow(
            "Output format:",
            field_with_help(self._format, "Output format", _MSA_PLOT_FORMAT_HELP),
        )
        form.addRow(
            "Quality:",
            field_with_help(self._dpi, "Quality", _MSA_PLOT_DPI_HELP),
        )
        form.addRow(
            "Sequence height:",
            field_with_help(self._seq_height, "Sequence height", _MSA_PLOT_HEIGHT_HELP),
        )
        form.addRow(
            "Sequence distance:",
            field_with_help(
                self._seq_distance, "Sequence distance", _MSA_PLOT_DISTANCE_HELP
            ),
        )
        self._msa_basic_form = form
        root_layout.addWidget(basic)

        advanced = QWidget(root)
        self._msa_advanced_panel = advanced
        advanced.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        advanced_form = QFormLayout(advanced)
        advanced_form.setContentsMargins(0, 4, 0, 0)
        advanced_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        advanced_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        advanced_form.setHorizontalSpacing(12)
        advanced_form.setVerticalSpacing(4)
        advanced_form.addRow(_form_section_separator())
        advanced_form.addRow(_section("Advanced Settings"))

        colors_host = QWidget(advanced)
        colors_layout = QVBoxLayout(colors_host)
        colors_layout.setContentsMargins(0, 0, 0, 0)
        colors_layout.setSpacing(8)
        colors_layout.addWidget(nuc_box)
        colors_layout.addWidget(aa_box)

        scroll = QScrollArea(advanced)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(colors_host)
        scroll.setMaximumHeight(self._COLORS_SCROLL_MAX_HEIGHT)
        self._colors_scroll = scroll
        advanced_form.addRow(scroll)
        root_layout.addWidget(advanced)
        root_layout.addStretch(1)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._advanced_settings_btn = QPushButton("Advanced Settings ▸", self)
        self._advanced_settings_btn.setCheckable(True)
        self._advanced_settings_btn.setChecked(False)
        self._advanced_settings_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._advanced_settings_btn.toggled.connect(self._on_advanced_settings_toggled)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore, self._advanced_settings_btn],
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(14)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(root, 0)
        layout.addWidget(self._button_box, 0)

        self._sync_msa_inputs(resize=False)
        self._set_msa_advanced_visible(False)
        self.resize(self._DIALOG_WIDTH, 420)
        QTimer.singleShot(0, lambda: self._resize_for_msa_advanced(False))

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._msa_compact_size = None
        self._sync_field_column()
        QTimer.singleShot(0, lambda: self._resize_for_msa_advanced(False))

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def input_source(self) -> str:
        if self._input_define.isChecked():
            return ALIGN_INPUT_DEFINE
        return ALIGN_INPUT_SELECTED

    def run_files(self, *, panel_files: list[Path] | None = None) -> list[Path]:
        from repeatviewer.backends.fastga_dotplot import expand_fastga_input_paths

        if self.input_source() == ALIGN_INPUT_DEFINE:
            return expand_fastga_input_paths(self._define_inputs.paths())
        return [Path(p) for p in (panel_files or [])]

    def _sync_msa_inputs(self, *, resize: bool = True) -> None:
        define = self.input_source() == ALIGN_INPUT_DEFINE
        form = self._msa_basic.layout()
        if isinstance(form, QFormLayout):
            form.setRowVisible(self._define_inputs, define)
        self._define_inputs.setVisible(define)
        if define:
            self._define_inputs.setMinimumHeight(0)
            self._define_inputs.setMaximumHeight(16777215)
        else:
            self._define_inputs.setFixedHeight(0)
        self._define_inputs.updateGeometry()
        self._msa_basic.updateGeometry()
        self._sync_field_column()
        if (
            resize
            and self.isVisible()
            and not self._advanced_settings_btn.isChecked()
        ):
            self._msa_compact_size = None
            self._resize_for_msa_advanced(False)

    def _sync_field_column(self) -> None:
        if not hasattr(self, "_inputs_radios_label"):
            return
        sync_form_field_column_spacers(
            self,
            [self._inputs_radios_label, self._inputs_define_label],
            [getattr(self, "_msa_basic_form", None)],
            ref_label=self._SETTINGS_LABEL_REF,
        )

    def _set_msa_advanced_visible(self, visible: bool) -> None:
        panel = self._msa_advanced_panel
        panel.setVisible(visible)
        if visible:
            panel.setMaximumHeight(16777215)
        else:
            panel.setMaximumHeight(0)
        panel.updateGeometry()
        root = self._msa_root.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        self._msa_root.updateGeometry()

    def _on_advanced_settings_toggled(self, checked: bool) -> None:
        self._advanced_settings_btn.setText(
            "Advanced Settings ▾" if checked else "Advanced Settings ▸"
        )
        if checked:
            self._msa_compact_size = self.size()
        self._set_msa_advanced_visible(checked)
        self._resize_for_msa_advanced(checked)

    def _dialog_chrome_height(self) -> int:
        lay = self.layout()
        margins = lay.contentsMargins() if lay is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 34
        spacing = lay.spacing() if lay is not None else 14
        buttons_h = self._button_box.sizeHint().height()
        return margin_h + spacing + buttons_h

    def _msa_compact_content_height(self) -> int:
        root = self._msa_root.layout()
        margins = root.contentsMargins() if root is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 20
        spacing = root.spacing() if root is not None else 4
        basic_h = max(
            self._msa_basic.minimumSizeHint().height(),
            self._msa_basic.sizeHint().height(),
        )
        return margin_h + basic_h + spacing

    def _msa_advanced_extra_height(self) -> int:
        """Header + scroll viewport — not the full unwrapped colour grids."""
        adv = self._msa_advanced_panel
        scroll_h = min(
            self._COLORS_SCROLL_MAX_HEIGHT,
            max(180, self._colors_scroll.sizeHint().height()),
        )
        # Separator + "Advanced Settings" title above the scroll.
        header = 36
        return header + scroll_h + 12

    def _resize_for_msa_advanced(self, expanded: bool) -> None:
        panel = self._msa_root
        root = panel.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        panel.updateGeometry()

        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        chrome = self._dialog_chrome_height()
        width = tool_params_dialog_width(
            define_inputs=self.input_source() == ALIGN_INPUT_DEFINE
        )
        self.setMinimumWidth(width)

        if not expanded:
            height = chrome + self._msa_compact_content_height()
            height = min(max(height, 220), max_h)
            self.setMinimumHeight(0)
            panel.setMinimumHeight(0)
            self.setMaximumHeight(height)
            self.resize(width, height)
            self.setMaximumHeight(16777215)
            self._msa_compact_size = QSize(width, height)
            return

        extra = self._msa_advanced_extra_height()
        base = self._msa_compact_size or self.size()
        height = min(base.height() + extra + 12, max_h)
        height = max(
            height,
            min(chrome + self._msa_compact_content_height() + extra, max_h),
        )
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def _restore_defaults(self) -> None:
        from repeatviewer.job_options import (
            DEFAULT_MSA_PLOT_DPI,
            DEFAULT_MSA_PLOT_SEQUENCE_DISTANCE_CM,
            DEFAULT_MSA_PLOT_SEQUENCE_HEIGHT_CM,
        )

        self._input_selected.setChecked(True)
        self._define_inputs.clear_paths()
        self._seq_type.setCurrentIndex(
            max(0, self._seq_type.findData(self._msa_seq_auto))
        )
        self._format.setCurrentIndex(
            max(0, self._format.findData(self._msa_fmt_png))
        )
        self._dpi.setValue(DEFAULT_MSA_PLOT_DPI)
        self._seq_height.setValue(DEFAULT_MSA_PLOT_SEQUENCE_HEIGHT_CM)
        self._seq_distance.setValue(DEFAULT_MSA_PLOT_SEQUENCE_DISTANCE_CM)
        self._nuc_grid.set_colors(dict(self._default_nuc))
        self._aa_grid.set_colors(dict(self._default_aa))
        self._sync_msa_inputs()

    def _on_file_input_field_activated(self) -> None:
        sender = self.sender()
        if isinstance(sender, TeMultiInputFilesWidget):
            self._focused_file_input = sender

    def apply_path_from_panel(self, path: Path) -> bool:
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        try:
            candidate = Path(path).expanduser().resolve()
        except OSError:
            return False
        if self._focused_file_input is not None:
            if self._focused_file_input.set_path_at_focus(candidate):
                return True
        return bool(self._define_inputs.accept_path(candidate))

    def input_files(self) -> list[Path]:
        """Legacy helper: defined paths only (prefer run_files + input_source)."""
        return self._define_inputs.paths()

    def options(self) -> MsaPlotJobOptions:
        nuc = self._nuc_grid.colors()
        if "U" in nuc:
            nuc["T"] = nuc["U"]
        return MsaPlotJobOptions(
            sequence_type=str(self._seq_type.currentData()),
            output_format=str(self._format.currentData()),
            dpi=int(self._dpi.value()),
            sequence_height_cm=float(self._seq_height.value()),
            sequence_distance_cm=float(self._seq_distance.value()),
            nucleotide_colors=nuc,
            amino_acid_colors=self._aa_grid.colors(),
        )


def _browse_file_line(
    parent: QWidget,
    line: QLineEdit,
    *,
    title: str,
    filt: str,
    placeholder: str = "",
) -> QWidget:
    return PathBrowseRow(
        parent,
        placeholder=placeholder or line.placeholderText(),
        browse_title=title,
        file_filter=filt,
        mode="file",
        line_edit=line,
    )


def _browse_dir_line(
    parent: QWidget,
    line: QLineEdit,
    *,
    title: str,
    placeholder: str = "",
) -> QWidget:
    return PathBrowseRow(
        parent,
        placeholder=placeholder or line.placeholderText(),
        browse_title=title,
        mode="directory",
        line_edit=line,
    )


class ConfusionMatrixBenchmarkDialog(QDialog):
    """Compare query RepeatMasker .out files to a reference library on the genome."""

    def __init__(
        self,
        *,
        initial: TeConfusionMatrixBenchmarkOptions | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Confusion matrix benchmarker")
        self.setMinimumWidth(620)
        opts = initial or TeConfusionMatrixBenchmarkOptions()

        self._tetrimmer_home = QLineEdit(opts.tetrimmer_home)
        self._tetrimmer_home.setPlaceholderText("TEtrimmer repository root folder")
        self._reference = QLineEdit(opts.reference_out)
        self._reference.setPlaceholderText("Reference RepeatMasker .out")
        self._less_class = QCheckBox("Simpler TE classes (merge LINE/SINE to Non-LTR)")
        self._less_class.setChecked(opts.less_classification)
        self._debug = QCheckBox("Keep intermediate files (debug)")
        self._debug.setChecked(opts.debug)

        form = QFormLayout()
        form.addRow(
            "TEtrimmer folder:",
            _browse_dir_line(
                self,
                self._tetrimmer_home,
                title="Select TEtrimmer repository root",
                placeholder="TEtrimmer repository root folder",
            ),
        )
        form.addRow(
            "Reference RM .out:",
            _browse_file_line(
                self,
                self._reference,
                title="Reference RepeatMasker .out",
                filt="RepeatMasker out (*.out);;All files (*)",
                placeholder="Reference RepeatMasker .out",
            ),
        )
        form.addRow("", self._less_class)
        form.addRow("", self._debug)

        buttons = _run_cancel_button_box(self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def options(self) -> TeConfusionMatrixBenchmarkOptions:
        return TeConfusionMatrixBenchmarkOptions(
            tetrimmer_home=self._tetrimmer_home.text().strip(),
            reference_out=self._reference.text().strip(),
            less_classification=self._less_class.isChecked(),
            debug=self._debug.isChecked(),
        )


class TeLibraryBenchmarkDialog(QDialog):
    """Benchmark a TE consensus library against a curated standard library."""

    def __init__(
        self,
        *,
        initial: TeLibraryBenchmarkOptions | None = None,
        query_file: Path | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("TE library benchmarker")
        self.setMinimumWidth(620)
        opts = initial or TeLibraryBenchmarkOptions()
        if query_file is not None:
            opts = TeLibraryBenchmarkOptions(
                query_fasta=str(query_file),
                standard_fasta=opts.standard_fasta,
                rm_out=opts.rm_out,
                run_repeatmasker=opts.run_repeatmasker,
                threads=opts.threads,
            )

        self._query = QLineEdit(opts.query_fasta)
        self._query.setPlaceholderText("Query TE library FASTA (.fa / .fasta)")
        self._standard = QLineEdit(opts.standard_fasta)
        self._standard.setPlaceholderText("Standard / curated library FASTA")
        self._rm_out = QLineEdit(opts.rm_out)
        self._rm_out.setPlaceholderText("Optional existing RepeatMasker .out")
        self._run_rm = QCheckBox("Run RepeatMasker (-lib standard query)")
        self._run_rm.setChecked(opts.run_repeatmasker)
        self._threads = QSpinBox()
        self._threads.setRange(1, 64)
        self._threads.setValue(max(1, int(opts.threads)))

        form = QFormLayout()
        form.addRow(
            "Query library FASTA:",
            _browse_file_line(
                self,
                self._query,
                title="Query TE library FASTA",
                filt="FASTA (*.fa *.fasta *.fas);;All files (*)",
                placeholder="Query TE library FASTA (.fa / .fasta)",
            ),
        )
        form.addRow(
            "Standard library FASTA:",
            _browse_file_line(
                self,
                self._standard,
                title="Standard TE library FASTA",
                filt="FASTA (*.fa *.fasta *.fas);;All files (*)",
                placeholder="Standard / curated library FASTA",
            ),
        )
        form.addRow(
            "Existing RM .out (optional):",
            _browse_file_line(
                self,
                self._rm_out,
                title="Existing lib-vs-lib RepeatMasker .out",
                filt="RepeatMasker out (*.out);;All files (*)",
                placeholder="Optional existing RepeatMasker .out",
            ),
        )
        form.addRow("", self._run_rm)
        form.addRow("RepeatMasker threads:", self._threads)

        buttons = _run_cancel_button_box(self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)
        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def options(self) -> TeLibraryBenchmarkOptions:
        return TeLibraryBenchmarkOptions(
            query_fasta=self._query.text().strip(),
            standard_fasta=self._standard.text().strip(),
            rm_out=self._rm_out.text().strip(),
            run_repeatmasker=self._run_rm.isChecked(),
            threads=int(self._threads.value()),
        )


# Backward-compatible alias (removed from menus; use TeDensityPlotDialog).
TeGenomeMapDefaultsDialog = TeDensityPlotDialog


class _ElidedPathLabel(QLabel):
    """One-line path label that shortens the middle of the path to fit.

    ``prefix`` and ``suffix`` always stay whole; only ``path`` is elided, and the
    full text is available on hover and via ``full_text()``.
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._prefix = ""
        self._path = ""
        self._suffix = ""
        self.setWordWrap(False)
        # Never widen the dialog for a long path; shrink the text instead.
        self.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Preferred)
        self.setMinimumWidth(1)

    def set_parts(self, prefix: str, path: str, suffix: str = "") -> None:
        self._prefix, self._path, self._suffix = prefix, path, suffix
        self.setToolTip(self.full_text())
        self._refit()

    def full_text(self) -> str:
        return f"{self._prefix}{self._path}{self._suffix}"

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._refit()

    def _refit(self) -> None:
        metrics = self.fontMetrics()
        fixed = metrics.horizontalAdvance(self._prefix + self._suffix)
        room = max(0, self.contentsRect().width() - fixed)
        middle = Qt.TextElideMode.ElideMiddle
        # Keep the last component (the run folder name) whole and shorten the
        # parent path in its middle; only if the name alone is too wide does
        # the whole path get shortened.
        cut = max(self._path.rfind("/"), self._path.rfind("\\")) + 1
        head, tail = self._path[:cut], self._path[cut:]
        tail_width = metrics.horizontalAdvance(tail)
        if metrics.horizontalAdvance(self._path) <= room:
            shown = self._path
        elif head and tail_width < room:
            shown = metrics.elidedText(head, middle, room - tail_width) + tail
        else:
            shown = metrics.elidedText(self._path, middle, room)
        super().setText(f"{self._prefix}{shown}{self._suffix}")


_RQ_LIBRARY_HELP = (
    "<b>Library</b>: the repeat library FASTA to quantify (monomers, "
    "satellites or TE consensus sequences).<br>"
    "Click a FASTA in the Files panel to fill this field, or use Browse."
)
_RQ_GENOME_HELP = (
    "<b>Genome</b>: the assembly FASTA (scaffolds or chromosomes) to search.<br>"
    "<b>Add another genome</b> quantifies the same library against more "
    "genomes; each genome gets its own report.<br>"
    "Click a FASTA in the Files panel to fill the field that last had focus."
)
_RQ_OUTPUT_HELP = (
    "<b>Auto</b> (empty): each genome's results go to a new folder next to "
    "that genome.<br>"
    "<b>Chosen folder</b>: the run folders are created inside it.<br>"
    "Every run gets its own timestamped folder, so nothing is overwritten.<br>"
    "Click this field, then a folder in the Folders panel, to fill it."
)
_RQ_THREADS_HELP = (
    "<b>Threads</b>: CPU threads for BLASTn (<code>-num_threads</code>) and "
    "RepeatMasker (<code>-pa</code>)."
)
_RQ_WORD_SIZE_HELP = (
    "<b>-word_size</b>: BLASTn seed length.<br>"
    "Lower is more sensitive but much slower on big, repetitive genomes "
    "(7 can hang; 11 is a good balance)."
)
_RQ_EVALUE_HELP = "<b>-evalue</b>: report BLASTn hits at or below this E-value."
_RQ_DIV_HELP = (
    "<b>-div</b>: RepeatMasker reports only copies less diverged than this "
    "percentage.<br>"
    "Higher captures older copies, but the result is noisier."
)
_RQ_SENSITIVITY_HELP = (
    "<b>Sensitivity</b>: RepeatMasker search speed.<br>"
    "<b>Sensitive (-s)</b>: slowest, finds the most copies.<br>"
    "<b>Default</b>: RepeatMasker's own default.<br>"
    "<b>Quick (-q)</b>: fastest on large genomes."
)
_RQ_PERC_IDENTITY_HELP = (
    "<b>-perc_identity</b>: minimum BLASTn percent identity per hit.<br>"
    "0 = no minimum."
)
_RQ_QCOV_HELP = (
    "<b>-qcov_hsp_perc</b>: minimum query (monomer) coverage per hit.<br>"
    "0 = no minimum."
)
_RQ_BLAST_INDEX_HELP = (
    "Repeat Quantify searches each genome's BLAST index. If the genome is "
    "already in Add/Manage databases (and the FASTA has not changed since), "
    "that index is reused and makeblastdb is skipped.<br>"
    "<b>Ticked</b>: a genome without a saved index is indexed once and added "
    "to Add/Manage databases, so the next run and the BLAST, TE-Aid and Extend "
    "tools reuse it.<br>"
    "<b>Unticked</b>: the index is temporary and rebuilt when the system "
    "clears its temp folder."
)


class RepeatQuantifyDialog(QDialog):
    """Repeat Quantify parameters (Tools → Repeat Quantify).

    Tabs: **BLASTn** | **RepeatMasker** | **Both** | **Options**, laid out like
    the BLAST and Mapping dialogs. The engine tabs share one page with
    Inputs / Outputs / Settings sections; each tab shows only its engine's
    rows. Restore Defaults and Advanced Settings sit bottom-left, Cancel /
    Save / Run bottom-right. Save stores the parameters without running.

    Both path fields start empty (nothing is guessed from the current folder,
    which could pick up a previous run's intermediate FASTAs). The dialog is
    shown non-modally, so the Folders and Files panels stay live: clicking a
    FASTA there drops it into whichever field last had focus
    (``apply_path_from_panel``), and Browse covers files elsewhere.

    The output folder is optional. Empty (Auto) writes each genome's results to
    a new ``<genome>_repeat_quant_<YYYY-MM-DD_HHMM>`` folder next to the genome;
    a chosen folder gets those new run folders inside it. It is deliberately not
    saved with the parameters, so every run starts on Auto.
    """

    _FASTA_FILTER = "FASTA (*.fa *.fasta *.fna *.fas);;All files (*)"
    # Mirrors _FASTA_FILTER: what a Files-panel click may drop into a field.
    _FASTA_SUFFIXES = frozenset({".fa", ".fasta", ".fna", ".fas"})
    _ENGINES = ("blastn", "repeatmasker", "both")
    _TAB_OPTIONS = 3
    _DIALOG_WIDTH = TOOL_PARAMS_DIALOG_WIDTH

    def __init__(
        self,
        *,
        start_dir: str = "",
        initial: RepeatQuantifyOptions | None = None,
        blast_databases: list | None = None,
        use_system_blast: bool = False,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Repeat Quantify parameters")
        self.setMinimumWidth(self._DIALOG_WIDTH)
        # Saved BLAST databases (Add/Manage databases), to show which genomes
        # already have an index Repeat Quantify will reuse.
        self._blast_databases = list(blast_databases or [])
        self._start_dir = start_dir
        opts = initial or RepeatQuantifyOptions()
        self._engine = opts.engine if opts.engine in self._ENGINES else "both"

        # --- Inputs ----------------------------------------------------------
        # Path fields rather than dropdowns, so the dialog can stay non-modal and
        # be filled by clicking files in the Files panel (see apply_path_from_panel).
        self._library_input = TeMultiInputFilesWidget(
            self,
            allow_multiple=False,
            start_dir=start_dir,
            placeholder="Library FASTA: click it in Files, or Browse",
            browse_title="Select repeat library FASTA",
            file_filter=self._FASTA_FILTER,
            help_title="Library",
            help_text=_RQ_LIBRARY_HELP,
        )
        self._genome_input = TeMultiInputFilesWidget(
            self,
            start_dir=start_dir,
            placeholder="Genome FASTA: click it in Files, or Browse",
            browse_title="Select genome FASTA",
            file_filter=self._FASTA_FILTER,
            add_button_label="Add another genome",
            add_button_tooltip=(
                "Quantify the same library against another genome. "
                "Each genome produces its own report."
            ),
            help_title="Genome",
            help_text=_RQ_GENOME_HELP,
        )
        self._focused_file_input: TeMultiInputFilesWidget | None = None
        for widget in (self._library_input, self._genome_input):
            widget.field_activated.connect(self._on_file_input_field_activated)
        self._db_status = QLabel(self)
        self._db_status.setWordWrap(True)
        self._db_status.setTextFormat(Qt.TextFormat.PlainText)
        self._db_status.setStyleSheet("color: palette(mid);")
        # A word-wrapped label asks for a narrow width; let it use the whole field.
        self._db_status.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred
        )

        # --- Settings / Advanced controls --------------------------------------
        self._threads = QSpinBox(self)
        self._threads.setRange(1, 256)
        self._word = QSpinBox(self)
        self._word.setRange(4, 24)
        self._evalue = QDoubleSpinBox(self)
        configure_report_evalue_spinbox(self._evalue)
        self._div = QSpinBox(self)
        self._div.setRange(1, 100)
        self._div.setSuffix(" %")
        self._sens = QComboBox(self)
        self._sens.addItem("Sensitive (-s)", "sensitive")
        self._sens.addItem("Default", "default")
        self._sens.addItem("Quick (-q)", "quick")
        self._minid = QDoubleSpinBox(self)
        self._minid.setRange(0.0, 100.0)
        self._minid.setDecimals(0)
        self._minid.setSuffix(" %")
        self._mincov = QDoubleSpinBox(self)
        self._mincov.setRange(0.0, 100.0)
        self._mincov.setDecimals(0)
        self._mincov.setSuffix(" %")
        self._save_index = QCheckBox("Save new indexes as BLAST databases", self)
        for field in (
            self._threads, self._word, self._evalue, self._div, self._sens,
            self._minid, self._mincov,
        ):
            field.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        for field in (
            self._threads, self._word, self._evalue, self._div,
            self._minid, self._mincov,
        ):
            field.setAlignment(Qt.AlignmentFlag.AlignRight)

        # --- engine page (one page; rows shown per engine tab) ----------------
        engine_page = QWidget(self)
        self._engine_page = engine_page
        engine_page.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )
        engine_root = QVBoxLayout(engine_page)
        engine_root.setContentsMargins(16, 10, 16, 10)
        engine_root.setSpacing(4)
        engine_root.setAlignment(Qt.AlignmentFlag.AlignTop)

        basic = QWidget(engine_page)
        self._engine_basic = basic
        basic.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        form = self._form(basic, grow=True)
        self._basic_form = form

        def _section(title: str, parent: QWidget) -> QLabel:
            return QLabel(f"<b>{title}</b>", parent)

        form.addRow(_section("Inputs", basic))
        form.addRow("Library:", self._library_input)
        form.addRow("Genome:", self._genome_input)
        self._db_status_label = QLabel("")
        form.addRow(self._db_status_label, self._db_status)

        form.addRow(_form_section_separator())
        form.addRow(_section("Outputs", basic))
        form.addRow("Output folder:", self._make_output_controls())

        form.addRow(_form_section_separator())
        form.addRow(_section("Settings", basic))
        form.addRow(
            "Threads:", field_with_help(self._threads, "Threads", _RQ_THREADS_HELP)
        )
        self._word_field = field_with_help(self._word, "-word_size", _RQ_WORD_SIZE_HELP)
        self._evalue_field = field_with_help(self._evalue, "-evalue", _RQ_EVALUE_HELP)
        self._div_field = field_with_help(self._div, "-div", _RQ_DIV_HELP)
        self._sens_field = field_with_help(
            self._sens, "Sensitivity", _RQ_SENSITIVITY_HELP
        )
        form.addRow("-word_size:", self._word_field)
        form.addRow("-evalue:", self._evalue_field)
        form.addRow("-div:", self._div_field)
        form.addRow("Sensitivity:", self._sens_field)
        engine_root.addWidget(basic)

        advanced = QWidget(engine_page)
        self._advanced_panel = advanced
        advanced.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        adv = self._form(advanced, grow=False)
        adv.setContentsMargins(0, 4, 0, 0)
        adv.addRow(_form_section_separator())
        adv.addRow(_section("Advanced Settings", advanced))
        self._minid_field = field_with_help(
            self._minid, "-perc_identity", _RQ_PERC_IDENTITY_HELP
        )
        self._mincov_field = field_with_help(
            self._mincov, "-qcov_hsp_perc", _RQ_QCOV_HELP
        )
        self._save_index_field = field_with_help(
            self._save_index, "BLAST index", _RQ_BLAST_INDEX_HELP
        )
        self._save_index_label = QLabel("")
        adv.addRow("-perc_identity:", self._minid_field)
        adv.addRow("-qcov_hsp_perc:", self._mincov_field)
        adv.addRow(self._save_index_label, self._save_index_field)
        self._advanced_form = adv
        engine_root.addWidget(advanced)

        # --- Options page ----------------------------------------------------
        options_page = QWidget(self)
        self._options_page = options_page
        options_page.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        ol = QVBoxLayout(options_page)
        ol.setContentsMargins(16, 14, 16, 14)
        ol.setSpacing(14)
        ol.setAlignment(Qt.AlignmentFlag.AlignTop)

        rc_box = QGroupBox("Right-click action", self)
        rc_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        rcl = QVBoxLayout(rc_box)
        rcl.setSpacing(6)
        rc_label = QLabel("Right-clicking the Repeat Quantify button:", self)
        rc_label.setWordWrap(True)
        rcl.addWidget(rc_label)
        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            (
                "Runs again on the last library and genomes, with saved settings",
                REPEAT_QUANTIFY_RIGHT_CLICK_RERUN,
            ),
        ):
            radio = QRadioButton(label, self)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
            rcl.addWidget(radio)
        ol.addWidget(rc_box)

        self._use_system_blast = QCheckBox(system_path_executable_label("BLAST+"), self)
        self._use_system_blast.setChecked(bool(use_system_blast))
        self._use_system_blast.setToolTip(
            "Use blastn and makeblastdb from PATH or the conda env in "
            "File → Settings → External tools instead of the bundled BLAST+. "
            "Shared with the BLAST window."
        )
        self._blast_path_status = make_executable_status_label(self)
        self._rm_path_status = make_executable_status_label(self)
        self._use_system_blast.toggled.connect(self._update_executable_status)
        exe_box = QGroupBox("Executables", self)
        exe_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        exl = QVBoxLayout(exe_box)
        exl.setContentsMargins(10, 8, 10, 8)
        exl.setSpacing(8)
        exl.addWidget(
            executable_checkbox_row(self, self._use_system_blast, self._blast_path_status)
        )
        exl.addWidget(self._rm_path_status)
        ol.addWidget(exe_box)

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="Repeat Quantify",
                checked=opts.show_process_window,
                close_when_finished=opts.close_process_window_when_finished,
            )
        )
        ol.addWidget(process_box)

        # --- tabs + stack ------------------------------------------------------
        self._tab_bar = QTabBar(self)
        self._tab_bar.setExpanding(False)
        self._tab_bar.setDocumentMode(True)
        for name in ("BLASTn", "RepeatMasker", "Both", "Options"):
            self._tab_bar.addTab(name)
        self._stack = QStackedWidget(self)
        self._stack.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        self._stack.addWidget(engine_page)
        self._stack.addWidget(options_page)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._advanced_settings_btn = QPushButton("Advanced Settings ▸", self)
        self._advanced_settings_btn.setCheckable(True)
        self._advanced_settings_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._advanced_settings_btn.toggled.connect(self._on_advanced_toggled)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore, self._advanced_settings_btn],
            on_accept=self._on_footer_accept,
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(8)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._tab_bar, 0)
        layout.addWidget(self._stack, 0)
        layout.addWidget(self._button_box, 0)

        self._apply_options(opts)
        self._genome_input.paths_changed.connect(self._update_db_status)
        self._genome_input.paths_changed.connect(self._refit_soon)
        self._save_index.toggled.connect(self._update_db_status)
        self._tab_bar.currentChanged.connect(self._on_tab_changed)
        self._set_advanced_visible(False)
        self._sync_field_column()
        self._on_tab_changed(self._tab_bar.currentIndex())
        self.resize(self._DIALOG_WIDTH, 420)

    # --- layout helpers ---------------------------------------------------------

    @staticmethod
    def _form(parent: QWidget, *, grow: bool) -> QFormLayout:
        form = QFormLayout(parent)
        form.setContentsMargins(0, 0, 0, 0)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
            if grow
            else QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)
        return form

    def _make_output_controls(self) -> QWidget:
        """Output folder field + Browse + ?, with a preview of the run folder below."""
        self._output_focused = False
        self._output_row = PathBrowseRow(
            self,
            placeholder="Auto: a new folder next to each genome",
            browse_title="Select output folder",
            mode="directory",
        )
        self._output_row.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self._output_path = self._output_row.line_edit()
        self._output_row.set_default_tooltip(
            "Leave empty to write each genome's results to a new folder next to "
            "that genome. Choose a folder to create the run folders inside it "
            "instead. Every run gets its own timestamped folder."
        )
        self._output_path.installEventFilter(self)
        self._output_row.text_changed.connect(lambda _t: self._update_output_preview())

        # One line, path shortened in the middle; full path on hover.
        self._output_preview = _ElidedPathLabel(self)
        self._output_preview.setStyleSheet("color: palette(mid);")

        host = QWidget(self)
        host.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        box = QVBoxLayout(host)
        box.setContentsMargins(0, 0, 0, 0)
        box.setSpacing(3)
        row = QHBoxLayout()
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        row.addWidget(self._output_row, stretch=1)
        row.addWidget(make_help_button(host, "Output folder", _RQ_OUTPUT_HELP))
        box.addLayout(row)
        box.addWidget(self._output_preview)

        self._genome_input.paths_changed.connect(self._update_output_preview)
        self._update_output_preview()
        return host

    def _sync_field_column(self) -> None:
        """Every label and spacer on the shared Settings field column."""
        sync_form_field_column_spacers(
            self,
            [self._db_status_label, self._save_index_label],
            [self._basic_form, self._advanced_form],
        )

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._sync_field_column()
        QTimer.singleShot(0, self._refit)

    def _refit_soon(self, *_args) -> None:
        if self.isVisible():
            QTimer.singleShot(0, self._refit)

    # --- tabs ---------------------------------------------------------------------

    def current_dialog_tab(self) -> str:
        """Active tab key (blastn / repeatmasker / both / options)."""
        index = self._tab_bar.currentIndex()
        if index == self._TAB_OPTIONS:
            return DIALOG_TAB_OPTIONS
        return self._ENGINES[index]

    def _on_tab_changed(self, index: int) -> None:
        on_options = index == self._TAB_OPTIONS
        self._stack.setCurrentIndex(1 if on_options else 0)
        if not on_options:
            self._engine = self._ENGINES[index]
        uses_blast = self._engine in ("blastn", "both")
        uses_rm = self._engine in ("repeatmasker", "both")
        form = self._basic_form
        form.setRowVisible(self._word_field, uses_blast)
        form.setRowVisible(self._evalue_field, uses_blast)
        form.setRowVisible(self._div_field, uses_rm)
        form.setRowVisible(self._sens_field, uses_rm)
        # Every Advanced row is a BLASTn option: no button on RepeatMasker / Options.
        show_advanced_btn = uses_blast and not on_options
        self._advanced_settings_btn.setVisible(show_advanced_btn)
        if not show_advanced_btn and self._advanced_settings_btn.isChecked():
            self._advanced_settings_btn.blockSignals(True)
            self._advanced_settings_btn.setChecked(False)
            self._advanced_settings_btn.setText("Advanced Settings ▸")
            self._advanced_settings_btn.blockSignals(False)
            self._set_advanced_visible(False)
        self._update_db_status()
        if on_options:
            self._update_executable_status()
        self._refit_soon()

    def _update_executable_status(self, *_args) -> None:
        """Refresh BLAST+ PATH / RepeatMasker lines on Options → Executables."""
        set_path_checkbox_status(
            self._blast_path_status,
            resolve_system_executable("blastn"),
            checked=self._use_system_blast.isChecked(),
            has_builtin=True,
        )
        set_named_executable_status(
            self._rm_path_status,
            "RepeatMasker",
            resolve_system_executable("RepeatMasker"),
            missing_text="not found on PATH",
        )
        self._refit_soon()

    # --- Advanced Settings and sizing (Mapping-dialog pattern) ------------------

    def _set_advanced_visible(self, visible: bool) -> None:
        panel = self._advanced_panel
        panel.setVisible(visible)
        panel.setMaximumHeight(16777215 if visible else 0)
        panel.updateGeometry()

    def _on_advanced_toggled(self, checked: bool) -> None:
        self._advanced_settings_btn.setText(
            "Advanced Settings ▾" if checked else "Advanced Settings ▸"
        )
        self._set_advanced_visible(checked)
        self._refit()

    def _dialog_chrome_height(self) -> int:
        lay = self.layout()
        margins = lay.contentsMargins()
        tab_h = max(24, self._tab_bar.sizeHint().height()) + 4
        return (
            margins.top() + margins.bottom() + lay.spacing() * 2
            + tab_h + self._button_box.sizeHint().height()
        )

    def _page_height(self, page: QWidget) -> int:
        page.ensurePolished()
        root = page.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        h = max(page.minimumSizeHint().height(), page.sizeHint().height())
        if root is not None:
            h = max(h, root.totalMinimumSize().height(), root.totalSizeHint().height())
        return h

    def _refit(self) -> None:
        """Size the window to the active page so no empty band is left."""
        page = self._stack.currentWidget()
        self._stack.setMaximumHeight(16777215)
        page.setMaximumHeight(16777215)
        content_h = max(1, self._page_height(page))
        self._stack.setMaximumHeight(content_h)
        page.setMaximumHeight(content_h)
        height = self._dialog_chrome_height() + content_h
        screen = self.screen().availableGeometry() if self.screen() else None
        if screen is not None:
            height = min(height, screen.height() - 40)
        self.setMinimumHeight(0)
        self.setMaximumHeight(height)
        self.resize(max(self.width(), self._DIALOG_WIDTH), height)
        self.setMaximumHeight(16777215)

    # --- output folder ------------------------------------------------------------

    def _update_output_preview(self, *_args) -> None:
        from repeatviewer.backends.repeat_quantify import (
            repeat_quant_output_dir,
            run_stamp,
        )

        stamp = run_stamp()
        base = self.output_dir()
        targets = self.targets()
        suffix = ""
        if targets:
            path = str(repeat_quant_output_dir(targets[0], stamp, base))
            if len(targets) > 1:
                suffix = f"  (+{len(targets) - 1} more)"  # one folder per genome
        else:
            name = f"<genome>_repeat_quant_{stamp}"
            path = str(base / name) if base else f"<genome folder>/{name}"
        self._output_preview.set_parts("→ ", path, suffix)
        if not self._output_path.text().strip() and targets:
            self._output_row.set_start_dir(str(targets[0].parent))

    def set_output_dir(self, folder: Path | None) -> None:
        """Set the output folder (None = Auto)."""
        if folder is None:
            self._output_path.clear()
        else:
            self._output_path.setText(str(Path(folder).expanduser().resolve()))

    def output_dir(self) -> Path | None:
        """Chosen output folder, or None for Auto (next to each genome)."""
        text = self._output_path.text().strip()
        return Path(text).expanduser() if text else None

    @staticmethod
    def run_stamp() -> str:
        """Timestamp for this run's folder names (taken at call time)."""
        from repeatviewer.backends.repeat_quantify import run_stamp

        return run_stamp()

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if obj is self._output_path and event.type() == QEvent.Type.FocusIn:
            self._output_focused = True
            self._focused_file_input = None
        return super().eventFilter(obj, event)

    # --- options ------------------------------------------------------------------

    def _apply_options(self, opts: RepeatQuantifyOptions) -> None:
        self._threads.setValue(max(1, int(opts.threads)))
        self._word.setValue(int(opts.blast_word_size))
        self._evalue.setValue(float(opts.blast_evalue))
        self._minid.setValue(float(opts.blast_min_identity))
        self._mincov.setValue(float(opts.blast_min_coverage))
        self._save_index.setChecked(bool(opts.blast_save_index))
        self._div.setValue(int(opts.rm_div))
        idx = self._sens.findData(opts.rm_sensitivity)
        self._sens.setCurrentIndex(max(0, idx))
        action = (
            opts.right_click_action
            if opts.right_click_action in self._right_click_buttons
            else REPEAT_QUANTIFY_RIGHT_CLICK_RERUN
        )
        self._right_click_buttons[action].setChecked(True)
        self._show_process_window.setChecked(bool(opts.show_process_window))
        self._close_process_window.setChecked(
            bool(opts.close_process_window_when_finished)
        )
        engine = opts.engine if opts.engine in self._ENGINES else "both"
        self._engine = engine
        self._tab_bar.setCurrentIndex(self._ENGINES.index(engine))

    def _restore_defaults(self) -> None:
        """Reset the parameters (not the chosen files) to their defaults."""
        self._apply_options(RepeatQuantifyOptions())

    def options(self) -> RepeatQuantifyOptions:
        right_click = next(
            (v for v, b in self._right_click_buttons.items() if b.isChecked()),
            REPEAT_QUANTIFY_RIGHT_CLICK_RERUN,
        )
        return RepeatQuantifyOptions(
            engine=self.engine(),
            blast_word_size=int(self._word.value()),
            blast_evalue=float(self._evalue.value()),
            blast_min_identity=float(self._minid.value()),
            blast_min_coverage=float(self._mincov.value()),
            rm_div=int(self._div.value()),
            rm_sensitivity=str(self._sens.currentData()),
            blast_save_index=self.save_blast_index(),
            threads=self.threads(),
            right_click_action=right_click,
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=self._close_process_window.isChecked(),
        )

    def threads(self) -> int:
        return int(self._threads.value())

    def use_system_blast(self) -> bool:
        return bool(self._use_system_blast.isChecked())

    def save_blast_index(self) -> bool:
        """Save a genome's new BLAST index to Add/Manage databases."""
        return bool(self._save_index.isChecked())

    def blast_database_status(self) -> list[str]:
        """One line per genome: whether its BLAST index is reused or built."""
        from repeatviewer.backends.database_prepare import saved_blast_database_for

        lines = []
        for target in self.targets():
            saved = saved_blast_database_for(target, self._blast_databases)
            if saved is not None:
                lines.append(f"{target.name}: uses saved BLAST database “{saved.name}”")
            elif self.save_blast_index():
                lines.append(
                    f"{target.name}: no saved BLAST database yet, so it is indexed "
                    "once and added to Add/Manage databases"
                )
            else:
                lines.append(
                    f"{target.name}: no saved BLAST database, so a temporary "
                    "index is built for this run"
                )
        return lines

    def _update_db_status(self, *_args) -> None:
        uses_blast = self.engine() in ("blastn", "both")
        lines = self.blast_database_status() if uses_blast else []
        self._db_status.setText("\n".join(lines))
        self._basic_form.setRowVisible(self._db_status, bool(lines))
        self._refit_soon()

    # --- run --------------------------------------------------------------------

    def should_run(self) -> bool:
        """True when the dialog was closed with Run (not Save or Cancel)."""
        return dialog_should_run(self)

    def _on_footer_accept(self) -> None:
        # Save only stores the parameters, so it needs no input files.
        if self._close_action == "run":
            self._try_accept()
        else:
            self.accept()

    def _on_file_input_field_activated(self) -> None:
        sender = self.sender()
        if isinstance(sender, TeMultiInputFilesWidget):
            self._focused_file_input = sender
            self._output_focused = False

    def apply_folder_from_panel(self, folder: Path) -> bool:
        """Folders-panel click: fill the output field if it last had focus."""
        if not self._output_focused or not Path(folder).is_dir():
            return False
        self.set_output_dir(Path(folder))
        return True

    def apply_path_from_panel(self, path: Path) -> bool:
        """Fill the last-focused path field from a Files-panel click."""
        if self._output_focused:
            return self.apply_folder_from_panel(path)
        if self._focused_file_input is None:
            return False
        if path.suffix.lower() not in self._FASTA_SUFFIXES:
            return False
        return self._focused_file_input.set_path_at_focus(path)

    def _try_accept(self) -> None:
        title = "Repeat Quantify"
        library = self.library()
        targets = self.targets()
        if library is None or not targets:
            QMessageBox.information(
                self, title,
                "Choose both a library FASTA and at least one genome FASTA "
                "(click a file in the Files panel, or use Browse).",
            )
            return
        missing = [p for p in [library, *targets] if not p.is_file()]
        if missing:
            QMessageBox.warning(
                self, title, "File not found:\n" + "\n".join(str(p) for p in missing[:5]),
            )
            return
        if library in targets:
            QMessageBox.warning(self, title, "Genome and library must be different files.")
            return
        if len(set(targets)) != len(targets):
            QMessageBox.warning(self, title, "The same genome is listed more than once.")
            return
        out = self.output_dir()
        if out is not None:
            if not out.is_absolute():
                QMessageBox.warning(
                    self, title, "Enter the full path of the output folder, or use Browse.",
                )
                return
            if out.exists() and not out.is_dir():
                QMessageBox.warning(
                    self, title, f"The output path is a file, not a folder:\n{out}",
                )
                return
        self.accept()

    def targets(self) -> list[Path]:
        """Genome FASTAs to quantify against, one report per genome."""
        return self._genome_input.paths()

    def target(self) -> Path | None:
        targets = self.targets()
        return targets[0] if targets else None

    def library(self) -> Path | None:
        paths = self._library_input.paths()
        return paths[0] if paths else None

    def engine(self) -> str:
        """Engine of the active tab (the last engine tab while on Options)."""
        return self._engine

    def params(self) -> dict:
        return self.options().params()
