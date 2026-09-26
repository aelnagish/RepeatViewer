"""Dialogs for tool-specific parameters (right-click / custom runs)."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import QEvent, QItemSelectionModel, QPoint, QSize, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QMouseEvent, QStandardItemModel
from PySide6.QtWidgets import (
    QAbstractItemView,
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
from repeatviewer.widgets.help_button import field_with_help, make_help_button  # noqa: E402


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


def make_process_output_box(
    parent: QWidget,
    *,
    tool_name: str,
    checked: bool = False,
) -> tuple[QGroupBox, QCheckBox]:
    """FastGA-style Options-tab group for the live process monitor window."""
    checkbox = QCheckBox(f"Show process window while {tool_name} runs", parent)
    checkbox.setChecked(bool(checked))
    checkbox.setToolTip(
        f"Open a window that first shows the {tool_name} command line, then "
        "streams progress and tool messages. Useful for long runs."
    )
    box = QGroupBox("Process output", parent)
    layout = QVBoxLayout(box)
    layout.addWidget(checkbox)
    return box, checkbox


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
from repeatviewer.backends.blast_genome_plot import format_genome_length_bp
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
    DEFAULT_ANNOTATION_DENSITY,
    DEFAULT_ANNOTATION_FEATURE_TYPES,
    CHROMOSOME_SORT_NATURAL,
    CHROMOSOME_SORT_SIZE_ASC,
    CHROMOSOME_SORT_SIZE_DESC,
    TeGenomeMapOptions,
    TeDivergencePlotOptions,
    MsaPlotJobOptions,
    RepeatQuantifyOptions,
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
    DEFAULT_CONSENSUS_AMBIGUOUS,
    DEFAULT_CONSENSUS_THRESHOLD,
    CONSENSUS_RIGHT_CLICK_RUN,
    CONSENSUS_RIGHT_CLICK_DIALOG,
    DEFAULT_EXTEND_LEFT,
    DEFAULT_EXTEND_RIGHT,
    EXTEND_RIGHT_CLICK_RUN,
    EXTEND_RIGHT_CLICK_DIALOG,
    FASTGA_RIGHT_CLICK_RUN,
    FASTGA_RIGHT_CLICK_DIALOG,
    DEFAULT_FLEXIDOT_LCS_SHADING,
    DEFAULT_FLEXIDOT_OUTPUT_FORMAT,
    DEFAULT_FLEXIDOT_WORDSIZE,
    DEFAULT_FLEXIDOT_MODES,
    FLEXIDOT_SEQUENCE_AA,
    FLEXIDOT_SEQUENCE_NUC,
    FLEXIDOT_RIGHT_CLICK_RUN,
    FLEXIDOT_RIGHT_CLICK_DIALOG,
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
    ConsensusOptions,
    normalize_consensus_ambiguous,
    configure_report_evalue_spinbox,
    DivergenceRowOptions,
    AlignOptions,
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
    ExtendOptions,
    FlexidotOptions,
    MAX_FLEXIDOT_WORDSIZE,
    MIN_FLEXIDOT_WORDSIZE,
    GapRowOptions,
    GappyColumnOptions,
    MsaCleanOptions,
    TeaAidOptions,
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
        title: str = "Alignment parameters",
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self._for_run = for_run
        opts = initial or AlignOptions()

        self._strategy = QComboBox(self)
        for label, value in (
            ("Auto (MAFFT chooses)", MAFFT_STRATEGY_AUTO),
            ("FFT-NS-1 (very fast)", MAFFT_STRATEGY_FFT_NS_1),
            ("FFT-NS-2 (fast)", MAFFT_STRATEGY_FFT_NS_2),
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
        self._threads.setToolTip("Number passed to MAFFT's --thread option.")

        self._retree = QSpinBox(self)
        self._retree.setRange(1, 100)
        self._retree.setValue(opts.retree)
        self._retree.setToolTip("Number of guide-tree rebuilds in Custom mode.")

        self._max_iterate = QSpinBox(self)
        self._max_iterate.setRange(0, 1000)
        self._max_iterate.setValue(opts.max_iterate)
        self._max_iterate.setToolTip("Maximum iterative-refinement cycles in Custom mode.")
        self._custom_retree = opts.retree
        self._custom_max_iterate = opts.max_iterate
        self._last_strategy: str | None = None

        self._gap_open = QDoubleSpinBox(self)
        self._gap_open.setRange(0.0, 20.0)
        self._gap_open.setDecimals(3)
        self._gap_open.setSingleStep(0.1)
        self._gap_open.setValue(opts.gap_open_penalty)
        self._gap_open.setToolTip("MAFFT --op gap-opening penalty; default 1.53.")

        self._offset = QDoubleSpinBox(self)
        self._offset.setRange(-10.0, 10.0)
        self._offset.setDecimals(3)
        self._offset.setSingleStep(0.05)
        self._offset.setValue(opts.offset)
        self._offset.setToolTip("MAFFT --ep offset value; default 0.0.")

        self._reorder = QCheckBox(
            "--reorder (group similar sequences in the output)", self
        )
        self._reorder.setChecked(opts.reorder)
        self._adjust_direction = QCheckBox(
            "--adjustdirectionaccurately (reverse-complement nucleotide sequences)",
            self,
        )
        self._adjust_direction.setChecked(opts.adjust_direction)
        self._adjust_direction.setToolTip("Uses MAFFT --adjustdirectionaccurately.")

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

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        form.addRow("MAFFT strategy:", self._strategy)
        form.addRow("Sequence type:", self._sequence_type)
        form.addRow("Threads:", self._threads)
        form.addRow("--retree (guide-tree rebuilds):", self._retree)
        form.addRow("--maxiterate (refinement cycles):", self._max_iterate)
        form.addRow("--op (gap-opening penalty):", self._gap_open)
        form.addRow("--ep (offset value):", self._offset)
        form.addRow("", self._reorder)
        form.addRow("", self._adjust_direction)
        self._mafft_combined = QCheckBox(
            "Combine all selected files into one MSA", self
        )
        self._mafft_combined.setChecked(opts.mafft_combined)
        self._mafft_combined.setToolTip(
            "Pool every sequence from the selected files into a single alignment, "
            "instead of aligning each file on its own."
        )
        form.addRow("", self._mafft_combined)

        mafft_tab = QWidget(self)
        mafft_tab.setLayout(form)

        self._muscle_mode = QComboBox(self)
        self._muscle_mode.addItem("Align (high accuracy)", MUSCLE_MODE_ALIGN)
        self._muscle_mode.addItem("Super5 (large datasets)", MUSCLE_MODE_SUPER5)
        muscle_mode_index = self._muscle_mode.findData(opts.muscle_mode)
        if muscle_mode_index >= 0:
            self._muscle_mode.setCurrentIndex(muscle_mode_index)
        self._muscle_mode.setFixedWidth(260)
        self._muscle_mode.setToolTip(
            "MUSCLE -align is the normal high-accuracy method; -super5 is intended "
            "for very large sequence sets."
        )
        self._muscle_threads = QSpinBox(self)
        self._muscle_threads.setRange(1, 64)
        self._muscle_threads.setValue(max(1, opts.muscle_num_threads))
        self._muscle_threads.setFixedWidth(260)
        self._muscle_threads.setToolTip("Number passed to MUSCLE's -threads option.")
        muscle_form = QFormLayout()
        muscle_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        muscle_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        muscle_form.setHorizontalSpacing(12)
        muscle_form.setVerticalSpacing(10)
        muscle_form.addRow("MUSCLE algorithm:", self._muscle_mode)
        muscle_form.addRow("Threads:", self._muscle_threads)
        self._muscle_combined = QCheckBox(
            "Combine all selected files into one MSA", self
        )
        self._muscle_combined.setChecked(opts.muscle_combined)
        self._muscle_combined.setToolTip(self._mafft_combined.toolTip())
        muscle_form.addRow("", self._muscle_combined)
        muscle_tab = QWidget(self)
        muscle_tab.setLayout(muscle_form)

        self._use_system_mafft = QCheckBox("Use MAFFT from the system PATH", self)
        self._use_system_mafft.setChecked(use_system_mafft)
        self._use_system_mafft.setToolTip(
            "Use the system mafft command instead of the bundled MAFFT."
        )
        self._use_system_muscle = QCheckBox("Use MUSCLE from the system PATH", self)
        self._use_system_muscle.setChecked(use_system_muscle)
        self._use_system_muscle.setToolTip(
            "Use the system muscle command instead of the bundled MUSCLE."
        )

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
        for label, value in (("MAFFT", ALIGNER_MAFFT), ("MUSCLE", ALIGNER_MUSCLE)):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_aligner == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
            rc_box_layout.addWidget(radio)
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[ALIGNER_MAFFT].setChecked(True)
        options_layout.addWidget(rc_box)

        exe_box = QGroupBox("Executables", self)
        exe_box_layout = QVBoxLayout(exe_box)
        exe_box_layout.addWidget(self._use_system_mafft)
        exe_box_layout.addWidget(self._use_system_muscle)
        options_layout.addWidget(exe_box)

        process_box, self._show_process_window = make_process_output_box(
            self,
            tool_name="Align",
            checked=bool(getattr(opts, "show_process_window", False)),
        )
        options_layout.addWidget(process_box)
        options_layout.addStretch(1)

        self._aligner_tabs = QTabWidget(self)
        self._aligner_tabs.addTab(mafft_tab, "MAFFT")
        self._aligner_tabs.addTab(muscle_tab, "MUSCLE")
        self._aligner_tabs.addTab(options_tab, "Options")
        self._options_tab_index = self._aligner_tabs.count() - 1

        self._close_action = "cancel"
        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore],
        )
        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(14)
        layout.addWidget(self._aligner_tabs)
        layout.addWidget(self._button_box)

        self._strategy.currentIndexChanged.connect(self._sync_controls)
        self._sequence_type.currentIndexChanged.connect(self._sync_controls)
        self._sync_controls()

    def should_run(self) -> bool:
        """True when the dialog was closed with Run (not Save or Cancel)."""
        return dialog_should_run(self)

    def selected_aligner(self) -> str:
        """The aligner of the active tab (MAFFT/MUSCLE); falls back to right-click."""
        index = self._aligner_tabs.currentIndex()
        if index == 0:
            return ALIGNER_MAFFT
        if index == 1:
            return ALIGNER_MUSCLE
        return self.right_click_aligner()

    def right_click_aligner(self) -> str:
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                return value
        return ALIGNER_MAFFT

    def combined_for(self, aligner: str) -> bool:
        """Whether the given aligner should combine into one MSA."""
        if aligner == ALIGNER_MUSCLE:
            return self._muscle_combined.isChecked()
        return self._mafft_combined.isChecked()

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
        return AlignOptions(
            combined_aligner=self.right_click_aligner(),
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
            mafft_combined=self._mafft_combined.isChecked(),
            muscle_combined=self._muscle_combined.isChecked(),
            right_click_aligner=self.right_click_aligner(),
            show_process_window=self._show_process_window.isChecked(),
        )

    def use_system_mafft(self) -> bool:
        return self._use_system_mafft.isChecked()

    def use_system_muscle(self) -> bool:
        return self._use_system_muscle.isChecked()

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
        self._mafft_combined.setChecked(defaults.mafft_combined)

        muscle_idx = self._muscle_mode.findData(defaults.muscle_mode)
        if muscle_idx >= 0:
            self._muscle_mode.setCurrentIndex(muscle_idx)
        self._muscle_threads.setValue(max(1, defaults.muscle_num_threads))
        self._muscle_combined.setChecked(defaults.muscle_combined)

        for value, radio in self._right_click_buttons.items():
            radio.setChecked(value == defaults.right_click_aligner)
        self._show_process_window.setChecked(defaults.show_process_window)
        # Default: bundled MAFFT / MUSCLE binaries, not system PATH.
        self._use_system_mafft.setChecked(False)
        self._use_system_muscle.setChecked(False)
        self._sync_controls()


class NewSequenceDialog(QDialog):
    """Type or paste a sequence and save it as a FASTA."""

    def __init__(self, *, parent=None) -> None:
        super().__init__(parent)
        from repeatviewer.backends.sequence_compose import MOLECULE_LABELS

        self.setWindowTitle("New Sequence")
        self.resize(720, 520)

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
        layout.addLayout(form)
        layout.addWidget(self._summary)
        layout.addWidget(self._buttons)

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

    def _sync(self) -> None:
        from repeatviewer.backends.sequence_compose import check_sequence

        check = check_sequence(self._text.toPlainText(), str(self._molecule.currentData()))
        parts = [f"{check.length:,} bp" if check.length else "empty"]
        if check.invalid:
            parts.append(f"unexpected characters: {check.invalid}")
        self._summary.setText("  ·  ".join(parts))
        ok = bool(check.length) and bool(self._name.text().strip())
        self._buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(ok)

    def values(self) -> dict:
        return {
            "name": self._name.text().strip(),
            "sequence": self._text.toPlainText(),
            "description": self._description.text().strip(),
            "organism": self._organism.text().strip(),
            "molecule": str(self._molecule.currentData()),
        }


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

    Two ways in: tick boxes (fine for a handful of monomers) or an ID list
    (the only workable option for thousands of reads). ``records`` may be
    ``None`` when the source is too large to list cheaply -- the tick-box tab
    is then replaced by a note.
    """

    #: Emitted with the parsed id list when the user presses Check.
    check_requested = Signal(object)

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
        self.resize(600, 520)
        self._source_name = source_name

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
        layout.addWidget(self._summary)
        layout.addWidget(self._buttons)

        self._tabs.currentChanged.connect(self._sync)
        # Opened with ids already in hand (from a hit table): start on that tab.
        if initial_ids.strip() or records is None:
            self._tabs.setCurrentIndex(1)
        self._sync()

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
        self._filter.setPlaceholderText("Search sequences…")
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
        top.addWidget(QLabel(f"Sequences in {self._source_name}:", page))
        top.addStretch(1)
        top.addWidget(select_all)
        top.addWidget(select_none)

        hint = QLabel(
            "Tick boxes to add or remove · Click a name to select · "
            "Shift+click for a range · ⌘/Ctrl+click to add or remove · "
            "drag to select several",
            page,
        )
        hint.setStyleSheet("color: palette(mid);")
        hint.setWordWrap(True)

        layout.addLayout(top)
        layout.addWidget(self._filter)
        layout.addWidget(self._list, stretch=1)
        layout.addWidget(hint)
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

    def single_file(self) -> bool:
        return self._single_file.isChecked()

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


class FastGaParamsDialog(QDialog):
    """FastGA whole-genome alignment parameters (BLAST-style form)."""

    def __init__(
        self,
        *,
        initial=None,
        file_count: int = 1,
        title: str = "FastGA parameters",
        parent=None,
    ) -> None:
        super().__init__(parent)
        from repeatviewer.backends.fastga_dotplot import FORMAT_LABELS
        from repeatviewer.job_options import FastGaOptions

        self.setWindowTitle(title)
        opts = initial or FastGaOptions()

        self._format = QComboBox(self)
        for label, value in FORMAT_LABELS:
            self._format.addItem(label, value)
        idx = self._format.findData(opts.output_format)
        if idx >= 0:
            self._format.setCurrentIndex(idx)
        self._format.setToolTip(
            "Output format: PAF and PSL are text files streamed by FastGA "
            "(-paf / -pafx / -pafm / -pafs / -pafS / -psl); -1 writes a binary "
            ".1aln (ONEcode) file."
        )

        self._threads = QSpinBox(self)
        self._threads.setRange(1, 128)
        self._threads.setValue(max(1, opts.threads))
        self._threads.setToolTip("FastGA -T: number of threads to use.")

        self._identity = QDoubleSpinBox(self)
        self._identity.setRange(0.5, 1.0)
        self._identity.setSingleStep(0.05)
        self._identity.setDecimals(2)
        self._identity.setValue(opts.min_identity)
        self._identity.setToolTip("FastGA -i: minimum alignment identity.")

        self._min_length = QSpinBox(self)
        # Below ~100 FastGA fails: an alignment must contain a 40 bp seed and
        # chains must cover at least 85 bases.
        self._min_length.setRange(100, 1_000_000)
        self._min_length.setSingleStep(50)
        self._min_length.setValue(max(100, opts.min_length))
        self._min_length.setToolTip(
            "FastGA -l: shortest alignment to report, in bp. Cannot go below 100 "
            "-- FastGA rejects limits shorter than its seed length.\n"
            "Raise this on satellite arrays, where short matches can swamp "
            "the output."
        )

        self._seed_freq = QSpinBox(self)
        self._seed_freq.setRange(1, 100_000)
        self._seed_freq.setSingleStep(50)
        self._seed_freq.setValue(max(1, opts.seed_freq))
        self._seed_freq.setToolTip(
            "FastGA -f: adaptive seed count cutoff. Seeds occurring more often "
            "than this are discarded as repetitive before alignment.\n"
            "FastGA's own default is 10, which removes satellite arrays "
            "entirely. Raise it to see repeat structure; higher values are slower."
        )

        self._chain_coverage = QSpinBox(self)
        self._chain_coverage.setRange(1, 1_000_000)
        self._chain_coverage.setValue(max(1, opts.chain_coverage))
        self._chain_coverage.setToolTip(
            "FastGA -c: minimum seed chain coverage (bases) in both genomes. "
            "FastGA default: 85."
        )

        self._chain_gap = QSpinBox(self)
        self._chain_gap.setRange(1, 10_000_000)
        self._chain_gap.setSingleStep(100)
        self._chain_gap.setValue(max(1, opts.chain_gap))
        self._chain_gap.setToolTip(
            "FastGA -s: maximum spacing between consecutive seeds when starting "
            "a new seed chain. FastGA default: 1000."
        )

        self._symmetric = QCheckBox("Symmetric seeding (-S)", self)
        self._symmetric.setChecked(opts.symmetric)
        self._symmetric.setToolTip(
            "FastGA -S: seed adaptamers from both genomes so the result does "
            "not depend on argument order. Useful for repeat structure; slower."
        )

        self._keep_index = QCheckBox("Keep genome index (-k)", self)
        self._keep_index.setChecked(opts.keep_index)
        self._keep_index.setToolTip(
            "FastGA -k: keep generated .1gdb / .gix files beside the FASTA so "
            "later runs skip the index build."
        )

        self._no_soft_mask = QCheckBox("Do not soft mask (-N)", self)
        self._no_soft_mask.setChecked(opts.no_soft_mask)
        self._no_soft_mask.setToolTip("FastGA -N: do not soft mask.")

        self._pairwise = QCheckBox("Do each pair comparison", self)
        self._pairwise.setChecked(bool(getattr(opts, "pairwise", False)))
        self._pairwise.setToolTip(
            "Off (default): each selected FASTA is compared against itself, "
            "one file at a time.\n"
            "On: compare every unique pair among the selected files and write "
            "all alignments into one file named "
            "<first_selected>_multi_FastGA.<ext>."
        )

        for field in (
            self._format,
            self._threads,
            self._identity,
            self._min_length,
            self._seed_freq,
            self._chain_coverage,
            self._chain_gap,
        ):
            field.setFixedWidth(280)

        self._file_count = int(file_count)
        self._mode = QLabel(self)
        self._mode.setWordWrap(True)
        self._mode.setStyleSheet("color: palette(mid);")
        self._pairwise.toggled.connect(self._update_mode_text)
        self._update_mode_text()

        params_tab = QWidget(self)
        form = QFormLayout(params_tab)
        form.setContentsMargins(16, 14, 16, 14)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        form.addRow(self._mode)
        form.addRow("Output format:", self._format)
        form.addRow("-T (threads):", self._threads)
        form.addRow("-i (min identity):", self._identity)
        form.addRow("-l (min length):", self._min_length)
        form.addRow("-f (seed cutoff):", self._seed_freq)
        form.addRow("-c (chain coverage):", self._chain_coverage)
        form.addRow("-s (chain gap):", self._chain_gap)
        form.addRow("", self._symmetric)
        form.addRow("", self._keep_index)
        form.addRow("", self._no_soft_mask)
        form.addRow("", self._pairwise)

        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("Run FastGA with saved parameters", FASTGA_RIGHT_CLICK_RUN),
            ("Open this dialog", FASTGA_RIGHT_CLICK_DIALOG),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[FASTGA_RIGHT_CLICK_RUN].setChecked(True)

        rc_box = QGroupBox("Right-click action", self)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(
            QLabel("Right-clicking the FastGA button does:", self)
        )
        for radio in self._right_click_buttons.values():
            rc_box_layout.addWidget(radio)

        process_box, self._show_process_window = make_process_output_box(
            self,
            tool_name="FastGA",
            checked=bool(getattr(opts, "show_process_window", False)),
        )

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)
        options_layout.addWidget(rc_box)
        options_layout.addWidget(process_box)
        options_layout.addStretch(1)

        self._tabs = QTabWidget(self)
        self._tabs.addTab(params_tab, "Parameters")
        self._tabs.addTab(options_tab, "Options")

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore],
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(14)
        layout.addWidget(self._tabs)
        layout.addWidget(self._button_box)

    def _update_mode_text(self) -> None:
        pairwise = self._pairwise.isChecked()
        n = self._file_count
        if pairwise:
            if n >= 2:
                pair_count = n * (n - 1) // 2
                text = (
                    f"{n} files selected: every unique pair will be compared "
                    f"({pair_count} pair(s)). Results are merged into one file "
                    "named after the first selection with _multi_FastGA."
                )
            else:
                text = ""
        elif n >= 1:
            text = (
                f"{n} file(s) selected: each will be compared against itself "
                "one by one (self-comparison)."
            )
        else:
            text = (
                "Self-comparison (default): select one or more FASTA files. "
                "Each is compared against itself, one by one. Enable "
                "Do each pair comparison to compare every unique pair among "
                "the selected files."
            )
        self._mode.setText(text)

    def _restore_defaults(self) -> None:
        from repeatviewer.job_options import FastGaOptions

        opts = FastGaOptions()
        idx = self._format.findData(opts.output_format)
        if idx >= 0:
            self._format.setCurrentIndex(idx)
        self._threads.setValue(max(1, opts.threads))
        self._identity.setValue(opts.min_identity)
        self._min_length.setValue(max(100, opts.min_length))
        self._seed_freq.setValue(opts.seed_freq)
        self._chain_coverage.setValue(opts.chain_coverage)
        self._chain_gap.setValue(opts.chain_gap)
        self._symmetric.setChecked(opts.symmetric)
        self._keep_index.setChecked(opts.keep_index)
        self._no_soft_mask.setChecked(opts.no_soft_mask)
        self._pairwise.setChecked(opts.pairwise)
        self._show_process_window.setChecked(opts.show_process_window)
        self._right_click_buttons[FASTGA_RIGHT_CLICK_RUN].setChecked(True)
        self._update_mode_text()

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def options(self):
        from repeatviewer.job_options import FastGaOptions

        right_click_action = FASTGA_RIGHT_CLICK_RUN
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click_action = value
                break

        return FastGaOptions(
            output_format=str(self._format.currentData()),
            threads=int(self._threads.value()),
            min_identity=float(self._identity.value()),
            min_length=int(self._min_length.value()),
            seed_freq=int(self._seed_freq.value()),
            chain_coverage=int(self._chain_coverage.value()),
            chain_gap=int(self._chain_gap.value()),
            symmetric=self._symmetric.isChecked(),
            keep_index=self._keep_index.isChecked(),
            no_soft_mask=self._no_soft_mask.isChecked(),
            pairwise=self._pairwise.isChecked(),
            temp_dir="",
            right_click_action=right_click_action,
            show_process_window=self._show_process_window.isChecked(),
        )


class GenomeBlastDatabaseRow(QWidget):
    """Genome combo matching Blast Database: current genome + nucleotide DBs + Manage."""

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        project: RepeatViewerProject | None = None,
        project_genome: Path | None = None,
        selected_path: str | None = None,
        manage_title: str = "Manage BLAST databases",
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
        self._combo.setMinimumWidth(260)
        self._combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self._combo.setMinimumContentsLength(28)
        self.refresh(selected_path=selected_path)

        manage = QPushButton("Add/Manage databases", self)
        manage.clicked.connect(self._manage_databases)

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        # Keep the button next to the combo (do not stretch the combo across the row).
        layout.addWidget(self._combo)
        layout.addWidget(manage)
        layout.addStretch(1)

    def selected_path(self) -> str | None:
        return self._combo.currentData() or None

    def reset_to_default(self) -> None:
        if self._combo.count() > 0:
            self._combo.setCurrentIndex(0)

    def refresh(self, *, selected_path: str | None = None) -> None:
        """Fill combo like Blast Database (current genome + nucleotide DBs)."""
        from repeatviewer.database_settings import reload_blast_databases_into_project

        target = selected_path
        if target is None and self._combo.count() > 0:
            target = self._combo.currentData()

        if self._project is not None:
            reload_blast_databases_into_project(self._project)

        self._combo.clear()
        if self._project_genome is not None:
            self._combo.addItem(
                f"Current genome ({self._project_genome.name})",
                None,
            )
        else:
            self._combo.addItem("Current genome (not configured)", None)

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

        if target:
            resolved = str(Path(str(target)).expanduser().resolve())
            index = self._combo.findData(resolved)
            if index < 0 and self._project_genome is not None:
                if resolved == str(self._project_genome):
                    index = 0
            if index < 0:
                self._combo.addItem(Path(resolved).name, resolved)
                index = self._combo.count() - 1
            self._combo.setCurrentIndex(index)
        else:
            self._combo.setCurrentIndex(0)

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
        dialog.finished.connect(self._on_manager_finished)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

    def _on_manager_finished(self, _result: int = 0) -> None:
        self._blast_manager_dialog = None
        self.refresh(selected_path=self._combo.currentData())


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

        self._left = QSpinBox()
        self._left.setRange(0, 500_000)
        self._left.setValue(opts.left)
        self._left.setSuffix(" bp")
        self._right = QSpinBox()
        self._right.setRange(0, 500_000)
        self._right.setValue(opts.right)
        self._right.setSuffix(" bp")

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
            ("Open this dialog", EXTEND_RIGHT_CLICK_DIALOG),
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

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)
        options_layout.addWidget(rc_box)
        options_layout.addStretch(1)

        self._tabs = QTabWidget(self)
        params_tab = QWidget(self)
        params_form = QFormLayout(params_tab)
        params_form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        params_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        params_form.setHorizontalSpacing(12)
        params_form.setVerticalSpacing(10)
        if self._pan_genome_names:
            params_form.addRow("Pan-genome:", self._pan_genome)
        else:
            params_form.addRow("Genome:", self._genome_row)
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
        layout.addWidget(self._tabs)
        layout.addWidget(self._button_box)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _restore_defaults(self) -> None:
        self._left.setValue(DEFAULT_EXTEND_LEFT)
        self._right.setValue(DEFAULT_EXTEND_RIGHT)
        self._genome_row.reset_to_default()
        self._right_click_buttons[EXTEND_RIGHT_CLICK_RUN].setChecked(True)

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
            left=self._left.value(),
            right=self._right.value(),
            pan_genome_name=pan_name,
            genome_path=genome_path,
            right_click_action=right_click_action,
        )


class TsdSearchParamsDialog(QDialog):
    def __init__(
        self,
        *,
        initial: TsdSearchOptions | None = None,
        title: str = "TSD search parameters",
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        opts = initial or TsdSearchOptions()

        self._left = QSpinBox()
        self._left.setRange(0, 10_000_000)
        self._left.setValue(opts.left_boundary)
        self._left.setToolTip("MSA column index for the left TE boundary.")

        self._right = QSpinBox()
        self._right.setRange(0, 10_000_000)
        self._right.setValue(opts.right_boundary)
        self._right.setToolTip("MSA column index for the right TE boundary.")

        self._radius = QSpinBox()
        self._radius.setRange(1, 10_000)
        self._radius.setValue(opts.nt_radius)
        self._radius.setSuffix(" nt")
        self._radius.setToolTip("Search window size on each side of the boundary.")

        self._max_sim = QDoubleSpinBox()
        self._max_sim.setRange(0.0, 1.0)
        self._max_sim.setDecimals(2)
        self._max_sim.setSingleStep(0.05)
        self._max_sim.setValue(opts.max_sim_threshold)

        self._nt_threshold = QDoubleSpinBox()
        self._nt_threshold.setRange(0.0, 1.0)
        self._nt_threshold.setDecimals(2)
        self._nt_threshold.setSingleStep(0.05)
        self._nt_threshold.setValue(opts.nt_threshold)

        self._cluster_eps = QDoubleSpinBox()
        self._cluster_eps.setRange(0.1, 10.0)
        self._cluster_eps.setDecimals(1)
        self._cluster_eps.setSingleStep(0.1)
        self._cluster_eps.setValue(opts.cluster_eps)

        form = QFormLayout()
        form.addRow("Left boundary (MSA column):", self._left)
        form.addRow("Right boundary (MSA column):", self._right)
        form.addRow("Search radius:", self._radius)
        form.addRow("Max TSD similarity:", self._max_sim)
        form.addRow("Poly-N threshold:", self._nt_threshold)
        form.addRow("Cluster distance (eps):", self._cluster_eps)

        self._button_box = make_cancel_save_run_bar(self)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(self._button_box)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def options(self) -> TsdSearchOptions:
        return TsdSearchOptions(
            left_boundary=self._left.value(),
            right_boundary=self._right.value(),
            nt_radius=self._radius.value(),
            max_sim_threshold=self._max_sim.value(),
            nt_threshold=self._nt_threshold.value(),
            cluster_eps=self._cluster_eps.value(),
        )


_KEEP_CURRENT = object()


class BlastParamsDialog(QDialog):
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

        # Always show the latest Manage-BLAST list (and current genome name).
        from repeatviewer.database_settings import reload_blast_databases_into_project

        reload_blast_databases_into_project(self._project)

        self._tabs = QTabWidget()
        self._blastn_database = QComboBox()
        self._blastp_database = QComboBox()
        for combo in (self._blastn_database, self._blastp_database):
            combo.setMinimumWidth(260)
            combo.setSizeAdjustPolicy(
                QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
            )
            combo.setMinimumContentsLength(28)

        self._evalue = QDoubleSpinBox()
        configure_report_evalue_spinbox(self._evalue)
        self._evalue.setValue(normalize_blast_evalue(opts.e_value))

        self._threads = QSpinBox()
        self._threads.setRange(1, 128)
        self._threads.setValue(max(1, opts.num_threads))
        self._threads.setToolTip("BLASTN -num_threads: CPU threads used for this search.")

        self._word_size = QSpinBox()
        self._word_size.setRange(4, 10_000)
        self._word_size.setValue(opts.word_size)
        self._word_size.setToolTip(
            "BLASTN -word_size: length of the exact seed match used to begin a search. "
            "Smaller values are more sensitive but slower."
        )

        self._gap_open = QSpinBox()
        self._gap_open.setRange(0, 1000)
        self._gap_open.setValue(opts.gap_open)
        self._gap_open.setToolTip("BLASTN -gapopen: cost for opening a gap.")

        self._gap_extend = QSpinBox()
        self._gap_extend.setRange(0, 1000)
        self._gap_extend.setValue(opts.gap_extend)
        self._gap_extend.setToolTip("BLASTN -gapextend: cost for extending an existing gap.")

        self._penalty = QSpinBox()
        self._penalty.setRange(-1000, 0)
        self._penalty.setValue(opts.penalty)
        self._penalty.setToolTip("BLASTN -penalty: score assigned to a nucleotide mismatch.")

        self._percent_identity = QDoubleSpinBox()
        self._percent_identity.setRange(0.0, 100.0)
        self._percent_identity.setDecimals(1)
        self._percent_identity.setSuffix(" %")
        self._percent_identity.setValue(opts.percent_identity)
        self._percent_identity.setToolTip(
            "BLASTN -perc_identity: discard HSPs below this percent identity; 0 disables the filter."
        )

        self._query_coverage = QDoubleSpinBox()
        self._query_coverage.setRange(0.0, 100.0)
        self._query_coverage.setDecimals(1)
        self._query_coverage.setSuffix(" %")
        self._query_coverage.setValue(opts.query_coverage_hsp)
        self._query_coverage.setToolTip(
            "BLASTN -qcov_hsp_perc: discard HSPs covering less than this percentage "
            "of the query; 0 disables the filter."
        )

        self._blastn_write_text = QCheckBox("Export BLAST text file (*.txt)")
        self._blastn_write_text.setChecked(opts.write_text)
        self._blastn_write_gff = QCheckBox("Export GFF file (*.gff)")
        self._blastn_write_gff.setChecked(opts.write_gff)

        self._blastp_evalue = QDoubleSpinBox()
        configure_report_evalue_spinbox(self._blastp_evalue)
        self._blastp_evalue.setValue(normalize_blast_evalue(opts.blastp_e_value))

        self._blastp_threads = QSpinBox()
        self._blastp_threads.setRange(1, 128)
        self._blastp_threads.setValue(max(1, opts.blastp_num_threads))
        self._blastp_threads.setToolTip("BLASTP -num_threads: CPU threads used for this search.")

        self._blastp_word_size = QSpinBox()
        self._blastp_word_size.setRange(2, 10_000)
        self._blastp_word_size.setValue(opts.blastp_word_size)
        self._blastp_word_size.setToolTip(
            "BLASTP -word_size: length of the protein seed used to begin a search."
        )

        self._blastp_gap_open = QSpinBox()
        self._blastp_gap_open.setRange(0, 1000)
        self._blastp_gap_open.setValue(opts.blastp_gap_open)
        self._blastp_gap_open.setToolTip("BLASTP -gapopen: cost for opening a gap.")

        self._blastp_gap_extend = QSpinBox()
        self._blastp_gap_extend.setRange(0, 1000)
        self._blastp_gap_extend.setValue(opts.blastp_gap_extend)
        self._blastp_gap_extend.setToolTip("BLASTP -gapextend: cost for extending a gap.")

        self._blastp_query_coverage = QDoubleSpinBox()
        self._blastp_query_coverage.setRange(0.0, 100.0)
        self._blastp_query_coverage.setDecimals(1)
        self._blastp_query_coverage.setSuffix(" %")
        self._blastp_query_coverage.setValue(opts.blastp_query_coverage_hsp)
        self._blastp_query_coverage.setToolTip(
            "BLASTP -qcov_hsp_perc: discard HSPs below this query coverage; 0 disables the filter."
        )
        self._blastp_write_text = QCheckBox("Export BLAST text file (*.txt)")
        self._blastp_write_text.setChecked(opts.blastp_write_text)
        self._blastp_write_gff = QCheckBox("Export GFF file (*.gff)")
        self._blastp_write_gff.setChecked(opts.blastp_write_gff)

        # --- blastx: nucleotide query translated against a protein database ---
        self._blastx_database = QComboBox()
        self._blastx_database.setMinimumWidth(260)
        self._blastx_database.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon
        )
        self._blastx_database.setMinimumContentsLength(28)

        self._blastx_evalue = QDoubleSpinBox()
        configure_report_evalue_spinbox(self._blastx_evalue)
        self._blastx_evalue.setValue(normalize_blast_evalue(opts.blastx_e_value))

        self._blastx_threads = QSpinBox()
        self._blastx_threads.setRange(1, 128)
        self._blastx_threads.setValue(max(1, opts.blastx_num_threads))
        self._blastx_threads.setToolTip("BLASTX -num_threads: CPU threads used for this search.")

        self._blastx_word_size = QSpinBox()
        self._blastx_word_size.setRange(2, 10_000)
        self._blastx_word_size.setValue(opts.blastx_word_size)
        self._blastx_word_size.setToolTip(
            "BLASTX -word_size: length of the protein seed used to begin a search."
        )

        self._blastx_gap_open = QSpinBox()
        self._blastx_gap_open.setRange(0, 1000)
        self._blastx_gap_open.setValue(opts.blastx_gap_open)
        self._blastx_gap_open.setToolTip("BLASTX -gapopen: cost for opening a gap.")

        self._blastx_gap_extend = QSpinBox()
        self._blastx_gap_extend.setRange(0, 1000)
        self._blastx_gap_extend.setValue(opts.blastx_gap_extend)
        self._blastx_gap_extend.setToolTip("BLASTX -gapextend: cost for extending a gap.")

        self._blastx_query_coverage = QDoubleSpinBox()
        self._blastx_query_coverage.setRange(0.0, 100.0)
        self._blastx_query_coverage.setDecimals(1)
        self._blastx_query_coverage.setSuffix(" %")
        self._blastx_query_coverage.setValue(opts.blastx_query_coverage_hsp)
        self._blastx_query_coverage.setToolTip(
            "BLASTX -qcov_hsp_perc: discard HSPs below this query coverage; 0 disables the filter."
        )
        self._blastx_write_text = QCheckBox("Export BLAST text file (*.txt)")
        self._blastx_write_text.setChecked(opts.blastx_write_text)
        self._blastx_write_gff = QCheckBox("Export GFF file (*.gff)")
        self._blastx_write_gff.setChecked(opts.blastx_write_gff)

        self._max_msa = QSpinBox()
        self._max_msa.setRange(10, 5000)
        self._max_msa.setValue(opts.max_msa_lines)
        self._max_msa_help = (
            "Set the maximum number of sequences to be included in a multiple sequence "
            "alignment. Default: 100"
        )
        self._max_msa.setToolTip(self._max_msa_help)

        self._top_msa = QSpinBox()
        self._top_msa.setRange(10, 5000)
        self._top_msa.setMaximum(self._max_msa.value())
        self._top_msa.setValue(opts.top_msa_lines)
        self._max_msa.valueChanged.connect(self._top_msa.setMaximum)
        self._top_msa_help = (
            "If the sequence number of the multiple sequence alignment (MSA) is greater "
            "than <max_msa_lines>, TEtrimmer first sorts sequences by length and chooses "
            "<top_msa_lines> sequences. It then randomly selects sequences from all remaining "
            "BLAST hits until <max_msa_lines> sequences are found for the multiple sequence "
            "alignment. Default: 100"
        )
        self._top_msa.setToolTip(self._top_msa_help)

        self._min_bit_score = QDoubleSpinBox()
        self._min_bit_score.setRange(0.0, 1_000_000.0)
        self._min_bit_score.setDecimals(1)
        self._min_bit_score.setSingleStep(1.0)
        self._min_bit_score.setValue(opts.minimum_bit_score)
        self._min_bit_score.setToolTip("Minimum bit score for hits written to GFF.")

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
            self._max_msa,
            self._top_msa,
            self._min_bit_score,
        ):
            field.setFixedWidth(260)
            field.setAlignment(Qt.AlignmentFlag.AlignRight)

        blastn_tab = QWidget()
        blastn_form = self._parameter_form(blastn_tab)
        blastn_form.addRow("Database:", self._database_row(self._blastn_database))
        blastn_form.addRow("-evalue:", self._evalue)
        blastn_form.addRow("Threads:", self._threads)
        blastn_form.addRow("-word_size:", self._word_size)
        blastn_form.addRow("-gapopen:", self._gap_open)
        blastn_form.addRow("-gapextend:", self._gap_extend)
        blastn_form.addRow("-penalty:", self._penalty)
        blastn_form.addRow("-perc_identity:", self._percent_identity)
        blastn_form.addRow("-qcov_hsp_perc:", self._query_coverage)
        blastn_form.addRow("", self._blastn_write_text)
        blastn_form.addRow("", self._blastn_write_gff)
        # The system-PATH toggle applies to all three programs, so it lives once
        # on the Options tab rather than being repeated on every program tab.
        self._use_system_blast = QCheckBox("Use BLAST+ from the system PATH")
        self._use_system_blast.setChecked(project.use_system_blast)
        self._use_system_blast.setToolTip(
            "Use BLAST+ programs installed on the system PATH instead of the bundled programs."
        )
        self._tabs.addTab(blastn_tab, "BLASTn")

        blastp_tab = QWidget()
        blastp_form = self._parameter_form(blastp_tab)
        blastp_form.addRow("Database:", self._database_row(self._blastp_database))
        blastp_form.addRow("-evalue:", self._blastp_evalue)
        blastp_form.addRow("Threads:", self._blastp_threads)
        blastp_form.addRow("-word_size:", self._blastp_word_size)
        blastp_form.addRow("-gapopen:", self._blastp_gap_open)
        blastp_form.addRow("-gapextend:", self._blastp_gap_extend)
        blastp_form.addRow("-qcov_hsp_perc:", self._blastp_query_coverage)
        blastp_form.addRow("", self._blastp_write_text)
        blastp_form.addRow("", self._blastp_write_gff)
        self._tabs.addTab(blastp_tab, "BLASTp")

        blastx_tab = QWidget()
        blastx_form = self._parameter_form(blastx_tab)
        blastx_form.addRow("Database:", self._database_row(self._blastx_database))
        blastx_form.addRow("-evalue:", self._blastx_evalue)
        blastx_form.addRow("Threads:", self._blastx_threads)
        blastx_form.addRow("-word_size:", self._blastx_word_size)
        blastx_form.addRow("-gapopen:", self._blastx_gap_open)
        blastx_form.addRow("-gapextend:", self._blastx_gap_extend)
        blastx_form.addRow("-qcov_hsp_perc:", self._blastx_query_coverage)
        blastx_form.addRow("", self._blastx_write_text)
        blastx_form.addRow("", self._blastx_write_gff)
        self._tabs.addTab(blastx_tab, "BLASTx")

        # --- Options: right-click program + shared executable setting ---
        options_tab = QWidget()
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)

        rc_box = QGroupBox("Right-click action", self)
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
            ("blastn", BLAST_PROGRAM_BLASTN),
            ("blastp", BLAST_PROGRAM_BLASTP),
            ("blastx", BLAST_PROGRAM_BLASTX),
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
        exe_box_layout = QVBoxLayout(exe_box)
        exe_box_layout.addWidget(self._use_system_blast)
        options_layout.addWidget(exe_box)

        process_box, self._show_process_window = make_process_output_box(
            self,
            tool_name="BLAST",
            checked=bool(getattr(opts, "show_process_window", False)),
        )
        options_layout.addWidget(process_box)

        options_layout.addStretch(1)
        self._tabs.addTab(options_tab, "Options")

        self._populate_database_combos(
            blastn_name=opts.blast_database_name,
            blastp_name=opts.blastp_database_name,
            blastx_name=opts.blastx_database_name,
        )
        program_tab = {
            BLAST_PROGRAM_BLASTP: 1,
            BLAST_PROGRAM_BLASTX: 2,
        }.get(opts.program, 0)
        self._tabs.setCurrentIndex(program_tab)

        self._post_box = QGroupBox("Post-processing")
        post_layout = QVBoxLayout(self._post_box)
        post_form = QFormLayout()
        post_form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        post_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        post_form.setHorizontalSpacing(12)
        post_form.setVerticalSpacing(10)
        post_form.addRow(
            "Max MSA lines:",
            self._field_with_help(self._max_msa, "Max MSA lines", self._max_msa_help),
        )
        post_form.addRow(
            "Top longest lines:",
            self._field_with_help(self._top_msa, "Top longest lines", self._top_msa_help),
        )
        post_form.addRow("Minimum bit score (GFF):", self._min_bit_score)
        post_layout.addLayout(post_form)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore],
        )

        layout = QVBoxLayout(self)
        layout.addWidget(self._tabs)
        layout.addWidget(self._post_box)
        layout.addWidget(self._button_box)
        self._tabs.currentChanged.connect(self._update_postprocessing_visibility)
        self._update_postprocessing_visibility(self._tabs.currentIndex())

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        # Refresh in case Manage / Add ran while this dialog was constructed,
        # and pick up a newly defined genome for the "Current genome" row.
        self._refresh_database_combos()
        # After dialog polish's adjustSize(), shrink to the active tab so the
        # wider Options page does not leave blank space on BLASTn/p/x.
        QTimer.singleShot(0, self._fit_to_content)

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
    def _parameter_form(parent: QWidget) -> QFormLayout:
        form = QFormLayout(parent)
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        return form

    def _database_row(self, combo: QComboBox) -> QWidget:
        manage = QPushButton("Add/Manage databases")
        manage.clicked.connect(self._manage_databases)
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.addWidget(combo, 1)
        layout.addWidget(manage)
        return row

    def _field_with_help(self, field: QWidget, title: str, help_text: str) -> QWidget:
        return field_with_help(field, title, help_text)

    def _update_postprocessing_visibility(self, tab_index: int) -> None:
        self._post_box.setVisible(tab_index == 0)
        if self.isVisible():
            QTimer.singleShot(0, self._fit_to_content)

    def _fit_to_content(self) -> None:
        """Size width/height to the active tab (not the widest tab page)."""
        page = self._tabs.currentWidget()
        if page is None or self.layout() is None:
            return
        page.ensurePolished()
        page.adjustSize()
        page_hint = page.sizeHint()
        post_w = self._post_box.sizeHint().width() if self._post_box.isVisible() else 0
        post_h = self._post_box.sizeHint().height() if self._post_box.isVisible() else 0
        tab_chrome = max(24, self._tabs.tabBar().sizeHint().height()) + 8
        buttons_h = self._button_box.sizeHint().height()
        buttons_w = self._button_box.sizeHint().width()
        margins = self.layout().contentsMargins()
        spacing = max(0, self.layout().spacing())
        gaps = 2 if self._post_box.isVisible() else 1
        width = max(
            420,
            page_hint.width() + margins.left() + margins.right() + 24,
            post_w + margins.left() + margins.right() + 24,
            buttons_w + margins.left() + margins.right(),
        )
        height = (
            page_hint.height()
            + tab_chrome
            + post_h
            + buttons_h
            + margins.top()
            + margins.bottom()
            + spacing * gaps
            + 4
        )
        self.resize(width, height)

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
            if self._project.has_genome_configured() and self._project.genome_file is not None:
                combo.addItem(f"Current genome ({self._project.genome_file.name})", None)
            else:
                combo.addItem("Current genome (not configured)", None)
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

        # Non-modal manager (and Add dialog) so Files-panel clicks can fill FASTA.
        existing = getattr(self, "_blast_manager_dialog", None)
        if existing is not None and existing.isVisible():
            existing.raise_()
            existing.activateWindow()
            return
        dialog = BlastDatabaseManagerDialog(self._project, parent=self.window())
        self._blast_manager_dialog = dialog
        dialog.finished.connect(self._on_blast_manager_finished)
        dialog.show()
        dialog.raise_()
        dialog.activateWindow()

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
        self._blastn_write_gff.setChecked(defaults.write_gff)
        self._blastp_evalue.setValue(defaults.blastp_e_value)
        self._blastp_threads.setValue(defaults.blastp_num_threads)
        self._blastp_word_size.setValue(defaults.blastp_word_size)
        self._blastp_gap_open.setValue(defaults.blastp_gap_open)
        self._blastp_gap_extend.setValue(defaults.blastp_gap_extend)
        self._blastp_query_coverage.setValue(defaults.blastp_query_coverage_hsp)
        self._blastp_write_text.setChecked(defaults.blastp_write_text)
        self._blastp_write_gff.setChecked(defaults.blastp_write_gff)
        self._blastx_evalue.setValue(defaults.blastx_e_value)
        self._blastx_threads.setValue(defaults.blastx_num_threads)
        self._blastx_word_size.setValue(defaults.blastx_word_size)
        self._blastx_gap_open.setValue(defaults.blastx_gap_open)
        self._blastx_gap_extend.setValue(defaults.blastx_gap_extend)
        self._blastx_query_coverage.setValue(defaults.blastx_query_coverage_hsp)
        self._blastx_write_text.setChecked(defaults.blastx_write_text)
        self._blastx_write_gff.setChecked(defaults.blastx_write_gff)
        self._max_msa.setValue(defaults.max_msa_lines)
        self._top_msa.setValue(defaults.top_msa_lines)
        self._min_bit_score.setValue(defaults.minimum_bit_score)
        self._blastn_database.setCurrentIndex(0)
        self._blastp_database.setCurrentIndex(0)
        self._blastx_database.setCurrentIndex(0)
        self._use_system_blast.setChecked(False)
        self._show_process_window.setChecked(defaults.show_process_window)
        self._right_click_buttons[BLAST_PROGRAM_BLASTN].setChecked(True)

    def use_system_blast(self) -> bool:
        return self._use_system_blast.isChecked()

    def options(self) -> BlastOptions:
        program = {
            0: BLAST_PROGRAM_BLASTN,
            1: BLAST_PROGRAM_BLASTP,
            2: BLAST_PROGRAM_BLASTX,
        }.get(self._tabs.currentIndex(), self._initial_program)
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
            write_gff=self._blastn_write_gff.isChecked(),
            minimum_bit_score=self._min_bit_score.value(),
            blast_database_name=self._blastn_database.currentData(),
            blastp_e_value=normalize_blast_evalue(self._blastp_evalue.value()),
            blastp_num_threads=self._blastp_threads.value(),
            blastp_word_size=self._blastp_word_size.value(),
            blastp_gap_open=self._blastp_gap_open.value(),
            blastp_gap_extend=self._blastp_gap_extend.value(),
            blastp_query_coverage_hsp=self._blastp_query_coverage.value(),
            blastp_database_name=self._blastp_database.currentData(),
            blastp_write_text=self._blastp_write_text.isChecked(),
            blastp_write_gff=self._blastp_write_gff.isChecked(),
            blastx_e_value=normalize_blast_evalue(self._blastx_evalue.value()),
            blastx_num_threads=self._blastx_threads.value(),
            blastx_word_size=self._blastx_word_size.value(),
            blastx_gap_open=self._blastx_gap_open.value(),
            blastx_gap_extend=self._blastx_gap_extend.value(),
            blastx_query_coverage_hsp=self._blastx_query_coverage.value(),
            blastx_database_name=self._blastx_database.currentData(),
            blastx_write_text=self._blastx_write_text.isChecked(),
            blastx_write_gff=self._blastx_write_gff.isChecked(),
            right_click_program=self._selected_right_click_program(),
            show_process_window=self._show_process_window.isChecked(),
        )

    def _selected_right_click_program(self) -> str:
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                return value
        return BLAST_PROGRAM_BLASTN


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

        self._threshold = QDoubleSpinBox()
        self._threshold.setRange(0.0, 1.0)
        self._threshold.setDecimals(2)
        self._threshold.setSingleStep(0.05)
        self._threshold.setValue(opts.threshold)
        self._threshold.setFixedWidth(260)
        self._threshold.setToolTip(
            "Minimum fraction of sequences that must share a base for it to "
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
        self._ambiguous.setFixedWidth(260)
        self._ambiguous.setToolTip(
            "Character used when no base reaches the majority threshold."
        )

        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("Run Consensus with saved parameters", CONSENSUS_RIGHT_CLICK_RUN),
            ("Open this dialog", CONSENSUS_RIGHT_CLICK_DIALOG),
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
        params_form.addRow("Majority threshold:", self._threshold)
        params_form.addRow("Ambiguous base:", self._ambiguous)
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
        page.ensurePolished()
        page.adjustSize()
        page_hint = page.sizeHint()
        tab_bar_h = max(24, self._tabs.tabBar().sizeHint().height())
        # Frame/padding around the page inside QTabWidget.
        tab_chrome = tab_bar_h + 8
        buttons_h = self._button_box.sizeHint().height()
        margins = self.layout().contentsMargins()
        spacing = self.layout().spacing()
        height = (
            page_hint.height()
            + tab_chrome
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
        self.setFixedSize(width, height)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _restore_defaults(self) -> None:
        self._threshold.setValue(DEFAULT_CONSENSUS_THRESHOLD)
        amb_index = self._ambiguous.findText(DEFAULT_CONSENSUS_AMBIGUOUS)
        self._ambiguous.setCurrentIndex(amb_index if amb_index >= 0 else 0)
        self._right_click_buttons[CONSENSUS_RIGHT_CLICK_RUN].setChecked(True)

    def options(self) -> ConsensusOptions:
        right_click_action = CONSENSUS_RIGHT_CLICK_RUN
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click_action = value
                break
        return ConsensusOptions(
            threshold=self._threshold.value(),
            ambiguous=normalize_consensus_ambiguous(self._ambiguous.currentText()),
            right_click_action=right_click_action,
        )


class FlexidotParamsDialog(QDialog):
    """Flexidot parameters.

    Left-click the toolbar button opens this dialog. The Parameters tab holds
    matching and output settings; the Options tab sets what right-clicking the
    toolbar button does.
    """

    def __init__(
        self,
        *,
        initial: FlexidotOptions | None = None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Flexidot parameters")
        opts = initial or FlexidotOptions()
        self._for_run = for_run
        self._run_after = False
        # These supported settings remain persisted but are intentionally not
        # exposed in this compact dialog.
        self._lcs_shading = opts.lcs_shading
        self._only_vs_first_value = opts.only_vs_first
        self._wordsize = QSpinBox()
        self._wordsize.setRange(MIN_FLEXIDOT_WORDSIZE, MAX_FLEXIDOT_WORDSIZE)
        self._wordsize.setValue(
            max(MIN_FLEXIDOT_WORDSIZE, min(MAX_FLEXIDOT_WORDSIZE, opts.wordsize))
        )
        self._wordsize.setToolTip(
            f"K-mer length (-k), {MIN_FLEXIDOT_WORDSIZE}–{MAX_FLEXIDOT_WORDSIZE}. "
            "Larger values are faster but less sensitive."
        )

        self._sequence_type = QComboBox()
        self._sequence_type.addItem("DNA / nucleotide (-t nuc)", FLEXIDOT_SEQUENCE_NUC)
        self._sequence_type.addItem("Protein (-t aa)", FLEXIDOT_SEQUENCE_AA)
        seq_type = normalize_flexidot_sequence_type(opts.sequence_type)
        seq_idx = self._sequence_type.findData(seq_type)
        self._sequence_type.setCurrentIndex(seq_idx if seq_idx >= 0 else 0)
        self._sequence_type.setToolTip(
            "Flexidot -t: sequence alphabet. DNA (nuc) is the default; "
            "choose Protein (aa) for amino-acid FASTA such as .faa files."
        )
        self._sequence_type.currentIndexChanged.connect(self._on_sequence_type_changed)

        self._use_gff = QCheckBox("Show GFF features on dotplot (-g)")
        self._use_gff.setChecked(opts.use_gff)
        self._use_gff.setToolTip(
            "Shade the dotplot with companion annotation GFF next to each FASTA "
            "(e.g. Name.RV.gff3 from the Annotation tab). "
            "BLAST hit tables (*_b.gff) are not used. "
            "If no companion GFF is found, Flexidot runs without shading."
        )

        self._output_format = QComboBox()
        # HTML (interactive Go-to) is the default and needs the forked Flexidot.
        # If the stock build is installed it automatically falls back to PNG at
        # run time, so it is always safe to offer/select here.
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
        self._output_format.setToolTip(
            f"Flexidot output format (-f). HTML gives interactive dotplots with "
            f"Go-to → Annotation inside {APP_DISPLAY_NAME} (forked Flexidot); it "
            "falls back to PNG automatically if that build isn't installed. "
            "PNG/SVG/PDF write a static image."
        )

        self._substitutions = QSpinBox()
        self._substitutions.setRange(0, 20)
        self._substitutions.setValue(opts.substitution_count)
        self._substitutions.setToolTip(
            "Flexidot -S: substitutions allowed inside each k-mer window. "
            "Zero uses the much faster exact-match engine."
        )
        self._wobble = QCheckBox("Expand ambiguous IUPAC bases (-w)")
        self._wobble.setChecked(opts.wobble_conversion)
        self._wobble.setToolTip(
            "Treat ambiguity codes such as R, Y, and N as their possible bases. "
            "This can substantially increase computation time."
        )
        self._max_n = QDoubleSpinBox()
        self._max_n.setRange(0.0, 100.0)
        self._max_n.setDecimals(1)
        self._max_n.setSuffix(" %")
        self._max_n.setValue(opts.max_n_percentage)
        self._max_n.setToolTip(
            "Flexidot --max_n: maximum N percentage allowed in a k-mer when "
            "ambiguous-base expansion is enabled."
        )
        self._wobble.toggled.connect(self._max_n.setEnabled)
        self._max_n.setEnabled(self._wobble.isChecked())

        self._include_reverse = QCheckBox(
            "Calculate reverse-complement matches (disable: -r / --norev)"
        )
        self._include_reverse.setChecked(opts.include_reverse)
        self._include_reverse.setToolTip(
            "Enabled by default. Clear it to pass Flexidot -r / --norev and "
            "omit reverse-complement matches."
        )
        self._line_width = QDoubleSpinBox()
        self._line_width.setRange(0.1, 20.0)
        self._line_width.setDecimals(1)
        self._line_width.setSingleStep(0.2)
        self._line_width.setValue(opts.line_width)
        self._line_width.setToolTip(
            "Flexidot -A / --line_width: thickness of forward and "
            "reverse-complement match lines in the generated dotplot."
        )

        for field in (
            self._wordsize,
            self._sequence_type,
            self._output_format,
            self._substitutions,
            self._max_n,
            self._line_width,
        ):
            field.setFixedWidth(260)

        self._length_scaling = QCheckBox(
            "Scale plots by sequence length (-L, pairwise / all-vs-all)"
        )
        self._length_scaling.setChecked(opts.length_scaling)
        self._length_scaling.setToolTip(
            "Scale each dotplot to the sequence lengths instead of a square "
            "panel. Only affects pairwise and all-vs-all comparisons."
        )
        self._combined = QCheckBox("Combine all selected files into one dotplot")
        self._combined.setChecked(opts.combined)
        self._combined.setToolTip(
            "Draw a single dotplot comparing all selected files. When unticked, "
            "each selected file is plotted separately."
        )

        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("Run Flexidot with saved parameters", FLEXIDOT_RIGHT_CLICK_RUN),
            ("Open this dialog", FLEXIDOT_RIGHT_CLICK_DIALOG),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[FLEXIDOT_RIGHT_CLICK_RUN].setChecked(True)

        rc_box = QGroupBox("Right-click action", self)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(
            QLabel("Right-clicking the Flexidot button does:", self)
        )
        for radio in self._right_click_buttons.values():
            rc_box_layout.addWidget(radio)

        process_box, self._show_process_window = make_process_output_box(
            self,
            tool_name="Flexidot",
            checked=bool(getattr(opts, "show_process_window", False)),
        )

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)
        options_layout.addWidget(rc_box)
        options_layout.addWidget(process_box)
        options_layout.addStretch(1)

        self._tabs = QTabWidget(self)
        params_tab = QWidget(self)
        params_form = QFormLayout(params_tab)
        params_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        params_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        params_form.setHorizontalSpacing(12)
        params_form.setVerticalSpacing(10)
        params_form.addRow("Word size (k-mer, -k):", self._wordsize)
        params_form.addRow("Sequence type (-t):", self._sequence_type)
        params_form.addRow("Output format (-f):", self._output_format)
        params_form.addRow("Allowed substitutions (-S):", self._substitutions)
        params_form.addRow("Maximum N content (--max_n):", self._max_n)
        params_form.addRow(
            "Match line thickness (-A / --line_width):", self._line_width
        )
        params_form.addRow("", self._wobble)
        params_form.addRow("", self._include_reverse)
        params_form.addRow("", self._length_scaling)
        params_form.addRow("", self._use_gff)
        params_form.addRow("", self._combined)
        self._tabs.addTab(params_tab, "Parameters")
        self._tabs.addTab(options_tab, "Options")

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore],
        )

        layout = QVBoxLayout(self)
        layout.addWidget(self._tabs)
        layout.addWidget(self._button_box)
        self._on_sequence_type_changed()

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _on_sequence_type_changed(self, *_args) -> None:
        is_protein = (
            str(self._sequence_type.currentData()) == FLEXIDOT_SEQUENCE_AA
        )
        # Reverse-complement matching only applies to DNA.
        self._include_reverse.setEnabled(not is_protein)
        if is_protein:
            self._include_reverse.setChecked(False)

    def _restore_defaults(self) -> None:
        defaults = FlexidotOptions()
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
        self._combined.setChecked(defaults.combined)
        self._show_process_window.setChecked(defaults.show_process_window)
        self._right_click_buttons[FLEXIDOT_RIGHT_CLICK_RUN].setChecked(True)
        self._on_sequence_type_changed()

    def options(self) -> FlexidotOptions:
        fmt = self._output_format.currentData()
        right_click_action = FLEXIDOT_RIGHT_CLICK_RUN
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click_action = value
                break
        seq_type = normalize_flexidot_sequence_type(self._sequence_type.currentData())
        include_reverse = self._include_reverse.isChecked()
        if seq_type == FLEXIDOT_SEQUENCE_AA:
            include_reverse = False
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
            combined=self._combined.isChecked(),
            # Comparison mode is automatic: self for one sequence,
            # pairwise/all-vs-all for multiple.
            modes=DEFAULT_FLEXIDOT_MODES,
            right_click_action=right_click_action,
            show_process_window=self._show_process_window.isChecked(),
        )


_TE_ANNOTATION_PARAMS_WIDTH = 580
_TE_ANNOTATION_OPTIONS_SIZE = (640, 560)


def _section_title_label(text: str, parent: QWidget | None = None) -> QLabel:
    label = QLabel(text, parent)
    font = label.font()
    font.setBold(True)
    label.setFont(font)
    return label


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
    """Batch TE annotation: Parameters select modules; Options holds settings/colors."""

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
        )

        super().__init__(parent)
        self.setWindowTitle("Sequence Annotation")
        opts = initial or TeAnnotationJobOptions()
        self._for_run = for_run
        self._run_after = False

        self._TeAnnotationOptions = TeAnnotationOptions
        self._SimilarityDatabaseRun = SimilarityDatabaseRun
        self._similarity_database_color = similarity_database_color
        self._AutoAnnotationMode = AutoAnnotationMode
        self._auto_all_updating = False
        self._similarity_all_updating = False

        selected_modes = set(opts.auto_modes)

        auto_section = QWidget(self)
        auto_outer = QVBoxLayout(auto_section)
        auto_outer.setContentsMargins(0, 0, 0, 0)
        auto_outer.setSpacing(6)
        auto_outer.addWidget(_section_title_label("Auto annotation", auto_section))

        auto_checks_row = QHBoxLayout()
        auto_checks_row.setSpacing(10)
        self._auto_ann_checks: dict = {}
        for label, mode in (
            ("LTR", AutoAnnotationMode.LTR),
            ("TIR", AutoAnnotationMode.TIR),
            ("Tandem", AutoAnnotationMode.TANDEM),
            ("ORF", AutoAnnotationMode.ORF),
        ):
            check = QCheckBox(label, auto_section)
            check.setChecked(mode.value in selected_modes)
            check.toggled.connect(self._on_auto_check_toggled)
            self._auto_ann_checks[mode] = check
            auto_checks_row.addWidget(check)

        self._auto_ann_all_check = QCheckBox("All", auto_section)
        self._auto_ann_all_check.setToolTip("Select or clear all auto-annotation detectors")
        self._auto_ann_all_check.toggled.connect(self._on_auto_all_toggled)
        auto_checks_row.addWidget(self._auto_ann_all_check)
        auto_checks_row.addStretch(1)
        auto_outer.addLayout(auto_checks_row)
        self._sync_auto_all_check()

        sim_section = QWidget(self)
        sim_outer = QVBoxLayout(sim_section)
        sim_outer.setContentsMargins(0, 0, 0, 0)
        sim_outer.setSpacing(6)
        sim_outer.addWidget(_section_title_label("Similarity annotation", sim_section))

        self._similarity_ann_all_check = QCheckBox("All", sim_section)
        self._similarity_ann_all_check.setToolTip("Select or clear all similarity databases")
        self._similarity_ann_all_check.toggled.connect(self._on_similarity_all_toggled)

        self._similarity_db_rows: list[tuple] = []
        self._db_flow = _CheckboxFlowWidget(sim_section)

        databases = load_similarity_annotation_databases()
        select_all_dbs = opts.similarity_database_names is None
        selected_dbs = set(opts.similarity_database_names or ())
        if not databases:
            empty = QLabel(
                "No databases defined. Add them in File → Define database.",
                sim_section,
            )
            empty.setWordWrap(True)
            sim_outer.addWidget(empty)
            self._similarity_ann_all_check.hide()
            self._db_flow.hide()
        else:
            flow_items: list[QWidget] = []
            for index, entry in enumerate(databases):
                check = self._make_similarity_db_check(
                    entry,
                    index,
                    self._db_flow,
                    checked=select_all_dbs or entry.name in selected_dbs,
                )
                self._similarity_db_rows.append((entry, check))
                flow_items.append(check)
            flow_items.append(self._similarity_ann_all_check)
            self._db_flow.set_items(flow_items)
            sim_outer.addWidget(self._db_flow)
        self._sync_similarity_all_check()

        params_tab = QWidget(self)
        params_layout = QVBoxLayout(params_tab)
        params_layout.setContentsMargins(12, 10, 12, 10)
        params_layout.setSpacing(12)
        params_layout.addWidget(auto_section)
        params_layout.addWidget(sim_section)

        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("Run Sequence annotation with saved selection", TE_ANNOTATION_RIGHT_CLICK_RUN),
            ("Open this dialog", TE_ANNOTATION_RIGHT_CLICK_DIALOG),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[TE_ANNOTATION_RIGHT_CLICK_RUN].setChecked(True)

        rc_box = QGroupBox("Right-click action", self)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(
            QLabel("Right-clicking the Sequence Annotation button does:", self)
        )
        for radio in self._right_click_buttons.values():
            rc_box_layout.addWidget(radio)

        self._auto_settings_panel = AutoAnnotationSettingsPanel(
            load_auto_annotation_settings(), self
        )
        self._similarity_settings_panel = SimilarityAnnotationSettingsPanel(
            load_similarity_annotation_settings(), self
        )
        self._colors_panel = AnnotationElementColorsPanel(self)

        auto_box = QGroupBox("Auto annotation settings", self)
        auto_box_layout = QVBoxLayout(auto_box)
        auto_box_layout.addWidget(self._auto_settings_panel)

        sim_box = QGroupBox("Similarity annotation settings", self)
        sim_box_layout = QVBoxLayout(sim_box)
        sim_box_layout.addWidget(self._similarity_settings_panel)

        color_box = QGroupBox("Element colors", self)
        color_box_layout = QVBoxLayout(color_box)
        color_box_layout.addWidget(self._colors_panel)

        options_body = QWidget(self)
        options_layout = QVBoxLayout(options_body)
        options_layout.setContentsMargins(12, 10, 12, 10)
        options_layout.setSpacing(12)
        options_layout.addWidget(rc_box)
        options_layout.addWidget(auto_box)
        options_layout.addWidget(sim_box)
        options_layout.addWidget(color_box)
        options_layout.addStretch(1)

        options_scroll = QScrollArea(self)
        options_scroll.setWidgetResizable(True)
        options_scroll.setFrameShape(QFrame.Shape.NoFrame)
        options_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        options_scroll.setWidget(options_body)

        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        self._tabs.addTab(params_tab, "Parameters")
        self._tabs.addTab(options_scroll, "Options")
        self._params_tab = params_tab

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
        QTimer.singleShot(0, self._fit_to_current_tab)

    def _on_tab_changed(self, tab_index: int) -> None:
        self._fit_to_current_tab()

    def _fit_to_current_tab(self) -> None:
        """Size the dialog to the active tab (Parameters compact; Options scrollable)."""
        margins = self.layout().contentsMargins()
        spacing = self.layout().spacing()
        tab_bar_h = max(24, self._tabs.tabBar().sizeHint().height())
        tab_chrome = tab_bar_h + 8
        buttons_h = self._button_box.sizeHint().height()
        chrome = (
            tab_chrome
            + buttons_h
            + margins.top()
            + margins.bottom()
            + spacing
            + 4
        )
        if self._tabs.currentIndex() == 0:
            page = self._params_tab
            page.ensurePolished()
            if hasattr(self, "_db_flow") and self._db_flow.isVisible():
                self._db_flow._last_flow_width = -1
                self._db_flow.setFixedWidth(
                    max(200, _TE_ANNOTATION_PARAMS_WIDTH - 40)
                )
                self._db_flow._reflow()
            page.adjustSize()
            hint = page.sizeHint()
            width = max(
                _TE_ANNOTATION_PARAMS_WIDTH,
                hint.width() + margins.left() + margins.right() + 24,
                self._button_box.sizeHint().width() + margins.left() + margins.right(),
            )
            height = hint.height() + chrome
            self.setFixedSize(width, height)
            return
        opt_w, opt_h = _TE_ANNOTATION_OPTIONS_SIZE
        self.setFixedSize(opt_w, opt_h)

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

    def _on_auto_check_toggled(self, _checked: bool) -> None:
        if self._auto_all_updating:
            return
        self._sync_auto_all_check()

    def _sync_auto_all_check(self) -> None:
        checks = list(self._auto_ann_checks.values())
        if not checks:
            return
        all_checked = all(check.isChecked() for check in checks)
        self._auto_all_updating = True
        self._auto_ann_all_check.setChecked(all_checked)
        self._auto_all_updating = False

    def _on_auto_all_toggled(self, checked: bool) -> None:
        if self._auto_all_updating:
            return
        self._auto_all_updating = True
        for check in self._auto_ann_checks.values():
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
        from repeatviewer.widgets.sequence_map_viewer import (
            AutoAnnotationSettingsPanel,
            SimilarityAnnotationSettingsPanel,
        )

        new_auto = AutoAnnotationSettingsPanel(AutoAnnotationSettings(), parent_auto)
        new_sim = SimilarityAnnotationSettingsPanel(
            SimilarityAnnotationSettings(), parent_sim
        )
        if parent_auto is not None:
            layout = parent_auto.layout()
            if layout is not None:
                layout.replaceWidget(self._auto_settings_panel, new_auto)
        self._auto_settings_panel.deleteLater()
        self._auto_settings_panel = new_auto
        if parent_sim is not None:
            layout = parent_sim.layout()
            if layout is not None:
                layout.replaceWidget(self._similarity_settings_panel, new_sim)
        self._similarity_settings_panel.deleteLater()
        self._similarity_settings_panel = new_sim
        self._colors_panel.restore_defaults()
        self._right_click_buttons[TE_ANNOTATION_RIGHT_CLICK_RUN].setChecked(True)

    def apply_settings_panels(self) -> None:
        """Persist Options-tab auto/similarity settings and colors."""
        from repeatviewer.backends.auto_annotation_settings import (
            save_auto_annotation_settings,
        )
        from repeatviewer.backends.similarity_annotation_settings import (
            save_similarity_annotation_settings,
        )

        save_auto_annotation_settings(self._auto_settings_panel.settings())
        save_similarity_annotation_settings(self._similarity_settings_panel.settings())
        self._colors_panel.save()
        self._refresh_similarity_db_checkbox_colors()

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
        from repeatviewer.backends.auto_annotation_settings import OrfSettings as _OrfSettings

        super().__init__(parent)
        self.setWindowTitle("ORF Prediction settings")
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

        def _bp_row(spin: QSpinBox) -> QWidget:
            row = QWidget(self)
            layout = QHBoxLayout(row)
            layout.setContentsMargins(0, 0, 0, 0)
            layout.setSpacing(6)
            layout.addWidget(spin)
            layout.addWidget(QLabel("bp", row))
            layout.addStretch()
            return row

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        form.addRow("Min ORF length:", _bp_row(self._min_len))
        form.addRow("Max ORF length:", _bp_row(self._max_len))
        form.addRow("Add other start codons:", self._alt_starts)

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
        defaults = self._OrfSettings()
        self._min_len.setValue(defaults.min_length)
        self._max_len.setValue(defaults.max_length)
        self._alt_starts.setText(defaults.alternative_start_codons)

    def options(self):
        return self._OrfSettings(
            min_length=self._min_len.value(),
            max_length=max(self._max_len.value(), self._min_len.value()),
            alternative_start_codons=self._alt_starts.text().strip(),
        )


class TeaAidParamsDialog(QDialog):
    """TE-Aid parameters.

    Left-click the toolbar button opens this dialog. The Parameters tab holds
    BLAST settings and hit colors; the Options tab sets right-click behavior.
    """

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
        opts = initial or TeaAidOptions()
        self._for_run = for_run
        self._run_after = False
        # genome_history kept for call-site compatibility; Genome list uses BLAST DBs.
        _ = genome_history

        high_pct, mid_pct = clamp_teaid_hit_len_pct(
            opts.hit_len_high_pct, opts.hit_len_mid_pct
        )

        # Genome: same chrome as Blast / Extend — current genome + nucleotide DBs.
        self._genome_row = GenomeBlastDatabaseRow(
            self,
            project=project,
            project_genome=project_genome,
            selected_path=opts.genome_path,
            manage_title="TE-Aid parameters",
        )

        self._evalue = QDoubleSpinBox()
        configure_report_evalue_spinbox(self._evalue)
        self._evalue.setValue(
            normalize_blast_evalue(opts.e_value, default=DEFAULT_TEAID_EVALUE)
        )
        self._evalue.setFixedWidth(120)
        self._evalue.setToolTip(
            "Reporting E-value for the TE-Aid genome BLAST search only."
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
        self._word_size.setFixedWidth(120)
        self._word_size.setToolTip(
            "BLASTN word size (-word_size) for the TE-Aid self-BLAST dot plot only (4–25)."
        )

        self._line_transparency = QDoubleSpinBox()
        self._line_transparency.setRange(0.0, 100.0)
        self._line_transparency.setDecimals(0)
        self._line_transparency.setSingleStep(5.0)
        self._line_transparency.setSuffix(" %")
        self._line_transparency.setFixedWidth(120)
        self._line_transparency.setValue(
            clamp_teaid_hit_line_transparency(opts.hit_line_transparency) * 100.0
        )
        self._line_transparency.setToolTip(
            "Transparency of all genome BLAST hit lines (black, mid, and full-length). "
            "0% = fully opaque; 100% = invisible."
        )

        self._hit_high_pct = QDoubleSpinBox()
        self._hit_high_pct.setRange(0.0, 100.0)
        self._hit_high_pct.setDecimals(1)
        self._hit_high_pct.setSingleStep(1.0)
        self._hit_high_pct.setSuffix(" %")
        self._hit_high_pct.setValue(high_pct)
        self._hit_high_pct.setFixedWidth(100)
        self._hit_high_pct.setToolTip(
            "Hits at least this long (percent of the input / consensus length) "
            "use the upper color."
        )
        self._hit_color_high = opts.hit_color_high or DEFAULT_TEAID_HIT_COLOR_HIGH
        self._hit_color_high_btn = QPushButton()
        self._hit_color_high_btn.setFixedSize(56, 22)
        self._hit_color_high_btn.setToolTip("Color for hits longer than the value above")
        self._hit_color_high_btn.clicked.connect(self._pick_hit_color_high)
        self._update_hit_color_button(self._hit_color_high_btn, self._hit_color_high)

        self._hit_mid_pct = QDoubleSpinBox()
        self._hit_mid_pct.setRange(0.0, 100.0)
        self._hit_mid_pct.setDecimals(1)
        self._hit_mid_pct.setSingleStep(1.0)
        self._hit_mid_pct.setSuffix(" %")
        self._hit_mid_pct.setValue(mid_pct)
        self._hit_mid_pct.setFixedWidth(100)
        self._hit_mid_pct.setToolTip(
            "Hits at least this long but shorter than the upper threshold use the "
            "mid color. Cannot be greater than or equal to the upper threshold."
        )
        self._hit_color_mid = opts.hit_color_mid or DEFAULT_TEAID_HIT_COLOR_MID
        self._hit_color_mid_btn = QPushButton()
        self._hit_color_mid_btn.setFixedSize(56, 22)
        self._hit_color_mid_btn.setToolTip(
            "Color for hits longer than this value but below the upper threshold"
        )
        self._hit_color_mid_btn.clicked.connect(self._pick_hit_color_mid)
        self._update_hit_color_button(self._hit_color_mid_btn, self._hit_color_mid)

        self._hit_high_pct.valueChanged.connect(self._on_hit_high_changed)
        self._hit_mid_pct.valueChanged.connect(self._on_hit_mid_changed)

        high_row = QWidget(self)
        high_layout = QHBoxLayout(high_row)
        high_layout.setContentsMargins(0, 0, 0, 0)
        high_layout.setSpacing(8)
        high_layout.addWidget(self._hit_color_high_btn)
        high_layout.addWidget(QLabel("Hit longer than:", self))
        high_layout.addWidget(self._hit_high_pct)
        high_layout.addWidget(QLabel("of sequence length", self))
        high_layout.addStretch(1)

        mid_row = QWidget(self)
        mid_layout = QHBoxLayout(mid_row)
        mid_layout.setContentsMargins(0, 0, 0, 0)
        mid_layout.setSpacing(8)
        mid_layout.addWidget(self._hit_color_mid_btn)
        mid_layout.addWidget(QLabel("Hit longer than:", self))
        mid_layout.addWidget(self._hit_mid_pct)
        mid_layout.addWidget(QLabel("but below the value above", self))
        mid_layout.addStretch(1)

        colors_box = QGroupBox("Genome BLAST hit colors", self)
        colors_layout = QVBoxLayout(colors_box)
        colors_layout.setContentsMargins(10, 8, 10, 8)
        colors_layout.setSpacing(8)
        colors_layout.addWidget(high_row)
        colors_layout.addWidget(mid_row)

        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("Run TE-Aid with saved parameters", TEAID_RIGHT_CLICK_RUN),
            ("Open this dialog", TEAID_RIGHT_CLICK_DIALOG),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[TEAID_RIGHT_CLICK_RUN].setChecked(True)

        rc_box = QGroupBox("Right-click action", self)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(QLabel("Right-clicking the TE-Aid button does:", self))
        for radio in self._right_click_buttons.values():
            rc_box_layout.addWidget(radio)

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(12, 10, 12, 10)
        options_layout.setSpacing(8)
        options_layout.addWidget(rc_box)

        self._x_axis_label = QLineEdit(self)
        self._x_axis_label.setText(
            (opts.x_axis_label or "").strip() or DEFAULT_TEAID_X_AXIS_LABEL
        )
        self._x_axis_label.setPlaceholderText(DEFAULT_TEAID_X_AXIS_LABEL)
        self._x_axis_label.setToolTip(
            "Axis title for the consensus position in the TE-Aid plots."
        )
        axis_box = QGroupBox("Plot labels", self)
        axis_form = QFormLayout(axis_box)
        axis_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        axis_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        axis_form.addRow(
            "Consensus axis label:",
            field_with_help(
                self._x_axis_label,
                "Consensus axis label",
                "Title of the consensus-position axis in the divergence, "
                "coverage and self dot plots (default: "
                f"\"{DEFAULT_TEAID_X_AXIS_LABEL}\"). Change it when the input "
                "is not a TE, e.g. \"satDNA monomer (bp)\". "
                "Leave empty to use the default.",
                parent=self,
            ),
        )
        options_layout.addWidget(axis_box)

        self._tabs = QTabWidget(self)
        self._tabs.setDocumentMode(True)
        params_tab = QWidget(self)
        params_layout = QVBoxLayout(params_tab)
        params_layout.setContentsMargins(12, 10, 12, 10)
        params_layout.setSpacing(10)
        params_form = QFormLayout()
        params_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        params_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        params_form.setHorizontalSpacing(12)
        params_form.setVerticalSpacing(8)
        params_form.addRow("Genome:", self._genome_row)
        params_form.addRow("Genome BLAST E-value:", self._evalue)
        params_form.addRow("Self BLAST word size:", self._word_size)
        params_form.addRow("Line transparency:", self._line_transparency)
        params_layout.addLayout(params_form)
        params_layout.addWidget(colors_box)
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
        layout.addWidget(self._tabs)
        layout.addWidget(self._button_box)
        self._tabs.currentChanged.connect(self._on_tab_changed)
        self._on_tab_changed(self._tabs.currentIndex())

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        QTimer.singleShot(0, self._fit_to_current_tab)

    def _on_tab_changed(self, tab_index: int) -> None:
        self._fit_to_current_tab()

    def _fit_to_current_tab(self) -> None:
        page = self._tabs.currentWidget()
        if page is None:
            return
        page.ensurePolished()
        page.adjustSize()
        page_hint = page.sizeHint()
        tab_bar_h = max(24, self._tabs.tabBar().sizeHint().height())
        tab_chrome = tab_bar_h + 8
        buttons_h = self._button_box.sizeHint().height()
        margins = self.layout().contentsMargins()
        spacing = self.layout().spacing()
        height = (
            page_hint.height()
            + tab_chrome
            + buttons_h
            + margins.top()
            + margins.bottom()
            + spacing
            + 4
        )
        width = max(
            600,
            page_hint.width() + margins.left() + margins.right() + 40,
            self._button_box.sizeHint().width() + margins.left() + margins.right(),
        )
        self.setFixedSize(width, height)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    @staticmethod
    def _update_hit_color_button(button: QPushButton, color: str) -> None:
        hex_color = QColor(color).name(QColor.NameFormat.HexRgb)
        button.setStyleSheet(f"background-color: {hex_color};")

    def _pick_hit_color_high(self) -> None:
        chosen = QColorDialog.getColor(
            QColor(self._hit_color_high), self, "Hit color (upper band)"
        )
        if not chosen.isValid():
            return
        self._hit_color_high = chosen.name(QColor.NameFormat.HexRgb)
        self._update_hit_color_button(self._hit_color_high_btn, self._hit_color_high)

    def _pick_hit_color_mid(self) -> None:
        chosen = QColorDialog.getColor(
            QColor(self._hit_color_mid), self, "Hit color (mid band)"
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
        self._genome_row.reset_to_default()
        self._evalue.setValue(DEFAULT_TEAID_EVALUE)
        self._word_size.setValue(DEFAULT_TEAID_SELF_BLAST_WORD_SIZE)
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
        self._right_click_buttons[TEAID_RIGHT_CLICK_RUN].setChecked(True)
        self._x_axis_label.setText(DEFAULT_TEAID_X_AXIS_LABEL)

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
        )


class MsaCleanParamsDialog(QDialog):
    """MSA clean parameters.

    Opened from the MSA Clean drop menu / Tools menu. Parameters holds the
    Combined pipeline and step thresholds; Options sets right-click behavior
    when wired on the toolbar.
    """

    def __init__(
        self,
        *,
        initial: MsaCleanOptions | None = None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("MSA clean parameters")
        opts = initial or MsaCleanOptions()
        self._for_run = for_run
        self._run_after = False

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

        self._method_order = QListWidget()
        self._method_order.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self._method_order.setDefaultDropAction(Qt.DropAction.MoveAction)
        self._method_order.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        self._method_order.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
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
        # Ensure every known method appears (in case order was incomplete).
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

        combined_help = (
            "Combined MSA clean runs the checked cleaning steps in list order "
            "on the selected alignment.\n\n"
            "How to use it:\n"
            "• Check a step to include it in the pipeline.\n"
            "• Uncheck a step to skip it.\n"
            "• Drag steps up or down to change the order they run.\n"
            "• At least one step must stay checked.\n\n"
            "Each step uses the parameters in its section below. "
            "Use the ? next to a step title for details on that method."
        )
        combined_box = QGroupBox(self)
        combined_layout = QVBoxLayout(combined_box)
        combined_layout.setContentsMargins(10, 8, 10, 8)
        combined_layout.setSpacing(6)
        combined_header = QWidget(combined_box)
        combined_header_layout = QHBoxLayout(combined_header)
        combined_header_layout.setContentsMargins(0, 0, 0, 0)
        combined_header_layout.setSpacing(6)
        combined_title = QLabel("Combined MSA clean", combined_header)
        title_font = combined_title.font()
        title_font.setBold(True)
        combined_title.setFont(title_font)
        combined_header_layout.addWidget(combined_title, 1)
        combined_header_layout.addWidget(
            self._msa_help_button("Combined MSA clean", combined_help),
            0,
            Qt.AlignmentFlag.AlignRight,
        )
        combined_layout.addWidget(combined_header)
        hint = QLabel("Check steps to include. Drag to set order.", combined_box)
        hint.setWordWrap(True)
        combined_layout.addWidget(hint)
        combined_layout.addWidget(self._method_order)

        spin_width = 90

        div_box = QGroupBox(self)
        self._div_threshold = QDoubleSpinBox()
        self._div_threshold.setRange(0.0, 1.0)
        self._div_threshold.setDecimals(2)
        self._div_threshold.setSingleStep(0.05)
        self._div_threshold.setValue(opts.divergence.threshold)
        self._div_threshold.setFixedWidth(spin_width)

        self._div_window = QSpinBox()
        self._div_window.setRange(1, 10_000)
        self._div_window.setValue(opts.divergence.window_size)
        self._div_window.setFixedWidth(spin_width)

        self._fill_labeled_grid(
            div_box,
            [
                ("Divergence threshold:", self._div_threshold),
                ("Window size (bp):", self._div_window),
            ],
            heading="Clean rows with divergence",
            help_text=msa_clean_method_description("msa_clean_div"),
        )

        gap_box = QGroupBox(self)
        self._row_gap_threshold = QDoubleSpinBox()
        self._row_gap_threshold.setRange(0.0, 1.0)
        self._row_gap_threshold.setDecimals(3)
        self._row_gap_threshold.setSingleStep(0.01)
        self._row_gap_threshold.setValue(opts.gap_row.gap_threshold)
        self._row_gap_threshold.setFixedWidth(spin_width)

        self._row_gap_window = QSpinBox()
        self._row_gap_window.setRange(1, 10_000)
        self._row_gap_window.setValue(opts.gap_row.window_size)
        self._row_gap_window.setFixedWidth(spin_width)

        self._fill_labeled_grid(
            gap_box,
            [
                ("Gap proportion threshold:", self._row_gap_threshold),
                ("Window size (bp):", self._row_gap_window),
            ],
            heading="Clean rows with gaps",
            help_text=msa_clean_method_description("msa_clean_gap"),
        )

        gappy_box = QGroupBox(self)
        self._gap_threshold = QDoubleSpinBox()
        self._gap_threshold.setRange(0.0, 1.0)
        self._gap_threshold.setDecimals(2)
        self._gap_threshold.setSingleStep(0.05)
        self._gap_threshold.setValue(opts.gappy.gap_threshold)
        self._gap_threshold.setFixedWidth(spin_width)

        self._simi_gap = QDoubleSpinBox()
        self._simi_gap.setRange(0.0, 1.0)
        self._simi_gap.setDecimals(2)
        self._simi_gap.setSingleStep(0.05)
        self._simi_gap.setValue(opts.gappy.simi_check_gap_thr)
        self._simi_gap.setFixedWidth(spin_width)

        self._similarity = QDoubleSpinBox()
        self._similarity.setRange(0.0, 1.0)
        self._similarity.setDecimals(2)
        self._similarity.setSingleStep(0.05)
        self._similarity.setValue(opts.gappy.similarity_thr)
        self._similarity.setFixedWidth(spin_width)

        self._min_nt = QSpinBox()
        self._min_nt.setRange(1, 50)
        self._min_nt.setValue(opts.gappy.min_nucleotide)
        self._min_nt.setFixedWidth(spin_width)

        self._fill_labeled_grid(
            gappy_box,
            [
                ("Gap threshold:", self._gap_threshold),
                ("Similarity check gap threshold:", self._simi_gap),
                ("Similarity threshold:", self._similarity),
                ("Min nucleotides per column:", self._min_nt),
            ],
            heading="Remove gappy columns",
            help_text=msa_clean_method_description("msa_clean_gappy"),
        )

        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for label, value in (
            ("Run Combined MSA clean", MSA_CLEAN_RIGHT_CLICK_COMBINED),
            ("Run Remove gappy columns", MSA_CLEAN_RIGHT_CLICK_GAPPY),
            ("Run Clean rows with divergence", MSA_CLEAN_RIGHT_CLICK_DIV),
            ("Run Clean rows with gaps", MSA_CLEAN_RIGHT_CLICK_GAP),
            ("Open this dialog", MSA_CLEAN_RIGHT_CLICK_DIALOG),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[MSA_CLEAN_RIGHT_CLICK_COMBINED].setChecked(True)

        rc_box = QGroupBox("Right-click action", self)
        rc_box_layout = QVBoxLayout(rc_box)
        rc_box_layout.setSpacing(6)
        rc_box_layout.addWidget(
            QLabel("Right-clicking the MSA Clean button does:", self)
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
        params_tab = QWidget(self)
        params_layout = QVBoxLayout(params_tab)
        params_layout.setContentsMargins(12, 10, 12, 10)
        params_layout.setSpacing(10)
        params_layout.addWidget(combined_box)
        params_layout.addWidget(div_box)
        params_layout.addWidget(gap_box)
        params_layout.addWidget(gappy_box)
        self._tabs.addTab(params_tab, "Parameters")
        self._tabs.addTab(options_tab, "Options")

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore],
            on_accept=self._on_accept,
            include_run=False,
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
        self._fit_method_list_height()
        QTimer.singleShot(0, self._fit_to_current_tab)

    def _on_tab_changed(self, tab_index: int) -> None:
        self._fit_to_current_tab()

    def _fit_to_current_tab(self) -> None:
        page = self._tabs.currentWidget()
        if page is None:
            return
        page.ensurePolished()
        page.adjustSize()
        page_hint = page.sizeHint()
        tab_bar_h = max(24, self._tabs.tabBar().sizeHint().height())
        tab_chrome = tab_bar_h + 8
        buttons_h = self._button_box.sizeHint().height()
        margins = self.layout().contentsMargins()
        spacing = self.layout().spacing()
        height = (
            page_hint.height()
            + tab_chrome
            + buttons_h
            + margins.top()
            + margins.bottom()
            + spacing
            + 4
        )
        width = max(
            480,
            page_hint.width() + margins.left() + margins.right() + 24,
            self._button_box.sizeHint().width() + margins.left() + margins.right(),
        )
        self.setFixedSize(width, height)

    def _on_accept(self) -> None:
        if not self._checked_method_ids():
            QMessageBox.warning(
                self,
                "MSA clean",
                "Enable at least one step in Combined MSA clean.",
            )
            self._tabs.setCurrentIndex(0)
            return
        self.accept()

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _msa_help_button(self, title: str, help_text: str) -> QToolButton:
        return make_help_button(self, title, help_text)

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

    def _fill_labeled_grid(
        self,
        box,
        rows,
        *,
        heading: str | None = None,
        help_text: str = "",
    ) -> None:
        """Lay out (label, widget) rows in a group box with a real label column.

        Uses an explicit QGridLayout rather than QFormLayout: the form's label
        column was collapsing to nothing here, clipping the descriptions to a
        stray character. A grid gives column 0 its natural width and stretches
        column 1, so every label is always fully visible.
        """
        outer = QVBoxLayout(box)
        outer.setContentsMargins(10, 8, 10, 8)
        outer.setSpacing(8)
        if heading:
            header = QWidget(box)
            header_layout = QHBoxLayout(header)
            header_layout.setContentsMargins(0, 0, 0, 0)
            header_layout.setSpacing(6)
            title = QLabel(heading, header)
            font = title.font()
            font.setBold(True)
            title.setFont(font)
            header_layout.addWidget(title, 1)
            if help_text:
                header_layout.addWidget(
                    self._msa_help_button(heading, help_text),
                    0,
                    Qt.AlignmentFlag.AlignRight,
                )
            outer.addWidget(header)

        grid = QGridLayout()
        grid.setHorizontalSpacing(12)
        grid.setVerticalSpacing(8)
        grid.setColumnStretch(0, 0)
        grid.setColumnStretch(1, 1)
        for row, (text, widget) in enumerate(rows):
            label = QLabel(text)
            grid.addWidget(
                label,
                row,
                0,
                Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
            )
            grid.addWidget(widget, row, 1)
        outer.addLayout(grid)

    def _fit_method_list_height(self) -> None:
        """Keep the Combined list only as tall as its rows, with a little padding."""
        count = max(1, self._method_order.count())
        row_h = self._method_order.sizeHintForRow(0)
        if row_h <= 0:
            row_h = 22
        frame = 2 * self._method_order.frameWidth()
        # Extra space above the first row and below the last so the list
        # does not feel flush against the frame (matches stylesheet padding).
        padding = 16
        self._method_order.setFixedHeight(row_h * count + frame + padding)

    def _on_method_item_changed(self, item: QListWidgetItem) -> None:
        # Keep at least one step checked.
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

        self._rebuild_check_list()
        if self._preselected is None and self._chrom_sizes:
            self._select_min_size_mb(2.0)

    def _select_by_min_size_spinbox(self) -> None:
        self._select_min_size_mb(self._min_size_mb.value())

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
            selected = {self._chrom_name_from_checkbox(cb) for cb in self._check_boxes if cb.isChecked()}
            if not selected and self._preselected is not None:
                selected = set(self._preselected)
        while self._checks_layout.count():
            item = self._checks_layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()
        self._check_boxes.clear()

        ordered = self._sorted_names(self._sort_mode)
        if not ordered:
            self._checks_layout.addWidget(QLabel("No items available", self._checks_container))
            self._checks_layout.addStretch(1)
            return

        for name in ordered:
            cb = QCheckBox(self._display_label(name), self._checks_container)
            cb.setProperty("chrom_name", name)
            cb.setChecked(name in selected)
            self._check_boxes.append(cb)
            self._checks_layout.addWidget(cb)
        self._checks_layout.addStretch(1)
        self._apply_filter(self._filter.text())

    def set_chromosome_sort(self, mode: str) -> None:
        if mode not in (
            CHROMOSOME_SORT_NATURAL,
            CHROMOSOME_SORT_SIZE_DESC,
            CHROMOSOME_SORT_SIZE_ASC,
        ):
            mode = CHROMOSOME_SORT_NATURAL
        self._sort_mode = mode
        self._rebuild_check_list()

    def _apply_filter(self, text: str) -> None:
        needle = text.strip().lower()
        for cb in self._check_boxes:
            cb.setVisible(not needle or needle in cb.text().lower())

    def selected_names(self) -> list[str]:
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
        for cb in self._check_boxes:
            if cb.isVisible():
                cb.setChecked(True)

    def _uncheck_all(self) -> None:
        for cb in self._check_boxes:
            if cb.isVisible():
                cb.setChecked(False)

    def _select_min_size_mb(self, min_mb: float) -> None:
        min_bp = int(float(min_mb) * 1_000_000)
        selected = {
            name
            for name in self._all_chrom_names
            if self._chrom_sizes.get(name, 0) >= min_bp
        }
        self._rebuild_check_list(selected if selected else None)

    def name_count(self) -> int:
        return len(self._all_chrom_names)


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
    ) -> None:
        super().__init__(parent)
        self._row = PathBrowseRow(
            self,
            initial=initial,
            placeholder=placeholder,
            browse_title=browse_title,
            file_filter=file_filter,
            start_dir=start_dir,
            mode="file",
        )
        self._path = self._row.line_edit()
        self._row.activated.connect(lambda _w: self.activated.emit(self))
        self._row.path_changed.connect(self.path_changed.emit)

        remove_btn = QPushButton("Remove", self)
        remove_btn.setVisible(removable)
        remove_btn.clicked.connect(self.remove_requested.emit)
        self._remove_btn = remove_btn

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        layout.addWidget(self._row, stretch=1)
        layout.addWidget(remove_btn)

    def set_removable(self, removable: bool) -> None:
        self._remove_btn.setVisible(removable)

    def set_path(self, path: Path) -> None:
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
        add_button_label: str = "Add another TE file",
        add_button_tooltip: str = (
            "Add another RepeatMasker .out or BLAST hit table. "
            "Track order in the plot matches file order: first row = top track."
        ),
        allow_multiple: bool = True,
        start_dir: str = "",
    ) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
        self._placeholder = placeholder
        self._browse_title = browse_title
        self._file_filter = file_filter
        self._allow_multiple = allow_multiple
        self._start_dir = start_dir

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
        row = _TeInputFileRow(
            self,
            initial=path,
            removable=len(self._rows) >= 1,
            placeholder=self._placeholder,
            browse_title=self._browse_title,
            file_filter=self._file_filter,
            start_dir=self._start_dir,
        )
        row.remove_requested.connect(lambda r=row: self._remove_row(r))
        row.activated.connect(self._on_row_activated)
        row.path_changed.connect(self._emit_paths_changed)
        self._rows.append(row)
        self._rows_layout.addWidget(row)
        self._refresh_remove_buttons()
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
    """Optional GFF/GTF/BED/VCF annotation path entry."""

    def __init__(self, parent: QWidget | None = None, *, initial: Path | None = None) -> None:
        super().__init__(parent)
        self._row = PathBrowseRow(
            self,
            initial=initial,
            placeholder="Optional GFF/GTF/BED/VCF annotation",
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
    """Thin horizontal rule between major dialog sections."""
    host = QWidget()
    host.setFixedHeight(9)
    layout = QVBoxLayout(host)
    layout.setContentsMargins(0, 2, 0, 2)
    line = QFrame()
    line.setFrameShape(QFrame.Shape.HLine)
    line.setFrameShadow(QFrame.Shadow.Plain)
    layout.addWidget(line)
    return host


class AnnotationFeatureTypesWidget(QWidget):
    """Checklist of GFF/GTF feature types (or BED intervals)."""

    _DEFAULT_TYPES = ("gene", "mRNA", "CDS", "exon", "transcript")

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._list = QListWidget(self)
        self._list.setMaximumHeight(92)
        self._list.setMinimumHeight(72)
        self._list.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        self._list.setToolTip(
            "Feature types from column 3 of GFF/GTF files. BED files use a single "
            "'interval' type. VCF files use SNP, INDEL, MNP, SV, and other variant classes."
        )
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._list)
        self._populate_default()

    def _populate_default(self) -> None:
        checked = {t.lower() for t in DEFAULT_ANNOTATION_FEATURE_TYPES}
        self._set_types(list(self._DEFAULT_TYPES), checked=checked)

    def _set_types(self, types: list[str], *, checked: set[str] | None = None) -> None:
        prev = checked if checked is not None else {t.lower() for t in self.selected_types()}
        self._list.clear()
        for feature_type in types:
            item = QListWidgetItem(feature_type)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            state = (
                Qt.CheckState.Checked
                if feature_type.lower() in prev
                else Qt.CheckState.Unchecked
            )
            item.setCheckState(state)
            self._list.addItem(item)

    def selected_types(self) -> list[str]:
        out: list[str] = []
        for index in range(self._list.count()):
            item = self._list.item(index)
            if item is not None and item.checkState() == Qt.CheckState.Checked:
                out.append(item.text())
        return out

    def set_checked_types(self, types: list[str]) -> None:
        want = {t.lower() for t in types if t.strip()}
        for index in range(self._list.count()):
            item = self._list.item(index)
            if item is None:
                continue
            item.setCheckState(
                Qt.CheckState.Checked
                if item.text().lower() in want
                else Qt.CheckState.Unchecked
            )

    def refresh_from_paths(self, paths: list[Path]) -> None:
        from repeatviewer.backends.gff_features import list_feature_types_in_annotation

        merged: set[str] = set()
        for raw in paths:
            path = raw.expanduser()
            if not path.is_file():
                continue
            merged.update(list_feature_types_in_annotation(path))
        if not merged:
            self._populate_default()
            return
        if merged == {"interval"}:
            self._set_types(["interval"], checked={"interval"})
            return
        merged_lower = {t.lower() for t in merged}
        if merged_lower.issubset({"snp", "indel", "mnp", "sv", "variant", "other"}):
            default = [
                t for t in sorted(merged, key=str.lower)
                if t.lower() in {"snp", "indel", "mnp", "sv"}
            ]
            if not default:
                default = sorted(merged, key=str.lower)
            prev = {t.lower() for t in self.selected_types()} or {t.lower() for t in default}
            self._set_types(sorted(merged, key=str.lower), checked=prev)
            return
        prev = {t.lower() for t in self.selected_types()} or {"gene"}
        self._set_types(sorted(merged, key=str.lower), checked=prev)


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


class TeDensityPlotDialog(QDialog):
    """TE density / genome coverage plot settings and input hit table."""

    def __init__(
        self,
        *,
        initial: TeGenomeMapOptions | None = None,
        chromosome_names: list[str] | None = None,
        chromosome_sizes: dict[str, int] | None = None,
        input_file: Path | None = None,
        consensus_input_file: Path | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Visualization")
        self.setMinimumWidth(720)
        self._opts = initial or TeGenomeMapOptions()

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
        )
        self._gene_input = TeMultiInputFilesWidget(
            self,
            placeholder="Optional GFF/GTF/BED/VCF annotation",
            browse_title="Select annotation file",
            file_filter="Annotations (*.gff *.gff3 *.gtf *.bed *.vcf *.vcf.gz);;All files (*)",
            add_button_label="Add another annotation file",
            add_button_tooltip=(
                "Add another GFF, GTF, BED, or VCF annotation file. "
                "Each file is drawn on its own track below TE tracks."
            ),
        )
        self._methylation_input = TeMultiInputFilesWidget(
            self,
            placeholder="Optional Bismark coverage (.cov / .cov.gz)",
            browse_title="Select Bismark coverage file",
            file_filter=(
                "Bismark coverage (*.cov *.cov.gz *bismark*.cov*);;All files (*)"
            ),
            add_button_label="Add another methylation file",
            add_button_tooltip=(
                "Add another Bismark coverage file. Each file adds four tracks per "
                "chromosome: total (genome), covered (.cov lines), covered fraction, "
                "and methylation level (0–1)."
            ),
        )
        self._focused_file_input: TeMultiInputFilesWidget | None = None
        for widget in (
            self._masking_input,
            self._consensus_input,
            self._gene_input,
            self._methylation_input,
        ):
            widget.field_activated.connect(self._on_file_input_field_activated)

        self._bin_kb = QSpinBox(self)
        self._bin_kb.setRange(5, 500)
        self._bin_kb.setValue(max(5, self._opts.bin_size_bp // 1000))
        self._bin_kb.setSuffix(" kb")
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
        self._methylation_coverage_min.setToolTip(
            "Reserved for future filtering. Covered tracks count every line in the "
            "Bismark .cov file whose position falls in each window."
        )
        self._split_strand = QCheckBox("Separate forward / reverse strands", self)
        self._split_strand.setChecked(bool(self._opts.split_strand_tracks))
        self._split_strand.setToolTip(
            "Draw separate forward- and reverse-strand tracks per chromosome for "
            "RepeatMasker .out TE hits and for GFF/GTF/BED/VCF annotation features."
        )
        self._collapse_overlaps = QCheckBox("Collapse overlapped TEs", self)
        self._collapse_overlaps.setChecked(bool(self._opts.collapse_overlaps))
        self._collapse_overlaps.setToolTip(
            "For total TE length per window, count each base pair at most once when TEs "
            "overlap. When color is by order, superfamily, or family, overlaps are collapsed "
            "separately within each TE type; with % identity coloring, all TEs in a window "
            "are collapsed together. With separate forward/reverse tracks, collapse applies "
            "within each strand track only."
        )
        self._show_te_diff = QCheckBox("Show TE difference track (2 TE sources)", self)
        self._show_te_diff.setChecked(bool(self._opts.show_te_diff_track))
        self._show_meth_diff = QCheckBox(
            "Show methylation difference track (2 methylation files)", self
        )
        self._show_meth_diff.setChecked(bool(self._opts.show_meth_diff_track))
        self._show_meth_diff.setToolTip(
            "When two Bismark coverage files are provided, add a diff track per "
            "chromosome comparing methylation level (first file above axis, second below)."
        )
        self._show_te_diff.setToolTip(
            "When exactly two TE hit-table sources are selected (masking files and/or "
            "consensus BLAST tables), add a per-chromosome diff track "
            "(exclusive TE bp: above = first source only, below = second only). "
            "Track order follows input file order top to bottom. "
            "This step is optional and can take extra time on large genomes."
        )

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
        self._chrom_picker = ChromosomeMultiSelectWidget(
            chromosome_names or [],
            chromosome_sizes=chromosome_sizes,
            preselected=self._opts.chromosomes_only,
            chromosome_sort=self._opts.chromosome_sort,
            parent=self,
        )
        sort_idx = self._chromosome_sort.findData(self._opts.chromosome_sort)
        if sort_idx >= 0:
            self._chromosome_sort.blockSignals(True)
            self._chromosome_sort.setCurrentIndex(sort_idx)
            self._chromosome_sort.blockSignals(False)
        self._chromosome_sort.currentIndexChanged.connect(self._on_chromosome_sort_changed)

        form_host = QWidget(self)
        form = QFormLayout()
        form_host.setLayout(form)

        form.addRow("TE masking file:", self._masking_input)
        form.addRow("TE consensus file:", self._consensus_input)
        form.addRow("Y axis content:", self._density)
        form.addRow("Color by:", self._color_by)
        form.addRow("", self._split_strand)
        form.addRow("", self._collapse_overlaps)
        form.addRow("", self._show_te_diff)
        form.addRow(_form_section_separator())
        form.addRow("Other annotation file:", self._gene_input)
        form.addRow("Display content:", self._annotation_features)
        form.addRow("Y axis content:", self._annotation_density)
        form.addRow(_form_section_separator())
        form.addRow("DNA methylation file:", self._methylation_input)
        form.addRow("", self._show_meth_diff)
        form.addRow("Methylation coverage min:", self._methylation_coverage_min)
        form.addRow(_form_section_separator())
        form.addRow("Window size:", self._bin_kb)
        form.addRow("Chromosome order:", self._chromosome_sort)
        self._chrom_picker.setToolTip(
            "Each chromosome is one row in the map. Many small scaffolds are slower "
            "than a few large ones — set a minimum size (Mb) and click Confirm."
        )
        form.addRow("Chromosomes to show:", self._chrom_picker)

        self._gene_input.paths_changed.connect(self._refresh_annotation_feature_types)
        self._refresh_annotation_feature_types()

        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(form_host)
        scroll.setMinimumHeight(320)
        scroll.setMaximumHeight(560)

        buttons = _run_cancel_button_box(self)
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addWidget(scroll)
        layout.addWidget(buttons)
        self.setLayout(layout)

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

    def gene_input_file(self) -> Path | None:
        paths = self.gene_input_files()
        return paths[0] if paths else None

    def options(self) -> TeGenomeMapOptions:
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


class MsaPlotDialog(QDialog):
    """MSA overview plot (nucleotide or protein) with custom residue colors."""

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
        self.setMinimumWidth(640)
        self._opts = initial or MsaPlotJobOptions()
        self._default_nuc = dict(DEFAULT_NUCLEOTIDE_COLORS)
        self._default_aa = dict(DEFAULT_AMINO_ACID_COLORS)

        self._input = TeMultiInputFilesWidget(
            self,
            initial=[input_file] if input_file is not None else None,
            allow_multiple=False,
            placeholder="Aligned FASTA (.fa / .fasta / .faa)",
            browse_title="Select MSA FASTA",
            file_filter="FASTA (*.fa *.fasta *.fas *.fna *.faa);;All files (*)",
        )
        self._focused_file_input: TeMultiInputFilesWidget | None = None
        self._input.field_activated.connect(self._on_file_input_field_activated)

        self._seq_type = QComboBox(self)
        self._seq_type.addItem("Auto-detect", MSA_PLOT_SEQ_AUTO)
        self._seq_type.addItem("Nucleotide (DNA/RNA)", MSA_PLOT_SEQ_NUCLEOTIDE)
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
        self._dpi.setToolTip("Image resolution / quality (higher = larger file).")

        self._seq_height = QDoubleSpinBox(self)
        self._seq_height.setRange(0.1, 5.0)
        self._seq_height.setDecimals(2)
        self._seq_height.setSingleStep(0.05)
        self._seq_height.setValue(float(self._opts.sequence_height_cm))
        self._seq_height.setSuffix(" cm")
        self._seq_height.setToolTip(
            "Figure height reserved for each sequence row (centimeters). "
            "Increase to make rows taller."
        )

        self._seq_distance = QDoubleSpinBox(self)
        self._seq_distance.setRange(0.0, 2.0)
        self._seq_distance.setDecimals(2)
        self._seq_distance.setSingleStep(0.01)
        self._seq_distance.setValue(float(self._opts.sequence_distance_cm))
        self._seq_distance.setSuffix(" cm")
        self._seq_distance.setToolTip(
            "White distance between consecutive sequence rows (centimeters). "
            "Set to 0 to remove spacing."
        )

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
        nuc_presets.addWidget(nuc_preset_label)
        nuc_tips = {
            "Default": "Standard A/C/G/T colours",
            "Color-blind": "Wong / Paul Tol palette (same as Annotation → nucleotide Color-blind)",
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

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.addRow("MSA FASTA:", self._input)
        form.addRow("Sequence type:", self._seq_type)
        form.addRow("Output format:", self._format)
        form.addRow("Quality:", self._dpi)
        form.addRow("Sequence height:", self._seq_height)
        form.addRow("Sequence distance:", self._seq_distance)

        self._colors_toggle = QPushButton("▶ Color", self)
        self._colors_toggle.setCheckable(True)
        self._colors_toggle.setChecked(False)
        self._colors_toggle.setCursor(Qt.CursorShape.PointingHandCursor)
        self._colors_toggle.setToolTip("Show or hide nucleotide / amino-acid colour settings")
        self._colors_toggle.toggled.connect(self._on_colors_section_toggled)

        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        colors_host = QWidget()
        colors_layout = QVBoxLayout(colors_host)
        colors_layout.setContentsMargins(0, 0, 0, 0)
        colors_layout.addWidget(nuc_box)
        colors_layout.addWidget(aa_box)
        scroll.setWidget(colors_host)
        scroll.setMinimumHeight(260)
        scroll.setVisible(False)
        self._colors_panel = scroll

        buttons = make_cancel_save_run_bar(self)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(20, 18, 20, 16)
        layout.setSpacing(12)
        layout.addLayout(form)
        layout.addWidget(self._colors_toggle, alignment=Qt.AlignmentFlag.AlignLeft)
        layout.addWidget(self._colors_panel, stretch=1)
        layout.addWidget(buttons)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _on_colors_section_toggled(self, checked: bool) -> None:
        self._colors_panel.setVisible(checked)
        self._colors_toggle.setText("▼ Color" if checked else "▶ Color")
        self.adjustSize()

    def _on_file_input_field_activated(self) -> None:
        sender = self.sender()
        if isinstance(sender, TeMultiInputFilesWidget):
            self._focused_file_input = sender

    def apply_path_from_panel(self, path: Path) -> bool:
        if self._focused_file_input is None:
            return self._input.set_path_at_focus(path)
        return self._focused_file_input.set_path_at_focus(path)

    def input_files(self) -> list[Path]:
        return self._input.paths()

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


class RepeatQuantifyDialog(QDialog):
    """Choose a library FASTA, one or more genome FASTAs, and the engine.

    Laid out like the BLAST / Align dialogs: right-aligned labels, fixed-width
    parameter fields, one tab per engine, and Restore Defaults + Cancel / Save /
    Run. Save stores the parameters without running.

    Both path fields start empty (nothing is guessed from the current folder,
    which could pick up a previous run's intermediate FASTAs). The
    dialog is shown non-modally, so the Folders and Files panels stay live:
    clicking a FASTA there drops it into whichever field last had focus
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
    _FIELD_WIDTH = 260  # same as the BLAST / Align parameter fields

    def __init__(
        self,
        *,
        start_dir: str = "",
        initial: RepeatQuantifyOptions | None = None,
        blast_databases: list | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Repeat quantify")
        # Saved BLAST databases (Add/Manage databases), to show which genomes
        # already have an index Repeat Quantify will reuse.
        self._blast_databases = list(blast_databases or [])
        # Wide enough for the labels plus a readable path, Browse and Remove;
        # any wider only leaves empty space beside the fixed-width parameters.
        self.setMinimumWidth(680)
        self._start_dir = start_dir
        opts = initial or RepeatQuantifyOptions()

        # Path fields rather than dropdowns, so the dialog can stay non-modal and
        # be filled by clicking files in the Files panel (see apply_path_from_panel).
        self._library_input = TeMultiInputFilesWidget(
            self,
            allow_multiple=False,
            start_dir=start_dir,
            placeholder="TE library FASTA (monomers / satellites)",
            browse_title="Select TE library FASTA",
            file_filter=self._FASTA_FILTER,
        )
        self._genome_input = TeMultiInputFilesWidget(
            self,
            start_dir=start_dir,
            placeholder="Genome FASTA (scaffolds / assembly)",
            browse_title="Select genome FASTA",
            file_filter=self._FASTA_FILTER,
            add_button_label="Add another genome",
            add_button_tooltip=(
                "Quantify the same library against another genome. "
                "Each genome produces its own report."
            ),
        )
        self._focused_file_input: TeMultiInputFilesWidget | None = None
        for widget in (self._library_input, self._genome_input):
            widget.field_activated.connect(self._on_file_input_field_activated)

        inputs_form = self._form()
        inputs_form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        inputs_form.setContentsMargins(0, 0, 0, 0)
        inputs_form.addRow("Library (monomers / satellites):", self._library_input)
        inputs_form.addRow("Genome (scaffolds / assembly):", self._genome_input)
        self._db_status = QLabel(self)
        self._db_status.setWordWrap(True)
        self._db_status.setTextFormat(Qt.TextFormat.PlainText)
        # A word-wrapped label asks for a narrow width; let it use the whole field.
        self._db_status.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred
        )
        self._db_status.setVisible(False)
        inputs_form.addRow("", self._db_status)
        inputs_form.addRow("Output folder:", self._make_output_controls())

        # Tabs are the engine choice; each tab holds that engine's parameters.
        # "Both" repeats the controls, kept in sync with the single-engine tabs.
        self._tabs = QTabWidget(self)
        self._blast_a = self._make_blast_controls()
        self._rm_b = self._make_rm_controls()
        self._blast_c = self._make_blast_controls()
        self._rm_d = self._make_rm_controls()
        self._link_controls(self._blast_a, self._blast_c)
        self._link_controls(self._rm_b, self._rm_d)
        self._tabs.addTab(
            self._params_tab(self._blast_rows(self._blast_a)), "BLASTn"
        )
        self._tabs.addTab(
            self._params_tab(self._rm_rows(self._rm_b)), "RepeatMasker"
        )
        both_rows = (
            [("BLASTn", None)]
            + self._blast_rows(self._blast_c)
            + [("RepeatMasker", None)]
            + self._rm_rows(self._rm_d)
        )
        self._tabs.addTab(self._params_tab(both_rows), "Both")
        self._apply_options(opts)
        self._genome_input.paths_changed.connect(self._update_db_status)
        self._tabs.currentChanged.connect(self._update_db_status)
        self._blast_a["save"].toggled.connect(self._update_db_status)
        self._update_db_status()

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore],
            on_accept=self._on_footer_accept,
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(14)
        layout.addLayout(inputs_form)
        layout.addWidget(self._tabs)
        layout.addWidget(self._button_box)

    def _make_output_controls(self) -> QWidget:
        """Output folder field + Browse, with a preview of the run folder below."""
        self._output_focused = False
        self._output_row = PathBrowseRow(
            self,
            placeholder="Auto — a new folder next to each genome",
            browse_title="Select output folder",
            mode="directory",
        )
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
        box = QVBoxLayout(host)
        box.setContentsMargins(0, 0, 0, 0)
        box.setSpacing(3)
        box.addWidget(self._output_row)
        box.addWidget(self._output_preview)

        self._genome_input.paths_changed.connect(self._update_output_preview)
        self._update_output_preview()
        return host

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

    @staticmethod
    def _form(parent: QWidget | None = None) -> QFormLayout:
        form = QFormLayout(parent) if parent is not None else QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        return form

    def _make_blast_controls(self) -> dict:
        word = QSpinBox(self)
        word.setRange(4, 24)
        word.setToolTip(
            "BLASTn -word_size. Lower is more sensitive but much slower on "
            "big/repetitive genomes (7 can hang; 11 is a good balance)."
        )
        evalue = QDoubleSpinBox(self)
        configure_report_evalue_spinbox(evalue)
        evalue.setToolTip("BLASTn -evalue: report hits at or below this E-value.")
        minid = QDoubleSpinBox(self)
        minid.setRange(0.0, 100.0)
        minid.setDecimals(0)
        minid.setSuffix(" %")
        minid.setToolTip("BLASTn -perc_identity. 0 = no minimum.")
        mincov = QDoubleSpinBox(self)
        mincov.setRange(0.0, 100.0)
        mincov.setDecimals(0)
        mincov.setSuffix(" %")
        mincov.setToolTip(
            "BLASTn -qcov_hsp_perc: minimum query (monomer) coverage per hit. "
            "0 = no minimum."
        )
        for w in (word, evalue, minid, mincov):
            w.setFixedWidth(self._FIELD_WIDTH)
            w.setAlignment(Qt.AlignmentFlag.AlignRight)
        save = QCheckBox("Save new indexes as BLAST databases", self)
        save.setToolTip(
            "Keep each genome's BLAST index in Add/Manage databases so later "
            "runs and other tools reuse it."
        )
        return {
            "word": word, "evalue": evalue, "minid": minid, "mincov": mincov,
            "save": save,
        }

    def _make_rm_controls(self) -> dict:
        div = QSpinBox(self)
        div.setRange(1, 100)
        div.setSuffix(" %")
        div.setFixedWidth(self._FIELD_WIDTH)
        div.setAlignment(Qt.AlignmentFlag.AlignRight)
        div.setToolTip(
            "RepeatMasker -div: report only copies less than this % diverged "
            "(higher captures older copies, but noisier)."
        )
        sens = QComboBox(self)
        sens.addItem("Sensitive (-s)", "sensitive")
        sens.addItem("Default", "default")
        sens.addItem("Quick (-q)", "quick")
        sens.setFixedWidth(self._FIELD_WIDTH)
        sens.setToolTip(
            "RepeatMasker search sensitivity. Sensitive is slowest; Quick is "
            "fastest on large genomes."
        )
        return {"div": div, "sens": sens}

    @staticmethod
    def _link_controls(first: dict, second: dict) -> None:
        """Mirror each control of ``first`` into ``second`` and back."""
        for key, a in first.items():
            b = second[key]
            if isinstance(a, QComboBox):
                a.currentIndexChanged.connect(b.setCurrentIndex)
                b.currentIndexChanged.connect(a.setCurrentIndex)
            elif isinstance(a, QCheckBox):
                a.toggled.connect(b.setChecked)
                b.toggled.connect(a.setChecked)
            else:
                a.valueChanged.connect(b.setValue)
                b.valueChanged.connect(a.setValue)

    @staticmethod
    def _blast_rows(c: dict) -> list:
        return [
            ("-word_size:", c["word"]),
            ("-evalue:", c["evalue"]),
            ("-perc_identity:", c["minid"]),
            ("-qcov_hsp_perc:", c["mincov"]),
            (
                "BLAST index:",
                field_with_help(
                    c["save"],
                    "BLAST index",
                    "Repeat Quantify searches each genome's BLAST index. If the "
                    "genome is already in Add/Manage databases (and the FASTA has "
                    "not changed since), that index is reused and makeblastdb is "
                    "skipped.\n\nWhen this box is ticked, a genome without a saved "
                    "index is indexed once and added to Add/Manage databases, so "
                    "the next run, and the BLAST, TE-Aid and Extend tools, reuse "
                    "it. When unticked, the index is temporary and rebuilt when "
                    "the system clears its temp folder.",
                ),
            ),
        ]

    @staticmethod
    def _rm_rows(c: dict) -> list:
        return [
            ("-div (divergence cutoff):", c["div"]),
            ("Sensitivity:", c["sens"]),
        ]

    def _params_tab(self, rows: list) -> QWidget:
        """Build a tab from (label, widget) rows; a None widget = bold header."""
        tab = QWidget()
        form = self._form(tab)
        form.setContentsMargins(16, 14, 16, 14)
        for label, widget in rows:
            if widget is None:
                header = QLabel(label, self)
                font = header.font()
                font.setBold(True)
                header.setFont(font)
                form.addRow(header)
            else:
                form.addRow(label, widget)
        return tab

    def _apply_options(self, opts: RepeatQuantifyOptions) -> None:
        # Setting the single-engine controls updates the "Both" copies too.
        self._blast_a["word"].setValue(int(opts.blast_word_size))
        self._blast_a["evalue"].setValue(float(opts.blast_evalue))
        self._blast_a["minid"].setValue(float(opts.blast_min_identity))
        self._blast_a["mincov"].setValue(float(opts.blast_min_coverage))
        self._blast_a["save"].setChecked(bool(opts.blast_save_index))
        self._rm_b["div"].setValue(int(opts.rm_div))
        idx = self._rm_b["sens"].findData(opts.rm_sensitivity)
        self._rm_b["sens"].setCurrentIndex(max(0, idx))
        engine = opts.engine if opts.engine in self._ENGINES else "both"
        self._tabs.setCurrentIndex(self._ENGINES.index(engine))

    def _restore_defaults(self) -> None:
        """Reset the parameters (not the chosen files) to their defaults."""
        self._apply_options(RepeatQuantifyOptions())

    def options(self) -> RepeatQuantifyOptions:
        return RepeatQuantifyOptions(
            engine=self.engine(),
            blast_word_size=int(self._blast_a["word"].value()),
            blast_evalue=float(self._blast_a["evalue"].value()),
            blast_min_identity=float(self._blast_a["minid"].value()),
            blast_min_coverage=float(self._blast_a["mincov"].value()),
            rm_div=int(self._rm_b["div"].value()),
            rm_sensitivity=str(self._rm_b["sens"].currentData()),
            blast_save_index=self.save_blast_index(),
        )

    def save_blast_index(self) -> bool:
        """Save a genome's new BLAST index to Add/Manage databases."""
        return bool(self._blast_a["save"].isChecked())

    def blast_database_status(self) -> list[str]:
        """One line per genome: whether its BLAST index is reused or built."""
        from repeatviewer.backends.database_prepare import saved_blast_database_for

        lines = []
        for target in self.targets():
            saved = saved_blast_database_for(target, self._blast_databases)
            if saved is not None:
                lines.append(f"{target.name}: uses saved BLAST database \u201c{saved.name}\u201d")
            elif self.save_blast_index():
                lines.append(
                    f"{target.name}: no saved BLAST database yet \u2014 indexed once "
                    "and added to Add/Manage databases"
                )
            else:
                lines.append(
                    f"{target.name}: no saved BLAST database \u2014 temporary index for this run"
                )
        return lines

    def _update_db_status(self, *_args) -> None:
        uses_blast = self.engine() in ("blastn", "both")
        lines = self.blast_database_status() if uses_blast else []
        self._db_status.setText("\n".join(lines))
        self._db_status.setVisible(bool(lines))

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
        library = self.library()
        targets = self.targets()
        if library is None or not targets:
            QMessageBox.information(
                self, "Repeat quantify",
                "Choose both a library FASTA and at least one genome FASTA "
                "(click a file in the Files panel, or use Browse).",
            )
            return
        missing = [p for p in [library, *targets] if not p.is_file()]
        if missing:
            QMessageBox.warning(
                self, "Repeat quantify",
                "File not found:\n" + "\n".join(str(p) for p in missing[:5]),
            )
            return
        if library in targets:
            QMessageBox.warning(
                self, "Repeat quantify",
                "Genome and library must be different files.",
            )
            return
        if len(set(targets)) != len(targets):
            QMessageBox.warning(
                self, "Repeat quantify",
                "The same genome is listed more than once.",
            )
            return
        out = self.output_dir()
        if out is not None:
            if not out.is_absolute():
                QMessageBox.warning(
                    self, "Repeat quantify",
                    "Enter the full path of the output folder, or use Browse.",
                )
                return
            if out.exists() and not out.is_dir():
                QMessageBox.warning(
                    self, "Repeat quantify",
                    f"The output path is a file, not a folder:\n{out}",
                )
                return
        self.accept()

    def targets(self) -> list[Path]:
        """Genome FASTAs to quantify against — one report per genome."""
        return self._genome_input.paths()

    def target(self) -> Path | None:
        targets = self.targets()
        return targets[0] if targets else None

    def library(self) -> Path | None:
        paths = self._library_input.paths()
        return paths[0] if paths else None

    def engine(self) -> str:
        return ("blastn", "repeatmasker", "both")[self._tabs.currentIndex()]

    def params(self) -> dict:
        return self.options().params()
