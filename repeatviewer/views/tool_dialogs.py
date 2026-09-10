"""Dialogs for tool-specific parameters (right-click / custom runs)."""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from PySide6.QtCore import QEvent, Qt, QTimer, Signal
from PySide6.QtGui import QColor, QStandardItemModel
from PySide6.QtWidgets import (
    QAbstractItemView,
    QCheckBox,
    QColorDialog,
    QComboBox,
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
    DEFAULT_FLEXIDOT_LCS_SHADING,
    DEFAULT_FLEXIDOT_OUTPUT_FORMAT,
    DEFAULT_FLEXIDOT_WORDSIZE,
    DEFAULT_FLEXIDOT_MODES,
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
    DEFAULT_TEAID_EVALUE,
    DEFAULT_TEAID_X_AXIS_LABEL,
    DEFAULT_TEAID_SELF_BLAST_WORD_SIZE,
    DEFAULT_TEAID_HIT_LEN_HIGH_PCT,
    DEFAULT_TEAID_HIT_LEN_MID_PCT,
    DEFAULT_TEAID_HIT_COLOR_HIGH,
    DEFAULT_TEAID_HIT_COLOR_MID,
    MAX_TEAID_SELF_BLAST_WORD_SIZE,
    MIN_TEAID_SELF_BLAST_WORD_SIZE,
    TEAID_RIGHT_CLICK_RUN,
    TEAID_RIGHT_CLICK_DIALOG,
    clamp_teaid_hit_len_pct,
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
)
from repeatviewer.plugins.custom_command import (
    CommandVarDef,
    UserCommandPluginSpec,
    get_user_command_plugin_spec,
    load_user_command_plugin_specs,
    normalize_var_mode,
    plugin_id_from_label,
    template_variables,
    user_command_plugin_ids,
)


class _CommandVarRow(QWidget):
    remove_requested = Signal(object)
    output_toggled = Signal(object, bool)

    def __init__(self, initial: CommandVarDef | None = None, *, is_output: bool = False, parent=None) -> None:
        super().__init__(parent)
        var = initial or CommandVarDef(name="", mode="file", value="")
        self._name = QLineEdit(var.name, self)
        self._name.setPlaceholderText("variable_name")
        self._mode = QComboBox(self)
        self._mode.addItem("File", "file")
        self._mode.addItem("Folder", "folder")
        self._mode.addItem("Text", "text")
        self._mode.addItem("Number", "number")
        self._mode.setMinimumWidth(76)
        mode_key = normalize_var_mode(var.mode)
        index = self._mode.findData(mode_key)
        if index >= 0:
            self._mode.setCurrentIndex(index)
        self._default = QLineEdit(var.value, self)
        self._default.setPlaceholderText("Default value (optional)")
        self._is_output = QCheckBox("Output", self)
        self._is_output.setChecked(is_output)
        self._is_output.setMinimumWidth(self._is_output.fontMetrics().horizontalAdvance("Output") + 28)
        self._is_output.toggled.connect(self._on_output_toggled)
        self._is_output.toggled.connect(self._refresh_output_fields)
        remove_btn = QPushButton("Remove", self)
        remove_btn.clicked.connect(lambda: self.remove_requested.emit(self))
        self._remove_btn = remove_btn
        remove_btn.setMinimumWidth(72)

        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(8)
        row.addWidget(self._name, 1)
        row.addWidget(self._mode, 0)
        row.addWidget(self._default, 1)
        row.addWidget(self._is_output, 0)
        row.addWidget(remove_btn, 0)
        self._refresh_output_fields()

    def _on_output_toggled(self, checked: bool) -> None:
        if checked:
            self.output_toggled.emit(self, True)

    def set_output_checked(self, checked: bool) -> None:
        self._is_output.blockSignals(True)
        self._is_output.setChecked(checked)
        self._is_output.blockSignals(False)
        self._refresh_output_fields()

    def _refresh_output_fields(self) -> None:
        is_output = self._is_output.isChecked()
        if is_output:
            mode = normalize_var_mode(str(self._mode.currentData()))
            if mode not in {"file", "folder"}:
                index = self._mode.findData("file")
                if index >= 0:
                    self._mode.setCurrentIndex(index)

    def definition(self) -> tuple[CommandVarDef | None, bool]:
        name = self._name.text().strip()
        if not name:
            return None, False
        is_output = self._is_output.isChecked()
        var = CommandVarDef(
            name=name,
            mode=str(self._mode.currentData()),
            value=self._default.text().strip(),
        )
        return var, is_output

    def mark_removable(self, removable: bool) -> None:
        self._remove_btn.setVisible(removable)


class CommandPluginDesignerDialog(QDialog):
    """Design a user-defined command plugin from a command template."""

    def __init__(
        self,
        *,
        initial: UserCommandPluginSpec | None = None,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Design your plugin")
        self.setMinimumWidth(760)
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
        self._label.setPlaceholderText("Plugin name (must be unique)")
        self._tooltip = QTextEdit(self)
        self._tooltip.setPlaceholderText("Short description (optional)")
        self._tooltip.setPlainText(spec.tooltip)
        self._tooltip.setFixedHeight(52)
        self._tooltip.setLineWrapMode(QTextEdit.LineWrapMode.WidgetWidth)
        self._command = QTextEdit(self)
        self._command.setFixedHeight(120)
        self._command.setPlaceholderText(
            "Example:\nFastGA ${genome1} ${genome2} -1:${output}"
        )
        self._command.setPlainText(spec.command_template)
        self._rows_host = QWidget(self)
        self._rows_layout = QVBoxLayout(self._rows_host)
        self._rows_layout.setContentsMargins(0, 0, 0, 0)
        self._rows_layout.setSpacing(6)
        self._rows: list[_CommandVarRow] = []
        out_set = set(spec.output_vars)
        for var in spec.var_defs:
            self._add_var_row(var, is_output=var.name in out_set)
        if not self._rows:
            self._add_var_row()

        detect_btn = QPushButton("Detect variables from command", self)
        detect_btn.clicked.connect(self._detect_variables)
        add_var_btn = QPushButton("Add variable", self)
        add_var_btn.clicked.connect(lambda: self._add_var_row())

        form = QFormLayout()
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        form.addRow("Plugin name:", self._label)
        form.addRow("Description:", self._tooltip)
        form.addRow("Command template:", self._command)
        btn_row = QWidget(self)
        btn_layout = QHBoxLayout(btn_row)
        btn_layout.setContentsMargins(0, 0, 0, 0)
        btn_layout.addWidget(detect_btn)
        btn_layout.addStretch(1)
        form.addRow("", btn_row)
        form.addRow("Define variables:", self._rows_host)
        add_row = QWidget(self)
        add_layout = QHBoxLayout(add_row)
        add_layout.setContentsMargins(0, 0, 0, 0)
        add_layout.addWidget(add_var_btn)
        add_layout.addStretch(1)
        form.addRow("", add_row)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addLayout(form)
        layout.addWidget(buttons)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        QTimer.singleShot(0, self._label.setFocus)

    def _add_var_row(self, var: CommandVarDef | None = None, *, is_output: bool = False) -> None:
        if not hasattr(self, "_rows_host"):
            return
        row = _CommandVarRow(var, is_output=is_output, parent=self._rows_host)
        row.remove_requested.connect(self._remove_var_row)
        row.output_toggled.connect(self._on_output_toggled)
        self._rows.append(row)
        self._rows_layout.addWidget(row)
        self._refresh_remove_buttons()

    def _remove_var_row(self, row: _CommandVarRow) -> None:
        if row not in self._rows:
            return
        if len(self._rows) <= 1:
            row._name.clear()
            row._default.clear()
            row._is_output.setChecked(False)
            row._mode.setCurrentIndex(row._mode.findData("file"))
            row._refresh_output_fields()
            return
        self._rows.remove(row)
        self._rows_layout.removeWidget(row)
        row.setParent(None)
        row.deleteLater()
        self._refresh_remove_buttons()

    def _refresh_remove_buttons(self) -> None:
        removable = len(self._rows) > 1
        for row in self._rows:
            row.mark_removable(removable)

    def _on_output_toggled(self, source_row: _CommandVarRow, checked: bool) -> None:
        if not checked:
            return
        for row in self._rows:
            if row is not source_row:
                row.set_output_checked(False)

    def _detect_variables(self) -> None:
        current_defs = {}
        outputs = set()
        for row in self._rows:
            var, is_output = row.definition()
            if var is None:
                continue
            current_defs[var.name] = var
            if is_output:
                outputs.add(var.name)
        vars_in_cmd = template_variables(self._command.toPlainText())
        for name in vars_in_cmd:
            if name not in current_defs:
                current_defs[name] = CommandVarDef(name=name, mode="file", value="")
        self._rebuilding_rows = True
        try:
            for row in list(self._rows):
                self._rows_layout.removeWidget(row)
                row.setParent(None)
                row.deleteLater()
            self._rows = []
            for name in sorted(current_defs):
                self._add_var_row(current_defs[name], is_output=name in outputs)
        finally:
            self._rebuilding_rows = False
        if not self._rows:
            self._add_var_row()

    def spec(self) -> UserCommandPluginSpec:
        outputs: list[str] = []
        var_defs: list[CommandVarDef] = []
        for row in self._rows:
            var, is_output = row.definition()
            if var is None:
                continue
            var_defs.append(var)
            if is_output:
                outputs.append(var.name)
        return UserCommandPluginSpec(
            id=plugin_id_from_label(self._label.text().strip()),
            label=self._label.text().strip(),
            command_template=self._command.toPlainText().strip(),
            output_vars=outputs[:1],
            var_defs=var_defs,
            fasta_only=self._initial_fasta_only,
            tooltip=self._tooltip.toPlainText().strip(),
        )

    def accept(self) -> None:
        spec = self.spec()
        if not spec.label:
            QMessageBox.warning(self, "Design your plugin", "Plugin name is required.")
            return
        label_key = spec.label.strip().lower()
        for existing in load_user_command_plugin_specs():
            if existing.id == self._initial_id:
                continue
            if existing.label.strip().lower() == label_key:
                QMessageBox.warning(
                    self,
                    "Design your plugin",
                    f"A plugin named “{spec.label}” already exists. Choose a different name.",
                )
                return
        if not spec.command_template:
            QMessageBox.warning(self, "Design your plugin", "Command template is required.")
            return
        output_count = 0
        for row in self._rows:
            var, is_output = row.definition()
            if var is None:
                continue
            if is_output:
                output_count += 1
                if normalize_var_mode(var.mode) not in {"file", "folder"}:
                    QMessageBox.warning(
                        self,
                        "Design your plugin",
                        f"Output variable “{var.name}” must use type File or Folder.",
                    )
                    return
        if output_count > 1:
            QMessageBox.warning(
                self,
                "Design your plugin",
                "Only one variable can be marked as Output.",
            )
            return
        super().accept()


class PluginManageDialog(QDialog):
    """List installed plugins; add, edit, or delete user-designed command plugins."""

    def __init__(
        self,
        plugin_registry,
        *,
        on_save,
        on_delete,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Manage plugins")
        self.setMinimumWidth(720)
        self.setMinimumHeight(420)
        self._registry = plugin_registry
        self._on_save = on_save
        self._on_delete = on_delete

        self._list = QListWidget(self)
        self._list.currentItemChanged.connect(self._update_buttons)

        self._add_btn = QPushButton("Add…", self)
        self._edit_btn = QPushButton("Edit…", self)
        self._delete_btn = QPushButton("Delete…", self)
        self._add_btn.clicked.connect(self._add_plugin)
        self._edit_btn.clicked.connect(self._edit_plugin)
        self._delete_btn.clicked.connect(self._delete_plugin)
        self._list.itemDoubleClicked.connect(lambda _item: self._edit_plugin())

        btn_row = QHBoxLayout()
        btn_row.addWidget(self._add_btn)
        btn_row.addWidget(self._edit_btn)
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
        custom_ids = user_command_plugin_ids()
        plugins = self._registry.all_plugins()
        if not plugins:
            item = QListWidgetItem("(No plugins loaded)")
            item.setFlags(Qt.ItemFlag.NoItemFlags)
            self._list.addItem(item)
            return
        for plugin in plugins:
            kind = "Custom" if plugin.id in custom_ids else "Installed"
            tooltip = (plugin.tooltip or "").strip()
            text = f"{plugin.label} — {kind}"
            item = QListWidgetItem(text)
            item.setData(Qt.ItemDataRole.UserRole, plugin.id)
            item.setData(Qt.ItemDataRole.UserRole + 1, plugin.id in custom_ids)
            if tooltip:
                item.setToolTip(tooltip)
            self._list.addItem(item)

    def _selected(self) -> tuple[str, bool] | None:
        item = self._list.currentItem()
        if item is None:
            return None
        plugin_id = item.data(Qt.ItemDataRole.UserRole)
        is_custom = bool(item.data(Qt.ItemDataRole.UserRole + 1))
        if not plugin_id:
            return None
        return str(plugin_id), is_custom

    def _update_buttons(self) -> None:
        selected = self._selected()
        has_custom = selected is not None and selected[1]
        self._edit_btn.setEnabled(has_custom)
        self._delete_btn.setEnabled(has_custom)

    def _add_plugin(self) -> None:
        dialog = CommandPluginDesignerDialog(parent=self)
        if not exec_modal_dialog(dialog):
            return
        self._on_save(dialog.spec())
        self._populate()

    def _edit_plugin(self) -> None:
        selected = self._selected()
        if selected is None or not selected[1]:
            return
        plugin_id = selected[0]
        spec = get_user_command_plugin_spec(plugin_id)
        if spec is None:
            QMessageBox.warning(self, "Manage plugins", f"Could not load plugin {plugin_id!r}.")
            return
        dialog = CommandPluginDesignerDialog(initial=spec, parent=self)
        if not exec_modal_dialog(dialog):
            return
        self._on_save(dialog.spec(), previous_id=plugin_id)
        self._populate()

    def _delete_plugin(self) -> None:
        selected = self._selected()
        if selected is None or not selected[1]:
            return
        plugin_id = selected[0]
        spec = get_user_command_plugin_spec(plugin_id)
        label = spec.label if spec is not None else plugin_id
        answer = QMessageBox.question(
            self,
            "Manage plugins",
            f"Delete custom plugin “{label}” ({plugin_id})?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        if not self._on_delete(plugin_id):
            QMessageBox.warning(self, "Manage plugins", f"Could not delete plugin {plugin_id!r}.")
            return
        self._populate()
        self._update_buttons()


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


class ExtractSequencesDialog(QDialog):
    """Pick records out of a multi-FASTA."""

    def __init__(self, *, records: list[tuple[str, str, int]], source_name: str, parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Extract Sequences from File")
        self.resize(560, 460)

        self._list = QListWidget(self)
        self._list.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)
        for rec_id, description, length in records:
            label = f"{rec_id}    ({length:,} bp)"
            if description:
                label += f"   {description}"
            item = QListWidgetItem(label, self._list)
            item.setData(Qt.ItemDataRole.UserRole, rec_id)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Unchecked)

        self._single_file = QCheckBox("Write all selected into one file", self)

        select_all = QPushButton("Select all", self)
        select_none = QPushButton("Select none", self)
        select_all.clicked.connect(lambda: self._set_all(Qt.CheckState.Checked))
        select_none.clicked.connect(lambda: self._set_all(Qt.CheckState.Unchecked))

        top = QHBoxLayout()
        top.addWidget(QLabel(f"Sequences in {source_name}:", self))
        top.addStretch(1)
        top.addWidget(select_all)
        top.addWidget(select_none)

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
        layout.addLayout(top)
        layout.addWidget(self._list, stretch=1)
        layout.addWidget(self._single_file)
        layout.addWidget(self._summary)
        layout.addWidget(self._buttons)

        self._list.itemChanged.connect(self._sync)
        self._sync()

    def _set_all(self, state) -> None:
        for index in range(self._list.count()):
            self._list.item(index).setCheckState(state)

    def _sync(self) -> None:
        chosen = self.selected_ids()
        self._summary.setText(
            f"{len(chosen)} of {self._list.count()} selected"
            if chosen else "Nothing selected"
        )
        self._buttons.button(QDialogButtonBox.StandardButton.Ok).setEnabled(bool(chosen))

    def selected_ids(self) -> set[str]:
        return {
            self._list.item(i).data(Qt.ItemDataRole.UserRole)
            for i in range(self._list.count())
            if self._list.item(i).checkState() == Qt.CheckState.Checked
        }

    def single_file(self) -> bool:
        return self._single_file.isChecked()


class FastGaParamsDialog(QDialog):
    """FastGA whole-genome alignment parameters."""

    def __init__(
        self,
        *,
        initial=None,
        file_count: int = 1,
        title: str = "FastGA",
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

        self._threads = QSpinBox(self)
        self._threads.setRange(1, 128)
        self._threads.setValue(max(1, opts.threads))

        self._identity = QDoubleSpinBox(self)
        self._identity.setRange(0.5, 1.0)
        self._identity.setSingleStep(0.05)
        self._identity.setDecimals(2)
        self._identity.setValue(opts.min_identity)

        self._min_length = QSpinBox(self)
        # Below ~100 FastGA fails: an alignment must contain a 40 bp seed and
        # chains must cover at least 85 bases.
        self._min_length.setRange(100, 1_000_000)
        self._min_length.setSingleStep(50)
        self._min_length.setValue(max(100, opts.min_length))
        self._min_length.setToolTip(
            "Shortest alignment to report, in bp. Cannot go below 100 -- "
            "FastGA rejects limits shorter than its seed length.\n"
            "Raise this on satellite arrays, where short matches can swamp "
            "the output."
        )

        self._seed_freq = QSpinBox(self)
        self._seed_freq.setRange(1, 100_000)
        self._seed_freq.setSingleStep(50)
        self._seed_freq.setValue(getattr(opts, "seed_freq", 200))
        self._seed_freq.setToolTip(
            "FastGA's -f: seeds occurring more often than this are discarded as "
            "'repetitive' before alignment.\n"
            "FastGA's own default is 10, which removes satellite arrays "
            "entirely, since every monomer k-mer occurs far more than 10 times.\n"
            "Raise it to see repeat structure; higher values are slower."
        )

        self._symmetric = QCheckBox("Symmetric seeding (-S)", self)
        self._symmetric.setChecked(opts.symmetric)
        self._symmetric.setToolTip(
            "Uses the seeds of both genomes so the result does not depend on "
            "argument order. Recommended when the repeat structure is the "
            "point; slower, and not needed when looking for synteny."
        )

        self._keep_index = QCheckBox("Keep genome index for reuse (-k)", self)
        self._keep_index.setChecked(opts.keep_index)
        self._keep_index.setToolTip(
            "Keeps the .1gdb / .gidx files beside the FASTA so later runs skip "
            "the index build, which is the slow part of a first run."
        )

        for field in (
            self._format,
            self._threads,
            self._identity,
            self._min_length,
            self._seed_freq,
        ):
            field.setFixedWidth(240)

        mode = QLabel(self)
        mode.setWordWrap(True)
        mode.setStyleSheet("color: palette(mid);")
        if file_count >= 2:
            mode.setText(
                "Two files selected: they will be compared against each other."
            )
        else:
            mode.setText(
                "One file selected: it will be compared against itself, which "
                "finds repetitive regions. Self matches are excluded."
            )

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        form.addRow("Output format:", self._format)
        form.addRow("Min identity:", self._identity)
        form.addRow("Min alignment length:", self._min_length)
        form.addRow("Repeat cutoff (-f):", self._seed_freq)
        form.addRow("Threads:", self._threads)
        form.addRow("", self._symmetric)
        form.addRow("", self._keep_index)

        self._button_box = make_cancel_save_run_bar(self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(14)
        layout.addWidget(mode)
        layout.addLayout(form)
        layout.addWidget(self._button_box)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def options(self):
        from repeatviewer.job_options import FastGaOptions

        return FastGaOptions(
            output_format=str(self._format.currentData()),
            threads=int(self._threads.value()),
            min_identity=float(self._identity.value()),
            min_length=int(self._min_length.value()),
            seed_freq=int(self._seed_freq.value()),
            symmetric=self._symmetric.isChecked(),
            keep_index=self._keep_index.isChecked(),
        )


class ExtendParamsDialog(QDialog):
    """Extend MSA ends.

    Left-click the toolbar button opens this dialog. The Parameters tab holds
    the distances and genome; the Options tab sets what right-clicking the
    toolbar button does.
    """

    _FASTA_FILTER = "FASTA (*.fa *.fasta *.fna *.fas);;All files (*)"

    def __init__(
        self,
        *,
        initial: ExtendOptions | None = None,
        title: str = "MSA Extend parameters",
        pan_genome_names: list[str] | None = None,
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

        # Genome selector (single-genome mode only): empty means use the
        # project genome from Define Genome.
        self._genome_combo = QComboBox()
        self._genome_combo.setFixedWidth(320)
        self._genome_combo.addItem("Default (project genome)", "")
        project_genome_str = (
            str(project_genome.expanduser().resolve())
            if project_genome is not None and project_genome.is_file()
            else ""
        )
        if project_genome_str:
            self._genome_combo.addItem(project_genome.name, project_genome_str)
        for path in genome_history or []:
            if not path.is_file():
                continue
            resolved = str(path.expanduser().resolve())
            if resolved == project_genome_str:
                continue
            self._genome_combo.addItem(path.name, resolved)
        if opts.genome_path:
            index = self._genome_combo.findData(opts.genome_path)
            if index >= 0:
                self._genome_combo.setCurrentIndex(index)
            else:
                self._genome_combo.addItem(Path(opts.genome_path).name, opts.genome_path)
                self._genome_combo.setCurrentIndex(self._genome_combo.count() - 1)

        genome_browse = QPushButton("Browse…", self)
        genome_browse.setFixedWidth(88)
        genome_browse.clicked.connect(self._browse_genome)

        genome_row = QWidget(self)
        genome_layout = QHBoxLayout(genome_row)
        genome_layout.setContentsMargins(0, 0, 0, 0)
        genome_layout.setSpacing(6)
        genome_layout.addWidget(self._genome_combo, stretch=1)
        genome_layout.addWidget(genome_browse)

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
            params_form.addRow("Genome:", genome_row)
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
        self._genome_combo.setCurrentIndex(0)
        self._right_click_buttons[EXTEND_RIGHT_CLICK_RUN].setChecked(True)

    def _browse_genome(self) -> None:
        chosen, _ = QFileDialog.getOpenFileName(
            self,
            "Select genome FASTA",
            str(Path.home()),
            self._FASTA_FILTER,
        )
        if not chosen:
            return
        resolved = str(Path(chosen).expanduser().resolve())
        index = self._genome_combo.findData(resolved)
        if index >= 0:
            self._genome_combo.setCurrentIndex(index)
        else:
            self._genome_combo.addItem(Path(chosen).name, resolved)
            self._genome_combo.setCurrentIndex(self._genome_combo.count() - 1)

    def options(self) -> ExtendOptions:
        pan_name = None
        if self._pan_genome_names:
            pan_name = self._pan_genome.currentText().strip() or None
        genome_path = None
        if not self._pan_genome_names:
            genome_path = self._genome_combo.currentData() or None
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
        self.setWindowTitle("Blast parameters")
        self._project = project
        self._for_run = for_run
        opts = initial or BlastOptions()
        self._initial_program = opts.program

        self._tabs = QTabWidget()
        self._blastn_database = QComboBox()
        self._blastp_database = QComboBox()
        for combo in (self._blastn_database, self._blastp_database):
            combo.setMinimumWidth(220)

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
        self._blastx_database.setMinimumWidth(220)

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
            self._max_msa,
            self._top_msa,
            self._min_bit_score,
        ):
            field.setFixedWidth(150)
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
        rc_box_layout.addWidget(
            QLabel("Right-clicking the BLAST button runs, with saved settings:", self)
        )
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
        manage = QPushButton("Add/Manage databases…")
        manage.clicked.connect(self._manage_databases)
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        layout.addWidget(combo, 1)
        layout.addWidget(manage)
        return row

    def _field_with_help(self, field: QWidget, title: str, help_text: str) -> QWidget:
        button = QToolButton()
        button.setText("?")
        button.setToolTip(help_text)
        button.setAutoRaise(False)
        button.setFixedSize(22, 22)
        button.setStyleSheet(
            "QToolButton {"
            "  border: 1px solid palette(mid);"
            "  border-radius: 10px;"
            "  background: palette(base);"
            "  font-weight: 600;"
            "  padding: 0px;"
            "}"
            "QToolButton:hover { background: palette(alternate-base); }"
        )
        button.clicked.connect(
            lambda checked=False, t=title, text=help_text: QMessageBox.information(
                self, t, text
            )
        )
        row = QWidget()
        layout = QHBoxLayout(row)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        layout.addWidget(field)
        layout.addWidget(button)
        return row

    def _update_postprocessing_visibility(self, tab_index: int) -> None:
        self._post_box.setVisible(tab_index == 0)
        if self.isVisible():
            self.adjustSize()

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

        dialog = BlastDatabaseManagerDialog(self._project, parent=self)
        exec_modal_dialog(dialog)
        self._populate_database_combos()

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

        self._use_gff = QCheckBox("Show GFF features on dotplot (-g)")
        self._use_gff.setChecked(opts.use_gff)
        self._use_gff.setToolTip(
            "Shade the dotplot with companion GFF features when available."
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

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)
        options_layout.addWidget(rc_box)
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

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _restore_defaults(self) -> None:
        defaults = FlexidotOptions()
        self._wordsize.setValue(defaults.wordsize)
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
        self._right_click_buttons[FLEXIDOT_RIGHT_CLICK_RUN].setChecked(True)

    def options(self) -> FlexidotOptions:
        fmt = self._output_format.currentData()
        right_click_action = FLEXIDOT_RIGHT_CLICK_RUN
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click_action = value
                break
        return FlexidotOptions(
            wordsize=self._wordsize.value(),
            lcs_shading=self._lcs_shading,
            use_gff=self._use_gff.isChecked(),
            output_format=normalize_flexidot_output_format(
                str(fmt or FLEXIDOT_OUTPUT_HTML)
            ),
            substitution_count=self._substitutions.value(),
            wobble_conversion=self._wobble.isChecked(),
            max_n_percentage=self._max_n.value(),
            include_reverse=self._include_reverse.isChecked(),
            only_vs_first=self._only_vs_first_value,
            line_width=self._line_width.value(),
            length_scaling=self._length_scaling.isChecked(),
            combined=self._combined.isChecked(),
            # Comparison mode is automatic: self for one sequence,
            # pairwise/all-vs-all for multiple.
            modes=DEFAULT_FLEXIDOT_MODES,
            right_click_action=right_click_action,
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
    """Batch TE annotation: Parameters select modules; Options holds settings/colors.

    Left-click the toolbar button opens this dialog. Setting / Color buttons stay
    on Parameters but jump to the Options tab, where the same content is expanded.
    """

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
        self.setWindowTitle("TE annotation")
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

        auto_header = QHBoxLayout()
        auto_title = _section_title_label("Auto annotation", auto_section)
        auto_header.addWidget(auto_title, 0, Qt.AlignmentFlag.AlignVCenter)
        auto_settings_btn = QPushButton("Setting", auto_section)
        auto_settings_btn.setToolTip(
            "Open the Options tab for LTR, TIR, tandem, and ORF detection settings"
        )
        auto_settings_btn.clicked.connect(self._show_auto_settings_on_options)
        auto_header.addWidget(auto_settings_btn, 0, Qt.AlignmentFlag.AlignVCenter)
        auto_color_btn = QPushButton("Color", auto_section)
        auto_color_btn.setToolTip(
            "Open the Options tab for feature-type and database colors"
        )
        auto_color_btn.clicked.connect(self._show_colors_on_options)
        auto_header.addWidget(auto_color_btn, 0, Qt.AlignmentFlag.AlignVCenter)
        auto_header.addStretch(1)
        auto_outer.addLayout(auto_header)

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

        sim_header = QHBoxLayout()
        sim_title = _section_title_label("Similarity annotation", sim_section)
        sim_header.addWidget(sim_title, 0, Qt.AlignmentFlag.AlignVCenter)
        sim_settings_btn = QPushButton("Setting", sim_section)
        sim_settings_btn.setToolTip(
            "Open the Options tab for similarity search parameters"
        )
        sim_settings_btn.clicked.connect(self._show_similarity_settings_on_options)
        sim_header.addWidget(sim_settings_btn, 0, Qt.AlignmentFlag.AlignVCenter)
        sim_header.addStretch(1)
        sim_outer.addLayout(sim_header)

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
            ("Run TE annotation with saved selection", TE_ANNOTATION_RIGHT_CLICK_RUN),
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
            QLabel("Right-clicking the TE annotation button does:", self)
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
        self._auto_settings_anchor = auto_box

        sim_box = QGroupBox("Similarity annotation settings", self)
        sim_box_layout = QVBoxLayout(sim_box)
        sim_box_layout.addWidget(self._similarity_settings_panel)
        self._similarity_settings_anchor = sim_box

        color_box = QGroupBox("Element colors", self)
        color_box_layout = QVBoxLayout(color_box)
        color_box_layout.addWidget(self._colors_panel)
        self._colors_anchor = color_box

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
        self._options_scroll = options_scroll

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

    def _show_options_section(self, widget: QWidget) -> None:
        self._tabs.setCurrentIndex(1)
        QTimer.singleShot(0, lambda: self._options_scroll.ensureWidgetVisible(widget))

    def _show_auto_settings_on_options(self) -> None:
        self._show_options_section(self._auto_settings_anchor)

    def _show_similarity_settings_on_options(self) -> None:
        self._show_options_section(self._similarity_settings_anchor)

    def _show_colors_on_options(self) -> None:
        self._show_options_section(self._colors_anchor)

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
    """TEAid parameters.

    Left-click the toolbar button opens this dialog. The Parameters tab holds
    BLAST settings and hit colors; the Options tab sets right-click behavior.
    """

    def __init__(
        self,
        *,
        initial: TeaAidOptions | None = None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("TEAid parameters")
        opts = initial or TeaAidOptions()
        self._for_run = for_run
        self._run_after = False

        high_pct, mid_pct = clamp_teaid_hit_len_pct(
            opts.hit_len_high_pct, opts.hit_len_mid_pct
        )

        self._evalue = QDoubleSpinBox()
        configure_report_evalue_spinbox(self._evalue)
        self._evalue.setValue(
            normalize_blast_evalue(opts.e_value, default=DEFAULT_TEAID_EVALUE)
        )
        self._evalue.setFixedWidth(120)
        self._evalue.setToolTip(
            "Reporting E-value for the TEAid genome BLAST search only."
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
            "BLASTN word size (-word_size) for the TEAid self-BLAST dot plot only (4–25)."
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
            ("Run TEAid with saved parameters", TEAID_RIGHT_CLICK_RUN),
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
        rc_box_layout.addWidget(QLabel("Right-clicking the TEAid button does:", self))
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
        params_form = QFormLayout()
        params_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        params_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint
        )
        params_form.setHorizontalSpacing(12)
        params_form.setVerticalSpacing(8)
        params_form.addRow("Genome BLAST E-value:", self._evalue)
        params_form.addRow("Self BLAST word size:", self._word_size)
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
        self._right_click_buttons[TEAID_RIGHT_CLICK_RUN].setChecked(True)

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
            x_axis_label=DEFAULT_TEAID_X_AXIS_LABEL,
            hit_len_high_pct=high_pct,
            hit_len_mid_pct=mid_pct,
            hit_color_high=self._hit_color_high or DEFAULT_TEAID_HIT_COLOR_HIGH,
            hit_color_mid=self._hit_color_mid or DEFAULT_TEAID_HIT_COLOR_MID,
            right_click_action=right_click_action,
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
            "QListWidget { padding-top: 4px; padding-bottom: 4px; }"
        )
        self._method_order.setToolTip(
            "Check steps to include in Combined MSA clean. Drag to set execution order."
        )
        for method_id in opts.method_order:
            item = QListWidgetItem(self._method_labels[method_id])
            item.setData(Qt.ItemDataRole.UserRole, method_id)
            item.setFlags(
                item.flags()
                | Qt.ItemFlag.ItemIsUserCheckable
                | Qt.ItemFlag.ItemIsEnabled
                | Qt.ItemFlag.ItemIsSelectable
                | Qt.ItemFlag.ItemIsDragEnabled
            )
            item.setCheckState(
                Qt.CheckState.Checked
                if enabled_methods.get(method_id, False)
                else Qt.CheckState.Unchecked
            )
            self._method_order.addItem(item)
        # Ensure every known method appears (in case order was incomplete).
        listed = {str(self._method_order.item(i).data(Qt.ItemDataRole.UserRole))
                  for i in range(self._method_order.count())}
        for method_id in DEFAULT_MSA_CLEAN_ORDER:
            if method_id in listed:
                continue
            item = QListWidgetItem(self._method_labels[method_id])
            item.setData(Qt.ItemDataRole.UserRole, method_id)
            item.setFlags(
                item.flags()
                | Qt.ItemFlag.ItemIsUserCheckable
                | Qt.ItemFlag.ItemIsEnabled
                | Qt.ItemFlag.ItemIsSelectable
                | Qt.ItemFlag.ItemIsDragEnabled
            )
            item.setCheckState(
                Qt.CheckState.Checked
                if enabled_methods.get(method_id, False)
                else Qt.CheckState.Unchecked
            )
            self._method_order.addItem(item)
        self._method_order.itemChanged.connect(self._on_method_item_changed)
        self._fit_method_list_height()

        combined_box = QGroupBox("Combined MSA clean", self)
        combined_layout = QVBoxLayout(combined_box)
        combined_layout.setContentsMargins(10, 8, 10, 8)
        combined_layout.setSpacing(6)
        hint = QLabel("Check steps to include. Drag to set order.", combined_box)
        hint.setWordWrap(True)
        combined_layout.addWidget(hint)
        combined_layout.addWidget(self._method_order)

        spin_width = 90

        div_box = QGroupBox("Clean rows with divergence")
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
        )

        gap_box = QGroupBox("Clean rows with gaps")
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
        )

        gappy_box = QGroupBox("Remove gappy columns")
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

    @staticmethod
    def _fill_labeled_grid(box, rows) -> None:
        """Lay out (label, widget) rows in a group box with a real label column.

        Uses an explicit QGridLayout rather than QFormLayout: the form's label
        column was collapsing to nothing here, clipping the descriptions to a
        stray character. A grid gives column 0 its natural width and stretches
        column 1, so every label is always fully visible.
        """
        grid = QGridLayout(box)
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
            item = QListWidgetItem(self._method_labels[method_id])
            item.setData(Qt.ItemDataRole.UserRole, method_id)
            item.setFlags(
                item.flags()
                | Qt.ItemFlag.ItemIsUserCheckable
                | Qt.ItemFlag.ItemIsEnabled
                | Qt.ItemFlag.ItemIsSelectable
                | Qt.ItemFlag.ItemIsDragEnabled
            )
            item.setCheckState(
                Qt.CheckState.Checked
                if enabled.get(method_id, False)
                else Qt.CheckState.Unchecked
            )
            self._method_order.addItem(item)
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
        file_filter: str = "TE hit tables (*.out *_blast.txt);;RepeatMasker (*.out);;All files (*)",
    ) -> None:
        super().__init__(parent)
        self._browse_title = browse_title
        self._file_filter = file_filter
        self._path = QLineEdit(self)
        self._path.setPlaceholderText(placeholder)
        self._path.setMinimumWidth(280)
        self._path.installEventFilter(self)
        if initial is not None:
            self._path.setText(str(initial.expanduser().resolve()))

        self._path.editingFinished.connect(self._emit_path_changed)

        browse_btn = QPushButton("Browse", self)
        browse_btn.clicked.connect(self._browse)

        remove_btn = QPushButton("Remove", self)
        remove_btn.setVisible(removable)
        remove_btn.clicked.connect(self.remove_requested.emit)
        self._remove_btn = remove_btn

        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        layout.addWidget(self._path, stretch=1)
        layout.addWidget(browse_btn)
        layout.addWidget(remove_btn)

    def eventFilter(self, obj, event) -> bool:  # noqa: N802
        if obj is self._path and event.type() == QEvent.Type.FocusIn:
            self.activated.emit(self)
        return super().eventFilter(obj, event)

    def set_removable(self, removable: bool) -> None:
        self._remove_btn.setVisible(removable)

    def set_path(self, path: Path) -> None:
        self._path.setText(str(path.expanduser().resolve()))
        self._emit_path_changed()

    def _emit_path_changed(self) -> None:
        self.path_changed.emit()

    def path(self) -> Path | None:
        text = self._path.text().strip()
        if not text:
            return None
        return Path(text).expanduser()

    def _browse(self) -> None:
        start = self._path.text().strip() or str(Path.home())
        picked, _ = QFileDialog.getOpenFileName(
            self,
            self._browse_title,
            start,
            self._file_filter,
        )
        if picked:
            self._path.setText(picked)
            self._emit_path_changed()


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
        file_filter: str = "TE hit tables (*.out *_blast.txt);;RepeatMasker (*.out);;All files (*)",
        add_button_label: str = "Add another TE file",
        add_button_tooltip: str = (
            "Add another RepeatMasker .out or BLAST hit table. "
            "Track order in the plot matches file order: first row = top track."
        ),
    ) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Minimum)
        self._placeholder = placeholder
        self._browse_title = browse_title
        self._file_filter = file_filter

        self._rows_layout = QVBoxLayout()
        self._rows_layout.setContentsMargins(0, 0, 0, 0)
        self._rows_layout.setSpacing(2)
        self._rows: list[_TeInputFileRow] = []
        self._active_row: _TeInputFileRow | None = None

        self._add_btn = QPushButton(add_button_label, self)
        self._add_btn.setToolTip(add_button_tooltip)
        self._add_btn.clicked.connect(self._on_add_clicked)

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
        removable = len(self._rows) > 1
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
        self._path = QLineEdit(self)
        self._path.setPlaceholderText("Optional GFF/GTF/BED/VCF annotation")
        self._path.setMinimumWidth(280)
        if initial is not None:
            self._path.setText(str(initial.expanduser().resolve()))
        browse_btn = QPushButton("Browse", self)
        browse_btn.clicked.connect(self._browse)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(6)
        layout.addWidget(self._path, stretch=1)
        layout.addWidget(browse_btn)

    def set_path(self, path: Path | None) -> None:
        if path is None:
            self._path.clear()
        else:
            self._path.setText(str(path.expanduser().resolve()))

    def path(self) -> Path | None:
        text = self._path.text().strip()
        if not text:
            return None
        return Path(text).expanduser()

    def _browse(self) -> None:
        start = self._path.text().strip() or str(Path.home())
        picked, _ = QFileDialog.getOpenFileName(
            self,
            "Select annotation file",
            start,
            "Annotations (*.gff *.gff3 *.gtf *.bed *.vcf *.vcf.gz);;All files (*)",
        )
        if picked:
            self._path.setText(picked)


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


def _browse_file_line(
    parent: QWidget,
    line: QLineEdit,
    *,
    title: str,
    filt: str,
) -> QWidget:
    row = QWidget(parent)
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

    btn = QPushButton("Browse", row)
    btn.clicked.connect(browse)
    h.addWidget(btn)
    return row


def _browse_dir_line(parent: QWidget, line: QLineEdit, *, title: str) -> QWidget:
    row = QWidget(parent)
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

    btn = QPushButton("Browse", row)
    btn.clicked.connect(browse)
    h.addWidget(btn)
    return row


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
        self._reference = QLineEdit(opts.reference_out)
        self._less_class = QCheckBox("Simpler TE classes (merge LINE/SINE to Non-LTR)")
        self._less_class.setChecked(opts.less_classification)
        self._debug = QCheckBox("Keep intermediate files (debug)")
        self._debug.setChecked(opts.debug)

        form = QFormLayout()
        form.addRow(
            "TEtrimmer folder:",
            _browse_dir_line(self, self._tetrimmer_home, title="Select TEtrimmer repository root"),
        )
        form.addRow(
            "Reference RM .out:",
            _browse_file_line(
                self,
                self._reference,
                title="Reference RepeatMasker .out",
                filt="RepeatMasker out (*.out);;All files (*)",
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
        self._standard = QLineEdit(opts.standard_fasta)
        self._rm_out = QLineEdit(opts.rm_out)
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
            ),
        )
        form.addRow(
            "Standard library FASTA:",
            _browse_file_line(
                self,
                self._standard,
                title="Standard TE library FASTA",
                filt="FASTA (*.fa *.fasta *.fas);;All files (*)",
            ),
        )
        form.addRow(
            "Existing RM .out (optional):",
            _browse_file_line(
                self,
                self._rm_out,
                title="Existing lib-vs-lib RepeatMasker .out",
                filt="RepeatMasker out (*.out);;All files (*)",
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


class RepeatQuantifyDialog(QDialog):
    """Choose a target (database) and a library FASTA, plus the engine.

    Both fields are pre-populated with the FASTA files already available (the
    current selection / folder) so they can be picked directly, and each has a
    Browse button for files elsewhere.
    """

    _FASTA_FILTER = "FASTA (*.fa *.fasta *.fna *.fas);;All files (*)"

    def __init__(self, *, candidates, start_dir: str = "", parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Repeat quantify")
        self._start_dir = start_dir
        cand = [Path(p) for p in candidates]

        self._genome_combo = QComboBox()
        self._library_combo = QComboBox()
        for combo in (self._genome_combo, self._library_combo):
            # Fixed width: compact, symmetric, and never collapses/jumps.
            combo.setFixedWidth(320)
            for path in cand:
                combo.addItem(path.name, str(path.expanduser().resolve()))
            combo.currentIndexChanged.connect(lambda _i, c=combo: self._sync_tooltip(c))
            self._sync_tooltip(combo)
        # Sensible default: library defaults to the 2nd candidate if distinct.
        if self._library_combo.count() > 1:
            self._library_combo.setCurrentIndex(1)

        genome_browse = QPushButton("Browse…")
        library_browse = QPushButton("Browse…")
        for browse in (genome_browse, library_browse):
            browse.setFixedWidth(88)  # short + symmetric
        genome_browse.clicked.connect(lambda: self._browse(self._genome_combo))
        library_browse.clicked.connect(lambda: self._browse(self._library_combo))

        inputs_box = QGroupBox("Inputs", self)
        form = QFormLayout(inputs_box)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(10)
        form.addRow(
            "Library (monomers / satellites):",
            self._field_row(self._library_combo, library_browse),
        )
        form.addRow(
            "Genome (scaffolds / assembly):",
            self._field_row(self._genome_combo, genome_browse),
        )

        # Tabs are the engine choice; each tab holds that engine's parameters.
        self._tabs = QTabWidget(self)
        self._blast_a = self._make_blast_controls()
        self._rm_b = self._make_rm_controls()
        self._blast_c = self._make_blast_controls()
        self._rm_d = self._make_rm_controls()
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
        self._tabs.setCurrentIndex(2)  # default: Both

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.button(QDialogButtonBox.StandardButton.Ok).setText("Run")
        buttons.accepted.connect(self._try_accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(16, 16, 16, 16)
        layout.setSpacing(14)
        layout.addWidget(inputs_box)
        layout.addWidget(self._tabs)
        layout.addWidget(buttons)

    def _make_blast_controls(self) -> dict:
        word = QSpinBox(self)
        word.setRange(4, 24)
        word.setValue(11)
        word.setToolTip(
            "BLASTn -word_size. Lower is more sensitive but much slower on "
            "big/repetitive genomes (7 can hang; 11 is a good balance)."
        )
        evalue = QDoubleSpinBox(self)
        configure_report_evalue_spinbox(evalue)
        evalue.setValue(1e-5)
        minid = QDoubleSpinBox(self)
        minid.setRange(0.0, 100.0)
        minid.setDecimals(0)
        minid.setSuffix(" %")
        minid.setValue(0)
        minid.setToolTip("BLASTn -perc_identity. 0 = no minimum.")
        mincov = QDoubleSpinBox(self)
        mincov.setRange(0.0, 100.0)
        mincov.setDecimals(0)
        mincov.setSuffix(" %")
        mincov.setValue(0)
        mincov.setToolTip(
            "BLASTn -qcov_hsp_perc: minimum query (monomer) coverage per hit. "
            "0 = no minimum."
        )
        for w in (word, evalue, minid, mincov):
            w.setFixedWidth(150)
        return {"word": word, "evalue": evalue, "minid": minid, "mincov": mincov}

    def _make_rm_controls(self) -> dict:
        div = QSpinBox(self)
        div.setRange(1, 100)
        div.setValue(40)
        div.setSuffix(" %")
        div.setFixedWidth(150)
        div.setToolTip(
            "RepeatMasker -div: report only copies less than this % diverged "
            "(higher captures older copies, but noisier)."
        )
        sens = QComboBox(self)
        sens.addItem("Sensitive (-s)", "sensitive")
        sens.addItem("Default", "default")
        sens.addItem("Quick (-q)", "quick")
        sens.setFixedWidth(150)
        sens.setToolTip(
            "RepeatMasker search sensitivity. Sensitive is slowest; Quick is "
            "fastest on large genomes."
        )
        return {"div": div, "sens": sens}

    @staticmethod
    def _blast_rows(c: dict) -> list:
        return [
            ("Word size:", c["word"]),
            ("E-value:", c["evalue"]),
            ("Min identity:", c["minid"]),
            ("Min coverage:", c["mincov"]),
        ]

    @staticmethod
    def _rm_rows(c: dict) -> list:
        return [
            ("Divergence cutoff (-div):", c["div"]),
            ("Sensitivity:", c["sens"]),
        ]

    def _params_tab(self, rows: list) -> QWidget:
        """Build a tab from (label, widget) rows; a None widget = bold header."""
        tab = QWidget()
        form = QFormLayout(tab)
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(8)
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

    @staticmethod
    def _field_row(combo: QComboBox, browse: QPushButton) -> QWidget:
        row = QWidget()
        hbox = QHBoxLayout(row)
        hbox.setContentsMargins(0, 0, 0, 0)
        hbox.setSpacing(8)
        hbox.addWidget(combo)
        hbox.addWidget(browse)
        hbox.addStretch(1)
        return row

    @staticmethod
    def _sync_tooltip(combo: QComboBox) -> None:
        data = combo.currentData()
        combo.setToolTip(str(data) if data else "")

    def _browse(self, combo: QComboBox) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Select FASTA", self._start_dir, self._FASTA_FILTER
        )
        if not path:
            return
        resolved = str(Path(path).expanduser().resolve())
        index = combo.findData(resolved)
        if index < 0:
            combo.addItem(Path(path).name, resolved)
            index = combo.count() - 1
        combo.setCurrentIndex(index)

    def _try_accept(self) -> None:
        target = self.target()
        library = self.library()
        if target is None or library is None:
            QMessageBox.information(
                self, "Repeat quantify",
                "Choose both a target FASTA and a library FASTA "
                "(pick from the list or Browse).",
            )
            return
        if target == library:
            QMessageBox.warning(
                self, "Repeat quantify",
                "Target and library must be different files.",
            )
            return
        self.accept()

    def target(self) -> Path | None:
        data = self._genome_combo.currentData()
        return Path(data) if data else None

    def library(self) -> Path | None:
        data = self._library_combo.currentData()
        return Path(data) if data else None

    def engine(self) -> str:
        return ("blastn", "repeatmasker", "both")[self._tabs.currentIndex()]

    def params(self) -> dict:
        idx = self._tabs.currentIndex()
        blast = {0: self._blast_a, 2: self._blast_c}.get(idx)
        rm = {1: self._rm_b, 2: self._rm_d}.get(idx)
        return {
            "blast_word_size": blast["word"].value() if blast else 11,
            "blast_evalue": blast["evalue"].value() if blast else 1e-5,
            "blast_min_identity": blast["minid"].value() if blast else 0.0,
            "blast_min_coverage": blast["mincov"].value() if blast else 0.0,
            "rm_div": rm["div"].value() if rm else 40,
            "rm_sensitivity": rm["sens"].currentData() if rm else "sensitive",
        }
