"""samtools parameters (Tools → samtools), BLAST-style Inputs / Settings / Options."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSize, QTimer, Qt
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QDialog,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QRadioButton,
    QSizePolicy,
    QStackedWidget,
    QTabBar,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.samtools_tools import (
    is_alignment_path,
    is_windows,
    missing_backend_message,
    pysam_available,
    resolve_samtools_executable,
)
from repeatviewer.job_options import (
    ALIGN_INPUT_DEFINE,
    ALIGN_INPUT_SELECTED,
    DIALOG_TAB_OPTIONS,
    DEFAULT_SAMTOOLS_CONVERT_FORMAT,
    SAMTOOLS_CONVERT_BAM,
    SAMTOOLS_CONVERT_CRAM,
    SAMTOOLS_CONVERT_SAM,
    SAMTOOLS_DIALOG_TABS,
    SAMTOOLS_OP_CONVERT,
    SAMTOOLS_OP_INDEX,
    SAMTOOLS_OP_SORT,
    SamtoolsOptions,
    _normalize_samtools_convert_format,
    opening_dialog_tab,
)
from repeatviewer.tool_citations import CITATION_SAMTOOLS
from repeatviewer.views.tool_dialogs import (
    TOOL_PARAMS_DIALOG_WIDTH,
    TOOL_PARAMS_FIELD_WIDTH,
    TeMultiInputFilesWidget,
    dialog_should_run,
    make_cancel_save_run_bar,
    make_citation_box,
    make_process_output_box,
    sync_form_field_column_spacers,
    tool_params_dialog_width,
    _form_section_separator,
)
from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox
from repeatviewer.widgets.help_button import field_with_help
from repeatviewer.widgets.spin_boxes import NoWheelSpinBox as QSpinBox

_ALIGNMENT_FILTER = (
    "Alignments (*.sam *.bam *.cram *.sam.gz);;"
    "SAM (*.sam *.sam.gz);;BAM (*.bam);;CRAM (*.cram);;All files (*)"
)

_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): Samtools uses the alignment files "
    "selected in the Files panel when you click Run.<br><br>"
    "<b>Define inputs</b>: choose different alignment files for this run only "
    "(click in Files/Folders or Browse). Paths are not remembered next time.<br><br>"
    "<b>Input format</b>: SAM, BAM, or CRAM "
    "(<code>.sam</code>, <code>.bam</code>, <code>.cram</code>, "
    "<code>.sam.gz</code>). Index needs a coordinate-sorted BAM/CRAM."
)

_SORT_HELP = (
    "Coordinate-sort alignments into a BAM file "
    "(<code>samtools sort</code> / pysam). Output is written beside the input as "
    "<code>*.sorted.bam</code>."
)

_INDEX_HELP = (
    "Build a BAM/CRAM index (<code>.bai</code> / <code>.crai</code>) beside the "
    "input. The file must already be coordinate-sorted."
)

_CONVERT_HELP = (
    "Rewrite format without sorting (<code>samtools view</code> / pysam). "
    "Input is auto-detected (SAM/BAM/CRAM). Choose a different BAM, SAM, or "
    "CRAM output — files already in that format are skipped. "
    "Lighter than Sort; BAM/CRAM from unsorted input usually still need Sort "
    "before Index. CRAM may need a matching reference in the header."
)

_CONVERT_OUTPUT_HELP = (
    "Output format for Convert. samtools rewrites among SAM, BAM, and CRAM "
    "(no coordinate sorting). Inputs that are already this format are skipped. "
    "Default is BAM."
)

_THREADS_HELP = (
    "Threads for Sort / Index / Convert (<code>-@</code> / pysam). "
    "The same value is shared across the three tabs."
)

_INDEX_AFTER_HELP = (
    "After sorting, also build a <code>.bai</code> so the light genome browser "
    "can open the BAM window-by-window."
)


class SamtoolsParamsDialog(QDialog):
    """Tools → samtools: Sort / Index / Convert (+ Options)."""

    _DIALOG_WIDTH = TOOL_PARAMS_DIALOG_WIDTH
    # Shared with BLAST / Mapping so Inputs radios line up with Settings fields.

    def __init__(
        self,
        *,
        initial: SamtoolsOptions | None = None,
        for_run: bool = True,
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Samtools parameters")
        self.setMinimumWidth(self._DIALOG_WIDTH)
        opts = initial or SamtoolsOptions()
        self._for_run = for_run
        self._compact_size: QSize | None = None
        self._saved_operation = (
            opts.operation
            if opts.operation
            in (SAMTOOLS_OP_SORT, SAMTOOLS_OP_INDEX, SAMTOOLS_OP_CONVERT)
            else SAMTOOLS_OP_SORT
        )

        self._input_group = QButtonGroup(self)
        self._input_selected = QRadioButton("Use selected sequences", self)
        self._input_define = QRadioButton("Define inputs", self)
        self._input_group.addButton(self._input_selected)
        self._input_group.addButton(self._input_define)
        if opts.input_source == ALIGN_INPUT_DEFINE:
            self._input_define.setChecked(True)
        else:
            self._input_selected.setChecked(True)

        self._define_inputs = TeMultiInputFilesWidget(
            self,
            initial=None,
            start_empty=False,
            placeholder="Click in Files/Folders or Browse…",
            browse_title="Choose a SAM/BAM/CRAM file",
            file_filter=_ALIGNMENT_FILTER,
            add_button_label="Add another input",
            add_button_tooltip=(
                "Add another SAM/BAM/CRAM path. Click a file in the Files panel "
                "to fill a row, or use Browse."
            ),
            path_mode="file",
        )
        self._input_selected.toggled.connect(self._sync_inputs)
        self._input_define.toggled.connect(self._sync_inputs)
        self._define_inputs.paths_changed.connect(lambda: self._sync_inputs())

        threads = max(1, opts.num_threads)
        self._threads_sort = self._make_threads_spin(threads)
        self._threads_index = self._make_threads_spin(threads)
        self._threads_convert = self._make_threads_spin(threads)
        self._threads_spins = (
            self._threads_sort,
            self._threads_index,
            self._threads_convert,
        )
        self._wire_synced_spins(self._threads_spins)

        self._index_after = QCheckBox("Index after sort", self)
        self._index_after.setChecked(bool(opts.index_after_sort))
        self._index_after.setToolTip(_INDEX_AFTER_HELP)

        self._convert_output = QComboBox(self)
        self._convert_output.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        self._convert_output.setToolTip(_CONVERT_OUTPUT_HELP)
        for label, value in (
            ("BAM", SAMTOOLS_CONVERT_BAM),
            ("SAM", SAMTOOLS_CONVERT_SAM),
            ("CRAM", SAMTOOLS_CONVERT_CRAM),
        ):
            self._convert_output.addItem(label, value)
        convert_fmt = _normalize_samtools_convert_format(opts.convert_format)
        fmt_index = self._convert_output.findData(convert_fmt)
        self._convert_output.setCurrentIndex(
            fmt_index if fmt_index >= 0 else 0
        )

        self._inputs_panel = QWidget(self)
        self._inputs_panel.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )
        inputs_form = QFormLayout(self._inputs_panel)
        inputs_form.setContentsMargins(16, 10, 16, 0)
        inputs_form.setVerticalSpacing(4)
        inputs_form.setHorizontalSpacing(12)
        inputs_form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        inputs_form.setFieldGrowthPolicy(
            QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow
        )
        self._inputs_form = inputs_form
        inputs_form.addRow(QLabel("<b>Inputs</b>", self._inputs_panel))
        radios = QWidget(self._inputs_panel)
        radios_lay = QVBoxLayout(radios)
        radios_lay.setContentsMargins(0, 0, 0, 0)
        radios_lay.setSpacing(2)
        radios_lay.addWidget(
            field_with_help(
                self._input_selected, "Samtools inputs", _INPUTS_HELP, stretch=True
            )
        )
        # One Inputs "?" only (covers Use selected + Define + format).
        radios_lay.addWidget(self._input_define)
        self._inputs_radios_label = QLabel("")
        self._inputs_define_label = QLabel("")
        inputs_form.addRow(self._inputs_radios_label, radios)
        inputs_form.addRow(self._inputs_define_label, self._define_inputs)
        inputs_form.addRow(_form_section_separator())

        self._settings_forms: list[QFormLayout] = []
        self._index_after.setToolTip(_INDEX_AFTER_HELP)
        self._sort_tab = self._build_settings_tab(
            tooltip=_SORT_HELP,
            extra_rows=[
                ("Threads:", self._threads_sort),
                ("", self._index_after),
            ],
        )
        self._index_tab = self._build_settings_tab(
            tooltip=_INDEX_HELP,
            note="Requires a coordinate-sorted BAM or CRAM.",
            extra_rows=[("Threads:", self._threads_index)],
        )
        convert_input = QLabel("Auto (SAM / BAM / CRAM)", self)
        convert_input.setObjectName("ToolHintLabel")
        convert_input.setToolTip(
            "Input format is detected from each file "
            "(.sam, .bam, .cram, .sam.gz)."
        )
        self._convert_tab = self._build_settings_tab(
            tooltip=_CONVERT_HELP,
            note="Written beside the input (no sorting).",
            extra_rows=[
                ("Threads:", self._threads_convert),
                ("Input:", convert_input),
                ("Output:", self._convert_output),
            ],
        )

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="Samtools",
                checked=bool(opts.show_process_window),
                close_when_finished=bool(opts.close_process_window_when_finished),
            )
        )

        options_tab = QWidget(self)
        options_layout = QVBoxLayout(options_tab)
        options_layout.setContentsMargins(16, 14, 16, 14)
        options_layout.setSpacing(14)
        options_layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        options_layout.addWidget(process_box)
        options_layout.addWidget(make_citation_box(self, CITATION_SAMTOOLS))
        self._options_tab = options_tab

        self._tabs = QTabBar(self)
        self._tabs.setExpanding(False)
        self._tabs.setDocumentMode(True)
        self._tabs.setDrawBase(False)
        self._tabs.addTab("Sort")
        self._tabs.addTab("Index")
        self._tabs.addTab("Convert")
        self._tabs.addTab("Options")
        self._options_tab_index = 3

        self._stack = QStackedWidget(self)
        self._stack.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._stack.addWidget(self._sort_tab)
        self._stack.addWidget(self._index_tab)
        self._stack.addWidget(self._convert_tab)
        self._stack.addWidget(self._options_tab)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(
            self, left_widgets=[restore], include_run=True
        )
        self._run_check = self._validate_for_run

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(0)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._tabs, 0)
        layout.addWidget(self._inputs_panel, 0)
        layout.addWidget(self._stack, 0)
        layout.addSpacing(8)
        layout.addWidget(self._button_box, 0)

        open_tab = opening_dialog_tab(
            opts.last_dialog_tab,
            allowed=SAMTOOLS_DIALOG_TABS,
            first_tab=opts.operation
            if opts.operation
            in (SAMTOOLS_OP_SORT, SAMTOOLS_OP_INDEX, SAMTOOLS_OP_CONVERT)
            else SAMTOOLS_OP_SORT,
        )
        tab_index = {
            SAMTOOLS_OP_SORT: 0,
            SAMTOOLS_OP_INDEX: 1,
            SAMTOOLS_OP_CONVERT: 2,
            DIALOG_TAB_OPTIONS: 3,
        }.get(open_tab, 0)
        self._tabs.setCurrentIndex(tab_index)
        self._stack.setCurrentIndex(tab_index)
        self._tabs.currentChanged.connect(self._on_tab_changed)
        self._on_tab_changed(tab_index)
        self._sync_inputs(resize=False)
        self._sync_field_column()
        width = tool_params_dialog_width(
            define_inputs=self.input_source() == ALIGN_INPUT_DEFINE
        )
        self.resize(width, 360)

    def prepare_for_show(self) -> None:
        self._sync_field_column()
        self._fit_to_content()

    def sizeHint(self) -> QSize:  # noqa: N802
        # Prefer the last fitted size so polish/adjustSize cannot re-expand
        # to the tallest stack page after Options → Index.
        if self._compact_size is not None:
            return self._compact_size
        return super().sizeHint()

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        if self._compact_size is not None:
            return QSize(self._compact_size.width(), 0)
        return super().minimumSizeHint()

    def _make_threads_spin(self, value: int) -> QSpinBox:
        spin = QSpinBox(self)
        spin.setRange(1, 128)
        spin.setValue(max(1, value))
        spin.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        spin.setToolTip(_THREADS_HELP)
        return spin

    @staticmethod
    def _wire_synced_spins(spins: tuple[QSpinBox, ...]) -> None:
        def _sync(source: QSpinBox) -> None:
            value = source.value()
            for spin in spins:
                if spin is source or spin.value() == value:
                    continue
                spin.blockSignals(True)
                spin.setValue(value)
                spin.blockSignals(False)

        for spin in spins:
            spin.valueChanged.connect(lambda _v, s=spin: _sync(s))

    def _build_settings_tab(
        self,
        *,
        tooltip: str = "",
        note: str = "",
        extra_rows: list[tuple[str, QWidget]] | None = None,
    ) -> QWidget:
        tab = QWidget(self)
        tab.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        root = QVBoxLayout(tab)
        # Same left/right inset as Inputs so Threads lines up with Define inputs.
        root.setContentsMargins(16, 0, 16, 10)
        root.setSpacing(4)
        root.setAlignment(Qt.AlignmentFlag.AlignTop)
        basic = QWidget(tab)
        basic.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        form = QFormLayout(basic)
        form.setContentsMargins(0, 0, 0, 0)
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)
        # BLAST-style: Settings heading only — no "?" after the title.
        settings_label = QLabel("<b>Settings</b>", basic)
        if tooltip:
            settings_label.setToolTip(tooltip)
        form.addRow(settings_label)
        for label, widget in extra_rows or []:
            form.addRow(label, widget)
        # Index / Convert: one-line hint under Threads (field column).
        if note:
            note_lab = QLabel(note, basic)
            note_lab.setWordWrap(False)
            note_lab.setObjectName("ToolHintLabel")
            note_lab.setToolTip(note)
            form.addRow("", note_lab)
        self._settings_forms.append(form)
        root.addWidget(basic)
        return tab

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._sync_field_column()

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def current_dialog_tab(self) -> str:
        idx = self._tabs.currentIndex()
        return {
            0: SAMTOOLS_OP_SORT,
            1: SAMTOOLS_OP_INDEX,
            2: SAMTOOLS_OP_CONVERT,
            3: DIALOG_TAB_OPTIONS,
        }.get(idx, SAMTOOLS_OP_SORT)

    def current_operation(self) -> str:
        tab = self.current_dialog_tab()
        if tab in (SAMTOOLS_OP_SORT, SAMTOOLS_OP_INDEX, SAMTOOLS_OP_CONVERT):
            self._saved_operation = tab
            return tab
        return self._saved_operation

    def input_source(self) -> str:
        if self._input_define.isChecked():
            return ALIGN_INPUT_DEFINE
        return ALIGN_INPUT_SELECTED

    def run_files(self, *, panel_files: list[Path] | None = None) -> list[Path]:
        if self.input_source() == ALIGN_INPUT_DEFINE:
            paths = [Path(p) for p in self._define_inputs.paths() if str(p).strip()]
        else:
            paths = list(panel_files or [])
        return [p for p in paths if p.is_file() and is_alignment_path(p)]

    def apply_path_from_panel(self, path: Path) -> bool:
        if self.input_source() != ALIGN_INPUT_DEFINE:
            return False
        if self._tabs.currentIndex() == self._options_tab_index:
            return False
        path = Path(path)
        if not is_alignment_path(path):
            return False
        return self._define_inputs.accept_path(path)

    def apply_folder_from_panel(self, folder: Path) -> bool:
        del folder
        return False

    def options(self) -> SamtoolsOptions:
        op = self.current_operation()
        convert_fmt = self._convert_output.currentData()
        return SamtoolsOptions(
            operation=op,
            num_threads=max(1, int(self._threads_sort.value())),
            index_after_sort=self._index_after.isChecked(),
            convert_format=_normalize_samtools_convert_format(
                convert_fmt or DEFAULT_SAMTOOLS_CONVERT_FORMAT
            ),
            input_source=self.input_source(),
            input_paths=[str(p) for p in self._define_inputs.paths()],
            last_dialog_tab=self.current_dialog_tab(),
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=self._close_process_window.isChecked(),
            # Backend chooses pysam (macOS/Linux) or PATH Win-NGS samtools automatically.
            use_system_samtools=True,
            samtools_path="",
        )

    def _restore_defaults(self) -> None:
        defaults = SamtoolsOptions()
        for spin in self._threads_spins:
            spin.blockSignals(True)
            spin.setValue(defaults.num_threads)
            spin.blockSignals(False)
        self._index_after.setChecked(defaults.index_after_sort)
        fmt_index = self._convert_output.findData(defaults.convert_format)
        self._convert_output.setCurrentIndex(fmt_index if fmt_index >= 0 else 0)
        self._show_process_window.setChecked(defaults.show_process_window)
        self._close_process_window.setChecked(
            defaults.close_process_window_when_finished
        )
        self._input_selected.setChecked(True)
        self._fit_to_content()

    def _on_tab_changed(self, index: int) -> None:
        self._stack.setCurrentIndex(index)
        on_options = index == self._options_tab_index
        self._inputs_panel.setVisible(not on_options)
        self._compact_size = None
        # Refit after visibility/layout settles so Options → Index shrinks.
        if self.isVisible():
            QTimer.singleShot(0, self._fit_to_content)
        else:
            self._fit_to_content()

    def _sync_inputs(self, *, resize: bool = True) -> None:
        define = self.input_source() == ALIGN_INPUT_DEFINE
        # Collapse the form row fully so height tracks content (BLAST-style).
        self._inputs_form.setRowVisible(self._define_inputs, define)
        self._inputs_define_label.setVisible(define)
        self._define_inputs.setVisible(define)
        if define:
            self._define_inputs.setMinimumHeight(0)
            self._define_inputs.setMaximumHeight(16777215)
        else:
            self._define_inputs.setFixedHeight(0)
        self._define_inputs.updateGeometry()
        self._inputs_panel.updateGeometry()
        self._sync_field_column()
        if resize and self.isVisible():
            self._compact_size = None
            QTimer.singleShot(0, self._fit_to_content)

    def _sync_field_column(self) -> None:
        """Indent Inputs radios / Define inputs to the Settings field column."""
        sync_form_field_column_spacers(
            self,
            [self._inputs_radios_label, self._inputs_define_label],
            list(getattr(self, "_settings_forms", [])),
        )

    def _validate_for_run(self, _dialog: QDialog) -> str | None:
        exe = resolve_samtools_executable(explicit_path=None, use_system=True)
        if is_windows() and not exe:
            return missing_backend_message(exe=exe)
        if not is_windows() and not pysam_available() and not exe:
            return missing_backend_message(exe=exe)
        return None

    def _dialog_chrome_height(self) -> int:
        margins = self.layout().contentsMargins() if self.layout() else None
        top = margins.top() if margins is not None else 18
        bottom = margins.bottom() if margins is not None else 16
        tabs_h = max(24, self._tabs.sizeHint().height())
        buttons_h = self._button_box.sizeHint().height()
        return top + bottom + tabs_h + buttons_h + 8

    def _inputs_height(self) -> int:
        if not self._inputs_panel.isVisible():
            return 0
        self._inputs_panel.ensurePolished()
        self._inputs_panel.updateGeometry()
        return max(
            0,
            self._inputs_panel.minimumSizeHint().height(),
            self._inputs_panel.sizeHint().height(),
        )

    def _prepare_stack_pages_for_measure(self) -> QWidget | None:
        """Collapse inactive stack pages so their sizeHint cannot leave a blank band."""
        page = self._stack.currentWidget()
        for i in range(self._stack.count()):
            widget = self._stack.widget(i)
            if widget is None:
                continue
            if widget is page:
                widget.setMinimumHeight(0)
                widget.setMaximumHeight(16777215)
                widget.setSizePolicy(
                    QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
                )
            else:
                widget.setFixedHeight(0)
                widget.setSizePolicy(
                    QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Ignored
                )
        return page

    def _stack_page_height(self) -> int:
        page = self._prepare_stack_pages_for_measure()
        if page is None:
            return 0
        page.ensurePolished()
        root = page.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        page.adjustSize()
        page.updateGeometry()

        if page is self._options_tab:
            content_h = max(page.minimumSizeHint().height(), page.sizeHint().height())
            if root is not None:
                content_h = max(
                    content_h,
                    root.totalMinimumSize().height(),
                    root.totalSizeHint().height(),
                )
            return max(0, content_h)

        # Measure Settings form only — page sizeHint can keep leftover height
        # after visiting Options.
        margins = root.contentsMargins() if root is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 10
        basic = root.itemAt(0).widget() if root is not None and root.count() else None
        if basic is None:
            return max(0, page.sizeHint().height())
        basic.ensurePolished()
        basic_h = max(basic.minimumSizeHint().height(), basic.sizeHint().height())
        return max(0, margin_h + basic_h)

    def _pin_stack_height(self, height: int | None) -> None:
        page = self._stack.currentWidget()
        if height is None:
            self._stack.setMinimumHeight(0)
            self._stack.setMaximumHeight(16777215)
            if page is not None:
                page.setMinimumHeight(0)
                page.setMaximumHeight(16777215)
            return
        pinned = max(0, int(height))
        self._stack.setFixedHeight(pinned)
        if page is not None:
            page.setMinimumHeight(0)
            page.setMaximumHeight(pinned)

    def _fit_to_content(self) -> None:
        """Widen for Define inputs; shrink height to the visible tab content."""
        self._sync_field_column()
        width = tool_params_dialog_width(
            define_inputs=(
                self._inputs_panel.isVisible()
                and self.input_source() == ALIGN_INPUT_DEFINE
            )
        )
        self.setMinimumWidth(width)

        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000

        page_h = self._stack_page_height()
        inputs_h = self._inputs_height()
        height = self._dialog_chrome_height() + inputs_h + page_h
        height = min(max(height, 200), max_h)
        self._pin_stack_height(page_h)
        self.setMinimumHeight(0)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)
        self._compact_size = QSize(width, height)
