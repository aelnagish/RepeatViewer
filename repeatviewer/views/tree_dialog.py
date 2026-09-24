"""Tree dialog (Tools → Tree): input check, Alignment / Tree / Support and
rooting / Options tabs, and the Cancel / Save / Run footer.

Follows the RepeatViewer tool-dialog rules: Cancel, Save, Run in that order;
an Options tab with the right-click action, the process window and "use from
PATH"; the input is the file selection from the file panel.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QButtonGroup,
    QCheckBox,
    QDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QListWidget,
    QListWidgetItem,
    QPushButton,
    QRadioButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import tree_builder as tb
from repeatviewer.job_options import (
    TREE_ALIGN_ALWAYS,
    TREE_ALIGN_AUTO,
    TREE_ALIGN_NEVER,
    TREE_ENGINE_FASTME,
    TREE_ENGINE_FASTTREE,
    TREE_ENGINE_IQTREE,
    TREE_MODEL_FINDER,
    TREE_RIGHT_CLICK_DIALOG,
    TREE_RIGHT_CLICK_RUN,
    TREE_ROOT_MIDPOINT,
    TREE_ROOT_NONE,
    TREE_ROOT_OUTGROUP,
    TreeOptions,
)
from repeatviewer.views.tool_dialogs import (
    QComboBox,
    QDoubleSpinBox,
    QSpinBox,
    dialog_should_run,
    make_cancel_save_run_bar,
    make_process_output_box,
)

_FIELD_WIDTH = 300


def _form(parent: QWidget) -> QFormLayout:
    form = QFormLayout(parent)
    form.setContentsMargins(16, 14, 16, 14)
    form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
    form.setHorizontalSpacing(12)
    form.setVerticalSpacing(10)
    return form


def _note(text: str, parent: QWidget) -> QLabel:
    label = QLabel(text, parent)
    label.setWordWrap(True)
    label.setStyleSheet("color: palette(mid);")
    return label


class TreeParamsDialog(QDialog):
    """Build a phylogenetic tree from the selected FASTA file(s)."""

    def __init__(
        self,
        *,
        initial: TreeOptions | None = None,
        files: list[Path] | None = None,
        title: str = "Tree",
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        opts = initial or TreeOptions()
        self._files = [Path(p) for p in (files or [])]
        self._outgroup_saved = list(opts.outgroup)

        # --- input check strip ------------------------------------------------
        self._check = QLabel(self)
        self._check.setWordWrap(True)
        self._check.setTextFormat(Qt.TextFormat.RichText)
        self._check.setStyleSheet(
            "QLabel { background: palette(alternate-base); border-radius: 6px; padding: 8px 10px; }"
        )

        # --- Alignment tab ----------------------------------------------------
        self._align_mode = QComboBox(self)
        for label, value in (
            ("Align if not already aligned", TREE_ALIGN_AUTO),
            ("Always realign", TREE_ALIGN_ALWAYS),
            ("Use the input as an alignment", TREE_ALIGN_NEVER),
        ):
            self._align_mode.addItem(label, value)
        self._set_combo(self._align_mode, opts.align_mode)
        self._align_mode.setToolTip(
            "Sequences of equal length are treated as an existing alignment.\n"
            "'Use the input as an alignment' refuses unequal lengths instead of padding them."
        )
        self._strategy = QComboBox(self)
        for label, value in tb.MAFFT_STRATEGY_LABELS:
            self._strategy.addItem(label, value)
        self._set_combo(self._strategy, opts.mafft_strategy)
        self._adjust_direction = QCheckBox("Adjust direction (reverse-complement copies "
                                           "from the other strand)", self)
        self._adjust_direction.setChecked(opts.adjust_direction)
        self._adjust_direction.setToolTip(
            "MAFFT --adjustdirectionaccurately. Keep on for TE and satellite copies taken "
            "from both strands; otherwise they look extremely divergent and misplace in the tree."
        )
        self._min_cov = QSpinBox(self)
        self._min_cov.setRange(0, 100)
        self._min_cov.setSuffix(" % of median length")
        self._min_cov.setValue(round(opts.min_coverage * 100))
        self._min_cov.setToolTip(
            "Exclude fragments shorter than this share of the median sequence length "
            "(e.g. 5'-truncated copies). 0 turns the filter off. Excluded sequences are "
            "listed in the report and saved to <input>_tree.excluded.fa."
        )
        self._min_cov.valueChanged.connect(self._update_check)
        self._gap_open = QDoubleSpinBox(self)
        self._gap_open.setRange(0.0, 20.0)
        self._gap_open.setDecimals(2)
        self._gap_open.setSingleStep(0.1)
        self._gap_open.setValue(opts.mafft_gap_open)
        self._gap_open.setToolTip("MAFFT --op (gap opening penalty). MAFFT default: 1.53.")
        self._offset = QDoubleSpinBox(self)
        self._offset.setRange(0.0, 10.0)
        self._offset.setDecimals(3)
        self._offset.setSingleStep(0.05)
        self._offset.setValue(opts.mafft_offset)
        self._offset.setToolTip("MAFFT --ep (offset, acts like a gap extension penalty). Default: 0.")
        for w in (self._align_mode, self._strategy, self._min_cov):
            w.setFixedWidth(_FIELD_WIDTH)
        for w in (self._gap_open, self._offset):
            w.setFixedWidth(120)

        align_tab = QWidget(self)
        af = _form(align_tab)
        af.addRow("Alignment:", self._align_mode)
        af.addRow("MAFFT strategy:", self._strategy)
        af.addRow("", self._adjust_direction)
        af.addRow("Minimum length:", self._min_cov)
        more = QGroupBox("More options", align_tab)
        mf = QFormLayout(more)
        mf.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        mf.addRow("--op (gap opening):", self._gap_open)
        mf.addRow("--ep (offset):", self._offset)
        more.setCheckable(True)
        more.setChecked(abs(opts.mafft_gap_open - 1.53) > 1e-9 or opts.mafft_offset != 0.0)
        self._more_form = mf
        more.toggled.connect(self._show_more)
        self._show_more(more.isChecked())
        af.addRow(more)
        self._align_mode.currentIndexChanged.connect(self._sync_enabled)

        # --- Tree tab ---------------------------------------------------------
        self._engine_group = QButtonGroup(self)
        self._engine_buttons: dict[str, QRadioButton] = {}
        engine_box = QWidget(self)
        eb = QVBoxLayout(engine_box)
        eb.setContentsMargins(0, 0, 0, 0)
        eb.setSpacing(4)
        for value in (TREE_ENGINE_IQTREE, TREE_ENGINE_FASTME, TREE_ENGINE_FASTTREE):
            label = tb.ENGINE_LABELS[value]
            if value == TREE_ENGINE_IQTREE:
                label += "  — recommended"
            radio = QRadioButton(label, self)
            radio.setChecked(opts.engine == value)
            self._engine_group.addButton(radio)
            self._engine_buttons[value] = radio
            eb.addWidget(radio)
            radio.toggled.connect(self._sync_enabled)
        if not any(b.isChecked() for b in self._engine_buttons.values()):
            self._engine_buttons[TREE_ENGINE_IQTREE].setChecked(True)
        self._engine_status = _note("", self)

        self._model = QComboBox(self)
        self._model.setEditable(True)
        self._model.addItem("ModelFinder, best model by BIC", TREE_MODEL_FINDER)
        for m in ("GTR+F+G4", "TN+F+G4", "HKY+F+G4", "K2P+G4", "LG+G4"):
            self._model.addItem(m, m)
        if opts.model.upper() == TREE_MODEL_FINDER:
            self._model.setCurrentIndex(0)
        else:
            self._model.setCurrentText(opts.model)
        self._model.setToolTip(
            "IQ-TREE -m. Keep ModelFinder (MFP) unless you must fix a model; the chosen "
            "model is written into the report and methods text."
        )
        self._runs = QSpinBox(self)
        self._runs.setRange(1, 100)
        self._runs.setValue(opts.runs)
        self._runs.setToolTip(
            "IQ-TREE --runs: independent tree searches; the best is kept and the report "
            "says how many found the same topology."
        )
        self._seed = QSpinBox(self)
        self._seed.setRange(1, 2_000_000_000)
        self._seed.setValue(opts.seed)
        self._seed.setToolTip("Random seed, recorded in the report so the run can be repeated exactly.")
        self._threads = QSpinBox(self)
        self._threads.setRange(1, 256)
        self._threads.setValue(max(1, opts.threads))
        for w in (self._model,):
            w.setFixedWidth(_FIELD_WIDTH)
        for w in (self._runs, self._seed, self._threads):
            w.setFixedWidth(120)

        tree_tab = QWidget(self)
        tf = _form(tree_tab)
        tf.addRow("Method:", engine_box)
        tf.addRow("", self._engine_status)
        tf.addRow("Substitution model:", self._model)
        tf.addRow("Independent searches:", self._runs)
        tf.addRow("Random seed:", self._seed)
        tf.addRow("Threads:", self._threads)

        # --- Support and rooting tab -------------------------------------------
        self._alrt = QSpinBox(self)
        self._alrt.setRange(0, 100_000)
        self._alrt.setSingleStep(1000)
        self._alrt.setValue(opts.alrt)
        self._alrt.setToolTip("IQ-TREE --alrt: SH-aLRT replicates (0 = off). Trust ≥ 80.")
        self._ufboot = QSpinBox(self)
        self._ufboot.setRange(0, 100_000)
        self._ufboot.setSingleStep(1000)
        self._ufboot.setValue(opts.ufboot)
        self._ufboot.setToolTip("IQ-TREE -B: ultrafast bootstrap replicates (0 = off, "
                                "minimum 1000). Trust ≥ 95 — not the classic 70.")
        self._bnni = QCheckBox("Reduce UFBoot overestimation (-bnni)", self)
        self._bnni.setChecked(opts.bnni)
        self._dist_boot = QSpinBox(self)
        self._dist_boot.setRange(0, 10_000)
        self._dist_boot.setSingleStep(100)
        self._dist_boot.setValue(opts.distance_bootstrap)
        self._dist_boot.setToolTip("FastME bootstrap replicates (0 = off). Trust ≥ 70.")
        for w in (self._alrt, self._ufboot, self._dist_boot):
            w.setFixedWidth(120)

        self._root_group = QButtonGroup(self)
        self._root_buttons: dict[str, QRadioButton] = {}
        root_box = QWidget(self)
        rb = QVBoxLayout(root_box)
        rb.setContentsMargins(0, 0, 0, 0)
        rb.setSpacing(4)
        for label, value in (
            ("Outgroup (recommended when one exists)", TREE_ROOT_OUTGROUP),
            ("Midpoint (assumes similar rates)", TREE_ROOT_MIDPOINT),
            ("None (unrooted)", TREE_ROOT_NONE),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.root_mode == value)
            self._root_group.addButton(radio)
            self._root_buttons[value] = radio
            rb.addWidget(radio)
            radio.toggled.connect(self._sync_enabled)
        if not any(b.isChecked() for b in self._root_buttons.values()):
            self._root_buttons[TREE_ROOT_NONE].setChecked(True)
        self._outgroup_filter = QLineEdit(self)
        self._outgroup_filter.setPlaceholderText("Filter sequence names")
        self._outgroup_filter.textChanged.connect(self._filter_outgroups)
        self._outgroup_list = QListWidget(self)
        self._outgroup_list.setMinimumHeight(120)
        self._outgroup_list.setFixedWidth(_FIELD_WIDTH)
        self._outgroup_filter.setFixedWidth(_FIELD_WIDTH)

        support_tab = QWidget(self)
        sf = _form(support_tab)
        sf.addRow("SH-aLRT replicates:", self._alrt)
        sf.addRow("UFBoot2 replicates:", self._ufboot)
        sf.addRow("", self._bnni)
        sf.addRow("Bootstrap (FastME):", self._dist_boot)
        sf.addRow("Root:", root_box)
        sf.addRow("Outgroup:", self._outgroup_filter)
        sf.addRow("", self._outgroup_list)

        # --- Options tab ---------------------------------------------------------
        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        rc_box = QGroupBox("Right-click action", self)
        rcl = QVBoxLayout(rc_box)
        rcl.setSpacing(6)
        rcl.addWidget(QLabel("Right-clicking the Tree button does:", self))
        for label, value in (
            ("Build the tree with saved parameters", TREE_RIGHT_CLICK_RUN),
            ("Open this dialog", TREE_RIGHT_CLICK_DIALOG),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
            rcl.addWidget(radio)
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[TREE_RIGHT_CLICK_RUN].setChecked(True)
        process_box, self._show_process_window = make_process_output_box(
            self, tool_name="the tree build", checked=opts.show_process_window
        )
        self._use_system_mafft = QCheckBox("Use MAFFT from the system PATH", self)
        self._use_system_mafft.setChecked(opts.use_system_mafft)
        self._use_system_mafft.setToolTip("Off: RepeatViewer's bundled MAFFT.")
        self._use_system_path = QCheckBox(
            "Use IQ-TREE / FastME / FastTree from the system PATH only", self
        )
        self._use_system_path.setChecked(opts.use_system_path)
        self._use_system_path.setToolTip(
            "Off: look in the conda env set in File → Setting → External tools first, "
            "then on PATH. On: PATH only."
        )
        self._use_system_path.toggled.connect(self._sync_enabled)
        exe_box = QGroupBox("Executables", self)
        exl = QVBoxLayout(exe_box)
        exl.addWidget(self._use_system_mafft)
        exl.addWidget(self._use_system_path)
        options_tab = QWidget(self)
        ol = QVBoxLayout(options_tab)
        ol.setContentsMargins(16, 14, 16, 14)
        ol.setSpacing(14)
        ol.addWidget(rc_box)
        ol.addWidget(process_box)
        ol.addWidget(exe_box)
        ol.addStretch(1)

        self._tabs = QTabWidget(self)
        self._tabs.addTab(align_tab, "Alignment")
        self._tabs.addTab(tree_tab, "Tree")
        self._tabs.addTab(support_tab, "Support and rooting")
        self._tabs.addTab(options_tab, "Options")
        self._tabs.setCurrentIndex(1)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(self, left_widgets=[restore])

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(12)
        layout.addWidget(self._check)
        layout.addWidget(self._tabs)
        layout.addWidget(self._button_box)

        self._summary = tb.input_summary(self._files, opts.min_coverage) if self._files else None
        self._fill_outgroups()
        self._update_check()
        self._sync_enabled()

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def _set_combo(combo, value) -> None:
        idx = combo.findData(value)
        if idx >= 0:
            combo.setCurrentIndex(idx)

    def _engine(self) -> str:
        for value, radio in self._engine_buttons.items():
            if radio.isChecked():
                return value
        return TREE_ENGINE_IQTREE

    def _root_mode(self) -> str:
        for value, radio in self._root_buttons.items():
            if radio.isChecked():
                return value
        return TREE_ROOT_NONE

    def _show_more(self, on: bool) -> None:
        """Gap penalties stay hidden unless asked for (MAFFT's defaults are tested)."""
        for w in (self._gap_open, self._offset):
            w.setVisible(on)
            label = self._more_form.labelForField(w)
            if label is not None:
                label.setVisible(on)

    def _fill_outgroups(self) -> None:
        self._outgroup_list.clear()
        names = self._summary.names if self._summary else list(self._outgroup_saved)
        wanted = set(self._outgroup_saved)
        for name in names:
            item = QListWidgetItem(name, self._outgroup_list)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked if name in wanted else Qt.CheckState.Unchecked)

    def _filter_outgroups(self, text: str) -> None:
        needle = text.strip().lower()
        for i in range(self._outgroup_list.count()):
            item = self._outgroup_list.item(i)
            item.setHidden(bool(needle) and needle not in item.text().lower())

    def _checked_outgroups(self) -> list[str]:
        out = []
        for i in range(self._outgroup_list.count()):
            item = self._outgroup_list.item(i)
            if item.checkState() == Qt.CheckState.Checked:
                out.append(item.text())
        return out

    def _update_check(self) -> None:
        s = self._summary
        if not self._files:
            self._check.setText("No FASTA file selected. Select one or more FASTA files in "
                                "the file panel; several files are combined into one tree.")
            return
        if s is None or s.error:
            self._check.setText(f"<b>Input</b>: could not read the selection ({(s.error if s else '')}).")
            return
        frag = tb.fragment_count(s.lengths, self._min_cov.value() / 100.0)
        names = ", ".join(p.name for p in self._files[:3]) + ("…" if len(self._files) > 3 else "")
        parts = [f"<b>{s.n_sequences}</b> sequences", "DNA" if s.dna else "protein"]
        if s.aligned:
            parts.append(f"aligned, <b>{s.length}</b> columns, "
                         f"<b>{s.informative_sites}</b> parsimony-informative")
        else:
            parts.append("unaligned — MAFFT will align them")
        if frag:
            parts.append(f"<span style='color:#b25d00'><b>{frag}</b> fragment(s) below "
                         f"{self._min_cov.value()}% will be excluded</span>")
        warn = ""
        if s.n_sequences < tb.MIN_TAXA:
            warn = f"<br><span style='color:#b00020'>A tree needs at least {tb.MIN_TAXA} sequences.</span>"
        elif s.aligned and s.informative_sites is not None and s.informative_sites < 10:
            warn = ("<br><span style='color:#b25d00'>Very few informative sites: expect "
                    "poorly supported branches.</span>")
        self._check.setText(f"<b>Input</b> · {names}<br>" + " · ".join(parts) + warn)

    def _sync_enabled(self, *_args) -> None:
        engine = self._engine()
        iq = engine == TREE_ENGINE_IQTREE
        for w in (self._model, self._runs, self._alrt, self._ufboot, self._bnni):
            w.setEnabled(iq)
        self._dist_boot.setEnabled(engine == TREE_ENGINE_FASTME)
        outgroup = self._root_mode() == TREE_ROOT_OUTGROUP
        self._outgroup_list.setEnabled(outgroup)
        self._outgroup_filter.setEnabled(outgroup)
        aligning = self._align_mode.currentData() != TREE_ALIGN_NEVER
        for w in (self._strategy, self._adjust_direction):
            w.setEnabled(aligning)
        exe = tb.engine_executable(engine, system_path_only=self._use_system_path.isChecked())
        if exe:
            self._engine_status.setText(f"Found: {exe}")
            self._engine_status.setStyleSheet("color: palette(mid);")
        else:
            self._engine_status.setText(
                "Not found. " + tb.engine_install_hint(engine).replace("\n", " ")
            )
            self._engine_status.setStyleSheet("color: #b00020;")

    def _restore_defaults(self) -> None:
        d = TreeOptions()
        self._set_combo(self._align_mode, d.align_mode)
        self._set_combo(self._strategy, d.mafft_strategy)
        self._adjust_direction.setChecked(d.adjust_direction)
        self._min_cov.setValue(round(d.min_coverage * 100))
        self._gap_open.setValue(d.mafft_gap_open)
        self._offset.setValue(d.mafft_offset)
        self._engine_buttons[d.engine].setChecked(True)
        self._model.setCurrentIndex(0)
        self._runs.setValue(d.runs)
        self._seed.setValue(d.seed)
        self._threads.setValue(max(1, d.threads))
        self._alrt.setValue(d.alrt)
        self._ufboot.setValue(d.ufboot)
        self._bnni.setChecked(d.bnni)
        self._dist_boot.setValue(d.distance_bootstrap)
        self._root_buttons[d.root_mode].setChecked(True)
        self._right_click_buttons[TREE_RIGHT_CLICK_RUN].setChecked(True)
        self._show_process_window.setChecked(d.show_process_window)
        self._use_system_mafft.setChecked(d.use_system_mafft)
        self._use_system_path.setChecked(d.use_system_path)
        self._sync_enabled()

    # ------------------------------------------------------------------ results

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def options(self) -> TreeOptions:
        model_idx = self._model.currentIndex()
        text = self._model.currentText().strip()
        if model_idx == 0 and text == self._model.itemText(0):
            model = TREE_MODEL_FINDER
        else:
            model = text or TREE_MODEL_FINDER
        right_click = TREE_RIGHT_CLICK_RUN
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click = value
        ufboot = self._ufboot.value()
        return TreeOptions(
            align_mode=str(self._align_mode.currentData()),
            mafft_strategy=str(self._strategy.currentData()),
            adjust_direction=self._adjust_direction.isChecked(),
            mafft_gap_open=self._gap_open.value(),
            mafft_offset=self._offset.value(),
            min_coverage=self._min_cov.value() / 100.0,
            engine=self._engine(),
            model=model,
            runs=self._runs.value(),
            seed=self._seed.value(),
            threads=self._threads.value(),
            ufboot=0 if ufboot == 0 else max(1000, ufboot),
            alrt=self._alrt.value(),
            bnni=self._bnni.isChecked(),
            distance_bootstrap=self._dist_boot.value(),
            root_mode=self._root_mode(),
            outgroup=self._checked_outgroups() if self._summary else list(self._outgroup_saved),
            right_click_action=right_click,
            show_process_window=self._show_process_window.isChecked(),
            use_system_mafft=self._use_system_mafft.isChecked(),
            use_system_path=self._use_system_path.isChecked(),
        )
