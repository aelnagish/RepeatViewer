"""Tree dialog (Tools → Tree): the choices that matter up front, the rest under
"More options", and the Options tab; Cancel / Save / Run footer.

Up front: Sequences (as they are / RT domain), Alignment (MAFFT / already
aligned), monomer rotation, Method, Branch support, Root. More options:
skip fragments, model, methylation check, searches, seed, threads. The MAFFT
method is picked automatically (L-INS-i up to 200 sequences, else FFT-NS-2)
and copies from the other strand are always turned round; special alignments
are made with the Align tool and used here as "Already aligned".

Follows the RepeatViewer tool-dialog rules: Cancel, Save, Run in that order;
an Options tab with the right-click action, the process window and "use from
PATH"; the input is the file selection from the file panel; help uses the
round "?" button of the other tool dialogs. Rows that do not apply are hidden.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QSize, Qt
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
    QStackedWidget,
    QTabWidget,
    QToolButton,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import tree_builder as tb
from repeatviewer.job_options import (
    TREE_ALIGN_ALWAYS,
    TREE_ALIGN_NEVER,
    TREE_COVERAGE_LONGEST,
    TREE_COVERAGE_MEDIAN,
    TREE_ENGINE_FASTME,
    TREE_ENGINE_FASTTREE,
    TREE_ENGINE_IQTREE,
    TREE_MODEL_FINDER,
    TREE_RIGHT_CLICK_DIALOG,
    TREE_RIGHT_CLICK_RUN,
    TREE_ROOT_MIDPOINT,
    TREE_ROOT_NONE,
    TREE_ROOT_OUTGROUP,
    TREE_SEQ_AS_IS,
    TREE_SEQ_RT,
    TreeOptions,
)
from repeatviewer.views.tool_dialogs import (
    QComboBox,
    QSpinBox,
    dialog_should_run,
    make_cancel_save_run_bar,
    make_process_output_box,
)
from repeatviewer.widgets.help_button import make_help_button
from repeatviewer.widgets.path_browse import PathBrowseRow

_FIELD_WIDTH = 300
_NUMBER_WIDTH = 110

ALIGN_INFO = (
    "<b>Align with MAFFT</b>: RepeatViewer aligns the sequences. It picks the method: "
    "L-INS-i (most accurate) up to 200 sequences, FFT-NS-2 (fast) above that. DNA copies "
    "from the other strand are turned round first.<br>"
    "<b>Already aligned, use as is</b>: the input is an alignment (all rows the same "
    "length). For a special alignment, make it with the Align tool and use it here."
)
ROTATE_INFO = (
    "Satellite monomers cut from an array start at arbitrary positions. This turns every "
    "monomer to the start and strand of the median-length monomer before aligning. "
    "Monomers too different to place are left as they are and listed in the report. "
    "DNA only."
)
METHOD_INFO = (
    "<b>Maximum likelihood (IQ-TREE)</b>: model-based; the most reliable tree, with "
    "SH-aLRT and UFBoot support. Slowest.<br>"
    "<b>Minimum evolution (FastME)</b>: distance-based (TN93); fast, with a "
    "conservative bootstrap.<br>"
    "<b>Approximate maximum likelihood (FastTree)</b>: for very large sets "
    "(more than ~1,000 sequences)."
)
SUPPORT_INFO = (
    "How much to trust each branch.<br>"
    "<b>SH-aLRT + UFBoot</b> (IQ-TREE): trust a branch at SH-aLRT ≥ 80 and UFBoot ≥ 95 "
    "(not the classic 70). UFBoot uses -bnni against overestimation.<br>"
    "<b>Bootstrap</b> (FastME): trust ≥ 70.<br>"
    "<b>FastTree</b>: its own local support from 0 to 1, always computed."
)
ROOT_INFO = (
    "<b>Unrooted</b>: no root; the tree shows relationships only.<br>"
    "<b>Midpoint</b>: root halfway along the longest path; assumes similar rates.<br>"
    "<b>Outgroup</b>: root on sequences known to lie outside the rest (best when one "
    "exists). Choose them with <b>Outgroup…</b>.<br>"
    "You can also reroot later in the tree viewer."
)
FRAG_INFO = (
    "Leave out fragments shorter than this share of the median length (or of the longest "
    "sequence), e.g. 5′-truncated TE copies; they would sit on long, wrong branches. "
    "0 = keep all. Left-out sequences are listed in the report."
)
SEARCH_INFO = (
    "<b>Searches</b> (IQ-TREE): independent tree searches; the tree with the best "
    "likelihood is kept, and the report says how many found the same tree.<br>"
    "<b>Random seed</b>: recorded in the report so the run can be repeated exactly."
)
MASK_INFO = (
    "<b>Methylation check</b>: builds a second tree without the alignment columns where any "
    "sequence has a C or G in a CG (or CHG) context, and compares the two trees in the report.<br>"
    "Methylated C mutates to T much more often than other bases change, so unrelated copies can "
    "share the same C→T by chance. If the well-supported groups stay the same without these "
    "sites, the tree does not depend on them.<br>"
    "<b>CG</b>: CG sites only. <b>CG and CHG</b>: also CHG sites (plants).<br>"
    "DNA only. Takes about twice as long."
)

MODEL_INFO = (
    "<b>ModelFinder</b>: tests the models and selects the best by BIC "
    "(Bayesian information criterion).<br>"
    "<b>GTR</b>: all six substitution types at their own rate.<br>"
    "<b>TN (TN93)</b>: two transition rates (A↔G, C↔T), one transversion rate.<br>"
    "<b>HKY</b>: one transition and one transversion rate, unequal base frequencies.<br>"
    "<b>K2P (K80)</b>: as HKY, with equal base frequencies.<br>"
    "<b>LG</b>: protein model.<br>"
    "<b>+F</b> empirical base frequencies · <b>+G4</b> rate variation among sites "
    "(4 gamma categories)."
)

SEQ_INFO = (
    "<b>As they are</b>: the sequences are aligned as given (DNA or protein).<br>"
    "<b>RT domain from DNA</b>: for LTR retrotransposons and LINEs. RepeatViewer finds the "
    "reverse transcriptase (RT) domain in the six reading frames of every copy and builds "
    "the tree from these protein pieces. The RT is the most conserved part of the element, "
    "so copies of different lineages can be compared. Reference RTs from REXdb are added, "
    "every copy gets its closest REXdb lineage, and the tips are coloured by lineage.<br>"
    "Copies without an RT (solo LTRs, deleted copies) are left out and listed in the report."
)

REXDB_INFO = (
    "REXdb: the reference database of plant retrotransposon protein domains "
    "(Neumann et al. 2019, Mobile DNA). It is not shipped with RepeatViewer; download it "
    "once from github.com/repeatexplorer/rexdb:<br>"
    "<b>Viridiplantae_v4.0.fasta</b> and <b>Viridiplantae_v4.0.classification</b>, "
    "saved in the same folder. Choose the .fasta here (or click it in Files).<br>"
    "The first run builds search profiles from it (about 10 s); later runs reuse them.<br>"
    "<b>References per lineage</b>: REXdb RTs added to the tree for each lineage of the "
    "RT types found, so every copy sits next to named lineages. 0 = none."
)


IQ_SUPPORT = (            # data "alrt,ufboot" (Qt keeps strings, not tuples)
    ("SH-aLRT + UFBoot, 1000 each", "1000,1000"),
    ("UFBoot only, 1000", "0,1000"),
    ("SH-aLRT only, 1000", "1000,0"),
    ("None", "0,0"),
)
FASTME_SUPPORT = (
    ("Bootstrap, 100 replicates", 100),
    ("Bootstrap, 1000 replicates", 1000),
    ("None", 0),
)


class _FitTabs(QTabWidget):
    """Tabs as tall as the tab on show (not the tallest tab), as wide as the widest."""

    def _fit(self, page_hint: QSize) -> QSize:
        full = super().sizeHint()
        biggest = QSize(0, 0)
        for i in range(self.count()):
            biggest = biggest.expandedTo(self.widget(i).sizeHint())
        return QSize(full.width(), full.height() - biggest.height() + page_hint.height())

    def sizeHint(self) -> QSize:
        page = self.currentWidget()
        return self._fit(page.sizeHint()) if page is not None else super().sizeHint()

    def minimumSizeHint(self) -> QSize:
        return self.sizeHint()


def _row(parent: QWidget, *widgets: QWidget, stretch: bool = True) -> QWidget:
    row = QWidget(parent)
    lay = QHBoxLayout(row)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(6)
    for w in widgets:
        lay.addWidget(w)
    if stretch:
        lay.addStretch(1)
    return row


def _note(text: str, parent: QWidget) -> QLabel:
    label = QLabel(text, parent)
    label.setWordWrap(True)
    label.setStyleSheet("color: palette(mid);")
    return label


def _help(parent: QWidget, title: str, text: str):
    return make_help_button(parent, title, text)


class OutgroupDialog(QDialog):
    """Pick outgroup sequences (filter + check list); Cancel / Save."""

    def __init__(self, names: list[str], chosen: list[str], parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Outgroup")
        self.resize(360, 420)
        self._filter = QLineEdit(self)
        self._filter.setPlaceholderText("Filter sequence names")
        self._filter.textChanged.connect(self._apply_filter)
        self._list = QListWidget(self)
        wanted = set(chosen)
        for name in names:
            item = QListWidgetItem(name, self._list)
            item.setFlags(item.flags() | Qt.ItemFlag.ItemIsUserCheckable)
            item.setCheckState(Qt.CheckState.Checked if name in wanted else Qt.CheckState.Unchecked)
        bar = make_cancel_save_run_bar(self, include_run=False)
        lay = QVBoxLayout(self)
        lay.addWidget(QLabel("Sequences known to lie outside the rest:", self))
        lay.addWidget(self._filter)
        lay.addWidget(self._list, 1)
        lay.addWidget(bar)

    def _apply_filter(self, text: str) -> None:
        needle = text.strip().lower()
        for i in range(self._list.count()):
            item = self._list.item(i)
            item.setHidden(bool(needle) and needle not in item.text().lower())

    def chosen(self) -> list[str]:
        return [self._list.item(i).text() for i in range(self._list.count())
                if self._list.item(i).checkState() == Qt.CheckState.Checked]


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
        self._initial = opts
        self._files = [Path(p) for p in (files or [])]
        self._outgroup = list(opts.outgroup)
        self._summary = tb.input_summary(self._files, opts.min_coverage) if self._files else None

        # --- up front -----------------------------------------------------------
        self._seq_mode = QComboBox(self)
        self._seq_mode.addItem("As they are", TREE_SEQ_AS_IS)
        self._seq_mode.addItem("RT domain from DNA (REXdb)", TREE_SEQ_RT)
        self._set_combo(self._seq_mode, opts.sequence_mode)
        self._rexdb = PathBrowseRow(
            self,
            initial=opts.rexdb_fasta or None,
            placeholder="REXdb protein FASTA: click it in Files, or Browse",
            browse_title="Choose the REXdb protein FASTA",
            file_filter="FASTA (*.fasta *.fa *.faa);;All files (*)",
        )
        self._rexdb.set_default_tooltip(
            "Viridiplantae_v4.0.fasta, with Viridiplantae_v4.0.classification next to it.")
        self._rexdb_refs = QSpinBox(self)
        self._rexdb_refs.setRange(0, 50)
        self._rexdb_refs.setValue(opts.rexdb_references)
        self._rexdb_refs.setSuffix(" per lineage")
        self._rexdb_refs.setFixedWidth(140)
        self._rexdb_status = _note("", self)
        self._rexdb_status.setWordWrap(False)

        self._align_mode = QComboBox(self)
        self._align_mode.addItem("Align with MAFFT", TREE_ALIGN_ALWAYS)
        self._align_mode.addItem("Already aligned, use as is", TREE_ALIGN_NEVER)
        self._align_mode.setCurrentIndex(1 if opts.align_mode == TREE_ALIGN_NEVER else 0)

        self._rotate = QCheckBox("Rotate satellite monomers to a common start", self)
        self._rotate.setChecked(opts.rotate_monomers)

        self._engine_group = QButtonGroup(self)
        self._engine_buttons: dict[str, QRadioButton] = {}
        engine_box = QWidget(self)
        eb = QVBoxLayout(engine_box)
        eb.setContentsMargins(0, 0, 0, 0)
        eb.setSpacing(6)
        for value in (TREE_ENGINE_IQTREE, TREE_ENGINE_FASTME, TREE_ENGINE_FASTTREE):
            radio = QRadioButton(tb.ENGINE_LABELS[value], self)
            radio.setChecked(opts.engine == value)
            self._engine_group.addButton(radio)
            self._engine_buttons[value] = radio
            eb.addWidget(radio)
        if not any(b.isChecked() for b in self._engine_buttons.values()):
            self._engine_buttons[TREE_ENGINE_IQTREE].setChecked(True)

        self._support_iq = self._support_combo(
            IQ_SUPPORT, f"{opts.alrt},{opts.ufboot}",
            lambda v: "SH-aLRT {}, UFBoot {}".format(*v.split(",")))
        self._support_fme = self._support_combo(FASTME_SUPPORT, opts.distance_bootstrap,
                                                lambda v: f"Bootstrap, {v} replicates")
        self._support_ft = QLabel("Local support (built in)", self)
        self._support_stack = QStackedWidget(self)
        for w in (self._support_iq, self._support_fme, self._support_ft):
            self._support_stack.addWidget(w)
        self._support_stack.setFixedWidth(_FIELD_WIDTH)
        self._support_stack.setFixedHeight(self._support_iq.sizeHint().height())

        self._root_mode = QComboBox(self)
        for label, value in (("Unrooted", TREE_ROOT_NONE), ("Midpoint", TREE_ROOT_MIDPOINT),
                             ("Outgroup", TREE_ROOT_OUTGROUP)):
            self._root_mode.addItem(label, value)
        self._set_combo(self._root_mode, opts.root_mode)
        self._outgroup_btn = QPushButton("Outgroup…", self)
        self._outgroup_btn.clicked.connect(self._choose_outgroup)

        # --- more options ---------------------------------------------------------
        self._min_cov = QSpinBox(self)
        self._min_cov.setRange(0, 100)
        self._min_cov.setSuffix(" %")
        self._min_cov.setValue(round(opts.min_coverage * 100))
        self._min_cov.setFixedWidth(90)
        self._cov_ref = QComboBox(self)
        self._cov_ref.addItem("of the median length", TREE_COVERAGE_MEDIAN)
        self._cov_ref.addItem("of the longest sequence", TREE_COVERAGE_LONGEST)
        self._set_combo(self._cov_ref, opts.coverage_reference)
        self._cov_ref.setFixedWidth(_FIELD_WIDTH - 96)
        self._frag_hint = _note("", self)
        self._frag_hint.setWordWrap(False)
        self._model = QComboBox(self)
        self._model.setEditable(True)
        self._model.addItem("ModelFinder, best model by BIC", TREE_MODEL_FINDER)
        for m in ("GTR+F+G4", "TN+F+G4", "HKY+F+G4", "K2P+G4", "LG+G4"):
            self._model.addItem(m, m)
        if opts.model.upper() == TREE_MODEL_FINDER:
            self._model.setCurrentIndex(0)
        else:
            self._model.setCurrentText(opts.model)
        self._mask = QComboBox(self)
        for label, value in (("Off", ""), ("CG", "CG"), ("CG and CHG", "CG_CHG")):
            self._mask.addItem(label, value)
        self._set_combo(self._mask, opts.methylation_mask)
        self._runs = QSpinBox(self)
        self._runs.setRange(1, 100)
        self._runs.setValue(opts.runs)
        self._seed = QSpinBox(self)
        self._seed.setRange(1, 2_000_000_000)
        self._seed.setValue(opts.seed)
        self._threads = QSpinBox(self)
        self._threads.setRange(1, 256)
        self._threads.setValue(max(1, opts.threads))
        self._threads.setToolTip("Most threads to use. IQ-TREE picks the best number up to this.")
        for w in (self._runs, self._seed, self._threads):
            w.setFixedWidth(_NUMBER_WIDTH)
        for w in (self._seq_mode, self._align_mode, self._model, self._mask):
            w.setFixedWidth(_FIELD_WIDTH)
        self._root_mode.setFixedWidth(_FIELD_WIDTH - 136)
        self._outgroup_btn.setFixedWidth(130)

        self._more_toggle = QToolButton(self)
        self._more_toggle.setText("More options")
        self._more_toggle.setCheckable(True)
        self._more_toggle.setChecked(False)
        self._more_toggle.setToolButtonStyle(Qt.ToolButtonStyle.ToolButtonTextBesideIcon)
        self._more_toggle.setAutoRaise(True)
        more_font = self._more_toggle.font()
        more_font.setPointSizeF(more_font.pointSizeF() + 1.0)
        more_font.setBold(True)
        self._more_toggle.setFont(more_font)
        self._more_toggle.setStyleSheet("QToolButton { border: none; }")
        self._more_toggle.toggled.connect(self._show_more)

        tree_tab = QWidget(self)
        f = QFormLayout(tree_tab)
        f.setContentsMargins(16, 14, 16, 14)
        f.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        f.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        f.setHorizontalSpacing(12)
        f.setVerticalSpacing(9)
        self._form = f
        f.addRow("Sequences:", _row(self, self._seq_mode, _help(self, "Sequences", SEQ_INFO)))
        self._rexdb_row = _row(self, self._rexdb, _help(self, "REXdb", REXDB_INFO), stretch=False)
        self._rexdb_row.setFixedWidth(_FIELD_WIDTH + 30)
        f.addRow("REXdb:", self._rexdb_row)
        f.addRow("References:", self._rexdb_refs)
        f.addRow("", self._rexdb_status)
        f.addRow("Alignment:", _row(self, self._align_mode, _help(self, "Alignment", ALIGN_INFO)))
        f.addRow("Monomers:", _row(self, self._rotate, _help(self, "Rotate monomers", ROTATE_INFO)))
        method_help = _help(self, "Method", METHOD_INFO)
        method_row = QWidget(self)
        mr = QHBoxLayout(method_row)
        mr.setContentsMargins(0, 0, 0, 0)
        mr.addWidget(engine_box)
        mr.addWidget(method_help, 0, Qt.AlignmentFlag.AlignTop)
        mr.addStretch(1)
        f.addRow("Method:", method_row)
        f.addRow("Branch support:", _row(self, self._support_stack, _help(self, "Branch support", SUPPORT_INFO)))
        f.addRow("Root:", _row(self, self._root_mode, self._outgroup_btn, _help(self, "Root", ROOT_INFO)))
        f.addRow(self._more_toggle)
        self._more_rows = [
            _row(self, self._min_cov, self._cov_ref, _help(self, "Skip fragments", FRAG_INFO)),
            self._frag_hint,
            _row(self, self._model, _help(self, "Substitution model", MODEL_INFO)),
            _row(self, self._mask, _help(self, "Methylation check", MASK_INFO)),
            _row(self, self._runs, _help(self, "Searches and seed", SEARCH_INFO)),
            self._seed,
            self._threads,
        ]
        for label, w in zip(("Skip fragments:", "", "Model:", "Methylation check:",
                             "Searches:", "Random seed:", "Threads:"), self._more_rows):
            f.addRow(label, w)

        # --- Options tab ------------------------------------------------------------
        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        rc_box = QGroupBox("Right-click on Tree", self)
        rcl = QVBoxLayout(rc_box)
        rcl.setSpacing(6)
        for label, value in (("Build the tree with saved settings", TREE_RIGHT_CLICK_RUN),
                             ("Open this window", TREE_RIGHT_CLICK_DIALOG)):
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
        self._use_system_path = QCheckBox("Use programs from the system PATH", self)
        self._use_system_path.setChecked(opts.use_system_path or opts.use_system_mafft)
        self._use_system_path.setToolTip(
            "Off: RepeatViewer's MAFFT, and IQ-TREE / FastME / FastTree from the conda env "
            "set in File → Settings → External tools (then PATH). On: all from PATH.")
        self._use_system_path.toggled.connect(self._sync)
        self._found = _note("", self)
        self._found.setWordWrap(False)
        exe_box = QGroupBox("Programs", self)
        exl = QVBoxLayout(exe_box)
        exl.addWidget(self._use_system_path)
        exl.addWidget(self._found)
        options_tab = QWidget(self)
        ol = QVBoxLayout(options_tab)
        ol.setContentsMargins(16, 14, 16, 14)
        ol.setSpacing(14)
        ol.addWidget(rc_box)
        ol.addWidget(process_box)
        ol.addWidget(exe_box)
        ol.addStretch(1)

        self._tabs = _FitTabs(self)
        self._tabs.addTab(tree_tab, "Tree")
        self._tabs.addTab(options_tab, "Options")
        # the window fits the tab on show, not the largest tab: no empty rows or columns
        self._tabs.currentChanged.connect(lambda _i: (self._tabs.updateGeometry(), self.adjustSize()))

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(self, left_widgets=[restore])
        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(12)
        layout.addWidget(self._tabs)
        layout.addWidget(self._button_box)
        # The window is always exactly as large as what it shows: no empty space,
        # and it shrinks back when More options closes.
        layout.setSizeConstraint(QVBoxLayout.SizeConstraint.SetFixedSize)

        # one label column for every row, shown or hidden, so the width never jumps
        labels = [f.itemAt(i, QFormLayout.ItemRole.LabelRole) for i in range(f.rowCount())]
        labels = [it.widget() for it in labels if it is not None and it.widget() is not None]
        width = max((lab.sizeHint().width() for lab in labels), default=0)
        for lab in labels:
            lab.setMinimumWidth(width)
            lab.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        for w in (self._seq_mode, self._align_mode, self._root_mode, self._cov_ref):
            w.currentIndexChanged.connect(self._sync)
        for radio in self._engine_buttons.values():
            radio.toggled.connect(self._sync)
        self._rexdb.path_changed.connect(self._sync)
        self._min_cov.valueChanged.connect(self._update_check)
        self._show_more(False)
        self._sync()

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def _set_combo(combo, value) -> None:
        idx = combo.findData(value)
        if idx >= 0:
            combo.setCurrentIndex(idx)

    def _support_combo(self, items, current, custom_label) -> QComboBox:
        box = QComboBox(self)
        for label, value in items:
            box.addItem(label, value)
        idx = box.findData(current)
        if idx < 0:                       # a setting made before: keep it
            box.insertItem(0, custom_label(current), current)
            idx = 0
        box.setCurrentIndex(idx)
        box.setFixedWidth(_FIELD_WIDTH)
        return box

    def _engine(self) -> str:
        for value, radio in self._engine_buttons.items():
            if radio.isChecked():
                return value
        return TREE_ENGINE_IQTREE

    def _protein(self) -> bool:
        s = self._summary
        return bool(s and not s.error and not s.dna)

    def _set_row(self, field: QWidget, visible: bool) -> None:
        self._form.setRowVisible(field, visible)

    def _show_more(self, on: bool) -> None:
        self._more_toggle.setArrowType(Qt.ArrowType.DownArrow if on else Qt.ArrowType.RightArrow)
        self._sync()

    def _choose_outgroup(self) -> None:
        names = self._summary.names if (self._summary and not self._summary.error) else list(self._outgroup)
        dlg = OutgroupDialog(names, self._outgroup, self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self._outgroup = dlg.chosen()
            self._sync()

    def _update_check(self, *_args) -> None:
        """One grey line under Skip fragments: how many sequences it leaves out."""
        s = self._summary
        show = (self._more_toggle.isChecked() and bool(self._files) and s is not None
                and not s.error and bool(s.lengths)
                and self._seq_mode.currentData() != TREE_SEQ_RT)
        if show:
            frag = tb.fragment_count(s.lengths, self._min_cov.value() / 100.0,
                                     str(self._cov_ref.currentData()))
            self._frag_hint.setText(
                f"{frag} of {s.n_sequences} sequences fall below this and will be left out."
                if frag else f"None of the {s.n_sequences} sequences fall below this.")
        self._set_row(self._frag_hint, show)

    def _sync(self, *_args) -> None:
        """Show only what applies to the current choices."""
        engine = self._engine()
        iq = engine == TREE_ENGINE_IQTREE
        protein = self._protein()
        self._seq_mode.setEnabled(not protein)
        rt = self._seq_mode.currentData() == TREE_SEQ_RT and not protein
        for w in (self._rexdb_row, self._rexdb_refs, self._rexdb_status):
            self._set_row(w, rt)
        if rt:
            from repeatviewer.backends.rt_domains import looks_like_rexdb

            path = self._rexdb.path()
            if path is None:
                msg, ok = "Choose REXdb to use this mode.", False
            elif looks_like_rexdb(path):
                msg, ok = "Found, with its .classification table.", True
            else:
                msg, ok = (f"{path.name}: not found, or {path.with_suffix('.classification').name} "
                           "is not next to it."), False
            self._rexdb_status.setText(msg)
            self._rexdb_status.setStyleSheet("color: palette(mid);" if ok else "color: #b00020;")
        aligning = self._align_mode.currentData() != TREE_ALIGN_NEVER
        self._set_row(self._rotate.parentWidget(), aligning and not protein and not rt)
        self._support_stack.setCurrentIndex(
            {TREE_ENGINE_IQTREE: 0, TREE_ENGINE_FASTME: 1}.get(engine, 2))
        outgroup = self._root_mode.currentData() == TREE_ROOT_OUTGROUP
        self._outgroup_btn.setEnabled(outgroup)
        n = len(self._outgroup)
        self._outgroup_btn.setText(f"Outgroup ({n})…" if n else "Outgroup…")
        self._outgroup_btn.setToolTip(", ".join(self._outgroup) if n else
                                      "Choose the sequences that lie outside the rest.")
        system = self._use_system_path.isChecked()
        for value, radio in self._engine_buttons.items():
            found = tb.engine_executable(value, system_path_only=system)
            radio.setText(tb.ENGINE_LABELS[value] + ("" if found else " (not found)"))
            radio.setToolTip(f"Found: {found}" if found else
                             "Not found. " + tb.engine_install_hint(value).replace("\n", " "))
        more = self._more_toggle.isChecked()
        visible = {
            0: more, 2: more and iq, 3: more and not protein and not rt, 4: more and iq,
            5: more, 6: more,
        }
        for i, w in enumerate(self._more_rows):
            if i in visible:
                self._set_row(w, visible[i])
        self._update_check()
        self._update_found()

    def _update_found(self) -> None:
        from repeatviewer.backends.sequence_aligner import mafft_command

        system = self._use_system_path.isChecked()
        lines, paths = [], []
        for engine, name in ((TREE_ENGINE_IQTREE, "IQ-TREE"), (TREE_ENGINE_FASTME, "FastME"),
                             (TREE_ENGINE_FASTTREE, "FastTree")):
            exe = tb.engine_executable(engine, system_path_only=system)
            lines.append(f"{name}: {'found' if exe else 'not found'}")
            paths.append(f"{name}: {exe or 'not found'}")
        cmd = mafft_command(use_system_mafft=system)
        lines.append(f"MAFFT: {'found' if cmd else 'not found'}")
        paths.append(f"MAFFT: {cmd[0][0] if cmd else 'not found'}")
        self._found.setText("\n".join(lines))
        self._found.setToolTip("\n".join(paths))

    def _restore_defaults(self) -> None:
        d = TreeOptions()
        self._set_combo(self._seq_mode, d.sequence_mode)
        self._rexdb_refs.setValue(d.rexdb_references)
        self._align_mode.setCurrentIndex(1 if d.align_mode == TREE_ALIGN_NEVER else 0)
        self._rotate.setChecked(d.rotate_monomers)
        self._set_combo(self._mask, d.methylation_mask)
        self._min_cov.setValue(round(d.min_coverage * 100))
        self._set_combo(self._cov_ref, d.coverage_reference)
        self._engine_buttons[d.engine].setChecked(True)
        self._model.setCurrentIndex(0)
        self._runs.setValue(d.runs)
        self._seed.setValue(d.seed)
        self._threads.setValue(max(1, d.threads))
        self._set_combo(self._support_iq, f"{d.alrt},{d.ufboot}")
        self._set_combo(self._support_fme, d.distance_bootstrap)
        self._set_combo(self._root_mode, d.root_mode)
        self._right_click_buttons[TREE_RIGHT_CLICK_RUN].setChecked(True)
        self._show_process_window.setChecked(d.show_process_window)
        self._use_system_path.setChecked(d.use_system_path or d.use_system_mafft)
        self._sync()

    # ------------------------------------------------------------------ results

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def apply_path_from_panel(self, path: Path) -> bool:
        """Files-panel click on a REXdb FASTA fills the REXdb field."""
        from repeatviewer.backends.rt_domains import looks_like_rexdb

        resolved = Path(path).expanduser()
        if not looks_like_rexdb(resolved):
            return False
        self._rexdb.set_path(resolved)
        return True

    def options(self) -> TreeOptions:
        text = self._model.currentText().strip()
        if self._model.currentIndex() == 0 and text == self._model.itemText(0):
            model = TREE_MODEL_FINDER
        else:
            model = text or TREE_MODEL_FINDER
        right_click = TREE_RIGHT_CLICK_RUN
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click = value
        alrt, ufboot = (int(x) for x in str(self._support_iq.currentData()).split(","))
        protein = self._protein()
        rt = self._seq_mode.currentData() == TREE_SEQ_RT and not protein
        system = self._use_system_path.isChecked()
        d = TreeOptions()
        return TreeOptions(
            align_mode=str(self._align_mode.currentData()),
            mafft_strategy=tb.MAFFT_CHOOSE,
            adjust_direction=True,
            rotate_monomers=self._rotate.isChecked(),
            sequence_mode=TREE_SEQ_AS_IS if protein else str(self._seq_mode.currentData()),
            rexdb_fasta=str(self._rexdb.path() or ""),
            rexdb_references=self._rexdb_refs.value(),
            methylation_mask="" if (protein or rt) else str(self._mask.currentData() or ""),
            mafft_gap_open=d.mafft_gap_open,
            mafft_offset=d.mafft_offset,
            min_coverage=self._min_cov.value() / 100.0,
            coverage_reference=str(self._cov_ref.currentData()),
            engine=self._engine(),
            model=model,
            runs=self._runs.value(),
            seed=self._seed.value(),
            threads=self._threads.value(),
            ufboot=0 if ufboot == 0 else max(1000, ufboot),
            alrt=alrt,
            bnni=True,
            distance_bootstrap=int(self._support_fme.currentData()),
            root_mode=str(self._root_mode.currentData()),
            outgroup=list(self._outgroup),
            right_click_action=right_click,
            show_process_window=self._show_process_window.isChecked(),
            use_system_mafft=system,
            use_system_path=system,
        )
