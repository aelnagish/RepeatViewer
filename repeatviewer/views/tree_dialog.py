"""Tree dialog (Tools → Tree), Align / BLAST–style function window.

Shared above the tabs: Mode (alignment) and Settings (monomers, root, skip
fragments). Tabs: IQ-TREE | FastME | FastTree | Options. Each engine tab has
Settings; IQ-TREE also has Advanced Settings (searches, seed, --bnni, --fast).
Save / Run take the engine from the open tab (on Options, from the last engine
tab). Input is used as given (DNA or protein); RT-domain trees are started from
Tools → Tree → RT-domain tree.

The dialog opens immediately; sequence lengths for Skip fragments / Outgroup are
loaded in the background so large FASTA files do not block the window.

Options tab: right-click engine, Executables (PATH checkboxes + grey status),
then Process output. REXdb is not configured here (RT-domain tree picks it when needed).
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import QThread, QTimer, Qt, Signal
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
    QSizePolicy,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import tree_builder as tb
from repeatviewer.backends import tree_portable_tools as tpt
from repeatviewer.job_options import (
    DIALOG_TAB_OPTIONS,
    TREE_ALIGN_ALWAYS,
    TREE_ALIGN_NEVER,
    TREE_COVERAGE_LONGEST,
    TREE_COVERAGE_MEDIAN,
    TREE_DIALOG_TABS,
    TREE_ENGINE_FASTME,
    TREE_ENGINE_FASTTREE,
    TREE_ENGINE_IQTREE,
    TREE_MODEL_FINDER,
    TREE_RIGHT_CLICK_RUN,
    TREE_ROOT_MIDPOINT,
    TREE_ROOT_NONE,
    TREE_ROOT_OUTGROUP,
    TREE_SEQ_AS_IS,
    TreeOptions,
)
from repeatviewer.views.tool_dialogs import (
    TOOL_PARAMS_DIALOG_WIDTH,
    TOOL_PARAMS_FIELD_WIDTH,
    QComboBox,
    QSpinBox,
    _form_section_separator,
    dialog_should_run,
    executable_checkbox_row,
    make_cancel_save_run_bar,
    make_executable_status_label,
    make_process_output_box,
    resolve_system_executable,
    set_path_checkbox_status,
    sync_form_field_column_spacers,
    system_path_executable_label,
)
from repeatviewer.widgets.help_button import field_with_help

ENGINE_TABS = (TREE_ENGINE_IQTREE, TREE_ENGINE_FASTME, TREE_ENGINE_FASTTREE)
_ENGINE_TITLES = {
    TREE_ENGINE_IQTREE: "IQ-TREE",
    TREE_ENGINE_FASTME: "FastME",
    TREE_ENGINE_FASTTREE: "FastTree",
}

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
IQ_SUPPORT_INFO = (
    "How much to trust each branch.<br>"
    "<b>SH-aLRT + UFBoot</b>: trust a branch at SH-aLRT ≥ 80 and UFBoot ≥ 95 "
    "(not the classic 70). UFBoot uses -B; SH-aLRT uses --alrt. "
    "UFBoot trees can be refined with --bnni (Advanced Settings)."
)
FASTME_SUPPORT_INFO = (
    "How much to trust each branch.<br>"
    "<b>Bootstrap</b>: trust ≥ 70. FastME is distance-based (TN93); fast, with a "
    "conservative bootstrap."
)
FASTTREE_SUPPORT_INFO = (
    "FastTree is approximate maximum likelihood, meant for very large sets "
    "(more than ~1,000 sequences). Its own local support, from 0 to 1, is always "
    "computed."
)
ROOT_INFO = (
    "<b>Unrooted</b>: no root; the tree shows relationships only.<br>"
    "<b>Midpoint</b>: root halfway along the longest path; assumes similar rates.<br>"
    "<b>Outgroup</b>: root on sequences known to lie outside the rest (best when one "
    "exists). Choose them with <b>Outgroup…</b> (IQ-TREE -o).<br>"
    "You can also reroot later in the tree viewer."
)
FRAG_INFO = (
    "Leave out fragments shorter than this share of the median length (or of the longest "
    "sequence), e.g. 5′-truncated TE copies; they would sit on long, wrong branches. "
    "0 = keep all. Left-out sequences are listed in the report."
)
SEARCH_INFO = (
    "<b>Searches</b> (IQ-TREE --runs equivalent): independent tree searches; the tree "
    "with the best likelihood is kept, and the report says how many found the same "
    "tree.<br>"
    "<b>Random seed</b> (--seed): recorded in the report so the run can be repeated "
    "exactly."
)
BNNI_INFO = (
    "IQ-TREE <b>--bnni</b>: optimize UFBoot trees by NNI on the bootstrap alignment. "
    "Reduces overestimation of support; recommended whenever UFBoot is used."
)
FAST_INFO = (
    "IQ-TREE <b>--fast</b>: faster tree search resembling FastTree. Use for a quick "
    "look; leave off for a publishable tree."
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
    "<b>ModelFinder</b> (-m MFP): tests models and selects the best by BIC.<br>"
    "<b>GTR</b>: all six substitution types at their own rate.<br>"
    "<b>TN (TN93)</b>: two transition rates (A↔G, C↔T), one transversion rate.<br>"
    "<b>HKY</b>: one transition and one transversion rate, unequal base frequencies.<br>"
    "<b>K2P (K80)</b>: as HKY, with equal base frequencies.<br>"
    "<b>TIM / JC</b>: further DNA models.<br>"
    "<b>LG / WAG / JTT</b>: protein models.<br>"
    "<b>+F</b> empirical frequencies · <b>+G4</b> rate variation among sites "
    "(4 gamma categories)."
)

IQ_SUPPORT = (  # data "alrt,ufboot" (Qt keeps strings, not tuples)
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


def _section(title: str, parent: QWidget | None = None) -> QLabel:
    return QLabel(f"<b>{title}</b>", parent)


def _form(parent: QWidget, *, margins: tuple[int, int, int, int]) -> QFormLayout:
    f = QFormLayout(parent)
    f.setContentsMargins(*margins)
    f.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    f.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
    f.setHorizontalSpacing(12)
    f.setVerticalSpacing(4)
    return f


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
            item.setCheckState(
                Qt.CheckState.Checked if name in wanted else Qt.CheckState.Unchecked
            )
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
        return [
            self._list.item(i).text()
            for i in range(self._list.count())
            if self._list.item(i).checkState() == Qt.CheckState.Checked
        ]


class _SummaryLoader(QThread):
    """Read FASTA lengths off the GUI thread so the dialog can open at once."""

    ready = Signal(object)

    def __init__(self, files: list[Path], min_coverage: float, parent=None) -> None:
        super().__init__(parent)
        self._files = list(files)
        self._min_coverage = float(min_coverage)

    def run(self) -> None:
        if not self._files:
            self.ready.emit(None)
            return
        try:
            self.ready.emit(tb.input_summary(self._files, self._min_coverage))
        except Exception as exc:  # noqa: BLE001
            self.ready.emit(
                tb.InputSummary(0, [], True, False, None, None, 0, error=str(exc))
            )


class TreeParamsDialog(QDialog):
    """Build a phylogenetic tree from the selected FASTA file(s)."""

    _DIALOG_WIDTH = TOOL_PARAMS_DIALOG_WIDTH

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
        self.setMinimumWidth(self._DIALOG_WIDTH)
        opts = initial or TreeOptions()
        self._initial = opts
        self._files = [Path(p) for p in (files or [])]
        self._outgroup = list(opts.outgroup)
        self._summary = None
        self._summary_loader: _SummaryLoader | None = None
        self._found_timer: QTimer | None = None
        self._iq_compact_size = None
        # REXdb is configured when running RT-domain tree, not on this form.
        self._rexdb_fasta = str(opts.rexdb_fasta or "")
        self._rexdb_references = int(opts.rexdb_references)

        self._build_shared(opts)
        self._build_engine_tabs(opts)
        self._build_options_tab(opts)

        self._tabs = QTabWidget(self)
        self._tabs.tabBar().setExpanding(False)
        self._tabs.tabBar().setDocumentMode(True)
        self._tabs.setDocumentMode(True)
        self._tabs.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._tab_index: dict[str, int] = {}
        for engine, page in zip(ENGINE_TABS, (self._iq_tab, self._fme_tab, self._ft_tab)):
            self._tab_index[engine] = self._tabs.addTab(page, _ENGINE_TITLES[engine])
        self._options_index = self._tabs.addTab(self._options_tab, "Options")
        self._last_engine = (
            opts.engine if opts.engine in ENGINE_TABS else TREE_ENGINE_IQTREE
        )
        open_tab = getattr(opts, "last_dialog_tab", self._last_engine)
        if open_tab == DIALOG_TAB_OPTIONS:
            self._tabs.setCurrentIndex(self._options_index)
        elif open_tab in self._tab_index:
            self._last_engine = open_tab
            self._tabs.setCurrentIndex(self._tab_index[open_tab])
        else:
            self._tabs.setCurrentIndex(self._tab_index[self._last_engine])
        self._tabs.currentChanged.connect(self._tab_changed)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._advanced_settings_btn = QPushButton("Advanced Settings ▸", self)
        self._advanced_settings_btn.setCheckable(True)
        self._advanced_settings_btn.setChecked(False)
        self._advanced_settings_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._advanced_settings_btn.toggled.connect(self._on_advanced_toggled)
        self._button_box = make_cancel_save_run_bar(
            self, left_widgets=[restore, self._advanced_settings_btn]
        )
        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(8)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._shared, 0)
        layout.addWidget(self._tabs, 0)
        layout.addWidget(self._button_box, 0)

        for w in (self._align_mode, self._root_mode, self._cov_ref):
            w.currentIndexChanged.connect(self._sync)
        self._min_cov.valueChanged.connect(self._update_check)
        self._tabs.currentChanged.connect(self._on_tab_for_found_and_fit)
        self._set_advanced_visible(False)
        self._update_advanced_settings_button()
        self._sync_field_columns()
        self._sync()
        QTimer.singleShot(0, self._start_summary_load)
        QTimer.singleShot(0, self._fit_to_content)
        self.resize(self._DIALOG_WIDTH, 420)

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._iq_compact_size = None
        self._sync_field_columns()
        QTimer.singleShot(0, self._fit_to_content)

    def _sync_field_columns(self) -> None:
        sync_form_field_column_spacers(
            self,
            [
                self._frag_hint_label,
                self._bnni_label,
                self._fast_label,
                self._ft_support_note_label,
            ],
            [
                self._top_form,
                self._iq_form,
                self._iq_adv_form,
                self._fme_form,
                self._ft_form,
            ],
        )

    # ------------------------------------------------------------------ building

    def _build_shared(self, opts: TreeOptions) -> None:
        self._align_mode = QComboBox(self)
        self._align_mode.addItem("Align with MAFFT", TREE_ALIGN_ALWAYS)
        self._align_mode.addItem("Already aligned, use as is", TREE_ALIGN_NEVER)
        self._align_mode.setCurrentIndex(1 if opts.align_mode == TREE_ALIGN_NEVER else 0)
        self._align_mode.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)

        self._rotate = QCheckBox("Rotate satellite monomers to a common start", self)
        self._rotate.setChecked(opts.rotate_monomers)

        self._root_mode = QComboBox(self)
        for label, value in (
            ("Unrooted", TREE_ROOT_NONE),
            ("Midpoint", TREE_ROOT_MIDPOINT),
            ("Outgroup", TREE_ROOT_OUTGROUP),
        ):
            self._root_mode.addItem(label, value)
        self._set_combo(self._root_mode, opts.root_mode)
        self._root_mode.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH - 136)
        self._outgroup_btn = QPushButton("Outgroup…", self)
        self._outgroup_btn.setFixedWidth(130)
        self._outgroup_btn.clicked.connect(self._choose_outgroup)

        self._min_cov = QSpinBox(self)
        self._min_cov.setRange(0, 100)
        self._min_cov.setSuffix(" %")
        self._min_cov.setValue(round(opts.min_coverage * 100))
        self._min_cov.setFixedWidth(90)
        self._cov_ref = QComboBox(self)
        self._cov_ref.addItem("of the median length", TREE_COVERAGE_MEDIAN)
        self._cov_ref.addItem("of the longest sequence", TREE_COVERAGE_LONGEST)
        self._set_combo(self._cov_ref, opts.coverage_reference)
        self._cov_ref.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH - 96)
        self._frag_hint = _note("Counting sequences…", self)
        self._frag_hint.setWordWrap(False)
        self._frag_hint_label = QLabel("")

        self._shared = QWidget(self)
        self._shared.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        f = _form(self._shared, margins=(16, 4, 16, 4))
        self._top_form = f
        f.addRow(_section("Mode", self._shared))
        f.addRow(
            "Alignment:",
            field_with_help(self._align_mode, "Alignment", ALIGN_INFO),
        )
        f.addRow(_form_section_separator())
        f.addRow(_section("Settings", self._shared))
        f.addRow(
            "Monomers:",
            field_with_help(self._rotate, "Rotate monomers", ROTATE_INFO),
        )
        f.addRow(
            "Root:",
            field_with_help(
                _row(self, self._root_mode, self._outgroup_btn, stretch=False),
                "Root",
                ROOT_INFO,
            ),
        )
        f.addRow(
            "Skip fragments:",
            field_with_help(
                _row(self, self._min_cov, self._cov_ref, stretch=False),
                "Skip fragments",
                FRAG_INFO,
            ),
        )
        f.addRow(self._frag_hint_label, self._frag_hint)
        self._set_row(f, self._frag_hint, bool(self._files))

    def _number_spin(self, lo: int, hi: int, value: int) -> QSpinBox:
        box = QSpinBox(self)
        box.setRange(lo, hi)
        box.setValue(value)
        box.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        return box

    def _build_engine_tabs(self, opts: TreeOptions) -> None:
        threads = max(1, opts.threads)
        self._support_iq = self._support_combo(
            IQ_SUPPORT,
            f"{opts.alrt},{opts.ufboot}",
            lambda v: "SH-aLRT {}, UFBoot {}".format(*v.split(",")),
        )
        self._model = QComboBox(self)
        self._model.setEditable(True)
        self._model.addItem("ModelFinder, best model by BIC", TREE_MODEL_FINDER)
        for m in (
            "GTR+F+G4",
            "TIM+F+G4",
            "TN+F+G4",
            "HKY+F+G4",
            "K2P+G4",
            "JC+G4",
            "LG+G4",
            "WAG+G4",
            "JTT+G4",
        ):
            self._model.addItem(m, m)
        if opts.model.upper() == TREE_MODEL_FINDER:
            self._model.setCurrentIndex(0)
        else:
            self._model.setCurrentText(opts.model)
        self._model.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        self._mask = QComboBox(self)
        for label, value in (("Off", ""), ("CG", "CG"), ("CG and CHG", "CG_CHG")):
            self._mask.addItem(label, value)
        self._set_combo(self._mask, opts.methylation_mask)
        self._mask.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        self._runs = self._number_spin(1, 100, opts.runs)
        self._seed = self._number_spin(1, 2_000_000_000, opts.seed)
        self._bnni = QCheckBox("Optimize UFBoot trees by NNI (--bnni)", self)
        self._bnni.setChecked(bool(opts.bnni))
        self._iqtree_fast = QCheckBox("Fast search (--fast)", self)
        self._iqtree_fast.setChecked(bool(getattr(opts, "iqtree_fast", False)))
        self._bnni_label = QLabel("")
        self._fast_label = QLabel("")
        self._threads = self._number_spin(1, 256, threads)
        self._threads.setToolTip(
            "Most threads to use. IQ-TREE picks the best number up to this (-T AUTO)."
        )
        self._support_fme = self._support_combo(
            FASTME_SUPPORT,
            opts.distance_bootstrap,
            lambda v: f"Bootstrap, {v} replicates",
        )
        self._threads_fme = self._number_spin(1, 256, threads)
        self._threads_fme.setToolTip("Most threads to use.")
        self._support_ft = QLabel("Local support (built in)", self)
        self._threads_ft = self._number_spin(1, 256, threads)
        self._threads_ft.setToolTip("Most threads to use.")
        self._thread_boxes = (self._threads, self._threads_fme, self._threads_ft)
        for box in self._thread_boxes:
            box.valueChanged.connect(self._threads_changed)

        # --- IQ-TREE tab ----------------------------------------------------
        self._iq_tab = QWidget(self)
        iq_root = QVBoxLayout(self._iq_tab)
        iq_root.setContentsMargins(16, 10, 16, 10)
        iq_root.setSpacing(4)
        iq_root.setAlignment(Qt.AlignmentFlag.AlignTop)

        iq_basic = QWidget(self._iq_tab)
        f = self._iq_form = _form(iq_basic, margins=(0, 0, 0, 0))
        f.addRow(_section("Settings", iq_basic))
        f.addRow(
            "Branch support:",
            field_with_help(self._support_iq, "Branch support", IQ_SUPPORT_INFO),
        )
        f.addRow(
            "Model:",
            field_with_help(self._model, "Substitution model", MODEL_INFO),
        )
        f.addRow(
            "Methylation check:",
            field_with_help(self._mask, "Methylation check", MASK_INFO),
        )
        f.addRow("Threads:", self._threads)
        iq_root.addWidget(iq_basic)

        self._iq_advanced_panel = QWidget(self._iq_tab)
        self._iq_advanced_panel.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        af = self._iq_adv_form = _form(self._iq_advanced_panel, margins=(0, 4, 0, 0))
        af.addRow(_form_section_separator())
        af.addRow(_section("Advanced Settings", self._iq_advanced_panel))
        af.addRow(
            "Searches:",
            field_with_help(self._runs, "Searches and seed", SEARCH_INFO),
        )
        af.addRow("Random seed:", self._seed)
        af.addRow(
            self._bnni_label,
            field_with_help(self._bnni, "UFBoot NNI", BNNI_INFO, stretch=True),
        )
        af.addRow(
            self._fast_label,
            field_with_help(self._iqtree_fast, "Fast search", FAST_INFO, stretch=True),
        )
        iq_root.addWidget(self._iq_advanced_panel)

        # --- FastME tab -----------------------------------------------------
        self._fme_tab = QWidget(self)
        fme_root = QVBoxLayout(self._fme_tab)
        fme_root.setContentsMargins(16, 10, 16, 10)
        fme_root.setSpacing(4)
        fme_root.setAlignment(Qt.AlignmentFlag.AlignTop)
        fme_basic = QWidget(self._fme_tab)
        f = self._fme_form = _form(fme_basic, margins=(0, 0, 0, 0))
        f.addRow(_section("Settings", fme_basic))
        f.addRow(
            "Branch support:",
            field_with_help(self._support_fme, "Branch support", FASTME_SUPPORT_INFO),
        )
        f.addRow("Threads:", self._threads_fme)
        fme_root.addWidget(fme_basic)

        # --- FastTree tab ---------------------------------------------------
        self._ft_tab = QWidget(self)
        ft_root = QVBoxLayout(self._ft_tab)
        ft_root.setContentsMargins(16, 10, 16, 10)
        ft_root.setSpacing(4)
        ft_root.setAlignment(Qt.AlignmentFlag.AlignTop)
        ft_basic = QWidget(self._ft_tab)
        f = self._ft_form = _form(ft_basic, margins=(0, 0, 0, 0))
        f.addRow(_section("Settings", ft_basic))
        f.addRow(
            "Branch support:",
            field_with_help(self._support_ft, "FastTree", FASTTREE_SUPPORT_INFO),
        )
        self._ft_support_note_label = QLabel("")
        f.addRow(
            self._ft_support_note_label,
            _note(
                "Local support values (0 to 1) are built in and always computed.",
                self._ft_tab,
            ),
        )
        f.addRow("Threads:", self._threads_ft)
        ft_root.addWidget(ft_basic)

        for page in (self._iq_tab, self._fme_tab, self._ft_tab):
            page.setSizePolicy(
                QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
            )

    def _build_options_tab(self, opts: TreeOptions) -> None:
        rc_box = QGroupBox("Right-click action", self)
        rc_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        rcl = QVBoxLayout(rc_box)
        rcl.setSpacing(6)
        rc_label = QLabel(
            "Right-clicking the Tree button builds the tree, with saved settings, using:",
            self,
        )
        rc_label.setWordWrap(True)
        rcl.addWidget(rc_label)
        self._right_click_group = QButtonGroup(self)
        self._right_click_buttons: dict[str, QRadioButton] = {}
        for value in ENGINE_TABS:
            radio = QRadioButton(_ENGINE_TITLES[value], self)
            radio.setChecked(opts.right_click_engine == value)
            self._right_click_group.addButton(radio)
            self._right_click_buttons[value] = radio
            rcl.addWidget(radio)
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[TREE_ENGINE_IQTREE].setChecked(True)

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="the tree build",
                checked=opts.show_process_window,
                close_when_finished=bool(
                    getattr(opts, "close_process_window_when_finished", False)
                ),
            )
        )

        self._system_boxes: dict[str, QCheckBox] = {}
        self._found: dict[str, QLabel] = {}
        exe_box = QGroupBox("Executables", self)
        exe_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        exl = QVBoxLayout(exe_box)
        exl.setContentsMargins(10, 8, 10, 8)
        exl.setSpacing(8)
        legacy = opts.use_system_path
        states = {
            TREE_ENGINE_IQTREE: opts.use_system_iqtree or legacy,
            TREE_ENGINE_FASTME: opts.use_system_fastme or legacy,
            TREE_ENGINE_FASTTREE: opts.use_system_fasttree or legacy,
            "mafft": opts.use_system_mafft,
        }
        names = {**_ENGINE_TITLES, "mafft": "MAFFT"}
        tooltips = {
            TREE_ENGINE_IQTREE: (
                f"Off: IQ-TREE {tpt.IQTREE_VERSION} bundled with RepeatViewer "
                "(macOS universal / Linux universal / Windows 64-bit). "
                "On: iqtree from PATH or the defined conda environment."
            ),
            TREE_ENGINE_FASTME: (
                "Off: FastME from the conda env in File → Settings → External tools, "
                "then PATH. On: PATH / conda only. FastME has no official multi-platform binary."
            ),
            TREE_ENGINE_FASTTREE: (
                "Off: FastTree from the conda env in File → Settings → External tools, "
                "then PATH. On: PATH / conda only. FastTree has Windows/Linux builds online but "
                "no macOS binary, so it is not bundled."
            ),
            "mafft": (
                "Off: the MAFFT bundled with RepeatViewer. "
                "On: MAFFT from PATH or the defined conda environment."
            ),
        }
        for key in (*ENGINE_TABS, "mafft"):
            box = QCheckBox(system_path_executable_label(names[key]), self)
            box.setChecked(bool(states[key]))
            box.setToolTip(tooltips[key])
            box.toggled.connect(self._schedule_found_update)
            status = make_executable_status_label(self)
            self._system_boxes[key] = box
            self._found[key] = status
            exl.addWidget(executable_checkbox_row(self, box, status))
        note = _note(
            f"IQ-TREE {tpt.IQTREE_VERSION} is included with RepeatViewer. "
            "FastME and FastTree need conda or your PATH.",
            self,
        )
        exl.addWidget(note)

        self._options_tab = QWidget(self)
        self._options_tab.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        ol = QVBoxLayout(self._options_tab)
        ol.setContentsMargins(16, 14, 16, 14)
        ol.setSpacing(14)
        ol.setAlignment(Qt.AlignmentFlag.AlignTop)
        ol.addWidget(rc_box)
        ol.addWidget(exe_box)
        ol.addWidget(process_box)
        ol.addStretch(1)

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
        if idx < 0:
            box.insertItem(0, custom_label(current), current)
            idx = 0
        box.setCurrentIndex(idx)
        box.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        return box

    def _engine(self) -> str:
        return self._last_engine

    def current_dialog_tab(self) -> str:
        """Active Tree window tab key (iqtree / fastme / fasttree / options)."""
        if self._tabs.currentIndex() == self._options_index:
            return DIALOG_TAB_OPTIONS
        return (
            self._last_engine
            if self._last_engine in TREE_DIALOG_TABS
            else TREE_ENGINE_IQTREE
        )

    def _tab_changed(self, index: int) -> None:
        for engine, i in self._tab_index.items():
            if i == index:
                self._last_engine = engine

    def _on_tab_for_found_and_fit(self, *_args) -> None:
        self._update_advanced_settings_button()
        self._schedule_found_update()
        if self.isVisible():
            QTimer.singleShot(0, self._fit_to_content)

    def _update_advanced_settings_button(self) -> None:
        on_iq = self._tabs.currentIndex() == self._tab_index[TREE_ENGINE_IQTREE]
        on_options = self._tabs.currentIndex() == self._options_index
        self._advanced_settings_btn.setVisible(on_iq and not on_options)
        if not on_iq and self._advanced_settings_btn.isChecked():
            self._advanced_settings_btn.blockSignals(True)
            self._advanced_settings_btn.setChecked(False)
            self._advanced_settings_btn.setText("Advanced Settings ▸")
            self._advanced_settings_btn.blockSignals(False)
            self._set_advanced_visible(False)

    def _set_advanced_visible(self, visible: bool) -> None:
        panel = self._iq_advanced_panel
        panel.setVisible(visible)
        if visible:
            panel.setMaximumHeight(16777215)
        else:
            panel.setMaximumHeight(0)
        panel.updateGeometry()
        root = self._iq_tab.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        self._iq_tab.updateGeometry()

    def _on_advanced_toggled(self, checked: bool) -> None:
        self._advanced_settings_btn.setText(
            "Advanced Settings ▾" if checked else "Advanced Settings ▸"
        )
        if checked:
            self._iq_compact_size = self.size()
        self._set_advanced_visible(checked)
        self._resize_for_iq_advanced(checked)

    def _threads_changed(self, value: int) -> None:
        for box in self._thread_boxes:
            if box.value() != value:
                box.blockSignals(True)
                box.setValue(value)
                box.blockSignals(False)

    def _protein(self) -> bool:
        s = self._summary
        return bool(s and not s.error and not s.dna)

    @staticmethod
    def _set_row(form: QFormLayout, field: QWidget, visible: bool) -> None:
        form.setRowVisible(field, visible)

    def _choose_outgroup(self) -> None:
        names = (
            self._summary.names
            if (self._summary and not self._summary.error)
            else list(self._outgroup)
        )
        dlg = OutgroupDialog(names, self._outgroup, self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self._outgroup = dlg.chosen()
            self._sync()

    def _start_summary_load(self) -> None:
        if self._summary_loader is not None:
            return
        if not self._files:
            self._frag_hint.setText("")
            self._set_row(self._top_form, self._frag_hint, False)
            return
        self._frag_hint.setText("Counting sequences…")
        self._set_row(self._top_form, self._frag_hint, True)
        loader = _SummaryLoader(self._files, self._min_cov.value() / 100.0, self)
        self._summary_loader = loader
        loader.ready.connect(self._on_summary_ready)
        loader.start()

    def _on_summary_ready(self, summary) -> None:
        self._summary_loader = None
        self._summary = summary
        self._sync()
        if self.isVisible():
            QTimer.singleShot(0, self._fit_to_content)

    def _update_check(self, *_args) -> None:
        s = self._summary
        show = bool(self._files) and s is not None and not s.error and bool(s.lengths)
        if show:
            frag = tb.fragment_count(
                s.lengths,
                self._min_cov.value() / 100.0,
                str(self._cov_ref.currentData()),
            )
            self._frag_hint.setText(
                f"{frag} of {s.n_sequences} sequences fall below this and will be left out."
                if frag
                else f"None of the {s.n_sequences} sequences fall below this."
            )
        elif self._files and s is None:
            self._frag_hint.setText("Counting sequences…")
            show = True
        elif self._files and s is not None and s.error:
            self._frag_hint.setText(s.error)
            show = True
        self._set_row(self._top_form, self._frag_hint, show)

    def _sync(self, *_args) -> None:
        protein = self._protein()
        aligning = self._align_mode.currentData() != TREE_ALIGN_NEVER
        self._set_row(
            self._top_form, self._rotate.parentWidget(), aligning and not protein
        )
        self._set_row(self._iq_form, self._mask.parentWidget(), not protein)
        outgroup = self._root_mode.currentData() == TREE_ROOT_OUTGROUP
        self._outgroup_btn.setEnabled(outgroup)
        n = len(self._outgroup)
        self._outgroup_btn.setText(f"Outgroup ({n})…" if n else "Outgroup…")
        self._outgroup_btn.setToolTip(
            ", ".join(self._outgroup)
            if n
            else "Choose the sequences that lie outside the rest."
        )
        self._update_check()
        self._schedule_found_update()

    def _schedule_found_update(self, *_args) -> None:
        if getattr(self, "_tabs", None) is None:
            return
        if self._tabs.currentIndex() != self._options_index:
            return
        if self._found_timer is None:
            self._found_timer = QTimer(self)
            self._found_timer.setSingleShot(True)
            self._found_timer.timeout.connect(self._update_found)
        self._found_timer.start(0)

    def _update_found(self) -> None:
        if getattr(self, "_tabs", None) is None:
            return
        if self._tabs.currentIndex() != self._options_index:
            return
        for engine in ENGINE_TABS:
            use_path = self._system_boxes[engine].isChecked()
            exe = (
                tb.engine_status_executable(engine, system_path_only=True)
                if use_path
                else None
            )
            set_path_checkbox_status(
                self._found[engine],
                str(exe) if exe else None,
                checked=use_path,
                has_builtin=(engine == TREE_ENGINE_IQTREE),
            )
            missing = use_path and exe is None and engine != TREE_ENGINE_IQTREE
            # IQ-TREE still has a bundled fallback when PATH is missing.
            if use_path and exe is None and engine == TREE_ENGINE_IQTREE:
                missing = False
            self._right_click_buttons[engine].setText(
                _ENGINE_TITLES[engine] + (" (not found)" if missing else "")
            )
        use_mafft = self._system_boxes["mafft"].isChecked()
        set_path_checkbox_status(
            self._found["mafft"],
            resolve_system_executable("mafft") if use_mafft else None,
            checked=use_mafft,
            has_builtin=True,
        )
        if self.isVisible():
            QTimer.singleShot(0, self._fit_to_content)

    def _dialog_chrome_height(self) -> int:
        lay = self.layout()
        margins = lay.contentsMargins() if lay is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 34
        spacing = lay.spacing() if lay is not None else 8
        shared_h = max(
            self._shared.minimumSizeHint().height(),
            self._shared.sizeHint().height(),
        )
        tab_bar_h = self._tabs.tabBar().sizeHint().height() + 8
        buttons_h = self._button_box.sizeHint().height()
        # shared + tabs + buttons → two spacings
        return margin_h + shared_h + tab_bar_h + buttons_h + (2 * spacing)

    def _iq_compact_content_height(self) -> int:
        root = self._iq_tab.layout()
        margins = root.contentsMargins() if root is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 20
        spacing = root.spacing() if root is not None else 4
        basic = self._iq_form.parentWidget()
        basic_h = max(basic.minimumSizeHint().height(), basic.sizeHint().height())
        return margin_h + basic_h + spacing

    def _resize_for_iq_advanced(self, expanded: bool) -> None:
        if self._tabs.currentIndex() != self._tab_index[TREE_ENGINE_IQTREE]:
            return
        if not expanded and self._iq_compact_size is not None:
            self.resize(self._iq_compact_size)
            self._iq_compact_size = None
            return
        self._fit_to_content()

    def _fit_to_content(self) -> None:
        """Shrink / grow to the active tab’s real height (no blank band)."""
        page = self._tabs.currentWidget()
        if page is None:
            return
        # Clear prior pins so newly shown PATH status lines can grow sizeHint.
        self._tabs.setMaximumHeight(16777215)
        page.setMaximumHeight(16777215)
        self._shared.ensurePolished()
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
        # When IQ Advanced is collapsed, measure basic only (panel maxHeight 0
        # can still leave a stale hint on some platforms).
        if (
            page is self._iq_tab
            and not self._iq_advanced_panel.isVisible()
        ):
            content_h = self._iq_compact_content_height()
        chrome = self._dialog_chrome_height()
        height = chrome + content_h
        width = max(self.width(), self._DIALOG_WIDTH)
        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        height = min(max(height, 220), max_h)
        self.setMinimumHeight(0)
        self._tabs.setMinimumHeight(0)
        page.setMinimumHeight(0)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)

    def done(self, result: int) -> None:  # noqa: N802
        loader = self._summary_loader
        if loader is not None and loader.isRunning():
            loader.ready.disconnect(self._on_summary_ready)
            loader.setParent(None)
            self._summary_loader = None
        super().done(result)

    # ------------------------------------------------------------------ defaults / results

    def _restore_defaults(self) -> None:
        d = TreeOptions()
        self._rexdb_fasta = str(d.rexdb_fasta or "")
        self._rexdb_references = int(d.rexdb_references)
        self._align_mode.setCurrentIndex(1 if d.align_mode == TREE_ALIGN_NEVER else 0)
        self._rotate.setChecked(d.rotate_monomers)
        self._set_combo(self._mask, d.methylation_mask)
        self._min_cov.setValue(round(d.min_coverage * 100))
        self._set_combo(self._cov_ref, d.coverage_reference)
        self._model.setCurrentIndex(0)
        self._runs.setValue(d.runs)
        self._seed.setValue(d.seed)
        self._bnni.setChecked(d.bnni)
        self._iqtree_fast.setChecked(d.iqtree_fast)
        self._threads.setValue(max(1, d.threads))
        self._set_combo(self._support_iq, f"{d.alrt},{d.ufboot}")
        self._set_combo(self._support_fme, d.distance_bootstrap)
        self._set_combo(self._root_mode, d.root_mode)
        self._right_click_buttons[d.right_click_engine].setChecked(True)
        self._show_process_window.setChecked(d.show_process_window)
        self._close_process_window.setChecked(d.close_process_window_when_finished)
        for box in self._system_boxes.values():
            box.setChecked(False)
        self._last_engine = d.engine
        self._tabs.setCurrentIndex(self._tab_index[d.engine])
        if self._advanced_settings_btn.isChecked():
            self._advanced_settings_btn.setChecked(False)
        self._sync()
        self._update_found()
        self._update_advanced_settings_button()
        QTimer.singleShot(0, self._fit_to_content)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def apply_path_from_panel(self, path: Path) -> bool:
        """Files-panel click on a REXdb FASTA still stores it for RT-domain runs."""
        from repeatviewer.backends.rt_domains import looks_like_rexdb

        resolved = Path(path).expanduser()
        if not looks_like_rexdb(resolved):
            return False
        self._rexdb_fasta = str(resolved)
        return True

    def options(self) -> TreeOptions:
        text = self._model.currentText().strip()
        if self._model.currentIndex() == 0 and text == self._model.itemText(0):
            model = TREE_MODEL_FINDER
        else:
            model = text or TREE_MODEL_FINDER
        right_click = TREE_ENGINE_IQTREE
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click = value
        alrt, ufboot = (int(x) for x in str(self._support_iq.currentData()).split(","))
        protein = self._protein()
        d = TreeOptions()
        engine = self._engine()
        threads = {
            TREE_ENGINE_IQTREE: self._threads,
            TREE_ENGINE_FASTME: self._threads_fme,
            TREE_ENGINE_FASTTREE: self._threads_ft,
        }[engine].value()
        sysbox = self._system_boxes
        return TreeOptions(
            align_mode=str(self._align_mode.currentData()),
            mafft_strategy=tb.MAFFT_CHOOSE,
            adjust_direction=True,
            rotate_monomers=self._rotate.isChecked(),
            sequence_mode=TREE_SEQ_AS_IS,
            rexdb_fasta=self._rexdb_fasta,
            rexdb_references=self._rexdb_references,
            methylation_mask="" if protein else str(self._mask.currentData() or ""),
            mafft_gap_open=d.mafft_gap_open,
            mafft_offset=d.mafft_offset,
            min_coverage=self._min_cov.value() / 100.0,
            coverage_reference=str(self._cov_ref.currentData()),
            engine=engine,
            last_dialog_tab=self.current_dialog_tab(),
            model=model,
            runs=self._runs.value(),
            seed=self._seed.value(),
            threads=max(1, threads),
            ufboot=0 if ufboot == 0 else max(1000, ufboot),
            alrt=alrt,
            bnni=self._bnni.isChecked(),
            iqtree_fast=self._iqtree_fast.isChecked(),
            distance_bootstrap=int(self._support_fme.currentData()),
            root_mode=str(self._root_mode.currentData()),
            outgroup=list(self._outgroup),
            right_click_action=TREE_RIGHT_CLICK_RUN,
            right_click_engine=right_click,
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=(
                self._close_process_window.isChecked()
            ),
            use_system_mafft=sysbox["mafft"].isChecked(),
            use_system_iqtree=sysbox[TREE_ENGINE_IQTREE].isChecked(),
            use_system_fastme=sysbox[TREE_ENGINE_FASTME].isChecked(),
            use_system_fasttree=sysbox[TREE_ENGINE_FASTTREE].isChecked(),
            use_system_path=False,
        )
