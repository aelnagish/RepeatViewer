"""Tree dialog (Tools → Tree), laid out like the BLAST parameters dialog.

Above the tabs: Alignment, monomer rotation, Root and Skip fragments. Tabs:
IQ-TREE | FastME | FastTree | Options. Save / Run take the engine from the open
tab (on Options, from the last engine tab). Input is used as given (DNA or protein);
RT-domain trees are started from Tools → Tree → RT-domain tree.

The dialog opens immediately; sequence lengths for Skip fragments / Outgroup are
loaded in the background so large FASTA files do not block the window.

Options tab: right-click engine, process window, per-program PATH toggles
(IQ-TREE and MAFFT are bundled; tick PATH to use your own), and optional REXdb
path for RT-domain trees.
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
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import tree_builder as tb
from repeatviewer.backends import tree_portable_tools as tpt
from repeatviewer.job_options import (
    TREE_ALIGN_ALWAYS,
    TREE_ALIGN_NEVER,
    TREE_COVERAGE_LONGEST,
    TREE_COVERAGE_MEDIAN,
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
    "(not the classic 70). UFBoot uses -bnni against overestimation."
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
    "exists). Choose them with <b>Outgroup…</b>.<br>"
    "You can also reroot later in the tree viewer."
)
FRAG_INFO = (
    "Leave out fragments shorter than this share of the median length (or of the longest "
    "sequence), e.g. 5′-truncated TE copies; they would sit on long, wrong branches. "
    "0 = keep all. Left-out sequences are listed in the report."
)
REXDB_INFO = (
    "Protein FASTA from REXdb (github.com/repeatexplorer/rexdb), with its "
    "<code>.classification</code> table in the same folder. Used only for "
    "<b>Tools → Tree → RT-domain tree</b>."
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


def _form(parent: QWidget, *, margins: tuple[int, int, int, int]) -> QFormLayout:
    f = QFormLayout(parent)
    f.setContentsMargins(*margins)
    f.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    f.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
    f.setHorizontalSpacing(12)
    f.setVerticalSpacing(9)
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
            self.ready.emit(tb.InputSummary(
                0, [], True, False, None, None, 0, error=str(exc)))


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
        self._summary = None
        self._summary_loader: _SummaryLoader | None = None
        self._found_timer: QTimer | None = None

        self._build_shared(opts)
        self._build_engine_tabs(opts)
        self._build_options_tab(opts)

        self._tabs = QTabWidget(self)
        self._tab_index: dict[str, int] = {}
        for engine, page in zip(ENGINE_TABS, (self._iq_tab, self._fme_tab, self._ft_tab)):
            self._tab_index[engine] = self._tabs.addTab(page, _ENGINE_TITLES[engine])
        self._options_index = self._tabs.addTab(self._options_tab, "Options")
        self._last_engine = opts.engine if opts.engine in ENGINE_TABS else TREE_ENGINE_IQTREE
        self._tabs.setCurrentIndex(self._tab_index[self._last_engine])
        self._tabs.currentChanged.connect(self._tab_changed)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(self, left_widgets=[restore])
        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(12)
        layout.addWidget(self._shared)
        layout.addWidget(self._tabs)
        layout.addWidget(self._button_box)

        # one label column across the shared rows and the engine tabs, so nothing jumps
        forms = (self._top_form, self._iq_form, self._fme_form, self._ft_form)
        labels = []
        for f in forms:
            for i in range(f.rowCount()):
                item = f.itemAt(i, QFormLayout.ItemRole.LabelRole)
                if item is not None and item.widget() is not None:
                    labels.append(item.widget())
        width = max((lab.sizeHint().width() for lab in labels), default=0)
        for lab in labels:
            lab.setMinimumWidth(width)
            lab.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)

        for w in (self._align_mode, self._root_mode, self._cov_ref):
            w.currentIndexChanged.connect(self._sync)
        self._min_cov.valueChanged.connect(self._update_check)
        self._tabs.currentChanged.connect(self._maybe_update_found)
        self._sync()
        # Large FASTAs: do not block opening — load lengths in the background.
        QTimer.singleShot(0, self._start_summary_load)

    # ------------------------------------------------------------------ building

    def _build_shared(self, opts: TreeOptions) -> None:
        self._align_mode = QComboBox(self)
        self._align_mode.addItem("Align with MAFFT", TREE_ALIGN_ALWAYS)
        self._align_mode.addItem("Already aligned, use as is", TREE_ALIGN_NEVER)
        self._align_mode.setCurrentIndex(1 if opts.align_mode == TREE_ALIGN_NEVER else 0)

        self._rotate = QCheckBox("Rotate satellite monomers to a common start", self)
        self._rotate.setChecked(opts.rotate_monomers)

        self._root_mode = QComboBox(self)
        for label, value in (("Unrooted", TREE_ROOT_NONE), ("Midpoint", TREE_ROOT_MIDPOINT),
                             ("Outgroup", TREE_ROOT_OUTGROUP)):
            self._root_mode.addItem(label, value)
        self._set_combo(self._root_mode, opts.root_mode)
        self._outgroup_btn = QPushButton("Outgroup…", self)
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
        self._cov_ref.setFixedWidth(_FIELD_WIDTH - 96)
        self._frag_hint = _note("Counting sequences…", self)
        self._frag_hint.setWordWrap(False)

        self._align_mode.setFixedWidth(_FIELD_WIDTH)
        self._root_mode.setFixedWidth(_FIELD_WIDTH - 136)
        self._outgroup_btn.setFixedWidth(130)

        self._shared = QWidget(self)
        f = _form(self._shared, margins=(16, 4, 16, 4))
        self._top_form = f
        f.addRow("Alignment:", _row(self, self._align_mode, _help(self, "Alignment", ALIGN_INFO)))
        f.addRow("Monomers:", _row(self, self._rotate, _help(self, "Rotate monomers", ROTATE_INFO)))
        f.addRow("Root:", _row(self, self._root_mode, self._outgroup_btn,
                               _help(self, "Root", ROOT_INFO)))
        f.addRow("Skip fragments:",
                 _row(self, self._min_cov, self._cov_ref, _help(self, "Skip fragments", FRAG_INFO)))
        f.addRow("", self._frag_hint)
        # Hide the hint until lengths are known (or there is nothing to count).
        self._set_row(f, self._frag_hint, bool(self._files))

    def _number_spin(self, lo: int, hi: int, value: int) -> QSpinBox:
        box = QSpinBox(self)
        box.setRange(lo, hi)
        box.setValue(value)
        box.setFixedWidth(_NUMBER_WIDTH)
        return box

    def _build_engine_tabs(self, opts: TreeOptions) -> None:
        threads = max(1, opts.threads)
        # IQ-TREE
        self._support_iq = self._support_combo(
            IQ_SUPPORT, f"{opts.alrt},{opts.ufboot}",
            lambda v: "SH-aLRT {}, UFBoot {}".format(*v.split(",")))
        self._model = QComboBox(self)
        self._model.setEditable(True)
        self._model.addItem("ModelFinder, best model by BIC", TREE_MODEL_FINDER)
        for m in ("GTR+F+G4", "TN+F+G4", "HKY+F+G4", "K2P+G4", "LG+G4"):
            self._model.addItem(m, m)
        if opts.model.upper() == TREE_MODEL_FINDER:
            self._model.setCurrentIndex(0)
        else:
            self._model.setCurrentText(opts.model)
        self._model.setFixedWidth(_FIELD_WIDTH)
        self._mask = QComboBox(self)
        for label, value in (("Off", ""), ("CG", "CG"), ("CG and CHG", "CG_CHG")):
            self._mask.addItem(label, value)
        self._set_combo(self._mask, opts.methylation_mask)
        self._mask.setFixedWidth(_FIELD_WIDTH)
        self._runs = self._number_spin(1, 100, opts.runs)
        self._seed = self._number_spin(1, 2_000_000_000, opts.seed)
        self._threads = self._number_spin(1, 256, threads)
        self._threads.setToolTip("Most threads to use. IQ-TREE picks the best number up to this.")
        # FastME
        self._support_fme = self._support_combo(FASTME_SUPPORT, opts.distance_bootstrap,
                                                lambda v: f"Bootstrap, {v} replicates")
        self._threads_fme = self._number_spin(1, 256, threads)
        self._threads_fme.setToolTip("Most threads to use.")
        # FastTree
        self._support_ft = QLabel("Local support (built in)", self)
        self._threads_ft = self._number_spin(1, 256, threads)
        self._threads_ft.setToolTip("Most threads to use.")
        self._thread_boxes = (self._threads, self._threads_fme, self._threads_ft)
        for box in self._thread_boxes:
            box.valueChanged.connect(self._threads_changed)

        margins = (16, 14, 16, 14)
        self._iq_tab = QWidget(self)
        f = self._iq_form = _form(self._iq_tab, margins=margins)
        f.addRow("Branch support:",
                 _row(self, self._support_iq, _help(self, "Branch support", IQ_SUPPORT_INFO)))
        f.addRow("Model:", _row(self, self._model, _help(self, "Substitution model", MODEL_INFO)))
        f.addRow("Methylation check:",
                 _row(self, self._mask, _help(self, "Methylation check", MASK_INFO)))
        f.addRow("Searches:", _row(self, self._runs, _help(self, "Searches and seed", SEARCH_INFO)))
        f.addRow("Random seed:", self._seed)
        f.addRow("Threads:", self._threads)

        self._fme_tab = QWidget(self)
        f = self._fme_form = _form(self._fme_tab, margins=margins)
        f.addRow("Branch support:",
                 _row(self, self._support_fme, _help(self, "Branch support", FASTME_SUPPORT_INFO)))
        f.addRow("Threads:", self._threads_fme)

        self._ft_tab = QWidget(self)
        f = self._ft_form = _form(self._ft_tab, margins=margins)
        f.addRow("Branch support:",
                 _row(self, self._support_ft, _help(self, "FastTree", FASTTREE_SUPPORT_INFO)))
        f.addRow("", _note("Local support values (0 to 1) are built in and always computed.",
                           self._ft_tab))
        f.addRow("Threads:", self._threads_ft)
        for page in (self._iq_tab, self._fme_tab, self._ft_tab):
            page.setMinimumWidth(_FIELD_WIDTH + 200)

    def _build_options_tab(self, opts: TreeOptions) -> None:
        # right-click program (like BLAST's right-click program)
        rc_box = QGroupBox("Right-click action", self)
        rcl = QVBoxLayout(rc_box)
        rcl.setSpacing(6)
        rc_label = QLabel(
            "Right-clicking the Tree button builds the tree, with saved settings, using:", self)
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
                self, tool_name="the tree build", checked=opts.show_process_window,
                close_when_finished=bool(
                    getattr(opts, "close_process_window_when_finished", False)
                ),
            )
        )

        # programs: bundled IQ-TREE / MAFFT, or PATH / conda for the rest
        self._system_boxes: dict[str, QCheckBox] = {}
        self._found: dict[str, QLabel] = {}
        exe_box = QGroupBox("Programs", self)
        exl = QVBoxLayout(exe_box)
        exl.setSpacing(6)
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
                "On: iqtree from your system PATH."
            ),
            TREE_ENGINE_FASTME: (
                "Off: FastME from the conda env in File → Settings → External tools, "
                "then PATH. On: PATH only. FastME has no official multi-platform binary."
            ),
            TREE_ENGINE_FASTTREE: (
                "Off: FastTree from the conda env in File → Settings → External tools, "
                "then PATH. On: PATH only. FastTree has Windows/Linux builds online but "
                "no macOS binary, so it is not bundled."
            ),
            "mafft": "Off: the MAFFT bundled with RepeatViewer. On: MAFFT from the system PATH.",
        }
        for key in (*ENGINE_TABS, "mafft"):
            box = QCheckBox(f"Use {names[key]} from the system PATH", self)
            box.setChecked(bool(states[key]))
            box.setToolTip(tooltips[key])
            box.toggled.connect(self._schedule_found_update)
            label = _note("", self)
            label.setWordWrap(False)
            self._system_boxes[key] = box
            self._found[key] = label
            exl.addWidget(_row(self, box, label))
        note = _note(
            f"IQ-TREE {tpt.IQTREE_VERSION} is included with RepeatViewer. "
            "FastME and FastTree need conda or your PATH.",
            self,
        )
        exl.addWidget(note)

        # REXdb: only for Tools → Tree → RT-domain tree (not shown on the main Tree form)
        rex_box = QGroupBox("REXdb (RT-domain tree)", self)
        rxf = _form(rex_box, margins=(10, 10, 10, 10))
        self._rexdb = PathBrowseRow(
            self,
            initial=opts.rexdb_fasta or None,
            placeholder="REXdb protein FASTA: click it in Files, or Browse",
            browse_title="Choose the REXdb protein FASTA",
            file_filter="FASTA (*.fasta *.fa *.faa);;All files (*)",
        )
        self._rexdb.set_default_tooltip(
            "Viridiplantae_v4.0.fasta, with Viridiplantae_v4.0.classification next to it.")
        self._rexdb.path_changed.connect(self._update_rexdb_status)
        self._rexdb_refs = QSpinBox(self)
        self._rexdb_refs.setRange(0, 50)
        self._rexdb_refs.setValue(opts.rexdb_references)
        self._rexdb_refs.setSuffix(" per lineage")
        self._rexdb_refs.setFixedWidth(140)
        self._rexdb_status = _note("", self)
        self._rexdb_status.setWordWrap(False)
        self._rexdb_row = _row(self, self._rexdb, _help(self, "REXdb", REXDB_INFO), stretch=False)
        rxf.addRow("REXdb:", self._rexdb_row)
        rxf.addRow("References:", self._rexdb_refs)
        rxf.addRow("", self._rexdb_status)

        self._options_tab = QWidget(self)
        ol = QVBoxLayout(self._options_tab)
        ol.setContentsMargins(16, 14, 16, 14)
        ol.setSpacing(14)
        ol.addWidget(rc_box)
        ol.addWidget(process_box)
        ol.addWidget(exe_box)
        ol.addWidget(rex_box)
        ol.addStretch(1)
        self._update_rexdb_status()

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
        """Engine of the open engine tab (on Options: the last engine tab used)."""
        return self._last_engine

    def _tab_changed(self, index: int) -> None:
        for engine, i in self._tab_index.items():
            if i == index:
                self._last_engine = engine

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
        names = self._summary.names if (self._summary and not self._summary.error) else list(self._outgroup)
        dlg = OutgroupDialog(names, self._outgroup, self)
        if dlg.exec() == QDialog.DialogCode.Accepted:
            self._outgroup = dlg.chosen()
            self._sync()

    def _start_summary_load(self) -> None:
        """Load FASTA lengths off the GUI thread so opening stays instant."""
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

    def _update_check(self, *_args) -> None:
        """One grey line under Skip fragments: how many sequences it leaves out."""
        s = self._summary
        show = bool(self._files) and s is not None and not s.error and bool(s.lengths)
        if show:
            frag = tb.fragment_count(s.lengths, self._min_cov.value() / 100.0,
                                     str(self._cov_ref.currentData()))
            self._frag_hint.setText(
                f"{frag} of {s.n_sequences} sequences fall below this and will be left out."
                if frag else f"None of the {s.n_sequences} sequences fall below this.")
        elif self._files and s is None:
            self._frag_hint.setText("Counting sequences…")
            show = True
        elif self._files and s is not None and s.error:
            self._frag_hint.setText(s.error)
            show = True
        self._set_row(self._top_form, self._frag_hint, show)

    def _update_rexdb_status(self, *_args) -> None:
        from repeatviewer.backends.rt_domains import looks_like_rexdb

        path = self._rexdb.path()
        if path is None:
            msg, ok = "Optional. Needed for Tools → Tree → RT-domain tree.", True
        elif looks_like_rexdb(path):
            msg, ok = "Found, with its .classification table.", True
        else:
            msg, ok = (f"{path.name}: not found, or {path.with_suffix('.classification').name} "
                       "is not next to it."), False
        self._rexdb_status.setText(msg)
        self._rexdb_status.setStyleSheet("color: palette(mid);" if ok else "color: #b00020;")

    def _sync(self, *_args) -> None:
        """Show only what applies to the current choices."""
        protein = self._protein()
        aligning = self._align_mode.currentData() != TREE_ALIGN_NEVER
        self._set_row(self._top_form, self._rotate.parentWidget(), aligning and not protein)
        self._set_row(self._iq_form, self._mask.parentWidget(), not protein)
        outgroup = self._root_mode.currentData() == TREE_ROOT_OUTGROUP
        self._outgroup_btn.setEnabled(outgroup)
        n = len(self._outgroup)
        self._outgroup_btn.setText(f"Outgroup ({n})…" if n else "Outgroup…")
        self._outgroup_btn.setToolTip(", ".join(self._outgroup) if n else
                                      "Choose the sequences that lie outside the rest.")
        self._update_check()
        self._schedule_found_update()

    def _schedule_found_update(self, *_args) -> None:
        """Update found labels after the Options tab has painted (avoids lag on click)."""
        if getattr(self, "_tabs", None) is None:
            return
        if self._tabs.currentIndex() != self._options_index:
            return
        if self._found_timer is None:
            self._found_timer = QTimer(self)
            self._found_timer.setSingleShot(True)
            self._found_timer.timeout.connect(self._update_found)
        self._found_timer.start(0)

    def _maybe_update_found(self, *_args) -> None:
        self._schedule_found_update()

    def _update_found(self) -> None:
        from repeatviewer.backends.sequence_aligner import mafft_command
        from repeatviewer.backends.tree_portable_tools import bundled_engine_executable

        if getattr(self, "_tabs", None) is None:
            return
        if self._tabs.currentIndex() != self._options_index:
            return
        for engine in ENGINE_TABS:
            use_path = self._system_boxes[engine].isChecked()
            if use_path:
                exe = tb.engine_status_executable(engine, system_path_only=True)
                label = "found" if exe else "not found"
                tip = str(exe) if exe else tb.engine_install_hint(engine)
                missing = exe is None
            else:
                exe = bundled_engine_executable(engine)
                if exe:
                    label, tip, missing = "bundled", str(exe), False
                else:
                    # FastME / FastTree: no vendored binary; runtime still tries conda then PATH.
                    # Do not probe conda here — that made Options lag on every open.
                    label, tip, missing = "conda / PATH", tb.engine_install_hint(engine), False
            self._found[engine].setText(label)
            self._found[engine].setToolTip(tip)
            self._right_click_buttons[engine].setText(
                _ENGINE_TITLES[engine] + (" (not found)" if missing else ""))
        cmd = mafft_command(use_system_mafft=self._system_boxes["mafft"].isChecked())
        if cmd and not self._system_boxes["mafft"].isChecked():
            self._found["mafft"].setText("bundled")
        else:
            self._found["mafft"].setText("found" if cmd else "not found")
        self._found["mafft"].setToolTip(cmd[0][0] if cmd else "MAFFT was not found.")

    def done(self, result: int) -> None:  # noqa: N802
        """Detach a still-running summary loader so the dialog can close."""
        loader = self._summary_loader
        if loader is not None and loader.isRunning():
            loader.ready.disconnect(self._on_summary_ready)
            loader.setParent(None)
            self._summary_loader = None
        super().done(result)

    # ------------------------------------------------------------------ defaults / results

    def _restore_defaults(self) -> None:
        d = TreeOptions()
        self._rexdb.set_path(Path(d.rexdb_fasta) if d.rexdb_fasta else None)
        self._rexdb_refs.setValue(d.rexdb_references)
        self._align_mode.setCurrentIndex(1 if d.align_mode == TREE_ALIGN_NEVER else 0)
        self._rotate.setChecked(d.rotate_monomers)
        self._set_combo(self._mask, d.methylation_mask)
        self._min_cov.setValue(round(d.min_coverage * 100))
        self._set_combo(self._cov_ref, d.coverage_reference)
        self._model.setCurrentIndex(0)
        self._runs.setValue(d.runs)
        self._seed.setValue(d.seed)
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
        self._sync()

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
        right_click = TREE_ENGINE_IQTREE
        for value, radio in self._right_click_buttons.items():
            if radio.isChecked():
                right_click = value
        alrt, ufboot = (int(x) for x in str(self._support_iq.currentData()).split(","))
        protein = self._protein()
        d = TreeOptions()
        engine = self._engine()
        threads = {TREE_ENGINE_IQTREE: self._threads, TREE_ENGINE_FASTME: self._threads_fme,
                   TREE_ENGINE_FASTTREE: self._threads_ft}[engine].value()
        sysbox = self._system_boxes
        return TreeOptions(
            align_mode=str(self._align_mode.currentData()),
            mafft_strategy=tb.MAFFT_CHOOSE,
            adjust_direction=True,
            rotate_monomers=self._rotate.isChecked(),
            sequence_mode=TREE_SEQ_AS_IS,
            rexdb_fasta=str(self._rexdb.path() or ""),
            rexdb_references=self._rexdb_refs.value(),
            methylation_mask="" if protein else str(self._mask.currentData() or ""),
            mafft_gap_open=d.mafft_gap_open,
            mafft_offset=d.mafft_offset,
            min_coverage=self._min_cov.value() / 100.0,
            coverage_reference=str(self._cov_ref.currentData()),
            engine=engine,
            model=model,
            runs=self._runs.value(),
            seed=self._seed.value(),
            threads=max(1, threads),
            ufboot=0 if ufboot == 0 else max(1000, ufboot),
            alrt=alrt,
            bnni=True,
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
