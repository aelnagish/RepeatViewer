"""Primer windows (Tools → Primers).

* :class:`PrimerDesignDialog`  – Primers | Regions | Specificity | Options
* :class:`MapPrimersDialog`    – Map Primers | Options (Geneious' Map Primers)
* :class:`PrimerDatabaseDialog` – the saved primers, with Tm / hairpin / dimer checks
* :class:`PrimerPickerDialog`  – choose a saved primer for a field

All follow the tool-window rules: Cancel / Save / Run, an Options tab, file
fields filled by a click in the Files panel, "?" help next to the item.
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QFont, QFontDatabase
from PySide6.QtWidgets import (
    QAbstractItemView,
    QButtonGroup,
    QCheckBox,
    QDialog,
    QFileDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QHeaderView,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPlainTextEdit,
    QPushButton,
    QRadioButton,
    QSizePolicy,
    QSplitter,
    QStackedWidget,
    QTabBar,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import design as dz
from repeatviewer.backends.primers import jobs as pj
from repeatviewer.backends.primers.oligos import (
    KIND_LABELS,
    PRIMER_KINDS,
    Primer,
    PrimerDatabase,
    parse_primer_text,
    write_primers,
)
from repeatviewer.backends.primers.thermo import (
    SALT_METHOD_LABELS,
    SALT_METHODS,
    TM_METHOD_LABELS,
    TM_METHODS,
    PrimerError,
    ThermoSettings,
    analyze_oligo,
    heterodimer,
    primer3_available,
    PRIMER3_MISSING,
)
from repeatviewer.job_options import (
    DIALOG_TAB_OPTIONS,
    PRIMER_DIALOG_TAB_PRIMERS,
    PRIMER_DIALOG_TAB_REGIONS,
    PRIMER_DIALOG_TAB_SPECIFICITY,
    PRIMER_DIALOG_TABS,
    PRIMER_RIGHT_CLICK_DIALOG,
    PRIMER_RIGHT_CLICK_RUN,
    PrimerOptions,
    opening_dialog_tab,
)
from repeatviewer.views.tool_dialogs import (
    TOOL_PARAMS_DIALOG_WIDTH_WIDE,
    TOOL_PARAMS_FIELD_WIDTH,
    _form_section_separator,
    dialog_should_run,
    make_cancel_save_run_bar,
    make_process_output_box,
    sync_form_field_column_spacers,
)
from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox
from repeatviewer.widgets.help_button import field_with_help
from repeatviewer.widgets.path_browse import PathBrowseRow
from repeatviewer.widgets.spin_boxes import NoWheelDoubleSpinBox as QDoubleSpinBox
from repeatviewer.widgets.spin_boxes import NoWheelSpinBox as QSpinBox

_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".fsa", ".ffn", ".seq", ".aln", ".mfa"}
_FASTA_FILTER = "FASTA (*.fa *.fasta *.fna *.fas *.fa.gz *.fasta.gz);;All files (*)"
_SMALL = 78          # width of one number box in a Min / Opt / Max row

# 5' tails offered next to each primer (lower case = does not bind).
TAIL_PRESETS: tuple[tuple[str, str], ...] = (
    ("None", ""),
    ("M13 forward (−20)", "gtaaaacgacggccagt"),
    ("M13 reverse", "caggaaacagctatgac"),
    ("T7 promoter", "taatacgactcactataggg"),
    ("SP6 promoter", "atttaggtgacactatag"),
    ("T3 promoter", "aattaaccctcactaaaggg"),
    ("EcoRI site (GAATTC)", "cgcggaattc"),
    ("BamHI site (GGATCC)", "cgcgggatcc"),
    ("HindIII site (AAGCTT)", "cgcgaagctt"),
    ("XhoI site (CTCGAG)", "cgcgctcgag"),
    ("XbaI site (TCTAGA)", "cgcgtctaga"),
    ("NcoI site (CCATGG)", "cgcgccatgg"),
    ("NotI site (GCGGCCGC)", "ataagaatgcggccgc"),
    ("Custom", None),
)


def _is_fasta(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return Path(name).suffix in _FASTA_SUFFIXES


def _mono_font(size: int | None = None) -> QFont:
    font = QFontDatabase.systemFont(QFontDatabase.SystemFont.FixedFont)
    if size:
        font.setPointSize(size)
    return font


def _set_combo(combo, value) -> None:
    idx = combo.findData(value)
    if idx >= 0:
        combo.setCurrentIndex(idx)


# --------------------------------------------------------------------------
# Help texts
# --------------------------------------------------------------------------

MODE_HELP = (
    "<b>Standard PCR</b>: a primer pair by Primer3, kept out of repeats (Regions tab).<br>"
    "<b>qPCR (copy number)</b>: short amplicons, matched Tm, optional hydrolysis probe; add a "
    "single-copy reference gene to get a pair for ΔCt copy-number estimates.<br>"
    "<b>Bisulfite PCR (BSP)</b>: primers on the bisulfite-converted strand with no CpG inside, "
    "so methylated and unmethylated DNA amplify alike.<br>"
    "<b>Methylation-specific PCR (MSP)</b>: an M pair and a U pair at the same place, their 3′ "
    "ends on CpG sites.<br>"
    "<b>Conserved in a family (MSA)</b>: on an alignment of copies, primers that catch (nearly) "
    "all of them, degenerate where needed.<br>"
    "<b>Subfamily-specific (MSA)</b>: primers that catch the sequences whose names match "
    "<i>Target sequences</i> and not the rest.<br>"
    "<b>Satellite ladder PCR</b>: outward-facing primers in one monomer. A tandem array then "
    "gives a ladder of 1×, 2×, 3× … monomer bands.<br>"
    "<b>FISH oligo probes</b>: 30–45 nt oligos found in every monomer / copy, without hairpins, "
    "ready to order with a 5′ label.<br>"
    "<b>TE insertion</b>: a flanking pair (empty site) and a primer inside the TE (occupied "
    "site) for presence / absence genotyping."
)
TEMPLATE_HELP = (
    "The FASTA the primers are designed on. Click a file in the Files panel to change it. "
    "Every sequence in the file is a template (up to 20); for the MSA modes the file must be "
    "an alignment, for ladder / FISH a set of monomers or copies (aligned or not)."
)
COVERAGE_HELP = (
    "A primer counts as binding a sequence when it has at most <b>Mismatches tolerated</b> "
    "mismatches and its last three 3′ bases match. <b>Min coverage</b> is the share of "
    "sequences a primer (and the pair) must bind."
)
DEGENERACY_HELP = (
    "Where the copies differ, a position may use an IUPAC code (R = A/G, Y = C/T …) covering "
    "every base seen in at least <b>Degenerate base ≥</b> of the sequences. <b>Max "
    "degeneracy</b> caps how many plain primers the mix contains (e.g. 16 = up to four "
    "two-fold positions)."
)
TARGET_HELP = (
    "Names (or part of them) of the sequences to amplify, e.g. <code>TEfamA</code> or a "
    "regular expression such as <code>^CroSat2_</code>. All other sequences in the alignment "
    "must not be amplified."
)
TE_HELP = (
    "Where the TE sits in the template (1-based). Leave 0 to use the longest repeat feature "
    "annotated on the sequence. Primers are designed on the flanks (empty site) and inside "
    "the first 300 bp of the TE (occupied site)."
)
MASK_HELP = (
    "Primers are kept out of these regions; the product may still span them. This keeps "
    "primers out of TEs and satellites, the commonest reason primers fail in repeat-rich "
    "genomes."
)
MISPRIME_HELP = (
    "A FASTA of repeats (your TE library, satellite monomers …). Primer3 rejects primers that "
    "would prime on any of them (its mispriming library)."
)
GENOME_HELP = (
    "Every primer of every pair is searched in this FASTA (both strands) and the products "
    "each pair would make are counted, like Primer-BLAST but local. On large files the "
    "search is seeded by the 3′ end: the last <b>3′ seed</b> bases must match exactly."
)
DROP_HELP = (
    "Drop pairs that make more products than this in the genome. 1 = only pairs with a "
    "single product. Untick for repeat primers, where many products are the point."
)
TAILS_HELP = (
    "A 5′ tail (adapter, promoter, restriction site with extra bases) added to the primer. "
    "Tails are written in lower case and do not count for the template Tm."
)
TM_HELP = (
    "Thermodynamic model and reaction conditions for every Tm and ΔG (Primer3)."
)
PRODUCT_HELP = (
    "Size range of the PCR product. For ladder PCR this is the smallest band; for TE "
    "insertion the empty-site band."
)


# --------------------------------------------------------------------------
# Small building blocks
# --------------------------------------------------------------------------


class RangeEdit(QWidget):
    """Min [–] Opt [–] Max number boxes on one row."""

    changed = Signal()

    def __init__(self, lo: float, hi: float, *, decimals: int = 0, opt: bool = True,
                 step: float = 1.0, suffix: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(4)
        self._boxes = []
        labels = ("Min", "Opt", "Max") if opt else ("Min", "Max")
        for i, lab in enumerate(labels):
            if decimals:
                box = QDoubleSpinBox(self)
                box.setDecimals(decimals)
            else:
                box = QSpinBox(self)
            box.setRange(lo, hi)
            box.setSingleStep(step)
            if suffix:
                box.setSuffix(suffix)
            box.setFixedWidth(_SMALL)
            box.setToolTip(lab)
            box.valueChanged.connect(lambda _v: self.changed.emit())
            if i:
                dash = QLabel(lab.lower() if opt else "to", self)
                dash.setEnabled(False)
                row.addWidget(dash)
            row.addWidget(box)
            self._boxes.append(box)
        self.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)

    def values(self) -> tuple:
        vals = tuple(b.value() for b in self._boxes)
        return vals

    def set_values(self, *values) -> None:
        for box, v in zip(self._boxes, values):
            box.blockSignals(True)
            box.setValue(v)
            box.blockSignals(False)

    def boxes(self) -> list:
        return list(self._boxes)


class ThermoBox:
    """Formula / salt correction / concentrations (Geneious' Tm Calculation).

    Not a widget of its own: :meth:`add_to` puts its rows into a dialog form
    so the labels line up with every other row.
    """

    def __init__(self, parent: QWidget) -> None:
        self._parent = parent
        self.formula = QComboBox(parent)
        for m in TM_METHODS:
            self.formula.addItem(TM_METHOD_LABELS[m], m)
        self.salt = QComboBox(parent)
        for m in SALT_METHODS:
            self.salt.addItem(SALT_METHOD_LABELS[m], m)

        def conc(lo, hi, step, suffix, decimals=1):
            b = QDoubleSpinBox(parent)
            b.setRange(lo, hi)
            b.setDecimals(decimals)
            b.setSingleStep(step)
            b.setSuffix(suffix)
            b.setFixedWidth(104)
            return b

        self.mono = conc(0, 2000, 5, " mM")
        self.di = conc(0, 200, 0.5, " mM", 2)
        self.dntp = conc(0, 50, 0.1, " mM", 2)
        self.oligo = conc(0.1, 100000, 10, " nM")
        self.set_settings(ThermoSettings())

    def add_to(self, form: "_Form") -> list[QWidget]:
        fields = [form.row("Formula:", self.formula, "Tm calculation", TM_HELP),
                  form.row("Salt correction:", self.salt)]

        def pair(a_label, a, b_label, b):
            w = QWidget(self._parent)
            h = QHBoxLayout(w)
            h.setContentsMargins(0, 0, 0, 0)
            h.setSpacing(6)
            h.addWidget(a)
            lab = QLabel(b_label, w)
            lab.setFixedWidth(52)
            lab.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            h.addWidget(lab)
            h.addWidget(b)
            h.addStretch(1)
            return form.row(a_label, w, fixed=False)

        fields.append(pair("Monovalent:", self.mono, "Oligo:", self.oligo))
        fields.append(pair("Divalent:", self.di, "dNTPs:", self.dntp))
        return fields

    def settings(self) -> ThermoSettings:
        return ThermoSettings(
            tm_method=str(self.formula.currentData()), salt_method=str(self.salt.currentData()),
            monovalent_mM=self.mono.value(), divalent_mM=self.di.value(),
            dntp_mM=self.dntp.value(), oligo_nM=self.oligo.value())

    def set_settings(self, t: ThermoSettings) -> None:
        _set_combo(self.formula, t.tm_method)
        _set_combo(self.salt, t.salt_method)
        self.mono.setValue(t.monovalent_mM)
        self.di.setValue(t.divalent_mM)
        self.dntp.setValue(t.dntp_mM)
        self.oligo.setValue(t.oligo_nM)


class TailEdit(QWidget):
    """Preset combo + editable tail sequence."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.combo = QComboBox(self)
        for label, seq in TAIL_PRESETS:
            self.combo.addItem(label, seq)
        self.combo.setFixedWidth(170)
        self.edit = QLineEdit(self)
        self.edit.setPlaceholderText("5′ tail (lower case)")
        self.edit.setFont(_mono_font())
        self.edit.setClearButtonEnabled(True)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        row.addWidget(self.combo)
        row.addWidget(self.edit, 1)
        self.combo.currentIndexChanged.connect(self._preset)
        self.edit.textEdited.connect(self._typed)

    def _preset(self, _i: int) -> None:
        seq = self.combo.currentData()
        if seq is None:
            return
        self.edit.setText(seq)

    def _typed(self, text: str) -> None:
        clean = "".join(ch for ch in text if ch.isalpha()).lower()
        if clean != text:
            pos = self.edit.cursorPosition()
            self.edit.setText(clean)
            self.edit.setCursorPosition(min(pos, len(clean)))
        match = next((i for i, (_l, s) in enumerate(TAIL_PRESETS) if s == clean), None)
        self.combo.blockSignals(True)
        self.combo.setCurrentIndex(match if match is not None else len(TAIL_PRESETS) - 1)
        self.combo.blockSignals(False)

    def tail(self) -> str:
        return "".join(ch for ch in self.edit.text() if ch.isalpha()).lower()

    def set_tail(self, tail: str) -> None:
        self.edit.setText(tail or "")
        self._typed(tail or "")


class _Form:
    """QFormLayout helper: bold sections, separators, help, checkbox rows."""

    def __init__(self, host: QWidget) -> None:
        self.host = host
        self.form = QFormLayout(host)
        self.form.setContentsMargins(0, 0, 0, 0)
        self.form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        self.form.setHorizontalSpacing(12)
        self.form.setVerticalSpacing(4)
        self.spacers: list[QLabel] = []
        self._first = True

    def section(self, title: str) -> list[QWidget]:
        out = []
        if not self._first:
            sep = _form_section_separator()
            self.form.addRow(sep)
            out.append(sep)
        self._first = False
        lab = QLabel(f"<b>{title}</b>", self.host)
        self.form.addRow(lab)
        out.append(lab)
        return out

    def row(self, label: str, widget: QWidget, help_title: str | None = None,
            help_text: str | None = None, *, fixed: bool = True) -> QWidget:
        if fixed and isinstance(widget, (QComboBox, QSpinBox, QDoubleSpinBox)):
            widget.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        field = field_with_help(widget, help_title or label.rstrip(":"), help_text) \
            if help_text else widget
        self.form.addRow(label, field)
        return field

    def check(self, box: QWidget, help_title: str | None = None,
              help_text: str | None = None) -> QWidget:
        spacer = QLabel("", self.host)
        self.spacers.append(spacer)
        field = field_with_help(box, help_title or "", help_text) if help_text else box
        self.form.addRow(spacer, field)
        return field

    def labelled(self, box: QWidget, widget: QWidget, help_title: str | None = None,
                 help_text: str | None = None, *, box_width: int = 150) -> QWidget:
        """Checkbox + its field in the field column (checkboxes line up, text on the
        right of the box as everywhere in RepeatViewer)."""
        host = QWidget(self.host)
        row = QHBoxLayout(host)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        box.setFixedWidth(box_width)
        row.addWidget(box)
        if widget.sizePolicy().horizontalPolicy() == QSizePolicy.Policy.Fixed:
            row.addWidget(widget)
            host.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        else:
            row.addWidget(widget, 1)
        return self.check(host, help_title, help_text)

    def widget_row(self, widget: QWidget) -> QWidget:
        spacer = QLabel("", self.host)
        self.spacers.append(spacer)
        self.form.addRow(spacer, widget)
        return widget

    def set_visible(self, field: QWidget, visible: bool) -> None:
        self.form.setRowVisible(field, visible)


def _path_row(parent: QWidget, placeholder: str, title: str, *, mode: str = "file",
              file_filter: str = _FASTA_FILTER) -> PathBrowseRow:
    row = PathBrowseRow(parent, placeholder=placeholder, browse_title=title, mode=mode,
                        file_filter=file_filter)
    row.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
    row.set_default_tooltip("Click a file in the Files panel to fill this path, or use Browse.")
    return row


def _page(parent: QWidget) -> tuple[QWidget, QVBoxLayout]:
    page = QWidget(parent)
    page.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
    lay = QVBoxLayout(page)
    lay.setContentsMargins(16, 10, 16, 10)
    lay.setSpacing(4)
    lay.setAlignment(Qt.AlignmentFlag.AlignTop)
    return page, lay


def _engine_box(parent: QWidget) -> QGroupBox:
    """Options tab: which Primer3 the windows use (like the Executables box)."""
    from repeatviewer.backends.primers.thermo import primer3_version

    box = QGroupBox("Engine", parent)
    box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
    lay = QVBoxLayout(box)
    lay.setContentsMargins(12, 10, 12, 10)
    version = primer3_version()
    text = (f"Primer3 {version} (primer3-py), thermodynamic alignment for hairpins and dimers"
            if version else "primer3-py is not installed")
    label = QLabel(text, box)
    label.setEnabled(False)
    label.setWordWrap(True)
    lay.addWidget(label)
    return box


class _TabbedToolDialog(QDialog):
    """Tab bar + pages + Cancel / Save / Run, sized to the visible page."""

    def __init__(self, title: str, width: int, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setModal(False)
        self.setWindowModality(Qt.WindowModality.NonModal)
        self.setWindowFlag(Qt.WindowType.WindowStaysOnTopHint, False)
        self._width = width
        self._tab_bar = QTabBar(self)
        self._tab_bar.setExpanding(False)
        self._tab_bar.setDocumentMode(True)
        self._stack = QStackedWidget(self)
        self._stack.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        self._forms: list[_Form] = []
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(22, 18, 22, 16)
        self._layout.setSpacing(8)
        self._layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        self._layout.addWidget(self._tab_bar, 0)
        self._layout.addWidget(self._stack, 0)
        self._tab_bar.currentChanged.connect(self._on_tab)

    def _add_page(self, title: str, page: QWidget) -> int:
        self._tab_bar.addTab(title)
        self._stack.addWidget(page)
        return self._tab_bar.count() - 1

    def _finish(self, bar: QWidget) -> None:
        self._layout.addWidget(bar, 0)
        self.resize(self._width, 400)

    def _on_tab(self, index: int) -> None:
        self._stack.setCurrentIndex(index)
        for i in range(self._stack.count()):
            w = self._stack.widget(i)
            w.setSizePolicy(QSizePolicy.Policy.Preferred,
                            QSizePolicy.Policy.Maximum if i == index else QSizePolicy.Policy.Ignored)
        if self.isVisible():
            QTimer.singleShot(0, self._fit)

    def _sync_columns(self) -> None:
        spacers = [s for f in self._forms for s in f.spacers]
        sync_form_field_column_spacers(self, spacers, [f.form for f in self._forms])

    def _fit(self) -> None:
        self._sync_columns()
        page = self._stack.currentWidget()
        if page is None:
            return
        page.adjustSize()
        self._stack.setFixedHeight(page.sizeHint().height())
        self.adjustSize()
        self.resize(self._width, self.sizeHint().height())

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._on_tab(self._tab_bar.currentIndex())
        QTimer.singleShot(0, self._fit)
        QTimer.singleShot(0, self._clear_path_focus)

    def _clear_path_focus(self) -> None:
        focused = self.focusWidget()
        if focused is not None and focused.objectName() == "PathField":
            self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)
            self.setFocus(Qt.FocusReason.OtherFocusReason)

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def _focused_path_row(self, rows: list[PathBrowseRow]) -> PathBrowseRow | None:
        focused = self.focusWidget()
        for row in rows:
            if row.isVisible() and focused is row.line_edit():
                return row
        return None


# --------------------------------------------------------------------------
# Design dialog
# --------------------------------------------------------------------------


class PrimerDesignDialog(_TabbedToolDialog):
    """Design primers: mode, conditions, repeat-aware regions, genome check, options."""

    _TAB_PRIMERS, _TAB_REGIONS, _TAB_SPECIFICITY, _TAB_OPTIONS = range(4)

    def __init__(self, *, initial: PrimerOptions, files: list[Path] | None = None,
                 project_genome: Path | None = None, parent: QWidget | None = None) -> None:
        super().__init__("Primer Design", TOOL_PARAMS_DIALOG_WIDTH_WIDE, parent)
        self._opts = initial
        self._project_genome = project_genome
        s = initial.design_settings()
        self._templates: list[Path] = [Path(p) for p in (files or []) if _is_fasta(Path(p))]
        self._template_features: list[str] = []

        # ---------------- Primers tab
        page, lay = _page(self)
        body = QWidget(page)
        f = _Form(body)
        self._forms.append(f)
        f.section("Inputs")
        self._template_label = QLabel(self)
        self._template_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._template_label.setWordWrap(True)
        f.row("Template:", self._template_label, "Template", TEMPLATE_HELP, fixed=False)
        self._mode = QComboBox(self)
        for m in dz.DESIGN_MODES:
            self._mode.addItem(dz.MODE_LABELS[m], m)
        _set_combo(self._mode, s.mode)
        f.row("Mode:", self._mode, "Design mode", MODE_HELP)

        # mode-specific settings
        self._mode_rows: dict[str, list[QWidget]] = {m: [] for m in dz.DESIGN_MODES}

        def mrow(modes, label, widget, help_title=None, help_text=None, fixed=True):
            field = f.row(label, widget, help_title, help_text, fixed=fixed)
            for m in modes:
                self._mode_rows[m].append(field)
            return field

        def mcheck(modes, box, help_title=None, help_text=None):
            field = f.check(box, help_title, help_text)
            for m in modes:
                self._mode_rows[m].append(field)
            return field

        self._probe = QCheckBox("Hydrolysis probe (TaqMan)", self)
        self._probe.setChecked(s.pick_probe)
        mcheck((dz.MODE_QPCR, dz.MODE_PCR), self._probe, "Probe",
               "Also pick an internal hydrolysis probe (Tm about 8–10 °C above the primers).")
        self._reference = _path_row(self, "Single-copy reference gene FASTA (optional)",
                                    "Choose the reference gene FASTA")
        self._reference.set_path(s.reference_fasta or None, emit=False)
        mrow((dz.MODE_QPCR,), "Reference gene:", self._reference, "Reference gene",
             "A single-copy gene (FASTA). A pair with the same Tm window is designed on it, "
             "for copy number = 2^(Ct reference − Ct repeat) × reference copies.", fixed=False)
        self._bs_strand = QComboBox(self)
        self._bs_strand.addItem("Top strand", "top")
        self._bs_strand.addItem("Bottom strand", "bottom")
        _set_combo(self._bs_strand, s.bisulfite_strand)
        mrow((dz.MODE_BSP, dz.MODE_MSP), "Converted strand:", self._bs_strand, "Strand",
             "After bisulfite treatment the two strands are no longer complementary: choose "
             "which converted strand the primers amplify.")
        self._min_conv = QSpinBox(self)
        self._min_conv.setRange(0, 15)
        self._min_conv.setValue(s.min_converted)
        mrow((dz.MODE_BSP,), "Min converted C:", self._min_conv, "Converted C",
             "Each primer must cover at least this many non-CpG C (T after conversion), so "
             "unconverted DNA is not amplified.")
        self._target_pattern = QLineEdit(s.target_pattern, self)
        self._target_pattern.setPlaceholderText("e.g. TEfamA or ^CroSat2_")
        self._target_pattern.setClearButtonEnabled(True)
        mrow((dz.MODE_SPECIFIC,), "Target sequences:", self._target_pattern, "Target sequences",
             TARGET_HELP, fixed=False)
        self._min_cov = QSpinBox(self)
        self._min_cov.setRange(10, 100)
        self._min_cov.setSuffix(" %")
        self._min_cov.setValue(int(round(s.min_coverage * 100)))
        mrow((dz.MODE_CONSERVED, dz.MODE_SPECIFIC, dz.MODE_FISH), "Min coverage:", self._min_cov,
             "Coverage", COVERAGE_HELP)
        self._cov_mm = QSpinBox(self)
        self._cov_mm.setRange(0, 8)
        self._cov_mm.setValue(s.coverage_mismatches)
        mrow((dz.MODE_CONSERVED, dz.MODE_SPECIFIC, dz.MODE_LADDER, dz.MODE_FISH),
             "Mismatches tolerated:", self._cov_mm, "Coverage", COVERAGE_HELP)
        self._max_off = QSpinBox(self)
        self._max_off.setRange(0, 100)
        self._max_off.setSuffix(" %")
        self._max_off.setValue(int(round(s.max_off_target * 100)))
        mrow((dz.MODE_SPECIFIC,), "Max off-target:", self._max_off, "Off-target",
             "Largest share of the non-target sequences a pair may still amplify.")
        self._max_deg = QSpinBox(self)
        self._max_deg.setRange(1, 4096)
        self._max_deg.setValue(s.max_degeneracy)
        mrow((dz.MODE_CONSERVED, dz.MODE_SPECIFIC), "Max degeneracy:", self._max_deg,
             "Degenerate primers", DEGENERACY_HELP)
        self._deg_freq = QSpinBox(self)
        self._deg_freq.setRange(1, 50)
        self._deg_freq.setSuffix(" %")
        self._deg_freq.setValue(int(round(s.degenerate_min_freq * 100)))
        mrow((dz.MODE_CONSERVED, dz.MODE_SPECIFIC), "Degenerate base ≥:", self._deg_freq,
             "Degenerate primers", DEGENERACY_HELP)
        self._fish_label = QComboBox(self)
        for lab in dz.FISH_LABELS:
            self._fish_label.addItem(lab, lab)
        _set_combo(self._fish_label, s.fish_label)
        mrow((dz.MODE_FISH,), "5′ label:", self._fish_label, "Label",
             "Written into the order line of each probe (5′-label-sequence-3′).")
        self._fish_avoid = _path_row(self, "Other repeats to avoid (FASTA, optional)",
                                     "Choose repeats the probes must not match")
        self._fish_avoid.set_path(s.fish_avoid_fasta or None, emit=False)
        mrow((dz.MODE_FISH,), "Avoid:", self._fish_avoid, "Cross-hybridization",
             "Probes matching any of these sequences (with the mismatches tolerated) are "
             "ranked last (e.g. the other satellite families).", fixed=False)
        te_row = QWidget(self)
        th = QHBoxLayout(te_row)
        th.setContentsMargins(0, 0, 0, 0)
        th.setSpacing(4)
        self._te_start = QSpinBox(self)
        self._te_end = QSpinBox(self)
        for b, v in ((self._te_start, s.te_start), (self._te_end, s.te_end)):
            b.setRange(0, 2_000_000_000)
            b.setFixedWidth(110)
            b.setSpecialValueText("auto")
            b.setValue(v)
        th.addWidget(self._te_start)
        to = QLabel("to", self)
        to.setEnabled(False)
        th.addWidget(to)
        th.addWidget(self._te_end)
        te_row.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._te_feature = QComboBox(self)
        self._te_feature.setEditable(True)
        self._te_feature.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self._te_feature.lineEdit().setPlaceholderText("Longest annotated repeat")
        self._te_feature.setEditText(s.te_feature)
        mrow((dz.MODE_INSERTION,), "TE feature:", self._te_feature, "TE feature",
             "An annotated feature of the template to use as the TE (its name). Ignored when "
             "TE position is set.")
        mrow((dz.MODE_INSERTION,), "TE position:", te_row, "TE position", TE_HELP, fixed=False)
        self._flank = QSpinBox(self)
        self._flank.setRange(50, 5000)
        self._flank.setSuffix(" bp")
        self._flank.setValue(s.flank_max)
        mrow((dz.MODE_INSERTION,), "Max flank distance:", self._flank, "Flank",
             "How far from the TE ends the flanking primers may sit.")

        f.section("Primer Conditions")
        self._size = RangeEdit(5, 60)
        self._size.set_values(s.size_min, s.size_opt, s.size_max)
        f.row("Primer size:", self._size, "Primer size", "Length in bases (min, optimum, max).",
              fixed=False)
        self._tm = RangeEdit(20, 95, decimals=1, step=0.5)
        self._tm.set_values(s.tm_min, s.tm_opt, s.tm_max)
        f.row("Tm (°C):", self._tm, "Tm", "Melting temperature of the binding part.",
              fixed=False)
        self._gc = RangeEdit(0, 100, decimals=0, opt=False, suffix=" %")
        self._gc.set_values(s.gc_min, s.gc_max)
        f.row("GC:", self._gc, fixed=False)
        self._product = RangeEdit(10, 100_000, opt=False)
        for b in self._product.boxes():
            b.setFixedWidth(96)
            b.setSuffix(" bp")
        self._product.set_values(s.product_min, s.product_max)
        f.row("Product size:", self._product, "Product size", PRODUCT_HELP, fixed=False)
        self._tm_diff = QDoubleSpinBox(self)
        self._tm_diff.setRange(0, 20)
        self._tm_diff.setDecimals(1)
        self._tm_diff.setSuffix(" °C")
        self._tm_diff.setValue(s.max_tm_diff)
        f.row("Max Tm difference:", self._tm_diff)
        self._num = QSpinBox(self)
        self._num.setRange(1, 100)
        self._num.setValue(s.num_return)
        f.row("Pairs to return:", self._num)

        # Advanced
        self._adv_fields: list[QWidget] = []
        self._adv_fields += f.section("Advanced Settings")
        self._poly = QSpinBox(self)
        self._poly.setRange(1, 20)
        self._poly.setValue(s.max_poly_x)
        self._adv_fields.append(f.row("Max poly-X:", self._poly, "Poly-X",
                                      "Longest run of one base allowed in a primer."))
        self._clamp = QSpinBox(self)
        self._clamp.setRange(0, 5)
        self._clamp.setValue(s.gc_clamp)
        self._adv_fields.append(f.row("GC clamp:", self._clamp, "GC clamp",
                                      "G or C bases required at the 3′ end."))
        self._end_stab = QDoubleSpinBox(self)
        self._end_stab.setRange(0, 30)
        self._end_stab.setDecimals(1)
        self._end_stab.setValue(s.max_end_stability)
        self._adv_fields.append(f.row("Max 3′ stability:", self._end_stab, "3′ stability",
                                      "Largest −ΔG (kcal/mol) of the last five 3′ bases "
                                      "(Primer3 default 9). Lower = less mispriming."))
        self._probe_size = RangeEdit(10, 60)
        self._probe_size.set_values(s.probe_size_min, s.probe_size_opt, s.probe_size_max)
        self._adv_fields.append(f.row("Probe size:", self._probe_size, fixed=False))
        self._probe_tm = RangeEdit(30, 95, decimals=1, step=0.5)
        self._probe_tm.set_values(s.probe_tm_min, s.probe_tm_opt, s.probe_tm_max)
        self._adv_fields.append(f.row("Probe Tm (°C):", self._probe_tm, fixed=False))
        self._fwd_tail = TailEdit(self)
        self._fwd_tail.set_tail(s.forward_tail)
        self._adv_fields.append(f.row("Forward 5′ tail:", self._fwd_tail, "5′ tails", TAILS_HELP,
                                      fixed=False))
        self._rev_tail = TailEdit(self)
        self._rev_tail.set_tail(s.reverse_tail)
        self._adv_fields.append(f.row("Reverse 5′ tail:", self._rev_tail, fixed=False))
        self._adv_fields += f.section("Tm Calculation")
        self._thermo = ThermoBox(self)
        self._thermo.set_settings(s.thermo)
        self._adv_fields += self._thermo.add_to(f)
        lay.addWidget(body)
        self._primers_form = f
        self._add_page("Primers", page)

        # ---------------- Regions tab
        page, lay = _page(self)
        body = QWidget(page)
        g = _Form(body)
        self._forms.append(g)
        g.section("Regions")
        self._use_target = QCheckBox("Target region (the product must span it)", self)
        self._use_target.setChecked(bool(s.target_start and s.target_end))
        g.check(self._use_target)
        self._target = RangeEdit(1, 2_000_000_000, opt=False)
        for b in self._target.boxes():
            b.setFixedWidth(110)
        self._target.set_values(s.target_start or 1, s.target_end or 1)
        g.row("Target:", self._target, fixed=False)
        self._use_included = QCheckBox("Included region (primers only inside it)", self)
        self._use_included.setChecked(bool(s.included_start and s.included_end))
        g.check(self._use_included)
        self._included = RangeEdit(1, 2_000_000_000, opt=False)
        for b in self._included.boxes():
            b.setFixedWidth(110)
        self._included.set_values(s.included_start or 1, s.included_end or 1)
        g.row("Included:", self._included, fixed=False)
        self._use_target.toggled.connect(self._target.setEnabled)
        self._use_included.toggled.connect(self._included.setEnabled)
        self._target.setEnabled(self._use_target.isChecked())
        self._included.setEnabled(self._use_included.isChecked())
        g.section("Repeat Masking")
        self._mask_lower = QCheckBox("Avoid soft-masked (lower-case) bases", self)
        self._mask_lower.setChecked(s.mask_lowercase)
        g.check(self._mask_lower, "Repeat masking", MASK_HELP)
        self._mask_feats = QCheckBox("Avoid repeats annotated on the sequence", self)
        self._mask_feats.setChecked(s.mask_features)
        g.check(self._mask_feats)
        self._mask_types = QLineEdit(s.mask_feature_types, self)
        self._mask_types.setPlaceholderText("all except ORF / domain / primers")
        self._mask_types.setClearButtonEnabled(True)
        g.row("Feature types:", self._mask_types, "Feature types",
              "Comma-separated feature types to avoid (e.g. <code>LTR, TIR, tandem repeat</code>). "
              "Empty = every annotated feature except ORFs, domains, genes and primers.",
              fixed=False)
        self._mask_tandem = QCheckBox("Avoid tandem repeats (found with pytrf)", self)
        self._mask_tandem.setChecked(s.mask_tandem)
        g.check(self._mask_tandem)
        self._mask_file = _path_row(self, "RepeatMasker .out, GFF or BED (optional)",
                                    "Choose a repeat annotation file",
                                    file_filter="Repeat annotation (*.out *.gff *.gff3 *.bed);;"
                                                "All files (*)")
        self._mask_file.set_path(s.mask_file or None, emit=False)
        g.row("Repeat file:", self._mask_file, "Repeat file",
              "Repeats in this file (rows whose sequence name matches the template) are "
              "avoided too.", fixed=False)
        self._misprime = _path_row(self, "Repeat library FASTA (optional)",
                                   "Choose a mispriming library")
        self._misprime.set_path(s.mispriming_library or None, emit=False)
        g.row("Mispriming library:", self._misprime, "Mispriming library", MISPRIME_HELP,
              fixed=False)
        lay.addWidget(body)
        self._add_page("Regions", page)

        # ---------------- Specificity tab
        page, lay = _page(self)
        body = QWidget(page)
        h = _Form(body)
        self._forms.append(h)
        h.section("Genome Specificity")
        self._check_genome = QCheckBox("Check every pair against a genome", self)
        self._check_genome.setChecked(s.check_genome)
        h.check(self._check_genome, "Genome check", GENOME_HELP)
        self._genome = _path_row(self, "Genome FASTA: click it in Files, or Browse",
                                 "Choose the genome FASTA")
        genome = s.genome_fasta or (str(project_genome) if project_genome else "")
        self._genome.set_path(genome or None, emit=False)
        h.row("Genome:", self._genome, fixed=False)
        self._g_mm = QSpinBox(self)
        self._g_mm.setRange(0, 6)
        self._g_mm.setValue(s.genome_max_mismatches)
        h.row("Max mismatches:", self._g_mm, "Mismatches",
              "Mismatches a primer may have and still count as binding.")
        self._g_exact = QSpinBox(self)
        self._g_exact.setRange(0, 15)
        self._g_exact.setValue(s.genome_three_prime_exact)
        h.row("3′ bases exact:", self._g_exact, "3′ end",
              "The last bases at the 3′ end must match: mismatches there stop the polymerase.")
        self._g_seed = QSpinBox(self)
        self._g_seed.setRange(6, 12)
        self._g_seed.setValue(s.genome_seed)
        h.row("3′ seed:", self._g_seed, "Seed", GENOME_HELP)
        self._g_maxprod = QSpinBox(self)
        self._g_maxprod.setRange(50, 100_000)
        self._g_maxprod.setSuffix(" bp")
        self._g_maxprod.setValue(s.genome_max_product)
        h.row("Max product size:", self._g_maxprod, "Product size",
              "Products longer than this are not counted (they would not amplify well).")
        drop_row = QWidget(self)
        dr = QHBoxLayout(drop_row)
        dr.setContentsMargins(0, 0, 0, 0)
        dr.setSpacing(6)
        self._drop = QCheckBox("Drop pairs with more than", self)
        self._drop.setChecked(s.drop_unspecific)
        self._drop_n = QSpinBox(self)
        self._drop_n.setRange(1, 1_000_000)
        self._drop_n.setValue(max(1, s.genome_max_products))
        self._drop_n.setFixedWidth(_SMALL + 10)
        dr.addWidget(self._drop)
        dr.addWidget(self._drop_n)
        dr.addWidget(QLabel("product(s)", self))
        drop_row.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        h.check(drop_row, "Drop unspecific pairs", DROP_HELP)
        self._check_genome.toggled.connect(self._sync_genome)
        self._drop.toggled.connect(self._drop_n.setEnabled)
        lay.addWidget(body)
        self._add_page("Specificity", page)

        # ---------------- Options tab
        page, lay = _page(self)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.setSpacing(14)
        rc_box = QGroupBox("Right-click on Primers", self)
        rcl = QVBoxLayout(rc_box)
        self._rc_run = QRadioButton("Design with saved settings", self)
        self._rc_dialog = QRadioButton("Open this window", self)
        grp = QButtonGroup(self)
        for b in (self._rc_run, self._rc_dialog):
            grp.addButton(b)
            rcl.addWidget(b)
        (self._rc_dialog if initial.right_click_action == PRIMER_RIGHT_CLICK_DIALOG
         else self._rc_run).setChecked(True)
        res_box = QGroupBox("Results", self)
        rl = QVBoxLayout(res_box)
        self._annotate = QCheckBox("Add the best pair to the sequence annotations", self)
        self._annotate.setChecked(initial.annotate_design)
        self._to_db = QCheckBox("Add the designed primers to the Primer Database", self)
        self._to_db.setChecked(initial.save_design_to_database)
        rl.addWidget(self._annotate)
        rl.addWidget(self._to_db)
        db_row = QWidget(self)
        dbl = QHBoxLayout(db_row)
        dbl.setContentsMargins(0, 4, 0, 0)
        dbl.setSpacing(6)
        dbl.addWidget(QLabel("Primer Database:", self))
        self._db_path = _path_row(self, "Default (RepeatViewer settings folder)",
                                  "Choose the Primer Database FASTA")
        self._db_path.set_path(initial.database_path or None, emit=False)
        dbl.addWidget(self._db_path, 1)
        rl.addWidget(db_row)
        process_box, self._show_process, self._close_process = make_process_output_box(
            self, tool_name="primer design", checked=initial.show_process_window,
            close_when_finished=initial.close_process_window_when_finished)
        lay.addWidget(rc_box)
        lay.addWidget(res_box)
        lay.addWidget(_engine_box(self))
        lay.addWidget(process_box)
        self._add_page("Options", page)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._adv_btn = QPushButton("Advanced Settings ▸", self)
        self._adv_btn.setCheckable(True)
        self._adv_btn.toggled.connect(self._toggle_advanced)
        self._finish(make_cancel_save_run_bar(self, left_widgets=[restore, self._adv_btn]))
        self._mode.currentIndexChanged.connect(self._on_mode)
        open_tab = opening_dialog_tab(initial.last_dialog_tab, allowed=PRIMER_DIALOG_TABS,
                                      first_tab=PRIMER_DIALOG_TAB_PRIMERS)
        self._tab_bar.setCurrentIndex({PRIMER_DIALOG_TAB_PRIMERS: 0, PRIMER_DIALOG_TAB_REGIONS: 1,
                                       PRIMER_DIALOG_TAB_SPECIFICITY: 2}.get(open_tab, 0))
        self._apply_mode_rows()
        self._toggle_advanced(False)
        self._sync_genome()
        self._drop_n.setEnabled(self._drop.isChecked())
        self._refresh_template_label()
        self._load_template_features()

    # ------------------------------------------------------------- behaviour

    def _all_path_rows(self) -> list[PathBrowseRow]:
        return [self._reference, self._fish_avoid, self._mask_file, self._misprime,
                self._genome, self._db_path]

    def _on_tab(self, index: int) -> None:
        if hasattr(self, "_adv_btn"):
            self._adv_btn.setVisible(index == self._TAB_PRIMERS)
        super()._on_tab(index)

    def _toggle_advanced(self, on: bool) -> None:
        self._adv_btn.setText("Advanced Settings ▾" if on else "Advanced Settings ▸")
        for w in self._adv_fields:
            self._primers_form.set_visible(w, on)
        if self.isVisible():
            QTimer.singleShot(0, self._fit)

    def _sync_genome(self, *_a) -> None:
        on = self._check_genome.isChecked()
        for w in (self._genome, self._g_mm, self._g_exact, self._g_seed, self._g_maxprod,
                  self._drop, self._drop_n):
            w.setEnabled(on)
        if on:
            self._drop_n.setEnabled(self._drop.isChecked())

    def _on_mode(self, _i: int) -> None:
        mode = str(self._mode.currentData())
        self._apply_preset(mode)
        self._apply_mode_rows()

    def _apply_mode_rows(self) -> None:
        mode = str(self._mode.currentData())
        shown = set(id(w) for w in self._mode_rows.get(mode, []))
        every = {id(w): w for rows in self._mode_rows.values() for w in rows}
        for key, w in every.items():
            self._primers_form.set_visible(w, key in shown)
        # probe sizes only matter with a probe
        self._refresh_template_label()
        if self.isVisible():
            QTimer.singleShot(0, self._fit)

    def _apply_preset(self, mode: str) -> None:
        p = dz.mode_preset(mode)
        self._size.set_values(p.get("size_min", 18), p.get("size_opt", 20), p.get("size_max", 25))
        self._tm.set_values(p.get("tm_min", 57.0), p.get("tm_opt", 60.0), p.get("tm_max", 63.0))
        self._gc.set_values(p.get("gc_min", 40.0), p.get("gc_max", 60.0))
        self._product.set_values(p.get("product_min", 100), p.get("product_max", 1000))
        self._tm_diff.setValue(p.get("max_tm_diff", 2.0))
        if "max_poly_x" in p:
            self._poly.setValue(p["max_poly_x"])
        if "gc_clamp" in p:
            self._clamp.setValue(p["gc_clamp"])
        if "coverage_mismatches" in p:
            self._cov_mm.setValue(p["coverage_mismatches"])
        if "min_coverage" in p:
            self._min_cov.setValue(int(round(p["min_coverage"] * 100)))
        if mode == dz.MODE_QPCR:
            self._probe.setChecked(True)

    def _restore_defaults(self) -> None:
        mode = str(self._mode.currentData())
        d = dz.DesignSettings(mode=mode, **dz.mode_preset(mode))
        self._load_settings(d)

    def _load_settings(self, s: dz.DesignSettings) -> None:
        _set_combo(self._mode, s.mode)
        self._size.set_values(s.size_min, s.size_opt, s.size_max)
        self._tm.set_values(s.tm_min, s.tm_opt, s.tm_max)
        self._gc.set_values(s.gc_min, s.gc_max)
        self._product.set_values(s.product_min, s.product_max)
        self._tm_diff.setValue(s.max_tm_diff)
        self._num.setValue(s.num_return)
        self._poly.setValue(s.max_poly_x)
        self._clamp.setValue(s.gc_clamp)
        self._end_stab.setValue(s.max_end_stability)
        self._probe.setChecked(s.pick_probe)
        self._probe_size.set_values(s.probe_size_min, s.probe_size_opt, s.probe_size_max)
        self._probe_tm.set_values(s.probe_tm_min, s.probe_tm_opt, s.probe_tm_max)
        self._fwd_tail.set_tail(s.forward_tail)
        self._rev_tail.set_tail(s.reverse_tail)
        self._thermo.set_settings(s.thermo)
        self._mask_lower.setChecked(s.mask_lowercase)
        self._mask_feats.setChecked(s.mask_features)
        self._mask_types.setText(s.mask_feature_types)
        self._mask_tandem.setChecked(s.mask_tandem)
        self._g_mm.setValue(s.genome_max_mismatches)
        self._g_exact.setValue(s.genome_three_prime_exact)
        self._g_seed.setValue(s.genome_seed)
        self._g_maxprod.setValue(s.genome_max_product)
        self._drop.setChecked(s.drop_unspecific)
        self._drop_n.setValue(max(1, s.genome_max_products))
        self._min_cov.setValue(int(round(s.min_coverage * 100)))
        self._cov_mm.setValue(s.coverage_mismatches)
        self._max_off.setValue(int(round(s.max_off_target * 100)))
        self._max_deg.setValue(s.max_degeneracy)
        self._deg_freq.setValue(int(round(s.degenerate_min_freq * 100)))
        self._min_conv.setValue(s.min_converted)
        self._flank.setValue(s.flank_max)
        self._apply_mode_rows()

    # ------------------------------------------------------------- inputs

    def _refresh_template_label(self) -> None:
        if not self._templates:
            self._template_label.setText(
                "<span style='color:palette(mid)'>Click a FASTA in the Files panel</span>")
            return
        names = ", ".join(p.name for p in self._templates[:3])
        more = f" +{len(self._templates) - 3}" if len(self._templates) > 3 else ""
        self._template_label.setText(f"<b>{names}</b>{more}")
        self._template_label.setToolTip("\n".join(str(p) for p in self._templates))

    def template_files(self) -> list[Path]:
        return list(self._templates)

    def set_templates(self, files: list[Path]) -> None:
        self._templates = [Path(p) for p in files if _is_fasta(Path(p))]
        self._refresh_template_label()
        self._load_template_features()

    def _load_template_features(self) -> None:
        """Annotated features of the template, offered as TE candidates."""
        names: list[str] = []
        if self._templates:
            feats = pj.record_features(self._templates[0])
            skip = {"primer_bind", "pcr product", "probe", "orf", "domain"}
            for flist in feats.values():
                for feat in sorted(flist, key=lambda x: x.start - x.end):
                    if str(feat.feature_type).lower() not in skip and feat.name not in names:
                        names.append(feat.name)
        text = self._te_feature.currentText()
        self._te_feature.blockSignals(True)
        self._te_feature.clear()
        self._te_feature.addItems(names[:200])
        self._te_feature.setCurrentIndex(-1)
        self._te_feature.setEditText(text)
        self._te_feature.blockSignals(False)

    def apply_path_from_panel(self, path: Path) -> bool:
        """Files-panel click: fills the path field being edited, else the template."""
        resolved = Path(path).expanduser()
        if not resolved.is_file():
            return False
        row = self._focused_path_row(self._all_path_rows())
        if row is not None:
            row.set_path(resolved)
            return True
        if not _is_fasta(resolved):
            return False
        self.set_templates([resolved])
        return True

    def current_dialog_tab(self) -> str:
        return {0: PRIMER_DIALOG_TAB_PRIMERS, 1: PRIMER_DIALOG_TAB_REGIONS,
                2: PRIMER_DIALOG_TAB_SPECIFICITY}.get(self._tab_bar.currentIndex(),
                                                      DIALOG_TAB_OPTIONS)

    def design_settings(self) -> dz.DesignSettings:
        smin, sopt, smax = self._size.values()
        tmin, topt, tmax = self._tm.values()
        gmin, gmax = self._gc.values()
        pmin, pmax = self._product.values()
        psmin, psopt, psmax = self._probe_size.values()
        ptmin, ptopt, ptmax = self._probe_tm.values()
        tgt = self._target.values() if self._use_target.isChecked() else (0, 0)
        inc = self._included.values() if self._use_included.isChecked() else (0, 0)

        def path(row: PathBrowseRow) -> str:
            p = row.path()
            return str(p) if p else ""

        return dz.DesignSettings(
            mode=str(self._mode.currentData()), size_min=smin, size_opt=sopt, size_max=smax,
            tm_min=tmin, tm_opt=topt, tm_max=tmax, gc_min=gmin, gc_max=gmax,
            max_tm_diff=self._tm_diff.value(), product_min=pmin, product_max=pmax,
            num_return=self._num.value(), max_poly_x=self._poly.value(),
            gc_clamp=self._clamp.value(), max_end_stability=self._end_stab.value(),
            pick_probe=self._probe.isChecked() and str(self._mode.currentData()) in
            (dz.MODE_PCR, dz.MODE_QPCR),
            probe_size_min=psmin, probe_size_opt=psopt, probe_size_max=psmax,
            probe_tm_min=ptmin, probe_tm_opt=ptopt, probe_tm_max=ptmax,
            thermo=self._thermo.settings(),
            target_start=tgt[0], target_end=tgt[1], included_start=inc[0], included_end=inc[1],
            mask_lowercase=self._mask_lower.isChecked(), mask_features=self._mask_feats.isChecked(),
            mask_feature_types=self._mask_types.text().strip(),
            mask_tandem=self._mask_tandem.isChecked(), mask_file=path(self._mask_file),
            mispriming_library=path(self._misprime),
            check_genome=self._check_genome.isChecked(), genome_fasta=path(self._genome),
            genome_max_mismatches=self._g_mm.value(),
            genome_three_prime_exact=self._g_exact.value(), genome_seed=self._g_seed.value(),
            genome_max_product=self._g_maxprod.value(),
            genome_max_products=self._drop_n.value(), drop_unspecific=self._drop.isChecked(),
            forward_tail=self._fwd_tail.tail(), reverse_tail=self._rev_tail.tail(),
            reference_fasta=path(self._reference),
            bisulfite_strand=str(self._bs_strand.currentData()),
            min_converted=self._min_conv.value(),
            min_coverage=self._min_cov.value() / 100.0,
            max_degeneracy=self._max_deg.value(),
            degenerate_min_freq=self._deg_freq.value() / 100.0,
            coverage_mismatches=self._cov_mm.value(),
            target_pattern=self._target_pattern.text().strip(),
            max_off_target=self._max_off.value() / 100.0,
            fish_label=str(self._fish_label.currentData()),
            fish_avoid_fasta=path(self._fish_avoid),
            te_start=self._te_start.value(), te_end=self._te_end.value(),
            te_feature=self._te_feature.currentText().strip(),
            flank_max=self._flank.value(),
        )

    def options(self) -> PrimerOptions:
        db = self._db_path.path()
        return replace(
            self._opts,
            design=self.design_settings().to_dict(),
            database_path=str(db) if db else "",
            right_click_action=(PRIMER_RIGHT_CLICK_DIALOG if self._rc_dialog.isChecked()
                                else PRIMER_RIGHT_CLICK_RUN),
            annotate_design=self._annotate.isChecked(),
            save_design_to_database=self._to_db.isChecked(),
            show_process_window=self._show_process.isChecked(),
            close_process_window_when_finished=self._close_process.isChecked(),
            last_dialog_tab=self.current_dialog_tab(),
        )


# --------------------------------------------------------------------------
# Map Primers dialog
# --------------------------------------------------------------------------


class PrimerPickerDialog(QDialog):
    """Pick a saved primer (filterable list)."""

    def __init__(self, db: PrimerDatabase, *, kind: str | None = None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Choose a Primer")
        self.setMinimumSize(560, 380)
        self._db = db
        self._chosen: Primer | None = None
        self._filter = QLineEdit(self)
        self._filter.setPlaceholderText("Filter by name, sequence or note")
        self._filter.setClearButtonEnabled(True)
        self._table = QTableWidget(0, 4, self)
        self._table.setHorizontalHeaderLabels(["Name", "Kind", "Sequence 5′→3′", "Note"])
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.SingleSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.verticalHeader().setVisible(False)
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.setAlternatingRowColors(True)
        self._table.doubleClicked.connect(lambda _i: self._accept())
        self._filter.textChanged.connect(self._fill)
        cancel = QPushButton("Cancel", self)
        choose = QPushButton("Choose", self)
        choose.setDefault(True)
        cancel.clicked.connect(self.reject)
        choose.clicked.connect(self._accept)
        bar = QHBoxLayout()
        bar.addStretch(1)
        bar.addWidget(cancel)
        bar.addWidget(choose)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.addWidget(self._filter)
        lay.addWidget(self._table, 1)
        lay.addLayout(bar)
        self._kind = kind
        self._fill()

    def _fill(self, *_a) -> None:
        text = self._filter.text().strip().lower()
        rows = [p for p in self._db.primers
                if not text or text in p.name.lower() or text in p.sequence.lower()
                or text in p.note.lower()]
        if self._kind:
            rows.sort(key=lambda p: (p.kind != self._kind, p.name.lower()))
        self._rows = rows
        self._table.setRowCount(len(rows))
        mono = _mono_font()
        for i, p in enumerate(rows):
            for j, val in enumerate((p.name, KIND_LABELS.get(p.kind, p.kind), p.sequence, p.note)):
                item = QTableWidgetItem(val)
                if j == 2:
                    item.setFont(mono)
                self._table.setItem(i, j, item)
        self._table.resizeColumnsToContents()
        if rows:
            self._table.selectRow(0)

    def _accept(self) -> None:
        row = self._table.currentRow()
        if 0 <= row < len(self._rows):
            self._chosen = self._rows[row]
            self.accept()

    def chosen(self) -> Primer | None:
        return self._chosen


class _PrimerField(QWidget):
    """Editable combo (saved names or a typed sequence) + Choose…"""

    def __init__(self, kind: str, db_getter, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._kind = kind
        self._db_getter = db_getter
        self.combo = QComboBox(self)
        self.combo.setEditable(True)
        self.combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self.combo.lineEdit().setPlaceholderText("Saved primer name, or type 5′→3′ sequence")
        self.combo.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.choose = QPushButton("Choose…", self)
        self.choose.clicked.connect(self._pick)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(6)
        row.addWidget(self.combo, 1)
        row.addWidget(self.choose)

    def reload(self) -> None:
        text = self.text()
        db = self._db_getter()
        self.combo.blockSignals(True)
        self.combo.clear()
        if db is not None:
            for p in sorted(db.primers, key=lambda x: (x.kind != self._kind, x.name.lower())):
                self.combo.addItem(p.name, p.name)
        self.combo.setCurrentIndex(-1)
        self.combo.setEditText(text)
        self.combo.blockSignals(False)

    def text(self) -> str:
        return self.combo.currentText().strip()

    def set_text(self, text: str) -> None:
        self.combo.setEditText(text or "")

    def _pick(self) -> None:
        db = self._db_getter()
        if db is None or not db.primers:
            QMessageBox.information(self, "Primer Database",
                                    "The Primer Database is empty. Add primers in Tools → "
                                    "Primers → Primer Database, or type a sequence here.")
            return
        dlg = PrimerPickerDialog(db, kind=self._kind, parent=self.window())
        if dlg.exec() == int(QDialog.DialogCode.Accepted) and dlg.chosen() is not None:
            self.set_text(dlg.chosen().name)


class MapPrimersDialog(_TabbedToolDialog):
    """Test primers against the selected sequences (Geneious' Map Primers)."""

    _TAB_MAP, _TAB_OPTIONS = range(2)

    def __init__(self, *, initial: PrimerOptions, files: list[Path] | None = None,
                 project_genome: Path | None = None, parent: QWidget | None = None) -> None:
        super().__init__("Map Primers", TOOL_PARAMS_DIALOG_WIDTH_WIDE, parent)
        self._opts = initial
        ms = initial.map_settings()
        self._files = [Path(p) for p in (files or []) if _is_fasta(Path(p))]
        self._db: PrimerDatabase | None = None

        page, lay = _page(self)
        body = QWidget(page)
        f = _Form(body)
        self._forms.append(f)
        f.section("Inputs")
        self._seq_label = QLabel(self)
        self._seq_label.setWordWrap(True)
        f.row("Sequences:", self._seq_label, "Sequences",
              "The FASTA files selected in the Files panel when you click Run (click a file "
              "to choose it). Each sequence is tested.", fixed=False)
        self._test_on = QComboBox(self)
        self._test_on.addItem("Sequences", pj.TEST_ON_SEQUENCES)
        self._test_on.addItem("Consensus (alignments)", pj.TEST_ON_CONSENSUS)
        _set_combo(self._test_on, ms.test_on)
        f.row("Test primers on:", self._test_on, "Test primers on",
              "For an alignment, test the primers on its consensus instead of on every row.")
        self._cons_threshold = QComboBox(self)
        for label, value in (("0 % – Majority", 0.0), ("25 %", 0.25), ("50 %", 0.5),
                             ("75 %", 0.75), ("95 %", 0.95), ("100 % – Identical", 1.0)):
            self._cons_threshold.addItem(label, value)
        _set_combo(self._cons_threshold, ms.consensus_threshold)
        self._cons_field = f.row("Consensus threshold:", self._cons_threshold, "Consensus",
                                 "Bases that together reach the threshold are written as one "
                                 "IUPAC code; 0 % writes the most common base.")
        self._cons_gaps = QCheckBox("Ignore gaps", self)
        self._cons_gaps.setChecked(ms.consensus_ignore_gaps)
        self._cons_gaps_field = f.check(self._cons_gaps)
        self._test_on.currentIndexChanged.connect(self._sync)

        f.section("Primers")
        self._specific = QRadioButton("Test specific primers against the selected sequences",
                                      self)
        self._database = QRadioButton("Search for saved primers that match the selected "
                                      "sequences", self)
        grp = QButtonGroup(self)
        grp.addButton(self._specific)
        grp.addButton(self._database)
        (self._database if ms.mode == pj.MAP_DATABASE else self._specific).setChecked(True)
        f.check(self._specific)
        self._use_f = QCheckBox("Forward primer", self)
        self._use_r = QCheckBox("Reverse primer", self)
        self._use_p = QCheckBox("Probe", self)
        self._f = _PrimerField("forward", self._database_obj, self)
        self._r = _PrimerField("reverse", self._database_obj, self)
        self._p = _PrimerField("probe", self._database_obj, self)
        self._spec_rows = []
        for box, field, on, text in ((self._use_f, self._f, ms.use_forward, ms.forward),
                                     (self._use_r, self._r, ms.use_reverse, ms.reverse),
                                     (self._use_p, self._p, ms.use_probe, ms.probe)):
            box.setChecked(on)
            field.set_text(text)
            field.setEnabled(on)
            box.toggled.connect(field.setEnabled)
            self._spec_rows.append(f.labelled(box, field))
        f.check(self._database)
        src = QWidget(self)
        sl = QHBoxLayout(src)
        sl.setContentsMargins(0, 0, 0, 0)
        sl.setSpacing(6)
        self._src = QComboBox(self)
        self._src.addItem("Primer Database", "db")
        self._src.addItem("FASTA file", "file")
        self._src.setFixedWidth(150)
        self._src_file = _path_row(self, "Primers FASTA", "Choose a primers FASTA")
        if ms.primer_file:
            self._src.setCurrentIndex(1)
            self._src_file.set_path(ms.primer_file, emit=False)
        sl.addWidget(self._src)
        sl.addWidget(self._src_file, 1)
        self._src_field = f.row("Primers from:", src, fixed=False)
        kinds = QWidget(self)
        kl = QHBoxLayout(kinds)
        kl.setContentsMargins(0, 0, 0, 0)
        kl.setSpacing(14)
        self._s_f = QCheckBox("Forward", self)
        self._s_r = QCheckBox("Reverse", self)
        self._s_p = QCheckBox("Probe", self)
        self._pairs_only = QCheckBox("Pairs only", self)
        for b, v in ((self._s_f, ms.search_forward), (self._s_r, ms.search_reverse),
                     (self._s_p, ms.search_probe), (self._pairs_only, ms.pairs_only)):
            b.setChecked(v)
            kl.addWidget(b)
        kinds.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._kinds_field = f.row("Search for:", kinds, "Search for",
                                  "Which saved primers to test. <b>Pairs only</b> keeps primers "
                                  "that form a product with another saved primer.", fixed=False)
        self._specific.toggled.connect(self._sync)
        self._src.currentIndexChanged.connect(self._sync)

        f.section("Products")

        def region(check_text, on, a, b):
            box = QCheckBox(check_text, self)
            box.setChecked(on)
            rng = RangeEdit(1, 2_000_000_000, opt=False)
            for x in rng.boxes():
                x.setFixedWidth(110)
            rng.set_values(a, b)
            rng.setEnabled(on)
            box.toggled.connect(rng.setEnabled)
            f.labelled(box, rng)
            return box, rng

        self._use_inc, self._inc = region("Included region", ms.use_included,
                                          ms.included_start, ms.included_end or 1)
        self._use_tgt, self._tgt = region("Target region", ms.use_target, ms.target_start,
                                          ms.target_end)
        self._use_prod, self._prod = region("Product size", ms.use_product_range,
                                            ms.product_min, ms.product_max)
        self._use_opt = QCheckBox("Optimal size", self)
        self._use_opt.setChecked(ms.use_optimal)
        self._opt = QSpinBox(self)
        self._opt.setRange(1, 1_000_000)
        self._opt.setValue(max(1, ms.optimal_size or 500))
        self._opt.setFixedWidth(110)
        self._opt.setSuffix(" bp")
        self._opt.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._opt.setEnabled(ms.use_optimal)
        self._use_opt.toggled.connect(self._opt.setEnabled)
        f.labelled(self._use_opt, self._opt, "Optimal product size",
                   "Sort products by how close they are to this size.")

        f.section("Mismatches")
        mm_row = QWidget(self)
        ml = QHBoxLayout(mm_row)
        ml.setContentsMargins(0, 0, 0, 0)
        ml.setSpacing(6)
        self._mm = QSpinBox(self)
        self._mm.setRange(0, 10)
        self._mm.setValue(ms.max_mismatches)
        self._mm.setFixedWidth(_SMALL)
        self._exact = QSpinBox(self)
        self._exact.setRange(0, 15)
        self._exact.setValue(ms.three_prime_exact)
        self._exact.setFixedWidth(_SMALL)
        ml.addWidget(self._mm)
        lab = QLabel("last", self)
        lab.setEnabled(False)
        ml.addWidget(lab)
        ml.addWidget(self._exact)
        lab2 = QLabel("3′ bases must match", self)
        lab2.setEnabled(False)
        ml.addWidget(lab2)
        mm_row.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        f.row("Maximum mismatches:", mm_row, "Mismatches",
              "Mismatches a primer may have at a binding site; mismatches in the last 3′ "
              "bases stop extension, so those must match.", fixed=False)

        self._thermo_fields = f.section("Tm Calculation")
        self._thermo = ThermoBox(self)
        self._thermo.set_settings(ms.thermo)
        self._thermo_fields += self._thermo.add_to(f)

        self._genome_fields = f.section("Genome")
        self._check_genome = QCheckBox("Also count binding sites and products in a genome", self)
        self._check_genome.setChecked(ms.check_genome)
        self._genome_fields.append(f.check(self._check_genome, "Genome", GENOME_HELP))
        self._genome = _path_row(self, "Genome FASTA: click it in Files, or Browse",
                                 "Choose the genome FASTA")
        g = ms.genome_fasta or (str(project_genome) if project_genome else "")
        self._genome.set_path(g or None, emit=False)
        self._genome.setEnabled(ms.check_genome)
        self._check_genome.toggled.connect(self._genome.setEnabled)
        self._genome_fields.append(f.row("Genome:", self._genome, fixed=False))
        lay.addWidget(body)
        self._form_main = f
        self._add_page("Map Primers", page)

        # Options tab
        page, lay = _page(self)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.setSpacing(14)
        rc_box = QGroupBox("Right-click on Map Primers", self)
        rcl = QVBoxLayout(rc_box)
        self._rc_run = QRadioButton("Map with saved settings", self)
        self._rc_dialog = QRadioButton("Open this window", self)
        rgrp = QButtonGroup(self)
        for b in (self._rc_run, self._rc_dialog):
            rgrp.addButton(b)
            rcl.addWidget(b)
        (self._rc_dialog if initial.map_right_click_action == PRIMER_RIGHT_CLICK_DIALOG
         else self._rc_run).setChecked(True)
        lay.addWidget(rc_box)
        res_box = QGroupBox("Results", self)
        rl = QVBoxLayout(res_box)
        self._annotate = QCheckBox("Add binding sites and products to the sequence annotations",
                                   self)
        self._annotate.setChecked(ms.annotate)
        rl.addWidget(self._annotate)
        db_row = QWidget(self)
        dbl = QHBoxLayout(db_row)
        dbl.setContentsMargins(0, 4, 0, 0)
        dbl.setSpacing(6)
        dbl.addWidget(QLabel("Primer Database:", self))
        self._db_path = _path_row(self, "Default (RepeatViewer settings folder)",
                                  "Choose the Primer Database FASTA")
        self._db_path.set_path(initial.database_path or None, emit=False)
        self._db_path.line_edit().editingFinished.connect(self._reload_db)
        dbl.addWidget(self._db_path, 1)
        rl.addWidget(db_row)
        process_box, self._show_process, self._close_process = make_process_output_box(
            self, tool_name="primer mapping", checked=initial.show_process_window,
            close_when_finished=initial.close_process_window_when_finished)
        lay.addWidget(res_box)
        lay.addWidget(_engine_box(self))
        lay.addWidget(process_box)
        self._add_page("Options", page)

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._adv_btn = QPushButton("Advanced Settings ▸", self)
        self._adv_btn.setCheckable(True)
        self._adv_btn.toggled.connect(self._toggle_advanced)
        self._finish(make_cancel_save_run_bar(self, left_widgets=[restore, self._adv_btn]))
        self._reload_db()
        self._toggle_advanced(False)
        self._sync()
        self._refresh_seq_label()

    def _database_obj(self) -> PrimerDatabase | None:
        return self._db

    def _reload_db(self) -> None:
        p = self._db_path.path()
        try:
            self._db = PrimerDatabase(p)
        except PrimerError:
            self._db = PrimerDatabase.__new__(PrimerDatabase)
            self._db.path = p
            self._db.primers = []
        for fld in (self._f, self._r, self._p):
            fld.reload()

    def _on_tab(self, index: int) -> None:
        if hasattr(self, "_adv_btn"):
            self._adv_btn.setVisible(index == self._TAB_MAP)
        super()._on_tab(index)

    def _toggle_advanced(self, on: bool) -> None:
        self._adv_btn.setText("Advanced Settings ▾" if on else "Advanced Settings ▸")
        for w in self._thermo_fields + self._genome_fields:
            self._form_main.set_visible(w, on)
        if self.isVisible():
            QTimer.singleShot(0, self._fit)

    def _sync(self, *_a) -> None:
        spec = self._specific.isChecked()
        for b, fld in ((self._use_f, self._f), (self._use_r, self._r), (self._use_p, self._p)):
            b.setEnabled(spec)
            fld.setEnabled(spec and b.isChecked())
        for w in (self._src, self._s_f, self._s_r, self._s_p, self._pairs_only):
            w.setEnabled(not spec)
        self._src_file.setEnabled(not spec and self._src.currentData() == "file")
        cons = self._test_on.currentData() == pj.TEST_ON_CONSENSUS
        self._form_main.set_visible(self._cons_field, cons)
        self._form_main.set_visible(self._cons_gaps_field, cons)
        if self.isVisible():
            QTimer.singleShot(0, self._fit)

    def _restore_defaults(self) -> None:
        d = pj.MapSettings()
        self._mm.setValue(d.max_mismatches)
        self._exact.setValue(d.three_prime_exact)
        self._thermo.set_settings(d.thermo)
        self._use_prod.setChecked(d.use_product_range)
        self._prod.set_values(d.product_min, d.product_max)
        self._use_inc.setChecked(False)
        self._use_tgt.setChecked(False)
        self._use_opt.setChecked(False)
        _set_combo(self._test_on, d.test_on)

    def _refresh_seq_label(self) -> None:
        if not self._files:
            self._seq_label.setText("<span style='color:palette(mid)'>The files selected in the "
                                    "Files panel</span>")
            return
        names = ", ".join(p.name for p in self._files[:3])
        more = f" +{len(self._files) - 3}" if len(self._files) > 3 else ""
        self._seq_label.setText(f"<b>{names}</b>{more}")

    def set_files(self, files: list[Path]) -> None:
        self._files = [Path(p) for p in files if _is_fasta(Path(p))]
        self._refresh_seq_label()

    def run_files(self) -> list[Path]:
        return list(self._files)

    def apply_path_from_panel(self, path: Path) -> bool:
        resolved = Path(path).expanduser()
        if not resolved.is_file():
            return False
        row = self._focused_path_row([self._src_file, self._genome, self._db_path])
        if row is not None:
            row.set_path(resolved)
            if row is self._db_path:
                self._reload_db()
            return True
        if not _is_fasta(resolved):
            return False
        self.set_files([resolved])
        return True

    def current_dialog_tab(self) -> str:
        return "map" if self._tab_bar.currentIndex() == self._TAB_MAP else DIALOG_TAB_OPTIONS

    def map_settings(self) -> pj.MapSettings:
        inc = self._inc.values()
        tgt = self._tgt.values()
        prod = self._prod.values()
        g = self._genome.path()
        src = self._src_file.path() if self._src.currentData() == "file" else None
        return pj.MapSettings(
            mode=pj.MAP_SPECIFIC if self._specific.isChecked() else pj.MAP_DATABASE,
            use_forward=self._use_f.isChecked(), forward=self._f.text(),
            use_reverse=self._use_r.isChecked(), reverse=self._r.text(),
            use_probe=self._use_p.isChecked(), probe=self._p.text(),
            primer_file=str(src) if src else "",
            search_forward=self._s_f.isChecked(), search_reverse=self._s_r.isChecked(),
            search_probe=self._s_p.isChecked(), pairs_only=self._pairs_only.isChecked(),
            use_included=self._use_inc.isChecked(), included_start=inc[0], included_end=inc[1],
            use_target=self._use_tgt.isChecked(), target_start=tgt[0], target_end=tgt[1],
            use_product_range=self._use_prod.isChecked(), product_min=prod[0],
            product_max=prod[1], use_optimal=self._use_opt.isChecked(),
            optimal_size=self._opt.value(), max_mismatches=self._mm.value(),
            three_prime_exact=self._exact.value(), thermo=self._thermo.settings(),
            test_on=str(self._test_on.currentData()),
            consensus_threshold=float(self._cons_threshold.currentData()),
            consensus_ignore_gaps=self._cons_gaps.isChecked(),
            check_genome=self._check_genome.isChecked(), genome_fasta=str(g) if g else "",
            annotate=self._annotate.isChecked(),
        )

    def options(self) -> PrimerOptions:
        db = self._db_path.path()
        return replace(
            self._opts, map=self.map_settings().to_dict(),
            database_path=str(db) if db else "",
            map_right_click_action=(PRIMER_RIGHT_CLICK_DIALOG if self._rc_dialog.isChecked()
                                    else PRIMER_RIGHT_CLICK_RUN),
            show_process_window=self._show_process.isChecked(),
            close_process_window_when_finished=self._close_process.isChecked(),
        )


# --------------------------------------------------------------------------
# Primer Database
# --------------------------------------------------------------------------


class PrimerDatabaseDialog(QDialog):
    """The saved primers: add, paste, import, export, delete, with live checks.

    Selecting one primer shows its Tm, GC, hairpin and self-dimer; selecting
    two shows the cross-dimer between them.
    """

    databaseChanged = Signal()

    COLS = ("Name", "Kind", "Sequence 5′→3′", "Length", "Tm °C", "GC %", "Note")

    def __init__(self, path: Path | None = None, *, thermo: ThermoSettings | None = None,
                 parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Primer Database")
        self.setModal(False)
        self.setWindowModality(Qt.WindowModality.NonModal)
        self.resize(940, 560)
        self._db = PrimerDatabase(path)
        self._thermo = thermo or ThermoSettings()
        self._stats_cache: dict[str, object] = {}

        self._filter = QLineEdit(self)
        self._filter.setPlaceholderText("Filter by name, sequence or note")
        self._filter.setClearButtonEnabled(True)
        self._filter.textChanged.connect(self._fill)
        self._path_label = QLabel(self)
        self._path_label.setEnabled(False)
        self._path_label.setText(str(self._db.path))
        self._path_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)

        self._table = QTableWidget(0, len(self.COLS), self)
        self._table.setHorizontalHeaderLabels(list(self.COLS))
        self._table.setSelectionBehavior(QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self._table.setEditTriggers(QAbstractItemView.EditTrigger.DoubleClicked
                                    | QAbstractItemView.EditTrigger.EditKeyPressed)
        self._table.verticalHeader().setVisible(False)
        self._table.verticalHeader().setDefaultSectionSize(22)
        self._table.setAlternatingRowColors(True)
        self._table.setSortingEnabled(True)
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        hh.setStretchLastSection(True)
        self._table.itemSelectionChanged.connect(self._show_details)
        self._table.itemChanged.connect(self._edited)
        self._table.horizontalHeader().setSortIndicator(0, Qt.SortOrder.AscendingOrder)

        self._details = QPlainTextEdit(self)
        self._details.setReadOnly(True)
        self._details.setFont(_mono_font())
        self._details.setPlaceholderText("Select a primer to see its Tm, GC, hairpin and "
                                         "self-dimer.\nSelect two to see the dimer between them.")
        split = QSplitter(Qt.Orientation.Vertical, self)
        split.addWidget(self._table)
        split.addWidget(self._details)
        split.setStretchFactor(0, 3)
        split.setStretchFactor(1, 2)

        add = QPushButton("Add…", self)
        add.setToolTip("Add primers: paste FASTA, name<TAB>sequence lines, or bare sequences")
        add.clicked.connect(self._add)
        imp = QPushButton("Import FASTA…", self)
        imp.clicked.connect(self._import)
        exp = QPushButton("Export…", self)
        exp.clicked.connect(self._export)
        self._delete = QPushButton("Delete", self)
        self._delete.clicked.connect(self._delete_selected)
        close = QPushButton("Close", self)
        close.setDefault(True)
        close.clicked.connect(self.accept)
        bar = QHBoxLayout()
        bar.setSpacing(8)
        for b in (add, imp, exp, self._delete):
            bar.addWidget(b)
        bar.addStretch(1)
        bar.addWidget(close)
        top = QHBoxLayout()
        top.addWidget(self._filter, 1)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.setSpacing(8)
        lay.addLayout(top)
        lay.addWidget(split, 1)
        lay.addWidget(self._path_label)
        lay.addLayout(bar)
        self._fill()

    # --------------------------------------------------------------- table

    def database(self) -> PrimerDatabase:
        return self._db

    def _stats(self, p: Primer):
        key = (p.sequence, self._thermo.to_dict().__repr__())
        if key not in self._stats_cache:
            try:
                self._stats_cache[key] = analyze_oligo(p.sequence, self._thermo)
            except Exception:  # noqa: BLE001
                self._stats_cache[key] = None
        return self._stats_cache[key]

    def _fill(self, *_a) -> None:
        text = self._filter.text().strip().lower()
        rows = [p for p in self._db.primers
                if not text or text in p.name.lower() or text in p.sequence.lower()
                or text in p.note.lower()]
        self._table.blockSignals(True)
        self._table.setSortingEnabled(False)
        self._table.setRowCount(len(rows))
        mono = _mono_font()
        for i, p in enumerate(rows):
            st = self._stats(p)
            vals = (p.name, KIND_LABELS.get(p.kind, p.kind), p.sequence, len(p.binding),
                    round(st.tm, 1) if st else "", round(st.gc) if st else "", p.note)
            for j, v in enumerate(vals):
                item = QTableWidgetItem()
                item.setData(Qt.ItemDataRole.DisplayRole, v)
                item.setData(Qt.ItemDataRole.UserRole, p.name)
                if j == 2:
                    item.setFont(mono)
                if j not in (0, 6):
                    item.setFlags(item.flags() & ~Qt.ItemFlag.ItemIsEditable)
                if j in (3, 4, 5):
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                          | Qt.AlignmentFlag.AlignVCenter)
                if st and st.warnings and j == 2:
                    item.setToolTip("; ".join(st.warnings))
                self._table.setItem(i, j, item)
        self._table.setSortingEnabled(True)
        self._table.sortByColumn(self._table.horizontalHeader().sortIndicatorSection(),
                                 self._table.horizontalHeader().sortIndicatorOrder())
        self._table.blockSignals(False)
        self._table.resizeColumnsToContents()
        self._table.setColumnWidth(2, max(220, min(380, self._table.columnWidth(2))))
        self._delete.setEnabled(bool(rows))
        self._show_details()

    def _selected(self) -> list[Primer]:
        names = []
        for idx in self._table.selectionModel().selectedRows():
            item = self._table.item(idx.row(), 0)
            if item is not None:
                names.append(item.data(Qt.ItemDataRole.UserRole))
        return [p for n in names if (p := self._db.get(n)) is not None]

    def _show_details(self) -> None:
        sel = self._selected()
        if not sel:
            self._details.clear()
            return
        lines: list[str] = []
        for p in sel[:2]:
            st = self._stats(p)
            if st is None:
                lines.append(f"{p.name}: cannot analyse {p.sequence}")
                continue
            lines.append(f"{p.name}  ({KIND_LABELS.get(p.kind, p.kind)})")
            lines.append(f"  5′-{p.sequence}-3′")
            lines.append(f"  {st.binding_length} nt binding"
                         + (f" + {len(p.tail)} nt tail" if p.tail else "")
                         + f" · Tm {st.tm:.1f} °C"
                         + (f" ({st.tm_min:.1f}–{st.tm_max:.1f})" if st.degeneracy > 1 else "")
                         + f" · GC {st.gc:.0f} % · 3′ ΔG {st.end_dg:.1f}")
            if st.warnings:
                lines.append("  ⚠ " + "; ".join(st.warnings))
            for s in (st.hairpin, st.self_dimer):
                if s.found:
                    tm = f", Tm {s.tm:.1f} °C" if s.tm > 0 else ""
                    lines.append(f"  {s.kind}: ΔG {s.dg:.1f} kcal/mol{tm}")
                    lines.extend("    " + x for x in s.lines)
                else:
                    lines.append(f"  {s.kind}: none")
            lines.append("")
        if len(sel) >= 2:
            d = heterodimer(sel[0].sequence, sel[1].sequence, self._thermo)
            lines.append(f"{sel[0].name} × {sel[1].name}: "
                         + (f"cross-dimer ΔG {d.dg:.1f} kcal/mol" if d.found else "no dimer"))
            lines.extend("    " + x for x in d.lines)
            ta, tb = self._stats(sel[0]), self._stats(sel[1])
            if ta and tb:
                lines.append(f"    Tm difference {abs(ta.tm - tb.tm):.1f} °C")
        self._details.setPlainText("\n".join(lines))

    def _edited(self, item: QTableWidgetItem) -> None:
        old = item.data(Qt.ItemDataRole.UserRole)
        p = self._db.get(old)
        if p is None:
            return
        if item.column() == 0:
            new = item.text().strip()
            if not new or (new != old and self._db.get(new) is not None):
                QTimer.singleShot(0, self._fill)
                return
            p.name = new
        elif item.column() == 6:
            p.note = item.text().strip()
        self._db.save()
        self.databaseChanged.emit()
        QTimer.singleShot(0, self._fill)

    # --------------------------------------------------------------- actions

    def add_primers(self, primers: list[Primer]) -> int:
        added = self._db.add(primers)
        if added:
            self._db.save()
            self.databaseChanged.emit()
            self._fill()
        return len(added)

    def _add(self) -> None:
        dlg = QDialog(self)
        dlg.setWindowTitle("Add Primers")
        dlg.resize(520, 300)
        text = QPlainTextEdit(dlg)
        text.setFont(_mono_font())
        text.setPlaceholderText(
            "CroSat2_F\tGACTTGCAAGTTCGAGAAC\nCroSat2_R\tGTTGCATCCACTTAGCTTAG\n\n"
            "or FASTA, or one sequence per line.\n"
            "Lower-case 5′ letters are a tail (e.g. caggaaacagctatgacGTTGCATCC…).")
        kind = QComboBox(dlg)
        kind.addItem("Kind from the name (…_F / …_R / …_P)", "")
        for k in PRIMER_KINDS:
            kind.addItem(KIND_LABELS[k], k)
        cancel = QPushButton("Cancel", dlg)
        save = QPushButton("Save", dlg)
        save.setDefault(True)
        cancel.clicked.connect(dlg.reject)
        save.clicked.connect(dlg.accept)
        bar = QHBoxLayout()
        bar.addWidget(kind)
        bar.addStretch(1)
        bar.addWidget(cancel)
        bar.addWidget(save)
        lay = QVBoxLayout(dlg)
        lay.addWidget(text, 1)
        lay.addLayout(bar)
        if dlg.exec() != int(QDialog.DialogCode.Accepted):
            return
        try:
            primers = parse_primer_text(text.toPlainText())
        except PrimerError as exc:
            QMessageBox.warning(self, "Add Primers", str(exc))
            return
        forced = kind.currentData()
        if forced:
            for p in primers:
                p.kind = forced
        n = self.add_primers(primers)
        if n < len(primers):
            QMessageBox.information(self, "Add Primers",
                                    f"Added {n} of {len(primers)}; the others were already "
                                    "in the database (same sequence).")

    def _import(self) -> None:
        path, _ = QFileDialog.getOpenFileName(self, "Import primers", str(Path.home()),
                                              "Primers (*.fa *.fasta *.txt *.tsv);;All files (*)")
        if not path:
            return
        from repeatviewer.backends.primers.oligos import load_primers

        try:
            primers = load_primers(Path(path))
        except PrimerError as exc:
            QMessageBox.warning(self, "Import primers", str(exc))
            return
        n = self.add_primers(primers)
        QMessageBox.information(self, "Import primers", f"Added {n} primer(s) from "
                                                        f"{Path(path).name}.")

    def _export(self) -> None:
        sel = self._selected() or list(self._db.primers)
        if not sel:
            return
        path, _ = QFileDialog.getSaveFileName(self, "Export primers",
                                              str(Path.home() / "primers.fa"),
                                              "FASTA (*.fa);;Order sheet (*.tsv)")
        if not path:
            return
        p = Path(path)
        if p.suffix.lower() == ".tsv":
            lines = ["name\tsequence\tlength\ttm\tgc\tnote"]
            for pr in sel:
                st = self._stats(pr)
                lines.append(f"{pr.name}\t{pr.sequence}\t{len(pr.sequence)}\t"
                             f"{st.tm:.1f}\t{st.gc:.0f}\t{pr.note}" if st else
                             f"{pr.name}\t{pr.sequence}\t{len(pr.sequence)}\t\t\t{pr.note}")
            p.write_text("\n".join(lines) + "\n", encoding="utf-8")
        else:
            write_primers(p, sel)

    def _delete_selected(self) -> None:
        sel = self._selected()
        if not sel:
            return
        names = ", ".join(p.name for p in sel[:5]) + (" …" if len(sel) > 5 else "")
        ok = QMessageBox.question(self, "Delete primers",
                                  f"Delete {len(sel)} primer(s) from the database?\n{names}")
        if ok != QMessageBox.StandardButton.Yes:
            return
        self._db.remove([p.name for p in sel])
        self._db.save()
        self.databaseChanged.emit()
        self._fill()


def primer3_missing_message(parent: QWidget | None) -> bool:
    """Show the install hint when primer3-py is missing; True when it is."""
    if primer3_available():
        return False
    QMessageBox.warning(parent, "Primers", PRIMER3_MISSING)
    return True


__all__ = [
    "MapPrimersDialog",
    "PrimerDatabaseDialog",
    "PrimerDesignDialog",
    "PrimerPickerDialog",
    "primer3_missing_message",
]
