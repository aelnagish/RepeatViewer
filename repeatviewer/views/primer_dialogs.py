"""Primer windows (Tools → Primers).

* :class:`PrimerDesignDialog`  – Primers | Regions | Specificity | Options
* :class:`MapPrimersDialog`    – Map Primers | Options (Geneious' Map Primers)
* :class:`PrimerDatabaseDialog` – the saved primers, with Tm / hairpin / dimer checks
* :class:`PrimerPickerDialog`  – choose a saved primer for a field

All follow the tool-window rules: Cancel / Save / Run, an Options tab, file
fields filled by a click in the Files panel, "?" help next to the item.

Developer: Abdullah Yusuf
"""

from __future__ import annotations

from dataclasses import replace
from pathlib import Path

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QFont, QFontDatabase, QStandardItem, QStandardItemModel
from PySide6.QtWidgets import (
    QAbstractItemView,
    QButtonGroup,
    QCheckBox,
    QDialog,
    QFormLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QMessageBox,
    QPushButton,
    QRadioButton,
    QSizePolicy,
    QStackedWidget,
    QScrollArea,
    QFrame,
    QApplication,
    QStyledItemDelegate,
    QStyleOptionViewItem,
    QTabBar,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.primers import design as dz
from repeatviewer.backends.primers import jobs as pj
from repeatviewer.backends.primers import painting as pt
from repeatviewer.backends.primers.oligos import (
    KIND_LABELS,
    Primer,
    PrimerDatabase,
)
from repeatviewer.backends.primers.thermo import (
    SALT_METHOD_LABELS,
    SALT_METHODS,
    TM_METHOD_LABELS,
    TM_METHODS,
    PrimerError,
    ThermoSettings,
    primer3_available,
    PRIMER3_MISSING,
)
from repeatviewer.job_options import (
    DIALOG_TAB_OPTIONS,
    PRIMER_DIALOG_TAB_CHORUS2,
    PRIMER_DIALOG_TAB_PRIMERS,
    PRIMER_DIALOG_TAB_SELECTION,
    PRIMER_DIALOG_TABS,
    PRIMER_RIGHT_CLICK_RUN,
    PrimerOptions,
    opening_dialog_tab,
)
from repeatviewer.views.tool_dialogs import (
    TOOL_PARAMS_DIALOG_WIDTH_WIDE,
    _form_section_separator,
    settings_form_label_column_width,
    dialog_should_run,
    make_cancel_save_run_bar,
    make_process_output_box,
    sync_form_field_column_spacers,
)
from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox
from repeatviewer.widgets.path_browse import PathBrowseRow
from repeatviewer.widgets.spin_boxes import NoWheelDoubleSpinBox as QDoubleSpinBox
from repeatviewer.widgets.spin_boxes import NoWheelSpinBox as QSpinBox

_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".fsa", ".ffn", ".seq", ".aln", ".mfa"}
_FASTA_FILTER = "FASTA (*.fa *.fasta *.fna *.fas *.fa.gz *.fasta.gz);;All files (*)"
# The form grid (same in Design and Map): every number box is BOX_W wide; a
# Min / Opt / Max row is box, word, box, word, box; a two-value row (box "to"
# box) uses the first two of those columns; everything else (lists, paths,
# text) is FIELD_W wide, so all rows end on one line, and every "?" sits in
# one column after it (rows without help keep the space).
BOX_W = 96
JOIN_W = 34
GRID_GAP = 6
FIELD_W = 3 * BOX_W + 2 * JOIN_W + 4 * GRID_GAP

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


_INDENT_ROLE = Qt.ItemDataRole.UserRole + 7


class _IndentDelegate(QStyledItemDelegate):
    """Indents the popup rows of the grouped Mode list (the closed box shows no indent)."""

    _STEP = 14

    def paint(self, painter, option, index) -> None:
        level = int(index.data(_INDENT_ROLE) or 0)
        if level:
            moved = QStyleOptionViewItem(option)
            moved.rect = option.rect.adjusted(self._STEP * level, 0, 0, 0)
            option = moved
        super().paint(painter, option, index)

    def sizeHint(self, option, index):  # noqa: N802
        size = super().sizeHint(option, index)
        size.setWidth(size.width() + self._STEP * int(index.data(_INDENT_ROLE) or 0))
        return size


def _fill_grouped_modes(combo) -> None:
    """Mode list grouped as Family > Group > modes; the headers cannot be chosen."""
    model = QStandardItemModel(combo)
    for family, groups in dz.MODE_FAMILIES:
        head = QStandardItem(family)
        font = QFont(combo.font())
        font.setBold(True)
        head.setFont(font)
        head.setFlags(Qt.ItemFlag.NoItemFlags)
        head.setData(0, _INDENT_ROLE)
        model.appendRow(head)
        for group, modes in groups:
            sub = QStandardItem(group)
            sfont = QFont(combo.font())
            sfont.setItalic(True)
            sub.setFont(sfont)
            sub.setFlags(Qt.ItemFlag.NoItemFlags)
            sub.setData(1, _INDENT_ROLE)
            model.appendRow(sub)
            for mode in modes:
                item = QStandardItem(dz.MODE_LABELS[mode])
                item.setData(mode, Qt.ItemDataRole.UserRole)
                item.setData(2, _INDENT_ROLE)
                item.setToolTip(f"{family} / {group}")
                model.appendRow(item)
    combo.setModel(model)
    combo.setItemDelegate(_IndentDelegate(combo))
    first = next(i for i in range(model.rowCount())
                 if model.item(i).data(Qt.ItemDataRole.UserRole))
    combo.setCurrentIndex(first)


# --------------------------------------------------------------------------
# Help texts
# --------------------------------------------------------------------------

MODE_HELP = (
    "<b>PCR primers</b><br>"
    "<b>Standard PCR</b>: a primer pair by Primer3 on the template.<br>"
    "<b>qPCR</b>: short amplicons, matched Tm, optional hydrolysis probe; add a "
    "single-copy reference gene to get a reference pair as well.<br>"
    "<b>Ladder PCR</b> (satellites): outward-facing primers in one monomer. A tandem array "
    "then gives a ladder of 1×, 2×, 3× … monomer bands.<br><br>"
    "<b>FISH probes</b><br>"
    "<b>Satellite repeat oligos</b>: 30–45 nt oligos found in every monomer / copy, without "
    "hairpins, ready to order with a 5′ label.<br>"
    "<b>Chromosome painting</b>: single-copy oligo pools along whole chromosomes, made by "
    "Chorus2."
)
INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): design on the FASTA selected in the Files "
    "panel. Clicking another file in the Files panel changes it.<br>"
    "<b>Define inputs</b>: choose the template FASTA for this run here. Click a file in the "
    "Files panel to fill the field, or use Browse."
)
PAINT_SOURCE_HELP = (
    "<b>Run Chorus2 here</b>: RepeatViewer runs Chorus2 from its conda env on this "
    "computer. A plant genome of a few Gb takes hours and needs a lot of memory.<br>"
    "<b>Import results from the server</b>: run Chorus2 on a Linux server (Copy Server "
    "Command), copy back its BED and choose it here. Selection and the order file are then "
    "made here."
)
PAINT_GENOME_HELP = (
    "The whole genome FASTA. Chorus2 checks every probe against it, so probes are unique in "
    "the genome. When importing, it is optional and only gives the chromosome lengths."
)
PAINT_INPUT_HELP = (
    "Leave empty to design along the whole genome. Give one chromosome or region (FASTA, "
    "with the same sequence names as the genome) to paint only that part."
)
PAINT_BED_HELP = (
    "The probe BED written by Chorus2 (<code>chorus2_out/&lt;input&gt;.bed</code>) or by "
    "ChorusNGSfilter (with a k-mer score column)."
)
PAINT_NGS_HELP = (
    "Optional. Genomic reads (WGS, a few x coverage) give each probe a k-mer score. Probes "
    "missing from the reads (deleted or mis-assembled) or too frequent (repeats missed by "
    "the assembly) are dropped. Needs jellyfish.<br>"
    "<b>Import</b>: tick it to use the score column of a BED made by ChorusNGSfilter; the "
    "quantile window is applied here."
)
PAINT_SERVER_HELP = (
    "Copies the Chorus2 command lines with these settings. On the server: put the genome "
    "(and reads) in one folder, <code>conda activate chorus2</code>, paste. Then copy "
    "<code>chorus2_out/*.bed</code> back and use <b>Import results from the server</b>."
)
PAINT_AMOUNT_HELP = (
    "<b>All probes</b>: every probe that passes the spacing.<br>"
    "<b>Density</b>: about this many probes per kb, evenly along each chromosome (painting "
    "usually needs 0.1 to 1 per kb).<br>"
    "<b>Pool size</b>: about this many probes in total, spread evenly (for example the size "
    "of an oligo pool you order)."
)
PAINT_REGIONS_HELP = (
    "Keep probes only here. Separate items with a comma and a space or with a semicolon: "
    "<code>chr1, chr2:1-5000000</code>. A name alone means the whole sequence. Empty = "
    "everything."
)
PAINT_ADAPTER_HELP = (
    "Added to every oligo in the order file, for amplifying the pool (for example a T7 "
    "promoter or the Chorus2 5′ primer). The BED and FASTA keep the bare probe."
)
MAP_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): test the primers on the FASTA files selected "
    "in the Files panel. Clicking another file changes it.<br>"
    "<b>Define inputs</b>: choose the FASTA for this run here. Click a file in the Files "
    "panel to fill the field, or use Browse."
)
TEMPLATE_HELP = (
    "The FASTA the primers are designed on. Every sequence in the file is a template (up to "
    "20).<br><b>Ladder / Satellite oligos</b>: a set of monomers or copies, aligned or not."
)
COVERAGE_HELP = (
    "A primer counts as binding a sequence when it has at most <b>Mismatches tolerated</b> "
    "mismatches and its last three 3′ bases match. <b>Min coverage</b> is the share of "
    "sequences a probe must bind."
)
GENOME_HELP = (
    "Every primer is also searched in this FASTA (both strands, with the mismatches above) "
    "and the products each pair would make are counted, like Primer-BLAST but local."
)
TAILS_HELP = (
    "A 5′ tail (adapter, promoter, restriction site with extra bases) added to the primer. "
    "Tails are written in lower case and do not count for the template Tm."
)
TM_HELP = (
    "Thermodynamic model and reaction conditions for every Tm and ΔG (Primer3)."
)
PRODUCT_HELP = (
    "Size range of the PCR product. For ladder PCR this is the smallest band."
)


# --------------------------------------------------------------------------
# Small building blocks
# --------------------------------------------------------------------------


class RangeEdit(QWidget):
    """Min [–] opt [–] max number boxes on the form grid (or Min [–] to [–] Max)."""

    changed = Signal()

    def __init__(self, lo: float, hi: float, *, decimals: int = 0, opt: bool = True,
                 step: float = 1.0, suffix: str = "", parent: QWidget | None = None) -> None:
        super().__init__(parent)
        row = QHBoxLayout(self)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(GRID_GAP)
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
            box.setFixedWidth(BOX_W)
            box.setToolTip(lab)
            box.valueChanged.connect(lambda _v: self.changed.emit())
            if i:
                word = QLabel(lab.lower() if opt else "to", self)
                word.setEnabled(False)
                word.setFixedWidth(JOIN_W)
                word.setAlignment(Qt.AlignmentFlag.AlignCenter)
                row.addWidget(word)
            row.addWidget(box)
            self._boxes.append(box)
        n = len(labels)
        self.setFixedWidth(n * BOX_W + (n - 1) * (JOIN_W + 2 * GRID_GAP))
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
            return b

        self.mono = conc(0, 2000, 5, " mM")
        self.di = conc(0, 200, 0.5, " mM", 2)
        self.dntp = conc(0, 50, 0.1, " mM", 2)
        self.oligo = conc(0.1, 100000, 10, " nM")
        self.set_settings(ThermoSettings())

    def add_to(self, form: "_Form") -> list[QWidget]:
        return [form.row("Formula:", self.formula, "Tm calculation", TM_HELP),
                form.row("Salt correction:", self.salt),
                form.row("Monovalent:", self.mono),
                form.row("Divalent:", self.di),
                form.row("dNTPs:", self.dntp),
                form.row("Oligo:", self.oligo)]

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


class TailEdit:
    """Preset list + editable tail sequence, as two rows of the form grid."""

    def __init__(self, parent: QWidget) -> None:
        self.combo = QComboBox(parent)
        for label, seq in TAIL_PRESETS:
            self.combo.addItem(label, seq)
        self.edit = QLineEdit(parent)
        self.edit.setPlaceholderText("5′ tail (lower case)")
        self.edit.setFont(_mono_font())
        self.edit.setClearButtonEnabled(True)
        self.combo.currentIndexChanged.connect(self._preset)
        self.edit.textEdited.connect(self._typed)

    def add_to(self, form: "_Form", label: str, help_title: str | None = None,
               help_text: str | None = None) -> list[QWidget]:
        return [form.row(label, self.combo, help_title, help_text),
                form.row("", self.edit)]

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
    """Form on the grid: right-aligned labels, fields of BOX_W or FIELD_W, one "?" column.

    Every field sits in a cell of the same width, so the right edges and the
    "?" buttons line up on every row of every tab.
    """

    def __init__(self, host: QWidget) -> None:
        self.host = host
        self.form = QFormLayout(host)
        self.form.setContentsMargins(0, 0, 0, 0)
        self.form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
        self.form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        self.form.setHorizontalSpacing(12)
        self.form.setVerticalSpacing(6)
        self.spacers: list[QWidget] = []
        self._first = True

    @staticmethod
    def _fit(widget: QWidget) -> None:
        if isinstance(widget, (QSpinBox, QDoubleSpinBox)):
            widget.setFixedWidth(BOX_W)
        elif isinstance(widget, RangeEdit):
            pass
        else:
            widget.setFixedWidth(FIELD_W)

    def _cell(self, widget: QWidget, help_title: str | None, help_text: str | None) -> QWidget:
        from repeatviewer.widgets.help_button import make_help_button

        cell = QWidget(self.host)
        row = QHBoxLayout(cell)
        row.setContentsMargins(0, 0, 0, 0)
        row.setSpacing(GRID_GAP)
        area = QWidget(cell)
        area.setFixedWidth(FIELD_W)
        al = QHBoxLayout(area)
        al.setContentsMargins(0, 0, 0, 0)
        al.setSpacing(0)
        al.addWidget(widget, 0, Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter)
        al.addStretch(1)
        row.addWidget(area)
        button = make_help_button(cell, help_title or "", help_text or "")
        if not help_text:
            # Keeps the "?" column: an invisible button that still takes its space.
            policy = button.sizePolicy()
            policy.setRetainSizeWhenHidden(True)
            button.setSizePolicy(policy)
            button.hide()
        row.addWidget(button, 0, Qt.AlignmentFlag.AlignVCenter)
        cell.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        cell.setProperty("gridField", widget)
        return cell

    def _label(self, text: str) -> QWidget:
        lab = QLabel(text, self.host)
        if not text:
            self.spacers.append(lab)
        return lab

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
        del fixed   # every field follows the grid
        self._fit(widget)
        cell = self._cell(widget, help_title or label.rstrip(":"), help_text)
        self.form.addRow(self._label(label), cell)
        return cell

    def check(self, box: QWidget, help_title: str | None = None,
              help_text: str | None = None) -> QWidget:
        cell = self._cell(box, help_title or "", help_text)
        self.form.addRow(self._label(""), cell)
        return cell

    def switch_row(self, box: QCheckBox, widget: QWidget, help_title: str | None = None,
                   help_text: str | None = None) -> QWidget:
        """A row whose label is a checkbox (``☐ Included region:``) turning the field on.

        The checkbox (text right of the box, rule 16) is right-aligned in the label
        column, so its field starts where every other field starts."""
        host = QWidget(self.host)
        hl = QHBoxLayout(host)
        hl.setContentsMargins(0, 0, 0, 0)
        hl.setSpacing(0)
        hl.addStretch(1)
        hl.addWidget(box)
        self.spacers.append(host)
        self._fit(widget)
        cell = self._cell(widget, help_title or box.text().rstrip(":"), help_text)
        self.form.addRow(host, cell)
        return cell

    def widget_row(self, widget: QWidget) -> QWidget:
        widget.setFixedWidth(FIELD_W)
        cell = self._cell(widget, None, None)
        self.form.addRow(self._label(""), cell)
        return cell

    def set_visible(self, field: QWidget, visible: bool) -> None:
        self.form.setRowVisible(field, visible)


def _path_row(parent: QWidget, placeholder: str, title: str, *, mode: str = "file",
              file_filter: str = _FASTA_FILTER) -> PathBrowseRow:
    row = PathBrowseRow(parent, placeholder=placeholder, browse_title=title, mode=mode,
                        file_filter=file_filter)
    row.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
    row.line_edit().setMinimumWidth(0)
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
        # Each page scrolls when it is taller than the screen (Advanced Settings
        # on a laptop), so the buttons always stay visible.
        scroll = QScrollArea(self)
        scroll.setWidgetResizable(True)
        scroll.setFrameShape(QFrame.Shape.NoFrame)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        scroll.setWidget(page)
        scroll.viewport().setAutoFillBackground(False)
        page.setAutoFillBackground(False)
        self._tab_bar.addTab(title)
        self._stack.addWidget(scroll)
        return self._tab_bar.count() - 1

    def _finish(self, bar: QWidget) -> None:
        self._layout.addWidget(bar, 0)
        self._width = self._grid_width()
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

    def _grid_width(self) -> int:
        """Window width that fits the form grid exactly: label column, field, "?"."""
        from repeatviewer.widgets.help_button import make_help_button

        probe = make_help_button(None, "", "")
        help_w = probe.sizeHint().width()
        probe.deleteLater()
        label_w = settings_form_label_column_width(self)
        margins = self._layout.contentsMargins()
        # room for the page's scroll bar, so a long page never covers the "?" column
        bar = self.style().pixelMetric(self.style().PixelMetric.PM_ScrollBarExtent)
        return (label_w + 12 + FIELD_W + GRID_GAP + help_w + 2 * 16 + bar
                + margins.left() + margins.right())

    def _fit(self) -> None:
        self._sync_columns()
        self._width = self._grid_width()
        scroll = self._stack.currentWidget()
        if scroll is None:
            return
        page = scroll.widget() if isinstance(scroll, QScrollArea) else scroll
        page.adjustSize()
        wanted = page.sizeHint().height()
        # Room on this screen for the page: screen height minus the title bar,
        # tab bar, button bar and margins.
        screen = self.screen() or QApplication.primaryScreen()
        avail = screen.availableGeometry() if screen is not None else None
        if avail is not None:
            chrome = (self._tab_bar.sizeHint().height() + 2 * self._layout.spacing()
                      + self._layout.contentsMargins().top()
                      + self._layout.contentsMargins().bottom() + 44 + 40)
            room = max(240, avail.height() - chrome)
        else:
            room = wanted
        self._stack.setFixedHeight(min(wanted, room))
        self.adjustSize()
        self.resize(self._width, self.sizeHint().height())
        if avail is not None:
            # Keep the whole window, buttons included, on the screen.
            frame = self.frameGeometry()
            top = frame.top()
            if frame.bottom() > avail.bottom():
                top = max(avail.top(), avail.bottom() - frame.height())
            if top != frame.top():
                self.move(frame.left(), top)

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

    _TAB_PRIMERS, _TAB_CHORUS2, _TAB_SELECTION, _TAB_OPTIONS = range(4)
    _TAB_NAMES = (PRIMER_DIALOG_TAB_PRIMERS, PRIMER_DIALOG_TAB_CHORUS2,
                  PRIMER_DIALOG_TAB_SELECTION, DIALOG_TAB_OPTIONS)

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
        self._input_group = QButtonGroup(self)
        self._input_selected = QRadioButton("Use selected sequences", self)
        self._input_define = QRadioButton("Define inputs", self)
        self._input_group.addButton(self._input_selected)
        self._input_group.addButton(self._input_define)
        (self._input_selected if self._templates else self._input_define).setChecked(True)
        f.check(self._input_selected, "Inputs", INPUTS_HELP)
        self._template_label = QLabel(self)
        self._template_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._template_label.setWordWrap(True)
        self._template_label.setIndent(22)  # under the radio text, not its circle
        self._template_label_field = f.widget_row(self._template_label)
        f.check(self._input_define)
        self._template_row = _path_row(self, "Template FASTA: click it in Files, or Browse",
                                       "Choose the template FASTA")
        self._template_field = f.row("Template:", self._template_row, "Template", TEMPLATE_HELP,
                                     fixed=False)
        self._input_selected.toggled.connect(self._sync_inputs)
        self._template_row.path_changed.connect(lambda *_a: self._on_template_changed())

        f.section("Settings")
        self._mode = QComboBox(self)
        _fill_grouped_modes(self._mode)
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
             "as the reference for relative quantification (ΔCt).")
        self._min_cov = QSpinBox(self)
        self._min_cov.setRange(10, 100)
        self._min_cov.setSuffix(" %")
        self._min_cov.setValue(int(round(s.min_coverage * 100)))
        mrow((dz.MODE_FISH,), "Min coverage:", self._min_cov, "Coverage", COVERAGE_HELP)
        self._cov_mm = QSpinBox(self)
        self._cov_mm.setRange(0, 8)
        self._cov_mm.setValue(s.coverage_mismatches)
        mrow((dz.MODE_LADDER, dz.MODE_FISH), "Mismatches tolerated:", self._cov_mm, "Coverage",
             COVERAGE_HELP)
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
             "ranked last (e.g. the other satellite families).")

        self._conditions_title = f.section("Primer Conditions")[-1]
        self._size = RangeEdit(5, 60)
        self._size.set_values(s.size_min, s.size_opt, s.size_max)
        self._size_field = f.row("Primer size:", self._size, "Primer size",
                                 "Length in bases (min, optimum, max).")
        self._tm = RangeEdit(20, 95, decimals=1, step=0.5)
        self._tm.set_values(s.tm_min, s.tm_opt, s.tm_max)
        f.row("Tm (°C):", self._tm, "Tm", "Melting temperature of the binding part.")
        self._gc = RangeEdit(0, 100, decimals=0, opt=False, suffix=" %")
        self._gc.set_values(s.gc_min, s.gc_max)
        f.row("GC:", self._gc)
        self._product = RangeEdit(10, 100_000, opt=False, suffix=" bp")
        self._product.set_values(s.product_min, s.product_max)
        self._product_field = f.row("Product size:", self._product, "Product size", PRODUCT_HELP)
        self._tm_diff = QDoubleSpinBox(self)
        self._tm_diff.setRange(0, 20)
        self._tm_diff.setDecimals(1)
        self._tm_diff.setSuffix(" °C")
        self._tm_diff.setValue(s.max_tm_diff)
        self._tm_diff_field = f.row("Max Tm difference:", self._tm_diff)
        self._num = QSpinBox(self)
        self._num.setRange(1, 100)
        self._num.setValue(s.num_return)
        self._num_field = f.row("Pairs to return:", self._num)

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
        self._adv_fields.append(f.row("Probe size:", self._probe_size))
        self._probe_tm = RangeEdit(30, 95, decimals=1, step=0.5)
        self._probe_tm.set_values(s.probe_tm_min, s.probe_tm_opt, s.probe_tm_max)
        self._adv_fields.append(f.row("Probe Tm (°C):", self._probe_tm))
        self._fwd_tail = TailEdit(self)
        self._fwd_tail.set_tail(s.forward_tail)
        self._adv_fields += self._fwd_tail.add_to(f, "Forward 5′ tail:", "5′ tails", TAILS_HELP)
        self._rev_tail = TailEdit(self)
        self._rev_tail.set_tail(s.reverse_tail)
        self._adv_fields += self._rev_tail.add_to(f, "Reverse 5′ tail:")
        self._adv_fields += f.section("Tm Calculation")
        self._thermo = ThermoBox(self)
        self._thermo.set_settings(s.thermo)
        self._adv_fields += self._thermo.add_to(f)
        lay.addWidget(body, 0, Qt.AlignmentFlag.AlignHCenter)
        self._primers_form = f
        self._add_page("Primers", page)

        self._build_painting_pages(initial.painting_settings(), s, project_genome)

        # ---------------- Options tab
        page, lay = _page(self)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.setSpacing(14)
        rc_box = QGroupBox("Right-click on Primers", self)
        rcl = QVBoxLayout(rc_box)
        self._rc_run = QRadioButton("Design with saved settings", self)
        self._rc_run.setChecked(True)
        grp = QButtonGroup(self)
        grp.addButton(self._rc_run)
        rcl.addWidget(self._rc_run)
        res_box = QGroupBox("Results", self)
        rl = QVBoxLayout(res_box)
        self._annotate = QCheckBox("Ask which pairs to add to the sequence after the run", self)
        self._annotate.setChecked(initial.annotate_design)
        self._to_db = QCheckBox("Add the designed primers to the Primer Database", self)
        self._to_db.setChecked(initial.save_design_to_database)
        rl.addWidget(self._annotate)
        rl.addWidget(self._to_db)
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
        self._paint_mode.currentIndexChanged.connect(self._on_paint_mode)
        self._apply_tab_family(select=False)
        open_tab = opening_dialog_tab(initial.last_dialog_tab, allowed=PRIMER_DIALOG_TABS,
                                      first_tab=PRIMER_DIALOG_TAB_PRIMERS)
        index = self._TAB_NAMES.index(open_tab) if open_tab in self._TAB_NAMES else 0
        if not self._tab_bar.isTabVisible(index):
            index = self._TAB_CHORUS2 if self._painting() else self._TAB_PRIMERS
        self._tab_bar.setCurrentIndex(index)
        self._apply_mode_rows()
        self._toggle_advanced(False)
        self._sync_inputs()

    # ------------------------------------------------------------- behaviour

    def _sync_inputs(self, *_a) -> None:
        define = self._input_define.isChecked()
        self._primers_form.set_visible(self._template_label_field, not define)
        self._primers_form.set_visible(self._template_field, define)
        self._refresh_template_label()
        if self.isVisible():
            QTimer.singleShot(0, self._fit)

    def _on_template_changed(self) -> None:
        pass

    def _all_path_rows(self) -> list[PathBrowseRow]:
        return [self._template_row, self._reference, self._fish_avoid, self._pt_genome,
                self._pt_input, self._pt_bed, self._pt_fastq]

    def _on_tab(self, index: int) -> None:
        if hasattr(self, "_adv_btn"):
            self._adv_btn.setVisible(index == self._TAB_PRIMERS)
        super()._on_tab(index)

    # ------------------------------------------------------------- painting

    def _painting(self) -> bool:
        return self._mode.currentData() == dz.MODE_PAINTING

    def is_painting(self) -> bool:
        """Chromosome painting (Chorus2) is chosen: Run starts a painting job."""
        return self._painting()

    def _apply_tab_family(self, *, select: bool) -> None:
        """Painting shows Chorus2 | Selection | Options; the rest Primers | Options."""
        paint = self._painting()
        was = self._tab_bar.currentIndex()
        first = self._TAB_CHORUS2 if paint else self._TAB_PRIMERS
        # Show the new family first, so hiding the old one never jumps to Options.
        self._tab_bar.setTabVisible(first, True)
        if was != self._TAB_OPTIONS or not select:
            self._tab_bar.setCurrentIndex(first)
        self._tab_bar.setTabVisible(self._TAB_PRIMERS, not paint)
        for i in (self._TAB_CHORUS2, self._TAB_SELECTION):
            self._tab_bar.setTabVisible(i, paint)
        if was == self._TAB_OPTIONS and select:
            self._tab_bar.setCurrentIndex(self._TAB_OPTIONS)

    def _on_paint_mode(self, index: int) -> None:
        if self._mode.currentIndex() != index:
            self._mode.setCurrentIndex(index)

    def _build_painting_pages(self, ps: pt.PaintingSettings, s: dz.DesignSettings,
                              project_genome: Path | None) -> None:
        # ---------------- Chorus2 tab
        page, lay = _page(self)
        body = QWidget(page)
        c = _Form(body)
        self._forms.append(c)
        c.section("Inputs")
        grp = QButtonGroup(self)
        self._pt_run = QRadioButton("Run Chorus2 here", self)
        self._pt_import = QRadioButton("Import results from the server", self)
        grp.addButton(self._pt_run)
        grp.addButton(self._pt_import)
        (self._pt_import if ps.source == pt.SOURCE_IMPORT else self._pt_run).setChecked(True)
        c.check(self._pt_run, "Where Chorus2 runs", PAINT_SOURCE_HELP)
        c.check(self._pt_import)
        self._pt_genome = _path_row(self, "Genome FASTA: click it in Files, or Browse",
                                    "Choose the genome FASTA")
        genome = ps.genome_fasta or (str(project_genome) if project_genome else "")
        self._pt_genome.set_path(genome or None, emit=False)
        self._pt_genome_field = c.row("Genome:", self._pt_genome, "Genome", PAINT_GENOME_HELP,
                                      fixed=False)
        self._pt_input = _path_row(self, "Optional: one chromosome or region (FASTA)",
                                   "Choose the chromosome or region FASTA")
        self._pt_input.set_path(ps.input_fasta or None, emit=False)
        self._pt_input_field = c.row("Paint only:", self._pt_input, "Paint only",
                                     PAINT_INPUT_HELP, fixed=False)
        self._pt_bed = _path_row(self, "Chorus2 BED: click it in Files, or Browse",
                                 "Choose the Chorus2 probe BED",
                                 file_filter="BED (*.bed);;All files (*)")
        self._pt_bed.set_path(ps.import_bed or None, emit=False)
        self._pt_bed_field = c.row("Chorus2 BED:", self._pt_bed, "Chorus2 BED", PAINT_BED_HELP,
                                   fixed=False)
        c.section("Settings")
        self._paint_mode = QComboBox(self)
        self._paint_mode.setModel(self._mode.model())
        self._paint_mode.setItemDelegate(_IndentDelegate(self._paint_mode))
        self._paint_mode.setCurrentIndex(self._mode.currentIndex())
        c.row("Mode:", self._paint_mode, "Design mode", MODE_HELP)

        self._pt_chorus_fields = c.section("Chorus2")

        def spin(lo, hi, value, suffix=""):
            b = QSpinBox(self)
            b.setRange(lo, hi)
            b.setValue(int(value))
            if suffix:
                b.setSuffix(suffix)
            return b

        def dspin(lo, hi, value, suffix="", step=1.0):
            b = QDoubleSpinBox(self)
            b.setRange(lo, hi)
            b.setDecimals(0 if step >= 1 else 2)
            b.setSingleStep(step)
            b.setValue(float(value))
            if suffix:
                b.setSuffix(suffix)
            return b

        self._pt_len = spin(30, 100, ps.length, " nt")
        self._pt_hom = dspin(50, 100, ps.homology, " %")
        self._pt_dtm = dspin(0, 37, ps.dtm, " °C")
        self._pt_step = spin(1, 100, ps.step, " bp")
        self._pt_ploidy = spin(1, 12, ps.ploidy)
        self._pt_threads = spin(1, 256, ps.threads)
        for label, w, title, text in (
                ("Probe length:", self._pt_len, "Probe length",
                 "Oligo length (Chorus2 <code>-l</code>, default 45 nt)."),
                ("Max homology:", self._pt_hom, "Homology",
                 "A probe is dropped when another genome position matches more than this "
                 "share of it (<code>--homology</code>, default 75 %). Lower = more unique."),
                ("Min dTm:", self._pt_dtm, "dTm",
                 "Hybrid Tm minus hairpin Tm must be at least this (<code>-d</code>, default "
                 "10 °C)."),
                ("k-mer step:", self._pt_step, "Step",
                 "Sliding-window step of the k-mer search (<code>--step</code>, default 5 bp). "
                 "Smaller = more candidates, slower."),
                ("Ploidy:", self._pt_ploidy, "Ploidy",
                 "Of the genome (<code>--ploidy</code>, default 2)."),
                ("Threads:", self._pt_threads, "Threads", "CPU threads for Chorus2.")):
            self._pt_chorus_fields.append(c.row(label, w, title, text))
        self._pt_env = QLineEdit(ps.conda_env, self)
        self._pt_env.setPlaceholderText(pt.DEFAULT_ENV)
        self._pt_env_field = c.row("Conda env:", self._pt_env, "Chorus2 env",
                                            "The conda env (name or path) with Chorus2, bwa "
                                            "and jellyfish.<br><br>" + pt.install_help(rich=True),
                                            fixed=False)
        self._pt_status = QLabel(self)
        self._pt_status.setWordWrap(True)
        self._pt_status.setMinimumWidth(420)
        self._pt_status.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._pt_status_field = c.widget_row(self._pt_status)
        self._pt_env.editingFinished.connect(self._refresh_chorus_status)

        self._pt_ngs_fields = c.section("NGS Filter")
        self._pt_ngs = QCheckBox("Score probes with genomic reads (ChorusNGSfilter)", self)
        self._pt_ngs.setToolTip("With an imported BED: use its k-mer score column")
        self._pt_ngs.setChecked(ps.use_ngs)
        self._pt_ngs_fields.append(c.check(self._pt_ngs, "NGS filter", PAINT_NGS_HELP))
        self._pt_fastq = _path_row(self, "Reads FASTQ (.fq.gz); several: separate with commas",
                                   "Choose the FASTQ reads",
                                   file_filter="FASTQ (*.fq *.fastq *.fq.gz *.fastq.gz);;"
                                               "All files (*)")
        self._pt_fastq.set_path(ps.ngs_fastq or None, emit=False)
        self._pt_ngs_fields.append(c.row("Reads:", self._pt_fastq, fixed=False))
        self._pt_kmer = spin(13, 31, ps.ngs_kmer)
        self._pt_ngs_fields.append(c.row("k-mer:", self._pt_kmer, "k-mer",
                                         "k-mer length for counting the reads (default 17)."))
        self._pt_q = RangeEdit(0, 100, opt=False, suffix=" %")
        self._pt_q.set_values(round(ps.ngs_low_quantile * 100), round(ps.ngs_high_quantile * 100))
        self._pt_ngs_fields.append(c.row("Keep score quantiles:", self._pt_q, "Score window",
                                         "Probes with a k-mer score in the lowest or highest "
                                         "quantiles are dropped: too low = not in the reads, "
                                         "too high = repetitive (default 10 % to 90 %).",
                                         fixed=False))
        srv = QPushButton("Copy Server Command", self)
        srv.setToolTip("Copy the Chorus2 command for a Linux server to the clipboard")
        srv.clicked.connect(self._copy_server_command)
        srv.setSizePolicy(QSizePolicy.Policy.Fixed, QSizePolicy.Policy.Fixed)
        self._pt_server_field = c.check(srv, "Run on a server", PAINT_SERVER_HELP)
        lay.addWidget(body, 0, Qt.AlignmentFlag.AlignHCenter)
        self._paint_form = c
        self._add_page("Chorus2", page)

        # ---------------- Selection tab
        page, lay = _page(self)
        body = QWidget(page)
        g = _Form(body)
        self._forms.append(g)
        g.section("Pool")
        self._pt_spacing = spin(0, 100_000, ps.spacing, " bp")
        g.row("Min gap between probes:", self._pt_spacing, "Spacing",
              "Neighbouring probes are at least this far apart (default 25 bp).")
        self._pt_strand = QComboBox(self)
        self._pt_strand.addItem("Both (alternate + / −)", pt.STRAND_BOTH)
        self._pt_strand.addItem("Plus strand only", pt.STRAND_PLUS)
        _set_combo(self._pt_strand, ps.strand)
        g.row("Strand:", self._pt_strand, "Strand",
              "<b>Both</b> (Chorus2 default): every second probe is reverse-complemented, so "
              "neighbouring probes do not anneal to each other.<br><b>Plus</b>: all from one "
              "strand.")
        self._pt_target = QComboBox(self)
        self._pt_target.addItem("All probes", pt.TARGET_ALL)
        self._pt_target.addItem("Density (probes per kb)", pt.TARGET_DENSITY)
        self._pt_target.addItem("Pool size (number of probes)", pt.TARGET_POOL)
        _set_combo(self._pt_target, ps.target)
        g.row("Amount:", self._pt_target, "Amount", PAINT_AMOUNT_HELP)
        self._pt_density = dspin(0.01, 100, ps.density_per_kb, " per kb", step=0.1)
        self._pt_density_field = g.row("Density:", self._pt_density)
        self._pt_pool = spin(1, 10_000_000, ps.pool_size)
        self._pt_pool_field = g.row("Pool size:", self._pt_pool)
        self._pt_regions = QLineEdit(ps.regions, self)
        self._pt_regions.setPlaceholderText("All; or e.g. chr1, chr2:1-5000000")
        self._pt_regions.setClearButtonEnabled(True)
        g.row("Regions:", self._pt_regions, "Regions", PAINT_REGIONS_HELP, fixed=False)
        g.section("Order")
        self._pt_pool_name = QLineEdit(ps.pool_name, self)
        g.row("Pool name:", self._pt_pool_name, "Pool name",
              "Used in every oligo name and in the order file.", fixed=False)
        self._pt_p5 = QLineEdit(ps.primer5, self)
        self._pt_p5.setPlaceholderText("None")
        g.row("5′ adapter:", self._pt_p5, "Adapters", PAINT_ADAPTER_HELP, fixed=False)
        self._pt_p3 = QLineEdit(ps.primer3, self)
        self._pt_p3.setPlaceholderText("None")
        g.row("3′ adapter:", self._pt_p3, fixed=False)
        lay.addWidget(body, 0, Qt.AlignmentFlag.AlignHCenter)
        self._select_form = g
        self._add_page("Selection", page)

        for w in (self._pt_run, self._pt_ngs):
            w.toggled.connect(self._sync_painting)
        self._pt_target.currentIndexChanged.connect(self._sync_painting)
        self._sync_painting()
        self._refresh_chorus_status()

    def _sync_painting(self, *_a) -> None:
        run = self._pt_run.isChecked()
        c = self._paint_form
        c.set_visible(self._pt_input_field, run)
        c.set_visible(self._pt_bed_field, not run)
        # Import: the Chorus2 values still feed the server command and the methods text;
        # only the local env is irrelevant.
        c.set_visible(self._pt_env_field, run)
        c.set_visible(self._pt_status_field, run)
        ngs = self._pt_ngs.isChecked()
        for w in (self._pt_fastq, self._pt_kmer, self._pt_q):
            w.setEnabled(ngs)
        target = self._pt_target.currentData()
        self._select_form.set_visible(self._pt_density_field, target == pt.TARGET_DENSITY)
        self._select_form.set_visible(self._pt_pool_field, target == pt.TARGET_POOL)
        label = c.form.labelForField(self._pt_genome_field)
        if label is not None:
            label.setText("Genome:" if run else "Genome (optional):")
        if self.isVisible():
            QTimer.singleShot(0, self._fit)

    def _refresh_chorus_status(self) -> None:
        inst = pt.find_chorus2(self._pt_env.text().strip() or pt.DEFAULT_ENV)
        if inst.ok:
            self._pt_status.setText(f"<span style='color:palette(mid)'>✓ {inst.message}</span>")
        else:
            self._pt_status.setText(
                f"<span style='color:#b45309'>{inst.message}</span> "
                "<span style='color:palette(mid)'>Install it (see ?), or run it on a server "
                "and import the BED.</span>")

    def _copy_server_command(self) -> None:
        from PySide6.QtWidgets import QApplication

        ps = self.painting_settings()
        if not ps.genome_fasta:
            QMessageBox.information(self, "Chromosome painting",
                                    "Choose the genome FASTA first.")
            return
        try:
            text = pt.server_command_text(ps)
        except PrimerError as exc:
            QMessageBox.information(self, "Chromosome painting", str(exc))
            return
        QApplication.clipboard().setText(text)
        QMessageBox.information(self, "Chromosome painting",
                                "Copied to the clipboard:\n\n" + text)

    def _load_painting(self, ps: pt.PaintingSettings) -> None:
        (self._pt_import if ps.source == pt.SOURCE_IMPORT else self._pt_run).setChecked(True)
        self._pt_genome.set_path(ps.genome_fasta or None, emit=False)
        self._pt_input.set_path(ps.input_fasta or None, emit=False)
        self._pt_bed.set_path(ps.import_bed or None, emit=False)
        self._pt_len.setValue(ps.length)
        self._pt_hom.setValue(ps.homology)
        self._pt_dtm.setValue(ps.dtm)
        self._pt_step.setValue(ps.step)
        self._pt_ploidy.setValue(ps.ploidy)
        self._pt_threads.setValue(ps.threads)
        self._pt_env.setText(ps.conda_env)
        self._pt_ngs.setChecked(ps.use_ngs)
        self._pt_fastq.set_path(ps.ngs_fastq or None, emit=False)
        self._pt_kmer.setValue(ps.ngs_kmer)
        self._pt_q.set_values(round(ps.ngs_low_quantile * 100), round(ps.ngs_high_quantile * 100))
        self._pt_spacing.setValue(ps.spacing)
        _set_combo(self._pt_strand, ps.strand)
        _set_combo(self._pt_target, ps.target)
        self._pt_density.setValue(ps.density_per_kb)
        self._pt_pool.setValue(ps.pool_size)
        self._pt_regions.setText(ps.regions)
        self._pt_pool_name.setText(ps.pool_name)
        self._pt_p5.setText(ps.primer5)
        self._pt_p3.setText(ps.primer3)
        self._sync_painting()

    def painting_settings(self) -> pt.PaintingSettings:
        def path(row: PathBrowseRow) -> str:
            p = row.path()
            return str(p) if p else ""

        lo, hi = self._pt_q.values()
        fastq = self._pt_fastq.line_edit().text().strip()
        return pt.PaintingSettings(
            source=pt.SOURCE_IMPORT if self._pt_import.isChecked() else pt.SOURCE_RUN,
            genome_fasta=path(self._pt_genome), input_fasta=path(self._pt_input),
            import_bed=path(self._pt_bed),
            conda_env=self._pt_env.text().strip() or pt.DEFAULT_ENV,
            threads=self._pt_threads.value(), length=self._pt_len.value(),
            homology=self._pt_hom.value(), dtm=self._pt_dtm.value(), step=self._pt_step.value(),
            ploidy=self._pt_ploidy.value(), use_ngs=self._pt_ngs.isChecked(),
            ngs_fastq=fastq, ngs_kmer=self._pt_kmer.value(),
            ngs_low_quantile=lo / 100.0, ngs_high_quantile=hi / 100.0,
            spacing=self._pt_spacing.value(), strand=str(self._pt_strand.currentData()),
            target=str(self._pt_target.currentData()),
            density_per_kb=self._pt_density.value(), pool_size=self._pt_pool.value(),
            regions=self._pt_regions.text().strip(), primer5=self._pt_p5.text().strip(),
            primer3=self._pt_p3.text().strip(),
            pool_name=self._pt_pool_name.text().strip() or "Pool1",
        )

    def _toggle_advanced(self, on: bool) -> None:
        self._adv_btn.setText("Advanced Settings ▾" if on else "Advanced Settings ▸")
        for w in self._adv_fields:
            self._primers_form.set_visible(w, on)
        if self.isVisible():
            QTimer.singleShot(0, self._fit)

    def _on_mode(self, _i: int) -> None:
        mode = str(self._mode.currentData())
        if self._paint_mode.currentIndex() != self._mode.currentIndex():
            self._paint_mode.blockSignals(True)
            self._paint_mode.setCurrentIndex(self._mode.currentIndex())
            self._paint_mode.blockSignals(False)
        if mode != dz.MODE_PAINTING:
            self._apply_preset(mode)
        self._apply_mode_rows()
        self._apply_tab_family(select=True)

    def _apply_mode_rows(self) -> None:
        mode = str(self._mode.currentData())
        shown = set(id(w) for w in self._mode_rows.get(mode, []))
        every = {id(w): w for rows in self._mode_rows.values() for w in rows}
        for key, w in every.items():
            self._primers_form.set_visible(w, key in shown)
        # FISH designs single probes: no product, no pair Tm difference.
        probes = mode == dz.MODE_FISH
        self._primers_form.set_visible(self._product_field, not probes)
        self._primers_form.set_visible(self._tm_diff_field, not probes)
        self._conditions_title.setText(
            "<b>Probe Conditions</b>" if probes else "<b>Primer Conditions</b>")
        form = self._primers_form.form
        for field, text in ((self._size_field, "Probe size:" if probes else "Primer size:"),
                            (self._num_field, "Probes to return:" if probes else "Pairs to return:")):
            label = form.labelForField(field)
            if label is not None:
                label.setText(text)
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
        if mode == dz.MODE_PAINTING:
            keep = self.painting_settings()
            self._load_painting(pt.PaintingSettings(
                source=keep.source, genome_fasta=keep.genome_fasta,
                input_fasta=keep.input_fasta, import_bed=keep.import_bed,
                ngs_fastq=keep.ngs_fastq, conda_env=keep.conda_env))
            return
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
        self._min_cov.setValue(int(round(s.min_coverage * 100)))
        self._cov_mm.setValue(s.coverage_mismatches)
        self._apply_mode_rows()

    # ------------------------------------------------------------- inputs

    def _refresh_template_label(self) -> None:
        if not self._templates:
            self._template_label.setText(
                "<span style='color:palette(mid)'>No FASTA selected in the Files panel</span>")
            self._template_label.setToolTip("")
            return
        names = ", ".join(p.name for p in self._templates[:3])
        more = f" +{len(self._templates) - 3}" if len(self._templates) > 3 else ""
        self._template_label.setText(f"<b>{names}</b>{more}")
        self._template_label.setToolTip("\n".join(str(p) for p in self._templates))

    def uses_selected(self) -> bool:
        """Use selected sequences (the Files-panel selection) is on."""
        return self._input_selected.isChecked()

    def template_files(self) -> list[Path]:
        """The FASTA files to design on (Define inputs, else the Files-panel selection)."""
        if self._input_define.isChecked():
            p = self._template_row.path()
            return [Path(p)] if p else []
        return list(self._templates)

    def set_templates(self, files: list[Path]) -> None:
        """Files-panel selection (used with Use selected sequences)."""
        self._templates = [Path(p) for p in files if _is_fasta(Path(p))]
        self._refresh_template_label()

    def apply_path_from_panel(self, path: Path) -> bool:
        """Files-panel click: fills the path field being edited, else the template
        (the Template field with Define inputs, the selection otherwise)."""
        resolved = Path(path).expanduser()
        if not resolved.is_file():
            return False
        row = self._focused_path_row(self._all_path_rows())
        if row is not None:
            row.set_path(resolved)
            return True
        if self._painting():
            # Painting: a BED fills the import field, a FASTA the genome.
            if self._pt_import.isChecked() and resolved.suffix.lower() == ".bed":
                self._pt_bed.set_path(resolved)
                return True
            if _is_fasta(resolved):
                self._pt_genome.set_path(resolved)
                return True
            return False
        if not _is_fasta(resolved):
            return False
        if self._input_define.isChecked():
            self._template_row.set_path(resolved)
        else:
            self.set_templates([resolved])
        return True

    def current_dialog_tab(self) -> str:
        i = self._tab_bar.currentIndex()
        return self._TAB_NAMES[i] if 0 <= i < len(self._TAB_NAMES) else DIALOG_TAB_OPTIONS

    def design_settings(self) -> dz.DesignSettings:
        smin, sopt, smax = self._size.values()
        tmin, topt, tmax = self._tm.values()
        gmin, gmax = self._gc.values()
        pmin, pmax = self._product.values()
        psmin, psopt, psmax = self._probe_size.values()
        ptmin, ptopt, ptmax = self._probe_tm.values()

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
            forward_tail=self._fwd_tail.tail(), reverse_tail=self._rev_tail.tail(),
            reference_fasta=path(self._reference),
            min_coverage=self._min_cov.value() / 100.0,
            coverage_mismatches=self._cov_mm.value(),
            fish_label=str(self._fish_label.currentData()),
            fish_avoid_fasta=path(self._fish_avoid),
        )

    def options(self) -> PrimerOptions:
        return replace(
            self._opts,
            design=self.design_settings().to_dict(),
            painting=self.painting_settings().to_dict(),
            right_click_action=PRIMER_RIGHT_CLICK_RUN,
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
        # Long saved names must not widen the Map Primers window.
        self.combo.setSizeAdjustPolicy(
            QComboBox.SizeAdjustPolicy.AdjustToMinimumContentsLengthWithIcon)
        self.combo.setMinimumContentsLength(14)
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
        self._input_group = QButtonGroup(self)
        self._input_selected = QRadioButton("Use selected sequences", self)
        self._input_define = QRadioButton("Define inputs", self)
        self._input_group.addButton(self._input_selected)
        self._input_group.addButton(self._input_define)
        (self._input_selected if self._files else self._input_define).setChecked(True)
        f.check(self._input_selected, "Inputs", MAP_INPUTS_HELP)
        self._seq_label = QLabel(self)
        self._seq_label.setWordWrap(True)
        self._seq_label.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        self._seq_label.setIndent(22)  # under the radio text, not its circle
        self._seq_label_field = f.widget_row(self._seq_label)
        f.check(self._input_define)
        self._seq_row = _path_row(self, "Sequences FASTA: click it in Files, or Browse",
                                  "Choose the FASTA to test the primers on")
        self._seq_field = f.row("Sequences:", self._seq_row, "Sequences",
                                "The FASTA the primers are tested on. Each sequence in it is "
                                "searched on both strands.", fixed=False)
        self._input_selected.toggled.connect(self._sync_inputs)
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
        self._use_f = QCheckBox("Forward:", self)
        self._use_r = QCheckBox("Reverse:", self)
        self._use_p = QCheckBox("Probe:", self)
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
            self._spec_rows.append(f.switch_row(box, field))
        f.check(self._database)
        self._src = QComboBox(self)
        self._src.addItem("Primer Database", "db")
        self._src.addItem("FASTA file", "file")
        self._src_field = f.row("Primers from:", self._src)
        self._src_file = _path_row(self, "Primers FASTA: click it in Files, or Browse",
                                   "Choose a primers FASTA")
        if ms.primer_file:
            self._src.setCurrentIndex(1)
            self._src_file.set_path(ms.primer_file, emit=False)
        self._src_file_field = f.row("Primers FASTA:", self._src_file)
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
        kl.addStretch(1)
        self._kinds_field = f.row("Search for:", kinds, "Search for",
                                  "Which saved primers to test. <b>Pairs only</b> keeps primers "
                                  "that form a product with another saved primer.")
        # From the Database window: only the primers selected there.
        self._only_names: list[str] = list(ms.primer_names or [])
        self._only_label = QLabel(self)
        self._only_label.setWordWrap(True)
        self._only_field = f.row("Testing:", self._only_label)
        self._only_all = QPushButton("Use All Saved Primers", self)
        self._only_all.clicked.connect(self._clear_only)
        self._only_all_field = f.check(self._only_all)
        self._specific.toggled.connect(self._sync)
        self._src.currentIndexChanged.connect(self._sync)

        f.section("Products")

        def region(check_text, on, a, b, suffix=""):
            box = QCheckBox(check_text, self)
            box.setChecked(on)
            rng = RangeEdit(1, 2_000_000_000, opt=False, suffix=suffix)
            rng.set_values(a, b)
            rng.setEnabled(on)
            box.toggled.connect(rng.setEnabled)
            f.switch_row(box, rng)
            return box, rng

        self._use_inc, self._inc = region("Included region:", ms.use_included,
                                          ms.included_start, ms.included_end or 1)
        self._use_tgt, self._tgt = region("Target region:", ms.use_target, ms.target_start,
                                          ms.target_end)
        self._use_prod, self._prod = region("Product size:", ms.use_product_range,
                                            ms.product_min, ms.product_max, " bp")
        self._use_opt = QCheckBox("Optimal size:", self)
        self._use_opt.setChecked(ms.use_optimal)
        self._opt = QSpinBox(self)
        self._opt.setRange(1, 1_000_000)
        self._opt.setValue(max(1, ms.optimal_size or 500))
        self._opt.setSuffix(" bp")
        self._opt.setEnabled(ms.use_optimal)
        self._use_opt.toggled.connect(self._opt.setEnabled)
        f.switch_row(self._use_opt, self._opt, "Optimal product size",
                     "Sort products by how close they are to this size.")

        f.section("Mismatches")
        self._mm = QSpinBox(self)
        self._mm.setRange(0, 10)
        self._mm.setValue(ms.max_mismatches)
        f.row("Max mismatches:", self._mm, "Mismatches",
              "Mismatches a primer may have at a binding site.")
        self._exact = QSpinBox(self)
        self._exact.setRange(0, 15)
        self._exact.setValue(ms.three_prime_exact)
        f.row("3′ bases exact:", self._exact, "3′ end",
              "The last bases at the 3′ end must match: mismatches there stop extension.")

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
        self._genome_fields.append(f.row("Genome:", self._genome))
        lay.addWidget(body, 0, Qt.AlignmentFlag.AlignHCenter)
        self._form_main = f
        self._add_page("Map Primers", page)

        # Options tab
        page, lay = _page(self)
        lay.setContentsMargins(16, 14, 16, 14)
        lay.setSpacing(14)
        rc_box = QGroupBox("Right-click on Map Primers", self)
        rcl = QVBoxLayout(rc_box)
        self._rc_run = QRadioButton("Map with saved settings", self)
        self._rc_run.setChecked(True)
        rgrp = QButtonGroup(self)
        rgrp.addButton(self._rc_run)
        rcl.addWidget(self._rc_run)
        lay.addWidget(rc_box)
        res_box = QGroupBox("Results", self)
        rl = QVBoxLayout(res_box)
        self._annotate = QCheckBox("Add binding sites and products to the sequence annotations",
                                   self)
        self._annotate.setChecked(ms.annotate)
        rl.addWidget(self._annotate)
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
        self._sync_inputs()

    def _sync_inputs(self, *_a) -> None:
        define = self._input_define.isChecked()
        self._form_main.set_visible(self._seq_label_field, not define)
        self._form_main.set_visible(self._seq_field, define)
        self._refresh_seq_label()
        if self.isVisible():
            QTimer.singleShot(0, self._fit)

    def _database_obj(self) -> PrimerDatabase | None:
        return self._db

    def _reload_db(self) -> None:
        p = self._opts.database_path or None
        try:
            from repeatviewer.backends.primers.inventory import open_database

            self._db = open_database(p)
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

    def _clear_only(self) -> None:
        self._only_names = []
        self._sync()

    def _sync(self, *_a) -> None:
        spec = self._specific.isChecked()
        from_file = self._src.currentData() == "file"
        only = bool(self._only_names) and not spec and not from_file
        form = self._form_main
        form.set_visible(self._only_field, only)
        form.set_visible(self._only_all_field, only)
        if only:
            names = ", ".join(self._only_names[:6]) + (" …" if len(self._only_names) > 6 else "")
            self._only_label.setText(f"<b>{len(self._only_names)} selected in the Primer "
                                     f"Database</b>: {names}")
        form.set_visible(self._kinds_field, not only)
        for b, fld in ((self._use_f, self._f), (self._use_r, self._r), (self._use_p, self._p)):
            b.setEnabled(spec)
            fld.setEnabled(spec and b.isChecked())
        for w in (self._src, self._s_f, self._s_r, self._s_p, self._pairs_only):
            w.setEnabled(not spec)
        form.set_visible(self._src_file_field, from_file)
        self._src_file.setEnabled(not spec)
        cons = self._test_on.currentData() == pj.TEST_ON_CONSENSUS
        form.set_visible(self._cons_field, cons)
        form.set_visible(self._cons_gaps_field, cons)
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
            self._seq_label.setText("<span style='color:palette(mid)'>No FASTA selected in the "
                                    "Files panel</span>")
            return
        names = ", ".join(p.name for p in self._files[:3])
        more = f" +{len(self._files) - 3}" if len(self._files) > 3 else ""
        self._seq_label.setText(f"<b>{names}</b>{more}")

    def set_files(self, files: list[Path]) -> None:
        self._files = [Path(p) for p in files if _is_fasta(Path(p))]
        self._refresh_seq_label()

    def uses_selected(self) -> bool:
        return self._input_selected.isChecked()

    def run_files(self) -> list[Path]:
        """FASTA files to test on (Define inputs, else the Files-panel selection)."""
        if self._input_define.isChecked():
            p = self._seq_row.path()
            return [Path(p)] if p else []
        return list(self._files)

    def apply_path_from_panel(self, path: Path) -> bool:
        resolved = Path(path).expanduser()
        if not resolved.is_file():
            return False
        row = self._focused_path_row([self._seq_row, self._src_file, self._genome])
        if row is not None:
            row.set_path(resolved)
            return True
        if not _is_fasta(resolved):
            return False
        if self._input_define.isChecked():
            self._seq_row.set_path(resolved)
        else:
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
            primer_names=list(self._only_names) if not self._specific.isChecked() else [],
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
        return replace(
            self._opts, map=self.map_settings().to_dict(),
            map_right_click_action=PRIMER_RIGHT_CLICK_RUN,
            show_process_window=self._show_process.isChecked(),
            close_process_window_when_finished=self._close_process.isChecked(),
        )


# --------------------------------------------------------------------------
# Primer Database
# --------------------------------------------------------------------------


def primer3_missing_message(parent: QWidget | None) -> bool:
    """Show the install hint when primer3-py is missing; True when it is."""
    if primer3_available():
        return False
    QMessageBox.warning(parent, "Primers", PRIMER3_MISSING)
    return True


# The Database window lives in its own module (it imports helpers from here).
from repeatviewer.views.primer_database_dialog import PrimerDatabaseDialog  # noqa: E402

__all__ = [
    "MapPrimersDialog",
    "PrimerDatabaseDialog",
    "PrimerDesignDialog",
    "PrimerPickerDialog",
    "primer3_missing_message",
]
