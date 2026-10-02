"""Map to Reference dialog (Tools → Map to Reference).

Mapping tab: the reference (click it in Files, or Browse), the mapper, the
sensitivity, what to do with sequences that fit in several places, and the
annotation check. Advanced: the limits a placement must pass. Output: consensus
and the files to write. Options: right-click action, process window, programs
from PATH. Cancel / Save / Run, as in every RepeatViewer tool dialog.

The sequences to map are the files selected in the file panel when the dialog
opens; a click in Files while it is open fills the Reference field.
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
    QPushButton,
    QRadioButton,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import map_reference as mr
from repeatviewer.job_options import (
    MAP_CONSENSUS_THRESHOLDS,
    MAP_ENGINE_AUTO,
    MAP_ENGINE_BOWTIE2,
    MAP_ENGINE_MAFFT,
    MAP_ENGINE_MINIMAP2,
    MAP_MINIMAP2_PRESETS,
    MAP_MULTI_ALL,
    MAP_MULTI_BEST,
    MAP_MULTI_NONE,
    MAP_NO_COVERAGE_N,
    MAP_NO_COVERAGE_REFERENCE,
    MAP_NO_COVERAGE_SKIP,
    MAP_RIGHT_CLICK_DIALOG,
    MAP_RIGHT_CLICK_RUN,
    MAP_SENS_FAST,
    MAP_SENS_HIGH,
    MAP_SENS_MEDIUM,
    MapOptions,
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

_FIELD_WIDTH = 320
_NUMBER_WIDTH = 110
_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".faa", ".fsa", ".ffn", ".seq"}
_SCAN_LIMIT_BYTES = 30_000_000     # count sequences for the summary line up to this size

REFERENCE_INFO = (
    "The sequence(s) to map to: a satellite monomer or array, a TE consensus, a BAC or a "
    "genome. Click it in the Files panel, or use Browse.<br>"
    "Annotations of the reference (its <b>.RV.gff3</b>, made in the sequence viewer or by "
    "Sequence Annotation) are shown above the mapped sequences and checked in each of them."
)
MAPPER_INFO = (
    "<b>Auto</b>: MAFFT for a few sequences on one short reference (up to 500 sequences, "
    "50 kb), Bowtie2 for many short reads (under 400 bp), otherwise minimap2.<br>"
    "<b>MAFFT add fragments</b> (bundled): the most accurate placement of a few sequences, "
    "e.g. satellite variants on an annotated monomer array. One reference sequence only.<br>"
    "<b>minimap2</b>: long reads (ONT, HiFi) or assembled contigs against a genome.<br>"
    "<b>Bowtie2</b>: short reads (Illumina). Bisulfite and EM-seq reads need a "
    "conversion-aware mapper such as Bismark and are not handled here."
)
SENSITIVITY_INFO = (
    "<b>Fast</b>: the mapper's defaults.<br>"
    "<b>Medium</b>: smaller seeds, finds more diverged copies (MAFFT: --multipair; "
    "minimap2: -k13 -w7; Bowtie2: --sensitive-local).<br>"
    "<b>High</b>: the most sensitive; for old, diverged repeat copies. Slower "
    "(minimap2: -k11 -w5; Bowtie2: --very-sensitive-local)."
)
MULTI_INFO = (
    "In a satellite array or a TE family a sequence often fits equally well in several "
    "places.<br>"
    "<b>Best place only</b>: one place, the mapper's best (ties chosen at random).<br>"
    "<b>Every equally good place</b>: the sequence is shown in all of them "
    "(minimap2: score within 95 % of the best; Bowtie2: up to 20 places).<br>"
    "<b>Do not place</b>: such sequences are listed with the unused ones as ambiguous.<br>"
    "MAFFT maps to a single sequence, so this does not apply."
)
LIFT_INFO = (
    "Every feature of the reference (subrepeats, motifs, domains…) is looked up in every "
    "mapped sequence: <b>intact</b>, <b>changed</b> (with the number of differing bases), "
    "<b>partial</b> (only part covered), or <b>lost</b> (half or more deleted).<br>"
    "The Features tab of the result shows this as a table, and the placed sequences and the "
    "consensus are written with the annotations moved onto them."
)
LIMITS_INFO = (
    "A placement is kept only when it passes all limits; the others are listed with the "
    "unused sequences, with the reason.<br>"
    "<b>Minimum identity</b>: share of aligned columns that match.<br>"
    "<b>Minimum aligned length</b>: bases of the sequence placed on the reference.<br>"
    "<b>Minimum mapping quality</b> (minimap2, Bowtie2): 0 keeps all."
)
CONSENSUS_INFO = (
    "<b>Majority</b>: the most common base at each position (ties as IUPAC codes).<br>"
    "<b>50 % … 100 %</b>: the fewest bases that together reach this share of the sequences "
    "at the position, written as one IUPAC code (e.g. R = A or G).<br>"
    "A deletion wins only when it is the most common state; an insertion is included when "
    "at least half the sequences (or the threshold) have it."
)
PRESET_INFO = (
    "<b>Auto</b>: sr for short reads, asm20 for long assembled contigs, map-ont otherwise "
    "(works for ONT and HiFi). Choose map-hifi for HiFi reads, asm5 / asm10 / asm20 for "
    "assemblies at < 1 %, < 5 %, < 20 % divergence."
)
DISPLAY_INFO = (
    "How many placements the contig view draws per reference sequence. Coverage, consensus "
    "and variants always use all of them; above this number a random sample is drawn."
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


def _form(parent: QWidget) -> QFormLayout:
    f = QFormLayout(parent)
    f.setContentsMargins(16, 14, 16, 14)
    f.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    f.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
    f.setHorizontalSpacing(12)
    f.setVerticalSpacing(9)
    return f


def _is_fasta(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return Path(name).suffix in _FASTA_SUFFIXES


def input_summary_text(files: list[Path], reference: Path | None) -> str:
    """The grey line under the reference: what will be mapped to what."""
    if not files:
        return "Select the FASTA / FASTQ files to map in the Files panel first."
    names = ", ".join(p.name for p in files[:3]) + (" …" if len(files) > 3 else "")
    try:
        size = sum(p.stat().st_size for p in files)
    except OSError:
        size = _SCAN_LIMIT_BYTES + 1
    if size <= _SCAN_LIMIT_BYTES:
        try:
            n = mr.scan_inputs(files, keep_names=0).n
            what = f"{n:,} sequence{'s' if n != 1 else ''} ({names})"
        except OSError:
            what = names
    else:
        what = f"{len(files)} file{'s' if len(files) != 1 else ''} ({names})"
    target = reference.name if reference else "the reference"
    return f"{what} will be mapped to {target}."


class MapReferenceDialog(QDialog):
    """Map the selected sequences or reads to a reference."""

    def __init__(
        self,
        *,
        initial: MapOptions | None = None,
        files: list[Path] | None = None,
        title: str = "Map to Reference",
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        opts = initial or MapOptions()
        self._files = [Path(p) for p in (files or [])]
        self._ref_info_cache: dict[str, tuple] = {}

        # --- Mapping tab -------------------------------------------------------
        self._reference = PathBrowseRow(
            self,
            initial=opts.reference_fasta or None,
            placeholder="Reference FASTA: click it in Files, or Browse",
            browse_title="Choose the reference FASTA",
            file_filter="FASTA (*.fa *.fasta *.fna *.fas *.fa.gz *.fasta.gz);;All files (*)",
        )
        self._reference.setFixedWidth(_FIELD_WIDTH)
        self._summary = _note("", self)
        self._summary.setFixedWidth(_FIELD_WIDTH + 30)
        self._engine = QComboBox(self)
        for value in (MAP_ENGINE_AUTO, MAP_ENGINE_MAFFT, MAP_ENGINE_MINIMAP2, MAP_ENGINE_BOWTIE2):
            self._engine.addItem(mr.ENGINE_LABELS[value], value)
        self._set_combo(self._engine, opts.engine)
        self._engine_note = _note("", self)
        self._engine_note.setFixedWidth(_FIELD_WIDTH + 30)
        self._sensitivity = QComboBox(self)
        for label, value in (("Fast", MAP_SENS_FAST), ("Medium", MAP_SENS_MEDIUM),
                             ("High, for diverged repeat copies", MAP_SENS_HIGH)):
            self._sensitivity.addItem(label, value)
        self._set_combo(self._sensitivity, opts.sensitivity)
        self._multi = QComboBox(self)
        for label, value in (("Best place only", MAP_MULTI_BEST),
                             ("Every equally good place", MAP_MULTI_ALL),
                             ("Do not place (list as ambiguous)", MAP_MULTI_NONE)):
            self._multi.addItem(label, value)
        self._set_combo(self._multi, opts.multiple_hits)
        self._lift = QCheckBox("Check the reference annotations in every sequence", self)
        self._lift.setChecked(opts.lift_annotations)
        for w in (self._engine, self._sensitivity, self._multi):
            w.setFixedWidth(_FIELD_WIDTH)

        mapping_tab = QWidget(self)
        f = _form(mapping_tab)
        f.addRow("Reference:", _row(self, self._reference,
                                    make_help_button(self, "Reference", REFERENCE_INFO)))
        f.addRow("", self._summary)
        f.addRow("Mapper:", _row(self, self._engine,
                                 make_help_button(self, "Mapper", MAPPER_INFO)))
        f.addRow("", self._engine_note)
        f.addRow("Sensitivity:", _row(self, self._sensitivity,
                                      make_help_button(self, "Sensitivity", SENSITIVITY_INFO)))
        f.addRow("Several places:", _row(self, self._multi,
                                         make_help_button(self, "Several equally good places",
                                                          MULTI_INFO)))
        f.addRow("Annotations:", _row(self, self._lift,
                                      make_help_button(self, "Annotation check", LIFT_INFO)))

        # --- Advanced tab --------------------------------------------------------
        self._min_identity = QSpinBox(self)
        self._min_identity.setRange(0, 100)
        self._min_identity.setSuffix(" %")
        self._min_identity.setValue(round(opts.min_identity * 100))
        self._min_aligned = QSpinBox(self)
        self._min_aligned.setRange(1, 1_000_000)
        self._min_aligned.setSuffix(" bp")
        self._min_aligned.setValue(opts.min_aligned_bp)
        self._min_mapq = QSpinBox(self)
        self._min_mapq.setRange(0, 255)
        self._min_mapq.setValue(opts.min_mapq)
        self._preset = QComboBox(self)
        for value in MAP_MINIMAP2_PRESETS:
            self._preset.addItem("Auto" if value == "auto" else value, value)
        self._set_combo(self._preset, opts.minimap2_preset)
        self._preset.setFixedWidth(160)
        self._threads = QSpinBox(self)
        self._threads.setRange(1, 256)
        self._threads.setValue(max(1, opts.threads))
        for w in (self._min_identity, self._min_aligned, self._min_mapq, self._threads):
            w.setFixedWidth(_NUMBER_WIDTH)
        advanced_tab = QWidget(self)
        a = _form(advanced_tab)
        a.addRow("Minimum identity:", _row(self, self._min_identity,
                                           make_help_button(self, "Limits", LIMITS_INFO)))
        a.addRow("Minimum aligned length:", self._min_aligned)
        a.addRow("Minimum mapping quality:", self._min_mapq)
        a.addRow("minimap2 preset:", _row(self, self._preset,
                                          make_help_button(self, "minimap2 preset", PRESET_INFO)))
        a.addRow("Threads:", self._threads)
        self._advanced_form = a

        # --- Output tab -------------------------------------------------------------
        self._threshold = QComboBox(self)
        for t in MAP_CONSENSUS_THRESHOLDS:
            label = "Majority (most common base)" if t == 0 else (
                "100 % (identical)" if t == 1 else f"{t:.0%}")
            self._threshold.addItem(label, t)
        self._set_combo(self._threshold, opts.consensus_threshold)
        self._no_cov = QComboBox(self)
        for label, value in (("N", MAP_NO_COVERAGE_N), ("Reference base", MAP_NO_COVERAGE_REFERENCE),
                             ("Leave out", MAP_NO_COVERAGE_SKIP)):
            self._no_cov.addItem(label, value)
        self._set_combo(self._no_cov, opts.no_coverage)
        self._min_var = QSpinBox(self)
        self._min_var.setRange(0, 100)
        self._min_var.setSuffix(" %")
        self._min_var.setSpecialValueText("Every difference")
        self._min_var.setValue(round(opts.min_variant_freq * 100))
        self._min_var.setFixedWidth(150)
        self._display = QSpinBox(self)
        self._display.setRange(50, 100_000)
        self._display.setSingleStep(500)
        self._display.setValue(opts.max_display_reads)
        self._display.setFixedWidth(_NUMBER_WIDTH)
        self._save_consensus = QCheckBox("Consensus", self)
        self._save_consensus.setChecked(opts.save_consensus)
        self._save_mapped = QCheckBox("Placed sequences (turned to the reference strand)", self)
        self._save_mapped.setChecked(opts.save_mapped)
        self._save_unused = QCheckBox("Unused sequences, with the reason", self)
        self._save_unused.setChecked(opts.save_unused)
        self._save_variants = QCheckBox("Variants table", self)
        self._save_variants.setChecked(opts.save_variants)
        for w in (self._threshold, self._no_cov):
            w.setFixedWidth(220)
        saves = QWidget(self)
        sl = QVBoxLayout(saves)
        sl.setContentsMargins(0, 0, 0, 0)
        sl.setSpacing(4)
        for w in (self._save_consensus, self._save_mapped, self._save_unused, self._save_variants):
            sl.addWidget(w)
        output_tab = QWidget(self)
        o = _form(output_tab)
        o.addRow("Consensus:", _row(self, self._threshold,
                                    make_help_button(self, "Consensus", CONSENSUS_INFO)))
        o.addRow("No coverage:", self._no_cov)
        o.addRow("Variants from:", self._min_var)
        o.addRow("Sequences drawn:", _row(self, self._display,
                                          make_help_button(self, "Sequences drawn", DISPLAY_INFO)))
        o.addRow("Save:", saves)
        o.addRow("", _note("Always written: the result (.rvmap, opens in the contig view), the "
                           "report, the reads and features tables and the aligner output, in "
                           "the <input>_map_details folder.", self))

        # --- Options tab ---------------------------------------------------------------
        self._right_click_buttons: dict[str, QRadioButton] = {}
        group = QButtonGroup(self)
        rc_box = QGroupBox("Right-click on Map to Reference", self)
        rcl = QVBoxLayout(rc_box)
        for label, value in (("Map with saved settings", MAP_RIGHT_CLICK_RUN),
                             ("Open this window", MAP_RIGHT_CLICK_DIALOG)):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            group.addButton(radio)
            self._right_click_buttons[value] = radio
            rcl.addWidget(radio)
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[MAP_RIGHT_CLICK_RUN].setChecked(True)
        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self, tool_name="the mapping", checked=opts.show_process_window,
                close_when_finished=bool(
                    getattr(opts, "close_process_window_when_finished", False)
                ),
            )
        )
        self._use_system_path = QCheckBox("Use programs from the system PATH", self)
        self._use_system_path.setChecked(opts.use_system_path)
        self._use_system_path.setToolTip(
            "Off: RepeatViewer's bundled MAFFT and minimap2; Bowtie2 from the conda env "
            "in File → Settings → External tools (then PATH). On: all from PATH.")
        self._found = _note("", self)
        self._found.setWordWrap(False)
        exe_box = QGroupBox("Programs", self)
        exl = QVBoxLayout(exe_box)
        exl.addWidget(self._use_system_path)
        exl.addWidget(self._found)
        from repeatviewer.backends.minimap2_tools import MINIMAP2_VERSION
        exl.addWidget(_note(
            f"minimap2 {MINIMAP2_VERSION} is included with RepeatViewer "
            "(Windows / Linux x86_64 / macOS). Bowtie2 still needs conda or PATH.",
            self))
        options_tab = QWidget(self)
        ol = QVBoxLayout(options_tab)
        ol.setContentsMargins(16, 14, 16, 14)
        ol.setSpacing(14)
        ol.addWidget(rc_box)
        ol.addWidget(process_box)
        ol.addWidget(exe_box)
        ol.addStretch(1)

        self._tabs = QTabWidget(self)
        self._tabs.addTab(mapping_tab, "Mapping")
        self._tabs.addTab(advanced_tab, "Advanced")
        self._tabs.addTab(output_tab, "Output")
        self._tabs.addTab(options_tab, "Options")

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._button_box = make_cancel_save_run_bar(self, left_widgets=[restore])
        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(12)
        layout.addWidget(self._tabs)
        layout.addWidget(self._button_box)

        self._reference.path_changed.connect(self._sync)
        self._engine.currentIndexChanged.connect(self._sync)
        self._use_system_path.toggled.connect(self._sync)
        self._sync()

    # ------------------------------------------------------------------ helpers

    @staticmethod
    def _set_combo(combo, value) -> None:
        idx = combo.findData(value)
        if idx >= 0:
            combo.setCurrentIndex(idx)

    def _reference_info(self, path: Path | None) -> tuple[dict[str, int], str]:
        """(record lengths, error) for the reference; cached per path."""
        if path is None:
            return {}, ""
        key = str(path)
        if key not in self._ref_info_cache:
            try:
                if not path.is_file():
                    raise OSError("not found")
                if path.stat().st_size > _SCAN_LIMIT_BYTES:
                    self._ref_info_cache[key] = ({}, "")
                else:
                    self._ref_info_cache[key] = (mr.reference_lengths(path), "")
            except OSError:
                self._ref_info_cache[key] = ({}, f"{path.name}: not found.")
        return self._ref_info_cache[key]

    def _auto_choice(self) -> str | None:
        lengths, _err = self._reference_info(self._reference.path())
        if not lengths or not self._files:
            return None
        try:
            size = sum(p.stat().st_size for p in self._files)
        except OSError:
            return None
        if size > _SCAN_LIMIT_BYTES:
            return None
        return mr.choose_engine(lengths, mr.scan_inputs(self._files, keep_names=0))

    def _sync(self, *_args) -> None:
        ref = self._reference.path()
        _lengths, err = self._reference_info(ref)
        self._summary.setText(err or input_summary_text(self._files, ref))
        self._summary.setStyleSheet("color: #b00020;" if err else "color: palette(mid);")
        system = self._use_system_path.isChecked()
        for i in range(self._engine.count()):
            value = self._engine.itemData(i)
            label = mr.ENGINE_LABELS[value]
            if value != MAP_ENGINE_AUTO and not mr.engine_available(value, system_path_only=system):
                label += " (not found)"
            self._engine.setItemText(i, label)
        engine = self._engine.currentData()
        note = ""
        if engine == MAP_ENGINE_AUTO:
            choice = self._auto_choice()
            if choice:
                note = f"Auto will use {mr.ENGINE_NAMES[choice]} for this input."
                if not mr.engine_available(choice, system_path_only=system):
                    note += " It was not found: " + mr.engine_install_hint(choice).split(":\n")[-1]
        elif engine == MAP_ENGINE_MAFFT and len(_lengths) > 1:
            note = (f"MAFFT maps to one sequence; {ref.name if ref else 'the reference'} has "
                    f"{len(_lengths)}. Choose minimap2 or Bowtie2.")
        elif engine != MAP_ENGINE_AUTO and not mr.engine_available(engine, system_path_only=system):
            note = "Not found. " + mr.engine_install_hint(engine).replace("\n", " ")
        self._engine_note.setText(note)
        self._engine_note.setVisible(bool(note))
        mafft = engine == MAP_ENGINE_MAFFT
        self._multi.setEnabled(not mafft)
        self._min_mapq.setEnabled(not mafft)
        self._preset.setEnabled(engine in (MAP_ENGINE_AUTO, MAP_ENGINE_MINIMAP2))
        self._update_found()

    def _update_found(self) -> None:
        from repeatviewer.backends.minimap2_tools import bundled_minimap2_executable
        from repeatviewer.backends.sequence_aligner import mafft_command

        system = self._use_system_path.isChecked()
        lines, paths = [], []
        cmd = mafft_command(use_system_mafft=system)
        if cmd and not system:
            mafft_label = "bundled"
        else:
            mafft_label = "found" if cmd else "not found"
        lines.append(f"MAFFT: {mafft_label}")
        paths.append(f"MAFFT: {cmd[0][0] if cmd else 'not found'}")
        if not system and bundled_minimap2_executable():
            exe = bundled_minimap2_executable()
            lines.append(f"{mr.ENGINE_NAMES[MAP_ENGINE_MINIMAP2]}: bundled")
            paths.append(f"{mr.ENGINE_NAMES[MAP_ENGINE_MINIMAP2]}: {exe}")
        else:
            exe = mr.engine_executable(MAP_ENGINE_MINIMAP2, system_path_only=system)
            lines.append(f"{mr.ENGINE_NAMES[MAP_ENGINE_MINIMAP2]}: "
                         f"{'found' if exe else 'not found'}")
            paths.append(f"{mr.ENGINE_NAMES[MAP_ENGINE_MINIMAP2]}: {exe or 'not found'}")
        if system:
            exe = mr.engine_executable(MAP_ENGINE_BOWTIE2, system_path_only=True)
            lines.append(f"{mr.ENGINE_NAMES[MAP_ENGINE_BOWTIE2]}: "
                         f"{'found' if exe else 'not found'}")
            paths.append(f"{mr.ENGINE_NAMES[MAP_ENGINE_BOWTIE2]}: {exe or 'not found'}")
        else:
            # Soft status without conda probe (runtime still searches conda then PATH).
            lines.append(f"{mr.ENGINE_NAMES[MAP_ENGINE_BOWTIE2]}: conda / PATH")
            paths.append(f"{mr.ENGINE_NAMES[MAP_ENGINE_BOWTIE2]}: "
                         + mr.engine_install_hint(MAP_ENGINE_BOWTIE2).replace("\n", " "))
        self._found.setText("\n".join(lines))
        self._found.setToolTip("\n".join(paths))

    def _restore_defaults(self) -> None:
        d = MapOptions()
        self._set_combo(self._engine, d.engine)
        self._set_combo(self._sensitivity, d.sensitivity)
        self._set_combo(self._multi, d.multiple_hits)
        self._lift.setChecked(d.lift_annotations)
        self._min_identity.setValue(round(d.min_identity * 100))
        self._min_aligned.setValue(d.min_aligned_bp)
        self._min_mapq.setValue(d.min_mapq)
        self._set_combo(self._preset, d.minimap2_preset)
        self._threads.setValue(max(1, d.threads))
        self._set_combo(self._threshold, d.consensus_threshold)
        self._set_combo(self._no_cov, d.no_coverage)
        self._min_var.setValue(round(d.min_variant_freq * 100))
        self._display.setValue(d.max_display_reads)
        for box, value in ((self._save_consensus, d.save_consensus),
                           (self._save_mapped, d.save_mapped),
                           (self._save_unused, d.save_unused),
                           (self._save_variants, d.save_variants)):
            box.setChecked(value)
        self._right_click_buttons[d.right_click_action].setChecked(True)
        self._show_process_window.setChecked(d.show_process_window)
        self._close_process_window.setChecked(d.close_process_window_when_finished)
        self._use_system_path.setChecked(d.use_system_path)
        self._sync()

    # ------------------------------------------------------------------ results

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def files(self) -> list[Path]:
        return list(self._files)

    def apply_path_from_panel(self, path: Path) -> bool:
        """A FASTA clicked in the Files panel becomes the reference."""
        resolved = Path(path).expanduser()
        if not resolved.is_file() or not _is_fasta(resolved):
            return False
        self._reference.set_path(resolved)
        return True

    def options(self) -> MapOptions:
        right_click = next((v for v, b in self._right_click_buttons.items() if b.isChecked()),
                           MAP_RIGHT_CLICK_RUN)
        ref = self._reference.path()
        return MapOptions(
            reference_fasta=str(ref) if ref else "",
            engine=str(self._engine.currentData()),
            sensitivity=str(self._sensitivity.currentData()),
            multiple_hits=str(self._multi.currentData()),
            minimap2_preset=str(self._preset.currentData()),
            min_identity=self._min_identity.value() / 100.0,
            min_aligned_bp=self._min_aligned.value(),
            min_mapq=self._min_mapq.value(),
            threads=self._threads.value(),
            consensus_threshold=float(self._threshold.currentData()),
            no_coverage=str(self._no_cov.currentData()),
            min_variant_freq=self._min_var.value() / 100.0,
            save_consensus=self._save_consensus.isChecked(),
            save_mapped=self._save_mapped.isChecked(),
            save_unused=self._save_unused.isChecked(),
            save_variants=self._save_variants.isChecked(),
            lift_annotations=self._lift.isChecked(),
            max_display_reads=self._display.value(),
            right_click_action=right_click,
            show_process_window=self._show_process_window.isChecked(),
            close_process_window_when_finished=(
                self._close_process_window.isChecked()
            ),
            use_system_path=self._use_system_path.isChecked(),
        )
