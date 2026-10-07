"""Mapping dialog (Tools → Mapping).

Tabs: **Minimap2** | **HISAT2** | **Options**. Each engine tab shows that
mapper's real CLI settings; RepeatViewer-only output options live under
Advanced. Queries come from the Files panel or Define inputs; a Files click
fills Reference.
"""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, QTimer
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
    QSizePolicy,
    QStackedWidget,
    QTabBar,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends import map_reference as mr
from repeatviewer.job_options import (
    ALIGN_INPUT_DEFINE,
    ALIGN_INPUT_SELECTED,
    DIALOG_TAB_OPTIONS,
    MAP_CONSENSUS_THRESHOLDS,
    MAP_DIALOG_TABS,
    MAP_ENGINE_AUTO,
    MAP_ENGINE_HISAT2,
    MAP_ENGINE_MINIMAP2,
    MAP_HISAT2_DNA,
    MAP_HISAT2_RNA,
    MAP_HISAT2_READS_PAIRED,
    MAP_HISAT2_READS_SINGLE,
    MAP_HISAT2_STRAND_F,
    MAP_HISAT2_STRAND_FR,
    MAP_HISAT2_STRAND_R,
    MAP_HISAT2_STRAND_RF,
    MAP_HISAT2_STRAND_UNSTRANDED,
    MAP_MINIMAP2_OUT_PAF,
    MAP_MINIMAP2_OUT_SAM,
    MAP_MINIMAP2_PRESET_KW,
    MAP_MINIMAP2_PRESET_LABELS,
    MAP_MINIMAP2_PRESETS,
    MAP_MULTI_ALL,
    MAP_MULTI_BEST,
    MAP_MULTI_NONE,
    MAP_NO_COVERAGE_N,
    MAP_NO_COVERAGE_REFERENCE,
    MAP_NO_COVERAGE_SKIP,
    MAP_RIGHT_CLICK_DIALOG,
    MAP_RIGHT_CLICK_RUN,
    MapOptions,
    opening_dialog_tab,
    persistable_dialog_tab,
)
from repeatviewer.views.tool_dialogs import (
    TOOL_PARAMS_DIALOG_WIDTH,
    TOOL_PARAMS_FIELD_WIDTH,
    TeMultiInputFilesWidget,
    _form_section_separator,
    dialog_should_run,
    executable_checkbox_row,
    make_cancel_save_run_bar,
    make_executable_status_label,
    make_process_output_box,
    resolve_system_executable,
    set_named_executable_status,
    set_path_checkbox_status,
    sync_form_field_column_spacers,
    tool_params_dialog_width,
    system_path_executable_label,
)
from repeatviewer.widgets.combo_boxes import NoWheelComboBox as QComboBox
from repeatviewer.widgets.help_button import HELP_BUTTON_SIZE, field_with_help
from repeatviewer.widgets.path_browse import PathBrowseRow
from repeatviewer.widgets.spin_boxes import NoWheelSpinBox as QSpinBox

_FASTA_SUFFIXES = {".fa", ".fasta", ".fna", ".fas", ".faa", ".fsa", ".ffn", ".seq"}

_MAP_INPUTS_HELP = (
    "<b>Use selected sequences</b> (default): map the FASTA / FASTQ files selected in "
    "the Files panel when you click Run.<br><br>"
    "<b>Define inputs</b>: choose different query files for this run only."
)
PRESET_INFO = (
    "minimap2 <b>-x</b> preset (applied before other options). Pick the preset that "
    "matches your data — there is no Auto mode.<br><br>"
    "<b>Long genomic reads</b><br>"
    "<b>map-ont</b>: Oxford Nanopore genomic reads (ordinary minimizers; default).<br>"
    "<b>map-pb</b>: PacBio CLR; uses homopolymer-compressed minimizers.<br>"
    "<b>map-hifi</b>: PacBio HiFi / CCS genomic reads.<br>"
    "<b>lr:hq</b>: Nanopore Q20 / other high-accuracy long genomic reads.<br>"
    "<b>map-iclr</b>: Illumina Complete Long Reads.<br><br>"
    "<b>RNA / cDNA (spliced)</b><br>"
    "<b>splice</b>: long mRNA / cDNA vs genome; strand unknown (Nanopore cDNA).<br>"
    "<b>splice:hq</b>: accurate long RNA (PacBio Iso-seq / Kinnex).<br>"
    "<b>splice:sr</b>: short RNA-seq (Illumina), spliced.<br><br>"
    "<b>Assemblies &amp; short DNA</b><br>"
    "<b>asm5 / asm10 / asm20</b>: assembly vs reference (~0.1 / 1 / 5% divergence).<br>"
    "<b>sr</b>: short genomic reads (no splicing)."
)
OUTPUT_FMT_INFO = (
    "<b>PAF</b> (default): text alignment. CIGAR (<code>-c</code>) and CS "
    "(<code>--cs=long</code>) are under Advanced Settings; the contig view needs CS.<br><br>"
    "<b>SAM</b>: <code>-a</code> with <code>--MD</code> for the contig view."
)
CIGAR_INFO = (
    "PAF only: add a CIGAR string with <code>-c</code> (base-level alignment). "
    "SAM always includes CIGAR. Leave on unless you want approximate mapping only."
)
CS_INFO = (
    "PAF only: add the long <code>cs</code> tag (<code>--cs=long</code>), which encodes "
    "mismatches and indels. The contig view needs this for PAF. SAM uses "
    "<code>--MD</code> instead."
)
KMER_INFO = (
    "minimap2 <b>-k</b> / <b>-w</b>. These change with <b>Map mode</b> (<code>-x</code>) "
    "(e.g. map-ont 15/10, splice 15/5, lr:hq 19/19, sr 21/11). Changing Map mode "
    "fills these numbers; you can still edit them. For noisy Nanopore Direct "
    "RNA-seq, try a smaller <b>-k</b> (e.g. 14) with <b>splice</b>."
)
MULTI_MM_INFO = (
    "When a query maps to more than one place on the reference, minimap2 ranks "
    "hits: the best is <b>primary</b>; others are <b>secondary</b> "
    "(e.g. near-identical copies).<br><br>"
    "<b>Primary only</b> (default): <code>--secondary=no</code> — keep only the "
    "best hit. One place per sequence.<br><br>"
    "<b>Up to 50 places (-N 50)</b>: <code>--secondary=yes -N 50 -p 0.95</code>. "
    "<code>-N</code> is how many secondary alignments to report (minimap2’s own "
    "default is 5; RepeatViewer allows up to 50). <code>-p 0.95</code> keeps a "
    "secondary only if its score is at least 95% of its primary.<br><br>"
    "<b>Do not place ambiguous</b>: after mapping, drop low-MAPQ primaries "
    "(unreliable best hit) instead of placing them."
)
MULTI_HS_INFO = (
    "HISAT2 multi-mapping (<code>-k</code>).<br><br>"
    "<b>Best place only</b> (default): <code>-k 1</code> — one primary alignment.<br><br>"
    "<b>Up to 20 places (-k 20)</b>: report up to 20 distinct alignments and mark "
    "extras as secondary (<code>--secondary</code>).<br><br>"
    "<b>Do not place ambiguous</b>: after mapping, drop low-MAPQ primaries."
)
HISAT2_LIBRARY_INFO = (
    "HISAT2 alignment mode.<br><br>"
    "<b>RNA-seq (spliced)</b> (default): allow splice junctions — for short RNA-seq "
    "or cDNA vs a genome.<br><br>"
    "<b>DNA-seq</b>: <code>--no-spliced-alignment</code> — genomic short reads "
    "(Illumina DNA) with no intron gaps."
)
HISAT2_READS_INFO = (
    "How selected FASTQ files are passed to HISAT2.<br><br>"
    "<b>Single-end</b> (default): <code>-U</code> — each file is unpaired "
    "(or several files are merged).<br><br>"
    "<b>Paired-end (Illumina R1/R2)</b>: <code>-1</code> / <code>-2</code>. "
    "Select both mates; names must differ by <code>R1</code>/<code>R2</code> "
    "(or <code>_1</code>/<code>_2</code>), e.g. "
    "<code>sample_R1.fastq.gz</code> + <code>sample_R2.fastq.gz</code>. "
    "Several samples are fine if every R1 has a matching R2.<br><br>"
    "Strandness FR/RF is the RNA kit protocol — it does <b>not</b> pair files."
)
HISAT2_STRAND_INFO = (
    "HISAT2 <code>--rna-strandness</code> (RNA-seq only; ignored for DNA-seq).<br><br>"
    "This is the library strand protocol, <b>not</b> how files are paired.<br><br>"
    "<b>Unstranded</b> (default): no strand constraint.<br>"
    "<b>F / R</b>: single-end forward or reverse (use with Reads → Single-end).<br>"
    "<b>FR / RF</b>: paired-end strand protocols (use with Reads → Paired-end; "
    "Illumina stranded kits are often RF)."
)
HISAT2_DTA_INFO = (
    "HISAT2 <code>--dta</code> (downstream transcriptome assembly). Longer anchors "
    "for splice discovery — recommended when feeding alignments to StringTie. "
    "RNA-seq mode only."
)
HISAT2_INTRON_INFO = (
    "HISAT2 <code>--min-intronlen</code> / <code>--max-intronlen</code> (RNA-seq only). "
    "Defaults 20 / 500000. Narrow the range for compact genomes or to avoid "
    "implausibly long introns."
)
HISAT2_SOFTCLIP_INFO = (
    "HISAT2 <code>--no-softclip</code>: disallow soft-clipping at read ends. "
    "Off by default."
)
_RIGHT_CLICK_ENGINE_HELP = (
    "When right-click runs the tool: <b>Auto</b> picks minimap2 or HISAT2 from the "
    "input; otherwise the saved engine from the last Minimap2 / HISAT2 tab you used."
)


def _is_fasta(path: Path) -> bool:
    name = path.name.lower()
    if name.endswith(".gz"):
        name = name[:-3]
    return Path(name).suffix in _FASTA_SUFFIXES


class MapReferenceDialog(QDialog):
    """Map selected sequences or reads to a reference (Align-style engine tabs)."""

    _DIALOG_WIDTH = TOOL_PARAMS_DIALOG_WIDTH
    _TAB_MINIMAP2 = 0
    _TAB_HISAT2 = 1
    _TAB_OPTIONS = 2
    # Longest Settings label — Inputs radios / Advanced checkboxes share this column.
    _SETTINGS_LABEL_REF = "Secondaries:"

    def __init__(
        self,
        *,
        initial: MapOptions | None = None,
        files: list[Path] | None = None,
        project_genome: Path | None = None,
        title: str = "Mapping",
        parent=None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setMinimumWidth(self._DIALOG_WIDTH)
        opts = initial or MapOptions()
        self._panel_files = [Path(p) for p in (files or [])]
        self._engine_compact_size = None
        self._project_genome = (
            project_genome.expanduser().resolve()
            if project_genome is not None and Path(project_genome).expanduser().is_file()
            else None
        )
        # Saved reference if set; otherwise the project genome from File → Define genome.
        initial_ref = str(opts.reference_fasta or "").strip()
        if not initial_ref and self._project_genome is not None:
            initial_ref = str(self._project_genome)

        # --- shared controls -------------------------------------------------
        self._reference = PathBrowseRow(
            self,
            initial=initial_ref or None,
            placeholder="Reference FASTA: click it in Files, or Browse",
            browse_title="Choose the reference FASTA",
            file_filter="FASTA (*.fa *.fasta *.fna *.fas *.fa.gz *.fasta.gz);;All files (*)",
        )
        # Expand with the dialog (same fill as Define-inputs / FastGA reference).
        self._reference.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        self._reference.line_edit().setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )

        self._input_group = QButtonGroup(self)
        self._input_selected = QRadioButton("Use selected sequences", self)
        self._input_define = QRadioButton("Define inputs", self)
        self._input_group.addButton(self._input_selected)
        self._input_group.addButton(self._input_define)
        self._input_selected.setChecked(True)

        def _make_map_define_inputs() -> TeMultiInputFilesWidget:
            return TeMultiInputFilesWidget(
                self,
                initial=None,
                allow_multiple=True,
                start_empty=False,
                placeholder="Query FASTA / FASTQ (.fa / .fq / .gz)",
                browse_title="Choose query files",
                file_filter=(
                    "Sequence (*.fa *.fasta *.fna *.fas *.fq *.fastq "
                    "*.fa.gz *.fasta.gz *.fq.gz *.fastq.gz);;All files (*)"
                ),
                path_mode="file",
                add_button_label="Add another input",
                add_button_tooltip=(
                    "Add another FASTA / FASTQ file. Click a file in the panel "
                    "to fill a row, or use Browse."
                ),
            )

        # Independent Define-inputs lists for Minimap2 and HISAT2.
        self._mm_define_inputs = _make_map_define_inputs()
        self._hs_define_inputs = _make_map_define_inputs()
        self._define_inputs = self._mm_define_inputs
        self._input_selected.toggled.connect(self._sync_inputs)
        self._input_define.toggled.connect(self._sync_inputs)
        self._mm_define_inputs.paths_changed.connect(self._on_define_paths_changed)
        self._hs_define_inputs.paths_changed.connect(self._on_define_paths_changed)

        self._hisat2_library = QComboBox(self)
        for label, value in (
            ("RNA-seq (spliced)", MAP_HISAT2_RNA),
            ("DNA-seq (--no-spliced-alignment)", MAP_HISAT2_DNA),
        ):
            self._hisat2_library.addItem(label, value)
        self._set_combo(
            self._hisat2_library,
            getattr(opts, "hisat2_library", MAP_HISAT2_RNA),
        )
        self._hisat2_library.currentIndexChanged.connect(self._sync_hisat2_library)

        self._hisat2_reads = QComboBox(self)
        for label, value in (
            ("Single-end (-U)", MAP_HISAT2_READS_SINGLE),
            ("Paired-end Illumina R1/R2 (-1/-2)", MAP_HISAT2_READS_PAIRED),
        ):
            self._hisat2_reads.addItem(label, value)
        self._set_combo(
            self._hisat2_reads,
            getattr(opts, "hisat2_reads", MAP_HISAT2_READS_SINGLE),
        )
        self._hisat2_reads.currentIndexChanged.connect(self._sync_hisat2_reads)

        self._hisat2_strand = QComboBox(self)
        for label, value in (
            ("Unstranded", MAP_HISAT2_STRAND_UNSTRANDED),
            ("F (single-end forward)", MAP_HISAT2_STRAND_F),
            ("R (single-end reverse)", MAP_HISAT2_STRAND_R),
            ("FR (paired)", MAP_HISAT2_STRAND_FR),
            ("RF (paired)", MAP_HISAT2_STRAND_RF),
        ):
            self._hisat2_strand.addItem(label, value)
        self._set_combo(
            self._hisat2_strand,
            getattr(opts, "hisat2_strandness", MAP_HISAT2_STRAND_UNSTRANDED),
        )

        self._multi_mm = QComboBox(self)
        for label, value in (
            ("Primary only", MAP_MULTI_BEST),
            ("Up to 50 places (-N 50)", MAP_MULTI_ALL),
            ("Do not place ambiguous", MAP_MULTI_NONE),
        ):
            self._multi_mm.addItem(label, value)
        self._set_combo(self._multi_mm, opts.multiple_hits)

        self._multi_hs = QComboBox(self)
        for label, value in (
            ("Best place only", MAP_MULTI_BEST),
            ("Up to 20 places (-k 20)", MAP_MULTI_ALL),
            ("Do not place ambiguous", MAP_MULTI_NONE),
        ):
            self._multi_hs.addItem(label, value)
        self._set_combo(self._multi_hs, opts.multiple_hits)

        def _sync_multi(source: QComboBox, target: QComboBox) -> None:
            target.blockSignals(True)
            self._set_combo(target, source.currentData())
            target.blockSignals(False)

        self._multi_mm.currentIndexChanged.connect(
            lambda _i: _sync_multi(self._multi_mm, self._multi_hs)
        )
        self._multi_hs.currentIndexChanged.connect(
            lambda _i: _sync_multi(self._multi_hs, self._multi_mm)
        )

        self._hisat2_dta = QCheckBox("Tailor for StringTie (--dta)", self)
        self._hisat2_dta.setChecked(bool(getattr(opts, "hisat2_dta", False)))
        self._hisat2_no_softclip = QCheckBox("Disallow soft-clipping (--no-softclip)", self)
        self._hisat2_no_softclip.setChecked(
            bool(getattr(opts, "hisat2_no_softclip", False))
        )
        self._hisat2_min_intron = QSpinBox(self)
        self._hisat2_min_intron.setRange(1, 1_000_000)
        self._hisat2_min_intron.setValue(int(getattr(opts, "hisat2_min_intron", 20) or 20))
        self._hisat2_max_intron = QSpinBox(self)
        self._hisat2_max_intron.setRange(1, 10_000_000)
        self._hisat2_max_intron.setValue(
            int(getattr(opts, "hisat2_max_intron", 500_000) or 500_000)
        )

        self._preset = QComboBox(self)
        for value in MAP_MINIMAP2_PRESETS:
            label = MAP_MINIMAP2_PRESET_LABELS.get(value, value)
            self._preset.addItem(label, value)
        preset = opts.minimap2_preset if opts.minimap2_preset in MAP_MINIMAP2_PRESETS else "map-ont"
        self._set_combo(self._preset, preset)
        self._preset.currentIndexChanged.connect(self._on_preset_changed)

        self._mm_format = QComboBox(self)
        self._mm_format.addItem("PAF", MAP_MINIMAP2_OUT_PAF)
        self._mm_format.addItem("SAM", MAP_MINIMAP2_OUT_SAM)
        self._set_combo(
            self._mm_format,
            getattr(opts, "minimap2_output", MAP_MINIMAP2_OUT_PAF),
        )
        self._mm_cigar = QCheckBox("Include CIGAR", self)
        self._mm_cigar.setChecked(bool(getattr(opts, "minimap2_include_cigar", True)))
        self._mm_cs = QCheckBox("Include CS tag", self)
        self._mm_cs.setChecked(bool(getattr(opts, "minimap2_include_cs", True)))
        self._mm_format.currentIndexChanged.connect(self._update_mm_output_extras)
        self._update_mm_output_extras()

        k0, w0 = self._kw_for_preset(preset, opts.minimap2_kmer, opts.minimap2_window)
        self._mm_kmer = QSpinBox(self)
        self._mm_kmer.setRange(1, 28)
        self._mm_kmer.setValue(k0)
        self._mm_window = QSpinBox(self)
        self._mm_window.setRange(1, 1000)
        self._mm_window.setValue(w0)

        self._lift = QCheckBox("Check the reference annotations in every sequence", self)
        self._lift.setChecked(opts.lift_annotations)

        self._threshold = QComboBox(self)
        for t in MAP_CONSENSUS_THRESHOLDS:
            label = (
                "Majority (most common base)"
                if t == 0
                else ("100 % (identical)" if t == 1 else f"{t:.0%}")
            )
            self._threshold.addItem(label, t)
        self._set_combo(self._threshold, opts.consensus_threshold)

        self._no_cov = QComboBox(self)
        for label, value in (
            ("N", MAP_NO_COVERAGE_N),
            ("Reference base", MAP_NO_COVERAGE_REFERENCE),
            ("Leave out", MAP_NO_COVERAGE_SKIP),
        ):
            self._no_cov.addItem(label, value)
        self._set_combo(self._no_cov, opts.no_coverage)

        self._min_var = QSpinBox(self)
        self._min_var.setRange(0, 100)
        self._min_var.setSuffix(" %")
        self._min_var.setSpecialValueText("Every difference")
        self._min_var.setValue(round(opts.min_variant_freq * 100))

        self._display = QSpinBox(self)
        self._display.setRange(50, 100_000)
        self._display.setSingleStep(500)
        self._display.setValue(opts.max_display_reads)

        self._save_consensus = QCheckBox("Consensus", self)
        self._save_consensus.setChecked(opts.save_consensus)
        self._save_mapped = QCheckBox("Placed sequences (turned to the reference strand)", self)
        self._save_mapped.setChecked(opts.save_mapped)
        self._save_unused = QCheckBox("Unused sequences, with the reason", self)
        self._save_unused.setChecked(opts.save_unused)
        self._save_variants = QCheckBox("Variants table", self)
        self._save_variants.setChecked(opts.save_variants)

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
        self._threads = QSpinBox(self)
        self._threads.setRange(1, 256)
        self._threads.setValue(max(1, opts.threads))

        for field in (
            self._hisat2_library,
            self._hisat2_reads,
            self._hisat2_strand,
            self._multi_mm,
            self._multi_hs,
            self._threshold,
            self._no_cov,
            self._preset,
            self._mm_format,
            self._threads,
            self._mm_kmer,
            self._mm_window,
            self._hisat2_min_intron,
            self._hisat2_max_intron,
        ):
            field.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        for field in (
            self._min_identity,
            self._min_aligned,
            self._min_mapq,
            self._min_var,
            self._display,
        ):
            field.setFixedWidth(120)

        # --- engine page (shared radios/Reference; per-engine Define inputs) -
        engine_page = QWidget(self)
        self._engine_page = engine_page
        engine_page.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )
        engine_root = QVBoxLayout(engine_page)
        engine_root.setContentsMargins(16, 10, 16, 10)
        engine_root.setSpacing(4)
        engine_root.setAlignment(Qt.AlignmentFlag.AlignTop)

        basic = QWidget(engine_page)
        self._engine_basic = basic
        basic.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum)
        form = QFormLayout(basic)
        form.setContentsMargins(0, 0, 0, 0)
        form.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
        form.setHorizontalSpacing(12)
        form.setVerticalSpacing(4)
        self._basic_form = form

        def _section(title: str) -> QLabel:
            return QLabel(f"<b>{title}</b>", engine_page)

        form.addRow(_section("Inputs"))
        inputs_box = QWidget(basic)
        inputs_layout = QVBoxLayout(inputs_box)
        inputs_layout.setContentsMargins(0, 0, 0, 0)
        inputs_layout.setSpacing(2)
        inputs_layout.addWidget(
            field_with_help(
                self._input_selected, "Map inputs", _MAP_INPUTS_HELP, stretch=True
            )
        )
        inputs_layout.addWidget(self._input_define)
        # Empty labels sized to Settings so radio dots line up with Map mode.
        self._inputs_radios_label = QLabel("")
        self._inputs_define_label = QLabel("")
        self._define_inputs_host = QWidget(basic)
        self._define_inputs_host.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Maximum
        )
        define_host_layout = QVBoxLayout(self._define_inputs_host)
        define_host_layout.setContentsMargins(0, 0, 0, 0)
        define_host_layout.setSpacing(0)
        define_host_layout.addWidget(self._mm_define_inputs)
        define_host_layout.addWidget(self._hs_define_inputs)
        form.addRow(self._inputs_radios_label, inputs_box)
        form.addRow(self._inputs_define_label, self._define_inputs_host)
        # Trailing help-sized spacer matches Define-inputs rows (Browse X align).
        self._reference_field = QWidget(basic)
        self._reference_field.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed
        )
        ref_row = QHBoxLayout(self._reference_field)
        ref_row.setContentsMargins(0, 0, 0, 0)
        ref_row.setSpacing(6)
        ref_row.addWidget(self._reference, stretch=1)
        ref_help_spacer = QWidget(self._reference_field)
        ref_help_spacer.setFixedSize(HELP_BUTTON_SIZE, HELP_BUTTON_SIZE)
        ref_row.addWidget(ref_help_spacer)
        form.addRow("Reference:", self._reference_field)

        form.addRow(_form_section_separator())
        form.addRow(_section("Settings"))
        # Minimap2-only settings (real minimap2 CLI).
        self._preset_field = field_with_help(self._preset, "Map mode (-x)", PRESET_INFO)
        form.addRow("Map mode:", self._preset_field)
        self._multi_mm_field = field_with_help(
            self._multi_mm, "Secondaries (-N)", MULTI_MM_INFO
        )
        # HISAT2-only settings (real HISAT2 CLI).
        self._hisat2_library_field = field_with_help(
            self._hisat2_library, "HISAT2 library", HISAT2_LIBRARY_INFO
        )
        form.addRow("Library:", self._hisat2_library_field)
        self._hisat2_reads_field = field_with_help(
            self._hisat2_reads, "HISAT2 reads", HISAT2_READS_INFO
        )
        form.addRow("Reads:", self._hisat2_reads_field)
        self._hisat2_strand_field = field_with_help(
            self._hisat2_strand, "RNA strandness", HISAT2_STRAND_INFO
        )
        form.addRow("Strandness:", self._hisat2_strand_field)
        self._multi_hs_field = field_with_help(
            self._multi_hs, "Multi-mapping", MULTI_HS_INFO
        )
        form.addRow("Secondaries:", self._multi_mm_field)
        form.addRow("Report:", self._multi_hs_field)
        form.addRow("Threads:", self._threads)

        # Minimap2 Output (FastGA-style Format + optional tags).
        self._mm_output_sep = _form_section_separator()
        form.addRow(self._mm_output_sep)
        self._mm_output_heading = _section("Output")
        form.addRow(self._mm_output_heading)
        self._mm_format_field = field_with_help(
            self._mm_format, "Output format", OUTPUT_FMT_INFO
        )
        form.addRow("Format:", self._mm_format_field)
        engine_root.addWidget(basic)

        advanced = QWidget(engine_page)
        self._advanced_panel = advanced
        advanced.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        adv = QFormLayout(advanced)
        adv.setContentsMargins(0, 4, 0, 0)
        adv.setLabelAlignment(
            Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter
        )
        adv.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.FieldsStayAtSizeHint)
        adv.setHorizontalSpacing(12)
        adv.setVerticalSpacing(4)
        adv.addRow(_form_section_separator())
        adv.addRow(_section("Advanced Settings"))
        self._kmer_field = field_with_help(self._mm_kmer, "minimap2 -k", KMER_INFO)
        self._window_field = field_with_help(self._mm_window, "minimap2 -w", KMER_INFO)
        adv.addRow("K-mer (-k):", self._kmer_field)
        adv.addRow("Window (-w):", self._window_field)
        # CIGAR / CS: needed for contig view defaults, rarely toggled → Advanced.
        self._mm_cigar_field = field_with_help(
            self._mm_cigar, "Include CIGAR", CIGAR_INFO
        )
        self._mm_cs_field = field_with_help(self._mm_cs, "Include CS tag", CS_INFO)
        self._mm_cigar_label = QLabel("")
        self._mm_cs_label = QLabel("")
        adv.addRow(self._mm_cigar_label, self._mm_cigar_field)
        adv.addRow(self._mm_cs_label, self._mm_cs_field)
        # HISAT2 Advanced (RNA-seq / softclip).
        self._hisat2_dta_field = field_with_help(
            self._hisat2_dta, "HISAT2 --dta", HISAT2_DTA_INFO
        )
        self._hisat2_dta_label = QLabel("")
        adv.addRow(self._hisat2_dta_label, self._hisat2_dta_field)
        self._hisat2_min_intron_field = field_with_help(
            self._hisat2_min_intron, "min intron", HISAT2_INTRON_INFO
        )
        self._hisat2_max_intron_field = field_with_help(
            self._hisat2_max_intron, "max intron", HISAT2_INTRON_INFO
        )
        adv.addRow("Min intron:", self._hisat2_min_intron_field)
        adv.addRow("Max intron:", self._hisat2_max_intron_field)
        self._hisat2_softclip_field = field_with_help(
            self._hisat2_no_softclip, "HISAT2 --no-softclip", HISAT2_SOFTCLIP_INFO
        )
        self._hisat2_softclip_label = QLabel("")
        adv.addRow(self._hisat2_softclip_label, self._hisat2_softclip_field)
        self._advanced_form = adv
        engine_root.addWidget(advanced)
        self._sync_field_column()
        # No stretch: dialog height follows content (avoids a tall empty window).

        # --- Options page ----------------------------------------------------
        self._right_click_buttons: dict[str, QRadioButton] = {}
        group = QButtonGroup(self)
        rc_box = QGroupBox("Right-click on Mapping", self)
        rcl = QVBoxLayout(rc_box)
        for label, value in (
            ("Map with saved settings", MAP_RIGHT_CLICK_RUN),
            ("Open this window", MAP_RIGHT_CLICK_DIALOG),
        ):
            radio = QRadioButton(label, self)
            radio.setChecked(opts.right_click_action == value)
            group.addButton(radio)
            self._right_click_buttons[value] = radio
            rcl.addWidget(radio)
        if not any(b.isChecked() for b in self._right_click_buttons.values()):
            self._right_click_buttons[MAP_RIGHT_CLICK_RUN].setChecked(True)

        self._rc_engine = QComboBox(self)
        self._rc_engine.addItem(mr.ENGINE_LABELS[MAP_ENGINE_AUTO], MAP_ENGINE_AUTO)
        self._rc_engine.addItem("minimap2", MAP_ENGINE_MINIMAP2)
        self._rc_engine.addItem("HISAT2", MAP_ENGINE_HISAT2)
        rc_engine = opts.engine if opts.engine in (
            MAP_ENGINE_AUTO,
            MAP_ENGINE_MINIMAP2,
            MAP_ENGINE_HISAT2,
        ) else MAP_ENGINE_AUTO
        self._set_combo(self._rc_engine, rc_engine)
        self._rc_engine.setFixedWidth(TOOL_PARAMS_FIELD_WIDTH)
        rc_engine_row = QWidget(self)
        rcel = QHBoxLayout(rc_engine_row)
        rcel.setContentsMargins(0, 0, 0, 0)
        rcel.addWidget(QLabel("Mapper for right-click / Save:", self))
        rcel.addWidget(
            field_with_help(
                self._rc_engine, "Right-click mapper", _RIGHT_CLICK_ENGINE_HELP
            )
        )
        rcel.addStretch(1)
        rcl.addWidget(rc_engine_row)

        self._use_system_path = QCheckBox(
            system_path_executable_label("Minimap2"), self
        )
        self._use_system_path.setChecked(opts.use_system_path)
        self._use_system_path.setToolTip(
            "Use minimap2 from PATH or the conda env in File → Settings → External tools "
            "instead of the bundled Minimap2."
        )
        self._minimap2_path_status = make_executable_status_label(self)
        self._hisat2_path_status = make_executable_status_label(self)
        self._use_system_path.toggled.connect(self._update_executable_status)
        exe_box = QGroupBox("Executables", self)
        exe_box.setSizePolicy(QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum)
        exl = QVBoxLayout(exe_box)
        exl.setContentsMargins(10, 8, 10, 8)
        exl.setSpacing(8)
        exl.addWidget(
            executable_checkbox_row(
                self, self._use_system_path, self._minimap2_path_status
            )
        )
        exl.addWidget(self._hisat2_path_status)

        process_box, self._show_process_window, self._close_process_window = (
            make_process_output_box(
                self,
                tool_name="the mapping",
                checked=opts.show_process_window,
                close_when_finished=bool(
                    getattr(opts, "close_process_window_when_finished", True)
                ),
            )
        )

        # RepeatViewer result options keep saved defaults (not shown — not mapper flags).
        for w in (
            self._min_identity,
            self._min_aligned,
            self._min_mapq,
            self._lift,
            self._threshold,
            self._no_cov,
            self._min_var,
            self._display,
            self._save_consensus,
            self._save_mapped,
            self._save_unused,
            self._save_variants,
        ):
            w.hide()

        options_page = QWidget(self)
        self._options_page = options_page
        options_page.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        ol = QVBoxLayout(options_page)
        ol.setContentsMargins(16, 14, 16, 14)
        ol.setSpacing(14)
        ol.setAlignment(Qt.AlignmentFlag.AlignTop)
        ol.addWidget(rc_box)
        ol.addWidget(exe_box)
        ol.addWidget(process_box)

        self._tab_bar = QTabBar(self)
        self._tab_bar.setExpanding(False)
        self._tab_bar.setDocumentMode(True)
        self._tab_bar.addTab("Minimap2")
        self._tab_bar.addTab("HISAT2")
        self._tab_bar.addTab("Options")
        self._stack = QStackedWidget(self)
        self._stack.setSizePolicy(
            QSizePolicy.Policy.Preferred, QSizePolicy.Policy.Maximum
        )
        self._stack.addWidget(engine_page)
        self._stack.addWidget(options_page)

        open_tab = opening_dialog_tab(
            getattr(opts, "last_dialog_tab", MAP_ENGINE_MINIMAP2),
            allowed=MAP_DIALOG_TABS,
            first_tab=(
                opts.engine
                if opts.engine in (MAP_ENGINE_MINIMAP2, MAP_ENGINE_HISAT2)
                else MAP_ENGINE_MINIMAP2
            ),
        )
        self._tab_bar.setCurrentIndex(
            {
                MAP_ENGINE_MINIMAP2: self._TAB_MINIMAP2,
                MAP_ENGINE_HISAT2: self._TAB_HISAT2,
                DIALOG_TAB_OPTIONS: self._TAB_OPTIONS,
            }.get(open_tab, self._TAB_MINIMAP2)
        )

        restore = QPushButton("Restore Defaults", self)
        restore.clicked.connect(self._restore_defaults)
        self._advanced_settings_btn = QPushButton("Advanced Settings ▸", self)
        self._advanced_settings_btn.setCheckable(True)
        self._advanced_settings_btn.setChecked(False)
        self._advanced_settings_btn.setCursor(Qt.CursorShape.PointingHandCursor)
        self._advanced_settings_btn.toggled.connect(self._on_advanced_toggled)
        self._button_box = make_cancel_save_run_bar(
            self,
            left_widgets=[restore, self._advanced_settings_btn],
        )

        layout = QVBoxLayout(self)
        layout.setContentsMargins(22, 18, 22, 16)
        layout.setSpacing(8)
        layout.setAlignment(Qt.AlignmentFlag.AlignTop)
        layout.addWidget(self._tab_bar, 0)
        layout.addWidget(self._stack, 0)
        layout.addWidget(self._button_box, 0)

        self._tab_bar.currentChanged.connect(self._on_tab_changed)
        self._set_advanced_visible(False)
        self._sync_inputs(resize=False)
        self._on_tab_changed(self._tab_bar.currentIndex())
        self.resize(self._DIALOG_WIDTH, 360)
        QTimer.singleShot(0, lambda: self._resize_for_advanced(False))

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._engine_compact_size = None
        self._sync_field_column()
        if self._tab_bar.currentIndex() == self._TAB_OPTIONS:
            self._update_executable_status()
        QTimer.singleShot(0, lambda: self._resize_for_advanced(False))

    def _sync_field_column(self) -> None:
        """Indent Inputs radios / Advanced checkboxes to the Settings field column."""
        if not hasattr(self, "_inputs_radios_label"):
            return
        sync_form_field_column_spacers(
            self,
            [
                self._inputs_radios_label,
                self._inputs_define_label,
                self._mm_cigar_label,
                self._mm_cs_label,
                self._hisat2_dta_label,
                self._hisat2_softclip_label,
            ],
            [self._basic_form, getattr(self, "_advanced_form", None)],
            ref_label=self._SETTINGS_LABEL_REF,
        )

    @staticmethod
    def _set_combo(combo: QComboBox, value) -> None:
        idx = combo.findData(value)
        if idx >= 0:
            combo.setCurrentIndex(idx)

    def _active_engine(self) -> str:
        if self._tab_bar.currentIndex() == self._TAB_HISAT2:
            return MAP_ENGINE_HISAT2
        if self._tab_bar.currentIndex() == self._TAB_MINIMAP2:
            return MAP_ENGINE_MINIMAP2
        # Options tab: keep last engine radio / combo
        return str(self._rc_engine.currentData() or MAP_ENGINE_AUTO)

    def current_dialog_tab(self) -> str:
        """Active Map window tab key (minimap2 / hisat2 / options)."""
        index = self._tab_bar.currentIndex()
        if index == self._TAB_HISAT2:
            return MAP_ENGINE_HISAT2
        if index == self._TAB_OPTIONS:
            return DIALOG_TAB_OPTIONS
        return MAP_ENGINE_MINIMAP2

    def _on_tab_changed(self, index: int) -> None:
        on_options = index == self._TAB_OPTIONS
        on_mm = index == self._TAB_MINIMAP2
        on_hs = index == self._TAB_HISAT2
        self._stack.setCurrentIndex(1 if on_options else 0)
        self._apply_engine_row_visibility(on_mm=on_mm, on_hs=on_hs)
        # Keep Options mapper in sync when leaving an engine tab.
        if index in (self._TAB_MINIMAP2, self._TAB_HISAT2):
            self._set_combo(self._rc_engine, self._active_engine())
        if on_options:
            self._advanced_settings_btn.setVisible(False)
            if self._advanced_settings_btn.isChecked():
                self._advanced_settings_btn.blockSignals(True)
                self._advanced_settings_btn.setChecked(False)
                self._advanced_settings_btn.setText("Advanced Settings ▸")
                self._advanced_settings_btn.blockSignals(False)
                self._set_advanced_visible(False)
            self._update_executable_status()
            if self.isVisible():
                self._fit_options()
            return
        # Swap which engine's Define-inputs list is shown.
        self._sync_inputs(resize=False)
        if on_hs:
            # Defer so the HISAT2 tab paints before any availability probe / reminder.
            QTimer.singleShot(0, self._remind_hisat2_if_needed)
        if self.isVisible():
            self._resize_for_advanced(self._advanced_settings_btn.isChecked())

    def _remind_hisat2_if_needed(self) -> None:
        """Show install/bind instructions when HISAT2 is not configured (FastGA style)."""
        if self._tab_bar.currentIndex() != self._TAB_HISAT2:
            return
        if mr.engine_available(MAP_ENGINE_HISAT2):
            return
        mr.show_hisat2_not_configured_dialog(self)

    def _update_executable_status(self, *_args) -> None:
        """Refresh Minimap2 PATH / HISAT2 lines on Options → Executables."""
        if not hasattr(self, "_minimap2_path_status"):
            return
        set_path_checkbox_status(
            self._minimap2_path_status,
            resolve_system_executable("minimap2"),
            checked=self._use_system_path.isChecked(),
            has_builtin=True,
        )
        set_named_executable_status(
            self._hisat2_path_status,
            "HISAT2",
            resolve_system_executable("hisat2"),
            missing_text="not configured",
        )
        # Grow / shrink the dialog for the new status line — do not squeeze groups.
        if (
            self.isVisible()
            and self._tab_bar.currentIndex() == self._TAB_OPTIONS
        ):
            QTimer.singleShot(0, self._fit_options)

    def _sync_hisat2_library(self, _index: int = 0) -> None:
        rna = self._hisat2_library.currentData() == MAP_HISAT2_RNA
        self._hisat2_strand.setEnabled(rna)
        self._hisat2_dta.setEnabled(rna)
        self._hisat2_min_intron.setEnabled(rna)
        self._hisat2_max_intron.setEnabled(rna)
        self._sync_hisat2_reads()

    def _sync_hisat2_reads(self, _index: int = 0) -> None:
        """Paired-end uses FR/RF strandness; single-end uses F/R."""
        paired = self._hisat2_reads.currentData() == MAP_HISAT2_READS_PAIRED
        current = self._hisat2_strand.currentData()
        for i in range(self._hisat2_strand.count()):
            value = self._hisat2_strand.itemData(i)
            if value in (MAP_HISAT2_STRAND_F, MAP_HISAT2_STRAND_R):
                self._hisat2_strand.model().item(i).setEnabled(not paired)
            elif value in (MAP_HISAT2_STRAND_FR, MAP_HISAT2_STRAND_RF):
                self._hisat2_strand.model().item(i).setEnabled(paired)
        if paired and current in (MAP_HISAT2_STRAND_F, MAP_HISAT2_STRAND_R):
            self._set_combo(self._hisat2_strand, MAP_HISAT2_STRAND_UNSTRANDED)
        if not paired and current in (MAP_HISAT2_STRAND_FR, MAP_HISAT2_STRAND_RF):
            self._set_combo(self._hisat2_strand, MAP_HISAT2_STRAND_UNSTRANDED)

    def _apply_engine_row_visibility(self, *, on_mm: bool, on_hs: bool) -> None:
        """Show Minimap2 vs HISAT2 Settings / Advanced / Output for the active tab."""
        basic = getattr(self, "_basic_form", None)
        adv = getattr(self, "_advanced_form", None)
        on_engine = on_mm or on_hs
        if basic is not None:
            basic.setRowVisible(self._preset_field, on_mm)
            basic.setRowVisible(self._multi_mm_field, on_mm)
            basic.setRowVisible(self._mm_output_sep, on_mm)
            basic.setRowVisible(self._mm_output_heading, on_mm)
            basic.setRowVisible(self._mm_format_field, on_mm)
            basic.setRowVisible(self._hisat2_library_field, on_hs)
            basic.setRowVisible(self._hisat2_reads_field, on_hs)
            basic.setRowVisible(self._hisat2_strand_field, on_hs)
            basic.setRowVisible(self._multi_hs_field, on_hs)
        if adv is not None:
            self._advanced_settings_btn.setVisible(
                on_engine and self._tab_bar.currentIndex() != self._TAB_OPTIONS
            )
            adv.setRowVisible(self._kmer_field, on_mm)
            adv.setRowVisible(self._window_field, on_mm)
            adv.setRowVisible(self._mm_cigar_field, on_mm)
            adv.setRowVisible(self._mm_cs_field, on_mm)
            adv.setRowVisible(self._hisat2_dta_field, on_hs)
            adv.setRowVisible(self._hisat2_min_intron_field, on_hs)
            adv.setRowVisible(self._hisat2_max_intron_field, on_hs)
            adv.setRowVisible(self._hisat2_softclip_field, on_hs)
            if on_hs:
                self._sync_hisat2_library()
            if on_engine and self._advanced_settings_btn.isChecked():
                # Keep advanced open when switching engines; just refresh rows.
                self._set_advanced_visible(True)

    def _kw_for_preset(
        self, preset: str, saved_k: int | None = None, saved_w: int | None = None
    ) -> tuple[int, int]:
        """Concrete -k/-w for the chosen ``-x`` preset."""
        table = MAP_MINIMAP2_PRESET_KW.get(preset, (15, 10))
        k = int(saved_k) if saved_k and int(saved_k) > 0 else table[0]
        w = int(saved_w) if saved_w and int(saved_w) > 0 else table[1]
        return k, w

    def _on_preset_changed(self, _index: int = 0) -> None:
        preset = str(self._preset.currentData() or "map-ont")
        k, w = self._kw_for_preset(preset)
        self._mm_kmer.blockSignals(True)
        self._mm_window.blockSignals(True)
        self._mm_kmer.setValue(k)
        self._mm_window.setValue(w)
        self._mm_kmer.blockSignals(False)
        self._mm_window.blockSignals(False)

    def _update_mm_output_extras(self, _index: int = 0) -> None:
        paf = str(self._mm_format.currentData()) == MAP_MINIMAP2_OUT_PAF
        self._mm_cigar.setEnabled(paf)
        self._mm_cs.setEnabled(paf)
        if not paf:
            self._mm_cigar.blockSignals(True)
            self._mm_cs.blockSignals(True)
            self._mm_cigar.setChecked(False)
            self._mm_cs.setChecked(False)
            self._mm_cigar.blockSignals(False)
            self._mm_cs.blockSignals(False)
        elif not self._mm_cigar.isChecked() and not self._mm_cs.isChecked():
            # Sensible PAF defaults for the contig view.
            self._mm_cigar.setChecked(True)
            self._mm_cs.setChecked(True)
    def _active_define_inputs(self) -> TeMultiInputFilesWidget:
        if self._tab_bar.currentIndex() == self._TAB_HISAT2:
            return self._hs_define_inputs
        return self._mm_define_inputs

    def _on_define_paths_changed(self) -> None:
        # Ignore path changes from the hidden engine's widget.
        if self.sender() is not None and self.sender() is not self._active_define_inputs():
            return
        if (
            self.isVisible()
            and self._tab_bar.currentIndex() != self._TAB_OPTIONS
            and not self._advanced_settings_btn.isChecked()
        ):
            self._engine_compact_size = None
            QTimer.singleShot(0, lambda: self._resize_for_advanced(False))

    def _sync_inputs(self, *, resize: bool = True) -> None:
        define = self._input_define.isChecked()
        active = self._active_define_inputs()
        self._define_inputs = active
        form = self._engine_basic.layout()
        if isinstance(form, QFormLayout):
            form.setRowVisible(self._define_inputs_host, define)
        self._define_inputs_host.setVisible(define)
        for widget in (self._mm_define_inputs, self._hs_define_inputs):
            show = define and widget is active
            widget.setVisible(show)
            if show:
                widget.setMinimumHeight(0)
                widget.setMaximumHeight(16777215)
            else:
                widget.setFixedHeight(0)
            widget.updateGeometry()
        if define:
            self._define_inputs_host.setMinimumHeight(0)
            self._define_inputs_host.setMaximumHeight(16777215)
        else:
            self._define_inputs_host.setFixedHeight(0)
        self._define_inputs_host.updateGeometry()
        self._engine_basic.updateGeometry()
        self._sync_field_column()
        if (
            resize
            and self.isVisible()
            and self._tab_bar.currentIndex() != self._TAB_OPTIONS
            and not self._advanced_settings_btn.isChecked()
        ):
            self._engine_compact_size = None
            QTimer.singleShot(0, lambda: self._resize_for_advanced(False))

    def _set_advanced_visible(self, visible: bool) -> None:
        panel = self._advanced_panel
        panel.setVisible(visible)
        if visible:
            panel.setMaximumHeight(16777215)
        else:
            panel.setMaximumHeight(0)
        panel.updateGeometry()

    def _on_advanced_toggled(self, checked: bool) -> None:
        self._advanced_settings_btn.setText(
            "Advanced Settings ▾" if checked else "Advanced Settings ▸"
        )
        if checked:
            self._engine_compact_size = self.size()
        self._set_advanced_visible(checked)
        self._resize_for_advanced(checked)

    def _dialog_chrome_height(self) -> int:
        lay = self.layout()
        margins = lay.contentsMargins() if lay is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 34
        spacing = lay.spacing() if lay is not None else 8
        # tab bar + stack + buttons → 2 gaps
        tab_h = max(24, self._tab_bar.sizeHint().height()) + 4
        buttons_h = self._button_box.sizeHint().height()
        return margin_h + (spacing * 2) + tab_h + buttons_h

    def _engine_page_height(self, *, expanded: bool) -> int:
        """Height of the engine stack page (basic ± Advanced), not Options."""
        root = self._engine_page.layout()
        margins = root.contentsMargins() if root is not None else None
        margin_h = (margins.top() + margins.bottom()) if margins is not None else 20
        spacing = root.spacing() if root is not None else 4
        self._engine_basic.ensurePolished()
        basic_h = max(
            self._engine_basic.minimumSizeHint().height(),
            self._engine_basic.sizeHint().height(),
        )
        if not expanded:
            return margin_h + basic_h + spacing
        self._advanced_panel.ensurePolished()
        extra = max(0, self._advanced_panel.sizeHint().height())
        return margin_h + basic_h + spacing + extra

    def _pin_stack_height(self, content_h: int) -> None:
        """Constrain the stack to the active page so Options cannot leave a gap."""
        h = max(1, content_h)
        self._stack.setMinimumHeight(0)
        self._stack.setMaximumHeight(h)
        page = self._stack.currentWidget()
        if page is not None:
            page.setMinimumHeight(0)
            page.setMaximumHeight(h)

    def _resize_for_advanced(self, expanded: bool) -> None:
        if self._tab_bar.currentIndex() == self._TAB_OPTIONS:
            self._fit_options()
            return
        # Clear prior height pins so newly added Define-input rows can grow.
        self._stack.setMaximumHeight(16777215)
        self._engine_page.setMaximumHeight(16777215)
        self._define_inputs_host.setMaximumHeight(16777215)
        active = self._active_define_inputs()
        if self._input_define.isChecked():
            active.setMinimumHeight(0)
            active.setMaximumHeight(16777215)
        self._engine_page.ensurePolished()
        self._engine_basic.ensurePolished()
        self._engine_basic.adjustSize()
        content_h = self._engine_page_height(expanded=expanded)
        chrome = self._dialog_chrome_height()
        width = tool_params_dialog_width(
            define_inputs=self._input_define.isChecked()
        )
        self.setMinimumWidth(width)
        screen = self.screen().availableGeometry() if self.screen() else None
        max_h = screen.height() - 40 if screen is not None else 10_000
        height = min(max(chrome + content_h, 200), max_h)
        self._pin_stack_height(content_h)
        self.setMinimumHeight(0)
        self.setMaximumHeight(height)
        self.resize(width, height)
        self.setMaximumHeight(16777215)
        if not expanded:
            self._engine_compact_size = self.size()

    def _fit_options(self) -> None:
        page = self._options_page
        # Clear a previous pin so newly shown status lines can increase sizeHint.
        self._stack.setMaximumHeight(16777215)
        page.setMaximumHeight(16777215)
        page.ensurePolished()
        page.adjustSize()
        root = page.layout()
        if root is not None:
            root.invalidate()
            root.activate()
        content_h = page.minimumSizeHint().expandedTo(page.sizeHint()).height()
        if root is not None:
            content_h = max(
                content_h,
                root.totalMinimumSize().height(),
                root.totalSizeHint().height(),
            )
        height = self._dialog_chrome_height() + content_h
        self._pin_stack_height(content_h)
        self.setMinimumHeight(0)
        self.setMaximumHeight(height)
        self.resize(self._DIALOG_WIDTH, height)
        self.setMaximumHeight(16777215)

    def _restore_defaults(self) -> None:
        d = MapOptions()
        self._input_selected.setChecked(True)
        self._mm_define_inputs.clear_paths()
        self._hs_define_inputs.clear_paths()
        self._set_combo(self._hisat2_library, d.hisat2_library)
        self._set_combo(self._hisat2_reads, d.hisat2_reads)
        self._set_combo(self._hisat2_strand, d.hisat2_strandness)
        self._hisat2_dta.setChecked(d.hisat2_dta)
        self._hisat2_no_softclip.setChecked(d.hisat2_no_softclip)
        self._hisat2_min_intron.setValue(d.hisat2_min_intron)
        self._hisat2_max_intron.setValue(d.hisat2_max_intron)
        self._set_combo(self._multi_mm, d.multiple_hits)
        self._set_combo(self._multi_hs, d.multiple_hits)
        self._lift.setChecked(d.lift_annotations)
        self._min_identity.setValue(round(d.min_identity * 100))
        self._min_aligned.setValue(d.min_aligned_bp)
        self._min_mapq.setValue(d.min_mapq)
        self._set_combo(self._preset, d.minimap2_preset)
        self._set_combo(self._mm_format, d.minimap2_output)
        self._mm_cigar.setChecked(d.minimap2_include_cigar)
        self._mm_cs.setChecked(d.minimap2_include_cs)
        self._update_mm_output_extras()
        self._on_preset_changed()
        self._threads.setValue(max(1, d.threads))
        self._set_combo(self._threshold, d.consensus_threshold)
        self._set_combo(self._no_cov, d.no_coverage)
        self._min_var.setValue(round(d.min_variant_freq * 100))
        self._display.setValue(d.max_display_reads)
        for box, value in (
            (self._save_consensus, d.save_consensus),
            (self._save_mapped, d.save_mapped),
            (self._save_unused, d.save_unused),
            (self._save_variants, d.save_variants),
        ):
            box.setChecked(value)
        self._right_click_buttons[d.right_click_action].setChecked(True)
        self._set_combo(self._rc_engine, MAP_ENGINE_AUTO)
        self._show_process_window.setChecked(d.show_process_window)
        self._close_process_window.setChecked(d.close_process_window_when_finished)
        self._use_system_path.setChecked(d.use_system_path)
        self._update_executable_status()
        self._tab_bar.setCurrentIndex(self._TAB_MINIMAP2)
        self._sync_hisat2_library()
        self._sync_inputs()

    def should_run(self) -> bool:
        return dialog_should_run(self)

    def input_source(self) -> str:
        if self._input_define.isChecked():
            return ALIGN_INPUT_DEFINE
        return ALIGN_INPUT_SELECTED

    def run_files(self, *, panel_files: list[Path] | None = None) -> list[Path]:
        if self.input_source() == ALIGN_INPUT_DEFINE:
            return self._active_define_inputs().paths()
        return [Path(p) for p in (panel_files if panel_files is not None else self._panel_files)]

    def files(self) -> list[Path]:
        return self.run_files()

    def apply_path_from_panel(self, path: Path) -> bool:
        """FASTA click: fills Define inputs when that mode is on, else Reference."""
        resolved = Path(path).expanduser()
        if not resolved.is_file():
            return False
        if self.input_source() == ALIGN_INPUT_DEFINE:
            return bool(self._active_define_inputs().accept_path(resolved))
        if not _is_fasta(resolved):
            return False
        self._reference.set_path(resolved)
        return True

    def options(self) -> MapOptions:
        right_click = next(
            (v for v, b in self._right_click_buttons.items() if b.isChecked()),
            MAP_RIGHT_CLICK_RUN,
        )
        ref = self._reference.path()
        # Engine: current engine tab when on Minimap2/HISAT2; Options combo otherwise.
        if self._tab_bar.currentIndex() in (self._TAB_MINIMAP2, self._TAB_HISAT2):
            engine = self._active_engine()
            self._set_combo(self._rc_engine, engine)
        else:
            engine = str(self._rc_engine.currentData() or MAP_ENGINE_AUTO)
        multi = (
            self._multi_hs.currentData()
            if self._tab_bar.currentIndex() == self._TAB_HISAT2
            else self._multi_mm.currentData()
        )
        return MapOptions(
            reference_fasta=str(ref) if ref else "",
            engine=engine,
            last_dialog_tab=persistable_dialog_tab(
                self.current_dialog_tab(),
                first_tab=MAP_ENGINE_MINIMAP2,
            ),
            multiple_hits=str(multi),
            hisat2_library=str(self._hisat2_library.currentData() or MAP_HISAT2_RNA),
            hisat2_reads=str(self._hisat2_reads.currentData() or MAP_HISAT2_READS_SINGLE),
            hisat2_strandness=str(self._hisat2_strand.currentData() or ""),
            hisat2_dta=bool(self._hisat2_dta.isChecked()),
            hisat2_min_intron=int(self._hisat2_min_intron.value()),
            hisat2_max_intron=int(self._hisat2_max_intron.value()),
            hisat2_no_softclip=bool(self._hisat2_no_softclip.isChecked()),
            minimap2_preset=str(self._preset.currentData()),
            minimap2_output=str(self._mm_format.currentData() or MAP_MINIMAP2_OUT_PAF),
            minimap2_include_cigar=bool(self._mm_cigar.isChecked()),
            minimap2_include_cs=bool(self._mm_cs.isChecked()),
            minimap2_kmer=int(self._mm_kmer.value()),
            minimap2_window=int(self._mm_window.value()),
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
