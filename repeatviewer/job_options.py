"""Per-job parameters passed from the UI into curation backends."""

from __future__ import annotations

import json
from dataclasses import dataclass, field


def _to_bool(value, default: bool = False) -> bool:
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


DEFAULT_TOOL_THREADS = 4


def _default_tool_threads() -> int:
    """Default Threads in tool dialogs (BLAST, MAFFT, TE-Aid, …)."""
    return DEFAULT_TOOL_THREADS


DEFAULT_EXTEND_LEFT = 1000
DEFAULT_EXTEND_RIGHT = 1000
# Preset lengths for the MSA Extend left/right editable dropdowns.
EXTEND_LENGTH_PRESETS: tuple[int, ...] = (
    0,
    100,
    200,
    300,
    500,
    1000,
    2000,
    5000,
    10_000,
)
EXTEND_RIGHT_CLICK_RUN = "run"
EXTEND_RIGHT_CLICK_DIALOG = "dialog"
EXTEND_RIGHT_CLICK_ACTIONS = (EXTEND_RIGHT_CLICK_RUN, EXTEND_RIGHT_CLICK_DIALOG)
# Default reporting E-value (shown as 0.0001 with four decimal places in the UI).
DEFAULT_REPORT_EVALUE = 1e-4
DEFAULT_BLAST_EVALUE = DEFAULT_REPORT_EVALUE
MIN_BLAST_EVALUE = 1e-4
DEFAULT_BLAST_THREADS = 1
DEFAULT_BLAST_MAX_MSA_LINES = 100
DEFAULT_BLAST_TOP_MSA_LINES = 100
DEFAULT_BLAST_WRITE_GFF = False
DEFAULT_BLAST_WRITE_TEXT = True
DEFAULT_BLAST_MINIMUM_BIT_SCORE = 0.0
DEFAULT_BLAST_WORD_SIZE = 11
DEFAULT_BLAST_GAP_OPEN = 5
DEFAULT_BLAST_GAP_EXTEND = 2
DEFAULT_BLAST_PENALTY = -3
DEFAULT_BLAST_PERCENT_IDENTITY = 0.0
DEFAULT_BLAST_QUERY_COVERAGE = 0.0
BLAST_PROGRAM_BLASTN = "blastn"
BLAST_PROGRAM_BLASTP = "blastp"
BLAST_PROGRAM_BLASTX = "blastx"
BLAST_PROGRAM_TBLASTN = "tblastn"
BLAST_PROGRAM_TBLASTX = "tblastx"
BLAST_PROGRAM_ONLINE = "online"
BLAST_PROGRAMS = (
    BLAST_PROGRAM_BLASTN,
    BLAST_PROGRAM_BLASTP,
    BLAST_PROGRAM_BLASTX,
    BLAST_PROGRAM_ONLINE,
)
# Last-open tab in BLAST parameters (program tabs + Options).
DIALOG_TAB_OPTIONS = "options"
BLAST_DIALOG_TABS = BLAST_PROGRAMS + (DIALOG_TAB_OPTIONS,)

# Process monitor: show while the tool runs; auto-close on success by default.
DEFAULT_SHOW_PROCESS_WINDOW = True
DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED = True


def opening_dialog_tab(
    saved: str | None,
    *,
    allowed: tuple[str, ...],
    first_tab: str,
) -> str:
    """Tab to select when opening a params dialog.

    Never reopen on Options — always land on a tool/Parameters tab. ``first_tab``
    is used when the saved value is Options or unknown.
    """
    if saved and saved != DIALOG_TAB_OPTIONS and saved in allowed:
        return saved
    return first_tab


def persistable_dialog_tab(current: str, *, first_tab: str) -> str:
    """Value to store for ``last_dialog_tab`` (never Options)."""
    if current == DIALOG_TAB_OPTIONS:
        return first_tab
    return current


# Local programs that search a protein BLAST database.
BLAST_PROGRAMS_PROTEIN_DB = (BLAST_PROGRAM_BLASTP, BLAST_PROGRAM_BLASTX)
# NCBI online BLAST suite (matches NCBI web program tabs).
BLAST_ONLINE_ALGORITHMS = (
    BLAST_PROGRAM_BLASTN,
    BLAST_PROGRAM_BLASTP,
    BLAST_PROGRAM_BLASTX,
    BLAST_PROGRAM_TBLASTN,
    BLAST_PROGRAM_TBLASTX,
)
DEFAULT_BLAST_ONLINE_PROGRAM = BLAST_PROGRAM_BLASTN
DEFAULT_BLAST_ONLINE_DATABASE = "core_nt"
DEFAULT_BLAST_ONLINE_EVALUE = DEFAULT_BLAST_EVALUE
DEFAULT_BLAST_ONLINE_HITLIST = 100
DEFAULT_BLAST_ONLINE_WRITE_TEXT = True
DEFAULT_BLASTP_EVALUE = DEFAULT_BLAST_EVALUE
DEFAULT_BLASTP_WORD_SIZE = 3
DEFAULT_BLASTP_GAP_OPEN = 11
DEFAULT_BLASTP_GAP_EXTEND = 1
DEFAULT_BLASTP_QUERY_COVERAGE = 0.0
DEFAULT_BLASTP_WRITE_TEXT = True
DEFAULT_BLASTP_WRITE_GFF = False

DEFAULT_BLASTX_EVALUE = DEFAULT_BLAST_EVALUE
DEFAULT_BLASTX_WORD_SIZE = 3
DEFAULT_BLASTX_GAP_OPEN = 11
DEFAULT_BLASTX_GAP_EXTEND = 1
DEFAULT_BLASTX_QUERY_COVERAGE = 0.0
DEFAULT_BLASTX_WRITE_TEXT = True
DEFAULT_BLASTX_WRITE_GFF = False

DEFAULT_TE_GENOME_MAP_BIN_BP = 100_000
TE_GENOME_MAP_DENSITY_HIT_COUNT = "hit_count"
TE_GENOME_MAP_DENSITY_TE_LENGTH = "te_length"
DEFAULT_TE_GENOME_MAP_DENSITY = TE_GENOME_MAP_DENSITY_TE_LENGTH
TE_GENOME_MAP_COLOR_IDENTITY = "identity"
TE_GENOME_MAP_COLOR_ORDER = "order"
TE_GENOME_MAP_COLOR_SUPERFAMILY = "superfamily"
TE_GENOME_MAP_COLOR_FAMILY = "family"
DEFAULT_TE_GENOME_MAP_COLOR_BY = TE_GENOME_MAP_COLOR_ORDER
DEFAULT_TE_GENOME_MAP_SPLIT_STRAND = False
DEFAULT_TE_GENOME_MAP_COLLAPSE_OVERLAPS = False
DEFAULT_TE_GENOME_MAP_SHOW_DIFF_TRACK = False
DEFAULT_TE_GENOME_MAP_SHOW_METH_DIFF_TRACK = False
DEFAULT_TE_GENOME_MAP_WINDOW_REF_LINE = True
CHROMOSOME_SORT_NATURAL = "natural"
CHROMOSOME_SORT_SIZE_DESC = "size_desc"
CHROMOSOME_SORT_SIZE_ASC = "size_asc"
DEFAULT_CHROMOSOME_SORT = CHROMOSOME_SORT_NATURAL
DEFAULT_TE_GENOME_MAP_METHYLATION_COVERAGE_MIN = 5
DEFAULT_ANNOTATION_DENSITY = TE_GENOME_MAP_DENSITY_TE_LENGTH
DEFAULT_ANNOTATION_FEATURE_TYPES = ("gene",)
TE_GENOME_MAP_DIALOG_TAB_TE = "te"
TE_GENOME_MAP_DIALOG_TAB_ANNOTATION = "annotation"
TE_GENOME_MAP_DIALOG_TAB_METHYLATION = "methylation"
# Legacy saved tab from the single-Parameters layout.
TE_GENOME_MAP_DIALOG_TAB_PARAMETERS = TE_GENOME_MAP_DIALOG_TAB_TE
TE_GENOME_MAP_DIALOG_TABS = (
    TE_GENOME_MAP_DIALOG_TAB_TE,
    TE_GENOME_MAP_DIALOG_TAB_ANNOTATION,
    TE_GENOME_MAP_DIALOG_TAB_METHYLATION,
    DIALOG_TAB_OPTIONS,
)


@dataclass
class TeGenomeMapOptions:
    """Genome density plot: bin window, density metric, and chromosome filter."""

    bin_size_bp: int = DEFAULT_TE_GENOME_MAP_BIN_BP
    density_mode: str = DEFAULT_TE_GENOME_MAP_DENSITY
    color_by: str = DEFAULT_TE_GENOME_MAP_COLOR_BY
    methylation_coverage_min: int = DEFAULT_TE_GENOME_MAP_METHYLATION_COVERAGE_MIN
    split_strand_tracks: bool = DEFAULT_TE_GENOME_MAP_SPLIT_STRAND
    collapse_overlaps: bool = DEFAULT_TE_GENOME_MAP_COLLAPSE_OVERLAPS
    show_te_diff_track: bool = DEFAULT_TE_GENOME_MAP_SHOW_DIFF_TRACK
    show_meth_diff_track: bool = DEFAULT_TE_GENOME_MAP_SHOW_METH_DIFF_TRACK
    show_window_ref_line: bool = DEFAULT_TE_GENOME_MAP_WINDOW_REF_LINE
    annotation_density_mode: str = DEFAULT_ANNOTATION_DENSITY
    annotation_feature_types: list[str] | None = None
    annotation_feature_types_by_file: dict[str, list[str]] | None = None
    chromosomes_only: list[str] | None = None
    superfamilies_only: list[str] | None = None
    chromosome_sort: str = DEFAULT_CHROMOSOME_SORT
    # When True, Run uses every track tab that has files. When False, Run uses
    # only the currently selected track tab (TE / Gene / Methylation).
    plot_all_tracks_together: bool = False
    # Reference genome FASTA for the plot. Empty means the project current genome.
    genome_path: str = ""
    show_process_window: bool = DEFAULT_SHOW_PROCESS_WINDOW
    close_process_window_when_finished: bool = DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED
    last_dialog_tab: str = TE_GENOME_MAP_DIALOG_TAB_TE


DEFAULT_TE_DIVERGENCE_BIN_WIDTH = 0.005


@dataclass
class TeDivergencePlotOptions:
    """TE divergence histogram: x-axis bin width on 0–1 divergence scale."""

    bin_width: float = DEFAULT_TE_DIVERGENCE_BIN_WIDTH


DEFAULT_TE_BENCHMARK_THREADS = 4


@dataclass
class TeConfusionMatrixBenchmarkOptions:
    """Confusion-matrix benchmark (query vs reference RepeatMasker .out on genome)."""

    tetrimmer_home: str = ""
    reference_out: str = ""
    less_classification: bool = False
    debug: bool = False


@dataclass
class TeLibraryBenchmarkOptions:
    """TE consensus-library benchmark (Perfect / Good / Present)."""

    query_fasta: str = ""
    standard_fasta: str = ""
    rm_out: str = ""
    run_repeatmasker: bool = True
    threads: int = DEFAULT_TE_BENCHMARK_THREADS


@dataclass
class ExtendOptions:
    left: int = DEFAULT_EXTEND_LEFT
    right: int = DEFAULT_EXTEND_RIGHT
    pan_genome_name: str | None = None
    # Single-genome mode: optional explicit FASTA to extend against. Empty means
    # use the project genome from Define Genome.
    genome_path: str | None = None
    # Toolbar right-click behaviour: run with saved settings, or open dialog.
    right_click_action: str = EXTEND_RIGHT_CLICK_RUN
    show_process_window: bool = DEFAULT_SHOW_PROCESS_WINDOW
    close_process_window_when_finished: bool = DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED


DEFAULT_TSD_LEFT_BOUNDARY = 0
DEFAULT_TSD_RIGHT_BOUNDARY = 0
DEFAULT_TSD_NT_RADIUS = 100
DEFAULT_TSD_MAX_SIM_THRESHOLD = 0.3
DEFAULT_TSD_NT_THRESHOLD = 0.95
DEFAULT_TSD_CLUSTER_EPS = 1.0


@dataclass
class TsdSearchOptions:
    """TSD search boundaries and filtering (MSA column positions)."""

    left_boundary: int = DEFAULT_TSD_LEFT_BOUNDARY
    right_boundary: int = DEFAULT_TSD_RIGHT_BOUNDARY
    nt_radius: int = DEFAULT_TSD_NT_RADIUS
    max_sim_threshold: float = DEFAULT_TSD_MAX_SIM_THRESHOLD
    nt_threshold: float = DEFAULT_TSD_NT_THRESHOLD
    cluster_eps: float = DEFAULT_TSD_CLUSTER_EPS


DEFAULT_ALIGN_THREADS = 5
ALIGN_MODE_COMBINED = "combined"
MAFFT_STRATEGY_AUTO = "auto"
MAFFT_STRATEGY_FFT_NS_1 = "fft_ns_1"
MAFFT_STRATEGY_FFT_NS_2 = "fft_ns_2"
MAFFT_STRATEGY_LINSI = "linsi"
MAFFT_STRATEGY_GINSI = "ginsi"
MAFFT_STRATEGY_EINSI = "einsi"
MAFFT_STRATEGY_CUSTOM = "custom"
MAFFT_STRATEGIES = (
    MAFFT_STRATEGY_AUTO,
    MAFFT_STRATEGY_FFT_NS_1,
    MAFFT_STRATEGY_FFT_NS_2,
    MAFFT_STRATEGY_LINSI,
    MAFFT_STRATEGY_GINSI,
    MAFFT_STRATEGY_EINSI,
    MAFFT_STRATEGY_CUSTOM,
)
MAFFT_SEQUENCE_AUTO = "auto"
MAFFT_SEQUENCE_NUCLEOTIDE = "nucleotide"
MAFFT_SEQUENCE_PROTEIN = "protein"
MAFFT_SEQUENCE_TYPES = (
    MAFFT_SEQUENCE_AUTO,
    MAFFT_SEQUENCE_NUCLEOTIDE,
    MAFFT_SEQUENCE_PROTEIN,
)
ALIGNER_MAFFT = "mafft"
ALIGNER_MUSCLE = "muscle"
ALIGNERS = (ALIGNER_MAFFT, ALIGNER_MUSCLE)
# FastGA sits in the Align window and can be the Align button's right-click
# action, but it is not an MSA aligner: never a combined / workflow aligner.
ALIGNER_FASTGA = "fastga"
RIGHT_CLICK_ALIGNERS = (ALIGNER_MAFFT, ALIGNER_MUSCLE, ALIGNER_FASTGA)
ALIGN_DIALOG_TABS = RIGHT_CLICK_ALIGNERS + (DIALOG_TAB_OPTIONS,)
MUSCLE_MODE_ALIGN = "align"
MUSCLE_MODE_SUPER5 = "super5"
MUSCLE_MODES = (MUSCLE_MODE_ALIGN, MUSCLE_MODE_SUPER5)
MUSCLE_PERM_NONE = "none"
MUSCLE_PERM_ABC = "abc"
MUSCLE_PERM_ACB = "acb"
MUSCLE_PERM_BCA = "bca"
MUSCLE_PERMS = (MUSCLE_PERM_NONE, MUSCLE_PERM_ABC, MUSCLE_PERM_ACB, MUSCLE_PERM_BCA)

# Where Align / FastGA query files come from (Inputs radios).
ALIGN_INPUT_SELECTED = "selected"  # Files-panel selection
ALIGN_INPUT_DEFINE = "define"      # explicit file/folder paths
ALIGN_INPUT_SOURCES = (ALIGN_INPUT_SELECTED, ALIGN_INPUT_DEFINE)
# FastGA keeps these aliases (same values).
FASTGA_INPUT_SELECTED = ALIGN_INPUT_SELECTED
FASTGA_INPUT_DEFINE = ALIGN_INPUT_DEFINE
FASTGA_INPUT_SOURCES = ALIGN_INPUT_SOURCES


@dataclass
class AlignOptions:
    """Combined-MSA aligner plus MAFFT and MUSCLE command-line parameters."""

    combined_aligner: str = ALIGNER_MAFFT
    num_threads: int = field(default_factory=_default_tool_threads)
    strategy: str = MAFFT_STRATEGY_FFT_NS_1
    sequence_type: str = MAFFT_SEQUENCE_AUTO
    retree: int = 1
    max_iterate: int = 0
    gap_open_penalty: float = 1.53
    offset: float = 0.0
    reorder: bool = False
    adjust_direction: bool = False
    muscle_mode: str = MUSCLE_MODE_ALIGN
    muscle_num_threads: int = field(default_factory=_default_tool_threads)
    muscle_perm: str = MUSCLE_PERM_NONE
    muscle_perturb: int = 0
    muscle_input_source: str = ALIGN_INPUT_SELECTED
    muscle_input_paths: list[str] = field(default_factory=list)
    mafft_input_source: str = ALIGN_INPUT_SELECTED
    mafft_input_paths: list[str] = field(default_factory=list)
    # Per-aligner "combine all selected files into one MSA" choice.
    mafft_combined: bool = False
    muscle_combined: bool = False
    # Which Align window tab was open last (mafft / muscle / fastga / options).
    last_dialog_tab: str = ALIGNER_MAFFT
    # Which aligner a right-click on the Align toolbar button runs.
    right_click_aligner: str = ALIGNER_MAFFT
    show_process_window: bool = DEFAULT_SHOW_PROCESS_WINDOW
    close_process_window_when_finished: bool = DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED


# FastGA output formats (CLI streaming flags, except 1aln which uses -1:file).
FASTGA_FORMAT_PAF = "paf"
FASTGA_FORMAT_PAFX = "pafx"
FASTGA_FORMAT_PAFM = "pafm"
FASTGA_FORMAT_PAFS = "pafs"
FASTGA_FORMAT_PAFS_LONG = "pafS"
FASTGA_FORMAT_PSL = "psl"
FASTGA_FORMAT_ALN = "1aln"
FASTGA_OUTPUT_FORMATS = (
    FASTGA_FORMAT_PAF,
    FASTGA_FORMAT_PAFX,
    FASTGA_FORMAT_PAFM,
    FASTGA_FORMAT_PAFS,
    FASTGA_FORMAT_PAFS_LONG,
    FASTGA_FORMAT_PSL,
    FASTGA_FORMAT_ALN,
)

FASTGA_RIGHT_CLICK_RUN = "run"
FASTGA_RIGHT_CLICK_DIALOG = "dialog"
FASTGA_RIGHT_CLICK_ACTIONS = (FASTGA_RIGHT_CLICK_RUN, FASTGA_RIGHT_CLICK_DIALOG)

DEFAULT_FASTGA_MIN_IDENTITY = 0.7
DEFAULT_FASTGA_MIN_LENGTH = 100
DEFAULT_FASTGA_SEED_FREQ = 10  # FastGA's own -f default
DEFAULT_FASTGA_CHAIN_COVERAGE = 85
DEFAULT_FASTGA_CHAIN_GAP = 1000
DEFAULT_FASTGA_SYMMETRIC = False  # FastGA default: off; enable in Advanced for repeats
DEFAULT_FASTGA_KEEP_INDEX = True
DEFAULT_FASTGA_SOFT_MASK = False
DEFAULT_FASTGA_SHOW_PROCESS_WINDOW = DEFAULT_SHOW_PROCESS_WINDOW

# What FastGA compares (FastGA tab → Mode).
FASTGA_MODE_SELF = "self"            # each selected file against itself
FASTGA_MODE_PAIRS = "pairs"          # every pair among the selected → one file
FASTGA_MODE_REFERENCE = "reference"  # each selected file against one reference
FASTGA_MODES = (FASTGA_MODE_SELF, FASTGA_MODE_PAIRS, FASTGA_MODE_REFERENCE)


@dataclass
class FastGaOptions:
    """FastGA whole-genome alignment.

    ``mode``: self (default) compares each selected FASTA with itself; pairs
    compares every unique pair and writes one combined file; reference runs
    each selected FASTA against ``reference_fasta`` (one file per query, or
    one combined file with ``reference_combined``).
    ``soft_mask``: seed only in uppercase bases; sent as -M or -N depending on
    the installed FastGA version.
    """

    output_format: str = FASTGA_FORMAT_PAF
    threads: int = field(default_factory=_default_tool_threads)
    min_identity: float = DEFAULT_FASTGA_MIN_IDENTITY
    min_length: int = DEFAULT_FASTGA_MIN_LENGTH
    seed_freq: int = DEFAULT_FASTGA_SEED_FREQ
    chain_coverage: int = DEFAULT_FASTGA_CHAIN_COVERAGE
    chain_gap: int = DEFAULT_FASTGA_CHAIN_GAP
    symmetric: bool = DEFAULT_FASTGA_SYMMETRIC
    keep_index: bool = DEFAULT_FASTGA_KEEP_INDEX
    soft_mask: bool = DEFAULT_FASTGA_SOFT_MASK
    mode: str = FASTGA_MODE_SELF
    reference_fasta: str = ""
    reference_combined: bool = False
    temp_dir: str = ""
    input_source: str = FASTGA_INPUT_SELECTED
    # Absolute FASTA files and/or folders (expanded at run time). Used when
    # input_source is FASTGA_INPUT_DEFINE.
    input_paths: list[str] = field(default_factory=list)
    right_click_action: str = FASTGA_RIGHT_CLICK_RUN
    show_process_window: bool = DEFAULT_FASTGA_SHOW_PROCESS_WINDOW
    close_process_window_when_finished: bool = DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED
    # Ephemeral: absolute path pairs chosen in the UI (None = every pair). Not saved.
    selected_pairs: list[tuple[str, str]] | None = None


def load_fastga_options() -> FastGaOptions:
    v = _load_settings_group("fastga")
    fmt = str(v.get("output_format", FASTGA_FORMAT_PAF))
    if fmt not in FASTGA_OUTPUT_FORMATS:
        fmt = FASTGA_FORMAT_PAF
    right_click_action = str(v.get("right_click_action", FASTGA_RIGHT_CLICK_RUN))
    if right_click_action not in FASTGA_RIGHT_CLICK_ACTIONS:
        right_click_action = FASTGA_RIGHT_CLICK_RUN
    mode = str(v.get("mode", "") or "")
    if mode not in FASTGA_MODES:
        # Settings from before the Mode choice: "Do each pair comparison".
        mode = FASTGA_MODE_PAIRS if _to_bool(v.get("pairwise"), False) else FASTGA_MODE_SELF
    input_source = str(v.get("input_source", FASTGA_INPUT_SELECTED) or "")
    if input_source not in FASTGA_INPUT_SOURCES:
        input_source = FASTGA_INPUT_SELECTED
    raw_inputs = str(v.get("input_paths", "") or "")
    input_paths = [p.strip() for p in raw_inputs.split("\n") if p.strip()]
    return FastGaOptions(
        output_format=fmt,
        threads=max(1, int(v.get("threads", _default_tool_threads()))),
        min_identity=max(
            0.0,
            min(1.0, float(v.get("min_identity", DEFAULT_FASTGA_MIN_IDENTITY))),
        ),
        min_length=max(1, int(v.get("min_length", DEFAULT_FASTGA_MIN_LENGTH))),
        seed_freq=max(
            1,
            DEFAULT_FASTGA_SEED_FREQ
            if int(v.get("seed_freq", DEFAULT_FASTGA_SEED_FREQ)) == 200
            else int(v.get("seed_freq", DEFAULT_FASTGA_SEED_FREQ)),
        ),  # 200 was RepeatViewer's old default; FastGA's -f default is 10
        chain_coverage=max(
            1, int(v.get("chain_coverage", DEFAULT_FASTGA_CHAIN_COVERAGE))
        ),
        chain_gap=max(1, int(v.get("chain_gap", DEFAULT_FASTGA_CHAIN_GAP))),
        symmetric=_to_bool(v.get("symmetric"), DEFAULT_FASTGA_SYMMETRIC),
        keep_index=_to_bool(v.get("keep_index"), DEFAULT_FASTGA_KEEP_INDEX),
        soft_mask=_to_bool(v.get("soft_mask"), DEFAULT_FASTGA_SOFT_MASK),
        mode=mode,
        reference_fasta=str(v.get("reference_fasta", "") or ""),
        reference_combined=_to_bool(v.get("reference_combined"), False),
        temp_dir=str(v.get("temp_dir", "") or ""),
        input_source=input_source,
        input_paths=input_paths,
        right_click_action=right_click_action,
        show_process_window=_to_bool(
            v.get("show_process_window"), DEFAULT_FASTGA_SHOW_PROCESS_WINDOW
        ),
        close_process_window_when_finished=_to_bool(
            v.get("close_process_window_when_finished"), DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED
        ),
    )


def save_fastga_options(opts: FastGaOptions) -> None:
    _save_settings_group(
        "fastga",
        {
            "output_format": opts.output_format,
            "threads": opts.threads,
            "min_identity": opts.min_identity,
            "min_length": opts.min_length,
            "seed_freq": opts.seed_freq,
            "chain_coverage": opts.chain_coverage,
            "chain_gap": opts.chain_gap,
            "symmetric": opts.symmetric,
            "keep_index": opts.keep_index,
            "soft_mask": opts.soft_mask,
            "mode": opts.mode,
            # Kept so an older RepeatViewer reading these settings still works.
            "pairwise": opts.mode == FASTGA_MODE_PAIRS,
            "reference_fasta": opts.reference_fasta,
            "reference_combined": opts.reference_combined,
            "temp_dir": opts.temp_dir,
            "input_source": opts.input_source,
            "input_paths": "\n".join(opts.input_paths),
            "right_click_action": opts.right_click_action,
            "show_process_window": opts.show_process_window,
            "close_process_window_when_finished": opts.close_process_window_when_finished,
        },
    )


@dataclass
class BlastOptions:
    program: str = BLAST_PROGRAM_BLASTN
    e_value: float = DEFAULT_BLAST_EVALUE
    num_threads: int = field(default_factory=_default_tool_threads)
    word_size: int = DEFAULT_BLAST_WORD_SIZE
    gap_open: int = DEFAULT_BLAST_GAP_OPEN
    gap_extend: int = DEFAULT_BLAST_GAP_EXTEND
    penalty: int = DEFAULT_BLAST_PENALTY
    percent_identity: float = DEFAULT_BLAST_PERCENT_IDENTITY
    query_coverage_hsp: float = DEFAULT_BLAST_QUERY_COVERAGE
    max_msa_lines: int = DEFAULT_BLAST_MAX_MSA_LINES
    top_msa_lines: int = DEFAULT_BLAST_TOP_MSA_LINES
    write_text: bool = DEFAULT_BLAST_WRITE_TEXT
    write_gff: bool = DEFAULT_BLAST_WRITE_GFF
    minimum_bit_score: float = DEFAULT_BLAST_MINIMUM_BIT_SCORE
    blast_database_name: str | None = None
    blastp_e_value: float = DEFAULT_BLASTP_EVALUE
    blastp_word_size: int = DEFAULT_BLASTP_WORD_SIZE
    blastp_gap_open: int = DEFAULT_BLASTP_GAP_OPEN
    blastp_gap_extend: int = DEFAULT_BLASTP_GAP_EXTEND
    blastp_query_coverage_hsp: float = DEFAULT_BLASTP_QUERY_COVERAGE
    blastp_database_name: str | None = None
    blastp_write_text: bool = DEFAULT_BLASTP_WRITE_TEXT
    blastp_write_gff: bool = DEFAULT_BLASTP_WRITE_GFF
    blastp_num_threads: int = field(default_factory=_default_tool_threads)
    blastx_e_value: float = DEFAULT_BLASTX_EVALUE
    blastx_word_size: int = DEFAULT_BLASTX_WORD_SIZE
    blastx_gap_open: int = DEFAULT_BLASTX_GAP_OPEN
    blastx_gap_extend: int = DEFAULT_BLASTX_GAP_EXTEND
    blastx_query_coverage_hsp: float = DEFAULT_BLASTX_QUERY_COVERAGE
    blastx_database_name: str | None = None
    blastx_write_text: bool = DEFAULT_BLASTX_WRITE_TEXT
    blastx_write_gff: bool = DEFAULT_BLASTX_WRITE_GFF
    blastx_num_threads: int = field(default_factory=_default_tool_threads)
    # NCBI online BLAST (Common URL API) — no local database / genome required.
    online_program: str = DEFAULT_BLAST_ONLINE_PROGRAM
    online_database: str = DEFAULT_BLAST_ONLINE_DATABASE
    online_e_value: float = DEFAULT_BLAST_ONLINE_EVALUE
    online_hitlist_size: int = DEFAULT_BLAST_ONLINE_HITLIST
    online_write_text: bool = DEFAULT_BLAST_ONLINE_WRITE_TEXT
    # Query Inputs (Align-style); Define paths are session-only.
    input_source: str = ALIGN_INPUT_SELECTED
    input_paths: list[str] = field(default_factory=list)
    # Which BLAST window tab was open last (blastn / … / online / options).
    last_dialog_tab: str = BLAST_PROGRAM_BLASTN
    # Which program a right-click on the BLAST toolbar button runs.
    right_click_program: str = BLAST_PROGRAM_BLASTN
    show_process_window: bool = DEFAULT_SHOW_PROCESS_WINDOW
    close_process_window_when_finished: bool = DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED

    @property
    def active_database_name(self) -> str | None:
        if self.program == BLAST_PROGRAM_ONLINE:
            return None
        if self.program == BLAST_PROGRAM_BLASTP:
            return self.blastp_database_name
        if self.program == BLAST_PROGRAM_BLASTX:
            return self.blastx_database_name
        return self.blast_database_name

    @property
    def is_online(self) -> bool:
        return self.program == BLAST_PROGRAM_ONLINE


DEFAULT_HMMER_EVALUE = DEFAULT_REPORT_EVALUE
DEFAULT_HMMER_DOMAIN_EVALUE = DEFAULT_REPORT_EVALUE
DEFAULT_HMMER_MAX_HITS = 100
DEFAULT_HMMER_MINIMUM_BIT_SCORE = 0.0
DEFAULT_HMMER_WRITE_TEXT = True
# 0 = no per-sequence ORF cap (protein HMM path).
DEFAULT_HMMER_MAX_ORFS = 0
HMMER_DIALOG_TAB_PARAMETERS = "parameters"
HMMER_DIALOG_TABS = (HMMER_DIALOG_TAB_PARAMETERS, DIALOG_TAB_OPTIONS)
HMMER_RIGHT_CLICK_RUN = "run"
HMMER_RIGHT_CLICK_DIALOG = "dialog"
HMMER_RIGHT_CLICK_ACTIONS = (HMMER_RIGHT_CLICK_RUN, HMMER_RIGHT_CLICK_DIALOG)


@dataclass
class HmmerOptions:
    """Tools → HMMER search against a project HMM database (pyhmmer)."""

    database_name: str = ""
    e_value: float = DEFAULT_HMMER_EVALUE
    domain_e_value: float = DEFAULT_HMMER_DOMAIN_EVALUE
    num_threads: int = field(default_factory=_default_tool_threads)
    max_hits: int = DEFAULT_HMMER_MAX_HITS
    write_text: bool = DEFAULT_HMMER_WRITE_TEXT
    minimum_bit_score: float = DEFAULT_HMMER_MINIMUM_BIT_SCORE
    max_orfs_per_sequence: int = DEFAULT_HMMER_MAX_ORFS
    input_source: str = ALIGN_INPUT_SELECTED
    input_paths: list[str] = field(default_factory=list)
    last_dialog_tab: str = HMMER_DIALOG_TAB_PARAMETERS
    show_process_window: bool = DEFAULT_SHOW_PROCESS_WINDOW
    close_process_window_when_finished: bool = DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED
    right_click_action: str = HMMER_RIGHT_CLICK_RUN


DEFAULT_CONSENSUS_THRESHOLD = 0.8
DEFAULT_CONSENSUS_AMBIGUOUS = "N"
CONSENSUS_AMBIGUOUS_CHOICES = ("N", "X")
CONSENSUS_RIGHT_CLICK_RUN = "run"
CONSENSUS_RIGHT_CLICK_DIALOG = "dialog"
CONSENSUS_RIGHT_CLICK_ACTIONS = (CONSENSUS_RIGHT_CLICK_RUN, CONSENSUS_RIGHT_CLICK_DIALOG)


@dataclass
class ConsensusOptions:
    threshold: float = DEFAULT_CONSENSUS_THRESHOLD
    ambiguous: str = DEFAULT_CONSENSUS_AMBIGUOUS
    sequence_type: str = MAFFT_SEQUENCE_AUTO
    right_click_action: str = CONSENSUS_RIGHT_CLICK_RUN


def normalize_consensus_sequence_type(raw: object) -> str:
    """Return ``auto``, ``nucleotide``, or ``protein`` for consensus generation."""
    value = str(raw or MAFFT_SEQUENCE_AUTO).strip().lower()
    if value in {"nuc", "nucleotide", "dna", "rna"}:
        return MAFFT_SEQUENCE_NUCLEOTIDE
    if value in {"aa", "protein", "amino", "aminoacid", "amino_acid"}:
        return MAFFT_SEQUENCE_PROTEIN
    if value in MAFFT_SEQUENCE_TYPES:
        return value
    return MAFFT_SEQUENCE_AUTO


def normalize_consensus_ambiguous(raw: object) -> str:
    """Return a supported ambiguous base (``N`` or ``X``)."""
    amb = str(raw or DEFAULT_CONSENSUS_AMBIGUOUS).upper()[:1]
    return amb if amb in CONSENSUS_AMBIGUOUS_CHOICES else DEFAULT_CONSENSUS_AMBIGUOUS


DEFAULT_FLEXIDOT_WORDSIZE = 5
DEFAULT_FLEXIDOT_LCS_SHADING = True
DEFAULT_FLEXIDOT_USE_GFF = False
DEFAULT_FLEXIDOT_OUTPUT_FORMAT = "png"
DEFAULT_FLEXIDOT_SUBSTITUTIONS = 0
DEFAULT_FLEXIDOT_WOBBLE_CONVERSION = False
DEFAULT_FLEXIDOT_MAX_N_PERCENTAGE = 10.0
DEFAULT_FLEXIDOT_INCLUDE_REVERSE = True
DEFAULT_FLEXIDOT_ONLY_VS_FIRST = False
DEFAULT_FLEXIDOT_LINE_WIDTH = 1.0
DEFAULT_FLEXIDOT_LENGTH_SCALING = False
DEFAULT_FLEXIDOT_COMBINED = True
DEFAULT_FLEXIDOT_SEQUENCE_TYPE = "nuc"
FLEXIDOT_SEQUENCE_NUC = "nuc"
FLEXIDOT_SEQUENCE_AA = "aa"
FLEXIDOT_SEQUENCE_TYPES = (FLEXIDOT_SEQUENCE_NUC, FLEXIDOT_SEQUENCE_AA)
# Flexidot comparison modes (-m): "0"=self, "1"=pairwise, "2"=all-vs-all.
# Empty tuple means automatic via ``auto_flexidot_mode`` (self / pairwise /
# all-vs-all from sequence count). Not exposed in the Flexidot dialog.
FLEXIDOT_MODE_CHOICES = ("0", "1", "2")
DEFAULT_FLEXIDOT_MODES: tuple[str, ...] = ()
FLEXIDOT_RIGHT_CLICK_RUN = "run"
FLEXIDOT_RIGHT_CLICK_DIALOG = "dialog"
FLEXIDOT_RIGHT_CLICK_ACTIONS = (FLEXIDOT_RIGHT_CLICK_RUN, FLEXIDOT_RIGHT_CLICK_DIALOG)


def normalize_flexidot_sequence_type(raw: object) -> str:
    """Return ``nuc`` or ``aa`` (Flexidot ``-t``); default DNA."""
    value = str(raw or DEFAULT_FLEXIDOT_SEQUENCE_TYPE).strip().lower()
    if value in {"aa", "protein", "amino", "aminoacid", "amino_acid"}:
        return FLEXIDOT_SEQUENCE_AA
    return FLEXIDOT_SEQUENCE_NUC


def auto_flexidot_mode(sequence_count: int) -> str:
    """Pick Flexidot ``-m`` from sequence count: self, pairwise, or all-vs-all."""
    if sequence_count <= 1:
        return "0"
    if sequence_count == 2:
        return "1"
    return "2"


def parse_flexidot_modes(raw: object) -> tuple[str, ...]:
    """Normalise a stored/typed modes value into an ordered ``0/1/2`` tuple."""
    tokens = {t.strip() for t in str(raw or "").replace(" ", "").split(",")}
    return tuple(m for m in FLEXIDOT_MODE_CHOICES if m in tokens)
FLEXIDOT_OUTPUT_HTML = "html"
FLEXIDOT_OUTPUT_PNG = "png"
FLEXIDOT_OUTPUT_PDF = "pdf"
FLEXIDOT_OUTPUT_SVG = "svg"
FLEXIDOT_OUTPUT_FORMATS = (
    FLEXIDOT_OUTPUT_HTML,
    FLEXIDOT_OUTPUT_PNG,
    FLEXIDOT_OUTPUT_PDF,
    FLEXIDOT_OUTPUT_SVG,
)
MIN_FLEXIDOT_WORDSIZE = 3
MAX_FLEXIDOT_WORDSIZE = 100


def normalize_flexidot_output_format(value: str) -> str:
    """Return one of html / png / pdf / svg for Flexidot ``-f``."""
    v = str(value).strip().lower()
    if v in FLEXIDOT_OUTPUT_FORMATS:
        return v
    return DEFAULT_FLEXIDOT_OUTPUT_FORMAT


@dataclass
class FlexidotOptions:
    """Saved matching, display, and output parameters for Flexidot."""

    wordsize: int = DEFAULT_FLEXIDOT_WORDSIZE
    # Flexidot -t: ``nuc`` (DNA, default) or ``aa`` (protein).
    sequence_type: str = DEFAULT_FLEXIDOT_SEQUENCE_TYPE
    lcs_shading: bool = DEFAULT_FLEXIDOT_LCS_SHADING
    use_gff: bool = DEFAULT_FLEXIDOT_USE_GFF
    output_format: str = DEFAULT_FLEXIDOT_OUTPUT_FORMAT
    substitution_count: int = DEFAULT_FLEXIDOT_SUBSTITUTIONS
    wobble_conversion: bool = DEFAULT_FLEXIDOT_WOBBLE_CONVERSION
    max_n_percentage: float = DEFAULT_FLEXIDOT_MAX_N_PERCENTAGE
    include_reverse: bool = DEFAULT_FLEXIDOT_INCLUDE_REVERSE
    only_vs_first: bool = DEFAULT_FLEXIDOT_ONLY_VS_FIRST
    line_width: float = DEFAULT_FLEXIDOT_LINE_WIDTH
    # Scale pairwise/all-vs-all plots by sequence length (flexidot -L).
    length_scaling: bool = DEFAULT_FLEXIDOT_LENGTH_SCALING
    # Combine all selected files into one dotplot; also the right-click mode.
    # When False, each selected file is plotted separately.
    combined: bool = DEFAULT_FLEXIDOT_COMBINED
    # Explicit flexidot comparison modes (-m); empty tuple = auto (see above).
    modes: tuple[str, ...] = DEFAULT_FLEXIDOT_MODES
    right_click_action: str = FLEXIDOT_RIGHT_CLICK_RUN
    show_process_window: bool = DEFAULT_SHOW_PROCESS_WINDOW
    close_process_window_when_finished: bool = DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED


TE_ANNOTATION_RIGHT_CLICK_RUN = "run"
TE_ANNOTATION_RIGHT_CLICK_DIALOG = "dialog"
TE_ANNOTATION_RIGHT_CLICK_ACTIONS = (
    TE_ANNOTATION_RIGHT_CLICK_RUN,
    TE_ANNOTATION_RIGHT_CLICK_DIALOG,
)
DEFAULT_TE_ANNOTATION_AUTO_MODES = ("LTR", "TIR", "tandem repeat", "ORF")


@dataclass
class TeAnnotationJobOptions:
    """Saved module selection and toolbar right-click behaviour for TE annotation."""

    # AutoAnnotationMode.value strings; empty means none selected.
    auto_modes: tuple[str, ...] = DEFAULT_TE_ANNOTATION_AUTO_MODES
    # Selected similarity database names. None means “all currently defined”.
    similarity_database_names: tuple[str, ...] | None = None
    right_click_action: str = TE_ANNOTATION_RIGHT_CLICK_RUN
    show_process_window: bool = DEFAULT_SHOW_PROCESS_WINDOW
    close_process_window_when_finished: bool = DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED


DEFAULT_TEAID_EVALUE = DEFAULT_REPORT_EVALUE
DEFAULT_TEAID_THREADS = DEFAULT_TOOL_THREADS
DEFAULT_TEAID_SELF_BLAST_WORD_SIZE = 6
# Text written on the consensus (x-)axis of the TE-Aid panels. Configurable so
# non-TE users (e.g. tandem repeats) can relabel it.
DEFAULT_TEAID_X_AXIS_LABEL = "Input sequence (bp)"
MIN_TEAID_SELF_BLAST_WORD_SIZE = 4
MAX_TEAID_SELF_BLAST_WORD_SIZE = 25
# Genome-BLAST hit length bands (% of consensus / input length).
DEFAULT_TEAID_HIT_LEN_HIGH_PCT = 80.0
DEFAULT_TEAID_HIT_LEN_MID_PCT = 50.0
DEFAULT_TEAID_HIT_COLOR_HIGH = "#ff0000"
DEFAULT_TEAID_HIT_COLOR_MID = "#1f77b4"
DEFAULT_TEAID_HIT_LINE_TRANSPARENCY = 0.2
TEAID_RIGHT_CLICK_RUN = "run"
TEAID_RIGHT_CLICK_DIALOG = "dialog"
TEAID_RIGHT_CLICK_ACTIONS = (TEAID_RIGHT_CLICK_RUN, TEAID_RIGHT_CLICK_DIALOG)


def clamp_teaid_self_blast_word_size(requested: int) -> int:
    """Clamp TE-Aid self-BLAST ``-word_size`` to 4–25."""
    return max(
        MIN_TEAID_SELF_BLAST_WORD_SIZE,
        min(MAX_TEAID_SELF_BLAST_WORD_SIZE, int(requested)),
    )


def clamp_teaid_hit_len_pct(high: float, mid: float) -> tuple[float, float]:
    """Keep high/mid in 0–100 and ensure mid < high."""
    high_pct = max(0.0, min(100.0, float(high)))
    mid_pct = max(0.0, min(100.0, float(mid)))
    if mid_pct >= high_pct:
        mid_pct = max(0.0, high_pct - 1.0) if high_pct > 0 else 0.0
    return high_pct, mid_pct


def clamp_teaid_hit_line_transparency(value: float) -> float:
    """Keep hit-line transparency in 0–1 (0 = opaque, 1 = invisible)."""
    return max(0.0, min(1.0, float(value)))


def configure_report_evalue_spinbox(spin) -> None:
    """Four decimal places (e.g. 0.0001), minimum 0.0001 for BLAST+."""
    spin.setDecimals(4)
    spin.setRange(MIN_BLAST_EVALUE, 1.0)
    spin.setSingleStep(0.0001)


def normalize_blast_evalue(e_value: float, *, default: float = DEFAULT_BLAST_EVALUE) -> float:
    """BLAST+ requires a strictly positive -evalue."""
    try:
        value = float(e_value)
    except (TypeError, ValueError):
        return default
    if value <= 0.0:
        return default
    return max(value, MIN_BLAST_EVALUE)


@dataclass
class TeaAidOptions:
    e_value: float = DEFAULT_TEAID_EVALUE
    self_blast_word_size: int = DEFAULT_TEAID_SELF_BLAST_WORD_SIZE
    num_threads: int = field(default_factory=_default_tool_threads)
    x_axis_label: str = DEFAULT_TEAID_X_AXIS_LABEL
    # Genome BLAST hit coloring by span as % of consensus length.
    hit_len_high_pct: float = DEFAULT_TEAID_HIT_LEN_HIGH_PCT
    hit_len_mid_pct: float = DEFAULT_TEAID_HIT_LEN_MID_PCT
    hit_color_high: str = DEFAULT_TEAID_HIT_COLOR_HIGH
    hit_color_mid: str = DEFAULT_TEAID_HIT_COLOR_MID
    # 0 = opaque, 1 = fully transparent; applies to all genome BLAST hit lines.
    hit_line_transparency: float = DEFAULT_TEAID_HIT_LINE_TRANSPARENCY
    right_click_action: str = TEAID_RIGHT_CLICK_RUN
    # Single-genome mode: empty means use the project genome from Define Genome.
    genome_path: str | None = None
    # TE-Aid Mode: pan-genome run uses pan_genomes (name, path) list.
    use_pan_genome: bool = False
    pan_genomes: list[tuple[str, str]] = field(default_factory=list)
    show_process_window: bool = DEFAULT_SHOW_PROCESS_WINDOW
    close_process_window_when_finished: bool = DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED


# remove_gaps_with_similarity_check defaults (TEtrimmerGUI)
DEFAULT_GAP_COL_GAP_THRESHOLD = 0.8
DEFAULT_GAP_COL_SIMI_GAP_THR = 0.4
DEFAULT_GAP_COL_SIMILARITY_THR = 0.7
DEFAULT_GAP_COL_MIN_NUCLEOTIDE = 5


@dataclass
class GappyColumnOptions:
    gap_threshold: float = DEFAULT_GAP_COL_GAP_THRESHOLD
    simi_check_gap_thr: float = DEFAULT_GAP_COL_SIMI_GAP_THR
    similarity_thr: float = DEFAULT_GAP_COL_SIMILARITY_THR
    min_nucleotide: int = DEFAULT_GAP_COL_MIN_NUCLEOTIDE


DEFAULT_DIV_THRESHOLD = 0.8
DEFAULT_DIV_WINDOW_SIZE = 20
DEFAULT_GAP_ROW_GAP_THRESHOLD = 0.05
DEFAULT_GAP_ROW_WINDOW_SIZE = 300


@dataclass
class DivergenceRowOptions:
    """Parameters for clean row with divergence (crop_end_div)."""

    threshold: float = DEFAULT_DIV_THRESHOLD
    window_size: int = DEFAULT_DIV_WINDOW_SIZE


@dataclass
class GapRowOptions:
    """Parameters for clean row with gap (crop_end_gap)."""

    gap_threshold: float = DEFAULT_GAP_ROW_GAP_THRESHOLD
    window_size: int = DEFAULT_GAP_ROW_WINDOW_SIZE


DEFAULT_MSA_CLEAN_METHOD = "msa_clean_gappy"
DEFAULT_MSA_CLEAN_ORDER = (
    "msa_clean_gappy",
    "msa_clean_div",
    "msa_clean_gap",
)
MSA_CLEAN_RIGHT_CLICK_COMBINED = "msa_clean"
MSA_CLEAN_RIGHT_CLICK_GAPPY = "msa_clean_gappy"
MSA_CLEAN_RIGHT_CLICK_DIV = "msa_clean_div"
MSA_CLEAN_RIGHT_CLICK_GAP = "msa_clean_gap"
MSA_CLEAN_RIGHT_CLICK_DIALOG = "dialog"
MSA_CLEAN_RIGHT_CLICK_ACTIONS = (
    MSA_CLEAN_RIGHT_CLICK_COMBINED,
    MSA_CLEAN_RIGHT_CLICK_GAPPY,
    MSA_CLEAN_RIGHT_CLICK_DIV,
    MSA_CLEAN_RIGHT_CLICK_GAP,
    MSA_CLEAN_RIGHT_CLICK_DIALOG,
)

MSA_CLEAN_METHOD_DESCRIPTIONS: dict[str, str] = {
    "msa_clean_gappy": (
        "Remove gappy columns: if a single column in the multiple sequence "
        "alignment has a gap proportion larger than the gap threshold and the "
        "proportion of the most common nucleotide in this column is less than "
        "the similarity threshold, this column will be removed.\n\n"
        "Example: if a column has less than 40% gaps (similarity-check gap "
        "threshold) and the portion of T (or any other nucleotide) is greater "
        "than 70% (similarity threshold), this column will be kept.\n\n"
        "Columns with gap proportion above the gap threshold are removed "
        "directly. Columns with fewer than the minimum nucleotides are also "
        "removed."
    ),
    "msa_clean_div": (
        "Clean rows with divergence converts each nucleotide in the multiple "
        "sequence alignment into a proportion value. It iteratively chooses a "
        "sliding window from each end of each sequence of the MSA and sums the "
        "proportion numbers in this window. Cropping continues until the sum of "
        "proportions is larger than the divergence threshold. Cropped "
        "nucleotides are converted to '-'.\n\n"
        "Window size (bp) is the sliding window used for this end-cropping "
        "process."
    ),
    "msa_clean_gap": (
        "Clean rows with gaps iteratively chooses a sliding window from each "
        "end of each sequence of the MSA and calculates the gap proportion in "
        "this window. Cropping continues until the gap proportion is smaller "
        "than the gap proportion threshold. Cropped nucleotides are converted "
        "to '-'.\n\n"
        "Window size (bp) is the sliding window used for this end-cropping "
        "process."
    ),
}

def msa_clean_method_description(method_id: str) -> str:
    return MSA_CLEAN_METHOD_DESCRIPTIONS.get(
        method_id,
        "",
    )


@dataclass
class MsaCleanOptions:
    """All MSA clean tool parameters (gappy columns + row cleaning)."""

    use_divergence: bool = True
    use_gap_row: bool = False
    use_gappy_columns: bool = True
    method_order: tuple[str, ...] = DEFAULT_MSA_CLEAN_ORDER
    gappy: GappyColumnOptions = field(default_factory=GappyColumnOptions)
    divergence: DivergenceRowOptions = field(default_factory=DivergenceRowOptions)
    gap_row: GapRowOptions = field(default_factory=GapRowOptions)
    right_click_action: str = MSA_CLEAN_RIGHT_CLICK_COMBINED
    show_process_window: bool = DEFAULT_SHOW_PROCESS_WINDOW
    close_process_window_when_finished: bool = DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED

    def __post_init__(self) -> None:
        valid = [method for method in self.method_order if method in DEFAULT_MSA_CLEAN_ORDER]
        valid.extend(method for method in DEFAULT_MSA_CLEAN_ORDER if method not in valid)
        self.method_order = tuple(valid)
        if not (self.use_divergence or self.use_gap_row or self.use_gappy_columns):
            self.use_gappy_columns = True
        if self.right_click_action not in MSA_CLEAN_RIGHT_CLICK_ACTIONS:
            self.right_click_action = MSA_CLEAN_RIGHT_CLICK_COMBINED


def load_msa_clean_options() -> MsaCleanOptions:
    from PySide6.QtCore import QSettings

    s = QSettings()
    s.beginGroup("msa_clean")
    if s.contains("use_divergence"):
        use_divergence = _to_bool(s.value("use_divergence"), False)
        use_gap_row = _to_bool(s.value("use_gap_row"), False)
        use_gappy_columns = _to_bool(s.value("use_gappy_columns"), True)
    elif s.contains("default_method"):
        # Migrate the former single-method toolbar setting into the equivalent
        # configurable pipeline.
        method = str(s.value("default_method", DEFAULT_MSA_CLEAN_METHOD))
        use_divergence = method in {"msa_clean_div", "msa_clean_gap", "msa_clean_row_col"}
        use_gap_row = method == "msa_clean_gap"
        use_gappy_columns = method in {"msa_clean_gappy", "msa_clean_row_col"}
    else:
        defaults = MsaCleanOptions()
        use_divergence = defaults.use_divergence
        use_gap_row = defaults.use_gap_row
        use_gappy_columns = defaults.use_gappy_columns
    stored_order = s.value("method_order", list(DEFAULT_MSA_CLEAN_ORDER))
    if isinstance(stored_order, str):
        stored_order = [stored_order]
    right_click_action = str(
        s.value("right_click_action", MSA_CLEAN_RIGHT_CLICK_COMBINED)
    )
    if right_click_action not in MSA_CLEAN_RIGHT_CLICK_ACTIONS:
        right_click_action = MSA_CLEAN_RIGHT_CLICK_COMBINED
    show_process_window = _to_bool(
        s.value("show_process_window"), DEFAULT_SHOW_PROCESS_WINDOW
    )
    close_process_window_when_finished = _to_bool(
        s.value("close_process_window_when_finished"),
        DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED,
    )
    opts = MsaCleanOptions(
        use_divergence=use_divergence,
        use_gap_row=use_gap_row,
        use_gappy_columns=use_gappy_columns,
        method_order=tuple(stored_order),
        gappy=GappyColumnOptions(
            gap_threshold=float(s.value("gappy_gap_threshold", DEFAULT_GAP_COL_GAP_THRESHOLD)),
            simi_check_gap_thr=float(
                s.value("gappy_simi_check_gap_thr", DEFAULT_GAP_COL_SIMI_GAP_THR)
            ),
            similarity_thr=float(
                s.value("gappy_similarity_thr", DEFAULT_GAP_COL_SIMILARITY_THR)
            ),
            min_nucleotide=int(s.value("gappy_min_nucleotide", DEFAULT_GAP_COL_MIN_NUCLEOTIDE)),
        ),
        divergence=DivergenceRowOptions(
            threshold=float(s.value("div_threshold", DEFAULT_DIV_THRESHOLD)),
            window_size=int(s.value("div_window_size", DEFAULT_DIV_WINDOW_SIZE)),
        ),
        gap_row=GapRowOptions(
            gap_threshold=float(s.value("gap_row_gap_threshold", DEFAULT_GAP_ROW_GAP_THRESHOLD)),
            window_size=int(s.value("gap_row_window_size", DEFAULT_GAP_ROW_WINDOW_SIZE)),
        ),
        right_click_action=right_click_action,
        show_process_window=show_process_window,
        close_process_window_when_finished=close_process_window_when_finished,
    )
    s.endGroup()
    return opts


def _load_settings_group(group: str) -> dict:
    from PySide6.QtCore import QSettings

    s = QSettings()
    s.beginGroup(group)
    keys = s.allKeys()
    out = {k: s.value(k) for k in keys}
    s.endGroup()
    return out


def _save_settings_group(group: str, values: dict) -> None:
    from PySide6.QtCore import QSettings

    s = QSettings()
    s.beginGroup(group)
    for key, value in values.items():
        s.setValue(key, value)
    s.endGroup()


def load_blast_options() -> BlastOptions:
    v = _load_settings_group("blast")
    program = str(v.get("program", BLAST_PROGRAM_BLASTN))
    if program not in BLAST_PROGRAMS:
        program = BLAST_PROGRAM_BLASTN
    return BlastOptions(
        program=program,
        e_value=normalize_blast_evalue(float(v.get("e_value", DEFAULT_BLAST_EVALUE))),
        num_threads=max(1, int(v.get("num_threads", _default_tool_threads()))),
        word_size=max(4, int(v.get("word_size", DEFAULT_BLAST_WORD_SIZE))),
        gap_open=max(0, int(v.get("gap_open", DEFAULT_BLAST_GAP_OPEN))),
        gap_extend=max(0, int(v.get("gap_extend", DEFAULT_BLAST_GAP_EXTEND))),
        penalty=min(0, int(v.get("penalty", DEFAULT_BLAST_PENALTY))),
        percent_identity=max(
            0.0,
            min(100.0, float(v.get("percent_identity", DEFAULT_BLAST_PERCENT_IDENTITY))),
        ),
        query_coverage_hsp=max(
            0.0,
            min(100.0, float(v.get("query_coverage_hsp", DEFAULT_BLAST_QUERY_COVERAGE))),
        ),
        max_msa_lines=int(v.get("max_msa_lines", DEFAULT_BLAST_MAX_MSA_LINES)),
        top_msa_lines=int(v.get("top_msa_lines", DEFAULT_BLAST_TOP_MSA_LINES)),
        write_text=bool(v.get("write_text", DEFAULT_BLAST_WRITE_TEXT)),
        write_gff=bool(v.get("write_gff", DEFAULT_BLAST_WRITE_GFF)),
        minimum_bit_score=float(
            v.get("minimum_bit_score", DEFAULT_BLAST_MINIMUM_BIT_SCORE)
        ),
        blast_database_name=(str(v["blast_database_name"]) if v.get("blast_database_name") else None),
        blastp_e_value=normalize_blast_evalue(
            float(v.get("blastp_e_value", DEFAULT_BLASTP_EVALUE))
        ),
        blastp_word_size=max(2, int(v.get("blastp_word_size", DEFAULT_BLASTP_WORD_SIZE))),
        blastp_gap_open=max(0, int(v.get("blastp_gap_open", DEFAULT_BLASTP_GAP_OPEN))),
        blastp_gap_extend=max(
            0, int(v.get("blastp_gap_extend", DEFAULT_BLASTP_GAP_EXTEND))
        ),
        blastp_query_coverage_hsp=max(
            0.0,
            min(
                100.0,
                float(v.get("blastp_query_coverage_hsp", DEFAULT_BLASTP_QUERY_COVERAGE)),
            ),
        ),
        blastp_database_name=(
            str(v["blastp_database_name"]) if v.get("blastp_database_name") else None
        ),
        blastp_write_text=bool(v.get("blastp_write_text", DEFAULT_BLASTP_WRITE_TEXT)),
        blastp_write_gff=bool(v.get("blastp_write_gff", DEFAULT_BLASTP_WRITE_GFF)),
        blastp_num_threads=max(
            1, int(v.get("blastp_num_threads", _default_tool_threads()))
        ),
        blastx_e_value=normalize_blast_evalue(
            float(v.get("blastx_e_value", DEFAULT_BLASTX_EVALUE))
        ),
        blastx_word_size=max(2, int(v.get("blastx_word_size", DEFAULT_BLASTX_WORD_SIZE))),
        blastx_gap_open=max(0, int(v.get("blastx_gap_open", DEFAULT_BLASTX_GAP_OPEN))),
        blastx_gap_extend=max(
            0, int(v.get("blastx_gap_extend", DEFAULT_BLASTX_GAP_EXTEND))
        ),
        blastx_query_coverage_hsp=max(
            0.0,
            min(
                100.0,
                float(v.get("blastx_query_coverage_hsp", DEFAULT_BLASTX_QUERY_COVERAGE)),
            ),
        ),
        blastx_database_name=(
            str(v["blastx_database_name"]) if v.get("blastx_database_name") else None
        ),
        blastx_write_text=bool(v.get("blastx_write_text", DEFAULT_BLASTX_WRITE_TEXT)),
        blastx_write_gff=bool(v.get("blastx_write_gff", DEFAULT_BLASTX_WRITE_GFF)),
        blastx_num_threads=max(
            1, int(v.get("blastx_num_threads", _default_tool_threads()))
        ),
        online_program=(
            str(v.get("online_program", DEFAULT_BLAST_ONLINE_PROGRAM))
            if str(v.get("online_program", DEFAULT_BLAST_ONLINE_PROGRAM))
            in BLAST_ONLINE_ALGORITHMS
            else DEFAULT_BLAST_ONLINE_PROGRAM
        ),
        online_database=str(
            v.get("online_database", DEFAULT_BLAST_ONLINE_DATABASE) or DEFAULT_BLAST_ONLINE_DATABASE
        ).strip()
        or DEFAULT_BLAST_ONLINE_DATABASE,
        online_e_value=normalize_blast_evalue(
            float(v.get("online_e_value", DEFAULT_BLAST_ONLINE_EVALUE))
        ),
        online_hitlist_size=max(
            1, min(500, int(v.get("online_hitlist_size", DEFAULT_BLAST_ONLINE_HITLIST)))
        ),
        online_write_text=bool(
            v.get("online_write_text", DEFAULT_BLAST_ONLINE_WRITE_TEXT)
        ),
        input_source=ALIGN_INPUT_SELECTED,
        input_paths=[],
        last_dialog_tab=(
            str(v.get("last_dialog_tab", program))
            if str(v.get("last_dialog_tab", program)) in BLAST_DIALOG_TABS
            else (
                program
                if program in BLAST_DIALOG_TABS
                else BLAST_PROGRAM_BLASTN
            )
        ),
        right_click_program=(
            str(v.get("right_click_program", BLAST_PROGRAM_BLASTN))
            if str(v.get("right_click_program", BLAST_PROGRAM_BLASTN)) in BLAST_PROGRAMS
            else BLAST_PROGRAM_BLASTN
        ),
        show_process_window=_to_bool(v.get("show_process_window"), DEFAULT_SHOW_PROCESS_WINDOW),
        close_process_window_when_finished=_to_bool(
            v.get("close_process_window_when_finished"), DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED
        ),
    )


def load_hmmer_options() -> HmmerOptions:
    v = _load_settings_group("hmmer")
    tab = str(v.get("last_dialog_tab", HMMER_DIALOG_TAB_PARAMETERS))
    if tab not in HMMER_DIALOG_TABS:
        tab = HMMER_DIALOG_TAB_PARAMETERS
    return HmmerOptions(
        database_name=str(v.get("database_name", "") or "").strip(),
        e_value=normalize_blast_evalue(
            float(v.get("e_value", DEFAULT_HMMER_EVALUE)),
            default=DEFAULT_HMMER_EVALUE,
        ),
        domain_e_value=normalize_blast_evalue(
            float(v.get("domain_e_value", DEFAULT_HMMER_DOMAIN_EVALUE)),
            default=DEFAULT_HMMER_DOMAIN_EVALUE,
        ),
        num_threads=max(1, int(v.get("num_threads", _default_tool_threads()))),
        max_hits=max(1, int(v.get("max_hits", DEFAULT_HMMER_MAX_HITS))),
        write_text=_to_bool(v.get("write_text"), DEFAULT_HMMER_WRITE_TEXT),
        minimum_bit_score=max(
            0.0, float(v.get("minimum_bit_score", DEFAULT_HMMER_MINIMUM_BIT_SCORE))
        ),
        max_orfs_per_sequence=max(
            0, int(v.get("max_orfs_per_sequence", DEFAULT_HMMER_MAX_ORFS))
        ),
        input_source=ALIGN_INPUT_SELECTED,
        input_paths=[],
        last_dialog_tab=tab,
        show_process_window=_to_bool(v.get("show_process_window"), DEFAULT_SHOW_PROCESS_WINDOW),
        close_process_window_when_finished=_to_bool(
            v.get("close_process_window_when_finished"), DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED
        ),
        right_click_action=_normalize_hmmer_right_click(
            v.get("right_click_action", HMMER_RIGHT_CLICK_RUN)
        ),
    )


def _normalize_hmmer_right_click(raw: object) -> str:
    value = str(raw or HMMER_RIGHT_CLICK_RUN).strip()
    if value in HMMER_RIGHT_CLICK_ACTIONS:
        return value
    return HMMER_RIGHT_CLICK_RUN


def save_hmmer_options(opts: HmmerOptions) -> None:
    _save_settings_group(
        "hmmer",
        {
            "database_name": opts.database_name or "",
            "e_value": opts.e_value,
            "domain_e_value": opts.domain_e_value,
            "num_threads": opts.num_threads,
            "max_hits": opts.max_hits,
            "write_text": opts.write_text,
            "minimum_bit_score": opts.minimum_bit_score,
            "max_orfs_per_sequence": opts.max_orfs_per_sequence,
            "last_dialog_tab": (
                opts.last_dialog_tab
                if opts.last_dialog_tab in HMMER_DIALOG_TABS
                else HMMER_DIALOG_TAB_PARAMETERS
            ),
            "show_process_window": opts.show_process_window,
            "close_process_window_when_finished": (
                opts.close_process_window_when_finished
            ),
            "right_click_action": (
                opts.right_click_action
                if opts.right_click_action in HMMER_RIGHT_CLICK_ACTIONS
                else HMMER_RIGHT_CLICK_RUN
            ),
        },
    )


def save_blast_options(opts: BlastOptions) -> None:
    _save_settings_group(
        "blast",
        {
            "program": opts.program,
            "e_value": opts.e_value,
            "num_threads": opts.num_threads,
            "word_size": opts.word_size,
            "gap_open": opts.gap_open,
            "gap_extend": opts.gap_extend,
            "penalty": opts.penalty,
            "percent_identity": opts.percent_identity,
            "query_coverage_hsp": opts.query_coverage_hsp,
            "max_msa_lines": opts.max_msa_lines,
            "top_msa_lines": opts.top_msa_lines,
            "write_text": opts.write_text,
            "write_gff": opts.write_gff,
            "minimum_bit_score": opts.minimum_bit_score,
            "blast_database_name": opts.blast_database_name or "",
            "blastp_e_value": opts.blastp_e_value,
            "blastp_word_size": opts.blastp_word_size,
            "blastp_gap_open": opts.blastp_gap_open,
            "blastp_gap_extend": opts.blastp_gap_extend,
            "blastp_query_coverage_hsp": opts.blastp_query_coverage_hsp,
            "blastp_database_name": opts.blastp_database_name or "",
            "blastp_write_text": opts.blastp_write_text,
            "blastp_write_gff": opts.blastp_write_gff,
            "blastp_num_threads": opts.blastp_num_threads,
            "blastx_e_value": opts.blastx_e_value,
            "blastx_word_size": opts.blastx_word_size,
            "blastx_gap_open": opts.blastx_gap_open,
            "blastx_gap_extend": opts.blastx_gap_extend,
            "blastx_query_coverage_hsp": opts.blastx_query_coverage_hsp,
            "blastx_database_name": opts.blastx_database_name or "",
            "blastx_write_text": opts.blastx_write_text,
            "blastx_write_gff": opts.blastx_write_gff,
            "blastx_num_threads": opts.blastx_num_threads,
            "online_program": opts.online_program,
            "online_database": opts.online_database,
            "online_e_value": opts.online_e_value,
            "online_hitlist_size": opts.online_hitlist_size,
            "online_write_text": opts.online_write_text,
            "last_dialog_tab": (
                opts.last_dialog_tab
                if opts.last_dialog_tab in BLAST_DIALOG_TABS
                else opts.program
            ),
            "right_click_program": opts.right_click_program,
            "show_process_window": opts.show_process_window,
            "close_process_window_when_finished": (
                opts.close_process_window_when_finished
            ),
        },
    )


def load_extend_options() -> ExtendOptions:
    v = _load_settings_group("extend")
    pan_name = v.get("pan_genome_name")
    genome_path = v.get("genome_path")
    right_click_action = str(v.get("right_click_action", EXTEND_RIGHT_CLICK_RUN))
    if right_click_action not in EXTEND_RIGHT_CLICK_ACTIONS:
        right_click_action = EXTEND_RIGHT_CLICK_RUN
    return ExtendOptions(
        left=int(v.get("left", DEFAULT_EXTEND_LEFT)),
        right=int(v.get("right", DEFAULT_EXTEND_RIGHT)),
        pan_genome_name=str(pan_name).strip() if pan_name else None,
        genome_path=str(genome_path).strip() if genome_path else None,
        right_click_action=right_click_action,
        show_process_window=_to_bool(
            v.get("show_process_window"), DEFAULT_SHOW_PROCESS_WINDOW
        ),
        close_process_window_when_finished=_to_bool(
            v.get("close_process_window_when_finished"),
            DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED,
        ),
    )


def save_extend_options(opts: ExtendOptions) -> None:
    _save_settings_group(
        "extend",
        {
            "left": opts.left,
            "right": opts.right,
            "pan_genome_name": opts.pan_genome_name or "",
            "genome_path": opts.genome_path or "",
            "right_click_action": opts.right_click_action,
            "show_process_window": opts.show_process_window,
            "close_process_window_when_finished": (
                opts.close_process_window_when_finished
            ),
        },
    )


def load_tsd_search_options() -> TsdSearchOptions:
    v = _load_settings_group("tsd_search")
    return TsdSearchOptions(
        left_boundary=int(v.get("left_boundary", DEFAULT_TSD_LEFT_BOUNDARY)),
        right_boundary=int(v.get("right_boundary", DEFAULT_TSD_RIGHT_BOUNDARY)),
        nt_radius=max(1, int(v.get("nt_radius", DEFAULT_TSD_NT_RADIUS))),
        max_sim_threshold=float(
            v.get("max_sim_threshold", DEFAULT_TSD_MAX_SIM_THRESHOLD)
        ),
        nt_threshold=float(v.get("nt_threshold", DEFAULT_TSD_NT_THRESHOLD)),
        cluster_eps=float(v.get("cluster_eps", DEFAULT_TSD_CLUSTER_EPS)),
    )


def save_tsd_search_options(opts: TsdSearchOptions) -> None:
    _save_settings_group(
        "tsd_search",
        {
            "left_boundary": opts.left_boundary,
            "right_boundary": opts.right_boundary,
            "nt_radius": opts.nt_radius,
            "max_sim_threshold": opts.max_sim_threshold,
            "nt_threshold": opts.nt_threshold,
            "cluster_eps": opts.cluster_eps,
        },
    )


def load_teaid_options() -> TeaAidOptions:
    v = _load_settings_group("teaid")
    word_size = int(v.get("self_blast_word_size", DEFAULT_TEAID_SELF_BLAST_WORD_SIZE))
    word_size = max(
        MIN_TEAID_SELF_BLAST_WORD_SIZE,
        min(MAX_TEAID_SELF_BLAST_WORD_SIZE, word_size),
    )
    high_pct, mid_pct = clamp_teaid_hit_len_pct(
        float(v.get("hit_len_high_pct", DEFAULT_TEAID_HIT_LEN_HIGH_PCT)),
        float(v.get("hit_len_mid_pct", DEFAULT_TEAID_HIT_LEN_MID_PCT)),
    )
    right_click_action = str(v.get("right_click_action", TEAID_RIGHT_CLICK_RUN))
    if right_click_action not in TEAID_RIGHT_CLICK_ACTIONS:
        right_click_action = TEAID_RIGHT_CLICK_RUN
    genome_path = v.get("genome_path")
    pan_raw = v.get("pan_genomes", [])
    pan_genomes: list[tuple[str, str]] = []
    if isinstance(pan_raw, list):
        for item in pan_raw:
            if isinstance(item, (list, tuple)) and len(item) >= 2:
                name = str(item[0]).strip()
                path = str(item[1]).strip()
            elif isinstance(item, dict):
                name = str(item.get("name", "")).strip()
                path = str(item.get("path", "")).strip()
            else:
                continue
            if name and path:
                pan_genomes.append((name, path))
    return TeaAidOptions(
        e_value=normalize_blast_evalue(
            float(v.get("e_value", DEFAULT_TEAID_EVALUE)),
            default=DEFAULT_TEAID_EVALUE,
        ),
        self_blast_word_size=word_size,
        num_threads=max(1, int(v.get("num_threads", _default_tool_threads()))),
        x_axis_label=(
            str(v.get("x_axis_label", DEFAULT_TEAID_X_AXIS_LABEL)).strip()
            or DEFAULT_TEAID_X_AXIS_LABEL
        ),
        hit_len_high_pct=high_pct,
        hit_len_mid_pct=mid_pct,
        hit_color_high=(
            str(v.get("hit_color_high", DEFAULT_TEAID_HIT_COLOR_HIGH)).strip()
            or DEFAULT_TEAID_HIT_COLOR_HIGH
        ),
        hit_color_mid=(
            str(v.get("hit_color_mid", DEFAULT_TEAID_HIT_COLOR_MID)).strip()
            or DEFAULT_TEAID_HIT_COLOR_MID
        ),
        hit_line_transparency=_load_teaid_line_transparency(v),
        right_click_action=right_click_action,
        genome_path=str(genome_path).strip() if genome_path else None,
        use_pan_genome=_to_bool(v.get("use_pan_genome"), False),
        pan_genomes=pan_genomes,
        show_process_window=_to_bool(v.get("show_process_window"), DEFAULT_SHOW_PROCESS_WINDOW),
        close_process_window_when_finished=_to_bool(
            v.get("close_process_window_when_finished"), DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED
        ),
    )


def _load_teaid_line_transparency(values: dict) -> float:
    """Load line transparency; migrate older per-band opacity settings if needed."""
    if "hit_line_transparency" in values:
        return clamp_teaid_hit_line_transparency(
            float(values.get("hit_line_transparency", DEFAULT_TEAID_HIT_LINE_TRANSPARENCY))
        )
    # Older builds stored opacity (1 = opaque). Convert the stronger of the two.
    legacy_high = values.get("hit_opacity_high")
    legacy_mid = values.get("hit_opacity_mid")
    if legacy_high is None and legacy_mid is None:
        return DEFAULT_TEAID_HIT_LINE_TRANSPARENCY
    opacities = []
    for raw in (legacy_high, legacy_mid):
        if raw is None:
            continue
        try:
            opacities.append(max(0.0, min(1.0, float(raw))))
        except (TypeError, ValueError):
            continue
    if not opacities:
        return DEFAULT_TEAID_HIT_LINE_TRANSPARENCY
    return clamp_teaid_hit_line_transparency(1.0 - min(opacities))


def save_teaid_options(opts: TeaAidOptions) -> None:
    high_pct, mid_pct = clamp_teaid_hit_len_pct(
        opts.hit_len_high_pct, opts.hit_len_mid_pct
    )
    _save_settings_group(
        "teaid",
        {
            "e_value": opts.e_value,
            "self_blast_word_size": opts.self_blast_word_size,
            "num_threads": max(1, int(opts.num_threads)),
            "x_axis_label": opts.x_axis_label,
            "hit_len_high_pct": high_pct,
            "hit_len_mid_pct": mid_pct,
            "hit_color_high": opts.hit_color_high,
            "hit_color_mid": opts.hit_color_mid,
            "hit_line_transparency": clamp_teaid_hit_line_transparency(
                opts.hit_line_transparency
            ),
            "right_click_action": (
                opts.right_click_action
                if opts.right_click_action in TEAID_RIGHT_CLICK_ACTIONS
                else TEAID_RIGHT_CLICK_RUN
            ),
            "genome_path": opts.genome_path or "",
            "use_pan_genome": bool(opts.use_pan_genome),
            "pan_genomes": [
                {"name": name, "path": path}
                for name, path in (opts.pan_genomes or [])
                if str(name).strip() and str(path).strip()
            ],
            "show_process_window": opts.show_process_window,
            "close_process_window_when_finished": (
                opts.close_process_window_when_finished
            ),
        },
    )


def load_align_options() -> AlignOptions:
    v = _load_settings_group("align")
    strategy = str(v.get("strategy", MAFFT_STRATEGY_FFT_NS_1))
    if strategy not in MAFFT_STRATEGIES:
        strategy = MAFFT_STRATEGY_FFT_NS_1
    sequence_type = str(v.get("sequence_type", MAFFT_SEQUENCE_AUTO))
    if sequence_type not in MAFFT_SEQUENCE_TYPES:
        sequence_type = MAFFT_SEQUENCE_AUTO
    combined_aligner = str(
        v.get("combined_aligner", v.get("default_aligner", ALIGNER_MAFFT))
    )
    if combined_aligner not in ALIGNERS:
        combined_aligner = ALIGNER_MAFFT
    muscle_mode = str(v.get("muscle_mode", MUSCLE_MODE_ALIGN))
    if muscle_mode not in MUSCLE_MODES:
        muscle_mode = MUSCLE_MODE_ALIGN
    muscle_perm = str(v.get("muscle_perm", MUSCLE_PERM_NONE))
    if muscle_perm not in MUSCLE_PERMS:
        muscle_perm = MUSCLE_PERM_NONE
    muscle_input_source = str(v.get("muscle_input_source", ALIGN_INPUT_SELECTED) or "")
    if muscle_input_source not in ALIGN_INPUT_SOURCES:
        muscle_input_source = ALIGN_INPUT_SELECTED
    raw_muscle_paths = v.get("muscle_input_paths", []) or []
    if isinstance(raw_muscle_paths, str):
        muscle_input_paths = [p.strip() for p in raw_muscle_paths.split("\n") if p.strip()]
    else:
        muscle_input_paths = [str(p).strip() for p in raw_muscle_paths if str(p).strip()]
    mafft_input_source = str(v.get("mafft_input_source", ALIGN_INPUT_SELECTED) or "")
    if mafft_input_source not in ALIGN_INPUT_SOURCES:
        mafft_input_source = ALIGN_INPUT_SELECTED
    raw_mafft_paths = v.get("mafft_input_paths", []) or []
    if isinstance(raw_mafft_paths, str):
        mafft_input_paths = [p.strip() for p in raw_mafft_paths.split("\n") if p.strip()]
    else:
        mafft_input_paths = [str(p).strip() for p in raw_mafft_paths if str(p).strip()]
    right_click_aligner = str(v.get("right_click_aligner", ALIGNER_MAFFT))
    if right_click_aligner not in RIGHT_CLICK_ALIGNERS:
        right_click_aligner = ALIGNER_MAFFT

    def _flag(key: str) -> bool:
        return str(v.get(key, "false")).lower() in {"1", "true", "yes"}

    return AlignOptions(
        combined_aligner=combined_aligner,
        num_threads=max(1, int(v.get("num_threads", _default_tool_threads()))),
        strategy=strategy,
        sequence_type=sequence_type,
        retree=max(1, int(v.get("retree", 1))),
        max_iterate=max(0, int(v.get("max_iterate", 0))),
        gap_open_penalty=float(v.get("gap_open_penalty", 1.53)),
        offset=float(v.get("offset", 0.0)),
        reorder=_flag("reorder"),
        adjust_direction=_flag("adjust_direction"),
        muscle_mode=muscle_mode,
        muscle_num_threads=max(
            1, int(v.get("muscle_num_threads", _default_tool_threads()))
        ),
        muscle_perm=muscle_perm,
        muscle_perturb=max(0, int(v.get("muscle_perturb", 0))),
        muscle_input_source=muscle_input_source,
        muscle_input_paths=muscle_input_paths,
        mafft_input_source=mafft_input_source,
        mafft_input_paths=mafft_input_paths,
        mafft_combined=_flag("mafft_combined"),
        muscle_combined=_flag("muscle_combined"),
        last_dialog_tab=(
            str(v.get("last_dialog_tab", ALIGNER_MAFFT))
            if str(v.get("last_dialog_tab", ALIGNER_MAFFT)) in ALIGN_DIALOG_TABS
            else ALIGNER_MAFFT
        ),
        right_click_aligner=right_click_aligner,
        show_process_window=_to_bool(
            v.get("show_process_window"), DEFAULT_SHOW_PROCESS_WINDOW
        ),
        close_process_window_when_finished=_to_bool(
            v.get("close_process_window_when_finished"),
            DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED,
        ),
    )


def save_align_options(opts: AlignOptions) -> None:
    _save_settings_group(
        "align",
        {
            "combined_aligner": opts.combined_aligner,
            "num_threads": opts.num_threads,
            "strategy": opts.strategy,
            "sequence_type": opts.sequence_type,
            "retree": opts.retree,
            "max_iterate": opts.max_iterate,
            "gap_open_penalty": opts.gap_open_penalty,
            "offset": opts.offset,
            "reorder": opts.reorder,
            "adjust_direction": opts.adjust_direction,
            "muscle_mode": opts.muscle_mode,
            "muscle_num_threads": opts.muscle_num_threads,
            "muscle_perm": opts.muscle_perm,
            "muscle_perturb": opts.muscle_perturb,
            "muscle_input_source": opts.muscle_input_source,
            "muscle_input_paths": list(opts.muscle_input_paths or []),
            "mafft_input_source": opts.mafft_input_source,
            "mafft_input_paths": list(opts.mafft_input_paths or []),
            "mafft_combined": opts.mafft_combined,
            "muscle_combined": opts.muscle_combined,
            "last_dialog_tab": (
                opts.last_dialog_tab
                if opts.last_dialog_tab in ALIGN_DIALOG_TABS
                else ALIGNER_MAFFT
            ),
            "right_click_aligner": opts.right_click_aligner,
            "show_process_window": opts.show_process_window,
            "close_process_window_when_finished": (
                opts.close_process_window_when_finished
            ),
        },
    )


def _parse_chromosomes_only_setting(raw) -> list[str] | None:
    if raw is None or raw == "":
        return None
    if isinstance(raw, (list, tuple)):
        names = [str(x).strip() for x in raw if str(x).strip()]
        return names or None
    text = str(raw).strip()
    if not text:
        return None
    names = [part.strip() for part in text.replace(",", "\n").splitlines() if part.strip()]
    return names or None


def _parse_name_list_setting(raw) -> list[str] | None:
    if raw is None or raw == "":
        return None
    if raw == "__NONE__":
        return []
    if isinstance(raw, (list, tuple)):
        values = [str(x).strip() for x in raw if str(x).strip()]
        return values or None
    text = str(raw).strip()
    if not text:
        return None
    values = [part.strip() for part in text.replace(",", "\n").splitlines() if part.strip()]
    return values or None


def _serialize_annotation_feature_types_by_file(
    by_file: dict[str, list[str]] | None,
) -> str:
    if not by_file:
        return ""
    payload: dict[str, list[str]] = {}
    for raw_path, types in by_file.items():
        path = str(raw_path).strip()
        if not path:
            continue
        cleaned = [str(t).strip() for t in types if str(t).strip()]
        if cleaned:
            payload[path] = cleaned
    if not payload:
        return ""
    return json.dumps(payload, sort_keys=True)


def _parse_annotation_feature_types_by_file(raw: object) -> dict[str, list[str]] | None:
    if raw is None:
        return None
    if isinstance(raw, dict):
        out: dict[str, list[str]] = {}
        for key, value in raw.items():
            path = str(key).strip()
            if not path:
                continue
            if isinstance(value, str):
                types = [p.strip() for p in value.replace(",", "\n").splitlines() if p.strip()]
            elif isinstance(value, (list, tuple)):
                types = [str(t).strip() for t in value if str(t).strip()]
            else:
                continue
            if types:
                out[path] = types
        return out or None
    if isinstance(raw, str) and raw.strip():
        try:
            parsed = json.loads(raw)
        except json.JSONDecodeError:
            return None
        return _parse_annotation_feature_types_by_file(parsed)
    return None


def load_te_genome_map_options() -> TeGenomeMapOptions:
    v = _load_settings_group("te_genome_map")
    mode = str(v.get("density_mode", DEFAULT_TE_GENOME_MAP_DENSITY))
    if mode not in (TE_GENOME_MAP_DENSITY_HIT_COUNT, TE_GENOME_MAP_DENSITY_TE_LENGTH):
        mode = DEFAULT_TE_GENOME_MAP_DENSITY
    color_by = str(v.get("color_by", DEFAULT_TE_GENOME_MAP_COLOR_BY))
    if color_by not in (
        TE_GENOME_MAP_COLOR_IDENTITY,
        TE_GENOME_MAP_COLOR_ORDER,
        TE_GENOME_MAP_COLOR_SUPERFAMILY,
        TE_GENOME_MAP_COLOR_FAMILY,
    ):
        color_by = DEFAULT_TE_GENOME_MAP_COLOR_BY
    split_strand = str(
        v.get("split_strand_tracks", str(DEFAULT_TE_GENOME_MAP_SPLIT_STRAND).lower())
    ).lower() in {"1", "true", "yes"}
    show_window_ref = str(
        v.get("show_window_ref_line", str(DEFAULT_TE_GENOME_MAP_WINDOW_REF_LINE).lower())
    ).lower() in {"1", "true", "yes"}
    collapse_overlaps = str(
        v.get("collapse_overlaps", str(DEFAULT_TE_GENOME_MAP_COLLAPSE_OVERLAPS).lower())
    ).lower() in {"1", "true", "yes"}
    show_te_diff = str(
        v.get("show_te_diff_track", str(DEFAULT_TE_GENOME_MAP_SHOW_DIFF_TRACK).lower())
    ).lower() in {"1", "true", "yes"}
    show_meth_diff = str(
        v.get(
            "show_meth_diff_track",
            str(DEFAULT_TE_GENOME_MAP_SHOW_METH_DIFF_TRACK).lower(),
        )
    ).lower() in {"1", "true", "yes"}
    ann_mode = str(v.get("annotation_density_mode", DEFAULT_ANNOTATION_DENSITY))
    if ann_mode not in (TE_GENOME_MAP_DENSITY_HIT_COUNT, TE_GENOME_MAP_DENSITY_TE_LENGTH):
        ann_mode = DEFAULT_ANNOTATION_DENSITY
    ann_types_raw = v.get("annotation_feature_types")
    ann_types: list[str] | None = None
    if isinstance(ann_types_raw, str) and ann_types_raw.strip():
        ann_types = [p.strip() for p in ann_types_raw.replace(",", "\n").splitlines() if p.strip()]
    ann_by_file = _parse_annotation_feature_types_by_file(
        v.get("annotation_feature_types_by_file"),
    )
    chrom_sort = str(v.get("chromosome_sort", DEFAULT_CHROMOSOME_SORT))
    if chrom_sort not in (
        CHROMOSOME_SORT_NATURAL,
        CHROMOSOME_SORT_SIZE_DESC,
        CHROMOSOME_SORT_SIZE_ASC,
    ):
        chrom_sort = DEFAULT_CHROMOSOME_SORT
    tab = str(v.get("last_dialog_tab", TE_GENOME_MAP_DIALOG_TAB_TE))
    if tab == "parameters":
        tab = TE_GENOME_MAP_DIALOG_TAB_TE
    if tab not in TE_GENOME_MAP_DIALOG_TABS:
        tab = TE_GENOME_MAP_DIALOG_TAB_TE
    if "plot_all_tracks_together" in v:
        plot_all = _to_bool(v.get("plot_all_tracks_together"), False)
    elif any(
        key in v
        for key in (
            "include_te_tracks",
            "include_annotation_tracks",
            "include_methylation_tracks",
        )
    ):
        plot_all = (
            _to_bool(v.get("include_te_tracks"), True)
            and _to_bool(v.get("include_annotation_tracks"), True)
            and _to_bool(v.get("include_methylation_tracks"), True)
        )
    else:
        plot_all = False
    return TeGenomeMapOptions(
        bin_size_bp=max(5_000, int(v.get("bin_size_bp", DEFAULT_TE_GENOME_MAP_BIN_BP))),
        density_mode=mode,
        color_by=color_by,
        methylation_coverage_min=max(
            1, int(v.get("methylation_coverage_min", DEFAULT_TE_GENOME_MAP_METHYLATION_COVERAGE_MIN))
        ),
        split_strand_tracks=split_strand,
        collapse_overlaps=collapse_overlaps,
        show_te_diff_track=show_te_diff,
        show_meth_diff_track=show_meth_diff,
        show_window_ref_line=show_window_ref,
        annotation_density_mode=ann_mode,
        annotation_feature_types=ann_types,
        annotation_feature_types_by_file=ann_by_file,
        chromosomes_only=_parse_chromosomes_only_setting(v.get("chromosomes_only")),
        superfamilies_only=_parse_name_list_setting(v.get("superfamilies_only")),
        chromosome_sort=chrom_sort,
        plot_all_tracks_together=plot_all,
        genome_path=str(v.get("genome_path", "") or "").strip(),
        show_process_window=_to_bool(
            v.get("show_process_window"), DEFAULT_SHOW_PROCESS_WINDOW
        ),
        close_process_window_when_finished=_to_bool(
            v.get("close_process_window_when_finished"),
            DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED,
        ),
        last_dialog_tab=tab,
    )


def save_te_genome_map_options(opts: TeGenomeMapOptions) -> None:
    chrom_value = ""
    if opts.chromosomes_only:
        chrom_value = "\n".join(opts.chromosomes_only)
    sf_value = ""
    if opts.superfamilies_only == []:
        sf_value = "__NONE__"
    elif opts.superfamilies_only:
        sf_value = "\n".join(opts.superfamilies_only)
    tab = opts.last_dialog_tab
    if tab == "parameters":
        tab = TE_GENOME_MAP_DIALOG_TAB_TE
    if tab not in TE_GENOME_MAP_DIALOG_TABS:
        tab = TE_GENOME_MAP_DIALOG_TAB_TE
    _save_settings_group(
        "te_genome_map",
        {
            "bin_size_bp": opts.bin_size_bp,
            "density_mode": opts.density_mode,
            "color_by": opts.color_by,
            "methylation_coverage_min": opts.methylation_coverage_min,
            "split_strand_tracks": opts.split_strand_tracks,
            "collapse_overlaps": opts.collapse_overlaps,
            "show_te_diff_track": opts.show_te_diff_track,
            "show_meth_diff_track": opts.show_meth_diff_track,
            "show_window_ref_line": opts.show_window_ref_line,
            "annotation_density_mode": opts.annotation_density_mode,
            "annotation_feature_types": (
                ",".join(opts.annotation_feature_types) if opts.annotation_feature_types else ""
            ),
            "annotation_feature_types_by_file": _serialize_annotation_feature_types_by_file(
                opts.annotation_feature_types_by_file,
            ),
            "chromosomes_only": chrom_value,
            "superfamilies_only": sf_value,
            "chromosome_sort": opts.chromosome_sort,
            "plot_all_tracks_together": opts.plot_all_tracks_together,
            "genome_path": opts.genome_path,
            "show_process_window": opts.show_process_window,
            "close_process_window_when_finished": opts.close_process_window_when_finished,
            "last_dialog_tab": tab,
        },
    )


def load_te_divergence_plot_options() -> TeDivergencePlotOptions:
    v = _load_settings_group("te_divergence_plot")
    return TeDivergencePlotOptions(
        bin_width=max(0.001, float(v.get("bin_width", DEFAULT_TE_DIVERGENCE_BIN_WIDTH))),
    )


def save_te_divergence_plot_options(opts: TeDivergencePlotOptions) -> None:
    _save_settings_group(
        "te_divergence_plot",
        {"bin_width": opts.bin_width},
    )


DEFAULT_MSA_PLOT_DPI = 100
# Per-sequence row height in centimeters (matplotlib figsize uses inches).
DEFAULT_MSA_PLOT_SEQUENCE_HEIGHT_CM = 0.5
# White distance between sequence rows (cm).
DEFAULT_MSA_PLOT_SEQUENCE_DISTANCE_CM = 0.01
# Bump when built-in nucleotide/amino defaults change so old QSettings palettes refresh.
MSA_PLOT_COLOR_SCHEME_VERSION = 3


@dataclass
class MsaPlotJobOptions:
    """MSA overview plot (Visualization → MSA plot)."""

    sequence_type: str = "auto"
    output_format: str = "png"
    dpi: int = DEFAULT_MSA_PLOT_DPI
    sequence_height_cm: float = DEFAULT_MSA_PLOT_SEQUENCE_HEIGHT_CM
    sequence_distance_cm: float = DEFAULT_MSA_PLOT_SEQUENCE_DISTANCE_CM
    nucleotide_colors: dict[str, str] = field(default_factory=dict)
    amino_acid_colors: dict[str, str] = field(default_factory=dict)


def _parse_color_map_setting(raw: object) -> dict[str, str]:
    if not raw:
        return {}
    text = str(raw).strip()
    if not text:
        return {}
    try:
        import json

        decoded = json.loads(text)
    except (json.JSONDecodeError, TypeError):
        return {}
    if not isinstance(decoded, dict):
        return {}
    colors: dict[str, str] = {}
    for key, value in decoded.items():
        color = str(value).strip()
        if key and color.startswith("#") and len(color) >= 4:
            colors[str(key).upper()] = color
    return colors


def load_msa_plot_options() -> MsaPlotJobOptions:
    from repeatviewer.backends.msa_plot import (
        DEFAULT_AMINO_ACID_COLORS,
        DEFAULT_NUCLEOTIDE_COLORS,
        MSA_PLOT_FORMAT_PDF,
        MSA_PLOT_FORMAT_PNG,
        MSA_PLOT_SEQ_AUTO,
        MSA_PLOT_SEQ_NUCLEOTIDE,
        MSA_PLOT_SEQ_PROTEIN,
    )

    v = _load_settings_group("msa_plot")
    seq_type = str(v.get("sequence_type", MSA_PLOT_SEQ_AUTO)).lower()
    if seq_type not in {
        MSA_PLOT_SEQ_AUTO,
        MSA_PLOT_SEQ_NUCLEOTIDE,
        MSA_PLOT_SEQ_PROTEIN,
    }:
        seq_type = MSA_PLOT_SEQ_AUTO
    fmt = str(v.get("output_format", MSA_PLOT_FORMAT_PNG)).lower()
    if fmt not in {MSA_PLOT_FORMAT_PNG, MSA_PLOT_FORMAT_PDF}:
        fmt = MSA_PLOT_FORMAT_PNG

    try:
        scheme_version = int(v.get("color_scheme_version", 0))
    except (TypeError, ValueError):
        scheme_version = 0

    if scheme_version < MSA_PLOT_COLOR_SCHEME_VERSION:
        # Refresh built-in palettes after Default / Color-blind updates.
        nuc = dict(DEFAULT_NUCLEOTIDE_COLORS)
        aa = dict(DEFAULT_AMINO_ACID_COLORS)
    else:
        nuc = dict(DEFAULT_NUCLEOTIDE_COLORS)
        nuc.update(_parse_color_map_setting(v.get("nucleotide_colors")))
        aa = dict(DEFAULT_AMINO_ACID_COLORS)
        aa.update(_parse_color_map_setting(v.get("amino_acid_colors")))
    # DNA/RNA: keep T and U identical in the UI and plot.
    if "U" in nuc:
        nuc["T"] = nuc["U"]

    # Prefer cm; migrate legacy inch values (sequence_height without _cm).
    if "sequence_height_cm" in v:
        height_cm = float(v.get("sequence_height_cm", DEFAULT_MSA_PLOT_SEQUENCE_HEIGHT_CM))
    elif "sequence_height" in v:
        height_cm = float(v.get("sequence_height")) * 2.54
    else:
        height_cm = DEFAULT_MSA_PLOT_SEQUENCE_HEIGHT_CM

    return MsaPlotJobOptions(
        sequence_type=seq_type,
        output_format=fmt,
        dpi=max(50, min(600, int(v.get("dpi", DEFAULT_MSA_PLOT_DPI)))),
        sequence_height_cm=max(0.1, min(5.0, height_cm)),
        sequence_distance_cm=max(
            0.0,
            min(
                2.0,
                float(
                    v.get(
                        "sequence_distance_cm",
                        v.get("sequence_gap_cm", DEFAULT_MSA_PLOT_SEQUENCE_DISTANCE_CM),
                    )
                ),
            ),
        ),
        nucleotide_colors=nuc,
        amino_acid_colors=aa,
    )


def save_msa_plot_options(opts: MsaPlotJobOptions) -> None:
    import json

    nuc = dict(opts.nucleotide_colors)
    if "U" in nuc:
        nuc["T"] = nuc["U"]
    _save_settings_group(
        "msa_plot",
        {
            "sequence_type": opts.sequence_type,
            "output_format": opts.output_format,
            "dpi": int(opts.dpi),
            "sequence_height_cm": float(opts.sequence_height_cm),
            "sequence_distance_cm": float(opts.sequence_distance_cm),
            "color_scheme_version": MSA_PLOT_COLOR_SCHEME_VERSION,
            "nucleotide_colors": json.dumps(nuc, sort_keys=True),
            "amino_acid_colors": json.dumps(opts.amino_acid_colors, sort_keys=True),
        },
    )


REPEAT_QUANTIFY_ENGINES = ("blastn", "repeatmasker", "both")
REPEAT_QUANTIFY_SENSITIVITIES = ("sensitive", "default", "quick")


@dataclass
class RepeatQuantifyOptions:
    """Repeat Quantify engine choice and per-engine parameters."""

    engine: str = "both"
    blast_word_size: int = 11
    # The dialog's E-value box has four decimals (min 0.0001), so 1e-5 was
    # never actually selectable there.
    blast_evalue: float = DEFAULT_BLAST_EVALUE
    blast_min_identity: float = 0.0
    blast_min_coverage: float = 0.0
    rm_div: int = 40
    rm_sensitivity: str = "sensitive"
    # Save a genome's BLAST index as a named database (Add/Manage databases)
    # so later runs and other tools reuse it. Not a quantify_library keyword.
    blast_save_index: bool = True

    def params(self) -> dict:
        """Keyword arguments for ``quantify_library`` (everything but engine)."""
        return {
            "blast_word_size": self.blast_word_size,
            "blast_evalue": self.blast_evalue,
            "blast_min_identity": self.blast_min_identity,
            "blast_min_coverage": self.blast_min_coverage,
            "rm_div": self.rm_div,
            "rm_sensitivity": self.rm_sensitivity,
        }


def _clamped(value, cast, default, low, high):
    try:
        out = cast(value)
    except (TypeError, ValueError):
        return default
    return min(max(out, low), high)


def load_repeat_quantify_options() -> RepeatQuantifyOptions:
    v = _load_settings_group("repeat_quantify")
    d = RepeatQuantifyOptions()
    engine = str(v.get("engine", d.engine))
    sensitivity = str(v.get("rm_sensitivity", d.rm_sensitivity))
    return RepeatQuantifyOptions(
        engine=engine if engine in REPEAT_QUANTIFY_ENGINES else d.engine,
        blast_word_size=_clamped(
            v.get("blast_word_size", d.blast_word_size), int, d.blast_word_size, 4, 24
        ),
        blast_evalue=_clamped(
            v.get("blast_evalue", d.blast_evalue), float, d.blast_evalue,
            MIN_BLAST_EVALUE, 1.0,
        ),
        blast_min_identity=_clamped(
            v.get("blast_min_identity", d.blast_min_identity), float,
            d.blast_min_identity, 0.0, 100.0,
        ),
        blast_min_coverage=_clamped(
            v.get("blast_min_coverage", d.blast_min_coverage), float,
            d.blast_min_coverage, 0.0, 100.0,
        ),
        rm_div=_clamped(v.get("rm_div", d.rm_div), int, d.rm_div, 1, 100),
        rm_sensitivity=(
            sensitivity if sensitivity in REPEAT_QUANTIFY_SENSITIVITIES
            else d.rm_sensitivity
        ),
        blast_save_index=_setting_bool(
            v.get("blast_save_index", d.blast_save_index), d.blast_save_index
        ),
    )


def _setting_bool(value, default: bool) -> bool:
    """QSettings hands booleans back as "true"/"false" strings on some platforms."""
    if isinstance(value, bool):
        return value
    text = str(value).strip().lower()
    if text in {"true", "1", "yes"}:
        return True
    if text in {"false", "0", "no"}:
        return False
    return default


def save_repeat_quantify_options(opts: RepeatQuantifyOptions) -> None:
    _save_settings_group(
        "repeat_quantify",
        {
            "engine": opts.engine,
            **opts.params(),
            "blast_save_index": bool(opts.blast_save_index),
        },
    )


def load_consensus_options() -> ConsensusOptions:
    v = _load_settings_group("consensus")
    right_click_action = str(v.get("right_click_action", CONSENSUS_RIGHT_CLICK_RUN))
    if right_click_action not in CONSENSUS_RIGHT_CLICK_ACTIONS:
        right_click_action = CONSENSUS_RIGHT_CLICK_RUN
    return ConsensusOptions(
        threshold=float(v.get("threshold", DEFAULT_CONSENSUS_THRESHOLD)),
        ambiguous=normalize_consensus_ambiguous(
            v.get("ambiguous", DEFAULT_CONSENSUS_AMBIGUOUS)
        ),
        sequence_type=normalize_consensus_sequence_type(
            v.get("sequence_type", MAFFT_SEQUENCE_AUTO)
        ),
        right_click_action=right_click_action,
    )


def save_consensus_options(opts: ConsensusOptions) -> None:
    _save_settings_group(
        "consensus",
        {
            "threshold": opts.threshold,
            "ambiguous": normalize_consensus_ambiguous(opts.ambiguous),
            "sequence_type": normalize_consensus_sequence_type(opts.sequence_type),
            "right_click_action": (
                opts.right_click_action
                if opts.right_click_action in CONSENSUS_RIGHT_CLICK_ACTIONS
                else CONSENSUS_RIGHT_CLICK_RUN
            ),
        },
    )


def load_flexidot_options() -> FlexidotOptions:
    v = _load_settings_group("flexidot")
    wordsize = int(v.get("wordsize", DEFAULT_FLEXIDOT_WORDSIZE))
    wordsize = max(MIN_FLEXIDOT_WORDSIZE, min(MAX_FLEXIDOT_WORDSIZE, wordsize))
    right_click_action = str(v.get("right_click_action", FLEXIDOT_RIGHT_CLICK_RUN))
    if right_click_action not in FLEXIDOT_RIGHT_CLICK_ACTIONS:
        right_click_action = FLEXIDOT_RIGHT_CLICK_RUN
    return FlexidotOptions(
        wordsize=wordsize,
        sequence_type=normalize_flexidot_sequence_type(
            v.get("sequence_type", DEFAULT_FLEXIDOT_SEQUENCE_TYPE)
        ),
        lcs_shading=str(v.get("lcs_shading", "true")).lower() in {"1", "true", "yes"},
        use_gff=str(v.get("use_gff", "false")).lower() in {"1", "true", "yes"},
        output_format=normalize_flexidot_output_format(
            str(v.get("output_format", DEFAULT_FLEXIDOT_OUTPUT_FORMAT))
        ),
        substitution_count=max(0, min(20, int(v.get("substitution_count", 0)))),
        wobble_conversion=str(v.get("wobble_conversion", "false")).lower()
        in {"1", "true", "yes"},
        max_n_percentage=max(0.0, min(100.0, float(v.get("max_n_percentage", 10.0)))),
        include_reverse=str(v.get("include_reverse", "true")).lower()
        in {"1", "true", "yes"},
        only_vs_first=str(v.get("only_vs_first", "false")).lower()
        in {"1", "true", "yes"},
        line_width=max(0.1, min(20.0, float(v.get("line_width", 1.0)))),
        length_scaling=str(v.get("length_scaling", "false")).lower()
        in {"1", "true", "yes"},
        combined=str(v.get("combined", "true")).lower() in {"1", "true", "yes"},
        # Comparison mode is always automatic (not exposed in the dialog).
        modes=DEFAULT_FLEXIDOT_MODES,
        right_click_action=right_click_action,
        show_process_window=str(v.get("show_process_window", "true" if DEFAULT_SHOW_PROCESS_WINDOW else "false")).lower()
        in {"1", "true", "yes"},
        close_process_window_when_finished=str(
            v.get("close_process_window_when_finished", "true" if DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED else "false")
        ).lower()
        in {"1", "true", "yes"},
    )


def save_flexidot_options(opts: FlexidotOptions) -> None:
    _save_settings_group(
        "flexidot",
        {
            "wordsize": opts.wordsize,
            "sequence_type": normalize_flexidot_sequence_type(opts.sequence_type),
            "lcs_shading": opts.lcs_shading,
            "use_gff": opts.use_gff,
            "output_format": normalize_flexidot_output_format(opts.output_format),
            "substitution_count": opts.substitution_count,
            "wobble_conversion": opts.wobble_conversion,
            "max_n_percentage": opts.max_n_percentage,
            "include_reverse": opts.include_reverse,
            "only_vs_first": opts.only_vs_first,
            "line_width": opts.line_width,
            "length_scaling": opts.length_scaling,
            "combined": opts.combined,
            "modes": ",".join(opts.modes),
            "right_click_action": (
                opts.right_click_action
                if opts.right_click_action in FLEXIDOT_RIGHT_CLICK_ACTIONS
                else FLEXIDOT_RIGHT_CLICK_RUN
            ),
            "show_process_window": opts.show_process_window,
            "close_process_window_when_finished": (
                opts.close_process_window_when_finished
            ),
        },
    )


def load_te_annotation_options() -> TeAnnotationJobOptions:
    from repeatviewer.backends.auto_annotation_settings import load_annotation_run_options

    v = _load_settings_group("te_annotation")
    right_click_action = str(v.get("right_click_action", TE_ANNOTATION_RIGHT_CLICK_RUN))
    if right_click_action not in TE_ANNOTATION_RIGHT_CLICK_ACTIONS:
        right_click_action = TE_ANNOTATION_RIGHT_CLICK_RUN
    raw_modes = str(v.get("auto_modes", "")).strip()
    if "auto_modes" not in v:
        auto_modes = DEFAULT_TE_ANNOTATION_AUTO_MODES
    elif not raw_modes:
        auto_modes = ()
    else:
        auto_modes = tuple(part.strip() for part in raw_modes.split(",") if part.strip())
    if "similarity_database_names" not in v:
        similarity_database_names: tuple[str, ...] | None = None
    else:
        raw_dbs = str(v.get("similarity_database_names", "")).strip()
        similarity_database_names = tuple(
            part.strip() for part in raw_dbs.split(",") if part.strip()
        )
    _threads, show, close = load_annotation_run_options()
    return TeAnnotationJobOptions(
        auto_modes=auto_modes,
        similarity_database_names=similarity_database_names,
        right_click_action=right_click_action,
        show_process_window=show,
        close_process_window_when_finished=close,
    )


def save_te_annotation_options(opts: TeAnnotationJobOptions) -> None:
    from PySide6.QtCore import QSettings

    from repeatviewer.backends.auto_annotation_settings import (
        load_annotation_run_options,
        save_annotation_run_options,
    )

    s = QSettings()
    s.beginGroup("te_annotation")
    s.setValue("auto_modes", ",".join(opts.auto_modes))
    s.setValue(
        "right_click_action",
        (
            opts.right_click_action
            if opts.right_click_action in TE_ANNOTATION_RIGHT_CLICK_ACTIONS
            else TE_ANNOTATION_RIGHT_CLICK_RUN
        ),
    )
    if opts.similarity_database_names is None:
        s.remove("similarity_database_names")
    else:
        s.setValue(
            "similarity_database_names",
            ",".join(opts.similarity_database_names),
        )
    s.endGroup()
    threads, _show, _close = load_annotation_run_options()
    save_annotation_run_options(
        num_threads=threads,
        show_process_window=opts.show_process_window,
        close_process_window_when_finished=opts.close_process_window_when_finished,
    )


def load_te_confusion_matrix_benchmark_options() -> TeConfusionMatrixBenchmarkOptions:
    v = _load_settings_group("te_confusion_benchmark")
    home = str(v.get("tetrimmer_home", "")).strip()
    if not home:
        from repeatviewer.backends.te_benchmark import default_tetrimmer_home

        detected = default_tetrimmer_home()
        if detected is not None:
            home = str(detected)
    return TeConfusionMatrixBenchmarkOptions(
        tetrimmer_home=home,
        reference_out=str(v.get("reference_out", "")).strip(),
        less_classification=str(v.get("less_classification", "false")).lower()
        in {"1", "true", "yes"},
        debug=str(v.get("debug", "false")).lower() in {"1", "true", "yes"},
    )


def save_te_confusion_matrix_benchmark_options(opts: TeConfusionMatrixBenchmarkOptions) -> None:
    _save_settings_group(
        "te_confusion_benchmark",
        {
            "tetrimmer_home": opts.tetrimmer_home,
            "reference_out": opts.reference_out,
            "less_classification": opts.less_classification,
            "debug": opts.debug,
        },
    )


def load_te_library_benchmark_options() -> TeLibraryBenchmarkOptions:
    v = _load_settings_group("te_library_benchmark")
    return TeLibraryBenchmarkOptions(
        query_fasta=str(v.get("query_fasta", "")).strip(),
        standard_fasta=str(v.get("standard_fasta", "")).strip(),
        rm_out=str(v.get("rm_out", "")).strip(),
        run_repeatmasker=str(v.get("run_repeatmasker", "true")).lower() in {"1", "true", "yes"},
        threads=max(1, int(v.get("threads", DEFAULT_TE_BENCHMARK_THREADS))),
    )


def save_te_library_benchmark_options(opts: TeLibraryBenchmarkOptions) -> None:
    _save_settings_group(
        "te_library_benchmark",
        {
            "query_fasta": opts.query_fasta,
            "standard_fasta": opts.standard_fasta,
            "rm_out": opts.rm_out,
            "run_repeatmasker": opts.run_repeatmasker,
            "threads": opts.threads,
        },
    )


def save_msa_clean_options(opts: MsaCleanOptions) -> None:
    from PySide6.QtCore import QSettings

    s = QSettings()
    s.beginGroup("msa_clean")
    s.remove("default_method")
    s.setValue("use_divergence", opts.use_divergence)
    s.setValue("use_gap_row", opts.use_gap_row)
    s.setValue("use_gappy_columns", opts.use_gappy_columns)
    s.setValue("method_order", list(opts.method_order))
    s.setValue(
        "right_click_action",
        (
            opts.right_click_action
            if opts.right_click_action in MSA_CLEAN_RIGHT_CLICK_ACTIONS
            else MSA_CLEAN_RIGHT_CLICK_COMBINED
        ),
    )
    s.setValue("show_process_window", opts.show_process_window)
    s.setValue(
        "close_process_window_when_finished",
        opts.close_process_window_when_finished,
    )
    g = opts.gappy
    s.setValue("gappy_gap_threshold", g.gap_threshold)
    s.setValue("gappy_simi_check_gap_thr", g.simi_check_gap_thr)
    s.setValue("gappy_similarity_thr", g.similarity_thr)
    s.setValue("gappy_min_nucleotide", g.min_nucleotide)
    d = opts.divergence
    s.setValue("div_threshold", d.threshold)
    s.setValue("div_window_size", d.window_size)
    gr = opts.gap_row
    s.setValue("gap_row_gap_threshold", gr.gap_threshold)
    s.setValue("gap_row_window_size", gr.window_size)
    s.endGroup()


# --------------------------------------------------------------------------
# Tree (phylogenetic tree builder)
# --------------------------------------------------------------------------

TREE_ENGINE_IQTREE = "iqtree"
TREE_ENGINE_FASTME = "fastme"
TREE_ENGINE_FASTTREE = "fasttree"
TREE_ENGINES = (TREE_ENGINE_IQTREE, TREE_ENGINE_FASTME, TREE_ENGINE_FASTTREE)
TREE_DIALOG_TABS = TREE_ENGINES + (DIALOG_TAB_OPTIONS,)

TREE_ALIGN_AUTO = "auto"      # align unless the input is already an alignment
TREE_ALIGN_ALWAYS = "always"
TREE_ALIGN_NEVER = "never"    # input must already be aligned
TREE_ALIGN_MODES = (TREE_ALIGN_AUTO, TREE_ALIGN_ALWAYS, TREE_ALIGN_NEVER)

TREE_COVERAGE_MEDIAN = "median"    # length filter relative to the median length
TREE_COVERAGE_LONGEST = "longest"  # ... or to the longest sequence
TREE_COVERAGE_REFERENCES = (TREE_COVERAGE_MEDIAN, TREE_COVERAGE_LONGEST)

TREE_ROOT_NONE = "none"
TREE_ROOT_OUTGROUP = "outgroup"
TREE_ROOT_MIDPOINT = "midpoint"
TREE_ROOT_MODES = (TREE_ROOT_NONE, TREE_ROOT_OUTGROUP, TREE_ROOT_MIDPOINT)

TREE_MODEL_FINDER = "MFP"      # IQ-TREE ModelFinder, best model by BIC

TREE_RIGHT_CLICK_RUN = "run"
TREE_RIGHT_CLICK_DIALOG = "dialog"
TREE_RIGHT_CLICK_ACTIONS = (TREE_RIGHT_CLICK_RUN, TREE_RIGHT_CLICK_DIALOG)


TREE_SEQ_AS_IS = "as_is"
TREE_SEQ_RT = "rt_domain"
TREE_SEQ_MODES = (TREE_SEQ_AS_IS, TREE_SEQ_RT)


@dataclass
class TreeOptions:
    """Phylogenetic tree building (Tools → Tree).

    Defaults follow current practice for a publishable tree: MAFFT L-INS-i with
    strand adjustment, IQ-TREE with ModelFinder (BIC), several independent
    searches, SH-aLRT + UFBoot2 (-bnni) branch support and a fixed seed.
    """

    # Alignment
    align_mode: str = TREE_ALIGN_AUTO
    mafft_strategy: str = "choose"     # L-INS-i up to 200 sequences, else FFT-NS-2
    adjust_direction: bool = True
    rotate_monomers: bool = False      # satellite monomers: common start before aligning
    # Sequences: TREE_SEQ_AS_IS, or TREE_SEQ_RT = RT domain found in the DNA (REXdb)
    sequence_mode: str = TREE_SEQ_AS_IS
    rexdb_fasta: str = ""              # user's REXdb protein FASTA (.classification beside it)
    rexdb_references: int = 3          # REXdb RTs added per lineage (0 = none)
    mafft_gap_open: float = 1.53
    mafft_offset: float = 0.0
    min_coverage: float = 0.7          # fraction of the reference length; 0 = off
    coverage_reference: str = TREE_COVERAGE_MEDIAN
    # Tree
    engine: str = TREE_ENGINE_IQTREE
    # Which Tree window tab was open last (iqtree / fastme / fasttree / options).
    last_dialog_tab: str = TREE_ENGINE_IQTREE
    model: str = TREE_MODEL_FINDER
    runs: int = 5
    seed: int = 12345
    threads: int = field(default_factory=_default_tool_threads)
    # Support and rooting
    ufboot: int = 1000
    alrt: int = 1000
    bnni: bool = True                   # IQ-TREE --bnni with UFBoot
    iqtree_fast: bool = False           # IQ-TREE --fast (FastTree-like search)
    distance_bootstrap: int = 100      # FastME bootstrap replicates
    root_mode: str = TREE_ROOT_NONE
    outgroup: list[str] = field(default_factory=list)
    # Options tab
    right_click_action: str = TREE_RIGHT_CLICK_RUN
    show_process_window: bool = DEFAULT_SHOW_PROCESS_WINDOW
    close_process_window_when_finished: bool = DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED
    right_click_engine: str = TREE_ENGINE_IQTREE   # program run by right-click
    use_system_mafft: bool = False
    use_system_iqtree: bool = False    # per program: system PATH only
    use_system_fastme: bool = False
    use_system_fasttree: bool = False
    use_system_path: bool = False      # old single switch: counts for every tree program
    # Methylation check: also build a tree without CG / CHG columns and compare.
    # Set in the Tree tab (saved) or for one run by Tools → Tree → masked run.
    methylation_mask: str = ""         # "", "CG" or "CG_CHG"


def load_tree_options() -> TreeOptions:
    v = _load_settings_group("tree")
    d = TreeOptions()

    def _pick(key: str, allowed: tuple[str, ...], default: str) -> str:
        value = str(v.get(key, default))
        return value if value in allowed else default

    def _num(key: str, default, cast, lo, hi):
        try:
            return max(lo, min(hi, cast(v.get(key, default))))
        except (TypeError, ValueError):
            return default

    outgroup = v.get("outgroup") or []
    if isinstance(outgroup, str):
        outgroup = [x for x in outgroup.split(",") if x]
    engine = _pick("engine", TREE_ENGINES, d.engine)
    legacy_system = _to_bool(v.get("use_system_path"), False)
    return TreeOptions(
        align_mode=_pick("align_mode", TREE_ALIGN_MODES, d.align_mode),
        mafft_strategy=str(v.get("mafft_strategy", d.mafft_strategy)),
        adjust_direction=_to_bool(v.get("adjust_direction"), d.adjust_direction),
        rotate_monomers=_to_bool(v.get("rotate_monomers"), d.rotate_monomers),
        methylation_mask=_pick("methylation_mask", ("", "CG", "CG_CHG"), d.methylation_mask),
        sequence_mode=_pick("sequence_mode", TREE_SEQ_MODES, d.sequence_mode),
        rexdb_fasta=str(v.get("rexdb_fasta", d.rexdb_fasta) or ""),
        rexdb_references=_num("rexdb_references", d.rexdb_references, int, 0, 50),
        mafft_gap_open=_num("mafft_gap_open", d.mafft_gap_open, float, 0.0, 20.0),
        mafft_offset=_num("mafft_offset", d.mafft_offset, float, 0.0, 10.0),
        min_coverage=_num("min_coverage", d.min_coverage, float, 0.0, 1.0),
        coverage_reference=_pick("coverage_reference", TREE_COVERAGE_REFERENCES,
                                 d.coverage_reference),
        engine=engine,
        last_dialog_tab=_pick("last_dialog_tab", TREE_DIALOG_TABS, engine),
        model=str(v.get("model", d.model) or d.model).strip(),
        runs=_num("runs", d.runs, int, 1, 100),
        seed=_num("seed", d.seed, int, 1, 2_000_000_000),
        threads=_num("threads", _default_tool_threads(), int, 1, 256),
        ufboot=_num("ufboot", d.ufboot, int, 0, 100_000),
        alrt=_num("alrt", d.alrt, int, 0, 100_000),
        bnni=_to_bool(v.get("bnni"), d.bnni),
        iqtree_fast=_to_bool(v.get("iqtree_fast"), d.iqtree_fast),
        distance_bootstrap=_num("distance_bootstrap", d.distance_bootstrap, int, 0, 10_000),
        root_mode=_pick("root_mode", TREE_ROOT_MODES, d.root_mode),
        outgroup=[str(x) for x in outgroup],
        right_click_action=_pick(
            "right_click_action", TREE_RIGHT_CLICK_ACTIONS, d.right_click_action
        ),
        show_process_window=_to_bool(v.get("show_process_window"), d.show_process_window),
        close_process_window_when_finished=_to_bool(
            v.get("close_process_window_when_finished"),
            d.close_process_window_when_finished,
        ),
        right_click_engine=_pick("right_click_engine", TREE_ENGINES, engine),
        use_system_mafft=_to_bool(v.get("use_system_mafft"), d.use_system_mafft),
        # settings saved before the per-program boxes: the one switch meant all three
        use_system_iqtree=_to_bool(v.get("use_system_iqtree"), legacy_system),
        use_system_fastme=_to_bool(v.get("use_system_fastme"), legacy_system),
        use_system_fasttree=_to_bool(v.get("use_system_fasttree"), legacy_system),
        use_system_path=False,
    )


def save_tree_options(opts: TreeOptions) -> None:
    _save_settings_group(
        "tree",
        {
            "align_mode": opts.align_mode,
            "mafft_strategy": opts.mafft_strategy,
            "adjust_direction": opts.adjust_direction,
            "rotate_monomers": opts.rotate_monomers,
            "methylation_mask": opts.methylation_mask,
            "sequence_mode": opts.sequence_mode,
            "rexdb_fasta": opts.rexdb_fasta,
            "rexdb_references": opts.rexdb_references,
            "mafft_gap_open": opts.mafft_gap_open,
            "mafft_offset": opts.mafft_offset,
            "min_coverage": opts.min_coverage,
            "coverage_reference": opts.coverage_reference,
            "engine": opts.engine,
            "last_dialog_tab": (
                opts.last_dialog_tab
                if opts.last_dialog_tab in TREE_DIALOG_TABS
                else opts.engine
            ),
            "model": opts.model,
            "runs": opts.runs,
            "seed": opts.seed,
            "threads": opts.threads,
            "ufboot": opts.ufboot,
            "alrt": opts.alrt,
            "bnni": opts.bnni,
            "iqtree_fast": opts.iqtree_fast,
            "distance_bootstrap": opts.distance_bootstrap,
            "root_mode": opts.root_mode,
            "outgroup": list(opts.outgroup),
            "right_click_action": opts.right_click_action,
            "show_process_window": opts.show_process_window,
            "close_process_window_when_finished": (
                opts.close_process_window_when_finished
            ),
            "right_click_engine": opts.right_click_engine,
            "use_system_mafft": opts.use_system_mafft,
            "use_system_iqtree": opts.use_system_iqtree or opts.use_system_path,
            "use_system_fastme": opts.use_system_fastme or opts.use_system_path,
            "use_system_fasttree": opts.use_system_fasttree or opts.use_system_path,
            "use_system_path": False,
        },
    )


# --------------------------------------------------------------------------
# Map to Reference (Tools → Map to Reference)
# --------------------------------------------------------------------------

MAP_ENGINE_AUTO = "auto"
MAP_ENGINE_MAFFT = "mafft"  # legacy saved value; migrated to auto on load
MAP_ENGINE_MINIMAP2 = "minimap2"
MAP_ENGINE_HISAT2 = "hisat2"
MAP_ENGINE_BOWTIE2 = "bowtie2"  # legacy; migrated to hisat2 on load
MAP_ENGINES = (MAP_ENGINE_AUTO, MAP_ENGINE_MINIMAP2, MAP_ENGINE_HISAT2)
MAP_DIALOG_TABS = (MAP_ENGINE_MINIMAP2, MAP_ENGINE_HISAT2, DIALOG_TAB_OPTIONS)

# Kept for old saved MapOptions; HISAT2 no longer uses Bowtie2 sensitivity presets.
MAP_SENS_FAST = "fast"
MAP_SENS_MEDIUM = "medium"
MAP_SENS_HIGH = "high"
MAP_SENSITIVITIES = (MAP_SENS_FAST, MAP_SENS_MEDIUM, MAP_SENS_HIGH)

# HISAT2 library / alignment mode (Settings).
MAP_HISAT2_RNA = "rna"  # spliced (default)
MAP_HISAT2_DNA = "dna"  # --no-spliced-alignment
MAP_HISAT2_LIBRARIES = (MAP_HISAT2_RNA, MAP_HISAT2_DNA)

# HISAT2 read layout (Settings): single-end (-U) vs Illumina paired-end (-1/-2).
MAP_HISAT2_READS_SINGLE = "single"
MAP_HISAT2_READS_PAIRED = "paired"
MAP_HISAT2_READS = (MAP_HISAT2_READS_SINGLE, MAP_HISAT2_READS_PAIRED)

# HISAT2 --rna-strandness (empty = unstranded / default).
MAP_HISAT2_STRAND_UNSTRANDED = ""
MAP_HISAT2_STRAND_F = "F"
MAP_HISAT2_STRAND_R = "R"
MAP_HISAT2_STRAND_FR = "FR"
MAP_HISAT2_STRAND_RF = "RF"
MAP_HISAT2_STRANDS = (
    MAP_HISAT2_STRAND_UNSTRANDED,
    MAP_HISAT2_STRAND_F,
    MAP_HISAT2_STRAND_R,
    MAP_HISAT2_STRAND_FR,
    MAP_HISAT2_STRAND_RF,
)

# A sequence that fits equally well in several places (common in satellite arrays).
MAP_MULTI_BEST = "best"    # one place (the mapper's best; ties at random)
MAP_MULTI_ALL = "all"      # every equally good place
MAP_MULTI_NONE = "none"    # none: listed as ambiguous with the unused sequences
MAP_MULTI_MODES = (MAP_MULTI_BEST, MAP_MULTI_ALL, MAP_MULTI_NONE)

MAP_MINIMAP2_PRESET_AUTO = "auto"  # legacy saved value; migrated to map-ont on load
MAP_MINIMAP2_PRESETS = (
    "map-ont",
    "map-hifi",
    "map-pb",
    "map-iclr",
    "lr:hq",
    "asm5",
    "asm10",
    "asm20",
    "sr",
    "splice",
    "splice:hq",
    "splice:sr",
)

# Short combo labels (fit TOOL_PARAMS_FIELD_WIDTH); details live in the ? help.
MAP_MINIMAP2_PRESET_LABELS: dict[str, str] = {
    "map-ont": "map-ont (ONT genomic)",
    "map-hifi": "map-hifi (PacBio HiFi)",
    "map-pb": "map-pb (PacBio CLR)",
    "map-iclr": "map-iclr (Illumina CLR)",
    "lr:hq": "lr:hq (Nanopore Q20)",
    "asm5": "asm5 (~0.1% divergence)",
    "asm10": "asm10 (~1% divergence)",
    "asm20": "asm20 (~5% divergence)",
    "sr": "sr (short genomic)",
    "splice": "splice (long RNA/cDNA)",
    "splice:hq": "splice:hq (Iso-seq / HQ RNA)",
    "splice:sr": "splice:sr (short RNA-seq)",
}

# minimap2 -k / -w for each -x preset (from minimap2.1). map-ont uses the
# tool defaults (15 / 10). Auto resolves at run time — see auto_minimap2_preset.
MAP_MINIMAP2_PRESET_KW: dict[str, tuple[int, int]] = {
    "map-ont": (15, 10),
    "map-pb": (19, 10),       # -H -k19
    "map-hifi": (19, 19),     # via lr:hq
    "lr:hq": (19, 19),
    "map-iclr": (19, 10),
    "asm5": (19, 19),
    "asm10": (19, 19),
    "asm20": (19, 10),
    "sr": (21, 11),
    "splice": (15, 5),        # -k15 -w5 --splice …
    "splice:hq": (15, 5),     # -xsplice + HQ scoring
    "splice:sr": (15, 5),     # -xsplice:hq + short-read flags
}

# minimap2 alignment file written next to the .rvmap (PAF is the viewer default).
MAP_MINIMAP2_OUT_PAF = "paf"
MAP_MINIMAP2_OUT_SAM = "sam"
MAP_MINIMAP2_OUTPUTS = (MAP_MINIMAP2_OUT_PAF, MAP_MINIMAP2_OUT_SAM)

# Consensus threshold, as a fraction of the sequences covering a position.
# 0 = majority (the most common base); higher values give IUPAC codes.
MAP_CONSENSUS_THRESHOLDS = (0.0, 0.5, 0.75, 0.95, 1.0)
MAP_NO_COVERAGE_N = "N"
MAP_NO_COVERAGE_REFERENCE = "reference"
MAP_NO_COVERAGE_SKIP = "skip"
MAP_NO_COVERAGE_MODES = (MAP_NO_COVERAGE_N, MAP_NO_COVERAGE_REFERENCE, MAP_NO_COVERAGE_SKIP)

MAP_RIGHT_CLICK_RUN = "run"
MAP_RIGHT_CLICK_DIALOG = "dialog"
MAP_RIGHT_CLICK_ACTIONS = (MAP_RIGHT_CLICK_RUN, MAP_RIGHT_CLICK_DIALOG)


@dataclass
class MapOptions:
    """Map sequences or reads to a reference (Tools → Mapping).

    minimap2 maps long reads or assemblies to a genome; HISAT2 maps short
    DNA or RNA-seq reads. ``auto`` picks between them by input.
    """

    reference_fasta: str = ""
    engine: str = MAP_ENGINE_AUTO
    # Which Map window tab was open last (minimap2 / hisat2 / options).
    last_dialog_tab: str = MAP_ENGINE_MINIMAP2
    sensitivity: str = MAP_SENS_MEDIUM  # legacy Bowtie2; unused by HISAT2
    multiple_hits: str = MAP_MULTI_BEST
    minimap2_preset: str = "map-ont"
    minimap2_output: str = MAP_MINIMAP2_OUT_PAF
    minimap2_include_cigar: bool = True   # PAF -c
    minimap2_include_cs: bool = True      # PAF --cs=long (needed for contig view)
    # minimap2 -k / -w (concrete numbers; updated from the chosen -x preset)
    minimap2_kmer: int = 15
    minimap2_window: int = 10
    # HISAT2 Settings / Advanced
    hisat2_library: str = MAP_HISAT2_RNA
    hisat2_reads: str = MAP_HISAT2_READS_SINGLE
    hisat2_strandness: str = MAP_HISAT2_STRAND_UNSTRANDED
    hisat2_dta: bool = False              # --dta (StringTie-oriented)
    hisat2_min_intron: int = 20           # --min-intronlen
    hisat2_max_intron: int = 500_000      # --max-intronlen
    hisat2_no_softclip: bool = False      # --no-softclip
    # Advanced: a mapped sequence is kept only when it passes all of these.
    min_identity: float = 0.70         # fraction of aligned columns that match
    min_aligned_bp: int = 25
    min_mapq: int = 0                  # minimap2 / HISAT2
    threads: int = field(default_factory=_default_tool_threads)
    # Output (RepeatViewer contig view / sidecars)
    consensus_threshold: float = 0.0   # 0 = majority
    no_coverage: str = MAP_NO_COVERAGE_N
    min_variant_freq: float = 0.0      # variants table: 0 = every difference
    save_consensus: bool = True
    save_mapped: bool = True
    save_unused: bool = True
    save_variants: bool = True
    lift_annotations: bool = True
    max_display_reads: int = 2000
    # Options tab
    right_click_action: str = MAP_RIGHT_CLICK_RUN
    show_process_window: bool = DEFAULT_SHOW_PROCESS_WINDOW
    close_process_window_when_finished: bool = DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED
    use_system_path: bool = False


def load_map_options() -> MapOptions:
    v = _load_settings_group("map_reference")
    d = MapOptions()

    def _pick(key: str, allowed: tuple, default):
        value = v.get(key, default)
        return value if value in allowed else default

    def _num(key: str, default, cast, lo, hi):
        try:
            return max(lo, min(hi, cast(v.get(key, default))))
        except (TypeError, ValueError):
            return default

    threshold = _num("consensus_threshold", d.consensus_threshold, float, 0.0, 1.0)
    engine = v.get("engine", d.engine)
    if engine == MAP_ENGINE_BOWTIE2:
        engine = MAP_ENGINE_HISAT2
    elif engine == MAP_ENGINE_MAFFT or engine not in MAP_ENGINES:
        engine = MAP_ENGINE_AUTO
    mm_preset = v.get("minimap2_preset", d.minimap2_preset)
    if mm_preset == MAP_MINIMAP2_PRESET_AUTO or mm_preset not in MAP_MINIMAP2_PRESETS:
        mm_preset = "map-ont"
    last_tab = v.get("last_dialog_tab", None)
    if last_tab not in MAP_DIALOG_TABS:
        last_tab = (
            engine
            if engine in (MAP_ENGINE_MINIMAP2, MAP_ENGINE_HISAT2)
            else MAP_ENGINE_MINIMAP2
        )
    return MapOptions(
        reference_fasta=str(v.get("reference_fasta", d.reference_fasta) or ""),
        engine=engine,
        last_dialog_tab=last_tab,
        sensitivity=_pick("sensitivity", MAP_SENSITIVITIES, d.sensitivity),
        multiple_hits=_pick("multiple_hits", MAP_MULTI_MODES, d.multiple_hits),
        minimap2_preset=mm_preset,
        minimap2_output=_pick("minimap2_output", MAP_MINIMAP2_OUTPUTS, d.minimap2_output),
        minimap2_include_cigar=_to_bool(
            v.get("minimap2_include_cigar"), d.minimap2_include_cigar
        ),
        minimap2_include_cs=_to_bool(v.get("minimap2_include_cs"), d.minimap2_include_cs),
        minimap2_kmer=_num("minimap2_kmer", d.minimap2_kmer, int, 1, 28),
        minimap2_window=_num("minimap2_window", d.minimap2_window, int, 1, 1000),
        hisat2_library=_pick("hisat2_library", MAP_HISAT2_LIBRARIES, d.hisat2_library),
        hisat2_reads=_pick("hisat2_reads", MAP_HISAT2_READS, d.hisat2_reads),
        hisat2_strandness=_pick(
            "hisat2_strandness", MAP_HISAT2_STRANDS, d.hisat2_strandness
        ),
        hisat2_dta=_to_bool(v.get("hisat2_dta"), d.hisat2_dta),
        hisat2_min_intron=_num("hisat2_min_intron", d.hisat2_min_intron, int, 1, 1_000_000),
        hisat2_max_intron=_num(
            "hisat2_max_intron", d.hisat2_max_intron, int, 1, 10_000_000
        ),
        hisat2_no_softclip=_to_bool(v.get("hisat2_no_softclip"), d.hisat2_no_softclip),
        min_identity=_num("min_identity", d.min_identity, float, 0.0, 1.0),
        min_aligned_bp=_num("min_aligned_bp", d.min_aligned_bp, int, 1, 1_000_000),
        min_mapq=_num("min_mapq", d.min_mapq, int, 0, 255),
        threads=_num("threads", _default_tool_threads(), int, 1, 256),
        consensus_threshold=threshold if threshold in MAP_CONSENSUS_THRESHOLDS else 0.0,
        no_coverage=_pick("no_coverage", MAP_NO_COVERAGE_MODES, d.no_coverage),
        min_variant_freq=_num("min_variant_freq", d.min_variant_freq, float, 0.0, 1.0),
        save_consensus=_to_bool(v.get("save_consensus"), d.save_consensus),
        save_mapped=_to_bool(v.get("save_mapped"), d.save_mapped),
        save_unused=_to_bool(v.get("save_unused"), d.save_unused),
        save_variants=_to_bool(v.get("save_variants"), d.save_variants),
        lift_annotations=_to_bool(v.get("lift_annotations"), d.lift_annotations),
        max_display_reads=_num("max_display_reads", d.max_display_reads, int, 50, 100_000),
        right_click_action=_pick("right_click_action", MAP_RIGHT_CLICK_ACTIONS,
                                 d.right_click_action),
        show_process_window=_to_bool(v.get("show_process_window"), d.show_process_window),
        close_process_window_when_finished=_to_bool(
            v.get("close_process_window_when_finished"),
            d.close_process_window_when_finished,
        ),
        use_system_path=_to_bool(v.get("use_system_path"), d.use_system_path),
    )


def save_map_options(opts: MapOptions) -> None:
    _save_settings_group(
        "map_reference",
        {
            "reference_fasta": opts.reference_fasta,
            "engine": opts.engine,
            "last_dialog_tab": (
                opts.last_dialog_tab
                if opts.last_dialog_tab in MAP_DIALOG_TABS
                else MAP_ENGINE_MINIMAP2
            ),
            "sensitivity": opts.sensitivity,
            "multiple_hits": opts.multiple_hits,
            "minimap2_preset": opts.minimap2_preset,
            "minimap2_output": opts.minimap2_output,
            "minimap2_include_cigar": opts.minimap2_include_cigar,
            "minimap2_include_cs": opts.minimap2_include_cs,
            "minimap2_kmer": opts.minimap2_kmer,
            "minimap2_window": opts.minimap2_window,
            "hisat2_library": opts.hisat2_library,
            "hisat2_reads": opts.hisat2_reads,
            "hisat2_strandness": opts.hisat2_strandness,
            "hisat2_dta": opts.hisat2_dta,
            "hisat2_min_intron": opts.hisat2_min_intron,
            "hisat2_max_intron": opts.hisat2_max_intron,
            "hisat2_no_softclip": opts.hisat2_no_softclip,
            "min_identity": opts.min_identity,
            "min_aligned_bp": opts.min_aligned_bp,
            "min_mapq": opts.min_mapq,
            "threads": opts.threads,
            "consensus_threshold": opts.consensus_threshold,
            "no_coverage": opts.no_coverage,
            "min_variant_freq": opts.min_variant_freq,
            "save_consensus": opts.save_consensus,
            "save_mapped": opts.save_mapped,
            "save_unused": opts.save_unused,
            "save_variants": opts.save_variants,
            "lift_annotations": opts.lift_annotations,
            "max_display_reads": opts.max_display_reads,
            "right_click_action": opts.right_click_action,
            "show_process_window": opts.show_process_window,
            "close_process_window_when_finished": (
                opts.close_process_window_when_finished
            ),
            "use_system_path": opts.use_system_path,
        },
    )
