"""Per-job parameters passed from the UI into curation backends."""

from __future__ import annotations

import json
from dataclasses import dataclass, field


def _to_bool(value, default: bool = False) -> bool:
    if value is None:
        return default
    return str(value).strip().lower() in {"1", "true", "yes", "on"}


def _default_tool_threads() -> int:
    """Live default thread count: the Preferences 'Threads for tools' value."""
    try:
        from repeatviewer.backends.thread_budget import load_batch_workers

        return load_batch_workers()
    except Exception:  # noqa: BLE001
        return 4


DEFAULT_EXTEND_LEFT = 1000
DEFAULT_EXTEND_RIGHT = 1000
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
DEFAULT_BLAST_WRITE_GFF = True
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
BLAST_PROGRAMS = (BLAST_PROGRAM_BLASTN, BLAST_PROGRAM_BLASTP, BLAST_PROGRAM_BLASTX)
# blastx translates a nucleotide query and searches a protein database, so it
# shares blastp's protein-style scoring but takes DNA in.
BLAST_PROGRAMS_PROTEIN_DB = (BLAST_PROGRAM_BLASTP, BLAST_PROGRAM_BLASTX)
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


@dataclass
class TeGenomeMapOptions:
    """TE density plot: bin window, density metric, and chromosome filter."""

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
MUSCLE_MODE_ALIGN = "align"
MUSCLE_MODE_SUPER5 = "super5"
MUSCLE_MODES = (MUSCLE_MODE_ALIGN, MUSCLE_MODE_SUPER5)


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
    # Per-aligner "combine all selected files into one MSA" choice.
    mafft_combined: bool = False
    muscle_combined: bool = False
    # Which aligner a right-click on the Align toolbar button runs.
    right_click_aligner: str = ALIGNER_MAFFT
    show_process_window: bool = False


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
DEFAULT_FASTGA_SEED_FREQ = 200
DEFAULT_FASTGA_CHAIN_COVERAGE = 85
DEFAULT_FASTGA_CHAIN_GAP = 1000
DEFAULT_FASTGA_SYMMETRIC = True
DEFAULT_FASTGA_KEEP_INDEX = True
DEFAULT_FASTGA_NO_SOFT_MASK = False
DEFAULT_FASTGA_SHOW_PROCESS_WINDOW = False
DEFAULT_FASTGA_PAIRWISE = False


@dataclass
class FastGaOptions:
    """FastGA whole-genome alignment.

    By default each selected FASTA is self-compared one by one. Pairwise mode
    ("Do each pair comparison") compares every unique pair among the selected
    files and writes one combined alignment file.
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
    no_soft_mask: bool = DEFAULT_FASTGA_NO_SOFT_MASK
    pairwise: bool = DEFAULT_FASTGA_PAIRWISE
    temp_dir: str = ""
    right_click_action: str = FASTGA_RIGHT_CLICK_RUN
    show_process_window: bool = DEFAULT_FASTGA_SHOW_PROCESS_WINDOW


def load_fastga_options() -> FastGaOptions:
    v = _load_settings_group("fastga")
    fmt = str(v.get("output_format", FASTGA_FORMAT_PAF))
    if fmt not in FASTGA_OUTPUT_FORMATS:
        fmt = FASTGA_FORMAT_PAF
    right_click_action = str(v.get("right_click_action", FASTGA_RIGHT_CLICK_RUN))
    if right_click_action not in FASTGA_RIGHT_CLICK_ACTIONS:
        right_click_action = FASTGA_RIGHT_CLICK_RUN
    return FastGaOptions(
        output_format=fmt,
        threads=max(1, int(v.get("threads", _default_tool_threads()))),
        min_identity=max(
            0.0,
            min(1.0, float(v.get("min_identity", DEFAULT_FASTGA_MIN_IDENTITY))),
        ),
        min_length=max(1, int(v.get("min_length", DEFAULT_FASTGA_MIN_LENGTH))),
        seed_freq=max(1, int(v.get("seed_freq", DEFAULT_FASTGA_SEED_FREQ))),
        chain_coverage=max(
            1, int(v.get("chain_coverage", DEFAULT_FASTGA_CHAIN_COVERAGE))
        ),
        chain_gap=max(1, int(v.get("chain_gap", DEFAULT_FASTGA_CHAIN_GAP))),
        symmetric=_to_bool(v.get("symmetric"), DEFAULT_FASTGA_SYMMETRIC),
        keep_index=_to_bool(v.get("keep_index"), DEFAULT_FASTGA_KEEP_INDEX),
        no_soft_mask=_to_bool(v.get("no_soft_mask"), DEFAULT_FASTGA_NO_SOFT_MASK),
        pairwise=_to_bool(v.get("pairwise"), DEFAULT_FASTGA_PAIRWISE),
        temp_dir=str(v.get("temp_dir", "") or ""),
        right_click_action=right_click_action,
        show_process_window=_to_bool(
            v.get("show_process_window"), DEFAULT_FASTGA_SHOW_PROCESS_WINDOW
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
            "no_soft_mask": opts.no_soft_mask,
            "pairwise": opts.pairwise,
            "temp_dir": opts.temp_dir,
            "right_click_action": opts.right_click_action,
            "show_process_window": opts.show_process_window,
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
    # Which program a right-click on the BLAST toolbar button runs.
    right_click_program: str = BLAST_PROGRAM_BLASTN
    show_process_window: bool = False

    @property
    def active_database_name(self) -> str | None:
        if self.program == BLAST_PROGRAM_BLASTP:
            return self.blastp_database_name
        if self.program == BLAST_PROGRAM_BLASTX:
            return self.blastx_database_name
        return self.blast_database_name


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
    right_click_action: str = CONSENSUS_RIGHT_CLICK_RUN


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
    show_process_window: bool = False


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


DEFAULT_TEAID_EVALUE = DEFAULT_REPORT_EVALUE
DEFAULT_TEAID_THREADS = 1
DEFAULT_TEAID_SELF_BLAST_WORD_SIZE = 6
# Text written on the consensus (x-)axis of the TE-Aid panels. Configurable so
# non-TE users (e.g. tandem repeats) can relabel it.
DEFAULT_TEAID_X_AXIS_LABEL = "TE consensus (bp)"
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
    num_threads: int = DEFAULT_TEAID_THREADS
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
        right_click_program=(
            str(v.get("right_click_program", BLAST_PROGRAM_BLASTN))
            if str(v.get("right_click_program", BLAST_PROGRAM_BLASTN)) in BLAST_PROGRAMS
            else BLAST_PROGRAM_BLASTN
        ),
        show_process_window=_to_bool(v.get("show_process_window"), False),
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
            "right_click_program": opts.right_click_program,
            "show_process_window": opts.show_process_window,
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
    return TeaAidOptions(
        e_value=normalize_blast_evalue(
            float(v.get("e_value", DEFAULT_TEAID_EVALUE)),
            default=DEFAULT_TEAID_EVALUE,
        ),
        self_blast_word_size=word_size,
        num_threads=max(1, int(v.get("num_threads", DEFAULT_TEAID_THREADS))),
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
    right_click_aligner = str(v.get("right_click_aligner", ALIGNER_MAFFT))
    if right_click_aligner not in ALIGNERS:
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
        mafft_combined=_flag("mafft_combined"),
        muscle_combined=_flag("muscle_combined"),
        right_click_aligner=right_click_aligner,
        show_process_window=_flag("show_process_window"),
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
            "mafft_combined": opts.mafft_combined,
            "muscle_combined": opts.muscle_combined,
            "right_click_aligner": opts.right_click_aligner,
            "show_process_window": opts.show_process_window,
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
        right_click_action=right_click_action,
    )


def save_consensus_options(opts: ConsensusOptions) -> None:
    _save_settings_group(
        "consensus",
        {
            "threshold": opts.threshold,
            "ambiguous": normalize_consensus_ambiguous(opts.ambiguous),
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
        show_process_window=str(v.get("show_process_window", "false")).lower()
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
        },
    )


def load_te_annotation_options() -> TeAnnotationJobOptions:
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
    return TeAnnotationJobOptions(
        auto_modes=auto_modes,
        similarity_database_names=similarity_database_names,
        right_click_action=right_click_action,
    )


def save_te_annotation_options(opts: TeAnnotationJobOptions) -> None:
    from PySide6.QtCore import QSettings

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
    model: str = TREE_MODEL_FINDER
    runs: int = 5
    seed: int = 12345
    threads: int = field(default_factory=_default_tool_threads)
    # Support and rooting
    ufboot: int = 1000
    alrt: int = 1000
    bnni: bool = True
    distance_bootstrap: int = 100      # FastME bootstrap replicates
    root_mode: str = TREE_ROOT_NONE
    outgroup: list[str] = field(default_factory=list)
    # Options tab
    right_click_action: str = TREE_RIGHT_CLICK_RUN
    show_process_window: bool = False
    use_system_mafft: bool = False
    use_system_path: bool = False      # tree programs: system PATH only
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
        engine=_pick("engine", TREE_ENGINES, d.engine),
        model=str(v.get("model", d.model) or d.model).strip(),
        runs=_num("runs", d.runs, int, 1, 100),
        seed=_num("seed", d.seed, int, 1, 2_000_000_000),
        threads=_num("threads", _default_tool_threads(), int, 1, 256),
        ufboot=_num("ufboot", d.ufboot, int, 0, 100_000),
        alrt=_num("alrt", d.alrt, int, 0, 100_000),
        bnni=_to_bool(v.get("bnni"), d.bnni),
        distance_bootstrap=_num("distance_bootstrap", d.distance_bootstrap, int, 0, 10_000),
        root_mode=_pick("root_mode", TREE_ROOT_MODES, d.root_mode),
        outgroup=[str(x) for x in outgroup],
        right_click_action=_pick(
            "right_click_action", TREE_RIGHT_CLICK_ACTIONS, d.right_click_action
        ),
        show_process_window=_to_bool(v.get("show_process_window"), d.show_process_window),
        use_system_mafft=_to_bool(v.get("use_system_mafft"), d.use_system_mafft),
        use_system_path=_to_bool(v.get("use_system_path"), d.use_system_path),
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
            "model": opts.model,
            "runs": opts.runs,
            "seed": opts.seed,
            "threads": opts.threads,
            "ufboot": opts.ufboot,
            "alrt": opts.alrt,
            "bnni": opts.bnni,
            "distance_bootstrap": opts.distance_bootstrap,
            "root_mode": opts.root_mode,
            "outgroup": list(opts.outgroup),
            "right_click_action": opts.right_click_action,
            "show_process_window": opts.show_process_window,
            "use_system_mafft": opts.use_system_mafft,
            "use_system_path": opts.use_system_path,
        },
    )
