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
        from teviewer.backends.thread_budget import load_batch_workers

        return load_batch_workers()
    except Exception:  # noqa: BLE001
        return 4


DEFAULT_EXTEND_LEFT = 1000
DEFAULT_EXTEND_RIGHT = 1000
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


@dataclass
class TreeOptions:
    """Phylogenetic tree building.

    ``engine`` "auto" prefers IQ-TREE/FastTree when installed and falls back to
    the built-in Biopython distance method. ``model`` applies only to the
    distance method; the maximum-likelihood engines choose their own model.
    """

    engine: str = "auto"
    model: str = "jc69"
    method: str = "nj"
    outgroup: str = ""
    align_first: bool = True
    render_image: bool = True


@dataclass
class FastGaOptions:
    """FastGA whole-genome alignment.

    One selected file is compared against itself (finds repetitive regions);
    two files are compared against each other.
    """

    output_format: str = "paf"
    threads: int = field(default_factory=_default_tool_threads)
    min_identity: float = 0.7
    min_length: int = 100
    seed_freq: int = 200
    symmetric: bool = True
    keep_index: bool = True


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

    @property
    def active_database_name(self) -> str | None:
        if self.program == BLAST_PROGRAM_BLASTP:
            return self.blastp_database_name
        if self.program == BLAST_PROGRAM_BLASTX:
            return self.blastx_database_name
        return self.blast_database_name


DEFAULT_CONSENSUS_THRESHOLD = 0.8
DEFAULT_CONSENSUS_AMBIGUOUS = "N"


@dataclass
class ConsensusOptions:
    threshold: float = DEFAULT_CONSENSUS_THRESHOLD
    ambiguous: str = DEFAULT_CONSENSUS_AMBIGUOUS


DEFAULT_FLEXIDOT_WORDSIZE = 3
DEFAULT_FLEXIDOT_LCS_SHADING = True
DEFAULT_FLEXIDOT_USE_GFF = False
DEFAULT_FLEXIDOT_OUTPUT_FORMAT = "html"
DEFAULT_FLEXIDOT_SUBSTITUTIONS = 0
DEFAULT_FLEXIDOT_WOBBLE_CONVERSION = False
DEFAULT_FLEXIDOT_MAX_N_PERCENTAGE = 10.0
DEFAULT_FLEXIDOT_INCLUDE_REVERSE = True
DEFAULT_FLEXIDOT_ONLY_VS_FIRST = False
DEFAULT_FLEXIDOT_LINE_WIDTH = 1.0
DEFAULT_FLEXIDOT_LENGTH_SCALING = False
DEFAULT_FLEXIDOT_COMBINED = True
# Also run each selected file separately (in addition to / instead of combined).
DEFAULT_FLEXIDOT_RUN_SEPARATELY = False
# Flexidot comparison modes (-m): "0"=self, "1"=pairwise, "2"=all-vs-all.
# Empty tuple means "auto": self for a single file, pairwise/all-vs-all when
# combining (matching the historical behaviour).
FLEXIDOT_MODE_CHOICES = ("0", "1", "2")
DEFAULT_FLEXIDOT_MODES: tuple[str, ...] = ()


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
    """Return one of html / png / pdf / svg for Flexidot ``-f`` (default html)."""
    v = str(value).strip().lower()
    if v in FLEXIDOT_OUTPUT_FORMATS:
        return v
    return FLEXIDOT_OUTPUT_HTML


@dataclass
class FlexidotOptions:
    """Saved matching, display, and output parameters for Flexidot."""

    wordsize: int = DEFAULT_FLEXIDOT_WORDSIZE
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
    combined: bool = DEFAULT_FLEXIDOT_COMBINED
    # Also run each file separately; may be combined with ``combined`` (both).
    run_separately: bool = DEFAULT_FLEXIDOT_RUN_SEPARATELY
    # Explicit flexidot comparison modes (-m); empty tuple = auto (see above).
    modes: tuple[str, ...] = DEFAULT_FLEXIDOT_MODES


DEFAULT_TEAID_EVALUE = DEFAULT_REPORT_EVALUE
DEFAULT_TEAID_THREADS = 1
DEFAULT_TEAID_SELF_BLAST_WORD_SIZE = 11
MIN_TEAID_SELF_BLAST_WORD_SIZE = 4
MAX_TEAID_SELF_BLAST_WORD_SIZE = 25


def clamp_teaid_self_blast_word_size(requested: int) -> int:
    """Clamp TEAid self-BLAST ``-word_size`` to 4–25."""
    return max(
        MIN_TEAID_SELF_BLAST_WORD_SIZE,
        min(MAX_TEAID_SELF_BLAST_WORD_SIZE, int(requested)),
    )


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

MSA_CLEAN_METHOD_DESCRIPTIONS: dict[str, str] = {
    "msa_clean_gappy": "Remove columns that are mostly gaps (uses the gappy-column parameters below).",
    "msa_clean_div": (
        "Trim alignment ends where rows diverge from the consensus "
        "(uses the clean row with divergence parameters below)."
    ),
    "msa_clean_gap": (
        "Trim alignment ends where rows are mostly gaps "
        "(uses the clean row with gap parameters below)."
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

    def __post_init__(self) -> None:
        valid = [method for method in self.method_order if method in DEFAULT_MSA_CLEAN_ORDER]
        valid.extend(method for method in DEFAULT_MSA_CLEAN_ORDER if method not in valid)
        self.method_order = tuple(valid)
        if not (self.use_divergence or self.use_gap_row or self.use_gappy_columns):
            self.use_gappy_columns = True


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
        },
    )


def load_extend_options() -> ExtendOptions:
    v = _load_settings_group("extend")
    pan_name = v.get("pan_genome_name")
    return ExtendOptions(
        left=int(v.get("left", DEFAULT_EXTEND_LEFT)),
        right=int(v.get("right", DEFAULT_EXTEND_RIGHT)),
        pan_genome_name=str(pan_name).strip() if pan_name else None,
    )


def save_extend_options(opts: ExtendOptions) -> None:
    _save_settings_group(
        "extend",
        {
            "left": opts.left,
            "right": opts.right,
            "pan_genome_name": opts.pan_genome_name or "",
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
    return TeaAidOptions(
        e_value=normalize_blast_evalue(
            float(v.get("e_value", DEFAULT_TEAID_EVALUE)),
            default=DEFAULT_TEAID_EVALUE,
        ),
        self_blast_word_size=word_size,
        num_threads=max(1, int(v.get("num_threads", DEFAULT_TEAID_THREADS))),
    )


def save_teaid_options(opts: TeaAidOptions) -> None:
    _save_settings_group(
        "teaid",
        {
            "e_value": opts.e_value,
            "self_blast_word_size": opts.self_blast_word_size,
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


def load_consensus_options() -> ConsensusOptions:
    v = _load_settings_group("consensus")
    amb = str(v.get("ambiguous", DEFAULT_CONSENSUS_AMBIGUOUS)).upper()[:1] or "N"
    return ConsensusOptions(
        threshold=float(v.get("threshold", DEFAULT_CONSENSUS_THRESHOLD)),
        ambiguous=amb,
    )


def save_consensus_options(opts: ConsensusOptions) -> None:
    _save_settings_group(
        "consensus",
        {"threshold": opts.threshold, "ambiguous": opts.ambiguous},
    )


def load_flexidot_options() -> FlexidotOptions:
    v = _load_settings_group("flexidot")
    wordsize = int(v.get("wordsize", DEFAULT_FLEXIDOT_WORDSIZE))
    wordsize = max(MIN_FLEXIDOT_WORDSIZE, min(MAX_FLEXIDOT_WORDSIZE, wordsize))
    return FlexidotOptions(
        wordsize=wordsize,
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
        run_separately=str(v.get("run_separately", "false")).lower()
        in {"1", "true", "yes"},
        modes=parse_flexidot_modes(v.get("modes", "")),
    )


def save_flexidot_options(opts: FlexidotOptions) -> None:
    _save_settings_group(
        "flexidot",
        {
            "wordsize": opts.wordsize,
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
            "run_separately": opts.run_separately,
            "modes": ",".join(opts.modes),
        },
    )


def load_te_confusion_matrix_benchmark_options() -> TeConfusionMatrixBenchmarkOptions:
    v = _load_settings_group("te_confusion_benchmark")
    home = str(v.get("tetrimmer_home", "")).strip()
    if not home:
        from teviewer.backends.te_benchmark import default_tetrimmer_home

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
