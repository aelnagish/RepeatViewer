import logging
import os
import platform
import shutil
import tempfile
import traceback
import warnings
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from Bio import BiopythonDeprecationWarning, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from plotly.subplots import make_subplots

from repeatviewer.backends.gui_functions import (
    _blast_path_has_unsafe_chars,
    blast,
    blast_database_basename,
    check_database,
)
from repeatviewer.job_options import (
    DEFAULT_TE_DIVERGENCE_BIN_WIDTH,
    DEFAULT_TEAID_SELF_BLAST_WORD_SIZE,
    MAFFT_SEQUENCE_AUTO,
    clamp_teaid_self_blast_word_size,
    normalize_consensus_ambiguous,
    normalize_consensus_sequence_type,
)

# Suppress all deprecation warnings
warnings.filterwarnings('ignore', category=BiopythonDeprecationWarning)

"""
The BLAST, coverate, and self-dot plot codes are derived from TE-Aid. https://github.com/clemgoub/TE-Aid.git

Developer: Jiangzhao Qian
"""

os_type = platform.system()

TEAID_DIVERGENCE_HIST_MAX = 0.6
DEFAULT_TEAID_X_AXIS_LABEL = "Input sequence (bp)"
TEAID_PAN_PANEL_PX_MIN = 260
TEAID_PAN_PANEL_PX_MAX = 380
TEAID_PAN_X_TICKANGLE = -45


class TeaidInputError(ValueError):
    """Input FASTA is not suitable for TE-Aid consensus generation."""


def _pan_blast_subplot_title(genome_name: str, blast_plot_title: str) -> str:
    """Two-line pan-genome BLAST panel title: genome + fragments, then full length."""
    marker = '; full length: '
    if marker in blast_plot_title:
        fragments_part, full_length_part = blast_plot_title.split(marker, 1)
        return f'{genome_name}: {fragments_part}<br>full length: {full_length_part}'
    return f'{genome_name}: {blast_plot_title}'


def _pan_genome_panel_px(num_genomes: int) -> int:
    """Per-column panel width: never below the 5-genome baseline; expand when n grows."""
    if num_genomes <= 0:
        return TEAID_PAN_PANEL_PX_MIN
    if num_genomes <= 5:
        return max(TEAID_PAN_PANEL_PX_MIN, min(TEAID_PAN_PANEL_PX_MAX, int(1200 / num_genomes)))
    return TEAID_PAN_PANEL_PX_MIN


def _resolve_consensus_sequence_type(records, sequence_type: str) -> str:
    """Resolve auto-detect to nucleotide or protein from alignment content."""
    resolved = normalize_consensus_sequence_type(sequence_type)
    if resolved != MAFFT_SEQUENCE_AUTO:
        return resolved
    from repeatviewer.backends.msa_plot import detect_sequence_type

    seqs = [list(str(rec.seq).upper()) for rec in records]
    return detect_sequence_type(np.asarray(seqs, dtype="<U1"))


def _alignment_dumb_consensus(
    sequences: list[str],
    *,
    threshold: float = 0.8,
    ambiguous: str = "N",
) -> str:
    """Majority consensus matching Biopython's removed ``SummaryInfo.dumb_consensus``.

    Per column: among non-gap residues, if a single letter's share of non-gap
    letters meets ``threshold``, emit that letter; otherwise emit ``ambiguous``.
    All-gap columns use ``ambiguous``. Ties use ``ambiguous``.
    """
    if not sequences:
        return ""
    width = len(sequences[0])
    if any(len(seq) != width for seq in sequences):
        raise ValueError("MSA rows must share one length")
    try:
        thresh = float(threshold)
    except (TypeError, ValueError):
        thresh = 0.8
    amb = str(ambiguous or "N")[:1] or "N"
    out: list[str] = []
    for col in range(width):
        counts: dict[str, int] = {}
        num_atoms = 0
        for seq in sequences:
            ch = seq[col]
            if ch in "-.":
                continue
            up = ch.upper()
            counts[up] = counts.get(up, 0) + 1
            num_atoms += 1
        if num_atoms <= 0:
            out.append(amb)
            continue
        max_size = max(counts.values())
        winners = [letter for letter, n in counts.items() if n == max_size]
        if len(winners) == 1 and max_size / num_atoms >= thresh:
            out.append(winners[0])
        else:
            out.append(amb)
    return "".join(out)


def con_generater(
    input_file,
    output_dir,
    threshold=0.8,
    ambiguous="N",
    sequence_type=MAFFT_SEQUENCE_AUTO,
):
    """
    Build a consensus sequence for TE-Aid and the Consensus tool.

    Accepts a single-sequence FASTA or a multiple sequence alignment (all rows
    the same length). Unaligned multi-FASTA files are rejected with a clear error.
    """
    records = list(SeqIO.parse(input_file, "fasta"))
    if not records:
        raise TeaidInputError(
            f"No sequences found in {os.path.basename(input_file)}."
        )

    amb = normalize_consensus_ambiguous(ambiguous)
    # Normalize for API compatibility; majority consensus needs no alphabet.
    normalize_consensus_sequence_type(sequence_type)

    if len(records) == 1:
        consensus_seq = Seq(str(records[0].seq).replace("-", "")).upper()
    else:
        lengths = {len(rec.seq) for rec in records}
        if len(lengths) > 1:
            min_len = min(lengths)
            max_len = max(lengths)
            raise TeaidInputError(
                f"{os.path.basename(input_file)} contains {len(records)} sequences "
                f"with different lengths ({min_len}–{max_len}). "
                "Provide a multiple sequence alignment (all rows the same length). "
                "Align the sequences first (Align menu) or use a single-sequence FASTA."
            )

        # Own majority consensus: Biopython removed SummaryInfo.dumb_consensus.
        consensus_seq = Seq(
            _alignment_dumb_consensus(
                [str(rec.seq) for rec in records],
                threshold=threshold,
                ambiguous=amb,
            )
        )

    consensus_seq_len = len(consensus_seq)
    consensus_seq_obj = SeqRecord(
        consensus_seq,
        id="consensus",
        description=f"threshold={threshold};amb={amb}",
    )

    output_file = os.path.join(output_dir, f'{os.path.basename(input_file)}_{threshold}_co.fa')
    with open(output_file, 'w') as file:
        SeqIO.write(consensus_seq_obj, file, 'fasta')

    return output_file, consensus_seq_len


def _nice_consensus_tick_step(max_len: int, *, target_ticks: int = 8) -> int:
    """Choose a 1/2/5×10ⁿ step so the axis has about ``target_ticks`` labels."""
    max_len = max(1, int(max_len))
    target = max(3, int(target_ticks))
    raw = max(1.0, float(max_len) / float(target))
    exp = int(np.floor(np.log10(raw)))
    base = 10.0 ** exp
    frac = raw / base
    if frac <= 1.0:
        nice = 1.0
    elif frac <= 2.0:
        nice = 2.0
    elif frac <= 5.0:
        nice = 5.0
    else:
        nice = 10.0
    return max(1, int(nice * base))


def _format_consensus_tick_label(pos: int) -> str:
    """Compact axis label (``15k``) so Plotly SI prefixes do not crowd the axis."""
    value = max(0, int(pos))
    if value >= 10_000 and value % 1000 == 0:
        return f"{value // 1000}k"
    if value >= 10_000:
        text = f"{value / 1000:.1f}k"
        return text.replace(".0k", "k")
    return str(value)


def _consensus_position_ticks(max_len: int) -> list[int]:
    """1-based TE consensus tick values; always includes 1, never the exact sequence end.

    Step grows with length (nice 1/2/5×10ⁿ) so long consensus axes stay readable
    instead of packing a tick every 1 kb on a 20–30 kb TE.
    """
    max_len = max(1, int(max_len))
    step = _nice_consensus_tick_step(max_len)
    ticks = [1]
    pos = step
    while pos < max_len:
        if pos != 1:
            ticks.append(pos)
        pos += step
    return ticks


def _consensus_axis_tick_kwargs(max_len: int, *, tickangle: int = 0) -> dict:
    ticks = _consensus_position_ticks(max_len)
    return {
        'tickmode': 'array',
        'tickvals': ticks,
        'ticktext': [_format_consensus_tick_label(t) for t in ticks],
        'tickangle': tickangle,
    }


#####################################################################################################
# Code block: genome blast
#####################################################################################################
def teaid_genome_blast(input_file, genome_file, output_dir, e_value, num_threads=1,
                       use_system_blast=False, output_label=None, cancel=None):
    # input_file: fasta file contain only one sequence
    # Used global variable: os_type

    if genome_file is None:
        logging.error("WARNING: Genome file not found. TE-Aid can't perform genome blast.")
        return False

    # Check if provided genome file exist
    elif not os.path.isfile(genome_file):
        logging.error(
            "WARNING: You provided genome file not found. TE-Aid can't perform genome blast."
        )
        return False

    # Check genome blast database availability
    blast_database = check_database(
        genome_file, os_type=os_type, use_system_blast=use_system_blast, cancel=cancel
    )

    # Check if makeblastdb is correctly installed
    if blast_database == 'makeblastdb_not_found':
        logging.error(
            "Error: makeblastdb command not found. TE-Aid can't perform genome blast."
        )
        return False

    # Check if error happened
    elif blast_database == 'makeblastdb_got_error':
        logging.error(
            "Error: BLAST database can't be established for the genome. TE-Aid can't perform genome blast."
        )
        return False

    # Run blast command
    genome_blast_out, _ = blast(
        input_file,
        blast_database,
        output_dir,
        e_value=e_value,
        os_type=os_type,
        num_threads=num_threads,
        use_system_blast=use_system_blast,
        output_label=output_label,
        cancel=cancel,
    )

    if genome_blast_out == 'blastn_not_found':
        logging.error("Error: BLAST command not found. TE-Aid can't perform genome blast.")
        return False

    elif genome_blast_out == 'blastn_got_error':
        logging.error(
            "TE-Aid can't perform genome blast. Refer to terminal for more information."
        )
        return False

    elif genome_blast_out == 'blast_n_zero':
        logging.warning('Genome blast hit number is 0 for this sequence.')
        return 'blast_n_zero'

    elif not genome_blast_out:
        logging.error(
            "TE-Aid can't perform genome blast. Refer to terminal for more information."
        )
        return False

    return genome_blast_out


#####################################################################################################
# Code block: self blast
#####################################################################################################
def self_blast(input_file, output_dir, use_system_blast=False, word_size=11, cancel=None):
    # input_file: fasta file contain only one sequence
    # Used global variable: os_type

    # NCBI BLAST mishandles paths with spaces; keep the work tree space-free.
    work_parent = output_dir
    if _blast_path_has_unsafe_chars(os.path.abspath(str(output_dir))):
        work_parent = tempfile.gettempdir()
    work_dir = tempfile.mkdtemp(prefix=".repeatviewer_selfblast_", dir=work_parent)

    query_path = os.path.abspath(str(input_file))
    if _blast_path_has_unsafe_chars(query_path):
        safe_query = os.path.join(work_dir, blast_database_basename(query_path))
        try:
            os.link(query_path, safe_query)
        except OSError:
            try:
                os.symlink(query_path, safe_query)
            except OSError:
                shutil.copy2(query_path, safe_query)
        query_path = safe_query

    # Build blast database
    blast_database = check_database(
        query_path,
        work_dir,
        os_type=os_type,
        use_system_blast=use_system_blast,
        cancel=cancel,
    )

    # Check if makeblastdb is correctly installed
    if blast_database == 'makeblastdb_not_found':
        logging.error(
            "Error: makeblastdb command not found. TE-Aid can't perform self blast."
        )
        return False, work_dir

    # Check if error happened
    elif blast_database == 'makeblastdb_got_error':
        logging.error(
            "Error: BLAST database can't be established for the genome. TE-Aid can't perform self blast."
        )
        return False, work_dir

    self_blast_out, _ = blast(
        query_path,
        blast_database,
        work_dir,
        e_value=0.05,
        self_blast=True,
        word_size=word_size,
        os_type=os_type,
        use_system_blast=use_system_blast,
        cancel=cancel,
    )

    if self_blast_out == 'blastn_not_found':
        logging.error("Error: BLAST command not found. TE-Aid can't perform self blast.")
        return False, work_dir

    elif self_blast_out == 'blastn_got_error':
        logging.error(
            "TE-Aid can't perform self blast. Refer to terminal for more information."
        )
        return False, work_dir

    elif self_blast_out == 'blast_n_zero':
        logging.warning('Self blast hit number is 0 for this sequence.')
        return 'blast_n_zero', work_dir

    elif not self_blast_out:
        logging.error(
            "TE-Aid can't perform self blast. Refer to terminal for more information."
        )
        return False, work_dir

    return self_blast_out, work_dir


#####################################################################################################
# Code block: Empty plot
#####################################################################################################
# Define function to plot an empty plot, which will be used when blast or rpstblastn result is zero or get an error.
def empty_plot(seq_len, width_n=600, height_n=600, custom_text=None):
    # Initialize an empty figure
    fig = go.Figure()

    # Add optional custom text at the center of the plot
    if custom_text:
        fig.add_trace(
            go.Scatter(
                x=[int(seq_len) / 2],  # Centered in the middle of the x-axis range
                y=[5],  # Centered in the middle of the y-axis range
                text=[custom_text],
                mode='text',
                textfont={'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
                textposition='middle center',
                showlegend=False,
                hoverinfo='none',
            )
        )

    # Update layout to match the style of other plots
    fig.update_layout(
        xaxis={
            'title': '',  # No title for empty plot
            'showgrid': False,
            'zeroline': False,
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
            'range': [1, int(seq_len)],  # BLAST positions are 1-based
            **_consensus_axis_tick_kwargs(int(seq_len)),
            'ticks': '',  # Hide tick marks
            'showticklabels': True,  # Hide tick labels
        },
        yaxis={
            'title': '',  # No title for empty plot
            'showgrid': False,
            'zeroline': False,
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
            'range': [0, 10],
            'ticks': '',
            'showticklabels': True,
        },
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin={'l': 60, 'r': 60, 't': 60, 'b': 60},
        font={'family': 'Arial, sans-serif', 'size': 12, 'color': 'black'},
        width=int(width_n),
        height=int(height_n),
    )

    return fig


#####################################################################################################
# Code block: genome blast plot
#####################################################################################################
# Function to create hover text
def create_hovertext(row):
    divergence = 100 - row['V3']
    return f'Start: {row["V7"]}<br>End: {row["V8"]}<br>Divergence: {divergence:.2f}%'


def _plotly_rgba(color, opacity=1.0) -> str:
    """Convert a color name/hex to ``rgba(r,g,b,a)`` for Plotly line opacity."""
    opacity = max(0.0, min(1.0, float(opacity)))
    text = str(color or "").strip()
    if not text:
        return f"rgba(0, 0, 0, {opacity})"
    if text.startswith("rgba(") or text.startswith("rgb("):
        return text
    named = {
        "red": (255, 0, 0),
        "blue": (0, 0, 255),
        "green": (0, 128, 0),
        "black": (0, 0, 0),
        "white": (255, 255, 255),
        "orange": (255, 165, 0),
        "purple": (128, 0, 128),
        "cyan": (0, 255, 255),
        "magenta": (255, 0, 255),
        "yellow": (255, 255, 0),
        "gray": (128, 128, 128),
        "grey": (128, 128, 128),
    }
    lower = text.lower()
    if lower in named:
        r, g, b = named[lower]
        return f"rgba({r}, {g}, {b}, {opacity})"
    if text.startswith("#"):
        hex_body = text[1:]
        if len(hex_body) == 3:
            hex_body = "".join(ch * 2 for ch in hex_body)
        if len(hex_body) >= 6:
            try:
                r = int(hex_body[0:2], 16)
                g = int(hex_body[2:4], 16)
                b = int(hex_body[4:6], 16)
                return f"rgba({r}, {g}, {b}, {opacity})"
            except ValueError:
                pass
    try:
        from matplotlib.colors import to_rgb

        r, g, b = to_rgb(text)
        return (
            f"rgba({int(round(r * 255))}, {int(round(g * 255))}, "
            f"{int(round(b * 255))}, {opacity})"
        )
    except Exception:
        return text


def blast_plot(
    df,
    cons_len,
    full_len_thr=0.8,
    mid_len_thr=0.5,
    full_color='red',
    mid_color='#1f77b4',
    line_transparency=0.0,
):
    # Length bands as fractions of consensus / input length.
    hit_len = (df['V7'] - df['V8']).abs()
    full = df[hit_len >= full_len_thr * cons_len]
    mid = df[
        (hit_len >= mid_len_thr * cons_len) & (hit_len < full_len_thr * cons_len)
    ]
    short = df[hit_len < mid_len_thr * cons_len]

    # Calculate the max number of divergence for the y-axis range
    y_max = max(100 - df['V3'])

    # Create a figure
    fig = go.Figure()
    # Transparency 0 = opaque, 1 = invisible. Apply to every hit line color.
    opacity = max(0.0, min(1.0, 1.0 - float(line_transparency)))
    short_line_color = _plotly_rgba('black', opacity)
    mid_line_color = _plotly_rgba(mid_color, opacity)
    full_line_color = _plotly_rgba(full_color, opacity)

    # Short hits (below mid threshold) — black only, no opaque underlay.
    for i in short.index:
        fig.add_trace(
            go.Scatter(
                x=[short.loc[i, 'V7'], short.loc[i, 'V8']],
                y=[100 - short.loc[i, 'V3'], 100 - short.loc[i, 'V3']],
                mode='lines',
                line={'color': short_line_color, 'width': 1},
                hoverinfo='text',
                text=create_hovertext(short.loc[i]),
                showlegend=False,
            )
        )

    # Mid-length hits (between mid and full thresholds)
    for i in mid.index:
        fig.add_trace(
            go.Scatter(
                x=[mid.loc[i, 'V7'], mid.loc[i, 'V8']],
                y=[100 - mid.loc[i, 'V3'], 100 - mid.loc[i, 'V3']],
                mode='lines',
                line={'color': mid_line_color, 'width': 2},
                hoverinfo='text',
                text=create_hovertext(mid.loc[i]),
                showlegend=False,
            )
        )

    # Full-length hits (at or above the upper threshold)
    for i in full.index:
        fig.add_trace(
            go.Scatter(
                x=[full.loc[i, 'V7'], full.loc[i, 'V8']],
                y=[100 - full.loc[i, 'V3'], 100 - full.loc[i, 'V3']],
                mode='lines',
                line={'color': full_line_color, 'width': 3},
                hoverinfo='text',
                text=create_hovertext(full.loc[i]),
                showlegend=False,
            )
        )

    # Define plot title
    full_bp = int(full_len_thr * cons_len)
    mid_bp = int(mid_len_thr * cons_len)
    plot_title = (
        f'Fragments number: {len(df)}; '
        f'full length: {len(full)} (>={full_bp}bp); '
        f'mid length: {len(mid)} (>={mid_bp}bp)'
    )

    # Update the layout to add titles and labels with adjusted axis ranges
    fig.update_layout(
        title=plot_title,
        xaxis_title=DEFAULT_TEAID_X_AXIS_LABEL,
        yaxis_title='Divergence to consensus (%)',
        xaxis={
            'title_font': {'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
            'showgrid': False,
            'range': [1, cons_len],
            **_consensus_axis_tick_kwargs(cons_len),
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
        },
        yaxis={
            'showgrid': False,
            'title_font': {'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
            'range': [-(y_max / 30), y_max + (y_max / 30)],
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
        },
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin={
            'l': 60, 'r': 60, 't': 60, 'b': 60
        },  # Adjust margins to ensure frame visibility
        font={'family': 'Arial, sans-serif', 'size': 12, 'color': 'black'},
        width=600,
        height=600,
    )
    # fig.show()
    return fig, y_max, plot_title


#####################################################################################################
# Code block: coverage plot
#####################################################################################################
def coverage_plot(df, cons_len):
    # Initialize coverage array
    coverage = np.zeros((len(df), cons_len))

    # Populate coverage array
    for i in range(len(df)):
        start = df.iloc[i]['V7'] - 1  # Convert to zero-indexed
        end = df.iloc[i]['V8']
        coverage[i, start:end] = 1  # Set covered positions to 1

    # Sum coverage over all rows to get total coverage at each position
    total_coverage = np.sum(coverage, axis=0)
    max_coverage = np.max(total_coverage)

    # Create a figure for plotting
    fig = go.Figure()

    # Add the coverage line
    fig.add_trace(
        go.Scatter(
            x=np.arange(1, cons_len + 1),  # X-axis from 1 to cons_len
            y=total_coverage,
            mode='lines',
            line={'color': 'black', 'width': 3},
            showlegend=False,
        )
    )

    # Update layout to add titles and labels
    fig.update_layout(
        title='TE Consensus Genomic Coverage Plot',
        xaxis_title=DEFAULT_TEAID_X_AXIS_LABEL,
        yaxis_title='Coverage (bp)',
        xaxis={
            'title_font': {'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
            'showgrid': False,
            'range': [1, cons_len],
            **_consensus_axis_tick_kwargs(cons_len),
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
        },
        yaxis={
            'showgrid': False,
            'title_font': {'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
            'range': [-(max_coverage / 30), max_coverage + (max_coverage / 30)],
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
        },
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin={
            'l': 60, 'r': 60, 't': 60, 'b': 60
        },  # Adjust margins to ensure frame visibility
        font={'family': 'Arial, sans-serif', 'size': 12, 'color': 'black'},
        width=600,
        height=600,
    )
    # fig.show()
    return fig, max_coverage


#####################################################################################################
# Code block: divergence histogram (TE divergence plot logic, single TE)
#####################################################################################################
def empty_divergence_histogram(custom_text=None):
    fig = go.Figure()
    y_max = 10
    y_center = (y_max * 1.08) / 2
    if custom_text:
        fig.add_trace(
            go.Scatter(
                x=[TEAID_DIVERGENCE_HIST_MAX / 2],
                y=[y_center],
                text=[custom_text],
                mode='text',
                textfont={'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
                textposition='middle center',
                showlegend=False,
                hoverinfo='none',
            )
        )
    fig.update_layout(
        title='Divergence Plot',
        xaxis_title='Divergence',
        yaxis_title='Total length (bp)',
        xaxis={
            'title_font': {'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
            'showgrid': True,
            'gridcolor': '#f0f0f0',
            'range': [0, TEAID_DIVERGENCE_HIST_MAX],
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
        },
        yaxis={
            'title_font': {'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
            'showgrid': True,
            'gridcolor': '#f4f4f4',
            'range': [0, y_max],
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
        },
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin={'l': 60, 'r': 60, 't': 60, 'b': 60},
        font={'family': 'Arial, sans-serif', 'size': 12, 'color': 'black'},
        width=600,
        height=600,
    )
    return fig, y_max


def divergence_histogram_plot(
    df,
    *,
    bin_width: float = DEFAULT_TE_DIVERGENCE_BIN_WIDTH,
):
    """
    Histogram of genome BLAST hit divergence (1 - identity), weighted by alignment length.

    Matches ``write_blast_divergence_plot_html`` in ``blast_genome_plot.py``.
    """
    bw = max(0.001, float(bin_width))
    div_max = TEAID_DIVERGENCE_HIST_MAX
    n_bins = int(np.ceil(div_max / bw))
    edges = np.linspace(0.0, div_max, n_bins + 1)
    centers = (edges[:-1] + edges[1:]) / 2.0

    div = ((100.0 - df['V3'].astype(float)) / 100.0).clip(lower=0.0, upper=div_max)
    bp = df['V4'].astype(float).clip(lower=0.0)
    bp_per_bin, _ = np.histogram(
        div.to_numpy(dtype=float),
        bins=edges,
        weights=bp.to_numpy(dtype=float),
    )
    active = bp_per_bin > 0
    if not np.any(active):
        return empty_divergence_histogram('No aligned length in divergence bins')

    active_idx = np.where(active)[0]
    plot_centers = centers[active]
    plot_bp = bp_per_bin[active]
    plot_edges_lo = edges[active_idx]
    plot_edges_hi = edges[active_idx + 1]
    y_max = max(1.0, float(plot_bp.max()))

    hover = [
        (
            f'Divergence: {lo:.3f} - {hi:.3f}<br>'
            f'Total length: {int(val):,} bp'
        )
        for lo, hi, val in zip(plot_edges_lo, plot_edges_hi, plot_bp)
    ]

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=plot_centers,
            y=plot_bp,
            width=[bw] * len(plot_centers),
            marker={'color': 'black', 'line': {'width': 0}},
            hovertext=hover,
            hoverinfo='text',
            showlegend=False,
        )
    )
    fig.update_layout(
        title='Divergence Plot',
        xaxis_title='Divergence',
        yaxis_title='Total length (bp)',
        xaxis={
            'title_font': {'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
            'showgrid': True,
            'gridcolor': '#f0f0f0',
            'range': [0, TEAID_DIVERGENCE_HIST_MAX],
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
        },
        yaxis={
            'title_font': {'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
            'showgrid': True,
            'gridcolor': '#f4f4f4',
            'range': [0, y_max * 1.08],
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
        },
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin={'l': 60, 'r': 60, 't': 60, 'b': 60},
        font={'family': 'Arial, sans-serif', 'size': 12, 'color': 'black'},
        width=600,
        height=600,
        bargap=0,
    )
    return fig, y_max


#####################################################################################################
# Code block: self dot plot
#####################################################################################################
def dot_plot(df):
    fig = go.Figure()

    # Loop through data and add segments
    for i in range(len(df)):
        segment = df.iloc[i]
        if segment['V10'] > segment['V9']:
            color = 'black'
        else:
            color = '#009E73'  # Use different color for TIR
        fig.add_trace(
            go.Scatter(
                x=[segment['V7'], segment['V8']],
                y=[segment['V9'], segment['V10']],
                mode='lines',
                line={'color': color, 'width': 3},
                hoverinfo='text',
                text=f'x left: {segment["V7"]}, y left: {segment["V9"]}<br>x right: {segment["V8"]}, y right: {segment["V10"]}',
                showlegend=False,
            )
        )

    axis_max = int(
        max(
            df['V7'].max(),
            df['V8'].max(),
            df['V9'].max(),
            df['V10'].max(),
        )
    )
    axis_ticks = _consensus_axis_tick_kwargs(axis_max)
    fig.update_layout(
        title='Self dot plot',
        xaxis_title=DEFAULT_TEAID_X_AXIS_LABEL,
        yaxis_title=DEFAULT_TEAID_X_AXIS_LABEL,
        xaxis={
            'title_font': {'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
            'showgrid': False,
            'range': [1, axis_max],
            **axis_ticks,
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
        },
        yaxis={
            'showgrid': False,
            'title_font': {'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
            'range': [1, axis_max],
            **axis_ticks,
            'linewidth': 1,
            'linecolor': 'black',
            'mirror': True,
        },
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin={
            'l': 60, 'r': 60, 't': 60, 'b': 60
        },  # Adjust margins to ensure frame visibility
        font={'family': 'Arial, sans-serif', 'size': 12, 'color': 'black'},
        width=600,
        height=600,
    )
    # fig.show()
    return fig


def _self_dotplot_axis_max(df: pd.DataFrame, *, fallback: int = 1) -> int:
    if df.empty:
        return max(1, int(fallback))
    return int(
        max(
            df["V7"].max(),
            df["V8"].max(),
            df["V9"].max(),
            df["V10"].max(),
            1,
        )
    )


def render_self_dotplot_png(
    df: pd.DataFrame,
    output_png: os.PathLike[str] | str,
    *,
    title: str = "Self dot plot",
) -> str:
    """Render TE-Aid self-BLAST hits to a static PNG (no WebEngine)."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_path = os.path.abspath(os.path.expanduser(str(output_png)))
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    axis_max = _self_dotplot_axis_max(df)

    fig, ax = plt.subplots(figsize=(6.4, 6.4), dpi=120)
    if not df.empty:
        for _index, segment in df.iterrows():
            color = "#000000" if segment["V10"] > segment["V9"] else "#009E73"
            ax.plot(
                [segment["V7"], segment["V8"]],
                [segment["V9"], segment["V10"]],
                color=color,
                linewidth=1.2,
                solid_capstyle="round",
            )
    ax.set_xlim(1, axis_max)
    ax.set_ylim(1, axis_max)
    ax.set_aspect("equal", adjustable="box")
    ax.set_title(title)
    ax.set_xlabel(DEFAULT_TEAID_X_AXIS_LABEL)
    ax.set_ylabel(DEFAULT_TEAID_X_AXIS_LABEL)
    ax.tick_params(direction="out")
    for spine in ax.spines.values():
        spine.set_visible(True)
    fig.tight_layout()
    fig.savefig(output_path, facecolor="white")
    plt.close(fig)
    return output_path


def _render_self_dotplot_message_png(
    message: str,
    seq_len: int,
    output_png: os.PathLike[str] | str,
) -> str:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    output_path = os.path.abspath(os.path.expanduser(str(output_png)))
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    seq_len = max(1, int(seq_len))
    fig, ax = plt.subplots(figsize=(6.4, 6.4), dpi=120)
    ax.set_xlim(1, seq_len)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.text(
        seq_len / 2,
        5,
        message,
        ha="center",
        va="center",
        fontsize=12,
        wrap=True,
    )
    fig.tight_layout()
    fig.savefig(output_path, facecolor="white")
    plt.close(fig)
    return output_path


def _preview_self_dotplot_paths(
    input_file: os.PathLike[str] | str,
    output_png: os.PathLike[str] | str,
    *,
    word_size: int,
    use_system_blast: bool,
    work_dir: os.PathLike[str] | str | None,
) -> tuple[str, int, str, str]:
    input_path = os.path.abspath(os.path.expanduser(str(input_file)))
    output_path = os.path.abspath(os.path.expanduser(str(output_png)))
    if not os.path.isfile(input_path):
        raise FileNotFoundError(f"FASTA not found: {input_path}")

    records = list(SeqIO.parse(input_path, "fasta"))
    if not records:
        raise TeaidInputError(f"No sequences found in {os.path.basename(input_path)}.")
    cons_len = len(str(records[0].seq).replace("-", ""))
    if cons_len < 1:
        raise TeaidInputError(f"No sequence in {os.path.basename(input_path)}.")

    scratch_dir = (
        os.path.abspath(os.path.expanduser(str(work_dir)))
        if work_dir is not None
        else os.path.dirname(output_path)
    )
    os.makedirs(scratch_dir, exist_ok=True)
    return input_path, cons_len, output_path, scratch_dir


def generate_preview_self_dotplot_png(
    input_file: os.PathLike[str] | str,
    output_png: os.PathLike[str] | str,
    *,
    word_size: int = 11,
    use_system_blast: bool = False,
    work_dir: os.PathLike[str] | str | None = None,
) -> str:
    """Build a TE-Aid self-BLAST dotplot PNG (used by TE-Aid / export helpers)."""
    input_path, cons_len, output_path, scratch_dir = _preview_self_dotplot_paths(
        input_file,
        output_png,
        word_size=word_size,
        use_system_blast=use_system_blast,
        work_dir=work_dir,
    )
    dotplot_blast, blast_work_dir = self_blast(
        input_path,
        scratch_dir,
        use_system_blast=use_system_blast,
        word_size=word_size,
    )
    try:
        if dotplot_blast == "blast_n_zero":
            return _render_self_dotplot_message_png(
                "No self BLAST hits found for this sequence.",
                cons_len,
                output_path,
            )
        if dotplot_blast:
            df_dotplot_blast = pd.read_csv(dotplot_blast)
            return render_self_dotplot_png(df_dotplot_blast, output_path)
        return _render_self_dotplot_message_png(
            "An error occurred. Please check the terminal for details.",
            cons_len,
            output_path,
        )
    finally:
        shutil.rmtree(blast_work_dir, ignore_errors=True)


#####################################################################################################
# Code block: TEtrimmer GUI plotter, integrate all plots
#####################################################################################################
_AXIS_TITLE_FONT = {'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'}
_SUBPLOT_TITLE_FONT = {'size': 13, 'family': 'Arial, sans-serif', 'color': 'black'}


def _build_genome_teaid_panels(
    con_seq,
    cons_len,
    genome_file,
    output_dir,
    e_value,
    *,
    num_threads=5,
    use_system_blast=False,
    output_label=None,
    hit_len_high_pct=80.0,
    hit_len_mid_pct=50.0,
    hit_color_high='red',
    hit_color_mid='#1f77b4',
    hit_line_transparency=0.0,
    cancel=None,
):
    """BLAST hits, coverage, and divergence panels for one genome."""
    genome_blast_out = teaid_genome_blast(
        con_seq,
        genome_file,
        output_dir,
        e_value,
        num_threads=num_threads,
        use_system_blast=use_system_blast,
        output_label=output_label,
        cancel=cancel,
    )

    if genome_blast_out:
        if genome_blast_out == 'blast_n_zero':
            fig_blast = empty_plot(cons_len, custom_text='No genome BLAST hits found')
            divergence_max = 10
            blast_plot_title = 'BLAST hits Plot'
            fig_coverage = empty_plot(
                cons_len, custom_text='No genome BLAST hits found'
            )
            coverage_max = 10
            fig_div_hist, div_hist_max = empty_divergence_histogram(
                'No genome BLAST hits found'
            )
        else:
            df_genome_blast_out = pd.read_csv(genome_blast_out)
            full_thr = max(0.0, min(1.0, float(hit_len_high_pct) / 100.0))
            mid_thr = max(0.0, min(full_thr, float(hit_len_mid_pct) / 100.0))
            if mid_thr >= full_thr and full_thr > 0:
                mid_thr = max(0.0, full_thr - 0.01)
            fig_blast, divergence_max, blast_plot_title = blast_plot(
                df_genome_blast_out,
                cons_len,
                full_len_thr=full_thr,
                mid_len_thr=mid_thr,
                full_color=hit_color_high,
                mid_color=hit_color_mid,
                line_transparency=hit_line_transparency,
            )
            fig_coverage, coverage_max = coverage_plot(df_genome_blast_out, cons_len)
            fig_div_hist, div_hist_max = divergence_histogram_plot(df_genome_blast_out)
    else:
        error_text = (
            'An error occurred. Please check the terminal for details.<br> '
            'Please provide a genome file if it has not been specified.'
        )
        fig_blast = empty_plot(cons_len, custom_text=error_text)
        divergence_max = 10
        blast_plot_title = 'BLAST hit Plot'
        fig_coverage = empty_plot(cons_len, custom_text=error_text)
        coverage_max = 10
        fig_div_hist, div_hist_max = empty_divergence_histogram(
            'An error occurred. Please check the terminal for details.'
        )

    return {
        'fig_blast': fig_blast,
        'fig_coverage': fig_coverage,
        'fig_div_hist': fig_div_hist,
        'divergence_max': divergence_max,
        'coverage_max': coverage_max,
        'div_hist_max': div_hist_max,
        'blast_plot_title': blast_plot_title,
    }


def _build_self_dot_figure(
    con_seq,
    output_dir,
    cons_len,
    *,
    use_system_blast=False,
    self_blast_word_size=DEFAULT_TEAID_SELF_BLAST_WORD_SIZE,
    cancel=None,
):
    self_word_size = clamp_teaid_self_blast_word_size(self_blast_word_size)
    dotplot_blast, work_dir = self_blast(
        con_seq,
        output_dir,
        use_system_blast=use_system_blast,
        word_size=self_word_size,
        cancel=cancel,
    )
    try:
        if dotplot_blast == 'blast_n_zero':
            return empty_plot(
                cons_len,
                custom_text='No self BLAST hits found for this sequence.',
            )
        if dotplot_blast:
            df_dotplot_blast = pd.read_csv(dotplot_blast)
            return dot_plot(df_dotplot_blast)
        return empty_plot(
            cons_len,
            custom_text='An error occurred. Please check the terminal for details.',
        )
    finally:
        shutil.rmtree(work_dir, ignore_errors=True)


def _configure_blast_subplot(
    fig,
    row,
    col,
    cons_len,
    divergence_max,
    consensus_ticks,
    *,
    show_x_title=True,
    show_y_title=True,
    x_tickangle=0,
    x_axis_label=DEFAULT_TEAID_X_AXIS_LABEL,
):
    fig.update_xaxes(
        title_text=x_axis_label if show_x_title else '',
        showgrid=False,
        showticklabels=True,
        title_font=_AXIS_TITLE_FONT,
        range=[1, cons_len],
        **_consensus_axis_tick_kwargs(cons_len, tickangle=x_tickangle),
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=row,
        col=col,
    )
    fig.update_yaxes(
        title_text='Divergence to consensus (%)' if show_y_title else '',
        showgrid=False,
        showticklabels=True,
        title_font=_AXIS_TITLE_FONT,
        range=[-(divergence_max / 30), divergence_max + (divergence_max / 30)],
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=row,
        col=col,
    )


def _configure_coverage_subplot(
    fig,
    row,
    col,
    cons_len,
    coverage_max,
    consensus_ticks,
    *,
    show_x_title=True,
    show_y_title=True,
    x_tickangle=0,
    x_axis_label=DEFAULT_TEAID_X_AXIS_LABEL,
):
    fig.update_xaxes(
        title_text=x_axis_label if show_x_title else '',
        showgrid=False,
        showticklabels=True,
        title_font=_AXIS_TITLE_FONT,
        range=[1, cons_len],
        **_consensus_axis_tick_kwargs(cons_len, tickangle=x_tickangle),
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=row,
        col=col,
    )
    fig.update_yaxes(
        title_text='Coverage (bp)' if show_y_title else '',
        showgrid=False,
        showticklabels=True,
        title_font=_AXIS_TITLE_FONT,
        range=[-(coverage_max / 30), coverage_max + (coverage_max / 30)],
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=row,
        col=col,
    )


def _configure_self_dot_subplot(
    fig,
    row,
    col,
    cons_len,
    *,
    show_x_title=True,
    show_y_title=True,
    x_tickangle=0,
    x_axis_label=DEFAULT_TEAID_X_AXIS_LABEL,
):
    fig.update_xaxes(
        title_text=x_axis_label if show_x_title else '',
        showgrid=False,
        showticklabels=True,
        title_font=_AXIS_TITLE_FONT,
        range=[1, cons_len],
        **_consensus_axis_tick_kwargs(cons_len, tickangle=x_tickangle),
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=row,
        col=col,
    )
    fig.update_yaxes(
        title_text=x_axis_label if show_y_title else '',
        showgrid=False,
        showticklabels=True,
        title_font=_AXIS_TITLE_FONT,
        range=[1, cons_len],
        **_consensus_axis_tick_kwargs(cons_len, tickangle=0),
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=row,
        col=col,
    )


def _configure_divergence_subplot(
    fig,
    row,
    col,
    div_hist_max,
    *,
    show_x_title=True,
    show_y_title=True,
    x_tickangle=0,
):
    fig.update_xaxes(
        title_text='Divergence' if show_x_title else '',
        showgrid=False,
        showticklabels=True,
        title_font=_AXIS_TITLE_FONT,
        range=[0, TEAID_DIVERGENCE_HIST_MAX],
        tickangle=x_tickangle,
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=row,
        col=col,
    )
    fig.update_yaxes(
        title_text='Total length (bp)' if show_y_title else '',
        showgrid=False,
        showticklabels=True,
        title_font=_AXIS_TITLE_FONT,
        range=[0, div_hist_max * 1.08],
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=row,
        col=col,
    )


def _consensus_axis_center(cons_len):
    return (1 + int(cons_len)) / 2


def _divergence_axis_center(divergence_max):
    lo = -(divergence_max / 30)
    hi = divergence_max + (divergence_max / 30)
    return (lo + hi) / 2


def _coverage_axis_center(coverage_max):
    lo = -(coverage_max / 30)
    hi = coverage_max + (coverage_max / 30)
    return (lo + hi) / 2


def _div_hist_axis_center(div_hist_max):
    return (div_hist_max * 1.08) / 2


def _add_subplot_traces(fig, source_fig, row, col, *, text_x=None, text_y=None):
    """Copy traces into a subplot, recentering text-mode empty-state labels."""
    for trace in source_fig.data:
        if (
            text_x is not None
            and text_y is not None
            and getattr(trace, 'mode', None) == 'text'
        ):
            textfont = trace.textfont
            if not textfont:
                textfont = {'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'}
            fig.add_trace(
                go.Scatter(
                    x=[text_x],
                    y=[text_y],
                    text=list(trace.text) if trace.text is not None else [''],
                    mode='text',
                    textfont=textfont,
                    textposition='middle center',
                    showlegend=False,
                    hoverinfo='none',
                ),
                row=row,
                col=col,
            )
        else:
            fig.add_trace(trace, row=row, col=col)


def _assemble_teaid_2x2_figure(
    plot_main_title,
    cons_len,
    panels,
    fig_dot,
    *,
    x_axis_label=DEFAULT_TEAID_X_AXIS_LABEL,
):
    divergence_max = panels['divergence_max']
    coverage_max = panels['coverage_max']
    div_hist_max = panels['div_hist_max']
    blast_plot_title = panels['blast_plot_title']

    panel_px = 380
    gap_between_row = 90
    top_margin = 100
    bottom_margin = 70
    side_margin = 70
    plot_block = 2 * panel_px + gap_between_row
    total_fig_height = plot_block + top_margin + bottom_margin
    total_fig_width = total_fig_height
    vertical_spacing = gap_between_row / plot_block
    consensus_ticks = _consensus_axis_tick_kwargs(cons_len)
    consensus_center = _consensus_axis_center(cons_len)

    fig = make_subplots(
        rows=2,
        cols=2,
        subplot_titles=(
            blast_plot_title,
            'BLAST Coverage Plot',
            'Self Dot Plot',
            'Divergence Plot',
        ),
        vertical_spacing=vertical_spacing,
        horizontal_spacing=0.10,
        row_heights=[1, 1],
        column_widths=[1, 1],
    )

    _add_subplot_traces(
        fig,
        panels['fig_blast'],
        1,
        1,
        text_x=consensus_center,
        text_y=_divergence_axis_center(divergence_max),
    )
    _add_subplot_traces(
        fig,
        panels['fig_coverage'],
        1,
        2,
        text_x=consensus_center,
        text_y=_coverage_axis_center(coverage_max),
    )
    _add_subplot_traces(
        fig,
        fig_dot,
        2,
        1,
        text_x=consensus_center,
        text_y=consensus_center,
    )
    _add_subplot_traces(
        fig,
        panels['fig_div_hist'],
        2,
        2,
        text_x=TEAID_DIVERGENCE_HIST_MAX / 2,
        text_y=_div_hist_axis_center(div_hist_max),
    )

    fig.update_annotations(font=_SUBPLOT_TITLE_FONT)
    _configure_blast_subplot(
        fig, 1, 1, cons_len, divergence_max, consensus_ticks, x_axis_label=x_axis_label
    )
    _configure_coverage_subplot(
        fig, 1, 2, cons_len, coverage_max, consensus_ticks, x_axis_label=x_axis_label
    )
    _configure_self_dot_subplot(fig, 2, 1, cons_len, x_axis_label=x_axis_label)
    _configure_divergence_subplot(fig, 2, 2, div_hist_max)

    fig.update_layout(
        height=total_fig_height,
        width=total_fig_width,
        title_text=plot_main_title,
        title_font={'size': 26, 'family': 'Arial, sans-serif', 'color': 'black'},
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin={'t': top_margin, 'b': bottom_margin, 'l': side_margin, 'r': side_margin},
    )
    return fig


def _assemble_teaid_pan_genome_figure(
    plot_main_title,
    cons_len,
    genome_panels,
    fig_dot,
    *,
    x_axis_label=DEFAULT_TEAID_X_AXIS_LABEL,
):
    """
    Pan-genome layout: four rows (BLAST hits, coverage, divergence, self dot),
    one column per genome. Self dot sits in the bottom-left cell only.
    """
    n = len(genome_panels)
    if n < 1:
        raise ValueError('Pan-genome TE-Aid requires at least one genome.')

    specs = [[{} for _ in range(n)] for _ in range(4)]
    specs[3] = [{}] + [None] * (n - 1)

    subplot_titles: list[str] = []
    for name, panels in genome_panels:
        subplot_titles.append(
            _pan_blast_subplot_title(name, panels['blast_plot_title'])
        )
    for name, _panels in genome_panels:
        subplot_titles.append(f'{name}:<br>BLAST Coverage Plot')
    for name, _panels in genome_panels:
        subplot_titles.append(f'{name}:<br>Divergence Plot')
    subplot_titles.append('Self Dot Plot')

    panel_px = _pan_genome_panel_px(n)
    row_gap = 24
    top_margin = 100
    bottom_margin = 90
    side_margin = 70
    plot_block = 4 * panel_px + 3 * row_gap
    vertical_spacing = row_gap / (panel_px + row_gap)
    total_fig_width = n * panel_px + side_margin * 2
    total_fig_height = plot_block + top_margin + bottom_margin
    consensus_ticks = _consensus_axis_tick_kwargs(cons_len)
    consensus_center = _consensus_axis_center(cons_len)

    divergence_max = max(panels['divergence_max'] for _name, panels in genome_panels)
    coverage_max = max(panels['coverage_max'] for _name, panels in genome_panels)
    div_hist_max = max(panels['div_hist_max'] for _name, panels in genome_panels)

    blast_text_y = _divergence_axis_center(divergence_max)
    coverage_text_y = _coverage_axis_center(coverage_max)
    div_text_x = TEAID_DIVERGENCE_HIST_MAX / 2
    div_text_y = _div_hist_axis_center(div_hist_max)

    fig = make_subplots(
        rows=4,
        cols=n,
        specs=specs,
        subplot_titles=tuple(subplot_titles),
        vertical_spacing=vertical_spacing,
        horizontal_spacing=0.04,
        row_heights=[1, 1, 1, 1],
        column_widths=[1] * n,
    )

    for col_idx, (_name, panels) in enumerate(genome_panels, start=1):
        left_col = col_idx == 1
        _add_subplot_traces(
            fig,
            panels['fig_blast'],
            1,
            col_idx,
            text_x=consensus_center,
            text_y=blast_text_y,
        )
        _add_subplot_traces(
            fig,
            panels['fig_coverage'],
            2,
            col_idx,
            text_x=consensus_center,
            text_y=coverage_text_y,
        )
        _add_subplot_traces(
            fig,
            panels['fig_div_hist'],
            3,
            col_idx,
            text_x=div_text_x,
            text_y=div_text_y,
        )
        _configure_blast_subplot(
            fig,
            1,
            col_idx,
            cons_len,
            divergence_max,
            consensus_ticks,
            show_x_title=False,
            show_y_title=left_col,
            x_tickangle=TEAID_PAN_X_TICKANGLE,
        )
        _configure_coverage_subplot(
            fig,
            2,
            col_idx,
            cons_len,
            coverage_max,
            consensus_ticks,
            show_x_title=False,
            show_y_title=left_col,
            x_tickangle=TEAID_PAN_X_TICKANGLE,
        )
        _configure_divergence_subplot(
            fig,
            3,
            col_idx,
            div_hist_max,
            show_x_title=left_col,
            show_y_title=left_col,
            x_tickangle=TEAID_PAN_X_TICKANGLE,
        )

    _add_subplot_traces(
        fig,
        fig_dot,
        4,
        1,
        text_x=consensus_center,
        text_y=consensus_center,
    )
    _configure_self_dot_subplot(
        fig, 4, 1, cons_len, x_tickangle=TEAID_PAN_X_TICKANGLE, x_axis_label=x_axis_label
    )

    fig.update_annotations(font={'size': 12, 'family': 'Arial, sans-serif', 'color': 'black'})
    genome_names = ', '.join(name for name, _panels in genome_panels)
    fig.update_layout(
        height=total_fig_height,
        width=total_fig_width,
        title_text=f'{plot_main_title} | Pan-genome ({genome_names})',
        title_font={'size': 24, 'family': 'Arial, sans-serif', 'color': 'black'},
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin={'t': top_margin, 'b': bottom_margin, 'l': side_margin, 'r': side_margin},
        autosize=False,
    )
    return fig


def _write_teaid_html(fig, output_dir, input_file, *, output_name_tag=None):
    from repeatviewer.pan_genome import genome_output_suffix

    tag = genome_output_suffix(output_name_tag) if output_name_tag else ''
    output_html_path = os.path.join(
        output_dir, f'{os.path.basename(input_file)}{tag}_TEAid.html'
    )
    pio.write_html(
        fig,
        file=output_html_path,
        auto_open=False,
        config={'responsive': False},
    )
    return output_html_path


def teaid_pan_genome_plotter(
    input_file,
    output_dir,
    genome_targets,
    e_value=1e-40,
    num_threads=5,
    use_system_blast=False,
    self_blast_word_size=DEFAULT_TEAID_SELF_BLAST_WORD_SIZE,
    x_axis_label=DEFAULT_TEAID_X_AXIS_LABEL,
    hit_len_high_pct=80.0,
    hit_len_mid_pct=50.0,
    hit_color_high='red',
    hit_color_mid='#1f77b4',
    hit_line_transparency=0.0,
    progress=None,
):
    """
    Combined TE-Aid report for pan-genome mode.

    ``genome_targets`` is a list of ``(short_name, genome_fasta_path)`` pairs.
    """
    from repeatviewer.backends.job_phases import TEAID_PHASES, emit_phase

    base = os.path.basename(input_file)
    logging.info(
        'TE-Aid pan-genome is running for:%s (%d genomes)',
        base,
        len(genome_targets),
    )

    try:
        emit_phase(progress, 0, TEAID_PHASES, f"{base}: preparing consensus…")
        con_seq, cons_len = con_generater(input_file, output_dir, threshold=0.3)
    except TeaidInputError:
        raise
    except Exception:
        logging.error(
            'An error for consensus sequence generation: \n %s',
            traceback.format_exc(),
        )
        return False

    plot_main_title = f'TE: {os.path.basename(input_file)} | {cons_len} bp'

    emit_phase(progress, 1, TEAID_PHASES, f"{base}: genome BLAST panels…")
    genome_panels = []
    for name, genome_file in genome_targets:
        panels = _build_genome_teaid_panels(
            con_seq,
            cons_len,
            genome_file,
            output_dir,
            e_value,
            num_threads=num_threads,
            use_system_blast=use_system_blast,
            output_label=name,
            hit_len_high_pct=hit_len_high_pct,
            hit_len_mid_pct=hit_len_mid_pct,
            hit_color_high=hit_color_high,
            hit_color_mid=hit_color_mid,
            hit_line_transparency=hit_line_transparency,
        )
        genome_panels.append((name, panels))

    emit_phase(progress, 2, TEAID_PHASES, f"{base}: self-dotplot…")
    fig_dot = _build_self_dot_figure(
        con_seq,
        output_dir,
        cons_len,
        use_system_blast=use_system_blast,
        self_blast_word_size=self_blast_word_size,
    )
    emit_phase(progress, 3, TEAID_PHASES, f"{base}: assembling TE-Aid figure…")
    fig = _assemble_teaid_pan_genome_figure(
        plot_main_title,
        cons_len,
        genome_panels,
        fig_dot,
        x_axis_label=x_axis_label,
    )
    emit_phase(progress, 4, TEAID_PHASES, f"{base}: writing HTML…")
    _write_teaid_html(fig, output_dir, input_file)
    emit_phase(progress, TEAID_PHASES, TEAID_PHASES, f"{base}: finished")
    logging.info('TE-Aid pan-genome finished: %s.', os.path.basename(input_file))
    return True


def teaid_plotter(
    input_file,
    output_dir,
    genome_file,
    current_win,
    e_value=1e-40,
    num_threads=5,
    use_system_blast=False,
    self_blast_word_size=DEFAULT_TEAID_SELF_BLAST_WORD_SIZE,
    conduct_cdd=False,
    prepared_cdd=None,
    output_name_tag=None,
    x_axis_label=DEFAULT_TEAID_X_AXIS_LABEL,
    hit_len_high_pct=80.0,
    hit_len_mid_pct=50.0,
    hit_color_high='red',
    hit_color_mid='#1f77b4',
    hit_line_transparency=0.0,
    cancel=None,
    progress=None,
):
    # Used global variable: os_type

    # output_dir refers to the temp folders

    from repeatviewer.backends.job_phases import TEAID_PHASES, emit_phase

    run_succeed = True
    base = os.path.basename(input_file)

    logging.info(f'TE-Aid is running for:{base}')

    try:
        emit_phase(progress, 0, TEAID_PHASES, f"{base}: preparing consensus…")
        # Generate consensus sequence
        con_seq, cons_len = con_generater(input_file, output_dir, threshold=0.3)

    except TeaidInputError as e:
        logging.error("TE-Aid input: %s", e)
        raise
    except Exception as e:
        logging.error(
            f'An error for consensus sequence generation: \n {traceback.format_exc()}'
        )
        logging.error(
            'Consensus sequence generation failed: %s',
            e,
        )
        return False

    # Define main plot title
    plot_main_title = f'TE: {os.path.basename(input_file)} | {cons_len} bp'

    emit_phase(progress, 1, TEAID_PHASES, f"{base}: genome BLAST panels…")
    panels = _build_genome_teaid_panels(
        con_seq,
        cons_len,
        genome_file,
        output_dir,
        e_value,
        num_threads=num_threads,
        use_system_blast=use_system_blast,
        output_label=output_name_tag,
        hit_len_high_pct=hit_len_high_pct,
        hit_len_mid_pct=hit_len_mid_pct,
        hit_color_high=hit_color_high,
        hit_color_mid=hit_color_mid,
        hit_line_transparency=hit_line_transparency,
        cancel=cancel,
    )
    emit_phase(progress, 2, TEAID_PHASES, f"{base}: self-dotplot…")
    fig_dot = _build_self_dot_figure(
        con_seq,
        output_dir,
        cons_len,
        use_system_blast=use_system_blast,
        self_blast_word_size=self_blast_word_size,
        cancel=cancel,
    )
    emit_phase(progress, 3, TEAID_PHASES, f"{base}: assembling TE-Aid figure…")
    fig = _assemble_teaid_2x2_figure(
        plot_main_title, cons_len, panels, fig_dot, x_axis_label=x_axis_label
    )
    emit_phase(progress, 4, TEAID_PHASES, f"{base}: writing HTML…")
    _write_teaid_html(fig, output_dir, input_file, output_name_tag=output_name_tag)
    emit_phase(progress, TEAID_PHASES, TEAID_PHASES, f"{base}: finished")

    logging.info(f'TE-Aid finished: {os.path.basename(input_file)}.')
    # fig.show()
    return run_succeed
