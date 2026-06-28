import logging
import os
import platform
import traceback
import warnings
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import plotly.io as pio
from Bio import AlignIO, BiopythonDeprecationWarning, SeqIO
from Bio.Seq import Seq
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord
from plotly.subplots import make_subplots

from teviewer.backends.gui_functions import blast, check_database
from teviewer.job_options import (
    DEFAULT_TE_DIVERGENCE_BIN_WIDTH,
    clamp_teaid_self_blast_word_size,
)

# Suppress all deprecation warnings
warnings.filterwarnings('ignore', category=BiopythonDeprecationWarning)

"""
The BLAST, coverate, and self-dot plot codes are derived from TE-Aid. https://github.com/clemgoub/TE-Aid.git

Developer: Jiangzhao Qian
"""

os_type = platform.system()

TEAID_DIVERGENCE_HIST_MAX = 0.6


def _consensus_position_ticks(max_len: int) -> list[int]:
    """1-based TE consensus tick values; always includes 1, never the exact sequence end."""
    max_len = max(1, int(max_len))
    if max_len <= 100:
        step = 10
    elif max_len <= 250:
        step = 25
    elif max_len <= 500:
        step = 50
    elif max_len <= 1000:
        step = 100
    elif max_len <= 2500:
        step = 500
    else:
        step = 1000
    ticks = [1]
    pos = step
    while pos < max_len:
        ticks.append(pos)
        pos += step
    return ticks


def _consensus_axis_tick_kwargs(max_len: int) -> dict:
    return {
        'tickmode': 'array',
        'tickvals': _consensus_position_ticks(max_len),
        'tickangle': 0,
    }


def con_generater(input_file, output_dir, threshold=0.8, ambiguous='N'):
    """
    Mimic a 'dumb consensus' without using AlignInfo.dumb_consensus.
    Counts only A/C/G/T, ignores gaps and other symbols.
    Returns a Bio.Seq.Seq (uppercased), same as the upper function.
    """
    alignment = AlignIO.read(input_file, "fasta")
    n_cols = alignment.get_alignment_length()

    consensus_chars = []
    for i in range(n_cols):
        col = alignment[:, i].upper()
        bases = [b for b in col if b in {"A", "C", "G", "T"}]

        if not bases:  # no A/C/G/T at this column
            consensus_chars.append(ambiguous)
            continue

        a = bases.count("A")
        c = bases.count("C")
        g = bases.count("G")
        t = bases.count("T")
        counts = {"A": a, "C": c, "G": g, "T": t}

        winner = max(counts, key=counts.get)
        freq = counts[winner] / len(bases)

        consensus_chars.append(winner if freq >= threshold else ambiguous)

    consensus_seq = Seq("".join(consensus_chars)).upper()

    consensus_seq_len = len(consensus_seq)

    consensus_seq_obj = SeqRecord(consensus_seq, id="consensus", description=f"threshold={threshold};amb={ambiguous}")

    # Write consensus sequence to a FASTA file
    output_file = os.path.join(output_dir, f'{os.path.basename(input_file)}_{threshold}_co.fa')
    with open(output_file, 'w') as file:
        SeqIO.write(consensus_seq_obj, file, 'fasta')

    # Return a Bio.Seq.Seq in uppercase (to match the upper function’s output type)
    return output_file, consensus_seq_len


#####################################################################################################
# Code block: genome blast
#####################################################################################################
def teaid_genome_blast(input_file, genome_file, output_dir, e_value, num_threads=1,
                       use_system_blast=False):
    # input_file: fasta file contain only one sequence
    # Used global variable: os_type

    if genome_file is None:
        logging.error("WARNING: Genome file not found. TEAid can't perform genome blast.")
        return False

    # Check if provided genome file exist
    elif not os.path.isfile(genome_file):
        logging.error(
            "WARNING: You provided genome file not found. TEAid can't perform genome blast."
        )
        return False

    # Check genome blast database availability
    blast_database = check_database(genome_file, os_type=os_type, use_system_blast = use_system_blast)

    # Check if makeblastdb is correctly installed
    if blast_database == 'makeblastdb_not_found':
        logging.error(
            "Error: makeblastdb command not found. TEAid can't perform genome blast."
        )
        return False

    # Check if error happened
    elif blast_database == 'makeblastdb_got_error':
        logging.error(
            "Error: BLAST database can't be established for the genome. TEAid can't perform genome blast."
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
        use_system_blast=use_system_blast
    )

    if genome_blast_out == 'blastn_not_found':
        logging.error("Error: BLAST command not found. TEAid can't perform genome blast.")
        return False

    elif genome_blast_out == 'blastn_got_error':
        logging.error(
            "TEAid can't perform genome blast. Refer to terminal for more information."
        )
        return False

    elif genome_blast_out == 'blast_n_zero':
        logging.warning('Genome blast hit number is 0 for this sequence.')
        return 'blast_n_zero'

    elif not genome_blast_out:
        logging.error(
            "TEAid can't perform genome blast. Refer to terminal for more information."
        )
        return False

    return genome_blast_out


#####################################################################################################
# Code block: self blast
#####################################################################################################
def self_blast(input_file, output_dir, use_system_blast=False, word_size=11):
    # input_file: fasta file contain only one sequence
    # Used global variable: os_type

    input_n = os.path.basename(input_file)

    output_dir_self_blast = os.path.join(output_dir, f'{input_n}_self_blast')

    os.makedirs(output_dir_self_blast, exist_ok=True)

    # Build blast database
    blast_database = check_database(input_file, output_dir_self_blast, os_type=os_type,
                                    use_system_blast=use_system_blast)

    # Check if makeblastdb is correctly installed
    if blast_database == 'makeblastdb_not_found':
        logging.error(
            "Error: makeblastdb command not found. TEAid can't perform self blast."
        )
        return False

    # Check if error happened
    elif blast_database == 'makeblastdb_got_error':
        logging.error(
            "Error: BLAST database can't be established for the genome. TEAid can't perform self blast."
        )
        return False

    self_blast_out, _ = blast(
        input_file,
        blast_database,
        output_dir_self_blast,
        e_value=0.05,
        self_blast=True,
        word_size=word_size,
        os_type=os_type,
        use_system_blast=use_system_blast
    )

    if self_blast_out == 'blastn_not_found':
        logging.error("Error: BLAST command not found. TEAid can't perform self blast.")
        return False

    elif self_blast_out == 'blastn_got_error':
        logging.error(
            "TEAid can't perform self blast. Refer to terminal for more information."
        )
        return False

    elif self_blast_out == 'blast_n_zero':
        logging.warning('Self blast hit number is 0 for this sequence.')
        return 'blast_n_zero'

    elif not self_blast_out:
        logging.error(
            "TEAid can't perform self blast. Refer to terminal for more information."
        )
        return False

    return self_blast_out


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


def blast_plot(df, cons_len, full_len_thr=0.8):
    # Calculate full matches based on the threshold and consensus length
    full = df[(abs(df['V7'] - df['V8']) >= full_len_thr * cons_len)]

    # Calculate the max number of divergence for the y-axis range
    y_max = max(100 - df['V3'])

    # Create a figure
    fig = go.Figure()

    # Add segments for all blast hits
    for i in df.index:
        fig.add_trace(
            go.Scatter(
                x=[df.loc[i, 'V7'], df.loc[i, 'V8']],
                y=[100 - df.loc[i, 'V3'], 100 - df.loc[i, 'V3']],
                mode='lines',
                line={'color': 'black'},
                hoverinfo='text',
                text=create_hovertext(df.loc[i]),
                showlegend=False,
            )
        )

    # Add segments for full matches
    for i in full.index:
        fig.add_trace(
            go.Scatter(
                x=[full.loc[i, 'V7'], full.loc[i, 'V8']],
                y=[100 - full.loc[i, 'V3'], 100 - full.loc[i, 'V3']],
                mode='lines',
                line={'color': 'red', 'width': 3},
                hoverinfo='text',
                text=create_hovertext(full.loc[i]),
                showlegend=False,
            )
        )

    # Define plot title
    plot_title = f'Fragments number: {len(df)}; full length: {len(full)} (>={int(full_len_thr * cons_len)}bp)'

    # Update the layout to add titles and labels with adjusted axis ranges
    fig.update_layout(
        title=plot_title,
        xaxis_title='TE consensus (bp)',
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
        xaxis_title='Position along TE consensus (bp)',
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
    if custom_text:
        fig.add_trace(
            go.Scatter(
                x=[0.5],
                y=[5],
                text=[custom_text],
                mode='text',
                textfont={'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
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
            'range': [0, 10],
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
    return fig, 10


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
        xaxis_title='TE consensus self dotplot (bp)',
        yaxis_title='TE consensus self dotplot (bp)',
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


#####################################################################################################
# Code block: rpstblasn plot
#####################################################################################################
def plot_rpsblast_hits(input_file, cons_len):
    """
    Create a TE consensus hit map plot from RPS-BLAST results with score-based coloring and hit directions.

    Parameters:
    - input_file: rpstblastn result file
    - cons_len: Integer representing the length of the TE consensus sequence.
    """

    #########################################################
    # Code block: Format dataframe
    #########################################################

    df_rps = pd.read_csv(input_file, sep='\t', header=None)
    # Assign column names
    df_rps.columns = [
        'qseqid',
        'sseqid',
        'pident',
        'length',
        'mismatch',
        'gapopen',
        'qstart',
        'qend',
        'sstart',
        'send',
        'evalue',
        'score',
        'description',
    ]

    # Convert positions to integers if they're not already
    df_rps['qstart'] = df_rps['qstart'].astype(int)
    df_rps['qend'] = df_rps['qend'].astype(int)

    # Convert Score to integers
    df_rps['score'] = df_rps['score'].astype(int)

    # Add 'start' and 'end' columns to handle reverse hits correctly
    df_rps['start'] = df_rps[['qstart', 'qend']].min(axis=1)
    df_rps['end'] = df_rps[['qstart', 'qend']].max(axis=1)

    # Function to assign tracks to hits to prevent overlaps
    def assign_tracks(df_rps_assign):
        tracks = []  # List to keep track of end positions in each track
        df_rps_assign_sorted = df_rps_assign.sort_values(
            by=['start', 'end']
        ).reset_index(drop=True)
        track_assignment = []  # list records which track each hit is assigned to

        for _idx, hit in df_rps_assign_sorted.iterrows():
            assigned = False
            for track_num, track_end in enumerate(tracks):
                if hit['start'] > track_end + (
                    cons_len / 20
                ):  # Adjusted gap to prevent overlaps
                    # Assign hit to this track
                    tracks[track_num] = hit['end']
                    track_assignment.append(track_num)
                    assigned = True
                    break
            if not assigned:
                # Create new track
                tracks.append(hit['end'])
                track_assignment.append(len(tracks) - 1)

        df_rps_assign_sorted['track'] = track_assignment
        return df_rps_assign_sorted

    # Assign tracks to hits
    df_rps_with_tracks = assign_tracks(df_rps)
    num_tracks = df_rps_with_tracks['track'].max() + 1

    # Adjust 'arrow_width' accordingly
    arrow_width = 0.12  # Keep arrow width proportional to track spacing

    #########################################################
    # Code block: Define arrow color and figure legend
    #########################################################
    score_ranges = [
        {
            'min_score': float('-inf'),
            'max_score': -1000,
            'color': 'white',
            'label': 'Alignment Score',
        },  # Show text only
        {'min_score': -1000, 'max_score': 40, 'color': 'black', 'label': '< 40'},
        {'min_score': 40, 'max_score': 50, 'color': 'blue', 'label': '40 - 50'},
        {'min_score': 50, 'max_score': 80, 'color': 'green', 'label': '50 - 80'},
        {'min_score': 80, 'max_score': 200, 'color': 'pink', 'label': '80 - 200'},
        {'min_score': 200, 'max_score': float('inf'), 'color': 'red', 'label': '≥ 200'},
    ]

    # Function to get score range and the color
    def get_score_range(score):
        for sr in score_ranges:
            if sr['min_score'] <= score < sr['max_score']:
                return sr
        return score_ranges[-1]  # Default to last score range

    #########################################################
    # Code block: Plotting
    #########################################################

    def wrap_text(text, max_length=70):
        """
        Efficiently wrap a long string into multiple lines with a specified maximum length.
        """
        words = text.split()
        wrapped_lines = []
        current_line = []
        current_length = 0  # Running total of the current line's length

        for word in words:
            word_length = len(word)

            # Check if adding the next word exceeds the max length
            if current_length + word_length + len(current_line) > max_length:
                wrapped_lines.append(' '.join(current_line))
                current_line = []
                current_length = 0
            current_line.append(word)
            current_length += word_length

        # Append any remaining words
        if current_line:
            wrapped_lines.append(' '.join(current_line))

        return '<br>'.join(wrapped_lines)

    # Create figure
    fig = go.Figure()

    # Add figure legend
    """
    for sr in score_ranges:
        fig.add_trace(go.Scatter(
            x=[None],
            y=[None],
            mode='markers',
            marker=dict(
                size=15,
                color=sr['color'],
                symbol='square'
            ),
            name=sr['label'],
            showlegend=True
        ))
    """
    #########################################################
    # Code block: Start to loop each line of dataframe
    #########################################################
    # Initialize list to hold annotations
    annotations = []

    for _idx, hit in df_rps_with_tracks.iterrows():
        x_start = hit['qstart']
        x_end = hit['qend']
        y = hit['track']

        # Determine direction
        if hit['qstart'] <= hit['qend']:
            direction = 'forward'
            plot_x_start = x_start
            plot_x_end = x_end
        else:
            direction = 'reverse'
            plot_x_start = x_end
            plot_x_end = x_start

        # Get score range and color
        score_range = get_score_range(hit['score'])
        color = score_range['color']

        # Hover text with detailed information
        wrapped_description = wrap_text(hit['description'], max_length=70)

        trace_hovertext = (
            f'Query Start: {hit["qstart"]}<br>'
            f'Query End: {hit["qend"]}<br>'
            f'Percent Identity: {hit["pident"]:.2f}<br>'
            f'E-value: {hit["evalue"]:.2e}<br>'
            f'Score: {hit["score"]}<br>'
            f'Description:<br>{wrapped_description}'
        )

        # Brief info to display under the bar (using the second element)
        description_parts = hit['description'].split(',')
        if len(description_parts) > 1:
            brief_info = description_parts[
                1
            ].strip()  # Take the second part of the description
        else:
            brief_info = description_parts[0].strip()

        # Adjust arrowhead length
        arrow_head_length = (plot_x_end - plot_x_start) * 0.25  # 25% of the bar length

        if direction == 'forward':
            # Coordinates for forward arrow
            x_coords = [
                plot_x_start,
                plot_x_end - arrow_head_length,
                plot_x_end - arrow_head_length,
                plot_x_end,
                plot_x_end - arrow_head_length,
                plot_x_end - arrow_head_length,
                plot_x_start,
                plot_x_start,
            ]
            y_coords = [
                y + arrow_width,
                y + arrow_width,
                y + arrow_width * 1.5,
                y,
                y - arrow_width * 1.5,
                y - arrow_width,
                y - arrow_width,
                y + arrow_width,
            ]
        else:
            # Coordinates for reverse arrow
            x_coords = [
                plot_x_end,
                plot_x_start + arrow_head_length,
                plot_x_start + arrow_head_length,
                plot_x_start,
                plot_x_start + arrow_head_length,
                plot_x_start + arrow_head_length,
                plot_x_end,
                plot_x_end,
            ]
            y_coords = [
                y + arrow_width,
                y + arrow_width,
                y + arrow_width * 1.5,
                y,
                y - arrow_width * 1.5,
                y - arrow_width,
                y - arrow_width,
                y + arrow_width,
            ]

        # Add the arrow as a filled scatter plot
        fig.add_trace(
            go.Scatter(
                x=x_coords,
                y=y_coords,
                mode='lines+markers',
                marker={'size': 0.1, 'opacity': 0},  # Invisible markers
                fill='toself',
                fillcolor=color,
                line={'color': 'black'},
                hoverinfo='text',
                hovertext=trace_hovertext,
                showlegend=False,  # Keep legend entries from arrows
            )
        )

        # Adjust short description position
        x_mid = (plot_x_start + plot_x_end) / 2
        annotation_y = y - 0.52  # Position annotation above the arrow

        # Add annotation for brief info under the bar
        annotations.append(
            {
                'x': x_mid,
                'y': annotation_y,
                'text': brief_info,
                'showarrow': False,
                'font': {'size': 10},
                'xanchor': 'center',
                'yanchor': 'bottom',  # Adjusted for placement
            }
        )

    # Add annotations to the figure
    fig.update_layout(annotations=annotations)

    # Update layout with adjusted y-axis range
    fig.update_layout(
        title=f'{os.path.basename(input_file)}',
        # title_x=0.5,  # Center the title
        xaxis_title='TE consensus (bp)',
        yaxis_title='TE Protein Domains',
        font={'family': 'Arial, sans-serif', 'size': 14, 'color': 'black'},
        xaxis={
            'range': [1, cons_len],
            'side': 'bottom',
            'showgrid': False,
            'zeroline': False,
            'showticklabels': True,
            'linecolor': 'black',
            'mirror': True,
            'linewidth': 2,
        },
        yaxis={
            # visible=False,
            # range=[max_y, min_y],
            'autorange': True,
            'fixedrange': False,
            'showticklabels': False,
            'linecolor': 'black',
            'mirror': True,
            'linewidth': 2,
        },
        plot_bgcolor='white',
        width=1600,
        # height=600,
        height=(165 + (num_tracks - 1) * 40),  # Adjusted figure height
        margin={'l': 60, 'r': 60, 't': 60, 'b': 60},  # Increased top margin for legend
        # Customize hover label appearance at the layout level
        hoverlabel={
            'bgcolor': 'lightgray',  # White background
            'font': {
                'size': 12,  # Font size
                'color': 'black',  # Black text
            },
            'bordercolor': 'black',  # black border
        },
    )

    return fig, num_tracks


#####################################################################################################
# Code block: TEtrimmer GUI plotter, integrate all plots
#####################################################################################################
def teaid_plotter(
    input_file,
    output_dir,
    genome_file,
    current_win,
    e_value=1e-40,
    num_threads=5,
    use_system_blast=False,
    self_blast_word_size=11,
    conduct_cdd=False,
    prepared_cdd=None,
):
    # Used global variable: os_type

    # output_dir refers to the temp folders

    run_succeed = True

    logging.info(f'TEAid is running for:{os.path.basename(input_file)}')

    try:
        # Generate consensus sequence
        con_seq, cons_len = con_generater(input_file, output_dir, threshold=0.3)

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

    # Perform genome blast
    genome_blast_out = teaid_genome_blast(
        con_seq, genome_file, output_dir, e_value, num_threads=num_threads, use_system_blast=use_system_blast
    )

    if genome_blast_out:
        # Check if blast hit number is zero
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

            # Genome blast plot, coverage, and divergence histogram
            fig_blast, divergence_max, blast_plot_title = blast_plot(
                df_genome_blast_out, cons_len
            )

            fig_coverage, coverage_max = coverage_plot(df_genome_blast_out, cons_len)
            fig_div_hist, div_hist_max = divergence_histogram_plot(df_genome_blast_out)
    else:
        fig_blast = empty_plot(
            cons_len,
            custom_text='An error occurred. Please check the terminal for details.<br> '
            'Please provide a genome file if it has not been specified.',
        )
        divergence_max = 10
        blast_plot_title = 'BLAST hit Plot'

        fig_coverage = empty_plot(
            cons_len,
            custom_text='An error occurred. Please check the terminal for details.<br> '
            'Please provide a genome file if it has not been specified.',
        )
        coverage_max = 10
        fig_div_hist, div_hist_max = empty_divergence_histogram(
            'An error occurred. Please check the terminal for details.'
        )

    # Perform self blast (word size clamped 4–20; lowered on short consensus)
    self_word_size = clamp_teaid_self_blast_word_size(self_blast_word_size)
    dotplot_blast = self_blast(
        con_seq,
        output_dir,
        use_system_blast=use_system_blast,
        word_size=self_word_size,
    )

    if dotplot_blast == 'blast_n_zero':
        fig_dot = empty_plot(
            cons_len,
            custom_text='No self BLAST hits found for this sequence.',
        )
    elif dotplot_blast:
        # self blast dot plot
        df_dotplot_blast = pd.read_csv(dotplot_blast)
        fig_dot = dot_plot(df_dotplot_blast)
    else:
        fig_dot = empty_plot(
            cons_len,
            custom_text='An error occurred. Please check the terminal for details.',
        )

    # 2×2 layout: equal square panels with extra space between rows.
    panel_px = 380
    gap_between_row = 90
    top_margin = 100
    bottom_margin = 70
    side_margin = 70
    plot_block = 2 * panel_px + gap_between_row
    total_fig_height = plot_block + top_margin + bottom_margin
    total_fig_width = total_fig_height
    vertical_spacing = gap_between_row / plot_block

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

    for trace in fig_blast.data:
        fig.add_trace(trace, row=1, col=1)
    for trace in fig_coverage.data:
        fig.add_trace(trace, row=1, col=2)
    for trace in fig_dot.data:
        fig.add_trace(trace, row=2, col=1)
    for trace in fig_div_hist.data:
        fig.add_trace(trace, row=2, col=2)

    fig.update_annotations(
        font={'size': 13, 'family': 'Arial, sans-serif', 'color': 'black'},
    )

    consensus_ticks = _consensus_axis_tick_kwargs(cons_len)
    # Configuring subplot 1
    fig.update_xaxes(
        title_text='TE consensus (bp)',
        showgrid=False,
        showticklabels=True,
        title_font={'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
        range=[1, cons_len],
        **consensus_ticks,
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=1,
        col=1,
    )
    fig.update_yaxes(
        title_text='Divergence to consensus (%)',
        showgrid=False,
        showticklabels=True,
        title_font={'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
        range=[-(divergence_max / 30), divergence_max + (divergence_max / 30)],
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=1,
        col=1,
    )

    # Configuring subplot 2 (coverage)
    fig.update_xaxes(
        title_text='TE consensus (bp)',
        showgrid=False,
        showticklabels=True,
        title_font={'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
        range=[1, cons_len],
        **consensus_ticks,
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=1,
        col=2,
    )
    fig.update_yaxes(
        title_text='Coverage (bp)',
        showgrid=False,
        showticklabels=True,
        title_font={'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
        range=[-(coverage_max / 30), coverage_max + (coverage_max / 30)],
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=1,
        col=2,
    )

    # Configuring subplot 3 (self dot plot: 1-based coordinates)
    fig.update_xaxes(
        title_text='TE consensus (bp)',
        showgrid=False,
        showticklabels=True,
        title_font={'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
        range=[1, cons_len],
        **consensus_ticks,
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=2,
        col=1,
    )
    fig.update_yaxes(
        title_text='TE consensus (bp)',
        showgrid=False,
        showticklabels=True,
        title_font={'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
        range=[1, cons_len],
        tickmode='array',
        tickvals=_consensus_position_ticks(cons_len),
        tickangle=0,
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=2,
        col=1,
    )

    # Configuring subplot 4 (divergence plot)
    fig.update_xaxes(
        title_text='Divergence',
        showgrid=True,
        gridcolor='#f0f0f0',
        showticklabels=True,
        title_font={'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
        range=[0, TEAID_DIVERGENCE_HIST_MAX],
        tickangle=0,
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=2,
        col=2,
    )
    fig.update_yaxes(
        title_text='Total length (bp)',
        showgrid=True,
        gridcolor='#f4f4f4',
        showticklabels=True,
        title_font={'size': 16, 'family': 'Arial, sans-serif', 'color': 'black'},
        range=[0, div_hist_max * 1.08],
        linewidth=1,
        linecolor='black',
        mirror=True,
        row=2,
        col=2,
    )

    fig.update_layout(
        height=total_fig_height,
        width=total_fig_width,
        title_text=plot_main_title,
        title_font={'size': 26, 'family': 'Arial, sans-serif', 'color': 'black'},
        plot_bgcolor='white',
        paper_bgcolor='white',
        margin={'t': top_margin, 'b': bottom_margin, 'l': side_margin, 'r': side_margin},
    )
    # Save the figure as an HTML file
    output_html_path = os.path.join(
        output_dir, f'{os.path.basename(input_file)}_TEAid.html'
    )
    pio.write_html(fig, file=output_html_path, auto_open=False)

    logging.info(f'TEAid finished: {os.path.basename(input_file)}.')
    # fig.show()
    return run_succeed
