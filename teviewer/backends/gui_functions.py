import logging
import os
import random
import shutil
import subprocess
import sys
import threading
import traceback

from teviewer.backends.job_cancel import JobCancelledError, run_subprocess_with_cancel
from teviewer.external_process import run_logged_process
import pandas as pd
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

special_character = {
    '/': '__',
    '\\': '__',
    ' ': '_',
    '|': '_',
    '<': '_',
    ':': '_',
    '"': '_',
    '*': '_',
    '?': '_',
    '#': '_',
}

_blast_db_locks: dict[str, threading.Lock] = {}
_blast_db_locks_guard = threading.Lock()

# Cache genome FASTA records to avoid re-reading huge genome files
# on every Extend/Blast extraction call.
_GENOME_FASTA_CACHE: dict[str, tuple[float, dict[str, SeqRecord]]] = {}
_GENOME_FASTA_CACHE_LOCK = threading.Lock()

def check_cmd_in_path(cmd):
    cmd_path = shutil.which(cmd)
    if cmd_path:
        #logging.info(f"{cmd} is available at: {cmd_path}")
        return True
    else:
        #logging.error(f"{cmd} is not available in the system's PATH")
        # Raise error if tool is not available
        #raise FileNotFoundError
        return False


def get_original_file_path():
    """Package directory containing bundled BLAST binaries (teviewer/backends/)."""
    if hasattr(sys, '_MEIPASS'):
        return os.path.join(sys._MEIPASS, 'teviewer', 'backends')
    return os.path.dirname(os.path.abspath(__file__))


def sanitize_name(name):
    for char, replacement in special_character.items():
        name = name.replace(char, replacement)
    return name


def _blast_db_lock(key: str) -> threading.Lock:
    with _blast_db_locks_guard:
        if key not in _blast_db_locks:
            _blast_db_locks[key] = threading.Lock()
        return _blast_db_locks[key]


def blast_database_basename(genome_file: str | os.PathLike) -> str:
    """Safe BLAST database basename (``#`` and path separators break BLAST index lookup)."""
    return sanitize_name(os.path.basename(str(genome_file)))


_FLEXIDOT_UNSAFE_ID_CHARS = frozenset("/\\")


def prepare_flexidot_input_fasta(input_file: str | os.PathLike, work_dir: str | os.PathLike) -> tuple[str, str | None]:
    """
    Flexidot names PNGs after sequence IDs; ``/`` in headers becomes invalid paths.

    Returns ``(path_to_run, temp_path_or_none)``. Caller should delete the temp file when done.
    """
    input_path = os.path.abspath(str(input_file))
    work = os.path.abspath(str(work_dir))
    os.makedirs(work, exist_ok=True)

    records = list(SeqIO.parse(input_path, "fasta"))
    if not records:
        raise RuntimeError(f"No sequences in {os.path.basename(input_path)}.")

    needs_copy = any(any(ch in rec.id for ch in _FLEXIDOT_UNSAFE_ID_CHARS) for rec in records)
    if not needs_copy:
        return input_path, None

    seen: dict[str, int] = {}
    for record in records:
        safe = sanitize_name(record.id)
        count = seen.get(safe, 0)
        seen[safe] = count + 1
        if count:
            safe = f"{safe}_{count + 1}"
        record.id = safe
        record.description = safe

    stem = sanitize_name(os.path.splitext(os.path.basename(input_path))[0]) or "flexidot_input"
    temp_path = os.path.join(work, f".teviewer_flexidot_{stem}.fa")
    SeqIO.write(records, temp_path, "fasta")
    return temp_path, temp_path


def _makeblastdb_executable(
    blast_dir: str,
    os_type: str,
    *,
    use_system_blast: bool,
) -> str:
    if use_system_blast and check_cmd_in_path("makeblastdb"):
        return "makeblastdb"
    if os_type == "Linux":
        return os.path.join(blast_dir, "makeblastdb_linux")
    if os_type == "Darwin":
        return os.path.join(blast_dir, "makeblastdb_mac")
    return os.path.join(blast_dir, "makeblastdb.exe")


def check_database(
        genome_file,
        output_dir=None,
        os_type='Darwin',
        use_system_blast=False
):
    #script_dir = get_original_file_path()
    #blast_dir = os.path.join(script_dir, 'blast')

    script_dir = get_original_file_path()

    blast_dir = os.path.join(script_dir, "blast")
    # Define database path (sanitize basename: ``#`` breaks BLAST index files)
    genome_n = blast_database_basename(genome_file)
    if output_dir is None:
        database_path = os.path.join(os.path.dirname(str(genome_file)), genome_n)
    else:
        database_path = os.path.join(output_dir, genome_n)

    database_file = f'{database_path}.nin'
    lock_key = os.path.abspath(database_path)

    with _blast_db_lock(lock_key):
        if not os.path.isfile(database_file):
            try:

                makeblastdb_cmd = [
                    _makeblastdb_executable(
                        blast_dir, os_type, use_system_blast=use_system_blast
                    ),
                    "-in",
                    str(genome_file),
                    "-dbtype",
                    "nucl",
                    "-out",
                    str(database_path),
                ]

                run_logged_process(
                    makeblastdb_cmd,
                    tool="makeblastdb",
                    inputs=[genome_file],
                    cwd=os.path.dirname(database_path),
                    check=True,
                    output_path=database_path,
                )

            except FileNotFoundError:
                logging.error(
                    f"'makeblastdb' command not found. Please ensure 'makeblastdb' is installed correctly.\n"
                    f'{traceback.format_exc()}'
                )
                return 'makeblastdb_not_found'

            except subprocess.CalledProcessError as e:
                logging.error(
                    f'makeblastdb failed with exit code {e.returncode} \n'
                    f'{traceback.format_exc()}'
                )
                logging.error(e.stderr)
                return 'makeblastdb_got_error'
        return database_path


def blast(
    input_file,
    blast_database,
    blast_out_dir,
    e_value=1e-40,
    bed_file=False,
    self_blast=False,
    word_size=11,
    gap_open=5,
    gap_extend=2,
    penalty=-3,
    percent_identity=0.0,
    query_coverage_hsp=0.0,
    program="blastn",
    os_type='Darwin',
    num_threads=1,
    use_system_blast=False,
    output_label=None,
):
    from teviewer.job_options import normalize_blast_evalue

    program = str(program).lower()
    if program not in {"blastn", "blastp", "blastx"}:
        raise ValueError(f"Unsupported BLAST program: {program}")
    e_value = normalize_blast_evalue(e_value)
    script_dir = get_original_file_path()
    blast_dir = os.path.join(script_dir, 'blast')

    base_name = os.path.basename(input_file)
    label_suffix = f"_{output_label}" if output_label else ""

    output_tag = {"blastn": "blast", "blastp": "blastp", "blastx": "blastx"}[program]
    # Define BLAST output without header.
    blast_out_file = os.path.join(
        blast_out_dir, f'{base_name}{label_suffix}_{output_tag}_no_header.b'
    )
    # Define BLAST output with header.
    blast_out_file_header = os.path.join(
        blast_out_dir, f'{base_name}{label_suffix}_{output_tag}.txt'
    )

    # Use the requested system executable when enabled, otherwise the bundled copy.
    if use_system_blast and check_cmd_in_path(program):
        blast_executable = program
    else:
        if os_type == "Linux":
            blast_executable = os.path.join(blast_dir, f"{program}_linux")
        elif os_type == "Darwin":
            blast_executable = os.path.join(blast_dir, f"{program}_mac")
        else:
            blast_executable = os.path.join(blast_dir, f"{program}.exe")

    # -outfmt 6 use "\t" as deliminator. -outfmt 10 use "," as deliminator
    blast_cmd = [
        str(blast_executable),
        "-query", input_file,
        "-db", blast_database,
        "-max_target_seqs", str(10000),
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore",
        "-evalue", str(e_value),
        "-out", blast_out_file,
        "-num_threads", str(num_threads)
    ]

    if self_blast:
        # Add more parameters for self-blast
        blast_cmd.extend([
            "-word_size", str(max(4, min(25, int(word_size)))),
            "-gapopen", "5",
            "-gapextend", "2",
            "-reward", "2",
            "-penalty", "-3"
        ])
    elif program == "blastn":
        blast_cmd.extend([
            "-word_size", str(max(4, int(word_size))),
            "-gapopen", str(max(0, int(gap_open))),
            "-gapextend", str(max(0, int(gap_extend))),
            "-penalty", str(min(0, int(penalty))),
            "-perc_identity", str(max(0.0, min(100.0, float(percent_identity)))),
            "-qcov_hsp_perc", str(max(0.0, min(100.0, float(query_coverage_hsp)))),
        ])
    else:
        blast_cmd.extend([
            "-word_size", str(max(2, int(word_size))),
            "-gapopen", str(max(0, int(gap_open))),
            "-gapextend", str(max(0, int(gap_extend))),
            "-qcov_hsp_perc", str(max(0.0, min(100.0, float(query_coverage_hsp)))),
        ])

    try:
        # Run the BLAST command
        run_logged_process(
            blast_cmd,
            tool=program,
            inputs=[input_file, blast_database],
            cwd=blast_out_dir,
            check=True,
            output_path=blast_out_file,
        )

    except FileNotFoundError:
        logging.error(
            f"'{program}' command not found. Please ensure '{program}' is installed correctly.\n"
            f'{traceback.format_exc()}'
        )
        return ('blastn_not_found' if program == 'blastn' else 'blastp_not_found'), False

    except subprocess.CalledProcessError as e:
        logging.error(f'An error occurred during BLAST: \n {traceback.format_exc()}')
        logging.error(f'\nBLAST failed with error code {e.returncode}')
        # Print the error generated from the BLAST itself
        logging.error(e.stderr)
        return ('blastn_got_error' if program == 'blastn' else 'blastp_got_error'), False

    # Check if the blast hit number is 0
    if os.path.isfile(blast_out_file) and os.path.getsize(blast_out_file) == 0:
        return ('blast_n_zero' if program == 'blastn' else 'blastp_zero_hits'), False

    if not os.path.isfile(blast_out_file):
        logging.error(f"BLAST produced no output file: {blast_out_file}")
        return ('blastn_got_error' if program == 'blastn' else 'blastp_got_error'), False

    try:
        # Always produce one canonical BLAST table with header.
        header = [
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
            'sstrand',
            'evalue',
            'bitscore',
        ]
        # TE-Aid expects V1–V13 column names in a comma-separated table (TEtrimmer format).
        if not bed_file:
            col_n = [
                "V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13",
            ]
            blast_out_teaid_file = f"{blast_out_file}_TEAid.b"
            blast_data_teaid_df = pd.read_csv(
                blast_out_file, delimiter="\t", header=None, names=col_n
            )
            blast_data_teaid_df.to_csv(blast_out_teaid_file, sep=",", index=False)
            try:
                os.remove(blast_out_file)
            except OSError:
                pass
            return blast_out_teaid_file, False

        blast_out_file_df = pd.read_csv(
            blast_out_file, delimiter="\t", header=None, names=header
        )
        blast_out_file_df.to_csv(blast_out_file_header, sep="\t", index=False)
        if program in ("blastp", "blastx"):
            # Protein hits: no strand/BED, just the tabular text output.
            try:
                os.remove(blast_out_file)
            except OSError:
                pass
            return False, blast_out_file_header

        blast_out_bed_file = f'{blast_out_file}.bed'
        blast_data_bed_df = pd.read_csv(blast_out_file, delimiter='\t', header=None)

        # Define a function to determine the strand and adjust columns accordingly
        def process_row(row):
            if 'plus' in row[10]:
                return [row[1], row[8], row[9], row[0], 0, '+']
            else:
                return [row[1], row[9], row[8], row[0], 0, '-']

        # Apply the function to each row and create a new DataFrame
        blast_data_bed_df_processed = blast_data_bed_df.apply(
            process_row, axis=1, result_type='expand'
        )

        # Drop duplicate rows
        blast_data_bed_df_processed_unique = (
            blast_data_bed_df_processed.drop_duplicates()
        )

        # Save bed file
        blast_data_bed_df_processed_unique.to_csv(
            blast_out_bed_file, sep='\t', index=False, header=False
        )

        try:
            os.remove(blast_out_file)
        except OSError:
            pass
        return blast_out_bed_file, blast_out_file_header

    except Exception:
        logging.error(
            f'An error occurred while processing BLAST output files:\n{traceback.format_exc()}'
        )
        return False, False


def select_top_longest_lines(lines, n):
    return sorted(lines, key=lambda line: int(line[2]) - int(line[1]), reverse=True)[:n]


def select_random_lines(n, remaining_lines):
    random.shuffle(remaining_lines)
    return remaining_lines[:n]


def process_bed_lines(
    input_file, output_dir, max_lines=100, top_longest_lines_count=70
):
    """
    Select desired line numbers from bed file.

    :param output_dir: str, the absolute path of output folder directory.
    :param max_lines: num default=100, the maximum line number for bed file to keep.
    :param top_longest_lines_count: num default=100 smaller or equal to max_lines.
    When the bed file line number is excess than the max_lines, sort bed file by sequence length and choose
    "top_longest_lines_count" lines. When this number is smaller than the max_lines, randomly choose the rest
    number of lines from the bed file.

    :return: the absolute selected bed file path.
    """
    # check if top_longest_lines_count is equal or smaller than max_lines.
    if top_longest_lines_count > max_lines:
        top_longest_lines_count = max_lines

    with open(input_file, 'r') as file:
        lines = [line.strip().split('\t') for line in file]

    bed_out_filter_file = os.path.join(
        output_dir, f'{os.path.basename(input_file)}_f.bed'
    )

    if len(lines) > max_lines:
        top_longest_lines = select_top_longest_lines(lines, top_longest_lines_count)

        # Eliminate selected lines
        remaining_lines = [line for line in lines if line not in top_longest_lines]

        # Randomly choose lines
        random_lines = select_random_lines(
            max_lines - top_longest_lines_count, remaining_lines
        )
        selected_lines = top_longest_lines + random_lines

        # Write the selected lines to the output file.
        with open(bed_out_filter_file, 'w') as file:
            for line in selected_lines:
                file.write('\t'.join(line) + '\n')
    else:
        # Copy the original file to the output directory with the new name to keep the file name consistency.
        shutil.copy(input_file, bed_out_filter_file)

    return bed_out_filter_file


# Generate bed file based on fasta header，this is used for the extend function
# The fasta header could look like
"""
> scaffold_1:343622-349068(+)
> scaffold_1:346171-347467(+)
> scaffold_1:346385-347665(+)
> scaffold_1:346200-347196(+)
"""


def fasta_header_to_bed(input_file, output_file):
    unique_lines = set()

    with open(input_file, 'r') as fasta, open(output_file, 'w') as bed:
        for line in fasta:
            if line.startswith('>'):
                header = line.strip().lstrip('>')  # Remove '>'
                parts = header.rsplit(':', 2)  # Split from the right side
                scaffold = parts[0]
                range_strand = parts[1].split('(')
                range_part = range_strand[0]
                strand = range_strand[1].split(')')[0]
                start, end = range_part.split('-')
                bed_line = f'{scaffold}\t{start}\t{end}\tTEtrimmer\t0\t{strand}\n'

                # Add the line to the set if it's not already present
                if bed_line not in unique_lines:
                    bed.write(bed_line)
                    unique_lines.add(bed_line)

    return output_file


# Avoid to use bedtools slop and getfasta to make it compatible with Windows system
def extend_bed_regions(
    bed_file, left_extension, right_extension, chrom_sizes, output_bed
):
    with open(bed_file, 'r') as infile, open(output_bed, 'w') as outfile:
        for line in infile:
            parts = line.strip().split()
            chrom, start, end, strand = parts[0], int(parts[1]), int(parts[2]), parts[5]

            # Adjust extensions based on strand
            if strand == '+':
                new_start = max(0, start - left_extension)
                new_end = min(chrom_sizes[chrom], end + right_extension)
            else:  # For anti sense strand
                new_start = max(
                    0, start - right_extension
                )  # Extend "start" less, as it's the "end"
                new_end = min(
                    chrom_sizes[chrom], end + left_extension
                )  # Extend "end" more, as it's the "start"

            # Write the extended region to the output BED file
            outfile.write(f'{chrom}\t{new_start}\t{new_end}\tTEtrimmer\t0\t{strand}\n')
    return output_bed


def read_bed(bed_file):
    with open(bed_file, 'r') as file:
        for line in file:
            parts = line.strip().split('\t')
            if len(parts) < 6:
                continue
            yield {
                'chrom': parts[0],
                'start': int(parts[1]),
                'end': int(parts[2]),
                'strand': parts[5],
            }


def _load_genome_records_cached(genome_fasta: str) -> dict[str, SeqRecord]:
    """Load genome FASTA into memory once per file mtime."""
    key = os.path.abspath(genome_fasta)
    try:
        mtime = os.path.getmtime(key)
    except OSError:
        mtime = -1.0
    with _GENOME_FASTA_CACHE_LOCK:
        cached = _GENOME_FASTA_CACHE.get(key)
        if cached is not None and cached[0] == mtime:
            return cached[1]
        records = SeqIO.to_dict(SeqIO.parse(genome_fasta, 'fasta'))
        _GENOME_FASTA_CACHE[key] = (mtime, records)
        return records


def extract_fasta_from_bed(genome_fasta, bed_file, output_fasta):
    genome = _load_genome_records_cached(genome_fasta)
    with open(output_fasta, 'w') as out_fasta:
        for region in read_bed(bed_file):
            chrom = region['chrom']
            start = region['start']
            end = region['end']
            strand = region['strand']
            if chrom in genome:
                sequence = genome[chrom].seq[start:end]
                if strand == '-':
                    sequence = (
                        sequence.reverse_complement()
                    )  # Reverse complement if on negative strand
                seq_record = SeqRecord(
                    sequence, id=f'{chrom}:{start}-{end}({strand})', description=''
                )
                SeqIO.write(seq_record, out_fasta, 'fasta')
    return output_fasta


def flexidot_executable() -> str | None:
    """Resolve the Flexidot CLI (PATH, then common pip user ``bin`` locations)."""
    import glob

    found = shutil.which("flexidot")
    if found:
        return found
    candidates: list[str] = []
    for pattern in (
        os.path.expanduser("~/Library/Python/*/bin/flexidot"),
        os.path.expanduser("~/.local/bin/flexidot"),
    ):
        candidates.extend(glob.glob(pattern))
    for base in (getattr(sys, "prefix", ""), os.path.expanduser("~/.local")):
        if base:
            candidates.append(os.path.join(base, "bin", "flexidot"))
    try:
        import site

        user_base = site.getuserbase()
        if user_base:
            candidates.append(os.path.join(user_base, "bin", "flexidot"))
    except Exception:  # noqa: BLE001
        pass
    for candidate in candidates:
        if os.path.isfile(candidate) and os.access(candidate, os.X_OK):
            return candidate
    return None


_flexidot_help_cache: dict[str, str] = {}


def _flexidot_help_text(exe: str) -> str:
    if exe in _flexidot_help_cache:
        return _flexidot_help_cache[exe]
    text = ""
    try:
        cmd = [exe, "-h"]
        result = run_logged_process(cmd, tool="Flexidot", timeout=20)
        text = ((result.stdout or "") + (result.stderr or "")).lower()
    except Exception:  # noqa: BLE001
        text = ""
    _flexidot_help_cache[exe] = text
    return text


def flexidot_supports_html(exe: str | None = None) -> bool:
    """True if the installed Flexidot lists ``html`` among the ``-f`` choices.

    The stock PyPI Flexidot only offers png/pdf/svg; the forked build adds html
    (interactive dotplot with Go-to). We look for html as a filetype choice (in
    the ``{...}`` list) rather than anywhere in the text, to avoid false matches.
    """
    if exe is None:
        exe = flexidot_executable()
    if not exe:
        return False
    help_text = _flexidot_help_text(exe)
    return any(tok in help_text for tok in ("{html", "html}", ",html", "html,"))


def flexidot_supports_teviewer_integration(exe: str | None = None) -> bool:
    """True if the forked Flexidot exposes the ``--teviewer-integration`` flag."""
    if exe is None:
        exe = flexidot_executable()
    if not exe:
        return False
    return "teviewer-integration" in _flexidot_help_text(exe)


def effective_flexidot_filetype(requested: str) -> str:
    """Map a requested Flexidot output format to one the installed build supports.

    ``html`` requires the fork; if it is unavailable, fall back to ``png`` so the
    dotplot still renders (and output-file discovery uses the right extension).
    """
    ft = str(requested).strip().lower()
    if ft not in {"html", "png", "pdf", "svg"}:
        ft = "html"
    if ft == "html" and not flexidot_supports_html():
        return "png"
    return ft


def run_flexidot(
    input_file,
    output_dir,
    wordsize=3,
    mode='0',
    *,
    filetype: str = "html",
    lcs_shading: bool = True,
    use_gff: bool = False,
    substitution_count: int = 0,
    wobble_conversion: bool = False,
    max_n_percentage: float = 10.0,
    include_reverse: bool = True,
    only_vs_first: bool = False,
    line_width: float = 1.0,
    gff_source_fastas: list[str | os.PathLike] | None = None,
    gff_output_stem: str | None = None,
    cancel=None,
):
    """
    Generates dotplots using FlexiDot.

    :param input_file: str, path to the input FASTA file
    :param output_dir: str, directory where results will be saved
    :param wordsize: int, k-mer length (-k)
    :param mode: str, plotting mode ('0'=self, '1'=paired, '2'=poly)
    :param filetype: str, Flexidot output format (-f): 'html' (interactive,
        default), 'png', 'pdf', or 'svg'
    :param lcs_shading: bool, add -x (LCS shading) for multi-sequence plots
    :param use_gff: bool, pass companion GFF files to Flexidot (-g / -G)
    :param gff_source_fastas: FASTA paths used to locate companion GFF3 files
    :param gff_output_stem: basename for generated ``*_flexidot.gff.colors.config``
    """
    out_dir = os.path.abspath(str(output_dir))
    os.makedirs(out_dir, exist_ok=True)
    run_path, temp_path = prepare_flexidot_input_fasta(input_file, out_dir)
    original_name = os.path.basename(os.path.abspath(str(input_file)))

    exe = flexidot_executable()
    if not exe:
        raise FileNotFoundError(
            "Flexidot command not found. Install with: pip install flexidot "
            "(ensure the Python scripts directory is on PATH, e.g. "
            "~/Library/Python/3.x/bin on macOS)."
        )

    logging.info("FlexiDot is running for %s...", original_name)

    flexidot_cmd = [
        exe,
        "-i",
        run_path,
        "-m",
        str(mode),
        "-k",
        str(wordsize),
        "-f",
        filetype,
        "--outdir",
        out_dir,
    ]
    if filetype == "html" and flexidot_supports_teviewer_integration(exe):
        # Show "Go to <seq> <start>-<end>" buttons in the dotplot that jump the
        # Annotation view; only available in the forked Flexidot build.
        flexidot_cmd.append("--teviewer-integration")
    if lcs_shading and str(mode) in {"1", "2"}:
        flexidot_cmd.append("-x")
    if substitution_count > 0:
        flexidot_cmd.extend(["-S", str(int(substitution_count))])
    if wobble_conversion:
        flexidot_cmd.extend(["-w", "--max_n", str(float(max_n_percentage))])
    if not include_reverse:
        flexidot_cmd.append("-r")
    if only_vs_first and str(mode) == "1":
        flexidot_cmd.append("-O")
    flexidot_cmd.extend(["-A", str(float(line_width))])

    temp_gff_paths: list[str] = []
    if use_gff:
        from pathlib import Path

        from teviewer.backends.flexidot_gff import prepare_flexidot_gff_inputs

        source_fastas = [
            Path(str(path)).expanduser().resolve()
            for path in (gff_source_fastas or [input_file])
        ]
        stem = gff_output_stem or sanitize_name(os.path.splitext(original_name)[0]) or "flexidot"
        prepared = prepare_flexidot_gff_inputs(
            source_fastas,
            Path(out_dir),
            output_stem=stem,
        )
        if prepared is None:
            logging.warning(
                "Flexidot GFF shading enabled but no companion GFF found for %s",
                original_name,
            )
        else:
            gff_paths, color_config, temp_gff_files = prepared
            flexidot_cmd.extend(["-g", *gff_paths])
            flexidot_cmd.extend(["-G", color_config])
            temp_gff_paths = [str(path) for path in temp_gff_files]

    env = os.environ.copy()
    exe_dir = os.path.dirname(exe)
    if exe_dir:
        env["PATH"] = exe_dir + os.pathsep + env.get("PATH", "")

    try:
        completed = run_subprocess_with_cancel(
            flexidot_cmd,
            cancel=cancel,
            env=env,
            check=True,
        )
        if completed.stderr.strip():
            logging.debug("FlexiDot stderr: %s", completed.stderr.strip())
        logging.info("FlexiDot finished successfully.")

    except JobCancelledError:
        raise
    except subprocess.CalledProcessError as e:
        detail = (e.stderr or e.stdout or str(e)).strip()
        logging.error("FlexiDot failed: %s", detail)
        hint = ""
        if temp_path is None and "/" in original_name:
            hint = (
                "\n\nIf the FASTA header contains '/', Flexidot cannot write output files; "
                "TEviewer should have sanitized IDs automatically — please report this case."
            )
        raise RuntimeError(
            f"Flexidot failed for {original_name} (exit {e.returncode}).\n\n{detail}{hint}"
        ) from e
    finally:
        if temp_path and os.path.isfile(temp_path):
            try:
                os.unlink(temp_path)
            except OSError:
                pass
        for gff_path in temp_gff_paths:
            if not os.path.basename(gff_path).startswith(".teviewer_flexidot_"):
                continue
            if os.path.isfile(gff_path):
                try:
                    os.unlink(gff_path)
                except OSError:
                    pass

    return out_dir


def align_sequences(input_file, output_file, num_thread=5):
    """
    Aligns FASTA sequences using MAFFT

    :param input_file: str, input file path
    :param output_dir: str, output file directory
    :return: absolute path for multiple sequence alignment file
    """
    input_file_n = os.path.basename(input_file)

    logging.info("Mafft is running......")

    # Construct the command as a list of strings
    mafft_cmd = ['mafft', 
                 '--nuc', 
                 '--retree', '1', 
                 '--thread', str(num_thread),
                 input_file
                 ]

    try:
        # Execute the command
        result = run_logged_process(
            mafft_cmd,
            tool="MAFFT",
            inputs=[input_file],
            cwd=os.path.dirname(os.path.abspath(output_file)),
            check=True,
            stdout_is_result=True,
            output_path=output_file,
        )

    except FileNotFoundError as e:
        logging.error(
            "'mafft' command not found. Please ensure 'mafft' is installed correctly."
        )
        raise Exception from e

    except subprocess.CalledProcessError as e:
        error_message = e.stderr
        logging.error(f'\nMAFFT failed for {input_file_n} with error code {e.returncode}\n{e.stdout}\n{e.stderr}\n')

        if (
            'Killed' in error_message
            and 'disttbfast' in error_message
            and 'memopt' in error_message
        ):
            logging.error(
                'It seems insufficient RAM was available for MAFFT multiple sequence alignment. Please assign more '
                'RAM or reduce the thread number to solve this problem.\n'
            )
        raise Exception from e

    # Write the output to the file
    with open(output_file, 'w') as f:
        f.write(result.stdout)

    logging.info("Mafft is finished.")

    return output_file


def check_cdd_index_files(directory):
    """
    Check if the cdd index files exist in the provided directory.
    :param directory: str, the directory to search for the cdd index files.
    :return: boolean, 'True' if all cdd index files exist, 'False' otherwise.
    """
    # Define the list of cdd index files
    cdd_files = [
        'cdd_profile.24.freq',
        'cdd_profile.24.aux',
        'cdd_profile.24.rps',
        'cdd_profile.24.loo',
        'cdd_profile.24.ptf',
        'cdd_profile.24.pot',
        'cdd_profile.24.pdb',
        'cdd_profile.24.pos',
    ]

    # Check if all index files exist, store missing files in a list
    missing_files = []
    for file in cdd_files:
        if not os.path.isfile(os.path.join(directory, file)):
            missing_files.append(file)

    # Return True if all files exist, False otherwise
    if not missing_files:
        return True
    else:
        return False
