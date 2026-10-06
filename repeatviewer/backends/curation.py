"""TE curation tools (ported from TEtrimmerGUI, standalone in RepeatViewer)."""

from __future__ import annotations

import logging
import os
import platform
import re
import shutil
import tempfile
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Iterable

import plotly.io as pio
from Bio import SeqIO

from repeatviewer.backends.blast_gff import blast_gff_output_path, write_blast_gff
from repeatviewer.tool_output_paths import write_genome_reference_link
from repeatviewer.backends.flexidot_gff import (
    persistent_flexidot_color_config_path,
    persistent_flexidot_legend_path,
    resolve_gff_paths_for_fastas,
    resolve_user_flexidot_color_config,
    write_flexidot_color_config,
    write_flexidot_gff_legend_png,
    _gff_feature_display_labels,
)
from repeatviewer.backends.crop_end_divergence import crop_end_div
from repeatviewer.backends.crop_end_gap import crop_end_gap
from repeatviewer.backends.gui_functions import (
    sanitize_name,
    blast as run_blast,
    check_database,
    effective_flexidot_filetype,
    extend_bed_regions,
    extract_fasta_from_bed,
    fasta_header_to_bed,
    process_bed_lines,
    run_flexidot,
)
from repeatviewer.backends.remove_gap import remove_gaps_with_similarity_check
from repeatviewer.backends.sequence_operations import (
    write_predicted_orfs_fasta,
    write_reverse_complement_fasta,
    write_translated_orfs_fasta,
)
from repeatviewer.backends.teaid_plotter import (
    TeaidInputError,
    con_generater,
    teaid_pan_genome_plotter,
    teaid_plotter,
)
from repeatviewer.backends.job_cancel import JobCancelledError, JobCancellation, check_cancelled

logger = logging.getLogger(__name__)
from repeatviewer.backends.thread_budget import analysis_threads_per_job
from repeatviewer.backends.tsd_search import TSD_identifier
from repeatviewer.job_options import (
    BLAST_PROGRAM_ONLINE,
    BlastOptions,
    ConsensusOptions,
    ExtendOptions,
    FlexidotOptions,
    HmmerOptions,
    MsaCleanOptions,
    TeaAidOptions,
    TsdSearchOptions,
    auto_flexidot_mode,
)
from repeatviewer.pan_genome import genome_output_suffix
from repeatviewer.backends.te_annotation_batch import TeAnnotationOptions, run_te_annotation_on_fasta
from repeatviewer.backends.separate_fasta import count_fasta_sequences
from repeatviewer.project import BlastDatabase, RepeatViewerProject
from repeatviewer.backends.database_prepare import blast_database_is_ready

ProgressCallback = Callable[[int, str], None]
# completed_in_file, total_in_file, label — completed=-1 means label-only update.
SequenceProgressCallback = Callable[[int, int, str], None]

_PER_SEQUENCE_FUNCS = frozenset({"te_annotation"})


def _fasta_sequence_units(path: Path) -> int:
    try:
        return max(1, int(count_fasta_sequences(path)))
    except Exception:  # noqa: BLE001
        return 1


_DISPLAY_NAMES = {
    "teaid": "TE-Aid",
    "blast": "BLAST",
    "hmmer": "HMMER",
    "extend": "MSA Extend",
    "msa_clean": "MSA Clean",
    "msa_clean_gappy": "MSA remove gappy columns",
    "msa_clean_div": "Clean row with divergence",
    "msa_clean_gap": "Clean row with gap",
    "consensus_generation": "Consensus",
    "flexidot": "Flexidot",
    "te_annotation": "Sequence Annotation",
    "reverse_complement": "Reverse complement",
    "orf_prediction": "ORF prediction",
    "dna_translation": "DNA translation",
}


def _snapshot_paths(root: Path) -> set[Path]:
    root = root.expanduser().resolve()
    if not root.is_dir():
        return set()
    return {p.resolve() for p in root.rglob("*")} | {root}


_flexidot_dir_locks: dict[Path, threading.Lock] = {}
_flexidot_locks_guard = threading.Lock()


def _flexidot_dir_lock(folder: Path) -> threading.Lock:
    key = folder.expanduser().resolve()
    with _flexidot_locks_guard:
        if key not in _flexidot_dir_locks:
            _flexidot_dir_locks[key] = threading.Lock()
        return _flexidot_dir_locks[key]


def cleanup_teaid_side_files(
    output_dir: Path,
    input_file: Path,
    keep_html: Path,
) -> None:
    """
    Remove intermediate TE-Aid files for one input only.

    Keeps the input FASTA, this run's HTML report, and Blast-tool sidecars
    (``*_b.txt``, ``*_b.gff``, ``*_b.txt.fa``, etc.). Safe when several TE-Aid
    jobs run in parallel in the same folder.
    """
    output_dir = output_dir.expanduser().resolve()
    input_file = input_file.expanduser().resolve()
    keep_html = keep_html.expanduser().resolve()
    input_name = input_file.name
    consensus_prefix = f"{input_name}_0.3_co.fa"

    to_remove: list[Path] = []
    if not output_dir.is_dir():
        return

    for entry in output_dir.iterdir():
        resolved = entry.resolve()
        if resolved in {keep_html, input_file}:
            continue
        name = entry.name
        if name.endswith("_TEAid.html"):
            continue
        if entry.is_dir() and name.endswith("_self_blast"):
            prefix = name[: -len("_self_blast")]
            if prefix == input_name or prefix.startswith(f"{input_name}_"):
                to_remove.append(resolved)
            continue
        if not (entry.is_file() or entry.is_symlink()):
            continue
        # Consensus intermediates from this TE-Aid run only — never wipe Blast
        # outputs that share the query name prefix ({input}_b.txt / .gff / .txt.fa).
        if name.startswith(consensus_prefix):
            to_remove.append(resolved)

    for path in sorted(to_remove, key=lambda p: len(p.parts), reverse=True):
        try:
            if path.is_dir():
                shutil.rmtree(path, ignore_errors=True)
            else:
                path.unlink(missing_ok=True)
        except OSError:
            pass


def _output_paths_in_dir(folder: Path, ext: str) -> set[Path]:
    folder = folder.expanduser().resolve()
    return {p.resolve() for p in folder.glob(f"*.{ext}") if p.is_file()}


def _png_paths_in_dir(folder: Path) -> set[Path]:
    return _output_paths_in_dir(folder, "png")


def flexidot_output_stem(seq_files: list[Path], genome_label: str = "") -> str:
    """Output basename: ``{first.name}_flexidot`` or ``{first.name}_flexidot_multi``."""
    from repeatviewer.tool_output_paths import kept_input_name

    base = kept_input_name(seq_files[0].expanduser().resolve())
    suffix = genome_output_suffix(genome_label)
    if len(seq_files) >= 2:
        return f"{base}_flexidot_multi{suffix}"
    return f"{base}_flexidot{suffix}"


def is_flexidot_html_file(path: Path) -> bool:
    """True for a Flexidot interactive dotplot HTML output (``*_flexidot*.html``)."""
    name = path.name.lower()
    if not name.endswith((".html", ".htm")):
        return False
    return "_flexidot" in name


_FLEXIDOT_EXTRA_FILE_RE = re.compile(
    r".+_flexidot(?:_multi)?.*_(\d+)\.(?:png|html)$",
    re.IGNORECASE,
)


def is_flexidot_extra_plot_file(path: Path) -> bool:
    """True for secondary Flexidot outputs (``*_flexidot_2.png``/``.html``, etc.)."""
    match = _FLEXIDOT_EXTRA_FILE_RE.match(path.name)
    return match is not None and int(match.group(1)) >= 2


def _remove_flexidot_extra_outputs(output_dir: Path, base: str, ext: str) -> None:
    """Delete numbered Flexidot extras left from older runs."""
    output_dir = output_dir.expanduser().resolve()
    for index in range(2, 20):
        extra = output_dir / f"{base}_{index}.{ext}"
        try:
            extra.unlink(missing_ok=True)
        except OSError:
            pass


def _genome_work_items(
    project: RepeatViewerProject,
    func_name: str,
    extend: ExtendOptions | None = None,
    teaid: TeaAidOptions | None = None,
) -> list[tuple[str, Path]]:
    """Return ``(label, fasta_path)`` pairs for one batch job."""
    if func_name == "teaid" and teaid is not None and teaid.use_pan_genome:
        entries = [
            (str(name).strip(), Path(path).expanduser().resolve())
            for name, path in (teaid.pan_genomes or [])
            if str(name).strip() and str(path).strip()
        ]
        if len(entries) < 2:
            raise RuntimeError(
                "TE-Aid pangenome mode needs at least two genomes "
                "(Mode → Pangenome)."
            )
        missing = [name for name, path in entries if not path.is_file()]
        if missing:
            raise RuntimeError(
                "TE-Aid pangenome genome file not found: " + ", ".join(missing)
            )
        return entries

    if not project.is_pan_genome_mode:
        override_path = None
        if func_name == "extend" and extend and extend.genome_path:
            override_path = extend.genome_path
        elif func_name == "teaid" and teaid and teaid.genome_path:
            override_path = teaid.genome_path
        if override_path:
            path = Path(override_path).expanduser().resolve()
            if path.is_file():
                return [("", path)]
            raise RuntimeError(f"Genome file not found: {path}")
        return [
            (label, path)
            for label, path in project.genome_targets()
        ]

    if func_name == "extend":
        extend_opts = extend or ExtendOptions()
        name = extend_opts.pan_genome_name
        if not name and project.pan_genomes:
            name = project.pan_genomes[0].name
        if not name:
            raise RuntimeError("Pan-genome mode: select a genome name for Extend.")
        entry = project.pan_genome_by_name(name)
        if entry is None:
            raise RuntimeError(f"Unknown pan-genome name: {name}")
        return [(entry.name, entry.path.expanduser().resolve())]

    return [
        (entry.name, entry.path.expanduser().resolve())
        for entry in project.pan_genomes
    ]


def _find_blast_database(
    project: RepeatViewerProject,
    name: str,
    *,
    molecule_type: str | None = None,
) -> BlastDatabase:
    for db in project.blast_databases:
        if db.name == name:
            if not blast_database_is_ready(db):
                raise RuntimeError(
                    f"BLAST database '{name}' is missing its index files "
                    f"(expected at {db.db_prefix_path}). Rebuild it in "
                    "Add/Manage databases."
                )
            if molecule_type is not None and db.molecule_type != molecule_type:
                expected = "protein" if molecule_type == "prot" else "nucleotide"
                raise RuntimeError(f"BLAST database '{name}' is not a {expected} database.")
            return db
    raise RuntimeError(f"BLAST database '{name}' is no longer registered.")


def _prepare_genome_resources(
    func_name: str,
    genome_items: list[tuple[str, Path]],
    project: RepeatViewerProject,
    os_type: str,
    cancel: JobCancellation | None = None,
) -> tuple[dict[Path, str], dict[Path, dict[str, int]]]:
    db_cache: dict[Path, str] = {}
    chrom_cache: dict[Path, dict[str, int]] = {}
    for _label, genome_path in genome_items:
        if func_name == "blast":
            genome_db = check_database(
                str(genome_path),
                os_type=os_type,
                use_system_blast=project.use_system_blast,
                cancel=cancel,
            )
            if genome_db in {"makeblastdb_not_found", "makeblastdb_got_error", None, False}:
                raise RuntimeError(
                    f"Cannot prepare genome BLAST database for {genome_path.name}: {genome_db}"
                )
            db_cache[genome_path] = genome_db
        if func_name == "extend":
            chrom_cache[genome_path] = genome_chrom_sizes(genome_path)
    return db_cache, chrom_cache


def _flexidot_outputs_for_run(
    output_dir: Path,
    before: set[Path],
    seq_files: list[Path],
    ext: str,
) -> list[Path]:
    """New Flexidot output files of ``ext`` from one run (optionally filtered to one input stem)."""
    after = _output_paths_in_dir(output_dir, ext)
    new_files = sorted(
        after - before,
        key=lambda p: p.stat().st_size,
        reverse=True,
    )
    if len(seq_files) == 1:
        # Flexidot may name intermediates from stem or full basename.
        markers = {seq_files[0].stem, seq_files[0].name}
        owned = [
            p
            for p in new_files
            if any(marker and marker in p.name for marker in markers)
        ]
        if owned:
            return owned
    return new_files


def cleanup_flexidot_intermediates(output_dir: Path, keep_paths: Iterable[Path]) -> None:
    """Remove Flexidot scratch files after plots and sidecars are finalized."""
    output_dir = output_dir.expanduser().resolve()
    keep = {p.expanduser().resolve() for p in keep_paths}
    if not output_dir.is_dir():
        return
    for entry in output_dir.iterdir():
        name = entry.name
        if not (
            name.startswith("flexidot_output")
            or name.startswith(".repeatviewer_flexidot_")
        ):
            continue
        resolved = entry.resolve()
        if resolved in keep:
            continue
        try:
            if entry.is_dir():
                shutil.rmtree(resolved, ignore_errors=True)
            else:
                entry.unlink(missing_ok=True)
        except OSError:
            pass


def cleanup_blast_scratch_files(
    output_dir: Path,
    seq_path: Path,
    *,
    genome_label: str = "",
) -> None:
    """Remove BLAST intermediates that are not needed after a successful run."""
    output_dir = output_dir.expanduser().resolve()
    label_suffix = genome_output_suffix(genome_label)
    base_name = seq_path.name
    from repeatviewer.backends.blast_output_names import LOCAL_BLAST_PROGRAM_TAGS

    legacy_tags = ("b", "blast", "blastp", "blastx")
    for tag in (*LOCAL_BLAST_PROGRAM_TAGS.values(), *legacy_tags):
        for name in (
            f"{base_name}{label_suffix}_{tag}_no_header.b",
            f"{base_name}{label_suffix}_{tag}_no_header.b.bed",
        ):
            try:
                (output_dir / name).unlink(missing_ok=True)
            except OSError:
                pass


def _new_flexidot_legend_pngs(output_dir: Path, before: set[Path]) -> list[Path]:
    """Flexidot-native GFF legend images created during a run."""
    after = _png_paths_in_dir(output_dir)
    return sorted(
        (
            path
            for path in (after - before)
            if path.is_file() and "legend" in path.name.lower() and "gff" in path.name.lower()
        ),
        key=lambda p: p.stat().st_size,
        reverse=True,
    )


def _finalize_flexidot_gff_sidecars(
    output_dir: Path,
    output_stem: str,
    before: set[Path],
    *,
    gff_paths: list[Path] | None = None,
) -> tuple[list[str], list[Path]]:
    """Keep editable color config and a legend PNG for GFF-shaded Flexidot plots."""
    output_dir = output_dir.expanduser().resolve()
    created: list[str] = []
    keep: list[Path] = []
    gff_paths = [path.expanduser().resolve() for path in (gff_paths or []) if path]

    persist_config = persistent_flexidot_color_config_path(output_dir, output_stem)
    legend_dest = persistent_flexidot_legend_path(output_dir, output_stem)
    temp_config = output_dir / f".repeatviewer_flexidot_{output_stem}.gff.colors.config"
    user_config = resolve_user_flexidot_color_config(output_dir, output_stem)

    if temp_config.is_file():
        shutil.copy2(temp_config, persist_config)
        try:
            temp_config.unlink(missing_ok=True)
        except OSError:
            pass
    elif user_config is not None and user_config.resolve() != persist_config.resolve():
        shutil.copy2(user_config, persist_config)
    elif gff_paths and not persist_config.is_file():
        write_flexidot_color_config(gff_paths, persist_config)

    if persist_config.is_file():
        keep.append(persist_config)
        created.append(str(persist_config.resolve()))
        legacy_name = f"{output_stem}_flexidot.gff.colors.config"
        if legacy_name != persist_config.name:
            try:
                (output_dir / legacy_name).unlink(missing_ok=True)
            except OSError:
                pass
        legacy_legend = f"{output_stem}_flexidot_legend.png"
        if legacy_legend != legend_dest.name:
            try:
                (output_dir / legacy_legend).unlink(missing_ok=True)
            except OSError:
                pass

    native_legends = _new_flexidot_legend_pngs(output_dir, before)
    if native_legends:
        if legend_dest.exists():
            legend_dest.unlink()
        native_legends[0].rename(legend_dest)
    elif persist_config.is_file():
        try:
            labels = _gff_feature_display_labels(gff_paths)
            write_flexidot_gff_legend_png(
                persist_config,
                legend_dest,
                display_labels=labels,
            )
        except (OSError, ValueError) as exc:
            logger.warning("Flexidot legend PNG skipped: %s", exc)

    if legend_dest.is_file():
        keep.append(legend_dest)
        created.append(str(legend_dest.resolve()))
    return created, keep


def finalize_flexidot_outputs(
    output_dir: Path,
    seq_files: list[Path],
    before: set[Path],
    *,
    ext: str = "html",
    use_gff: bool = False,
    gff_source_fastas: list[Path] | None = None,
) -> list[str]:
    """Keep the primary Flexidot output as ``{name}_flexidot.{ext}``; drop extra panels.

    ``before`` must be the snapshot taken with the same ``ext`` (e.g. via
    ``_output_paths_in_dir(out_dir, ext)``) right before the run.
    """
    output_dir = output_dir.expanduser().resolve()
    new_files = _flexidot_outputs_for_run(output_dir, before, seq_files, ext)
    if not new_files:
        return []

    base = flexidot_output_stem(seq_files)
    _remove_flexidot_extra_outputs(output_dir, base, ext)

    src = new_files[0]
    from repeatviewer.backends.output_paths import next_numbered_path

    dest = next_numbered_path(output_dir / f"{base}.{ext}")
    if not src.is_file():
        cleanup_flexidot_intermediates(output_dir, [])
        return []

    if src.resolve() != dest.resolve():
        src.rename(dest)
    for extra in new_files[1:]:
        try:
            if extra.is_file():
                extra.unlink(missing_ok=True)
        except OSError:
            pass

    keep_paths = [dest]
    outputs = [str(dest.resolve())]
    if use_gff:
        gff_paths = resolve_gff_paths_for_fastas(gff_source_fastas or seq_files)
        # Only keep color-config / legend when a companion GFF was actually used.
        # Otherwise those files look like Flexidot outputs but have no meaning.
        if gff_paths:
            sidecar_paths, sidecar_keep = _finalize_flexidot_gff_sidecars(
                output_dir,
                base,
                before,
                gff_paths=gff_paths,
            )
            outputs.extend(sidecar_paths)
            keep_paths.extend(sidecar_keep)

    cleanup_flexidot_intermediates(output_dir, keep_paths)
    return outputs


def effective_batch_workers(requested: int, num_sequences: int) -> int:
    """Clamp parallel worker count to sensible limits."""
    if num_sequences <= 1 or requested <= 1:
        return 1
    cpus = os.cpu_count() or 4
    return max(1, min(requested, num_sequences, cpus))


def genome_chrom_sizes(genome_path: Path) -> dict[str, int]:
    chrom_sizes = {}
    for rec in SeqIO.parse(str(genome_path), "fasta"):
        chrom_sizes[rec.id] = len(rec.seq)
    if not chrom_sizes:
        raise RuntimeError("Genome FASTA has no sequences.")
    return chrom_sizes


def run_flexidot_combined(
    seq_files: list[Path],
    *,
    flexidot_options: FlexidotOptions | None = None,
    cancel: JobCancellation | None = None,
    genome_label: str = "",
) -> list[str]:
    """Flexidot for 2+ inputs: mode 1 (paired) or 2 (poly)."""
    if len(seq_files) < 2:
        raise ValueError("Combined Flexidot requires at least two sequences.")
    flex_opts = flexidot_options or FlexidotOptions()
    filetype = effective_flexidot_filetype(flex_opts.output_format)
    out_dir = seq_files[0].parent.resolve()
    combined = out_dir / f"_flexidot_{len(seq_files)}seq.fa"
    records = []
    for path in seq_files:
        records.extend(SeqIO.parse(str(path.expanduser().resolve()), "fasta"))
    if not records:
        raise RuntimeError("No sequences found in selected FASTA files.")
    # Explicit modes win; otherwise self / pairwise / all-vs-all from count.
    mode: str | list[str] = (
        list(flex_opts.modes) if flex_opts.modes else auto_flexidot_mode(len(records))
    )
    seen: dict[str, int] = {}
    for record in records:
        safe = sanitize_name(record.id)
        count = seen.get(safe, 0)
        seen[safe] = count + 1
        if count:
            safe = f"{safe}_{count + 1}"
        record.id = safe
        record.description = safe
    SeqIO.write(records, str(combined), "fasta")
    with _flexidot_dir_lock(out_dir):
        before = _output_paths_in_dir(out_dir, filetype)
        try:
            run_flexidot(
                str(combined),
                str(out_dir),
                wordsize=flex_opts.wordsize,
                mode=mode,
                filetype=filetype,
                sequence_type=flex_opts.sequence_type,
                lcs_shading=flex_opts.lcs_shading,
                use_gff=flex_opts.use_gff,
                substitution_count=flex_opts.substitution_count,
                wobble_conversion=flex_opts.wobble_conversion,
                max_n_percentage=flex_opts.max_n_percentage,
                include_reverse=flex_opts.include_reverse,
                only_vs_first=flex_opts.only_vs_first,
                line_width=flex_opts.line_width,
                length_scaling=flex_opts.length_scaling,
                gff_source_fastas=[str(path) for path in seq_files],
                gff_output_stem=flexidot_output_stem(seq_files, genome_label),
                cancel=cancel,
            )
        finally:
            if combined.is_file():
                combined.unlink(missing_ok=True)
        return finalize_flexidot_outputs(
            out_dir,
            seq_files,
            before,
            ext=filetype,
            use_gff=flex_opts.use_gff,
            gff_source_fastas=seq_files,
        )


def run_on_sequence(
    func_name: str,
    seq_path: Path,
    project: RepeatViewerProject,
    *,
    genome_path: Path | None = None,
    genome_label: str = "",
    genome_db: str | None = None,
    chrom_sizes: dict[str, int] | None = None,
    flexidot_options: FlexidotOptions | None = None,
    flexidot_mode: str = "0",
    extend: ExtendOptions | None = None,
    blast_options: BlastOptions | None = None,
    hmmer_options: HmmerOptions | None = None,
    teaid_options: TeaAidOptions | None = None,
    msa_clean_options: MsaCleanOptions | None = None,
    consensus_options: ConsensusOptions | None = None,
    te_annotation_options: TeAnnotationOptions | None = None,
    tsd_search_options: TsdSearchOptions | None = None,
    orf_prediction_settings=None,
    cancel: JobCancellation | None = None,
    sequence_progress: SequenceProgressCallback | None = None,
) -> list[str]:
    """Run one tool on one FASTA. Returns list of created output paths (strings)."""
    seq_path = seq_path.expanduser().resolve()
    out_dir = seq_path.parent
    created: list[str] = []
    # Windows has no os.uname(); the old fallback defaulted to "Darwin" and
    # launched blastn_mac → WinError 193. Always use platform.system().
    os_type = platform.system()
    active_genome = genome_path or project.genome_file
    label_suffix = genome_output_suffix(genome_label)

    if func_name == "blast":
        blast_opts = blast_options or BlastOptions()
        if blast_opts.program == BLAST_PROGRAM_ONLINE:
            from repeatviewer.backends.ncbi_blast_online import (
                NcbiBlastOnlineError,
                run_ncbi_online_blast,
                strip_online_blast_from_annotations,
            )

            try:
                outputs = run_ncbi_online_blast(
                    seq_path,
                    out_dir,
                    program=blast_opts.online_program,
                    database=blast_opts.online_database,
                    evalue=blast_opts.online_e_value,
                    hitlist_size=blast_opts.online_hitlist_size,
                    write_text=bool(blast_opts.online_write_text),
                    cancel=cancel,
                    progress=sequence_progress,
                )
            except NcbiBlastOnlineError as exc:
                raise RuntimeError(str(exc)) from exc
            for path in outputs:
                created.append(str(Path(path).resolve()))
            try:
                for ann in strip_online_blast_from_annotations(seq_path):
                    created.append(str(ann.resolve()))
            except Exception as exc:  # noqa: BLE001
                raise RuntimeError(
                    f"Online BLAST finished but clearing Annotation overlay failed: {exc}"
                ) from exc
            return created
        if genome_db is None:
            raise RuntimeError("BLAST database not prepared.")
        program = blast_opts.program
        # blastp and blastx both search a protein database and give hits in
        # protein coordinates (no strand / BED), differing only in the query.
        is_protein = program in ("blastp", "blastx")
        if program == "blastx":
            p_evalue = blast_opts.blastx_e_value
            p_word = blast_opts.blastx_word_size
            p_gap_open = blast_opts.blastx_gap_open
            p_gap_extend = blast_opts.blastx_gap_extend
            p_qcov = blast_opts.blastx_query_coverage_hsp
            p_threads = blast_opts.blastx_num_threads
            p_write_text = blast_opts.blastx_write_text
            p_write_gff = blast_opts.blastx_write_gff
        else:  # blastp
            p_evalue = blast_opts.blastp_e_value
            p_word = blast_opts.blastp_word_size
            p_gap_open = blast_opts.blastp_gap_open
            p_gap_extend = blast_opts.blastp_gap_extend
            p_qcov = blast_opts.blastp_query_coverage_hsp
            p_threads = blast_opts.blastp_num_threads
            p_write_text = blast_opts.blastp_write_text
            p_write_gff = blast_opts.blastp_write_gff

        blast_bed, blast_txt = run_blast(
            str(seq_path),
            str(genome_db),
            str(out_dir),
            e_value=(p_evalue if is_protein else blast_opts.e_value),
            word_size=(p_word if is_protein else blast_opts.word_size),
            gap_open=(p_gap_open if is_protein else blast_opts.gap_open),
            gap_extend=(p_gap_extend if is_protein else blast_opts.gap_extend),
            penalty=blast_opts.penalty,
            percent_identity=(0.0 if is_protein else blast_opts.percent_identity),
            query_coverage_hsp=(p_qcov if is_protein else blast_opts.query_coverage_hsp),
            program=program,
            bed_file=True,
            os_type=os_type,
            num_threads=(p_threads if is_protein else blast_opts.num_threads),
            use_system_blast=project.use_system_blast,
            output_label=genome_label or None,
            cancel=cancel,
            progress=sequence_progress,
        )
        if is_protein:
            if not blast_txt:
                raise RuntimeError(f"{program.upper()} failed: {blast_bed}")
            blast_txt_path = Path(str(blast_txt)).resolve()
            if p_write_gff:
                gff_path = write_blast_gff(blast_txt_path, protein_coordinates=True)
                created.append(str(gff_path.resolve()))
            else:
                blast_gff_output_path(blast_txt_path).unlink(missing_ok=True)
            if p_write_text:
                created.append(str(blast_txt_path))
            else:
                blast_txt_path.unlink(missing_ok=True)
            return created
        if not blast_bed:
            raise RuntimeError(f"Blast failed: {blast_bed}")
        max_lines = blast_opts.max_msa_lines
        top_lines = blast_opts.top_msa_lines
        blast_bed_filtered = process_bed_lines(
            str(blast_bed),
            str(out_dir),
            max_lines=max_lines,
            top_longest_lines_count=top_lines,
        )
        blast_fa = f"{blast_txt}.fa"
        extract_fasta_from_bed(
            str(active_genome),
            str(blast_bed_filtered),
            blast_fa,
        )
        # Expose selected Blast exports so the UI can highlight them in green.
        if blast_txt:
            blast_txt_path = Path(str(blast_txt)).resolve()
            if blast_opts.write_gff:
                gff_path = write_blast_gff(
                    blast_txt_path,
                    minimum_bit_score=blast_opts.minimum_bit_score,
                )
                created.append(str(gff_path.resolve()))
            else:
                blast_gff_output_path(blast_txt_path).unlink(missing_ok=True)
            # Density / genome-map jobs force write_text=True when they need
            # the hit table; honor the Export checkbox for ordinary Blast.
            if blast_opts.write_text:
                created.append(str(blast_txt_path))
            else:
                blast_txt_path.unlink(missing_ok=True)
        # Remove temporary BED intermediates.
        for tmp in (Path(str(blast_bed)), Path(str(blast_bed_filtered))):
            try:
                tmp.unlink(missing_ok=True)
            except OSError:
                pass
        created.append(blast_fa)
        if active_genome is not None and Path(active_genome).is_file():
            ref_link = write_genome_reference_link(
                out_dir,
                Path(active_genome),
                genome_label=genome_label,
            )
            created.append(str(ref_link))
        cleanup_blast_scratch_files(out_dir, seq_path, genome_label=genome_label)

    elif func_name == "extend":
        if chrom_sizes is None:
            raise RuntimeError("Genome chromosome sizes not prepared.")
        if active_genome is None:
            raise RuntimeError("Genome FASTA is required for Extend.")
        extend_opts = extend or ExtendOptions()
        left = extend_opts.left
        right = extend_opts.right
        from repeatviewer.tool_output_paths import kept_input_name

        base = kept_input_name(seq_path)
        tag = f"{label_suffix}_L{left}_R{right}"
        in_bed = out_dir / f"{base}{tag}_hdr.bed"
        ext_bed = out_dir / f"{base}{tag}.bed"
        ext_fa = out_dir / f"{base}{tag}.fa"
        fasta_header_to_bed(str(seq_path), str(in_bed))
        extend_bed_regions(str(in_bed), left, right, chrom_sizes, str(ext_bed))
        extract_fasta_from_bed(str(active_genome), str(ext_bed), str(ext_fa))
        for bed_path in (in_bed, ext_bed):
            try:
                bed_path.unlink(missing_ok=True)
            except OSError:
                pass
        created.append(str(ext_fa))

    elif func_name == "teaid":
        from repeatviewer.backends.output_paths import next_numbered_path

        teaid_opts = teaid_options or TeaAidOptions()
        html_path = next_numbered_path(
            out_dir / f"{seq_path.name}{label_suffix}_TEAid.html"
        )
        teaid_work_dir = out_dir
        scratch_dir: Path | None = None
        pan_genome_run = bool(
            teaid_opts.use_pan_genome and len(teaid_opts.pan_genomes or []) >= 2
        )
        if pan_genome_run:
            scratch_dir = out_dir / f"_teaid_{seq_path.stem}_pan"
            scratch_dir.mkdir(parents=True, exist_ok=True)
            teaid_work_dir = scratch_dir
            html_path = next_numbered_path(out_dir / f"{seq_path.name}_TEAid.html")
        elif genome_label:
            scratch_dir = out_dir / f"_teaid_{seq_path.stem}{label_suffix}"
            scratch_dir.mkdir(parents=True, exist_ok=True)
            teaid_work_dir = scratch_dir
        _patch_plotly_no_auto_open()
        try:
            if pan_genome_run:
                genome_targets = [
                    (str(name).strip(), str(Path(path).expanduser().resolve()))
                    for name, path in teaid_opts.pan_genomes
                    if str(name).strip() and str(path).strip()
                ]
                run_ok = teaid_pan_genome_plotter(
                    str(seq_path),
                    str(teaid_work_dir),
                    genome_targets,
                    e_value=teaid_opts.e_value,
                    num_threads=max(1, int(teaid_opts.num_threads)),
                    use_system_blast=project.use_system_blast,
                    self_blast_word_size=teaid_opts.self_blast_word_size,
                    x_axis_label=teaid_opts.x_axis_label,
                    hit_len_high_pct=teaid_opts.hit_len_high_pct,
                    hit_len_mid_pct=teaid_opts.hit_len_mid_pct,
                    hit_color_high=teaid_opts.hit_color_high,
                    hit_color_mid=teaid_opts.hit_color_mid,
                    hit_line_transparency=teaid_opts.hit_line_transparency,
                    progress=sequence_progress,
                )
            else:
                run_ok = teaid_plotter(
                    str(seq_path),
                    str(teaid_work_dir),
                    str(active_genome) if active_genome else None,
                    None,
                    e_value=teaid_opts.e_value,
                    num_threads=max(1, int(teaid_opts.num_threads)),
                    use_system_blast=project.use_system_blast,
                    self_blast_word_size=teaid_opts.self_blast_word_size,
                    output_name_tag=genome_label or None,
                    x_axis_label=teaid_opts.x_axis_label,
                    hit_len_high_pct=teaid_opts.hit_len_high_pct,
                    hit_len_mid_pct=teaid_opts.hit_len_mid_pct,
                    hit_color_high=teaid_opts.hit_color_high,
                    hit_color_mid=teaid_opts.hit_color_mid,
                    hit_line_transparency=teaid_opts.hit_line_transparency,
                    cancel=cancel,
                    progress=sequence_progress,
                )
        except TeaidInputError as exc:
            if scratch_dir is not None:
                shutil.rmtree(scratch_dir, ignore_errors=True)
            raise RuntimeError(str(exc)) from exc
        if not run_ok:
            if scratch_dir is not None:
                shutil.rmtree(scratch_dir, ignore_errors=True)
            raise RuntimeError("TE-Aid returned failure.")
        scratch_name = f"{seq_path.name}_TEAid.html" if pan_genome_run else (
            f"{seq_path.name}{label_suffix}_TEAid.html"
        )
        scratch_html = teaid_work_dir / scratch_name
        if scratch_dir is not None:
            if not scratch_html.is_file():
                shutil.rmtree(scratch_dir, ignore_errors=True)
                raise RuntimeError(f"TE-Aid HTML not created: {scratch_html}")
            if html_path.exists():
                html_path.unlink()
            shutil.move(str(scratch_html), str(html_path))
            cleanup_teaid_side_files(scratch_dir, seq_path, scratch_html)
            shutil.rmtree(scratch_dir, ignore_errors=True)
        else:
            if not html_path.is_file():
                raise RuntimeError(f"TE-Aid HTML not created: {html_path}")
            cleanup_teaid_side_files(out_dir, seq_path, html_path)
        created.append(str(html_path.resolve()))

    elif func_name == "msa_clean_gappy":
        from repeatviewer.backends.output_paths import next_numbered_path

        msa_opts = msa_clean_options or MsaCleanOptions()
        g = msa_opts.gappy
        gcol_out = str(next_numbered_path(Path(f"{seq_path}_gap.fa")))
        remove_gaps_with_similarity_check(
            str(seq_path),
            gcol_out,
            gap_threshold=g.gap_threshold,
            simi_check_gap_thr=g.simi_check_gap_thr,
            similarity_thr=g.similarity_thr,
            min_nucleotide=g.min_nucleotide,
        )
        created.append(gcol_out)

    elif func_name == "msa_clean":
        from repeatviewer.backends.output_paths import next_numbered_path

        msa_opts = msa_clean_options or MsaCleanOptions()
        div = msa_opts.divergence
        gap = msa_opts.gap_row
        g = msa_opts.gappy
        enabled_methods = [
            method
            for method in msa_opts.method_order
            if (
                (method == "msa_clean_div" and msa_opts.use_divergence)
                or (method == "msa_clean_gap" and msa_opts.use_gap_row)
                or (method == "msa_clean_gappy" and msa_opts.use_gappy_columns)
            )
        ]
        if not enabled_methods:
            raise ValueError("Enable at least one MSA cleaning method.")

        combined_suffix = {
            "msa_clean_gappy": "cgap",
            "msa_clean_div": "crowdiv",
            "msa_clean_gap": "crowgap",
        }
        final_output = str(
            next_numbered_path(
                Path(f"{seq_path}_{combined_suffix[enabled_methods[-1]]}.fa")
            )
        )
        with tempfile.TemporaryDirectory(prefix="repeatviewer_msa_clean_") as temp_dir:
            current = str(seq_path)
            for index, method in enumerate(enabled_methods):
                is_last = index == len(enabled_methods) - 1
                output = (
                    final_output
                    if is_last
                    else str(Path(temp_dir) / f"stage_{index + 1}.fa")
                )
                if method == "msa_clean_div":
                    crop_end_div(
                        current,
                        output,
                        threshold=div.threshold,
                        window_size=div.window_size,
                    )
                elif method == "msa_clean_gap":
                    crop_end_gap(
                        current,
                        output,
                        gap_threshold=gap.gap_threshold,
                        window_size=gap.window_size,
                    )
                else:
                    remove_gaps_with_similarity_check(
                        current,
                        output,
                        gap_threshold=g.gap_threshold,
                        simi_check_gap_thr=g.simi_check_gap_thr,
                        similarity_thr=g.similarity_thr,
                        min_nucleotide=g.min_nucleotide,
                    )
                current = output
        created.append(final_output)

    elif func_name in {"msa_clean_div", "msa_clean_gap"}:
        from repeatviewer.backends.output_paths import next_numbered_path

        msa_opts = msa_clean_options or MsaCleanOptions()
        div = msa_opts.divergence
        gap = msa_opts.gap_row
        div_out = str(next_numbered_path(Path(f"{seq_path}_rowdiv.fa")))
        gap_out = str(next_numbered_path(Path(f"{seq_path}_rowgap.fa")))
        if func_name == "msa_clean_div":
            crop_end_div(
                str(seq_path),
                div_out,
                threshold=div.threshold,
                window_size=div.window_size,
            )
            created.append(div_out)
        else:
            crop_end_gap(
                str(seq_path),
                gap_out,
                gap_threshold=gap.gap_threshold,
                window_size=gap.window_size,
            )
            created.append(gap_out)

    elif func_name == "consensus_generation":
        cons_opts = consensus_options or ConsensusOptions()
        cons_out, _ = con_generater(
            str(seq_path),
            str(out_dir),
            threshold=cons_opts.threshold,
            ambiguous=cons_opts.ambiguous,
            sequence_type=cons_opts.sequence_type,
        )
        created.append(cons_out)

    elif func_name == "flexidot":
        from repeatviewer.backends.job_phases import FLEXIDOT_PHASES, emit_phase

        flex_opts = flexidot_options or FlexidotOptions()
        seq_count = sum(1 for _ in SeqIO.parse(str(seq_path), "fasta"))
        if seq_count < 1:
            raise RuntimeError(f"No sequences in {seq_path.name} (Flexidot needs at least one).")
        out_path = Path(out_dir).resolve()
        flexidot_filetype = effective_flexidot_filetype(flex_opts.output_format)
        if flex_opts.modes:
            mode: str | list[str] = list(flex_opts.modes)
        else:
            mode = auto_flexidot_mode(seq_count)
        emit_phase(
            sequence_progress, 0, FLEXIDOT_PHASES, f"{seq_path.name}: preparing Flexidot…"
        )
        with _flexidot_dir_lock(out_path):
            before = _output_paths_in_dir(out_path, flexidot_filetype)
            emit_phase(
                sequence_progress,
                1,
                FLEXIDOT_PHASES,
                f"{seq_path.name}: running Flexidot…",
            )
            run_flexidot(
                str(seq_path),
                str(out_path),
                wordsize=flex_opts.wordsize,
                mode=mode,
                filetype=flexidot_filetype,
                sequence_type=flex_opts.sequence_type,
                lcs_shading=flex_opts.lcs_shading,
                use_gff=flex_opts.use_gff,
                substitution_count=flex_opts.substitution_count,
                wobble_conversion=flex_opts.wobble_conversion,
                max_n_percentage=flex_opts.max_n_percentage,
                include_reverse=flex_opts.include_reverse,
                only_vs_first=flex_opts.only_vs_first,
                line_width=flex_opts.line_width,
                length_scaling=flex_opts.length_scaling,
                gff_source_fastas=[str(seq_path)],
                gff_output_stem=flexidot_output_stem([seq_path.resolve()], genome_label),
                cancel=cancel,
            )
            emit_phase(
                sequence_progress,
                2,
                FLEXIDOT_PHASES,
                f"{seq_path.name}: writing Flexidot outputs…",
            )
            created.extend(
                finalize_flexidot_outputs(
                    out_path,
                    [seq_path.resolve()],
                    before,
                    ext=flexidot_filetype,
                    use_gff=flex_opts.use_gff,
                    gff_source_fastas=[seq_path.resolve()],
                )
            )
        emit_phase(
            sequence_progress, FLEXIDOT_PHASES, FLEXIDOT_PHASES, f"{seq_path.name}: finished"
        )

    elif func_name == "te_annotation":
        te_opts = te_annotation_options or TeAnnotationOptions()
        outs = run_te_annotation_on_fasta(
            seq_path,
            te_opts,
            cancel=cancel,
            progress=sequence_progress,
        )
        created.extend(str(p) for p in outs)

    elif func_name == "reverse_complement":
        created.append(str(write_reverse_complement_fasta(seq_path)))

    elif func_name == "orf_prediction":
        created.append(
            str(write_predicted_orfs_fasta(seq_path, settings=orf_prediction_settings))
        )

    elif func_name == "dna_translation":
        created.append(
            str(write_translated_orfs_fasta(seq_path, settings=orf_prediction_settings))
        )

    elif func_name == "tsd_search":
        tsd_opts = tsd_search_options or TsdSearchOptions()
        if tsd_opts.left_boundary <= 0 or tsd_opts.right_boundary <= 0:
            raise RuntimeError(
                "Set left and right MSA boundary columns before running TSD search."
            )
        if tsd_opts.left_boundary >= tsd_opts.right_boundary:
            raise RuntimeError("Left boundary must be less than the right boundary.")
        from repeatviewer.tool_output_paths import kept_input_name

        out_file = out_dir / f"{kept_input_name(seq_path)}.tsd_search.txt"
        TSD_identifier(
            tsd_opts.left_boundary,
            tsd_opts.right_boundary,
            tsd_opts.nt_radius,
            seq_path,
            out_file,
            max_sim_threshold=tsd_opts.max_sim_threshold,
            nt_threshold=tsd_opts.nt_threshold,
            cluster_eps=tsd_opts.cluster_eps,
        )
        created.append(str(out_file.resolve()))

    elif func_name == "hmmer":
        from repeatviewer.backends.hmmer_tool import resolve_hmm_database, run_hmmer_search
        from repeatviewer.database_settings import load_user_databases_from_settings

        hmmer_opts = hmmer_options or HmmerOptions()
        db = resolve_hmm_database(
            load_user_databases_from_settings(),
            hmmer_opts.database_name,
        )
        outs = run_hmmer_search(seq_path, db, options=hmmer_opts, cancel=cancel)
        created.extend(str(p) for p in outs)

    else:
        raise ValueError(f"Unknown function: {func_name}")

    return created


def run_batch(
    func_name: str,
    seq_files: Iterable[Path],
    project: RepeatViewerProject,
    progress_cb: ProgressCallback,
    *,
    flexidot_options: FlexidotOptions | None = None,
    flexidot_mode: str = "0",
    extend: ExtendOptions | None = None,
    blast_options: BlastOptions | None = None,
    hmmer_options: HmmerOptions | None = None,
    teaid_options: TeaAidOptions | None = None,
    msa_clean_options: MsaCleanOptions | None = None,
    consensus_options: ConsensusOptions | None = None,
    te_annotation_options: TeAnnotationOptions | None = None,
    tsd_search_options: TsdSearchOptions | None = None,
    orf_prediction_settings=None,
    max_workers: int = 4,
    flexidot_combined: bool = True,
    cancel: JobCancellation | None = None,
    on_outputs_ready: Callable[[list[str]], None] | None = None,
) -> dict:
    seq_list = list(seq_files)
    failures: list[str] = []
    created_outputs: list[str] = []
    source_files: list[str] = []
    os_type = platform.system()

    genome_dependent = func_name in {"blast", "extend", "teaid", "flexidot"}
    genome_items: list[tuple[str, Path]] = []
    db_cache: dict[Path, str] = {}
    chrom_cache: dict[Path, dict[str, int]] = {}

    selected_blast_db_name = blast_options.active_database_name if blast_options else None
    if func_name == "blast" and blast_options and blast_options.program == BLAST_PROGRAM_ONLINE:
        # NCBI web BLAST — no local genome / database preparation.
        genome_dependent = False
        genome_items = []
        db_cache = {}
        chrom_cache = {}
    elif func_name == "blast" and selected_blast_db_name:
        expected_molecule = (
            "prot"
            if blast_options and blast_options.program in ("blastp", "blastx")
            else "nucl"
        )
        db = _find_blast_database(
            project,
            selected_blast_db_name,
            molecule_type=expected_molecule,
        )
        genome_items = [(db.name, db.source_fasta)]
        db_cache = {db.source_fasta: str(db.db_prefix_path)}
    elif func_name == "blast" and blast_options and blast_options.program == "blastp":
        raise RuntimeError("Select a protein database in the blastp parameter tab.")
    elif genome_dependent:
        teaid_pan = (
            func_name == "teaid"
            and teaid_options is not None
            and bool(teaid_options.use_pan_genome)
        )
        if not teaid_pan and not project.has_genome_configured():
            raise RuntimeError("Genome FASTA is required for this tool.")
        genome_items = _genome_work_items(
            project, func_name, extend, teaid=teaid_options
        )
        db_cache, chrom_cache = _prepare_genome_resources(
            func_name, genome_items, project, os_type, cancel=cancel
        )

    display_name = _DISPLAY_NAMES.get(func_name, func_name)

    if func_name == "flexidot" and len(seq_list) > 1 and flexidot_combined:
        from repeatviewer.backends.job_phases import FLEXIDOT_PHASES

        targets = genome_items or [("", Path(project.genome_file or ""))]
        total_steps = max(1, len(targets) * FLEXIDOT_PHASES)
        progress_cb(0, f"{display_name}: {len(seq_list)} sequences (comparison)")
        check_cancelled(cancel)
        completed = 0
        for genome_label, genome_path in targets:
            check_cancelled(cancel)
            label = f" ({genome_label})" if genome_label else ""
            try:
                progress_cb(
                    completed,
                    f"{display_name}: preparing comparison{label}",
                )
                progress_cb(
                    completed + 1,
                    f"{display_name}: running Flexidot{label}",
                )
                outs = run_flexidot_combined(
                    seq_list,
                    flexidot_options=flexidot_options,
                    cancel=cancel,
                    genome_label=genome_label,
                )
                progress_cb(
                    completed + 2,
                    f"{display_name}: writing outputs{label}",
                )
                created_outputs.extend(outs)
                for seq_path in seq_list:
                    source_files.append(str(seq_path.resolve()))
            except JobCancelledError:
                raise
            except Exception as exc:  # noqa: BLE001
                failures.append(f"batch{label}: {exc}")
            completed += FLEXIDOT_PHASES
            progress_cb(completed, f"{display_name}: finished{label}")
        return {
            "func_name": func_name,
            "total": total_steps,
            "failures": failures,
            "created_outputs": created_outputs,
            "source_files": source_files,
        }

    work_items: list[tuple[Path, str, Path | None]] = []
    if genome_dependent and genome_items:
        teaid_pan = (
            func_name == "teaid"
            and teaid_options is not None
            and teaid_options.use_pan_genome
        )
        if teaid_pan:
            for seq_path in seq_list:
                work_items.append((seq_path, "", None))
        else:
            for seq_path in seq_list:
                for genome_label, genome_path in genome_items:
                    work_items.append((seq_path, genome_label, genome_path))
    else:
        for seq_path in seq_list:
            work_items.append((seq_path, "", None))

    per_sequence = func_name in _PER_SEQUENCE_FUNCS
    # Online BLAST needs phase updates (submit / poll / download) even though
    # each FASTA file is one work unit.
    online_blast = (
        func_name == "blast"
        and blast_options is not None
        and blast_options.program == BLAST_PROGRAM_ONLINE
    )
    from repeatviewer.backends.job_phases import tool_units_for_fasta

    phased_tools = {"blast", "flexidot", "teaid"}
    pass_sequence_progress = per_sequence or online_blast or func_name in phased_tools
    item_units = [
        (
            _fasta_sequence_units(seq_path)
            if per_sequence
            else tool_units_for_fasta(
                func_name, seq_path, online_blast=online_blast
            )
        )
        for seq_path, _genome_label, _genome_path in work_items
    ]
    total_units = max(1, sum(item_units))

    workers = 1 if func_name == "te_annotation" else effective_batch_workers(
        max_workers, len(work_items)
    )
    # NCBI asks clients not to hammer the service with parallel web BLAST jobs.
    if online_blast:
        workers = 1
    progress_lock = threading.Lock()
    completed_count = 0
    from repeatviewer.external_process import current_process_log_sink, process_log_sink

    log_sink = current_process_log_sink()

    def _run_single(
        seq_path: Path,
        genome_label: str,
        genome_path: Path | None,
        units: int,
    ) -> tuple[Path, str, list[str], str | None, int, int]:
        """Returns path, label, outs, error, units, sequences_reported."""
        last_reported = 0

        def _sequence_progress(completed: int, total: int, label: str) -> None:
            nonlocal last_reported, completed_count
            display = f"{display_name}: {label}"
            with progress_lock:
                if completed < 0:
                    progress_cb(completed_count, display)
                    return
                done = max(0, min(int(completed), int(units)))
                delta = done - last_reported
                if delta > 0:
                    last_reported = done
                    completed_count += delta
                progress_cb(completed_count, display)

        try:
            with process_log_sink(log_sink):
                check_cancelled(cancel)
                outs = run_on_sequence(
                    func_name,
                    seq_path,
                    project,
                    genome_path=genome_path,
                    genome_label=genome_label,
                    genome_db=db_cache.get(genome_path) if genome_path else None,
                    chrom_sizes=chrom_cache.get(genome_path) if genome_path else None,
                    flexidot_options=flexidot_options,
                    flexidot_mode=(
                        ",".join(flexidot_options.modes)
                        if func_name == "flexidot"
                        and flexidot_options is not None
                        and flexidot_options.modes
                        else "0"
                        if func_name == "flexidot"
                        else flexidot_mode
                    ),
                    extend=extend,
                    blast_options=blast_options,
                    hmmer_options=hmmer_options,
                    teaid_options=teaid_options,
                    msa_clean_options=msa_clean_options,
                    consensus_options=consensus_options,
                    te_annotation_options=te_annotation_options,
                    tsd_search_options=tsd_search_options,
                    orf_prediction_settings=orf_prediction_settings,
                    cancel=cancel,
                    sequence_progress=(
                        _sequence_progress if pass_sequence_progress else None
                    ),
                )
            return seq_path, genome_label, outs, None, units, last_reported
        except JobCancelledError:
            raise
        except Exception as exc:  # noqa: BLE001
            return seq_path, genome_label, [], str(exc), units, last_reported

    def _record_finished(
        seq_path: Path,
        genome_label: str,
        outs: list[str],
        error: str | None,
        units: int,
        reported: int,
    ) -> None:
        nonlocal completed_count
        label = f" ({genome_label})" if genome_label else ""
        with progress_lock:
            # Phase / per-sequence progress already advanced completed_count;
            # only add any remaining units for this work item.
            if pass_sequence_progress:
                remaining = max(0, int(units) - int(reported))
                if remaining:
                    completed_count += remaining
            else:
                completed_count += 1
            progress_cb(
                completed_count,
                f"{display_name}: {seq_path.name}{label}",
            )
        if error:
            failures.append(f"{seq_path.name}{label}: {error}")
        else:
            created_outputs.extend(outs)
            source_files.append(str(seq_path.resolve()))
            if on_outputs_ready is not None and outs:
                try:
                    on_outputs_ready(list(outs))
                except Exception:  # noqa: BLE001 — UI refresh must not abort the batch
                    pass

    def _cancelled_result() -> dict:
        return {
            "func_name": func_name,
            "total": total_units,
            "failures": failures,
            "created_outputs": created_outputs,
            "source_files": source_files,
            "cancelled": True,
        }

    try:
        if workers == 1:
            for (seq_path, genome_label, genome_path), units in zip(
                work_items, item_units, strict=True
            ):
                check_cancelled(cancel)
                label = f" ({genome_label})" if genome_label else ""
                progress_cb(
                    completed_count,
                    f"{display_name}: {seq_path.name}{label}",
                )
                path, glabel, outs, error, units_out, reported = _run_single(
                    seq_path, genome_label, genome_path, units
                )
                _record_finished(path, glabel, outs, error, units_out, reported)
        else:
            with ThreadPoolExecutor(max_workers=workers) as pool:
                futures = [
                    pool.submit(_run_single, seq_path, genome_label, genome_path, units)
                    for (seq_path, genome_label, genome_path), units in zip(
                        work_items, item_units, strict=True
                    )
                ]
                for future in as_completed(futures):
                    try:
                        path, label, outs, error, units_out, reported = future.result()
                    except JobCancelledError:
                        for pending in futures:
                            pending.cancel()
                        # Keep outputs from work items that already finished.
                        raise
                    _record_finished(path, label, outs, error, units_out, reported)
    except JobCancelledError:
        # Return whatever finished so the UI can show partial Online BLAST (etc.).
        return _cancelled_result()

    progress_cb(total_units, f"{display_name}: finished")

    return {
        "func_name": func_name,
        "total": total_units,
        "failures": failures,
        "created_outputs": created_outputs,
        "source_files": source_files,
    }


def _patch_plotly_no_auto_open() -> None:
    _orig = pio.write_html

    def _write_html(*args, **kwargs):
        kwargs.setdefault("auto_open", False)
        return _orig(*args, **kwargs)

    pio.write_html = _write_html
