"""TE curation tools (ported from TEtrimmerGUI, standalone in TEviewer)."""

from __future__ import annotations

import os
import shutil
import threading
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Callable, Iterable

import plotly.io as pio
from Bio import SeqIO

from teviewer.backends.crop_end_divergence import crop_end_div
from teviewer.backends.crop_end_gap import crop_end_gap
from teviewer.backends.gui_functions import (
    sanitize_name,
    blast as run_blast,
    check_database,
    extend_bed_regions,
    extract_fasta_from_bed,
    fasta_header_to_bed,
    process_bed_lines,
    run_flexidot,
)
from teviewer.backends.remove_gap import remove_gaps_with_similarity_check
from teviewer.backends.teaid_plotter import con_generater, teaid_plotter
from teviewer.backends.job_cancel import JobCancelledError, JobCancellation, check_cancelled
from teviewer.backends.thread_budget import analysis_threads_per_job
from teviewer.job_options import (
    BlastOptions,
    ConsensusOptions,
    ExtendOptions,
    FlexidotOptions,
    MsaCleanOptions,
    TeaAidOptions,
)
from teviewer.backends.te_annotation_batch import TeAnnotationOptions, run_te_annotation_on_fasta
from teviewer.project import TEviewerProject

ProgressCallback = Callable[[int, str], None]

_DISPLAY_NAMES = {
    "teaid": "TEAid",
    "blast": "Blast",
    "extend": "Extend",
    "msa_clean": "MSA clean",
    "msa_clean_gappy": "MSA remove gappy columns",
    "msa_clean_div": "Clean row with divergence",
    "msa_clean_gap": "Clean row with gap",
    "msa_clean_row_col": "Clean row and column",
    "consensus_generation": "Consensus",
    "flexidot": "Flexidot",
    "te_annotation": "TE annotation",
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
    Remove intermediate TEAid files for one input only.

    Keeps the input FASTA and this run's HTML report. Safe when several TEAid
    jobs run in parallel in the same folder (does not delete other sequences'
    HTML reports or their source files).
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
        if name == f"{input_name}_self_blast" and entry.is_dir():
            to_remove.append(resolved)
            continue
        if not (entry.is_file() or entry.is_symlink()):
            continue
        if name.startswith(consensus_prefix):
            to_remove.append(resolved)
            continue
        if name.startswith(input_name) and name != input_name:
            to_remove.append(resolved)

    for path in sorted(to_remove, key=lambda p: len(p.parts), reverse=True):
        try:
            if path.is_dir():
                shutil.rmtree(path, ignore_errors=True)
            else:
                path.unlink(missing_ok=True)
        except OSError:
            pass


def _png_paths_in_dir(folder: Path) -> set[Path]:
    folder = folder.expanduser().resolve()
    return {p.resolve() for p in folder.glob("*.png") if p.is_file()}


def flexidot_output_stem(seq_files: list[Path]) -> str:
    """Output basename: ``{first_stem}_flexidot`` or ``{first_stem}_flexidot_multi``."""
    stem = seq_files[0].expanduser().resolve().stem
    if len(seq_files) >= 2:
        return f"{stem}_flexidot_multi"
    return f"{stem}_flexidot"


def _flexidot_pngs_for_run(
    output_dir: Path,
    before: set[Path],
    seq_files: list[Path],
) -> list[Path]:
    """New PNGs from one Flexidot run (optionally filtered to one input stem)."""
    after = _png_paths_in_dir(output_dir)
    new_files = sorted(
        after - before,
        key=lambda p: p.stat().st_size,
        reverse=True,
    )
    if len(seq_files) == 1:
        stem = seq_files[0].stem
        owned = [p for p in new_files if stem in p.name]
        if owned:
            return owned
    return new_files


def cleanup_flexidot_intermediates(output_dir: Path, keep_paths: Iterable[Path]) -> None:
    """Remove Flexidot scratch files (``flexidot_output-*``) after plots are renamed."""
    output_dir = output_dir.expanduser().resolve()
    keep = {p.expanduser().resolve() for p in keep_paths}
    if not output_dir.is_dir():
        return
    for entry in output_dir.iterdir():
        if not entry.name.startswith("flexidot_output"):
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


def finalize_flexidot_outputs(
    output_dir: Path,
    seq_files: list[Path],
    before: set[Path],
) -> list[str]:
    """Rename new Flexidot PNGs to ``{stem}_flexidot[.png|_N.png]``."""
    output_dir = output_dir.expanduser().resolve()
    new_files = _flexidot_pngs_for_run(output_dir, before, seq_files)
    if not new_files:
        return []

    base = flexidot_output_stem(seq_files)
    created: list[str] = []
    dest_paths: list[Path] = []
    for index, src in enumerate(new_files):
        dest = output_dir / (f"{base}.png" if index == 0 else f"{base}_{index + 1}.png")
        if dest.exists():
            dest.unlink()
        if not src.is_file():
            continue
        src.rename(dest)
        dest_paths.append(dest)
        created.append(str(dest.resolve()))
    cleanup_flexidot_intermediates(output_dir, dest_paths)
    return created


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
) -> list[str]:
    """Flexidot for 2+ inputs: mode 1 (paired) or 2 (poly)."""
    if len(seq_files) < 2:
        raise ValueError("Combined Flexidot requires at least two sequences.")
    flex_opts = flexidot_options or FlexidotOptions()
    mode = "1" if len(seq_files) == 2 else "2"
    out_dir = seq_files[0].parent.resolve()
    combined = out_dir / f"_flexidot_{len(seq_files)}seq.fa"
    records = []
    for path in seq_files:
        records.extend(SeqIO.parse(str(path.expanduser().resolve()), "fasta"))
    if not records:
        raise RuntimeError("No sequences found in selected FASTA files.")
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
        before = _png_paths_in_dir(out_dir)
        try:
            run_flexidot(
                str(combined),
                str(out_dir),
                wordsize=flex_opts.wordsize,
                mode=mode,
                lcs_shading=flex_opts.lcs_shading,
                use_gff=flex_opts.use_gff,
                gff_source_fastas=[str(path) for path in seq_files],
                gff_output_stem=flexidot_output_stem(seq_files),
                cancel=cancel,
            )
        finally:
            if combined.is_file():
                combined.unlink(missing_ok=True)
        return finalize_flexidot_outputs(out_dir, seq_files, before)


def run_on_sequence(
    func_name: str,
    seq_path: Path,
    project: TEviewerProject,
    *,
    genome_db: str | None = None,
    chrom_sizes: dict[str, int] | None = None,
    flexidot_options: FlexidotOptions | None = None,
    flexidot_mode: str = "0",
    extend: ExtendOptions | None = None,
    blast_options: BlastOptions | None = None,
    teaid_options: TeaAidOptions | None = None,
    msa_clean_options: MsaCleanOptions | None = None,
    consensus_options: ConsensusOptions | None = None,
    te_annotation_options: TeAnnotationOptions | None = None,
    cancel: JobCancellation | None = None,
) -> list[str]:
    """Run one tool on one FASTA. Returns list of created output paths (strings)."""
    seq_path = seq_path.expanduser().resolve()
    out_dir = seq_path.parent
    created: list[str] = []
    os_type = os.uname().sysname if hasattr(os, "uname") else "Darwin"

    if func_name == "blast":
        if genome_db is None:
            raise RuntimeError("Genome BLAST database not prepared.")
        blast_opts = blast_options or BlastOptions()
        blast_bed, blast_txt = run_blast(
            str(seq_path),
            str(genome_db),
            str(out_dir),
            e_value=blast_opts.e_value,
            bed_file=True,
            os_type=os_type,
            num_threads=analysis_threads_per_job(project.batch_workers),
            use_system_blast=project.use_system_blast,
        )
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
            str(project.genome_file),
            str(blast_bed_filtered),
            blast_fa,
        )
        # Expose both Blast outputs so the UI highlights all newly generated files.
        if blast_txt:
            created.append(str(Path(str(blast_txt)).resolve()))
        # Keep *_blast.txt for downstream plots; remove temporary BED intermediates.
        for tmp in (Path(str(blast_bed)), Path(str(blast_bed_filtered))):
            try:
                tmp.unlink(missing_ok=True)
            except OSError:
                pass
        created.append(blast_fa)

    elif func_name == "extend":
        if chrom_sizes is None:
            raise RuntimeError("Genome chromosome sizes not prepared.")
        extend_opts = extend or ExtendOptions()
        left = extend_opts.left
        right = extend_opts.right
        stem = seq_path.stem
        tag = f"_L{left}_R{right}"
        in_bed = out_dir / f"{stem}{tag}_hdr.bed"
        ext_bed = out_dir / f"{stem}{tag}.bed"
        ext_fa = out_dir / f"{stem}{tag}.fa"
        fasta_header_to_bed(str(seq_path), str(in_bed))
        extend_bed_regions(str(in_bed), left, right, chrom_sizes, str(ext_bed))
        extract_fasta_from_bed(str(project.genome_file), str(ext_bed), str(ext_fa))
        for bed_path in (in_bed, ext_bed):
            try:
                bed_path.unlink(missing_ok=True)
            except OSError:
                pass
        created.append(str(ext_fa))

    elif func_name == "teaid":
        teaid_opts = teaid_options or TeaAidOptions()
        out_path = out_dir
        html_path = out_path / f"{seq_path.name}_TEAid.html"
        _patch_plotly_no_auto_open()
        run_ok = teaid_plotter(
            str(seq_path),
            str(out_path),
            str(project.genome_file) if project.genome_file else None,
            None,
            e_value=teaid_opts.e_value,
            num_threads=analysis_threads_per_job(project.batch_workers),
            use_system_blast=project.use_system_blast,
            self_blast_word_size=teaid_opts.self_blast_word_size,
        )
        if not run_ok:
            raise RuntimeError("TEAid returned failure.")
        if not html_path.is_file():
            raise RuntimeError(f"TEAid HTML not created: {html_path}")
        cleanup_teaid_side_files(out_path, seq_path, html_path)
        created.append(str(html_path.resolve()))

    elif func_name == "msa_clean_gappy":
        msa_opts = msa_clean_options or MsaCleanOptions()
        g = msa_opts.gappy
        gcol_out = f"{seq_path}_gcol.fa"
        remove_gaps_with_similarity_check(
            str(seq_path),
            gcol_out,
            gap_threshold=g.gap_threshold,
            simi_check_gap_thr=g.simi_check_gap_thr,
            similarity_thr=g.similarity_thr,
            min_nucleotide=g.min_nucleotide,
        )
        created.append(gcol_out)

    elif func_name == "msa_clean_row_col":
        msa_opts = msa_clean_options or MsaCleanOptions()
        div = msa_opts.divergence
        g = msa_opts.gappy
        div_out = f"{seq_path}_cd.fa"
        clean_out = f"{seq_path}_clean.fa"
        crop_end_div(
            str(seq_path),
            div_out,
            threshold=div.threshold,
            window_size=div.window_size,
        )
        remove_gaps_with_similarity_check(
            div_out,
            clean_out,
            gap_threshold=g.gap_threshold,
            simi_check_gap_thr=g.simi_check_gap_thr,
            similarity_thr=g.similarity_thr,
            min_nucleotide=g.min_nucleotide,
        )
        created.append(clean_out)

    elif func_name in {"msa_clean_div", "msa_clean_gap"}:
        msa_opts = msa_clean_options or MsaCleanOptions()
        div = msa_opts.divergence
        gap = msa_opts.gap_row
        div_out = f"{seq_path}_cd.fa"
        gap_out = f"{seq_path}_cg.fa"
        crop_end_div(
            str(seq_path),
            div_out,
            threshold=div.threshold,
            window_size=div.window_size,
        )
        if func_name == "msa_clean_div":
            created.append(div_out)
        else:
            crop_end_gap(
                div_out,
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
        )
        created.append(cons_out)

    elif func_name == "flexidot":
        flex_opts = flexidot_options or FlexidotOptions()
        seq_count = sum(1 for _ in SeqIO.parse(str(seq_path), "fasta"))
        if seq_count < 1:
            raise RuntimeError(f"No sequences in {seq_path.name} (Flexidot needs at least one).")
        out_path = Path(out_dir).resolve()
        with _flexidot_dir_lock(out_path):
            before = _png_paths_in_dir(out_path)
            run_flexidot(
                str(seq_path),
                str(out_path),
                wordsize=flex_opts.wordsize,
                mode=flexidot_mode,
                lcs_shading=flex_opts.lcs_shading,
                use_gff=flex_opts.use_gff,
                gff_source_fastas=[str(seq_path)],
                gff_output_stem=flexidot_output_stem([seq_path.resolve()]),
                cancel=cancel,
            )
            created.extend(
                finalize_flexidot_outputs(out_path, [seq_path.resolve()], before)
            )

    elif func_name == "te_annotation":
        te_opts = te_annotation_options or TeAnnotationOptions()
        outs = run_te_annotation_on_fasta(
            seq_path,
            te_opts,
            cancel=cancel,
        )
        created.extend(str(p) for p in outs)

    else:
        raise ValueError(f"Unknown function: {func_name}")

    return created


def run_batch(
    func_name: str,
    seq_files: Iterable[Path],
    project: TEviewerProject,
    progress_cb: ProgressCallback,
    *,
    flexidot_options: FlexidotOptions | None = None,
    flexidot_mode: str = "0",
    extend: ExtendOptions | None = None,
    blast_options: BlastOptions | None = None,
    teaid_options: TeaAidOptions | None = None,
    msa_clean_options: MsaCleanOptions | None = None,
    consensus_options: ConsensusOptions | None = None,
    te_annotation_options: TeAnnotationOptions | None = None,
    max_workers: int = 4,
    flexidot_combined: bool = True,
    cancel: JobCancellation | None = None,
) -> dict:
    seq_list = list(seq_files)
    failures: list[str] = []
    created_outputs: list[str] = []
    source_files: list[str] = []
    os_type = os.uname().sysname if hasattr(os, "uname") else "Darwin"

    genome_db = None
    chrom_sizes = None

    if func_name == "blast":
        if not project.genome_file:
            raise RuntimeError("Genome FASTA is required for Blast.")
        genome_db = check_database(
            str(project.genome_file),
            os_type=os_type,
            use_system_blast=project.use_system_blast,
        )
        if genome_db in {"makeblastdb_not_found", "makeblastdb_got_error", None, False}:
            raise RuntimeError(f"Cannot prepare genome BLAST database: {genome_db}")

    if func_name == "extend":
        if not project.genome_file:
            raise RuntimeError("Genome FASTA is required for Extend.")
        chrom_sizes = genome_chrom_sizes(project.genome_file)

    display_name = _DISPLAY_NAMES.get(func_name, func_name)

    if func_name == "flexidot" and len(seq_list) > 1 and flexidot_combined:
        progress_cb(0, f"{display_name}: {len(seq_list)} sequences (comparison)")
        check_cancelled(cancel)
        try:
            outs = run_flexidot_combined(
                seq_list,
                flexidot_options=flexidot_options,
                cancel=cancel,
            )
            created_outputs.extend(outs)
            for seq_path in seq_list:
                source_files.append(str(seq_path.resolve()))
        except JobCancelledError:
            raise
        except Exception as exc:  # noqa: BLE001
            failures.append(f"batch: {exc}")
        progress_cb(len(seq_list), f"{display_name}: finished")
        return {
            "func_name": func_name,
            "total": len(seq_list),
            "failures": failures,
            "created_outputs": created_outputs,
            "source_files": source_files,
        }

    workers = 1 if func_name == "te_annotation" else effective_batch_workers(max_workers, len(seq_list))
    progress_lock = threading.Lock()
    completed_count = 0

    def _run_single(seq_path: Path) -> tuple[Path, list[str], str | None]:
        try:
            check_cancelled(cancel)
            outs = run_on_sequence(
                func_name,
                seq_path,
                project,
                genome_db=genome_db,
                chrom_sizes=chrom_sizes,
                flexidot_options=flexidot_options,
                flexidot_mode="0" if func_name == "flexidot" else flexidot_mode,
                extend=extend,
                blast_options=blast_options,
                teaid_options=teaid_options,
                msa_clean_options=msa_clean_options,
                consensus_options=consensus_options,
                te_annotation_options=te_annotation_options,
                cancel=cancel,
            )
            return seq_path, outs, None
        except JobCancelledError:
            raise
        except Exception as exc:  # noqa: BLE001
            return seq_path, [], str(exc)

    def _record_finished(seq_path: Path, outs: list[str], error: str | None) -> None:
        nonlocal completed_count
        with progress_lock:
            completed_count += 1
            progress_cb(completed_count - 1, f"{display_name}: {seq_path.name}")
        if error:
            failures.append(f"{seq_path.name}: {error}")
        else:
            created_outputs.extend(outs)
            source_files.append(str(seq_path.resolve()))

    if workers == 1:
        for seq_path in seq_list:
            check_cancelled(cancel)
            path, outs, error = _run_single(seq_path)
            _record_finished(path, outs, error)
    else:
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(_run_single, seq_path) for seq_path in seq_list]
            for future in as_completed(futures):
                try:
                    path, outs, error = future.result()
                except JobCancelledError:
                    for pending in futures:
                        pending.cancel()
                    raise
                _record_finished(path, outs, error)

    progress_cb(len(seq_list), f"{display_name}: finished")

    return {
        "func_name": func_name,
        "total": len(seq_list),
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
