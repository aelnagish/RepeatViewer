"""Build and run a RepeatMasker command; return created output paths."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from repeatviewer.backends.job_cancel import (
    JobCancellation,
    JobCancelledError,
    run_subprocess_with_cancel,
)
from repeatviewer.backends.repeat_quantify import _sanitize_fasta_ids
from repeatviewer.backends.te_library_benchmark import repeatmasker_executable


def build_repeatmasker_cmd(
    *,
    exe: str,
    query_fasta: Path,
    output_dir: Path,
    options: dict[str, Any],
) -> list[str]:
    """Assemble ``RepeatMasker [options] query -dir out`` from dialog options."""
    pa = max(1, int(options.get("parallel", 1)))
    cmd: list[str] = [exe, "-pa", str(pa)]

    sensitivity = str(options.get("sensitivity", "default"))
    if sensitivity == "sensitive":
        cmd.append("-s")
    elif sensitivity == "quick":
        cmd.append("-q")
    elif sensitivity == "rush":
        cmd.append("-qq")

    engine = str(options.get("engine", "") or "").strip()
    if engine:
        cmd.extend(["-e", engine])

    if options.get("use_div", True):
        cmd.extend(["-div", str(int(options.get("div", 40)))])

    lib_mode = str(options.get("library_mode", "lib"))
    if lib_mode == "species":
        species = str(options.get("species", "")).strip()
        if not species:
            raise ValueError("Species name is required for -species mode.")
        cmd.extend(["-species", species])
    else:
        lib = options.get("library_fasta")
        if not lib:
            raise ValueError("Library FASTA is required for -lib mode.")
        lib_path = Path(lib).expanduser().resolve()
        if not lib_path.is_file():
            raise FileNotFoundError(f"Library FASTA not found: {lib_path}")
        cmd.extend(["-lib", str(lib_path)])
        cutoff = options.get("cutoff")
        if cutoff is not None and int(cutoff) > 0:
            cmd.extend(["-cutoff", str(int(cutoff))])

    if options.get("nolow"):
        cmd.append("-nolow")
    if options.get("norna"):
        cmd.append("-norna")
    if options.get("no_is"):
        cmd.append("-no_is")
    if options.get("noint"):
        cmd.append("-noint")
    if options.get("gff"):
        cmd.append("-gff")
    if options.get("xsmall"):
        cmd.append("-xsmall")
    if options.get("html"):
        cmd.append("-html")
    if options.get("alignments"):
        cmd.append("-a")
    if options.get("uncurated"):
        cmd.append("-uncurated")

    cmd.extend([str(query_fasta), "-dir", str(output_dir)])
    return cmd


def run_repeatmasker_on_query(
    query_fasta: Path,
    *,
    options: dict[str, Any],
    output_dir: Path | None = None,
    cancel: JobCancellation | None = None,
) -> list[str]:
    """Run RepeatMasker on one query FASTA; return created artifact paths."""
    exe = repeatmasker_executable()
    if not exe:
        raise RuntimeError(
            "RepeatMasker is not on PATH. Install it (e.g. conda install -c bioconda "
            "repeatmasker) or set your conda env under File → Setting → External tools."
        )

    query_fasta = query_fasta.expanduser().resolve()
    if not query_fasta.is_file():
        raise FileNotFoundError(f"Query FASTA not found: {query_fasta}")

    work = (
        output_dir.expanduser().resolve()
        if output_dir is not None
        else query_fasta.parent / f"{query_fasta.stem}_RepeatMasker"
    )
    work.mkdir(parents=True, exist_ok=True)

    # RepeatMasker FastaDB rejects ids longer than 50 characters.
    short_query = work / f"{query_fasta.stem}_query_short.fa"
    _sanitize_fasta_ids(query_fasta, short_query, "q")

    run_opts = dict(options)
    if run_opts.get("library_mode", "lib") == "lib" and run_opts.get("library_fasta"):
        lib_src = Path(run_opts["library_fasta"]).expanduser().resolve()
        short_lib = work / "lib_short.fa"
        _sanitize_fasta_ids(lib_src, short_lib, "L")
        run_opts["library_fasta"] = str(short_lib)

    cmd = build_repeatmasker_cmd(
        exe=exe,
        query_fasta=short_query,
        output_dir=work,
        options=run_opts,
    )
    try:
        proc = run_subprocess_with_cancel(
            cmd,
            cancel=cancel,
            cwd=work,
            check=False,
        )
    except JobCancelledError:
        raise
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip()
        low = detail.lower()
        hint = ""
        if "trf" in low and ("not found" in low or "no such" in low):
            hint = (
                "\n\nRepeatMasker needs TRF on PATH "
                "(conda install -c bioconda trf)."
            )
        elif "rmblast" in low or "makeblastdb" in low or "search engine" in low:
            hint = (
                "\n\nSearch engine (rmblast) may be missing "
                "(conda install -c bioconda rmblast)."
            )
        elif "configure" in low or "no repeat" in low or "libraries" in low:
            hint = (
                "\n\nRepeatMasker may need configuration, or -species is unknown "
                "to your FamDB install."
            )
        raise RuntimeError(
            f"RepeatMasker failed (exit {proc.returncode}).\n\n{detail[-2000:]}{hint}"
        )

    created: list[str] = []
    preferred = work / f"{short_query.name}.out"
    if preferred.is_file():
        created.append(str(preferred.resolve()))
    else:
        outs = sorted(work.glob("*.out"))
        if not outs:
            raise RuntimeError("RepeatMasker finished but no .out file was found.")
        created.append(str(outs[0].resolve()))

    for pattern in ("*.gff", "*.tbl", "*.masked", "*.align", "*.html"):
        for path in sorted(work.glob(pattern)):
            resolved = str(path.resolve())
            if resolved not in created:
                created.append(resolved)
    return created
