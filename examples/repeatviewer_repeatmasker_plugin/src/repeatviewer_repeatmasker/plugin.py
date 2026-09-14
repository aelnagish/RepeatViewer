"""RepeatViewer plugin: run RepeatMasker on the defined/chosen genome."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QWidget

from repeatviewer.backends.repeat_quantify import repeatmasker_executable_available
from repeatviewer.plugins.api import RepeatViewerPlugin
from repeatviewer.project import RepeatViewerProject

from repeatviewer_repeatmasker.dialog import RepeatMaskerParamsDialog
from repeatviewer_repeatmasker.runner import run_repeatmasker_on_query


class RepeatMaskerPlugin(RepeatViewerPlugin):
    id = "repeatmasker"
    label = "RepeatMasker"
    tooltip = (
        "Mask a genome FASTA with RepeatMasker using the selected TE library (-lib) "
        "or a -species / clade database"
    )
    fasta_only = True
    # One RM job can be heavy; keep sequential by default.
    supports_parallel = False

    def prepare_job(
        self,
        parent: QWidget,
        files: list[Path],
        project: RepeatViewerProject,
    ) -> bool:
        dialog = RepeatMaskerParamsDialog(
            parent,
            selected_files=files,
            genome_file=project.genome_file,
        )
        if dialog.exec() != dialog.DialogCode.Accepted:
            return False
        self._job_options = dialog.options()
        return True

    def validate(
        self,
        files: list[Path],
        project: RepeatViewerProject,
    ) -> str | None:
        if not files:
            return "Select at least one TE library FASTA file."
        if not repeatmasker_executable_available():
            return (
                "RepeatMasker is not on PATH. Install it "
                "(e.g. conda install -c bioconda repeatmasker) or configure "
                "File → Setting → External tools."
            )
        opts = getattr(self, "_job_options", {}) or {}
        genome_raw = str(opts.get("genome_fasta") or "").strip()
        genome = Path(genome_raw).expanduser() if genome_raw else project.genome_file
        if genome is None or not Path(genome).is_file():
            return "Choose a genome FASTA in the RepeatMasker dialog (or File → Define genome)."
        if opts.get("library_mode", "lib") == "lib" and len(files) < 1:
            return "Select a TE library FASTA for -lib."
        return None

    def run(
        self,
        file_path: Path,
        project: RepeatViewerProject,
        **job_options,
    ) -> list[str]:
        cancel = job_options.pop("cancel", None)
        opts = dict(job_options)
        lib_path = file_path.expanduser().resolve()

        genome_raw = str(opts.get("genome_fasta") or "").strip()
        if genome_raw:
            genome = Path(genome_raw).expanduser().resolve()
        elif project.genome_file is not None:
            genome = project.genome_file.expanduser().resolve()
        else:
            raise RuntimeError("Genome FASTA is missing.")
        if not genome.is_file():
            raise RuntimeError(f"Genome FASTA not found: {genome}")

        run_opts = dict(opts)
        run_opts["target"] = "genome"
        if run_opts.get("library_mode", "lib") == "lib":
            # Selected panel file is always the -lib TE library.
            run_opts["library_fasta"] = str(lib_path)
            out_dir = lib_path.parent / f"{lib_path.stem}_vs_{genome.stem}_RepeatMasker"
        else:
            out_dir = genome.parent / f"{genome.stem}_RepeatMasker"

        return run_repeatmasker_on_query(
            genome,
            options=run_opts,
            output_dir=out_dir,
            cancel=cancel,
        )
