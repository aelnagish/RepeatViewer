"""RepeatViewer plugin: run RepeatMasker on selected FASTA or the defined genome."""

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
        "Run RepeatMasker on selected FASTA (-lib or -species), "
        "or mask the defined genome with the selected library"
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
        dialog = RepeatMaskerParamsDialog(parent)
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
            return "Select at least one FASTA file."
        if not repeatmasker_executable_available():
            return (
                "RepeatMasker is not on PATH. Install it "
                "(e.g. conda install -c bioconda repeatmasker) or configure "
                "File → Setting → External tools."
            )
        opts = getattr(self, "_job_options", {}) or {}
        if opts.get("target") == "genome":
            if not project.genome_file or not project.genome_file.is_file():
                return "Define a genome FASTA first (File → Define genome)."
            if len(files) != 1:
                return "Select exactly one library FASTA when masking the defined genome."
        return None

    def run(
        self,
        file_path: Path,
        project: RepeatViewerProject,
        **job_options,
    ) -> list[str]:
        opts = dict(job_options)
        path = file_path.expanduser().resolve()

        if opts.get("target") == "genome":
            genome = project.genome_file
            if genome is None or not genome.is_file():
                raise RuntimeError("Defined genome FASTA is missing.")
            run_opts = dict(opts)
            run_opts["library_mode"] = "lib"
            run_opts["library_fasta"] = str(path)
            out_dir = path.parent / f"{path.stem}_vs_{genome.stem}_RepeatMasker"
            return run_repeatmasker_on_query(
                genome.expanduser().resolve(),
                options=run_opts,
                output_dir=out_dir,
            )

        out_dir = path.parent / f"{path.stem}_RepeatMasker"
        return run_repeatmasker_on_query(path, options=opts, output_dir=out_dir)
