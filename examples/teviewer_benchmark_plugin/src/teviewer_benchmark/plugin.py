"""TEviewer plugin wrapping TEtrimmer TE_anno_benchmarker (Sankey / confusion matrix)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtWidgets import QWidget

from teviewer.plugins.api import TEviewerPlugin
from teviewer.project import TEviewerProject

from teviewer_benchmark.dialog import BenchmarkParamsDialog, SETTINGS_GROUP
from teviewer_benchmark.tetrimmer_runner import (
    benchmark_outputs,
    default_tetrimmer_home,
    load_calculate_overlap,
    require_bedtools,
)


class TEBenchmarkPlugin(TEviewerPlugin):
    id = "te_benchmark"
    label = "TE benchmark"
    tooltip = (
        "Compare RepeatMasker .out files to a reference library "
        "(TEtrimmer TE_anno_benchmarker: Sankey, heatmap, confusion matrix)"
    )
    fasta_only = False
    file_suffixes = {".out"}
    supports_parallel = False

    def prepare_job(
        self,
        parent: QWidget,
        files: list[Path],
        project: TEviewerProject,
    ) -> bool:
        dialog = BenchmarkParamsDialog(parent)
        if dialog.exec() != dialog.DialogCode.Accepted:
            return False
        home = dialog.tetrimmer_home() or default_tetrimmer_home()
        reference = dialog.reference_out()
        if home is None:
            from PySide6.QtWidgets import QMessageBox

            QMessageBox.warning(
                parent,
                "TE benchmarking",
                "Set the TEtrimmer repository folder in the dialog, or set "
                "environment variable TETRIMMER_HOME to the repo root.",
            )
            return False
        if reference is None or not reference.is_file():
            from PySide6.QtWidgets import QMessageBox

            QMessageBox.warning(
                parent,
                "TE benchmarking",
                "Choose a reference RepeatMasker .out file (curated TE library).",
            )
            return False
        self._job_options = {
            "tetrimmer_home": home,
            "reference": reference,
            "less_classification": dialog.less_classification(),
            "debug": dialog.debug(),
        }
        return True

    def validate(
        self,
        files: list[Path],
        project: TEviewerProject,
    ) -> str | None:
        if not files:
            return "Select one or more RepeatMasker .out files (query)."
        if not project.genome_file or not project.genome_file.is_file():
            return "Set genome FASTA in File → Project settings (used for genome length)."
        try:
            require_bedtools()
        except RuntimeError as exc:
            return str(exc)
        home = default_tetrimmer_home()
        from PySide6.QtCore import QSettings

        s = QSettings()
        s.beginGroup(SETTINGS_GROUP)
        saved_home = s.value("tetrimmer_home", "")
        s.endGroup()
        if not home and not saved_home:
            return (
                "Set TETRIMMER_HOME to your TEtrimmer repo, or configure the path "
                "in the dialog when you run TE benchmark."
            )
        return None

    def run(
        self,
        file_path: Path,
        project: TEviewerProject,
        **job_options,
    ) -> list[str]:
        query = file_path.expanduser().resolve()
        if not query.is_file():
            raise FileNotFoundError(query)

        home = Path(job_options["tetrimmer_home"]).expanduser().resolve()
        reference = Path(job_options["reference"]).expanduser().resolve()
        out_dir = query.parent

        calculate_overlap = load_calculate_overlap(home)
        calculate_overlap(
            query=str(query),
            reference=str(reference),
            output_dir=str(out_dir),
            genome_dir=str(project.genome_file),
            genome_length=None,
            less_classification=bool(job_options.get("less_classification")),
            debug=bool(job_options.get("debug")),
            query2=None,
        )

        created = [str(p) for p in benchmark_outputs(query, out_dir)]
        if not created:
            raise RuntimeError("Benchmark finished but no expected output files were found.")
        return created
