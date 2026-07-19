from __future__ import annotations

from dataclasses import fields
from pathlib import Path

from PySide6.QtWidgets import QWidget

from teviewer.plugins.api import TEviewerPlugin
from teviewer.project import TEviewerProject

from teviewer_minimap2.dialog import Minimap2ParamsDialog
from teviewer_minimap2.executable import minimap2_executable
from teviewer_minimap2.options import Minimap2Options
from teviewer_minimap2.runner import run_minimap2


class Minimap2Plugin(TEviewerPlugin):
    id = "minimap2"
    label = "Minimap2"
    tooltip = "Map selected FASTA to a reference genome with minimap2"
    fasta_only = True
    supports_parallel = True
    picks_files_in_dialog = True

    def validate(self, files: list[Path], project: TEviewerProject) -> str | None:
        try:
            minimap2_executable()
        except (FileNotFoundError, PermissionError) as exc:
            return str(exc)
        if not files:
            return "Add at least one query FASTA file."
        reference = self._reference_path(project)
        if reference is None or not reference.is_file():
            return "Choose a reference genome FASTA."
        return None

    def prepare_job(
        self,
        parent: QWidget,
        files: list[Path],
        project: TEviewerProject,
    ) -> bool:
        dialog = Minimap2ParamsDialog(
            parent,
            project=project,
            initial_queries=files,
        )
        if dialog.exec() != dialog.DialogCode.Accepted:
            return False
        opts = dialog.options()
        self._job_options = opts.__dict__
        self._job_options["query_files"] = [str(path) for path in dialog.query_files()]
        return True

    def run(self, file_path: Path, project: TEviewerProject, **job_options) -> list[str]:
        allowed = {field.name for field in fields(Minimap2Options)}
        opts = Minimap2Options(
            **{key: value for key, value in job_options.items() if key in allowed}
        )
        reference = self._resolve_reference(project, opts)
        outputs = run_minimap2(file_path.resolve(), reference, opts)
        return [str(path.resolve()) for path in outputs]

    def _reference_path(self, project: TEviewerProject) -> Path | None:
        opts = getattr(self, "_job_options", {})
        raw = opts.get("reference_fasta") or project.genome_file
        if not raw:
            return None
        return Path(str(raw)).expanduser()

    def _resolve_reference(self, project: TEviewerProject, opts: Minimap2Options) -> Path:
        raw = opts.reference_fasta.strip() or (
            str(project.genome_file) if project.genome_file else ""
        )
        if not raw:
            raise RuntimeError("No reference FASTA configured.")
        reference = Path(raw).expanduser().resolve()
        if not reference.is_file():
            raise FileNotFoundError(reference)
        return reference
