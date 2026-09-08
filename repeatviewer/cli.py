"""Command-line entry for RepeatViewer."""

from __future__ import annotations

from repeatviewer.qt_runtime import configure_qt_runtime

configure_qt_runtime()

from pathlib import Path
from typing import Optional

import click

from repeatviewer import __version__
from repeatviewer.app import run_app
from repeatviewer.views.settings_dialog import SettingsDialog, save_workspace


@click.command()
@click.option(
    "--input-dir",
    "-i",
    type=click.Path(exists=True, file_okay=False, path_type=Path),
    default=None,
    help="Input workspace directory (e.g. TEtrimmer_for_proof_curation or any folder).",
)
@click.option(
    "--genome",
    "-g",
    type=click.Path(exists=True, dir_okay=False, path_type=Path),
    default=None,
    help="Genome FASTA.",
)
@click.version_option(__version__, prog_name="RepeatViewer")
def main(
    input_dir: Optional[Path],
    genome: Optional[Path],
) -> None:
    """
    Launch RepeatViewer (PySide6).

    Example:

        python RepeatViewer.py -i ./TEtrimmer_for_proof_curation -g genome.fa

    Tool outputs are written next to each input FASTA; a separate output directory
    is not required.
    """
    if input_dir is not None or genome is not None:
        project = SettingsDialog.load_project_from_settings()
        if input_dir is not None:
            project.input_dir = input_dir
        if genome is not None:
            project.genome_file = genome
        dialog = SettingsDialog(project)
        project = dialog.apply_to_project()
        roots = []
        if project.input_dir is not None:
            roots.append(project.input_dir)
        save_workspace(project, root_folders=roots, last_folder=project.input_dir)

    exit_code = run_app()
    raise SystemExit(exit_code)


if __name__ == "__main__":
    main()
