"""Command-line entry for RepeatViewer."""

from __future__ import annotations

import sys


def _maybe_run_flexidot_entrypoint() -> None:
    """Frozen-app helper: ``RepeatViewer.app --repeatviewer-flexidot …`` → Flexidot."""
    marker = "--repeatviewer-flexidot"
    if marker not in sys.argv:
        return
    # Drop the marker (and anything before it that is the app path already in argv[0]).
    idx = sys.argv.index(marker)
    sys.argv = [sys.argv[0], *sys.argv[idx + 1 :]]
    # Re-entering the .app binary would otherwise show a second Dock icon.
    from repeatviewer.ui.platform import hide_from_macos_dock

    hide_from_macos_dock()
    from flexidot.app import main as flexidot_main

    raise SystemExit(flexidot_main())


def _maybe_run_import_check() -> None:
    """Frozen-app diagnostic: ``--repeatviewer-import-check`` → print key deps and exit."""
    marker = "--repeatviewer-import-check"
    if marker not in sys.argv:
        return
    checks = (
        ("orfipy_core", "orfipy_core"),
        ("pyhmmer", "pyhmmer"),
        ("pyhmmer.plan7", "pyhmmer.plan7"),
        ("pyhmmer.hmmer", "pyhmmer.hmmer"),
        ("pytrf", "pytrf"),
    )
    lines: list[str] = []
    for label, modname in checks:
        try:
            mod = __import__(modname, fromlist=["*"])
            where = getattr(mod, "__file__", "?")
            lines.append(f"OK  {label}: {where}")
        except Exception as exc:  # noqa: BLE001
            lines.append(f"FAIL {label}: {type(exc).__name__}: {exc}")
    sys.stderr.write("\n".join(lines) + "\n")
    raise SystemExit(0 if all(line.startswith("OK") for line in lines) else 1)


_maybe_run_import_check()
_maybe_run_flexidot_entrypoint()

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
