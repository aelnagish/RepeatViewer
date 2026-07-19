"""Open FASTA in AliView or system viewer; MAFFT alignment (TEtrimmerGUI behaviour)."""

from __future__ import annotations

import logging
import os
import platform
import subprocess
import sys
from pathlib import Path

from teviewer.backends.gui_functions import get_original_file_path
from teviewer.backends.sequence_aligner import run_alignment

OS_TYPE = platform.system()


def aliview_paths() -> tuple[Path | None, Path | None]:
    """
    Return (launcher_script, jar_path) for bundled AliView under teviewer/backends/aliview/.
    """
    base = Path(get_original_file_path()) / "aliview"
    launcher = base / "aliview"
    jar_path = base / "aliview.jar"
    if OS_TYPE == "Windows":
        jar_path = base / "aliview.jar"
        if jar_path.is_file():
            return None, jar_path
        return None, None
    if launcher.is_file():
        return launcher, jar_path if jar_path.is_file() else None
    return None, None


def open_fasta(path: Path, *, use_system_viewer: bool = False) -> None:
    """Open a FASTA/MSA file in AliView or the OS default application."""
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()
    if suffix not in {".fa", ".fasta", ".fas"}:
        raise ValueError(f"Not a FASTA file: {path}")

    if use_system_viewer:
        _open_with_system(path)
        return

    launcher, jar_path = aliview_paths()
    if OS_TYPE == "Windows":
        if jar_path:
            subprocess.Popen(
                ["java", "-jar", str(jar_path), str(path)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return
        raise FileNotFoundError("AliView JAR not found in teviewer/backends/aliview/")

    if launcher and launcher.is_file():
        subprocess.Popen(
            [str(launcher), str(path)],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return

    logging.warning("Bundled AliView not found; falling back to system default.")
    _open_with_system(path)


def _open_with_system(path: Path) -> None:
    if OS_TYPE == "Darwin":
        subprocess.Popen(["open", str(path)])
    elif OS_TYPE == "Windows":
        os.startfile(str(path))  # type: ignore[attr-defined]
    else:
        subprocess.Popen(["xdg-open", str(path)])


def run_mafft_alignment(input_fasta: Path, num_threads: int = 5) -> Path:
    """Run MAFFT (legacy alias). Prefer ``sequence_aligner.run_alignment``."""
    return run_alignment(input_fasta, backend="mafft", num_threads=num_threads)
