"""Runner for the external Flexidot-Classifier tool.

Flexidot-Classifier (https://github.com/tudipffmgt/Flexidot-Classifier) is a
computer-vision heuristic that reads Flexidot *self* dotplots (PDF/PNG) and
labels each panel as an LTR retrotransposon, tandem repeat, or satellite. We
run it as an *external* program in its own environment (the same isolation
model used for Flexidot itself) rather than importing it, so none of its code
or its heavy dependencies (opencv, PyMuPDF) live inside TEviewer.

The tool is a folder-driven script: it reads ``<cwd>/input`` and writes
``<cwd>/output`` (all paths are relative to the working directory). We honour
that contract by staging a temporary working directory per run.

------------------------------------------------------------------------------
TODO (calibration test): the classifier's crop margins and contour thresholds
are calibrated to the *authors'* Flexidot rendering (750x750 panels, fixed
37/15 px crops, 200 dpi). Before trusting results, generate a handful of self
dotplots from TEviewer, run this on them, and eyeball ``*_annotated.png`` to
confirm the boxes land on the right panels. If they are off, we either match
the authors' render settings when producing the plots or expose the Config
knobs (edge_threshold, crop_bottom_left/top_right, plots_per_page).

TODO (write-back option, deferred): right now we only surface the results
(annotated image + txt lists). The richer option is to map each detected plot
index back to its source sequence and tag the sequence with the suggested
class. That mapping is reliable only when the classifier's panel order matches
Flexidot's grid order; the robust way to get it is one-PNG-per-sequence input
(-m 0 -f png), where the index is just the filename. Revisit when we want
tagging.
------------------------------------------------------------------------------
"""

from __future__ import annotations

import os
import shutil
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path

from teviewer.external_process import run_logged_process

# Class labels the tool emits, in escalating-density order (see its README).
CLASSIFIER_CLASSES = ("ltr", "tandem", "satellite")

# Input files the tool accepts. It refuses a folder mixing the two.
_PDF_EXTS = {".pdf"}
_PNG_EXTS = {".png"}


def default_flexidot_classifier_home() -> Path | None:
    """Locate a Flexidot-Classifier clone (env var, then common locations).

    The marker file is ``main.py`` at the repository root.
    """
    env = os.environ.get("FLEXIDOT_CLASSIFIER_HOME", "").strip()
    if env:
        path = Path(env).expanduser()
        if (path / "main.py").is_file():
            return path
    candidates = [
        Path.home() / "Documents" / "Flexidot-Classifier",
        Path.home() / "Documents" / "GitHub" / "Flexidot-Classifier",
        Path.home() / "GitHub" / "Flexidot-Classifier",
        Path.home() / "Flexidot-Classifier",
    ]
    for path in candidates:
        if (path / "main.py").is_file():
            return path
    return None


def flexidot_classifier_python(home: Path) -> str:
    """Pick the interpreter that has the classifier's dependencies.

    Prefer a virtual environment created next to / inside the clone (that is
    where its opencv/PyMuPDF deps should live); fall back to the interpreter
    running TEviewer as a last resort.
    """
    home = home.expanduser().resolve()
    for rel in (".venv", "venv", "env"):
        for exe in ("bin/python", "bin/python3", "Scripts/python.exe"):
            candidate = home / rel / exe
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return str(candidate)
    # Last resort: TEviewer's own interpreter. Only works if the user installed
    # the classifier's requirements into it.
    return sys.executable


@dataclass
class ClassifierResult:
    """Outcome of one classifier run."""

    output_dir: Path
    counts: dict[str, int] = field(default_factory=dict)
    indices: dict[str, list[int]] = field(default_factory=dict)
    annotated_images: list[Path] = field(default_factory=list)
    txt_files: list[Path] = field(default_factory=list)

    @property
    def total_detected(self) -> int:
        return sum(self.counts.values())


def _classify_input_kind(plot_files: list[Path]) -> str:
    """Return 'pdf' or 'png'; raise if mixed or unsupported (tool's own rule)."""
    exts = {p.suffix.lower() for p in plot_files}
    if exts <= _PDF_EXTS:
        return "pdf"
    if exts <= _PNG_EXTS:
        return "png"
    if exts & _PDF_EXTS and exts & _PNG_EXTS:
        raise ValueError(
            "Flexidot-Classifier cannot mix PDF and PNG inputs; feed it one kind."
        )
    raise ValueError(
        f"Unsupported input for Flexidot-Classifier: {sorted(exts)} "
        "(expects .pdf or .png)."
    )


def run_flexidot_classifier(
    plot_files: list[Path],
    *,
    home: Path | None = None,
    python_exe: str | None = None,
    work_root: Path | None = None,
    keep_workdir: bool = True,
    timeout: float | None = 900.0,
) -> ClassifierResult:
    """Run Flexidot-Classifier on a set of Flexidot self-dotplot files.

    Parameters
    ----------
    plot_files:
        Flexidot output images (all ``.pdf`` OR all ``.png``). These are the
        self-mode (``-m 0``) plots the classifier is designed to read.
    home:
        Path to a Flexidot-Classifier clone; auto-detected if omitted.
    python_exe:
        Interpreter with the classifier's deps; auto-detected if omitted.
    work_root:
        Where to create the temporary ``input/``+``output/`` working dir.
    keep_workdir:
        Leave the working dir (with annotated images) in place for the caller
        to display. When False it is deleted after parsing.
    """
    plot_files = [Path(p).expanduser().resolve() for p in plot_files]
    if not plot_files:
        raise ValueError("No plot files supplied to Flexidot-Classifier.")
    for path in plot_files:
        if not path.is_file():
            raise FileNotFoundError(f"Plot file not found: {path}")
    kind = _classify_input_kind(plot_files)

    if home is None:
        home = default_flexidot_classifier_home()
    if home is None:
        raise RuntimeError(
            "Flexidot-Classifier not found. Clone "
            "https://github.com/tudipffmgt/Flexidot-Classifier and set "
            "FLEXIDOT_CLASSIFIER_HOME to its folder."
        )
    home = Path(home).expanduser().resolve()
    main_py = home / "main.py"
    if not main_py.is_file():
        raise RuntimeError(f"main.py not found in Flexidot-Classifier home: {home}")

    if python_exe is None:
        python_exe = flexidot_classifier_python(home)

    workdir = Path(
        tempfile.mkdtemp(prefix="flexidot_classifier_", dir=str(work_root) if work_root else None)
    )
    input_dir = workdir / "input"
    input_dir.mkdir(parents=True, exist_ok=True)
    for path in plot_files:
        shutil.copy2(path, input_dir / path.name)

    # The script keys off relative "input"/"output"; run it with cwd=workdir.
    run_logged_process(
        [python_exe, str(main_py)],
        tool="Flexidot-Classifier",
        inputs=[str(p) for p in plot_files],
        cwd=workdir,
        check=True,
        timeout=timeout,
    )

    result = _parse_output(workdir / "output" / "final")
    if not keep_workdir:
        shutil.rmtree(workdir, ignore_errors=True)
    return result


def _parse_output(final_dir: Path) -> ClassifierResult:
    """Read the tool's ``output/final`` folder into a ClassifierResult."""
    result = ClassifierResult(output_dir=final_dir)
    if not final_dir.is_dir():
        return result
    for cls in CLASSIFIER_CLASSES:
        txt = final_dir / f"{cls}.txt"
        indices: list[int] = []
        if txt.is_file():
            result.txt_files.append(txt)
            for line in txt.read_text().splitlines():
                line = line.strip()
                if not line:
                    continue
                try:
                    indices.append(int(line))
                except ValueError:
                    continue
        result.indices[cls] = indices
        result.counts[cls] = len(indices)
    result.annotated_images = sorted(final_dir.glob("*_annotated.png"))
    return result
