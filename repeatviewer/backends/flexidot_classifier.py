"""Runner for the external Flexidot-Classifier tool.

Flexidot-Classifier (https://github.com/tudipffmgt/Flexidot-Classifier) is a
computer-vision heuristic that reads Flexidot *self* dotplots (PDF/PNG) and
labels each panel as an LTR retrotransposon, tandem repeat, or satellite. We
run it as an *external* program in its own environment (the same isolation
model used for Flexidot itself) rather than importing it, so none of its code
or its heavy dependencies (opencv, PyMuPDF) live inside RepeatViewer.

The tool is a folder-driven script: it reads ``<cwd>/input`` and writes
``<cwd>/output`` (all paths are relative to the working directory). We honour
that contract by staging a temporary working directory per run.

------------------------------------------------------------------------------
TODO (calibration test): the classifier's crop margins and contour thresholds
are calibrated to the *authors'* Flexidot rendering (750x750 panels, fixed
37/15 px crops, 200 dpi). Before trusting results, generate a handful of self
dotplots from RepeatViewer, run this on them, and eyeball ``*_annotated.png`` to
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

from repeatviewer.external_process import run_logged_process

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
    running RepeatViewer as a last resort.
    """
    home = home.expanduser().resolve()
    for rel in (".venv", "venv", "env"):
        for exe in ("bin/python", "bin/python3", "Scripts/python.exe"):
            candidate = home / rel / exe
            if candidate.is_file() and os.access(candidate, os.X_OK):
                return str(candidate)
    # Last resort: RepeatViewer's own interpreter. Only works if the user installed
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
    _classify_input_kind(plot_files)  # validates input; raises on mixed/unsupported

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
    proc = run_logged_process(
        [python_exe, str(main_py)],
        tool="Flexidot-Classifier",
        inputs=[str(p) for p in plot_files],
        cwd=workdir,
        check=False,  # inspect ourselves so we can surface the real traceback
        timeout=timeout,
    )
    if proc.returncode != 0:
        detail = (proc.stderr or proc.stdout or "").strip()
        low = detail.lower()
        hint = ""
        if "modulenotfounderror" in low or "no module named" in low:
            hint = (
                "\n\nThe classifier's Python is missing its dependencies "
                "(numpy, opencv-python, Pillow, PyMuPDF). It ran with:\n"
                f"  {python_exe}\n"
                "Set up its own environment so it doesn't use RepeatViewer's Python:\n"
                f"  cd {home}\n"
                "  python3 -m venv .venv\n"
                "  ./.venv/bin/pip install -r requirements.txt\n"
                "(Or use the button's automatic install.)"
            )
        raise RuntimeError(
            f"Flexidot-Classifier failed (exit {proc.returncode}).\n\n"
            f"{detail[-2000:]}{hint}"
        )

    result = _parse_output(workdir / "output" / "final")
    if not keep_workdir:
        shutil.rmtree(workdir, ignore_errors=True)
    return result


INSTALL_REPO_URL = "https://github.com/tudipffmgt/Flexidot-Classifier.git"

# Copy-pasteable manual steps, shown as the fallback to auto-install.
MANUAL_INSTALL_COMMANDS = (
    "git clone https://github.com/tudipffmgt/Flexidot-Classifier.git"
    " ~/Documents/Flexidot-Classifier\n"
    "cd ~/Documents/Flexidot-Classifier\n"
    "python3.12 -m venv .venv\n"
    "./.venv/bin/pip install -r requirements.txt"
)


def default_install_dir() -> Path:
    """Where auto-install places the clone (matches the auto-detect path)."""
    return Path.home() / "Documents" / "Flexidot-Classifier"


def _python_for_venv() -> str | None:
    """Find a system Python (prefer 3.12+) to build the classifier's venv."""
    for name in ("python3.12", "python3.13", "python3.11", "python3.10", "python3"):
        found = shutil.which(name)
        if found:
            return found
    return None


def _venv_python(venv_dir: Path) -> Path:
    posix = venv_dir / "bin" / "python"
    win = venv_dir / "Scripts" / "python.exe"
    return posix if posix.is_file() else win


def install_flexidot_classifier(
    dest: Path | None = None,
    *,
    progress=None,
    cancel=None,
) -> Path:
    """Clone the classifier and provision its virtual environment.

    This only *fetches* the public repository (it does not redistribute the
    code), so it is not gated by the licence question that bundling would be.
    Requires ``git`` and a system Python (>=3.10) and network access.
    """

    def report(message: str) -> None:
        if progress is not None:
            progress(message)

    dest = Path(dest).expanduser() if dest else default_install_dir()

    git = shutil.which("git")
    if git is None:
        raise RuntimeError(
            "git was not found on PATH. Install git, or run the manual commands."
        )

    if (dest / "main.py").is_file():
        report(f"Repository already present at {dest}; skipping clone.")
    else:
        dest.parent.mkdir(parents=True, exist_ok=True)
        report(f"Cloning Flexidot-Classifier into {dest}…")
        run_logged_process(
            [git, "clone", INSTALL_REPO_URL, str(dest)],
            tool="Flexidot-Classifier setup",
            check=True,
            timeout=600,
        )

    venv_dir = dest / ".venv"
    venv_py = _venv_python(venv_dir)
    if not venv_py.is_file():
        system_py = _python_for_venv()
        if system_py is None:
            raise RuntimeError(
                "No system Python (>=3.10) found to create the environment."
            )
        report("Creating virtual environment (.venv)…")
        run_logged_process(
            [system_py, "-m", "venv", str(venv_dir)],
            tool="Flexidot-Classifier setup",
            check=True,
            timeout=300,
        )
        venv_py = _venv_python(venv_dir)
    if not venv_py.is_file():
        raise RuntimeError(f"Failed to create a virtual environment in {venv_dir}.")

    report("Installing dependencies (numpy, opencv, Pillow, PyMuPDF)…")
    run_logged_process(
        [str(venv_py), "-m", "pip", "install", "-r", str(dest / "requirements.txt")],
        tool="Flexidot-Classifier setup",
        check=True,
        timeout=1800,
    )
    report("Flexidot-Classifier setup complete.")
    return dest


def generate_self_dotplots(
    seq_files: list[Path],
    *,
    flexidot_options=None,
    filetype: str = "png",
    out_dir: Path | None = None,
    cancel=None,
) -> list[Path]:
    """Produce Flexidot self-mode (-m 0) plots for the classifier.

    Runs Flexidot per selected FASTA in ``png`` (or ``pdf``) — the formats the
    classifier can read, never HTML. Self-dotplots need no genome, so this
    calls ``run_flexidot`` directly rather than the genome-gated batch path.

    Each file is rendered into its *own fresh directory* so every produced plot
    is unambiguously new. (A before/after diff in the data folder silently
    misses plots because Flexidot overwrites a deterministic ``<seq>_flexidot``
    name that may already exist from an earlier run.) Returns the plot files.
    """
    from repeatviewer.backends.gui_functions import run_flexidot
    from repeatviewer.job_options import FlexidotOptions

    opts = flexidot_options or FlexidotOptions()
    ext = "pdf" if str(filetype).lower() == "pdf" else "png"
    root = Path(out_dir).expanduser().resolve() if out_dir else None
    if root is not None:
        root.mkdir(parents=True, exist_ok=True)
    created: list[Path] = []
    for index, seq_file in enumerate(seq_files):
        seq_file = Path(seq_file).expanduser().resolve()
        if root is not None:
            target = root / f"{index:03d}_{seq_file.stem}"
            target.mkdir(parents=True, exist_ok=True)
        else:
            target = Path(tempfile.mkdtemp(prefix="flexidot_selfplot_"))
        run_flexidot(
            str(seq_file),
            str(target),
            wordsize=opts.wordsize,
            mode="0",  # self-comparison — the classifier's expected input
            filetype=ext,
            lcs_shading=opts.lcs_shading,
            use_gff=opts.use_gff,
            substitution_count=opts.substitution_count,
            wobble_conversion=opts.wobble_conversion,
            max_n_percentage=opts.max_n_percentage,
            include_reverse=opts.include_reverse,
            only_vs_first=opts.only_vs_first,
            line_width=opts.line_width,
            length_scaling=False,  # -L is pairwise/all-vs-all only
            gff_source_fastas=[str(seq_file)],
            cancel=cancel,
        )
        # Fresh dir → every plot of this type is a product of this run.
        created.extend(sorted(target.glob(f"*.{ext}")))
    return created


def generate_and_classify(
    seq_files: list[Path],
    *,
    flexidot_options=None,
    home: Path | None = None,
    filetype: str = "png",
    work_root: Path | None = None,
    cancel=None,
) -> tuple[ClassifierResult, list[Path]]:
    """One-click: make self-dotplots for ``seq_files`` then classify them."""
    seq_files = [Path(p).expanduser().resolve() for p in seq_files]
    if not seq_files:
        raise ValueError("No sequence files selected to classify.")
    plots = generate_self_dotplots(
        seq_files,
        flexidot_options=flexidot_options,
        filetype=filetype,
        out_dir=work_root,
        cancel=cancel,
    )
    if not plots:
        raise RuntimeError(
            "Flexidot produced no PNG/PDF plots to classify. Check that the "
            "selected FASTA files contain sequences."
        )
    result = run_flexidot_classifier(plots, home=home)
    return result, plots


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
