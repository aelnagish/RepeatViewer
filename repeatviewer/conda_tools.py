"""Resolve optional conda env for external CLI tools (RepeatMasker, FastGA, …).

Users set an env **name** (``tetrimmer``) or an env **path**
(``~/miniconda3/envs/tetrimmer``) in Preferences. RepeatViewer prepends that
env's ``bin`` to ``PATH`` so GUI launches find conda-installed tools.
"""

from __future__ import annotations

import logging
import os
import shutil
import subprocess
from pathlib import Path

from PySide6.QtCore import QSettings

from repeatviewer.app_preferences import APP_SETTINGS_GROUP

logger = logging.getLogger(__name__)

_CONDA_PREF_KEY = "conda_tools_env"
_PATH_MARKER = "REPEATVIEWER_CONDA_BIN"
_original_path: str | None = None


def load_conda_tools_env() -> str:
    """Raw preference: env name, env prefix path, or empty."""
    s = QSettings()
    s.beginGroup(APP_SETTINGS_GROUP)
    value = s.value(_CONDA_PREF_KEY, "")
    s.endGroup()
    if value is None:
        return ""
    return str(value).strip()


def save_conda_tools_env(value: str) -> None:
    text = (value or "").strip()
    s = QSettings()
    s.beginGroup(APP_SETTINGS_GROUP)
    if text:
        s.setValue(_CONDA_PREF_KEY, text)
    else:
        s.remove(_CONDA_PREF_KEY)
    s.endGroup()


def _candidate_conda_roots() -> list[Path]:
    roots: list[Path] = []
    seen: set[str] = set()

    def add(path: Path | None) -> None:
        if path is None:
            return
        try:
            resolved = path.expanduser().resolve()
        except OSError:
            return
        key = str(resolved)
        if key in seen:
            return
        seen.add(key)
        roots.append(resolved)

    conda_exe = os.environ.get("CONDA_EXE") or shutil.which("conda") or shutil.which("mamba")
    if conda_exe:
        add(Path(conda_exe).resolve().parent.parent)

    for key in ("CONDA_PREFIX", "MAMBA_ROOT_PREFIX"):
        raw = os.environ.get(key)
        if raw:
            prefix = Path(raw).expanduser()
            # Active env: …/envs/name → root is parent of envs; base env → itself.
            if prefix.name != "envs" and prefix.parent.name == "envs":
                add(prefix.parent.parent)
            else:
                add(prefix)

    home = Path.home()
    for name in (
        "miniconda3",
        "anaconda3",
        "mambaforge",
        "miniforge3",
        "micromamba",
    ):
        add(home / name)
        add(home / "opt" / name)

    add(Path("/opt/homebrew/Caskroom/miniconda/base"))
    add(Path("/usr/local/Caskroom/miniconda/base"))
    add(Path("/opt/miniconda3"))
    add(Path("/opt/anaconda3"))

    # `conda info --base` when conda is on PATH.
    if conda_exe:
        try:
            result = subprocess.run(
                [conda_exe, "info", "--base"],
                capture_output=True,
                text=True,
                timeout=8,
                check=False,
            )
            if result.returncode == 0 and result.stdout.strip():
                add(Path(result.stdout.strip()))
        except (OSError, subprocess.SubprocessError):
            pass

    return roots


def resolve_conda_env_prefix(spec: str | None = None) -> Path | None:
    """Return the conda env prefix directory, or None if unresolved."""
    text = (spec if spec is not None else load_conda_tools_env()).strip()
    if not text:
        return None

    path = Path(text).expanduser()
    # Explicit path: env root, or …/bin, or …/envs/name.
    if path.is_absolute() or text.startswith(("~", "/", ".")):
        try:
            path = path.resolve()
        except OSError:
            return None
        if path.name == "bin" and path.is_dir():
            return path.parent
        if (path / "bin").is_dir():
            return path
        if path.is_dir() and path.parent.name == "envs":
            return path
        return None

    # Env name: search known conda roots.
    name = text
    for root in _candidate_conda_roots():
        for candidate in (
            root / "envs" / name,
            root / name,  # rare layouts
        ):
            if (candidate / "bin").is_dir():
                return candidate
    return None


def resolve_conda_env_bin(spec: str | None = None) -> Path | None:
    """Return ``…/envs/<name>/bin`` when the preference resolves."""
    prefix = resolve_conda_env_prefix(spec)
    if prefix is None:
        return None
    bin_dir = prefix / "bin"
    return bin_dir if bin_dir.is_dir() else None


def tool_search_path(spec: str | None = None) -> str:
    """PATH string with conda env ``bin`` prepended when configured."""
    current = os.environ.get("PATH", "")
    bin_dir = resolve_conda_env_bin(spec)
    if bin_dir is None:
        return current
    bin_text = str(bin_dir)
    parts = [p for p in current.split(os.pathsep) if p and p != bin_text]
    return os.pathsep.join([bin_text, *parts])


def which_tool(cmd: str, *, spec: str | None = None) -> str | None:
    """Like ``shutil.which``, preferring the configured conda env."""
    return shutil.which(cmd, path=tool_search_path(spec))


def apply_conda_tools_path(spec: str | None = None) -> Path | None:
    """Prepend the configured conda ``bin`` to ``os.environ['PATH']``.

    Returns the bin directory when applied, else None. Safe to call repeatedly.
    Restores the original PATH when the preference is cleared.
    """
    global _original_path
    if _original_path is None:
        _original_path = os.environ.get("PATH", "")

    base_parts = [p for p in _original_path.split(os.pathsep) if p]

    text = (spec if spec is not None else load_conda_tools_env()).strip()
    if not text:
        os.environ["PATH"] = os.pathsep.join(base_parts) if base_parts else _original_path
        os.environ.pop(_PATH_MARKER, None)
        return None

    bin_dir = resolve_conda_env_bin(text)
    if bin_dir is None:
        logger.warning("Conda tools env not found for %r", text)
        os.environ["PATH"] = os.pathsep.join(base_parts) if base_parts else _original_path
        os.environ.pop(_PATH_MARKER, None)
        return None

    bin_text = str(bin_dir)
    new_parts = [bin_text, *[p for p in base_parts if p != bin_text]]
    os.environ["PATH"] = os.pathsep.join(new_parts)
    os.environ[_PATH_MARKER] = bin_text
    logger.info("Using conda tools PATH prefix: %s", bin_text)
    return bin_dir


def describe_conda_tools_status(spec: str | None = None) -> str:
    """Short status string for the Preferences UI."""
    text = (spec if spec is not None else load_conda_tools_env()).strip()
    if not text:
        return "Not set — using the default system PATH only."
    bin_dir = resolve_conda_env_bin(text)
    if bin_dir is None:
        return f"Could not resolve “{text}”. Use an env name (e.g. tetrimmer) or a full path to the env folder."
    samples = []
    for name, label in (
        ("RepeatMasker", "RepeatMasker"),
        ("FastGA", "FastGA"),
        ("blastn", "BLASTn"),
    ):
        found = shutil.which(name, path=str(bin_dir))
        if found:
            samples.append(label)
    if samples:
        sample_text = "like " + ", ".join(samples)
    else:
        sample_text = "no common tools detected in bin/"
    return f"Using {bin_dir} ({sample_text})."
