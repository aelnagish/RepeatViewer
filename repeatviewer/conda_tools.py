"""Resolve optional conda env for external CLI tools (RepeatMasker, FastGA, …).

Users set an env **name** (``tetrimmer``) or an env **path** in Preferences.

* **macOS / Linux:** prepends that env's ``bin`` to ``PATH``.
* **Windows:** uses the same preference against **conda inside WSL2**. Linux
  binaries cannot sit on the Windows PATH, so external tools are launched via
  ``wsl`` with that env's ``bin`` prepended, and Windows paths in arguments are
  converted with ``wslpath``.
"""

from __future__ import annotations

import logging
import os
import re
import shlex
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Sequence

from PySide6.QtCore import QSettings

from repeatviewer.app_preferences import APP_SETTINGS_GROUP

logger = logging.getLogger(__name__)

_CONDA_PREF_KEY = "conda_tools_env"
_PATH_MARKER = "REPEATVIEWER_CONDA_BIN"
_WSL_BIN_ENV = "REPEATVIEWER_WSL_CONDA_BIN"
_original_path: str | None = None
_UNSET = object()  # sentinel: WSL bin not looked up yet
_wsl_bin_cache: str | None | object = _UNSET
_wsl_resolve_cache: dict[str, str | None] = {}
_wsl_which_cache: dict[tuple[str, str], str | None] = {}
_conda_roots_cache: list[Path] | None = None

# Login shells inherit Windows HOME/PATH and conda hook noise. Always start a
# clean bash and restore the Linux home from passwd (Windows HOME is often
# C:\Users\… which becomes C:Users… inside bash).
_WSL_HOME_PREAMBLE = (
    'linux_home=$(getent passwd "$(id -un)" | cut -d: -f6 2>/dev/null || true); '
    '[ -n "$linux_home" ] && export HOME="$linux_home"; '
)


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


def wsl_available() -> bool:
    if sys.platform != "win32":
        return False
    return shutil.which("wsl") is not None


def _normalize_linux_script(script: str) -> str:
    """Drop CR so Windows-sourced multiline bash is valid in WSL."""
    return script.replace("\r\n", "\n").replace("\r", "\n")


def _wsl_bash_argv(linux_command: str) -> list[str]:
    wrapped = _WSL_HOME_PREAMBLE + _normalize_linux_script(linux_command)
    return ["wsl", "-e", "bash", "--noprofile", "--norc", "-c", wrapped]


def _wsl_run(
    linux_command: str,
    *,
    timeout: float = 20,
) -> subprocess.CompletedProcess[str]:
    from repeatviewer.external_process import windows_hidden_process_kwargs

    return subprocess.run(
        _wsl_bash_argv(linux_command),
        capture_output=True,
        text=True,
        timeout=timeout,
        check=False,
        encoding="utf-8",
        errors="replace",
        **windows_hidden_process_kwargs(),
    )


def _wsl_first_line(result: subprocess.CompletedProcess[str]) -> str | None:
    """First nonempty stdout line. Do not require exit 0 — ``wsl``/bash -c
    sometimes still prints the path then exits 1 (``set -e``, trailing
    ``exit 1``, login-shell rc)."""
    for raw in (result.stdout or "").splitlines():
        line = raw.strip()
        if line:
            return line
    return None


def _clear_wsl_caches() -> None:
    global _wsl_bin_cache, _conda_roots_cache
    _wsl_bin_cache = _UNSET
    _wsl_resolve_cache.clear()
    _wsl_which_cache.clear()
    _conda_roots_cache = None


def _candidate_conda_roots() -> list[Path]:
    """Conda install roots used to resolve env names. Cached (``conda info --base`` is slow)."""
    global _conda_roots_cache
    if _conda_roots_cache is not None:
        return _conda_roots_cache

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

    # Prefer filesystem guesses first; only shell out if we still have no roots.
    # ``conda info --base`` can take seconds and was blocking the FastGA tab.
    if not roots and conda_exe:
        try:
            from repeatviewer.external_process import windows_hidden_process_kwargs

            result = subprocess.run(
                [conda_exe, "info", "--base"],
                capture_output=True,
                text=True,
                timeout=8,
                check=False,
                **windows_hidden_process_kwargs(),
            )
            if result.returncode == 0 and result.stdout.strip():
                add(Path(result.stdout.strip()))
        except (OSError, subprocess.SubprocessError):
            pass

    _conda_roots_cache = roots
    return roots


def resolve_conda_env_prefix(spec: str | None = None) -> Path | None:
    """Return a *native* conda env prefix (macOS/Linux/Windows native), or None."""
    text = (spec if spec is not None else load_conda_tools_env()).strip()
    if not text:
        return None

    path = Path(text).expanduser()
    if path.is_absolute() or text.startswith(("~", "/", ".", "\\\\")) or (
        len(text) > 2 and text[1] == ":"
    ):
        try:
            path = path.resolve()
        except OSError:
            return None
        if path.name == "bin" and path.is_dir():
            return path.parent
        if path.name == "Scripts" and path.is_dir():
            return path.parent
        if (path / "bin").is_dir() or (path / "Scripts").is_dir():
            return path
        if path.is_dir() and path.parent.name == "envs":
            return path
        return None

    name = text
    for root in _candidate_conda_roots():
        for candidate in (root / "envs" / name, root / name):
            if (candidate / "bin").is_dir() or (candidate / "Scripts").is_dir():
                return candidate
    return None


def resolve_conda_env_bin(spec: str | None = None) -> Path | None:
    """Native env ``bin`` / ``Scripts`` directory, or None."""
    prefix = resolve_conda_env_prefix(spec)
    if prefix is None:
        return None
    for sub in ("bin", "Scripts"):
        bin_dir = prefix / sub
        if bin_dir.is_dir():
            return bin_dir
    return None


def _unc_wsl_path_to_linux(text: str) -> str | None:
    """Map ``\\\\wsl$\\Distro\\home\\…`` / ``\\\\wsl.localhost\\…`` to a Linux path."""
    normalized = text.strip().replace("/", "\\")
    lower = normalized.lower()
    for prefix in (r"\\wsl$", r"\\wsl.localhost"):
        if not lower.startswith(prefix.lower()):
            continue
        # \\wsl$\Ubuntu\home\user\… → ['', '', 'wsl$', 'Ubuntu', 'home', …]
        parts = [p for p in normalized.split("\\") if p != ""]
        if len(parts) < 3:
            return None
        # parts[0]=wsl$/wsl.localhost, parts[1]=distro, parts[2:]=linux path
        return "/" + "/".join(parts[2:])
    return None


def resolve_wsl_conda_bin(spec: str | None = None) -> str | None:
    """Linux path to ``…/envs/<name>/bin`` inside WSL, or None."""
    if not wsl_available():
        return None
    text = (spec if spec is not None else load_conda_tools_env()).strip()
    if not text:
        return None
    if text in _wsl_resolve_cache:
        return _wsl_resolve_cache[text]
    found = _resolve_wsl_conda_bin_uncached(text)
    _wsl_resolve_cache[text] = found
    return found


def _resolve_wsl_conda_bin_uncached(text: str) -> str | None:
    unc_linux = _unc_wsl_path_to_linux(text)
    if unc_linux is not None:
        text = unc_linux

    # Explicit Linux path (or Windows path into \\wsl$\… converted by user as Linux).
    if text.startswith("/") or text.startswith("~"):
        if text.startswith("~"):
            # Expand ~ inside WSL (HOME already restored by _wsl_run).
            cmd = f'bin=$(eval echo {shlex.quote(text)}/bin); [ -d "$bin" ] && printf "%s\\n" "$bin"'
        else:
            cmd = (
                f'if [ -d {shlex.quote(text + "/bin")} ]; then printf "%s\\n" {shlex.quote(text + "/bin")}; '
                f'elif [ -d {shlex.quote(text)} ] && [ "$(basename {shlex.quote(text)})" = bin ]; '
                f"then printf '%s\\n' {shlex.quote(text)}; fi"
            )
        try:
            result = _wsl_run(cmd, timeout=15)
        except (OSError, subprocess.SubprocessError):
            return None
        return _wsl_first_line(result)

    # Env name: one-liner (no bash line-continuation), search common roots.
    # Also glob /home/*/miniconda3 in case Windows HOME leaked before preamble.
    name = shlex.quote(text)
    script = (
        f"name={name}; "
        'for root in "$HOME/miniconda3" "$HOME/anaconda3" "$HOME/mambaforge" '
        '"$HOME/miniforge3" "$HOME/micromamba" /opt/conda /opt/miniconda3 /opt/anaconda3 '
        "/home/*/miniconda3 /home/*/anaconda3 /home/*/mambaforge /home/*/miniforge3; do "
        'if [ -d "$root/envs/$name/bin" ]; then printf "%s\\n" "$root/envs/$name/bin"; exit 0; fi; '
        "done; "
        "if command -v conda >/dev/null 2>&1; then "
        "  base=$(conda info --base 2>/dev/null || true); "
        '  if [ -n "$base" ] && [ -d "$base/envs/$name/bin" ]; then '
        '    printf "%s\\n" "$base/envs/$name/bin"; exit 0; '
        "  fi; "
        "fi; "
        "exit 1"
    )
    try:
        result = _wsl_run(script, timeout=25)
    except (OSError, subprocess.SubprocessError):
        return None
    return _wsl_first_line(result)


def uses_wsl_conda(spec: str | None = None) -> bool:
    """True when Windows should run external tools through WSL conda."""
    if sys.platform != "win32" or not wsl_available():
        return False
    text = (spec if spec is not None else load_conda_tools_env()).strip()
    if not text:
        return False
    # Prefer WSL for Linux / \\wsl$ paths, or when native Windows conda resolve fails.
    if (
        text.startswith("/")
        or text.startswith("~")
        or _unc_wsl_path_to_linux(text) is not None
    ):
        return resolve_wsl_conda_bin(text) is not None
    if resolve_conda_env_bin(text) is not None:
        return False
    return resolve_wsl_conda_bin(text) is not None


def _active_wsl_bin() -> str | None:
    global _wsl_bin_cache
    if _wsl_bin_cache is not _UNSET:
        return _wsl_bin_cache  # type: ignore[return-value]
    cached = os.environ.get(_WSL_BIN_ENV, "").strip()
    if cached:
        _wsl_bin_cache = cached
        return cached
    bin_dir = resolve_wsl_conda_bin() if uses_wsl_conda() else None
    _wsl_bin_cache = bin_dir
    return bin_dir


def tool_search_path(spec: str | None = None) -> str:
    """PATH string with native conda env ``bin`` prepended when configured."""
    current = os.environ.get("PATH", "")
    if uses_wsl_conda(spec):
        return current
    bin_dir = resolve_conda_env_bin(spec)
    if bin_dir is None:
        return current
    bin_text = str(bin_dir)
    parts = [p for p in current.split(os.pathsep) if p and p != bin_text]
    return os.pathsep.join([bin_text, *parts])


def _wsl_which(cmd: str, *, bin_dir: str | None = None) -> str | None:
    bin_dir = bin_dir or _active_wsl_bin()
    if not bin_dir:
        return None
    key = (bin_dir, cmd)
    if key in _wsl_which_cache:
        return _wsl_which_cache[key]
    script = (
        f'export PATH={shlex.quote(bin_dir)}:$PATH; '
        f'command -v {shlex.quote(cmd)} 2>/dev/null'
    )
    try:
        result = _wsl_run(script, timeout=12)
    except (OSError, subprocess.SubprocessError):
        _wsl_which_cache[key] = None
        return None
    found = _wsl_first_line(result)
    _wsl_which_cache[key] = found
    return found


def which_tool(cmd: str, *, spec: str | None = None) -> str | None:
    """Like ``shutil.which``, preferring the configured conda env (native or WSL)."""
    if spec is None and sys.platform == "win32":
        active = _active_wsl_bin()
        if active:
            found = _wsl_which(cmd, bin_dir=active)
            return cmd if found else None
    if uses_wsl_conda(spec):
        found = _wsl_which(cmd, bin_dir=resolve_wsl_conda_bin(spec))
        # Return the bare command name; prepare_command() wraps via wsl.
        return cmd if found else None
    return shutil.which(cmd, path=tool_search_path(spec))


_WIN_PATH_RE = re.compile(r"^[A-Za-z]:[\\/]")
_EMBEDDED_WIN_PATH_RE = re.compile(r"[A-Za-z]:[\\/][^\s]*")


def _looks_like_windows_path(value: str) -> bool:
    if _WIN_PATH_RE.match(value):
        return True
    if value.startswith("\\\\"):
        return True
    return False


def _windows_drive_to_wsl(value: str) -> str | None:
    """Map ``C:\\Users\\…`` / ``C:/Users/…`` to ``/mnt/c/Users/…`` without calling wsl.

    ``wsl wslpath C:\\Users`` (no ``-e``) eats backslashes and becomes ``C:Users``.
    """
    if not _WIN_PATH_RE.match(value):
        return None
    drive = value[0].lower()
    rest = value[2:].replace("\\", "/").lstrip("/")
    return f"/mnt/{drive}/{rest}" if rest else f"/mnt/{drive}"


def _map_windows_path_to_linux(value: str) -> str | None:
    unc = _unc_wsl_path_to_linux(value)
    if unc is not None:
        return unc
    return _windows_drive_to_wsl(value)


def _to_wsl_path(value: str) -> str:
    mapped = _map_windows_path_to_linux(value)
    if mapped:
        return mapped
    if not _looks_like_windows_path(value):
        return value
    # Last resort: wslpath with -e so backslashes survive argv.
    posixish = value.replace("\\", "/")
    try:
        from repeatviewer.external_process import windows_hidden_process_kwargs

        result = subprocess.run(
            ["wsl", "-e", "wslpath", "-a", "-u", posixish],
            capture_output=True,
            text=True,
            timeout=8,
            check=False,
            encoding="utf-8",
            errors="replace",
            **windows_hidden_process_kwargs(),
        )
    except (OSError, subprocess.SubprocessError):
        return value
    line = _wsl_first_line(result)
    return line if line else value


def _rewrite_arg_for_wsl(value: str) -> str:
    """Convert a whole Windows path, or one glued into a flag (``-PC:\\…``)."""
    mapped = _map_windows_path_to_linux(value)
    if mapped:
        return mapped

    def _replace(match: re.Match[str]) -> str:
        converted = _map_windows_path_to_linux(match.group(0))
        return converted if converted else match.group(0)

    return _EMBEDDED_WIN_PATH_RE.sub(_replace, value)


def _conda_lib_dir(bin_dir: str) -> str:
    if bin_dir.endswith("/bin"):
        return bin_dir[:-3] + "lib"
    return bin_dir.rstrip("/") + "/lib"


def prepare_command(command: Sequence[object]) -> list[str]:
    """Rewrite an external command for WSL conda when needed; otherwise unchanged."""
    cmd = [str(part) for part in command]
    if not cmd or sys.platform != "win32" or not wsl_available():
        return cmd
    bin_dir = _active_wsl_bin()
    if not bin_dir and uses_wsl_conda():
        bin_dir = resolve_wsl_conda_bin()
    if not bin_dir:
        return cmd

    tool = Path(cmd[0]).name
    # Skip wrapping for Windows-native / bundled tools.
    if tool.lower().endswith((".exe", ".bat", ".cmd", ".ps1")):
        return cmd
    if not _wsl_which(tool, bin_dir=bin_dir):
        return cmd

    args = [_rewrite_arg_for_wsl(arg) for arg in cmd[1:]]
    inner = shlex.join([tool, *args])
    lib_dir = _conda_lib_dir(bin_dir)
    script = (
        f"export PATH={shlex.quote(bin_dir)}:$PATH; "
        f"export LD_LIBRARY_PATH={shlex.quote(lib_dir)}:${{LD_LIBRARY_PATH:-}}; "
        f"exec {inner}"
    )
    return _wsl_bash_argv(script)


def apply_conda_tools_path(spec: str | None = None) -> Path | str | None:
    """Apply the conda preference to the process environment.

    * Native (macOS/Linux/Windows conda): prepend ``bin``/``Scripts`` to ``PATH``.
    * Windows + WSL2 conda: set ``REPEATVIEWER_WSL_CONDA_BIN`` (no Linux path on PATH).

    Returns the bin directory (Path or Linux path str) when applied, else None.
    """
    global _original_path, _wsl_bin_cache
    if _original_path is None:
        _original_path = os.environ.get("PATH", "")
    _clear_wsl_caches()

    base_parts = [p for p in _original_path.split(os.pathsep) if p]
    text = (spec if spec is not None else load_conda_tools_env()).strip()

    # Always restore native PATH first.
    os.environ["PATH"] = os.pathsep.join(base_parts) if base_parts else _original_path
    os.environ.pop(_PATH_MARKER, None)
    os.environ.pop(_WSL_BIN_ENV, None)

    if not text:
        return None

    if sys.platform == "win32" and wsl_available():
        # Prefer WSL when native resolve fails or Linux / \\wsl$ path given.
        native = resolve_conda_env_bin(text)
        wsl_bin = resolve_wsl_conda_bin(text)
        prefer_wsl = (
            text.startswith(("/", "~"))
            or _unc_wsl_path_to_linux(text) is not None
            or native is None
        )
        if wsl_bin and prefer_wsl:
            os.environ[_WSL_BIN_ENV] = wsl_bin
            _wsl_bin_cache = wsl_bin
            logger.info("Using WSL2 conda tools PATH prefix: %s", wsl_bin)
            return wsl_bin
        if native is not None:
            bin_text = str(native)
            os.environ["PATH"] = os.pathsep.join(
                [bin_text, *[p for p in base_parts if p != bin_text]]
            )
            os.environ[_PATH_MARKER] = bin_text
            _wsl_bin_cache = None
            logger.info("Using native Windows conda tools PATH prefix: %s", bin_text)
            return native
        logger.warning("Conda tools env not found (native or WSL) for %r", text)
        return None

    bin_dir = resolve_conda_env_bin(text)
    if bin_dir is None:
        logger.warning("Conda tools env not found for %r", text)
        return None
    bin_text = str(bin_dir)
    os.environ["PATH"] = os.pathsep.join(
        [bin_text, *[p for p in base_parts if p != bin_text]]
    )
    os.environ[_PATH_MARKER] = bin_text
    logger.info("Using conda tools PATH prefix: %s", bin_text)
    return bin_dir


def describe_conda_tools_status(spec: str | None = None) -> str:
    """Status string for the Preferences UI (full paths; UI wraps/scrolls)."""
    text = (spec if spec is not None else load_conda_tools_env()).strip()
    if not text:
        if sys.platform == "win32":
            return (
                "Not set — using the default system PATH only. "
                "On Windows, set a WSL2 conda env name (e.g. tetrimmer) for RepeatMasker / FastGA."
            )
        return "Not set — using the default system PATH only."

    if sys.platform == "win32":
        if not wsl_available():
            native = resolve_conda_env_bin(text)
            if native is not None:
                return _describe_native_bin(native)
            return (
                "WSL2 not found. Install WSL2 with conda inside Ubuntu "
                "(or install a native Windows conda env), then set the env name."
            )
        wsl_bin = resolve_wsl_conda_bin(text)
        native = resolve_conda_env_bin(text)
        prefer_wsl = (
            text.startswith(("/", "~"))
            or _unc_wsl_path_to_linux(text) is not None
            or native is None
        )
        if wsl_bin and prefer_wsl:
            samples = []
            for name, label in (
                ("RepeatMasker", "RepeatMasker"),
                ("FastGA", "FastGA"),
                ("blastn", "BLASTn"),
            ):
                if _wsl_which(name, bin_dir=wsl_bin):
                    samples.append(label)
            sample_text = (
                "like " + ", ".join(samples) if samples else "no common tools detected in bin/"
            )
            return f"Using WSL2 {wsl_bin} ({sample_text})."
        if native is not None:
            return _describe_native_bin(native)
        return (
            f"Could not resolve “{text}” in WSL2. "
            "Use an env name (e.g. tetrimmer) or a Linux path like ~/miniconda3/envs/tetrimmer."
        )

    bin_dir = resolve_conda_env_bin(text)
    if bin_dir is None:
        return (
            f"Could not resolve “{text}”. "
            "Use an env name (e.g. tetrimmer) or a full path to the env folder."
        )
    return _describe_native_bin(bin_dir)


def _describe_native_bin(bin_dir: Path) -> str:
    samples = []
    for name, label in (
        ("RepeatMasker", "RepeatMasker"),
        ("FastGA", "FastGA"),
        ("blastn", "BLASTn"),
    ):
        if shutil.which(name, path=str(bin_dir)):
            samples.append(label)
    sample_text = (
        "like " + ", ".join(samples) if samples else "no common tools detected in bin/"
    )
    return f"Using {bin_dir} ({sample_text})."
