"""Open FASTA in AliView or system viewer; MAFFT alignment (TEtrimmerGUI behaviour)."""

from __future__ import annotations

import logging
import os
import platform
import shutil
import subprocess
import sys
from pathlib import Path

from repeatviewer.backends.gui_functions import get_original_file_path
from repeatviewer.backends.sequence_aligner import run_alignment
from repeatviewer.external_process import log_command

OS_TYPE = platform.system()

_ALIVIEW_MIN_MEM = "512M"
_ALIVIEW_MAX_MEM = "1024M"


def aliview_paths() -> tuple[Path | None, Path | None]:
    """
    Return (launcher_script, jar_path) for bundled AliView under repeatviewer/backends/aliview/.
    """
    base = Path(get_original_file_path()) / "aliview"
    launcher = base / "aliview"
    jar_path = base / "aliview.jar"
    if OS_TYPE == "Windows":
        if jar_path.is_file():
            return None, jar_path
        return None, None
    if launcher.is_file():
        return launcher, jar_path if jar_path.is_file() else None
    if jar_path.is_file():
        return None, jar_path
    return None, None


def _java_runtime_works(java_bin: Path) -> bool:
    """True if ``java_bin`` is a real JVM (not the macOS /usr/bin/java stub)."""
    try:
        completed = subprocess.run(
            [str(java_bin), "-version"],
            capture_output=True,
            timeout=8,
            check=False,
        )
    except (OSError, subprocess.TimeoutExpired):
        return False
    text = (
        (completed.stderr or b"") + (completed.stdout or b"")
    ).decode("utf-8", errors="replace")
    if "Unable to locate a Java Runtime" in text:
        return False
    return "version" in text.lower()


def find_java_executable() -> str | None:
    """Locate a usable ``java`` binary for launching AliView."""
    candidates: list[Path] = []

    java_home = (os.environ.get("JAVA_HOME") or "").strip()
    if java_home:
        candidates.append(Path(java_home) / "bin" / "java")

    if OS_TYPE == "Darwin":
        # Prefer full JDKs under JavaVirtualMachines over the browser
        # JavaAppletPlugin that ``java_home`` sometimes returns.
        jvm_root = Path("/Library/Java/JavaVirtualMachines")
        if jvm_root.is_dir():
            candidates.extend(sorted(jvm_root.glob("*/Contents/Home/bin/java")))
        for brew_java in (
            Path("/opt/homebrew/opt/openjdk/bin/java"),
            Path("/usr/local/opt/openjdk/bin/java"),
        ):
            candidates.append(brew_java)
        try:
            home = subprocess.check_output(
                ["/usr/libexec/java_home"],
                stderr=subprocess.DEVNULL,
                text=True,
                timeout=5,
            ).strip()
            if home and "JavaAppletPlugin" not in home:
                candidates.append(Path(home) / "bin" / "java")
        except (OSError, subprocess.CalledProcessError, subprocess.TimeoutExpired):
            pass

    which = shutil.which("java")
    if which:
        candidates.append(Path(which))

    seen: set[str] = set()
    for candidate in candidates:
        key = str(candidate)
        if key in seen:
            continue
        seen.add(key)
        if "JavaAppletPlugin" in key:
            continue
        if candidate.is_file() and _java_runtime_works(candidate):
            return str(candidate)
    return None


def open_fasta(path: Path, *, use_system_viewer: bool = False) -> None:
    """Open a FASTA/MSA file in AliView or the OS default application."""
    path = path.expanduser().resolve()
    if not path.is_file():
        raise FileNotFoundError(path)
    suffix = path.suffix.lower()
    if suffix not in {".fa", ".fasta", ".fas", ".fna", ".faa"}:
        raise ValueError(f"Not a FASTA file: {path}")

    if use_system_viewer:
        _open_with_system_aliview(path)
        return

    _open_with_bundled_aliview(path)


def _open_with_bundled_aliview(path: Path) -> None:
    """Launch the JAR shipped with RepeatViewer via a real Java runtime."""
    _launcher, jar_path = aliview_paths()
    if jar_path is None or not jar_path.is_file():
        raise FileNotFoundError(
            "Bundled AliView JAR not found in repeatviewer/backends/aliview/. "
            "Reinstall RepeatViewer, or enable Settings → Use system sequence viewer "
            "after installing AliView.app."
        )

    java = find_java_executable()
    if java is None:
        raise RuntimeError(
            "AliView needs Java, but no working Java runtime was found.\n\n"
            "Install a JDK (e.g. Temurin/OpenJDK), or set JAVA_HOME to your JDK, "
            "then try again.\n\n"
            "On this Mac, /usr/bin/java is only a stub and cannot run AliView."
        )

    command = [
        java,
        f"-Xmx{_ALIVIEW_MAX_MEM}",
        f"-Xms{_ALIVIEW_MIN_MEM}",
        "-jar",
        str(jar_path),
        str(path),
    ]
    log_command("AliView", command, inputs=[path])
    subprocess.Popen(
        command,
        stdout=subprocess.DEVNULL,
        stderr=subprocess.DEVNULL,
        start_new_session=True,
    )


def _mac_open_aliview_app(path: Path) -> bool:
    """Try ``open -a AliView`` (installed .app). Returns True on success."""
    for app_name in ("AliView", "aliview"):
        completed = subprocess.run(
            ["open", "-a", app_name, "--", str(path)],
            capture_output=True,
            check=False,
        )
        if completed.returncode == 0:
            log_command(
                "system-AliView",
                ["open", "-a", app_name, "--", str(path)],
                inputs=[path],
            )
            return True
    return False


def _open_with_system_aliview(path: Path) -> None:
    """Open with a user-installed AliView (not the generic file association)."""
    if OS_TYPE == "Darwin":
        if _mac_open_aliview_app(path):
            return
        # Fall through to PATH / clear error — never bare ``open file.fa``
        # (macOS has no default handler for .fa → kLSApplicationNotFoundErr).

    launcher = shutil.which("aliview")
    if launcher:
        command = [launcher, str(path)]
        log_command("system-AliView", command, inputs=[path])
        subprocess.Popen(
            command,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            start_new_session=True,
        )
        return

    if OS_TYPE == "Windows":
        # Last resort: let Windows file association try.
        os.startfile(str(path))  # type: ignore[attr-defined]
        return

    raise RuntimeError(
        "System AliView was not found.\n\n"
        "Install AliView.app (macOS) or put ``aliview`` on your PATH, "
        "or uncheck Settings → Use system sequence viewer to use the "
        "AliView bundled with RepeatViewer (requires Java)."
    )


def _open_with_system(path: Path) -> None:
    """Legacy helper: prefer AliView; avoid unclaimed .fa file associations."""
    _open_with_system_aliview(path)


def run_mafft_alignment(input_fasta: Path, num_threads: int = 5) -> Path:
    """Run MAFFT (legacy alias). Prefer ``sequence_aligner.run_alignment``."""
    return run_alignment(input_fasta, num_threads=num_threads)
