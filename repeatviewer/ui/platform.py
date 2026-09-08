"""OS integration helpers."""

from __future__ import annotations

import platform
import subprocess
from pathlib import Path
from repeatviewer.external_process import log_command


def reveal_in_file_manager(path: Path) -> None:
    """Show a file or folder in the system file manager."""
    path = path.expanduser().resolve()
    system = platform.system()
    try:
        if system == "Darwin":
            command = ["open", "-R", str(path)]
        elif system == "Windows":
            command = ["explorer", "/select,", str(path)]
        elif path.is_file():
            command = ["xdg-open", str(path.parent)]
        else:
            command = ["xdg-open", str(path)]
        log_command("file-manager", command, inputs=[path])
        subprocess.run(command, check=False)
    except OSError:
        pass
