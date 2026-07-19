"""OS integration helpers."""

from __future__ import annotations

import platform
import subprocess
from pathlib import Path


def reveal_in_file_manager(path: Path) -> None:
    """Show a file or folder in the system file manager."""
    path = path.expanduser().resolve()
    system = platform.system()
    try:
        if system == "Darwin":
            subprocess.run(["open", "-R", str(path)], check=False)
        elif system == "Windows":
            subprocess.run(["explorer", "/select,", str(path)], check=False)
        elif path.is_file():
            subprocess.run(["xdg-open", str(path.parent)], check=False)
        else:
            subprocess.run(["xdg-open", str(path)], check=False)
    except OSError:
        pass
