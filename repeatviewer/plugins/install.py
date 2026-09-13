"""Install / uninstall third-party RepeatViewer plugins via pip."""

from __future__ import annotations

import importlib
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from repeatviewer.plugins.loader import ENTRY_POINT_GROUP, _iter_entry_points


@dataclass(frozen=True)
class InstalledEntryPlugin:
    """One package that registers under ``repeatviewer.plugins``."""

    entry_name: str
    distribution: str
    module: str
    label: str = ""


def list_installed_entry_plugins() -> list[InstalledEntryPlugin]:
    """Return entry-point plugins currently installed in this Python env."""
    out: list[InstalledEntryPlugin] = []
    seen: set[tuple[str, str]] = set()
    for ep in _iter_entry_points():
        dist_name = ""
        dist = getattr(ep, "dist", None)
        if dist is not None:
            dist_name = str(getattr(dist, "name", "") or "").strip()
            if not dist_name:
                try:
                    dist_name = str(dist.metadata["Name"]).strip()
                except Exception:  # noqa: BLE001
                    dist_name = ""
        if not dist_name:
            dist_name = _distribution_for_entry(str(ep.name)) or ""
        if not dist_name:
            continue
        key = (str(ep.name), dist_name)
        if key in seen:
            continue
        seen.add(key)
        # Prefer a human label from the plugin class when import works.
        label = str(ep.name)
        try:
            loaded = ep.load()
            plugin = loaded() if isinstance(loaded, type) else loaded
            plugin_label = str(getattr(plugin, "label", "") or "").strip()
            if plugin_label:
                label = plugin_label
        except Exception:  # noqa: BLE001
            pass
        out.append(
            InstalledEntryPlugin(
                entry_name=str(ep.name),
                distribution=dist_name,
                module=str(getattr(ep, "value", "") or ""),
                label=label,
            )
        )
    out.sort(key=lambda item: item.label.lower())
    return out


def _distribution_for_entry(entry_name: str) -> str | None:
    try:
        from importlib.metadata import distributions
    except ImportError:  # pragma: no cover
        return None
    for dist in distributions():
        try:
            eps = dist.entry_points
            try:
                selected = list(eps.select(group=ENTRY_POINT_GROUP))
            except Exception:  # noqa: BLE001
                selected = [ep for ep in eps if getattr(ep, "group", "") == ENTRY_POINT_GROUP]
            for ep in selected:
                if ep.name == entry_name:
                    return str(dist.metadata["Name"])
        except Exception:  # noqa: BLE001
            continue
    return None


def run_pip(args: list[str], *, timeout_s: float = 600.0) -> tuple[int, str]:
    """Run ``python -m pip …``; return (returncode, combined output)."""
    cmd = [sys.executable, "-m", "pip", *args]
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout_s,
            check=False,
        )
    except subprocess.TimeoutExpired as exc:
        out = f"{exc.stdout or ''}{exc.stderr or ''}".strip()
        return 1, (out + "\nTimed out.").strip()
    except OSError as exc:
        return 1, str(exc)
    combined = "\n".join(
        part for part in (proc.stdout or "", proc.stderr or "") if part.strip()
    ).strip()
    return int(proc.returncode), combined


def install_plugin_path(path: Path) -> tuple[bool, str]:
    """Install a plugin package from a folder (editable) or a wheel/sdist file."""
    target = path.expanduser().resolve()
    if target.is_dir():
        if not (
            (target / "pyproject.toml").is_file()
            or (target / "setup.py").is_file()
            or (target / "setup.cfg").is_file()
        ):
            return (
                False,
                "Selected folder is not a Python package "
                "(need pyproject.toml or setup.py).",
            )
        code, output = run_pip(["install", "-e", str(target)])
    elif target.is_file() and target.suffix.lower() in {".whl", ".zip"}:
        code, output = run_pip(["install", str(target)])
    elif target.is_file() and target.name.endswith(".tar.gz"):
        code, output = run_pip(["install", str(target)])
    else:
        return False, "Choose a plugin folder or a .whl / source archive."
    if code != 0:
        return False, output or f"pip failed with exit code {code}."
    _refresh_import_caches()
    return True, output or "Installed."


def uninstall_plugin_distribution(distribution: str) -> tuple[bool, str]:
    """Uninstall a distribution that provides RepeatViewer plugins."""
    name = (distribution or "").strip()
    if not name:
        return False, "No package name."
    # Never uninstall the host app itself.
    if name.lower().replace("_", "-") in {"repeatviewer", "repeeatviewer"}:
        return False, f"Refusing to uninstall the host package “{name}”."
    code, output = run_pip(["uninstall", "-y", name])
    if code != 0:
        return False, output or f"pip failed with exit code {code}."
    _refresh_import_caches()
    return True, output or "Uninstalled."


def _refresh_import_caches() -> None:
    importlib.invalidate_caches()
    try:
        from importlib.metadata import entry_points

        entry_points()
    except Exception:  # noqa: BLE001
        pass
