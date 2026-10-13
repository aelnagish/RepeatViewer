"""Loaded plugin registry."""

from __future__ import annotations

import logging
from pathlib import Path

from repeatviewer.plugins.api import ProgressCallback, RepeatViewerPlugin
from repeatviewer.project import RepeatViewerProject

logger = logging.getLogger(__name__)

_FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".faa"}


class PluginRegistry:
    def __init__(self) -> None:
        self._plugins: dict[str, RepeatViewerPlugin] = {}

    def register(self, plugin: RepeatViewerPlugin) -> None:
        plugin_id = plugin.id.strip()
        if not plugin_id or not plugin_id.replace("_", "").isalnum():
            raise ValueError(f"Invalid plugin id: {plugin.id!r}")
        if plugin_id in self._plugins:
            raise ValueError(f"Duplicate plugin id: {plugin_id}")
        self._plugins[plugin_id] = plugin

    def register_or_replace(self, plugin: RepeatViewerPlugin) -> None:
        plugin_id = plugin.id.strip()
        if not plugin_id or not plugin_id.replace("_", "").isalnum():
            raise ValueError(f"Invalid plugin id: {plugin.id!r}")
        self._plugins[plugin_id] = plugin

    def unregister(self, plugin_id: str) -> bool:
        return self._plugins.pop(plugin_id, None) is not None

    def all_plugins(self) -> list[RepeatViewerPlugin]:
        return sorted(self._plugins.values(), key=lambda p: p.label.lower())

    def get(self, plugin_id: str) -> RepeatViewerPlugin | None:
        return self._plugins.get(plugin_id)

    def filter_files(self, plugin: RepeatViewerPlugin, paths: list[Path]) -> list[Path]:
        if plugin.file_suffixes is not None:
            allowed = {s.lower() for s in plugin.file_suffixes}
            return [p for p in paths if p.is_file() and p.suffix.lower() in allowed]
        if plugin.fasta_only:
            return [p for p in paths if p.suffix.lower() in _FASTA_SUFFIXES]
        return [p for p in paths if p.is_file()]

    def run_plugin(
        self,
        plugin_id: str,
        seq_files: list[Path],
        project: RepeatViewerProject,
        progress_cb: ProgressCallback,
        *,
        max_workers: int = 5,
        cancel=None,
    ) -> dict:
        plugin = self._plugins.get(plugin_id)
        if plugin is None:
            raise KeyError(f"Unknown plugin: {plugin_id}")
        return plugin.run_batch(
            seq_files,
            project,
            progress_cb,
            max_workers=max_workers,
            cancel=cancel,
        ).as_dict()


_registry: PluginRegistry | None = None


def get_registry() -> PluginRegistry:
    global _registry
    if _registry is None:
        _registry = PluginRegistry()
    return _registry


def reset_registry() -> None:
    global _registry
    _registry = PluginRegistry()
