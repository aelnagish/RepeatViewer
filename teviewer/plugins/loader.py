"""Discover plugins via Python entry points (``teviewer.plugins`` group)."""

from __future__ import annotations

import logging
import sys
from importlib.metadata import entry_points

from teviewer.plugins.custom_command import (
    UserCommandPlugin,
    load_user_command_plugin_specs,
)
from teviewer.plugins.api import TEviewerPlugin
from teviewer.plugins.registry import PluginRegistry, get_registry, reset_registry

logger = logging.getLogger(__name__)

ENTRY_POINT_GROUP = "teviewer.plugins"


def _iter_entry_points():
    if sys.version_info >= (3, 10):
        eps = entry_points()
        group = eps.select(group=ENTRY_POINT_GROUP)
        if group is not None:
            yield from group
        return
    for ep in entry_points().get(ENTRY_POINT_GROUP, []):
        yield ep


def load_plugins(registry: PluginRegistry | None = None) -> PluginRegistry:
    """
    Load all plugins registered under ``[project.entry-points."teviewer.plugins"]``.

    Third-party packages install their plugin class, for example::

        [project.entry-points."teviewer.plugins"]
        my_tool = "my_package.teviewer_plugin:MyPlugin"
    """
    reg = registry or get_registry()
    errors: list[str] = []

    for ep in _iter_entry_points():
        try:
            loaded = ep.load()
            plugin = loaded() if isinstance(loaded, type) else loaded
            if not isinstance(plugin, TEviewerPlugin):
                raise TypeError(
                    f"{ep.name} must return a TEviewerPlugin instance, got {type(plugin)}"
                )
            reg.register(plugin)
            logger.info("Loaded TEviewer plugin: %s (%s)", plugin.label, plugin.id)
        except Exception as exc:  # noqa: BLE001
            msg = f"{ep.name}: {exc}"
            errors.append(msg)
            logger.exception("Failed to load TEviewer plugin %s", ep.name)

    for spec in load_user_command_plugin_specs():
        try:
            reg.register(UserCommandPlugin(spec))
            logger.info("Loaded user command plugin: %s (%s)", spec.label, spec.id)
        except Exception as exc:  # noqa: BLE001
            msg = f"{spec.id}: {exc}"
            errors.append(msg)
            logger.exception("Failed to load user command plugin %s", spec.id)

    if errors:
        logger.warning("Plugin load errors:\n%s", "\n".join(errors))
    return reg


def reload_plugins() -> PluginRegistry:
    reset_registry()
    return load_plugins()
