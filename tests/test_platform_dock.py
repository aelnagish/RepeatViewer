"""macOS Dock helper for frozen Flexidot re-entries."""

from repeatviewer.ui.platform import hide_from_macos_dock


def test_hide_from_macos_dock_is_safe() -> None:
    # Must not raise on any platform (no-op off macOS; best-effort on macOS).
    hide_from_macos_dock()
