"""Shared pytest fixtures for RepeatViewer tests."""

from __future__ import annotations

import os

import pytest

# Headless-friendly Qt for dialog / widget tests in CI and sandboxes.
os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")


@pytest.fixture(scope="session")
def qapp():
    """Session-scoped QApplication for widget tests."""
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance()
    if app is None:
        app = QApplication([])
    return app
