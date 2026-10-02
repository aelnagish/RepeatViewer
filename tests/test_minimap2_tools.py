"""Bundled minimap2 2.31 lookup."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends import map_reference as mr
from repeatviewer.backends import minimap2_tools as m2
from repeatviewer.job_options import MAP_ENGINE_MINIMAP2


def test_bundled_minimap2_exists_for_each_platform():
    assert m2.bundled_minimap2_executable("Darwin", machine="arm64") is not None
    assert m2.bundled_minimap2_executable("Darwin", machine="x86_64") is not None
    assert Path(m2.bundled_minimap2_executable("Linux", machine="x86_64")).name == "minimap2"
    assert m2.bundled_minimap2_executable("Linux", machine="aarch64") is None
    win = m2.bundled_minimap2_executable("Windows", machine="x86_64")
    assert win is not None and win.lower().endswith("minimap2.exe")


def test_engine_prefers_bundled_minimap2(monkeypatch):
    bundled = m2.bundled_minimap2_executable()
    assert bundled
    monkeypatch.setattr(mr.shutil, "which", lambda name: "/usr/bin/fake-minimap2")
    assert mr.engine_executable(MAP_ENGINE_MINIMAP2) == bundled
    assert mr.engine_executable(MAP_ENGINE_MINIMAP2, system_path_only=True) == "/usr/bin/fake-minimap2"


def test_install_hint_mentions_bundle():
    hint = mr.engine_install_hint(MAP_ENGINE_MINIMAP2)
    assert "bundled" in hint.lower() and m2.MINIMAP2_VERSION in hint
