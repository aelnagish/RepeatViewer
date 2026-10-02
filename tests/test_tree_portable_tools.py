"""Bundled IQ-TREE lookup and PATH override."""

from __future__ import annotations

import os
import stat
from pathlib import Path

import pytest

from repeatviewer.backends import tree_builder as tb
from repeatviewer.backends import tree_portable_tools as tpt
from repeatviewer.job_options import TREE_ENGINE_FASTME, TREE_ENGINE_IQTREE, TreeOptions


def test_bundled_iqtree_exists_for_each_platform():
    assert tpt.bundled_iqtree_executable("Darwin") is not None
    assert Path(tpt.bundled_iqtree_executable("Darwin")).name == "iqtree3"
    assert tpt.bundled_iqtree_executable("Linux") is not None
    win = tpt.bundled_iqtree_executable("Windows", machine="x86_64")
    assert win is not None and win.lower().endswith("iqtree3.exe")
    assert tpt.bundled_iqtree_executable("Windows", machine="aarch64") is None


def test_bundled_preferred_over_path(monkeypatch):
    bundled = tpt.bundled_engine_executable(TREE_ENGINE_IQTREE)
    assert bundled
    monkeypatch.setattr(tb.shutil, "which", lambda name: "/usr/bin/fake-iqtree")
    assert tb.engine_executable(TREE_ENGINE_IQTREE) == bundled
    assert tb.engine_executable(TREE_ENGINE_IQTREE, system_path_only=True) == "/usr/bin/fake-iqtree"


def test_fastme_and_fasttree_are_not_bundled():
    assert tpt.bundled_engine_executable(TREE_ENGINE_FASTME) is None
    assert tpt.bundled_engine_executable("fasttree") is None
    assert tpt.missing_portable_engines() == []


def test_engine_uses_system_path_per_engine():
    o = TreeOptions(use_system_fastme=True)
    assert tb.engine_uses_system_path(o, TREE_ENGINE_FASTME)
    assert not tb.engine_uses_system_path(o, TREE_ENGINE_IQTREE)
    assert tb.engine_uses_system_path(TreeOptions(use_system_path=True), TREE_ENGINE_IQTREE)


def test_download_helper_removed():
    with pytest.raises(tpt.PortableToolError):
        tpt.download_portable_engine("iqtree")


@pytest.mark.skipif(os.name == "nt", reason="unix executable bit")
def test_legacy_portable_alias_points_at_bundle(tmp_path, monkeypatch):
    # Old name still resolves to the vendored IQ-TREE.
    assert tpt.portable_engine_executable(TREE_ENGINE_IQTREE) == tpt.bundled_iqtree_executable()
    assert "iqtree" in str(tpt.portable_tree_tools_dir())
