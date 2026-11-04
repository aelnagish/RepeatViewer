"""Cytogenetics workspace: marker file, summary, tab names, and (with Qt) the page."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import karyotype as ky
from repeatviewer.backends.cytogenetics import model as md
from repeatviewer.backends.cytogenetics import workspace as ws


def test_marker_roundtrip(tmp_path):
    marker = ws.ensure_workspace(tmp_path)
    assert marker == ws.workspace_path(tmp_path) and marker.is_file()
    assert marker.parent.name == ci.STORE_DIR_NAME and marker.name == "workspace.rvcyto"
    assert ws.is_workspace_path(marker) and not ws.is_workspace_path(tmp_path / "a.fa")
    assert ws.workspace_folder(marker) == tmp_path.resolve()
    before = marker.stat().st_mtime_ns
    assert ws.ensure_workspace(tmp_path) == marker and marker.stat().st_mtime_ns == before


def test_summary(tmp_path):
    assert ws.folder_summary(tmp_path) == ws.WorkspaceSummary()
    a = ci.CytoImage(id="a", file="a.png")
    b = ci.CytoImage(id="b", file="b.png")
    md.set_image_details(a, md.ImageDetails(channels=[md.Channel(id="c", file="c.tif")],
                                            signals=[md.Signal(id="s", probe="P", x=0, y=0)]))
    ci.save_images(tmp_path, [a, b])
    md.add_slide(tmp_path, "Slide 1")
    ky.save_karyotype(tmp_path, ky.new_karyotype("k", 8, 2))
    s = ws.folder_summary(tmp_path)
    assert (s.images, s.slides, s.multichannel, s.with_signals, s.karyotypes) == (2, 1, 1, 1, 1)
    assert ws.summary_lines(s)[0] == ("Images", "2")


needs_qt = pytest.mark.skipif(importlib.util.find_spec("PySide6") is None, reason="PySide6 not installed")


def test_normalize_tab():
    assert ws.normalize_tab("signals") == "Signals"
    assert ws.normalize_tab("nope") == "Library" and ws.normalize_tab("") == "Library"


@needs_qt
def test_preview_routes_marker():
    from repeatviewer.widgets import preview_tabs as pt

    marker = Path("/x/.repeatviewer_cytogenetics/workspace.rvcyto")
    assert pt._preview_kind(marker) == "cytogenetics"
    assert pt.is_previewable_file(marker)
    assert pt.preview_path_priority(marker) == 2


@needs_qt
def test_workspace_widget(qapp, tmp_path):
    from repeatviewer.widgets.cytogenetics_workspace import CytogeneticsWorkspace

    ci.save_images(tmp_path, [ci.CytoImage(id="a", file="a.png", targets=[ci.ImageTarget("CroSat2.fasta")])])
    w = CytogeneticsWorkspace(ws.ensure_workspace(tmp_path))
    assert w.folder() == tmp_path.resolve() and w.tab_name() == "Library"
    w.show_tab("figure")
    assert w.tab_name() == "Figure"
    w.set_sequence_filter("CroSat2.fasta")
    assert w.sequence_filter() == "CroSat2.fasta"
    ci.save_images(tmp_path, [])
    from repeatviewer.widgets.cyto_images_ui import notify_changed

    notify_changed(tmp_path)
