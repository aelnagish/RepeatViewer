"""Contig view: lane packing, difference marks, consensus on the reference, and
loading a real result into the panel."""

from __future__ import annotations

import json

import numpy as np

from repeatviewer.backends import map_reference as mr
from repeatviewer.widgets import contig_view as cv


def test_pack_lanes():
    lanes, n = cv.pack_lanes([(0, 10), (5, 15), (12, 20), (30, 40)], gap=2)
    assert n == 2
    assert lanes[0] == lanes[2] == 0 and lanes[1] == 1 and lanes[3] == 0


def test_difference_marks():
    pos, bases, runs = cv.difference_marks("ACTT--GA", 10, "N" * 10 + "ACGTACGA")
    assert pos == [12] and bases == ["T"]
    assert runs == [(14, 16)]


def test_consensus_on_reference_marks_deletions():
    m = np.array([0, 1, -1, 2])
    assert cv.consensus_on_reference("ACT", m, 4) == "AC-T"


def test_feature_cell_text():
    assert cv.feature_cell_text(mr.STATUS_CHANGED, 3) == "3 diff"
    assert cv.feature_cell_text(mr.STATUS_INTACT, 0) == "intact"


def _payload():
    ref = "ACGTACGTACGTACGTACGT"
    reads = [
        {"name": "a", "start": 0, "end": 10, "strand": "+", "row": "ACGTTCGTAC", "ins": [],
         "identity": 0.9, "mismatches": 1, "deleted": 0, "inserted": 0, "mapq": -1,
         "kind": "primary", "qlen": 10, "clip_left": 0, "placement": 1,
         "lift": [[0, "changed", 1, 0.8, 1, 5]]},
        {"name": "b", "start": 5, "end": 20, "strand": "-", "row": "CGTAC--TACGTACGT",
         "ins": [[8, "G"]], "identity": 0.85, "mismatches": 0, "deleted": 2, "inserted": 1,
         "mapq": -1, "kind": "primary", "qlen": 14, "clip_left": 0, "placement": 1,
         "lift": [[0, "partial", 0, 1.0, 1, 1]]},
    ]
    reads[1]["row"] = reads[1]["row"][:15]
    return {
        "format": mr.RVMAP_FORMAT, "version": 1, "title": "x mapped to ref",
        "reference": "", "queries": [], "engine": "mafft", "engine_version": "7.5",
        "summary": {"n_sequences": 3, "n_placed": 2, "n_placements": 2, "n_reversed": 1,
                    "n_unused": 1, "mean_identity": 0.875},
        "warnings": [], "details_dir": "x_map_details", "details": {},
        "unused": [["c", "no alignment found"]],
        "contigs": [{
            "id": "ref", "length": 20, "detailed": True, "n_reads": 2, "mean_identity": 0.875,
            "reference": ref, "consensus": ref, "ref_to_cons": mr.rle(np.arange(20)),
            "coverage": mr.rle(np.array([1] * 5 + [2] * 5 + [1] * 10)), "bin_size": 1,
            "features": [{"start": 1, "end": 5, "name": "motif", "feature_type": "motif",
                          "color": "#e0413a", "strand": "+", "track": 0, "parent": ""}],
            "feature_summary": [{"name": "motif", "type": "motif", "start": 1, "end": 5,
                                 "intact": 0, "changed": 1, "partial": 1, "lost": 0,
                                 "mean_identity": 0.8}],
            "variants": [[5, "A", "T", 1, 2, "", "motif"]],
            "reads": reads,
        }],
    }


def test_panel_loads_result_and_fills_tables(qapp, tmp_path):
    path = tmp_path / "x.fa_map.rvmap"
    path.write_text(json.dumps(_payload()))
    panel = cv.ContigViewPanel(path)
    panel.resize(900, 600)
    panel.show()
    qapp.processEvents()
    model = panel.current_model()
    assert model is not None and model.length == 20 and len(model.reads) == 2
    assert model.reads[1].label == "b (reversed)"
    assert model.reads[0].mism_pos == [4]
    assert panel._reads_view.model().rowCount() == 2
    assert panel._features_view.model().columnCount() == 3          # name + motif + hidden #
    assert panel._variants_view.model().rowCount() == 1
    assert panel._unused_view.model().rowCount() == 1
    panel.canvas.show_base(10, zoom_to_letters=True)
    assert panel.canvas.letters()
    panel.canvas.show_read(1)
    assert panel.canvas.selected == 1
    panel.grab()                                                     # paints without error
    panel.canvas.set_one_per_row(False)
    panel.grab()


def test_panel_reports_a_bad_file(qapp, tmp_path):
    path = tmp_path / "bad.rvmap"
    path.write_text("{}")
    panel = cv.ContigViewPanel(path)
    assert panel.current_model() is None
    assert "Cannot open" in panel._title.text()
