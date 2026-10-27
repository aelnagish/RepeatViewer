"""Contig view: lane packing, difference marks, consensus on the reference, and
loading a real result into the panel."""

from __future__ import annotations

import json

import numpy as np

from repeatviewer.backends import map_reference as mr
from repeatviewer.backends.sequence_features import SequenceFeature
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


def test_bases_match_is_iupac_aware_and_case_blind():
    assert cv.bases_match("A", "A") and cv.bases_match("a", "A")
    assert cv.bases_match("Y", "C") and cv.bases_match("Y", "t")
    assert cv.bases_match("N", "G")
    assert not cv.bases_match("Y", "A")
    assert not cv.bases_match("A", "G")


def test_difference_marks_ignore_ambiguity_codes_in_reference():
    pos, bases, runs = cv.difference_marks("ACTA", 0, "ACYG")
    assert pos == [3] and bases == ["A"] and runs == []


def test_difference_profile_counts_mismatches_and_deletions():
    reads = [
        cv.ReadView({}, "a", 0, 6, "+", 1.0, [2], ["T"], [(4, 6)]),
        cv.ReadView({}, "b", 0, 6, "-", 1.0, [2, 3], ["T", "A"], []),
    ]
    prof = cv.difference_profile(reads, 8)
    assert prof.tolist() == [0, 0, 2, 1, 1, 1, 0, 0]


def test_step_difference_wraps_both_ways():
    pos = [10, 20, 30]
    assert cv.step_difference(pos, None, 1) == 10
    assert cv.step_difference(pos, None, -1) == 30
    assert cv.step_difference(pos, 10, 1) == 20
    assert cv.step_difference(pos, 30, 1) == 10
    assert cv.step_difference(pos, 10, -1) == 30
    assert cv.step_difference(pos, 25, -1) == 20
    assert cv.step_difference([], 5, 1) is None


def test_difference_positions_from_variants():
    variants = [[50, "C", "A", 1, 1, "CG", ""], [52, "G", "C", 1, 1, "", ""],
                [50, "C", "T", 1, 1, "CG", ""]]
    assert cv.difference_positions(variants) == [49, 51]


def test_read_passes_filters():
    d = {"mapq": 30, "kind": "primary"}
    assert cv.read_passes(d, 0.95)
    assert not cv.read_passes(d, 0.85, min_identity=0.9)
    assert not cv.read_passes(d, 0.99, min_mapq=40)
    assert cv.read_passes({"mapq": -1}, 0.99, min_mapq=40)       # no MAPQ reported: kept
    assert not cv.read_passes({"kind": "secondary"}, 1.0, hide_secondary=True)
    assert cv.read_passes({"kind": "secondary"}, 1.0)


def test_sort_order_keys():
    reads = [
        cv.ReadView({}, "b", 50, 90, "-", 0.97),
        cv.ReadView({}, "a", 10, 60, "+", 0.90),
        cv.ReadView({}, "c", 30, 80, "+", 1.00),
    ]
    assert cv.sort_order(reads, cv.SORT_START) == [1, 2, 0]
    assert cv.sort_order(reads, cv.SORT_STRAND) == [1, 2, 0]
    assert cv.sort_order(reads, cv.SORT_IDENTITY) == [2, 0, 1]
    assert cv.sort_order(reads, cv.SORT_NAME) == [1, 0, 2]


def test_read_color_hex_modes():
    assert cv.read_color_hex("+", 1.0, 60, cv.COLOR_STRAND) != cv.read_color_hex(
        "-", 1.0, 60, cv.COLOR_STRAND)
    assert cv.read_color_hex("+", 1.0, None, cv.COLOR_IDENTITY) == "#5fae66"
    assert cv.read_color_hex("+", 0.5, None, cv.COLOR_IDENTITY) == "#d0574a"
    assert cv.read_color_hex("+", 1.0, -1, cv.COLOR_MAPQ) == "#cfcfcf"


def test_feature_status_text():
    assert cv.feature_status_text({"intact": 7}) == "intact"
    assert cv.feature_status_text({"changed": 2}) == "changed"
    assert cv.feature_status_text({"intact": 3, "partial": 4}) == "mixed"
    assert cv.feature_status_text({}) == ""


def _rv(name, start, row, strand="+", lift=None):
    r = cv.ReadView({"row": row, "lift": lift or []}, name, start, start + len(row), strand, 1.0)
    return r


def test_count_matrix_and_strand_coverage():
    reads = [_rv("a", 0, "ACGT"), _rv("b", 1, "CG-T", "-"), _rv("c", 2, "GT")]
    reads[2].hidden = True                                   # filtered out: not counted
    counts = cv.count_matrix(reads, 6)
    assert counts[0].tolist() == [1, 0, 0, 0, 0, 0]
    assert counts[1].tolist() == [0, 2, 0, 0, 0, 0]
    assert counts[3].tolist() == [0, 0, 0, 1, 0, 1]          # one T, one gap
    fwd, rev = cv.strand_coverage(reads, 6)
    assert fwd.tolist() == [1, 1, 1, 1, 0, 0]
    assert rev.tolist() == [0, 1, 1, 1, 1, 0]


def test_consensus_from_counts_thresholds_and_no_coverage():
    counts = np.zeros((4, 6), dtype=np.int32)
    counts[0, 0] = 3                       # A A A
    counts[1, 1], counts[1, 3] = 2, 2      # C/T tie
    counts[2, 5] = 3                       # deleted
    assert cv.consensus_from_counts(counts, 0.0, "?") == "AY-?"
    assert cv.consensus_from_counts(counts, 0.0, "N")[-1] == "N"


def test_feature_type_counts_order_and_colour():
    feats = [SequenceFeature(start=1, end=5, name="R1", feature_type="repeat", color="#111111"),
             SequenceFeature(start=2, end=3, name="M", feature_type="motif", color="#222222"),
             SequenceFeature(start=6, end=9, name="R2", feature_type="repeat", color="#333333")]
    assert cv.feature_type_counts(feats) == [("repeat", 2, "#111111"), ("motif", 1, "#222222")]


def test_contig_statistics_rows():
    reads = [_rv("a", 0, "ACGT"), _rv("b", 0, "ACTT", "-")]
    reads[1].mism_pos, reads[1].mism_base = [2], ["T"]
    m = cv.ContigModel(ref_id="r", length=6, detailed=True, reference="ACGTAA",
                       consensus="ACGT", cons_on_ref="ACGTNN",
                       coverage=np.array([2, 2, 2, 2, 0, 0]), bin_size=1, features=[],
                       reads=reads, variants=[[3, "G", "T", 1, 2, "", ""]],
                       feature_summary=[{"intact": 1, "lost": 1}], n_reads=2)
    m.diff_count = cv.difference_profile(reads, 6)
    m.diff_pos = cv.difference_positions(m.variants)
    rows = dict(cv.contig_statistics(m))
    assert rows["Forward / reverse"] == "1 / 1"
    assert rows["Reference covered"].startswith("4 bp")
    assert rows["Identical sites"].startswith("3 of 4")
    assert rows["Variant positions"] == "1"
    assert "1 intact" in rows["Features in sequences"] and "1 lost" in rows["Features in sequences"]
