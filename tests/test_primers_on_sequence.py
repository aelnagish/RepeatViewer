"""Primers on a sequence: placement, tails, sidecar, Primers tab of the Annotation view."""

from __future__ import annotations

import os

import pytest

os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")

from repeatviewer.backends.primers import on_sequence as osq  # noqa: E402
from repeatviewer.backends.primers.oligos import Primer  # noqa: E402

pytest.importorskip("primer3")

TEMPLATE = ("GGCATCGATCGATTACGGA" * 3
            + "AAGTTCGATCCTGCAAGTCGTAGACCTTGAACGGTTCAAGGCTTACCGTAATGGCAACTTGAGC"
            + "TTAGCGGCATACGATTTGCAGCAAGT" * 4)
FWD_WITH_TAIL = "ATGGTACCCCTGCTAGTCGTAGACCTTGAA"   # KpnI tail + 1 mismatch
REV = "ACTTGCTGCAAATCGTATGCC"


def test_guess_tail_and_place():
    seq = osq.guess_tail(FWD_WITH_TAIL, TEMPLATE)
    assert seq.startswith("atggtacc") and seq[8:] == FWD_WITH_TAIL[8:]
    place = osq.place_primer(seq, TEMPLATE)
    assert place.strand == "+" and place.tail == 8 and place.annealed == 22
    assert place.mismatches == 1
    assert TEMPLATE[place.mismatch_positions[0] - 1] == "A"
    rev = osq.place_primer(REV, TEMPLATE)
    assert rev.strand == "-" and rev.mismatches == 0


def test_fix_resize_snap():
    seq = osq.guess_tail(FWD_WITH_TAIL, TEMPLATE)
    place = osq.place_primer(seq, TEMPLATE)
    fixed = osq.fix_mismatches(seq, place)
    assert osq.place_primer(fixed, TEMPLATE).mismatches == 0
    longer = osq.resize_on_template(fixed, place, TEMPLATE, 2, 1)
    assert len(longer) == len(fixed) + 3
    assert osq.place_primer(longer, TEMPLATE).mismatches == 0
    snapped = osq.snap_to_tm(fixed, place, TEMPLATE, 55.0)
    assert snapped.endswith(fixed[-5:])


def test_tails_and_enzymes():
    assert osq.restriction_tail("KpnI", "AT") == "atggtacc"
    assert osq.with_tail("ACGTACGTACGTACGTAC", "gaattc").startswith("gaattc")
    assert "KpnI" in osq.enzymes_not_cutting(TEMPLATE)
    assert osq.enzyme_label("EcoRI").startswith("EcoRI")


def test_sidecar_and_features(tmp_path):
    fasta = tmp_path / "t.fa"
    fasta.write_text(">t\n" + TEMPLATE + "\n")
    prim = [Primer("F", osq.guess_tail(FWD_WITH_TAIL, TEMPLATE)), Primer("R", REV)]
    path = osq.save_sequence_primers(fasta, prim)
    assert path.name == "t.primers_RV.json"
    back = osq.load_sequence_primers(fasta)
    assert [p.sequence for p in back] == [p.sequence for p in prim]
    places = {p.name: osq.place_primer(p.sequence, TEMPLATE) for p in back}
    feats = osq.merge_primer_features([], back, places)
    assert len(feats) == 2 and {osq.feature_primer_name(f) for f in feats} == {"F", "R"}
    again = osq.merge_primer_features(feats, back[:1], places)
    assert [osq.feature_primer_name(f) for f in again] == ["F"]
    legacy = osq.primers_from_features(feats, TEMPLATE, [])
    assert {p.name for p in legacy} == {"F", "R"}
    osq.save_sequence_primers(fasta, [])
    assert not path.exists()
    from repeatviewer.backends.sequence_features import is_sequence_map_hidden_file
    assert is_sequence_map_hidden_file(path)


def test_primers_tab_in_annotation_view(tmp_path):
    from PySide6.QtWidgets import QApplication

    app = QApplication.instance() or QApplication([])
    from repeatviewer.widgets.sequence_map_viewer import SequenceMapPanel

    fasta = tmp_path / "t.fa"
    fasta.write_text(">t\n" + TEMPLATE + "\n")
    panel = SequenceMapPanel(fasta, autoload=True)
    panel.show_primers_tab()
    pp = panel.primers_panel()
    pp.add_primers([Primer("F", FWD_WITH_TAIL), Primer("R", REV)], tick=True)
    app.processEvents()
    assert osq.primers_sidecar_path(fasta).is_file()
    names = {osq.feature_primer_name(f) for f in panel.primers_host_features()
             if osq.is_primer_feature(f)}
    assert names == {"F", "R"}
    assert pp.placement("F").tail == 8
    assert len(pp.ticked_primers()) == 2
    pp.remove_primers(["R"])
    assert [p.name for p in pp.primers()] == ["F"]
    from repeatviewer.views.primer_edit_dialog import EditPrimerDialog

    dlg = EditPrimerDialog(pp.primers()[0], template=TEMPLATE)
    dlg._update()
    assert dlg._sites and "Binding site" in dlg.metrics.text()
    dlg._on_apply_tail()
    panel.deleteLater()
