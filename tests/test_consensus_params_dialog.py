"""Consensus parameters dialog: sequence type and ambiguous character."""

from __future__ import annotations

from PySide6.QtWidgets import QApplication, QLabel

from repeatviewer.job_options import (
    ConsensusOptions,
    MAFFT_SEQUENCE_PROTEIN,
)
from repeatviewer.views.tool_dialogs import ConsensusParamsDialog


def test_consensus_dialog_sequence_type_and_help(qapp: QApplication) -> None:
    del qapp
    dlg = ConsensusParamsDialog(initial=ConsensusOptions())
    dlg.show()
    assert any(lab.text() == "Sequence type:" for lab in dlg.findChildren(QLabel))
    assert dlg._sequence_type.currentText() == "Auto-detect"
    assert dlg._ambiguous_label.text() == "Ambiguous base:"


def test_consensus_protein_switches_ambiguous_label(qapp: QApplication) -> None:
    del qapp
    dlg = ConsensusParamsDialog(initial=ConsensusOptions())
    dlg.show()
    idx = dlg._sequence_type.findData(MAFFT_SEQUENCE_PROTEIN)
    dlg._sequence_type.setCurrentIndex(idx)
    assert dlg._ambiguous_label.text() == "Ambiguous residue:"
    assert dlg._ambiguous.currentText() == "X"
    opts = dlg.options()
    assert opts.sequence_type == MAFFT_SEQUENCE_PROTEIN
    assert opts.ambiguous == "X"
