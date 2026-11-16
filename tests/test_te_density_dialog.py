"""Genome density plot dialog: tab height follows content."""

from __future__ import annotations

from PySide6.QtWidgets import QApplication

from repeatviewer.job_options import (
    TE_GENOME_MAP_DIALOG_TAB_METHYLATION,
    TeGenomeMapOptions,
)
from repeatviewer.views.tool_dialogs import TeDensityPlotDialog


def test_te_density_methylation_height_follows_tab(qapp: QApplication) -> None:
    """Methylation must shrink vs TE — QTabWidget sizeHint alone stays on TE."""
    dlg = TeDensityPlotDialog(
        initial=TeGenomeMapOptions(
            last_dialog_tab=TE_GENOME_MAP_DIALOG_TAB_METHYLATION
        )
    )
    dlg.prepare_for_show()
    dlg.show()
    qapp.processEvents()
    meth_h = dlg.height()
    dlg._tabs.setCurrentIndex(0)  # TE
    qapp.processEvents()
    te_h = dlg.height()
    dlg._tabs.setCurrentIndex(2)  # Methylation
    qapp.processEvents()
    back_meth_h = dlg.height()
    assert meth_h < te_h
    assert back_meth_h < te_h
    assert back_meth_h <= meth_h + 20
    assert back_meth_h < 500
