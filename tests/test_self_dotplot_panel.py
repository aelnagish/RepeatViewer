"""Self dotplot panel refresh scheduling."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pytest
from PySide6.QtWidgets import QApplication

from teviewer.widgets.self_dotplot_panel import SelfDotplotPanel

_app = QApplication.instance() or QApplication([])


@pytest.fixture
def fasta_file(tmp_path: Path) -> Path:
    path = tmp_path / "reads.fasta"
    path.write_text(">seq1\nACGTACGT\n")
    return path


@pytest.fixture
def panel(fasta_file: Path) -> SelfDotplotPanel:
    widget = SelfDotplotPanel(fasta_file)
    yield widget
    widget.close()
    widget.deleteLater()


def test_output_paths_are_unique_per_token(panel: SelfDotplotPanel) -> None:
    p1 = panel._output_path_for_token(1)
    p2 = panel._output_path_for_token(2)
    assert p1 != p2
    assert p1.name.endswith("_1.png")
    assert p2.name.endswith("_2.png")


def test_refresh_while_busy_queues_rerun(panel: SelfDotplotPanel) -> None:
    panel._worker_busy = True
    panel.refresh()
    assert panel._rerun_pending is True


def test_begin_render_clears_rerun_pending(panel: SelfDotplotPanel) -> None:
    panel._rerun_pending = True
    with patch.object(panel, "_start_worker") as start_worker:
        panel._begin_render()
    assert panel._rerun_pending is False
    start_worker.assert_called_once()


def test_schedule_refresh_starts_debounce_timer(panel: SelfDotplotPanel) -> None:
    panel._debounce.stop()
    panel._schedule_refresh()
    assert panel._debounce.isActive()
