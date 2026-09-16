"""BLAST cancel wiring and process-monitor bridge (main-thread UI updates)."""

from __future__ import annotations

import inspect

from repeatviewer.backends import gui_functions
from repeatviewer.backends.curation import run_on_sequence
from repeatviewer.widgets.process_monitor_dialog import (
    ProcessMonitorBridge,
    ProcessMonitorDialog,
)


def test_blast_accepts_cancel_parameter() -> None:
    sig = inspect.signature(gui_functions.blast)
    assert "cancel" in sig.parameters


def test_check_database_accepts_cancel_parameter() -> None:
    sig = inspect.signature(gui_functions.check_database)
    assert "cancel" in sig.parameters


def test_run_on_sequence_passes_cancel_into_signature() -> None:
    sig = inspect.signature(run_on_sequence)
    assert "cancel" in sig.parameters


def test_process_monitor_bridge_buffers_lines_without_live_window(qapp) -> None:
    dialog = ProcessMonitorDialog("Align")
    buffer: list[str] = []
    bridge = ProcessMonitorBridge(
        dialog,
        buffer,
        live=False,
        title="Align",
        cancelled_message="__CANCELLED__",
    )
    bridge.on_line("cmd blastn …")
    bridge.on_progress(1, "step two")
    assert buffer == ["cmd blastn …", "step two"]
    # Failure with live=False opens the captured log.
    bridge.on_failed("boom")
    assert dialog.isVisible()
    dialog.close()
    bridge.deleteLater()
    dialog.deleteLater()


def test_process_monitor_bridge_cancelled_is_quiet(qapp) -> None:
    dialog = ProcessMonitorDialog("BLAST")
    buffer: list[str] = []
    bridge = ProcessMonitorBridge(
        dialog,
        buffer,
        live=True,
        title="BLAST",
        cancelled_message="__CANCELLED__",
    )
    dialog.prepare("BLAST")
    bridge.on_failed("__CANCELLED__")
    # Cancelled live jobs mark stopped, not a failure popup path.
    assert "Stopped" in dialog._status.text() or dialog._close_btn.isEnabled()
    dialog.close()
    bridge.deleteLater()
    dialog.deleteLater()
