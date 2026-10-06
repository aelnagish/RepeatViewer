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


def test_process_monitor_bridge_closes_when_finished(qapp) -> None:
    dialog = ProcessMonitorDialog("Flexidot")
    buffer: list[str] = []
    bridge = ProcessMonitorBridge(
        dialog,
        buffer,
        live=True,
        title="Flexidot",
        cancelled_message="__CANCELLED__",
        close_when_finished=True,
    )
    dialog.prepare("Flexidot")
    assert dialog.isVisible()
    bridge.on_finished()
    assert not dialog.isVisible()
    bridge.deleteLater()
    dialog.deleteLater()


def test_process_monitor_bridge_keeps_window_on_failure(qapp) -> None:
    dialog = ProcessMonitorDialog("BLAST")
    buffer: list[str] = []
    bridge = ProcessMonitorBridge(
        dialog,
        buffer,
        live=True,
        title="BLAST",
        cancelled_message="__CANCELLED__",
        close_when_finished=True,
    )
    dialog.prepare("BLAST")
    bridge.on_failed("boom")
    assert dialog.isVisible()
    dialog.close()
    bridge.deleteLater()
    dialog.deleteLater()


def test_close_other_process_monitors_keeps_current_tool(qapp) -> None:
    from repeatviewer.main_window import MainWindow

    class Host:
        def __init__(self) -> None:
            self._flexidot_process_dialog = ProcessMonitorDialog("Flexidot")
            self._blast_process_dialog = ProcessMonitorDialog("BLAST")
            self._flexidot_process_dialog.prepare("Flexidot")
            self._blast_process_dialog.prepare("BLAST")

    host = Host()
    assert host._flexidot_process_dialog.isVisible()
    assert host._blast_process_dialog.isVisible()
    MainWindow._close_other_process_monitors(host, "blast")
    assert not host._flexidot_process_dialog.isVisible()
    assert host._blast_process_dialog.isVisible()
    host._blast_process_dialog.close()
    host._flexidot_process_dialog.deleteLater()
    host._blast_process_dialog.deleteLater()


def test_dismiss_open_tool_dialogs_closes_others(qapp) -> None:
    from PySide6.QtWidgets import QDialog

    from repeatviewer.main_window import MainWindow

    class Host:
        def __init__(self) -> None:
            self._te_viz_input_dialog = QDialog()
            self._blast_params_dialog = QDialog()
            self._hmmer_params_dialog = QDialog()
            self._te_viz_input_dialog.show()
            self._blast_params_dialog.show()
            self._hmmer_params_dialog.show()

    host = Host()
    keep = host._blast_params_dialog
    MainWindow._dismiss_open_tool_dialogs(host, keep=keep)
    assert not host._te_viz_input_dialog.isVisible()
    assert not host._hmmer_params_dialog.isVisible()
    assert host._blast_params_dialog.isVisible()
    host._blast_params_dialog.close()
    host._te_viz_input_dialog.deleteLater()
    host._blast_params_dialog.deleteLater()
    host._hmmer_params_dialog.deleteLater()


def test_process_window_defaults_are_show_and_close() -> None:
    from repeatviewer.job_options import (
        BlastOptions,
        DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED,
        DEFAULT_SHOW_PROCESS_WINDOW,
        FlexidotOptions,
        HmmerOptions,
    )

    assert DEFAULT_SHOW_PROCESS_WINDOW is True
    assert DEFAULT_CLOSE_PROCESS_WINDOW_WHEN_FINISHED is True
    assert BlastOptions().show_process_window is True
    assert BlastOptions().close_process_window_when_finished is True
    assert FlexidotOptions().show_process_window is True
    assert FlexidotOptions().close_process_window_when_finished is True
    assert HmmerOptions().show_process_window is True
    assert HmmerOptions().close_process_window_when_finished is True
