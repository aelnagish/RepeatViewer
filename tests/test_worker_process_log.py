"""Worker routes external_process logs to log_line for process monitors."""

from __future__ import annotations

from PySide6.QtCore import QThread
from PySide6.QtWidgets import QApplication

from repeatviewer.external_process import emit_process_log, log_command
from repeatviewer.workers.base import Worker


def test_worker_emits_command_lines_on_log_line(qapp: QApplication) -> None:
    lines: list[str] = []

    def _job(*, progress=None):
        log_command("blastn", ["blastn", "-query", "q.fa", "-db", "genome"])
        emit_process_log("stdout sample")
        if progress is not None:
            progress("1/1 Pfam")
        return "ok"

    worker = Worker(_job, progress=True)
    worker.log_line.connect(lambda text: lines.append(str(text)))
    thread = QThread()
    worker.moveToThread(thread)
    finished: list[object] = []
    worker.finished.connect(lambda result: finished.append(result))
    worker.finished.connect(thread.quit)
    thread.started.connect(worker.run)
    thread.start()
    while thread.isRunning():
        qapp.processEvents()
        thread.wait(50)
    qapp.processEvents()

    assert finished == ["ok"]
    assert "Command:" in lines
    assert any("blastn -query q.fa -db genome" in line for line in lines)
    assert "stdout sample" in lines
