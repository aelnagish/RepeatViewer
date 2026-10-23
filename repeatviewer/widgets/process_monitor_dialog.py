"""Non-modal window that streams a long-running tool's command and progress."""

from __future__ import annotations

import re
from collections.abc import Mapping, Sequence

from PySide6.QtCore import QObject, Qt, Signal, Slot
from PySide6.QtGui import QFont
from PySide6.QtWidgets import (
    QDialog,
    QHBoxLayout,
    QLabel,
    QPlainTextEdit,
    QPushButton,
    QVBoxLayout,
    QWidget,
)


# Matches "4 hit(s)", "4 BLASTN hit(s)", "4 pyhmmer domain hit(s)", etc.
_HIT_COUNT_RE = re.compile(r"(\d+)\s+(?:\S+\s+)*hit\(s\)", re.IGNORECASE)


def format_annotation_run_summary(
    *,
    n_sequences: int,
    counts_by_type: Mapping[str, int] | None = None,
    preferred_types: Sequence[str] | None = None,
    total_new: int | None = None,
) -> list[str]:
    """Plain-language end-of-job lines for Feature / Similarity annotation.

    Example::

        Finished annotating 15 sequences.
        Found 20 new features (LTR 5, PBS 15).
    """
    n_sequences = max(0, int(n_sequences))
    counts = {
        str(key).strip(): max(0, int(value))
        for key, value in dict(counts_by_type or {}).items()
        if str(key).strip()
    }
    if total_new is None:
        total = sum(counts.values())
    else:
        total = max(0, int(total_new))

    lines: list[str] = []
    if n_sequences <= 0:
        lines.append("Finished annotation.")
    elif n_sequences == 1:
        lines.append("Finished annotating 1 sequence.")
    else:
        lines.append(f"Finished annotating {n_sequences} sequences.")

    preferred = [str(t).strip() for t in (preferred_types or ()) if str(t).strip()]
    ordered: list[str] = []
    seen: set[str] = set()
    for key in preferred:
        if key in seen:
            continue
        seen.add(key)
        ordered.append(key)
    for key in sorted(counts):
        if key in seen:
            continue
        if counts[key] <= 0 and key not in preferred:
            continue
        seen.add(key)
        ordered.append(key)

    # Drop preferred types that were never selected and have no hits.
    ordered = [key for key in ordered if counts.get(key, 0) > 0 or key in preferred]

    if total <= 0 and not ordered:
        lines.append("Found no new features.")
        return lines
    if total <= 0:
        # Selected types ran but none produced new features.
        zero_parts = [f"{key} 0" for key in ordered]
        if zero_parts:
            lines.append(
                "Found no new features"
                + (f" ({', '.join(zero_parts)})." if len(zero_parts) > 1 else ".")
            )
        else:
            lines.append("Found no new features.")
        return lines

    parts = [f"{key} {counts.get(key, 0)}" for key in ordered if counts.get(key, 0) > 0]
    # Include explicit zeros for selected types that found nothing when several ran.
    if len(preferred) > 1:
        parts = [f"{key} {counts.get(key, 0)}" for key in ordered]

    feat_word = "feature" if total == 1 else "features"
    if not parts:
        lines.append(f"Found {total} new {feat_word}.")
    elif len(parts) == 1 and counts.get(ordered[0], 0) == total:
        lines.append(f"Found {total} new {feat_word} ({parts[0]}).")
    else:
        lines.append(f"Found {total} new {feat_word} ({', '.join(parts)}).")
    return lines


def annotation_job_summary_lines(result: object) -> list[str]:
    """Lines to append when Feature / Similarity / Transfer jobs finish.

    Prefer the job ``message`` (module detail) and add a plain-language
    sequence / per-type feature summary when it can be derived.
    """
    if result is None:
        return []
    # Transfer annotation returns ``(AnnotationTransferResult, …)``.
    if isinstance(result, tuple) and result:
        result = result[0]

    lines: list[str] = []
    message = getattr(result, "message", None)
    if isinstance(message, str):
        for raw in message.splitlines():
            text = raw.strip()
            if text:
                lines.append(text)

    type_name = type(result).__name__

    if type_name == "TeAnnotationBatchResult":
        # Full plain summary is already streamed via emit_process_log in the batch.
        return []

    if type_name == "AutoAnnotationResult":
        features = getattr(result, "features", None) or ()
        mode = getattr(result, "mode", None)
        mode_name = str(getattr(mode, "value", mode) or "").strip()
        counts = {}
        for feat in features:
            key = str(getattr(feat, "feature_type", "") or "").strip() or "others"
            counts[key] = counts.get(key, 0) + 1
        preferred = [mode_name] if mode_name else list(counts)
        summary = format_annotation_run_summary(
            n_sequences=1,
            counts_by_type=counts,
            preferred_types=preferred,
            total_new=len(features),
        )
        for line in summary:
            if line not in lines:
                lines.append(line)
        return lines

    if type_name == "MultiAutoAnnotationResult":
        features = getattr(result, "features", None) or ()
        modes = getattr(result, "modes", None) or ()
        preferred = [
            str(getattr(mode, "value", mode) or "").strip()
            for mode in modes
            if str(getattr(mode, "value", mode) or "").strip()
        ]
        type_names = set(preferred)
        counts: dict[str, int] = {key: 0 for key in preferred}
        matched = 0
        for feat in features:
            key = str(getattr(feat, "feature_type", "") or "").strip()
            if type_names and key not in type_names:
                continue
            if not key:
                key = "others"
            counts[key] = counts.get(key, 0) + 1
            matched += 1
        summary = format_annotation_run_summary(
            n_sequences=1,
            counts_by_type=counts,
            preferred_types=preferred,
            total_new=matched if type_names else len(features),
        )
        for line in summary:
            if line not in lines:
                lines.append(line)
        return lines

    if type_name == "SimilarityAnnotationResult":
        total_hits = 0
        found = False
        for line in lines:
            match = _HIT_COUNT_RE.search(line)
            if match:
                total_hits += int(match.group(1))
                found = True
        if found:
            summary = format_annotation_run_summary(
                n_sequences=1,
                counts_by_type={"similarity hit": total_hits},
                preferred_types=("similarity hit",),
                total_new=total_hits,
            )
            for line in summary:
                if line not in lines:
                    lines.append(line)
        return lines

    if type_name == "AnnotationTransferResult":
        # Detailed transfer message is already streamed via progress.
        lines.clear()
        lifted = getattr(result, "lifted", None) or []
        identified = sum(
            1 for item in lifted if getattr(item, "feature", None) is not None
        )
        summary = format_annotation_run_summary(
            n_sequences=1,
            counts_by_type={"transferred": identified},
            preferred_types=("transferred",),
            total_new=identified,
        )
        return summary

    return lines


class ProcessMonitorDialog(QDialog):
    """Show a command line first, then append intermediate process output."""

    def __init__(self, title: str = "Process", parent=None) -> None:
        super().__init__(parent)
        self.setWindowTitle(title)
        self.setModal(False)
        self.setWindowFlag(Qt.WindowType.Window, True)
        self.resize(720, 420)

        self._status = QLabel("Running…", self)
        self._status.setStyleSheet("color: palette(mid);")

        self._log = QPlainTextEdit(self)
        self._log.setReadOnly(True)
        self._log.setLineWrapMode(QPlainTextEdit.LineWrapMode.WidgetWidth)
        font = QFont("Menlo")
        if not font.exactMatch():
            font = QFont("Consolas")
        if not font.exactMatch():
            font = QFont("Courier New")
        font.setStyleHint(QFont.StyleHint.Monospace)
        font.setPointSize(11)
        self._log.setFont(font)

        self._close_btn = QPushButton("Close", self)
        self._close_btn.setEnabled(False)
        self._close_btn.clicked.connect(self.close)

        footer = QHBoxLayout()
        footer.addWidget(self._status, 1)
        footer.addWidget(self._close_btn)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(12, 12, 12, 12)
        layout.setSpacing(8)
        layout.addWidget(self._log, 1)
        layout.addLayout(footer)

    def prepare(self, title: str) -> None:
        self.setWindowTitle(title)
        self._log.clear()
        self._status.setText("Running…")
        self._close_btn.setEnabled(False)
        self.show()
        self.raise_()
        self.activateWindow()

    def _is_scrolled_to_bottom(self) -> bool:
        bar = self._log.verticalScrollBar()
        # Near bottom (small slack) → keep following new output.
        return bar.maximum() <= 0 or bar.value() >= bar.maximum() - 8

    def append_line(self, text: str) -> None:
        if text is None:
            return
        line = str(text).rstrip("\r\n")
        bar = self._log.verticalScrollBar()
        stick = self._is_scrolled_to_bottom()
        saved = bar.value()
        self._log.appendPlainText(line)
        if stick:
            # Follow new lines only when the user was already at the bottom.
            bar.setValue(bar.maximum())
        else:
            # User scrolled up — keep that position (don't jump).
            bar.setValue(saved)

    def mark_finished(self, ok: bool, detail: str = "") -> None:
        if ok:
            self._status.setText("Finished.")
            self.append_line("")
            self.append_line("— finished —")
        else:
            message = detail.strip() or "failed"
            self._status.setText("Failed.")
            self.append_line("")
            self.append_line(f"— failed: {message} —")
        self._close_btn.setEnabled(True)

    def show_failure_log(
        self,
        title: str,
        lines: list[str],
        *,
        error: str = "",
    ) -> None:
        """Open (or reveal) this window with a captured command log after a failure."""
        self.setWindowTitle(title)
        self._log.clear()
        for line in lines:
            self.append_line(line)
        self.mark_finished(False, error)
        self.show()
        self.raise_()
        self.activateWindow()


class ProcessMonitorBridge(QObject):
    """Marshal worker log/finish signals onto the GUI thread.

    Worker signals connected to plain lambdas use a DirectConnection (no
    receiver QObject), so dialog.show() would run on the QThread and abort
    on macOS (NSWindow / NSPanel must be used from the main thread).
    Parent this bridge to the main window and connect with QueuedConnection.
    """

    def __init__(
        self,
        dialog: ProcessMonitorDialog,
        buffer: list[str],
        *,
        live: bool,
        title: str,
        cancelled_message: str,
        close_when_finished: bool = True,
        parent: QObject | None = None,
    ) -> None:
        super().__init__(parent)
        self._dialog = dialog
        self._buffer = buffer
        self._live = live
        self._title = title
        self._cancelled_message = cancelled_message
        self._close_when_finished = bool(close_when_finished)

    @Slot(str)
    def on_line(self, text: str) -> None:
        line = str(text)
        self._buffer.append(line)
        if self._live and self._dialog.isVisible():
            self._dialog.append_line(line)

    @Slot(int, str)
    def on_progress(self, _value: int, label: str) -> None:
        if label:
            self.on_line(label)

    @Slot()
    def on_finished(self) -> None:
        if self._live:
            self._dialog.mark_finished(True)
            if self._close_when_finished:
                self._dialog.close()

    @Slot(str)
    def on_failed(self, error: str) -> None:
        if error == self._cancelled_message:
            if self._live:
                self._dialog.mark_finished(False, "Stopped")
            return
        if self._live:
            self._dialog.mark_finished(False, str(error))
            self._dialog.raise_()
            self._dialog.activateWindow()
        else:
            self._dialog.show_failure_log(
                self._title, list(self._buffer), error=str(error)
            )


def attach_annotation_worker_process_monitor(
    host: QWidget,
    worker: QObject,
    *,
    title: str,
    enabled: bool,
    close_when_finished: bool = True,
    tool_key: str = "annotation",
) -> None:
    """Show a process window for Feature detection / Similarity ``Worker`` jobs.

    ``Worker`` emits ``progress(str)``, ``finished(object)``, ``error(str)`` —
    bridge those onto the same dialog used by BLAST / TE-Aid.
    """
    from repeatviewer.backends.job_cancel import CANCELLED_MESSAGE

    win = host.window()
    parent = win if isinstance(win, QWidget) else host
    # Walk to MainWindow when host is a floating Preview duplicate.
    cursor: QWidget | None = parent
    attach = None
    while cursor is not None:
        candidate = getattr(cursor, "_attach_process_monitor", None)
        if callable(candidate):
            attach = candidate
            parent = cursor
            break
        cursor = cursor.parentWidget()

    if callable(attach):
        # Prefer MainWindow so Stop / other-tool closing stays coordinated.
        class _Adapter(QObject):
            progress = Signal(int, str)
            finished = Signal()
            failed = Signal(str)
            log_line = Signal(str)

            def __init__(self) -> None:
                super().__init__(parent)
                queued = Qt.ConnectionType.QueuedConnection
                # BLAST / hmmscan command lines (process_log_sink → Worker.log_line).
                if hasattr(worker, "log_line"):
                    worker.log_line.connect(self.log_line, queued)
                worker.progress.connect(self._on_progress, queued)
                worker.finished.connect(self._on_finished, queued)
                worker.error.connect(self.failed, queued)

            @Slot(str)
            def _on_progress(self, message: str) -> None:
                text = str(message or "").strip()
                if text:
                    self.log_line.emit(text)
                    self.progress.emit(0, text)

            @Slot(object)
            def _on_finished(self, result: object) -> None:
                summary = annotation_job_summary_lines(result)
                if summary:
                    self.log_line.emit("")
                    for line in summary:
                        self.log_line.emit(line)
                self.finished.emit()

        adapter = _Adapter()
        # Keep adapter alive for the job lifetime.
        setattr(parent, f"_{tool_key}_worker_adapter", adapter)
        attach(
            adapter,
            title=title,
            enabled=bool(enabled),
            tool_key=tool_key,
            mirror_progress=False,
            close_when_finished=bool(close_when_finished),
        )
        return

    # Fallback when no MainWindow (tests / standalone preview).
    if not enabled:
        return
    dialog = ProcessMonitorDialog(title, parent=parent)
    dialog.prepare(title)
    bridge = ProcessMonitorBridge(
        dialog,
        [],
        live=True,
        title=title,
        cancelled_message=CANCELLED_MESSAGE,
        close_when_finished=bool(close_when_finished),
        parent=parent,
    )

    class _FinishRelay(QObject):
        @Slot(object)
        def on_finished(self, result: object) -> None:
            summary = annotation_job_summary_lines(result)
            if summary:
                bridge.on_line("")
                for line in summary:
                    bridge.on_line(line)
            bridge.on_finished()

    relay = _FinishRelay(parent)
    queued = Qt.ConnectionType.QueuedConnection
    if hasattr(worker, "log_line"):
        worker.log_line.connect(bridge.on_line, queued)
    worker.progress.connect(bridge.on_line, queued)
    worker.finished.connect(relay.on_finished, queued)
    worker.error.connect(bridge.on_failed, queued)
    setattr(parent, f"_{tool_key}_process_dialog", dialog)
    setattr(parent, f"_{tool_key}_process_bridge", bridge)
    setattr(parent, f"_{tool_key}_process_relay", relay)