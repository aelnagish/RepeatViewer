"""Self dotplot preview tab for single-sequence FASTA files (TEAid BLAST)."""

from __future__ import annotations

import tempfile
from pathlib import Path

from PySide6.QtCore import Qt, QThread, QTimer
from PySide6.QtWidgets import (
    QHBoxLayout,
    QLabel,
    QSizePolicy,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

from teviewer.backends.similarity_annotation import load_use_system_blast_from_settings
from teviewer.backends.teaid_plotter import generate_preview_self_dotplot_png
from teviewer.job_options import (
    DEFAULT_TEAID_SELF_BLAST_WORD_SIZE,
    MAX_TEAID_SELF_BLAST_WORD_SIZE,
    MIN_TEAID_SELF_BLAST_WORD_SIZE,
)
from teviewer.workers.base import Worker

_DEBOUNCE_MS = 400
_DEFAULT_WORD_SIZE = DEFAULT_TEAID_SELF_BLAST_WORD_SIZE


class SelfDotplotPanel(QWidget):
    """Self dotplot preview tab for single-sequence FASTA files (TEAid BLAST)."""

    def __init__(self, path: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._path = path.expanduser().resolve()
        self._render_token = 0
        self._output_png = self._output_path_for_token(0)
        self._thread: QThread | None = None
        self._worker: Worker | None = None
        self._worker_busy = False
        self._rerun_pending = False
        self._closed = False
        self._image_preview = None

        self._debounce = QTimer(self)
        self._debounce.setSingleShot(True)
        self._debounce.setInterval(_DEBOUNCE_MS)
        self._debounce.timeout.connect(self.refresh)

        self.setObjectName("SelfDotplotPanel")

        root = QHBoxLayout(self)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)

        self._image_host = QWidget(self)
        image_layout = QVBoxLayout(self._image_host)
        image_layout.setContentsMargins(0, 0, 0, 0)
        self._status_label = QLabel("Open the Dotplot tab to generate.")
        self._status_label.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._status_label.setWordWrap(True)
        image_layout.addWidget(self._status_label, stretch=1)
        self._image_host.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        root.addWidget(self._image_host, stretch=1)

        controls = QWidget(self)
        controls.setObjectName("DotplotControls")
        controls.setFixedWidth(108)
        controls_layout = QVBoxLayout(controls)
        controls_layout.setContentsMargins(8, 8, 8, 8)
        controls_layout.setSpacing(6)

        word_label = QLabel("Word size")
        word_label.setAlignment(Qt.AlignmentFlag.AlignLeft)
        controls_layout.addWidget(word_label)

        self._word_size = QSpinBox()
        self._word_size.setRange(MIN_TEAID_SELF_BLAST_WORD_SIZE, MAX_TEAID_SELF_BLAST_WORD_SIZE)
        self._word_size.setValue(_DEFAULT_WORD_SIZE)
        self._word_size.setToolTip(
            "Self-BLAST word size (TEAid). Same algorithm as the TEAid self dot plot."
        )
        self._word_size.setFixedWidth(72)
        controls_layout.addWidget(self._word_size)
        controls_layout.addStretch()

        root.addWidget(controls)

        self._word_size.valueChanged.connect(self._schedule_refresh)

    def _output_path_for_token(self, token: int) -> Path:
        stem = self._path.stem or "sequence"
        return Path(tempfile.gettempdir()) / f"teviewer_dotplot_{stem}_{token}.png"

    def _work_dir_for_token(self, token: int) -> Path:
        stem = self._path.stem or "sequence"
        return Path(tempfile.gettempdir()) / f"teviewer_dotplot_{stem}_{token}_work"

    def load_path(self, path: Path) -> None:
        self._path = path.expanduser().resolve()
        self._render_token += 1
        self._output_png = self._output_path_for_token(self._render_token)
        self.refresh()

    def _schedule_refresh(self, *_args) -> None:
        if self._closed:
            return
        self._debounce.start()

    def _set_busy(self, message: str) -> None:
        self._status_label.setText(message)
        self._status_label.show()
        if self._image_preview is not None:
            self._image_preview.hide()

    def _show_png(self, png_path: Path) -> None:
        from teviewer.widgets.preview_tabs import ScaledImagePreview

        if not png_path.is_file():
            self._set_busy("Dotplot image was not created.")
            return

        if self._image_preview is None:
            self._image_preview = ScaledImagePreview(self._image_host)
            layout = self._image_host.layout()
            if layout is not None:
                layout.insertWidget(0, self._image_preview, stretch=1)
        self._status_label.hide()
        self._image_preview.show()
        if not self._image_preview.set_image_file(png_path):
            self._set_busy(f"Could not load dotplot:\n{png_path}")
            return
        self._image_preview.refresh_scale()

    def refresh(self) -> None:
        if self._closed:
            return
        if self._worker_busy:
            self._rerun_pending = True
            return
        self._begin_render()

    def _begin_render(self) -> None:
        if self._closed:
            return
        self._rerun_pending = False
        self._render_token += 1
        token = self._render_token
        self._output_png = self._output_path_for_token(token)
        if not self._path.is_file():
            self._set_busy("File not found.")
            return
        self._set_busy("Running self-BLAST and building dotplot…")
        self._start_worker(token)

    def _start_worker(self, token: int) -> None:
        work_dir = self._work_dir_for_token(token)
        worker = Worker(
            generate_preview_self_dotplot_png,
            self._path,
            self._output_png,
            word_size=self._word_size.value(),
            use_system_blast=load_use_system_blast_from_settings(),
            work_dir=work_dir,
        )
        thread = QThread(self)
        worker.moveToThread(thread)
        thread.started.connect(worker.run)
        worker.signals.finished.connect(
            lambda result, t=token: self._on_dotplot_finished(t, result),
        )
        worker.signals.error.connect(
            lambda message, t=token: self._on_dotplot_error(t, message),
        )
        worker.signals.finished.connect(thread.quit)
        worker.signals.error.connect(thread.quit)
        thread.finished.connect(worker.deleteLater)
        thread.finished.connect(thread.deleteLater)
        thread.finished.connect(self._on_worker_thread_finished)
        self._worker = worker
        self._thread = thread
        self._worker_busy = True
        thread.start()

    def _on_worker_thread_finished(self) -> None:
        self._worker_busy = False
        self._thread = None
        self._worker = None
        if self._rerun_pending and not self._closed:
            QTimer.singleShot(0, self._begin_render)

    def _on_dotplot_finished(self, token: int, result: object) -> None:
        if self._closed or token != self._render_token:
            return
        png_path = Path(str(result))
        self._show_png(png_path)

    def _on_dotplot_error(self, token: int, message: str) -> None:
        if self._closed or token != self._render_token:
            return
        self._set_busy(f"Dotplot failed:\n{message}")

    def showEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().showEvent(event)
        if self._image_preview is not None and self._image_preview.isVisible():
            QTimer.singleShot(0, self._image_preview.refresh_scale)

    def closeEvent(self, event) -> None:  # noqa: N802 (Qt API)
        self._closed = True
        self._debounce.stop()
        self._rerun_pending = False
        super().closeEvent(event)

    def open_dotplot_tab(self) -> None:
        self.refresh()

    def apply_preview_zoom(self, percent: int, fit_mode) -> None:
        if self._image_preview is not None:
            self._image_preview.set_zoom_percent(percent)
            self._image_preview.set_fit_mode(fit_mode)
