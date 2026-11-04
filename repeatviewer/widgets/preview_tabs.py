"""Tabbed document preview (text, PDF, HTML, images, FASTA)."""

from __future__ import annotations

import gzip
import logging
import time
from pathlib import Path

from PySide6.QtCore import QEvent, QSize, QSizeF, Qt, QTimer, Signal, Slot
from PySide6.QtGui import (
    QFont,
    QKeyEvent,
    QKeySequence,
    QPainter,
    QPixmap,
    QShortcut,
    QTextCursor,
)
from PySide6.QtPdf import QPdfDocument
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QMenu,
    QMessageBox,
    QScrollArea,
    QSizePolicy,
    QStackedWidget,
    QTabBar,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.backends.bismark_cov import is_bismark_cov_file
from repeatviewer.backends.selection_preview import (
    is_selection_preview_path,
    selection_preview_sources,
)
from repeatviewer.backends.sequence_features import (
    format_multi_sequence_fasta_summary_text,
    format_single_sequence_summary_text,
    summarize_loaded_sequence,
    summarize_single_sequence_fasta,
)
from repeatviewer.browser_file_metadata import scan_fasta_headers
from repeatviewer.widgets.html_preview import _qobject_alive, create_html_preview_widget
from repeatviewer.widgets.sequence_map_viewer import SequenceMapPanel, _is_delete_key
from repeatviewer.widgets.preview_zoom import (
    MAX_PREVIEW_FONT_SIZE,
    MAX_PREVIEW_ZOOM_PERCENT,
    MIN_PREVIEW_FONT_SIZE,
    MIN_PREVIEW_ZOOM_PERCENT,
    DefaultFastaPreviewTab,
    PreviewFitMode,
    PreviewZoomKind,
    load_default_fasta_preview_tab,
    load_preview_fit_mode,
    load_show_sequence_summary,
    load_preview_font_size,
    load_preview_zoom_percent,
)

logger = logging.getLogger(__name__)

_FASTA_SUFFIXES = {".fa", ".fasta", ".fas", ".faa"}
_FASTA_PREVIEW_MAX_BYTES = 1024**3
# Show status-bar progress when opening files at least this large.
_PREVIEW_OPEN_PROGRESS_MIN_BYTES = 10 * 1024 * 1024
_FASTA_TOO_LARGE_MESSAGE = "The file size is too big to show"


def _count_fasta_sequences_fast(path: Path, *, stop_at: int | None = None) -> int:
    """Count FASTA records by header lines (no Biopython).

    Pass ``stop_at=2`` when only "one sequence or many?" matters -- the scan
    then stops at the second header instead of reading the whole file. Leave it
    ``None`` when the real total is displayed to the user.
    """
    return scan_fasta_headers(path, stop_at=stop_at)[0]


def _is_fasta_file(path: Path) -> bool:
    suffix = path.suffix.lower()
    return suffix in _FASTA_SUFFIXES or path.name.lower().endswith(".fa.gz")


def _is_fasta_too_large_for_preview(path: Path) -> bool:
    if not _is_fasta_file(path):
        return False
    try:
        return path.expanduser().resolve().stat().st_size > _FASTA_PREVIEW_MAX_BYTES
    except OSError:
        return False


def _should_show_preview_open_progress(path: Path) -> bool:
    try:
        return path.expanduser().resolve().stat().st_size >= _PREVIEW_OPEN_PROGRESS_MIN_BYTES
    except OSError:
        return False


def _is_single_fasta_file(path: Path) -> bool:
    if not _is_fasta_file(path):
        return False
    if _is_fasta_too_large_for_preview(path):
        return False
    # Only the one-vs-many answer is needed here, so stop at the second header.
    # On a multi-sequence library this returns after the first few KB instead of
    # reading the entire file on the UI thread.
    return _count_fasta_sequences_fast(path, stop_at=2) == 1


def _is_multi_fasta_file(path: Path) -> bool:
    """True when the FASTA has two or more records (overview Annotation UI)."""
    if not _is_fasta_file(path):
        return False
    if _is_fasta_too_large_for_preview(path):
        return False
    return _count_fasta_sequences_fast(path, stop_at=2) >= 2


def _is_sequence_map_fasta(path: Path) -> bool:
    """True when Annotation / Sequence Map tabs should be available."""
    if not _is_fasta_file(path):
        return False
    if _is_fasta_too_large_for_preview(path):
        return False
    return _count_fasta_sequences_fast(path, stop_at=1) >= 1


_TEXT_PREVIEW_SUFFIXES = _FASTA_SUFFIXES | {
    ".txt",
    ".text",
    ".tsv",
    ".csv",
    ".bed",
    ".vcf",
    ".gff",
    ".gff3",
    ".gtf",
    ".log",
    ".md",
    ".json",
    ".xml",
    ".out",
    ".sizes",
}
_TEXT_PREVIEW_MAX_CHARS = 500_000
_DEFAULT_TEXT_MAX_CHARS = 200_000
_EDITABLE_TEXT_MAX_CHARS = 2_000_000
_TRUNCATION_MARKER = "\n\n... (truncated preview"
_EMPTY_FILE_PREFIX = "(empty file)\n"
_IMAGE_PREVIEW_SUFFIXES = {".png", ".jpg", ".jpeg", ".svg", ".gif", ".bmp", ".webp"}
_HTML_PREVIEW_SUFFIXES = {".html", ".htm"}
_UNSUPPORTED_PREVIEW_SUFFIXES = {
    ".docx",
    ".doc",
    ".xlsx",
    ".xls",
    ".pptx",
    ".ppt",
    ".odt",
    ".ods",
    ".odp",
    ".zip",
    ".gz",
    ".bz2",
    ".7z",
    ".exe",
    ".dmg",
    ".deb",
    ".rpm",
    ".bam",
    ".sam",
    ".cram",
}


def _preview_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix == _CONTIG_PREVIEW_SUFFIX:
        return "contig"
    if suffix in _DOTPLOT_PREVIEW_SUFFIXES:
        return "dotplot"
    if _is_tree_path(path):
        return "tree"
    if suffix in _HTML_PREVIEW_SUFFIXES:
        return "html"
    if suffix in _FASTA_SUFFIXES or path.name.lower().endswith(".fa.gz"):
        if _is_sequence_map_fasta(path):
            return "fasta_annotation"
        return "fasta"
    if suffix == ".pdf":
        return "pdf"
    if suffix in _IMAGE_PREVIEW_SUFFIXES:
        return "image"
    return "text"


def _switching_html_preview(old_path: Path, new_path: Path) -> bool:
    old_kind = _preview_kind(old_path)
    new_kind = _preview_kind(new_path)
    html_like = {"html"}
    # Recreate the page when crossing HTML ↔ non-HTML so FASTA Text/Annotation
    # chrome and WebEngine hosts do not leave a blank leftover panel.
    return (old_kind in html_like) != (new_kind in html_like)


def _is_text_like_preview(path: Path) -> bool:
    suffix = path.suffix.lower()
    name = path.name.lower()
    if is_bismark_cov_file(path):
        return True
    from repeatviewer.backends.blast_output_names import is_blast_table_filename

    return suffix in _TEXT_PREVIEW_SUFFIXES or is_blast_table_filename(name)


def _format_byte_size(num_bytes: int) -> str:
    size = float(max(0, num_bytes))
    for unit in ("B", "KB", "MB", "GB", "TB"):
        if size < 1024 or unit == "TB":
            if unit == "B":
                return f"{int(size)} {unit}"
            return f"{size:.1f} {unit}"
        size /= 1024
    return f"{size:.1f} TB"


def _read_text_preview(path: Path, max_chars: int) -> tuple[str, bool]:
    """Return ``(text, truncated)`` for preview display."""
    resolved = path.expanduser().resolve()
    if resolved.name.lower().endswith(".gz"):
        with gzip.open(resolved, "rt", encoding="utf-8", errors="replace") as handle:
            text = handle.read(max_chars + 1)
        try:
            file_size = resolved.stat().st_size
            size_note = f"compressed size {_format_byte_size(file_size)}"
        except OSError:
            size_note = "file is large"
    else:
        with resolved.open(encoding="utf-8", errors="replace") as handle:
            text = handle.read(max_chars + 1)
        try:
            file_size = resolved.stat().st_size
            size_note = f"file size {_format_byte_size(file_size)}"
        except OSError:
            size_note = "file is large"
    if not text.strip():
        return f"(empty file)\n\n{resolved}", False
    if len(text) > max_chars:
        return (
            text[:max_chars]
            + "\n\n... (truncated preview — "
            + f"{size_note} is too big for full preview; "
            + f"showing the first {max_chars:,} characters)"
        ), True
    return text, False


def _write_text_file(path: Path, text: str) -> None:
    resolved = path.expanduser().resolve()
    if resolved.name.lower().endswith(".gz"):
        raise OSError("Saving compressed (.gz) text previews is not supported.")
    resolved.write_text(text, encoding="utf-8")


# Alignment files shown as a dotplot rather than a table of numbers.
_DOTPLOT_PREVIEW_SUFFIXES = frozenset({".paf", ".1aln"})
# Newick trees (.treefile, .nwk, .newick, .tre, .tree) open in the Tree tab.
from repeatviewer.backends.tree_layout import is_tree_file as _is_tree_path  # noqa: E402
# Map to Reference results (.rvmap) open in the contig view.
_CONTIG_PREVIEW_SUFFIX = ".rvmap"


def is_previewable_file(path: Path) -> bool:
    """Return True when RepeatViewer can show the file in the Preview panel."""
    resolved = path.expanduser()
    suffix = resolved.suffix.lower()
    name = resolved.name.lower()
    if suffix in _DOTPLOT_PREVIEW_SUFFIXES or _is_tree_path(resolved):
        return True
    if suffix == _CONTIG_PREVIEW_SUFFIX:
        return True
    if is_bismark_cov_file(resolved):
        return True
    if suffix in _UNSUPPORTED_PREVIEW_SUFFIXES:
        return False
    if suffix in _HTML_PREVIEW_SUFFIXES | {".pdf"} | _IMAGE_PREVIEW_SUFFIXES:
        return True
    from repeatviewer.backends.blast_output_names import is_blast_table_filename

    if suffix in _TEXT_PREVIEW_SUFFIXES or is_blast_table_filename(name):
        return True
    return False


def preview_path_priority(path: Path) -> int:
    """Lower values are preferred when auto-opening after a tool run."""
    name = path.name.lower()
    suffix = path.suffix.lower()
    if suffix in _HTML_PREVIEW_SUFFIXES:
        if "divergence" in name or "genome_map" in name:
            return 0
        return 1
    if suffix in _DOTPLOT_PREVIEW_SUFFIXES or _is_tree_path(path):
        return 2
    if suffix == _CONTIG_PREVIEW_SUFFIX:
        return 2
    if suffix == ".pdf":
        return 3
    if suffix in _IMAGE_PREVIEW_SUFFIXES:
        return 4
    from repeatviewer.backends.blast_output_names import is_blast_table_filename

    if suffix in _FASTA_SUFFIXES or is_blast_table_filename(name):
        return 5
    if is_bismark_cov_file(path):
        return 5
    if suffix in _TEXT_PREVIEW_SUFFIXES:
        return 6
    return 99


def _looks_like_binary_file(path: Path) -> bool:
    try:
        head = path.read_bytes()[:512]
    except OSError:
        return False
    if head.startswith(b"PK"):
        return True
    return b"\x00" in head


class ScaledImagePreview(QWidget):
    """PNG/JPEG/SVG preview with scroll bars (e.g. Flexidot / export figures).

    SVG is re-rendered with ``QSvgRenderer`` at the current zoom so it stays
    sharp (unlike ``QPixmap(svg)``, which rasterizes once then upscales).
    """

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self._original: QPixmap | None = None
        self._svg_renderer = None
        self._svg_size = QSize()
        self._zoom_percent = load_preview_zoom_percent("image")
        self._fit_mode = PreviewFitMode.FIT_WINDOW
        self._scale_pending = False
        self._last_physical_w = 0
        self._reset_scroll_on_refresh = True
        self._refresh_timer = QTimer(self)
        self._refresh_timer.setSingleShot(True)
        self._refresh_timer.setInterval(50)
        self._refresh_timer.timeout.connect(self._do_refresh_scale)

        self._scroll = QScrollArea(self)
        self._scroll.setObjectName("PreviewScrollArea")
        self._scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        self._scroll.setWidgetResizable(False)
        self._scroll.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop
        )
        self._scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )
        self._scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded
        )

        self._label = QLabel()
        self._label.setAlignment(
            Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignTop
        )
        self._scroll.setWidget(self._label)
        self._scroll.viewport().installEventFilter(self)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._scroll)

    def set_image_file(self, path: Path) -> bool:
        self._svg_renderer = None
        self._svg_size = QSize()
        if path.suffix.lower() == ".svg":
            try:
                from PySide6.QtSvg import QSvgRenderer
            except Exception:
                self._original = None
                self._label.clear()
                self._last_physical_w = 0
                return False
            renderer = QSvgRenderer(str(path))
            if not renderer.isValid():
                self._original = None
                self._label.clear()
                self._last_physical_w = 0
                return False
            size = renderer.defaultSize()
            if size.width() < 1 or size.height() < 1:
                vb = renderer.viewBox()
                size = QSize(max(1, int(vb.width())), max(1, int(vb.height())))
            self._svg_renderer = renderer
            self._svg_size = size
            # Aspect placeholder for fit math; pixels come from the renderer.
            self._original = QPixmap(size)
            self._original.fill(Qt.GlobalColor.white)
            self._last_physical_w = 0
            self.refresh_scale(immediate=True, reset_scroll=True)
            return True

        pix = QPixmap(str(path))
        if pix.isNull():
            self._original = None
            self._label.clear()
            self._last_physical_w = 0
            return False
        self._original = pix
        self._last_physical_w = 0
        self.refresh_scale(immediate=True, reset_scroll=True)
        return True

    def reload_image_file(self, path: Path) -> bool:
        """Re-read the file from disk (e.g. after Flexidot DPI post-processing)."""
        return self.set_image_file(path)

    def eventFilter(self, watched, event) -> bool:  # noqa: N802 (Qt API)
        try:
            if watched is self._scroll.viewport() and event.type() == QEvent.Type.Resize:
                self.refresh_scale(reset_scroll=False)
        except KeyboardInterrupt:
            # Ctrl+C during a heavy QPixmap scale must not tear down Qt's filter chain.
            pass
        return super().eventFilter(watched, event)

    def resizeEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().resizeEvent(event)
        self.refresh_scale(reset_scroll=False)

    def showEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().showEvent(event)
        # After PDF/HTML → image, scroll can be stale; keep zoom %, pin top-left.
        self.refresh_scale(immediate=True, reset_scroll=True)

    def set_zoom_percent(self, percent: int) -> None:
        # Keep zoom %; always pin to top-left so PDF → PNG cannot leave a blank view.
        self._zoom_percent = max(25, min(MAX_PREVIEW_ZOOM_PERCENT, int(percent)))
        self.refresh_scale(immediate=True, reset_scroll=True)

    def set_fit_mode(self, mode: PreviewFitMode | str) -> None:
        if isinstance(mode, str):
            try:
                mode = PreviewFitMode(mode)
            except ValueError:
                mode = PreviewFitMode.FIT_WINDOW
        self._fit_mode = mode
        self.refresh_scale(immediate=True, reset_scroll=True)

    def reset_scroll_to_top_left(self) -> None:
        """Pin scrollbars to the top-left corner of the image."""
        self._scroll.horizontalScrollBar().setValue(0)
        self._scroll.verticalScrollBar().setValue(0)

    def refresh_scale(
        self,
        *,
        immediate: bool = False,
        reset_scroll: bool | None = None,
    ) -> None:
        """Rescale the pixmap to the viewport.

        Resize storms debounce via a short timer so SmoothTransformation on large
        Flexidot images does not block the UI (and so Ctrl+C mid-scale is rare).
        """
        if reset_scroll is not None:
            self._reset_scroll_on_refresh = reset_scroll
        if immediate:
            self._refresh_timer.stop()
            self._do_refresh_scale()
            return
        self._refresh_timer.start()

    def _do_refresh_scale(self) -> None:
        if self._original is None or self._original.isNull():
            return
        viewport = self._scroll.viewport().size()
        if viewport.width() < 20 or viewport.height() < 20:
            self._scale_pending = True
            return
        self._scale_pending = False

        if self._svg_renderer is not None and self._svg_size.isValid():
            ow = max(1, self._svg_size.width())
            oh = max(1, self._svg_size.height())
        else:
            ow = self._original.width()
            oh = self._original.height()
        vw = viewport.width()
        vh = viewport.height()

        if self._fit_mode == PreviewFitMode.ORIGINAL:
            base_scale = 1.0
        elif self._fit_mode == PreviewFitMode.FIT_WIDTH:
            base_scale = vw / ow if ow > 0 else 1.0
        else:
            base_scale = min(vw / ow, vh / oh) if ow > 0 and oh > 0 else 1.0

        scale = base_scale * (self._zoom_percent / 100.0)
        target_w = max(1, int(ow * scale))
        target_h = max(1, int(oh * scale))
        dpr = max(1.0, self.devicePixelRatioF())
        physical_w = max(1, int(target_w * dpr))
        physical_h = max(1, int(target_h * dpr))
        if physical_w == self._last_physical_w:
            current = self._label.pixmap()
            if current is not None and not current.isNull():
                if self._reset_scroll_on_refresh:
                    self.reset_scroll_to_top_left()
                    self._reset_scroll_on_refresh = False
                return

        if self._svg_renderer is not None:
            # Re-rasterize from vectors at the display size (stays sharp on zoom).
            try:
                display = QPixmap(physical_w, physical_h)
                display.fill(Qt.GlobalColor.white)
                painter = QPainter(display)
                try:
                    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
                    painter.setRenderHint(
                        QPainter.RenderHint.SmoothPixmapTransform, True
                    )
                    self._svg_renderer.render(painter)
                finally:
                    painter.end()
            except KeyboardInterrupt:
                return
            display.setDevicePixelRatio(dpr)
        else:
            # Fast path while dragging/resizing large images; smooth once settled.
            transform = (
                Qt.TransformationMode.FastTransformation
                if max(ow, oh) >= 2500
                else Qt.TransformationMode.SmoothTransformation
            )
            try:
                display = self._original.scaledToWidth(physical_w, transform)
            except KeyboardInterrupt:
                return
            display.setDevicePixelRatio(dpr)

        self._last_physical_w = physical_w
        self._label.setPixmap(display)
        self._label.setFixedSize(target_w, target_h)
        if self._reset_scroll_on_refresh:
            self.reset_scroll_to_top_left()
            self._reset_scroll_on_refresh = False


class PdfPreviewWidget(QWidget):
    """PDF preview: fit-to-width pages, scaled by the PDF preview zoom percent."""

    _PAGE_GAP_PX = 10

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self._doc = QPdfDocument(self)
        self._render_generation = 0
        self._zoom_percent = load_preview_zoom_percent(PreviewZoomKind.PDF)
        self._pending_scroll_frac: tuple[float, float] | None = None

        self._scroll = QScrollArea(self)
        self._scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        self._scroll.setWidgetResizable(True)
        self._scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._scroll.viewport().installEventFilter(self)

        self._content = QWidget(self._scroll)
        self._layout = QVBoxLayout(self._content)
        self._layout.setContentsMargins(6, 6, 6, 6)
        self._layout.setSpacing(self._PAGE_GAP_PX)
        self._scroll.setWidget(self._content)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(self._scroll)

        self._render_timer = QTimer(self)
        self._render_timer.setSingleShot(True)
        self._render_timer.setInterval(100)
        self._render_timer.timeout.connect(self._render_pages)

        self._doc.statusChanged.connect(self._on_doc_status)
        if hasattr(self._doc, "pageCountChanged"):
            self._doc.pageCountChanged.connect(self._schedule_render)

    def set_zoom_percent(self, percent: int) -> None:
        hbar = self._scroll.horizontalScrollBar()
        vbar = self._scroll.verticalScrollBar()
        vw = max(1, self._scroll.viewport().width())
        vh = max(1, self._scroll.viewport().height())
        old_w = max(1, self._content.width())
        old_h = max(1, self._content.height())
        self._pending_scroll_frac = (
            (hbar.value() + vw / 2.0) / old_w,
            (vbar.value() + vh / 2.0) / old_h,
        )
        self._zoom_percent = max(
            MIN_PREVIEW_ZOOM_PERCENT,
            min(MAX_PREVIEW_ZOOM_PERCENT, int(percent)),
        )
        self._schedule_render()

    def eventFilter(self, watched, event) -> bool:  # noqa: N802 (Qt API)
        if watched is self._scroll.viewport() and event.type() == QEvent.Type.Resize:
            self._schedule_render()
        return super().eventFilter(watched, event)

    def resizeEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().resizeEvent(event)
        self._schedule_render()

    def _on_doc_status(self, status: QPdfDocument.Status) -> None:
        if status == QPdfDocument.Status.Ready:
            self._schedule_render()
            QTimer.singleShot(200, self._schedule_render)

    def _schedule_render(self) -> None:
        self._render_timer.start()

    def _clear_pages(self) -> None:
        while self._layout.count():
            item = self._layout.takeAt(0)
            widget = item.widget()
            if widget is not None:
                widget.deleteLater()

    @staticmethod
    def _page_size(doc: QPdfDocument, page: int) -> QSizeF:
        size = doc.pagePointSize(page)
        if size.width() <= 0 or size.height() <= 0:
            return QSizeF(612.0, 792.0)
        return size

    def _render_pages(self) -> None:
        if self._doc.status() != QPdfDocument.Status.Ready:
            return
        page_count = self._doc.pageCount()
        if page_count < 1:
            return

        viewport_w = max(200, self._scroll.viewport().width())
        zoom = max(0.25, self._zoom_percent / 100.0)
        # 100% = fit viewport width; other values scale that base width.
        logical_w = max(200, int(viewport_w * zoom))
        dpr = max(1.0, self.devicePixelRatioF())
        render_w = max(1, int(logical_w * dpr))

        generation = self._render_generation + 1
        self._render_generation = generation
        self._clear_pages()

        for page in range(page_count):
            if self._render_generation != generation:
                return
            page_size = self._page_size(self._doc, page)
            render_h = max(1, int(render_w * page_size.height() / page_size.width()))
            image = self._doc.render(page, QSize(render_w, render_h))
            label = QLabel(self._content)
            label.setAlignment(Qt.AlignmentFlag.AlignHCenter | Qt.AlignmentFlag.AlignTop)
            if image.isNull():
                label.setText(f"(Page {page + 1} could not be rendered)")
                label.setMinimumWidth(logical_w)
            else:
                pixmap = QPixmap.fromImage(image)
                pixmap.setDevicePixelRatio(dpr)
                label.setPixmap(pixmap)
                label.setFixedWidth(logical_w)
            self._layout.addWidget(label)

        self._content.adjustSize()
        if self._pending_scroll_frac is not None:
            frac_x, frac_y = self._pending_scroll_frac
            self._pending_scroll_frac = None
            hbar = self._scroll.horizontalScrollBar()
            vbar = self._scroll.verticalScrollBar()
            vw = max(1, self._scroll.viewport().width())
            vh = max(1, self._scroll.viewport().height())
            new_w = max(1, self._content.width())
            new_h = max(1, self._content.height())
            hbar.setValue(int(frac_x * new_w - vw / 2.0))
            vbar.setValue(int(frac_y * new_h - vh / 2.0))

    def load_file(self, path: Path) -> QPdfDocument.Error:
        self._clear_pages()
        resolved = path.expanduser().resolve()
        status = self._doc.load(str(resolved))
        if status == QPdfDocument.Error.None_:
            self._schedule_render()
        return status

    def refresh(self) -> None:
        """Re-render after tab show, resize, or zoom change."""
        self._schedule_render()


class FilePreviewPage(QWidget):
    """Single-file preview page inside a tab."""

    annotationsSaved = Signal(object)  # Path to a newly written annotation file
    sequenceExtracted = Signal(object)  # Path to a FASTA extracted from a selection
    openInAdditionalPreviewTab = Signal(object)  # Path to open as a new Preview tab
    dotplotRegionSelected = Signal(dict)  # payload from window.__flexidotLastClick
    fileWritten = Signal(object)  # Path of a file written from the preview (e.g. tree export)
    openFileRequested = Signal(object)  # Path a preview asks to open (contig view → Results)
    # Tree view → "Show in alignment": (alignment path, record ids, extract folder)
    alignmentRequested = Signal(object, object, object)
    # Overview right-click → "Show self-dotplot": (file path, sequence id), new tab
    selfDotplotTabRequested = Signal(object, str)
    blastRequested = Signal(object)  # FASTA a preview sends to BLAST (contig view right-click)
    # Source FASTAs of a temporary multi-file preview, to combine into a real file.
    saveSelectionAsFileRequested = Signal(list)

    def __init__(self, path: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.path = path.expanduser().resolve()
        self._single_fasta_mode = False
        self._multi_fasta_mode = False
        self._msa_aligned_mode = False  # multi-FASTA with equal lengths → residue rows
        self._annotation_panel: SequenceMapPanel | None = None
        self._annotation_overview = None  # MultiSequenceOverviewPanel, multi-FASTA only
        self._dotplot_panel = None  # built lazily on first view
        self._alignment_dotplot = None  # PAF/.1aln viewer, built on first use
        self._tree_panel = None  # Newick tree viewer, built on first use
        self._stack = QStackedWidget(self)
        self._stack.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self._placeholder = QLabel("Loading…")
        self._placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._base_text_font_pt = load_preview_font_size()
        self._text_truncated = False
        self._text_edit_mode = False
        self._text_dirty = False
        self._text_changed_at = 0.0  # time.monotonic() of the last typed change
        self._text_undone_at = 0.0  # time.monotonic() of the last undo of typing
        self._text_redo_steps = 0
        self._text = QTextEdit()
        self._text.setReadOnly(True)
        self._text.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self._text.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._text.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._text.setMinimumSize(200, 120)
        self._text.textChanged.connect(self._on_text_editor_changed)
        self._text_scroll = QScrollArea(self)
        self._text_scroll.setWidgetResizable(True)
        self._text_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        self._text_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._text_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._text_scroll.setWidget(self._text)

        self._fasta_text = QTextEdit()
        self._fasta_text.setReadOnly(True)
        self._fasta_text.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self._fasta_text.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._fasta_text.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._fasta_text.setMinimumSize(200, 120)
        self._fasta_text.textChanged.connect(self._on_text_editor_changed)
        # Right-click: the standard text menu plus Paste Reverse Complement.
        for editor in (self._text, self._fasta_text):
            editor.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)
            editor.customContextMenuRequested.connect(
                lambda pos, e=editor: self._show_text_context_menu(e, pos)
            )
        self._fasta_text_scroll = QScrollArea(self)
        self._fasta_text_scroll.setWidgetResizable(True)
        self._fasta_text_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        self._fasta_text_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._fasta_text_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._fasta_text_scroll.setWidget(self._fasta_text)

        self._pdf = PdfPreviewWidget(self)
        self._html = create_html_preview_widget(self)
        if hasattr(self._html, "dotplotRegionSelected"):
            self._html.dotplotRegionSelected.connect(self.dotplotRegionSelected.emit)
        self._image = ScaledImagePreview(self)
        # After _pdf/_html/_image exist — apply_text_font_size → preview_zoom_kind.
        self.apply_text_font_size(self._base_text_font_pt)
        self._fasta_stack = QStackedWidget(self)
        self._fasta_annotation_holder = QWidget(self)
        self._fasta_annotation_holder.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._fasta_dotplot_holder = QWidget(self)
        self._fasta_dotplot_holder.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._fasta_cyto_holder = QWidget(self)
        self._fasta_cyto_holder.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._cyto_panel = None  # built lazily on first view
        # Selection dragged on the dotplot, waiting to be applied to the
        # Annotation view. select_by_position() only retries for ~2 s while the
        # tab is hidden, so hold it here until Annotation is actually opened.
        self._pending_annotation_selection: tuple[int, int] | None = None
        # Prefer this Overview record when opening Dotplot (Self) from context menu.
        self._pending_dotplot_sequence_id: str | None = None
        # Annotation/overview data is loaded lazily on first Annotation tab open
        # so Text preview of large genomes does not wait on a full FASTA parse.
        self._annotation_data_loaded_for: Path | None = None
        # When Annotation is the active tab, skip filling the Text editor until
        # the user opens Text (setPlainText is costly when flipping many files).
        self._fasta_text_deferred: bool = False
        self._fasta_stack.addWidget(self._fasta_text_scroll)
        self._fasta_stack.addWidget(self._fasta_annotation_holder)
        self._fasta_stack.addWidget(self._fasta_dotplot_holder)
        self._fasta_stack.addWidget(self._fasta_cyto_holder)
        self._fasta_stack.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self._sequence_info_bar = QLabel(self)
        self._sequence_info_bar.setTextInteractionFlags(
            Qt.TextInteractionFlag.TextSelectableByMouse
        )
        self._sequence_info_bar.hide()

        for widget in (
            self._placeholder,
            self._text_scroll,
            self._pdf,
            self._html,
            self._image,
        ):
            widget.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
            self._stack.addWidget(widget)

        # One expanding host for FASTA vs other previews. Two stretch=1 siblings
        # used to leave a large blank region when mode switches left both claiming space.
        self._mode_stack = QStackedWidget(self)
        self._mode_stack.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._mode_stack.addWidget(self._stack)
        self._mode_stack.addWidget(self._fasta_stack)

        self._selection_banner = QFrame(self)
        self._selection_banner.setObjectName("selectionPreviewBanner")
        self._selection_banner.setStyleSheet(
            "QFrame#selectionPreviewBanner {"
            " background: rgba(255, 196, 0, 45);"
            " border-bottom: 1px solid palette(mid); }"
        )
        banner_layout = QHBoxLayout(self._selection_banner)
        banner_layout.setContentsMargins(8, 4, 8, 4)
        self._selection_banner_label = QLabel(self._selection_banner)
        self._selection_banner_label.setWordWrap(True)
        banner_layout.addWidget(self._selection_banner_label, stretch=1)
        self._selection_save_button = QPushButton(
            "Save as One File…", self._selection_banner
        )
        self._selection_save_button.clicked.connect(self._on_save_selection_as_file)
        banner_layout.addWidget(self._selection_save_button)
        self._selection_banner.hide()

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
        layout.addWidget(self._selection_banner)
        layout.addWidget(self._mode_stack, stretch=1)
        layout.addWidget(self._sequence_info_bar)
        self._load()
        self._apply_view_mode()
        self._install_annotation_delete_shortcuts()

    def _install_annotation_delete_shortcuts(self) -> None:
        for key in (QKeySequence.StandardKey.Delete, QKeySequence(Qt.Key.Key_Backspace)):
            shortcut = QShortcut(key, self)
            shortcut.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            shortcut.activated.connect(self._delete_annotation_features_if_active)

    def _delete_annotation_features_if_active(self) -> None:
        if not self.is_annotation_preview_active() or self._multi_fasta_mode:
            return
        if self._annotation_panel is None:
            return
        self._annotation_panel._map._sync_feature_selection_from_annot_track()
        self._annotation_panel._delete_feature(confirm_empty=False)

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802
        if (
            _is_delete_key(event.key())
            and self.is_annotation_preview_active()
            and not self._multi_fasta_mode
        ):
            if self._annotation_panel is not None:
                self._annotation_panel._delete_feature(confirm_empty=False)
                event.accept()
                return
        super().keyPressEvent(event)

    def _apply_view_mode(self) -> None:
        # Index 0 = general preview stack; index 1 = single-FASTA Text/Annotation.
        self._mode_stack.setCurrentIndex(1 if self._single_fasta_mode else 0)
        self._notify_preview_host()

    def _notify_preview_host(self) -> None:
        tab_host = self.parent()
        while tab_host is not None and not isinstance(tab_host, PreviewTabWidget):
            tab_host = tab_host.parent()
        if isinstance(tab_host, PreviewTabWidget):
            tab_host._sync_header_from_current_page()
            tab_host._sync_zoom_bar_enabled()

    def current_fasta_tab_index(self) -> int:
        if not self._single_fasta_mode:
            return -1
        return self._fasta_stack.currentIndex()

    def set_fasta_tab_index(self, index: int) -> None:
        if not self._single_fasta_mode or index < 0 or index >= self._fasta_stack.count():
            return
        if self._fasta_stack.currentIndex() == index:
            return
        self._fasta_stack.setCurrentIndex(index)
        self._on_fasta_tab_changed(index)

    def load_path(self, path: Path) -> None:
        """Show a different file in this page (reuse one preview tab)."""
        old_path = self.path
        self.path = path.expanduser().resolve()
        self._cyto_images_cache_path = None
        self._cyto_images_cache_count = None
        self._sequence_info_bar_cache_key = None
        if (
            _preview_kind(old_path) in {"html"}
            and _preview_kind(self.path) not in {"html"}
            and hasattr(self._html, "teardown")
        ):
            self._html.teardown()
        self._load()
        self._apply_view_mode()
        QTimer.singleShot(0, self._after_load)

    def is_single_fasta_preview(self) -> bool:
        return self._single_fasta_mode

    def is_multi_fasta_preview(self) -> bool:
        return bool(self._single_fasta_mode and self._multi_fasta_mode)

    def is_single_record_fasta_preview(self) -> bool:
        return bool(self._single_fasta_mode and not self._multi_fasta_mode)

    def is_text_preview_active(self) -> bool:
        if self._single_fasta_mode:
            return self._fasta_stack.currentIndex() == 0
        return self._stack.currentWidget() is self._text_scroll

    def is_annotation_preview_active(self) -> bool:
        return self._single_fasta_mode and self._fasta_stack.currentIndex() == 1

    def persist_annotation_view_for_session(self) -> None:
        """Write the current annotation viewport to QSettings for later reopen."""
        if not _is_sequence_map_fasta(self.path):
            return
        panel = self._annotation_panel
        if panel is None:
            return
        state = panel._map.view_state_for_session()
        if not state:
            return
        from repeatviewer.session_restore import save_annotation_view

        save_annotation_view(
            self.path,
            state,
            sequence_id=str(state.get("sequence_id") or ""),
        )

    def show_annotation_tab(self) -> None:
        """Switch this FASTA preview to the Annotation (or Overview) sub-tab."""
        if not self._single_fasta_mode:
            return
        self._activate_fasta_tab(DefaultFastaPreviewTab.ANNOTATION)

    def is_cyto_preview_active(self) -> bool:
        return self._single_fasta_mode and self._fasta_stack.currentIndex() == 3

    def has_cyto_images(self) -> bool:
        from repeatviewer.backends.cyto_images import image_count_for_sequence

        # Cache per path — header sync used to re-scan the folder on every
        # Text/Annotation/Dotplot click.
        cache_path = getattr(self, "_cyto_images_cache_path", None)
        cache_count = getattr(self, "_cyto_images_cache_count", None)
        try:
            resolved = self.path.expanduser().resolve()
        except OSError:
            resolved = self.path
        if cache_path == resolved and isinstance(cache_count, int):
            return cache_count > 0
        try:
            count = int(image_count_for_sequence(resolved))
        except Exception:  # noqa: BLE001
            count = 0
        self._cyto_images_cache_path = resolved
        self._cyto_images_cache_count = count
        return count > 0

    def _ensure_cyto_panel(self):
        from repeatviewer.widgets.cyto_images_ui import CytogeneticsPanel

        if self._cyto_panel is None:
            self._cyto_panel = CytogeneticsPanel(self.path, self)
            self._cyto_panel.showSourceRequested.connect(self._show_cyto_source)
            holder_layout = self._fasta_cyto_holder.layout()
            if holder_layout is None:
                holder_layout = QVBoxLayout(self._fasta_cyto_holder)
                holder_layout.setContentsMargins(0, 0, 0, 0)
                holder_layout.setSpacing(0)
            holder_layout.addWidget(self._cyto_panel, stretch=1)
        elif self._cyto_panel.path() != self.path:
            self._cyto_panel.set_path(self.path)
        return self._cyto_panel

    def _open_cyto_tab(self) -> None:
        self._ensure_cyto_panel().reload()

    def open_cyto_tab(self) -> None:
        """Switch this page to the Cytogenetics sub-tab."""
        if not self._single_fasta_mode:
            return
        self._ensure_cyto_panel()
        self.set_fasta_tab_index(3)
        self._notify_preview_host()

    def _show_cyto_source(self, ref) -> None:
        """Show source: Annotation tab, jumped to the probe (whole sequence if none)."""
        panel = self._ensure_annotation_panel()
        tab_host = self._preview_tab_host()
        if tab_host is not None and tab_host._annotation_zoom_bar is not None:
            panel.bind_annotation_zoom_bar(tab_host._annotation_zoom_bar)
        self.set_fasta_tab_index(1)
        self._notify_preview_host()
        if ref is not None:
            QTimer.singleShot(0, lambda r=ref: panel.show_feature(r))

    def is_dotplot_preview_active(self) -> bool:
        if self._single_fasta_mode:
            return self._fasta_stack.currentIndex() == 2
        alignment = getattr(self, "_alignment_dotplot", None)
        return alignment is not None and self._stack.currentWidget() is alignment


    def preview_zoom_kind(self):
        """Which persisted size setting applies to the current view.

        Text → font size (pt). HTML / PDF / image → zoom percent.
        """
        from repeatviewer.widgets.preview_zoom import PreviewZoomKind

        if (
            self.is_annotation_preview_active()
            or self.is_dotplot_preview_active()
            or self.is_tree_preview_active()
            or self.is_contig_preview_active()
        ):
            return PreviewZoomKind.TEXT
        if not self._single_fasta_mode:
            current = self._stack.currentWidget()
            pdf = getattr(self, "_pdf", None)
            html = getattr(self, "_html", None)
            image = getattr(self, "_image", None)
            text_scroll = getattr(self, "_text_scroll", None)
            if pdf is not None and current is pdf:
                return PreviewZoomKind.PDF
            if html is not None and current is html:
                return PreviewZoomKind.HTML
            if image is not None and current is image:
                return PreviewZoomKind.IMAGE
            if (text_scroll is not None and current is text_scroll) or self._uses_text_preview():
                return PreviewZoomKind.TEXT
        elif self.is_text_preview_active():
            return PreviewZoomKind.TEXT
        kind = _preview_kind(self.path)
        if kind == "html":
            return PreviewZoomKind.HTML
        if kind == "pdf":
            return PreviewZoomKind.PDF
        if kind == "image":
            return PreviewZoomKind.IMAGE
        return PreviewZoomKind.TEXT

    def _on_fasta_tab_changed(self, index: int) -> None:
        if not self._single_fasta_mode:
            return
        self._refresh_sequence_info_bar()
        if index == 0:
            QTimer.singleShot(0, self._open_text_tab)
        elif index == 1:
            QTimer.singleShot(0, self._open_annotation_tab)
        elif index == 2:
            QTimer.singleShot(0, self._open_dotplot_tab)
        elif index == 3:
            QTimer.singleShot(0, self._open_cyto_tab)
        tab_host = self.parent()
        while tab_host is not None and not isinstance(tab_host, PreviewTabWidget):
            tab_host = tab_host.parent()
        if isinstance(tab_host, PreviewTabWidget):
            tab_host._sync_zoom_bar_enabled()
            # _sync_zoom_bar_enabled already refreshes the FASTA tab strip.

    def _after_load(self) -> None:
        if self._single_fasta_mode:
            index = self._fasta_stack.currentIndex()
            if index == 0:
                self._open_text_tab()
            elif index == 1:
                # _load_single_fasta_preview already loaded + scheduled open;
                # only reopen when that path did not run (e.g. Text→Annotation).
                if not self._annotation_data_is_current():
                    self._open_annotation_tab()
            elif index == 2:
                self._open_dotplot_tab()
            elif index == 3:
                self._open_cyto_tab()
            self._notify_preview_host()
            return
        current = self._stack.currentWidget()
        if current is self._text_scroll:
            self._refresh_text_view()
        elif current is self._html:
            self._reload_html_if_needed()
        elif current is self._image:
            self._image.refresh_scale()
        elif current is self._pdf:
            self._pdf.refresh()

    def showEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().showEvent(event)
        if self._single_fasta_mode:
            index = self._fasta_stack.currentIndex()
            if index == 0:
                self.apply_text_font_size(load_preview_font_size())
                QTimer.singleShot(0, self._open_text_tab)
            elif index == 1:
                QTimer.singleShot(0, self._open_annotation_tab)
            return
        current = self._stack.currentWidget()
        if current is self._text_scroll:
            self.apply_text_font_size(load_preview_font_size())
            QTimer.singleShot(0, self._refresh_text_view)
            return
        if current in (self._html, self._image, self._pdf):
            self.apply_preview_zoom(
                load_preview_zoom_percent(self.preview_zoom_kind()),
                load_preview_fit_mode(),
            )
        if current is self._html:
            QTimer.singleShot(0, self._reload_html_if_needed)
        elif current is self._image:
            QTimer.singleShot(
                0,
                lambda: self._image.refresh_scale(immediate=True, reset_scroll=True),
            )
        elif current is self._pdf:
            QTimer.singleShot(0, self._pdf.refresh)

    def _reload_html_if_needed(self) -> None:
        if self.path.suffix.lower() in {".html", ".htm"}:
            self._schedule_html_load(self.path)

    def _schedule_html_load(self, path: Path, *, force: bool = False) -> None:
        self._stack.setCurrentWidget(self._html)
        resolved = path.expanduser().resolve()
        page = self

        def _load() -> None:
            if not _qobject_alive(page):
                return
            if page._stack.currentWidget() is not page._html:
                return
            page._html.load_html_file(resolved, force=force)

        QTimer.singleShot(0, _load)

    def apply_text_font_size(self, point_size: int) -> None:
        self._base_text_font_pt = max(
            MIN_PREVIEW_FONT_SIZE,
            min(MAX_PREVIEW_FONT_SIZE, int(point_size)),
        )
        self._apply_text_font()

    def _uses_text_preview(self) -> bool:
        return _is_text_like_preview(self.path)

    def _apply_text_font(self) -> None:
        size = self._base_text_font_pt
        mono = QFont("Courier New", size)
        if not mono.exactMatch():
            mono = QFont("Menlo", size)
        self._text.setFont(mono)
        self._fasta_text.setFont(mono)
        scroll = getattr(self, "_text_scroll", None)
        fasta_mode = getattr(self, "_single_fasta_mode", False)
        if scroll is not None and (
            self._stack.currentWidget() is scroll
            or (fasta_mode and self._fasta_stack.currentIndex() == 0)
        ):
            QTimer.singleShot(0, self._refresh_text_view)

    def apply_preview_zoom(self, percent: int, fit_mode: PreviewFitMode | str) -> None:
        self._image.set_zoom_percent(percent)
        self._image.set_fit_mode(fit_mode)
        if hasattr(self._html, "set_zoom_percent"):
            self._html.set_zoom_percent(percent)
        if hasattr(self._html, "set_fit_mode"):
            self._html.set_fit_mode(fit_mode)
        if hasattr(self._pdf, "set_zoom_percent"):
            self._pdf.set_zoom_percent(percent)
        if self._uses_text_preview() or self.is_text_preview_active():
            self._apply_text_font()

    def reload_from_disk(self) -> None:
        """Refresh the current view after the file on disk changed."""
        if self.has_unsaved_changes():
            # Never drop unsaved Edit-text or annotation changes on a silent
            # reload (same file clicked again, folder watcher, tab switch).
            # Leaving the file on purpose goes through
            # PreviewTabWidget.confirm_leave_page().
            return
        suffix = self.path.suffix.lower()
        if suffix in {".png", ".jpg", ".jpeg", ".svg", ".gif", ".bmp", ".webp"}:
            if self._stack.currentWidget() is self._image:
                self._image.reload_image_file(self.path)
            return
        if suffix in {".html", ".htm"}:
            if self._stack.currentWidget() is self._html:
                self._schedule_html_load(self.path, force=True)
            return
        if suffix == ".pdf":
            if self._stack.currentWidget() is self._pdf:
                self._pdf.load_file(self.path)
            return
        self._load()

    def release_shared_annotation_zoom_bar(self) -> None:
        """Return AnnotationZoomBar to the header before this page is destroyed.

        Multi-seq Overview reparents the shared bar into its toolbar. If the
        page is deleted without unbinding, Qt destroys the bar with the page
        and the next HTML/preview sync crashes (``already deleted``).
        """
        overview = getattr(self, "_annotation_overview", None)
        if overview is None:
            return
        unbind = getattr(overview, "_unbind_header_zoom_bar", None)
        if callable(unbind):
            try:
                unbind()
            except RuntimeError:
                pass

    def _set_single_fasta_mode(self, enabled: bool) -> None:
        if self._single_fasta_mode and not enabled:
            self.release_shared_annotation_zoom_bar()
        self._single_fasta_mode = enabled
        if not enabled:
            self._multi_fasta_mode = False
            self._msa_aligned_mode = False
        self._apply_view_mode()
        self._refresh_sequence_info_bar()

    def _refresh_sequence_info_bar(self) -> None:
        path = self.path
        if (
            self._single_fasta_mode
            and self.is_annotation_preview_active()
        ):
            self._sequence_info_bar.hide()
            if (
                not self._multi_fasta_mode
                and self._annotation_panel is not None
            ):
                self._annotation_panel.refresh_sequence_summary()
            return

        # Dotplot shows Length/GC or Sequence number in its Display rail.
        if self.is_dotplot_preview_active():
            self._sequence_info_bar.hide()
            self._sequence_info_bar_cache_key = None
            panel = getattr(self, "_dotplot_panel", None)
            if panel is None:
                panel = getattr(self, "_alignment_dotplot", None)
            if panel is not None and hasattr(panel, "refresh_sequence_summary"):
                panel.refresh_sequence_summary()
            return

        if not _is_fasta_file(path):
            self._sequence_info_bar.clear()
            self._sequence_info_bar.hide()
            self._sequence_info_bar_cache_key = None
            return

        if _is_fasta_too_large_for_preview(path):
            self._sequence_info_bar.hide()
            self._sequence_info_bar_cache_key = None
            return

        try:
            resolved = path.expanduser().resolve()
            mtime = resolved.stat().st_mtime_ns
        except OSError:
            resolved = path
            mtime = None
        show_summary = load_show_sequence_summary()
        cache_key = (str(resolved), mtime, show_summary, self._multi_fasta_mode)
        if (
            cache_key == getattr(self, "_sequence_info_bar_cache_key", None)
            and self._sequence_info_bar.isVisible()
            and self._sequence_info_bar.text()
        ):
            return

        # Prefer already-loaded Annotation data over re-reading the FASTA.
        if (
            not self._multi_fasta_mode
            and self._annotation_panel is not None
            and show_summary
        ):
            data = getattr(self._annotation_panel, "_loaded_data", None)
            if data is not None and getattr(data, "fasta_path", None) == resolved:
                summary = summarize_loaded_sequence(data.sequence_id, data.sequence)
                self._sequence_info_bar.setText(
                    format_single_sequence_summary_text(summary)
                )
                self._sequence_info_bar.show()
                self._sequence_info_bar_cache_key = cache_key
                return

        # Only need single-vs-multi for the label when summary is off or multi.
        sequence_count = _count_fasta_sequences_fast(
            path, stop_at=None if self._multi_fasta_mode else 2
        )
        if sequence_count <= 0:
            self._sequence_info_bar.hide()
            self._sequence_info_bar_cache_key = None
            return

        if sequence_count == 1 and show_summary:
            try:
                summary = summarize_single_sequence_fasta(path)
            except (OSError, ValueError):
                self._sequence_info_bar.hide()
                self._sequence_info_bar_cache_key = None
                return
            self._sequence_info_bar.setText(format_single_sequence_summary_text(summary))
        else:
            if self._multi_fasta_mode:
                sequence_count = _count_fasta_sequences_fast(path)
            self._sequence_info_bar.setText(
                format_multi_sequence_fasta_summary_text(sequence_count)
            )
        self._sequence_info_bar.show()
        self._sequence_info_bar_cache_key = cache_key

    def _ensure_annotation_overview(self):
        from repeatviewer.widgets.multi_sequence_overview import MultiSequenceOverviewPanel

        if self._annotation_overview is None:
            self._annotation_overview = MultiSequenceOverviewPanel(self.path, self)
            self._annotation_overview.annotationsSaved.connect(self.annotationsSaved.emit)
            self._annotation_overview.sequenceExtracted.connect(self.sequenceExtracted.emit)
            self._annotation_overview.openInAdditionalPreviewTab.connect(
                self.openInAdditionalPreviewTab.emit
            )
            self._annotation_overview.showSelfDotplotRequested.connect(
                self._on_overview_show_self_dotplot
            )
            self._annotation_overview.sequenceOrderChanged.connect(
                self._on_overview_sequence_order_changed
            )
            holder_layout = self._fasta_annotation_holder.layout()
            if holder_layout is None:
                holder_layout = QVBoxLayout(self._fasta_annotation_holder)
                holder_layout.setContentsMargins(0, 0, 0, 0)
                holder_layout.setSpacing(0)
            holder_layout.addWidget(self._annotation_overview, stretch=1)
        return self._annotation_overview

    def _show_annotation_content(self, *, multi: bool) -> None:
        """Show overview (multi) or sequence map (single) inside the Annotation tab."""
        tab_host = self._preview_tab_host()
        zoom_bar = (
            tab_host._annotation_zoom_bar
            if tab_host is not None
            else None
        )
        if multi:
            overview = self._ensure_annotation_overview()
            overview.set_annotation_editing_enabled(not self.is_selection_preview())
            overview.show()
            if self._annotation_panel is not None:
                self._annotation_panel.hide()
            if zoom_bar is not None:
                overview.bind_annotation_zoom_bar(zoom_bar)
        else:
            panel = self._ensure_annotation_panel()
            panel.show()
            if self._annotation_overview is not None:
                self._annotation_overview.hide()
            if zoom_bar is not None:
                zoom_bar.set_overview_mode(False)
                panel.bind_annotation_zoom_bar(zoom_bar)

    def _ensure_annotation_panel(self) -> SequenceMapPanel:
        if self._annotation_panel is None:
            self._annotation_panel = SequenceMapPanel(self.path, self, autoload=False)
            self._annotation_panel.annotationsSaved.connect(self.annotationsSaved.emit)
            self._annotation_panel.sequenceExtracted.connect(self.sequenceExtracted.emit)
            self._annotation_panel.activeSequenceChanged.connect(
                self._on_annotation_active_sequence_changed
            )
            holder_layout = self._fasta_annotation_holder.layout()
            if holder_layout is None:
                holder_layout = QVBoxLayout(self._fasta_annotation_holder)
                holder_layout.setContentsMargins(0, 0, 0, 0)
                holder_layout.setSpacing(0)
            holder_layout.addWidget(self._annotation_panel, stretch=1)
        return self._annotation_panel

    @Slot(str)
    def _on_annotation_active_sequence_changed(self, _seq_id: str) -> None:
        if self._dotplot_panel is None:
            return
        self._sync_dotplot_active_sequence(self._dotplot_panel)
        if self.is_dotplot_preview_active():
            self._dotplot_panel.ensure_loaded()

    def _show_alignment_dotplot(self, path: Path) -> None:
        """Plot a PAF/.1aln alignment file instead of dumping it as text."""
        from repeatviewer.widgets.dotplot_view import DotplotPanel

        if self._alignment_dotplot is None:
            self._alignment_dotplot = DotplotPanel(path, self)
            self._alignment_dotplot.setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
            )
            self._stack.addWidget(self._alignment_dotplot)
        self._stack.setCurrentWidget(self._alignment_dotplot)
        self._alignment_dotplot.load_alignment(path)

    def _show_tree(self, path: Path) -> None:
        """Draw a Newick tree in the Tree view instead of showing the text."""
        from repeatviewer.widgets.tree_view import TreePanel

        if self._tree_panel is None:
            self._tree_panel = TreePanel(None, self)
            self._tree_panel.setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
            )
            self._tree_panel.fileWritten.connect(self.fileWritten.emit)
            self._tree_panel.alignmentRequested.connect(self.alignmentRequested.emit)
            self._stack.addWidget(self._tree_panel)
        self._stack.setCurrentWidget(self._tree_panel)
        self._tree_panel.load_tree(path)

    def _show_contig(self, path: Path) -> None:
        """Draw a Map to Reference result (.rvmap) in the contig view."""
        from repeatviewer.widgets.contig_view import ContigViewPanel

        if getattr(self, "_contig_panel", None) is None:
            self._contig_panel = ContigViewPanel(None, self)
            self._contig_panel.setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding
            )
            self._contig_panel.openRequested.connect(self.openFileRequested.emit)
            self._contig_panel.blastRequested.connect(self.blastRequested.emit)
            self._contig_panel.fileWritten.connect(self.fileWritten.emit)
            self._stack.addWidget(self._contig_panel)
        self._stack.setCurrentWidget(self._contig_panel)
        self._contig_panel.load(path)

    def contig_panel(self):
        """The contig view of this page, if it was ever shown (else None)."""
        return getattr(self, "_contig_panel", None)

    def tree_panel(self):
        """The Tree view of this page, if it was ever shown (else None)."""
        return getattr(self, "_tree_panel", None)

    def is_tree_preview_active(self) -> bool:
        panel = getattr(self, "_tree_panel", None)
        return (
            not self._single_fasta_mode
            and panel is not None
            and self._stack.currentWidget() is panel
        )

    def is_contig_preview_active(self) -> bool:
        """Mapping result (.rvmap): the contig view zooms itself (−, +, wheel, pinch)."""
        panel = getattr(self, "_contig_panel", None)
        return (
            not self._single_fasta_mode
            and panel is not None
            and self._stack.currentWidget() is panel
        )

    def _ensure_dotplot_panel(self):
        """Create the dotplot panel on first use.

        Nothing is computed here -- DotplotPanel.ensure_loaded() does that when
        the tab is actually shown, so selecting a file never pays the cost.
        """
        if self._dotplot_panel is None:
            from repeatviewer.widgets.dotplot_view import DotplotPanel

            self._dotplot_panel = DotplotPanel(self.path, self)
            self._dotplot_panel.suggest_default_word_size()
            # Reuse the Flexidot Go-to chain: clicking a dot jumps the
            # Annotation tab to that region.
            # Selection is shared with the Annotation tab, but switching tabs
            # is left to the user -- clicking the plot never navigates.
            self._dotplot_panel.rangeSelected.connect(self._on_dotplot_range_selected)
            self._dotplot_panel.viewInAnnotationRequested.connect(
                self._on_dotplot_view_in_annotation
            )
            holder_layout = self._fasta_dotplot_holder.layout()
            if holder_layout is None:
                holder_layout = QVBoxLayout(self._fasta_dotplot_holder)
                holder_layout.setContentsMargins(0, 0, 0, 0)
                holder_layout.setSpacing(0)
            holder_layout.addWidget(self._dotplot_panel, stretch=1)
        else:
            self._dotplot_panel.set_path(self.path)
        return self._dotplot_panel

    def _on_overview_show_self_dotplot(self, sequence_id: str) -> None:
        """Overview context menu → Dotplot (Self) for the chosen sequence, in a new
        Preview tab so the alignment stays open in this one."""
        sid = str(sequence_id or "").strip()
        if not sid:
            return
        if getattr(self, "_dotplot_in_new_tab", False):     # set by the Preview tabs
            self.selfDotplotTabRequested.emit(self.path, sid)
            return
        self.show_self_dotplot(sid)

    def show_self_dotplot(self, sequence_id: str) -> None:
        """Switch this page to Dotplot (Self) for ``sequence_id``."""
        self._pending_dotplot_sequence_id = str(sequence_id)
        self._activate_fasta_tab(DefaultFastaPreviewTab.DOTPLOT)

    def _on_overview_sequence_order_changed(self) -> None:
        """Keep Dotplot Sequence list numbered like the Overview after Sorting."""
        if self._dotplot_panel is None:
            return
        self._sync_dotplot_active_sequence(self._dotplot_panel)

    def _overview_region_for_dotplot(self, panel) -> tuple[int, int] | None:
        """Overview region to mirror onto Dotplot, or ``None`` to skip highlight.

        Requires exactly one selected sequence (whole-row or regional) whose id
        matches the Dotplot's active sequence. Multi-sequence / full-column
        selections without a single matching row do not highlight the plot.
        """
        overview = self._annotation_overview
        if overview is None or panel is None:
            return None
        try:
            region = overview._bars.region_selection()
        except Exception:  # noqa: BLE001
            logger.exception("Could not read Overview selection")
            return None
        if region is None:
            return None
        dot_id = str(panel.sequence_id() or "").strip()
        if not dot_id:
            return None
        ids = overview.record_ids()
        regional = overview._bars.regional_rows()
        if regional is not None:
            if len(regional) != 1:
                return None
            row = next(iter(regional))
            if not (0 <= row < len(ids)) or ids[row] != dot_id:
                return None
            return int(region[0]), int(region[1])
        # Full-height column band: only when exactly one name-selected sequence
        # matches the Dotplot sequence.
        rows = overview._bars.selected_rows()
        if len(rows) != 1:
            return None
        row = int(rows[0])
        if not (0 <= row < len(ids)) or ids[row] != dot_id:
            return None
        return int(region[0]), int(region[1])

    def _open_dotplot_tab(self) -> None:
        self._capture_text_selection_into_pending()
        panel = self._ensure_dotplot_panel()
        self._sync_dotplot_active_sequence(panel)
        panel.ensure_loaded()
        # Carry the current Annotation/Overview selection across to the plot. A
        # pending selection that originated on the plot takes priority, so
        # switching back and forth does not overwrite what you just dragged.
        if self._pending_annotation_selection is not None:
            start, end = self._pending_annotation_selection
            panel.show_selection(int(start), int(end))
            return
        if self._multi_fasta_mode and self._annotation_overview is not None:
            mapped = self._overview_region_for_dotplot(panel)
            if mapped is not None:
                panel.show_selection(int(mapped[0]), int(mapped[1]))
            else:
                # Avoid leaving a highlight from another sequence or multi-select.
                panel.clear_selection()
            return
        if self._annotation_panel is not None:
            try:
                current = self._annotation_panel.selection_range()
            except Exception:  # noqa: BLE001
                logger.exception("Could not read Annotation selection")
                current = None
            if current:
                panel.show_selection(int(current[0]), int(current[1]))

    def _sync_dotplot_active_sequence(self, panel=None) -> None:
        """Point Dotplot (Self) at the active Annotation/Overview sequence."""
        if panel is None:
            panel = self._dotplot_panel
        if panel is None:
            return
        if self._multi_fasta_mode:
            overview = self._annotation_overview
            if overview is None:
                return
            try:
                ov_path = Path(overview._path).expanduser().resolve()
                page_path = self.path.expanduser().resolve()
            except OSError:
                return
            if ov_path != page_path:
                return
            record_ids = overview.record_ids_for_dotplot()
            if not record_ids:
                return
            pending = str(self._pending_dotplot_sequence_id or "").strip()
            self._pending_dotplot_sequence_id = None
            # Keep the Dotplot Sequence pick unless Overview explicitly asked
            # to open a self-dotplot (pending) or the plot has no sequence yet.
            current = str(panel.sequence_id() or "").strip()
            preferred = (
                pending
                or current
                or overview.selected_record_id()
                or record_ids[0]
            )
            # Ids only + lazy loader: avoid copying every MSA row into the combo.
            panel.set_record_choices(
                record_ids,
                ungap=bool(self._msa_aligned_mode),
                preferred_id=preferred,
                sequence_loader=overview.sequence_text_for_dotplot,
            )
            # Overview "Show self-dotplot" asks for that one sequence, not Grid.
            if pending:
                panel.open_single_for_sequence(pending)
            return
        if self._annotation_panel is None:
            return
        try:
            ann_path = self._annotation_panel.fasta_path().expanduser().resolve()
            page_path = self.path.expanduser().resolve()
        except OSError:
            return
        if ann_path != page_path:
            # Annotation still holds another file — do not override Dotplot.
            return
        data = getattr(self._annotation_panel, "_loaded_data", None)
        if data is None or not getattr(data, "sequence", ""):
            return
        panel.set_active_sequence(data.sequence_id, data.sequence)

    def _ensure_fasta_text_loaded(self) -> None:
        """Fill the Text editor if it was deferred while Annotation was shown."""
        if not self._fasta_text_deferred:
            return
        self._fasta_text_deferred = False
        path = self.path
        try:
            text, truncated = _read_text_preview(path, _TEXT_PREVIEW_MAX_CHARS)
        except OSError as exc:
            text = f"Could not read file:\n{path}\n{exc}"
            truncated = False
        self._text_truncated = truncated
        self._exit_text_edit_mode(clear_dirty=True)
        self._fasta_text.setPlainText(text)
        cursor = self._fasta_text.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.Start)
        self._fasta_text.setTextCursor(cursor)

    def _open_text_tab(self) -> None:
        self._ensure_fasta_text_loaded()
        self._refresh_text_view()
        selection = self._current_sequence_selection()
        if selection is not None:
            self._highlight_fasta_range(int(selection[0]), int(selection[1]))

    def _capture_text_selection_into_pending(self) -> None:
        """If Text has a sequence selection, store it for Annotation / Dotplot."""
        if not self._single_fasta_mode:
            return
        cursor = self._fasta_text.textCursor()
        if not cursor.hasSelection():
            return
        mapped = self._fasta_document_range_to_bp(
            self._fasta_text.toPlainText(),
            cursor.selectionStart(),
            cursor.selectionEnd(),
        )
        if mapped is not None:
            self._pending_annotation_selection = mapped

    def _current_sequence_selection(self) -> tuple[int, int] | None:
        if self._pending_annotation_selection is not None:
            return self._pending_annotation_selection
        if self._single_fasta_mode:
            cursor = self._fasta_text.textCursor()
            if cursor.hasSelection():
                mapped = self._fasta_document_range_to_bp(
                    self._fasta_text.toPlainText(),
                    cursor.selectionStart(),
                    cursor.selectionEnd(),
                )
                if mapped is not None:
                    return mapped
        if self._annotation_panel is not None:
            try:
                return self._annotation_panel.selection_range()
            except Exception:  # noqa: BLE001
                logger.exception("Could not read Annotation selection")
        if self._dotplot_panel is not None:
            coords = getattr(self._dotplot_panel, "_selection_xy", None)
            if coords is not None:
                x1, x2, _y1, _y2 = coords
                return (int(x1), int(x2)) if x1 <= x2 else (int(x2), int(x1))
        return None

    @staticmethod
    def _fasta_bp_to_document_index(text: str, bp: int) -> int:
        """Map a 1-based sequence position to a character index in FASTA text."""
        if bp < 1:
            return 0
        i = 0
        if text.startswith(">"):
            nl = text.find("\n")
            i = nl + 1 if nl >= 0 else 0
        seen = 0
        while i < len(text):
            ch = text[i]
            if ch in " \t\r\n":
                i += 1
                continue
            seen += 1
            if seen >= bp:
                return i
            i += 1
        return len(text)

    @staticmethod
    def _fasta_document_range_to_bp(
        text: str, start_idx: int, end_idx: int
    ) -> tuple[int, int] | None:
        """Map a document selection to a 1-based inclusive sequence range."""
        lo, hi = (start_idx, end_idx) if start_idx <= end_idx else (end_idx, start_idx)
        if lo >= hi or not text:
            return None
        i = 0
        if text.startswith(">"):
            nl = text.find("\n")
            i = nl + 1 if nl >= 0 else 0
        first_bp: int | None = None
        last_bp: int | None = None
        seen = 0
        while i < len(text):
            ch = text[i]
            if ch in " \t\r\n":
                i += 1
                continue
            seen += 1
            if lo <= i < hi:
                if first_bp is None:
                    first_bp = seen
                last_bp = seen
            elif i >= hi:
                break
            i += 1
        if first_bp is None or last_bp is None:
            return None
        return (first_bp, last_bp)

    def _highlight_fasta_range(self, start: int, end: int) -> None:
        """Select [start, end] (1-based inclusive) in the Text preview."""
        editor = self._fasta_text
        plain = editor.toPlainText()
        if not plain:
            return
        lo, hi = (start, end) if start <= end else (end, start)
        begin = self._fasta_bp_to_document_index(plain, lo)
        # End is inclusive: advance to the character after the last selected base.
        end_idx = self._fasta_bp_to_document_index(plain, hi)
        if end_idx < len(plain) and plain[end_idx] not in " \t\r\n":
            end_idx += 1
        cursor = editor.textCursor()
        cursor.setPosition(begin)
        cursor.setPosition(end_idx, QTextCursor.MoveMode.KeepAnchor)
        editor.setTextCursor(cursor)
        editor.ensureCursorVisible()

    def _on_dotplot_range_selected(
        self, x_start: int, x_end: int, y_start: int, y_end: int
    ) -> None:
        """Mirror a dotplot rectangle onto the Annotation view's selection.

        The x range is used: it is the region under the cursor horizontally,
        which is what "select this part of the sequence" means. The tab is not
        switched -- you stay on the plot and move across when you choose to.
        """
        del y_start, y_end
        self._pending_annotation_selection = (int(x_start), int(x_end))
        # If Annotation happens to be the visible tab already, apply it now;
        # otherwise it is applied when that tab is opened.
        if self._fasta_stack.currentIndex() == 1:
            self._apply_pending_annotation_selection()
        elif self._fasta_stack.currentIndex() == 0:
            self._highlight_fasta_range(int(x_start), int(x_end))

    def _on_dotplot_view_in_annotation(self, start: int, end: int) -> None:
        """Switch to Annotation and highlight the chosen x or y range."""
        self._pending_annotation_selection = (int(start), int(end))
        self._activate_fasta_tab(DefaultFastaPreviewTab.ANNOTATION)

    def _apply_pending_annotation_selection(self) -> None:
        pending = self._pending_annotation_selection
        if pending is None:
            return
        try:
            if self._multi_fasta_mode:
                overview = self._ensure_annotation_overview()
                if not overview.record_ids():
                    overview.load_path(self.path, msa_mode=self._msa_aligned_mode)
                seq_id = None
                if self._dotplot_panel is not None:
                    seq_id = self._dotplot_panel.sequence_id() or None
                overview.select_by_position(
                    int(pending[0]), int(pending[1]), sequence_id=seq_id
                )
                # Layout may settle after the tab switch; re-center the row.
                if seq_id:
                    row = overview._row_for_id(seq_id)
                    if row >= 0:
                        QTimer.singleShot(
                            50,
                            lambda r=row: overview._center_selection_vertically(
                                rows={r}
                            ),
                        )
            else:
                panel = self._ensure_annotation_panel()
                # Load only if this panel has no sequence yet: load_path() re-parses
                # the whole FASTA, so calling it per selection would be very slow.
                if getattr(panel, "_loaded_data", None) is None:
                    panel.load_path(self.path)
                panel.select_by_position(int(pending[0]), int(pending[1]))
        except Exception:  # noqa: BLE001
            # Do not swallow silently: a missing method here previously made
            # the whole dotplot/annotation sync a no-op with no visible error.
            logger.exception("Could not apply dotplot selection to Annotation")
            return
        self._pending_annotation_selection = None

    def bind_annotation_header_zoom(self, bar) -> None:
        self._ensure_annotation_panel().bind_annotation_zoom_bar(bar)

    def _preview_tab_host(self):
        tab_host = self.parent()
        while tab_host is not None:
            from repeatviewer.widgets.preview_tabs import PreviewTabWidget

            if isinstance(tab_host, PreviewTabWidget):
                return tab_host
            tab_host = tab_host.parent()
        return None

    def _activate_fasta_tab(self, tab: DefaultFastaPreviewTab) -> None:
        if tab == DefaultFastaPreviewTab.DOTPLOT:
            self.set_fasta_tab_index(2)
            self._refresh_sequence_info_bar()
            QTimer.singleShot(0, self._open_dotplot_tab)
            return
        if tab == DefaultFastaPreviewTab.ANNOTATION:
            if self._multi_fasta_mode:
                self._show_annotation_content(multi=True)
            else:
                self._ensure_annotation_panel()
                self._show_annotation_content(multi=False)
            self.set_fasta_tab_index(1)
            self._refresh_sequence_info_bar()
            QTimer.singleShot(0, self._open_annotation_tab)
            return
        self.set_fasta_tab_index(0)
        self._refresh_sequence_info_bar()
        QTimer.singleShot(0, self._open_text_tab)

    def _open_annotation_tab(self) -> None:
        self._capture_text_selection_into_pending()
        if self._multi_fasta_mode:
            self._show_annotation_content(multi=True)
            if self._fasta_stack.currentIndex() != 1:
                self.set_fasta_tab_index(1)
                return
            self._ensure_annotation_data_loaded(
                progress_detail="Loading annotations…"
            )
            self._refresh_sequence_info_bar()
            if self._pending_annotation_selection is not None:
                QTimer.singleShot(0, self._apply_pending_annotation_selection)
            return
        panel = self._ensure_annotation_panel()
        self._show_annotation_content(multi=False)
        tab_host = self._preview_tab_host()
        if tab_host is not None and tab_host._annotation_zoom_bar is not None:
            panel.bind_annotation_zoom_bar(tab_host._annotation_zoom_bar)
        if self._fasta_stack.currentIndex() != 1:
            self.set_fasta_tab_index(1)
            return
        self._ensure_annotation_data_loaded(progress_detail="Loading annotations…")
        self._refresh_sequence_info_bar()
        panel.open_annotation_tab()
        # Apply after open_annotation_tab()'s own deferred layout work has run,
        # otherwise its pending fit can land on top of the selection.
        if self._pending_annotation_selection is not None:
            QTimer.singleShot(0, self._apply_pending_annotation_selection)

    def contains_sequence_id(self, seq_id: str) -> bool:
        if not self._single_fasta_mode:
            return False
        from repeatviewer.backends.gui_functions import sanitize_name
        from repeatviewer.backends.gff_features import _seqid_matches
        from repeatviewer.backends.sequence_features import list_fasta_record_ids

        try:
            record_ids = list_fasta_record_ids(self.path)
        except (OSError, ValueError):
            return False
        wanted = str(seq_id or "").strip()
        for record_id in record_ids:
            if (
                record_id == wanted
                or sanitize_name(record_id) == wanted
                or _seqid_matches(record_id, wanted)
            ):
                return True
        return False

    def matches_sequence_id_for_go_to(self, seq_id: str) -> bool:
        """True if a Go-to click on ``seq_id`` should jump this tab's Annotation view."""
        return self._single_fasta_mode and self.contains_sequence_id(seq_id)

    # ------------------------------------------------ Edit > Find / Go To Position

    def supports_find(self) -> bool:
        """Edit > Find and Go To Position work on FASTA previews."""
        return bool(self._single_fasta_mode)

    def _show_annotation_for_find(self):
        """Switch to the Annotation tab with its data loaded; returns overview or panel."""
        if self._multi_fasta_mode:
            view = self._ensure_annotation_overview()
            self._show_annotation_content(multi=True)
        else:
            view = self._ensure_annotation_panel()
            tab_host = self._preview_tab_host()
            if tab_host is not None and tab_host._annotation_zoom_bar is not None:
                view.bind_annotation_zoom_bar(tab_host._annotation_zoom_bar)
        if self._fasta_stack.currentIndex() != 1:
            self.set_fasta_tab_index(1)
            self._notify_preview_host()
        self._ensure_annotation_data_loaded(progress_detail="Loading annotations…")
        if not self._multi_fasta_mode and view._map.data() is None:
            view.load_path(self.path, fit_after_load=False)
        return view

    def select_records(self, record_ids: list[str], *, extract_folder: object = None) -> bool:
        """Show the alignment / sequence view with these records selected and the
        first one scrolled into view. False if none of them is in this file."""
        view = self._show_annotation_for_find()
        if self._multi_fasta_mode:
            if extract_folder:
                view.set_extract_folder(Path(extract_folder))
            return view.select_records(record_ids) > 0
        return bool(record_ids)

    def annotation_features_for_find(self) -> dict | None:
        """Annotations as open in the viewer (unsaved edits included), or None."""
        if self._multi_fasta_mode:
            overview = self._annotation_overview
            if overview is not None and overview.record_ids():
                return overview.features_by_id()
            return None
        panel = self._annotation_panel
        data = getattr(panel, "_loaded_data", None) if panel is not None else None
        if data is None:
            return None
        features = {rid: list(feats) for rid, feats in data.features_by_id.items()}
        if data.sequence_id:
            features[data.sequence_id] = list(panel._map.features())
        return features

    def find_start_position(self) -> tuple[str, int] | None:
        """Where Find Next starts: the Annotation selection (record, position)."""
        if self._multi_fasta_mode:
            overview = self._annotation_overview
            rid = overview.selected_record_id() if overview is not None else None
            return (rid, 1) if rid else None
        panel = self._annotation_panel
        if panel is None or panel._map.data() is None:
            return None
        selection = panel.selection_range()
        return panel.active_sequence_id(), (selection[0] if selection else 1)

    def show_find_item(self, item, hits_by_id: dict[str, list[tuple[int, int]]]) -> bool:
        """Show one Find hit in the Annotation view. False if it cannot be shown."""
        view = self._show_annotation_for_find()
        kind = getattr(item, "kind", "") or ""
        whole = kind == "name"
        ann = [item] if kind == "annotation" else None
        if self._multi_fasta_mode:
            return view.show_find_hit(
                item.seq_id,
                item.start,
                item.end,
                hits_by_id,
                whole_sequence=whole,
                select_annotations=ann,
            )
        return view.show_find_hit(
            item.seq_id,
            item.start,
            item.end,
            hits_by_id.get(item.seq_id, []),
            whole_sequence=whole,
            select_annotations=ann,
        )

    def clear_find_hits(self) -> None:
        if self._annotation_overview is not None:
            self._annotation_overview.clear_find_hits()
        if self._annotation_panel is not None:
            self._annotation_panel.clear_find_hits()

    def go_to_region(self, seq_id: str | None, start: int, end: int) -> str | None:
        """Edit > Go To Position. Returns a message for the user when it cannot."""
        view = self._show_annotation_for_find()
        if self._multi_fasta_mode:
            rid = seq_id or view.selected_record_id()
            if not rid:
                return "This file has several sequences: type the name too, e.g. scaffold_5:1000-2000."
            row = view._row_for_id(rid)
            if row < 0:
                return f'There is no sequence "{rid}" in {self.path.name}.'
            rid = view.record_ids()[row]
            length = view.record_length(rid)
            hits = view._bars._find_hits
        else:
            rid = seq_id or view.active_sequence_id()
            if seq_id and not view.ensure_active_sequence(seq_id):
                return f'There is no sequence "{seq_id}" in {self.path.name}.'
            if view.is_too_long_for_display():
                return "This sequence is too long to show in the Annotation view."
            rid = view.active_sequence_id()
            length = view.active_sequence_length()
            hits = view._map._canvas.find_hits()
        if length and start > length:
            return f"Position {start:,} is past the end of {rid} ({length:,} bp)."
        end = min(end, length) if length else end
        view.show_find_hit(rid, start, end, hits)
        return None

    def save_find_items(self, items: list, name: str) -> int:
        """Save motif hits as annotations. Returns how many were saved."""
        from repeatviewer.backends.sequence_features import (
            REPEATVIEWER_INTERNAL_SOURCE,
            SequenceFeature,
        )

        view = self._show_annotation_for_find()
        by_id: dict[str, list[SequenceFeature]] = {}
        for item in items:
            by_id.setdefault(item.seq_id, []).append(
                SequenceFeature(
                    start=item.start,
                    end=item.end,
                    name=name,
                    feature_type="others",
                    color="#f59e0b",
                    strand=item.strand if item.strand in {"+", "-"} else "+",
                    source=REPEATVIEWER_INTERNAL_SOURCE,
                )
            )
        if not by_id:
            return 0
        if self._multi_fasta_mode:
            return sum(len(v) for v in by_id.values()) if view.add_features_and_save(by_id) else 0
        rid = view.active_sequence_id()
        features = by_id.get(rid) or next(iter(by_id.values()))
        return len(features) if view.add_features_and_save(features) else 0

    def jump_annotation_to_region(self, seq_id: str, start: int, end: int) -> None:
        """Switch to the Annotation sub-tab and jump to [start, end].

        Call only after this tab has been made the outer-current tab (see
        ``PreviewTabWidget.show_annotation_region_for_sequence``) - the sequence map
        widget needs to be visible for its zoom/scroll math to run immediately instead
        of falling back to its own visibility-polling retry loop.
        """
        if self._multi_fasta_mode:
            overview = self._ensure_annotation_overview()
            self._show_annotation_content(multi=True)
            if self._fasta_stack.currentIndex() != 1:
                self.set_fasta_tab_index(1)
            self._ensure_annotation_data_loaded(progress_detail="Loading annotations…")
            # Select matching row in the overview list.
            from repeatviewer.backends.gff_features import _seqid_matches

            wanted = str(seq_id or "").strip()
            for row in range(overview._list.count()):
                rid = overview._list.item(row).data(Qt.ItemDataRole.UserRole)
                if rid == wanted or _seqid_matches(str(rid), wanted):
                    overview._list.setCurrentRow(row)
                    overview._list.scrollToItem(overview._list.item(row))
                    break
            return
        panel = self._ensure_annotation_panel()
        tab_host = self._preview_tab_host()
        if tab_host is not None and tab_host._annotation_zoom_bar is not None:
            panel.bind_annotation_zoom_bar(tab_host._annotation_zoom_bar)
        if self._fasta_stack.currentIndex() != 1:
            self.set_fasta_tab_index(1)
        self._ensure_annotation_data_loaded(progress_detail="Loading annotations…")
        panel.ensure_active_sequence(seq_id)
        panel.jump_to_region(start, end)

    def refresh_nucleotide_colors(self) -> None:
        if self._annotation_panel is not None:
            self._annotation_panel.refresh_nucleotide_colors()

    def refresh_similarity_databases(self) -> None:
        if self._annotation_panel is not None:
            self._annotation_panel.refresh_similarity_databases()

    def _preview_progress_host(self):
        window = self.window()
        begin = getattr(window, "begin_preview_open_progress", None)
        if not callable(begin):
            return None
        return window

    def _begin_open_progress(self, path: Path, *, total: int = 0) -> bool:
        host = self._preview_progress_host()
        if host is None or not _should_show_preview_open_progress(path):
            return False
        began = host.begin_preview_open_progress(path.name, total=total)
        # begin may no-op while a background job owns the bar.
        return bool(began) if began is not None else bool(
            getattr(host, "_preview_open_progress", False)
        )

    def _update_open_progress(self, completed: int, detail: str = "") -> None:
        host = self._preview_progress_host()
        if host is None:
            return
        update = getattr(host, "update_preview_open_progress", None)
        if callable(update):
            update(completed, detail)

    def _end_open_progress(self) -> None:
        host = self._preview_progress_host()
        if host is None:
            return
        end = getattr(host, "end_preview_open_progress", None)
        if callable(end):
            end()

    def _annotation_data_is_current(self) -> bool:
        loaded_for = self._annotation_data_loaded_for
        if loaded_for is None:
            return False
        try:
            return loaded_for == self.path.expanduser().resolve()
        except OSError:
            return loaded_for == self.path

    def _ensure_annotation_data_loaded(self, *, progress_detail: str = "") -> None:
        """Load Annotation/overview data once per path (may take a while)."""
        if self._annotation_data_is_current():
            return
        path = self.path
        host = self._preview_progress_host()
        started_here = False
        if host is not None and not getattr(host, "_preview_open_progress", False):
            if progress_detail or _should_show_preview_open_progress(path):
                began = host.begin_preview_open_progress(path.name, total=0)
                started_here = bool(began) if began is not None else bool(
                    getattr(host, "_preview_open_progress", False)
                )
                if started_here and progress_detail:
                    self._update_open_progress(0, progress_detail)
        try:
            if self._multi_fasta_mode:
                self._ensure_annotation_overview().load_path(
                    path, msa_mode=self._msa_aligned_mode
                )
            else:
                self._ensure_annotation_panel().load_path(path)
            try:
                self._annotation_data_loaded_for = path.expanduser().resolve()
            except OSError:
                self._annotation_data_loaded_for = path
        finally:
            if started_here:
                self._end_open_progress()

    def _load_single_fasta_preview(self, path: Path) -> None:
        if _is_fasta_too_large_for_preview(path):
            self._set_single_fasta_mode(False)
            self._multi_fasta_mode = False
            self._fasta_text_deferred = False
            self._show_text(_FASTA_TOO_LARGE_MESSAGE)
            return
        # Stay on Annotation when the user is already flipping files there.
        stay_on_annotation = (
            self._single_fasta_mode and self._fasta_stack.currentIndex() == 1
        )
        show_progress = self._begin_open_progress(path, total=4)
        try:
            self._set_single_fasta_mode(True)
            self._multi_fasta_mode = _is_multi_fasta_file(path)
            if self._multi_fasta_mode:
                from repeatviewer.backends.sequence_features import is_aligned_fasta

                self._msa_aligned_mode = is_aligned_fasta(path)
            else:
                self._msa_aligned_mode = False
            self._annotation_data_loaded_for = None
            default_tab = load_default_fasta_preview_tab()
            if stay_on_annotation:
                default_tab = DefaultFastaPreviewTab.ANNOTATION
            open_annotation = default_tab == DefaultFastaPreviewTab.ANNOTATION
            # Selections belong to one sequence; never carry them to the next file.
            self._pending_annotation_selection = None
            self._exit_text_edit_mode(clear_dirty=True)
            if open_annotation:
                # Skip Text setPlainText — Annotation does not need it.
                self._fasta_text_deferred = True
                self._text_truncated = False
                if self._fasta_text.toPlainText():
                    self._fasta_text.clear()
            else:
                self._fasta_text_deferred = False
                if show_progress:
                    self._update_open_progress(0, "Reading text…")
                try:
                    text, truncated = _read_text_preview(path, _TEXT_PREVIEW_MAX_CHARS)
                except OSError as exc:
                    text = f"Could not read file:\n{path}\n{exc}"
                    truncated = False
                if show_progress:
                    self._update_open_progress(1, "Filling preview…")
                self._text_truncated = truncated
                self._fasta_text.setPlainText(text)
                cursor = self._fasta_text.textCursor()
                cursor.movePosition(QTextCursor.MoveOperation.Start)
                self._fasta_text.setTextCursor(cursor)
            if self._dotplot_panel is not None:
                # Must reset path/state — clear_selection alone leaves the
                # previous file's "too long" prompt and any active sequence.
                self._dotplot_panel.set_path(path)
            if self._annotation_panel is not None:
                # Drop stale loaded sequence so Dotplot sync cannot reinject
                # the previous (e.g. 50 Mbp) record into a small new file.
                # Switching FASTA also ends any live Transfer search/preview.
                self._annotation_panel.cancel_transfer_annotation()
                self._annotation_panel._loaded_data = None
                try:
                    self._annotation_panel._path = path.expanduser().resolve()
                except OSError:
                    self._annotation_panel._path = path
            if show_progress:
                self._update_open_progress(2, "Preparing tabs…")
            self._show_annotation_content(multi=self._multi_fasta_mode)
            if self._cyto_panel is not None:
                self._cyto_panel.set_path(path)
            # Only pay for Annotation parse when that tab is the open target.
            if open_annotation:
                if show_progress:
                    self._update_open_progress(3, "Loading annotations…")
                self._ensure_annotation_data_loaded()
            elif (
                default_tab == DefaultFastaPreviewTab.DOTPLOT
                and self._multi_fasta_mode
            ):
                # Dotplot needs Overview sequences for the Sequence picker.
                self._ensure_annotation_data_loaded()
            self._activate_fasta_tab(default_tab)
            if show_progress:
                self._update_open_progress(3, "Summarizing…")
            self._refresh_sequence_info_bar()
        finally:
            if show_progress:
                self._end_open_progress()

    def is_selection_preview(self) -> bool:
        """True for the temporary combined view of several selected files."""
        return is_selection_preview_path(self.path)

    def _sync_selection_preview_banner(self) -> None:
        if not self.is_selection_preview():
            self._selection_banner.hide()
            return
        count = len(selection_preview_sources(self.path))
        self._selection_banner_label.setText(
            f"Temporary preview of {count} selected files. "
            "Annotation editing is off, and nothing here is saved to your files."
        )
        self._selection_banner.show()

    def _on_save_selection_as_file(self) -> None:
        sources = [p for p in selection_preview_sources(self.path) if p.is_file()]
        if sources:
            self.saveSelectionAsFileRequested.emit(sources)

    def _load(self) -> None:
        path = self.path
        self._sync_selection_preview_banner()
        if _is_fasta_too_large_for_preview(path):
            self._set_single_fasta_mode(False)
            self._show_text(_FASTA_TOO_LARGE_MESSAGE)
            return
        if _is_sequence_map_fasta(path):
            self._load_single_fasta_preview(path)
            return
        self._set_single_fasta_mode(False)
        suffix = path.suffix.lower()
        if not is_previewable_file(path):
            self._show_text(
                f"Preview is not available for this file type ({suffix or 'unknown'}).\n\n"
                f"{path}\n\n"
                "Open it with an external application (for example Word or Excel)."
            )
            return
        if suffix in _DOTPLOT_PREVIEW_SUFFIXES:
            self._show_alignment_dotplot(path)
            return
        if _is_tree_path(path):
            self._show_tree(path)
            return
        if suffix == _CONTIG_PREVIEW_SUFFIX:
            self._show_contig(path)
            return
        if suffix == ".pdf":
            status = self._pdf.load_file(path)
            if status == QPdfDocument.Error.None_:
                self._stack.setCurrentWidget(self._pdf)
            else:
                self._show_text(f"Could not load PDF:\n{path}\nstatus={status}")
            return
        if suffix in {".html", ".htm"}:
            self._schedule_html_load(path)
            return
        if suffix in {".png", ".jpg", ".jpeg", ".svg", ".gif", ".bmp", ".webp"}:
            if not self._image.set_image_file(path):
                self._show_text(f"Could not load image:\n{path}")
                return
            self._stack.setCurrentWidget(self._image)
            return
        if _is_text_like_preview(path):
            show_progress = self._begin_open_progress(path, total=2)
            try:
                if show_progress:
                    self._update_open_progress(0, "Reading text…")
                try:
                    max_chars = (
                        _TEXT_PREVIEW_MAX_CHARS
                        if _is_text_like_preview(path)
                        else _DEFAULT_TEXT_MAX_CHARS
                    )
                    text, truncated = _read_text_preview(path, max_chars)
                    if show_progress:
                        self._update_open_progress(1, "Filling preview…")
                    self._show_text(text, truncated=truncated)
                except OSError as exc:
                    self._show_text(
                        f"Could not read file:\n{path}\n{exc}", truncated=False
                    )
            finally:
                if show_progress:
                    self._end_open_progress()
            return
        if _looks_like_binary_file(path):
            self._show_text(
                f"Preview is not available for this binary file.\n\n{path}\n\n"
                "Open it with an external application.",
                truncated=False,
            )
            return
        try:
            text, truncated = _read_text_preview(path, _DEFAULT_TEXT_MAX_CHARS)
            self._show_text(text, truncated=truncated)
        except OSError as exc:
            self._show_text(f"Could not read file:\n{path}\n{exc}", truncated=False)

    def _refresh_text_view(self) -> None:
        editor = self._fasta_text if self._single_fasta_mode else self._text
        # NoWrap previews should not force a document text width — that can leave
        # a tiny painted strip at the top of an otherwise blank preview panel.
        if editor.lineWrapMode() == QTextEdit.LineWrapMode.NoWrap:
            editor.document().setTextWidth(-1)
        else:
            editor.document().setTextWidth(max(1, editor.viewport().width()))
        editor.updateGeometry()

    def _active_text_editor(self) -> QTextEdit:
        return self._fasta_text if self._single_fasta_mode else self._text

    def is_html_preview_active(self) -> bool:
        if self._single_fasta_mode:
            return False
        return (
            self._stack.currentWidget() is self._html
            and self.path.suffix.lower() in {".html", ".htm"}
        )

    def suspend_html_preview(self) -> None:
        if hasattr(self._html, "suspend"):
            self._html.suspend()

    def resume_html_preview(self) -> None:
        if not self.is_html_preview_active():
            return
        if hasattr(self._html, "resume"):
            self._html.resume()
        self._reload_html_if_needed()

    def teardown_html_preview(self) -> None:
        if hasattr(self._html, "teardown"):
            self._html.teardown()

    def discard_msa_blast_coverage(self) -> None:
        """Drop MSA Cov plot + ``*.blastcov_RV.json`` when this preview closes."""
        overview = getattr(self, "_annotation_overview", None)
        if overview is None:
            return
        discard = getattr(overview, "discard_blast_coverage_session", None)
        if callable(discard):
            try:
                discard()
            except RuntimeError:
                pass

    def _active_text_editor(self) -> QTextEdit:
        return self._fasta_text if self._single_fasta_mode else self._text

    def can_edit_text(self) -> bool:
        if not self.is_text_preview_active():
            return False
        if self.is_selection_preview():
            return False
        if self.path.name.lower().endswith(".gz"):
            return False
        if not self.path.is_file():
            return False
        current = self._active_text_editor().toPlainText()
        if current.startswith("Could not read file:") or current.startswith(
            "Preview is not available"
        ):
            return False
        if current == _FASTA_TOO_LARGE_MESSAGE:
            return False
        return True

    def can_save_text(self) -> bool:
        return (
            self.can_edit_text()
            and self._text_edit_mode
            and self._text_dirty
            and not self._text_truncated
        )

    def begin_edit_text(self) -> bool:
        if not self.can_edit_text():
            return False
        editor = self._active_text_editor()
        if self._text_truncated or _TRUNCATION_MARKER in editor.toPlainText():
            try:
                text, truncated = _read_text_preview(self.path, _EDITABLE_TEXT_MAX_CHARS)
            except OSError as exc:
                QMessageBox.warning(
                    self,
                    "Edit text",
                    f"Could not load the full file for editing:\n{exc}",
                )
                return False
            if truncated:
                QMessageBox.warning(
                    self,
                    "Edit text",
                    "This file is too large to edit safely in Preview.\n"
                    "Open it in an external editor instead.",
                )
                return False
            self._text_truncated = False
            editor.blockSignals(True)
            editor.setPlainText(text)
            editor.blockSignals(False)
        elif editor.toPlainText().startswith(_EMPTY_FILE_PREFIX):
            editor.blockSignals(True)
            editor.setPlainText("")
            editor.blockSignals(False)
        editor.setReadOnly(False)
        editor.document().setModified(False)
        editor.setFocus(Qt.FocusReason.OtherFocusReason)
        self._text_edit_mode = True
        self._text_dirty = False
        return True

    def save_text(self) -> bool:
        if not self._text_edit_mode:
            if not self.begin_edit_text():
                return False
        if self._text_truncated:
            QMessageBox.warning(
                self,
                "Save text",
                "Cannot save a truncated preview. Use Edit text to load the full file first.",
            )
            return False
        editor = self._active_text_editor()
        try:
            _write_text_file(self.path, editor.toPlainText())
        except OSError as exc:
            QMessageBox.warning(self, "Save text", f"Could not save file:\n{exc}")
            return False
        self._text_dirty = False
        self._text_truncated = False
        editor.document().setModified(False)
        if self._single_fasta_mode:
            self._refresh_sequence_info_bar()
            if self._multi_fasta_mode:
                from repeatviewer.backends.sequence_features import is_aligned_fasta

                self._msa_aligned_mode = is_aligned_fasta(self.path)
                if self._annotation_overview is not None:
                    self._annotation_overview.load_path(
                        self.path, msa_mode=self._msa_aligned_mode
                    )
                    self._annotation_overview.set_annotation_editing_enabled(
                        not self.is_selection_preview()
                    )
                    try:
                        self._annotation_data_loaded_for = self.path.expanduser().resolve()
                    except OSError:
                        self._annotation_data_loaded_for = self.path
            elif self._annotation_panel is not None:
                self._annotation_panel.load_path(self.path)
                try:
                    self._annotation_data_loaded_for = self.path.expanduser().resolve()
                except OSError:
                    self._annotation_data_loaded_for = self.path
        return True

    def _exit_text_edit_mode(self, *, clear_dirty: bool = False) -> None:
        self._text.setReadOnly(True)
        self._fasta_text.setReadOnly(True)
        self._text_edit_mode = False
        if clear_dirty:
            self._text_dirty = False

    def _on_text_editor_changed(self) -> None:
        if not self._text_edit_mode:
            return
        self._text_dirty = True
        self._text_changed_at = time.monotonic()
        # Qt updates the undo/redo counts after this signal: look right after.
        QTimer.singleShot(0, self._note_text_undo)
        tab_host = self.parent()
        while tab_host is not None and not isinstance(tab_host, PreviewTabWidget):
            tab_host = tab_host.parent()
        if isinstance(tab_host, PreviewTabWidget):
            tab_host._sync_text_tool_actions()

    def has_unsaved_text(self) -> bool:
        """Edit text is on and the text differs from what was loaded or last saved.

        Undoing back to the saved text counts as no change (QTextDocument keeps
        its "modified" flag in step with the undo history).
        """
        if not (self._text_edit_mode and self._text_dirty):
            return False
        try:
            return self._active_text_editor().document().isModified()
        except RuntimeError:
            return True

    def _annotation_views(self) -> list:
        views = []
        for view in (self._annotation_panel, self._annotation_overview):
            if view is not None and hasattr(view, "has_unsaved_annotation_changes"):
                views.append(view)
        return views

    def has_unsaved_annotations(self) -> bool:
        """Annotation edits (added / edited / deleted features) not saved yet."""
        for view in self._annotation_views():
            try:
                if view.has_unsaved_annotation_changes():
                    return True
            except RuntimeError:
                continue
        return False

    def has_unsaved_changes(self) -> bool:
        return self.has_unsaved_text() or self.has_unsaved_annotations()

    def save_annotations_now(self) -> bool:
        """Save button behaviour for the leave-file prompt. True when nothing is left unsaved."""
        for view in self._annotation_views():
            try:
                if not view.has_unsaved_annotation_changes():
                    continue
                if hasattr(view, "save_annotations"):
                    view.save_annotations(silent=True)
                else:
                    view._save_overview_annotations()
                if view.has_unsaved_annotation_changes():
                    return False
            except RuntimeError:
                continue
        return True

    def discard_annotation_changes(self) -> None:
        for view in self._annotation_views():
            try:
                view.discard_annotation_changes()
            except RuntimeError:
                continue

    def can_undo_unsaved_text(self) -> bool:
        """Unsaved typing that Edit > Undo can take back."""
        if not self.has_unsaved_text():
            return False
        try:
            return self._active_text_editor().document().isUndoAvailable()
        except RuntimeError:
            return False

    def undo_text(self) -> None:
        editor = self._active_text_editor()
        editor.undo()
        editor.ensureCursorVisible()
        self._text_undone_at = time.monotonic()
        self._text_redo_steps = editor.document().availableRedoSteps()

    def _note_text_undo(self) -> None:
        """Remember when typing was last undone (by any route: menu, key, right-click)."""
        try:
            steps = self._active_text_editor().document().availableRedoSteps()
        except RuntimeError:
            return
        if steps > self._text_redo_steps:
            self._text_undone_at = time.monotonic()
        self._text_redo_steps = steps

    def _show_text_context_menu(self, editor: QTextEdit, pos) -> None:
        from repeatviewer.widgets.paste_reverse_complement import add_to_context_menu

        menu = editor.createStandardContextMenu()
        status = getattr(self.window(), "_set_status_message", None)
        add_to_context_menu(
            menu,
            editor,
            on_message=(lambda text: status(text, 5000)) if callable(status) else None,
        )
        menu.exec(editor.viewport().mapToGlobal(pos))
        menu.deleteLater()

    def can_redo_text(self) -> bool:
        """Undone typing that Edit > Redo can bring back."""
        if not self._text_edit_mode:
            return False
        try:
            return self._active_text_editor().document().isRedoAvailable()
        except RuntimeError:
            return False

    def redo_text(self) -> None:
        editor = self._active_text_editor()
        editor.redo()
        editor.ensureCursorVisible()
        self._text_redo_steps = editor.document().availableRedoSteps()

    def text_changed_at(self) -> float:
        return self._text_changed_at

    def text_undone_at(self) -> float:
        return self._text_undone_at

    def _show_text(self, text: str, *, truncated: bool = False) -> None:
        self._text_truncated = truncated
        self._exit_text_edit_mode(clear_dirty=True)
        self._text.setPlainText(text)
        cursor = self._text.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.Start)
        self._text.setTextCursor(cursor)
        self._stack.setCurrentWidget(self._text_scroll)
        self._text_scroll.show()
        self._text.show()
        QTimer.singleShot(0, self._refresh_text_view)


# Cap per-file tab width so a few long names do not stretch across the title bar.
MAX_PREVIEW_FILE_TAB_WIDTH = 200

# Tabs now persist until closed by hand, so an unbounded number of live
# WebEngine renderers (one OS process each) could pile up. Keep only the N
# most-recently-used pages loaded; older tabs stay in the bar but release
# their HTML preview and reload from disk when clicked again.
MAX_LIVE_PREVIEW_TABS = 8


class PreviewFileTabBar(QTabBar):
    """Title-bar file tabs: content-sized with a max width; scroll when too many."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setObjectName("PreviewFileTabBar")
        self.setElideMode(Qt.TextElideMode.ElideRight)
        self.setExpanding(False)
        self.setDrawBase(False)
        self.setUsesScrollButtons(True)
        self.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        self.setMinimumWidth(0)
        self.setContextMenuPolicy(Qt.ContextMenuPolicy.CustomContextMenu)

    def tabSizeHint(self, index: int) -> QSize:  # noqa: N802
        size = super().tabSizeHint(index)
        if size.width() > MAX_PREVIEW_FILE_TAB_WIDTH:
            size.setWidth(MAX_PREVIEW_FILE_TAB_WIDTH)
        return size

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        # Allow the title-bar layout to shrink this strip; otherwise tabs force
        # themselves under the zoom controls when many files are open.
        return QSize(48, super().minimumSizeHint().height())


class PreviewTabWidget(QTabWidget):
    """Tab bar for open files; reuses tabs when the same path is opened again."""

    annotationsSaved = Signal(object)  # Path to a newly written annotation file
    sequenceExtracted = Signal(object)  # Path to a FASTA extracted from a selection
    openInAdditionalPreviewTab = Signal(object)  # Path to open as a new Preview tab
    dotplotRegionSelected = Signal(dict)  # payload from window.__flexidotLastClick
    fileWritten = Signal(object)  # Path of a file written from a preview (e.g. tree export)
    openFileRequested = Signal(object)  # Path a preview asks to open (contig view → Results)
    blastRequested = Signal(object)  # FASTA a preview sends to BLAST (contig view right-click)
    # Source FASTAs of a temporary multi-file preview, to combine into a real file.
    saveSelectionAsFileRequested = Signal(list)
    # Active preview file path (``Path``) or ``None`` when no tab is open.
    currentPathChanged = Signal(object)

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        zoom_bar: QWidget | None = None,
        header_bar: QWidget | None = None,
        annotation_zoom_bar: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        from repeatviewer.widgets.preview_controls import PreviewHeaderBar, PreviewZoomBar

        self._zoom_bar = zoom_bar if zoom_bar is not None else PreviewZoomBar(self)
        self._header_bar = header_bar
        self._annotation_zoom_bar = annotation_zoom_bar
        if self._annotation_zoom_bar is None and isinstance(self._header_bar, PreviewHeaderBar):
            self._annotation_zoom_bar = self._header_bar._annotation_zoom_bar
        self._zoom_bar.zoom_changed.connect(self._on_preview_zoom_changed)
        self._zoom_bar.fit_mode_changed.connect(self._on_preview_fit_changed)
        if hasattr(self._zoom_bar, "font_size_changed"):
            self._zoom_bar.font_size_changed.connect(self._on_preview_font_size_changed)
        if hasattr(self._zoom_bar, "edit_text_requested"):
            self._zoom_bar.edit_text_requested.connect(self._on_edit_text_requested)
        if hasattr(self._zoom_bar, "save_text_requested"):
            self._zoom_bar.save_text_requested.connect(self._on_save_text_requested)
        if hasattr(self._zoom_bar, "duplicate_window_requested"):
            self._zoom_bar.duplicate_window_requested.connect(self.open_duplicate_window)
        if hasattr(self._zoom_bar, "split_window_requested"):
            self._zoom_bar.split_window_requested.connect(self.open_split_window)
        if self._annotation_zoom_bar is not None and hasattr(
            self._annotation_zoom_bar, "duplicate_window_requested"
        ):
            self._annotation_zoom_bar.duplicate_window_requested.connect(
                self.open_duplicate_window
            )
        if self._annotation_zoom_bar is not None and hasattr(
            self._annotation_zoom_bar, "split_window_requested"
        ):
            self._annotation_zoom_bar.split_window_requested.connect(
                self.open_split_window
            )
        self._save_text_shortcut = QShortcut(QKeySequence.StandardKey.Save, self)
        self._save_text_shortcut.activated.connect(self._on_save_text_requested)
        from repeatviewer.widgets.preview_controls import PreviewWheelZoomFilter

        self._wheel_zoom_filter = PreviewWheelZoomFilter(self, self._zoom_bar)
        if self._header_bar is not None:
            self._header_bar.fasta_tab_changed.connect(self._on_header_fasta_tab_changed)

        self.setTabBar(PreviewFileTabBar(self))
        self.tabBar().customContextMenuRequested.connect(
            lambda pos: self._show_tab_context_menu(self.tabBar(), pos)
        )
        self.setTabsClosable(True)
        self.setMovable(True)
        self.setDocumentMode(True)
        self.setUsesScrollButtons(True)
        self.tabCloseRequested.connect(self._close_tab)
        self._path_to_index: dict[str, int] = {}
        self._lru_order: list[int] = []  # id(FilePreviewPage), most recent first
        self._duplicate_windows: list[QWidget] = []
        self._external_tab_bar_host: QWidget | None = None
        self._title_tab_bar: PreviewFileTabBar | None = None
        self._tab_bar_sync_guard = False
        self.setMinimumHeight(120)

        self.currentChanged.connect(self._on_preview_tab_changed)
        from repeatviewer.widgets.cyto_images_ui import cyto_notifier

        cyto_notifier().changed.connect(self._on_cyto_images_changed)
        self._apply_zoom_to_current_tab()
        self._install_annotation_delete_shortcuts()
        self._sync_tab_bar_visible()
        self._last_emitted_current_path: Path | None = None

    def current_preview_path(self) -> Path | None:
        """Path of the file shown in the active Preview tab, if any."""
        widget = self.currentWidget()
        if isinstance(widget, FilePreviewPage):
            path = getattr(widget, "path", None)
            if isinstance(path, Path):
                return path
        return None

    def _emit_current_path_changed(self) -> None:
        path = self.current_preview_path()
        if path is not None:
            try:
                path = path.expanduser().resolve()
            except OSError:
                path = path.expanduser()
        if path == getattr(self, "_last_emitted_current_path", None):
            return
        self._last_emitted_current_path = path
        self.currentPathChanged.emit(path)

    def persist_session_view_state(self) -> None:
        """Save annotation viewport for the active preview tab (session restore)."""
        widget = self.currentWidget()
        if isinstance(widget, FilePreviewPage):
            widget.persist_annotation_view_for_session()

    def mount_file_tab_bar(self, host: QWidget) -> None:
        """Show open-file tabs in an external host (Preview dock title bar).

        Qt requires QTabWidget's tab bar to stay a child of the tab widget; reparenting
        it into the dock title bar leaves a dangling layout item and segfaults on show().
        Use a mirrored tab bar in the host and keep the real tab bar hidden.
        """
        self._external_tab_bar_host = host
        host.setMinimumWidth(0)
        host.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Preferred)
        self.setCornerWidget(None, Qt.Corner.TopRightCorner)

        internal = self.tabBar()
        internal.hide()

        if self._title_tab_bar is None:
            self._title_tab_bar = PreviewFileTabBar(host)
            host_layout = host.layout()
            if host_layout is not None:
                host_layout.addWidget(self._title_tab_bar)
            self._title_tab_bar.setSizePolicy(
                QSizePolicy.Policy.Expanding,
                QSizePolicy.Policy.Preferred,
            )
            self._title_tab_bar.setMinimumWidth(0)
            self._title_tab_bar.setExpanding(False)
            self._title_tab_bar.setUsesScrollButtons(True)
            self._title_tab_bar.setMovable(True)
            self._title_tab_bar.tabCloseRequested.connect(self._close_tab)
            self._title_tab_bar.currentChanged.connect(self._on_title_tab_current_changed)
            self._title_tab_bar.tabMoved.connect(self._on_title_tab_moved)
            self._title_tab_bar.customContextMenuRequested.connect(
                lambda pos: self._show_tab_context_menu(self._title_tab_bar, pos)
            )

        self._sync_tab_bar_visible()

    def _show_tab_context_menu(self, tab_bar: QTabBar, pos) -> None:
        """Right-click menu on a preview tab: Close / Close Others / Close All.

        Tabs persist until closed by hand, so this is the bulk-cleanup path.
        """
        index = tab_bar.tabAt(pos)
        if index < 0:
            return
        menu = QMenu(tab_bar)
        close_action = menu.addAction("Close Tab")
        others_action = menu.addAction("Close Other Tabs")
        all_action = menu.addAction("Close All Tabs")
        others_action.setEnabled(self.count() > 1)
        chosen = menu.exec(tab_bar.mapToGlobal(pos))
        if chosen is close_action:
            self._close_tab(index)
        elif chosen is others_action:
            self.close_other_tabs(index)
        elif chosen is all_action:
            self.close_all_tabs()

    def _on_title_tab_current_changed(self, index: int) -> None:
        if self._tab_bar_sync_guard or self._title_tab_bar is None:
            return
        if 0 <= index < self.count() and index != self.currentIndex():
            self.setCurrentIndex(index)

    def _on_title_tab_moved(self, from_index: int, to_index: int) -> None:
        if self._tab_bar_sync_guard or self._title_tab_bar is None:
            return
        if not (0 <= from_index < self.count() and 0 <= to_index < self.count()):
            return
        widget = self.widget(from_index)
        if widget is None:
            return
        text = self.tabText(from_index)
        self._tab_bar_sync_guard = True
        try:
            self.removeTab(from_index)
            self.insertTab(to_index, widget, text)
            self.setCurrentIndex(to_index)
            self._rebuild_path_index()
        finally:
            self._tab_bar_sync_guard = False
        self._sync_title_tab_bar()

    def _sync_title_tab_bar(self) -> None:
        title = self._title_tab_bar
        if title is None:
            return
        self._tab_bar_sync_guard = True
        try:
            while title.count() > self.count():
                title.removeTab(title.count() - 1)
            for index in range(title.count(), self.count()):
                title.insertTab(index, self.tabText(index))
            for index in range(self.count()):
                if title.tabText(index) != self.tabText(index):
                    title.setTabText(index, self.tabText(index))
            if title.currentIndex() != self.currentIndex():
                title.setCurrentIndex(self.currentIndex())
            title.setTabsClosable(self.tabsClosable())
        finally:
            self._tab_bar_sync_guard = False

    def _relayout_file_tabs(self) -> None:
        """Refresh title-bar tab sizing after open/close (capped width, no stretch)."""
        tab_bar = self._title_tab_bar if self._title_tab_bar is not None else self.tabBar()
        tab_bar.setExpanding(False)
        tab_bar.setUsesScrollButtons(True)
        tab_bar.setMinimumWidth(0)
        tab_bar.updateGeometry()
        host = self._external_tab_bar_host
        if host is not None:
            host.setMinimumWidth(0)
            host.updateGeometry()

    def _on_header_fasta_tab_changed(self, index: int) -> None:
        widget = self.currentWidget()
        if isinstance(widget, FilePreviewPage):
            widget.set_fasta_tab_index(index)

    def _sync_header_from_current_page(self) -> None:
        if self._header_bar is None:
            return
        widget = self.currentWidget()
        if isinstance(widget, FilePreviewPage) and widget.is_single_fasta_preview():
            # Cytogenetics tab only for sequences that have attached images.
            has_images = widget.has_cyto_images()
            multi = widget.is_multi_fasta_preview()
            if hasattr(self._header_bar, "set_fasta_tab_visible"):
                self._header_bar.set_fasta_tab_visible(2, True)  # Dotplot (Self)
                self._header_bar.set_fasta_tab_visible(3, has_images and not multi)
            if hasattr(self._header_bar, "set_fasta_tab_label"):
                # Same label for single-seq map and multi-seq overview.
                self._header_bar.set_fasta_tab_label(1, "Overview")
            if not has_images and widget.current_fasta_tab_index() == 3:
                widget.set_fasta_tab_index(1)  # last image removed: back to Annotation
                return
            self._header_bar.set_fasta_tabs_visible(True)
            self._header_bar.sync_fasta_tab(widget.current_fasta_tab_index())
        else:
            self._header_bar.set_fasta_tabs_visible(False)

    def _on_cyto_images_changed(self, _folder: str) -> None:
        widget = self.currentWidget()
        if isinstance(widget, FilePreviewPage):
            widget._cyto_images_cache_path = None
            widget._cyto_images_cache_count = None
        self._sync_header_from_current_page()

    def _sync_tab_bar_visible(self) -> None:
        """Show filename tabs only when more than one file is open."""
        multi = self.count() > 1
        if self._title_tab_bar is not None:
            self.tabBar().hide()
            self.setTabsClosable(multi)
            self._sync_title_tab_bar()
            self._title_tab_bar.setVisible(multi)
            self._title_tab_bar.setTabsClosable(multi)
            host = self._external_tab_bar_host
            if host is not None:
                host.setVisible(True)
        else:
            tab_bar = self.tabBar()
            tab_bar.setVisible(multi)
            self.setTabsClosable(multi)
        self._relayout_file_tabs()

    def _install_annotation_delete_shortcuts(self) -> None:
        for key in (QKeySequence.StandardKey.Delete, QKeySequence(Qt.Key.Key_Backspace)):
            shortcut = QShortcut(key, self)
            shortcut.setContext(Qt.ShortcutContext.WidgetWithChildrenShortcut)
            shortcut.activated.connect(self._delete_annotation_features_on_current_tab)

    def _delete_annotation_features_on_current_tab(self) -> None:
        widget = self.currentWidget()
        if isinstance(widget, FilePreviewPage):
            widget._delete_annotation_features_if_active()

    def _touch_lru(self, widget: FilePreviewPage) -> None:
        """Mark ``widget`` as the most recently used page."""
        wid = id(widget)
        try:
            self._lru_order.remove(wid)
        except ValueError:
            pass
        self._lru_order.insert(0, wid)

    def _forget_lru(self, widget: object) -> None:
        try:
            self._lru_order.remove(id(widget))
        except ValueError:
            pass

    def _enforce_live_tab_limit(self) -> None:
        """Release HTML previews on tabs beyond MAX_LIVE_PREVIEW_TABS.

        The tab itself stays in the bar; only its WebEngine renderer is torn
        down. Activating the tab again calls reload_from_disk(), which rebuilds
        the view. This keeps memory flat no matter how many tabs are open.
        """
        pages: dict[int, FilePreviewPage] = {}
        for tab_index in range(self.count()):
            widget = self.widget(tab_index)
            if isinstance(widget, FilePreviewPage):
                pages[id(widget)] = widget
        self._lru_order = [wid for wid in self._lru_order if wid in pages]
        for wid in self._lru_order[MAX_LIVE_PREVIEW_TABS:]:
            page = pages[wid]
            if getattr(page, "_preview_unloaded", False):
                continue
            try:
                page.teardown_html_preview()
            except Exception:  # noqa: BLE001 - never let eviction break tab switching
                logger.debug("Preview eviction failed for %s", page.path, exc_info=True)
                continue
            page._preview_unloaded = True

    def _on_preview_tab_changed(self, index: int) -> None:
        for tab_index in range(self.count()):
            widget = self.widget(tab_index)
            if not isinstance(widget, FilePreviewPage):
                continue
            if tab_index == index:
                self._touch_lru(widget)
                if getattr(widget, "_preview_unloaded", False):
                    widget._preview_unloaded = False
                    widget.reload_from_disk()
                widget.resume_html_preview()
            else:
                widget.suspend_html_preview()
        self._enforce_live_tab_limit()
        if self._title_tab_bar is not None and not self._tab_bar_sync_guard:
            self._tab_bar_sync_guard = True
            try:
                if self._title_tab_bar.currentIndex() != index:
                    self._title_tab_bar.setCurrentIndex(index)
            finally:
                self._tab_bar_sync_guard = False
        self._sync_zoom_bar_enabled(index)
        self._emit_current_path_changed()

    def teardown_all_html_previews(self) -> None:
        for index in range(self.count()):
            widget = self.widget(index)
            if isinstance(widget, FilePreviewPage):
                widget.teardown_html_preview()

    def _wire_preview_page(self, page: FilePreviewPage) -> None:
        page.annotationsSaved.connect(self.annotationsSaved.emit)
        page.sequenceExtracted.connect(self.sequenceExtracted.emit)
        page.openInAdditionalPreviewTab.connect(self.openInAdditionalPreviewTab.emit)
        page.dotplotRegionSelected.connect(self.dotplotRegionSelected.emit)
        page.fileWritten.connect(self.fileWritten.emit)
        page.openFileRequested.connect(self.openFileRequested.emit)
        page.alignmentRequested.connect(self.show_records_in_alignment)
        page.selfDotplotTabRequested.connect(self.open_self_dotplot_tab)
        page._dotplot_in_new_tab = True
        page.blastRequested.connect(self.blastRequested.emit)
        page.saveSelectionAsFileRequested.connect(self.saveSelectionAsFileRequested.emit)
        self._wheel_zoom_filter.register_page(page)

    def show_annotation_region_for_sequence(self, seq_id: str, start: int, end: int) -> bool:
        """Jump the Annotation view of whichever open tab holds ``seq_id`` to [start, end].

        Makes the matching tab current *before* jumping so its sequence map widget is
        already visible (zoom/scroll math needs real widget geometry, not the geometry
        of a still-hidden tab).
        """
        matched = False
        for index in range(self.count()):
            widget = self.widget(index)
            if not isinstance(widget, FilePreviewPage):
                continue
            if not widget.matches_sequence_id_for_go_to(seq_id):
                continue
            matched = True
            self.setCurrentIndex(index)
            self.setCurrentWidget(widget)
            widget.show()
            widget.setFocus(Qt.FocusReason.OtherFocusReason)
            widget.jump_annotation_to_region(seq_id, start, end)
        if matched:
            logger.debug("Go-to: jumped to seq=%r %d-%d", seq_id, start, end)
        return matched

    def _add_preview_page(self, page: FilePreviewPage, *, name: str) -> int:
        self._wire_preview_page(page)
        return self.addTab(page, name)

    @property
    def zoom_bar(self):
        return self._zoom_bar

    def _on_preview_zoom_changed(self, percent: int) -> None:
        self._apply_zoom_to_current_tab(percent=percent)

    def _on_preview_font_size_changed(self, point_size: int) -> None:
        widget = self.currentWidget()
        if isinstance(widget, FilePreviewPage):
            widget.apply_text_font_size(point_size)

    def _on_preview_fit_changed(self, fit_mode: str) -> None:
        self._apply_zoom_to_current_tab(fit_mode=fit_mode)

    def apply_preview_font_size(self, point_size: int | None = None) -> None:
        size = load_preview_font_size() if point_size is None else point_size
        for index in range(self.count()):
            widget = self.widget(index)
            if isinstance(widget, FilePreviewPage):
                widget.apply_text_font_size(size)

    def refresh_nucleotide_color_scheme(self) -> None:
        for index in range(self.count()):
            widget = self.widget(index)
            if isinstance(widget, FilePreviewPage):
                widget.refresh_nucleotide_colors()

    def refresh_similarity_databases(self) -> None:
        for index in range(self.count()):
            widget = self.widget(index)
            if isinstance(widget, FilePreviewPage):
                widget.refresh_similarity_databases()

    def _apply_zoom_to_current_tab(
        self,
        *,
        percent: int | None = None,
        fit_mode: str | None = None,
        sync_chrome: bool = True,
    ) -> None:
        from repeatviewer.widgets.preview_zoom import load_preview_fit_mode

        widget = self.currentWidget()
        if isinstance(widget, FilePreviewPage) and (
            widget.is_text_preview_active()
            or (
                widget._uses_text_preview()
                and not widget.is_annotation_preview_active()
                and not widget.is_dotplot_preview_active()
                and widget.preview_zoom_kind() == PreviewZoomKind.TEXT
            )
        ):
            if hasattr(self._zoom_bar, "current_font_size") and getattr(
                self._zoom_bar, "_font_spin", None
            ) is not None and self._zoom_bar._font_spin.isVisible():
                widget.apply_text_font_size(self._zoom_bar.current_font_size())
            else:
                widget.apply_text_font_size(load_preview_font_size())
            if sync_chrome:
                self._sync_zoom_bar_enabled()
            return

        if percent is None:
            if isinstance(widget, FilePreviewPage) and hasattr(self._zoom_bar, "current_zoom_kind"):
                percent = load_preview_zoom_percent(widget.preview_zoom_kind())
            else:
                percent = load_preview_zoom_percent()
            if self._zoom_bar.isVisible() and self._zoom_bar._spin.isVisible():
                percent = self._zoom_bar.current_zoom_percent()
        if fit_mode is None:
            fit_mode = load_preview_fit_mode().value
            if self._zoom_bar.isVisible():
                fit_mode = self._zoom_bar.current_fit_mode().value
        else:
            fit_mode = str(fit_mode)

        if isinstance(widget, FilePreviewPage):
            widget.apply_preview_zoom(percent, fit_mode)
        if sync_chrome:
            self._sync_zoom_bar_enabled()

    def _sync_zoom_bar_enabled(self, _index: int | None = None) -> None:
        from repeatviewer.widgets.preview_controls import PreviewHeaderBar

        widget = self.currentWidget()
        is_text = False
        is_pdf = False
        is_annotation = False
        is_dotplot = False
        is_cyto = False
        is_media = False
        zoom_kind = None
        if isinstance(widget, FilePreviewPage):
            is_annotation = widget.is_annotation_preview_active()
            is_multi_overview = is_annotation and widget.is_multi_fasta_preview()
            is_dotplot = (
                widget.is_dotplot_preview_active()
                or widget.is_tree_preview_active()
                or widget.is_contig_preview_active()
            )
            is_cyto = widget.is_cyto_preview_active()
            zoom_kind = widget.preview_zoom_kind()
            if is_annotation and not is_multi_overview:
                is_text = False
                is_pdf = False
            elif is_multi_overview:
                # Overview has no nucleotide zoom chrome.
                is_annotation = False
                is_text = False
                is_pdf = False
            elif widget.is_single_fasta_preview():
                is_text = widget.is_text_preview_active()
                is_pdf = False
            else:
                current = widget._stack.currentWidget()
                is_text = current is widget._text_scroll
                is_pdf = current is widget._pdf
                is_media = current in (widget._html, widget._image)
                is_dotplot = is_dotplot or (
                    getattr(widget, "_alignment_dotplot", None) is not None
                    and current is widget._alignment_dotplot
                )

        if isinstance(self._header_bar, PreviewHeaderBar):
            if is_cyto:
                # Cytogenetics: no text/annotation zoom chrome.
                self._header_bar.set_header_zoom_mode("hidden")
            elif (
                isinstance(widget, FilePreviewPage)
                and widget.is_annotation_preview_active()
                and widget.is_multi_fasta_preview()
            ):
                # Multi-seq Overview: −/+/Fit live in the overview strip.
                # Title-bar Actions stay for the main dock via mode "preview".
                # Duplicate windows hide Actions — use "tabs" so no empty gray
                # strip is reserved under Text/Overview.
                if getattr(self._zoom_bar, "_actions_visible", True):
                    self._header_bar.set_header_zoom_mode("preview")
                else:
                    self._header_bar.set_header_zoom_mode("tabs")
                    self._zoom_bar.hide()
            elif is_annotation:
                # Single-seq Annotation: Actions on the Preview title bar (like Dotplot).
                if getattr(self._zoom_bar, "_actions_visible", True):
                    self._header_bar.set_header_zoom_mode("preview")
                else:
                    self._header_bar.set_header_zoom_mode("tabs")
                    self._zoom_bar.hide()
            else:
                self._header_bar.set_header_zoom_mode("preview")
        else:
            self._zoom_bar.setVisible(not is_annotation)

        single_annotation = (
            isinstance(widget, FilePreviewPage)
            and widget.is_annotation_preview_active()
            and not widget.is_multi_fasta_preview()
        )
        if hasattr(self._zoom_bar, "set_chrome_mode"):
            if is_dotplot or single_annotation:
                # Actions only on the dock title (default preview tab).
                self._zoom_bar.set_chrome_mode("dotplot")
            elif is_text or (
                isinstance(widget, FilePreviewPage) and widget._uses_text_preview()
            ):
                self._zoom_bar.set_chrome_mode("text")
            else:
                # HTML, image, and PDF share % / Reset / default-tab Actions.
                self._zoom_bar.set_chrome_mode("media")

        kind_changed = False
        if (
            zoom_kind is not None
            and hasattr(self._zoom_bar, "set_zoom_kind")
            and not is_annotation
            and not is_dotplot
            and not is_cyto
            and not single_annotation
        ):
            previous = self._zoom_bar.current_zoom_kind()
            self._zoom_bar.set_zoom_kind(zoom_kind)
            kind_changed = previous != zoom_kind

        is_multi_fasta = isinstance(widget, FilePreviewPage) and (
            widget.is_multi_fasta_preview()
        )
        is_multi_overview_active = (
            isinstance(widget, FilePreviewPage)
            and widget.is_annotation_preview_active()
            and is_multi_fasta
        )
        show_preview_chrome = not is_cyto and not is_multi_overview_active
        for control in (
            self._zoom_bar._spin,
            getattr(self._zoom_bar, "_font_spin", None),
        ):
            if control is not None:
                control.setEnabled(show_preview_chrome)
        # Keep Actions enabled on multi-seq Overview (and Text) for default tab.
        actions_enabled = (
            show_preview_chrome
            or is_dotplot
            or single_annotation
            or is_multi_fasta
        ) and not is_cyto
        if hasattr(self._zoom_bar, "_actions_btn"):
            self._zoom_bar._actions_btn.setEnabled(actions_enabled)
            # Respect PreviewZoomBar.set_actions_visible(False) in duplicate windows.
            show_actions = getattr(self._zoom_bar, "_actions_visible", True)
            self._zoom_bar._actions_btn.setVisible(bool(show_actions))
        if is_multi_overview_active and hasattr(self._zoom_bar, "set_chrome_mode"):
            # Title-bar Actions only (font/scale spins stay off).
            self._zoom_bar.set_chrome_mode("dotplot")
        if hasattr(self._zoom_bar, "set_default_preview_tab_visible"):
            # FASTA previews: Text / Annotation|Overview [/ Dotplot for single].
            show_default_tab = isinstance(widget, FilePreviewPage) and (
                widget.is_single_fasta_preview() or is_multi_fasta
            )
            self._zoom_bar.set_default_preview_tab_visible(show_default_tab)
        if hasattr(self._zoom_bar, "set_multi_fasta_default_tabs"):
            self._zoom_bar.set_multi_fasta_default_tabs(is_multi_fasta)
        ann_bar = self._annotation_zoom_bar
        if ann_bar is not None:
            try:
                from shiboken6 import isValid

                if not isValid(ann_bar):
                    self._annotation_zoom_bar = None
                    ann_bar = None
            except Exception:
                self._annotation_zoom_bar = None
                ann_bar = None
        if (
            is_multi_overview_active
            and ann_bar is not None
            and hasattr(ann_bar, "set_multi_fasta_default_tabs")
        ):
            ann_bar.set_multi_fasta_default_tabs(True)
        self._sync_text_tool_actions()
        self._sync_header_from_current_page()
        if isinstance(widget, FilePreviewPage):
            self._wheel_zoom_filter.register_page(widget)
        if kind_changed and isinstance(widget, FilePreviewPage):
            if zoom_kind == PreviewZoomKind.TEXT:
                if hasattr(self._zoom_bar, "current_font_size"):
                    widget.apply_text_font_size(self._zoom_bar.current_font_size())
                else:
                    widget.apply_text_font_size(load_preview_font_size())
            else:
                widget.apply_preview_zoom(
                    self._zoom_bar.current_zoom_percent(),
                    self._zoom_bar.current_fit_mode(),
                )

    def _sync_text_tool_actions(self) -> None:
        if not hasattr(self._zoom_bar, "set_text_actions_enabled"):
            return
        widget = self.currentWidget()
        can_edit = False
        can_save = False
        if isinstance(widget, FilePreviewPage):
            can_edit = widget.can_edit_text()
            can_save = widget.can_save_text()
        self._zoom_bar.set_text_actions_enabled(can_edit=can_edit, can_save=can_save)

    def _on_edit_text_requested(self) -> None:
        widget = self.currentWidget()
        if not isinstance(widget, FilePreviewPage):
            return
        if widget.begin_edit_text():
            self._sync_text_tool_actions()

    def _on_save_text_requested(self) -> None:
        widget = self.currentWidget()
        if not isinstance(widget, FilePreviewPage):
            return
        if widget.save_text():
            self._sync_text_tool_actions()
            parent = self.window()
            refresh = getattr(parent, "_refresh_current_folder", None)
            if callable(refresh):
                refresh()

    def _replace_page_at(self, index: int, resolved: Path) -> None:
        """Show ``resolved`` in the existing tab ``index``.

        The page is recreated in place when the preview kind crosses the
        HTML <-> non-HTML boundary (see :func:`_switching_html_preview`);
        otherwise the existing page just loads the new path. Tabs other than
        ``index`` are never touched.
        """
        widget = self.widget(index)
        if not isinstance(widget, FilePreviewPage):
            return
        if str(widget.path) == str(resolved):
            widget.reload_from_disk()
        elif _switching_html_preview(widget.path, resolved):
            was_current = self.currentIndex() == index
            # Unbind before deleteLater — Overview may own the shared zoom bar.
            widget.discard_msa_blast_coverage()
            widget.release_shared_annotation_zoom_bar()
            self._forget_lru(widget)
            widget.teardown_html_preview()
            self.removeTab(index)
            widget.deleteLater()
            page = FilePreviewPage(resolved)
            self._wire_preview_page(page)
            self.insertTab(index, page, resolved.name)
            if was_current:
                self.setCurrentIndex(index)
                self.setCurrentWidget(page)
                page.show()
        else:
            widget.load_path(resolved)
        self.setTabText(index, resolved.name)
        self._rebuild_path_index()

    def unsaved_pages(self) -> list["FilePreviewPage"]:
        pages = []
        for index in range(self.count()):
            widget = self.widget(index)
            if isinstance(widget, FilePreviewPage) and widget.has_unsaved_changes():
                pages.append(widget)
        return pages

    def confirm_leave_page(self, page: object) -> bool:
        """Unsaved Edit-text or annotation changes: ask Save / Don't Save / Cancel.

        Returns True when it is fine to leave the page (saved, or the user chose
        not to save), False to stay (Cancel, or saving failed).
        """
        if not isinstance(page, FilePreviewPage):
            return True
        text_dirty = page.has_unsaved_text()
        annot_dirty = page.has_unsaved_annotations()
        if not (text_dirty or annot_dirty):
            return True
        index = self.indexOf(page)
        if index >= 0 and index != self.currentIndex():
            self.setCurrentIndex(index)
        box = QMessageBox(self.window())
        box.setIcon(QMessageBox.Icon.Warning)
        box.setWindowTitle("Unsaved Changes")
        box.setText(f'Save changes to "{page.path.name}"?')
        box.setStandardButtons(
            QMessageBox.StandardButton.Save
            | QMessageBox.StandardButton.Discard
            | QMessageBox.StandardButton.Cancel
        )
        box.setDefaultButton(QMessageBox.StandardButton.Save)
        box.setEscapeButton(QMessageBox.StandardButton.Cancel)
        box.exec()
        clicked = box.clickedButton()
        answer = (
            box.standardButton(clicked)
            if clicked is not None
            else QMessageBox.StandardButton.Cancel
        )
        if answer == QMessageBox.StandardButton.Save:
            if text_dirty and not page.save_text():
                return False
            if annot_dirty and not page.save_annotations_now():
                return False
            self._sync_text_tool_actions()
            return True
        if answer == QMessageBox.StandardButton.Discard:
            if text_dirty:
                page._exit_text_edit_mode(clear_dirty=True)
            if annot_dirty:
                page.discard_annotation_changes()
            self._sync_text_tool_actions()
            return True
        return False

    def confirm_close_pages(self, pages: list["FilePreviewPage"]) -> bool:
        """Ask for each page with unsaved edits; False as soon as one says Cancel."""
        for page in pages:
            if not self.confirm_leave_page(page):
                return False
        return True

    def open_duplicate_window(self) -> bool:
        """Open a floating Preview window for the current file (Actions menu)."""
        widget = self.currentWidget()
        if not isinstance(widget, FilePreviewPage):
            return False
        path = widget.path.expanduser().resolve()
        if not path.is_file():
            return False
        from repeatviewer.widgets.preview_duplicate_window import PreviewDuplicateWindow

        # Parent to the main window so app quit closes these; keep a Python
        # reference so WA_DeleteOnClose windows are not garbage-collected early.
        host = self.window()
        win = PreviewDuplicateWindow(widget, source_tabs=self, parent=host)
        return self._show_floating_preview_window(win)

    def open_split_window(self) -> bool:
        """Open a floating left/right split Preview of the current file."""
        widget = self.currentWidget()
        if not isinstance(widget, FilePreviewPage):
            return False
        path = widget.path.expanduser().resolve()
        if not path.is_file():
            return False
        from repeatviewer.widgets.preview_split_window import PreviewSplitWindow

        host = self.window()
        win = PreviewSplitWindow(widget, source_tabs=self, parent=host)
        return self._show_floating_preview_window(win)

    def _show_floating_preview_window(self, win: QWidget) -> bool:
        """Track, show, and raise a Duplicate / Split preview window."""
        self._duplicate_windows.append(win)

        def _forget(*_args, window=win) -> None:
            try:
                self._duplicate_windows.remove(window)
            except ValueError:
                pass

        win.destroyed.connect(_forget)
        win.show()
        win.raise_()
        win.activateWindow()
        return True

    def open_file(self, path: Path, *, replace: bool = False) -> bool:
        """Open one file. Returns False if the user chose Cancel at the save prompt.

        ``replace=True`` (single click in the file list) reuses the *active*
        tab: whatever tab you are currently on becomes the file you clicked.
        Every other tab keeps its content and stays open.

        ``replace=False`` (right-click "Open in Additional Tab") always appends
        a new tab and makes it active, leaving existing tabs untouched.

        Either way, a file that is already open somewhere just gets focused
        instead of being loaded a second time.
        """
        resolved = path.expanduser().resolve()
        key = str(resolved)

        existing = self._path_to_index.get(key)
        if existing is not None:
            if 0 <= existing < self.count():
                widget = self.widget(existing)
                if isinstance(widget, FilePreviewPage):
                    widget.reload_from_disk()
                self.setCurrentIndex(existing)
                self._apply_zoom_to_current_tab()
                self._sync_tab_bar_visible()
                self._emit_current_path_changed()
                return True
            del self._path_to_index[key]

        if replace and self.count():
            target = self.currentIndex()
            if not (0 <= target < self.count()):
                target = 0
            # Showing another file in this tab would drop unsaved edits.
            if not self.confirm_leave_page(self.widget(target)):
                return False
            self._replace_page_at(target, resolved)
            self.setCurrentIndex(target)
            self._apply_zoom_to_current_tab()
            self._sync_tab_bar_visible()
            # Same tab index: currentChanged may not fire after replace.
            self._emit_current_path_changed()
            return True

        page = FilePreviewPage(resolved)
        index = self._add_preview_page(page, name=resolved.name)
        self._path_to_index[key] = index
        self.setCurrentIndex(index)
        self.setCurrentWidget(page)
        page.show()
        self._apply_zoom_to_current_tab()
        self._sync_tab_bar_visible()
        self._emit_current_path_changed()
        return True

    def open_files_in_tabs(self, paths: list[Path]) -> None:
        """Open many files (e.g. right-click); keeps existing tabs and adds new ones."""
        for path in paths:
            self.open_file(path, replace=False)
        self._sync_tab_bar_visible()

    def open_self_dotplot_tab(self, path: object, sequence_id: str) -> FilePreviewPage | None:
        """A second tab on the same file, showing Dotplot (Self) for one sequence.
        It is an extra view: the file's own tab stays the one Files and reloads use."""
        resolved = Path(path).expanduser().resolve()
        if not resolved.is_file():
            return None
        page = FilePreviewPage(resolved)
        page._extra_view = True
        index = self._add_preview_page(page, name=f"Dotplot: {sequence_id}")
        self.setTabToolTip(index, f"{resolved.name}, Dotplot (Self) of {sequence_id}")
        self.setCurrentIndex(index)
        self.setCurrentWidget(page)
        page.show()
        # The new page queued its usual first view while loading; switch to the
        # dotplot after that has run, or it would land on Overview again.
        QTimer.singleShot(
            0, lambda pg=page, sid=sequence_id: _qobject_alive(pg) and pg.show_self_dotplot(sid)
        )
        self._apply_zoom_to_current_tab()
        self._sync_tab_bar_visible()
        return page

    def show_records_in_alignment(
        self, path: object, record_ids: object, extract_folder: object = None
    ) -> bool:
        """Open ``path`` in a tab of its own (or the tab it is already in), show its
        alignment and select ``record_ids`` (tree → "Show in alignment").
        ``extract_folder``: where that tab's Extract saves (the tree's folder)."""
        resolved = Path(path).expanduser().resolve()
        if not self.open_file(resolved, replace=False):
            return False
        index = self._path_to_index.get(str(resolved))
        if index is None or not (0 <= index < self.count()):
            return False
        widget = self.widget(index)
        if not isinstance(widget, FilePreviewPage):
            return False
        self.setCurrentIndex(index)
        self.setCurrentWidget(widget)
        widget.show()
        return widget.select_records(list(record_ids or []), extract_folder=extract_folder)

    def show_annotation_for_path(self, path: Path) -> None:
        """Open ``path`` in preview (if needed) and switch to its Annotation tab."""
        resolved = path.expanduser().resolve()
        self.open_file(resolved, replace=False)
        key = str(resolved)
        index = self._path_to_index.get(key)
        if index is None or not (0 <= index < self.count()):
            return
        widget = self.widget(index)
        if not isinstance(widget, FilePreviewPage):
            return
        self.setCurrentIndex(index)
        self.setCurrentWidget(widget)
        widget.show()
        widget.show_annotation_tab()

    def reload_paths(self, paths: list[Path]) -> None:
        """Reload open preview tabs whose paths match (e.g. after DPI change)."""
        keys = {str(p.expanduser().resolve()) for p in paths}
        for index in range(self.count()):
            widget = self.widget(index)
            if isinstance(widget, FilePreviewPage) and str(widget.path) in keys:
                widget.reload_from_disk()
        if self.currentWidget() and str(getattr(self.currentWidget(), "path", "")) in keys:
            self._apply_zoom_to_current_tab()

    def close_all_tabs(self) -> None:
        if not self.confirm_close_pages(self.unsaved_pages()):
            return
        while self.count():
            widget = self.widget(0)
            if isinstance(widget, FilePreviewPage):
                widget.discard_msa_blast_coverage()
                widget.release_shared_annotation_zoom_bar()
                widget.teardown_html_preview()
            self.removeTab(0)
            if widget is not None:
                self._forget_lru(widget)
                widget.deleteLater()
        self._path_to_index.clear()
        self._lru_order.clear()
        self._sync_tab_bar_visible()
        self._emit_current_path_changed()

    def close_other_tabs(self, keep_index: int) -> None:
        """Close every tab except ``keep_index``."""
        if not (0 <= keep_index < self.count()):
            return
        keep = self.widget(keep_index)
        others = [page for page in self.unsaved_pages() if page is not keep]
        if not self.confirm_close_pages(others):
            return
        for index in range(self.count() - 1, -1, -1):
            if self.widget(index) is not keep:
                self._close_tab(index, ask=False)

    def close_paths(self, paths: list[Path]) -> None:
        """Close preview tabs whose files were deleted."""
        keys = {str(p.expanduser().resolve()) for p in paths}
        for index in range(self.count() - 1, -1, -1):
            widget = self.widget(index)
            if isinstance(widget, FilePreviewPage) and str(widget.path) in keys:
                # The file was deleted or renamed: nothing left to save to.
                self._close_tab(index, ask=False)

    def close_under_folder(self, folder: Path) -> None:
        """Close preview tabs for files inside a deleted folder."""
        folder = folder.expanduser().resolve()
        for index in range(self.count() - 1, -1, -1):
            widget = self.widget(index)
            if not isinstance(widget, FilePreviewPage):
                continue
            try:
                if widget.path.resolve().is_relative_to(folder):
                    self._close_tab(index, ask=False)
            except ValueError:
                continue

    def _close_tab(self, index: int, *, ask: bool = True) -> None:
        widget = self.widget(index)
        if ask and not self.confirm_leave_page(widget):
            return
        if isinstance(widget, FilePreviewPage):
            widget.persist_annotation_view_for_session()
            key = str(widget.path)
            self._path_to_index.pop(key, None)
            widget.discard_msa_blast_coverage()
            widget.release_shared_annotation_zoom_bar()
            widget.teardown_html_preview()
        self.removeTab(index)
        if widget is not None:
            self._forget_lru(widget)
            widget.deleteLater()
        self._rebuild_path_index()
        self._sync_tab_bar_visible()
        self._emit_current_path_changed()

    def _rebuild_path_index(self) -> None:
        self._path_to_index.clear()
        for index in range(self.count()):
            widget = self.widget(index)
            if isinstance(widget, FilePreviewPage) and not getattr(widget, "_extra_view", False):
                self._path_to_index.setdefault(str(widget.path), index)
