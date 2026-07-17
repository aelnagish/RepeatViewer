"""Tabbed document preview (text, PDF, HTML, images, alignments)."""

from __future__ import annotations

import gzip
import logging
from collections.abc import Callable
from pathlib import Path

from PySide6.QtCore import QEvent, QSize, QSizeF, Qt, QTimer, Signal
from PySide6.QtGui import QFont, QKeyEvent, QKeySequence, QPixmap, QShortcut, QTextCursor
from PySide6.QtPdf import QPdfDocument
from PySide6.QtWidgets import (
    QLabel,
    QScrollArea,
    QSizePolicy,
    QStackedWidget,
    QTabBar,
    QTabWidget,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from teviewer.backends.bismark_cov import is_bismark_cov_file
from teviewer.backends.igv_viewer import IgvPreviewError, build_igv_view_url, is_alignment_preview_file
from teviewer.backends.sequence_features import (
    format_multi_sequence_fasta_summary_text,
    format_single_sequence_summary_text,
    load_single_sequence_fasta,
    summarize_single_sequence_fasta,
)
from teviewer.widgets.html_preview import _qobject_alive, create_html_preview_widget
from teviewer.widgets.sequence_map_viewer import SequenceMapPanel, _is_delete_key
from teviewer.widgets.preview_zoom import (
    MAX_PREVIEW_FONT_SIZE,
    MIN_PREVIEW_FONT_SIZE,
    DefaultFastaPreviewTab,
    PreviewFitMode,
    load_default_fasta_preview_tab,
    load_preview_fit_mode,
    load_preview_font_size,
    load_preview_zoom_percent,
)

logger = logging.getLogger(__name__)

_FASTA_SUFFIXES = {".fa", ".fasta", ".fas"}
_ALIGNMENT_PREVIEW_SUFFIXES = {".bam", ".sam"}
_FASTA_PREVIEW_MAX_BYTES = 1024**3
_FASTA_TOO_LARGE_MESSAGE = "The file size is too big to show"

_genome_path_provider: Callable[[], Path | None] | None = None


def set_genome_path_provider(provider: Callable[[], Path | None] | None) -> None:
    """Return project reference FASTA for igv.js alignment preview."""
    global _genome_path_provider
    _genome_path_provider = provider


def _genome_path_for_preview() -> Path | None:
    if _genome_path_provider is None:
        return None
    try:
        path = _genome_path_provider()
    except Exception:
        return None
    if path is None:
        return None
    try:
        resolved = path.expanduser().resolve()
    except OSError:
        return None
    return resolved if resolved.is_file() else None


def _count_fasta_sequences_fast(path: Path) -> int:
    """Count FASTA records by header lines (no Biopython)."""
    resolved = path.expanduser().resolve()
    if not resolved.is_file():
        return 0
    count = 0
    try:
        if resolved.name.lower().endswith(".gz"):
            with gzip.open(resolved, "rt", encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.startswith(">"):
                        count += 1
        else:
            with resolved.open(encoding="utf-8", errors="replace") as handle:
                for line in handle:
                    if line.startswith(">"):
                        count += 1
    except OSError:
        return 0
    return count


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


def _is_single_fasta_file(path: Path) -> bool:
    if not _is_fasta_file(path):
        return False
    if _is_fasta_too_large_for_preview(path):
        return False
    return _count_fasta_sequences_fast(path) == 1


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
}
_TEXT_PREVIEW_MAX_CHARS = 500_000
_DEFAULT_TEXT_MAX_CHARS = 200_000
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
}


def _preview_kind(path: Path) -> str:
    suffix = path.suffix.lower()
    if suffix in _ALIGNMENT_PREVIEW_SUFFIXES:
        return "igv"
    if suffix in _HTML_PREVIEW_SUFFIXES:
        return "html"
    if suffix in _FASTA_SUFFIXES or path.name.lower().endswith(".fa.gz"):
        if _is_single_fasta_file(path):
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
    html_like = {"html", "igv"}
    # Recreate the page when crossing HTML ↔ non-HTML so FASTA Text/Annotation
    # chrome and WebEngine hosts do not leave a blank leftover panel.
    return (old_kind in html_like) != (new_kind in html_like)


def _is_text_like_preview(path: Path) -> bool:
    suffix = path.suffix.lower()
    name = path.name.lower()
    if is_bismark_cov_file(path):
        return True
    return suffix in _TEXT_PREVIEW_SUFFIXES or name.endswith("blast.txt")


def _read_text_preview(path: Path, max_chars: int) -> str:
    resolved = path.expanduser().resolve()
    if resolved.name.lower().endswith(".gz"):
        with gzip.open(resolved, "rt", encoding="utf-8", errors="replace") as handle:
            text = handle.read(max_chars + 1)
    else:
        with resolved.open(encoding="utf-8", errors="replace") as handle:
            text = handle.read(max_chars + 1)
    if not text.strip():
        return f"(empty file)\n\n{resolved}"
    if len(text) > max_chars:
        return text[:max_chars] + "\n\n... (truncated preview)"
    return text


def is_previewable_file(path: Path) -> bool:
    """Return True when TEviewer can show the file in the Preview panel."""
    resolved = path.expanduser()
    suffix = resolved.suffix.lower()
    name = resolved.name.lower()
    if is_bismark_cov_file(resolved):
        return True
    if suffix in _UNSUPPORTED_PREVIEW_SUFFIXES:
        return False
    if suffix in _ALIGNMENT_PREVIEW_SUFFIXES:
        return True
    if suffix in _HTML_PREVIEW_SUFFIXES | {".pdf"} | _IMAGE_PREVIEW_SUFFIXES:
        return True
    if suffix in _TEXT_PREVIEW_SUFFIXES or name.endswith("blast.txt"):
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
    if suffix in _ALIGNMENT_PREVIEW_SUFFIXES:
        return 2
    if suffix == ".pdf":
        return 3
    if suffix in _IMAGE_PREVIEW_SUFFIXES:
        return 4
    if suffix in _FASTA_SUFFIXES or name.endswith("blast.txt"):
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
    """PNG/JPEG preview with scroll bars (e.g. Flexidot dotplots)."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self._original: QPixmap | None = None
        self._zoom_percent = load_preview_zoom_percent()
        self._fit_mode = PreviewFitMode.FIT_WINDOW
        self._scale_pending = False

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
        pix = QPixmap(str(path))
        if pix.isNull():
            self._original = None
            self._label.clear()
            return False
        self._original = pix
        self.refresh_scale()
        return True

    def reload_image_file(self, path: Path) -> bool:
        """Re-read the file from disk (e.g. after Flexidot DPI post-processing)."""
        return self.set_image_file(path)

    def eventFilter(self, watched, event) -> bool:  # noqa: N802 (Qt API)
        if watched is self._scroll.viewport() and event.type() == QEvent.Type.Resize:
            self.refresh_scale()
        return super().eventFilter(watched, event)

    def resizeEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().resizeEvent(event)
        self.refresh_scale()

    def showEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().showEvent(event)
        if self._scale_pending:
            QTimer.singleShot(0, self.refresh_scale)

    def set_zoom_percent(self, percent: int) -> None:
        self._zoom_percent = max(25, min(200, int(percent)))
        self.refresh_scale()

    def set_fit_mode(self, mode: PreviewFitMode | str) -> None:
        if isinstance(mode, str):
            try:
                mode = PreviewFitMode(mode)
            except ValueError:
                mode = PreviewFitMode.FIT_WINDOW
        self._fit_mode = mode
        self.refresh_scale()

    def refresh_scale(self) -> None:
        if self._original is None or self._original.isNull():
            return
        viewport = self._scroll.viewport().size()
        if viewport.width() < 20 or viewport.height() < 20:
            self._scale_pending = True
            return
        self._scale_pending = False

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
        dpr = max(1.0, self.devicePixelRatioF())
        physical_w = max(1, int(target_w * dpr))
        display = self._original.scaledToWidth(
            physical_w,
            Qt.TransformationMode.SmoothTransformation,
        )
        display.setDevicePixelRatio(dpr)

        self._label.setPixmap(display)
        self._label.setFixedSize(display.size())
        self._scroll.horizontalScrollBar().setValue(0)
        self._scroll.verticalScrollBar().setValue(0)


class PdfPreviewWidget(QWidget):
    """PDF preview: fit-to-width pages at device resolution (ignores plot zoom settings)."""

    _PAGE_GAP_PX = 10

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

        self._doc = QPdfDocument(self)
        self._render_generation = 0

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

        logical_w = max(200, self._scroll.viewport().width())
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

    def load_file(self, path: Path) -> QPdfDocument.Error:
        self._clear_pages()
        resolved = path.expanduser().resolve()
        status = self._doc.load(str(resolved))
        if status == QPdfDocument.Error.None_:
            self._schedule_render()
        return status

    def refresh(self) -> None:
        """Re-render fit-to-width after tab show or resize."""
        self._schedule_render()


class FilePreviewPage(QWidget):
    """Single-file preview page inside a tab."""

    annotationsSaved = Signal(object)  # Path to a newly written annotation file
    dotplotRegionSelected = Signal(dict)  # payload from window.__flexidotLastClick

    def __init__(self, path: Path, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self.path = path.expanduser().resolve()
        self._single_fasta_mode = False
        self._annotation_panel: SequenceMapPanel | None = None
        self._stack = QStackedWidget(self)
        self._stack.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self._placeholder = QLabel("Loading…")
        self._placeholder.setAlignment(Qt.AlignmentFlag.AlignCenter)
        self._base_text_font_pt = load_preview_font_size()
        self._text = QTextEdit()
        self._text.setReadOnly(True)
        self._text.setLineWrapMode(QTextEdit.LineWrapMode.NoWrap)
        self._text.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._text.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._text.setMinimumSize(200, 120)
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
        self._fasta_text_scroll = QScrollArea(self)
        self._fasta_text_scroll.setWidgetResizable(True)
        self._fasta_text_scroll.setFrameShape(QScrollArea.Shape.NoFrame)
        self._fasta_text_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._fasta_text_scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        self._fasta_text_scroll.setWidget(self._fasta_text)

        self.apply_text_font_size(self._base_text_font_pt)
        self._pdf = PdfPreviewWidget(self)
        self._html = create_html_preview_widget(self)
        if hasattr(self._html, "dotplotRegionSelected"):
            self._html.dotplotRegionSelected.connect(self.dotplotRegionSelected.emit)
        self._image = ScaledImagePreview(self)
        self._fasta_stack = QStackedWidget(self)
        self._fasta_annotation_holder = QWidget(self)
        self._fasta_annotation_holder.setSizePolicy(
            QSizePolicy.Policy.Expanding,
            QSizePolicy.Policy.Expanding,
        )
        self._fasta_stack.addWidget(self._fasta_text_scroll)
        self._fasta_stack.addWidget(self._fasta_annotation_holder)
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

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)
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
        if not self.is_annotation_preview_active():
            return
        if self._annotation_panel is None:
            return
        self._annotation_panel._map._sync_feature_selection_from_annot_track()
        self._annotation_panel._delete_feature(confirm_empty=False)

    def keyPressEvent(self, event: QKeyEvent) -> None:  # noqa: N802
        if _is_delete_key(event.key()) and self.is_annotation_preview_active():
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
        if (
            _preview_kind(old_path) in {"html", "igv"}
            and _preview_kind(self.path) not in {"html", "igv"}
            and hasattr(self._html, "teardown")
        ):
            self._html.teardown()
        self._load()
        self._apply_view_mode()
        QTimer.singleShot(0, self._after_load)

    def is_single_fasta_preview(self) -> bool:
        return self._single_fasta_mode

    def is_text_preview_active(self) -> bool:
        if self._single_fasta_mode:
            return self._fasta_stack.currentIndex() == 0
        return self._stack.currentWidget() is self._text_scroll

    def is_annotation_preview_active(self) -> bool:
        return self._single_fasta_mode and self._fasta_stack.currentIndex() == 1

    def _on_fasta_tab_changed(self, index: int) -> None:
        if not self._single_fasta_mode:
            return
        self._refresh_sequence_info_bar()
        if index == 0:
            QTimer.singleShot(0, self._refresh_text_view)
        elif index == 1:
            QTimer.singleShot(0, self._open_annotation_tab)
        tab_host = self.parent()
        while tab_host is not None and not isinstance(tab_host, PreviewTabWidget):
            tab_host = tab_host.parent()
        if isinstance(tab_host, PreviewTabWidget):
            tab_host._sync_zoom_bar_enabled()
            tab_host._sync_header_from_current_page()

    def _after_load(self) -> None:
        if self._single_fasta_mode:
            index = self._fasta_stack.currentIndex()
            if index == 0:
                self._refresh_text_view()
            elif index == 1:
                self._open_annotation_tab()
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
                self.apply_preview_zoom(
                    load_preview_zoom_percent(),
                    load_preview_fit_mode(),
                )
                QTimer.singleShot(0, self._refresh_text_view)
            elif index == 1:
                QTimer.singleShot(0, self._open_annotation_tab)
            return
        current = self._stack.currentWidget()
        if current in (self._html, self._image, self._text_scroll):
            self.apply_preview_zoom(
                load_preview_zoom_percent(),
                load_preview_fit_mode(),
            )
        if current is self._html:
            QTimer.singleShot(0, self._reload_html_if_needed)
        elif current is self._image:
            QTimer.singleShot(0, self._image.refresh_scale)
        elif current is self._text_scroll:
            QTimer.singleShot(0, self._refresh_text_view)
        elif current is self._pdf:
            QTimer.singleShot(0, self._pdf.refresh)

    def _reload_html_if_needed(self) -> None:
        if self.path.suffix.lower() in {".html", ".htm"}:
            self._schedule_html_load(self.path)
            return
        if is_alignment_preview_file(self.path):
            try:
                url = build_igv_view_url(self.path, _genome_path_for_preview())
            except IgvPreviewError:
                return
            self._schedule_igv_load(url)

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

    def _schedule_igv_load(self, url: str) -> None:
        self._stack.setCurrentWidget(self._html)
        page = self

        def _load() -> None:
            if not _qobject_alive(page):
                return
            if page._stack.currentWidget() is not page._html:
                return
            page._html.load_url(url)

        QTimer.singleShot(0, _load)

    def apply_text_font_size(self, point_size: int) -> None:
        self._base_text_font_pt = max(
            MIN_PREVIEW_FONT_SIZE,
            min(MAX_PREVIEW_FONT_SIZE, int(point_size)),
        )
        self._apply_text_zoom(load_preview_zoom_percent())

    def _uses_text_preview(self) -> bool:
        return _is_text_like_preview(self.path)

    def _apply_text_zoom(self, percent: int) -> None:
        scaled = max(
            MIN_PREVIEW_FONT_SIZE,
            min(
                MAX_PREVIEW_FONT_SIZE,
                int(round(self._base_text_font_pt * int(percent) / 100)),
            ),
        )
        mono = QFont("Courier New", scaled)
        if not mono.exactMatch():
            mono = QFont("Menlo", scaled)
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
        if self._uses_text_preview() or self.is_text_preview_active():
            self._apply_text_zoom(percent)

    def reload_from_disk(self) -> None:
        """Refresh the current view after the file on disk changed."""
        suffix = self.path.suffix.lower()
        if suffix in {".png", ".jpg", ".jpeg", ".svg", ".gif", ".bmp", ".webp"}:
            if self._stack.currentWidget() is self._image:
                self._image.reload_image_file(self.path)
            return
        if suffix in {".html", ".htm"}:
            if self._stack.currentWidget() is self._html:
                self._schedule_html_load(self.path, force=True)
            return
        if is_alignment_preview_file(self.path):
            if self._stack.currentWidget() is self._html:
                self._load()
            return
        if suffix == ".pdf":
            if self._stack.currentWidget() is self._pdf:
                self._pdf.load_file(self.path)
            return
        self._load()

    def _set_single_fasta_mode(self, enabled: bool) -> None:
        self._single_fasta_mode = enabled
        self._apply_view_mode()
        self._refresh_sequence_info_bar()

    def _refresh_sequence_info_bar(self) -> None:
        path = self.path
        if (
            self._single_fasta_mode
            and self.is_annotation_preview_active()
        ):
            self._sequence_info_bar.hide()
            if self._annotation_panel is not None:
                self._annotation_panel.refresh_sequence_summary()
            return

        if not _is_fasta_file(path):
            self._sequence_info_bar.clear()
            self._sequence_info_bar.hide()
            return

        if _is_fasta_too_large_for_preview(path):
            self._sequence_info_bar.hide()
            return

        sequence_count = _count_fasta_sequences_fast(path)
        if sequence_count <= 0:
            self._sequence_info_bar.hide()
            return

        if sequence_count == 1:
            try:
                summary = summarize_single_sequence_fasta(path)
            except (OSError, ValueError):
                self._sequence_info_bar.hide()
                return
            self._sequence_info_bar.setText(format_single_sequence_summary_text(summary))
        else:
            self._sequence_info_bar.setText(
                format_multi_sequence_fasta_summary_text(sequence_count)
            )
        self._sequence_info_bar.show()

    def _ensure_annotation_panel(self) -> SequenceMapPanel:
        if self._annotation_panel is None:
            self._annotation_panel = SequenceMapPanel(self.path, self, autoload=False)
            self._annotation_panel.annotationsSaved.connect(self.annotationsSaved.emit)
            holder_layout = self._fasta_annotation_holder.layout()
            if holder_layout is None:
                holder_layout = QVBoxLayout(self._fasta_annotation_holder)
                holder_layout.setContentsMargins(0, 0, 0, 0)
                holder_layout.setSpacing(0)
            holder_layout.addWidget(self._annotation_panel, stretch=1)
        return self._annotation_panel

    def bind_annotation_header_zoom(self, bar) -> None:
        self._ensure_annotation_panel().bind_annotation_zoom_bar(bar)

    def _preview_tab_host(self):
        tab_host = self.parent()
        while tab_host is not None:
            from teviewer.widgets.preview_tabs import PreviewTabWidget

            if isinstance(tab_host, PreviewTabWidget):
                return tab_host
            tab_host = tab_host.parent()
        return None

    def _activate_fasta_tab(self, tab: DefaultFastaPreviewTab) -> None:
        if tab == DefaultFastaPreviewTab.ANNOTATION:
            self._ensure_annotation_panel()
            self.set_fasta_tab_index(1)
            self._refresh_sequence_info_bar()
            QTimer.singleShot(0, self._open_annotation_tab)
            return
        self.set_fasta_tab_index(0)
        self._refresh_sequence_info_bar()
        QTimer.singleShot(0, self._refresh_text_view)

    def _open_annotation_tab(self) -> None:
        panel = self._ensure_annotation_panel()
        tab_host = self._preview_tab_host()
        if tab_host is not None and tab_host._annotation_zoom_bar is not None:
            panel.bind_annotation_zoom_bar(tab_host._annotation_zoom_bar)
        if self._fasta_stack.currentIndex() != 1:
            self.set_fasta_tab_index(1)
            return
        self._refresh_sequence_info_bar()
        panel.open_annotation_tab()

    def contains_sequence_id(self, seq_id: str) -> bool:
        if not self._single_fasta_mode:
            return False
        try:
            record_id, _ = load_single_sequence_fasta(self.path)
        except (OSError, ValueError):
            return False
        from teviewer.backends.gui_functions import sanitize_name

        return record_id == seq_id or sanitize_name(record_id) == seq_id

    def matches_sequence_id_for_go_to(self, seq_id: str) -> bool:
        """True if a Go-to click on ``seq_id`` should jump this tab's Annotation view."""
        return self._single_fasta_mode and self.contains_sequence_id(seq_id)

    def jump_annotation_to_region(self, seq_id: str, start: int, end: int) -> None:
        """Switch to the Annotation sub-tab and jump to [start, end].

        Call only after this tab has been made the outer-current tab (see
        ``PreviewTabWidget.show_annotation_region_for_sequence``) - the sequence map
        widget needs to be visible for its zoom/scroll math to run immediately instead
        of falling back to its own visibility-polling retry loop.
        """
        panel = self._ensure_annotation_panel()
        tab_host = self._preview_tab_host()
        if tab_host is not None and tab_host._annotation_zoom_bar is not None:
            panel.bind_annotation_zoom_bar(tab_host._annotation_zoom_bar)
        if self._fasta_stack.currentIndex() != 1:
            self.set_fasta_tab_index(1)
        panel.jump_to_region(start, end)

    def refresh_nucleotide_colors(self) -> None:
        if self._annotation_panel is not None:
            self._annotation_panel.refresh_nucleotide_colors()

    def refresh_similarity_databases(self) -> None:
        if self._annotation_panel is not None:
            self._annotation_panel.refresh_similarity_databases()

    def _load_single_fasta_preview(self, path: Path) -> None:
        if _is_fasta_too_large_for_preview(path):
            self._set_single_fasta_mode(False)
            self._show_text(_FASTA_TOO_LARGE_MESSAGE)
            return
        self._set_single_fasta_mode(True)
        try:
            max_chars = _TEXT_PREVIEW_MAX_CHARS
            text = _read_text_preview(path, max_chars)
        except OSError as exc:
            text = f"Could not read file:\n{path}\n{exc}"
        self._fasta_text.setPlainText(text)
        cursor = self._fasta_text.textCursor()
        cursor.movePosition(QTextCursor.MoveOperation.Start)
        self._fasta_text.setTextCursor(cursor)
        if self._annotation_panel is not None:
            self._annotation_panel.load_path(path)
        self._activate_fasta_tab(load_default_fasta_preview_tab())
        self._refresh_sequence_info_bar()

    def _load(self) -> None:
        path = self.path
        if _is_fasta_too_large_for_preview(path):
            self._set_single_fasta_mode(False)
            self._show_text(_FASTA_TOO_LARGE_MESSAGE)
            return
        if _is_single_fasta_file(path):
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
        if is_alignment_preview_file(path):
            try:
                url = build_igv_view_url(path, _genome_path_for_preview())
            except IgvPreviewError as exc:
                self._show_text(
                    "Alignment preview (igv.js)\n\n"
                    f"{path}\n\n"
                    f"{exc}"
                )
                return
            self._schedule_igv_load(url)
            return
        if suffix in {".png", ".jpg", ".jpeg", ".svg", ".gif", ".bmp", ".webp"}:
            if not self._image.set_image_file(path):
                self._show_text(f"Could not load image:\n{path}")
                return
            self._stack.setCurrentWidget(self._image)
            return
        if _is_text_like_preview(path):
            try:
                max_chars = (
                    _TEXT_PREVIEW_MAX_CHARS
                    if _is_text_like_preview(path)
                    else _DEFAULT_TEXT_MAX_CHARS
                )
                self._show_text(_read_text_preview(path, max_chars))
            except OSError as exc:
                self._show_text(f"Could not read file:\n{path}\n{exc}")
            return
        if _looks_like_binary_file(path):
            self._show_text(
                f"Preview is not available for this binary file.\n\n{path}\n\n"
                "Open it with an external application."
            )
            return
        try:
            self._show_text(_read_text_preview(path, _DEFAULT_TEXT_MAX_CHARS))
        except OSError as exc:
            self._show_text(f"Could not read file:\n{path}\n{exc}")

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

    def _show_text(self, text: str) -> None:
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
    dotplotRegionSelected = Signal(dict)  # payload from window.__flexidotLastClick

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        zoom_bar: QWidget | None = None,
        header_bar: QWidget | None = None,
        annotation_zoom_bar: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        from teviewer.widgets.preview_controls import PreviewHeaderBar, PreviewZoomBar

        self._zoom_bar = zoom_bar if zoom_bar is not None else PreviewZoomBar(self)
        self._header_bar = header_bar
        self._annotation_zoom_bar = annotation_zoom_bar
        if self._annotation_zoom_bar is None and isinstance(self._header_bar, PreviewHeaderBar):
            self._annotation_zoom_bar = self._header_bar._annotation_zoom_bar
        self._zoom_bar.zoom_changed.connect(self._on_preview_zoom_changed)
        self._zoom_bar.fit_mode_changed.connect(self._on_preview_fit_changed)
        from teviewer.widgets.preview_controls import PreviewWheelZoomFilter

        self._wheel_zoom_filter = PreviewWheelZoomFilter(self, self._zoom_bar)
        if self._header_bar is not None:
            self._header_bar.fasta_tab_changed.connect(self._on_header_fasta_tab_changed)

        self.setTabBar(PreviewFileTabBar(self))
        self.setTabsClosable(True)
        self.setMovable(True)
        self.setDocumentMode(True)
        self.setUsesScrollButtons(True)
        self.tabCloseRequested.connect(self._close_tab)
        self._path_to_index: dict[str, int] = {}
        self._external_tab_bar_host: QWidget | None = None
        self._title_tab_bar: PreviewFileTabBar | None = None
        self._tab_bar_sync_guard = False
        self.setMinimumHeight(120)

        self.currentChanged.connect(self._on_preview_tab_changed)
        self._apply_zoom_to_current_tab()
        self._install_annotation_delete_shortcuts()
        self._sync_tab_bar_visible()

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

        self._sync_tab_bar_visible()

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
            self._header_bar.set_fasta_tabs_visible(True)
            self._header_bar.sync_fasta_tab(widget.current_fasta_tab_index())
        else:
            self._header_bar.set_fasta_tabs_visible(False)

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

    def _on_preview_tab_changed(self, index: int) -> None:
        for tab_index in range(self.count()):
            widget = self.widget(tab_index)
            if not isinstance(widget, FilePreviewPage):
                continue
            if tab_index == index:
                widget.resume_html_preview()
            else:
                widget.suspend_html_preview()
        if self._title_tab_bar is not None and not self._tab_bar_sync_guard:
            self._tab_bar_sync_guard = True
            try:
                if self._title_tab_bar.currentIndex() != index:
                    self._title_tab_bar.setCurrentIndex(index)
            finally:
                self._tab_bar_sync_guard = False
        self._sync_zoom_bar_enabled(index)

    def teardown_all_html_previews(self) -> None:
        for index in range(self.count()):
            widget = self.widget(index)
            if isinstance(widget, FilePreviewPage):
                widget.teardown_html_preview()

    def _wire_preview_page(self, page: FilePreviewPage) -> None:
        page.annotationsSaved.connect(self.annotationsSaved.emit)
        page.dotplotRegionSelected.connect(self.dotplotRegionSelected.emit)
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
    ) -> None:
        from teviewer.widgets.preview_zoom import load_preview_fit_mode

        if percent is None:
            percent = load_preview_zoom_percent()
            if self._zoom_bar.isVisible():
                percent = self._zoom_bar.current_zoom_percent()
        if fit_mode is None:
            fit_mode = load_preview_fit_mode().value
            if self._zoom_bar.isVisible():
                fit_mode = self._zoom_bar.current_fit_mode().value
        else:
            fit_mode = str(fit_mode)

        widget = self.currentWidget()
        if isinstance(widget, FilePreviewPage):
            widget.apply_preview_zoom(percent, fit_mode)
        self._sync_zoom_bar_enabled()

    def _sync_zoom_bar_enabled(self, _index: int | None = None) -> None:
        from teviewer.widgets.preview_controls import PreviewHeaderBar

        widget = self.currentWidget()
        is_text = False
        is_media = False
        is_pdf = False
        is_annotation = False
        if isinstance(widget, FilePreviewPage):
            is_annotation = widget.is_annotation_preview_active()
            if is_annotation:
                is_text = False
                is_media = False
                is_pdf = False
            elif widget.is_single_fasta_preview():
                is_text = widget.is_text_preview_active()
                is_media = False
                is_pdf = False
            else:
                current = widget._stack.currentWidget()
                is_text = current is widget._text_scroll
                is_media = current in (widget._html, widget._image)
                is_pdf = current is widget._pdf
        show_controls = not is_pdf and not is_annotation
        if isinstance(self._header_bar, PreviewHeaderBar):
            # Zoom controls live in the Preview header for Text/HTML/images.
            # For Annotation mode, zoom/position/font live on the embedded map toolbar.
            self._header_bar.set_header_zoom_mode("preview" if show_controls else "hidden")
        else:
            self._zoom_bar.setVisible(show_controls)
        # Fit mode is configured in Preferences; title bar only has % / Reset.
        for control in (
            self._zoom_bar._spin,
            self._zoom_bar._reset_btn,
        ):
            control.setEnabled(show_controls)
        self._sync_header_from_current_page()
        if isinstance(widget, FilePreviewPage):
            self._wheel_zoom_filter.register_page(widget)

    def open_file(self, path: Path, *, replace: bool = False) -> None:
        """Open one file. replace=True uses a single reusable preview tab."""
        resolved = path.expanduser().resolve()
        key = str(resolved)

        if replace:
            if self.count() == 1:
                widget = self.widget(0)
                if isinstance(widget, FilePreviewPage):
                    if str(widget.path) == key:
                        widget.reload_from_disk()
                        self._sync_tab_bar_visible()
                    elif _switching_html_preview(widget.path, resolved):
                        widget.teardown_html_preview()
                        widget.deleteLater()
                        self.removeTab(0)
                        page = FilePreviewPage(resolved)
                        self._add_preview_page(page, name=resolved.name)
                        self._path_to_index = {key: 0}
                        self.setCurrentIndex(0)
                        self.setCurrentWidget(page)
                        page.show()
                        self._apply_zoom_to_current_tab()
                        self._sync_tab_bar_visible()
                        return
                    else:
                        widget.load_path(resolved)
                    self.setTabText(0, resolved.name)
                    self._path_to_index = {key: 0}
                    self.setCurrentIndex(0)
                    self._apply_zoom_to_current_tab()
                    self._sync_tab_bar_visible()
                    return
            self.close_all_tabs()
            page = FilePreviewPage(resolved)
            self._add_preview_page(page, name=resolved.name)
            self._path_to_index = {key: 0}
            self.setCurrentIndex(0)
            self.setCurrentWidget(page)
            page.show()
            self._apply_zoom_to_current_tab()
            self._sync_tab_bar_visible()
            return

        if key in self._path_to_index:
            index = self._path_to_index[key]
            if 0 <= index < self.count():
                widget = self.widget(index)
                if isinstance(widget, FilePreviewPage):
                    widget.reload_from_disk()
                self.setCurrentIndex(index)
                self._apply_zoom_to_current_tab()
                self._sync_tab_bar_visible()
                return
            del self._path_to_index[key]

        page = FilePreviewPage(resolved)
        index = self._add_preview_page(page, name=resolved.name)
        self._path_to_index[key] = index
        self.setCurrentIndex(index)
        self.setCurrentWidget(page)
        page.show()
        self._apply_zoom_to_current_tab()
        self._sync_tab_bar_visible()

    def open_files_in_tabs(self, paths: list[Path]) -> None:
        """Open many files (e.g. right-click); keeps existing tabs and adds new ones."""
        for path in paths:
            self.open_file(path, replace=False)
        self._sync_tab_bar_visible()

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
        while self.count():
            widget = self.widget(0)
            if isinstance(widget, FilePreviewPage):
                widget.teardown_html_preview()
            self.removeTab(0)
            if widget is not None:
                widget.deleteLater()
        self._path_to_index.clear()
        self._sync_tab_bar_visible()

    def close_paths(self, paths: list[Path]) -> None:
        """Close preview tabs whose files were deleted."""
        keys = {str(p.expanduser().resolve()) for p in paths}
        for index in range(self.count() - 1, -1, -1):
            widget = self.widget(index)
            if isinstance(widget, FilePreviewPage) and str(widget.path) in keys:
                self._close_tab(index)

    def close_under_folder(self, folder: Path) -> None:
        """Close preview tabs for files inside a deleted folder."""
        folder = folder.expanduser().resolve()
        for index in range(self.count() - 1, -1, -1):
            widget = self.widget(index)
            if not isinstance(widget, FilePreviewPage):
                continue
            try:
                if widget.path.resolve().is_relative_to(folder):
                    self._close_tab(index)
            except ValueError:
                continue

    def _close_tab(self, index: int) -> None:
        widget = self.widget(index)
        if isinstance(widget, FilePreviewPage):
            key = str(widget.path)
            self._path_to_index.pop(key, None)
            widget.teardown_html_preview()
        self.removeTab(index)
        if widget is not None:
            widget.deleteLater()
        self._rebuild_path_index()
        self._sync_tab_bar_visible()

    def _rebuild_path_index(self) -> None:
        self._path_to_index.clear()
        for index in range(self.count()):
            widget = self.widget(index)
            if isinstance(widget, FilePreviewPage):
                self._path_to_index[str(widget.path)] = index
