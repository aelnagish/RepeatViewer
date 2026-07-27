"""Embedded HTML preview (TEAid Plotly reports, genome maps)."""

from __future__ import annotations

import json
import logging
import os
import sys
from pathlib import Path

from PySide6.QtCore import QEvent, Qt, QTimer, QUrl, Signal
from PySide6.QtGui import QDesktopServices
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QMessageBox,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from teviewer.backends.blast_genome_plot import is_genome_map_html_file
from teviewer.backends.curation import is_flexidot_html_file
from teviewer.backends.flexidot_go_to import (
    FLEXIDOT_GO_TO_BRIDGE_JS,
    payload_from_console_message,
)
from teviewer.backends.html_preview_policy import (
    embedded_html_disabled_reason,
    flexidot_browser_fallback_reason,
    may_attempt_embedded_html_preview,
    should_use_webengine_for_html,
    use_embedded_html_preview,
)
from teviewer.backends.igv_viewer import (
    IgvPreviewError,
    add_track_to_igv_session,
    is_igv_view_url,
    parse_igv_session_id,
)
from teviewer.plot_preferences import save_view_state, sidecar_prefs_path, view_state_from_dict
from teviewer.widgets.preview_zoom import (
    PreviewFitMode,
    load_preview_fit_mode,
    load_preview_zoom_percent,
)

logger = logging.getLogger(__name__)

try:
    from shiboken6 import isValid as _qt_is_valid
except ImportError:  # pragma: no cover
    def _qt_is_valid(obj) -> bool:  # type: ignore[misc]
        return obj is not None


def _qobject_alive(obj) -> bool:
    if obj is None:
        return False
    try:
        return _qt_is_valid(obj)
    except RuntimeError:
        return False

# Nominal report size when content dimensions are not yet measured.
TEAID_NOMINAL_WIDTH = 1600
TEAID_NOMINAL_HEIGHT = 900

_CONTENT_SIZE_JS = """
(() => {
  const plot = document.querySelector('.plotly-graph-div');
  if (plot) {
    const r = plot.getBoundingClientRect();
    return [Math.ceil(r.width || plot.scrollWidth || 0), Math.ceil(r.height || plot.scrollHeight || 0)];
  }
  const de = document.documentElement;
  if (!de) {
    return [0, 0];
  }
  return [de.scrollWidth || 0, de.scrollHeight || 0];
})()
"""

# Do not import QWebEngineView at module load — it can block or crash app startup on macOS.


def _local_file_url(path: Path) -> QUrl:
    """Build a file:// URL that survives ``#`` and other special characters in paths."""
    return QUrl.fromLocalFile(str(path.expanduser().resolve()))


def _webengine_safe_html_path(path: Path) -> Path:
    """Return a load path without ``#`` when needed (WebEngine + large TEAid HTML).

    ``QUrl.fromLocalFile`` encodes ``#`` as ``%23``, but some macOS WebEngine builds
    still fail to render large Plotly reports from those URLs. Copying/linking to a
    temp name without ``#`` keeps ``file://`` loading (needed for multi‑MB TEAid
    HTML) instead of stuffing the whole file into ``setHtml``.
    """
    resolved = path.expanduser().resolve()
    if "#" not in resolved.name or not resolved.is_file():
        return resolved
    import hashlib
    import shutil
    import tempfile

    digest = hashlib.sha1(str(resolved).encode("utf-8")).hexdigest()[:12]
    safe_name = resolved.name.replace("#", "__")
    dest = Path(tempfile.gettempdir()) / f"teviewer_html_{digest}_{safe_name}"
    try:
        src_mtime = resolved.stat().st_mtime_ns
        if dest.is_file() and dest.stat().st_mtime_ns == src_mtime:
            return dest
        if dest.exists():
            dest.unlink()
        try:
            dest.hardlink_to(resolved)
        except OSError:
            try:
                dest.symlink_to(resolved)
            except OSError:
                shutil.copy2(resolved, dest)
        try:
            os.utime(dest, ns=(resolved.stat().st_atime_ns, src_mtime))
        except OSError:
            pass
        return dest
    except OSError:
        logger.debug("Could not create hash-safe HTML path for %s", resolved, exc_info=True)
        return resolved


def open_html_in_browser(path: Path) -> None:
    """Open an HTML report in the system browser."""
    QDesktopServices.openUrl(_local_file_url(path.expanduser().resolve()))


def create_html_preview_widget(parent: QWidget | None = None) -> QWidget:
    """Lightweight placeholder; real WebEngine loads only when opening HTML."""
    return LazyHtmlPreview(parent)


def _configure_web_settings(settings) -> None:
    import sys

    from PySide6.QtWebEngineCore import QWebEngineProfile, QWebEngineSettings

    QWebEngineProfile.defaultProfile().setSpellCheckEnabled(False)

    attrs = [
        QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls,
        QWebEngineSettings.WebAttribute.LocalContentCanAccessFileUrls,
        QWebEngineSettings.WebAttribute.JavascriptEnabled,
        QWebEngineSettings.WebAttribute.LocalStorageEnabled,
        QWebEngineSettings.WebAttribute.Accelerated2dCanvasEnabled,
    ]
    if hasattr(QWebEngineSettings.WebAttribute, "WebGLEnabled"):
        attrs.append(QWebEngineSettings.WebAttribute.WebGLEnabled)
    for attribute in attrs:
        settings.setAttribute(attribute, True)


class LazyHtmlPreview(QWidget):
    """Defers Qt WebEngine creation until the first HTML file is opened."""

    dotplotRegionSelected = Signal(dict)  # payload from window.__flexidotLastClick

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self._layout = QVBoxLayout(self)
        self._layout.setContentsMargins(0, 0, 0, 0)
        self._inner: QWidget | None = None
        self._zoom_percent = load_preview_zoom_percent()
        self._fit_mode = load_preview_fit_mode()
        self._status = QLabel("HTML preview loads when you open a .html file.")
        self._status.setWordWrap(True)
        self._layout.addWidget(self._status)

    def set_zoom_percent(self, percent: int) -> None:
        self._zoom_percent = max(25, min(200, int(percent)))
        if self._inner is not None and hasattr(self._inner, "set_zoom_percent"):
            if hasattr(self._inner, "is_view_active") and not self._inner.is_view_active():
                return
            self._inner.set_zoom_percent(self._zoom_percent)

    def set_fit_mode(self, mode: PreviewFitMode | str) -> None:
        if isinstance(mode, str):
            try:
                mode = PreviewFitMode(mode)
            except ValueError:
                mode = PreviewFitMode.FIT_WINDOW
        self._fit_mode = mode
        if self._inner is not None and hasattr(self._inner, "set_fit_mode"):
            if hasattr(self._inner, "is_view_active") and not self._inner.is_view_active():
                return
            self._inner.set_fit_mode(mode)

    def load_html_file(self, path: Path, *, force: bool = False) -> None:
        resolved = path.expanduser().resolve()
        if not force and not should_use_webengine_for_html(resolved):
            self._load_browser_fallback(resolved, auto_open=is_flexidot_html_file(resolved))
            return
        if not may_attempt_embedded_html_preview(force=force):
            self._load_browser_fallback(
                resolved,
                reason=embedded_html_disabled_reason(),
            )
            return
        if isinstance(self._inner, HtmlFallbackPreview):
            self._layout.removeWidget(self._inner)
            self._inner.deleteLater()
            self._inner = None
        if self._inner is None:
            self._inner = self._create_inner_viewer()
            if isinstance(self._inner, HtmlFallbackPreview):
                self._inner.load_html_file(resolved)
                self._status.hide()
                self._layout.addWidget(self._inner, stretch=1)
                return
            self._status.hide()
            self._layout.addWidget(self._inner, stretch=1)

        if hasattr(self._inner, "load_html_file"):
            self._inner.load_html_file(path, force=force)
        if self._inner is not None:
            self._status.hide()
        if hasattr(self._inner, "is_view_active") and not self._inner.is_view_active():
            return
        if hasattr(self._inner, "set_zoom_percent"):
            self._inner.set_zoom_percent(self._zoom_percent)
        if hasattr(self._inner, "set_fit_mode"):
            self._inner.set_fit_mode(self._fit_mode)

    def try_embedded_preview(self, path: Path) -> None:
        """Retry embedded Plotly/HTML preview after browser fallback."""
        from teviewer.backends.webengine_probe import clear_webengine_probe_cache

        clear_webengine_probe_cache()
        if self._inner is not None:
            if hasattr(self._inner, "teardown"):
                self._inner.teardown()
            self._layout.removeWidget(self._inner)
            self._inner.deleteLater()
            self._inner = None
        self.load_html_file(path, force=True)

    def _load_browser_fallback(
        self,
        path: Path,
        *,
        reason: str = "",
        auto_open: bool = False,
    ) -> None:
        if reason == "" and is_flexidot_html_file(path):
            reason = flexidot_browser_fallback_reason(path)
        if self._inner is not None:
            if hasattr(self._inner, "teardown"):
                self._inner.teardown()
            self._layout.removeWidget(self._inner)
            self._inner.deleteLater()
            self._inner = None
        self._inner = HtmlFallbackPreview(self, reason=reason or "HTML preview unavailable.")
        self._status.hide()
        self._layout.addWidget(self._inner, stretch=1)
        self._inner.load_html_file(path)
        self._inner.set_retry_handler(lambda p=path: self.try_embedded_preview(p))
        if auto_open:
            QTimer.singleShot(150, lambda p=path: open_html_in_browser(p))

    def load_url(self, url: str) -> None:
        """Load an http(s) page (for example igv.js served on localhost)."""
        if not use_embedded_html_preview():
            if self._inner is None:
                self._inner = HtmlFallbackPreview(
                    self,
                    reason=embedded_html_disabled_reason(),
                )
                self._status.hide()
                self._layout.addWidget(self._inner, stretch=1)
            if hasattr(self._inner, "load_url"):
                self._inner.load_url(url)
            return
        if self._inner is None:
            self._inner = self._create_inner_viewer()
            self._status.hide()
            self._layout.addWidget(self._inner, stretch=1)
        if hasattr(self._inner, "load_url"):
            self._inner.load_url(url)
        if hasattr(self._inner, "set_zoom_percent"):
            self._inner.set_zoom_percent(self._zoom_percent)
        if hasattr(self._inner, "set_fit_mode"):
            self._inner.set_fit_mode(self._fit_mode)

    def suspend(self) -> None:
        """Stop WebEngine work when this preview is hidden behind another tab."""
        if self._inner is not None and hasattr(self._inner, "suspend"):
            self._inner.suspend()

    def resume(self) -> None:
        """Restart JS polling stopped by :meth:`suspend` (no-op if never suspended)."""
        if self._inner is not None and hasattr(self._inner, "resume"):
            self._inner.resume()

    def teardown(self) -> None:
        """Release WebEngine before showing a non-HTML preview in the same tab."""
        if self._inner is None:
            return
        if hasattr(self._inner, "teardown"):
            self._inner.teardown()
        self._layout.removeWidget(self._inner)
        self._inner.setParent(None)
        self._inner.deleteLater()
        self._inner = None
        self._status.show()

    def _create_inner_viewer(self) -> QWidget:
        try:
            from PySide6.QtWebEngineWidgets import QWebEngineView  # noqa: F401

            viewer = HtmlWebEnginePreview(self)
            viewer.dotplotRegionSelected.connect(self.dotplotRegionSelected.emit)
            from teviewer.backends.webengine_probe import mark_webengine_safe

            mark_webengine_safe()
            return viewer
        except ImportError:
            return HtmlFallbackPreview(
                self,
                reason="PySide6-Addons not installed.\n\npip install 'PySide6-Addons>=6.6'",
            )
        except Exception as exc:  # noqa: BLE001
            logger.exception("WebEngine init failed")
            return HtmlFallbackPreview(self, reason=str(exc))


class HtmlWebEnginePreview(QWidget):
    """In-app viewer for Plotly HTML (TEAid, genome maps)."""

    dotplotRegionSelected = Signal(dict)  # payload from window.__flexidotLastClick

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        from PySide6.QtWebEngineWidgets import QWebEngineView

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(0)

        toolbar = QWidget(self)
        toolbar.setObjectName("htmlPreviewToolbar")
        toolbar.setFixedHeight(22)
        toolbar.setStyleSheet(
            "QWidget#htmlPreviewToolbar {"
            "  background: palette(midlight);"
            "  border-bottom: 1px solid palette(mid);"
            "}"
            "QWidget#htmlPreviewToolbar QPushButton {"
            "  min-height: 18px;"
            "  max-height: 18px;"
            "  padding: 0px 8px;"
            "  margin: 0px;"
            "}"
        )
        toolbar_layout = QHBoxLayout(toolbar)
        toolbar_layout.setContentsMargins(4, 1, 4, 1)
        toolbar_layout.setSpacing(4)
        self._add_igv_track_btn = QPushButton("Add track…")
        self._add_igv_track_btn.setToolTip(
            "Add GFF, GTF, BED, BAM, or SAM track to the igv.js viewer"
        )
        self._add_igv_track_btn.clicked.connect(self._add_igv_track)
        self._add_igv_track_btn.hide()
        toolbar_layout.addWidget(self._add_igv_track_btn)
        toolbar_layout.addStretch()
        open_browser = QPushButton("Open in browser")
        open_browser.clicked.connect(self._open_in_browser)
        toolbar_layout.addWidget(open_browser)
        layout.addWidget(toolbar)

        self._view = QWebEngineView(self)
        self._view.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self._view, stretch=1)

        _configure_web_settings(self._view.settings())
        self._view.loadFinished.connect(self._on_load_finished)
        self._view.page().printRequested.connect(self._save_print_as_pdf)
        self._view.page().profile().downloadRequested.connect(self._handle_download)
        self._connect_js_console_bridge()
        self._view.installEventFilter(self)
        self._current_path: Path | None = None
        self._current_mtime_ns: int | None = None
        self._current_url: str | None = None
        self._last_load_ok = False
        self._used_html_fallback = False
        self._zoom_percent = load_preview_zoom_percent()
        self._fit_mode = load_preview_fit_mode()
        self._content_w = TEAID_NOMINAL_WIDTH
        self._content_h = TEAID_NOMINAL_HEIGHT
        self._prefs_save_timer: QTimer | None = None
        self._dotplot_click_timer: QTimer | None = None
        self._last_dotplot_click_ts: float | None = None
        self._load_debounce_timer: QTimer | None = None
        self._pending_load_path: Path | None = None
        self._pending_load_force = False
        self._teardown_scheduled = False
        self._work_generation = 0

    def _export_start_path(self, suffix: str) -> str:
        """Suggested native export path based on the HTML currently being viewed."""
        if self._current_path is not None:
            return str(self._current_path.with_suffix(suffix))
        return str(Path.home() / f"teviewer_plot{suffix}")

    def _handle_download(self, download) -> None:
        """Route Plotly's in-page SVG download through a native Qt save dialog."""
        try:
            if download.page() is not self._view.page():
                return
        except RuntimeError:
            return
        suggested = download.downloadFileName() or "teviewer_plot.svg"
        start = Path(self._export_start_path(".svg")).with_name(suggested)
        destination, _ = QFileDialog.getSaveFileName(
            self,
            "Save plot as SVG",
            str(start),
            "SVG image (*.svg);;All files (*)",
        )
        if not destination:
            download.cancel()
            return
        output = Path(destination).expanduser()
        if output.suffix.lower() != ".svg":
            output = output.with_suffix(".svg")
        download.setDownloadDirectory(str(output.parent))
        download.setDownloadFileName(output.name)
        download.accept()

    def _save_print_as_pdf(self) -> None:
        """Handle window.print() from plot HTML using Qt WebEngine's PDF writer."""
        destination, _ = QFileDialog.getSaveFileName(
            self,
            "Save plot as PDF",
            self._export_start_path(".pdf"),
            "PDF document (*.pdf)",
        )
        if not destination or not self._view_active():
            return
        output = Path(destination).expanduser()
        if output.suffix.lower() != ".pdf":
            output = output.with_suffix(".pdf")
        self._view.page().printToPdf(str(output))

    def _connect_js_console_bridge(self) -> None:
        """Optional faster Go-to pickup when this Qt build exposes the page console signal."""
        if not self._view_active():
            return
        page = self._view.page()
        signal = getattr(page, "javaScriptConsoleMessage", None)
        if signal is not None and hasattr(signal, "connect"):
            signal.connect(self._on_js_console_message)

    def is_view_active(self) -> bool:
        return _qobject_alive(self._view)

    def _view_active(self) -> bool:
        return _qobject_alive(self._view) and _qobject_alive(self._view.page())

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        # Force a repaint when this tab becomes visible again after being hidden
        # behind another preview tab, in case Chromium's background-tab throttling
        # left it in a stale render state that would make clicks miss their target.
        if self._view_active():
            self._view.update()
            self._defer(0, self._apply_zoom)
            self._defer(100, self._apply_zoom)
        if (
            self._view_active()
            and self._current_path is not None
            and is_flexidot_html_file(self._current_path)
            and self._last_load_ok
        ):
            self._install_dotplot_click_bridge()

    def _on_js_console_message(self, _level, message, _line, _source) -> None:
        payload = payload_from_console_message(str(message))
        if payload is not None:
            self._emit_dotplot_payload(payload)

    def _emit_dotplot_payload(self, payload: dict) -> None:
        ts = payload.get("ts")
        if not isinstance(ts, (int, float)):
            return
        if self._last_dotplot_click_ts is not None and ts <= self._last_dotplot_click_ts:
            return
        self._last_dotplot_click_ts = ts
        logger.debug(
            "Flexidot Go-to bridge: payload x_seq=%r y_seq=%r ts=%s",
            payload.get("x_seq"),
            payload.get("y_seq"),
            ts,
        )
        self.dotplotRegionSelected.emit(payload)

    def _defer(self, delay_ms: int, callback) -> None:
        generation = self._work_generation

        def _run() -> None:
            if generation != self._work_generation or not self._view_active():
                return
            callback()

        QTimer.singleShot(delay_ms, _run)

    def _invalidate_pending_work(self) -> None:
        self._work_generation += 1
        if self._prefs_save_timer is not None:
            self._prefs_save_timer.stop()
        if self._dotplot_click_timer is not None:
            self._dotplot_click_timer.stop()

    def set_zoom_percent(self, percent: int) -> None:
        if not self._view_active():
            return
        self._zoom_percent = max(25, min(200, int(percent)))
        self._apply_zoom()

    def set_fit_mode(self, mode: PreviewFitMode | str) -> None:
        if not self._view_active():
            return
        if isinstance(mode, str):
            try:
                mode = PreviewFitMode(mode)
            except ValueError:
                mode = PreviewFitMode.FIT_WINDOW
        self._fit_mode = mode
        self._apply_zoom()

    def eventFilter(self, watched, event) -> bool:  # noqa: N802 (Qt API)
        if watched is self._view and event.type() == QEvent.Type.Resize:
            if self._view_active():
                self._apply_zoom()
        return super().eventFilter(watched, event)

    def resizeEvent(self, event) -> None:  # noqa: N802 (Qt API)
        super().resizeEvent(event)
        if self._view_active():
            self._apply_zoom()

    def _apply_zoom(self) -> None:
        if not self._view_active():
            return
        try:
            vp = self._view.size()
        except RuntimeError:
            return
        vw = int(vp.width())
        vh = int(vp.height())
        # Dock/tab layout often reports 0×0 before the first real resize. Fitting
        # against that zooms Plotly to ~5% and leaves a blank-looking preview.
        if vw < 40 or vh < 40:
            self._defer(50, self._apply_zoom)
            return
        cw = max(1, self._content_w)
        ch = max(1, self._content_h)

        if self._fit_mode == PreviewFitMode.ORIGINAL:
            base = 1.0
        elif self._fit_mode == PreviewFitMode.FIT_WIDTH:
            base = vw / cw
        else:
            base = min(vw / cw, vh / ch)

        factor = base * (self._zoom_percent / 100.0)
        factor = max(0.15, min(5.0, factor))
        try:
            self._view.setZoomFactor(factor)
        except RuntimeError:
            return

    def _query_content_size(self) -> None:
        if not self._view_active():
            return
        self._view.page().runJavaScript(_CONTENT_SIZE_JS, 0, self._on_content_size)

    def _on_content_size(self, size) -> None:
        if not self._view_active():
            return
        if not size or not isinstance(size, (list, tuple)) or len(size) < 2:
            return
        try:
            w = float(size[0])
            h = float(size[1])
        except (TypeError, ValueError):
            return
        if w > 50 and h > 50:
            self._content_w = int(w)
            self._content_h = int(h)
            self._apply_zoom()

    def load_html_file(self, path: Path, *, force: bool = False) -> None:
        if not self._view_active():
            return
        resolved = path.expanduser().resolve()
        self._pending_load_path = resolved
        self._pending_load_force = force
        if self._load_debounce_timer is None:
            self._load_debounce_timer = QTimer(self)
            self._load_debounce_timer.setSingleShot(True)
            self._load_debounce_timer.setInterval(120)
            self._load_debounce_timer.timeout.connect(self._flush_pending_html_load)
        self._load_debounce_timer.start()

    def _flush_pending_html_load(self) -> None:
        if not self._view_active() or self._pending_load_path is None:
            return
        self._load_html_file_now(self._pending_load_path, force=self._pending_load_force)
        self._pending_load_path = None

    def _load_html_file_now(self, path: Path, *, force: bool = False) -> None:
        if not self._view_active():
            return
        resolved = path.expanduser().resolve()
        try:
            mtime_ns = resolved.stat().st_mtime_ns if resolved.is_file() else None
        except OSError:
            mtime_ns = None
        if (
            not force
            and resolved == self._current_path
            and mtime_ns is not None
            and mtime_ns == self._current_mtime_ns
            and self._last_load_ok
            and not self._used_html_fallback
        ):
            return

        self._current_path = resolved
        self._current_mtime_ns = mtime_ns
        self._current_url = None
        if self._prefs_save_timer is not None:
            self._prefs_save_timer.stop()
        if self._dotplot_click_timer is not None:
            self._dotplot_click_timer.stop()
        self._last_dotplot_click_ts = None
        self._last_load_ok = False
        self._used_html_fallback = False
        self._content_w = TEAID_NOMINAL_WIDTH
        self._content_h = TEAID_NOMINAL_HEIGHT

        if not resolved.is_file():
            try:
                self._view.setHtml(f"<p>File not found: {resolved}</p>")
            except RuntimeError:
                return
            return

        # Prefer file:// (works for multi‑MB TEAid Plotly). For names with ``#``,
        # load via a temp path without ``#`` so WebEngine does not blank the page.
        load_path = _webengine_safe_html_path(resolved)
        try:
            self._view.load(_local_file_url(load_path))
        except RuntimeError:
            return
        self._sync_igv_toolbar()

    def _load_html_via_set_html(self, resolved: Path) -> None:
        # Only for small files — large TEAid reports must use file:// loading.
        try:
            size = resolved.stat().st_size
        except OSError:
            size = 0
        if size > 4_000_000:
            load_path = _webengine_safe_html_path(resolved)
            try:
                self._view.load(_local_file_url(load_path))
            except RuntimeError:
                return
            self._sync_igv_toolbar()
            return
        try:
            html = resolved.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            try:
                self._view.setHtml(f"<p>Could not read HTML file: {exc}</p>")
            except RuntimeError:
                return
            return
        base_url = _local_file_url(resolved.parent)
        try:
            self._view.setHtml(html, base_url)
        except RuntimeError:
            return
        self._sync_igv_toolbar()
        # setHtml() sometimes does not emit a reliable loadFinished on macOS.
        # Mark as OK and apply zoom now so the plot does not stay blank.
        self._last_load_ok = True
        self._apply_zoom()
        self._defer(50, self._query_content_size)
        self._defer(250, self._query_content_size)
        self._defer(600, self._query_content_size)

    def suspend(self) -> None:
        """Pause JS polling while hidden behind another tab.

        Deliberately does not navigate the page away (e.g. to about:blank): doing so
        used to blank the Flexidot dotplot whenever Go-to opened/focused the Annotation
        tab, and forced a full reload (re-parsing the whole inline Plotly JSON) on every
        switch back. Chromium already throttles hidden-tab rendering on its own.
        """
        self._invalidate_pending_work()
        if self._load_debounce_timer is not None:
            self._load_debounce_timer.stop()
        self._pending_load_path = None

    def resume(self) -> None:
        """Restart JS polling stopped by :meth:`suspend` without reloading the page."""
        if not self._view_active() or not self._last_load_ok or self._current_path is None:
            return
        if is_flexidot_html_file(self._current_path):
            self._install_dotplot_click_bridge()
        if is_genome_map_html_file(self._current_path):
            self._install_plot_prefs_bridge()

    def load_url(self, url: str) -> None:
        if not self._view_active():
            return
        self._current_path = None
        self._current_mtime_ns = None
        self._current_url = url
        if self._prefs_save_timer is not None:
            self._prefs_save_timer.stop()
        self._last_load_ok = False
        self._used_html_fallback = False
        self._content_w = TEAID_NOMINAL_WIDTH
        self._content_h = TEAID_NOMINAL_HEIGHT
        try:
            self._view.load(QUrl(url))
        except RuntimeError:
            return
        self._sync_igv_toolbar()

    def _sync_igv_toolbar(self) -> None:
        show = is_igv_view_url(self._current_url)
        self._add_igv_track_btn.setVisible(show)

    def _add_igv_track(self) -> None:
        session_id = parse_igv_session_id(self._current_url or "")
        if not session_id:
            return
        paths, _ = QFileDialog.getOpenFileNames(
            self,
            "Add track to igv.js",
            "",
            "Tracks (*.gff *.gff3 *.gtf *.bed *.bam *.sam);;All files (*)",
        )
        if not paths:
            return
        errors: list[str] = []
        for raw in paths:
            try:
                track = add_track_to_igv_session(session_id, Path(raw))
            except IgvPreviewError as exc:
                errors.append(f"{Path(raw).name}: {exc}")
                continue
            except Exception as exc:  # noqa: BLE001
                errors.append(f"{Path(raw).name}: {exc}")
                continue
            script = (
                "window.teviewerIgvLoadTrack("
                + json.dumps(track)
                + ").catch(function(err){"
                + 'alert("Add track failed:\\n" + err);'
                + "});"
            )
            if self._view_active():
                self._view.page().runJavaScript(script)
        if errors:
            QMessageBox.warning(
                self,
                "Add track",
                "Some tracks could not be added:\n\n" + "\n".join(errors),
            )

    def _on_load_finished(self, ok: bool) -> None:
        if not self._view_active():
            return
        if ok:
            self._last_load_ok = True
            self._install_plot_prefs_bridge()
            self._install_dotplot_click_bridge()
            self._apply_zoom()
            self._defer(50, self._query_content_size)
            self._defer(250, self._query_content_size)
            self._defer(600, self._query_content_size)
        if ok or self._used_html_fallback or self._current_path is None:
            return
        self._used_html_fallback = True
        resolved = self._current_path
        logger.warning(
            "HTML preview failed to load %s via file URL; falling back to setHtml()",
            resolved,
        )
        self._defer(0, lambda p=resolved: self._load_html_via_set_html(p) if self._view_active() else None)

    def _install_plot_prefs_bridge(self) -> None:
        if not self._view_active():
            return
        if self._current_path is None or not is_genome_map_html_file(self._current_path):
            return
        prefs_path = sidecar_prefs_path(self._current_path)
        bridge_js = f"""
        window.teviewerSavePlotPrefs = function(payload) {{
            window.__teviewerPendingPlotPrefs = payload;
        }};
        window.__teviewerPrefsPath = {json.dumps(str(prefs_path))};
        """
        self._view.page().runJavaScript(bridge_js)
        if self._prefs_save_timer is None:
            self._prefs_save_timer = QTimer(self)
            self._prefs_save_timer.setInterval(1500)
            self._prefs_save_timer.timeout.connect(self._poll_plot_prefs_save)
        self._prefs_save_timer.start()

    def _poll_plot_prefs_save(self) -> None:
        if not self._view_active():
            return
        if self._current_path is None or not is_genome_map_html_file(self._current_path):
            return

        def _on_payload(payload) -> None:
            if not self._view_active():
                return
            if not isinstance(payload, dict) or self._current_path is None:
                return
            try:
                save_view_state(self._current_path, view_state_from_dict(payload))
            except OSError as exc:
                logger.warning("Could not save plot preferences: %s", exc)

        self._view.page().runJavaScript(
            "window.teviewerCollectPlotPrefs ? window.teviewerCollectPlotPrefs() : null",
            0,
            _on_payload,
        )

    def _install_dotplot_click_bridge(self) -> None:
        """Poll a Flexidot dotplot HTML for match-line clicks (window.__flexidotLastClick)."""
        if not self._view_active():
            return
        if self._current_path is None or not is_flexidot_html_file(self._current_path):
            return
        logger.debug("Flexidot Go-to bridge: starting poll timer for %s", self._current_path)
        self._view.page().runJavaScript(FLEXIDOT_GO_TO_BRIDGE_JS)
        if self._dotplot_click_timer is None:
            self._dotplot_click_timer = QTimer(self)
            self._dotplot_click_timer.setInterval(200)
            self._dotplot_click_timer.timeout.connect(self._poll_dotplot_click)
        self._dotplot_click_timer.start()

    def _poll_dotplot_click(self) -> None:
        if not self._view_active():
            return
        if self._current_path is None or not is_flexidot_html_file(self._current_path):
            return

        def _on_payload(raw) -> None:
            if not self._view_active() or not isinstance(raw, str) or not raw:
                return
            try:
                payload = json.loads(raw)
            except (TypeError, ValueError):
                logger.debug("Dotplot click payload was not valid JSON: %r", raw)
                return
            if isinstance(payload, dict):
                self._emit_dotplot_payload(payload)

        self._view.page().runJavaScript(
            """
(() => {
  var pending = window.__teviewerFlexidotNotifyJson;
  if (pending) {
    window.__teviewerFlexidotNotifyJson = null;
    return pending;
  }
  return window.__flexidotLastClick ? JSON.stringify(window.__flexidotLastClick) : null;
})()
""",
            0,
            _on_payload,
        )

    def _open_in_browser(self) -> None:
        if self._current_url:
            QDesktopServices.openUrl(QUrl(self._current_url))
            return
        if self._current_path and self._current_path.is_file():
            QDesktopServices.openUrl(_local_file_url(self._current_path))

    def teardown(self) -> None:
        if self._teardown_scheduled:
            return
        self._teardown_scheduled = True
        self._invalidate_pending_work()
        if self._load_debounce_timer is not None:
            self._load_debounce_timer.stop()
        self._pending_load_path = None
        try:
            self._view.loadFinished.disconnect(self._on_load_finished)
        except (RuntimeError, TypeError):
            pass
        try:
            self._view.page().printRequested.disconnect(self._save_print_as_pdf)
        except (RuntimeError, TypeError):
            pass
        try:
            self._view.page().profile().downloadRequested.disconnect(self._handle_download)
        except (RuntimeError, TypeError):
            pass

        def _finish_teardown() -> None:
            if not _qobject_alive(self):
                return
            try:
                if self._view_active():
                    self._view.stop()
                    self._view.setUrl(QUrl("about:blank"))
            except RuntimeError:
                logger.debug("WebEngine teardown stop failed", exc_info=True)
            self._current_path = None
            self._current_mtime_ns = None
            self._last_load_ok = False
            self._teardown_scheduled = False

        QTimer.singleShot(0, _finish_teardown)


class HtmlFallbackPreview(QWidget):
    """Open HTML in the system browser when WebEngine is unavailable."""

    def __init__(
        self,
        parent: QWidget | None = None,
        *,
        reason: str = "",
    ) -> None:
        super().__init__(parent)
        layout = QVBoxLayout(self)
        self._label = QLabel(
            "HTML cannot be shown inside TEviewer.\n\n"
            f"{reason}\n\n"
            "Try:\n"
            "  pip install 'PySide6-Addons>=6.6'\n\n"
            "Or open the report in your browser:"
        )
        self._label.setWordWrap(True)
        self._try_btn = QPushButton("Try in-app preview")
        self._try_btn.clicked.connect(self._try_in_app)
        self._try_btn.setToolTip(
            "Attempt to show this report inside TEviewer (requires Qt WebEngine)."
        )
        self._open_btn = QPushButton("Open in browser")
        self._open_btn.clicked.connect(self._open_in_browser)
        btn_row = QHBoxLayout()
        btn_row.addWidget(self._try_btn)
        btn_row.addWidget(self._open_btn)
        btn_row.addStretch()
        layout.addWidget(self._label)
        layout.addLayout(btn_row)
        layout.addStretch()
        self._current_path: Path | None = None
        self._current_url: str | None = None
        self._retry_handler = None

    def set_retry_handler(self, handler) -> None:
        self._retry_handler = handler

    def _try_in_app(self) -> None:
        if self._retry_handler is not None:
            self._retry_handler()

    def load_html_file(self, path: Path) -> None:
        self._current_path = path.expanduser().resolve()
        self._current_url = None
        self._label.setText(
            "HTML report — open in browser\n\n"
            f"File:\n{self._current_path}"
        )

    def load_url(self, url: str) -> None:
        self._current_url = url
        self._current_path = None
        self._label.setText(
            "Alignment viewer — open in browser\n\n"
            f"URL:\n{url}"
        )

    def set_zoom_percent(self, percent: int) -> None:
        del percent

    def set_fit_mode(self, mode: PreviewFitMode | str) -> None:
        del mode

    def is_view_active(self) -> bool:
        return True

    def _open_in_browser(self) -> None:
        if self._current_url:
            QDesktopServices.openUrl(QUrl(self._current_url))
            return
        if self._current_path and self._current_path.is_file():
            QDesktopServices.openUrl(_local_file_url(self._current_path))
