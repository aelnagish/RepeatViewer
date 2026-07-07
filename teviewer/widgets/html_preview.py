"""Embedded HTML preview (TEAid Plotly reports, genome maps)."""

from __future__ import annotations

import json
import logging
import os
import sys
from pathlib import Path

from PySide6.QtCore import QEvent, Qt, QTimer, QUrl
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

_warmup_widget: QWidget | None = None

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


def use_embedded_html_preview() -> bool:
    """
    Return True when Plotly HTML should load inside Qt WebEngine.

    Set ``TEVIEWER_EMBED_HTML=0`` to open HTML only in the system browser.
    """
    flag = os.environ.get("TEVIEWER_EMBED_HTML", "").strip().lower()
    if flag in {"0", "false", "no", "off"}:
        return False
    return True


def create_html_preview_widget(parent: QWidget | None = None) -> QWidget:
    """Lightweight placeholder; real WebEngine loads only when opening HTML."""
    return LazyHtmlPreview(parent)


def warm_html_webengine(parent: QWidget | None = None) -> bool:
    """
    Start Chromium during idle time so the first HTML preview does not freeze the UI.

    Safe to call repeatedly; only the first call does work.
    """
    global _warmup_widget
    if _warmup_widget is not None:
        return True
    try:
        from PySide6.QtWebEngineWidgets import QWebEngineView

        holder = QWidget(parent)
        holder.setFixedSize(0, 0)
        holder.hide()
        view = QWebEngineView(holder)
        _configure_web_settings(view.settings())
        view.setHtml("<!DOCTYPE html><html><body></body></html>")
        _warmup_widget = holder
        logger.debug("WebEngine warmup complete")
        return True
    except Exception:
        logger.debug("WebEngine warmup skipped", exc_info=True)
        return False


def _configure_web_settings(settings) -> None:
    import sys

    from PySide6.QtWebEngineCore import QWebEngineProfile, QWebEngineSettings

    QWebEngineProfile.defaultProfile().setSpellCheckEnabled(False)

    attrs = [
        QWebEngineSettings.WebAttribute.LocalContentCanAccessRemoteUrls,
        QWebEngineSettings.WebAttribute.LocalContentCanAccessFileUrls,
        QWebEngineSettings.WebAttribute.JavascriptEnabled,
        QWebEngineSettings.WebAttribute.LocalStorageEnabled,
    ]
    if sys.platform != "darwin":
        attrs.append(QWebEngineSettings.WebAttribute.Accelerated2dCanvasEnabled)
        if hasattr(QWebEngineSettings.WebAttribute, "WebGLEnabled"):
            attrs.append(QWebEngineSettings.WebAttribute.WebGLEnabled)
    for attribute in attrs:
        settings.setAttribute(attribute, True)


class LazyHtmlPreview(QWidget):
    """Defers Qt WebEngine creation until the first HTML file is opened."""

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
        if not use_embedded_html_preview():
            if self._inner is None:
                self._inner = HtmlFallbackPreview(
                    self,
                    reason="Embedded HTML preview disabled (TEVIEWER_EMBED_HTML=0).",
                )
                self._status.hide()
                self._layout.addWidget(self._inner, stretch=1)
            self._inner.load_html_file(path)
            return
        if self._inner is None:
            self._inner = self._create_inner_viewer()
            self._status.hide()
            self._layout.addWidget(self._inner, stretch=1)

        if hasattr(self._inner, "load_html_file"):
            self._inner.load_html_file(path, force=force)
        if hasattr(self._inner, "is_view_active") and not self._inner.is_view_active():
            return
        if hasattr(self._inner, "set_zoom_percent"):
            self._inner.set_zoom_percent(self._zoom_percent)
        if hasattr(self._inner, "set_fit_mode"):
            self._inner.set_fit_mode(self._fit_mode)

    def load_url(self, url: str) -> None:
        """Load an http(s) page (for example igv.js served on localhost)."""
        if not use_embedded_html_preview():
            if self._inner is None:
                self._inner = HtmlFallbackPreview(
                    self,
                    reason="Embedded HTML preview disabled (TEVIEWER_EMBED_HTML=0).",
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

    def teardown(self) -> None:
        """Release WebEngine before showing a non-HTML preview in the same tab."""
        if self._inner is None:
            return
        if hasattr(self._inner, "teardown"):
            self._inner.teardown()
        self._layout.removeWidget(self._inner)
        self._inner.deleteLater()
        self._inner = None
        self._status.show()

    def _create_inner_viewer(self) -> QWidget:
        try:
            from PySide6.QtWebEngineWidgets import QWebEngineView  # noqa: F401

            return HtmlWebEnginePreview(self)
        except ImportError:
            return HtmlFallbackPreview(self, reason="PySide6-Addons not installed")
        except Exception as exc:  # noqa: BLE001
            logger.exception("WebEngine init failed")
            return HtmlFallbackPreview(self, reason=str(exc))


class HtmlWebEnginePreview(QWidget):
    """In-app viewer for Plotly HTML (TEAid, genome maps)."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        from PySide6.QtWebEngineWidgets import QWebEngineView

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)

        toolbar = QHBoxLayout()
        self._add_igv_track_btn = QPushButton("Add track…")
        self._add_igv_track_btn.setToolTip(
            "Add GFF, GTF, BED, BAM, or SAM track to the igv.js viewer"
        )
        self._add_igv_track_btn.clicked.connect(self._add_igv_track)
        self._add_igv_track_btn.hide()
        toolbar.addWidget(self._add_igv_track_btn)
        toolbar.addStretch()
        open_browser = QPushButton("Open in browser")
        open_browser.clicked.connect(self._open_in_browser)
        toolbar.addWidget(open_browser)
        layout.addLayout(toolbar)

        self._view = QWebEngineView(self)
        self._view.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        layout.addWidget(self._view, stretch=1)

        _configure_web_settings(self._view.settings())
        self._view.loadFinished.connect(self._on_load_finished)
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
        self._work_generation = 0

    def is_view_active(self) -> bool:
        return _qobject_alive(self._view)

    def _view_active(self) -> bool:
        return _qobject_alive(self._view) and _qobject_alive(self._view.page())

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
        vw = max(1, vp.width())
        vh = max(1, vp.height())
        cw = max(1, self._content_w)
        ch = max(1, self._content_h)

        if self._fit_mode == PreviewFitMode.ORIGINAL:
            base = 1.0
        elif self._fit_mode == PreviewFitMode.FIT_WIDTH:
            base = vw / cw
        else:
            base = min(vw / cw, vh / ch)

        factor = base * (self._zoom_percent / 100.0)
        factor = max(0.05, min(5.0, factor))
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

        try:
            self._view.load(_local_file_url(resolved))
        except RuntimeError:
            return
        self._sync_igv_toolbar()

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
        try:
            html = resolved.read_text(encoding="utf-8", errors="replace")
        except OSError as exc:
            self._view.setHtml(f"<p>Could not read HTML file: {exc}</p>")
            return
        base_url = _local_file_url(resolved.parent)

        def _load_fallback() -> None:
            if not self._view_active():
                return
            try:
                self._view.setHtml(html, base_url)
            except RuntimeError:
                return
            self._last_load_ok = True
            self._apply_zoom()
            self._defer(50, self._query_content_size)
            self._defer(250, self._query_content_size)
            self._defer(600, self._query_content_size)

        self._defer(0, _load_fallback)

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

    def _open_in_browser(self) -> None:
        if self._current_url:
            QDesktopServices.openUrl(QUrl(self._current_url))
            return
        if self._current_path and self._current_path.is_file():
            QDesktopServices.openUrl(_local_file_url(self._current_path))

    def teardown(self) -> None:
        self._invalidate_pending_work()
        try:
            self._view.loadFinished.disconnect(self._on_load_finished)
        except (RuntimeError, TypeError):
            pass
        try:
            self._view.stop()
            self._view.setUrl(QUrl("about:blank"))
        except RuntimeError:
            logger.debug("WebEngine teardown stop failed", exc_info=True)
        except Exception:
            logger.debug("WebEngine teardown stop failed", exc_info=True)
        self._current_path = None
        self._current_mtime_ns = None
        self._last_load_ok = False


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
        self._open_btn = QPushButton("Open in browser")
        self._open_btn.clicked.connect(self._open_in_browser)
        layout.addWidget(self._label)
        layout.addWidget(self._open_btn)
        layout.addStretch()
        self._current_path: Path | None = None
        self._current_url: str | None = None

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
