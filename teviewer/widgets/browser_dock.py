"""Helpers for browser dock widgets."""

from __future__ import annotations

import os

from PySide6.QtCore import QEvent, QObject, QPoint, QSize, Qt
from PySide6.QtGui import QMouseEvent
from PySide6.QtWidgets import QApplication, QDockWidget, QSizePolicy, QWidget


class ZeroDockTitleBar(QWidget):
    """Placeholder title bar so QDockWidget does not reserve a native title strip."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self.setFixedHeight(0)
        self.setSizePolicy(QSizePolicy.Policy.Ignored, QSizePolicy.Policy.Fixed)

    def sizeHint(self) -> QSize:  # noqa: N802
        return QSize(0, 0)

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return self.sizeHint()


class PreviewDockWidget(QDockWidget):
    """Preview dock that keeps panel content aligned with mouse hit-testing on macOS.

    QDockWidget can reserve ~24px for a title strip while painting the content
    widget from y=0, so clicks land below visible controls. Force the panel to
    fill the dock client rect and never install a custom title-bar widget.
    """

    def __init__(self, title: str, parent: QWidget | None = None) -> None:
        super().__init__(title, parent)
        self._layout_content_full_bleed()

    def setWidget(self, widget: QWidget | None) -> None:  # noqa: N802
        super().setWidget(widget)
        self._layout_content_full_bleed()

    def resizeEvent(self, event) -> None:  # noqa: N802
        super().resizeEvent(event)
        self._layout_content_full_bleed()

    def showEvent(self, event) -> None:  # noqa: N802
        super().showEvent(event)
        self._layout_content_full_bleed()

    def _layout_content_full_bleed(self) -> None:
        content = self.widget()
        if content is None:
            return
        content.setGeometry(self.rect())
        content.raise_()


class PreviewClickDebugFilter(QObject):
    """Log widget-under-cursor on click when ``TEVIEWER_DEBUG_CLICKS=1``."""

    def __init__(self, root: QWidget, status_callback) -> None:
        super().__init__(root)
        self._root = root
        self._status = status_callback
        root.installEventFilter(self)

    def eventFilter(self, watched: QObject, event: QEvent) -> bool:  # noqa: N802
        if watched is not self._root or event.type() != QEvent.Type.MouseButtonPress:
            return False
        if not isinstance(event, QMouseEvent):
            return False
        global_pos = event.globalPosition().toPoint()
        local = self._root.mapFromGlobal(global_pos)
        hit = QApplication.widgetAt(global_pos)
        names: list[str] = []
        widget = hit
        while widget is not None:
            name = widget.objectName() or widget.__class__.__name__
            geo = widget.geometry()
            top_left = widget.mapToGlobal(QPoint(0, 0))
            names.append(
                f"{name}@{top_left.y()}-{top_left.y() + geo.height()} "
                f"(local y={widget.mapFromGlobal(global_pos).y()})"
            )
            widget = widget.parentWidget()
            if widget is self._root.window():
                break
        self._status(
            f"click panel-local y={local.y()} · hit: {' → '.join(names[:4])}"
        )
        return False


def install_preview_click_debug(root: QWidget, status_callback) -> None:
    if os.environ.get("TEVIEWER_DEBUG_CLICKS", "").strip().lower() not in {
        "1",
        "true",
        "yes",
        "on",
    }:
        return
    PreviewClickDebugFilter(root, status_callback)


class _GlobalClickLogFilter(QObject):
    """Logs exactly what Qt delivers for every real mouse press, app-wide.

    Installed on QApplication so ``watched`` is the actual widget Qt has
    decided to deliver the event to (real OS dispatch, not a recomputation),
    which is what a synthetic/injected click cannot reproduce.
    """

    def __init__(self, app: QApplication, log_path: str) -> None:
        super().__init__(app)
        self._log_path = log_path
        with open(self._log_path, "a", encoding="utf-8") as fh:
            fh.write("--- TEviewer click log session start ---\n")

    def eventFilter(self, watched: QObject, event: QEvent) -> bool:  # noqa: N802
        if event.type() not in (QEvent.Type.MouseButtonPress, QEvent.Type.MouseButtonDblClick):
            return False
        if not isinstance(event, QMouseEvent):
            return False
        if not isinstance(watched, QWidget):
            return False

        global_pos = event.globalPosition().toPoint()
        local_in_watched = watched.mapFromGlobal(global_pos)
        watched_rect = watched.rect()
        inside = watched_rect.contains(local_in_watched)

        painted_hit = QApplication.widgetAt(global_pos)
        mismatch = painted_hit is not watched

        lines = [
            f"global={global_pos.x()},{global_pos.y()}",
            f"delivered_to={watched.__class__.__name__}#{watched.objectName()!r}",
            f"delivered_local={local_in_watched.x()},{local_in_watched.y()}",
            f"delivered_rect={watched_rect.width()}x{watched_rect.height()}",
            f"local_inside_rect={inside}",
            (
                f"widgetAt(global)={painted_hit.__class__.__name__}#{painted_hit.objectName()!r}"
                if painted_hit is not None
                else "widgetAt(global)=None"
            ),
            f"MISMATCH_delivered_vs_widgetAt={mismatch}",
        ]
        try:
            with open(self._log_path, "a", encoding="utf-8") as fh:
                fh.write(" | ".join(lines) + "\n")
        except OSError:
            pass
        return False


def install_global_click_log(app: QApplication) -> None:
    """Opt-in, app-wide click diagnostic: set TEVIEWER_CLICK_LOG=<path>.

    Unlike install_preview_click_debug, this sees every real click Qt
    dispatches anywhere in the app (not just ones landing on one root
    widget's own background), and it runs on genuine OS-delivered events.
    """
    log_path = os.environ.get("TEVIEWER_CLICK_LOG", "").strip()
    if not log_path:
        return
    filt = _GlobalClickLogFilter(app, log_path)
    app.installEventFilter(filt)
