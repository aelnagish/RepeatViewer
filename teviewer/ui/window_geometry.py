"""Keep top-level windows and floating docks on a visible screen."""

from __future__ import annotations

from PySide6.QtCore import QRect
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QWidget


def _screen_rects() -> list[QRect]:
    rects: list[QRect] = []
    for screen in QGuiApplication.screens():
        geo = screen.availableGeometry()
        if geo.isValid():
            rects.append(geo)
    return rects


def clamp_frame_to_screens(frame: QRect, *, min_visible: int = 120) -> QRect:
    """Return ``frame`` nudged onto a visible screen when it would sit off-screen."""
    if not frame.isValid():
        return frame

    screens = _screen_rects()
    if not screens:
        return frame

    best_rect = frame
    best_area = 0
    for available in screens:
        inter = frame.intersected(available)
        area = max(0, inter.width()) * max(0, inter.height())
        if area > best_area:
            best_area = area
            best_rect = available

    min_area = min_visible * min_visible
    if best_area >= min_area:
        available = best_rect
        width = min(frame.width(), available.width())
        height = min(frame.height(), available.height())
        x = frame.x()
        y = frame.y()
        if x < available.x():
            x = available.x()
        if y < available.y():
            y = available.y()
        if x + width > available.right() + 1:
            x = max(available.x(), available.right() + 1 - width)
        if y + height > available.bottom() + 1:
            y = max(available.y(), available.bottom() + 1 - height)
        return QRect(x, y, width, height)

    primary = QGuiApplication.primaryScreen()
    if primary is None:
        return frame
    available = primary.availableGeometry()
    width = min(max(frame.width(), 320), available.width())
    height = min(max(frame.height(), 240), available.height())
    x = available.x() + max(0, (available.width() - width) // 2)
    y = available.y() + max(0, (available.height() - height) // 2)
    return QRect(x, y, width, height)


def ensure_widget_on_screen(widget: QWidget | None, *, min_visible: int = 120) -> bool:
    """Move ``widget`` onto a visible screen if needed. Returns True when adjusted."""
    if widget is None:
        return False
    frame = widget.frameGeometry()
    clamped = clamp_frame_to_screens(frame, min_visible=min_visible)
    if clamped == frame:
        return False
    widget.setGeometry(clamped)
    return True


def ensure_floating_docks_on_screen(root: QWidget) -> None:
    """Clamp the main window and any floating dock panels to visible screens."""
    from PySide6.QtWidgets import QDockWidget

    ensure_widget_on_screen(root)
    for dock in root.findChildren(QDockWidget):
        if dock.isFloating() and dock.isVisible():
            ensure_widget_on_screen(dock)
