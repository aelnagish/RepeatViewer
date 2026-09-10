"""Keep top-level windows and floating docks on a visible screen."""

from __future__ import annotations

from PySide6.QtCore import QRect, QSize
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


def _decoration_size(widget: QWidget) -> QSize:
    """Estimate non-client (title bar / borders) size for a top-level window."""
    frame = widget.frameGeometry()
    geo = widget.geometry()
    if not frame.isValid() or not geo.isValid():
        return QSize(0, 0)
    return QSize(
        max(0, frame.width() - geo.width()),
        max(0, frame.height() - geo.height()),
    )


def fit_widget_to_available_screen(widget: QWidget | None) -> bool:
    """Resize/move ``widget`` so its *frame* fits in the screen's available area.

    On Windows, a restored or default geometry often exceeds ``availableGeometry``
    (taskbar / DPI), which hides the OS title bar (close) and status bar (Stop).
    ``setGeometry`` adjusts the client area; we account for frame decorations.
    """
    if widget is None:
        return False

    screen = widget.screen()
    if screen is None:
        screen = QGuiApplication.primaryScreen()
    if screen is None:
        return False

    available = screen.availableGeometry()
    if not available.isValid():
        return False

    decor = _decoration_size(widget)
    max_client_w = max(320, available.width() - decor.width())
    max_client_h = max(240, available.height() - decor.height())

    # Cap QWidget minimumSize so Windows can actually shrink to fit.
    min_size = widget.minimumSize()
    capped_min_w = min(min_size.width(), max_client_w) if min_size.width() > 0 else 0
    capped_min_h = min(min_size.height(), max_client_h) if min_size.height() > 0 else 0
    if capped_min_w != min_size.width() or capped_min_h != min_size.height():
        widget.setMinimumSize(capped_min_w, capped_min_h)

    client = widget.geometry()
    new_w = min(max(client.width(), capped_min_w or 320), max_client_w)
    new_h = min(max(client.height(), capped_min_h or 240), max_client_h)

    # Prefer keeping the current top-left if it already sits on this screen;
    # otherwise center within available.
    frame = widget.frameGeometry()
    x = frame.x()
    y = frame.y()
    if not available.contains(frame.topLeft()):
        x = available.x() + max(0, (available.width() - (new_w + decor.width())) // 2)
        y = available.y() + max(0, (available.height() - (new_h + decor.height())) // 2)
    if x < available.x():
        x = available.x()
    if y < available.y():
        y = available.y()
    if x + new_w + decor.width() > available.right() + 1:
        x = max(available.x(), available.right() + 1 - (new_w + decor.width()))
    if y + new_h + decor.height() > available.bottom() + 1:
        y = max(available.y(), available.bottom() + 1 - (new_h + decor.height()))

    changed = (
        new_w != client.width()
        or new_h != client.height()
        or x != frame.x()
        or y != frame.y()
    )
    if not changed:
        return False

    widget.resize(new_w, new_h)
    frame_after = widget.frameGeometry()
    dx = x - frame_after.x()
    dy = y - frame_after.y()
    if dx or dy:
        widget.move(widget.x() + dx, widget.y() + dy)
    return True


def ensure_widget_on_screen(widget: QWidget | None, *, min_visible: int = 120) -> bool:
    """Move ``widget`` onto a visible screen if needed. Returns True when adjusted."""
    if widget is None:
        return False
    # Prefer a full fit (keeps title/status bars visible) over a partial nudge.
    if fit_widget_to_available_screen(widget):
        return True
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
