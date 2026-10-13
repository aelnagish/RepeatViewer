"""Actions → Dispatch / Restore must bring Preview back under Files after a move."""

from __future__ import annotations

from PySide6.QtCore import QRect, Qt
from PySide6.QtWidgets import (
    QApplication,
    QDockWidget,
    QLabel,
    QMainWindow,
    QWidget,
)


def _bind_host_methods(host: QMainWindow) -> None:
    from repeatviewer.main_window import MainWindow

    host.is_preview_window_dispatched = MainWindow.is_preview_window_dispatched.__get__(
        host, type(host)
    )
    host.toggle_preview_window_dispatch = (
        MainWindow.toggle_preview_window_dispatch.__get__(host, type(host))
    )
    host._preview_dispatch_split_sizes = (
        MainWindow._preview_dispatch_split_sizes.__get__(host, type(host))
    )
    host._restore_preview_after_dispatch = (
        MainWindow._restore_preview_after_dispatch.__get__(host, type(host))
    )
    host._apply_preview_dispatch_split_sizes = (
        MainWindow._apply_preview_dispatch_split_sizes.__get__(host, type(host))
    )
    host._sync_preview_dispatch_action = (
        MainWindow._sync_preview_dispatch_action.__get__(host, type(host))
    )
    host._redock_browser_panel = MainWindow._redock_browser_panel.__get__(
        host, type(host)
    )
    host._ensure_preview_visible = MainWindow._ensure_preview_visible.__get__(
        host, type(host)
    )
    host._ensure_preview_dock_attached = (
        MainWindow._ensure_preview_dock_attached.__get__(host, type(host))
    )
    host._make_floating_dock_independent = (
        MainWindow._make_floating_dock_independent.__get__(host, type(host))
    )
    host._detach_floating_dock_transient = (
        MainWindow._detach_floating_dock_transient
    )
    host._on_browser_dock_floating_changed = (
        MainWindow._on_browser_dock_floating_changed.__get__(host, type(host))
    )


def _make_host() -> QMainWindow:
    class Host(QMainWindow):
        def __init__(self) -> None:
            super().__init__()
            self._preview_tabs = None
            self._preview_host = QWidget(self)
            self._preview_dispatch_restore_sizes = None

            self._dock_folders = QDockWidget("Folders", self)
            self._dock_folders.setObjectName("dock_folders")
            self._dock_folders.setWidget(QLabel("folders", self._dock_folders))

            self._dock_files = QDockWidget("Files", self)
            self._dock_files.setObjectName("dock_files")
            self._dock_files.setWidget(QLabel("files", self._dock_files))

            self._dock_preview = QDockWidget("Preview", self)
            self._dock_preview.setObjectName("dock_preview")
            self._dock_preview.setWidget(self._preview_host)

            self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self._dock_folders)
            self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self._dock_files)
            self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self._dock_preview)
            self.splitDockWidget(
                self._dock_folders, self._dock_files, Qt.Orientation.Horizontal
            )
            self.splitDockWidget(
                self._dock_files, self._dock_preview, Qt.Orientation.Vertical
            )
            self.resize(900, 700)

    host = Host()
    _bind_host_methods(host)
    return host


def test_actions_dispatch_skips_independent_float_and_restores_after_move(
    qapp: QApplication,
) -> None:
    host = _make_host()
    host.show()
    qapp.processEvents()

    preview = host._dock_preview
    files = host._dock_files
    assert not preview.isFloating()
    assert preview.isVisible()
    assert preview.height() > 50

    # Dispatch (Actions path): normal float, no frameless independent window.
    host.toggle_preview_window_dispatch()
    qapp.processEvents()
    assert preview.isFloating()
    assert bool(preview.property("_rv_actions_dispatch"))
    assert not bool(preview.windowFlags() & Qt.WindowType.FramelessWindowHint)

    # Deferred promote must not turn Actions Dispatch into Frameless.
    host._make_floating_dock_independent(preview)
    qapp.processEvents()
    assert preview.isFloating()
    assert not bool(preview.windowFlags() & Qt.WindowType.FramelessWindowHint)

    # Move the floating window (the case that previously broke Restore).
    geo = preview.geometry()
    preview.setGeometry(QRect(geo.x() + 120, geo.y() + 80, geo.width(), geo.height()))
    qapp.processEvents()

    # Restore must put Preview back under Files with visible height.
    host.toggle_preview_window_dispatch()
    qapp.processEvents()
    qapp.processEvents()

    assert not preview.isFloating()
    assert preview.isVisible()
    assert files.isVisible()
    assert preview not in host.tabifiedDockWidgets(files)
    assert preview.height() > 80
    assert not bool(preview.property("_rv_actions_dispatch"))
    assert not bool(preview.property("_rv_float_restoring"))

    host.close()


def test_restore_survives_stale_float_promotion_race(qapp: QApplication) -> None:
    """Clearing dispatch flags mid-restore must not re-promote Preview."""
    host = _make_host()
    host.show()
    qapp.processEvents()

    preview = host._dock_preview
    host.toggle_preview_window_dispatch()
    qapp.processEvents()
    geo = preview.geometry()
    preview.setGeometry(QRect(geo.x() + 80, geo.y() + 60, geo.width(), geo.height()))
    qapp.processEvents()

    # Simulate the old bug: flags cleared before unfloat, then promote runs.
    preview.setProperty("_rv_actions_dispatch", False)
    preview.setProperty("_rv_float_restoring", True)
    host._make_floating_dock_independent(preview)
    qapp.processEvents()
    assert not bool(preview.windowFlags() & Qt.WindowType.FramelessWindowHint)

    preview.setProperty("_rv_float_restoring", False)
    preview.setProperty("_rv_actions_dispatch", True)
    host._restore_preview_after_dispatch()
    qapp.processEvents()

    assert not preview.isFloating()
    assert preview.isVisible()
    assert preview.height() > 80
    assert not bool(preview.windowFlags() & Qt.WindowType.FramelessWindowHint)

    host.close()


def test_restore_after_legacy_frameless_independent_float(qapp: QApplication) -> None:
    """Restore must undo Window|Frameless promotion from older builds / title-bar float."""
    host = _make_host()
    host.show()
    qapp.processEvents()

    preview = host._dock_preview
    preview.setProperty("_rv_actions_dispatch", True)
    preview.setFloating(True)
    qapp.processEvents()
    preview.hide()
    preview.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.FramelessWindowHint)
    preview.setProperty("_rv_float_independent", True)
    preview.show()
    qapp.processEvents()
    geo = preview.geometry()
    preview.setGeometry(QRect(geo.x() + 90, geo.y() + 70, geo.width(), geo.height()))
    qapp.processEvents()

    host._restore_preview_after_dispatch()
    qapp.processEvents()
    qapp.processEvents()

    assert not preview.isFloating()
    assert preview.isVisible()
    assert preview.height() > 80
    assert not bool(preview.windowFlags() & Qt.WindowType.FramelessWindowHint)
    assert not bool(preview.property("_rv_float_independent"))

    host.close()
