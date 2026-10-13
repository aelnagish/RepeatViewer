"""Hiding a floating browser dock must leave it recoverable via View / Reset layout."""

from __future__ import annotations

from PySide6.QtCore import Qt
from PySide6.QtWidgets import QApplication, QDockWidget, QLabel, QMainWindow


def test_hide_floating_dock_redocks_before_hide(qapp: QApplication) -> None:
    from repeatviewer.main_window import MainWindow

    # Use a minimal stand-in: call the same helpers MainWindow uses.
    class Host(QMainWindow):
        def __init__(self) -> None:
            super().__init__()
            self._dock = QDockWidget("Preview", self)
            self._dock.setObjectName("dock_preview")
            self._dock.setWidget(QLabel("content", self._dock))
            self.addDockWidget(Qt.DockWidgetArea.LeftDockWidgetArea, self._dock)
            self._dock_folders = None
            self._dock_files = None
            self._dock_preview = self._dock

        hide_browser_dock = MainWindow.hide_browser_dock
        _redock_browser_panel = MainWindow._redock_browser_panel

    host = Host()
    host.resize(800, 600)
    host.show()
    qapp.processEvents()

    dock = host._dock
    dock.setFloating(True)
    dock.setWindowFlags(Qt.WindowType.Window | Qt.WindowType.FramelessWindowHint)
    dock.show()
    qapp.processEvents()
    assert dock.isFloating()

    host.hide_browser_dock(dock)
    qapp.processEvents()

    assert not dock.isVisible()
    assert not dock.isFloating()
    assert bool(dock.property("_rv_user_hidden"))
    # View toggle / Reset layout can show it again inside the main window.
    dock.setProperty("_rv_user_hidden", False)
    dock.show()
    qapp.processEvents()
    assert dock.isVisible()
    assert not dock.isFloating()

    host.close()
