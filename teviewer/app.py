"""QApplication bootstrap."""

from __future__ import annotations

import sys

from teviewer.qt_runtime import apply_macos_app_name, configure_qt_runtime

configure_qt_runtime()

from PySide6.QtCore import Qt
from PySide6.QtGui import QGuiApplication, QPixmap
from PySide6.QtWidgets import QApplication, QSplashScreen


def run_app(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv

    try:
        if sys.platform == "darwin":
            argv = list(argv)
            argv[0] = "TEviewer"

        # High-DPI on Qt6 is enabled by default; keep AA flag for older Qt builds.
        QApplication.setHighDpiScaleFactorRoundingPolicy(
            Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
        )

        app = QApplication(argv)
        app.setApplicationName("TEviewer")
        app.setApplicationDisplayName("TEviewer")
        QGuiApplication.setDesktopFileName("TEviewer")
        app.setOrganizationName("TEtrimmer")
        app.setOrganizationDomain("github.com/qjiangzhao/TEtrimmer")
        apply_macos_app_name()

        splash = QSplashScreen(QPixmap(480, 120))
        splash.showMessage(
            "Starting TEviewer…",
            Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignHCenter,
            Qt.GlobalColor.white,
        )
        splash.show()
        app.processEvents()

        from teviewer.main_window import MainWindow

        window = MainWindow()
        window.show()
        splash.finish(window)
        return app.exec()
    except Exception:
        import traceback

        traceback.print_exc()
        return 1
