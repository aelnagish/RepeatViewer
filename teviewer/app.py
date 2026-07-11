"""QApplication bootstrap."""

from __future__ import annotations

import sys

from teviewer.branding import APP_DISPLAY_NAME
from teviewer.qt_runtime import apply_macos_app_name, configure_qt_runtime

configure_qt_runtime()

from PySide6.QtCore import Qt
from PySide6.QtGui import QGuiApplication, QPixmap
from PySide6.QtWidgets import QApplication, QSplashScreen

from teviewer.widgets.app_icon import app_icon


def run_app(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv

    try:
        if sys.platform == "darwin":
            argv = list(argv)
            argv[0] = APP_DISPLAY_NAME
            apply_macos_app_name()

        # High-DPI on Qt6 is enabled by default; keep AA flag for older Qt builds.
        QApplication.setHighDpiScaleFactorRoundingPolicy(
            Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
        )

        app = QApplication(argv)
        app.setApplicationName(APP_DISPLAY_NAME)
        app.setApplicationDisplayName(APP_DISPLAY_NAME)
        QGuiApplication.setDesktopFileName(APP_DISPLAY_NAME)
        app.setOrganizationName("TEtrimmer")
        app.setOrganizationDomain("github.com/qjiangzhao/TEtrimmer")
        icon = app_icon()
        if not icon.isNull():
            app.setWindowIcon(icon)

        splash = QSplashScreen(QPixmap(480, 120))
        splash.showMessage(
            f"Starting {APP_DISPLAY_NAME}…",
            Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignHCenter,
            Qt.GlobalColor.white,
        )
        splash.show()
        app.processEvents()

        from teviewer.main_window import MainWindow

        window = MainWindow()
        if not icon.isNull():
            window.setWindowIcon(icon)
        window.show()
        splash.finish(window)
        return app.exec()
    except Exception:
        import traceback

        traceback.print_exc()
        return 1
