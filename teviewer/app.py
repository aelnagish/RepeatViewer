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

        # Organisation / application identity must be set before any QSettings
        # read, so preferences resolve to the same store the dialog writes to.
        app.setApplicationName(APP_DISPLAY_NAME)
        app.setApplicationDisplayName(APP_DISPLAY_NAME)
        QGuiApplication.setDesktopFileName(APP_DISPLAY_NAME)
        app.setOrganizationName("TEtrimmer")
        app.setOrganizationDomain("github.com/qjiangzhao/TEtrimmer")

        # Consistent cross-platform look: Fusion + pyqtdarktheme (auto light/dark),
        # with the accents in teviewer.ui.theme layered on top per-window.
        from teviewer.app_preferences import apply_log_level, setup_file_logging
        from teviewer.ui.theme_mode import apply_theme

        apply_theme()
        setup_file_logging()
        apply_log_level()

        from teviewer.views.dialog_polish import install_dialog_polisher

        install_dialog_polisher(app)

        import logging

        from teviewer import __version__

        logging.getLogger(__name__).info("%s %s starting", APP_DISPLAY_NAME, __version__)
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
