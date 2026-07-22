"""QApplication bootstrap."""

from __future__ import annotations

import faulthandler
import logging
import sys

from teviewer.branding import APP_DISPLAY_NAME
from teviewer.qt_runtime import apply_macos_app_name, configure_qt_runtime

configure_qt_runtime()

from PySide6.QtCore import Qt
from PySide6.QtGui import QGuiApplication, QPalette, QPixmap
from PySide6.QtWidgets import QApplication, QSplashScreen

from teviewer.widgets.app_icon import app_icon

_FAULT_HANDLER_FILE = None


def run_app(argv: list[str] | None = None) -> int:
    if argv is None:
        argv = sys.argv
    write_exception_report_fn = None

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

        from teviewer.qt_runtime import install_macos_stderr_fd_filter

        install_macos_stderr_fd_filter()

        from teviewer.widgets.browser_dock import install_global_click_log

        install_global_click_log(app)

        # Organisation / application identity must be set before any QSettings
        # read, so preferences resolve to the same store the dialog writes to.
        app.setApplicationName(APP_DISPLAY_NAME)
        app.setApplicationDisplayName(APP_DISPLAY_NAME)
        QGuiApplication.setDesktopFileName(APP_DISPLAY_NAME)
        app.setOrganizationName("TEtrimmer")
        app.setOrganizationDomain("github.com/qjiangzhao/TEtrimmer")

        # Consistent cross-platform look: Fusion + light/dark QPalette
        # (teviewer.ui.theme accents layered on top per-window).
        from teviewer.app_preferences import (
            apply_log_level,
            crash_reports_dir,
            install_exception_logging_hooks,
            log_file_path,
            setup_file_logging,
            write_exception_report,
        )
        write_exception_report_fn = write_exception_report
        from teviewer.ui.theme_mode import apply_theme

        apply_theme()
        setup_file_logging()
        apply_log_level()
        install_exception_logging_hooks()

        global _FAULT_HANDLER_FILE
        try:
            crash_dir = crash_reports_dir()
            crash_dir.mkdir(parents=True, exist_ok=True)
            fault_log = crash_dir / "fatal-faults.log"
            _FAULT_HANDLER_FILE = fault_log.open("a", encoding="utf-8")
            faulthandler.enable(file=_FAULT_HANDLER_FILE, all_threads=True)
        except OSError:
            _FAULT_HANDLER_FILE = None

        from teviewer.views.dialog_polish import install_dialog_polisher

        install_dialog_polisher(app)

        from teviewer import __version__

        logging.getLogger(__name__).info("%s %s starting", APP_DISPLAY_NAME, __version__)
        logging.getLogger(__name__).info("Diagnostics log: %s", log_file_path())
        icon = app_icon()
        if not icon.isNull():
            app.setWindowIcon(icon)

        splash_pix = QPixmap(480, 120)
        splash_pix.fill(app.palette().color(QPalette.ColorRole.Window))
        splash = QSplashScreen(splash_pix)
        splash.showMessage(
            f"Starting {APP_DISPLAY_NAME}…",
            Qt.AlignmentFlag.AlignBottom | Qt.AlignmentFlag.AlignHCenter,
            app.palette().color(QPalette.ColorRole.WindowText),
        )
        splash.show()
        app.processEvents()

        from teviewer.main_window import MainWindow

        window = MainWindow()
        if not icon.isNull():
            window.setWindowIcon(icon)
        window.show()
        app.processEvents()
        splash.finish(window)
        return app.exec()
    except Exception:
        import traceback

        try:
            exc_type, exc_value, exc_traceback = sys.exc_info()
            if (
                write_exception_report_fn is not None
                and exc_type is not None
                and exc_value is not None
                and exc_traceback is not None
            ):
                report = write_exception_report_fn(
                    title="startup-exception",
                    exc_type=exc_type,
                    exc_value=exc_value,
                    exc_traceback=exc_traceback,
                    context="run_app startup",
                )
                if report is not None:
                    logging.getLogger(__name__).error(
                        "Startup exception report written to %s", report
                    )
            logging.getLogger(__name__).exception("Application startup failed")
        except Exception:
            pass
        traceback.print_exc()
        return 1
