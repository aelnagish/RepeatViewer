"""Application theme (light/dark) built on Fusion + pyqtdarktheme.

The look is layered:
  1. Qt "Fusion" style  -> consistent chrome on macOS / Windows / Linux
  2. pyqtdarktheme       -> modern flat palette + base stylesheet, auto light/dark
  3. MainWindow QSS      -> the accents in ``teviewer.ui.theme`` on top of that

``pyqtdarktheme`` is an optional dependency: if it is missing we still fall
back to plain Fusion so the app always launches.
"""

from __future__ import annotations

from PySide6.QtCore import Qt, QSettings
from PySide6.QtGui import QGuiApplication
from PySide6.QtWidgets import QApplication

THEME_SETTINGS_GROUP = "appearance"
DEFAULT_THEME_MODE = "dark"
THEME_MODES = ("auto", "light", "dark")
THEME_LABELS = (("Auto (follow system)", "auto"), ("Light", "light"), ("Dark", "dark"))

DEFAULT_ACCENT_COLOR = "#3d6ea8"
DEFAULT_CORNER_SHAPE = "rounded"
CORNER_SHAPES = ("rounded", "sharp")
CORNER_LABELS = (("Rounded", "rounded"), ("Sharp", "sharp"))


def load_theme_mode() -> str:
    s = QSettings()
    s.beginGroup(THEME_SETTINGS_GROUP)
    value = str(s.value("theme_mode", DEFAULT_THEME_MODE))
    s.endGroup()
    return value if value in THEME_MODES else DEFAULT_THEME_MODE


def save_theme_mode(mode: str) -> None:
    if mode not in THEME_MODES:
        mode = DEFAULT_THEME_MODE
    s = QSettings()
    s.beginGroup(THEME_SETTINGS_GROUP)
    s.setValue("theme_mode", mode)
    s.endGroup()


def _valid_hex(value: str) -> bool:
    if not isinstance(value, str) or not value.startswith("#") or len(value) != 7:
        return False
    try:
        int(value[1:], 16)
    except ValueError:
        return False
    return True


def load_accent_color() -> str:
    s = QSettings()
    s.beginGroup(THEME_SETTINGS_GROUP)
    value = str(s.value("accent_color", DEFAULT_ACCENT_COLOR))
    s.endGroup()
    return value if _valid_hex(value) else DEFAULT_ACCENT_COLOR


def save_accent_color(value: str) -> None:
    if not _valid_hex(value):
        value = DEFAULT_ACCENT_COLOR
    s = QSettings()
    s.beginGroup(THEME_SETTINGS_GROUP)
    s.setValue("accent_color", value)
    s.endGroup()


def load_corner_shape() -> str:
    s = QSettings()
    s.beginGroup(THEME_SETTINGS_GROUP)
    value = str(s.value("corner_shape", DEFAULT_CORNER_SHAPE))
    s.endGroup()
    return value if value in CORNER_SHAPES else DEFAULT_CORNER_SHAPE


def save_corner_shape(value: str) -> None:
    if value not in CORNER_SHAPES:
        value = DEFAULT_CORNER_SHAPE
    s = QSettings()
    s.beginGroup(THEME_SETTINGS_GROUP)
    s.setValue("corner_shape", value)
    s.endGroup()


def apply_theme(mode: str | None = None) -> None:
    """Apply the theme to the current QApplication.

    Safe to call at startup and again whenever the setting changes.
    """
    if mode is None:
        mode = load_theme_mode()
    if mode not in THEME_MODES:
        mode = DEFAULT_THEME_MODE

    app = QApplication.instance()
    if app is None:
        return

    # 1) Set the platform colour-scheme FIRST (macOS title bar / native frames).
    #    Doing this before qdarktheme paints means qdarktheme's palette is applied
    #    last and always wins for the widget area, while the native chrome still
    #    follows light/dark correctly.
    _apply_color_scheme(mode)

    try:
        import qdarktheme
    except ImportError:
        # pyqtdarktheme not installed: still give consistent cross-platform chrome.
        app.setStyle("Fusion")
        return

    _patch_qdarktheme_icon_engine()

    # 2) setup_theme() sets the palette + base stylesheet (Fusion-based proxy
    #    style). We resolve "auto" to a concrete dark/light ourselves so
    #    qdarktheme never installs its macOS palette-change listener, which
    #    re-enters the app event filters and crashes on some Qt/macOS builds.
    concrete = _resolve_os_theme() if mode == "auto" else mode
    qdarktheme.setup_theme(
        concrete,
        corner_shape=load_corner_shape(),
        custom_colors={"primary": load_accent_color()},
    )

    # 3) Force already-shown widgets to re-read the new palette immediately so a
    #    live switch repaints without needing a restart.
    _repolish(app)


def _patch_qdarktheme_icon_engine() -> None:
    """Work around a pyqtdarktheme bug: its SvgIconEngine.pixmap() builds a
    QPainter on a zero-size pixmap (during layout/resize) and never ends it,
    spamming 'QPainter not active' errors and rendering some icons — notably the
    QToolBar overflow chevron — empty. Guard zero sizes and close the painter.
    """
    try:
        from qdarktheme._icon.icon_engine import SvgIconEngine
    except Exception:
        return
    if getattr(SvgIconEngine, "_teviewer_patched", False):
        return

    from PySide6.QtCore import QPoint, QRect, Qt
    from PySide6.QtGui import QImage, QPainter, QPixmap

    def pixmap(self, size, mode, state):
        min_size = min(size.width(), size.height())
        if min_size <= 0:
            return QPixmap()
        size.setWidth(min_size)
        size.setHeight(min_size)
        image = QImage(size, QImage.Format.Format_ARGB32)
        image.fill(Qt.GlobalColor.transparent)
        result = QPixmap.fromImage(image)
        painter = QPainter(result)
        try:
            self.paint(painter, QRect(QPoint(0, 0), size), mode, state)
        finally:
            painter.end()
        return result

    SvgIconEngine.pixmap = pixmap
    SvgIconEngine._teviewer_patched = True


def _resolve_os_theme() -> str:
    """Return 'dark' or 'light' for the current OS appearance (for 'auto')."""
    hints = QGuiApplication.styleHints()
    scheme = hints.colorScheme() if hasattr(hints, "colorScheme") else None
    if scheme == Qt.ColorScheme.Dark:
        return "dark"
    if scheme == Qt.ColorScheme.Light:
        return "light"
    try:
        import darkdetect

        return "dark" if darkdetect.isDark() else "light"
    except Exception:
        return "dark"


def _apply_color_scheme(mode: str) -> None:
    """Drive Qt's colour-scheme hint (Qt 6.8+) so native chrome follows the theme."""
    hints = QGuiApplication.styleHints()
    if not hasattr(hints, "setColorScheme"):
        return
    if mode == "light":
        hints.setColorScheme(Qt.ColorScheme.Light)
    elif mode == "dark":
        hints.setColorScheme(Qt.ColorScheme.Dark)
    else:  # auto -> follow the operating system
        hints.unsetColorScheme()


def _repolish(app: QApplication) -> None:
    from PySide6.QtWidgets import QWidget

    for widget in app.allWidgets():
        try:
            style = widget.style()
            style.unpolish(widget)
            style.polish(widget)
            # Call the base QWidget.update() explicitly: some custom widgets
            # shadow update() with a different signature (e.g. progress panels).
            QWidget.update(widget)
        except Exception:
            continue
