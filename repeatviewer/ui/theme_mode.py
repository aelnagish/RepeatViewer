"""Application theme (light/dark) built on Fusion + QPalette.

Default look (what forks should see):
  1. Qt "Fusion" style  -> consistent chrome on macOS / Windows / Linux
  2. Fusion QPalette    -> light/dark colours
  3. MainWindow QSS     -> accents in ``repeatviewer.ui.theme`` on top of that

Optional flat theme: set env ``REPEATVIEWER_USE_QDARKTHEME=1`` and install
``pyqtdarktheme-fork`` — not used by default so installs match across machines.

Corners are always rounded (no sharp option).
"""

from __future__ import annotations

import os
import sys

from PySide6.QtCore import QEvent, QObject, Qt, QSettings
from PySide6.QtGui import QColor, QGuiApplication, QPalette
from PySide6.QtWidgets import (
    QAbstractItemView,
    QApplication,
    QComboBox,
    QMenu,
    QProxyStyle,
    QStyle,
    QStyledItemDelegate,
    QStyleFactory,
    QWidget,
)

# Hover → submenu open delay (ms). Fusion's default feels sluggish for Tools.
_MENU_SUBMENU_POPUP_DELAY_MS = 50

THEME_SETTINGS_GROUP = "appearance"
DEFAULT_THEME_MODE = "auto"
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
    """Always rounded — sharp corners are no longer offered."""
    return DEFAULT_CORNER_SHAPE


def save_corner_shape(value: str) -> None:
    s = QSettings()
    s.beginGroup(THEME_SETTINGS_GROUP)
    s.setValue("corner_shape", DEFAULT_CORNER_SHAPE)
    s.endGroup()
    del value


_APP_STYLE = None


# Popups drawn by the stylesheet with rounded corners (theme.py). Their window must be
# transparent or the corners stay square behind the rounded border.
_ROUNDED_POPUP_CLASSES = ("QComboBoxPrivateContainer", "QTipLabel")
_ROUNDED_POPUP_NAMES = ("previewScaleSliderPopup",)
POPUP_RADIUS = 8

# Completer lists are stand-alone popup views the app stylesheet can't single out,
# so they get this sheet directly. The rounded panel behind them is painted by
# _POPUP_PAINTER; the list itself is transparent on top.
_POPUP_LIST_CSS = (
    "QAbstractItemView { background: transparent; border: none; padding: 5px;"
    " outline: 0; }"
    "QAbstractItemView::item { padding: 3px 8px; border-radius: 5px; min-height: 18px; }"
    "QAbstractItemView::item:selected { background: palette(highlight);"
    " color: palette(highlighted-text); }"
)


def paint_rounded_popup(widget) -> None:
    """Background for a hand-made popup frame: rounded, base colour, hairline border."""
    from PySide6.QtCore import QRectF
    from PySide6.QtGui import QPainter, QPen

    p = QPainter(widget)
    p.setRenderHint(QPainter.RenderHint.Antialiasing)
    # app palette: a popup list's own stylesheet ('transparent') overrides its palette
    pal = QApplication.palette()
    p.setPen(QPen(pal.color(pal.ColorRole.Midlight), 1))
    p.setBrush(pal.color(pal.ColorRole.Base))
    p.drawRoundedRect(QRectF(widget.rect()).adjusted(0.5, 0.5, -0.5, -0.5),
                      POPUP_RADIUS, POPUP_RADIUS)
    p.end()


class _PopupPainter(QObject):
    """Paints the rounded panel for popup windows whose own painting is square: the
    combo-box list container (Fusion draws a square frame there) and completer lists."""

    def eventFilter(self, obj, event) -> bool:  # noqa: N802 (Qt API)
        if event.type() == QEvent.Type.Paint and isinstance(obj, QWidget):
            paint_rounded_popup(obj)
            return True
        return False


_POPUP_PAINTER: _PopupPainter | None = None


def _popup_painter() -> _PopupPainter:
    """One shared painter, owned by the application so Qt deletes it with the app
    (a module-level QObject outliving QApplication crashed Python at exit)."""
    global _POPUP_PAINTER
    try:
        from shiboken6 import isValid

        alive = _POPUP_PAINTER is not None and isValid(_POPUP_PAINTER)
    except Exception:  # noqa: BLE001
        alive = _POPUP_PAINTER is not None
    if not alive:
        _POPUP_PAINTER = _PopupPainter(QApplication.instance())
    return _POPUP_PAINTER


_IS_WINDOWS = sys.platform == "win32"
_WINDOWS_POPUP_FLAGS = (
    Qt.WindowType.FramelessWindowHint | Qt.WindowType.NoDropShadowWindowHint
)


def _windows_transparent_popup(widget: QWidget) -> None:
    """Windows only: a translucent popup must also be frameless and without the
    system drop shadow. Otherwise Windows fills the area outside the rounded
    panel (and the class drop shadow) black: the dark corners on menus, combo
    lists and tooltips. Set before the native window exists."""
    if (widget.windowFlags() & _WINDOWS_POPUP_FLAGS) == _WINDOWS_POPUP_FLAGS:
        return
    if widget.property("_rv_win_popup_flags"):
        return  # re-entered from setWindowFlags (re-polish)
    widget.setProperty("_rv_win_popup_flags", True)
    try:
        widget.setWindowFlags(widget.windowFlags() | _WINDOWS_POPUP_FLAGS)
    finally:
        widget.setProperty("_rv_win_popup_flags", False)


class _MenuStyle(QProxyStyle):
    """Fusion, with every popup (menus, combo-box lists, tooltips, completer lists, the
    preview zoom slider) on a transparent window, so the stylesheet's rounded corners
    show. One look for all popups across the app."""

    def styleHint(self, hint, option=None, widget=None, returnData=None):  # noqa: N802
        if hint == QStyle.StyleHint.SH_Menu_SubMenuPopupDelay:
            return _MENU_SUBMENU_POPUP_DELAY_MS
        return super().styleHint(hint, option, widget, returnData)

    def polish(self, target):  # Qt override: QWidget, QPalette or QApplication
        if isinstance(target, QWidget):
            rounded = (
                isinstance(target, QMenu)
                or type(target).__name__ in _ROUNDED_POPUP_CLASSES
                or target.metaObject().className() in _ROUNDED_POPUP_CLASSES
                or target.objectName() in _ROUNDED_POPUP_NAMES
            )
            if (not rounded and isinstance(target, QAbstractItemView)
                    and target.isWindow()
                    and target.windowType() == Qt.WindowType.Popup):
                # QCompleter list: a stand-alone popup view
                rounded = True
                if not target.styleSheet():
                    target.setStyleSheet(_POPUP_LIST_CSS)
                target.installEventFilter(_popup_painter())
            if isinstance(target, QComboBox):
                # The internal combo delegates ignore the stylesheet's rounded item
                # highlight; the standard one follows it. Custom delegates are kept.
                name = target.itemDelegate().metaObject().className()
                if name in ("QComboMenuDelegate", "QComboBoxDelegate"):
                    target.setItemDelegate(QStyledItemDelegate(target))
            if target.metaObject().className() == "QComboBoxPrivateContainer":
                target.installEventFilter(_popup_painter())
            if rounded and _IS_WINDOWS and target.isWindow():
                _windows_transparent_popup(target)
            if rounded and not target.testAttribute(
                    Qt.WidgetAttribute.WA_TranslucentBackground):
                target.setAttribute(Qt.WidgetAttribute.WA_TranslucentBackground, True)
        return super().polish(target)


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

    fusion = QStyleFactory.create("Fusion")
    if fusion is not None:
        global _APP_STYLE
        _APP_STYLE = _MenuStyle(fusion)
        app.setStyle(_APP_STYLE)

    _apply_color_scheme(mode)
    concrete = _resolve_os_theme() if mode == "auto" else mode

    # Default: Fusion palette (stable look for all platforms / forks).
    # Opt-in flat theme only when explicitly requested.
    use_qdark = os.environ.get("REPEATVIEWER_USE_QDARKTHEME", "").strip().lower() in {
        "1",
        "true",
        "yes",
        "on",
    }
    if use_qdark:
        try:
            import qdarktheme
        except ImportError:
            use_qdark = False
        else:
            _patch_qdarktheme_icon_engine()
            qdarktheme.setup_theme(
                concrete,
                corner_shape="rounded",
                custom_colors={"primary": load_accent_color()},
            )
            _repolish(app)
            return

    _apply_fusion_palette(app, concrete)
    _repolish(app)


def _apply_fusion_palette(app: QApplication, concrete: str) -> None:
    """Dark/light Fusion palette when pyqtdarktheme is not installed."""
    palette = QPalette()
    accent = QColor(load_accent_color())
    if concrete == "dark":
        window = QColor(45, 45, 45)
        base = QColor(35, 35, 35)
        alt = QColor(53, 53, 53)
        text = QColor(220, 220, 220)
        disabled = QColor(127, 127, 127)
        mid = QColor(66, 66, 66)
        highlight = accent
        highlighted_text = QColor(255, 255, 255)
        button = QColor(53, 53, 53)
        link = accent.lighter(120)
        placeholder = QColor(255, 255, 255)
    else:
        window = QColor(240, 240, 240)
        base = QColor(255, 255, 255)
        alt = QColor(247, 247, 247)
        text = QColor(33, 33, 33)
        disabled = QColor(160, 160, 160)
        mid = QColor(200, 200, 200)
        highlight = accent
        highlighted_text = QColor(255, 255, 255)
        button = QColor(240, 240, 240)
        link = accent.darker(110)
        placeholder = QColor(140, 140, 140)

    palette.setColor(QPalette.ColorRole.Window, window)
    palette.setColor(QPalette.ColorRole.WindowText, text)
    palette.setColor(QPalette.ColorRole.Base, base)
    palette.setColor(QPalette.ColorRole.AlternateBase, alt)
    palette.setColor(QPalette.ColorRole.ToolTipBase, base)
    palette.setColor(QPalette.ColorRole.ToolTipText, text)
    palette.setColor(QPalette.ColorRole.Text, text)
    palette.setColor(QPalette.ColorRole.Button, button)
    palette.setColor(QPalette.ColorRole.ButtonText, text)
    palette.setColor(QPalette.ColorRole.BrightText, QColor(255, 255, 255))
    palette.setColor(QPalette.ColorRole.Link, link)
    palette.setColor(QPalette.ColorRole.Highlight, highlight)
    palette.setColor(QPalette.ColorRole.HighlightedText, highlighted_text)
    palette.setColor(QPalette.ColorRole.PlaceholderText, placeholder)
    palette.setColor(QPalette.ColorRole.Mid, mid)
    palette.setColor(QPalette.ColorRole.Midlight, alt)
    palette.setColor(QPalette.ColorRole.Dark, mid)
    palette.setColor(QPalette.ColorRole.Light, window.lighter(110))
    palette.setColor(QPalette.ColorRole.Shadow, QColor(20, 20, 20))
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.Text, disabled)
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.WindowText, disabled)
    palette.setColor(QPalette.ColorGroup.Disabled, QPalette.ColorRole.ButtonText, disabled)
    app.setPalette(palette)


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
    if getattr(SvgIconEngine, "_repeatviewer_patched", False):
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
    SvgIconEngine._repeatviewer_patched = True


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
