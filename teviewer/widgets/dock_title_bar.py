"""Custom QDockWidget title bars that host panel controls (filter, zoom, …)."""

from __future__ import annotations

from PySide6.QtCore import Qt, QSize, QTimer, QRect
from PySide6.QtGui import QMouseEvent
from teviewer.ui.window_geometry import clamp_frame_to_screens, ensure_widget_on_screen
from PySide6.QtWidgets import (
    QAbstractButton,
    QComboBox,
    QDockWidget,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QSizePolicy,
    QSlider,
    QSpinBox,
    QTabBar,
    QToolButton,
    QWidget,
)

# Compact title row; sizeHint must match so dock content is not laid out underneath.
BROWSER_DOCK_TITLE_HEIGHT = 24


class BrowserDockTitleBar(QWidget):
    """Title row: label · leading · stretch · trailing · close.

    Empty regions ignore mouse presses so QDockWidget can still drag the panel.
    Interactive children (line edits, sliders, buttons, tab bars) keep their own events.
    """

    def __init__(
        self,
        dock: QDockWidget,
        title: str,
        parent: QWidget | None = None,
        *,
        leading: QWidget | None = None,
        trailing: QWidget | None = None,
        leading_stretch: int = 0,
    ) -> None:
        super().__init__(parent)
        self._dock = dock
        self.setObjectName("BrowserDockTitleBar")
        self.setAttribute(Qt.WidgetAttribute.WA_StyledBackground, True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Fixed)
        self.setFixedHeight(BROWSER_DOCK_TITLE_HEIGHT)
        self.setMinimumHeight(BROWSER_DOCK_TITLE_HEIGHT)
        self.setMaximumHeight(BROWSER_DOCK_TITLE_HEIGHT)

        layout = QHBoxLayout(self)
        # 3px above/below so the filter's own bottom border is not clipped.
        layout.setContentsMargins(6, 3, 2, 3)
        layout.setSpacing(4)

        self._title = QLabel(title, self)
        self._title.setObjectName("BrowserDockTitleLabel")
        self._title.setToolTip(
            "Double-click to detach this panel to a separate window "
            "(uses a second screen when available)"
        )
        layout.addWidget(self._title, 0, Qt.AlignmentFlag.AlignVCenter)

        self._leading = leading
        if leading is not None:
            leading.setParent(self)
            layout.addWidget(
                leading,
                max(0, int(leading_stretch)),
                Qt.AlignmentFlag.AlignVCenter,
            )

        # When leading already expands (e.g. Preview file tabs), do not add another
        # stretch or it leaves a dead gap before Fit window / chrome.
        if int(leading_stretch) <= 0:
            layout.addStretch(1)

        if trailing is not None:
            trailing.setParent(self)
            # Prevent the dock from treating control clicks as title-bar drags.
            trailing.setAttribute(Qt.WidgetAttribute.WA_MacShowFocusRect, False)
            layout.addWidget(trailing, 0, Qt.AlignmentFlag.AlignVCenter)

        self._close_btn = self._make_chrome_button(
            tip="Hide panel",
            text="✕",
        )
        self._close_btn.clicked.connect(dock.close)
        layout.addWidget(self._close_btn, 0, Qt.AlignmentFlag.AlignVCenter)

        # Keep zoom / chrome above overflowing file tabs so the 6th tab cannot
        # paint on top of Reset / close.
        if trailing is not None:
            trailing.raise_()
        self._close_btn.raise_()

    def sizeHint(self) -> QSize:  # noqa: N802
        return QSize(max(super().sizeHint().width(), 80), BROWSER_DOCK_TITLE_HEIGHT)

    def minimumSizeHint(self) -> QSize:  # noqa: N802
        return QSize(0, BROWSER_DOCK_TITLE_HEIGHT)

    def _make_chrome_button(self, *, tip: str, text: str = "") -> QToolButton:
        btn = QToolButton(self)
        btn.setObjectName("BrowserDockChromeButton")
        btn.setAutoRaise(True)
        btn.setFocusPolicy(Qt.FocusPolicy.NoFocus)
        btn.setToolTip(tip)
        btn.setText(text)
        btn.setFixedSize(QSize(18, 18))
        return btn

    @staticmethod
    def _is_interactive(widget: QWidget | None) -> bool:
        current = widget
        while current is not None:
            if isinstance(
                current,
                (QSlider, QSpinBox, QLineEdit, QComboBox, QAbstractButton, QToolButton, QTabBar),
            ):
                return True
            name = current.objectName()
            if name in {
                "PreviewDockZoomHost",
                "PreviewDockFileTabs",
                "FilesDockFilterHost",
                "FilesDockNavHost",
                "FilesDockTrailing",
                "FoldersDockNavHost",
                "FoldersDockTrailing",
            }:
                return True
            current = current.parentWidget()
            if isinstance(current, BrowserDockTitleBar):
                break
        return False

    def mousePressEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
        if self._is_interactive(self.childAt(pos)):
            # Keep the event with the control; do not start a dock drag.
            event.accept()
            return
        # Empty title-bar space: let QDockWidget start a panel drag.
        event.ignore()

    def mouseDoubleClickEvent(self, event: QMouseEvent) -> None:  # noqa: N802
        pos = event.position().toPoint() if hasattr(event, "position") else event.pos()
        if self._is_interactive(self.childAt(pos)):
            event.accept()
            return
        self._toggle_dock_float()
        event.accept()

    def _toggle_dock_float(self) -> None:
        dock = self._dock
        if dock.isFloating():
            dock.setFloating(False)
            return
        size = dock.size()
        dock.setFloating(True)
        dock.resize(max(size.width(), 320), max(size.height(), 240))
        dock.show()
        QTimer.singleShot(0, lambda: self._place_floating_dock(dock))

    @staticmethod
    def _place_floating_dock(dock: QDockWidget) -> None:
        from PySide6.QtGui import QGuiApplication

        screens = QGuiApplication.screens()
        screen = screens[1] if len(screens) > 1 else QGuiApplication.primaryScreen()
        if screen is None:
            return
        geo = screen.availableGeometry()
        width = min(max(dock.width(), 320), geo.width())
        height = min(max(dock.height(), 240), geo.height())
        x = geo.x() + max(0, (geo.width() - width) // 2)
        y = geo.y() + max(0, (geo.height() - height) // 2)
        frame = clamp_frame_to_screens(QRect(x, y, width, height))
        window = dock.window()
        handle = window.windowHandle() if window is not None else None
        if handle is not None:
            handle.setScreen(screen)
        dock.setGeometry(frame)
        ensure_widget_on_screen(dock)
