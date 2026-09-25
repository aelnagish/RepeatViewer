"""Start-up banner (splash screen), in the spirit of Geneious.

Shows the RepeatViewer artwork, the version, what is loading and one tip;
"Next tip" shows another. A click anywhere else closes it. It stays up at
least ``MIN_SHOW_MS`` so the tip can be read, and can be turned off in
Settings → Startup. Tips come from ``assets/icons/splash/tips.txt``.
"""

from __future__ import annotations

import random
import time
from pathlib import Path

from PySide6.QtCore import QRectF, Qt, QTimer
from PySide6.QtGui import QColor, QFont, QFontMetricsF, QGuiApplication, QPainter, QPixmap
from PySide6.QtSvg import QSvgRenderer
from PySide6.QtWidgets import QApplication, QSplashScreen, QWidget

SPLASH_DIR = Path(__file__).resolve().parent.parent / "assets" / "icons" / "splash"
BANNER = SPLASH_DIR / "banner.svg"
TIPS = SPLASH_DIR / "tips.txt"
WIDTH, HEIGHT = 640, 330
MARGIN = 36
MIN_SHOW_MS = 1500
_FALLBACK_TIP = "Every \"?\" button explains the option next to it."


def short_version(version: str) -> str:
    """Major.minor for the banner: "1.0.0" → "1.0", "1.2.3" → "1.2"."""
    parts = version.split(".")
    return ".".join(parts[:2]) if len(parts) >= 2 else version


def load_tips(path: Path = TIPS) -> list[str]:
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except OSError:
        return [_FALLBACK_TIP]
    tips = [ln.strip() for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
    return tips or [_FALLBACK_TIP]


def _banner_pixmap() -> QPixmap:
    screen = QGuiApplication.primaryScreen()
    ratio = screen.devicePixelRatio() if screen is not None else 1.0
    pix = QPixmap(int(WIDTH * ratio), int(HEIGHT * ratio))
    pix.fill(QColor("#1560b0"))
    renderer = QSvgRenderer(str(BANNER))
    if renderer.isValid():
        p = QPainter(pix)
        renderer.render(p)
        p.end()
    pix.setDevicePixelRatio(ratio)
    return pix


class RepeatViewerSplash(QSplashScreen):
    def __init__(self, version: str, tips: list[str] | None = None,
                 *, seed: int | None = None) -> None:
        super().__init__(_banner_pixmap())
        self._version = short_version(version)
        self._tips = tips or load_tips()
        self._rng = random.Random(seed)
        self._tip = self._rng.randrange(len(self._tips))
        self._status = ""
        self._shown_at = time.monotonic()
        self._next_rect = QRectF()

    # ------------------------------------------------------------ state

    @property
    def tip(self) -> str:
        return self._tips[self._tip]

    def next_tip(self) -> None:
        self._tip = (self._tip + 1) % len(self._tips)
        self.repaint()

    def set_status(self, text: str) -> None:
        """Show what is loading; paints at once (the event loop may be busy)."""
        self._status = text
        self.repaint()
        QApplication.processEvents()

    def showEvent(self, event) -> None:  # noqa: N802 (Qt API)
        self._shown_at = time.monotonic()
        super().showEvent(event)

    def finish_after_minimum(self, window: QWidget, min_ms: int = MIN_SHOW_MS) -> None:
        """Close when ``window`` is up, but not before ``min_ms`` on screen."""
        left = int(min_ms - (time.monotonic() - self._shown_at) * 1000)
        if left <= 0:
            self.finish(window)
        else:
            QTimer.singleShot(left, lambda: self.finish(window))

    # ------------------------------------------------------------ input

    def mousePressEvent(self, event) -> None:  # noqa: N802 (Qt API)
        if self._next_rect.contains(event.position()):
            self.next_tip()
            event.accept()
            return
        super().mousePressEvent(event)      # QSplashScreen: a click closes it

    # ------------------------------------------------------------ painting

    def drawContents(self, painter: QPainter) -> None:  # noqa: N802 (Qt API)
        painter.setRenderHint(QPainter.RenderHint.TextAntialiasing, True)
        white = QColor("#ffffff")
        soft = QColor(255, 255, 255, 190)
        base = QFont(self.font())

        version = QFont(base)
        version.setPointSizeF(44)
        version.setWeight(QFont.Weight.Light)
        painter.setFont(version)
        painter.setPen(white)
        painter.drawText(QRectF(0, 30, WIDTH - MARGIN, 70),
                         int(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignTop),
                         self._version)

        # "↻ Next tip" (clickable), then the tip in up to two lines
        head = QFont(base)
        head.setPointSizeF(10.5)
        head.setBold(True)
        painter.setFont(head)
        painter.setPen(white)
        text = "↻  Next tip"
        fm = QFontMetricsF(head)
        self._next_rect = QRectF(MARGIN - 4, 244, fm.horizontalAdvance(text) + 8, fm.height() + 6)
        painter.drawText(QRectF(MARGIN, 247, 300, fm.height()), int(Qt.AlignmentFlag.AlignLeft), text)
        body = QFont(base)
        body.setPointSizeF(10.5)
        painter.setFont(body)
        painter.drawText(QRectF(MARGIN, 268, WIDTH - 2 * MARGIN, 36),
                         int(Qt.AlignmentFlag.AlignLeft | Qt.TextFlag.TextWordWrap), self.tip)

        small = QFont(base)
        small.setPointSizeF(9)
        painter.setFont(small)
        painter.setPen(soft)
        painter.drawText(QRectF(MARGIN, HEIGHT - 24, WIDTH - 2 * MARGIN, 16),
                         int(Qt.AlignmentFlag.AlignLeft), self._status)
