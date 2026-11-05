"""Look of the Cytogenetics workspace, taken from the agreed mockup (colours, fields, buttons, tabs).

One stylesheet is set on the workspace; it reaches every tab and the dialogs opened from it.
"""

from __future__ import annotations

import tempfile
from pathlib import Path

INK = "#1C1C1E"
SOFT = "#45484D"
MUTED = "#5A5D63"
LINE = "#D0D0D5"
LINE_LIGHT = "#E0E0E4"
FIELD_LINE = "#C4C4CA"
CHROME = "#F4F4F6"
BAND = "#E9E9EB"
ASIDE = "#F7F7F9"
PAGE_BG = "#DADBDF"
BLUE = "#1F4E9A"
BLUE_HOVER = "#163A75"
SELECTED = "#E3EBF8"
TICK = "#2F6FDB"

_CHECK_SVG = ('<svg xmlns="http://www.w3.org/2000/svg" width="16" height="16" viewBox="0 0 16 16">'
              '<path d="M3.6 8.4l2.9 2.9 5.9-6.6" fill="none" stroke="#ffffff" stroke-width="2.1" '
              'stroke-linecap="round" stroke-linejoin="round"/></svg>')
_CHEVRON_SVG = ('<svg xmlns="http://www.w3.org/2000/svg" width="10" height="10" viewBox="0 0 10 10">'
                '<path d="M2 3.5l3 3 3-3" fill="none" stroke="#45484D" stroke-width="1.5" '
                'stroke-linecap="round" stroke-linejoin="round"/></svg>')


def _icon(name: str, svg: str) -> str:
    folder = Path(tempfile.gettempdir()) / "repeatviewer_cyto_icons"
    path = folder / name
    try:
        folder.mkdir(parents=True, exist_ok=True)
        if not path.is_file() or path.read_text(encoding="utf-8") != svg:
            path.write_text(svg, encoding="utf-8")
    except OSError:
        return ""
    return path.as_posix()


def workspace_css() -> str:
    """Only backgrounds of the side columns and a few text colours: every control stays the app's native one."""
    return """
#cytoAside { background: palette(window); }
QTabWidget#cytoTabs { background: palette(window); }
#cytoMain { background: palette(base); }
#cytoPage { background: palette(base); }
QLabel#cytoMuted, QLabel#cytoCrumb { color: gray; }
QPushButton#cytoLink { background: transparent; border: none; color: palette(link); padding: 0px; text-align: left; }
QPushButton#cytoLink:hover { text-decoration: underline; }
QScrollArea { border: none; background: transparent; }
"""


def app_slider(slider) -> None:
    """The app's one slider style (UI rule 12): blue track left of a round blue handle."""
    from repeatviewer.widgets.preview_controls import apply_annotation_zoom_slider_style

    slider.setObjectName("annotationZoomSlider")
    apply_annotation_zoom_slider_style(slider)


def link_button(text: str, parent=None):
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QPushButton

    b = QPushButton(text, parent)
    b.setObjectName("cytoLink")
    b.setFlat(True)
    b.setCursor(Qt.CursorShape.PointingHandCursor)
    return b


def muted(label):
    label.setObjectName("cytoMuted")
    return label


def clear_layout(layout) -> None:
    """Empty a layout at once: widgets are hidden and unparented now, deleted later; nested layouts too.

    Only calling deleteLater() leaves the old widgets painted until the event loop gets to them,
    and widgets inside nested layouts were never removed at all (they showed as ghost text).
    """
    if layout is None:
        return
    while layout.count():
        item = layout.takeAt(0)
        widget = item.widget()
        if widget is not None:
            widget.hide()
            widget.setParent(None)
            widget.deleteLater()
            continue
        child = item.layout()
        if child is not None:
            clear_layout(child)
            child.setParent(None)
            child.deleteLater()


BOX_CSS = "QFrame#cytoBox { border: 1px solid palette(mid); border-radius: 8px; background: palette(base); }"


def rounded_box(inner, parent=None, margin: int = 2):
    """Put a list, table or canvas in the app's round-edged box (the same frame as the Files and Preview panels)."""
    from PySide6.QtWidgets import QFrame, QVBoxLayout

    box = QFrame(parent)
    box.setObjectName("cytoBox")
    box.setStyleSheet(BOX_CSS)
    lay = QVBoxLayout(box)
    lay.setContentsMargins(margin, margin, margin, margin)
    lay.setSpacing(0)
    if isinstance(inner, QFrame):
        inner.setFrameShape(QFrame.Shape.NoFrame)
    lay.addWidget(inner)
    return box


def form_layout(parent=None):
    """The form of the app's tool windows: labels right-aligned in one column, fields growing."""
    from PySide6.QtCore import Qt
    from PySide6.QtWidgets import QFormLayout

    form = QFormLayout(parent)
    form.setContentsMargins(0, 0, 0, 0)
    form.setLabelAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
    form.setFieldGrowthPolicy(QFormLayout.FieldGrowthPolicy.ExpandingFieldsGrow)
    form.setHorizontalSpacing(12)
    form.setVerticalSpacing(6)
    return form


def section_label(text: str, parent=None):
    from PySide6.QtWidgets import QLabel

    return QLabel(f"<b>{text}</b>", parent)


def fit_side_scroll(scroll) -> None:
    """Keep a side column at least as wide as its content, so nothing is cut off at the right.

    The panels have no sideways scroll bar; their content width depends on the platform's fonts (macOS is wider
    than Linux), so the minimum width is measured whenever the content changes, not fixed in pixels.
    """
    from PySide6.QtCore import QEvent, QObject, QTimer

    class _Fit(QObject):
        def __init__(self, area) -> None:
            super().__init__(area)
            self._area = area
            self._queued = False
            if area.widget() is not None:
                area.widget().installEventFilter(self)
            self.queue()

        def queue(self) -> None:
            if not self._queued:
                self._queued = True
                QTimer.singleShot(0, self.fit)

        def eventFilter(self, obj, event) -> bool:  # noqa: N802
            if event.type() in (QEvent.Type.LayoutRequest, QEvent.Type.Show, QEvent.Type.FontChange):
                self.queue()
            return False

        def fit(self) -> None:
            self._queued = False
            area = self._area
            body = area.widget()
            if body is None:
                return
            need = (body.minimumSizeHint().width() + area.verticalScrollBar().sizeHint().width()
                    + 2 * area.frameWidth() + 2)
            if area.minimumWidth() < need:
                area.setMinimumWidth(need)
            w = area
            for _ in range(3):  # a stacked widget or frame around the column must not cap it either
                if w is None:
                    break
                if w.maximumWidth() < need:
                    w.setMaximumWidth(need)
                w = w.parentWidget()

    scroll._fit_side = _Fit(scroll)
