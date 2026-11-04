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
    check = _icon("check.svg", _CHECK_SVG)
    chevron = _icon("chevron.svg", _CHEVRON_SVG)
    check_rule = f"image: url({check});" if check else ""
    chevron_rule = f"image: url({chevron});" if chevron else ""
    return f"""
#cytoWorkspace {{ background: {BAND}; }}
#cytoChrome {{ background: {CHROME}; border-bottom: 1px solid {LINE}; }}
#cytoWorkspace QLabel {{ color: {INK}; background: transparent; }}
#cytoAside {{ background: {ASIDE}; }}
#cytoMain {{ background: #FFFFFF; }}
#cytoPage {{ background: {PAGE_BG}; }}
QLabel#cytoMuted {{ color: {MUTED}; }}
QLabel#cytoCrumb {{ color: {MUTED}; }}

QTabWidget::pane {{ border: none; border-top: 1px solid {LINE}; top: -1px; }}
QTabWidget::tab-bar {{ left: 12px; }}
QTabBar {{ background: {BAND}; }}
QTabBar::tab {{ padding: 8px 16px; margin-top: 6px; margin-right: 2px; color: {SOFT}; background: transparent;
               border: 1px solid transparent; border-bottom: none;
               border-top-left-radius: 6px; border-top-right-radius: 6px; }}
QTabBar::tab:selected {{ background: #FFFFFF; color: {INK}; font-weight: 600; border-color: {LINE}; }}
QTabBar::tab:hover:!selected {{ color: {INK}; }}

QPushButton {{ background: #FFFFFF; color: {INK}; border: 1px solid {FIELD_LINE}; border-radius: 6px; padding: 5px 12px; }}
QPushButton:hover {{ background: #F3F5F9; }}
QPushButton:pressed {{ background: #E8EBF2; }}
QPushButton:checked {{ background: {SELECTED}; border-color: {BLUE}; }}
QPushButton:disabled {{ color: #A3A6AB; border-color: #DADADF; background: #FAFAFB; }}
QPushButton::menu-indicator {{ width: 0px; }}
QPushButton#cytoPrimary {{ background: {BLUE}; color: #FFFFFF; border: 1px solid {BLUE}; padding: 6px 14px; }}
QPushButton#cytoPrimary:hover {{ background: {BLUE_HOVER}; border-color: {BLUE_HOVER}; }}
QPushButton#cytoLink {{ background: transparent; border: none; color: {BLUE}; padding: 0px; text-align: left; }}
QPushButton#cytoLink:hover {{ text-decoration: underline; }}

QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox, QTextEdit, QPlainTextEdit {{
    background: #FFFFFF; color: {INK}; border: 1px solid {FIELD_LINE}; border-radius: 6px;
    selection-background-color: {SELECTED}; selection-color: {INK}; }}
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{ padding: 4px 8px; min-height: 18px; }}
QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus, QTextEdit:focus {{ border-color: {BLUE}; }}
QComboBox::drop-down {{ border: none; width: 22px; }}
QComboBox::down-arrow {{ {chevron_rule} width: 10px; height: 10px; }}
QComboBox QAbstractItemView {{ background: #FFFFFF; color: {INK}; border: 1px solid {LINE};
    selection-background-color: {SELECTED}; selection-color: {INK}; outline: 0; }}

QCheckBox {{ color: {INK}; spacing: 8px; background: transparent; }}
QCheckBox::indicator {{ width: 15px; height: 15px; border: 1px solid {FIELD_LINE}; border-radius: 4px; background: #FFFFFF; }}
QCheckBox::indicator:checked {{ background: {TICK}; border-color: {TICK}; {check_rule} }}
QCheckBox::indicator:disabled {{ background: #F0F0F2; }}

QSlider::groove:horizontal {{ height: 4px; background: {LINE}; border-radius: 2px; }}
QSlider::sub-page:horizontal {{ background: {TICK}; border-radius: 2px; }}
QSlider::handle:horizontal {{ background: {TICK}; width: 16px; height: 16px; margin: -6px 0; border-radius: 8px; }}

QSplitter::handle {{ background: {LINE_LIGHT}; }}
QScrollArea {{ border: none; background: transparent; }}

QTableWidget {{ background: #FFFFFF; color: {INK}; gridline-color: {LINE_LIGHT}; border: 1px solid {LINE_LIGHT};
               selection-background-color: {SELECTED}; selection-color: {INK}; }}
QHeaderView::section {{ background: {ASIDE}; color: {SOFT}; border: none; border-bottom: 1px solid {LINE}; padding: 5px 6px; }}
QProgressBar {{ border: 1px solid {LINE}; border-radius: 4px; background: #FFFFFF; height: 8px; text-align: center; }}
QProgressBar::chunk {{ background: {TICK}; border-radius: 3px; }}
""" + spin_css()


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


_CHEVRON_UP_SVG = ('<svg xmlns="http://www.w3.org/2000/svg" width="8" height="8" viewBox="0 0 10 10">'
                   '<path d="M2 6.5l3-3 3 3" fill="none" stroke="#45484D" stroke-width="1.5" '
                   'stroke-linecap="round" stroke-linejoin="round"/></svg>')


def spin_css() -> str:
    up = _icon("chevron_up.svg", _CHEVRON_UP_SVG)
    down = _icon("chevron.svg", _CHEVRON_SVG)
    return f"""
QAbstractSpinBox {{ padding-right: 20px; }}
QAbstractSpinBox::up-button {{ subcontrol-origin: border; subcontrol-position: top right; width: 18px; border: none;
    border-left: 1px solid #E0E0E4; border-top-right-radius: 6px; background: #FFFFFF; }}
QAbstractSpinBox::down-button {{ subcontrol-origin: border; subcontrol-position: bottom right; width: 18px; border: none;
    border-left: 1px solid #E0E0E4; border-bottom-right-radius: 6px; background: #FFFFFF; }}
QAbstractSpinBox::up-arrow {{ {f"image: url({up});" if up else ""} width: 8px; height: 8px; }}
QAbstractSpinBox::down-arrow {{ {f"image: url({down});" if down else ""} width: 8px; height: 8px; }}
"""
