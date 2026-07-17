"""Application-wide Qt stylesheet fragments."""

from __future__ import annotations

from pathlib import Path


def _icon_url(name: str) -> str:
    """Absolute filesystem path for Qt stylesheet ``image: url(...)``."""
    path = Path(__file__).resolve().parent.parent / "assets" / "icons" / name
    # Qt stylesheets expect a plain path (not file://); quote for spaces.
    return str(path).replace("\\", "/")



def application_stylesheet(
    *,
    toolbar_font_pt: int,
    browser_panel_font_pt: int = 12,
    accent: str = "#3d6ea8",
) -> str:
    size = toolbar_font_pt
    size_lg = toolbar_font_pt + 14
    panel_size = browser_panel_font_pt
    radius = 8
    radius_sm = 6
    spin_up = _icon_url("spinbox-up.svg")
    spin_down = _icon_url("spinbox-down.svg")
    return f"""
    QMainWindow {{
        background: palette(window);
    }}
    /* Thin resize handles — no wide gutter between panels. */
    QMainWindow::separator {{
        width: 3px;
        height: 3px;
        background: palette(window);
        margin: 0px;
        padding: 0px;
    }}
    QMainWindow::separator:hover {{
        background: {accent};
    }}
    QWidget#BrowserCentral {{
        background: palette(window);
    }}
    QDockWidget {{
        font-weight: 600;
        border: 1px solid palette(mid);
        border-radius: {radius}px;
        background: palette(base);
        margin: 2px;
        titlebar-close-icon: none;
        titlebar-normal-icon: none;
    }}
    QDockWidget::title {{
        background: palette(midlight);
        padding: 5px 8px;
        border-bottom: 1px solid palette(mid);
        border-top-left-radius: {radius}px;
        border-top-right-radius: {radius}px;
    }}
    QWidget#BrowserDockTitleBar {{
        background: palette(midlight);
        border: none;
        border-bottom: 1px solid palette(mid);
        border-top-left-radius: {radius}px;
        border-top-right-radius: {radius}px;
        min-height: 24px;
        max-height: 24px;
    }}
    QWidget#FoldersPanel,
    QWidget#FilesPanel,
    QWidget#PreviewHost {{
        background: palette(base);
        border: none;
        border-bottom-left-radius: {radius}px;
        border-bottom-right-radius: {radius}px;
    }}
    QLabel#BrowserDockTitleLabel {{
        font-weight: 600;
        font-size: {panel_size}pt;
    }}
    QToolButton#BrowserDockChromeButton {{
        padding: 0px;
        margin: 0px;
        border: 1px solid transparent;
        background: transparent;
        color: palette(text);
        border-radius: {radius_sm}px;
        min-width: 18px;
        max-width: 18px;
        min-height: 18px;
        max-height: 18px;
    }}
    QToolButton#BrowserDockChromeButton:hover {{
        background: palette(alternate-base);
        border-color: palette(mid);
    }}
    QLineEdit#FileFilter {{
        border: 1px solid palette(mid);
        padding: 0px 6px;
        min-height: 16px;
        max-height: 18px;
        border-radius: {radius_sm}px;
        background: palette(base);
    }}
    QLineEdit#FileFilter:focus {{
        border-color: {accent};
    }}
    QWidget#FolderNavigationBar QToolButton {{
        border: none;
        background: transparent;
        border-radius: {radius_sm}px;
        color: palette(text);
        font-size: 11px;
    }}
    QWidget#FolderNavigationBar QToolButton:hover {{
        background: palette(alternate-base);
    }}
    QWidget#FolderNavigationBar QToolButton:disabled {{
        color: palette(mid);
    }}
    QTreeWidget, QListWidget, QTableWidget {{
        border: none;
        border-radius: 0px;
        outline: none;
        font-size: {panel_size}pt;
        background: palette(base);
    }}
    QHeaderView::section {{
        border: none;
        border-right: 1px solid palette(mid);
        border-bottom: 1px solid palette(mid);
        padding: 1px 4px;
        background: palette(midlight);
        border-radius: 0px;
        min-height: 16px;
        max-height: 18px;
    }}
    QTableWidget#FileList QHeaderView::section {{
        padding: 0px 4px;
        min-height: 16px;
        max-height: 18px;
    }}
    QTreeWidget::item, QListWidget::item {{
        padding: 0px 3px;
    }}
    QTreeWidget#FolderTree {{
        show-decoration-selected: 0;
        border: none;
        padding-top: 2px;
        background: palette(base);
    }}
    QTreeWidget#FolderTree::item {{
        padding: 1px 3px;
        margin: 0px;
        min-height: 18px;
    }}
    QTreeWidget#FolderTree::branch {{
        background: transparent;
        border: none;
        width: 0px;
        margin: 0px;
        image: none;
    }}
    QTreeWidget#FolderTree::branch:has-children:!has-siblings,
    QTreeWidget#FolderTree::branch:closed:has-children:has-siblings,
    QTreeWidget#FolderTree::branch:open:has-children:has-siblings {{
        image: none;
        border: none;
        background: transparent;
    }}
    QTreeWidget#FolderTree::item:selected {{
        background: transparent;
        color: white;
    }}
    QTreeWidget#FolderTree::item:hover:!selected {{
        background: transparent;
        color: palette(text);
    }}
    QListWidget::item:selected {{
        background: {accent};
        color: white;
        border-radius: {radius_sm}px;
    }}
    QListWidget::item:hover:!selected {{
        background: palette(alternate-base);
        border-radius: {radius_sm}px;
    }}
    QPushButton {{
        padding: 4px 10px;
        border-radius: {radius_sm}px;
    }}
    QComboBox {{
        border-radius: {radius_sm}px;
        padding: 2px 6px;
    }}
    QComboBox::drop-down {{
        border: none;
        border-top-right-radius: {radius_sm}px;
        border-bottom-right-radius: {radius_sm}px;
    }}
    QLineEdit, QSpinBox, QDoubleSpinBox {{
        border-radius: {radius_sm}px;
        background: palette(base);
    }}
    QSpinBox, QDoubleSpinBox {{
        border: 1px solid palette(mid);
        padding: 0px 4px;
        min-height: 18px;
        max-height: 20px;
    }}
    QSpinBox::up-button, QSpinBox::down-button, QDoubleSpinBox::up-button, QDoubleSpinBox::down-button {{
        border: none;
        background: palette(base);
        width: 14px;
        margin: 0px;
    }}
    QSpinBox::up-button, QDoubleSpinBox::up-button {{
        border-left: 1px solid palette(mid);
        border-top-right-radius: {radius_sm}px;
    }}
    QSpinBox::down-button, QDoubleSpinBox::down-button {{
        border-left: 1px solid palette(mid);
        border-bottom-right-radius: {radius_sm}px;
    }}
    QSpinBox::up-arrow, QDoubleSpinBox::up-arrow {{
        image: url({spin_up});
        width: 8px;
        height: 5px;
    }}
    QSpinBox::down-arrow, QDoubleSpinBox::down-arrow {{
        image: url({spin_down});
        width: 8px;
        height: 5px;
    }}
    QSpinBox::up-button:hover, QSpinBox::down-button:hover, QDoubleSpinBox::up-button:hover, QDoubleSpinBox::down-button:hover {{
        background: palette(alternate-base);
    }}
    QWidget#PreviewHeaderBar {{
        min-height: 24px;
        max-height: 24px;
        background: palette(base);
    }}
    QWidget#FastaTabStrip {{
        min-height: 22px;
        max-height: 22px;
    }}
    QPushButton#FastaPreviewTab {{
        padding: 0px 8px;
        margin: 0px;
        min-height: 20px;
        max-height: 22px;
        border: none;
        border-radius: 0px;
        background: transparent;
    }}
    QPushButton#FastaPreviewTab:checked {{
        font-weight: 600;
        border-bottom: 2px solid {accent};
    }}
    QFrame#PreviewHeaderDivider {{
        max-height: 1px;
        min-height: 1px;
        background: palette(mid);
    }}
    QWidget#AnnotationMapToolbar {{
        min-height: 26px;
        max-height: 28px;
    }}
    QWidget#AnnotationMapToolbar QSpinBox {{
        min-height: 18px;
        max-height: 18px;
    }}
    QWidget#AnnotationMapToolbar QPushButton {{
        padding: 1px 8px;
        min-height: 18px;
        max-height: 20px;
    }}
    QWidget#AnnotationMapToolbar QLabel {{
        padding: 0px;
    }}
    QTabWidget::pane {{
        border: none;
        border-radius: 0px;
        top: -1px;
    }}
    QTabBar::tab {{
        padding: 6px 12px;
        border-top-left-radius: {radius_sm}px;
        border-top-right-radius: {radius_sm}px;
    }}
    QTabBar::tab:selected {{
        font-weight: 600;
        border-bottom: 2px solid {accent};
    }}
    QTabBar#PreviewFileTabBar {{
        background: transparent;
        qproperty-drawBase: 0;
        qproperty-expanding: false;
        qproperty-usesScrollButtons: true;
        min-width: 0px;
    }}
    QTabBar#PreviewFileTabBar::tab {{
        padding: 0px 4px;
        margin: 0px 1px;
        border-radius: {radius_sm}px;
        border-bottom: none;
        min-width: 28px;
        max-width: 200px;
        max-height: 18px;
    }}
    QTabBar#PreviewFileTabBar::scroller {{
        width: 16px;
    }}
    QWidget#PreviewDockFileTabs {{
        min-width: 0px;
        background: transparent;
    }}
    QWidget#PreviewDockZoomHost {{
        background: transparent;
    }}
    QWidget#PreviewZoomBar {{
        background: transparent;
    }}
    QPushButton#PreviewZoomReset {{
        border: 1px solid palette(mid);
        background: palette(base);
        padding: 0px 8px;
        min-height: 18px;
        max-height: 18px;
        border-radius: {radius_sm}px;
    }}
    QPushButton#PreviewZoomReset:hover {{
        background: palette(alternate-base);
    }}
    QPushButton#PreviewZoomReset:pressed {{
        background: palette(mid);
    }}
    QTabBar#PreviewFileTabBar::tab:selected {{
        background: palette(highlight);
        color: palette(highlighted-text);
        font-weight: 600;
        border-bottom: none;
    }}
    QTabBar#PreviewFileTabBar::close-button {{
        margin: 1px;
    }}
    QWidget#PreviewTabs::pane {{
        border: none;
        top: 0px;
        margin: 0px;
        padding: 0px;
    }}
    QToolBar#FunctionsToolBar {{
        spacing: 1px;
        padding: 2px 6px;
        border-bottom: 1px solid palette(mid);
        background: palette(base);
    }}
    QToolBar#FunctionsToolBar QToolButton {{
        font-size: {size}pt;
        font-weight: 600;
        padding: 2px 8px;
        margin: 0px;
        border-radius: {radius_sm}px;
        text-align: center;
    }}
    QToolButton#ToolbarOverflowButton {{
        font-size: {size_lg}pt;
        font-weight: 700;
        min-width: 22px;
        padding: 0px 8px;
        margin-left: 2px;
        border-radius: {radius_sm}px;
    }}
    QToolButton#ToolbarOverflowButton::menu-indicator {{
        image: none;
        width: 0px;
    }}
    QToolButton#FileColumnsButton::menu-indicator {{
        image: none;
        width: 0px;
    }}
    QToolBar#FunctionsToolBar QToolButton:hover {{
        background: palette(alternate-base);
    }}
    QToolBar#FunctionsToolBar QToolButton:pressed {{
        background: palette(mid);
    }}
    QStatusBar {{
        border-top: 1px solid palette(mid);
    }}
    /* Scrollbars: thin rounded thumb, no visible arrows/buttons. */
    QScrollBar:vertical {{
        border: none;
        background: transparent;
        margin: 0px;
        width: 8px;
    }}
    QScrollBar:horizontal {{
        border: none;
        background: transparent;
        margin: 0px;
        height: 8px;
    }}
    QScrollBar::sub-page:horizontal, QScrollBar::add-page:horizontal {{
        background: transparent;
    }}
    QScrollBar::sub-page:vertical, QScrollBar::add-page:vertical {{
        background: transparent;
    }}
    QScrollBar::handle:horizontal {{
        background: palette(mid);
        border-radius: 4px;
        min-width: 16px;
        height: 6px;
        margin: 1px 2px;
    }}
    QScrollBar::handle:vertical {{
        background: palette(mid);
        border-radius: 4px;
        min-height: 16px;
        width: 6px;
        margin: 2px 1px;
    }}
    QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
        width: 0px;
        height: 0px;
        border: none;
    }}
    QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
        width: 0px;
        height: 0px;
        border: none;
    }}
    QScrollBar::up-arrow:horizontal, QScrollBar::down-arrow:horizontal {{
        width: 0px;
        height: 0px;
    }}
    QScrollBar::left-arrow:vertical, QScrollBar::right-arrow:vertical {{
        width: 0px;
        height: 0px;
    }}
    QSpinBox#previewZoomSpin {{
        min-width: 64px;
        max-width: 80px;
        padding-right: 2px;
    }}
    QSlider::groove:horizontal {{
        border: 1px solid palette(mid);
        background: palette(base);
        border-radius: 3px;
        height: 6px;
    }}
    QSlider::handle:horizontal {{
        background: {accent};
        border: 1px solid {accent};
        width: 14px;
        height: 14px;
        margin: -5px 0;
        border-radius: 7px;
    }}
    QSlider#annotationZoomSlider::groove:horizontal {{
        border: 1px solid palette(mid);
        background: palette(base);
        height: 4px;
        border-radius: 2px;
    }}
    QSlider#annotationZoomSlider::sub-page:horizontal {{
        background: {accent};
        border: 1px solid {accent};
        height: 4px;
        border-radius: 2px;
    }}
    QSlider#annotationZoomSlider::handle:horizontal {{
        background: {accent};
        border: 1px solid palette(text);
        width: 10px;
        height: 10px;
        margin: -4px 0;
        border-radius: 5px;
    }}
    """
