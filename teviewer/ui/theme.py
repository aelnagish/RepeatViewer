"""Application-wide Qt stylesheet fragments."""

from __future__ import annotations


def application_stylesheet(
    *,
    toolbar_font_pt: int,
    browser_panel_font_pt: int = 12,
    accent: str = "#3d6ea8",
) -> str:
    size = toolbar_font_pt
    size_lg = toolbar_font_pt + 14
    panel_size = browser_panel_font_pt
    return f"""
    QMainWindow {{
        background: palette(window);
    }}
    QDockWidget {{
        font-weight: 600;
    }}
    QDockWidget::title {{
        background: palette(midlight);
        padding: 5px 8px;
        border-bottom: 1px solid palette(mid);
    }}
    QTreeWidget, QListWidget {{
        border: 1px solid palette(mid);
        border-radius: 4px;
        outline: none;
        font-size: {panel_size}pt;
    }}
    QTreeWidget::item, QListWidget::item {{
        padding: 0px 3px;
    }}
    QTreeWidget#FolderTree {{
        show-decoration-selected: 0;
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
    QTreeWidget#FolderTree::item:selected,
    QTreeWidget#FolderTree::item:hover {{
        background: transparent;
        color: palette(text);
    }}
    QListWidget::item:selected {{
        background: {accent};
        color: white;
    }}
    QListWidget::item:hover:!selected {{
        background: palette(alternate-base);
    }}
    QLineEdit#FileFilter {{
        border: 1px solid palette(mid);
        border-radius: 4px;
        padding: 2px 6px;
        min-height: 0px;
        max-height: 26px;
    }}
    QLineEdit#FileFilter:focus {{
        border-color: {accent};
    }}
    QComboBox#FileSortCombo {{
        border: 1px solid palette(mid);
        border-radius: 4px;
        padding: 2px 6px;
        min-height: 0px;
        max-height: 26px;
        combobox-popup: 0;
    }}
    QComboBox#FileSortCombo::drop-down {{
        width: 16px;
        border: none;
    }}
    QComboBox#FileSortCombo QAbstractItemView {{
        border: 1px solid palette(mid);
        selection-background-color: {accent};
        selection-color: white;
        outline: none;
        min-width: 118px;
    }}
    QComboBox#FileSortCombo QAbstractItemView::item {{
        min-height: 22px;
        padding: 2px 6px;
    }}
    QPushButton#FileRefreshButton {{
        padding: 2px 10px;
        min-height: 0px;
        max-height: 26px;
    }}
    QPushButton {{
        padding: 4px 10px;
        border-radius: 4px;
    }}
    QTabWidget::pane {{
        border: 1px solid palette(mid);
        top: -1px;
    }}
    QTabBar::tab {{
        padding: 6px 12px;
    }}
    QTabBar::tab:selected {{
        font-weight: 600;
        border-bottom: 2px solid {accent};
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
        border-radius: 4px;
        text-align: center;
    }}
    QToolButton#ToolbarOverflowButton {{
        font-size: {size_lg}pt;
        font-weight: 700;
        min-width: 22px;
        padding: 0px 8px;
        margin-left: 2px;
    }}
    QToolButton#ToolbarOverflowButton::menu-indicator {{
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
    """
