"""Application-wide Qt stylesheet fragments."""

from __future__ import annotations


def application_stylesheet(
    *,
    toolbar_font_pt: int,
    browser_panel_font_pt: int = 12,
) -> str:
    size = toolbar_font_pt
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
    QTreeWidget::item:selected, QListWidget::item:selected {{
        background: #3d6ea8;
        color: white;
    }}
    QTreeWidget::item:hover:!selected, QListWidget::item:hover:!selected {{
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
        border-color: #3d6ea8;
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
        selection-background-color: #3d6ea8;
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
    QPushButton:hover {{
        background: palette(alternate-base);
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
        border-bottom: 2px solid #3d6ea8;
    }}
    QToolBar#FunctionsToolBar {{
        spacing: 2px;
        padding: 0px 6px;
        border-bottom: 1px solid palette(mid);
        background: palette(base);
    }}
    QToolBar#FunctionsToolBar QToolButton {{
        font-size: {size}pt;
        font-weight: 600;
        padding: 0px 8px;
        margin: 0px;
        border-radius: 3px;
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
