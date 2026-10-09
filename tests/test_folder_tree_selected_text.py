"""Selected folder rows must use high-contrast text over custom folder colors."""

from __future__ import annotations

from PySide6.QtCore import QRect, QSize, Qt
from PySide6.QtGui import QBrush, QColor, QImage, QPainter, QPalette
from PySide6.QtWidgets import QStyle, QStyleOptionViewItem, QTreeWidget, QTreeWidgetItem

from repeatviewer.ui.theme import application_stylesheet
from repeatviewer.widgets.folder_tree_delegate import (
    FolderTreeItemDelegate,
    _contrasting_on,
)


def _paint_selected_row(
    *,
    label: str,
    folder_color: QColor,
    with_app_stylesheet: bool,
) -> tuple[QImage, QColor, QColor]:
    tree = QTreeWidget()
    tree.setObjectName("FolderTree")
    tree.setHeaderHidden(True)
    if with_app_stylesheet:
        tree.setStyleSheet(
            application_stylesheet(toolbar_font_pt=12, accent="#3b82f6")
        )
    item = QTreeWidgetItem([label])
    item.setData(0, Qt.ItemDataRole.UserRole, "/tmp/colored-folder")
    item.setForeground(0, QBrush(folder_color))
    tree.addTopLevelItem(item)
    tree.setCurrentItem(item)
    item.setSelected(True)
    tree.resize(320, 60)
    tree.setItemDelegate(FolderTreeItemDelegate(tree))

    img = QImage(320, 40, QImage.Format.Format_ARGB32)
    img.fill(QColor("#ffffff"))
    painter = QPainter(img)
    index = tree.indexFromItem(item)
    opt = QStyleOptionViewItem()
    opt.initFrom(tree)
    opt.rect = QRect(0, 0, 300, 24)
    opt.state = QStyle.StateFlag.State_Enabled | QStyle.StateFlag.State_Selected
    opt.widget = tree
    opt.palette = tree.palette()
    opt.features = QStyleOptionViewItem.ViewItemFeature.HasDisplay
    opt.text = label
    opt.decorationSize = QSize(16, 16)
    tree.itemDelegate().paint(painter, opt, index)
    painter.end()

    highlight = tree.palette().color(QPalette.ColorRole.Highlight)
    expected = _contrasting_on(highlight)
    tree.close()
    return img, expected, folder_color


def _count_near(img: QImage, target: QColor, *, skip_highlight: QColor | None) -> int:
    count = 0
    for y in range(2, 22):
        for x in range(8, 200):
            c = QColor(img.pixel(x, y))
            if c.alpha() < 200:
                continue
            if skip_highlight is not None:
                # Ignore selection-bar blues so text pixels dominate the count.
                if (
                    abs(c.red() - skip_highlight.red()) < 45
                    and abs(c.green() - skip_highlight.green()) < 45
                    and abs(c.blue() - skip_highlight.blue()) < 45
                ):
                    continue
                if c.blue() > c.red() + 40 and c.blue() > 100 and c.lightness() < 190:
                    continue
            if (
                abs(c.red() - target.red()) < 55
                and abs(c.green() - target.green()) < 55
                and abs(c.blue() - target.blue()) < 55
            ):
                count += 1
    return count


def test_selected_colored_folder_text_is_contrasting(qapp):
    custom = QColor("#1a7a4c")  # dark green — unreadable on default blue selection
    img, expected, _ = _paint_selected_row(
        label="GreenFolderName",
        folder_color=custom,
        with_app_stylesheet=False,
    )
    highlight = QColor("#308cc6")  # approximate; count uses expected vs custom
    # Prefer counting against expected contrast color vs the custom folder color.
    contrast_pixels = _count_near(img, expected, skip_highlight=None)
    custom_pixels = _count_near(img, custom, skip_highlight=highlight)
    assert contrast_pixels > custom_pixels
    assert custom_pixels == 0


def test_selected_blue_folder_text_stays_readable_with_stylesheet(qapp):
    # Same hue family as the selection bar — the worst case for contrast.
    custom = QColor("#2563eb")
    img, expected, _ = _paint_selected_row(
        label="BlueFolderName",
        folder_color=custom,
        with_app_stylesheet=True,
    )
    contrast_pixels = _count_near(img, expected, skip_highlight=custom)
    custom_pixels = _count_near(img, custom, skip_highlight=expected)
    assert contrast_pixels > 0
    assert contrast_pixels > custom_pixels
