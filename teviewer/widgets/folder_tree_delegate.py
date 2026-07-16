"""Paint folder-tree selection/hover from the icon onward (not the indent gutter)."""

from __future__ import annotations

from PySide6.QtCore import Qt, QRect, QModelIndex
from PySide6.QtGui import QPainter, QPalette
from PySide6.QtWidgets import (
    QApplication,
    QStyle,
    QStyleOptionViewItem,
    QStyledItemDelegate,
    QTreeWidget,
)



class FolderTreeItemDelegate(QStyledItemDelegate):
    """Clip row highlight so indented subfolders do not show blue before the chevron."""

    def _icon_left(self, opt: QStyleOptionViewItem, widget) -> int:
        style = widget.style() if widget is not None else QApplication.style()
        deco_rect = style.subElementRect(
            QStyle.SubElement.SE_ItemViewItemDecoration,
            opt,
            widget,
        )
        if deco_rect.isValid() and deco_rect.width() > 0:
            return deco_rect.left()
        return opt.rect.left()

    def _highlight_rect(self, opt: QStyleOptionViewItem, widget) -> QRect:
        rect = opt.rect
        rect.setLeft(self._icon_left(opt, widget))
        return rect

    def _current_index(self, widget, index: QModelIndex) -> bool:
        tree = widget
        if not isinstance(tree, QTreeWidget):
            tree = widget.parent() if widget is not None else None
        if not isinstance(tree, QTreeWidget):
            return False
        current = tree.currentItem()
        if current is None:
            return False
        return index == tree.indexFromItem(current)

    def _readable_text_brush(self, opt: QStyleOptionViewItem):
        text = opt.palette.color(QPalette.ColorRole.Text)
        if not text.isValid() or text.alpha() == 0:
            text = opt.palette.color(QPalette.ColorRole.WindowText)
        return text

    def paint(self, painter: QPainter, option: QStyleOptionViewItem, index) -> None:  # noqa: N802
        opt = QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        widget = opt.widget

        is_current = self._current_index(widget, index)
        is_hover = bool(opt.state & QStyle.StateFlag.State_MouseOver)

        # Prevent Qt/macOS from painting range/multi-select chrome on other rows.
        opt.state &= ~QStyle.StateFlag.State_Selected
        opt.state &= ~QStyle.StateFlag.State_HasFocus
        opt.state &= ~QStyle.StateFlag.State_MouseOver

        if is_current or is_hover:
            bg_rect = self._highlight_rect(opt, widget)
            painter.setPen(Qt.PenStyle.NoPen)
            if is_current:
                # Selection colour follows the Preferences accent (qdarktheme maps
                # the accent onto the palette Highlight role).
                painter.setBrush(opt.palette.color(QPalette.ColorRole.Highlight))
                painter.drawRoundedRect(bg_rect, 4, 4)
                sel_text = opt.palette.color(QPalette.ColorRole.HighlightedText)
                opt.palette.setBrush(QPalette.ColorRole.Text, sel_text)
                opt.palette.setBrush(QPalette.ColorRole.WindowText, sel_text)
            else:
                painter.setBrush(opt.palette.color(QPalette.ColorRole.AlternateBase))
                painter.drawRoundedRect(bg_rect, 4, 4)
                text = self._readable_text_brush(opt)
                opt.palette.setBrush(QPalette.ColorRole.Text, text)
                opt.palette.setBrush(QPalette.ColorRole.WindowText, text)

        super().paint(painter, opt, index)
