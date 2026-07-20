"""Paint folder-tree selection/hover from the icon onward (not the indent gutter)."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, QRect, QModelIndex
from PySide6.QtGui import QColor, QPainter, QPalette
from PySide6.QtWidgets import (
    QApplication,
    QStyle,
    QStyleOptionViewItem,
    QStyledItemDelegate,
    QTreeWidget,
)


def _contrasting_on(background: QColor) -> QColor:
    """White or near-black text, whichever reads on ``background``."""
    if not background.isValid():
        return QColor(Qt.GlobalColor.white)
    # Rec. 709 relative luminance.
    lum = (
        0.2126 * background.redF()
        + 0.7152 * background.greenF()
        + 0.0722 * background.blueF()
    )
    return QColor(Qt.GlobalColor.white) if lum < 0.62 else QColor(20, 20, 20)


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
        if current is not None and index == tree.indexFromItem(current):
            return True

        active_key = getattr(tree, "_active_folder_key", None)
        if not active_key:
            return False
        raw = index.data(Qt.ItemDataRole.UserRole)
        if not raw:
            return False
        try:
            item_key = str(Path(str(raw)).expanduser().resolve())
        except OSError:
            item_key = str(raw)
        return item_key == active_key

    def _readable_text_brush(self, opt: QStyleOptionViewItem):
        text = opt.palette.color(QPalette.ColorRole.Text)
        if not text.isValid() or text.alpha() == 0:
            text = opt.palette.color(QPalette.ColorRole.WindowText)
        return text

    def paint(self, painter: QPainter, option: QStyleOptionViewItem, index) -> None:  # noqa: N802
        opt = QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        widget = opt.widget

        # Read selection from the original option before we clear Qt's chrome flags.
        is_selected = bool(option.state & QStyle.StateFlag.State_Selected)
        is_current = self._current_index(widget, index)
        is_hover = bool(opt.state & QStyle.StateFlag.State_MouseOver)

        # Draw our own clipped highlight (icon → text), not the full-row macOS
        # selection chrome — but still honor ExtendedSelection multi-select.
        opt.state &= ~QStyle.StateFlag.State_Selected
        opt.state &= ~QStyle.StateFlag.State_HasFocus
        opt.state &= ~QStyle.StateFlag.State_MouseOver

        if is_selected or is_current or is_hover:
            bg_rect = self._highlight_rect(opt, widget)
            painter.setPen(Qt.PenStyle.NoPen)
            if is_selected or is_current:
                # Selection colour follows the Preferences accent (qdarktheme maps
                # the accent onto the palette Highlight role).
                highlight = opt.palette.color(QPalette.ColorRole.Highlight)
                painter.setBrush(highlight)
                painter.drawRoundedRect(bg_rect, 4, 4)
                # Prefer HighlightedText when it contrasts; otherwise force white/black.
                # (Custom accent + qdarktheme often leaves HighlightedText dark.)
                sel_text = opt.palette.color(QPalette.ColorRole.HighlightedText)
                if abs(sel_text.lightness() - highlight.lightness()) < 80:
                    sel_text = _contrasting_on(highlight)
                opt.palette.setBrush(QPalette.ColorRole.Text, sel_text)
                opt.palette.setBrush(QPalette.ColorRole.WindowText, sel_text)
                opt.palette.setBrush(QPalette.ColorRole.HighlightedText, sel_text)
                # Re-enable Selected so folder chevrons use the white Selected
                # icon pixmap (Normal is dark grey and disappears on blue).
                opt.state |= QStyle.StateFlag.State_Selected
            else:
                painter.setBrush(opt.palette.color(QPalette.ColorRole.AlternateBase))
                painter.drawRoundedRect(bg_rect, 4, 4)
                text = self._readable_text_brush(opt)
                opt.palette.setBrush(QPalette.ColorRole.Text, text)
                opt.palette.setBrush(QPalette.ColorRole.WindowText, text)

        super().paint(painter, opt, index)
