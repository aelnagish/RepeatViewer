"""Paint folder-tree selection/hover on the full row without shifting labels."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt, QRect, QModelIndex
from PySide6.QtGui import QColor, QIcon, QPainter, QPalette, QPen
from PySide6.QtWidgets import (
    QApplication,
    QStyle,
    QStyleOptionViewItem,
    QStyledItemDelegate,
    QTreeWidget,
)


_FOLDER_HOVER_BG_LIGHT = QColor(220, 235, 255)
_FOLDER_HOVER_BG_DARK = QColor(45, 65, 95)


def _folder_hover_background(palette: QPalette) -> QColor:
    window = palette.color(QPalette.ColorRole.Window)
    if window.isValid() and window.lightness() < 128:
        return _FOLDER_HOVER_BG_DARK
    return _FOLDER_HOVER_BG_LIGHT


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


def _guide_pen(palette: QPalette) -> QPen:
    guide = palette.color(QPalette.ColorRole.Mid)
    if not guide.isValid() or guide.alpha() == 0:
        guide = palette.color(QPalette.ColorRole.Dark)
    if not guide.isValid() or guide.alpha() == 0:
        guide = QColor(120, 120, 120)
    guide.setAlpha(160)
    return QPen(guide, 1)


def _index_depth(index: QModelIndex) -> int:
    depth = 0
    parent = index.parent()
    while parent.isValid():
        depth += 1
        parent = parent.parent()
    return depth


def _is_last_child(index: QModelIndex) -> bool:
    parent = index.parent()
    model = index.model()
    if model is None:
        return True
    if not parent.isValid():
        return index.row() >= model.rowCount(QModelIndex()) - 1
    return index.row() >= model.rowCount(parent) - 1


def _ancestor_at_depth(index: QModelIndex, target_depth: int) -> QModelIndex:
    """Return the ancestor of ``index`` whose depth equals ``target_depth``."""
    current_depth = _index_depth(index)
    node = index
    while current_depth > target_depth and node.isValid():
        node = node.parent()
        current_depth -= 1
    return node


def paint_folder_indent_guides(
    painter: QPainter,
    *,
    tree: QTreeWidget,
    index: QModelIndex,
    top: int,
    bottom: int,
    content_left: int,
    chevron_center_x: int | None = None,
) -> None:
    """Draw vertical indent guides centered on folder chevrons.

    The last child at any level gets a short horizontal end-cap (└).
    """
    depth = _index_depth(index)
    if depth <= 0:
        return
    indent = max(1, int(tree.indentation()))
    if chevron_center_x is not None:
        # Offset from this row's content left to the chevron tip / icon center.
        chevron_offset = chevron_center_x - content_left
    else:
        chevron_offset = max(1, int(tree.iconSize().width()) // 2)

    painter.save()
    painter.setClipping(False)
    painter.setPen(_guide_pen(tree.palette()))

    for column in range(depth):
        # Column N lines up with the chevron of the ancestor at depth N.
        steps_from_content = depth - column
        x = content_left - steps_from_content * indent + chevron_offset
        if x >= content_left:
            continue

        ancestor = _ancestor_at_depth(index, column + 1)
        last = _is_last_child(ancestor) if ancestor.isValid() else True
        innermost = column == depth - 1

        if innermost:
            if last:
                mid = (top + bottom) // 2
                painter.drawLine(x, top, x, mid)
                stub = max(6, min(indent // 2, 10))
                painter.drawLine(x, mid, x + stub, mid)
            else:
                painter.drawLine(x, top, x, bottom)
        elif not last:
            # Keep the stem only while that ancestor still has siblings below.
            painter.drawLine(x, top, x, bottom)

    painter.restore()


class FolderTreeItemDelegate(QStyledItemDelegate):
    """Row highlight for select/hover; label position stays aligned with siblings."""

    def _tree_widget(self, widget) -> QTreeWidget | None:
        if isinstance(widget, QTreeWidget):
            return widget
        parent = widget.parent() if widget is not None else None
        return parent if isinstance(parent, QTreeWidget) else None

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
        """Full row width so select/hover does not shift label vs siblings."""
        del widget
        return opt.rect.adjusted(1, 1, -1, -1)

    def _current_index(self, widget, index: QModelIndex) -> bool:
        tree = self._tree_widget(widget)
        if tree is None:
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

    def _is_drop_target(self, widget, index: QModelIndex) -> bool:
        tree = self._tree_widget(widget)
        if tree is None or not hasattr(tree, "drop_target_item"):
            return False
        target = tree.drop_target_item()
        if target is None:
            return False
        return index == tree.indexFromItem(target)

    def paint(self, painter: QPainter, option: QStyleOptionViewItem, index) -> None:  # noqa: N802
        opt = QStyleOptionViewItem(option)
        self.initStyleOption(opt, index)
        widget = opt.widget

        # Read selection from the original option before we clear Qt's chrome flags.
        is_selected = bool(option.state & QStyle.StateFlag.State_Selected)
        is_current = self._current_index(widget, index)
        is_hover = bool(opt.state & QStyle.StateFlag.State_MouseOver)
        is_drop_target = self._is_drop_target(widget, index)

        # Custom row background only — do not let Qt shift text/decoration on select.
        opt.state &= ~QStyle.StateFlag.State_Selected
        opt.state &= ~QStyle.StateFlag.State_HasFocus
        opt.state &= ~QStyle.StateFlag.State_MouseOver
        opt.showDecorationSelected = False

        tree = self._tree_widget(widget)
        if tree is not None:
            deco_left = self._icon_left(opt, widget)
            style = widget.style() if widget is not None else QApplication.style()
            deco = style.subElementRect(
                QStyle.SubElement.SE_ItemViewItemDecoration,
                opt,
                widget,
            )
            if deco.isValid() and deco.width() > 0:
                chevron_cx = (deco.left() + deco.right()) // 2
            else:
                chevron_cx = deco_left + max(1, int(tree.iconSize().width()) // 2)
            paint_folder_indent_guides(
                painter,
                tree=tree,
                index=index,
                top=opt.rect.top(),
                bottom=opt.rect.bottom(),
                content_left=opt.rect.left(),
                chevron_center_x=chevron_cx,
            )

        paint_selected_icon = False
        if is_drop_target:
            # Same light-blue hover as mouse-over (not selection / not green).
            bg_rect = self._highlight_rect(opt, widget)
            painter.setPen(Qt.PenStyle.NoPen)
            painter.setBrush(_folder_hover_background(opt.palette))
            painter.drawRoundedRect(bg_rect, 4, 4)
            text = self._readable_text_brush(opt)
            opt.palette.setBrush(QPalette.ColorRole.Text, text)
            opt.palette.setBrush(QPalette.ColorRole.WindowText, text)
        elif is_selected or is_current or is_hover:
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
                # Always override user folder colors (orange/green/blue…) so text
                # stays readable on the blue selection bar — same as Files.
                sel_text = opt.palette.color(QPalette.ColorRole.HighlightedText)
                if abs(sel_text.lightness() - highlight.lightness()) < 80:
                    sel_text = _contrasting_on(highlight)
                for group in (
                    QPalette.ColorGroup.Active,
                    QPalette.ColorGroup.Inactive,
                    QPalette.ColorGroup.Disabled,
                ):
                    opt.palette.setBrush(group, QPalette.ColorRole.Text, sel_text)
                    opt.palette.setBrush(group, QPalette.ColorRole.WindowText, sel_text)
                    opt.palette.setBrush(
                        group, QPalette.ColorRole.HighlightedText, sel_text
                    )
                # State_Selected was cleared so Qt would keep Normal (grey) chevrons;
                # paint the Selected (white) pixmap so the arrow matches the text.
                paint_selected_icon = not opt.icon.isNull()
            else:
                painter.setBrush(_folder_hover_background(opt.palette))
                painter.drawRoundedRect(bg_rect, 4, 4)
                text = self._readable_text_brush(opt)
                opt.palette.setBrush(QPalette.ColorRole.Text, text)
                opt.palette.setBrush(QPalette.ColorRole.WindowText, text)

        force_sel_text: QColor | None = None
        if is_selected or is_current:
            # Re-read from the palette we just wrote (may be HighlightedText or
            # the luminance fallback). Used for an explicit text draw below.
            force_sel_text = opt.palette.color(
                QPalette.ColorGroup.Active, QPalette.ColorRole.Text
            )

        if paint_selected_icon:
            style = widget.style() if widget is not None else QApplication.style()
            deco = style.subElementRect(
                QStyle.SubElement.SE_ItemViewItemDecoration,
                opt,
                widget,
            )
            if deco.isValid() and deco.width() > 0:
                opt.icon.paint(
                    painter,
                    deco,
                    Qt.AlignmentFlag.AlignCenter,
                    QIcon.Mode.Selected,
                    QIcon.State.Off,
                )
                opt.icon = QIcon()

        # Draw via the style directly — ``super().paint`` re-runs initStyleOption
        # and would restore the item's custom ForegroundRole color over sel_text.
        # When selected, also strip ``opt.text`` and draw the label ourselves so
        # QStyleSheetStyle cannot re-apply ForegroundRole (blue/green folders on
        # a blue selection bar would otherwise stay unreadable).
        style = widget.style() if widget is not None else QApplication.style()
        label = opt.text
        if force_sel_text is not None and label:
            opt.text = ""
        style.drawControl(QStyle.ControlElement.CE_ItemViewItem, opt, painter, widget)
        if force_sel_text is not None and label:
            text_rect = style.subElementRect(
                QStyle.SubElement.SE_ItemViewItemText,
                opt,
                widget,
            )
            if not text_rect.isValid() or text_rect.width() <= 0:
                text_rect = opt.rect.adjusted(2, 0, -2, 0)
            style.drawItemText(
                painter,
                text_rect,
                int(opt.displayAlignment or Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter),
                opt.palette,
                bool(opt.state & QStyle.StateFlag.State_Enabled),
                label,
                QPalette.ColorRole.Text,
            )
