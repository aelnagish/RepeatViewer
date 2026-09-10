"""Toolbar appearance, visibility, and button order."""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtCore import Qt
from PySide6.QtGui import QColor, QDrag, QDropEvent, QPainter, QPen
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QListView,
    QMenu,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.spin_boxes import NoWheelSpinBox as QSpinBox


if TYPE_CHECKING:
    from repeatviewer.plugins.registry import PluginRegistry

TOOLBAR_SETTINGS_GROUP = "toolbar"
DEFAULT_TOOLBAR_FONT_SIZE = 12
MAX_TOOLBAR_FONT_SIZE = 20
MIN_TOOLBAR_FONT_SIZE = 8

# Removed from toolbar customization (kept elsewhere in the app if needed).
TOOLBAR_CUSTOMIZE_HIDDEN = frozenset({"te_divergence_plot", "_sep_plugins"})

# Add-button menus mirror the main menu bar (toolbar-capable items only).
FILE_TOOLBAR_IDS: tuple[str, ...] = (
    "define_genome",
)

SEQUENCE_TOOLBAR_IDS: tuple[str, ...] = (
    "new_sequence",
    "group_sequences",
    "extract_sequences",
    "reverse_complement",
    "orf_prediction",
    "translate",
)

TOOLS_TOOLBAR_IDS: tuple[str, ...] = (
    "aliview",
    "blast",
    "align",
    "fastga",
    "repeat_quantify",
    "extend",
    "tsd_search",
    "teaid",
    "msa_clean",
    "consensus",
    "flexidot",
    "te_annotation",
    "te_genome_map",
    "te_benchmark",
)

# Built-in tools known to the customize dialog.
ALL_BUILTIN_TOOLBAR_IDS: tuple[str, ...] = (
    *FILE_TOOLBAR_IDS,
    *SEQUENCE_TOOLBAR_IDS,
    *TOOLS_TOOLBAR_IDS,
    "_sep1",
    "plugin:minimap2",
)

# Default visible toolbar (Settings → Toolbar → Reset to default).
DEFAULT_TOOLBAR_ORDER: list[str] = [
    "blast",
    "align",
    "flexidot",
    "teaid",
    "msa_clean",
    "extend",
    "consensus",
    "repeat_quantify",
    "te_genome_map",
]

TOOLBAR_LABELS: dict[str, str] = {
    "define_genome": "Define genome",
    "new_sequence": "New Sequence…",
    "group_sequences": "Combine Sequences to One File",
    "extract_sequences": "Extract Sequences from File",
    "reverse_complement": "Reverse Complement",
    "orf_prediction": "ORF Prediction",
    "translate": "Translate",
    "aliview": "AliView",
    "_sep1": "— separator —",
    "align": "Align",
    "fastga": "FastGA",
    "repeat_quantify": "Repeat Quantify",
    "blast": "BLAST",
    "extend": "MSA Extend",
    "teaid": "TEAid",
    "msa_clean": "MSA Clean",
    "consensus": "Consensus",
    "flexidot": "Flexidot",
    "te_annotation": "TE annotation",
    "te_genome_map": "Visualization",
    "te_benchmark": "Benchmarking",
    "te_divergence_plot": "TE divergence",
    "tsd_search": "TSD search",
    "plugin:minimap2": "Minimap2",
}


def is_toolbar_separator(tool_id: str) -> bool:
    return tool_id.startswith("_sep")


def toolbar_item_label(tool_id: str, labels: dict[str, str] | None = None) -> str:
    lookup = labels or TOOLBAR_LABELS
    if is_toolbar_separator(tool_id):
        return "— separator —"
    return lookup.get(tool_id, tool_id)


def next_separator_tool_id(existing_ids: set[str]) -> str:
    if "_sep1" not in existing_ids:
        return "_sep1"
    index = 2
    while f"_sep{index}" in existing_ids:
        index += 1
    return f"_sep{index}"


def plugin_tool_id(plugin_id: str) -> str:
    return f"plugin:{plugin_id}"


def toolbar_labels_for(registry: PluginRegistry | None = None) -> dict[str, str]:
    labels = dict(TOOLBAR_LABELS)
    if registry is not None:
        for plugin in registry.all_plugins():
            labels[plugin_tool_id(plugin.id)] = plugin.label
    return labels


def normalize_toolbar_order(
    order: list[str],
    *,
    known: set[str] | None = None,
) -> list[str]:
    """Drop removed ids and dedupe while preserving order."""
    cleaned: list[str] = []
    seen: set[str] = set()
    for tool_id in order:
        if tool_id in TOOLBAR_CUSTOMIZE_HIDDEN:
            continue
        if tool_id in seen:
            continue
        if known is not None and tool_id not in known and not is_toolbar_separator(tool_id):
            continue
        cleaned.append(tool_id)
        seen.add(tool_id)
    return cleaned


def default_toolbar_order_for(registry: PluginRegistry | None = None) -> list[str]:
    del registry  # Default set is fixed; plugins are added via Add button.
    return normalize_toolbar_order(list(DEFAULT_TOOLBAR_ORDER))


def build_toolbar_catalog(
    visible_order: list[str],
    registry: PluginRegistry | None = None,
) -> list[tuple[str, bool]]:
    """Toolbar customize list: only tools already on / listed for the toolbar."""
    labels = toolbar_labels_for(registry)
    catalog: list[tuple[str, bool]] = []
    seen: set[str] = set()
    for tool_id in normalize_toolbar_order(visible_order, known=set(labels)):
        if tool_id in seen:
            continue
        if tool_id not in labels and not is_toolbar_separator(tool_id):
            continue
        catalog.append((tool_id, True))
        seen.add(tool_id)
    return catalog


def load_toolbar_order(registry: PluginRegistry | None = None) -> list[str]:
    from PySide6.QtCore import QSettings

    labels = toolbar_labels_for(registry)
    default_order = default_toolbar_order_for(registry)
    s = QSettings()
    s.beginGroup(TOOLBAR_SETTINGS_GROUP)
    raw = s.value("order")
    s.endGroup()
    if not raw:
        return list(default_order)
    if isinstance(raw, (list, tuple)):
        ids = [str(item) for item in raw]
    else:
        ids = [part.strip() for part in str(raw).split(",") if part.strip()]
    known = set(labels)
    visible = normalize_toolbar_order(
        [tool_id for tool_id in ids if tool_id in known or is_toolbar_separator(tool_id)],
        known=known,
    )
    while visible and is_toolbar_separator(visible[0]):
        visible.pop(0)
    return normalize_toolbar_order(visible, known=known)


def save_toolbar_order(
    order: list[str],
    registry: PluginRegistry | None = None,
) -> None:
    from PySide6.QtCore import QSettings

    known = set(toolbar_labels_for(registry))
    cleaned = normalize_toolbar_order(order, known=known)
    s = QSettings()
    s.beginGroup(TOOLBAR_SETTINGS_GROUP)
    s.setValue("order", cleaned)
    s.endGroup()


def load_toolbar_font_size() -> int:
    from PySide6.QtCore import QSettings

    s = QSettings()
    s.beginGroup(TOOLBAR_SETTINGS_GROUP)
    size = int(s.value("font_size", DEFAULT_TOOLBAR_FONT_SIZE))
    s.endGroup()
    return max(MIN_TOOLBAR_FONT_SIZE, min(MAX_TOOLBAR_FONT_SIZE, size))


def save_toolbar_font_size(size: int) -> None:
    from PySide6.QtCore import QSettings

    s = QSettings()
    s.beginGroup(TOOLBAR_SETTINGS_GROUP)
    s.setValue("font_size", max(MIN_TOOLBAR_FONT_SIZE, min(MAX_TOOLBAR_FONT_SIZE, int(size))))
    s.endGroup()


class _ReorderOnlyListWidget(QListWidget):
    """List that only reorders rows — items cannot be freely stacked on each other."""

    def __init__(self, parent: QWidget | None = None) -> None:
        super().__init__(parent)
        self._drag_rows: list[int] = []
        self._drop_line_y: int | None = None
        self.setViewMode(QListView.ViewMode.ListMode)
        self.setMovement(QListView.Movement.Static)
        self.setFlow(QListView.Flow.TopToBottom)
        self.setWrapping(False)
        self.setSelectionMode(QAbstractItemView.SelectionMode.ExtendedSelection)
        self.setDragEnabled(True)
        self.setAcceptDrops(True)
        self.setDropIndicatorShown(False)  # custom insertion line in paintEvent
        self.setDragDropOverwriteMode(False)
        self.setDefaultDropAction(Qt.DropAction.CopyAction)
        # DragDrop (not InternalMove): we reorder ourselves in dropEvent.
        # InternalMove makes Qt delete the dragged rows after drop, which made
        # items disappear when combined with a custom dropEvent.
        self.setDragDropMode(QAbstractItemView.DragDropMode.DragDrop)

    def startDrag(self, supported_actions) -> None:  # noqa: N802
        selected = self.selectedItems()
        if not selected:
            return
        self._drag_rows = sorted({self.row(item) for item in selected if self.row(item) >= 0})
        indexes = self.selectedIndexes()
        if not indexes:
            return
        mime = self.model().mimeData(indexes)
        if mime is None:
            return
        drag = QDrag(self)
        drag.setMimeData(mime)
        first = selected[0]
        rect = self.visualItemRect(first)
        if rect.isValid():
            drag.setPixmap(self.viewport().grab(rect))
            drag.setHotSpot(rect.center() - rect.topLeft())
        # Prefer Copy so Qt never removes source rows after the drag finishes.
        drag.exec(
            Qt.DropAction.CopyAction | Qt.DropAction.MoveAction,
            Qt.DropAction.CopyAction,
        )
        self._drag_rows = []
        self._clear_drop_indicator()

    def _insertion_y_for_pos(self, pos) -> int:
        """Y coordinate of the insertion line for the pointer position."""
        target_item = self.itemAt(pos)
        if target_item is None:
            if self.count() == 0:
                return 2
            last = self.item(self.count() - 1)
            return self.visualItemRect(last).bottom()
        rect = self.visualItemRect(target_item)
        if pos.y() > rect.center().y():
            return rect.bottom()
        return rect.top()

    def _set_drop_indicator(self, y: int | None) -> None:
        if self._drop_line_y == y:
            return
        self._drop_line_y = y
        self.viewport().update()

    def _clear_drop_indicator(self) -> None:
        self._set_drop_indicator(None)

    def dropEvent(self, event: QDropEvent) -> None:  # noqa: N802
        if event.source() is not self:
            event.ignore()
            self._clear_drop_indicator()
            return

        rows = list(self._drag_rows)
        if not rows:
            rows = sorted({self.row(item) for item in self.selectedItems() if self.row(item) >= 0})
        if not rows:
            event.ignore()
            self._clear_drop_indicator()
            return

        entries: list[tuple[str, object]] = []
        for index in range(self.count()):
            item = self.item(index)
            if item is None:
                continue
            entries.append((item.text(), item.data(Qt.ItemDataRole.UserRole)))

        moving = [entries[row] for row in rows if 0 <= row < len(entries)]
        if not moving:
            event.ignore()
            self._clear_drop_indicator()
            return
        remaining = [entry for index, entry in enumerate(entries) if index not in set(rows)]

        pos = event.position().toPoint()
        target_item = self.itemAt(pos)
        if target_item is None:
            drop_row = len(remaining)
        else:
            target_row = self.row(target_item)
            # Map target row into the "remaining" list (rows already removed).
            drop_row = sum(1 for row in range(target_row) if row not in set(rows))
            rect = self.visualItemRect(target_item)
            if pos.y() > rect.center().y():
                drop_row += 1
            drop_row = max(0, min(drop_row, len(remaining)))

        new_order = remaining[:drop_row] + moving + remaining[drop_row:]
        self.clear()
        selected_ids: set[object] = {data for _text, data in moving}
        first_selected: QListWidgetItem | None = None
        for text, data in new_order:
            item = QListWidgetItem(text)
            item.setData(Qt.ItemDataRole.UserRole, data)
            item.setFlags(
                Qt.ItemFlag.ItemIsEnabled
                | Qt.ItemFlag.ItemIsSelectable
                | Qt.ItemFlag.ItemIsDragEnabled
                | Qt.ItemFlag.ItemIsDropEnabled
            )
            self.addItem(item)
            if data in selected_ids:
                item.setSelected(True)
                if first_selected is None:
                    first_selected = item
        if first_selected is not None:
            self.setCurrentItem(first_selected)

        self._clear_drop_indicator()
        event.setDropAction(Qt.DropAction.CopyAction)
        event.accept()

    def dragMoveEvent(self, event) -> None:  # noqa: N802
        if event.source() is self:
            self._set_drop_indicator(self._insertion_y_for_pos(event.position().toPoint()))
            event.setDropAction(Qt.DropAction.CopyAction)
            event.accept()
            return
        self._clear_drop_indicator()
        event.ignore()

    def dragEnterEvent(self, event) -> None:  # noqa: N802
        if event.source() is self:
            self._set_drop_indicator(self._insertion_y_for_pos(event.position().toPoint()))
            event.setDropAction(Qt.DropAction.CopyAction)
            event.accept()
            return
        event.ignore()

    def dragLeaveEvent(self, event) -> None:  # noqa: N802
        self._clear_drop_indicator()
        super().dragLeaveEvent(event)

    def paintEvent(self, event) -> None:  # noqa: N802
        super().paintEvent(event)
        if self._drop_line_y is None:
            return
        painter = QPainter(self.viewport())
        painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        y = self._drop_line_y
        left = 6
        right = self.viewport().width() - 6
        pen = QPen(QColor("#0a84ff"), 2)
        pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        painter.setPen(pen)
        painter.drawLine(left + 4, y, right - 4, y)
        # Small end caps so the insertion point is easy to see.
        painter.setBrush(QColor("#0a84ff"))
        painter.setPen(Qt.PenStyle.NoPen)
        painter.drawEllipse(left - 1, y - 3, 6, 6)
        painter.drawEllipse(right - 5, y - 3, 6, 6)
        painter.end()

class ToolbarCustomizeWidget(QWidget):
    """Inline panel to choose toolbar buttons and font size."""

    def __init__(
        self,
        catalog: list[tuple[str, bool]],
        parent: QWidget | None = None,
        *,
        registry: PluginRegistry | None = None,
    ) -> None:
        super().__init__(parent)
        self._catalog = list(catalog)
        self._registry = registry
        self._labels = toolbar_labels_for(registry)
        self._default_order = default_toolbar_order_for(registry)

        hint = QLabel("Drag to change button order.")
        hint.setWordWrap(True)

        self._font_size = QSpinBox()
        self._font_size.setRange(MIN_TOOLBAR_FONT_SIZE, MAX_TOOLBAR_FONT_SIZE)
        self._font_size.setValue(load_toolbar_font_size())
        self._font_size.setSuffix(" pt")
        self._font_size.setFixedWidth(88)
        self._font_size.setToolTip("Font size for function toolbar button labels.")

        font_row = QHBoxLayout()
        font_row.setContentsMargins(0, 0, 0, 0)
        font_row.setSpacing(8)
        font_label = QLabel("Toolbar font size:")
        font_row.addWidget(font_label)
        font_row.addWidget(self._font_size)
        font_row.addStretch()

        self._list = _ReorderOnlyListWidget()
        self._list.setMinimumHeight(280)
        self._populate(self._catalog)

        button_row = QHBoxLayout()
        button_row.setContentsMargins(0, 4, 0, 0)
        button_row.setSpacing(8)

        self._add_btn = QPushButton("Add button…")
        self._add_btn.setToolTip(
            "Add File, Sequence, Tools, or Plugins items, or a separator"
        )
        self._add_menu = QMenu(self)
        self._add_btn.clicked.connect(self._show_add_menu)
        button_row.addWidget(self._add_btn)

        remove_btn = QPushButton("Remove")
        remove_btn.setToolTip("Remove the selected toolbar button(s)")
        remove_btn.clicked.connect(self._remove_selected)
        button_row.addWidget(remove_btn)

        reset_btn = QPushButton("Reset to default")
        reset_btn.clicked.connect(self._reset_default)
        button_row.addStretch()
        button_row.addWidget(reset_btn)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 12)
        layout.setSpacing(10)
        layout.addWidget(hint)
        layout.addLayout(font_row)
        layout.addWidget(self._list, stretch=1)
        layout.addSpacing(6)
        layout.addLayout(button_row)

    def _current_tool_ids(self) -> set[str]:
        return set(self._list_tool_ids())

    def _list_tool_ids(self) -> list[str]:
        ids: list[str] = []
        for row in range(self._list.count()):
            item = self._list.item(row)
            if item is None:
                continue
            tool_id = item.data(Qt.ItemDataRole.UserRole)
            if tool_id:
                ids.append(str(tool_id))
        return ids

    def _append_list_item(self, tool_id: str, *, visible: bool = True) -> QListWidgetItem:
        del visible
        item = QListWidgetItem(toolbar_item_label(tool_id, self._labels))
        item.setData(Qt.ItemDataRole.UserRole, tool_id)
        item.setFlags(
            Qt.ItemFlag.ItemIsEnabled
            | Qt.ItemFlag.ItemIsSelectable
            | Qt.ItemFlag.ItemIsDragEnabled
        )
        self._list.addItem(item)
        return item

    def _show_add_menu(self) -> None:
        self._populate_add_menu()
        self._add_menu.exec(self._add_btn.mapToGlobal(self._add_btn.rect().bottomLeft()))

    def _add_category_menu(
        self,
        title: str,
        tool_ids: tuple[str, ...],
        current_ids: set[str],
    ) -> None:
        menu = self._add_menu.addMenu(title)
        added = 0
        for tool_id in tool_ids:
            if tool_id in current_ids or tool_id in TOOLBAR_CUSTOMIZE_HIDDEN:
                continue
            menu.addAction(
                toolbar_item_label(tool_id, self._labels),
                lambda checked=False, tid=tool_id: self._add_tool(tid),
            )
            added += 1
        if added == 0:
            empty = menu.addAction(f"(all {title} tools listed)")
            empty.setEnabled(False)

    def _populate_add_menu(self) -> None:
        self._add_menu.clear()
        current_ids = self._current_tool_ids()

        self._add_category_menu("File", FILE_TOOLBAR_IDS, current_ids)
        self._add_category_menu("Sequence", SEQUENCE_TOOLBAR_IDS, current_ids)
        self._add_category_menu("Tools", TOOLS_TOOLBAR_IDS, current_ids)

        plugins_menu = self._add_menu.addMenu("Plugins")
        plugins_added = 0
        if self._registry is not None:
            for plugin in self._registry.all_plugins():
                tool_id = plugin_tool_id(plugin.id)
                if tool_id in current_ids:
                    continue
                plugins_menu.addAction(
                    plugin.label,
                    lambda checked=False, tid=tool_id: self._add_tool(tid),
                )
                plugins_added += 1
        if plugins_added == 0:
            empty = plugins_menu.addAction("(all plugins listed)")
            empty.setEnabled(False)

        self._add_menu.addSeparator()
        self._add_menu.addAction(
            "— separator —",
            lambda: self._add_tool(next_separator_tool_id(current_ids)),
        )

    def _add_tool(self, tool_id: str) -> None:
        if tool_id in TOOLBAR_CUSTOMIZE_HIDDEN:
            return
        if is_toolbar_separator(tool_id):
            tool_id = next_separator_tool_id(self._current_tool_ids())
        if tool_id in self._current_tool_ids():
            return
        item = self._append_list_item(tool_id, visible=True)
        self._list.setCurrentItem(item)
        self._list.scrollToItem(item)

    def _remove_selected(self) -> None:
        selected = self._list.selectedItems()
        if not selected:
            return
        rows = sorted(
            (self._list.row(item) for item in selected),
            reverse=True,
        )
        for row in rows:
            if row >= 0:
                self._list.takeItem(row)

    def _populate(self, catalog: list[tuple[str, bool]]) -> None:
        self._list.clear()
        ordered_ids = normalize_toolbar_order(
            [tool_id for tool_id, _visible in catalog if tool_id not in TOOLBAR_CUSTOMIZE_HIDDEN],
            known=set(self._labels),
        )
        for tool_id in ordered_ids:
            if tool_id not in self._labels and not is_toolbar_separator(tool_id):
                continue
            self._append_list_item(tool_id, visible=True)

    def _reset_default(self) -> None:
        self._font_size.setValue(DEFAULT_TOOLBAR_FONT_SIZE)
        self._populate([(tool_id, True) for tool_id in self._default_order])

    def has_visible_tools(self) -> bool:
        return bool(self.result_order())

    def result_order(self) -> list[str]:
        return normalize_toolbar_order(self._list_tool_ids(), known=set(self._labels))

    def result_font_size(self) -> int:
        return self._font_size.value()


class ToolbarCustomizeDialog(QDialog):
    """Standalone dialog wrapping :class:`ToolbarCustomizeWidget`."""

    def __init__(
        self,
        catalog: list[tuple[str, bool]],
        parent: QWidget | None = None,
        *,
        registry: PluginRegistry | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Customize function toolbar")
        self.setMinimumWidth(420)
        self._panel = ToolbarCustomizeWidget(catalog, self, registry=registry)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.addWidget(self._panel)
        layout.addWidget(buttons)

    def result_order(self) -> list[str]:
        return self._panel.result_order()

    def result_font_size(self) -> int:
        return self._panel.result_font_size()

    def accept(self) -> None:
        if not self._panel.has_visible_tools():
            QMessageBox.warning(
                self,
                "Customize function toolbar",
                "Show at least one toolbar button.",
            )
            return
        super().accept()
