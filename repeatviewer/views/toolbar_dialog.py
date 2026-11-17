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
    QMessageBox,
    QPushButton,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)

from repeatviewer.widgets.spin_boxes import NoWheelSpinBox as QSpinBox


if TYPE_CHECKING:
    from repeatviewer.plugins.registry import PluginRegistry

TOOLBAR_SETTINGS_GROUP = "toolbar"
# One-shot flag: the "drop Repeat Quantify from the old default" migration in
# load_toolbar_order() must run at most once, or it silently deletes the button
# again every launch for anyone who deliberately adds it back.
REPEAT_QUANTIFY_MIGRATION_KEY = "repeat_quantify_migrated_v1"
# One-shot: users still on the previous curated default (MSA Clean / Extend /
# Consensus) move to Mapping / Tree / Primers without wiping custom toolbars.
DEFAULT_TOOLBAR_V2_MIGRATION_KEY = "default_toolbar_v2_migrated"
DEFAULT_TOOLBAR_FONT_SIZE = 12
MAX_TOOLBAR_FONT_SIZE = 20
MIN_TOOLBAR_FONT_SIZE = 8
# List height for the default set (~8 rows); scrolls when more buttons are added.
_TOOLBAR_LIST_ROW_PAD = 8
_TOOLBAR_LIST_ROW_FALLBACK = 26

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
    "ncbi_download",
    "aliview",
    "blast",
    "hmmer",
    "align",
    "fastga",
    "tree",
    "map_reference",
    "primers",
    "repeat_quantify",
    "extend",
    "teaid",
    "msa_clean",
    "consensus",
    "flexidot",
    "te_annotation",
    "te_genome_map",
)

# Built-in tools known to the customize dialog.
ALL_BUILTIN_TOOLBAR_IDS: tuple[str, ...] = (
    *FILE_TOOLBAR_IDS,
    *SEQUENCE_TOOLBAR_IDS,
    *TOOLS_TOOLBAR_IDS,
    "_sep1",
)

# Default visible toolbar (Settings → Toolbar → Reset to default).
DEFAULT_TOOLBAR_ORDER: list[str] = [
    "blast",
    "align",
    "flexidot",
    "teaid",
    "map_reference",
    "tree",
    "primers",
    "te_genome_map",
]

# Previous factory default; used only for the one-shot v2 migration.
_PREVIOUS_DEFAULT_TOOLBAR_ORDER: tuple[str, ...] = (
    "blast",
    "align",
    "flexidot",
    "teaid",
    "msa_clean",
    "extend",
    "consensus",
    "te_genome_map",
)

TOOLBAR_LABELS: dict[str, str] = {
    "define_genome": "Manage Database",
    "new_sequence": "New Sequence",
    "group_sequences": "Combine Sequences to One File",
    "extract_sequences": "Extract Sequences from File",
    "reverse_complement": "Reverse Complement",
    "orf_prediction": "ORF Prediction",
    "translate": "Translate",
    "aliview": "AliView",
    "ncbi_download": "NCBI Download",
    "_sep1": "— separator —",
    "align": "Align",
    "fastga": "FastGA",
    "tree": "Tree",
    "map_reference": "Mapping",
    "primers": "Primers and Probes",
    "repeat_quantify": "Repeat Quantify",
    "blast": "BLAST",
    "hmmer": "HMMER",
    "extend": "MSA Extend",
    "teaid": "TE-Aid",
    "msa_clean": "MSA Clean",
    "consensus": "Consensus",
    "flexidot": "Dotplot",
    "te_annotation": "Sequence Annotation",
    "te_genome_map": "Visualization",
    "te_divergence_plot": "TE divergence",
    "workflows_manage": "Design and manage workflows",
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


def workflow_tool_id(workflow_id: str) -> str:
    return f"workflow:{workflow_id}"


def toolbar_labels_for(registry: PluginRegistry | None = None) -> dict[str, str]:
    labels = dict(TOOLBAR_LABELS)
    if registry is not None:
        for plugin in registry.all_plugins():
            labels[plugin_tool_id(plugin.id)] = plugin.label
    try:
        from repeatviewer.workflows import load_workflow_specs

        for spec in load_workflow_specs():
            labels[workflow_tool_id(spec.id)] = spec.name
    except Exception:  # noqa: BLE001
        pass
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


def _settings_flag_done(key: str) -> bool:
    """True once a one-shot toolbar migration has run (and mark it done)."""
    from PySide6.QtCore import QSettings

    s = QSettings()
    s.beginGroup(TOOLBAR_SETTINGS_GROUP)
    raw = s.value(key)
    if isinstance(raw, str):
        done = raw.strip().lower() in {"1", "true", "yes", "on"}
    else:
        done = bool(raw)
    if not done:
        s.setValue(key, True)
    s.endGroup()
    return done


def _repeat_quantify_migration_done() -> bool:
    """True once the one-shot Repeat Quantify migration has run (and mark it done)."""
    return _settings_flag_done(REPEAT_QUANTIFY_MIGRATION_KEY)


def _default_toolbar_v2_migration_done() -> bool:
    """True once the curated-default v2 migration has run (and mark it done)."""
    return _settings_flag_done(DEFAULT_TOOLBAR_V2_MIGRATION_KEY)


# Old ids that became one button: Primers + Map Primers -> Primers and Probes.
# A rename only (the old ids no longer exist), so it is safe on every load.
_RENAMED_TOOLBAR_IDS = {"primer_design": "primers", "primer_map": "primers"}


def _rename_old_toolbar_ids(ids: list[str]) -> list[str]:
    out: list[str] = []
    for tool_id in ids:
        new = _RENAMED_TOOLBAR_IDS.get(tool_id, tool_id)
        if new in out and not is_toolbar_separator(new):
            continue
        out.append(new)
    return out


def load_toolbar_order(registry: PluginRegistry | None = None) -> list[str]:
    from PySide6.QtCore import QSettings

    labels = toolbar_labels_for(registry)
    default_order = default_toolbar_order_for(registry)
    migration_done = _repeat_quantify_migration_done()
    v2_done = _default_toolbar_v2_migration_done()
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
    ids = _rename_old_toolbar_ids(ids)
    known = set(labels)
    visible = normalize_toolbar_order(
        [tool_id for tool_id in ids if tool_id in known or is_toolbar_separator(tool_id)],
        known=known,
    )
    while visible and is_toolbar_separator(visible[0]):
        visible.pop(0)
    visible = normalize_toolbar_order(visible, known=known)
    # Migrate users still on the previous default (which included Repeat Quantify).
    # One-shot only: after this has run once, Repeat Quantify is a normal toolbar
    # button the user is free to keep.
    if not migration_done and "repeat_quantify" in visible:
        stripped = [tool_id for tool_id in visible if tool_id != "repeat_quantify"]
        if stripped == list(default_order):
            save_toolbar_order(stripped, registry)
            return stripped
    # Users still on the previous curated default → new Mapping/Tree/Primers set.
    if not v2_done and visible == list(_PREVIOUS_DEFAULT_TOOLBAR_ORDER):
        save_toolbar_order(default_order, registry)
        return list(default_order)
    return visible


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

class _AddToolbarButtonsDialog(QDialog):
    """Pick one or more tools to add to the toolbar in one step."""

    def __init__(
        self,
        *,
        available: list[tuple[str, list[tuple[str, str]]]],
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Add toolbar buttons")
        self.setMinimumWidth(420)
        self.setMinimumHeight(360)

        hint = QLabel("Check one or more functions, then click Add.")
        hint.setWordWrap(True)

        self._tree = QTreeWidget(self)
        self._tree.setHeaderHidden(True)
        self._tree.setRootIsDecorated(True)
        self._tree.setUniformRowHeights(True)
        self._tree.setSelectionMode(QAbstractItemView.SelectionMode.NoSelection)

        for category, entries in available:
            if not entries:
                continue
            parent_item = QTreeWidgetItem([category])
            parent_item.setFlags(
                Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsUserCheckable
            )
            parent_item.setCheckState(0, Qt.CheckState.Unchecked)
            for tool_id, label in entries:
                child = QTreeWidgetItem([label])
                child.setFlags(
                    Qt.ItemFlag.ItemIsEnabled | Qt.ItemFlag.ItemIsUserCheckable
                )
                child.setCheckState(0, Qt.CheckState.Unchecked)
                child.setData(0, Qt.ItemDataRole.UserRole, tool_id)
                parent_item.addChild(child)
            self._tree.addTopLevelItem(parent_item)
            parent_item.setExpanded(True)

        self._tree.itemChanged.connect(self._on_item_changed)

        select_row = QHBoxLayout()
        select_row.setContentsMargins(0, 0, 0, 0)
        select_all = QPushButton("Select all")
        select_all.clicked.connect(lambda: self._set_all_checked(True))
        select_none = QPushButton("Select none")
        select_none.clicked.connect(lambda: self._set_all_checked(False))
        select_row.addWidget(select_all)
        select_row.addWidget(select_none)
        select_row.addStretch()

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        ok = buttons.button(QDialogButtonBox.StandardButton.Ok)
        if ok is not None:
            ok.setText("Add")
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(14, 14, 14, 14)
        layout.setSpacing(10)
        layout.addWidget(hint)
        layout.addWidget(self._tree, stretch=1)
        layout.addLayout(select_row)
        layout.addWidget(buttons)

    def _on_item_changed(self, item: QTreeWidgetItem, column: int) -> None:
        del column
        if item.childCount() <= 0:
            return
        state = item.checkState(0)
        self._tree.blockSignals(True)
        try:
            for index in range(item.childCount()):
                child = item.child(index)
                if child is not None:
                    child.setCheckState(0, state)
        finally:
            self._tree.blockSignals(False)

    def _set_all_checked(self, checked: bool) -> None:
        state = Qt.CheckState.Checked if checked else Qt.CheckState.Unchecked
        self._tree.blockSignals(True)
        try:
            for i in range(self._tree.topLevelItemCount()):
                parent = self._tree.topLevelItem(i)
                if parent is None:
                    continue
                parent.setCheckState(0, state)
                for j in range(parent.childCount()):
                    child = parent.child(j)
                    if child is not None:
                        child.setCheckState(0, state)
        finally:
            self._tree.blockSignals(False)

    def selected_tool_ids(self) -> list[str]:
        picked: list[str] = []
        for i in range(self._tree.topLevelItemCount()):
            parent = self._tree.topLevelItem(i)
            if parent is None:
                continue
            for j in range(parent.childCount()):
                child = parent.child(j)
                if child is None or child.checkState(0) != Qt.CheckState.Checked:
                    continue
                tool_id = child.data(0, Qt.ItemDataRole.UserRole)
                if tool_id:
                    picked.append(str(tool_id))
        return picked


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
        self._populate(self._catalog)

        button_row = QHBoxLayout()
        button_row.setContentsMargins(0, 4, 0, 0)
        button_row.setSpacing(8)

        self._add_btn = QPushButton("Add buttons…")
        self._add_btn.setToolTip(
            "Add one or more File, Sequence, Tools, Plugins, or Workflows items, "
            "or a separator"
        )
        self._add_btn.clicked.connect(self._show_add_dialog)
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
        # Fixed to the default set height so Preferences is not padded with
        # empty list space; extra buttons scroll inside the list.
        layout.addWidget(self._list, stretch=0)
        layout.addSpacing(6)
        layout.addLayout(button_row)
        layout.addStretch(1)

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

    def _available_add_categories(
        self,
    ) -> list[tuple[str, list[tuple[str, str]]]]:
        current_ids = self._current_tool_ids()
        categories: list[tuple[str, list[tuple[str, str]]]] = []

        def _entries(tool_ids: tuple[str, ...]) -> list[tuple[str, str]]:
            out: list[tuple[str, str]] = []
            for tool_id in tool_ids:
                if tool_id in current_ids or tool_id in TOOLBAR_CUSTOMIZE_HIDDEN:
                    continue
                out.append((tool_id, toolbar_item_label(tool_id, self._labels)))
            return out

        for title, ids in (
            ("File", FILE_TOOLBAR_IDS),
            ("Sequence", SEQUENCE_TOOLBAR_IDS),
            ("Tools", TOOLS_TOOLBAR_IDS),
        ):
            entries = _entries(ids)
            if entries:
                categories.append((title, entries))

        plugins: list[tuple[str, str]] = []
        wrappers: list[tuple[str, str]] = []
        if self._registry is not None:
            from repeatviewer.plugins.custom_command import UserCommandPlugin

            for plugin in self._registry.all_plugins():
                tool_id = plugin_tool_id(plugin.id)
                if tool_id in current_ids:
                    continue
                entry = (tool_id, plugin.label)
                if isinstance(plugin, UserCommandPlugin):
                    wrappers.append(entry)
                else:
                    plugins.append(entry)
        if plugins:
            categories.append(("Plugins", plugins))
        if wrappers:
            categories.append(("Command wrappers", wrappers))

        workflows: list[tuple[str, str]] = []
        if "workflows_manage" not in current_ids:
            workflows.append(
                (
                    "workflows_manage",
                    toolbar_item_label("workflows_manage", self._labels),
                )
            )
        try:
            from repeatviewer.workflows import load_workflow_specs

            for spec in load_workflow_specs():
                tool_id = workflow_tool_id(spec.id)
                if tool_id in current_ids:
                    continue
                workflows.append((tool_id, spec.name))
        except Exception:  # noqa: BLE001
            pass
        if workflows:
            categories.append(("Workflows", workflows))

        categories.append(
            ("Layout", [("__separator__", "— separator —")])
        )
        return categories

    def _show_add_dialog(self) -> None:
        available = self._available_add_categories()
        # Layout/separator alone still counts as something to add.
        has_tools = any(
            title != "Layout" and entries for title, entries in available
        )
        if not has_tools and not any(title == "Layout" for title, _ in available):
            QMessageBox.information(
                self,
                "Add toolbar buttons",
                "All available functions are already on the toolbar.",
            )
            return
        dialog = _AddToolbarButtonsDialog(available=available, parent=self)
        if dialog.exec() != QDialog.DialogCode.Accepted:
            return
        last_item: QListWidgetItem | None = None
        for tool_id in dialog.selected_tool_ids():
            last_item = self._add_tool(tool_id) or last_item
        if last_item is not None:
            self._list.setCurrentItem(last_item)
            self._list.scrollToItem(last_item)

    def _add_tool(self, tool_id: str) -> QListWidgetItem | None:
        if tool_id in TOOLBAR_CUSTOMIZE_HIDDEN:
            return None
        if tool_id == "__separator__" or is_toolbar_separator(tool_id):
            tool_id = next_separator_tool_id(self._current_tool_ids())
        if tool_id in self._current_tool_ids():
            return None
        return self._append_list_item(tool_id, visible=True)

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

    def _fit_list_height(self) -> None:
        """Size the list to the default toolbar set (scrolls when longer)."""
        rows = max(len(DEFAULT_TOOLBAR_ORDER), 1)
        hinted = self._list.sizeHintForRow(0) if self._list.count() else 0
        row_h = max(
            hinted,
            self._list.fontMetrics().height() + 10,
            _TOOLBAR_LIST_ROW_FALLBACK,
        )
        height = rows * row_h + _TOOLBAR_LIST_ROW_PAD
        self._list.setMinimumHeight(height)
        self._list.setMaximumHeight(height)

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
        self._fit_list_height()

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
