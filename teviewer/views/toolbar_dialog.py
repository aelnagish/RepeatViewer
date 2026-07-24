"""Toolbar appearance, visibility, and button order."""

from __future__ import annotations

from typing import TYPE_CHECKING

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QAbstractItemView,
    QDialog,
    QDialogButtonBox,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMenu,
    QMessageBox,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
    from teviewer.plugins.registry import PluginRegistry

TOOLBAR_SETTINGS_GROUP = "toolbar"
DEFAULT_TOOLBAR_FONT_SIZE = 12
MAX_TOOLBAR_FONT_SIZE = 20
MIN_TOOLBAR_FONT_SIZE = 8

# Always first on the toolbar / in Settings → Toolbar.
FIXED_TOOLBAR_PREFIX: tuple[str, ...] = ("blast", "align")

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
    "tree",
    "fastga",
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
    "define_genome",
    "tree",
    "fastga",
    "msa_clean",
    "flexidot",
    "teaid",
    "te_genome_map",
    "aliview",
    "extend",
    "consensus",
    "te_benchmark",
    "plugin:minimap2",
]

TOOLBAR_LABELS: dict[str, str] = {
    "define_genome": "Define genome",
    "new_sequence": "New Sequence…",
    "group_sequences": "Group Sequences into a List…",
    "extract_sequences": "Extract Sequences from List…",
    "reverse_complement": "Reverse Complement",
    "orf_prediction": "ORF Prediction",
    "translate": "Translate",
    "aliview": "AliView",
    "_sep1": "— separator —",
    "align": "Align",
    "tree": "Tree",
    "fastga": "FastGA",
    "blast": "BLAST",
    "extend": "MSA Extend",
    "teaid": "TEAid",
    "msa_clean": "MSA clean",
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
    """Drop removed ids, dedupe, and pin BLAST / Align at the front."""
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
    fixed = [tool_id for tool_id in FIXED_TOOLBAR_PREFIX if tool_id in cleaned]
    rest = [tool_id for tool_id in cleaned if tool_id not in FIXED_TOOLBAR_PREFIX]
    return fixed + rest


def default_toolbar_order_for(registry: PluginRegistry | None = None) -> list[str]:
    order = list(DEFAULT_TOOLBAR_ORDER)
    if registry is None:
        return normalize_toolbar_order(order)
    plugin_tool_ids = [plugin_tool_id(p.id) for p in registry.all_plugins()]
    known = set(order)
    for tool_id in plugin_tool_ids:
        if tool_id not in known:
            order.append(tool_id)
            known.add(tool_id)
    return normalize_toolbar_order(order)


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
    if "define_genome" not in visible:
        insert_at = len([t for t in FIXED_TOOLBAR_PREFIX if t in visible])
        visible.insert(insert_at, "define_genome")
    for tool_id in default_order:
        if tool_id not in visible and tool_id not in TOOLBAR_CUSTOMIZE_HIDDEN:
            visible.append(tool_id)
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
        self._enforcing_order = False

        hint = QLabel(
            "Drag to change button order. BLAST and Align stay at positions 1 and 2."
        )
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

        self._list = QListWidget()
        self._list.setMinimumHeight(280)
        self._list.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self._list.setDefaultDropAction(Qt.DropAction.MoveAction)
        self._list.model().rowsMoved.connect(self._enforce_fixed_prefix)
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
        remove_btn.setToolTip("Remove the selected toolbar button")
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
        ids: set[str] = set()
        for row in range(self._list.count()):
            item = self._list.item(row)
            if item is None:
                continue
            tool_id = item.data(Qt.ItemDataRole.UserRole)
            if tool_id:
                ids.add(str(tool_id))
        return ids

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
            | Qt.ItemFlag.ItemIsDropEnabled
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
            if tool_id in FIXED_TOOLBAR_PREFIX:
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
        self._enforce_fixed_prefix()
        self._list.setCurrentItem(item)
        self._list.scrollToItem(item)

    def _remove_selected(self) -> None:
        item = self._list.currentItem()
        if item is None:
            return
        tool_id = item.data(Qt.ItemDataRole.UserRole)
        if tool_id in FIXED_TOOLBAR_PREFIX:
            QMessageBox.information(
                self,
                "Toolbar",
                "BLAST and Align stay on the toolbar and cannot be removed.",
            )
            return
        row = self._list.row(item)
        self._list.takeItem(row)

    def _populate(self, catalog: list[tuple[str, bool]]) -> None:
        self._list.clear()
        ordered_ids = normalize_toolbar_order(
            [tool_id for tool_id, _visible in catalog if tool_id not in TOOLBAR_CUSTOMIZE_HIDDEN],
            known=set(self._labels),
        )
        for index, tool_id in enumerate(FIXED_TOOLBAR_PREFIX):
            if tool_id not in ordered_ids and tool_id in self._labels:
                ordered_ids.insert(index, tool_id)
        for tool_id in ordered_ids:
            if tool_id not in self._labels and not is_toolbar_separator(tool_id):
                continue
            self._append_list_item(tool_id, visible=True)

    def _enforce_fixed_prefix(self, *_args) -> None:
        if self._enforcing_order:
            return
        entries = self._list_tool_ids()
        if not entries:
            return
        ordered_ids = normalize_toolbar_order(entries, known=set(self._labels))
        for tool_id in FIXED_TOOLBAR_PREFIX:
            if tool_id not in ordered_ids and tool_id in self._labels:
                insert_at = list(FIXED_TOOLBAR_PREFIX).index(tool_id)
                ordered_ids.insert(insert_at, tool_id)
        if entries == ordered_ids:
            return
        self._enforcing_order = True
        try:
            self._populate([(tool_id, True) for tool_id in ordered_ids])
        finally:
            self._enforcing_order = False

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
