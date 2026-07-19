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
TOOLBAR_CUSTOMIZE_HIDDEN = frozenset({"te_divergence_plot"})

# Sequence menu tools that can be added to the toolbar (menu order).
SEQUENCE_TOOLBAR_IDS: tuple[str, ...] = (
    "aliview",
    "blast",
    "align",
    "tree",
    "extend",
    "tsd_search",
    "teaid",
    "msa_clean",
    "consensus",
    "flexidot",
    "te_annotation",
    "te_genome_map",
    "te_divergence_plot",
    "te_benchmark",
)


# Built-in tools known to the customize dialog.
ALL_BUILTIN_TOOLBAR_IDS: tuple[str, ...] = (
    "define_genome",
    *SEQUENCE_TOOLBAR_IDS,
    "_sep1",
    "_sep_plugins",
    "plugin:minimap2",
)

# Default visible toolbar (Settings → Toolbar → Reset to default).
DEFAULT_TOOLBAR_ORDER: list[str] = [
    "define_genome",
    "blast",
    "align",
    "tree",
    "msa_clean",
    "flexidot",
    "teaid",
    "te_genome_map",
    "aliview",
    "extend",
    "consensus",
    "te_benchmark",
    "_sep_plugins",
    "plugin:minimap2",
    "te_divergence_plot",
]

TOOLBAR_LABELS: dict[str, str] = {
    "define_genome": "Define genome",
    "aliview": "AliView",
    "_sep1": "— separator —",
    "align": "Align (MAFFT)",
    "tree": "Tree",
    "blast": "BLAST",
    "extend": "Extend",
    "teaid": "TEAid",
    "msa_clean": "MSA clean",
    "consensus": "Consensus",
    "flexidot": "Flexidot",
    "te_annotation": "TE annotation",
    "te_genome_map": "TE visualization",
    "te_benchmark": "Benchmarking",
    "te_divergence_plot": "TE divergence",
    "tsd_search": "TSD search",
    "_sep_plugins": "— plugins —",
    "plugin:minimap2": "Minimap2",
}


def is_toolbar_separator(tool_id: str) -> bool:
    return tool_id.startswith("_sep")


def toolbar_item_label(tool_id: str, labels: dict[str, str] | None = None) -> str:
    lookup = labels or TOOLBAR_LABELS
    if tool_id.startswith("_sep") and tool_id != "_sep_plugins":
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


def default_toolbar_order_for(registry: PluginRegistry | None = None) -> list[str]:
    order = list(DEFAULT_TOOLBAR_ORDER)
    if registry is None:
        return order
    plugin_tool_ids = [plugin_tool_id(p.id) for p in registry.all_plugins()]
    if not plugin_tool_ids:
        return [
            tool_id
            for tool_id in order
            if tool_id != "_sep_plugins" and not tool_id.startswith("plugin:")
        ]
    known = set(order)
    if "_sep_plugins" not in known:
        order.append("_sep_plugins")
        known.add("_sep_plugins")
    for tool_id in plugin_tool_ids:
        if tool_id not in known:
            order.append(tool_id)
            known.add(tool_id)
    return order


def build_toolbar_catalog(
    visible_order: list[str],
    registry: PluginRegistry | None = None,
) -> list[tuple[str, bool]]:
    """Toolbar customize list: only tools already on / listed for the toolbar."""
    labels = toolbar_labels_for(registry)
    catalog: list[tuple[str, bool]] = []
    seen: set[str] = set()
    for tool_id in visible_order:
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
    visible = [
        tool_id
        for tool_id in ids
        if tool_id in known or is_toolbar_separator(tool_id)
    ]
    while visible and visible[0].startswith("_sep"):
        visible.pop(0)
    if "define_genome" not in visible:
        if "blast" in visible:
            visible.insert(visible.index("blast"), "define_genome")
        else:
            visible.insert(0, "define_genome")
    for tool_id in default_order:
        if tool_id not in visible:
            visible.append(tool_id)
    return visible


def save_toolbar_order(
    order: list[str],
    registry: PluginRegistry | None = None,
) -> None:
    from PySide6.QtCore import QSettings

    known = set(toolbar_labels_for(registry))
    cleaned = [
        tool_id
        for tool_id in order
        if tool_id in known or is_toolbar_separator(tool_id)
    ]
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

        self._list = QListWidget()
        self._list.setMinimumHeight(280)
        self._list.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self._list.setDefaultDropAction(Qt.DropAction.MoveAction)
        self._populate(self._catalog)

        button_row = QHBoxLayout()
        button_row.setContentsMargins(0, 4, 0, 0)
        button_row.setSpacing(8)

        self._add_btn = QPushButton("Add button…")
        self._add_btn.setToolTip(
            "Add Sequence menu tools, plugins, or a separator to this list"
        )
        self._add_menu = QMenu(self)
        self._add_btn.clicked.connect(self._show_add_menu)
        button_row.addWidget(self._add_btn)

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

    def _append_list_item(self, tool_id: str, *, visible: bool = True) -> QListWidgetItem:
        item = QListWidgetItem(toolbar_item_label(tool_id, self._labels))
        item.setData(Qt.ItemDataRole.UserRole, tool_id)
        item.setFlags(
            item.flags()
            | Qt.ItemFlag.ItemIsUserCheckable
            | Qt.ItemFlag.ItemIsEnabled
            | Qt.ItemFlag.ItemIsSelectable
        )
        item.setCheckState(
            Qt.CheckState.Checked if visible else Qt.CheckState.Unchecked
        )
        self._list.addItem(item)
        return item

    def _show_add_menu(self) -> None:
        self._populate_add_menu()
        self._add_menu.exec(self._add_btn.mapToGlobal(self._add_btn.rect().bottomLeft()))

    def _populate_add_menu(self) -> None:
        self._add_menu.clear()
        current_ids = self._current_tool_ids()

        sequence_menu = self._add_menu.addMenu("Sequence")
        sequence_added = 0
        for tool_id in SEQUENCE_TOOLBAR_IDS:
            if tool_id in current_ids:
                continue
            label = toolbar_item_label(tool_id, self._labels)
            sequence_menu.addAction(
                label,
                lambda checked=False, tid=tool_id: self._add_tool(tid),
            )
            sequence_added += 1
        if sequence_added == 0:
            empty = sequence_menu.addAction("(all Sequence tools listed)")
            empty.setEnabled(False)

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

        if "_sep_plugins" not in current_ids:
            self._add_menu.addAction(
                toolbar_item_label("_sep_plugins", self._labels),
                lambda: self._add_tool("_sep_plugins"),
            )

        self._add_menu.addSeparator()
        self._add_menu.addAction(
            "— separator —",
            lambda: self._add_tool(next_separator_tool_id(current_ids)),
        )

    def _add_tool(self, tool_id: str) -> None:
        if is_toolbar_separator(tool_id) and tool_id != "_sep_plugins":
            tool_id = next_separator_tool_id(self._current_tool_ids())
        if tool_id in self._current_tool_ids():
            return
        item = self._append_list_item(tool_id, visible=True)
        self._list.setCurrentItem(item)
        self._list.scrollToItem(item)

    def _populate(self, catalog: list[tuple[str, bool]]) -> None:
        self._list.clear()
        for tool_id, visible in catalog:
            if tool_id not in self._labels and not is_toolbar_separator(tool_id):
                continue
            self._append_list_item(tool_id, visible=visible)

    def _reset_default(self) -> None:
        self._font_size.setValue(DEFAULT_TOOLBAR_FONT_SIZE)
        self._populate(
            [
                (tool_id, True)
                for tool_id in self._default_order
            ]
        )

    def has_visible_tools(self) -> bool:
        return bool(self.result_order())

    def result_order(self) -> list[str]:
        order: list[str] = []
        for row in range(self._list.count()):
            item = self._list.item(row)
            if item is None or item.checkState() != Qt.CheckState.Checked:
                continue
            tool_id = item.data(Qt.ItemDataRole.UserRole)
            if tool_id and (str(tool_id) in self._labels or is_toolbar_separator(str(tool_id))):
                order.append(str(tool_id))
        return order

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
