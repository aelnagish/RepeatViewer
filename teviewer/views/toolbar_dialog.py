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
    QMessageBox,
    QPushButton,
    QSpinBox,
    QVBoxLayout,
    QWidget,
)

if TYPE_CHECKING:
    from teviewer.plugins.registry import PluginRegistry

TOOLBAR_SETTINGS_GROUP = "toolbar"
DEFAULT_TOOLBAR_FONT_SIZE = 13
TOOLBAR_CUSTOMIZE_HIDDEN = frozenset({"te_divergence_plot"})

DEFAULT_TOOLBAR_ORDER: list[str] = [
    "aliview",
    "_sep1",
    "align",
    "blast",
    "extend",
    "teaid",
    "msa_clean",
    "consensus",
    "flexidot",
    "te_genome_map",
    "te_benchmark",
    "_sep_plugins",
    "plugin:minimap2",
]

TOOLBAR_LABELS: dict[str, str] = {
    "aliview": "AliView",
    "_sep1": "— separator —",
    "align": "Align (MAFFT)",
    "blast": "Blast",
    "extend": "Extend",
    "teaid": "TEAid",
    "msa_clean": "MSA clean",
    "consensus": "Consensus",
    "flexidot": "Flexidot",
    "te_genome_map": "TE visualization",
    "te_benchmark": "Benchmarking",
    "te_divergence_plot": "TE divergence",
    "_sep_plugins": "— plugins —",
    "plugin:minimap2": "Minimap2",
}


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
    """Full toolbar catalog for the customize dialog: (tool_id, visible)."""
    labels = toolbar_labels_for(registry)
    default_order = default_toolbar_order_for(registry)
    known = set(labels)
    visible_set = set(visible_order)
    catalog: list[tuple[str, bool]] = [
        (tool_id, True)
        for tool_id in visible_order
        if tool_id in known and tool_id not in TOOLBAR_CUSTOMIZE_HIDDEN
    ]
    for tool_id in default_order:
        if tool_id not in visible_set and tool_id not in TOOLBAR_CUSTOMIZE_HIDDEN:
            catalog.append((tool_id, False))
    for tool_id in labels:
        if tool_id in TOOLBAR_CUSTOMIZE_HIDDEN:
            continue
        if tool_id not in {entry[0] for entry in catalog}:
            catalog.append((tool_id, False))
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
        if tool_id in known and tool_id not in TOOLBAR_CUSTOMIZE_HIDDEN
    ]
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
        if tool_id in known and tool_id not in TOOLBAR_CUSTOMIZE_HIDDEN
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
    return max(8, min(16, size))


def save_toolbar_font_size(size: int) -> None:
    from PySide6.QtCore import QSettings

    s = QSettings()
    s.beginGroup(TOOLBAR_SETTINGS_GROUP)
    s.setValue("font_size", max(8, min(16, int(size))))
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
        self._font_size.setRange(8, 16)
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

        reset_btn = QPushButton("Reset to default")
        reset_btn.clicked.connect(self._reset_default)

        reset_row = QHBoxLayout()
        reset_row.setContentsMargins(0, 4, 0, 0)
        reset_row.addStretch()
        reset_row.addWidget(reset_btn)

        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 12)
        layout.setSpacing(10)
        layout.addWidget(hint)
        layout.addLayout(font_row)
        layout.addWidget(self._list, stretch=1)
        layout.addSpacing(6)
        layout.addLayout(reset_row)

    def _populate(self, catalog: list[tuple[str, bool]]) -> None:
        self._list.clear()
        for tool_id, visible in catalog:
            if tool_id not in self._labels or tool_id in TOOLBAR_CUSTOMIZE_HIDDEN:
                continue
            item = QListWidgetItem(self._labels[tool_id])
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

    def _reset_default(self) -> None:
        self._font_size.setValue(DEFAULT_TOOLBAR_FONT_SIZE)
        self._populate(
            [
                (tool_id, True)
                for tool_id in self._default_order
                if tool_id not in TOOLBAR_CUSTOMIZE_HIDDEN
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
            if tool_id and str(tool_id) not in TOOLBAR_CUSTOMIZE_HIDDEN:
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
