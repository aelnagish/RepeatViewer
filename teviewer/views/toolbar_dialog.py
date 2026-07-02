"""Toolbar appearance, visibility, and button order."""

from __future__ import annotations

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
    QVBoxLayout,
    QWidget,
)

TOOLBAR_SETTINGS_GROUP = "toolbar"
DEFAULT_TOOLBAR_FONT_SIZE = 13

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
}


def build_toolbar_catalog(visible_order: list[str]) -> list[tuple[str, bool]]:
    """Full toolbar catalog for the customize dialog: (tool_id, visible)."""
    known = set(TOOLBAR_LABELS)
    visible_set = set(visible_order)
    catalog: list[tuple[str, bool]] = [
        (tool_id, True) for tool_id in visible_order if tool_id in known
    ]
    for tool_id in DEFAULT_TOOLBAR_ORDER:
        if tool_id not in visible_set:
            catalog.append((tool_id, False))
    for tool_id in TOOLBAR_LABELS:
        if tool_id not in {entry[0] for entry in catalog}:
            catalog.append((tool_id, False))
    return catalog


def load_toolbar_order() -> list[str]:
    from PySide6.QtCore import QSettings

    s = QSettings()
    s.beginGroup(TOOLBAR_SETTINGS_GROUP)
    raw = s.value("order")
    s.endGroup()
    if not raw:
        return list(DEFAULT_TOOLBAR_ORDER)
    if isinstance(raw, (list, tuple)):
        ids = [str(item) for item in raw]
    else:
        ids = [part.strip() for part in str(raw).split(",") if part.strip()]
    known = set(TOOLBAR_LABELS)
    visible = [tool_id for tool_id in ids if tool_id in known]
    for tool_id in DEFAULT_TOOLBAR_ORDER:
        if tool_id not in visible:
            visible.append(tool_id)
    return visible


def save_toolbar_order(order: list[str]) -> None:
    from PySide6.QtCore import QSettings

    known = set(TOOLBAR_LABELS)
    cleaned = [tool_id for tool_id in order if tool_id in known]
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


class ToolbarCustomizeDialog(QDialog):
    """Choose which function-toolbar buttons are shown and drag to reorder them."""

    def __init__(
        self,
        catalog: list[tuple[str, bool]],
        parent: QWidget | None = None,
    ) -> None:
        super().__init__(parent)
        self.setWindowTitle("Customize function toolbar")
        self.setMinimumWidth(420)
        self._catalog = list(catalog)

        hint = QLabel(
            "Check items to show them on the toolbar. Drag rows to change button order. "
            "You can also drag the toolbar grip to move or float it."
        )
        hint.setWordWrap(True)

        self._list = QListWidget()
        self._list.setDragDropMode(QAbstractItemView.DragDropMode.InternalMove)
        self._list.setDefaultDropAction(Qt.DropAction.MoveAction)
        self._populate(self._catalog)

        reset_btn = QPushButton("Reset to default")
        reset_btn.clicked.connect(self._reset_default)

        buttons = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Ok | QDialogButtonBox.StandardButton.Cancel
        )
        buttons.accepted.connect(self.accept)
        buttons.rejected.connect(self.reject)

        row = QHBoxLayout()
        row.addStretch()
        row.addWidget(reset_btn)

        layout = QVBoxLayout(self)
        layout.addWidget(hint)
        layout.addWidget(self._list)
        layout.addLayout(row)
        layout.addWidget(buttons)

    def _populate(self, catalog: list[tuple[str, bool]]) -> None:
        self._list.clear()
        for tool_id, visible in catalog:
            if tool_id not in TOOLBAR_LABELS:
                continue
            item = QListWidgetItem(TOOLBAR_LABELS[tool_id])
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
        self._populate([(tool_id, True) for tool_id in DEFAULT_TOOLBAR_ORDER])

    def result_order(self) -> list[str]:
        order: list[str] = []
        for row in range(self._list.count()):
            item = self._list.item(row)
            if item is None or item.checkState() != Qt.CheckState.Checked:
                continue
            tool_id = item.data(Qt.ItemDataRole.UserRole)
            if tool_id:
                order.append(str(tool_id))
        return order

    def accept(self) -> None:
        if not self.result_order():
            QMessageBox.warning(
                self,
                "Customize function toolbar",
                "Show at least one toolbar button.",
            )
            return
        super().accept()
