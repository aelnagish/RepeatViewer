"""Dialog to install / uninstall pip-based RepeatViewer plugins."""

from __future__ import annotations

from pathlib import Path

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QListWidget,
    QListWidgetItem,
    QMessageBox,
    QPushButton,
    QVBoxLayout,
)

from repeatviewer.plugins.install import (
    InstalledEntryPlugin,
    install_plugin_path,
    list_installed_entry_plugins,
    uninstall_plugin_distribution,
)

class InstalledPluginManageDialog(QDialog):
    """List installed entry-point plugins; add from folder/wheel or delete."""

    def __init__(self, *, parent=None, on_changed=None) -> None:
        super().__init__(parent)
        self.setWindowTitle("Manage plugin")
        self.setMinimumWidth(560)
        self.setMinimumHeight(360)
        self._on_changed = on_changed
        self._items: list[InstalledEntryPlugin] = []

        hint = QLabel(
            "Install a plugin package (folder with pyproject.toml, or a .whl). "
            "Plugins register via the repeatviewer.plugins entry point.",
            self,
        )
        hint.setWordWrap(True)
        hint.setStyleSheet("color: #666666;")

        self._list = QListWidget(self)
        self._list.currentItemChanged.connect(self._update_buttons)

        self._add_btn = QPushButton("Add", self)
        self._delete_btn = QPushButton("Delete", self)
        self._add_btn.clicked.connect(self._add_plugin)
        self._delete_btn.clicked.connect(self._delete_plugin)

        btn_row = QHBoxLayout()
        btn_row.addWidget(self._add_btn)
        btn_row.addWidget(self._delete_btn)
        btn_row.addStretch(1)

        buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Close)
        buttons.rejected.connect(self.reject)
        close_btn = buttons.button(QDialogButtonBox.StandardButton.Close)
        if close_btn is not None:
            close_btn.clicked.connect(self.accept)

        layout = QVBoxLayout(self)
        layout.addWidget(hint)
        layout.addWidget(self._list, 1)
        layout.addLayout(btn_row)
        layout.addWidget(buttons)

        self._populate()
        self._update_buttons()

    def _populate(self) -> None:
        self._list.clear()
        self._items = list_installed_entry_plugins()
        if not self._items:
            item = QListWidgetItem("(No installed plugins — click Add)")
            item.setFlags(Qt.ItemFlag.NoItemFlags)
            self._list.addItem(item)
            return
        for index, info in enumerate(self._items, start=1):
            item = QListWidgetItem(f"{index}. {info.label}  ({info.distribution})")
            item.setData(Qt.ItemDataRole.UserRole, index - 1)
            item.setToolTip(
                f"Entry point: {info.entry_name}\n"
                f"Package: {info.distribution}\n"
                f"{info.module}"
            )
            self._list.addItem(item)

    def _selected(self) -> InstalledEntryPlugin | None:
        item = self._list.currentItem()
        if item is None:
            return None
        index = item.data(Qt.ItemDataRole.UserRole)
        if index is None:
            return None
        try:
            return self._items[int(index)]
        except (TypeError, ValueError, IndexError):
            return None

    def _update_buttons(self) -> None:
        self._delete_btn.setEnabled(self._selected() is not None)

    def _notify_changed(self) -> None:
        if self._on_changed is not None:
            self._on_changed()

    def _add_plugin(self) -> None:
        box = QMessageBox(self)
        box.setWindowTitle("Manage plugin")
        box.setText("Install a plugin from:")
        folder_btn = box.addButton("Folder", QMessageBox.ButtonRole.AcceptRole)
        wheel_btn = box.addButton("Wheel file", QMessageBox.ButtonRole.AcceptRole)
        box.addButton(QMessageBox.StandardButton.Cancel)
        box.exec()
        clicked = box.clickedButton()
        if clicked is folder_btn:
            folder = QFileDialog.getExistingDirectory(
                self,
                "Select plugin package folder",
                str(Path.home()),
            )
            if not folder:
                return
            target = Path(folder)
        elif clicked is wheel_btn:
            path_str, _ = QFileDialog.getOpenFileName(
                self,
                "Select plugin wheel or archive",
                str(Path.home()),
                "Python packages (*.whl *.tar.gz *.zip);;All files (*)",
            )
            if not path_str:
                return
            target = Path(path_str)
        else:
            return

        self.setEnabled(False)
        try:
            ok, message = install_plugin_path(target)
        finally:
            self.setEnabled(True)
        if not ok:
            QMessageBox.warning(self, "Manage plugin", message)
            return
        self._populate()
        self._update_buttons()
        self._notify_changed()
        detail = message.strip()
        info = "Plugin installed. It should appear under Tools → Plugin (Beta)."
        if detail:
            info = f"{info}\n\n{detail[-1500:]}"
        QMessageBox.information(self, "Manage plugin", info)

    def _delete_plugin(self) -> None:
        info = self._selected()
        if info is None:
            return
        answer = QMessageBox.question(
            self,
            "Manage plugin",
            f"Uninstall package “{info.distribution}” "
            f"(entry point “{info.entry_name}”)?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No,
        )
        if answer != QMessageBox.StandardButton.Yes:
            return
        self.setEnabled(False)
        try:
            ok, message = uninstall_plugin_distribution(info.distribution)
        finally:
            self.setEnabled(True)
        if not ok:
            QMessageBox.warning(self, "Manage plugin", message)
            return
        self._populate()
        self._update_buttons()
        self._notify_changed()
        QMessageBox.information(
            self,
            "Manage plugin",
            f"Uninstalled “{info.distribution}”.",
        )
