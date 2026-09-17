"""Folders panel: Cmd/Ctrl+click and Shift+click keep a multi-selection.

Opening a folder refreshes the Files panel, which calls back into the tree
(select_folder_item / set_active_folder). Those used to reset the selection
while Qt was still handling the click, so Cmd+click emptied it.
"""

from pathlib import Path

import pytest
from PySide6.QtCore import Qt
from PySide6.QtTest import QTest
from PySide6.QtWidgets import QTreeWidgetItem

mw = pytest.importorskip("repeatviewer.main_window")

M = Qt.KeyboardModifier
NAMES = ("a", "b", "c", "d")


@pytest.fixture
def tree(qapp, tmp_path):
    tree = mw.FolderTreeWidget()
    root = QTreeWidgetItem([tmp_path.name])
    root.setData(0, mw.ROLE_PATH, str(tmp_path))
    tree.addTopLevelItem(root)
    for name in NAMES:
        (tmp_path / name).mkdir()
        child = QTreeWidgetItem([name])
        child.setData(0, mw.ROLE_PATH, str(tmp_path / name))
        root.addChild(child)
    root.setExpanded(True)

    # What MainWindow does when a folder is opened.
    def on_folder_selected(folder: Path) -> None:
        item = tree.find_item_for_folder(folder)
        tree.select_folder_item(item)
        tree.set_active_folder(folder)

    tree.folder_selected.connect(on_folder_selected)
    tree.opened = []
    tree.folder_selected.connect(lambda f: tree.opened.append(f.name))
    tree.resize(300, 300)
    tree.show()
    qapp.processEvents()
    yield tree
    tree.close()


def _selected(tree) -> list[str]:
    return sorted(item.text(0) for item in tree.selectedItems())


def _click(tree, name, modifier=M.NoModifier):
    item = next(
        tree.topLevelItem(0).child(i)
        for i in range(tree.topLevelItem(0).childCount())
        if tree.topLevelItem(0).child(i).text(0) == name
    )
    pos = tree.visualItemRect(item).center()
    QTest.mouseClick(tree.viewport(), Qt.MouseButton.LeftButton, modifier, pos)
    QTest.qWait(10)  # let the deferred "unsuppress" timers run, as in the app
    return _selected(tree)


def test_cmd_click_adds_and_removes_folders(tree):
    assert _click(tree, "a") == ["a"]
    assert _click(tree, "b", M.ControlModifier) == ["a", "b"]
    assert _click(tree, "d", M.ControlModifier) == ["a", "b", "d"]
    assert _click(tree, "b", M.ControlModifier) == ["a", "d"]
    assert tree.opened[-1] == "b"  # Files panel still follows the clicked folder


def test_plain_click_resets_and_shift_click_extends(tree):
    _click(tree, "a")
    _click(tree, "c", M.ControlModifier)
    assert _click(tree, "b") == ["b"]
    assert _click(tree, "d", M.ShiftModifier) == ["b", "c", "d"]


def test_arrow_keys_still_move_a_single_selection(tree):
    _click(tree, "a")
    QTest.keyClick(tree, Qt.Key.Key_Down)
    QTest.qWait(10)
    assert _selected(tree) == ["b"]
    assert tree.opened[-1] == "b"
    QTest.keyClick(tree, Qt.Key.Key_Down, M.ShiftModifier)
    QTest.qWait(10)
    assert _selected(tree) == ["b", "c"]


def test_programmatic_open_selects_exactly_that_folder(tree, tmp_path):
    _click(tree, "a")
    _click(tree, "c", M.ControlModifier)
    tree.folder_selected.emit(tmp_path / "d")  # e.g. history / go-to folder
    assert _selected(tree) == ["d"]


def test_multi_selection_offers_deleting_all(tree, monkeypatch):
    labels = []

    class RecordingMenu(mw.QMenu):
        def exec(self, *args, **kwargs):
            labels.extend(a.text() for a in self.actions() if a.text())

    monkeypatch.setattr(mw, "QMenu", RecordingMenu)
    _click(tree, "a")
    _click(tree, "c", M.ControlModifier)
    item_c = tree.topLevelItem(0).child(2)
    tree._show_context_menu(tree.visualItemRect(item_c).center())
    assert "Delete 2 folders" in labels
