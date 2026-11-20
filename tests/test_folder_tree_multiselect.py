"""Folders panel: Cmd/Ctrl+click and Shift+click keep a multi-selection.

Opening a folder refreshes the Files panel, which calls back into the tree
(select_folder_item / set_active_folder). Those used to reset the selection
while Qt was still handling the click, so Cmd+click emptied it.
"""

from pathlib import Path

import pytest
from PySide6.QtCore import Qt, QPoint
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


def test_folder_context_menu_includes_reveal(tree, monkeypatch):
    labels = []

    class RecordingMenu(mw.QMenu):
        def exec(self, *args, **kwargs):
            labels.extend(a.text() for a in self.actions() if a.text())

    monkeypatch.setattr(mw, "QMenu", RecordingMenu)
    _click(tree, "a")
    item_a = tree.topLevelItem(0).child(0)
    tree._show_context_menu(tree.visualItemRect(item_a).center())
    assert "Reveal in Finder" in labels or "Show in File Manager" in labels
    assert "Copy Path" in labels


def test_top_level_roots_can_move_up_down(qapp, tmp_path, monkeypatch):
    """Main workspace folders can change position (Move up / Move down)."""
    tree = mw.FolderTreeWidget()
    ordered: list[list[str]] = []
    tree.root_order_changed.connect(
        lambda: ordered.append(
            [tree.topLevelItem(i).text(0) for i in range(tree.topLevelItemCount())]
        )
    )
    for name in ("first", "second", "third"):
        (tmp_path / name).mkdir()
        item = QTreeWidgetItem([name])
        item.setData(0, mw.ROLE_PATH, str(tmp_path / name))
        tree.addTopLevelItem(item)
    tree.resize(300, 300)
    tree.show()
    qapp.processEvents()

    second = tree.topLevelItem(1)
    tree._move_top_level_item(second, -1)
    assert [tree.topLevelItem(i).text(0) for i in range(3)] == [
        "second",
        "first",
        "third",
    ]
    assert ordered[-1] == ["second", "first", "third"]

    labels = []

    class RecordingMenu(mw.QMenu):
        def exec(self, *args, **kwargs):
            labels.extend(a.text() for a in self.actions() if a.text())

    monkeypatch.setattr(mw, "QMenu", RecordingMenu)
    tree.clearSelection()
    tree.topLevelItem(0).setSelected(True)
    tree._show_context_menu(tree.visualItemRect(tree.topLevelItem(0)).center())
    assert "Move up" in labels and "Move down" in labels
    tree.close()


def test_top_level_drag_reorder_keeps_both_roots(qapp, tmp_path):
    """Dragging a main folder below another must not make it disappear."""

    class _FakeDrop:
        def __init__(self, pos: QPoint) -> None:
            self._pos = pos
            self._action = Qt.DropAction.IgnoreAction
            self._accepted = False

        def position(self):
            return type("P", (), {"toPoint": lambda _s=None: self._pos})()

        def setDropAction(self, action) -> None:  # noqa: N802
            self._action = action

        def accept(self) -> None:
            self._accepted = True

        def ignore(self) -> None:
            self._accepted = False

    tree = mw.FolderTreeWidget()
    for name in ("TEviewer", "database"):
        (tmp_path / name).mkdir()
        item = QTreeWidgetItem([name])
        item.setData(0, mw.ROLE_PATH, str(tmp_path / name))
        tree.addTopLevelItem(item)
    # Nested child must not count as a reorderable root.
    nested = QTreeWidgetItem(["AnnoFlexi"])
    nested.setData(0, mw.ROLE_PATH, str(tmp_path / "TEviewer" / "AnnoFlexi"))
    nested.setFlags(nested.flags() & ~Qt.ItemFlag.ItemIsDragEnabled)
    tree.topLevelItem(0).addChild(nested)
    tree.topLevelItem(0).setExpanded(True)
    tree.resize(320, 360)
    tree.show()
    qapp.processEvents()

    te = tree.topLevelItem(0)
    db = tree.topLevelItem(1)
    # Nested selection alone is not a top-level drag source.
    tree.clearSelection()
    nested.setSelected(True)
    assert tree._selected_top_level_items() == []

    tree.clearSelection()
    te.setSelected(True)
    # Drop in the lower half of the database block → after database.
    drop_pos = tree.visualItemRect(db).center()
    drop_pos.setY(tree.visualItemRect(db).bottom() - 1)
    tree._set_drop_line_y(tree._insertion_y_for_top_level_drop(drop_pos))
    assert tree._drop_line_y is not None

    assert tree._reorder_top_level_drop(_FakeDrop(drop_pos)) is True
    assert tree.topLevelItemCount() == 2
    assert [tree.topLevelItem(i).text(0) for i in range(2)] == ["database", "TEviewer"]
    # Children stay with the moved root; nothing left selected after the move.
    moved = tree.topLevelItem(1)
    assert moved is not None and moved.childCount() == 1
    assert moved.child(0).text(0) == "AnnoFlexi"
    assert tree.selectedItems() == []
    tree.close()


def test_move_up_keeps_active_files_folder(qapp, tmp_path):
    """Move up/down must not steal the Files-panel highlight from another root."""
    tree = mw.FolderTreeWidget()
    items = {}
    for name in ("data", "RepeatMasker_output", "database"):
        (tmp_path / name).mkdir()
        item = QTreeWidgetItem([name])
        item.setData(0, mw.ROLE_PATH, str(tmp_path / name))
        tree.addTopLevelItem(item)
        items[name] = item
    tree.resize(320, 360)
    tree.show()
    qapp.processEvents()

    data_path = tmp_path / "data"
    tree.set_active_folder(data_path)
    assert tree._active_folder_key == str(data_path.resolve())

    # Move RepeatMasker up while browsing data — data stays the active folder.
    tree._move_top_level_item(items["RepeatMasker_output"], -1)
    assert [tree.topLevelItem(i).text(0) for i in range(3)] == [
        "RepeatMasker_output",
        "data",
        "database",
    ]
    assert tree._active_folder_key == str(data_path.resolve())
    assert tree.selectedItems() == []
    assert tree.currentItem() is items["data"]
    tree.close()


def _branch_item(tree, tmp_path, name: str):
    parent = tmp_path / name
    parent.mkdir()
    (parent / "child").mkdir()
    item = QTreeWidgetItem([name])
    item.setData(0, mw.ROLE_PATH, str(parent))
    item.setChildIndicatorPolicy(
        QTreeWidgetItem.ChildIndicatorPolicy.ShowIndicator
    )
    placeholder = QTreeWidgetItem([""])
    placeholder.setData(0, mw.ROLE_LAZY_PLACEHOLDER, True)
    item.addChild(placeholder)
    tree.topLevelItem(0).addChild(item)
    return item


def test_double_click_toggles_folder_expand_collapse(tree, tmp_path):
    """Label double-click expands, then collapses on the next double-click."""
    item = _branch_item(tree, tmp_path, "nested")
    item.setExpanded(False)
    tree.resize(300, 300)
    tree.show()

    label_pos = tree.visualItemRect(item).center()
    QTest.mouseDClick(tree.viewport(), Qt.MouseButton.LeftButton, M.NoModifier, label_pos)
    QTest.qWait(10)
    assert item.isExpanded() is True
    QTest.mouseDClick(tree.viewport(), Qt.MouseButton.LeftButton, M.NoModifier, label_pos)
    QTest.qWait(10)
    assert item.isExpanded() is False


def test_mixed_chevron_label_double_click_stays_expanded(tree, tmp_path):
    """Chevron press then label double-click must not collapse (old race)."""
    from PySide6.QtGui import QMouseEvent
    from PySide6.QtCore import QPointF, QEvent

    item = _branch_item(tree, tmp_path, "mixed")
    item.setExpanded(False)
    tree.resize(300, 300)
    tree.show()

    rect = tree.visualItemRect(item)
    chevron_pos = QPoint(rect.left() + 4, rect.center().y())
    label_pos = rect.center()

    def _press(pos: QPoint) -> None:
        ev = QMouseEvent(
            QEvent.Type.MouseButtonPress,
            QPointF(pos),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
        tree.mousePressEvent(ev)

    def _dbl(pos: QPoint) -> None:
        ev = QMouseEvent(
            QEvent.Type.MouseButtonDblClick,
            QPointF(pos),
            Qt.MouseButton.LeftButton,
            Qt.MouseButton.LeftButton,
            Qt.KeyboardModifier.NoModifier,
        )
        tree.mouseDoubleClickEvent(ev)

    _press(chevron_pos)  # expands via chevron
    assert item.isExpanded() is True
    _press(label_pos)  # second press of a double-click (must not toggle)
    _dbl(label_pos)  # keep chevron result; do not collapse
    assert item.isExpanded() is True
