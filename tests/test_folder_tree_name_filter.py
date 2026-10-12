"""Folders panel name filter matches folder names and contained file names."""

from pathlib import Path

import pytest
from PySide6.QtWidgets import QTreeWidgetItem

mw = pytest.importorskip("repeatviewer.main_window")


def _add_folder(parent: QTreeWidgetItem, path: Path) -> QTreeWidgetItem:
    path.mkdir(parents=True, exist_ok=True)
    item = QTreeWidgetItem([path.name])
    item.setData(0, mw.ROLE_PATH, str(path))
    parent.addChild(item)
    return item


@pytest.fixture
def filtered_tree(qapp, tmp_path):
    tree = mw.FolderTreeWidget()
    root = QTreeWidgetItem([tmp_path.name])
    root.setData(0, mw.ROLE_PATH, str(tmp_path))
    tree.addTopLevelItem(root)

    keep = _add_folder(root, tmp_path / "keep_me")
    (tmp_path / "keep_me" / "target_seq.fa").write_text(">a\nAT\n")
    other = _add_folder(root, tmp_path / "other")
    (tmp_path / "other" / "readme.txt").write_text("nope")
    nested_parent = _add_folder(root, tmp_path / "nested_parent")
    nested = _add_folder(nested_parent, tmp_path / "nested_parent" / "deep")
    (tmp_path / "nested_parent" / "deep" / "target_hit.fasta").write_text(">b\nGC\n")

    root.setExpanded(True)
    nested_parent.setExpanded(True)
    tree.resize(300, 300)
    tree.show()
    qapp.processEvents()
    yield tree, {
        "root": root,
        "keep": keep,
        "other": other,
        "nested_parent": nested_parent,
        "nested": nested,
    }
    tree.close()


def _visible(items: dict[str, QTreeWidgetItem]) -> dict[str, bool]:
    return {name: not item.isHidden() for name, item in items.items()}


def test_name_filter_matches_file_names(filtered_tree):
    tree, items = filtered_tree
    matching_files = {
        Path(str(items["keep"].data(0, mw.ROLE_PATH))),
        Path(str(items["nested"].data(0, mw.ROLE_PATH))),
    }

    tree.apply_name_filter(
        "target",
        folder_has_matching_file=lambda folder: folder in matching_files,
    )

    shown = _visible(items)
    assert shown["keep"] is True
    assert shown["nested"] is True
    assert shown["nested_parent"] is True  # ancestor of file match
    assert shown["other"] is False
    assert shown["root"] is True
    assert items["nested_parent"].isExpanded()


def test_name_filter_still_matches_folder_names(filtered_tree):
    tree, items = filtered_tree
    tree.apply_name_filter("other", folder_has_matching_file=lambda _folder: False)
    shown = _visible(items)
    assert shown["other"] is True
    assert shown["keep"] is False
    assert shown["nested"] is False
    assert shown["nested_parent"] is False


def test_name_filter_clear_shows_all(filtered_tree):
    tree, items = filtered_tree
    tree.apply_name_filter("target", folder_has_matching_file=lambda _folder: True)
    tree.apply_name_filter("")
    assert all(not item.isHidden() for item in items.values())


def test_folder_has_file_name_matching_ignores_hidden(qapp, tmp_path):
    """MainWindow helper: match visible files only, not blast sidecars / dotfiles."""
    window = mw.MainWindow.__new__(mw.MainWindow)
    folder = tmp_path / "lib"
    folder.mkdir()
    (folder / "HitMe.fa").write_text(">x\nA\n")
    (folder / ".secret.fa").write_text(">y\nA\n")
    (folder / "db.nhr").write_bytes(b"x")

    assert window._folder_has_file_name_matching(folder, "hitme") is True
    assert window._folder_has_file_name_matching(folder, "secret") is False
    assert window._folder_has_file_name_matching(folder, "nhr") is False
    assert window._folder_has_file_name_matching(folder, "missing") is False
