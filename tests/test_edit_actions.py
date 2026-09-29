"""Edit menu helpers: the Undo label says what Undo will restore."""

from __future__ import annotations

from repeatviewer.edit_actions import browser_undo_label


def test_nothing_to_undo():
    assert browser_undo_label(None) == "Undo"
    assert browser_undo_label({}) == "Undo"
    assert browser_undo_label({"kind": "something-new"}) == "Undo"


def test_deleted_files():
    assert browser_undo_label({"kind": "files", "items": [1]}) == "Undo Delete 1 File"
    assert browser_undo_label({"kind": "files", "items": [1, 2, 3]}) == "Undo Delete 3 Files"


def test_deleted_folders():
    assert browser_undo_label({"kind": "delete_folders", "items": [1]}) == "Undo Delete 1 Folder"
    assert browser_undo_label({"kind": "delete_folders", "items": [1, 2]}) == "Undo Delete 2 Folders"


def test_folders_removed_from_list():
    assert (
        browser_undo_label({"kind": "hide_folders", "folders": ["a"]})
        == "Undo Remove 1 Folder from List"
    )
    assert (
        browser_undo_label({"kind": "hide_folders", "folders": ["a", "b"]})
        == "Undo Remove 2 Folders from List"
    )
