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


def test_labels_for_paste_move_duplicate_rename():
    assert browser_undo_label({"kind": "transfer", "action": "Paste", "count": 2}) == "Undo Paste 2 Files"
    assert browser_undo_label({"kind": "transfer", "action": "Move", "count": 1}) == "Undo Move 1 File"
    assert browser_undo_label({"kind": "transfer", "action": "Duplicate", "count": 3}) == "Undo Duplicate 3 Files"
    assert browser_undo_label({"kind": "created_folders", "items": [1]}) == "Undo Duplicate 1 Folder"
    assert browser_undo_label({"kind": "rename_file"}) == "Undo Rename File"
    assert browser_undo_label({"kind": "rename_folder"}) == "Undo Rename Folder"


def test_fingerprint_notices_changes_but_not_finder_files(tmp_path):
    import os

    from repeatviewer.edit_actions import path_fingerprint

    folder = tmp_path / "copy"
    (folder / "sub").mkdir(parents=True)
    (folder / "a.fa").write_text(">a\nACGT\n")
    (folder / "sub" / "b.fa").write_text(">b\nAAAA\n")
    before = path_fingerprint(folder)
    assert before is not None and before == path_fingerprint(folder)

    (folder / ".DS_Store").write_bytes(b"finder")          # Finder opened it
    assert path_fingerprint(folder) == before

    (folder / "sub" / "b.fa").write_text(">b\nAAAAC\n")    # edited
    assert path_fingerprint(folder) != before

    single = tmp_path / "x.fa"
    single.write_text(">x\nA\n")
    fp = path_fingerprint(single)
    os.utime(single, ns=(1, 1))
    assert path_fingerprint(single) != fp
    assert path_fingerprint(tmp_path / "missing") is None
