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


def test_batch_rename_and_redo_labels():
    from repeatviewer.edit_actions import browser_redo_label

    assert browser_undo_label({"kind": "batch_rename", "count": 12}) == "Undo Rename 12 Files"
    assert browser_redo_label(None) == "Redo"
    assert browser_redo_label({"kind": "something-new"}) == "Redo"
    assert browser_redo_label({"kind": "files", "count": 2}) == "Redo Delete 2 Files"
    assert browser_redo_label({"kind": "transfer", "action": "Paste", "count": 1}) == "Redo Paste 1 File"
    assert browser_redo_label({"kind": "hide_folders", "count": 1}) == "Redo Remove 1 Folder from List"
    assert browser_redo_label({"kind": "created_folders", "count": 2}) == "Redo Duplicate 2 Folders"
    assert browser_redo_label({"kind": "rename_file", "count": 1}) == "Redo Rename File"
    assert browser_redo_label({"kind": "batch_rename", "count": 3}) == "Redo Rename 3 Files"


def test_reverse_complement_text_plain_and_wrapped():
    from repeatviewer.edit_actions import reverse_complement_text

    assert reverse_complement_text("ACGTN") == ("NACGT", 5)
    assert reverse_complement_text("acgRY\n") == ("RYcgt\n", 5)
    # Line width of the original is kept.
    assert reverse_complement_text("AAAACCCC\nGGGGTT\n") == ("AACCCCGG\nGGTTTT\n", 14)
    # Numbers and spaces (GenBank style) are ignored, gaps kept.
    assert reverse_complement_text("  1 acgt acgt\n  9 aa") == ("ttacgtac\ngt", 10)
    assert reverse_complement_text("AC-GT") == ("AC-GT", 4)
    # RNA stays RNA.
    assert reverse_complement_text("UUAG") == ("CUAA", 4)


def test_reverse_complement_text_fasta_and_rejects():
    from repeatviewer.edit_actions import reverse_complement_text

    text = ">s1 desc\nAAAC\nGG\n>s2\nUUAG"
    assert reverse_complement_text(text) == (
        ">s1 desc reverse_complement\nCCGT\nTT\n>s2 reverse_complement\nCUAA",
        10,
    )
    assert reverse_complement_text("MKLVE") is None      # protein
    assert reverse_complement_text("   ") is None
    assert reverse_complement_text(">only header") is None
