"""Folder Duplicate: naming, size check, copy, Stop and error clean-up."""

from __future__ import annotations

import os
from pathlib import Path

import pytest

from repeatviewer.folder_duplicate import (
    CopyCancelled,
    copy_folder,
    format_size,
    measure_folders,
    outermost_folders,
    unique_folder_copy_path,
)


def _make_tree(root: Path) -> Path:
    src = root / "run2.fasta_out"
    (src / "sub" / "deep").mkdir(parents=True)
    (src / "a.fa").write_text(">a\nACGT\n")
    (src / "sub" / "b.fa").write_text(">b\n" + "A" * 1000 + "\n")
    (src / "sub" / "deep" / "c.txt").write_text("x" * 10)
    return src


def test_copy_name_keeps_the_whole_folder_name(tmp_path):
    src = _make_tree(tmp_path)
    assert unique_folder_copy_path(src).name == "run2.fasta_out_copy"
    (tmp_path / "run2.fasta_out_copy").mkdir()
    assert unique_folder_copy_path(src).name == "run2.fasta_out_copy2"
    (tmp_path / "run2.fasta_out_copy2").mkdir()
    assert unique_folder_copy_path(src).name == "run2.fasta_out_copy3"


def test_nested_selection_is_copied_once(tmp_path):
    src = _make_tree(tmp_path)
    other = tmp_path / "other"
    other.mkdir()
    picked = outermost_folders([src / "sub", src, other, src])
    assert picked == [src.resolve(), other.resolve()]


def test_measure_counts_files_and_stops_early(tmp_path):
    src = _make_tree(tmp_path)
    size = measure_folders([src])
    assert size.complete and size.files == 3
    assert size.bytes == sum(p.stat().st_size for p in src.rglob("*") if p.is_file())
    early = measure_folders([src], stop_above_bytes=100)
    assert not early.complete and early.bytes > 100
    limited = measure_folders([src], entry_limit=2)
    assert not limited.complete


def test_measure_does_not_follow_symlinks(tmp_path):
    src = _make_tree(tmp_path)
    big = tmp_path / "big.bin"
    big.write_bytes(b"0" * 50_000)
    os.symlink(big, src / "link.bin")
    os.symlink(tmp_path, src / "loop")  # would recurse forever if followed
    size = measure_folders([src])
    assert size.complete and size.files == 3 and size.bytes < 50_000


def test_copy_folder_copies_everything_and_reports_progress(tmp_path):
    src = _make_tree(tmp_path)
    os.symlink("a.fa", src / "link.fa")
    dest = unique_folder_copy_path(src)
    seen = []
    copy_folder(src, dest, progress=seen.append)
    assert (dest / "sub" / "deep" / "c.txt").read_text() == "x" * 10
    assert (dest / "link.fa").is_symlink()          # links stay links
    assert sum(seen) == measure_folders([src]).bytes
    assert (src / "a.fa").exists()                  # source untouched


def test_stop_removes_the_partial_copy(tmp_path):
    src = _make_tree(tmp_path)
    dest = unique_folder_copy_path(src)
    calls = {"n": 0}

    def cancelled():
        calls["n"] += 1
        return calls["n"] > 1  # stop after the first file

    with pytest.raises(CopyCancelled):
        copy_folder(src, dest, cancelled=cancelled)
    assert not dest.exists()
    assert (src / "sub" / "b.fa").exists()


def test_existing_destination_is_never_overwritten(tmp_path):
    src = _make_tree(tmp_path)
    dest = tmp_path / "taken"
    dest.mkdir()
    (dest / "keep.txt").write_text("mine")
    with pytest.raises(FileExistsError):
        copy_folder(src, dest)
    assert (dest / "keep.txt").read_text() == "mine"


def test_format_size():
    assert format_size(512) == "512 bytes"
    assert format_size(2048) == "2.0 KB"
    assert format_size(5 * 1024**2) == "5.0 MB"
    assert format_size(3 * 1024**3) == "3.0 GB"
    assert format_size(2048 * 1024**3) == "2048.0 GB"
