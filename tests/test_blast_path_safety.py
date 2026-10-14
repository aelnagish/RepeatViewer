"""BLAST path helpers must not treat path separators as unsafe."""

from __future__ import annotations

import os
from pathlib import Path

from repeatviewer.backends.gui_functions import (
    _blast_path_has_unsafe_chars,
    _blast_safe_fasta_for_makeblastdb,
)


def test_absolute_unix_path_without_spaces_is_safe() -> None:
    assert not _blast_path_has_unsafe_chars("/var/folders/xy/tmp/seq.fa")
    assert not _blast_path_has_unsafe_chars("/tmp/repeatviewer_msa_annot_abc/ungapped.fa")


def test_space_and_hash_are_unsafe() -> None:
    assert _blast_path_has_unsafe_chars("/tmp/my project/seq.fa")
    assert _blast_path_has_unsafe_chars("/tmp/seq#1.fa")


def test_parallel_ungapped_links_do_not_collide(tmp_path: Path) -> None:
    """Two files both named ungapped.fa must get distinct makeblastdb -in links."""
    spaced = tmp_path / "project dir"
    spaced.mkdir()
    a_dir = spaced / "a"
    b_dir = spaced / "b"
    a_dir.mkdir()
    b_dir.mkdir()
    a = a_dir / "ungapped.fa"
    b = b_dir / "ungapped.fa"
    a.write_text(">a\nACGT\n", encoding="utf-8")
    b.write_text(">b\nTTTT\n", encoding="utf-8")

    in_a, cleanup_a = _blast_safe_fasta_for_makeblastdb(str(a))
    in_b, cleanup_b = _blast_safe_fasta_for_makeblastdb(str(b))
    try:
        assert in_a != in_b
        assert os.path.exists(in_a)
        assert os.path.exists(in_b)
        assert Path(in_a).read_text(encoding="utf-8") == a.read_text(encoding="utf-8")
        assert Path(in_b).read_text(encoding="utf-8") == b.read_text(encoding="utf-8")
    finally:
        for path in (cleanup_a, cleanup_b):
            if path and os.path.lexists(path):
                os.remove(path)
