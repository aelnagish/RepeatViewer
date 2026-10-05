"""Tests for shared tool-output path helpers."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.output_paths import next_numbered_path


def test_next_numbered_path_inserts_before_suffix(tmp_path: Path) -> None:
    path = tmp_path / "reads.fa_FastGA.paf"
    assert next_numbered_path(path) == path
    path.write_text("x")
    assert next_numbered_path(path).name == "reads.fa_FastGA_2.paf"
    (tmp_path / "reads.fa_FastGA_2.paf").write_text("y")
    assert next_numbered_path(path).name == "reads.fa_FastGA_3.paf"
