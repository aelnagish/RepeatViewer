from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.safe_cli_paths import (
    basename_has_spaces,
    find_names_with_spaces,
    rename_spaces_in_basename,
    rename_spaced_paths,
    spaced_basename_target,
)


def test_basename_has_spaces() -> None:
    assert basename_has_spaces("foo bar.fa") is True
    assert basename_has_spaces("foo_bar.fa") is False
    assert basename_has_spaces(Path("/tmp/a b/c.fa")) is False


def test_spaced_basename_target_replaces_spaces(tmp_path: Path) -> None:
    path = tmp_path / "my seq.fa"
    path.write_text(">a\nACGT\n")
    target = spaced_basename_target(path)
    assert target.name == "my_seq.fa"
    assert target.parent == path.parent


def test_spaced_basename_target_avoids_collision(tmp_path: Path) -> None:
    existing = tmp_path / "my_seq.fa"
    existing.write_text(">a\nACGT\n")
    path = tmp_path / "my seq.fa"
    path.write_text(">b\nTTTT\n")
    target = spaced_basename_target(path)
    assert target.name == "my_seq_1.fa"


def test_rename_spaces_in_basename(tmp_path: Path) -> None:
    path = tmp_path / "TE family.fa"
    path.write_text(">a\nACGT\n")
    final, renamed = rename_spaces_in_basename(path)
    assert renamed is True
    assert final.name == "TE_family.fa"
    assert final.exists()
    assert not path.exists()


def test_find_and_rename_spaced_paths(tmp_path: Path) -> None:
    spaced_dir = tmp_path / "TEtrimmer output"
    spaced_dir.mkdir()
    fasta = spaced_dir / "LTR element.fa"
    fasta.write_text(">a\nACGT\n")
    (spaced_dir / "ok.fa").write_text(">b\nGGGG\n")

    hits, truncated = find_names_with_spaces(tmp_path)
    assert truncated is False
    names = {h.path.name for h in hits}
    assert "TEtrimmer output" in names
    assert "LTR element.fa" in names
    assert "ok.fa" not in names

    pairs, errors = rename_spaced_paths([h.path for h in hits])
    assert errors == []
    assert pairs
    assert not (tmp_path / "TEtrimmer output").exists()
    assert (tmp_path / "TEtrimmer_output" / "LTR_element.fa").exists()
