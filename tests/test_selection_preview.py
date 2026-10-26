"""Temporary combined FASTA used to preview several selected files together."""

from __future__ import annotations

from pathlib import Path

import pytest

from repeatviewer.backends import selection_preview as sp


@pytest.fixture()
def cache_root(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    root = tmp_path / "cache" / "selection_preview"
    monkeypatch.setattr(sp, "selection_preview_root", lambda: root)
    return root


def _fasta(folder: Path, name: str, records: list[tuple[str, str]]) -> Path:
    path = folder / name
    path.write_text("".join(f">{rid}\n{seq}\n" for rid, seq in records), encoding="utf-8")
    return path


def _ids(path: Path) -> list[str]:
    return [line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")]


def test_combined_file_keeps_selection_order(tmp_path: Path, cache_root: Path) -> None:
    a = _fasta(tmp_path, "CroSat3.fasta", [("CroSat3", "ACGTACGT")])
    b = _fasta(tmp_path, "CroSat2.fasta", [("CroSat2", "GGGGCCCC")])
    out = sp.build_selection_preview([a, b])
    assert _ids(out) == ["CroSat3", "CroSat2"]
    assert out.parent.parent == cache_root
    assert sp.is_selection_preview_path(out)
    assert sp.selection_preview_sources(out) == [a.resolve(), b.resolve()]


def test_nothing_is_written_next_to_the_sources(tmp_path: Path, cache_root: Path) -> None:
    folder = tmp_path / "data"
    folder.mkdir()
    a = _fasta(folder, "a.fa", [("a", "ACGT")])
    b = _fasta(folder, "b.fa", [("b", "TTTT")])
    before = sorted(p.name for p in folder.iterdir())
    sp.build_selection_preview([a, b])
    after = sorted(p.name for p in folder.iterdir())
    assert after == before


def test_unchanged_selection_reuses_the_cached_file(tmp_path: Path, cache_root: Path) -> None:
    a = _fasta(tmp_path, "a.fa", [("a", "ACGT")])
    b = _fasta(tmp_path, "b.fa", [("b", "TTTT")])
    first = sp.build_selection_preview([a, b])
    mtime = first.stat().st_mtime_ns
    second = sp.build_selection_preview([a, b])
    assert second == first
    assert second.stat().st_mtime_ns == mtime


def test_edited_source_gives_a_new_combined_file(tmp_path: Path, cache_root: Path) -> None:
    a = _fasta(tmp_path, "a.fa", [("a", "ACGT")])
    b = _fasta(tmp_path, "b.fa", [("b", "TTTT")])
    first = sp.build_selection_preview([a, b])
    _fasta(tmp_path, "a.fa", [("a", "ACGTACGTAA")])
    second = sp.build_selection_preview([a, b])
    assert second != first
    assert "ACGTACGTAA" in second.read_text()


def test_duplicate_ids_stay_unique(tmp_path: Path, cache_root: Path) -> None:
    a = _fasta(tmp_path, "a.fa", [("seq", "ACGT")])
    b = _fasta(tmp_path, "b.fa", [("seq", "TTTT")])
    ids = _ids(sp.build_selection_preview([a, b]))
    assert len(ids) == len(set(ids)) == 2


def test_one_file_is_not_a_selection_preview(tmp_path: Path, cache_root: Path) -> None:
    a = _fasta(tmp_path, "a.fa", [("a", "ACGT")])
    with pytest.raises(ValueError):
        sp.build_selection_preview([a])
    assert not sp.is_selection_preview_path(a)
    assert not sp.is_selection_preview_path(None)


def test_limit_message(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    a = _fasta(tmp_path, "a.fa", [("a", "ACGT" * 100)])
    b = _fasta(tmp_path, "b.fa", [("b", "TTTT")])
    assert sp.selection_preview_limit_message([a, b]) is None
    monkeypatch.setattr(sp, "SELECTION_PREVIEW_MAX_FILES", 1)
    assert "at most 1" in (sp.selection_preview_limit_message([a, b]) or "")
    monkeypatch.setattr(sp, "SELECTION_PREVIEW_MAX_FILES", 100)
    monkeypatch.setattr(sp, "SELECTION_PREVIEW_MAX_BYTES", 10)
    assert "too large" in (sp.selection_preview_limit_message([a, b]) or "")
