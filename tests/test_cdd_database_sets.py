"""Recognize full CDD and NCBI subset / LE packages."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.database_prepare import (
    cdd_database_present,
    cdd_index_is_ready,
    describe_database_problem,
    find_cdd_pn,
    find_cdd_preformatted_prefix,
    prepared_cdd_prefix,
    suggest_cdd_database_name,
)
from repeatviewer.views.define_database_dialog import resolve_panel_database_path


def test_source_pfam_pn_is_recognized(tmp_path: Path) -> None:
    (tmp_path / "Pfam.pn").write_text("pfam00001.smp\n", encoding="utf-8")
    assert cdd_database_present(tmp_path)
    assert find_cdd_pn(tmp_path) == (tmp_path / "Pfam.pn").resolve()
    assert suggest_cdd_database_name(tmp_path) == "Pfam"
    assert not cdd_index_is_ready(tmp_path)
    assert prepared_cdd_prefix(tmp_path) is None
    assert describe_database_problem("cdd", tmp_path) is None


def test_source_cdd_pn_preferred_over_subset(tmp_path: Path) -> None:
    (tmp_path / "Pfam.pn").write_text("x\n", encoding="utf-8")
    (tmp_path / "Cdd.pn").write_text("y\n", encoding="utf-8")
    assert find_cdd_pn(tmp_path).name == "Cdd.pn"
    assert suggest_cdd_database_name(tmp_path) == "Cdd"


def test_preformatted_pfam_pal_is_ready(tmp_path: Path) -> None:
    (tmp_path / "Pfam.pal").write_text('DBLIST "Pfam.00"\n', encoding="utf-8")
    assert cdd_database_present(tmp_path)
    assert cdd_index_is_ready(tmp_path)
    assert find_cdd_preformatted_prefix(tmp_path) == (tmp_path / "Pfam").resolve()
    assert prepared_cdd_prefix(tmp_path) == (tmp_path / "Pfam").resolve()
    assert suggest_cdd_database_name(tmp_path) == "Pfam"
    assert describe_database_problem("cdd", tmp_path) is None


def test_empty_folder_reports_expected_markers(tmp_path: Path) -> None:
    problem = describe_database_problem("cdd", tmp_path)
    assert problem is not None
    assert "Cdd.pn" in problem
    assert "Pfam.pal" in problem or "Cdd.pal" in problem


def test_resolve_panel_path_accepts_pfam_pal(tmp_path: Path) -> None:
    pal = tmp_path / "Pfam.pal"
    pal.write_text("x\n", encoding="utf-8")
    assert resolve_panel_database_path(pal) == ("cdd", tmp_path.resolve())


def test_resolve_panel_path_accepts_smart_pn(tmp_path: Path) -> None:
    pn = tmp_path / "Smart.pn"
    pn.write_text("x\n", encoding="utf-8")
    assert resolve_panel_database_path(pn) == ("cdd", tmp_path.resolve())
