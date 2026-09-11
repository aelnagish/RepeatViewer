"""Destination + name validation for the create-folder dialog.

Pure-function tests: no QApplication needed, so they run in any CI without a
display. The dialog widget only renders what ``validate_folder_name`` decides.
"""

from __future__ import annotations

from datetime import datetime
from pathlib import Path

from repeatviewer.views.create_folder_dialog import (
    common_prefix_name,
    expand_name_tokens,
    unique_folder_name,
    validate_folder_name,
)


def test_plain_name_is_accepted(tmp_path: Path) -> None:
    ok, message, suggestion = validate_folder_name(tmp_path, "run3")
    assert ok
    assert message == ""
    assert suggestion == ""


def test_empty_name_is_rejected(tmp_path: Path) -> None:
    for raw in ("", "   "):
        ok, message, _ = validate_folder_name(tmp_path, raw)
        assert not ok
        assert "Enter a folder name" in message


def test_existing_name_offers_the_next_free_name(tmp_path: Path) -> None:
    (tmp_path / "run3").mkdir()
    ok, message, suggestion = validate_folder_name(tmp_path, "run3")
    assert not ok
    assert "already exists" in message
    assert suggestion == "run3 2"

    (tmp_path / "run3 2").mkdir()
    _, _, suggestion = validate_folder_name(tmp_path, "run3")
    assert suggestion == "run3 3"


def test_reserved_and_illegal_names_are_rejected(tmp_path: Path) -> None:
    for raw in (".", "..", "a:b", "a\\b"):
        ok, _, _ = validate_folder_name(tmp_path, raw)
        assert not ok, raw


def test_absolute_paths_are_rejected(tmp_path: Path) -> None:
    for raw in ("/etc/passwd", "~/Desktop"):
        ok, message, _ = validate_folder_name(tmp_path, raw)
        assert not ok
        assert "absolute path" in message


def test_nested_name_is_allowed_and_announced(tmp_path: Path) -> None:
    ok, message, _ = validate_folder_name(tmp_path, "run3/monomers")
    assert ok
    assert "2 nested folders" in message


def test_empty_segment_is_rejected(tmp_path: Path) -> None:
    ok, message, _ = validate_folder_name(tmp_path, "run3//monomers")
    assert not ok
    assert "empty path segment" in message


def test_no_destination_is_rejected() -> None:
    ok, message, _ = validate_folder_name(None, "run3")
    assert not ok
    assert "Choose a folder" in message


def test_date_token_expands(tmp_path: Path) -> None:
    stamp = datetime(2026, 9, 11, 14, 30)
    assert expand_name_tokens("run_{date}", now=stamp) == "run_2026-09-11"
    assert expand_name_tokens("run_{date}_{time}", now=stamp) == "run_2026-09-11_1430"
    ok, message, _ = validate_folder_name(tmp_path, "run_{date}")
    assert ok
    assert "Creates" in message


def test_unique_folder_name_counts_up(tmp_path: Path) -> None:
    assert unique_folder_name(tmp_path, "x") == "x"
    (tmp_path / "x").mkdir()
    assert unique_folder_name(tmp_path, "x") == "x 2"


def test_common_prefix_name_groups_a_family() -> None:
    paths = [
        Path("CroSat2_2.features.gff3"),
        Path("CroSat2_2_flexidot.gff.colors.config"),
        Path("CroSat2_2_flexidot.html"),
        Path("CroSat2_2_flexidot.png"),
    ]
    assert common_prefix_name(paths) == "CroSat2_2"


def test_common_prefix_name_declines_weak_prefixes() -> None:
    assert common_prefix_name([Path("alpha.fa"), Path("beta.fa")]) == ""
    assert common_prefix_name([]) == ""
    assert common_prefix_name([Path("only_one.fa")]) == "only_one"
