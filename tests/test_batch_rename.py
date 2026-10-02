from __future__ import annotations

from pathlib import Path

import pytest

from repeatviewer.batch_rename import (
    RenameOptions,
    apply_moves,
    check_moves,
    name_problem,
    new_name_for,
    plan_batch_rename,
    split_extension,
)


def _touch(folder: Path, *names: str) -> list[Path]:
    paths = []
    for name in names:
        path = folder / name
        path.write_text(name)
        paths.append(path)
    return paths


@pytest.mark.parametrize(
    ("name", "expected"),
    [
        ("sample.fa", ("sample", ".fa")),
        ("sample.fa.gz", ("sample", ".fa.gz")),
        ("reads.fastq.gz", ("reads", ".fastq.gz")),
        ("sample.v2.fa", ("sample.v2", ".fa")),
        ("x.fa_b.txt.fa", ("x.fa_b.txt", ".fa")),
        (".hidden", (".hidden", "")),
        ("README", ("README", "")),
        ("odd.name with space", ("odd.name with space", "")),
    ],
)
def test_split_extension(name, expected):
    assert split_extension(name) == expected


def test_new_name_find_replace_prefix_suffix():
    opts = RenameOptions(find="contig", replace="Chr", prefix="Cs_", suffix="_v2")
    assert new_name_for("contig_1.fa", opts) == "Cs_Chr_1_v2.fa"
    assert new_name_for("CONTIG_2.fa.gz", opts) == "Cs_Chr_2_v2.fa.gz"


def test_new_name_match_case_and_extension_untouched():
    opts = RenameOptions(find="fa", replace="X", match_case=True)
    assert new_name_for("fa_FA.fa", opts) == "X_FA.fa"


def test_new_name_without_keep_extension():
    opts = RenameOptions(find=".fasta", replace=".fa", keep_extension=False, suffix="")
    assert new_name_for("a.fasta", opts) == "a.fa"


def test_new_name_regex_groups():
    opts = RenameOptions(find=r"^seq_(\d+)$", replace=r"CroSat_\1", regex=True)
    assert new_name_for("seq_12.fa", opts) == "CroSat_12.fa"
    assert new_name_for("other_12.fa", opts) == "other_12.fa"


def test_plan_reports_broken_regex(tmp_path):
    paths = _touch(tmp_path, "a.fa")
    plan = plan_batch_rename(paths, RenameOptions(find="(", regex=True))
    assert plan.error.startswith("Regular expression")
    assert not plan.ok


def test_name_problem():
    assert name_problem("") != ""
    assert name_problem("a/b.fa") != ""
    assert name_problem("ok.fa") == ""
    assert name_problem("bad?.fa", windows=True) != ""
    assert name_problem("x" * 300) != ""


def test_plan_detects_existing_and_duplicate_names(tmp_path):
    a, b, c = _touch(tmp_path, "a_1.fa", "b_1.fa", "keep.fa")
    _touch(tmp_path, "a_2.fa")
    plan = plan_batch_rename(
        [a, b, c], RenameOptions(find="_1", replace="_2"), case_insensitive=False
    )
    rows = {row.source.name: row for row in plan.rows}
    assert rows["a_1.fa"].problem == "a file with this name already exists"
    assert rows["b_1.fa"].problem == ""
    assert rows["keep.fa"].changed is False
    assert not plan.ok

    plan = plan_batch_rename(
        [a, b], RenameOptions(find=r"^[ab]", replace="x", regex=True), case_insensitive=False
    )
    assert all(row.problem == "same new name as another file" for row in plan.rows)


def test_plan_allows_swaps_and_chains(tmp_path):
    a, b = _touch(tmp_path, "a.fa", "b.fa")
    # a -> b and b -> c in one go: b is free once it moved away.
    opts = RenameOptions(find=r"^b$", replace="c", regex=True)
    plan = plan_batch_rename([b], opts, case_insensitive=False)
    assert plan.ok
    plan = plan_batch_rename(
        [a, b],
        RenameOptions(find=r"^(a|b)$", replace=r"\1x", regex=True),
        case_insensitive=False,
    )
    assert plan.ok


def test_plan_case_only_rename_on_case_insensitive_disk(tmp_path):
    (a,) = _touch(tmp_path, "sample.fa")
    plan = plan_batch_rename(
        [a],
        RenameOptions(find="sample", replace="Sample", match_case=True),
        case_insensitive=True,
    )
    assert plan.ok


def test_plan_clash_with_kept_name_cascades(tmp_path):
    a, b, c = _touch(tmp_path, "a.fa", "b.fa", "c.fa")
    # a and b both want "x"; neither moves, so c must not take "a".
    opts = RenameOptions(find=r"^(a|b|c)$", replace="x", regex=True)
    plan = plan_batch_rename([a, b], opts, case_insensitive=False)
    assert all(row.problem for row in plan.rows)


def test_plan_companions_follow_and_are_checked(tmp_path):
    fa, side = _touch(tmp_path, "s.fa", "s.seqmap_RV.json")
    _touch(tmp_path, "t.seqmap_RV.json")

    def companion_files(path: Path) -> list[Path]:
        if path.suffix == ".fa":
            json = path.parent / f"{path.stem}.seqmap_RV.json"
            return [json] if json.exists() else []
        return []

    def companion_renames(old: Path, new: Path) -> list[tuple[Path, Path]]:
        return [
            (src, new.parent / f"{new.stem}.seqmap_RV.json") for src in companion_files(old)
        ]

    plan = plan_batch_rename(
        [fa, side],
        RenameOptions(find="s", replace="u"),
        companion_files=companion_files,
        companion_renames=companion_renames,
        case_insensitive=False,
    )
    assert [row.source.name for row in plan.rows] == ["s.fa"]
    assert plan.ok
    assert (side, tmp_path / "u.seqmap_RV.json") in plan.moves

    plan = plan_batch_rename(
        [fa],
        RenameOptions(find="s", replace="t"),
        companion_files=companion_files,
        companion_renames=companion_renames,
        case_insensitive=False,
    )
    assert "annotation file t.seqmap_RV.json already exists" in plan.rows[0].problem


def test_apply_moves_swap_and_rollback(tmp_path):
    a, b = _touch(tmp_path, "a.fa", "b.fa")
    assert apply_moves([(a, b), (b, a)]) == ""
    assert a.read_text() == "b.fa" and b.read_text() == "a.fa"

    c = tmp_path / "c.fa"
    blocker = tmp_path / "taken.fa"
    blocker.write_text("x")
    error = apply_moves([(a, c), (b, blocker)])
    assert "taken.fa" in error
    # Nothing changed, nothing overwritten, no temporary files left behind.
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a.fa", "b.fa", "taken.fa"]
    assert blocker.read_text() == "x"


def test_check_moves(tmp_path):
    a, b = _touch(tmp_path, "a.fa", "b.fa")
    assert check_moves([(a, b), (b, a)], case_insensitive=False) == []
    assert check_moves([(a, b)], case_insensitive=False) == ["b.fa: this name is in use"]
    assert check_moves([(tmp_path / "gone.fa", tmp_path / "z.fa")], case_insensitive=False)


def test_prefix_chain_renames_in_two_steps(tmp_path):
    a, pa = _touch(tmp_path, "a.fa", "pa.fa")
    plan = plan_batch_rename([a, pa], RenameOptions(prefix="p"), case_insensitive=False)
    assert plan.ok
    assert apply_moves(plan.moves) == ""
    assert (tmp_path / "pa.fa").read_text() == "a.fa"
    assert (tmp_path / "ppa.fa").read_text() == "pa.fa"
