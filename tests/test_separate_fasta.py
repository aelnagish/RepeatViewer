from pathlib import Path

from repeatviewer.backends.separate_fasta import (
    sanitize_name,
    separate_sequences,
    separation_output_dir,
)


def test_sanitize_name_replaces_hash_with_underscore() -> None:
    assert sanitize_name("Ama5#Class_I__SINE") == "Ama5_Class_I__SINE"
    assert "#" not in sanitize_name("foo#bar#baz")


def test_separate_sequences_sanitizes_hash_in_filenames(tmp_path: Path) -> None:
    src = tmp_path / "multi.fa"
    src.write_text(
        ">Ama5-lic_Beta#Class_I__SINE\nACGT\n>other/name\nTTTT\n",
        encoding="utf-8",
    )
    out = tmp_path / "Separated_TE_sequence"
    created, n = separate_sequences(src, out, write_mapping=True)
    assert n == 2
    names = {p.name for p in created}
    assert "Ama5-lic_Beta_Class_I__SINE.fa" in names
    assert "other__name.fa" in names
    assert all("#" not in p.name for p in created)


def test_separate_sequences_subfolder_sanitizes_hash(tmp_path: Path) -> None:
    src = tmp_path / "multi.fa"
    src.write_text(">TE#LTR/Gypsy\nAAAA\n", encoding="utf-8")
    out = tmp_path / "out"
    created, n = separate_sequences(
        src, out, write_mapping=False, subfolder_per_sequence=True
    )
    assert n == 1
    assert created[0].name == "TE_LTR__Gypsy.fa"
    assert created[0].parent.name == "TE_LTR__Gypsy"
    assert "#" not in str(created[0])


def test_separate_sequences_preserves_faa_extension(tmp_path: Path) -> None:
    src = tmp_path / "proteins.faa"
    src.write_text(">prot1\nMKTAYI\n>prot2\nACDEFG\n", encoding="utf-8")
    out = tmp_path / "sep"
    created, n = separate_sequences(src, out, write_mapping=False)
    assert n == 2
    assert {p.suffix for p in created} == {".faa"}
    assert {p.name for p in created} == {"prot1.faa", "prot2.faa"}
    src = tmp_path / "my file.fa"
    src.write_text(">a\nAAAA\n", encoding="utf-8")
    out = separation_output_dir(src)
    assert out.name == "my_file_sep"


def test_separation_output_dir_conflict_adds_number(tmp_path: Path) -> None:
    src = tmp_path / "multi.fa"
    src.write_text(">a\nAAAA\n", encoding="utf-8")
    (tmp_path / "multi_sep").mkdir()
    out = separation_output_dir(src)
    assert out.name == "multi_sep_2"
