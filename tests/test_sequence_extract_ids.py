"""Extracting FASTA records named in an ID list (Extract Sequences, ID tab)."""

from __future__ import annotations

import gzip
from pathlib import Path

import pytest

from repeatviewer.backends.sequence_compose import (
    SequenceComposeError,
    count_matching_ids,
    extract_records_by_ids,
    parse_id_list,
    read_id_file,
)

FASTA = (
    ">read_1 length=4\nACGT\n"
    ">read_2\nTTTTGGGG\n"
    ">read_3 some description\nCCCC\n"
    ">read_2\nAAAA\n"  # duplicate header: legitimate in BLAST hit FASTAs
)


def _fasta(tmp_path: Path, name: str = "reads.fasta") -> Path:
    path = tmp_path / name
    path.write_text(FASTA, encoding="utf-8")
    return path


def _headers(path: Path) -> list[str]:
    return [
        line[1:].strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.startswith(">")
    ]


# --- parse_id_list ---------------------------------------------------------


def test_parse_id_list_ignores_blanks_comments_and_markers() -> None:
    parsed = parse_id_list("\n".join(["read_1", "", "# a comment", ">read_2", '"read_3"']))
    assert parsed.ids == ["read_1", "read_2", "read_3"]
    assert parsed.ignored == 2


def test_parse_id_list_handles_crlf_and_duplicates() -> None:
    parsed = parse_id_list("read_1\r\nread_2\r\nread_1\r\n")
    assert parsed.ids == ["read_1", "read_2"]
    assert parsed.duplicates == 1


def test_parse_id_list_takes_first_token_of_a_header() -> None:
    # BLAST reports the first whitespace token, so a pasted header still matches.
    assert parse_id_list(">read_3 some description").ids == ["read_3"]


def test_parse_id_list_reads_a_chosen_column_of_a_blast_table() -> None:
    table = (
        "qseqid\tsseqid\tpident\n"
        "monomer_1\tscaffold_5\t98.1\n"
        "monomer_1\tscaffold_7\t95.0\n"
        "monomer_2\tscaffold_5\t97.3\n"
    )
    subjects = parse_id_list(table, column=2)
    assert subjects.ids == ["scaffold_5", "scaffold_7"]
    assert subjects.duplicates == 1
    queries = parse_id_list(table, column=1)
    assert queries.ids == ["monomer_1", "monomer_2"]


def test_parse_id_list_skips_rows_without_the_wanted_column() -> None:
    parsed = parse_id_list("a\tb\nlonely\n", column=2)
    assert parsed.ids == ["b"]
    assert parsed.ignored == 1


def test_read_id_file(tmp_path: Path) -> None:
    ids = tmp_path / "ids.txt"
    ids.write_text("read_1\nread_3\n", encoding="utf-8")
    assert read_id_file(ids).ids == ["read_1", "read_3"]


# --- extraction ------------------------------------------------------------


def test_extract_writes_one_combined_file_in_fasta_order(tmp_path: Path) -> None:
    src = _fasta(tmp_path)
    summary = extract_records_by_ids(src, ["read_3", "read_1"], tmp_path)

    assert len(summary.written) == 1
    out = summary.written[0]
    assert out.name == "reads.fasta_extr.fa"
    # File order, not ID-list order.
    assert _headers(out) == ["read_1 length=4", "read_3 some description"]
    assert summary.matched == 2
    assert summary.missing == []


def test_extract_reports_missing_ids_instead_of_raising(tmp_path: Path) -> None:
    src = _fasta(tmp_path)
    summary = extract_records_by_ids(src, ["read_1", "absent_a", "absent_b"], tmp_path)

    assert summary.matched == 1
    assert summary.found_ids == 1
    assert summary.requested_ids == 3
    assert summary.missing == ["absent_a", "absent_b"]
    assert summary.missing_path is not None
    assert summary.missing_path.read_text(encoding="utf-8").split() == [
        "absent_a",
        "absent_b",
    ]


def test_extract_emits_every_record_sharing_a_duplicated_header(tmp_path: Path) -> None:
    src = _fasta(tmp_path)
    summary = extract_records_by_ids(src, ["read_2"], tmp_path)
    assert summary.matched == 2
    assert _headers(summary.written[0]) == ["read_2", "read_2"]


def test_extract_exclude_mode_keeps_everything_else(tmp_path: Path) -> None:
    src = _fasta(tmp_path)
    summary = extract_records_by_ids(src, ["read_2"], tmp_path, exclude=True)
    assert _headers(summary.written[0]) == ["read_1 length=4", "read_3 some description"]
    assert summary.matched == 2


def test_extract_one_file_per_record(tmp_path: Path) -> None:
    src = _fasta(tmp_path)
    out_dir = tmp_path / "out"
    summary = extract_records_by_ids(
        src, ["read_1", "read_3"], out_dir, single_file=False
    )
    names = sorted(p.name for p in summary.written)
    assert names == ["read_1.fasta", "read_3.fasta"]
    assert _headers(out_dir / "read_1.fasta") == ["read_1 length=4"]


def test_extract_reads_gzipped_input(tmp_path: Path) -> None:
    src = tmp_path / "reads.fasta.gz"
    with gzip.open(src, "wt", encoding="utf-8") as handle:
        handle.write(FASTA)
    summary = extract_records_by_ids(src, ["read_3"], tmp_path)
    assert summary.matched == 1
    assert _headers(summary.written[0]) == ["read_3 some description"]


def test_extract_sequence_content_is_preserved(tmp_path: Path) -> None:
    src = tmp_path / "wrapped.fasta"
    src.write_text(">long\n" + "ACGT\n" * 40 + ">short\nTT\n", encoding="utf-8")
    summary = extract_records_by_ids(src, ["long"], tmp_path)
    body = summary.written[0].read_text(encoding="utf-8").splitlines()[1:]
    assert "".join(body) == "ACGT" * 40


def test_extract_without_ids_is_an_error(tmp_path: Path) -> None:
    src = _fasta(tmp_path)
    with pytest.raises(SequenceComposeError):
        extract_records_by_ids(src, [], tmp_path)


def test_extract_cancelled_leaves_no_partial_output(tmp_path: Path) -> None:
    src = tmp_path / "many.fasta"
    src.write_text(
        "".join(f">read_{i}\nACGT\n" for i in range(6000)), encoding="utf-8"
    )
    summary = extract_records_by_ids(
        src, [f"read_{i}" for i in range(6000)], tmp_path, cancel=lambda: True
    )
    assert summary.cancelled is True
    assert summary.written == []
    assert not (tmp_path / "many.fasta_extr.fa").exists()


# --- the Check button ------------------------------------------------------


def test_count_matching_ids_reports_found_and_missing(tmp_path: Path) -> None:
    src = _fasta(tmp_path)
    summary = count_matching_ids(src, ["read_1", "read_2", "nope"])
    assert summary.found_ids == 2
    assert summary.requested_ids == 3
    assert summary.missing == ["nope"]
    assert summary.records_read == 4


def test_count_matching_ids_writes_nothing(tmp_path: Path) -> None:
    src = _fasta(tmp_path)
    before = sorted(p.name for p in tmp_path.iterdir())
    count_matching_ids(src, ["read_1"])
    assert sorted(p.name for p in tmp_path.iterdir()) == before
