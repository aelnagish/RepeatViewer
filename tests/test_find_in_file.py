"""Edit > Find over a whole FASTA file: names, annotations, motif, order, limits, Stop."""

from __future__ import annotations

import gzip
from dataclasses import dataclass

import pytest

from repeatviewer.backends.find_in_file import (
    FindRequest,
    find_in_fasta,
    first_index_after,
    iter_fasta_records,
    next_index,
)
from repeatviewer.backends.motif_search import FindCancelled, FindOptions

FASTA = (
    ">CroSat1 monomer one\nACGTAC\nGTAAAA\n\n"
    ">CroSat2\r\nTTTTACGTTT\r\n"
    ">CroSat1 duplicate header\nGGGGGG\n"
)


@pytest.fixture
def fasta(tmp_path):
    path = tmp_path / "monomers.fasta"
    path.write_text(FASTA)
    return path


def test_record_names_match_the_annotation_views(fasta):
    from repeatviewer.backends.sequence_features import fasta_record_lengths

    records = list(iter_fasta_records(fasta))
    assert [(rid, len(seq)) for rid, seq in records] == list(fasta_record_lengths(fasta).items())
    assert [rid for rid, _ in records] == ["CroSat1", "CroSat2", "CroSat1#2"]
    assert records[0][1] == "ACGTACGTAAAA"


def test_gzipped_fasta(tmp_path):
    path = tmp_path / "x.fa.gz"
    with gzip.open(path, "wt") as handle:
        handle.write(">a\nACGT\n")
    assert list(iter_fasta_records(path)) == [("a", "ACGT")]


def test_motif_hits_in_file_order_with_record_ids(fasta):
    request = FindRequest("ACGT", FindOptions(both_strands=False), in_names=False, in_annotations=False)
    result = find_in_fasta(fasta, request)
    assert [(i.seq_id, i.start, i.end) for i in result.items] == [
        ("CroSat1", 1, 4), ("CroSat1", 5, 8), ("CroSat2", 5, 8)
    ]
    assert result.record_ids == ["CroSat1", "CroSat2", "CroSat1#2"]


def test_sequence_names_and_case(fasta):
    request = FindRequest("crosat1", in_sequence=False, in_annotations=False)
    assert [i.seq_id for i in find_in_fasta(fasta, request).items] == ["CroSat1", "CroSat1#2"]
    strict = FindRequest("crosat1", in_sequence=False, in_annotations=False, case_sensitive=True)
    assert find_in_fasta(fasta, strict).items == []


@dataclass
class Feat:
    start: int
    end: int
    name: str
    feature_type: str = "others"
    strand: str = "+"


def test_annotation_names_from_the_open_view(fasta):
    features = {"CroSat2": [Feat(2, 7, "LTR_copia", "LTR", "-")], "CroSat1": [Feat(1, 3, "probe")]}
    request = FindRequest("copia", in_sequence=False, in_names=False)
    items = find_in_fasta(fasta, request, features_by_id=features).items
    assert [(i.seq_id, i.start, i.end, i.strand, i.kind, i.label) for i in items] == [
        ("CroSat2", 2, 7, "-", "annotation", "LTR_copia")
    ]
    by_type = find_in_fasta(fasta, FindRequest("ltr", in_sequence=False, in_names=False), features_by_id=features)
    assert len(by_type.items) == 1


def test_words_that_are_not_dna_still_search_names(fasta):
    result = find_in_fasta(fasta, FindRequest("CroSat2"))
    assert result.sequence_skipped and [i.kind for i in result.items] == ["name"]


def test_limit_progress_and_stop(fasta):
    ticks = []
    capped = find_in_fasta(
        fasta,
        FindRequest("A", FindOptions(both_strands=False, max_hits=3), in_names=False, in_annotations=False),
        progress=lambda done, rid: ticks.append((done, rid)),
    )
    assert len(capped.items) == 3 and capped.truncated
    assert ticks[0] == (1, "CroSat1")
    with pytest.raises(FindCancelled):
        find_in_fasta(fasta, FindRequest("ACGT"), cancelled=lambda: True)


def test_next_previous_and_first_after():
    from repeatviewer.backends.find_in_file import FoundItem

    items = [FoundItem("a", 5, 8), FoundItem("a", 20, 23), FoundItem("b", 1, 4)]
    assert next_index(items, -1, +1) == 0 and next_index(items, -1, -1) == 2
    assert next_index(items, 2, +1) == 0 and next_index(items, 0, -1) == 2
    assert first_index_after(items, "a", 10, ["a", "b"]) == 1
    assert first_index_after(items, "b", 2, ["a", "b"]) == 0  # wraps to the first hit
    assert next_index([], 0, 1) == -1
