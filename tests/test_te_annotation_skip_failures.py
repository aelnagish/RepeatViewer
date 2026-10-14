"""Batch feature detection should skip failing sequences and finish the rest."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.auto_annotation import AutoAnnotationMode
from repeatviewer.backends.te_annotation_batch import (
    TeAnnotationBatchResult,
    TeAnnotationOptions,
    run_te_annotation_on_fasta,
)


def test_skipped_reminder_text() -> None:
    result = TeAnnotationBatchResult(
        features_by_id={},
        record_ids=["a", "b"],
        lengths_by_id={"a": 10, "b": 10},
        skipped_sequences=(
            ("seq_bad", "makeblastdb failed while building the self-BLAST database."),
            ("seq_empty", "Sequence is empty"),
        ),
    )
    headline, body = result.skipped_reminder_parts(kind="Feature detection")
    assert "2 sequence(s) were skipped" in headline
    assert "seq_bad" in body
    assert "makeblastdb" in body
    assert "Other sequences were annotated normally." in body
    assert "makeblastdb" not in headline
    text = result.skipped_reminder_text(kind="Feature detection")
    assert "2 sequence(s) were skipped" in text
    assert TeAnnotationBatchResult(
        features_by_id={},
        record_ids=[],
        lengths_by_id={},
    ).skipped_reminder_text() == ""


def test_batch_skips_failed_sequence_and_continues(tmp_path: Path, monkeypatch) -> None:
    fasta = tmp_path / "copies.fa"
    fasta.write_text(
        ">good\n" + ("ACGT" * 50) + "\n"
        ">bad\n" + ("ACGT" * 40) + "\n"
        ">also_good\n" + ("TTTT" * 50) + "\n",
        encoding="utf-8",
    )

    def _fake_annotate(record_fasta, options, **kwargs):
        del options, kwargs
        # Temp files are named seq_{index}_{id}.fa
        if "_bad" in Path(record_fasta).name or Path(record_fasta).name == "bad.fa":
            raise RuntimeError(
                "makeblastdb failed while building the self-BLAST database."
            )
        return []

    monkeypatch.setattr(
        "repeatviewer.backends.te_annotation_batch._annotate_one_record_fasta",
        _fake_annotate,
    )
    # Force sequential so progress ordering is irrelevant.
    monkeypatch.setattr(
        "repeatviewer.backends.te_annotation_batch._annotation_parallel_sequences",
        lambda _n: 1,
    )

    result = run_te_annotation_on_fasta(
        fasta,
        TeAnnotationOptions(auto_modes=(AutoAnnotationMode.LTR,)),
        persist=False,
    )
    assert set(result.record_ids) == {"good", "bad", "also_good"}
    assert len(result.skipped_sequences) == 1
    skipped_id, reason = result.skipped_sequences[0]
    assert skipped_id == "bad"
    assert "makeblastdb" in reason
    # Failed sequence keeps prior (empty) features; others still annotated.
    assert result.features_by_id["bad"] == []
    assert "bad" in result.features_by_id


def test_empty_sequence_skips_self_blast_without_makeblastdb(
    tmp_path: Path, monkeypatch
) -> None:
    from repeatviewer.backends.terminal_repeat_detect import auto_annotate_terminal_repeats

    empty = tmp_path / "empty.fa"
    empty.write_text(">empty\n\n", encoding="utf-8")

    def _boom(*_args, **_kwargs):
        raise AssertionError("run_self_blast should not be called for empty sequence")

    monkeypatch.setattr(
        "repeatviewer.backends.terminal_repeat_detect.run_self_blast",
        _boom,
    )
    result = auto_annotate_terminal_repeats(empty, repeat_mode="LTR")
    assert result.features == ()
    assert result.repeat_type is None
    assert "empty" in result.message.lower()
