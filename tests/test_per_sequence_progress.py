"""Per-sequence progress units for BLAST and Sequence Annotation."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.auto_annotation import AutoAnnotationMode
from repeatviewer.backends.curation import _PER_SEQUENCE_FUNCS, _fasta_sequence_units
from repeatviewer.backends.te_annotation_batch import (
    TeAnnotationOptions,
    run_te_annotation_on_fasta,
)
from repeatviewer.widgets.task_progress import TaskProgressPanel


def test_fasta_sequence_units_counts_records(tmp_path: Path) -> None:
    fasta = tmp_path / "multi.fa"
    fasta.write_text(
        ">a\nACGT\n>b\nACGT\n>c\nACGT\n",
        encoding="utf-8",
    )
    assert _fasta_sequence_units(fasta) == 3
    assert "te_annotation" in _PER_SEQUENCE_FUNCS
    assert "blast" not in _PER_SEQUENCE_FUNCS


def test_task_progress_panel_shows_fraction(qapp) -> None:
    panel = TaskProgressPanel()
    panel.start("Running BLAST…", 5)
    assert panel._bar.maximum() == 5
    panel.update(2, "lib.fa: sequence 2/4 (TE2)")
    assert panel._bar.value() == 2
    assert panel._bar.format() == "%v / %m"
    panel.deleteLater()


def test_te_annotation_emits_per_sequence_progress(tmp_path: Path, monkeypatch) -> None:
    fasta = tmp_path / "lib.fa"
    fasta.write_text(
        ">TE1\n" + ("A" * 60) + "\n"
        ">TE2\n" + ("C" * 60) + "\n",
        encoding="utf-8",
    )
    events: list[tuple[int, int, str]] = []

    monkeypatch.setattr(
        "repeatviewer.backends.te_annotation_batch._annotate_one_record_fasta",
        lambda *_args, **_kwargs: ([], {}, {}, {}),
    )
    monkeypatch.setattr(
        "repeatviewer.backends.te_annotation_batch.save_annotations_multi",
        lambda path, **kwargs: (tmp_path / "out.seqmap_RV.json", None),
    )
    monkeypatch.setattr(
        "repeatviewer.backends.te_annotation_batch.load_auto_annotation_settings",
        lambda: object(),
    )
    monkeypatch.setattr(
        "repeatviewer.backends.te_annotation_batch.load_similarity_annotation_settings",
        lambda: object(),
    )
    monkeypatch.setattr(
        "repeatviewer.backends.te_annotation_batch.load_use_system_blast_from_settings",
        lambda: False,
    )

    def progress(completed: int, total: int, label: str) -> None:
        events.append((completed, total, label))

    options = TeAnnotationOptions(auto_modes=(AutoAnnotationMode.ORF,))
    run_te_annotation_on_fasta(fasta, options, progress=progress)

    steps = [e for e in events if e[0] >= 0]
    totals = {t for _c, t, _label in steps}
    assert totals == {2}
    assert {c for c, _t, _label in steps} >= {0, 1, 2}
