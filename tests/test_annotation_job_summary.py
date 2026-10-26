"""Process-monitor end-of-job feature count summaries."""

from __future__ import annotations

from dataclasses import dataclass, field

from repeatviewer.widgets.process_monitor_dialog import (
    annotation_job_summary_lines,
    format_annotation_run_summary,
)


@dataclass(frozen=True)
class _Feat:
    feature_type: str


@dataclass(frozen=True)
class AutoAnnotationResult:
    features: tuple
    mode: object
    message: str


@dataclass(frozen=True)
class MultiAutoAnnotationResult:
    features: tuple
    modes: tuple
    message: str
    found_by_type: dict = field(default_factory=dict)
    overlap_by_type: dict = field(default_factory=dict)
    new_by_type: dict = field(default_factory=dict)


@dataclass(frozen=True)
class SimilarityAnnotationResult:
    features: tuple
    message: str


@dataclass
class _Lifted:
    feature: object | None


@dataclass
class AnnotationTransferResult:
    lifted: list
    message: str = ""


@dataclass
class TeAnnotationBatchResult:
    features_by_id: dict
    record_ids: list
    identified_count: int = 0
    identified_by_type: dict | None = None
    annotated_sequence_count: int = 0


def test_format_summary_plain_sentences_and_per_type() -> None:
    lines = format_annotation_run_summary(
        n_sequences=15,
        counts_by_type={"LTR": 5, "PBS": 15},
        preferred_types=("LTR", "TIR", "PBS"),
    )
    assert lines[0] == "Finished annotating 15 sequences."
    assert lines[1] == "Found 20 new features (LTR 5, TIR 0, PBS 15)."


def test_format_summary_found_overlap_and_new() -> None:
    lines = format_annotation_run_summary(
        n_sequences=15,
        preferred_types=("ORF", "PBS"),
        found_by_type={"ORF": 49, "PBS": 0},
        overlap_by_type={"ORF": 49},
        counts_by_type={"ORF": 0, "PBS": 0},
        total_new=0,
    )
    assert lines == [
        "Finished annotating 15 sequences.",
        "Found 49 features (ORF 49, PBS 0): 49 overlap existing (ORF 49), 0 new.",
    ]


def test_format_summary_mixed_new_and_overlap() -> None:
    lines = format_annotation_run_summary(
        n_sequences=1,
        preferred_types=("ORF", "PBS"),
        found_by_type={"ORF": 10, "PBS": 3},
        overlap_by_type={"ORF": 7},
        counts_by_type={"ORF": 3, "PBS": 3},
        total_new=6,
    )
    assert lines[1] == (
        "Found 13 features (ORF 10, PBS 3): "
        "7 overlap existing (ORF 7), 6 new (ORF 3, PBS 3)."
    )


def test_format_summary_overlap_breakdown_by_kind() -> None:
    lines = format_annotation_run_summary(
        n_sequences=15,
        preferred_types=("LTR", "TIR", "tandem repeat", "ORF"),
        found_by_type={
            "LTR": 40,
            "TIR": 130,
            "tandem repeat": 233,
            "ORF": 52,
        },
        overlap_by_type={"LTR": 8, "ORF": 24},
        counts_by_type={
            "LTR": 32,
            "TIR": 130,
            "tandem repeat": 233,
            "ORF": 28,
        },
        total_new=423,
    )
    assert lines[1] == (
        "Found 455 features (LTR 40, TIR 130, tandem repeat 233, ORF 52): "
        "32 overlap existing (LTR 8, ORF 24), "
        "423 new (LTR 32, TIR 130, tandem repeat 233, ORF 28)."
    )


def test_format_summary_zero_features() -> None:
    lines = format_annotation_run_summary(
        n_sequences=15,
        counts_by_type={},
        preferred_types=("PBS",),
        total_new=0,
    )
    assert lines == [
        "Finished annotating 15 sequences.",
        "Found no new features.",
    ]


def test_format_summary_single_type() -> None:
    lines = format_annotation_run_summary(
        n_sequences=1,
        counts_by_type={"PBS": 1},
        preferred_types=("PBS",),
    )
    assert lines == [
        "Finished annotating 1 sequence.",
        "Found 1 new feature (PBS 1).",
    ]


def test_auto_annotation_summary_includes_count() -> None:
    result = AutoAnnotationResult(
        features=(_Feat("ORF"), _Feat("ORF"), _Feat("ORF")),
        mode="ORF",
        message="Detected 3 ORF(s) with orfipy.",
    )
    lines = annotation_job_summary_lines(result)
    assert "Detected 3 ORF(s) with orfipy." in lines
    assert "Finished annotating 1 sequence." in lines
    assert "Found 3 features (ORF 3): 0 overlap existing, 3 new (ORF 3)." in lines


def test_multi_auto_counts_only_selected_modes() -> None:
    class Mode:
        def __init__(self, value: str) -> None:
            self.value = value

    result = MultiAutoAnnotationResult(
        features=(
            _Feat("ORF"),
            _Feat("ORF"),
            _Feat("domain"),
            _Feat("LTR"),
        ),
        modes=(Mode("ORF"), Mode("LTR")),
        message="Detected 2 ORF(s).\nLTR: 1 pair.",
        found_by_type={"ORF": 2, "LTR": 1},
        overlap_by_type={"ORF": 2},
        new_by_type={"LTR": 1},
    )
    lines = annotation_job_summary_lines(result)
    assert "Finished annotating 1 sequence." in lines
    assert (
        "Found 3 features (ORF 2, LTR 1): "
        "2 overlap existing (ORF 2), 1 new (LTR 1)."
    ) in lines


def test_similarity_sums_hit_counts() -> None:
    result = SimilarityAnnotationResult(
        features=(),
        message=(
            "Pfam: 4 pyhmmer domain hit(s) (E-value ≤ 1e-5).\n"
            "Dfam: 2 BLASTN hit(s) (E-value ≤ 1e-5)."
        ),
    )
    lines = annotation_job_summary_lines(result)
    assert "Found 6 new features (similarity hit 6)." in lines


def test_transfer_counts_placed_features() -> None:
    result = AnnotationTransferResult(
        lifted=[_Lifted("f"), _Lifted(None), _Lifted("g")],
        message="placed 2 of 3.",
    )
    lines = annotation_job_summary_lines(result)
    # Message is streamed live during transfer; summary is the count only.
    assert lines == [
        "Finished annotating 1 sequence.",
        "Found 2 new features (transferred 2).",
    ]


def test_transfer_tuple_payload() -> None:
    result = (
        AnnotationTransferResult(lifted=[_Lifted("f")], message="ok"),
        "tid",
        "seq",
        None,
    )
    lines = annotation_job_summary_lines(result)
    assert "Found 1 new feature (transferred 1)." in lines


def test_batch_summary_already_streamed() -> None:
    result = TeAnnotationBatchResult(
        features_by_id={"a": [1, 2], "b": [3]},
        record_ids=["a", "b"],
        identified_count=5,
        annotated_sequence_count=2,
    )
    # Batch emits the plain summary during the job; do not repeat it here.
    assert annotation_job_summary_lines(result) == []
