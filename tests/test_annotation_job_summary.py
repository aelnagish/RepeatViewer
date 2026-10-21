"""Process-monitor end-of-job feature count summaries."""

from __future__ import annotations

from dataclasses import dataclass

from repeatviewer.widgets.process_monitor_dialog import annotation_job_summary_lines


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


def test_auto_annotation_summary_includes_count() -> None:
    result = AutoAnnotationResult(
        features=(_Feat("ORF"), _Feat("ORF"), _Feat("ORF")),
        mode="ORF",
        message="Detected 3 ORF(s) with orfipy.",
    )
    lines = annotation_job_summary_lines(result)
    assert "Detected 3 ORF(s) with orfipy." in lines
    assert "Identified features: 3" in lines


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
    )
    lines = annotation_job_summary_lines(result)
    assert "Identified features: 3" in lines


def test_similarity_sums_hit_counts() -> None:
    result = SimilarityAnnotationResult(
        features=(),
        message=(
            "Pfam: 4 pyhmmer domain hit(s) (E-value ≤ 1e-5).\n"
            "Dfam: 2 BLASTN hit(s) (E-value ≤ 1e-5)."
        ),
    )
    lines = annotation_job_summary_lines(result)
    assert "Identified features: 6" in lines


def test_transfer_counts_placed_features() -> None:
    result = AnnotationTransferResult(
        lifted=[_Lifted("f"), _Lifted(None), _Lifted("g")],
        message="placed 2 of 3.",
    )
    lines = annotation_job_summary_lines(result)
    # Message is streamed live during transfer; summary is the count only.
    assert lines == ["Identified features: 2"]


def test_transfer_tuple_payload() -> None:
    result = (
        AnnotationTransferResult(lifted=[_Lifted("f")], message="ok"),
        "tid",
        "seq",
        None,
    )
    assert "Identified features: 1" in annotation_job_summary_lines(result)


def test_batch_adds_sequence_count_only() -> None:
    result = TeAnnotationBatchResult(
        features_by_id={"a": [1, 2], "b": [3]},
        record_ids=["a", "b"],
        identified_count=5,
    )
    lines = annotation_job_summary_lines(result)
    assert lines == ["Annotated sequences: 2"]
