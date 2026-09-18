"""Detect open reading frames with orfipy for sequence-map auto annotation."""

from __future__ import annotations

from dataclasses import dataclass

from repeatviewer.backends.auto_annotation_settings import OrfSettings
from repeatviewer.backends.sequence_features import SequenceFeature

try:
    import orfipy_core
except ImportError:  # pragma: no cover - optional at import time
    orfipy_core = None  # type: ignore[assignment]


@dataclass(frozen=True)
class OrfDetectResult:
    features: tuple[SequenceFeature, ...]
    message: str


def orfipy_available() -> bool:
    global orfipy_core
    if orfipy_core is not None:
        return True
    try:
        import orfipy_core as _orfipy_core
    except ImportError:
        return False
    orfipy_core = _orfipy_core
    return True


def _require_orfipy() -> None:
    if orfipy_core is None:
        raise RuntimeError(
            "orfipy is not installed. Install RepeatViewer with ORF detection support:\n"
            "  pip install orfipy"
        )


def parse_start_codons(settings: OrfSettings | None = None) -> list[str]:
    """Build orfipy start-codon list: ATG plus optional alternatives."""
    cfg = settings or OrfSettings()
    codons = ["ATG"]
    for token in cfg.alternative_start_codons.replace(",", " ").split():
        codon = token.strip().upper()
        if len(codon) == 3 and codon not in codons:
            codons.append(codon)
    return codons


def _orfipy_to_feature(
    start: int,
    stop: int,
    strand: str,
    *,
    index: int,
) -> SequenceFeature | None:
    """Convert orfipy 0-based half-open coords to a 1-based inclusive feature."""
    feat_start = int(start) + 1
    feat_end = int(stop)
    if feat_end < feat_start:
        return None
    feat_strand = "-" if str(strand).strip() == "-" else "+"
    return SequenceFeature(
        start=feat_start,
        end=feat_end,
        name=f"ORF{index}",
        feature_type="ORF",
        strand=feat_strand,
    )


def detect_orf_features(
    sequence_id: str,
    sequence: str,
    *,
    settings: OrfSettings | None = None,
) -> OrfDetectResult:
    """Run orfipy on a single sequence and return ORF interval features."""
    _require_orfipy()
    cfg = settings or OrfSettings()
    seq = str(sequence).upper()
    if not seq:
        return OrfDetectResult((), "Sequence is empty.")

    min_len = max(1, int(cfg.min_length))
    max_len = max(min_len, int(cfg.max_length))
    starts = parse_start_codons(cfg)

    features: list[SequenceFeature] = []
    index = 0
    for start, stop, strand, _description in orfipy_core.orfs(
        seq,
        name=sequence_id or "Seq",
        minlen=min_len,
        maxlen=max_len,
        strand="b",
        starts=starts,
    ):
        index += 1
        feat = _orfipy_to_feature(start, stop, strand, index=index)
        if feat is not None:
            features.append(feat)

    if not features:
        start_label = ", ".join(starts)
        return OrfDetectResult(
            (),
            f"No ORFs found (length {min_len}–{max_len} bp, starts: {start_label}).",
        )

    start_label = ", ".join(starts)
    return OrfDetectResult(
        tuple(features),
        f"Detected {len(features)} ORF(s) ({min_len}–{max_len} bp, starts: {start_label}).",
    )
