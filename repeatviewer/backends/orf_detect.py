"""Detect open reading frames with orfipy for sequence-map auto annotation."""

from __future__ import annotations

import re
from dataclasses import dataclass

from repeatviewer.backends.auto_annotation_settings import (
    DEFAULT_ORF_FRAMES,
    OrfSettings,
    normalize_orf_frames,
)
from repeatviewer.backends.sequence_features import SequenceFeature

try:
    import orfipy_core
except ImportError:  # pragma: no cover - optional at import time
    orfipy_core = None  # type: ignore[assignment]

_ORF_FRAME_RE = re.compile(r"ORF_frame=(-?\d+)")


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


def _orfipy_strand_arg(frames: tuple[int, ...]) -> str:
    """Map selected frames to orfipy's strand flag (f / r / b)."""
    has_fwd = any(frame > 0 for frame in frames)
    has_rev = any(frame < 0 for frame in frames)
    if has_fwd and has_rev:
        return "b"
    if has_rev:
        return "r"
    return "f"


def _parse_orf_frame(description: str, start: int, strand: str) -> int:
    """Read orfipy's ORF_frame, with a start-position fallback."""
    match = _ORF_FRAME_RE.search(str(description or ""))
    if match:
        try:
            return int(match.group(1))
        except ValueError:
            pass
    # Fallback when description lacks ORF_frame (should be rare).
    offset = (int(start) % 3) + 1
    return -offset if str(strand).strip() == "-" else offset


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
    frames = normalize_orf_frames(cfg.frames)
    frame_set = set(frames)
    strand_arg = _orfipy_strand_arg(frames)

    features: list[SequenceFeature] = []
    index = 0
    for start, stop, strand, description in orfipy_core.orfs(
        seq,
        name=sequence_id or "Seq",
        minlen=min_len,
        maxlen=max_len,
        strand=strand_arg,
        starts=starts,
    ):
        frame = _parse_orf_frame(description, start, strand)
        if frame not in frame_set:
            continue
        index += 1
        feat = _orfipy_to_feature(start, stop, strand, index=index)
        if feat is not None:
            features.append(feat)

    frame_label = ",".join(str(frame) for frame in frames)
    start_label = ", ".join(starts)
    if not features:
        return OrfDetectResult(
            (),
            f"No ORFs found with orfipy (length {min_len}–{max_len} bp, "
            f"frames {frame_label}, starts: {start_label}).",
        )

    return OrfDetectResult(
        tuple(features),
        f"Detected {len(features)} ORF(s) with orfipy "
        f"({min_len}–{max_len} bp, frames {frame_label}, starts: {start_label}).",
    )


__all__ = [
    "DEFAULT_ORF_FRAMES",
    "OrfDetectResult",
    "detect_orf_features",
    "orfipy_available",
    "parse_start_codons",
]
