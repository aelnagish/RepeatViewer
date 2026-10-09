"""Map features between gapped MSA columns and ungapped sequence coordinates.

Detectors / BLAST / transfer should run on gap-free DNA/protein. Features are
then projected back onto the alignment (``-`` / ``.`` columns) for MSA display.
"""

from __future__ import annotations

from dataclasses import replace

from repeatviewer.backends.sequence_features import SequenceFeature

# Same gap alphabet as motif Find ("Allow gaps in sequence").
GAP_CHARS = frozenset("-.")


def sequence_has_gaps(sequence: str) -> bool:
    return any(ch in GAP_CHARS for ch in sequence)


def ungap_sequence(sequence: str) -> tuple[str, list[int] | None]:
    """Return ungapped text and 0-based gapped index of each kept residue.

    When there are no gaps, returns ``(sequence, None)``.
    """
    if not sequence_has_gaps(sequence):
        return sequence, None
    keep = [i for i, ch in enumerate(sequence) if ch not in GAP_CHARS]
    return "".join(sequence[i] for i in keep), keep


def ungapped_to_gapped_1based(ungapped_pos: int, keep: list[int]) -> int:
    """Map a 1-based ungapped position to a 1-based gapped (alignment) position."""
    if not keep:
        return max(1, int(ungapped_pos))
    idx = max(0, min(len(keep) - 1, int(ungapped_pos) - 1))
    return int(keep[idx]) + 1


def gapped_to_ungapped_1based(gapped_pos: int, keep: list[int]) -> int | None:
    """Map a 1-based gapped position onto the nearest kept ungapped residue.

    Returns ``None`` when ``keep`` is empty.
    """
    if not keep:
        return None
    target = max(0, int(gapped_pos) - 1)
    # First kept column at or after target; else last kept before.
    for i, g0 in enumerate(keep):
        if g0 >= target:
            return i + 1
    return len(keep)


def project_feature_to_gapped(
    feat: SequenceFeature, keep: list[int] | None
) -> SequenceFeature:
    """Project an ungapped-coordinate feature onto MSA columns."""
    if keep is None:
        return feat
    start = ungapped_to_gapped_1based(feat.start, keep)
    end = ungapped_to_gapped_1based(feat.end, keep)
    if end < start:
        start, end = end, start
    return replace(feat, start=start, end=end)


def project_feature_to_ungapped(
    feat: SequenceFeature, keep: list[int] | None
) -> SequenceFeature | None:
    """Project a gapped-coordinate feature onto the ungapped sequence."""
    if keep is None:
        return feat
    start = gapped_to_ungapped_1based(feat.start, keep)
    end = gapped_to_ungapped_1based(feat.end, keep)
    if start is None or end is None:
        return None
    if end < start:
        start, end = end, start
    return replace(feat, start=int(start), end=int(end))


def project_features_to_gapped(
    features: list[SequenceFeature], keep: list[int] | None
) -> list[SequenceFeature]:
    if keep is None:
        return list(features)
    return [project_feature_to_gapped(f, keep) for f in features]


def project_features_to_ungapped(
    features: list[SequenceFeature], keep: list[int] | None
) -> list[SequenceFeature]:
    if keep is None:
        return list(features)
    out: list[SequenceFeature] = []
    for feat in features:
        mapped = project_feature_to_ungapped(feat, keep)
        if mapped is not None:
            out.append(mapped)
    return out
