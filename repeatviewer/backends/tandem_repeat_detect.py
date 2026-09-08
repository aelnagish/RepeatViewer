"""Detect tandem repeats with pytrf (STR + GTR + ATR) for sequence-map auto annotation."""

from __future__ import annotations

from dataclasses import dataclass

from repeatviewer.backends.sequence_features import SequenceFeature
from repeatviewer.backends.auto_annotation_settings import TandemRepeatSettings

try:
    import pytrf
except ImportError:  # pragma: no cover - optional at import time
    pytrf = None  # type: ignore[assignment]


@dataclass(frozen=True)
class TandemRepeatResult:
    features: tuple[SequenceFeature, ...]
    message: str


def _require_pytrf():
    if pytrf is None:
        raise RuntimeError(
            "pytrf is not installed. Install RepeatViewer with tandem-repeat support:\n"
            "  pip install pytrf"
        )


def _intervals_overlap(a_start: int, a_end: int, b_start: int, b_end: int) -> bool:
    return a_start <= b_end and b_start <= a_end


def _locus_span(start: int, end: int) -> int:
    return int(end) - int(start) + 1


def _overlaps_accepted_loci(
    start: int,
    end: int,
    accepted: list[tuple[int, int]],
) -> bool:
    """True when a prior locus overlaps and has the same total span (bp length)."""
    span = _locus_span(start, end)
    for lo, hi in accepted:
        if _locus_span(lo, hi) != span:
            continue
        if _intervals_overlap(start, end, lo, hi):
            return True
    return False


def _atr_min_identity(value: float) -> float:
    """Convert UI fraction (0–1) to pytrf percent scale (0–100)."""
    value = float(value)
    if value <= 1.0:
        return value * 100.0
    return value


def _parse_pytrf_locus(row: tuple, *, prefix: str, locus_index: int) -> list[SequenceFeature]:
    """
    Expand one pytrf locus tuple into per-unit arrow features.

    Tuple layout from pytrf ``as_list()``:
    (name, start, end, motif, motif_len, repeat_count, total_length)
    """
    if len(row) < 6:
        return []

    start = int(row[1])
    end = int(row[2])
    motif = str(row[3])
    unit_len = int(row[4]) if int(row[4]) > 0 else max(1, len(motif))
    repeat_count = int(float(row[5]))
    if repeat_count <= 0 or unit_len <= 0:
        return []

    expected_end = start + repeat_count * unit_len - 1
    if expected_end > end:
        repeat_count = max(1, (end - start + 1) // unit_len)

    features: list[SequenceFeature] = []
    locus_label = f"{prefix}{locus_index}"
    for unit_idx in range(repeat_count):
        unit_start = start + unit_idx * unit_len
        unit_end = unit_start + unit_len - 1
        if unit_end > end:
            break
        features.append(
            SequenceFeature(
                start=unit_start,
                end=unit_end,
                name=f"{locus_label}_u{unit_idx + 1}",
                feature_type="tandem repeat",
                strand="+",
            )
        )
    return features


def _append_non_overlapping_loci(
    rows: list[tuple],
    *,
    prefix: str,
    accepted_intervals: list[tuple[int, int]],
    locus_counter: int,
) -> tuple[list[SequenceFeature], int, int]:
    """Expand pytrf rows into unit features, skipping duplicate overlapping loci."""
    features: list[SequenceFeature] = []
    locus_count = 0
    for row in rows:
        if len(row) < 4:
            continue
        start = int(row[1])
        end = int(row[2])
        if _overlaps_accepted_loci(start, end, accepted_intervals):
            continue
        accepted_intervals.append((start, end))
        locus_counter += 1
        locus_count += 1
        features.extend(_parse_pytrf_locus(row, prefix=prefix, locus_index=locus_counter))
    return features, locus_counter, locus_count


def _run_atr_finder(name: str, seq: str, cfg: TandemRepeatSettings):
    return pytrf.ATRFinder(
        name,
        seq,
        min_motif=max(1, int(cfg.atr_min_motif)),
        max_motif=max(1, int(cfg.atr_max_motif)),
        min_seedrep=max(1, int(cfg.atr_min_seed_repeat)),
        min_seedlen=max(1, int(cfg.atr_min_seed_length)),
        min_identity=_atr_min_identity(cfg.atr_min_identity),
        max_errors=max(0, int(cfg.atr_max_errors)),
        max_extend=max(1, int(cfg.atr_max_extend)),
    )


def detect_tandem_repeat_features(
    sequence_id: str,
    sequence: str,
    *,
    settings: TandemRepeatSettings | None = None,
) -> TandemRepeatResult:
    """Run pytrf STR/GTR and optional ATR; return one arrow per repeat unit."""
    _require_pytrf()
    cfg = settings or TandemRepeatSettings()

    name = sequence_id or "sequence"
    seq = str(sequence).upper()
    if not seq:
        return TandemRepeatResult(features=(), message="Sequence is empty.")

    ssrs = pytrf.STRFinder(
        name,
        seq,
        100,
        cfg.min_repeat_2,
        cfg.min_repeat_3,
        cfg.min_repeat_4,
        cfg.min_repeat_5,
        cfg.min_repeat_6,
    )
    gtrs = pytrf.GTRFinder(
        name,
        seq,
        min_motif=7,  # motifs 1–6 bp are handled by STRFinder
        max_motif=100,
        min_repeat=cfg.min_repeat_6,
        min_length=10,
    )

    features: list[SequenceFeature] = []
    accepted_intervals: list[tuple[int, int]] = []
    locus_counter = 0
    ssr_n = 0
    gtr_n = 0
    atr_n = 0

    # ATR runs first so long imperfect blocks are not discarded by exact STR/GTR hits.
    if cfg.enable_approximate:
        atrs = _run_atr_finder(name, seq, cfg)
        added, locus_counter, atr_n = _append_non_overlapping_loci(
            atrs.as_list(),
            prefix="ATR",
            accepted_intervals=accepted_intervals,
            locus_counter=locus_counter,
        )
        features.extend(added)

    added, locus_counter, ssr_n = _append_non_overlapping_loci(
        ssrs.as_list(),
        prefix="SSR",
        accepted_intervals=accepted_intervals,
        locus_counter=locus_counter,
    )
    features.extend(added)

    added, locus_counter, gtr_n = _append_non_overlapping_loci(
        gtrs.as_list(),
        prefix="GTR",
        accepted_intervals=accepted_intervals,
        locus_counter=locus_counter,
    )
    features.extend(added)

    if not features:
        mode = "pytrf STR/GTR"
        if cfg.enable_approximate:
            mode += "/ATR"
        return TandemRepeatResult(
            features=(),
            message=f"No tandem repeats were found ({mode}).",
        )

    parts: list[str] = []
    if cfg.enable_approximate:
        parts.append(f"{atr_n} ATR")
    parts.extend([f"{ssr_n} STR", f"{gtr_n} GTR"])
    return TandemRepeatResult(
        features=tuple(features),
        message=(
            f"Detected {len(features)} tandem repeat unit(s) "
            f"({', '.join(parts)} locus/loci; same-length overlaps skipped)."
        ),
    )


def merge_tandem_repeat_features(
    existing: list[SequenceFeature],
    detected: list[SequenceFeature] | tuple[SequenceFeature, ...],
) -> list[SequenceFeature]:
    """Replace existing tandem-repeat features with newly detected units."""
    kept = [f for f in existing if f.feature_type != "tandem repeat"]
    return kept + list(detected)
