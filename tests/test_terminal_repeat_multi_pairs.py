"""LTR/TIR detection can keep multiple nested pairs (outermost first)."""

from __future__ import annotations

import pandas as pd

from repeatviewer.backends.auto_annotation_settings import TerminalRepeatSettings
from repeatviewer.backends.terminal_repeat_detect import (
    detect_terminal_repeat_features,
)


def _ltr_row(
    *,
    qstart: int,
    qend: int,
    sstart: int,
    send: int,
    pident: float = 95.0,
) -> dict:
    return {
        "qseqid": "seq",
        "sseqid": "seq",
        "pident": pident,
        "length": qend - qstart + 1,
        "mismatch": 0,
        "gapopen": 0,
        "qstart": qstart,
        "qend": qend,
        "sstart": sstart,
        "send": send,
        "sstrand": "plus",
        "evalue": 0.0,
        "bitscore": 100.0,
    }


def test_detects_up_to_three_ltr_pairs_outermost_first() -> None:
    # Nested LTRs on a 2000 bp sequence: outer, middle, inner.
    record_len = 2000
    df = pd.DataFrame(
        [
            _ltr_row(qstart=10, qend=110, sstart=1890, send=1990),  # outer
            _ltr_row(qstart=200, qend=300, sstart=1700, send=1800),  # middle
            _ltr_row(qstart=400, qend=500, sstart=1500, send=1600),  # inner
            _ltr_row(qstart=11, qend=111, sstart=1891, send=1991),  # near-dup of outer
        ]
    )
    result = detect_terminal_repeat_features(
        df,
        record_len,
        repeat_mode="LTR",
        ltr_settings=TerminalRepeatSettings(
            min_length=50,
            max_length=2000,
            min_pair_distance=50,
            max_distance_to_end=500,
            similarity_threshold=0.5,
            max_pairs=3,
        ),
    )
    assert result.repeat_type == "LTR"
    names = [feat.name for feat in result.features]
    assert names[0].startswith("LTR5")
    assert names[1].startswith("LTR3")
    assert any(n.startswith("LTR5_2") for n in names)
    assert any(n.startswith("LTR3_2") for n in names)
    assert any(n.startswith("LTR5_3") for n in names)
    assert any(n.startswith("LTR3_3") for n in names)
    # Near-duplicate outer should not create a 4th pair.
    assert sum(1 for n in names if n.startswith("LTR5")) == 3


def test_max_pairs_one_keeps_only_outermost() -> None:
    record_len = 2000
    df = pd.DataFrame(
        [
            _ltr_row(qstart=10, qend=110, sstart=1890, send=1990),
            _ltr_row(qstart=200, qend=300, sstart=1700, send=1800),
        ]
    )
    result = detect_terminal_repeat_features(
        df,
        record_len,
        repeat_mode="LTR",
        ltr_settings=TerminalRepeatSettings(
            min_length=50,
            max_length=2000,
            min_pair_distance=50,
            max_distance_to_end=500,
            similarity_threshold=0.5,
            max_pairs=1,
        ),
    )
    names = [feat.name for feat in result.features]
    assert len(result.features) == 2
    assert names[0].startswith("LTR5")
    assert names[1].startswith("LTR3")
    assert not any("_2" in n for n in names)
