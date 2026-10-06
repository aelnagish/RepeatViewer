"""Logical progress phases for long-running analysis tools.

The status bar reserves one extra step for “Opening results…”, so a single
work unit used to show ``0 / 2`` for the entire BLAST run. These phase counts
split each file (or file×genome) into stages that match the real pipeline.
"""

from __future__ import annotations

from pathlib import Path
from typing import Callable

# Local BLAST: prepare → run → write outputs (single-sequence files).
BLAST_PHASES = 3
# Online BLAST: connect/submit → wait → download → HTML → table.
ONLINE_BLAST_PHASES = 5
# Align (MAFFT / MUSCLE): prepare → align → finish.
ALIGN_PHASES = 3
# Flexidot: prepare → run → write plots.
FLEXIDOT_PHASES = 3
# TE-Aid: consensus → genome BLAST → self-dot → assemble → write HTML.
TEAID_PHASES = 5
# Mapping: start → read/index → map → collect/consensus → write.
MAP_PHASES = 5

ProgressFn = Callable[[int, int, str], None]


def emit_phase(
    progress: ProgressFn | None,
    step: int,
    total: int,
    detail: str,
) -> None:
    """Report ``step`` of ``total`` (``step`` may be 0 … total)."""
    if progress is None:
        return
    progress(max(0, int(step)), max(1, int(total)), str(detail))


def blast_units_for_fasta(
    seq_path: Path,
    *,
    online: bool = False,
    sequence_count: int | None = None,
) -> int:
    """Progress units for one BLAST query file."""
    if online:
        return ONLINE_BLAST_PHASES
    if sequence_count is None:
        from repeatviewer.backends.curation import _fasta_sequence_units

        sequence_count = _fasta_sequence_units(seq_path)
    count = max(1, int(sequence_count))
    # Multi-sequence files advance once per record (existing BLAST logic).
    if count > 1:
        return count
    return BLAST_PHASES


def tool_units_for_fasta(func_name: str, seq_path: Path, *, online_blast: bool = False) -> int:
    """Progress units for one work item of ``func_name``."""
    if func_name == "blast":
        return blast_units_for_fasta(seq_path, online=online_blast)
    if func_name == "flexidot":
        return FLEXIDOT_PHASES
    if func_name == "teaid":
        return TEAID_PHASES
    if func_name == "te_annotation":
        from repeatviewer.backends.curation import _fasta_sequence_units

        return _fasta_sequence_units(seq_path)
    return 1
