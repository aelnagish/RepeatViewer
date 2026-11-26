"""Detect LTR and TIR terminal repeats via self-BLAST (ported from TEtrimmer)."""

from __future__ import annotations

import logging
import os
import platform
import shutil
import tempfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd
from Bio import SeqIO

from repeatviewer.backends.gui_functions import (
    check_cmd_in_path,
    get_original_file_path,
)
from repeatviewer.backends.job_cancel import JobCancelledError, JobCancellation, check_cancelled, run_subprocess_with_cancel
from repeatviewer.backends.auto_annotation_settings import (
    DEFAULT_LTR_MIN_LENGTH,
    DEFAULT_LTR_SIMILARITY_THRESHOLD,
    DEFAULT_TIR_MIN_LENGTH,
    DEFAULT_TIR_SIMILARITY_THRESHOLD,
    TerminalRepeatSettings,
    normalize_max_terminal_pairs,
)
from repeatviewer.backends.sequence_features import SequenceFeature

MIN_LTR_LEN = DEFAULT_LTR_MIN_LENGTH
MIN_TIR_LEN = DEFAULT_TIR_MIN_LENGTH
DEFAULT_TERMINAL_REPEAT_WORD_SIZE = 11
MIN_TERMINAL_REPEAT_WORD_SIZE = 4
# Collapse near-identical BLAST rows when selecting multiple pairs.
_PAIR_DEDUP_TOLERANCE_BP = 25


def word_size_for_terminal_repeat(min_repeat_length: int) -> int:
    """
    Self-BLAST ``-word_size`` for LTR or TIR detection.

    Default 11. When the module's minimum repeat length is below 11 bp, use that
    value instead (floored at 4). LTR and TIR settings are independent.
    """
    min_repeat_length = int(min_repeat_length)
    if min_repeat_length < DEFAULT_TERMINAL_REPEAT_WORD_SIZE:
        return max(MIN_TERMINAL_REPEAT_WORD_SIZE, min_repeat_length)
    return DEFAULT_TERMINAL_REPEAT_WORD_SIZE


@dataclass(frozen=True)
class TerminalRepeatResult:
    """Outcome of terminal-repeat auto annotation."""

    features: tuple[SequenceFeature, ...]
    repeat_type: str | None  # "LTR", "TIR", or None
    message: str


def _record_length(fasta_path: Path) -> int:
    record = SeqIO.read(str(fasta_path), "fasta")
    return len(record.seq)


def _load_self_blast_table(blast_table: Path) -> pd.DataFrame:
    df = pd.read_csv(blast_table, sep="\t", header=0)
    required = ("qstart", "qend", "sstart", "send")
    missing = [col for col in required if col not in df.columns]
    if missing:
        raise ValueError(f"Self-BLAST table missing columns: {', '.join(missing)}")
    for col in required:
        df[col] = pd.to_numeric(df[col], errors="coerce")
    df = df.dropna(subset=list(required))
    if df.empty:
        return df
    return df.astype({col: int for col in required})


def _subject_inverted(df: pd.DataFrame) -> pd.Series:
    """True when the subject alignment is on the minus strand."""
    inverted = df["sstart"] > df["send"]
    if "sstrand" in df.columns:
        strand = df["sstrand"].astype(str).str.strip().str.lower()
        inverted = inverted | strand.isin({"minus", "-", "c", "reverse", "rev"})
    return inverted


def _prepare_blast_hits(df: pd.DataFrame) -> pd.DataFrame:
    """
    Normalize BLAST coordinates and derive orientation flags.

    BLAST+ often reports sstart < send together with sstrand=minus for inverted
    terminal-repeat hits. TEtrimmer's original filter used sstart > send only.
    """
    hits = df.copy()
    hits["qstart"], hits["qend"] = (
        hits[["qstart", "qend"]].min(axis=1),
        hits[["qstart", "qend"]].max(axis=1),
    )
    hits["s_lo"] = hits[["sstart", "send"]].min(axis=1)
    hits["s_hi"] = hits[["sstart", "send"]].max(axis=1)
    hits["s_inverted"] = _subject_inverted(hits)
    hits["s_forward"] = ~hits["s_inverted"]
    return hits


def _filter_blast_hits(
    hits: pd.DataFrame,
    *,
    min_len: int,
    max_len: int,
    similarity_threshold: float,
) -> pd.DataFrame:
    filtered = hits.copy()
    span = filtered["qend"] - filtered["qstart"]
    filtered = filtered[(span >= min_len) & (span <= max_len)]
    if "pident" in filtered.columns and similarity_threshold > 0:
        min_pident = similarity_threshold * 100.0
        filtered = filtered[pd.to_numeric(filtered["pident"], errors="coerce") >= min_pident]
    return filtered


def _copies_non_overlapping(hit: pd.Series) -> bool:
    """Left copy ``[qstart, qend]`` must end before right copy ``[s_lo, s_hi]``."""
    return int(hit["qend"]) < int(hit["s_lo"])


def _pair_distance_ok(hit: pd.Series, *, min_pair_distance: int) -> bool:
    # Always forbid overlapping terminal copies; optional min gap on top.
    if not _copies_non_overlapping(hit):
        return False
    if min_pair_distance <= 0:
        return True
    q3 = int(hit["qend"])
    s5 = int(hit["s_lo"])
    return s5 - q3 >= min_pair_distance


def _terminal_ends_ok(hit: pd.Series, record_len: int, *, max_distance_to_end: int) -> bool:
    """
    Require the 5′ copy near the sequence start and the 3′ copy near the end.

    ``max_distance_to_end`` is the maximum allowed gap (bp) from each terminus:
    left copy start within that many bp of position 1, right copy end within
    that many bp of ``record_len``.
    """
    allowed = max(0, int(max_distance_to_end))
    left_gap = max(0, int(hit["qstart"]) - 1)
    right_gap = max(0, record_len - int(hit["s_hi"]))
    return left_gap <= allowed and right_gap <= allowed


def _pair_coords(hit: pd.Series) -> tuple[int, int, int, int]:
    return (
        int(hit["qstart"]),
        int(hit["qend"]),
        int(hit["s_lo"]),
        int(hit["s_hi"]),
    )


def _pairs_nearly_identical(
    a: pd.Series,
    b: pd.Series,
    *,
    tol: int = _PAIR_DEDUP_TOLERANCE_BP,
) -> bool:
    return all(
        abs(x - y) <= tol for x, y in zip(_pair_coords(a), _pair_coords(b), strict=True)
    )


def _select_terminal_pairs(
    candidates: pd.DataFrame,
    record_len: int,
    *,
    cfg: TerminalRepeatSettings,
    score_column: str,
) -> list[pd.Series]:
    """Keep up to ``cfg.max_pairs`` outermost pairs that pass distance filters."""
    if candidates.empty:
        return []
    ranked = candidates.copy()
    ranked["score"] = ranked[score_column]
    ranked = ranked.sort_values("score", ascending=False)
    max_pairs = normalize_max_terminal_pairs(cfg.max_pairs)
    selected: list[pd.Series] = []
    for _, hit in ranked.iterrows():
        if not _terminal_ends_ok(
            hit,
            record_len,
            max_distance_to_end=cfg.max_distance_to_end,
        ):
            continue
        if not _pair_distance_ok(hit, min_pair_distance=cfg.min_pair_distance):
            continue
        if any(_pairs_nearly_identical(hit, kept) for kept in selected):
            continue
        selected.append(hit)
        if len(selected) >= max_pairs:
            break
    return selected


def _check_terminal_repeat(
    df: pd.DataFrame,
    record_len: int,
    *,
    repeat_mode: str = "both",
    ltr_settings: TerminalRepeatSettings | None = None,
    tir_settings: TerminalRepeatSettings | None = None,
) -> tuple[str | None, list[pd.Series]]:
    """
    Port of TEtrimmer ``check_terminal_repeat`` hit selection.

    Returns ``(repeat_type, hits)`` where ``hits`` are outermost-first pairs
    (up to ``max_pairs``). Unlike TEtrimmer, TIR is still evaluated when
    LTR-like hits exist but fail the terminal-end proximity checks.
    """
    if df.empty:
        return None, []

    ltr_cfg = ltr_settings or TerminalRepeatSettings(
        min_length=MIN_LTR_LEN,
        similarity_threshold=DEFAULT_LTR_SIMILARITY_THRESHOLD,
    )
    tir_cfg = tir_settings or TerminalRepeatSettings(
        min_length=MIN_TIR_LEN,
        similarity_threshold=DEFAULT_TIR_SIMILARITY_THRESHOLD,
    )

    hits = _prepare_blast_hits(df)

    mode = repeat_mode.strip().upper()
    try_ltr = mode in {"BOTH", "LTR"}
    try_tir = mode in {"BOTH", "TIR"}

    if try_ltr:
        df_ltr = hits[
            (hits["qstart"] != hits["sstart"])
            & (hits["qend"] != hits["send"])
            & hits["s_forward"]
            & (hits["sstart"] > hits["qstart"])
            # Left and right LTR must not overlap (right starts after left ends).
            & (hits["sstart"] > hits["qend"])
        ].copy()
        df_ltr = _filter_blast_hits(
            df_ltr,
            min_len=ltr_cfg.min_length,
            max_len=ltr_cfg.max_length,
            similarity_threshold=ltr_cfg.similarity_threshold,
        )
        if not df_ltr.empty:
            # Outermost = largest span from 5′ start to 3′ end.
            df_ltr["score"] = df_ltr["send"] - df_ltr["qstart"]
            ltr_hits = _select_terminal_pairs(
                df_ltr,
                record_len,
                cfg=ltr_cfg,
                score_column="score",
            )
            if ltr_hits:
                return "LTR", ltr_hits

    if not try_tir:
        return None, []

    df_tir = hits[
        (hits["qstart"] != hits["sstart"])
        & (hits["qend"] != hits["send"])
        & hits["s_inverted"]
        & (hits["s_lo"] > hits["qstart"])
        # Terminal copies must not overlap.
        & (hits["s_lo"] > hits["qend"])
    ].copy()
    df_tir = _filter_blast_hits(
        df_tir,
        min_len=tir_cfg.min_length,
        max_len=tir_cfg.max_length,
        similarity_threshold=tir_cfg.similarity_threshold,
    )
    if df_tir.empty:
        return None, []

    # TEtrimmer uses sstart - qstart; with normalized coords use the far subject end.
    df_tir["score"] = df_tir["s_hi"] - df_tir["qstart"]
    tir_hits = _select_terminal_pairs(
        df_tir,
        record_len,
        cfg=tir_cfg,
        score_column="score",
    )
    if tir_hits:
        return "TIR", tir_hits
    return None, []


def _terminal_repeat_identity_suffix(hit: pd.Series) -> str:
    if "pident" not in hit.index:
        return ""
    try:
        pident = float(hit["pident"])
    except (TypeError, ValueError):
        return ""
    if pident <= 0:
        return ""
    if abs(pident - round(pident)) < 0.05:
        return f" ({int(round(pident))}%)"
    return f" ({pident:.1f}%)"


def _terminal_repeat_base_name(prefix: str, pair_index: int) -> str:
    """``LTR5``, ``LTR5_2``, ``LTR5_3`` (1-based pair index)."""
    if pair_index <= 1:
        return prefix
    return f"{prefix}_{pair_index}"


def _terminal_repeat_name(base: str, hit: pd.Series) -> str:
    suffix = _terminal_repeat_identity_suffix(hit)
    return f"{base}{suffix}" if suffix else base


def _features_from_hit(
    repeat_type: str,
    hit: pd.Series,
    *,
    pair_index: int = 1,
) -> list[SequenceFeature]:
    q5 = int(hit["qstart"])
    q3 = int(hit["qend"])
    s_lo = int(hit["s_lo"])
    s_hi = int(hit["s_hi"])

    if repeat_type == "LTR":
        return [
            SequenceFeature(
                start=q5,
                end=q3,
                name=_terminal_repeat_name(
                    _terminal_repeat_base_name("LTR5", pair_index), hit
                ),
                feature_type="LTR",
            ),
            SequenceFeature(
                start=s_lo,
                end=s_hi,
                name=_terminal_repeat_name(
                    _terminal_repeat_base_name("LTR3", pair_index), hit
                ),
                feature_type="LTR",
            ),
        ]

    if repeat_type == "TIR":
        return [
            SequenceFeature(
                start=q5,
                end=q3,
                name=_terminal_repeat_name(
                    _terminal_repeat_base_name("TIR5", pair_index), hit
                ),
                feature_type="TIR",
                strand="+",
            ),
            SequenceFeature(
                start=s_lo,
                end=s_hi,
                name=_terminal_repeat_name(
                    _terminal_repeat_base_name("TIR3", pair_index), hit
                ),
                feature_type="TIR",
                strand="-",
            ),
        ]

    return []


def detect_terminal_repeat_features(
    df: pd.DataFrame,
    record_len: int,
    *,
    repeat_mode: str = "both",
    ltr_settings: TerminalRepeatSettings | None = None,
    tir_settings: TerminalRepeatSettings | None = None,
) -> TerminalRepeatResult:
    """Parse a self-BLAST table and return LTR or TIR annotation features."""
    mode = repeat_mode.strip().upper()
    if df.empty:
        return TerminalRepeatResult(
            features=(),
            repeat_type=None,
            message="Self-BLAST returned no hits (sequence may contain too many N/X).",
        )

    repeat_type, hit_rows = _check_terminal_repeat(
        df,
        record_len,
        repeat_mode=mode,
        ltr_settings=ltr_settings,
        tir_settings=tir_settings,
    )
    if repeat_type is None or not hit_rows:
        hits = _prepare_blast_hits(df)
        if mode == "TIR":
            tir_candidates = int(
                (
                    (hits["qend"] - hits["qstart"] >= MIN_TIR_LEN)
                    & hits["s_inverted"]
                ).sum()
            )
            if tir_candidates > 0:
                message = (
                    "No TIR passed terminal-end checks "
                    f"({tir_candidates} inverted hit(s) were too short or not terminal)."
                )
            else:
                message = "No TIR terminal repeats were identified."
        elif mode == "LTR":
            message = "No LTR terminal repeats were identified."
        else:
            tir_candidates = int(
                (
                    (hits["qend"] - hits["qstart"] >= MIN_TIR_LEN)
                    & hits["s_inverted"]
                ).sum()
            )
            if tir_candidates > 0:
                message = (
                    "No LTR or TIR passed terminal-end checks "
                    f"({tir_candidates} inverted hit(s) were too short or not terminal)."
                )
            else:
                message = "No LTR or TIR terminal repeats were identified."
        return TerminalRepeatResult(features=(), repeat_type=None, message=message)

    features: list[SequenceFeature] = []
    for index, hit in enumerate(hit_rows, start=1):
        features.extend(_features_from_hit(repeat_type, hit, pair_index=index))
    if mode == "LTR" and repeat_type != "LTR":
        return TerminalRepeatResult(features=(), repeat_type=None, message="No LTR terminal repeats were identified.")
    if mode == "TIR" and repeat_type != "TIR":
        return TerminalRepeatResult(features=(), repeat_type=None, message="No TIR terminal repeats were identified.")
    n_pairs = len(hit_rows)
    pair_label = "pair" if n_pairs == 1 else "pairs"
    return TerminalRepeatResult(
        features=tuple(features),
        repeat_type=repeat_type,
        message=(
            f"Detected {repeat_type} repeats "
            f"({n_pairs} {pair_label}, {len(features)} features)."
        ),
    )


def _blastn_executable(blast_dir: Path, os_type: str, *, use_system_blast: bool) -> str:
    if use_system_blast and check_cmd_in_path("blastn"):
        return "blastn"
    if os_type == "Linux":
        return str(blast_dir / "blastn_linux")
    if os_type == "Darwin":
        return str(blast_dir / "blastn_mac")
    return str(blast_dir / "blastn.exe")


def _run_cancellable_self_blastn(
    fasta_path: Path,
    database_path: str | None,
    work_dir: Path,
    *,
    cancel: JobCancellation | None = None,
    use_system_blast: bool = False,
    word_size: int = DEFAULT_TERMINAL_REPEAT_WORD_SIZE,
    num_threads: int = 1,
) -> Path:
    os_type = platform.system()
    script_dir = Path(get_original_file_path())
    blast_dir = script_dir / "blast"
    blast_raw = work_dir / f"{fasta_path.name}_no_header.b"
    from repeatviewer.backends.blast_output_names import local_blast_table_basename

    blast_header = work_dir / local_blast_table_basename(fasta_path.name, "blastn")
    blastn = _blastn_executable(blast_dir, os_type, use_system_blast=use_system_blast)
    effective_ws = max(MIN_TERMINAL_REPEAT_WORD_SIZE, int(word_size))
    cmd = [
        blastn,
        "-query",
        str(fasta_path),
        # No database: ``-subject`` aligns the sequence against itself directly.
        *(["-db", database_path] if database_path else ["-subject", str(fasta_path)]),
        "-max_target_seqs",
        "10000",
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore",
        "-evalue",
        "0.05",
        "-out",
        str(blast_raw),
        # blastn ignores (and warns about) -num_threads with -subject.
        *(["-num_threads", str(max(1, int(num_threads)))] if database_path else []),
        "-word_size",
        str(effective_ws),
        "-gapopen",
        "5",
        "-gapextend",
        "2",
        "-reward",
        "2",
        "-penalty",
        "-3",
    ]
    run_subprocess_with_cancel(cmd, cancel=cancel)
    check_cancelled(cancel)
    header = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "mismatch",
        "gapopen",
        "qstart",
        "qend",
        "sstart",
        "send",
        "sstrand",
        "evalue",
        "bitscore",
    ]
    if not blast_raw.is_file() or blast_raw.stat().st_size == 0:
        pd.DataFrame(columns=header).to_csv(blast_header, sep="\t", index=False)
        try:
            blast_raw.unlink(missing_ok=True)
        except OSError:
            pass
        return blast_header

    df = pd.read_csv(blast_raw, delimiter="\t", header=None, names=header)
    df.to_csv(blast_header, sep="\t", index=False)
    try:
        blast_raw.unlink()
    except OSError:
        pass
    return blast_header


def run_self_blast(
    fasta_path: Path,
    work_dir: Path,
    *,
    cancel: JobCancellation | None = None,
    use_system_blast: bool = False,
    word_size: int = DEFAULT_TERMINAL_REPEAT_WORD_SIZE,
    num_threads: int = 1,
) -> Path:
    """Self-BLAST one sequence; return the headered BLAST table path.

    Uses ``blastn -subject`` instead of building a database with makeblastdb:
    the hits are identical for a single sequence, and it saves one process
    launch per sequence (about 3x faster on macOS, where launches are slow).
    ``num_threads`` is ignored because blastn does not thread ``-subject`` runs.
    """
    del num_threads
    check_cancelled(cancel)
    fasta_path = fasta_path.expanduser().resolve()
    work_dir.mkdir(parents=True, exist_ok=True)
    return _run_cancellable_self_blastn(
        fasta_path,
        None,
        work_dir,
        cancel=cancel,
        use_system_blast=use_system_blast,
        word_size=word_size,
    )


def auto_annotate_terminal_repeats(
    fasta_path: Path,
    *,
    cancel: JobCancellation | None = None,
    use_system_blast: bool = False,
    work_dir: Path | None = None,
    cleanup: bool = True,
    repeat_mode: str = "both",
    ltr_settings: TerminalRepeatSettings | None = None,
    tir_settings: TerminalRepeatSettings | None = None,
    num_threads: int = 1,
    blast_cache: dict[int, pd.DataFrame] | None = None,
    shared_word_size: int | None = None,
) -> TerminalRepeatResult:
    """Run self-BLAST on a single-sequence FASTA and detect LTR or TIR features.

    ``blast_cache`` (word size -> hit table) lets a caller that runs LTR and
    TIR on the same sequence reuse one self-BLAST. ``shared_word_size`` makes
    both modes use that word size so the cached table can be shared.
    """
    fasta_path = fasta_path.expanduser().resolve()
    try:
        record_len = _record_length(fasta_path)
    except Exception as exc:
        raise RuntimeError(
            f"Could not read sequence for terminal-repeat detection: {exc}"
        ) from exc
    if record_len <= 0:
        return TerminalRepeatResult(
            features=(),
            repeat_type=None,
            message=(
                "Sequence is empty (or only gaps after MSA ungapping); "
                "skipped LTR/TIR self-BLAST."
            ),
        )
    ltr_cfg = ltr_settings or TerminalRepeatSettings(
        min_length=MIN_LTR_LEN,
        similarity_threshold=DEFAULT_LTR_SIMILARITY_THRESHOLD,
    )
    tir_cfg = tir_settings or TerminalRepeatSettings(
        min_length=MIN_TIR_LEN,
        similarity_threshold=DEFAULT_TIR_SIMILARITY_THRESHOLD,
    )
    ltr_word_size = word_size_for_terminal_repeat(ltr_cfg.min_length)
    tir_word_size = word_size_for_terminal_repeat(tir_cfg.min_length)
    mode = repeat_mode.strip().upper()

    temp_dir: tempfile.TemporaryDirectory[str] | None = None
    if work_dir is None:
        temp_dir = tempfile.TemporaryDirectory(prefix="repeatviewer_autoann_")
        work_dir = Path(temp_dir.name)
    else:
        work_dir = work_dir.expanduser().resolve()
        work_dir.mkdir(parents=True, exist_ok=True)

    blast_threads = max(1, int(num_threads))

    def _self_blast_df(subdir: str, word_size: int) -> pd.DataFrame:
        if shared_word_size is not None:
            word_size = int(shared_word_size)
        if blast_cache is not None and word_size in blast_cache:
            return blast_cache[word_size].copy()
        table = run_self_blast(
            fasta_path,
            work_dir / subdir if subdir else work_dir,
            cancel=cancel,
            use_system_blast=use_system_blast,
            word_size=word_size,
            num_threads=blast_threads,
        )
        frame = _load_self_blast_table(table)
        if blast_cache is not None:
            blast_cache[word_size] = frame.copy()
        return frame

    try:
        if mode == "LTR":
            df = _self_blast_df("", ltr_word_size)
            return detect_terminal_repeat_features(
                df,
                record_len,
                repeat_mode="LTR",
                ltr_settings=ltr_cfg,
                tir_settings=tir_cfg,
            )

        if mode == "TIR":
            df = _self_blast_df("tir_self_blast", tir_word_size)
            return detect_terminal_repeat_features(
                df,
                record_len,
                repeat_mode="TIR",
                ltr_settings=ltr_cfg,
                tir_settings=tir_cfg,
            )

        # BOTH: LTR and TIR word sizes are independent — separate self-BLAST when they differ.
        if ltr_word_size == tir_word_size:
            blast_table = run_self_blast(
                fasta_path,
                work_dir,
                cancel=cancel,
                use_system_blast=use_system_blast,
                word_size=ltr_word_size,
                num_threads=blast_threads,
            )
            df = _load_self_blast_table(blast_table)
            return detect_terminal_repeat_features(
                df,
                record_len,
                repeat_mode="BOTH",
                ltr_settings=ltr_cfg,
                tir_settings=tir_cfg,
            )

        ltr_table = run_self_blast(
            fasta_path,
            work_dir / "ltr_self_blast",
            cancel=cancel,
            use_system_blast=use_system_blast,
            word_size=ltr_word_size,
            num_threads=blast_threads,
        )
        ltr_df = _load_self_blast_table(ltr_table)
        ltr_result = detect_terminal_repeat_features(
            ltr_df,
            record_len,
            repeat_mode="LTR",
            ltr_settings=ltr_cfg,
            tir_settings=tir_cfg,
        )
        if ltr_result.repeat_type == "LTR":
            return ltr_result

        tir_table = run_self_blast(
            fasta_path,
            work_dir / "tir_self_blast",
            cancel=cancel,
            use_system_blast=use_system_blast,
            word_size=tir_word_size,
            num_threads=blast_threads,
        )
        tir_df = _load_self_blast_table(tir_table)
        return detect_terminal_repeat_features(
            tir_df,
            record_len,
            repeat_mode="TIR",
            ltr_settings=ltr_cfg,
            tir_settings=tir_cfg,
        )
    except JobCancelledError:
        raise
    except Exception:
        logging.exception("Terminal repeat auto annotation failed for %s", fasta_path.name)
        raise
    finally:
        if cleanup and temp_dir is not None:
            temp_dir.cleanup()
        elif cleanup and work_dir is not None and work_dir.name.startswith("repeatviewer_autoann_"):
            shutil.rmtree(work_dir, ignore_errors=True)


def merge_terminal_repeat_features(
    existing: list[SequenceFeature],
    detected: list[SequenceFeature] | tuple[SequenceFeature, ...],
) -> list[SequenceFeature]:
    """Replace existing LTR/TIR features with newly detected terminal repeats."""
    kept = [f for f in existing if f.feature_type not in {"LTR", "TIR"}]
    return kept + list(detected)
