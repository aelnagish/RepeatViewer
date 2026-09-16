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
    check_database,
    get_original_file_path,
)
from repeatviewer.backends.job_cancel import JobCancelledError, JobCancellation, check_cancelled, run_subprocess_with_cancel
from repeatviewer.backends.auto_annotation_settings import (
    DEFAULT_LTR_SIMILARITY_THRESHOLD,
    DEFAULT_TIR_SIMILARITY_THRESHOLD,
    TerminalRepeatSettings,
)
from repeatviewer.backends.sequence_features import SequenceFeature

TIR_ADJ = 50
LTR_ADJ = 100
MIN_LTR_LEN = 100
MIN_TIR_LEN = 20
DEFAULT_TERMINAL_REPEAT_WORD_SIZE = 11
MIN_TERMINAL_REPEAT_WORD_SIZE = 4


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


def _pair_distance_ok(hit: pd.Series, *, min_pair_distance: int) -> bool:
    if min_pair_distance <= 0:
        return True
    q3 = int(hit["qend"])
    s5 = int(hit["s_lo"])
    return s5 - q3 >= min_pair_distance


def _check_terminal_repeat(
    df: pd.DataFrame,
    record_len: int,
    *,
    repeat_mode: str = "both",
    tir_adj: int = TIR_ADJ,
    ltr_adj: int = LTR_ADJ,
    ltr_settings: TerminalRepeatSettings | None = None,
    tir_settings: TerminalRepeatSettings | None = None,
) -> tuple[str | None, pd.Series | None]:
    """
    Port of TEtrimmer ``check_terminal_repeat`` hit selection.

    Returns (repeat_type, winning_row) where repeat_type is ``LTR``, ``TIR``, or None.
    Unlike TEtrimmer, TIR is still evaluated when LTR-like hits exist but fail the
    terminal-span threshold.
    """
    if df.empty:
        return None, None

    ltr_cfg = ltr_settings or TerminalRepeatSettings(
        min_length=MIN_LTR_LEN,
        similarity_threshold=DEFAULT_LTR_SIMILARITY_THRESHOLD,
    )
    tir_cfg = tir_settings or TerminalRepeatSettings(
        min_length=MIN_TIR_LEN,
        similarity_threshold=DEFAULT_TIR_SIMILARITY_THRESHOLD,
    )

    hits = _prepare_blast_hits(df)
    span_min_ltr = max(0, record_len - ltr_adj)
    span_min_tir = max(0, record_len - tir_adj)

    mode = repeat_mode.strip().upper()
    try_ltr = mode in {"BOTH", "LTR"}
    try_tir = mode in {"BOTH", "TIR"}

    if try_ltr:
        df_ltr = hits[
            (hits["qstart"] != hits["sstart"])
            & (hits["qend"] != hits["send"])
            & hits["s_forward"]
            & (hits["sstart"] > hits["qstart"])
        ].copy()
        df_ltr = _filter_blast_hits(
            df_ltr,
            min_len=ltr_cfg.min_length,
            max_len=ltr_cfg.max_length,
            similarity_threshold=ltr_cfg.similarity_threshold,
        )
        if not df_ltr.empty:
            df_ltr["score"] = df_ltr["send"] - df_ltr["qstart"]
            ltr_hit = df_ltr.loc[df_ltr["score"].idxmax()]
            if (
                abs(int(ltr_hit["send"]) - int(ltr_hit["qstart"])) >= span_min_ltr
                and _pair_distance_ok(ltr_hit, min_pair_distance=ltr_cfg.min_pair_distance)
            ):
                return "LTR", ltr_hit

    if not try_tir:
        return None, None

    df_tir = hits[
        (hits["qstart"] != hits["sstart"])
        & (hits["qend"] != hits["send"])
        & hits["s_inverted"]
        & (hits["s_lo"] > hits["qstart"])
    ].copy()
    df_tir = _filter_blast_hits(
        df_tir,
        min_len=tir_cfg.min_length,
        max_len=tir_cfg.max_length,
        similarity_threshold=tir_cfg.similarity_threshold,
    )
    if df_tir.empty:
        return None, None

    # TEtrimmer uses sstart - qstart; with normalized coords use the far subject end.
    df_tir["score"] = df_tir["s_hi"] - df_tir["qstart"]
    tir_hit = df_tir.loc[df_tir["score"].idxmax()]
    if (
        abs(int(tir_hit["s_hi"]) - int(tir_hit["qstart"])) >= span_min_tir
        and _pair_distance_ok(tir_hit, min_pair_distance=tir_cfg.min_pair_distance)
    ):
        return "TIR", tir_hit
    return None, None


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


def _terminal_repeat_name(base: str, hit: pd.Series) -> str:
    suffix = _terminal_repeat_identity_suffix(hit)
    return f"{base}{suffix}" if suffix else base


def _features_from_hit(repeat_type: str, hit: pd.Series) -> list[SequenceFeature]:
    q5 = int(hit["qstart"])
    q3 = int(hit["qend"])
    s_lo = int(hit["s_lo"])
    s_hi = int(hit["s_hi"])

    if repeat_type == "LTR":
        return [
            SequenceFeature(
                start=q5,
                end=q3,
                name=_terminal_repeat_name("LTR5", hit),
                feature_type="LTR",
            ),
            SequenceFeature(
                start=s_lo,
                end=s_hi,
                name=_terminal_repeat_name("LTR3", hit),
                feature_type="LTR",
            ),
        ]

    if repeat_type == "TIR":
        return [
            SequenceFeature(
                start=q5,
                end=q3,
                name=_terminal_repeat_name("TIR5", hit),
                feature_type="TIR",
                strand="+",
            ),
            SequenceFeature(
                start=s_lo,
                end=s_hi,
                name=_terminal_repeat_name("TIR3", hit),
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

    repeat_type, hit = _check_terminal_repeat(
        df,
        record_len,
        repeat_mode=mode,
        ltr_settings=ltr_settings,
        tir_settings=tir_settings,
    )
    if repeat_type is None or hit is None:
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
                    "No TIR passed terminal-span checks "
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
                    "No LTR or TIR passed TEtrimmer terminal-span checks "
                    f"({tir_candidates} inverted hit(s) were too short or not terminal)."
                )
            else:
                message = "No LTR or TIR terminal repeats were identified."
        return TerminalRepeatResult(features=(), repeat_type=None, message=message)

    features = _features_from_hit(repeat_type, hit)
    if mode == "LTR" and repeat_type != "LTR":
        return TerminalRepeatResult(features=(), repeat_type=None, message="No LTR terminal repeats were identified.")
    if mode == "TIR" and repeat_type != "TIR":
        return TerminalRepeatResult(features=(), repeat_type=None, message="No TIR terminal repeats were identified.")
    return TerminalRepeatResult(
        features=tuple(features),
        repeat_type=repeat_type,
        message=f"Detected {repeat_type} repeats ({len(features)} features).",
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
    database_path: str,
    work_dir: Path,
    *,
    cancel: JobCancellation | None = None,
    use_system_blast: bool = False,
    word_size: int = DEFAULT_TERMINAL_REPEAT_WORD_SIZE,
) -> Path:
    os_type = platform.system()
    script_dir = Path(get_original_file_path())
    blast_dir = script_dir / "blast"
    blast_raw = work_dir / f"{fasta_path.name}_no_header.b"
    blast_header = work_dir / f"{fasta_path.name}_blast.txt"
    blastn = _blastn_executable(blast_dir, os_type, use_system_blast=use_system_blast)
    effective_ws = max(MIN_TERMINAL_REPEAT_WORD_SIZE, int(word_size))
    cmd = [
        blastn,
        "-query",
        str(fasta_path),
        "-db",
        database_path,
        "-max_target_seqs",
        "10000",
        "-outfmt",
        "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue bitscore",
        "-evalue",
        "0.05",
        "-out",
        str(blast_raw),
        "-num_threads",
        "1",
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
) -> Path:
    """Build a database and run self-BLAST; return the headered BLAST table path."""
    check_cancelled(cancel)
    os_type = platform.system()
    fasta_path = fasta_path.expanduser().resolve()
    work_dir.mkdir(parents=True, exist_ok=True)

    database = check_database(
        str(fasta_path),
        str(work_dir),
        os_type=os_type,
        use_system_blast=use_system_blast,
        cancel=cancel,
    )
    if database == "makeblastdb_not_found":
        raise RuntimeError(
            "makeblastdb was not found. Install BLAST+ or use the bundled binaries."
        )
    if database == "makeblastdb_got_error":
        raise RuntimeError("makeblastdb failed while building the self-BLAST database.")

    check_cancelled(cancel)
    return _run_cancellable_self_blastn(
        fasta_path,
        database,
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
) -> TerminalRepeatResult:
    """Run self-BLAST on a single-sequence FASTA and detect LTR or TIR features."""
    fasta_path = fasta_path.expanduser().resolve()
    record_len = _record_length(fasta_path)
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

    try:
        if mode == "LTR":
            blast_table = run_self_blast(
                fasta_path,
                work_dir,
                cancel=cancel,
                use_system_blast=use_system_blast,
                word_size=ltr_word_size,
            )
            df = _load_self_blast_table(blast_table)
            return detect_terminal_repeat_features(
                df,
                record_len,
                repeat_mode="LTR",
                ltr_settings=ltr_cfg,
                tir_settings=tir_cfg,
            )

        if mode == "TIR":
            blast_table = run_self_blast(
                fasta_path,
                work_dir / "tir_self_blast",
                cancel=cancel,
                use_system_blast=use_system_blast,
                word_size=tir_word_size,
            )
            df = _load_self_blast_table(blast_table)
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
