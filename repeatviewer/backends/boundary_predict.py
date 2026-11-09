"""MSA boundary prediction ported from TEtrimmer ``boundarycrop.final_MSA``.

The open MSA is uncleaned, so it is first cleaned exactly like TEtrimmer
(``CleanAndSelectColumn(threshold=0.01).clean_column`` then
``remove_gaps_with_similarity_check``). Boundaries are then defined on the
cleaned MSA by:

* MSA identity  — ``CropEnd`` → ``crop_end_and_remove_gap`` → ``DefineBoundary``
* BLAST coverage — consensus vs the project genome (``GenomeBlastCoverage``)
* LTR / TIR      — self-BLAST of the cropped 0.7 consensus
* poly-A         — LINE / SINE only

TEtrimmer's MSA extension step is skipped: the open MSA is used as is. All
positions are mapped back to 0-based columns of the open MSA (end exclusive).

Developer: Jiangzhao Qian
"""

from __future__ import annotations

import os
import platform
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Sequence

import numpy as np
import pandas as pd

from repeatviewer.backends.job_cancel import JobCancelledError

os_type = platform.system()

# TEtrimmer "conserved" preset (config.json) values used by final_MSA.
CROP_END_THR = 0.7
CROP_END_WIN = 40
CROP_END_GAP_THR = 0.1
CROP_END_GAP_WIN = 250
POLY_PATTERNS = "A"
POLY_LEN = 10

# Rows TEtrimmer adds to the RepeatViewer export (its input sequence and its
# output consensus); they were not part of TEtrimmer's boundary MSA.
TETRIMMER_ROW_PREFIXES = ("TEtrimmer_input_", "TEtrimmer_output_")



def te_type_from_name(name: str | Path) -> str:
    """"LINE" / "SINE" when the MSA file name says so, else "Other".

    Like TEtrimmer, only LINE and SINE change the boundary path (poly-A, and
    for LINE the gap end-crop); every other type runs the default path.
    """
    stem = Path(str(name)).name
    if "LINE" in stem:
        return "LINE"
    if "SINE" in stem:
        return "SINE"
    return "Other"

_NUCL = ("a", "c", "g", "t")


@dataclass
class BoundaryPrediction:
    """Boundary columns on the open MSA: 0-based start, exclusive end."""

    msa_identity: tuple[int, int] | None = None
    blast_coverage: tuple[int, int] | None = None
    terminal_repeat: tuple[int, int] | None = None
    terminal_repeat_type: str | None = None  # "LTR" / "TIR"
    te_type: str | None = None  # path used: "LINE" / "SINE" / "Other"
    poly_a_end: int | None = None
    final: tuple[int, int] | None = None
    evaluation: str | None = None
    sequence_count: int = 0
    cleaned_columns: int = 0
    full_length_hits: int = 0
    notes: list[str] = field(default_factory=list)


# TEtrimmer / RepeatViewer sidecar: ``{msa_stem}.boundaries_RV.json``.
_BOUNDARIES_SUFFIX = ".boundaries_RV.json"
_BOUNDARIES_FORMAT_VERSION = 1


def boundaries_sidecar_path(fasta_path: Path | str) -> Path:
    """Hidden guides next to the MSA: ``{stem}.boundaries_RV.json``."""
    path = Path(fasta_path).expanduser().resolve()
    return path.parent / f"{path.stem}{_BOUNDARIES_SUFFIX}"


def is_boundaries_sidecar(path: Path | str) -> bool:
    return Path(path).name.endswith(_BOUNDARIES_SUFFIX)


def _parse_span(value) -> tuple[int, int] | None:
    if not isinstance(value, dict):
        return None
    try:
        start = int(value["start"])
        end = int(value["end"])
    except (KeyError, TypeError, ValueError):
        return None
    return start, end


def load_boundary_prediction_json(
    fasta_path: Path | str,
    *,
    msa_width: int | None = None,
) -> BoundaryPrediction | None:
    """Load TEtrimmer method boundaries beside ``fasta_path``, or ``None``.

    When ``msa_width`` is set, reject sidecars whose ``msa_width`` disagrees.
    """
    import json

    path = boundaries_sidecar_path(fasta_path)
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError, TypeError, ValueError):
        return None
    if not isinstance(payload, dict):
        return None
    stored_width = payload.get("msa_width")
    if msa_width is not None and stored_width is not None:
        try:
            if int(stored_width) != int(msa_width):
                return None
        except (TypeError, ValueError):
            return None
    term = payload.get("terminal_repeat")
    term_span = _parse_span(term)
    term_type = None
    if isinstance(term, dict) and term_span is not None:
        raw_type = term.get("type")
        term_type = str(raw_type) if raw_type else None
    hits = payload.get("full_length_hits")
    try:
        full_hits = int(hits) if hits is not None else 0
    except (TypeError, ValueError):
        full_hits = 0
    evaluation = payload.get("evaluation")
    return BoundaryPrediction(
        msa_identity=_parse_span(payload.get("msa_identity")),
        blast_coverage=_parse_span(payload.get("blast_coverage")),
        terminal_repeat=term_span,
        terminal_repeat_type=term_type,
        final=_parse_span(payload.get("final")),
        evaluation=str(evaluation) if evaluation else None,
        full_length_hits=full_hits,
        notes=["Loaded from TEtrimmer boundaries_RV.json"],
    )


#####################################################################################################
# Code block: alignment helpers (TEtrimmer reads/writes FASTA between steps; keep rows in memory)
#####################################################################################################
def _to_matrix(sequences: Sequence[str]) -> np.ndarray:
    return np.array([list(s) for s in sequences], dtype="<U1")


def _to_rows(matrix: np.ndarray) -> list[str]:
    return ["".join(row) for row in matrix]


def _check_cancel(cancel) -> None:
    if cancel is not None and bool(getattr(cancel, "requested", False)):
        raise JobCancelledError("Boundary prediction stopped by user.")


#####################################################################################################
# Code block: MSA cleaning (MSAcluster.CleanAndSelectColumn, functions.remove_gaps_with_similarity_check)
#####################################################################################################
def clean_column(matrix: np.ndarray, threshold: float = 0.01) -> np.ndarray:
    """``CleanAndSelectColumn.clean_column``: rare nucleotides (< threshold) → gap.

    Proportions count lowercase a/c/g/t only; the replacement test is
    case-sensitive on lowercase letters, as in TEtrimmer.
    """
    out = matrix.copy()
    if out.size == 0:
        return out
    lower = np.char.lower(matrix)
    counts = np.stack([(lower == n).sum(axis=0) for n in _NUCL])  # (4, cols)
    total = counts.sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        props = np.where(total > 0, counts / np.maximum(total, 1), 0.0)
    for k, n in enumerate(_NUCL):
        rare = (matrix == n) & (props[k] < threshold)[None, :]
        out[rare] = "-"
    return out


def _calc_conservation(col: np.ndarray) -> float:
    lower = np.char.lower(col)
    counts = [int((lower == n).sum()) for n in _NUCL]
    total = sum(counts)
    return max(counts) / total if total > 0 else 0


def remove_gaps_with_similarity_check(
    matrix: np.ndarray,
    gap_threshold: float = 0.8,
    simi_check_gap_thre: float = 0.4,
    similarity_threshold: float = 0.7,
    min_nucleotide: int = 5,
) -> tuple[np.ndarray, dict[int, int]] | tuple[None, None]:
    """Keep columns by gap fraction / conservation; ``None`` under 50 columns.

    ``column_mapping`` maps kept index → input index, plus one extra key
    ``len(keep)`` → last + 1 (TEtrimmer's exclusive-end convention).
    """
    n_rows, n_cols = matrix.shape
    keep_list: list[int] = []
    column_mapping: dict[int, int] = {}
    gap_counts = (matrix == "-").sum(axis=0)
    for col_idx in range(n_cols):
        gap_count = int(gap_counts[col_idx])
        gap_fraction = gap_count / n_rows
        nt_count = n_rows - gap_count
        if nt_count < min_nucleotide:
            continue
        if gap_fraction <= gap_threshold:
            if simi_check_gap_thre > gap_fraction:
                keep_list.append(col_idx)
                column_mapping[len(keep_list) - 1] = col_idx
            elif simi_check_gap_thre <= gap_fraction:
                if _calc_conservation(matrix[:, col_idx]) < similarity_threshold:
                    continue
                keep_list.append(col_idx)
                column_mapping[len(keep_list) - 1] = col_idx
    if len(keep_list) < 50:
        return None, None
    column_mapping[len(keep_list)] = column_mapping[len(keep_list) - 1] + 1
    return matrix[:, keep_list], column_mapping


#####################################################################################################
# Code block: consensus (functions.con_generater / con_generater_no_file)
#####################################################################################################
def consensus(matrix: np.ndarray, threshold: float = 0.8, ambiguous: str = "N") -> str:
    upper = np.char.upper(matrix)
    order = ("A", "C", "G", "T")
    counts = np.stack([(upper == b).sum(axis=0) for b in order])  # (4, cols)
    total = counts.sum(axis=0)
    # max(counts, key=counts.get) keeps the first of tied bases in A, C, G, T order.
    winner = counts.argmax(axis=0)
    out = []
    for i in range(matrix.shape[1]):
        n = int(total[i])
        if n == 0:
            out.append(ambiguous)
            continue
        freq = int(counts[winner[i], i]) / n
        out.append(order[winner[i]] if freq >= threshold else ambiguous)
    return "".join(out)


#####################################################################################################
# Code block: CropEnd (boundaryclass.CropEnd)
#####################################################################################################
def _first_window_over(
    values: np.ndarray, window: int, threshold: float, starts: np.ndarray, width: int
) -> int | None:
    """First of ``starts`` whose window ``values[s:s+width]`` sums above threshold.

    Uses prefix sums to skip clear cases, then the same ``ndarray.sum`` pandas
    uses for borderline windows so float results match TEtrimmer exactly.
    """
    if starts.size == 0:
        return None
    csum = np.concatenate(([0.0], np.cumsum(values)))
    approx = csum[starts + width] - csum[starts]
    tol = 1e-6
    candidates = np.nonzero(approx > threshold - tol)[0]
    for idx in candidates:
        s = int(starts[idx])
        if approx[idx] > threshold + tol or values[s : s + width].sum() > threshold:
            return s
    return None


def crop_end(
    matrix: np.ndarray,
    threshold: float = 16,
    window_size: int = 20,
    crop_l: bool = True,
    crop_r: bool = True,
) -> np.ndarray:
    n_rows, n_cols = matrix.shape
    lower = np.char.lower(matrix)
    counts = np.stack([(lower == n).sum(axis=0) for n in _NUCL])
    total = counts.sum(axis=0)
    # Columns with fewer than 5 nucleotides give proportion 0.
    with np.errstate(divide="ignore", invalid="ignore"):
        props = np.where(total >= 5, counts / np.maximum(total, 1), 0.0)
    prop_mat = np.full((n_rows, n_cols), np.nan)
    for k, n in enumerate(_NUCL):
        hit = matrix == n  # case-sensitive lookup, as in TEtrimmer
        prop_mat[hit] = np.broadcast_to(props[k], (n_rows, n_cols))[hit]
    # DataFrame.round(2) is numpy round; Series.sum skips NaN (zero-filled sum).
    prop_mat = np.round(prop_mat, 2)
    values_all = np.where(np.isnan(prop_mat), 0.0, prop_mat)

    out = matrix.copy()
    w = int(window_size)
    for r in range(n_rows):
        values = values_all[r]
        start = 0
        end = 0
        if crop_l:
            found = _first_window_over(
                values, w, threshold, np.arange(0, max(0, n_cols - w + 1)), w
            )
            if found is not None:
                start = found
        if crop_r:
            # Window row[i - w : i + 1] holds w + 1 values; scanned right to left.
            idx = np.arange(n_cols - 1, w - 1, -1)
            found = _first_window_over(values, w, threshold, idx - w, w + 1)
            if found is not None:
                end = found + w + 1
        else:
            end = n_cols
        # No passing window keeps end at 0, which gaps the whole row.
        out[r, :start] = "-"
        out[r, end:] = "-"
    return out


def crop_end_by_gap(
    matrix: np.ndarray, gap_threshold: float = 0.05, window_size: int = 300
) -> np.ndarray:
    """``boundaryclass.CropEndByGap`` (both ends)."""
    n_rows, n_cols = matrix.shape
    w = int(window_size)
    out = matrix.copy()
    for r in range(n_rows):
        gaps = np.concatenate(([0], np.cumsum(matrix[r] == "-")))
        start = 0
        end = 0
        for i in range(n_cols - w + 1):
            if (gaps[i + w] - gaps[i]) / w <= gap_threshold:
                start = i
                break
        for i in range(n_cols - 1, w - 2, -1):
            if (gaps[i + 1] - gaps[i - w + 1]) / w <= gap_threshold:
                end = i + 1 if i + 1 > start else start
                break
        out[r, :start] = "-"
        out[r, end:] = "-"
    return out


#####################################################################################################
# Code block: DefineBoundary (boundaryclass.DefineBoundary, polyA_position=None)
#####################################################################################################
def define_boundary(matrix: np.ndarray, check_window: int = 4, max_X: float = 0) -> tuple[int, int]:
    con = consensus(matrix, threshold=0.8, ambiguous="N")
    n = len(con)
    start_post = None
    for i in range(n):
        if i + check_window <= n:
            if con[i : i + check_window].count("N") <= max_X * check_window:
                start_post = i
                break
    if start_post is None:
        start_post = n
    end_post = None
    for i in range(n - 1, -1, -1):
        i = i + 1
        if i - check_window >= 0:
            if con[i - check_window : i].count("N") <= max_X * check_window:
                end_post = i
                break
    if end_post is None:
        end_post = 1
    return start_post, end_post


#####################################################################################################
# Code block: poly-A (functions.find_poly_a_end_position)
#####################################################################################################
def find_poly_a_end_position(sequence: str, poly_patterns: str = "A", min_length: int = 10):
    sequence = sequence.upper()
    seq_length = len(sequence)
    seq_mid_position = seq_length - int(seq_length // 2.5)
    patterns = [p.strip().upper() for p in poly_patterns.split(",") if p.strip()]
    for pattern in patterns:
        pattern_len = len(pattern)
        if pattern_len == 1:
            poly_length = 0
            interruption_used = False
            for i in range(seq_length - 1, seq_mid_position - 1, -1):
                if sequence[i] == pattern:
                    poly_length += 1
                elif not interruption_used and 5 <= poly_length < min_length:
                    interruption_used = True
                    poly_length += 1
                else:
                    if poly_length >= min_length:
                        return i + poly_length
                    poly_length = 0
                    interruption_used = False
        else:
            poly_length = 0
            i = seq_length - pattern_len
            while i >= seq_mid_position:
                if sequence[i : i + pattern_len] == pattern:
                    poly_length += 1
                    if poly_length >= min_length:
                        return i + pattern_len
                    i -= pattern_len
                else:
                    poly_length = 0
                    i -= 1
    return None


def _clip_con(start: int, stop: int, seq: str) -> str:
    return seq[max(0, start) : min(len(seq), stop)]


#####################################################################################################
# Code block: BLAST runners (TEtrimmer blast() / check_terminal_repeat command lines)
#####################################################################################################
def _blast_executable(name: str, use_system_blast: bool) -> str:
    from repeatviewer.backends.gui_functions import check_cmd_in_path, get_original_file_path

    if use_system_blast and check_cmd_in_path(name):
        return name
    blast_dir = os.path.join(get_original_file_path(), "blast")
    if os_type == "Linux":
        return os.path.join(blast_dir, f"{name}_linux")
    if os_type == "Darwin":
        return os.path.join(blast_dir, f"{name}_mac")
    return os.path.join(blast_dir, f"{name}.exe")


def _write_fasta(path: Path, seq_id: str, seq: str) -> None:
    path.write_text(f">{seq_id}\n{seq}\n", encoding="utf-8")


def _run(cmd: list[str], tool: str, inputs, cwd: Path, cancel) -> None:
    from repeatviewer.external_process import run_logged_process

    run_logged_process(cmd, tool=tool, inputs=inputs, cwd=str(cwd), check=True, cancel=cancel)


class GenomeBlastCoverage:
    """``boundaryclass.GenomeBlastCoverage`` against a prepared BLAST database."""

    def __init__(self, con_seq: str, database: str, work_dir: Path, *, num_threads=1,
                 use_system_blast=False, cancel=None, label="con"):
        self.con_seq = con_seq
        self.database = database
        self.work_dir = work_dir
        self.num_threads = num_threads
        self.use_system_blast = use_system_blast
        self.cancel = cancel
        self.label = label
        self.coverage_list: list[int] | None = None
        self.start_point = None
        self.end_point = None
        self.blast_df: pd.DataFrame | None = None

    def calculate_blast_coverage(self) -> bool:
        query = self.work_dir / f"{self.label}.fa"
        _write_fasta(query, self.label, self.con_seq)
        out = self.work_dir / f"{self.label}.b"
        # blast(min_length=50, blast_qcov_hsp_perc=1, evalue='1e-8'); default task (megablast).
        _run(
            [
                _blast_executable("blastn", self.use_system_blast),
                "-max_target_seqs", "10000",
                "-query", str(query),
                "-db", self.database,
                "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send sstrand evalue qcovhsp",
                "-evalue", "1e-8",
                "-qcov_hsp_perc", "1",
                "-num_threads", str(max(1, int(self.num_threads))),
                "-out", str(out),
            ],
            "blastn", [query], self.work_dir, self.cancel,
        )
        rows = []
        for line in out.read_text(encoding="utf-8").splitlines():
            parts = line.split("\t")
            if len(parts) < 11 or not (int(parts[3]) > 50):
                continue
            # .b.bed columns used later: alnlen ($4), qstart ($7), qend ($8).
            rows.append((int(parts[3]), int(parts[6]), int(parts[7])))
        if not rows:
            return False
        query_length = len(self.con_seq)
        self.blast_df = pd.DataFrame(rows, columns=["alnlen", "qstart", "qend"])
        self.blast_df["weight"] = 1
        top_10_indices = self.blast_df.nlargest(10, "alnlen").index
        self.blast_df.loc[top_10_indices, "weight"] = 3
        qstart = self.blast_df["qstart"].clip(lower=1, upper=query_length).astype(int).to_numpy()
        qend = self.blast_df["qend"].clip(lower=1, upper=query_length).astype(int).to_numpy()
        weights = self.blast_df["weight"].astype(int).to_numpy()
        diff = np.zeros(query_length + 1, dtype=int)
        for s, e, w in zip(qstart, qend, weights):
            diff[min(s, e) - 1] += w
            diff[max(s, e)] -= w
        self.coverage_list = np.cumsum(diff[:-1]).tolist()
        return True

    def find_boundary_blast_coverage(self, window_size=20, threshold=100):
        cov = self.coverage_list or []
        cov_len = len(cov)
        for i in range(cov_len - window_size + 1):
            window = cov[i : i + window_size]
            if sum(window) > threshold:
                for j, k in enumerate(window):
                    if k >= 5:
                        self.start_point = i + j
                        break
                if self.start_point is not None:
                    break
        for i in range(cov_len - 1, window_size - 2, -1):
            window = cov[i - window_size + 1 : i + 1]
            if sum(window) > threshold:
                for j in range(window_size - 1, -1, -1):
                    if window[j] >= 5:
                        self.end_point = i - (window_size - 1 - j) + 1
                        break
                if self.end_point is not None:
                    break
        # Trim low-coverage tails at a >= 5x jump.
        if self.start_point is not None:
            start_limit = min(self.start_point + 20, cov_len - 1)
            for i in range(self.start_point, start_limit):
                if cov[i + 1] >= max(1, cov[i]) * 5:
                    self.start_point = i + 1
                    break
        if self.end_point is not None:
            end_limit = max(1, self.end_point - 10)
            for i in range(self.end_point - 1, end_limit - 1, -1):
                if cov[i - 1] >= max(1, cov[i]) * 5:
                    self.end_point = i
                    break
        full_length_hit_count = 0
        lower_percent_hit_count = 0
        if self.start_point is not None and self.end_point is not None and self.blast_df is not None:
            seq_start = self.start_point + 1
            seq_end = self.end_point
            full_length_buffer = 20
            lower_buffer = max(20, int(abs(seq_end - seq_start) * 0.15))
            df = self.blast_df
            full_length_hit_count = len(
                df[(df["qstart"] <= seq_start + full_length_buffer) & (df["qend"] >= seq_end - full_length_buffer)]
            )
            lower_percent_hit_count = len(
                df[(df["qstart"] <= seq_start + lower_buffer) & (df["qend"] >= seq_end - lower_buffer)]
            )
        return self.start_point, self.end_point, full_length_hit_count, lower_percent_hit_count


def check_terminal_repeat(con_seq: str, work_dir: Path, *, use_system_blast=False, cancel=None,
                          TIR_adj=5000, LTR_adj=5000):
    """``functions.check_terminal_repeat`` with max_pairs=1.

    Returns ``(LTR_boundary, TIR_boundary)``; each is ``[start0, end]`` or None.
    """
    record_len = len(con_seq)
    query = work_dir / "tem_con.fa"
    _write_fasta(query, "tem_con", con_seq)
    db = work_dir / "Tem_blast_database"
    _run(
        [_blast_executable("makeblastdb", use_system_blast), "-in", str(query),
         "-dbtype", "nucl", "-out", str(db)],
        "makeblastdb", [query], work_dir, cancel,
    )
    out = work_dir / "tem_blast_out.b"
    _run(
        [_blast_executable("blastn", use_system_blast), "-query", str(query), "-db", str(db),
         "-outfmt", "6 qseqid qstart qend sstart send pident",
         "-evalue", "0.05", "-word_size", "11", "-gapopen", "5", "-gapextend", "2",
         "-reward", "2", "-penalty", "-3", "-out", str(out)],
        "blastn", [query], work_dir, cancel,
    )
    text = out.read_text(encoding="utf-8")
    if text.strip() == "":
        return None, None
    rows = [line.split("\t") for line in text.splitlines() if line.strip()]
    df = pd.DataFrame(
        [[r[0], int(r[1]), int(r[2]), int(r[3]), int(r[4]), float(r[5])] for r in rows]
    )
    if df.empty:
        return None, None

    def _select(df_hits, score, adj):
        ranked = df_hits.copy()
        ranked["_span_score"] = score
        ranked = ranked.sort_values("_span_score", ascending=False, kind="mergesort")
        for _, row in ranked.iterrows():
            if abs(float(row["_span_score"])) < (record_len - adj):
                continue
            return row  # max_pairs = 1: the outermost passing hit
        return None

    df_LTR = df[
        (df.iloc[:, 2] - df.iloc[:, 1] >= 100)
        & (df.iloc[:, 1] != df.iloc[:, 3])
        & (df.iloc[:, 2] != df.iloc[:, 4])
        & (df.iloc[:, 3] < df.iloc[:, 4])
        & (df.iloc[:, 3] > df.iloc[:, 1])
        # Left LTR [qstart,qend] and right LTR [sstart,send] must not overlap.
        & (df.iloc[:, 3] > df.iloc[:, 2])
    ].copy()
    if not df_LTR.empty:
        hit = _select(df_LTR, df_LTR.iloc[:, 4] - df_LTR.iloc[:, 1], LTR_adj)
        # A direct-repeat HSP blocks TIR, even when none pass the span test.
        if hit is None:
            return None, None
        return [int(hit.iloc[1]) - 1, int(hit.iloc[4])], None

    df_TIR = df[
        (df.iloc[:, 2] - df.iloc[:, 1] >= 20)
        & (df.iloc[:, 1] != df.iloc[:, 3])
        & (df.iloc[:, 2] != df.iloc[:, 4])
        & (df.iloc[:, 3] > df.iloc[:, 4])
        & (df.iloc[:, 4] > df.iloc[:, 1])
        # Left TIR must end before the right (inverted) copy starts.
        & (df.iloc[:, 4] > df.iloc[:, 2])
    ].copy()
    if not df_TIR.empty:
        hit = _select(df_TIR, df_TIR.iloc[:, 3] - df_TIR.iloc[:, 1], TIR_adj)
        if hit is not None:
            return None, [int(hit.iloc[1]) - 1, int(hit.iloc[3])]
    return None, None


#####################################################################################################
# Code block: final_MSA boundary definition
#####################################################################################################
def boundary_rows(record_ids: Sequence[str], sequences: Sequence[str]) -> list[int]:
    """Row indices used for prediction (drops TEtrimmer's input/output rows)."""
    return [
        i for i, rid in enumerate(record_ids)
        if not str(rid).startswith(TETRIMMER_ROW_PREFIXES)
    ]


def clean_msa(matrix: np.ndarray):
    """Lowercase (MAFFT case), column clean, gap removal → (gap_sim, column_map)."""
    lowered = np.char.lower(matrix)
    cleaned = clean_column(lowered, threshold=0.01)
    return remove_gaps_with_similarity_check(cleaned)


def predict_boundaries(
    record_ids: Sequence[str],
    sequences: Sequence[str],
    *,
    genome_path: Path | str | None,
    te_type: str = "Other",
    num_threads: int = 1,
    use_system_blast: bool = False,
    cancel=None,
    progress: Callable[[int, str], None] | None = None,
) -> BoundaryPrediction:
    """Run TEtrimmer's ``final_MSA`` boundary logic on the open MSA."""
    result = BoundaryPrediction()

    def _step(n: int, text: str) -> None:
        _check_cancel(cancel)
        if progress is not None:
            progress(n, text)

    rows = boundary_rows(record_ids, sequences)
    if len(rows) < 5:
        raise ValueError("Boundary prediction needs at least 5 aligned sequences.")
    seqs = [str(sequences[i]) for i in rows]
    width = len(seqs[0])
    if any(len(s) != width for s in seqs):
        raise ValueError("MSA rows must share one length.")
    if len(rows) < len(sequences):
        result.notes.append(
            f"Ignored {len(sequences) - len(rows)} TEtrimmer input/output row(s)."
        )
    result.sequence_count = len(rows)
    result.te_type = te_type
    is_line = te_type.upper().startswith("LINE")
    is_sine = te_type.upper().startswith("SINE")

    _step(1, "cleaning MSA columns")
    gap_sim, column_mapping_initial = clean_msa(_to_matrix(seqs))
    if gap_sim is None:
        raise ValueError("Fewer than 50 columns remain after TEtrimmer column cleaning.")
    result.cleaned_columns = gap_sim.shape[1]

    def to_raw_start(pos: int) -> int:
        return int(column_mapping_initial[pos])

    def to_raw_end(pos: int) -> int:
        # Exclusive end on gap_sim → exclusive end on the open MSA.
        return int(column_mapping_initial[pos - 1]) + 1 if pos > 0 else to_raw_start(0)

    def to_raw(span):
        return (to_raw_start(span[0]), to_raw_end(span[1]))

    con_08 = consensus(gap_sim, threshold=0.8)
    con_07 = consensus(gap_sim, threshold=0.7)
    con_055 = consensus(gap_sim, threshold=0.55)
    gap_sim_cp = crop_end(gap_sim, threshold=20, window_size=40)
    cp_con_070 = consensus(gap_sim_cp, threshold=0.7)

    _step(2, "MSA identity")
    if is_line:
        gap_sim_cp = crop_end_by_gap(gap_sim_cp, CROP_END_GAP_THR, CROP_END_GAP_WIN)
    cropped = crop_end(gap_sim_cp, threshold=CROP_END_THR * CROP_END_WIN, window_size=CROP_END_WIN)
    cropped_gs, column_mapping_MSA = remove_gaps_with_similarity_check(cropped)
    if cropped_gs is None:
        raise ValueError("Fewer than 50 columns remain after end cropping.")
    s_db, e_db = define_boundary(cropped_gs, check_window=4, max_X=0)
    start_posit_MSA = column_mapping_MSA[s_db]
    end_posit_MSA = column_mapping_MSA[e_db]
    result.msa_identity = to_raw((start_posit_MSA, end_posit_MSA))

    _step(3, "LTR / TIR self-BLAST")
    with tempfile.TemporaryDirectory(prefix="rv_boundary_") as tmp:
        work = Path(tmp)
        ltr, tir = check_terminal_repeat(
            cp_con_070, work, use_system_blast=use_system_blast, cancel=cancel
        )
        left_rep = right_rep = None
        if ltr is not None:
            left_rep, right_rep = ltr
            result.terminal_repeat_type = "LTR"
        elif tir is not None:
            left_rep, right_rep = tir
            result.terminal_repeat_type = "TIR"
        if left_rep is not None and right_rep is not None:
            result.terminal_repeat = to_raw((left_rep, right_rep))

        poly_a = None
        if is_line or is_sine:
            poly_a = find_poly_a_end_position(con_08, POLY_PATTERNS, POLY_LEN)
            beyond_ns = 1
            if poly_a is not None:
                beyond_ns = _clip_con(poly_a, poly_a + 100, con_08).count("N")
            if beyond_ns < 50:
                poly_a = None
        if poly_a is not None:
            result.poly_a_end = to_raw_end(poly_a)

        _step(4, "BLAST coverage vs genome")
        genome = Path(genome_path) if genome_path else None
        if genome is None or not genome.is_file():
            result.notes.append("No project genome: BLAST coverage and final boundary skipped.")
            return result
        from repeatviewer.backends.gui_functions import check_database

        database = check_database(
            str(genome), os_type=os_type, use_system_blast=use_system_blast, cancel=cancel
        )
        if database in ("makeblastdb_not_found", "makeblastdb_got_error"):
            raise RuntimeError("Could not build the genome BLAST database.")

        def _coverage(con_seq: str, label: str):
            obj = GenomeBlastCoverage(
                con_seq, str(database), work, num_threads=num_threads,
                use_system_blast=use_system_blast, cancel=cancel, label=label,
            )
            if not obj.calculate_blast_coverage():
                return obj, None
            return obj, obj.find_boundary_blast_coverage()

        initial_con = cp_con_070 if is_line else con_07
        cov_obj, found = _coverage(initial_con, "con_07")
        if found is None:
            result.notes.append("No genome BLAST hits for the consensus.")
            return result
        start_cov, end_cov, full_hits, lower_hits = found
        if start_cov is None or end_cov is None or start_cov >= end_cov:
            result.notes.append("BLAST coverage did not define a boundary.")
            return result
        if not is_line:
            if start_cov - start_posit_MSA >= 30 or end_posit_MSA - end_cov >= 30:
                _step(5, "BLAST coverage (0.55 consensus)")
                cov_obj, found = _coverage(con_055, "con_055")
                if found is None:
                    result.notes.append("No genome BLAST hits for the 0.55 consensus.")
                    return result
                start_cov, end_cov, full_hits, lower_hits = found
        if start_cov is None or end_cov is None or start_cov >= end_cov:
            result.notes.append("BLAST coverage did not define a boundary.")
            return result

    result.blast_coverage = to_raw((start_cov, end_cov))
    result.full_length_hits = int(full_hits)

    final_start = start_cov
    final_end = poly_a if poly_a is not None else end_cov
    result.final = to_raw((final_start, final_end))

    evaluation = "Need_check"
    diff_start = abs(start_posit_MSA - start_cov)
    diff_end = abs(end_posit_MSA - end_cov)
    end_boundary_buffer = len(con_08) - 40
    msa_seq_n = len(rows)
    if (
        diff_start <= 15 and diff_end <= 15 and full_hits >= 3
        and start_cov >= 40 and end_cov <= end_boundary_buffer
    ):
        evaluation = "Perfect"
    elif (
        start_cov >= 40 and end_cov <= end_boundary_buffer
        and (
            full_hits >= 2
            or (
                (left_rep is not None and (left_rep - start_cov <= 40))
                and (right_rep is not None and (right_rep - end_cov <= 40))
            )
        )
    ):
        evaluation = "Good"
    elif msa_seq_n >= 20 and start_cov >= 40 and end_cov <= end_boundary_buffer and lower_hits >= 2:
        evaluation = "Reco_check"
    result.evaluation = evaluation
    return result
