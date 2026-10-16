"""Map features between gapped MSA columns and ungapped sequence coordinates.

Detectors / BLAST / transfer should run on gap-free DNA/protein. Features are
then projected back onto the alignment (``-`` / ``.`` columns) for MSA display.
"""

from __future__ import annotations

from collections.abc import Sequence
from dataclasses import dataclass, replace
from pathlib import Path

from repeatviewer.backends.sequence_features import SequenceFeature

# Same gap alphabet as motif Find ("Allow gaps in sequence").
GAP_CHARS = frozenset("-.")

# MSA residue paint modes (Display → Color).
# Schemes apply at overview and letter zoom. Overview cells the scheme does
# not colour use the MSA light-grey bar colour.
MSA_COLOR_ALL_BACKGROUNDS = "all_backgrounds"
MSA_COLOR_LETTERS_ONLY = "letters_only"
# Fill letters whose column frequency is below the divergence threshold.
MSA_COLOR_NO_BG_WHEN_CONSERVED = "no_bg_when_conserved"
# Fill letters whose column frequency meets the divergence threshold.
MSA_COLOR_BG_WHEN_CONSERVED = "bg_when_conserved"
# Fill the reference fully; other rows only where they differ from it.
MSA_COLOR_DIFF_FROM_SELECTED = "diff_from_selected"
# Retired id kept for settings migration only.
_MSA_COLOR_CONSERVATION_GREY_LEGACY = "conservation_grey"

MSA_COLOR_MODE_CHOICES: tuple[tuple[str, str], ...] = (
    (MSA_COLOR_LETTERS_ONLY, "No color"),
    (MSA_COLOR_ALL_BACKGROUNDS, "Color all"),
    (MSA_COLOR_BG_WHEN_CONSERVED, "Color conserved"),
    (MSA_COLOR_NO_BG_WHEN_CONSERVED, "Color divergent"),
    (MSA_COLOR_DIFF_FROM_SELECTED, "Color vs reference"),
)

MSA_COLOR_MODE_IDS = frozenset(mode_id for mode_id, _label in MSA_COLOR_MODE_CHOICES)

# Fraction of non-gap letters in a column required for a residue to count as conserved.
DEFAULT_MSA_DIVERGENCE_THRESHOLD = 0.4


def normalize_msa_color_mode(mode: str | None) -> str:
    text = str(mode or "").strip()
    if text == _MSA_COLOR_CONSERVATION_GREY_LEGACY:
        return MSA_COLOR_ALL_BACKGROUNDS
    if text in MSA_COLOR_MODE_IDS:
        return text
    return MSA_COLOR_ALL_BACKGROUNDS


def normalize_msa_divergence_threshold(value: float | int | None) -> float:
    try:
        return max(0.0, min(1.0, float(value)))
    except (TypeError, ValueError):
        return DEFAULT_MSA_DIVERGENCE_THRESHOLD


def column_residue_frequencies(sequences: list[str]) -> list[dict[str, float]]:
    """Per-column frequency of each non-gap residue (upper-case keys).

    Frequencies sum to 1.0 over non-gap letters. All-gap columns yield ``{}``.
    """
    if not sequences:
        return []
    width = len(sequences[0])
    if any(len(s) != width for s in sequences):
        raise ValueError("MSA rows must share one length")
    out: list[dict[str, float]] = []
    n_rows = len(sequences)
    for col in range(width):
        counts: dict[str, int] = {}
        total = 0
        for row in range(n_rows):
            ch = sequences[row][col]
            if ch in GAP_CHARS:
                continue
            up = ch.upper()
            counts[up] = counts.get(up, 0) + 1
            total += 1
        if total <= 0:
            out.append({})
        else:
            out.append({letter: count / total for letter, count in counts.items()})
    return out


def column_conservation_mask(
    sequences: list[str],
    *,
    threshold: float = 1.0,
) -> list[bool]:
    """True when at least one non-gap residue meets ``threshold`` in the column.

    ``threshold=1.0`` matches the old unanimous-column definition.
    """
    thresh = normalize_msa_divergence_threshold(threshold)
    freqs = column_residue_frequencies(sequences)
    return [
        any(freq >= thresh for freq in col_freqs.values()) for col_freqs in freqs
    ]


def letter_is_conserved(
    column_freqs: dict[str, float] | None,
    residue: str,
    *,
    threshold: float,
) -> bool:
    """Whether ``residue`` meets the column frequency threshold."""
    ch = str(residue or "").upper()
    if not ch or ch in GAP_CHARS:
        return False
    freqs = column_freqs or {}
    return float(freqs.get(ch, 0.0)) >= normalize_msa_divergence_threshold(threshold)


def residue_differs_from_reference(residue: str, reference: str | None) -> bool:
    """True when ``residue`` differs from the selected standard at one column.

    Matching letters (case-insensitive) do not highlight. Missing standard
    never highlights. A gap on the compared row highlights when the standard
    has a residue; a letter vs a gap on the standard also counts. Matching
    gaps (both gap) do not highlight.
    """
    if reference is None:
        return False
    ch = str(residue or "")
    ref = str(reference or "")
    ch_gap = (not ch) or ch in GAP_CHARS
    ref_gap = (not ref) or ref in GAP_CHARS
    if ch_gap and ref_gap:
        return False
    if ch_gap or ref_gap:
        return True
    return ch.upper() != ref.upper()


def msa_residue_fill_background(
    mode: str | None,
    *,
    residue: str,
    column_freqs: dict[str, float] | None = None,
    threshold: float = DEFAULT_MSA_DIVERGENCE_THRESHOLD,
    column_conserved: bool | None = None,
    reference_residue: str | None = None,
    is_reference_row: bool = False,
) -> bool:
    """Whether an MSA cell should paint a residue background fill.

    Gaps never take a residue colour fill except in
    :data:`MSA_COLOR_DIFF_FROM_SELECTED` (caller may still draw ``-``).
    ``column_conserved`` is accepted for older callers and ignored when
    ``column_freqs`` is provided.

    For :data:`MSA_COLOR_DIFF_FROM_SELECTED` (Color vs reference):
    - ``is_reference_row=True`` fills every non-gap residue (full reference).
    - other rows fill when ``residue`` differs from ``reference_residue``,
      including gaps where the reference has a letter.
    """
    ch = str(residue or "")
    normalized = normalize_msa_color_mode(mode)
    if normalized == MSA_COLOR_DIFF_FROM_SELECTED:
        if ch in GAP_CHARS:
            if is_reference_row:
                return False
            return residue_differs_from_reference(ch, reference_residue)
        if is_reference_row:
            return True
        return residue_differs_from_reference(ch, reference_residue)
    if ch in GAP_CHARS:
        return False
    if normalized == MSA_COLOR_ALL_BACKGROUNDS:
        return True
    if normalized == MSA_COLOR_LETTERS_ONLY:
        return False
    if column_freqs is not None:
        conserved = letter_is_conserved(
            column_freqs, ch, threshold=threshold
        )
    elif column_conserved is not None:
        conserved = bool(column_conserved)
    else:
        conserved = False
    if normalized == MSA_COLOR_NO_BG_WHEN_CONSERVED:
        # Divergent: fill letters that do not meet the threshold.
        return not conserved
    if normalized == MSA_COLOR_BG_WHEN_CONSERVED:
        return conserved
    return True


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


@dataclass
class MsaAnnotationWorkspace:
    """Temporary ungapped FASTA for detectors; project features back when done."""

    work_fasta: Path
    base_features: tuple[SequenceFeature, ...]
    keep: list[int] | None
    _tmp: object | None = None

    def project_features(
        self, features: list[SequenceFeature] | tuple[SequenceFeature, ...]
    ) -> list[SequenceFeature]:
        return project_features_to_gapped(list(features), self.keep)

    def close(self) -> None:
        tmp = self._tmp
        self._tmp = None
        if tmp is not None:
            close = getattr(tmp, "cleanup", None)
            if callable(close):
                close()


def prepare_msa_annotation_workspace(
    fasta_path: Path,
    base_features: Sequence[SequenceFeature] | None = None,
) -> MsaAnnotationWorkspace:
    """If ``fasta_path`` has gaps, write an ungapped temp FASTA for detectors.

    Existing feature coordinates are mapped into ungapped space. Callers must
    :meth:`MsaAnnotationWorkspace.project_features` on detector output and
    :meth:`MsaAnnotationWorkspace.close` when finished.
    """
    import tempfile

    from Bio import SeqIO
    from Bio.Seq import Seq
    from Bio.SeqRecord import SeqRecord

    resolved = Path(fasta_path).expanduser().resolve()
    base = tuple(base_features or ())
    try:
        records = list(SeqIO.parse(str(resolved), "fasta"))
    except (OSError, ValueError):
        records = []
    if not records:
        return MsaAnnotationWorkspace(
            work_fasta=resolved, base_features=base, keep=None
        )
    raw_seq = str(records[0].seq)
    ungapped, keep = ungap_sequence(raw_seq)
    if keep is None:
        return MsaAnnotationWorkspace(
            work_fasta=resolved, base_features=base, keep=None
        )
    tmp = tempfile.TemporaryDirectory(prefix="repeatviewer_msa_annot_")
    rid = str(records[0].id or "seq")
    # Unique basename avoids parallel makeblastdb races on a shared
    # ``ungapped.fa`` link name under the BLAST-safe temp directory.
    safe_rid = "".join(ch if ch.isalnum() or ch in "-_" else "_" for ch in rid)[:40]
    out = Path(tmp.name) / f"ungapped_{safe_rid or 'seq'}.fa"
    SeqIO.write(
        SeqRecord(Seq(ungapped), id=rid, name=rid, description=rid),
        str(out),
        "fasta",
    )
    return MsaAnnotationWorkspace(
        work_fasta=out,
        base_features=tuple(project_features_to_ungapped(list(base), keep)),
        keep=keep,
        _tmp=tmp,
    )


def remove_alignment_columns(
    sequences: list[str],
    start_1based: int,
    end_1based: int,
) -> list[str]:
    """Drop inclusive 1-based columns ``[start, end]`` from every MSA row."""
    if not sequences:
        return []
    width = len(sequences[0])
    if any(len(s) != width for s in sequences):
        raise ValueError("MSA rows must share one length")
    a = max(1, min(width, int(start_1based)))
    b = max(1, min(width, int(end_1based)))
    if b < a:
        a, b = b, a
    # 0-based slice ends exclusive
    lo, hi = a - 1, b
    return [s[:lo] + s[hi:] for s in sequences]


def splice_rows_region(
    sequences: list[str],
    start_1based: int,
    end_1based: int,
    row_indices: set[int] | list[int],
) -> list[str]:
    """Remove ``[start, end]`` only from ``row_indices``, then pad with gaps.

    Other rows are unchanged. All rows are right-padded with ``-`` so the MSA
    stays equal-length.
    """
    if not sequences:
        return []
    width = len(sequences[0])
    if any(len(s) != width for s in sequences):
        raise ValueError("MSA rows must share one length")
    a = max(1, min(width, int(start_1based)))
    b = max(1, min(width, int(end_1based)))
    if b < a:
        a, b = b, a
    lo, hi = a - 1, b
    rows = {int(i) for i in row_indices}
    out: list[str] = []
    for i, seq in enumerate(sequences):
        if i in rows:
            out.append(seq[:lo] + seq[hi:])
        else:
            out.append(seq)
    max_len = max((len(s) for s in out), default=0)
    if max_len <= 0:
        return out
    return [s + ("-" * (max_len - len(s))) for s in out]


# Residues accepted when typing into an MSA selection (AliView-style edit).
_MSA_NUCLEOTIDE_EDIT_CHARS = frozenset("ACGTUNIRYWSMKBDHVX-.")
_MSA_PROTEIN_EDIT_CHARS = frozenset("ACDEFGHIKLMNPQRSTVWYBXZJUO*-.")


def normalize_msa_edit_residue(residue: str, *, protein: bool = False) -> str | None:
    """Return a single editable MSA character, or ``None`` if invalid."""
    ch = str(residue or "")
    if len(ch) != 1:
        return None
    if ch in "-.":
        return ch
    up = ch.upper()
    allowed = _MSA_PROTEIN_EDIT_CHARS if protein else _MSA_NUCLEOTIDE_EDIT_CHARS
    if up in allowed:
        return up
    return None


def fill_rows_region(
    sequences: list[str],
    start_1based: int,
    end_1based: int,
    row_indices: set[int] | list[int] | None = None,
    *,
    fill: str = "-",
) -> list[str]:
    """Replace ``[start, end]`` with ``fill`` on the given rows (or all rows)."""
    if not sequences:
        return []
    width = len(sequences[0])
    if any(len(s) != width for s in sequences):
        raise ValueError("MSA rows must share one length")
    ch = str(fill or "")
    if len(ch) != 1:
        raise ValueError("fill must be a single character")
    a = max(1, min(width, int(start_1based)))
    b = max(1, min(width, int(end_1based)))
    if b < a:
        a, b = b, a
    lo, hi = a - 1, b
    block = ch * (hi - lo)
    if row_indices is None:
        rows = set(range(len(sequences)))
    else:
        rows = {int(i) for i in row_indices}
    out: list[str] = []
    for i, seq in enumerate(sequences):
        if i in rows:
            out.append(seq[:lo] + block + seq[hi:])
        else:
            out.append(seq)
    return out


def mask_rows_region(
    sequences: list[str],
    start_1based: int,
    end_1based: int,
    row_indices: set[int] | list[int] | None = None,
) -> list[str]:
    """Replace ``[start, end]`` with ``-`` on the given rows (or all rows)."""
    return fill_rows_region(
        sequences,
        start_1based,
        end_1based,
        row_indices,
        fill="-",
    )


def region_has_residue(
    sequence: str,
    start_1based: int,
    end_1based: int,
) -> bool:
    """True when ``sequence[start:end]`` (1-based inclusive) has a non-gap."""
    if not sequence:
        return False
    a = int(start_1based)
    b = int(end_1based)
    if b < a:
        a, b = b, a
    lo = max(0, a - 1)
    hi = min(len(sequence), b)
    for i in range(lo, hi):
        if sequence[i] not in GAP_CHARS:
            return True
    return False


def drop_features_overlapping_region(
    features: list[SequenceFeature],
    start_1based: int,
    end_1based: int,
    sequence: str | None = None,
) -> list[SequenceFeature]:
    """Drop features that overlap inclusive ``[start, end]``.

    When ``sequence`` is set, a feature is dropped only if the overlap with
    the edited span contains at least one non-gap residue. Overlaps that fall
    only on ``-`` / ``.`` columns are kept so deleting gap columns does not
    strip MSA annotations.
    """
    a = int(start_1based)
    b = int(end_1based)
    if b < a:
        a, b = b, a
    kept: list[SequenceFeature] = []
    for feat in features:
        fs, fe = int(feat.start), int(feat.end)
        if fe < fs:
            fs, fe = fe, fs
        if fe < a or fs > b:
            kept.append(feat)
            continue
        if sequence is None:
            # No sequence context → any overlap removes the feature.
            continue
        inter_a = max(fs, a)
        inter_b = min(fe, b)
        if region_has_residue(sequence, inter_a, inter_b):
            continue
        kept.append(feat)
    return kept


def remap_interval_after_column_delete(
    start: int,
    end: int,
    del_start: int,
    del_end: int,
) -> tuple[int, int] | None:
    """Map a 1-based inclusive interval after deleting ``[del_start, del_end]``.

    Returns ``None`` when the interval is entirely removed.
    """
    a = int(start)
    b = int(end)
    if b < a:
        a, b = b, a
    d0 = int(del_start)
    d1 = int(del_end)
    if d1 < d0:
        d0, d1 = d1, d0
    n = d1 - d0 + 1
    if b < d0:
        return a, b
    if a > d1:
        return a - n, b - n

    def map_start(p: int) -> int:
        if p < d0:
            return p
        if p > d1:
            return p - n
        return d0

    def map_end(p: int) -> int:
        if p < d0:
            return p
        if p > d1:
            return p - n
        return d0 - 1

    ns, ne = map_start(a), map_end(b)
    if ne < ns:
        return None
    return ns, ne


def remap_feature_after_column_delete(
    feat: SequenceFeature,
    del_start: int,
    del_end: int,
) -> SequenceFeature | None:
    """Shift/shrink a feature after MSA columns are deleted; drop if gone."""
    mapped = remap_interval_after_column_delete(
        feat.start, feat.end, del_start, del_end
    )
    if mapped is None:
        return None
    return replace(feat, start=mapped[0], end=mapped[1])


def remap_features_after_column_delete(
    features: list[SequenceFeature],
    del_start: int,
    del_end: int,
) -> list[SequenceFeature]:
    out: list[SequenceFeature] = []
    for feat in features:
        mapped = remap_feature_after_column_delete(feat, del_start, del_end)
        if mapped is not None:
            out.append(mapped)
    return out
