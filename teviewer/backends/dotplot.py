"""Fast self-comparison dotplots (the "Dotplot (Self)" preview tab).

Approach is the same as the EMBOSS ``dottup`` tool Geneious uses: index every
word (k-mer) of the sequence in a hash table, then read matching positions
straight out of it. That is roughly O(n + matches) instead of the O(n^2) of
comparing every position against every other.

The important second step is **collapsing runs into segments**. A 150 bp repeat
match is one diagonal line, not 140 separate dots. Storing segments rather than
points cuts the data by orders of magnitude and is what keeps panning and
zooming smooth -- the canvas redraws a few thousand lines, never millions of
points. Nothing is recomputed while zooming; only a change of word size or the
reverse-complement setting triggers a new pass.

This module is deliberately Qt-free so it can be tested and profiled on its own.
It does not touch the Flexidot or TEAid code paths.
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path

# Sequences longer than this are refused outright: a self dotplot is O(matches),
# and for long repetitive sequences the match count explodes well before the
# picture stays readable. The preview shows a message instead of hanging.
MAX_DOTPLOT_SEQUENCE_BYTES = 5_000_000

# Upper bound on raw k-mer matches before we give up and ask for a larger word.
# Satellite arrays are the pathological case: every repeat copy shares k-mers
# with every other copy, so matches grow with the square of the copy number.
DEFAULT_MAX_MATCHES = 4_000_000

MIN_WORD_SIZE = 4
MAX_WORD_SIZE = 64
DEFAULT_WORD_SIZE = 10

STRAND_FORWARD = "forward"
STRAND_REVERSE = "reverse"

_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


class DotplotError(RuntimeError):
    """Raised when a dotplot cannot be computed (too large, too many matches)."""


@dataclass(frozen=True)
class DotplotSegment:
    """One diagonal run of matching positions.

    ``x``/``y`` are the start coordinates in sequence space and ``length`` is
    how many bases the run covers. Forward segments run down-right (slope +1);
    reverse segments run up-right (slope -1), which is how inverted repeats --
    TIRs and hairpins -- show up.
    """

    x: int
    y: int
    length: int
    strand: str = STRAND_FORWARD

    @property
    def end_x(self) -> int:
        return self.x + self.length

    @property
    def end_y(self) -> int:
        if self.strand == STRAND_REVERSE:
            return self.y - self.length
        return self.y + self.length


@dataclass
class DotplotResult:
    sequence_id: str
    sequence_length: int
    word_size: int
    segments: list[DotplotSegment] = field(default_factory=list)
    include_reverse: bool = False
    raw_match_count: int = 0
    truncated: bool = False
    note: str = ""

    @property
    def segment_count(self) -> int:
        return len(self.segments)


def reverse_complement(sequence: str) -> str:
    return sequence.translate(_COMPLEMENT)[::-1]


def suggest_word_size(sequence_length: int) -> int:
    """A word size that keeps the match count sane for this length.

    Random sequence produces roughly n^2 / 4^k spurious matches, so the word
    has to grow with the sequence. These thresholds keep a random sequence
    under a few hundred thousand background hits.
    """
    if sequence_length <= 2_000:
        return 7
    if sequence_length <= 20_000:
        return 10
    if sequence_length <= 100_000:
        return 14
    if sequence_length <= 1_000_000:
        return 18
    return 22


def _build_word_index(sequence: str, word_size: int) -> dict[str, list[int]]:
    """Map each k-mer to the positions where it starts.

    Words containing anything other than A/C/G/T are skipped: runs of N would
    otherwise match each other everywhere and swamp the plot.
    """
    index: dict[str, list[int]] = defaultdict(list)
    limit = len(sequence) - word_size + 1
    valid = set("ACGT")
    for i in range(limit):
        word = sequence[i : i + word_size]
        # Cheap rejection first; set ops on every word would dominate runtime.
        if not valid.issuperset(word):
            continue
        index[word].append(i)
    return index


def _merge_runs(starts: list[int], word_size: int) -> list[tuple[int, int]]:
    """Merge overlapping/adjacent word hits on one diagonal into runs.

    Input is the list of x-positions of matching words on a single diagonal.
    Returns (start, length) pairs covering them.
    """
    if not starts:
        return []
    starts.sort()
    runs: list[tuple[int, int]] = []
    run_start = starts[0]
    run_end = run_start + word_size  # exclusive
    for pos in starts[1:]:
        if pos <= run_end:  # overlapping or exactly adjacent -> same line
            run_end = max(run_end, pos + word_size)
        else:
            runs.append((run_start, run_end - run_start))
            run_start = pos
            run_end = pos + word_size
    runs.append((run_start, run_end - run_start))
    return runs


def compute_self_dotplot(
    sequence: str,
    *,
    sequence_id: str = "",
    word_size: int = DEFAULT_WORD_SIZE,
    include_reverse: bool = True,
    include_main_diagonal: bool = False,
    max_matches: int = DEFAULT_MAX_MATCHES,
) -> DotplotResult:
    """Self-comparison dotplot as a list of diagonal segments.

    Only one triangle is produced for the forward strand: a self dotplot is
    symmetric, so the viewer mirrors segments across the main diagonal rather
    than storing both halves. ``include_main_diagonal`` adds the full-length
    self match, which the viewer normally draws itself.
    """
    seq = sequence.upper()
    n = len(seq)
    if n == 0:
        raise DotplotError("The sequence is empty.")
    if n > MAX_DOTPLOT_SEQUENCE_BYTES:
        raise DotplotError(
            f"Sequence is {n:,} bp. Self dotplots are limited to "
            f"{MAX_DOTPLOT_SEQUENCE_BYTES:,} bp because the number of matches "
            "grows faster than the sequence length."
        )
    word_size = max(MIN_WORD_SIZE, min(int(word_size), MAX_WORD_SIZE))
    if word_size > n:
        raise DotplotError(
            f"Word size ({word_size}) is longer than the sequence ({n} bp)."
        )

    result = DotplotResult(
        sequence_id=sequence_id,
        sequence_length=n,
        word_size=word_size,
        include_reverse=include_reverse,
    )

    index = _build_word_index(seq, word_size)
    raw = 0

    # ---- forward strand: group hits by diagonal (y - x), then merge runs ----
    diagonals: dict[int, list[int]] = defaultdict(list)
    for positions in index.values():
        count = len(positions)
        if count < 2:
            continue
        for a_i in range(count - 1):
            a = positions[a_i]
            for b in positions[a_i + 1 :]:
                diagonals[b - a].append(a)
                raw += 1
        if raw > max_matches:
            raise DotplotError(
                f"Too many matches (over {max_matches:,}) at word size "
                f"{word_size}. This usually means a highly repetitive sequence. "
                f"Try a word size of {suggest_word_size(n)} or larger."
            )

    for offset, starts in diagonals.items():
        for start, length in _merge_runs(starts, word_size):
            result.segments.append(
                DotplotSegment(x=start, y=start + offset, length=length)
            )

    if include_main_diagonal:
        result.segments.append(DotplotSegment(x=0, y=0, length=n))

    # ---- reverse strand: inverted repeats (TIRs, hairpins) ----
    if include_reverse:
        rc = reverse_complement(seq)
        rc_index = _build_word_index(rc, word_size)
        anti: dict[int, list[int]] = defaultdict(list)
        for word, positions in index.items():
            rc_positions = rc_index.get(word)
            if not rc_positions:
                continue
            for i in positions:
                for j in rc_positions:
                    # rc[j:j+k] is the reverse complement of seq[n-j-k : n-j]
                    y = n - j - word_size
                    anti[i + y].append(i)
                    raw += 1
            if raw > max_matches:
                result.truncated = True
                result.note = (
                    "Reverse-strand matches were truncated; try a larger word size."
                )
                break

        for key, starts in anti.items():
            for start, length in _merge_runs(starts, word_size):
                # Anti-diagonal: y decreases as x increases.
                result.segments.append(
                    DotplotSegment(
                        x=start,
                        y=key - start,
                        length=length,
                        strand=STRAND_REVERSE,
                    )
                )

    result.raw_match_count = raw
    return result


def compute_self_dotplot_for_fasta(
    path: Path,
    *,
    word_size: int | None = None,
    include_reverse: bool = True,
    max_matches: int = DEFAULT_MAX_MATCHES,
) -> DotplotResult:
    """Read a single-record FASTA and compute its self dotplot.

    Raises DotplotError with a readable message rather than letting Biopython
    exceptions reach the UI.
    """
    from teviewer.backends.sequence_features import load_single_sequence_fasta

    try:
        sequence_id, sequence = load_single_sequence_fasta(path)
    except (OSError, ValueError) as exc:
        raise DotplotError(str(exc)) from exc

    resolved_word = word_size or suggest_word_size(len(sequence))
    return compute_self_dotplot(
        sequence,
        sequence_id=sequence_id,
        word_size=resolved_word,
        include_reverse=include_reverse,
        max_matches=max_matches,
    )
