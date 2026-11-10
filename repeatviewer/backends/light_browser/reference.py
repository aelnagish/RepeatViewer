"""Pick the reference FASTA a SAM / BAM / GFF file is drawn against.

The project's genome (or pan-genomes) comes first, then the genomes listed in
Manage Database, then any FASTA the user picked for this file before. A
reference *matches* a track when every sequence name in the file is found in it
(``chr1`` / ``1`` / ``scaffold_1`` forms count) and known lengths agree
(BAM/SAM ``@SQ`` must equal the FASTA; annotation spans may not run past it).

Pure Python, no Qt: tested in ``tests/test_light_browser.py``.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

from repeatviewer.backends.light_browser import genome_tracks as gt


MIN_MATCH_FRACTION = 0.5


@dataclass
class ReferenceMatch:
    path: Path
    label: str              # "Project genome", "Pan-genome Col0", "Database: TAIR10"…
    matched: int            # track sequences found in the reference
    total: int              # track sequences
    lengths: dict[str, int]  # reference sequence name -> length
    # Matched names where both the file and the reference know a length, and how
    # many of those lengths disagree (BAM/SAM @SQ, or max feature end vs FASTA).
    length_checked: int = 0
    length_mismatches: int = 0

    @property
    def fraction(self) -> float:
        return self.matched / self.total if self.total else 0.0

    @property
    def lengths_agree(self) -> bool:
        return self.length_mismatches == 0

    @property
    def ok(self) -> bool:
        # Every sequence name in the file must be found in the reference (lenient
        # forms like scaffold_1 ~ Chr1 are fine), and known lengths must agree.
        # A few coincidental name hits are not enough (see MIN_MATCH_FRACTION for
        # ranking only when comparing partial candidates).
        return (
            self.matched > 0
            and self.matched == self.total
            and self.lengths_agree
        )


def reference_candidates(project=None, remembered: Path | None = None) -> list[tuple[str, Path]]:
    """``(label, fasta)`` to try, best first, without duplicates or missing files."""
    out: list[tuple[str, Path]] = []
    seen: set[str] = set()

    def add(label: str, path) -> None:
        if path is None:
            return
        p = Path(path).expanduser()
        try:
            key = str(p.resolve())
        except OSError:
            key = str(p)
        if key in seen or not p.is_file() or p.name.lower().endswith(".gz"):
            return
        seen.add(key)
        out.append((label, p))

    if remembered is not None:
        add("Chosen for this file", remembered)
    if project is not None:
        try:
            targets = project.genome_targets()
        except Exception:  # noqa: BLE001 - an odd project must not break the preview
            targets = []
        for name, path in targets:
            add(f"Pan-genome {name}" if name else "Project genome", path)
        for g in getattr(project, "managed_genomes", None) or []:
            add(f"Database: {g.name}", g.path)
    return out


def _norm(name: str) -> str:
    parts = name.split()
    s = parts[0] if parts else name
    s = re.sub(r"^(chr|chromosome|scaffold)[_\-]?", "", s, flags=re.IGNORECASE)
    return re.sub(r"\.\d+$", "", s).lower()


class NameMatcher:
    """``genome_tracks.match_chrom`` for many names at once (one dict, not a scan per
    name: genomes with tens of thousands of scaffolds stay fast)."""

    def __init__(self, names) -> None:
        self.names = set(names)
        by_norm: dict[str, list[str]] = {}
        for n in self.names:
            by_norm.setdefault(_norm(n), []).append(n)
        self._norm = {k: v[0] for k, v in by_norm.items() if len(v) == 1}

    def __call__(self, wanted: str) -> str | None:
        if wanted in self.names:
            return wanted
        return self._norm.get(_norm(wanted))


def track_sequence_names(track) -> list[str]:
    try:
        return list(track.chrom_names())
    except Exception:  # noqa: BLE001
        return []


def track_sequence_lengths(track) -> tuple[dict[str, int], frozenset[str]]:
    """Lengths for the track's sequences, plus which came from a real header.

    Prefers BAM/SAM ``@SQ`` / BAM index lengths (authoritative). Otherwise uses
    the furthest feature / read end (estimate only: a mismatch is then only
    when that end runs past the reference chromosome).
    """
    names = track_sequence_names(track)
    out: dict[str, int] = {n: 0 for n in names}
    header: dict[str, int] = dict(getattr(track, "_lengths", None) or {})
    if not header and str(getattr(track, "path", "")).lower().endswith(gt.SAM_SUFFIXES):
        header = sam_header_lengths(track.path)
    authoritative = frozenset(n for n, length in header.items() if int(length or 0) > 0)
    for n in names:
        length = int(header.get(n, 0) or 0)
        if not length:
            ix = track._chroms.get(n) if hasattr(track, "_chroms") else None
            if ix is not None and len(ix):
                length = int(ix.ends.max())
        out[n] = max(0, length)
    return out, authoritative


def match_reference(
    fasta: Path,
    names: list[str],
    label: str = "",
    *,
    track_lengths: dict[str, int] | None = None,
    authoritative_lengths: frozenset[str] | set[str] | None = None,
) -> ReferenceMatch:
    """How many of ``names`` the FASTA holds (builds its .fai on first use).

    When ``track_lengths`` is given, matched chromosomes are length-checked:
    authoritative (BAM/SAM) lengths must equal the reference; feature-span
    estimates only fail when they run past the reference end.
    """
    index = gt.fasta_index(fasta)
    lengths = {k: e.length for k, e in index.items()}
    find = NameMatcher(lengths)
    matched = 0
    length_checked = 0
    length_mismatches = 0
    file_lengths = track_lengths or {}
    auth = frozenset(authoritative_lengths or ())
    for n in names:
        hit = find(n)
        if hit is None:
            continue
        matched += 1
        file_len = int(file_lengths.get(n, 0) or 0)
        ref_len = int(lengths.get(hit, 0) or 0)
        if file_len <= 0 or ref_len <= 0:
            continue
        length_checked += 1
        if n in auth:
            if file_len != ref_len:
                length_mismatches += 1
        elif file_len > ref_len:
            length_mismatches += 1
    return ReferenceMatch(
        Path(fasta),
        label or Path(fasta).name,
        matched,
        len(names),
        lengths,
        length_checked=length_checked,
        length_mismatches=length_mismatches,
    )


def choose_reference(
    names: list[str],
    candidates: list[tuple[str, Path]],
    *,
    track_lengths: dict[str, int] | None = None,
    authoritative_lengths: frozenset[str] | set[str] | None = None,
) -> tuple[ReferenceMatch | None, list[ReferenceMatch]]:
    """The candidate holding the largest share of ``names`` (at least
    ``MIN_MATCH_FRACTION``; earlier candidates win ties, a full match stops the
    search), and every candidate tried."""
    tried: list[ReferenceMatch] = []
    best: ReferenceMatch | None = None
    for label, path in candidates:
        try:
            m = match_reference(
                path,
                names,
                label,
                track_lengths=track_lengths,
                authoritative_lengths=authoritative_lengths,
            )
        except (OSError, gt.TrackError):
            continue
        tried.append(m)
        if m.ok and (best is None or m.fraction > best.fraction):
            best = m
        if best is not None and best.fraction >= 1.0 and best.lengths_agree:
            break
    return best, tried


def sam_header_lengths(path: Path) -> dict[str, int]:
    """``@SQ`` names and lengths of a SAM file (empty when it has no header)."""
    out: dict[str, int] = {}
    try:
        with gt._open_text(Path(path)) as fh:
            for line in fh:
                if not line.startswith("@"):
                    break
                if line.startswith("@SQ"):
                    tags = dict(f.split(":", 1) for f in line.rstrip("\n").split("\t")[1:]
                                if ":" in f)
                    if "SN" in tags:
                        try:
                            out[tags["SN"]] = int(tags.get("LN", "0"))
                        except ValueError:
                            out[tags["SN"]] = 0
    except OSError:
        pass
    return out


def sequence_lengths(track, reference: ReferenceMatch | None) -> dict[str, int]:
    """Track sequence name -> length, from the reference when it holds the name, then
    the file's own header (BAM / SAM), then the furthest item end."""
    names = track_sequence_names(track)
    out: dict[str, int] = {}
    header: dict[str, int] = dict(getattr(track, "_lengths", None) or {})
    if not header and str(track.path).lower().endswith(gt.SAM_SUFFIXES):
        header = sam_header_lengths(track.path)
    find = NameMatcher(reference.lengths) if reference is not None else None
    for n in names:
        length = 0
        if find is not None:
            hit = find(n)
            if hit is not None:
                length = reference.lengths.get(hit, 0)
        if not length:
            length = header.get(n, 0)
        if not length:
            ix = track._chroms.get(n) if hasattr(track, "_chroms") else None
            if ix is not None and len(ix):
                length = int(ix.ends.max())
        out[n] = max(1, length)
    return out


def genome_order(track, reference: ReferenceMatch | None) -> list[str]:
    """The track's sequences in the genome's own order: the reference FASTA's record
    order when it holds them, else the file's header (BAM / SAM ``@SQ``), else the
    order they first appear in the file. Never sorted by name."""
    names = track_sequence_names(track)
    header: list[str] = list(getattr(track, "_lengths", None) or [])
    if not header and str(track.path).lower().endswith(gt.SAM_SUFFIXES):
        header = list(sam_header_lengths(track.path))
    rank_in_header = {n: i for i, n in enumerate(header)}
    rank_in_file = {n: i for i, n in enumerate(names)}
    ref_rank: dict[str, int] = {}
    if reference is not None:
        find = NameMatcher(reference.lengths)
        order = {n: i for i, n in enumerate(reference.lengths)}
        for n in names:
            hit = find(n)
            if hit is not None:
                ref_rank[n] = order[hit]
    big = len(names) + len(header) + 1

    def key(n: str):
        if n in ref_rank:
            return (0, ref_rank[n])
        if n in rank_in_header:
            return (1, rank_in_header[n])
        return (2, rank_in_file.get(n, big))

    return sorted(names, key=key)
