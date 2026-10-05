"""BLAST output filename tags (local bn/bp/bx; online obn/obp/obx/…)."""

from __future__ import annotations

from pathlib import Path

# Local BLAST: ``{query}_{tag}.txt`` (e.g. reads.fa_bn.txt).
LOCAL_BLAST_PROGRAM_TAGS: dict[str, str] = {
    "blastn": "bn",
    "blastp": "bp",
    "blastx": "bx",
}

# Online NCBI BLAST: ``{query}_{tag}.txt`` / ``.html``.
ONLINE_BLAST_PROGRAM_TAGS: dict[str, str] = {
    "blastn": "obn",
    "blastp": "obp",
    "blastx": "obx",
    "tblastn": "obtn",
    "tblastx": "obtx",
}

# Legacy local tags (still discovered, not written for new runs).
_LEGACY_LOCAL_BLAST_TAGS: dict[str, str] = {
    "blastn": "b",
    "blastp": "blastp",
    "blastx": "blastx",
}

_LEGACY_ONLINE_TAG = "ob"


def local_blast_run_tag(program: str) -> str:
    key = str(program or "blastn").lower()
    try:
        return LOCAL_BLAST_PROGRAM_TAGS[key]
    except KeyError as exc:
        raise ValueError(f"Unsupported local BLAST program: {program}") from exc


def online_blast_run_tag(program: str) -> str:
    key = str(program or "blastn").lower()
    try:
        return ONLINE_BLAST_PROGRAM_TAGS[key]
    except KeyError as exc:
        raise ValueError(f"Unsupported online BLAST program: {program}") from exc


def _table_suffixes_for_tag(tag: str) -> str:
    return f"_{tag}.txt"


def all_blast_table_txt_suffixes() -> tuple[str, ...]:
    """All ``*_*.txt`` endings recognized as BLAST hit tables (longest first)."""
    tags: list[str] = []
    tags.extend(ONLINE_BLAST_PROGRAM_TAGS.values())
    tags.append(_LEGACY_ONLINE_TAG)
    tags.append("online")
    tags.extend(LOCAL_BLAST_PROGRAM_TAGS.values())
    tags.extend(_LEGACY_LOCAL_BLAST_TAGS.values())
    tags.append("blast")
    # Longer tags first so ``_obtn.txt`` wins over ``_ob.txt``.
    tags = sorted(set(tags), key=len, reverse=True)
    return tuple(_table_suffixes_for_tag(tag) for tag in tags)


def local_blast_table_basename(base_name: str, program: str, *, label_suffix: str = "") -> str:
    """Hit table filename (with extension) next to the query FASTA."""
    tag = local_blast_run_tag(program)
    return f"{base_name}{label_suffix}_{tag}.txt"


def online_blast_table_basename(query_basename: str, program: str) -> str:
    tag = online_blast_run_tag(program)
    return f"{query_basename}_{tag}.txt"


def online_blast_html_basename(query_basename: str, program: str) -> str:
    tag = online_blast_run_tag(program)
    return f"{query_basename}_{tag}.html"


def blast_table_source_tag(filename: str) -> str:
    """Short tag for TE plot suffixes (``bn``, ``b``, ``RM``, …)."""
    name = filename.lower()
    if name.endswith(".out"):
        return "RM"
    stripped = strip_blast_table_suffix(filename)
    if stripped is None:
        return "b"
    rest = filename[len(stripped) :]
    if rest.lower().endswith(".txt"):
        rest = rest[:-4]
    return rest.lstrip("_") or "b"


def strip_blast_table_suffix(name: str) -> str | None:
    """If ``name`` is a BLAST table, return the query basename; else ``None``."""
    lower = name.lower()
    for suffix in all_blast_table_txt_suffixes():
        if lower.endswith(suffix):
            return name[: -len(suffix)]
    return None


def blast_table_candidate_names(query_basename: str, *, stem_only: bool = False) -> list[str]:
    """Possible hit-table filenames for a query (canonical tags first)."""
    names: list[str] = []
    seen: set[str] = set()

    def _add(name: str) -> None:
        if name not in seen:
            seen.add(name)
            names.append(name)

    for tag in LOCAL_BLAST_PROGRAM_TAGS.values():
        _add(f"{query_basename}_{tag}.txt")
    for tag in _LEGACY_LOCAL_BLAST_TAGS.values():
        _add(f"{query_basename}_{tag}.txt")
    _add(f"{query_basename}_blast.txt")
    for tag in ONLINE_BLAST_PROGRAM_TAGS.values():
        _add(f"{query_basename}_{tag}.txt")
    _add(f"{query_basename}_{_LEGACY_ONLINE_TAG}.txt")
    _add(f"{query_basename}_online.txt")

    if stem_only and "." in query_basename:
        stem = Path(query_basename).stem
        for tag in (*LOCAL_BLAST_PROGRAM_TAGS.values(), *_LEGACY_LOCAL_BLAST_TAGS.values()):
            _add(f"{stem}_{tag}.txt")
        _add(f"{stem}_blast.txt")
        for tag in ONLINE_BLAST_PROGRAM_TAGS.values():
            _add(f"{stem}_{tag}.txt")
        _add(f"{stem}_{_LEGACY_ONLINE_TAG}.txt")
        _add(f"{stem}_online.txt")
    return names


def local_blast_txt_path(query_fasta: Path, program: str = "blastn") -> Path:
    query_fasta = query_fasta.expanduser().resolve()
    return query_fasta.parent / local_blast_table_basename(query_fasta.name, program)


def online_blast_txt_path(query_fasta: Path, program: str = "blastn") -> Path:
    query_fasta = query_fasta.expanduser().resolve()
    return query_fasta.parent / online_blast_table_basename(query_fasta.name, program)


def online_blast_html_path(query_fasta: Path, program: str = "blastn") -> Path:
    query_fasta = query_fasta.expanduser().resolve()
    return query_fasta.parent / online_blast_html_basename(query_fasta.name, program)


def local_blast_gff_suffixes() -> tuple[str, ...]:
    tags = list(LOCAL_BLAST_PROGRAM_TAGS.values()) + list(
        _LEGACY_LOCAL_BLAST_TAGS.values()
    )
    tags.append("blast")
    return tuple(f"_{tag}.gff" for tag in sorted(set(tags), key=len, reverse=True))


def local_blast_hit_fasta_suffixes() -> tuple[str, ...]:
    return tuple(f"{s}.fa" for s in all_blast_table_txt_suffixes())


def is_blast_table_filename(name: str) -> bool:
    return strip_blast_table_suffix(name) is not None


def online_blast_html_candidate_names(query_basename: str) -> list[str]:
    names: list[str] = []
    seen: set[str] = set()

    def _add(name: str) -> None:
        if name not in seen:
            seen.add(name)
            names.append(name)

    for tag in ONLINE_BLAST_PROGRAM_TAGS.values():
        _add(f"{query_basename}_{tag}.html")
    _add(f"{query_basename}_{_LEGACY_ONLINE_TAG}.html")
    _add(f"{query_basename}_online.html")
    return names


def online_blast_table_candidate_names(query_basename: str) -> list[str]:
    names: list[str] = []
    seen: set[str] = set()

    def _add(name: str) -> None:
        if name not in seen:
            seen.add(name)
            names.append(name)

    for tag in ONLINE_BLAST_PROGRAM_TAGS.values():
        _add(f"{query_basename}_{tag}.txt")
    _add(f"{query_basename}_{_LEGACY_ONLINE_TAG}.txt")
    _add(f"{query_basename}_online.txt")
    return names


def is_online_blast_html_filename(name: str) -> bool:
    lower = name.lower()
    if not lower.endswith(".html"):
        return False
    for tag in (*ONLINE_BLAST_PROGRAM_TAGS.values(), _LEGACY_ONLINE_TAG, "online"):
        if lower.endswith(f"_{tag}.html"):
            return True
    return False


def is_online_blast_table_filename(name: str) -> bool:
    lower = name.lower()
    if not lower.endswith(".txt"):
        return False
    for tag in (*ONLINE_BLAST_PROGRAM_TAGS.values(), _LEGACY_ONLINE_TAG, "online"):
        if lower.endswith(f"_{tag}.txt"):
            return True
    return False


def find_online_blast_table_path(query_fasta: Path) -> Path | None:
    query_fasta = query_fasta.expanduser().resolve()
    for name in online_blast_table_candidate_names(query_fasta.name):
        path = query_fasta.parent / name
        if path.is_file():
            return path
    return None


def find_online_blast_html_path(query_fasta: Path) -> Path | None:
    query_fasta = query_fasta.expanduser().resolve()
    for name in online_blast_html_candidate_names(query_fasta.name):
        path = query_fasta.parent / name
        if path.is_file():
            return path
    return None
