"""Create a new sequence, or group several FASTA files into one.

Two small operations that mirror Geneious's Sequence menu:

* **New sequence** -- type or paste bases and save them as a FASTA. Pasted
  text usually carries line numbers, spaces or FASTA headers from wherever it
  came from, so it is cleaned before being written.
* **Group sequences into a list** -- concatenate several FASTA files into one
  multi-record file, which is the counterpart of the existing "Separate FASTA"
  action (Geneious calls that "extract sequences from list").
"""

from __future__ import annotations

import gzip
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable, Iterable

from Bio import SeqIO

MOLECULE_DNA = "DNA"
MOLECULE_RNA = "RNA"
MOLECULE_PROTEIN = "Protein"

MOLECULE_LABELS = (
    ("Nucleotide (DNA)", MOLECULE_DNA),
    ("Nucleotide (RNA)", MOLECULE_RNA),
    ("Protein", MOLECULE_PROTEIN),
)

_DNA_OK = set("ACGTUNRYKMSWBDHVacgtunrykmswbdhv-*.")
_PROTEIN_OK = set("ACDEFGHIKLMNPQRSTVWYXBZJUO*-.acdefghiklmnpqrstvwyxbzjuo")

_FASTA_LINE_WIDTH = 60


class SequenceComposeError(RuntimeError):
    """The sequence could not be created or the files could not be grouped."""


@dataclass(frozen=True)
class SequenceCheck:
    cleaned: str
    length: int
    invalid: str  # characters that do not belong to the chosen molecule type


def clean_pasted_sequence(text: str) -> str:
    """Strip everything that is not sequence from pasted text.

    Handles the usual cases: FASTA headers, digits from numbered alignment
    output, whitespace and line breaks. Gaps and ambiguity codes are kept --
    they may be meaningful, and are reported as invalid only if they do not
    suit the chosen molecule type.
    """
    lines = []
    for raw in (text or "").splitlines():
        line = raw.strip()
        if not line or line.startswith((">", ";")):
            continue  # FASTA header or comment
        lines.append(line)
    joined = "".join(lines)
    # Digits and whitespace appear in numbered/formatted sequence output.
    return re.sub(r"[\s\d]", "", joined)


def parse_pasted_header(text: str) -> tuple[str, str]:
    """(name, description) from the first FASTA header in pasted text.

    Pasting a record usually brings its ``>`` line along, which carries the
    name the sequence already has -- better to offer it than make the user
    retype it. The first whitespace-separated token becomes the name and the
    remainder the description, matching how FASTA is normally read.
    """
    for raw in (text or "").splitlines():
        line = raw.strip()
        if not line.startswith(">"):
            continue
        header = line[1:].strip()
        if not header:
            return "", ""
        parts = header.split(None, 1)
        name = parts[0]
        description = parts[1].strip() if len(parts) > 1 else ""
        return name, description
    return "", ""


def check_sequence(text: str, molecule: str = MOLECULE_DNA) -> SequenceCheck:
    """Clean the text and report any characters that do not fit ``molecule``."""
    cleaned = clean_pasted_sequence(text)
    allowed = _PROTEIN_OK if molecule == MOLECULE_PROTEIN else _DNA_OK
    invalid = sorted({ch for ch in cleaned if ch not in allowed})
    return SequenceCheck(
        cleaned=cleaned, length=len(cleaned), invalid="".join(invalid)
    )


def _wrap(sequence: str, width: int = _FASTA_LINE_WIDTH) -> str:
    return "\n".join(
        sequence[i : i + width] for i in range(0, len(sequence), width)
    ) or ""


def _unique_path(folder: Path, stem: str, suffix: str = ".fasta") -> Path:
    """A path that does not exist yet, adding _1, _2 ... if needed."""
    candidate = folder / f"{stem}{suffix}"
    index = 1
    while candidate.exists():
        candidate = folder / f"{stem}_{index}{suffix}"
        index += 1
    return candidate


def safe_file_stem(name: str) -> str:
    """Filesystem-safe stem; the record keeps the name as typed."""
    cleaned = "".join(ch if ch.isalnum() or ch in "-_." else "_" for ch in name)
    return cleaned.strip("._") or "sequence"


def extract_source_suffix(source: Path) -> str:
    """FASTA suffix for extracted files (always ``.fa``, or ``.faa`` for protein)."""
    name = source.name.lower()
    if name.endswith(".faa") or name.endswith(".faa.gz") or name.endswith(".faa.bgz"):
        return ".faa"
    return ".fa"


def proposed_extract_path(
    source: Path,
    record_ids: list[str],
    out_dir: Path | None = None,
    *,
    single_file: bool = True,
) -> Path:
    """Proposed output path for an extract (may already exist).

    - One sequence → ``{sanitize_name(id)}{source_suffix}`` (Separate Sequences rule).
    - Several sequences into one file → ``{source.name}_extr.fa``.
    """
    from repeatviewer.backends.separate_fasta import sanitize_name

    source = source.expanduser().resolve()
    folder = (out_dir or source.parent).expanduser().resolve()
    ids = [str(i) for i in record_ids if str(i).strip()]
    if single_file and len(ids) == 1:
        return folder / f"{sanitize_name(ids[0])}{extract_source_suffix(source)}"
    if single_file:
        return folder / f"{source.name}_extr.fa"
    # One-file-per-sequence caller should use proposed_single_record_path.
    return folder / f"{source.name}_extr.fa"


def msa_extract_name_stem(source: Path, record_id: str) -> str:
    """``{original file stem}_{sequence name}`` for an MSA extract.

    The file stem keeps characters such as ``#`` that appear in the alignment
    name. The sequence id is sanitized for the filesystem.
    """
    from repeatviewer.backends.separate_fasta import sanitize_name

    file_part = Path(source).stem
    for char, replacement in (
        ("/", "__"),
        ("\\", "__"),
        (" ", "_"),
        ("|", "_"),
        ("<", "_"),
        (":", "_"),
        ('"', "_"),
        ("*", "_"),
        ("?", "_"),
    ):
        file_part = file_part.replace(char, replacement)
    file_part = file_part.strip("._") or "sequence"
    seq_part = sanitize_name(str(record_id)).strip("._") or "sequence"
    return f"{file_part}_{seq_part}"


def proposed_msa_extract_path(
    source: Path,
    record_id: str,
    out_dir: Path | None = None,
) -> Path:
    """MSA display extract: ``{original_file_name}_{sequence_name}{suffix}``.

    The multi-sequence (unaligned) overview keeps :func:`proposed_extract_path`.
    """
    source = Path(source).expanduser().resolve()
    folder = (out_dir or source.parent).expanduser().resolve()
    return folder / (
        f"{msa_extract_name_stem(source, record_id)}{extract_source_suffix(source)}"
    )


def proposed_single_record_path(source: Path, record_id: str, out_dir: Path | None = None) -> Path:
    """Path for one extracted sequence, using Separate Sequences sanitization."""
    from repeatviewer.backends.separate_fasta import sanitize_name

    source = source.expanduser().resolve()
    folder = (out_dir or source.parent).expanduser().resolve()
    return folder / f"{sanitize_name(str(record_id))}{extract_source_suffix(source)}"


def write_new_sequence(
    folder: Path,
    *,
    name: str,
    sequence: str,
    description: str = "",
    organism: str = "",
    molecule: str = MOLECULE_DNA,
    overwrite_path: Path | None = None,
) -> Path:
    """Write a single-record FASTA and return its path."""
    check = check_sequence(sequence, molecule)
    if not check.cleaned:
        raise SequenceComposeError("The sequence is empty.")
    if not name.strip():
        raise SequenceComposeError("Enter a name for the sequence.")

    header_bits = [name.strip()]
    if description.strip():
        header_bits.append(description.strip())
    if organism.strip():
        header_bits.append(f"[organism={organism.strip()}]")
    header = " ".join(header_bits)

    out_path = overwrite_path or _unique_path(folder, safe_file_stem(name))
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8") as handle:
            handle.write(f">{header}\n")
            handle.write(_wrap(check.cleaned) + "\n")
    except OSError as exc:
        raise SequenceComposeError(f"Could not write {out_path.name}: {exc}") from exc
    return out_path


def combine_annotation_json(sources: list[Path], dest_fasta: Path) -> None:
    """Merge JSON annotations from source FASTAs into the ``dest_fasta`` sidecar.

    Records are matched by order, so ``dest_fasta`` must hold the source
    records in the same order (as written by :func:`group_sequences`).
    """
    from repeatviewer.backends.sequence_features import (
        fasta_record_lengths,
        load_features_by_id_from_json,
        save_annotations_multi,
    )

    dest_lengths = fasta_record_lengths(dest_fasta)
    dest_ids = list(dest_lengths.keys())
    merged: dict[str, list] = {rid: [] for rid in dest_ids}
    dest_index = 0
    for src in sources:
        src_lengths = fasta_record_lengths(src)
        src_ids = list(src_lengths.keys())
        src_feats = load_features_by_id_from_json(src, src_ids)
        for sid in src_ids:
            if dest_index >= len(dest_ids):
                break
            feats = src_feats.get(sid) or []
            if feats:
                merged[dest_ids[dest_index]].extend(feats)
            dest_index += 1
    if any(merged.values()):
        save_annotations_multi(
            dest_fasta,
            features_by_id=merged,
            lengths_by_id=dest_lengths,
            record_ids=dest_ids,
        )


def group_sequences(
    files: list[Path],
    out_path: Path,
    *,
    prefix_with_filename: bool = False,
) -> tuple[Path, int]:
    """Concatenate FASTA files into one multi-record file.

    Duplicate record ids are disambiguated rather than silently written twice,
    since most downstream tools (including BLAST) require
    unique names.
    """
    if len(files) < 2:
        raise SequenceComposeError("Select at least two FASTA files to group.")

    seen: dict[str, int] = {}
    written = 0
    try:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8") as handle:
            for path in files:
                try:
                    records = list(SeqIO.parse(str(path), "fasta"))
                except (OSError, ValueError) as exc:
                    raise SequenceComposeError(
                        f"Could not read {path.name}: {exc}"
                    ) from exc
                if not records:
                    continue
                for record in records:
                    base = str(record.id) or path.stem
                    if prefix_with_filename:
                        base = f"{path.stem}_{base}"
                    if base in seen:
                        seen[base] += 1
                        name = f"{base}_{seen[base]}"
                    else:
                        seen[base] = 0
                        name = base
                    description = str(record.description or "")
                    extra = description[len(str(record.id)):].strip()
                    header = f"{name} {extra}".strip()
                    handle.write(f">{header}\n")
                    handle.write(_wrap(str(record.seq)) + "\n")
                    written += 1
    except SequenceComposeError:
        raise
    except OSError as exc:
        raise SequenceComposeError(f"Could not write {out_path.name}: {exc}") from exc

    if written == 0:
        raise SequenceComposeError("No sequences were found in the selected files.")
    return out_path, written


def list_records(path: Path) -> list[tuple[str, str, int]]:
    """(id, description, length) for every record, for a pick list.

    Reads ids and lengths without holding the sequences, so a large
    multi-FASTA can be listed cheaply.
    """
    entries: list[tuple[str, str, int]] = []
    try:
        for record in SeqIO.parse(str(path), "fasta"):
            rec_id = str(record.id)
            description = str(record.description or "")
            extra = description[len(rec_id):].strip()
            entries.append((rec_id, extra, len(record.seq)))
    except (OSError, ValueError) as exc:
        raise SequenceComposeError(f"Could not read {path.name}: {exc}") from exc
    if not entries:
        raise SequenceComposeError(f"No sequences found in {path.name}.")
    return entries


def extract_records(
    path: Path,
    wanted: set[str],
    out_dir: Path,
    *,
    single_file: bool = False,
    combined_name: str = "",
    output_path: Path | None = None,
    output_paths_by_id: dict[str, Path] | None = None,
) -> list[Path]:
    """Pull chosen records out of a multi-FASTA.

    Unlike "Separate FASTA", which splits everything into its own directory,
    this takes only the records asked for and writes them beside the source --
    either one file each, or all of them into one new multi-FASTA.

    When ``output_path`` is set (combined / single-record extract), that path is
    written directly (caller handles overwrite/rename). When
    ``output_paths_by_id`` is set for multi-file extract, those paths are used
    (caller handles conflicts). Otherwise a unique path is chosen automatically.
    """
    from repeatviewer.backends.separate_fasta import sanitize_name

    if not wanted:
        raise SequenceComposeError("Select at least one sequence to extract.")

    path = path.expanduser().resolve()
    try:
        records = [
            r for r in SeqIO.parse(str(path), "fasta") if str(r.id) in wanted
        ]
    except (OSError, ValueError) as exc:
        raise SequenceComposeError(f"Could not read {path.name}: {exc}") from exc
    if not records:
        raise SequenceComposeError("None of the selected sequences were found.")

    out_dir = out_dir.expanduser().resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    suffix = extract_source_suffix(path)
    path_map = {
        str(k): Path(v).expanduser().resolve()
        for k, v in (output_paths_by_id or {}).items()
    }

    if single_file:
        if output_path is not None:
            out_path = output_path.expanduser().resolve()
            out_path.parent.mkdir(parents=True, exist_ok=True)
        elif len(records) == 1:
            proposed = proposed_extract_path(
                path, [str(records[0].id)], out_dir, single_file=True
            )
            out_path = (
                proposed
                if not proposed.exists()
                else _unique_path(out_dir, sanitize_name(str(records[0].id)), suffix)
            )
        else:
            if combined_name:
                stem = safe_file_stem(combined_name)
                out_path = _unique_path(out_dir, stem, ".fa")
            else:
                proposed = proposed_extract_path(
                    path, [str(r.id) for r in records], out_dir, single_file=True
                )
                if proposed.exists():
                    # ``name_extr.fa`` → stem ``name_extr`` when unique-ifying.
                    stem = proposed.name[: -len(".fa")] if proposed.name.endswith(".fa") else proposed.stem
                    out_path = _unique_path(out_dir, stem, ".fa")
                else:
                    out_path = proposed
        with out_path.open("w", encoding="utf-8") as handle:
            for record in records:
                handle.write(f">{record.description or record.id}\n")
                handle.write(_wrap(str(record.seq)) + "\n")
        return [out_path]

    for record in records:
        rid = str(record.id)
        if rid in path_map:
            out_path = path_map[rid]
            out_path.parent.mkdir(parents=True, exist_ok=True)
        elif output_path is not None and len(records) == 1:
            out_path = output_path.expanduser().resolve()
        else:
            stem = sanitize_name(rid)
            proposed = out_dir / f"{stem}{suffix}"
            out_path = (
                proposed
                if not proposed.exists()
                else _unique_path(out_dir, stem, suffix)
            )
        with out_path.open("w", encoding="utf-8") as handle:
            handle.write(f">{record.description or record.id}\n")
            handle.write(_wrap(str(record.seq)) + "\n")
        written.append(out_path)
    return written


def default_group_name(files: list[Path]) -> str:
    """A sensible name for the grouped file, based on the selection."""
    from repeatviewer.tool_output_paths import kept_input_name

    if not files:
        return "grouped_sequences"
    if len(files) == 2:
        return f"{kept_input_name(files[0])}_and_{kept_input_name(files[1])}_comb"
    return f"{kept_input_name(files[0])}_and_{len(files) - 1}_comb"


# ---------------------------------------------------------------------------
# Extract records named in an ID list
#
# The pick-list dialog above works record by record, which suits a handful of
# monomers. Selecting a few thousand reads out of a read set needs the other
# shape: hand the tool a list of names. Everything here streams -- one pass
# over the FASTA, memory proportional to the ID list rather than the file --
# so a multi-GB (optionally gzipped) read set is fine.
# ---------------------------------------------------------------------------

_GZIP_SUFFIXES = {".gz", ".bgz"}

# First column of a RepeatViewer BLAST table; skipped so a hit table can be
# pasted straight in without deleting its header row first.
_TABLE_HEADER_TOKENS = {"qseqid", "sseqid", "query", "subject", "id", "name"}

_ID_PROGRESS_INTERVAL = 2000


@dataclass(frozen=True)
class IdListParse:
    """Result of reading an ID list: the names, plus what was thrown away."""

    ids: list[str] = field(default_factory=list)
    ignored: int = 0  # blank lines, comments, rows without the wanted column
    duplicates: int = 0  # repeated names, collapsed to one


@dataclass(frozen=True)
class ExtractSummary:
    """What an ID-driven extraction did."""

    written: list[Path] = field(default_factory=list)
    matched: int = 0  # records written
    found_ids: int = 0  # distinct requested ids present in the FASTA
    requested_ids: int = 0
    records_read: int = 0
    missing: list[str] = field(default_factory=list)
    missing_path: Path | None = None
    cancelled: bool = False


def _split_fields(line: str) -> list[str]:
    """Tab-separated when tabs are present, otherwise any whitespace."""
    return line.split("\t") if "\t" in line else line.split()


def _clean_id_token(value: str) -> str:
    value = value.strip().strip("'\"").strip()
    if value.startswith(">"):
        value = value[1:].strip()
    if not value:
        return ""
    # BLAST reports the first whitespace token of a FASTA header, and Biopython
    # calls the same token ``record.id`` -- match on that so a list harvested
    # from a hit table needs no massaging.
    return value.split()[0]


def parse_id_lines(lines: Iterable[str], *, column: int = 1) -> IdListParse:
    """Read sequence names out of arbitrary pasted or loaded text.

    Absorbs the usual mess: blank lines, ``#`` comments, a stray ``>``,
    quotes, CRLF line endings, and whole tables where the names sit in some
    column other than the first. A leading BLAST-style header row is dropped.
    """
    column = max(1, int(column))
    seen: dict[str, None] = {}
    ignored = 0
    duplicates = 0
    first = True

    for raw in lines:
        line = raw.strip()
        if not line or line.startswith("#") or line.startswith(";"):
            ignored += 1
            continue
        fields = _split_fields(line)
        if len(fields) < column:
            ignored += 1
            continue
        value = _clean_id_token(fields[column - 1])
        if not value:
            ignored += 1
            continue
        if first and value.lower() in _TABLE_HEADER_TOKENS:
            first = False
            ignored += 1
            continue
        first = False
        if value in seen:
            duplicates += 1
            continue
        seen[value] = None

    return IdListParse(ids=list(seen), ignored=ignored, duplicates=duplicates)


def parse_id_list(text: str, *, column: int = 1) -> IdListParse:
    """``parse_id_lines`` for a block of pasted text."""
    return parse_id_lines((text or "").splitlines(), column=column)


def read_id_file(path: Path, *, column: int = 1) -> IdListParse:
    """``parse_id_lines`` for a file, read line by line (tables can be large)."""
    path = Path(path).expanduser().resolve()
    try:
        with path.open("r", encoding="utf-8", errors="replace") as handle:
            return parse_id_lines(handle, column=column)
    except OSError as exc:
        raise SequenceComposeError(f"Could not read {path.name}: {exc}") from exc


def _open_fasta_text(path: Path):
    """Open a FASTA for reading, transparently handling ``.gz``."""
    if path.suffix.lower() in _GZIP_SUFFIXES:
        return gzip.open(path, "rt", encoding="utf-8", errors="replace")
    return path.open("r", encoding="utf-8", errors="replace")


def _progress_total(path: Path) -> int:
    """Bytes to read, or 0 when it cannot be known (compressed input)."""
    if path.suffix.lower() in _GZIP_SUFFIXES:
        return 0
    try:
        return path.stat().st_size
    except OSError:
        return 0


def _record_id(header: str) -> str:
    header = header.strip()
    if not header:
        return ""
    return header.split()[0]


def _cancelled(cancel: Callable[[], bool] | None) -> bool:
    return bool(cancel is not None and cancel())


def count_matching_ids(
    source: Path,
    ids: Iterable[str],
    *,
    progress: Callable[[int, int], None] | None = None,
    cancel: Callable[[], bool] | None = None,
) -> ExtractSummary:
    """Check an ID list against a FASTA without writing anything.

    Reads headers only, so it is much cheaper than a full extraction -- enough
    to catch the usual mistake of an ID list from one assembly being used
    against the FASTA of another.
    """
    source = Path(source).expanduser().resolve()
    wanted = list(dict.fromkeys(ids))
    if not wanted:
        raise SequenceComposeError("No sequence IDs were given.")

    seen: set[str] = set()
    wanted_set = set(wanted)
    records = 0
    processed = 0
    total = _progress_total(source)

    try:
        with _open_fasta_text(source) as handle:
            for line in handle:
                processed += len(line)
                if not line.startswith(">"):
                    continue
                records += 1
                rec_id = _record_id(line[1:])
                if rec_id in wanted_set:
                    seen.add(rec_id)
                if records % _ID_PROGRESS_INTERVAL == 0:
                    if _cancelled(cancel):
                        return ExtractSummary(
                            requested_ids=len(wanted),
                            records_read=records,
                            cancelled=True,
                        )
                    if progress is not None:
                        progress(processed, total)
    except (OSError, EOFError) as exc:
        raise SequenceComposeError(f"Could not read {source.name}: {exc}") from exc

    return ExtractSummary(
        matched=0,
        found_ids=len(seen),
        requested_ids=len(wanted),
        records_read=records,
        missing=[i for i in wanted if i not in seen],
    )


def missing_ids_path(out_dir: Path, stem: str) -> Path:
    return Path(out_dir) / f"{safe_file_stem(stem)}_missing_ids.txt"


def extract_records_by_ids(
    source: Path,
    ids: Iterable[str],
    out_dir: Path,
    *,
    single_file: bool = True,
    exclude: bool = False,
    combined_name: str = "",
    output_path: Path | None = None,
    output_paths_by_id: dict[str, Path] | None = None,
    write_missing: bool = True,
    progress: Callable[[int, int], None] | None = None,
    cancel: Callable[[], bool] | None = None,
) -> ExtractSummary:
    """Write the records named in ``ids`` (or everything but them).

    One streaming pass, so the FASTA is never held in memory and the input may
    be gzipped. Records whose ids are absent are reported rather than raised
    on, and every record matching a duplicated header is written -- BLAST hit
    FASTAs legitimately repeat ``chrom:start-end(strand)`` headers.
    """
    source = Path(source).expanduser().resolve()
    out_dir = Path(out_dir).expanduser().resolve()
    wanted = list(dict.fromkeys(ids))
    if not wanted:
        raise SequenceComposeError("No sequence IDs were given.")
    wanted_set = set(wanted)

    out_dir.mkdir(parents=True, exist_ok=True)
    from repeatviewer.backends.separate_fasta import sanitize_name

    forced_path = output_path.expanduser().resolve() if output_path is not None else None
    path_map = {
        str(k): Path(v).expanduser().resolve()
        for k, v in (output_paths_by_id or {}).items()
    }
    if combined_name:
        stem = safe_file_stem(combined_name)
        combined_suffix = ".fa"
    elif single_file and len(wanted) == 1 and not exclude:
        stem = sanitize_name(wanted[0])
        combined_suffix = extract_source_suffix(source)
    else:
        # Multi-record combined file: ``{source.name}_extr.fa``.
        stem = f"{source.name}_extr"
        combined_suffix = ".fa"

    written: list[Path] = []
    seen: set[str] = set()
    matched = 0
    records = 0
    processed = 0
    total = _progress_total(source)
    combined_handle = None
    combined_path: Path | None = None

    def _emit(header: str, sequence: str) -> None:
        nonlocal matched, combined_handle, combined_path
        if single_file:
            if combined_handle is None:
                if forced_path is not None:
                    combined_path = forced_path
                    combined_path.parent.mkdir(parents=True, exist_ok=True)
                else:
                    combined_path = _unique_path(out_dir, stem, combined_suffix)
                combined_handle = combined_path.open("w", encoding="utf-8")
            handle = combined_handle
            handle.write(f">{header}\n")
            handle.write(_wrap(sequence) + "\n")
        else:
            rec_id = _record_id(header)
            if rec_id in path_map:
                rec_path = path_map[rec_id]
                rec_path.parent.mkdir(parents=True, exist_ok=True)
            else:
                rec_stem = sanitize_name(rec_id)
                rec_path = _unique_path(
                    out_dir, rec_stem, extract_source_suffix(source)
                )
            with rec_path.open("w", encoding="utf-8") as record_handle:
                record_handle.write(f">{header}\n")
                record_handle.write(_wrap(sequence) + "\n")
            written.append(rec_path)
        matched += 1

    def _discard_partial() -> None:
        nonlocal combined_handle
        if combined_handle is not None:
            combined_handle.close()
            combined_handle = None
        if combined_path is not None:
            combined_path.unlink(missing_ok=True)
        for path in written:
            path.unlink(missing_ok=True)

    try:
        with _open_fasta_text(source) as handle:
            header: str | None = None
            chunks: list[str] = []

            def _flush() -> None:
                nonlocal header, chunks
                if header is None:
                    return
                rec_id = _record_id(header)
                hit = rec_id in wanted_set
                if hit:
                    seen.add(rec_id)
                if hit != exclude:
                    _emit(header, "".join(chunks))
                header = None
                chunks = []

            for line in handle:
                processed += len(line)
                if line.startswith(">"):
                    _flush()
                    records += 1
                    header = line[1:].strip()
                    chunks = []
                    if records % _ID_PROGRESS_INTERVAL == 0:
                        if _cancelled(cancel):
                            _discard_partial()
                            return ExtractSummary(
                                requested_ids=len(wanted),
                                records_read=records,
                                cancelled=True,
                            )
                        if progress is not None:
                            progress(processed, total)
                else:
                    chunks.append(line.strip())
            _flush()
    except (OSError, EOFError) as exc:
        _discard_partial()
        raise SequenceComposeError(f"Could not read {source.name}: {exc}") from exc
    finally:
        if combined_handle is not None:
            combined_handle.close()

    if combined_path is not None:
        written.append(combined_path)

    missing = [i for i in wanted if i not in seen]
    missing_path: Path | None = None
    if missing and write_missing and written:
        missing_path = missing_ids_path(out_dir, stem)
        try:
            missing_path.write_text(
                "\n".join(missing) + "\n", encoding="utf-8"
            )
        except OSError:
            missing_path = None

    return ExtractSummary(
        written=written,
        matched=matched,
        found_ids=len(seen),
        requested_ids=len(wanted),
        records_read=records,
        missing=missing,
        missing_path=missing_path,
    )
