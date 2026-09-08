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

import re
from dataclasses import dataclass
from pathlib import Path

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


def group_sequences(
    files: list[Path],
    out_path: Path,
    *,
    prefix_with_filename: bool = False,
) -> tuple[Path, int]:
    """Concatenate FASTA files into one multi-record file.

    Duplicate record ids are disambiguated rather than silently written twice,
    since most downstream tools (including BLAST and tree building) require
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
) -> list[Path]:
    """Pull chosen records out of a multi-FASTA.

    Unlike "Separate FASTA", which splits everything into its own directory,
    this takes only the records asked for and writes them beside the source --
    either one file each, or all of them into one new multi-FASTA.
    """
    if not wanted:
        raise SequenceComposeError("Select at least one sequence to extract.")

    try:
        records = [
            r for r in SeqIO.parse(str(path), "fasta") if str(r.id) in wanted
        ]
    except (OSError, ValueError) as exc:
        raise SequenceComposeError(f"Could not read {path.name}: {exc}") from exc
    if not records:
        raise SequenceComposeError("None of the selected sequences were found.")

    out_dir.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    if single_file:
        stem = safe_file_stem(combined_name or f"{path.stem}_extracted")
        out_path = _unique_path(out_dir, stem)
        with out_path.open("w", encoding="utf-8") as handle:
            for record in records:
                handle.write(f">{record.description or record.id}\n")
                handle.write(_wrap(str(record.seq)) + "\n")
        return [out_path]

    for record in records:
        out_path = _unique_path(out_dir, safe_file_stem(str(record.id)))
        with out_path.open("w", encoding="utf-8") as handle:
            handle.write(f">{record.description or record.id}\n")
            handle.write(_wrap(str(record.seq)) + "\n")
        written.append(out_path)
    return written


def default_group_name(files: list[Path]) -> str:
    """A sensible name for the grouped file, based on the selection."""
    if not files:
        return "grouped_sequences"
    if len(files) == 2:
        return f"{files[0].stem}_and_{files[1].stem}"
    return f"{files[0].stem}_and_{len(files) - 1}_more"
