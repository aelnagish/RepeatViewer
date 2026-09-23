"""Split multi-sequence FASTA files into one file per sequence (from TEtrimmerGUI)."""

from __future__ import annotations

import logging
from pathlib import Path

from Bio import SeqIO

logger = logging.getLogger(__name__)

SPECIAL_CHARACTER = {
    "/": "__",
    "\\": "__",
    " ": "_",
    "|": "_",
    "<": "_",
    ":": "_",
    '"': "_",
    "*": "_",
    "?": "_",
    "#": "_",
}


def sanitize_name(name: str) -> str:
    for char, replacement in SPECIAL_CHARACTER.items():
        name = name.replace(char, replacement)
    return name


def count_fasta_sequences(path: Path) -> int:
    path = path.expanduser().resolve()
    if not path.is_file():
        return 0
    return sum(1 for _ in SeqIO.parse(str(path), "fasta"))


def separation_output_dir(input_file: Path) -> Path:
    """
    Choose an output folder for a multi-FASTA split.

    Naming rule:
    - Always based on the input file name (stem), e.g. ``myfile_sep``.
    - If the folder already exists, append ``_2``, ``_3``, ... to avoid conflicts.
    """
    input_file = input_file.expanduser().resolve()
    parent = input_file.parent
    base = f"{sanitize_name(input_file.stem)}_sep"

    candidate = parent / base
    if not candidate.exists():
        return candidate

    index = 2
    while True:
        candidate = parent / f"{base}_{index}"
        if not candidate.exists():
            return candidate
        index += 1


def _mapping_path(input_file: Path, output_dir: Path) -> Path:
    parent = input_file.parent
    return parent / f"Sequence_name_mapping_{sanitize_name(input_file.stem)}.txt"


def separate_sequences(
    input_file: Path,
    output_dir: Path | None = None,
    *,
    write_mapping: bool = True,
    subfolder_per_sequence: bool = False,
) -> tuple[list[Path], int]:
    """
    Write one FASTA file per sequence (TEtrimmerGUI ``separate_sequences`` behaviour).

    When ``subfolder_per_sequence`` is True, each sequence is written to
    ``{output_dir}/{sanitized_id}/{sanitized_id}{ext}`` where ``ext`` matches
    the input suffix (e.g. ``.faa`` for protein FASTA).

    Returns (list of created file paths, number of sequences written).
    """
    input_file = input_file.expanduser().resolve()
    if not input_file.is_file():
        raise FileNotFoundError(input_file)

    if output_dir is None:
        output_dir = separation_output_dir(input_file)
    output_dir = output_dir.expanduser().resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    mapping_path = _mapping_path(input_file, output_dir)
    created: list[Path] = []
    id_list: list[str] = []

    mapping_handle = None
    if write_mapping:
        mapping_handle = mapping_path.open("w", encoding="utf-8")
        mapping_handle.write("original_input_seq_name\tRepeatViewer_modified_seq_name\n")

    try:
        with input_file.open("r", encoding="utf-8", errors="replace") as fasta_file:
            for record in SeqIO.parse(fasta_file, "fasta"):
                original_id = record.id
                sanitized_id = sanitize_name(original_id)
                header_id = original_id

                if sanitized_id in id_list:
                    id_list.append(sanitized_id)
                    count = id_list.count(sanitized_id)
                    sanitized_id = f"{sanitized_id}_{count}"
                    header_id = f"{original_id}_{count}"
                    logger.warning("Duplicate name %s renamed to %s", original_id, sanitized_id)
                else:
                    id_list.append(sanitized_id)

                if mapping_handle is not None:
                    mapping_handle.write(f"{original_id}\t{sanitized_id}\n")

                ext = input_file.suffix.lower() or ".fa"
                if subfolder_per_sequence:
                    seq_dir = output_dir / sanitized_id
                    seq_dir.mkdir(parents=True, exist_ok=True)
                    output_filename = seq_dir / f"{sanitized_id}{ext}"
                else:
                    output_filename = output_dir / f"{sanitized_id}{ext}"
                # Keep biological header naming (including '/') in FASTA output.
                # Only the file name itself is sanitized.
                record.id = header_id
                record.description = ""

                with output_filename.open("w", encoding="utf-8") as output_file:
                    SeqIO.write(record, output_file, "fasta")
                created.append(output_filename)
    finally:
        if mapping_handle is not None:
            mapping_handle.close()

    return created, len(created)
