"""NCBI genetic code tables for sequence translation (transl_table).

Tables follow https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi
and are resolved via Biopython's CodonTable (same IDs as GenBank).
"""

from __future__ import annotations

from Bio.Data import CodonTable
from Bio.Data.IUPACData import protein_letters_1to3

# One-letter → three-letter (NCBI-style Ter for stop).
_AA_1TO3: dict[str, str] = {
    **{k: v for k, v in protein_letters_1to3.items()},
    "*": "Ter",
    "X": "Xaa",
    "B": "Asx",
    "Z": "Glx",
    "J": "Xle",
    "U": "Sec",
    "O": "Pyl",
}


def genetic_code_choices() -> list[tuple[int, str]]:
    """Return ``(table_id, label)`` for every unambiguous DNA genetic code."""
    out: list[tuple[int, str]] = []
    for table_id in sorted(CodonTable.unambiguous_dna_by_id):
        table = CodonTable.unambiguous_dna_by_id[table_id]
        name = table.names[0] if table.names else f"Table {table_id}"
        out.append((int(table_id), f"{table_id}. {name}"))
    return out


def normalize_genetic_code_id(table_id: object, *, default: int = 1) -> int:
    try:
        value = int(table_id)  # type: ignore[arg-type]
    except (TypeError, ValueError):
        return default
    if value in CodonTable.unambiguous_dna_by_id:
        return value
    return default


def amino_acid_display(residue: str, *, three_letter: bool) -> str:
    """Format a translated residue for the sequence track."""
    letter = (residue or "X")[:1].upper()
    if not three_letter:
        return letter
    return _AA_1TO3.get(letter, letter if letter else "Xaa")


def translate_dna_codons(dna: str, *, table: int = 1) -> str:
    """Translate a codon-aligned DNA fragment with an NCBI genetic code table."""
    from Bio.Seq import Seq

    text = (dna or "").upper().replace("U", "T")
    rem = len(text) % 3
    if rem:
        text = text[: len(text) - rem]
    if not text:
        return ""
    table_id = normalize_genetic_code_id(table)
    try:
        return str(Seq(text).translate(table=table_id, to_stop=False))
    except Exception:  # noqa: BLE001
        return ""
