"""Annotation sequence molecule detection and amino-acid colors."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.sequence_features import detect_sequence_is_protein
from repeatviewer.widgets.preview_zoom import (
    amino_acid_color_map,
    default_amino_acid_backgrounds,
    load_custom_amino_acid_backgrounds,
    save_custom_amino_acid_backgrounds,
)


def test_detect_protein_from_faa_suffix(tmp_path: Path) -> None:
    path = tmp_path / "orf.faa"
    # Nucleotide-like letters only — suffix still forces protein.
    assert detect_sequence_is_protein("ACGTACGT", path=path) is True


def test_detect_protein_from_composition() -> None:
    protein = "MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVVHS"
    dna = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
    assert detect_sequence_is_protein(protein) is True
    assert detect_sequence_is_protein(dna) is False


def test_amino_acid_color_map_defaults_to_clustal() -> None:
    from repeatviewer.backends.msa_plot import AMINO_ACID_COLORS_CLUSTAL

    save_custom_amino_acid_backgrounds({})
    defaults = default_amino_acid_backgrounds()
    assert defaults["A"] == AMINO_ACID_COLORS_CLUSTAL["A"]
    palette = amino_acid_color_map()
    assert palette["A"][0] == AMINO_ACID_COLORS_CLUSTAL["A"]
    assert palette["A"][1] in {"#111111", "#ffffff"}


def test_custom_amino_acid_backgrounds_roundtrip() -> None:
    save_custom_amino_acid_backgrounds({"A": "#123456", "C": "#abcdef"})
    try:
        loaded = load_custom_amino_acid_backgrounds()
        assert loaded["A"].upper() == "#123456"
        assert loaded["C"].upper() == "#ABCDEF"
        palette = amino_acid_color_map()
        assert palette["A"][0].upper() == "#123456"
    finally:
        save_custom_amino_acid_backgrounds({})
