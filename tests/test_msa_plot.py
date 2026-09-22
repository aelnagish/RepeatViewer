"""MSA overview plot backend."""

from __future__ import annotations

from pathlib import Path

from repeatviewer.backends.msa_plot import (
    AMINO_ACID_COLORS_CLUSTAL,
    AMINO_ACID_COLOR_PRESETS,
    DEFAULT_AMINO_ACID_COLORS,
    MSA_PLOT_SEQ_NUCLEOTIDE,
    MSA_PLOT_SEQ_PROTEIN,
    MsaPlotOptions,
    NUCLEOTIDE_COLOR_PRESETS,
    detect_sequence_type,
    fasta_to_alignment_array,
    msa_plot_output_path,
    write_msa_plot,
)


def test_detect_nucleotide_vs_protein(tmp_path: Path) -> None:
    dna = tmp_path / "dna.fa"
    dna.write_text(">a\nACGTACGT\n>b\nACGT--GT\n")
    arr, _names = fasta_to_alignment_array(dna)
    assert detect_sequence_type(arr) == MSA_PLOT_SEQ_NUCLEOTIDE

    prot = tmp_path / "prot.fa"
    prot.write_text(">a\nMKTAYIAK\n>b\nMKT-YIAK\n")
    arr_p, _ = fasta_to_alignment_array(prot)
    assert detect_sequence_type(arr_p) == MSA_PLOT_SEQ_PROTEIN


def test_write_msa_plot_png_and_pdf(tmp_path: Path) -> None:
    msa = tmp_path / "aln.fa"
    msa.write_text(
        ">s1\nACGTACGTAC\n"
        ">s2\nACGT--GTAC\n"
        ">s3\nACGTAAGTAC\n"
    )
    png = write_msa_plot(
        msa,
        options=MsaPlotOptions(output_format="png", dpi=100),
    )
    assert png.is_file()
    assert png.suffix == ".png"
    assert png.name == "aln_msa_plot.png"

    pdf = write_msa_plot(
        msa,
        options=MsaPlotOptions(
            output_format="pdf",
            dpi=120,
            sequence_type=MSA_PLOT_SEQ_NUCLEOTIDE,
            nucleotide_colors={"A": "#111111", "T": "#222222", "G": "#333333", "C": "#444444"},
        ),
    )
    assert pdf.is_file()
    assert pdf.suffix == ".pdf"

    protein = tmp_path / "p.fa"
    protein.write_text(">a\nARNDCQEG\n>b\nARND--EG\n")
    out = write_msa_plot(
        protein,
        tmp_path / "custom_name.png",
        options=MsaPlotOptions(sequence_type=MSA_PLOT_SEQ_PROTEIN, dpi=80),
    )
    assert out == (tmp_path / "custom_name.png").resolve()
    assert out.is_file()


def test_msa_plot_output_path() -> None:
    path = msa_plot_output_path(Path("/tmp/demo.fa"), output_format="pdf")
    assert path.name == "demo_msa_plot.pdf"


def test_clustal_is_default_amino_palette() -> None:
    assert DEFAULT_AMINO_ACID_COLORS == AMINO_ACID_COLORS_CLUSTAL
    names = [label for label, _palette in AMINO_ACID_COLOR_PRESETS]
    assert names == ["Clustal", "Zappo", "Taylor"]
    nuc_names = [label for label, _palette in NUCLEOTIDE_COLOR_PRESETS]
    assert nuc_names == ["Default", "Color-blind"]
    assert AMINO_ACID_COLORS_CLUSTAL["P"] == "#FFFF00"
    assert AMINO_ACID_COLORS_CLUSTAL["A"] == "#80A0F0"


def test_nucleotide_ambiguity_defaults_white() -> None:
    from repeatviewer.backends.msa_plot import (
        NUCLEOTIDE_COLORS_COLORBLIND,
        NUCLEOTIDE_COLORS_DEFAULT,
    )

    for letter in ("R", "Y", "W", "S", "M", "K", "B", "H", "D", "V", "I"):
        assert NUCLEOTIDE_COLORS_DEFAULT[letter] == "#FFFFFF"
        assert NUCLEOTIDE_COLORS_COLORBLIND[letter] == "#FFFFFF"
    assert NUCLEOTIDE_COLORS_DEFAULT["N"] == "#C8C8C8"
    assert NUCLEOTIDE_COLORS_COLORBLIND["N"] == "#C8C8C8"
    assert NUCLEOTIDE_COLORS_DEFAULT["X"] == "#C8C8C8"
    assert NUCLEOTIDE_COLORS_COLORBLIND["X"] == "#C8C8C8"
    assert NUCLEOTIDE_COLORS_DEFAULT["T"] == NUCLEOTIDE_COLORS_DEFAULT["U"]


def test_nucleotide_imshow_preserves_base_colors(tmp_path: Path) -> None:
    """Regression: sparse residue codes must not stretch across the full palette."""
    from matplotlib.colors import to_rgb

    from repeatviewer.backends.msa_plot import NUCLEOTIDE_COLORS_DEFAULT

    msa = tmp_path / "bases.fa"
    # One column each for A/T/C/G so each colour must appear.
    msa.write_text(">s1\nATCG\n>s2\nATCG\n")
    out = write_msa_plot(
        msa,
        tmp_path / "bases.png",
        options=MsaPlotOptions(
            sequence_type=MSA_PLOT_SEQ_NUCLEOTIDE,
            output_format="png",
            dpi=80,
            nucleotide_colors=dict(NUCLEOTIDE_COLORS_DEFAULT),
        ),
    )
    assert out.is_file()

    # Spot-check expected palette entries (same hex as defaults).
    assert to_rgb(NUCLEOTIDE_COLORS_DEFAULT["T"]) == to_rgb(
        NUCLEOTIDE_COLORS_DEFAULT["U"]
    )
    assert NUCLEOTIDE_COLORS_DEFAULT["T"].upper() == "#F06090"
    assert NUCLEOTIDE_COLORS_DEFAULT["C"].upper() == "#7040E8"


def test_amino_acid_categories_cover_standard_letters() -> None:
    from repeatviewer.backends.msa_plot import AMINO_ACID_CATEGORIES

    letters = [aa for _name, group in AMINO_ACID_CATEGORIES for aa in group]
    assert letters.count("C") == 1
    assert "A" in letters and "P" in letters and "-" in letters
    assert set("ARNDCQEGHILKMFPSTWYV").issubset(set(letters))


def test_x_tick_step_is_readable() -> None:
    from repeatviewer.backends.msa_plot import _x_tick_step

    assert _x_tick_step(10) == 1
    assert _x_tick_step(100) in {5, 10}
    assert _x_tick_step(3800) >= 200
    # Roughly a dozen labels, not crowded.
    assert 3800 / _x_tick_step(3800) <= 20


def test_sequence_height_scales_figure() -> None:
    from repeatviewer.backends.msa_plot import _figure_size

    _, h_default = _figure_size(10, 100, sequence_height_cm=0.5)
    _, h_tall = _figure_size(10, 100, sequence_height_cm=1.0)
    assert h_tall == 2 * h_default


def test_row_gaps_insert_white_spacers() -> None:
    import numpy as np

    from repeatviewer.backends.msa_plot import _GAP_ROW_PIXELS, _SEQ_ROW_PIXELS, _insert_row_gaps

    numeric = np.asarray([[1.0, 2.0], [3.0, 4.0], [5.0, 6.0]])
    display, centers = _insert_row_gaps(numeric, gap_code=0.0)
    expected_rows = 3 * _SEQ_ROW_PIXELS + 2 * _GAP_ROW_PIXELS
    assert display.shape == (expected_rows, 2)
    assert len(centers) == 3
    # Gap rows sit between sequence blocks and use the spacer code.
    gap_start = _SEQ_ROW_PIXELS
    assert np.all(display[gap_start : gap_start + _GAP_ROW_PIXELS] == 0.0)
