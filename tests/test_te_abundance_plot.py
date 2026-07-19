from pathlib import Path

import pandas as pd

from teviewer.backends.te_abundance_plot import (
    TE_ABUNDANCE_PERCENT,
    _merged_interval_bp,
    family_occupied_bp,
    write_te_abundance_plot_html,
)


def test_repeatmasker_query_coordinates_are_inclusive() -> None:
    assert _merged_interval_bp([(97218, 97398)]) == 181


def test_family_coverage_merges_only_within_same_strand() -> None:
    hits = pd.DataFrame(
        [
            {"family": "A", "chrom": "chr1", "sstrand": "+", "start": 1, "end": 10},
            {"family": "A", "chrom": "chr1", "sstrand": "+", "start": 5, "end": 15},
            {"family": "A", "chrom": "chr1", "sstrand": "-", "start": 5, "end": 15},
            {"family": "B", "chrom": "chr1", "sstrand": "+", "start": 20, "end": 24},
        ]
    )

    abundance = family_occupied_bp(hits)
    assert abundance["A"] == 26  # 15 bp forward + 11 bp reverse
    assert abundance["B"] == 5


def test_abundance_html_defaults_to_ranked_percentage_curve(tmp_path: Path) -> None:
    source = tmp_path / "sample.out"
    source.write_text(
        "   SW  perc perc perc  query       position in query     matching repeat\n"
        "score  div. del. ins.  sequence     begin end (left)   repeat class/family\n"
        "  100  3.0  0.0  0.0  chr1 1 10 (90) + Family_A LINE/L1 1 10 (0) 1\n"
        "  100  3.0  0.0  0.0  chr1 5 15 (85) C Family_A LINE/L1 1 11 (0) 2\n"
        "  100  3.0  0.0  0.0  chr1 30 34 (66) + Family_B DNA/hAT 1 5 (0) 3\n",
        encoding="utf-8",
    )

    output = write_te_abundance_plot_html(source, metric=TE_ABUNDANCE_PERCENT, top_n=100)
    page = output.read_text(encoding="utf-8")
    assert "Cumulative occupied genome (%)" in page
    assert "Family_A" in page
    assert "Family_B" in page
    assert "TE families (ranked)" not in page
    assert "LINE/L1" not in page
    assert "Ranked TE families" in page
    assert "Occupied bp" in page
    assert "Genome %" in page
    assert "overflow-wrap:anywhere" in page
    assert "text-overflow:ellipsis" not in page
    assert 'id="te-list-drag-handle"' in page
    assert 'id="te-list-toggle"' in page
    assert "resize:both" in page
    assert ">21<" in page  # Family_A: 10 bp forward + 11 bp reverse.
    assert ">21.00<" in page  # 21 occupied bp out of a 100 bp query.
    assert "%{customdata[2]:.2f}%" in page
    assert '"plot_bgcolor":"#ffffff"' in page
    assert '"showgrid":false' in page
    assert '"showline":true' in page


def test_multiple_repeatmasker_outputs_get_distinct_curves(tmp_path: Path) -> None:
    header = (
        "   SW  perc perc perc  query       position in query     matching repeat\n"
        "score  div. del. ins.  sequence     begin end (left)   repeat class/family\n"
    )
    first = tmp_path / "genome_a.out"
    second = tmp_path / "genome_b.out"
    first.write_text(
        header + " 100 3.0 0.0 0.0 chr1 1 10 (90) + Family_A LINE/L1 1 10 (0) 1\n",
        encoding="utf-8",
    )
    second.write_text(
        header + " 100 3.0 0.0 0.0 chr1 1 20 (80) + Family_B DNA/hAT 1 20 (0) 1\n",
        encoding="utf-8",
    )

    output = write_te_abundance_plot_html([first, second])
    page = output.read_text(encoding="utf-8")
    assert output.name == "genome_a_comparison_te_abundance_plot.html"
    assert "genome_a" in page
    assert "genome_b" in page
    assert "#3b82f6" in page
    assert "#ef4444" in page
    assert "RepeatMasker outputs" in page
    assert 'id="te-sf-save-svg"' in page
    assert 'id="te-sf-save-pdf"' in page
