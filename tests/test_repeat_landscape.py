"""Repeat landscape: CpG-adjusted Kimura from .align, no double counting, bins.

The landscape used to bin the raw ``.out`` "perc div." column and summed hit
lengths, so overlapping hits were counted twice and the axis was not the
Kimura metric published landscapes use.
"""

from pathlib import Path

import pytest

from repeatviewer.backends.repeat_quantify import (
    QuantifyResult,
    _aggregate_rm_rows,
    _bin_start,
    _landscape_from_rows,
    _parse_repeatmasker_out,
    _write_landscape_tsv,
    _write_summary_tsv,
    align_path_for_out,
    apply_kimura_to_rows,
    jukes_cantor_percent,
    parse_repeatmasker_align_kimura,
    write_landscape_html,
    write_summary_html,
)

OUT_TEXT = """\
   SW   perc perc perc  query  position in query      matching repeat    position in repeat
score   div. del. ins.  sequence  begin end (left)   repeat  class/family  begin end (left)  ID

  1200   8.2  0.0  0.0  q0     1   150 (850) +  L0  Unspecified   1 150   (0)  1
  1100  12.0  0.0  0.0  q0   101   250 (750) C  L0  Unspecified (0) 150     1  2
   900  20.0  1.0  0.0  q0   400   549 (451) +  L1  Unspecified   1 150   (0)  3
"""

# Minimal .align in RepeatMasker -a layout (forward and complement headers).
ALIGN_TEXT = """\
1200 8.20 0.00 0.00 q0 1 150 (850) L0#Unspecified 1 150 (0) 1

  q0                     1 ACGTACGTAC 10
                           i
  L0#Unspecifie          1 ACGTACGTAC 10

Matrix = 20p43g.matrix
Kimura (with divCpGMod) = 9.10
Transitions / transversions = 1.40 (7/5)
Gap_init rate = 0.00 (0 / 149), avg. gap size = 0.00 (0 / 0)

1100 12.00 0.00 0.00 q0 101 250 (750) C L0#Unspecified (0) 150 1 2

  q0                   101 ACGTACGTAC 110

C L0#Unspecifie        150 ACGTACGTAC 141

Matrix = 20p43g.matrix
Kimura (with divCpGMod) = 13.55
Transitions / transversions = 1.20 (6/5)
Gap_init rate = 0.00 (0 / 149), avg. gap size = 0.00 (0 / 0)

"""


def _write_run(tmp_path: Path, *, with_align: bool = True) -> Path:
    out = tmp_path / "scaf_rm.fa.out"
    out.write_text(OUT_TEXT)
    if with_align:
        (tmp_path / "scaf_rm.fa.align").write_text(ALIGN_TEXT)
    return out


def test_align_path_found_next_to_out(tmp_path):
    out = _write_run(tmp_path)
    assert align_path_for_out(out) == tmp_path / "scaf_rm.fa.align"


def test_align_path_missing_returns_none(tmp_path):
    assert align_path_for_out(_write_run(tmp_path, with_align=False)) is None


def test_parse_align_forward_and_complement(tmp_path):
    _write_run(tmp_path)
    kimura = parse_repeatmasker_align_kimura(tmp_path / "scaf_rm.fa.align")
    assert kimura == {
        ("q0", 1, 150, "L0"): 9.10,
        ("q0", 101, 250, "L0"): 13.55,
    }


def test_unmatched_hits_fall_back_to_jukes_cantor(tmp_path):
    out = _write_run(tmp_path)
    rows = _parse_repeatmasker_out(out)
    kimura = parse_repeatmasker_align_kimura(tmp_path / "scaf_rm.fa.align")
    unmatched = apply_kimura_to_rows(rows, kimura)
    assert unmatched == 1
    by_start = {r["start"]: r["landscape_div"] for r in rows}
    assert by_start[1] == pytest.approx(9.10)
    assert by_start[101] == pytest.approx(13.55)
    assert by_start[400] == pytest.approx(jukes_cantor_percent(20.0))


def test_jukes_cantor_values():
    assert jukes_cantor_percent(0.0) == 0.0
    assert jukes_cantor_percent(10.0) == pytest.approx(10.73, abs=0.01)
    assert jukes_cantor_percent(30.0) == pytest.approx(38.31, abs=0.01)
    assert jukes_cantor_percent(80.0) == 50.0  # saturated -> cap


def test_overlap_counted_once_and_given_to_least_diverged(tmp_path):
    rows = _parse_repeatmasker_out(_write_run(tmp_path))
    kimura = parse_repeatmasker_align_kimura(tmp_path / "scaf_rm.fa.align")
    apply_kimura_to_rows(rows, kimura)
    landscape = _landscape_from_rows(rows, div_key="landscape_div")
    l0 = {b: bp for fam, b, bp in landscape if fam == "L0"}
    # Hit 1 (1-150, 9.1 %) owns the 101-150 overlap; hit 2 keeps 151-250.
    assert l0 == {9.0: 150, 13.0: 100}


def test_landscape_totals_match_summary_covered_bp(tmp_path):
    rows = _parse_repeatmasker_out(_write_run(tmp_path))
    families = {f.family: f.covered_bp for f in _aggregate_rm_rows(rows, 1000)}
    totals: dict[str, int] = {}
    for fam, _b, bp in _landscape_from_rows(rows):
        totals[fam] = totals.get(fam, 0) + bp
    assert totals == families


def test_nested_and_disjoint_hits():
    rows = [
        {"family": "F", "scaffold": "s", "start": 1, "end": 100, "divergence": 20.0},
        {"family": "F", "scaffold": "s", "start": 30, "end": 40, "divergence": 5.0},
        {"family": "F", "scaffold": "s", "start": 500, "end": 509, "divergence": 20.0},
        {"family": "F", "scaffold": "t", "start": 30, "end": 40, "divergence": 5.0},
    ]
    got = dict((b, bp) for _f, b, bp in _landscape_from_rows(rows))
    # s: 11 bp at 5 %, 89 + 10 bp at 20 %; t: 11 bp at 5 %.
    assert got == {5.0: 22, 20.0: 99}


def test_bin_start_edges():
    assert _bin_start(8.99, 1.0, 50.0) == 8.0
    assert _bin_start(9.0, 1.0, 50.0) == 9.0
    assert _bin_start(-1.0, 1.0, 50.0) == 0.0
    assert _bin_start(50.0, 1.0, 50.0) == 49.0  # ceiling stays in range
    assert _bin_start(0.3, 0.1, 50.0) == pytest.approx(0.3)


def test_tsv_records_bin_end_and_metric(tmp_path):
    path = _write_landscape_tsv(
        [("L0", 9.0, 150)], tmp_path / "l.tsv", bin_width=1.0, metric="kimura"
    )
    assert path.read_text().splitlines() == [
        "family\tdivergence_bin\tbp\tdivergence_bin_end\tmetric",
        "L0\t9.0\t150\t10.0\tkimura_cpg_adjusted",
    ]


@pytest.mark.parametrize(
    "metric, label",
    [
        ("kimura", "Kimura divergence, CpG-adjusted (%)"),
        ("raw", "Divergence from consensus, raw (%)"),
    ],
)
def test_html_centres_bars_and_labels_metric(tmp_path, metric, label):
    result = QuantifyResult(
        engine="repeatmasker",
        target_length=1000,
        landscape=[("L0", 9.0, 150), ("L0", 13.0, 100)],
        landscape_metric=metric,
    )
    html = write_landscape_html(result, tmp_path / "l.html").read_text()
    assert label in html
    assert '"x":[9.5,13.5]' in html
    # Plotly JSON-escapes the en dash in hover labels.
    assert "9–10%" in html or "9\\u201310%" in html


def _kimura_rows(tmp_path: Path) -> list[dict]:
    rows = _parse_repeatmasker_out(_write_run(tmp_path))
    kimura = parse_repeatmasker_align_kimura(tmp_path / "scaf_rm.fa.align")
    apply_kimura_to_rows(rows, kimura)
    return rows


def test_table_divergence_is_kimura_identity_stays_raw(tmp_path):
    fams = {f.family: f for f in _aggregate_rm_rows(_kimura_rows(tmp_path), 1000)}
    l0 = fams["L0"]
    # Length-weighted over two 150 bp hits.
    assert l0.mean_divergence == pytest.approx((9.10 + 13.55) / 2, abs=0.01)
    assert l0.mean_identity == pytest.approx(100 - (8.2 + 12.0) / 2, abs=0.01)
    assert fams["L1"].mean_divergence == pytest.approx(
        jukes_cantor_percent(20.0), abs=0.01
    )


def test_table_without_align_keeps_raw(tmp_path):
    rows = _parse_repeatmasker_out(_write_run(tmp_path, with_align=False))
    l0 = {f.family: f for f in _aggregate_rm_rows(rows, 1000)}["L0"]
    assert l0.mean_divergence == pytest.approx(10.1, abs=0.01)
    assert l0.mean_identity == pytest.approx(89.9, abs=0.01)


@pytest.mark.parametrize(
    "engine, metric, expected",
    [
        ("repeatmasker", "kimura", "kimura_cpg_adjusted"),
        ("repeatmasker", "raw", "raw_perc_div"),
        ("blastn", "raw", ""),
    ],
)
def test_summary_tsv_records_divergence_metric(tmp_path, engine, metric, expected):
    result = QuantifyResult(
        engine=engine,
        target_length=1000,
        families=_aggregate_rm_rows(_kimura_rows(tmp_path), 1000),
        landscape_metric=metric,
    )
    lines = _write_summary_tsv(result, tmp_path / "s.tsv").read_text().splitlines()
    header = lines[0].split("\t")
    assert header[-1] == "divergence_metric"
    assert all(line.split("\t")[-1] == expected for line in lines[1:])
    assert all(len(line.split("\t")) == len(header) for line in lines[1:])


@pytest.mark.parametrize(
    "metric, label", [("kimura", "mean_divergence (Kimura)"), ("raw", "mean_divergence (raw)")]
)
def test_summary_html_labels_divergence_metric(tmp_path, metric, label):
    result = QuantifyResult(
        engine="repeatmasker",
        target_length=1000,
        families=_aggregate_rm_rows(_kimura_rows(tmp_path), 1000),
        landscape_metric=metric,
    )
    html = write_summary_html([result], tmp_path / "s.html").read_text()
    assert f"<th>{label}</th>" in html
