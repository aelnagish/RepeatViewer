"""Trait table for the Tree tab (pure): parsing, arrays, methylation, colours."""

from __future__ import annotations

import gzip

import pytest

from repeatviewer.backends import tree_traits as tt


def test_paths_follow_rule_one(tmp_path):
    tree = tmp_path / "sat.fa_tree.treefile"
    assert tt.traits_path(tree).name == "sat.fa_tree.treefile_traits.tsv"
    details = tmp_path / "sat.fa_tree_details"
    details.mkdir()
    (details / "sat.fa_tree.names.tsv").write_text("tree_name\toriginal_header\tsource_file\n")
    assert tt.names_table_path(tree) == details / "sat.fa_tree.names.tsv"
    assert tt.names_table_path(tmp_path / "sat.fa_tree_masked.treefile") == details / "sat.fa_tree.names.tsv"
    assert tt.names_table_path(tmp_path / "other.nwk") is None


@pytest.mark.parametrize("header,expected", [
    ("Chr1:12345-12500(+)", ("Chr1", 12345, 12500, "+")),
    ("Chr1:12,345-12,500 some description", ("Chr1", 12345, 12500, "")),
    ("scaffold_1:900-100(-)", ("scaffold_1", 100, 900, "-")),
    ("Chr1_12345-12500_+", ("Chr1", 12345, 12500, "+")),     # tree-safe tip name
    ("scaffold_1_12345-12500", ("scaffold_1", 12345, 12500, "")),
])
def test_parse_location(header, expected):
    loc = tt.parse_location(header)
    assert (loc.chrom, loc.start, loc.end, loc.strand) == expected


def test_parse_location_rejects_plain_names():
    assert tt.parse_location("CroSat3_monomer_7") is None


def test_arrays_group_nearby_copies():
    locs = {
        "a": tt.Location("Chr1", 1000, 1156), "b": tt.Location("Chr1", 1160, 1316),
        "c": tt.Location("Chr1", 50_000, 50_156), "d": tt.Location("Chr2", 1000, 1156),
    }
    arrays = tt.assign_arrays(locs, gap=10_000)
    assert arrays["a"] == arrays["b"] == "Chr1:1.0-1.3 kb"
    assert arrays["c"] != arrays["a"] and arrays["d"].startswith("Chr2")


def test_bismark_methylation_over_intervals(tmp_path):
    cov = tmp_path / "x_bismark.cov.gz"
    with gzip.open(cov, "wt") as fh:
        fh.write("Chr1\t100\t100\t100.0\t4\t0\n")   # inside a
        fh.write("Chr1\t150\t150\t0.0\t0\t4\n")     # inside a
        fh.write("Chr1\t900\t900\t50.0\t1\t1\n")    # inside b
        fh.write("Chr9\t100\t100\t100.0\t9\t0\n")   # other chromosome
    locs = {"a": tt.Location("Chr1", 90, 200), "b": tt.Location("Chr1", 800, 1000),
            "c": tt.Location("Chr1", 5000, 6000)}
    assert tt.methylation_format(cov) == "bismark"
    levels = tt.interval_methylation(cov, locs, chunk=2)
    assert levels == {"a": 50.0, "b": 50.0, "c": None}


def test_bedmethyl_with_space_separated_counts(tmp_path):
    bed = tmp_path / "x.bedmethyl"
    row = "Chr1\t{p}\t{q}\t{code}\t10\t+\t{p}\t{q}\t255,0,0\t10 {pct} {m} {u} 0 0 0 0 0\n"
    bed.write_text(row.format(p=100, q=101, code="m", pct="80.0", m=8, u=2)
                   + row.format(p=100, q=101, code="h", pct="10.0", m=1, u=9)   # 5hmC: ignored
                   + row.format(p=150, q=151, code="m", pct="20.0", m=2, u=8))
    assert tt.methylation_format(bed) == "bedmethyl"
    levels = tt.interval_methylation(bed, {"a": tt.Location("Chr1", 90, 200)})
    assert levels["a"] == pytest.approx(50.0)


def test_build_table_from_a_tree_run(tmp_path):
    tree = tmp_path / "sat.fa_tree.treefile"
    details = tmp_path / "sat.fa_tree_details"
    details.mkdir()
    (details / "sat.fa_tree.names.tsv").write_text(
        "tree_name\toriginal_header\tsource_file\n"
        "Chr1_100-256_+\tChr1:100-256(+)\tsat.fa\n"
        "Chr1_260-416_-\tChr1:260-416(-)\tsat.fa\n"
        "Chr2_9000-9156_+\tChr2:9000-9156(+)\tsat.fa\n")
    tips = ["Chr1_100-256_+", "Chr1_260-416_-", "Chr2_9000-9156_+"]
    table = tt.build_trait_table(tree, tips)
    assert table.columns == ["header", "chromosome", "start", "end", "strand", "array"]
    assert table.rows["Chr1_260-416_-"]["strand"] == "-"
    assert table.rows["Chr1_100-256_+"]["array"] == table.rows["Chr1_260-416_-"]["array"]
    out = tt.write_trait_table(tt.traits_path(tree), table)
    back = tt.read_trait_table(out)
    assert back.rows == table.rows and back.columns == table.columns


def test_rebuild_keeps_columns_the_user_added(tmp_path):
    tips = ["Chr1:1-100", "Chr1:200-300"]
    user = tt.TraitTable(["chromosome", "karyotype"], {t: {"karyotype": "B"} for t in tips})
    table = tt.build_trait_table(tmp_path / "x.nwk", tips, existing=user)
    assert "karyotype" in table.columns and table.rows[tips[0]]["karyotype"] == "B"


def test_several_source_files_give_source_and_fish(tmp_path):
    tree = tmp_path / "multi.fa_multi_tree.treefile"
    details = tmp_path / "multi.fa_multi_tree_details"
    details.mkdir()
    (details / "multi.fa_multi_tree.names.tsv").write_text(
        "tree_name\toriginal_header\tsource_file\nA\tA\tCroSat1.fa\nB\tB\tCroSat2.fa\n")
    table = tt.build_trait_table(tree, ["A", "B"], fish_files={"CroSat1.fa"})
    assert table.rows["A"]["fish"] == "yes" and table.rows["B"]["fish"] == "no"
    assert table.rows["B"]["source_file"] == "CroSat2.fa"


def test_colour_scales():
    table = tt.TraitTable(["chrom", "meth"], {
        "a": {"chrom": "Chr1", "meth": "10"}, "b": {"chrom": "Chr1", "meth": "90"},
        "c": {"chrom": "Chr2", "meth": "NA"}, "d": {"chrom": "", "meth": "50"}})
    cat = tt.colour_scale(table, "chrom", ["a", "b", "c", "d"])
    assert not cat.numeric and cat.categories[0] == ("Chr1", tt.CATEGORICAL[0])
    assert cat.colors["a"] == cat.colors["b"] != cat.colors["c"] and "d" not in cat.colors
    num = tt.colour_scale(table, "meth", ["a", "b", "c", "d"])
    assert num.numeric and (num.low, num.high) == (10.0, 90.0) and num.n_missing == 1
    assert num.colors["a"] == tt.viridis(0.0) and num.colors["b"] == tt.viridis(1.0)


def test_many_categories_fold_into_other():
    tips = [f"t{i}" for i in range(25)]
    table = tt.TraitTable(["x"], {t: {"x": t} for t in tips})
    scale = tt.colour_scale(table, "x", tips)
    assert len(scale.categories) == len(tt.CATEGORICAL) + 1
    assert scale.categories[-1] == ("8 more values", tt.OTHER_COLOR)


def test_array_labels_stay_informative():
    assert tt._span_label("Chr1", 12_400_000, 12_401_000) == "Chr1:12.40 Mb (1.0 kb)"
    assert tt._span_label("Chr1", 12_400_000, 12_910_000) == "Chr1:12.40-12.91 Mb"
    assert tt._span_label("Chr1", 1_000, 1_316) == "Chr1:1.0-1.3 kb"
