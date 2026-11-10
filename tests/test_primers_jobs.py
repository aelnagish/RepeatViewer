"""Tools > Primers jobs: result files, annotations, Map Primers, Primer Database."""

from __future__ import annotations

import json

import pytest

pytest.importorskip("primer3")

from primer_testdata import make_data  # noqa: E402

from repeatviewer.backends.primers import design as dz  # noqa: E402
from repeatviewer.backends.primers import jobs as pj  # noqa: E402
from repeatviewer.backends.primers.oligos import Primer, PrimerDatabase  # noqa: E402
from repeatviewer.backends.primers.thermo import PrimerError  # noqa: E402


@pytest.fixture()
def data(tmp_path):
    return make_data(tmp_path / "d", genome_flank=20_000)


def _s(mode: str, **kw) -> dz.DesignSettings:
    return dz.DesignSettings(mode=mode, **{**dz.mode_preset(mode), **kw})


def _features(fasta):
    from repeatviewer.backends.sequence_features import load_sequence_map_data

    return load_sequence_map_data(fasta).features_by_id


def test_design_writes_report_next_to_input(data):
    main, files, payload = pj.run_design([data.gene], _s(dz.MODE_PCR))
    assert main == data.gene.parent / "gene.fa_primers.html"
    details = data.gene.parent / "gene.fa_primers_details"
    names = {f.name for f in files}
    assert {"gene.fa_primers.html", "gene.fa_primers_primers.fa", "gene.fa_primers_pairs.tsv",
            "gene.fa_primers_gel.svg", "gene.fa_primers.json"} <= names
    assert all(f.parent in (data.gene.parent, details) for f in files)
    html = main.read_text(encoding="utf-8")
    assert "<h2>Primers</h2>" in html and "Virtual gel" in html and "Methods text" in html
    result = pj.load_result(details / "gene.fa_primers.json")
    assert result["format"] == pj.RESULT_FORMAT and len(result["pairs"]) == 5
    assert "—" not in html


def test_second_run_numbers_or_reuses_confirmed_path(data):
    main1, _f, _p = pj.run_design([data.gene], _s(dz.MODE_PCR))
    main2, _f, _p = pj.run_design([data.gene], _s(dz.MODE_PCR))
    assert main2.name == "gene.fa_primers_2.html"
    # Overwrite confirmed by the user: same name, same details folder.
    main3, _f, _p = pj.run_design([data.gene], _s(dz.MODE_PCR), report_path=main1)
    assert main3 == main1
    renamed = data.gene.parent / "my_pcr.html"
    main4, files4, _p = pj.run_design([data.gene], _s(dz.MODE_PCR), report_path=renamed)
    assert main4 == renamed and (data.gene.parent / "my_pcr_details").is_dir()
    assert pj.default_report_path(data.gene, pj.DESIGN_TAG).name == "gene.fa_primers.html"


def test_design_annotates_input_and_replaces_old_primers(data):
    pj.run_design([data.gene], _s(dz.MODE_PCR))
    feats = _features(data.gene)["gene1"]
    primers = [f for f in feats if f.source == pj.ANNOTATION_SOURCE]
    assert {f.feature_type for f in primers} == {pj.FEATURE_PRIMER, pj.FEATURE_PRODUCT}
    assert len([f for f in primers if f.feature_type == pj.FEATURE_PRIMER]) == 2
    # A second run replaces, never stacks, the primer annotations.
    pj.run_design([data.gene], _s(dz.MODE_PCR))
    again = [f for f in _features(data.gene)["gene1"] if f.source == pj.ANNOTATION_SOURCE]
    assert len(again) == len(primers)


def test_design_without_annotation_leaves_input_alone(data):
    pj.run_design([data.gene], _s(dz.MODE_PCR), annotate=False)
    assert not (data.gene.parent / "gene.seqmap_RV.json").exists()


@pytest.mark.parametrize("mode,attr,kw", [
    (dz.MODE_QPCR, "gene", {"pick_probe": True}),
    (dz.MODE_LADDER, "monomers", {}),
    (dz.MODE_FISH, "monomers", {}),
])
def test_every_mode_runs_end_to_end(data, mode, attr, kw):
    path = getattr(data, attr)
    main, files, payload = pj.run_design([path], _s(mode, **kw))
    assert main.is_file() and payload["pairs"], payload.get("explain")
    assert payload["mode"] == mode


def test_qpcr_reference_gene_pairs(data):
    s = _s(dz.MODE_QPCR, pick_probe=False, reference_fasta=str(data.refgene))
    _main, _files, payload = pj.run_design([data.gene], s)
    groups = {p["group"] for p in payload["pairs"]}
    assert "reference" in groups and "" in groups


def test_save_pairs_to_database(data, tmp_path):
    _main, _files, payload = pj.run_design([data.gene], _s(dz.MODE_PCR))
    db_path = tmp_path / "db.fa"
    added = pj.save_pairs_to_database(payload, db_path)
    assert len(added) == 10
    again = pj.save_pairs_to_database(payload, db_path)
    assert again == []  # identical sequences are not added twice
    kinds = {p.kind for p in PrimerDatabase(db_path).primers}
    assert kinds == {"forward", "reverse"}


# --------------------------------------------------------------------------
# Map Primers
# --------------------------------------------------------------------------


def test_map_typed_pair_finds_the_product(data):
    ms = pj.MapSettings(use_forward=True, forward=data.gene_f, use_reverse=True,
                        reverse=data.gene_r, max_mismatches=0)
    main, files, payload = pj.run_map([data.gene], ms, database_path=data.root / "none.fa")
    assert main.name == "gene.fa_primer_map.html"
    # Forward = gene[100:121], reverse site = gene[1400:1421] (left 800 + repeat 400 + 200).
    assert [(p["start"], p["end"], p["size"]) for p in payload["products"]] == [(101, 1421, 1321)]
    prod = payload["products"][0]
    assert prod["forward"] == "Forward primer" and prod["reverse"] == "Reverse primer"
    assert data.gene_seq.upper()[prod["start"] - 1:prod["start"] - 1 + 21] == data.gene_f
    assert payload["counts"]["sites"] == 2
    feats = [f for f in _features(data.gene)["gene1"] if f.source == pj.ANNOTATION_SOURCE]
    assert {f.feature_type for f in feats} == {pj.FEATURE_PRIMER, pj.FEATURE_PRODUCT}


def test_map_saved_primers_from_database_and_pairs_only(data, tmp_path):
    db = PrimerDatabase(tmp_path / "db.fa")
    db.add([Primer("Gene_F", data.gene_f), Primer("Gene_R", data.gene_r),
            Primer("Sat_F", data.sat_f), Primer("Sat_R", data.sat_r)])
    db.save()
    ms = pj.MapSettings(mode=pj.MAP_DATABASE, pairs_only=True, max_mismatches=1)
    _main, _files, payload = pj.run_map([data.gene], ms, database_path=db.path)
    used = {s["primer"] for s in payload["sites"]}
    assert used == {"Gene_F", "Gene_R"}  # satellite primers do not bind gene.fa
    _main, _files, arr = pj.run_map([data.array], ms, database_path=db.path)
    sizes = sorted({p["size"] for p in arr["products"]})
    # Ladder on the tandem array: steps of one monomer.
    assert len(sizes) >= 3 and all((b - a) % len(data.mono) in (0, 1, len(data.mono) - 1)
                                   for a, b in zip(sizes, sizes[1:]))


def test_map_on_consensus_of_alignment(data):
    names = [line[1:].split()[0] for line in data.te_aln.read_text().splitlines()
             if line.startswith(">")]
    assert len(names) == 12
    from repeatviewer.backends.primers import binding as bd

    cons = bd.consensus([s for _n, s in bd.iter_fasta(data.te_aln)])
    ms = pj.MapSettings(use_forward=True, forward=cons[100:122], use_reverse=False,
                        test_on=pj.TEST_ON_CONSENSUS, max_mismatches=0)
    _main, _files, payload = pj.run_map([data.te_aln], ms, database_path=data.root / "x.fa")
    assert [t["id"] for t in payload["templates"]] == ["te_family_aln.fa consensus"]
    assert payload["counts"]["sites"] == 1
    # The consensus is not a file: no annotations written.
    assert not (data.te_aln.parent / "te_family_aln.seqmap_RV.json").exists()


def test_map_genome_counts(data):
    ms = pj.MapSettings(use_forward=True, forward=data.sat_f, use_reverse=True, reverse=data.sat_r,
                        check_genome=True, genome_fasta=str(data.genome), max_mismatches=3)
    _main, _files, payload = pj.run_map([data.array], ms, database_path=data.root / "x.fa")
    assert payload["genome"]["products"] > 20


def test_map_errors_are_clear(data):
    with pytest.raises(PrimerError, match="choose a saved primer or type"):
        pj.run_map([data.gene], pj.MapSettings(use_forward=True, forward=""),
                   database_path=data.root / "x.fa")
    with pytest.raises(PrimerError, match="neither a saved primer nor a sequence"):
        pj.run_map([data.gene], pj.MapSettings(use_forward=True, forward="Nope_F"),
                   database_path=data.root / "x.fa")


def test_result_json_is_plain_json(data):
    main, _files, _payload = pj.run_design([data.gene], _s(dz.MODE_PCR))
    raw = (data.gene.parent / "gene.fa_primers_details" / "gene.fa_primers.json").read_text()
    assert json.loads(raw)["kind"] == "design"
