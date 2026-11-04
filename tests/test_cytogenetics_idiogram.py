"""Idiogram logic: pooling FISH signals, assembly arrays, comparison, settings."""

from __future__ import annotations

import importlib.util

import numpy as np
import pytest

from repeatviewer.backends import cyto_images as ci
from repeatviewer.backends.cytogenetics import idiogram as idg
from repeatviewer.backends.cytogenetics import karyotype as kt
from repeatviewer.backends.cytogenetics import model as md
from repeatviewer.backends.cytogenetics import signals_session as ss
from repeatviewer.backends.cytogenetics.importer import ImportRequest, run_request

tifffile = pytest.importorskip("tifffile")


def _folder(tmp_path, cells=2):
    src = tmp_path / "scope"
    src.mkdir()
    for n in range(1, cells + 1):
        for dye in ("DAPI", "FITC"):
            tifffile.imwrite(str(src / f"s1_c0{n}_{dye}.tif"), np.full((120, 120), 100, dtype=np.uint16))
    work = tmp_path / "data"
    work.mkdir()
    run_request(ImportRequest(folder=work, source=src, slide_name="S1", probes_by_dye={"FITC": "CroSat2"}))
    k = kt.new_karyotype("Cs", 2, 2, species="Crocus sativus")
    k.chromosomes[0].centromere_index, k.chromosomes[0].length = 40, 10
    k.chromosomes[1].centromere_index, k.chromosomes[1].length = 50, 5
    k.chromosomes[0].assembly_seq, k.chromosomes[1].assembly_seq = "Chr1", "Chr2"
    kt.save_karyotype(work, k)
    return work, k


def _mark_cells(work, *, both_homologues=True):
    for index in range(len(ci.load_images(work))):
        s = ss.SignalsSession(work, ci.load_images(work)[index])
        for label, x in (("1a", 20), ("1b", 60)):
            m = s.add_mark(label, x, 5)
            s.set_axis(m.id, (x, 0), (x, 100))  # centromere at y = 40
            if label == "1a" or both_homologues:
                s.add_signal("CroSat2", x, 99, refine=False)  # q terminal
        s.save()


def test_signal_fraction_roundtrip():
    sig = md.Signal(id="s", probe="P", x=0, y=0, arm="q", rel_distance=0.5)
    assert idg.signal_fraction(sig, 0.4) == pytest.approx(0.7)
    assert idg.fraction_to_arm(0.7, 0.4) == ("q", pytest.approx(0.5))
    assert idg.signal_fraction(md.Signal(id="c", probe="P", x=0, y=0, arm="cen"), 0.3) == 0.3
    assert idg.signal_fraction(md.Signal(id="n", probe="P", x=0, y=0), 0.3) is None
    p = md.Signal(id="p", probe="P", x=0, y=0, arm="p", position="terminal")
    assert idg.signal_fraction(p, 0.4) == pytest.approx(0.4 - 0.97 * 0.4)


def test_collect_fish_pools_homologues(tmp_path):
    work, k = _folder(tmp_path)
    _mark_cells(work, both_homologues=False)
    data = idg.collect_fish(work, k)
    assert data.cells == 2 and data.homologues == 4 and data.unplaced == 0
    [band] = data.bands
    assert (band.chromosome, band.probe, band.arm, band.position) == ("1", "CroSat2", "q", "terminal")
    assert (band.seen, band.total, band.cells) == (2, 4, 2) and band.label == "2/4"
    assert band.fraction == pytest.approx(0.99, abs=0.01)
    assert idg.visible_bands(data, set(), 0.5) == [band]
    assert idg.visible_bands(data, set(), 0.75) == []
    assert idg.visible_bands(data, {"CroSat2"}, 0.0) == []
    per = idg.collect_fish(work, k, per_homologue=True)
    assert [(b.homologue, b.seen, b.total) for b in per.bands] == [("a", 2, 2)]


def test_merge_hits_and_silico():
    arrays = idg.merge_hits([(1000, 2000), (2500, 4000), (900_000, 901_000)], 1_000_000, gap=1000)
    assert arrays == [(1000, 4000, 2502, 2), (900_000, 901_000, 1001, 1)]
    k = kt.new_karyotype("Cs", 1, 2)
    k.chromosomes[0].assembly_seq = "Chr1"
    hits = {"CroSat2": {"Chr1": [(980_000, 990_000), (100, 200)]}}
    out = idg.collect_silico(k, {"Chr1": 1_000_000}, hits, min_array_kb=5)
    assert len(out) == 1 and out[0].start == pytest.approx(0.98, abs=0.001)
    flipped = idg.collect_silico(k, {"Chr1": 1_000_000}, hits, min_array_kb=5, flips={"1": True})
    assert flipped[0].end == pytest.approx(0.02, abs=0.001)


def test_compare_statuses(tmp_path):
    work, k = _folder(tmp_path)
    _mark_cells(work)
    data = idg.collect_fish(work, k)
    settings = idg.IdiogramSettings()
    lengths = {"Chr1": 1_000_000, "Chr2": 500_000}
    agree = [idg.SilicoArray("1", "CroSat2", 0.95, 1.0, 20000, 4)]
    rows = idg.build_rows(k, settings, data, agree, lengths)
    assert [r.status for r in rows] == [idg.STATUS_AGREES, idg.STATUS_NONE]
    far = [idg.SilicoArray("1", "CroSat2", 0.0, 0.05, 20000, 4)]
    assert idg.build_rows(k, settings, data, far, lengths)[0].status == idg.STATUS_DIFFERS
    assert idg.build_rows(k, settings, data, [], lengths)[0].status == idg.STATUS_FISH
    only = [idg.SilicoArray("2", "CroSat2", 0.4, 0.5, 20000, 4)]
    assert idg.build_rows(k, settings, data, only, lengths)[1].status == idg.STATUS_ASM
    assert idg.build_rows(k, settings, data, [], {})[0].status == idg.STATUS_UNLINKED
    tsv = idg.comparison_tsv(rows)
    assert tsv.splitlines()[0].startswith("chromosome\tsequence\tstatus") and "in situ\tCroSat2" in tsv


def test_layout(tmp_path):
    work, k = _folder(tmp_path)
    s = idg.IdiogramSettings(order="length")
    assert [c.name for c in idg.ordered_chromosomes(k, s, {})] == ["1", "2"]
    s.order, s.custom_order = "custom", ["2", "1"]
    assert [c.name for c in idg.ordered_chromosomes(k, s, {})] == ["2", "1"]
    rows = idg.build_rows(k, idg.IdiogramSettings(homologues=True), idg.FishData([], 0, 0, 0, {}), [], {})
    assert [r.name for r in rows] == ["1a", "1b", "2a", "2b"]
    assert idg.relative_heights(rows[::2], "karyotype") == [1.0, 0.5]
    assert idg.relative_heights(rows[::2], "same") == [1.0, 1.0]
    assert idg.title_text(k, ["CroSat2", "45S rDNA"]) == ("CroSat2 and 45S rDNA on the", "Crocus sativus")


def test_move_band(tmp_path):
    work, k = _folder(tmp_path)
    _mark_cells(work)
    [band] = idg.collect_fish(work, k).bands
    assert idg.move_band(work, band, 0.4, 0.4) == 4
    [moved] = idg.collect_fish(work, k).bands
    assert moved.arm == "cen" and moved.fraction == pytest.approx(0.4)
    sig = md.image_details(ci.load_images(work)[0]).signals[0]
    assert sig.manual


def test_settings_roundtrip(tmp_path):
    s = idg.IdiogramSettings(karyotype="Cs", hidden_probes=["X"], min_share=0.0, genome="g.fa", order="custom",
                             custom_order=["2", "1"], lengths="same", homologues=True, flips={"1": True})
    s.source("CroSat2").match = "CroSat2_mono"
    idg.save_settings(tmp_path, s)
    again = idg.load_settings(tmp_path)
    assert again.to_dict() == s.to_dict()
    assert idg.IdiogramSettings.from_dict({"order": "bad", "lengths": 3}).order == "karyotype"


def test_features_from_file(tmp_path):
    bed = tmp_path / "hits.bed"
    bed.write_text("Chr1\t100\t20000\tCroSat2_array\nChr1\t5000\t6000\tother\n", encoding="utf-8")
    src = idg.ProbeSource(kind="file", file=str(bed))
    hits = idg.features_for_probe(tmp_path / "none.fa", ["Chr1"], "CroSat2", src, cache={})
    assert len(hits["Chr1"]) == 1
    assert idg.features_for_probe(tmp_path / "none.fa", ["Chr1"], "CroSat2", idg.ProbeSource(kind="none")) == {}


def test_genome_candidates(tmp_path):
    (tmp_path / "b.fa").write_text(">Chr1\nACGT\n")
    (tmp_path / "a.fasta").write_text(">Chr2\nACGTAC\n")
    (tmp_path / "x.txt").write_text("no")
    assert [p.name for p in idg.genome_candidates(tmp_path)] == ["a.fasta", "b.fa"]
    assert idg.resolve_genome(tmp_path, "b.fa") == tmp_path / "b.fa"
    if importlib.util.find_spec("Bio") is not None:  # Biopython is in the app, not always in a bare test env
        assert idg.sequence_lengths(tmp_path / "a.fasta") == {"Chr2": 6}
